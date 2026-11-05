"""Read-only checks for unusable ligand coordinates before engine dispatch.

These checks do not certify conformer quality, valence or stereochemistry, and
never regenerate, filter, rename or rewrite a molecule.
"""
from __future__ import annotations

import math
from pathlib import Path


def audit_sdf_geometry(path: Path) -> dict:
    """Inspect every SDF record, including ordinary libraries without state tags."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ValueError("RDKit is required for ligand coordinate checks") from exc

    report = {
        "success": False,
        "source_file": str(path),
        "scope": "basic coordinates only; not a stereochemistry or pose-quality certificate",
        "record_count": 0,
        "invalid_record_count": 0,
        "invalid_records": [],
    }
    with Path(path).open("rb") as handle:
        # Do not sanitize or remove hydrogens: preserve supplied atom indices
        # and leave chemical-state validation to its existing contract.
        supplier = Chem.ForwardSDMolSupplier(
            handle, removeHs=False, sanitize=False, strictParsing=True,
        )
        for index, molecule in enumerate(supplier):
            report["record_count"] += 1
            name = None
            reason = None
            if molecule is None:
                reason = "unparseable_record"
            else:
                if molecule.HasProp("_Name"):
                    name = molecule.GetProp("_Name")
                if molecule.GetNumAtoms() == 0:
                    reason = "no_atoms"
                elif molecule.GetNumConformers() == 0:
                    reason = "missing_coordinates"
                else:
                    conformer = molecule.GetConformer()
                    positions = [
                        tuple(conformer.GetAtomPosition(i))
                        for i in range(molecule.GetNumAtoms())
                    ]
                    if not all(math.isfinite(v) for xyz in positions for v in xyz):
                        reason = "nonfinite_coordinates"
                    elif len(positions) > 1 and all(
                        math.dist(positions[0], xyz) <= 0.001
                        for xyz in positions[1:]
                    ):
                        reason = "collapsed_coordinates"
                    elif any(
                        math.dist(
                            positions[bond.GetBeginAtomIdx()],
                            positions[bond.GetEndAtomIdx()],
                        ) <= 0.001
                        for bond in molecule.GetBonds()
                    ):
                        reason = "zero_length_bond"
            if reason:
                report["invalid_records"].append({
                    "source_molecule_index": index,
                    "record_number": index + 1,
                    "name": name,
                    "reason": reason,
                })
    report["invalid_record_count"] = len(report["invalid_records"])
    report["success"] = report["record_count"] > 0 and not report["invalid_records"]
    if not report["record_count"]:
        report["error"] = "empty_library"
    return report


def geometry_failure_message(report: dict) -> str:
    records = report["invalid_records"]
    details = "; ".join(
        f"record {row['record_number']} ({row['name'] or 'unnamed'}): {row['reason']}"
        for row in records[:10]
    )
    if len(records) > 10:
        details += f"; and {len(records) - 10} more (see geometry_qc)"
    return (
        f"Invalid ligand coordinates in {len(records)}/{report['record_count']} "
        f"SDF records: {details or report.get('error', 'unknown error')}. "
        "Regenerate failed records from the original declared SMILES and verify "
        "stereochemistry; do not infer it from collapsed coordinates."
    )

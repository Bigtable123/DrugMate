"""Per-pose ProLIF evaluation for the structure-evaluator.

The canonical ProLIF implementation lives in pocket-analysis.  This adapter
reuses it for every docked SDF record and exposes the atom-aware interaction
observations expected by the Step 4 blueprint matcher.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path


BIO_ROOT = Path(__file__).resolve().parents[2]
POCKET_SCRIPTS = BIO_ROOT / "pocket-analysis" / "scripts"
if str(POCKET_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(POCKET_SCRIPTS))

from interaction_observations import (  # noqa: E402
    PROLIF_TO_UNIFIED,
    normalize_interaction_observations,
)
from prolif_utils import run_prolif_fingerprint  # noqa: E402


def run_prolif_per_pose(
    receptor_file: Path,
    docked_sdf: Path,
    *,
    chain: str = "A",
    vicinity_cutoff: float = 6.5,
    interaction_types: list[str] | None = None,
    reference_ligand_file: Path | None = None,
) -> list[dict]:
    """Run ProLIF independently for every record in ``docked_sdf``.

    A per-pose temporary SDF is intentional: the shared Step 1 helper analyses
    one ligand at a time, and isolating records prevents a malformed pose from
    shifting identities for every following result.
    """
    from rdkit import Chem

    receptor_path = Path(receptor_file)
    supplier = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)
    results: list[dict] = []
    reference_bits, reference_error = None, None
    if reference_ligand_file:
        try:
            import hashlib
            raw_reference = run_prolif_fingerprint(
                receptor_path, Path(reference_ligand_file),
                interaction_types=interaction_types, vicinity_cutoff=float(vicinity_cutoff))
            if raw_reference.get("success") is False or raw_reference.get("validation_failed"):
                raise ValueError(_raw_error(raw_reference))
            observations = normalize_interaction_observations(None, raw_reference, default_chain=chain)
            reference_bits = residue_fingerprint(observations.get("observations") or [])
            reference_hash = hashlib.sha256(Path(reference_ligand_file).read_bytes()).hexdigest()
        except Exception as exc:
            reference_bits = None
            reference_error = str(exc)
    with tempfile.TemporaryDirectory(prefix="sbdd-prolif-") as temp_dir:
        temp_root = Path(temp_dir)
        for pose_index, mol in enumerate(supplier):
            mol_name = (
                mol.GetProp("_Name").strip()
                if mol is not None and mol.HasProp("_Name") and mol.GetProp("_Name").strip()
                else f"pose_{pose_index}"
            )
            if mol is None:
                results.append({
                    "success": False,
                    "pose_index": pose_index,
                    "mol_name": mol_name,
                    "error": "RDKit failed to parse pose",
                    "num_interactions": None,
                    "interactions": [],
                })
                continue

            pose_file = temp_root / f"pose_{pose_index:06d}.sdf"
            writer = Chem.SDWriter(str(pose_file))
            try:
                writer.write(mol)
            finally:
                writer.close()

            try:
                raw = run_prolif_fingerprint(
                    receptor_path,
                    pose_file,
                    interaction_types=interaction_types,
                    vicinity_cutoff=float(vicinity_cutoff),
                )
                normalized = normalize_interaction_observations(
                    None,
                    raw,
                    default_chain=chain,
                )
                interactions = [
                    item
                    for item in normalized.get("observations", [])
                    if item.get("tool") == "ProLIF"
                ]
                failed = bool(raw.get("validation_failed") or raw.get("success") is False)
                bits = residue_fingerprint(interactions)
                reference_union = set(bits) | set(reference_bits or [])
                results.append({
                    "success": not failed,
                    "pose_index": pose_index,
                    "mol_name": mol_name,
                    "num_interactions": len(interactions),
                    "interactions": interactions,
                    "interaction_observations": normalized,
                    "evaluator_capabilities": _capabilities(raw, interactions),
                    "raw_summary": raw.get("summary") or {},
                    "tool_version": raw.get("tool_version"),
                    "diagnostics": raw.get("diagnostics"),
                    "residue_level_bits": bits,
                    "tanimoto_to_reference": (
                        len(set(bits) & set(reference_bits)) / len(reference_union)
                        if not failed and reference_bits and reference_union else None
                    ),
                    "reference_comparison": {
                        "status": "failed" if reference_error else "available" if reference_bits else "not_applicable",
                        "reason": reference_error or (None if reference_bits else "no usable reference interaction fingerprint"),
                        "reference_sha256": reference_hash if reference_bits else None,
                        "feature_definition": "protein chain/residue/insertion code plus interaction subtype",
                    },
                    **({"error": _raw_error(raw)} if failed else {}),
                })
            except Exception as exc:
                results.append({
                    "success": False,
                    "pose_index": pose_index,
                    "mol_name": mol_name,
                    "error": str(exc),
                    "num_interactions": None,
                    "interactions": [],
                })
    return results


def residue_fingerprint(interactions: list[dict]) -> list[str]:
    """Residue/type bits independent of ligand atom indices and input order."""
    import json

    bits = set()
    for item in interactions:
        residue = item.get("residue") or item.get("residue_record")
        if not isinstance(residue, dict):
            continue
        bits.add(json.dumps([
            residue.get("chain"), residue.get("resnum"), residue.get("icode") or "",
            residue.get("resname"), item.get("interaction_subtype") or item.get("interaction_type"),
        ], separators=(",", ":")))
    return sorted(bits)


def _capabilities(raw: dict, interactions: list[dict]) -> dict:
    """Build a conservative, machine-readable ProLIF capability manifest."""
    raw_types = (
        (raw.get("capabilities") or {}).get("enabled_pairwise_interaction_types")
        or []
    )
    interaction_types = sorted({
        PROLIF_TO_UNIFIED.get(str(value), str(value))
        for value in raw_types
        if value
    })
    fields_by_type: dict[str, set[str]] = {}
    for observation in interactions:
        interaction_type = str(observation.get("interaction_type") or "")
        if not interaction_type:
            continue
        fields = fields_by_type.setdefault(interaction_type, set())
        for field in (
            "ligand_role",
            "protein_role",
            "protein_atoms",
            "ligand_atoms",
            "mediator",
        ):
            if observation.get(field) not in (None, "", [], {}):
                fields.add(field)
        for key, value in (observation.get("geometry") or {}).items():
            if value not in (None, "", [], {}):
                fields.add(f"geometry.{key}")
    return {
        "source": "prolif",
        "interaction_types": interaction_types,
        "fields": [],
        "fields_by_type": {
            key: sorted(values)
            for key, values in sorted(fields_by_type.items())
        },
    }


def _raw_error(raw: dict) -> str:
    diagnostics = raw.get("diagnostics") or {}
    return str(
        diagnostics.get("error")
        or diagnostics.get("detail")
        or raw.get("error")
        or "ProLIF pose evaluation failed"
    )

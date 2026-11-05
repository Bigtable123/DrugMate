"""PDBe Arpeggio wrapper with strict mmCIF and chemistry-quality checks."""
from __future__ import annotations

import json
import math
import subprocess
import tempfile
from pathlib import Path
from typing import Any


_PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "ASX", "GLX", "SEC", "PYL", "MSE",
}

_METAL_RESNAMES = {
    "LI", "NA", "K", "RB", "CS", "BE", "MG", "CA", "SR", "BA",
    "AL", "GA", "IN", "SN", "PB", "BI", "SC", "TI", "V", "CR",
    "MN", "FE", "CO", "NI", "CU", "ZN", "Y", "ZR", "NB", "MO",
    "TC", "RU", "RH", "PD", "AG", "CD", "HF", "TA", "W", "RE",
    "OS", "IR", "PT", "AU", "HG", "LA", "CE", "PR", "ND", "SM",
    "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
}

_CONTACT_TYPES = {
    "hbond": "HBond",
    "weak_hbond": "WeakHBond",
    "xbond": "HalogenBond",
    "ionic": "SaltBridge",
    "metal": "MetalComplex",
    "metal_complex": "MetalComplex",
    "hydrophobic": "Hydrophobic",
    "vdw": "VdW",
}

_CHEMISTRY_WARNING_MARKERS = (
    "failed to kekulize",
    "perceivebondorders",
    "failed to perceive bond",
    "aromatic bond",
)


def _version(distribution: str) -> str | None:
    try:
        from importlib.metadata import version

        return version(distribution)
    except Exception:
        return None


def _clean_icode(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text in {".", "?", "\x00"} else text


def _endpoint_identity(endpoint: dict) -> tuple[str, int | None, str, str]:
    try:
        number = int(endpoint.get("auth_seq_id"))
    except (TypeError, ValueError):
        number = None
    return (
        str(endpoint.get("auth_asym_id") or ""),
        number,
        _clean_icode(endpoint.get("pdbx_PDB_ins_code")),
        str(endpoint.get("label_comp_id") or "").upper(),
    )


def _atom_from_endpoint(endpoint: dict, namespace: str, *, protein: bool) -> dict:
    name = str(endpoint.get("auth_atom_id") or "").strip()
    element = ""
    for char in name:
        if char.isalpha():
            element += char
        elif element:
            break
    # Protein atom names such as CA/CB encode topology, not a two-letter
    # element. Ligand atom names may genuinely start with CL/BR. Atom name is
    # still the primary identity; this field is only a conservative hint.
    if protein and element:
        element = element[:1]
    elif len(element) > 1 and element[:2].upper() not in {
        "CL", "BR", "SI", "SE", "NA", "MG", "AL", "CA", "FE", "ZN",
    }:
        element = element[:1]
    return {
        "name": name,
        "element": element.upper(),
        "namespace": namespace,
    }


def _sdf_centroid(ligand_file: Path) -> tuple[float, float, float] | None:
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(ligand_file), removeHs=False, sanitize=False)
        molecule = next((mol for mol in supplier if mol is not None), None)
        if molecule is None or molecule.GetNumConformers() == 0:
            return None
        conformer = molecule.GetConformer()
        coords = [conformer.GetAtomPosition(i) for i in range(molecule.GetNumAtoms())]
        if not coords:
            return None
        return tuple(sum(getattr(pos, axis) for pos in coords) / len(coords) for axis in ("x", "y", "z"))
    except Exception:
        return None


def _mmcif_ligand_candidates(
    complex_file: Path,
    ligand_resname: str,
) -> tuple[list[dict], dict[str, str], list[str]]:
    import gemmi

    document = gemmi.cif.read_file(str(complex_file))
    block = document.sole_block()
    comp_ids = list(block.find_values("_chem_comp.id"))
    comp_types = list(block.find_values("_chem_comp.type"))
    type_by_comp = {
        str(comp_id).upper(): str(comp_type)
        for comp_id, comp_type in zip(comp_ids, comp_types)
    }
    warnings = []
    if not comp_ids:
        warnings.append("mmCIF is missing the _chem_comp category")
    ligand_type = type_by_comp.get(ligand_resname.upper(), "")
    if not ligand_type or ligand_type.strip() in {".", "?"}:
        warnings.append(f"mmCIF has no usable _chem_comp.type for ligand {ligand_resname}")

    structure = gemmi.read_structure(str(complex_file))
    candidates = []
    if len(structure) == 0:
        return candidates, type_by_comp, warnings + ["mmCIF has no coordinate model"]
    for chain in structure[0]:
        for residue in chain:
            if residue.name.upper() != ligand_resname.upper():
                continue
            coords = [(atom.pos.x, atom.pos.y, atom.pos.z) for atom in residue if atom.element.name != "H"]
            if not coords:
                continue
            candidates.append({
                "chain": chain.name,
                "resnum": int(residue.seqid.num),
                "icode": _clean_icode(residue.seqid.icode),
                "resname": residue.name.upper(),
                "centroid": [sum(point[i] for point in coords) / len(coords) for i in range(3)],
            })
    if not candidates:
        warnings.append(f"ligand {ligand_resname} was not found in mmCIF coordinates")
    return candidates, type_by_comp, warnings


def _select_ligand_identity(
    candidates: list[dict],
    ligand_file: Path,
    preferred_chain: str | None,
) -> dict | None:
    if not candidates:
        return None
    pool = [item for item in candidates if item["chain"] == preferred_chain] or candidates
    target = _sdf_centroid(ligand_file)
    if target is None:
        return sorted(pool, key=lambda item: (item["chain"], item["resnum"], item["icode"]))[0]
    return min(
        pool,
        key=lambda item: math.dist(item["centroid"], target),
    )


def _is_protein_endpoint(endpoint: dict) -> bool:
    comp_type = str(endpoint.get("label_comp_type") or "").upper()
    resname = str(endpoint.get("label_comp_id") or "").upper()
    return comp_type == "P" or resname in _PROTEIN_RESNAMES


def _is_metal_endpoint(endpoint: dict) -> bool:
    return str(endpoint.get("label_comp_id") or "").upper() in _METAL_RESNAMES


def parse_arpeggio_records(
    records: list[dict],
    ligand_identity: dict,
    *,
    blueprint_eligible: bool,
) -> tuple[list[dict], dict]:
    """Convert Arpeggio JSON records to the shared profiler detail contract."""
    ligand_key = (
        ligand_identity["chain"], int(ligand_identity["resnum"]),
        _clean_icode(ligand_identity.get("icode")), ligand_identity["resname"].upper(),
    )
    details = []
    qc_counts = {"clash": 0, "vdw_clash": 0, "covalent": 0, "proximal": 0}
    atom_atom_records = 0
    direct_records = 0
    for record in records:
        if not isinstance(record, dict) or record.get("interacting_entities") != "INTER":
            continue
        bgn = record.get("bgn") or {}
        end = record.get("end") or {}
        bgn_is_ligand = _endpoint_identity(bgn) == ligand_key
        end_is_ligand = _endpoint_identity(end) == ligand_key
        if bgn_is_ligand == end_is_ligand:
            continue
        ligand = bgn if bgn_is_ligand else end
        protein = end if bgn_is_ligand else bgn
        protein_endpoint = _is_protein_endpoint(protein)
        metal_endpoint = _is_metal_endpoint(protein)
        if not protein_endpoint and not metal_endpoint:
            continue
        direct_records += 1
        is_atom_atom = record.get("type") == "atom-atom"
        atom_atom_records += int(is_atom_atom)
        contacts = {str(item).lower() for item in (record.get("contact") or [])}
        for qc_name in ("clash", "vdw_clash", "covalent"):
            qc_counts[qc_name] += int(qc_name in contacts)
        # The published fraction is proximal atom-atom contacts divided by all
        # direct atom-atom contacts.  Group/plane records must enter neither side.
        qc_counts["proximal"] += int(is_atom_atom and "proximal" in contacts)
        for contact in sorted(contacts):
            raw_type = _CONTACT_TYPES.get(contact)
            if raw_type is None:
                continue
            weak = contact == "weak_hbond"
            if metal_endpoint and raw_type != "MetalComplex":
                continue
            detail = {
                "type": raw_type,
                "residue_record": {
                    "resname": str(protein.get("label_comp_id") or "").upper(),
                    "chain": str(protein.get("auth_asym_id") or ""),
                    "resnum": protein.get("auth_seq_id"),
                    "icode": _clean_icode(protein.get("pdbx_PDB_ins_code")),
                },
                "protein_atoms": [
                    _atom_from_endpoint(
                        protein, "arpeggio_auth_atom", protein=protein_endpoint,
                    )
                ],
                "ligand_atoms": [
                    _atom_from_endpoint(
                        ligand, "arpeggio_auth_atom", protein=False,
                    )
                ],
                "geometry": {"distance": record.get("distance")},
                "mediator": "metal" if raw_type == "MetalComplex" else "direct",
                "ligand_role": "metal_coordinator" if metal_endpoint else None,
                "protein_role": "metal_center" if metal_endpoint else None,
                "blueprint_eligible": bool(blueprint_eligible and not weak),
                "quality_flags": ["weak_hbond"] if weak else [],
                "raw_metadata": {
                    "arpeggio_record_type": record.get("type"),
                    "contact_labels": sorted(contacts),
                    "interacting_entities": record.get("interacting_entities"),
                    "target_entity_class": (
                        "metal" if metal_endpoint else "protein"
                    ),
                },
            }
            details.append(detail)
    return details, {"direct_ligand_protein_records": direct_records, "atom_atom_records": atom_atom_records, **qc_counts}


def run_arpeggio_interactions(
    complex_file: Path,
    ligand_file: Path,
    ligand_resname: str,
    chain: str = "A",
    *,
    timeout_seconds: int = 300,
) -> dict:
    if complex_file.suffix.lower() not in {".cif", ".mmcif"}:
        return {
            "success": False,
            "validation_failed": True,
            "not_evaluable": True,
            "error": "PDBe Arpeggio 1.4.4 requires a source complex in mmCIF format",
            "tool_version": _version("pdbe-arpeggio"),
        }

    candidates, type_by_comp, input_warnings = _mmcif_ligand_candidates(
        complex_file, ligand_resname,
    )
    ligand_identity = _select_ligand_identity(candidates, ligand_file, chain)
    if ligand_identity is None:
        return {
            "success": False,
            "validation_failed": True,
            "not_evaluable": True,
            "error": "; ".join(input_warnings) or "ligand identity could not be resolved in mmCIF",
            "tool_version": _version("pdbe-arpeggio"),
        }

    selection_residue = f"{ligand_identity['resnum']}{ligand_identity.get('icode') or ''}"
    selection = f"/{ligand_identity['chain']}/{selection_residue}/"
    with tempfile.TemporaryDirectory(prefix="arpeggio_step1_") as tmp:
        output_dir = Path(tmp)
        command = [
            "pdbe-arpeggio",
            str(complex_file), "-s", selection, "-o", str(output_dir), "-m",
        ]
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
        output_file = output_dir / f"{complex_file.stem}.json"
        if completed.returncode != 0 or not output_file.is_file():
            detail = (completed.stderr or completed.stdout or "Arpeggio produced no JSON").strip()
            raise RuntimeError(f"pdbe-arpeggio failed with code {completed.returncode}: {detail}")
        records = json.loads(output_file.read_text(encoding="utf-8"))

    combined_log = "\n".join([completed.stdout or "", completed.stderr or ""])
    backend_chemistry_warnings = [
        marker for marker in _CHEMISTRY_WARNING_MARKERS
        if marker in combined_log.lower()
    ]
    from rdkit import Chem

    mol = next(
        (m for m in Chem.SDMolSupplier(str(ligand_file), removeHs=False) if m),
        None,
    )
    if mol is None:
        raise ValueError(f"RDKit failed to parse Arpeggio ligand: {ligand_file}")
    plane_count = sum(
        record.get("type") in {"atom-plane", "plane-plane", "group-group"}
        for record in records
    )
    chemistry_warnings = list(backend_chemistry_warnings)
    # An aromatic ligand need not make a pi/plane contact in a valid pose.  The
    # absence of such a record is interaction evidence, not a chemistry-parser
    # failure; only explicit backend/input chemistry warnings degrade the result.
    quality_flags = sorted(set(input_warnings + chemistry_warnings))
    quality = "degraded" if quality_flags else "valid"
    details, qc = parse_arpeggio_records(
        records,
        ligand_identity,
        blueprint_eligible=quality == "valid",
    )
    for detail in details:
        detail["quality_flags"] = sorted(set(detail.get("quality_flags") or []) | set(quality_flags))

    try:
        from openbabel import openbabel

        openbabel_version = openbabel.OBReleaseVersion()
    except Exception:
        openbabel_version = None
    return {
        "success": True,
        "tool": "Arpeggio",
        "tool_version": _version("pdbe-arpeggio"),
        "openbabel_version": openbabel_version,
        "quality": quality,
        "quality_flags": quality_flags,
        "blueprint_eligible": quality == "valid",
        "ligand_identity": ligand_identity,
        "selection": selection,
        "chemical_representation": {
            "source_complex": str(complex_file),
            "format": "mmCIF",
            "ligand_chem_comp_type": type_by_comp.get(ligand_resname.upper()),
            "bond_order_warning_detected": bool(backend_chemistry_warnings),
        },
        "raw_record_count": len(records),
        "raw_contact_records": records,
        "plane_or_group_record_count": plane_count,
        "num_interactions_total": len(details),
        "interactions_detail": details,
        "qc_counts": qc,
        "capabilities": {
            "atom_identity": True,
            "distance": True,
            "directional_roles": False,
            "angles": False,
            "structured_json": True,
        },
        "runtime": {
            "command": command,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        },
    }

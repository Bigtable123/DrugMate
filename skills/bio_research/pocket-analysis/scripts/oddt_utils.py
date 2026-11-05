"""ODDT 0.7 residue IFP/SPLIF wrapper for the isolated legacy environment."""
from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any


_CHANNELS = (
    ("Hydrophobic", None, None),
    ("PiStacking", None, None),
    ("EdgeToFace", None, None),
    ("HBond", "acceptor", "donor"),       # protein is donor
    ("HBond", "donor", "acceptor"),       # protein is acceptor
    ("SaltBridge", "anion", "cation"),    # protein is positive
    ("SaltBridge", "cation", "anion"),    # protein is negative
    (None, None, None),                       # upstream metal-index bug: ignore
)


def _version(distribution: str) -> str | None:
    try:
        from importlib.metadata import version

        return version(distribution)
    except Exception:
        return None


def _clean_icode(value: Any) -> str:
    text = str(value or "").replace("\x00", "").strip()
    return "" if text in {".", "?"} else text


def _residue_record(residue: Any) -> dict:
    ob_residue = residue.OBResidue
    return {
        "resname": str(ob_residue.GetName() or residue.name).upper(),
        "resnum": int(ob_residue.GetNum()),
        "chain": str(ob_residue.GetChain() or residue.chain),
        "icode": _clean_icode(ob_residue.GetInsertionCode()),
    }


def _first_molecule(path: Path, fmt: str):
    import oddt

    try:
        return next(oddt.toolkit.readfile(fmt, str(path)))
    except StopIteration as exc:
        raise ValueError(f"ODDT could not parse {fmt} file: {path}") from exc


def _finite_count(array) -> int:
    import numpy as np

    return int(np.isfinite(array).sum())


def run_oddt_fingerprints(
    protein_file: Path,
    ligand_file: Path,
    ligand_resname: str,
    *,
    strict: bool = True,
    splif_depth: int = 1,
    splif_size: int = 4096,
    splif_distance_cutoff: float = 4.5,
) -> dict:
    # ODDT chooses its backend at import time.  Pin the tested backend before
    # importing it rather than inheriting an arbitrary parent-process value.
    os.environ["ODDT_TOOLKIT"] = "ob"
    import numpy as np
    import oddt
    from oddt import fingerprints

    protein = _first_molecule(protein_file, "pdb")
    ligand_format = ligand_file.suffix.lower().lstrip(".")
    if ligand_format == "sd":
        ligand_format = "sdf"
    ligand = _first_molecule(ligand_file, ligand_format)
    protein.protein = True

    ifp = fingerprints.InteractionFingerprint(ligand, protein, strict=bool(strict))
    residues = list(protein.residues)
    expected_size = len(residues) * len(_CHANNELS)
    if int(ifp.size) != expected_size:
        raise RuntimeError(
            f"ODDT IFP length {ifp.size} does not match {len(residues)} residues x 8 channels"
        )
    matrix = np.asarray(ifp, dtype=np.uint8).reshape(len(residues), len(_CHANNELS))
    details = []
    active_residues = []
    for residue_index, counts in enumerate(matrix):
        if not np.any(counts):
            continue
        record = _residue_record(residues[residue_index])
        channel_counts = {}
        for channel_index, count_value in enumerate(counts.tolist()):
            count = int(count_value)
            if not count:
                continue
            raw_type, ligand_role, protein_role = _CHANNELS[channel_index]
            channel_counts[str(channel_index)] = count
            if raw_type is None:
                continue
            details.append({
                "type": raw_type,
                "residue_record": record,
                "ligand_role": ligand_role,
                "protein_role": protein_role,
                "count": count,
                "blueprint_eligible": False,
                "quality_flags": ["residue_level_count_without_atom_geometry"],
                "raw_metadata": {
                    "oddt_internal_residue_index": residue_index,
                    "ifp_channel": channel_index,
                    "contact_count": count,
                },
            })
        active_residues.append({
            "residue_record": record,
            "internal_residue_index": residue_index,
            "channel_counts": channel_counts,
        })

    simple_ifp = fingerprints.SimpleInteractionFingerprint(
        ligand, protein, strict=bool(strict),
    )
    splif = fingerprints.SPLIF(
        ligand,
        protein,
        depth=int(splif_depth),
        size=int(splif_size),
        distance_cutoff=float(splif_distance_cutoff),
    )
    hashes = [int(value) for value in splif["hash"].tolist()]
    try:
        from openbabel import openbabel

        openbabel_runtime_version = openbabel.OBReleaseVersion()
    except Exception:
        openbabel_runtime_version = None
    return {
        "success": True,
        "tool": "ODDT",
        "tool_version": _version("oddt"),
        "backend": getattr(oddt.toolkit, "__name__", str(oddt.toolkit)),
        "numpy_version": np.__version__,
        "openbabel_package_version": _version("openbabel"),
        "openbabel_runtime_version": openbabel_runtime_version,
        "num_interactions_total": sum(
            int(detail.get("count") or 0) for detail in details
        ),
        "interactions_detail": details,
        "interaction_fingerprint": {
            "kind": "residue_contact_counts",
            "strict": bool(strict),
            "shape": list(matrix.shape),
            "channel_order": [
                "hydrophobic", "pi_face_to_face", "pi_edge_to_face",
                "hbond_protein_donor", "hbond_protein_acceptor",
                "salt_protein_positive", "salt_protein_negative",
                "metal_ignored_upstream_index_bug",
            ],
            "active_residues": active_residues,
            "nonzero_channel_count": int(np.count_nonzero(matrix)),
            "contact_count_sum_excluding_metal": int(matrix[:, :7].sum()),
        },
        "simple_interaction_fingerprint": {
            "shape": list(np.asarray(simple_ifp).shape),
            "values": [int(value) for value in np.asarray(simple_ifp).tolist()],
        },
        "splif": {
            "depth": int(splif_depth),
            "size": int(splif_size),
            "distance_cutoff": float(splif_distance_cutoff),
            "record_count": int(len(splif)),
            "unique_hash_count": len(set(hashes)),
            "hashes": hashes,
            "finite_ligand_coordinate_values": _finite_count(splif["ligand_coords"]),
            "finite_protein_coordinate_values": _finite_count(splif["protein_coords"]),
            "interpretation": "relative similarity only; not an absolute hard gate",
        },
        "quality_flags": [
            "residue_level_counts_without_atom_geometry",
            "metal_channel_ignored_due_upstream_mapping_bug",
            "splif_similarity_is_relative_not_absolute",
        ],
        "blueprint_eligible": False,
        "capabilities": {
            "residue_level_ifp": True,
            "atom_identity": False,
            "geometry": False,
            "splif_relative_similarity": True,
            "hard_gate": False,
        },
        "chemical_representation": {
            "protein_file": str(protein_file),
            "ligand_file": str(ligand_file),
            "ligand_resname": ligand_resname,
            "hydrogens_added_by_wrapper": False,
        },
    }

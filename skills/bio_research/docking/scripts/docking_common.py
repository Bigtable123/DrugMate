"""Engine-neutral geometry helpers for Step 4 docking."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from gnina_utils import normalize_docking_box
from pose_identity import LIGAND_LINEAGE_PROPERTIES, build_pose_identity


_IDENTITY_ALIASES = {
    "source": ("source", "generator", "origin"),
    "microstate_id": ("microstate_id", "protonation_state_id", "tautomer_id"),
    "variant_id": ("variant_id", "stereoisomer_id"),
    "input_compound_id": ("input_compound_id", "parent_compound_id"),
}

POSE_IDENTITY_BASIS = "sdf_molblock_roundtrip_v1"


def reference_ligand_box(reference_ligand_file: Path, autobox_add: float) -> dict:
    """Build an explicit box matching GNINA autobox padding semantics."""
    from rdkit import Chem

    padding = float(autobox_add)
    if not math.isfinite(padding) or padding < 0:
        raise ValueError("autobox_add must be a finite number >= 0")

    supplier = Chem.SDMolSupplier(
        str(reference_ligand_file), removeHs=False, strictParsing=False,
    )
    ligand = next((mol for mol in supplier if mol is not None), None)
    if ligand is None or ligand.GetNumConformers() == 0:
        raise ValueError("reference ligand has no parseable 3D conformer")

    atom_indices = [
        atom.GetIdx() for atom in ligand.GetAtoms() if atom.GetAtomicNum() > 1
    ] or list(range(ligand.GetNumAtoms()))
    conformer = ligand.GetConformer()
    if not conformer.Is3D():
        raise ValueError("reference ligand conformer is not marked as 3D")
    axes = [[], [], []]
    for atom_index in atom_indices:
        position = conformer.GetAtomPosition(atom_index)
        for axis, value in enumerate((position.x, position.y, position.z)):
            number = float(value)
            if not math.isfinite(number):
                raise ValueError("reference ligand contains non-finite coordinates")
            axes[axis].append(number)

    center = [(min(values) + max(values)) / 2.0 for values in axes]
    # GNINA's autobox_add is added on both sides of every dimension.
    size = [
        max((max(values) - min(values)) + 2.0 * padding, 1.0)
        for values in axes
    ]
    return normalize_docking_box({
        "center": center,
        "size": size,
        "unit": "angstrom",
        "source": "reference_ligand_autobox",
        "autobox_add": padding,
    })


def compound_name(mol: Any, index: int) -> str:
    """Return a stable display/source name for one input molecule."""
    for prop in ("compound_id", "CompoundID", "compound_name", "_Name"):
        if mol.HasProp(prop):
            value = mol.GetProp(prop).strip()
            if value:
                return value
    return f"compound_{index:06d}"


def identity_source(mol: Any, compound_id: str) -> dict:
    """Match the structure-evaluator source-identity contract exactly."""
    identity = {"compound_id": str(compound_id)}
    for target, properties in _IDENTITY_ALIASES.items():
        for prop in properties:
            if mol.HasProp(prop):
                value = mol.GetProp(prop).strip()
                if value:
                    identity[target] = value
                    break
    return identity


def copy_identity_properties(source: Any, destination: Any) -> None:
    """Copy only source-lineage properties used by downstream pose identity."""
    properties = {
        prop
        for aliases in _IDENTITY_ALIASES.values()
        for prop in aliases
    }
    properties.update(LIGAND_LINEAGE_PROPERTIES)
    for prop in sorted(properties):
        if source.HasProp(prop):
            value = source.GetProp(prop).strip()
            if value:
                destination.SetProp(prop, value)


def annotate_canonical_pose(
    mol: Any,
    *,
    compound_id: str,
    pose_index: int,
    source_molecule_index: int,
    engine: str,
    engine_pose_index: int,
    engine_pose_rank: int,
) -> dict:
    """Attach a canonical pose identity shared by SDF and score records.

    The explicit ``source_identity`` property makes both downstream identity
    readers recompute the same identifier.  A failure is intentionally fatal:
    publishing an arbitrary fallback ID would permit score/coordinate mixing.
    """
    mol.SetProp("_Name", str(compound_id))
    mol.SetProp("compound_id", str(compound_id))
    mol.SetProp("docking_engine", str(engine))
    mol.SetIntProp("pose_index", int(pose_index))
    mol.SetIntProp("source_molecule_index", int(source_molecule_index))
    mol.SetIntProp("engine_pose_index", int(engine_pose_index))
    mol.SetIntProp("engine_pose_rank", int(engine_pose_rank))
    source_identity = identity_source(mol, str(compound_id))
    # Pose coordinates are serialized by SDWriter at the mol-block format's
    # precision.  Hashing the in-memory coordinates directly can therefore
    # produce an ID that no reader can recompute from the published SDF.  Build
    # the identity from an in-memory mol-block round trip, matching the bytes a
    # downstream RDKit reader will see.
    try:
        from rdkit import Chem

        serialized = Chem.MolFromMolBlock(
            Chem.MolToMolBlock(mol),
            removeHs=False,
            sanitize=False,
            strictParsing=False,
        )
    except Exception as exc:
        raise ValueError(
            f"could not round-trip docked pose through SDF mol block: {exc}"
        ) from exc
    if serialized is None:
        raise ValueError("could not round-trip docked pose through SDF mol block")
    identity = build_pose_identity(
        serialized,
        source_identity=source_identity,
        lineage_metadata={
            "operation": "docking_output",
            "source_pose_index": int(pose_index),
        },
    )
    chemical = identity["chemical_identity"]
    mol.SetProp("identity_schema_version", identity["schema_version"])
    mol.SetProp("pose_id", identity["pose_id"])
    mol.SetProp("chemical_id", chemical["chemical_id"])
    mol.SetProp("chemical_identity_id", chemical["identity_id"])
    mol.SetProp("source_identity_id", chemical["source_identity_id"])
    mol.SetProp(
        "source_identity",
        json.dumps(source_identity, sort_keys=True, separators=(",", ":")),
    )
    mol.SetProp("coordinate_signature", identity["coordinates"]["signature"])
    mol.SetProp("pose_identity_basis", POSE_IDENTITY_BASIS)
    mol.SetIntProp("coordinate_revision", 0)
    return {
        "pose_id": identity["pose_id"],
        "chemical_id": chemical["chemical_id"],
        "chemical_identity_id": chemical["identity_id"],
        "source_identity_id": chemical["source_identity_id"],
        "coordinate_signature": identity["coordinates"]["signature"],
        "pose_identity_basis": POSE_IDENTITY_BASIS,
    }

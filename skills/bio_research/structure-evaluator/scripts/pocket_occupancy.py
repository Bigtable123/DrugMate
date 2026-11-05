"""Pocket residue-contact fractions for docked poses and subpocket definitions."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import extract_resnum as _extract_resnum


def compute_pocket_occupancy(
    receptor_file: Path,
    mol,
    subpockets: list[dict],
    chain: str = "A",
    distance_cutoff: float = 4.0,
    *,
    subpocket_coords: list[list[list[float]]] | None = None,
) -> dict[str, float | None]:
    """Compute the heavy-atom contact fraction for each residue-property region.

    For each subpocket, returns the fraction of residue heavy atoms within
    distance_cutoff of any ligand heavy atom.  An unresolved coordinate set is
    returned as ``None`` rather than being conflated with measured zero contact.
    """
    ligand_coords = _get_ligand_coords(mol)
    if not ligand_coords:
        return {
            **{sp["id"]: None for sp in subpockets},
            "overall": None,
        }

    results = {}
    total_contacted = 0
    total_atoms = 0

    for index, sp in enumerate(subpockets):
        residues = sp.get("residues", [])
        sp_coords = (
            subpocket_coords[index] if subpocket_coords is not None
            else _get_pocket_atom_coords(receptor_file, residues, chain)
        )
        if not sp_coords:
            # A region that cannot be resolved against the receptor is missing
            # evidence, not a measured zero-contact region.
            results[sp["id"]] = None
            continue

        contacted = _count_contacted_atoms(sp_coords, ligand_coords, distance_cutoff)
        occupancy = contacted / len(sp_coords)
        # Preserve the measured fraction for thresholding.  Rounding here can
        # move a pose across the registered 0.05/0.20 gate boundaries (for
        # example 0.049 -> 0.05 or 0.196 -> 0.20).  Presentation layers may
        # round a copy, but the evaluator must consume the raw value.
        results[sp["id"]] = occupancy
        total_contacted += contacted
        total_atoms += len(sp_coords)

    # Keep a genuine measured zero when atoms were evaluated, but never turn an
    # empty denominator into zero: downstream gates must see it as NA.
    results["overall"] = (
        total_contacted / total_atoms if total_atoms > 0 else None
    )
    return results


def load_subpockets(pocket_features_file: Path) -> list[dict]:
    """Load subpocket definitions from pocket_features.json."""
    with open(pocket_features_file) as f:
        data = json.load(f)
    features = data if "key_subpockets" in data else data.get("pocket_features", {})
    return features.get("key_subpockets", [])


def _get_ligand_coords(mol) -> list[list[float]]:
    """Extract heavy atom coordinates from an RDKit mol."""
    from rdkit import Chem

    if mol is None:
        raise ValueError("cannot compute pocket occupancy for an unparsed ligand")
    mol_noH = Chem.RemoveHs(mol)
    conf = mol_noH.GetConformer()
    return [[conf.GetAtomPosition(i).x,
             conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z]
            for i in range(mol_noH.GetNumAtoms())]


def _count_contacted_atoms(
    pocket_coords: list[list[float]],
    ligand_coords: list[list[float]],
    cutoff: float,
) -> int:
    """Count pocket atoms within cutoff of any ligand atom."""
    import numpy as np

    pocket_arr = np.array(pocket_coords)
    ligand_arr = np.array(ligand_coords)
    contacted = 0
    for patom in pocket_arr:
        dists = np.linalg.norm(ligand_arr - patom, axis=1)
        if np.min(dists) <= cutoff:
            contacted += 1
    return contacted


def _get_pocket_atom_coords(
    protein_file: Path, pocket_residues: list[str], chain: str
) -> list[list[float]]:
    """Extract heavy atom coordinates for pocket residues from PDB."""
    return _select_pocket_coords(_read_pocket_atoms(protein_file, chain), pocket_residues)


def prepare_subpocket_coords(
    protein_file: Path, subpockets: list[dict], chain: str = "A",
) -> list[list[list[float]]]:
    """Read receptor coordinates once for all subpockets and all batch poses."""
    atoms = _read_pocket_atoms(protein_file, chain)
    return [_select_pocket_coords(atoms, sp.get("residues", [])) for sp in subpockets]


def _select_pocket_coords(atoms, pocket_residues):
    target_resnums = {_extract_resnum(res) for res in pocket_residues}
    target_resnums.discard(None)
    return [coords for resseq, coords in atoms if resseq in target_resnums]


def _read_pocket_atoms(protein_file: Path, chain: str):
    atoms = []

    with open(protein_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            line_chain = line[21].strip()
            if chain and line_chain != chain:
                continue
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((resseq, [x, y, z]))
            except (ValueError, IndexError):
                continue
    return atoms


# _extract_resnum is now imported from shared.residue_utils

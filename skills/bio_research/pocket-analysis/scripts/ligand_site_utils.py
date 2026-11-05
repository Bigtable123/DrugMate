"""Co-crystal ligand-defined binding-site utilities.

The observed ligand is the primary pocket definition when explicitly selected.
This module deliberately has no dependency on fpocket/P2Rank so it remains a
deterministic internal tool and can also supply a docking box directly.
"""
from __future__ import annotations

import sys
import math
from pathlib import Path
from typing import Iterable

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import (  # noqa: E402
    deduplicate_residue_records,
    legacy_residue_label,
    make_residue_record,
)


def build_box_site(protein_file: Path, docking_box: dict, chain: str | None = None) -> dict:
    """Characterize residues with a heavy atom inside the supplied box.

    This describes overlap with the selected chain, not evidence of a cavity
    or binding. Box geometry is never expanded to the residue envelope.
    """
    if chain is not None and not isinstance(chain, str):
        raise ValueError("chain must be a string or omitted")
    box_error = "provided box requires finite center, positive size, and angstrom units"
    if not isinstance(docking_box, dict) or docking_box.get("unit", "angstrom") != "angstrom":
        raise ValueError(box_error)
    for field in ("center", "size"):
        values = docking_box.get(field)
        if not isinstance(values, (list, tuple)) or len(values) != 3:
            raise ValueError(box_error)
        if any(
            isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or (field == "size" and value <= 0)
            for value in values
        ):
            raise ValueError(box_error)
    requested_chain = chain.strip() if chain is not None else None
    center, size = docking_box["center"], docking_box["size"]
    hits_by_chain = {}
    invalid_chains = set()
    with protein_file.open() as handle:
        for line in handle:
            if not line.startswith("ATOM") or _is_hydrogen_pdb_line(line):
                continue
            record = _residue_record_from_pdb_line(line)
            if record is None:
                continue
            hits = hits_by_chain.setdefault(record["chain"], [])
            coord = _pdb_coord(line)
            if coord is None or not all(math.isfinite(value) for value in coord):
                invalid_chains.add(record["chain"])
                continue
            if all(abs(coord[i] - center[i]) <= size[i] / 2 for i in range(3)):
                hits.append(record)
    counts = {key: len(value) for key, value in sorted(hits_by_chain.items())}
    overlap = [key for key, count in counts.items() if count]
    details = ", ".join(f"{key or '<blank>'}={count}" for key, count in counts.items()) or "none"
    if requested_chain:
        chain = requested_chain
        if chain in invalid_chains:
            raise ValueError(f"selected receptor chain {chain!r} contains invalid coordinates")
        if not counts.get(chain):
            raise ValueError(
                f"selected chain {chain!r} has no heavy atoms inside the provided box; "
                f"per-chain overlap: {details}. Review context.chain against the receptor "
                "and fixed box; the explicit chain was not changed."
            )
    else:
        if invalid_chains:
            raise ValueError("cannot infer chain because receptor chains contain invalid coordinates")
        if len(overlap) != 1 or not overlap[0]:
            raise ValueError(
                "provided box requires an explicit reviewed chain unless exactly one named "
                f"receptor chain overlaps; per-chain overlap: {details}. "
                "Review the target/assembly and set context.chain; do not assume A or choose "
                "the chain with most atoms."
            )
        chain = overlap[0]
    hits = hits_by_chain[chain]
    records = deduplicate_residue_records(hits)
    return {
        "candidate_id": "P001",
        "center": list(center), "box_center": list(center), "size": list(size),
        "residue_records": records,
        "residues": [legacy_residue_label(record) for record in records],
        "source": "user_provided_box",
        "selection_method": "selected_chain_heavy_atoms_inside_provided_box",
        "heavy_atom_count_inside_box": len(hits),
        "chain": chain,
        "chain_selection": {
            "requested_chain": requested_chain,
            "selected_chain": chain,
            "source": "explicit" if requested_chain else "unique_box_overlap",
            "heavy_atom_counts_by_chain": counts,
        },
        "detector_support": [], "support_count": 0,
    }


def build_ligand_site(
    protein_file: Path | str,
    ligand_file: Path | str,
    chain: str | Iterable[str] | None = "A",
    residue_cutoff: float = 5.0,
    box_padding: float = 4.0,
    min_box_size: float = 12.0,
) -> dict:
    """Build a ligand-defined pocket and docking box.

    Parameters
    ----------
    protein_file:
        Prepared receptor PDB.
    ligand_file:
        Co-crystal ligand in SDF or PDB format.  Only heavy atoms define the
        center, shell and bounds.
    chain:
        Protein chain filter.  A string selects one chain, an iterable selects
        several, and ``None``/``""`` includes every chain.
    residue_cutoff:
        A protein residue is in the shell when any heavy atom lies within this
        distance (angstrom) of any ligand heavy atom.
    box_padding:
        Padding added on both sides of each ligand-coordinate extent.
    min_box_size:
        Minimum length of each docking-box axis, preventing degenerate boxes
        for very small or nearly planar ligands.
    """
    protein_path = Path(protein_file)
    ligand_path = Path(ligand_file)
    if not protein_path.exists():
        raise FileNotFoundError(f"Protein file not found: {protein_path}")
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_path}")
    if residue_cutoff <= 0:
        raise ValueError("residue_cutoff must be positive")
    if box_padding < 0:
        raise ValueError("box_padding must be non-negative")
    if min_box_size < 0:
        raise ValueError("min_box_size must be non-negative")

    ligand_coords = load_ligand_heavy_atom_coords(ligand_path)
    if not ligand_coords:
        raise ValueError(f"No ligand heavy-atom coordinates found in {ligand_path}")

    box_center, size, bounds = box_from_coordinates(ligand_coords, padding=box_padding)
    size = [round(max(axis_size, float(min_box_size)), 4) for axis_size in size]
    centroid = [
        round(sum(point[axis] for point in ligand_coords) / len(ligand_coords), 4)
        for axis in range(3)
    ]
    chain_filter = _normalize_chain_filter(chain)
    cutoff_sq = float(residue_cutoff) ** 2

    residue_hits: list[dict] = []
    with protein_path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            record = _residue_record_from_pdb_line(line)
            if record is None:
                continue
            if chain_filter is not None and record["chain"] not in chain_filter:
                continue
            coord = _pdb_coord(line)
            if coord is None or _is_hydrogen_pdb_line(line):
                continue
            if any(_distance_sq(coord, ligand_coord) <= cutoff_sq for ligand_coord in ligand_coords):
                residue_hits.append(record)

    residue_records = deduplicate_residue_records(residue_hits)
    residues = [legacy_residue_label(record) for record in residue_records]
    return {
        "center": centroid,
        "centroid": centroid,
        "box_center": box_center,
        "size": size,
        "bounds": bounds,
        "residues": residues,
        "residue_records": residue_records,
        "source": "co_crystal_ligand",
        "selection_method": "ligand_heavy_atom_shell",
        "residue_cutoff": float(residue_cutoff),
        "box_padding": float(box_padding),
        "min_box_size": float(min_box_size),
        "ligand_heavy_atom_count": len(ligand_coords),
    }


def load_ligand_heavy_atom_coords(ligand_file: Path | str) -> list[list[float]]:
    """Read heavy-atom coordinates from an SDF or PDB ligand file."""
    path = Path(ligand_file)
    suffix = path.suffix.lower()
    if suffix in {".pdb", ".ent", ".pdbqt"}:
        return _load_pdb_heavy_atom_coords(path)
    if suffix in {".sdf", ".mol"}:
        return _load_sdf_with_rdkit(path)
    # A structure-prep artifact should normally be SDF, but accepting PDB-like
    # text makes the helper robust to extension-less temporary files.
    coords = _load_pdb_heavy_atom_coords(path)
    if coords:
        return coords
    return _load_sdf_with_rdkit(path)


def box_from_coordinates(
    coords: Iterable[Iterable[float]], padding: float = 4.0,
) -> tuple[list[float], list[float], dict[str, list[float]]]:
    """Return ``(center, size, bounds)`` for coordinates plus padding."""
    points = [[float(value) for value in point] for point in coords]
    if not points or any(len(point) != 3 for point in points):
        raise ValueError("coords must contain at least one xyz point")
    minima = [min(point[axis] for point in points) for axis in range(3)]
    maxima = [max(point[axis] for point in points) for axis in range(3)]
    center = [round((minima[axis] + maxima[axis]) / 2.0, 4) for axis in range(3)]
    size = [round(maxima[axis] - minima[axis] + 2.0 * padding, 4) for axis in range(3)]
    bounds = {
        "min": [round(value, 4) for value in minima],
        "max": [round(value, 4) for value in maxima],
    }
    return center, size, bounds


def _normalize_chain_filter(chain: str | Iterable[str] | None) -> set[str] | None:
    if chain is None or chain == "":
        return None
    if isinstance(chain, str):
        return {chain.strip()}
    chains = {str(item).strip() for item in chain}
    return chains or None


def _load_pdb_heavy_atom_coords(path: Path) -> list[list[float]]:
    coords: list[list[float]] = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")) or _is_hydrogen_pdb_line(line):
            continue
        coord = _pdb_coord(line)
        if coord is not None:
            coords.append(coord)
    return coords


def _load_sdf_with_rdkit(path: Path) -> list[list[float]]:
    from rdkit import Chem

    try:
        if path.suffix.lower() == ".mol":
            molecules = [Chem.MolFromMolFile(str(path), removeHs=False)]
        else:
            molecules = Chem.SDMolSupplier(str(path), removeHs=False)
        mol = next((item for item in molecules if item is not None), None)
        if mol is None or mol.GetNumConformers() == 0:
            return []
        conformer = mol.GetConformer()
        coords = []
        for index, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() == 1:
                continue
            position = conformer.GetAtomPosition(index)
            coords.append([float(position.x), float(position.y), float(position.z)])
        return coords
    except Exception:
        return []
def _residue_record_from_pdb_line(line: str) -> dict | None:
    try:
        return make_residue_record(
            line[17:20].strip(), int(line[22:26].strip()),
            line[21:22].strip(), line[26:27].strip(),
        )
    except (ValueError, IndexError):
        return None


def _pdb_coord(line: str) -> list[float] | None:
    try:
        return [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    except (ValueError, IndexError):
        return None


def _is_hydrogen_pdb_line(line: str) -> bool:
    atom_name = line[12:16].strip() if len(line) >= 16 else ""
    element = line[76:78].strip() if len(line) >= 78 else ""
    element = element or atom_name.lstrip("0123456789")[:1]
    return element.upper() in {"H", "D", "T"}


def _distance_sq(left: list[float], right: list[float]) -> float:
    return sum((left[index] - right[index]) ** 2 for index in range(3))

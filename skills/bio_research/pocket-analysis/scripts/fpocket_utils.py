"""fpocket utilities: invoke fpocket CLI and parse output files."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import (  # noqa: E402
    deduplicate_residue_records,
    legacy_residue_label,
    make_residue_record,
)


def run_fpocket(
    protein_file: Path,
    output_dir: Path,
    min_radius: float = 3.0,
    max_radius: float = 6.0,
    min_spheres: int = 30,
) -> Path:
    """Run fpocket on a PDB file.

    Args:
        protein_file: Path to input PDB
        output_dir: Directory for fpocket output
        min_radius: Minimum alpha sphere radius (-m)
        max_radius: Maximum alpha sphere radius (-M)
        min_spheres: Minimum number of alpha spheres per pocket (-i)

    Returns:
        Path to fpocket output directory
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    work_pdb = output_dir / protein_file.name
    shutil.copy2(protein_file, work_pdb)

    cmd = [
        "fpocket",
        "-f", str(work_pdb),
        "-m", str(min_radius),
        "-M", str(max_radius),
        "-i", str(min_spheres),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
    if result.returncode != 0:
        raise RuntimeError(f"fpocket failed: {result.stderr}")

    return output_dir


def parse_fpocket_output(output_dir: Path, protein_stem: str) -> dict:
    """Parse fpocket output files to extract pocket information.

    Args:
        output_dir: fpocket output directory
        protein_stem: stem of the protein file name (without extension)

    Returns:
        Dictionary with pocket detection results
    """
    fpocket_dir = output_dir / f"{protein_stem}_out"
    if not fpocket_dir.exists():
        for d in output_dir.iterdir():
            if d.is_dir() and d.name.endswith("_out"):
                fpocket_dir = d
                break

    if not fpocket_dir.exists():
        raise FileNotFoundError(f"fpocket output directory not found in {output_dir}")

    info_file = fpocket_dir / f"{protein_stem}_info.txt"
    if not info_file.exists():
        for f in fpocket_dir.glob("*_info.txt"):
            info_file = f
            break

    if not info_file.exists():
        raise FileNotFoundError(f"fpocket info file not found in {fpocket_dir}")
    pockets = _parse_info_file(info_file)

    pockets_dir = fpocket_dir / "pockets"
    if pockets_dir.exists():
        for pocket in pockets:
            rank = pocket["rank"]
            atm_file = pockets_dir / f"pocket{rank}_atm.pdb"
            if atm_file.exists():
                records = _extract_residue_records_from_pdb(atm_file)
                pocket["residue_records"] = records
                pocket["residues"] = [legacy_residue_label(record) for record in records]

            vert_file = pockets_dir / f"pocket{rank}_vert.pqr"
            if vert_file.exists():
                pocket["center"] = _compute_center_from_pqr(vert_file)
                pocket["support_file"] = str(vert_file.resolve())

    pockets.sort(key=lambda pocket: pocket["rank"])
    top_pocket = pockets[0] if pockets else None

    return {
        "num_pockets": len(pockets),
        "top_pocket": top_pocket,
        "all_pockets": pockets[:10],
    }


def _parse_info_file(info_file: Path) -> list[dict]:
    """Parse fpocket _info.txt file to extract pocket properties."""
    pockets = []
    current_pocket = None

    with open(info_file) as f:
        for line in f:
            line = line.strip()

            pocket_match = re.match(r"Pocket\s+(\d+)\s*:", line)
            if pocket_match:
                if current_pocket:
                    pockets.append(current_pocket)
                current_pocket = {"rank": int(pocket_match.group(1))}
                continue

            if current_pocket is None:
                continue

            if re.match(r"^Score\s*:", line, flags=re.IGNORECASE):
                val = _extract_float(line)
                if val is not None:
                    current_pocket["score"] = val
            elif re.match(r"^Druggability Score\s*:", line, flags=re.IGNORECASE):
                val = _extract_float(line)
                if val is not None:
                    current_pocket["druggability_score"] = val
            # fpocket 4 writes physical cavity volume as ``Volume : ...``;
            # older variants sometimes annotate it with ``MC``.  Anchor the
            # property name so ``Volume score:`` is not mistaken for volume.
            elif re.match(r"^Volume\s*:", line, flags=re.IGNORECASE):
                val = _extract_float(line)
                if val is not None:
                    current_pocket["volume"] = val
            elif "Number of Alpha Spheres" in line:
                val = _extract_float(line)
                if val is not None:
                    current_pocket["num_alpha_spheres"] = int(val)

    if current_pocket:
        pockets.append(current_pocket)

    return pockets


def _extract_float(line: str) -> Optional[float]:
    """Extract the last float from a line like 'Property : 3.14'."""
    match = re.search(r":\s*([-+]?\d*\.?\d+)", line)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            pass
    numbers = re.findall(r"[-+]?\d*\.?\d+", line)
    if numbers:
        try:
            return float(numbers[-1])
        except ValueError:
            pass
    return None


def _extract_residues_from_pdb(atm_file: Path) -> list[str]:
    """Extract legacy residue labels from a pocket atom PDB file."""
    return [
        legacy_residue_label(record)
        for record in _extract_residue_records_from_pdb(atm_file)
    ]


def _extract_residue_records_from_pdb(atm_file: Path) -> list[dict]:
    """Extract chain/insertion-code-safe residues from a pocket atom PDB."""
    residues = []
    with open(atm_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resseq = line[22:26].strip()
                chain = line[21].strip()
                icode = line[26].strip() if len(line) > 26 else ""
                if resname and resseq:
                    residues.append(make_residue_record(resname, resseq, chain, icode))
    return deduplicate_residue_records(residues)


def _compute_center_from_pqr(vert_file: Path) -> list[float] | None:
    """Compute geometric center from alpha sphere PQR file."""
    coords = []
    with open(vert_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue

    if not coords:
        return None

    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    # Preserve native geometry for the strict pocket matching thresholds.
    return [cx, cy, cz]

"""P2Rank utilities: invoke prank CLI and parse output CSV files."""
from __future__ import annotations

import csv
import math
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
    parse_residue_record,
)


def run_p2rank(protein_file: Path, output_dir: Path) -> Path:
    """Run P2Rank pocket prediction.

    Args:
        protein_file: Path to input PDB
        output_dir: Directory for P2Rank output

    Returns:
        Path to P2Rank output directory
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "prank",
        "predict",
        "-f", str(protein_file.resolve()),
        "-o", str(output_dir.resolve()),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"P2Rank failed: {result.stderr}")

    return output_dir


def parse_p2rank_output(output_dir: Path, protein_filename: str) -> dict:
    """Parse P2Rank predictions CSV.

    Args:
        output_dir: P2Rank output directory
        protein_filename: Original protein filename (e.g., 'protein.pdb')

    Returns:
        Dictionary with pocket prediction results
    """
    predictions_file = _find_predictions_file(output_dir, protein_filename)
    if predictions_file is None:
        raise FileNotFoundError(
            f"P2Rank predictions CSV not found in {output_dir} for {protein_filename}"
        )

    residues_file = _find_residues_file(output_dir, protein_filename)
    residue_details = {}
    if residues_file:
        residue_details = _parse_residues_csv(residues_file)

    pockets = _parse_predictions_csv(predictions_file, residue_details)
    top_pocket = pockets[0] if pockets else None

    return {
        "num_pockets": len(pockets),
        "top_pocket": top_pocket,
        "all_pockets": pockets[:10],
        "residue_details": residue_details if residue_details else None,
    }


def _find_predictions_file(output_dir: Path, protein_filename: str) -> Optional[Path]:
    """Find the predictions CSV in P2Rank output."""
    stem = Path(protein_filename).stem

    candidates = [
        output_dir / f"{protein_filename}_predictions.csv",
        output_dir / f"{stem}.pdb_predictions.csv",
    ]
    for c in candidates:
        if c.exists():
            return c

    for f in output_dir.rglob("*predictions.csv"):
        return f

    return None


def _find_residues_file(output_dir: Path, protein_filename: str) -> Optional[Path]:
    """Find the residues CSV in P2Rank output."""
    for f in output_dir.rglob("*residues.csv"):
        return f
    return None


def _parse_predictions_csv(
    csv_file: Path, residue_details: Optional[dict] = None
) -> list[dict]:
    """Parse P2Rank predictions CSV file."""
    pockets = []

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)

        headers = reader.fieldnames
        if not headers:
            return []

        clean_headers = [h.strip() for h in headers]

        for row in reader:
            clean_row = {k.strip(): v.strip() if v else "" for k, v in row.items()}

            raw_residues = clean_row.get(
                "residue_ids", clean_row.get("residues", "")
            )
            residue_records = _parse_residue_records(raw_residues, residue_details)
            center = _parse_center(clean_row)
            pocket = {
                "rank": _safe_int(clean_row.get("rank", "")),
                "name": clean_row.get("name", ""),
                "score": _safe_float(clean_row.get("score", "")),
                "probability": _safe_float(clean_row.get("prob", clean_row.get("probability", ""))),
                "center": center,
                "residues": [legacy_residue_label(record) for record in residue_records],
                "residue_records": residue_records,
            }
            if clean_row.get("surf_atom_ids"):
                serials = [int(value) for value in clean_row["surf_atom_ids"].split()]
                if len(serials) != int(clean_row["surf_atoms"]):
                    raise ValueError("P2Rank surface atom count mismatch")
                pocket["support_atom_serials"] = serials
            pockets.append(pocket)

    pockets.sort(key=lambda p: p.get("rank") or 999)
    return pockets


def _parse_residues_csv(csv_file: Path) -> dict:
    """Parse P2Rank residues CSV for per-residue pocket assignments.

    The residues CSV uses `residue_label` = chain-prefixed id (e.g. "A_LEU83")
    and `residue_name` = "LEU". The predictions CSV `residue_ids` column contains
    the same chain-prefixed tokens (e.g. "A_LEU83"), so the map must be keyed by
    that token verbatim. `residue_id` is normalized to the fpocket-style "LEU83"
    (name+number) so downstream consensus set math compares the same residue as
    one member rather than "LEU83" vs bare "83".
    """
    residue_map = {}

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            clean_row = {k.strip(): v.strip() if v else "" for k, v in row.items()}
            chain = clean_row.get("chain", "")
            res_label = clean_row.get("residue_label", clean_row.get("residue", ""))
            res_name = clean_row.get("residue_name", "")
            res_number = clean_row.get("residue_number", "")
            pocket = clean_row.get("pocket", clean_row.get("pocket_name", ""))
            score = _safe_float(clean_row.get("score", clean_row.get("pocket_score", "0")))
            if not res_label:
                continue
            # P2Rank versions differ here: predictions may refer to bare labels
            # ("10") or chain-prefixed labels ("A_10"). Index both forms.
            source_ids = [res_label]
            if chain and "_" not in res_label:
                source_ids.append(f"{chain}_{res_label}")
            # Prefer name+number ("LEU83"); fall back to stripping the chain
            # prefix off the label when name/number columns are absent.
            label_without_chain = _strip_chain(res_label)
            record = parse_residue_record(
                label_without_chain,
                default_chain=chain,
                default_resname=res_name,
            )
            if record is None and res_number:
                record = make_residue_record(res_name, res_number, chain)
            if record is None:
                continue
            residue_id = legacy_residue_label(record)
            detail = {
                "chain": chain,
                "residue_name": res_name,
                "residue_label": res_label,
                "residue_id": residue_id,
                "residue_record": record,
                "pocket": pocket,
                "score": score,
            }
            for source_id in source_ids:
                residue_map[source_id] = detail

    return residue_map


def _strip_chain(label: str) -> str:
    """Drop a leading chain prefix from a residue label: 'A_LEU83' -> 'LEU83'."""
    return label.split("_", 1)[1] if "_" in label else label


def _parse_residue_ids(
    residue_str: str, residue_details: Optional[dict] = None
) -> list[str]:
    """Parse residue IDs from old (A_LEU83) and new (A_83) P2Rank formats."""
    return [
        legacy_residue_label(record)
        for record in _parse_residue_records(residue_str, residue_details)
    ]


def _parse_residue_records(
    residue_str: str, residue_details: Optional[dict] = None
) -> list[dict]:
    """Parse P2Rank residue ids without discarding chain or insertion code."""
    if not residue_str:
        return []

    residues: list[dict] = []
    parts = residue_str.replace('"', "").replace("'", "").split()
    for part in parts:
        part = part.strip().rstrip(",")
        detail = (residue_details or {}).get(part, {})
        record = detail.get("residue_record")
        if record is None:
            record = parse_residue_record(part)
        if record is not None:
            residues.append(record)
    return deduplicate_residue_records(residues)


def _safe_float(value: str) -> float:
    """Safely convert string to float."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0


def _parse_center(row: dict[str, str]) -> list[float] | None:
    """Return a real xyz center; never manufacture an origin for missing data."""
    raw = [row.get("center_x", ""), row.get("center_y", ""), row.get("center_z", "")]
    if any(value in (None, "") for value in raw):
        return None
    try:
        center = [float(value) for value in raw]
    except (TypeError, ValueError):
        return None
    return center if all(math.isfinite(value) for value in center) else None


def _safe_int(value: str) -> int:
    """Safely convert string to int."""
    try:
        return int(value)
    except (ValueError, TypeError):
        return 0

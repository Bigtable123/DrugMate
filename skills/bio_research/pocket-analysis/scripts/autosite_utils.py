"""AutoSite native cluster adapter, matching the supplied benchmark CSV parser."""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
import shutil
import tempfile

from tool_runtime import run_tool


def parse_autosite_output(raw_dir: Path) -> dict:
    """Preserve native cluster IDs/ranks and use exported point centroids."""
    raw_dir = Path(raw_dir).resolve()
    summary = raw_dir / "receptor_summary.csv"
    pockets, seen = [], set()
    with summary.open(newline="") as handle:
        reader = csv.reader(handle)
        if not next(reader, None):
            raise ValueError("AutoSite summary has no header")
        for row in reader:
            if not row or not any(cell.strip() for cell in row):
                continue
            rank = int(row[0])
            if rank < 1 or rank in seen or len(row) < 2:
                raise ValueError("AutoSite has invalid or duplicate native cluster ranks")
            seen.add(rank)
            properties = [float(value) for value in row[1:]]
            if not all(math.isfinite(value) for value in properties):
                raise ValueError("AutoSite summary contains nonfinite values")
            points_file = raw_dir / f"receptor_cl_{rank:03d}.pdb"
            points = []
            for line in points_file.read_text().splitlines():
                if line.startswith(("ATOM  ", "HETATM")):
                    point = [float(line[start:start + 8]) for start in (30, 38, 46)]
                    if not all(math.isfinite(value) for value in point):
                        raise ValueError("AutoSite cluster contains nonfinite coordinates")
                    points.append(point)
            if not points:
                raise ValueError(f"AutoSite native cluster {rank} is empty")
            pockets.append({
                "rank": rank, "pocket_id": str(rank),
                "center": [math.fsum(p[axis] for p in points) / len(points) for axis in range(3)],
                "center_definition": "centroid of exported AutoSite cluster points",
                "score": properties[-1], "native_properties": properties,
                "support_file": str(points_file),
                "support_sha256": hashlib.sha256(points_file.read_bytes()).hexdigest(),
                "support_point_count": len(points),
            })
    if 1 not in seen:
        raise ValueError("AutoSite output has no native rank-1 cluster")
    pockets.sort(key=lambda pocket: pocket["rank"])
    return {
        "success": True, "all_pockets": pockets, "top_pocket": pockets[0],
        "native_summary": str(summary),
        "native_summary_sha256": hashlib.sha256(summary.read_bytes()).hexdigest(),
    }


def run_autosite(protein_file: Path, output_dir: Path, *, runtime: dict, timeout: int = 1200) -> dict:
    """Prepare a private receptor copy and run native defaults in ADFRsuite."""
    protein_file = Path(protein_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix="attempt_", dir=output_dir))
    receptor = attempt / "receptor.pdb"
    shutil.copyfile(protein_file, receptor)
    commands = runtime.get("commands") or {}
    records = []
    for label, command, args in (
        ("prepare", commands.get("prepare_receptor", "prepare_receptor"),
         ["-r", "receptor.pdb", "-o", "receptor.pdbqt", "-A", "hydrogens"]),
        ("autosite", commands.get("engine", runtime.get("command", "autosite")),
         ["-r", "receptor.pdbqt", "-o", "raw"]),
    ):
        result = run_tool(runtime, args, command=command, cwd=attempt, timeout=timeout, check=False)
        for stream in ("stdout", "stderr"):
            (attempt / f"{label}.{stream}.log").write_text(result.get(stream, ""))
        records.append({key: value for key, value in result.items() if key not in {"stdout", "stderr"}})
        if not result.get("success"):
            raise RuntimeError(f"AutoSite {label} failed; inspect {attempt}")
        if label == "prepare":
            prepared = attempt / "receptor.pdbqt"
            if not prepared.is_file() or not any(
                line.startswith(("ATOM  ", "HETATM")) for line in prepared.read_text().splitlines()
            ):
                raise ValueError("AutoSite preparation produced no receptor atoms")
    result = parse_autosite_output(attempt / "raw")
    result.update({
        "input": {"protein_file": str(protein_file),
                  "protein_sha256": hashlib.sha256(receptor.read_bytes()).hexdigest()},
        "commands": records, "output_dir": str(attempt),
        "parameters": {"detector_parameters": "AutoSite native defaults", "timeout_seconds": timeout},
    })
    return result

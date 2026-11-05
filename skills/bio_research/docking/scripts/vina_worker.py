#!/usr/bin/env python3
"""Run prepared PDBQT ligands through the AutoDock Vina Python API.

This worker intentionally imports only Vina and the Python standard library so
it can run in the small dedicated ``vina`` Conda environment.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
from pathlib import Path


def run_manifest(manifest: dict) -> dict:
    """Dock every prepared ligand and retain per-pose native energy terms."""
    from vina import Vina

    receptor = Path(manifest["receptor_pdbqt"])
    if not receptor.is_file():
        raise FileNotFoundError(f"receptor PDBQT not found: {receptor}")

    params = dict(manifest.get("params") or {})
    scoring_function = str(params.get("scoring_function", "vina"))
    if scoring_function not in {"vina", "vinardo"}:
        raise ValueError("scoring_function must be vina or vinardo")

    box = manifest.get("docking_box") or {}
    center = _finite_triplet(box.get("center"), "docking_box.center", positive=False)
    size = _finite_triplet(box.get("size"), "docking_box.size", positive=True)
    output_dir = Path(manifest["native_output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    def new_engine():
        instance = Vina(
            sf_name=scoring_function,
            cpu=int(params.get("cpu", 1)),
            seed=int(params.get("seed", 42)),
            no_refine=bool(params.get("no_refine", False)),
            verbosity=int(params.get("verbosity", 0)),
        )
        instance.set_receptor(str(receptor))
        # Computing maps before loading any ligand creates the complete atom-type
        # map set, which can then be reused safely across a heterogeneous library.
        instance.compute_vina_maps(
            center=center,
            box_size=size,
            spacing=float(params.get("spacing", 0.375)),
            force_even_voxels=bool(params.get("force_even_voxels", False)),
        )
        return instance

    engine = new_engine()

    entries: list[dict] = []
    ligands = manifest.get("ligands") or []
    for ligand_index, ligand in enumerate(ligands):
        result = dict(ligand)
        input_pdbqt = Path(ligand["input_pdbqt"])
        source_index = int(ligand["source_molecule_index"])
        variant_index = int(ligand.get("preparation_variant_index", 0))
        output_pdbqt = output_dir / (
            f"ligand_{source_index:06d}_{variant_index:03d}_out.pdbqt"
        )
        try:
            engine.set_ligand_from_file(str(input_pdbqt))
            engine.dock(
                exhaustiveness=int(params.get("exhaustiveness", 16)),
                n_poses=int(params.get("num_modes", 5)),
                min_rmsd=float(params.get("min_rmsd", 1.0)),
                max_evals=int(params.get("max_evals", 0)),
            )
            energies = engine.energies(
                n_poses=int(params.get("num_modes", 5)),
                energy_range=float(params.get("energy_range", 3.0)),
            )
            energy_rows = [
                [float(value) for value in row]
                for row in energies.tolist()
            ]
            if not energy_rows or not all(
                row and all(math.isfinite(value) for value in row)
                for row in energy_rows
            ):
                raise RuntimeError("Vina returned no finite pose energies")
            engine.write_poses(
                str(output_pdbqt),
                n_poses=int(params.get("num_modes", 5)),
                energy_range=float(params.get("energy_range", 3.0)),
                overwrite=True,
            )
            if not output_pdbqt.is_file() or output_pdbqt.stat().st_size == 0:
                raise RuntimeError("Vina did not write a pose PDBQT")
            result.update({
                "success": True,
                "output_pdbqt": str(output_pdbqt),
                "energies": energy_rows,
            })
        except Exception as exc:
            result.update({
                "success": False,
                "error": str(exc),
                "output_pdbqt": str(output_pdbqt),
                "energies": [],
            })
            # Some Vina failures leave the object in an uncertain ligand state.
            # Rebuild it before attempting the next molecule, while retaining
            # the first failure in the result manifest.
            if ligand_index + 1 < len(ligands):
                engine = new_engine()
        entries.append(result)
        status_path = manifest.get("status_path")
        if status_path:
            _write_json_atomic(
                Path(status_path),
                {
                    "success": False,
                    "complete": False,
                    "engine": "autodock_vina",
                    "expected_preparations": len(ligands),
                    "completed_preparations": len(entries),
                    "ligands": entries,
                },
            )

    successful = sum(1 for entry in entries if entry.get("success"))
    return {
        "success": bool(entries) and successful == len(entries),
        "engine": "autodock_vina",
        "engine_version": _vina_version(),
        "scoring_function": scoring_function,
        "docking_box": {"center": center, "size": size},
        "expected_molecules": int(manifest.get("expected_molecules", len(entries))),
        "expected_preparations": len(entries),
        "successful_preparations": successful,
        "ligands": entries,
        "params": params,
    }


def _vina_version() -> str:
    try:
        return importlib.metadata.version("vina")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _finite_triplet(value, name: str, *, positive: bool) -> list[float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly three numbers")
    triplet = [float(item) for item in value]
    if not all(math.isfinite(item) for item in triplet):
        raise ValueError(f"{name} contains a non-finite number")
    if positive and not all(item > 0 for item in triplet):
        raise ValueError(f"{name} values must be > 0")
    return triplet


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    partial.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    result_path = Path(args.result)
    try:
        manifest = json.loads(Path(args.manifest).read_text())
        result = run_manifest(manifest)
    except Exception as exc:
        result = {"success": False, "error": str(exc), "ligands": []}
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    raise SystemExit(0 if result.get("success") else 1)


if __name__ == "__main__":
    main()

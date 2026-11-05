#!/usr/bin/env python3
"""Prepare and export AutoDock Vina ligands in the validated Meeko runtime."""
from __future__ import annotations

import argparse
import io
import json
import math
from pathlib import Path


_RESERVED_PROPERTIES = {
    "meeko", "vina_affinity", "vina_inter", "vina_intra", "vina_torsions",
    "vina_intra_best_pose", "docking_engine", "score_direction", "pose_index",
    "source_molecule_index", "preparation_variant_index", "engine_pose_index",
    "engine_pose_rank", "heavy_atom_count", "ligand_efficiency",
    "score_gap_1_2",
}


def prepare_ligands(input_sdf: Path, output_dir: Path) -> dict:
    """Convert a multi-molecule 3D SDF into mapped per-ligand PDBQT files."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    output_dir.mkdir(parents=True, exist_ok=True)
    supplier = Chem.SDMolSupplier(
        str(input_sdf), removeHs=False, strictParsing=False,
    )
    ligands: list[dict] = []
    rejects: list[dict] = []
    prepared_sources: set[int] = set()
    total_records = 0

    for source_index, molecule in enumerate(supplier):
        total_records += 1
        if molecule is None:
            rejects.append({
                "source_molecule_index": source_index,
                "reason": "RDKit could not parse the SDF record",
            })
            continue
        compound_id = _compound_id(molecule, source_index)
        try:
            _validate_3d(molecule)
            source_properties = {
                name: molecule.GetProp(name)
                for name in molecule.GetPropNames(
                    includePrivate=False, includeComputed=False,
                )
            }
            molecule_h = Chem.AddHs(molecule, addCoords=True)
            setups = MoleculePreparation(add_index_map=True).prepare(molecule_h)
            if not setups:
                raise ValueError("Meeko produced no ligand preparation")
            for variant_index, setup in enumerate(setups):
                pdbqt_text, success, error = PDBQTWriterLegacy.write_string(
                    setup, add_index_map=True,
                )
                if not success:
                    raise ValueError(error or "Meeko PDBQT writer failed")
                pdbqt_path = output_dir / (
                    f"ligand_{source_index:06d}_{variant_index:03d}.pdbqt"
                )
                pdbqt_path.write_text(pdbqt_text)
                ligands.append({
                    "source_molecule_index": source_index,
                    "preparation_variant_index": variant_index,
                    "compound_id": compound_id,
                    "input_pdbqt": str(pdbqt_path),
                    "source_properties": source_properties,
                })
            prepared_sources.add(source_index)
        except Exception as exc:
            rejects.append({
                "source_molecule_index": source_index,
                "compound_id": compound_id,
                "reason": str(exc),
            })

    return {
        "success": (
            total_records > 0
            and not rejects
            and len(prepared_sources) == total_records
            and bool(ligands)
        ),
        "input_sdf": str(input_sdf),
        "total_records": total_records,
        "prepared_molecules": len(prepared_sources),
        "prepared_variants": len(ligands),
        "ligands": ligands,
        "rejects": rejects,
    }


def export_poses(run_result: dict, output_sdf: Path) -> dict:
    """Convert all Vina PDBQT poses to one lineage-preserving canonical SDF."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    partial_sdf = output_sdf.with_suffix(".partial.sdf")
    writer = Chem.SDWriter(str(partial_sdf))
    scores: list[dict] = []
    failures: list[dict] = []
    successful_molecules: set[int] = set()
    successful_preparations: set[tuple[int, int]] = set()
    try:
        for entry in run_result.get("ligands") or []:
            source_index = int(entry.get("source_molecule_index", -1))
            variant_index = int(entry.get("preparation_variant_index", 0))
            compound_id = str(entry.get("compound_id") or f"MOL_{source_index:06d}")
            if not entry.get("success"):
                failures.append({
                    "source_molecule_index": source_index,
                    "compound_id": compound_id,
                    "reason": entry.get("error", "Vina docking failed"),
                })
                continue
            try:
                pdbqt_molecule = PDBQTMolecule.from_file(entry["output_pdbqt"])
                sdf_text, conversion_failures = (
                    RDKitMolCreate.write_sd_string(pdbqt_molecule)
                )
                if conversion_failures:
                    raise ValueError(
                        "Meeko could not reconstruct molecule indices "
                        f"{conversion_failures}"
                    )
                molecules = list(Chem.ForwardSDMolSupplier(
                    io.BytesIO(sdf_text.encode()),
                    removeHs=False,
                    strictParsing=False,
                ))
                energies = entry.get("energies") or []
                if len(molecules) != len(energies):
                    raise ValueError(
                        "Meeko pose count does not match Vina energy count "
                        f"({len(molecules)} != {len(energies)})"
                    )
                affinity_values = [
                    float(row[0])
                    for row in energies
                    if isinstance(row, (list, tuple)) and row
                ]
                score_gap = None
                if (
                    len(affinity_values) >= 2
                    and all(math.isfinite(value) for value in affinity_values)
                ):
                    ordered_affinities = sorted(affinity_values)
                    score_gap = ordered_affinities[1] - ordered_affinities[0]
                file_pose_count = 0
                for engine_pose_index, (molecule, energy_row) in enumerate(
                    zip(molecules, energies, strict=True)
                ):
                    if molecule is None:
                        raise ValueError(
                            f"Meeko returned an empty pose at index {engine_pose_index}"
                        )
                    values = [float(value) for value in energy_row]
                    if not values or not all(math.isfinite(value) for value in values):
                        raise ValueError(
                            f"non-finite Vina energy row at pose {engine_pose_index}"
                        )
                    affinity = values[0]
                    heavy_atom_count = int(molecule.GetNumHeavyAtoms())
                    if heavy_atom_count <= 0:
                        raise ValueError(
                            f"pose {engine_pose_index} has no heavy atoms"
                        )
                    ligand_efficiency = -affinity / heavy_atom_count
                    engine_pose_rank = engine_pose_index + 1
                    meeko_metadata = _meeko_metadata(molecule)
                    exported_affinity = meeko_metadata.get("free_energy")
                    if exported_affinity is not None and not math.isclose(
                        float(exported_affinity), affinity, abs_tol=0.02,
                    ):
                        raise ValueError(
                            "Meeko/Vina affinity mismatch at pose "
                            f"{engine_pose_index}: {exported_affinity} != {affinity}"
                        )

                    for name, value in (entry.get("source_properties") or {}).items():
                        target = f"source_{name}" if name in _RESERVED_PROPERTIES else name
                        molecule.SetProp(str(target), str(value))
                    pose_index = len(scores)
                    molecule.SetProp("_Name", compound_id)
                    molecule.SetProp("docking_engine", "autodock_vina")
                    molecule.SetProp("score_direction", "lower")
                    molecule.SetDoubleProp("vina_affinity", affinity)
                    molecule.SetIntProp("heavy_atom_count", heavy_atom_count)
                    molecule.SetDoubleProp(
                        "ligand_efficiency", ligand_efficiency,
                    )
                    molecule.SetIntProp("engine_pose_rank", engine_pose_rank)
                    if score_gap is not None:
                        molecule.SetDoubleProp("score_gap_1_2", score_gap)
                    molecule.SetIntProp("pose_index", pose_index)
                    molecule.SetIntProp("source_molecule_index", source_index)
                    molecule.SetIntProp(
                        "preparation_variant_index", variant_index,
                    )
                    molecule.SetIntProp("engine_pose_index", engine_pose_index)
                    component_names = (
                        "vina_inter", "vina_intra", "vina_torsions",
                        "vina_intra_best_pose",
                    )
                    for component_name, value in zip(
                        component_names, values[1:], strict=False,
                    ):
                        molecule.SetDoubleProp(component_name, value)
                    writer.write(molecule)

                    score = {
                        "pose_index": pose_index,
                        "source_molecule_index": source_index,
                        "preparation_variant_index": variant_index,
                        "engine_pose_index": engine_pose_index,
                        "engine_pose_rank": engine_pose_rank,
                        "compound_id": compound_id,
                        "mol_name": compound_id,
                        "engine": "autodock_vina",
                        "vina_affinity": affinity,
                        "heavy_atom_count": heavy_atom_count,
                        "ligand_efficiency": ligand_efficiency,
                    }
                    if score_gap is not None:
                        score["score_gap_1_2"] = score_gap
                    for component_name, value in zip(
                        component_names, values[1:], strict=False,
                    ):
                        score[component_name] = value
                    scores.append(score)
                    file_pose_count += 1
                if file_pose_count:
                    successful_molecules.add(source_index)
                    successful_preparations.add((source_index, variant_index))
            except Exception as exc:
                failures.append({
                    "source_molecule_index": source_index,
                    "compound_id": compound_id,
                    "reason": str(exc),
                })
    finally:
        writer.close()

    expected = int(run_result.get("expected_molecules", 0))
    expected_preparations = int(
        run_result.get("expected_preparations", expected)
    )
    diagnostics = {
        "expected_molecules": expected,
        "successful_molecules": len(successful_molecules),
        "expected_preparations": expected_preparations,
        "successful_preparations": len(successful_preparations),
        "failures": failures,
    }
    if not scores:
        return {
            "success": False,
            "error": (
                "AutoDock Vina exported no scored pose: "
                f"0/{expected} molecules and "
                f"0/{expected_preparations} preparations"
            ),
            "scores": [],
            "partial_docked_sdf": str(partial_sdf),
            "diagnostics": diagnostics,
        }
    # A molecule or preparation variant Vina could not export is that candidate's
    # missing value, not the batch's failure.  Returning here discarded every pose
    # the run did produce, so one unexportable ligand ended docking for the whole
    # cohort -- while the executable policy scopes a single tool fault to the single
    # candidate and the executable policy removes a missing value from numerator and
    # denominator instead.  The shortfall stays itemized in ``diagnostics``
    # (``failures`` already names each one) and is declared in
    # ``operational_errors`` for the failure registry.
    diagnostics["complete"] = (
        not failures
        and len(successful_molecules) == expected
        and len(successful_preparations) == expected_preparations
    )

    partial_sdf.replace(output_sdf)
    merged = {
        "success": True,
        "scores": scores,
        "diagnostics": diagnostics,
        "docked_sdf": str(output_sdf),
    }
    if not diagnostics["complete"]:
        merged["operational_errors"] = [{
            "tool": "autodock_vina",
            "reason": (
                "AutoDock Vina output is incomplete: "
                f"{len(successful_molecules)}/{expected} molecules and "
                f"{len(successful_preparations)}/{expected_preparations} "
                "preparations exported"
            ),
        }]
    return merged


def _compound_id(molecule, source_index: int) -> str:
    if molecule.HasProp("_Name") and molecule.GetProp("_Name").strip():
        return molecule.GetProp("_Name").strip()
    for prop in ("compound_id", "mol_name", "ID", "id"):
        if molecule.HasProp(prop) and molecule.GetProp(prop).strip():
            return molecule.GetProp(prop).strip()
    return f"MOL_{source_index:06d}"


def _validate_3d(molecule) -> None:
    if molecule.GetNumConformers() == 0:
        raise ValueError("molecule has no conformer")
    conformer = molecule.GetConformer()
    if not conformer.Is3D():
        raise ValueError("molecule conformer is not marked as 3D")
    for atom_index in range(molecule.GetNumAtoms()):
        position = conformer.GetAtomPosition(atom_index)
        if not all(math.isfinite(value) for value in (position.x, position.y, position.z)):
            raise ValueError("molecule contains non-finite coordinates")


def _meeko_metadata(molecule) -> dict:
    if not molecule.HasProp("meeko"):
        return {}
    try:
        value = json.loads(molecule.GetProp("meeko"))
    except (TypeError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--input-sdf", required=True)
    prepare_parser.add_argument("--output-dir", required=True)
    prepare_parser.add_argument("--result", required=True)

    export_parser = subparsers.add_parser("export")
    export_parser.add_argument("--run-result", required=True)
    export_parser.add_argument("--output-sdf", required=True)
    export_parser.add_argument("--result", required=True)

    args = parser.parse_args()
    result_path = Path(args.result)
    try:
        if args.command == "prepare":
            result = prepare_ligands(
                Path(args.input_sdf), Path(args.output_dir),
            )
        else:
            run_result = json.loads(Path(args.run_result).read_text())
            result = export_poses(run_result, Path(args.output_sdf))
    except Exception as exc:
        result = {"success": False, "error": str(exc)}
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    raise SystemExit(0 if result.get("success") else 1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Isolated OpenMM simulation and MDAnalysis trajectory worker."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _simulate(payload: dict) -> dict:
    import openmm
    from openmm import app, unit

    prmtop_path = Path(payload["prmtop"]).resolve()
    inpcrd_path = Path(payload["inpcrd"]).resolve()
    output_dir = Path(payload["output_dir"]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in (prmtop_path, inpcrd_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    prmtop = app.AmberPrmtopFile(str(prmtop_path))
    inpcrd = app.AmberInpcrdFile(str(inpcrd_path))
    if inpcrd.boxVectors is None:
        raise RuntimeError("short MD requires periodic box vectors in Amber inpcrd")

    temperature = float(payload["temperature_kelvin"]) * unit.kelvin
    pressure = float(payload["pressure_bar"]) * unit.bar
    timestep_fs = float(payload["timestep_femtoseconds"])
    if not 0.0 < timestep_fs <= 2.0:
        raise ValueError(
            "short MD without hydrogen-mass repartitioning requires a "
            "timestep in (0, 2] fs"
        )
    friction = float(payload["friction_per_picosecond"]) / unit.picosecond
    timestep = timestep_fs * unit.femtoseconds
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    seed = int(payload.get("seed", 42))
    barostat = openmm.MonteCarloBarostat(pressure, temperature, 25)
    barostat.setRandomNumberSeed(seed)
    system.addForce(barostat)
    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, timestep)
    integrator.setRandomNumberSeed(seed)

    platform_name = str(payload.get("platform", "CUDA"))
    platform = openmm.Platform.getPlatformByName(platform_name)
    properties = {}
    if platform_name in {"CUDA", "OpenCL"}:
        properties["Precision"] = str(payload.get("precision", "mixed"))
        properties["DeviceIndex"] = str(payload.get("device_index", 0))
    simulation = app.Simulation(
        prmtop.topology, system, integrator, platform, properties
    )
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.context.setPositions(inpcrd.positions)

    starting_pdb = output_dir / "starting_docked_structure.pdb"
    with starting_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(
            prmtop.topology, inpcrd.positions, handle, keepIds=True
        )

    max_iterations = int(payload.get("minimization_max_iterations", 500))
    if max_iterations > 0:
        simulation.minimizeEnergy(maxIterations=max_iterations)
    simulation.context.setVelocitiesToTemperature(temperature, seed)

    equilibration_ps = float(payload.get("equilibration_picoseconds", 20.0))
    production_ps = float(payload.get("production_picoseconds", 100.0))
    report_ps = float(payload.get("report_interval_picoseconds", 2.0))
    steps_per_ps = 1000.0 / timestep_fs
    equilibration_steps = max(0, int(round(equilibration_ps * steps_per_ps)))
    production_steps = max(1, int(round(production_ps * steps_per_ps)))
    report_interval = max(1, int(round(report_ps * steps_per_ps)))
    if report_interval > production_steps:
        report_interval = production_steps

    if equilibration_steps:
        simulation.step(equilibration_steps)

    # The COM-drift and ligand-RMSD definitions use the equilibrated production
    # baseline, not the pre-minimization docked coordinates.  Keep both files:
    # the original remains useful provenance, while this snapshot is the metric
    # reference required by the executable policy.
    production_reference = output_dir / "production_reference.pdb"
    production_start_state = simulation.context.getState(getPositions=True)
    with production_reference.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(
            prmtop.topology,
            production_start_state.getPositions(),
            handle,
            keepIds=True,
        )

    # Reporter scheduling is based on the absolute simulation step. Reset the
    # production clock so arbitrary equilibration/report combinations yield an
    # exact, auditable production-frame count and a clear time origin.
    simulation.currentStep = 0
    simulation.context.setTime(0.0 * unit.picoseconds)

    trajectory = output_dir / "trajectory.dcd"
    state_csv = output_dir / "state.csv"
    simulation.reporters.append(
        app.DCDReporter(str(trajectory), report_interval, enforcePeriodicBox=False)
    )
    simulation.reporters.append(
        app.StateDataReporter(
            str(state_csv),
            report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            density=True,
            volume=True,
            separator=",",
        )
    )
    simulation.step(production_steps)

    checkpoint = output_dir / "final.chk"
    final_state = output_dir / "final_state.xml"
    final_pdb = output_dir / "final.pdb"
    simulation.saveCheckpoint(str(checkpoint))
    simulation.saveState(str(final_state))
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    potential_energy = float(
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    kinetic_energy = float(
        state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    if not all(
        math.isfinite(value) for value in (potential_energy, kinetic_energy)
    ):
        raise RuntimeError("OpenMM final energy is non-finite")
    with final_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(
            prmtop.topology, state.getPositions(), handle, keepIds=True
        )
    output_paths = (
        starting_pdb, production_reference, trajectory, state_csv, checkpoint,
        final_state, final_pdb,
    )
    for path in output_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"OpenMM did not create a non-empty output: {path}")

    expected_frames = production_steps // report_interval
    with state_csv.open(newline="", encoding="utf-8") as handle:
        state_rows = list(csv.DictReader(handle))
    if len(state_rows) != expected_frames or expected_frames < 1:
        raise RuntimeError(
            "OpenMM reporter frame count is inconsistent: "
            f"expected {expected_frames}, found {len(state_rows)}"
        )
    for row_index, row in enumerate(state_rows):
        for label, raw_value in row.items():
            try:
                value = float(raw_value)
            except (TypeError, ValueError) as exc:
                raise RuntimeError(
                    f"OpenMM state row {row_index} field {label!r} is non-numeric"
                ) from exc
            if not math.isfinite(value):
                raise RuntimeError(
                    f"OpenMM state row {row_index} field {label!r} is non-finite"
                )

    properties_reported = {
        name: platform.getPropertyValue(simulation.context, name)
        for name in platform.getPropertyNames()
    }
    return {
        "success": True,
        "action": "simulate",
        "openmm_version": str(openmm.__version__),
        "platform": platform_name,
        "platform_properties": properties_reported,
        "atom_count": int(system.getNumParticles()),
        "periodic": True,
        "equilibration_steps": equilibration_steps,
        "production_steps": production_steps,
        "report_interval_steps": report_interval,
        "frame_count_expected": expected_frames,
        "final_energy_kilojoule_per_mole": {
            "potential": potential_energy,
            "kinetic": kinetic_energy,
        },
        "files": {
            "starting_docked_pdb": str(starting_pdb),
            "production_reference_pdb": str(production_reference),
            "trajectory_dcd": str(trajectory),
            "state_csv": str(state_csv),
            "checkpoint": str(checkpoint),
            "final_state_xml": str(final_state),
            "final_pdb": str(final_pdb),
        },
        "file_sha256": {path.name: _sha256(path) for path in output_paths},
    }


def _statistics(values: list[float]) -> dict:
    import numpy as np

    data = np.asarray(values, dtype=float)
    if data.size == 0 or not bool(np.all(np.isfinite(data))):
        raise RuntimeError("trajectory metric is empty or non-finite")
    return {
        "mean": float(np.mean(data)),
        "median": float(np.median(data)),
        "p95": float(np.percentile(data, 95)),
        "stddev": float(np.std(data)),
        "min": float(np.min(data)),
        "max": float(np.max(data)),
        "final": float(data[-1]),
    }


def _graph_automorphisms(
    node_labels: list[tuple[str, str]],
    bonds: list[tuple[int, int]],
    *,
    max_permutations: int = 4096,
) -> list[tuple[int, ...]]:
    """Return exact labeled-graph automorphisms for symmetry-only atom swaps."""
    import networkx as nx

    graph = nx.Graph()
    for index, label in enumerate(node_labels):
        graph.add_node(index, label=tuple(label))
    graph.add_edges_from((int(left), int(right)) for left, right in bonds)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=lambda left, right: left.get("label") == right.get("label"),
    )
    permutations: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for mapping in matcher.isomorphisms_iter():
        permutation = tuple(int(mapping[index]) for index in range(len(node_labels)))
        if permutation in seen:
            continue
        seen.add(permutation)
        permutations.append(permutation)
        if len(permutations) > int(max_permutations):
            raise RuntimeError(
                "ligand symmetry exceeds the audited automorphism limit "
                f"({max_permutations})"
            )
    if not permutations:
        raise RuntimeError("could not construct a ligand symmetry mapping")
    return permutations


def _atom_element(atom) -> str:
    """Read an element from MDAnalysis, with a conservative Amber-name fallback."""
    try:
        element = str(atom.element or "").strip().upper()
    except Exception:
        element = ""
    if element:
        return element
    try:
        name = "".join(character for character in str(atom.name) if character.isalpha()).upper()
    except Exception:
        name = ""
    for two_letter in ("CL", "BR", "SI", "NA", "MG", "AL", "CA", "FE", "ZN"):
        if name.startswith(two_letter):
            return two_letter
    if name:
        return name[0]
    raise RuntimeError("ligand atom element is unavailable for symmetry mapping")


def _ligand_symmetry_permutations(ligand) -> list[tuple[int, ...]]:
    """Build symmetry maps from element, Amber type, and bonded topology."""
    global_to_local = {
        int(global_index): local_index
        for local_index, global_index in enumerate(ligand.indices)
    }
    labels = []
    for atom in ligand.atoms:
        try:
            atom_type = str(atom.type or "").strip().lower()
        except Exception:
            atom_type = ""
        labels.append((_atom_element(atom), atom_type))
    bonds = []
    try:
        bond_indices = ligand.bonds.indices
    except Exception as exc:
        raise RuntimeError(
            "ligand bonds are unavailable for symmetry-aware RMSD"
        ) from exc
    for left, right in bond_indices:
        left_local = global_to_local.get(int(left))
        right_local = global_to_local.get(int(right))
        if left_local is not None and right_local is not None:
            bonds.append((left_local, right_local))
    if len(labels) > 1 and not bonds:
        raise RuntimeError("ligand bonded topology is empty for symmetry-aware RMSD")
    return _graph_automorphisms(labels, bonds)


def _minimum_symmetry_rmsd(current, reference, permutations) -> float:
    """Minimum fixed-frame RMSD over atom automorphisms, without ligand refit."""
    import numpy as np

    current_array = np.asarray(current, dtype=float)
    reference_array = np.asarray(reference, dtype=float)
    if current_array.shape != reference_array.shape or current_array.ndim != 2:
        raise RuntimeError("ligand RMSD coordinate arrays are inconsistent")
    best = math.inf
    for permutation in permutations:
        mapped_reference = reference_array[np.asarray(permutation, dtype=int)]
        delta = current_array - mapped_reference
        value = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
        if value < best:
            best = value
    if not math.isfinite(best):
        raise RuntimeError("symmetry-aware ligand RMSD is non-finite")
    return best


def _analyze(payload: dict) -> dict:
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis import align
    from MDAnalysis.lib.distances import distance_array
    from MDAnalysis.transformations.nojump import NoJump

    prmtop = Path(payload["prmtop"]).resolve()
    trajectory = Path(payload["trajectory_dcd"]).resolve()
    reference_coordinates = Path(payload["reference_pdb"]).resolve()
    output_dir = Path(payload["output_dir"]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in (prmtop, trajectory, reference_coordinates):
        if not path.is_file():
            raise FileNotFoundError(path)

    alignment_selection = str(payload.get("alignment_selection", "protein and backbone"))
    ligand_selection = str(payload.get("ligand_selection", "resname LIG and not name H*"))
    contact_cutoff = float(payload.get("contact_cutoff_angstrom", 4.5))
    stride = int(payload.get("analysis_stride", 1))

    universe = mda.Universe(str(prmtop), str(trajectory))
    reference = mda.Universe(str(prmtop), str(reference_coordinates))
    universe.trajectory.add_transformations(NoJump())
    align_atoms = universe.select_atoms(alignment_selection)
    reference_align_atoms = reference.select_atoms(alignment_selection)
    ligand = universe.select_atoms(ligand_selection)
    reference_ligand = reference.select_atoms(ligand_selection)
    protein_heavy = universe.select_atoms("protein and not name H*")
    reference_protein_heavy = reference.select_atoms("protein and not name H*")
    if align_atoms.n_atoms < 3 or align_atoms.n_atoms != reference_align_atoms.n_atoms:
        raise RuntimeError("protein alignment selection is empty or inconsistent")
    if ligand.n_atoms < 1 or ligand.n_atoms != reference_ligand.n_atoms:
        raise RuntimeError("ligand analysis selection is empty or inconsistent")
    if protein_heavy.n_atoms < 1:
        raise RuntimeError("protein heavy-atom selection is empty")
    if protein_heavy.n_atoms != reference_protein_heavy.n_atoms:
        raise RuntimeError("protein heavy-atom reference selection is inconsistent")

    align.AlignTraj(
        universe, reference, select=alignment_selection, in_memory=True
    ).run()
    universe.trajectory[0]
    reference_ligand_positions = reference_ligand.positions.copy()
    reference_backbone_positions = reference_align_atoms.positions.copy()
    reference_center = reference_ligand.center_of_mass().copy()
    initial_distances = distance_array(
        reference_ligand.positions,
        reference_protein_heavy.positions,
        box=None,
    )
    contact_pairs = np.argwhere(initial_distances <= contact_cutoff)
    has_reference_contacts = bool(contact_pairs.size)
    pocket_indices = np.flatnonzero(np.any(initial_distances <= 6.0, axis=0))
    pocket = protein_heavy[pocket_indices]
    reference_pocket = reference_protein_heavy[pocket_indices]
    reference_com_distance = (
        float(np.linalg.norm(reference_ligand.center_of_mass() - reference_pocket.center_of_mass()))
        if pocket.n_atoms else None
    )
    symmetry_permutations = _ligand_symmetry_permutations(ligand)

    rows: list[dict] = []
    for ts in universe.trajectory[::stride]:
        backbone_delta = align_atoms.positions - reference_backbone_positions
        ligand_rmsd = _minimum_symmetry_rmsd(
            ligand.positions,
            reference_ligand_positions,
            symmetry_permutations,
        )
        backbone_rmsd = float(np.sqrt(np.mean(np.sum(backbone_delta * backbone_delta, axis=1))))
        center_displacement = float(np.linalg.norm(ligand.center_of_mass() - reference_center))
        pocket_com_drift = (
            abs(float(np.linalg.norm(ligand.center_of_mass() - pocket.center_of_mass())) - reference_com_distance)
            if reference_com_distance is not None else None
        )
        contact_retention = None
        if has_reference_contacts:
            pair_distances = distance_array(
                ligand.positions, protein_heavy.positions, box=None
            )[contact_pairs[:, 0], contact_pairs[:, 1]]
            contact_retention = float(np.mean(pair_distances <= contact_cutoff))
        values = (ligand_rmsd, backbone_rmsd, center_displacement)
        if contact_retention is not None:
            values = (*values, contact_retention)
        if pocket_com_drift is not None:
            values = (*values, pocket_com_drift)
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("MDAnalysis produced a non-finite trajectory metric")
        rows.append({
            "frame": int(ts.frame),
            "time_ps": float(ts.time),
            "ligand_heavy_atom_rmsd_angstrom": ligand_rmsd,
            "protein_backbone_rmsd_angstrom": backbone_rmsd,
            "ligand_center_displacement_angstrom": center_displacement,
            "ligand_pocket_com_distance_drift_angstrom": pocket_com_drift,
            "initial_contact_retention_fraction": contact_retention,
        })
    if len(rows) < 2:
        raise RuntimeError(
            "trajectory stability analysis requires at least two analyzed frames"
        )

    frame_csv = output_dir / "trajectory_metrics.csv"
    with frame_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    times = [row["time_ps"] for row in rows]
    ligand_rmsd_values = [row["ligand_heavy_atom_rmsd_angstrom"] for row in rows]
    slope = 0.0
    if len(rows) >= 2 and max(times) > min(times):
        slope = float(np.polyfit(times, ligand_rmsd_values, 1)[0])
    report_path = output_dir / "trajectory_stability.json"
    contact_metric = (
        {
            **_statistics([
                row["initial_contact_retention_fraction"] for row in rows
            ]),
            "applicability": "applicable",
        }
        if has_reference_contacts else {
            "applicability": "not_evaluable",
            "reason": (
                "the equilibrated production-start reference has no "
                "protein-ligand heavy-atom contacts "
                "at the configured cutoff"
            ),
        }
    )
    report = {
        "success": True,
        "action": "analyze",
        "status": "analyzed",
        "result_semantics": (
            "Short-trajectory evidence; the shared sbdd numeric rules assess symmetry-aware "
            "ligand RMSD median, pocket COM-distance drift P95 and mean initial-contact retention. "
            "These short simulations do not establish binding affinity or equilibrium stability."
        ),
        "mdanalysis_version": str(mda.__version__),
        "selections": {
            "alignment": alignment_selection,
            "ligand": ligand_selection,
            "contact_cutoff_angstrom": contact_cutoff,
            "initial_contact_pair_count": int(len(contact_pairs)),
            "pocket_heavy_atom_indices": pocket.indices.tolist(),
            "pocket_definition": "fixed protein heavy atoms within 6 A of reference ligand",
            "com_drift_definition": "absolute change in mass-weighted ligand-pocket COM distance",
            "reference_com_distance_angstrom": reference_com_distance,
            "ligand_rmsd_definition": (
                "minimum fixed-protein-frame RMSD over labeled-topology atom "
                "automorphisms; no ligand-only rotational/translational refit"
            ),
            "ligand_symmetry_automorphism_count": len(symmetry_permutations),
            "ligand_symmetry_node_labels": "element_and_amber_atom_type",
        },
        "reference": {
            "kind": "equilibrated_production_start_coordinates",
            "file": str(reference_coordinates),
            "sha256": _sha256(reference_coordinates),
        },
        "frame_count": len(rows),
        "metrics": {
            "ligand_heavy_atom_rmsd_angstrom": _statistics(ligand_rmsd_values),
            "protein_backbone_rmsd_angstrom": _statistics([
                row["protein_backbone_rmsd_angstrom"] for row in rows
            ]),
            "ligand_center_displacement_angstrom": _statistics([
                row["ligand_center_displacement_angstrom"] for row in rows
            ]),
            "initial_contact_retention_fraction": contact_metric,
            "ligand_pocket_com_distance_drift_angstrom": (
                _statistics([row["ligand_pocket_com_distance_drift_angstrom"] for row in rows])
                if reference_com_distance is not None else {
                    "applicability": "not_applicable", "reason": "no reference pocket atoms within 6 A"}
            ),
            "ligand_rmsd_drift_angstrom_per_ps": slope,
        },
        "pbc_handling": (
            "trajectory DCD stores unwrapped coordinates; NoJump is applied "
            "before protein-backbone alignment"
        ),
        "input_sha256": {
            "prmtop": _sha256(prmtop),
            "trajectory_dcd": _sha256(trajectory),
            "reference_pdb": _sha256(reference_coordinates),
        },
        "files": {"frame_metrics_csv": str(frame_csv), "report": str(report_path)},
        "file_sha256": {"frame_metrics_csv": _sha256(frame_csv)},
    }
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    # The report cannot recursively contain its own digest. Return the digest
    # alongside the on-disk report so the parent manifest can bind it.
    report["file_sha256"]["report"] = _sha256(report_path)
    return report


def run(payload: dict) -> dict:
    action = str(payload.get("action") or "")
    if action == "simulate":
        return _simulate(payload)
    if action == "analyze":
        return _analyze(payload)
    raise ValueError("action must be simulate or analyze")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    payload = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
    result = run(payload)
    output = Path(args.output_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

"""Audited short-MD orchestration across Amber topology and OpenMM runtimes."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from tool_runtime import ToolRuntimeError, run_tool


SCHEMA_VERSION = "openmm-short-md/v1"
WORKER = Path(__file__).resolve().with_name("openmm_md_worker.py")


class OpenMMShortMDError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _execution_record(result: Mapping[str, Any]) -> dict:
    return {
        "command": list(result.get("command") or []),
        "environment": result.get("environment"),
        "returncode": result.get("returncode"),
        "elapsed_seconds": result.get("elapsed_seconds"),
        "stdout_tail": str(result.get("stdout") or "")[-3000:],
        "stderr_tail": str(result.get("stderr") or "")[-3000:],
    }


def _run_worker(
    runtime: Mapping[str, Any],
    payload: dict,
    request_path: Path,
    response_path: Path,
    *,
    timeout: int,
) -> tuple[dict, dict]:
    _atomic_json(request_path, payload)
    try:
        execution = run_tool(
            dict(runtime),
            [WORKER, "--input-json", request_path, "--output-json", response_path],
            command="python",
            cwd=request_path.parent,
            timeout=int(timeout),
        )
    except ToolRuntimeError as exc:
        raise OpenMMShortMDError(str(exc)) from exc
    if not response_path.is_file():
        raise OpenMMShortMDError("isolated worker wrote no response JSON")
    response = json.loads(response_path.read_text(encoding="utf-8"))
    if response.get("success") is not True:
        raise OpenMMShortMDError("isolated worker returned success=false")
    return response, _execution_record(execution)


def run_openmm_short_md(
    candidate_topologies_manifest: str | Path,
    output_dir: str | Path,
    openmm_runtime: Mapping[str, Any],
    *,
    mdanalysis_runtime: Mapping[str, Any] | None = None,
    run_trajectory_analysis: bool = False,
    top_k: int = 3,
    platform: str = "CUDA",
    precision: str = "mixed",
    device_index: int = 0,
    temperature_kelvin: float = 300.0,
    pressure_bar: float = 1.0,
    timestep_femtoseconds: float = 2.0,
    friction_per_picosecond: float = 1.0,
    minimization_max_iterations: int = 500,
    equilibration_picoseconds: float = 20.0,
    production_picoseconds: float = 100.0,
    report_interval_picoseconds: float = 2.0,
    seed: int = 42,
    timeout_per_candidate_seconds: int = 86400,
    fail_fast: bool = False,
    alignment_selection: str = "protein and backbone",
    ligand_selection: str = "resname LIG and not name H*",
    contact_cutoff_angstrom: float = 4.5,
    analysis_stride: int = 1,
    analysis_timeout_per_candidate_seconds: int = 3600,
) -> dict:
    topology_path = Path(candidate_topologies_manifest).resolve()
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "short_md_manifest.json"
    stability_path = destination / "trajectory_stability_report.json"
    if not topology_path.is_file():
        raise OpenMMShortMDError("candidate topology manifest is missing")
    if not WORKER.is_file():
        raise OpenMMShortMDError(f"OpenMM worker is missing: {WORKER}")
    topology_manifest = json.loads(topology_path.read_text(encoding="utf-8"))
    if topology_manifest.get("schema_version") != "candidate-topologies/v1":
        raise OpenMMShortMDError("unsupported candidate topology manifest schema")
    raw_candidates = topology_manifest.get("candidates") or []
    if not isinstance(raw_candidates, list) or any(
        not isinstance(entry, dict) for entry in raw_candidates
    ):
        raise OpenMMShortMDError(
            "candidate topology manifest candidates must be an object list"
        )
    all_pose_ids = [str(entry.get("pose_id") or "") for entry in raw_candidates]
    if any(not value for value in all_pose_ids):
        raise OpenMMShortMDError(
            "candidate topology manifest contains a missing pose_id"
        )
    if len(set(all_pose_ids)) != len(all_pose_ids):
        raise OpenMMShortMDError(
            "candidate topology manifest contains duplicate pose_id"
        )
    candidates = raw_candidates[:int(top_k)]
    pose_ids = all_pose_ids[:int(top_k)]

    def manifest_file(value: Any) -> Path:
        path = Path(str(value or ""))
        if not path.is_absolute():
            path = topology_path.parent / path
        return path.resolve()

    results: list[dict] = []
    abort_reason: str | None = None
    for index, candidate in enumerate(candidates):
        pose_id = pose_ids[index]
        candidate_dir = destination / f"candidate_{index:05d}"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        record = {
            "pose_id": pose_id,
            "chemical_id": candidate.get("chemical_id"),
            "candidate_rank": index + 1,
            "upstream_topology_status": candidate.get("status"),
            "simulation_status": "not_evaluable",
            "simulation_reason": None,
            "analysis_status": (
                "not_evaluable" if run_trajectory_analysis else "not_requested"
            ),
            "analysis_reason": None,
            "simulation": None,
            "analysis": None,
            "simulation_input_sha256": None,
        }
        results.append(record)
        if abort_reason:
            record["simulation_reason"] = (
                "not run after fail_fast simulation failure: " + abort_reason
            )
            if run_trajectory_analysis:
                record["analysis_reason"] = "simulation was not run"
            continue
        if candidate.get("status") != "passed":
            record["simulation_reason"] = (
                "candidate topology is not evaluable: "
                + str(candidate.get("reason") or candidate.get("status"))
            )
            if run_trajectory_analysis:
                record["analysis_reason"] = "simulation was not run"
            continue
        files = dict(candidate.get("output_files") or {})
        declared_hashes = dict(candidate.get("output_sha256") or {})
        prmtop = manifest_file(files.get("solvated_complex_prmtop"))
        inpcrd = manifest_file(files.get("solvated_complex_inpcrd"))
        if not prmtop.is_file() or not inpcrd.is_file():
            record["simulation_status"] = "failed"
            record["simulation_reason"] = (
                "passed candidate lacks explicit-solvent Amber topology files"
            )
            if fail_fast:
                abort_reason = record["simulation_reason"]
            continue
        actual_hashes = {
            "solvated_complex_prmtop": _sha256(prmtop),
            "solvated_complex_inpcrd": _sha256(inpcrd),
        }
        record["simulation_input_sha256"] = actual_hashes
        hash_errors = [
            name
            for name, actual in actual_hashes.items()
            if not declared_hashes.get(name)
            or str(declared_hashes.get(name)) != actual
        ]
        if hash_errors:
            record["simulation_status"] = "failed"
            record["simulation_reason"] = (
                "candidate topology content hash is missing or mismatched for: "
                + ", ".join(hash_errors)
            )
            if run_trajectory_analysis:
                record["analysis_reason"] = "simulation input integrity failed"
            if fail_fast:
                abort_reason = record["simulation_reason"]
            continue
        try:
            simulation_payload = {
                "action": "simulate",
                "prmtop": str(prmtop.resolve()),
                "inpcrd": str(inpcrd.resolve()),
                "output_dir": str(candidate_dir / "simulation"),
                "platform": platform,
                "precision": precision,
                "device_index": int(device_index),
                "temperature_kelvin": float(temperature_kelvin),
                "pressure_bar": float(pressure_bar),
                "timestep_femtoseconds": float(timestep_femtoseconds),
                "friction_per_picosecond": float(friction_per_picosecond),
                "minimization_max_iterations": int(minimization_max_iterations),
                "equilibration_picoseconds": float(equilibration_picoseconds),
                "production_picoseconds": float(production_picoseconds),
                "report_interval_picoseconds": float(report_interval_picoseconds),
                "seed": int(seed) + index,
            }
            simulation, execution = _run_worker(
                openmm_runtime,
                simulation_payload,
                candidate_dir / "simulate_request.json",
                candidate_dir / "simulate_response.json",
                timeout=int(timeout_per_candidate_seconds),
            )
            record["simulation"] = {**simulation, "execution": execution}
            record["simulation_status"] = "passed"
        except Exception as exc:
            record["simulation_status"] = "failed"
            record["simulation_reason"] = str(exc)
            if run_trajectory_analysis:
                record["analysis_reason"] = "simulation failed"
            if fail_fast:
                abort_reason = str(exc)
            continue

        if run_trajectory_analysis:
            try:
                trajectory = Path(simulation["files"]["trajectory_dcd"])
                reference_pdb = Path(
                    simulation["files"]["production_reference_pdb"]
                )
                analysis_payload = {
                    "action": "analyze",
                    "prmtop": str(prmtop.resolve()),
                    "trajectory_dcd": str(trajectory.resolve()),
                    "reference_pdb": str(reference_pdb.resolve()),
                    "output_dir": str(candidate_dir / "analysis"),
                    "alignment_selection": alignment_selection,
                    "ligand_selection": ligand_selection,
                    "contact_cutoff_angstrom": float(contact_cutoff_angstrom),
                    "analysis_stride": int(analysis_stride),
                }
                analysis, analysis_execution = _run_worker(
                    mdanalysis_runtime or openmm_runtime,
                    analysis_payload,
                    candidate_dir / "analyze_request.json",
                    candidate_dir / "analyze_response.json",
                    timeout=int(analysis_timeout_per_candidate_seconds),
                )
                record["analysis"] = {**analysis, "execution": analysis_execution}
                record["analysis_status"] = "passed"
            except Exception as exc:
                # Keep a successful simulation successful; trajectory analysis
                # is a separate tool/result with separate failure attribution.
                record["analysis_status"] = "failed"
                record["analysis_reason"] = str(exc)

    def aggregate(statuses: list[str]) -> str:
        if not statuses:
            return "not_evaluable"
        if "failed" in statuses:
            return "failed"
        if all(value == "passed" for value in statuses):
            return "passed"
        if "passed" in statuses:
            return "partial"
        return "not_evaluable"

    simulation_status = aggregate([
        entry["simulation_status"] for entry in results
    ])
    analysis_status = (
        aggregate([entry["analysis_status"] for entry in results])
        if run_trajectory_analysis else "disabled"
    )
    counts = {
        state: sum(entry["simulation_status"] == state for entry in results)
        for state in ("passed", "failed", "not_evaluable")
    }
    counts["requested"] = len(candidates)
    analysis_counts = {
        state: sum(entry["analysis_status"] == state for entry in results)
        for state in ("passed", "failed", "not_evaluable", "not_requested")
    }
    analysis_counts["requested"] = len(candidates) if run_trajectory_analysis else 0
    for entry in results:
        # Backward-compatible aliases explicitly refer to simulation only.
        entry["status"] = entry["simulation_status"]
        entry["reason"] = entry["simulation_reason"]
    stability_report = {
        "schema_version": "mdanalysis-stability/v1",
        "success": analysis_status != "failed",
        "status": analysis_status,
        "result_semantics": (
            "Short-trajectory evidence; the shared sbdd numeric policy assesses "
            "RMSD median, pocket COM drift P95 and mean contact retention. "
            "This is not evidence of equilibrium binding stability."
        ),
        "candidates": [
            {
                "pose_id": entry["pose_id"],
                "status": entry["analysis_status"],
                "reason": entry["analysis_reason"],
                "analysis": entry["analysis"],
            }
            for entry in results
        ],
        "report_file": str(stability_path),
    }
    stability_report_sha256 = None
    if run_trajectory_analysis:
        _atomic_json(stability_path, stability_report)
        stability_report_sha256 = _sha256(stability_path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "success": simulation_status != "failed",
        "status": simulation_status,
        "simulation_status": simulation_status,
        "analysis_status": analysis_status,
        "result_semantics": (
            "Explicit-solvent periodic short MD for post-QC poses. This is a "
            "stability diagnostic, not an equilibrium free-energy calculation."
        ),
        "protocol": {
            "platform": platform,
            "precision": precision,
            "temperature_kelvin": float(temperature_kelvin),
            "pressure_bar": float(pressure_bar),
            "timestep_femtoseconds": float(timestep_femtoseconds),
            "friction_per_picosecond": float(friction_per_picosecond),
            "minimization_max_iterations": int(minimization_max_iterations),
            "equilibration_picoseconds": float(equilibration_picoseconds),
            "production_picoseconds": float(production_picoseconds),
            "report_interval_picoseconds": float(report_interval_picoseconds),
            "seed": int(seed),
            "device_index": int(device_index),
            "top_k": int(top_k),
            "fail_fast": bool(fail_fast),
            "trajectory_analysis_enabled": bool(run_trajectory_analysis),
            "analysis_stride": int(analysis_stride),
            "alignment_selection": alignment_selection,
            "ligand_selection": ligand_selection,
            "contact_cutoff_angstrom": float(contact_cutoff_angstrom),
            "trajectory_reference": "equilibrated_production_start_coordinates",
            "trajectory_pbc_handling": (
                "DCD written with enforcePeriodicBox=false; MDAnalysis applies "
                "NoJump before protein-backbone alignment"
            ),
            "worker_sha256": _sha256(WORKER),
            "candidate_topologies_sha256": _sha256(topology_path),
            "candidate_topology_protocol": topology_manifest.get("protocol"),
        },
        "counts": counts,
        "analysis_counts": analysis_counts,
        "source_candidate_count": len(raw_candidates),
        "candidate_topologies_manifest": str(topology_path),
        "candidates": results,
        "manifest_file": str(manifest_path),
        "trajectory_stability_report": (
            str(stability_path) if run_trajectory_analysis else None
        ),
        "trajectory_stability_report_sha256": stability_report_sha256,
    }
    _atomic_json(manifest_path, manifest)
    return manifest


__all__ = ["OpenMMShortMDError", "SCHEMA_VERSION", "run_openmm_short_md"]

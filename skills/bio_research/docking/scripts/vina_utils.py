"""AutoDock Vina adapter for the canonical Step 3 docking contract.

The installed Vina Python binding and the validated Meeko/RDKit installation
live in separate Conda environments.  This adapter keeps orchestration in
Drugmate while running preparation/export in ``reduce`` and docking in ``vina``.
"""
from __future__ import annotations

import json
import math
import tempfile
import time
from pathlib import Path

from docking_common import reference_ligand_box
from gnina_utils import normalize_docking_box
from tool_runtime import ToolRuntimeError, run_tool


SCRIPT_DIR = Path(__file__).resolve().parent
VINA_WORKER = SCRIPT_DIR / "vina_worker.py"
VINA_IO_WORKER = SCRIPT_DIR / "vina_io_worker.py"

DEFAULT_RUNTIME = {
    "kind": "composite",
    "components": {
        "engine": {
            "kind": "command",
            "environment": "vina",
            "command": "python",
        },
        "io": {
            "kind": "command",
            "environment": "reduce",
            "command": "python",
        },
    },
}


def run_autodock_vina_docking(
    *,
    receptor_pdbqt_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    runtime: dict | None = None,
    autobox_add: float = 4.0,
    scoring_function: str = "vina",
    exhaustiveness: int = 16,
    num_modes: int = 5,
    energy_range: float = 3.0,
    min_rmsd: float = 1.0,
    max_evals: int = 0,
    cpu: int = 1,
    seed: int = 42,
    spacing: float = 0.375,
    force_even_voxels: bool = False,
    no_refine: bool = False,
    timeout_seconds: int = 86400,
) -> dict:
    """Prepare, dock, export, and validate a multi-molecule Vina run."""
    if scoring_function not in {"vina", "vinardo"}:
        return {
            "success": False,
            "error": "AutoDock Vina scoring_function must be vina or vinardo",
        }
    if int(cpu) < 1:
        return {"success": False, "error": "AutoDock Vina cpu must be >= 1"}
    if float(spacing) <= 0:
        return {"success": False, "error": "AutoDock Vina spacing must be > 0"}
    if float(energy_range) <= 0:
        return {
            "success": False,
            "error": "AutoDock Vina energy_range must be > 0",
        }

    try:
        resolved_box = (
            reference_ligand_box(reference_ligand_file, autobox_add)
            if reference_ligand_file is not None
            else normalize_docking_box(docking_box)
        )
        engine_runtime, io_runtime = _runtime_components(runtime)
    except (TypeError, ValueError, ToolRuntimeError) as exc:
        return {"success": False, "error": str(exc)}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "vina_native"
    native_root.mkdir(parents=True, exist_ok=True)
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    prepared_dir = native_output_dir / "prepared_ligands"
    docked_native_dir = native_output_dir / "docked_pdbqt"
    docked_native_dir.mkdir(parents=True, exist_ok=True)
    docked_sdf = output_dir / "docked_poses.sdf"
    docked_sdf.unlink(missing_ok=True)
    log_file = output_dir / "autodock_vina.log"
    deadline = time.monotonic() + int(timeout_seconds)
    stage_records: list[dict] = []

    def execute_stage(
        label: str,
        stage_runtime: dict,
        args: list[str],
    ) -> dict:
        remaining = int(math.ceil(deadline - time.monotonic()))
        if remaining < 1:
            raise TimeoutError(
                f"AutoDock Vina total timeout exhausted before {label}"
            )
        process = run_tool(
            stage_runtime,
            args,
            timeout=remaining,
            check=False,
        )
        stage_records.append({"stage": label, **process})
        return process

    try:
        prepare_result_file = native_output_dir / "prepare_result.json"
        prepare_process = execute_stage(
            "ligand_preparation",
            io_runtime,
            [
                str(VINA_IO_WORKER),
                "prepare",
                "--input-sdf",
                str(ligand_file),
                "--output-dir",
                str(prepared_dir),
                "--result",
                str(prepare_result_file),
            ],
        )
        prepare_result = _read_result(prepare_result_file, "ligand preparation")
        if not prepare_process.get("success") or not prepare_result.get("success"):
            return _failure(
                prepare_result.get(
                    "error",
                    "AutoDock Vina ligand preparation rejected one or more records",
                ),
                log_file,
                stage_records,
                search_space=resolved_box,
                diagnostics=prepare_result,
                native_output_dir=native_output_dir,
            )

        engine_manifest = {
            "receptor_pdbqt": str(receptor_pdbqt_file),
            "docking_box": resolved_box,
            "native_output_dir": str(docked_native_dir),
            "status_path": str(native_output_dir / "engine_status.json"),
            "expected_molecules": int(prepare_result["total_records"]),
            "ligands": prepare_result["ligands"],
            "params": {
                "scoring_function": scoring_function,
                "exhaustiveness": int(exhaustiveness),
                "num_modes": int(num_modes),
                "energy_range": float(energy_range),
                "min_rmsd": float(min_rmsd),
                "max_evals": int(max_evals),
                "cpu": int(cpu),
                "seed": int(seed),
                "spacing": float(spacing),
                "force_even_voxels": bool(force_even_voxels),
                "no_refine": bool(no_refine),
                "verbosity": 0,
            },
        }
        engine_manifest_file = native_output_dir / "engine_manifest.json"
        engine_manifest_file.write_text(
            json.dumps(engine_manifest, indent=2, ensure_ascii=False)
        )
        engine_result_file = native_output_dir / "engine_result.json"
        engine_process = execute_stage(
            "vina_docking",
            engine_runtime,
            [
                str(VINA_WORKER),
                "--manifest",
                str(engine_manifest_file),
                "--result",
                str(engine_result_file),
            ],
        )
        engine_result = _read_result(engine_result_file, "Vina docking")
        if not engine_process.get("success") or not engine_result.get("success"):
            return _failure(
                engine_result.get(
                    "error",
                    "AutoDock Vina failed for one or more ligand preparations",
                ),
                log_file,
                stage_records,
                search_space=resolved_box,
                diagnostics=engine_result,
                native_output_dir=native_output_dir,
            )

        export_result_file = native_output_dir / "export_result.json"
        staged_sdf = native_output_dir / "docked_poses.sdf"
        export_process = execute_stage(
            "pose_export",
            io_runtime,
            [
                str(VINA_IO_WORKER),
                "export",
                "--run-result",
                str(engine_result_file),
                "--output-sdf",
                str(staged_sdf),
                "--result",
                str(export_result_file),
            ],
        )
        export_result = _read_result(export_result_file, "pose export")
        if not export_process.get("success") or not export_result.get("success"):
            return _failure(
                export_result.get(
                    "error",
                    "AutoDock Vina pose export failed validation",
                ),
                log_file,
                stage_records,
                search_space=resolved_box,
                diagnostics=export_result,
                native_output_dir=native_output_dir,
            )
        staged_sdf.replace(docked_sdf)
    except (OSError, TimeoutError, ToolRuntimeError, ValueError) as exc:
        return _failure(
            f"could not execute AutoDock Vina adapter: {exc}",
            log_file,
            stage_records,
            search_space=resolved_box,
            native_output_dir=native_output_dir,
        )

    _write_log(log_file, stage_records)
    redocking_validation = {
        "performed": False,
        "passed": None,
        "applicable": False,
        "not_applicable": True,
        "reason": "disabled_by_config",
        "gate": "not_applicable",
    }
    return {
        "success": True,
        "engine": "autodock_vina",
        "docked_sdf": str(docked_sdf),
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "num_poses": len(export_result["scores"]),
        "num_input_molecules": int(prepare_result["total_records"]),
        "num_preparation_variants": int(prepare_result["prepared_variants"]),
        "scores": export_result["scores"],
        "score_metadata": score_metadata(scoring_function),
        "redocking_validation": redocking_validation,
        "quality_gate_passed": None,
        "search_space": {
            "mode": (
                "reference_autobox"
                if reference_ligand_file is not None
                else "explicit_box"
            ),
            **resolved_box,
        },
        "software": {
            "engine": "AutoDock Vina Python API",
            "engine_version": engine_result.get("engine_version", "unknown"),
            "engine_runtime_environment": engine_runtime.get(
                "environment", "vina"
            ),
            "io_runtime_environment": io_runtime.get("environment", "reduce"),
            "ligand_io": "Meeko",
        },
        "commands": [_command_record(record) for record in stage_records],
        "export_diagnostics": export_result["diagnostics"],
        # Carried through so the failure registry can name the molecules and
        # preparation variants Vina could not export; empty for a complete batch.
        "operational_errors": export_result.get("operational_errors", []),
    }


def score_metadata(scoring_function: str = "vina") -> dict:
    """Machine-readable native AutoDock Vina score semantics."""
    return {
        "vina_affinity": {
            "meaning": f"AutoDock Vina {scoring_function} total energy",
            "unit": "kcal/mol",
            "better": "lower",
        },
        "ligand_efficiency": {
            "meaning": "Derived as -vina_affinity divided by heavy atom count",
            "unit": "kcal/mol/heavy atom",
            "better": "higher",
        },
        "engine_pose_rank": {
            "meaning": "One-based rank in the native output for one preparation",
            "unit": "rank",
            "better": "context_only",
        },
        "score_gap_1_2": {
            "meaning": "Second-best minus best affinity within one preparation",
            "unit": "kcal/mol",
            "better": "context_only",
        },
        "vina_inter": {
            "meaning": "Vina intermolecular energy component",
            "unit": "kcal/mol",
            "better": "lower",
        },
        "vina_intra": {
            "meaning": "Vina intramolecular energy component",
            "unit": "kcal/mol",
            "better": "context_only",
        },
        "vina_torsions": {
            "meaning": "Vina torsional energy component",
            "unit": "kcal/mol",
            "better": "context_only",
        },
        "vina_intra_best_pose": {
            "meaning": "Vina intramolecular energy of the best pose",
            "unit": "kcal/mol",
            "better": "context_only",
        },
    }


def _runtime_components(runtime: dict | None) -> tuple[dict, dict]:
    resolved = dict(runtime or DEFAULT_RUNTIME)
    if resolved.get("kind") != "composite":
        raise ToolRuntimeError(
            "AutoDock Vina runtime must be a composite engine + io declaration"
        )
    components = resolved.get("components")
    if not isinstance(components, dict):
        raise ToolRuntimeError(
            "AutoDock Vina composite runtime is missing components"
        )
    engine = components.get("engine")
    io = components.get("io")
    if not isinstance(engine, dict) or not isinstance(io, dict):
        raise ToolRuntimeError(
            "AutoDock Vina runtime requires engine and io components"
        )
    return dict(engine), dict(io)


def _read_result(path: Path, label: str) -> dict:
    if not path.is_file():
        raise ValueError(f"{label} did not write its result manifest")
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} wrote invalid result JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} result must be a JSON object")
    return payload


def _failure(
    error: str,
    log_file: Path,
    stage_records: list[dict],
    *,
    search_space: dict,
    native_output_dir: Path,
    diagnostics: dict | None = None,
) -> dict:
    _write_log(log_file, stage_records)
    return {
        "success": False,
        "error": str(error),
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "search_space": search_space,
        "commands": [_command_record(record) for record in stage_records],
        "diagnostics": diagnostics or {},
    }


def _write_log(log_file: Path, stage_records: list[dict]) -> None:
    sections = []
    for record in stage_records:
        sections.extend([
            f"[{record.get('stage', 'stage')}]",
            "command: "
            + " ".join(str(part) for part in (record.get("command") or [])),
            f"returncode: {record.get('returncode')}",
            "[stdout]",
            str(record.get("stdout") or ""),
            "[stderr]",
            str(record.get("stderr") or ""),
        ])
    log_file.write_text("\n".join(sections), encoding="utf-8")


def _command_record(record: dict) -> dict:
    return {
        key: record.get(key)
        for key in (
            "stage",
            "command",
            "returncode",
            "elapsed_seconds",
            "environment",
        )
    }

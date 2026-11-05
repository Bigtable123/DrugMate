"""Uni-Dock adapter for the canonical Step 3 docking contract.

Raw multi-molecule SDF files must not be sent directly to ``unidock
--gpu_batch``: without Uni-Dock's ``fragInfo``/``torsionInfo``/``atomInfo``
properties the ligands are treated as rigid.  This adapter therefore invokes
the installed ``unidocktools unidock_pipeline`` wrapper, then validates and
merges its per-ligand SDF outputs.
"""
from __future__ import annotations

import math
import re
import tempfile
from pathlib import Path

from docking_common import compound_name, reference_ligand_box
from gnina_utils import normalize_docking_box
from tool_runtime import ToolRuntimeError, run_tool


DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "unidock",
    "command": "unidocktools",
}

def build_unidock_pipeline_args(
    *,
    receptor_file: Path,
    ligand_file: Path,
    docking_box: dict,
    native_output_dir: Path,
    batch_size: int,
    scoring_function: str,
    search_mode: str,
    exhaustiveness: int,
    max_step: int,
    num_modes: int,
    refine_step: int,
    energy_range: float,
    seed: int,
) -> list[str]:
    """Return argv after ``unidocktools`` for one prepared pipeline run."""
    box = normalize_docking_box(docking_box)
    center = box["center"]
    size = box["size"]
    args = [
        "unidock_pipeline",
        "-r", str(receptor_file),
        "-l", str(ligand_file),
        "-cx", str(center[0]),
        "-cy", str(center[1]),
        "-cz", str(center[2]),
        "-sx", str(size[0]),
        "-sy", str(size[1]),
        "-sz", str(size[2]),
        "-sd", str(native_output_dir),
        "-bs", str(batch_size),
        "-sf", str(scoring_function),
        "-nm", str(num_modes),
        "-rs", str(refine_step),
        "-er", str(energy_range),
        "-topn", str(num_modes),
        "--seed", str(seed),
    ]
    # Uni-Dock Tools treats search_mode as an alternative to the explicit
    # exhaustiveness/max_step pair.
    if search_mode:
        args.extend(["-sm", str(search_mode)])
    else:
        args.extend(["-ex", str(exhaustiveness), "-ms", str(max_step)])
    return args


def run_unidock_docking(
    *,
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    runtime: dict | None = None,
    autobox_add: float = 4.0,
    batch_size: int = 18000,
    scoring_function: str = "vina",
    search_mode: str = "",
    exhaustiveness: int = 16,
    max_step: int = 20,
    num_modes: int = 5,
    refine_step: int = 3,
    energy_range: float = 3.0,
    device: int = 0,
    seed: int = 42,
    timeout_seconds: int = 86400,
) -> dict:
    """Run Uni-Dock and return the same engine-neutral fields as GNINA."""
    if scoring_function not in {"vina", "vinardo"}:
        return {
            "success": False,
            "error": "Uni-Dock adapter currently supports scoring_function vina or vinardo",
        }
    if search_mode not in {"", "fast", "balance", "detail"}:
        return {
            "success": False,
            "error": "Uni-Dock search_mode must be empty, fast, balance, or detail",
        }

    try:
        resolved_box = (
            reference_ligand_box(reference_ligand_file, autobox_add)
            if reference_ligand_file is not None
            else normalize_docking_box(docking_box)
        )
    except (TypeError, ValueError) as exc:
        return {"success": False, "error": str(exc)}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "unidock_native"
    native_root.mkdir(parents=True, exist_ok=True)
    # Never merge stale per-ligand files left by a prior failed/reconfigured run.
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    log_file = output_dir / "unidock.log"
    docked_sdf = output_dir / "docked_poses.sdf"

    args = build_unidock_pipeline_args(
        receptor_file=Path(receptor_file),
        ligand_file=Path(ligand_file),
        docking_box=resolved_box,
        native_output_dir=native_output_dir,
        batch_size=int(batch_size),
        scoring_function=scoring_function,
        search_mode=search_mode,
        exhaustiveness=int(exhaustiveness),
        max_step=int(max_step),
        num_modes=int(num_modes),
        refine_step=int(refine_step),
        energy_range=float(energy_range),
        seed=int(seed),
    )
    resolved_runtime = dict(runtime or DEFAULT_RUNTIME)
    try:
        process = run_tool(
            resolved_runtime,
            args,
            timeout=int(timeout_seconds),
            check=False,
            env={"CUDA_VISIBLE_DEVICES": str(device)},
        )
    except (OSError, ToolRuntimeError) as exc:
        return {
            "success": False,
            "error": f"could not execute Uni-Dock runtime: {exc}",
            "search_space": resolved_box,
        }

    log_file.write_text(
        (process.get("stdout") or "")
        + ("\n[stderr]\n" if process.get("stderr") else "")
        + (process.get("stderr") or ""),
        encoding="utf-8",
    )
    command_record = {
        key: process.get(key)
        for key in ("command", "returncode", "elapsed_seconds", "environment")
    }
    if not process.get("success"):
        return {
            "success": False,
            "error": f"Uni-Dock pipeline exited with code {process.get('returncode')}",
            "log_file": str(log_file),
            "command": command_record,
            "search_space": resolved_box,
        }

    try:
        input_compound_ids = _input_compound_ids(Path(ligand_file))
        expected_molecules = len(input_compound_ids)
        merge = _merge_unidock_outputs(
            native_output_dir=native_output_dir,
            output_sdf=docked_sdf,
            expected_molecules=expected_molecules,
            input_compound_ids=input_compound_ids,
        )
    except Exception as exc:
        return {
            "success": False,
            "error": f"could not validate Uni-Dock outputs: {exc}",
            "log_file": str(log_file),
            "command": command_record,
            "search_space": resolved_box,
        }
    if not merge["success"]:
        return {
            **merge,
            "log_file": str(log_file),
            "command": command_record,
            "search_space": resolved_box,
            "native_output_dir": str(native_output_dir),
        }

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
        "engine": "unidock",
        "docked_sdf": str(docked_sdf),
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "num_poses": len(merge["scores"]),
        "num_input_molecules": expected_molecules,
        "scores": merge["scores"],
        "score_metadata": score_metadata(scoring_function),
        "redocking_validation": redocking_validation,
        "quality_gate_passed": None,
        "search_space": {
            "mode": (
                "reference_autobox"
                if reference_ligand_file is not None else "explicit_box"
            ),
            **resolved_box,
        },
        "software": {
            "engine": "Uni-Dock",
            "engine_version": _unidock_version(resolved_runtime),
            "runtime_environment": resolved_runtime.get("environment", "unidock"),
            "frontend": "unidocktools unidock_pipeline",
        },
        "command": command_record,
        "merge_diagnostics": merge["diagnostics"],
        # Carried through so the failure registry can name the ligands Uni-Dock
        # could not score; empty for a complete batch.
        "operational_errors": merge.get("operational_errors", []),
    }


def _unidock_version(runtime: dict) -> str:
    """Query the actual Uni-Dock engine in the configured runtime."""
    try:
        result = run_tool(
            runtime,
            ["--version"],
            command="unidock",
            timeout=30,
            check=False,
        )
    except (OSError, ToolRuntimeError):
        return "unknown"
    output = "\n".join((result.get("stdout") or "", result.get("stderr") or ""))
    match = re.search(r"(?<!\w)v?(\d+(?:\.\d+){1,3}(?:[-+][\w.-]+)?)", output)
    return match.group(1) if match else "unknown"


def score_metadata(scoring_function: str = "vina") -> dict:
    """Machine-readable native Uni-Dock score semantics."""
    return {
        "docking_score": {
            "meaning": f"Uni-Dock {scoring_function} docking energy",
            "unit": "kcal/mol",
            "better": "lower",
        },
        "engine_pose_rank": {
            "meaning": "One-based rank in the native output for one ligand",
            "unit": "rank",
            "better": "context_only",
        },
    }


def _valid_input_molecule_count(ligand_file: Path) -> int:
    return len(_input_compound_ids(ligand_file))


def _input_compound_ids(ligand_file: Path) -> list[str]:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(
        str(ligand_file), removeHs=False, strictParsing=False,
    )
    identifiers = [
        compound_name(molecule, index)
        for index, molecule in enumerate(supplier)
        if molecule is not None
    ]
    if not identifiers:
        raise ValueError("input ligand SDF contains no parseable molecules")
    return identifiers


def _merge_unidock_outputs(
    *,
    native_output_dir: Path,
    output_sdf: Path,
    expected_molecules: int,
    input_compound_ids: list[str] | None = None,
) -> dict:
    """Validate per-ligand files and stream them into one canonical SDF."""
    from rdkit import Chem

    result_files = sorted(
        native_output_dir.glob("*.sdf"),
        key=lambda path: _natural_key(path.name),
    )
    partial_sdf = output_sdf.with_suffix(".partial.sdf")
    writer = Chem.SDWriter(str(partial_sdf))
    scores: list[dict] = []
    successful_files = 0
    successful_source_indices: set[int] = set()
    invalid_pose_records: list[dict] = []
    source_indices_by_id: dict[str, list[int]] = {}
    for index, identifier in enumerate(input_compound_ids or []):
        source_indices_by_id.setdefault(str(identifier), []).append(index)
    try:
        for fallback_source_index, result_file in enumerate(result_files):
            file_pose_count = 0
            supplier = Chem.SDMolSupplier(
                str(result_file), removeHs=False, strictParsing=False,
            )
            for engine_pose_index, molecule in enumerate(supplier):
                if molecule is None or not molecule.HasProp("docking_score"):
                    invalid_pose_records.append({
                        "file": str(result_file),
                        "engine_pose_index": engine_pose_index,
                        "reason": "missing molecule or docking_score",
                    })
                    continue
                try:
                    score = float(molecule.GetProp("docking_score"))
                except (TypeError, ValueError):
                    score = float("nan")
                if not math.isfinite(score):
                    invalid_pose_records.append({
                        "file": str(result_file),
                        "engine_pose_index": engine_pose_index,
                        "reason": "non-finite docking_score",
                    })
                    continue

                pose_index = len(scores)
                compound_id = (
                    molecule.GetProp("_Name").strip()
                    if molecule.HasProp("_Name") else result_file.stem
                ) or result_file.stem
                matches = source_indices_by_id.get(compound_id, [])
                source_molecule_index = (
                    matches[0] if len(matches) == 1 else fallback_source_index
                )
                molecule.SetProp("docking_engine", "unidock")
                molecule.SetProp("score_direction", "lower")
                molecule.SetIntProp("pose_index", pose_index)
                molecule.SetIntProp(
                    "source_molecule_index", source_molecule_index,
                )
                molecule.SetIntProp("engine_pose_index", engine_pose_index)
                molecule.SetIntProp("engine_pose_rank", engine_pose_index + 1)
                writer.write(molecule)
                scores.append({
                    "pose_index": pose_index,
                    "source_molecule_index": source_molecule_index,
                    "engine_pose_index": engine_pose_index,
                    "engine_pose_rank": engine_pose_index + 1,
                    "compound_id": compound_id,
                    "mol_name": compound_id,
                    "engine": "unidock",
                    "docking_score": score,
                })
                file_pose_count += 1
            if file_pose_count:
                successful_files += 1
                successful_source_indices.add(source_molecule_index)
    finally:
        writer.close()

    diagnostics = {
        "expected_molecules": expected_molecules,
        "result_sdf_files": len(result_files),
        "successful_molecule_files": successful_files,
        "successful_molecules": len(successful_source_indices),
        "invalid_pose_records": invalid_pose_records,
    }
    unscored_source_indices = sorted(
        set(range(expected_molecules)) - successful_source_indices
    )
    if unscored_source_indices:
        # A ligand Uni-Dock could not score is that ligand's missing value, not
        # the batch's failure: the executable policy scopes a single tool fault to the
        # single candidate it happened to, and the executable policy removes a missing value
        # from numerator and denominator rather than failing the run.  Returning
        # here threw away every pose the engine *did* produce, so one awkward
        # input ended docking for the whole cohort.  The shortfall stays visible
        # in ``diagnostics`` and is declared in ``operational_errors`` for the
        # failure registry (reflective-coordinator/SKILL.md:260-266); the
        # molecules that never scored simply carry no pose, which coverage
        # already reads as missing evidence.
        diagnostics["unscored_molecules"] = len(unscored_source_indices)
        diagnostics["unscored_source_molecule_indices"] = (
            unscored_source_indices
        )
        diagnostics["complete"] = False
    else:
        diagnostics["unscored_molecules"] = 0
        diagnostics["unscored_source_molecule_indices"] = []
        diagnostics["complete"] = True
    if not scores:
        return {
            "success": False,
            "error": "Uni-Dock produced no finite scored poses",
            "scores": [],
            "partial_docked_sdf": str(partial_sdf),
            "diagnostics": diagnostics,
        }

    partial_sdf.replace(output_sdf)
    merged = {
        "success": True,
        "scores": scores,
        "diagnostics": diagnostics,
    }
    if not diagnostics["complete"]:
        merged["operational_errors"] = [{
            "tool": "unidock",
            "reason": (
                "Uni-Dock output is incomplete: "
                f"{len(successful_source_indices)}/{expected_molecules} input molecules "
                "have at least one finite scored pose"
            ),
        }]
    return merged


def _natural_key(value: str) -> list[object]:
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    ]

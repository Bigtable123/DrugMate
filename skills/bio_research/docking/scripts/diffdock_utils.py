"""DiffDock adapter for the canonical Step 4 docking contract.

DiffDock is a blind generative pose model.  Its confidence output is not a
binding affinity and must never be combined numerically with Vina-family or
DOCK 6 energies.  This adapter isolates every run, splits a library into
one-ligand inputs, requires a complete ranked output for every input molecule,
and publishes canonical pose identities.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import tempfile
from pathlib import Path
from typing import Any

from deployment_config import deployment_path
from docking_common import (
    annotate_canonical_pose,
    compound_name,
    copy_identity_properties,
)
from tool_runtime import ToolRuntimeError, run_tool


DEFAULT_SOURCE_ROOT = deployment_path("${SBDD_DIFFDOCK_SOURCE_DIR}")
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "diffdock",
    "command": "python",
}
_CONFIDENCE_FILE = re.compile(
    r"^rank(?P<rank>\d+)_confidence"
    r"(?P<confidence>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\.sdf$"
)


def run_diffdock_docking(
    *,
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    runtime: dict | None = None,
    source_root: str | Path = DEFAULT_SOURCE_ROOT,
    model_dir: str | Path | None = None,
    confidence_model_dir: str | Path | None = None,
    checkpoint: str = "best_ema_inference_epoch_model.pt",
    confidence_checkpoint: str = "best_model_epoch75.pt",
    samples_per_complex: int = 5,
    batch_size: int = 10,
    inference_steps: int = 20,
    actual_steps: int | None = None,
    old_confidence_model: bool = True,
    device: int = 0,
    seed: int = 42,
    timeout_seconds: int = 86400,
) -> dict:
    """Run DiffDock and publish a complete canonical pose cohort."""
    source_root = Path(source_root).resolve()
    model_dir = Path(model_dir or source_root / "workdir/v1.1/score_model").resolve()
    confidence_model_dir = Path(
        confidence_model_dir
        or source_root / "workdir/v1.1/confidence_model"
    ).resolve()
    score_checkpoint = model_dir / checkpoint
    confidence_checkpoint_path = confidence_model_dir / confidence_checkpoint
    required = {
        "DiffDock inference entrypoint": source_root / "inference.py",
        "DiffDock score-model parameters": model_dir / "model_parameters.yml",
        "DiffDock score checkpoint": score_checkpoint,
        "DiffDock confidence-model parameters": (
            confidence_model_dir / "model_parameters.yml"
        ),
        "DiffDock confidence checkpoint": confidence_checkpoint_path,
    }
    missing = [label for label, path in required.items() if not path.is_file()]
    if missing:
        return {
            "success": False,
            "error": "missing required DiffDock resource(s): " + ", ".join(missing),
        }
    try:
        samples_per_complex = int(samples_per_complex)
        batch_size = int(batch_size)
        inference_steps = int(inference_steps)
        timeout_seconds = int(timeout_seconds)
        seed = int(seed)
        device = int(device)
        if samples_per_complex < 1 or batch_size < 1 or inference_steps < 1:
            raise ValueError(
                "samples_per_complex, batch_size, and inference_steps must be >= 1"
            )
        if actual_steps is not None:
            actual_steps = int(actual_steps)
            if actual_steps < 1 or actual_steps > inference_steps:
                raise ValueError("actual_steps must be between 1 and inference_steps")
    except (TypeError, ValueError) as exc:
        return {"success": False, "error": f"invalid DiffDock parameters: {exc}"}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "diffdock_native"
    native_root.mkdir(parents=True, exist_ok=True)
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    input_dir = native_output_dir / "inputs"
    result_dir = native_output_dir / "results"
    input_dir.mkdir()
    result_dir.mkdir()
    log_file = output_dir / "diffdock.log"
    docked_sdf = output_dir / "docked_poses.sdf"
    docked_sdf.unlink(missing_ok=True)

    try:
        molecules, input_rows = _prepare_inputs(
            ligand_file=Path(ligand_file),
            receptor_file=Path(receptor_file),
            input_dir=input_dir,
        )
    except Exception as exc:
        return {
            "success": False,
            "error": f"could not prepare DiffDock inputs: {exc}",
            "native_output_dir": str(native_output_dir),
        }

    csv_path = native_output_dir / "protein_ligand.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "complex_name",
                "protein_path",
                "ligand_description",
                "protein_sequence",
            ),
        )
        writer.writeheader()
        writer.writerows(input_rows)

    inference_config = {
        "model_dir": str(model_dir),
        "ckpt": checkpoint,
        "confidence_model_dir": str(confidence_model_dir),
        "confidence_ckpt": confidence_checkpoint,
        "samples_per_complex": samples_per_complex,
        "batch_size": batch_size,
        "inference_steps": inference_steps,
        "actual_steps": actual_steps,
        "old_score_model": False,
        "old_confidence_model": bool(old_confidence_model),
        "no_final_step_noise": True,
    }
    config_path = native_output_dir / "inference_config.json"
    config_path.write_text(
        json.dumps(inference_config, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    worker = Path(__file__).with_name("diffdock_worker.py").resolve()
    args = [
        str(worker),
        "--source-root", str(source_root),
        "--seed", str(seed),
        "--",
        "--config", str(config_path),
        "--protein_ligand_csv", str(csv_path),
        "--out_dir", str(result_dir),
        "--loglevel", "INFO",
    ]
    resolved_runtime = dict(runtime or DEFAULT_RUNTIME)
    try:
        process = run_tool(
            resolved_runtime,
            args,
            cwd=source_root,
            timeout=timeout_seconds,
            check=False,
            env={"CUDA_VISIBLE_DEVICES": str(device)},
        )
    except (OSError, ToolRuntimeError) as exc:
        return {
            "success": False,
            "error": f"could not execute DiffDock runtime: {exc}",
            "native_output_dir": str(native_output_dir),
        }
    _write_log(log_file, process)
    command_record = _command_record(process)
    if not process.get("success"):
        return {
            "success": False,
            "error": f"DiffDock exited with code {process.get('returncode')}",
            "log_file": str(log_file),
            "native_output_dir": str(native_output_dir),
            "command": command_record,
        }

    try:
        merge = _merge_diffdock_outputs(
            native_output_dir=result_dir,
            input_molecules=molecules,
            expected_poses_per_molecule=samples_per_complex,
            output_sdf=docked_sdf,
        )
    except Exception as exc:
        return {
            "success": False,
            "error": f"could not validate DiffDock outputs: {exc}",
            "log_file": str(log_file),
            "native_output_dir": str(native_output_dir),
            "command": command_record,
        }

    protocol_assets = {
        "score_checkpoint": _file_identity(score_checkpoint),
        "score_model_parameters": _file_identity(model_dir / "model_parameters.yml"),
        "confidence_checkpoint": _file_identity(confidence_checkpoint_path),
        "confidence_model_parameters": _file_identity(
            confidence_model_dir / "model_parameters.yml"
        ),
    }
    return {
        "success": True,
        "engine": "diffdock",
        "docked_sdf": str(docked_sdf),
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "num_poses": len(merge["scores"]),
        "num_input_molecules": len(molecules),
        "scores": merge["scores"],
        "score_metadata": score_metadata(),
        "redocking_validation": {
            "performed": False,
            "passed": None,
            "applicable": False,
            "not_applicable": True,
            "reason": "DiffDock route has no registered native-redocking protocol",
            "gate": "not_applicable",
        },
        "quality_gate_passed": None,
        "search_space": {
            "mode": "whole_receptor_blind_generative",
            "pocket_conditioned": False,
            "declared_pocket_context": docking_box,
            "reference_ligand_supplied": reference_ligand_file is not None,
        },
        "software": {
            "engine": "DiffDock",
            "engine_version": _source_revision(source_root),
            "runtime_environment": resolved_runtime.get("environment", "diffdock"),
            "frontend": "inference.py",
        },
        "protocol_assets": protocol_assets,
        "command": command_record,
        "merge_diagnostics": merge["diagnostics"],
        # Carried through so the failure registry can name the inputs and pose
        # files DiffDock could not deliver; empty for a complete batch.
        "operational_errors": merge.get("operational_errors", []),
    }


def _prepare_inputs(
    *,
    ligand_file: Path,
    receptor_file: Path,
    input_dir: Path,
) -> tuple[list[Any], list[dict[str, str]]]:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(
        str(ligand_file), removeHs=False, strictParsing=False,
    )
    molecules: list[Any] = []
    rows: list[dict[str, str]] = []
    invalid: list[int] = []
    for index, molecule in enumerate(supplier):
        if molecule is None:
            invalid.append(index)
            continue
        token = f"drugmate_src_{index:06d}"
        ligand_path = input_dir / f"{token}.sdf"
        writer = Chem.SDWriter(str(ligand_path))
        writer.write(molecule)
        writer.close()
        if not ligand_path.is_file() or ligand_path.stat().st_size == 0:
            raise ValueError(f"could not serialize ligand record {index}")
        molecules.append(molecule)
        rows.append({
            "complex_name": token,
            "protein_path": str(receptor_file.resolve()),
            "ligand_description": str(ligand_path.resolve()),
            "protein_sequence": "",
        })
    if invalid:
        raise ValueError(f"input ligand SDF has unparsable record(s): {invalid}")
    if not molecules:
        raise ValueError("input ligand SDF contains no parseable molecules")
    return molecules, rows


def _merge_diffdock_outputs(
    *,
    native_output_dir: Path,
    input_molecules: list[Any],
    expected_poses_per_molecule: int,
    output_sdf: Path,
) -> dict:
    from rdkit import Chem

    partial_sdf = output_sdf.with_suffix(".partial.sdf")
    partial_sdf.unlink(missing_ok=True)
    writer = Chem.SDWriter(str(partial_sdf))
    scores: list[dict] = []
    diagnostics = {
        "molecules": [],
        "invalid_records": [],
        "incomplete_molecules": [],
    }
    # Every gap below is recorded and stepped over rather than raised.  A pose
    # file DiffDock wrote badly, or an input it produced fewer ranks for, belongs
    # to that pose and that molecule: the executable policy scopes a single tool fault to
    # the single candidate, and the executable policy removes a missing value from numerator
    # and denominator instead of discarding the row.  Raising here threw away
    # every pose the run did produce, so one malformed output file ended docking
    # for the whole cohort.  Only a batch with no usable pose at all is a failure.
    try:
        for source_index, source_molecule in enumerate(input_molecules):
            token = f"drugmate_src_{source_index:06d}"
            complex_dir = native_output_dir / token
            ranked: dict[int, tuple[float, Path]] = {}
            if complex_dir.is_dir():
                for path in sorted(complex_dir.glob("*.sdf")):
                    match = _CONFIDENCE_FILE.fullmatch(path.name)
                    if match is None:
                        continue
                    rank = int(match.group("rank"))
                    confidence = float(match.group("confidence"))
                    if not math.isfinite(confidence):
                        diagnostics["invalid_records"].append({
                            "source_molecule_index": source_index,
                            "file": str(path),
                            "reason": "non-finite confidence",
                        })
                        continue
                    if rank in ranked:
                        diagnostics["invalid_records"].append({
                            "source_molecule_index": source_index,
                            "file": str(path),
                            "reason": f"duplicate DiffDock rank {rank}",
                        })
                        continue
                    ranked[rank] = (confidence, path)
            expected_ranks = set(range(1, expected_poses_per_molecule + 1))
            if set(ranked) != expected_ranks:
                diagnostics["incomplete_molecules"].append({
                    "source_molecule_index": source_index,
                    "expected_ranks": sorted(expected_ranks),
                    "found_ranks": sorted(ranked),
                })
            if not ranked:
                continue

            compound_id = compound_name(source_molecule, source_index)
            molecule_diagnostic = {
                "source_molecule_index": source_index,
                "compound_id": compound_id,
                "ranks": sorted(ranked),
            }
            diagnostics["molecules"].append(molecule_diagnostic)
            for rank in sorted(ranked):
                confidence, path = ranked[rank]
                parsed = [
                    molecule
                    for molecule in Chem.SDMolSupplier(
                        str(path), removeHs=False, strictParsing=False,
                    )
                    if molecule is not None
                ]
                if len(parsed) != 1:
                    diagnostics["invalid_records"].append({
                        "source_molecule_index": source_index,
                        "file": str(path),
                        "reason": (
                            "expected exactly one parseable pose; got "
                            f"{len(parsed)}"
                        ),
                    })
                    continue
                molecule = parsed[0]
                if molecule.GetNumConformers() != 1:
                    diagnostics["invalid_records"].append({
                        "source_molecule_index": source_index,
                        "file": str(path),
                        "reason": "does not contain exactly one conformer",
                    })
                    continue
                copy_identity_properties(source_molecule, molecule)
                pose_index = len(scores)
                molecule.SetDoubleProp("diffdock_confidence", confidence)
                molecule.SetProp("score_direction", "higher")
                identity = annotate_canonical_pose(
                    molecule,
                    compound_id=compound_id,
                    pose_index=pose_index,
                    source_molecule_index=source_index,
                    engine="diffdock",
                    engine_pose_index=rank - 1,
                    engine_pose_rank=rank,
                )
                writer.write(molecule)
                scores.append({
                    "pose_index": pose_index,
                    "source_molecule_index": source_index,
                    "engine_pose_index": rank - 1,
                    "engine_pose_rank": rank,
                    "compound_id": compound_id,
                    "mol_name": compound_id,
                    "engine": "diffdock",
                    "diffdock_confidence": confidence,
                    "score_reliability": "model_confidence_only",
                    **identity,
                })
    finally:
        writer.close()
    if not scores:
        raise ValueError("DiffDock produced no finite confidence-scored poses")
    diagnostics["complete"] = not (
        diagnostics["invalid_records"] or diagnostics["incomplete_molecules"]
    )
    partial_sdf.replace(output_sdf)
    merged = {"success": True, "scores": scores, "diagnostics": diagnostics}
    if not diagnostics["complete"]:
        merged["operational_errors"] = [{
            "tool": "diffdock",
            "reason": (
                "DiffDock output is incomplete: "
                f"{len(diagnostics['incomplete_molecules'])} input molecule(s) "
                "did not yield the expected pose ranks and "
                f"{len(diagnostics['invalid_records'])} pose file(s) were "
                "unusable"
            ),
        }]
    return merged


def score_metadata() -> dict:
    return {
        "diffdock_confidence": {
            "meaning": "DiffDock confidence-model output; not binding affinity",
            "unit": "model score",
            "better": "higher",
            "cross_engine_comparable": False,
        },
        "engine_pose_rank": {
            "meaning": "One-based rank assigned by DiffDock confidence",
            "unit": "rank",
            "better": "context_only",
        },
    }


def _write_log(path: Path, process: dict) -> None:
    path.write_text(
        (process.get("stdout") or "")
        + ("\n[stderr]\n" if process.get("stderr") else "")
        + (process.get("stderr") or ""),
        encoding="utf-8",
    )


def _command_record(process: dict) -> dict:
    return {
        key: process.get(key)
        for key in ("command", "returncode", "elapsed_seconds", "environment")
    }


def _file_identity(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": path.stat().st_size,
    }


def _source_revision(source_root: Path) -> str:
    git_dir = source_root / ".git"
    head = git_dir / "HEAD"
    try:
        value = head.read_text(encoding="utf-8").strip()
        if value.startswith("ref: "):
            value = (git_dir / value[5:]).read_text(encoding="utf-8").strip()
        return value[:40] if value else "unknown"
    except OSError:
        return "unknown"

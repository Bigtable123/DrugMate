"""Pocket2Mol structure-based molecule generation via local GPU execution.

Generates drug-like molecules by running Pocket2Mol locally on a GPU machine,
guided by pocket trimming from the interaction blueprint, then filters
results against blueprint constraints.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from deployment_config import deployment_path
from residue_utils import normalize_residue, extract_resnum
from tool_runtime import ToolRuntimeError, run_tool
from molecule_failures import MoleculeFailureTracker, is_rdkit_parse_error

_DEFAULT_TIMEOUT = 86400
_DEFAULT_POCKET2MOL_DIR = str(deployment_path("${SBDD_POCKET2MOL_SOURCE_DIR}"))
_DEFAULT_CONDA_ENV = "Pocket2Mol"
_DEFAULT_CHECKPOINT = str(deployment_path("${SBDD_POCKET2MOL_SOURCE_DIR}/pretrained_Pocket2Mol.pt"))
_DEFAULT_CHECKPOINT_SHA256 = (
    "e9d8152dfdf034dabf48e566904e9d7c6cc1b96843bdc4d0026d891678913911"
)
_POCKET_SHELL_RADIUS = 6.0
# Reporting bounds only. Molecular elimination belongs to the completed
# property/pose category folds, including for values outside these bounds.
_DRUGLIKE_FILTER = {"MW_min": 50, "MW_max": 1050, "LogP_min": -3, "LogP_max": 10, "RotBonds_max": 15}


def check_pocket2mol_available(
    pocket2mol_dir: str = _DEFAULT_POCKET2MOL_DIR,
    conda_env: str = _DEFAULT_CONDA_ENV,
    checkpoint_path: str = _DEFAULT_CHECKPOINT,
    conda_executable: str | None = None,
    checkpoint_sha256: str | None = None,
) -> bool:
    """Check if Pocket2Mol assets, conda env, and CUDA are ready locally."""
    return check_pocket2mol_status(
        pocket2mol_dir=pocket2mol_dir,
        conda_env=conda_env,
        checkpoint_path=checkpoint_path,
        conda_executable=conda_executable,
        checkpoint_sha256=checkpoint_sha256,
    )["available"]


def check_pocket2mol_status(
    pocket2mol_dir: str = _DEFAULT_POCKET2MOL_DIR,
    conda_env: str = _DEFAULT_CONDA_ENV,
    checkpoint_path: str = _DEFAULT_CHECKPOINT,
    conda_executable: str | None = None,
    checkpoint_sha256: str | None = None,
) -> dict:
    """Return a diagnostic readiness report for local Pocket2Mol execution."""
    p2m_dir = Path(pocket2mol_dir).expanduser()
    if not (p2m_dir / "sample_for_pdb.py").exists():
        return {
            "available": False,
            "reason": f"Pocket2Mol entry script not found: {p2m_dir / 'sample_for_pdb.py'}",
        }
    checkpoint = _resolve_path(checkpoint_path, p2m_dir)
    if not checkpoint.exists():
        return {"available": False, "reason": f"Pocket2Mol checkpoint not found: {checkpoint}"}
    actual_sha256 = _sha256(checkpoint)
    expected_sha256 = _expected_checkpoint_sha256(
        checkpoint, checkpoint_sha256
    )
    if expected_sha256 and actual_sha256 != expected_sha256:
        return {
            "available": False,
            "reason": (
                f"Pocket2Mol checkpoint checksum mismatch for {checkpoint}: "
                f"expected {expected_sha256}, got {actual_sha256}"
            ),
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": actual_sha256,
            "expected_checkpoint_sha256": expected_sha256,
            "checkpoint_checksum_verified": False,
        }
    conda_cmd = None
    try:
        probe_args = ["-c", "import torch; print(torch.cuda.is_available())"]
        if conda_executable is None and str(conda_env) == _DEFAULT_CONDA_ENV:
            # Use the same registered placement as actual generation. A local
            # prefix override must also apply to this readiness check.
            runtime = {"kind": "command", "environment": "pocket2mol", "command": "python"}
            args = probe_args
        else:
            conda_cmd = _find_conda_executable(conda_env, conda_executable)
            runtime = {"kind": "command", "environment": "current", "command": conda_cmd}
            args = ["run", *_conda_env_args(conda_env), "python", *probe_args]
        result = run_tool(runtime, args, timeout=30, check=False)
    except Exception as exc:
        return {"available": False, "reason": f"Pocket2Mol CUDA probe failed to start: {exc}"}

    stdout = str(result.get("stdout") or "").strip()
    stderr = str(result.get("stderr") or "").strip()
    if result.get("returncode") == 124:
        return {
            "available": False,
            "reason": "Pocket2Mol CUDA probe timed out after 30s",
        }
    if result.get("returncode") != 0:
        return {
            "available": False,
            "reason": (
                "Pocket2Mol CUDA probe exited with code "
                f"{result.get('returncode')}: {stderr or stdout}"
            ),
        }
    if "True" not in stdout:
        return {
            "available": False,
            "reason": f"Pocket2Mol conda env loaded but torch.cuda.is_available() was not True: {stdout}",
        }
    return {
        "available": True,
        "reason": "Pocket2Mol assets, conda env, and CUDA are ready",
        "conda_executable": conda_cmd,
        "conda_env": conda_env,
        "runtime_environment": result.get("environment"),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": actual_sha256,
        "expected_checkpoint_sha256": expected_sha256,
        "checkpoint_checksum_verified": bool(expected_sha256),
        "cuda_probe_stdout": stdout,
    }


def prepare_pocket_pdb(
    protein_pdb: Path,
    blueprint_residues: list[str],
    chain: str,
    output_dir: Path,
    shell_radius: float = _POCKET_SHELL_RADIUS,
) -> Path:
    """Trim protein PDB to blueprint residues + surrounding shell.

    Keeps only:
    1. Residues explicitly listed in blueprint (Required + Desired)
    2. Residues within shell_radius of any blueprint residue atom

    This biases Pocket2Mol generation toward key interaction sites.
    """
    target_resnums = set()
    for res in blueprint_residues:
        num = extract_resnum(res)
        if num is not None:
            target_resnums.add(num)

    if not target_resnums:
        # No blueprint residues — use full protein as pocket
        return protein_pdb

    # Pass 1: collect coordinates of blueprint residue atoms
    blueprint_coords = []
    all_lines = []
    with open(protein_pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            all_lines.append(line)
            line_chain = line[21].strip()
            if chain and line_chain != chain:
                continue
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            if resseq in target_resnums:
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    blueprint_coords.append((x, y, z))
                except (ValueError, IndexError):
                    pass

    if not blueprint_coords:
        return protein_pdb

    # Pass 2: select residues within shell_radius of any blueprint atom
    shell_resnums = set(target_resnums)
    for line in all_lines:
        line_chain = line[21].strip()
        if chain and line_chain != chain:
            continue
        try:
            resseq = int(line[22:26].strip())
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except (ValueError, IndexError):
            continue
        if resseq in shell_resnums:
            continue
        for bx, by, bz in blueprint_coords:
            dist_sq = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
            if dist_sq <= shell_radius ** 2:
                shell_resnums.add(resseq)
                break

    # Write trimmed PDB
    output_dir.mkdir(parents=True, exist_ok=True)
    descriptor, trimmed_name = tempfile.mkstemp(
        prefix="pocket2mol_pocket-",
        suffix=".pdb",
        dir=str(output_dir),
    )
    os.close(descriptor)
    trimmed_path = Path(trimmed_name)
    with open(protein_pdb) as fin, open(trimmed_path, "w") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                line_chain = line[21].strip()
                if chain and line_chain != chain:
                    continue
                try:
                    resseq = int(line[22:26].strip())
                except ValueError:
                    continue
                if resseq in shell_resnums:
                    fout.write(line)
            elif line.startswith("END"):
                fout.write(line)
        fout.write("END\n")

    return trimmed_path


def run_pocket2mol(
    pocket_pdb: Path,
    center: list[float] | tuple[float, float, float] | str,
    num_samples: int = 100,
    pocket2mol_dir: str = _DEFAULT_POCKET2MOL_DIR,
    conda_env: str = _DEFAULT_CONDA_ENV,
    checkpoint_path: str = _DEFAULT_CHECKPOINT,
    gpu_device: int = 0,
    taskset_cpu: int | None = 0,
    timeout_seconds: int = _DEFAULT_TIMEOUT,
    output_dir: Path | None = None,
    conda_executable: str | None = None,
    seed: int = 42,
    checkpoint_sha256: str | None = None,
    return_report: bool = False,
    failure_tracker: MoleculeFailureTracker | None = None,
) -> list[dict] | tuple[list[dict], dict]:
    """Run Pocket2Mol locally on GPU machine.

    Returns molecule dicts containing SMILES plus the unmodified generated
    MolBlock so downstream assembly can preserve the original 3D pose.
    """
    work_root = (
        Path(output_dir)
        if output_dir is not None
        else Path(tempfile.mkdtemp(prefix="p2m_local_"))
    )
    work_root.mkdir(parents=True, exist_ok=True)
    # A new run-* directory is the invocation boundary.  Pocket2Mol may create
    # nested output directories, but parsing never escapes this fresh directory.
    run_dir = Path(tempfile.mkdtemp(prefix="run-", dir=str(work_root)))
    raw_output_dir = run_dir / "output"
    p2m_dir = Path(pocket2mol_dir).expanduser()
    entry_script = p2m_dir / "sample_for_pdb.py"
    checkpoint = _resolve_path(checkpoint_path, p2m_dir)
    expected_sha256 = _expected_checkpoint_sha256(
        checkpoint, checkpoint_sha256
    )
    provenance = {
        "tool": "pocket2mol",
        "work_dir": str(run_dir),
        "output_dir": str(raw_output_dir),
        "seed": int(seed),
        "num_samples_requested": int(num_samples),
        "checkpoint": str(checkpoint),
        "expected_checkpoint_sha256": expected_sha256,
        "gpu_device": int(gpu_device),
        "taskset_cpu": taskset_cpu,
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    try:
        if not entry_script.is_file():
            raise RuntimeError(f"Pocket2Mol entry script not found: {entry_script}")
        if not checkpoint.is_file():
            raise RuntimeError(f"Pocket2Mol checkpoint not found: {checkpoint}")
        actual_sha256 = _sha256(checkpoint)
        provenance["checkpoint_sha256"] = actual_sha256
        provenance["checkpoint_checksum_verified"] = False
        if expected_sha256 and actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Pocket2Mol checkpoint checksum mismatch for {checkpoint}: "
                f"expected {expected_sha256}, got {actual_sha256}"
            )
        provenance["checkpoint_checksum_verified"] = bool(expected_sha256)

        pocket_copy = run_dir / "pocket.pdb"
        shutil.copy2(pocket_pdb, pocket_copy)
        provenance["pocket_pdb"] = str(pocket_copy)
        provenance["pocket_pdb_sha256"] = _sha256(pocket_copy)

        config_path = run_dir / "sample_for_pdb.yml"
        _write_pocket2mol_config(
            template_path=p2m_dir / "configs" / "sample_for_pdb.yml",
            output_path=config_path,
            checkpoint_path=checkpoint,
            num_samples=num_samples,
            seed=seed,
        )
        provenance["config"] = str(config_path)
        provenance["config_sha256"] = _sha256(config_path)

        sample_cmd = [
            "python", str(entry_script),
            "--pdb_path", str(pocket_copy),
            "--center", _format_center(center),
            "--config", str(config_path),
            "--device", "cuda",
            "--outdir", str(raw_output_dir),
        ]
        taskset_prefix = _taskset_prefix(taskset_cpu)
        if taskset_prefix:
            sample_cmd = taskset_prefix + sample_cmd

        execution = _execute_pocket2mol(
            sample_cmd=sample_cmd,
            conda_env=conda_env,
            conda_executable=conda_executable,
            cwd=p2m_dir,
            timeout_seconds=int(timeout_seconds),
            gpu_device=int(gpu_device),
        )
        execution = _write_execution_provenance(
            run_dir, execution, provenance=provenance
        )
        if not execution.get("success"):
            raise RuntimeError(
                execution.get("error")
                or "Pocket2Mol execution returned a non-zero status"
            )
    except Exception as exc:
        if not (run_dir / "execution.json").is_file():
            _write_execution_provenance(
                run_dir,
                {
                    "success": False,
                    "returncode": None,
                    "stdout": "",
                    "stderr": str(exc),
                    "error": str(exc),
                },
                provenance=provenance,
            )
        raise RuntimeError(
            f"Pocket2Mol invocation failed in {run_dir}: {exc}"
        ) from exc

    try:
        compounds, parse_diagnostics = _parse_generated_sdfs(
            raw_output_dir, return_diagnostics=True, failure_tracker=failure_tracker,
        )
    except Exception as exc:
        execution = _finalize_execution_output_validation(
            execution,
            success=False,
            diagnostics={"parse_error": str(exc)},
            error=f"Pocket2Mol output parsing failed: {exc}",
        )
        raise RuntimeError(
            f"Pocket2Mol output parsing failed in {run_dir}: {exc}"
        ) from exc
    execution = _finalize_execution_output_validation(
        execution,
        success=bool(compounds),
        diagnostics=parse_diagnostics,
        error=(
            None
            if compounds
            else "Pocket2Mol produced no valid sanitized 3D molecules"
        ),
    )
    report = {
        "tool": "pocket2mol",
        "success": bool(compounds),
        "num_samples_requested": int(num_samples),
        "num_valid": len(compounds),
        "work_dir": str(run_dir),
        "output_dir": str(raw_output_dir),
        "seed": int(seed),
        "checkpoint": str(checkpoint),
        "asset_sha256": provenance.get("checkpoint_sha256"),
        "diagnostics": parse_diagnostics,
        "execution": execution,
    }
    if return_report:
        return compounds, report
    return compounds


def filter_by_blueprint(
    generated_mols: list[dict],
    blueprint_file: str | None,
    pocket_center: list[float] | None,
    return_diagnostics: bool = False,
    failure_tracker: MoleculeFailureTracker | None = None,
) -> list[dict] | tuple[list[dict], dict]:
    """Annotate raw candidates; only unreadable structures are removed.

    Partial property or 2D pharmacophore checks cannot decide a category gate.
    Keep their marks for audit and let property-filter / structure-evaluator
    decide after collecting the complete planned category evidence.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from metric_scoring import PROPERTY_FILTER_METRIC_RULES, assess_record

    diagnostics = {
        "raw_generated": len(generated_mols), "parsed": 0,
        "druglike_pass": 0, "blueprint_pass": 0, "kept": 0,
        "blueprint_filter": "not_provided",
        "rejection_reasons": {"parse_failed": 0},
        "thresholds": dict(_DRUGLIKE_FILTER),
        "decision_role": "marks_only_until_category_fold",
    }
    blueprint = _load_blueprint_constraints(blueprint_file) if blueprint_file else None
    if blueprint_file:
        diagnostics["blueprint_filter"] = (
            "annotation_only" if blueprint else "not_loaded_or_empty"
        )
    metric_ids = {"rdkit.exact_mol_wt", "rdkit.clogp", "rdkit.rotatable_bonds"}
    rules = tuple(rule for rule in PROPERTY_FILTER_METRIC_RULES
                  if rule["metric_id"] in metric_ids)
    kept = []
    tracker = failure_tracker or MoleculeFailureTracker()
    for source in generated_mols:
        try:
            mol = Chem.MolFromSmiles(source.get("smiles", ""))
            if mol is None or mol.GetNumAtoms() == 0:
                diagnostics["rejection_reasons"]["parse_failed"] += 1
                tracker.record(source, "pocket2mol.annotations", "molecule_parse_failed",
                               "RDKit could not parse the generated SMILES")
                continue
            item = dict(source)
            props = {
                "MW": Descriptors.ExactMolWt(mol),
                "LogP": Descriptors.MolLogP(mol),
                "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            }
            item.update(props)
            item["HBD"] = rdMolDescriptors.CalcNumHBD(mol)
            item["HBA"] = rdMolDescriptors.CalcNumHBA(mol)
            marks = assess_record({"properties": props}, rules)
            item["generation_metric_labels"] = {
                mark["metric_id"]: mark["label"] for mark in marks
            }
            diagnostics["parsed"] += 1
            diagnostics["druglike_pass"] += int(
                all(mark["label"] != "fail" for mark in marks)
            )
            if blueprint:
                matched = _matches_blueprint_pharmacophore(item, blueprint)
                item["blueprint_filter"] = "passed" if matched else "flag"
                diagnostics["blueprint_pass"] += int(matched)
            kept.append(item)
        except Exception as exc:
            tracker.record(source, "pocket2mol.annotations", "generated_annotation_failed", str(exc),
                           exception=None if is_rdkit_parse_error(exc) else exc)
            continue

    diagnostics["kept"] = len(kept)
    return _filter_result(kept, diagnostics, return_diagnostics)


# ===== Internal helpers =====


def _filter_result(
    kept: list[dict],
    diagnostics: dict,
    return_diagnostics: bool,
) -> list[dict] | tuple[list[dict], dict]:
    if return_diagnostics:
        return kept, diagnostics
    return kept


def _resolve_path(path_value: str | Path, base_dir: Path) -> Path:
    path = Path(path_value).expanduser()
    if path.is_absolute():
        return path
    return base_dir / path


def _expected_checkpoint_sha256(
    checkpoint: Path,
    configured: str | None,
) -> str | None:
    if configured:
        return str(configured).strip().lower()
    try:
        if checkpoint.resolve() == Path(_DEFAULT_CHECKPOINT).resolve():
            return _DEFAULT_CHECKPOINT_SHA256
    except (OSError, RuntimeError):
        pass
    return None


def _find_conda_executable(conda_env: str, configured: str | None = None) -> str:
    """Find conda even when the interactive shell has not initialized it."""
    if configured:
        return str(Path(configured).expanduser())

    if os.environ.get("CONDA_EXE"):
        return str(Path(os.environ["CONDA_EXE"]).expanduser())

    conda = shutil.which("conda")
    if conda:
        return conda

    env_file = Path.home() / ".conda" / "environments.txt"
    try:
        for line in env_file.read_text().splitlines():
            env_path = Path(line.strip()).expanduser()
            if not env_path:
                continue
            if env_path.name == conda_env:
                try:
                    conda_candidate = env_path.parents[1] / "bin" / "conda"
                    return str(conda_candidate)
                except IndexError:
                    continue
    except Exception:
        pass

    return "conda"


def _conda_env_args(conda_env: str) -> list[str]:
    env_value = str(conda_env).strip()
    if not env_value:
        raise ValueError("Pocket2Mol conda environment name/path is empty")

    env_path = Path(env_value).expanduser()
    is_explicit_path = (
        env_path.is_absolute()
        or env_value.startswith("~")
        or os.sep in env_value
        or (os.altsep is not None and os.altsep in env_value)
    )
    if is_explicit_path:
        return ["-p", str(env_path)]
    return ["-n", env_value]


def _format_center(center: list[float] | tuple[float, float, float] | str) -> str:
    values = _coerce_center(center)
    center_str = f"{values[0]},{values[1]},{values[2]}"
    if center_str.startswith("-"):
        return f" {center_str}"
    return center_str


def _coerce_center(center: list[float] | tuple[float, float, float] | str) -> list[float]:
    if isinstance(center, str):
        parts = center.replace(";", ",").split(",")
        if len(parts) != 3:
            raise ValueError("Pocket2Mol center string must contain three comma-separated values")
        raw_values = parts
    elif isinstance(center, (list, tuple)) and len(center) == 3:
        raw_values = center
    else:
        raise ValueError("Pocket2Mol center must be a 3-element list/tuple or comma-separated string")

    try:
        values = [float(value) for value in raw_values]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Pocket2Mol center contains non-numeric values: {center}") from exc
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"Pocket2Mol center contains non-finite values: {center}")
    return values


def _taskset_prefix(taskset_cpu: int | None) -> list[str]:
    if taskset_cpu is None:
        return []
    try:
        cpu = int(taskset_cpu)
    except (TypeError, ValueError):
        return []
    if cpu < 0:
        return []
    taskset = shutil.which("taskset")
    if not taskset:
        return []
    return [taskset, "-c", str(cpu)]


def _write_pocket2mol_config(
    template_path: Path,
    output_path: Path,
    checkpoint_path: Path,
    num_samples: int,
    seed: int,
) -> None:
    """Write a per-run config with explicit checkpoint, sample count, and seed."""
    if not template_path.exists():
        raise RuntimeError(f"Pocket2Mol config template not found: {template_path}")

    checkpoint_value = str(checkpoint_path).replace('"', '\\"')
    lines = template_path.read_text().splitlines()
    wrote_checkpoint = False
    wrote_num_samples = False
    wrote_seed = False
    updated_lines = []

    for line in lines:
        stripped = line.strip()
        indent = line[: len(line) - len(line.lstrip())]
        if stripped.startswith("checkpoint:"):
            updated_lines.append(f'{indent}checkpoint: "{checkpoint_value}"')
            wrote_checkpoint = True
        elif stripped.startswith("num_samples:"):
            updated_lines.append(f"{indent}num_samples: {int(num_samples)}")
            wrote_num_samples = True
        elif stripped.startswith("seed:"):
            updated_lines.append(f"{indent}seed: {int(seed)}")
            wrote_seed = True
        else:
            updated_lines.append(line)

    if not wrote_checkpoint:
        _insert_yaml_field(
            updated_lines,
            section="model",
            field_line=f'  checkpoint: "{checkpoint_value}"',
        )
    missing_sample_fields = []
    if not wrote_seed:
        missing_sample_fields.append(f"  seed: {int(seed)}")
    if not wrote_num_samples:
        missing_sample_fields.append(f"  num_samples: {int(num_samples)}")
    for field_line in reversed(missing_sample_fields):
        _insert_yaml_field(
            updated_lines,
            section="sample",
            field_line=field_line,
        )

    output_path.write_text("\n".join(updated_lines) + "\n")


def _insert_yaml_field(
    lines: list[str],
    *,
    section: str,
    field_line: str,
) -> None:
    """Insert a field into a top-level mapping without duplicating the section."""
    for index, line in enumerate(lines):
        if line == f"{section}:":
            lines.insert(index + 1, field_line)
            return
    if lines and lines[-1]:
        lines.append("")
    lines.extend([f"{section}:", field_line])


def _parse_generated_sdfs(
    raw_dir: Path,
    *,
    return_diagnostics: bool = False,
    failure_tracker: MoleculeFailureTracker | None = None,
) -> list[dict] | tuple[list[dict], dict]:
    """Parse native 3D outputs and retain record-level rejection provenance."""
    from rdkit import Chem

    tracker = failure_tracker or MoleculeFailureTracker()
    start = len(tracker.failures)
    results = []
    diagnostics = dict.fromkeys(("sdf_files", "records", "invalid_or_unsanitized",
                                "missing_conformer", "not_marked_3d", "nonfinite_coordinates",
                                "multicomponent", "valid"), 0)
    sdf_files = sorted(Path(raw_dir).rglob("*.sdf"))
    diagnostics["sdf_files"] = len(sdf_files)
    for sdf_file in sdf_files:
        supplier = Chem.SDMolSupplier(str(sdf_file), sanitize=True, removeHs=False, strictParsing=True)
        for mol_index, mol in enumerate(supplier):
            record = {"source": "pocket2mol", "candidate_role": "generated",
                      "source_file": str(sdf_file.resolve()), "source_mol_index": mol_index}
            diagnostics["records"] += 1
            reason = None
            try:
                if mol is None:
                    reason = "invalid_or_unsanitized"
                else:
                    Chem.SanitizeMol(mol)
                    if mol.GetNumConformers() < 1:
                        reason = "missing_conformer"
                    elif not mol.GetConformer().Is3D():
                        reason = "not_marked_3d"
                    elif not all(math.isfinite(float(v)) for row in mol.GetConformer().GetPositions() for v in row):
                        reason = "nonfinite_coordinates"
                    elif len(Chem.GetMolFrags(mol, sanitizeFrags=False)) != 1:
                        reason = "multicomponent"
                    else:
                        record["smiles"] = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)
                        record["raw_pose_molblock"] = Chem.MolToMolBlock(mol)
                        record["coordinate_origin"] = "pocket2mol_raw"
                        if not record["smiles"]:
                            reason = "invalid_or_unsanitized"
            except Exception as exc:
                diagnostics["invalid_or_unsanitized"] += 1
                tracker.record(record, "pocket2mol.parse", "generated_parse_failed", str(exc),
                               exception=None if is_rdkit_parse_error(exc) else exc)
                continue
            if reason:
                diagnostics[reason] += 1
                tracker.record(record, "pocket2mol.parse", "generated_" + reason, reason)
                continue
            results.append(record)
    diagnostics["valid"] = len(results)
    if failure_tracker is not None:
        diagnostics["molecule_failures"] = tracker.failures[start:]
    if return_diagnostics:
        return results, diagnostics
    return results


def _execute_pocket2mol(
    *,
    sample_cmd: list[str],
    conda_env: str,
    conda_executable: str | None,
    cwd: Path,
    timeout_seconds: int,
    gpu_device: int,
) -> dict:
    """Execute through the shared runtime when using the registered environment."""
    execution_env = {"CUDA_VISIBLE_DEVICES": str(gpu_device)}
    if conda_executable is None and str(conda_env) == _DEFAULT_CONDA_ENV:
        runtime = {
            "kind": "command",
            "environment": "pocket2mol",
            "command": sample_cmd[0],
        }
        try:
            return run_tool(
                runtime,
                sample_cmd[1:],
                cwd=cwd,
                timeout=timeout_seconds,
                check=True,
                env=execution_env,
            )
        except ToolRuntimeError as exc:
            if exc.result:
                return {**exc.result, "error": str(exc)}
            raise

    conda_cmd = _find_conda_executable(conda_env, conda_executable)
    command = [
        conda_cmd,
        "run",
        *_conda_env_args(conda_env),
        *sample_cmd,
    ]
    try:
        record = run_tool(
            {
                "kind": "command",
                "environment": "current",
                "command": command[0],
            },
            command[1:],
            cwd=cwd,
            timeout=timeout_seconds,
            check=False,
            env=execution_env,
        )
        record["environment"] = str(conda_env)
        if not record["success"]:
            record.setdefault(
                "error",
                (
                    str(record.get("stderr") or "").strip()
                    or str(record.get("stdout") or "").strip()
                    or f"command exited with code {record.get('returncode')}"
                ),
            )
        return record
    except ToolRuntimeError as exc:
        if exc.result:
            return {
                **exc.result,
                "environment": str(conda_env),
                "error": str(exc),
            }
        raise


def _finalize_execution_output_validation(
    execution: dict,
    *,
    success: bool,
    diagnostics: dict,
    error: str | None,
) -> dict:
    """Bind the local execution record to parse/QC success, not exit code alone."""
    record_file = execution.get("record_file")
    record_path = Path(str(record_file)) if record_file else None
    record = dict(execution)
    if record_path is not None and record_path.is_file():
        try:
            loaded = json.loads(record_path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                record = loaded
        except (OSError, json.JSONDecodeError):
            pass
    record["command_success"] = bool(record.get("success"))
    record["output_validation"] = diagnostics
    record["success"] = bool(success)
    if error:
        record["error"] = error
    else:
        record.pop("error", None)
    record["finished_at"] = datetime.now(timezone.utc).isoformat()
    if record_path is not None:
        record_path.write_text(
            json.dumps(record, indent=2, ensure_ascii=False, default=str),
            encoding="utf-8",
        )
        record["record_file"] = str(record_path)
    return record


def _write_execution_provenance(
    run_dir: Path,
    result: dict,
    *,
    provenance: dict,
) -> dict:
    """Persist command/output logs and a compact execution record."""
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    stdout_path.write_text(_stream_text(result.get("stdout")), encoding="utf-8")
    stderr_path.write_text(_stream_text(result.get("stderr")), encoding="utf-8")

    command = result.get("command")
    command_path = run_dir / "command.json"
    command_path.write_text(
        json.dumps(command, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    record = {
        **provenance,
        **{
            key: result.get(key)
            for key in (
                "success",
                "command",
                "returncode",
                "elapsed_seconds",
                "environment",
                "error",
                "process_group_terminated",
            )
            if result.get(key) is not None
        },
        "command_file": str(command_path),
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    record_path = run_dir / "execution.json"
    record_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    record["record_file"] = str(record_path)
    return record


def _stream_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return str(value)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_blueprint_constraints(blueprint_file: str) -> list[dict]:
    """Load Required + Desired interactions from blueprint JSON."""
    try:
        bp = json.loads(Path(blueprint_file).read_text())
        constraints = []
        for item in bp.get("required", []):
            constraints.append({**item, "tier": "required"})
        for item in bp.get("desired", []):
            constraints.append({**item, "tier": "desired"})
        return constraints
    except Exception:
        return []


def _matches_blueprint_pharmacophore(mol_data: dict, constraints: list[dict]) -> bool:
    """Check if molecule has pharmacophore features matching blueprint constraints.

    Simple heuristic: for Required H-bond constraints, molecule must have
    corresponding HBD/HBA count. For Hydrophobic constraints, check LogP > 0.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(mol_data.get("smiles", ""))
    if mol is None:
        return False

    required_constraints = [c for c in constraints if c.get("tier") == "required"]
    if not required_constraints:
        return True  # No required constraints = pass

    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)

    # Count required interaction types
    required_hbond = sum(1 for c in required_constraints if "H-bond" in c.get("interaction_type", ""))
    required_hydrophobic = sum(1 for c in required_constraints if "Hydrophobic" in c.get("interaction_type", ""))

    # Molecule must have enough HBD+HBA to satisfy H-bond requirements
    if required_hbond > 0 and (hbd + hba) < required_hbond:
        return False

    # For hydrophobic requirements, molecule should have some hydrophobic character
    if required_hydrophobic > 0:
        logp = mol_data.get("LogP", 0)
        if logp < -1:
            return False

    return True

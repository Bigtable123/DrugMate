"""GNINA docking utilities.

GNINA is a molecular docking program with CNN scoring. It extends AutoDock Vina
with learned pose scoring (CNN pose score) and affinity prediction (CNN affinity).

Install: download binary from https://github.com/gnina/gnina/releases
"""
from __future__ import annotations

import json
import io
import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from time import monotonic, sleep

from pose_identity import sdf_molblock_sha256

logger = logging.getLogger(__name__)


def check_gnina_available() -> bool:
    return shutil.which("gnina") is not None


def normalize_docking_box(docking_box: dict) -> dict:
    """Validate and normalize an explicit GNINA search box.

    The box is an upstream scientific artifact, not a free-form GNINA config.
    Reject malformed, non-finite, and non-positive dimensions before spawning
    the docking process so an explicit-box run cannot silently search the wrong space.
    """
    if not isinstance(docking_box, dict):
        raise ValueError("docking_box must be a JSON object")

    def triplet(name: str, *, positive: bool = False) -> list[float]:
        raw = docking_box.get(name)
        if not isinstance(raw, (list, tuple)) or len(raw) != 3:
            raise ValueError(f"docking_box.{name} must contain exactly 3 numbers")
        values: list[float] = []
        for index, value in enumerate(raw):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"docking_box.{name}[{index}] must be a number")
            number = float(value)
            if not math.isfinite(number):
                raise ValueError(f"docking_box.{name}[{index}] must be finite")
            if positive and number <= 0:
                raise ValueError(f"docking_box.{name}[{index}] must be > 0")
            values.append(number)
        return values

    unit = docking_box.get("unit", docking_box.get("units", "angstrom"))
    if str(unit).strip().lower() not in {"angstrom", "angstroms", "a", "å"}:
        raise ValueError("docking_box unit must be angstrom")

    normalized = dict(docking_box)
    normalized["center"] = triplet("center")
    normalized["size"] = triplet("size", positive=True)
    normalized["unit"] = "angstrom"
    normalized.pop("units", None)
    return normalized


def load_docking_box(path: Path) -> dict:
    """Load a canonical docking-box JSON artifact and validate its schema."""
    path = Path(path)
    try:
        payload = json.loads(path.read_text())
    except OSError as exc:
        raise ValueError(f"could not read docking_box_file {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"docking_box_file is not valid JSON: {exc}") from exc
    return normalize_docking_box(payload)


def build_search_space_args(
    reference_ligand_file: Path | None,
    docking_box: dict | None,
    autobox_add: float,
) -> list[str]:
    """Build mutually-exclusive GNINA search-space arguments.

    A usable co-crystal ligand retains historical priority.  Only ligand-free
    runs use the explicit pocket-analysis box.
    """
    if reference_ligand_file is not None:
        return [
            "--autobox_ligand", str(reference_ligand_file),
            "--autobox_add", str(autobox_add),
        ]
    box = normalize_docking_box(docking_box)
    center = box["center"]
    size = box["size"]
    return [
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
    ]


def run_gnina_docking(
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    autobox_add: float = 4.0,
    exhaustiveness: int = 16,
    num_modes: int = 5,
    device: int = 0,
    seed: int = 42,
    pose_sort_order: str = "CNNscore",
    cnn_scoring: str = "rescore",
    validate_redocking: bool = True,
    redocking_rmsd_threshold: float = 2.0,
    redocking_rmsd_hard_threshold: float = 2.5,
    gpu_monitor_interval_seconds: int = 5,
    docking_box: dict | None = None,
    timeout_seconds: int = 86400,
    redocking_timeout_seconds: int = 7200,
) -> dict:
    """Run GNINA docking and parse results."""
    receptor_file = Path(receptor_file)
    ligand_file = Path(ligand_file)
    reference_ligand_file = (
        Path(reference_ligand_file) if reference_ligand_file is not None else None
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        redocking_rmsd_threshold = float(redocking_rmsd_threshold)
        redocking_rmsd_hard_threshold = float(
            redocking_rmsd_hard_threshold
        )
    except (TypeError, ValueError):
        redocking_rmsd_threshold = float("nan")
        redocking_rmsd_hard_threshold = float("nan")
    if (
        not math.isfinite(redocking_rmsd_threshold)
        or not math.isfinite(redocking_rmsd_hard_threshold)
        or redocking_rmsd_threshold < 0
        or redocking_rmsd_hard_threshold < redocking_rmsd_threshold
    ):
        return {
            "success": False,
            "error": (
                "redocking RMSD thresholds require 0 <= soft <= hard "
                f"(got soft={redocking_rmsd_threshold}, "
                f"hard={redocking_rmsd_hard_threshold})"
            ),
        }

    try:
        normalized_box = (
            None if reference_ligand_file is not None
            else normalize_docking_box(docking_box)
        )
        search_space_args = build_search_space_args(
            reference_ligand_file, normalized_box, autobox_add,
        )
    except ValueError as exc:
        return {"success": False, "error": str(exc)}

    docked_sdf = output_dir / "docked_poses.sdf"
    log_file = output_dir / "gnina.log"
    gpu_diagnostics_file = output_dir / "gpu_diagnostics.json"
    gpu_usage_file = output_dir / "gpu_usage_timeseries.json"

    cmd = [
        "gnina",
        "--receptor", str(receptor_file),
        "--ligand", str(ligand_file),
        *search_space_args,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--device", str(device),
        "--seed", str(seed),
        "--pose_sort_order", str(pose_sort_order),
        "--cnn_scoring", str(cnn_scoring),
        "--out", str(docked_sdf),
        "--log", str(log_file),
    ]

    gpu_diagnostics = _capture_gpu_diagnostics(device)
    gpu_diagnostics["gnina_command"] = cmd
    gpu_diagnostics_file.write_text(
        json.dumps(gpu_diagnostics, indent=2, ensure_ascii=False)
    )

    logger.info(f"Running GNINA: {' '.join(cmd)}")

    try:
        result = _run_gnina_with_gpu_monitor(
            cmd,
            device=device,
            gpu_usage_file=gpu_usage_file,
            monitor_interval_seconds=gpu_monitor_interval_seconds,
            timeout_seconds=timeout_seconds,
        )
    except TimeoutError:
        gpu_usage = _read_gpu_usage(gpu_usage_file)
        return {
            "success": False,
            "error": f"GNINA timed out after {timeout_seconds}s",
            "gpu_diagnostics": gpu_diagnostics,
            "gpu_diagnostics_file": str(gpu_diagnostics_file),
            "gpu_usage_file": str(gpu_usage_file),
            "gpu_usage_summary": gpu_usage.get("summary"),
        }
    except FileNotFoundError:
        gpu_usage = _read_gpu_usage(gpu_usage_file)
        return {
            "success": False,
            "error": "gnina binary not found",
            "gpu_diagnostics": gpu_diagnostics,
            "gpu_diagnostics_file": str(gpu_diagnostics_file),
            "gpu_usage_file": str(gpu_usage_file),
            "gpu_usage_summary": gpu_usage.get("summary"),
        }

    gpu_usage = _read_gpu_usage(gpu_usage_file)

    if result.returncode != 0:
        return {
            "success": False,
            "error": f"GNINA exited with code {result.returncode}",
            "stderr": result.stderr[:2000],
            "gpu_diagnostics": gpu_diagnostics,
            "gpu_diagnostics_file": str(gpu_diagnostics_file),
            "gpu_usage_file": str(gpu_usage_file),
            "gpu_usage_summary": gpu_usage.get("summary"),
        }

    if not docked_sdf.exists():
        return {
            "success": False,
            "error": "GNINA produced no output SDF",
            "gpu_diagnostics": gpu_diagnostics,
            "gpu_diagnostics_file": str(gpu_diagnostics_file),
            "gpu_usage_file": str(gpu_usage_file),
            "gpu_usage_summary": gpu_usage.get("summary"),
        }

    scores = _parse_gnina_sdf_scores(docked_sdf)
    if not scores:
        return {
            "success": False,
            "error": "GNINA output SDF contained no parseable scored poses",
            "docked_sdf": str(docked_sdf),
            "log_file": str(log_file),
            "gpu_diagnostics": gpu_diagnostics,
            "gpu_diagnostics_file": str(gpu_diagnostics_file),
            "gpu_usage_file": str(gpu_usage_file),
            "gpu_usage_summary": gpu_usage.get("summary"),
        }
    score_coverage = _gnina_score_coverage(
        docked_sdf,
        scores,
        cnn_scoring=cnn_scoring,
    )
    # Incomplete coverage means *some* poses lack a required field, not that the
    # batch failed: the poses that do carry all of them are perfectly usable, and
    # The executable policy removes a missing value from numerator and denominator instead of
    # discarding the row.  Returning ``success: False`` here threw away every
    # scored pose GNINA produced because one pose was short a field, so a single
    # awkward ligand ended docking for the whole cohort (the executable policy scopes a
    # tool fault to the candidate it happened to).  The shortfall is already
    # itemized in ``score_coverage`` and is declared in ``operational_errors`` for
    # the failure registry (reflective-coordinator/SKILL.md:260-266).
    operational_errors = []
    if not score_coverage["complete"]:
        operational_errors.append({
            "tool": "gnina",
            "reason": (
                "GNINA output SDF has incomplete score coverage: "
                f"{score_coverage['num_complete_poses']}/"
                f"{score_coverage['num_valid_poses']} valid poses contain all "
                "required fields "
                f"{', '.join(score_coverage['required_fields'])}"
            ),
        })

    redocking_validation = {
        "metric_id": "redocking.native_ligand_rmsd",
        "engine": "gnina",
        "provider": "GNINA",
        "performed": False,
        "passed": None,
        "applicable": False,
        "not_applicable": True,
        "reason": (
            "disabled_by_config" if not validate_redocking
            else "reference ligand unavailable"
        ),
        "gate": "not_applicable",
        "rmsd_threshold_angstrom": redocking_rmsd_threshold,
        "rmsd_soft_threshold_angstrom": redocking_rmsd_threshold,
        "rmsd_hard_threshold_angstrom": redocking_rmsd_hard_threshold,
    }
    if validate_redocking and reference_ligand_file is not None:
        redocking_validation = _run_reference_redocking(
            receptor_file=receptor_file,
            reference_ligand_file=reference_ligand_file,
            output_dir=output_dir,
            autobox_add=autobox_add,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            device=device,
            seed=seed,
            pose_sort_order=pose_sort_order,
            cnn_scoring=cnn_scoring,
            rmsd_soft_threshold=redocking_rmsd_threshold,
            rmsd_hard_threshold=redocking_rmsd_hard_threshold,
            timeout_seconds=redocking_timeout_seconds,
        )

    # Only a redocking that actually measured an RMSD carries a gate verdict.
    # Disabled, no-reference, crashed and atom-mapping failures all leave
    # ``passed`` as None, and absence of a measurement is not a protocol failure
    # -- keying this on ``performed`` alone charged a crashed run as
    # a hard gate failure.
    if redocking_validation.get("passed") is None:
        quality_gate_passed = None
    else:
        quality_gate_passed = redocking_validation.get("passed") is True

    return {
        "success": True,
        "docked_sdf": str(docked_sdf),
        "log_file": str(log_file),
        "num_poses": len(scores),
        "scores": scores,
        "score_metadata": score_metadata(),
        "score_coverage": score_coverage,
        "operational_errors": operational_errors,
        "redocking_validation": redocking_validation,
        "quality_gate_passed": quality_gate_passed,
        "search_space": (
            {
                "mode": "reference_autobox",
                "reference_ligand_file": str(reference_ligand_file),
                "autobox_add": autobox_add,
            }
            if reference_ligand_file is not None
            else {"mode": "explicit_box", **(normalized_box or {})}
        ),
        "software": {"gnina_version": _gnina_version()},
        "gpu_diagnostics": gpu_diagnostics,
        "gpu_diagnostics_file": str(gpu_diagnostics_file),
        "gpu_usage_file": str(gpu_usage_file),
        "gpu_usage_summary": gpu_usage.get("summary"),
    }


def _gnina_version() -> str | None:
    try:
        result = subprocess.run(["gnina", "--version"], capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return None
    value = (result.stdout or result.stderr).strip()
    return value.splitlines()[0] if value else None


def score_metadata() -> dict:
    """Machine-readable score semantics used by downstream ranking."""
    return {
        "minimizedAffinity": {
            "meaning": "Vina-like minimized energy",
            "unit": "kcal/mol",
            "better": "lower",
        },
        "CNNscore": {
            "meaning": "CNN pose-quality score",
            "unit": "probability-like score",
            "better": "higher",
        },
        "CNNaffinity": {
            "meaning": "CNN-predicted binding affinity on a pK scale",
            "unit": "pK",
            "better": "higher",
        },
        "intramol": {
            "meaning": "Optional GNINA ligand-internal score term",
            "unit": "GNINA score unit",
            "better": "context_only",
        },
        "CNNvariance": {
            "meaning": "Optional GNINA ensemble-model variance",
            "unit": "ensemble variance",
            "better": "context_only",
        },
        "engine_pose_rank": {
            "meaning": "One-based output rank under the configured pose sort order",
            "unit": "rank",
            "better": "context_only",
        },
    }


def _run_reference_redocking(
    receptor_file: Path,
    reference_ligand_file: Path,
    output_dir: Path,
    autobox_add: float,
    exhaustiveness: int,
    num_modes: int,
    device: int,
    seed: int,
    pose_sort_order: str,
    cnn_scoring: str,
    rmsd_soft_threshold: float,
    rmsd_hard_threshold: float,
    timeout_seconds: int,
) -> dict:
    """Redock the co-crystal ligand and compare the best pose with its crystal pose."""
    redocked_sdf = output_dir / "redocked_reference.sdf"
    # A same-directory rerun must never treat a previous reference pose as the
    # output of this subprocess when GNINA fails before opening --out.
    redocked_sdf.unlink(missing_ok=True)
    redocking_log = output_dir / "redocking_reference.log"
    cmd = [
        "gnina",
        "--receptor", str(receptor_file),
        "--ligand", str(reference_ligand_file),
        "--autobox_ligand", str(reference_ligand_file),
        "--autobox_add", str(autobox_add),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--device", str(device),
        "--seed", str(seed),
        "--pose_sort_order", str(pose_sort_order),
        "--cnn_scoring", str(cnn_scoring),
        "--out", str(redocked_sdf),
        "--log", str(redocking_log),
    ]
    base = {
        "metric_id": "redocking.native_ligand_rmsd",
        "engine": "gnina",
        "provider": "GNINA",
        "performed": True,
        "applicable": True,
        "not_applicable": False,
        "rmsd_threshold_angstrom": rmsd_soft_threshold,
        "rmsd_soft_threshold_angstrom": rmsd_soft_threshold,
        "rmsd_hard_threshold_angstrom": rmsd_hard_threshold,
        "redocked_sdf": str(redocked_sdf),
        "log_file": str(redocking_log),
        "command": cmd,
    }
    # A redocking run that crashed, timed out, exited non-zero, or yielded no
    # atom-compatible pose measured no RMSD at all.  Reporting those as
    # ``passed: False`` / ``gate: "reject"`` made them indistinguishable from a
    # genuine RMSD past the hard threshold, and the coordinator reads the gate
    # literally (reflective-coordinator/run.py:1407-1429) -- so a benchmark that
    # never produced a number triggered a protocol rollback. The executable policy drops a
    # missing value from numerator and denominator alike: no verdict, not a
    # failing one.  This is the same tri-state shape the disabled / no-reference
    # branch already publishes (see ``redocking_validation`` above), and it keeps
    # the attempt itself auditable via ``performed`` plus ``reason``.
    unmeasured = {
        **base,
        "passed": None,
        "measured": False,
        "applicable": False,
        "not_applicable": True,
        "gate": "not_applicable",
    }
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            **unmeasured,
            "reason": f"redocking failed within {timeout_seconds}s: {exc}",
        }
    if proc.returncode != 0 or not redocked_sdf.exists():
        return {
            **unmeasured,
            "reason": f"GNINA redocking exited with code {proc.returncode}",
            "stderr": proc.stderr[-2000:],
        }

    try:
        rmsd, pose_index = _best_reference_rmsd(reference_ligand_file, redocked_sdf)
    except Exception as exc:
        return {**unmeasured, "reason": f"RMSD calculation failed: {exc}"}
    if rmsd is None:
        return {
            **unmeasured,
            "reason": "no atom-compatible redocked pose for RMSD",
        }
    gate = _redocking_rmsd_gate(
        rmsd,
        soft_threshold=rmsd_soft_threshold,
        hard_threshold=rmsd_hard_threshold,
    )
    return {
        **base,
        "measured": True,
        # Keep the machine-readable scalar unrounded so the exported metric and
        # the gate computed above cannot disagree at 2.0/2.5 A.  The rounded
        # copy is presentation-only.
        "best_rmsd_angstrom": float(rmsd),
        "best_rmsd_angstrom_display": round(float(rmsd), 3),
        "best_rmsd_pose_index": pose_index,
        "gate": gate,
        "passed": gate == "keep",
        "reason": {
            "keep": "RMSD is within the soft threshold",
            "modify": "RMSD is between the soft and hard thresholds; rerun protocol",
            "reject": "RMSD exceeds the hard threshold; roll back protocol",
        }[gate],
    }


def _redocking_rmsd_gate(
    rmsd: float,
    *,
    soft_threshold: float = 2.0,
    hard_threshold: float = 2.5,
) -> str:
    """Return the protocol gate for a finite native-redocking RMSD."""
    value = float(rmsd)
    soft = float(soft_threshold)
    hard = float(hard_threshold)
    if not all(math.isfinite(item) for item in (value, soft, hard)):
        raise ValueError("redocking RMSD and thresholds must be finite")
    if soft < 0 or hard < soft:
        raise ValueError("redocking RMSD thresholds require 0 <= soft <= hard")
    if value > hard:
        return "reject"
    if value > soft:
        return "modify"
    return "keep"


def _best_reference_rmsd(reference_sdf: Path, redocked_sdf: Path) -> tuple[float | None, int | None]:
    """Return the smallest symmetry-aware heavy-atom RMSD across redocked poses."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    def parse_physical_record(block: str):
        framed = block if block.endswith(("\n", "\r")) else block + "\n"
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO((framed + "$$$$\n").encode("utf-8")),
            removeHs=False,
            strictParsing=False,
        )
        parsed = list(supplier)
        return parsed[0] if len(parsed) == 1 else None

    reference_blocks = _physical_sdf_blocks(reference_sdf)
    if len(reference_blocks) != 1:
        raise ValueError(
            "reference ligand SDF must contain exactly one physical record"
        )
    reference = parse_physical_record(reference_blocks[0])
    if reference is None:
        raise ValueError("reference ligand physical record could not be parsed")
    reference = Chem.RemoveHs(reference)

    best: tuple[float, int] | None = None
    # Enumerate delimiter-framed records, not successful supplier outputs.  A
    # malformed record must retain its physical slot; otherwise the returned
    # index can name different coordinates when the artifact is later sealed.
    for index, block in enumerate(_physical_sdf_blocks(redocked_sdf)):
        pose = parse_physical_record(block)
        if pose is None:
            continue
        pose = Chem.RemoveHs(pose)
        if pose.GetNumAtoms() != reference.GetNumAtoms():
            continue
        try:
            # Docking accuracy is measured in the receptor frame. Aligning the
            # ligand alone would erase an incorrect translation or rotation.
            rmsd = float(rdMolAlign.CalcRMS(pose, reference))
        except Exception:
            continue
        if best is None or rmsd < best[0]:
            best = (rmsd, index)
    return best if best is not None else (None, None)


def _run_gnina_with_gpu_monitor(
    cmd: list[str],
    device: int,
    gpu_usage_file: Path,
    monitor_interval_seconds: int,
    timeout_seconds: int,
) -> subprocess.CompletedProcess:
    """Run GNINA while sampling GPU usage with nvidia-smi."""
    interval = max(1, int(monitor_interval_seconds or 5))
    samples: list[dict] = []
    started_at = datetime.now(timezone.utc).isoformat()

    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stdout_file, tempfile.TemporaryFile(
        mode="w+t",
        encoding="utf-8",
    ) as stderr_file:
        proc = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, text=True)
        deadline = monotonic() + timeout_seconds
        next_sample_at = monotonic()
        timed_out = False

        try:
            while True:
                now = monotonic()
                if now >= next_sample_at:
                    samples.append(_capture_gpu_usage_sample(device, proc.pid))
                    _write_gpu_usage(gpu_usage_file, started_at, samples)
                    next_sample_at = now + interval

                if proc.poll() is not None:
                    break

                if now >= deadline:
                    timed_out = True
                    proc.kill()
                    proc.wait()
                    break

                sleep(min(1.0, max(0.1, next_sample_at - now)))
        finally:
            samples.append(_capture_gpu_usage_sample(device, proc.pid))
            _write_gpu_usage(gpu_usage_file, started_at, samples)

        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read()
        stderr = stderr_file.read()

    if timed_out:
        raise TimeoutError("GNINA timed out")

    return subprocess.CompletedProcess(cmd, proc.returncode, stdout, stderr)


def _capture_gpu_usage_sample(device: int, gnina_pid: int) -> dict:
    sample = {
        "captured_at": datetime.now(timezone.utc).isoformat(),
        "requested_device": device,
        "gnina_pid": gnina_pid,
        "nvidia_smi_available": shutil.which("nvidia-smi") is not None,
    }

    if not sample["nvidia_smi_available"]:
        sample["error"] = "nvidia-smi not found in PATH"
        return sample

    sample["gpu"] = _query_gpu_usage(device)
    sample["compute_processes"] = _query_gpu_processes()

    processes = sample.get("compute_processes") or []
    sample["gnina_process_seen"] = any(
        proc.get("pid") == gnina_pid or "gnina" in proc.get("process_name", "").lower()
        for proc in processes
    )
    return sample


def _query_gpu_usage(device: int) -> dict:
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                f"--id={device}",
                "--query-gpu=timestamp,index,name,utilization.gpu,memory.used,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"error": str(exc)}

    if result.returncode != 0:
        return {"error": result.stderr.strip(), "returncode": result.returncode}

    line = result.stdout.strip().splitlines()[0] if result.stdout.strip() else ""
    parts = [part.strip() for part in line.split(",")]
    if len(parts) < 6:
        return {"raw": line}

    return {
        "timestamp": parts[0],
        "index": _safe_int(parts[1]),
        "name": parts[2],
        "utilization_gpu_percent": _safe_int(parts[3]),
        "memory_used_mib": _safe_int(parts[4]),
        "memory_total_mib": _safe_int(parts[5]),
    }


def _query_gpu_processes() -> list[dict]:
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-compute-apps=pid,process_name,used_gpu_memory",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []

    if result.returncode != 0:
        return []

    processes = []
    for line in result.stdout.strip().splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) < 3:
            continue
        processes.append(
            {
                "pid": _safe_int(parts[0]),
                "process_name": parts[1],
                "used_gpu_memory_mib": _safe_int(parts[2]),
            }
        )
    return processes


def _write_gpu_usage(gpu_usage_file: Path, started_at: str, samples: list[dict]) -> None:
    payload = {
        "started_at": started_at,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(samples),
        "samples": samples,
        "summary": _summarize_gpu_usage(samples),
    }
    gpu_usage_file.write_text(json.dumps(payload, indent=2, ensure_ascii=False))


def _read_gpu_usage(gpu_usage_file: Path) -> dict:
    if not gpu_usage_file.exists():
        return {}
    try:
        return json.loads(gpu_usage_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _summarize_gpu_usage(samples: list[dict]) -> dict:
    utils = []
    memory = []
    process_memory = []
    process_seen = False

    for sample in samples:
        gpu = sample.get("gpu") or {}
        if isinstance(gpu.get("utilization_gpu_percent"), int):
            utils.append(gpu["utilization_gpu_percent"])
        if isinstance(gpu.get("memory_used_mib"), int):
            memory.append(gpu["memory_used_mib"])
        for proc in sample.get("compute_processes") or []:
            if isinstance(proc.get("used_gpu_memory_mib"), int):
                process_memory.append(proc["used_gpu_memory_mib"])
        process_seen = process_seen or bool(sample.get("gnina_process_seen"))

    return {
        "gnina_process_seen": process_seen,
        "max_gpu_util_percent": max(utils) if utils else None,
        "max_gpu_memory_mib": max(memory) if memory else None,
        "max_process_gpu_memory_mib": max(process_memory) if process_memory else None,
    }


def _safe_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _capture_gpu_diagnostics(device: int) -> dict:
    """Capture GPU visibility before GNINA starts."""
    diagnostics = {
        "captured_at": datetime.now(timezone.utc).isoformat(),
        "requested_device": device,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "nvidia_smi_available": shutil.which("nvidia-smi") is not None,
    }

    if not diagnostics["nvidia_smi_available"]:
        diagnostics["nvidia_smi_error"] = "nvidia-smi not found in PATH"
        return diagnostics

    try:
        result = subprocess.run(
            ["nvidia-smi"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        diagnostics["nvidia_smi_error"] = str(exc)
        return diagnostics

    diagnostics["nvidia_smi_returncode"] = result.returncode
    diagnostics["nvidia_smi_stdout"] = result.stdout
    diagnostics["nvidia_smi_stderr"] = result.stderr
    return diagnostics


def _physical_sdf_blocks(sdf_path: Path) -> list[str]:
    """Split SDF only on exact delimiter lines and preserve physical order."""
    blocks: list[str] = []
    current: list[str] = []
    with Path(sdf_path).open(errors="ignore", newline="") as handle:
        for line in handle:
            if line.rstrip("\r\n") == "$$$$":
                # Preserve empty/malformed physical records.  Dropping one
                # would left-shift the pose_index of every later score.
                blocks.append("".join(current))
                current = []
            else:
                current.append(line)
    if current and any(part.strip() for part in current):
        blocks.append("".join(current))
    return blocks


def _parse_gnina_sdf_scores(sdf_path: Path) -> list[dict]:
    """Parse scores per physical block so a bad row cannot shift later IDs."""
    scores = []
    try:
        from rdkit import Chem
        ranks_by_name: dict[str, int] = {}
        for i, block in enumerate(_physical_sdf_blocks(sdf_path)):
            framed = block if block.endswith(("\n", "\r")) else block + "\n"
            payload = (framed + "$$$$\n").encode("utf-8")
            supplier = Chem.ForwardSDMolSupplier(
                io.BytesIO(payload), removeHs=False, strictParsing=False,
            )
            molecules = list(supplier)
            if len(molecules) != 1 or molecules[0] is None:
                continue
            mol = molecules[0]
            mol_name = (
                mol.GetProp("_Name")
                if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
                else f"pose_{i}"
            )
            ranks_by_name[mol_name] = ranks_by_name.get(mol_name, 0) + 1
            entry = {
                "pose_index": i,
                "mol_name": mol_name,
                "engine": "gnina",
                "engine_pose_rank": ranks_by_name[mol_name],
                # GNINA embeds every score in the same physical SDF record as
                # its coordinates.  Preserve a content witness for that exact
                # mol block so downstream joins never have to trust row index.
                "sdf_molblock_sha256": sdf_molblock_sha256(block),
            }
            for prop in (
                "minimizedAffinity", "CNNscore", "CNNaffinity",
                "intramol", "CNNvariance",
            ):
                if mol.HasProp(prop):
                    try:
                        value = float(mol.GetProp(prop))
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(value):
                        entry[prop] = value
            if any(
                prop in entry
                for prop in ("minimizedAffinity", "CNNscore", "CNNaffinity")
            ):
                scores.append(entry)
    except ImportError:
        scores = _parse_gnina_sdf_scores_fallback(sdf_path)
    return scores


def _count_valid_sdf_poses(sdf_path: Path) -> int:
    """Count poses that downstream RDKit-based evaluation can actually read."""
    try:
        from rdkit import Chem
    except ImportError:
        return sum(
            1
            for block in Path(sdf_path).read_text(errors="ignore").split("$$$$")
            if any(line.strip() for line in block.splitlines())
        )
    count = 0
    for block in _physical_sdf_blocks(sdf_path):
        framed = block if block.endswith(("\n", "\r")) else block + "\n"
        payload = (framed + "$$$$\n").encode("utf-8")
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(payload), removeHs=False, strictParsing=False,
        )
        molecules = list(supplier)
        count += int(len(molecules) == 1 and molecules[0] is not None)
    return count


def _gnina_score_coverage(
    sdf_path: Path,
    scores: list[dict],
    *,
    cnn_scoring: str,
) -> dict:
    """Require the configured GNINA primary fields on every valid output pose."""
    required_fields = ["minimizedAffinity"]
    if str(cnn_scoring).casefold() != "none":
        required_fields.extend(["CNNscore", "CNNaffinity"])
    num_valid_poses = _count_valid_sdf_poses(sdf_path)
    num_physical_poses = len(_physical_sdf_blocks(sdf_path))
    missing_fields = []
    for score in scores:
        missing = [field for field in required_fields if field not in score]
        if missing:
            missing_fields.append({
                "pose_index": score.get("pose_index"),
                "mol_name": score.get("mol_name"),
                "missing_fields": missing,
            })
    num_unparsed_poses = max(0, num_valid_poses - len(scores))
    num_complete_poses = max(
        0,
        num_valid_poses - num_unparsed_poses - len(missing_fields),
    )
    return {
        "complete": (
            num_valid_poses == len(scores)
            and not missing_fields
        ),
        "required_fields": required_fields,
        "num_physical_poses": num_physical_poses,
        "num_invalid_poses": num_physical_poses - num_valid_poses,
        "physical_records_complete": num_physical_poses == num_valid_poses,
        "num_valid_poses": num_valid_poses,
        "num_parsed_score_records": len(scores),
        "num_complete_poses": num_complete_poses,
        "num_unparsed_poses": num_unparsed_poses,
        "poses_missing_required_fields": missing_fields,
    }


def _parse_gnina_sdf_scores_fallback(sdf_path: Path) -> list[dict]:
    """Parse SDF without RDKit by reading property blocks."""
    scores = []
    prop_re = re.compile(r"^>\s+<([^>]+)>")
    ranks_by_name: dict[str, int] = {}
    for pose_idx, block in enumerate(_physical_sdf_blocks(sdf_path)):
        lines = block.splitlines()
        if not any(line.strip() for line in lines):
            continue

        mol_name = lines[0].strip() if lines and lines[0].strip() else f"pose_{pose_idx}"
        ranks_by_name[mol_name] = ranks_by_name.get(mol_name, 0) + 1
        current = {
            "pose_index": pose_idx,
            # The SDF title is exactly the first line. Do not search for the
            # first non-empty line: a blank title would otherwise turn the
            # V2000 counts line into a fake compound name.
            "mol_name": mol_name,
            "engine": "gnina",
            "engine_pose_rank": ranks_by_name[mol_name],
            "sdf_molblock_sha256": sdf_molblock_sha256(block),
        }
        i = 0
        while i < len(lines):
            match = prop_re.match(lines[i])
            if not match:
                i += 1
                continue

            prop_name = match.group(1)
            value_line = lines[i + 1].strip() if i + 1 < len(lines) else ""
            if prop_name in (
                "minimizedAffinity", "CNNscore", "CNNaffinity",
                "intramol", "CNNvariance",
            ):
                try:
                    value = float(value_line)
                except ValueError:
                    value = None
                if value is not None and math.isfinite(value):
                    current[prop_name] = value
            i += 2

        if any(prop in current for prop in ("minimizedAffinity", "CNNscore", "CNNaffinity")):
            scores.append(current)
    return scores

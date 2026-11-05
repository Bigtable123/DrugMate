"""Drugmate-side orchestration for the isolated ChemFH worker."""
from __future__ import annotations

import json
from pathlib import Path

from deployment_config import deployment_path
import tool_runtime


WORKER = Path(__file__).resolve().with_name("chemfh_worker.py")
DEFAULT_SOURCE_DIR = deployment_path("${SBDD_CHEMFH_SOURCE_DIR}")
DEFAULT_CHECKPOINT = (
    DEFAULT_SOURCE_DIR
    / "static"
    / "media"
    / "chemprop"
    / "checkpoint"
    / "fold_0"
    / "model_0"
    / "model.pt"
)
DEFAULT_CHECKPOINT_SHA256 = (
    "9f5c82630bed2d787eaf16d99eeac500763748d54474818757207659433a8b6b"
)
DEFAULT_RULE_BUNDLE_SHA256 = (
    "f67d0cf228af15cd3ccf200034f0363a1495aaba0f25388e0393379cd1e55519"
)
DEFAULT_SOURCE_COMMIT = "f8c046fac9a1d331979025a104382d3d2b2614ca"
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "chemfh",
    "command": "python",
    "source_dir": str(DEFAULT_SOURCE_DIR),
    "source_commit": DEFAULT_SOURCE_COMMIT,
    "checkpoint": str(DEFAULT_CHECKPOINT),
    "checkpoint_sha256": DEFAULT_CHECKPOINT_SHA256,
    "rule_bundle_sha256": DEFAULT_RULE_BUNDLE_SHA256,
}


def check_chemfh_available(runtime: dict | None = None) -> bool:
    """Run ChemFH's source-aware functional probe in its registered runtime."""
    candidate = _runtime(runtime)
    candidate.setdefault("probe_timeout_seconds", 180)
    candidate.setdefault("probe_args", _worker_args(candidate, probe=True))
    return bool(tool_runtime.probe_tool(candidate, timeout=180).get("available"))


def run_chemfh_worker(
    candidates: list[dict],
    output_dir: Path,
    *,
    runtime: dict | None = None,
    timeout_seconds: int = 7200,
    seed: int = 42,
) -> dict:
    """Run one strict ChemFH batch and return its structured JSON output."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / "chemfh_worker_input.json"
    output_path = output_dir / "chemfh_worker_output.json"
    input_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "seed": int(seed),
                "candidates": candidates,
            },
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )

    candidate_runtime = _runtime(runtime)
    args = _worker_args(
        candidate_runtime,
        input_path=input_path,
        output_path=output_path,
        seed=seed,
    )
    source_dir = Path(candidate_runtime["source_dir"]).resolve()
    try:
        execution = tool_runtime.run_tool(
            candidate_runtime,
            args,
            cwd=source_dir,
            timeout=max(1, int(timeout_seconds)),
            check=False,
            env={
                "CUDA_VISIBLE_DEVICES": "-1",
                "DJANGO_SETTINGS_MODULE": "ChemFH.settings",
                "PYTHONDONTWRITEBYTECODE": "1",
            },
        )
    except (OSError, tool_runtime.ToolRuntimeError) as exc:
        execution = {
            "success": False,
            "error": f"could not execute ChemFH runtime: {exc}",
            "environment": candidate_runtime.get("environment"),
        }
    result = _read_worker_output(output_path, execution)
    result["files"] = {
        "input_json": str(input_path),
        "output_json": str(output_path),
    }
    return result


def _runtime(runtime: dict | None) -> dict:
    candidate = dict(DEFAULT_RUNTIME)
    candidate.update(runtime or {})
    candidate.setdefault("kind", "command")
    candidate.setdefault("environment", "chemfh")
    candidate.setdefault("command", "python")
    candidate.setdefault("source_dir", str(DEFAULT_SOURCE_DIR))
    candidate.setdefault("source_commit", DEFAULT_SOURCE_COMMIT)
    candidate.setdefault("checkpoint", str(DEFAULT_CHECKPOINT))
    candidate.setdefault("checkpoint_sha256", DEFAULT_CHECKPOINT_SHA256)
    candidate.setdefault("rule_bundle_sha256", DEFAULT_RULE_BUNDLE_SHA256)
    return candidate


def _worker_args(
    runtime: dict,
    *,
    input_path: Path | None = None,
    output_path: Path | None = None,
    seed: int = 42,
    probe: bool = False,
) -> list[str]:
    args = [
        "-B",
        str(WORKER),
        "--source-dir", str(runtime["source_dir"]),
        "--checkpoint", str(runtime["checkpoint"]),
        "--expected-checkpoint-sha256", str(runtime["checkpoint_sha256"]),
        "--expected-rule-bundle-sha256", str(runtime["rule_bundle_sha256"]),
        "--expected-source-commit", str(runtime["source_commit"]),
        "--seed", str(int(seed)),
    ]
    if probe:
        return [*args, "--probe"]
    if input_path is None or output_path is None:
        raise ValueError("ChemFH worker execution needs input and output paths")
    return [
        *args,
        "--input-file", str(input_path),
        "--output-file", str(output_path),
    ]


def _read_worker_output(output_path: Path, execution: dict) -> dict:
    if output_path.is_file():
        try:
            result = json.loads(output_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            result = {
                "success": False,
                "error": f"invalid ChemFH worker output: {exc}",
            }
    else:
        result = {
            "success": False,
            "error": execution.get("error") or "ChemFH worker wrote no output",
        }
    if not execution.get("success"):
        result["success"] = False
        result.setdefault(
            "error",
            execution.get("error")
            or str(execution.get("stderr") or "").strip()
            or str(execution.get("stdout") or "").strip()
            or "ChemFH worker failed",
        )
    result["execution"] = {
        "success": bool(execution.get("success")),
        "command": execution.get("command"),
        "returncode": execution.get("returncode"),
        "elapsed_seconds": execution.get("elapsed_seconds"),
        "environment": execution.get("environment"),
        "stdout_tail": str(execution.get("stdout") or "")[-4000:],
        "stderr_tail": str(execution.get("stderr") or "")[-4000:],
    }
    return result

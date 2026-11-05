"""Drugmate-side launcher for the isolated DeepPocket worker."""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

from deployment_config import deployment_path
from tool_runtime import run_tool


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BRIDGE = SCRIPT_DIR / "deeppocket_bridge.py"
DEFAULT_SOURCE_DIR = deployment_path("${SBDD_DEEPPOCKET_SOURCE_DIR}")
DEFAULT_CLASS_CHECKPOINT = (
    DEFAULT_SOURCE_DIR / "models" / "first_model_fold1_best_test_auc_85001.pth.tar"
)
DEFAULT_SEGMENTATION_CHECKPOINT = (
    DEFAULT_SOURCE_DIR / "models" / "seg0_best_test_IOU_91.pth.tar"
)
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "deeppocket",
    "command": "python",
    "bridge_path": str(DEFAULT_BRIDGE),
    "source_dir": str(DEFAULT_SOURCE_DIR),
    "class_checkpoint": str(DEFAULT_CLASS_CHECKPOINT),
    "segmentation_checkpoint": str(DEFAULT_SEGMENTATION_CHECKPOINT),
}


def _diagnostic_text(execution: dict) -> str:
    detail = (
        execution.get("error")
        or str(execution.get("stderr") or "").strip()
        or str(execution.get("stdout") or "").strip()
        or "DeepPocket worker failed without diagnostic output"
    )
    return " ".join(str(detail).split())[-4000:]


def run_deeppocket(
    protein_file: Path,
    output_json: Path,
    *,
    runtime: dict | None = None,
    params: dict | None = None,
) -> dict[str, Any]:
    """Execute DeepPocket and return its validated, normalized candidates."""
    runtime_spec = {**DEFAULT_RUNTIME, **dict(runtime or {})}
    settings = dict(params or {})
    bridge = runtime_spec.get("bridge_path") or str(DEFAULT_BRIDGE)
    source_dir = runtime_spec.get("source_dir") or str(DEFAULT_SOURCE_DIR)
    class_checkpoint = (
        runtime_spec.get("class_checkpoint") or str(DEFAULT_CLASS_CHECKPOINT)
    )
    segmentation_checkpoint = (
        runtime_spec.get("segmentation_checkpoint")
        or str(DEFAULT_SEGMENTATION_CHECKPOINT)
    )
    timeout = int(settings.get("timeout_seconds", 1200))
    raw_output_dir = output_json.parent / "deeppocket_output"

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.unlink(missing_ok=True)
    # A forced rerun may reuse a Step-1 directory.  Remove only this tool's
    # dedicated raw directory so a failed new run cannot expose stale masks or
    # segmented residues from an earlier successful invocation.
    if raw_output_dir.is_dir():
        shutil.rmtree(raw_output_dir)
    arguments = [
        bridge,
        "predict",
        "--source-dir", source_dir,
        "--class-checkpoint", class_checkpoint,
        "--segmentation-checkpoint", segmentation_checkpoint,
        "--protein-file", str(Path(protein_file).resolve()),
        "--raw-output-dir", str(raw_output_dir.resolve()),
        "--output-json", str(output_json.resolve()),
        "--max-candidates", str(settings.get("max_candidates", 3)),
        "--threshold", str(settings.get("segmentation_threshold", 0.5)),
        "--mask-distance", str(settings.get("mask_distance", 3.5)),
    ]
    execution = run_tool(
        runtime_spec,
        arguments,
        timeout=timeout,
        check=False,
    )

    payload: dict[str, Any] | None = None
    if output_json.is_file():
        try:
            candidate = json.loads(output_json.read_text(encoding="utf-8"))
            if isinstance(candidate, dict):
                payload = candidate
        except (OSError, json.JSONDecodeError):
            payload = None

    execution_summary = {
        "success": bool(execution.get("success")),
        "returncode": execution.get("returncode"),
        "elapsed_seconds": execution.get("elapsed_seconds"),
        "environment": execution.get("environment"),
        "command": execution.get("command"),
    }
    if payload is None:
        payload = {
            "schema_version": 1,
            "success": False,
            "mode": "predict",
            "error_code": "worker_output_missing_or_invalid",
            "error": _diagnostic_text(execution),
        }
    payload["worker_execution"] = execution_summary

    # Require agreement between process state and the parsed scientific
    # contract.  This prevents a stale success JSON from masking a failed run.
    if not execution.get("success"):
        payload["success"] = False
        payload.setdefault("error_code", "deeppocket_worker_failed")
        payload.setdefault("error", _diagnostic_text(execution))
    elif payload.get("success") is not True:
        payload["success"] = False
        payload.setdefault("error_code", "deeppocket_contract_failed")
        payload.setdefault(
            "error", "DeepPocket worker returned an unsuccessful scientific contract"
        )
    return payload

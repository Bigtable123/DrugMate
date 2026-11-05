"""Drugmate-side launcher for the isolated PocketMiner worker.

TensorFlow, MDTraj, and the upstream GVP source are intentionally never imported
here.  The registry-owned runtime executes ``pocketminer_bridge.py`` inside the
dedicated ``pocketminer`` Conda environment and returns its JSON contract.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from deployment_config import deployment_path
from tool_runtime import run_tool


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BRIDGE = SCRIPT_DIR / "pocketminer_bridge.py"
DEFAULT_SOURCE_DIR = deployment_path("${SBDD_POCKETMINER_SOURCE_DIR}")
DEFAULT_CHECKPOINT = DEFAULT_SOURCE_DIR / "models" / "pocketminer"
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "pocketminer",
    "command": "python",
    "bridge_path": str(DEFAULT_BRIDGE),
    "source_dir": str(DEFAULT_SOURCE_DIR),
    "checkpoint": str(DEFAULT_CHECKPOINT),
}


def _diagnostic_text(execution: dict) -> str:
    detail = (
        execution.get("error")
        or str(execution.get("stderr") or "").strip()
        or str(execution.get("stdout") or "").strip()
        or "PocketMiner worker failed without diagnostic output"
    )
    return " ".join(str(detail).split())[-4000:]


def run_pocketminer(
    protein_file: Path,
    chain: str,
    output_json: Path,
    *,
    runtime: dict | None = None,
    params: dict | None = None,
) -> dict[str, Any]:
    """Run strict PocketMiner inference and return the normalized JSON payload."""
    runtime_spec = {**DEFAULT_RUNTIME, **dict(runtime or {})}
    settings = dict(params or {})
    bridge = runtime_spec.get("bridge_path") or str(DEFAULT_BRIDGE)
    source_dir = runtime_spec.get("source_dir") or str(DEFAULT_SOURCE_DIR)
    checkpoint = runtime_spec.get("checkpoint") or str(DEFAULT_CHECKPOINT)
    timeout = int(settings.get("timeout_seconds", 600))

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.unlink(missing_ok=True)
    arguments = [
        bridge,
        "predict",
        "--source-dir", source_dir,
        "--checkpoint", checkpoint,
        "--protein-file", str(Path(protein_file).resolve()),
        "--chain", str(chain),
        "--output-json", str(output_json.resolve()),
        "--provisional-quantile", str(settings.get("provisional_quantile", 0.9)),
        "--cluster-ca-cutoff", str(settings.get("cluster_ca_cutoff", 8.0)),
        "--min-cluster-residues", str(settings.get("min_cluster_residues", 3)),
        "--max-hypotheses", str(settings.get("max_hypotheses", 5)),
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
    payload["execution"] = execution_summary

    # Both layers must agree.  A stale/forged success JSON cannot override a
    # non-zero worker exit, and an exit-zero worker must still assert success.
    if not execution.get("success"):
        payload["success"] = False
        payload.setdefault("error_code", "pocketminer_worker_failed")
        payload.setdefault("error", _diagnostic_text(execution))
    elif payload.get("success") is not True:
        payload["success"] = False
        payload.setdefault("error_code", "pocketminer_contract_failed")
        payload.setdefault("error", "PocketMiner worker returned an unsuccessful contract")
    return payload

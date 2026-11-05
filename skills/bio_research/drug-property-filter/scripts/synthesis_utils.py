"""Drugmate-side orchestration for the isolated AiZynthFinder worker."""
from __future__ import annotations

import json
from pathlib import Path

from deployment_config import deployment_path
import tool_runtime


WORKER = Path(__file__).resolve().with_name("aizynth_worker.py")
DEFAULT_CONFIG = deployment_path("${SBDD_AIZYNTHFINDER_CONFIG}")
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "aizynth",
    "command": "python",
    "config_path": str(DEFAULT_CONFIG),
}


def run_aizynth_worker(
    candidates: list[dict],
    output_dir: Path,
    *,
    runtime: dict | None = None,
    config_path: str | Path | None = None,
    timeout_seconds: int = 21600,
    search_time_limit_seconds: int = 60,
    iteration_limit: int = 100,
    max_routes: int = 5,
    return_first: bool = True,
    expansion_policies: list[str] | None = None,
    filter_policies: list[str] | None = None,
    stocks: list[str] | None = None,
) -> dict:
    """Plan a batch of targets and preserve one output row per candidate id."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / "aizynth_worker_input.json"
    output_path = output_dir / "aizynth_worker_output.json"
    candidate_runtime = dict(runtime or DEFAULT_RUNTIME)
    candidate_runtime.setdefault("kind", "command")
    candidate_runtime.setdefault("environment", "aizynth")
    candidate_runtime.setdefault("command", "python")
    resolved_config = Path(
        config_path
        or candidate_runtime.get("config_path")
        or DEFAULT_CONFIG
    ).expanduser().resolve()
    payload = {
        "schema_version": 1,
        "config_path": str(resolved_config),
        "search_time_limit_seconds": max(1, int(search_time_limit_seconds)),
        "iteration_limit": max(1, int(iteration_limit)),
        "max_routes": max(1, int(max_routes)),
        "return_first": bool(return_first),
        "expansion_policies": list(expansion_policies or ["uspto"]),
        "filter_policies": list(filter_policies or ["uspto"]),
        "stocks": list(stocks or ["zinc"]),
        "candidates": candidates,
    }
    input_path.write_text(
        json.dumps(
            payload,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )

    execution = tool_runtime.run_tool(
        candidate_runtime,
        [
            str(WORKER),
            "--input-file", str(input_path),
            "--output-file", str(output_path),
        ],
        timeout=max(1, int(timeout_seconds)),
        check=False,
    )
    result = _read_worker_output(output_path, execution)
    result["files"] = {
        "input_json": str(input_path),
        "output_json": str(output_path),
    }
    return result


def _read_worker_output(output_path: Path, execution: dict) -> dict:
    if output_path.is_file():
        try:
            result = json.loads(output_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            result = {
                "success": False,
                "error": f"invalid aizynthfinder worker output: {exc}",
            }
    else:
        result = {
            "success": False,
            "error": execution.get("error") or "aizynthfinder worker wrote no output",
        }
    if not execution.get("success"):
        result["success"] = False
        result.setdefault(
            "error",
            execution.get("error")
            or execution.get("stderr", "").strip()
            or execution.get("stdout", "").strip()
            or "aizynthfinder worker failed",
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

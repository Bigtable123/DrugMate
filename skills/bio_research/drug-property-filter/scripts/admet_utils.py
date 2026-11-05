"""Drugmate-side orchestration for the isolated ADMET-AI worker."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import tool_runtime


WORKER = Path(__file__).resolve().with_name("admet_worker.py")
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "admet_ai",
    "command": "python",
}


def check_admet_available(runtime: dict | None = None) -> bool:
    """Probe the registered isolated backend without importing it in Drugmate."""
    candidate = dict(runtime or DEFAULT_RUNTIME)
    candidate.setdefault("kind", "command")
    candidate.setdefault("command", "python")
    candidate.setdefault(
        "probe_code",
        (
            "import os;"
            "os.environ['CUDA_VISIBLE_DEVICES']='-1';"
            "from admet_ai import ADMETModel;"
            "ADMETModel()"
        ),
    )
    return bool(tool_runtime.probe_tool(candidate, timeout=60).get("available"))


def run_admet_worker(
    candidates: list[dict],
    output_dir: Path,
    *,
    runtime: dict | None = None,
    timeout_seconds: int = 7200,
) -> dict:
    """Run one strict batch prediction and return its structured JSON output."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / "admet_worker_input.json"
    output_path = output_dir / "admet_worker_output.json"
    input_path.write_text(
        json.dumps(
            {"schema_version": 1, "candidates": candidates},
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )

    candidate_runtime = dict(runtime or DEFAULT_RUNTIME)
    candidate_runtime.setdefault("kind", "command")
    candidate_runtime.setdefault("environment", "admet_ai")
    candidate_runtime.setdefault("command", "python")
    execution = tool_runtime.run_tool(
        candidate_runtime,
        [
            str(WORKER),
            "--input-file", str(input_path),
            "--output-file", str(output_path),
        ],
        timeout=max(1, int(timeout_seconds)),
        check=False,
        env={"CUDA_VISIBLE_DEVICES": "-1"},
    )
    result = _read_worker_output(output_path, execution, backend_id="admet_ai")
    result["files"] = {
        "input_json": str(input_path),
        "output_json": str(output_path),
    }
    return result


def predict_admet(smiles_list: list[str]) -> list[dict]:
    """Backward-compatible standalone helper using the isolated batch worker."""
    candidates = [
        {"candidate_id": f"candidate-{index}", "smiles": smiles}
        for index, smiles in enumerate(smiles_list)
    ]
    with tempfile.TemporaryDirectory(prefix="admet-ai-") as tmpdir:
        result = run_admet_worker(candidates, Path(tmpdir))
    if not result.get("success"):
        return [{"status": "error", "error": result.get("error")} for _ in candidates]
    by_id = {entry["candidate_id"]: entry for entry in result.get("results") or []}
    return [by_id.get(entry["candidate_id"], {
        "status": "error",
        "error": "ADMET-AI worker omitted candidate",
    }) for entry in candidates]


def _read_worker_output(
    output_path: Path,
    execution: dict,
    *,
    backend_id: str,
) -> dict:
    if output_path.is_file():
        try:
            result = json.loads(output_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            result = {
                "success": False,
                "error": f"invalid {backend_id} worker output: {exc}",
            }
    else:
        result = {
            "success": False,
            "error": execution.get("error") or f"{backend_id} worker wrote no output",
        }
    if not execution.get("success"):
        result["success"] = False
        result.setdefault(
            "error",
            execution.get("error")
            or execution.get("stderr", "").strip()
            or execution.get("stdout", "").strip()
            or f"{backend_id} worker failed",
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

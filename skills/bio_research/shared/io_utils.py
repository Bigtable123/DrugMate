"""Shared I/O utilities for SBDD skills."""
from __future__ import annotations

import argparse
import json
import warnings
import re
import sys
from pathlib import Path
from typing import Any


def skill_result_summary(result: dict) -> dict:
    """Return control fields; scientific rows remain in the stage report files.

    These fields serve artifact wiring, scientific outcomes and bundle checks.
    In-process run() callers still receive the full result; CLI stdout and the
    runner's result.json do not duplicate the scientific evidence.
    """
    fields = (
        "success", "error", "error_code", "error_details", "outcome", "partial",
        "files", "summary", "round_summary", "assessment_schema_version",
        "round_assessment", "decision", "decision_rationale", "decision_outcome", "goal_met",
        "valid_for_nomination", "quality_gate_status", "quality_gate_passed",
        "rollback_instructions", "metadata", "docking_box", "pocket_definition",
        "plots_status", "mode", "candidate_count", "steps_completed",
        "steps_skipped", "publication_status", "main_viewer_error", "pymol_cover_error",
        "structure_2d_status", "cards_status", "missing_sources", "asset_path_warning",
    )
    summary = {key: result[key] for key in fields if key in result}
    summary["result_schema"] = "sbdd-step-result/v1"
    # provided_box routing and legacy scientific outcomes need these scalars.
    # Full tool parameters remain in the stage report.
    if isinstance(result.get("parameters"), dict):
        summary["parameters"] = {
            key: value for key, value in result["parameters"].items()
            if value is None or isinstance(value, (str, bool, int, float))
        }
    assessment = result.get("round_assessment")
    if isinstance(assessment, dict):
        summary["round_assessment"] = {
            key: value for key, value in assessment.items()
            if key != "review_required_molecules"
        }
        if "review_required_molecules" in assessment:
            summary["round_assessment"]["review_required_molecule_count"] = len(
                assessment["review_required_molecules"] or []
            )
    if isinstance(result.get("round_summary"), dict):
        summary["round_summary"] = {
            key: value for key, value in result["round_summary"].items()
            if key != "round_assessment"
        }
    return summary


def print_skill_result(result: dict) -> None:
    """Emit one compact JSON object with incremental writes under python -u."""
    json.dump(skill_result_summary(result), sys.stdout, ensure_ascii=False,
              allow_nan=False, default=str)
    sys.stdout.write("\n")


def number(value, default: float) -> float:
    """Safe float coercion with fallback."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def slugify(value: Any) -> str:
    """Convert value to URL-safe slug."""
    text = str(value or "sbdd").strip().lower()
    text = re.sub(r"[^a-z0-9]+", "-", text).strip("-")
    return text or "sbdd"


def load_scores(path_value: str | None) -> list[dict]:
    """Load docking scores JSON file, returning list of score dicts."""
    if not path_value:
        return []
    path = Path(path_value)
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return []
    if isinstance(data, list):
        return data
    return data.get("scores", []) if isinstance(data, dict) else []


def results_by_pose_index(results: list[dict] | None) -> dict[int, dict]:
    """Index a list of result dicts by their pose_index field."""
    indexed = {}
    for position, result in enumerate(results or []):
        try:
            pose_index = int(result.get("pose_index", position))
        except (TypeError, ValueError):
            pose_index = position
        if pose_index in indexed:
            warnings.warn(
                f"Duplicate pose_index {pose_index}; keeping first result and ignoring duplicate",
                RuntimeWarning,
                stacklevel=2,
            )
            continue
        indexed[pose_index] = result
    return indexed


def load_smiles_map(manifest_path: str | None) -> dict[str, str]:
    """Load compound_id → SMILES mapping from library_manifest.csv."""
    if not manifest_path:
        return {}
    path = Path(manifest_path)
    if not path.exists():
        return {}
    import csv
    mapping = {}
    try:
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                cid = row.get("compound_id", "").strip()
                smi = row.get("smiles", "").strip()
                if cid and smi:
                    mapping[cid] = smi
    except (OSError, KeyError):
        pass
    return mapping


def tool_enabled(input_data: dict, tool_id: str, default: bool = True) -> bool:
    """Resolve one interface-v2 tool switch from a skill's input.

    The orchestrator passes a `tools` block: {tool_id: {"enabled": bool}}. When
    the block is absent (standalone skill invocation) the tool uses its default.
    A tool the caller explicitly disables is skipped without penalty.
    """
    tools = input_data.get("tools")
    if isinstance(tools, dict) and tool_id in tools:
        return bool((tools[tool_id] or {}).get("enabled", default))
    return default


def tool_switches(input_data: dict, tool_ids, default: bool = True) -> dict:
    """Resolve a set of interface-v2 tool switches to {tool_id: enabled}."""
    return {tid: tool_enabled(input_data, tid, default) for tid in tool_ids}


def tool_config(input_data: dict, tool_id: str) -> dict:
    """Return one resolved interface-v2 tool config (params + runtime)."""
    tools = input_data.get("tools")
    if isinstance(tools, dict) and isinstance(tools.get(tool_id), dict):
        return dict(tools[tool_id])
    return {}


def parse_skill_args(description: str) -> tuple[dict, Path]:
    """Unified CLI argument parsing for skill entry points.

    Returns:
        (input_data, output_dir)
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", type=str, help="JSON input string")
    parser.add_argument("--input-file", type=str, help="Path to JSON input file")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    args = parser.parse_args()

    if args.input:
        input_data = json.loads(args.input)
    elif args.input_file:
        with open(args.input_file) as f:
            input_data = json.load(f)
    else:
        parser.error("Must provide --input or --input-file")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return input_data, output_dir

"""SBDD pipeline config materialization — explicit selection + param defaults.

Pipeline:
    explicit user selection  →  merge parameter defaults  →  resolved plan

`resolve()` returns a fully-expanded run plan: every enabled step/tool listed with
all parameter defaults filled in, plus the ordered list of steps to execute.  Step
and tool enablement is never inferred from facts, profiles, or registry defaults:
the input config must state the complete boolean selection matrix.
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path
from typing import Any

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import pipeline_registry as R
import pipeline_timeouts as T


def default_template() -> dict:
    """Build a complete registry-recommendation template.

    This helper is useful for policy comparison and tests.  ``merge_user_config``
    never uses it: Agent-authored configs must provide every selection switch,
    and execution materialization has a separate parameter-only default base.
    """
    steps: dict[str, dict] = {}
    for step_id in R.steps_in_order():
        spec = R.STEPS[step_id]
        tools = {}
        for tool_id, tool in spec["tools"].items():
            tools[tool_id] = {
                "enabled": bool(tool["default_enabled"]),
                "params": _param_defaults(tool.get("params", {})),
            }
        steps[step_id] = {
            "enabled": bool(spec["default_enabled"]),
            "tools": tools,
            "params": _param_defaults(spec.get("params", {})),
        }
    return {
        "project": _param_defaults(R.PROJECT_SCHEMA),
        "context": _param_defaults(R.CONTEXT_SCHEMA),
        "steps": steps,
    }


def _materialization_defaults() -> dict:
    """Build schema defaults without any execution-selection switches.

    Registry ``default_enabled`` values are recommendation metadata.  Keeping
    them out of the merge base makes it structurally impossible for config
    materialization to obtain a step/tool choice from those recommendations.
    """
    steps: dict[str, dict] = {}
    for step_id in R.steps_in_order():
        spec = R.STEPS[step_id]
        steps[step_id] = {
            "params": _param_defaults(spec.get("params", {})),
            "tools": {
                tool_id: {
                    "params": _param_defaults(tool.get("params", {})),
                }
                for tool_id, tool in spec["tools"].items()
            },
        }
    return {
        "project": _param_defaults(R.PROJECT_SCHEMA),
        "context": _param_defaults(R.CONTEXT_SCHEMA),
        "steps": steps,
    }


def _param_defaults(param_schema: dict) -> dict:
    """Return {name: default} for every param that has a non-None default."""
    out = {}
    for name, spec in param_schema.items():
        if spec.get("default") is not None:
            out[name] = copy.deepcopy(spec["default"])
    return out


def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into a deep copy of base (override wins)."""
    merged = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def apply_step_range(config: dict, start_step: Any = None, stop_step: Any = None) -> dict:
    """Translate start_step/stop_step convenience bounds into per-step enabled flags.

    Steps whose registry order falls outside [start, stop] are disabled. Steps
    inside the range keep whatever `enabled` the merged config already has.
    """
    if start_step is None and stop_step is None:
        return config
    try:
        lo = float(start_step) if start_step is not None else float("-inf")
    except (TypeError, ValueError):
        raise ValueError(
            f"start_step must be a number (got {start_step!r})"
        )
    try:
        hi = float(stop_step) if stop_step is not None else float("inf")
    except (TypeError, ValueError):
        raise ValueError(
            f"stop_step must be a number (got {stop_step!r})"
        )
    out = copy.deepcopy(config)
    for step_id, spec in R.STEPS.items():
        order = spec["order"]
        if order < lo or order > hi:
            out.setdefault("steps", {}).setdefault(step_id, {})["enabled"] = False
    return out


def explicit_selection_errors(user_config: dict) -> list[str]:
    """Return omissions in the required explicit step/tool selection matrix."""
    errors: list[str] = []
    if not isinstance(user_config, dict):
        return ["pipeline config must be an object"]
    steps = user_config.get("steps")
    if not isinstance(steps, dict):
        return ["steps must be an object containing the complete explicit selection matrix"]

    for step_id, step_spec in R.STEPS.items():
        step = steps.get(step_id)
        if not isinstance(step, dict):
            errors.append(f"steps.{step_id} must be present as an object")
            continue
        if "enabled" not in step:
            errors.append(f"steps.{step_id}.enabled must be explicitly set")
        tools = step.get("tools")
        if not isinstance(tools, dict):
            errors.append(
                f"steps.{step_id}.tools must contain every registered tool"
            )
            continue
        for tool_id in step_spec["tools"]:
            tool = tools.get(tool_id)
            if not isinstance(tool, dict):
                errors.append(
                    f"steps.{step_id}.tools.{tool_id} must be present as an object"
                )
            elif "enabled" not in tool:
                errors.append(
                    f"steps.{step_id}.tools.{tool_id}.enabled must be explicitly set"
                )
    return errors


def merge_user_config(user_config: dict) -> dict:
    """Validate explicit selection, then fill only registered parameter defaults."""
    user = user_config if user_config is not None else {}
    selection_errors = explicit_selection_errors(user)
    if selection_errors:
        body = "\n".join(f"  - {error}" for error in selection_errors)
        raise ValueError(
            "pipeline config requires an explicit step/tool selection matrix:\n"
            f"{body}"
        )

    # The merge base contains no enable flags: only schema-owned
    # parameter/project/context defaults may survive.
    merged = deep_merge(_materialization_defaults(), user)

    start_step = user.get("start_step")
    stop_step = user.get("stop_step")
    merged = apply_step_range(merged, start_step, stop_step)
    return merged


def resolve(merged_config: dict) -> dict:
    """Build the fully-expanded run plan from a merged (and validated) config.

    Returns:
        {
          "project": {...}, "context": {...},
          "execution_order": ["structure_prep", ...],   # enabled only
          "steps": { step_id: {enabled, skill, order, best_effort,
                               consumes, produces, active_produces,
                               tools: {tool_id: {enabled, params}},
                               tools_enabled: [...], tools_disabled: [...],
                               params: {...}} }
        }
    """
    plan_steps: dict[str, dict] = {}
    execution_order: list[str] = []

    for step_id in R.steps_in_order():
        spec = R.STEPS[step_id]
        cfg_step = (merged_config.get("steps") or {})[step_id]
        enabled = bool(cfg_step["enabled"])

        tools_plan: dict[str, dict] = {}
        tools_enabled: list[str] = []
        tools_disabled: list[str] = []
        for tool_id, tool_spec in spec["tools"].items():
            cfg_tool = (cfg_step.get("tools") or {})[tool_id]
            t_enabled = bool(cfg_tool["enabled"])
            params = deep_merge(
                _param_defaults(tool_spec.get("params", {})),
                cfg_tool.get("params", {}),
            )
            tools_plan[tool_id] = {
                "enabled": t_enabled,
                "params": params,
                "runtime": copy.deepcopy(tool_spec.get("runtime")),
            }
            (tools_enabled if t_enabled else tools_disabled).append(tool_id)

        step_params = deep_merge(
            _param_defaults(spec.get("params", {})),
            cfg_step.get("params", {}),
        )

        conditional_produces = {
            artifact
            for tool_spec in spec["tools"].values()
            for artifact in tool_spec.get("produces", [])
        }
        enabled_tool_produces = {
            artifact
            for tool_id in tools_enabled
            for artifact in spec["tools"][tool_id].get("produces", [])
        }
        active_produces = [
            artifact
            for artifact in spec["produces"]
            if enabled and (
                artifact not in conditional_produces
                or artifact in enabled_tool_produces
            )
        ]
        if step_id == "pocket_analysis" and step_params.get("site_selection_strategy") == "provided_box":
            active_produces = [name for name in active_produces if name != "docking_box"]

        plan_steps[step_id] = {
            "enabled": enabled,
            "skill": spec["skill"],
            "runtime": copy.deepcopy(spec.get("runtime")),
            "order": spec["order"],
            "best_effort": bool(spec.get("best_effort", False)),
            "consumes": list(spec["consumes"]),
            "optional_consumes": list(spec.get("optional_consumes", [])),
            "produces": list(spec["produces"]),
            "active_produces": active_produces,
            "tools": tools_plan,
            "tools_enabled": tools_enabled,
            "tools_disabled": tools_disabled,
            "params": step_params,
            "timeout": T.timeout_budget_breakdown(step_id, merged_config),
        }
        plan_steps[step_id]["timeout"]["resolved_seconds"] = (
            T.resolve_step_timeout(step_id, merged_config)
            if enabled else None
        )
        if enabled:
            execution_order.append(step_id)

    return {
        "project": copy.deepcopy(merged_config.get("project", {})),
        "context": copy.deepcopy(merged_config.get("context", {})),
        "input_artifacts": copy.deepcopy(merged_config.get("input_artifacts", {})),
        "execution_order": execution_order,
        "steps": plan_steps,
    }

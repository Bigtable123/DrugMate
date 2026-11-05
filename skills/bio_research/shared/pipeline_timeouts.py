"""Deterministic outer timeout budgets for the SBDD pipeline.

Budgets cover active execution only: the host resource lease must be acquired
before starting a step clock. Summing backend hard caps remains a conservative
bound when individual runners parallelize work; headroom covers preparation
and merging without allowing smaller legacy defaults to terminate them first.
"""
from __future__ import annotations

from typing import Any

from constants import DEFAULT_STEP_TIMEOUTS


MAX_STEP_TIMEOUT_SECONDS = 30 * 24 * 60 * 60
FIXED_OVERHEAD_SECONDS = 1800
HEADROOM_RATIO = 0.10
_HEADROOM_PERCENT = 10

_POCKET_NATIVE_TOOLS = {"plip", "prolif", "arpeggio", "luna", "oddt"}


class TimeoutBudgetError(ValueError):
    """Raised when a resolved automatic timeout exceeds the safety ceiling."""


def _positive_int(value: Any, default: int = 0) -> int:
    if type(value) is not int:
        return default
    return value if value > 0 else default


def _step_config(config: dict, step_id: str) -> dict:
    steps = config.get("steps")
    if not isinstance(steps, dict):
        return {}
    step = steps.get(step_id)
    return step if isinstance(step, dict) else {}


def _enabled_tools(step: dict) -> dict[str, dict]:
    tools = step.get("tools")
    if not isinstance(tools, dict):
        return {}
    return {
        tool_id: tool
        for tool_id, tool in tools.items()
        if isinstance(tool, dict) and tool.get("enabled") is True
    }


def timeout_budget_breakdown(step_id: str, config: dict) -> dict[str, Any]:
    """Return the serial inner budget and its auditable component list."""
    step = _step_config(config, step_id)
    enabled = _enabled_tools(step)
    components: list[dict[str, Any]] = []

    # A tool-level timeout_seconds is an overall cap for that serial backend.
    for tool_id, tool in enabled.items():
        params = tool.get("params")
        params = params if isinstance(params, dict) else {}
        seconds = _positive_int(params.get("timeout_seconds"))
        if seconds:
            components.append({
                "component": f"tool:{tool_id}",
                "seconds": seconds,
                "calculation": "timeout_seconds",
            })

    step_params = step.get("params")
    step_params = step_params if isinstance(step_params, dict) else {}

    if step_id == "docking":
        # These existing wrappers have fixed serial caps. Moving them from
        # Step 0 must also move their budgets, so engine time remains available.
        for tool_id, seconds in {"reduce": 600, "pdb2pqr": 600, "meeko": 600,
                                 "ligand_microstates": 300}.items():
            if tool_id in enabled:
                components.append({"component": f"docking_prep:{tool_id}",
                                   "seconds": seconds, "calculation": "preparation wrapper hard cap"})

    if step_id == "pocket_analysis":
        native_timeout = _positive_int(
            step_params.get("native_step_timeout_seconds")
        )
        native_count = len(_POCKET_NATIVE_TOOLS & set(enabled))
        if native_timeout and native_count:
            components.append({
                "component": "native_interaction_profilers",
                "seconds": native_timeout * native_count,
                "calculation": (
                    f"native_step_timeout_seconds({native_timeout})"
                    f" * enabled_profilers({native_count})"
                ),
            })

    if step_id == "pose_evaluation":
        overall = _positive_int(step_params.get("overall_timeout_seconds"))
        if overall:
            components.append({
                "component": "core_pose_qc",
                "seconds": overall,
                "calculation": "overall_timeout_seconds",
            })

        topology = enabled.get("candidate_topology")
        if topology:
            params = topology.get("params") or {}
            if params.get("execution_mode") == "build":
                count = _positive_int(params.get("max_candidates"), 1)
                per_candidate = _positive_int(
                    params.get("timeout_per_candidate_seconds")
                )
                if per_candidate:
                    components.append({
                        "component": "candidate_topology",
                        "seconds": count * per_candidate,
                        "calculation": (
                            f"max_candidates({count})"
                            f" * timeout_per_candidate_seconds({per_candidate})"
                        ),
                    })

        mmgbsa = enabled.get("mmgbsa")
        if mmgbsa:
            params = mmgbsa.get("params") or {}
            if params.get("execution_mode") == "score":
                count = _positive_int(params.get("top_k"), 1)
                per_candidate = _positive_int(
                    params.get("timeout_per_candidate_seconds")
                )
                if per_candidate:
                    components.append({
                        "component": "mmgbsa",
                        "seconds": count * per_candidate,
                        "calculation": (
                            f"top_k({count})"
                            f" * timeout_per_candidate_seconds({per_candidate})"
                        ),
                    })

        ties = enabled.get("ties_rbfe")
        if ties:
            params = ties.get("params") or {}
            seconds = _positive_int(
                params.get("runtime_probe_timeout_seconds")
            )
            if seconds:
                components.append({
                    "component": "ties_rbfe",
                    "seconds": seconds,
                    "calculation": "runtime_probe_timeout_seconds",
                })

    serial_budget = sum(component["seconds"] for component in components)
    floor = int(DEFAULT_STEP_TIMEOUTS[step_id])
    calculated = 0
    if serial_budget:
        headroom = (
            serial_budget * _HEADROOM_PERCENT + 99
        ) // 100
        calculated = serial_budget + headroom + FIXED_OVERHEAD_SECONDS
    automatic = max(floor, calculated)
    return {
        "step": step_id,
        "clock": "active_execution_after_resource_admission",
        "includes_resource_wait": False,
        "floor_seconds": floor,
        "serial_budget_seconds": serial_budget,
        "headroom_ratio": HEADROOM_RATIO,
        "fixed_overhead_seconds": FIXED_OVERHEAD_SECONDS,
        "automatic_timeout_seconds": automatic,
        "components": components,
    }


def automatic_step_timeout(step_id: str, config: dict) -> int:
    """Resolve the safe automatic outer budget for one step."""
    budget = timeout_budget_breakdown(step_id, config)
    timeout = int(budget["automatic_timeout_seconds"])
    if timeout > MAX_STEP_TIMEOUT_SECONDS:
        raise TimeoutBudgetError(
            f"automatic timeout for {step_id} is {timeout}s, above the "
            f"{MAX_STEP_TIMEOUT_SECONDS}s (30 day) safety limit; split the "
            "library or reduce the candidate cohort"
        )
    return timeout


def resolve_step_timeout(step_id: str, config: dict) -> int:
    """Return an explicit validated hard cap or the automatic budget."""
    automatic = automatic_step_timeout(step_id, config)
    overrides = config.get("step_timeouts")
    if isinstance(overrides, dict) and step_id in overrides:
        explicit = overrides[step_id]
        if type(explicit) is not int or explicit < 1:
            raise TimeoutBudgetError(
                f"step_timeouts.{step_id} must be an integer >= 1"
            )
        if explicit > MAX_STEP_TIMEOUT_SECONDS:
            raise TimeoutBudgetError(
                f"step_timeouts.{step_id} must be <= "
                f"{MAX_STEP_TIMEOUT_SECONDS} seconds (30 days)"
            )
        if explicit < automatic:
            raise TimeoutBudgetError(
                f"step_timeouts.{step_id}={explicit} is below the required "
                f"{automatic}s budget for enabled serial tools"
            )
        return explicit
    return automatic

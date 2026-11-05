#!/usr/bin/env python3
"""Structured user-confirmation tool for the SBDD bio-research workflow.

The generic ``clarify`` tool intentionally remains flexible.  SBDD execution
gates need a narrower contract: fixed safe choices, a readable summary, a
typed decision, and a deterministic timeout policy.  This module reuses the
platform's existing clarify callback/TUI without changing the generic tool.
"""

from __future__ import annotations

import inspect
import json
import re
from typing import Any, Callable, NamedTuple, Optional

from tools.registry import registry, tool_error


SBDD_CONFIRM_TIMEOUT_SECONDS = 120
MAX_SUMMARY_CHARS = 6_000

_CLI_TIMEOUT_PREFIX = (
    "The user did not provide a response within the time limit."
)
_CLI_CANCEL_PREFIX = "The user cancelled."
_GATE_NAMES = frozenset({"plan", "step", "rewind", "final"})


class _GateSpec(NamedTuple):
    title: str
    question: str
    choices: tuple[str, str, str]


def _gate_spec(gate: str, rewind_target: Optional[str]) -> _GateSpec:
    if gate == "plan":
        return _GateSpec(
            "SBDD 实验方案确认",
            "是否按以上方案启动本轮 SBDD 实验？",
            ("确认并执行", "调整方案", "取消实验"),
        )
    if gate == "step":
        return _GateSpec(
            "SBDD 步骤结果确认",
            "是否接受本步骤结果并继续执行下一步？",
            ("继续下一步", "调整配置", "停止实验"),
        )
    if gate == "final":
        return _GateSpec(
            "SBDD 本轮结果确认",
            "是否接受本轮计算结果并进入总结或下一轮决策？",
            (
                "接受本轮结果并继续总结/迭代",
                "调整配置",
                "停止实验",
            ),
        )

    target = str(rewind_target or "").strip()
    return _GateSpec(
        "SBDD 重新执行确认",
        f"已有结果发生变化，是否确认从 {target} 重新执行？",
        (f"确认重新执行 {target}", "调整配置", "停止实验"),
    )


def _invalid_result(message: str, *, gate: Any = None) -> str:
    return tool_error(
        message,
        gate=gate,
        approved=False,
        decision=None,
        decision_source=None,
        decision_status="invalid_request",
        selected_choice=None,
        default_applied=False,
        timeout_seconds=SBDD_CONFIRM_TIMEOUT_SECONDS,
    )


def _callback_accepts_keyword(callback: Callable[..., Any], keyword: str) -> bool:
    """Return whether a callback can safely receive *keyword*."""
    try:
        parameters = inspect.signature(callback).parameters
    except (TypeError, ValueError):
        return False
    parameter = parameters.get(keyword)
    if parameter is not None and parameter.kind != inspect.Parameter.POSITIONAL_ONLY:
        return True
    return any(
        item.kind == inspect.Parameter.VAR_KEYWORD
        for item in parameters.values()
    )


def _invoke_callback(
    callback: Callable[..., Any],
    question: str,
    choices: list[str],
) -> Any:
    """Use the richer callback protocol when the platform supports it."""
    kwargs: dict[str, Any] = {}
    if _callback_accepts_keyword(callback, "timeout_seconds"):
        kwargs["timeout_seconds"] = SBDD_CONFIRM_TIMEOUT_SECONDS
    if _callback_accepts_keyword(callback, "return_status"):
        kwargs["return_status"] = True
    return callback(question, choices, **kwargs)


def _legacy_callback_status(raw_response: Any) -> tuple[str, str]:
    """Classify callbacks that still implement the legacy string protocol."""
    response = "" if raw_response is None else str(raw_response).strip()
    if response.startswith(_CLI_TIMEOUT_PREFIX):
        return "timeout", ""
    if re.fullmatch(r"\[clarify timed out after \d+(?:\.\d+)?s\]", response):
        return "timeout", ""
    if response.startswith(_CLI_CANCEL_PREFIX) or response == "[clarify cancelled]":
        return "cancelled", ""
    if response.startswith("[oneshot mode:"):
        return "unavailable", response
    if response == "[clarify prompt could not be delivered]":
        return "unavailable", response
    return ("selected" if response else "no_response"), response


def _callback_status(raw_response: Any) -> tuple[str, str]:
    if not isinstance(raw_response, dict):
        return _legacy_callback_status(raw_response)

    status = str(raw_response.get("status") or "").strip().lower()
    response = str(raw_response.get("response") or "").strip()
    if status in {
        "selected",
        "timeout",
        "cancelled",
        "unavailable",
        "delivery_failed",
        "callback_error",
        "no_response",
    }:
        return status, response
    return "callback_error", response


def _selected_choice(response: str, choices: list[str]) -> Optional[str]:
    if response in choices:
        return response
    return {
        str(index): choice
        for index, choice in enumerate(choices, start=1)
    }.get(response)


def sbdd_confirm_tool(
    gate: str,
    summary: str,
    rewind_target: Optional[str] = None,
    callback: Optional[Callable[..., Any]] = None,
) -> str:
    """Present one SBDD execution gate and return a structured decision."""
    if not isinstance(gate, str) or gate not in _GATE_NAMES:
        return _invalid_result(
            "gate must be one of: plan, step, rewind, final.",
            gate=gate,
        )
    if not isinstance(summary, str) or not summary.strip():
        return _invalid_result(
            "A non-empty natural-language summary is required.",
            gate=gate,
        )

    summary = summary.strip()
    if len(summary) > MAX_SUMMARY_CHARS:
        return _invalid_result(
            f"summary must not exceed {MAX_SUMMARY_CHARS} characters.",
            gate=gate,
        )

    if gate == "rewind":
        if not isinstance(rewind_target, str) or not rewind_target.strip():
            return _invalid_result(
                "rewind_target is required when gate='rewind'.",
                gate=gate,
            )
        rewind_target = rewind_target.strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,128}", rewind_target):
            return _invalid_result(
                "rewind_target must be a simple step identifier.",
                gate=gate,
            )
    elif rewind_target is not None and str(rewind_target).strip():
        return _invalid_result(
            "rewind_target is only valid when gate='rewind'.",
            gate=gate,
        )

    spec = _gate_spec(gate, rewind_target)
    choices = list(spec.choices)
    display_question = (
        f"{spec.title}\n\n"
        f"自然语言摘要：\n{summary}\n\n"
        f"{spec.question}\n"
        f"（{SBDD_CONFIRM_TIMEOUT_SECONDS} 秒内未响应时将默认选择："
        f"{choices[0]}）"
    )

    if callback is None:
        return tool_error(
            "SBDD confirmation is not available in this execution context.",
            gate=gate,
            summary=summary,
            question=spec.question,
            choices_offered=choices,
            user_response="",
            approved=False,
            decision=None,
            decision_source=None,
            decision_status="unavailable",
            selected_choice=None,
            default_applied=False,
            timeout_seconds=SBDD_CONFIRM_TIMEOUT_SECONDS,
        )

    try:
        raw_response = _invoke_callback(callback, display_question, choices)
    except Exception as exc:
        return tool_error(
            f"Failed to get SBDD confirmation: {exc}",
            gate=gate,
            summary=summary,
            question=spec.question,
            choices_offered=choices,
            user_response="",
            approved=False,
            decision=None,
            decision_source=None,
            decision_status="callback_error",
            selected_choice=None,
            default_applied=False,
            timeout_seconds=SBDD_CONFIRM_TIMEOUT_SECONDS,
        )

    callback_status, response = _callback_status(raw_response)
    default_applied = callback_status == "timeout"
    selected = choices[0] if default_applied else (
        _selected_choice(response, choices)
        if callback_status == "selected"
        else None
    )
    decision = None
    if selected is not None:
        decision = ("approve", "adjust", "stop")[choices.index(selected)]

    decision_source = None
    if decision is not None:
        decision_source = "timeout_default" if default_applied else "user"

    if default_applied:
        decision_status = "timeout_defaulted"
    elif callback_status == "selected" and selected is not None:
        decision_status = "selected"
    elif callback_status == "selected":
        decision_status = "no_match"
    else:
        decision_status = callback_status

    return json.dumps({
        "gate": gate,
        "summary": summary,
        "question": spec.question,
        "choices_offered": choices,
        "user_response": response,
        "approved": decision == "approve",
        "decision": decision,
        "decision_source": decision_source,
        "decision_status": decision_status,
        "selected_choice": selected,
        "default_applied": default_applied,
        "timeout_seconds": SBDD_CONFIRM_TIMEOUT_SECONDS,
    }, ensure_ascii=False)


def check_sbdd_confirm_requirements() -> bool:
    return True


SBDD_CONFIRM_SCHEMA = {
    "name": "sbdd_confirm",
    "description": (
        "Request a structured user decision at an SBDD workflow gate. The "
        "tool renders a required natural-language summary and fixed safe "
        "choices through the existing clarify UI. It waits 120 seconds; a "
        "real timeout selects the affirmative default, while cancellation, "
        "delivery failure, and callback errors do not approve execution. "
        "Use only for the bio_research SBDD plan, step, rewind, and final gates."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "gate": {
                "type": "string",
                "enum": ["plan", "step", "rewind", "final"],
                "description": "The SBDD confirmation boundary being shown.",
            },
            "summary": {
                "type": "string",
                "minLength": 1,
                "maxLength": MAX_SUMMARY_CHARS,
                "description": (
                    "A concise natural-language explanation for the user, not "
                    "raw config JSON. For plan gates summarize goals, tools, "
                    "steps, key cost parameters, and outputs. For step/final "
                    "gates summarize work completed, key results, artifacts, "
                    "warnings, and what would happen next."
                ),
            },
            "rewind_target": {
                "type": "string",
                "pattern": "^[A-Za-z0-9_.-]{1,128}$",
                "description": (
                    "Step identifier to rerun. Required only for gate='rewind'."
                ),
            },
        },
        "required": ["gate", "summary"],
        "additionalProperties": False,
    },
}


registry.register(
    name="sbdd_confirm",
    toolset="sbdd_confirmation",
    schema=SBDD_CONFIRM_SCHEMA,
    handler=lambda args, **kw: sbdd_confirm_tool(
        gate=args.get("gate", ""),
        summary=args.get("summary", ""),
        rewind_target=args.get("rewind_target"),
        callback=kw.get("callback"),
    ),
    check_fn=check_sbdd_confirm_requirements,
    emoji="🧬",
)

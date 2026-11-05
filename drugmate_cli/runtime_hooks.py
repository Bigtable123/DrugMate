"""Internal lifecycle hooks and tool policy dispatch.

Configured shell hooks use this dispatcher for per-session events and tool
approval gates. It does not scan directories, install packages, expose custom
commands, or load external Python extensions.
"""

from __future__ import annotations

import logging
import threading
from typing import Any, Callable, Dict, List, Optional, Set

from drugmate_cli.middleware import OBSERVER_SCHEMA_VERSION

logger = logging.getLogger(__name__)

VALID_HOOKS: Set[str] = {
    "pre_tool_call", "post_tool_call",
    "transform_terminal_output", "transform_tool_result", "transform_llm_output",
    "pre_llm_call", "post_llm_call", "pre_api_request", "post_api_request",
    "api_request_error", "on_session_start", "on_session_end",
    "on_session_finalize", "on_session_reset", "subagent_start", "subagent_stop",
    "pre_approval_request", "post_approval_response",
}


class HookManager:
    """Own callbacks registered by the built-in runtime."""

    def __init__(self) -> None:
        self._hooks: Dict[str, List[Callable]] = {}
        self._middleware: Dict[str, List[Callable]] = {}

    def invoke_hook(self, hook_name: str, **kwargs: Any) -> List[Any]:
        """Call all registered callbacks for *hook_name*.

        Each callback is wrapped in its own try/except so a misbehaving
        callback cannot break the core agent loop.

        Returns a list of non-``None`` return values from callbacks.

        For ``pre_llm_call``, callbacks may return a dict describing
        context to inject into the current turn's user message::

            {"context": "recalled text..."}
            "recalled text..."          # plain string, equivalent

        Context is ALWAYS injected into the user message, never the
        system prompt.  This preserves the prompt cache prefix — the
        system prompt stays identical across turns so cached tokens
        are reused.  All injected context is ephemeral — never
        persisted to session DB.
        """
        kwargs.setdefault("telemetry_schema_version", OBSERVER_SCHEMA_VERSION)
        callbacks = self._hooks.get(hook_name, [])
        results: List[Any] = []
        for cb in callbacks:
            try:
                ret = cb(**kwargs)
                if ret is not None:
                    results.append(ret)
            except Exception as exc:
                logger.warning(
                    "Hook '%s' callback %s raised: %s",
                    hook_name,
                    getattr(cb, "__name__", repr(cb)),
                    exc,
                )
        return results

    def has_hook(self, hook_name: str) -> bool:
        """Return True when at least one callback is registered for a hook."""
        return bool(self._hooks.get(hook_name))

    def has_middleware(self, kind: str) -> bool:
        """Return True when at least one callback is registered for middleware."""
        return bool(self._middleware.get(kind))

    def invoke_middleware(self, kind: str, **kwargs: Any) -> List[Any]:
        """Call registered middleware callbacks for *kind*.

        Each callback is isolated so one callback cannot break the base runtime
        path. Middleware that wants to change behavior must return the shape
        documented by the caller-specific contract.
        """
        callbacks = self._middleware.get(kind, [])
        results: List[Any] = []
        for cb in callbacks:
            try:
                ret = cb(**kwargs)
                if ret is not None:
                    results.append(ret)
            except Exception as exc:
                logger.warning(
                    "Middleware '%s' callback %s raised: %s",
                    kind,
                    getattr(cb, "__name__", repr(cb)),
                    exc,
                )
        return results


_manager = HookManager()


def get_hook_manager() -> HookManager:
    """Return the process-local internal callback dispatcher."""
    return _manager


def invoke_hook(hook_name: str, **kwargs: Any) -> List[Any]:
    """Invoke a lifecycle hook on the internal runtime dispatcher.

    Returns a list of non-``None`` return values from registered callbacks.
    """
    return get_hook_manager().invoke_hook(hook_name, **kwargs)


def invoke_middleware(kind: str, **kwargs: Any) -> List[Any]:
    """Invoke registered middleware callbacks.

    Returns a list of non-``None`` return values from middleware callbacks.
    """
    return get_hook_manager().invoke_middleware(kind, **kwargs)


def has_middleware(kind: str) -> bool:
    """Return True when middleware callbacks are registered for ``kind``."""
    manager = get_hook_manager()
    method = getattr(manager, "has_middleware", None)
    if callable(method):
        return bool(method(kind))
    return bool(getattr(manager, "_middleware", {}).get(kind))


def has_hook(hook_name: str) -> bool:
    """Return True when a hook has registered callbacks."""
    return get_hook_manager().has_hook(hook_name)


_thread_tool_whitelist = threading.local()


def set_thread_tool_whitelist(
    allowed: Optional[Set[str]],
    deny_msg_fmt: str = "Tool '{tool_name}' denied: not in this thread's tool whitelist",
) -> None:
    _thread_tool_whitelist.allowed = allowed
    _thread_tool_whitelist.fmt = deny_msg_fmt


def clear_thread_tool_whitelist() -> None:
    _thread_tool_whitelist.allowed = None


def get_pre_tool_call_block_message(
    tool_name: str,
    args: Optional[Dict[str, Any]],
    task_id: str = "",
    session_id: str = "",
    tool_call_id: str = "",
    turn_id: str = "",
    api_request_id: str = "",
    middleware_trace: Optional[List[Dict[str, Any]]] = None,
) -> Optional[str]:
    """Check ``pre_tool_call`` hooks for a blocking directive.

    Configured hooks that enforce policy (rate limiting, security
    restrictions, approval workflows) can return::

        {"action": "block", "message": "Reason the tool was blocked"}

    from their ``pre_tool_call`` callback.  The first valid block
    directive wins.  Invalid or irrelevant hook return values are
    silently ignored so existing observer-only hooks are unaffected.
    """
    allowed = getattr(_thread_tool_whitelist, "allowed", None)
    if allowed is not None and tool_name not in allowed:
        fmt = getattr(_thread_tool_whitelist, "fmt", "Tool '{tool_name}' denied")
        return fmt.format(tool_name=tool_name)

    hook_results = invoke_hook(
        "pre_tool_call",
        tool_name=tool_name,
        args=args if isinstance(args, dict) else {},
        task_id=task_id,
        session_id=session_id,
        tool_call_id=tool_call_id,
        turn_id=turn_id,
        api_request_id=api_request_id,
        middleware_trace=list(middleware_trace or []),
    )

    for result in hook_results:
        if not isinstance(result, dict):
            continue
        if result.get("action") != "block":
            continue
        message = result.get("message")
        if isinstance(message, str) and message:
            return message

    return None



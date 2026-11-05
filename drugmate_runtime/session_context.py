"""Task-local session identity for CLI tools and delegated research work."""

from contextvars import ContextVar, Token
from typing import Any

# Sentinel to distinguish "never set in this context" from "explicitly set to empty".
# When a contextvar holds _UNSET, we fall back to os.environ (CLI compatibility).
# When it holds "" (after clear_session_vars resets it), we return "" — no fallback.
_UNSET: Any = object()

# ---------------------------------------------------------------------------
# Per-task session variables
# ---------------------------------------------------------------------------

_SESSION_PLATFORM: ContextVar = ContextVar("DRUGMATE_SESSION_PLATFORM", default=_UNSET)
_SESSION_CHAT_ID: ContextVar = ContextVar("DRUGMATE_SESSION_CHAT_ID", default=_UNSET)
_SESSION_CHAT_NAME: ContextVar = ContextVar("DRUGMATE_SESSION_CHAT_NAME", default=_UNSET)
_SESSION_THREAD_ID: ContextVar = ContextVar("DRUGMATE_SESSION_THREAD_ID", default=_UNSET)
_SESSION_USER_ID: ContextVar = ContextVar("DRUGMATE_SESSION_USER_ID", default=_UNSET)
_SESSION_USER_NAME: ContextVar = ContextVar("DRUGMATE_SESSION_USER_NAME", default=_UNSET)
_SESSION_KEY: ContextVar = ContextVar("DRUGMATE_SESSION_KEY", default=_UNSET)
_SESSION_ID: ContextVar = ContextVar("DRUGMATE_SESSION_ID", default=_UNSET)
# Stable owner of SBDD assets created during this logical top-level Session.
# Unlike ``_SESSION_ID``, this value must not rotate when context compression
# creates a continuation Session or when ``delegate_task`` runs a child agent.
_SBDD_OWNER_SESSION_ID: ContextVar = ContextVar(
    "DRUGMATE_SBDD_OWNER_SESSION_ID", default=_UNSET
)
# ID of the message that triggered the current turn. Used as a reply anchor
# so background-process notifications stay inside the originating Telegram
# private-chat topic (those lanes route only with thread id + reply anchor).
_SESSION_MESSAGE_ID: ContextVar = ContextVar("DRUGMATE_SESSION_MESSAGE_ID", default=_UNSET)

_VAR_MAP = {
    "DRUGMATE_SESSION_PLATFORM": _SESSION_PLATFORM,
    "DRUGMATE_SESSION_CHAT_ID": _SESSION_CHAT_ID,
    "DRUGMATE_SESSION_CHAT_NAME": _SESSION_CHAT_NAME,
    "DRUGMATE_SESSION_THREAD_ID": _SESSION_THREAD_ID,
    "DRUGMATE_SESSION_USER_ID": _SESSION_USER_ID,
    "DRUGMATE_SESSION_USER_NAME": _SESSION_USER_NAME,
    "DRUGMATE_SESSION_KEY": _SESSION_KEY,
    "DRUGMATE_SESSION_ID": _SESSION_ID,
    "DRUGMATE_SBDD_OWNER_SESSION_ID": _SBDD_OWNER_SESSION_ID,
    "DRUGMATE_SESSION_MESSAGE_ID": _SESSION_MESSAGE_ID,
}


def set_current_session_id(session_id: str) -> Token:
    """Synchronize ``DRUGMATE_SESSION_ID`` across ContextVar and ``os.environ``.

    Long-lived single-process entrypoints like the CLI can rotate sessions via
    ``/new``, ``/resume``, ``/branch``, or compression splits without
    reconstructing the entire agent. Tools still consult
    ``get_session_env("DRUGMATE_SESSION_ID")`` with an ``os.environ`` fallback,
    so both storage paths must move together when the active session changes.
    """
    import os

    os.environ["DRUGMATE_SESSION_ID"] = session_id
    return _SESSION_ID.set(session_id)


def reset_current_session_id(token: Token) -> None:
    """Restore a temporary ``DRUGMATE_SESSION_ID`` binding.

    ``delegate_task`` uses this around child construction/execution so a child
    agent cannot leave the parent task pointing at the child's Session.
    """
    import os

    _SESSION_ID.reset(token)
    restored = _SESSION_ID.get()
    if restored is _UNSET:
        os.environ.pop("DRUGMATE_SESSION_ID", None)
    else:
        os.environ["DRUGMATE_SESSION_ID"] = restored


def set_sbdd_owner_session_id(session_id: str) -> Token:
    """Temporarily bind the stable SBDD owner for this execution context."""
    return _SBDD_OWNER_SESSION_ID.set(session_id)


def reset_sbdd_owner_session_id(token: Token) -> None:
    """Restore a prior stable SBDD owner binding."""
    _SBDD_OWNER_SESSION_ID.reset(token)


def set_session_vars(
    platform: str = "",
    chat_id: str = "",
    chat_name: str = "",
    thread_id: str = "",
    user_id: str = "",
    user_name: str = "",
    session_key: str = "",
    session_id: str = "",
    message_id: str = "",
    cwd: str = "",
    sbdd_owner_session_id: str | None = None,
) -> list:
    """Set all session context variables and return reset tokens.

    Call ``clear_session_vars(tokens)`` in a ``finally`` block when the handler
    exits. Note ``clear_session_vars`` resets every var to ``""`` (to suppress
    the ``os.environ`` fallback) rather than restoring prior values — these
    helpers are not nestable/stack-safe, and the returned tokens are accepted
    only for API compatibility.

    ``cwd`` pins the logical working directory for this context.
    """
    current_owner = _SBDD_OWNER_SESSION_ID.get()
    if sbdd_owner_session_id is None:
        # Preserve an inherited logical owner for nested contexts. A fresh
        # top-level session context starts with _UNSET/"" and is initialized
        # from its real Session id instead.
        owner_session_id = (
            current_owner
            if current_owner is not _UNSET and current_owner
            else session_id
        )
    else:
        owner_session_id = sbdd_owner_session_id

    tokens = [
        _SESSION_PLATFORM.set(platform),
        _SESSION_CHAT_ID.set(chat_id),
        _SESSION_CHAT_NAME.set(chat_name),
        _SESSION_THREAD_ID.set(thread_id),
        _SESSION_USER_ID.set(user_id),
        _SESSION_USER_NAME.set(user_name),
        _SESSION_KEY.set(session_key),
        _SESSION_ID.set(session_id),
        _SBDD_OWNER_SESSION_ID.set(owner_session_id),
        _SESSION_MESSAGE_ID.set(message_id),
    ]
    try:
        from agent.runtime_cwd import set_session_cwd

        set_session_cwd(cwd)
    except Exception:
        pass
    return tokens


def clear_session_vars(tokens: list) -> None:
    """Mark session context variables as explicitly cleared.

    Sets all variables to ``""`` so that ``get_session_env`` returns an empty
    string instead of falling back to (potentially stale) ``os.environ``
    values.  The *tokens* argument is accepted for API compatibility with
    callers that saved the return value of ``set_session_vars``, but the
    actual clearing uses ``var.set("")`` rather than ``var.reset(token)``
    to ensure the "explicitly cleared" state is distinguishable from
    "never set" (which holds the ``_UNSET`` sentinel).
    """
    for var in (
        _SESSION_PLATFORM,
        _SESSION_CHAT_ID,
        _SESSION_CHAT_NAME,
        _SESSION_THREAD_ID,
        _SESSION_USER_ID,
        _SESSION_USER_NAME,
        _SESSION_KEY,
        _SESSION_ID,
        _SBDD_OWNER_SESSION_ID,
        _SESSION_MESSAGE_ID,
    ):
        var.set("")
    try:
        from agent.runtime_cwd import clear_session_cwd

        clear_session_cwd()
    except Exception:
        pass


def get_session_env(name: str, default: str = "") -> str:
    """Read a session context variable by its legacy ``DRUGMATE_SESSION_*`` name.

    Drop-in replacement for ``os.getenv("DRUGMATE_SESSION_*", default)``.

    Resolution order:
    1. Context variable (set for task-local access).
       If the variable was explicitly set (even to ``""``) via
       ``set_session_vars`` or ``clear_session_vars``, that value is
       returned — **no fallback to os.environ**.
    2. ``os.environ`` (only when the context variable was never set in
       this context — i.e. CLI processes that
       don't use ``set_session_vars`` at all).
    3. *default*
    """
    import os

    var = _VAR_MAP.get(name)
    if var is not None:
        value = var.get()
        if value is not _UNSET:
            return value
    # Fall back to os.environ for CLI compatibility
    return os.getenv(name, default)

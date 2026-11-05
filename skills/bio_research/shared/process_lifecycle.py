"""Owned subprocess execution and cancellation for scientific workers.

Each command inherits a unique scope marker.  On Linux this also identifies
descendants that started another session or outlived their immediate parent;
cleanup never uses executable names or other users' resource processes.
"""
from __future__ import annotations

import contextlib
import ctypes
import os
import signal
import subprocess
import threading
import time
import uuid
from pathlib import Path


_SCOPE_ENV = "SBDD_PROCESS_SCOPES"
_subreaper_lock = threading.Lock()
_subreaper_users = 0
_subreaper_previous = 0
_signal_scope_active = False


class PipelineCancelled(BaseException):
    """Cancellation must escape backend handlers that catch Exception."""

    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"pipeline cancelled by {signal.Signals(signum).name}")


@contextlib.contextmanager
def _adopt_orphans():
    """Reap our exact owned descendants instead of leaving zombies to PID 1."""
    global _subreaper_users, _subreaper_previous
    libc = None
    if Path("/proc/self/stat").exists():
        try:
            libc = ctypes.CDLL(None, use_errno=True)
            with _subreaper_lock:
                if _subreaper_users == 0:
                    previous = ctypes.c_int()
                    if libc.prctl(37, ctypes.byref(previous), 0, 0, 0) != 0:
                        raise OSError("PR_GET_CHILD_SUBREAPER failed")
                    _subreaper_previous = previous.value
                    if libc.prctl(36, 1, 0, 0, 0) != 0:
                        raise OSError("PR_SET_CHILD_SUBREAPER failed")
                _subreaper_users += 1
        except (AttributeError, OSError):
            libc = None
    try:
        yield
    finally:
        if libc is not None:
            with _subreaper_lock:
                _subreaper_users -= 1
                if _subreaper_users == 0:
                    libc.prctl(36, _subreaper_previous, 0, 0, 0)


@contextlib.contextmanager
def cancellation_scope():
    """Turn terminal cancellation signals into a catchable unwind on main."""
    global _signal_scope_active
    if threading.current_thread() is not threading.main_thread() or _signal_scope_active:
        yield
        return
    previous = {}
    cancelling = False

    def cancel(signum, _frame):
        nonlocal cancelling
        if not cancelling:
            cancelling = True
            raise PipelineCancelled(signum)

    for signum in (signal.SIGINT, signal.SIGTERM, getattr(signal, "SIGHUP", None)):
        if signum is not None:
            previous[signum] = signal.signal(signum, cancel)
    _signal_scope_active = True
    try:
        yield
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)
        _signal_scope_active = False


def _process_info(pid: int) -> tuple[int, str, str] | None:
    try:
        # comm can contain spaces and parentheses; fields after its final ')' are fixed.
        fields = (Path("/proc") / str(pid) / "stat").read_text().rsplit(")", 1)[1].split()
        return int(fields[1]), fields[19], fields[0]  # ppid, starttime, state
    except (OSError, ValueError, IndexError):
        return None


def _owned_members(scope: str, root_pid: int, root_identity) -> dict[int, tuple[int, str, str]]:
    members = {}
    all_processes = {}
    try:
        entries = list(Path("/proc").iterdir())
    except OSError:
        return members
    for entry in entries:
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        info = _process_info(pid)
        if info is None:
            continue
        all_processes[pid] = info
        try:
            variables = (entry / "environ").read_bytes().split(b"\0")
            marker = next((value.split(b"=", 1)[1] for value in variables
                           if value.startswith((_SCOPE_ENV + "=").encode())), b"")
            if scope.encode() in marker.split(b":"):
                members[pid] = info
        except OSError:
            pass
    # Also catch children that deliberately replaced their environment.
    if (root_identity is not None and root_pid in all_processes
            and all_processes[root_pid][1] == root_identity[1]):
        members[root_pid] = all_processes[root_pid]
    while True:
        descendants = {pid: info for pid, info in all_processes.items()
                       if pid not in members and info[0] in members}
        if not descendants:
            return members
        members.update(descendants)


def _signal_owned(members: dict, signum: int, *, exclude: int | None = None):
    # Descendants first; direct parents get a chance to reap their own children.
    def depth(pid):
        seen = set()
        while pid in members and pid not in seen:
            seen.add(pid)
            pid = members[pid][0]
        return len(seen)

    for pid in sorted(members, key=depth, reverse=True):
        if pid == exclude:
            continue
        current = _process_info(pid)
        if current is None or current[1] != members[pid][1] or current[2] == "Z":
            continue
        try:
            os.kill(pid, signum)
        except (ProcessLookupError, PermissionError):
            pass


def _reap_owned(process, members):
    process.poll()  # Let Popen retain its own child's exit status.
    for pid, original in members.items():
        if pid == process.pid:
            continue
        current = _process_info(pid)
        if current is not None and current[1] == original[1] and current[2] == "Z":
            try:
                os.waitpid(pid, os.WNOHANG)
            except (ChildProcessError, OSError):
                pass


def _has_owned_orphans(scope: str) -> bool:
    """Inspect only our adopted children on normal exit, not all host PIDs."""
    try:
        tasks = list(Path("/proc/self/task").iterdir())
    except OSError:
        return False
    for task in tasks:
        try:
            children = (task / "children").read_text().split()
        except OSError:
            continue
        for pid in children:
            try:
                variables = (Path("/proc") / pid / "environ").read_bytes().split(b"\0")
                marker = next((value.split(b"=", 1)[1] for value in variables
                               if value.startswith((_SCOPE_ENV + "=").encode())), b"")
                if scope.encode() in marker.split(b":"):
                    return True
            except OSError:
                pass
    return False


def _stop_command(process, scope: str, root_identity, grace_seconds: float = 2.0):
    members = _owned_members(scope, process.pid, root_identity)
    if os.name == "posix":
        # The command always owns this new group. Also enumerate separately
        # above because scientific launchers may call setsid() themselves.
        owns_group = process.poll() is None
        for pid, identity in members.items():
            try:
                current = _process_info(pid)
                owns_group |= (current is not None and current[1] == identity[1]
                               and os.getpgid(pid) == process.pid)
            except ProcessLookupError:
                pass
        if owns_group:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    else:
        process.terminate()
    _signal_owned(members, signal.SIGTERM)
    deadline = time.monotonic() + grace_seconds
    while time.monotonic() < deadline:
        _reap_owned(process, members)
        if process.poll() is not None and all(
            (info := _process_info(pid)) is None or info[1] != original[1]
            for pid, original in members.items()
        ):
            break
        time.sleep(0.05)
    members.update(_owned_members(scope, process.pid, root_identity))
    _signal_owned(members, signal.SIGKILL)
    if process.poll() is None:
        process.kill()
    try:
        stdout, stderr = process.communicate(timeout=2)
    except subprocess.TimeoutExpired:
        # A daemon that scrubbed its environment and escaped ancestry may
        # retain a pipe; never let cancellation hang on draining that pipe.
        for stream in (process.stdout, process.stderr):
            if stream is not None:
                stream.close()
        process.wait(timeout=2)
        stdout, stderr = "", ""
    for _ in range(20):
        _reap_owned(process, members)
        if all(_process_info(pid) is None for pid in members):
            break
        time.sleep(0.01)
    return stdout, stderr


def run_command(argv, *, timeout, env=None, cwd=None, cancel_event=None):
    """Run argv with captured text output; clean the whole owned tree on abort."""
    if cancel_event is not None and cancel_event.is_set():
        raise PipelineCancelled(signal.SIGTERM)
    scope = uuid.uuid4().hex
    process_env = dict(os.environ if env is None else env)
    inherited = process_env.get(_SCOPE_ENV, "")
    process_env[_SCOPE_ENV] = f"{inherited}:{scope}" if inherited else scope
    with _adopt_orphans(), cancellation_scope():
        process = subprocess.Popen(
            argv, cwd=str(cwd) if cwd is not None else None, env=process_env,
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, start_new_session=os.name == "posix",
        )
        root_identity = None
        try:
            root_identity = _process_info(process.pid)
            if cancel_event is None:
                stdout, stderr = process.communicate(timeout=timeout)
            else:
                deadline = time.monotonic() + timeout if timeout is not None else None
                while True:
                    if cancel_event.is_set():
                        raise PipelineCancelled(signal.SIGTERM)
                    remaining = deadline - time.monotonic() if deadline is not None else None
                    try:
                        stdout, stderr = process.communicate(
                            timeout=min(0.1, max(0, remaining)) if remaining is not None else 0.1
                        )
                        break
                    except subprocess.TimeoutExpired:
                        if deadline is not None and time.monotonic() >= deadline:
                            raise
        except subprocess.TimeoutExpired as exc:
            stdout, stderr = _stop_command(process, scope, root_identity)
            raise subprocess.TimeoutExpired(argv, timeout, output=stdout or exc.output,
                                            stderr=stderr or exc.stderr) from None
        except BaseException as exc:
            stdout, stderr = _stop_command(process, scope, root_identity)
            exc.stdout, exc.stderr = stdout, stderr
            raise
        if _has_owned_orphans(scope):
            _stop_command(process, scope, root_identity)
        return subprocess.CompletedProcess(argv, process.returncode, stdout, stderr)

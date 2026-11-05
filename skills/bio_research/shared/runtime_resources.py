"""Host-wide cooperative SBDD admission with a separate queue budget."""
from __future__ import annotations

import contextlib
import fcntl
import hashlib
import json
import math
import os
import socket
import time
from pathlib import Path


DEFAULT_RESOURCE_SLOTS = 3


class ResourceUnavailable(RuntimeError):
    def __init__(self, message: str, *, waited_seconds: float = 0):
        super().__init__(message)
        self.result = {"success": False, "status": "resource_unavailable",
                       "retryable": True, "error": message,
                       "resource_wait_seconds": round(waited_seconds, 3)}


def validate_resource_options(options: dict) -> dict:
    """Validate the lease's internal options independently from scientific caps."""
    if not isinstance(options, dict):
        raise ValueError("runtime_resources must be an object")
    unknown = set(options) - {"slots", "wait_timeout_seconds", "poll_interval_seconds"}
    if unknown:
        raise ValueError(f"unknown runtime_resources keys: {', '.join(sorted(unknown))}")
    slots = options.get("slots", DEFAULT_RESOURCE_SLOTS)
    wait = options.get("wait_timeout_seconds", 0)
    poll = options.get("poll_interval_seconds", 1)
    if type(slots) is not int or not 1 <= slots <= 64:
        raise ValueError("runtime_resources.slots must be an integer from 1 to 64")
    if type(wait) is not int or not 0 <= wait <= 30 * 86400:
        raise ValueError("runtime_resources.wait_timeout_seconds must be an integer from 0 to 2592000")
    for name, value, minimum in (("wait_timeout_seconds", wait, 0),
                                 ("poll_interval_seconds", poll, 0.05)):
        if type(value) not in (int, float) or not math.isfinite(value) or value < minimum:
            raise ValueError(f"runtime_resources.{name} must be a finite number >= {minimum}")
    return {"slots": slots, "wait_timeout_seconds": wait, "poll_interval_seconds": poll}


def _open_slot(path: Path):
    flags = os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(path, flags | os.O_CREAT | os.O_EXCL, 0o666)
    except FileExistsError:
        fd = os.open(path, flags)
    else:
        # A shared data root can serve different Unix users; umask must not
        # accidentally create a lock other participants cannot acquire.
        os.fchmod(fd, 0o666)
    return os.fdopen(fd, "r+")


def _check_capacity(directory: Path, slots: int):
    # Persist the host policy so a caller cannot bypass an existing capacity
    # by choosing a larger local slot count. Operators may change the file
    # while no experiment is running; it is never removed during release.
    with _open_slot(directory / "capacity.lock") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        contents = handle.read().strip()
        if contents:
            try:
                capacity = json.loads(contents)["slots"]
            except (ValueError, KeyError, TypeError) as exc:
                raise ResourceUnavailable("invalid shared SBDD capacity.lock policy") from exc
            if capacity != slots:
                raise ResourceUnavailable(
                    f"shared SBDD resource capacity is {capacity}, requested {slots}; "
                    "all users must use the same runtime_resources.slots setting"
                )
        else:
            json.dump({"slots": slots}, handle)
            handle.flush()


@contextlib.contextmanager
def resource_lease(resource_dir: Path, *, wait_timeout_seconds=0, slots=DEFAULT_RESOURCE_SLOTS,
                   poll_interval_seconds=1, on_wait=None):
    """Acquire one host slot before starting any execution deadline.

    All participants must use the same data root and slot count. Lock files
    stay in place: unlinking a released lock would allow two inode owners.
    """
    options = validate_resource_options({"slots": slots,
                                        "wait_timeout_seconds": wait_timeout_seconds,
                                        "poll_interval_seconds": poll_interval_seconds})
    host = hashlib.sha256(socket.gethostname().encode()).hexdigest()[:16]
    directory = Path(resource_dir) / host
    try:
        if not directory.exists():
            directory.mkdir(parents=True, exist_ok=True, mode=0o777)
            directory.chmod(0o777)
        _check_capacity(directory, options["slots"])
    except OSError as exc:
        raise ResourceUnavailable(f"cannot access SBDD resource locks: {exc}") from exc
    started = time.monotonic()
    handle = None
    try:
        while handle is None:
            for slot in range(options["slots"]):
                try:
                    candidate = _open_slot(directory / f"slot-{slot}.lock")
                    try:
                        fcntl.flock(candidate.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    except BlockingIOError:
                        candidate.close()
                        continue
                    handle = candidate
                except OSError as exc:
                    raise ResourceUnavailable(f"cannot acquire SBDD resource lock: {exc}",
                                              waited_seconds=time.monotonic() - started) from exc
                break
            elapsed = time.monotonic() - started
            if handle is not None:
                break
            if elapsed >= options["wait_timeout_seconds"]:
                raise ResourceUnavailable(
                    "No SBDD execution slot is currently available; retry later or increase "
                    "runtime_resources.wait_timeout_seconds (queue time is separate from execution).",
                    waited_seconds=elapsed,
                )
            if on_wait is not None:
                on_wait({"status": "waiting_for_resources", "wait_seconds": round(elapsed, 3),
                         "wait_timeout_seconds": options["wait_timeout_seconds"]})
            time.sleep(min(options["poll_interval_seconds"],
                           options["wait_timeout_seconds"] - elapsed))
        record = {"slot": slot, "wait_seconds": round(time.monotonic() - started, 3),
                  "resource_directory": str(directory), "pid": os.getpid(),
                  "hostname": socket.gethostname()}
        handle.seek(0)
        handle.truncate()
        json.dump(record, handle)
        handle.flush()
        yield record
    finally:
        if handle is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()

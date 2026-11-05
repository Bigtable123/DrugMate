"""Bounded profiler execution and run-scoped I/O for pose evaluation."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import tempfile
from time import monotonic


def cached_file(cache: dict | None, path: Path, kind: object, loader):
    """Reuse parsed immutable inputs within one run, invalidating changed files."""
    if cache is None:
        return loader(path)
    stat = path.stat()
    key = (kind, str(path.resolve()))
    version = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    cached = cache.get(key)
    if cached is None or cached[0] != version:
        cached = (version, loader(path))
        cache[key] = cached
    return cached[1]


def profiler_worker_count(value=2) -> int:
    """Use a small CPU-affinity-bound budget; no nested pose/batch pools."""
    requested = int(value)
    if requested < 1:
        raise ValueError("profiler_workers must be at least 1")
    try:
        available = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        available = os.cpu_count() or 1
    return min(requested, 4, max(1, available))


def run_profiler_tasks(tasks: dict, workers: int, *, timings: dict | None = None) -> dict:
    """Return results/errors in submission order, irrespective of completion order."""
    def invoke(name, task):
        started = monotonic()
        try:
            return task(), None
        except Exception as exc:
            return None, exc
        finally:
            if timings is not None:
                timings[name] = round(monotonic() - started, 3)

    if workers == 1 or len(tasks) < 2:
        return {name: invoke(name, task) for name, task in tasks.items()}
    with ThreadPoolExecutor(max_workers=min(workers, len(tasks)), thread_name_prefix="pose-qc") as pool:
        pending = {name: pool.submit(invoke, name, task) for name, task in tasks.items()}
        return {name: future.result() for name, future in pending.items()}


def write_json(path: Path, payload, *, additional_paths=()) -> None:
    """Serialize compact JSON once and atomically publish each requested copy."""
    serialized = json.dumps(payload, ensure_ascii=False, default=str, separators=(",", ":")) + "\n"
    for target in (path, *additional_paths):
        target = Path(target)
        fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(serialized)
            os.replace(temporary, target)
        finally:
            Path(temporary).unlink(missing_ok=True)

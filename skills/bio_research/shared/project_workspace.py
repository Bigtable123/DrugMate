"""Request-bound round allocation and validation of initialized run directories."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from output_layout import write_json

INITIALIZATION = "project_initialization.json"


def _load(path: Path) -> dict:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected an object in {path}")
    return value


def _initialize(directory: Path, *, request_id: str | None, naming: dict,
                recover: bool = False) -> dict:
    record = directory / INITIALIZATION
    if recover and record.is_file():
        existing = _load(record)
        if (existing.get("request_id") != request_id or existing.get("naming") != naming
                or existing.get("project_dir") != str(directory)):
            raise ValueError("request binding does not match the initialization record")
        return existing
    if directory.exists():
        # A request binding is persisted before mkdir. Recover only an empty
        # allocation interrupted before its initialization record was committed.
        allowed = {directory / "inputs", directory / "inputs/normalization"}
        if not recover or any(p not in allowed or not p.is_dir() for p in directory.rglob("*")):
            raise ValueError(f"refusing to initialize existing directory: {directory}")
    directory.mkdir(parents=True, exist_ok=recover)
    inputs = directory / "inputs"
    normalization = inputs / "normalization"
    normalization.mkdir(parents=True, exist_ok=recover)
    result = {
        "success": True, "status": "project_initialized",
        "pipeline_invoked": False, "pipeline_complete": False,
        "project_dir": str(directory), "round_dir": str(directory),
        "project_root": str(directory.parent),
        "inputs_dir": str(inputs), "normalization_dir": str(normalization),
        "initialization_file": str(record), "request_id": request_id,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "drugmate_session_id": (os.environ.get("DRUGMATE_SBDD_OWNER_SESSION_ID")
                              or os.environ.get("DRUGMATE_SESSION_ID") or None),
        "naming": naming,
    }
    write_json(record, result)
    return result


def reserve(directory: Path, *, data_root: str | None, request_id: str | None,
            naming: dict, explicit_directory: bool = False) -> dict:
    """Reuse only a stable request identity; never infer identity from a PDB/slug."""
    if request_id is None:
        return _initialize(directory, request_id=None, naming=naming)
    if not isinstance(request_id, str) or not request_id.strip() or len(request_id) > 256:
        raise ValueError("request-id must be a nonempty string of at most 256 characters")
    if not data_root:
        raise ValueError("request-id requires an explicit project.data_root")
    request_id = request_id.strip()
    root = Path(data_root).expanduser().resolve()
    if not directory.is_relative_to(root):
        raise ValueError("request-bound round directory must be inside project.data_root")
    key = hashlib.sha256(f"{request_id}\n{naming['round_number']}".encode()).hexdigest()
    index = root / "_requests" / f"{key}.json"
    index.parent.mkdir(parents=True, exist_ok=True)
    with index.with_suffix(".lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        recovering = index.is_file()
        if recovering:
            binding = _load(index)
            if binding.get("request_id") != request_id or binding.get("naming") != naming:
                raise ValueError("request-id is already bound to different project naming; reuse its round_dir")
            bound_dir = Path(binding["round_dir"]).resolve()
            if not bound_dir.is_relative_to(root):
                raise ValueError("request binding points outside project.data_root")
            if explicit_directory and directory != bound_dir:
                raise ValueError(f"request-id is already bound to {bound_dir}")
            directory = bound_dir
        else:
            if directory.exists():
                raise ValueError(f"refusing to claim existing directory: {directory}")
            # Commit identity first so a crash cannot allocate a second UUID.
            write_json(index, {"request_id": request_id, "naming": naming,
                               "round_dir": str(directory)})
        return _initialize(directory, request_id=request_id, naming=naming, recover=recovering)


def validate_run_directory(directory: Path, round_number: int) -> None:
    """Reject a project root mistaken for its initialized round, before writes."""
    record = directory / INITIALIZATION
    if record.is_file():
        initialized = _load(record)
        if Path(initialized["project_dir"]).resolve() != directory:
            raise ValueError("initialized round directory does not match project_dir")
        declared_round = (initialized.get("naming") or {}).get("round_number")
        if declared_round is not None and declared_round != round_number:
            raise ValueError("project.round_number does not match the initialized round")
        return
    children = sorted(directory.glob(f"round*/{INITIALIZATION}")) if directory.is_dir() else []
    if not children:
        return  # Legacy layouts without an initialization record remain valid.
    checkpoint = directory / "pipeline_checkpoint.json"
    if checkpoint.is_file() and _load(checkpoint).get("project_dir") == str(directory):
        return  # Preserve already-running historical layouts; never migrate them.
    rounds = ", ".join(str(path.parent) for path in children)
    raise ValueError(f"project_dir is a project root, not an initialized round; use --project-dir {rounds}")

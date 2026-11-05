#!/usr/bin/env python3
"""Build an isolated reader-facing view for an existing SBDD round.

The legacy ``step*`` working tree is treated as immutable source material.  A
new publication is built in a sibling staging directory, verified, and renamed
into place only after the source tree's aggregate digest is unchanged.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any


SKILL_DIR = Path(__file__).resolve().parents[1]
COORDINATOR_RUN = SKILL_DIR / "run.py"
SHARED_DIR = SKILL_DIR.parent / "shared"
if str(SHARED_DIR) not in sys.path:
    sys.path.insert(0, str(SHARED_DIR))

import output_layout as ol  # noqa: E402
import partial_publication  # noqa: E402


def _load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ValueError(f"required round artifact is missing: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in required round artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"required round artifact must contain an object: {path}")
    return payload


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.resolve().relative_to(directory.resolve())
    except ValueError:
        return False
    return True


def _source_digest(round_dir: Path, excluded: Path) -> dict[str, Any]:
    """Hash path, size and bytes for every source file outside ``excluded``."""
    digest = hashlib.sha256()
    count = 0
    total_bytes = 0
    for path in sorted(round_dir.rglob("*")):
        if not path.is_file() or _is_within(path, excluded):
            continue
        relative = path.relative_to(round_dir).as_posix()
        file_digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1 << 20), b""):
                file_digest.update(block)
        size = path.stat().st_size
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(size).encode("ascii"))
        digest.update(b"\0")
        digest.update(file_digest.digest())
        count += 1
        total_bytes += size
    return {
        "files": count,
        "bytes": total_bytes,
        "sha256": digest.hexdigest(),
    }


def _load_coordinator_module():
    spec = importlib.util.spec_from_file_location(
        "sbdd_republish_coordinator", COORDINATOR_RUN,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load coordinator module: {COORDINATOR_RUN}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _verify_publication(
    publication_dir: Path,
    *,
    expected_session_id: str | None = None,
) -> dict[str, Any]:
    try:
        return ol.verify_publication(
            publication_dir,
            expected_session_id=expected_session_id,
            require_pipeline_final=True,
        )
    except ValueError as exc:
        raise RuntimeError(str(exc)) from exc


def republish(round_dir: Path, publication_dir: Path, *, partial: bool = False) -> dict[str, Any]:
    requested_publication = Path(publication_dir).expanduser()
    if requested_publication.is_symlink():
        raise ValueError(
            f"publication target must not be a symlink: {requested_publication}",
        )
    round_dir = Path(round_dir).expanduser().resolve()
    publication_dir = requested_publication.resolve()
    if not round_dir.is_dir():
        raise ValueError(f"round directory does not exist: {round_dir}")
    if publication_dir == round_dir:
        raise ValueError("publication directory must be separate from the working round")
    if publication_dir.exists():
        raise ValueError(
            f"publication directory already exists; refusing to overwrite: {publication_dir}"
        )

    lock_handle = None
    lock_path = round_dir / ".pipeline.lock"
    if not lock_path.is_file():
        raise ValueError(
            f"round has no pipeline lock; cannot prove it is idle: {lock_path}",
        )
    try:
        try:
            lock_handle = lock_path.open("r")
        except FileNotFoundError as exc:
            raise ValueError(
                f"round pipeline lock disappeared before acquisition: {lock_path}",
            ) from exc
        try:
            fcntl.flock(
                lock_handle.fileno(),
                fcntl.LOCK_EX | fcntl.LOCK_NB,
            )
        except BlockingIOError as exc:
            raise ValueError(
                f"round is busy; cannot republish while pipeline holds {lock_path}",
            ) from exc

        # Failed/limited rounds have no coordinator cohort. Build a diagnostic
        # snapshot from their durable manifest without changing source files.
        if partial:
            if publication_dir.exists():
                raise ValueError(f"publication directory already exists: {publication_dir}")
            snapshot = _load_json(round_dir / "manifest.json")
            # Older runner manifests omitted errors, even though result.json
            # and the checkpoint retained them. Enrich only failed steps;
            # never discover/promote pending results from the working tree.
            checkpoint_path = round_dir / "pipeline_checkpoint.json"
            checkpoint = _load_json(checkpoint_path) if checkpoint_path.is_file() else {}
            for step_id, info in (snapshot.get("steps") or {}).items():
                if info.get("status") not in {"failed", "failed_best_effort"}:
                    continue
                result_path = info.get("result_json")
                if result_path and Path(result_path).is_file():
                    try:
                        recorded = _load_json(Path(result_path))
                        if recorded.get("success") is False:
                            info["error"] = info.get("error") or recorded.get("error")
                    except ValueError as exc:
                        info["error"] = info.get("error") or str(exc)
                if checkpoint.get("current_step") == step_id:
                    info["error"] = info.get("error") or checkpoint.get("error")
            if checkpoint.get("error") and snapshot.get("status") != "complete":
                snapshot["error"] = snapshot.get("error") or checkpoint["error"]
            written = partial_publication.publish_snapshot(
                round_dir, snapshot, publication_dir=publication_dir,
            )
            return {
                "success": True, "source_round": str(round_dir),
                "publication_dir": str(publication_dir), **written,
                **ol.verify_publication(publication_dir, expected_session_id=round_dir.name),
            }

        source_before = _source_digest(round_dir, publication_dir)
        report = _load_json(round_dir / "step6" / "coordinator_report.json")
        eval_path = round_dir / "step4" / "evaluation_report.json"
        prop_path = round_dir / "step5" / "property_filter_report.json"
        docking_path = round_dir / "step3" / "docking_report.json"
        eval_report = _load_json(eval_path)
        prop_report = _load_json(prop_path)
        docking_report = _load_json(docking_path)
        source_loaded = _source_digest(round_dir, publication_dir)
        if source_loaded != source_before:
            raise RuntimeError("working tree changed while source reports were loaded")

        publication_dir.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=f".{publication_dir.name}.staging-",
            dir=str(publication_dir.parent),
        ))
        try:
            coordinator = _load_coordinator_module()
            coordinator._publish_output_layout(
                round_dir / "step6",
                report,
                eval_file=str(eval_path),
                prop_file=str(prop_path),
                docking_file=str(docking_path),
                eval_report=eval_report,
                prop_report=prop_report,
                docking_report=docking_report,
                publication_root=staging,
                publication_mode="pipeline_final",
            )
            verification = _verify_publication(
                staging,
                expected_session_id=round_dir.name,
            )
            source_after = _source_digest(round_dir, staging)
            if source_after != source_before:
                raise RuntimeError(
                    "working-tree digest changed while publication was built",
                )
            staging.rename(publication_dir)
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise
    finally:
        if lock_handle is not None:
            try:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
            finally:
                lock_handle.close()

    return {
        "success": True,
        "source_round": str(round_dir),
        "publication_dir": str(publication_dir),
        "source_integrity": source_before,
        **verification,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a non-destructive published view of an SBDD round.",
    )
    parser.add_argument("--round-dir", required=True, type=Path)
    parser.add_argument(
        "--partial", action="store_true",
        help="Publish existing files, step states and errors without requiring coordination/reporting.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Default: <round-dir>/published",
    )
    args = parser.parse_args()
    output_dir = args.output_dir or args.round_dir / "published"
    try:
        result = republish(args.round_dir, output_dir, partial=args.partial)
    except (OSError, RuntimeError, ValueError) as exc:
        print(json.dumps({"success": False, "error": str(exc)}, ensure_ascii=False))
        return 1
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

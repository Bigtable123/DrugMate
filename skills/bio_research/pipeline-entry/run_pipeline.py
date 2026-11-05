#!/usr/bin/env python3
"""Deterministic SBDD pipeline orchestrator (registry/artifact-driven).

Config is a three-layer model (step → tool → param). This script:
  1. verifies the explicit step/tool matrix and fills registered param defaults,
  2. validates it hard (VR1..VR10),
  3. writes a fully-expanded resolved_run_plan.json,
  4. runs each enabled step in registry order, wiring inter-step artifacts
     automatically (users never pass file paths),
  5. persists an atomic step cursor so interrupt mode can return control to the
     Agent for user confirmation without re-running completed work,
  6. keeps the agent out of the fragile parts: project dirs, per-step input JSON,
     timeouts, logs, manifests.

pipeline_adapters translates the resolved plan into each skill's input schema,
including its explicit tool switches and wired artifact paths.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile
import uuid
from contextlib import ExitStack
from datetime import datetime
from pathlib import Path
from typing import Any

BIO_ROOT = Path(__file__).resolve().parents[1]
SHARED_DIR = BIO_ROOT / "shared"
sys.path.insert(0, str(SHARED_DIR))

import pipeline_registry as registry
from deployment_config import deployment_path
import pipeline_config as pconfig
import pipeline_validation as pvalidation
import pipeline_adapters as adapters
import pipeline_timeouts
import tool_runtime
import output_layout as publication_layout
import partial_publication
import campaign
import project_workspace
from process_lifecycle import PipelineCancelled, cancellation_scope, run_command
from runtime_resources import DEFAULT_RESOURCE_SLOTS, ResourceUnavailable, resource_lease
from pipeline_artifacts import ArtifactStore, LEGACY_STEP_DIR


CHECKPOINT_SCHEMA_VERSION = 2
CHECKPOINT_FILENAME = "pipeline_checkpoint.json"
LOCK_FILENAME = ".pipeline.lock"
RUN_CONFIG_FILENAME = "pipeline_config.json"
SCIENTIFIC_OUTCOME_FILENAME = "scientific_outcome.json"
PUBLISHED_DIRNAME = "published"
DECISION_OUTCOMES = {"iterate", "success", "max_rounds", "no_viable", "invalid"}
# Session ids are opaque to the runner; the length cap mirrors the gateway's
# manifest validation so a corrupt provenance value never reaches a manifest.
_MAX_SESSION_ID_LENGTH = 256


class PipelineError(RuntimeError):
    def __init__(self, step: str, message: str, result: dict | None = None):
        super().__init__(message)
        self.step = step
        self.result = result or {}


def _quarantine_canonical_publication(
    project_dir: Path,
    *,
    reason: str,
) -> Path | None:
    """Atomically remove an obsolete reader view from its canonical path.

    A failed or rewound invocation must not leave ``published/`` looking like
    the output of that invocation.  The old derived view is cheap to rebuild,
    but retaining it outside the canonical path is useful for audit and manual
    recovery.  Renaming within the round directory is atomic on the supported
    filesystems.
    """
    publication = project_dir / PUBLISHED_DIRNAME
    if not publication.exists() and not publication.is_symlink():
        return None

    archive_root = project_dir / ".publication_archive"
    if archive_root.is_symlink():
        raise RuntimeError(
            f"publication archive root must not be a symlink: {archive_root}"
        )
    archive_root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    archive = archive_root / f"published-{stamp}-{uuid.uuid4().hex}"
    publication.rename(archive)
    audit_path = archive_root / f"{archive.name}.json"
    _atomic_write_json(audit_path, {
        "quarantined_at": datetime.now().isoformat(timespec="seconds"),
        "reason": reason,
        "canonical_path": str(publication),
        "archive_path": str(archive),
    })
    return archive


def _verified_publication(project_dir: Path) -> Path | None:
    """Return the current round's publication only when its manifest is exact."""
    published_dir = project_dir / PUBLISHED_DIRNAME
    try:
        publication_layout.verify_publication(
            published_dir,
            expected_session_id=project_dir.name,
            require_pipeline_final=True,
        )
    except (OSError, ValueError):
        return None
    return published_dir


def _result_may_expose_publication(result: dict, project_dir: Path) -> Path | None:
    """Partial evidence is viewable independently of execution success."""
    if result.get("publication_mode") == "pipeline_partial":
        try:
            verification = publication_layout.verify_publication(
                project_dir / PUBLISHED_DIRNAME, expected_session_id=project_dir.name,
            )
        except (OSError, ValueError):
            return None
        if verification["publication_mode"] == "pipeline_partial":
            return project_dir / PUBLISHED_DIRNAME
        return None
    if (
        result.get("success") is not True
        or result.get("pipeline_complete") is not True
    ):
        return None
    if not {"coordination", "reporting"}.issubset(
        set(result.get("steps_completed") or ())
    ):
        return None
    return _verified_publication(project_dir)


def _load_json_object(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read publication source {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"publication source must contain an object: {path}")
    return payload


def _canonical_artifact_path(store: ArtifactStore, name: str) -> Path:
    """Resolve and revalidate the active canonical path for final publication."""
    value = store.get(name)
    if not isinstance(value, (str, os.PathLike)) or not str(value).strip():
        raise ValueError(
            f"required canonical publication artifact is unavailable: {name}"
        )
    path = Path(os.fspath(value)).expanduser().resolve()
    if not path.is_file():
        raise ValueError(
            f"required canonical publication artifact is missing: {name} ({path})"
        )
    return path


def _load_coordinator_publisher():
    """Load the existing coordinator publication adapter without importing it globally."""
    path = BIO_ROOT / "reflective-coordinator" / "run.py"
    module_name = "_sbdd_pipeline_final_publication"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load publication adapter: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _finalize_pipeline_publication(
    project_dir: Path,
    plan: dict,
    results: dict[str, dict],
    store: ArtifactStore,
) -> Path | None:
    """Seal the reader view only after both final pipeline results are durable."""
    required_steps = {"coordination", "reporting"}
    scheduled_steps = set(plan.get("execution_order") or ())
    scheduled_final_steps = required_steps & scheduled_steps
    if not scheduled_final_steps:
        return None
    if scheduled_final_steps != required_steps:
        missing = sorted(required_steps - scheduled_final_steps)
        raise RuntimeError(
            "final publication boundary is incomplete; missing scheduled step(s): "
            + ", ".join(missing)
        )
    incomplete = [
        step_id
        for step_id in sorted(required_steps)
        if (results.get(step_id) or {}).get("success") is not True
    ]
    if incomplete:
        raise RuntimeError(
            "final publication requires successful final steps: "
            + ", ".join(incomplete)
        )

    for step_id in ("coordination", "reporting"):
        inputs_match, reasons = store.verify_result_inputs(
            step_id, results[step_id],
        )
        if not inputs_match:
            raise RuntimeError(
                f"{step_id} input lineage changed before publication: "
                + "; ".join(reasons)
            )

    current = _verified_publication(project_dir)
    if current is not None:
        return current

    coordinator_path = _canonical_artifact_path(store, "coordinator_report")
    eval_path = _canonical_artifact_path(store, "evaluation_report")
    property_path = _canonical_artifact_path(store, "property_report")
    docking_path = _canonical_artifact_path(store, "docking_report")
    # The reporting step is part of the durable finalization boundary.  Its
    # bundle manifest is independently validated by the artifact extractor;
    # require the concrete primary files too so a successful-but-empty result
    # cannot unlock publication.
    for name in (
        "report_html",
        "report_data",
        "report_assets_manifest",
        "report_artifact_manifest",
    ):
        _canonical_artifact_path(store, name)
    coordinator = _load_coordinator_publisher()
    coordinator._publish_output_layout(
        project_dir / LEGACY_STEP_DIR["coordination"],
        _load_json_object(coordinator_path),
        eval_file=str(eval_path),
        prop_file=str(property_path),
        docking_file=str(docking_path),
        eval_report=_load_json_object(eval_path),
        prop_report=_load_json_object(property_path),
        docking_report=_load_json_object(docking_path),
        publication_mode="pipeline_final",
    )
    publication = _verified_publication(project_dir)
    if publication is None:
        raise RuntimeError("final publication did not pass strict verification")
    return publication


def run(config: dict) -> dict:
    with cancellation_scope():
        return _run(config)


def _run(config: dict) -> dict:
    if not isinstance(config, dict):
        return {
            "success": False,
            "status": "config_invalid",
            "pipeline_complete": False,
            "requires_confirmation": False,
            "error": "pipeline config must be an object",
            "config_errors": ["pipeline config must be an object"],
        }
    try:
        project_dir = _project_dir(config)
        project_workspace.validate_run_directory(
            project_dir, _run_setting(config, "round_number", 1),
        )
    except (ValueError, OSError, KeyError) as exc:
        return {
            "success": False, "status": "config_invalid", "pipeline_complete": False,
            "requires_confirmation": False, "error": str(exc), "config_errors": [str(exc)],
        }
    project_dir.mkdir(parents=True, exist_ok=True)
    lock_path = project_dir / LOCK_FILENAME
    lock_handle = lock_path.open("a+")
    resources = ExitStack()
    try:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return {
                "success": False,
                "status": "project_busy",
                "pipeline_complete": False,
                "requires_confirmation": False,
                "project_dir": str(project_dir),
                "error": (
                    "another pipeline invocation is already using this "
                    "project_dir"
                ),
                "lock_file": str(lock_path),
            }
        lock_handle.seek(0)
        lock_handle.truncate()
        json.dump(
            {
                "pid": os.getpid(),
                "acquired_at": datetime.now().isoformat(timespec="seconds"),
            },
            lock_handle,
        )
        lock_handle.flush()
        os.fsync(lock_handle.fileno())
        campaign_errors = campaign.validate_campaign_config(config, runtime=True)
        if campaign_errors:
            return {
                "success": False, "status": "config_invalid",
                "pipeline_complete": False, "requires_confirmation": False,
                "project_dir": str(project_dir),
                "error": "campaign plan validation failed",
                "config_errors": campaign_errors,
            }
        try:
            campaign_run = campaign.reserve_campaign_run(config, project_dir)
        except (ValueError, OSError) as exc:
            return {
                "success": False, "status": "config_invalid",
                "pipeline_complete": False, "requires_confirmation": False,
                "project_dir": str(project_dir),
                "error": f"campaign reservation failed: {exc}",
                "config_errors": [str(exc)],
            }
        config_path = project_dir / RUN_CONFIG_FILENAME
        manifest_path = project_dir / "manifest.json"
        previous_manifest = None
        publication_archives: list[str] = []
        diagnostic_published = False
        try:
            previous_manifest = manifest_path.read_bytes() if manifest_path.is_file() else None
            _atomic_write_json(config_path, config)
            result = _run_locked(
                config,
                project_dir,
                publication_archives=publication_archives,
                resource_stack=resources,
            )
        except BaseException as exc:
            # Preserve the last validated prefix from this invocation only.
            # An exception during config/checkpoint loading must not publish
            # successful results from an older invocation as current work.
            try:
                fresh = manifest_path.is_file() and manifest_path.read_bytes() != previous_manifest
            except OSError:
                fresh = False
            interrupted = isinstance(exc, (PipelineCancelled, KeyboardInterrupt, SystemExit))
            status = "interrupted" if interrupted else "failed"
            if isinstance(exc, ResourceUnavailable):
                status = "resource_unavailable"
            _publish_diagnostic_view(
                project_dir, status=status,
                error=f"{type(exc).__name__}: {exc}", preserve_steps=fresh,
            )
            diagnostic_published = True
            if not isinstance(exc, (PipelineCancelled, KeyboardInterrupt, ResourceUnavailable)):
                try:
                    campaign.complete_campaign_run(config, project_dir, {
                        "success": False, "status": status,
                        "pipeline_complete": False, "error": str(exc),
                    })
                except (ValueError, OSError) as state_exc:
                    if hasattr(exc, "add_note"):
                        exc.add_note(f"Campaign state also failed to persist: {state_exc}")
                raise
            result = {
                "success": False, "status": status, "pipeline_complete": False,
                "requires_confirmation": False, "project_dir": str(project_dir),
                "manifest_json": str(manifest_path), "error": str(exc),
                "retryable": isinstance(exc, ResourceUnavailable),
            }
            if isinstance(exc, ResourceUnavailable):
                result.update(exc.result)
        if result.get("success") is not True and not diagnostic_published:
            _publish_diagnostic_view(
                project_dir, status=result.get("status") or "failed",
                error=result.get("error") or "pipeline invocation failed",
                preserve_steps=result.get("status") not in {"config_invalid", "checkpoint_invalid"},
            )
        if publication_archives:
            result["quarantined_publications"] = list(publication_archives)
        result.setdefault("pipeline_config", str(config_path))
        result.setdefault("working_dir", str(project_dir))
        manifest = _read_publication_manifest(project_dir)
        for key in ("publication_mode", "publication_error"):
            if key in manifest:
                result[key] = manifest[key]
        published_dir = _result_may_expose_publication(result, project_dir)
        if published_dir is not None:
            result.setdefault("published_dir", str(published_dir))
            result.setdefault("user_output_dir", str(published_dir))
            report = published_dir / "07_final_report" / "final_report.md"
            if report.is_file():
                result["execution_report"] = str(report)
        if campaign_run is not None:
            result["campaign"] = {
                key: campaign_run[key] for key in ("run_id", "stage", "attempt")
            }
            result["campaign"]["state_file"] = config["project"]["campaign_state_file"]
            try:
                campaign.complete_campaign_run(config, project_dir, result)
            except (ValueError, OSError) as exc:
                # Completed scientific evidence remains on disk; a persistence
                # error must not authorize the next stage or trigger recompute.
                result["campaign_state_error"] = str(exc)
                result["success"] = False
                result["status"] = "campaign_state_failed"
                result["requires_confirmation"] = False
        return result
    finally:
        resources.close()
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        finally:
            lock_handle.close()


def _context_artifact_values(context: dict) -> dict[str, Any]:
    """Map explicit context files onto the canonical artifact vocabulary."""
    return {
        artifact: context.get(context_key)
        for artifact, context_key in registry.CONTEXT_ARTIFACTS.items()
        if context.get(context_key)
    }


def _inject_explicit_artifacts(store: ArtifactStore, merged: dict) -> None:
    """Seed context and input artifacts after rejecting ambiguous ownership."""
    context_values = _context_artifact_values(merged.get("context") or {})
    input_values = merged.get("input_artifacts") or {}
    conflicts = sorted(set(context_values) & {
        name for name, value in input_values.items() if value
    })
    if conflicts:
        raise ValueError(
            "canonical artifacts have both context and input_artifacts sources: "
            + ", ".join(conflicts)
        )
    store.inject(
        context_values,
        producer="context",
        source="config.context",
    )
    store.inject(input_values, source="config.input_artifacts")


def _source_is_within(source: Any, directory: Path | None) -> bool:
    """Return whether a provenance source is a file below ``directory``."""
    if not source or directory is None:
        return False
    try:
        Path(source).expanduser().resolve().relative_to(directory.resolve())
    except (OSError, RuntimeError, TypeError, ValueError):
        return False
    return True


def _scheduled_invalidations(
    step_id: str,
    step_plan: dict,
    store: ArtifactStore,
    *,
    upstream_project_dir: Path | None,
) -> list[str]:
    """Select outputs owned by this execution rather than external seeds."""
    declared = list(step_plan["produces"])
    active = set(step_plan.get("active_produces", declared))
    invalidations: list[str] = []
    for artifact in declared:
        if artifact in active:
            invalidations.append(artifact)
            continue
        provenance = store.provenance(artifact)
        if provenance.get("producer") != step_id:
            continue
        if _source_is_within(provenance.get("source"), upstream_project_dir):
            continue
        invalidations.append(artifact)
    return invalidations


def _run_locked(
    config: dict,
    project_dir: Path,
    *,
    publication_archives: list[str] | None = None,
    resource_stack: ExitStack | None = None,
) -> dict:
    logs_dir = project_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    project = config.get("project") or {}
    upstream_raw = project.get("upstream_project_dir") or config.get("upstream_project_dir")
    upstream_project_dir = (
        Path(upstream_raw).expanduser().resolve() if upstream_raw else None
    )

    # --- resolve + validate config ---
    try:
        merged = pconfig.merge_user_config(config)
        preload = _preload_artifacts(project_dir, upstream_project_dir, merged)
        _inject_explicit_artifacts(preload, merged)
        _write_artifact_lineage(logs_dir, preload)
        pvalidation.validate(merged, available_artifacts=preload.names())
    except (pvalidation.ConfigValidationError, ValueError) as exc:
        _append_error(logs_dir, f"config validation failed: {exc}")
        errors = exc.errors if hasattr(exc, "errors") else [str(exc)]
        return {
            "success": False,
            "status": "config_invalid",
            "pipeline_complete": False,
            "requires_confirmation": False,
            "project_dir": str(project_dir),
            "error": str(exc),
            "config_errors": errors,
            "artifact_lineage": str(logs_dir / "artifact_lineage.json"),
        }
    plan = pconfig.resolve(merged)
    producer_config_fingerprints = _producer_config_fingerprints(merged)
    checkpoint_path = project_dir / CHECKPOINT_FILENAME
    try:
        resume_state = _restore_checkpoint(
            checkpoint_path=checkpoint_path,
            project_dir=project_dir,
            plan=plan,
            merged=merged,
            store=preload,
            producer_config_fingerprints=producer_config_fingerprints,
        )
    except ValueError as exc:
        message = f"pipeline checkpoint is invalid: {exc}"
        _append_error(logs_dir, message)
        return {
            "success": False,
            "status": "checkpoint_invalid",
            "pipeline_complete": False,
            "requires_confirmation": False,
            "project_dir": str(project_dir),
            "error": message,
            "checkpoint_json": str(checkpoint_path),
            "resolved_run_plan": (
                str(project_dir / "resolved_run_plan.json")
                if (project_dir / "resolved_run_plan.json").exists()
                else None
            ),
            "artifact_lineage": str(logs_dir / "artifact_lineage.json"),
        }

    # Record the explicit-or-derived provenance in the plan before it is serialized, so
    # ``resolved_run_plan.json`` shows which screening-mode default the run
    # picked instead of leaving it to be re-derived from the step matrix.
    plan.setdefault("context", {})["candidate_source"] = _candidate_source(
        plan, preload
    )
    (project_dir / "resolved_run_plan.json").write_text(
        json.dumps(plan, indent=2, ensure_ascii=False, default=str)
    )

    results: dict[str, dict] = resume_state["results"]
    processed_steps: list[str] = resume_state["processed_steps"]
    processed_result_fingerprints: dict[str, str] = resume_state[
        "processed_result_fingerprints"
    ]
    resume_rewinds: list[str] = resume_state["rewind_reasons"]
    pending_steps = plan["execution_order"][len(processed_steps):]
    if resume_rewinds:
        archived = _quarantine_canonical_publication(
            project_dir,
            reason="pipeline has rewound work",
        )
        if archived is not None and publication_archives is not None:
            publication_archives.append(str(archived))
    execution_mode = str(
        (plan.get("project") or {}).get("execution_mode", "interrupt")
    )
    _write_artifact_lineage(logs_dir, preload)

    # A changed config/result/artifact can invalidate an already-confirmed
    # prefix. In interrupt mode, persist the safe prefix and return before
    # executing the rewind target: approval to run the old "next step" must
    # never be stretched into approval to redo an earlier, possibly expensive
    # step.
    if execution_mode == "interrupt" and resume_rewinds and pending_steps:
        rewind_target = pending_steps[0]
        _write_checkpoint(
            checkpoint_path,
            state="awaiting_step_confirmation",
            project_dir=project_dir,
            plan=plan,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            results=results,
            producer_config_fingerprints=producer_config_fingerprints,
            current_step=None,
            next_step=rewind_target,
        )
        _write_manifest(
            project_dir,
            f"awaiting_confirmation_for_rewind_to_{rewind_target}",
            plan,
            results,
        )
        return _rewind_pause_result(
            project_dir=project_dir,
            plan=plan,
            results=results,
            rewind_target=rewind_target,
            checkpoint_path=checkpoint_path,
            resume_rewinds=resume_rewinds,
        )

    # Runtime preflight is repeated for the work that remains, but already
    # checkpointed steps are not re-probed or re-executed merely because the
    # Agent paused to ask the user for confirmation.
    _write_checkpoint(
        checkpoint_path,
        state="preflight",
        project_dir=project_dir,
        plan=plan,
        processed_steps=processed_steps,
        processed_result_fingerprints=processed_result_fingerprints,
        results=results,
        producer_config_fingerprints=producer_config_fingerprints,
        current_step=None,
        next_step=pending_steps[0] if pending_steps else None,
    )
    # Publish the Session-owned project before runtime probing begins. Some
    # probes are slow enough that waiting until the first step starts leaves a
    # newly opened Session with no bindable project in the workspace UI.
    _write_manifest(project_dir, "preflight", plan, results)
    if pending_steps and resource_stack is not None:
        settings = merged.get("runtime_resources") or {}
        data_root = _run_setting(merged, "data_root", None) or os.environ.get("SBDD_DATA_ROOT")
        resource_dir = Path(settings.get("lock_directory") or (
            Path(data_root).expanduser() / "_resources" if data_root
            else deployment_path("${SBDD_DATA_ROOT}/_resources")
        )).expanduser().resolve()
        wait_published = False

        def waiting(details: dict) -> None:
            nonlocal wait_published
            _atomic_write_json(project_dir / "resource_status.json", {
                **details, "status": "waiting_resources",
            })
            if not wait_published:
                _write_manifest(project_dir, "waiting_resources", plan, results)
                wait_published = True

        lease = resource_stack.enter_context(resource_lease(
            resource_dir,
            wait_timeout_seconds=settings.get("wait_timeout_seconds", 0),
            slots=settings.get("slots", DEFAULT_RESOURCE_SLOTS), on_wait=waiting,
        ))
        _atomic_write_json(project_dir / "resource_status.json", {"status": "acquired", **lease})
        if wait_published:
            _write_manifest(project_dir, "preflight", plan, results)
    runtime_plan = copy.deepcopy(plan)
    runtime_plan["execution_order"] = list(pending_steps)
    fresh_runtime_report = tool_runtime.probe_plan(runtime_plan)
    runtime_report_path = project_dir / "tool_runtime_report.json"
    runtime_report = _merge_runtime_report(
        plan=plan,
        fresh_report=fresh_runtime_report,
        existing_path=runtime_report_path,
        processed_steps=processed_steps,
    )
    _atomic_write_json(runtime_report_path, runtime_report)
    if not runtime_report["success"]:
        missing = ", ".join(runtime_report["missing"])
        message = f"enabled tool runtime preflight failed: {missing}"
        _append_error(logs_dir, message)
        _write_checkpoint(
            checkpoint_path,
            state="preflight",
            project_dir=project_dir,
            plan=plan,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            results=results,
            producer_config_fingerprints=producer_config_fingerprints,
            current_step=None,
            next_step=pending_steps[0] if pending_steps else None,
            error=message,
        )
        _write_manifest(project_dir, "runtime_preflight_failed", plan, results)
        return {
            "success": False,
            "status": "runtime_preflight_failed",
            "pipeline_complete": False,
            "requires_confirmation": False,
            "project_dir": str(project_dir),
            "error": message,
            "runtime_errors": runtime_report["errors"],
            "processed_steps": list(processed_steps),
            "steps_completed": _successful_steps(plan, results),
            "tool_runtime_report": str(runtime_report_path),
            "checkpoint_json": str(checkpoint_path),
            "resolved_run_plan": str(project_dir / "resolved_run_plan.json"),
            "artifact_lineage": str(logs_dir / "artifact_lineage.json"),
        }

    for step_id in plan["execution_order"]:
        (project_dir / LEGACY_STEP_DIR[step_id]).mkdir(exist_ok=True)
    _write_env_paths(project_dir, logs_dir, merged)

    store = preload
    context = dict(plan["context"])
    # Level 2 of the gate contract needs to know where the molecules came from
    #.  Every step that folds metric marks accepts an explicit
    # `screening_mode`; when the Agent leaves it unset, this derived source is
    # what selects the default -- V3 for freshly generated molecules, V2 for a
    # sampled library, V1 (marking only) for a user-supplied full set.
    context["candidate_source"] = _candidate_source(plan, store)
    docking_tools = ((plan.get("steps") or {}).get("docking") or {}).get("tools") or {}
    docking_prep_params = (docking_tools.get("docking_prep") or {}).get("params") or {}
    # Retain the adapter metadata key while taking its value from the planned
    # final protonation context, rather than the geometry-stage diagnostic.
    structure_prep_ph = docking_prep_params.get("ph")
    round_number = int(_run_setting(config, "round_number", 1))
    resume_existing = bool(_run_setting(config, "resume_existing", True))
    _write_manifest(
        project_dir,
        "resumed" if processed_steps else "initialized",
        plan,
        results,
    )

    if not pending_steps:
        publication_failure = _finalize_publication_or_failure(
            project_dir=project_dir,
            plan=plan,
            results=results,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            producer_config_fingerprints=producer_config_fingerprints,
            checkpoint_path=checkpoint_path,
            store=store,
        )
        if publication_failure is not None:
            return publication_failure
        _write_checkpoint(
            checkpoint_path,
            state="complete",
            project_dir=project_dir,
            plan=plan,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            results=results,
            producer_config_fingerprints=producer_config_fingerprints,
            current_step=None,
            next_step=None,
        )
        _write_manifest(project_dir, "complete", plan, results)
        return _final_result(
            project_dir=project_dir,
            plan=plan,
            results=results,
            runtime_report_path=runtime_report_path,
            checkpoint_path=checkpoint_path,
            execution_mode=execution_mode,
            resume_rewinds=resume_rewinds,
        )

    try:
        pending_index = 0
        while pending_index < len(pending_steps):
            step_id = pending_steps[pending_index]
            step_plan = plan["steps"][step_id]
            _write_checkpoint(
                checkpoint_path,
                state="running",
                project_dir=project_dir,
                plan=plan,
                processed_steps=processed_steps,
                processed_result_fingerprints=processed_result_fingerprints,
                results=results,
                producer_config_fingerprints=producer_config_fingerprints,
                current_step=step_id,
                next_step=step_id,
            )
            # A step that is about to run owns a fresh generation of all of its
            # outputs. Clear any preloaded result for this producer and every
            # downstream derivative before wiring inputs. A fingerprint-matched
            # cached result may immediately republish them below.
            # ``produces`` is the registry capability envelope. Conditional
            # tool outputs are only replaced when their producer is enabled in
            # this resolved plan; injected artifacts otherwise remain valid.
            for artifact in _scheduled_invalidations(
                step_id,
                step_plan,
                store,
                upstream_project_dir=upstream_project_dir,
            ):
                store.invalidate(
                    artifact,
                    reason=f"{step_id} scheduled for a fresh output generation",
                )
            _write_artifact_lineage(logs_dir, store)
            _write_manifest(project_dir, f"running_{step_id}", plan, results)
            wired_names = list(step_plan["consumes"])
            expected_optional_inputs = {
                name
                for name in step_plan.get("optional_consumes", [])
                if _use_optional_artifact(step_id, name, context, step_plan.get("params"))
            }
            wired_names.extend(sorted(expected_optional_inputs))
            wired = {
                name: store.get(name)
                for name in wired_names
                if store.has(name)
            }
            artifact_inputs = {
                name: {
                    key: value
                    for key, value in store.provenance(name).items()
                    if key in {"producer", "source", "fingerprint", "kind", "path"}
                }
                for name in wired
            }
            runtime = {
                "project_dir": str(project_dir),
                "upstream_project_dir": str(upstream_project_dir) if upstream_project_dir else None,
                "round_number": round_number,
                "inferred_target_name": _infer_target_name(results.get("structure_prep", {})),
                "structure_prep_ph": structure_prep_ph,
            }
            step_input = adapters.to_legacy_input(step_id, wired, step_plan, context, runtime)
            timeout = _step_timeout(step_id, merged)
            producer_config_fingerprint = producer_config_fingerprints[step_id]

            try:
                result = _run_step(
                    step_id=step_id, skill=step_plan["skill"], step_input=step_input,
                    project_dir=project_dir, logs_dir=logs_dir, timeout_seconds=timeout,
                    resume_existing=resume_existing,
                    python_executable=_run_setting(config, "python_executable", None),
                    skill_runtime=step_plan.get("runtime"),
                )
                inputs_match, input_reasons = store.verify_result_inputs(
                    step_id,
                    {"artifact_inputs": artifact_inputs},
                    expected_optional_inputs=expected_optional_inputs,
                )
                if not inputs_match:
                    message = (
                        "canonical artifact inputs changed while the step was "
                        f"running: {'; '.join(input_reasons)}"
                    )
                    result.update({
                        "success": False,
                        "error": message,
                        "artifact_inputs": artifact_inputs,
                        "producer_config_fingerprint": (
                            producer_config_fingerprint
                        ),
                    })
                    result_file = (
                        project_dir
                        / LEGACY_STEP_DIR[step_id]
                        / "result.json"
                    )
                    _atomic_write_json(result_file, result)
                    _write_artifact_lineage(logs_dir, store)
                    raise PipelineError(step_id, message, result)
            except PipelineError as step_exc:
                if step_plan["best_effort"]:
                    _append_error(logs_dir, f"{step_id} (best-effort) failed but pipeline continues: {step_exc}")
                    result = dict(step_exc.result or {})
                    result.update({
                        "success": False,
                        "step": step_id,
                        "error": str(step_exc),
                        "note": "best-effort step failure does not block pipeline success",
                        "artifact_inputs": artifact_inputs,
                        "producer_config_fingerprint": producer_config_fingerprint,
                    })
                    results[step_id] = result
                    result_file = project_dir / LEGACY_STEP_DIR[step_id] / "result.json"
                    processed_result_fingerprints[step_id] = _atomic_write_json(
                        result_file,
                        result,
                    )
                    _write_artifact_lineage(logs_dir, store)
                    _write_manifest(project_dir, f"{step_id}_failed", plan, results)
                else:
                    raise
            else:
                results[step_id] = result
                result["artifact_inputs"] = artifact_inputs
                result["producer_config_fingerprint"] = producer_config_fingerprint
                result_file = project_dir / LEGACY_STEP_DIR[step_id] / "result.json"
                produced = store.update_from(
                    step_id, result, source=str(result_file)
                )
                result["artifact_outputs"] = {
                    name: _artifact_fingerprint_record(store.provenance(name))
                    for name in produced
                }
                processed_result_fingerprints[step_id] = _atomic_write_json(
                    result_file,
                    result,
                )
                _write_artifact_lineage(logs_dir, store)
                _write_manifest(project_dir, f"{step_id}_complete", plan, results)

            processed_steps.append(step_id)
            next_step = (
                pending_steps[pending_index + 1]
                if pending_index + 1 < len(pending_steps)
                else None
            )

            if execution_mode == "interrupt" and next_step is not None:
                _write_checkpoint(
                    checkpoint_path,
                    state="awaiting_step_confirmation",
                    project_dir=project_dir,
                    plan=plan,
                    processed_steps=processed_steps,
                    processed_result_fingerprints=processed_result_fingerprints,
                    results=results,
                    producer_config_fingerprints=producer_config_fingerprints,
                    current_step=None,
                    next_step=next_step,
                )
                _write_manifest(
                    project_dir,
                    f"awaiting_confirmation_after_{step_id}",
                    plan,
                    results,
                )
                return _pause_result(
                    project_dir=project_dir,
                    plan=plan,
                    results=results,
                    completed_step=step_id,
                    next_step=next_step,
                    runtime_report_path=runtime_report_path,
                    checkpoint_path=checkpoint_path,
                    resume_rewinds=resume_rewinds,
                )

            pending_index += 1

        publication_failure = _finalize_publication_or_failure(
            project_dir=project_dir,
            plan=plan,
            results=results,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            producer_config_fingerprints=producer_config_fingerprints,
            checkpoint_path=checkpoint_path,
            store=store,
        )
        if publication_failure is not None:
            return publication_failure
        _write_checkpoint(
            checkpoint_path,
            state="complete",
            project_dir=project_dir,
            plan=plan,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            results=results,
            producer_config_fingerprints=producer_config_fingerprints,
            current_step=None,
            next_step=None,
        )
        _write_manifest(project_dir, "complete", plan, results)
        return _final_result(
            project_dir=project_dir,
            plan=plan,
            results=results,
            runtime_report_path=runtime_report_path,
            checkpoint_path=checkpoint_path,
            execution_mode=execution_mode,
            resume_rewinds=resume_rewinds,
        )
    except PipelineError as exc:
        _append_error(logs_dir, f"Step {exc.step} failed: {exc}")
        failure_result = dict(exc.result or {})
        failure_result.update({
            "success": False,
            "step": exc.step,
            "error": str(exc),
        })
        results[exc.step] = failure_result
        _write_checkpoint(
            checkpoint_path,
            state="failed",
            project_dir=project_dir,
            plan=plan,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            results=results,
            producer_config_fingerprints=producer_config_fingerprints,
            current_step=exc.step,
            next_step=exc.step,
            error=str(exc),
        )
        _write_manifest(project_dir, f"failed_{exc.step}", plan, results)
        return {
            "success": False,
            "status": "failed",
            "pipeline_complete": False,
            "requires_confirmation": False,
            "project_dir": str(project_dir), "failed_step": exc.step,
            "error": str(exc), "step_result": failure_result,
            "processed_steps": list(processed_steps),
            "steps_completed": _successful_steps(plan, results),
            "checkpoint_json": str(checkpoint_path),
            "manifest_json": str(project_dir / "manifest.json"),
            "manifest_md": str(project_dir / "manifest.md"),
            "artifact_lineage": str(logs_dir / "artifact_lineage.json"),
        }
# __ORCH_HELPERS__


def _restore_checkpoint(
    *,
    checkpoint_path: Path,
    project_dir: Path,
    plan: dict,
    merged: dict,
    store: ArtifactStore,
    producer_config_fingerprints: dict[str, str],
) -> dict:
    """Restore the longest still-valid processed-step prefix.

    The checkpoint is an orchestration cursor, not a replacement for artifact
    lineage. Every skipped step is revalidated against its current producer
    fingerprint plus the v2 input/output SHA-256 contracts. If a user changes
    an already-processed step (or an artifact is edited), continuation safely
    rewinds to the first affected step; changes confined to future steps remain
    allowed.
    """
    empty = {
        "processed_steps": [],
        "processed_result_fingerprints": {},
        "results": {},
        "rewind_reasons": [],
    }
    if not checkpoint_path.exists():
        return empty

    try:
        checkpoint = json.loads(checkpoint_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {checkpoint_path.name}: {exc}") from exc
    if not isinstance(checkpoint, dict):
        raise ValueError("checkpoint root must be a JSON object")
    if checkpoint.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError(
            "unsupported checkpoint schema_version "
            f"{checkpoint.get('schema_version')!r}"
        )
    if checkpoint.get("state") not in {
        "preflight", "running", "awaiting_step_confirmation", "failed", "complete",
    }:
        raise ValueError(f"unsupported checkpoint state {checkpoint.get('state')!r}")
    recorded_project_dir = checkpoint.get("project_dir")
    if recorded_project_dir != str(project_dir):
        raise ValueError(
            "checkpoint project_dir does not match the requested project directory"
        )

    recorded_steps = checkpoint.get("processed_steps")
    recorded_fingerprints = checkpoint.get("processed_step_fingerprints")
    recorded_result_fingerprints = checkpoint.get("processed_result_fingerprints")
    if not isinstance(recorded_steps, list) or not all(
        isinstance(step_id, str) for step_id in recorded_steps
    ):
        raise ValueError("processed_steps must be a list of step ids")
    if len(recorded_steps) != len(set(recorded_steps)):
        raise ValueError("processed_steps contains duplicates")
    if not isinstance(recorded_fingerprints, dict):
        raise ValueError("processed_step_fingerprints must be an object")
    if not isinstance(recorded_result_fingerprints, dict):
        raise ValueError("processed_result_fingerprints must be an object")

    execution_order = list(plan.get("execution_order") or [])
    reusable: list[str] = []
    reusable_result_fingerprints: dict[str, str] = {}
    results: dict[str, dict] = {}
    rewind_reasons: list[str] = []

    for index, step_id in enumerate(recorded_steps):
        if index >= len(execution_order) or execution_order[index] != step_id:
            expected = execution_order[index] if index < len(execution_order) else None
            rewind_reasons.append(
                f"checkpoint step {step_id!r} no longer matches execution-order "
                f"position {index} ({expected!r})"
            )
            break
        current_fingerprint = producer_config_fingerprints.get(step_id)
        if recorded_fingerprints.get(step_id) != current_fingerprint:
            rewind_reasons.append(
                f"{step_id} producer configuration or implementation changed"
            )
            break

        result_file = project_dir / LEGACY_STEP_DIR[step_id] / "result.json"
        try:
            result_payload = result_file.read_bytes()
        except OSError as exc:
            rewind_reasons.append(
                f"{step_id} result cannot be read: {exc}"
            )
            break
        active_result_fingerprint = _result_payload_fingerprint(result_payload)
        if recorded_result_fingerprints.get(step_id) != active_result_fingerprint:
            rewind_reasons.append(
                f"{step_id} result generation changed after its checkpoint"
            )
            break
        try:
            result = json.loads(result_payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            rewind_reasons.append(
                f"{step_id} result cannot be restored: {exc}"
            )
            break
        if not isinstance(result, dict):
            rewind_reasons.append(f"{step_id} result is not a JSON object")
            break
        if result.get("producer_config_fingerprint") != current_fingerprint:
            rewind_reasons.append(
                f"{step_id} result producer fingerprint does not match checkpoint"
            )
            break

        if result.get("success") is not True:
            if not plan["steps"][step_id].get("best_effort"):
                rewind_reasons.append(
                    f"{step_id} checkpoint result is not successful"
                )
                break
        else:
            result_ok, reasons = _checkpoint_result_is_current(
                step_id=step_id,
                result=result,
                merged=merged,
                store=store,
                producer_config_fingerprint=current_fingerprint,
            )
            if not result_ok:
                rewind_reasons.extend(
                    f"{step_id}: {reason}" for reason in reasons
                )
                break

        reusable.append(step_id)
        reusable_result_fingerprints[step_id] = active_result_fingerprint
        results[step_id] = result

    # A rewind makes stale current/downstream outputs ineligible immediately;
    # the normal execution loop will republish the first affected producer.
    if len(reusable) < len(recorded_steps):
        for step_id in execution_order[len(reusable):]:
            store.invalidate_producer(
                step_id,
                reason="interactive checkpoint rewound before this step",
            )

    return {
        "processed_steps": reusable,
        "processed_result_fingerprints": reusable_result_fingerprints,
        "results": results,
        "rewind_reasons": rewind_reasons,
    }


def _checkpoint_result_is_current(
    *,
    step_id: str,
    result: dict,
    merged: dict,
    store: ArtifactStore,
    producer_config_fingerprint: str | None,
) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if result.get("producer_config_fingerprint") != producer_config_fingerprint:
        reasons.append("result producer fingerprint does not match current plan")
    inputs_match, input_reasons = store.verify_result_inputs(
        step_id,
        result,
        expected_optional_inputs=_expected_optional_inputs(step_id, merged),
    )
    outputs_match, output_reasons = store.verify_result_outputs(step_id, result)
    reasons.extend(input_reasons)
    reasons.extend(output_reasons)

    if not reasons and inputs_match and outputs_match:
        # Re-publish in checkpoint order after config.input_artifacts has been
        # injected. The generic preload scan necessarily runs before injection,
        # so a valid derived result may have been conservatively rejected there
        # even though its recorded external input is now the active generation.
        store.update_from(
            step_id,
            result,
            source=f"checkpoint:{step_id}",
        )
        for name, recorded in (result.get("artifact_outputs") or {}).items():
            active = store.provenance(name)
            if active.get("producer") != step_id:
                reasons.append(
                    f"output {name!r} is not active for producer {step_id}"
                )
                continue
            if isinstance(recorded, dict) and (
                recorded.get("fingerprint") != active.get("fingerprint")
            ):
                reasons.append(f"output {name!r} active fingerprint changed")
    return inputs_match and outputs_match and not reasons, reasons


def _plan_fingerprint(plan: dict) -> str:
    """Hash the resolved plan while excluding UI-only continuation controls."""
    payload = copy.deepcopy(plan)
    project = payload.get("project") or {}
    for key in ("execution_mode", "project_dir", "resume_existing"):
        project.pop(key, None)
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _result_payload_fingerprint(payload: bytes) -> str:
    """Fingerprint the same immutable byte buffer used for JSON parsing."""
    digest = hashlib.sha256()
    digest.update(payload)
    return f"sha256:{digest.hexdigest()}"


def _merge_runtime_report(
    *,
    plan: dict,
    fresh_report: dict,
    existing_path: Path,
    processed_steps: list[str],
) -> dict:
    """Keep one current-plan runtime record across interrupt invocations."""
    existing_report: dict = {}
    if existing_path.exists():
        try:
            candidate = json.loads(existing_path.read_text())
            if isinstance(candidate, dict):
                existing_report = candidate
        except (OSError, json.JSONDecodeError):
            existing_report = {}

    def keyed(entries: Any) -> dict[tuple[str, str], dict]:
        out: dict[tuple[str, str], dict] = {}
        for entry in entries or []:
            if not isinstance(entry, dict):
                continue
            step_id = entry.get("step")
            tool_id = entry.get("tool")
            if isinstance(step_id, str) and isinstance(tool_id, str):
                out[(step_id, tool_id)] = entry
        return out

    existing = keyed(existing_report.get("tools"))
    fresh = keyed(fresh_report.get("tools"))
    desired: list[tuple[str, str]] = []
    for step_id in plan.get("execution_order") or []:
        for tool_id, tool in (plan["steps"][step_id].get("tools") or {}).items():
            if tool.get("enabled") and tool.get("runtime"):
                desired.append((step_id, tool_id))

    def required_for_preflight(key: tuple[str, str]) -> bool:
        runtime = plan["steps"][key[0]]["tools"][key[1]].get("runtime") or {}
        return runtime.get("required_for_preflight") is not False

    tools: list[dict] = []
    missing_evidence: list[str] = []
    for key in desired:
        entry = fresh.get(key) or existing.get(key)
        if entry is None:
            missing_evidence.append(f"{key[0]}.{key[1]}")
        else:
            tools.append(entry)

    missing = [
        f"{entry['step']}.{entry['tool']}"
        for entry in tools
        if not entry.get("available")
        and required_for_preflight((entry["step"], entry["tool"]))
    ]
    optional_missing = [
        f"{entry['step']}.{entry['tool']}"
        for entry in tools
        if not entry.get("available")
        and not required_for_preflight((entry["step"], entry["tool"]))
    ]
    for name in missing_evidence:
        step_id, tool_id = name.split(".", 1)
        target = missing if required_for_preflight((step_id, tool_id)) else optional_missing
        target.append(name)
    for tool_name in fresh_report.get("missing") or []:
        if tool_name not in missing:
            missing.append(tool_name)
    errors = [
        {
            "tool": f"{entry['step']}.{entry['tool']}",
            "detail": entry.get("detail"),
        }
        for entry in tools
        if not entry.get("available")
    ]
    errors.extend(
        {
            "tool": tool_name,
            "detail": "runtime evidence is missing from the cumulative report",
        }
        for tool_name in missing_evidence
    )
    for fresh_error in fresh_report.get("errors") or []:
        if fresh_error not in errors:
            errors.append(fresh_error)
    return {
        "success": not missing,
        "tools": tools,
        "missing": missing,
        "optional_missing": optional_missing,
        "errors": errors,
        "conda_executable": (
            fresh_report.get("conda_executable")
            or existing_report.get("conda_executable")
        ),
        "processed_steps_reused": list(processed_steps),
        "last_probe_steps": sorted(
            {entry["step"] for entry in fresh.values()}
        ),
        "updated_at": datetime.now().isoformat(timespec="seconds"),
    }


def _write_checkpoint(
    path: Path,
    *,
    state: str,
    project_dir: Path,
    plan: dict,
    processed_steps: list[str],
    processed_result_fingerprints: dict[str, str],
    results: dict[str, dict],
    producer_config_fingerprints: dict[str, str],
    current_step: str | None,
    next_step: str | None,
    error: str | None = None,
) -> None:
    scientific = _persist_scientific_outcome(project_dir, results)
    payload = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "state": state,
        "project_dir": str(project_dir),
        "plan_fingerprint": _plan_fingerprint(plan),
        "execution_order": list(plan.get("execution_order") or []),
        "processed_steps": list(processed_steps),
        "successful_steps": _successful_steps(plan, results),
        "processed_step_fingerprints": {
            step_id: producer_config_fingerprints[step_id]
            for step_id in processed_steps
        },
        "processed_result_fingerprints": {
            step_id: processed_result_fingerprints[step_id]
            for step_id in processed_steps
        },
        "current_step": current_step,
        "last_processed_step": processed_steps[-1] if processed_steps else None,
        "next_step": next_step,
        "error": error,
        **scientific,
        "updated_at": datetime.now().isoformat(timespec="seconds"),
    }
    _atomic_write_json(path, payload)


def _atomic_write_json(path: Path, payload: dict) -> str:
    """Publish JSON atomically and return the exact written-byte fingerprint."""
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(
        payload,
        indent=2,
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
        try:
            directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        except (AttributeError, OSError):
            directory_fd = None
        if directory_fd is not None:
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return _result_payload_fingerprint(encoded)


def _successful_steps(plan: dict, results: dict[str, dict]) -> list[str]:
    return [
        step_id
        for step_id in plan.get("execution_order", [])
        if (results.get(step_id) or {}).get("success") is True
    ]


def _scientific_outcome(results: dict[str, dict]) -> dict:
    """Normalize new and legacy coordinator conclusions without changing evidence.

    ``success`` elsewhere in this runner is execution success.  This helper
    deliberately keeps the independent nomination/scientific conclusion in a
    separate contract so an invalid protocol cannot look like a successful
    scientific nomination merely because every subprocess exited cleanly.
    """
    report = results.get("coordination")
    if not isinstance(report, dict) or report.get("success") is False:
        return {
            "scientific_status": "not_evaluated",
            "decision_outcome": None,
            "goal_met": None,
            "valid_for_nomination": None,
            "rollback_instructions": None,
        }

    if report.get("assessment_schema_version") == "round-assessment/v1":
        assessment = report.get("round_assessment")
        if not isinstance(assessment, dict):
            assessment = {}
        status = assessment.get("status")
        if status not in {"valid", "review_required", "invalid", "not_evaluated"}:
            status = "not_evaluated"
        valid = (
            status in {"valid", "review_required"}
            and assessment.get("valid_for_nomination") is True
        )
        return {
            "assessment_schema_version": "round-assessment/v1",
            "round_assessment": copy.deepcopy(assessment),
            "scientific_status": status,
            "decision_outcome": None,
            "goal_met": None,
            "valid_for_nomination": valid,
            "rollback_instructions": (
                copy.deepcopy(report.get("rollback_instructions"))
                if status == "invalid" else None
            ),
        }

    explicit_outcome = report.get("decision_outcome")
    outcome = (
        explicit_outcome
        if isinstance(explicit_outcome, str)
        and explicit_outcome in DECISION_OUTCOMES
        else _legacy_decision_outcome(report)
    )

    explicit_goal_met = report.get("goal_met")
    if isinstance(explicit_goal_met, bool):
        goal_met = explicit_goal_met
    elif outcome is None:
        goal_met = None
    else:
        goal_met = outcome == "success"

    explicit_valid = report.get("valid_for_nomination")
    if outcome != "success" or goal_met is False:
        valid_for_nomination = False if outcome is not None else None
    elif isinstance(explicit_valid, bool):
        valid_for_nomination = explicit_valid
    else:
        valid_for_nomination = _legacy_nomination_valid(report)

    rollback = None
    if outcome == "invalid":
        explicit_rollback = report.get("rollback_instructions")
        rollback = (
            copy.deepcopy(explicit_rollback)
            if isinstance(explicit_rollback, dict)
            else _legacy_rollback_instructions(report)
        )

    return {
        "scientific_status": outcome or "not_evaluated",
        "decision_outcome": outcome,
        "goal_met": goal_met,
        "valid_for_nomination": valid_for_nomination,
        "rollback_instructions": rollback,
    }


def _legacy_decision_outcome(report: dict) -> str | None:
    """Map the pre-structured coordinator decision onto the current enum."""
    decision = report.get("decision")
    if decision == "terminate_invalid":
        return "invalid"
    if decision == "iterate":
        return "iterate"
    if decision != "terminate":
        return None

    rationale = str(report.get("decision_rationale") or "").lower()
    if "maximum round" in rationale or "max round" in rationale:
        return "max_rounds"
    if (
        "no molecules amenable" in rationale
        or "no viable" in rationale
        or "no candidate" in rationale
    ):
        return "no_viable"

    summary = report.get("round_summary") or {}
    candidate_keep = summary.get("candidate_keep")
    candidate_modify = summary.get("candidate_modify")
    if (
        isinstance(candidate_keep, (int, float))
        and isinstance(candidate_modify, (int, float))
        and candidate_keep <= 0
        and candidate_modify <= 0
    ):
        return "no_viable"
    return "success"


def _legacy_nomination_valid(report: dict) -> bool:
    """Conservatively infer nomination validity from an old terminate report."""
    summary = report.get("round_summary") or {}
    parameters = report.get("parameters") or {}
    keep = summary.get("candidate_keep")
    total = summary.get("candidate_total")
    pass_rate = summary.get("candidate_pass_rate")
    if not isinstance(keep, (int, float)):
        keep = sum(
            1
            for record in report.get("candidate_molecules") or []
            if isinstance(record, dict)
            and record.get("classification") == "keep"
            and record.get("ranking_blocked") is not True
            and not str(record.get("nomination_status") or "").startswith("blocked")
        )
    if not isinstance(total, (int, float)):
        total = len(report.get("candidate_molecules") or [])
    if not isinstance(pass_rate, (int, float)):
        pass_rate = keep / total if total else 0.0

    min_candidates = parameters.get("min_candidates", 1)
    pass_rate_threshold = parameters.get("pass_rate_threshold", 0.0)
    if not isinstance(min_candidates, (int, float)):
        min_candidates = 1
    if not isinstance(pass_rate_threshold, (int, float)):
        pass_rate_threshold = 0.0

    require_quality_gates = parameters.get("require_quality_gates", True)
    quality_status = (
        report.get("quality_gate_status")
        or summary.get("quality_gate_status")
    )
    quality_ok = not require_quality_gates or quality_status == "keep"
    return bool(
        quality_ok
        and keep >= min_candidates
        and pass_rate >= pass_rate_threshold
    )


def _legacy_rollback_instructions(report: dict) -> dict:
    """Build a non-destructive rollback marker for legacy invalid reports."""
    summary = report.get("round_summary") or {}
    invalid_gates = list(summary.get("failed_quality_gates") or [])
    if not invalid_gates:
        for name, gate in (report.get("quality_gates") or {}).items():
            if not isinstance(gate, dict):
                continue
            status = gate.get("gate") or gate.get("status")
            if status == "reject" or gate.get("passed") is False:
                invalid_gates.append(str(name))
    return {
        "required": True,
        "scope": "protocol",
        "action": "rollback_invalid_protocol",
        "invalid_quality_gates": list(dict.fromkeys(map(str, invalid_gates))),
        "preserve_molecule_evidence": True,
    }


def _persist_scientific_outcome(
    project_dir: Path,
    results: dict[str, dict],
) -> dict:
    """Persist the current scientific conclusion without deleting artifacts."""
    outcome = _scientific_outcome(results)
    marker_path = project_dir / SCIENTIFIC_OUTCOME_FILENAME
    _atomic_write_json(
        marker_path,
        {
            "schema_version": 1,
            **outcome,
            "artifact_policy": "preserve",
            "invalidation_mode": "logical_marker_only",
            "non_destructive": True,
            "updated_at": datetime.now().isoformat(timespec="seconds"),
        },
    )
    return {
        **outcome,
        "scientific_outcome_json": str(marker_path),
    }


def _rewind_pause_result(
    *,
    project_dir: Path,
    plan: dict,
    results: dict[str, dict],
    rewind_target: str,
    checkpoint_path: Path,
    resume_rewinds: list[str],
) -> dict:
    last_completed = (
        next(
            (
                step_id
                for step_id in reversed(plan.get("execution_order") or [])
                if step_id in results
            ),
            None,
        )
    )
    scientific = _persist_scientific_outcome(project_dir, results)
    return {
        "success": True,
        "status": "awaiting_step_confirmation",
        "confirmation_reason": "checkpoint_rewind",
        "pipeline_complete": False,
        "requires_confirmation": True,
        "execution_mode": "interrupt",
        "project_dir": str(project_dir),
        "last_completed_step": last_completed,
        "last_step_success": (
            bool((results.get(last_completed) or {}).get("success"))
            if last_completed else None
        ),
        "next_step": rewind_target,
        "rewind_target": rewind_target,
        "processed_steps": [
            step_id
            for step_id in plan.get("execution_order", [])
            if step_id in results
        ],
        "steps_completed": _successful_steps(plan, results),
        "result_json": (
            str(project_dir / LEGACY_STEP_DIR[last_completed] / "result.json")
            if last_completed else None
        ),
        "checkpoint_json": str(checkpoint_path),
        "resume_project_dir": str(project_dir),
        "resume_rewinds": list(resume_rewinds),
        "resolved_run_plan": str(project_dir / "resolved_run_plan.json"),
        "tool_runtime_report": (
            str(project_dir / "tool_runtime_report.json")
            if (project_dir / "tool_runtime_report.json").exists()
            else None
        ),
        "manifest_json": str(project_dir / "manifest.json"),
        "manifest_md": str(project_dir / "manifest.md"),
        "artifact_lineage": str(project_dir / "logs" / "artifact_lineage.json"),
        **scientific,
    }


def _pause_result(
    *,
    project_dir: Path,
    plan: dict,
    results: dict[str, dict],
    completed_step: str,
    next_step: str,
    runtime_report_path: Path,
    checkpoint_path: Path,
    resume_rewinds: list[str],
) -> dict:
    result_file = project_dir / LEGACY_STEP_DIR[completed_step] / "result.json"
    scientific = _persist_scientific_outcome(project_dir, results)
    return {
        "success": True,
        "status": "awaiting_step_confirmation",
        "pipeline_complete": False,
        "requires_confirmation": True,
        "execution_mode": "interrupt",
        "project_dir": str(project_dir),
        "last_completed_step": completed_step,
        "last_step_success": bool((results.get(completed_step) or {}).get("success")),
        "next_step": next_step,
        "processed_steps": [
            step_id for step_id in plan["execution_order"] if step_id in results
        ],
        "steps_completed": _successful_steps(plan, results),
        "result_json": str(result_file),
        "checkpoint_json": str(checkpoint_path),
        "resume_project_dir": str(project_dir),
        "resume_rewinds": list(resume_rewinds),
        "resolved_run_plan": str(project_dir / "resolved_run_plan.json"),
        "tool_runtime_report": str(runtime_report_path),
        "manifest_json": str(project_dir / "manifest.json"),
        "manifest_md": str(project_dir / "manifest.md"),
        "artifact_lineage": str(project_dir / "logs" / "artifact_lineage.json"),
        **scientific,
    }


def _final_result(
    *,
    project_dir: Path,
    plan: dict,
    results: dict[str, dict],
    runtime_report_path: Path,
    checkpoint_path: Path,
    execution_mode: str,
    resume_rewinds: list[str],
) -> dict:
    last_step = (
        (plan.get("execution_order") or [None])[-1]
        if plan.get("execution_order")
        else None
    )
    interrupt = execution_mode == "interrupt" and last_step is not None
    scientific = _persist_scientific_outcome(project_dir, results)
    return {
        "success": True,
        "status": "awaiting_final_confirmation" if interrupt else "complete",
        "pipeline_complete": True,
        "requires_confirmation": interrupt,
        "execution_mode": execution_mode,
        "project_dir": str(project_dir),
        "last_completed_step": last_step,
        "last_step_success": (
            bool((results.get(last_step) or {}).get("success"))
            if last_step else None
        ),
        "result_json": (
            str(project_dir / LEGACY_STEP_DIR[last_step] / "result.json")
            if last_step else None
        ),
        "next_step": None,
        "processed_steps": [
            step_id for step_id in plan.get("execution_order", []) if step_id in results
        ],
        "steps_completed": _successful_steps(plan, results),
        "checkpoint_json": str(checkpoint_path),
        "resume_rewinds": list(resume_rewinds),
        "resolved_run_plan": str(project_dir / "resolved_run_plan.json"),
        "tool_runtime_report": str(runtime_report_path),
        "manifest_json": str(project_dir / "manifest.json"),
        "manifest_md": str(project_dir / "manifest.md"),
        "artifact_lineage": str(project_dir / "logs" / "artifact_lineage.json"),
        "round_report": _round_report_path(project_dir),
        "final_report": _project_final_report_path(project_dir),
        **scientific,
    }


def _finalize_publication_or_failure(
    *,
    project_dir: Path,
    plan: dict,
    results: dict[str, dict],
    processed_steps: list[str],
    processed_result_fingerprints: dict[str, str],
    producer_config_fingerprints: dict[str, str],
    checkpoint_path: Path,
    store: ArtifactStore,
) -> dict | None:
    """Finalize the deliverable, preserving completed scientific step results."""
    try:
        _finalize_pipeline_publication(project_dir, plan, results, store)
    except (OSError, RuntimeError, ValueError) as exc:
        message = f"final publication failed: {exc}"
        _append_error(project_dir / "logs", message)
        _write_checkpoint(
            checkpoint_path,
            state="failed",
            project_dir=project_dir,
            plan=plan,
            processed_steps=processed_steps,
            processed_result_fingerprints=processed_result_fingerprints,
            results=results,
            producer_config_fingerprints=producer_config_fingerprints,
            current_step=None,
            next_step=None,
            error=message,
        )
        _write_manifest(project_dir, "publication_failed", plan, results)
        scientific = _persist_scientific_outcome(project_dir, results)
        return {
            "success": False,
            "status": "publication_failed",
            "pipeline_complete": False,
            "requires_confirmation": False,
            "project_dir": str(project_dir),
            "failed_step": "publication",
            "error": message,
            "processed_steps": list(processed_steps),
            "steps_completed": _successful_steps(plan, results),
            "checkpoint_json": str(checkpoint_path),
            "manifest_json": str(project_dir / "manifest.json"),
            "manifest_md": str(project_dir / "manifest.md"),
            "artifact_lineage": str(
                project_dir / "logs" / "artifact_lineage.json"
            ),
            **scientific,
        }
    return None


def _preload_artifacts(
    project_dir: Path,
    upstream_project_dir: Path | None,
    merged: dict,
) -> ArtifactStore:
    """Scan this project dir and the upstream round for already-produced artifacts.

    Lets a resume / round-N run (steps disabled via start_step) satisfy the DAG
    from prior results instead of re-running upstream steps.
    """
    store = ArtifactStore()
    search_dirs = [d for d in (upstream_project_dir, project_dir) if d is not None]
    try:
        producer_fingerprints = _producer_config_fingerprints(merged)
    except (AttributeError, KeyError, TypeError, ValueError):
        # Validation below owns malformed-config diagnostics. Fail closed for
        # cache preload without letting fingerprint resolution mask them.
        producer_fingerprints = {
            step_id: "<INVALID_CONFIG>"
            for step_id in LEGACY_STEP_DIR
        }
    for base in search_dirs:
        for step_id, dir_name in LEGACY_STEP_DIR.items():
            result_file = base / dir_name / "result.json"
            if not result_file.exists():
                continue
            if not _preload_step_authorized(step_id, merged):
                store.reject_preload(
                    step_id,
                    source=str(result_file),
                    reason=(
                        "producer step is disabled without a start_step "
                        "resume boundary"
                    ),
                )
                continue
            try:
                result = json.loads(result_file.read_text())
            except (OSError, json.JSONDecodeError):
                store.reject_preload(
                    step_id,
                    source=str(result_file),
                    reason="newer result.json is unreadable or invalid JSON",
                )
                continue
            if not isinstance(result, dict) or result.get("success") is not True:
                store.reject_preload(
                    step_id,
                    source=str(result_file),
                    reason="newer result.json is not a successful step result",
                )
                continue
            store.preload_from(
                step_id,
                result,
                source=str(result_file),
                expected_producer_fingerprint=producer_fingerprints[step_id],
                expected_optional_inputs=_expected_optional_inputs(
                    step_id,
                    merged,
                ),
            )
    return store


def _preload_step_authorized(step_id: str, merged: dict) -> bool:
    """Allow current-cache preload or an explicit start-step resume boundary."""
    step = (merged.get("steps") or {}).get(step_id) or {}
    if not isinstance(step, dict):
        return False
    if step.get("enabled", False):
        return True
    start_step = merged.get("start_step")
    if start_step is None:
        return False
    try:
        return registry.STEPS[step_id]["order"] < float(start_step)
    except (TypeError, ValueError):
        return False


def _producer_config_fingerprints(merged: dict) -> dict[str, str]:
    """Hash each adapter-resolved scientific config without artifact inputs.

    Deliberately excludes only a provenance-*derived* ``candidate_source``.  It
    cannot be included: this fingerprint is what authorizes artifact preload,
    and the provenance it reads lives in the store that preload populates, so
    hashing that derived value would be self-referential.  A user-declared
    ``context.candidate_source`` is already present in ``merged`` and therefore
    is covered by ``_screening_policy`` just like an explicit ``screening_mode``;
    changing either invalidates the affected producer cache.
    """
    disabled = pconfig.deep_merge(merged, {
        "steps": {
            step_id: {"enabled": False}
            for step_id in registry.STEPS
        },
    })
    config_plan = pconfig.resolve(disabled)
    project_config = merged.get("project") or {}
    fingerprints: dict[str, str] = {}
    for step_id in LEGACY_STEP_DIR:
        spec = registry.STEPS[step_id]
        entrypoint = BIO_ROOT / spec["skill"] / "run.py"
        stable_runtime = {
            "project_dir": "<PROJECT_DIR>",
            "upstream_project_dir": project_config.get(
                "upstream_project_dir"
            ),
            "round_number": project_config.get("round_number", 1),
            "inferred_target_name": None,
        }
        producer_input = adapters.to_legacy_input(
            step_id,
            {},
            config_plan["steps"][step_id],
            config_plan.get("context") or {},
            stable_runtime,
        )
        payload = {
            "producer_contract_version": 2,
            "step_id": step_id,
            "skill": spec["skill"],
            "skill_entrypoint": str(entrypoint),
            "implementation_fingerprint": _implementation_fingerprint(
                spec["skill"],
                producer_input,
            ),
            "producer_input_without_artifacts": producer_input,
        }
        payload["local_file_signatures"] = _local_file_signatures({
            "producer_input_without_artifacts": producer_input,
            "skill_entrypoint": str(entrypoint),
        })
        canonical = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            default=str,
        ).encode("utf-8")
        fingerprints[step_id] = hashlib.sha256(canonical).hexdigest()
    return fingerprints


def _producer_config_fingerprint(
    step_id: str,
    merged: dict,
) -> str:
    """Return one adapter-scoped producer hash."""
    return _producer_config_fingerprints(merged)[step_id]


def _artifact_fingerprint_record(provenance: dict) -> dict:
    """Keep the stable subset needed to verify a published artifact."""
    return {
        key: value
        for key, value in provenance.items()
        if key in {
            "producer",
            "source",
            "fingerprint",
            "kind",
            "path",
            "size",
        }
    }


def _run_step(
    step_id: str,
    skill: str,
    step_input: dict,
    project_dir: Path,
    logs_dir: Path,
    timeout_seconds: int,
    resume_existing: bool,
    python_executable: str | None,
    skill_runtime: dict | None,
) -> dict:
    step_dir = project_dir / LEGACY_STEP_DIR[step_id]
    step_dir.mkdir(parents=True, exist_ok=True)
    result_file = step_dir / "result.json"
    input_fingerprint = _step_input_fingerprint(
        step_id=step_id,
        skill=skill,
        step_input=step_input,
        timeout_seconds=timeout_seconds,
        python_executable=python_executable,
        skill_runtime=skill_runtime,
    )
    cache_reuse_allowed = not (
        step_id == "pocket_analysis"
        and bool(
            (((step_input.get("tools") or {}).get("pocketminer") or {}).get("enabled"))
        )
    )
    if resume_existing and cache_reuse_allowed and result_file.exists():
        try:
            existing = json.loads(result_file.read_text())
        except (OSError, json.JSONDecodeError):
            existing = {}
        if (
            existing.get("success")
            and existing.get("input_fingerprint") == input_fingerprint
            and existing.get("cache_contract_version") == 2
        ):
            outputs_match, _ = ArtifactStore().verify_result_outputs(
                step_id, existing
            )
            if outputs_match:
                return existing

    input_file = logs_dir / f"{step_id}_{skill}_input.json"
    stdout_file = logs_dir / f"{step_id}_{skill}.stdout.log"
    stderr_file = logs_dir / f"{step_id}_{skill}.stderr.log"
    input_file.write_text(json.dumps(step_input, indent=2, ensure_ascii=False))

    step_args = [
        str(BIO_ROOT / skill / "run.py"), "--input-file", str(input_file),
        "--output-dir", str(step_dir),
    ]
    process_env = None
    if python_executable:
        cmd = [python_executable, *step_args]
    elif skill_runtime:
        step_runtime = {
            "environment": skill_runtime.get("environment", "current"),
            "command": "python",
        }
        cmd = tool_runtime.build_command(step_runtime, step_args)
        process_env = tool_runtime.build_process_env(
            step_runtime, cmd, overrides=None,
        )
    else:
        cmd = [sys.executable, *step_args]
    # Compare with the pre-launch inventory only on a failed attempt. Old
    # files retained after a rewind are not evidence from this execution.
    before_files = _step_file_inventory(step_dir)
    research_execution = {
        "schema_version": 1,
        "execution_id": uuid.uuid4().hex,
        "project_dir": str(project_dir.resolve()),
        "step_id": step_id,
        "skill": skill,
        "input_fingerprint": input_fingerprint,
        "step_input": copy.deepcopy(step_input),
        "input_file_signatures": _local_file_signatures(step_input),
    }

    def research_provenance() -> dict:
        # Retain the identity of this actual execution on cache reuse. Only
        # JSON reports created/changed by this attempt may supply new measured
        # evidence; stale reports left by previous attempts are excluded.
        current_files = _step_file_inventory(step_dir)
        changed_json = [
            path for path, identity in current_files.items()
            if before_files.get(path) != identity
            and Path(path).suffix.lower() == ".json"
            and Path(path) != result_file
        ]
        return {**research_execution,
                "evidence_file_signatures": _local_file_signatures(changed_json)}

    def failure_result(error: str, stdout: Any = "", stderr: Any = "") -> dict:
        def as_text(value: Any) -> str:
            return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else (value or "")

        stdout_file.write_text(as_text(stdout))
        stderr_file.write_text(as_text(stderr))
        payload = _parse_step_stdout(as_text(stdout)) or {}
        current = _step_file_inventory(step_dir)
        payload.update({
            "success": False, "error": error, "command": cmd,
            "timeout_seconds": timeout_seconds,
            "input_fingerprint": input_fingerprint, "cache_contract_version": 2,
            "research_provenance": research_provenance(),
        })
        files = payload.get("files")
        if not isinstance(files, dict):
            files = {}
        files["partial_attempt_files"] = [
            path for path, identity in current.items()
            if before_files.get(path) != identity and Path(path) != result_file
        ]
        payload["files"] = files
        _atomic_write_json(result_file, payload)
        return payload

    try:
        proc = run_command(cmd, timeout=timeout_seconds, env=process_env)
    except subprocess.TimeoutExpired as exc:
        result = failure_result(f"{skill} timed out after {timeout_seconds}s", exc.stdout, exc.stderr)
        raise PipelineError(step_id, result["error"], result)
    except (PipelineCancelled, KeyboardInterrupt) as exc:
        result = failure_result(str(exc) or "experiment interrupted",
                                getattr(exc, "stdout", ""), getattr(exc, "stderr", ""))
        # This result is known to belong to the current attempt. Do not
        # rediscover an arbitrary old result.json from an exception handler.
        manifest = _read_publication_manifest(project_dir)
        info = manifest.get("steps", {}).get(step_id)
        if isinstance(info, dict):
            info.update(status="interrupted", error=result["error"],
                        files=result["files"], result_json=str(result_file))
            _atomic_write_json(project_dir / "manifest.json", manifest)
        raise
    except OSError as exc:
        result = failure_result(f"cannot execute {skill}: {exc}")
        raise PipelineError(step_id, result["error"], result) from exc

    stdout_file.write_text(proc.stdout or "")
    stderr_file.write_text(proc.stderr or "")
    result = _parse_step_stdout(proc.stdout)
    if result is None:
        message = (
            f"{skill} produced an invalid or incomplete JSON result "
            f"(exit code {proc.returncode}); expected one complete object with boolean success"
        )
        result = failure_result(message, proc.stdout, proc.stderr)
        raise PipelineError(step_id, message, result)
    result["command"] = cmd
    result["timeout_seconds"] = timeout_seconds
    result["input_fingerprint"] = input_fingerprint
    result["cache_contract_version"] = 2
    result["cache_reuse_policy"] = (
        "disabled_for_enabled_pocketminer"
        if not cache_reuse_allowed
        else "resolved_input_local_file_and_output_sha256"
    )
    result["research_provenance"] = research_provenance()
    result_file.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=str))

    if proc.returncode != 0 or not result.get("success"):
        message = result.get("error") or (
            f"{skill} exited with code {proc.returncode}"
            if proc.returncode != 0
            else f"{skill} reported success=false (exit code 0)"
        )
        result = failure_result(message, proc.stdout, proc.stderr)
        raise PipelineError(step_id, message, result)
    return result


def _step_file_inventory(step_dir: Path) -> dict[str, tuple[int, int, int]]:
    inventory = {}
    for path in step_dir.rglob("*"):
        try:
            if path.is_file() and not path.is_symlink():
                stat = path.stat()
                inventory[str(path)] = (stat.st_ino, stat.st_mtime_ns, stat.st_size)
        except OSError:
            continue
    return inventory


def _step_input_fingerprint(
    *,
    step_id: str,
    skill: str,
    step_input: dict,
    timeout_seconds: int,
    python_executable: str | None,
    skill_runtime: dict | None,
) -> str:
    """Hash the resolved execution contract used to authorize step-result reuse.

    The adapter input contains every resolved tool switch, parameter, runtime,
    context value, and wired artifact path. Old results without this fingerprint
    are intentionally not reusable: otherwise toggling PocketMiner (or any
    other tool) could silently return an output generated under a different
    scientific configuration.
    """
    payload = {
        "cache_contract_version": 2,
        "step_id": step_id,
        "skill": skill,
        "skill_entrypoint": str(BIO_ROOT / skill / "run.py"),
        "implementation_fingerprint": _implementation_fingerprint(
            skill,
            step_input,
        ),
        "step_input": step_input,
        "timeout_seconds": int(timeout_seconds),
        "python_executable": python_executable,
        "skill_runtime": skill_runtime or {},
    }
    payload["local_file_signatures"] = _local_file_signatures(payload)
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


_SKILL_IMPLEMENTATION_DEPENDENCIES = {
    "docking": (BIO_ROOT / "structure-prep" / "scripts",),
    "pocket-analysis": (BIO_ROOT / "pipeline-entry" / "scripts",),
    "structure-evaluator": (
        BIO_ROOT / "pocket-analysis" / "scripts",
    ),
}
_IMPLEMENTATION_SUFFIXES = {".py", ".sh", ".yaml", ".yml", ".toml"}


def _implementation_fingerprint(
    skill: str,
    execution_payload: Any | None = None,
) -> str:
    """Hash local/transitive skill code plus declared external source trees."""
    roots = {
        BIO_ROOT / skill,
        SHARED_DIR,
        *_SKILL_IMPLEMENTATION_DEPENDENCIES.get(skill, ()),
        *_declared_source_directories(execution_payload),
    }
    digest = hashlib.sha256()
    for root in sorted(roots, key=lambda item: str(item)):
        digest.update(f"root:{root}".encode("utf-8"))
        digest.update(b"\0")
        if not root.is_dir():
            digest.update(b"<missing-directory>\0")
            continue
        source_files = sorted(
            (
                path
                for path in root.rglob("*")
                if path.is_file()
                and path.suffix.lower() in _IMPLEMENTATION_SUFFIXES
            ),
            key=lambda item: str(item),
        )
        for path in source_files:
            try:
                label = path.relative_to(root)
            except ValueError:
                label = path
            digest.update(str(label).encode("utf-8"))
            digest.update(b"\0")
            try:
                with path.open("rb") as handle:
                    for chunk in iter(
                        lambda: handle.read(1024 * 1024),
                        b"",
                    ):
                        digest.update(chunk)
            except OSError as exc:
                digest.update(
                    f"<unreadable:{type(exc).__name__}>".encode("utf-8")
                )
            digest.update(b"\0")
    return digest.hexdigest()


def _declared_source_directories(value: Any) -> set[Path]:
    """Collect enabled-tool runtime source trees without scanning user data."""
    directories: set[Path] = set()
    tools = value.get("tools") if isinstance(value, dict) else None
    if not isinstance(tools, dict):
        return directories
    for tool in tools.values():
        if not isinstance(tool, dict) or not tool.get("enabled"):
            continue
        runtime = tool.get("runtime")
        if not isinstance(runtime, dict):
            continue
        source_dir = runtime.get("source_dir")
        if isinstance(source_dir, str) and source_dir:
            directories.add(Path(source_dir).expanduser())
    return directories


def _use_optional_artifact(
    step_id: str,
    artifact: str,
    context: dict,
    step_params: dict | None = None,
) -> bool:
    """Resolve optional wiring from the scientific input actually in use."""
    if step_id == "pocket_analysis" and artifact == "docking_box":
        return (step_params or {}).get("site_selection_strategy") == "provided_box"
    if (
        step_id == "molecule_generation"
        and artifact == "inferred_target_name"
    ):
        return not bool(context.get("target_name"))
    return False


def _expected_optional_inputs(
    step_id: str,
    config: dict,
) -> set[str]:
    """Return optional lineage inputs required by the current run config."""
    context = config.get("context") or {}
    return {
        name
        for name in registry.STEPS[step_id].get("optional_consumes", [])
        if _use_optional_artifact(
            step_id, name, context,
            ((config.get("steps") or {}).get(step_id) or {}).get("params"),
        )
    }


def _local_file_signatures(value: Any) -> dict[str, dict[str, Any]]:
    """Hash local files (including TensorFlow checkpoint prefixes) in a payload."""
    discovered: dict[str, dict[str, Any]] = {}
    strings: list[str] = []

    def collect(item: Any) -> None:
        if isinstance(item, dict):
            for key in sorted(item):
                collect(item[key])
        elif isinstance(item, (list, tuple)):
            for child in item:
                collect(child)
        elif isinstance(item, str) and item:
            strings.append(item)

    def add_file(path: Path) -> None:
        try:
            resolved = path.expanduser().resolve(strict=True)
            if not resolved.is_file() or str(resolved) in discovered:
                return
        except (OSError, RuntimeError, ValueError):
            return
        digest = hashlib.sha256()
        try:
            with resolved.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            stat = resolved.stat()
        except OSError:
            return
        discovered[str(resolved)] = {
            "sha256": digest.hexdigest(),
            "size_bytes": stat.st_size,
        }

    collect(value)
    for text in strings:
        # The payload also contains prose, SMILES and other non-path strings.
        # Probe defensively without restricting legitimate relative filenames.
        try:
            path = Path(text).expanduser()
        except (OSError, RuntimeError, ValueError):
            continue
        add_file(path)
        # TensorFlow object checkpoints are addressed by a non-file prefix.
        index_path = Path(f"{path}.index")
        try:
            has_index = index_path.is_file()
        except (OSError, ValueError):
            has_index = False
        if has_index:
            add_file(index_path)
            try:
                data_files = sorted(path.parent.glob(f"{path.name}.data-*"))
            except OSError:
                data_files = []
            for data_file in data_files:
                add_file(data_file)
    return dict(sorted(discovered.items()))


def _parse_step_stdout(stdout: str) -> dict | None:
    """Accept one complete result, optionally preceded by non-JSON log lines.

    Never salvage an inner object from a truncated report: it may look valid
    while belonging to an entirely different schema (e.g. a metric result).
    """
    text = (stdout or "").strip()
    if not text:
        return None
    # Legacy tools may log before their JSON. Once a JSON document starts, its
    # entire suffix must parse; later braces are not alternate result starts.
    if not text.startswith(("{", "[")):
        lines = text.splitlines(keepends=True)
        start = next((i for i, line in enumerate(lines)
                      if line.lstrip().startswith(("{", "["))), None)
        if start is None:
            return None
        text = "".join(lines[start:]).strip()
    try:
        result = json.loads(text)
    except (json.JSONDecodeError, ValueError):
        return None
    if not isinstance(result, dict) or type(result.get("success")) is not bool:
        return None
    if "files" in result and not isinstance(result["files"], dict):
        return None
    if result.get("result_schema", "sbdd-step-result/v1") != "sbdd-step-result/v1":
        return None
    return result


_GENERATOR_PROVENANCE_PRODUCERS = frozenset({"molecule_generation"})
_OPERATOR_PROVENANCE_PRODUCERS = frozenset({"input_artifacts", "context"})


def _candidate_source(plan: dict, store: ArtifactStore | None = None) -> str:
    """Name where this round's molecules came from, for screening-mode defaults.

    The executable policy ties the level-2 elimination thresholds to provenance: a
    generative round can afford the strict mode because it can simply generate
    more, a sampled library gets the standard mode, and a set the user supplied
    in full must only be marked, never thinned.

    Provenance has to come from the artifact store, not from the step-enable
    matrix.  An iteration round starts at step 4 with the previous round's
    generator output preloaded (ITERATIVE_GUIDE.md:183,196,212-213), and
    ``apply_step_range`` disables ``molecule_generation`` to do it.  Reading only
    ``steps`` therefore classified every iteration round as
    ``user_supplied_full`` -- V1, marking only, zero elimination for the whole
    round -- when the molecules were in fact machine-generated and the policy
    asks for V3.  The store records who produced each artifact, so a library
    carried over from ``molecule_generation`` is distinguishable from one the
    user injected through ``input_artifacts``.
    """
    explicit = str(
        ((plan.get("context") or {}).get("candidate_source") or "")
    ).strip().lower()
    if explicit:
        allowed = {
            "molecule_generation", "sampled_library", "user_supplied_full",
        }
        if explicit not in allowed:
            raise ValueError(
                "context.candidate_source must be one of "
                + ", ".join(sorted(allowed))
            )
        return explicit

    steps = plan.get("steps") or {}
    generation = steps.get("molecule_generation") or {}
    if generation.get("enabled"):
        return "molecule_generation"

    if store is not None:
        # ``library`` first: it is the original set, so it answers "generated or
        # handed over?" directly.  ``prefiltered_library`` only says that
        # ``property_prefilter`` ran, which is true of either origin, so it can
        # never claim generation -- it falls through to the neutral V2 default.
        for artifact in ("library", "prefiltered_library"):
            producer = str(
                (store.provenance(artifact) or {}).get("producer") or ""
            )
            if not producer:
                continue
            if producer in _GENERATOR_PROVENANCE_PRODUCERS:
                return "molecule_generation"
            if producer in _OPERATOR_PROVENANCE_PRODUCERS:
                return "user_supplied_full"
            return "sampled_library"

    injected = plan.get("input_artifacts") or {}
    if isinstance(injected, dict) and (
        injected.get("library") or injected.get("prefiltered_library")
    ):
        return "user_supplied_full"
    return "sampled_library"


def _run_setting(config: dict, key: str, default):
    """Read a run/project setting from config.project first, then top-level.

    The `project` block is the documented home for these; top-level is accepted
    as a fallback for backward compatibility and CLI overrides (--project-dir).
    """
    project = config.get("project") or {}
    if key in project and project[key] is not None:
        return project[key]
    if key in config and config[key] is not None:
        return config[key]
    return default


def _project_dir(config: dict) -> Path:
    explicit = _run_setting(config, "project_dir", None)
    if explicit:
        return Path(explicit).expanduser().resolve()

    round_number = int(_run_setting(config, "round_number", 1))
    upstream = _run_setting(config, "upstream_project_dir", None)
    if upstream:
        return (Path(upstream).expanduser().resolve().parent / f"round{round_number}")

    slug = _slugify(_run_setting(config, "project_slug", None)
                    or (config.get("context") or {}).get("pdb_id")
                    or (config.get("context") or {}).get("target_name")
                    or "sbdd")
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_id = uuid.uuid4().hex[:8]
    data_root = _run_setting(config, "data_root", None)
    if data_root:
        root = Path(data_root).expanduser()
    else:
        root = deployment_path("${SBDD_DATA_ROOT}")
    return (
        root / f"{timestamp}-{run_id}-{slug}" / f"round{round_number}"
    ).resolve()


def initialize_project(config: dict, *, request_id: str | None = None) -> dict:
    """Reserve a fresh round before input conversion, without running a plan.

    Only project naming facts are needed here. The authoritative scientific
    configuration is assembled after normalization and passed to ``run`` with
    the returned directory as a CLI override.
    """
    try:
        if not isinstance(config, dict):
            raise ValueError("project initialization input must be an object")
        if not isinstance(config.get("project", {}), dict):
            raise ValueError("project must be an object")
        round_number = _run_setting(config, "round_number", 1)
        if isinstance(round_number, bool) or not isinstance(round_number, int) or round_number < 1:
            raise ValueError("project.round_number must be a positive integer")
        project_dir = _project_dir(config)
        naming = {
            "round_number": round_number,
            "project_slug": _slugify(_run_setting(config, "project_slug", None)
                                     or (config.get("context") or {}).get("pdb_id")
                                     or (config.get("context") or {}).get("target_name") or "sbdd"),
            "upstream_project_dir": _run_setting(config, "upstream_project_dir", None),
        }
        return project_workspace.reserve(
            project_dir, request_id=request_id, naming=naming,
            data_root=_run_setting(config, "data_root", None),
            explicit_directory=bool(_run_setting(config, "project_dir", None)),
        )
    except (OSError, ValueError, TypeError) as exc:
        return {
            "success": False,
            "status": "project_initialization_failed",
            "pipeline_invoked": False,
            "pipeline_complete": False,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _slugify(value: Any) -> str:
    from io_utils import slugify
    return slugify(value)


def _step_timeout(step_id: str, config: dict) -> int:
    return pipeline_timeouts.resolve_step_timeout(step_id, config)


def _infer_target_name(step0_result: dict) -> str | None:
    metadata = (step0_result or {}).get("metadata") or {}
    for chain in metadata.get("chains") or []:
        if chain.get("description"):
            return chain["description"]
    return None


def _write_env_paths(project_dir: Path, logs_dir: Path, config: dict) -> None:
    pose_config = (
        (config.get("steps") or {}).get("pose_evaluation") or {}
    )
    pose_timeout: int | str = (
        _step_timeout("pose_evaluation", config)
        if pose_config.get("enabled", False)
        else "DISABLED"
    )
    payload = [
        f"BIO_SKILLS_ROOT={BIO_ROOT}",
        f"PROJECT_DIR={project_dir}",
        f"PYTHON={_run_setting(config, 'python_executable', None) or sys.executable}",
        f"STEP4_TIMEOUT_SECONDS={pose_timeout}",
    ]
    (logs_dir / "env_paths.txt").write_text("\n".join(payload) + "\n")


def _recorded_session_owners(provenance: dict) -> list[str]:
    """Return the Session ids already recorded on a Round, primary first."""
    candidates: list = [provenance.get("drugmate_session_id")]
    recorded = provenance.get("drugmate_session_ids")
    if isinstance(recorded, list):
        candidates.extend(recorded)
    owners: list[str] = []
    for value in candidates:
        session_id = value.strip() if isinstance(value, str) else ""
        if (
            session_id
            and len(session_id) <= _MAX_SESSION_ID_LENGTH
            and session_id not in owners
        ):
            owners.append(session_id)
    return owners


def _manifest_execution_provenance(project_dir: Path) -> dict:
    """Return durable runtime provenance without consulting pipeline config."""
    manifest_path = project_dir / "manifest.json"
    provenance: dict = {}
    if manifest_path.exists():
        try:
            previous_manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError):
            previous_manifest = None
        if isinstance(previous_manifest, dict):
            previous_provenance = previous_manifest.get("execution_provenance")
            if isinstance(previous_provenance, dict):
                provenance = copy.deepcopy(previous_provenance)

    # Keep every Session that has driven this Round, in order, and leave the
    # first one in place as the primary owner. Rewriting the primary when a
    # later Session resumes would move the Round out of the panel of the
    # Session that created it, which is what made panels disappear.
    owners = _recorded_session_owners(provenance)

    # Delegated agents have their own execution Session id, but SBDD assets
    # belong to the top-level logical Session that initiated the work. The
    # terminal bridge injects this stable owner ContextVar into subprocesses.
    # Keep DRUGMATE_SESSION_ID as a compatibility fallback for direct CLI runs.
    drugmate_session_id = (
        os.environ.get("DRUGMATE_SBDD_OWNER_SESSION_ID", "").strip()
        or os.environ.get("DRUGMATE_SESSION_ID", "").strip()
    )
    if drugmate_session_id and drugmate_session_id not in owners:
        owners.append(drugmate_session_id)

    if owners:
        provenance["drugmate_session_id"] = owners[0]
        provenance["drugmate_session_ids"] = owners
    else:
        provenance.pop("drugmate_session_id", None)
        provenance.pop("drugmate_session_ids", None)
    return provenance


def _write_manifest(project_dir: Path, status: str, plan: dict, results: dict) -> None:
    scientific = _persist_scientific_outcome(project_dir, results)
    execution_provenance = _manifest_execution_provenance(project_dir)
    manifest = {
        "project_dir": str(project_dir), "status": status,
        **scientific,
        "updated_at": datetime.now().isoformat(timespec="seconds"), "steps": {},
    }
    coordination = results.get("coordination") or {}
    reporting = results.get("reporting") or {}
    published_dir = (
        _verified_publication(project_dir)
        if status == "complete"
        and coordination.get("success") is True
        and reporting.get("success") is True else None
    )
    if published_dir is not None:
        manifest["published_dir"] = str(published_dir)
        manifest["publication_mode"] = "pipeline_final"
    if execution_provenance:
        manifest["execution_provenance"] = execution_provenance
    for step_id, step_plan in plan["steps"].items():
        step_dir = project_dir / LEGACY_STEP_DIR[step_id]
        result = results.get(step_id)
        if result and result.get("success"):
            state = "done"
        elif result is not None:
            state = "failed_best_effort" if step_plan["best_effort"] else "failed"
        elif not step_plan["enabled"]:
            state = "disabled"
        else:
            state = "pending"
        manifest["steps"][step_id] = {
            "skill": step_plan["skill"], "order": step_plan["order"],
            "dir": LEGACY_STEP_DIR[step_id], "status": state,
            "result_json": str(step_dir / "result.json"),
            "files": (result or {}).get("files", {}),
            "error": (result or {}).get("error"),
        }
        if status == f"running_{step_id}" and result is None:
            manifest["steps"][step_id]["status"] = "running"
    _atomic_write_json(project_dir / "manifest.json", manifest)
    if published_dir is None:
        _refresh_partial_view(project_dir, manifest)
    _atomic_write_json(project_dir / "manifest.json", manifest)
    _write_manifest_md(project_dir, manifest)


def _read_publication_manifest(project_dir: Path) -> dict:
    """Read publication metadata without making an error path fail again."""
    try:
        return _load_json_object(project_dir / "manifest.json")
    except ValueError:
        return {}


def _refresh_partial_view(project_dir: Path, manifest: dict) -> None:
    for key in ("published_dir", "publication_mode", "publication_error"):
        manifest.pop(key, None)
    try:
        partial_publication.publish_snapshot(project_dir, manifest)
    except Exception as exc:
        # Publishing diagnostics must not mask the scientific failure or turn
        # a successful step into a failed one. Retire an obsolete reader view.
        manifest["publication_error"] = f"partial publication failed: {exc}"
        try:
            _quarantine_canonical_publication(project_dir, reason=manifest["publication_error"])
            _append_error(project_dir / "logs", manifest["publication_error"])
        except (OSError, RuntimeError):
            pass
    else:
        manifest["published_dir"] = str(project_dir / PUBLISHED_DIRNAME)
        manifest["publication_mode"] = "pipeline_partial"


def _publish_diagnostic_view(
    project_dir: Path, *, status: str, error: str, preserve_steps: bool,
) -> None:
    """Also cover early validation errors and unexpected Python exceptions."""
    manifest = _read_publication_manifest(project_dir) if preserve_steps else {}
    if not manifest:
        manifest = {
            "project_dir": str(project_dir), "steps": {},
            "scientific_status": "not_evaluated", "goal_met": None,
            "valid_for_nomination": None,
            "scientific_outcome_json": str(project_dir / SCIENTIFIC_OUTCOME_FILENAME),
        }
    if status == "failed" and str(manifest.get("status", "")).startswith("failed_"):
        status = manifest["status"]
    manifest.update(status=status, error=error, updated_at=datetime.now().isoformat(timespec="seconds"))
    for info in manifest.get("steps", {}).values():
        if info.get("status") == "running":
            info.update(status="interrupted" if status == "interrupted" else "failed",
                        error=error, files={}, result_json=None)
    try:
        _refresh_partial_view(project_dir, manifest)
        _atomic_write_json(project_dir / "manifest.json", manifest)
        _write_manifest_md(project_dir, manifest)
    except OSError:
        # The original exception/status remains authoritative if the disk is
        # unavailable, including while trying to save the publication error.
        pass


def _write_manifest_md(project_dir: Path, manifest: dict) -> None:
    lines = [
        "# Project Manifest", "",
        f"- Execution status: {manifest['status']}",
        f"- Scientific status: `{manifest['scientific_status']}`",
        f"- Goal met: `{manifest['goal_met']}`",
        f"- Valid for nomination: `{manifest['valid_for_nomination']}`",
        f"- Scientific outcome marker: `{manifest['scientific_outcome_json']}`",
        f"- Updated: {manifest['updated_at']}",
        f"- Project directory: `{manifest['project_dir']}`", "",
        "| Order | Step | Skill | Status | Result |",
        "|---:|---|---|---|---|",
    ]
    for step_id, info in sorted(manifest["steps"].items(), key=lambda kv: kv[1]["order"]):
        lines.append(f"| {info['order']} | {step_id} | {info['skill']} | {info['status']} | `{info['result_json']}` |")
    lines.append("")
    (project_dir / "manifest.md").write_text("\n".join(lines))


def _write_artifact_lineage(logs_dir: Path, store: ArtifactStore) -> Path:
    """Persist active artifact provenance and every invalidation decision."""
    path = logs_dir / "artifact_lineage.json"
    existing_history: list[dict] = []
    if path.exists():
        try:
            existing = json.loads(path.read_text())
            if isinstance(existing, dict) and isinstance(
                existing.get("history"), list
            ):
                existing_history = [
                    event
                    for event in existing["history"]
                    if isinstance(event, dict)
                ]
        except (OSError, json.JSONDecodeError):
            existing_history = []

    current_history = store.history()
    persisted_count = int(
        getattr(store, "_pipeline_lineage_persisted_count", 0)
    )
    new_events = current_history[persisted_count:]
    cumulative = list(existing_history)
    for local_event in new_events:
        event = dict(local_event)
        event["event"] = len(cumulative) + 1
        cumulative.append(event)
    setattr(
        store,
        "_pipeline_lineage_persisted_count",
        len(current_history),
    )
    _atomic_write_json(
        path,
        {
            "active": store.provenance(),
            "history": cumulative,
        },
    )
    return path


def _append_error(logs_dir: Path, message: str) -> None:
    with (logs_dir / "error_log.md").open("a") as handle:
        handle.write(f"- {datetime.now().isoformat(timespec='seconds')} {message}\n")


def _round_report_path(project_dir: Path) -> str | None:
    candidate = project_dir / LEGACY_STEP_DIR["coordination"] / "round_report.md"
    return str(candidate) if candidate.exists() else None


def _project_final_report_path(project_dir: Path) -> str | None:
    project_root = project_dir.parent if re.fullmatch(r"round\d+", project_dir.name) else project_dir
    candidate = project_root / "final_report.md"
    return str(candidate) if candidate.exists() else None


def _archive_cli_input_config(result: dict, archived_config: dict) -> None:
    """Keep the invocation's ``--project-dir`` out of the round config.

    The CLI override selects a reserved or existing checkpoint directory; it is
    not part of the Agent's scientific config and would point a copied round-N+1
    config back at round N. ``run`` has already created the round, so this
    replaces only that config snapshot with the unmodified input-file payload.
    """
    project_dir = result.get("project_dir") if isinstance(result, dict) else None
    if not project_dir or not result.get("pipeline_config"):
        return
    config_path = Path(str(project_dir)) / RUN_CONFIG_FILENAME
    _atomic_write_json(config_path, archived_config)
    result["pipeline_config"] = str(config_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="SBDD pipeline orchestrator")
    parser.add_argument("--input", type=str, help="JSON input string")
    parser.add_argument("--input-file", type=str, help="Path to JSON input file")
    parser.add_argument("--project-dir", type=str, help="Override project_dir")
    parser.add_argument("--request-id", type=str,
                        help="Stable initialization request identity; reuse on retries (init only)")
    parser.add_argument(
        "--init-project-dir", action="store_true",
        help="Reserve a fresh round and its input directories; do not validate or execute a plan",
    )
    args = parser.parse_args()
    if args.request_id is not None and not args.init_project_dir:
        parser.error("--request-id is only for --init-project-dir; execute using its returned --project-dir")

    if args.input:
        config = json.loads(args.input)
    elif args.input_file:
        config = json.loads(Path(args.input_file).read_text())
    else:
        parser.error("Must provide --input or --input-file")

    archived_config = copy.deepcopy(config)
    if args.project_dir:
        config.setdefault("project", {})["project_dir"] = args.project_dir

    result = initialize_project(config, request_id=args.request_id) if args.init_project_dir else run(config)
    if args.project_dir and not args.init_project_dir:
        _archive_cli_input_config(result, archived_config)
    print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
    sys.exit(0 if result.get("success") else 1)


if __name__ == "__main__":
    main()

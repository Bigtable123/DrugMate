"""Persistence and hard bounds for Agent-authored R1–R3 campaign plans.

This module deliberately does not choose experiments, rank conditions, or decide
whether a scientific objective has been met. Configs without a campaign retain
the standalone runner contract. All paths in the campaign contract are absolute.
"""
from __future__ import annotations

import copy
import fcntl
import hashlib
import json
import math
import os
import tempfile
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

STAGES = ("R1", "R2", "R3")
AUTONOMY_LIMITS = {"A0": "R3", "A1": "R2", "A2": "R3", "A3": "R3"}
LATE_STAGE_TOOLS = {"openmm_md", "mdanalysis", "mmgbsa", "openfe_rbfe", "ties_rbfe"}
PAUSED_STATUSES = {
    "awaiting_step_confirmation", "awaiting_final_confirmation", "interrupted",
    "waiting_resources", "resource_unavailable", "awaiting_redocking",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _integer(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 1


def _text(value) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _absolute(value) -> bool:
    return _text(value) and Path(value).expanduser().is_absolute()


def _load(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return data


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


@contextmanager
def _state_lock(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_name(path.name + ".lock").open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def config_hash(config: dict) -> str:
    """Bind scientific inputs while allowing execution-mode/resume CLI changes."""
    value = copy.deepcopy(config)
    for block in (value, value.get("project", {})):
        for key in ("project_dir", "resume_existing", "execution_mode"):
            block.pop(key, None)
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False, allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _input_hashes(config: dict) -> dict:
    """Freeze caller-supplied input content, not merely its mutable filename."""
    hashes = {}
    for section in ("input_artifacts", "context"):
        for key, value in (config.get(section) or {}).items():
            if section == "context" and not key.endswith("_file"):
                continue
            if not isinstance(value, str):
                continue
            path = Path(value).expanduser()
            if path.is_file():
                hashes[f"{section}.{key}"] = {"path": str(path.resolve()), "sha256": _file_hash(path)}
    return hashes


def initialize_campaign(path: Path, *, autonomy: str = "A1",
                        authorized_through: str | None = None, max_runs: int = 6,
                        authorization_ref: str | None = None) -> dict:
    """Record session authorization supplied by the Agent; never infer approval.

    Without a reference this creates a proposed intent, which cannot execute.
    The skill must check that a supplied reference corresponds to actual user
    authorization; a local JSON file cannot authenticate conversational consent.
    """
    if autonomy not in AUTONOMY_LIMITS:
        raise ValueError("autonomy must be A0, A1, A2, or A3")
    limit = authorized_through or ("R1" if autonomy == "A0" else AUTONOMY_LIMITS[autonomy])
    if limit not in STAGES or STAGES.index(limit) > STAGES.index(AUTONOMY_LIMITS[autonomy]):
        raise ValueError(f"{autonomy} cannot authorize through {limit}")
    if not _integer(max_runs):
        raise ValueError("max_runs must be a positive integer")
    state = {"schema_version": 1, "autonomy": autonomy, "authorized_through": limit,
             "max_runs": max_runs, "authorization_ref": authorization_ref,
             "authorization_status": "recorded" if _text(authorization_ref) else "proposed",
             "created_at": _now(), "runs": []}
    path = Path(path).expanduser().resolve()
    with _state_lock(path):
        if path.exists():
            raise ValueError(f"campaign state already exists: {path}; refusing to erase history")
        _atomic_json(path, state)
    return state


def _state_errors(state: dict) -> list[str]:
    errors = []
    autonomy = state.get("autonomy")
    stage = state.get("authorized_through")
    if state.get("schema_version") != 1:
        errors.append("unsupported campaign state schema_version")
    if autonomy not in AUTONOMY_LIMITS or stage not in STAGES:
        errors.append("campaign requires valid autonomy and authorized_through")
    elif STAGES.index(stage) > STAGES.index(AUTONOMY_LIMITS[autonomy]):
        errors.append(f"{autonomy} cannot authorize through {stage}")
    if not _integer(state.get("max_runs")):
        errors.append("campaign max_runs must be a positive integer")
    if not isinstance(state.get("runs"), list) or any(not isinstance(r, dict) for r in state.get("runs", [])):
        errors.append("campaign runs must be a list of objects")
    else:
        seen = set()
        for run in state["runs"]:
            key = (run.get("stage"), run.get("attempt"))
            if run.get("stage") not in STAGES or not _integer(run.get("attempt")):
                errors.append("campaign history contains an invalid stage/attempt")
                continue
            elif key in seen:
                errors.append("campaign history contains duplicate stage/attempt")
            seen.add(key)
            if not _absolute(run.get("project_dir")) or not _text(run.get("config_hash")):
                errors.append("campaign history requires project_dir and config_hash")
            if run.get("run_id") != f"{key[0]}-a{key[1]}" or run.get("status") not in {"running", "paused", "complete", "failed"}:
                errors.append("campaign history contains an invalid run_id/status")
    return errors


def authorize_campaign(path: Path, *, authorization_ref: str, autonomy: str | None = None,
                       authorized_through: str | None = None, max_runs: int | None = None) -> dict:
    """Record a newly supplied session authorization without deleting run history."""
    if not _text(authorization_ref):
        raise ValueError("authorization_ref must refer to actual user authorization")
    path = Path(path).expanduser().resolve()
    with _state_lock(path):
        state = _load(path)
        previous = {key: state.get(key) for key in
                    ("autonomy", "authorized_through", "max_runs", "authorization_ref")}
        for key, value in (("autonomy", autonomy), ("authorized_through", authorized_through),
                           ("max_runs", max_runs)):
            if value is not None:
                state[key] = value
        errors = _state_errors(state)
        if not errors and state["max_runs"] < len(state["runs"]):
            errors.append("max_runs cannot be lower than already reserved runs")
        if errors:
            raise ValueError("; ".join(errors))
        state.setdefault("authorization_history", []).append({"previous": previous, "recorded_at": _now()})
        state.update(authorization_ref=authorization_ref, authorization_status="recorded")
        _atomic_json(path, state)
        return state


def confirm_final(path: Path, *, run_id: str, authorization_ref: str) -> dict:
    """Record final user confirmation of a finished interrupt-mode run, without rerunning."""
    if not _text(authorization_ref):
        raise ValueError("final confirmation requires a user authorization_ref")
    path = Path(path).expanduser().resolve()
    with _state_lock(path):
        state = _load(path)
        errors = _state_errors(state)
        if errors:
            raise ValueError("; ".join(errors))
        run = state["runs"][-1] if state["runs"] else {}
        if run.get("run_id") != run_id:
            raise ValueError("final confirmation must target the latest campaign run")
        if run.get("status") == "complete" and run.get("final_confirmation"):
            return copy.deepcopy(run)
        if run.get("status") != "paused" or run.get("last_status") != "awaiting_final_confirmation":
            raise ValueError("run is not awaiting final confirmation")
        project_dir = Path(run["project_dir"]).resolve()
        artifact = project_dir / "campaign_result.json"
        record = _load(artifact)
        result = record.get("result") or {}
        if result.get("success") is not True or result.get("pipeline_complete") is not True or result.get("status") != "awaiting_final_confirmation":
            raise ValueError("final confirmation requires a successfully completed pipeline result")
        checkpoint_path = Path(result.get("checkpoint_json") or "").resolve()
        if not checkpoint_path.is_relative_to(project_dir) or not checkpoint_path.is_file():
            raise ValueError("final confirmation requires this run's checkpoint")
        checkpoint = _load(checkpoint_path)
        if (checkpoint.get("state") != "complete" or checkpoint.get("next_step") is not None or
                not {"coordination", "reporting"}.issubset(set(checkpoint.get("successful_steps") or []))):
            raise ValueError("checkpoint does not verify completed coordination and reporting")
        confirmation = {"authorization_ref": authorization_ref, "confirmed_at": _now()}
        record["final_confirmation"] = confirmation
        run.update(status="complete", finished_at=_now(), final_confirmation=confirmation)
        _atomic_json(artifact, record)
        _atomic_json(path, state)
        return copy.deepcopy(run)


def abandon_run(path: Path, *, run_id: str, reason: str) -> dict:
    """Close a paused/crashed attempt after proving no runner owns its directory.

    This does not terminate processes. A live runner must first be stopped using
    its normal cancellation path; lock acquisition is deliberately nonblocking.
    """
    if not _text(reason):
        raise ValueError("abandon requires a reason")
    path = Path(path).expanduser().resolve()
    with _state_lock(path):
        state = _load(path)
        errors = _state_errors(state)
        if errors:
            raise ValueError("; ".join(errors))
        run = state["runs"][-1] if state["runs"] else {}
        if run.get("run_id") != run_id or run.get("status") not in {"running", "paused"}:
            raise ValueError("abandon must target the latest unfinished campaign run")
        project_dir = Path(run["project_dir"])
        with (project_dir / ".pipeline.lock").open("a+") as handle:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError("runner is still active; stop it normally before abandoning this attempt") from exc
            try:
                artifact = project_dir / "campaign_result.json"
                record = _load(artifact) if artifact.is_file() else {
                    "schema_version": 1, "run_id": run_id,
                    "result": {"success": False, "pipeline_complete": False, "status": "abandoned"},
                }
                termination = {"reason": reason, "abandoned_at": _now(), "previous_status": run["status"]}
                record["termination"] = termination
                run.update(status="failed", last_status="abandoned", finished_at=_now(),
                           termination=termination, result_file=str(artifact))
                _atomic_json(artifact, record)
                _atomic_json(path, state)
                return copy.deepcopy(run)
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _basic_errors(config: dict) -> list[str]:
    if not isinstance(config, dict) or (config.get("project") is not None and not isinstance(config.get("project"), dict)):
        return ["campaign config and project must be objects"]
    project = config.get("project") or {}
    if not project.get("campaign_state_file"):
        if any(project.get(k) is not None for k in ("campaign_stage", "campaign_attempt", "round_decision_file")):
            return ["campaign fields require project.campaign_state_file"]
        return []
    errors = []
    for field in ("campaign_state_file", "round_decision_file"):
        if not _absolute(project.get(field)):
            errors.append(f"project.{field} must be an absolute path")
    if project.get("campaign_stage") not in STAGES:
        errors.append("project.campaign_stage must be R1, R2, or R3")
    if not _integer(project.get("campaign_attempt")):
        errors.append("project.campaign_attempt must be a positive integer")
    if not isinstance(config.get("steps"), dict):
        return errors + ["campaign steps must be an object"]
    for step in ("coordination", "reporting"):
        if (config.get("steps", {}).get(step) or {}).get("enabled") is not True:
            errors.append(f"campaign requires steps.{step}.enabled=true")
    start, stop = _start(config), config.get("stop_step")
    if (isinstance(start, bool) or not isinstance(start, (float, int)) or
            not math.isfinite(start) or start > 7):
        errors.append("campaign start_step must include coordination (order 7)")
    if stop is not None and (isinstance(stop, bool) or not isinstance(stop, (float, int)) or
                             not math.isfinite(stop) or stop < 8):
        errors.append("campaign stop_step must include reporting (order 8)")
    for step in (config.get("steps") or {}).values():
        if not isinstance(step, dict) or step.get("enabled") is not True:
            continue
        for tool_id, tool in (step.get("tools") or {}).items():
            if tool_id in LATE_STAGE_TOOLS and isinstance(tool, dict) and tool.get("enabled") is True:
                errors.append(f"{tool_id} requires a later-stage workflow; this campaign supports R1–R3")
    return errors


def _start(config: dict):
    value = config.get("start_step", (config.get("project") or {}).get("start_step", 0))
    return 0 if value is None else value


def _seed_values(config: dict) -> set:
    seeds = set()
    context_seed = (config.get("context") or {}).get("seed")
    docking = (config.get("steps") or {}).get("docking") or {}
    if docking.get("enabled"):
        # Use the registry's actual parameter defaults for selected engines.
        try:
            import pipeline_registry as registry
            specs = registry.STEPS["docking"]["tools"]
        except ImportError:
            specs = {}
        for name, tool in (docking.get("tools") or {}).items():
            if not isinstance(tool, dict) or tool.get("enabled") is not True:
                continue
            seed = (tool.get("params") or {}).get("seed", ((specs.get(name) or {}).get("params", {}).get("seed") or {}).get("default"))
            if isinstance(seed, int) and not isinstance(seed, bool):
                seeds.add(seed)
    if not seeds and isinstance(context_seed, int) and not isinstance(context_seed, bool):
        seeds.add(context_seed)
    return seeds


def _plan_errors(config: dict, plan: dict, *, bound: bool) -> list[str]:
    errors = []
    project = config["project"]
    stage = project["campaign_stage"]
    if plan.get("schema_version") != 1:
        errors.append("round decision requires schema_version=1")
    if plan.get("stage") != stage or plan.get("attempt") != project["campaign_attempt"]:
        errors.append("round decision stage/attempt do not match config")
    if bound and plan.get("config_hash") != config_hash(config):
        errors.append("round decision config_hash does not match config; regenerate the reviewed plan")
    if bound and plan.get("input_artifact_hashes") != _input_hashes(config):
        errors.append("round decision input artifact content changed; regenerate and review the plan")
    if not _text(plan.get("reason")):
        errors.append("round decision requires a reason")
    for field in ("evidence_refs", "parameter_changes", "candidate_ids"):
        if not isinstance(plan.get(field), list):
            errors.append(f"round decision {field} must be a list")
    for change in plan.get("parameter_changes", []) if isinstance(plan.get("parameter_changes"), list) else []:
        if not isinstance(change, dict) or not _text(change.get("path")) or "before" not in change or "after" not in change:
            errors.append("parameter_changes entries require path, before, and after")
    start = plan.get("start_step")
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start) or start != _start(config):
        errors.append("round decision start_step must match config start_step (default 0)")
    if not isinstance(plan.get("molecules_changed"), bool):
        errors.append("round decision molecules_changed must be a boolean")
    purposes = {"R1": {"initial_screening", "protocol_repair"},
                "R2": {"replication", "parameter_improvement", "protocol_repair"},
                "R3": {"molecular_optimization", "protocol_repair"}}
    if plan.get("purpose") not in purposes[stage]:
        errors.append(f"invalid purpose for {stage}")
    injected = config.get("input_artifacts") or {}
    generation = (config.get("steps") or {}).get("molecule_generation") or {}
    new_r3_input = stage == "R3" and (generation.get("enabled") is True or bool(injected.get("library")))
    if stage == "R3" and injected.get("prefiltered_library"):
        errors.append("R3 cannot inject prefiltered_library; inject library and run property_prefilter")
    if plan.get("molecules_changed") or new_r3_input:
        if not isinstance(start, (int, float)) or start > 3:
            errors.append("changed/new molecules must restart at property_prefilter (order 3) or earlier")
        prefilter = (config.get("steps") or {}).get("property_prefilter") or {}
        if prefilter.get("enabled") is not True:
            errors.append("changed/new molecules require property_prefilter enabled")
    if stage == "R2":
        design = plan.get("replicate_design")
        if not isinstance(design, dict):
            errors.append("R2 requires replicate_design with replicate_id, seed and rerun_steps")
        else:
            if not _text(design.get("replicate_id")):
                errors.append("R2 replicate_design requires replicate_id")
            seed = design.get("seed")
            if isinstance(seed, bool) or not isinstance(seed, int) or seed not in _seed_values(config):
                errors.append("R2 replicate_design.seed must match an enabled docking tool seed or context.seed")
            rerun = design.get("rerun_steps")
            if not isinstance(rerun, list) or not rerun or any(not isinstance(s, str) or (config.get("steps", {}).get(s) or {}).get("enabled") is not True for s in rerun):
                errors.append("R2 replicate_design.rerun_steps must name enabled steps")
    return errors


def _history_errors(config: dict, plan: dict, state: dict, project_dir: Path | None,
                    *, enforce_authorization: bool = True) -> list[str]:
    errors = _state_errors(state)
    if errors:
        return errors
    project = config["project"]
    stage, attempt = project["campaign_stage"], project["campaign_attempt"]
    if enforce_authorization and (state.get("authorization_status") != "recorded" or not _text(state.get("authorization_ref"))):
        errors.append("campaign authorization is proposed; record an existing user authorization reference before execution")
    if enforce_authorization and STAGES.index(stage) > STAGES.index(state["authorized_through"]):
        errors.append(f"{stage} exceeds authorized_through={state['authorized_through']}")
    if enforce_authorization and state["autonomy"] == "A0" and not _text(plan.get("authorization_ref")):
        errors.append("A0 requires a per-run user authorization_ref in the round decision")
    runs = state["runs"]
    existing = next((r for r in runs if r.get("stage") == stage and r.get("attempt") == attempt), None)
    if existing:
        if runs[-1] is not existing:
            errors.append("cannot resume a superseded campaign run")
        if existing.get("last_status") == "abandoned":
            errors.append("cannot resume an abandoned campaign run; create a new attempt")
        if existing.get("status") == "paused" and existing.get("last_status") == "awaiting_final_confirmation":
            errors.append("run awaits final confirmation; use campaign.py confirm-final without rerunning tools")
        if existing["config_hash"] != config_hash(config):
            errors.append("stage/attempt already reserved with a different config")
        if project_dir is not None and existing["project_dir"] != str(project_dir):
            errors.append("stage/attempt already reserved in a different project_dir")
        previous = runs[runs.index(existing) - 1] if runs.index(existing) else None
    else:
        if enforce_authorization and len(runs) >= state["max_runs"]:
            errors.append("campaign max_runs budget exhausted")
        if project_dir is not None and any(r["project_dir"] == str(project_dir) for r in runs):
            errors.append("project_dir already belongs to another campaign run")
        previous = runs[-1] if runs else None
        if (enforce_authorization and state["autonomy"] == "A0" and previous and
                plan.get("authorization_ref") == previous.get("authorization_ref")):
            errors.append("A0 requires a new per-run authorization_ref")
        if previous is None:
            if stage != "R1" or attempt != 1:
                errors.append("campaign must begin with R1 attempt 1")
        else:
            if previous.get("status") not in {"complete", "failed"}:
                errors.append("previous campaign run is unfinished; resume it before planning another run")
            if stage == previous["stage"]:
                if attempt != previous["attempt"] + 1:
                    errors.append("same-stage retry must increment attempt by exactly one")
            elif STAGES.index(stage) != STAGES.index(previous["stage"]) + 1 or attempt != 1:
                errors.append("campaign may only advance one stage, starting at attempt 1")
            elif previous.get("status") != "complete":
                errors.append("a failed run cannot advance the campaign stage; repair within the same stage")
    refs = plan.get("evidence_refs") or []
    referenced = set()
    for ref in refs:
        if not isinstance(ref, dict) or not _absolute(ref.get("path")) or not _text(ref.get("run_id")):
            errors.append("evidence_refs entries require run_id and absolute path")
            continue
        run = next((r for r in runs if r.get("run_id") == ref["run_id"]), None)
        if run is None or run is existing or run.get("status") not in {"complete", "failed"}:
            errors.append(f"evidence references an unavailable prior run: {ref['run_id']}")
            continue
        path = Path(ref["path"]).resolve()
        if not path.is_relative_to(Path(run["project_dir"]).resolve()) or not path.is_file():
            errors.append(f"evidence must be an existing artifact within {run['run_id']}: {path}")
            continue
        if ref.get("sha256") != _file_hash(path):
            errors.append(f"evidence hash mismatch: {path}")
        referenced.add(run["run_id"])
    if previous and previous["run_id"] not in referenced:
        errors.append("round decision must reference evidence from the immediately preceding run")
    if previous and previous.get("protocol_rejected") and plan.get("purpose") != "protocol_repair":
        errors.append("previous protocol rejection requires an explicit protocol_repair plan")
    if previous and previous.get("protocol_rejected") and stage != previous["stage"]:
        errors.append("a rejected protocol must be repaired within the same campaign stage before advancing")
    if stage == "R2" and previous and plan.get("purpose") == "replication":
        design = plan.get("replicate_design") or {}
        if design.get("seed") in previous.get("seed_values", []) and not _text(design.get("independent_replicate_reason")):
            errors.append("unchanged replication seed requires independent_replicate_reason; round numbering is not independence")
        stop = config.get("stop_step")
        if "docking" not in design.get("rerun_steps", []) or _start(config) > 4 or (stop is not None and stop < 4):
            errors.append("replication must actually rerun docking (order 4)")
        if not existing and project.get("resume_existing", config.get("resume_existing", True)) is not False:
            errors.append("a new replication run requires project.resume_existing=false")
    return errors


def validate_campaign_config(config: dict, *, runtime: bool = False) -> list[str]:
    """Validate static fields; runtime also checks bound plan, history and consent."""
    try:
        errors = _basic_errors(config)
        if errors:
            return errors
        project = config.get("project") or {}
        if not project.get("campaign_state_file") or not runtime:
            return errors
        plan = _load(Path(project["round_decision_file"]).expanduser())
        errors = _plan_errors(config, plan, bound=True)
        if errors:
            return errors
        state = _load(Path(project["campaign_state_file"]).expanduser())
        explicit_dir = project.get("project_dir") or config.get("project_dir")
        return _history_errors(config, plan, state, Path(explicit_dir).expanduser().resolve() if explicit_dir else None)
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        return [f"campaign validation failed: {exc}"]


def bind_round_decision(config: dict, decision: dict) -> dict:
    """Validate an Agent-authored decision and bind it to config/artifact hashes."""
    errors = _basic_errors(config)
    if not (config.get("project") or {}).get("campaign_state_file"):
        errors.append("planning requires project.campaign_state_file")
    if errors:
        raise ValueError("; ".join(errors))
    plan = copy.deepcopy(decision)
    errors = _plan_errors(config, plan, bound=False)
    if errors:
        raise ValueError("; ".join(errors))
    for ref in plan["evidence_refs"]:
        if not isinstance(ref, dict) or not _absolute(ref.get("path")):
            raise ValueError("evidence_refs entries require absolute paths")
        ref["sha256"] = _file_hash(Path(ref["path"]).expanduser())
    plan["config_hash"] = config_hash(config)
    plan["input_artifact_hashes"] = _input_hashes(config)
    plan["created_at"] = _now()
    state = _load(Path(config["project"]["campaign_state_file"]).expanduser())
    explicit_dir = config["project"].get("project_dir") or config.get("project_dir")
    directory = Path(explicit_dir).expanduser().resolve() if explicit_dir else None
    errors = _history_errors(config, plan, state, directory, enforce_authorization=False)
    if errors:
        raise ValueError("; ".join(errors))
    blockers = _history_errors(config, plan, state, directory)
    plan["execution_authorized"] = not blockers
    plan["execution_blockers"] = blockers
    return plan


def reserve_campaign_run(config: dict, project_dir: Path) -> dict | None:
    """Reserve one budget unit under a file lock, idempotently on resume."""
    project = config.get("project") or {}
    if not project.get("campaign_state_file"):
        return None
    path = Path(project["campaign_state_file"]).expanduser().resolve()
    project_dir = Path(project_dir).expanduser().resolve()
    with _state_lock(path):
        errors = _basic_errors(config)
        if errors:
            raise ValueError("; ".join(errors))
        state = _load(path)
        plan = _load(Path(project["round_decision_file"]).expanduser())
        errors = _plan_errors(config, plan, bound=True)
        if not errors:
            errors = _history_errors(config, plan, state, project_dir)
        if errors:
            raise ValueError("; ".join(errors))
        run_id = f"{project['campaign_stage']}-a{project['campaign_attempt']}"
        existing = next((r for r in state["runs"] if r["run_id"] == run_id), None)
        archive = project_dir / "round_decision.json"
        snapshot = project_dir / "profile_snapshot.json"
        if existing:
            if not archive.is_file() or _load(archive) != plan:
                raise ValueError("reserved round decision was changed or its archive is missing")
            if existing["status"] in {"failed", "paused"}:
                _resume_execution(path, state, existing, project_dir)
            return copy.deepcopy(existing)
        if snapshot.exists() and _load(snapshot).get("config_hash") != config_hash(config):
            raise ValueError("project_dir contains a different frozen campaign config")
        _atomic_json(archive, plan)
        _atomic_json(snapshot, {"schema_version": 1, "run_id": run_id,
                               "config_hash": config_hash(config), "config": config,
                               "input_artifact_hashes": plan["input_artifact_hashes"],
                               "frozen_at": _now()})
        run = {"run_id": run_id, "stage": project["campaign_stage"],
               "attempt": project["campaign_attempt"], "status": "running",
               "project_dir": str(project_dir), "config_hash": config_hash(config),
               "round_decision_file": str(archive), "profile_snapshot": str(snapshot),
               "authorization_ref": plan.get("authorization_ref") or state.get("authorization_ref"),
               "seed_values": sorted(_seed_values(config)), "started_at": _now()}
        state["runs"].append(run)
        _atomic_json(path, state)
        return copy.deepcopy(run)


def _resume_execution(state_path: Path, state: dict, run: dict, project_dir: Path) -> None:
    """Called under the campaign lock after all resume checks have passed."""
    result_path = project_dir / "campaign_result.json"
    previous = _load(result_path) if result_path.is_file() else None
    history_path = project_dir / "execution_history" / f"{run['run_id']}-{uuid.uuid4().hex}.json"
    now = _now()
    _atomic_json(history_path, {
        "schema_version": 1, "archived_at": now,
        "run": copy.deepcopy(run), "campaign_result": previous,
    })
    run.setdefault("execution_history", []).append(str(history_path))
    run.update(status="running", last_status="running", updated_at=now,
               resumed_at=now, result_file=str(result_path))
    run.pop("finished_at", None)
    # A current result must not continue advertising the preceding failure.
    # Its immutable archive survives even if either atomic write is interrupted.
    _atomic_json(result_path, {
        "schema_version": 1, "run_id": run["run_id"],
        "protocol_rejected": run.get("protocol_rejected", False),
        "previous_execution_file": str(history_path), "updated_at": now,
        "result": {
            "success": False, "status": "running", "pipeline_complete": False,
            "requires_confirmation": False, "project_dir": str(project_dir),
            "checkpoint_json": str(project_dir / "pipeline_checkpoint.json"),
            "manifest_json": str(project_dir / "manifest.json"),
        },
    })
    _atomic_json(state_path, state)


def _protocol_rejected(result: dict) -> bool:
    assessment = result.get("round_assessment") or {}
    gates = assessment.get("protocol_gates") or result.get("protocol_gates") or {}
    return (assessment.get("protocol_status") == "reject" or
            assessment.get("status") == "invalid" or
            result.get("scientific_status") == "invalid" or
            result.get("scientific_outcome") == "invalid" or
            (isinstance(gates, dict) and gates.get("overall_status") == "reject"))


def complete_campaign_run(config: dict, project_dir: Path, result: dict) -> None:
    """Persist completion/failure; confirmations and transient pauses stay open."""
    project = config.get("project") or {}
    if not project.get("campaign_state_file"):
        return
    path = Path(project["campaign_state_file"]).expanduser().resolve()
    project_dir = Path(project_dir).expanduser().resolve()
    with _state_lock(path):
        state = _load(path)
        errors = _state_errors(state)
        if errors:
            raise ValueError("; ".join(errors))
        run = next((r for r in state["runs"] if r["project_dir"] == str(project_dir)), None)
        if run is None or run["config_hash"] != config_hash(config):
            raise ValueError("cannot complete an unreserved or differently configured campaign run")
        if run is not state["runs"][-1]:
            raise ValueError("cannot modify a superseded campaign run")
        status = result.get("status")
        paused = status in PAUSED_STATUSES or result.get("requires_confirmation") is True
        complete = result.get("pipeline_complete") is True and result.get("success") is True and status == "complete"
        # Once complete, accidental repeat callbacks cannot downgrade its evidence.
        if run["status"] == "complete":
            return
        run["status"] = "complete" if complete else "paused" if paused else "failed"
        run["last_status"] = status
        run["updated_at"] = _now()
        run["protocol_rejected"] = run.get("protocol_rejected", False) or _protocol_rejected(result)
        if not paused:
            run["finished_at"] = _now()
        else:
            run.pop("finished_at", None)
        artifact = project_dir / "campaign_result.json"
        _atomic_json(artifact, {"schema_version": 1, "run_id": run["run_id"],
                                "protocol_rejected": run["protocol_rejected"], "result": result})
        run["result_file"] = str(artifact)
        _atomic_json(path, state)

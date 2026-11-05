"""Pure-stdlib aggregator: reads Step 0-6 (and optional cross-round) artifacts
and returns a single report_data dict handed to the HTML template and renderers.

No RDKit / matplotlib / PyMOL imports here — this module always runs.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from pose_scores import (
    CANDIDATE_RANKING_PREFERENCE,
    LEGACY_SCORE_SPECS,
    collect_pose_scores,
    find_score,
    resolve_primary_metric,
)
from metric_scoring import SCORING_SCHEMA_VERSION
from candidate_review import build_candidate_reviews
from docking_artifact_contract import (
    REDOCKING_SEALED_FIELDS,
    evaluation_docking_lineage_error,
    redocking_measurement_error,
)


_COHORT_COMPONENT_FIELDS = (
    "protein_state_id", "pocket_id", "round_id", "protocol_profile_id",
    "metric_catalog_version",
)

_DOCKING_SCORE_METRIC_IDS = frozenset(
    str(spec["metric_id"]) for spec in LEGACY_SCORE_SPECS.values()
)
_DOCKING_SCORE_METHODS = frozenset(
    str(spec["method"]).casefold() for spec in LEGACY_SCORE_SPECS.values()
)

_AIZYNTH_ATTEMPTED_STATUSES = frozenset({
    "solved", "not_found_within_budget", "error",
})
_AIZYNTH_NOT_ATTEMPTED_STATUSES = frozenset({
    "disabled", "not_eligible", "not_selected_top_k",
})
_AIZYNTH_STATUS_ORDER = (
    "solved",
    "not_found_within_budget",
    "error",
    "not_selected_top_k",
    "not_eligible",
    "disabled",
    "missing",
)
_AIZYNTH_DEFAULT_REASONS = {
    "solved": "solved_route_found",
    "not_found_within_budget": "no_solved_route_within_search_budget",
    "error": "tool_execution_failed_without_recorded_detail",
    "not_selected_top_k": "outside_configured_top_k",
    "not_eligible": "excluded_by_pre_synthesis_eligibility_gate",
    "disabled": "disabled_by_config",
    "missing": "synthesis_record_missing",
}
_AIZYNTH_STATUS_MEANINGS = {
    "solved": "Attempt completed and at least one solved route was found.",
    "not_found_within_budget": (
        "Attempt completed without a solved route in the configured budget; "
        "this is a search outcome, not a tool failure or automatic reject."
    ),
    "error": (
        "The selected tool call failed; metrics are N/A and protocol review "
        "is required, not automatic molecule rejection."
    ),
    "not_selected_top_k": (
        "Not attempted because the row was outside the configured compute "
        "slice; synthesis metrics are N/A."
    ),
    "not_eligible": (
        "Not attempted because a pre-synthesis eligibility rule excluded the "
        "row; this status is not a tool failure."
    ),
    "disabled": (
        "Not attempted because AiZynthFinder was disabled; synthesis metrics "
        "are N/A."
    ),
    "missing": (
        "No usable synthesis status was recorded; do not infer an attempt, "
        "failure, or molecule verdict."
    ),
}


def _is_docking_pose_score(score: dict[str, Any]) -> bool:
    """Identify evidence that may only be sourced from committed Step 3.

    Older Step 4/6 artifacts duplicate docking values, sometimes without a
    reliable ``stage`` tag.  Accepting those copies after the Step 3 artifact
    contract failed would bypass the hash/identity checks performed by
    :func:`_read_step3` and could resurrect a score belonging to another pose.
    """
    stage = str(score.get("stage") or "").strip().casefold()
    metric_id = str(score.get("metric_id") or "").strip()
    method = str(score.get("method") or "").strip().casefold()
    return (
        stage == "docking"
        or metric_id in _DOCKING_SCORE_METRIC_IDS
        or method in _DOCKING_SCORE_METHODS
    )


# ---------- top-level entrypoints ----------

def aggregate_round(
    project_dir: Path,
    round_number: int | None = None,
    top_n: int = 10,
) -> dict[str, Any]:
    """Aggregate one round's step0..step6 artifacts.

    Missing files are recorded in `missing_sources[]` but never raise — the
    report should still render with whatever is present.
    """
    project_dir = Path(project_dir).expanduser().resolve()
    resolved_round = round_number if round_number is not None else _infer_round_number(project_dir)

    data: dict[str, Any] = {
        "mode": "round",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "project_dir": str(project_dir),
        "project_root": str(project_dir.parent),
        "round_number": resolved_round,
        "missing_sources": [],
        "step0": _read_step0(project_dir),
        "step1": _read_step1(project_dir),
        "step3": _read_step3(project_dir),
        "step4": _read_step4(project_dir),
        "step5": _read_step5(project_dir),
        "step6": _read_step6(project_dir),
    }
    _enforce_step3_step4_lineage(data)

    data["header"] = _build_header(data)
    data["candidates"] = _build_candidates(data, top_n=top_n)
    step6 = data["step6"]
    data["candidate_reviews"] = build_candidate_reviews(
        step6.get("candidate_ranking_json") or []
        if step6.get("available") and step6.get("success") else [],
        source_round=resolved_round,
    )
    data["synthesis_summary"] = _build_synthesis_summary(data.get("step5"))
    # ``identity_index`` is a transient join accelerator.  Several of its maps
    # use tuple keys, so exposing it in the returned report payload makes the
    # otherwise JSON-native sidecar fail during serialization.  Candidate
    # enrichment is complete at this point; do not leak the internal index into
    # the public report_data contract.
    step5 = data.get("step5")
    if isinstance(step5, dict):
        step5.pop("identity_index", None)
    data["candidate_ranking_contract"] = _candidate_ranking_contract(top_n)
    data["quality_gates"] = _build_quality_gates(data)
    data["screening"] = _build_screening(data)
    data["missing_sources"] = _collect_missing(data)
    return data


def aggregate_final(
    project_root: Path,
    top_n: int = 10,
) -> dict[str, Any]:
    """Aggregate all rounds under {project_root}/round*/ for the final report."""
    project_root = Path(project_root).expanduser().resolve()
    rounds = _discover_rounds(project_root)

    per_round: list[dict[str, Any]] = []
    for round_number, round_dir in rounds:
        per_round.append(aggregate_round(round_dir, round_number=round_number, top_n=top_n))

    candidate_source_round, candidate_source_reason = _pick_candidate_source_round(per_round)

    data: dict[str, Any] = {
        "mode": "final",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "project_dir": str(project_root),
        "project_root": str(project_root),
        "round_number": None,
        "rounds": per_round,
        "candidate_source_round": candidate_source_round,
        "candidate_source_reason": candidate_source_reason,
        "trend": _build_trend(per_round),
    }

    # Candidate rows bind only to a nomination-valid source round.  When none
    # exists, retain the last round's target/QC context for diagnosis but do not
    # silently promote its invalid or incomplete candidates.
    src = _find_round_by_number(per_round, candidate_source_round)
    context_src = src or (per_round[-1] if per_round else None)
    review_src = src or next(
        (r for r in reversed(per_round) if (r.get("candidate_reviews") or {}).get("groups")),
        None,
    )
    data["candidate_reviews"] = (
        review_src["candidate_reviews"] if review_src is not None
        else build_candidate_reviews([], source_round=None)
    )
    if context_src is not None:
        data["header"] = _build_header(
            context_src, override_round_label="final"
        )
        data["candidates"] = (
            _build_candidates(src, top_n=top_n, nomination_only=True)
            if src is not None else []
        )
        data["candidate_ranking_contract"] = _candidate_ranking_contract(top_n)
        data["quality_gates"] = (
            context_src.get("quality_gates")
            or _build_quality_gates(context_src)
        )
        data["screening"] = (
            context_src.get("screening") or _build_screening(context_src)
        )
        data["synthesis_summary"] = context_src.get("synthesis_summary") or (
            _build_synthesis_summary(context_src.get("step5"))
        )
        data["step0"] = context_src.get("step0", {})
        data["step1"] = context_src.get("step1", {})
    else:
        data["header"] = {"note": "no rounds discovered"}
        data["candidates"] = []
        data["candidate_ranking_contract"] = _candidate_ranking_contract(top_n)
        data["quality_gates"] = []
        data["screening"] = {"available": False}
        data["synthesis_summary"] = {
            "available": False,
            "reason": "no rounds discovered",
        }
        data["step0"] = {}
        data["step1"] = {}

    data["missing_sources"] = _collect_missing(data)
    return data


# ---------- step readers (silent on missing) ----------

def _read_step0(project_dir: Path) -> dict[str, Any]:
    step_dir = project_dir / "step0"
    if not step_dir.exists():
        return {"available": False, "reason": "step0 directory missing"}
    result_json = _safe_json(step_dir / "result.json")
    files = (result_json.get("files") if isinstance(result_json, dict) else {}) or {}
    prepared_pdb = _first_existing(
        files.get("prepared_pdb"),
        files.get("receptor_file"),
        _first_matching(step_dir, "*_receptor_prepared.pdb"),
    )
    ligand_sdf = _first_existing(
        files.get("ligand_sdf"),
        _first_matching(step_dir, "*_ligand_*.sdf"),
        _first_matching(step_dir, "*ligand*.sdf"),
    )
    prepared_noh = _first_existing(
        files.get("prepared_noH_pdb"),
        _first_matching(step_dir, "*_noH.pdb"),
    )
    return {
        "available": bool(prepared_pdb),
        "prepared_pdb": prepared_pdb,
        "prepared_noH_pdb": prepared_noh,
        "reference_ligand_sdf": ligand_sdf,
        "metadata": (result_json or {}).get("metadata") if isinstance(result_json, dict) else None,
        "result_json": _existing_or_none(step_dir / "result.json"),
    }


def _read_step1(project_dir: Path) -> dict[str, Any]:
    step_dir = project_dir / "step1"
    if not step_dir.exists():
        return {"available": False, "reason": "step1 directory missing"}
    pocket_report = _safe_json(step_dir / "pocket_analysis_report.json")
    blueprint = _safe_json(step_dir / "interaction_blueprint.json")
    pocket_features = _safe_json(step_dir / "pocket_features.json")
    interaction_fp = _safe_json(step_dir / "interaction_fingerprint.json")

    consensus_center = None
    consensus_residues: list[str] = []
    if isinstance(pocket_report, dict):
        detection = pocket_report.get("pocket_detection") or {}
        consensus = detection.get("consensus_pocket") or {}
        consensus_center = consensus.get("center")
        consensus_residues = list(consensus.get("residues") or [])

    blueprint_validation = (
        (blueprint.get("validation") or {}) if isinstance(blueprint, dict) else {}
    )
    blueprint_not_applicable = bool(
        blueprint_validation.get("not_applicable") is True
        or blueprint_validation.get("applicable") is False
        or (isinstance(blueprint, dict) and blueprint.get("not_applicable") is True)
        or (isinstance(blueprint, dict) and blueprint.get("applicable") is False)
    )
    blueprint_valid = blueprint_validation.get("blueprint_valid")
    if blueprint_valid is None:
        blueprint_valid = not bool(blueprint_validation.get("invalid"))
    required_hard_gate_applicable = blueprint_validation.get(
        "required_hard_gate_applicable"
    )
    if required_hard_gate_applicable is None:
        required_hard_gate_applicable = blueprint_validation.get(
            "hard_gate_applicable"
        )
    if required_hard_gate_applicable is None:
        required_hard_gate_applicable = blueprint_validation.get(
            "valid_for_hard_gating"
        )
    if required_hard_gate_applicable is None:
        required_hard_gate_applicable = bool(
            (blueprint or {}).get("required") if isinstance(blueprint, dict) else []
        )
    desired_metric_applicable = blueprint_validation.get(
        "desired_metric_applicable"
    )
    if desired_metric_applicable is None:
        desired_metric_applicable = bool(
            not blueprint_not_applicable
            and blueprint_valid
            and isinstance(blueprint, dict)
            and blueprint.get("desired")
        )
    molecule_hard_gate_applicable = blueprint_validation.get(
        "molecule_hard_gate_applicable"
    )
    if molecule_hard_gate_applicable is None:
        molecule_hard_gate_applicable = bool(
            required_hard_gate_applicable or desired_metric_applicable
        )

    return {
        "available": bool(pocket_report),
        "pocket_analysis_report": pocket_report,
        "interaction_blueprint": blueprint,
        "pocket_features": pocket_features,
        "interaction_fingerprint": interaction_fp,
        "consensus_pocket_center": consensus_center,
        "consensus_pocket_residues": consensus_residues,
        # valid_for_hard_gating is nested under interaction_blueprint.validation.
        "valid_for_hard_gating": (
            blueprint_validation.get("valid_for_hard_gating")
        ),
        "blueprint_valid": None if blueprint_not_applicable else bool(blueprint_valid),
        "blueprint_not_applicable": blueprint_not_applicable,
        "required_hard_gate_applicable": bool(required_hard_gate_applicable),
        "desired_metric_applicable": bool(desired_metric_applicable),
        "molecule_hard_gate_applicable": bool(molecule_hard_gate_applicable),
    }


def _artifact_identity_matches(path: Path, expected: object) -> bool:
    if not isinstance(expected, dict) or not path.is_file():
        return False
    expected_digest = expected.get("sha256")
    expected_size = expected.get("size_bytes")
    if not isinstance(expected_digest, str) or not expected_digest:
        return False
    try:
        if int(expected_size) != path.stat().st_size:
            return False
    except (TypeError, ValueError, OSError):
        return False
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                digest.update(chunk)
    except OSError:
        return False
    return digest.hexdigest() == expected_digest


def _count_sdf_record_delimiters(path: Path) -> int | None:
    try:
        with path.open(encoding="utf-8", errors="replace", newline="") as handle:
            return sum(line.rstrip("\r\n") == "$$$$" for line in handle)
    except OSError:
        return None


def _physical_sdf_record_bytes(path: Path) -> list[bytes] | None:
    """Read delimiter-framed records without compressing malformed entries."""
    blocks: list[bytes] = []
    current: list[bytes] = []
    saw_delimiter = False
    try:
        with path.open("rb") as handle:
            for line in handle:
                if line.rstrip(b"\r\n") == b"$$$$":
                    saw_delimiter = True
                    blocks.append(b"".join(current))
                    current = []
                else:
                    current.append(line)
    except OSError:
        return None
    if not saw_delimiter or b"".join(current).strip():
        return None
    return blocks


def _redocking_identity_matches(path: Path, expected: object) -> bool:
    """Verify both the redocking file and its claimed best physical record."""
    if not _artifact_identity_matches(path, expected):
        return False
    if not isinstance(expected, dict):
        return False
    blocks = _physical_sdf_record_bytes(path)
    index = _to_pose_index(expected.get("best_rmsd_pose_index"))
    expected_records = _to_pose_index(expected.get("physical_sdf_records"))
    expected_record_digest = expected.get("best_pose_sdf_record_sha256")
    if (
        blocks is None
        or expected_records != len(blocks)
        or index is None
        or index < 0
        or index >= len(blocks)
        or not isinstance(expected_record_digest, str)
        or not expected_record_digest
    ):
        return False
    return hashlib.sha256(blocks[index]).hexdigest() == expected_record_digest


def _unavailable_redocking_validation(
    redocking: object,
    *,
    reason: str,
) -> dict[str, Any]:
    """Preserve attempt diagnostics but remove every unverified RMSD verdict."""
    sanitized = dict(redocking) if isinstance(redocking, dict) else {}
    reported_reason = sanitized.get("reason")
    if reported_reason not in (None, "") and reported_reason != reason:
        sanitized["reported_reason"] = reported_reason
    sanitized.update({
        "measured": False,
        "applicable": False,
        "not_applicable": True,
        "passed": None,
        "gate": "not_applicable",
        "best_rmsd_angstrom": None,
        "best_rmsd_angstrom_display": None,
        "best_rmsd_pose_index": None,
        "artifact_available": False,
        "reason": reason,
    })
    sanitized.pop("redocked_sdf_artifact", None)
    return sanitized


def _read_step3(project_dir: Path) -> dict[str, Any]:
    step_dir = project_dir / "step3"
    if not step_dir.exists():
        return {"available": False, "reason": "step3 directory missing"}
    docking_report = _safe_json(step_dir / "docking_report.json")
    contract = (
        docking_report.get("docked_pose_artifact_contract")
        if isinstance(docking_report, dict) else None
    )
    if not (
        isinstance(docking_report, dict)
        and docking_report.get("success") is True
        and isinstance(contract, dict)
        and contract.get("status") == "passed"
        and contract.get("committed") is True
    ):
        return {
            "available": False,
            "reason": "docking artifact set is not committed and verified",
            "docking_report": docking_report,
            "docking_scores": [],
            "docked_poses_sdf": None,
            "redocked_reference_sdf": None,
            "redocking_validation": None,
        }

    files = docking_report.get("files") or {}
    docked_path = Path(str(files.get("docked_poses_sdf") or ""))
    scores_path = Path(str(files.get("docking_scores_json") or ""))
    expected_sdf = contract.get("docked_sdf") or {}
    expected_scores = contract.get("docking_scores_json") or {}
    seal_payload: dict[str, Any] = {
        "docked_sdf_sha256": expected_sdf.get("sha256"),
        "docking_scores_sha256": expected_scores.get("sha256"),
        "physical_sdf_records": contract.get("physical_sdf_records"),
        "score_records": contract.get("score_records"),
    }
    try:
        artifact_set_schema = int(
            contract.get("artifact_set_schema_version", 1)
        )
    except (TypeError, ValueError):
        artifact_set_schema = -1
    if artifact_set_schema in {2, 3, 4}:
        expected_redocked = contract.get("redocked_reference_sdf") or {}
        redocked_present = contract.get("redocked_reference_present") is True
        seal_payload.update({
            "redocked_reference_present": redocked_present,
            "redocked_reference_sdf_sha256": (
                expected_redocked.get("sha256") if redocked_present else None
            ),
            "redocked_reference_sdf_size_bytes": (
                expected_redocked.get("size_bytes") if redocked_present else None
            ),
            "redocked_reference_physical_sdf_records": (
                expected_redocked.get("physical_sdf_records")
                if redocked_present else None
            ),
            "redocked_reference_best_rmsd_pose_index": (
                expected_redocked.get("best_rmsd_pose_index")
                if redocked_present else None
            ),
            "redocked_reference_best_pose_sdf_record_sha256": (
                expected_redocked.get("best_pose_sdf_record_sha256")
                if redocked_present else None
            ),
        })
        if artifact_set_schema in {3, 4}:
            sealed_measurement = contract.get("redocking_measurement")
            if not isinstance(sealed_measurement, dict):
                return {
                    "available": False,
                    "reason": "docking artifact set lacks its redocking measurement seal",
                    "docking_report": docking_report,
                    "docking_scores": [],
                    "docked_poses_sdf": None,
                    "redocked_reference_sdf": None,
                    "redocking_validation": None,
                }
            seal_payload["redocking_measurement"] = sealed_measurement
        if artifact_set_schema == 4:
            evaluation_receptor = contract.get("evaluation_receptor")
            protocol_fingerprint = contract.get(
                "score_protocol_fingerprint_sha256"
            )
            reported_fingerprint = (
                (docking_report.get("score_protocol") or {}).get(
                    "fingerprint_sha256"
                )
                if isinstance(docking_report.get("score_protocol"), dict)
                else None
            )
            try:
                protocol_fingerprint_valid = bool(
                    isinstance(protocol_fingerprint, str)
                    and len(protocol_fingerprint) == 64
                    and int(protocol_fingerprint, 16) >= 0
                )
            except ValueError:
                protocol_fingerprint_valid = False
            if not (
                isinstance(evaluation_receptor, dict)
                and protocol_fingerprint_valid
                and protocol_fingerprint == reported_fingerprint
                and _artifact_identity_matches(
                    Path(str(evaluation_receptor.get("path") or "")),
                    evaluation_receptor,
                )
            ):
                return {
                    "available": False,
                    "reason": (
                        "docking artifact set lacks a verified evaluation "
                        "receptor or score-protocol identity"
                    ),
                    "docking_report": docking_report,
                    "docking_scores": [],
                    "docked_poses_sdf": None,
                    "redocked_reference_sdf": None,
                    "redocking_validation": None,
                }
            seal_payload.update({
                "evaluation_receptor_sha256": evaluation_receptor.get(
                    "sha256"
                ),
                "evaluation_receptor_size_bytes": evaluation_receptor.get(
                    "size_bytes"
                ),
                "score_protocol_fingerprint_sha256": protocol_fingerprint,
            })
    elif artifact_set_schema != 1:
        return {
            "available": False,
            "reason": "unsupported docking artifact-set schema",
            "docking_report": docking_report,
            "docking_scores": [],
            "docked_poses_sdf": None,
            "redocked_reference_sdf": None,
            "redocking_validation": None,
        }
    expected_artifact_set_id = "sha256:" + hashlib.sha256(
        json.dumps(
            seal_payload, sort_keys=True, separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    if not (
        _artifact_identity_matches(docked_path, expected_sdf)
        and _artifact_identity_matches(scores_path, expected_scores)
        and contract.get("artifact_set_id") == expected_artifact_set_id
    ):
        return {
            "available": False,
            "reason": "docking artifact hash/size does not match commit marker",
            "docking_report": docking_report,
            "docking_scores": [],
            "docked_poses_sdf": None,
            "redocked_reference_sdf": None,
            "redocking_validation": None,
        }
    docking_scores = _safe_json(scores_path)
    scores_list = (
        docking_scores if isinstance(docking_scores, list)
        else (docking_scores or {}).get("scores")
        if isinstance(docking_scores, dict) else None
    )
    physical_records = _count_sdf_record_delimiters(docked_path)
    if not (
        isinstance(scores_list, list)
        and physical_records == contract.get("physical_sdf_records")
        and len(scores_list) == contract.get("score_records")
    ):
        return {
            "available": False,
            "reason": "docking artifact record counts do not match commit marker",
            "docking_report": docking_report,
            "docking_scores": [],
            "docked_poses_sdf": None,
            "redocked_reference_sdf": None,
            "redocking_validation": None,
        }

    docked_sdf = str(docked_path)
    redocking = (docking_report or {}).get("redocking_validation") if isinstance(docking_report, dict) else None
    if not redocking and isinstance(docking_report, dict):
        # Some earlier runs stored it flat
        redocking = {
            "rmsd": docking_report.get("redocking_rmsd"),
            "threshold": docking_report.get("redocking_rmsd_threshold"),
            "passed": docking_report.get("redocking_passed"),
        }
    redocked_sdf = None
    measurement_seal_valid = True
    redocking_evidence_status = "not_applicable"
    if artifact_set_schema in {3, 4}:
        sealed_measurement = contract.get("redocking_measurement") or {}
        reported_measurement = {
            field: redocking.get(field)
            if isinstance(redocking, dict) else None
            for field in REDOCKING_SEALED_FIELDS
        }
        measurement_seal_valid = (
            set(sealed_measurement) == set(reported_measurement)
            and json.dumps(
                sealed_measurement, sort_keys=True, separators=(",", ":"),
            ) == json.dumps(
                reported_measurement, sort_keys=True, separators=(",", ":"),
            )
        )
        if not measurement_seal_valid:
            redocking = _unavailable_redocking_validation(
                redocking,
                reason=(
                    "redocking measurement fields do not match the committed "
                    "docking report seal"
                ),
            )
            redocking_evidence_status = "invalid"
        else:
            semantic_error = redocking_measurement_error(
                sealed_measurement,
                redocked_reference_present=redocked_present,
            )
            if semantic_error is not None:
                measurement_seal_valid = False
                redocking = _unavailable_redocking_validation(
                    redocking,
                    reason=(
                        "invalid committed redocking measurement: "
                        f"{semantic_error}"
                    ),
                )
                redocking_evidence_status = "invalid"

    if artifact_set_schema in {3, 4} and measurement_seal_valid:
        redocked_present = contract.get("redocked_reference_present") is True
        if redocked_present:
            redocked_path_value = files.get("redocked_reference_sdf")
            contract_identity = contract.get("redocked_reference_sdf")
            validation_identity = (
                redocking.get("redocked_sdf_artifact")
                if isinstance(redocking, dict) else None
            )
            valid_measurement = bool(
                isinstance(redocking, dict)
                and redocking.get("performed") is True
                and redocking.get("measured") is True
                and redocking.get("artifact_available") is True
            )
            try:
                rmsd = float(
                    redocking.get("best_rmsd_angstrom")
                    if isinstance(redocking, dict) else None
                )
                valid_measurement = valid_measurement and math.isfinite(rmsd)
            except (TypeError, ValueError):
                valid_measurement = False
            if redocked_path_value not in (None, ""):
                redocked_path = Path(str(redocked_path_value))
                expected_path = (
                    contract_identity.get("path")
                    if isinstance(contract_identity, dict) else None
                )
                validation_path = (
                    validation_identity.get("path")
                    if isinstance(validation_identity, dict) else None
                )
                try:
                    paths_agree = bool(
                        expected_path not in (None, "")
                        and validation_path not in (None, "")
                        and redocked_path.resolve()
                        == Path(str(expected_path)).resolve()
                        == Path(str(validation_path)).resolve()
                    )
                except OSError:
                    paths_agree = False
                valid_measurement = bool(
                    valid_measurement
                    and paths_agree
                    and contract_identity == validation_identity
                    and _redocking_identity_matches(
                        redocked_path, contract_identity
                    )
                    and _to_pose_index(redocking.get("best_rmsd_pose_index"))
                    == _to_pose_index(
                        contract_identity.get("best_rmsd_pose_index")
                    )
                )
                if valid_measurement:
                    redocked_sdf = str(redocked_path)
            if redocked_sdf is None:
                redocking = _unavailable_redocking_validation(
                    redocking,
                    reason=(
                        "redocking artifact or best-pose identity does not "
                        "match the committed docking report"
                    ),
                )
                redocking_evidence_status = "invalid"
            else:
                redocking_evidence_status = "verified"
        elif isinstance(redocking, dict):
            redocking = _unavailable_redocking_validation(
                redocking,
                reason=(
                    str(redocking.get("reason"))
                    if redocking.get("measured") is not True
                    and redocking.get("reason") not in (None, "")
                    else "docking artifact set does not commit a redocking measurement"
                ),
            )
    elif (
        measurement_seal_valid
        and isinstance(redocking, dict)
        and (
            redocking.get("performed") is True
            or redocking.get("measured") is True
            or redocking.get("passed") is not None
            or redocking.get("best_rmsd_angstrom") is not None
            or redocking.get("rmsd") is not None
        )
    ):
        # Schema 1/2 reports did not commit the measurement fields. Preserve
        # the attempt for audit, but do not expose its RMSD/gate as evidence.
        redocking = _unavailable_redocking_validation(
            redocking,
            reason=(
                "legacy redocking artifact lacks a sealed measurement "
                "identity"
            ),
        )
        redocking_evidence_status = "invalid"

    return {
        "available": True,
        "artifact_set_id": contract.get("artifact_set_id"),
        "docking_report": docking_report,
        "docking_scores": scores_list or [],
        "docked_poses_sdf": docked_sdf,
        "redocked_reference_sdf": redocked_sdf,
        "redocking_validation": redocking,
        "redocking_evidence_status": redocking_evidence_status,
    }


def _read_step4(project_dir: Path) -> dict[str, Any]:
    step_dir = project_dir / "step4"
    if not step_dir.exists():
        return {"available": False, "reason": "step4 directory missing"}
    evaluation_report = _safe_json(step_dir / "evaluation_report.json")
    has_committed_evaluation = bool(
        isinstance(evaluation_report, dict)
        and any(
            isinstance(evaluation_report.get(field), list)
            for field in ("evaluation", "molecules", "results")
        )
    )
    if not has_committed_evaluation:
        return {
            "available": False,
            "reason": "step4 has no committed per-pose evaluation artifact",
            "evaluation_report": evaluation_report,
            "molecules": [],
            "per_molecule": {},
            "per_molecule_entries": {},
            "selected_poses_sdf": None,
            "selection_manifest": [],
        }
    selected_sdf = _existing_or_none(step_dir / "selected_poses.sdf")
    selection_manifest = _safe_csv(step_dir / "selection_manifest.csv")

    molecules: list[dict] = []
    if isinstance(evaluation_report, dict):
        # structure-evaluator writes the per-pose list under "evaluation".
        # Accept "molecules"/"results" too for forward/backward compatibility.
        molecules = (
            evaluation_report.get("evaluation")
            or evaluation_report.get("molecules")
            or evaluation_report.get("results")
            or []
        )

    per_mol = {}
    per_mol_entries: dict[str, list[dict]] = {}
    for m in molecules:
        name = m.get("mol_name") or m.get("name")
        if name:
            name = str(name)
            # Keep the legacy single-entry lookup deterministic, but retain all
            # poses for the identity-aware join used by _build_candidates.
            per_mol.setdefault(name, m)
            per_mol_entries.setdefault(name, []).append(m)

    return {
        "available": bool(evaluation_report or selected_sdf),
        "evaluation_report": evaluation_report,
        "selected_poses_sdf": selected_sdf,
        "selection_manifest": selection_manifest,
        "molecules": molecules,
        "per_molecule": per_mol,
        "per_molecule_entries": per_mol_entries,
    }


def _enforce_step3_step4_lineage(data: dict[str, Any]) -> None:
    """Fail closed when Step 4 was produced from a different docking run."""
    step3 = data.get("step3") or {}
    step4 = data.get("step4") or {}
    if step3.get("available") is not True or step4.get("available") is not True:
        return
    error = evaluation_docking_lineage_error(
        step4.get("evaluation_report"),
        step3.get("docking_report"),
    )
    if error is None:
        step4["docking_lineage_status"] = "verified"
        return
    step4["docking_lineage_status"] = "invalid"
    step4["docking_lineage_error"] = error


def _read_step5(project_dir: Path) -> dict[str, Any]:
    step_dir = project_dir / "step5"
    if not step_dir.exists():
        return {"available": False, "reason": "step5 directory missing"}
    property_report = _safe_json(step_dir / "property_filter_report.json")
    results = (property_report or {}).get("results") or [] if isinstance(property_report, dict) else []
    return {
        "available": bool(property_report),
        "property_filter_report": property_report,
        "results": results,
        "identity_index": _index_records_by_identity(
            results,
            source_cohort_ids=_declared_cohort_ids(property_report),
        ),
    }


def _normalize_synthesis_audit(value: Any) -> dict[str, Any]:
    """Return an explicit AiZynthFinder attempt/status/reason contract.

    Current Step 5 records carry ``attempted`` directly.  Historical reports
    predate that bit, so their attempt state is inferred only for the finite
    set of statuses whose execution semantics are defined.  Unknown statuses
    remain unknown rather than being guessed into either failure or rejection.
    """
    synthesis = dict(value) if isinstance(value, Mapping) else {}
    status = str(synthesis.get("status") or "missing").strip().casefold()
    if not status:
        status = "missing"

    expected_attempted: bool | None = None
    if status in _AIZYNTH_ATTEMPTED_STATUSES:
        expected_attempted = True
    elif status in _AIZYNTH_NOT_ATTEMPTED_STATUSES:
        expected_attempted = False

    recorded_attempted = synthesis.get("attempted")
    if isinstance(recorded_attempted, bool):
        attempted = recorded_attempted
        attempted_source = "explicit"
    elif expected_attempted is not None:
        attempted = expected_attempted
        attempted_source = "inferred_from_legacy_status"
    else:
        attempted = None
        attempted_source = "unknown"

    contract_consistent = not (
        isinstance(recorded_attempted, bool)
        and expected_attempted is not None
        and recorded_attempted is not expected_attempted
    )
    reason = synthesis.get("reason") or synthesis.get("error")
    if reason in (None, ""):
        reason = _AIZYNTH_DEFAULT_REASONS.get(
            status, "unrecognized_synthesis_status"
        )

    synthesis.update({
        "status": status,
        "attempted": attempted,
        "attempted_source": attempted_source,
        "reason": str(reason),
        "contract_consistent": contract_consistent,
        "interpretation": _AIZYNTH_STATUS_MEANINGS.get(
            status,
            "Unrecognized status; do not infer an attempt, failure, or "
            "molecule verdict.",
        ),
    })
    if not contract_consistent:
        synthesis["contract_warning"] = (
            f"attempted={recorded_attempted!r} conflicts with status={status!r}"
        )
    return synthesis


def _build_synthesis_summary(step5: Any) -> dict[str, Any]:
    """Summarize execution states independently from molecule verdict counts."""
    if not isinstance(step5, Mapping) or not step5.get("available"):
        reason = (
            step5.get("reason")
            if isinstance(step5, Mapping)
            else None
        )
        return {
            "available": False,
            "reason": reason or "Step 5 property report unavailable",
        }

    results = [
        row for row in (step5.get("results") or [])
        if isinstance(row, Mapping)
    ]
    audits = [
        _normalize_synthesis_audit(row.get("synthesis")) for row in results
    ]
    status_counts: dict[str, int] = {}
    attempted_by_status: dict[str, set[bool | None]] = {}
    for audit in audits:
        status = str(audit["status"])
        status_counts[status] = status_counts.get(status, 0) + 1
        attempted_by_status.setdefault(status, set()).add(audit["attempted"])

    ordered_statuses = [
        status for status in _AIZYNTH_STATUS_ORDER if status in status_counts
    ]
    ordered_statuses.extend(sorted(set(status_counts) - set(ordered_statuses)))
    status_rows = []
    for status in ordered_statuses:
        observed = attempted_by_status[status]
        attempted: bool | str | None
        if len(observed) == 1:
            attempted = next(iter(observed))
        else:
            attempted = "mixed"
        status_rows.append({
            "status": status,
            "count": status_counts[status],
            "attempted": attempted,
            "meaning": _AIZYNTH_STATUS_MEANINGS.get(
                status,
                "Unrecognized status; no execution or verdict inference is safe.",
            ),
        })

    return {
        "available": True,
        "evaluated_records": len(results),
        "attempted": sum(audit["attempted"] is True for audit in audits),
        "not_attempted": sum(audit["attempted"] is False for audit in audits),
        "attempt_unknown": sum(audit["attempted"] is None for audit in audits),
        "operational_errors": status_counts.get("error", 0),
        "contract_inconsistent": sum(
            audit["contract_consistent"] is False for audit in audits
        ),
        "status_counts": status_counts,
        "status_rows": status_rows,
        "verdict_independence": (
            "AiZynthFinder execution status is not a molecule classification. "
            "Not-attempted and error rows must not be counted as automatic "
            "molecule rejects."
        ),
    }


def _read_step6(project_dir: Path) -> dict[str, Any]:
    step_dir = project_dir / "step6"
    if not step_dir.exists():
        return {"available": False, "reason": "step6 directory missing"}
    coordinator_report = _safe_json(step_dir / "coordinator_report.json")
    round_summary = _safe_json(step_dir / "round_summary.json")
    round_history = _safe_json(step_dir / "round_history.json")
    candidate_ranking_json = _safe_json(step_dir / "candidate_ranking.json")

    step_result = _safe_json(step_dir / "result.json")
    success = bool((step_result or {}).get("success")) if isinstance(step_result, dict) else bool(coordinator_report)

    decision = None
    decision_outcome = None
    goal_met = None
    valid_for_nomination = None
    quality_gate_passed = None
    diagnosis = None
    round_assessment = None
    assessment_schema_version = None
    if isinstance(coordinator_report, dict):
        decision = coordinator_report.get("decision")
        embedded_summary = coordinator_report.get("round_summary") or round_summary or {}
        decision_outcome = coordinator_report.get(
            "decision_outcome", embedded_summary.get("decision_outcome")
        )
        goal_met = coordinator_report.get("goal_met", embedded_summary.get("goal_met"))
        explicit_nomination_validity = coordinator_report.get(
            "valid_for_nomination", embedded_summary.get("valid_for_nomination")
        )
        if isinstance(explicit_nomination_validity, bool):
            valid_for_nomination = explicit_nomination_validity
        elif decision == "terminate_invalid":
            # Legacy coordinator reports predate the explicit field, but their
            # hard-invalid decision must never be a final candidate source.
            valid_for_nomination = False
        # quality_gate_passed lives inside round_summary, not at the top level.
        quality_gate_passed = (
            coordinator_report.get("quality_gate_passed")
            if coordinator_report.get("quality_gate_passed") is not None
            else embedded_summary.get("quality_gate_passed")
        )
        diagnosis = ((coordinator_report.get("next_round_instructions") or {}).get("diagnosis")
                     if coordinator_report.get("next_round_instructions") else None)
        assessment_schema_version = coordinator_report.get("assessment_schema_version")
        if assessment_schema_version == "round-assessment/v1":
            round_assessment = coordinator_report.get("round_assessment") or {}
            valid_for_nomination = (
                round_assessment.get("valid_for_nomination") is True
                and round_assessment.get("status") in {"valid", "review_required"}
            )
            decision_outcome = None
            goal_met = None
            diagnosis = (coordinator_report.get("diagnostics") or {}).get("diagnosis")

    return {
        "available": bool(coordinator_report),
        "success": success,
        "coordinator_report": coordinator_report,
        "round_summary": round_summary,
        "round_history": round_history,
        "candidate_ranking_json": candidate_ranking_json,
        "decision": decision,
        "decision_outcome": decision_outcome,
        "goal_met": goal_met,
        "valid_for_nomination": valid_for_nomination,
        "quality_gate_passed": quality_gate_passed,
        "diagnosis": diagnosis,
        "assessment_schema_version": assessment_schema_version,
        "round_assessment": round_assessment,
    }


# ---------- derived: header / candidates / gates / screening / trend ----------

def _build_header(data: dict, override_round_label: str | None = None) -> dict:
    step0 = data.get("step0") or {}
    step1 = data.get("step1") or {}
    step6 = data.get("step6") or {}
    coord = step6.get("coordinator_report") or {}
    metadata = step0.get("metadata") or {}
    # Prefer the round_summary embedded in coordinator_report.json (always
    # complete) over the standalone round_summary.json side file.
    round_summary = coord.get("round_summary") or step6.get("round_summary") or {}

    return {
        "title": _guess_title(metadata, coord),
        "target_name": metadata.get("target_name") or coord.get("target_name"),
        "pdb_id": metadata.get("pdb_id") or coord.get("pdb_id"),
        "chembl_target_id": metadata.get("chembl_target_id") or coord.get("chembl_target_id"),
        "round_label": override_round_label or (f"round{data.get('round_number')}" if data.get("round_number") else "round"),
        "assessment_schema_version": step6.get("assessment_schema_version"),
        "round_assessment": step6.get("round_assessment"),
        "assessment_status": (step6.get("round_assessment") or {}).get("status"),
        "decision": step6.get("decision"),
        "decision_outcome": step6.get("decision_outcome"),
        "goal_met": step6.get("goal_met"),
        "valid_for_nomination": step6.get("valid_for_nomination"),
        "decision_rationale": coord.get("decision_rationale"),
        "quality_gate_passed": step6.get("quality_gate_passed"),
        "primary_bottleneck": (
            round_summary.get("primary_bottleneck")
            or ((step6.get("diagnosis") or {}).get("bottleneck_type") if step6.get("diagnosis") else None)
        ),
        # round_summary keys are total_evaluated / candidate_total / keep / modify / reject.
        "num_screened": round_summary.get("total_evaluated") or round_summary.get("total_molecules") or round_summary.get("num_total"),
        "candidate_pool_size": round_summary.get("candidate_total") or round_summary.get("candidate_pool_size"),
        "num_keep": round_summary.get("keep") if round_summary.get("keep") is not None else round_summary.get("num_keep"),
        "num_modify": round_summary.get("modify") if round_summary.get("modify") is not None else round_summary.get("num_modify"),
        "num_reject": round_summary.get("reject") if round_summary.get("reject") is not None else round_summary.get("num_reject"),
        "candidate_intrinsic_keep": round_summary.get("candidate_intrinsic_keep"),
        "candidate_quality_gate_blocked_keep": round_summary.get("candidate_quality_gate_blocked_keep"),
    }


def _build_candidates(
    data: dict,
    top_n: int,
    *,
    nomination_only: bool = False,
) -> list[dict]:
    """Top-N per cohort from Step 6 (DEC_ already excluded upstream).

    Enrich each row with step3 scores, step4 evaluation diagnostics, and step5
    properties.  OverallRank is cohort-local, so rows are partitioned before
    sorting/truncation and are never compared across cohorts.
    """
    if top_n <= 0:
        return []

    step6 = data.get("step6") or {}
    # candidate_ranking.json is authoritative for every numeric field, so it also
    # drives the row loop; the retired CSV only ever duplicated it.
    ranking_json = step6.get("candidate_ranking_json") or []
    ranking_rows = ranking_json
    step3 = data.get("step3") or {}
    step4 = data.get("step4") or {}
    step5 = data.get("step5") or {}
    step4_lineage_error = step4.get("docking_lineage_error")

    step3_index = _index_records_by_identity(
        step3.get("docking_scores") or [],
        source_cohort_ids=_declared_cohort_ids(step3.get("docking_report")),
    )
    step4_index = _index_records_by_identity(
        step4.get("molecules")
        or [
            item
            for values in (step4.get("per_molecule_entries") or {}).values()
            for item in values
        ],
        source_cohort_ids=_declared_cohort_ids(step4.get("evaluation_report")),
    )
    step5_index = step5.get("identity_index") or _index_records_by_identity(
        step5.get("results") or [],
        source_cohort_ids=_declared_cohort_ids(step5.get("property_filter_report")),
    )

    # candidate_ranking.json is authoritative for numeric fields when present.
    ranking_json_entries: list[dict] = []
    if isinstance(ranking_json, list):
        for entry in ranking_json:
            if isinstance(entry, dict):
                ranking_json_entries.append(entry)
    elif isinstance(ranking_json, dict):
        for entry in ranking_json.get("candidates", []) or []:
            if isinstance(entry, dict):
                ranking_json_entries.append(entry)
    ranking_json_index = _index_records_by_identity(
        ranking_json_entries,
        source_cohort_ids=_declared_cohort_ids(step6.get("coordinator_report")),
    )

    candidates: list[dict] = []
    for row in ranking_rows:
        mol_name = row.get("mol_name") or row.get("name") or row.get("compound_id")
        if not mol_name:
            continue
        if str(mol_name).startswith("DEC_"):
            # Defensive: coordinator should have excluded these; skip if it didn't.
            continue

        mol_name = str(mol_name)
        row_pose_index = _first_non_none(
            _to_pose_index(row.get("pose_index")),
            _to_pose_index(row.get("selected_pose_index")),
        )
        pose_id = row.get("pose_id")
        row_cohort_id = row.get("cohort_id")
        json_entry = _lookup_identity(
            ranking_json_index, mol_name, row_pose_index, pose_id=pose_id,
            cohort_id=row_cohort_id,
            cohort_identity=row,
        )
        candidate_cohort_identity = (
            json_entry if _cohort_tuple_valid(json_entry) else row
        )
        cohort_id = str(
            json_entry.get("cohort_id")
            or row_cohort_id
            or "cohort:unspecified"
        )
        classification = (
            json_entry.get("classification")
            or row.get("classification")
            # ``classification`` is the published spelling of the level-2 fold
            # verdict; ``metric_gate`` is the same value under its machine name
            #.  Older artifacts carry only the latter.
            or json_entry.get("metric_gate")
            or row.get("metric_gate")
        )
        nomination_status = (
            json_entry.get("nomination_status")
            or row.get("nomination_status")
            or "eligible"
        )
        if nomination_only:
            # Nomination requires an affirmative ``keep`` (the executable policy carries
            # only keep molecules forward).  A blank or absent verdict is an
            # evidence gap, not a keep: the row stays in the full report but may
            # not be presented as a nominated candidate.
            if str(classification or "").strip().casefold() != "keep":
                continue
            # Only a real nomination block excludes a row.  A coverage
            # annotation such as ``coverage_review`` says the evidence is thin
            # (the executable policy publishes that rank as provisional), not that the
            # molecule is disqualified.
            if str(nomination_status).startswith("blocked"):
                continue
        pose_id = pose_id or json_entry.get("pose_id")
        json_pose_index = _first_non_none(
            _to_pose_index(json_entry.get("pose_index")),
            _to_pose_index(json_entry.get("selected_pose_index")),
        )
        pose_index = _first_non_none(
            json_pose_index, row_pose_index
        ) if pose_id else _first_non_none(row_pose_index, json_pose_index)
        eval_entry = _lookup_identity(
            step4_index,
            mol_name,
            pose_index,
            pose_id=pose_id,
            cohort_id=cohort_id,
            cohort_identity=candidate_cohort_identity,
        )
        if step4_lineage_error:
            # Keep the raw report under step4 for audit, but do not project its
            # pose QC (or Step 6 copies of it) into a candidate assembled with
            # another Step 3 artifact set.
            eval_entry = {}
        pose_id = pose_id or eval_entry.get("pose_id")
        eval_pose_index = _first_non_none(
            _to_pose_index(eval_entry.get("pose_index")),
            _to_pose_index(eval_entry.get("selected_pose_index")),
        )
        pose_index = _first_non_none(
            eval_pose_index, pose_index
        ) if pose_id else _first_non_none(pose_index, eval_pose_index)
        source_pose_index = _first_non_none(
            _to_pose_index(row.get("source_pose_index")),
            _to_pose_index(json_entry.get("source_pose_index")),
            _to_pose_index(eval_entry.get("source_pose_index")),
        )
        chemical_identity_id = (
            row.get("chemical_identity_id")
            or json_entry.get("chemical_identity_id")
            or (
                (eval_entry.get("chemical_identity") or {}).get("identity_id")
            )
            or eval_entry.get("chemical_identity_id")
        )
        sdf_molblock_witness = (
            row.get("sdf_molblock_sha256")
            or json_entry.get("sdf_molblock_sha256")
            or eval_entry.get("sdf_molblock_sha256")
        )
        prop_entry = _lookup_identity(
            step5_index,
            mol_name,
            pose_index,
            pose_id=pose_id,
            chemical_identity_id=chemical_identity_id,
            source_pose_index=source_pose_index,
            cohort_id=cohort_id,
            cohort_identity=candidate_cohort_identity,
        )

        docking_scores = _lookup_identity(
            step3_index,
            mol_name,
            pose_index,
            pose_id=pose_id,
            sdf_molblock_sha256=sdf_molblock_witness,
            cohort_id=cohort_id,
            cohort_identity=candidate_cohort_identity,
        )
        docking_lineage_error = (
            str(step4_lineage_error) if step4_lineage_error else None
        )
        expected_docking_set_id = step3.get("artifact_set_id")
        if (
            docking_lineage_error is None
            and expected_docking_set_id not in (None, "")
        ):
            candidate_docking_set_id = (
                row.get("docking_artifact_set_id")
                or json_entry.get("docking_artifact_set_id")
            )
            if candidate_docking_set_id in (None, ""):
                docking_lineage_error = (
                    "Step 6 candidate lacks verified Step 3 artifact-set "
                    "lineage"
                )
            elif str(candidate_docking_set_id) != str(
                expected_docking_set_id
            ):
                docking_lineage_error = (
                    "Step 6 candidate belongs to a different Step 3 "
                    "docking artifact set"
                )
        if docking_lineage_error is not None:
            docking_scores = {}
        # Step 4 and Step 6 historically embedded copies of Step 3 scores.
        # Those copies are diagnostics only: once the committed Step 3 artifact
        # is unavailable or fails this identity join, they must not become an
        # alternate source of docking evidence.  Non-docking rescoring/QC
        # metrics remain usable from their owning downstream stage.
        step3_pose_scores = collect_pose_scores(docking_scores)
        downstream_pose_scores = (
            [] if step4_lineage_error else [
                score
                for score in collect_pose_scores(
                    json_entry,
                    eval_entry,
                    {
                        **row,
                        "mol_name": mol_name,
                        "pose_index": pose_index,
                    },
                )
                if not _is_docking_pose_score(score)
            ]
        )
        pose_score_contract = collect_pose_scores({
            "pose_scores": [*step3_pose_scores, *downstream_pose_scores],
        })
        docking_evidence_status = (
            "identity_verified" if docking_scores
            else "invalid" if docking_lineage_error
            else "unavailable"
        )
        docking_evidence_reason = None
        if not docking_scores:
            docking_evidence_reason = (
                docking_lineage_error
                or
                step3.get("reason")
                or "no identity-matched record in committed Step 3 scores"
            )
            if not str(nomination_status).startswith("blocked"):
                nomination_status = "blocked_docking_evidence"
            if nomination_only:
                continue
        primary_metric_id = (
            json_entry.get("primary_score_metric_id")
            or eval_entry.get("primary_score_metric_id")
            or row.get("primary_score_metric_id")
        )
        if not primary_metric_id:
            primary_metric_id = resolve_primary_metric(
                [pose_score_contract],
                preferred_metric_ids=CANDIDATE_RANKING_PREFERENCE,
            )
        primary_score = find_score(pose_score_contract, primary_metric_id)

        smiles = (
            row.get("smiles")
            or json_entry.get("smiles")
            or eval_entry.get("smiles")
            or prop_entry.get("smiles")
        )
        declared_schema_version = _first_non_none(
            _to_int(json_entry.get("scoring_schema_version")),
            _to_int(row.get("scoring_schema_version")),
        )
        declared_ranking_schema = str(
            json_entry.get("ranking_schema")
            or row.get("ranking_schema")
            or ""
        )
        explicitly_legacy = (
            declared_schema_version == 1
            or declared_ranking_schema.startswith("legacy")
        )
        modern_scoring = (
            (
                declared_schema_version is not None
                and declared_schema_version >= 2
            )
            or (
                not explicitly_legacy
                and (
                    "overall_rank" in json_entry
                    or "overall_score" in json_entry
                    or "overall_rank" in row
                    or "overall_score" in row
                )
            )
        )
        overall_rank = _first_non_none(
            _to_rank(json_entry.get("overall_rank")),
            _to_rank(row.get("overall_rank")),
        )
        # Schema-v2 explicit nulls are authoritative.  Never backfill one from
        # rank/weighted_score_rank written by the former flat-score contract.
        rank = overall_rank if modern_scoring else _first_non_none(
            _to_rank(json_entry.get("rank")),
            _to_rank(row.get("rank")),
        )
        normalization_pool = str(
            json_entry.get("normalization_pool")
            or row.get("normalization_pool")
            or "nomination_candidates"
        )
        ranking_schema = (
            "hierarchical_v2"
            if modern_scoring
            else (declared_ranking_schema or "legacy")
        )
        step_scores = json_entry.get("step_scores") or {}

        candidates.append({
            "rank": rank,
            "overall_rank": overall_rank,
            "overall_score": _first_float(
                json_entry.get("overall_score"),
                row.get("overall_score"),
            ),
            "rank_status": (
                json_entry.get("rank_status")
                or row.get("rank_status")
                or row.get("coverage_status")
            ),
            "metric_gate": (
                json_entry.get("metric_gate") or row.get("metric_gate")
            ),
            "gate_status": json_entry.get("gate_status"),
            "coverage_status": (
                json_entry.get("coverage_status")
                or row.get("coverage_status")
            ),
            "metric_assessments": json_entry.get("metric_assessments") or [],
            "planned_metric_coverage_incomplete": json_entry.get(
                "planned_metric_coverage_incomplete"
            ),
            "planned_metric_evidence_incomplete": json_entry.get(
                "planned_metric_evidence_incomplete"
            ),
            "planned_but_failed_metric_ids": json_entry.get(
                "planned_but_failed_metric_ids"
            ),
            "scoring_schema_version": (
                (
                    declared_schema_version or SCORING_SCHEMA_VERSION
                    if modern_scoring else declared_schema_version
                )
            ),
            "ranking_schema": ranking_schema,
            "cohort_id": cohort_id,
            "protein_state_id": (
                json_entry.get("protein_state_id")
                or row.get("protein_state_id")
            ),
            "pocket_id": (
                json_entry.get("pocket_id") or row.get("pocket_id")
            ),
            "round_id": (
                json_entry.get("round_id") or row.get("round_id")
            ),
            "protocol_profile_id": (
                json_entry.get("protocol_profile_id")
                or row.get("protocol_profile_id")
            ),
            "metric_catalog_version": (
                json_entry.get("metric_catalog_version")
                or row.get("metric_catalog_version")
            ),
            "cohort_identity_complete": _first_non_none(
                json_entry.get("cohort_identity_complete"),
                row.get("cohort_identity_complete"),
            ),
            "cohort_identity_missing_fields": (
                json_entry.get("cohort_identity_missing_fields")
                if json_entry.get("cohort_identity_missing_fields") is not None
                else row.get("cohort_identity_missing_fields")
            ),
            "cohort_identity_derived_fields": (
                json_entry.get("cohort_identity_derived_fields")
                if json_entry.get("cohort_identity_derived_fields") is not None
                else row.get("cohort_identity_derived_fields")
            ),
            "cohort_provenance_sources": (
                json_entry.get("cohort_provenance_sources")
                or row.get("cohort_provenance_sources")
                or {}
            ),
            "normalization_pool": normalization_pool,
            "rank_partition_id": (
                json_entry.get("rank_partition_id")
                or row.get("rank_partition_id")
                or f"{cohort_id}|normalization_pool={normalization_pool}"
            ),
            "rank_scope": (
                json_entry.get("rank_scope")
                or row.get("rank_scope")
                or "within_cohort_and_normalization_pool"
            ),
            "category_weighted_coverage": _first_float(
                json_entry.get("category_weighted_coverage"),
                row.get("category_weighted_coverage"),
            ),
            "category_scores": json_entry.get("category_scores") or {},
            "step_scores": step_scores,
            **{
                f"{step_id}_{field}": _first_float(
                    json_entry.get(f"{step_id}_{field}"),
                    row.get(f"{step_id}_{field}"),
                    (step_scores.get(step_id) or {}).get(
                        "rank_avg" if field == "rank" else "score"
                    ),
                )
                for step_id in ("step4", "step5", "step6")
                for field in ("score", "rank")
            },
            "mol_name": str(mol_name),
            "smiles": smiles,
            "pose_index": pose_index,
            "source_pose_index": source_pose_index,
            "pose_id": pose_id,
            "chemical_identity_id": chemical_identity_id,
            "docking_artifact_set_id": (
                expected_docking_set_id
                if docking_evidence_status == "identity_verified"
                else None
            ),
            "selection_rank": _to_int(eval_entry.get("selection_rank")),
            "classification": classification,
            "intrinsic_classification": (
                json_entry.get("intrinsic_classification")
                or row.get("intrinsic_classification")
                or json_entry.get("classification")
                or row.get("classification")
            ),
            "nomination_status": nomination_status,
            "docking_evidence_status": docking_evidence_status,
            "docking_evidence_reason": docking_evidence_reason,
            "verdict": (
                None if step4_lineage_error
                else eval_entry.get("verdict")
                or row.get("eval_verdict")
                or row.get("verdict")
            ),
            "flags": (
                [] if step4_lineage_error
                else eval_entry.get("flags") or row.get("flags") or []
            ),
            "diagnosis": (
                [] if step4_lineage_error
                else eval_entry.get("diagnosis") or []
            ),
            "modification_hint": (
                None if step4_lineage_error
                else eval_entry.get("modification_hint")
                or row.get("modification_hint")
            ),
            "recommendation": (
                None if step4_lineage_error
                else eval_entry.get("recommendation")
                or row.get("recommendation")
            ),
            "docking_scores": {
                "CNNscore": _first_float(
                    docking_scores.get("CNNscore"),
                ),
                "CNNaffinity": _first_float(
                    docking_scores.get("CNNaffinity"),
                ),
                "minimizedAffinity": _first_float(
                    docking_scores.get("minimizedAffinity"),
                ),
            },
            "pose_scores": pose_score_contract,
            "primary_score_metric_id": primary_metric_id,
            "primary_score": primary_score,
            "blueprint_match": eval_entry.get("blueprint_match") or {},
            "posebusters": eval_entry.get("posebusters") or {},
            "posecheck": eval_entry.get("posecheck") or {},
            "pocket_contact_fraction": eval_entry.get("pocket_contact_fraction") or {},
            "properties": (prop_entry.get("properties") if isinstance(prop_entry, dict) else None) or {},
            "lipinski": prop_entry.get("lipinski") or {},
            "pains": prop_entry.get("pains") or {},
            "admet": prop_entry.get("admet"),
            "chemfh": prop_entry.get("chemfh"),
            "synthesis": _normalize_synthesis_audit(
                prop_entry.get("synthesis")
            ),
            "nomination_eligibility": prop_entry.get("nomination_eligibility"),
            "property_verdict": prop_entry.get("verdict"),
            "property_flags": prop_entry.get("flags") or [],
        })

    return _partition_and_truncate_candidates(candidates, top_n=top_n)


def _candidate_ranking_contract(top_n: int) -> dict:
    return {
        "schema_version": 2,
        "rank_source": "overall_rank",
        "rank_type": "float_or_null",
        "partition_fields": ["cohort_id", "normalization_pool"],
        "cross_cohort_comparison": False,
        "top_n_scope": "per_partition",
        "top_n_per_partition": top_n,
        "mixed_schema_policy": "legacy_rows_unranked_when_schema2_present",
    }


def _partition_and_truncate_candidates(
    candidates: list[dict],
    *,
    top_n: int,
) -> list[dict]:
    """Sort and truncate within independent cohort/pool partitions only."""
    has_schema2 = any(
        candidate.get("ranking_schema") == "hierarchical_v2"
        for candidate in candidates
    )
    has_legacy = any(
        candidate.get("ranking_schema") != "hierarchical_v2"
        for candidate in candidates
    )
    if has_schema2 and has_legacy:
        for candidate in candidates:
            if candidate.get("ranking_schema") == "hierarchical_v2":
                continue
            candidate["rank"] = None
            candidate["overall_rank"] = None
            candidate["ranking_schema"] = "legacy_unranked_in_mixed_schema"
            candidate["rank_status"] = "legacy_schema_unranked_in_mixed_output"

    partitions: dict[tuple[str, str], list[dict]] = {}
    for candidate in candidates:
        partition_key = (
            str(candidate.get("cohort_id") or "cohort:unspecified"),
            str(candidate.get("normalization_pool") or "nomination_candidates"),
        )
        partitions.setdefault(partition_key, []).append(candidate)

    selected: list[dict] = []
    for partition_key in sorted(partitions):
        partition = sorted(
            partitions[partition_key],
            key=_candidate_partition_sort_key,
        )
        if not has_schema2:
            # Legacy artifacts sometimes omitted a rank.  Preserve their old
            # presentation fallback, but scope it to this cohort rather than
            # fabricating one global sequence across unrelated cohorts.
            for position, candidate in enumerate(partition, start=1):
                if candidate.get("rank") is None:
                    candidate["rank"] = position
        selected.extend(partition[:top_n])
    return selected


def _candidate_partition_sort_key(candidate: dict) -> tuple:
    rank = _to_rank(candidate.get("rank"))
    schema_order = (
        0 if candidate.get("ranking_schema") == "hierarchical_v2" else 1
    )
    pose_index = _to_int(candidate.get("pose_index"))
    return (
        schema_order,
        rank is None,
        rank if rank is not None else float("inf"),
        str(candidate.get("mol_name") or ""),
        pose_index if pose_index is not None else -1,
    )


def _build_quality_gates(data: dict) -> list[dict]:
    step1 = data.get("step1") or {}
    step3 = data.get("step3") or {}
    step4 = data.get("step4") or {}
    step6 = data.get("step6") or {}

    pocket_report = step1.get("pocket_analysis_report") or {}
    pocket_validation = pocket_report.get("validation") or {}
    blueprint = step1.get("interaction_blueprint") or {}
    # valid_for_hard_gating is nested under interaction_blueprint.validation.
    blueprint_validation = blueprint.get("validation") if isinstance(blueprint, dict) else {}
    blueprint_validation = blueprint_validation or {}
    blueprint_not_applicable = bool(step1.get("blueprint_not_applicable"))
    blueprint_valid = step1.get("blueprint_valid")
    if blueprint_valid is None and not blueprint_not_applicable:
        blueprint_valid = blueprint_validation.get("blueprint_valid")
    if blueprint_valid is None and not blueprint_not_applicable:
        blueprint_valid = not bool(blueprint_validation.get("invalid"))
    redocking = step3.get("redocking_validation") or {}
    redocking_evidence_status = step3.get("redocking_evidence_status")
    coord = step6.get("coordinator_report") or {}
    coord_gates = coord.get("quality_gates") or {}
    screening_metrics = coord.get("screening_validation") or {}

    # Redocking: a check that produced no RMSD is not a failure — surface it as
    # not-applicable rather than a red gate.  That covers validate_redocking=False
    # and a missing reference ligand (``performed: False``) as well as a run that
    # started and then crashed, timed out, or found no atom-compatible pose: the
    # engine reports those with the tri-state ``not_applicable``/``applicable``
    # flags while ``performed`` stays True, because the attempt itself is worth
    # auditing (docking/scripts/gnina_utils.py ``unmeasured``).
    redocking_invalid = bool(
        redocking_evidence_status == "invalid"
        or step3.get("available") is False
    )
    if redocking_invalid and redocking_evidence_status in (None, ""):
        redocking_evidence_status = "invalid"
    redocking_not_applicable = bool(
        not redocking_invalid
        and (
            redocking.get("performed") is False
            or redocking.get("not_applicable") is True
            or redocking.get("applicable") is False
            or str(redocking.get("gate") or "").casefold()
            == "not_applicable"
        )
    )
    redocking_passed = None if redocking_not_applicable else redocking.get("passed")

    gates = [
        {
            "name": "Redocking RMSD",
            "passed": redocking_passed,
            "not_applicable": redocking_not_applicable,
            "evidence_status": redocking_evidence_status,
            "protocol_blocker": redocking_invalid,
            "reason": redocking.get("reason") or step3.get("reason"),
            "value": redocking.get("best_rmsd_angstrom", redocking.get("rmsd")),
            "threshold": redocking.get("rmsd_threshold_angstrom", redocking.get("threshold")),
            "unit": "Å",
            "source": "step3/docking_report.json",
        },
        {
            "name": "Pocket validation",
            "passed": pocket_validation.get("passed"),
            "value": pocket_validation.get("known_site_recall", pocket_validation.get("recall")),
            "threshold": pocket_validation.get("recall_threshold", 0.8),
            "unit": "recall",
            "source": "step1/pocket_analysis_report.json",
        },
        {
            "name": "Blueprint validity",
            "passed": None if blueprint_not_applicable else bool(blueprint_valid),
            "not_applicable": blueprint_not_applicable,
            "value": None if blueprint_not_applicable else bool(blueprint_valid),
            "threshold": True,
            "unit": "boolean",
            "required_constraints": (
                len((blueprint or {}).get("required") or [])
                if isinstance(blueprint, dict) else None
            ),
            "required_hard_gate_applicable": bool(
                step1.get("required_hard_gate_applicable")
            ),
            "desired_metric_applicable": bool(
                step1.get("desired_metric_applicable")
            ),
            "molecule_hard_gate_applicable": bool(
                step1.get("molecule_hard_gate_applicable")
            ),
            "source": "step1/interaction_blueprint.json",
        },
        {
            "name": "Screening AUROC",
            "passed": (coord_gates.get("screening_benchmark") or {}).get("passed"),
            "value": screening_metrics.get("AUROC") if screening_metrics.get("AUROC") is not None else screening_metrics.get("auroc"),
            "threshold": ((coord_gates.get("screening_benchmark") or {}).get("criteria") or {}).get("min_AUROC"),
            "unit": "AUROC",
            "source": "step6/coordinator_report.json",
        },
        {
            "name": "Coordinator quality gate",
            "passed": step6.get("quality_gate_passed"),
            "value": coord.get("quality_gate_status"),
            "threshold": None,
            "unit": None,
            "source": "step6/coordinator_report.json",
        },
    ]
    if step4.get("docking_lineage_status") in {"verified", "invalid"}:
        lineage_verified = step4["docking_lineage_status"] == "verified"
        gates.insert(0, {
            "name": "Step 3 → Step 4 artifact lineage",
            "passed": lineage_verified,
            "not_applicable": False,
            "evidence_status": step4["docking_lineage_status"],
            "protocol_blocker": not lineage_verified,
            "reason": step4.get("docking_lineage_error"),
            "value": (
                step3.get("artifact_set_id") if lineage_verified else None
            ),
            "threshold": "exact artifact-set/receptor/protocol match",
            "unit": None,
            "source": "step3/docking_report.json + step4/evaluation_report.json",
        })
    return gates


def _build_screening(data: dict) -> dict:
    coord = (data.get("step6") or {}).get("coordinator_report") or {}
    # The coordinator writes the retrospective benchmark under
    # "screening_validation" (ranking_engine.screening_metrics output). Accept a
    # couple of legacy aliases. Subkeys are AUROC / AUROC_bootstrap_95ci /
    # EF1%/EF5%/EF10% / num_actives / num_decoys — mind the exact casing.
    bench = (
        coord.get("screening_validation")
        or coord.get("screening_benchmark")
        or coord.get("retrospective_screening")
        or {}
    )
    if not isinstance(bench, dict) or not bench.get("available"):
        return {"available": False, "reason": "no available screening_validation in coordinator_report.json"}
    average_precision = bench.get("average_precision")
    if average_precision is None and bench.get("PR_AUC_method") == "tie_aware_average_precision":
        # A generic PR-AUC may mean trapezoidal integration, so only use the
        # older rounded alias when its recorded method explicitly identifies AP.
        average_precision = bench.get("PR_AUC")
    return {
        "available": True,
        "auroc": bench.get("AUROC") if bench.get("AUROC") is not None else bench.get("auroc"),
        "auroc_ci": bench.get("AUROC_bootstrap_95ci") or bench.get("auroc_ci") or bench.get("bootstrap_ci"),
        "average_precision": average_precision,
        "average_precision_method": bench.get("average_precision_method") or bench.get("PR_AUC_method"),
        "PR_AUC": bench.get("PR_AUC"),
        "PR_AUC_method": bench.get("PR_AUC_method"),
        "PR_AUC_bootstrap_95ci": bench.get("PR_AUC_bootstrap_95ci"),
        **{
            key: bench.get(key)
            for key in (
                "top20_size", "top20_actives", "top20_hit_ratio",
                "top20_observed_tie_policy", "top20_boundary_tie",
                "top20_boundary_tie_size", "top20_boundary_slots",
                "top20_expected_actives", "top20_expected_hit_ratio",
                "top20_expectation_method", "top20_active_count_range",
                "top20_hit_ratio_range",
            )
        },
        "ef": {
            "1%": bench.get("EF1%", bench.get("ef_1", bench.get("EF1"))),
            "5%": bench.get("EF5%", bench.get("ef_5", bench.get("EF5"))),
            "10%": bench.get("EF10%", bench.get("ef_10", bench.get("EF10"))),
        },
        "num_actives": bench.get("num_actives"),
        "num_decoys": bench.get("num_decoys"),
        # The screening pass/fail lives in the quality_gates.screening_benchmark
        # entry, not in the metrics dict itself.
        "passed": ((coord.get("quality_gates") or {}).get("screening_benchmark") or {}).get("passed"),
        "roc_curve": bench.get("roc_curve"),  # {"fpr": [...], "tpr": [...]} if present
    }


def _build_trend(per_round: list[dict]) -> list[dict]:
    rows = []
    for r in per_round:
        step6 = r.get("step6") or {}
        summary = step6.get("round_summary") or {}
        diagnosis = step6.get("diagnosis") or {}
        rows.append({
            "round_number": r.get("round_number"),
            "decision": step6.get("decision"),
            "assessment_status": (step6.get("round_assessment") or {}).get("status"),
            "quality_gate_passed": step6.get("quality_gate_passed"),
            "total_molecules": summary.get("total_evaluated") or summary.get("total_molecules") or summary.get("num_total"),
            "candidate_pool_size": summary.get("candidate_total") or summary.get("candidate_pool_size"),
            "num_keep": summary.get("keep") if summary.get("keep") is not None else summary.get("num_keep"),
            "num_modify": summary.get("modify") if summary.get("modify") is not None else summary.get("num_modify"),
            "num_reject": summary.get("reject") if summary.get("reject") is not None else summary.get("num_reject"),
            "candidate_pass_rate": summary.get("candidate_pass_rate"),
            "pass_rate": summary.get("pass_rate"),
            "primary_bottleneck": summary.get("primary_bottleneck") or diagnosis.get("bottleneck_type"),
            "suggested_strategy": diagnosis.get("suggested_strategy"),
        })
    return rows


# ---------- helpers ----------

def _safe_json(path: Path) -> Any:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _safe_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return [dict(row) for row in reader]
    except OSError:
        return []


def _existing_or_none(path: Path) -> str | None:
    return str(path) if path.exists() else None


def _first_existing(*candidates) -> str | None:
    for c in candidates:
        if c and Path(c).exists():
            return str(c)
    return None


def _first_matching(directory: Path, pattern: str) -> str | None:
    matches = sorted(directory.glob(pattern))
    return str(matches[0]) if matches else None


def _to_int(value) -> int | None:
    try:
        if value is None or value == "":
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _to_pose_index(value) -> int | None:
    """Parse an artifact index without bools or lossy float truncation."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped and stripped.lstrip("-").isdigit():
            return int(stripped)
    return None


def _to_rank(value) -> float | None:
    try:
        if value is None or value == "" or isinstance(value, bool):
            return None
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _to_float(value) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _first_non_none(*args):
    for a in args:
        if a is not None:
            return a
    return None


def _first_float(*values) -> float | None:
    """Return the first parseable value without treating 0.0 as missing."""
    for value in values:
        parsed = _to_float(value)
        if parsed is not None:
            return parsed
    return None


def _declared_cohort_ids(payload: Any) -> set[str]:
    """Return cohort IDs explicitly declared by one source/report.

    These declarations are used only to support cohort-less legacy rows.  A
    legacy row may be inferred into a candidate cohort when the complete source
    declares exactly that one cohort; an absent or multi-cohort declaration is
    not sufficient evidence.
    """
    if not isinstance(payload, dict):
        return set()

    cohort_ids: set[str] = set()

    def add_manifest(item: Any) -> None:
        if isinstance(item, dict) and _cohort_tuple_valid(item):
            cohort_ids.add(str(item["cohort_id"]))

    scientific_cohort = payload.get("scientific_cohort")
    add_manifest(scientific_cohort)
    for cohort in payload.get("scientific_cohorts") or []:
        add_manifest(cohort)
    return cohort_ids


def _cohort_tuple_valid(item: Any) -> bool:
    if not isinstance(item, dict):
        return False
    if item.get("cohort_manifest_valid") is False:
        return False
    if any(
        item.get(field) in (None, "")
        for field in ("cohort_id", *_COHORT_COMPONENT_FIELDS)
    ):
        return False
    canonical = "|".join(
        f"{field}={item[field]}" for field in _COHORT_COMPONENT_FIELDS
    )
    cohort_id = str(item["cohort_id"])
    return cohort_id == canonical or (
        cohort_id.startswith("declared_cohort_id=")
        and cohort_id.endswith(f"|{canonical}")
    )


def _cohort_tuples_compatible(left: Any, right: Any) -> bool:
    return bool(
        _cohort_tuple_valid(left)
        and _cohort_tuple_valid(right)
        and str(left.get("cohort_id")) == str(right.get("cohort_id"))
        and all(
            str(left.get(field)) == str(right.get(field))
            for field in _COHORT_COMPONENT_FIELDS
        )
    )


def _index_records_by_identity(
    records: list[dict],
    *,
    source_cohort_ids: set[str] | list[str] | tuple[str, ...] | None = None,
    _include_cohortless_index: bool = True,
) -> dict[str, Any]:
    by_pose_id: dict[str, list[dict]] = {}
    by_chemical_identity: dict[str, list[dict]] = {}
    by_identity: dict[tuple[str, int], list[dict]] = {}
    by_name: dict[str, list[dict]] = {}
    by_pose: dict[int, list[dict]] = {}
    by_source_pose: dict[int, list[dict]] = {}
    by_cohort_pose_id: dict[tuple[str, str], list[dict]] = {}
    by_cohort_chemical_identity: dict[tuple[str, str], list[dict]] = {}
    by_cohort_identity: dict[tuple[str, str, int], list[dict]] = {}
    by_cohort_name: dict[tuple[str, str], list[dict]] = {}
    by_cohort_pose: dict[tuple[str, int], list[dict]] = {}
    by_cohort_source_pose: dict[tuple[str, int], list[dict]] = {}
    indexed_records: list[dict] = []
    cohortless_records: list[dict] = []
    declared_source_cohort_ids = {
        str(value)
        for value in (source_cohort_ids or [])
        if value not in (None, "")
    }
    observed_cohort_ids: set[str] = set()
    for record in records or []:
        if not isinstance(record, dict):
            continue
        indexed_records.append(record)
        name = record.get("mol_name") or record.get("name") or record.get("compound_id")
        raw_pose_indices = [
            record.get("pose_index"), record.get("selected_pose_index")
        ]
        declared_pose_indices = [
            value for value in raw_pose_indices if value not in (None, "")
        ]
        parsed_pose_indices = [
            _to_pose_index(value) for value in declared_pose_indices
        ]
        weak_pose_metadata_valid = bool(
            all(value is not None for value in parsed_pose_indices)
            and len(set(parsed_pose_indices)) <= 1
        )
        pose_index = (
            parsed_pose_indices[0]
            if weak_pose_metadata_valid and parsed_pose_indices
            else None
        )
        cohort_id = record.get("cohort_id")
        cohort_id = (
            str(cohort_id) if cohort_id not in (None, "") else None
        )
        if cohort_id is None:
            cohortless_records.append(record)
        else:
            observed_cohort_ids.add(cohort_id)
        if name not in (None, "") and weak_pose_metadata_valid:
            name = str(name)
            by_name.setdefault(name, []).append(record)
            if pose_index is not None:
                by_identity.setdefault((name, pose_index), []).append(record)
            if cohort_id is not None:
                by_cohort_name.setdefault((cohort_id, name), []).append(record)
                if pose_index is not None:
                    by_cohort_identity.setdefault(
                        (cohort_id, name, pose_index), []
                    ).append(record)
        pose_id = record.get("pose_id")
        if pose_id not in (None, ""):
            by_pose_id.setdefault(str(pose_id), []).append(record)
            if cohort_id is not None:
                by_cohort_pose_id.setdefault(
                    (cohort_id, str(pose_id)), []
                ).append(record)
        for alternate_pose_id in record.get("pose_ids") or []:
            if alternate_pose_id not in (None, ""):
                by_pose_id.setdefault(str(alternate_pose_id), []).append(record)
                if cohort_id is not None:
                    by_cohort_pose_id.setdefault(
                        (cohort_id, str(alternate_pose_id)), []
                    ).append(record)
        chemical_identity_id = (
            record.get("chemical_identity_id")
            or ((record.get("chemical_identity") or {}).get("identity_id"))
        )
        if chemical_identity_id not in (None, ""):
            by_chemical_identity.setdefault(
                str(chemical_identity_id), []
            ).append(record)
            if cohort_id is not None:
                by_cohort_chemical_identity.setdefault(
                    (cohort_id, str(chemical_identity_id)), []
                ).append(record)
        source_pose_indices = list(record.get("source_pose_indices") or [])
        source_pose_indices.append(record.get("source_pose_index"))
        normalized_source_indices = (
            {
                source_index
                for source_index in (
                    _to_pose_index(source_pose_index)
                    for source_pose_index in source_pose_indices
                )
                if source_index is not None
            }
            if weak_pose_metadata_valid
            else set()
        )
        for source_index in normalized_source_indices:
            by_source_pose.setdefault(source_index, []).append(record)
            if cohort_id is not None:
                by_cohort_source_pose.setdefault(
                    (cohort_id, source_index), []
                ).append(record)
        if pose_index is not None:
            by_pose.setdefault(pose_index, []).append(record)
            if cohort_id is not None:
                by_cohort_pose.setdefault(
                    (cohort_id, pose_index), []
                ).append(record)
    index: dict[str, Any] = {
        "by_pose_id": by_pose_id,
        "by_chemical_identity": by_chemical_identity,
        "by_identity": by_identity,
        "by_name": by_name,
        "by_pose": by_pose,
        "by_source_pose": by_source_pose,
        "by_cohort_pose_id": by_cohort_pose_id,
        "by_cohort_chemical_identity": by_cohort_chemical_identity,
        "by_cohort_identity": by_cohort_identity,
        "by_cohort_name": by_cohort_name,
        "by_cohort_pose": by_cohort_pose,
        "by_cohort_source_pose": by_cohort_source_pose,
        "source_cohort_ids": sorted(declared_source_cohort_ids),
        "observed_cohort_ids": sorted(observed_cohort_ids),
        "legacy_inference_cohort_id": (
            next(iter(declared_source_cohort_ids))
            if (
                len(declared_source_cohort_ids) == 1
                and observed_cohort_ids.issubset(declared_source_cohort_ids)
            )
            else None
        ),
        "records": indexed_records,
    }
    if _include_cohortless_index:
        index["cohortless_index"] = _index_records_by_identity(
            cohortless_records,
            source_cohort_ids=declared_source_cohort_ids,
            _include_cohortless_index=False,
        )
    return index


def _lookup_identity(
    index: dict[str, Any],
    mol_name: str,
    pose_index: int | None,
    *,
    pose_id: str | None = None,
    chemical_identity_id: str | None = None,
    sdf_molblock_sha256: str | None = None,
    source_pose_index: int | None = None,
    cohort_id: str | None = None,
    cohort_identity: dict | None = None,
) -> dict:
    """Prefer chemical/pose identities, then unambiguous exact fallbacks."""
    inferred_cohortless = False

    def accepted(candidate: Any) -> dict:
        if not isinstance(candidate, dict):
            return {}
        candidate_witness = candidate.get("sdf_molblock_sha256")
        if (
            candidate_witness not in (None, "")
            or sdf_molblock_sha256 not in (None, "")
        ) and (
            candidate_witness in (None, "")
            or sdf_molblock_sha256 in (None, "")
            or str(candidate_witness) != str(sdf_molblock_sha256)
        ):
            return {}
        if inferred_cohortless or not _cohort_tuple_valid(cohort_identity):
            return candidate
        return (
            candidate
            if _cohort_tuples_compatible(cohort_identity, candidate)
            else {}
        )

    def weakly_accepted(
        candidate: Any,
        *,
        require_supplied_stable_ids: bool = False,
    ) -> dict:
        """Validate labels/IDs before accepting an integer or name fallback."""
        candidate = accepted(candidate)
        if not candidate:
            return {}
        candidate_name = (
            candidate.get("mol_name")
            or candidate.get("name")
            or candidate.get("compound_id")
        )
        if mol_name not in (None, "") and (
            candidate_name in (None, "")
            or str(candidate_name) != str(mol_name)
        ):
            return {}
        alternate_pose_ids = candidate.get("pose_ids")
        if isinstance(alternate_pose_ids, (list, tuple, set)):
            alternates = list(alternate_pose_ids)
        elif alternate_pose_ids in (None, ""):
            alternates = []
        else:
            alternates = [alternate_pose_ids]
        candidate_pose_ids = {
            str(value)
            for value in [
                candidate.get("pose_id"),
                *alternates,
            ]
            if value not in (None, "")
        }
        if pose_id not in (None, ""):
            if candidate_pose_ids:
                if str(pose_id) not in candidate_pose_ids:
                    return {}
            elif require_supplied_stable_ids:
                return {}
        candidate_chemical_identity_id = (
            candidate.get("chemical_identity_id")
            or (
                (candidate.get("chemical_identity") or {}).get("identity_id")
                if isinstance(candidate.get("chemical_identity"), dict)
                else None
            )
        )
        if (
            chemical_identity_id not in (None, "")
            and (
                (
                    require_supplied_stable_ids
                    and candidate_chemical_identity_id in (None, "")
                )
                or (
                    candidate_chemical_identity_id not in (None, "")
                    and str(candidate_chemical_identity_id)
                    != str(chemical_identity_id)
                )
            )
        ):
            return {}
        return candidate

    def strongly_accepted(candidate: Any, *, matched_by: str) -> dict:
        """Accept a stable-ID hit only when every explicit ID is compatible.

        A chemical identity identifies a molecule rather than one concrete
        pose.  Therefore, when the query also supplies ``pose_id``, a hit
        reached through the chemical namespace must carry and match that pose
        ID.  A direct pose-ID hit may tolerate a missing legacy chemical ID,
        but never an explicit conflicting one.
        """
        candidate = accepted(candidate)
        if not candidate:
            return {}
        alternate_pose_ids = candidate.get("pose_ids")
        if isinstance(alternate_pose_ids, (list, tuple, set)):
            alternates = list(alternate_pose_ids)
        elif alternate_pose_ids in (None, ""):
            alternates = []
        else:
            alternates = [alternate_pose_ids]
        candidate_pose_ids = {
            str(value)
            for value in [candidate.get("pose_id"), *alternates]
            if value not in (None, "")
        }
        if pose_id not in (None, ""):
            if candidate_pose_ids and str(pose_id) not in candidate_pose_ids:
                return {}
            if matched_by != "pose_id" and not candidate_pose_ids:
                return {}
        candidate_chemical_identity_id = (
            candidate.get("chemical_identity_id")
            or (
                (candidate.get("chemical_identity") or {}).get("identity_id")
                if isinstance(candidate.get("chemical_identity"), dict)
                else None
            )
        )
        if (
            chemical_identity_id not in (None, "")
            and candidate_chemical_identity_id not in (None, "")
            and str(candidate_chemical_identity_id)
            != str(chemical_identity_id)
        ):
            return {}
        return candidate

    cohort_id = str(cohort_id) if cohort_id not in (None, "") else None
    if cohort_id == "cohort:unspecified":
        # This is a presentation sentinel created by _build_candidates, not a
        # provenance assertion.  Treat it as missing so entirely cohort-less
        # legacy artifacts can still join by a proven pose/name identity.
        cohort_id = None
    if cohort_id is not None:
        if chemical_identity_id not in (None, ""):
            matches = index.get("by_cohort_chemical_identity", {}).get(
                (cohort_id, str(chemical_identity_id)), []
            )
            if len(matches) == 1:
                candidate = strongly_accepted(
                    matches[0], matched_by="chemical_identity_id"
                )
                if candidate:
                    return candidate
        if pose_id not in (None, ""):
            exact_poses = index.get("by_cohort_pose_id", {}).get(
                (cohort_id, str(pose_id)), []
            )
            if len(exact_poses) == 1:
                candidate = strongly_accepted(
                    exact_poses[0], matched_by="pose_id"
                )
                if candidate:
                    return candidate
        if source_pose_index is not None:
            matches = index.get("by_cohort_source_pose", {}).get(
                (cohort_id, int(source_pose_index)), []
            )
            if len(matches) == 1:
                candidate = weakly_accepted(matches[0])
                if candidate:
                    return candidate
        if pose_index is not None:
            matches = index.get("by_cohort_identity", {}).get(
                (cohort_id, mol_name, pose_index), []
            )
            if len(matches) == 1:
                candidate = weakly_accepted(matches[0])
                if candidate:
                    return candidate
            matches = index.get("by_cohort_pose", {}).get(
                (cohort_id, pose_index), []
            )
            if len(matches) == 1:
                candidate = weakly_accepted(matches[0])
                if candidate:
                    return candidate
        matches = index.get("by_cohort_name", {}).get(
            (cohort_id, mol_name), []
        )
        if len(matches) == 1:
            candidate = weakly_accepted(
                matches[0], require_supplied_stable_ids=True
            )
            if candidate:
                return candidate
        # A candidate with a cohort must never fall back to a row carrying a
        # different cohort.  Cohort-less legacy rows are eligible only when
        # the complete source/report proves it contains exactly this cohort.
        if index.get("legacy_inference_cohort_id") != cohort_id:
            return {}
        cohortless_index = index.get("cohortless_index")
        if not isinstance(cohortless_index, dict):
            return {}
        index = cohortless_index
        inferred_cohortless = True
    if chemical_identity_id not in (None, ""):
        chemical_matches = index.get("by_chemical_identity", {}).get(
            str(chemical_identity_id), []
        )
        if len(chemical_matches) == 1:
            candidate = strongly_accepted(
                chemical_matches[0], matched_by="chemical_identity_id"
            )
            if candidate:
                return candidate
    if pose_id not in (None, ""):
        exact_poses = index.get("by_pose_id", {}).get(str(pose_id), [])
        if len(exact_poses) == 1:
            candidate = strongly_accepted(
                exact_poses[0], matched_by="pose_id"
            )
            if candidate:
                return candidate
    if source_pose_index is not None:
        source_matches = index.get("by_source_pose", {}).get(
            int(source_pose_index), []
        )
        if len(source_matches) == 1:
            candidate = weakly_accepted(source_matches[0])
            if candidate:
                return candidate
    if pose_index is not None:
        exact_matches = index.get("by_identity", {}).get(
            (mol_name, pose_index), []
        )
        if len(exact_matches) == 1:
            candidate = weakly_accepted(exact_matches[0])
            if candidate:
                return candidate
        pose_matches = index.get("by_pose", {}).get(pose_index, [])
        if len(pose_matches) == 1:
            candidate = weakly_accepted(pose_matches[0])
            if candidate:
                return candidate
    name_matches = index.get("by_name", {}).get(mol_name, [])
    if len(name_matches) == 1:
        candidate = weakly_accepted(
            name_matches[0], require_supplied_stable_ids=True
        )
        if candidate:
            return candidate
    return {}


def _index_docking_scores_by_pose(scores: list[dict]) -> dict[int, dict]:
    """Legacy helper retained for callers outside this module.

    Internal report joins use ``_index_records_by_identity`` so duplicate
    compound names cannot silently overwrite one another.
    """
    indexed: dict[int, dict] = {}
    for s in scores or []:
        if not isinstance(s, dict):
            continue
        pi = _to_pose_index(s.get("pose_index"))
        if pi is None:
            continue
        if pi not in indexed:
            indexed[pi] = s
    return indexed


def _infer_round_number(project_dir: Path) -> int | None:
    match = re.fullmatch(r"round(\d+)", project_dir.name)
    return int(match.group(1)) if match else None


def _discover_rounds(project_root: Path) -> list[tuple[int, Path]]:
    rounds: list[tuple[int, Path]] = []
    for child in sorted(project_root.iterdir()):
        if not child.is_dir():
            continue
        match = re.fullmatch(r"round(\d+)", child.name)
        if match:
            rounds.append((int(match.group(1)), child))
    rounds.sort(key=lambda x: x[0])
    return rounds


def _pick_candidate_source_round(per_round: list[dict]) -> tuple[int | None, str]:
    """Choose the last round with eligible evidence, independent of campaign stop.

    Legacy reports without an explicit nomination field retain the old behavior,
    except that ``terminate_invalid`` is always excluded.  New structured reports
    fail closed when they explicitly declare ``valid_for_nomination=false``.
    """
    for r in reversed(per_round):
        step6 = r.get("step6") or {}
        report = step6.get("coordinator_report") or {}
        if not step6.get("success") or not report:
            continue
        if report.get("assessment_schema_version") == "round-assessment/v1":
            assessment = report.get("round_assessment") or {}
            if (
                assessment.get("status") in {"valid", "review_required"}
                and assessment.get("valid_for_nomination") is True
            ):
                return r.get("round_number"), "last nomination-valid round assessment"
            continue
        explicit = step6.get("valid_for_nomination")
        if explicit is None and isinstance(report, dict):
            explicit = report.get("valid_for_nomination")
        outcome = step6.get("decision_outcome")
        if outcome is None and isinstance(report, dict):
            outcome = report.get("decision_outcome")
        if outcome is not None or explicit is not None:
            # Once either structured field is present, fail closed unless the
            # complete success contract is internally consistent.
            if outcome != "success" or explicit is not True:
                continue
        elif (
            step6.get("decision") == "terminate_invalid"
            or report.get("decision") == "terminate_invalid"
        ):
            continue
        return r.get("round_number"), "last nomination-valid successful step6"
    if per_round:
        return None, (
            "no round completed a nomination-valid step6; candidate nomination "
            "is suppressed"
        )
    return None, "no rounds discovered"


def _find_round_by_number(per_round: list[dict], number: int | None) -> dict | None:
    if number is None:
        return None
    for r in per_round:
        if r.get("round_number") == number:
            return r
    return None


def _guess_title(metadata: dict, coord: dict) -> str:
    target = (metadata or {}).get("target_name") or (coord or {}).get("target_name")
    pdb = (metadata or {}).get("pdb_id") or (coord or {}).get("pdb_id")
    if target and pdb:
        return f"SBDD report — {target} ({pdb})"
    if target:
        return f"SBDD report — {target}"
    if pdb:
        return f"SBDD report — {pdb}"
    return "SBDD report"


def _collect_missing(data: dict) -> list[dict]:
    missing = list(data.get("missing_sources") or [])
    for step_key in ("step0", "step1", "step3", "step4", "step5", "step6"):
        step = data.get(step_key)
        if isinstance(step, dict) and step.get("available") is False:
            missing.append({"step": step_key, "reason": step.get("reason", "unavailable")})
    return missing

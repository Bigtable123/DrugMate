"""Verify the committed Step 3 docking artifact set before consumption.

The docking report is the commit marker for its SDF, score table, and optional
native-redocking pose.  Consumers must not trust duplicated score/RMSD fields
when that marker or any sealed file identity fails verification.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping


REDOCKING_SEALED_FIELDS = (
    "performed",
    "measured",
    "applicable",
    "not_applicable",
    "passed",
    "gate",
    "best_rmsd_angstrom",
    "best_rmsd_pose_index",
    "rmsd_threshold_angstrom",
    "rmsd_soft_threshold_angstrom",
    "rmsd_hard_threshold_angstrom",
)


def redocking_measurement_error(
    measurement: object,
    *,
    redocked_reference_present: bool,
) -> str | None:
    """Return why a sealed native-redocking state is internally impossible.

    A completed measurement is evidence only when the committed coordinates,
    RMSD, thresholds, and verdict describe one coherent state.  Conversely, a
    tool invocation that produced no RMSD is a valid N/A state: ``performed``
    records the attempt, while ``measured=False`` prevents it from becoming a
    protocol pass/fail verdict.
    """
    if not (
        isinstance(measurement, Mapping)
        and set(measurement) == set(REDOCKING_SEALED_FIELDS)
    ):
        return "redocking measurement seal has an incomplete field set"

    performed = measurement.get("performed")
    measured = measurement.get("measured")
    applicable = measurement.get("applicable")
    not_applicable = measurement.get("not_applicable")
    passed = measurement.get("passed")
    gate = measurement.get("gate")
    rmsd_value = measurement.get("best_rmsd_angstrom")
    pose_index = measurement.get("best_rmsd_pose_index")

    if measured is True:
        if not redocked_reference_present:
            return "measured redocking has no committed redocked SDF"
        if performed is not True:
            return "measured redocking is not marked as performed"
        if applicable is not True or not_applicable is not False:
            return "measured redocking has inconsistent applicability flags"
        if gate not in {"keep", "modify", "reject"}:
            return "measured redocking has no valid three-tier gate"
        if not isinstance(passed, bool) or passed is not (gate == "keep"):
            return "redocking passed flag disagrees with its gate"

        numeric: dict[str, float] = {}
        for field in (
            "best_rmsd_angstrom",
            "rmsd_threshold_angstrom",
            "rmsd_soft_threshold_angstrom",
            "rmsd_hard_threshold_angstrom",
        ):
            value = measurement.get(field)
            if isinstance(value, bool):
                return f"redocking {field} is not a finite number"
            try:
                parsed = float(value)
            except (TypeError, ValueError):
                return f"redocking {field} is not a finite number"
            if not math.isfinite(parsed) or parsed < 0:
                return f"redocking {field} is not a non-negative finite number"
            numeric[field] = parsed

        soft = numeric["rmsd_soft_threshold_angstrom"]
        hard = numeric["rmsd_hard_threshold_angstrom"]
        nominal = numeric["rmsd_threshold_angstrom"]
        if soft > hard:
            return "redocking soft threshold exceeds its hard threshold"
        if nominal != soft:
            return "redocking nominal threshold disagrees with its soft threshold"
        rmsd = numeric["best_rmsd_angstrom"]
        expected_gate = (
            "keep" if rmsd <= soft
            else "modify" if rmsd <= hard
            else "reject"
        )
        if gate != expected_gate:
            return "redocking gate disagrees with RMSD and sealed thresholds"
        index = _canonical_index(pose_index)
        if index is None or index < 0:
            return "measured redocking has no valid best-pose index"
        return None

    # No measurement means no molecular/protocol verdict.  ``performed=True``
    # is intentionally allowed here: it distinguishes an attempted tool run
    # from a disabled/not-selected one without converting failure into reject.
    if redocked_reference_present:
        return "committed redocked SDF is paired with no RMSD measurement"
    if measured not in {None, False}:
        return "redocking measured flag is not tri-state"
    if performed not in {None, False, True}:
        return "redocking performed flag is not tri-state"
    if performed is True and measured is not False:
        return "attempted redocking does not explicitly record measured=false"
    if applicable not in {None, False} or not_applicable not in {None, True}:
        return "unmeasured redocking has inconsistent applicability flags"
    if passed is not None:
        return "unmeasured redocking carries a pass/fail verdict"
    if gate not in {None, "not_applicable"}:
        return "unmeasured redocking carries a three-tier gate verdict"
    if rmsd_value is not None or pose_index is not None:
        return "unmeasured redocking carries RMSD or pose-index evidence"
    return None


def _invalid(report: object, reason: str) -> dict[str, Any]:
    return {
        "available": False,
        "reason": reason,
        "docking_report": report if isinstance(report, dict) else {},
        "docking_scores": [],
        "docked_poses_sdf": None,
        "redocked_reference_sdf": None,
        "redocking_validation": None,
        "redocking_evidence_status": "invalid",
    }


def _read_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _artifact_identity_matches(path: Path, expected: object) -> bool:
    if not isinstance(expected, Mapping) or not path.is_file():
        return False
    expected_digest = expected.get("sha256")
    expected_size = expected.get("size_bytes")
    if not isinstance(expected_digest, str) or len(expected_digest) != 64:
        return False
    try:
        int(expected_digest, 16)
        if isinstance(expected_size, bool) or int(expected_size) != path.stat().st_size:
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
    return digest.hexdigest() == expected_digest.casefold()


def _physical_sdf_records(path: Path) -> list[bytes] | None:
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


def _canonical_index(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped and stripped.lstrip("-").isdigit():
            return int(stripped)
    return None


def _unavailable_redocking(redocking: object, reason: str) -> dict[str, Any]:
    sanitized = dict(redocking) if isinstance(redocking, Mapping) else {}
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


def evaluation_docking_lineage_error(
    evaluation_report: object,
    docking_report: object,
    *,
    minimum_schema: int = 4,
) -> str | None:
    """Return why a Step 4 report is not bound to this Step 3 set.

    The ligand pose hash alone is insufficient because the same coordinates
    can be evaluated against a different receptor.  Schema 4 seals the
    evaluation receptor and score protocol in addition to the pose/score
    artifacts; consumers must compare the complete lineage before merging any
    Step 3 score with Step 4 QC evidence.
    """
    lineage = (
        evaluation_report.get("source_docking_artifact_set")
        if isinstance(evaluation_report, Mapping) else None
    )
    contract = (
        docking_report.get("docked_pose_artifact_contract")
        if isinstance(docking_report, Mapping) else None
    )
    if not isinstance(contract, Mapping):
        return "committed docking report has no artifact-set contract"
    try:
        contract_schema = int(contract.get("artifact_set_schema_version"))
    except (TypeError, ValueError):
        return "committed docking report has an invalid artifact-set schema"
    if contract_schema < minimum_schema:
        return (
            "Step 3 artifact set does not seal the evaluation receptor and "
            "score protocol; rerun docking"
        )
    if not isinstance(lineage, Mapping) or lineage.get("status") != "verified":
        return (
            "evaluation report does not carry verified Step 3 artifact-set "
            "lineage; rerun structure-evaluator with docking_report_file"
        )
    if lineage.get("artifact_set_id") != contract.get("artifact_set_id"):
        return (
            "evaluation report belongs to a different Step 3 docking "
            "artifact set"
        )
    try:
        lineage_schema = int(lineage.get("artifact_set_schema_version"))
    except (TypeError, ValueError):
        return "evaluation report has an invalid Step 3 artifact-set schema"
    if lineage_schema != contract_schema:
        return (
            "evaluation report Step 3 artifact-set schema does not match "
            "docking_report_file"
        )

    for field in (
        "docked_sdf", "docking_scores_json", "evaluation_receptor",
    ):
        observed = lineage.get(field)
        expected = contract.get(field)
        if not isinstance(observed, Mapping) or not isinstance(expected, Mapping):
            return f"evaluation report lacks the sealed {field} identity"
        try:
            identity_matches = bool(
                observed.get("sha256") == expected.get("sha256")
                and not isinstance(observed.get("size_bytes"), bool)
                and not isinstance(expected.get("size_bytes"), bool)
                and int(observed.get("size_bytes"))
                == int(expected.get("size_bytes"))
            )
        except (TypeError, ValueError):
            identity_matches = False
        if not identity_matches:
            return (
                f"evaluation report {field} identity does not match the "
                "committed Step 3 artifact set"
            )
    if lineage.get("score_protocol_fingerprint_sha256") != contract.get(
        "score_protocol_fingerprint_sha256"
    ):
        return (
            "evaluation report score-protocol identity does not match the "
            "committed Step 3 artifact set"
        )
    return None


def load_committed_docking_artifacts(report_file: Path) -> dict[str, Any]:
    """Load Step 3 only when its report and every declared artifact agree."""
    report_file = Path(report_file)
    if not report_file.is_file():
        return _invalid({}, f"docking report missing: {report_file}")
    report = _read_json(report_file)
    contract = (
        report.get("docked_pose_artifact_contract")
        if isinstance(report, dict) else None
    )
    if not (
        isinstance(report, dict)
        and report.get("success") is True
        and isinstance(contract, dict)
        and contract.get("status") == "passed"
        and contract.get("committed") is True
    ):
        return _invalid(
            report, "docking artifact set is not committed and verified"
        )

    files = report.get("files") or {}
    docked_path = Path(str(files.get("docked_poses_sdf") or ""))
    scores_path = Path(str(files.get("docking_scores_json") or ""))
    expected_sdf = contract.get("docked_sdf") or {}
    expected_scores = contract.get("docking_scores_json") or {}
    try:
        schema = int(contract.get("artifact_set_schema_version", 1))
    except (TypeError, ValueError):
        schema = -1
    if schema not in {1, 2, 3, 4}:
        return _invalid(report, "unsupported docking artifact-set schema")

    seal_payload: dict[str, Any] = {
        "docked_sdf_sha256": expected_sdf.get("sha256"),
        "docking_scores_sha256": expected_scores.get("sha256"),
        "physical_sdf_records": contract.get("physical_sdf_records"),
        "score_records": contract.get("score_records"),
    }
    expected_redocked = contract.get("redocked_reference_sdf") or {}
    redocked_present = contract.get("redocked_reference_present") is True
    if schema >= 2:
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
    sealed_measurement = contract.get("redocking_measurement")
    if schema == 3:
        if not (
            isinstance(sealed_measurement, dict)
            and set(sealed_measurement) == set(REDOCKING_SEALED_FIELDS)
        ):
            return _invalid(
                report, "docking artifact set lacks its redocking measurement seal"
            )
        seal_payload["redocking_measurement"] = sealed_measurement
    elif schema == 4:
        if not (
            isinstance(sealed_measurement, dict)
            and set(sealed_measurement) == set(REDOCKING_SEALED_FIELDS)
        ):
            return _invalid(
                report, "docking artifact set lacks its redocking measurement seal"
            )
        evaluation_receptor = contract.get("evaluation_receptor")
        protocol_fingerprint = contract.get(
            "score_protocol_fingerprint_sha256"
        )
        reported_fingerprint = (
            (report.get("score_protocol") or {}).get("fingerprint_sha256")
            if isinstance(report.get("score_protocol"), Mapping) else None
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
            isinstance(evaluation_receptor, Mapping)
            and protocol_fingerprint_valid
            and protocol_fingerprint == reported_fingerprint
            and _artifact_identity_matches(
                Path(str(evaluation_receptor.get("path") or "")),
                evaluation_receptor,
            )
        ):
            return _invalid(
                report,
                "docking artifact set lacks a verified evaluation receptor "
                "or score-protocol identity",
            )
        seal_payload.update({
            "redocking_measurement": sealed_measurement,
            "evaluation_receptor_sha256": evaluation_receptor.get("sha256"),
            "evaluation_receptor_size_bytes": evaluation_receptor.get(
                "size_bytes"
            ),
            "score_protocol_fingerprint_sha256": protocol_fingerprint,
        })
    expected_set_id = "sha256:" + hashlib.sha256(
        json.dumps(
            seal_payload, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    if not (
        _artifact_identity_matches(docked_path, expected_sdf)
        and _artifact_identity_matches(scores_path, expected_scores)
        and contract.get("artifact_set_id") == expected_set_id
    ):
        return _invalid(
            report, "docking artifact hash/size does not match commit marker"
        )

    scores_payload = _read_json(scores_path)
    scores = (
        scores_payload
        if isinstance(scores_payload, list)
        else scores_payload.get("scores")
        if isinstance(scores_payload, dict)
        else None
    )
    physical_records = _physical_sdf_records(docked_path)
    if not (
        isinstance(scores, list)
        and physical_records is not None
        and len(physical_records) == contract.get("physical_sdf_records")
        and len(scores) == contract.get("score_records")
    ):
        return _invalid(
            report, "docking artifact record counts do not match commit marker"
        )

    redocking = report.get("redocking_validation")
    redocked_sdf: str | None = None
    redocking_status = "not_applicable"
    if schema in {3, 4}:
        reported_measurement = {
            field: redocking.get(field)
            if isinstance(redocking, dict) else None
            for field in REDOCKING_SEALED_FIELDS
        }
        if json.dumps(
            reported_measurement, sort_keys=True, separators=(",", ":")
        ) != json.dumps(
            sealed_measurement, sort_keys=True, separators=(",", ":")
        ):
            redocking = _unavailable_redocking(
                redocking,
                "redocking measurement fields do not match the committed "
                "docking report seal",
            )
            redocking_status = "invalid"

    if schema in {3, 4} and redocking_status != "invalid":
        semantic_error = redocking_measurement_error(
            sealed_measurement,
            redocked_reference_present=redocked_present,
        )
        if semantic_error is not None:
            redocking = _unavailable_redocking(
                redocking,
                f"invalid committed redocking measurement: {semantic_error}",
            )
            redocking_status = "invalid"

    if schema in {3, 4} and redocked_present and redocking_status != "invalid":
        validation_identity = (
            redocking.get("redocked_sdf_artifact")
            if isinstance(redocking, dict) else None
        )
        redocked_path_value = files.get("redocked_reference_sdf")
        valid = bool(
            isinstance(redocking, dict)
            and redocking.get("performed") is True
            and redocking.get("measured") is True
            and redocking.get("artifact_available") is True
            and isinstance(expected_redocked, dict)
            and validation_identity == expected_redocked
            and redocked_path_value not in (None, "")
        )
        try:
            rmsd = float(redocking.get("best_rmsd_angstrom"))
            valid = valid and math.isfinite(rmsd)
        except (AttributeError, TypeError, ValueError):
            valid = False
        redocked_path = Path(str(redocked_path_value or ""))
        try:
            valid = bool(
                valid
                and redocked_path.resolve()
                == Path(str(expected_redocked.get("path") or "")).resolve()
                == Path(str(validation_identity.get("path") or "")).resolve()
            )
        except (AttributeError, OSError):
            valid = False
        blocks = _physical_sdf_records(redocked_path) if valid else None
        best_index = _canonical_index(expected_redocked.get(
            "best_rmsd_pose_index"
        ))
        expected_count = _canonical_index(expected_redocked.get(
            "physical_sdf_records"
        ))
        expected_record_digest = expected_redocked.get(
            "best_pose_sdf_record_sha256"
        )
        valid = bool(
            valid
            and _artifact_identity_matches(redocked_path, expected_redocked)
            and blocks is not None
            and expected_count == len(blocks)
            and best_index is not None
            and 0 <= best_index < len(blocks)
            and isinstance(expected_record_digest, str)
            and hashlib.sha256(blocks[best_index]).hexdigest()
            == expected_record_digest
            and _canonical_index(redocking.get("best_rmsd_pose_index"))
            == best_index
        )
        if valid:
            redocked_sdf = str(redocked_path)
            redocking_status = "verified"
        else:
            redocking = _unavailable_redocking(
                redocking,
                "redocking artifact or best-pose identity does not match the "
                "committed docking report",
            )
            redocking_status = "invalid"
    elif (
        redocking_status != "invalid"
        and schema in {1, 2}
        and isinstance(redocking, dict)
        and (
            redocking.get("performed") is True
            or redocking.get("measured") is True
            or redocking.get("passed") is not None
            or redocking.get("best_rmsd_angstrom") is not None
            or redocking.get("rmsd") is not None
        )
    ):
        redocking = _unavailable_redocking(
            redocking,
            "docking artifact set does not commit a verifiable redocking "
            "measurement",
        )
        redocking_status = "invalid"
    elif isinstance(redocking, dict):
        redocking = _unavailable_redocking(
            redocking,
            str(redocking.get("reason") or "redocking was not performed"),
        )

    return {
        "available": True,
        "docking_report": report,
        "docking_scores": scores,
        "docked_poses_sdf": str(docked_path),
        "redocked_reference_sdf": redocked_sdf,
        "redocking_validation": redocking,
        "redocking_evidence_status": redocking_status,
        "evaluation_receptor_file": (
            str((contract.get("evaluation_receptor") or {}).get("path"))
            if schema >= 4 else None
        ),
    }

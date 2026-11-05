"""Match pose interactions against legacy and schema-v2 interaction blueprints.

Required constraints are hard constraints: every independent group must pass.
Detailed schema-v2 constraints are evaluated conservatively with three states:
``matched``, ``missing``, and ``unmatchable``.  The last state is important when
the pose evaluator (currently PLIP) cannot observe a ProLIF/Arpeggio-specific
atom, role, geometry, or interaction class; lack of capability is not evidence
that the pose lacks the interaction.
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path
from typing import Any

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from constants import (
    BLUEPRINT_ANGLE_TOLERANCE_DEGREES,
    BLUEPRINT_CONFIDENCE_WEIGHTS,
    BLUEPRINT_DESIRED_DEFAULT_WEIGHT,
    BLUEPRINT_DISTANCE_TOLERANCE_ANGSTROM,
    normalize_interaction_type,
)
from residue_utils import (
    legacy_residue_label,
    make_residue_record,
    normalize_residue as _shared_normalize_residue,
    parse_residue_record,
)


# The Step-5 evaluator currently obtains pose interactions from PLIP.  This
# declaration describes the wrapper below, not every capability of PLIP itself.
# Callers can provide a richer manifest through ``evaluator_capabilities``.
DEFAULT_PLIP_CAPABILITIES = {
    "source": "plip",
    "interaction_types": {
        "H-bond",
        "Hydrophobic",
        "Pi-stacking",
        "Pi-cation",
        "Salt-bridge",
        "Water-bridge",
        "Halogen-bond",
        "Metal-complex",
    },
    "fields_by_type": {
        "H-bond": {
            "ligand_role", "protein_role", "protein_atoms", "ligand_atoms", "mediator",
            "geometry.distance", "geometry.heavy_atom_distance",
            "geometry.acceptor_hydrogen_distance",
            "geometry.angle", "geometry.donor_angle",
            "geometry.distance_unit", "geometry.angle_unit",
            "geometry.distance_metric",
        },
        "Hydrophobic": {
            "protein_atoms", "ligand_atoms", "mediator", "geometry.distance",
            "geometry.distance_unit", "geometry.distance_metric",
        },
        "Pi-stacking": {
            "protein_atoms", "ligand_atoms", "ligand_role", "protein_role", "mediator",
            "geometry.distance", "geometry.angle", "geometry.distance_unit",
            "geometry.angle_unit", "geometry.distance_metric",
        },
        "Pi-cation": {
            "protein_atoms", "ligand_atoms", "ligand_role", "protein_role", "mediator",
            "geometry.distance", "geometry.distance_unit", "geometry.distance_metric",
        },
        "Salt-bridge": {
            "protein_atoms", "ligand_atoms", "ligand_role", "protein_role", "mediator",
            "geometry.distance", "geometry.distance_unit", "geometry.distance_metric",
        },
        "Water-bridge": {
            "protein_atoms", "ligand_atoms", "ligand_role", "protein_role", "mediator",
            "geometry.distance", "geometry.donor_angle", "geometry.water_angle",
            "geometry.distance_unit", "geometry.angle_unit", "geometry.distance_metric",
        },
        "Halogen-bond": {
            "protein_atoms", "ligand_atoms", "ligand_role", "protein_role", "mediator",
            "geometry.distance", "geometry.donor_angle", "geometry.acceptor_angle",
            "geometry.distance_unit", "geometry.angle_unit", "geometry.distance_metric",
        },
        "Metal-complex": {
            "protein_atoms", "ligand_atoms", "ligand_role", "protein_role", "mediator",
            "geometry.distance", "geometry.distance_unit", "geometry.distance_metric",
        },
    },
}

_DETAIL_FIELDS = (
    "ligand_role",
    "protein_role",
    "protein_atoms",
    "ligand_atoms",
    "mediator",
    "geometry",
    "distance_range",
    "angle_range",
)


def load_blueprint(blueprint_file: Path) -> dict:
    """Load interaction blueprint JSON."""
    with open(blueprint_file) as handle:
        return json.load(handle)


def match_pose_to_blueprint(
    pose_interactions: list[dict],
    blueprint: dict,
    evaluator_capabilities: dict | None = None,
) -> dict:
    """Compare one pose to a blueprint without converting unknowns to misses.

    Legacy blueprints retain residue/type matching.  For ``schema_version >= 2``
    any atom, role, geometry, or mediator fields present on a constraint become
    part of the match contract.
    """
    validation = blueprint.get("validation") or {}
    schema_version = _schema_version(blueprint)
    strict_details = schema_version >= 2 or bool(
        (blueprint.get("match_policy") or {}).get("strict_details")
    )
    blueprint_valid = validation.get("blueprint_valid")
    if blueprint_valid is None:
        blueprint_valid = not bool(validation.get("invalid"))
    blueprint_applicable = not bool(
        validation.get("not_applicable") is True
        or validation.get("applicable") is False
        or blueprint.get("not_applicable") is True
        or blueprint.get("applicable") is False
    )

    capabilities = _normalize_capabilities(evaluator_capabilities)
    required = [item for item in blueprint.get("required", []) if isinstance(item, dict)]
    desired = [item for item in blueprint.get("desired", []) if isinstance(item, dict)]
    optional = [item for item in blueprint.get("optional", []) if isinstance(item, dict)]

    required_groups = _required_groups(blueprint, required)
    required_results_by_key: dict[str, dict] = {}
    for group in required_groups:
        for member in group["members"]:
            key = _constraint_key(member)
            if key not in required_results_by_key:
                required_results_by_key[key] = _evaluate_constraint(
                    member,
                    pose_interactions,
                    capabilities,
                    strict_details=strict_details,
                )

    group_results = []
    for group in required_groups:
        members = [required_results_by_key[_constraint_key(item)] for item in group["members"]]
        group_results.append(_evaluate_group(group, members))

    required_results = list(required_results_by_key.values())
    desired_results = [
        _evaluate_constraint(item, pose_interactions, capabilities, strict_details=strict_details)
        for item in desired
    ]
    optional_results = [
        _evaluate_constraint(item, pose_interactions, capabilities, strict_details=strict_details)
        for item in optional
    ]

    required_hits = [r for r in required_results if r["status"] == "matched"]
    required_misses = [r for r in required_results if r["status"] == "missing"]
    required_unresolved = [r for r in required_results if r["status"] == "unmatchable"]
    desired_hits = [r for r in desired_results if r["status"] == "matched"]
    desired_misses = [r for r in desired_results if r["status"] == "missing"]
    desired_unresolved = [r for r in desired_results if r["status"] == "unmatchable"]
    optional_hits = [r for r in optional_results if r["status"] == "matched"]

    required_total = len(required_results)
    required_hit_rate = len(required_hits) / required_total if required_total else None
    desired_score = _desired_score(desired_results)
    declared_desired_metric = validation.get("desired_metric_applicable")
    desired_metric_applicable = bool(
        blueprint_applicable
        and blueprint_valid
        and desired
        and (
            bool(declared_desired_metric)
            if declared_desired_metric is not None else True
        )
    )
    optional_metric_applicable = bool(
        blueprint_applicable and blueprint_valid and optional
    )

    declared_hard_gate = _declared_hard_gate(validation, bool(required_groups))
    if not blueprint_applicable:
        verdict = "not_applicable"
        hard_gate_applicable = False
        hard_gate_status = "not_applicable"
    elif not blueprint_valid:
        verdict = "unresolved"
        hard_gate_applicable = False
        hard_gate_status = "invalid_blueprint"
    elif not required_groups or not declared_hard_gate:
        verdict = "not_applicable"
        hard_gate_applicable = False
        hard_gate_status = "not_applicable"
    elif any(group["status"] == "failed" for group in group_results):
        verdict = "fail"
        hard_gate_applicable = True
        hard_gate_status = "failed"
    elif any(group["status"] == "unresolved" for group in group_results):
        verdict = "unresolved"
        hard_gate_applicable = False
        hard_gate_status = "unresolved"
    else:
        verdict = "pass"
        hard_gate_applicable = True
        hard_gate_status = "passed"

    desired_ranking_gate_applicable = bool(
        desired_metric_applicable
        and desired_score["ranking_score"] is not None
    )
    desired_coverage_protocol_gate_applicable = bool(
        desired_metric_applicable
        and desired_score["evaluable_coverage"] is not None
    )

    return {
        "schema_version": schema_version,
        "blueprint_valid": bool(blueprint_valid),
        "blueprint_applicable": blueprint_applicable,
        "hard_gate_declared": declared_hard_gate,
        "hard_gate_applicable": hard_gate_applicable,
        # ``hard_gate_*`` is retained for schema compatibility and refers only
        # to the curated Required-group verdict.  Desired has its own numeric
        # 0.10/0.50 gate in the metric table, so a valid desired-only blueprint
        # must not be described as if every Blueprint gate were disabled.
        "required_hard_gate_declared": declared_hard_gate,
        "required_hard_gate_applicable": hard_gate_applicable,
        "desired_metric_applicable": desired_metric_applicable,
        "desired_ranking_gate_applicable": desired_ranking_gate_applicable,
        "desired_coverage_protocol_gate_applicable": (
            desired_coverage_protocol_gate_applicable
        ),
        "molecule_hard_gate_applicable": bool(
            hard_gate_applicable or desired_ranking_gate_applicable
        ),
        "hard_gate_status": hard_gate_status,
        "required_logic": "all_groups",
        "required_group_results": group_results,
        "required_hit_rate": required_hit_rate,
        "required_hit_rate_display": (
            round(required_hit_rate, 3)
            if required_hit_rate is not None else None
        ),
        "required_hits": len(required_hits),
        "required_total": required_total,
        "required_missed": [r["constraint"] for r in required_misses],
        "required_missed_labels": [_entry_label(r["constraint"]) for r in required_misses],
        "required_unresolved": [r["constraint"] for r in required_unresolved],
        "required_unresolved_labels": [_entry_label(r["constraint"]) for r in required_unresolved],
        "required_details": required_results,
        "desired_hits": len(desired_hits),
        "desired_total": len(desired),
        "desired_hit_rate": (
            desired_score["hit_rate"] if desired_metric_applicable else None
        ),
        "desired_hit_rate_display": (
            desired_score["hit_rate_display"]
            if desired_metric_applicable else None
        ),
        "desired_weight_hit": desired_score["weight_hit"],
        "desired_weight_hit_display": desired_score["weight_hit_display"],
        "desired_weight_total": desired_score["weight_total"],
        "desired_weight_total_display": desired_score["weight_total_display"],
        "desired_weighted_score": (
            desired_score["weighted_score"] if desired_metric_applicable else None
        ),
        "desired_weighted_score_display": (
            desired_score["weighted_score_display"]
            if desired_metric_applicable else None
        ),
        "desired_evaluable_count": desired_score["evaluable_count"],
        "desired_unmatchable_count": desired_score["unmatchable_count"],
        "desired_evaluable_coverage": (
            desired_score["evaluable_coverage"]
            if desired_metric_applicable else None
        ),
        "desired_evaluable_coverage_display": (
            desired_score["evaluable_coverage_display"]
            if desired_metric_applicable else None
        ),
        "desired_ranking_score": (
            desired_score["ranking_score"] if desired_metric_applicable else None
        ),
        "desired_ranking_score_display": (
            desired_score["ranking_score_display"]
            if desired_metric_applicable else None
        ),
        "desired_missed": [r["constraint"] for r in desired_misses],
        "desired_unresolved": [r["constraint"] for r in desired_unresolved],
        "desired_details": desired_results,
        "optional_metric_applicable": optional_metric_applicable,
        "optional_hits": (
            len(optional_hits) if optional_metric_applicable else None
        ),
        "optional_total": len(optional),
        "optional_details": optional_results,
        "evaluator_capabilities": _json_capabilities(capabilities),
        "verdict": verdict,
    }


def unresolved_match_result(blueprint: dict, reason: str) -> dict:
    """Build a per-pose result for unavailable interaction evidence."""
    required = [item for item in blueprint.get("required", []) if isinstance(item, dict)]
    groups = _required_groups(blueprint, required)
    grouped_required = {
        _constraint_key(member): member
        for group in groups
        for member in group["members"]
    }
    validation = blueprint.get("validation") or {}
    declared = _declared_hard_gate(validation, bool(groups))
    blueprint_valid = bool(validation.get("blueprint_valid", True))
    blueprint_applicable = not bool(
        validation.get("not_applicable") is True
        or validation.get("applicable") is False
        or blueprint.get("not_applicable") is True
        or blueprint.get("applicable") is False
    )
    desired = [
        item for item in blueprint.get("desired", []) if isinstance(item, dict)
    ]
    declared_desired_metric = validation.get("desired_metric_applicable")
    desired_metric_applicable = bool(
        blueprint_applicable
        and blueprint_valid
        and desired
        and (
            bool(declared_desired_metric)
            if declared_desired_metric is not None else True
        )
    )
    return {
        "schema_version": _schema_version(blueprint),
        "blueprint_valid": blueprint_valid,
        "blueprint_applicable": blueprint_applicable,
        "hard_gate_declared": declared,
        "hard_gate_applicable": False,
        "required_hard_gate_declared": declared,
        "required_hard_gate_applicable": False,
        "desired_metric_applicable": desired_metric_applicable,
        "desired_ranking_gate_applicable": False,
        "desired_coverage_protocol_gate_applicable": desired_metric_applicable,
        "molecule_hard_gate_applicable": False,
        "hard_gate_status": "unresolved" if declared else "not_applicable",
        "required_logic": "all_groups",
        "required_group_results": [],
        "required_hit_rate": None,
        "required_hit_rate_display": None,
        "required_hits": 0,
        "required_total": len(grouped_required),
        "required_missed": [],
        "required_missed_labels": [],
        "required_unresolved": list(grouped_required.values()),
        "required_unresolved_labels": [
            _entry_label(item) for item in grouped_required.values()
        ],
        "desired_hits": 0,
        "desired_total": len(desired),
        "desired_hit_rate": None,
        "desired_hit_rate_display": None,
        "desired_weight_hit": 0.0,
        "desired_weight_hit_display": 0.0,
        "desired_weight_total": 0.0,
        "desired_weight_total_display": 0.0,
        "desired_weighted_score": None,
        "desired_weighted_score_display": None,
        "desired_evaluable_count": 0,
        "desired_unmatchable_count": len(desired),
        "desired_evaluable_coverage": (
            0.0 if desired_metric_applicable else None
        ),
        "desired_evaluable_coverage_display": (
            0.0 if desired_metric_applicable else None
        ),
        "desired_ranking_score": None,
        "desired_ranking_score_display": None,
        "optional_hits": None,
        "optional_metric_applicable": False,
        "optional_total": len(blueprint.get("optional", [])),
        "verdict": "unresolved" if declared else "not_applicable",
        "reason": reason,
    }


def _schema_version(blueprint: dict) -> int:
    value = blueprint.get("schema_version", 1)
    try:
        return int(str(value).split(".", 1)[0])
    except (TypeError, ValueError):
        return 1


def _declared_hard_gate(validation: dict, has_required_groups: bool) -> bool:
    required = validation.get("required_hard_gate_applicable")
    if required is not None:
        return bool(required) and has_required_groups
    explicit = validation.get("hard_gate_applicable")
    if explicit is not None:
        return bool(explicit) and has_required_groups
    legacy = validation.get("valid_for_hard_gating")
    if legacy is not None:
        return bool(legacy) and has_required_groups
    return has_required_groups


def _required_groups(blueprint: dict, required: list[dict]) -> list[dict]:
    """Resolve inline or ID-referenced required logic groups.

    Supported locations and aliases intentionally make the matcher tolerant of
    early schema-v2 files while retaining a single canonical output shape.
    """
    constraints = blueprint.get("constraints") or {}
    raw_groups = (
        blueprint.get("required_groups")
        or constraints.get("required_groups")
        or constraints.get("logic_groups")
        or []
    )
    by_id = {
        str(item.get("constraint_id") or item.get("id")): item
        for item in required
        if item.get("constraint_id") or item.get("id")
    }
    groups: list[dict] = []
    referenced: set[str] = set()
    for index, raw in enumerate(raw_groups):
        if not isinstance(raw, dict):
            continue
        members_raw = (
            raw.get("members")
            or raw.get("constraints")
            or raw.get("interactions")
            or raw.get("member_ids")
            or raw.get("constraint_ids")
            or raw.get("member_constraint_ids")
            or []
        )
        members = []
        for member in members_raw:
            if isinstance(member, dict):
                members.append(member)
            elif str(member) in by_id:
                members.append(by_id[str(member)])
        if not members:
            continue
        logic = str(raw.get("logic") or raw.get("operator") or "all_of").lower()
        logic = {"all": "all_of", "any": "any_of", "at_least_n": "at_least_k"}.get(
            logic, logic,
        )
        if logic not in {"all_of", "any_of", "at_least_k"}:
            logic = "all_of"
        k = raw.get(
            "at_least_k",
            raw.get("at_least_n", raw.get("k", raw.get("min_required"))),
        )
        if logic == "all_of":
            k = len(members)
        elif logic == "any_of":
            k = 1
        else:
            try:
                k = int(k)
            except (TypeError, ValueError):
                k = 1
            k = max(1, min(k, len(members)))
        groups.append({
            "group_id": str(raw.get("group_id") or raw.get("id") or f"required_group_{index + 1}"),
            "logic": logic,
            "at_least_k": k,
            "members": members,
        })
        referenced.update(_constraint_key(member) for member in members)

    # A required item omitted from explicit groups remains an independent hard
    # all-of constraint.  This prevents an accidental group omission weakening
    # the legacy required list.
    for index, item in enumerate(required):
        key = _constraint_key(item)
        if key in referenced:
            continue
        groups.append({
            "group_id": str(item.get("group_id") or f"required_{index + 1}"),
            "logic": "all_of",
            "at_least_k": 1,
            "members": [item],
        })
    return groups


def _evaluate_group(group: dict, members: list[dict]) -> dict:
    hits = sum(item["status"] == "matched" for item in members)
    unresolved = sum(item["status"] == "unmatchable" for item in members)
    missing = sum(item["status"] == "missing" for item in members)
    k = int(group["at_least_k"])
    if hits >= k:
        status = "passed"
    elif hits + unresolved >= k:
        status = "unresolved"
    else:
        status = "failed"
    return {
        "group_id": group["group_id"],
        "logic": group["logic"],
        "at_least_k": k,
        "status": status,
        "matched": hits,
        "missing": missing,
        "unmatchable": unresolved,
        "total": len(members),
        "member_results": members,
    }


def _evaluate_constraint(
    constraint: dict,
    pose_interactions: list[dict],
    capabilities: dict,
    *,
    strict_details: bool,
) -> dict:
    expected_type = _normalize_type(constraint.get("interaction_type", constraint.get("type", "")))
    if not _capability_supports_type(capabilities, expected_type):
        return _constraint_result(
            constraint, "unmatchable", "interaction_type_not_supported",
            unsupported_fields=["interaction_type"],
        )

    expected_key = _normalize_blueprint_entry(constraint)
    candidates = [
        item for item in pose_interactions
        if _base_interaction_matches(expected_key, item)
    ]
    if not candidates:
        if _icode_unobservable(constraint, pose_interactions, expected_type):
            return _constraint_result(
                constraint, "unmatchable", "residue_insertion_code_not_observable",
                unsupported_fields=["residue.icode"],
            )
        return _constraint_result(constraint, "missing", "no_residue_type_match")

    if not strict_details or not any(constraint.get(field) not in (None, "", [], {}) for field in _DETAIL_FIELDS):
        return _constraint_result(constraint, "matched", "residue_type_match", matched=candidates[0])

    mismatch_comparisons: list[dict] = []
    unsupported: set[str] = set()
    for candidate in candidates:
        comparison = _compare_details(constraint, candidate, capabilities, expected_type)
        if comparison["status"] == "matched":
            return _constraint_result(constraint, "matched", "strict_match", matched=candidate)
        if comparison["status"] == "missing":
            mismatch_comparisons.append(comparison)
        else:
            unsupported.update(comparison["unsupported_fields"])

    # A fully comparable candidate with a value mismatch is affirmative negative
    # evidence. Report the closest candidate rather than unioning unrelated
    # mismatches from every interaction at the same residue.
    if mismatch_comparisons:
        closest = min(
            mismatch_comparisons,
            key=lambda item: (
                len(set(item["mismatches"])),
                len(set(item["unsupported_fields"])),
                tuple(sorted(set(item["mismatches"]))),
                tuple(sorted(set(item["unsupported_fields"]))),
            ),
        )
        return _constraint_result(
            constraint, "missing", "detail_mismatch",
            mismatches=sorted(set(closest["mismatches"])),
        )
    return _constraint_result(
        constraint, "unmatchable", "required_detail_not_observable",
        unsupported_fields=sorted(unsupported),
    )


def _compare_details(
    expected: dict,
    actual: dict,
    capabilities: dict,
    interaction_type: str,
) -> dict:
    mismatches: list[str] = []
    unsupported: list[str] = []

    for field in ("ligand_role", "protein_role", "mediator"):
        expected_value = expected.get(field)
        if expected_value in (None, "", [], {}):
            continue
        actual_value = _field_value(actual, field)
        if actual_value in (None, "", [], {}):
            unsupported.append(field)
        elif not _values_match(expected_value, actual_value):
            mismatches.append(field)

    for field in ("protein_atoms", "ligand_atoms"):
        expected_value = expected.get(field)
        if expected_value in (None, "", [], {}):
            continue
        actual_value = _field_value(actual, field)
        if actual_value in (None, "", [], {}):
            unsupported.append(field)
        elif field == "protein_atoms" and not _protein_atoms_match(
            expected_value, actual_value, expected, actual,
        ):
            mismatches.append(field)
        elif field == "ligand_atoms" and not _ligand_atoms_match(
            expected_value, actual_value,
        ):
            mismatches.append(field)

    geometry = dict(expected.get("geometry") or {})
    if expected.get("distance_range") not in (None, "", [], {}):
        geometry["distance_range"] = expected["distance_range"]
    if expected.get("angle_range") not in (None, "", [], {}):
        geometry["angle_range"] = expected["angle_range"]
    if geometry not in (None, "", [], {}):
        geometry_result = _compare_geometry(geometry, actual, capabilities, interaction_type)
        mismatches.extend(geometry_result["mismatches"])
        unsupported.extend(geometry_result["unsupported_fields"])

    if mismatches:
        return {"status": "missing", "mismatches": mismatches, "unsupported_fields": unsupported}
    if unsupported:
        return {"status": "unmatchable", "mismatches": [], "unsupported_fields": unsupported}
    return {"status": "matched", "mismatches": [], "unsupported_fields": []}


def _compare_geometry(
    expected_geometry: Any,
    actual: dict,
    capabilities: dict,
    interaction_type: str,
) -> dict:
    expected = expected_geometry if isinstance(expected_geometry, dict) else {"distance": expected_geometry}
    actual_geometry = actual.get("geometry") if isinstance(actual.get("geometry"), dict) else {}
    mismatches: list[str] = []
    unsupported: list[str] = []
    measures = _geometry_measures(expected)

    for metadata_field in ("distance_metric", "distance_unit", "angle_unit"):
        expected_value = expected.get(metadata_field)
        if expected_value in (None, ""):
            continue
        if metadata_field == "angle_unit" and not any("angle" in name for name in measures):
            continue
        if metadata_field == "distance_unit" and not any("angle" not in name for name in measures):
            continue
        actual_value = actual_geometry.get(metadata_field, actual.get(metadata_field))
        field_path = f"geometry.{metadata_field}"
        if actual_value in (None, "") or not _capability_supports_field(
            capabilities, interaction_type, field_path,
        ):
            unsupported.append(field_path)
        elif _normal_text(expected_value) != _normal_text(actual_value):
            # A metric/unit mismatch makes the values incomparable rather than
            # proving that the interaction is absent.
            unsupported.append(field_path)

    for measure in measures:
        tolerance = (
            BLUEPRINT_ANGLE_TOLERANCE_DEGREES
            if "angle" in measure
            else BLUEPRINT_DISTANCE_TOLERANCE_ANGSTROM
        )
        specification = _geometry_specification(expected, measure)
        if specification is None:
            continue
        actual_value = actual_geometry.get(measure, actual.get(measure))
        if actual_value is None or not _capability_supports_field(capabilities, interaction_type, f"geometry.{measure}"):
            unsupported.append(f"geometry.{measure}")
            continue
        if not _geometry_matches(specification, actual_value, tolerance):
            mismatches.append(f"geometry.{measure}")
    return {"mismatches": mismatches, "unsupported_fields": unsupported}


def _geometry_measures(geometry: dict) -> list[str]:
    """Return every numeric/range measure instead of silently ignoring new keys."""
    metadata = {"distance_metric", "metric", "distance_unit", "angle_unit"}
    measures = set()
    for key in geometry:
        if key in metadata:
            continue
        if key.startswith("min_"):
            measures.add(key[4:])
        elif key.startswith("max_"):
            measures.add(key[4:])
        elif key.endswith("_min"):
            measures.add(key[:-4])
        elif key.endswith("_max"):
            measures.add(key[:-4])
        elif key.endswith("_range"):
            measures.add(key[:-6])
        else:
            measures.add(key)
    return sorted(measures)


def _geometry_specification(geometry: dict, measure: str) -> Any:
    if measure in geometry:
        return geometry[measure]
    low = geometry.get(f"{measure}_min", geometry.get(f"min_{measure}"))
    high = geometry.get(f"{measure}_max", geometry.get(f"max_{measure}"))
    if low is not None or high is not None:
        return {"min": low, "max": high}
    if f"{measure}_range" in geometry:
        return geometry[f"{measure}_range"]
    return None


def _geometry_matches(specification: Any, actual_value: Any, tolerance: float) -> bool:
    try:
        actual = float(actual_value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(actual):
        return False
    if isinstance(specification, dict):
        low = specification.get("min")
        high = specification.get("max")
        target = specification.get("value")
        local_tolerance = specification.get("tolerance", tolerance)
        if low is not None and actual < float(low):
            return False
        if high is not None and actual > float(high):
            return False
        if target is not None and abs(actual - float(target)) > float(local_tolerance):
            return False
        return True
    if isinstance(specification, (list, tuple)) and len(specification) == 2:
        return float(specification[0]) <= actual <= float(specification[1])
    try:
        return abs(actual - float(specification)) <= tolerance
    except (TypeError, ValueError):
        return False


def _desired_score(results: list[dict]) -> dict:
    evaluable = [item for item in results if item["status"] != "unmatchable"]
    all_weight = sum(_constraint_weight(item["constraint"]) for item in results)
    weight_total = sum(_constraint_weight(item["constraint"]) for item in evaluable)
    weight_hit = sum(
        _constraint_weight(item["constraint"])
        for item in evaluable
        if item["status"] == "matched"
    )
    hit_rate = (
        sum(item["status"] == "matched" for item in evaluable) / len(evaluable)
        if evaluable else None
    )
    coverage = weight_total / all_weight if all_weight else None
    weighted_score = weight_hit / weight_total if weight_total else None
    ranking_score = weight_hit / all_weight if all_weight else None
    return {
        "weight_hit": weight_hit,
        "weight_hit_display": round(weight_hit, 3),
        "weight_total": weight_total,
        "weight_total_display": round(weight_total, 3),
        "weighted_score": weighted_score,
        "weighted_score_display": (
            round(weighted_score, 3) if weighted_score is not None else None
        ),
        "hit_rate": hit_rate,
        "hit_rate_display": round(hit_rate, 3) if hit_rate is not None else None,
        "evaluable_count": len(evaluable),
        "unmatchable_count": len(results) - len(evaluable),
        "evaluable_coverage": coverage,
        "evaluable_coverage_display": (
            round(coverage, 3) if coverage is not None else None
        ),
        # Ranking must account for profiler coverage.  Keep the conditional
        # weighted score for interpretation, but never let 1/1 evaluable look
        # equivalent to 10/10 when nine constraints are unobservable.
        "ranking_score": ranking_score,
        "ranking_score_display": (
            round(ranking_score, 3) if ranking_score is not None else None
        ),
    }


def _constraint_weight(constraint: dict) -> float:
    if constraint.get("weight") is not None:
        try:
            return max(0.0, float(constraint["weight"]))
        except (TypeError, ValueError):
            pass
    confidence = str(constraint.get("confidence") or "").lower()
    return BLUEPRINT_CONFIDENCE_WEIGHTS.get(confidence, BLUEPRINT_DESIRED_DEFAULT_WEIGHT)


def _constraint_result(
    constraint: dict,
    status: str,
    reason: str,
    *,
    matched: dict | None = None,
    mismatches: list[str] | None = None,
    unsupported_fields: list[str] | None = None,
) -> dict:
    result = {
        "constraint_id": constraint.get("constraint_id") or constraint.get("id"),
        "label": _entry_label(constraint),
        "status": status,
        "reason": reason,
        "constraint": constraint,
    }
    if matched is not None:
        result["matched_interaction"] = matched
    if mismatches:
        result["mismatches"] = mismatches
    if unsupported_fields:
        result["unsupported_fields"] = unsupported_fields
    return result


def _normalize_capabilities(capabilities: dict | None) -> dict:
    source = capabilities or DEFAULT_PLIP_CAPABILITIES
    return {
        "source": source.get("source", "unknown"),
        "interaction_types": {
            _normalize_type(value) for value in source.get("interaction_types", [])
        },
        "fields": set(source.get("fields", [])),
        "fields_by_type": {
            _normalize_type(key): set(value)
            for key, value in (source.get("fields_by_type") or {}).items()
        },
    }


def _json_capabilities(capabilities: dict) -> dict:
    return {
        "source": capabilities["source"],
        "interaction_types": sorted(capabilities["interaction_types"]),
        "fields": sorted(capabilities["fields"]),
        "fields_by_type": {
            key: sorted(value) for key, value in sorted(capabilities["fields_by_type"].items())
        },
    }


def _capability_supports_type(capabilities: dict, interaction_type: str) -> bool:
    supported = capabilities.get("interaction_types") or set()
    return not supported or interaction_type in supported


def _capability_supports_field(capabilities: dict, interaction_type: str, field: str) -> bool:
    if field in capabilities.get("fields", set()):
        return True
    return field in capabilities.get("fields_by_type", {}).get(interaction_type, set())


def _base_interaction_matches(expected_key: tuple[str, str, str], actual: dict) -> bool:
    chain, residue, interaction_type = expected_key
    actual_chain, actual_residue = _normalize_residue_fields(
        actual.get("residue", ""), actual.get("chain", ""), actual.get("icode", ""),
    )
    if chain and chain != actual_chain:
        return False
    return residue == actual_residue and interaction_type == _normalize_type(
        actual.get("interaction_type", actual.get("type", ""))
    )


def _icode_unobservable(
    constraint: dict, pose_interactions: list[dict], expected_type: str,
) -> bool:
    """Distinguish a missing iCode capability from an affirmative mismatch."""
    expected = parse_residue_record(
        constraint.get("residue", ""),
        default_chain=str(constraint.get("chain") or "").strip(),
    )
    if expected is None:
        return False
    expected_icode = str(constraint.get("icode") or expected.get("icode") or "").strip()
    if not expected_icode:
        return False
    expected_chain = str(expected.get("chain") or "").strip()
    expected_name = str(expected.get("resname") or "").upper()
    expected_number = str(expected.get("resnum") or "")
    found_without_icode = False
    for actual in pose_interactions:
        if _normalize_type(actual.get("interaction_type", actual.get("type", ""))) != expected_type:
            continue
        record = parse_residue_record(
            actual.get("residue", ""),
            default_chain=str(actual.get("chain") or "").strip(),
        )
        if record is None:
            continue
        if actual.get("icode") not in (None, ""):
            record = make_residue_record(
                record.get("resname", ""), record.get("resnum", ""),
                record.get("chain", ""), actual.get("icode", ""),
            )
        if (
            str(record.get("chain") or "").strip() == expected_chain
            and str(record.get("resname") or "").upper() == expected_name
            and str(record.get("resnum") or "") == expected_number
        ):
            actual_icode = str(record.get("icode") or "").strip()
            if actual_icode:
                return False
            found_without_icode = True
    return found_without_icode


def _field_value(interaction: dict, field: str) -> Any:
    aliases = {
        "protein_atoms": ("protein_atoms", "protein_atom", "protein_atom_names"),
        "ligand_atoms": ("ligand_atoms", "ligand_atom", "ligand_atom_indices", "ligand_atom_names"),
        "ligand_role": ("ligand_role",),
        "protein_role": ("protein_role",),
        "mediator": ("mediator", "interaction_mediator"),
    }
    for alias in aliases.get(field, (field,)):
        if interaction.get(alias) not in (None, "", [], {}):
            return interaction[alias]
    return None


def _values_match(expected: Any, actual: Any) -> bool:
    if isinstance(expected, (list, tuple, set)) or isinstance(actual, (list, tuple, set)):
        expected_values = {_normal_text(value) for value in _as_list(expected)}
        actual_values = {_normal_text(value) for value in _as_list(actual)}
        return expected_values == actual_values
    return _normal_text(expected) == _normal_text(actual)


def _protein_atoms_match(
    expected_value: Any,
    actual_value: Any,
    expected_interaction: dict,
    actual_interaction: dict,
) -> bool:
    """Match transferable protein atom constraints against profiler evidence.

    Blueprint protein atoms are receptor atom constraints, usually stable atom
    names such as ``N`` or ``NZ``.  A profiler may include additional atoms that
    participate in the same feature (for example the donor hydrogen), so every
    expected atom must be present but the actual set need not be exactly equal.
    Tool-local indices and explicit canonical IDs remain strict identities.
    """
    expected_atoms = _as_list(expected_value)
    if not all(
        _protein_atom_matches_interaction_residue(atom, expected_interaction)
        for atom in expected_atoms
    ):
        return False
    actual_atoms = [
        atom for atom in _as_list(actual_value)
        if _protein_atom_matches_interaction_residue(atom, actual_interaction)
    ]
    return bool(actual_atoms) and all(
        any(_atom_identity_matches(expected, actual) for actual in actual_atoms)
        for expected in expected_atoms
    )


def _ligand_atoms_match(expected_value: Any, actual_value: Any) -> bool:
    """Match ligand atom records exactly while honoring every declared identity."""
    expected_atoms = _as_list(expected_value)
    actual_atoms = _as_list(actual_value)
    if len(expected_atoms) != len(actual_atoms):
        return False
    if not all(
        _atom_record_has_consistent_canonical_name(atom)
        for atom in [*expected_atoms, *actual_atoms]
    ):
        return False
    return _atom_lists_have_exact_matching(expected_atoms, actual_atoms)


def _atom_lists_have_exact_matching(expected_atoms: list[Any], actual_atoms: list[Any]) -> bool:
    if not expected_atoms:
        return not actual_atoms
    expected = expected_atoms[0]
    for index, actual in enumerate(actual_atoms):
        if not _atom_identity_matches(expected, actual):
            continue
        if _atom_lists_have_exact_matching(
            expected_atoms[1:], actual_atoms[:index] + actual_atoms[index + 1:],
        ):
            return True
    return False


def _atom_identity_matches(expected: Any, actual: Any) -> bool:
    if isinstance(expected, dict):
        checks: list[bool] = []
        canonical_id = expected.get("canonical_id")
        if canonical_id not in (None, ""):
            checks.append(
                isinstance(actual, dict)
                and _normal_text(canonical_id) == _normal_text(actual.get("canonical_id"))
            )

        expected_name = expected.get("name") or expected.get("atom_name")
        if expected_name not in (None, ""):
            checks.append(_normal_text(expected_name) == _atom_name(actual))

        expected_id = expected.get("id")
        if expected_id not in (None, ""):
            checks.append(
                isinstance(actual, dict)
                and _normal_text(expected_id) == _normal_text(actual.get("id"))
            )

        has_explicit_identity = any(
            value not in (None, "")
            for value in (canonical_id, expected_name, expected_id, expected.get("index"))
        )
        if has_explicit_identity and expected.get("namespace") not in (None, ""):
            checks.append(
                isinstance(actual, dict)
                and _normal_text(expected["namespace"])
                == _normal_text(actual.get("namespace"))
            )

        if expected.get("index") is not None:
            actual_index = actual.get("index") if isinstance(actual, dict) else None
            expected_namespace = _normal_text(expected.get("namespace")) or "unspecified"
            actual_namespace = (
                _normal_text(actual.get("namespace")) or "unspecified"
                if isinstance(actual, dict)
                else "unspecified"
            )
            checks.append(
                actual_index is not None
                and expected_namespace == actual_namespace
                and _normal_text(expected["index"]) == _normal_text(actual_index)
            )

        if checks:
            return all(checks)

    if isinstance(expected, str):
        return _normal_text(expected) == _atom_name(actual)
    return _canonical_atom_set([expected]) == _canonical_atom_set([actual])


def _atom_name(atom: Any) -> str:
    if not isinstance(atom, dict):
        return _normal_text(atom) if isinstance(atom, str) else ""
    name = atom.get("name") or atom.get("atom_name")
    if name not in (None, ""):
        return _normal_text(name)
    canonical_id = atom.get("canonical_id")
    if canonical_id not in (None, "") and ":" in str(canonical_id):
        return _normal_text(str(canonical_id).rsplit(":", 1)[-1])
    return ""


def _protein_atom_matches_interaction_residue(atom: Any, interaction: dict) -> bool:
    """Reject atom metadata that contradicts its interaction residue."""
    if not isinstance(atom, dict):
        return True
    if not _atom_record_has_consistent_canonical_name(atom):
        return False
    context = parse_residue_record(
        interaction.get("residue", ""),
        default_chain=str(interaction.get("chain") or "").strip(),
    )
    if context is None:
        return True
    if interaction.get("icode") not in (None, ""):
        context = make_residue_record(
            context.get("resname", ""), context.get("resnum", ""),
            context.get("chain", ""), interaction["icode"],
        )

    explicit_fields = {
        "chain": atom.get("chain", atom.get("chain_id")),
        "resname": atom.get("resname", atom.get("res_name")),
        "resnum": atom.get("resnum", atom.get("res_seq")),
        "icode": atom.get("icode", atom.get("i_code")),
    }
    for field, value in explicit_fields.items():
        if value is None or (value == "" and field != "icode"):
            continue
        if not _residue_field_matches(field, value, context.get(field)):
            return False

    canonical_id = atom.get("canonical_id")
    if canonical_id not in (None, ""):
        canonical_context = _residue_from_atom_canonical_id(canonical_id)
        if canonical_context is None:
            return False
        for field in ("chain", "resname", "resnum", "icode"):
            value = canonical_context.get(field)
            if value in (None, "") and field != "icode":
                continue
            if not _residue_field_matches(field, value, context.get(field)):
                return False
    return True


def _atom_record_has_consistent_canonical_name(atom: Any) -> bool:
    if not isinstance(atom, dict) or atom.get("canonical_id") in (None, ""):
        return True
    canonical_name = _canonical_atom_name(atom["canonical_id"])
    if not canonical_name:
        return False
    explicit_name = atom.get("name") or atom.get("atom_name")
    return (
        explicit_name in (None, "")
        or _normal_text(explicit_name) == canonical_name
    )


def _canonical_atom_name(value: Any) -> str:
    if value in (None, "") or ":" not in str(value):
        return ""
    return _normal_text(str(value).rsplit(":", 1)[-1])


def _residue_from_atom_canonical_id(value: Any) -> dict | None:
    if value in (None, "") or ":" not in str(value):
        return None
    residue_token = str(value).rsplit(":", 1)[0]
    if ":" in residue_token:
        chain, residue = residue_token.split(":", 1)
        numbered = re.fullmatch(r"(-?\d+)([A-Za-z]?)", residue.strip())
        if numbered:
            return make_residue_record(
                "", numbered.group(1), chain, numbered.group(2),
            )
        return parse_residue_record(residue, default_chain=chain)
    return parse_residue_record(residue_token)


def _residue_field_matches(field: str, first: Any, second: Any) -> bool:
    if field == "resnum":
        try:
            return int(first) == int(second)
        except (TypeError, ValueError):
            pass
    return _normal_text(first) == _normal_text(second)


def _canonical_atom_set(value: Any) -> set[str]:
    atoms = _as_list(value)
    normalized = set()
    for atom in atoms:
        if isinstance(atom, dict):
            namespace = _normal_text(atom.get("namespace")) or "unspecified"
            if atom.get("canonical_id") not in (None, ""):
                atom_value = f"canonical:{_normal_text(atom['canonical_id'])}"
            elif atom.get("name") or atom.get("atom_name"):
                atom_value = f"name:{_normal_text(atom.get('name') or atom.get('atom_name'))}"
            elif atom.get("index") is not None:
                atom_value = f"{namespace}:index:{_normal_text(atom['index'])}"
            else:
                atom_value = f"{namespace}:id:{_normal_text(atom.get('id'))}"
        else:
            atom_value = atom
        normalized.add(_normal_text(atom_value))
    return normalized


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def _normal_text(value: Any) -> str:
    raw = "" if value is None else str(value)
    return re.sub(r"[\s_-]+", "", raw.strip().lower())


def _constraint_key(entry: dict) -> str:
    explicit = entry.get("constraint_id") or entry.get("id")
    if explicit:
        return f"id:{explicit}"
    return json.dumps(entry, sort_keys=True, default=str, separators=(",", ":"))


def _normalize_blueprint_entry(entry: dict) -> tuple[str, str, str]:
    chain, residue = _normalize_residue_fields(
        entry.get("residue", ""), entry.get("chain", ""), entry.get("icode", ""),
    )
    interaction_type = _normalize_type(entry.get("interaction_type", entry.get("type", "")))
    return chain, residue, interaction_type


def _normalize_residue(res: str) -> str:
    """Normalize residue ID while retaining an insertion code when present."""
    _, residue = _normalize_residue_fields(res)
    return residue


def _normalize_residue_fields(
    residue: str, chain: str = "", icode: str = "",
) -> tuple[str, str]:
    """Return chain plus a legacy label that still contains the insertion code."""
    record = parse_residue_record(residue, default_chain=str(chain or "").strip())
    if record is None:
        return str(chain or _extract_chain(residue)).strip(), _shared_normalize_residue(residue)
    if icode not in (None, ""):
        record = make_residue_record(
            record.get("resname", ""), record.get("resnum", ""),
            record.get("chain", ""), icode,
        )
    return str(record.get("chain", "")).strip(), legacy_residue_label(record)


def _extract_chain(res: str) -> str:
    if "." in str(res):
        parts = str(res).split(".")
        return parts[1].strip() if len(parts) > 1 else ""
    if "_" in str(res):
        return str(res).split("_", 1)[0].strip()
    return ""


def _normalize_type(interaction_type: str) -> str:
    return normalize_interaction_type(str(interaction_type or ""))


def _entry_label(entry: dict) -> str:
    chain, residue = _normalize_residue_fields(
        entry.get("residue", ""), entry.get("chain", ""), entry.get("icode", ""),
    )
    residue = residue or "?"
    if chain:
        residue = f"{chain}:{residue}"
    return f"{residue}:{entry.get('interaction_type', entry.get('type', '?'))}"

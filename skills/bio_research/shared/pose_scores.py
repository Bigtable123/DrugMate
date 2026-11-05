"""Canonical pose-score records and engine-neutral primary-score selection.

The pipeline historically stored scorer values as flat engine-specific keys.
This module keeps those files readable while exposing a single JSON-friendly
contract to Step 4 and later consumers.  It deliberately selects one comparable
metric for a ranking cohort; it never combines scores with different units.
"""
from __future__ import annotations

import math
from collections import Counter
from typing import Any, Iterable


SCHEMA_VERSION = 1

# Flat fields emitted by the existing docking adapters.  The canonical metric
# id includes the provider so similarly named values from different methods do
# not silently become comparable.
LEGACY_SCORE_SPECS: dict[str, dict[str, str]] = {
    "minimizedAffinity": {
        "metric_id": "gnina.minimizedAffinity",
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "gnina",
    },
    "CNNscore": {
        "metric_id": "gnina.CNNscore",
        "unit": "probability-like score",
        "direction": "higher",
        "method": "gnina",
    },
    "CNNaffinity": {
        "metric_id": "gnina.CNNaffinity",
        "unit": "pK",
        "direction": "higher",
        "method": "gnina",
    },
    "intramol": {
        "metric_id": "gnina.intramol",
        "unit": "GNINA score unit",
        "direction": "context_only",
        "method": "gnina",
    },
    "CNNvariance": {
        "metric_id": "gnina.CNNvariance",
        "unit": "ensemble variance",
        "direction": "context_only",
        "method": "gnina",
    },
    "docking_score": {
        "metric_id": "unidock.docking_score",
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "unidock",
    },
    "vina_affinity": {
        "metric_id": "autodock_vina.vina_affinity",
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "autodock_vina",
    },
    "vina_inter": {
        "metric_id": "autodock_vina.vina_inter",
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "autodock_vina",
    },
    "vina_intra": {
        "metric_id": "autodock_vina.vina_intra",
        "unit": "kcal/mol",
        "direction": "context_only",
        "method": "autodock_vina",
    },
    "vina_torsions": {
        "metric_id": "autodock_vina.vina_torsions",
        "unit": "kcal/mol",
        "direction": "context_only",
        "method": "autodock_vina",
    },
    "vina_intra_best_pose": {
        "metric_id": "autodock_vina.vina_intra_best_pose",
        "unit": "kcal/mol",
        "direction": "context_only",
        "method": "autodock_vina",
    },
    "diffdock_confidence": {
        "metric_id": "diffdock.confidence",
        "unit": "DiffDock confidence",
        "direction": "higher",
        "method": "diffdock",
    },
    "dock6_grid_score": {
        "metric_id": "dock6.grid_score",
        "unit": "DOCK 6 grid score",
        "direction": "lower",
        "method": "dock6",
    },
    "rdock_score": {
        "metric_id": "rdock.score",
        "unit": "rDock score",
        "direction": "lower",
        "method": "rdock",
    },
    "autodock4_binding_energy": {
        "metric_id": "autodock4.binding_energy",
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "autodock4",
    },
    "adfr_affinity": {
        "metric_id": "adfr.affinity",
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "adfr",
    },
    "ligand_efficiency": {
        "metric_id": "autodock_vina.ligand_efficiency",
        "unit": "kcal/mol/heavy atom",
        "direction": "higher",
        "method": "autodock_vina",
    },
    "score_gap_1_2": {
        "metric_id": "autodock_vina.score_gap_1_2",
        "unit": "kcal/mol",
        "direction": "context_only",
        "method": "autodock_vina",
    },
    "engine_pose_rank": {
        "metric_id": "autodock_vina.engine_pose_rank",
        "unit": "rank",
        "direction": "context_only",
        "method": "autodock_vina",
    },
}

POSE_SELECTION_PREFERENCE = (
    "gnina.CNNscore",
    "gnina.CNNaffinity",
    "gnina.minimizedAffinity",
    "unidock.docking_score",
    "autodock_vina.vina_affinity",
    "autodock4.binding_energy",
    "rdock.score",
    "adfr.affinity",
    "diffdock.confidence",
    "dock6.grid_score",
)

CANDIDATE_RANKING_PREFERENCE = (
    "gnina.CNNaffinity",
    "gnina.CNNscore",
    "gnina.minimizedAffinity",
    "unidock.docking_score",
    "autodock_vina.vina_affinity",
    "autodock4.binding_energy",
    "rdock.score",
    "adfr.affinity",
    "diffdock.confidence",
    "dock6.grid_score",
)

_UNUSABLE_RELIABILITY = {
    "failed",
    "invalid",
    "not_evaluated",
    "unreliable",
    "unsupported",
}
_RELIABILITY_ORDER = {
    "calibrated": 0,
    "high": 1,
    "supported": 2,
    "medium": 3,
    "unassessed": 4,
    "unknown": 5,
    "low": 6,
}


def collect_pose_scores(*records: Any) -> list[dict]:
    """Normalize and de-duplicate score evidence from one or more records.

    Structured ``pose_scores`` entries take precedence over legacy flat keys.
    All returned entries contain the required contract fields and a pose
    identifier (``pose_id`` or ``pose_index``).
    """
    structured: list[dict] = []
    legacy: list[dict] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        structured.extend(_structured_scores(record))
        legacy.extend(_legacy_scores(record))

        nested = record.get("docking_scores")
        if isinstance(nested, dict) and nested is not record:
            structured.extend(_structured_scores(nested, parent=record))
            legacy.extend(_legacy_scores(nested, parent=record))

    result: list[dict] = []
    seen: set[tuple] = set()
    for item in [*structured, *legacy]:
        key = (
            item["metric_id"],
            item.get("pose_id"),
            item.get("pose_index"),
            item["method"],
            item.get("method_version"),
            item["stage"],
            item.get("scoring_function"),
            item.get("protocol_id"),
            item.get("protocol_identity_complete"),
        )
        if key in seen:
            continue
        seen.add(key)
        result.append(item)
    return result


def resolve_primary_metric(
    score_groups: Iterable[Iterable[dict]],
    *,
    requested: str | None = None,
    preferred_metric_ids: Iterable[str] = (),
) -> str | None:
    """Choose one comparable metric id for an entire ranking cohort.

    An explicit request is authoritative.  Otherwise a structured score marked
    ``primary`` wins, followed by the caller's compatibility preference and
    finally the reliable metric with greatest cohort coverage.  A metric whose
    unit or direction changes within the cohort is rejected.
    """
    groups = [list(group) for group in score_groups]
    flat = [score for group in groups for score in group if isinstance(score, dict)]
    if not flat:
        return None

    if requested:
        matched = _matching_metric_ids(flat, requested)
        return _best_compatible_metric(flat, matched)

    declared = [
        str(score.get("metric_id"))
        for score in flat
        if score.get("primary") is True and is_rankable_score(score)
    ]
    if declared:
        metric = _best_compatible_metric(flat, set(declared))
        if metric:
            return metric

    for requested_metric in preferred_metric_ids:
        matched = _matching_metric_ids(flat, requested_metric)
        metric = _best_compatible_metric(flat, matched)
        if metric:
            return metric

    candidates = {
        str(score.get("metric_id"))
        for score in flat
        if is_rankable_score(score)
    }
    return _best_compatible_metric(flat, candidates)


def find_score(scores: Iterable[dict], metric_id: str | None) -> dict | None:
    """Return a rankable score matching a canonical id or legacy alias."""
    if not metric_id:
        return None
    candidates = [
        score for score in scores
        if is_rankable_score(score) and metric_matches(score.get("metric_id"), metric_id)
    ]
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda item: (
            _RELIABILITY_ORDER.get(str(item.get("reliability")), 50),
            str(item.get("method_version") or ""),
        ),
    )


def score_sort_key(score: dict | None) -> tuple[int, float]:
    """Ascending sort component with missing scores always last."""
    if not is_rankable_score(score):
        return (1, 0.0)
    value = float(score["value"])
    if score["direction"] == "higher":
        value = -value
    return (0, value)


def metric_matches(actual: Any, requested: Any) -> bool:
    """Match canonical ids while accepting historical unqualified field names."""
    actual_text = str(actual or "").strip()
    requested_text = str(requested or "").strip()
    if not actual_text or not requested_text:
        return False
    if actual_text.casefold() == requested_text.casefold():
        return True
    return actual_text.rsplit(".", 1)[-1].casefold() == requested_text.casefold()


def is_rankable_score(score: Any) -> bool:
    if not isinstance(score, dict):
        return False
    if score.get("direction") not in {"higher", "lower"}:
        return False
    if score.get("applicability") != "applicable":
        return False
    if str(score.get("reliability", "")).casefold() in _UNUSABLE_RELIABILITY:
        return False
    value = score.get("value")
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _structured_scores(record: dict, parent: dict | None = None) -> list[dict]:
    raw_scores = record.get("pose_scores")
    if not isinstance(raw_scores, list):
        return []
    context = {**(parent or {}), **record}
    normalized = []
    for raw in raw_scores:
        if not isinstance(raw, dict):
            continue
        item = _normalize_structured_score(raw, context)
        if item is not None:
            normalized.append(item)
    return normalized


def _normalize_structured_score(raw: dict, context: dict) -> dict | None:
    metric_id = str(raw.get("metric_id") or raw.get("name") or "").strip()
    if not metric_id:
        return None
    method = str(
        raw.get("method")
        or context.get("method")
        or context.get("engine")
        or metric_id.split(".", 1)[0]
        or "unknown"
    ).strip()
    if "." not in metric_id and method and method != "unknown":
        metric_id = f"{method}.{metric_id}"

    direction = _normalize_direction(
        raw.get("direction", raw.get("better", "context_only")),
    )
    pose_id, pose_index = _pose_identity(raw, context)
    item = {
        "schema_version": SCHEMA_VERSION,
        "metric_id": metric_id,
        "value": _finite_number(raw.get("value")),
        "unit": str(raw.get("unit") or "dimensionless"),
        "direction": direction,
        "method": method or "unknown",
        "method_version": str(
            raw.get("method_version")
            or context.get("method_version")
            or context.get("engine_version")
            or "unknown"
        ),
        "stage": str(raw.get("stage") or context.get("stage") or "pose_evaluation"),
        "applicability": _normalize_applicability(
            raw.get("applicability", raw.get("applicable", "applicable")),
        ),
        "reliability": _normalize_reliability(
            raw.get("reliability", "unassessed"),
        ),
    }
    _set_pose_identity(item, pose_id, pose_index)
    if raw.get("primary") is True or raw.get("is_primary") is True:
        item["primary"] = True
    if raw.get("uncertainty") is not None:
        item["uncertainty"] = raw["uncertainty"]
    for field in ("scoring_function", "protocol_id"):
        value = raw.get(field, context.get(field))
        if value not in (None, ""):
            item[field] = str(value)
    protocol_complete = raw.get(
        "protocol_identity_complete",
        context.get("protocol_identity_complete"),
    )
    if isinstance(protocol_complete, bool):
        item["protocol_identity_complete"] = protocol_complete
    return item


def _legacy_scores(record: dict, parent: dict | None = None) -> list[dict]:
    context = {**(parent or {}), **record}
    pose_id, pose_index = _pose_identity(record, parent or {})
    result = []
    for field, spec in LEGACY_SCORE_SPECS.items():
        if field not in record:
            continue
        value = _finite_number(record.get(field))
        if value is None:
            continue
        method = _legacy_method(spec["method"], context)
        metric_id = spec["metric_id"]
        if field == "docking_score" and method != "unidock":
            metric_id = f"{method}.docking_score"
        elif field == "engine_pose_rank":
            metric_id = f"{method}.engine_pose_rank"
        item = {
            "schema_version": SCHEMA_VERSION,
            "metric_id": metric_id,
            "value": value,
            "unit": spec["unit"],
            "direction": spec["direction"],
            "method": method,
            "method_version": str(
                context.get("method_version")
                or context.get("engine_version")
                or context.get(f"{method}_version")
                or "unknown"
            ),
            "stage": str(context.get("stage") or "docking"),
            "applicability": "applicable",
            "reliability": str(context.get("score_reliability") or "supported"),
        }
        _set_pose_identity(item, pose_id, pose_index)
        for provenance_field in ("scoring_function", "protocol_id"):
            provenance_value = context.get(provenance_field)
            if provenance_value not in (None, ""):
                item[provenance_field] = str(provenance_value)
        protocol_complete = context.get("protocol_identity_complete")
        if isinstance(protocol_complete, bool):
            item["protocol_identity_complete"] = protocol_complete
        result.append(item)
    return result


def _legacy_method(default: str, context: dict) -> str:
    engine = str(context.get("engine") or context.get("method") or "").strip()
    if not engine:
        return default
    aliases = {
        "vina": "autodock_vina",
        "autodock-vina": "autodock_vina",
        "uni-dock": "unidock",
    }
    return aliases.get(engine.casefold(), engine.casefold())


def _pose_identity(record: dict, context: dict) -> tuple[str | None, int | None]:
    pose_id = record.get("pose_id")
    if pose_id in (None, ""):
        pose_id = record.get("selected_pose_id")
    if pose_id in (None, ""):
        pose_id = context.get("pose_id") or context.get("selected_pose_id")

    pose_index = _integer(
        record.get(
            "pose_index",
            record.get("selected_pose_index"),
        ),
    )
    if pose_index is None:
        pose_index = _integer(
            context.get(
                "pose_index",
                context.get("selected_pose_index"),
            ),
        )

    if pose_id in (None, "") and pose_index is None:
        name = (
            record.get("mol_name")
            or record.get("compound_id")
            or context.get("mol_name")
            or context.get("compound_id")
        )
        if name not in (None, ""):
            pose_id = str(name)
    return (str(pose_id) if pose_id not in (None, "") else None, pose_index)


def _set_pose_identity(item: dict, pose_id: str | None, pose_index: int | None) -> None:
    if pose_id is not None:
        item["pose_id"] = pose_id
    if pose_index is not None:
        item["pose_index"] = pose_index
    if pose_id is None and pose_index is None:
        # A score without pose lineage is not safe to join or rank.  Retain it
        # as evidence but make the missing identity explicit and non-rankable.
        item["pose_id"] = "unknown"
        item["applicability"] = "not_evaluable"
        item["reliability"] = "invalid"


def _normalize_direction(value: Any) -> str:
    text = str(value or "").strip().casefold()
    aliases = {
        "maximize": "higher",
        "maximum": "higher",
        "high": "higher",
        "higher_is_better": "higher",
        "minimize": "lower",
        "minimum": "lower",
        "low": "lower",
        "lower_is_better": "lower",
        "none": "context_only",
    }
    normalized = aliases.get(text, text)
    return normalized if normalized in {"higher", "lower", "context_only"} else "context_only"


def _normalize_applicability(value: Any) -> str:
    if value is True:
        return "applicable"
    if value is False:
        return "not_applicable"
    text = str(value or "").strip().casefold()
    aliases = {
        "yes": "applicable",
        "true": "applicable",
        "no": "not_applicable",
        "false": "not_applicable",
        "n/a": "not_applicable",
        "na": "not_applicable",
        "unsupported": "not_evaluable",
    }
    return aliases.get(text, text or "unknown")


def _normalize_reliability(value: Any) -> str:
    text = str(value or "").strip().casefold()
    return text or "unknown"


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _integer(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) and number.is_integer() else None


def _matching_metric_ids(scores: list[dict], requested: str) -> set[str]:
    return {
        str(score.get("metric_id"))
        for score in scores
        if metric_matches(score.get("metric_id"), requested)
    }


def _best_compatible_metric(scores: list[dict], candidates: set[str]) -> str | None:
    ranked: list[tuple[int, int, str]] = []
    for metric_id in candidates:
        metric_scores = [
            score for score in scores
            if score.get("metric_id") == metric_id and is_rankable_score(score)
        ]
        if not metric_scores:
            continue
        signatures = {
            (
                score.get("unit"),
                score.get("direction"),
                score.get("method"),
                score.get("method_version"),
                score.get("protocol_id")
                or (
                    f"scoring_function:{score.get('scoring_function')}"
                    if score.get("scoring_function") not in (None, "")
                    else None
                ),
            )
            for score in metric_scores
        }
        if len(signatures) != 1:
            continue
        coverage = len(metric_scores)
        reliability = Counter(
            str(score.get("reliability") or "unknown") for score in metric_scores
        )
        reliability_penalty = sum(
            count * _RELIABILITY_ORDER.get(label, 50)
            for label, count in reliability.items()
        )
        ranked.append((-coverage, reliability_penalty, metric_id))
    return min(ranked)[2] if ranked else None

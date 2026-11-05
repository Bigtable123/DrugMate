"""Compact KEEP/MODIFY review tables; never release or recompute a rank."""
from __future__ import annotations

import math

from pose_identity import benchmark_label


def build_candidate_reviews(records: list[dict], *, source_round: int | None) -> dict:
    """Select five scored rows per class inside each scientific cohort/pool."""
    partitions: dict[tuple[str, str], list[dict]] = {}
    for record in records:
        if not isinstance(record, dict) or benchmark_label(record) == "DEC":
            continue
        if record.get("classification") not in {"keep", "modify"}:
            continue
        key = (
            str(record.get("cohort_id") or "cohort:unspecified"),
            str(record.get("normalization_pool") or "nomination_candidates"),
        )
        partitions.setdefault(key, []).append(record)

    groups = []
    for (cohort_id, pool), members in sorted(partitions.items()):
        for classification in ("keep", "modify"):
            classified = [r for r in members if r.get("classification") == classification]
            scored = [r for r in classified if _finite_score(r.get("overall_score")) is not None]
            scored.sort(key=lambda r: (
                -float(r["overall_score"]),
                str(r.get("mol_name") or ""),
                str(r.get("pose_id") or ""),
            ))
            groups.append({
                "cohort_id": cohort_id,
                "normalization_pool": pool,
                "classification": classification,
                "scored_count": len(scored),
                "unscored_count": len(classified) - len(scored),
                "rows": [_review_row(r, position) for position, r in enumerate(scored[:5], 1)],
            })
    return {
        "review_only": True,
        "source_round": source_round,
        "source": "step6/candidate_ranking.json",
        "order_by": "existing_overall_score_descending",
        "tie_break_fields": ["mol_name", "pose_id"],
        "position_field": "review_position",
        "top_n_per_classification": 5,
        "groups": groups,
    }


def _finite_score(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    return score if math.isfinite(score) else None


def _review_row(record: dict, position: int) -> dict:
    molecular_findings = []
    missing_metrics = set(record.get("planned_but_failed_metric_ids") or [])
    for category_id, category in sorted((record.get("category_scores") or {}).items()):
        if not isinstance(category, dict):
            continue
        status = category.get("category_verdict")
        failures = category.get("fail_metric_ids") or []
        flags = category.get("flag_metric_ids") or []
        if status in {"modify", "reject"} or failures or flags:
            detail = [f"{category_id}: {status or 'unreported'}"]
            if failures:
                detail.append("fail=" + ", ".join(map(str, failures)))
            if flags:
                detail.append("flag=" + ", ".join(map(str, flags)))
            molecular_findings.append("; ".join(detail))
        missing_metrics.update(category.get("missing_metric_ids") or [])
    evidence_gaps = [f"missing metric: {metric}" for metric in sorted(missing_metrics)]
    evidence_gaps.extend(
        f"missing step: {step}" for step in (record.get("missing_active_steps") or [])
    )
    return {
        "review_position": position,
        **{key: record.get(key) for key in (
            "mol_name", "pose_id", "smiles", "classification",
            "overall_score", "overall_rank", "rank_status", "nomination_status",
        )},
        "molecular_findings": molecular_findings,
        "evidence_gaps": evidence_gaps,
        "protocol_blockers": list(record.get("ranking_blocked_by") or []),
    }

"""Compound-level pose selection and score-aware candidate ranking."""
from __future__ import annotations

import sys
from collections import defaultdict
from itertools import groupby
import math
import random
from pathlib import Path

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from io_utils import number as _number
from pipeline_registry import DOCKING_PRIMARY_ENGINES
from pose_identity import benchmark_label
from pose_scores import (
    CANDIDATE_RANKING_PREFERENCE,
    LEGACY_SCORE_SPECS,
    POSE_SELECTION_PREFERENCE,
    collect_pose_scores,
    find_score,
    metric_matches,
    resolve_primary_metric,
    score_sort_key,
)


EVAL_ORDER = {
    "pass": 0,
    "warn": 1,
    "fail": 2,
    "not_evaluated": 3,
    "not_evaluable": 3,
    "unknown": 4,
}
# Sort order, not severity.  ``not_evaluated`` is listed last to match this
# file's own ``EVAL_ORDER`` convention (unlabeled states sort after ``fail``): a
# molecule with no evidence has proven nothing, so it must not outrank a
# ``modify`` that was actually measured.  This is deliberately not the same as
# the keep < not_evaluated < modify < reject severity used to *decide* a label.
# Making the entry explicit changes nothing -- the lookup already defaulted
# unknown labels to 3 -- it only records the intent.
CLASS_ORDER = {"keep": 0, "modify": 1, "reject": 2, "not_evaluated": 3}
BLUEPRINT_ORDER = {
    "pass": 0,
    "not_applicable": 1,
    "unknown": 1,
    "unresolved": 2,
    "fail": 3,
}


def select_compound_poses(
    eval_results: list[dict],
    docking_scores: list[dict] | None = None,
    primary_metric_id: str | None = None,
    authoritative_pose_ids: set[str] | None = None,
    eligible_verdicts: set[str] | frozenset[str] | None = None,
    authoritative_pose_records: dict[str, dict] | None = None,
) -> list[dict]:
    """Keep the evaluator's representative and its exact-pose property evidence.

    A supplied manifest is authoritative: downstream property calculations have
    already consumed those poses. Reports without a selection contract may use
    the best eligible Step-5 StepScore, or the QC-first fallback without scores.
    """
    all_poses = prepare_pose_records(eval_results, docking_scores)
    requested_metric = primary_metric_id or _configured_primary_metric(all_poses)
    primary_metric = resolve_primary_metric(
        (_entry_scores(pose) for pose in all_poses),
        requested=requested_metric,
        preferred_metric_ids=POSE_SELECTION_PREFERENCE,
    )

    authoritative_ids = (
        {str(value) for value in authoritative_pose_ids if str(value)}
        if authoritative_pose_ids is not None
        else None
    )
    if authoritative_ids is not None and eligible_verdicts is None:
        raise ValueError(
            "authoritative pose selection requires its declared "
            "eligible_verdicts policy"
        )
    normalized_eligible_verdicts = (
        {
            str(value).strip()
            for value in eligible_verdicts or ()
            if str(value).strip()
        }
        if authoritative_ids is not None
        else {"pass", "warn", "not_evaluated"}
    )
    if authoritative_ids is not None:
        available_ids = {
            str(pose.get("pose_id"))
            for pose in all_poses
            if pose.get("pose_id") not in (None, "")
        }
        unknown = sorted(authoritative_ids - available_ids)
        if unknown:
            raise ValueError(
                "authoritative selected pose_id is absent from evaluation: "
                + ", ".join(unknown)
            )
        record_ids = set(authoritative_pose_records or {})
        if (
            authoritative_pose_records is not None
            and record_ids != authoritative_ids
        ):
            raise ValueError(
                "authoritative pose metadata and selected pose_id set disagree"
            )

    grouped: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for pose in all_poses:
        pose_index = _score_pose_index(pose)
        fallback_name = (
            f"pose_{pose_index}"
            if pose_index is not None and pose_index >= 0
            else "pose_invalid_index"
        )
        mol_name = str(pose.get("mol_name") or fallback_name)
        # The evaluator creates its post-QC manifest one molecule label at a
        # time.  Re-grouping a modern manifest by another identity domain here
        # can merge two evaluator groups or split one and manufacture an
        # omission.  Use the exact evaluator grouping whenever its contract is
        # present; retain the richer legacy key only for reports predating it.
        group_key = (
            _evaluator_selection_group_key(pose, mol_name)
            if authoritative_ids is not None
            else _compound_group_key(pose, mol_name)
        )
        grouped[group_key].append(pose)

    selected = []
    for _group_key, poses in grouped.items():
        authoritative_matches = (
            [
                pose for pose in poses
                if str(pose.get("pose_id") or "") in authoritative_ids
            ]
            if authoritative_ids is not None
            else []
        )
        eligible = [
            pose for pose in poses
            if str(pose.get("verdict") or "")
            in normalized_eligible_verdicts
        ]
        if len(authoritative_matches) > 1:
            raise ValueError(
                "authoritative selection contains multiple poses for one "
                f"compound group: {_group_key!r}"
            )
        # StepScore selection is only for reports without an authoritative pose.
        step5_pool = [pose for pose in eligible if _step5_pose_is_rankable(pose)]
        if authoritative_matches:
            best = authoritative_matches[0]
            if best not in eligible:
                raise ValueError(
                    "authoritative selection names a structure-ineligible "
                    f"pose: {best.get('pose_id')!r}"
                )
            selection_method = "structure_evaluator_selection_manifest"
        elif authoritative_ids is not None and eligible:
            raise ValueError(
                "authoritative selection omits an eligible compound group: "
                f"{_group_key!r}"
            )
        elif step5_pool:
            best = min(
                step5_pool,
                key=lambda pose: _step5_selection_key(pose, primary_metric),
            )
            selection_method = "coordinator_step5_stepscore"
        else:
            best = min(
                poses,
                key=lambda pose: _pose_selection_key(pose, primary_metric),
            )
            selection_method = (
                "structure_evaluator_diagnostic_pose"
                if authoritative_ids is not None
                else "structure_evaluator_qc_fallback"
            )
        mol_name = str(best.get("mol_name") or f"pose_{best.get('pose_index', 0)}")
        primary_score = find_score(_entry_scores(best), primary_metric)
        verdict_counts = {
            verdict: sum(
                str(pose.get("verdict") or pose.get("eval_verdict") or "unknown")
                == verdict
                for pose in poses
            )
            for verdict in ("pass", "warn", "fail", "not_evaluated")
        }
        unknown_count = len(poses) - sum(verdict_counts.values())
        if unknown_count:
            verdict_counts["unknown"] = unknown_count
        selected_entry = {
            **best,
            "mol_name": mol_name,
            "num_poses_considered": len(poses),
            "pose_verdicts": verdict_counts,
            "primary_score_metric_id": primary_metric,
            "primary_score": primary_score,
            "representative_pose_selection": selection_method,
            "selected_step5_score": _step5_score(best),
            "selected_step5_rank_status": best.get("rank_status"),
        }
        authoritative_record = (authoritative_pose_records or {}).get(
            str(best.get("pose_id") or "")
        )
        if authoritative_record is not None:
            selected_entry.update({
                "selection_rank": authoritative_record["selection_rank"],
                "post_qc_verdict": authoritative_record["post_qc_verdict"],
                "authoritative_pose_selection": dict(authoritative_record),
            })
        selected.append(selected_entry)
    return selected


def _evaluator_selection_group_key(
    entry: dict,
    mol_name: str,
) -> tuple[str, str, str, str]:
    """Reproduce the evaluator's molecule-label selection partition."""
    return (
        str(entry.get("cohort_id") or "cohort:unspecified"),
        str(
            entry.get("normalization_pool")
            or (
                "benchmark_decoys"
                if benchmark_label(entry) == "DEC"
                else "nomination_candidates"
            )
        ),
        "mol_name",
        mol_name,
    )


def prepare_pose_records(
    eval_results: list[dict],
    docking_scores: list[dict] | None = None,
    *,
    docking_evidence_unavailable_reason: str | None = None,
) -> list[dict]:
    """Attach canonical docking evidence while preserving every input pose.

    Passing ``docking_scores`` (including an empty list) declares it to be the
    verified Step 3 source of truth.  In that mode, duplicated docking values
    embedded in evaluator artifacts are removed before the join, so a missing
    or invalid Step 3 contract cannot be bypassed by a downstream copy.
    """
    authoritative_docking = docking_scores is not None
    scores_by_identity, scores_by_pose, scores_by_stable_identity = (
        _index_score_records(docking_scores or [])
    )
    prepared: list[dict] = []
    for entry in eval_results:
        merged = dict(entry)
        if authoritative_docking:
            _remove_docking_evidence(merged)
        pose_index = _score_pose_index(entry)
        if pose_index is None or pose_index < 0:
            _discard_invalid_pose_index_docking_evidence(
                merged,
                raw_pose_index=entry.get("pose_index"),
            )
            prepared.append(merged)
            continue
        mol_name = str(entry.get("mol_name") or f"pose_{pose_index}")
        identity_conflict: dict | None = None
        join_issue_reason = "pose_index_label_or_stable_identity_conflict"
        score_record = _score_record_by_stable_identity(
            entry,
            pose_index,
            scores_by_stable_identity,
        )
        if score_record is None:
            exact_matches = scores_by_identity.get((mol_name, pose_index), [])
            if len(exact_matches) == 1 and _stable_identities_compatible(
                entry, exact_matches[0]
            ):
                score_record = exact_matches[0]
            elif exact_matches:
                identity_conflict = exact_matches[0]
                if len(exact_matches) > 1:
                    join_issue_reason = "ambiguous_duplicate_score_identity"
        if score_record is None and identity_conflict is None:
            pose_records = scores_by_pose.get(pose_index, [])
            if len(pose_records) == 1:
                candidate = pose_records[0]
                # A bare integer is only a safe legacy key when it does not
                # contradict an explicit molecule label.  This used to accept
                # e.g. M1K000381's unique score row for M1K000761 solely because
                # both happened to carry pose_index=1140.
                if (
                    _record_names_compatible(entry, candidate)
                    and _stable_identities_compatible(entry, candidate)
                ):
                    score_record = candidate
                else:
                    identity_conflict = candidate
            elif pose_records:
                identity_conflict = pose_records[0]
                join_issue_reason = "ambiguous_duplicate_pose_index"
        if score_record is None and identity_conflict is not None:
            _discard_conflicted_docking_evidence(
                merged,
                expected=entry,
                observed=identity_conflict,
                pose_index=pose_index,
                reason=join_issue_reason,
            )
            prepared.append(merged)
            continue
        if score_record is None and authoritative_docking:
            reason = (
                docking_evidence_unavailable_reason
                or "no committed docking score matches this evaluated pose"
            )
            merged["docking_evidence_status"] = "unavailable"
            merged["docking_score_join"] = {
                "status": "unavailable",
                "reason": reason,
                "pose_index": pose_index,
                "expected_mol_name": mol_name,
            }
            metric_context = merged.get("metric_context")
            metric_context = (
                dict(metric_context) if isinstance(metric_context, dict) else {}
            )
            blockers = list(metric_context.get("protocol_blockers") or [])
            if "docking_evidence_unavailable" not in blockers:
                blockers.append("docking_evidence_unavailable")
            metric_context["protocol_blockers"] = blockers
            merged["metric_context"] = metric_context
        merged_scores = _merge_score_payload(merged, score_record)
        if merged_scores:
            merged["docking_scores"] = merged_scores
        pose_scores = collect_pose_scores(merged, score_record or {})
        if pose_scores:
            merged["pose_scores"] = pose_scores
            merged.setdefault("docking_scores", {})["pose_scores"] = pose_scores
        if score_record is not None and authoritative_docking:
            merged["docking_evidence_status"] = "verified"
            merged["docking_score_join"] = {
                "status": "verified",
                "pose_index": pose_index,
                "mol_name": mol_name,
            }
        prepared.append(merged)
    return prepared


def _compound_group_key(
    entry: dict,
    mol_name: str,
) -> tuple[str, str, str, str]:
    """Group by scientific cohort, normalization pool, and identity."""
    cohort_id = str(entry.get("cohort_id") or "cohort:unspecified")
    normalization_pool = str(
        entry.get("normalization_pool")
        or (
            "benchmark_decoys"
            if benchmark_label(entry) == "DEC"
            else "nomination_candidates"
        )
    )
    chemical_identity = entry.get("chemical_identity")
    chemical_identity_id = (
        chemical_identity.get("identity_id")
        if isinstance(chemical_identity, dict)
        else None
    ) or entry.get("chemical_identity_id")
    if chemical_identity_id not in (None, ""):
        return (
            cohort_id,
            normalization_pool,
            "chemical_identity_id",
            str(chemical_identity_id),
        )
    chemical_id = entry.get("chemical_id")
    if chemical_id not in (None, ""):
        return cohort_id, normalization_pool, "chemical_id", str(chemical_id)
    return cohort_id, normalization_pool, "mol_name", mol_name


def rank_candidates(
    records: list[dict],
    primary_metric_id: str | None = None,
) -> list[dict]:
    """Rank candidates for physical output.

    When a hierarchical overall score is present, refresh its released ranks
    after final nomination blocks, preserving cohort partitions and average ties.
    Rejected, protocol-blocked, and insufficient-coverage records remain in the
    audit output but receive no rank. An all-legacy artifact retains its
    historical gate-aware ordering; legacy rows mixed into schema 2 are
    unranked so old weighted totals cannot contaminate OverallRank.
    """
    for record in records:
        scores = _entry_scores(record)
        if scores:
            record["pose_scores"] = scores

    requested_metric = primary_metric_id or _configured_primary_metric(records)
    primary_metric = resolve_primary_metric(
        (_entry_scores(record) for record in records),
        requested=requested_metric,
        preferred_metric_ids=CANDIDATE_RANKING_PREFERENCE,
    )
    for record in records:
        record["primary_score_metric_id"] = primary_metric
        record["primary_score"] = find_score(_entry_scores(record), primary_metric)

    use_hierarchical_scoring = any(_is_schema2_record(record) for record in records)
    if use_hierarchical_scoring:
        # Coverage and nomination blocks may have been added after scoring.
        # Re-rank only previously released scores; a null OverallRank must not
        # be resurrected from an old weighted-score alias.
        pools = defaultdict(list)
        for record in records:
            if not _is_schema2_record(record):
                continue
            score = _finite_rank(record.get("overall_score"))
            if (
                score is not None
                and _finite_rank(record.get("overall_rank")) is not None
                and record.get("rank_status") in {"formal", "provisional"}
                and record.get("classification") != "reject"
                and record.get("gate_status") != "reject"
                and not _ranking_is_blocked(record)
            ):
                pools[_ranking_partition_key(record)].append(record)
            record["overall_rank"] = record["weighted_score_rank"] = None
        for pool in pools.values():
            ordered = sorted(pool, key=lambda row: float(row["overall_score"]), reverse=True)
            position = 1
            for _, tied in groupby(ordered, key=lambda row: float(row["overall_score"])):
                tied = list(tied)
                rank = position + (len(tied) - 1) / 2
                for record in tied:
                    record["overall_rank"] = record["weighted_score_rank"] = rank
                position += len(tied)

        # A flat physical file may contain several independently normalized
        # cohorts/pools.  Keep each partition contiguous before applying its
        # local OverallRank; a local rank of 1 in cohort B is never compared to
        # rank 2 in cohort A.  Legacy weighted_* rows in a mixed artifact are
        # audit-only so their old totals cannot enter the schema-v2 ordering.
        ranked = sorted(
            records,
            key=lambda record: (
                *_ranking_partition_key(record),
                0 if _is_schema2_record(record) else 1,
                *_candidate_ranking_key(
                    record,
                    primary_metric,
                    use_weighted_total=_is_schema2_record(record),
                ),
            ),
        )
        for record in ranked:
            cohort_id, normalization_pool = _ranking_partition_key(record)
            record["cohort_id"] = cohort_id
            record["normalization_pool"] = normalization_pool
            record["rank_partition_id"] = (
                f"{cohort_id}|normalization_pool={normalization_pool}"
            )
            record["rank_scope"] = "within_cohort_and_normalization_pool"
            if _is_schema2_record(record):
                # Presence of schema-v2 is authoritative.  In particular, do
                # not resurrect a null OverallRank from weighted_score_rank.
                record["rank"] = _finite_rank(record.get("overall_rank"))
                record["ranking_schema"] = "hierarchical_v2"
            else:
                record["rank"] = None
                record["ranking_schema"] = "legacy_unranked_in_mixed_schema"
                record["rank_status"] = "legacy_schema_unranked_in_mixed_output"
    else:
        ranked = sorted(
            records,
            key=lambda record: _candidate_ranking_key(
                record,
                primary_metric,
                use_weighted_total=False,
            ),
        )
        trusted_rank = 0
        for record in ranked:
            if _ranking_is_blocked(record):
                record["rank"] = None
                continue
            trusted_rank += 1
            record["rank"] = trusted_rank
    return ranked


def screening_metrics(
    records: list[dict],
    bootstrap_samples: int = 500,
    seed: int = 42,
    primary_metric_id: str | None = None,
) -> dict:
    """Compute retrospective metrics from ACT_/DEC_ labels under one protocol.

    Top-20 observations retain input order within exact score ties. A tie across
    the cutoff also reports the expected hit count under uniform tie ordering
    and the feasible hit range; the observed ordering follows the supplied rows.
    """
    labeled_records = [
        record
        for record in records
        if benchmark_label(record) is not None
    ]
    expected_actives = sum(
        benchmark_label(record) == "ACT"
        for record in labeled_records
    )
    expected_decoys = len(labeled_records) - expected_actives
    coverage_metadata = {
        "num_expected_actives": expected_actives,
        "num_expected_decoys": expected_decoys,
        "num_expected_records": len(labeled_records),
    }
    score_groups = [_docking_scores(record) for record in labeled_records]
    configured_metric = primary_metric_id or _configured_primary_metric(records)
    requested_metric = (
        configured_metric
        if configured_metric
        and any(
            metric_matches(score.get("metric_id"), configured_metric)
            for group in score_groups
            for score in group
        )
        else None
    )
    primary_metric = _full_coverage_docking_metric(
        score_groups,
        requested_metric=requested_metric,
    )
    if not expected_actives or not expected_decoys:
        return {
            "available": False,
            "reason": "both ACT_ and DEC_ labels are required",
            "score": primary_metric,
            "num_scored_records": 0,
            **coverage_metadata,
        }
    if primary_metric is None:
        return {
            "available": False,
            "reason": (
                "one common protocol-compatible docking score is required "
                "for every ACT_/DEC_ record"
            ),
            "score": None,
            "num_scored_records": 0,
            **coverage_metadata,
        }
    labeled = []
    direction = None
    selected_scores = []
    missing_score_records = []
    for record in labeled_records:
        name = str(record.get("mol_name", ""))
        score = find_score(_docking_scores(record), primary_metric)
        if score is None:
            missing_score_records.append(name)
            continue
        selected_scores.append(score)
        direction = score["direction"]
        rank_value = float(score["value"])
        if direction == "lower":
            rank_value = -rank_value
        labeled.append((1 if benchmark_label(record) == "ACT" else 0, rank_value))
    if missing_score_records:
        return {
            "available": False,
            "reason": (
                f"{len(missing_score_records)} labeled records are missing "
                f"the pre-registered docking score {primary_metric}"
            ),
            "score": primary_metric,
            "num_scored_records": len(selected_scores),
            "missing_score_records": missing_score_records[:50],
            **coverage_metadata,
        }
    protocol_ids = {
        str(score.get("protocol_id"))
        for score in selected_scores
        if score.get("protocol_id") not in (None, "")
    }
    incomplete_protocol = any(
        score.get("protocol_identity_complete") is not True
        for score in selected_scores
    )
    if len(protocol_ids) != 1 or incomplete_protocol:
        return {
            "available": False,
            "reason": (
                "every ACT_/DEC_ score requires one complete shared docking "
                "protocol identity"
            ),
            "score": primary_metric,
            "num_scored_records": len(selected_scores),
            "docking_protocol_ids": sorted(protocol_ids),
            **coverage_metadata,
        }
    positives = [item for item in labeled if item[0] == 1]
    negatives = [item for item in labeled if item[0] == 0]
    if not positives or not negatives:
        return {
            "available": False,
            "reason": "ACT_/DEC_ labels and one common rankable primary score are required",
            "score": primary_metric,
            "num_scored_records": len(selected_scores),
            **coverage_metadata,
        }

    auc = _auc(labeled)
    pr_auc = _average_precision(labeled)
    rng = random.Random(seed)
    boot_auc = []
    boot_pr_auc = []
    for _ in range(bootstrap_samples):
        sample = [rng.choice(positives) for _ in positives] + [rng.choice(negatives) for _ in negatives]
        boot_auc.append(_auc(sample))
        boot_pr_auc.append(_average_precision(sample))
    auc_ci = _bootstrap_interval(boot_auc)
    pr_auc_ci = _bootstrap_interval(boot_pr_auc)

    ranked = sorted(labeled, key=lambda item: item[1], reverse=True)
    active_fraction = len(positives) / len(labeled)
    enrichment = {}
    for percent in (1, 5, 10):
        n = max(1, math.ceil(len(ranked) * percent / 100))
        found = sum(label for label, _ in ranked[:n])
        enrichment[f"EF{percent}%"] = round((found / n) / active_fraction, 3)
        enrichment[f"top{percent}%_actives"] = found
        enrichment[f"top{percent}%_decoy_fraction"] = round(
            (n - found) / n, 3,
        )
        enrichment[f"top{percent}%_size"] = n
    return {
        "available": True,
        "num_actives": len(positives),
        "num_decoys": len(negatives),
        "num_scored_records": len(selected_scores),
        "docking_protocol_id": next(iter(protocol_ids)),
        **coverage_metadata,
        "AUROC": round(auc, 3),
        "AUROC_bootstrap_95ci": auc_ci,
        "PR_AUC": round(pr_auc, 3),
        "PR_AUC_method": "tie_aware_average_precision",
        "PR_AUC_bootstrap_95ci": pr_auc_ci,
        "average_precision": pr_auc,
        "average_precision_method": "tie_aware_average_precision",
        "active_prevalence": round(active_fraction, 6),
        "PR_AUC_relative_to_prevalence": round(
            pr_auc / active_fraction, 3,
        ),
        "score": primary_metric,
        "score_direction": direction,
        **_top20_metrics(labeled),
        **enrichment,
    }


def _top20_metrics(items: list[tuple[int, float]]) -> dict:
    """Report observed and tie-aware Top-20 hits without inventing a tie rank.

    Scores must already be oriented so larger is better. Stable input order
    defines the observed list only. Expected hits assume every ordering within
    the cutoff tie is equally likely; bounds require no such assumption.
    """
    ranked = sorted(items, key=lambda item: item[1], reverse=True)
    size = min(20, len(ranked))
    observed = sum(label for label, _ in ranked[:size])
    boundary_score = ranked[size - 1][1] if size else None
    above = [label for label, score in ranked if score > boundary_score] if size else []
    boundary = [label for label, score in ranked if score == boundary_score] if size else []
    slots = size - len(above)
    boundary_hits = sum(boundary)
    expected = sum(above) + slots * boundary_hits / len(boundary) if boundary else 0.0
    minimum = sum(above) + max(0, slots - (len(boundary) - boundary_hits))
    maximum = sum(above) + min(slots, boundary_hits)
    return {
        "top20_size": size,
        "top20_actives": observed,
        "top20_hit_ratio": observed / size if size else None,
        "top20_observed_tie_policy": "stable_input_order_within_equal_scores",
        "top20_boundary_tie": bool(boundary and slots < len(boundary)),
        "top20_boundary_tie_size": len(boundary),
        "top20_boundary_slots": slots,
        "top20_expected_actives": expected,
        "top20_expected_hit_ratio": expected / size if size else None,
        "top20_expectation_method": "uniform_order_within_cutoff_score_tie",
        "top20_active_count_range": [minimum, maximum],
        "top20_hit_ratio_range": [minimum / size, maximum / size] if size else None,
    }


def _auc(items: list[tuple[int, float]]) -> float:
    positives = [score for label, score in items if label == 1]
    negatives = [score for label, score in items if label == 0]
    wins = sum((p > n) + 0.5 * (p == n) for p in positives for n in negatives)
    return wins / (len(positives) * len(negatives))


def _average_precision(items: list[tuple[int, float]]) -> float:
    """Tie-aware average precision, reported as the retrospective PR-AUC."""
    positives = sum(label for label, _score in items)
    if positives <= 0:
        return 0.0
    ranked = sorted(items, key=lambda item: item[1], reverse=True)
    true_positives = 0
    observed = 0
    area = 0.0
    index = 0
    while index < len(ranked):
        score = ranked[index][1]
        group_positives = 0
        group_size = 0
        while index < len(ranked) and ranked[index][1] == score:
            group_positives += ranked[index][0]
            group_size += 1
            index += 1
        true_positives += group_positives
        observed += group_size
        if group_positives:
            area += (group_positives / positives) * (
                true_positives / observed
            )
    return area


def _bootstrap_interval(values: list[float]) -> list[float] | None:
    if not values:
        return None
    ordered = sorted(values)
    low = ordered[max(0, int(0.025 * len(ordered)) - 1)]
    high = ordered[min(len(ordered) - 1, int(0.975 * len(ordered)))]
    return [round(low, 3), round(high, 3)]


def _step5_selection_key(entry: dict, primary_metric_id: str | None = None) -> tuple:
    """Representative-pose ordering by Step-5 StepScore.

    ``eligible`` already applied the Step-6 gate (only survivor verdicts reach
    here), so among the survivors the representative is the pose with the best
    Step-5 StepScore.  The remaining terms only break exact StepScore ties and
    reproduce the QC-first fallback ordering so behaviour is continuous when two
    poses score identically.
    """
    verdict_order = {"pass": 0, "warn": 1, "not_evaluated": 2, "fail": 3}
    step5 = _step5_score(entry)
    blueprint = entry.get("blueprint_match") or {}
    desired_score = _finite_rank(
        blueprint.get("desired_ranking_score")
        if isinstance(blueprint, dict)
        else None
    )
    primary_score = find_score(_entry_scores(entry), primary_metric_id)
    source_pose_index = _score_pose_index({
        "pose_index": entry.get("source_pose_index", entry.get("pose_index"))
    })
    return (
        -(step5 if step5 is not None else 0.0),
        verdict_order.get(str(entry.get("verdict") or ""), 4),
        -(desired_score if desired_score is not None else 0.0),
        score_sort_key(primary_score),
        (
            source_pose_index
            if source_pose_index is not None and source_pose_index >= 0
            else math.inf
        ),
    )


def _pose_selection_key(entry: dict, primary_metric_id: str | None = None) -> tuple:
    """Compatibility form of the structure evaluator's QC-first ordering."""
    verdict_order = {
        "pass": 0,
        "warn": 1,
        "not_evaluated": 2,
        "fail": 3,
    }
    blueprint = entry.get("blueprint_match") or {}
    desired_score = _finite_rank(
        blueprint.get("desired_ranking_score")
        if isinstance(blueprint, dict)
        else None
    )
    primary_score = find_score(_entry_scores(entry), primary_metric_id)
    source_pose_index = _score_pose_index({
        "pose_index": entry.get("source_pose_index", entry.get("pose_index"))
    })
    return (
        verdict_order.get(str(entry.get("verdict") or ""), 4),
        -(desired_score if desired_score is not None else 0.0),
        score_sort_key(primary_score),
        (
            source_pose_index
            if source_pose_index is not None and source_pose_index >= 0
            else math.inf
        ),
    )


def _step5_score(entry: dict) -> float | None:
    steps = entry.get("step_scores")
    value = (
        (steps.get("step5") or {}).get("score")
        if isinstance(steps, dict)
        else None
    )
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _step5_pose_is_rankable(entry: dict) -> bool:
    """Return whether Step 5 supplied a usable score for pose selection."""
    step5_gate = str(
        entry.get("gate_status")
        or entry.get("metric_gate")
        or "unknown"
    ).casefold()
    return (
        _step5_score(entry) is not None
        and step5_gate not in {"reject", "fail", "failed"}
        and not _ranking_is_blocked(entry)
    )


def _candidate_ranking_key(
    entry: dict,
    primary_metric_id: str | None = None,
    *,
    use_weighted_total: bool = False,
) -> tuple:
    primary_score = find_score(_entry_scores(entry), primary_metric_id)
    legacy_key = (
        CLASS_ORDER.get(entry.get("classification", "reject"), 3),
        EVAL_ORDER.get(entry.get("eval_verdict", "unknown"), 4),
        _blueprint_gate_order(entry),
        -_desired_blueprint_score(entry),
        score_sort_key(primary_score),
        entry.get("mol_name", ""),
    )
    if not use_weighted_total:
        blocked_key = 1 if _ranking_is_blocked(entry) else 0
        return (blocked_key, *legacy_key)

    weighted_total = _finite_weighted_total(entry)
    overall_rank = _finite_rank(entry.get("overall_rank"))
    rankable = (
        weighted_total is not None
        and overall_rank is not None
        and not _ranking_is_blocked(entry)
        and str(entry.get("classification") or "").casefold() != "reject"
    )
    weighted_key = (
        (0, overall_rank)
        if rankable
        else (1, *legacy_key)
    )
    chemical_identity = entry.get("chemical_identity")
    chemical_identity_id = (
        chemical_identity.get("identity_id")
        if isinstance(chemical_identity, dict)
        else None
    ) or entry.get("chemical_identity_id") or entry.get("chemical_id") or ""
    pose_index = entry.get("selected_pose_index")
    if pose_index is None:
        pose_index = entry.get("pose_index")
    stable_pose_index = _score_pose_index({"pose_index": pose_index})
    if stable_pose_index is None or stable_pose_index < 0:
        stable_pose_index = math.inf
    return (
        *weighted_key,
        str(entry.get("mol_name") or ""),
        str(chemical_identity_id),
        stable_pose_index,
    )


def _ranking_is_blocked(entry: dict) -> bool:
    """Return whether a record may receive a trusted nomination rank."""
    if entry.get("ranking_blocked") is True:
        return True
    status = str(entry.get("nomination_status") or "eligible").casefold()
    return status.startswith("blocked")


def _is_schema2_record(entry: dict) -> bool:
    """Return whether hierarchical fields, including explicit nulls, govern.

    Explicit legacy provenance is authoritative.  CSV rows commonly contain
    empty ``overall_*`` columns, and transitional artifacts may retain stale
    values in those columns; neither may promote a declared schema-1 row into
    schema 2.
    """
    schema_version = entry.get("scoring_schema_version")
    try:
        numeric_schema_version = (
            None if isinstance(schema_version, bool) else float(schema_version)
        )
    except (TypeError, ValueError):
        numeric_schema_version = None
    ranking_schema = str(entry.get("ranking_schema") or "").strip().casefold()
    if numeric_schema_version == 1.0 or ranking_schema.startswith("legacy"):
        return False
    return (
        (numeric_schema_version is not None and numeric_schema_version >= 2.0)
        or "overall_score" in entry
        or "overall_rank" in entry
    )


def _ranking_partition_key(entry: dict) -> tuple[str, str]:
    """Return the non-comparable scientific cohort and normalization pool."""
    cohort_id = str(entry.get("cohort_id") or "cohort:unspecified")
    normalization_pool = entry.get("normalization_pool")
    if normalization_pool in (None, ""):
        normalization_pool = (
            "benchmark_decoys"
            if benchmark_label(entry) == "DEC"
            else "nomination_candidates"
        )
    return cohort_id, str(normalization_pool)


def _finite_rank(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _finite_weighted_total(entry: dict) -> float | None:
    # Presence of the schema-v2 field is authoritative even when it is null;
    # falling back to a stale legacy total would accidentally rank an
    # ineligible or insufficient-coverage record.
    value = (
        entry.get("overall_score")
        if "overall_score" in entry
        else entry.get("weighted_total_score")
    )
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _desired_blueprint_score(entry: dict) -> float:
    """Score desired interactions without confusing them with hard gates."""
    match = entry.get("blueprint_match") or {}
    value = match.get("desired_ranking_score")
    if value is None:
        value = match.get("desired_weighted_score")
    return _number(value, default=-1.0)


def _blueprint_gate_order(entry: dict) -> int:
    match = entry.get("blueprint_match") or {}
    verdict = match.get("verdict")
    if not verdict:
        verdict = {
            "passed": "pass",
            "failed": "fail",
            "not_applicable": "not_applicable",
            "unresolved": "unresolved",
        }.get(str(match.get("hard_gate_status")), "unknown")
    verdict = str(verdict)
    return BLUEPRINT_ORDER.get(verdict, BLUEPRINT_ORDER["unknown"])


def _entry_scores(entry: dict) -> list[dict]:
    return collect_pose_scores(entry)


def _docking_scores(entry: dict) -> list[dict]:
    """Limit retrospective screening validation to native docking engines."""
    result = []
    native_methods = {method.casefold() for method in DOCKING_PRIMARY_ENGINES}
    for score in _entry_scores(entry):
        method = str(score.get("method") or "").casefold()
        if method == "vina":
            method = "autodock_vina"
        if method in native_methods:
            result.append(score)
    return result


def _full_coverage_docking_metric(
    score_groups: list[list[dict]],
    *,
    requested_metric: str | None,
) -> str | None:
    """Resolve one compatible metric, preferring complete labeled coverage."""
    if requested_metric:
        return resolve_primary_metric(score_groups, requested=requested_metric)
    candidate_requests = [
        *CANDIDATE_RANKING_PREFERENCE,
        *sorted({
            str(score.get("metric_id"))
            for group in score_groups
            for score in group
            if score.get("metric_id") not in (None, "")
        }),
    ]
    for requested in dict.fromkeys(candidate_requests):
        metric_id = resolve_primary_metric(score_groups, requested=requested)
        if metric_id and all(
            find_score(group, metric_id) is not None for group in score_groups
        ):
            return metric_id
    return None


def _configured_primary_metric(records: list[dict]) -> str | None:
    configured = []
    for record in records:
        value = record.get("primary_score_metric_id")
        if value in (None, ""):
            docking_scores = record.get("docking_scores")
            if isinstance(docking_scores, dict):
                value = docking_scores.get("primary_score_metric_id")
        if value not in (None, ""):
            configured.append(str(value))
    return configured[0] if configured and len(set(configured)) == 1 else None


def _index_score_records(
    records: list[dict],
) -> tuple[
    dict[tuple[str, int], list[dict]],
    dict[int, list[dict]],
    dict[tuple[str, str], list[dict]],
]:
    by_identity: dict[tuple[str, int], list[dict]] = defaultdict(list)
    by_pose: dict[int, list[dict]] = defaultdict(list)
    by_stable_identity: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for record in records:
        if not isinstance(record, dict):
            continue
        pose_index = _score_pose_index(record)
        name = record.get("mol_name") or record.get("compound_id")
        if name not in (None, "") and pose_index is not None:
            by_identity[(str(name), pose_index)].append(record)
        if pose_index is not None:
            by_pose[pose_index].append(record)
        for field in ("pose_id", "chemical_identity_id", "chemical_id"):
            for value in _stable_identity_values(record, field):
                by_stable_identity[(field, value)].append(record)
    return by_identity, by_pose, by_stable_identity


def _score_record_by_stable_identity(
    entry: dict,
    pose_index: int,
    by_stable_identity: dict[tuple[str, str], list[dict]],
) -> dict | None:
    """Resolve an aliased score row by a shared stable identity when unique."""
    for field in ("pose_id", "chemical_identity_id", "chemical_id"):
        for value in _stable_identity_values(entry, field):
            matches = list(by_stable_identity.get((field, value), []))
            if field != "pose_id":
                # Chemical IDs identify a molecule, not one of its poses.
                matches = [
                    record for record in matches
                    if _score_pose_index(record) == pose_index
                ]
            matches = [
                record for record in matches
                if _stable_identities_compatible(entry, record)
            ]
            if len(matches) == 1:
                return matches[0]
    return None


def _score_pose_index(record: dict) -> int | None:
    value = record.get("pose_index")
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped and stripped.lstrip("-").isdigit():
            return int(stripped)
    return None


def _discard_invalid_pose_index_docking_evidence(
    record: dict,
    *,
    raw_pose_index: object,
) -> None:
    """Preserve an invalid evaluation row without letting it acquire a score."""
    _remove_docking_evidence(record)
    record["docking_evidence_status"] = "invalid"
    record["docking_score_join"] = {
        "status": "unavailable",
        "reason": "invalid_evaluation_pose_index",
        "pose_index": raw_pose_index,
        "expected_mol_name": _record_name(record),
    }
    metric_context = record.get("metric_context")
    metric_context = (
        dict(metric_context) if isinstance(metric_context, dict) else {}
    )
    blockers = list(metric_context.get("protocol_blockers") or [])
    for blocker in (
        "invalid_evaluation_pose_index",
        "docking_evidence_unavailable",
    ):
        if blocker not in blockers:
            blockers.append(blocker)
    metric_context["protocol_blockers"] = blockers
    record["metric_context"] = metric_context


def _record_name(record: dict) -> str | None:
    value = record.get("mol_name") or record.get("compound_id")
    return str(value) if value not in (None, "") else None


def _record_names_compatible(left: dict, right: dict) -> bool:
    left_name = _record_name(left)
    right_name = _record_name(right)
    # This helper guards the weak integer fallback.  If either side lacks a
    # name there is no evidence beyond pose_index, so fail closed.
    return bool(left_name and right_name and left_name == right_name)


def _stable_identity_values(record: dict, field: str) -> tuple[str, ...]:
    if field == "pose_id":
        alternate_pose_ids = record.get("pose_ids")
        if isinstance(alternate_pose_ids, (list, tuple, set)):
            alternates = list(alternate_pose_ids)
        elif alternate_pose_ids in (None, ""):
            alternates = []
        else:
            alternates = [alternate_pose_ids]
        raw_values = [record.get("pose_id"), *alternates]
    else:
        chemical_identity = record.get("chemical_identity")
        nested = chemical_identity if isinstance(chemical_identity, dict) else {}
        nested_field = "identity_id" if field == "chemical_identity_id" else field
        raw_values = [record.get(field), nested.get(nested_field)]
    return tuple(dict.fromkeys(
        str(value) for value in raw_values if value not in (None, "")
    ))


def _stable_identities_compatible(left: dict, right: dict) -> bool:
    """Require coordinate witnesses when either side declares one."""
    left_witness = left.get("sdf_molblock_sha256")
    right_witness = right.get("sdf_molblock_sha256")
    if left_witness not in (None, "") or right_witness not in (None, ""):
        if (
            left_witness in (None, "")
            or right_witness in (None, "")
            or str(left_witness) != str(right_witness)
        ):
            return False
    for field in ("pose_id", "chemical_identity_id", "chemical_id"):
        left_values = set(_stable_identity_values(left, field))
        right_values = set(_stable_identity_values(right, field))
        if left_values and right_values and left_values.isdisjoint(right_values):
            return False
    return True


def _discard_conflicted_docking_evidence(
    record: dict,
    *,
    expected: dict,
    observed: dict,
    pose_index: int,
    reason: str = "pose_index_label_or_stable_identity_conflict",
) -> None:
    """Remove mislabeled docking evidence and turn the conflict into a blocker."""
    _remove_docking_evidence(record)

    expected_name = _record_name(expected)
    observed_name = _record_name(observed)
    record["docking_evidence_status"] = "invalid"
    record["docking_score_join"] = {
        "status": "identity_mismatch",
        "reason": reason,
        "pose_index": pose_index,
        "expected_mol_name": expected_name,
        "observed_mol_name": observed_name,
    }
    metric_context = record.get("metric_context")
    metric_context = dict(metric_context) if isinstance(metric_context, dict) else {}
    blockers = list(metric_context.get("protocol_blockers") or [])
    if "docking_score_identity_mismatch" not in blockers:
        blockers.append("docking_score_identity_mismatch")
    metric_context["protocol_blockers"] = blockers
    record["metric_context"] = metric_context


def _remove_docking_evidence(record: dict) -> None:
    """Delete all known flat and structured copies of Step 3 score evidence."""
    record.pop("docking_scores", None)
    for field in LEGACY_SCORE_SPECS:
        record.pop(field, None)

    docking_metric_ids = {
        str(spec.get("metric_id"))
        for spec in LEGACY_SCORE_SPECS.values()
        if spec.get("metric_id") not in (None, "")
    }
    docking_methods = {
        str(method).casefold() for method in DOCKING_PRIMARY_ENGINES
    } | {"vina"}

    def is_docking_score(value: object) -> bool:
        return bool(
            isinstance(value, dict)
            and (
                str(value.get("stage") or "").lower() == "docking"
                or str(value.get("metric_id") or "") in docking_metric_ids
                or str(value.get("method") or "").casefold()
                in docking_methods
            )
        )

    pose_scores = record.get("pose_scores")
    if isinstance(pose_scores, list):
        retained = [score for score in pose_scores if not is_docking_score(score)]
        if retained:
            record["pose_scores"] = retained
        else:
            record.pop("pose_scores", None)
    if (
        is_docking_score(record.get("primary_score"))
        or str(record.get("primary_score_metric_id") or "")
        in docking_metric_ids
    ):
        record.pop("primary_score", None)
        record.pop("primary_score_metric_id", None)

def _merge_score_payload(entry: dict, score_record: dict | None) -> dict:
    existing = entry.get("docking_scores")
    payload = dict(existing) if isinstance(existing, dict) else {}
    if isinstance(score_record, dict):
        # Preserve the matched source's complete score context. A field whitelist
        # lost protocol completeness and created a second, incomplete legacy
        # score during normalization. False and unknown values stay unchanged.
        payload.update(score_record)
    return payload

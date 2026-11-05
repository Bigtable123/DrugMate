"""Publish one round's scored cohort as the reader-facing output tree.

The tree is an isolated *view*: the step working directories keep owning the raw
tool output, and every table here references them by publication-relative path instead of
copying a pose or a trajectory twice (the contract defines ``source_map.json`` as
"tool_run_id, the metric columns, and the relative path", which only makes sense
if the path points at the original).

Three conventions cover pipeline outputs the directory demo does not name:

* An empty cell means "not measured".  Never ``0``, never the string ``NA`` --
  a number in a cell is always a measurement.
* ``*_table.tsv`` is the machine-readable flat table for its layer;
  ``*_report.tsv`` is the reviewer-facing one, carrying the score, its rank and
  the evidence behind it rather than every raw column.
* A tool with neither a measured metric nor execution evidence gets no
  directory. A failed attempted run still gets an audit summary, even when all
  of its scalar metrics are blank.

Round and backtracking artifacts are deliberately absent: that logic does not
exist in the pipeline yet, so publishing a placeholder for it would be a lie.
"""

from __future__ import annotations

import fcntl
import json
import os
import shutil
import sys
import tempfile
import uuid
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

_SHARED = Path(__file__).resolve().parents[2] / "shared"
if str(_SHARED) not in sys.path:
    sys.path.insert(0, str(_SHARED))

import output_layout as ol
from metric_scoring import (  # noqa: E402  (path bootstrap must precede)
    CATEGORY_CONFIG,
    POSE_PROPERTY_METRIC_RULES,
    STEP_IDS,
)

NA = ol.NA

# ---------------------------------------------------------------------------
# The rule catalog, indexed the three ways the tables read it.
# ---------------------------------------------------------------------------

STEP_DIR_KEYS = frozenset(STEP_IDS)

# The legacy runner directories are an implementation detail.  The publication
# groups them by the scientific workflow a reader expects to follow.
_PIPELINE_STEP_NAMES: dict[str, str] = {
    "step0": "structure_prep",
    "step1": "pocket_analysis",
    "step2": "molecule_generation",
    "step2_5": "property_prefilter",
    "step3": "docking",
    "step4": "pose_evaluation",
    "step5": "property_filter",
    "step6": "coordination",
    "step7": "reporting",
}

_WORKFLOW_STAGE_SOURCES: dict[str, tuple[str, ...]] = {
    ol.STRUCTURE_PREP_DIR: ("step0",),
    ol.POCKET_ANALYSIS_DIR: ("step1",),
    ol.MOLECULE_GENERATION_DIR: ("step2",),
    # The pre-docking property pass is published as its own phase below.  The
    # main property scoring view is derived from the post-pose property pass.
    ol.SM_PROPERTY_DIR: ("step5",),
    ol.DOCKING_DIR: ("step3",),
    ol.POST_DOCKING_DIR: ("step4",),
    ol.FINAL_REPORT_DIR: ("step6", "step7"),
}

_CONTROL_INPUT_FILES: tuple[str, ...] = (
    "pipeline_config.json",
    "resolved_run_plan.json",
    "tool_runtime_report.json",
)

_FINALIZATION_RESULT_FILES: tuple[str, ...] = (
    "step6/result.json",
    "step7/result.json",
)

_PUBLICATION_MODES = frozenset({"standalone", "pipeline_final"})

_ARTIFACT_INDEX_COLUMNS: tuple[str, ...] = (
    "pipeline_step_id",
    "declared_by",
    "artifact_name",
    "kind",
    "producer",
    "path",
    "exists",
    "bytes",
    "sha256",
)

_PATH_SUFFIXES = frozenset({
    ".cif", ".csv", ".html", ".json", ".jsonl", ".log", ".md",
    ".mol", ".mol2", ".pdb", ".pdbqt", ".sdf", ".smi", ".tsv",
    ".txt", ".xml", ".xtc", ".yaml", ".yml", ".zst",
})


def _rules_by_category(
    rules: Sequence[Mapping[str, Any]],
) -> dict[str, list[Mapping[str, Any]]]:
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for rule in rules:
        code = str(rule.get("category_code") or "")
        if code in CATEGORY_CONFIG:
            grouped.setdefault(code, []).append(rule)
    return grouped


def _rules_by_tool(
    rules: Sequence[Mapping[str, Any]],
) -> dict[tuple[str, str], list[Mapping[str, Any]]]:
    """Group by physical tool-run directory, not planner switch.

    Several planner switches (for example Lipinski, QED and SA) deliberately
    share the ``rdkit_descriptors`` archive.  Grouping on the switches would
    write that path repeatedly and leave only the last metric in ``tool_table``.
    """
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for rule in rules:
        step_id = str(rule.get("step_id") or "")
        if step_id in STEP_DIR_KEYS:
            tool_run_id = ol.tool_dir_name(ol.rule_tool_id(rule))
            grouped.setdefault((step_id, tool_run_id), []).append(rule)
    return grouped


def _assessments(record: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    """``metric_id`` -> its assessment for one record."""
    return {
        str(item.get("metric_id")): item
        for item in record.get("metric_assessments") or ()
        if item.get("metric_id")
    }


def _measured(assessment: Mapping[str, Any] | None) -> bool:
    """True when the metric was actually observed for this row."""
    if not assessment:
        return False
    return assessment.get("availability_status") == "available"


def _measured_metric_ids(
    records: Sequence[Mapping[str, Any]], metric_ids: Sequence[str],
) -> list[str]:
    """The subset some row actually measured, in the order given.

    The contract drops the directory of a tool that did not run; a column of blanks is
    the same claim at table scale, so it is dropped for the same reason.  The
    final ``all_metrics_table.tsv`` keeps the full catalog width, so nothing
    disappears from the run -- only from the layer that has nothing to say.
    """
    per_row = [_assessments(record) for record in records]
    return [
        metric_id for metric_id in metric_ids
        if any(_measured(assessments.get(metric_id)) for assessments in per_row)
    ]


# ---------------------------------------------------------------------------
# Column resolution.  A published column name is a contract with the reader, so
# it never changes to match whatever the producing stage happened to call it.
# ---------------------------------------------------------------------------

#: Published name -> the record keys that may carry it, best first.  A dotted
#: key walks into a nested block.
_ALIASES: dict[str, tuple[str, ...]] = {
    "ligand_id": ("ligand_id", "chemical_identity_id", "mol_name"),
    "molecule_state_id": ("molecule_state_id", "chemical_identity_id", "mol_name"),
    # The column promises a canonical SMILES, so the RDKit-canonicalised string
    # computed for pose identity comes first; ``smiles`` is the library's raw
    # input string and only stands in when identity could not be built.
    "canonical_SMILES": (
        "canonical_SMILES",
        "chemical_identity.canonical_isomeric_smiles",
        "smiles",
    ),
    "engine": ("engine", "docking_engine"),
    "profile_id": ("profile_id", "protocol_profile_id"),
    "pose_file": ("pose_file", "pose_path"),
    "artifact_path": ("artifact_path", "pose_file", "pose_path"),
    "planned_weighted_coverage_overall": ("planned_category_weighted_coverage",),
    "confidence_grade": ("coverage_status",),
    # ``uncertainty_grade`` deliberately has no alias.  It asks how precise the
    # numbers are (MM-GBSA single-frame endpoints, missing SEM, absent
    # convergence evidence); ``coverage_status`` answers how many metrics ran.
    # Pointing both at coverage shipped one value under two headings and implied
    # two independent judgements, so this stays empty until a real producer
    # exists.
    "missing_required_metrics": ("ranking_blocked_by",),
    "reason_codes": ("reason_codes", "ranking_blocked_by"),
    "notes": ("modification_hint",),
    "fallback_used": ("fallback_used",),
}


def _walk(source: Mapping[str, Any], key: str) -> Any:
    """Resolve a possibly dotted key against one mapping."""
    if "." not in key:
        return source.get(key)
    value: Any = source
    for part in key.split("."):
        if not isinstance(value, Mapping):
            return None
        value = value.get(part)
    return value


def _lookup(
    column: str,
    record: Mapping[str, Any],
    identity: Mapping[str, Any],
) -> Any:
    """Row-level value if the record has one, else the run-level identity."""
    keys = _ALIASES.get(column, (column,))
    for source in (record, identity):
        for key in keys:
            value = _walk(source, key)
            if value not in (None, "", [], {}):
                return value
    return NA


# ---------------------------------------------------------------------------
# Root tables.
# ---------------------------------------------------------------------------


#: Severity used only to pick which of a record's verdict fields to publish.
#: Matches the coordinator's consumer table (keep < not_evaluated < modify <
#: reject) rather than ``metric_scoring._combine_gate_statuses``, where
#: ``not_evaluated`` ranks *below* ``keep``.
_DECISION_SEVERITY = {"keep": 0, "not_evaluated": 1, "modify": 2, "reject": 3}


def _decision(record: Mapping[str, Any]) -> str:
    """Publish the molecule-level verdict, not just the per-category fold.

    ``metric_gate`` alone was the wrong source.  A structural reject
    (``valid_structure is False``, ``nomination_eligibility == "rejected"``) never
    enters the fold, so it leaves ``metric_gate`` unset and the decision column
    printed ``NA`` for a molecule the run had already rejected.  ``gate_status``
    alone is wrong in the other direction: it is combined under
    ``metric_scoring``'s convention where ``not_evaluated`` ranks below ``keep``,
    so a molecule with no metric evidence at all surfaces there as ``keep``.

    So all three are read and the most severe known label wins: ``classification``
    is the combined molecule verdict and normally decides, while ``metric_gate``
    and ``gate_status`` cover records written before it existed and guarantee that
    neither a structural nor a metric reject can be dropped.  ``not_evaluated``
    has no verdict of its own in the executable policy, and ``ol.decision`` maps it to
    ``NA`` -- which is the honest answer for a molecule nothing could measure, and
    keeps it out of ``KEEP`` (the executable policy forbids turning that into a rejection).
    """
    labels = [
        str(record.get(key) or "").strip().casefold()
        for key in ("classification", "metric_gate", "gate_status")
    ]
    known = [label for label in labels if label in _DECISION_SEVERITY]
    worst = max(known, key=_DECISION_SEVERITY.__getitem__) if known else None
    return ol.decision(worst, record.get("rank_status"))


def _final_row(
    record: Mapping[str, Any], identity: Mapping[str, Any],
) -> dict[str, Any]:
    row = {
        column: _lookup(column, record, identity)
        for column in ol.FINAL_REPORT_COLUMNS
    }
    row["decision"] = _decision(record)
    return row


def _full_row(
    record: Mapping[str, Any],
    identity: Mapping[str, Any],
    unimplemented: frozenset[str] = frozenset(),
) -> dict[str, Any]:
    row = {
        column: _lookup(column, record, identity)
        for column in ol.FINAL_REPORT_FULL_COLUMNS
    }
    row["decision"] = _decision(record)
    steps = record.get("step_scores") or {}
    for step_id in STEP_IDS:
        step = steps.get(step_id) or {}
        row[f"{step_id}_score"] = step.get("score", NA)
        row[f"{step_id}_rank"] = step.get("rank_avg", NA)
        row[f"{step_id}_coverage"] = step.get("category_weighted_coverage", NA)
    categories = record.get("category_scores") or {}
    for code in ol.CATEGORY_ORDER:
        category = categories.get(code) or {}
        row[f"{code}_score"] = category.get("score", NA)
        row[f"{code}_rank"] = category.get("rank_avg", NA)
        row[f"{code}_status"] = _category_status(category, unimplemented, code)
        row[f"{code}_coverage"] = category.get("planned_weighted_coverage", NA)
    return row


#: A gate-only category has no score to report, so its status carries the gate
#: verdict instead.  Prefixed, because this column also carries evidence-coverage
#: labels for the scored categories and the two must not be confused.
_GATE_STATUS: dict[str, str] = {
    "keep": "gate_pass",
    "modify": "gate_modify",
    "reject": "gate_reject",
}


def _category_status(
    category: Mapping[str, Any], unimplemented: frozenset[str] = frozenset(),
    code: str = "",
) -> str:
    """Why a category has the score it has -- or why it has none.

    The executable policy gives ALERT/BENCH/IDQC/TOPO/CHEM/DQC/FP the form ``gate/status``:
    their catalog weight is zero, so ``A_ic`` is zero, so they can never have a
    score and the gate verdict is the whole of their result.  Reporting them as
    ``not_planned`` claimed a QC step had been skipped when it had in fact run
    and passed, which is worse than reporting nothing.
    """
    if code and code in unimplemented:
        return "unimplemented"
    if not category:
        return "not_planned"
    if category.get("score") is None:
        gate = _GATE_STATUS.get(str(category.get("gate_status") or ""))
        if gate:
            return gate
        return "no_evidence" if category.get("planned") else "not_planned"
    if category.get("consequential_missing_metric_ids"):
        return "partial"
    return "complete"


#: Every flat table repeats this stem so any two of them can be joined.
JOIN_COLUMNS: tuple[str, ...] = (
    "row_id", "mol_name", "pose_id", "cohort_id",
)

ALL_METRICS_STEM: tuple[str, ...] = JOIN_COLUMNS + (
    "decision", "overall_score", "overall_rank",
    "aizynthfinder.status", "aizynthfinder.attempted",
    "aizynthfinder.reason",
)


def _ordered_metric_ids(rules: Sequence[Mapping[str, Any]]) -> list[str]:
    """Catalog order: step, then the category's directory order, then the id.

    A reader scrolling right through ``all_metrics_table.tsv`` therefore walks
    the same path as a reader clicking down through the directories.
    """
    def sort_key(rule: Mapping[str, Any]) -> tuple[Any, ...]:
        code = str(rule.get("category_code") or "")
        step_id = str(rule.get("step_id") or "")
        return (
            list(STEP_IDS).index(step_id) if step_id in STEP_DIR_KEYS else len(STEP_IDS),
            ol.CATEGORY_DIR.get(code, "zz"),
            str(rule.get("metric_id")),
        )

    return [str(rule["metric_id"]) for rule in sorted(rules, key=sort_key)]


def _stem(record: Mapping[str, Any], step_id: str | None = None) -> dict[str, Any]:
    return {
        "row_id": ol.row_id(record, step_id),
        "mol_name": record.get("mol_name", NA),
        "pose_id": record.get("pose_id", NA),
        "cohort_id": record.get("cohort_id", NA),
    }


def _synthesis_audit_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    """Expose why AiZynthFinder scalar cells are blank for this row."""
    property_record = record.get("property")
    source = property_record if isinstance(property_record, Mapping) else record
    synthesis = source.get("synthesis")
    if not isinstance(synthesis, Mapping):
        if "property" in record and property_record is None:
            property_join = record.get("property_join")
            reason = (
                property_join.get("reason")
                if isinstance(property_join, Mapping)
                else None
            )
            return {
                "aizynthfinder.status": "no_property_result",
                "aizynthfinder.attempted": NA,
                "aizynthfinder.reason": reason or "no property-filter row",
            }
        return {
            "aizynthfinder.status": NA,
            "aizynthfinder.attempted": NA,
            "aizynthfinder.reason": NA,
        }

    status = str(synthesis.get("status") or "missing")
    attempted = synthesis.get("attempted")
    if not isinstance(attempted, bool):
        if status in {"solved", "not_found_within_budget", "error"}:
            attempted = True
        elif status in {"disabled", "not_eligible", "not_selected_top_k"}:
            attempted = False
        else:
            attempted = NA
    return {
        "aizynthfinder.status": status,
        "aizynthfinder.attempted": attempted,
        "aizynthfinder.reason": (
            synthesis.get("reason") or synthesis.get("error") or NA
        ),
    }


def _has_tool_execution_evidence(
    records: Sequence[Mapping[str, Any]], tool_ids: Sequence[str],
) -> bool:
    """Whether an otherwise metric-empty tool run still needs publication."""
    if "aizynthfinder" not in tool_ids:
        return False
    for record in records:
        audit = _synthesis_audit_fields(record)
        if audit["aizynthfinder.attempted"] is True:
            return True
        if audit["aizynthfinder.status"] in {
            "solved", "not_found_within_budget", "error",
        }:
            return True
    return False


def _all_metrics_row(
    record: Mapping[str, Any], metric_ids: Sequence[str],
) -> dict[str, Any]:
    assessments = _assessments(record)
    row = _stem(record)
    row["decision"] = _decision(record)
    row["overall_score"] = record.get("overall_score", NA)
    row["overall_rank"] = record.get("overall_rank", NA)
    row.update(_synthesis_audit_fields(record))
    for metric_id in metric_ids:
        assessment = assessments.get(metric_id)
        row[metric_id] = (
            assessment.get("value") if _measured(assessment) else NA
        )
        row[f"{metric_id}__label"] = (
            assessment.get("label", NA) if _measured(assessment) else NA
        )
        row[f"{metric_id}__gate"] = (
            assessment.get("gate", NA) if _measured(assessment) else NA
        )
    return row


# ---------------------------------------------------------------------------
# Step layer.
# ---------------------------------------------------------------------------


def _step_table(
    records: Sequence[Mapping[str, Any]],
    step_id: str,
    metric_ids: Sequence[str],
) -> tuple[tuple[str, ...], list[dict[str, Any]]]:
    columns = JOIN_COLUMNS + (
        f"{step_id}_score", f"{step_id}_rank", f"{step_id}_coverage",
    ) + tuple(metric_ids)
    rows = []
    for record in records:
        assessments = _assessments(record)
        step = (record.get("step_scores") or {}).get(step_id) or {}
        row = _stem(record, step_id)
        row[f"{step_id}_score"] = step.get("score", NA)
        row[f"{step_id}_rank"] = step.get("rank_avg", NA)
        row[f"{step_id}_coverage"] = step.get("category_weighted_coverage", NA)
        for metric_id in metric_ids:
            assessment = assessments.get(metric_id)
            row[metric_id] = (
                assessment.get("value") if _measured(assessment) else NA
            )
        rows.append(row)
    return columns, rows


def _step_report_rows(
    records: Sequence[Mapping[str, Any]], step_id: str,
) -> list[dict[str, Any]]:
    rows = []
    for record in records:
        step = (record.get("step_scores") or {}).get(step_id) or {}
        rows.append({
            "row_id": ol.row_id(record, step_id),
            "step_id": step_id,
            "step_score": step.get("score", NA),
            "step_rank": step.get("rank_avg", NA),
            "category_weighted_coverage": step.get(
                "category_weighted_coverage", NA,
            ),
            "planned_category_weighted_coverage": step.get(
                "planned_category_weighted_coverage", NA,
            ),
            "catalog_category_weighted_coverage": step.get(
                "catalog_category_weighted_coverage", NA,
            ),
            "available_parent_weight": step.get("available_parent_weight", NA),
            "reference_parent_weight": step.get("reference_parent_weight", NA),
            "available_categories": step.get("available_categories", ()),
            "missing_reference_categories": step.get(
                "missing_reference_categories", (),
            ),
            "rank_status": record.get("rank_status", NA),
            "gate_status": record.get("gate_status", NA),
        })
    return rows


# ---------------------------------------------------------------------------
# Category layer.
# ---------------------------------------------------------------------------

#: Per-metric drill-down columns inside a category table.  ``__label`` is the
#: canonical pass/flag/fail mark; ``__gate`` remains as a legacy
#: keep/modify/reject alias for readers of older publication trees.
_METRIC_SUFFIXES: tuple[str, ...] = (
    "", "__label", "__gate", "__percentile",
)


def _category_table(
    records: Sequence[Mapping[str, Any]],
    code: str,
    metric_ids: Sequence[str],
    step_id: str,
) -> tuple[tuple[str, ...], list[dict[str, Any]]]:
    """The category's flat table: its score, then the evidence under it.

    Each metric contributes its raw value, the gate that value tripped and the
    cohort percentile the score actually used, because a reviewer asking "why is
    this category 42?" needs all three and no fourth thing.
    """
    columns = JOIN_COLUMNS + (
        "category_score", "category_rank", "category_status", "category_coverage",
    ) + tuple(
        f"{metric_id}{suffix}"
        for metric_id in metric_ids for suffix in _METRIC_SUFFIXES
    )
    rows = []
    for record in records:
        assessments = _assessments(record)
        category = (record.get("category_scores") or {}).get(code) or {}
        row = _stem(record, step_id)
        row["category_score"] = category.get("score", NA)
        row["category_rank"] = category.get("rank_avg", NA)
        row["category_status"] = _category_status(category, frozenset(), "")
        row["category_coverage"] = category.get("planned_weighted_coverage", NA)
        for metric_id in metric_ids:
            assessment = assessments.get(metric_id)
            if _measured(assessment):
                row[metric_id] = assessment.get("value")
                row[f"{metric_id}__label"] = assessment.get("label", NA)
                row[f"{metric_id}__gate"] = assessment.get("gate", NA)
                row[f"{metric_id}__percentile"] = assessment.get(
                    "rank_percentage", NA,
                )
            else:
                for suffix in _METRIC_SUFFIXES:
                    row[f"{metric_id}{suffix}"] = NA
        rows.append(row)
    return columns, rows


def _category_report_rows(
    records: Sequence[Mapping[str, Any]], code: str,
) -> list[dict[str, Any]]:
    config = CATEGORY_CONFIG[code]
    step_id = config["step_id"]
    rows = []
    for record in records:
        category = (record.get("category_scores") or {}).get(code) or {}
        rows.append({
            "row_id": ol.row_id(record, step_id),
            "category_code": code,
            "category_name": ol.CATEGORY_DIR[code],
            "step_id": step_id,
            "parent_weight": config["parent_weight"],
            "score": category.get("score", NA),
            "rank": category.get("rank_avg", NA),
            "status": _category_status(category, frozenset(), ""),
            "coverage": category.get("planned_weighted_coverage", NA),
            "available_abs_weight": category.get("available_abs_weight", NA),
            "planned_abs_weight": category.get("planned_abs_weight", NA),
            "catalog_abs_weight": category.get(
                "catalog_abs_weight", config["catalog_abs_weight"],
            ),
            "catalog_implementation_fraction": category.get(
                "catalog_implementation_fraction", NA,
            ),
            "evidence_strength": category.get("evidence_strength", NA),
            "available_metric_ids": category.get("available_metric_ids", ()),
            "consequential_missing_metric_ids": category.get(
                "consequential_missing_metric_ids", (),
            ),
            "audit_only_missing_metric_ids": category.get(
                "audit_only_missing_metric_ids", (),
            ),
        })
    return rows


# ---------------------------------------------------------------------------
# Tool layer.  The raw output stays where the tool wrote it; what is published
# here is the tool's metric columns plus a pointer back to that output.
# ---------------------------------------------------------------------------


def _tool_table(
    records: Sequence[Mapping[str, Any]],
    metric_ids: Sequence[str],
    step_id: str,
) -> tuple[tuple[str, ...], list[dict[str, Any]]]:
    columns = JOIN_COLUMNS + tuple(
        f"{metric_id}{suffix}"
        for metric_id in metric_ids
        for suffix in ("", "__label", "__gate")
    )
    rows = []
    for record in records:
        assessments = _assessments(record)
        row = _stem(record, step_id)
        for metric_id in metric_ids:
            assessment = assessments.get(metric_id)
            measured = _measured(assessment)
            row[metric_id] = assessment.get("value") if measured else NA
            row[f"{metric_id}__label"] = (
                assessment.get("label", NA) if measured else NA
            )
            row[f"{metric_id}__gate"] = (
                assessment.get("gate", NA) if measured else NA
            )
        rows.append(row)
    return columns, rows


def _tool_summary(
    records: Sequence[Mapping[str, Any]],
    rules: Sequence[Mapping[str, Any]],
    tool_run_id: str,
    tool_ids: Sequence[str],
    step_id: str,
) -> dict[str, Any]:
    metrics = []
    for rule in rules:
        metric_id = str(rule["metric_id"])
        availability_statuses = Counter()
        for record in records:
            assessment = _assessments(record).get(metric_id)
            if not assessment:
                availability_statuses["missing_assessment"] += 1
                continue
            status = str(assessment.get("availability_status") or "unknown")
            availability_statuses[status] += 1
        measured = availability_statuses.get("available", 0)
        missing = len(records) - measured
        not_planned = availability_statuses.get("not_planned", 0)
        planned_missing = availability_statuses.get("planned_missing", 0)
        metrics.append({
            "metric_id": metric_id,
            "output_metric": rule.get("output_metric"),
            "category_code": rule.get("category_code"),
            "weight": rule.get("weight"),
            "signed_weight": rule.get("signed_weight"),
            "rows_measured": measured,
            "rows_unmeasured": missing,
            "rows_missing": missing,
            "rows_missing_deprecated_alias_of": "rows_unmeasured",
            "rows_not_planned": not_planned,
            "rows_planned_missing": planned_missing,
            "rows_other_unavailable": missing - not_planned - planned_missing,
            "availability_status_counts": dict(sorted(availability_statuses.items())),
            "rows_missing_semantics": (
                "coverage absence, not an execution-failure count; inspect "
                "availability_status_counts and execution_summary"
            ),
        })
    logical_ids = sorted(set(tool_ids))
    tool_names = sorted({
        str(rule.get("tool")) for rule in rules if rule.get("tool")
    })
    summary = {
        "tool_run_id": tool_run_id,
        "tool_id": logical_ids[0] if len(logical_ids) == 1 else tool_run_id,
        "tool_ids": logical_ids,
        "tool_name": tool_names[0] if len(tool_names) == 1 else tool_run_id,
        "tool_names": tool_names,
        "step_id": step_id,
        "rows": len(records),
        "row_id_field": ol.ROW_ID_FIELD[step_id],
        "metrics": sorted(metrics, key=lambda item: item["metric_id"]),
    }
    if "aizynthfinder" in logical_ids:
        attempts = Counter()
        statuses = Counter()
        reasons = Counter()
        search_completed = 0
        operational_errors = 0
        contract_inconsistent = 0
        for record in records:
            audit = _synthesis_audit_fields(record)
            attempted = audit["aizynthfinder.attempted"]
            if attempted is True:
                attempts["attempted"] += 1
            elif attempted is False:
                attempts["not_attempted"] += 1
            else:
                attempts["unknown"] += 1
            raw_status = audit["aizynthfinder.status"]
            status = str(raw_status).strip() if raw_status is not None else ""
            status = status or "unknown"
            statuses[status] += 1
            raw_reason = audit["aizynthfinder.reason"]
            reason = str(raw_reason).strip() if raw_reason is not None else ""
            reasons[reason or "unspecified"] += 1
            expected_attempted = None
            if status in {"solved", "not_found_within_budget", "error"}:
                expected_attempted = True
            elif status in {
                "disabled", "not_eligible", "not_selected_top_k",
            }:
                expected_attempted = False
            if (
                expected_attempted is not None
                and attempted is not expected_attempted
            ):
                contract_inconsistent += 1
            if attempted is True and status in {
                "solved", "not_found_within_budget",
            }:
                search_completed += 1
            if attempted is True and status == "error":
                operational_errors += 1
        summary["execution_summary"] = {
            "rows_attempted": attempts["attempted"],
            "rows_not_attempted": attempts["not_attempted"],
            "rows_attempt_status_unknown": attempts["unknown"],
            "rows_search_completed": search_completed,
            "rows_operational_errors": operational_errors,
            "rows_status_error": statuses["error"],
            "contract_inconsistent_rows": contract_inconsistent,
            "status_counts": dict(sorted(statuses.items())),
            "reason_counts": dict(sorted(reasons.items())),
            "failure_count_semantics": (
                "only consistent attempted=true,status=error rows are "
                "operational failures; rows_unmeasured, deprecated "
                "rows_missing, and molecular reject counts are not tool failures"
            ),
        }
    return summary


# ---------------------------------------------------------------------------
# metric_catalog_and_scoring.md -- frozen alongside the tables it explains.
# ---------------------------------------------------------------------------

_FORMULAS = """\
## Scoring

Gating runs on the raw value, before any normalization.  Level 1
marks one metric: a value outside a hard bound is `fail`; inside a hard bound
but outside a soft bound is `flag`; inside both soft bounds is `pass`.  All four
bounds are inclusive.  A missing value is `not_applicable` -- it is never
treated as `0` and never enters the fail-ratio denominator.

Level 2 is the only molecule decision.  Within each parent category calculate
`count(fail)` and `ratio(fail) = fail / (fail + flag + pass)`.  V1 never
eliminates; V2 rejects at `count(fail) >= 3` or `ratio(fail) >= 0.50`; V3
rejects at `count(fail) >= 2` or `ratio(fail) >= 0.30`. Otherwise every
evaluated category needs strictly more than half its labels to be `pass` for
`keep`; other labeled cases give `modify`.

    RP_im       = (n - RANK.AVG(x_im, desc)) / (n - 1)      # 0.5 when n = 1
    Score_ic    = 100 x (Raw_ic + B-_ic) / A_ic             # NA when A_ic = 0
    StepScore_is= 100 x sum(W'_c x RP_ic) / sum(W'_c)       # available c only
    OverallScore_i = sum(APW_is x StepScore_is) / sum(APW_is)

`A_ic` is the available absolute weight of category `c`, `B-_ic` the absolute
weight of its "lower is better" metrics, and `APW_is` the parent weight step `s`
could actually use.  A metric with no value leaves both sides of every ratio,
so a missing measurement shrinks the denominator instead of scoring zero.

## Release status

`category_weighted_coverage` divides the available parent weight by the parent
weight this round planned, so a category nobody planned leaves both
sides of the ratio rather than charging every molecule for it.  At or above 0.50 a
run is `formal`, from 0.10 to 0.50 `provisional`, below that `NA` plus a
backtrack request.  Coverage never rejects a molecule or a pose.
`catalog_category_weighted_coverage` reports the same numerator over the full
section 3 catalog; it is a pipeline-local audit ratio that says how much of the
published catalog this profile reaches and never decides a band.
Independently of either ratio, any active step with no score at all caps the
status at `provisional`: high coverage elsewhere cannot stand in for a step that
was never measured.

## Decision

`REJECT` and `MODIFY` come straight from the gate.  A row that cleared every
gate is `KEEP` when its release status is `formal` and `RESERVE` otherwise --
held back for thin evidence rather than for anything it did wrong.

## Category status

`<CODE>_status` answers two different questions, because two kinds of category
exist.  A scored category reports how complete its evidence was: `complete`,
`partial` (a weighted metric is missing), `no_evidence`, `not_planned`, or
`unimplemented` (the specification lists this category but no tool in this build
can compute it).  A gate-only category carries no child weight, so `A_ic` is zero
and it can never have a score; it reports the worst gate its members tripped
instead -- `gate_pass`, `gate_modify` or `gate_reject`.

## Deliberately empty columns

`uncertainty_grade` has no producer yet.  It asks how precise the numbers are,
which is not what `coverage_status` answers, so it stays blank rather than
repeating the coverage judgement under a second heading.
"""


def _catalog_markdown(
    rules: Sequence[Mapping[str, Any]],
    metadata: Mapping[str, Any],
    profile: Mapping[str, Any],
) -> str:
    lines = [
        "# Metric catalog and scoring profile",
        "",
        "Frozen with the tables in this directory: the scores here were produced",
        "by exactly these rules.  Re-running with a different catalog produces a",
        "different directory, never an edit to this one.",
        "",
        "## Profile",
        "",
        "| field | value |",
        "| --- | --- |",
    ]
    for key, value in sorted(profile.items()):
        lines.append(f"| {key} | {ol.format_cell(value)} |")
    lines += ["", _FORMULAS, "", "## Categories", "",
              "| code | step | directory | parent weight | catalog weight |",
              "| --- | --- | --- | --- | --- |"]
    for code in ol.CATEGORY_ORDER:
        config = CATEGORY_CONFIG[code]
        lines.append(
            f"| {code} | {config['step_id']} | {ol.CATEGORY_DIR[code]} "
            f"| {ol.format_cell(config['parent_weight'])} "
            f"| {ol.format_cell(config['catalog_abs_weight'])} |"
        )
    unimplemented = metadata.get("unimplemented_reference_categories") or []
    if unimplemented:
        lines += [
            "",
            "Referenced by section 3 but with no executable rule in this build, "
            "so never scored: " + ", ".join(sorted(unimplemented)) + ".",
        ]
    lines += ["", "## Metrics", "",
              "| metric_id | tool | category | step | weight | direction "
              "| hard_lower | soft_lower | soft_upper | hard_upper | transform |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for metric_id in _ordered_metric_ids(rules):
        rule = next(item for item in rules if str(item["metric_id"]) == metric_id)
        lines.append("| " + " | ".join(ol.format_cell(cell) for cell in (
            metric_id, rule.get("tool"), rule.get("category_code"),
            rule.get("step_id"), rule.get("signed_weight"), rule.get("direction"),
            rule.get("hard_lower"), rule.get("soft_lower"),
            rule.get("soft_upper"), rule.get("hard_upper"),
            rule.get("rank_transform"),
        )) + " |")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The entry point.
# ---------------------------------------------------------------------------


def _ordered_records(
    records: Sequence[Mapping[str, Any]],
) -> list[Mapping[str, Any]]:
    """Best first, grouped by partition; unranked rows trail their own group.

    ``overall_rank`` is local to one ``(cohort_id, normalization_pool)``
    partition, so sorting on it alone interleaves partitions and can put a
    benchmark decoy that ranked 1st in the decoy pool above the best nomination
    candidate. The partition leads the key, with ``nomination_candidates`` first,
    so the top row of the table is always a real candidate.
    """
    def key(item: tuple[int, Mapping[str, Any]]) -> tuple[Any, ...]:
        index, record = item
        rank = record.get("overall_rank")
        pool = str(record.get("normalization_pool") or "")
        return (
            0 if pool in ("", "nomination_candidates") else 1,
            pool,
            str(record.get("cohort_id") or ""),
            0 if rank is not None else 1,
            float(rank) if rank is not None else 0.0,
            index,
        )

    return [record for _index, record in sorted(enumerate(records), key=key)]


def _active_steps(metadata: Mapping[str, Any]) -> tuple[str, ...]:
    declared = metadata.get("active_steps") or STEP_IDS
    return tuple(step for step in STEP_IDS if step in set(declared))


#: Scorer facts that are identical in every normalization pool, so they belong to
#: the run rather than to a cohort.
_RUN_LEVEL_SCORING_FIELDS: tuple[str, ...] = (
    "schema_version", "active_steps", "active_steps_source",
    "coverage_release_thresholds", "reference_step_parent_weight",
    "unimplemented_reference_categories", "normalization_method",
    "normalization_scope", "family_weight_caps", "category_config",
)


def flatten_scoring_metadata(metadata: Mapping[str, Any]) -> dict[str, Any]:
    """Accept either a raw ``score_records`` metadata block or the coordinator's.

    The coordinator scores each cohort and each normalization pool separately, so
    its envelope keeps only the fields it owns and files the scorer's own metadata
    under ``cohorts[cohort][pool]``.  The run-level facts there are identical
    across pools by construction, so lifting them costs nothing and lets one
    publisher serve both shapes.
    """
    flat = dict(metadata)
    for cohort in (metadata.get("cohorts") or {}).values():
        if not isinstance(cohort, Mapping):
            continue
        for pool in cohort.values():
            if not isinstance(pool, Mapping):
                continue
            for key in _RUN_LEVEL_SCORING_FIELDS:
                if flat.get(key) in (None, [], {}, ()) and key in pool:
                    flat[key] = pool[key]
    return flat


def _default_profile(metadata: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "scoring_schema_version": metadata.get("schema_version"),
        "metric_rules_version": metadata.get("metric_rules_version"),
        "normalization_method": metadata.get("normalization_method"),
        "normalization_scope": metadata.get("normalization_scope"),
        "active_steps": metadata.get("active_steps"),
        "active_steps_source": metadata.get("active_steps_source"),
        "coverage_release_thresholds": metadata.get(
            "coverage_release_thresholds",
        ),
        "layout_contract": ol.LAYOUT_CONTRACT,
    }


# ---------------------------------------------------------------------------
# Workflow-stage evidence indexes.  These make the early pipeline stages
# discoverable without pretending they produced candidate-level score tables.
# Raw artifacts remain in the resumable working tree and are referenced only.
# ---------------------------------------------------------------------------


def _read_json_object(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _require_completed_pipeline_sources(working_root: Path) -> None:
    """Refuse to seal a pipeline view before coordination and reporting end."""
    incomplete = []
    for relative in _FINALIZATION_RESULT_FILES:
        result = _read_json_object(working_root / relative)
        if not isinstance(result, dict):
            incomplete.append(f"{relative}: missing_or_invalid")
        elif result.get("success") is not True:
            incomplete.append(f"{relative}: not_successful")
    if incomplete:
        raise ValueError(
            "pipeline publication requires completed coordination and reporting: "
            + "; ".join(incomplete)
        )


def _resolve_source_path(
    value: Any,
    *,
    working_root: Path,
    result_path: Path | None = None,
) -> Path | None:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        return None
    raw = str(value).strip()
    if "://" in raw:
        return None
    candidate = Path(raw).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()
    choices = []
    if result_path is not None:
        choices.append((result_path.parent / candidate).resolve())
    choices.append((working_root / candidate).resolve())
    for choice in choices:
        if choice.exists():
            return choice
    return choices[-1]


def _artifact_row(
    *,
    layout: ol.RunLayout,
    pipeline_step_id: str,
    declared_by: str,
    artifact_name: str,
    path: Path,
    kind: str = "path",
    producer: str = "",
) -> dict[str, Any]:
    resolved_kind = "directory" if path.is_dir() else kind
    digest = (
        ol.directory_digest(path)
        if resolved_kind == "directory"
        else ol.file_digest(path)
    )
    return {
        "pipeline_step_id": pipeline_step_id,
        "declared_by": declared_by,
        "artifact_name": artifact_name,
        "kind": resolved_kind,
        "producer": producer,
        "path": layout.relative(path),
        "exists": digest.get("exists", False),
        "bytes": digest.get("bytes", NA),
        "sha256": digest.get("sha256", NA),
    }


def _validate_recorded_path_identity(
    *,
    section: str,
    artifact_name: str,
    record: Mapping[str, Any],
    path: Path,
) -> None:
    """Reject a path whose live bytes differ from the step's lineage record."""
    if path.is_symlink() or not path.is_file():
        raise ValueError(
            f"{section} artifact {artifact_name!r} is missing or not a file: {path}"
        )
    fingerprint = str(record.get("fingerprint") or "")
    if not fingerprint.startswith("sha256:") or len(fingerprint) != 71:
        raise ValueError(
            f"{section} artifact {artifact_name!r} has no verifiable SHA-256 fingerprint"
        )
    try:
        int(fingerprint.removeprefix("sha256:"), 16)
    except ValueError as exc:
        raise ValueError(
            f"{section} artifact {artifact_name!r} has an invalid SHA-256 fingerprint"
        ) from exc
    live = ol.file_digest(path)
    if fingerprint.removeprefix("sha256:") != live.get("sha256"):
        raise ValueError(
            f"{section} artifact {artifact_name!r} fingerprint changed"
        )
    recorded_size = record.get("size")
    alternate_size = record.get("size_bytes")
    if recorded_size is not None and alternate_size is not None:
        try:
            sizes_agree = int(recorded_size) == int(alternate_size)
        except (TypeError, ValueError):
            sizes_agree = False
        if not sizes_agree:
            raise ValueError(
                f"{section} artifact {artifact_name!r} has conflicting recorded sizes"
            )
    size_value = recorded_size if recorded_size is not None else alternate_size
    if size_value is not None:
        try:
            expected_size = int(size_value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{section} artifact {artifact_name!r} has an invalid recorded size"
            ) from exc
        if isinstance(size_value, bool) or expected_size < 0:
            raise ValueError(
                f"{section} artifact {artifact_name!r} has an invalid recorded size"
            )
        if expected_size != live.get("bytes"):
            raise ValueError(
                f"{section} artifact {artifact_name!r} size changed"
            )


def _looks_like_file_reference(value: str) -> bool:
    candidate = value.strip()
    if not candidate or "://" in candidate:
        return False
    return (
        "/" in candidate
        or "\\" in candidate
        or Path(candidate).suffix.lower() in _PATH_SUFFIXES
    )


def _declared_file_leaves(
    value: Any,
    prefix: str = "",
) -> Iterable[tuple[str, str]]:
    if isinstance(value, Mapping):
        for key, child in value.items():
            name = f"{prefix}.{key}" if prefix else str(key)
            yield from _declared_file_leaves(child, name)
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            name = f"{prefix}[{index}]"
            yield from _declared_file_leaves(child, name)
    elif isinstance(value, (str, Path)) and _looks_like_file_reference(str(value)):
        yield prefix or "file", str(value)


def _report_bundle_artifact_rows(
    *,
    layout: ol.RunLayout,
    working_root: Path,
    result_path: Path,
    result: Mapping[str, Any],
    seen_paths: set[str],
) -> list[dict[str, Any]]:
    """Validate and index every member sealed by the reporting result.

    Report HTML depends on plots and shared JavaScript that are not individual
    top-level ``files`` entries.  Re-hashing only ``result.json`` and the HTML
    would therefore miss a damaged image or asset.  The report-owned manifest
    is treated as an index, not as trust: every member is independently checked
    here and later included in the pipeline source snapshot.
    """
    files = result.get("files") or {}
    if not isinstance(files, Mapping):
        raise ValueError("successful reporting result has no files object")
    manifest_path = _resolve_source_path(
        files.get("report_artifact_manifest_json"),
        working_root=working_root,
        result_path=result_path,
    )
    if manifest_path is None or not manifest_path.is_file():
        raise ValueError(
            "successful reporting result has no report artifact manifest"
        )
    manifest = _read_json_object(manifest_path)
    if not isinstance(manifest, Mapping) or manifest.get("schema_version") != 1:
        raise ValueError("report artifact manifest is missing or invalid")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ValueError("report artifact manifest has no artifacts")

    required = {
        "report_html": files.get("report_html"),
        "report_data": files.get("report_data_json"),
        "assets_manifest": files.get("assets_manifest_json"),
    }
    if files.get("pymol_cover_png"):
        required["pymol_cover"] = files.get("pymol_cover_png")
    expected_paths: dict[str, Path] = {}
    for role, raw_path in required.items():
        path = _resolve_source_path(
            raw_path,
            working_root=working_root,
            result_path=result_path,
        )
        if path is None:
            raise ValueError(f"reporting result is missing required file {role}")
        expected_paths[role] = path

    assets_manifest = expected_paths["assets_manifest"]
    assets_payload = _read_json_object(assets_manifest)
    assets = assets_payload.get("assets") if isinstance(assets_payload, Mapping) else None
    if not isinstance(assets, Mapping):
        raise ValueError("report assets manifest is missing or invalid")
    for filename, record in sorted(assets.items()):
        if not isinstance(record, Mapping):
            raise ValueError(
                f"report assets manifest record is invalid: {filename!r}"
            )
        if record.get("present") is not True:
            continue
        role = f"asset:{filename}"
        path = _resolve_source_path(
            record.get("path"),
            working_root=working_root,
            result_path=assets_manifest,
        )
        if path is None or role in expected_paths:
            raise ValueError(f"report asset declaration is ambiguous: {filename!r}")
        expected_paths[role] = path

    plots_status = result.get("plots_status") or {}
    figures = plots_status.get("figures") if isinstance(plots_status, Mapping) else None
    if figures is not None and not isinstance(figures, Mapping):
        raise ValueError("report plots_status.figures is invalid")
    for name, raw_path in sorted((figures or {}).items()):
        if not raw_path:
            continue
        role = f"plot:{name}"
        path = _resolve_source_path(
            raw_path,
            working_root=working_root,
            result_path=manifest_path,
        )
        if path is None or role in expected_paths:
            raise ValueError(f"report plot declaration is ambiguous: {name!r}")
        expected_paths[role] = path

    rows: list[dict[str, Any]] = []
    roles: dict[str, Path] = {}
    member_paths: set[Path] = set()
    for index, record in enumerate(artifacts):
        if not isinstance(record, Mapping):
            raise ValueError(
                f"report artifact manifest member {index} is invalid"
            )
        role = str(record.get("role") or "").strip()
        path = _resolve_source_path(
            record.get("path"),
            working_root=working_root,
            result_path=manifest_path,
        )
        if not role or path is None or role in roles or path in member_paths:
            raise ValueError(
                f"report artifact manifest member {index} is ambiguous"
            )
        row = _artifact_row(
            layout=layout,
            pipeline_step_id="reporting",
            declared_by="report_artifact_manifest",
            artifact_name=role,
            path=path,
            producer="reporting",
        )
        try:
            expected_bytes = int(record.get("size_bytes"))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"report artifact manifest member {role!r} has invalid size"
            ) from exc
        expected_sha256 = str(record.get("sha256") or "")
        if (
            row.get("exists") is not True
            or row.get("bytes") != expected_bytes
            or row.get("sha256") != expected_sha256
        ):
            raise ValueError(
                f"report artifact manifest member changed: {role}"
            )
        roles[role] = path
        member_paths.add(path)
        if row["path"] not in seen_paths:
            rows.append(row)
            seen_paths.add(str(row["path"]))

    if set(roles) != set(expected_paths):
        missing = sorted(set(expected_paths) - set(roles))
        extra = sorted(set(roles) - set(expected_paths))
        raise ValueError(
            "report artifact manifest membership disagrees with report declarations: "
            f"missing={missing}; extra={extra}"
        )
    for role, expected_path in expected_paths.items():
        if roles.get(role) != expected_path:
            raise ValueError(
                f"report artifact manifest does not seal result file {role}"
            )
    return rows


def _runtime_step_evidence(
    layout: ol.RunLayout,
    working_root: Path | None,
    legacy_step: str,
    *,
    strict_lineage: bool = False,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pipeline_step_id = _PIPELINE_STEP_NAMES[legacy_step]
    if working_root is None:
        return ({
            "pipeline_step_id": pipeline_step_id,
            "working_directory": None,
            "result_json": None,
            "status": "source_root_not_supplied",
        }, [])

    result_path = working_root / legacy_step / "result.json"
    relative_result = layout.relative(result_path)
    if not result_path.is_file():
        return ({
            "pipeline_step_id": pipeline_step_id,
            "working_directory": layout.relative(result_path.parent),
            "result_json": relative_result,
            "status": "not_run",
        }, [])

    rows = [_artifact_row(
        layout=layout,
        pipeline_step_id=pipeline_step_id,
        declared_by="result_record",
        artifact_name="result_json",
        path=result_path,
        producer=pipeline_step_id,
    )]
    result = _read_json_object(result_path)
    if result is None:
        return ({
            "pipeline_step_id": pipeline_step_id,
            "working_directory": layout.relative(result_path.parent),
            "result_json": relative_result,
            "status": "invalid_result_json",
        }, rows)

    success = result.get("success")
    status = "complete" if success is True else "failed" if success is False else "unknown"
    seen_paths = {rows[0]["path"]}
    for section in ("artifact_inputs", "artifact_outputs"):
        artifacts = result.get(section) or {}
        if not isinstance(artifacts, Mapping):
            if strict_lineage:
                raise ValueError(
                    f"successful {pipeline_step_id} result has invalid {section}"
                )
            continue
        for artifact_name, raw_record in sorted(artifacts.items()):
            if not isinstance(raw_record, Mapping):
                if strict_lineage:
                    raise ValueError(
                        f"{pipeline_step_id} {section} record {artifact_name!r} is invalid"
                    )
                continue
            kind = str(raw_record.get("kind") or "path")
            if kind == "inline":
                continue
            if not raw_record.get("path"):
                if strict_lineage:
                    raise ValueError(
                        f"{pipeline_step_id} {section} path is missing for "
                        f"{artifact_name!r}"
                    )
                continue
            path = _resolve_source_path(
                raw_record.get("path"),
                working_root=working_root,
                result_path=result_path,
            )
            if path is None:
                if strict_lineage:
                    raise ValueError(
                        f"{pipeline_step_id} {section} path is invalid for "
                        f"{artifact_name!r}"
                    )
                continue
            if strict_lineage:
                _validate_recorded_path_identity(
                    section=f"{pipeline_step_id} {section}",
                    artifact_name=str(artifact_name),
                    record=raw_record,
                    path=path,
                )
            relative = layout.relative(path)
            if relative in seen_paths:
                continue
            rows.append(_artifact_row(
                layout=layout,
                pipeline_step_id=pipeline_step_id,
                declared_by=section,
                artifact_name=str(artifact_name),
                path=path,
                kind=kind,
                producer=str(raw_record.get("producer") or ""),
            ))
            seen_paths.add(relative)

    files = result.get("files") or {}
    if isinstance(files, Mapping):
        visible_files = {
            key: value for key, value in files.items()
            if not str(key).startswith("layout_")
            and key not in {"published_dir", "user_output_dir"}
        }
        for artifact_name, raw_path in _declared_file_leaves(visible_files):
            path = _resolve_source_path(
                raw_path,
                working_root=working_root,
                result_path=result_path,
            )
            if path is None or (not path.is_file() and not path.is_dir()):
                if strict_lineage:
                    raise ValueError(
                        f"{pipeline_step_id} declared file {artifact_name!r} "
                        f"is missing: {raw_path}"
                    )
                continue
            relative = layout.relative(path)
            if relative in seen_paths:
                continue
            rows.append(_artifact_row(
                layout=layout,
                pipeline_step_id=pipeline_step_id,
                declared_by="files",
                artifact_name=artifact_name,
                path=path,
                producer=pipeline_step_id,
            ))
            seen_paths.add(relative)

    if pipeline_step_id == "reporting" and success is True:
        rows.extend(_report_bundle_artifact_rows(
            layout=layout,
            working_root=working_root,
            result_path=result_path,
            result=result,
            seen_paths=seen_paths,
        ))

    return ({
        "pipeline_step_id": pipeline_step_id,
        "working_directory": layout.relative(result_path.parent),
        "result_json": relative_result,
        "status": status,
        "artifact_inputs": len(result.get("artifact_inputs") or {}),
        "artifact_outputs": len(result.get("artifact_outputs") or {}),
        "indexed_files": len(rows),
    }, rows)


def _stage_status(source_results: Sequence[Mapping[str, Any]]) -> str:
    statuses = {str(item.get("status") or "unknown") for item in source_results}
    if statuses == {"complete"}:
        return "complete"
    if statuses and statuses <= {"not_run", "source_root_not_supplied"}:
        return "not_run"
    if "failed" in statuses or "invalid_result_json" in statuses:
        return "failed" if len(statuses) == 1 else "partial"
    return "partial" if statuses else "not_run"


def _publish_stage_evidence(
    layout: ol.RunLayout,
    *,
    directory: Path,
    publication_stage_id: str,
    publication_stage_name: str,
    legacy_steps: Sequence[str],
    working_root: Path | None,
    strict_lineage: bool = False,
    manifest_name: str = ol.STEP_MANIFEST,
    manifest_extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    source_results = []
    rows = []
    for legacy_step in legacy_steps:
        summary, step_rows = _runtime_step_evidence(
            layout,
            working_root,
            legacy_step,
            strict_lineage=strict_lineage,
        )
        source_results.append(summary)
        rows.extend(step_rows)
    cited = [
        str(item["result_json"])
        for item in source_results if item.get("result_json")
    ]
    index_path = layout.publish_tsv(
        directory / ol.ARTIFACT_INDEX,
        _ARTIFACT_INDEX_COLUMNS,
        rows,
        role="artifact_index",
        derived_from=cited,
    )
    payload = {
        "manifest_type": "publication_stage",
        "path_base": "published_root",
        "publication_stage_id": publication_stage_id,
        "publication_stage_name": publication_stage_name,
        "status": _stage_status(source_results),
        "generated_at": layout.generated_at,
        "raw_artifact_policy": "reference_only",
        "artifact_index": layout.relative(index_path),
        "pipeline_step_ids": [
            _PIPELINE_STEP_NAMES[legacy_step] for legacy_step in legacy_steps
        ],
        "source_results": source_results,
        **dict(manifest_extra or {}),
    }
    layout.publish_json(
        directory / manifest_name,
        payload,
        role="stage_manifest",
        derived_from=[layout.relative(index_path)],
    )
    return payload


def _publish_workflow_evidence(
    layout: ol.RunLayout,
    working_root: Path | None,
    *,
    strict_lineage: bool = False,
) -> dict[str, dict[str, Any]]:
    stage_manifests: dict[str, dict[str, Any]] = {}

    input_rows = []
    if working_root is not None:
        for name in _CONTROL_INPUT_FILES:
            path = working_root / name
            if path.is_file():
                input_rows.append(_artifact_row(
                    layout=layout,
                    pipeline_step_id="pipeline",
                    declared_by="run_control",
                    artifact_name=Path(name).stem,
                    path=path,
                    producer="pipeline_entry",
                ))
            elif strict_lineage:
                raise ValueError(
                    f"pipeline-final control input is missing: {path}"
                )
    input_index = layout.publish_tsv(
        layout.path(ol.INPUTS_DIR, ol.ARTIFACT_INDEX),
        _ARTIFACT_INDEX_COLUMNS,
        input_rows,
        role="artifact_index",
        derived_from=[str(row["path"]) for row in input_rows],
    )
    input_payload: dict[str, Any] = {
        "manifest_type": "publication_inputs",
        "path_base": "published_root",
        "publication_stage_id": "00",
        "publication_stage_name": "inputs",
        "status": "complete" if input_rows else "no_local_control_files",
        "generated_at": layout.generated_at,
        "raw_artifact_policy": "reference_only",
        "artifact_index": layout.relative(input_index),
        "source_files": [str(row["path"]) for row in input_rows],
    }
    if working_root is not None:
        config = _read_json_object(working_root / "pipeline_config.json") or {}
        input_payload["project"] = config.get("project") or {}
        input_payload["context"] = config.get("context") or {}
    layout.publish_json(
        layout.path(ol.INPUTS_DIR, ol.INPUT_MANIFEST),
        input_payload,
        role="input_manifest",
        derived_from=[layout.relative(input_index)],
    )
    stage_manifests[ol.INPUTS_DIR] = input_payload

    for directory_name, legacy_steps in _WORKFLOW_STAGE_SOURCES.items():
        stage_id, stage_name = directory_name.split("_", 1)
        stage_manifests[directory_name] = _publish_stage_evidence(
            layout,
            directory=layout.path(directory_name),
            publication_stage_id=stage_id,
            publication_stage_name=stage_name,
            legacy_steps=legacy_steps,
            working_root=working_root,
            strict_lineage=strict_lineage,
        )

    _publish_stage_evidence(
        layout,
        directory=layout.path(ol.SM_PROPERTY_DIR, "00_prefilter"),
        publication_stage_id="04",
        publication_stage_name="sm_property_prefilter",
        legacy_steps=("step2_5",),
        working_root=working_root,
        strict_lineage=strict_lineage,
        manifest_name="phase_manifest.json",
        manifest_extra={
            "phase_id": "property_prefilter",
            "phase_scope": "pre_docking_molecule_gate",
            "separate_from": "post_pose_property_scoring",
        },
    )
    stage_manifests[ol.SM_PROPERTY_DIR]["phases"] = [{
        "phase_id": "property_prefilter",
        "manifest": f"{ol.SM_PROPERTY_DIR}/00_prefilter/phase_manifest.json",
    }]
    layout.publish_json(
        layout.path(ol.SM_PROPERTY_DIR, ol.STEP_MANIFEST),
        stage_manifests[ol.SM_PROPERTY_DIR],
        role="stage_manifest",
        derived_from=[
            str(stage_manifests[ol.SM_PROPERTY_DIR]["artifact_index"]),
            f"{ol.SM_PROPERTY_DIR}/00_prefilter/phase_manifest.json",
        ],
    )
    return stage_manifests


def _pipeline_source_snapshot(
    layout: ol.RunLayout,
    working_root: Path,
) -> dict[str, dict[str, Any]]:
    """Fingerprint local source files so verification detects stale views."""
    rows: list[dict[str, Any]] = []
    for name in _CONTROL_INPUT_FILES:
        path = working_root / name
        if path.is_file():
            rows.append(_artifact_row(
                layout=layout,
                pipeline_step_id="pipeline",
                declared_by="run_control",
                artifact_name=Path(name).stem,
                path=path,
                producer="pipeline_entry",
            ))
    for legacy_step in _PIPELINE_STEP_NAMES:
        _summary, step_rows = _runtime_step_evidence(
            layout,
            working_root,
            legacy_step,
            strict_lineage=True,
        )
        rows.extend(step_rows)

    # Seal every local file cited by the runtime evidence, including injected
    # inputs and artifacts reused from an upstream round.  A final view is stale
    # if any generation that fed it changes, regardless of directory ownership.
    snapshot: dict[str, dict[str, Any]] = {}
    for row in rows:
        raw_relative = str(row.get("path") or "")
        if not raw_relative:
            continue
        if row.get("exists") is not True:
            raise ValueError(
                f"local publication source is missing: {raw_relative}"
            )
        identity = {
            "kind": row.get("kind") or "path",
            "bytes": row.get("bytes"),
            "sha256": row.get("sha256"),
        }
        previous = snapshot.get(raw_relative)
        if previous is not None and previous != identity:
            raise ValueError(
                f"publication source changed while indexing: {raw_relative}"
            )
        snapshot[raw_relative] = identity

    required = {
        layout.relative(working_root / relative)
        for relative in _FINALIZATION_RESULT_FILES
    }
    missing = sorted(required - set(snapshot))
    if missing:
        raise ValueError(
            "pipeline finalization results are absent from source snapshot: "
            + ", ".join(missing)
        )
    return {key: snapshot[key] for key in sorted(snapshot)}


def publish_run(
    root: Path | str,
    *,
    records: Sequence[Mapping[str, Any]],
    metadata: Mapping[str, Any],
    rules: Sequence[Mapping[str, Any]] = POSE_PROPERTY_METRIC_RULES,
    identity: Mapping[str, Any] | None = None,
    profile: Mapping[str, Any] | None = None,
    step_sources: Mapping[str, Sequence[Any]] | None = None,
    shared_assets: Mapping[str, Any] | None = None,
    inputs: Mapping[str, Any] | None = None,
    failures: Iterable[Mapping[str, Any]] = (),
    session_id: str = "",
    working_root: Path | str | None = None,
    publication_mode: str = "standalone",
) -> dict[str, str]:
    """Write the whole tree and return ``role -> path`` for what was written."""
    if publication_mode not in _PUBLICATION_MODES:
        raise ValueError(
            f"invalid publication_mode {publication_mode!r}; expected "
            "'standalone' or 'pipeline_final'"
        )
    layout = ol.RunLayout(root=Path(root), session_id=session_id)
    resolved_working_root = (
        Path(working_root).expanduser().resolve()
        if working_root is not None else None
    )
    if publication_mode == "pipeline_final" and resolved_working_root is None:
        raise ValueError("pipeline_final publication requires working_root")
    if publication_mode == "pipeline_final":
        _require_completed_pipeline_sources(resolved_working_root)
    metadata = flatten_scoring_metadata(metadata)
    rows = _ordered_records(records)
    steps = _active_steps(metadata)
    by_category = _rules_by_category(rules)
    by_tool = _rules_by_tool(rules)
    step_sources = dict(step_sources or {})
    identity = {
        "session_id": layout.session_id,
        "generated_at": layout.generated_at,
        "scoring_schema_version": metadata.get("schema_version"),
        "metric_catalog_version": metadata.get("metric_rules_version"),
        **dict(identity or {}),
    }
    written: dict[str, str] = {}

    def note(role: str, path: Path) -> None:
        written[role] = str(path)

    workflow_manifests = _publish_workflow_evidence(
        layout,
        resolved_working_root,
        strict_lineage=publication_mode == "pipeline_final",
    )
    source_snapshot = (
        _pipeline_source_snapshot(layout, resolved_working_root)
        if publication_mode == "pipeline_final" else {}
    )

    # -- tool layer first: the tables above it cite these paths -------------
    tool_paths: dict[tuple[str, str], str] = {}
    tool_runs_by_step: dict[str, set[str]] = {step_id: set() for step_id in steps}
    for (step_id, tool_run_id), tool_rules in sorted(by_tool.items()):
        if step_id not in steps:
            continue
        logical_tool_ids = sorted({ol.rule_tool_id(rule) for rule in tool_rules})
        metric_ids = _measured_metric_ids(rows, _ordered_metric_ids(tool_rules))
        if not metric_ids and not _has_tool_execution_evidence(
            rows, logical_tool_ids,
        ):
            # Neither a measurement nor an attempted execution needs an audit.
            continue
        columns, tool_rows = _tool_table(rows, metric_ids, step_id)
        directory = layout.tool_run_dir(step_id, tool_run_id)
        path = layout.publish_tsv(
            directory / ol.TOOL_TABLE, columns, tool_rows, role="tool_table",
        )
        for tool_id in logical_tool_ids:
            tool_paths[(step_id, tool_id)] = layout.relative(path)
        tool_runs_by_step[step_id].add(tool_run_id)
        layout.publish_json(
            directory / ol.INTERNAL_DIR / ol.TOOL_SUMMARY,
            _tool_summary(
                rows, tool_rules, tool_run_id, logical_tool_ids, step_id,
            ),
            role="tool_summary",
        )
        layout.publish_json(
            directory / ol.INTERNAL_DIR / ol.PROVENANCE,
            {
                "tool_run_id": tool_run_id,
                "tool_id": (
                    logical_tool_ids[0] if len(logical_tool_ids) == 1
                    else tool_run_id
                ),
                "tool_ids": logical_tool_ids,
                "step_id": step_id,
                "generated_at": layout.generated_at,
                "read_from": [
                    layout.relative(source)
                    for source in step_sources.get(step_id, ())
                ],
                "raw_output_retained_by": "the step working directory",
            },
            role="provenance",
        )
        layout.write_internal_index(step_id, tool_run_id)

    # -- category layer -----------------------------------------------------
    category_paths: dict[str, list[str]] = {step_id: [] for step_id in steps}
    for code in ol.CATEGORY_ORDER:
        category_rules = by_category.get(code)
        step_id = CATEGORY_CONFIG[code]["step_id"]
        if not category_rules or step_id not in steps:
            continue
        metric_ids = _measured_metric_ids(
            rows, _ordered_metric_ids(category_rules),
        )
        columns, category_rows = _category_table(rows, code, metric_ids, step_id)
        measured = set(metric_ids)
        tool_columns: dict[str, list[str]] = {}
        for rule in category_rules:
            if str(rule["metric_id"]) not in measured:
                continue
            tool_columns.setdefault(ol.rule_tool_id(rule), []).append(
                str(rule["metric_id"]),
            )
        cited = sorted({
            tool_paths[(step_id, tool_id)]
            for tool_id in sorted(tool_columns)
            if (step_id, tool_id) in tool_paths
        })
        directory = layout.category_dir(code)
        table = layout.publish_tsv(
            directory / ol.CATEGORY_TABLE, columns, category_rows,
            role="category_table", derived_from=cited,
        )
        layout.publish_tsv(
            directory / ol.CATEGORY_REPORT,
            ol.CATEGORY_REPORT_COLUMNS, _category_report_rows(rows, code),
            role="category_report", derived_from=[layout.relative(table)],
        )
        layout.write_source_map(
            code,
            metric_columns={
                tool_id: sorted(columns_for)
                for tool_id, columns_for in tool_columns.items()
            },
            artifacts={
                tool_id: [
                    source for source in step_sources.get(step_id, ())
                ]
                for tool_id in tool_columns
            },
        )
        category_paths[step_id].append(layout.relative(table))

    # -- step layer ---------------------------------------------------------
    step_paths: list[str] = []
    for step_id in steps:
        step_rules = [
            rule for rule in rules if str(rule.get("step_id")) == step_id
        ]
        metric_ids = _measured_metric_ids(rows, _ordered_metric_ids(step_rules))
        columns, step_rows = _step_table(rows, step_id, metric_ids)
        directory = layout.step_dir(step_id)
        table = layout.publish_tsv(
            directory / ol.STEP_TABLE, columns, step_rows,
            role="step_table", derived_from=category_paths.get(step_id, []),
        )
        layout.publish_tsv(
            directory / ol.STEP_REPORT,
            ol.STEP_REPORT_COLUMNS, _step_report_rows(rows, step_id),
            role="step_report", derived_from=[layout.relative(table)],
        )
        workflow_manifest = workflow_manifests.get(
            ol.STEP_DIR[step_id],
            {},
        )
        layout.publish_json(
            directory / ol.STEP_MANIFEST,
            {
                **workflow_manifest,
                "step_id": step_id,
                "scoring_step_id": step_id,
                "step_dir": ol.STEP_DIR[step_id],
                "generated_at": layout.generated_at,
                "row_id_field": ol.ROW_ID_FIELD[step_id],
                "rows": len(rows),
                "categories": [
                    ol.CATEGORY_DIR[code]
                    for code in ol.STEP_CATEGORY_ORDER[step_id]
                    if code in by_category
                ],
                "tool_runs": sorted(tool_runs_by_step.get(step_id, set())),
                "reference_parent_weight": (
                    metadata.get("reference_step_parent_weight") or {}
                ).get(step_id),
                "read_from": [
                    layout.relative(source)
                    for source in step_sources.get(step_id, ())
                ],
            },
            role="step_manifest",
        )
        step_paths.append(layout.relative(table))
        note(f"{step_id}_table", table)

    # -- the two step sidecars the output contract singles out -------------
    if "step6" in steps:
        layout.publish_json(
            layout.step_dir("step6") / ol.INPUT_REFS,
            {
                "step_id": "step6",
                "generated_at": layout.generated_at,
                "upstream_step": "step5" if "step5" in steps else None,
                "upstream_table": (
                    layout.relative(layout.step_dir("step5") / ol.STEP_TABLE)
                    if "step5" in steps else None
                ),
                "rows": [
                    {
                        "row_id": ol.row_id(record, "step6"),
                        "pose_id": record.get("pose_id"),
                        "selected_pose_index": record.get("selected_pose_index"),
                        "source_pose_index": record.get("source_pose_index"),
                        "pose_file": (
                            layout.relative(_lookup("pose_file", record, identity))
                            if _lookup("pose_file", record, identity) else None
                        ),
                    }
                    for record in rows
                ],
            },
            role="input_refs",
        )
    if shared_assets and "step5" in steps:
        layout.publish_json(
            layout.shared_assets_dir("step5") / ol.ASSET_MANIFEST,
            {
                "generated_at": layout.generated_at,
                "note": (
                    "References, not copies: the prepared receptor, pocket and "
                    "ligand files stay where the preparing step wrote them."
                ),
                **{
                    key: [layout.relative(item) for item in value]
                    if isinstance(value, (list, tuple))
                    else layout.relative(value)
                    if isinstance(value, (str, Path)) and value else value
                    for key, value in shared_assets.items()
                },
            },
            role="asset_manifest",
        )

    # -- final report layer -------------------------------------------------
    final_dir = layout.final_dir()
    metric_ids = _ordered_metric_ids([
        rule for rule in rules if str(rule.get("step_id")) in set(steps)
    ])
    note("all_metrics_table", layout.publish_tsv(
        final_dir / ol.ALL_METRICS_TABLE,
        ALL_METRICS_STEM + tuple(
            f"{metric_id}{suffix}"
            for metric_id in metric_ids
            for suffix in ("", "__label", "__gate")
        ),
        [_all_metrics_row(record, metric_ids) for record in rows],
        role="all_metrics_table", derived_from=step_paths,
    ))
    unimplemented = frozenset(
        metadata.get("unimplemented_reference_categories") or ()
    )
    note("final_report_full", layout.publish_tsv(
        final_dir / ol.FINAL_REPORT_FULL,
        ol.FINAL_REPORT_FULL_COLUMNS,
        [_full_row(record, identity, unimplemented) for record in rows],
        role="final_report_full", derived_from=step_paths,
    ))
    note("final_report", layout.publish_tsv(
        final_dir / ol.FINAL_REPORT,
        ol.FINAL_REPORT_COLUMNS,
        [_final_row(record, identity) for record in rows],
        role="final_report",
        derived_from=[layout.relative(final_dir / ol.FINAL_REPORT_FULL)],
    ))
    resolved_profile = {**_default_profile(metadata), **dict(profile or {})}
    note("metric_catalog", layout.publish_text(
        final_dir / ol.METRIC_CATALOG,
        _catalog_markdown(rules, metadata, resolved_profile),
        role="metric_catalog",
    ))

    for failure in failures:
        layout.register_failure(
            stage=str(failure.get("stage") or "unknown"),
            tool_id=failure.get("tool_id"),
            reason=str(failure.get("reason") or ""),
            detail=failure.get("detail"),
            row_ids=failure.get("row_ids") or (),
        )
    registry = final_dir / ol.FAILURE_REGISTRY
    if registry.exists():
        layout.record(registry, role="failure_registry")
        note("failure_registry", registry)

    manifest, lineage = layout.finalize(
        profile=resolved_profile,
        inputs=dict(inputs or {}),
        extra={
            "row_count": len(rows),
            "active_steps": list(steps),
            "publication_mode": publication_mode,
            "source_root": (
                layout.relative(resolved_working_root)
                if resolved_working_root is not None else None
            ),
            "source_snapshot": source_snapshot,
            "finalization": {
                "status": (
                    "complete"
                    if publication_mode == "pipeline_final" else "not_applicable"
                ),
                "required_result_files": list(_FINALIZATION_RESULT_FILES),
            },
            "publication_stages": list(ol.PUBLICATION_STAGE_DIRS),
            "runtime_stage_map": {
                directory: [
                    _PIPELINE_STEP_NAMES[legacy_step]
                    for legacy_step in legacy_steps
                ]
                for directory, legacy_steps in _WORKFLOW_STAGE_SOURCES.items()
            } | {
                ol.INPUTS_DIR: ["pipeline_config", "resolved_run_plan"],
                ol.SM_PROPERTY_DIR: [
                    "property_prefilter",
                    "property_filter",
                ],
            },
        },
    )
    note("run_manifest", manifest)
    note("lineage_graph", lineage)
    return written


def _owned_swap_directories(parent: Path, pattern: str) -> list[Path]:
    """Return publisher-owned swap directories, rejecting link tricks."""
    directories = []
    for path in sorted(parent.glob(pattern)):
        if path.is_symlink() or not path.is_dir():
            raise ValueError(f"unexpected publication swap entry: {path}")
        directories.append(path)
    return directories


def publish_run_transactionally(
    root: Path | str,
    *,
    renderer: Callable[..., dict[str, str]] = publish_run,
    **publish_kwargs: Any,
) -> dict[str, str]:
    """Build a fresh publication, verify it, then recoverably swap it in.

    A resume must never write incrementally over an older view: a tool removed
    from the new plan, or a failure that no longer exists, would otherwise stay
    on disk while disappearing from the manifest.  Staging is a sibling of the
    final directory so publication-relative ``../stepN`` references survive the
    rename unchanged.  POSIX cannot atomically replace a non-empty directory;
    the short two-rename swap is serialized on the parent directory, rolls back
    Python-level failures, and restores an orphan backup on the next invocation
    after a process-level crash.
    """
    requested_target = Path(root).expanduser()
    if requested_target.is_symlink():
        raise ValueError(
            f"publication target must not be a symlink: {requested_target}",
        )
    target = requested_target.resolve()
    session_id = str(publish_kwargs.get("session_id") or "").strip()
    if not session_id:
        raise ValueError("transactional publication requires an explicit session_id")
    target.parent.mkdir(parents=True, exist_ok=True)
    directory_fd = os.open(
        target.parent,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0),
    )
    staging: Path | None = None
    try:
        fcntl.flock(directory_fd, fcntl.LOCK_EX)
        if requested_target.is_symlink():
            raise ValueError(
                f"publication target must not be a symlink: {requested_target}",
            )
        backups = _owned_swap_directories(
            target.parent, f".{target.name}.backup-*",
        )
        orphan_staging = _owned_swap_directories(
            target.parent, f".{target.name}.staging-*",
        )
        if not target.exists() and backups:
            if len(backups) != 1:
                raise RuntimeError(
                    f"cannot choose among orphan publication backups: {backups}",
                )
            recovered = backups.pop()
            recovered.rename(target)
        for orphan in orphan_staging:
            shutil.rmtree(orphan)
        if target.exists() and not target.is_dir():
            raise ValueError(f"publication target is not a directory: {target}")

        staging = Path(tempfile.mkdtemp(
            prefix=f".{target.name}.staging-",
            dir=str(target.parent),
        ))
        staged_written = renderer(
            staging,
            **publish_kwargs,
        )
        ol.verify_publication(
            staging,
            expected_session_id=session_id,
            require_pipeline_final=(
                publish_kwargs.get("publication_mode") == "pipeline_final"
            ),
        )
        backup: Path | None = None
        try:
            if target.exists():
                backup = target.with_name(
                    f".{target.name}.backup-{uuid.uuid4().hex}",
                )
                target.rename(backup)
            staging.rename(target)
        except BaseException:
            if backup is not None and backup.exists() and not target.exists():
                backup.rename(target)
            raise
        written = {
            role: str(target / Path(path).resolve().relative_to(staging))
            for role, path in staged_written.items()
        }
        for old_view in [*backups, backup]:
            if old_view is None or not old_view.exists():
                continue
            # The prior publication is a derived view; the immutable step tree
            # remains its recovery source.  Failure to clean a backup must not
            # invalidate the newly verified view.
            try:
                shutil.rmtree(old_view)
            except OSError:
                pass
        return written
    finally:
        try:
            if staging is not None and staging.exists():
                shutil.rmtree(staging)
        finally:
            try:
                fcntl.flock(directory_fd, fcntl.LOCK_UN)
            finally:
                os.close(directory_fd)

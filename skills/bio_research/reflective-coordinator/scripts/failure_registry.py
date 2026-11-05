"""Harvest one round's "why is this empty" ledger for ``failure_registry.jsonl``.

The output-layout contract puts a registry in ``07_final_report/`` of the
published tree.  The tables give it its meaning: an
empty cell there is honest but mute -- it says a value is absent without saying
whether the tool crashed, was switched off, or was never part of this build.
Every entry here answers that question for one gap, and ``detail.kind`` says
whether the reason was a fault or a choice.

Three sources, because each is blind where the others see:

* The producing skills' own ``steps_skipped`` / ``operational_errors``.  A tool
  that failed for the whole cohort leaves no per-row trace at all: with no row
  holding a value, the scorer cannot know the metric was ever planned, so every
  assessment reads ``not_planned`` -- indistinguishable from a metric this build
  never intended to compute.  Only the producer knows, and only here does it say.
* ``planned_missing`` assessments.  The complement: a tool that answered for some
  rows and not others, which the producer usually does not notice.
* The coordinator's own blocks -- quality gates that could not run or did not
  pass, and per-row ranking blocks.

Nothing here is inferred from an absence alone.  Every entry names a declaration
or a computed status, and the metric columns it claims are blank are checked
against the published rows rather than assumed from the reason text.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


_SHARED = Path(__file__).resolve().parents[2] / "shared"
if str(_SHARED) not in sys.path:
    sys.path.insert(0, str(_SHARED))

import output_layout as ol  # noqa: E402
from metric_scoring import POSE_PROPERTY_METRIC_RULES, STEP_IDS  # noqa: E402


#: Listing every affected row would turn one failure into a wall of ids.
ROW_ID_LIMIT = 25

#: Upstream skips name a sub-step; the catalog names the tool that owns the
#: metric.  Spelled out rather than fuzzy-matched: ``mmgbsa`` and ``amber_mmpbsa``
#: share no token, and a near-miss would file the failure against a tool whose
#: columns are fine.
DECLARED_ALIASES: dict[str, str] = {
    "plip_interactions": "plip",
    "prolif_interactions": "prolif",
    "mmgbsa": "amber_mmpbsa",
    "structural_alerts": "pains",
    "posecheck_native_calibration": "posecheck",
    "rescore_rtmscore": "rtmscore",
    "openmm_md": "short_md",
}

#: Wording a producer uses when a step was switched off rather than broken.
DELIBERATE_MARKERS: tuple[str, ...] = (
    "disabled",
    "not requested",
    "not configured",
    "not enabled",
    "no candidates",
    "nothing to",
    "skipped by",
)

#: Where the coordinator's own findings sit, after the three step stages.
_STAGE_ORDER: dict[str, int] = {step_id: index for index, step_id in enumerate(STEP_IDS)}
_STAGE_ORDER["protocol"] = len(STEP_IDS)
_STAGE_ORDER["ranking"] = len(STEP_IDS) + 1


# ---------------------------------------------------------------------------
# Routing a declared name onto a step and a tool directory.
# ---------------------------------------------------------------------------


def _tool_index(rules: Sequence[Mapping[str, Any]]) -> dict[str, tuple[str, str]]:
    """Alias -> ``(step_id, tool_id)`` for every tool that owns a metric column."""
    index: dict[str, tuple[str, str]] = {}
    for rule in rules:
        step_id = str(rule.get("step_id") or "")
        if step_id not in _STAGE_ORDER:
            continue
        tool_id = ol.rule_tool_id(rule)
        aliases = (tool_id, ol.slug_tool_id(tool_id), ol.slug_tool_id(rule.get("tool")))
        for alias in aliases:
            if alias:
                index.setdefault(alias, (step_id, tool_id))
    return index


def _metrics_by_tool(
    rules: Sequence[Mapping[str, Any]],
) -> dict[tuple[str, str], list[str]]:
    """``(step_id, tool_id)`` -> the metric columns that tool is the sole source of."""
    owned: dict[tuple[str, str], list[str]] = {}
    for rule in rules:
        key = (str(rule.get("step_id") or ""), ol.rule_tool_id(rule))
        owned.setdefault(key, []).append(str(rule.get("metric_id")))
    return owned


def _route(
    name: Any, default_step: str, index: Mapping[str, tuple[str, str]]
) -> tuple[str, str, bool]:
    """``(step_id, tool_id, owns_published_columns)`` for a name a producer declared."""
    slug = ol.slug_tool_id(name)
    target = index.get(DECLARED_ALIASES.get(slug, slug)) or index.get(slug)
    if target:
        return target[0], target[1], True
    return default_step, slug or "unknown", False


def _deliberate(reason: Any) -> bool:
    """Whether a producer's wording describes a choice instead of a fault."""
    text = str(reason or "").lower()
    return any(marker in text for marker in DELIBERATE_MARKERS)

# ---------------------------------------------------------------------------
# One entry.
# ---------------------------------------------------------------------------


def _entry(
    stage: str,
    tool_id: str | None,
    reason: str,
    kind: str,
    *,
    row_ids: Sequence[str] = (),
    row_id_limit: int = ROW_ID_LIMIT,
    **detail: Any,
) -> dict[str, Any]:
    """Shape one ledger line the way ``RunLayout.register_failure`` expects it."""
    rows = [str(item) for item in row_ids if item]
    payload: dict[str, Any] = {"kind": kind}
    payload.update({key: value for key, value in detail.items() if value != []})
    if rows:
        payload["affected_rows"] = len(rows)
    if len(rows) > row_id_limit:
        payload["row_ids_truncated"] = True
        rows = rows[:row_id_limit]
    return {
        "stage": stage,
        "tool_id": tool_id,
        "reason": str(reason),
        "detail": payload,
        "row_ids": rows,
    }


def _measured_metric_ids(records: Sequence[Mapping[str, Any]]) -> set[str]:
    """Metric ids at least one published row holds a value for."""
    return {
        str(assessment.get("metric_id"))
        for record in records
        for assessment in record.get("metric_assessments") or ()
        if assessment.get("availability_status") == "available"
    }

# ---------------------------------------------------------------------------
# Source 1: what the producing skill said about itself.
# ---------------------------------------------------------------------------


def _columns(
    key: tuple[str, str],
    owned: Mapping[tuple[str, str], Sequence[str]],
    measured: set[str],
    routed: bool,
) -> dict[str, Any]:
    """The metric columns a named tool owns, split by whether a value arrived.

    Checked against the published rows: a producer that reports a skip but whose
    numbers are all present is worth seeing, and so is the opposite.
    """
    if not routed:
        return {"catalog_tool": False}
    metric_ids = sorted(owned.get(key, ()))
    return {
        "blank_metric_ids": [item for item in metric_ids if item not in measured],
        "populated_metric_ids": [item for item in metric_ids if item in measured],
    }


def _declared(
    report: Mapping[str, Any],
    *,
    default_step: str,
    source: str,
    index: Mapping[str, tuple[str, str]],
    owned: Mapping[tuple[str, str], Sequence[str]],
    measured: set[str],
    row_id_limit: int,
) -> list[dict[str, Any]]:
    """Harvest ``steps_skipped`` / ``operational_errors`` / ``error`` from one report."""
    entries: list[dict[str, Any]] = []
    for item in report.get("steps_skipped") or ():
        name = item.get("step") if isinstance(item, Mapping) else item
        reason = item.get("reason") if isinstance(item, Mapping) else ""
        step_id, tool_id, routed = _route(name, default_step, index)
        entries.append(_entry(
            step_id,
            tool_id,
            reason or f"{name} was skipped without a stated reason",
            "step_disabled" if _deliberate(reason) else "step_skipped",
            row_id_limit=row_id_limit,
            source=source,
            declared_as=str(name),
            **_columns((step_id, tool_id), owned, measured, routed),
        ))
    for item in report.get("operational_errors") or ():
        name = item.get("tool") if isinstance(item, Mapping) else item
        reason = item.get("reason") if isinstance(item, Mapping) else ""
        row_ids = (
            item.get("row_ids")
            or item.get("affected_molecules")
            or ()
        ) if isinstance(item, Mapping) else ()
        step_id, tool_id, routed = _route(name, default_step, index)
        entries.append(_entry(
            step_id,
            tool_id,
            reason or f"{name} reported an error without a stated reason",
            "tool_error",
            row_ids=row_ids,
            row_id_limit=row_id_limit,
            source=source,
            declared_as=str(name),
            **_columns((step_id, tool_id), owned, measured, routed),
        ))
    summary_error = str(report.get("error") or "").strip()
    if summary_error and not (report.get("operational_errors") or ()):
        entries.append(_entry(
            default_step,
            None,
            summary_error,
            "stage_error",
            row_id_limit=row_id_limit,
            source=source,
        ))
    return entries

# ---------------------------------------------------------------------------
# Source 2: a report the run counted on but never got.
# ---------------------------------------------------------------------------


def _step_report(
    report: Mapping[str, Any] | None,
    path: Any,
    *,
    step_id: str,
    payload_key: str,
    source: str,
    row_id_limit: int,
) -> list[dict[str, Any]]:
    """Flag a declared input that was absent or carried no rows.

    A property report that was never declared is not a failure -- the step simply
    is not active.  One that *was* declared and then could not be opened is,
    because coverage still holds the step accountable while the caller sees no
    error: every metric it owns silently reads ``not_planned``.
    """
    if not path:
        return []
    if not Path(str(path)).exists():
        return [_entry(
            step_id,
            None,
            f"{source} was declared at {path} but no file is there",
            "missing_step_report",
            row_id_limit=row_id_limit,
            source=source,
            declared_path=str(path),
        )]
    if not (report or {}).get(payload_key):
        return [_entry(
            step_id,
            None,
            f"{source} carries no {payload_key} rows",
            "empty_step_report",
            row_id_limit=row_id_limit,
            source=source,
            declared_path=str(path),
        )]
    return []

# ---------------------------------------------------------------------------
# Source 3: a tool that answered for some rows and not others.
# ---------------------------------------------------------------------------


def _partial_coverage(
    records: Sequence[Mapping[str, Any]],
    rules: Sequence[Mapping[str, Any]],
    *,
    row_id_limit: int,
) -> list[dict[str, Any]]:
    """Group ``planned_missing`` assessments per tool.

    ``planned_missing`` is the scorer's own verdict: the metric was planned for
    this cohort and this row has no value for it.  The producer rarely notices a
    partial miss, so this is the only place it surfaces.
    """
    by_metric = {str(rule.get("metric_id")): rule for rule in rules}
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    for record in records:
        for assessment in record.get("metric_assessments") or ():
            if assessment.get("availability_status") != "planned_missing":
                continue
            metric_id = str(assessment.get("metric_id"))
            rule = by_metric.get(metric_id)
            step_id = str((rule or assessment).get("step_id") or "")
            tool_id = (
                ol.rule_tool_id(rule) if rule
                else ol.slug_tool_id(assessment.get("tool")) or "unknown"
            )
            bucket = grouped.setdefault(
                (step_id, tool_id), {"metrics": {}, "rows": {}},
            )
            bucket["metrics"][metric_id] = bucket["metrics"].get(metric_id, 0) + 1
            bucket["rows"][ol.row_id(record, step_id)] = None
    entries: list[dict[str, Any]] = []
    for (step_id, tool_id), bucket in sorted(grouped.items()):
        metrics = bucket["metrics"]
        rows = [item for item in bucket["rows"] if item]
        entries.append(_entry(
            step_id,
            tool_id,
            f"{len(metrics)} planned metric(s) produced no value for "
            f"{len(rows)} of {len(records)} rows",
            "partial_metric_coverage",
            row_ids=rows,
            row_id_limit=row_id_limit,
            missing_rows_by_metric=dict(sorted(metrics.items())),
            rows_total=len(records),
        ))
    return entries

# ---------------------------------------------------------------------------
# Sources 4 and 5: the coordinator's own blocks.
# ---------------------------------------------------------------------------


def _blocked_rows(records: Sequence[Mapping[str, Any]]) -> dict[str, list[str]]:
    """Blocker string -> the published row ids it kept out of the ranking."""
    grouped: dict[str, list[str]] = {}
    for record in records:
        for blocker in record.get("ranking_blocked_by") or ():
            grouped.setdefault(str(blocker), []).append(ol.row_id(record))
    return grouped


def _quality_gates(
    report: Mapping[str, Any],
    blocked: Mapping[str, Sequence[str]],
    *,
    row_id_limit: int,
) -> list[dict[str, Any]]:
    """One entry per gate that did not pass or could not be evaluated.

    ``round_summary`` already holds the aggregated verdicts, so the reasons a
    reader sees here are the ones the decision was actually made on.
    """
    gates = report.get("quality_gates") or {}
    summary = report.get("round_summary") or {}
    entries: list[dict[str, Any]] = []
    for key, kind, fallback in (
        ("failed_quality_gates", "quality_gate_failed", "did not pass"),
        (
            "not_applicable_quality_gates",
            "quality_gate_not_applicable",
            "could not be evaluated: no valid reference for it exists",
        ),
    ):
        for name in summary.get(key) or ():
            gate = gates.get(name) or {}
            entries.append(_entry(
                "protocol",
                None,
                str(gate.get("reason") or f"quality gate {name} {fallback}"),
                kind,
                row_id_limit=row_id_limit,
                gate=str(name),
                gate_status=str(gate.get("gate") or ""),
                blocked_rows=len(blocked.get(f"quality_gate:{name}", ())),
            ))
    return entries


def _row_blocks(
    blocked: Mapping[str, Sequence[str]], *, row_id_limit: int
) -> list[dict[str, Any]]:
    """Per-row ranking blocks, grouped by blocker.

    Gate-derived blockers are skipped: ``_quality_gates`` already carries them,
    with the reason text this one does not have.
    """
    return [
        _entry(
            "ranking",
            None,
            f"ranking blocked by {blocker}",
            "row_blocked",
            row_ids=rows,
            row_id_limit=row_id_limit,
            blocker=blocker,
        )
        for blocker, rows in sorted(blocked.items())
        if not blocker.startswith("quality_gate:")
    ]

# ---------------------------------------------------------------------------
# The harvest.
# ---------------------------------------------------------------------------


def _sort_key(entry: Mapping[str, Any]) -> tuple[Any, ...]:
    stage = str(entry.get("stage") or "")
    return (
        _STAGE_ORDER.get(stage, len(_STAGE_ORDER)),
        stage,
        str((entry.get("detail") or {}).get("kind") or ""),
        str(entry.get("tool_id") or ""),
        str(entry.get("reason") or ""),
    )


def collect_failures(
    report: Mapping[str, Any],
    *,
    eval_report: Mapping[str, Any] | None = None,
    prop_report: Mapping[str, Any] | None = None,
    docking_report: Mapping[str, Any] | None = None,
    eval_file: Any = None,
    prop_file: Any = None,
    docking_file: Any = None,
    rules: Sequence[Mapping[str, Any]] | None = None,
    row_id_limit: int = ROW_ID_LIMIT,
) -> list[dict[str, Any]]:
    """Every reason this round has a gap, ready for ``publish_run(failures=...)``.

    Returns ``[]`` for a clean round, which is what keeps ``failure_registry.jsonl``
    off disk entirely: the file's presence is itself the signal.
    """
    rules = list(rules or POSE_PROPERTY_METRIC_RULES)
    records = list(report.get("molecules") or ())
    index = _tool_index(rules)
    owned = _metrics_by_tool(rules)
    measured = _measured_metric_ids(records)
    blocked = _blocked_rows(records)

    entries: list[dict[str, Any]] = []
    entries += _step_report(
        eval_report, eval_file, step_id="step6", payload_key="evaluation",
        source="evaluation_report", row_id_limit=row_id_limit,
    )
    entries += _step_report(
        prop_report, prop_file, step_id="step4", payload_key="results",
        source="property_filter_report", row_id_limit=row_id_limit,
    )
    for source_report, default_step, source in (
        (docking_report, "step5", "docking_report"),
        (eval_report, "step6", "evaluation_report"),
        (prop_report, "step4", "property_filter_report"),
    ):
        entries += _declared(
            source_report or {},
            default_step=default_step,
            source=source,
            index=index,
            owned=owned,
            measured=measured,
            row_id_limit=row_id_limit,
        )
    entries += _partial_coverage(records, rules, row_id_limit=row_id_limit)
    entries += _quality_gates(report, blocked, row_id_limit=row_id_limit)
    entries += _row_blocks(blocked, row_id_limit=row_id_limit)
    entries.sort(key=_sort_key)
    return entries

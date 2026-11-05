"""HTML template — pure stdlib. Jinja2 optional.

Reads the fully-aggregated `report_data` dict + asset manifest + plot status
and produces a single HTML string. Every renderer already handled its own
graceful-skip; this file just decides what sections to show and stitches
strings together.
"""
from __future__ import annotations

import html as _html
import json
from pathlib import Path
from typing import Any


def check_jinja_available() -> bool:
    try:
        import jinja2  # noqa: F401
        return True
    except ImportError:
        return False


def render_html(
    report_data: dict[str, Any],
    asset_manifest: dict[str, Any],
    plots_status: dict[str, Any],
    pymol_cover_rel: str | None,
    main_viewer_html: str | None,
    html_output_path: Path,
    project_root: Path,
) -> str:
    """Produce a self-contained HTML string.

    We use plain string building — Jinja2 is intentionally not used to keep
    the skill's hard-dependency surface at zero. `check_jinja_available` is
    exposed for callers that want to log the choice.
    """
    header = report_data.get("header") or {}
    candidates = report_data.get("candidates") or []
    gates = report_data.get("quality_gates") or []
    screening = report_data.get("screening") or {}
    trend = report_data.get("trend") or []
    missing = report_data.get("missing_sources") or []
    figures = (plots_status or {}).get("figures") or {}
    plot_errors = (plots_status or {}).get("errors") or {}

    parts: list[str] = []
    parts.append(_head(header, project_root, html_output_path, asset_manifest))
    parts.append('<body>')
    parts.append('<div class="report-container">')

    parts.append(_render_header(header, gates, report_data.get("mode"), pymol_cover_rel))
    parts.append(_render_trend(trend) if report_data.get("mode") == "final" and trend else "")
    parts.append(_render_main_viewer(main_viewer_html))
    parts.append(_render_synthesis_audit(report_data.get("synthesis_summary") or {}))
    parts.append(_render_candidates(candidates))
    parts.append(_render_candidate_reviews(report_data.get("candidate_reviews") or {}))
    parts.append(_render_qc_panels(figures, plot_errors))
    parts.append(_render_screening(screening, figures))
    parts.append(_render_footer(report_data, asset_manifest, plots_status, missing))

    parts.append('</div>')
    parts.append(_body_scripts(project_root, html_output_path, asset_manifest))
    parts.append('</body></html>')
    return "\n".join(p for p in parts if p)


# ---------- <head> + shared CSS ----------

def _head(header: dict, project_root: Path, html_path: Path, asset_manifest: dict) -> str:
    title = _e(header.get("title") or "SBDD Report")
    molstar_css = _relpath(project_root / "assets" / "molstar.css", html_path.parent)
    has_molstar_css = ((asset_manifest.get("assets") or {}).get("molstar.css") or {}).get("present")
    css_link = f'<link rel="stylesheet" href="{_e(molstar_css)}">' if has_molstar_css else ""
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
{css_link}
<style>
{_CSS}
</style>
</head>"""


_CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 0; background: #f6f7f9; color: #24292f; }
.report-container { max-width: 1280px; margin: 0 auto; padding: 24px; }
.section { background: #fff; border: 1px solid #e5e7eb; border-radius: 8px; padding: 20px 24px; margin-bottom: 20px; box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
.section h2 { margin: 0 0 12px 0; font-size: 18px; border-bottom: 1px solid #eef0f3; padding-bottom: 8px; }
.section h3 { margin: 12px 0 8px 0; font-size: 14px; color: #6b7280; text-transform: uppercase; letter-spacing: 0.03em; }
.header-row { display: flex; gap: 16px; align-items: baseline; flex-wrap: wrap; }
.header-title { font-size: 22px; font-weight: 600; }
.header-sub { font-size: 13px; color: #6b7280; }
.chip { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 12px; margin-right: 6px; margin-bottom: 4px; background: #eef0f3; color: #24292f; }
.chip-pass { background: #dcfce7; color: #166534; }
.chip-warn { background: #fef3c7; color: #92400e; }
.chip-fail { background: #fee2e2; color: #991b1b; }
.chip-unknown { background: #e5e7eb; color: #6b7280; }
.gate-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 10px; margin-top: 10px; }
.gate { padding: 10px 12px; border-radius: 6px; border: 1px solid #e5e7eb; background: #fafbfc; }
.gate-name { font-size: 12px; color: #6b7280; margin-bottom: 4px; }
.gate-value { font-size: 15px; font-weight: 600; }
.gate-source { font-size: 11px; color: #9ca3af; margin-top: 4px; font-family: SFMono-Regular, Menlo, monospace; }
.card { border: 1px solid #e5e7eb; border-radius: 8px; padding: 14px; margin-bottom: 14px; background: #fff; }
.card-header { display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 10px; }
.card-title { font-size: 16px; font-weight: 600; }
.card-grid { display: grid; grid-template-columns: 380px 300px 1fr; gap: 16px; align-items: start; }
@media (max-width: 900px) { .card-grid { grid-template-columns: 1fr; } }
.viewer-3dmol { border: 1px solid #e5e7eb; border-radius: 6px; background: #fff; }
.viewer-placeholder { background: #fafbfc; }
.viewer-error { padding: 8px; font-size: 12px; color: #b45309; text-align: center; }
.molstar-main-viewer { border: 1px solid #e5e7eb; border-radius: 6px; overflow: hidden; }
.structure-2d { text-align: center; }
.structure-2d img { max-width: 100%; height: auto; border: 1px solid #eef0f3; border-radius: 4px; background: #fff; }
.structure-2d .smiles { font-family: SFMono-Regular, Menlo, monospace; font-size: 11px; color: #6b7280; margin-top: 6px; word-break: break-all; }
.metric-table { font-size: 13px; width: 100%; border-collapse: collapse; }
.metric-table td { padding: 3px 6px; }
.metric-table td.k { color: #6b7280; width: 55%; }
.metric-table td.v { text-align: right; font-family: SFMono-Regular, Menlo, monospace; font-weight: 500; }
.diagnosis-list { margin: 8px 0 0 0; padding-left: 18px; font-size: 12px; color: #4b5563; }
.diagnosis-list li { margin-bottom: 2px; }
.plots-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; }
.plot img { width: 100%; height: auto; border: 1px solid #e5e7eb; border-radius: 4px; background: #fff; }
.plot-caption { font-size: 12px; color: #6b7280; margin-top: 4px; text-align: center; }
.trend-table { width: 100%; border-collapse: collapse; font-size: 13px; }
.trend-table th, .trend-table td { text-align: left; padding: 6px 8px; border-bottom: 1px solid #eef0f3; }
.trend-table th { background: #f9fafb; color: #6b7280; font-weight: 500; }
.footer { font-size: 12px; color: #9ca3af; padding: 12px 0; }
.missing-list { font-family: SFMono-Regular, Menlo, monospace; font-size: 11px; color: #b45309; background: #fffbeb; border-left: 3px solid #f59e0b; padding: 8px 12px; margin-top: 8px; }
.pymol-cover { display: block; max-width: 100%; margin: 12px 0; border-radius: 6px; border: 1px solid #e5e7eb; }
"""


# ---------- sections ----------

def _render_header(header: dict, gates: list[dict], mode: str | None, pymol_cover_rel: str | None) -> str:
    decision = header.get("decision") or "—"
    outcome = header.get("decision_outcome") or "unknown"
    goal = header.get("goal_met")
    nomination = header.get("valid_for_nomination")
    outcome_class = _chip_class_for_outcome(header)
    goal_label = "met" if goal is True else "not met" if goal is False else "unknown"
    nomination_label = "eligible" if nomination is True else "not eligible" if nomination is False else "unknown"

    chips = [
        f'<span class="chip">{_e(header.get("target_name") or "target?")}</span>',
        f'<span class="chip">PDB: {_e(header.get("pdb_id") or "?")}</span>',
        f'<span class="chip">{_e(header.get("round_label") or "round?")}</span>',
    ]
    if header.get("assessment_schema_version") == "round-assessment/v1":
        assessment = header.get("round_assessment") or {}
        chips.extend([
            f'<span class="chip {outcome_class}">assessment: {_e(assessment.get("status") or "unknown")}</span>',
            '<span class="chip">campaign: Agent review</span>',
        ])
    else:
        chips.extend([
            f'<span class="chip">decision: {_e(decision)}</span>',
            f'<span class="chip {outcome_class}">outcome: {_e(outcome)}</span>',
            f'<span class="chip">goal: {goal_label}</span>',
        ])
    chips.append(f'<span class="chip">nomination: {nomination_label}</span>')
    if header.get("primary_bottleneck"):
        chips.append(f'<span class="chip chip-warn">bottleneck: {_e(header["primary_bottleneck"])}</span>')

    counts = []
    for key, label in (("num_screened", "screened"), ("candidate_pool_size", "candidates"), ("num_keep", "keep"),
                      ("num_modify", "modify"), ("num_reject", "reject")):
        val = header.get(key)
        if val is not None:
            counts.append(f'<span class="chip">{label}: {_e(str(val))}</span>')
    intrinsic_keep = header.get("candidate_intrinsic_keep")
    blocked_keep = header.get("candidate_quality_gate_blocked_keep")
    if intrinsic_keep is not None:
        counts.append(f'<span class="chip chip-pass">intrinsic keep: {_e(str(intrinsic_keep))}</span>')
    if blocked_keep:
        counts.append(f'<span class="chip chip-warn">blocked keep: {_e(str(blocked_keep))}</span>')

    gate_cards = "".join(_render_gate_card(g) for g in gates)

    pymol_block = ""
    if pymol_cover_rel:
        pymol_block = f'<img class="pymol-cover" src="{_e(pymol_cover_rel)}" alt="PyMOL cover">'

    mode_label = "final report" if mode == "final" else "per-round report"
    rationale = header.get("decision_rationale")
    rationale_html = f'<p class="header-sub">{_e(rationale)}</p>' if rationale else ""

    return f"""<div class="section">
<div class="header-row">
  <div class="header-title">{_e(header.get("title") or "SBDD Report")}</div>
  <div class="header-sub">{_e(mode_label)}</div>
</div>
<div style="margin-top: 10px;">{"".join(chips)}</div>
{rationale_html}
<div style="margin-top: 6px;">{"".join(counts)}</div>
<h3>quality gates</h3>
<div class="gate-grid">{gate_cards}</div>
{pymol_block}
</div>"""


def _render_gate_card(gate: dict) -> str:
    passed = gate.get("passed")
    not_applicable = gate.get("not_applicable") is True
    css = "chip-unknown" if not_applicable else _chip_class_for_bool(passed)
    value = gate.get("value")
    threshold = gate.get("threshold")
    unit = gate.get("unit")
    value_str = _fmt_value(value, unit)
    threshold_str = _fmt_value(threshold, unit) if threshold is not None else "—"
    word = "n/a" if not_applicable else _pass_word(passed)
    return f"""<div class="gate">
<div class="gate-name">{_e(gate.get("name") or "gate")}</div>
<div class="gate-value">{value_str} <span class="chip {css}" style="font-size:10px;">{_e(word)}</span></div>
<div class="gate-source">threshold: {threshold_str} · {_e(gate.get("source") or "")}</div>
</div>"""


def _render_trend(trend: list[dict]) -> str:
    if not trend:
        return ""
    rows = []
    for row in trend:
        rows.append(
            "<tr>"
            f"<td>{_e(str(row.get('round_number', '')))}</td>"
            f"<td>{_e(str(row.get('assessment_status') or row.get('decision', '')))}</td>"
            f"<td>{_e(_fmt_pct(row.get('candidate_pass_rate')))}</td>"
            f"<td>{_e(_fmt_pct(row.get('pass_rate')))}</td>"
            f"<td>{_e(str(row.get('num_keep', '')))}</td>"
            f"<td>{_e(str(row.get('num_modify', '')))}</td>"
            f"<td>{_e(str(row.get('num_reject', '')))}</td>"
            f"<td>{_e(str(row.get('primary_bottleneck', '')))}</td>"
            f"<td>{_e(str(row.get('suggested_strategy', '')))}</td>"
            "</tr>"
        )
    return f"""<div class="section">
<h2>Cross-round trend</h2>
<table class="trend-table">
<thead><tr><th>round</th><th>assessment / legacy decision</th><th>cand. pass rate</th><th>pass rate</th><th>keep</th><th>modify</th><th>reject</th><th>bottleneck</th><th>strategy</th></tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>
</div>"""


def _render_main_viewer(main_viewer_html: str | None) -> str:
    if not main_viewer_html:
        return ""
    return f"""<div class="section">
<h2>Pocket viewer (Mol*)</h2>
<div>{main_viewer_html}</div>
</div>"""


def _render_candidates(candidates: list[dict]) -> str:
    if not candidates:
        return """<div class="section"><h2>Top candidates by cohort</h2><div class="header-sub">No candidates in candidate_ranking.json (DEC_* excluded per coordinator rules).</div></div>"""

    partitions: dict[tuple[str, str], list[dict]] = {}
    for c in candidates:
        key = (
            str(c.get("cohort_id") or "cohort:unspecified"),
            str(c.get("normalization_pool") or "nomination_candidates"),
        )
        partitions.setdefault(key, []).append(c)
    groups = []
    for (cohort_id, normalization_pool), partition in sorted(partitions.items()):
        cards = "".join(_render_candidate_card(c) for c in partition)
        groups.append(
            f'<h3>Cohort <code>{_e(cohort_id)}</code> · '
            f'<code>{_e(normalization_pool)}</code></h3>'
            f'<div class="header-sub">Showing {len(partition)} candidate(s); '
            "ranks are local to this cohort and pool.</div>"
            f"{cards}"
        )
    return f"""<div class="section">
<h2>Top candidates by cohort</h2>
<div class="header-sub">Source: <code>step6/candidate_ranking.json</code> (DEC_* excluded upstream). Overall ranks are not compared across cohorts.</div>
{"".join(groups)}
</div>"""


def _render_candidate_reviews(reviews: dict) -> str:
    groups = reviews.get("groups") or []
    if not groups:
        return ""
    sections = []
    for group in groups:
        rows = []
        for row in group.get("rows") or []:
            evidence = list(row.get("evidence_gaps") or []) + list(row.get("protocol_blockers") or [])
            cells = [
                row.get("review_position"), row.get("mol_name"),
                row.get("overall_score"), row.get("overall_rank"),
                "; ".join(row.get("molecular_findings") or []) or "No flagged category evidence recorded",
                "; ".join(evidence) or "No gaps recorded",
                f"{row.get('rank_status') or 'unreported'} / {row.get('nomination_status') or 'unreported'}",
            ]
            rows.append("<tr>" + "".join(
                f"<td>{_e(str(value)) if value is not None else '—'}</td>" for value in cells
            ) + "</tr>")
        if not rows:
            rows.append('<tr><td colspan="7">No candidates with an existing OverallScore.</td></tr>')
        sections.append(
            f"<h3>{_e(str(group['classification']).upper())} Top 1–5 · "
            f"cohort <code>{_e(group['cohort_id'])}</code> · {_e(group['normalization_pool'])}</h3>"
            f'<div class="header-sub">{group["scored_count"]} scored; '
            f'{group["unscored_count"]} without a score (excluded from this review order).</div>'
            '<div style="overflow-x:auto"><table><thead><tr>'
            '<th>Review #</th><th>Molecule</th><th>Existing OverallScore</th><th>OverallRank</th>'
            '<th>Molecular category findings</th><th>Evidence gaps / protocol blockers</th>'
            '<th>Rank status / nomination status</th></tr></thead><tbody>'
            + "".join(rows) + "</tbody></table></div>"
        )
    return (
        '<div class="section"><h2>KEEP and MODIFY Top 1–5 for review</h2>'
        f'<div class="header-sub">Source round: {_e(str(reviews.get("source_round")))}. '
        'Review only: ordered by existing OverallScore within each cohort, pool and classification. '
        'Review # is a display position, not OverallRank. Scores are not recalculated; '
        'blocked rows remain blocked and MODIFY rows are not nominations.</div>'
        + "".join(sections) + "</div>"
    )


def _render_candidate_card(c: dict) -> str:
    mol_name = _e(c.get("mol_name") or "?")
    verdict = c.get("verdict") or "unknown"
    verdict_css = _chip_class_for_verdict(verdict)
    property_verdict = c.get("property_verdict") or "unknown"
    prop_css = _chip_class_for_verdict(property_verdict)
    cohort_id = str(c.get("cohort_id") or "cohort:unspecified")

    flags = c.get("flags") or []
    prop_flags = c.get("property_flags") or []
    all_flags = list(flags) + [f"prop:{f}" for f in prop_flags]
    if c.get("nomination_status") == "blocked_by_quality_gate":
        all_flags.append("protocol-blocked")
    flag_chips = "".join(f'<span class="chip chip-warn">{_e(f)}</span>' for f in all_flags) or '<span class="chip chip-pass">no flags</span>'

    # 3D viewer html was already inlined by mol_renderer_3d
    viewer_html = c.get("viewer_3d_html") or _placeholder_html("no 3D pose available")

    # 2D structure
    if c.get("structure_2d_png"):
        structure_2d = f'<img src="{_e(c["structure_2d_png"])}" alt="{mol_name}">'
    else:
        err = c.get("structure_2d_error") or "no 2D image"
        structure_2d = f'<div class="viewer-error">{_e(err)}</div>'
    smiles = c.get("smiles") or ""
    structure_2d_block = f'<div class="structure-2d">{structure_2d}<div class="smiles">{_e(smiles)}</div></div>'

    metrics = _render_candidate_metrics(c)
    diagnosis = c.get("diagnosis") or []
    diag_html = ""
    if diagnosis:
        items = "".join(f"<li>{_e(str(d))}</li>" for d in diagnosis[:10])
        diag_html = f'<ul class="diagnosis-list">{items}</ul>'
    hint = c.get("modification_hint")
    hint_html = f'<div class="header-sub" style="margin-top: 8px;"><strong>hint:</strong> {_e(str(hint))}</div>' if hint else ""

    rank_value = c.get("rank")
    rank = _e(str(rank_value)) if rank_value is not None else "unranked"
    return f"""<div class="card">
<div class="card-header">
  <div class="card-title">#{rank} · {mol_name}
    <span class="chip">cohort: {_e(cohort_id)}</span>
    <span class="chip {verdict_css}">eval: {_e(str(verdict))}</span>
    <span class="chip {prop_css}">prop: {_e(str(property_verdict))}</span>
    {_intrinsic_chip(c)}
  </div>
  <div>{flag_chips}</div>
</div>
<div class="card-grid">
  <div>{viewer_html}</div>
  {structure_2d_block}
  <div>{metrics}{diag_html}{hint_html}</div>
</div>
</div>"""


def _intrinsic_chip(candidate: dict) -> str:
    intrinsic = candidate.get("intrinsic_classification")
    classification = candidate.get("classification")
    if not intrinsic or intrinsic == classification:
        return ""
    css = _chip_class_for_verdict(str(intrinsic))
    return f'<span class="chip {css}">intrinsic: {_e(str(intrinsic))}</span>'


def _render_candidate_metrics(c: dict) -> str:
    scores = c.get("docking_scores") or {}
    primary_score = c.get("primary_score") or {}
    bp = c.get("blueprint_match") or {}
    pc = c.get("posecheck") or {}
    pb = c.get("posebusters") or {}
    props = c.get("properties") or {}
    lip = c.get("lipinski") or {}
    admet = c.get("admet") or {}
    chemfh = c.get("chemfh") or {}

    rows: list[tuple[str, str]] = []

    def _add(label: str, value):
        if value is None or value == "":
            return
        rows.append((label, _fmt_number(value)))

    primary_metric = (
        c.get("primary_score_metric_id")
        or primary_score.get("metric_id")
    )
    if primary_score.get("value") not in (None, ""):
        unit = primary_score.get("unit")
        direction = primary_score.get("direction")
        qualifiers = ", ".join(
            str(value) for value in (unit, direction) if value not in (None, "")
        )
        label = f"Primary · {primary_metric or 'score'}"
        if qualifiers:
            label += f" ({qualifiers})"
        _add(label, primary_score.get("value"))

    _add("CNNscore", scores.get("CNNscore"))
    _add("CNNaffinity", scores.get("CNNaffinity"))
    _add("minimizedAffinity", scores.get("minimizedAffinity"))
    _add("blueprint required hit rate", bp.get("required_hit_rate"))
    if bp.get("required_hits") is not None and bp.get("required_total") is not None:
        rows.append(("blueprint hits", f'{bp.get("required_hits")}/{bp.get("required_total")}'))
    _add("clash count (heavy)", pc.get("steric_clash_count"))
    _add("strain kcal/mol", pc.get("strain_energy_kcal"))
    _add("PoseBusters violations", pb.get("num_violations"))
    _add("QED", props.get("QED") or props.get("qed"))
    _add("SA", props.get("SA_score") or props.get("sa_score"))
    _add("MW", props.get("MW"))
    _add("LogP", props.get("LogP"))
    _add(
        "Lipinski violations",
        lip.get("violations")
        if lip.get("violations") is not None
        else lip.get("num_violations"),
    )
    if isinstance(admet, dict):
        admet_summary = admet.get("summary") or {}
        _add(
            "hERG risk",
            ((admet_summary.get("hERG") or {}).get("risk"))
            or admet.get("hERG_risk")
            or admet.get("hERG"),
        )
        _add(
            "CYP3A4 inhibitor",
            ((admet_summary.get("CYP3A4_Veith") or {}).get("inhibitor"))
            or admet.get("CYP3A4_inhibitor"),
        )
        _add(
            "Ames",
            ((admet_summary.get("AMES") or {}).get("call"))
            or admet.get("Ames")
            or admet.get("AMES"),
        )
        _add(
            "oral bioavailability",
            ((admet_summary.get("Bioavailability_Ma") or {}).get("call"))
            or admet.get("oral_bioavailability"),
        )
    if isinstance(chemfh, dict):
        chemfh_summary = chemfh.get("summary") or {}
        _add("ChemFH status", chemfh.get("status"))
        _add(
            "ChemFH positive endpoints",
            chemfh_summary.get("positive_endpoint_count"),
        )
        positive_names = chemfh_summary.get("positive_endpoints") or []
        if positive_names:
            rows.append((
                "ChemFH positive calls",
                ", ".join(str(value) for value in positive_names),
            ))
        _add("ChemFH rule hits (evidence)", chemfh_summary.get("rule_hit_count"))
    synthesis = c.get("synthesis") or {}
    if isinstance(synthesis, dict):
        attempted = synthesis.get("attempted")
        rows.append((
            "AiZynthFinder attempted",
            _fmt_number(attempted) if isinstance(attempted, bool) else "unknown",
        ))
        _add("AiZynthFinder status", synthesis.get("status") or "missing")
        _add("AiZynthFinder reason", synthesis.get("reason") or "not_recorded")
        if synthesis.get("attempted_source") not in (None, "explicit"):
            _add(
                "AiZynthFinder attempted source",
                synthesis.get("attempted_source"),
            )
        _add("AiZynthFinder route count", synthesis.get("route_count"))
        _add("AiZynthFinder meaning", synthesis.get("interpretation"))
        _add(
            "AiZynthFinder contract warning",
            synthesis.get("contract_warning"),
        )

    trs = "".join(f'<tr><td class="k">{_e(k)}</td><td class="v">{_e(v)}</td></tr>' for k, v in rows)
    return f'<table class="metric-table"><tbody>{trs}</tbody></table>'


def _render_synthesis_audit(summary: dict) -> str:
    """Render AiZynthFinder execution counts apart from molecule verdicts."""
    if not summary.get("available"):
        return f"""<div class="section">
<h2>AiZynthFinder execution audit</h2>
<div class="header-sub">Unavailable: {_e(summary.get("reason") or "no Step 5 synthesis evidence")}. No attempt, failure, or molecule verdict is inferred.</div>
</div>"""

    chips = [
        f'<span class="chip">Step 5 rows: {_e(summary.get("evaluated_records"))}</span>',
        f'<span class="chip chip-pass">attempted: {_e(summary.get("attempted"))}</span>',
        f'<span class="chip">not attempted / N/A: {_e(summary.get("not_attempted"))}</span>',
        f'<span class="chip chip-warn">attempt unknown: {_e(summary.get("attempt_unknown"))}</span>',
        f'<span class="chip chip-warn">operational errors: {_e(summary.get("operational_errors"))}</span>',
    ]
    if summary.get("contract_inconsistent"):
        chips.append(
            '<span class="chip chip-fail">attempt/status conflicts: '
            f'{_e(summary.get("contract_inconsistent"))}</span>'
        )

    table_rows = []
    for row in summary.get("status_rows") or []:
        attempted = row.get("attempted")
        attempted_text = (
            _fmt_number(attempted)
            if isinstance(attempted, bool)
            else str(attempted) if attempted is not None
            else "unknown"
        )
        table_rows.append(
            "<tr>"
            f'<td>{_e(row.get("status"))}</td>'
            f'<td>{_e(attempted_text)}</td>'
            f'<td>{_e(row.get("count"))}</td>'
            f'<td>{_e(row.get("meaning"))}</td>'
            "</tr>"
        )
    table = "".join(table_rows) or (
        '<tr><td colspan="4">No per-row synthesis records.</td></tr>'
    )
    return f"""<div class="section">
<h2>AiZynthFinder execution audit</h2>
<div>{"".join(chips)}</div>
<div class="header-sub" style="margin: 8px 0;">{_e(summary.get("verdict_independence"))} <code>not_selected_top_k</code>, <code>not_eligible</code>, and <code>disabled</code> mean <code>attempted=false</code>. <code>error</code> is an operational failure with N/A metrics. <code>not_found_within_budget</code> is an attempted search outcome, not a tool failure.</div>
<table class="trend-table">
<thead><tr><th>Status</th><th>Attempted</th><th>Count</th><th>Interpretation</th></tr></thead>
<tbody>{table}</tbody>
</table>
</div>"""


def _render_qc_panels(figures: dict, plot_errors: dict) -> str:
    plot_order = [
        ("verdict_bar", "Evaluator verdicts"),
        ("blueprint_hit_hist", "Blueprint hit rate"),
        ("score_scatter", "Docking scores"),
        ("qed_sa_scatter", "QED vs SA"),
        ("redocking_gauge", "Redocking RMSD"),
        ("trend_pass_rate", "Cross-round trend"),
    ]
    cells = []
    for key, caption in plot_order:
        rel = figures.get(key)
        if not rel:
            continue
        cells.append(
            f'<div class="plot"><img src="{_e(rel)}" alt="{_e(caption)}">'
            f'<div class="plot-caption">{_e(caption)}</div></div>'
        )
    if not cells:
        return f"""<div class="section">
<h2>QC & druglikeness</h2>
<div class="header-sub">Plots disabled or matplotlib unavailable.{" errors: " + _e(json.dumps(plot_errors)) if plot_errors else ""}</div>
</div>"""
    return f"""<div class="section">
<h2>QC & druglikeness</h2>
<div class="plots-grid">{"".join(cells)}</div>
</div>"""


def _render_screening(screening: dict, figures: dict) -> str:
    if not screening.get("available"):
        return ""

    roc = figures.get("roc_curve")
    ef = figures.get("ef_bar")
    auroc = screening.get("auroc")
    ci = screening.get("auroc_ci")
    passed = screening.get("passed")

    stat_bits = [
        f'<span class="chip {_chip_class_for_bool(passed)}">screening gate: {_e(_pass_word(passed))}</span>',
        f'<span class="chip">AUROC: {_e(_fmt_number(auroc))}</span>',
    ]
    if ci:
        stat_bits.append(f'<span class="chip">95% CI: {_e(_fmt_ci(ci))}</span>')
    for pct in ("1%", "5%", "10%"):
        v = (screening.get("ef") or {}).get(pct)
        if v is not None:
            stat_bits.append(f'<span class="chip">EF{pct}: {_e(_fmt_number(v))}</span>')

    metric_rows = []
    average_precision = screening.get("average_precision")
    if average_precision is not None:
        metric_rows.append(
            f'<tr><th>Average precision (AP)</th><td>{_e(_fmt_pct(average_precision, 2))}</td></tr>'
        )
    if screening.get("top20_size") is not None:
        metric_rows.extend([
            f'<tr><th>Top-20 compounds evaluated</th><td>{_e(_fmt_number(screening.get("top20_size")))}</td></tr>',
            f'<tr><th>Top-20 active hits (observed order)</th><td>{_e(_fmt_number(screening.get("top20_actives")))}</td></tr>',
            f'<tr><th>Top-20 hit ratio (observed order)</th><td>{_e(_fmt_pct(screening.get("top20_hit_ratio"), 2))}</td></tr>',
        ])
    tie_note = ""
    if screening.get("top20_boundary_tie"):
        metric_rows.extend([
            f'<tr><th>Top-20 active hits (tie-aware expectation)</th><td>{_e(_fmt_number(screening.get("top20_expected_actives")))}</td></tr>',
            f'<tr><th>Top-20 hit ratio (tie-aware expectation)</th><td>{_e(_fmt_pct(screening.get("top20_expected_hit_ratio"), 2))}</td></tr>',
            f'<tr><th>Possible Top-20 active hits across the cutoff tie</th><td>{_e(_fmt_ci(screening.get("top20_active_count_range")))}</td></tr>',
        ])
        tie_note = (
            '<p class="header-sub">Equal scores cross the Top-20 cutoff. Observed hits '
            'retain input order; expected hits assume uniform ordering within that tie. '
            'The range covers every possible ordering. This is not a reconstruction '
            "of a unique ordering within tied scores.</p>"
        )
    metrics_table = (
        '<table class="metric-table"><tbody>' + "".join(metric_rows) + '</tbody></table>'
        if metric_rows else ""
    )

    plots = []
    if roc:
        plots.append(f'<div class="plot"><img src="{_e(roc)}" alt="ROC"><div class="plot-caption">ROC curve</div></div>')
    if ef:
        plots.append(f'<div class="plot"><img src="{_e(ef)}" alt="EF"><div class="plot-caption">Enrichment factors</div></div>')

    return f"""<div class="section">
<h2>Retrospective screening (ACT_/DEC_)</h2>
<div class="header-sub">These metrics validate the docking + evaluation protocol; DEC_* records are excluded from candidate nomination.</div>
<div style="margin-top: 8px;">{"".join(stat_bits)}</div>
{metrics_table}
{tie_note}
<div class="plots-grid" style="margin-top: 12px;">{"".join(plots)}</div>
</div>"""


def _render_footer(report_data: dict, asset_manifest: dict, plots_status: dict, missing: list) -> str:
    generated_at = _e(report_data.get("generated_at") or "")
    asset_bits = []
    for name, meta in ((asset_manifest.get("assets") or {}).items()):
        status = "present" if meta.get("present") else "missing"
        asset_bits.append(f"{name}={status}")
    asset_line = ", ".join(asset_bits) or "no assets tracked"

    plot_avail = "yes" if plots_status.get("available") else f"no ({plots_status.get('reason')})"

    missing_block = ""
    if missing:
        items = []
        for m in missing:
            step = m.get("step", "?")
            reason = m.get("reason", "")
            items.append(f"<li>{_e(f'{step}: {reason}')}</li>")
        rows = "".join(items)
        missing_block = f'<div class="missing-list"><strong>missing sources</strong><ul>{rows}</ul></div>'

    return f"""<div class="section footer">
<div>Generated at {generated_at}</div>
<div>Assets: {_e(asset_line)}</div>
<div>Plots: {_e(plot_avail)}</div>
{missing_block}
</div>"""


def _body_scripts(project_root: Path, html_path: Path, asset_manifest: dict) -> str:
    """Load 3Dmol.js and Mol* bundles once at the end of body."""
    has_3dmol = ((asset_manifest.get("assets") or {}).get("3Dmol-min.js") or {}).get("present")
    has_molstar = ((asset_manifest.get("assets") or {}).get("molstar.js") or {}).get("present")
    scripts = []
    if has_3dmol:
        rel = _relpath(project_root / "assets" / "3Dmol-min.js", html_path.parent)
        scripts.append(f'<script src="{_e(rel)}"></script>')
    if has_molstar:
        rel = _relpath(project_root / "assets" / "molstar.js", html_path.parent)
        scripts.append(f'<script src="{_e(rel)}"></script>')
    return "\n".join(scripts)


# ---------- utility ----------

def _e(value: Any) -> str:
    return _html.escape(str(value) if value is not None else "", quote=True)


def _relpath(target: Path, start: Path) -> str:
    import os
    try:
        return os.path.relpath(str(Path(target).resolve()), start=str(Path(start).resolve()))
    except ValueError:
        return str(target)


def _fmt_number(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int,)):
        return str(value)
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if abs(f) >= 100:
        return f"{f:.1f}"
    if abs(f) >= 1:
        return f"{f:.2f}"
    return f"{f:.3f}"


def _fmt_value(value, unit: str | None) -> str:
    if value is None:
        return "—"
    n = _fmt_number(value)
    return f"{n} {unit}" if unit and unit not in n else n


def _fmt_pct(value, decimals: int = 1) -> str:
    if value is None:
        return "—"
    try:
        return f"{float(value) * 100:.{decimals}f}%"
    except (TypeError, ValueError):
        return str(value)


def _fmt_ci(ci) -> str:
    if isinstance(ci, (list, tuple)) and len(ci) == 2:
        return f"{_fmt_number(ci[0])} – {_fmt_number(ci[1])}"
    return str(ci)


def _chip_class_for_verdict(verdict: str | None) -> str:
    if not verdict:
        return "chip-unknown"
    v = str(verdict).lower()
    if v in ("pass", "keep"):
        return "chip-pass"
    if v in ("warn", "flag", "modify"):
        return "chip-warn"
    if v in ("fail", "reject"):
        return "chip-fail"
    return "chip-unknown"


def _chip_class_for_bool(value) -> str:
    if value is True:
        return "chip-pass"
    if value is False:
        return "chip-fail"
    return "chip-unknown"


def _chip_class_for_outcome(header: dict) -> str:
    """Protocol success or legacy termination alone cannot establish success."""
    if header.get("assessment_schema_version") == "round-assessment/v1":
        status = (header.get("round_assessment") or {}).get("status")
        return {
            "valid": "chip-pass", "review_required": "chip-warn",
            "invalid": "chip-fail", "not_evaluated": "chip-unknown",
        }.get(status, "chip-unknown")
    outcome = str(header.get("decision_outcome") or "").lower()
    if outcome == "invalid" or header.get("decision") == "terminate_invalid":
        return "chip-fail"
    if outcome == "success":
        return (
            "chip-pass"
            if header.get("goal_met") is True
            and header.get("valid_for_nomination") is True
            and header.get("quality_gate_passed") is True
            else "chip-warn"
        )
    if outcome in ("iterate", "max_rounds", "no_viable") or header.get("decision") == "iterate":
        return "chip-warn"
    if header.get("goal_met") is False or header.get("valid_for_nomination") is False:
        return "chip-warn"
    return "chip-unknown"


def _pass_word(value) -> str:
    if value is True:
        return "passed"
    if value is False:
        return "failed"
    return "unknown"


def _placeholder_html(message: str) -> str:
    return f'<div class="viewer-3dmol viewer-placeholder"><div class="viewer-error">{_e(message)}</div></div>'

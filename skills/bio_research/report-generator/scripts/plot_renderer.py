"""QC / druglikeness / ROC-EF plots via matplotlib.

All charts are optional — graceful-skip when matplotlib is unavailable.
Each writes a PNG to `{output_dir}/figures/` and returns the relative href
that the HTML template embeds via `<img src="figures/...">`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def check_matplotlib_available() -> bool:
    try:
        import matplotlib  # noqa: F401
        return True
    except ImportError:
        return False


def render_all(
    report_data: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    """Render every enabled chart. Returns a status dict:

    {
      "available": bool,
      "reason": str | None,
      "figures": {name: relative_path or None},
      "errors": {name: error_message},
    }
    """
    figures: dict[str, str | None] = {}
    errors: dict[str, str] = {}

    if not check_matplotlib_available():
        return {"available": False, "reason": "matplotlib not installed", "figures": {}, "errors": {}}

    import matplotlib
    matplotlib.use("Agg")  # noqa: E402 — must be set before pyplot import
    import matplotlib.pyplot as plt  # noqa: E402

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    _wrap(plt, figures, errors, figures_dir, "verdict_bar", lambda: _plot_verdict_bar(report_data, figures_dir))
    _wrap(plt, figures, errors, figures_dir, "qed_sa_scatter", lambda: _plot_qed_sa_scatter(report_data, figures_dir))
    _wrap(plt, figures, errors, figures_dir, "score_scatter", lambda: _plot_score_scatter(report_data, figures_dir))
    _wrap(plt, figures, errors, figures_dir, "blueprint_hit_hist", lambda: _plot_blueprint_hit_histogram(report_data, figures_dir))
    _wrap(plt, figures, errors, figures_dir, "redocking_gauge", lambda: _plot_redocking_gauge(report_data, figures_dir))

    # Cross-round trend + ROC only meaningful in some modes/data
    if report_data.get("mode") == "final" and report_data.get("trend"):
        _wrap(plt, figures, errors, figures_dir, "trend_pass_rate", lambda: _plot_trend(report_data, figures_dir))
    if (report_data.get("screening") or {}).get("available"):
        _wrap(plt, figures, errors, figures_dir, "roc_curve", lambda: _plot_roc(report_data, figures_dir))
        _wrap(plt, figures, errors, figures_dir, "ef_bar", lambda: _plot_ef_bars(report_data, figures_dir))

    return {"available": True, "figures": figures, "errors": errors}


# ---------- individual plots ----------

def _plot_verdict_bar(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    step4 = report_data.get("step4") or {}
    evaluation = step4.get("evaluation_report") or {}
    summary = evaluation.get("summary") or {}
    if not summary:
        return None

    labels = ["pass", "warn", "fail", "not_evaluated"]
    values = [int(summary.get(k, 0) or 0) for k in labels]
    if sum(values) == 0:
        return None

    fig, ax = plt.subplots(figsize=(5, 3), dpi=140)
    colors = ["#2ca02c", "#ff7f0e", "#d62728", "#7f7f7f"]
    ax.bar(labels, values, color=colors)
    ax.set_ylabel("count")
    ax.set_title("Evaluator verdicts")
    for i, v in enumerate(values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    _tight_save(fig, figures_dir / "verdict_bar.png")
    return "figures/verdict_bar.png"


def _plot_qed_sa_scatter(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    step5 = report_data.get("step5") or {}
    results = (step5.get("property_filter_report") or {}).get("results") or []
    if not results:
        return None

    xs, ys, colors, labels = [], [], [], []
    verdict_color = {"pass": "#2ca02c", "flag": "#ff7f0e", "fail": "#d62728"}
    for r in results:
        props = (r.get("properties") or {})
        qed = props.get("QED") or props.get("qed")
        sa = props.get("SA_score") or props.get("sa_score") or props.get("SA")
        if qed is None or sa is None:
            continue
        try:
            xs.append(float(qed))
            ys.append(float(sa))
        except (TypeError, ValueError):
            continue
        colors.append(verdict_color.get(r.get("verdict"), "#7f7f7f"))
        labels.append(str(r.get("mol_name") or ""))

    if not xs:
        return None

    fig, ax = plt.subplots(figsize=(5.2, 4), dpi=140)
    ax.scatter(xs, ys, c=colors, s=32, alpha=0.75, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("QED (higher is better)")
    ax.set_ylabel("SA score (lower is better)")
    ax.set_title("Drug-likeness landscape")
    ax.axvline(0.4, linestyle="--", color="gray", alpha=0.6, label="QED 0.4")
    ax.axhline(5.0, linestyle="--", color="gray", alpha=0.6, label="SA 5.0")
    ax.legend(fontsize=8, loc="best")
    _tight_save(fig, figures_dir / "qed_sa_scatter.png")
    return "figures/qed_sa_scatter.png"


def _plot_score_scatter(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    scores = (report_data.get("step3") or {}).get("docking_scores") or []
    if not scores:
        return None

    step4 = report_data.get("step4") or {}
    per_mol = step4.get("per_molecule") or {}

    xs, ys, colors = [], [], []
    verdict_color = {"pass": "#2ca02c", "warn": "#ff7f0e", "fail": "#d62728"}
    for s in scores:
        try:
            aff = float(s.get("CNNaffinity"))
            energy = float(s.get("minimizedAffinity"))
        except (TypeError, ValueError):
            continue
        mol_name = s.get("mol_name") or s.get("name")
        verdict = (per_mol.get(str(mol_name)) or {}).get("verdict") if mol_name else None
        color = verdict_color.get(verdict, "#1f77b4")
        # De-emphasize DEC_
        if isinstance(mol_name, str) and mol_name.startswith("DEC_"):
            color = "#c7c7c7"
        xs.append(aff)
        ys.append(energy)
        colors.append(color)

    if not xs:
        return None

    fig, ax = plt.subplots(figsize=(5.2, 4), dpi=140)
    ax.scatter(xs, ys, c=colors, s=22, alpha=0.7, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("CNNaffinity (higher = stronger predicted binding)")
    ax.set_ylabel("minimizedAffinity kcal/mol (lower is better)")
    ax.set_title("GNINA score landscape")
    _tight_save(fig, figures_dir / "score_scatter.png")
    return "figures/score_scatter.png"


def _plot_blueprint_hit_histogram(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    step4 = report_data.get("step4") or {}
    evaluation = step4.get("evaluation_report") or {}
    molecules = (
        evaluation.get("evaluation")
        or evaluation.get("molecules")
        or evaluation.get("results")
        or []
    )
    if not molecules:
        return None

    rates_by_verdict: dict[str, list[float]] = {"pass": [], "warn": [], "fail": []}
    for m in molecules:
        bp = m.get("blueprint_match") or {}
        rate = bp.get("required_hit_rate")
        verdict = m.get("verdict")
        if rate is None or verdict not in rates_by_verdict:
            continue
        try:
            rates_by_verdict[verdict].append(float(rate))
        except (TypeError, ValueError):
            continue

    if not any(rates_by_verdict.values()):
        return None

    fig, ax = plt.subplots(figsize=(5.2, 3.5), dpi=140)
    colors = {"pass": "#2ca02c", "warn": "#ff7f0e", "fail": "#d62728"}
    bins = [i / 10 for i in range(0, 11)]
    for verdict in ("pass", "warn", "fail"):
        if rates_by_verdict[verdict]:
            ax.hist(
                rates_by_verdict[verdict],
                bins=bins,
                alpha=0.6,
                label=f"{verdict} (n={len(rates_by_verdict[verdict])})",
                color=colors[verdict],
                edgecolor="black",
                linewidth=0.4,
            )
    ax.axvline(0.6, linestyle="--", color="gray", alpha=0.7, label="pass threshold 0.6")
    ax.axvline(0.4, linestyle=":", color="gray", alpha=0.7, label="warn threshold 0.4")
    ax.set_xlabel("blueprint required_hit_rate")
    ax.set_ylabel("count")
    ax.set_title("Blueprint match distribution")
    ax.legend(fontsize=8)
    _tight_save(fig, figures_dir / "blueprint_hit_hist.png")
    return "figures/blueprint_hit_hist.png"


def _plot_redocking_gauge(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    redocking = ((report_data.get("step3") or {}).get("redocking_validation")) or {}
    # A redocking check that was not performed has no RMSD to chart.
    if redocking.get("performed") is False:
        return None
    rmsd_raw = redocking.get("best_rmsd_angstrom", redocking.get("rmsd"))
    try:
        rmsd = float(rmsd_raw)
    except (TypeError, ValueError):
        return None
    try:
        threshold = float(redocking.get("rmsd_threshold_angstrom", redocking.get("threshold", 2.0)))
    except (TypeError, ValueError):
        threshold = 2.0

    fig, ax = plt.subplots(figsize=(5.2, 1.8), dpi=140)
    ax.barh(["redocking RMSD"], [rmsd], color=("#2ca02c" if rmsd <= threshold else "#d62728"))
    ax.axvline(threshold, linestyle="--", color="gray", alpha=0.7, label=f"threshold {threshold} Å")
    ax.set_xlim(0, max(threshold * 2, rmsd * 1.3, 3.0))
    ax.set_xlabel("RMSD (Å) — lower is better")
    ax.set_title("Native redocking")
    ax.legend(fontsize=8, loc="lower right")
    _tight_save(fig, figures_dir / "redocking_gauge.png")
    return "figures/redocking_gauge.png"


def _plot_trend(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    trend = report_data.get("trend") or []
    if len(trend) < 1:
        return None

    xs, pass_rate, cand_rate, keep = [], [], [], []
    bottlenecks: list[tuple[int, str]] = []
    for row in trend:
        rn = row.get("round_number")
        if rn is None:
            continue
        xs.append(rn)
        pr = row.get("pass_rate")
        cr = row.get("candidate_pass_rate")
        try:
            pass_rate.append(float(pr) if pr is not None else None)
        except (TypeError, ValueError):
            pass_rate.append(None)
        try:
            cand_rate.append(float(cr) if cr is not None else None)
        except (TypeError, ValueError):
            cand_rate.append(None)
        try:
            keep.append(int(row.get("num_keep") or 0))
        except (TypeError, ValueError):
            keep.append(0)
        if row.get("primary_bottleneck"):
            bottlenecks.append((rn, str(row["primary_bottleneck"])))

    if not xs:
        return None

    fig, ax1 = plt.subplots(figsize=(6.5, 3.8), dpi=140)
    if any(v is not None for v in pass_rate):
        ax1.plot(xs, [v if v is not None else 0 for v in pass_rate], marker="o", label="pass_rate", color="#1f77b4")
    if any(v is not None for v in cand_rate):
        ax1.plot(xs, [v if v is not None else 0 for v in cand_rate], marker="s", label="candidate_pass_rate", color="#2ca02c")
    ax1.set_xlabel("Round")
    ax1.set_ylabel("rate")
    ax1.set_ylim(0, 1.05)
    ax1.legend(loc="upper left", fontsize=8)

    ax2 = ax1.twinx()
    ax2.bar(xs, keep, alpha=0.2, color="#7f7f7f", label="num_keep")
    ax2.set_ylabel("num_keep")

    for rn, b in bottlenecks:
        ax1.annotate(b, xy=(rn, 0.05), fontsize=7, rotation=25, color="#d62728")

    ax1.set_title("Cross-round trend")
    _tight_save(fig, figures_dir / "trend_pass_rate.png")
    return "figures/trend_pass_rate.png"


def _plot_roc(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    screening = report_data.get("screening") or {}
    roc = screening.get("roc_curve") or {}
    fpr = roc.get("fpr") or []
    tpr = roc.get("tpr") or []
    if not fpr or not tpr or len(fpr) != len(tpr):
        return None

    fig, ax = plt.subplots(figsize=(4.5, 4.5), dpi=140)
    ax.plot(fpr, tpr, color="#1f77b4", label=f"AUROC = {screening.get('auroc')}")
    ax.plot([0, 1], [0, 1], "--", color="gray", alpha=0.6)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Retrospective screening ROC")
    ax.legend(fontsize=8)
    _tight_save(fig, figures_dir / "roc_curve.png")
    return "figures/roc_curve.png"


def _plot_ef_bars(report_data: dict, figures_dir: Path) -> str | None:
    import matplotlib.pyplot as plt

    ef = (report_data.get("screening") or {}).get("ef") or {}
    labels, values = [], []
    for pct in ("1%", "5%", "10%"):
        v = ef.get(pct)
        try:
            if v is not None:
                labels.append(f"EF{pct}")
                values.append(float(v))
        except (TypeError, ValueError):
            continue
    if not labels:
        return None

    fig, ax = plt.subplots(figsize=(4, 3), dpi=140)
    ax.bar(labels, values, color="#2ca02c")
    for i, v in enumerate(values):
        ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Enrichment factor")
    ax.set_title("Enrichment (higher is better)")
    _tight_save(fig, figures_dir / "ef_bar.png")
    return "figures/ef_bar.png"


# ---------- utility ----------

def _wrap(plt, figures: dict, errors: dict, figures_dir: Path, name: str, fn) -> None:
    try:
        rel = fn()
        figures[name] = rel
    except Exception as e:  # noqa: BLE001 — plots must never crash the report
        figures[name] = None
        errors[name] = str(e)
    finally:
        plt.close("all")


def _tight_save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")

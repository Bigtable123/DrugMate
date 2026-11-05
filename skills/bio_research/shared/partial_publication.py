"""Publish durable execution snapshots without requiring a scored cohort.

The runner manifest determines which results belong to this invocation. Never
discover scientific results by scanning step directories: they may contain an
invalidated generation after a rewind. Raw files are references with observed
size/existence, not scientifically validated or source-sealed deliverables.
"""
from __future__ import annotations

import importlib.util
import sys
from html import escape
from pathlib import Path
from typing import Any, Iterable, Mapping
from urllib.parse import quote

import output_layout as ol


STAGE_SOURCES = {
    ol.STRUCTURE_PREP_DIR: ("structure_prep",),
    ol.POCKET_ANALYSIS_DIR: ("pocket_analysis",),
    ol.MOLECULE_GENERATION_DIR: ("molecule_generation",),
    ol.SM_PROPERTY_DIR: ("property_prefilter", "property_filter"),
    ol.DOCKING_DIR: ("docking",),
    ol.POST_DOCKING_DIR: ("pose_evaluation",),
    ol.FINAL_REPORT_DIR: ("coordination", "reporting"),
}
INDEX_COLUMNS = (
    "pipeline_step_id", "artifact_name", "path", "kind", "exists", "bytes", "error",
)
TERMINAL_STATES = {
    "done", "failed", "failed_best_effort", "interrupted", "cancelled", "resource_unavailable",
}


def _write_execution_report(
    layout: ol.RunLayout,
    snapshot: Mapping[str, Any],
    step_artifacts: Mapping[str, list[dict]],
) -> dict[str, str]:
    """Write an offline report even when no molecule reached scoring.

    Reuse the artifact observations already made for the indexes: rendering a
    report must not reread large scientific files or run another report skill.
    """
    steps = snapshot.get("steps") or {}
    rows = []
    for directory, step_ids in STAGE_SOURCES.items():
        for step_id in step_ids:
            info = steps.get(step_id) or {}
            rows.append({
                "step_id": step_id, "status": info.get("status", "not_run"),
                "error": info.get("error"),
                "artifact_index": f"{directory}/{ol.ARTIFACT_INDEX}",
                "artifacts": step_artifacts.get(step_id, []),
            })
    completed = [row["step_id"] for row in rows if row["status"] == "done"]
    failed = [row["step_id"] for row in rows if row["status"] in TERMINAL_STATES - {"done"}]
    unfinished = [row["step_id"] for row in rows if row["status"] not in TERMINAL_STATES | {"disabled"}]
    artifacts = [item for row in rows for item in row["artifacts"]]
    available = [item for item in artifacts if item["exists"] and item["kind"] == "file"]
    report = {
        "report_kind": "execution_summary", "publication_mode": "pipeline_partial",
        "finalized": False, "generated_at": layout.generated_at,
        "execution_status": snapshot.get("status", "unknown"), "error": snapshot.get("error"),
        "scientific_status": snapshot.get("scientific_status", "not_evaluated"),
        "source_integrity": "unverified_references", "path_base": "published_root",
        "completed_steps": completed, "failed_steps": failed, "unfinished_steps": unfinished,
        "available_file_count": len({item["path"] for item in available}),
        "missing_artifacts": [item for item in artifacts if not item["exists"]],
        "steps": rows,
    }

    def md(value: Any) -> str:
        text = escape(str(value or "")).replace("\n", " ").replace("\r", " ")
        for character in ("\\", "|", "[", "]", "`", "*", "_"):
            text = text.replace(character, "\\" + character)
        return text

    def href(path: str) -> str:
        # Index paths are relative to published/, reports live one level down.
        return quote("../" + path, safe="/")

    summary = (
        f"{len(completed)} stages completed, {len(failed)} failed or interrupted, "
        f"and {len(unfinished)} incomplete or not run; {report['available_file_count']} existing files retained."
    )
    caveat = (
        "This report summarizes the results obtained and the reason execution stopped. "
        "Full scoring and final source verification are incomplete, so this report does not establish "
        "experiment success or eligibility for candidate nomination. File links point to the working "
        "directory; retain it alongside published/."
    )
    markdown = [
        "# Experiment report (partial results)", "",
        f"Execution status: {md(report['execution_status'])}.", "", summary, "", caveat, "",
    ]
    html = [
        '<!doctype html><html lang="en"><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        '<title>Experiment report (partial results)</title>',
        '<style>body{font:16px/1.6 system-ui,sans-serif;max-width:1100px;margin:2em auto;'
        'padding:0 1em}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ddd;'
        'padding:.5em;text-align:left;overflow-wrap:anywhere}code{overflow-wrap:anywhere}'
        '.notice{background:#fff6df;padding:1em}</style><body>',
        '<h1>Experiment report (partial results)</h1>',
        f"<p>Execution status: <code>{escape(str(report['execution_status']))}</code></p>",
        f"<p>{summary}</p><p class=\"notice\">{caveat}</p>",
    ]
    if report["error"]:
        markdown.extend(["## Reason execution stopped", "", md(report["error"]), ""])
        html.append(f"<h2>Reason execution stopped</h2><p>{escape(str(report['error']))}</p>")
    if not available:
        message = (
            "No registered stage artifacts are currently available. Inputs, logs and failure "
            "diagnostics can still be inspected through the indexes below."
        )
        markdown.extend([message, ""])
        html.append(f"<p>{message}</p>")
    markdown.extend([
        "[Log index](../00_inputs/log_index.tsv) · [Inputs and run control files](../00_inputs/artifact_index.tsv)",
        "", "## Stage results", "", "| Stage | Status | Existing files | Error |", "|---|---|---:|---|",
    ])
    html.append('<p><a href="../00_inputs/log_index.tsv">Log index</a> · '
                '<a href="../00_inputs/artifact_index.tsv">Inputs and run control files</a></p>'
                '<h2>Stage results</h2><table><tr><th>Stage</th><th>Status</th><th>Existing files</th><th>Error</th></tr>')
    for row in rows:
        count = len({item["path"] for item in row["artifacts"] if item["exists"] and item["kind"] == "file"})
        link = href(row["artifact_index"])
        markdown.append(f"| [{row['step_id']}]({link}) | {md(row['status'])} | {count} | {md(row['error'])} |")
        html.append(f'<tr><td><a href="{link}">{row["step_id"]}</a></td>'
                    f'<td>{escape(str(row["status"]))}</td><td>{count}</td>'
                    f'<td>{escape(str(row["error"] or ""))}</td></tr>')
    markdown.extend(["", "## Registered artifacts", "", "| Stage | Artifact | File status |", "|---|---|---|"])
    html.append('</table><h2>Registered artifacts</h2><table><tr><th>Stage</th><th>Artifact</th><th>File status</th></tr>')
    for item in artifacts:
        observed = "Missing"
        if item["exists"]:
            observed = "Directory" if item["kind"] == "directory" else f"{item['bytes']} bytes"
        name = f"{item['artifact_name']} ({item['path']})"
        markdown.append(f"| {item['pipeline_step_id']} | [{md(name)}]({href(item['path'])}) | {observed} |")
        html.append(f'<tr><td>{item["pipeline_step_id"]}</td>'
                    f'<td><a href="{href(item["path"])}">{escape(name)}</a></td><td>{observed}</td></tr>')
    html.append("</table></body></html>\n")
    written = {}
    for name, content, role in (
        ("final_report.md", "\n".join(markdown) + "\n", "execution_report_markdown"),
        ("final_report.html", "\n".join(html), "execution_report_html"),
    ):
        written[role] = str(layout.publish_text(layout.final_dir() / name, content, role=role))
    written["execution_report_data"] = str(layout.publish_json(
        layout.final_dir() / "final_report_data.json", report, role="execution_report_data",
        derived_from=[f"{directory}/{ol.ARTIFACT_INDEX}" for directory in STAGE_SOURCES],
    ))
    return written


def _file_leaves(value: Any, prefix: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(value, Mapping):
        for key, child in value.items():
            if str(key).startswith("layout_") or key in {"published_dir", "user_output_dir"}:
                continue
            yield from _file_leaves(child, f"{prefix}.{key}" if prefix else str(key))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            yield from _file_leaves(child, f"{prefix}[{index}]")
    elif isinstance(value, (str, Path)) and str(value).strip() and "://" not in str(value):
        yield prefix, str(value)


def render_snapshot(
    root: Path,
    *,
    working_root: Path,
    snapshot: Mapping[str, Any],
    session_id: str,
    publication_mode: str = "pipeline_partial",
) -> dict[str, str]:
    if publication_mode != "pipeline_partial":
        raise ValueError("execution snapshots must use pipeline_partial mode")
    working_root = Path(working_root).resolve()
    layout = ol.RunLayout(root=root, session_id=session_id)
    steps = snapshot.get("steps") or {}
    step_artifacts: dict[str, list[dict]] = {}

    def artifact(step_id: str, name: str, value: str | Path) -> dict:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = working_root / path
        row = {
            "pipeline_step_id": step_id, "artifact_name": name,
            "path": layout.relative(path), "kind": "file", "exists": False,
            "bytes": "", "error": "",
        }
        try:
            row["exists"] = path.exists()
            if path.is_dir():
                row["kind"] = "directory"
            elif path.is_file():
                row["bytes"] = path.stat().st_size
        except OSError as exc:
            row["error"] = str(exc)
        return row

    def stage(directory: str, step_ids: tuple[str, ...]) -> None:
        rows, sources = [], []
        for step_id in step_ids:
            info = steps.get(step_id) or {}
            state = info.get("status", "not_run")
            source = {"pipeline_step_id": step_id, "status": state, "error": info.get("error")}
            # Pending/disabled steps may still have old files on disk. Their
            # mere existence must never promote them to current results.
            if step_id in step_artifacts:
                rows.extend(step_artifacts[step_id])
                sources.append(source)
                continue
            first_row = len(rows)
            if state in TERMINAL_STATES:
                if info.get("result_json"):
                    rows.append(artifact(step_id, "result_json", info["result_json"]))
                rows.extend(artifact(step_id, name, path) for name, path in _file_leaves(info.get("files")))
                skill = info.get("skill")
                if skill:
                    for suffix in (".stdout.log", ".stderr.log", "_input.json"):
                        path = working_root / "logs" / f"{step_id}_{skill}{suffix}"
                        if path.is_file():
                            rows.append(artifact(step_id, path.name, path))
            step_artifacts[step_id] = rows[first_row:]
            sources.append(source)
        index = layout.publish_tsv(
            layout.path(directory, ol.ARTIFACT_INDEX), INDEX_COLUMNS, rows,
            role="artifact_index", derived_from=[row["path"] for row in rows],
        )
        states = {source["status"] for source in sources}
        state = next(iter(states)) if len(states) == 1 else "partial"
        layout.publish_json(
            layout.path(directory, ol.STEP_MANIFEST),
            {
                "publication_mode": publication_mode, "path_base": "published_root",
                "status": state, "source_results": sources,
                "artifact_index": layout.relative(index), "raw_artifact_policy": "reference_only",
            }, role="stage_manifest", derived_from=[layout.relative(index)],
        )

    controls = []
    for name in ("pipeline_config.json", "resolved_run_plan.json", "pipeline_checkpoint.json",
                 "manifest.json", "scientific_outcome.json", "tool_runtime_report.json"):
        path = working_root / name
        if path.is_file():
            controls.append(artifact("pipeline", name, path))
    layout.publish_tsv(
        layout.path(ol.INPUTS_DIR, ol.ARTIFACT_INDEX), INDEX_COLUMNS, controls,
        role="artifact_index", derived_from=[row["path"] for row in controls],
    )
    layout.publish_json(
        layout.path(ol.INPUTS_DIR, ol.INPUT_MANIFEST),
        {"publication_mode": publication_mode, "path_base": "published_root",
         "raw_artifact_policy": "reference_only", "execution_status": snapshot.get("status")},
        role="input_manifest",
    )
    logs = [
        artifact("pipeline", path.name, path)
        for path in sorted((working_root / "logs").glob("*"))
        if path.is_file() and path.suffix in {".log", ".md"}
    ]
    layout.publish_tsv(
        layout.path(ol.INPUTS_DIR, "log_index.tsv"), INDEX_COLUMNS, logs,
        role="log_index", derived_from=[row["path"] for row in logs],
    )
    for directory, step_ids in STAGE_SOURCES.items():
        stage(directory, step_ids)
    # Keep the two property populations separate, as in the full publication.
    stage(f"{ol.SM_PROPERTY_DIR}/00_prefilter", ("property_prefilter",))

    status = {
        "publication_mode": publication_mode, "finalized": False,
        "generated_at": layout.generated_at, "execution_status": snapshot.get("status"),
        "error": snapshot.get("error"),
        "scientific_status": snapshot.get("scientific_status", "not_evaluated"),
        "source_integrity": "unverified_references",
        "steps": {
            name: {key: info.get(key) for key in ("status", "error")}
            for name, info in steps.items()
        },
    }
    status_path = layout.publish_json(
        layout.final_dir() / "publication_status.json", status, role="execution_snapshot",
    )
    overview = [
        "# Partial experiment results", "",
        f"Execution status: `{snapshot.get('status', 'unknown')}`.", "",
        "This snapshot contains existing artifacts and error details, not a complete deliverable. "
        "The published view does not change the experiment's success or failure status.",
        "File indexes record existence and size at snapshot time. Original files remain in the "
        "working directory and have not passed final source-integrity verification.",
        "Index paths are relative to published/. Old artifacts from steps that were not run or "
        "were invalidated by a rewind are not listed.", "",
        "[Experiment report (HTML)](final_report.html) · [Experiment report (Markdown)](final_report.md)", "",
        "[Log index](../00_inputs/log_index.tsv) · [Inputs and run control files](../00_inputs/artifact_index.tsv)", "",
        "| Stage | Status | Error |", "|---|---|---|",
    ]
    if snapshot.get("error"):
        overview[3:3] = [str(snapshot["error"]), ""]
    for directory, step_ids in STAGE_SOURCES.items():
        for step_id in step_ids:
            info = steps.get(step_id) or {}
            error = str(info.get("error") or "").replace("|", "\\|").replace("\n", " ")
            overview.append(
                f"| [{step_id}](../{directory}/artifact_index.tsv) | {info.get('status', 'not_run')} | {error} |"
            )
    readme = layout.publish_text(
        layout.final_dir() / "README.md", "\n".join(overview) + "\n", role="execution_overview",
    )
    reports = _write_execution_report(layout, snapshot, step_artifacts)
    manifest, _ = layout.finalize(extra={
        "publication_mode": publication_mode, "execution_status": snapshot.get("status"),
        "source_root": layout.relative(working_root),
        "source_integrity": "unverified_references", "finalization": {"status": "partial"},
    })
    return {"run_manifest": str(manifest), "publication_status": str(status_path), "overview": str(readme), **reports}


def publish_snapshot(
    working_root: Path, snapshot: Mapping[str, Any], *, publication_dir: Path | None = None,
) -> dict[str, str]:
    """Use the same verified directory swap as the full cohort publisher."""
    module_name = "_sbdd_snapshot_publication_writer"
    writer = sys.modules.get(module_name)
    if writer is None:
        path = Path(__file__).resolve().parents[1] / "reflective-coordinator/scripts/publish_layout.py"
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot load publication writer: {path}")
        writer = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = writer
        spec.loader.exec_module(writer)
    return writer.publish_run_transactionally(
        publication_dir or Path(working_root) / ol.PUBLISHED_DIR, renderer=render_snapshot,
        working_root=working_root, snapshot=snapshot, session_id=Path(working_root).name,
        publication_mode="pipeline_partial",
    )

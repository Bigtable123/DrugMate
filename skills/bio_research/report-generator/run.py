#!/usr/bin/env python3
"""Report Generator Skill — turn Step 0-6 artifacts into a self-contained HTML.

Orchestrates: aggregator → asset_manager → mol_renderer_2d → mol_renderer_3d →
plot_renderer → pymol_renderer (optional) → html_template.

Every renderer graceful-skips on missing dependencies; the aggregator and HTML
writer are pure stdlib and always run.

Usage:
    python run.py --input '{"mode":"round","project_dir":"..."}'  --output-dir ...
    python run.py --input-file params.json --output-dir ...
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from aggregator import aggregate_round, aggregate_final  # noqa: E402
from asset_manager import ensure_assets, asset_available  # noqa: E402
import mol_renderer_2d  # noqa: E402
import mol_renderer_3d  # noqa: E402
import plot_renderer  # noqa: E402
import pymol_renderer  # noqa: E402
import html_template  # noqa: E402


SCIENTIFIC_COHORT_REPORT_FIELDS = (
    "cohort_id",
    "protein_state_id",
    "pocket_id",
    "round_id",
    "protocol_profile_id",
    "metric_catalog_version",
    "cohort_identity_complete",
    "cohort_identity_missing_fields",
    "cohort_identity_derived_fields",
    "cohort_provenance_sources",
)


def _artifact_identity(path: Path, *, role: str) -> dict[str, Any]:
    """Return an immutable identity for one report-owned file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "role": role,
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _write_report_artifact_manifest(
    *,
    output_dir: Path,
    project_root: Path,
    primary_files: dict[str, Path],
    asset_manifest: dict[str, Any],
    plots_status: dict[str, Any],
    pymol_cover_rel: str | None,
) -> Path:
    """Seal every file referenced by the generated HTML/report payload.

    The pipeline fingerprints this manifest and independently revalidates each
    member before reusing a cached reporting result.  This prevents an edited
    plot, viewer notice or shared JavaScript asset from hiding behind an unchanged
    ``report.html``/``report_data.json`` pair.
    """
    candidates: list[tuple[str, Path]] = list(primary_files.items())
    for section, role in (("assets", "asset"), ("licenses", "license")):
        for filename, record in sorted((asset_manifest.get(section) or {}).items()):
            if isinstance(record, dict) and record.get("present") and record.get("path"):
                asset_path = Path(str(record["path"])).expanduser().resolve()
                if not asset_path.is_file():
                    raise RuntimeError(
                        f"report {role} disappeared before sealing: {asset_path}"
                    )
                current = _artifact_identity(
                    asset_path, role=f"{role}:{filename}",
                )
                if (
                    record.get("size_bytes") != current["size_bytes"]
                    or record.get("sha256") != current["sha256"]
                ):
                    raise RuntimeError(
                        f"report {role} changed after assets manifest creation: {asset_path}"
                    )
                candidates.append((f"{role}:{filename}", asset_path))
    for name, value in sorted((plots_status.get("figures") or {}).items()):
        if value:
            candidates.append((f"plot:{name}", output_dir / str(value)))
    if pymol_cover_rel:
        candidates.append(("pymol_cover", output_dir / pymol_cover_rel))

    artifacts: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for role, raw_path in candidates:
        path = raw_path.expanduser().resolve()
        if not path.is_file():
            raise RuntimeError(
                f"report declared artifact is missing ({role}): {path}"
            )
        if path in seen:
            continue
        seen.add(path)
        artifacts.append(_artifact_identity(path, role=role))

    manifest_path = output_dir / "report_artifacts.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "project_root": str(project_root.resolve()),
                "output_dir": str(output_dir.resolve()),
                "artifacts": artifacts,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return manifest_path


def run(input_data: dict, output_dir: Path) -> dict:
    mode = input_data.get("mode", "round")
    top_n = int(input_data.get("top_n", 10))
    # Interface-v2 tool switches (rdkit_2d / mol_3d / molstar_main / plots /
    # pymol_cover). Absent `tools` ⇒ fall back to the legacy enable_* flags so
    # standalone invocations keep working.
    from io_utils import tool_enabled
    enable_pymol_cover = tool_enabled(
        input_data, "pymol_cover", default=bool(input_data.get("enable_pymol_cover", False)))
    enable_plots = tool_enabled(
        input_data, "plots", default=bool(input_data.get("enable_plots", True)))
    enable_mol_main_viewer = tool_enabled(
        input_data, "molstar_main", default=bool(input_data.get("enable_mol_main_viewer", True)))
    enable_2d = tool_enabled(input_data, "rdkit_2d", default=True)
    enable_mol_3d = tool_enabled(input_data, "mol_3d", default=True)
    assets_source_dir = input_data.get("assets_source_dir")
    override_title = input_data.get("title")

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. project_root resolution + aggregation ----
    if mode == "round":
        project_dir = input_data.get("project_dir")
        if not project_dir:
            return {"success": False, "error": "project_dir required in round mode"}
        project_dir_path = Path(project_dir).expanduser().resolve()
        project_root_path = Path(input_data.get("project_root") or project_dir_path.parent).expanduser().resolve()
        report_data = aggregate_round(
            project_dir_path,
            round_number=input_data.get("round_number"),
            top_n=top_n,
        )
        html_filename = "report.html"
        data_filename = "report_data.json"
    elif mode == "final":
        project_root = input_data.get("project_root")
        if not project_root:
            return {"success": False, "error": "project_root required in final mode"}
        project_root_path = Path(project_root).expanduser().resolve()
        report_data = aggregate_final(project_root_path, top_n=top_n)
        html_filename = "final_report.html"
        data_filename = "final_report_data.json"
    else:
        return {"success": False, "error": f"unknown mode: {mode!r} (expected 'round' or 'final')"}

    if override_title and isinstance(report_data.get("header"), dict):
        report_data["header"]["title"] = override_title
    scientific_cohort = {
        field: input_data.get(field)
        for field in SCIENTIFIC_COHORT_REPORT_FIELDS
        if input_data.get(field) not in (None, "")
    }
    input_cohorts = input_data.get("scientific_cohorts")
    scientific_cohorts = [
        dict(item) for item in input_cohorts
        if isinstance(item, dict)
    ] if isinstance(input_cohorts, list) else []
    if scientific_cohorts:
        report_data["scientific_cohorts"] = scientific_cohorts
    if len(scientific_cohorts) == 1 and not scientific_cohort:
        scientific_cohort = dict(scientific_cohorts[0])
    if scientific_cohort and len(scientific_cohorts) <= 1:
        # Keep run-level provenance visible even when no molecule survives to
        # a candidate card (and therefore no candidate row can carry it).
        report_data["scientific_cohort"] = scientific_cohort

    html_output_path = output_dir / html_filename

    # Each report owns an immutable asset snapshot.  Project-global assets made
    # an older round unverifiable as soon as a later round rewrote the shared
    # manifest (or replaced a JavaScript bundle).
    report_asset_root = output_dir.resolve()
    asset_path_warning: str | None = None

    # ---- 2. assets ----
    bundled_assets_dir = SKILL_DIR / "assets_src"
    asset_manifest = ensure_assets(
        report_asset_root,
        source_dir=Path(assets_source_dir).expanduser().resolve() if assets_source_dir else None,
        bundled_asset_dir=bundled_assets_dir if bundled_assets_dir.exists() else None,
    )

    # ---- 3. 2D structures ----
    if enable_2d:
        structure_2d_status = mol_renderer_2d.render_candidates(report_data.get("candidates") or [])
    else:
        structure_2d_status = {"available": False, "reason": "disabled_by_config"}

    # ---- 4. 3D per-molecule cards + main Mol* viewer ----
    step0 = report_data.get("step0") or {}
    step1 = report_data.get("step1") or {}
    step3 = report_data.get("step3") or {}
    step4 = report_data.get("step4") or {}
    receptor_pdb = step0.get("prepared_pdb")
    reference_ligand_sdf = step0.get("reference_ligand_sdf")
    selected_poses_sdf = step4.get("selected_poses_sdf")
    docked_poses_sdf = step3.get("docked_poses_sdf")
    pocket_residues = step1.get("consensus_pocket_residues") or []

    cards_status = {"rendered": 0, "failed": 0, "total": 0}
    if not enable_mol_3d:
        for c in report_data.get("candidates") or []:
            c["viewer_3d_html"] = (
                '<div class="viewer-3dmol viewer-placeholder">'
                '<div class="viewer-error">3D cards disabled</div></div>'
            )
            c["viewer_3d_error"] = "mol_3d disabled_by_config"
            c["viewer_3d_has_pose"] = False
        cards_status = {"rendered": 0, "failed": 0, "total": 0, "reason": "disabled_by_config"}
    elif asset_available(asset_manifest, "3Dmol-min.js"):
        cards_status = mol_renderer_3d.render_candidate_cards(
            candidates=report_data.get("candidates") or [],
            receptor_pdb_path=receptor_pdb,
            selected_sdf_path=selected_poses_sdf,
            docked_sdf_path=docked_poses_sdf,
            pocket_residues=pocket_residues,
        )
    else:
        # Attach placeholder HTML per candidate so template can still render.
        for c in report_data.get("candidates") or []:
            c["viewer_3d_html"] = (
                '<div class="viewer-3dmol viewer-placeholder">'
                '<div class="viewer-error">3Dmol.js asset missing</div></div>'
            )
            c["viewer_3d_error"] = "3Dmol-min.js not in assets/"
            c["viewer_3d_has_pose"] = False

    main_viewer_html: str | None = None
    main_viewer_error: str | None = None
    if enable_mol_main_viewer:
        if asset_available(asset_manifest, "molstar.js"):
            main_viewer_html, main_viewer_error = mol_renderer_3d.render_molstar_main(
                receptor_pdb_path=receptor_pdb,
                reference_ligand_sdf_path=reference_ligand_sdf,
                html_output_path=html_output_path,
            )
        else:
            main_viewer_error = "molstar.js asset missing"

    # ---- 5. optional PyMOL static cover ----
    pymol_cover_rel: str | None = None
    pymol_cover_error: str | None = None
    if enable_pymol_cover:
        pymol_cover_rel, pymol_cover_error = pymol_renderer.render_cover(
            receptor_pdb_path=receptor_pdb,
            reference_ligand_sdf_path=reference_ligand_sdf,
            pocket_residues=pocket_residues,
            output_dir=output_dir,
        )

    # ---- 6. optional matplotlib plots ----
    plots_status: dict[str, Any]
    if enable_plots:
        plots_status = plot_renderer.render_all(report_data, output_dir=output_dir)
    else:
        plots_status = {"available": False, "reason": "plots disabled by input"}

    # ---- 7. write HTML + report_data.json ----
    html_text = html_template.render_html(
        report_data=report_data,
        asset_manifest=asset_manifest,
        plots_status=plots_status,
        pymol_cover_rel=pymol_cover_rel,
        main_viewer_html=main_viewer_html,
        html_output_path=html_output_path,
        project_root=report_asset_root,
    )
    html_output_path.write_text(html_text, encoding="utf-8")

    # Lookup indexes share records in memory, but JSON would serialize every
    # alias again. Keep canonical evidence once and stream the machine sidecar.
    report_data_dump = _slim_for_json(report_data)
    with (output_dir / data_filename).open("w", encoding="utf-8") as stream:
        json.dump(
            report_data_dump, stream, ensure_ascii=False,
            separators=(",", ":"), default=str,
        )

    report_data_path = output_dir / data_filename
    assets_manifest_path = report_asset_root / "assets" / "assets_manifest.json"
    primary_files = {
        "report_html": html_output_path,
        "report_data": report_data_path,
    }
    if assets_manifest_path.is_file():
        primary_files["assets_manifest"] = assets_manifest_path
    report_artifact_manifest = _write_report_artifact_manifest(
        output_dir=output_dir,
        project_root=project_root_path,
        primary_files=primary_files,
        asset_manifest=asset_manifest,
        plots_status=plots_status,
        pymol_cover_rel=pymol_cover_rel,
    )

    files = {
        "report_html": str(html_output_path),
        "report_data_json": str(report_data_path),
        "report_artifact_manifest_json": str(report_artifact_manifest),
    }
    if assets_manifest_path.is_file():
        files["assets_manifest_json"] = str(assets_manifest_path)
    legacy_data_key = f"{Path(data_filename).stem}_json"
    files.setdefault(legacy_data_key, str(report_data_path))
    if pymol_cover_rel:
        files["pymol_cover_png"] = str(output_dir / pymol_cover_rel)

    result = {
        "success": True,
        "mode": mode,
        "project_root": str(project_root_path),
        "html_output": str(html_output_path),
        "assets": asset_manifest,
        "structure_2d_status": structure_2d_status,
        "cards_status": cards_status,
        "main_viewer_error": main_viewer_error,
        "pymol_cover_error": pymol_cover_error,
        "plots_status": plots_status,
        "candidate_count": len(report_data.get("candidates") or []),
        "missing_sources": report_data.get("missing_sources") or [],
        "asset_path_warning": asset_path_warning,
        "files": files,
    }
    if scientific_cohorts:
        result["scientific_cohorts"] = scientific_cohorts
    if scientific_cohort and len(scientific_cohorts) <= 1:
        result["scientific_cohort"] = scientific_cohort
        result.update(scientific_cohort)
    return result


def _slim_for_json(report_data: dict) -> dict:
    """Omit redundant lookup views without mutating rendering data or evidence."""
    slim = dict(report_data)
    storage = report_data.get("report_storage") or {}
    omitted: dict[str, str] = dict(storage.get("omitted_derived_fields") or {})
    step4 = report_data.get("step4")
    if isinstance(step4, dict):
        slim["step4"] = dict(step4)
        evaluation = step4.get("evaluation_report")
        if isinstance(evaluation, dict):
            # Match the reader's current and legacy per-pose array spellings.
            fields = [
                key for key in ("evaluation", "molecules", "results")
                if isinstance(evaluation.get(key), list)
            ]
            field = next((key for key in fields if evaluation[key]), None)
            field = field or (fields[0] if fields else None)
            if field is not None:
                rows = evaluation[field]
                per_molecule = {}
                per_molecule_entries: dict[str, list] = {}
                for row in rows:
                    name = row.get("mol_name") or row.get("name")
                    if name:
                        name = str(name)
                        per_molecule.setdefault(name, row)
                        per_molecule_entries.setdefault(name, []).append(row)
                for key, derived in (
                    ("molecules", rows),
                    ("per_molecule", per_molecule),
                    ("per_molecule_entries", per_molecule_entries),
                ):
                    # Retain nonstandard/legacy views if they contain evidence
                    # absent from the canonical report; never discard by name.
                    if key in step4 and step4[key] == derived:
                        slim["step4"].pop(key)
                        omitted[f"step4.{key}"] = f"step4.evaluation_report.{field}"
    step5 = report_data.get("step5")
    if isinstance(step5, dict):
        slim["step5"] = dict(step5)
        properties = step5.get("property_filter_report")
        if (
            isinstance(properties, dict)
            and isinstance(properties.get("results"), list)
            and "results" in step5
            and step5["results"] == properties["results"]
        ):
            slim["step5"].pop("results")
            omitted["step5.results"] = "step5.property_filter_report.results"
    slim["report_storage"] = {
        "schema_version": 1,
        "omitted_derived_fields": omitted,
    }
    slim_candidates = []
    for c in report_data.get("candidates") or []:
        c2 = dict(c)
        c2.pop("viewer_3d_html", None)
        # keep an indicator whether a 2D image was produced without embedding the base64 blob
        if c2.get("structure_2d_png"):
            c2["structure_2d_png"] = "(embedded in report.html)"
        slim_candidates.append(c2)
    slim["candidates"] = slim_candidates
    if slim.get("mode") == "final":
        slim_rounds = []
        for r in slim.get("rounds") or []:
            slim_rounds.append(_slim_for_json(r))
        slim["rounds"] = slim_rounds
    return slim


def main():
    from io_utils import parse_skill_args, print_skill_result
    input_data, output_dir = parse_skill_args("SBDD Report Generator")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result.get("success") else 1)


if __name__ == "__main__":
    main()

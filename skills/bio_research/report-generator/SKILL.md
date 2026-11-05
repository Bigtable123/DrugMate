---
name: report-generator
description: Assemble a self-contained HTML report for one SBDD round (or the whole project) from Step 0-6 artifacts. Renders per-molecule cards with 3D pose viewers (3Dmol.js), 2D structures (RDKit), and metric badges; a top Mol* main viewer for the pocket + reference ligand; optional QC/druglikeness charts (matplotlib) and ROC/EF panels when ACT_/DEC_ labels are present. All heavy renderers (RDKit, matplotlib, PyMOL) graceful-skip when unavailable; the aggregation + HTML layer are pure stdlib and always run.
metadata:
  skill-author: SBDD Agent
  version: "1.0"
  dependencies: []
  optional-dependencies: [rdkit, matplotlib, jinja2, pymol]
  external-tools: []
  python: ">=3.10"
---

# Report Generator

Turn Step 0-6 (and optionally cross-round) artifacts into a self-contained HTML report intended for a human reviewer or a downstream frontend to iframe. This is Step 7 of the SBDD pipeline; it does **not** feed the Agent's iteration decision — the Agent still reads `coordinator_report.json`.

## When to Use

- **Per-round mode**: called automatically as Step 7 by [pipeline-entry/run_pipeline.py](../pipeline-entry/run_pipeline.py) after Step 6. Renders `round{N}/step7/report.html` from that round's artifacts.
- **Final mode**: called once by the Agent after all rounds finish, before writing `final_report.md`. Renders `{project_dir}/final_report.html` aggregating every round's `round_history.json` / `candidate_ranking.json`.

The report does not affect Agent iteration; it can be regenerated at any time from the same artifacts.

## Input Schema

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "mode": "round",
  "project_dir": "<SBDD_DATA_ROOT>/20260626-120000-cdk2-4kd1/round1",
  "project_root": "<SBDD_DATA_ROOT>/20260626-120000-cdk2-4kd1",
  "round_number": 1,
  "top_n": 10,
  "include_reference_ligand": true,
  "enable_pymol_cover": false,
  "enable_plots": true,
  "enable_mol_main_viewer": true,
  "assets_source_dir": null,
  "title": null
}
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `mode` | string | Yes | - | `"round"` (single round) or `"final"` (cross-round) |
| `project_dir` | string | Required when `mode=round` | - | Path to `round{N}/` directory holding step0..step6 |
| `project_root` | string | Required when `mode=final` | - | Path to `{data_root}/{timestamp}-{slug}/`; report scans all `round*/` subdirectories |
| `round_number` | int | No | Parsed from `project_dir` | Overrides the round label shown in the header |
| `top_n` | int | No | 10 | Number of non-decoy candidates to render per scientific cohort/pool (source: `candidate_ranking.json`) |
| `include_reference_ligand` | bool | No | true | Overlay reference ligand from Step 0 in the main Mol* viewer |
| `enable_pymol_cover` | bool | No | false | Try to render a static PyMOL cover PNG (graceful-skip when PyMOL unavailable) |
| `enable_plots` | bool | No | true | Render QC/druglikeness/ROC-EF charts with matplotlib (graceful-skip) |
| `enable_mol_main_viewer` | bool | No | true | Include the Mol* top viewer (pocket + reference ligand + optional top pose overlay) |
| `assets_source_dir` | string | No | Bundled in `report-generator/assets_src/` if present, otherwise use CDN placeholder text | Directory containing `3Dmol-min.js`, `molstar.js`, `molstar.css` to copy into `{project_root}/assets/` |
| `title` | string | No | Auto-generated from target/PDB metadata | Overrides HTML `<title>` |

Files are referenced from `{project_root}/assets/` via relative paths, so a single asset copy serves every round's report.

## Data Sources

The aggregator (pure stdlib) reads only files that already exist under `step0..step6`; nothing is recomputed.

**Per-round mode reads:**

| Source | Consumed for |
|---|---|
| `step0/*_receptor_prepared.pdb`, `step0/*_ligand_*.sdf` | Receptor + reference ligand for 3D viewers |
| `step1/pocket_analysis_report.json`, `interaction_blueprint.json`, `pocket_features.json` | Pocket header + hard-gate badges + blueprint required list |
| `step3/docking_report.json`, `docking_scores.json`, `docked_poses.sdf` | Redocking evidence + native per-pose score records (GNINA/Uni-Dock/Vina) |
| `step4/evaluation_report.json`, `selected_poses.sdf`, `selection_manifest.csv` | Per-molecule verdict, PoseBusters/PoseCheck/PLIP diagnostics, selected pose blocks |
| `step5/property_filter_report.json`, `synthesis_report.json` | Identity-aware QED/SA/Lipinski/alerts/ADMET plus top-K retrosynthesis status per candidate; executable public property metrics are also consumed by Step 6 gating/ranking |
| `step6/coordinator_report.json`, `candidate_ranking.json`, `round_summary.json` | Cohort-partitioned Top-N non-decoy lists (DEC_ excluded here — see reflective-coordinator rules), local OverallRank fields, header badges, ROC/EF panel |

**Final mode additionally reads:**

- `{project_root}/round*/step6/round_summary.json` and `round_history.json` → cross-round trend table.
- For `round-assessment/v1`, candidate source resolves to the latest successfully
  assessed round with `round_assessment.valid_for_nomination=true` and status
  `valid` or `review_required`. Eligible keep candidates remain visible while
  the Agent continues the campaign; campaign completion is not a requirement.
  Invalid, not-evaluated, or nomination-blocked rounds remain diagnostic context.
  Final cards still require an eligible keep record, so modify/reject/blocked
  records cannot displace a nomination card.
- Legacy reports retain their existing convergence contract: structured
  `decision_outcome=success` plus `valid_for_nomination=true`, or legacy fallback
  excluding `terminate_invalid`.
- New headers and cross-round trends display round assessment status. They do
  not fabricate campaign success, convergence, or a stop recommendation.

## Output Files

**Per-round mode** (`--output-dir` typically `{project_dir}/step7/`):

| File | Description |
|------|-------------|
| `report.html` | Self-contained HTML for that round |
| `report_data.json` | Compact JSON retaining canonical evidence, candidate ranking contracts and cohort lineage; redundant rendering indexes are omitted (see [storage contract](references/layout.md#json-storage-contract)) |
| `figures/*.png` | matplotlib PNG figures (only when matplotlib available) |
| `structures/mol_*.html` | Per-molecule 3Dmol.js viewer HTML fragments referenced by the report |

**Final mode** (`--output-dir` typically `{project_root}/`):

| File | Description |
|------|-------------|
| `final_report.html` | Cross-round HTML report |
| `final_report_data.json` | Compact cross-round JSON with the same storage and cohort-local ranking contracts |
| `figures/`, `structures/` | Same layout as per-round mode |

**Shared under `{project_root}/assets/`:**

| File | Description |
|------|-------------|
| `3Dmol-min.js` | Cards' 3D viewer (small — kept inline via `<script src="../../assets/3Dmol-min.js">`) |
| `molstar.js`, `molstar.css` | Main viewer bundle |
| `assets_manifest.json` | Which asset came from where + checksum + timestamp |

If `assets_source_dir` is not provided and no bundled copy exists, `assets_manifest.json` records the missing files and the HTML falls back to inline placeholder text (report still renders, 3D viewers show an "asset missing" message with instructions).

## HTML Layout

```
┌──────────────────────────────────────────────────────────────┐
│ Header (lightweight overview)                                │
│  • Target / PDB / Round N / Decision badge                   │
│  • Quality gate chips: redocking, pocket, blueprint, AUROC   │
│  • Cross-round mini-table (final mode only)                  │
├──────────────────────────────────────────────────────────────┤
│ Main viewer (Mol*, optional)                                 │
│  Pocket residues + reference ligand + toggleable Top-N poses │
├──────────────────────────────────────────────────────────────┤
│ AiZynthFinder execution audit                                │
│  attempted / not attempted / error counts + status meanings  │
├──────────────────────────────────────────────────────────────┤
│ Candidate cards (Top-N from candidate_ranking.json)          │
│  ┌──────────────────────────────────────────────────────┐    │
│  │ mol_name  verdict badge  flags                       │    │
│  │ ┌─── 3Dmol.js pose ───┐  ┌─ RDKit 2D ─┐  ┌ metrics ┐ │    │
│  │ │ receptor + pose     │  │ SMILES pic │  │ scores  │ │    │
│  │ │ (blueprint hits)    │  │            │  │ QED/SA  │ │    │
│  │ └─────────────────────┘  └────────────┘  │ Lipinski│ │    │
│  │ diagnosis, modification_hint             │ PAINS   │ │    │
│  │                                          │ ADMET   │ │    │
│  │                                          └─────────┘ │    │
│  └──────────────────────────────────────────────────────┘    │
├──────────────────────────────────────────────────────────────┤
│ QC & druglikeness panels (optional, matplotlib)              │
│  redocking RMSD · pass/warn/fail histogram · QED-SA scatter  │
├──────────────────────────────────────────────────────────────┤
│ Retrospective screening (optional)                           │
│  Rendered only when ACT_/DEC_ counts pass min thresholds:    │
│  AUROC curve + bootstrap CI + EF1/5/10 bars                  │
├──────────────────────────────────────────────────────────────┤
│ Footer: skips log, asset manifest, generator version         │
└──────────────────────────────────────────────────────────────┘
```

## Candidate Source Rules

- Top-N candidates are selected independently inside each `cohort_id +
  normalization_pool` partition in `step6/candidate_ranking.json`. Each card
  preserves the coordinator-provided float/null local OverallRank instead of
  being globally renumbered; ranks from different cohorts are never compared.
  The report propagates the cohort components and labels every card with its
  cohort. That file already excludes `DEC_*` and `reject` rows per the
  coordinator's rules. Final mode displays only eligible keep rows; round mode
  may show the other gate tiers it contains (`modify`, protocol-blocked keep)
  for diagnosis. A rejected molecule is never a candidate card — read it from
  `coordinator_report.json.molecules` or the published root tables instead.
- In mixed transitional artifacts, schema-2 `overall_score`/`overall_rank`
  fields are authoritative even when null. Legacy flat `weighted_*` rows are
  retained only as unranked audit rows and cannot enter the schema-2 ordering.
- Never source candidates from `docked_poses.sdf` or `evaluation_report.json.molecules` — those include decoys.
- In both round and final mode, also show separate **KEEP Top 1–5** and
  **MODIFY Top 1–5** review tables for each cohort/pool. These compact tables
  come from all scored rows in `candidate_ranking.json`, independently of the
  card `top_n` limit. Sort by the existing `overall_score`, never recalculate
  scores within a class, and exclude decoys and rows without a finite score.
  Equal scores use molecule name then pose ID as deterministic tie-breaks.
  Keep the original nullable `overall_rank`; label 1–5 as `review_position`,
  not a released OverallRank. `report_data.json.candidate_reviews` carries
  `review_only=true`, `source_round`, and the same grouped rows.
- Review tables include blocked rows with their original status and separate
  molecular category findings from missing evidence and protocol blockers.
  They do not nominate MODIFY or blocked KEEP records. Final review tables use
  the nomination source round when one exists; otherwise they use the latest
  readable, successfully assessed round, with its source round stated.
- Each rendered candidate binds to the exact selected pose: `mol_name`, `pose_index`, `source_pose_index`, `selection_rank`.

## Final Markdown Handoff

After rendering `final_report.html`, the Agent writes `final_report.md`. The
Markdown must also include the grouped KEEP/MODIFY review tables from
`final_report_data.json.candidate_reviews`, preserving the source round, review
label, original scores/ranks, and nomination blockers. The synthesis section
must use `final_report_data.json.synthesis_summary`
and each candidate's normalized `synthesis` object; it must explicitly print
`attempted`, `status`, and `reason`. Never infer execution from a blank metric
cell or from the overall keep/modify/reject totals.

- `not_selected_top_k`, `not_eligible`, and `disabled` mean
  `attempted=false`; their synthesis metrics are N/A, not failed.
- `error` means `attempted=true` but the tool call failed. It is an operational
  failure requiring protocol review, not an automatic molecular reject.
- `not_found_within_budget` means `attempted=true` and the search completed
  without a solved route in the configured budget. It is a search outcome, not
  a tool failure, and must not be rewritten as `reject`.
- Keep AiZynthFinder execution counts separate from molecular verdict counts.
  A statement such as "998 rejected" must never be presented as "998
  AiZynthFinder failures" unless 998 rows actually have
  `attempted=true,status=error`.

## Graceful-Skip Contract

Every renderer records into `steps_skipped[]` with `{"step": ..., "reason": ...}` when it can't run:

| Renderer | Skip trigger | HTML fallback |
|---|---|---|
| RDKit 2D structure | `rdkit` unimportable | Card shows SMILES as monospace text |
| matplotlib plots | `matplotlib` unimportable | Section replaced with a "plots disabled" note |
| PyMOL cover | `pymol` unimportable OR `enable_pymol_cover=false` | Section hidden |
| Mol* main viewer | Asset missing / `enable_mol_main_viewer=false` | Section hidden, top-viewer container reads "main viewer disabled" |
| 3Dmol.js card | Asset missing OR no pose SDF | Card shows PDB path + pose text with a note |
| Jinja2 template | `jinja2` unimportable | Fall back to f-string template (identical layout, less flexible) |

The **aggregator and HTML writer are pure stdlib** — `report.html` and `report_data.json` are always produced.

## Usage

```bash
# Per-round (called automatically by pipeline-entry as Step 7)
python report-generator/run.py \
  --input "$(python3 -c 'import json, os; print(json.dumps({"mode": "round", "project_dir": os.path.join(os.environ["SBDD_DATA_ROOT"], "20260626-120000-cdk2-4kd1/round1"), "top_n": 10}))')" \
  --output-dir "${SBDD_DATA_ROOT}/20260626-120000-cdk2-4kd1/round1/step7/"

# Final (Agent invokes after all rounds finish)
python report-generator/run.py \
  --input "$(python3 -c 'import json, os; print(json.dumps({"mode": "final", "project_root": os.path.join(os.environ["SBDD_DATA_ROOT"], "20260626-120000-cdk2-4kd1"), "top_n": 10}))')" \
  --output-dir "${SBDD_DATA_ROOT}/20260626-120000-cdk2-4kd1/"
```

## Integration with pipeline-entry

`STEP_SKILLS` is extended:

```python
STEP_SKILLS = {
    0: "structure-prep",
    1: "pocket-analysis",
    2: "molecule-generator",
    3: "docking",
    4: "structure-evaluator",
    5: "drug-property-filter",
    6: "reflective-coordinator",
    7: "report-generator",
}
```

`_build_step_input` for step 7 constructs:

```json
{
  "mode": "round",
  "project_dir": "<project_dir>",
  "project_root": "<project_dir.parent>",
  "round_number": "<round_number>",
  "top_n": "<config.report_top_n or 10>",
  "enable_pymol_cover": "<config.report_enable_pymol or false>",
  "enable_plots": "<config.report_enable_plots or true>",
  "enable_mol_main_viewer": "<config.report_enable_mol or true>",
  "assets_source_dir": "<config.report_assets_source_dir>"
}
```

`stop_step` default rises from 6 to 7. Existing runs that pass `stop_step: 6` are unaffected. Report generation failure never blocks pipeline success flags for step0-6 — Step 7 is treated as best-effort (failure logged but overall `success: true` if step6 passed and Step 7 wrote at least `report_data.json`).

## References

- `references/layout.md` — HTML section-by-section field-source table
- `references/assets.md` — how to obtain 3Dmol.js and Mol* offline bundles and place them under `assets_source_dir`

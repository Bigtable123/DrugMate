# Report Layout — section-by-section source mapping

Every value the HTML shows comes from Step 0-6 artifacts already produced by
[pipeline-entry](../../pipeline-entry/). Nothing is recomputed. All aggregation
lives in [scripts/aggregator.py](../scripts/aggregator.py).

## Header (lightweight overview)

| Field | Source key | File |
|---|---|---|
| `title` | `metadata.target_name` + `metadata.pdb_id` | `step0/result.json` |
| `target_name` | `metadata.target_name` | `step0/result.json` |
| `pdb_id` | `metadata.pdb_id` | `step0/result.json` |
| `chembl_target_id` | `metadata.chembl_target_id` | `step0/result.json` |
| `round_label` | derived from `project_dir` name (e.g. `round2`) | — |
| `assessment_status` (legacy: `decision`) | `round_assessment.status` (legacy: `decision`) | `step6/coordinator_report.json` |
| `quality_gate_passed` | `quality_gate_passed` | `step6/coordinator_report.json` |
| `primary_bottleneck` | `diagnostics.diagnosis.bottleneck_type` (legacy: `next_round_instructions`) | `step6/coordinator_report.json` |
| `num_screened` / `num_keep` / `num_modify` / `num_reject` | matching keys | `step6/round_summary.json` |

Header always renders — missing values are shown as `—`.

## Quality gates

| Gate | Passed | Value / Threshold | File |
|---|---|---|---|
| Redocking RMSD | `redocking_validation.passed` | `.rmsd` / `.threshold` | `step3/docking_report.json` |
| Pocket validation | `validation.passed` | `.recall` / `.recall_threshold` | `step1/pocket_analysis_report.json` |
| Blueprint hard-gating | `valid_for_hard_gating` | `len(required)` / 1 | `step1/interaction_blueprint.json` |
| Screening AUROC | `screening_benchmark.passed` | `.auroc` / `.min_auroc` | `step6/coordinator_report.json` |
| Round assessment | `quality_gate_passed` | `round_assessment.status` | `step6/coordinator_report.json` |

## Cross-round trend (final mode only)

`aggregate_final()` walks every `round*/step6/round_summary.json` and
`coordinator_report.json` to build one row per round. `candidate_source_round`
is the last round whose Step 6 executed successfully and whose coordinator
report does not explicitly block nomination. With the structured contract this
means `round_assessment.valid_for_nomination == true` and assessment status
`valid` or `review_required`, independent of campaign completion. Invalid and
not-evaluated rounds supply diagnostic context without candidate cards. Legacy
reports retain the old outcome/nomination contract.

## Mol* main viewer

Loaded via relative URL from the HTML file to:

- `receptor_pdb_path` = `step0/*_receptor_prepared.pdb`
- `reference_ligand_sdf_path` = `step0/*_ligand_*.sdf`

Molstar bundle: `../../assets/molstar.js` (per-round HTML at
`round{N}/step7/report.html`) or `assets/molstar.js` (final report at
`{project_root}/final_report.html`).

Skipped when `enable_mol_main_viewer=false` or asset absent — placeholder text
shown instead of the viewer.

## Per-candidate cards

Source of truth: `step6/candidate_ranking.json` (DEC_* already excluded by
[reflective-coordinator](../../reflective-coordinator/SKILL.md)).
Cards are grouped by `cohort_id + normalization_pool`, sorted by the local
float/null `overall_rank`, and truncated to `top_n` independently per group.
Ranks from different cohorts are never compared or globally renumbered.

For each candidate:

| Card field | Source |
|---|---|
| `mol_name`, `pose_index`, `smiles`, local `rank`/`overall_rank` | `candidate_ranking.json` (authoritative for every numeric/null field) |
| `cohort_id`, cohort components, normalization pool, rank scope | `candidate_ranking.json`, with CSV fallbacks |
| category/step/overall scores, coverage, rank status | `candidate_ranking.json`; schema-2 nulls never fall back to legacy `weighted_*` values |
| verdict, flags, diagnosis, modification_hint | `step4/evaluation_report.json.molecules[].{verdict,flags,diagnosis,modification_hint}` |
| primary score metric/value/unit/direction | `step6/candidate_ranking.json` + pose-bound `pose_scores`; legacy GNINA fields remain visible when present |
| blueprint_match | `step4/evaluation_report.json.molecules[].blueprint_match` |
| posebusters | `step4/evaluation_report.json.molecules[].posebusters` |
| posecheck | `step4/evaluation_report.json.molecules[].posecheck` |
| pocket_contact_fraction | `step4/evaluation_report.json.molecules[].pocket_contact_fraction` |
| QED / SA / MW / LogP | `step5/property_filter_report.json.results[].properties` |
| lipinski | `step5/property_filter_report.json.results[].lipinski` |
| pains | `step5/property_filter_report.json.results[].pains` |
| admet | `step5/property_filter_report.json.results[].admet` |
| AiZynthFinder attempted/status/reason | normalized from `step5/property_filter_report.json.results[].synthesis`; legacy statuses infer only the defined attempt bit and are labelled as inferred |

## KEEP and MODIFY Top 1–5 review tables

`candidate_reviews` is a compact view of the full Step 6 ranking artifact,
independent of candidate-card `top_n`. Each cohort/pool has separate KEEP and
MODIFY groups sorted by existing OverallScore. No score is recomputed and a
missing OverallRank stays null. `review_position` is only a display position;
rows without a finite score are counted but do not enter score Top 5.

Show `source_round` and `review_only=true` in final and round reports. Category
verdicts and observed flag/fail metric IDs are molecular findings; missing
metric IDs and protocol blockers are separate columns. All nomination status
values remain unchanged, including on blocked rows. The final report may use
its latest readable, successfully assessed round for review when there is no
nomination-valid source; its `candidates` list remains empty in that case.

## AiZynthFinder execution audit

The report summarizes every Step 5 row independently of the coordinator's
molecule verdict totals. It prints the attempted, not-attempted, unknown, and
operational-error counts plus one row per status. `not_selected_top_k`,
`not_eligible`, and `disabled` are not executions; `error` is an operational
failure with N/A metrics; and `not_found_within_budget` is a completed search
outcome rather than a tool failure. The same normalized fields are retained in
`report_data.json` for the Agent's final Markdown handoff.

Per-card 3D view (3Dmol.js) inlines the receptor PDB text + the pose block extracted from
`step4/selected_poses.sdf` (or `step3/docked_poses.sdf` if selected is missing) as JSON strings.
Pose extraction matches by `_Name` first, then by index.

## QC / druglikeness plots (optional)

`plot_renderer.render_all()` produces PNGs into `figures/`:

- `verdict_bar.png` ← `step4/evaluation_report.json.summary`
- `blueprint_hit_hist.png` ← `step4/evaluation_report.json.molecules[].blueprint_match.required_hit_rate` grouped by verdict
- `score_scatter.png` ← `step3/docking_scores.json` × `step4` verdict, DEC_ desaturated
- `qed_sa_scatter.png` ← `step5/property_filter_report.json.results[].{properties.QED, properties.SA_score, verdict}`
- `redocking_gauge.png` ← `step3/docking_report.json.redocking_validation`
- `trend_pass_rate.png` (final mode) ← `trend` table

## Retrospective screening (optional, ACT_/DEC_ only)

Rendered when `step6/coordinator_report.json.screening_validation.available` is
true. The legacy `screening_benchmark` and `retrospective_screening` block names
remain supported.

- AUROC and bootstrap 95% CI from `AUROC` / `AUROC_bootstrap_95ci`; lowercase
  legacy aliases remain accepted.
- EF1/5/10 from `EF1%`, `EF5%`, and `EF10%`, retaining older aliases.
- Average precision from unrounded `average_precision`, displayed as a percentage.
  The rounded `PR_AUC` alias is used only when `PR_AUC_method` explicitly says
  `tie_aware_average_precision`; arbitrary PR-AUC values are not relabeled as AP.
- Top-20 evaluated count, active hits and hit ratio from `top20_size`,
  `top20_actives` and `top20_hit_ratio`. With fewer than 20 evaluable labeled
  compounds, the denominator is the actual count. Observed equal-score order is
  stable input order.
- If `top20_boundary_tie` is true, show the expected hit count and ratio under
  uniform ordering within the cutoff tie, plus the feasible hit-count range.
  These are uncertainty summaries for this report; they do not change the
  observed ranking order. Exact AP and all Top-20 tie metadata remain in
  the aggregated JSON without presentation rounding.
- ROC curve from `roc_curve.{fpr,tpr}` when the coordinator recorded it.

## Footer

- Generated-at timestamp
- Asset presence line (per-file)
- Plot availability line + errors
- Missing sources list — every step directory that was absent or unreadable

## JSON storage contract

HTML rendering uses the complete in-memory aggregate. The JSON sidecar keeps
canonical scientific evidence, scores, ranks, candidate reviews, and cohort/pose
identities unchanged, but omits duplicate derived lookup views:

| Omitted view | Canonical records retained in JSON |
|---|---|
| `step4.molecules` | `step4.evaluation_report.evaluation` |
| `step4.per_molecule` | Same records; first row per `mol_name` (legacy fallback: `name`) |
| `step4.per_molecule_entries` | Same records; all rows grouped by molecule name |
| `step5.results` | `step5.property_filter_report.results` |

Legacy evaluation arrays named `molecules` or `results` are preserved under their
original key. A view is removed only if it matches the canonical records; unique
legacy evidence is retained. Do not use molecule names to merge scientific
records across poses or cohorts.

`report_storage.schema_version=1` identifies this storage contract.
`report_storage.omitted_derived_fields` maps each removed view to the actual
canonical array path within the same JSON. Older sidecars have no storage
metadata and may contain both views; consumers should prefer the canonical
arrays, falling back to the legacy views when no canonical array exists.

The rule also applies to each entry in final-mode `rounds`. JSON is written
incrementally without indentation; embedded image/viewer payloads are still
excluded. No external evidence files are required to replace the omitted views,
and the input aggregate and original experiment artifacts are not modified.

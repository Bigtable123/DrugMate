---
name: reflective-coordinator
description: Assess and rank SBDD evidence for research campaigns.
metadata:
  skill-author: SBDD Agent
  version: "1.0"
  dependencies: []
  optional-dependencies: []
  external-tools: []
  python: ">=3.10"
---

# Reflective Coordinator

Aggregate upstream results into a deterministic round assessment, rankings, and diagnostic evidence. The Agent owns R1–R3 campaign planning and the choice to continue, pause, or stop.

## When to Use

- After structure evaluation (Step 4) and drug property filter (Step 5) to assess round evidence
- When the pipeline needs to iterate: modify candidates and re-run docking/evaluation
- To systematically identify why candidates are failing and what to change

## Input Schema

```json
{
  "evaluation_report_file": "path/to/evaluation_report.json",
  "property_filter_report_file": "path/to/property_filter_report.json",
  "docking_scores_file": "path/to/docking_scores.json",
  "docking_report_file": "path/to/docking_report.json",
  "pocket_analysis_report_file": "path/to/pocket_analysis_report.json",
  "require_quality_gates": true,
  "min_validation_actives": 20,
  "min_validation_decoys": 50,
  "min_validation_auroc": 0.70,
  "min_validation_auroc_hard": 0.50,
  "min_validation_ef1": 3.0,
  "min_validation_ef1_hard": 1.0,
  "round_number": 1,
  "normalization_method": "rank_avg_percentile_v1",
  "screening_mode": "V2",
  "candidate_source": "sampled_library",
  "metric_weight_overrides": {},
  "category_weight_overrides": {}
}
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `evaluation_report_file` | string | Yes | - | Path to structure evaluator report |
| `property_filter_report_file` | string | No | - | Path to drug property filter report |
| `docking_scores_file` | string | No | - | Pose-bound native/structured scores used for exact-pose selection and ranking |
| `docking_report_file` | string | Recommended | - | Supplies the native-redocking quality gate |
| `pocket_analysis_report_file` | string | Recommended | - | Supplies pocket and blueprint quality gates |
| `require_quality_gates` | bool | No | true | Enforce protocol evidence validity and nomination blocking; reject marks the round invalid |
| `min_validation_actives` | int | No | 20 | Minimum labeled actives for a retrospective screening benchmark |
| `min_validation_decoys` | int | No | 50 | Minimum labeled decoys for a retrospective screening benchmark |
| `min_validation_auroc` | float | No | 0.70 | Soft keep boundary for compound-level AUROC; the hard reject boundary is 0.50 |
| `min_validation_auroc_hard` | float | No | 0.50 | Hard reject boundary for compound-level AUROC |
| `min_validation_ef1` | float | No | 3.0 | Soft keep boundary for retrospective EF1% |
| `min_validation_ef1_hard` | float | No | 1.0 | Hard reject boundary for retrospective EF1% |
| `round_number` | int | No | 1 | Current optimization round |
| `max_rounds` | int | No | 3 | Deprecated compatibility input, recorded but ignored; the Agent owns campaign budget |
| `pass_rate_threshold` | float | No | 0.3 | Deprecated compatibility input, recorded but ignored |
| `min_candidates` | int | No | 3 | Deprecated compatibility input, recorded but ignored |
| `normalization_method` | string | No | `rank_avg_percentile_v1` | Descending RANK.AVG percentile at metric and category levels inside one scientific cohort |
| `screening_mode` | string | No | - | `V1` marking / `V2` standard / `V3` strict level-2 fold. Unset means it is derived from `candidate_source`, falling back to `V2` |
| `candidate_source` | string | No | - | Provenance used to derive `screening_mode`: `molecule_generation` → `V3`, `sampled_library` → `V2`, `user_supplied_full` → `V1` |
| `metric_weight_overrides` | object | No | `{}` | Session-local `metric_id → 0..5` magnitude; preserves catalog sign/direction and every hard gate |
| `category_weight_overrides` | object | No | `{}` | Session-local category-code parent weights from `0..5`; hard gates remain unchanged |

## Round Assessment Contract

`assessment_schema_version: "round-assessment/v1"` identifies the new contract.
`round_assessment.status` describes evidence validity:

| Status | Meaning |
|---|---|
| `valid` | Evaluated candidate evidence and required protocol gates are complete; keep, modify, and reject may coexist. |
| `review_required` | Required protocol evidence or non-rejected candidate evidence needs repair/review. |
| `invalid` | A required protocol gate rejected the experiment; preserve evidence and publish non-destructive rollback instructions. |
| `not_evaluated` | No nomination candidate has an evaluated molecular classification. |

`valid_for_nomination` is true when at least one trusted eligible keep exists
and the required protocol gates pass. There is no 3-candidate or 30% success
threshold. A partial candidate evidence gap may make the round
`review_required` while other candidates remain eligible. All-reject rounds
can have valid evidence and no eligible candidates; that never terminates the
campaign automatically.

`diagnostics` contains advisory bottlenecks, affected molecules, and chemical
modification hints. Protocol findings have empty chemical targets. These hints
never choose the next round or execute a tool. The Agent writes its campaign
decision separately, using the round facts and the authorized budget.

Compatibility fields are `decision: "agent_review"`, `decision_outcome: null`,
`goal_met: null`, and `next_round_instructions: null`. Do not interpret the
compatibility decision as convergence or a stop instruction. `max_rounds`,
`pass_rate_threshold`, and `min_candidates` remain accepted historical inputs
but are ignored for assessment and nomination.

The structure evaluator's `selection_manifest.csv` (or embedded authoritative
selection) fixes the representative pose before property filtering. Preserve
that exact `pose_id`, coordinates and lineage during coordination, including
when a backup pose has a higher StepScore. Reject an invalid manifest instead
of silently replacing its pose. Property evidence remains joined by exact
canonical pose ID; changing the representative requires regenerating its
dependent property evidence first.

For legacy reports without an authoritative selection, group poses identity-first by
`chemical_identity_id`, then `chemical_id`, and use `mol_name` only as a legacy
fallback. Partition by `protein_state_id + pocket_id + round_id +
protocol_profile_id + metric_catalog_version`, with nomination candidates and
benchmark decoys in separate normalization pools. Mark and score every pose
using Step 5 categories only. Among poses that survived the Step 5 and Step 6
folds, select the highest Step 5 StepScore. A `warn` pose with a higher score
beats a `pass` pose with a lower score; QC and raw-score terms only break exact
StepScore ties. If Step 6 rejects the best pose, select the next eligible pose.
When Step 5 scores are unavailable, explicitly report the existing QC fallback.
Preserve identity, lineage, and evidence from the exact selected record.

Stage 04 prepares molecular structures and conformers without a separate
candidate-state validator. Nomination follows observed identity, pose QC,
property evidence, and the campaign protocol. Historical state-ensemble,
population-prior, and state-nomination flags do not create independent gates.
Source, compound, optional microstate, conformer, and pose identity remain
traceable through the ordinary ligand/pose lineage.

When every selected and evaluated pose explicitly fails, retain the molecule's diagnostic record
with `all_poses_failed` and its QC reasons, exclude it from candidates, and
continue coordination for the remaining molecules. No docking retry is required.
Historical `rerun_docking` requests do not block coordination; classify using
the actual QC evidence. Missing or failed tool evidence is not pose exhaustion.

After joining that pose to its property row, mark every executable raw metric
`pass`/`flag`/`fail`, then fold the marks per parent category under the active
`screening_mode` (`V1` marking / `V2` standard, the default / `V3` strict) to get
the molecule's keep/modify/reject. `screening_mode` is a step parameter; left
unset it is derived from `candidate_source`, so a generative round defaults to
`V3`, a sampled library to `V2`, and a user-supplied full set to `V1`.
Fold-rejected records never enter a percentile pool and receive no category,
step, overall score, or rank. For each finite metric among the surviving records,
calculate descending RANK.AVG percentile; `n=1` gives `0.5`, ties share average
rank, and missing values leave only that candidate's numerator and denominator.
Explicit target-range metrics first use the hard/soft trapezoid desirability.
Weight-zero metrics still mark and report their percentile but contribute zero.
Redocking RMSD and the retrospective screening benchmark are the deliberate
exception: they are protocol gates that invalidate a run rather than a molecule,
so they stay single-level and no category fold applies.

Within each category, preserve signed catalog weights and compute a missing-
aware 0–100 score. Rank category scores again by descending RANK.AVG, then use
parent-category weights to produce StepScore and OverallScore. Vina-family
energies retain a total effective child-weight cap of 5. Planned parent-category
coverage — `Σ_available W'_c / Σ_planned W'_c`, the categories this round
actually planned — controls release: `>=0.50` formal, `[0.10,0.50)` provisional,
and `<0.10` no OverallScore/Rank plus a backtrack request. A category a round
never planned (no enabled
tool produced a rule for it) leaves numerator and denominator together and
cannot demote anything. Coverage affects `rank_status`, `coverage_gate`, and
nomination; it never rejects a molecule or a pose and never rewrites the
molecule's own classification, which reads
`metric_gate`. OverallRank is an average-tie descending rank over
only the score-released records the fold did not reject, inside the same cohort
and
candidate/decoy normalization pool. `rank` mirrors OverallRank; rejected or
blocked records stay in audit outputs with `rank`, `overall_rank` and
`weighted_score_rank` all null. Apply final nomination blocks before publishing
RANK.AVG over the remaining released scores, so blocked rows cannot consume
rank positions. Physical artifacts group
these non-comparable partitions instead of constructing a single global list. If
a transitional artifact mixes schema 2 rows with legacy flat `weighted_*` rows,
the legacy rows remain audit-only with `rank=null`; their old totals never enter
an OverallRank ordering.

The actual cohort partition key always includes protein state, pocket, round,
protocol profile, and metric-catalog version. A supplied `cohort_id` is kept as
declared provenance and cannot override those components. Each finite StepScore
also receives a local descending RANK.AVG, exposed both under
`step_scores.<step>.rank_avg` and as flat `step4_score/step4_rank` through
`step6_score/step6_rank` fields.

The report publishes three coverage layers separately, per category and per
step: `category_weighted_coverage` (the planned ratio, the one that gates),
`catalog_category_weighted_coverage` (the same numerator over the full static
catalog, information only — it never penalizes a profile for tools this round
did not plan), and metric-level planned coverage. Catalog coverage uses original
table weights; session overrides and family-cap scaling do not change it.
`planned_metric_coverage_incomplete` follows the [metric policy](../METRIC_POLICY.md)
per-category repair floors: PHY and ABS below 0.60, and TOX below 0.50. One
missing endpoint inside a well-covered category is reported
(`planned_metric_evidence_incomplete`, `planned_but_failed_metric_ids`, the
failure registry) without demoting the molecule. A whole active step with no
evidence caps release at `provisional`. Deprecated `weighted_*` fields are
aliases of the new overall fields for one schema transition, not the former flat
composite.

A metric no molecule in the cohort resolved is only dropped from
`Σ_planned` when the upstream step declared the tool skipped
(`steps_skipped`/`operational_errors`, see the contract below). Absent that
declaration it stays planned-but-failed: coverage falls, and the metric is named
in the failure registry rather than silently scoring 1.0.

`ACT_*`/`DEC_*` labels have two different meanings:

- `ACT_*` and `DEC_*` together are used for retrospective screening validation, including AUROC and enrichment factors.
- `DEC_*` records are benchmark decoys and must not be nominated as candidates.
- `candidate_molecules`, `candidate_ranking.json`, and the round report candidate table exclude `DEC_*`, and also exclude `reject` rows — a molecule that did not advance has no rank.
- `molecules` in `coordinator_report.json` remains the full evaluated set for auditability and benchmark statistics, rejects included with `rank=null`. The tables in `published/07_final_report/` mirror it row-for-row.
- Candidate and decoy normalization profiles are fit separately. Decoys also
  cannot drive bottleneck detection, modification targets, or failed-scaffold
  exclusions.

Benchmarking is opt-in, including standalone calls. Enable
`tools.screening_benchmark.enabled=true` only when the input contains both
verified classes; see `../pipeline-entry/references/benchmark-labels.md`.
If explicitly enabled but either class is absent, the benchmark gate is
`modify` and blocks nomination when quality gates are required. Missing labels
must not silently remove a requested protocol gate.

Molecule-local docking failures remain in `molecules` with
`classification=not_evaluated`, null scores/ranks, and their ledger reason.
They never enter normalization or candidate ranking. Labeled failed molecules
remain in the expected benchmark cohort: quarantine must not improve AUROC/EF
by dropping difficult active/decoy records from the denominator.

## Molecule Classification

| Classification | Criteria |
|----------------|----------|
| **keep** | Current shared `metric_gate` is keep: every evaluated category has strictly more than 50% `pass` labels, and no category trips its reject threshold. |
| **modify** | Current fold has a flag/fail but no category trips its mode. |
| **reject** | Current category fold trips, or an explicit structural/identity failure is carried into the shared gate. |

Final `classification` equals the current `metric_gate`. Legacy and intrinsic
labels remain audit fields and cannot override this result; nomination and rank
publication are independently blocked by missing evidence or protocol errors.
The numeric V2/V3 count OR ratio formula has neither a minimum fail-count
condition on its ratio branch nor an observed-coverage guard. No applicable
labels gives `not_evaluated`. V1 only preserves structural/identity rejections.
Geometry, SASA, MD, Arpeggio and fingerprint evidence must survive the pose-to-
molecule join; weighted scalar rules are normalized with the rest of the cohort.

Experiment quality-gate results do not rewrite an intrinsic molecule
classification. Protocol `modify` blocks nomination and marks evidence for
review; protocol `reject` marks the round `invalid` and publishes non-destructive
`rollback_instructions`. The Agent plans any recovery. Protocol/data-review
diagnostics keep `target_molecules=[]`; impacted IDs are audit-only
`affected_molecules`. Blueprint protocol failures expose `ranking_blocked` and
`ranking_blocked_by` without inventing a chemical defect.

## Protocol Quality Classification

For `site_selection_strategy=provided_box`, read the same Step 1
`validation.passed` field. Its gate covers fixed-box geometry and selected-chain
overlap, with the box/receptor hashes retained in `details.input_evidence`.
Detector consensus is not required; a pass does not establish binding-site or
pose correctness. Missing/failed reports retain the ordinary modify/reject
semantics below.

Protocol gates use the same inclusive three-tier boundary convention but have
weight 0 and never enter a molecule's weighted total:

| Protocol metric | keep | modify | reject |
|-----------------|------|--------|--------|
| Native-ligand redocking RMSD | `<= 2.0 Å` | `(2.0, 2.5] Å` | `> 2.5 Å` |
| Retrospective AUROC | `>= 0.70` | `[0.50, 0.70)` | `< 0.50` |
| Retrospective EF1% | `>= 3` | `[1, 3)` | `< 1` |
| Benchmark sample sufficiency | Meets both configured active and decoy minima | Below either minimum | — |

A report that was not supplied, or a path that does not exist on disk, is
`modify` with `passed: null` — not `reject`. The protocol may be fine; the
coordinator simply has no verdict to read. That produces a review finding for the missing evidence, without invalidating
a run whose only defect is a missing input path. A
`passed` flag that is present and literally `false` is still `reject`.

Aggregate these gates with reject dominance. A reference-dependent gate that
was not performed because no valid reference exists is `not_applicable`, not a
pass or failure. If every required quality gate is `not_applicable`, the
aggregate remains `not_applicable` and requests applicable evidence instead of
silently becoming keep. Report PR-AUC together with active prevalence and its relative
lift (`PR-AUC / prevalence`), and report the top-ranked decoy fraction as audit
evidence; these audit signals have no universal hard cutoff and weight 0.

## Bottleneck Classification

| Bottleneck | Type | Implication |
|------------|------|-------------|
| `steric_clashes` | scaffold-level | Need different chemical series |
| `posebusters_violations` | scaffold-level | Need different chemical series |
| `blueprint_mismatch` | modification-amenable | Add functional groups for confirmed missing constraints |
| `blueprint_unresolved` | protocol/capability | Re-evaluate with a profiler capable of the required atom/role/geometry field |
| `high_strain` | modification-amenable | Simplify scaffold |
| `low_qed` | modification-amenable | Improve drug-likeness |
| `pains` | modification-amenable | Replace reactive substructure |
| `lipinski` | modification-amenable | Reduce size/lipophilicity |

## Output Files

| File | Description |
|------|-------------|
| `coordinator_report.json` | Full report: cohort identity, metric policy, raw gates, RANK.AVG percentiles, category/step/overall scores, three coverage layers, rank status, molecule classifications, and nomination validity |
| `round_summary.json` | Aggregate statistics for this round |
| `candidate_ranking.json` | Non-decoy `keep`+`modify` candidates grouped by cohort/pool and ordered by local released OverallRank, with exact selected-pose metrics; unranked audit rows follow within their partition. Rejects never enter it — gate first, then score |
| `round_report.md` | Deterministic per-round non-decoy candidate report with a separate Top-25 `keep`+`modify` table per cohort/pool |

The on-disk full report stores round history through `files.round_history_json`
rather than embedding it again. Its `candidate_molecules` contains compact
identity/rank/classification/nomination summaries; read
`files.candidate_ranking_json` for complete candidate evidence. Nested property
copies of `metric_assessments`, `category_scores`, `step_scores`,
`property_category_scores`, and `property_step_scores` are replaced by
`property.scoring_report_file`, which references the upstream property report.
The coordinator's top-level current-round scoring and raw property evidence
remain available; this storage change does not recalculate or remove them.

### The published tree

The runner also maintains a `pipeline_partial` view before full finalization,
including on failure or when coordination/reporting are disabled. This view
contains stage/file/log indexes and `07_final_report/README.md` plus
`publication_status.json`; it does not contain the scored tables below.
Check `publication_mode` before treating a published path as a final deliverable.

The coordinator keeps the resumable `step0..step7` working tree unchanged and
publishes the round's reader-facing tree under `{round_dir}/published/`.  The
tracked contract and migration procedure are in
[references/output-layout.md](references/output-layout.md).

| Path | Description |
|------|-------------|
| `published/00_inputs/` through `03_molecule_generation/` | Real run inputs and early-stage artifacts, indexed by relative path/hash without fabricated score tables |
| `published/04_sm_property/` | Pre-docking property evidence plus the scorer's Step 4 tables/categories/tools |
| `published/05_docking/` | Docking evidence plus the scorer's Step 5 tables/categories/tools |
| `published/06_post_docking/` | Pose-evaluation evidence plus the scorer's Step 6 tables/categories/tools |
| `published/07_final_report/final_report.tsv` | The 19 user-facing columns specified in [output layout](references/output-layout.md), best row first — the entry point |
| `published/07_final_report/final_report_full.tsv` | Every identity, gate, step, category and ranking column for the same rows |
| `published/07_final_report/all_metrics_table.tsv` | Full catalog width with each metric's value plus canonical `__label` (`pass`/`flag`/`fail`) and compatibility `__gate` columns |
| `published/07_final_report/metric_catalog_and_scoring.md` | The rules and the scoring profile these tables were produced with, frozen |
| `published/07_final_report/run_manifest.json`, `lineage_graph.json` | What was written, its digest, and which table each was derived from |
| `published/07_final_report/failure_registry.jsonl` | Written only when something failed; a missing value stays visible as a missing value |
| `published/0{4,5,6}_*/01_categories/*/` | `category_table.tsv`, `category_report.tsv`, `source_map.json` per category |
| `published/0{4,5,6}_*/02_tool_runs/*/` | `tool_table.tsv` plus `_internal/` for every physical tool run that produced a value this round |

Three conventions the contract leaves open:

- An empty cell means "not measured" — never `0`, never the string `NA`.
- A column of blanks is dropped from a step, category or tool table for the same
  reason a tool that did not run gets no directory. `all_metrics_table.tsv`
  keeps the full catalog width, so nothing disappears from the run.
- Paths in a manifest or `source_map.json` are publication-relative references. Poses and
  trajectories are never copied; the step directories keep owning them.
- The eight top-level names are reader-facing scientific stages.  Manifests map
  them to runner sources (`step0`, `step1`, `step2`, `step2_5` + `step5`,
  `step3`, `step4`, `step6` + `step7`) while preserving the scorer's internal
  `step4`/`step5`/`step6` IDs and column names.

### Why a cell is empty

`failure_registry.jsonl` is the ledger that turns a blank cell into a reason. One
JSON line per gap, carrying `stage`, `tool_id`, `reason`, the affected `row_ids`,
and a `detail.kind` that says whether the gap was a fault or a choice:
`tool_error`, `stage_error`, `step_skipped`, `step_disabled`,
`missing_step_report`, `empty_step_report`, `partial_metric_coverage`,
`quality_gate_failed`, `quality_gate_not_applicable`, `row_blocked`. Tool entries
also list which metric columns were left blank, checked against the published
rows rather than assumed. No entry is inferred from an absence alone, and a clean
round writes no file at all.

Two of those kinds depend on upstream skills, and only they can report them: a
tool that failed for the whole cohort leaves no per-row trace, so with no row
holding a value the scorer cannot tell that gap apart from a metric this build
never planned. Producing skills must therefore keep declaring
`steps_skipped: [{step, reason}]` and `operational_errors: [{tool, reason}]` in
their reports; the coordinator routes each declared name onto the tool that owns
the metric and files it here.

When ACT_/DEC_ labels are present, also report AUROC with a bootstrap 95%
interval, EF1%/EF5%/EF10%, PR-AUC, active prevalence, PR-AUC relative lift, and
the top-ranked decoy fraction. The selected native docking metric must cover
every labeled ACT_/DEC_ record, and all selected scores must carry the same
complete protocol identity. Missing score coverage or missing/mixed protocol
identity is a protocol `modify`; the coordinator must not compute a favorable
benchmark from the remaining subset. Insufficient benchmark size is `modify`;
AUROC or EF1% below its hard boundary is a protocol `reject`. None of these
cases may be presented as validated evidence.

## Round Report Logic

`round_report.md` is written by `reflective-coordinator/run.py::_write_round_report`; it is a per-round deterministic report, not the project-level final conclusion. The project-level `final_report.md` is written only after the Agent has reviewed all rounds. The round report must keep validation evidence separate from hit nomination:

1. **Executive Summary**: round assessment, quality-gate status, total screened molecules, non-decoy nomination pool size, decoy keep count, and primary bottleneck.
2. **Retrospective Screening Validation**: AUROC, bootstrap CI, benchmark size, enrichment factors, PR-AUC with prevalence/relative lift, and top-ranked decoy fraction from the full ACT_/DEC_ set.
3. **Quality Gates**: native redocking, pocket validation, blueprint validation, and screening benchmark.
4. **Candidate Pool Notes**: explicitly state that `DEC_*` decoys are excluded from nomination, and summarize high-clash/high-strain risk among keep candidates.
5. **Top Ranked Candidates (keep + modify)**: separate Top-25 tables per cohort/pool over the non-decoy molecules that advanced — `keep` and `modify`, which is the pool the scores were computed over — using the exact selected pose. Rejects are excluded; a separate section lists protocol-blocked keep candidates. Never truncate one cohort because another cohort filled a global quota. Include `mol_name`, SMILES, local OverallRank, primary metric/value, exact required-group status, weighted desired-Blueprint score, steric clash count, strain energy, QED, SA, and risk flags. Keep legacy GNINA columns only in machine-readable CSV compatibility fields.
6. **KEEP and MODIFY Top 1–5 for Review**: add separate tables per cohort/pool
   and classification from all non-decoy candidate rows, including blocked
   rows. Order by existing OverallScore without recalculation; only finite
   scores enter Top 5. Label the sequence `Review #`, preserve nullable
   OverallRank and nomination status, and state the source round. Separate
   molecular category findings from missing evidence/protocol blockers. These
   tables are review-only and never nominate MODIFY or blocked KEEP records.
7. **Recommended Next Steps**: concise follow-up actions, especially re-docking/visual inspection for high-ranking molecules with clashes or strain.

Read QED and SA from `property.properties.QED` and `property.properties.SA_score` when present; fall back to legacy lowercase keys only for old reports.

## Usage

```bash
# Standard usage after evaluation + property filter
python skills/bio_research/reflective-coordinator/run.py \
  --input '{"evaluation_report_file": "step4/evaluation_report.json", "property_filter_report_file": "step5/property_filter_report.json", "round_number": 1}' \
  --output-dir ./step6/

# Without property filter (evaluation only)
python skills/bio_research/reflective-coordinator/run.py \
  --input '{"evaluation_report_file": "step4/evaluation_report.json", "round_number": 2}' \
  --output-dir ./step6_r2/
```

## Agent Campaign Review

1. Read `round_assessment`, the raw quality gates, candidate evidence, and
   `diagnostics` after each completed round.
2. Treat protocol or missing-evidence problems as protocol work; do not convert
   them to chemical modification targets.
3. The Agent chooses the R1–R3 transition, candidate retention, parameter changes,
   earliest affected step, and continue/pause/stop action under the campaign budget.
4. Reuse historical artifacts only within their lineage/profile contract. Keep
   raw round history for audit; never compare different-profile percentiles directly.
5. Eligible candidates can be inspected while the campaign continues. Round
   assessment validity is independent of campaign completion.

## Dependencies

None — deterministic Python evidence assessment. Reads JSON files produced by upstream skills.

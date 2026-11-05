# SBDD reader-facing output layout

This is the tracked contract implemented by `scripts/publish_layout.py`.  It is
the publication subset of the larger SBDD directory proposal: execution and
publication are separate concerns.

## Roots

For a pipeline round at `<project>/roundN/`:

```text
roundN/
├── step0/ ... step7/       # resumable working/audit tree
├── logs/                   # runner logs
└── published/              # reader-facing tree described below
```

The runner must keep returning `project_dir`/`working_dir` for checkpoint
resume.  `published_dir` (also returned as `user_output_dir`) is the path a UI
or a human should open. It is exposed for both `pipeline_partial` evidence
snapshots and `pipeline_final` deliverables. Callers must inspect the publication
mode and execution/scientific status, not infer success from the path. Only
`pipeline_final` requires successful coordination/reporting and sealed sources.

## Partial execution snapshots

The runner refreshes `published/` at step boundaries, confirmation pauses,
limited-plan completion and caught failures. The eight stage directories each
contain a manifest and `artifact_index.tsv`; `00_inputs/log_index.tsv` points to
existing logs. `07_final_report/README.md` and `publication_status.json` show
execution status and errors. `final_report.html`, `final_report.md` and
`final_report_data.json` in that directory provide an offline execution report
with completed/failed/unfinished stages, artifact links and missing-file flags.
This report is available even if validation fails before the first experiment
stage; it remains explicitly partial and never fabricates ranking or scores.

The current runner manifest selects visible results. Pending, disabled and
rewound steps do not expose old files merely because their directories exist.
Partial indexes record existence/size with `source_integrity=unverified_references`;
only the generated view has exact file-set/hash verification. Raw files remain
in the working tree and may change on resume. Full publication retains its
strict source verification and replaces the partial view after finalization.

`republish_existing_round.py --round-dir <roundN> --partial` builds a missing
partial view for an idle historical round without changing its source files.

## Full publication tree

```text
published/
├── 00_inputs/
│   ├── input_manifest.json
│   └── artifact_index.tsv
├── 01_structure_prep/
│   ├── step_manifest.json
│   └── artifact_index.tsv
├── 02_pocket_analysis/
│   ├── step_manifest.json
│   └── artifact_index.tsv
├── 03_molecule_generation/
│   ├── step_manifest.json
│   └── artifact_index.tsv
├── 04_sm_property/
│   ├── step_table.tsv
│   ├── step_report.tsv
│   ├── step_manifest.json
│   ├── artifact_index.tsv              # post-pose property source files
│   ├── 00_prefilter/
│   │   ├── phase_manifest.json
│   │   └── artifact_index.tsv          # pre-docking property source files
│   ├── 01_categories/<category>/
│   │   ├── category_table.tsv
│   │   ├── category_report.tsv
│   │   └── source_map.json
│   └── 02_tool_runs/<tool>/
│       ├── tool_table.tsv
│       └── _internal/
│           ├── tool_summary.json
│           ├── internal_index.tsv
│           └── provenance.json
├── 05_docking/                         # same three scoring layers
├── 06_post_docking/                    # same three scoring layers
└── 07_final_report/
    ├── all_metrics_table.tsv
    ├── final_report.tsv
    ├── final_report_full.tsv
    ├── metric_catalog_and_scoring.md
    ├── run_manifest.json
    ├── lineage_graph.json
    ├── failure_registry.jsonl          # only when gaps exist
    ├── step_manifest.json
    └── artifact_index.tsv
```

The mandatory public names are `final_report.tsv`, `final_report_full.tsv`,
`metric_catalog_and_scoring.md`, `*_table.tsv`, and `*_report.tsv`. A tool-run
directory is created when at least one row has a measured metric owned by that
physical tool run or when an attempted execution must be audited despite
producing no scalar value. A tool with neither measurements nor execution
evidence is omitted.

Each `tool_summary.json` separates metric coverage from execution. Per-metric
`rows_unmeasured` is the number of blank measurements and is split into
`rows_not_planned`, `rows_planned_missing`, and other unavailable rows, with the
full `availability_status_counts` retained. The legacy `rows_missing` field is
a deprecated alias of `rows_unmeasured`; neither may be presented as a
tool-failure count. AiZynthFinder summaries additionally publish an
`execution_summary` derived from row-level `attempted/status/reason`, including
status/reason counts and contract inconsistencies. Only a consistent
`attempted=true,status=error` contributes to `rows_operational_errors`.

Several planner switches may share one physical archive.  For example,
Lipinski, QED, and SA metrics are published once under `rdkit_descriptors/`;
their columns must be combined into one `tool_table.tsv`, and the step manifest
must list that directory only once.

The eight names above are publication stages, not the runner's historical step
numbers.  Their source mapping is explicit in every stage manifest:

| Publication stage | Resumable working source |
|---|---|
| `00_inputs` | pipeline config and resolved run plan |
| `01_structure_prep` | `step0` (`structure_prep`) |
| `02_pocket_analysis` | `step1` (`pocket_analysis`) |
| `03_molecule_generation` | `step2` (`molecule_generation`) |
| `04_sm_property` | `step2_5` prefilter plus `step5` full property filter |
| `05_docking` | `step3` (`docking`) |
| `06_post_docking` | `step4` (`pose_evaluation`) |
| `07_final_report` | `step6` coordination plus `step7` reporting when present |

The scorer's `step4`, `step5`, and `step6` identifiers remain unchanged inside
tables and manifests: they are scoring-policy IDs, not execution-order labels.
The two property passes also remain distinct because they operate on different
populations; `00_prefilter/` is an evidence index, never merged into the
post-pose score table.

## References and identity

- Working artifacts are not moved or renamed.  `artifact_index.tsv`,
  `source_map.json`, `input_refs.json`, asset manifests, provenance, and
  run-manifest inputs use POSIX paths relative to `published/`, including
  `../stepN/...` where needed.  JSON manifests declare
  `path_base: "published_root"`.
- Large poses, trajectories, and topologies remain owned by the working tree.
  The publication is therefore a view, not a standalone archive.
- Tables join on stable IDs: `molecule_state_id` for the small-molecule property
  stage and `pose_id` for docking/post-docking.  File order is never an identity.
- One empty cell means not measured.  It is never written as `0`, `NA`, or a
  sentinel.  `failure_registry.jsonl` explains declared gaps.
- Metric drill-down tables expose both the canonical `__label`
  (`pass`/`flag`/`fail`) and the compatibility `__gate`
  (`keep`/`modify`/`reject`), so a single metric mark cannot be mistaken for a
  molecule-level decision.
- `07_final_report/run_manifest.json` records the layout version and the size
  and SHA-256 of every generated publication artifact.
  `07_final_report/lineage_graph.json` contains unique source-to-table edges.
- Every normal publication is built in a fresh sibling staging directory,
  verified, and then recoverably swapped into place under a directory lock.
  Resume never writes incrementally over an older tree, so removed tools and
  resolved failures cannot linger; an interrupted two-rename swap restores its
  orphan backup on the next invocation.
- A publication declares `standalone`, `pipeline_partial` or `pipeline_final` mode.
  `pipeline_final` is sealed only after successful `step6/result.json` and
  `step7/result.json`; its source snapshot covers every cited working artifact,
  including externally injected inputs and the report bundle's plots/assets.
  Verification re-hashes those live sources, so a changed source invalidates
  the reader view instead of silently serving stale evidence.

## Known boundary

`checkpoints/round_snapshot_manifest.json`, standalone native archives,
per-tool command/config sidecars, and materialized protein/pocket/ligand assets
are not currently produced by this publisher.  Early-stage `artifact_index.tsv`
files describe only real source files; candidate-level score tables are not
fabricated for stages that do not produce them.
The legacy working tree remains the authoritative source for those artifacts.

Use `scripts/republish_existing_round.py --round-dir <roundN>` to create an
isolated publication for an older mixed-layout round.  It never overwrites the
source or an existing publication, requires and exclusively acquires the
runner's existing `.pipeline.lock`, verifies the source aggregate SHA-256 around
both report loading and publication, and verifies the exact artifact set
registered in the new run manifest.

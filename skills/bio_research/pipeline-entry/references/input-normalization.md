# Pre-config input normalization

`pipeline-entry` remains the SBDD workflow entry point. When a user supplies a
local receptor or a fixed screening library, inspect and normalize those files
before assembling the authoritative pipeline config. Executable behavior lives
in `../../input-normalization/`; this reference defines the orchestration and
config boundary.

## Canonical boundary

| User input | Pre-config action | Pipeline-facing value |
|---|---|---|
| receptor PDB | validate/copy | PDB path |
| receptor MOL2 | Open Babel format conversion plus atom/coordinate QC | PDB path |
| fixed SDF library | check each record; log and skip unusable molecules | usable SDF path on success |
| fixed SMI/SMILES library | Open Babel `-h --gen3d fast` per molecule, bounded parallel workers and a per-molecule timeout; log and skip failures | usable SDF path on success |
| provided box JSON (`center`/`size` or `*_angstrom`) | schema and finite-coordinate validation | canonical docking-box JSON |

The canonical format is an internal interoperability contract, not a claim that
PDB/SDF are the only scientifically valid external formats.

SMILES conversion defaults to at most 8 workers, capped by available CPUs and
record count, and 60 seconds per molecule. `library.conversion_timeout_seconds`
accepts a positive integer override; omitted or `null` uses 60 seconds, never
unlimited time. A timeout records `conversion_timeout`, skips that molecule and
allows the remaining molecules to continue. All library coordinate generation
uses `--gen3d fast`; receptor format conversion preserves existing coordinates
and does not generate new ones. Final SDF records remain in source order.
Existing SDF input is checked without coordinate regeneration. Explicit other
library formats retain the whole-file conversion fallback: the timeout applies
to that file and fails normalization, rather than skipping an individual record.

The adapter probes native fast-mode coordinate generation before processing a
library. If it detects the Open Babel 3.1.1 coordinate-writeback defect, it builds
and caches `runtime/openbabel_fast_compat` in the normalization output directory.
This preserves the fast algorithm, including distance-geometry fallback, and
fixes writeback without changing the global Conda environment. Only that path
requires a C++17 compiler, matching Open Babel development headers/library and
Eigen 3 headers. Missing dependencies or failed probes fail the normalization
run, rather than marking library molecules as unusable.

## Agent procedure

1. Read `../../input-normalization/SKILL.md` completely.
2. Decide only facts required by normalization. In particular,
   `receptor.preparation_state=raw|prepared` must come from the user/protocol;
   never infer it from a `.pdb` or `.mol2` suffix.
3. Before writing new experiment files, reserve the pipeline round with
   `../run_pipeline.py --init-project-dir --input "$(python3 -c 'import json, os; print(json.dumps({"project": {"data_root": os.environ["SBDD_DATA_ROOT"], "project_slug": "request-short-name", "round_number": 1}}))')"`.
   This accepts directory naming facts only; no scientific plan is validated or
   executed. Use a known target or request label for the slug without inferring
   protein identity. Save the returned `project_dir` for this attempt, request
   JSON/config drafts in `inputs_dir`, and invoke
   `../../input-normalization/run.py --output-dir <normalization_dir>` with the
   returned path (`<project_dir>/inputs/normalization`). Existing directories
   cannot be re-initialized; reuse a previously reserved directory directly.
   Start conversion with `terminal(background=true, cancel_on_interrupt=true,
   notify_on_complete=true, pty=false)` and an `exec` Python command. Save its
   `session_id`; use `process(action="poll", session_id=...)` or
   `process(action="wait", session_id=..., timeout=30)`. Every wait explicitly
   uses 1–30 seconds. Read `<normalization_dir>/normalization_progress.json`
   with `read_file` for completed/total counts, skipped records and active
   molecules, and report these facts while the task runs. Do not use
   `while kill -0 PID; do sleep 30; done` or another unbounded shell wait loop,
   `nohup`, or a trailing `&`. A tool wait timeout does not cancel the library;
   the normalization script owns each molecule's timeout and continuation.
4. Read the command diagnostics, `library.skipped_records` and `library.geometry_qc` in
   `normalization_report.json`. Merge `config_patch.json` only on `success=true`.
   `library.conversion_execution.command_backend` identifies native fast mode
   (`openbabel_native_fast`) or the writeback fix (`openbabel_fast_writeback_compat`).
   Keep `runtime/` probe files and, when used, `fast_compat_build.json` plus the
   compiled executable with the other normalization outputs.
   A bad molecule is recorded and skipped even if Open Babel returned 0. The
   remaining usable molecules continue without replacement sampling. Preserve
   skipped entries as preprocessing failures with null scores in the final
   report, and state input, retained and skipped counts. The unfiltered SDF is
   retained separately. Empty usable libraries and infrastructure errors still
   fail; a failed retry removes a stale config patch.
   No manifest, ID assignment, state enumeration or separate approval is added.
5. Assemble and validate the full pipeline config only after successful
   normalization, then follow the existing plan authorization process. Execute
   `../run_pipeline.py --input-file <full-config.json> --project-dir <project_dir>`
   with the reserved path. Keep this CLI override out of the saved config so
   future attempts cannot accidentally target the old round. Conversion files
   and pipeline outputs now share one round directory from their first write.
   Do not store source artifacts inside `published/`, which is rebuilt by the
   runner. Existing `_staging` inputs remain valid and must not be moved while
   referenced by an older experiment.

## Config mapping

- `raw` receptor: set `context.receptor_file` to the normalized PDB, disable
  Step 0 intake, and let policy choose the remaining preparation tools.
- `prepared` receptor: inject `input_artifacts.prepared_receptor_noH` and disable
  Step 0 plus every Step 0 tool.
- fixed normalized library: inject `input_artifacts.library`, then disable
  molecule generation plus every generator tool. External libraries use ordinary
  3D SDFs; normalization does not enumerate chemical states or create state seals.
- provided docking box: inject `input_artifacts.docking_box`, enable pocket
  analysis with `params.site_selection_strategy=provided_box`, and select
  `docking.params.search_space_source=docking_box`. Step 1 uses the unchanged
  supplied box to characterize receptor overlap and emits the report consumed by
  coordination in the same experiment. Detector/profiler tools remain optional
  explicit choices; no auxiliary experiment or predicted replacement box is needed.

Property prefiltering is not disabled automatically: whether an uploaded fixed
library may be filtered before docking is a scientific/protocol decision. A
benchmark requiring every compound to be screened keeps property prefiltering
off. Conversion diagnostics are in the Open Babel stdout/stderr; docking
accounting is recorded by the docking step.

## MAPK1 smoke example

For an oracle public task containing `receptor/receptor.mol2`,
`library/library.smi`, `library/compound_manifest.jsonl`, and
`oracle_docking_box.json`, use a request shaped like:

Replace `<INPUTS_DIR>` with the actual input directory before saving this JSON;
JSON does not expand shell variables.

```json
{
  "receptor": {
    "path": "<INPUTS_DIR>/receptor/receptor.mol2",
    "format": "auto",
    "preparation_state": "prepared"
  },
  "library": {
    "path": "<INPUTS_DIR>/library/library.smi",
    "format": "auto",
    "workers": 8,
    "conversion_timeout_seconds": 60
  },
  "docking_box": {
    "path": "<INPUTS_DIR>/oracle_docking_box.json"
  }
}
```

Use only files made visible by the selected public campaign. An oracle/reference
docking box belongs only to a campaign that explicitly exposes it. The MAPK1
`predicted-pocket/public` campaign does not contain this file and must run an
enabled pocket predictor; the separate `oracle-pocket/public` control exposes
it and uses `provided_box` analysis without pocket prediction. Never read a private/hidden box or activity
label while solving either campaign.

For the predicted campaign, omit `docking_box` from the normalization request.
Use `site_selection_strategy=detector_consensus` with DeepPocket, fpocket,
P2Rank, and pyKVFinder enabled for Top-3 local fusion. Declare pocket memory mode
and compatible query/index inputs; fixed DeepPocket#1 rank fusion is the fallback
when memory is unavailable or explicitly disabled. When UniDock is the declared engine, set the explicit
Step 1 parameter `detector_box_padding=1.0`; the MAPK1 smoke receptor otherwise
produces a default-padded detector box larger than UniDock's current GPU grid
buffer. This parameter changes only the envelope around the same selected
residues and is recorded in the resolved plan. Treat detector quality and box
size as scientific/runtime facts in Gate A, not as hidden benchmark hints.

## Complete fixed-library ranking

SMILES names and conversion output are handled directly by Open Babel.
Normalization reports input/retained/skipped counts and skipped record names/indices;
it does not assign IDs, rewrite valid records or generate an input ledger.

Input normalization does not rank molecules. After docking, join
`docking_scores.json` to `docking_molecule_ledger.json` by `compound_id` (or the
unambiguous `mol_name` alias):

- a scored UniDock pose becomes `status=success`, `score=docking_score`,
  `engine=unidock`;
- multiple poses remain separate records; the benchmark normalizer selects the
  best score for the declared direction;
- every ledger row with no finite pose remains one `status=failed` record with
  a null score;
- normalization-skipped compounds also remain `status=failed` with a null score;
  recover these from `normalization_report.json`, since they never entered the
  docking ledger;
- use the benchmark input manifest when assembling the final ranking; it is not
  an input to library conversion.

Then use the benchmark's public `normalize_agent_run.py` with
`--score-direction minimize`. It validates the public manifest ID set, assigns
continuous ranks, puts failed compounds after scored compounds, and writes the
required `compound_ranking.csv`. For MAPK1, require exactly 1001 unique rows
before submission. The public contract is 10 labelled actives, 990 labelled
inactives, and one unlabelled native sentinel—not 11 labelled actives. Never
inspect private labels while producing the ranking.

Prepared receptor normalization preserves the original structure and water coordinates,
produces a separate noH derivative plus prep manifest, and still requires final
`docking.docking_prep` before engine execution. See the preparation contract.

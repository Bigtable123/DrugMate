---
name: input-normalization
description: Prepare external receptors, libraries and docking boxes.
---

# SBDD Input Normalization Skill

Normalize external receptors, fixed molecular libraries and docking boxes before
`pipeline-entry` assembles the authoritative config. This adapter prepares input
files; pocket selection, molecule generation and ranking belong to other steps.

## When to Use

Use when the user supplies a local receptor, fixed SMI/SMILES/SDF library, or
box JSON. Read `../pipeline-entry/references/input-normalization.md` for wiring.

## Prerequisites

Open Babel and RDKit must be available in the registered `drugmate` scientific
runtime. Deployment settings select its installed Conda name or prefix.
RDKit reads the output for basic coordinate checks; Open Babel owns conversion.
Before converting a library, the adapter checks native `--gen3d fast` with a
small ethane probe. Some Open Babel 3.1.1 builds fail to write the generated
coordinates back. Only when this defect is detected, the adapter builds a local
compatibility executable, requiring a C++17 compiler, matching Open Babel
development headers/library and Eigen 3 headers. It keeps the fast algorithm and
its distance-geometry fallback, and fixes coordinate writeback only. The global
Conda installation is unchanged. Missing build dependencies or a failed runtime
probe fail normalization before processing the library; they are not bad-molecule
failures to skip.

## How to Run

For a new pipeline run, first use `terminal` to reserve its output directory:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/run_pipeline.py" \
  --init-project-dir \
  --request-id '<stable ID for this request; reuse on initialization retries>' \
  --input "$(python3 -c 'import json, os; print(json.dumps({"project": {"data_root": os.environ["SBDD_DATA_ROOT"], "project_slug": "mapk1-screen", "round_number": 1}}))')"
```

This returns `round_dir` (`project_dir` for compatibility), `project_root`,
and `normalization_dir` without starting scientific
steps. Save the normalization request as `<project_dir>/inputs/input_normalization.json`,
then start conversion with
`terminal(background=true, cancel_on_interrupt=true, notify_on_complete=true, pty=false)`:

```bash
exec "${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/input-normalization/run.py" \
  --input-file "<project_dir>/inputs/input_normalization.json" \
  --output-dir "<normalization_dir>"
```

Replace placeholders with the returned paths. Assemble the final config from
`config_patch.json`, and execute the pipeline with `--project-dir <project_dir>`
so normalization and pipeline outputs stay in the same run. For an already
allocated round, use its `inputs/normalization/` directly. A standalone conversion
may still use any explicit `--output-dir`. See
`../pipeline-entry/references/input-normalization.md` for the full handoff.
Use the returned round directory, never its parent `project_root`. Preserve the
request ID across initialization retries; obtain it from the request/message ID
or generate one token once. When continuing an existing task after interruption,
recover its initialization record and round directory instead of allocating a
new project. Distinct requests for the same PDB remain separate experiments.

SMI records use isolated Open Babel processes in a bounded worker pool. The
adapter resolves the registered environment once per library (at most one
`conda run`), then calls the checked native or compatibility executable directly
for every record. The native command has this form:

```bash
"${SBDD_ENV_DRUGMATE_PREFIX}/bin/obabel" \
  -ismi "${SBDD_DATA_ROOT}/current_request/inputs/one_record.smi" -osdf \
  -O "${SBDD_DATA_ROOT}/current_request/inputs/one_record.sdf" -h --gen3d fast
```

The adapter selects the actual executable; do not hardcode the example path or
bypass its functional check. Defaults are at most 8 workers, capped by available
CPUs and record count, with one OpenMP/BLAS thread per worker. Results are merged
in original input order even when conversions finish out of order. A failed
record does not stop later molecules. Do not create a separate shell/Conda loop
around individual molecules.

Save the returned `session_id`. Use `process(action="poll", session_id=...)`
for status, or `process(action="wait", session_id=..., timeout=30)` for a
bounded wait; every wait must explicitly set a timeout from 1 to 30 seconds.
After waiting, read `<normalization_dir>/normalization_progress.json` with
`read_file`: `status`, `completed_records` / `total_records`, `retained_records`,
`skipped_records`, and `running_records` show progress. Running entries include
the molecule name, source index, record number and elapsed seconds. Statuses
distinguish `starting`, `preparing_environment`, `converting`, `completed`,
`failed` and `cancelled`.
Use these progress facts in user updates and rely on the completion notification
for the final result. A wait timeout only ends that observation window; it does
not mean conversion failed and must not trigger killing the whole library.

Do not use `while kill -0 PID; do sleep 30; done`, any other unbounded shell
PID-wait loop, `nohup`, or a trailing shell `&`. Do not wrap conversion in
`execute_code`: its outer time limit can interrupt the managed task. The script
enforces each molecule's timeout and skips it without stopping later molecules.

## Quick Reference

Provide at least one of `receptor`, `library`, or `docking_box`:

Replace `<INPUTS_DIR>` with the actual input directory before saving this JSON;
JSON does not expand shell variables.

```json
{
  "receptor": {
    "path": "<INPUTS_DIR>/receptor.mol2",
    "format": "auto",
    "preparation_state": "prepared",
    "coordinate_tolerance": 0.002
  },
  "library": {
    "path": "<INPUTS_DIR>/library.smi",
    "format": "auto",
    "workers": 8,
    "conversion_timeout_seconds": 60
  },
  "docking_box": {"path": "<INPUTS_DIR>/oracle_docking_box.json"}
}
```

- Receptors: PDB or MOL2. `preparation_state=raw|prepared` is a user/protocol fact;
  ask only when it is unclear and changes the plan. A suffix does not imply it.
- Libraries: `path` and optional `format` (default `auto`). Auto treats `.sdf`
  and `.sd` as SDF and other suffixes, including `.txt`, as SMILES. An explicit
  format is passed to Open Babel without content-based cross-checking.
- SMI input runs `obabel -ismi ... -osdf -O ... -h --gen3d fast` per molecule.
  SDF input is checked record by record. Failed molecules are logged and skipped;
  remaining records keep their original names, order and SDF bytes. There is no
  automatic regeneration, replacement sampling, property filtering or state
  enumeration; old preparation/manifest options have no effect.
- SMILES execution options: `workers` is an integer from 1 to 64 (default 8;
  capped by available CPUs and record count). Use 1 for serial execution.
  `conversion_timeout_seconds` is a positive integer per molecule; omitted or
  `null` means 60 seconds. Unlimited molecule conversion is not supported.
  On timeout the script stops that molecule, records `conversion_timeout` and
  continues with the rest. All library conversions that generate coordinates
  use `--gen3d fast`; receptor conversion preserves existing coordinates and
  does not generate new ones.
- Explicit non-SMILES, non-SDF library formats use the existing whole-file
  conversion fallback with `--gen3d fast`. Its 60-second default (or positive
  override) applies to the whole file; timing out fails normalization instead
  of attributing a timeout to individual records. Existing SDF input is checked
  without regenerating coordinates.
- Boxes: `center`/`size` or `center_angstrom`/`size_angstrom`, each with three
  finite numbers and strictly positive sizes. Supplied boxes are preserved.

## Procedure

1. Run normalization into the requested output directory. Full conversion output
   (or the supplied SDF) is preserved in `library/library_unfiltered.sdf`.
   Only usable records are written to `library/library_3d.sdf`.
   While running, monitor `normalization_progress.json` and the managed process
   with bounded tool calls. Progress follows completed work even when a slower
   earlier molecule is still running; final SDF records retain source order.
2. Read `library.skipped_records`, the input/retained/skipped record counts,
   `library.conversion_runs` (per-SMILES command diagnostics),
   `library.conversion_execution` (worker count, environment setup, command
   backend and total conversion wall time), and
   `library.geometry_qc` in `normalization_report.json`. Conversion timeouts,
   conversion failures, empty/unreadable/incomplete records, missing or non-finite coordinates,
   collapsed coordinates and zero-length bonds skip the affected molecule.
   The report identifies it by original source index, one-based record number
   and name (plus source line number for SMILES). `retained_source_indices` maps
   output order back to the original input. Geometry QC describes the retained
   library, not a full conformer-quality or stereochemistry assessment.
   Continue when at least one usable molecule remains. An empty usable library
   or a file/runtime setup error still fails normalization.
3. Only on `success=true`, merge `config_patch.json`: fixed libraries map to `input_artifacts.library`
   with molecule generation disabled. Raw receptors map to `context.receptor_file`
   with intake disabled; prepared receptors map to
   `input_artifacts.prepared_receptor_noH` with Step 0 disabled. Supplied boxes
   map to `input_artifacts.docking_box`, enable pocket analysis with
   `site_selection_strategy=provided_box`, and set
   `docking.params.search_space_source=docking_box`. Step 1 characterizes the
   selected chain inside the supplied box and produces its validation report in
   the same round. The patch leaves detector/profiler tools off; select any
   supplemental tools explicitly. Do not create an auxiliary pocket run.
   When a prepared receptor and box are both supplied, the patch sets
   `context.chain` only if exactly one named receptor chain overlaps the box.
   Read `docking_box.chain_selection` for the per-chain heavy-atom counts.
   Preserve any chain explicitly selected by the user when merging the patch;
   validation checks that selection instead of replacing it. Multiple overlapping
   chains, no overlap, or missing chain identifiers leave chain selection for
   review and add warnings/plan requirements without failing format conversion.
   Never assume chain A or select whichever chain has the most atoms. Raw
   receptors require chain review after preparation.
4. Carry skipped molecules into the final report as preprocessing failures with
   null scores. State how many were supplied, retained and skipped. Do not
   regenerate or replace skipped molecules just to restore the requested count.

## Pitfalls

Open Babel can exit successfully while individual molecules have no usable
coordinates, including with `--gen3d fast`. The command exit code
and the message "molecules converted" do not certify coordinate generation.

Coordinate failures are recorded and skipped by default. If the user later
requests regeneration, start from the original declared SMILES and
verify its specified stereochemistry before reuse. Re-embedding the collapsed
SDF can lose double-bond stereo or assign different tetrahedral stereocenters.
Do not silently drop failed molecules or change the molecule under its old ID.

Names remain supplied by Open Babel; normalization does not assign IDs, require
manifests, enumerate states or create state seals. There is no separate
conversion approval step. Direct SDF inputs also receive the coordinate check,
and docking preparation repeats it before invoking receptor tools or an engine.
Property filtering and docking preparation assign stable IDs when ligand IDs
are missing, preserving the normalized source file and explicit identities.

For a retrospective benchmark, prepare verified `ACT_<original_id>` and
`DEC_<original_id>` names **before** normalization and preserve the source/ID
mapping; see `../pipeline-entry/references/benchmark-labels.md`. Numeric IDs do
not imply activity. The adapter preserves these names but never invents labels.
If either labeled class is unavailable, disable `coordination.screening_benchmark`.

## Verification

Library outputs are `library/library_3d.sdf`, `library/library_unfiltered.sdf`
and `normalization_report.json`, plus `config_patch.json` only on success.
`normalization_progress.json` records live conversion progress; use the final
report and `success` value to decide whether inputs are ready for the pipeline.
Individual molecule failures produce warnings and skipped-record entries while
the usable library continues. An all-failed retry removes a previous config
patch. Valid SDF records are not rewritten.
Library coordinate-generation commands directly invoke the selected native or
compatibility executable with `--gen3d fast`, without a per-record Conda launcher.
`library.conversion_execution.command_backend` is `openbabel_native_fast` or
`openbabel_fast_writeback_compat`. Probe files stay in `runtime/`; the compatibility
path also keeps `runtime/openbabel_fast_compat` and `runtime/fast_compat_build.json`.
The first compatible build is cached in that output directory and reused while
its source and dependencies match. Timeout entries identify
the skipped molecule as `conversion_timeout`. Temporary per-record directories
are removed after conversion; the merged SDF files and diagnostics remain in
`inputs/normalization/`.

Receptor conversion checks atom count, order, finite coordinates and drift
(default tolerance 0.002 Å). MOL2 bond types and partial charges cannot all be
represented in PDB; the report records this limitation and the source checksum.
Prepared receptors retain source/water coordinates and receive a separate noH
copy and prep manifest. Final `docking.docking_prep` still applies; see
`../pipeline-entry/references/preparation-contract.md`.

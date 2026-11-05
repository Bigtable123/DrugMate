# Evaluated research memory

Use this protocol for ordinary memory-enabled research. Frozen benchmark runs
and memory-off ablations must neither read nor write this store. Record the scope
(`research` or `frozen`) in `memory_consultation.md`; it is not an unregistered
pipeline config field. The `--mode frozen` commands return disabled without
opening query inputs or the store.

The store lives under `SBDD_DATA_ROOT/research_memory`, outside the source checkout
and frozen bundles. Each immutable observation retains the actual prepared
inputs, encoder identity, tool settings, evaluated measurements and hashed source
references. Repeated completed experiments append observations to the same
prepared-input profile; exact re-imports are idempotent. Different prepared input
bytes or encoder versions form distinct profiles, even if the chemical name is
the same. This avoids deduplicating experiments by approximate embedding equality.

## Append after evaluation

After a run reaches a terminal state and its evaluation artifacts are written,
create a small import record under that run's output directory and call:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" append-research-run \
  --mode research --record "${SBDD_DATA_ROOT}/REQUEST/round1/memory_observation.json"
```

The record has `schema_version=1`, the absolute round directory as `run_id`,
`case_id`, `group_id`, `task` (`docking`, `pocket_prediction`, or
`molecule_generation`), and these fields:

| Field | Actual source |
| --- | --- |
| `inputs` | `receptor_file`, compatible `query_npz`, bound `query_metadata`; also `ligand_sdf` for docking/generation |
| `execution` | A `{path, sha256, pointer}` reference to `/research_provenance` in the relevant stage's `result.json`; a list permits linked preparation/evaluation receipts from the same round |
| `evidence.status` | A reference to the terminal `/status` in that round's `manifest.json` |
| `evidence.tool` | A reference to the actual tool ID or a nonempty list of distinct tool IDs in a generated report |
| `evidence.parameters` | A reference to the actual parameter object in that report, e.g. `/parameters` |
| `measurements` | One or more references to bounded named measurement blocks, with original values, units, availability and failure status retained |

Every reference is an object with exactly `path`, `sha256`, and `pointer`.
Use the file's actual SHA-256 and an RFC 6901 JSON pointer. Do not rewrite a report,
insert expected scores, or manufacture a receipt to make ingestion succeed.
The runner captures input signatures before execution and the JSON reports
produced by that execution in `research_provenance`. The importer verifies these
bindings and rejects evidence from another run or different inputs. Cached results
retain their original execution identity and are not new experiments.

Use the reports produced by the actual stage; their schemas are different:

| Task | Tool reference | Parameters and measured evidence |
| --- | --- | --- |
| Docking | `docking_report.json` → `/engine` | Report `/parameters` and the evaluated molecule's own report/ledger row |
| Pocket prediction | `pocket_analysis_report.json` → `/pocket_detection/consensus_pocket/active_tools` | Report `/parameters`, `/validation` and bounded detector-result blocks; absent reference DCA remains unavailable |
| Molecule generation | `source_cohort_report.json` → `/cohorts/0/source_tools` for the actual cohort | `library_report.json` or stage `result.json` → `/parameters`; cohort report `/sources` retains per-source results and terminal states |

For multiple generation cohorts, import each relevant cohort with its own tool
reference and evaluated evidence. A tools list records attempted sources, not
proof that every source succeeded. Do not substitute an invented common `/engine`
field for the actual pocket or generation report structure. A run may be terminal
as `failed_<step_id>`; retain the failed measurements and reasons without turning
them into successful labels. Intermediate paused stages are not terminal runs.

Use linked receipts when an evaluation stage consumes prepared inputs from an
earlier stage. An encoded hydrogen-free receptor and a protonated docking receptor
are distinct files; retain their real preparation lineage instead of assigning
the same hash. If old runs lack execution receipts or the current inputs have no
compatible encoding, keep their reports and record ingestion as unavailable.
Do not rerun scientific tools merely to make a memory record.

### Evaluated molecules from a docking library

For a multi-molecule docking input, append one profile observation for each
evaluated molecule with compatible query features. Extract its original SDF
record byte for byte, including the `$$$$` delimiter, into this round's inputs;
set `inputs.ligand_sdf` to that single-record file. Keep the original source
coordinates and chemistry; a newly selected best pose is not a replacement
retrieval input. Protein encoding can be reused only with its verified receptor,
task and encoder binding.

Add the optional top-level `ligand_source` object to the import record:

```json
{
  "ligand_source": {
    "path": "<captured input library's absolute path>",
    "sha256": "<SHA-256 of that library>",
    "record_index": 0,
    "identity": {
      "path": "<round docking directory>/docking_prep/microstate_identity_map.json",
      "sha256": "<SHA-256 of the original preparation map>",
      "pointer": "/records/0"
    }
  }
}
```

Replace `0` with the molecule's zero-based source-library index in both fields.
The importer verifies that the library and preparation map belong to the same
docking execution and that the map's `preserve_declared` policy retains the
source identities. Each selected measurement must retain that record's
`compound_id` and, when present, its `source_molecule_index`. Select its row in
the actual docking ledger's `/molecules` or evaluation report's `/evaluation`;
find the row by identity, since report ordering need not equal source ordering.
Aggregate library scores and other compounds' rows cannot become that molecule's
observations. Include the evaluation execution receipt when using its report.

The stored observation and later recall retain this source-record binding.
Matching prepared input and encoder identities append to the existing profile;
new identities create profiles. Re-importing a cached execution does not append
another experiment, even if unrelated manifest metadata changed. Different real
executions remain separate observations.

Only real, terminal observations may be appended. Explicit failed or partial
evaluations remain diagnostic evidence, not successful outcomes. A best-pose
redocking RMSD is not native Top-1 RMSD; a single generator run is not a winner
over untested generators; a predicted site without a reference does not have
measured DCA. Missing measurements stay missing.

## Recall before the next plan

When compatible current query inputs are ready, retrieve research evidence
before finalizing the next tool matrix. For docking, for example:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" recall-research \
  --mode research --task docking \
  --receptor-file "${SBDD_DATA_ROOT}/REQUEST/inputs/receptor.pdb" \
  --ligand-sdf "${SBDD_DATA_ROOT}/REQUEST/inputs/ligand.sdf" \
  --query-npz "${SBDD_DATA_ROOT}/REQUEST/inputs/embeddings.npz" \
  --query-metadata "${SBDD_DATA_ROOT}/REQUEST/inputs/query_metadata.json" \
  --exclude-run "${SBDD_DATA_ROOT}/REQUEST/round2" \
  > "${SBDD_DATA_ROOT}/REQUEST/round2/inputs/research_memory_recall.json"
```

Use the actual allocated paths. Pocket queries omit `--ligand-sdf` and use
`--task pocket_prediction`; generation uses `--task molecule_generation` with
its conditioning ligand. The default retrieval is three ligand-similar profiles
with protein reweighting for docking/generation, or five distinct protein groups
for pocket evidence. An explicit `-k` is a diagnostic override, not a change to
the fixed selector settings.

Prior completed observations of the same profile are useful in research mode;
exclude the current run and any disallowed cases/groups explicitly. Frozen
evaluation keeps its separate case/group exclusion and no-write policy. Never
use held-out answers or another ablation arm's observations as research memory.

Read the returned protocol, evidence and failure records before adopting advice.
Record whether a tool choice follows the frozen selector, user requirements or
new research evidence. Preserve the original numerical recommendation and the
reason for any planner override. Do not average different protocols or allow
many repeats of one profile to become independent nearest neighbors.

Research recall returns similarities and evaluated observations for the planner.
It does not silently modify the frozen selectors' numerical posterior estimates
or generate new cross-tool winner labels. Promoting observations to a new
benchmark snapshot remains a separate comparative evaluation and index-building
operation with an explicit new version and data split.

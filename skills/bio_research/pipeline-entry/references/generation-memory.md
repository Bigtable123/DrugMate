# Memory for Molecule Generation and Design

This is the generation branch of the [unified memory entry point](benchmark-memory.md).
Before pipeline execution, `ligand_protein_generation_reweight_v1` uses one valid
reference ligand and the current protein encoding to recommend tools separately
for five metrics. It neither runs generation nor edits the config, and does not
replace the chemical-modification contract of `molecule-modifier`.

## Queries and decisions

1. Confirm that the current task requires generation. If the user supplies a
   complete library, requests only known-active retrieval, or explicitly excludes
   generation, record this branch as `not_applicable`. Memory advice must not
   replace the user's task.
2. Prepare one clearly valid reference-ligand SDF, the current receptor,
   compatible `embeddings.npz`, and binding metadata. Without a reference ligand,
   the ligand-similarity policy cannot run. Record the limitation and plan from
   current facts and the catalog; do not substitute the candidate library's first
   molecule or a hidden evaluation answer.
3. Check generation tools' current inputs, runtimes, models, and resources.
   Explicitly pass the executable set through `--available-tools`. Registry
   registration only establishes that an adapter exists, not that the current
   environment can run it. DiffLinker requires valid user-supplied fragments;
   ResGen and MolCRAFT require a reference ligand.
4. Invoke the selector and record neighbors, weights, per-metric scores,
   recommendations, actual tools, and override reasons. Without a priority
   metric, use `default_recommendation` (the property-filter-pass-rate
   recommendation) as the minimal generation route. For one explicitly specified
   metric, adopt its recommendation. If multiple metrics disagree, the Agent
   selects a route according to authorized priorities and the total molecule
   budget, or divides the budget for independent comparisons. Do not combine
   metrics into one score or automatically enable every recommended tool.
5. Apply the actual choice to the complete explicit tool matrix and revalidate
   the config. Retain existing post-generation quality checks, including property
   assessment and PoseBusters/PoseCheck. Historical winner scores do not establish
   activity.
6. In ordinary memory-enabled research, append the completed run's verified
   observations with `append-research-run` and recall compatible prior runs with
   `recall-research --task molecule_generation` before the next plan, following
   the [research-memory protocol](research-memory.md). Record any planner override
   separately from the frozen selector's recommendation. Benchmark evaluations
   and memory-off ablations skip this research store entirely.

Overall `success/status` checks the user's priority metric, or
`gate_pass_rate_unique` if none is specified. Other metrics still return their
own `selected/fallback/unavailable` status. An unavailable supplementary
recommendation does not block a satisfied primary objective. Before adoption,
check the selected metric's tool and input requirements rather than reading only
the overall status.

Use the input directory already prepared for this request, for example:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" select-generation \
  --selector-index "${SBDD_MEMORY_ROOT}/generation_memory/ligand_protein_v1.json" \
  --ligand-sdf "${SBDD_DATA_ROOT}/current_request/inputs/reference_ligand.sdf" \
  --query-npz "${SBDD_DATA_ROOT}/current_request/inputs/embeddings.npz" \
  --query-metadata "${SBDD_DATA_ROOT}/current_request/inputs/query.json" \
  --receptor-file "${SBDD_DATA_ROOT}/current_request/inputs/receptor.pdb" \
  --available-tools pocket2mol tamgen diffsbdd \
  > "${SBDD_DATA_ROOT}/current_request/inputs/generation_memory_decision.json"
```

`select-generation` does not accept `--config`. To prioritize metrics, repeat
options such as `--priority-metric median_qed` or
`--priority-metric scaffold_entropy_n50`. Using the default profile does not mean
that all five metrics are user priorities.

Existing frozen inputs may use the generation release's `query.json` format:
`task_type` is `molecule_generation`, with `identity`, `files`, and `file_hashes`
binding the NPZ, reference ligand, receptor, and `embedding_metadata.json` to the
current inputs. For new encodings, use `bind-query-metadata` below to verify
original completion records and produce binding metadata. A filename, a manually
added current hash, or replacing an old encoding's hash with a new receptor's
hash is not a valid binding. The encoding model, dtype, windowing, pooling, and
pLDDT mask must be compatible with the index. Queries read only protein vectors,
not pocket features.

If no encoding exists, configure compatible local model and Foldseek paths
according to the [encoder guide](../../../../memory/saprot_embedding/README.md),
then run `run_batch.py` in a compatible environment. Create a manifest containing
only the current query, for example `inputs/generation_embedding_manifest.csv`:

```bash
python3 - <<'PY_MANIFEST'
import csv
import os
from pathlib import Path
inputs = Path(os.environ["SBDD_DATA_ROOT"]) / "current_request/inputs"
with (inputs / "generation_embedding_manifest.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["task_type", "dataset", "target", "protein_path", "plddt_mask", "plddt_threshold"])
    writer.writerow(["molecule_generation", "current_request", "target1", str(inputs / "receptor.pdb"), "False", ""])
PY_MANIFEST
```

```bash
python "${SBDD_MEMORY_ROOT}/saprot_embedding/run_batch.py" \
  --manifest "${SBDD_DATA_ROOT}/current_request/inputs/generation_embedding_manifest.csv" \
  --output "${SBDD_DATA_ROOT}/current_request/inputs/encoding" \
  --task-type all --dtype float32 --storage-dtype float16
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" bind-query-metadata \
  --encoding-dir "${SBDD_DATA_ROOT}/current_request/inputs/encoding/molecule_generation/current_request/target1" \
  --receptor-file "${SBDD_DATA_ROOT}/current_request/inputs/receptor.pdb" \
  --output "${SBDD_DATA_ROOT}/current_request/inputs/query_metadata.json"
```

Here, `--task-type all` must receive only that one-row manifest, never the default
full-library manifest. The actual task remains `molecule_generation`; do not
change it to docking or fabricate a cohort. This task generates a whole-protein
encoding and requires no reference box. Pass `embeddings.npz` from that same
encoding directory and the generated `query_metadata.json` to
`select-generation`. Preserve original provenance records such as `prepared.json`
and `complete.json`. If binding fails, retain the error and record the fallback;
do not bypass provenance checks. This input encoding does not require running
any generation or docking tool.

When query case/group identity is known, use `--query-case dataset/target` and
`--query-group ID`; `--exclude-case` and `--exclude-group` may also be repeated.
The case in metadata is excluded automatically. If it belongs to memory, its
entire protein group is also excluded. Refit descriptor scales and priors after
exclusion. Do not claim leave-protein-group-out evaluation when the group is
unknown.

## Fixed five-metric policy

Ligands use a custom radius 0–3, 2048-bit circular graph fingerprint and an
18-dimensional descriptor vector, **not Morgan**. Ligand similarity
`0.8 × Tanimoto + 0.2 × descriptor RBF` retrieves the top-3 neighbors. The
normalized 1280-dimensional protein vector provides only percentile-based
multiplicative reweighting, with `beta=0.10`. Neighbor weights are not normalized.
Historical labels identify each case's winning tool per metric (all tied winners
receive 1); they are neither docking RMSD labels nor individual-molecule
success/failure labels.

All five metrics use `prior_strength=10.0` and `switch_margin=0.03`.
`--priority-metric` chooses which objective's
recommendation the Agent should use and which recommendations must be available
for overall success. It does not change retrieval, historical scores, priors,
or switching thresholds. The frozen feature and label indexes are unchanged.

Compute prior-weighted winner scores separately for each metric. Switch from
Pocket2Mol only when the advantage is **strictly greater than 0.03**; ties retain
Pocket2Mol. Scores are uncalibrated, and five recommendations
do not imply that one tool achieves all five outcomes simultaneously. Joint
QED–SA pass rate is not used. Property pass rate is the target-level proportion
`(KEEP + MODIFY) / valid_unique_candidates` among valid unique molecules.

For invalid inputs/indexes, fewer than three neighbors after exclusion, or an
unavailable recommended tool, emit the specific fallback reason and use
Pocket2Mol if available. If it is also unavailable, return `unavailable` rather
than claiming an executable plan. Record the original `proposed_tool` separately
from the actual `selected_tool`.

## Offline index and data isolation

Queries use the supplied `memory/generation_memory/ligand_protein_v1.json` index
by default. If the frozen release changes, build a separate index offline on the
administrator side, never during formal ablation queries:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" build-generation-index \
  --benchmark-bundle "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-generation-memory/memory" \
  --property-winners "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-generation-memory/evaluation_v2/target_best_tools.csv" \
  --groups "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-generation-memory/audit/similarity_groups.csv" \
  --output "${SBDD_DATA_ROOT}/generation_memory/ligand_protein_v1.json"
```

First extract `drugmate-generation-memory.zip` and set
`SBDD_MEMORY_BUNDLES_ROOT` as described in the
[unified memory protocol](benchmark-memory.md#downloads-and-deployment-location).
`--benchmark-bundle` points to `memory/` inside the extracted directory, not the
extracted directory itself. The two CSVs are in its sibling `evaluation_v2/` and
`audit/` directories. `--output` must be outside the frozen package. Keep the
original supplied index rather than overwriting it with a rebuild. These argument
paths are not pipeline config fields.

The release identity is `agent_memory_generation_v2_20261018`. Configure
`SBDD_MEMORY_BUNDLES_ROOT` for the source bundle. Preserve its input hashes and
record the version of any separately rebuilt index.

The incoming download already lacked 49 memory and 12 validation `cavity.as`
files. The current generation index builder and bundled reader do not read them,
but cavity-dependent preprocessing or full original-file validation still requires
recovering them from the original source. Do not alter manifests or source hashes,
or skip inputs actually required at runtime, because of this observation.

Frozen memory cards already contain winners for the first four metrics. The fifth
metric and protein groups come from administrator-side CSVs. The builder creates
an allowlist from memory-card case IDs before accepting matching fields with
`split=memory`, validating input hashes and recording source versions. It does
not retain validation/final labels. A few historical multi-record ligands retain
the released study's record 0 feature convention. V3000 ligands retain the study's
dative bond → order 1 graph-feature convention. These are historical compatibility
rules only: they neither change the online single-reference-ligand contract nor
rewrite pipeline molecules.

`select-generation` reads only the frozen derived index and current valid inputs. Do not
substitute the old `memory/read_memory.py recommend` command, which uses different
multimodal weights, or invoke parameter-search or whole-split recomputation
scripts. Formal evaluations isolate administrator-side evaluation/audit data,
other arms' outputs, and hidden labels. Keep source observations in this round's
reports; research-mode ingestion copies verified evidence into the separate
store rather than writing it back to frozen memory. Record input/index
versions, query status, exclusion rules, recommendations, and actual decisions in
`memory_consultation.md`, and present important choices at Gate A.

This is the evaluation-snapshot policy. For user research, follow the separate
[append and recall protocol](research-memory.md): retain real evaluated runs,
including repeats, with their input/protocol/encoder identities and evidence.
Retrieved observations guide the planner; they do not automatically replace the
frozen winner labels. Do not treat one generator run as a winner over untested
tools. Missing comparison labels remain unknown. Offline rebuilding of the
supplied snapshot alone does not make new user results available to retrieval.

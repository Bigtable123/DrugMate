# Docking Memory: Redocking, Pose Generation, and Independent Objectives

Follow the [unified memory protocol](benchmark-memory.md), then complete the
docking branch described here. Memory provides planning advice only. The Agent
still explicitly configures the actual engine according to inputs, availability,
and the user's requirements.

## Task boundaries

The redocking selector is `top1_rmsd_lt_2a`, with GNINA as the default,
prior strength 2 and a strict switching margin of 0.10. Historical labels refer
to the released observation's published Top-1 pose, not the best RMSD among
several generated poses. Preserve its published ranking; do not claim a native
rank where the original memory card records none.
The runtime's best-pose redocking calibration gate is a separate protocol check;
its pass/fail result cannot be imported as a native-top-1 memory label.

Select `--task pose_generation` for the **rDock** baseline. This route
uses the same ligand retrieval and protein reweighting equations and defaults to
published Top-1 `posebusters_valid` labels from released docking observations. The
choice of that objective and reuse of prior strength 2 / switching margin 0.10
are implementation assumptions recorded with the selection. This policy does
not establish experimentally measured performance for a new query.

The output records `task`, `default_tool`, `metric`, and `task_policy`, including
these assumptions. The pose route reports `policy_version=pose_generation_inferred_v1`;
its `index_policy_version` retains the frozen index's version. Other explicitly
selected metrics are marked `objective_scope=explicit_extension`. In particular,
the continuous RMSD loss cap and its switching rule are additional policies, not
the default binary reliability equations. Follow the selected engine's real
cavity, receptor and ligand contracts, and enable only one primary docking engine.

## Modes and minimal commands

Set the mode in this round's existing pipeline config (this is only a field
fragment):

```json
{"context":{"docking_memory_mode":"ligand_protein"}}
```

| Mode | Docking planning behavior |
|---|---|
| `ligand_protein` (default) | Top-3 ligand retrieval, protein-similarity reweighting, and deterministic tool scoring |
| `ligand_only` | The same ligand retrieval and scoring, with a protein-reweighting coefficient of 0 |
| `legacy` | Original protein-cosine neighbors and lightweight cards; supports `--task redocking` only |
| `off` | Reads no docking query features, index, or historical library; returns the selected task's GNINA or rDock baseline |

After setting the config mode to `off`, run:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" select-docking --task redocking --config "${SBDD_DATA_ROOT}/current_request/config.json"
```

For a modern enabled mode, run after preparing inputs and before docking,
replacing placeholders with actual paths:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" select-docking --task redocking --selector-index "${SBDD_MEMORY_ROOT}/docking_memory/ligand_protein_v2.json" --config "${SBDD_DATA_ROOT}/current_request/config.json" --ligand-sdf "${SBDD_DATA_ROOT}/current_request/query_ligand.sdf" --query-npz "${SBDD_DATA_ROOT}/current_request/embeddings.npz" --query-metadata "${SBDD_DATA_ROOT}/current_request/metadata.json" --receptor-file "${SBDD_DATA_ROOT}/current_request/receptor.pdb" --available-tools gnina autodock_vina
```

For memory-guided pose generation, use the same verified inputs with the explicit
task, for example:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" select-docking \
  --task pose_generation \
  --config "${SBDD_DATA_ROOT}/REQUEST/config.json" \
  --selector-index "${SBDD_MEMORY_ROOT}/docking_memory/ligand_protein_v2.json" \
  --ligand-sdf "${SBDD_DATA_ROOT}/REQUEST/inputs/ligand.sdf" \
  --query-npz "${SBDD_DATA_ROOT}/REQUEST/inputs/embeddings.npz" \
  --query-metadata "${SBDD_DATA_ROOT}/REQUEST/inputs/query_metadata.json" \
  --receptor-file "${SBDD_DATA_ROOT}/REQUEST/inputs/receptor.pdb" \
  --available-tools rdock gnina autodock_vina
```

To inspect its memory-off baseline without reading those inputs:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" select-docking \
  --task pose_generation --mode off --available-tools rdock
```

Candidate tools in this example must already meet the current runtime, input
protocol, and resource requirements. Adjust `--available-tools` to the actually
available set. Omitting it means only that candidates are registered, not that
they are installed or satisfy this round's input requirements. `ligand_only`
requires no protein NPZ/metadata. If metadata is supplied, only case/group
identity is extracted to keep both modes' exclusions consistent; the NPZ is not
read and protein encoding is not checked. Standalone inspection may use
`--mode off|legacy|ligand_only|ligand_protein`; formal planning reads the mode from
`--config`. Supplying both with conflicting values is an error. The default is
`ligand_protein`; explicit null and unknown modes are invalid.

The mode is an Agent planning switch. The command only outputs a decision; it
does not edit the config or execute docking. Changing modes does not cause the
runner to query memory or replace an engine. The Agent uses the result and current
constraints to update the complete explicit tool matrix and `search_space_source`,
runs `validate_config.py`, and presents the final plan at Gate A. A user-specified
engine takes precedence; record why the actual choice overrides a recommendation.
`--task` is a selector argument recorded in its decision, not a new pipeline
context field.

## Choose one objective in advance

`--task redocking` is the default and resolves an omitted `--metric` to
`top1_rmsd_lt_2a`; `--task pose_generation` resolves it to `posebusters_valid`.
`--metric <metric>` overrides that objective without changing the task's baseline.
Fix the task and metric before running. Do not choose them using query ground truth or interpret independent
recommendations as simultaneously achieving all five best outcomes in one run.

| `--metric` | Meaning | Selection direction |
|---|---|---|
| `top1_rmsd_lt_2a` | Rank-1 pose RMSD strictly below 2 Å | Higher is better |
| `posebusters_valid` | Rank-1 pose passes the released version's PoseBusters checks | Higher is better |
| `joint_success` | Both Top-1 RMSD <2 Å and PB-valid | Higher is better |
| `top1_rmsd_a` | Continuous RMSD of the rank-1 pose | Lower is better |
| `median_rmsd_a` | Median RMSD over all valid poses within each case/tool | Lower is better |

For the two continuous RMSD objectives, priors and scoring use only
`min(value, 10 Å)` for historical RMSD. The index and neighbor audit retain raw
values. A within-case median is neither the full dataset's median Top-1 RMSD nor
the best pose. PB-valid and geometric reproduction measure different properties;
neither substitutes for the other.

## Modern read-only index

Queries use the supplied `memory/docking_memory/ligand_protein_v2.json` index
by default. Rebuild a separate copy only when explicitly updating the memory
release; do not rescan experiments for every query:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" build-docking-index \
  --benchmark-bundle "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-docking-memory" \
  --output "${SBDD_DATA_ROOT}/docking_memory/ligand_protein_v2.json"
```

First extract `drugmate-docking-memory.zip` and set
`SBDD_MEMORY_BUNDLES_ROOT` according to the
[unified memory protocol](benchmark-memory.md#downloads-and-deployment-location).
`--benchmark-bundle` points to the extracted directory itself, which contains the
historical records at `memory/agent_memory.jsonl`. Do not append an extra
`memory/` layer to the bundle path. Write outputs outside the frozen package and
do not overwrite the read-only index distributed with the source.

The builder extracts only case/group IDs, valid input-ligand features, protein
vectors, and published labels for the five objectives from the frozen memory
partition. It caches feature version, provenance, scales, and priors. It does not
load final labels, search parameters, rerun scientific tools, or modify the
original benchmark. To use another index, pass `--selector-index <index.json>` to
`select-docking` and hold that version fixed throughout the ablation.

Online selection reads only the cache and current valid inputs. A query accepts
exactly one ligand. Features use the study's custom radius 0–3, 2048-bit graph
fingerprint and 18-dimensional descriptors, not RDKit Morgan. Geometric
descriptors depend on the input conformation, which must be prepared consistently
with memory and independently of the true bound pose. Do not read the current
query's RMSD or true pose. `site.json` is not used for tool scoring. Protein mode
uses compatible 1280-dimensional `protein_embedding_normalized`, with metadata
for encoder and provenance checks. `--receptor-file` may verify receptor binding.
For newly generated raw SaProt outputs, first call `bind-query-metadata` as
specified in [Binding current encoding provenance](benchmark-memory.md#binding-current-encoding-provenance),
then pass the separate binding file as `--query-metadata`. Do not manually insert
hashes into or overwrite original encoder metadata. Encoding completion records
and the preparation recipe provide the binding evidence.

Both selection tasks consume the existing **docking encoder protocol**. Bound
metadata keeps `source_record.task_type="docking"`; the decision's separate
`task` selects the planning objective and baseline. Do not relabel encoder
metadata as `pose_generation` to make a query appear compatible.

The fixed policy is `ligand_protein_reweight_v1`: ligand similarity is 0.8
fingerprint Tanimoto + 0.2 standardized-descriptor kernel, retrieving the top-3.
Weight by ligand rank using `1/log2(rank+2)`; protein mode further multiplies by
`exp(0.10 * (protein_percentile - 0.5))`. Weights are not normalized to sum to 1.
For binary metrics, smooth that metric's labels with a tool-global success-rate
prior of strength 2.0. Switch tools only when the success score exceeds the
task's baseline score by strictly more than 0.10. For continuous RMSD, use a global mean-clipped-loss
prior of the same strength, switching only when estimated loss is strictly more
than 0.10 Å below the baseline's. Both types retain the task's baseline on ties.
The baseline is GNINA for redocking and rDock for pose generation; the latter's
guard is the explicit inference described above. `ligand_only` removes
only the multiplicative protein correction. Output the `metric`, direction,
units, tool scores, and neighbor weights for inspection. Binary scores are not
probability-calibrated.

An old `ligand_protein_v1.json` may still be passed explicitly, but supports only
the original Top-1 <2 Å objective. Choosing another objective with an old index
returns `unsupported_metric_in_legacy_index`; do not claim that retrieval for
that objective succeeded. Rebuild a v2 index as above. The default-objective
algorithm is unchanged between old and new indexes. Historical AutoDock-GPU has
no current adapter, and current ADFR has no labels in this release. Neither may
borrow another engine's results.

For invalid inputs/indexes, fewer than three valid neighbors, or an unavailable
recommended tool, record the specific reason and fall back to the task's baseline
if available. If that baseline is unavailable, return `unavailable`; a pose task
must not silently adopt GNINA as its fallback. Missing labels
must not be treated as failures. If a valid single-ligand query cannot be obtained,
or the task involves a multi-molecule screening library, record the applicability
limitation and plan according to the user's objective and current tool contracts.
Do not use the library's first record to represent the entire library or describe
a single-ligand model as a library-level policy.

## Ablation and data boundaries

In four-arm experiments, hold the task, metric, query, memory/index version, input preparation,
available tool set, seed, protocol, and evaluation definitions fixed; change only
the memory mode. Use isolated sessions and result directories for each arm. An
`off` session must not first read another arm's docking history cards, tool
results, recommendations, or consultation records, and must not add docking
recall afterward. Keep pocket and molecular-design workflows identical across
all four arms according to the unified protocol. This switch does not claim to
disable all Agent memory. Maintain the declared mode across rounds rather than
recovering docking history from an old session.

When query identity is known, pass `--query-case dataset/target`; when the protein
group is known, pass `--query-group <group_id>`. `--exclude-case` /
`--exclude-group` may be repeated. Evaluation within memory must exclude the
query case and its protein group. Do not claim leave-group-out evaluation when
the group is unknown. Descriptor scales and priors for formal cross-validation
must also be fitted only on training candidates remaining after exclusion.
Formal split validation still requires restricting the Agent's file access to
final labels and other arms' results.

Keep new results in this round's reports and artifact lineage. For ordinary
memory-enabled research, use the [append and recall protocol](research-memory.md)
to make evaluated observations available to later plans without changing the
frozen selector index. Build new supervised snapshots separately after obtaining
comparable evaluation labels; predicted scores cannot replace RMSD labels.

In `memory_consultation.md`, record the task, objective metric, mode, baseline, task-policy assumptions, index/feature
versions, query provenance, exclusions, neighbor similarities/weights, tool
scores, recommendations and actual tools, and fallback/override reasons. `off`
records only its disabled state and actual baseline, without loading historical
material to fill out the record. Raw JSON outputs may remain in the same request
directory. Include important choices in the Gate A plan.

## Original legacy inspection workflow

The following applies to `docking_memory_mode=legacy`. The original pocket-task
evidence reader may still be queried when needed, but does not replace the unified
protocol's pocket policy. Legacy mode accepts the query arguments below through
`select-docking --config <config.json>`, internally reusing original recall.
Direct `recall --task docking` is only for explicitly chosen legacy inspection
or reader checks and must not be added to modern or `off` experiments. Legacy
returns evidence cards, not deterministic tool recommendations.

The reader queries the fixed read-only bundle at
`${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-docking-memory/`. Configure that bundle
root before using legacy queries. Its release identity is
`agent_memory_split_v1_20261012`. The reader reuses the
bundle's `DecisionMemory` for source-hash verification and exact cosine retrieval.
Lightweight `planning` cards are derived views of original evidence, bound to
original-card context/revision, source files, and supplementary timing records.
Original cards, vectors, raw evidence, and the lightweight planning library remain
frozen.

This download retains 487 docking and 407 pocket_prediction contexts. `assets/`,
`derived/log_excerpts/`, `historical_references/`, databases, and timing
supplements are evidence sources for the original reader. Do not remove them as
ordinary logs; pocket selection also needs historical DCA and protein evidence. The
package does not include a generated `planning/` sidecar; building that view still
requires the existing reader and its data.

New project observations remain in this round's reports and artifact lineage.
In research mode, verified observations can also be appended to the separate
[research memory](research-memory.md) and recalled before later plans. Neither
that store nor legacy queries automatically rewrites frozen selector scores or
the pipeline config.

### Submit a legacy query

Direct `recall` requires `--task`, either `docking` or `pocket_prediction`.
`select-docking --mode legacy` accepts only `--task redocking` (the default),
internally querying docking evidence. It does not support `--task pose_generation`.
Choose one query method:

- Current inputs: supply both `--query-npz` and `--query-metadata`.
  `--receptor-file` can verify the encoding's content binding to the current
  source receptor.
- Bundled case: supply both `--dataset` and `--target`, for reader checks or
  queries of already encoded cases.

Replace the following paths with actual inputs:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" recall --task docking --query-npz "${SBDD_DATA_ROOT}/current_request/embeddings.npz" --query-metadata "${SBDD_DATA_ROOT}/current_request/metadata.json" --receptor-file "${SBDD_DATA_ROOT}/current_request/receptor.pdb"
```

Example of a bundled-case check:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" recall --task docking --dataset astex_diverse --target 1G9V_RQ3 -k 3
```

Bundled queries and NPZ queries with explicit dataset/target identity in metadata
automatically exclude the current case. Repeat `--exclude-case dataset/target`
when needed. Cases in the same `group_id` are not independent experimental
replicates. Formal split validation also requires restricting the Agent's file
access so it does not read isolated new observations.

If no compatible current encoding exists, reuse
`${SBDD_MEMORY_ROOT}/saprot_embedding/run_batch.py`, configuring compatible local
model and Foldseek paths through the [encoder guide](../../../../memory/saprot_embedding/README.md).
Use a `--manifest` containing only current inputs and an `--output` for this
request. Match the library's model, `--dtype float32`, windowing, storage
precision, and explicit pLDDT mask. Docking encoding requires this request's valid
known box; unknown-site pocket prediction must not introduce a reference ligand
or reference box. Do not use the default full-library manifest or rerun the
entire library or experimental tools just to query history.

If the encoding environment is unavailable, `recall --task docking` (or
`pocket_prediction`) may still be called to record `not_requested` and the specific
reason in the consultation, then plan normally from current inputs. Missing
fields in a submitted query, incompatible encoding, or read failures return
`error`, `success=false`, and a nonzero exit code. Neither `not_requested` nor
`error` means "no similar cases." When metadata supplies NPZ or source-receptor
hashes, the reader verifies them.
`sources.benchmark.query.npz_hash_checked / receptor_hash_checked` indicates
whether binding verification was completed. If not, inspect the current encoding
records; identical filenames do not establish identical provenance.

### Assess whether legacy neighbors are usable

Retrieval ranks compatible encodings by cosine similarity and has **no minimum
similarity threshold**; low-similarity cards are still returned. A returned card
identifies a nearest neighbor within the library. Similarity is not a success
rate, a best-tool recommendation for the current protein, or a guarantee that
parameters transfer.

Both task types in the current release have candidate libraries with valid
compatible encodings. `no_matches` means only that no candidates remain after
exclusion. Failure to find a compatible encoder is an `error`, not an empty
retrieval result.

Before adoption, the Agent checks against current inputs:

- Whether structure, author chain, conformation, species/mutation/construct, and
  current pocket have a defensible basis for comparison.
- Site and coordinate frame, site/search strategy, ligand mechanism, protonation,
  and water/metal/cofactor conditions.
- Tool version, scoring and preparation protocol, valid input conditions, input
  cohort, and scope of quality comparisons.
- Whether original evidence supports the cited observation and the current
  catalog still supports the relevant tools and parameters.

Observations under unknown or different conditions form hypotheses to be tested,
not established conclusions. If none transfer, all cases may be excluded; plan
normally from current facts and the catalog and record exclusion reasons. Do
not directly copy old boxes, residue numbers, paths, checkpoints, or tool-switch
matrices. Historical observations do not automatically become
`known_required_interactions` or new gates. Artifact reuse follows the
[artifact-lineage contract](artifact-reuse.md); do not rank scores from different
protocols together.

### Read lightweight cards and details

The default `--benchmark-view planning` returns the complete lightweight `tools`
table, task conditions, valid input features, protocol differences, and detail
locations. Read each row according to `tool_columns`. `values` contains original
metric values in `metric_columns` order; `states` supplies
`[status, optional reason]` by column index, otherwise inheriting `states['*']`
or `common.metric_status`. Preserve `false`, `0`, `null`, and reasons for missing
values.

`run` orders execution status, parsing status, and return code. Null row values
for `run / flags / signals` inherit the corresponding `common` value; row
`protocol / cost` overrides the corresponding `common` field. Scope is
`extra.scope` or `common.scope`. A `seconds` key is an index into `time_columns`,
and its value is `[value, status, optional reason]`, in seconds. Historical
hardware was not controlled. Do not double-count combined timings and sub-stages;
a stage's elapsed time is not this round's total budget.

RMSD, DCA, and DCC are in Å. A Top-K/best-RMSD hit does not establish a PoseBusters
pass. `quality_cohorts` entries `[included, reason]` describe eligibility for the
source benchmark's quality comparison. Distinguish exclusion from a cohort,
absence of evaluation, execution failure, and scientific-metric failure.

The default returns at most 3 cards (`-k`), with a 48000-character cap on the
complete JSON. `--max-chars` may set 4000–60000. Tool rows are never truncated
individually. If output is too large, omit entire cards from lowest to highest
similarity. `sources.benchmark.retrieved / output_cards / omitted_cards`
distinguishes retrieved, emitted, and omitted counts. Empty `cards` with
`retrieved>0` or `omitted_cards>0` reflects output-budget omission, not
`no_matches`. When `summary_incomplete=true`, requery relevant cases with smaller
`-k` or larger `--max-chars`. Do not infer the absence of counterexamples from
omitted cards.

`sources.benchmark.planning_source=materialized` means existing lightweight cards
were read; `dynamic` means projection from verified original cards without
rerunning tools. Hash mismatch, stale provenance, or corruption in an existing
lightweight card is an error, not grounds for silent fallback.
`--benchmark-view summary` provides a more detailed summary with the same
retrieval ranking and exclusions.

For exact parameters, follow the card's `details.script / args` to call the
original bundle reader's `show --view execution --tool <tool_id>`. For raw
evidence, call `evidence` with the artifact ID and selector. Preserve
`external_only` / `source_unavailable` states rather than bypassing bundle
boundaries to fetch source material. Expand only needed tools and evidence; do
not put original cards, complete comparison outputs, or vectors into the
conversation.

### Record legacy tradeoffs

In this request directory's `memory_consultation.md`, record the queried task and
input provenance, query status, benchmark release, context/revision and
observation ID, similarities, evidence locations, applicability checks, and
reasons for adoption, further testing, or exclusion. Record errors and whole-card
omissions honestly. Explain important choices in the natural-language Gate A
plan, then assemble and validate the config against the current catalog. This
round's results, cross-round conflicts, and new-parameter checks reference the
project's reports and artifact lineage directly.

## Evaluation data boundaries

Each objective is selected independently. Oracle uses final ground truth and
must never be an online input. Online calls do not run evaluation or plotting
scripts or read final labels or per-query evaluation tables. Index builders
read only the frozen memory partition.

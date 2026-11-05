# Unified SBDD Memory: Pre-run Consultation Protocol

For pocket selection, docking and molecule generation/design, assess memory applicability separately
**before** finalizing the tool matrix and passing Gate A. Query only the branches
needed for this round. Pocket selection for an unknown site uses protein-memory
reliability fusion when compatible assets are available, with fixed DeepPocket#1
rank fusion only as an explicit unavailable-memory/off fallback. Record each domain's policy and
execution order together, without mixing their labels or scores.

## Entry point and responsibilities

The unified CLI is `pipeline-entry/scripts/benchmark_memory.py`, invoked through
`terminal` with `${SBDD_PYTHON}` from the configured scientific environment.
Frozen selectors perform read-only retrieval and deterministic inference;
`append-research-run` and `recall-research` manage separate evaluated research
observations under the [research-memory protocol](research-memory.md). The Agent
maps results to the current catalog's tools, parameters, and artifacts, then
validates them. The runner neither reads consultation records nor changes tools
automatically. Pocket Step 1 directly invokes its selector using registered
`pocket_memory_*` parameters. Do not invent unregistered `context.memory_*` fields.

The modern selector indexes and encoder entry points shipped with the source are
under the repository's `memory/`; see the [memory asset guide](../../../../memory/README.md).
The runtime uses the distributed `memory/` directory by default. Set
`SBDD_MEMORY_ROOT` to override it, or set it explicitly for these shell examples
after setting `BIO_SKILLS_ROOT`:

```bash
export SBDD_MEMORY_ROOT="${BIO_SKILLS_ROOT}/../../memory"
```

An explicit `--selector-index` takes precedence over `SBDD_MEMORY_ROOT`: use
`$SBDD_MEMORY_ROOT/docking_memory/ligand_protein_v2.json` for docking and
`$SBDD_MEMORY_ROOT/generation_memory/ligand_protein_v1.json` for generation.
The supplied indexes are read-only. Builders write separate derived indexes
under `SBDD_DATA_ROOT`.
Current query encodings, logs, and decision artifacts remain under `SBDD_DATA_ROOT`.
Legacy recall/planning cards still require the original reader and data in the
docking bundle described below.

### Downloads and deployment location

The repository already contains `memory/`. Extract both large bundles into its
`bundles/` subdirectory. From the repository root:

```bash
mkdir -p memory/bundles
unzip "${SBDD_DOWNLOADS_DIR:?Set the download directory}/drugmate-generation-memory.zip" -d memory/bundles
unzip "${SBDD_DOWNLOADS_DIR:?Set the download directory}/drugmate-docking-memory.zip" -d memory/bundles
export SBDD_MEMORY_BUNDLES_ROOT="$PWD/memory/bundles"
```

Each archive contains one outer directory matching its name. After extraction:

```text
memory/bundles/
  drugmate-generation-memory/
    memory/                  # generation --benchmark-bundle
    evaluation_v2/
    audit/
  drugmate-docking-memory/    # docking --benchmark-bundle
    memory/agent_memory.jsonl
    code/bundle_reader.py
    inputs/
```

`SBDD_MEMORY_BUNDLES_ROOT` defaults to `bundles/` under the resolved
`SBDD_MEMORY_ROOT`. Set it explicitly before launch to use external storage
containing the same two bundle directories. `SBDD_MEMORY_ROOT` selects the small
indexes and encoder scripts; `SBDD_MEMORY_BUNDLES_ROOT` selects the source
bundles.

`memory/bundles/` is excluded from Git and Python packages. Experiment outputs
still go under `SBDD_DATA_ROOT`. Online queries that use only the modern derived
indexes do not require these large downloads. Preserve frozen scientific
measurements, labels, input hashes and release identities when relocating
bundles. See the
[memory deployment guide](../../../../docs/memory-and-paths.md#download-and-placement)
for full paths and bundle scope.

| Domain | Applicable tasks and inputs | Selection entry point | Defaults and outputs |
|---|---|---|---|
| Redocking | One valid query ligand; protein mode also requires compatible encoding and provenance binding | `select-docking --task redocking`; see [docking memory](docking-memory.md) | Default: GNINA and published Top-1 RMSD <2 Å; one objective and engine per plan |
| Pose generation | The same docking input/encoder contract, with a pose-generation planning objective | `select-docking --task pose_generation` | Baseline: rDock; published Top-1 PoseBusters validity and alpha 2 / margin 0.10 are recorded implementation assumptions |
| Molecule generation/design | New candidates are needed, with one valid conditioning ligand and protein encoding | `select-generation`; see [generation memory](generation-memory.md) | Pocket2Mol baseline; separate recommendations for five metrics, without combining them into one score |
| Pocket generation/selection | Unknown site and current hydrogen-free receptor; memory needs a compatible protein encoding and separate pocket index | `plan-pocket-tools`, then main Step 1 `detector_consensus`; see [pocket memory](pocket-memory.md); `select-pocket` is for independent inspection | Five protein-group neighbors supply a repeated detector recommendation list from DeepPocket, fpocket, P2Rank, AutoSite and pyKVFinder. Execute its unique tools plus the required DeepPocket anchor; fuse using historical DCA reliability and the guarded P2Rank switch. Only unavailable/off memory uses fixed DeepPocket#1 four-detector rank fusion. |

Fix the metric, mode, and budget according to the user's objective before reading
the branch-specific policy. A user-specified engine, generator, site, or box takes
precedence over memory advice. A complete screening library, retrieval of known
actives only, or evaluation of existing poses does not justify additional molecule
generation. For purely de novo tasks without a valid conditioning ligand, mark
the generation similarity model `not_applicable` and plan under current tool
contracts. Generator selection does not replace `molecule-modifier` reaction or
scaffold constraints; DiffLinker still requires actual fragments/anchors.

Specify the docking task before querying: `redocking` is the default, while
`pose_generation` selects rDock as its guarded baseline. Omitted `--metric`
resolves to Top-1 RMSD <2 Å or published Top-1 PoseBusters validity, respectively.
Explicit metric overrides are supported extensions. Preserve returned
`task_policy` assumptions and use the [documented task scope](docking-memory.md#task-boundaries).
Both routes retain `task_type="docking"` in encoder metadata; a planning task
does not authorize rewriting provenance.

## Pre-run sequence

1. Initialize this round's directories through the entry point and normalize
   inputs. All query files, derived indexes, precomputations, logs, and consultation
   artifacts remain under the configured `${SBDD_DATA_ROOT}`. Frozen benchmarks
   serve only as read-only sources.
2. Create a three-row domain record in the returned
   `inputs_dir/memory_consultation.md`, stating the objective, mode/metric,
   required valid inputs, and current status. Before the config is complete,
   queries may use branch CLI arguments, but the final config and consultation
   modes must agree.
3. Decide the pocket strategy first, then resolve generation and docking
   selection inputs. Preserve an existing user box, reviewed `known_site`, or
   explicit `ligand_site`. For an unknown site, configure main Step 1
   `detector_consensus`. Run `plan-pocket-tools` with the bound receptor/query
   assets and pocket index before selecting the detector switches. Preserve its
   five-entry `recommended_multiset`, including repeats; explicitly enable the
   unique `execution_tools` once each, including the default DeepPocket anchor.
   The registered methods are DeepPocket, fpocket, P2Rank, AutoSite and pyKVFinder;
   enable AutoSite only when the plan selects it and its runtime is available.
   A selected detector's runtime failure requires an explicit repair or revised
   plan; it must not be mislabeled as unavailable memory to trigger the fallback.
   Step 1 performs native Top-3 point5 normalization and memory-guided fusion.
   Supply identical memory inputs through the registered `pocket_memory_*`
   parameters. Only memory-off/unavailable plans use all four fallback detectors
   (DeepPocket, fpocket, P2Rank and pyKVFinder). Prepare missing assets when possible; if memory is unavailable,
   record the reason before accepting `auto` fallback. `required` prevents fallback.
   Do not run extra generation or docking, or read the current query's ground
   truth, merely for pre-run consultation.
4. For pocket/generation/docking, reuse encodings only when provenance, task type, and
   encoder recipe all match. A shared 1280-dimensional shape does not justify
   cross-task reuse. If an encoding is missing, prepare only the current inputs;
   do not use a full-library manifest or rerun historical experiments. Ligand
   descriptor conformations must be independent of the true bound pose. Do not
   use the first molecule of a multi-molecule library to represent the entire
   library.
5. Once inputs are ready, invoke the relevant selector and save the raw JSON and
   source version. Check runtime, tool inputs, hardware, and budget. Include only
   tools that actually meet the requirements in `--available-tools`; registration
   does not establish executability. Build an index only on first setup or when
   explicitly updating the release. Do not invoke tuning, plotting, or evaluation
   scripts online.
   In ordinary memory-enabled research, also run `recall-research` for each
   applicable domain before choosing the final plan. Read its measured outcomes,
   actual protocols and failures; record whether this evidence changes the
   snapshot recommendation. An empty store is a recorded first-use condition,
   not a reason to omit future recall. Frozen evaluation and memory-off branches
   must bypass this research store.
6. Record differences between recommendations and actual choices. For generation
   without a priority metric, use the property-pass-rate branch as the primary
   selection. Use a user-specified metric when given. When recommendations
   conflict, state priorities or allocate multiple sources within the total
   budget; do not claim simultaneous optimality on all five metrics. Each docking
   config has exactly one primary engine. Recommendations for multiple objectives
   inform tradeoffs or separate experiments, not simultaneous activation of
   mutually exclusive engines.
7. Apply decisions to the complete explicit switch matrix, site/search strategy,
   parameters, and input artifacts. For an unknown site, use main Step 1
   `detector_consensus` plus `search_space_source=docking_box`, preserving the
   fused center and the selected anchor's point5 residues. Do not recompute residues through
   `provided_box`. Run `validate_config.py`, inspect the resolved plan, and check
   existing authorization scope at Gate A before executing the main experiment.
8. After evaluation reaches a terminal state, follow
   [research memory](research-memory.md) to append the actual evaluated
   observations when research memory is enabled. Record the returned profile
   and observation IDs, or the exact reason ingestion is unavailable. Preserve
   earlier observations; do not mark an append complete merely because the
   report was written. The next plan must query the updated store through step 5.

If current upstream artifacts are still needed, record `pending_inputs` and
perform necessary authorized preparation, but **do not start the downstream main
experiment with unresolved selection dependencies**. If required inputs truly
cannot be obtained, record why, explicitly choose an applicable baseline or user
plan, and revalidate. Missing inputs, read errors, and inapplicability are not
"no similar cases"; a selector's baseline fallback does not establish successful
memory inference. Resuming the same frozen plan after an interruption requires
no new consultation. Requery affected domains if inputs, metrics, mode, available
tools, or index change.

## Binding current encoding provenance

This section applies to all memory branches that require protein features.
Pocket memory uses `run_batch.py --task-type pocket_prediction` with the current
receptor only. Reuse compatible encodings of this
round's valid inputs. Otherwise, invoke the local SaProt encoder with only the
current-input manifest and this round's output directory. Docking `run_batch.py`
uses `--task-type docking`. General molecule-generation inputs use a one-row
manifest with `task_type=molecule_generation` and `--task-type all`, preserving
the actual task identity rather than relabeling another task's encoding. Set the
model, float32 inference, storage precision, windowing, and pLDDT policy to the
branch's frozen recipe. Valid generation-asset outputs from
`run_target_assets.py` may also be reused.

Raw `metadata.json` generally lacks the two top-level binding hashes required by
the query interface. After successful encoding, retain `prepared.json`,
`complete.json`, `metadata.json`, and `embeddings.npz`. Use the actual encoding
case directory to create a separate binding file, replacing the example paths
with this round's actual inputs:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" bind-query-metadata --encoding-dir "${SBDD_DATA_ROOT}/current_request/inputs/embeddings/docking/current_request/target1" --receptor-file "${SBDD_DATA_ROOT}/current_request/inputs/receptor.pdb" --output "${SBDD_DATA_ROOT}/current_request/inputs/query_metadata.json"
```

This command checks the original file SHA-256 values in the completion record,
recomputes preparation and encoding fingerprints in the actual encoder format,
and checks the current receptor against the protein SHA-256 in the original
preparation recipe before emitting `embedding_sha256` / `source_input_sha256`
and verification provenance. Generation-asset format additionally checks pocket
preparation and its completed artifacts. The output must be a new file. Do not
overwrite original encoding files or manually insert a "current receptor hash"
as a substitute for encoding provenance. Missing original binding evidence must
produce an explicit failure; re-encode this round's valid inputs or record a
justified fallback.

Subsequent selectors use this separate file for `--query-metadata`, while
`--query-npz` still points to the original NPZ. Binding establishes auditable
input provenance only: it does not modify task, encoder, or model parameters,
run scientific tools, or automatically establish compatibility with a memory
branch's encoding recipe. Branch queries still check recipes and numerical data.

## A consultation record

Keep `memory_consultation.md` and branch JSON outputs in the same request's
`inputs_dir`. Later rounds reference upstream records and document their own
tradeoffs without overwriting earlier records. Include at least:

- This round's objective; each domain's `selected / fallback / disabled /
  not_applicable / pending_inputs / error` status and specific reason. Use
  `no_matches` only when an actual query yields no candidates.
- Policy/index/feature versions, paths, and hashes; provenance of query ligand,
  receptor, encoding, and candidate sources, plus binding-check results; known
  case/group identity, exclusion rules, actually available tools, and query
  command.
- Per-metric recommendations, neighbors/weights/scores or their JSON locations;
  the actual pocket policy name, selected anchor, switch guard, family support, fusion
  members, and independent box-sizing rule. Uncalibrated scores are not success
  probabilities.
- Final engine/generation sources/pocket and config fields; adoption, override,
  or fallback reasons; the total multi-source candidate budget; unresolved
  limitations and tasks not covered by current evidence.

Recording these facts requires no additional database, persistent service,
extra approval, or new pipeline step. New results belong in this project's
reports and [artifact lineage](artifact-reuse.md), not automatically in frozen
supervised memory. Predicted scores are not new labels. In ordinary research mode,
follow the [research-memory append and recall protocol](research-memory.md) after
evaluation and before the next plan. The separate store preserves repeated runs
and exposes measured evidence to the planner, while frozen selector posteriors
remain tied to their declared snapshot. Benchmark evaluation and memory-off
ablations must neither read nor write this research store.

## Data and ablation boundaries

Docking labels measure quality for individual objectives; generation labels
identify the historical winning tool for each metric. Pocket memory uses raw
native-top-1 DCA to derive strict `<4 Å` labels. RMSD/PoseBusters, generation properties, and pocket DCA
cannot substitute for one another. Evaluation within memory must exclude both
the query case and its entire protein group. Fit scales and priors only on the
remaining training subset. Do not claim leave-group-out evaluation when the group
is unknown.

`context.docking_memory_mode=off` **disables only docking memory**, preserving the
original four-arm ablation meaning. State separately in the consultation record
whether generation memory is enabled, and declare the separate
`steps.pocket_analysis.params.pocket_memory_mode`. If an experiment declares
"all memory disabled," set pocket memory to `off`, avoid generation/docking
historical queries, and explicitly choose baselines.
Before reading history, fix isolated sessions, directories, and inputs for the
ablation. `off` must not reconstruct docking history from another arm's records.
Old `recall` is only for explicit legacy evidence inspection or reader checks,
not an additional hidden input to modern branches.

Do not use final/validation labels, Oracle results, research tuning outputs, or
the query's reference pose/box for online selection. Mixed-source CSVs required
for generation-index construction must first be filtered by frozen memory IDs
and split, emitting only memory records. Formula and artifact checks do not by themselves
establish scientific benchmark performance.

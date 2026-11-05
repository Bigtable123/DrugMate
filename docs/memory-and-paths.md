# Memory resources and deployment

DrugMATE has three separate memory stores: conversation state under
`DRUGMATE_HOME`, frozen scientific selector indexes, and evaluated research
observations under `SBDD_DATA_ROOT`. Starting the CLI and using conversation
memory do not require scientific models or benchmark downloads.

## Included resources

The [`memory/` inventory](../memory/README.md) lists the supplied indexes,
versions, case counts, sizes and file hashes:

| Resource | Included path | Purpose |
| --- | --- | --- |
| Docking index | `memory/docking_memory/ligand_protein_v2.json` | Ligand-conditioned docking recommendations |
| Generation index | `memory/generation_memory/ligand_protein_v1.json` | Five independent generation objectives |
| Query encoders | `memory/saprot_embedding/` | Prepare and encode current protein inputs |

Docking and generation selectors default to their distributed indexes and read
those files without opening a large source bundle. Use `--selector-index` to
select another index explicitly. `SBDD_MEMORY_ROOT` overrides the root used for
default index paths; when unset, it resolves to the distributed `memory/`
directory. An explicit `--selector-index` takes precedence.

The supplied indexes are frozen scientific snapshots. Their `memory_version`,
features, labels and source hashes identify that snapshot; they are not rewritten
for each deployment. Source paths are provenance references relative to
`SBDD_MEMORY_BUNDLES_ROOT`, not instructions to open files during an indexed
query. Keep rebuilt indexes and all new query inputs, encodings, bindings and
decisions outside the checkout.

The pocket selector needs a separate `pocket_protein_v1.json`, built from the
docking/pocket download. Neither distributed index contains its required DCA
labels. The commands below write this index under `SBDD_DATA_ROOT`.

## Configure paths

The default output root is `$HOME/drugmate-runs`. External sources and models
use `$HOME/drugmate-tools`, and benchmark resources use
`$HOME/drugmate-benchmarks`. Set deployment variables when your installations
live elsewhere; a configured path does not download files or start remote work.
Run the CLI on the machine that has the required scientific resources.

From the checkout root, with the scientific Conda environment active:

```bash
export BIO_SKILLS_ROOT="$PWD/skills/bio_research"
export SBDD_MEMORY_ROOT="$PWD/memory"
export SBDD_DATA_ROOT="$HOME/drugmate-runs"
export SBDD_SOURCE_ROOT="$HOME/drugmate-tools"
export SBDD_EXTERNAL_ROOT="$HOME/drugmate-tools"
export SBDD_BENCHMARK_ROOT="$HOME/drugmate-benchmarks"
export SBDD_ENV_DRUGMATE_PREFIX="$CONDA_PREFIX"
export SBDD_PYTHON="$CONDA_PREFIX/bin/python"

conda activate drugmate
python drugmate --cli
```

`SBDD_PYTHON` selects the scientific orchestration interpreter; individual tools
still need their registered environments, source trees and models. Do not set
both `SBDD_ENV_<ID>_PREFIX` and `SBDD_ENV_<ID>_NAME` for the same environment.
See the [deployment reference](../skills/bio_research/pipeline-entry/references/deployment.md)
for per-tool overrides and the [local deployment file](../skills/bio_research/pipeline-entry/references/deployment.md#local-deployment-file)
for persistent settings. Keep model credentials in local DrugMATE settings.

Shell examples expand environment variables. JSON configuration files must
contain actual resolved paths: JSON does not expand `$HOME` or `$SBDD_DATA_ROOT`.
Keep `SBDD_DATA_ROOT` outside the checkout and outside the downloaded bundles;
do not use the repository's `memory/` directory as an experiment output root.

## Download and placement

Use the Google Drive [download links in the main README](../README.md#research-memory).
The two archives are needed for index rebuilding and their original memory
readers; modern docking/generation queries can use the included indexes alone.
Downloads do not include installed scientific tools or encoder model weights.

From the checkout root, adjust `SBDD_DOWNLOADS_DIR` to your download directory:

```bash
export SBDD_DOWNLOADS_DIR="$HOME/Downloads"
export SBDD_MEMORY_BUNDLES_ROOT="$PWD/memory/bundles"
mkdir -p "$SBDD_MEMORY_BUNDLES_ROOT"
unzip "$SBDD_DOWNLOADS_DIR/drugmate-generation-memory.zip" -d "$SBDD_MEMORY_BUNDLES_ROOT"
unzip "$SBDD_DOWNLOADS_DIR/drugmate-docking-memory.zip" -d "$SBDD_MEMORY_BUNDLES_ROOT"
```

Each bundle directory must sit directly under `SBDD_MEMORY_BUNDLES_ROOT`:

```text
memory/bundles/
  drugmate-generation-memory/
    memory/
      bundle_manifest.json
      memory/agent_memory.jsonl
      memory/memory_index.sqlite
      inputs/molecule_generation/
      read_memory.py
    evaluation_v2/target_best_tools.csv
    audit/similarity_groups.csv
    distribution/
  drugmate-docking-memory/
    memory/agent_memory.jsonl
    memory/memory_index.sqlite
    inputs/docking/
    inputs/pocket_prediction/
    code/bundle_reader.py
    distribution/
```

The tree shows key entrypoints, not every required file. Preserve the complete
extracted directories. Generation has two internal `memory/` levels in its
record path; docking has one. Do not add an extra wrapper directory.

An external disk is also suitable: set `SBDD_MEMORY_BUNDLES_ROOT` to the absolute
parent of the two extracted directories before starting the CLI. This variable
takes precedence for locating these bundles; when unset, it defaults to
`bundles/` under the resolved memory root. It is separate from the output root
and the index directory. Save it in local deployment settings if needed.
`memory/bundles/` is excluded from Git, source distributions and wheels.

A supplied ZIP checksum verifies only that exact archive and is not a runtime
input. Original manifests, input bytes and content hashes inside each bundle
must remain unchanged. Do not rewrite paths inside frozen records, filter
`agent_memory.jsonl`, trim its fields or remove referenced evidence. SQLite
indexes can depend on JSONL byte offsets. Evidence assets, log excerpts and
runtime supplements are part of the scientific inputs, not disposable caches.

## Package boundaries and known missing assets

The generation download contains its memory reader, case inputs, retrieval data
and supporting evaluation/split tables. Its `distribution/` inventories describe
the supplied files and omissions. It lacks 49 memory and 12 held-out `cavity.as`
files. The current generation index builder and memory reader do not open them;
cavity-dependent preparation and full release validation still require those
files. Availability for indexed selection does not establish complete experiment
reproducibility.

The docking download contains memory cases for docking and pocket prediction,
plus their reader, inputs and referenced evidence. It does not include held-out
answer sets. Some evidence is explicitly external or unavailable; the reader
reports that status rather than fetching missing sources. Retain the release's
read-only access boundary.

Legacy `recall`, evidence drill-down and planning cards require the complete
compatible docking reader and its data. A `planning/` sidecar or a copied
`bundle_reader.py` alone is not sufficient. The generation bundle's
`memory/read_memory.py` is a different reader. Consult each bundle's
`distribution/` records before attempting operations beyond indexed selection.

## Query encodings and evaluation isolation

Current protein queries need a compatible local SaProt/Foldseek environment and
model assets, or compatible precomputed encodings. Follow the
[encoder guide](../memory/saprot_embedding/README.md). Model/tokenizer identity,
dtype, storage precision, preprocessing, windowing, pooling and confidence-mask
policy must match the index contract; a 1280-dimensional vector alone does not
prove compatibility.

Keep the encoder's original `prepared.json`, `complete.json`, `metadata.json`
and `embeddings.npz`, plus any declared preparation inputs. Use
`benchmark_memory.py bind-query-metadata` to verify and bind them to the current
receptor. It does not compute embeddings. Never edit recorded hashes or encoder
identities to make an old encoding appear compatible with changed inputs.

During held-out evaluation, expose only the permitted read-only memory and
current query. Evaluation answers and administrative tables must not enter
online retrieval, priors or threshold selection. The generation builder uses
its supporting CSVs only for allowlisted memory cases and `split=memory` rows.
Keep these offline builder inputs separate from agent-visible query inputs.
Exclude known query cases and protein groups as required by the task protocol.

## Rebuilding indexes

Set the deployment variables above and use a compatible scientific environment.
These commands rebuild derived indexes without running docking, generation or
pocket detectors. Outputs must remain outside the frozen bundles and checkout:

```bash
"$SBDD_PYTHON" "$BIO_SKILLS_ROOT/pipeline-entry/scripts/benchmark_memory.py" \
  build-docking-index \
  --benchmark-bundle "$SBDD_MEMORY_BUNDLES_ROOT/drugmate-docking-memory" \
  --output "$SBDD_DATA_ROOT/docking_memory/ligand_protein_v2.json"

"$SBDD_PYTHON" "$BIO_SKILLS_ROOT/pipeline-entry/scripts/benchmark_memory.py" \
  build-generation-index \
  --benchmark-bundle "$SBDD_MEMORY_BUNDLES_ROOT/drugmate-generation-memory/memory" \
  --property-winners "$SBDD_MEMORY_BUNDLES_ROOT/drugmate-generation-memory/evaluation_v2/target_best_tools.csv" \
  --groups "$SBDD_MEMORY_BUNDLES_ROOT/drugmate-generation-memory/audit/similarity_groups.csv" \
  --output "$SBDD_DATA_ROOT/generation_memory/ligand_protein_v1.json"

"$SBDD_PYTHON" "$BIO_SKILLS_ROOT/pipeline-entry/scripts/benchmark_memory.py" \
  build-pocket-index \
  --benchmark-bundle "$SBDD_MEMORY_BUNDLES_ROOT/drugmate-docking-memory" \
  --output "$SBDD_DATA_ROOT/docking_memory/pocket_protein_v1.json"
```

Pass a rebuilt index through `--selector-index`; rebuilding it does not replace
the distributed default. Preserve its generated provenance and version alongside
the supplied snapshot. Do not remove historical cases to shrink an index:
changing the memory population changes the selector's results.

## What each mode requires

| Mode | Required resources and behavior |
| --- | --- |
| Conversation memory | Included runtime; private notes and session SQLite state under `DRUGMATE_HOME` |
| Docking `ligand_protein` | Docking index, one query ligand, compatible current protein encoding and bound provenance |
| Docking `ligand_only` | Docking index and one query ligand; no current protein encoding required |
| Docking `legacy` / explicit `recall` | Complete compatible reader release and permitted query inputs; the legacy selector supports redocking only |
| Docking `off` | No docking history or query-feature reads; selected task's GNINA/redocking or rDock/pose-generation baseline must be available |
| Generation memory | Generation index, one conditioning reference ligand, compatible current protein encoding and bound provenance |
| Pocket memory `auto` / `required` | Current receptor, separate pocket index and compatible bound protein encoding; plan detectors before execution and supply their current predictions for fusion |
| Pocket memory `off` | Current receptor and four-detector evidence; fixed DeepPocket#1 rank fusion without historical inputs |

Pocket `auto` can report an unavailable-memory fallback; `required` rejects that
fallback. Corrupt evidence and invalid detector outputs remain errors. Check
actual tool inputs, environments and resources before adopting any recommendation.
The selectors provide planning advice; the agent explicitly configures and
validates the complete tool matrix. Follow the
[docking](../skills/bio_research/pipeline-entry/references/docking-memory.md),
[generation](../skills/bio_research/pipeline-entry/references/generation-memory.md)
and [pocket](../skills/bio_research/pipeline-entry/references/pocket-memory.md)
protocols for the full contracts.

`context.docking_memory_mode=off` disables only docking memory. An all-memory-off
run must also disable pocket memory, skip generation/history queries and skip
the research store. Record missing inputs, failed retrieval and baseline
fallbacks explicitly; do not report them as successful retrieval with no hits.
Follow the existing Gate A/Gate B requirements before changing a study arm or
proceeding with a fallback.

Ordinary memory-enabled research appends verified evaluated runs to the separate
store under `SBDD_DATA_ROOT/research_memory` and recalls them before later
planning. Frozen evaluations and memory-off runs neither read nor write this
store. It preserves independent observations without overwriting frozen labels
or assigning measured winners to untested tools. See the
[research-memory protocol](../skills/bio_research/pipeline-entry/references/research-memory.md).

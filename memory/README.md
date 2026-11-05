# Scientific memory assets

This directory contains the supplied frozen selector indexes and query-encoding
entrypoints. Ordinary DrugMATE conversation memory uses `DRUGMATE_HOME` and does
not depend on these files.

```text
memory/
  docking_memory/ligand_protein_v2.json
  generation_memory/ligand_protein_v1.json
  saprot_embedding/
    run_batch.py
    run_target_assets.py
    README.md
  bundles/                         # Separate downloads; ignored by Git
    drugmate-generation-memory/
    drugmate-docking-memory/
```

## Included indexes

| Index | Schema | Memory cases | Bytes |
| --- | --- | ---: | ---: |
| `docking_memory/ligand_protein_v2.json` | `sbdd_docking_selector_index_v2` | 487 | 14,210,940 |
| `generation_memory/ligand_protein_v1.json` | `generation_selector/1` | 141 | 4,076,404 |

SHA-256 of the supplied indexes:

```text
docking_memory/ligand_protein_v2.json
ac6fbc0c642d596d8ccc10d8126333075fba96f297e3683a8e0aa0393c6392d5

generation_memory/ligand_protein_v1.json
6c80168f54bd019cbc8e58ef797174853c9d47d17b7cc46084e8f0e5aecccc57
```

The docking v2 index includes Top-1 RMSD, PoseBusters, joint-success and
RMSD objectives. The generation index uses its own `v1` format. Both indexes
declare the `memory` split and retain release identities, source hashes,
scientific measurements and labels. Source locations are relative to the
configured bundle parent.

Neither supplied index is the pocket selector's DCA index. Build
`pocket_protein_v1.json` separately from the docking/pocket bundle into an output
directory outside this checkout, then pass its path explicitly; see the
[pocket memory protocol](../skills/bio_research/pipeline-entry/references/pocket-memory.md#prepare-the-pocket-index-and-current-query).
The builder has been checked against the real downloaded bundle: 407 pocket
cases and 267 protein groups. Detector routing, fusion and implementation choices
are documented in the
[pocket memory protocol](../skills/bio_research/pipeline-entry/references/pocket-memory.md).

The selectors do not append new experiments to these frozen files. Use the
[research-memory protocol](../skills/bio_research/pipeline-entry/references/research-memory.md)
to retain evaluated runs in a separate store under `SBDD_DATA_ROOT` and retrieve
their evidence in later planning. New prepared-input identity creates a profile;
another run for the same identity appends evidence. Identity uses input hashes,
task, feature and encoder versions, not floating point vector equality.
A 1280-dimensional vector here represents the protein; ligand features are a
fingerprint and descriptors. An untested tool's quality remains unknown, so a
single new result cannot establish the generation selector's comparative winners.

## Use the local copies

Selectors use the distributed `memory/` directory by default. Set
`SBDD_MEMORY_ROOT` through the
[deployment configuration](../skills/bio_research/pipeline-entry/references/deployment.md#local-deployment-file)
to use another installation. An explicit `--selector-index` overrides that root;
for example:

| Command | Index argument |
| --- | --- |
| `select-docking` | `--selector-index "$SBDD_MEMORY_ROOT/docking_memory/ligand_protein_v2.json"` |
| `select-generation` | `--selector-index "$SBDD_MEMORY_ROOT/generation_memory/ligand_protein_v1.json"` |

Follow the complete [docking](../skills/bio_research/pipeline-entry/references/docking-memory.md)
or [generation](../skills/bio_research/pipeline-entry/references/generation-memory.md)
query protocol for ligand, receptor, encoding, metadata and authorization
requirements. The index loaders read their JSON directly without reopening
the source bundles.

Keep every new query input, embedding, metadata binding, selection result and
rebuilt index under `SBDD_DATA_ROOT`, outside this checkout. Do not point
`SBDD_DATA_ROOT` at this directory. Use explicit index arguments instead of
relying on an unreviewed index location.

## What remains external

The [encoder entrypoints](saprot_embedding/README.md) are included, but the
compatible scientific environments, model weights and preparation tools must
be provisioned separately. Current queries require a compatible protein
encoding with its original provenance, or compatible precomputed query assets;
the historical indexes do not supply new query encodings.

Large benchmark bundles are distributed separately through the Google Drive
[download links in the main README](../README.md#research-memory). Extract
`drugmate-generation-memory.zip` and `drugmate-docking-memory.zip` into
this directory's `bundles/` subdirectory. The resulting bundle directories must
be `bundles/drugmate-generation-memory/` and `bundles/drugmate-docking-memory/`;
do not add another `memory/` layer. From the checkout root, set
`SBDD_MEMORY_BUNDLES_ROOT="$PWD/memory/bundles"` before launching the CLI.
`SBDD_MEMORY_BUNDLES_ROOT` defaults to the `bundles/` directory under the
resolved `SBDD_MEMORY_ROOT`; set it explicitly for external storage. See the
[placement guide](../docs/memory-and-paths.md#download-and-placement) for full
commands and external-disk placement. Downloaded bundles are excluded from
Git and Python packages; they do not replace the tracked indexes above.
The supplied generation package contains its own `memory/read_memory.py`,
SQLite index, retrieval arrays and case inputs. The supplied docking package
contains the original `code/bundle_reader.py`, its supporting modules and data
for legacy recall, evidence drill-down and planning-card construction. It retains
894 historical contexts: 487 docking and 407 pocket-prediction cases. Its
evidence assets, log excerpts, databases and runtime supplement remain part of
that frozen release. A generated `planning/` sidecar is not included.
The generation reader does not replace the docking reader. The two included
indexes support modern selectors without either large source bundle.

If ZIP checksum files are supplied, they support optional transfer verification
and are not runtime inputs; see the
[placement guide](../docs/memory-and-paths.md#download-and-placement). Checksums
for the earlier `.tar.gz` archives do not apply to these ZIP files.
In the prepared generation bundle, 49 memory and 12 validation `cavity.as`
files were already absent. The current generation selector and bundled reader
do not open them, but cavity-dependent preprocessing and full original-release
validation still require those files. The package records their paths in
`distribution/missing-before-cleanup.tsv`.

See [memory deployment details](../docs/memory-and-paths.md) for encoder identity
rules, external inputs and index rebuilding.

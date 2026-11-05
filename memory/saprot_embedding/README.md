# SaProt query encoders

Two query-encoding entrypoints are available:

- `run_batch.py`: protein inputs described by a CSV; docking and pocket tasks
  have explicit filters, while a generation query uses its task-labelled row
  with `--task-type all` as specified by the scientific protocol.
- `run_target_assets.py`: generation target-asset directories; reuses
  `run_batch.py` for protein preparation and encoding.

Run them with a compatible scientific Python environment.
They require NumPy, Biopython, PyTorch and Transformers, a compatible local
SaProt model directory and the Foldseek executable. Target-assets encoding also
uses SciPy and, when cavity files are present, rDock's `rbcavity` and data files.
The scripts use POSIX file locking. These dependencies and model weights are not
included or installed by DrugMATE's application environment. Model loading uses
`local_files_only=True` and does not download weights.

## Deployment settings

Supply paths explicitly through CLI arguments or environment variables. CLI
arguments take precedence.

| CLI argument | Environment variable | Meaning |
| --- | --- | --- |
| `--manifest` | `SAPROT_MANIFEST` | Input CSV for `run_batch.py` |
| `--assets` | `SAPROT_TARGET_ASSETS` | Target-asset root for `run_target_assets.py` |
| `--output` | `SAPROT_OUTPUT` | Encoding output root under the allocated research output directory |
| `--model` | `SAPROT_MODEL` | Compatible local SaProt model, tokenizer and weight files |
| `--foldseek` | `SAPROT_FOLDSEEK` | Existing Foldseek executable |
| `--rbcavity` | `SAPROT_RBCAVITY` | Native rDock cavity reader; required for cavity-bearing inputs |
| `--rbt-root` | `SAPROT_RBT_ROOT` | Existing rDock data directory; required with the native reader |

The local model is not needed for `--prepare-only`; Foldseek still runs during
preparation, and cavity inputs still require the native rDock check. Store
outputs outside this source directory and outside frozen memory releases.

From the checkout root, with these environment variables set and the scientific
environment activated, the retained entrypoints are:

```bash
# Configure CPU library thread limits for the selected environment.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"

python -u memory/saprot_embedding/run_batch.py --task-type docking
python -u memory/saprot_embedding/run_target_assets.py
```

These are separate operations: invoke the one required by the reviewed task.
Follow the scientific skill's input review, output placement and authorization
contracts before running them.

## Inputs and selection

The batch CSV must contain `task_type`, `dataset`, `target` and `protein_path`.
Docking rows additionally require `known_pocket_box_path`. Preserve supplied
protein/receptor hashes, validation flags and explicit pLDDT policy fields;
the script checks those bindings during preparation. Path values inside the CSV
are literal filesystem paths, not shell-variable templates.

`--task-type` selects `docking`, `pocket_prediction` or `all`; `--dataset` and
`--target` can be repeated. `--limit` bounds the selected CSV rows. The
target-assets encoder supports `--target` and expects the established
`ADB_*`, `PB_*` or `PL_*` directory names, with `metadata.json` and `protein.pdb`.
If `cavity.as` exists, `protein.mol2` is also required. A missing cavity is
recorded explicitly; it is not synthesized from a reference ligand.

## Encoding contract and outputs

Both entrypoints preserve the supplied encoding implementation: separate
TER/chain fragments, residue mappings, confidence policy, Foldseek conversion,
overlapping windows, float32 pooling and output bindings. Current defaults are
`--window-size 1022`, `--window-overlap 128`, `--storage-dtype float16` and
`--dtype auto` (`float16` on CUDA, `float32` on CPU). Generation cavity selection
uses `--cavity-distance 4.0` and checks the parsed grid against native rDock.
Match the frozen index's model hashes, dtype, window, pooling and preparation
contract when producing a query; a matching vector dimension is insufficient.

Each selected case writes `prepared.json`, `embeddings.npz`, `metadata.json`
and `complete.json` beneath `OUTPUT/task_type/dataset/target/`. Generation
outputs also include `pocket_prepared.json`, plus cavity geometry and pocket
arrays when available. Resume skips a case only when its configuration
fingerprint, completion marker and output hashes match; `--overwrite` requests
recomputation. Newly generated source hashes reflect this release's scripts;
keep frozen index contents and original provenance unchanged.

Bind compatible query assets using the scientific pipeline's
[`benchmark_memory.py`](../../skills/bio_research/pipeline-entry/scripts/benchmark_memory.py)
`bind-query-metadata` operation before selection. The scripts themselves do not
rebuild the frozen docking or generation memory indexes.


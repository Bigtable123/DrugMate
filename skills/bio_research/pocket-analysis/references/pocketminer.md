# PocketMiner: isolated installation and smoke test

PocketMiner is an experimental cryptic-pocket predictor from the official [`Mickdub/gvp`](https://github.com/Mickdub/gvp) repository's `pocket_pred` branch. Keep it in a **separate Conda environment**: its TensorFlow-era dependency stack is incompatible with the broader `drugmate` environment, and the repository's official YAML contains a machine-specific `prefix` that must not be copied into Drugmate.

This guide pins commit `187062df3c94127e991669768009141a08fd5d8b`. Do not run from a floating branch in a reproducible experiment.

## 1. Clone and pin the source

```bash
: "${SBDD_POCKETMINER_SOURCE_DIR:?Set the PocketMiner source directory}"
export PM_ROOT="${SBDD_POCKETMINER_SOURCE_DIR}"
mkdir -p "$(dirname "$PM_ROOT")"

git clone --branch pocket_pred --single-branch \
  https://github.com/Mickdub/gvp.git "$PM_ROOT"
git -C "$PM_ROOT" checkout --detach \
  187062df3c94127e991669768009141a08fd5d8b
git -C "$PM_ROOT" rev-parse HEAD
```

The last command must print exactly:

```text
187062df3c94127e991669768009141a08fd5d8b
```

## 2. Create a dedicated environment

The upstream README reports TensorFlow `2.9.1` as tested. Python 3.9 is a practical compatible baseline. Create the environment explicitly instead of applying the upstream YAML's local `prefix`:

```bash
conda create -n pocketminer --override-channels -c conda-forge -y \
  python=3.9 pip numpy=1.23.5 scipy=1.9.3 pandas=1.5.3 \
  tqdm=4.64.1 mdtraj=1.9.7

conda run -n pocketminer python -m pip install \
  "tensorflow==2.9.1"
conda run -n pocketminer python -m pip check
```

This is a modern compatibility environment, not a byte-for-byte recreation of
the original Python 3.7.6 / TensorFlow 2.1.0 stack.  Upstream also reports
TensorFlow 2.6.2 and 2.9.1 as tested.  Let pip resolve TensorFlow's compatible
`protobuf` dependency, then use `pip check` above to detect an inconsistent
environment rather than upgrading `protobuf` independently.

Do **not** install these packages into `drugmate`. It is safe for the current shell
to have Drugmate activated as long as every PocketMiner command still uses
`conda run -n pocketminer`. If you choose to derive an environment from the
official `pocketminer.yml`, first make a copy outside the repository and remove
its `prefix:` line; still give the environment a distinct name.

## 3. Import smoke test

```bash
conda run -n pocketminer python -c \
  "import tensorflow as tf, mdtraj, numpy, scipy, pandas; print('tensorflow', tf.__version__); print('PocketMiner imports OK')"
```

Expected: TensorFlow reports `2.9.1`, all imports complete, and the command exits with status 0. Prediction does not require a GPU; TensorFlow CUDA warnings on a CPU-only host are not by themselves a failed smoke test.

Import success alone is insufficient.  Incompatible TensorFlow/Keras releases can
silently leave model variables at random initialization because upstream
`util.load_checkpoint()` does not inspect the restore status.

## 4. Run a deterministic checkpoint smoke test

The README mentions an ACE2 demo, but that file is absent from this pinned
checkout.  Use a structure that is actually bundled with the PocketMiner dataset
and the parameterized `xtal_predict_UI.py` entry point.  Run it in two independent
processes: identical input must produce identical text output.

```bash
export PM_TEST_PDB="$PM_ROOT/data/pm-dataset/apo-structures/1RRGA_clean_h.pdb"
PM_SMOKE_DIR="$(mktemp -d "${SBDD_DATA_ROOT:?Set the experiment output root}/pocketminer-smoke.XXXXXX")"
mkdir "$PM_SMOKE_DIR/run1" "$PM_SMOKE_DIR/run2"

conda run -n pocketminer python "$PM_ROOT/src/xtal_predict_UI.py" \
  "$PM_ROOT/models/pocketminer" "$PM_TEST_PDB" "$PM_SMOKE_DIR/run1"
conda run -n pocketminer python "$PM_ROOT/src/xtal_predict_UI.py" \
  "$PM_ROOT/models/pocketminer" "$PM_TEST_PDB" "$PM_SMOKE_DIR/run2"

test -s "$PM_SMOKE_DIR/run1/predictions.txt"
test -s "$PM_SMOKE_DIR/run2/predictions.txt"
cmp "$PM_SMOKE_DIR/run1/predictions.txt" "$PM_SMOKE_DIR/run2/predictions.txt"
```

Any difference is a hard failure, even when both files contain finite values in
the interval 0–1.  The pipeline adapter must additionally build the model before
restore and call TensorFlow's `assert_existing_objects_matched()`; the two-run
test is a smoke check, not a replacement for strict checkpoint validation.

The checkpoint argument is the prefix `$PM_ROOT/models/pocketminer`, without
`.index` and without treating it as a directory.  The UI script supports one
PDB per invocation and writes a fixed `predictions.txt`, so every run needs its
own output directory.

Keep the pinned source tree and record its commit, environment export, input structure, checkpoint, and output paths with every experiment. Run the script on a small known structure before scheduling a batch.

## Integration status

PocketMiner is connected to Step 1 as an isolated auxiliary tool. The Agent
decides whether to enable it through
`../../pipeline-entry/TOOL_SELECTION_POLICY.md`, and config assembly records
that decision as `steps.pocket_analysis.tools.pocketminer.enabled=true|false`.
The runner does not enable it from structure facts. When selected, the resolved
runtime invokes `pocketminer_bridge.py` in the dedicated `pocketminer` Conda
environment and performs all of the following before scores can influence
candidate ordering:

- verifies the pinned repository commit and checkpoint files;
- rejects tracked source changes and verifies the pinned checkpoint index/data SHA-256 digests;
- builds the complete model and requires
  `assert_existing_objects_matched()` to pass;
- maps one requested PDB chain by `(chain, resSeq, insertion code)` and requires
  N/CA/C/O for every encoded residue;
- checks prediction count, finite values, the sigmoid range, and repeated
  `train=False` inference determinism.

Run the exact pipeline preflight directly with:

```bash
conda run -n pocketminer python \
  "${BIO_SKILLS_ROOT}/pocket-analysis/scripts/pocketminer_bridge.py" \
  probe --source-dir "${SBDD_POCKETMINER_SOURCE_DIR}" \
  --checkpoint "${SBDD_POCKETMINER_SOURCE_DIR}/models/pocketminer"
```

A successful process exit alone does not establish checkpoint compatibility.
The strict probe must match every model object and produce repeatable scores;
otherwise the pipeline reports the tool unavailable.

Even after the probe passes, PocketMiner remains residue-level evidence. It does
not join the cavity-detector `one_of`, increment detector support, report volume,
or independently publish a canonical pocket/docking box. Its provisional
high-score residue clusters are hypotheses for visualization or conformational
sampling only.

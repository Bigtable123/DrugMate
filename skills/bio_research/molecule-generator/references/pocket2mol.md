# Pocket2Mol Reference

## Overview

Pocket2Mol is an equivariant graph neural network model for 3D structure-based molecule generation.
It generates drug-like molecules that fit into protein binding pockets.

- GitHub: https://github.com/pengxingang/pocket2mol
- Paper: Peng et al., ICML 2022

## Installation

```bash
# Create dedicated conda environment
conda env create -f env_cuda113.yml
# Or manually:
conda create -n Pocket2Mol python=3.8.12
conda activate Pocket2Mol
conda install pytorch==1.10.1 cudatoolkit=11.3 -c pytorch -c conda-forge
conda install pyg -c pyg
conda install -c conda-forge rdkit biopython pyyaml easydict python-lmdb
```

Required model checkpoint for this workspace:
`${SBDD_POCKET2MOL_SOURCE_DIR}/pretrained_Pocket2Mol.pt`.
Its registry-pinned SHA256 is
`e9d8152dfdf034dabf48e566904e9d7c6cc1b96843bdc4d0026d891678913911`;
readiness and execution both reject a mismatched asset.

## Which Config to Use

Use `configs/sample_for_pdb.yml` for this SBDD pipeline. It is paired with
`sample_for_pdb.py`, which accepts a prepared protein PDB plus an explicit
pocket center.

Do not use `configs/sample.yml` for this pipeline. It is paired with `sample.py`
and the bundled CrossDocked test-set workflow (`--data_id`), where pocket data
comes from `data.crossdocked_pocket10` instead of a user-provided PDB.

## CLI Usage

```bash
# Generate molecules for a custom PDB pocket
python sample_for_pdb.py \
    --pdb_path ./pocket.pdb \
    --center "32.0,28.0,36.0" \
    --config ./configs/sample_for_pdb.yml \
    --outdir ./output/

# IMPORTANT: if the first coordinate is negative, add a leading space:
python sample_for_pdb.py \
    --pdb_path ./pocket.pdb \
    --center " -1.5,28.0,36.0"
```

## Input Format

| Parameter | Type | Description |
|-----------|------|-------------|
| `--pdb_path` | string | Path to protein pocket PDB file |
| `--center` | string | Bounding box center as "x,y,z" (space-prefix if negative first value) |
| `--bbox_size` | float | Bounding box side length (default: 23 Å) |
| `--config` | string | YAML config containing checkpoint path, `sample.num_samples`, and `sample.seed` |
| `--outdir` | string | Output directory for generated SDF files |

## Output Format

Generated molecules are saved as SDF files in the output directory.
Each SDF entry contains 3D coordinates of atoms placed within the pocket.

## Performance Tips

- **Use `taskset -c 0`** to restrict to one CPU core during sampling — significantly faster.
- This skill enables that by default through `pocket2mol_taskset_cpu: 0`; set it to `-1` to disable.
- GPU acceleration: set `CUDA_VISIBLE_DEVICES` appropriately.
- Typical speed: ~100 molecules in 5–10 minutes on a single GPU.
- Memory: ~4 GB GPU memory for standard pockets.

## Integration Notes (Local Execution)

This skill runs Pocket2Mol locally via `conda run` because:

1. Requires PyTorch + PyG + CUDA (heavy dependencies)
2. Requires GPU for practical speed
3. Dependencies conflict with the main SBDD skill RDKit-only environment

Execution flow:

1. Allocate a fresh `pocket2mol_work/run-*` invocation directory
2. Copy the pocket PDB into that invocation and write a config with the registered checkpoint, requested sample count, and pipeline seed
3. Run `taskset -c 0 python sample_for_pdb.py` via `conda run -n Pocket2Mol`
4. Save `command.json`, `stdout.log`, `stderr.log`, and `execution.json`; the
   execution record distinguishes command success from parsed-output QC success
5. Parse SDF files only below that invocation's fresh `output/` directory
6. Retain only sanitized, single-component records with a conformer marked as 3D and finite coordinates; zero valid records make the invocation unsuccessful
7. Preserve compatible raw coordinates in `library.sdf`; generate ETKDG coordinates only when the standardized graph no longer matches or no raw pose exists

The default Pocket2Mol timeout is 24 hours (`86400` seconds). Each invocation is
isolated even when it reuses the same Step 2 output directory; a failed or empty
run cannot consume SDF files from an earlier invocation.

## Known Issues

| Issue | Workaround |
|-------|-----------|
| No macOS support | Run on Linux with CUDA |
| Model checkpoint missing or mismatched | Restore the registry-pinned checkpoint; readiness and execution both verify its SHA-256 |
| Slow with multiple CPUs | Use `taskset -c 0` for single-core |
| Large bbox generates too many atoms | Keep default 23 Å or reduce |
| Some generated SMILES invalid | Filter with RDKit post-generation |

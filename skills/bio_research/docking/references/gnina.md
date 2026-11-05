# GNINA Reference

## Overview

GNINA is a molecular docking program based on AutoDock Vina with integrated CNN scoring functions for pose prediction and binding affinity estimation.

- GitHub: https://github.com/gnina/gnina
- Paper: McNutt et al., J Cheminform 2021

## Installation

```bash
# Download pre-built binary (Linux, CUDA optional)
wget https://github.com/gnina/gnina/releases/download/v1.1/gnina
chmod +x gnina
install -m 755 gnina "${SBDD_ENV_DRUGMATE_PREFIX:?Set the scientific environment prefix}/bin/gnina"

# macOS: no official binary, use Docker
docker pull gnina/gnina:latest
```

## Key CLI Options

```bash
gnina \
  --receptor protein.pdb \       # Receptor file (PDB/PDBQT)
  --ligand library.sdf \         # Ligand(s) (SDF/MOL2/PDBQT)
  --autobox_ligand ref.sdf \     # Define box from reference ligand
  --autobox_add 4 \              # Padding around autobox (Angstroms)
  --exhaustiveness 16 \          # Search effort (default 8, recommend 16+)
  --num_modes 20 \                # Default poses per ligand; configurable
  --cnn_scoring rescore \        # CNN scoring mode: none|rescore|refinement|all
  --out docked.sdf \             # Output file
  --log gnina.log                # Log file with scores
```

## Scoring Modes

| Mode | Speed | Accuracy | Use Case |
|------|-------|----------|----------|
| `--cnn_scoring none` | Fast | Low | Vina-only, large libraries |
| `--cnn_scoring rescore` | Medium | High | Default: Vina search + CNN rescore |
| `--cnn_scoring refinement` | Slow | Highest | CNN-guided refinement after Vina |

## Output SDF Properties

Each pose in output SDF has these properties:
- `minimizedAffinity`: Vina energy after minimization (kcal/mol)
- `CNNscore`: Pose quality score (0-1)
- `CNNaffinity`: Predicted binding affinity on a pK scale (higher is better; not kcal/mol)

When supplied by the selected GNINA build/model, the adapter also preserves
finite `intramol` and `CNNvariance` values. Canonical score records include a
one-based `engine_pose_rank` derived from output order within each molecule
name. These three fields are audit context, not weighted evidence.

GNINA source sorts `CNNaffinity` in descending order and empirical Energy in ascending order. The GNINA 1.0 paper describes the affinity model output as pK.

## Executable screening policy

Use the canonical `minimizedAffinity` field for the Vina-like energy; do not
substitute an undocumented `affinity` alias. The entries below are level-1
per-metric marks, not molecule decisions: the shared scorer folds the resulting
`pass`/`flag`/`fail` marks by category before assigning `keep`/`modify`/`reject`.
Apply the fallback gates only to finite values from one target/pocket cohort run
with the same GNINA version, CNN model/configuration, ligand preparation, box,
and search settings.

| Metric | Pass mark | Flag mark | Fail mark | Suggested signed weight |
|---|---|---|---|---:|
| `minimizedAffinity` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -5 |
| `CNNscore` | `>= 0.50` | `[0.20, 0.50)` | `< 0.20` | +4 |
| `CNNaffinity` | `>= 6` pK | `[2, 6)` | `< 2` | +3 |

Normalize applicable non-zero-weight metrics with the same declared method for
the entire cohort, orient lower-is-better energy before weighting, compute
`sum(normalized_value * weight)`, and rank by descending total. Missing,
non-finite, or configuration-inapplicable values are NA and do not enter the
normalization or total. In particular, CNN fields are normally NA when
`--cnn_scoring none` is selected.

GNINA is Vina-derived. Its `minimizedAffinity` must not be added to Uni-Dock
`docking_score` or AutoDock Vina `vina_affinity` for the same molecule. When
cross-engine evidence is available, choose one predeclared Vina-family energy
slot. `CNNscore` and `CNNaffinity` remain distinct, namespaced GNINA evidence.

`intramol`, model-ensemble `CNNvariance`, pose rank, and any score-gap summary
have weight 0. Rank `<=20` is inside the declared export window; a larger rank
is rejected from the pose-export set, while the parent molecule is mapped to
modify so an alternate pose or rerun can rescue it. Preserve the
other fields for audit when a validated engine configuration
produces them. The adapter preserves finite `intramol` and `CNNvariance` values
when present and publishes `engine_pose_rank`; absent optional properties are
NA rather than zero or failure. Because current GNINA lineage groups rank by
molecule name, duplicate names belonging to different inputs make that audit
rank ambiguous. A small score gap may mark pose ordering as uncertain, but
cannot reject a molecule by itself.

## Native-redocking policy

The implemented redocking route is GNINA-specific and reports the lowest
symmetry-aware heavy-atom RMSD found across all generated reference-ligand
modes. The default soft/hard parameters are `2.0` and `2.5` Å respectively.
The validation command replays the main run's `num_modes`, `pose_sort_order`,
`cnn_scoring`, exhaustiveness, device, and seed so the RMSD gate evaluates the
same sampling/scoring protocol.
Apply it as a target×pocket×protocol gate:

| Best redocking RMSD | Protocol action |
|---|---|
| `<= 2.0 Å` | keep |
| `(2.0, 2.5] Å` | modify/recalibrate |
| `> 2.5 Å` | reject and roll back the protocol |

This is best-of-generated-modes RMSD, not the RMSD of the score-top-1 pose and
not a top-N active/decoy success rate. It never rejects an individual candidate
molecule. Disabled or reference-free redocking is NA/not-applicable, not a
pass and not a molecular score.

## Performance Tips

- Use `--exhaustiveness 16-32` for production runs
- `num_modes` requests up to `20` per molecule, capped to 20/5/1 by cohort size; users can request fewer in the
  experiment config. Report actual output counts against the requested budget.
- For large libraries (>1000 mols): consider `--cnn_scoring none` first pass, then rescore top hits
- GPU acceleration: significant speedup with CUDA-enabled build

## Known Issues

| Issue | Workaround |
|-------|-----------|
| No macOS binary | Use Docker or Linux VM |
| PDBQT requirement dropped in 1.1 | Use PDB/SDF directly |
| Memory with large SDF | Split into batches of ~100 mols |
| Hydrogen handling | GNINA adds H internally; input can be with or without H |

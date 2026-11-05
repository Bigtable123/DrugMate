# Uni-Dock adapter notes

The Step 3 adapter targets the open-source Uni-Dock engine through
`unidocktools unidock_pipeline` in the dedicated `unidock` Conda environment.
Runtime preflight requires both the `unidocktools` front end and the underlying
`unidock` executable. The pipeline orchestrator automatically gives the outer
Step 3 process at least the enabled engine timeout plus 300 seconds; an explicit
top-level `step_timeouts.docking` value still takes precedence.

## Why the Tools wrapper is required

The raw engine accepts SDF in GPU-batch mode, but a normal SDF does not contain
Uni-Dock's `fragInfo`, `torsionInfo`, and `atomInfo` properties. Without those
properties the engine reports that fragment information is absent and performs
rigid-ligand docking. The Tools wrapper constructs the flexible topology before
launching the engine.

## Native output

The wrapper writes one SDF per input molecule and a native
lower-is-better `docking_score`. The adapter treats the SDF files, not
`results.csv`, as the source of truth because the wrapper may overwrite the CSV
between batches. It verifies that every parseable input molecule has at least
one finite scored pose before publishing canonical outputs.

## Executable screening policy

Treat `docking_score` as a lower-is-better Vina/Vinardo-family energy. The table
contains level-1 per-metric marks only; `pass`/`flag`/`fail` are folded with the
other category metrics before a molecule receives `keep`/`modify`/`reject`.
Apply the fallback gate only to finite values from one target/pocket cohort with
the same Uni-Dock build/runtime, scoring function, receptor and ligand
preparation, box, and material search settings.

| Metric | Pass mark | Flag mark | Fail mark | Suggested signed weight |
|---|---|---|---|---:|
| `docking_score` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -3 |

Normalize applicable scores with one declared method for the whole cohort,
orient lower-is-better values before weighting, compute
`sum(normalized_value * weight)`, and rank by descending total. Missing or
non-finite values are NA and do not enter normalization or the weighted total;
the adapter separately rejects an incomplete run whose input molecules do not
have finite poses.

The adapter publishes `engine_pose_rank = engine_pose_index + 1` as the
audit-only local mode rank. If two finite modes exist for a source molecule,
`gap_1_2 = score_2 - score_1` can be derived after sorting `docking_score`; a
gap below `0.5` kcal/mol marks ordering uncertainty and requests
review/modification, but never rejects the molecule. Rank and gap have weight
0; the current Uni-Dock adapter does not publish the derived gap.

The default request is `num_modes=20`, capped to 20/5/1 by cohort size; users may request fewer in the
experiment config.

Uni-Dock `docking_score`, GNINA `minimizedAffinity`, and AutoDock Vina
`vina_affinity` are homologous Vina-family evidence. Select one predeclared
primary family energy and do not sum them for the same molecule. Native
redocking is not implemented for the Uni-Dock adapter; its redocking status is
NA/not-applicable rather than pass. The current GNINA redocking route measures
best-of-generated-modes RMSD and is not a top-N active/decoy success rate.

## Current scale boundary

The wrapper creates many temporary and result files and prepares topology with
a multiprocessing pool. The present adapter is suitable for validated
shard-sized runs. Before claiming monolithic 100k–1M molecule support, add a
shared sharding/resume/failure-manifest layer and enforce CPU/GPU resource
limits per shard.

## Supported modes

- receptor: prepared PDB
- ligand library: 3D SDF
- scoring: `vina` or `vinardo`
- search: explicit `exhaustiveness` + `max_step`, or one `search_mode`
- box: explicit pocket box, or a reference-ligand bounding box with GNINA-like
  padding on both sides of each axis

AD4 scoring is intentionally not exposed by this adapter because its map/grid
dependency has not been validated in the registered runtime.

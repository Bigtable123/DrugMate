# AutoDock Vina adapter notes

The Step 4 registry ID is `autodock_vina`. This reference describes the classic
CPU adapter after it has been explicitly selected. Engine choice belongs only
to `../../pipeline-entry/TOOL_SELECTION_POLICY.md`; the runner never chooses
GNINA or Vina from artifact availability.

## Validated configuration contract

Step 4 must enable exactly one of `gnina`, `unidock`, `autodock_vina`,
`diffdock`, `dock6`, `autodock4`, `rdock`, and `adfr`.
An explicitly selected AutoDock Vina plan requires all of the following:

- enable `docking.docking_prep` with its required justifications, one supported
  final receptor backend, an explicit water policy and declared ligand states;
- enable `docking.meeko` in Step 4 to produce `receptor_pdbqt` from the current
  final prepared receptor;
- disable every other primary engine;
- disable `redocking_validation`, which is still a GNINA-specific protocol;
- set `docking.params.search_space_source` explicitly to `reference_ligand` or
  `docking_box` and supply the corresponding validated artifact.

The validator rejects an incomplete combination before execution. Step 4 runs
final receptor preparation and then Meeko before invoking Vina. The resulting
PDBQT is the engine input; supplying an older `receptor_pdbqt_file` does not
bypass this preparation. The adapter does not enable Meeko or infer a
search-space source.

## Runtime and data flow

The adapter deliberately separates two registered runtimes:

1. the `reduce` environment uses RDKit and Meeko to validate each 3D SDF
   record, prepare ligand PDBQT variants, and convert native docked PDBQT poses
   back to a lineage-preserving SDF;
2. the `vina` environment uses the AutoDock Vina Python API for CPU docking.

The engine computes the full Vina/Vinardo map set before loading ligands, then
reuses those maps across the library. Ligand preparation requires every input
record to succeed, and the engine stage requires every preparation variant to
succeed. A failure in either stage stops the batch before canonical pose export;
native files and diagnostics remain available for review.

Pose export handles per-record conversion failures separately. If at least one
usable scored pose remains, it publishes those poses, lists failed records in
`export_diagnostics.failures`, and records incomplete coverage in
`operational_errors`. Zero usable exported poses fails the batch. Exported poses
still pass the common docking identity, geometry and score-artifact checks before
downstream publication.

Supported current scope:

- rigid receptor PDBQT;
- 3D multi-record ligand SDF;
- Vina or Vinardo scoring;
- reference-ligand autobox or explicit normalized box;
- CPU execution controlled by `cpu`;
- per-run `exhaustiveness`, `num_modes`, `energy_range`, `min_rmsd`,
  `max_evals`, `spacing`, `force_even_voxels`, `no_refine`, `seed`, and
  `timeout_seconds`.

AD4 maps/scoring, flexible receptor docking, hydrated docking, GPU execution,
and engine-neutral native-ligand redocking are not exposed by this adapter.

## Native outputs and score semantics

Canonical poses are written to `docked_poses.sdf`; audited preparation and
native PDBQT files remain under `vina_native/`. Each pose retains its source
molecule, preparation variant, and engine pose index.

| Field | Meaning | Direction |
|---|---|---|
| `vina_affinity` | Vina/Vinardo total energy in kcal/mol | lower is better |
| `heavy_atom_count` | Heavy-atom count used by the derived efficiency | context only |
| `ligand_efficiency` | `-vina_affinity / heavy_atom_count` | higher is better |
| `engine_pose_rank` | One-based mode rank within one preparation | context only |
| `score_gap_1_2` | Second-best minus best affinity within one preparation | context only |
| `vina_inter` | intermolecular energy component | lower is better |
| `vina_intra` | intramolecular component | context only |
| `vina_torsions` | torsional component | context only |
| `vina_intra_best_pose` | best-pose intramolecular component | context only |

Do not compare or average these raw values with GNINA CNN scores or Uni-Dock
scores.

## Executable screening policy

Apply these fallback gates only to finite values from one target/pocket cohort
with the same Vina version, scoring function, receptor and ligand preparation,
box, and material search settings. These are level-1 per-metric `pass`/`flag`/
`fail` marks. The category fold, rather than an individual mark, determines the
molecule-level `keep`/`modify`/`reject` result.

| Metric | Pass mark | Flag mark | Fail mark | Suggested signed weight |
|---|---|---|---|---:|
| `vina_affinity` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -5 |
| ligand efficiency, `-vina_affinity / HeavyAtomCount` | `>= 0.25` | `[0.05, 0.25)` | `< 0.05` | +2 |

Ligand efficiency is derived, not a native Vina field. The adapter calculates
and publishes both `heavy_atom_count` and `ligand_efficiency` from the exact
lineage-matched exported pose. A missing score, missing heavy-atom count,
non-finite value, or a preparation with fewer than the poses needed for a
derived statistic is NA; do not replace it with zero or include it in scoring.

The adapter publishes `engine_pose_rank = engine_pose_index + 1` as an
audit-only rank within one Meeko preparation variant; the global `pose_index`
is not that rank. Rank `<=20` is inside the declared export window; a larger
rank is rejected from the pose-export set, while the parent molecule is mapped
to modify so an alternate pose or rerun can rescue it. When
two finite modes exist the adapter also publishes
`score_gap_1_2 = score_2 - score_1` within the same
`(source_molecule_index, preparation_variant_index)` group. A gap `>= 0.5`
kcal/mol needs no warning; a smaller gap marks the pose ordering as uncertain
and requests review/modification, but never rejects the molecule. Pose rank and
score gap both have weight 0.

The default request is `num_modes=20`, capped to 20/5/1 by cohort size; users may request fewer in the
experiment config.

Vina RMSD lower/upper bounds are not part of the current canonical score JSON.
Do not infer them from the energy array. If they are later parsed from native
PDBQT output, keep them as weight-0 mode-diversity diagnostics. Likewise,
`vina_inter`, `vina_intra`, `vina_torsions`, and `vina_intra_best_pose` have
weight 0 because `vina_affinity` already contains the native total-energy
evidence.

Normalize applicable non-zero-weight metrics with one declared method for the
whole cohort, orient lower-is-better energy before weighting, calculate
`sum(normalized_value * weight)`, and rank by descending total. AutoDock Vina,
Uni-Dock, and GNINA minimized energy are one Vina-family evidence class: choose
one predeclared primary family energy and never stack their weights for the same
molecule.

Native-ligand redocking remains unavailable for this adapter. The current
redocking implementation is GNINA best-of-generated-modes validation, not a
Vina adapter result and not top-N active/decoy success.

## Scientific role and downstream use

AutoDock Vina is useful as a reproducible classic CPU baseline and for users
who explicitly need the Vina Python implementation. It is not a strongly
orthogonal validation vote for GNINA or Uni-Dock: GNINA is Vina-derived and
Uni-Dock also implements the Vina/Vinardo method family. Agreement among them
therefore supports within-family reproducibility, not independent confirmation
of a binding pose.

Downstream pose-score normalization recognizes `vina_affinity` under its
AutoDock Vina namespace. Preserve its native identity and units; it must never
be copied into fake GNINA CNN fields.

## Official sources

- AutoDock Vina repository: https://github.com/ccsb-scripps/AutoDock-Vina
- AutoDock Vina Python API tutorial:
  https://autodock-vina.readthedocs.io/en/latest/docking_python.html
- Official basic docking and Meeko preparation tutorial:
  https://autodock-vina.readthedocs.io/en/latest/docking_basic.html
- AutoDock Vina 1.2 paper:
  https://doi.org/10.1021/acs.jcim.1c00203
- Original AutoDock Vina paper:
  https://doi.org/10.1002/jcc.21334
- Meeko repository: https://github.com/forlilab/Meeko

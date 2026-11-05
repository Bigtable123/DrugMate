# DiffDock adapter

## Execution contract

When explicitly selected with
`docking.params.search_space_source=whole_receptor`, DiffDock consumes the
prepared receptor PDB and each ligand record separately, performs blind
whole-receptor inference, and writes ranked native predictions under an
isolated `diffdock_native/run_*` directory. A Step 1 box or reference ligand is
retained as protocol context but does not condition this local DiffDock v1.1
implementation. The validator rejects another source choice; the adapter does
not infer or rewrite it.

The registered runtime should use the prefix environment
`${SBDD_ENV_DIFFDOCK_PREFIX}`, command `python`, and source tree
`${SBDD_DIFFDOCK_SOURCE_DIR}`. Both local checkpoints and both
`model_parameters.yml` files are mandatory; the adapter does not download
models during a pipeline run.

The default request is `samples_per_complex=20`, capped to 20/5/1 by cohort size; users may request fewer in
the experiment config. Report short or incomplete output against that request.

## Output contract

For every input molecule, the adapter requires exactly the requested
`samples_per_complex` files named `rankN_confidenceS.sdf`, with a finite native
confidence and contiguous ranks starting at one. DiffDock can exit zero after
logging a failed or skipped complex, so process status alone is never accepted
as success. Any missing molecule/rank, duplicate rank, unparsable SDF, missing
3D conformer, non-finite confidence, or canonical identity failure rejects the
whole run.

Each published pose contains:

- global `pose_index` and canonical `pose_id`;
- `source_molecule_index`, `engine_pose_index`, and `engine_pose_rank`;
- `docking_engine=diffdock` and `score_direction=higher`;
- the native `diffdock_confidence` value.

## Score semantics

`diffdock_confidence` is confidence-model evidence about pose plausibility. It
is not a binding affinity, has no kcal/mol or pK unit, and must not be compared
or summed with GNINA, Uni-Dock, AutoDock Vina, DOCK 6, RTMScore, or physics-based
energies. Rank only within a cohort with the same receptor, ligand preparation,
model/checkpoint identities, sampling parameters, and protocol ID.

For the `AFF` category this is level-1 audit evidence: a finite confidence is
recorded as `pass`, while missing/non-finite evidence is `not_applicable` and an
availability blocker. It has no universal absolute cutoff and suggested signed
weight `0`; category folding and V1/V2/V3, not this field alone, determine the
molecule classification.

## Reproducibility

The launcher seeds Python, NumPy, and Torch before importing DiffDock and pins
`CUDA_VISIBLE_DEVICES`. GPU kernels can still have nondeterministic behavior;
the score protocol therefore records seed, sampling parameters, source revision,
and SHA-256 identities for both model checkpoints and parameter files.

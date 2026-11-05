# RTMScore / RescoreRTMScore contract

RTMScore is an opt-in canonical Step 5 pose-rescoring backend. The executable
runtime uses the `${SBDD_ENV_RTMSCORE_PREFIX}` Conda prefix and must pass its
registered model and dependency checks before execution.

## Execution order and identity

1. PoseBusters/PoseCheck and interaction evaluation run first.
2. Every `pass`/`warn` pose in the evaluated cohort is copied to an isolated SDF
   with a deterministic upstream token while retaining its canonical `pose_id`.
3. The receptor is cropped by whole residues around the exact ligand cohort.
4. Each selected pinned model runs in the isolated RTMScore environment.
5. Every model must return exactly one finite score for every input token. A
   missing, duplicate, extra, or unknown identifier rejects the complete join;
   duplicate canonical `pose_id` values are rejected before model inference.
6. The parent joins only by canonical `pose_id`, validates all records before
   committing any update, and then reselects one pose per compound.
7. Candidate topologies, short MD, and MM-GBSA are built only after this
   optional reselection.

`rtmscore` runs the real model ensemble and publishes
`rtmscore.ensemble_score` (`direction=higher`). `rescore_rtmscore` is the
pipeline-internal control that optionally makes this metric the primary
post-QC rank; it is not a second model or Conda environment.

## Scientific semantics

- The value is a machine-learning pose score in an arbitrary RTMScore unit,
  not an experimental affinity or free energy.
- Compare it only within the recorded receptor, pocket crop, model ensemble,
  and protocol cohort.
- Never average or add it to Vina/GNINA/DOCK6 energies.
- Ensemble standard deviation is model-disagreement evidence, not a calibrated
  experimental uncertainty.
- The manifest records input/model hashes, parameters, resolved devices, and
  package versions. It deliberately marks the transitive protocol closure as
  incomplete because every imported upstream source/environment file is not
  content-addressed.
- Failure also publishes a registered manifest with the completed-model
  provenance captured so far; partial model scores are never rankable, and an
  enabled `rescore_rtmscore` is reported as dependency `not_evaluable`.

The portable default is `num_workers=0`; increase it only where forked data
loader workers are known to function. `device=auto` uses CUDA when exposed and
otherwise runs on CPU.

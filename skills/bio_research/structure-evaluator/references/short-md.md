# Short OpenMM MD and MDAnalysis contract

`openmm_md` and `mdanalysis` are independently selected canonical Step 5 tools
that share the registered `openmm` Conda environment. Selection belongs to
`../../pipeline-entry/TOOL_SELECTION_POLICY.md`; this reference defines their
execution contract. The runtime preflight creates a real CUDA Context; merely
importing OpenMM is not considered sufficient. The registered adapter currently
exposes only `platform=CUDA`; CPU/OpenCL require a separately validated runtime
contract instead of bypassing this preflight.

## Required topology protocol

- `candidate_topology` must be enabled in `build` mode with `fail_fast=false`.
- The currently validated protocol is ff14SB protein + GAFF2/AM1-BCC ligand,
  neutralized TIP3P explicit solvent in a periodic octahedral box.
- Generic metal, covalent, disconnected, radical, unsupported-element, or
  unparameterized cofactor systems remain `not_evaluable`; there is no silent
  generic fallback.
- Candidate topology capacity must cover the requested MD `top_k`.
- Every passed candidate must provide matching SHA-256 values for its solvated
  prmtop/inpcrd; the adapter re-hashes both files immediately before launch.

The short-MD cohort is bound before applicability filtering. If original rank
1 is not evaluable and `top_k=1`, rank 2 is not substituted. Every original
top-K rank remains visible in the manifest.

## Simulation protocol

- PME, 1.0 nm nonbonded cutoff, HBonds constraints, rigid TIP3P water,
  LangevinMiddle integration, and Monte Carlo pressure coupling.
- No hydrogen-mass repartitioning is used, so the timestep is capped at 2 fs.
- Integrator, velocities, and barostat receive the recorded deterministic seed.
- The production step/time origin is reset after equilibration, so reporter
  scheduling and expected-frame counts remain exact for non-default intervals.
- The exact solvated docked coordinates are saved before minimization and used
  as the analysis reference.
- State rows and final energies must be finite; expected reporter row counts and
  every material output hash are recorded.

This is a short pose-stability diagnostic. It is not an equilibrated production
ensemble and does not estimate binding free energy.

## Trajectory analysis

When `mdanalysis` is enabled, the trajectory uses no-jump periodic-boundary
reconstruction and protein-backbone alignment against the saved equilibrated
production-start reference. The separate pre-minimization docked structure is
retained as provenance but is not the COM-drift/RMSD baseline. It reports ligand
heavy-atom RMSD minimized over element/Amber-type topology automorphisms (with
no ligand-only coordinate refit), protein-backbone RMSD, ligand
center displacement, initial-contact retention, and ligand-RMSD drift. At least
two frames must remain after `analysis_stride`.

If the equilibrated production-start reference has no heavy-atom contact at the configured cutoff,
contact retention alone is `not_evaluable`; RMSD/displacement analysis remains
valid and this is not reported as an MDAnalysis runtime failure. Input and
analysis-output hashes, the no-jump/alignment policy, and the aggregate stability
report hash are retained in provenance.

The shared numeric bounds are applied as screening heuristics when
this optional tool is enabled. They are not universal experimental cutoffs;
the trajectory remains descriptive evidence for same-protocol comparison.

OpenMM simulation and MDAnalysis have separate candidate and aggregate states.
An analysis failure never rewrites a successful simulation as failed; the
trajectory remains available for diagnosis or reanalysis.

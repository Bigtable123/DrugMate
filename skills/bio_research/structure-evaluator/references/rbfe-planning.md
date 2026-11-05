# RBFE planning backends

`openfe_rbfe` and `ties_rbfe` are optional planning backends. When explicitly
enabled, they consume only the final post-QC `selected_poses.sdf`: at most one
accepted pose per compound, with canonical chemical/source/pose identity and
coordinate lineage intact.

They are planning tools in the current workflow. Neither backend runs
production alchemical simulations, and neither may report or imply a numerical
ΔΔG.

## Selection and assembly boundary

`../../pipeline-entry/TOOL_SELECTION_POLICY.md` is the only guide for choosing
OpenFE, TIES, both, or neither. `QUERY_TO_CONFIG.md` only records that decision
in the complete matrix. There is no `auto/openfe/ties/both` routing profile and
neither the validator nor runner chooses a backend or falls back to another
one. The validator checks cohort, chemistry, dependency, and artifact
requirements; the runner executes each accepted enabled backend. Each enabled
tool keeps its registered `preflight` parameter default unless the explicit
config requests `plan`.

The generic per-candidate Amber topology route is also selected explicitly by
enabling `pose_evaluation.candidate_topology` when MM-GBSA or short MD needs it. OpenFE
constructs its own OpenFF/OpenMM systems, while TIES requires pair-specific
hybrid/alchemical inputs; the ordinary `candidate_topology` artifact is not a
substitute for either.

## Shared hard gates

- At least two distinct, usable post-QC chemical identities are required.
- Canonical `pose_id`, `chemical_id`, `source_identity_id`, coordinates, and
  lineage are recomputed and checked before a node enters a plan.
- Each chemical contributes only its one final post-QC pose.
- Generic metal-site and covalent-complex routes are `not_evaluable`.
- Non-water receptor HETATM records require an explicit, reviewed
  parameterization route; the generic planners do not silently drop them.
- Disconnected ligands, radicals, unsupported elements, invalid 3D records, and
  identity mismatches cannot enter the network.
- TIES requires one formal-charge cohort and never permits a charge-changing
  edge in this route. OpenFE charge-changing mappings are disabled by default
  and require explicit expert acknowledgement plus a reviewed correction
  protocol.

## OpenFE

`preflight` loads the protein and ligand components in the registered `openfe`
Conda environment. `plan` additionally creates a ligand mapping network and
GraphML using Kartograf or LOMAP with a minimal-spanning,
minimal-redundant, or radial strategy. Mapping score, mapped-heavy-atom
fraction, network connectivity, input hashes, versions, and command provenance
are recorded in `openfe_rbfe/openfe_rbfe_manifest.json`.

`minimal_spanning` is the bounded planning default and has no cycle-closure
redundancy. Before any future execution campaign, review the mappings and
normally re-plan with `minimal_redundant` plus an explicit edge/resource budget.

This is deliberately a mapping-level plan. It does not call `openfe quickrun`,
does not generate a completed complex/solvent thermodynamic cycle, and always
records:

```json
{
  "simulation_executed": false,
  "free_energy_estimates_produced": false,
  "simulation_ready": false
}
```

The default contract requires explicit ligand hydrogens and assigned
stereochemistry. Relax either only after a reviewed upstream microstate policy.

## TIES

The runtime probe requires both roles:

- the `ties` builder API (`Pair`, `Config`, and `Pair.superimpose`) plus its
  CLI for pair mapping/hybrid input generation;
- the real `TIES_MD.TIES.TIES` runner API, importable analysis module,
  ABI-compatible OpenMM, and `ties_md`/`ties_ana` CLIs for later simulation
  and analysis;
- AmberTools `antechamber`, `parmchk2`, `tleap`, and `sqm` for builder-side
  ligand charge assignment and hybrid input preparation.

These components remain atomic in the registered `ties` environment. On the
validated host, `drugmate` uses Python 3.11 while `ties` uses Python 3.10, and
each environment has its own compiled OpenMM extension. Never splice the
`drugmate` OpenMM package into TIES with `PYTHONPATH`: TIES_MD must use the
ABI-compatible OpenMM installed beside it. The `drugmate` OpenMM installation is
the selected runtime for a future standalone short-MD backend; the current
pipeline prepares candidate topologies for that route but does not yet launch
short MD.

`simulation_platform=auto` performs a real tiny-context OpenMM calculation in
CUDA, OpenCL, then CPU order and records the first functional TIES_MD platform.
An importable OpenMM package is not sufficient evidence that CUDA or any other
platform can execute. An explicit platform value probes only that platform and
never silently falls back. OpenMM's Reference platform is intentionally not a
TIES execution route.

`preflight` checks runtime, identity, and applicability. `plan` writes a
deterministic star network with pending stages for TIES hybrid setup, complex
and ligand ensemble simulations, and thermodynamic-cycle analysis. The output
is `ties_rbfe/ties_rbfe_plan.json`.

The current wrapper does not guess an unverified TIES API or manufacture hybrid
topologies when the registered environment is incomplete. It returns
`not_evaluable/ties_runtime_unavailable` with probe evidence instead.

## Evidence semantics

Backend and candidate/transformation records use:

- `passed`: the declared preflight or planning contract completed;
- `not_evaluable`: the chemistry, cohort, or runtime is outside the declared
  route;
- `failed`: an expected operation was attempted but its contract failed.

A passed plan is not a binding-affinity result. Do not insert a plan, mapping
score, or future pairwise ΔΔG into the pose-level `pose_scores` contract.
Future analyzed RBFE needs a separate edge/compound-level artifact containing
both thermodynamic legs, repeat-level results, uncertainty, convergence,
cycle-closure diagnostics, method/version, and provenance. OpenFE and TIES
outputs remain separate even when both are requested; never average them merely
because their final units match.

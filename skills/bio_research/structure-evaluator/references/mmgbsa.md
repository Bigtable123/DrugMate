# Amber MM-GBSA escalation

This reference defines execution and scientific semantics after `mmgbsa` has
been selected explicitly.

## Selection and validation boundary

`../../pipeline-entry/TOOL_SELECTION_POLICY.md` decides whether
`candidate_topology` and `mmgbsa` are enabled; `QUERY_TO_CONFIG.md` only
assembles the explicit decision. The validator requires
`candidate_topology.execution_mode=build` and
`candidate_topology.max_candidates >= mmgbsa.top_k`. It does not require the
Step 0 reference-complex `amber_topology`: execution consumes only
candidate-specific topologies built from the post-QC survivor SDF. The runner
does not add this dependency or enable either tool implicitly.

## Data and execution contract

`mmgbsa` reads `candidate_topologies.json`, never the original docked SDF. Join
every result by canonical `pose_id`. Validate the three Amber topology files
and coordinate files before execution and record their hashes.

The registered environment id is `mmgbsa`, mapped on this host to the Conda
environment named `gmxMMPBSA`. The validated implementation is AmberTools
`MMPBSA.py`; the separately distributed `gmx_MMPBSA` package and CLI are not
installed and must not be claimed in provenance.

`preflight` validates the runtime and topology bundle without calculating an
energy. `score` runs `MMPBSA.py` separately for up to `top_k` eligible poses and
writes `mmgbsa/mmgbsa_manifest.json`.

## Scientific semantics

The current backend uses one `complex.inpcrd` snapshot per pose with one
consistent GB model and salt concentration. It does not:

- sample an MD ensemble;
- estimate configurational entropy;
- minimize or alter the post-QC coordinates;
- produce an experimental or absolute binding free energy;
- combine its value numerically with GNINA/Vina scores or RBFE results.

Expose the primary value as
`amber_mmpbsa.delta_total_no_entropy`, in kcal/mol with lower-is-better
direction and `low` reliability. It is a within-target, within-protocol
endpoint reranking signal. When `use_as_primary_metric=true`, use it as the
single primary metric for the scored cohort; retain all other scores as
independent evidence.

Do not add an implicit minimization step. A coordinate-changing relaxation
would require a child `pose_id`/coordinate revision and renewed pose QC.
Future trajectory-based MM-GBSA should be a separate execution mode with frame
counts and uncertainty, not silently substituted for this single-frame
contract.

## Status handling

- `passed`: at least one requested score was produced and no operational
  failure occurred.
- `not_evaluable`: no scientifically eligible topology exists; the manifest is
  still a successful diagnostic artifact.
- `failed`: runtime missing, manifest/identity corruption, topology/coordinate
  inconsistency, timeout, non-zero `MMPBSA.py` exit, or unparseable/non-finite
  output. An explicitly enabled MM-GBSA operational failure fails Step 4.

OpenMM is not required for this backend. Its presence in `drugmate`, `openfe`, or
`ties` does not make it available inside the isolated MM-GBSA environment.

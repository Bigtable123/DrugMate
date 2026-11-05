# Receptor Preparation and Pre-Docking Preparation

This workflow follows the preparation-stage split of Stage 02 / Stage 06 in the
delivery specification. The current registry, generated catalog, and validator
govern configuration. Do not copy a tool matrix directly from an old experiment.

## Structure Preparation: Establish Only the Geometric Baseline

Even when the user supplies a PDB ID, first verify the target entity, author
chain, assembly, model, and reference-ligand instance. `context.chain` has no
default. For example, CDK6 is chain B in 2EUF; do not select A merely because an
example uses A. Retain partners required to form the pocket. For complex
assemblies, first supply a reviewed local structure. Explicitly choose one model
from a multi-model input before execution; the runner does not choose for you.

| Setting | Current default/requirement |
|---|---|
| `structure_prep.pdbfixer.add_missing_atoms` | `true`; only for justified repair of standard heavy atoms |
| `add_missing_residues` | `false`; enabling requires `repair_reason`; only internal short gaps of at most 3 residues are accepted, while terminal and long gaps are rejected |
| `add_hydrogens` | Must be `false` |
| `seed` | `0`, for reproducible heavy-atom repair |
| `structure_prep.propka` | Disabled by default; explicitly enable when needed, such as for pH-sensitive questions, and provide a diagnostic pH. Its report does not establish that protonation states have been rewritten |
| `structure_prep.reduce` | Disabled by default; enable only when pocket-detection/interaction tools need explicit H, and supply `params.reason` |

The structure stage produces `prepared_receptor_noH`, `prep_manifest`, basic QC,
the original structure, and `water_inventory`. The inventory retains original
water coordinates; the receptor's geometric baseline does not preselect the
final waters. Optional Reduce produces a separate `provisional_receptor` and a
`protonation_manifest` marked provisional. It neither overwrites the hydrogen-free
baseline nor becomes the final docking receptor.

For ProLIF fingerprints, have `pocket_analysis.prolif` consume the hydrogenated
`provisional_receptor`: enable `structure_prep.reduce`, set `build=true`, and
provide the interaction-analysis rationale, or inject a reviewed hydrogenated
provisional structure. Without a reference ligand, ProLIF records `not_evaluable`
and skips analysis. Without suitable protein hydrogens, it records an analysis
failure. Neither condition blocks other evidence or pocket results. Water/ligand
hydrogens or a `prepared_with_H` label do not establish that the receptor has
been hydrogenated.

Checks before and after preparation are the existing basic file/atom-record
checks. **Enhanced valence, clash, and reconstructed-geometry QC has not yet been
implemented.** Do not present `passed=true` as proof of reliable geometry.
Suspicious added-atom results still require review. Do not hide this implementation
limitation with `allow_bad_residues=true`.

## Step 2: Standardize the Library Without Prejudging Nomination

The current `molecule_generation` corresponds to Stage 04 in the delivery
specification. It acquires sources, standardizes structures, deduplicates by the
exact standardized graph, preserves compatible native coordinates or generates
ETKDG conformers, and delivers `library` and `library_manifest`. Step 2 does not
require pH microstate enumeration or randomly replace a structure with another
protonation state or stereoisomer. All libraries are handled using actual molecular structures and
source/conformer identities.

Preserve source, role, parent-molecule, and conformer identities.
`microstate_id` may be empty when states have not been expanded. Check actual
structures, consistency between declared stereochemistry and coordinates, usable
finite 3D coordinates, and artifact counts. Isolate invalid records and preserve
their failure sources and reasons. Multiple possible states or a missing state
abundance model are review information, not automatic prohibitions on candidate
ranking or nomination. `DEC_*` entries are never nominated; final eligibility
still depends on subsequent docking, pose evaluation, properties, and this
round's protocol evidence.

The old state scripts prepare chemical representations and check consistency;
they are not drug-property filters. Step 2 retains necessary parsing,
standardization, identity, role, and usable-3D checks. Lipinski, QED, SA, PAINS,
and other drug-property checks remain in subsequent `property_prefilter` /
`property_filter`, following this round's tool configuration.

## Docking: Complete docking_prep Before Invoking an Engine

The order is structural geometry preparation → pocket selection → ligand
library/prefiltering → **docking.docking_prep** → selected docking engine.
`docking_prep` is a required substage within docking, not a new top-level step.
Its tool choices, preflight across Conda runtimes, parameters, input SHA, results,
and timeouts enter the current plan, reports, and cache fingerprints.

The agent must explicitly populate `steps.docking.tools.docking_prep`:

- `enabled=true`.
- `params.reason`: justification for the selected final protonation/flip backend.
- `params.ph`: required when using PDB2PQR or ligand-microstate diagnostics;
  the delivery specification does not prescribe a universal pH of 7.4.
- `params.water_policy`: `dry` or `selected`, with `water_reason`.
- `params.ligand_state_policy="preserve_declared"` and `ligand_state_reason`:
  explain why the ligand representation after prefiltering is retained and its
  associated uncertainty. Ordinary non-enumerated libraries may also use this
  route. Preserve their identities and record order; retaining a representation
  does not establish the dominant microstate.

Select exactly one final-receptor backend:

- `docking.reduce`: **preferred for ordinary noncovalent protein receptors when
  protonation-state assignment at a target pH is not required**. Explicitly set
  `enabled=true`, `build=true`, and `flip=true`, and pass runtime preflight.
  It adds hydrogens and optimizes His/Asn/Gln flips; **do not claim that it applies
  a target pH**.
- `docking.pdb2pqr`: choose when protonation states must be assigned at a target
  pH, for a protonation-sensitive protocol, or when the user explicitly specifies
  PDB2PQR + PROPKA. Explicitly set `enabled=true` and an experimentally justified
  `docking_prep.params.ph`, and retain `optimize_hbond_network=true`. It adds
  hydrogens, assigns protonation states by pH, and optimizes hydrogen bonds and
  His/Asn/Gln flips. This backend calls PROPKA internally; there is no need to
  additionally enable the diagnostic `structure_prep.propka` for this purpose.
- `docking.protoss`: the current adapter can submit only a PDB ID and cannot
  preserve the prepared-complex identity, so validation rejects this selection.
  It can be enabled only after prepared-complex upload is implemented. Do not
  automatically download a different original structure to replace the input,
  or automatically switch backends.

These priorities guide the agent when assembling new configs. Preserve
user-specified backends and existing experimental protocols. If retained metals,
cofactors, or nonstandard residues are present, first check tool and force-field
support for those components. Do not delete necessary components to use a
preferred backend. If a protocol requires a target pH but lacks a justified
value, obtain that justification first. Do not universally fill in 7.4 or switch
to a route without target-pH assignment merely because Reduce is available.
The config must remain explicit with exactly one backend enabled. If the chosen
backend fails, report the reason; do not silently switch tools and continue.

The `selected` water policy requires `docking.water_selection` and the original
`water_inventory`. With a reference ligand, use the existing
distance/occupancy/B-factor rules. Without one, explicit water residues are
required. Preparation fails if the backend loses selected waters. For `dry`,
state the reason explicitly and do not run the water selector.

For Vina, enable `docking.meeko` to make PDBQT after the final receptor is
generated. The GNINA PDB/SDF route does not enable Meeko. UniDock's native adapter
handles its required format conversion. Do not inject an old `receptor_pdbqt`
to bypass this round's final-receptor preparation.

`docking.ligand_microstates` is an optional **reference-ligand** microstate
diagnostic. It must not silently replace the candidate library or crystal
reference conformation. The delivery specification places justified candidate
microstate expansion in docking_prep, but the current adapter implements only
`preserve_declared`; automatic candidate-library enumeration is not registered.
When a protocol requires additional states, document the input-preparation
rationale, use the existing library-injection route to pass through prefiltering
again, and preserve the microstate-to-parent mapping. Old state seals are no
longer read, and old state labels are not nomination rules. Actual structures,
sources, and coordinates are still checked as usual. This implementation does
not automatically infer state abundances; non-unique states alone do not veto
nomination. Amber topology remains part of per-candidate MM preparation in
`pose_evaluation.candidate_topology`. Do not prepare a generic topology for the
crystal reference ligand in the structure stage. Keep `docking.openbabel=false`
when the current route has no Open Babel consumer.

Before docking, produce separate `prepared_receptor`, optional
`receptor_pdbqt`, `receptor_state_manifest`, `water_report`,
`docking_species_library`, `microstate_identity_map`, and `docking_prep_qc`.
Subsequent pose evaluation consumes this final receptor. Stop on preparation
failure without invoking an engine. GNINA defaults to `exhaustiveness=16` and
`autobox_add=4.0`. Limit requested poses per molecule according to the number
of molecules in this round's `docking_species_library`: at most 20 for ≤50,
at most 5 for 51–1000, and at most 1 for >1000. Explicit configuration may further
tighten the cap.

PoseCheck directly reuses this round's hydrogenated `prepared_receptor` from
docking_prep. The ordinary route retains Reduce's hydrogenation and flip results;
a PDB2PQR protocol retains that backend's results. Do not switch backends during
evaluation. The local adapter uses RDKit to read existing atoms and coordinates,
preserving explicit H and negative residue numbers. It no longer calls Hydride
for hydrogen addition/optimization or converts the receptor through ProLIF.
If no readable atoms, coordinates, or explicit H are present, record a PoseCheck
preparation error. Do not substitute the hydrogen-free baseline or provisional
receptor. All-atom clashes remain diagnostic, while heavy-atom clashes continue
to govern gating. Ligand-strain reference-conformer sampling parameters are
unchanged.

The pocket-stage Blueprint remains pre-preparation evidence, with its manifest
requiring a final interaction review. No capability to automatically generate
`interaction_blueprint_final` has been added. Evaluate interactions of the final
receptor according to pose-evaluation tools and evidence rules; a provisional
state is not evidence of the final state.

## Entrypoints, External Inputs, and Migration of Old Experiments

1. Reassemble a complete explicit tool matrix from the current catalog, including
   all newly added docking-preparation tools.
2. Remove the old `water_selection/pdb2pqr/protoss/meeko/
   ligand_microstates/openbabel/amber_topology` keys from `structure_prep.tools`;
   assign them to their new stages as described above.
3. Set an explicit target chain, preparation rationale, and water/ligand-state
   policies; run `validate_config.py`.
4. After input-normalization, an external `prepared` receptor retains its original
   file and produces separate `prepared_receptor_noH`, water inventory, and prep
   manifest. This declares prepared geometry; final docking_prep is still
   required. For evaluation of existing poses only, the final-receptor input
   remains `prepared_receptor`.
5. Run in a new round directory without modifying existing experimental artifacts.
   Code and cross-stage dependencies enter the fingerprint. Old preparation
   results cannot be reused as artifacts of the new contract merely because
   their files exist.

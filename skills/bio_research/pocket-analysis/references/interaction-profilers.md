# Interaction profiler execution and runtime placement

Use this reference when changing the Step 1 interaction-profiler stack. Keep
pocket detection and interaction analysis in the existing single Step 1. The
Agent selects profilers only through
`../../pipeline-entry/TOOL_SELECTION_POLICY.md`; this file defines evidence and
execution behavior after those switches are explicit.

## Evidence policy

- Record PLIP and ProLIF agreement as an audit observation. ProLIF never
  contributes to Blueprint constraints or scientific gates. Eligible PLIP,
  LUNA and chemistry-checked Arpeggio observations may supply automatic
  Desired or Optional evidence.
- Treat profiler output as an observation, not a biological requirement.
- Normalize protein/ligand atoms or chemical features, direction, mediator and
  interaction-specific geometry before comparing tools.
- Count `not_evaluable` and `tool_failed` separately from a true negative.
- Compile a hard constraint only from reviewed evidence with a source,
  rationale and binding-mode scope. Multi-tool agreement alone is not causal
  evidence.
- Apply the same normalized matcher, or an explicit capability matrix, to
  reference structures and docked poses.
- Keep downstream Blueprint metric bands aligned with the executable policy: desired
  ranking `<0.10`/`[0.10,0.50)`/`>=0.50` maps to level-1
  `fail`/`flag`/`pass`; desired evaluable coverage has only a `0.80` soft lower
  bound, with `<0.30` stopping trusted ranking and requesting re-evaluation.
  Coverage alone never rejects a molecule or pose.

## Runtime matrix

| Tool | Runtime placement | Pipeline role | Notes |
|---|---|---|---|
| PLIP | `drugmate` | Core static profiler and pose evaluator | Keep in an isolated subprocess; retain atom IDs, direction and geometry. |
| ProLIF | `drugmate` | Record-only static/multi-pose/ensemble profiler | Read rich `Fingerprint.ifp` metadata rather than only the boolean DataFrame. |
| PDBe Arpeggio 1.4.4 | `drugmate` | Conditional atom-level adjudicator | Invoke by CLI subprocess. Accept only the Step 0 source complex mmCIF and consume JSON. Do not treat its OpenBabel-based output as an independent strict-consensus vote. |
| LUNA 0.14.0 | dedicated `luna` environment | Opt-in atom-aware interactions and hashed/count IFP | Invoke the core Python API with `out_pse=False`; do not use the PyMOL-dependent CLI for headless Step 1. |
| ODDT 0.7 | dedicated `oddt` environment | Legacy residue IFP/SPLIF comparison | Use the tested OpenBabel backend with Python 3.9 and NumPy 1.26. Interpret IFP values as counts and SPLIF only as a relative similarity/ranking feature. |

LigPlot+ is deliberately excluded: do not install it, register a tool ID, emit
an artifact for it, or add it to the tool matrix. Its GUI/report role does not improve
the structured machine-readable evidence contract used by this pipeline.

## Blueprint eligibility

- Keep ProLIF observations and agreement groups in the audit artifact with
  `blueprint_eligible=false`. They do not create `candidate_required`, Desired,
  Optional, or Required constraints.
- Permit atom-aware LUNA observations and chemistry-QC-passing Arpeggio
  observations to support Desired constraints when the selected consensus
  pocket contains the residue. Do not promote them to Required.
- Mark LUNA weak H-bonds, degraded Arpeggio observations and every ODDT record
  `blueprint_eligible=false`. Preserve them in the evidence artifact.
- Ignore ODDT's metal channel for interpretation because the upstream
  direction/index semantics are unreliable in this supported legacy version.

## Failure semantics

- Require explicitly enabled gate-bearing tools to pass pipeline runtime
  preflight. ProLIF is record-only: report a missing runtime under
  `optional_missing` and continue Step 1; its worker result records the fault.
- Record a ProLIF worker failure in `profiler_failures` and the failure
  registry without invalidating a sound Step 1 contract. An enabled PLIP
  worker failure currently invalidates it. `fail_on_profiler_error` remains
  a compatibility field and does not change this behavior.
- After preflight passes, record structure-specific Arpeggio/LUNA/ODDT worker
  failures in `supplemental_profiler_failures`. They currently invalidate the
  Step 1 contract; `fail_on_supplemental_profiler_error` is a compatibility
  field and does not change this behavior.
- Report a missing source mmCIF for Arpeggio as `not_evaluable`, not
  `tool_failed`. Report a profiler returning no matching contact as
  `not_detected`; never collapse either state into an empty successful result.

## Execution prerequisites

An optional profiler selected by policy must satisfy the relevant prerequisites
below before scientific execution:

1. Resolve dependencies without replacing Drugmate' RDKit, OpenBabel, NumPy or
   MDAnalysis packages.
2. Probe imports/CLI and record exact versions in `tool_runtime_report.json`.
3. Run a chemically complete reference-complex golden set covering directional
   H-bonds, salt bridges, water bridges, halogens, metals and insertion codes.
4. Verify canonical observation conversion and disagreement diagnostics.
5. Verify that Step 5 can match every constraint type, or marks it
   `unmatchable`; absence of evaluator capability must never become a miss.

For Arpeggio, a PDB file renamed or naively converted to mmCIF is insufficient:
preserve ligand bond orders/formal charges in `_chem_comp` using CCD or another
reviewed chemical template. Capture stdout and stderr. Treat OpenBabel
kekulization/aromaticity warnings as degraded chemistry quality and exclude the
affected observations from Blueprint compilation.

For LUNA, use the registered `luna` Conda environment and its core API. Preserve
`InteractionType.params`, because JSON summaries alone omit required geometry.
Map ligand atom serial numbers to the 1-based SDF order. If the ligand file has
no explicit hydrogens, report ligand-donor H-bond capability as not evaluable.
Maintain a golden hash/fingerprint regression because hashing dependencies can
change across environments.

For ODDT, set `ODDT_TOOLKIT=ob` before importing ODDT and map the returned
zero-based residue indices through the loaded protein residue list. Do not add
hydrogens in this wrapper. Do not expose SPLIF as a human-readable hard
interaction constraint; use it only for relative pose/series comparison.

Official sources:

- PLIP: <https://github.com/pharmai/plip>
- ProLIF: <https://prolif.readthedocs.io/en/latest/>
- Arpeggio: <https://github.com/PDBeurope/arpeggio>
- LUNA environment: <https://github.com/keiserlab/LUNA/blob/master/luna-env.yml>
- ODDT: <https://github.com/oddt/oddt>

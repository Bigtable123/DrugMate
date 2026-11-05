---
name: molecule-generator
description: Build a standardized 3D library for SBDD screening.
---

# Molecule Generator

Build the Step 2 candidate library from one or more registered sources. The
Agent chooses sources with `../pipeline-entry/TOOL_SELECTION_POLICY.md` and
writes every source switch explicitly through
`../pipeline-entry/QUERY_TO_CONFIG.md`. The validator rejects an empty or
incompatible source set; this skill and the runner only execute the accepted
switches. Infrastructure paths, Conda environments, checkpoints, and checksums
come only from `shared/pipeline_registry.py`.

## Step 2 boundary and nomination

Step 2 implements delivery Stage 04 `04-molecule-source-generation`: acquire
sources, standardize, deduplicate, prepare coordinates, and publish `library`
plus `library_manifest`. Its successful output is ready for property prefilter
and docking. Final nomination uses downstream pose, property, and campaign
evidence; molecular-state uniqueness is not a nomination prerequisite.

Standardization does not enumerate campaign-pH protomers, tautomers, or
stereoisomers. Preserve the supplied standardized representation and its declared
stereochemistry; record unresolved chemistry as review context. Do not replace
it with a randomly sampled state, require a population prior, or reject a valid
molecule merely because another state may exist. Libraries use the molecular
graph and source/conformer identity.

Step 2 performs parse/standardization, role/identity, and usable-3D
checks. Lipinski, QED, SA, PAINS, and related property assessment run in
`property_prefilter` and `property_filter` according to their registered plans.

When a campaign requires additional states, use an explicitly justified input
and an implemented docking-preparation route. The current `docking_prep` adapter
supports `preserve_declared`; its optional `ligand_microstates` tool diagnoses
the reference ligand only. It does not enumerate the candidate library. Supply
additional prepared candidate states through the existing library-injection
route and retain their parent mapping; do not invent a candidate-enumeration
tool or claim a dominant state. Existing input files are handled as molecular
records with source and coordinate lineage, without interpreting old state
tags as nomination rules or requiring a separate state manifest.

Keep `DEC_*` records benchmark-only throughout the pipeline. They never enter
candidate nomination or prospective ranking. Preparation success and favorable
generator scores do not establish binding activity.

## Selection boundary

Before choosing generation sources, complete the
[unified memory consultation](../pipeline-entry/references/benchmark-memory.md)
and applicable [generation selection](../pipeline-entry/references/generation-memory.md).
Use the five independent metric recommendations and record the actual source
choice and total candidate budget. An existing full library or a task without
a legal conditioning ligand does not become a generation task just to enable
memory. A recommendation cannot supply missing DiffLinker fragments, replace
the user's requested source, or override downstream property/pose gates.

Tool recommendations, user-intent mapping, and the no-preference decision live
only in the generated `TOOL_SELECTION_POLICY.md`. This skill does not enable
ChEMBL, Pocket2Mol, or any other source because a preference is absent, and it
does not silently add every installed generator. Every registered source has an
explicit `enabled=true|false` in a full pipeline config.

Several explicitly enabled sources may run together. They run independently and
are merged afterward. One source producing zero valid molecules does not
discard successful outputs from other selected sources. If every selected
source fails or yields zero valid molecules, Step 2 fails with per-source
diagnostics.

These tools are candidate sources, but they are not equivalent:

- ChEMBL retrieves experimentally reported actives.
- Pocket2Mol, ResGen, DiffSBDD, TargetDiff, Lingo3DMol, TamGen, and MolCRAFT
  generate molecules conditioned on a protein pocket. TamGen emits SMILES;
  the others can provide native 3D coordinates.
- DiffLinker elaborates supplied fragments; it is not pocket-only de novo.
- REINVENT4 currently samples 2D SMILES without pocket conditioning. This
  installation does not include Mol2Mol/LibInvent/LinkInvent priors, so this
  Step 2 REINVENT4 adapter cannot perform lead modification.

### ChEMBL unavailable: Pocket2Mol recovery

For candidate discovery, follow the ChEMBL fallback policy in
`../pipeline-entry/TOOL_SELECTION_POLICY.md`: when ChEMBL retrieval fails,
use Pocket2Mol for the candidate budget that ChEMBL was meant to supply.
This is an Agent config change, not an implicit runner switch. Respect an
explicit ChEMBL-only/known-active requirement or exclusion of generation;
a generated library cannot fulfill a retrospective active-control request.

- ChEMBL target/activity requests wait up to 60 seconds per socket operation;
  this is not a one-minute limit for the full paginated retrieval. The separate
  status probe remains 10 seconds. A status-probe failure alone does not prove
  target/activity endpoints are unavailable: a supplied target ID can still
  be used, and name/UniProt resolution attempts the actual target endpoint.
- Target lookup/search and every activity page use an explicit `.json` path.
  Each API read makes at most three total attempts (initial + two retries),
  with 1s/2s backoff for HTTP 408/429/5xx, URL/connection errors, socket read
  timeouts, and incomplete HTTP bodies. Other 4xx and invalid JSON fail
  immediately. Failed attempts are logged; exhausted retries preserve the
  source error and follow the recovery policy below. The lightweight status
  probe remains a single request.
- On a retrieval timeout, connection error or HTTP service failure, preserve
  the source error and raw capture before changing the plan. Do not repeatedly
  restart the unchanged ChEMBL-only pipeline. Do not treat an incomplete page
  set as a complete active library or label generated molecules as `ACT_*`.
- If usable generated/local candidates already exist, continue with them;
  do not regenerate merely to replace a missing control. Otherwise set
  `steps.molecule_generation.tools.chembl.enabled=false` and
  `steps.molecule_generation.tools.pocket2mol.enabled=true`. For a ChEMBL-only
  plan this makes Pocket2Mol the sole candidate source. Keep standardization
  and conformer preparation enabled. When docking is planned, enable both
  `steps.pose_evaluation.tools.posebusters` and `.posecheck` and retain the
  remaining downstream QC.
- Set `steps.molecule_generation.tools.pocket2mol.params.num_samples` from
  the requested total candidate budget, including any existing Pocket2Mol
  allocation. For a 60-candidate task with no retained candidates, request
  60 Pocket2Mol samples, not the default 100 or 60 from each source. If there
  is no separate total budget, transfer the planned ChEMBL `max_actives` plus
  any already planned Pocket2Mol sample allocation. Report actual usable
  counts after deduplication/3D QC; do not loop to replenish rejected records.
- Check the prepared receptor, pocket center, registered Pocket2Mol runtime,
  checkpoint and CUDA prerequisites. Revalidate the complete config. A source
  change in an already started campaign needs a new same-stage attempt under
  `../pipeline-entry/references/campaign.md`; preserve the failed attempt and
  reuse verified Step 0/1 artifacts through the existing reuse contract. Plan
  for this one fallback within the authorized attempt/resource budget; do not
  ask again when existing authorization covers it. If Pocket2Mol also fails,
  retain diagnostics and stop instead of cycling back to ChEMBL.
- Generated candidates remain `GEN_*`. For a prospective, generated-only
  fallback without ACT/DEC evidence, set
  `steps.coordination.tools.screening_benchmark.enabled=false` and disable
  `dude_decoys`/`rdkit_decoys` in the revised molecule-generation config.
  Report AUROC/EF as unavailable and continue prospective screening; never
  report a passed benchmark. This does not waive an explicitly requested
  retrospective benchmark. Preserve a usable co-crystal reference ligand
  for pocket analysis and redocking; that reference is independent of
  downloading ChEMBL actives.

### Decoy acquisition

Validation decoys support retrospective ACT/DEC screening benchmarks; they are
not nomination candidates. Prefer the selected DUD-E download/cache source for
covered targets. Local RDKit BRICS generation is optional and defaults to off;
enable it only when the user explicitly requests local decoy generation.
Unavailable downloads or isolated invalid DUD-E molecules are not a reason to
enable RDKit automatically. Preserve individual failures under the isolation
contract below and continue with the usable downloaded records.

On download failure, record the reason and continue the candidate-only workflow;
do not retry the experiment just to obtain decoys or enable local generation to
satisfy an AUROC/EF gate. When creating or revising a plan, do not carry forward
an old `rdkit_decoys.enabled=true` unless the current user request explicitly
asks for local decoy generation.

If no validation decoys remain, a candidate-only library can still proceed to
docking. Report ACT/DEC AUROC and enrichment as unavailable, never as passed.
Downloading precomputed decoys avoids local recombination but still requires
the same molecular-state and coordinate preparation as other library inputs.

## Inputs

The orchestrator wires files automatically from Step 0/1. Agent-authored config
contains only context, tool switches, and scientific parameters.

| Input | Used by |
|---|---|
| `target_name` / `chembl_target_id` / `uniprot_id` | ChEMBL; at least one is required when selected |
| `pocket_pdb_file` | All pocket-conditioned generators |
| `pocket_center` | Pocket2Mol, ResGen, DiffSBDD, TargetDiff, Lingo3DMol, TamGen |
| `reference_ligand_file` | ResGen and MolCRAFT require it; DiffSBDD uses it when available |
| `pocket_definition_file` | Provenance and future protocol extensions |
| `interaction_blueprint_file` | Pocket2Mol pocket trimming and pharmacophore annotations |
| `fragment_file` | DiffLinker only; comes from `context.fragment_file` |
| `seed` | All stochastic generators |

Tool parameters and defaults are generated into
`../pipeline-entry/tool_catalog.json`; do not duplicate runtime paths in the
user config.

Do not copy a partial source snippet into a pipeline config. Assemble the full
registered step/tool matrix as specified by `QUERY_TO_CONFIG.md`, then let
`validate_config.py` validate and materialize it before execution.

## Execution contract

1. Execute each explicitly enabled candidate source.
2. Run registered runtime preflight before scientific execution.
3. Verify the fixed model/prior SHA-256 before loading it.
4. Run selected generators serially through isolated environment-local workers
   and a new per-invocation work directory, so retries cannot ingest stale
   SDF/MOL/CSV files or overwrite fixed upstream `output/` directories.
5. Parse actual outputs; exit code zero alone is not success.
6. Reject unparseable, disconnected, non-finite, or missing-3D outputs from 3D
   generators. REINVENT4 and TamGen outputs are validated as 2D SMILES.
7. Standardize all candidate/decoy SMILES.
8. Deduplicate source records within each role using the exact standardized
   isomeric graph, preserving charge, tautomer, isotope, and stereochemical
   distinctions. Prefer a generated record with compatible native 3D
   coordinates and retain `source_tools` /
   `also_generated_by`. Preserve source structures and identifiers in
   `parent_input_records_json` (SDF/CSV) and the report's `parent_deduplication`.
   ChEMBL measurements retain their existing censoring/aggregation semantics.
   Collisions between roles remain explicit errors; never relabel a validation
   decoy as an active or generated candidate.
9. Apply round exclusions (`exclude_scaffolds`, `max_heavy_atoms`, `max_mw`) to
   ChEMBL and every generated candidate, never to validation decoys.
10. Preserve native coordinates only when connectivity, charge, atom mapping,
    and declared stereochemistry still match the standardized molecule.
    Otherwise generate ETKDGv3 + MMFF/UFF coordinates.
11. Check parseable chemistry, finite usable 3D coordinates, stable IDs and
    source/role/parent mapping. Write library and metadata only for molecules
    actually written to SDF, and reconcile input, excluded, and output counts.

### Molecule failures and tool faults

Step 2 isolates individual invalid molecular inputs, incompatible declared
stereochemistry, and failed conformer generation. Record the input/source identity,
SMILES or SDF file/record index, operation, error code, and reason in
`molecule_failure_report.json`; update it as failures occur. Such records have
`status=error`, not an inactive label or a property-filter rejection. Other
parents continue through preparation and docking-library assembly.

Release valid records independently. A failed record does not invalidate valid
sibling structures; retain failed-state provenance and report partial coverage
when an explicit expanded input contains siblings. Never claim an incomplete
state set is exhaustive. Preserve coordinate failure evidence in the failure
audit. Step 2 fails only when no usable prospective candidate remains or a
required source/runtime or shared input-integrity failure prevents execution;
it does not count final nomination eligibility at this stage.

Record unexpected per-molecule exceptions with the operation, exception type,
and message. Repeated similar failures do not prove a global runtime failure
and do not discard usable molecules from other records or sources. Missing
required runtimes, failed source processes, shared input corruption, and file
publication errors remain explicit failures at their own boundaries. Keep the
failure report and actual source counts; never label a failed source as a
successful generator.

RDKit BRICS products are explicitly sanitized before descriptors and Morgan
fingerprints are calculated on that same molecular object. Fingerprint errors
are audited rather than treated as property/similarity rejection. Fragment
ordering and BRICS reagent ordering use the explicit local seed, independent
of the process-global random generator.

Successful partial processing reports `completed_with_molecule_errors` in the
molecule-failure summary. SDF, CSV and final counts describe only released
records; input counts and all exclusions remain in the audit. Quarantined
ACT/DEC records change benchmark availability and sample size: report those
losses and apply the existing downstream validation gates. Processing success
does not itself validate or nominate the surviving molecules.

TargetDiff is launched as `python -m scripts.sample_for_pocket` in the registered
`targetdiff-gpu` environment and receives a center/radius-cropped standard-residue
pocket. Do not use the separately installed CPU-only `targetdiff` environment.

When DiffSBDD has no reference ligand, its center-cropped pocket is written with
chain `A` and sequential integer residue IDs. This preserves coordinates while
avoiding the upstream selector's inability to represent blank chains or PDB
insertion codes.

ResGen is launched with the explicit validated checkpoint and a per-run sampling
config. Its installed CLI still requires a reference ligand SDF even when the
Step 1 center is supplied, so the validator rejects a selected ResGen plan when
that artifact is absent.

REINVENT4 uses only the fixed reviewed `reinvent_pubchem.prior`. The adapter
checks its SHA-256 and records the official Zenodo checksum verification. The
upstream model emits an internal metadata-hash warning under this Python version;
this is recorded, not hidden, and does not replace the external asset checksum.

The Vina/Meeko/PDBQT extras present in some upstream repositories are not called
by this skill. Step 2 performs molecule generation only; docking remains the
pipeline docking step.

## Mandatory downstream geometry QC

Pocket2Mol, ResGen, DiffSBDD, TargetDiff, DiffLinker, Lingo3DMol, and MolCRAFT
native coordinates may contain high strain, incorrect bond perception, or
shape-complementary but chemically implausible structures. Their raw docking
score must never be treated as validation.

When docking is enabled, validation requires both:

- `pose_evaluation.posebusters=true`
- `pose_evaluation.posecheck=true`

The config validator hard-fails if either is disabled while a native 3D
generator is selected. A generation-only run may stop at Step 2, but its
candidates remain unvalidated and must not be nominated.

## Outputs and identity

| Output | Contract |
|---|---|
| `library.sdf` | Screening-ready 3D records only |
| `library_manifest.csv` | Stable source/compound/conformer identity, optional microstate and parent mapping, cohort and role, standardized/original SMILES, descriptors, coordinate and multi-generator provenance |
| `library_report.json` | Per-source counts/failures, runtimes, model checksums, filtering and coordinate diagnostics |
| `molecule_failure_report.json` | Persisted per-record error ledger, repeated-program-error threshold, and final processing/tool status |
| `conformer_report.json` | Coordinate preparation counts, native/rebuilt provenance, and per-record failure evidence |
| `<tool>_work/run-*/execution.json` | Exact argv, environment ID, exit code, elapsed time, and log paths for one isolated invocation |

IDs define downstream role:

- `ACT_*`: ChEMBL active candidate
- `GEN_*`: generated candidate from any generator
- `DEC_*`: validation-only decoy; never nominate as a candidate

Important SDF tags include `source`, `candidate_role`, `source_tools`,
`also_generated_by`, `coordinate_origin`, and `raw_pose_preserved`.
Preserve `source_record_id -> compound_id -> microstate_id|null -> conformer_id`
with source/cohort identity and coordinate lineage. A null `microstate_id`
means that no state expansion occurred; it is not missing preparation evidence.
Derived protomer/tautomer/stereo hashes, ensemble indices, population priors,
and state-nomination blockers are not mandatory tags for this route.

Pocket2Mol retains every parseable, nonempty candidate for downstream category
assessment. MW/LogP/rotatable-bond thresholds and the preliminary pharmacophore
check produce annotations only. A single hard-bound breach or blueprint mismatch
does not delete a candidate before property-filter/structure-evaluator finish
the relevant category fold.

## Failure semantics

| Failure | Behavior |
|---|---|
| Enabled runtime/env/checkpoint fails global preflight | Hard fail before any scientific step |
| One selected source fails at generation time | Record source failure; continue other selected sources |
| A generator exits zero but writes no valid molecules | Treat as zero-output failure, not success |
| Individual molecule parsing/declared-stereo/conformer preparation fails | Isolate the failed record, save identity and reason, continue usable records |
| Multiple molecules raise the same per-record exception | Record each failure; retain valid records without an arbitrary batch-abort threshold |
| Every selected candidate source yields zero | Step 2 hard fail |
| DUD-E and RDKit decoys both unavailable | Candidate-only library; benchmark metrics are unavailable |
| Native pose graph changes during standardization | Reject native pose and use ETKDG fallback |

DiffLinker commonly produces a mixture of valid and invalid bond assignments
after Open Babel conversion; the adapter filters each record individually.

## References

- `references/chembl_api.md`
- `references/dude.md`
- `references/pocket2mol.md`
- `references/external_generators.md`

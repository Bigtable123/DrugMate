---
name: structure-evaluator
description: Validate, filter, and rerank docked SBDD poses. Use after docking to enforce stable pose identity, run PoseBusters and PoseCheck quality gates, use PLIP for an interaction blueprint and record ProLIF evidence for audit, run opt-in RTMScore rescoring, prepare Amber candidate topologies, run short OpenMM MD with MDAnalysis stability analysis, run AmberTools MM-GBSA endpoint rescoring, or create a planning-only OpenFE/TIES RBFE campaign.
---

# Structure Evaluator

Multi-dimensional quality assessment of docked poses using geometric validation, energetic analysis, and interaction blueprint matching.

For a full pipeline, `../pipeline-entry/TOOL_SELECTION_POLICY.md` is the only
tool-selection guide and `../pipeline-entry/QUERY_TO_CONFIG.md` only assembles
that explicit decision. The validator rejects incomplete/incompatible
combinations; this skill and the runner execute enabled tools without adding a
default evaluation route.

## When to Use

- After docking (canonical Step 4; legacy disk directory `step3`) to filter and rank docked poses
- When you need to validate pose quality beyond docking scores
- When you want to check if poses satisfy the required interaction pattern
- Before final candidate selection / lead optimization

## Input Schema

```json
{
  "receptor_file": "path/to/docking_prep/prepared_receptor.pdb",
  "docked_sdf": "path/to/docked_poses.sdf",
  "interaction_blueprint_file": "path/to/interaction_blueprint.json",
  "pocket_features_file": "path/to/pocket_features.json",
  "reference_ligand_file": "path/to/reference_ligand.sdf",
  "docking_scores_file": "path/to/docking_scores.json",
  "docking_report_file": "path/to/docking_report.json",
  "library_manifest_file": "path/to/library_manifest.csv",
  "top_n": 10,
  "selection_strategy": "ranked_per_compound",
  "max_poses": 200,
  "batch_size": 50,
  "resume": true,
  "overall_timeout_seconds": 5400,
  "posecheck_num_conformers": 10,
  "posecheck_thresholds": {
    "strain_warn_kcal": 8.0,
    "strain_fail_kcal": 60.0,
    "clash_warn_count": 0,
    "clash_fail_count": 2
  },
  "calibrate_posecheck_against_reference": false,
  "posecheck_native_baseline_limits": {
    "strain_max_kcal": 100.0,
    "clash_max_count": 3
  },
  "screening_mode": "V2",
  "candidate_source": "sampled_library",
  "chain": "A"
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `receptor_file` | string | Yes | Final explicit-H receptor PDB from docking preparation; reuse the same receptor used for docking |
| `docked_sdf` | string | Yes | Docked poses SDF (from docking skill or any source) |
| `interaction_blueprint_file` | string | No | Blueprint JSON from pocket-analysis (enables blueprint matching) |
| `pocket_features_file` | string | No | Pocket features JSON from pocket-analysis (enables pocket occupancy) |
| `reference_ligand_file` | string | No | Co-crystal ligand; used for PoseCheck calibration only when the explicit calibration flag is true |
| `docking_scores_file` | string | No | Per-pose GNINA scores; binds scores to the same pose diagnostics |
| `docking_report_file` | string | Pipeline: Yes | Step 3 commit marker sealing the docked SDF and score file. The orchestrated pipeline requires it; an external/standalone QC run may omit it, but its docking evidence cannot later be treated as verified for nomination. |
| `library_manifest_file` | string | No | library_manifest.csv from molecule-generator (enables SMILES in report) |
| `top_n` | integer | No | Number of top candidates to include in summary (default: 10) |
| `docking_attempt` | integer | No | Default `1`; provenance for historical or explicitly requested docking runs. Included in the batch fingerprint; does not trigger a retry. |
| `selection_strategy` | string | No | `ranked_per_compound` (default): evaluate docking-ranked poses up to the cohort quota (20 when ≤50 compounds, 5 when 51–1000, 1 above 1000). Explicit legacy strategies remain available. |
| `poses_per_compound` | integer | No | Additional per-compound cap for `ranked_per_compound` (default `20`), or shortlist size for `best_per_compound`; it cannot raise the adaptive quota. |
| `max_poses` | integer | No | Cap actual evaluations and report partial if selected poses remain pending. All input poses stay in the source SDF and manifest. |
| `batch_size` | integer | No | Batch size for resumable evaluation (default: 50) |
| `profiler_workers` | integer | No | Independent profiler concurrency (default: 2; bounded by CPU affinity and 4). Set 1 to serialize. Batches remain sequential. |
| `checkpoint_interval_seconds` | number | No | Minimum interval between full cumulative checkpoints (default: 300; 0 writes every batch). Every batch always saves its full resumable report; first/last batch, timeout, and final completion force a cumulative report. |
| `resume` | boolean | No | Reuse completed `batches/batch_*/evaluation_report.json` files (default: true) |
| `overall_timeout_seconds` | integer | No | Stop between batches after this budget and write a resumable partial report (default: 5400) |
| `posecheck_thresholds` | object | No | Absolute soft/hard marking bounds by default: strain 8/60 kcal/mol and heavy-atom clashes 0/2 |
| `posecheck_num_conformers` | integer | No | Reference conformers for PoseCheck strain (default: `10`), including native-ligand calibration. Pipeline config: `tools.posecheck.params.num_conformers`. |
| `calibrate_posecheck_against_reference` | boolean | No | Explicitly opt into adding a validated native-ligand baseline to the configured bounds (default: false) |
| `posecheck_native_baseline_limits` | object | No | Sanity limits for accepting native-ligand PoseCheck calibration; if exceeded, fall back to absolute thresholds |
| `screening_mode` | string | No | `V1` marking / `V2` standard / `V3` strict level-2 category fold. Unset means it is derived from `candidate_source`, falling back to `V2` |
| `candidate_source` | string | No | Provenance used to derive `screening_mode`: `molecule_generation` → `V3`, `sampled_library` → `V2`, `user_supplied_full` → `V1` |
| `chain` | string | No | Protein chain for pocket occupancy (default: "A") |

Pose selection preserves original docking identity. Final
`evaluation_report.json` uses the original docking `pose_index`, retains the
per-record `sdf_molblock_sha256` coordinate witness, and records `pose_id`,
chemical/source identity, `source_pose_index`, `batch_pose_index`, coordinate
lineage, and selection ranks. When `docking_report_file` is supplied, it must
seal the exact `docked_sdf` and `docking_scores_file`; the report then records
their `source_docking_artifact_set` identity for downstream verification.

PoseCheck reuses the final `prepared_receptor` from docking preparation, normally
prepared by Reduce; preserve the selected backend's result when the protocol uses
PDB2PQR. The local wrapper reads its existing atoms and coordinates through RDKit,
without rerunning Hydride or converting the receptor through ProLIF. An unreadable,
empty, coordinate-free or hydrogen-free receptor produces failed PoseCheck evidence;
complete receptor preparation before evaluation. Do not substitute the structure
stage's `prepared_receptor_noH` or `provisional_receptor`. See
`references/posecheck.md` for loading and clash-count semantics.

## Workflow

For the explicitly enabled tools, execute in this order:

PoseBusters, PoseCheck, PLIP and ProLIF may execute concurrently under
`profiler_workers`; join their results before blueprint matching and verdicts.
Keep this budget within the run's CPU allocation and do not add nested batch or
pose worker pools. Input SDF blocks, blueprint/pocket JSON, receptor subpocket
coordinates and content signatures are cached only within one invocation, with
file-change invalidation. Reports use compact JSON and atomic replacement;
resumable batch reports retain every pose even between cumulative checkpoints.
Batch reports record `timings_seconds` for PoseBusters, PoseCheck, PLIP, ProLIF,
pose geometry extensions, and batch wall time. Final reports sum the recorded
batch timings, including reused batches' original measurements; overlapping
profiler times must not be added to infer wall time.

1. Bind every SDF record to a canonical chemical identity, pose ID, coordinate
   revision, and source lineage. Reject records that cannot be identified.
   Keep `identity_qc` enabled for every orchestrated Step 5 run; config
   validation treats this integrity check as mandatory.
2. By default, choose each compound's docking-ranked poses before QC using
   `engine_pose_rank` (existing score order when rank is unavailable). Evaluate
   up to 20 poses per compound when there are at most 50 compounds, up to 5 for
   51–1000 compounds, and only 1 above 1000 compounds. Evaluate every selected
   pose directly; do not advance to another pose based on a top-1 QC failure.
   Tool errors and missing evidence are unresolved, not scientific failures or
   passes. Preserve the full input SDF and manifest; unevaluated poses are never
   candidates. `max_poses` limits actual evaluations. An exhausted budget with
   selected poses remaining produces a resumable partial report and cannot
   trigger an all-poses-failed classification. Use `all` explicitly when the
   research question requires the full pose ensemble.
3. Run PoseBusters with per-check registered marking rules. A registered,
   chemically critical failed check is a `fail` mark. An unregistered check has
   no threshold and remains `not_applicable`; an aggregate failure that names no
   registered check is retained as audit evidence but cannot invent a separate
   hard gate. Table-declared soft failures (for example internal energy or
   water/metal-context checks) mark `flag`. Run
   PoseCheck with 10 reference conformers, absolute strain 8/60 kcal/mol and
   heavy-clash 0/2 bounds; native-baseline calibration is explicit opt-in and
   uses the same conformer count. This count controls reference sampling, not
   how many docking poses are retained; fewer conformers reduce sampling
   coverage, so record the setting and keep it consistent within the cohort.
4. Run PLIP and ProLIF independently. PLIP supplies the blueprint gate and
   ranking metrics; ProLIF matches are recorded as an independent audit cross-check.
   Never union their hits into a synthetic pass. Record detector disagreement
   without changing the PLIP-based gate.
   PLIP must use the exact candidate binding site `LIG:L:1` written by its
   per-pose adapter. Receptor cofactors and modified residues remain in the
   complex but their own ligand sites are excluded. A missing candidate site is
   a tool failure/NA; a detected candidate site with no contacts is a real zero.
5. Keep all pose diagnostics in the report, but export at most one pose per
   compound. Prefer a `pass`/`warn` pose; if every assessed pose is inconclusive,
   retain one with `post_qc_verdict=not_evaluated` for downstream review, where it
   cannot be nominated. A failed pose or a backup pose never evaluated by QC
   does not reach `selected_poses.sdf`.
   Only when every selected pose has been evaluated and explicitly failed,
   `selection.pose_exhausted_compounds` records `route=reject_molecule`,
   `reason_code=all_poses_failed`, the evaluated count and individual QC reasons.
   Skip that molecule and continue without redocking it or the cohort. This
   describes pose failure under the current protocol, not biological inactivity.
   Tool faults, unreadable structures and unmeasured poses do not qualify.
6. Preserve all evaluated pose rows for coordinator audit. The coordinator marks
   the raw Step 5 metrics, folds them per category, then applies descending
   RANK.AVG metric/category percentiles,
   signed category weights, and parent weights. The exported selection manifest
   fixes the representative used by property filtering and coordination;
   coordination must not replace it after property evidence has been computed.
   Do not average raw scores with unlike methods, units, directions,
   stages, or applicability.
7. If `rtmscore` is explicitly enabled, score every pass/warn pose in the
   evaluated cohort with the pinned RTMScore ensemble. `rescore_rtmscore` may
   make that one signal the primary within-cohort ranking metric; only after
   this join does the evaluator select one pose per compound. Never label the
   score as an experimental affinity or mix it arithmetically with docking
   energy. Read `references/rtmscore.md` first.
8. If `candidate_topology` is explicitly enabled, prepare candidate-specific
   Amber inputs only from the post-QC survivor SDF.
9. If `openmm_md` is explicitly enabled, request neutralized ff14SB/GAFF2,
   TIP3P explicit-solvent periodic topologies and run only the bound original
   top-K cohort. Re-hash each solvated prmtop/inpcrd before launch. Use the
   CUDA-preflighted protocol and at most a 2 fs timestep (no HMR). This is a short
   stability diagnostic, not an equilibrium or free-energy calculation.
10. If `mdanalysis` is explicitly enabled with `openmm_md`, analyze the same
    environment's trajectories relative to the saved equilibrated
    production-start coordinates after protein alignment and no-jump PBC
    reconstruction. Minimize ligand RMSD over exact element/Amber-type topology
    automorphisms without a ligand-only refit. Report ligand RMSD, center
    displacement, protein-backbone RMSD, and starting-contact retention. The
    shared numeric bounds are one-sided screening heuristics, not a
    claim of equilibrium stability or a universal experimental cutoff.
11. If `mmgbsa` is explicitly enabled, consume only the built
   `candidate_topologies.json` records. Keep the current single-frame,
   no-entropy result labeled as a low-reliability endpoint reranking signal;
   never present it as an absolute binding free energy or average it with a
   docking score. Read `references/mmgbsa.md` for the enabled tool's contract.
12. If an RBFE backend is explicitly enabled, consume only those post-QC
   survivors. Run OpenFE/TIES preflight or planning after QC, keep the backend
   independent of `candidate_topology`, and never interpret a passed plan as a
   computed ΔΔG.

When GNINA is the explicitly selected Step 4 engine, reuse its pose-bound CNN
fields as the cheap score gate; Step 5 does not rerun GNINA. Uni-Dock/Vina runs
retain and consume their own native metrics through the same `pose_scores`
contract.

## Verdict Logic

The pose verdict is the second level of the gate. Every registered scalar first
receives a level-1 `pass`/`flag`/`fail` mark, those marks are folded per parent
metric category under the active `screening_mode` (`V1` marking / `V2` standard,
the default / `V3` strict), and only that fold — published as `metric_gate` with
per-category detail in `metric_category_folds` — turns marks into a pose verdict:

| Verdict | Criteria |
|---------|----------|
| **pass** | Required evidence is available and no metric was marked `flag` or `fail`: geometry/strain/heavy-atom clashes pass and every required Blueprint group satisfies its `all_of`, `any_of`, or `at_least_k` logic |
| **warn** | At least one `flag` or `fail` mark exists but no category fold rejected, or the fold returned `modify`; a modify-only PoseBusters check fails; or a required constraint is unresolved |
| **fail** | A metric category tripped the active screening mode (`count(fail) >= 3` or `ratio(fail) >= 50%` under `V2`; `>= 2` or `>= 30%` under `V3`), or the pose has an unreadable identity |
| **not_evaluated** | No usable evaluation evidence is available; preserve `Score=NA` for audit without counting it as a pass or scientific failure. Unevaluated backup poses stay in the input SDF/manifest and never enter candidate nomination. |

A single hard-bound breach is therefore a `fail` *mark*, recorded under
`verdict_reasons.fail`, and it caps the pose at `warn` unless its category fold
rejects. The count and ratio branches are independent alternatives and neither
carries a minimum fail count, so a lone breach in a
category with few other labels is ratio 100% and does reject. No additional
observed-coverage guard modifies the numeric count OR ratio formula.
An unattributed PoseBusters aggregate failure is retained under
`verdict_reasons.not_applicable`; without a named registered check it cannot
participate in the category fold or reject a pose. `screening_mode` is a step parameter; left unset
it is derived from `candidate_source` — a generative round defaults to `V3`,
`sampled_library` to `V2`, `user_supplied_full` to `V1`, and an unknown provenance
to `V2`. Under `V1` no category ever rejects, so only an unreadable identity can
still produce `fail`.

Evidence counts as usable when it comes from any of the five interaction/geometry
profilers or from pocket occupancy, RTMScore, or MM-GBSA. A PoseBusters check
with no registered rule carries no threshold, so it is recorded as
`posebusters_unregistered:<check>` under `verdict_reasons.not_applicable` and is
never a fail — several checks report False purely because production screening
has no reference ligand. A PLIP row that
reported no interaction count is NA rather than 0, which would have hit
`hard_lower=2`. A
Blueprint verdict of `unresolved` — what an absent or crashed profiler produces —
scores as `None` and leaves both the numerator and the denominator,
rather than landing between the hard and soft bounds of a weight-5 metric.

PLIP `num_interactions` is a weight-zero plausibility metric: `0` or `1` marks
`fail` and `>=2` marks `pass`; no synthetic sub-1.36 Å counter is used. That mark
reaches a rejection only through its category fold.
`required_hit_rate <0.80` is a weight-zero soft MODIFY diagnostic and never
replaces exact required-group logic. Atom, chemical-role, geometry, and
mediator fields are matched when present. A field or interaction class outside
the current PLIP wrapper's capabilities is reported as
`unmatchable`/`unresolved`, not as a missing pose interaction. The supplied
metric profile gates `desired_ranking_score` at
0.10/0.50: values `<0.10` mark `fail`, `[0.10,0.50)` mark `flag`, and
`>=0.50` mark `pass`. `desired_evaluable_coverage` has `soft_lower=0.80` and
no hard lower; values below `0.80` mark `flag`, while `<0.30` additionally
stops trusted ranking and requests re-evaluation. Coverage never becomes a
molecule or pose rejection by itself.

The default absolute thresholds (from
`shared/constants.py :: DEFAULT_POSECHECK_THRESHOLDS`) are:

```json
{
  "strain_warn_kcal": 8.0,
  "strain_fail_kcal": 60.0,
  "clash_warn_count": 0,
  "clash_fail_count": 2
}
```

Set `calibrate_posecheck_against_reference=true` only when the run intentionally
uses native-relative thresholds. The report records the effective bounds and
calibration mode. See `../METRIC_POLICY.md` for all executable rules, weights,
unavailable table rows, and interpretation choices.

## Output Files

| File | Description |
|------|-------------|
| `evaluation_report.json` | Full per-molecule diagnostics, preserved PoseBusters checks/success ratio, Stage 07 pose/docking metric bounds and gates, protocol ranking blocks, and verdict; Stage 08 property/ADMET placeholders are excluded |
| `evaluation_summary.json` | Aggregate stats: pass/warn/fail counts |
| `evaluation_input_poses.sdf` | All ordered docking poses, including unevaluated backups, in the default strategy; pre-QC shortlist for legacy strategies |
| `evaluation_input_manifest.csv` | Identity and docking-score manifest for the evaluation input and backup pool |
| `selected_poses.sdf` | Final post-QC representative poses for downstream property filtering and coordination |
| `selection_manifest.csv` | Post-QC identity, source pose, verdict, and primary-score provenance |
| `batches/batch_*/evaluation_report.json` | Resumable per-batch reports |
| `candidate_topologies/candidate_topologies.json` | Optional candidate-specific topology manifest, produced only when explicitly enabled |
| `rtmscore/rtmscore_manifest.json` | Optional exact-pose RTMScore ensemble results and model/runtime provenance |
| `short_md/short_md_manifest.json` | Optional explicit-solvent OpenMM simulation manifest with separate simulation/analysis states |
| `short_md/trajectory_stability_report.json` | Optional MDAnalysis descriptive stability report |
| `mmgbsa/mmgbsa_manifest.json` | Optional pose-bound AmberTools MMPBSA.py manifest; current score is single-frame and excludes entropy |
| `openfe_rbfe/openfe_rbfe_manifest.json` | Optional OpenFE applicability/mapping-network manifest; contains no simulation result |
| `ties_rbfe/ties_rbfe_plan.json` | Optional TIES runtime/pair-stage plan; contains no simulation result |

Resume a batch only when its content fingerprint matches the receptor, pose
coordinates and IDs, relevant artifacts, full enabled-tool matrix, parameters,
runtime/model/checkpoint/version metadata, and implementation files. Recompute
instead of trusting an old report when any component changes.

### Per-Molecule Report Structure

The following JSON-with-comments illustration omits individual interaction
records; actual reports contain the complete records and use plain JSON.

```jsonc
{
  "pose_index": 0,
  "mol_name": "ACT_0001",
  "posebusters": {
    "all_passed": true,
    "violations": [],
    "num_violations": 0
  },
  "posecheck": {
    "steric_clash_count": 0,
    "steric_clash_count_all_atom": 2,
    "steric_clash_count_mode": "heavy_atom_only",
    "strain_energy_kcal": 2.3
  },
  "plip": {
    "num_interactions": 5,
    "interactions": [/* interaction records omitted */]
  },
  "blueprint_match": {
    "hard_gate_applicable": true,
    "hard_gate_status": "passed",
    "required_hits": 2,
    "required_total": 3,
    "required_group_results": [
      {"group_id": "hinge", "logic": "any_of", "status": "passed"},
      {"group_id": "catalytic_anchor", "logic": "all_of", "status": "passed"}
    ],
    "desired_hits": 2,
    "desired_weighted_score": 0.67,
    "verdict": "pass"
  },
  "pocket_contact_fraction": {
    "hydrophobic_region": 0.72,
    "polar_region": 0.45,
    "charged_region": 0.10,
    "aromatic_region": 0.83,
    "overall": 0.55
  },
  "diagnosis": [
    "GOOD: Strain energy 2.3 kcal/mol within threshold",
    "GOOD: all required Blueprint groups satisfied",
    "LOW_OCCUPANCY: charged_region (10%)"
  ],
  "recommendation": "modify",
  "modification_hint": "Extend scaffold into charged region near ASP145, GLU81",
  "verdict": "pass"
}
```

## Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| RDKit | SDF parsing, mol conversion | `pip install rdkit` |
| PoseBusters | Geometric validation | `pip install posebusters` |
| PoseCheck | Strain + clashes | `pip install posecheck` |
| PLIP | Interaction analysis | `pip install plip` |
| ProLIF + MDAnalysis | Independent interaction fingerprint | `pip install prolif MDAnalysis` |
| AmberTools | Optional candidate topologies | Use the registered `amber` Conda environment |
| RTMScore | Optional ML pose rescoring | Use the pinned `${SBDD_ENV_RTMSCORE_PREFIX}` prefix; CPU fallback is supported |
| OpenMM + MDAnalysis | Optional short MD and trajectory analysis | Both use the registered `openmm` Conda environment; its CUDA Context is functionally probed before execution |
| AmberTools MMPBSA.py | Optional single-frame MM-GBSA reranking | Use the registered `mmgbsa` environment (mapped to `gmxMMPBSA`) |
| OpenFE | Optional mapping-level RBFE planning | Use the registered `openfe` Conda environment |
| TIES + TIES_MD + OpenMM + AmberTools | Optional TIES pair/campaign planning | Builder API/CLI, runner/analysis API/CLI, an ABI-compatible and functionally tested OpenMM platform, antechamber, parmchk2, tleap, and sqm must all be present in the registered `ties` Conda environment; do not inject the Python 3.11 `drugmate` OpenMM build into this Python 3.10 runtime |

RDKit identity is part of the default integrity contract. Other tools report
explicit unavailable/unresolved evidence rather than manufacturing a pass.
See `references/candidate-topology.md` for the Amber topology execution
contract. See `references/rbfe-planning.md` for OpenFE or TIES execution and
applicability contracts.
See `references/mmgbsa.md` for the endpoint-reranking contract,
`references/rtmscore.md` for ML rescoring, and `references/short-md.md` for
OpenMM/MDAnalysis execution semantics.

## Usage

```bash
# Full evaluation with blueprint + pocket occupancy
python skills/bio_research/structure-evaluator/run.py \
  --input '{"receptor_file": "step3/docking_prep/prepared_receptor.pdb", "docked_sdf": "step3/docked_poses.sdf", "docking_scores_file": "step3/docking_scores.json", "docking_report_file": "step3/docking_report.json", "interaction_blueprint_file": "step1/interaction_blueprint.json", "pocket_features_file": "step1/pocket_features.json", "selection_strategy": "ranked_per_compound", "posecheck_num_conformers": 10, "batch_size": 50, "overall_timeout_seconds": 5400}' \
  --output-dir ./step4/

# Without blueprint (geometric + energetic checks only)
python skills/bio_research/structure-evaluator/run.py \
  --input '{"receptor_file": "docking_prep/prepared_receptor.pdb", "docked_sdf": "docked.sdf"}' \
  --output-dir ./eval_out/
```

## Error Handling

| Scenario | Behavior |
|----------|----------|
| No tools installed (only RDKit) | `success: true`; every pose is `not_evaluated`, carried with `Score=NA` and ranked behind every measured pose, and none is nominated |
| RDKit unavailable | `success: false` once, with an `operational_errors` entry `{"tool": "rdkit", "code": "rdkit_unavailable"}`. `identity_qc` cannot be disabled, so this is a run fault, never a per-pose verdict |
| PoseBusters unavailable | Skip, record in `steps_skipped` |
| PoseCheck unavailable | Skip, record in `steps_skipped` |
| PoseCheck receptor is unreadable, empty, has no coordinates or no explicit H | Record failed PoseCheck evidence/NA with a preparation error; supply the final prepared receptor before retrying |
| PLIP unavailable | Keep the Blueprint result as `unresolved`; do not call unsupported constraints missing |
| PLIP/ProLIF disagree | Preserve both profiler results; the Blueprint gate follows PLIP and records ProLIF disagreement for audit |
| Blueprint file missing | Skip matching, still run geometric checks |
| Individual pose fails parsing or PoseCheck execution | Preserve `success=false` and the per-pose error, mark that evidence `tool_failed`/NA, and continue others |
| RDKit ETKDG generates no reference conformer | Record `not_evaluable` with the reason; strain/TFD/RMSD are NA. Do not optimize an empty conformer set or fabricate zero-valued metrics. |
| Pocket occupancy has no ligand/region coordinates or an empty denominator | Emit `null` for the unavailable region/overall scalar; never convert missing evidence to zero |
| Resume fingerprint mismatch | Recompute that batch and record the rejected resume evidence |
| Candidate topology is inapplicable | Write `not_evaluable` with metal/covalent/cofactor evidence; do not run a generic fallback |
| RTMScore output omits/duplicates a pose ID | Reject the entire score join atomically; keep the evaluation cohort unmodified |
| A top-K topology is not evaluable | Preserve its original rank as `not_evaluable`; never backfill from rank K+1 |
| A passed topology hash is missing/mismatched | Do not launch OpenMM; mark that bound rank as an input-integrity failure |
| OpenMM simulation fails | Attribute the operational failure to `openmm_md`; dependent trajectory analysis is `not_evaluable` |
| MDAnalysis fails after successful MD | Keep simulation status `passed`, mark only `mdanalysis` failed, and preserve the trajectory |
| Docked reference has no contacts at the configured cutoff | Keep RMSD/displacement analysis; report only contact retention as `not_evaluable` |
| MM-GBSA is scientifically inapplicable | Write a successful `not_evaluable` manifest and do not fabricate a score |
| Explicit MM-GBSA runtime/output fails | Write the partial manifest and fail Step 5 operationally |

## Upstream / Downstream

### Upstream
| From | File | Usage |
|------|------|-------|
| docking preparation | `docking_prep/prepared_receptor.pdb` | Final receptor, including its prepared hydrogens, reused for clash checks |
| docking | `docked_poses.sdf` | Poses to evaluate |
| pocket-analysis | `interaction_blueprint.json` | Expected interaction pattern |

### Downstream
| To | File | Usage |
|----|------|-------|
| Ranking/Selection | `evaluation_report.json` | Filter candidates by verdict |
| Lead optimization | Top "pass" molecules | Starting points for optimization |
| MM-GBSA reranking | `mmgbsa_manifest.json` + pose-bound `pose_scores` | Explicit top-K endpoint-energy signal |
| Short MD / stability analysis | `short_md_manifest.json` + `trajectory_stability_report.json` | Explicit opt-in simulation and descriptive trajectory evidence |

## References

- `references/posebusters.md` — PoseBusters API, config options, check descriptions
- `references/posecheck.md` — PoseCheck API, strain energy calculation details
- `references/candidate-topology.md` — opt-in Amber topology contract and applicability boundaries
- `references/rtmscore.md` — pinned ensemble, exact pose join, and ranking semantics
- `references/short-md.md` — explicit-solvent OpenMM and MDAnalysis protocol contract
- `references/mmgbsa.md` — explicit AmberTools MMPBSA.py execution and single-frame scientific semantics
- `references/rbfe-planning.md` — OpenFE/TIES planning-only execution and applicability boundaries

## Optional metric extensions

Candidate-pose `rdkit_pose`, `freesasa`, `arpeggio`, `luna`, and `oddt` tools are
registered and off until explicitly selected. `prolif` adds reference fingerprint
Tanimoto for audit when a usable reference is supplied; it never gates or scores. `rtmscore` now contributes the
table's DYN weight +3; analyzed OpenMM trajectories supply three DYN metrics.
See [metric extension controls and definitions](../pipeline-entry/references/metric-extensions.md).
Do not infer RBFE free energies from successful OpenFE/TIES planning.

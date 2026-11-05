---
name: drug-property-filter
description: Filter SBDD candidates with property and risk models.
---

# Drug Property Filter

This skill has two registry uses:

- `property_prefilter` (registry order 3, disk `step2_5`): RDKit-only library
  screening before docking.
- `property_filter` (registry order 6, disk `step5`): RDKit → ADMET-AI +
  ChemFH → top-K AiZynthFinder after pose evaluation.

The main process runs in the registered `drugmate` scientific runtime selected
through deployment configuration. ADMET-AI, ChemFH, and AiZynthFinder are
isolated JSON workers in Conda `admet_ai`, prefix
`${SBDD_ENV_CHEMFH_PREFIX}`, and Conda `aizynth`; never import those backends into
Drugmate.

Whether either step and each of its tools is enabled is decided only through
`../pipeline-entry/TOOL_SELECTION_POLICY.md` and recorded explicitly through
`../pipeline-entry/QUERY_TO_CONFIG.md`. The validator rejects an incompatible
matrix; the runner and this skill do not choose between prefilter and final
filter behavior.

## Execution contracts

### Pre-docking prefilter

When `property_prefilter` and its RDKit tools are explicitly enabled, the
adapter performs this contract:

1. Calculate MW, LogP, TPSA, HBD/HBA, rotatable bonds, QED, Lipinski
   violations, and Ertl SA score with RDKit.
2. Collect every PAINS A/B/C, Brenk, and NIH alert, including matched atom
   indices.
3. Write `prefiltered_library.sdf`; keep `pass`, `flag`, and `fail` by default.
   A `fail` here is only the level-1 mark of one metric, so carrying it forward
   is deliberate: the keep/modify/reject decision belongs to the per-category
   fold.
4. Do not run ADMET-AI, ChemFH, pose deduplication, or retrosynthesis.

At library intake, preserve declared `compound_id` values (including existing
`CompoundID` / `compound_name` aliases). A unique title may supply a missing ID;
blank or repeated titles receive deterministic source-record IDs. The result
JSON, manifest CSV and filtered SDF carry the same `compound_id`, so removing or
reordering records does not renumber survivors. Original SDF titles and
parent/microstate/conformer metadata remain available. Docking preparation sets
the engine-input title to `compound_id` and retains its original value as
`input_mol_name`, because native engines also use titles to join scores. It uses
the same ID rule when prefiltering is disabled; duplicate explicit IDs require
unique species IDs from the caller and are never silently merged.
Verified `ACT_` / `DEC_` labels are read from `input_mol_name` after engine-title
normalization, so an explicit ID does not change a compound's benchmark role.

### Final property filter

When `property_filter` and the corresponding tools are explicitly enabled, the
adapter performs this contract:

1. Read Step 4's selected-pose SDF and preserve `chemical_identity_id`,
   `chemical_id`, `source_identity`, `pose_id`, parent/revision metadata,
   `source_pose_index`, `selection_rank`, and lineage JSON.
   If the existing post-QC SDF is empty, the pipeline passes
   `allow_empty_input=true`: write an empty successful report with
   `steps_skipped.reason_code=no_selected_poses` and invoke no property tools.
   Do not substitute rejected docking poses. Standalone empty inputs remain
   errors unless this behavior is explicitly enabled; missing or malformed
   files remain input errors.
2. When `dedup` is enabled, deduplicate repeated poses by `chemical_identity_id`; fall back to
   `chemical_id`, InChIKey, then canonical isomeric SMILES. The representative
   is the lowest `selection_rank`, then the earliest input record. A single
   identity mapped to conflicting graphs is a hard data-integrity failure.
3. Run the explicitly enabled RDKit first-pass tools. PAINS, Brenk, and NIH are
   all counts under the same relaxed band: `0` marks `pass`, `1-4` marks `flag`,
   and only `>4` marks `fail`. PAINS_A/B/C are summed into one
   family count first, because the catalog defines a single `pains_alert_count`
   metric. A single alert hit therefore never ends a molecule here. Covalent
   status adds review context but never rescues an alert by itself; only an
   explicit per-catalog hard-fail opt-out may drop a catalog to flags.
4. When `admet` is enabled, send all structurally valid unique candidates to ADMET-AI in one CPU batch.
   Join results only by the generated `candidate_id`; missing, extra, duplicate,
   or SMILES-mismatched outputs are rejected. A `complete` result must contain
   exactly the scorer's 41 registered native endpoints with finite values:
   probabilities must be in `[0,1]`, PPBR in `[0,100]`, and VDss, half-life,
   hepatocyte clearance, and microsome clearance must be nonnegative. Endpoint
   rename/addition/removal or any domain violation fails closed. ADMET-AI being
   switched off is not a finding about any molecule: it is recorded in
   `steps_skipped`, the categories it owns leave the planned coverage
   denominator, and no candidate is downgraded. An explicitly
   optional ADMET outage retains missing coverage without an independent
   nomination veto; required endpoints and conflicting required/N/A declarations
   still require evidence or plan repair.
5. When `chemfh` is enabled, send the same valid unique candidate batch to ChemFH on CPU. Require exact
   `candidate_id` + SMILES coverage and validate both complete and partial
   payloads endpoint-locally. Retain all seven model endpoints and 17 native
   rule results. Use the worker's native rule screening without SVG rendering;
   the pipeline only consumes SMARTS and atom-index evidence from these rules.
   Apply each valid probability's table gate at `0.30/0.90`;
   `>0.90` marks `fail`, while the provider's native positive call (`>0.5`)
   remains a conservative `review_only` floor — that is a prediction about the
   molecule, so it counts. Reactive compounds weigh `1`, promiscuous compounds
   and colloidal aggregators `0.5` each, and the four fluorescence/interference
   endpoints `0`. Map the provider's last ten rules to the tabulated
   `substructure_alert_count` (count hit rules, not SMARTS multiplicity) and use
   the numeric `soft=2, hard=10, weight=2` contract. Individual native rule hits
   remain evidence rows, and ChemFH's PAINS result does not add a second PAINS
   gate. A missing optional ChemFH run is an audit gap; a required run blocks
   nomination until rerun or a versioned, full-coverage, confidence-bearing
   equivalent fallback is recorded. An ordinary RDKit FilterCatalog pass is not
   silently equivalent to the seven ChemFH mechanisms.
6. Mark every available RDKit/alert/ADMET/ChemFH raw metric, then fold the marks
   per parent category under the active `screening_mode`. A candidate the fold
   rejects must not be submitted to retrosynthesis; a candidate carrying `fail`
   marks that no category fold rejected still may be.
7. When `aizynthfinder` is enabled, sort candidates the fold did not reject and
   that have complete ADMET by Step 4
   `selection_rank`, then stable input order. Send only the first `top_k`
   (default 25) to AiZynthFinder. Preserve `solved`,
   `not_found_within_budget`, and `error` as different states; a partial route
   is not a solved route.
8. Project the existing Step 5 decisions onto the shared three-tier contract:
   `rejected → reject`; `fail`, `flag`, and `review_only → modify`, otherwise
   `keep`. A bare `fail` verdict is capped at `modify` because it is a level-1
   mark; only the fold in `metric_gate` and a structural failure (unreadable
   chemistry, explicitly `rejected` eligibility) may reject. Emit both `gate` and
   `classification`; reject always dominates.
   A non-pass verdict that rests only on a coverage gap is exempt from that
   projection. `compute_verdict` reports the split in `rdkit_evidence`
   (`molecular_findings` versus `coverage_gaps`), and a missing sascorer, an
   unavailable FilterCatalog, a fully disabled RDKit pass or a raised RDKit pass
   goes to `protocol_review_reasons` only: it blocks automatic nomination without
   downgrading the chemistry, because the executable policy forbids penalising a molecule for
   a tool or coverage problem. Rows carrying no `rdkit_evidence` keep the older,
   more conservative reading.

AiZynthFinder uses the registered local configuration at
`${SBDD_AIZYNTHFINDER_CONFIG}`, default USPTO expansion/filter
policies, and the ZINC stock.

## Metric gate and weighted ranking

The gate has two levels. For each structurally valid record,
materialize all legacy/backend marks and then mark raw scalar values before
percentile scoring: outside a hard bound is `fail`, outside the soft interval but
within the hard bounds is `flag`, and inside the soft interval is `pass`. Missing
bound sides are infinite; a missing value is `not_applicable` and is not counted.

A mark is never a molecule decision. Once a parent category is complete, count
its marks and apply the active `screening_mode`:

| mode | rejects when |
| --- | --- |
| `V1` (marking) | never |
| `V2` (standard, default) | `count(fail) >= 3` or `ratio(fail) >= 50%` in one category |
| `V3` (strict) | `count(fail) >= 2` or `ratio(fail) >= 30%` in one category |

The count and ratio branches are independent alternatives and neither carries a
minimum fail count: one `fail` among few labeled
metrics is ratio 100% and rejects. Missing values are excluded from the ratio;
no extra observed-coverage guard modifies the numeric formula.

A tripped category → `reject`; otherwise every evaluated category must have
strictly more than half its labeled metrics marked `pass` for `keep`. Other
labeled cases → `modify`; nothing marked → `not_evaluated`. The current fold
also controls molecular nomination eligibility: an old single-metric RDKit or
ChemFH review mark cannot veto a `keep`, while protocol evidence gaps still
require review. `screening_mode` is a step parameter; left unset,
`candidate_source` selects it — `molecule_generation` →
`V3`, `sampled_library` → `V2`, `user_supplied_full` → `V1`. Records the fold
rejects are removed before every ranking pool and have null category/step/overall
score and rank.

Every finite metric among the surviving records uses descending RANK.AVG
percentile inside one scientific cohort. `DEC_*` benchmark decoys use a
separate pool and never change nomination-candidate percentiles; ties share
average rank and one valid value receives `0.5`. A declared `cohort_id` is a
lineage label, not a substitute for the protein/pocket/round/profile/catalog
components; those components remain in the grouping key to prevent accidental
cross-cohort deduplication or normalization.
ExactMolWt, cLogP, TPSA, formal charge, and heavy-atom count first use the
hard/soft trapezoid desirability. A weight-zero metric still marks and reports
its percentile but contributes zero. Missing evidence leaves only that
candidate's category numerator and denominator; it does not disable the metric
   for complete candidates. A required planned missing weight-zero gate/status
   blocks automatic nomination through coverage/protocol status even though it
   is absent from weighted coverage; an explicitly optional ChemFH outage remains
   an audit gap. Signed weights produce
0–100 PHY/ALERT/SYN/ABS/DIST/MET/EXCR/TOX category scores,
category scores are percentile-ranked again, and parent weights produce Step 4
and overall fields. Planned, catalog, and category coverage are distinct;
`>=0.50` releases formal rank, `[0.10,0.50)` provisional rank, and `<0.10`
leaves score/rank null and asks for a backtrack. Incomplete planned coverage is
recorded in `coverage_gate`/`rank_status` and never rewrites the chemical
classification or rejects a molecule/pose.
`property_weighted_*` fields are deprecated aliases of the hierarchical
overall fields during the schema transition. Invalid structures remain in the
audit result with schema-v2 category/step/overall score and rank fields present
but null; they never enter a percentile pool.

The stable executable subset is:

- RDKit ExactMolWt, cLogP, TPSA, HBD/HBA, rotatable bonds, ring counts,
  FractionCSP3, formal charge, heavy-atom count, and chiral-center count.
  Chiral centers mark `<=2` `pass` and `>2` `flag`, with no hard bound. cLogP and
  TPSA have weight `3`; the other descriptor weights are `0`.
- Lipinski RO5 violations (weight `0`) and QED (weight `3`; default
  `hard_lower=0.10`, `soft_lower=0.45`). RO5 follows the relaxed numeric cells
  `soft_upper=1`, `hard_upper=3`: `0-1` marks `pass`, `2-3` `flag`, `>3` `fail`.
  `lipinski_max_violations` (default `3`) configures only the hard bound.
- PAINS A/B/C, Brenk, and NIH alert counts (weight `1` each inside ALERT) while
  retaining alert names and matched atoms. All three share one relaxed band:
  `0` marks `pass`, `1-4` `flag`, `>4` `fail`, with PAINS_A/B/C summed into one
  family count. Only explicit catalog rescue drops a catalog to flags. RDKit does
  not expose the table's normalized Brenk severity, so a severity-specific
  single-hit hard rule is not inferred.
- ChemFH's seven native model probabilities use the exact ALERT bounds from the
  table: `<=0.30` marks `pass`, `0.30-0.90` `flag`, `>0.90` `fail`. Reactive
  compounds weigh `1`, promiscuous compounds and colloidal aggregators `0.5`
  each, and the four fluorescence/interference endpoints `0`. Uncertainty and the
  individual native structural-rule matches remain audit context; they feed only
  the aggregated `chemfh.substructure_alert_count` (`soft=2`, `hard=10`, weight
  `2`) and create no second PAINS gate.
- SA score (`<=6` `pass`, `6-9` `flag`, `>9` `fail`, weight `2`).
- ADMET-AI: 41 exact native endpoints from the provider schema plus derived
  CYP substrate/inhibition counts. These implement ABS, DIST, MET, and
  TOX gates and signed weights. Lipophilicity and VDss use target-range
  transforms; `admet.params.bbb_profile` defaults to `non_cns`, while explicit
  `cns` reverses both the BBB gate and score direction. The three EXCR endpoints have blank
  bounds and weight zero, remain audit-only, and do not activate EXCR scoring
  alone. Exact endpoint-level bounds and weights are published in
  `../pipeline-entry/TOOL_CATALOG.md`. CYP2C9/2D6/3A4 substrate probabilities
  use signed weights `-0.1/-0.1/-0.5`, penalizing higher
  probability under the implemented scoring formula. Provider
  columns outside this exact contract are retained under
  `unregistered_predictions` as audit evidence and never satisfy or expand the
  complete scoring payload. Planned weighted coverage below `0.60` for ABS or
  below `0.70` for TOX requests evidence repair; it never becomes a molecular
  reject by itself.
- AiZynthFinder route-found, state score, reaction count, precursor count, and
  an explicit all-precursors-in-stock scalar computed by equality of precursor
  and in-stock counts. State score is the only weighted route metric
  (`2`), is accepted only with exact `state score` analysis provenance, and is
  comparable only under one fixed version/config/search budget. Search failure
  is a `flag` mark; reaction count `<=5` is `pass` and `>5` is `flag`
  (`soft_upper=5`, with no hard-fail bound), so this metric alone never adds a
  SYN `fail`. Price and route-cost
  fields remain NA because the current worker has no such data source.

When AiZynthFinder is enabled, candidates outside its configured `top_k` keep
their gate label and are marked `synthesis.status="not_selected_top_k"`: a
compute budget is not a finding about a molecule. Preserve the distinction between
[scientific gates and rank publication](../METRIC_POLICY.md#gate-and-score-contract).
The [delivery-label mapping](../shared/output_layout.py) uses RESERVE with
`gate_status=PASS` for a gate-passing candidate whose rank is not formally released.
Retrospective `DEC_*` decoys never enter that `top_k` cohort:
they are marked `synthesis.status="not_eligible"` with reason `benchmark_decoy`,
keep their existing nomination eligibility, and are counted in
`excluded_benchmark_decoys`. Decoys exist only to score the pipeline, and the
retrospective benchmark ranks them on docking scores alone, so spending
retrosynthesis budget on them would both waste search time and push real
candidates out of the searched slice. State score contributes for
candidates where it is finite; missing values are removed per candidate and
reduce planned coverage.

ADMETlab 3.0, SwissADME, RAScore, SCScore, SYBA, ASKCOS, and pkCSM do
not exist in this pipeline and their rows must be skipped. The registered
`admet` tool is ADMET-AI, an additional backend rather than an alias. Never map
ADMETlab/SwissADME/pkCSM thresholds onto ADMET-AI endpoints by name. Only the
43 exact ADMET-AI rules above enter the scorer; its other provider outputs
remain evidence only.

ChemFH is registered as provider-specific frequent-hitter/interference
evidence, not as an ADMET-AI alias. Its seven table-matched probabilities are
ALERT marks at `0.30/0.90` with reactive compounds at weight `1`, promiscuous
compounds and colloidal aggregators at `0.5`, and the four
fluorescence/interference endpoints at `0`; uncertainties and native structural
rules remain audit context and feed only the aggregated
`chemfh.substructure_alert_count`. See `references/chemfh.md` for endpoint,
asset, and overlap policy.

## Interface-v2 input

After policy selection and validation, the pipeline adapter supplies paths,
thresholds, structural facts, and the fully materialized explicit `tools`
block:

```json
{
  "molecules_sdf": "step4/selected_poses.sdf",
  "qed_threshold": 0.45,
  "sa_threshold": 6.0,
  "lipinski_max_violations": 3,
  "is_covalent": false,
  "screening_mode": "V2",
  "candidate_source": "sampled_library",
  "tools": {
    "lipinski": {"enabled": true, "params": {"lipinski_max_violations": 3}},
    "qed": {"enabled": true, "params": {"qed_threshold": 0.45}},
    "sa_score": {"enabled": true, "params": {"sa_threshold": 6.0}},
    "pains": {
      "enabled": true,
      "params": {
        "catalogs": ["PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"],
        "hard_fail_catalogs": [
          "PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"
        ]
      }
    },
    "admet": {"enabled": true, "params": {"timeout_seconds": 7200}},
    "chemfh": {
      "enabled": true,
      "params": {"timeout_seconds": 7200, "seed": 42}
    },
    "dedup": {"enabled": true, "params": {}},
    "aizynthfinder": {
      "enabled": true,
      "params": {
        "top_k": 25,
        "search_time_limit_seconds": 60,
        "iteration_limit": 100,
        "max_routes": 5,
        "return_first": true
      }
    }
  }
}
```

Standalone calls may instead supply `molecules_smiles`. Runtime/environment
paths are infrastructure metadata from the registry and must not be added to
user configs.

## Outputs and failure semantics

| File | Contract |
|---|---|
| `property_filter_report.json` | Strict JSON with properties, all alerts, ADMET, ChemFH, compact synthesis assessment, identities, pose lineage, coverage, backend provenance, and per-record `gate`/`classification` |
| `property_filter_summary.json` | Input/valid/invalid/unique counts; pass/flag/fail and keep/modify/reject counts; ADMET, ChemFH, and synthesis coverage |
| `admet_ai/admet_worker_output.json` | Exact 41-endpoint scorer payload per complete candidate, domain-validation diagnostics, separate unregistered provider evidence, and provider schema provenance |
| `chemfh/chemfh_worker_output.json` | Strict candidate-level seven-endpoint predictions, uncertainties, 17 rule results, and pinned asset provenance |
| `synthesis_report.json` | Top-K selection policy plus full AiZynthFinder route trees and starting materials |
| `prefiltered_library.sdf` | Optional RDKit prefilter survivors |
| `property_prefilter_manifest.csv` | Optional keep/remove audit trail |

Zero records and batches with no RDKit-processable molecule return
`success=false`. Mixed valid/invalid batches continue and retain a failure row
for every invalid record. Output serialization rejects NaN and Infinity.

ADMET, ChemFH, or retrosynthesis execution failure does not invent favorable
evidence, and does not charge the outage to the molecules either. A backend that
was planned and failed as a whole is declared on the run: `steps_skipped` plus a
per-candidate `metric_context.operational_errors` entry naming the tool, which
keeps the metrics it owns inside `Σ_planned` so coverage falls and the failure
registry names it. Only findings about the molecule
itself — an RDKit property review, a positive ChemFH call, a retrosynthesis
search that ran and found no route — reach `review_only`.
The global runtime preflight normally catches a missing enabled environment or
model before Step 5 starts.

Only the stable RDKit/Lipinski/QED/PAINS-Brenk-NIH/SA/ChemFH-probability/
ADMET-AI/AiZynthFinder subset defined above enters property metric scoring.
Unregistered ADMET-AI values, ChemFH native structural rules, and dedup/backend
metadata remain evidence, and absent third-party tools never acquire fabricated
values or contributions. See `../METRIC_POLICY.md` for the exact executable
inventory and interpretation choices.

## References

- `references/chemfh.md` — endpoint calls, uncertainty, rule overlap, and pinned assets

## RAscore

The optional `rascore` tool runs the installed CPU XGB model in its own
`${SBDD_ENV_RASCORE_PREFIX}` runtime, joins by candidate ID, and enters SYN before
the category fold (hard/soft lower 0.2/0.5, weight +3). This predicts accessibility;
it does not return an experimental or computed synthesis route. Failed selected
predictions retain their errors and remain NA in the existing SYN scoring and
coverage calculation, without an additional nomination veto unless the resolved
plan explicitly requires `rascore` or `rascore.score`. Enable it explicitly in
either property stage; see [metric extension controls](../pipeline-entry/references/metric-extensions.md).

The adapter retains the original SMILES and candidate ID, and creates a separate
`model_smiles` for this installed chirality-independent ECFP6-counts model. Stereo
is removed only from the model copy; isotopes, charges, the source library, and
docking chemistry are retained. This bridges newer non-tetrahedral SMILES and
the older RAscore RDKit parser. Requests/results record both representations and
`model_representation`; joins validate both. Do not reuse this transformation
for a model with chirality-dependent features. Invalid inputs still produce
per-molecule failures and null scores.

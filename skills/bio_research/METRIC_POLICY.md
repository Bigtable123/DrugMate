# SBDD metric gate and weighted-ranking policy

This document describes the implemented evaluation policy: 142 candidate rules
and 3 protocol rules (145 total), using 20 metric categories. The
[executable metric tables](shared/metric_scoring.py) contain the current bounds
and weights, which are also exposed in the
[generated catalog](pipeline-entry/TOOL_CATALOG.md). Threshold and weight roles
can overlap.

Do not register placeholder computations to fill unavailable measurements. Optional
OpenFE/TIES planning outputs, for example, are not measured free energies.

## Gate and score contract

The gate has two levels. Level 1 marks one metric on one molecule; level 2 folds
those marks per parent metric category and is the only place a molecule is
decided.

**Level 1 — per-metric mark.** Gate raw values before normalization.

- `value < hard_lower` or `value > hard_upper` is `fail`.
- A value outside the soft interval, but on or inside the hard boundary, is
  `flag`.
- A value on or inside the soft interval is `pass`.
- A missing or non-finite value is `not_applicable` and is not counted.
- A missing lower or upper side is negative or positive infinity.
- A mark is never a molecule decision. `fail > flag > pass` only when marks on
  the *same* metric are reduced.

**Level 2 — per-category fold.** Once a parent category has been computed for a
molecule, count its marks: `count(fail)`, and
`ratio(fail) = count(fail)/(count(fail)+count(flag)+count(pass))`. The active
screening mode turns those two numbers into the molecule verdict:

| mode | name | eliminates | rejects when |
| --- | --- | --- | --- |
| `V1` | marking | no | never; marks are recorded only |
| `V2` | standard (default) | yes | `count(fail) >= 3` or `ratio(fail) >= 50%` in any one category |
| `V3` | strict | yes | `count(fail) >= 2` or `ratio(fail) >= 30%` in any one category |

- Any category trips the mode → `reject`.
- Otherwise, every evaluated category has `count(pass) / count(labeled) > 50%` → `keep`.
- Otherwise → `modify`; a tie at 50% is `modify`.
- Nothing was marked at all → `not_evaluated`; absence of evidence never
  becomes `keep`.
- The count and ratio branches are independent alternatives and neither is
  qualified by a minimum fail count. The rejection formula is
  `count(fail) >= N` **or** `ratio(fail) >= P%`, so
  one `fail` among a handful of labeled metrics is ratio 100% and does reject.
  Missing values are excluded from both numerator and denominator. No
  observed-coverage guard modifies this formula. Coverage only controls rank
  publication; NA itself never becomes a fail.
- The mode defaults from candidate provenance when it is not named explicitly:
  freshly generated molecules get `V3` because the run can
  simply generate more, a pipeline-sampled library gets `V2`, and a set the user
  supplied in full gets `V1` so it is marked but never thinned. An explicitly
  named mode always wins.
- Only structural defects reject outside this fold: an unreadable identity or
  parse failure, and a PoseBusters failure not attributable to any registered
  check. Every threshold breach is a level-1 mark.
- Protocol gates are the deliberate single-level exception. Redocking RMSD and
  the retrospective screening benchmark invalidate a *run*, not a molecule, so
  no category fold applies to them.
- Build one scientific cohort from `protein_state_id + pocket_id + round_id +
  protocol_profile_id + metric_catalog_version`. Nomination candidates and
  `DEC_*` benchmark decoys use separate percentile pools inside that cohort.
  A caller-provided `cohort_id` is retained as declared provenance but never
  replaces those five components in the actual partition key; reusing a label
  across different components therefore cannot merge incomparable records.
- Remove fold-rejected and protocol-blocked records before every percentile
  pool. Their metric/category/step/overall scores and ranks remain null.
- For each finite metric, calculate descending `RANK.AVG` and
  `RP=(n-rank_avg)/(n-1)`; one valid value receives `RP=0.5`. Ties retain their
  average rank. Weight-zero metrics still mark and expose RP, but contribute
  zero.
- Explicit target-range metrics first map the hard/soft interval to the
  implemented trapezoid desirability `U`, then rank `U` descending. The
  implemented set is ExactMolWt, cLogP, TPSA, formal charge, and heavy-atom
  count, plus ADMET-AI Lipophilicity and VDss. Hard gates still inspect the raw
  value first.
- Missing is not zero. A missing metric leaves that candidate's category
  numerator and denominator; it never disables the metric for candidates that
  do have a value. An explicit `not_applicable_metric_ids` declaration takes
  precedence over any stale scalar carried by the record and removes that
  metric from the candidate's plan as well as its score.
- Within category `c`, preserve the catalog sign with
  `Score_c=100*(sum(w*RP)+sum(max(-w,0)))/sum(abs(w))`, using only values
  available for that candidate. Rank each finite `Score_c` again with the same
  descending `RANK.AVG` rule.
- Compute StepScore and OverallScore from category percentiles and positive
  parent-category weights. Parent weights follow the overview table:
  `PHY=3, ALERT=1, SYN=2, ABS=4, DIST=3, MET=3, EXCR=2, TOX=5, AFF=5, DRAW=3,
  FLEX=3, INTRA=4, INTER=5, INTX=5, POCKET=5, DYN=3`; implemented categories
  absent from the selected profile are not planned. `IDQC`, `CHEM`, and `PBQ`
  carry parent weight `0`: they mark and gate but never enter a step score.
- Rank each finite StepScore descending within the same scientific cohort and
  normalization pool, using gate/protocol-eligible records only. The local
  average rank is `step_scores.<step>.rank_avg`; ties retain their RANK.AVG and
  a one-record step pool receives rank `1`. Step rank is descriptive and does
  not feed OverallScore. Final machine outputs also expose the flat aliases
  `step4_score/step4_rank` through `step6_score/step6_rank`.
- Report planned weighted coverage, catalog weighted coverage, and parent
  category weighted coverage separately. Catalog coverage always uses original
  catalog weights, not session overrides or family-cap scaling.
- Missing planned metrics remain NA and are reported explicitly. Non-zero
  metrics affect release only through the pipeline-local PHY/ABS/TOX
  planned-coverage floors and the parent `0.50/0.10` bands; a single ordinary
  missing endpoint does not invent a new threshold.
- A missing bounded weight-zero gate blocks nomination when the plan declares
  that evidence required. ChemFH follows the explicit `required_tools` /
  `required_metric_ids` contract: an optional failed run is an audit gap, while
  a required failed run soft-modifies protocol eligibility and triggers
  retry/fallback. Legacy rows with no required/optional distinction remain
  conservative.
- When `planned_metric_ids` is present it is authoritative for the planned
  coverage denominator; `required_metric_ids` is an additive required subset,
  and `not_applicable_metric_ids` is removed. A valid extra value may remain
  catalog/scoring evidence, but it is not retroactively called planned.
- Parent coverage `>=0.50` releases a formal rank; `0.10 <= coverage <0.50`
  releases a provisional rank; `<0.10` leaves OverallScore/Rank null and asks
  for a backtrack that refills the missing categories. These coverage outcomes
  are carried in `coverage_gate`/`rank_status` and do not rewrite the current
  chemical `classification`; coverage and tool problems must never reject a
  molecule or a pose. OverallRank is descending `RANK.AVG` over
  only score-released, fold-eligible records.
- Session-local metric and category overrides accept magnitudes `0..5`.
  Metric overrides preserve the original sign/direction and neither kind can
  disable a hard gate.
- `weighted_total_score`, `weighted_score_fraction`, and
  `weighted_score_rank` are deprecated aliases of OverallScore,
  `OverallScore/100`, and OverallRank during schema transition; `rank` mirrors
  OverallRank and remains null for ineligible records.
- Machine outputs carry both vocabularies during the schema transition: `label`
  / `molecule_label` hold the level-1 mark (`pass`/`flag`/`fail`), while `gate` /
  `molecule_gate` keep the legacy `keep`/`modify`/`reject` spelling for readers
  that predate the fold. Only `metric_gate` — the fold result — is a molecule
  decision.

## Docking candidate-score and benchmark policy

- Candidate-level docking metrics use the same two-level gate before scoring:
  values outside a hard bound are marked `fail`, values outside a soft bound but
  inside the hard bounds are marked `flag`, values inside the soft interval are
  marked `pass`, and the AFF/DRAW/INTER/... category folds decide the molecule.

- Percentile scoring and ranking are allowed only within one scientific and
  metric-protocol cohort.
  The protocol identity covers the engine/scoring function and version, the
  content identity of the receptor representation actually used by that engine,
  the ligand batch and optional preparation manifest, search box, and material
  search settings. File paths remain audit metadata rather than protocol-hash
  inputs. An unknown engine version or incomplete content identity blocks trusted
  ranking. Each active scalar forms its own descending RANK.AVG percentile
  pool; lower-is-better catalog metrics retain a negative signed weight.
- Raw Vina-family energy signals share an effective family-weight cap of `5` in
  the candidate composite. GNINA, Uni-Dock, and AutoDock Vina must not each gain
  a full independent contribution for homologous energy terms; the cap limits
  family weight, not the raw kcal/mol value.
- Records from different protocol identities must not receive one shared rank.
  A mixed-protocol batch blocks composite ranking until it is split into
  comparable cohorts or explicitly rerun under one protocol.
- Redocking and benchmark results are protocol gates. They can invalidate or
  request review/rerun of a docking protocol, but they do not reject candidate
  molecules and do not contribute to the candidate weighted score. Because they
  judge a run rather than a molecule, they stay single-level: their `fail` is
  read directly and never enters a category fold.
- A retrospective ACT_/DEC_ benchmark is valid only when the selected native
  docking metric covers every labeled record and every selected score carries
  one shared, complete protocol identity. Missing coverage or identity is a
  protocol `modify`, never a silently reduced benchmark subset.
- AutoDock Vina/GNINA local pose rank is an export gate, not a molecular score:
  rank `>20` drops that pose from the export window, maps the parent molecule
  to `modify`, and permits an in-window alternate pose or rerun. Like the
  protocol gates it bypasses the fold, because the subject of the decision is a
  pose rather than the molecule.
- The supplied table encodes lower-is-better importance as a negative suggested
  weight. The public policy exposes both the nonnegative magnitude and
  `signed_weight`; scoring uses the signed value.

The pipeline does not register or run smina, QuickVina, AutoDock4Zn, or covalent
AutoDock4/ADFR protocols. No covalent-geometry or pre-reaction-pose scoring rule
is registered. AutoSite is registered as a pocket detector for memory-guided
selection: it supplies native cluster ranks, support points, and centroids.
Its detector output does not add AutoSite scoring metrics or change the 145-rule
evaluation scope. POCKET scoring uses `pocket_occupancy`; unavailable metric
values must not be fabricated.
ADFR's *primary score* is registered; only its flexible-receptor sampling is
not, which the rule records as an applicability caveat rather than a missing
rule (see below). The registered rigid-receptor, non-covalent engines use these
primary-score rules:

- AutoDock4 estimated free energy of binding: `soft_upper=-7`,
  `hard_upper=-1`, weight `4` — a relaxed Vina-like fallback,
  used only as a coarse screen alongside within-batch
  rank.
- rDock total score: no absolute bound, lower-is-better within one protocol,
  weight `5`.
- ADFR affinity: `soft_upper=-6`, `hard_upper=0`, weight `4` — the AutoDockFR
  bounds configured here (signed weight `-4`: magnitude 4,
  lower-is-better). The registered runtime is rigid-receptor ADFR while that row
  asks for a fixed flexible-residue set, so the value is read as a same-protocol
  rank and a coarse sanity bound rather than a flexible-receptor result. The
  weight and bounds are not withheld on that basis: the applicability caveat is
  carried in the rule description, not by zeroing a catalog weight.
- DiffDock confidence and DOCK 6 grid score have no row in the supplied table.
  Both remain same-protocol audit signals with `weight=0` unless a session-local
  administrator/user override explicitly activates them.

Finite no-bound scores remain available for same-protocol ranking but their gate
label is `not_applicable`: without a declared bound they must not manufacture a
`pass` or dilute `ratio(fail)`. Missing/non-finite primary scores are blocked by
the adapter and coordinator contracts. A registered
engine with no canonical nonzero metric cannot by itself produce a formal
overall score without an explicit weight override; this is reported as
coverage, not hidden behind an invented catalog weight.

Vina `vina_inter` and `vina_intra_best_pose` are weight-zero gates with
`hard_upper=0`: a positive component energy is physically
implausible and earns a `fail` mark, but it carries no weight and, like any
other mark, only matters through the DRAW fold. `vina_torsions`, rDock
component scores, a Uni-Dock local `engine_pose_rank`, a hashed
`score_protocol` identity, plus lineage/completeness, search-box, runtime, and
GPU diagnostics have no executable scalar policy and remain audit evidence
rather than implicit weighted components.

The current redocking implementation is GNINA native-ligand redocking and
reports the best symmetry-aware RMSD over generated modes. It does not implement
Vina mode `RMSD l.b./u.b.`, top-N aggregation, independent-replicate
consistency, cross-docking RMSD, BEDROC, or a true cross-docking decoy-bias
measurement.

## Executable table subset

The current pipeline produces enough structured data to implement these
families now:

- GNINA `minimizedAffinity`, `CNNscore`, and `CNNaffinity`, plus weight-0
  `intramol`, `CNNvariance`, and local pose rank when present.
- Uni-Dock `docking_score` and AutoDock Vina `vina_affinity`; AutoDock Vina also
  emits derived heavy-atom ligand efficiency, local pose rank, and score gap.
- DiffDock confidence and DOCK 6 grid score are same-protocol audit-only signals
  without universal absolute cutoffs. rDock remains a signed-weight
  same-protocol score, while AutoDock4 and ADFR each carry the documented
  fallback energy bounds and a nonzero weight. All eight registered
  primary engines participate in mixed-engine/protocol and missing-primary-score
  blockers.
- GNINA best-of-generated-modes native-redocking RMSD and retrospective
  ACT_/DEC_ AUROC/EF1 protocol gates. PR-AUC/prevalence and decoy fractions are
  audit outputs without universal cutoffs.
- PoseBusters `dock` checks that are present in the installed configuration,
  plus an overall boolean success ratio.
- PoseCheck strain energy and heavy-atom steric-clash count.
- PLIP interaction count with a hard lower bound of `2` and weight `0`; the
  previous synthetic sub-1.36 Å counter is not an executable catalog rule.
- ProLIF interaction counts and reference-fingerprint similarity are audit-only;
  failed/error results are NA rather than zero-interaction observations.
- PLIP-based Blueprint verdict, required/desired diagnostics, desired ranking score,
  desired evaluable coverage, and optional hits. When a blueprint declares no
  optional constraints (`optional_total=0`), optional hits are NA and unplanned
  rather than a weighted zero.
- Pocket contact fractions for the four regions and overall occupancy.
- The current single-frame Amber MMPBSA
  `amber_mmpbsa.delta_total_no_entropy` signal, mapped to the table's weak
  MM-GBSA sanity gate with its low-reliability provenance retained.
- RTMScore `rtmscore.ensemble_score`: DYN weight +3; no absolute gate, same-protocol percentile only.
- OpenMM/MDAnalysis: symmetry-aware ligand heavy-atom RMSD median (soft/hard
  2/5 Å), minimized over element/Amber-type topology automorphisms after
  protein alignment without a ligand-only refit; P95
  mass-weighted ligand–pocket COM-distance drift (2/5 Å), and mean initial
  heavy-atom contact retention (hard/soft 0.2/0.5). The pocket is the fixed set
  of protein heavy atoms within 6 Å of the reference ligand; drift is the
  absolute change from the equilibrated production-start COM distance. The
  separate pre-minimization docked structure remains provenance only. P95,
  median, slope and per-frame CSV remain in the report. Unselected top-K poses
  are unplanned for these metrics.
  Short MD is a screen, not a demonstration of equilibrium binding stability.
- RDKit: MMFF94s docked strain relative to the lowest converged conformer among
  eight independently embedded ETKDGv3 samples; the docked pose is not admitted
  into the reference ensemble, and its heavy atoms stay fixed during hydrogen
  relaxation. TFD and symmetry-aware heavy-atom RMSD use that
  sampled reference (not a claimed global minimum). No TFD for torsionless
  molecules; metal/unsupported chemistry is NA. Parameters and seed are recorded.
- FreeSASA: heavy-atom RDKit van der Waals radii, 1.4 Å probe. Ligand-side buried
  surface is free ligand SASA minus its SASA in the receptor complex.
  `contact_surface_area` independently follows the table formula
  `0.5 * (free receptor SASA + free ligand SASA - complex SASA)`; each fraction
  divides its own area by free ligand SASA. The two correlated rows keep table
  weights +3/+1 and their declared bounds. Free/bound receptor and ligand SASA,
  complex SASA, and the absolute Å² values are retained for audit.
- Arpeggio: candidate-specific derived mmCIF, explicit atom/bond metadata,
  chemistry-quality checks, numeric clash/vdW-clash/proximal-fraction gates.
  Raw contact records are retained; degraded chemistry leaves metrics missing,
  never a zero clash count. The derived complex is labeled as such.
- ProLIF: residue/type fingerprint Tanimoto to the same-receptor reference,
  retained only for audit and independent Blueprint cross-check.
  No reference or an empty reference fingerprint is NA. The FP row has weight
  zero and no gate bounds. PLIP/ProLIF interaction-type
  counts, LUNA hashed fingerprints, and ODDT IFP/SPLIF remain audit data; raw
  unbounded records do not add artificial pass marks to the fail-ratio denominator.

### Property-filter stable subset

Step 3/6 first preserves the legacy RDKit and backend evidence, then applies
the same raw-value marking pass and cohort scorer to the stable scalars that
the current pipeline actually emits, and lets the PHY/ALERT/SYN/ABS/DIST/MET/
EXCR/TOX folds decide the molecule. The executable property subset is:

- RDKit descriptors: ExactMolWt, cLogP, TPSA, HBD, HBA, rotatable bonds, total
  and aromatic ring counts, FractionCSP3, formal charge, heavy-atom count, and
  chiral-center count. cLogP and TPSA have weight `3`; the remaining descriptor
  rules have weight `0` and therefore mark/audit without entering the total.
- Lipinski RO5 violation count (weight `0`) and RDKit QED (weight `3`). The
  default QED bounds are `hard_lower=0.10`, `soft_lower=0.45`.
- PAINS A/B/C, Brenk, and NIH alert counts, each weight `1` inside ALERT, with
  matched rule names and atoms retained. All three are counts under the same
  relaxed band: `0` is a `pass` mark, `1-4` a `flag` mark, and
  only `>4` a `fail` mark, so a single structural-alert hit no longer ends the
  molecule anywhere in the pipeline. RDKit reports PAINS_A/B/C as three separate
  catalogs behind the one `pains_alert_count` metric, so their hits are summed
  into a single family count before the bound is applied. Covalent status alone
  never rescues a hit; only an explicit per-catalog hard-fail opt-out drops a
  catalog to flags, while preserving the evidence.
  RDKit FilterCatalog does not publish the table's normalized Brenk severity,
  so the proposed "one severe hit" shortcut is not fabricated; the rule name
  and matched atoms remain available for review.
- ChemFH's seven exact model probabilities are ALERT gates with
  `soft_upper=0.30` and `hard_upper=0.90`. Reactive compounds carry weight `1`,
  promiscuous compounds and colloidal aggregators `0.5` each, and the four
  fluorescence/interference endpoints weight `0`. ChemFH's
  ten-rule substructure alert count is registered separately with
  `soft_upper=2`, `hard_upper=10`, weight `2`, matching the table's
  `none=0 / mild=2 / severe=10` encoding. Their uncertainty fields remain audit
  evidence; native PAINS overlap is not counted as a second structural-alert
  gate.
- RDKit SA score with `soft_upper=6`, `hard_upper=9`, weight `2`.
- ADMET-AI's 41 provider-native endpoints registered in the executable tables,
  plus deterministic `CYP Substrate count = n_sum(p<0.5)` and
  `CYP inhibition count = n_sum(p>0.7)`. They populate ABS, DIST, MET, EXCR,
  and TOX using the provider's exact endpoint IDs. Lipophilicity and VDss use
  the target-range transform; BBB defaults to the documented non-CNS profile,
  while explicit `admet.params.bbb_profile=cns` reverses both its hard/soft
  gate and signed scoring direction.
  The three provider-native excretion rows have blank bounds/weight zero and
  are audit-only, so they do not activate a weighted EXCR score by themselves.
  CYP2C9/2D6/3A4 substrate probabilities use the explicitly selected negative
  signed weights `-0.1/-0.1/-0.5`, so higher probability is penalized. This is
  the implemented signed scoring convention. The generated tool catalog
  carries every effective bound, signed
  weight, unit-oriented description, and category.
- AiZynthFinder route-found, state score, reaction count, precursor count, and
  an explicit `precursor_count == precursors_in_stock_count` completion scalar.
  Only state score has nonzero weight (`2`) and is
  emitted only when the route analysis provenance names exactly `state score`;
  it is comparable only within one fixed AiZynthFinder
  version/config/search-budget cohort. Route/search failures and overlong
  search failures are `flag` marks, not eliminations. A reaction
  count `<=5` is `pass` and `>5` is `flag` under the canonical numeric cells
  (`soft_upper=5`, with no hard-fail bound), so reaction count alone never adds
  a SYN `fail`. Price sum and route cost remain NA because the registered worker has
  no price/cost source.

Each applicable scalar uses descending RANK.AVG percentile inside the property
cohort, after the explicit target-range transform where required. Weight-zero
values may still change the mark but contribute exactly zero. Missing/NA values
leave only that candidate's numerator and denominator and lower its planned
coverage; complete candidates retain the metric. Step 4 category/overall fields
are authoritative. `property_weighted_total_score`,
`property_weighted_score_fraction`, and `property_weighted_score_rank` remain
deprecated aliases during schema transition. The more severe of the legacy and
scalar marks is materialized before percentile pools are built, and the legacy
RDKit screen's own `fail` verdict is capped at `modify` because it is a mark
rather than a fold result.

ADMETlab 3.0, SwissADME, RAScore, SCScore, SYBA, ASKCOS, and pkCSM are
not registered or executed by this pipeline, so their rows are skipped rather
than synthesized from similarly named fields. ADMET-AI is an additional
pipeline backend, not an alias for any of those providers. Only the 43 exact
ADMET-AI rows described above are registered; thresholds from ADMETlab,
SwissADME, or pkCSM must never be applied to its remaining endpoints by name
alone.
ChemFH is registered separately for frequent-hitter/interference evidence. Its
seven exact table-matched endpoint probabilities are executable ALERT gates:
`<=0.30` is a `pass` mark, `0.30-0.90` a `flag` mark, and `>0.90` a `fail`
mark; reactive compounds weigh `1`, promiscuous compounds and colloidal
aggregators `0.5` each, and the four fluorescence/interference endpoints `0`.
Their uncertainties remain audit context. The native structural rules
(including overlapping PAINS) are aggregated into the single
`chemfh.substructure_alert_count` gate that mirrors the table's
`none=0 / mild=2 / severe=10` encoding; individual rule matches make no further
gate or weighted contribution. Because exactly ten rules are mapped, its
`hard_upper=10` is only reachable after a catalog-version change, so in practice
it contributes flags and weighted rank rather than `fail` marks.

## Supplied table rows not executable in this pipeline

These rows are intentionally not fabricated when no corresponding candidate
scalar exists:

- every ADMETlab 3.0, SwissADME, SCScore, SYBA, ASKCOS, and pkCSM
  row, because those tools are absent from the registry/runtime;
- explicit RDKit canonical-SMILES/InChIKey/atom-map/stereo/tautomer/valence/
  formal-charge consistency booleans;
- candidate-topology conformer coverage and undefined-stereocenter count;
- PoseBusters molecular formula, molecular bonds, double-bond stereo,
  tetrahedral chirality, and benchmark RMSD under the installed `dock` config;
- PoseCheck close-contact minimum distance, interactions, interaction
  fingerprint, and redocked RMSD;
- candidate LUNA/ODDT normalized scalar gates (their per-pose contacts and
  fingerprints are retained as weight-zero audit artifacts);
- ProLIF trajectory interaction frequency;
- pocket centroid distance, volume occupancy,
  protrusion, unsatisfied polar atoms, Lawrence-Colman shape/electric
  complementarity, and a separate key-residue coverage scalar;
- GROMACS trajectory metrics (no GROMACS short-MD backend is registered);
- GenScore;
- MM-GBSA decomposition, SD/SEM, and block convergence;
- computed OpenFE/TIES free energies, uncertainties, overlap/convergence,
  hysteresis, cycle closure, and replica consistency; current OpenFE/TIES
  routes are planning/preflight only;
- PMX, which is not registered.

## Pipeline tools and signals missing from the supplied table

The table does not define thresholds/weights for these existing pipeline
capabilities, so they remain lineage, audit, preflight, or legacy primary-score
evidence rather than implicit weighted components:

- Step 4: `identity_qc`; PoseBusters molecule-loading and `no_radicals`
  checks; PoseCheck all-atom clash
  diagnostics; PLIP/ProLIF capability and disagreement evidence;
  candidate-topology chemistry/charge/topology checks; OpenFE mapping/LoMAP
  diagnostics; TIES planning diagnostics; MM-GBSA topology consistency and
  unweighted energy terms.
- Step 5: the ADMET-AI endpoints without an exact registered provider row,
  along with all ADMET-AI task/unit metadata; the ChemFH endpoints without an
  exact table row, every endpoint uncertainty, the individual rule matches
  behind `chemfh.substructure_alert_count`, and pinned assets; identity-aware
  `dedup`, backend completeness, and route-tree provenance. Registered
  ADMET-AI rows still never borrow ADMETlab/SwissADME/pkCSM semantics; dedup
  and provenance are not scalar weighted components.
- Upstream pocket analysis: ligand-site, fpocket, pyKVFinder, P2Rank,
  PocketMiner, blueprint generation, and reference-only
  Arpeggio/LUNA/ODDT/ProLIF profiling.

## Explicit interpretation choices

- Numeric cells reproduced in `shared/metric_scoring.py` and the generated
  catalog are authoritative where a prose sentence
  conflicts with them. For `pocket_occupancy.overall`, this means
  `hard_lower=0.05`, `soft_lower=0.20`, weight `3`.
- Where PoseBusters prose says reject but the numeric table supplies only a
  soft lower bound (organic/inorganic cofactor and water context), the numeric
  columns are used: failure defaults to a `flag` mark. A PoseBusters failure
  that cannot be attributed to any registered check stays a structural hard
  reject, because it is an unreadable result rather than a threshold breach.
- Blueprint desired-evaluable coverage uses only `soft_lower=0.80` from the
  numeric table. There is no separate hard lower bound at `0.30`.
- V2/V3 use the count OR ratio formula exactly, including a ratio reached by
  one fail. Missing planned values leave the fail ratio; the additional 0.50
  observed-coverage guard has been removed. Rank-release coverage remains
  separate. Parent category weights have not changed.
- Pocket region occupancy has weight zero but retains its numeric hard/soft
  bounds and contributes fail labels to its category. Zero score weight does
  not erase a numeric gate.
- Reference-dependent metrics remain not applicable when their reference is
  absent.
- A non-pass RDKit first-pass verdict is a molecular finding only when the screen
  actually marked this molecule. `compute_verdict` publishes the split in
  `rdkit_evidence`: `molecular_findings` (alerts present, Lipinski violations
  above the max, QED below threshold, SA above threshold) versus `coverage_gaps`
  (`no_check_enabled`, `alert_screening_incomplete`, `sa_score_unavailable`), and
  a raised RDKit pass is already reported as `property_tool_error`. Gaps reach
  `protocol_review_reasons` and cap nomination at `review_only`; they never enter
  `molecule_review_reasons`, never drive `legacy_gate` to `modify`, and are
  demoted to `not_evaluated` rather than `modify` in the coordinator's combined
  classification. Rows carrying no `rdkit_evidence` keep the previous reading.
- RO5 uses `soft_upper=1`, `hard_upper=3`:
  `0-1` is a `pass` mark, `2-3` a `flag` mark, `>3` a `fail` mark. It remains
  weight `0`. `lipinski_max_violations` configures the hard bound only and can
  never fall below the soft bound.
- Chiral-center count uses the numeric `soft_upper=2`, weight `0`, and no hard
  upper bound: `<=2` is a `pass` mark and `>2` a `flag` mark, so the count alone
  can never contribute a `fail` to the PHY fold.
- Vina ligand efficiency is emitted by the runtime as positive
  `-vina_affinity/HAC`; its equivalent marks are `<0.05` `fail`, `0.05-0.25`
  `flag`, `>=0.25` `pass`, with positive runtime signed weight `+2` equivalent
  to the table's negative raw-value weight.
- PoseCheck strain uses `soft_upper=8`, `hard_upper=60`, signed weight `-4`, and
  its steric-clash count `soft_upper=0`, `hard_upper=2`, signed weight `-5`.
  Both are INTRA marks: `>60 kcal/mol` and `>2` clashes constitute a *pose*
  failure, and the INTRA fold is what turns such marks into a molecule verdict.
- The current redocking adapter is GNINA-specific although the canonical row is
  named Vina. Its numeric `2.0/2.5 Å` protocol gate is retained and labeled with
  the actual provider; it is not presented as semantic equivalence.
- AFF and DYN use parent weights `5` and `3`.
- ALERT and SYN use parent weights `1` and `2`.
  ALERT is deliberately the lowest-weighted scoring category:
  structural alerts are context-dependent, so they earn marks and review notes
  rather than a large share of the composite.

## Evidence-plan ordering

Resolve planned/available/missing metric states before the molecular fold.
The inferred plan uses structurally evaluable cohort inputs and is then fixed
for this scoring pass; it must not depend on the survivor population it helps
classify. The same assessments determine both percentile eligibility and the
published metric gate. Rejected records receive no score and enter no ranking
pool. Protocol comparability is checked after molecular eligibility.

### Newly available switches and remaining high-cost endpoints

`pose_evaluation.tools` accepts `rdkit_pose`, `freesasa`, `arpeggio`, `luna`,
`oddt`, and the existing `prolif`, `rtmscore`, `openmm_md`, `mdanalysis`.
`property_prefilter.tools` and `property_filter.tools` accept `rascore`.
New switches default off and must be explicitly selected in a complete pipeline
plan; enabling a tool makes its produced scalars available to the shared rules.
For MD, enable `candidate_topology`, `openmm_md`, and `mdanalysis` together;
existing top-K and trajectory-length controls remain in force.

RAscore uses the installed XGB ChEMBL/ECFP model on CPU: SYN hard/soft lower
0.2/0.5, weight +3. It predicts accessibility and does not demonstrate a route.
A failed selected worker is planned-missing and blocks automatic nomination.
The shipped model is a ChEMBL-region proxy and its authors do not claim validated
generalization to generated/exotic chemistry, so generated molecules need an
applicability-domain warning and route-search confirmation rather than treating a
high RAScore as proof of synthesis.

Per-metric workers only emit labels. Category folds run after the category's
planned evidence reaches a terminal available/NA state; in particular SYN is not
folded before AiZynthFinder finishes. Final NA values leave both fail-ratio terms,
so denominator shrinkage may still make the final V2/V3 ratio reject. Pending
values are not final NA. Finite rows with no hard/soft bound (including RTMScore)
are score-only and do not add an automatic pass to that denominator.

Not added: OpenFE/TIES actual RBFE simulations or numerical free energies;
MD H-bond occupancy/interaction frequencies without predeclared key contacts;
ΔRMSF without an apo/reference trajectory; restraint occupancy without a
distance specification; composite short-MD stability without an explicit
condition set; multi-frame MMGBSA/SEM/convergence. OpenFE/TIES remain
preflight/plan only. See `audit/integration_20260910/INTEGRATION_REPORT.md`
for tested coverage and the remaining metric table.

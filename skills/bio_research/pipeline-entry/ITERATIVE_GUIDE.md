# SBDD Pipeline — Agent Decisions Across Rounds

Use this guide after the current round and its applicable final gate are complete.
First read the [execution contract](SKILL.md#user-confirmation-and-execution-modes-mandatory-contract)
and the [campaign state, planning, and authorization protocol](references/campaign.md).
`run_pipeline.py` executes one plan per invocation; the coordinator assesses the round,
and the Agent decides the next action.

## Read the Current Round's Evidence

Read `round_assessment`, quality gates, `candidate_ranking.json`, category
score/rank/status/coverage, and failure reasons from `step6/coordinator_report.json`
and associated outputs. New reports use `decision="agent_review"`,
`decision_outcome=null`, and `goal_met=null`. Do not wait for `iterate` before planning.
A KEEP fraction of 30%, a count of 3, or `steps.coordination.params.max_rounds` is
not a campaign termination condition. `valid_for_nomination` remains constrained
by deterministic round evidence; the Agent must not set it to true independently.

After checking the round's outcome, read this project's round reports, configs,
and artifact lineage directly. Preserve unresolved conflicts and the conditions
under which failed parameters were used. Before each new plan, review pocket,
generation, and docking under the [unified memory protocol](references/benchmark-memory.md).
Retain declared modes, metrics, presets, and index versions. Requery affected domains
when inputs, objectives, available tools, or versions change; cite unchanged results
only when their provenance is verified. `context.docking_memory_mode=off` reads no
docking history. Pocket memory has its own `pocket_memory_mode`; compatible history
drives reliability fusion, while unavailable/off memory uses fixed DeepPocket#1
rank fusion. Legacy recall does not replace this selector. Record suggestions,
actual choices, and override reasons in this round's `memory_consultation.md`.
Do not mix history from other ablation arms or overwrite earlier records.
Memory neither extends round budgets nor replaces the next Gate A.

Treat insufficient evidence, coverage gaps, and tool failures first as protocol/execution
issues to repair. Judge molecules under the [metric policy](../METRIC_POLICY.md),
metric catalog, and V1/V2/V3 formulas. Retain the 2026-09-10 interpretation: numerical
formulas take precedence; do not prematurely collapse pending values into NA;
final NA is excluded from failure-fraction denominators; weight 0 does not disable
an explicit numerical gate. Cross-round planning does not change scoring tables.

## R1 → R2: Replication and Condition Selection

R1 runs protein preparation, pocket analysis, molecule preparation/prefiltering,
docking, and post-docking evaluation, delivering candidates and failure/coverage
evidence. Explicitly choose R2's purpose from the current evidence:

| Observation | Agent's plan |
| --- | --- |
| Incorrect identity, units, receptor preparation, or protocol baseline; many invalid poses | Repair preparation/protocol in a new attempt at the current stage, not by advancing to R2 |
| Broadly poor Step 5 categories and fewer than 10 KEEP candidates | Inspect protein/receptor preparation and record category distributions and evidence; 10 is an investigation cue, not a molecule gate |
| Concentrated POCKET rejections; many candidates miss the same location | Inspect pocket/blueprint and search space; rerun from the earliest invalidated step |
| Usable R1 evidence; stability needs assessment | R2 `replication`: specify new seeds, replicate ID, fixed candidates, and steps that must actually rerun |
| A specific improvement hypothesis | R2 `parameter_improvement`: preferably change few conditions; record receptor/pocket/tool/parameter changes and reasons |

`round_number + 1` does not establish independent replication. Static preparation
may be reused under lineage rules, but docking/evaluation declared as replicated
must actually rerun. Record each tool's seed; changing only `context.seed` does
not establish that every tool's seed changed. Receptor/pocket/grid/constraint changes
invalidate affected downstream evidence. Keeping gating/scoring fixed differs from
changing an experimental protocol. New conditions create a new profile; do not
silently relax hard gates to increase KEEP.

R2 comparisons include shared chemical identities, KEEP/Top-K overlap, failure rates,
coverage, benchmarks, and key-category changes. Separate sources of change when
candidates or profiles change. Do not directly subtract or jointly rank percentile/overall
scores across profiles/cohorts. Numerical ranking changes require shared candidates
reevaluated under a common comparable protocol. Report fixed-protocol replication
stability separately from evidence after condition changes; higher R2 scores alone
do not establish improvement.

## R2 → R3: Candidate Accumulation and Local Optimization

1. Combine R1/R2 protocol validity, replication stability, coverage, and candidate
   quality to select the R3 receptor/pocket/protocol. Cite evidence and tradeoffs;
   do not compare only each round's maximum score.
2. Collect all KEEP candidates from completed, evidence-valid R1/R2 attempts.
   Deduplicate by stable chemical identity, retaining source runs, identities, and
   original verdicts. Protocol-invalid records await validation and are not trusted
   KEEP. Old scores from different profiles serve only as provenance.
3. Select MODIFY parents with concrete objectives/constraints, such as better
   interactions without added toxicity/structural alerts, scaffold preservation,
   synthesis feasibility, and limited modification size. Retain REJECT diagnostics;
   do not automatically recycle them. Use molecule-generator, molecule-modifier, or both.
4. Merge unchanged KEEP with generated/modified molecules. Preserve `parent_id`,
   `parent_smiles`, modification rationale, and provenance. Record exact duplicates
   and discarded duplicate sources. Preserve scaffold/diversity considerations;
   do not select only the highest docking score within one scaffold.
5. **Start the new library at `property_prefilter` (registry order 3) or earlier.**
   Recalculate properties after identity/structure changes, then run docking, pose
   evaluation, full property filtering, coordination, and reporting. Do not bypass
   prefiltering for new molecules with `start_step: 4`. Inject its `library_manifest`.
   When generated/modified molecules are included, explicitly use source
   `molecule_generation` (V3), not the external fixed-library V1 default.
6. Reevaluate under the selected common protocol. Report parent–child pairs, property
   changes, strengths/weaknesses, and Pareto-nondominated candidates. Do not combine
   each round's winner into one overall ranking. Mark uncomputed diversity, Pareto
   analysis, or gains as unevaluated; an output column does not prove analysis occurred.

`molecule-modifier` is a limited keyword-to-SMARTS rule engine. Its
`modification_instructions` use `target_molecules`, `per_molecule`, and `suggestions`,
written by the Agent from evidence; the coordinator need not command the next round.
A request such as adding a donor toward a residue does not establish 3D directional
design. New molecules still require full validation. R3 does not enable MD/FEP;
this implementation ends at R3 delivery and does not automatically enter R4.

## Repairs and Minimum Rerun Scope

| Change or issue | Starting point and candidate handling |
| --- | --- |
| Only weights change; raw evidence is complete | Rescore and publish under a new scoring profile; do not change structures |
| Missing tools/coverage | Run missing planned tools and dependencies; do not remove candidates to improve coverage |
| Ligand/protomer/tautomer/stereochemistry changes | `property_prefilter` order 3 or earlier, then all dependencies |
| Docking box, receptor state, constraints, or seed changes | Docking order 4 and all downstream steps; earlier if receptor preparation changes |
| Pocket/blueprint definition changes | `pocket_analysis` order 1 or earlier; recompute affected downstream steps |
| Protein sequence/structure changes | `structure_prep` order 0; new branch and cohort |
| One or more poses are invalid | Evaluate preselected poses under the entry-point policy; reject the molecule only when all selected poses explicitly fail, retain diagnostics, and continue other molecules without automatic redocking |
| Protocol modify/reject | Do not change candidate structures merely because of protocol problems; repair the protocol and retain original evidence and invalidity markers |

Native-redocking RMSD, AUROC/EF, and molecule-gate numerical rules remain in the
entry-point Quality Gates and current catalogs. Maximum attempts do not imply
protocol reject; a high KEEP count does not imply protocol pass.

## Planning, Execution, and Delivery

Copy the full config from the attempt with the selected protocol. Remove the old
`project_dir` override; retain the full matrix of unchanged upstream producers and
declare reuse through `start_step`. Increment physical-directory `round_number`
for a new attempt; record scientific stage/attempt under the campaign protocol.
A local repair increments the same stage's attempt; resuming a checkpoint does not.
Set `upstream_project_dir` to the actual selected upstream attempt. The evidence
list still references the immediately preceding run and adds the basis for selecting
historical conditions. If an unfinished attempt needs a scientific config change,
preserve its evidence and use campaign `abandon` before creating a new attempt;
do not overwrite frozen config.

Use the [campaign CLI](references/campaign.md) to freeze the Agent decision and
config, run `validate_config.py`, and inspect the resolved plan. Gate A checks
existing authorization. Continue within scope without asking every round.
`interrupt` still requires per-step/final Gate B; `auto` completes the authorized
current round and does not authorize the next stage.

At delivery, present candidates, protocol limitations, failures/coverage, and
comparability with history, then the next action's rationale, parameter changes,
and estimated cost. Pause after A0 delivery; A1 defaults to pausing after R2.
A2/A3 extend at most through R3 within authorized attempts and resource budgets.
When a repair exceeds authorization, preserve reviewable results and a next-round
draft, then wait for additional authorization. No reply is not authorization.

# R1–R3 Campaign: State, Plans, and Authorization

Use this protocol together with the [iteration decision guide](../ITERATIVE_GUIDE.md).
The agent determines scientific stages, replication/improvement hypotheses, and
next steps. Code validates plans, records state, and executes a single round.
This protocol does not change the scoring, gates, or 2026-09-10 numerical
conventions. Use the executable [metric policy](../../METRIC_POLICY.md),
implemented in `../../shared/metric_scoring.py`, as the current authority.

## Stages and Attempts

| Field | Meaning |
|---|---|
| `campaign_stage` | `R1`: rapid screening; `R2`: independent replication/condition improvement; `R3`: candidate consolidation and optimization |
| `campaign_attempt` | Attempt number within the current scientific stage, starting at 1; repairs within the same stage increment it |
| `round_number` | The runner's physical directory number; it does not identify a scientific stage |
| Campaign `run_id` | A stage/attempt identifier such as `R1-a1`, distinct from the random run token in the directory name |

Resuming the same interrupt checkpoint, waiting for resources, or recovering
from an execution error does not automatically advance the scientific stage.
For a new attempt requiring changed scientific inputs/parameters, preserve old
results and files and write a new plan. Do not overwrite archived evidence.
When all poses selected for a molecule in this run have been evaluated and
explicitly fail, mark `all_poses_failed` and skip the molecule. Do not automatically
redock or create a new attempt or R2 because of this result.

When resuming a recoverable failed/paused run, the runner changes the current
campaign status and `campaign_result.json` to running and clears the old
`finished_at`. The previous execution result and status are retained in
`execution_history/`. Stage/attempt identifiers and budget counts do not change.
Completed runs remain idempotent. Runs awaiting final confirmation, already
abandoned, or superseded by a new attempt cannot be rerun through ordinary resume.

If an unfinished attempt needs a change to its scientific config, first stop
its actual process and confirm exit. The existing runner lifecycle handles
process cancellation; the campaign helper does not kill it for you. Then close
the old attempt, increment the attempt within the same stage, and write a new plan:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign.py" abandon \
  --state "${SBDD_DATA_ROOT}/_staging/request-actual/campaign_state.json" \
  --run-id R1-a1 --reason 'Record the actual reason for noncompletion and the scientific parameters that must change'
```

`abandon` operates only on the latest unfinished attempt, checks that the runner
does not hold the directory lock, and preserves existing evidence. When tools
and inputs are unchanged and recovery is possible, continue from the original
checkpoint. Do not abandon merely to reset counting.

Advance R1→R2→R3 in order. If the protocol needs repair, retry within the current
stage first. Poor R1 evaluation results do not automatically mean R1 is complete
or R2 has begun. See the iteration guide for the evidence required to choose
better conditions in R2. A1 pauses for delivery after R2; prepare R3 only when
authorization permits it.

## Authorization and Budgets

| Autonomy level | Cross-round boundary in this implementation |
|---|---|
| `A0` | Pause after each delivery; the next attempt requires new explicit authorization |
| `A1` (default recommendation) | With explicit preauthorization, proceed automatically through R2, then pause |
| `A2` / `A3` | Complete at most R3 within an explicit budget; this implementation does not execute R4–R6 |

`project.execution_mode=interrupt/auto` controls only the current round's steps.
A1/A2/A3 combined with interrupt still require step/final confirmations. For
continuous execution within the authorized scope, the user selects auto.
Do not autonomously switch execution_mode because A2 is set, or treat auto as
cross-round authorization.

Before starting, explain the objective, allowed stages, planned tool/parameter
changes, attempt limit, available resources, and cost/duration limits. Existing
explicit user authorization may be referenced directly without asking again.
`authorization_ref` records the actual user statement or a traceable conversation
reference. Config defaults, silence, a timeout default, and the agent's own
recommendation are not explicit cross-round authorization. Without authorization,
a proposed state and reviewable plan may be created, but not executed.

Current code validates stage/attempt, authorization references, attempt limits,
and the specific config hash. The agent must still enforce the semantic scope
of permitted tools/parameters and total cost/GPU/duration limits against the
user's authorization. Do not describe `max_runs` or step timeouts as complete
cost accounting. For every new attempt, state estimates and cumulative usage.
If it exceeds the user's scope, preserve the draft and existing results and wait
for additional authorization. Every A0 attempt requires separate authorization.

## Create State

Keep all experiment files under `${SBDD_DATA_ROOT}`. A new experiment may initially
store state in the existing request's `_staging` directory; subsequent configs
continue referencing the same absolute path. Do not rebuild state and discard
runs. Replace the paths and authorization text below with actual values for this
request. `max-runs 6` is an example budget that must be supported by authorization;
the example does not itself grant permission for 6 computations.

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign.py" init \
  --state "${SBDD_DATA_ROOT}/_staging/request-actual/campaign_state.json" \
  --autonomy A1 --authorized-through R2 --max-runs 6 \
  --authorization-ref 'The actual explicit user authorization or conversation reference for this request'

"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign.py" status \
  --state "${SBDD_DATA_ROOT}/_staging/request-actual/campaign_state.json"
```

Omitting the authorization reference saves intent only. Do not reinitialize
existing state. After receiving new authorization, update it in place and retain
runs. The following command applies only when the user has actually extended
authorization through R3:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign.py" authorize \
  --state "${SBDD_DATA_ROOT}/_staging/request-actual/campaign_state.json" \
  --autonomy A2 --authorized-through R3 \
  --authorization-ref 'The actual user statement or conversation reference extending authorization through R3'
```

Increase `--max-runs` only after obtaining corresponding authorization for a
larger total attempt budget. Do not clear runs or modify counters to bypass it.

## Assemble the Config and Agent Decision

Add four campaign fields under project in the complete config. Do not treat
this fragment as a complete configuration:

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "project": {
    "data_root": "<SBDD_DATA_ROOT>",
    "round_number": 2,
    "campaign_state_file": "<SBDD_DATA_ROOT>/_staging/request-actual/campaign_state.json",
    "campaign_stage": "R2",
    "campaign_attempt": 1,
    "round_decision_file": "<SBDD_DATA_ROOT>/_staging/request-actual/R2-a1-plan.json"
  }
}
```

Write a separate agent-decision draft. The following is a structural example;
actual paths must exist, and candidate IDs must come from this experiment's
evidence. An initial R1 plan may have empty `evidence_refs`. Later plans must
reference at least evidence from the immediately preceding completed/failed run
directory. `campaign_result.json` may serve as evidence of execution facts.

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "schema_version": 1,
  "stage": "R2",
  "attempt": 1,
  "purpose": "replication",
  "reason": "Check result stability with independent sampling under the same candidates and protocol",
  "evidence_refs": [
    {"run_id": "R1-a1", "path": "<SBDD_DATA_ROOT>/project-actual/round1/campaign_result.json"}
  ],
  "parameter_changes": [
    {"path": "steps.docking.tools.gnina.params.seed", "before": 42, "after": 43,
     "reason": "Independent sampling"}
  ],
  "start_step": 4,
  "candidate_ids": ["actual-candidate-identity"],
  "molecules_changed": false,
  "replicate_design": {
    "replicate_id": "replicate-2",
    "seed": 43,
    "rerun_steps": ["docking", "pose_evaluation", "property_filter", "coordination", "reporting"]
  }
}
```

Use `initial_screening`, `replication`, `parameter_improvement`,
`protocol_repair`, or `molecular_optimization` for `purpose`. An A0 decision
also contains a new `authorization_ref` for that attempt. `start_step` must
match the config; if unset, execution begins at 0. `parameter_changes` records
actual changes, old/new values, and reasons. Use an empty list when nothing changes.

R2 requires `replicate_design`, and every rerun step must be enabled. Independent
replication must actually rerun docking. On first launch, set
`project.resume_existing=false`; resuming the same checkpoint after an interrupt
retains its processed prefix. The seed must match the actual docking-tool
parameters. Reusing a seed requires an auditable basis for an independent
replicate, stated in `independent_replicate_reason`. A changed round number is
not evidence of independence.

At any stage, `molecules_changed=true` requires execution from order 3 or earlier
with prefiltering enabled. R3 also mandates prefiltering when molecule_generation
is enabled or a library is injected, even if `molecules_changed=false` was set
incorrectly. R3 forbids direct injection of prefiltered_library. An R3 merged
library is a new input; even if it contains previous KEEP molecules, identity
and the complete downstream workflow must be checked again. Every campaign
plan enables coordination and reporting, and its execution scope must include
order 7/8 to retain a deliverable assessment for the current round.
MD/FEP tools are outside the R1–R3 campaign's execution scope.

## Freeze the Plan and Execute

Write the complete config and decision draft first, then use the helper to
generate a plan bound to the config hash. A reviewable plan can still be generated
when state is proposed, its stage is unauthorized, or the attempt budget is
exceeded. `execution_authorized=false` and `execution_blockers` record blockers
at generation time. They are snapshots, not replacements for runtime authorization
checks. After new authorization, the runner rereads current state. Without
authorization, it still refuses execution; successful plan generation is not
approval. `--output` must equal `project.round_decision_file` in the config:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign.py" plan \
  --config "${SBDD_DATA_ROOT}/_staging/request-actual/R2-a1-config.json" \
  --decision "${SBDD_DATA_ROOT}/_staging/request-actual/R2-a1-decision.json" \
  --output "${SBDD_DATA_ROOT}/_staging/request-actual/R2-a1-plan.json"

"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/validate_config.py" \
  --input-file "${SBDD_DATA_ROOT}/_staging/request-actual/R2-a1-config.json"
```

After checking that the current config, agent decision, resolved plan, and user
authorization agree, launch the runner with that same config under the entry
skill's lifecycle rules. Replan and validate after every config change. Do not
forge hashes manually or reuse an old decision to execute unreviewed new
parameters. The plan also binds content fingerprints of existing `input_artifacts`
files and file inputs in `context`. Replacing a file at the same path also
invalidates the plan and requires review and binding again. Do not ask again
round by round within an already authorized scope.

Each run directory retains `round_decision.json`, `profile_snapshot.json`, and
`campaign_result.json`; state records the run and execution result. The profile
snapshot supports reproducibility and explains differences. It does not replace
each tool's own protocol/cohort identity. Do not combine percentiles from old and
new profiles. To compare rankings, first align the protocol and candidate set.

## Collect, Compare, and Assemble R3 Inputs

Use `scripts/campaign_evidence.py` on existing structured evidence. Replace
example paths with actual completed directories; outputs still belong under
`${SBDD_DATA_ROOT}`:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign_evidence.py" collect \
  --round-dir "${SBDD_DATA_ROOT}/project-actual/round1" \
  --round-dir "${SBDD_DATA_ROOT}/project-actual/round2" \
  --output "${SBDD_DATA_ROOT}/project-actual/candidate_archive.json"

"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign_evidence.py" compare \
  --before "${SBDD_DATA_ROOT}/project-actual/round1" \
  --after "${SBDD_DATA_ROOT}/project-actual/round2" --top-k 20 \
  --output "${SBDD_DATA_ROOT}/project-actual/round_comparison.json"

"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign_evidence.py" assemble \
  --archive "${SBDD_DATA_ROOT}/project-actual/candidate_archive.json" \
  --modifications "${SBDD_DATA_ROOT}/project-actual/modified/modification_report.json" \
  --output-dir "${SBDD_DATA_ROOT}/project-actual/inputs/R3-a1" --seed 42 --max-candidates 100
```

`collect` merges sources and retains REJECT audit records. `assemble` preserves
KEEP and MODIFY parent molecules as comparison baselines, then adds specified
derivatives. `--modifications` may be repeated or omitted. Assembly regenerates
3D structures from SMILES and writes `library.sdf`, `library_manifest.csv`, and
`candidate_pool.json`. Exceeding `--max-candidates` causes rejection; historical
KEEP molecules are not silently truncated. If the budget is insufficient,
explicitly choose a smaller target scope or obtain budget authorization, then
record the choice and rationale again. This step only prepares inputs; it does
not establish molecular improvement. Inject the resulting `library` and
`library_manifest` together, retaining a producer bundle from the same source,
then begin prefiltering at order 3. When newly generated/modified molecules are
included, explicitly set `context.candidate_source="molecule_generation"` to
use V3, or explicitly set `screening_mode="V3"` in every relevant gating step.
Do not infer V1 as a complete fixed user library merely because the library is
an external SDF. When only replicating an unchanged fixed library, retain its
original source and mode.

The following fragment wires a complete R3 config; the full step/tool matrix
must still be retained:

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "start_step": 3,
  "context": {"candidate_source": "molecule_generation"},
  "input_artifacts": {
    "library": "<SBDD_DATA_ROOT>/project-actual/inputs/R3-a1/library.sdf",
    "library_manifest": "<SBDD_DATA_ROOT>/project-actual/inputs/R3-a1/library_manifest.csv"
  }
}
```

Benchmark/decoy entries do not enter the candidate pool. If the round needs
control data, prepare them separately under existing tool contracts and retain
actual active/decoy labels. Do not label generated molecules as known actives.
`compare` reports rank changes only when protein/pocket/protocol/catalog,
scoring configuration, and candidate sets are comparable. Raw metric changes
and Top-K overlap across protocols are observations, not total-score improvement
on a common scale.

Pareto objectives must be explicitly specified by the user. Pass
`--objectives <JSON-file>` to `compare`, with a mapping from `metric_id` to
`min`/`max`. Without explicit objectives, do not select them automatically or
claim Pareto optimization. Do not fill missing metrics or incompatible protocols
with fabricated 0 values to obtain a supposedly better result.

## Final Confirmation in Stepwise Execution

When interrupt returns `awaiting_final_confirmation`, first present results
under the entry skill's Gate B and obtain approval through
`sbdd_confirm(gate="final", summary="...")`. Do not run the runner again or switch
to auto to remove the pause. Unlock the campaign attempt using the approval
record that actually occurred:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/campaign.py" confirm-final \
  --state "${SBDD_DATA_ROOT}/_staging/request-actual/campaign_state.json" \
  --run-id R1-a1 --authorization-ref 'The actual approval record for this final gate'
```

The command checks the completed checkpoint and results; it does not create
user authorization. Without final-gate approval, preserve the paused state and
do not invent a reference. This confirmation unlocks only the current round's
delivery; it does not expand the budget for later stages.

## Results and Compatibility

The coordinator outputs only `round_assessment`, while retaining per-candidate
classification, ranking, protocol gates, coverage, diagnosis, and nomination
eligibility. `decision="agent_review"`, `decision_outcome=null`, and
`goal_met=null` mean that the cross-round outcome awaits the agent's judgment,
not that execution failed. Whether results are deliverable and whether to
continue are separate questions.

Old complete configs without campaign fields retain single-round operation.
Running one does not automatically give it R1/R2/R3 meaning. Existing
success/iterate/max_rounds labels in old reports are historical records only;
they cannot confer new authorization. To bring old results into a new campaign,
first preserve their original evidence and specify their stage, applicable
protocol, and authorization. Do not merely renumber an old directory to pretend
independent replication is complete.

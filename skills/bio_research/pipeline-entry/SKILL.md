---
name: pipeline-entry
description: Run auditable SBDD experiments and research campaigns.
platforms: [linux]
---

# SBDD Pipeline Entry

**Entry point**: `pipeline-entry/run_pipeline.py` (a config-driven, single-round executor).

Before first use or migration to another machine, read the
[deployment configuration](references/deployment.md) and check the tool environments
and output root resolved by the current process. That configuration supersedes the
fixed host paths in reference documents; all scientific gates, authorization,
artifact contracts, and output-boundary requirements remain unchanged. Legacy
paths in the commands below are default-deployment examples and must be replaced
with the actual configured values before execution.

The Agent assembles the config, drives single-round execution, and determines
plans across rounds from the evidence. For a new experiment, first read the
[R1–R3 campaign protocol](references/campaign.md) and establish records for stages,
attempts, and authorization. Legacy single-round configs can still run independently.
Read the corresponding stage's SKILL.md before invoking any tool. The coordinator
assesses only the current round; gates, scoring, pose retries, and artifact reuse
remain the responsibility of deterministic code.

Before assembling a new experiment or migrating an old config, you must read the
[receptor preparation contract](references/preparation-contract.md). The structure
stage produces only a hydrogen-free geometric baseline. Final protonation, water
selection, and Meeko run in `docking.docking_prep` after prefiltering.
`context.chain` has no schema default: verify the target entity and set it
explicitly before chain-specific work. Some legacy downstream adapters still
fall back to chain A when this field is omitted; do not rely on that fallback.
The entry instructions, registry, adapter, and validator jointly enforce this constraint.

`molecule_generation` (Step 2, corresponding to Stage 04 in the delivery layout)
only acquires sources, standardizes structures, deduplicates exact structures,
and prepares a 3D library. `microstate_id` may be
empty when states have not been expanded. Multiple possible microstates are
recorded only for review and must not become an additional nomination veto.
`DEC_*` entries are never nominated; final candidates are still determined by
subsequent scientific evidence and the protocol. See the preparation contract
above for the implementation boundary of candidate-state expansion; reference-ligand
diagnostics must not be treated as candidate-library enumeration. The legacy state
scripts prepare chemical representations, rather than screen drug properties.
Property checks such as Lipinski, QED, SA, and PAINS remain the responsibility of
`property_prefilter` and `property_filter`, according to the selected config.

If ChEMBL acquisition fails, handle candidate sourcing through the
[molecule-generator Pocket2Mol fallback procedure](../molecule-generator/SKILL.md#chembl-unavailable-pocket2mol-recovery).
Preserve error records, the total candidate budget, and valid upstream artifacts;
do not repeatedly restart the same ChEMBL-only config. The Agent explicitly
updates and validates the plan; the runner does not automatically rewrite tool switches.

## Skipping a Molecule When All Selected Poses Explicitly Fail

The default is `pose_evaluation.params.selection_strategy=ranked_per_compound`:
preselect each molecule's poses by docking rank and evaluate all selected poses.
For ≤50 molecules, select at most 20 poses per molecule; for 51–1000, at most 5;
for >1000, only 1. Do not wait for the top-1 QC result before selecting additional
poses. Keep the remaining poses as supplementary SDF/manifest artifacts; unevaluated
poses do not enter the candidates. Tool execution failures and NA remain unresolved
and must not be falsely reported as passes or as failure of all poses. `max_poses`
can further limit the number actually evaluated; report partial coverage when the
budget is exhausted. Other explicitly configured strategies retain their existing
behavior. PoseCheck's `tools.posecheck.params.num_conformers` defaults to `10`,
and reference-ligand calibration uses the same setting. It is independent of the
number of docking output poses. When `mdanalysis` or `mmgbsa` is enabled, it may
change the verdict during postprocessing. The default `ranked_per_compound`
strategy has already evaluated every preselected pose, so those results can be
used directly for postprocessing.

Only when every pose preselected for a molecule in this run has been evaluated
and explicitly marked `fail` does
`pose_evaluation.selection.pose_exhausted_compounds` record `route=reject_molecule`,
`reason_code=all_poses_failed`, the number evaluated, and the specific failure
reasons. Exclude that molecule from `selected_poses.sdf`, retain its diagnostics,
and continue processing the other molecules. Do not trigger redocking of either
the individual molecule or the entire batch. This record means only that all
poses evaluated under this protocol were unsuitable; it does not prove that the
molecule is inactive. Tool execution failures, parsing failures, NA, and poses
left unevaluated because of insufficient budget do not count as scientific failure
of all poses.

Continuous mode proceeds downstream; interrupt mode pauses at the normal next
step. Resume through the ordinary checkpoint without consuming a historical
`pose_retry_state.json` or rewinding docking because of that state. Evaluation of
externally imported poses does not need docking to be enabled merely to handle
molecules whose poses all fail.

## Runtime Environments: Separate the Main and Tool Environments

- **Orchestration entry point**: run the validator and runner with the drugmate
  scientific main environment selected by deployment configuration through
  `SBDD_PYTHON`.
- **Stage scripts**: the runner selects an interpreter according to the step's
  `runtime`; if none is declared, it uses the entry-point interpreter.
- **Individual tools**: invoke tools according to their `runtime` and the
  registry's `TOOL_ENVIRONMENTS`. Do not assume that all tools reside in `drugmate`,
  or infer environment names from tool names. Autonomous orchestration follows
  the same rule.

Environments are managed by the registry and deployment environment variables.
Do not put Conda paths in experiment configs. Do not set
`project.python_executable` or a same-named top-level override to unify stage
interpreters, mix tool dependencies into `drugmate`, or silently replace tool
environments. The runtime layer launches across environments; manual activation
for each tool is unnecessary.

If environment preflight fails, read the current round's
`tool_runtime_report.json` and
[runtime environment troubleshooting](references/runtime-environments.md) to
identify the affected tool environment.

## Experiment Lifecycle and Shared Resources

Before launch, read [cancellation, resource waits, and failure cleanup](references/execution-lifecycle.md).
Start the experiment with the local
`terminal(background=true, cancel_on_interrupt=true, notify_on_complete=true, pty=false)`.
Run the runner with `exec` and retain the returned `session_id`. Do not detach it
from management with `nohup`, shell `&`, or a fixed `timeout 8h`. An Agent
interruption cancels this task's experiment process tree. To stop it explicitly,
call `process(action="kill", session_id=...)`. Read status with
`process(action="poll", session_id=...)`, or wait for a bounded interval with
`process(action="wait", session_id=..., timeout=30)`. Every wait must explicitly
set a timeout of 1–30 seconds; do not omit it. Unbounded shell PID-wait loops such
as `while kill -0 PID; do sleep 30; done` are prohibited. A polling timeout means
only that execution has not finished; it is not grounds for killing the experiment.

Use the same managed background parameters for input conversion. Read completed
and total counts, skipped counts, and the molecules currently running from this
round's `inputs/normalization/normalization_progress.json` before updating the
user. The final authority is `normalization_report.json`. SMILES conversion uses
at most 8 concurrent workers by default, consistently uses `--gen3d fast`, and
defaults to 60 seconds per molecule. Omitting `library.conversion_timeout_seconds`
or setting it to `null` both select 60 seconds; a positive integer can override
that limit. The script marks timed-out molecules as `conversion_timeout` and
continues with the others. The Agent must not substitute a wait timeout for the
per-molecule timeout or stop the entire library because one `process` wait timed
out. A fast-mode capability probe runs before conversion. If it detects the
Open Babel 3.1.1 coordinate-writeback defect, the script compiles and caches a
compatibility executable under `inputs/normalization/runtime/` that fixes only
writeback, retaining the fast algorithm. Only this compatibility path requires a
C++17 compiler, matching Open Babel development headers/libraries, and Eigen 3.
It does not modify the global Conda environment. Probe or dependency errors fail
the entire normalization run and must not be skipped as bad molecules. Inspect
`library.conversion_execution.command_backend` and the runtime artifacts to
confirm the backend actually used.

By default, at most 3 experiments run concurrently on the same host using the
same resource directory. If resources are busy, return
`status=resource_unavailable`, `retryable=true` immediately and tell the user to
retry later. `runtime_resources.wait_timeout_seconds` can explicitly set a
separate queueing budget; queueing does not count against step execution limits.
While execution runs in the background, read the manifest/resource_status rather
than imposing an independent 8-hour wall-clock limit.

## Output Root and Config Storage (Mandatory)

From the first write onward, all experiment outputs on disk must reside under
the deployment-configured `SBDD_DATA_ROOT`.
This includes each round's config drafts, submitted configs, input-normalization
artifacts, preflight results, temporary files, caches, tool outputs, and logs.
Do not write them to the system temporary directory, the source checkout, or any other location and move
them afterward. Existing external inputs may be accessed read-only.

- Every config explicitly sets `project.data_root` to the absolute `data_root`
  returned by the deployment check.
- Before the first write, call the runner's `--init-project-dir` with minimal
  JSON specifying `project.data_root`, `project.project_slug`, and
  `project.round_number`. Also pass `--request-id '<stable ID for this request>'`:
  prefer the platform request/message ID. If none is available, generate a token
  once for this task, retain it in the session, and reuse it on initialization
  retries instead of generating a new token each time. The same data_root,
  request-id, and round_number bind to the same directory under a lock; different
  requests are not merged. This operation only allocates the current round's
  directory. It does not run steps or require a scientific config that has not yet
  been assembled. Use the known target or a short request name for `project_slug`;
  the name is not evidence of protein identity.
- The returned `round_dir` equals the compatibility field `project_dir` and is
  the execution directory. `project_root` is the parent containing the rounds
  and must not be passed to execution commands. The initialization record stores
  the request-id, creation time, and owner session. An initialization response
  only describes the directory binding; it does not establish whether the
  experiment is running or complete. When continuing the same task after an
  interruption, recover round_dir from the existing initialization record first.
  A new session explicitly continuing an old task also reuses that directory
  rather than allocating another. Do not automatically merge different tasks
  even if they use the same PDB.
- Save config drafts, sampled inputs, preflight logs, and similar files to the
  returned `inputs_dir`. Save input-conversion artifacts to the returned
  `normalization_dir` (that is, `project_dir/inputs/normalization/`). After the
  complete config has been assembled, validated, and authorized, execute with
  `--project-dir` pointing to the same directory. Do not write conversion
  artifacts directly into `published/`, which is rebuilt. Preallocating a
  directory or saving drafts does not constitute plan approval; Gate A below
  still applies.
- After acquiring the round-directory lock, the runner saves the submitted
  config to the returned `project_dir/pipeline_config.json`. It writes
  `project_dir/resolved_run_plan.json` after validation and checkpoint checks
  succeed. Here, `project_dir` is `roundN/`, not its parent project root.
- For subsequent attempts, copy the selected upstream `pipeline_config.json` to
  `configs/roundN.json` under the same project root, then modify, validate, and
  confirm it. The runner archives it again as the new
  `roundN/pipeline_config.json`. Never overwrite the previous round's config.
- For both initial execution and resume, use CLI `--project-dir` to specify the
  current round's directory returned by the runner. Do not put this override in
  a config that will later be copied to the next round. If an existing config
  contains `project.project_dir` or a top-level `project_dir`, remove it when
  creating the next round's config; otherwise it will still point to the old
  round.
- Check every output, working, temporary, cache, and log path. All must remain
  within the configured root after resolving relative paths and symbolic links.
  Do not move or delete legacy `_staging` inputs still referenced by experiment
  configs. Legacy initialization calls without a request-id still allocate a
  new directory and refuse to occupy an existing one; retries with the same
  request-id only return the original binding and do not overwrite inputs.
  Naming conflicts or explicit-directory conflicts produce errors. Before the
  first execution, the runner rejects a parent directory containing initialized
  rounds when it is passed as round_dir. Historical layouts with existing
  checkpoints remain compatible.

This is a path constraint the Agent must enforce, **not a global path sandbox
already implemented by the runner**. The code still permits other `data_root`
values and explicit directories. Set the intended output root explicitly in
each config and verify it through the deployment check.

## User Confirmation and Execution Modes (Mandatory Contract)

`project.execution_mode` controls whether execution pauses after each step
**once the current round's plan has been authorized**:

- `interrupt` (default): after each new step is persisted, the runner writes
  `pipeline_checkpoint.json` and returns control to the Agent. The Agent must
  request user confirmation before resuming.
- `auto`: automatically execute all remaining steps in the round from the
  current checkpoint.

Both modes must pass Gate A. `auto` only disables per-step confirmation; campaign
A0–A3 determine the scientific stage through which automatic progression is
authorized. These are independent controls. A0 pauses after each round's delivery;
A1 defaults to R2. In this implementation, A2/A3 extend at most through R3 and
must never be used as authorization to start R4–R6 or MD/FEP.

### Gate A: Plan Validation and Authorization Scope (All Modes)

1. Complete query → fact bundle → policy selection → config assembly in memory
   or in a config draft within the output directory above. Run
   `validate_config.py` until `valid=true`. An invalid plan must not execute.
   Single-round confirmation prompts in the validation report do not replace
   campaign authorization checks.
2. Do not present config JSON as the confirmation content. First provide a
   natural-language summary the user can understand: the round's objective and
   inputs; the explicit site/search strategy and its rationale; the ordered
   steps; each step's tools and their purposes; key time/cost parameters;
   expected artifacts and main risks; and the current `execution_mode`. A config
   path or JSON may be offered only as optional advanced detail.
3. First check existing user authorization. If the campaign plan passes the
   [authorization and budget checks](references/campaign.md), and its target,
   stage, tools, parameters, cost, and execution mode are all within explicit
   preauthorization, record that basis and execute without asking again whether
   to continue this round. If existing authorization does not cover the plan,
   call `sbdd_confirm(gate="plan", summary="...")` on its own and show the
   requested campaign scope. The tool provides fixed options for confirming and
   executing, adjusting the plan, or cancelling the experiment. Use the exact
   labels supplied by the tool; these English descriptions do not replace its
   choice strings. It must be the only tool call in that assistant
   tool-call batch; do not start the runner in the same batch.
4. Only `approved=true` constitutes approval. It can result only from the user
   explicitly selecting the approval option or from the tool applying its
   timeout default after 120 seconds. The latter explicitly returns
   `decision_status="timeout_defaulted"`, `default_applied=true`, and must not be
   described as an active user click. A request to adjust/stop, free text,
   cancellation, an empty response, failed delivery, callback error, or
   oneshot/headless execution **does not count as approval**. Do not start the
   experiment in those cases.
5. If the user requests changes, modify the proposed config, repeat the dry-run,
   and recheck authorization scope. Request additional authorization only when
   the stage, tool/parameter scope, or budget exceeds what was authorized;
   existing authorization remains valid. A timeout default can only apply the
   original tool's single-round approval contract and cannot substitute for
   explicit authorization of a budget across rounds. Remain paused when there
   is no reply at A0 or at an authorization boundary.

The exact gate labels are defined in
[`_gate_spec` in the confirmation tool](../../../tools/sbdd_clarify_tool.py).
Consult that definition if the tool is unavailable; do not invent or translate
the affirmative label required by the fallback below.

For Gate A or Gate B, if the current execution context does not expose
`sbdd_confirm`, or if the tool returns `decision_status="unavailable"` because
there is no interactive callback, present the same natural-language summary in
an ordinary message and end the current Agent turn. Wait for the user to reply
with the relevant option's **complete affirmative label**. This fallback has no
background timer; do not simulate timeout approval yourself. Do not import/call
`sbdd_clarify_tool.py` inside the pipeline subprocess, which has no platform callback.

### Gate B: Per-Step Confirmation (`interrupt` Only)

Each computational step's CLI stdout and `step*/result.json` use a compact
`result_schema="sbdd-step-result/v1"` control object. It retains `success`,
errors, `files`, a summary, and metadata needed for orchestration; it does not
embed scores and evidence for the entire molecule batch. Read the full contents
from report/ranking paths in `files`. Do not treat stdout or `result.json` as a
complete scientific report or print an entire report to the terminal. The
per-molecule review list in `round_assessment` remains in the coordinator report;
the control object contains only its count. The runner accepts only complete
JSON that satisfies the return contract. Truncated or malformed output is a
failure; never assemble a successful result from internal fragments.

The structured return from `run_pipeline.py` is a control contract. Do not check
only `success`:

- `status="awaiting_step_confirmation"`: the current step has been persisted,
  `pipeline_complete=false`, and `next_step` identifies the next step.
- `status="awaiting_final_confirmation"`: the round's last enabled step has
  been processed and `pipeline_complete=true`. Confirmation is still required
  before the round summary or next-round decision.
- `status="complete"`: in `auto` mode, the round's computation has finished
  automatically. `interrupt` also completes directly when there are no enabled
  steps, because there is no per-step gate. This means only execution complete,
  **not that the scientific objective has been achieved or molecules can be nominated**.
- **Only the two `awaiting_*` statuses above may enter Gate B**.
  `config_invalid`, `runtime_preflight_failed`, `project_busy`,
  `resource_unavailable`, `interrupted`, `failed`, `checkpoint_invalid`,
  `publication_failed`, and any unknown status fail closed. Do not resume them
  as confirmation points; report and handle the error first.

Likewise, `success=true` and `pipeline_complete=true` describe only successful
runner/step execution. When coordination has run, also read the consistent
scientific-outcome fields in the runner return, `pipeline_checkpoint.json`,
`manifest.json`, and `scientific_outcome.json`:

- `decision="agent_review"`, `decision_outcome=null`, `goal_met=null`: the
  coordinator no longer decides whether to continue, terminate, or declare the
  campaign objective achieved. Read the current round's facts in `round_assessment`.
- `valid_for_nomination`: indicates only whether the round's deterministic
  evidence permits candidate nomination, not whether the entire campaign has
  ended. The protocol and candidate eligibility do not change because the Agent
  chooses to continue or finish.
- `round_assessment`: molecule/pose classifications, protocol status, coverage,
  diagnostics, and nomination eligibility. When the protocol fails, follow
  `rollback_instructions` to preserve invalid results and lineage non-destructively.
- If coordination has not run, scientific assessment is `not_evaluated`.
  Historical outcomes from old reports may be shown, but legacy
  `iterate/success/max_rounds` must not be treated as control instructions for a
  new campaign.

Whenever `requires_confirmation=true`:

1. Read the returned `result_json` (if present), `manifest_json`, and key
   reports. Produce a natural-language summary instead of merely pasting paths
   or JSON: what the step did, which tools it called, key results/metrics,
   artifacts generated, warnings/failures/applicability limitations, and what
   `next_step` will do if execution continues.
2. For `awaiting_step_confirmation`, call
   `sbdd_confirm(gate="step", summary="...")` on its own. If
   `confirmation_reason="checkpoint_rewind"`, instead call
   `sbdd_confirm(gate="rewind", rewind_target="<rewind_target>", summary="...")`
   and explain `resume_rewinds` in the summary. The tool generates fixed options;
   resume only when it returns `approved=true`.
3. For `awaiting_final_confirmation`, call
   `sbdd_confirm(gate="final", summary="...")` on its own. Enter the round
   summary/iteration only after it returns `approved=true`.
4. **Only approval of `awaiting_step_confirmation` triggers another runner
   call**. Use the same config and the returned `project_dir`, for example by
   appending `--project-dir '<project_dir>'`. The persistent checkpoint resumes
   from the next unprocessed step, including PocketMiner, which normally
   prohibits cache reuse. Do not create a new project directory or infer the
   current position from repeated cache hits.
5. After `awaiting_final_confirmation` is approved, **do not call the runner
   again**. For a campaign, use the
   [confirm-final command](references/campaign.md#final-confirmation-in-stepwise-execution) to record
   this actual approval and unlock the completed attempt. In the legacy
   single-round workflow, the current conversation layer tracks approval.
   Then read coordinator/reporting results and proceed to the round summary or
   next-round decision.
6. For the tool's config-adjustment option, repeat the dry-run and return to Gate A
   to check authorization for the revised plan. A campaign's scientific config
   is frozen. If an unfinished attempt requires changes to scientific inputs
   or parameters, follow the
   [campaign protocol](references/campaign.md#stages-and-attempts) to preserve and end
   the old attempt and create a new attempt at the same stage. In legacy
   single-round mode, the runner can still rewind to the earliest affected step
   based on fingerprints; this does not permit overwriting a frozen campaign plan.
7. Interpret approval, timeout, and failure responses as in Gate A item 4.
   Without approval, do not call the runner again or enter the summary/next round.

Only at `awaiting_step_confirmation`, and only if the user explicitly asks to
“automatically finish the remaining steps,” may `execution_mode` change to
`auto`. After the change, repeat the dry-run and pass Gate A in full.
Do not switch modes at `awaiting_final_confirmation`; complete the final gate.
The Agent must not switch modes on its own. After approval, keep the same
`project_dir`; the runner preserves the confirmed prefix and automatically runs
the remaining steps. This UI/orchestration switch does not change scientific
step fingerprints.

## Query → Config (Required Before the First Round)

When the user supplies a fixed docking box, enable
`pocket_analysis.params.site_selection_strategy=provided_box` in the same round,
retaining `input_artifacts.docking_box` and
`docking.params.search_space_source=docking_box`. Step 1 analyzes the receptor
region corresponding to that box and produces reports needed by coordination;
it does not overwrite the original box. Explicitly select detectors and
interaction tools as needed. A separate auxiliary project is unnecessary for
supplementing the reports.

After receiving a user query, handle the responsibilities in the following
order. Do not duplicate individual tool-selection policies in this entry file:

1. Preallocate the round directory as described above. If the user supplies a
   local receptor, fixed library, or docking box, read
   `../input-normalization/SKILL.md` and the
   [input wiring instructions](references/input-normalization.md). Invoke
   `../input-normalization/run.py` in a managed background process and write its
   outputs to the returned `normalization_dir`. For SMI, run Open Babel
   `-h --gen3d fast` with bounded per-molecule concurrency, defaulting to 60
   seconds per molecule. Mark timeouts as `conversion_timeout` and skip those
   molecules; inspect SDF input per record. Read `normalization_progress.json`
   during execution. Record and skip a molecule whose conversion or coordinate
   check fails, and continue with the others; do not resample to replenish the
   count. Read command diagnostics, `library.skipped_records`, input/retained/
   skipped counts, and `library.geometry_qc` from `normalization_report.json`.
   Merge `config_patch.json` only if `success=true`. No usable molecules or an
   infrastructure error still constitutes failure. The final report retains
   preprocessing failure reasons and empty scores for skipped molecules.
   Preserve raw conversion outputs separately. Output usable records in source
   order without rewriting coordinates. The conversion layer does not fill in
   IDs, enumerate chemical states, or introduce conversion approval. Receptor
   format conversion preserves coordinates and does not add 3D-generation
   parameters; handle boxes according to their respective formats.
2. Assemble verifiable facts, the user's original intent, uncertainties, and
   external artifacts from the query, original files, and normalization report.
   Do not infer whether receptor chemical preparation is complete from a file
   extension, PDB ID, or artifact name. Connect external libraries through
   `library`; after conversion, proceed to property prefiltering or docking
   preparation according to the plan. Use “Unified Memory Consultation” below
   to determine the explicit pocket-site/memory-fusion strategy and check the
   applicability of generation and docking memory. Complete selection once
   valid inputs are ready and record what was actually adopted or why it was
   overridden. Do not leave unresolved memory-input dependencies for downstream
   execution to guess.
3. Read the generated `TOOL_SELECTION_POLICY.md` and select every step, tool,
   required parameter, `pocket_analysis.params.site_selection_strategy`, and
   `docking.params.search_space_source`. Read tool IDs, parameter schemas, and
   strategy choices only from `tool_catalog.json` / `TOOL_CATALOG.md`.
4. Assemble the config according to `QUERY_TO_CONFIG.md`, satisfying the
   complete switch matrix and strategy requirements below.
5. **Pre-run config validation (dry-run)**: run `validate_config.py` using the
   Quick Start command. Exit 0 means the static config is valid. On exit 1,
   correct the `errors` and rerun. Check that `resolved_plan` / `plan_summary`
   match the original selections. The dry-run does not validate actual upstream
   files, lineage, runtime, or final-publication completeness. `valid=true` does
   not guarantee successful runner execution. See
   [config validation instructions](QUERY_TO_CONFIG.md#pre-run-config-validation) for
   static assumptions across rounds.
6. For a new experiment, save the R1 plan and explicit authorization according
   to the [campaign protocol](references/campaign.md). Present input-normalization
   results, strategy, the full plan, and authorization scope under Gate A. After
   approval, pass the same config to `run_pipeline.py`.

## Unified Memory Consultation

Before assembling the main experiment, read the
[unified memory protocol](references/benchmark-memory.md). In the current round's
`inputs_dir/memory_consultation.md`, record applicability, strategy, query status,
version, recommendations, and the actual decision for pocket, generation, and
docking together. Record whether pocket memory was used, unavailable, disabled,
or bypassed by an explicit site, with the actual reason and decision provenance.
For generation/docking, record missing inputs as pending preparation, or explain
why the task is inapplicable. An unperformed query must not be reported as finding
no neighbors; do not defer the tool matrix decision until after the runner starts.

- [Docking](references/docking-memory.md): retain
  `context.docking_memory_mode` choices `ligand_protein` (default),
  `ligand_only`, `legacy`, and `off`. `select-docking --task redocking` defaults
  to Top-1 RMSD <2 Å and GNINA. `--task pose_generation` uses the rDock
  baseline with a published Top-1 PoseBusters objective and
  shared prior strength 2 / switching margin 0.10. Record its `task_policy`
  assumptions; this policy does not establish measured performance for a new query.
  Memory-off retains the selected task's baseline. Every plan records one task,
  objective and primary engine; see the
  [task boundaries](references/docking-memory.md#task-boundaries).
- [Molecule generation/design](references/generation-memory.md): recommendations
  are independent for five metrics, with `property_focused` as the default.
  Priority selects the objective; every objective retains prior strength 10 and
  switch margin 0.03. Do not
  force a query when a complete candidate library already exists or no valid
  conditioning ligand is available. Do not treat generation memory as a
  molecule-modifier reaction strategy.
- [Pocket generation/selection](references/pocket-memory.md): for an unknown
  site, use `detector_consensus` in the main Step 1. Prepare compatible protein
  encoding and a pocket index, then call `plan-pocket-tools` before finalizing
  the matrix. It retrieves five protein groups and records their best-detector
  multiset. Enable exactly its `execution_tools`: unique recommended detectors
  plus DeepPocket. All five detectors, including AutoSite, are registered.
  Normalize native Top-3 support to point5, estimate historical DCA reliability,
  apply all six P2Rank anchor guards, and fuse all matching pockets by
  reliability/native rank/family.
  Keep the selected anchor's point5 residues. Use fixed DeepPocket#1 rank fusion
  only when memory is unavailable or explicitly disabled; record the reason.
  Configure registered `pocket_memory_*` step parameters and prepare query assets
  before the run. The fallback plan enables DeepPocket, fpocket, P2Rank and
  pyKVFinder. Follow the pocket reference's selection and normalization rules.
  A fixed box/reviewed site still takes priority.

Invoke generation/docking selectors through `scripts/benchmark_memory.py`; they
do not automatically modify the config. Pocket fusion runs within the main
Step 1. If predictions from every planned detector already exist, use the same CLI's
`normalize-pocket-candidates` / `select-pocket` for independent review. Do not
preselect and then switch to `provided_box` to recalculate residues. The Agent
checks runtime, inputs, and budget; applies the complete switch matrix,
parameters, and site/search strategy; revalidates; and checks authorization
under Gate A. Explicit user choices take priority. The original benchmark stays
frozen. In ordinary memory-enabled research, use the
[research-memory protocol](references/research-memory.md) to append evaluated
results after each run and recall prior observations before the next plan.
Record the scope as research or frozen in consultation; benchmark evaluation and
memory-off ablations bypass this store. `context.docking_memory_mode=off`
disables docking memory only. Pocket memory is controlled separately by
`steps.pocket_analysis.params.pocket_memory_mode`. For ablation rules, see the
unified protocol. Do not add unregistered memory fields.

## Configuration Model (Three Layers)

Every registered step and every registered tool must explicitly declare
`enabled=true/false`. Registry defaults, scientific facts, and the runner do not
select switches on the Agent's behalf. Set every tool under an explicitly
disabled step to false. Reused producers before `start_step` retain their
original matrix; see the reuse boundaries below. Parameters without overrides
may be omitted and are deep-merged with schema defaults.

```
step   (whether to execute this step)    steps.<step_id>.enabled
tool   (which tools this step uses)     steps.<step_id>.tools.<tool_id>.enabled
param  (tool/step parameters)           steps.<step_id>.tools.<tool_id>.params.*  /  steps.<step_id>.params.*
```

### Tool Selection and Execution Boundaries

Two modes that change the scientific route must be explicitly set in step params:

- `pocket_analysis.params.site_selection_strategy`: `ligand_site` / `known_site` / `detector_consensus` / `provided_box`
- `docking.params.search_space_source`: `reference_ligand` / `docking_box` / `whole_receptor` / `precomputed_grid`

`context.scientific_facts.has_metal/is_covalent/protonation_sensitive` are
scientific inputs and validator evidence, not selectors. The orchestrator
automatically wires file paths between steps according to artifact dependencies.

### Top-Level Structure

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "project": {
    "project_slug": "cdk2-4kd1",
    "data_root": "<SBDD_DATA_ROOT>",
    "round_number": 1,
    "upstream_project_dir": null,
    "execution_mode": "interrupt"
  },
  "context": {
    "pdb_id": "4KD1",
    "target_name": "CDK2",
    "chembl_target_id": "CHEMBL301",
    "chain": "A",
    "seed": 42,
    "docking_memory_mode": "ligand_protein",
    "candidate_source": "sampled_library",
    "scientific_facts": {
      "has_metal": false,
      "is_covalent": false,
      "protonation_sensitive": false
    },
    "known_site_residues": ["Ile10","Val18","Ala31","Lys33","Glu81","Leu83"],
    "known_site_provenance": {
      "source": "reviewed kinase-site annotation",
      "method": "independent mechanism/SAR review",
      "independent_of_reference_structure": true,
      "review_status": "approved"
    },
    "known_required_interactions": [{
      "residue": "Leu83", "chain": "A", "interaction_type": "H-bond",
      "protein_atom_names": ["N"], "protein_role": "donor", "ligand_role": "acceptor",
      "rationale": "reviewed ATP-site hinge anchor",
      "source": "PDB:4KD1 plus reviewed CDK2 mechanism",
      "binding_mode_scope": "ATP-competitive CDK2 ligands",
      "review_status": "approved",
      "logic_group": "hinge_anchor", "logic_operator": "any_of"
    }]
  },
  "step_timeouts": {},
  "steps": {
    "<every_registered_step>": {
      "enabled": false,
      "params": {},
      "tools": {
        "<every_registered_tool_in_this_step>": {
          "enabled": false,
          "params": {}
        }
      }
    }
  }
}
```

This JSON illustrates only the config hierarchy. It is neither a submittable
config nor a tool-selection template. Expand every placeholder and populate it
according to actual facts and policy decisions. Do not copy the example's site
or interactions without establishing their applicability.

`context` collects identity/shared parameters such as
`chain`/`seed`/`chembl_target_id`/`known_*` instead of repeating them per step.
If an externally injected library was randomly sampled outside the pipeline,
explicitly set `candidate_source="sampled_library"` to enable automatic V2
selection. Use `user_supplied_full` (V1) for a small user library retained in
full, and `molecule_generation` (V3) for a library generated by this pipeline.
If gating steps explicitly specify `screening_mode`, that explicit mode takes priority.

`step_timeouts.<step_id>` is the step's hard execution limit after resources are
acquired, excluding queueing. If unset, the orchestrator sums tool timeouts
within the step, adds a 10% + 1800-second cleanup margin, and applies a
conservative minimum step budget suitable for large molecule batches. Parallel
tools retain a conservative budget. An explicit limit must not be below the
automatic budget or exceed 30 days; split exceptionally large tasks into shards.

## Quick Start

Execute the following commands in stages, not as an unattended block. Gates A/B
are handled in the Agent conversation layer. Only directory-naming information
is needed before conversion. Assemble the full config after conversion succeeds,
explicitly setting `project.data_root` to the deployment-configured root,
`project.round_number=1`, and the actual `project_slug`.

```bash
: "${BIO_SKILLS_ROOT:?First configure the absolute skill path for this checkout as described in the deployment guide}"
: "${SBDD_PYTHON:?First specify the absolute Python path of the configured scientific main environment}"

# Preallocate before the first write; replace the example project_slug with a short name for this request.
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/run_pipeline.py" \
  --init-project-dir \
  --request-id '<stable ID for this request; reuse it on initialization retries>' \
  --input "$(python3 -c 'import json, os; print(json.dumps({"project": {"data_root": os.environ["SBDD_DATA_ROOT"], "project_slug": "mapk1-screen", "round_number": 1}}))')"

# Read the preceding JSON response and use its actual absolute paths in subsequent terminal calls.
SBDD_PROJECT_DIR='<project_dir returned above>'
SBDD_REQUEST_DIR="${SBDD_PROJECT_DIR}/inputs"
CONFIG_FILE="${SBDD_REQUEST_DIR}/round1.json"

# If input conversion is needed, save its request JSON in inputs/ and set --output-dir to the returned normalization_dir.
# After successful conversion, merge config_patch.json and write the complete config to CONFIG_FILE.
# Save preflight logs and other request artifacts in SBDD_REQUEST_DIR as well.
# 0. Complete unified memory consultation and apply the relevant choices for all three domains; then dry-run query→config.
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/validate_config.py" --input-file "${CONFIG_FILE}"

# 0.5 All modes: present plan_summary; check campaign preauthorization and, if it does not cover the plan, confirm with sbdd_confirm(gate="plan").

# 1. After approval, execute the same config file that was just validated/confirmed; default interrupt advances only one new step per call.
# Submit this command to terminal(background=true, cancel_on_interrupt=true, notify_on_complete=true, pty=false).
exec "${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/run_pipeline.py" \
  --input-file "${CONFIG_FILE}" \
  --project-dir "${SBDD_PROJECT_DIR}"

# 2. In interrupt mode, on awaiting_step_confirmation, first summarize in natural language and call sbdd_confirm.
#    Only approved=true (user selection or timeout_defaulted) permits resume using the same returned directory:
exec "${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/run_pipeline.py" \
  --input-file "${CONFIG_FILE}" \
  --project-dir '<project_dir returned by the previous call>'
```

Shell variables are not guaranteed to persist across terminal calls. In later
calls, use the actual absolute paths generated the first time; do not initialize
another directory. If stdout/stderr from a call need to be saved, write them to
the same request directory or the current round's already-created `logs/`.
After the runner returns, its `pipeline_config` field points to the current
round's snapshot. The originally submitted config also remains within the
deployment-configured root.

## Pipeline Steps (Named Step IDs)

| order | step_id | skill | Main output artifacts |
|---:|---|---|---|
| 0 | `structure_prep` | structure-prep | prepared_receptor_noH, prep_manifest, source_complex_mmcif, reference_ligand, water_inventory |
| 1 | `pocket_analysis` | pocket-analysis | pocket_definition, docking_box, interaction_observations, interaction_blueprint, pocket_features, pocket_center; also deeppocket_predictions when explicitly enabled |
| 2 | `molecule_generation` | molecule-generator | library, library_manifest |
| 3 | `property_prefilter` | drug-property-filter | prefiltered_library |
| 4 | `docking` | docking | docked_poses, docking_scores, docking_report |
| 5 | `pose_evaluation` | structure-evaluator | selected_poses, evaluation_report; also rtmscore_manifest, candidate_topologies, short_md_manifest, trajectory_stability_report, mmgbsa_manifest, openfe_rbfe_manifest, ties_rbfe_manifest when explicitly enabled |
| 6 | `property_filter` | drug-property-filter | property_report, synthesis_report |
| 7 | `coordination` | reflective-coordinator | coordinator_report, round_report.md |
| 8 | `reporting` | report-generator | report.html, report_data.json, and resource/artifact manifests; final publication requires this step to succeed; it does not drive scientific iteration decisions |

> Disk directories retain legacy numbering for compatibility with the
> report/aggregator: registry order 5, `pose_evaluation`, writes to `step4`;
> order 6, `property_filter`, writes to `step5`; and `property_prefilter` writes
> to `step2_5`. The step_id is the canonical name in configs and manifests.

Full publication (`publication_mode=pipeline_final`) requires joint completion
of `coordination` and `reporting`. If a plan includes either one, it must include
both. At the final-publication check, a missing stage, a missing successful
result from either stage, or failed artifact validation causes the runner to
return `status="publication_failed"`, `success=false`, `pipeline_complete=false`.
Although the registry still marks reporting as `best_effort=true`, that flag
affects only exception handling within the step and does not bypass the
final-publication check. A run that executes only earlier experimental stages
with both coordination and reporting disabled can complete a partial pipeline
and produce a partial-results view with `publication_mode=pipeline_partial`.

The runner updates `published/` at stage start/end, per-step pauses, and
catchable failures. A partial view contains current stage status, an index of
existing files, log references, and `07_final_report/final_report.md`,
`final_report.html`, `final_report_data.json`, and `publication_status.json`.
A report is generated even if docking/Pose-QC or an earlier stage fails. It
lists the reason execution stopped, existing files, and unfinished stages;
unfinished metrics do not produce aggregate scores or rankings.
`published_dir` / `user_output_dir` indicate only that a publication view is
available, **not that the workflow succeeded or delivery is complete**. Always
read `publication_mode`, execution status, and the scientific outcome together.
Failures remain `success=false`; partial publication does not bypass Gate A/B
or the strict checks for full publication. If partial publication itself fails,
the runner returns `publication_error` while retaining the original experiment status.

## Iteration Logic

After completing this round's computation and any applicable final gate, read
`step6/coordinator_report.json`, reconcile the round report and artifact lineage,
and use the [iteration decision guide](ITERATIVE_GUIDE.md) to decide whether to
repair within the same stage, advance to the next stage, or deliver results and
pause. The coordinator's `round_assessment` is evidence, not an instruction
across rounds. The Agent records the decision, rationale, and input changes in
the campaign plan.

R1 performs rapid screening. R2 independently reproduces results/improves
conditions. R3 combines KEEP results from R1–R2 and optimizes MODIFY results.
The scientific `campaign_stage` is separate from `campaign_attempt` and the
physical directory number `round_number`: repairing R1 remains R1 and does not
become R2 just because the directory number increases. Stage budgets and
authorization are governed by the [campaign protocol](references/campaign.md).
Legacy `max_rounds`, a 30% pass rate, or a KEEP count of 3 no longer control the
campaign. After reaching R3 or an authorization boundary, deliver the existing
results and pause.

### Directory Structure

Each attempt uses a separate subdirectory. `run_pipeline.py --init-project-dir`
automatically names a new project
`{YYYYMMDD-HHMMSS}-{run_id}-{project_slug}`. The random `run_id` prevents
experiments launched in the same second from sharing a checkpoint. The Agent
sets only `project.project_slug`, rather than constructing a dated `project_dir`
manually. Campaign state and plan locations are explicitly referenced by
`project.campaign_state_file` and `project.round_decision_file`. They may be
saved in this request's `inputs/` directory. Preserve existing `_staging` references.

```
{data_root}/{YYYYMMDD-HHMMSS}-{run_id}-{project_slug}/
├── configs/         ← Submitted configs for subsequent rounds: roundN.json (created by the Agent)
├── round1/          ← First attempt; see campaign records for stage/attempt
│   ├── inputs/      ← Config drafts, sampled inputs, preflight logs, and request artifacts
│   │   └── normalization/ ← library/, receptor/, normalization_progress.json, normalization_report.json, config_patch.json
│   ├── project_initialization.json
│   ├── step0/ ... step7/ (also step2_5/ for property_prefilter)
│   ├── published/   ← Partial view updated with progress, or full view after completion
│   ├── pipeline_config.json
│   ├── resolved_run_plan.json
│   ├── round_decision.json
│   ├── profile_snapshot.json
│   ├── campaign_result.json
│   ├── manifest.json
│   └── logs/
├── round2/          ← Second attempt; not automatically R2
└── final_report.md  ← Final report across rounds, written after the Agent's loop ends
```

The runner continues to use `project_dir` / `working_dir` as checkpoint-resume
paths. After partial or full publication passes validation, it additionally
returns `published_dir` / `user_output_dir`. Users and UIs should prefer the
latter paths, using `publication_mode` to distinguish partial results from
complete delivery. `pipeline_final` requires successful coordination and
reporting, completion of the entire round, and successful final-publication
validation. For publication invalidation and archiving, see
[artifact reuse and recovery](references/artifact-reuse.md).

The fixed top-level order within `published/` is `00_inputs/`,
`01_structure_prep/`, `02_pocket_analysis/`, `03_molecule_generation/`,
`04_sm_property/`, `05_docking/`, `06_post_docking/`, and `07_final_report/`.
These are scientific stages from the user's perspective. Each directory's
manifest identifies its actual runner `step*` source; this presentation does
not permit changes to the internal scoring semantics or column names for
`step4/step5/step6`.

### Reuse Across Rounds and Execution Scope

Copy and store configs according to “Output Root and Config Storage” above.
Do not set `upstream_project_dir` for the first round. For later attempts, set
the actual upstream directory for the selected protocol and an incremented
`round_number`, reusing the same project root. Selecting the R1 protocol does
not require the R2 directory to be the upstream source. Campaign decisions must
still reference evidence from the immediately preceding attempt and the selected
historical conditions.

- `start_step` / `stop_step` limit execution by registry order (0–8); they do not
  replace the complete switch matrix.
- To reuse a disabled upstream producer, declare `start_step` and retain that
  producer's original step/tool/parameter matrix before the boundary. Scope
  filtering disables those steps during actual execution.
- Explicitly disable tools after `stop_step`. If upstream inputs or parameters
  change, rerun from the earliest affected stage. Never force reuse of old
  results that fail validation.
- Before injecting existing artifacts or investigating rejected cache reuse
  or a checkpoint rewind, read
  [artifact reuse and recovery](references/artifact-reuse.md) to avoid combining
  inconsistent upstream and downstream evidence.

For later-round config examples and mapping diagnostics to a rerun scope, see
the [iteration decision guide](ITERATIVE_GUIDE.md).

## Quality Gates

The AUROC/EF1% values below are default thresholds; the round's confirmed
config and structured gate results take precedence. Docking validation must
use artifacts that pass integrity checks. Missing or corrupted validation
evidence must not be interpreted as `not_applicable`; the current coordinator
marks it as `modify`, requiring protocol repair.

1. `pocket_analysis_report.json` → strategy-aware `validation.passed == true`.
   `ligand_site` checks reference-ligand coverage; `known_site` checks known-site
   recall; `detector_consensus` checks point5 from planned tools + memory-guided
   fusion (or explicit unavailable/off fallback) + box geometry; `provided_box` checks fixed-box
   geometry and overlap with target-chain heavy atoms, without treating the
   binding site as validated.
2. `interaction_blueprint.json` → `validation.blueprint_valid == true`. Only
   constraints whose provenance, scope, and review have been confirmed make the
   corresponding hard-constraint gate applicable. PLIP∩ProLIF agreement is for
   audit only; ProLIF does not generate Blueprint constraints.
   Arpeggio/LUNA/ODDT do not promote an observation to Required. If there are no
   reliable hard constraints, the status is `not_applicable`, not failure. If
   hard constraints exist, enforce explicit `all_of/any_of/at_least_k` groups;
   do not relax them using an overall hit rate.
3. `docking_report.json` → GNINA native redocking RMSD: `<=2.0 Å` keep,
   `(2.0,2.5] Å` modify, `>2.5 Å` reject. If redocking is not explicitly enabled,
   the search-space source is inapplicable, or the route does not use GNINA,
   the status is `not_applicable`, not failure.
4. Retrospective screening → AUROC: `>=0.70` keep, `[0.50,0.70)` modify,
   `<0.50` reject. EF1%: `>=3` keep, `[1,3)` modify, `<1` reject. If either
   actives or decoys are fewer than the configured minimum, the status is modify.
5. Benchmark audit → report PR-AUC, active prevalence, `PR-AUC/prevalence`
   relative lift, and top-ranked decoy fraction together. These items have
   weight 0; do not impose a falsely universal hard cutoff.
6. `coordinator_report.json` → `round_assessment` retains protocol modify/reject,
   diagnostics, and nomination-blocking reasons. Protocol modify requires
   repairing/extending the protocol, not modifying molecules on that basis.
   Protocol reject prohibits nomination; retain non-destructive
   `rollback_instructions`. The Agent plans whether to rerun within authorization.
   Protocol outcomes do not rewrite the molecule's own keep/modify/reject classification.

## Project Final Report

After a campaign ends or pauses at an authorization/budget boundary, the Agent
**must** write `final_report.md` in the project root. This is the current
project's final deliverable; each coordination run produces only `round_report.md`.
Distinguish stage completion, protocol validity, candidate nomination eligibility,
and experimental objectives that remain unverified. Do not describe completed
computation as an experimental hit.

Base the final report on existing structured data from **all rounds**. Read each
round's available `step6/round_report.md`, `step6/coordinator_report.json`,
`step6/round_summary.json`, and `step6/candidate_ranking.json`, plus key gating
data (the report JSON for each of `step1`/`step2_5`/`step3`/`step4`/`step5`).
If a round failed early, read its `published/07_final_report/final_report.md`
and manifest. Mark missing data as not executed/unavailable. Do not omit a
final report because Step6 is missing, or rerun experiments without authorization
merely to complete the report.

Include: an experiment overview; separation of validation and candidates
(`DEC_*` is used only for benchmarks such as AUROC/EF and never enters final
candidates); three-tier AUROC/EF1 status; audits of PR-AUC + prevalence/relative
lift and top decoy fraction; a per-round summary table; final candidate molecules
(non-decoy Top-N, with their source rounds); quality and risk conclusions; the
progression of iterative improvements; and failure analysis. Do not combine
benchmark conclusions and candidate recommendations into “N hits discovered.”
Describe recommendation tables only as “candidate nominations / keep candidates,”
and remind readers that molecules with high clashes/high strain require
redocking and manual inspection.

In addition to the formal nomination table, provide separate **KEEP Top 1–5**
and **MODIFY Top 1–5** review tables for each cohort/pool. Prefer
`final_report_data.json.candidate_reviews` (or the corresponding round's
`round_report.md` / `candidate_ranking.json`). Sort by the original
`overall_score` descending; do not recalculate scores within classifications.
Retain nullable `overall_rank`, and label 1–5 as “within-group review display
order.” Record source rounds. If no round permits valid nominations, review
tables from the most recent successfully evaluated round may be shown, but
clearly mark them as review only. Entries without scores do not enter the
score-based Top5; decoys never enter these tables. Separately list issues with
molecule classifications/category metrics and missing evidence/protocol blocks.
Do not describe MODIFY or blocked KEEP entries as nominated candidates. This
requirement applies only to subsequent reporting; do not rewrite existing
experiments to fill the tables.

Preserve each round's parameters, applicable conditions, failure reasons, and
evidence references in reports. Consult those reports and artifact lineage
directly for subsequent planning.

# Fixed-config baseline without an LLM controller

Entry: `../scripts/run_fixed_pipeline.py`. This reuses normalization,
`validate_config.py`, and `run_pipeline.py`. It does not instantiate Drugmate,
read skills through an LLM, call a chat model, or change scientific choices.

## What this ablation measures

The comparison is **LLM-controlled planning/adaptation versus a predeclared
fixed policy with the same computational backend**. GNINA CNNs and other
scientific ML models remain available to both arms. Removing the LLM does not
mean removing every learned scientific model.

Two experiments answer different questions:

| Control | What to freeze | Interpretation |
|---|---|---|
| Primary fixed-policy baseline | A protocol selected on development targets before held-out evaluation | Value of LLM planning/adaptation relative to that fixed protocol |
| Matched-plan replay | An agent-produced configuration, then run independently through the same backend | Execution consistency; not an independent no-LLM planning baseline |

If inputs, configuration, seeds and backend are identical, similar results are
expected. A difference then suggests runtime nondeterminism or an integration
problem, rather than evidence of LLM scientific reasoning.

Match input visibility, target splits, tool availability/versions, evaluator,
scientific thresholds, hardware, and compute budgets. If the agent gets several
rounds but the control gets one, report that difference. A later fixed-budget
control can predeclare a seed/parameter schedule; do not choose it using held-out
results. Gold poses or activity labels must be available only in the same roles
for both arms. A gold structure used for evaluation is not automatically allowed
as a ligand-preparation fallback or search-space oracle.

## Configuration

Supply an ordinary, complete pipeline configuration. The script deliberately
does not select a universal pocket detector, protonation backend, pH, docking
engine or filtering strategy. These are part of the baseline's scientific
protocol and must be chosen before the comparison.

Required project settings (fragment, not a complete config):

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "project": {
    "data_root": "<SBDD_DATA_ROOT>",
    "project_slug": "fixed-policy-baseline",
    "execution_mode": "auto",
    "round_number": 1
  }
}
```

Every step/tool must have an explicit boolean selection, as in the agent arm.
Do not use `default_template()` as an automatically executable scientific
protocol. For an external fixed library, disable generation and all its tools;
for a supplied box, disable pocket analysis and set `search_space_source` to
`docking_box`. Prepared receptors still require final docking preparation.

`project_dir` and `upstream_project_dir` overrides are rejected: each baseline
starts fresh, without accidentally consuming the agent arm's cached results.
The entry runs exactly one outer round and preserves the coordinator's scoring
and gates. New reports expose `round_assessment` with `decision="agent_review"`;
the baseline publishes that assessment without launching another round. Legacy
`iterate` outcomes remain historical results. The deprecated coordinator
`max_rounds` parameter no longer controls campaign progression.

## Commands

Preparation and validation (does not launch the pipeline):

```bash
"${SBDD_PYTHON}" \
  "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/run_fixed_pipeline.py" \
  --config "${SBDD_DATA_ROOT}/baseline_protocol/config.json" \
  --normalization-input "${SBDD_DATA_ROOT}/baseline_protocol/inputs.json" \
  --output-dir "${SBDD_DATA_ROOT}/baseline_prepare_target01"
```

The example paths are placeholders; supply real reviewed inputs. Optional
`--normalization-input` uses the existing normalization schema unchanged.
Preparation can run Open Babel for SMILES/MOL2 conversion and RDKit coordinate
checks; it is not a computation-free static dry run. A normalization patch may
wire files but may not silently change the selected steps/tools/search route.

After reviewing the fixed protocol, execute with a **new** output directory:

```bash
"${SBDD_PYTHON}" \
  "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/run_fixed_pipeline.py" \
  --config "${SBDD_DATA_ROOT}/baseline_prepare_target01/frozen_config.json" \
  --output-dir "${SBDD_DATA_ROOT}/baseline_execute_target01" \
  --execute
```

This reuses the prepared input files and frozen defaults. Do not pass
`--normalization-input` again unless independently repeating normalization is
part of the declared experiment. Keep the preparation directory available.

For a pre-reviewed raw-input protocol, one invocation with both
`--normalization-input` and `--execute` performs the chain end to end. Launching
this headless baseline is an operator action; the wrapper does not simulate a
`sbdd_confirm` response or change the interactive Agent's gate implementation.

Both arms retain existing deterministic stage/pose retry logic, output QC,
coordinator classifications and registered stage runtimes. The wrapper does
not automatically replace failed tools, repair inputs or omit failed targets.

## Outputs and failure behavior

```text
baseline_execute_target01/
  requested_config.json
  frozen_config.json
  config_validation_report.json
  resolved_run_plan.json
  baseline_audit.json
  baseline_report.md
  pipeline_result.json             # after runner returns
  normalization/                   # when requested
  round1/
    pipeline_config.json
    pipeline_checkpoint.json
    manifest.json
    step0/ ... step7/              # enabled stages, including step2_5
    logs/
    published/                    # subject to normal publication checks
```

The audit contains input SHA256s, frozen config SHA256, actual Python-source
hashes (including uncommitted edits), execution scope and scientific outcome.
It records zero calls by this controller to an LLM, not a network-level audit
of arbitrary third-party backends. Environment/tool evidence remains in the
normal runner reports. Local source files are hashed, not made immutable.

Missing inputs, invalid configs, normalization failures, runtime failures and
incomplete/confirmation states are not experiment successes. Existing output
directories are never overwritten. The Markdown report is deterministic and
links to the evidence; it does not invent a scientific interpretation or a
missing metric. If pose evaluation is disabled, the audit says so explicitly.
`success=true` means preparation succeeded, or (with `--execute`) the configured
pipeline completed; inspect `scientific_status` and `valid_for_nomination`
separately. A plan that intentionally stops at docking is not a full SBDD run.

## Recommended first comparison

Start with a declared fixed-library docking + pose-evaluation task on a small
development subset; keep generation out of that initial comparison in both
arms. Then freeze the protocol and evaluate held-out targets with paired seeds.
Use exactly the same evaluator for both arms, and report all targets including
failed ones. Summarize completion rate, all-pose/top-ranked pose quality (with
an explicitly defined ranking), time/GPU budget, retries, manual interventions,
and LLM token/time usage in the agent arm. Do not treat unmeasured RMSD as zero.

This entry also supports generation/filtering/reporting when enabled in the
supplied config. It does not add PoseBench RMSD integration to a backend that
lacks it: an independent benchmark evaluator must be shared by both arms.

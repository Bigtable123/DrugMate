# Query → Config: Configuration Assembly

This document only assembles confirmed user facts and the selection decision from
`TOOL_SELECTION_POLICY.md` into a valid config. It does not judge which tool is
better, turn holo/apo status, artifact availability, or scientific facts into tool
switches, or supply fallback strategies.

## File Responsibilities

| File | Sole responsibility |
| --- | --- |
| `QUERY_TO_CONFIG.md` | Facts + selection decision → config field assembly |
| `TOOL_SELECTION_POLICY.md` | User objectives/facts → step, tool, strategy, and parameter selection |
| `tool_catalog.json` | Machine-readable schemas, tools, parameters, dependencies, and artifacts |
| `TOOL_CATALOG.md` | Human/Agent-readable catalog view |
| `validate_config.py` | Hard validation and resolved-plan materialization; no tool selection or experiment execution |
| `SKILL.md` | User confirmation, execution, checkpoint, and iteration contracts |
| `references/input-normalization.md` | Original user files → canonical artifacts and config merge patch |
| `references/campaign.md` | Campaign authorization, scientific stages/attempts, and Agent decision records |
| `references/benchmark-memory.md` | Unified pre-run consultation and branches for pocket, generation, and docking |

`build_tool_catalog.py` generates `tool_catalog.json`, `TOOL_CATALOG.md`, and
`TOOL_SELECTION_POLICY.md` from registry and validator contracts. The Agent makes
all fact → tool decisions according to the policy.

## Preparation Parameters and Stage Mapping

Read the [receptor preparation contract](references/preparation-contract.md) first.
Set `structure_prep.pdbfixer.add_hydrogens=false`. Residue repair is disabled by
default; enabling it requires a short-gap repair rationale. Diagnostic PROPKA and
provisional Reduce are disabled by default. Put the final backend, pH rationale,
water strategy, and candidate-state rationale in
`steps.docking.tools.docking_prep.params`, and backend switches in
`steps.docking.tools`. Do not copy old Step 0 Meeko, final hydrogen-addition, or
water-selection switches into new configs. Set `context.chain` from the target
entity/author-chain mapping; the example A is not a default. Enabling
`pocket_analysis.prolif` requires `structure_prep.reduce` with `build=true` and a
`reason`, or an explicitly supplied, reviewed `input_artifacts.provisional_receptor`.
The hydrogen-free baseline does not meet ProLIF's input requirements.

## Assembly Inputs

Pose evaluation defaults to
`steps.pose_evaluation.params.selection_strategy="ranked_per_compound"`.
The maximum preselected poses per molecule depends on batch size: 20 for ≤50
molecules, 5 for 51–1000, and 1 for >1000. `poses_per_compound` may lower this limit;
explicit legacy strategies retain their behavior. Set PoseCheck's reference-conformer
count in `steps.pose_evaluation.tools.posecheck.params.num_conformers` (default `10`),
not docking's `num_modes`. `max_poses` is the actual evaluation budget and may produce
partial coverage. Explicitly choose `all` only to compare every pose. The default
strategy evaluates preselected poses, permitting later `mdanalysis` or `mmgbsa`
postprocessing to update their verdicts.

Obtain user facts and the selection decision before assembly. Local receptor,
library, or docking-box inputs require normalization records even if already PDB/3D
SDF. Store new campaigns' cross-round decisions and authorization separately;
do not insert explanations into tool configs:

1. Confirmed user/config facts, including project, context, and external artifacts.
2. The complete policy selection decision: every step/tool switch, explicit
   site/search strategy, parameter overrides, and rationale.
3. Successful `input-normalization` outputs: `normalization_report.json` records
   conversion, validation, the failure ledger, and provenance; `config_patch.json`
   supplies only canonical paths and producer-disable constraints.
4. `memory_consultation.md` covering all three domains and applicable selector
   outputs, or explicit non-applicability/disabled/fallback reasons. Unknown-site
   selection calls `plan-pocket-tools` before finalizing the matrix; Step 1
   verifies the exact `execution_tools` and fuses their native Top-3 candidates
   using historical reliability. Fixed DeepPocket#1 four-tool rank fusion is
   the unavailable-memory/off fallback. Declare memory mode and prepare its
   query encoding/index before the run; do not preselect a box.
5. Campaign scientific stage, attempt, user authorization, and Agent decision draft,
   saved as separate JSON under the [campaign protocol](references/campaign.md).
   Legacy single-round configs may omit them.

Facts and selection decisions are Agent working records. Only registered catalog
keys belong in config; intent, rationale, uncertainties, and rejected alternatives
must not masquerade as top-level fields.

## Top-Level Mapping

```text
Project/round/execution mode          → project.*
Scientific stage/attempt/decision     → project.campaign_stage / campaign_attempt / round_decision_file
Campaign state reference             → project.campaign_state_file
Input identity/shared knowledge      → context.*
Docking memory ablation mode         → context.docking_memory_mode
Metal/covalency/protonation facts     → context.scientific_facts.*
External canonical artifact          → input_artifacts.<artifact_name>
Step switch                          → steps.<step_id>.enabled
Tool switch                          → steps.<step_id>.tools.<tool_id>.enabled
Step parameter                       → steps.<step_id>.params.<param>
Tool parameter                       → steps.<step_id>.tools.<tool_id>.params.<param>
Resource queue budget/shared capacity → runtime_resources.*
```

See the [execution lifecycle](references/execution-lifecycle.md) for optional
`runtime_resources` fields, default immediate return when resources are busy, and
shared-directory constraints. Queue budgets differ from `step_timeouts`; larger
execution timeouts must not hide unavailable resources.

`context.docking_memory_mode` accepts `off`, `legacy`, `ligand_only`, or
`ligand_protein` (default). Save the mode explicitly during planning and read it
with `scripts/benchmark_memory.py select-docking --config <config.json>`.
If `--mode` is also supplied, it must match config. The mode controls only Agent
memory consultation/selection; it does not replace
`steps.docking.tools.<engine>.enabled`, and the runner does not switch engines from it.
Once inputs are ready, obtain recommendations, update the full explicit tool matrix,
validate, and present Gate A. User-specified engines take precedence. Record suggestions,
actual choices, override/fallback reasons, and index versions in
`memory_consultation.md`, not unregistered config fields. See the
[memory protocol](references/benchmark-memory.md) for commands and ablation boundaries.

Docking's `--task`/`--metric` and generation's `--priority-metric` are pre-run selector
arguments, recorded in consultation and query JSON. Map generation recommendations
to registered `steps.molecule_generation.tools`, preserving total candidate budgets
and source-input requirements. Pocket selection runs inside main Step 1: set
`steps.pocket_analysis.params.pocket_memory_mode` (`auto`, `off`, or `required`),
`pocket_memory_index`, `pocket_memory_query_npz`, and
`pocket_memory_query_metadata`. Known identities and exclusions use the registered
`pocket_memory_query_case/group` and `pocket_memory_exclude_cases/groups` fields.
Use `site_selection_strategy=detector_consensus` and downstream
`search_space_source=docking_box`, retaining the chosen anchor's residues.
Before Gate A, call `plan-pocket-tools` with the same query options and set
the detector switches to exactly its `execution_tools`; Step 1 rechecks this
plan before launching any detector.
Do not route a memory-exported box back through `provided_box`, which would
recompute the residue set. Preserve the original fusion decision for audit.

The canonical scientific-fact fields are:

```json
{
  "context": {
    "scientific_facts": {
      "has_metal": true,
      "is_covalent": false,
      "protonation_sensitive": true
    }
  }
}
```

Library provenance affects the default gating mode. For `input_artifacts.library`
randomly sampled outside the pipeline from a large user library, set
`context.candidate_source="sampled_library"` (V2). A small user library used in full
uses `user_supplied_full` (V1). Do not infer sampling from artifact producer names.
Explicit `steps.<gate_step>.params.screening_mode` always takes precedence.

Benchmark applicability requires actual labels, not `sampled_library` or directory
names. Follow the [label-input contract](references/benchmark-labels.md): assign
`ACT_<original-ID>` / `DEC_<original-ID>` before normalization from verifiable sources,
retaining original-ID mappings and provenance. Do not prefix unknown-activity molecules
with a class. Disable `coordination.screening_benchmark` if either class is absent.
Enabling it with an external SDF missing labels fails config validation. The coordinator
still checks generated/reused libraries; missing labels/scores mean protocol `modify`,
not silent gate omission.

Set scientific facts to `true`/`false` only with reliable evidence; omit unknowns.
These fields provide scientific inputs and validator capability/safety evidence,
not tool selection.

## Complete Explicit Switch Matrix

The final config must not rely on registry defaults to guess this round's tools:

- Explicit `enabled: true|false` for every registered step.
- Explicit `enabled: true|false` for every registered tool in each step.
- Normally set all tools under disabled steps to `false`, avoiding ignored opt-ins.
- The sole exception is an audited upstream producer before `start_step`: preserve
  its complete previous tool/parameter matrix for lineage/cache checks; range
  restrictions disable execution only.
- `params` may contain only policy-required overrides; other fields use schema defaults.

A step has this shape. Angle brackets are placeholders, not literal input:

```json
{
  "steps": {
    "<registered_step_id>": {
      "enabled": true,
      "params": {"<registered_step_param>": "<policy_value>"},
      "tools": {
        "<registered_tool_a>": {"enabled": true, "params": {}},
        "<registered_tool_b>": {"enabled": false, "params": {}}
      }
    }
  }
}
```

Expand every registered tool ID for the step; the example does not require only two tools.

## Explicit Site/Search Strategies

The policy selects both strategies that change the scientific path; assemble them unchanged:

```json
{
  "steps": {
    "pocket_analysis": {
      "params": {
        "site_selection_strategy": "<ligand_site|known_site|detector_consensus|provided_box>"
      }
    },
    "docking": {
      "params": {
        "search_space_source": "<reference_ligand|docking_box|whole_receptor|precomputed_grid>"
      }
    }
  }
}
```

This specifies field locations and choices, not when to choose them. Do not change
strategies during assembly based on reference-ligand presence. If artifact prerequisites
are unmet, revise the plan at the policy/validator layer; do not introduce an assembly
fallback. In later `start_step` rounds, excluded Step 1 retains its previous
`site_selection_strategy` for producer fingerprint/lineage only. It neither defines
newly injected docking search space nor executes in that round.

## Config Skeleton

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "project": {
    "project_slug": "example-project",
    "data_root": "<SBDD_DATA_ROOT>",
    "round_number": 1,
    "execution_mode": "interrupt"
  },
  "context": {
    "pdb_id": "<optional-pdb-id>",
    "receptor_file": "<optional-local-receptor>",
    "target_name": "<optional-target-name>",
    "chain": "A",
    "scientific_facts": {}
  },
  "steps": {
    "<expand-every-registered-step>": {
      "enabled": false,
      "params": {},
      "tools": {"<expand-every-registered-tool>": {"enabled": false, "params": {}}}
    }
  },
  "input_artifacts": {},
  "step_timeouts": {}
}
```

This is a structural illustration, not an executable config. Copy IDs, required
fields, types, and choices only from the current catalog. Do not add runtime/Conda
paths or manually construct inter-step artifact paths.

## Execution Range and Later Rounds

`project.round_number` controls physical directories, not R1/R2/R3 identity.
Campaigns use `project.campaign_stage` and `project.campaign_attempt`: R1 repair
retains R1 and increments the attempt; checkpoint resume does not. Fill campaign
fields together under the protocol; do not infer stages for legacy configs.

`start_step` / `stop_step` are execution-range conveniences, not replacements for
the full switch matrix. Resolved enabled states must match the range. Set
`project.upstream_project_dir` to the actual selected upstream attempt and modify
its full config; R3 selecting R1 conditions need not use R2 as upstream. Preserve
pre-`start_step` producers' original step/tool/parameter matrices. Materialization
sets their resolved enabled states to false, retaining raw matrices for lineage/cache
checks only. Unexecuted tools after `stop_step` do not receive this exception.

An upstream step may therefore retain `enabled=true` in raw JSON but be
`enabled=false` in `resolved_plan`: these represent raw lineage versus effective
execution, not competing choices. Inspect execution states in the dry-run's `resolved_plan`.

## Inject External Artifacts

Replace `<INPUTS_DIR>` with the actual input directory before saving this JSON;
JSON does not expand shell variables.

```json
{
  "input_artifacts": {
    "library": "<INPUTS_DIR>/modified_library.sdf"
  }
}
```

Allowed names are in catalog `injectable_artifacts`. Explicitly disable a same-name
producer or exclude it from execution; otherwise validation rejects the ambiguity
of both injecting and producing. Injection changes provenance, not tool/strategy switches.
Read [artifact reuse and recovery](references/artifact-reuse.md) before injection;
identify affected upstream/downstream evidence and do not combine new inputs with
invalidated old artifacts.

Use `input_artifacts.library` for a raw external library, starting at
`property_prefilter` when prefiltering remains necessary. Only completed, audited
prefiltered libraries use `input_artifacts.prefiltered_library`. Do not inject both
and leave the runner to guess.

New R3 generated/modified libraries must enter as `library` at `property_prefilter`
(order 3) or earlier with prefiltering enabled; `start_step: 4` must not bypass new
molecules' property/identity reevaluation. Inject the library's own `library_manifest`,
not an old producer's. Explicitly use
`context.candidate_source="molecule_generation"` or V3 at each relevant gate when
including generated/modified molecules, preventing external SDF from defaulting to
full user-library V1. Fixed-library replication retains original provenance and modes.
Reevaluate merged historical KEEP and new molecules under the selected common protocol;
do not reuse old profiles' scores directly. R2 independent replication requires new
seeds, replicate IDs, and actual rerun steps; incrementing `round_number` or reusing
cached output is not replication.

External libraries use SDF paths. Convert SMI/SMILES directly through Open Babel
under `references/input-normalization.md`; copy SDF unchanged. Conversion adds no
unique-ID, manifest, count, coordinate, or chemical-state acceptance requirements,
and does not fill or rename records. Convert and verify MOL2 receptors as PDB.
Only explicit user designation as prepared permits separate `prepared_receptor_noH`,
water inventory, and prep-manifest injection; otherwise use `context.receptor_file`
and Step 0 prep-only processing.

## Pre-Run Config Validation

Use `validate_config.py --input-file` from [Quick Start](SKILL.md#quick-start).
Store configs and validation logs under the configured `SBDD_DATA_ROOT`
according to the entry-point rules.

- Exit 0: `valid=true`. Compare each step's `resolved_plan` / `plan_summary` with
  the policy decision; confirm tool switches and strategies were not changed.
- Exit 1: use `errors` to repair facts, policy, or assembly and retry; do not bypass constraints.
- `project.upstream_project_dir` implies static upstream-artifact availability;
  `--assume-upstream` explicitly enables similar assumptions. Execution still checks
  files and reuse eligibility; static validity does not verify upstream evidence.

A valid dry-run does not authorize execution. Use the actual campaign CLI to create
a plan record matching the config, then check explicit preauthorization under Gate A.
Within scope, record the basis and execute; otherwise request additional authorization.
After config changes, regenerate the plan, validate, and recheck authorization;
do not continue with the old config hash.

## Regenerate Machine-Readable Materials

After registry or validator contracts change, run:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/build_tool_catalog.py"
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/build_tool_catalog.py" --check
```

The generator maintains `tool_catalog.json`, `TOOL_CATALOG.md`, and
`TOOL_SELECTION_POLICY.md` together. It does not execute selection or choose tools
on the Agent's behalf.

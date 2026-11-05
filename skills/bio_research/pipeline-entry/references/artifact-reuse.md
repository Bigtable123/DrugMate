# Artifact Reuse and Recovery

Read this reference when injecting external artifacts, investigating rejected
cross-round caches, rewinding checkpoints, or diagnosing invalid published
directories. Normal configuration and reuse boundaries are defined in the
[entry skill](../SKILL.md#reuse-across-rounds-and-execution-scope); configuration
fields are documented in [QUERY_TO_CONFIG.md](../QUERY_TO_CONFIG.md).

## Content and Implementation Jointly Determine Reuse Eligibility

Cache contract v2 records the producer's tools, parameters, relevant context,
implementation fingerprints for the skill, shared code, helpers from other
skills, and enabled external tools' `source_dir`, together with SHA-256 hashes
of `artifact_inputs` and `artifact_outputs`. Matching filenames do not establish
reuse eligibility. Old caches without v2 lineage are rejected and rerun. See
the current round's `logs/artifact_lineage.json` for detailed decisions.

The target name derived from metadata in Step 0 is also a formal artifact, but
it enters Step 2's actual dependency chain only when `context.target_name` is
not explicitly set. This avoids reruns caused by irrelevant metadata changes.
Each step also rechecks input content before publishing outputs, preventing
mixed results from being published after an input is rewritten during execution.

## Replacement and Invalidation Propagation

Replacing or injecting an upstream artifact, rerunning its producer, rewriting
an output file, or encountering a newer `result.json` that is corrupt or fails
validation invalidates the old producer bundle and its downstream derivatives.
Do not silently fall back to an older round.

An injection replaces the producer's entire previous lineage. Sibling outputs
not supplied in the same injection and their downstream derivatives also become
invalid. For example, after injecting a new `docking_box`, the old pocket
features, interaction blueprint, and pocket report are no longer reusable
evidence. Plans depending on them may fail static validation or actual preloading.

- Reuse the same artifact directly; repeated injection is unnecessary.
- When the pocket or input definition genuinely changes, rerun from a stage
  capable of producing consistent evidence, or explicitly disable consumers
  that no longer apply in the revised plan, then validate and confirm it again.
- Do not combine a new box with old pocket evidence as results of one experiment.

In `provided_box` mode, Step 1 consumes an external `docking_box` and publishes
only the corresponding pocket definition and analysis report. It neither takes
ownership of nor overwrites the box. The box's SHA-256 is recorded in lineage
as an actual input. Changing it still invalidates the old report and downstream
results, requiring the analysis to run again.

## Static Validation Versus Actual Recovery

When `upstream_project_dir` is configured or `--assume-upstream` is supplied,
`validate_config.py` optimistically assumes that upstream artifacts from disabled
steps exist; it does not read actual upstream lineage. The runner checks actual
artifacts, content, and producer configuration, and permits preloading only from
disabled producers before the `start_step` boundary. Merely disabling a step can
therefore pass a dry run while leaving execution without a usable artifact.

Normal cache reuse and the processed step prefix of an interactive checkpoint
are separate mechanisms. When PocketMiner is enabled, an ordinary pending
Step 1 cannot reuse its cache. However, resuming the same round after confirmation
can retain the verified completed prefix. A pause should not cause a new project
or repeated execution. Apply the entry skill's Gate B rewind rules when rewinding.

## Published Directory Integrity

`published/` is fully generated in a sibling staging directory and its manifest
and hashes are validated before it replaces the official directory. Every stage
start/end, pause, and catchable failure updates the `pipeline_partial` view with
the artifact index, status, and logs of currently valid stages. Failure does not
remove the existing results view or rewrite the experiment's failure status.
On rewind, the old published directory moves into the current round's
`.publication_archive/` for audit; the new view does not reference old artifacts
from invalidated steps. If building a partial publication fails, the old view is
also isolated and `publication_error` is recorded, preventing continued display
of stale results. Only `pipeline_final` must still meet the entry skill's defined
coordination/reporting completion boundary.

For historical failed or partially executed rounds, use
`reflective-coordinator/scripts/republish_existing_round.py --round-dir <roundN> --partial`
to build a view retrospectively. It must acquire the existing `.pipeline.lock`,
refuses to overwrite an existing publication, and does not rerun experiments or
rewrite original artifacts. A partial view validates only its own files; original
experiment files are indexed by existence and size, without claiming scientific
or provenance-integrity validation. `kill -9` and power loss cannot trigger a final
publication, so only the latest snapshot remains. A view cannot be guaranteed
when the disk is not writable.

When a protocol is invalid, apply non-destructive rollback/invalidation markers
and preserve molecules, poses, reports, and lineage. Artifact invalidation does
not authorize deletion of original experiment files.

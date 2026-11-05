# Cancellation, Resource Waits, and Failure Cleanup

## Managed Launch and Stop

Use the local `terminal` with `background=true`, `cancel_on_interrupt=true`,
`notify_on_complete=true`, and `pty=false`. Write the approved entry command as
`exec "${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/run_pipeline.py" --input-file <absolute-config-path> --project-dir <current-round-directory>`.
Save the returned `session_id`; use `process` to read progress or request a stop.
Do not detach the command with `nohup` or shell `&`. Input-normalization commands
must use the same background parameters and `exec`. This general cancellation
parameter requires a Drugmate version containing the corresponding runtime
change. If an older version does not support it, update the running Drugmate
first. Copying the skill alone does not establish that interrupting the agent
will automatically stop background tasks.

The runner catches SIGINT/SIGTERM/SIGHUP. Cancellation or a step timeout
terminates descendant processes belonging to this command, including tools that
start a separate session, then saves logs, a partial artifact index, and an
execution report. Do not kill processes in bulk by program name, or stop another
experiment because it occupies resources. Repeated stop signals do not interrupt
the first cleanup. SIGKILL and machine power loss cannot run cleanup; the latest
published snapshot remains readable, but its status may still be running.

Query status with `process(action="poll", session_id=...)`. To wait, call
`process(action="wait", session_id=..., timeout=30)`, explicitly limiting each
wait to 1–30 seconds. Do not use unbounded shell PID loops such as
`while kill -0 PID; do sleep 30; done`, or omit the `wait` timeout. After each wait,
read actual progress, handle new input, and update the user promptly. Do not
continue waiting blindly merely because a process is alive. A timeout returned
by `process` wait only ends that observation window; it does not justify killing
the whole library or experiment. Do not impose a fixed 8-hour kill on the entire
queue-plus-execution lifecycle through the shell, agent, or frontend. Execution
limits are managed by the runner's `step_timeouts` and individual tool budgets.
If a frontend has a hard total time limit, use the default immediate return when
resources are busy.

## Input-Conversion Progress and Per-Molecule Timeouts

SMILES conversion in `input-normalization/run.py` defaults to at most 8 concurrent
tasks, all invoking Open Babel `-h --gen3d fast`. The per-molecule budget is
controlled by `library.conversion_timeout_seconds` in the normalization request:
the default is 60 seconds, and omission or `null` still means 60 seconds. Supply
a positive integer to override it; unlimited conversion is unsupported. When a
molecule exceeds its budget, the script terminates that conversion, records
`conversion_timeout`, and processes the remaining molecules. Timed-out molecules
retain a reason and empty scores, just like other failed molecules. Receptor
format conversion preserves the original coordinates and does not add 3D
generation arguments.

Before converting molecules, a fast-feature probe checks native Open Babel.
If it detects a coordinate write-back defect, the script compiles and caches a
compatibility program on first use under the current round's
`inputs/normalization/runtime/`. This preparation requires a C++17 compiler,
matching Open Babel development headers/libraries, and Eigen 3. It does not
change the global environment or switch the fast algorithm. A working native
program requires no compilation. Check `library.conversion_execution.command_backend`,
the runtime probe files, and `fast_compat_build.json` to establish the actual
backend. Preparation errors are environment problems for the entire run. Do not
count a failed probe as a user molecule's failure or hide it by skipping library
molecules.

During execution, use `read_file` on the current round's
`inputs/normalization/normalization_progress.json`. Inspect `completed_records` /
`total_records`, `retained_records`, `skipped_records`, and `running_records`
(molecule names, source-record positions, and elapsed time), together with
`status` and `process` output, to update the user. Progress updates follow task
completion; the final SDF is still merged in original input order. If progress
stops updating for an extended period, inspect these records and diagnostics from
`process(action="log", session_id=...)` to locate and report the stall. Do not
create another PID-wait script or stop the entire library because one 30-second
observation window expires. On completion notification, read
`normalization_report.json` and merge `config_patch.json` only when `success=true`.
If every molecule fails, normalization still fails; do not start a pipeline
that depends on that library.

## Resource Policy

Optional top-level configuration:

```json
{
  "runtime_resources": {
    "wait_timeout_seconds": 0,
    "slots": 3
  }
}
```

- Default `wait_timeout_seconds=0`: when resources are busy, immediately return
  `resource_unavailable` with `retryable=true`, preserve the report, and tell
  the user that no execution resources are currently available and to retry later.
- To permit queuing, set a separate duration in seconds, such as `3600`. Resource
  waiting happens before preflight/step launch and does not consume the step
  timeout. Expired waits still return resource unavailability, which must not be
  recorded as scientific experiment failure.
- Default `slots=3` admits at most 3 concurrent experiments on the same host;
  bounded parallelism may still occur inside each experiment. Change slots
  consistently only after confirming host capacity. All participants must use
  the same capacity setting.
- The default lock directory is `project.data_root/_resources/<host-id>/`, shared
  across projects. If different entrypoints use different data roots, explicitly
  set the same absolute `runtime_resources.lock_directory`. Do not create a
  separate resource lock for each project or delete lock files in use. An
  administrator may change the `capacity.lock` capacity record only after all
  experiments have stopped.

This is cooperative concurrency limiting among experiments using this runner,
not a GPU/memory scheduler. It neither probes nor terminates other programs.
The shared directory must be readable and writable by participating accounts.
An experiment releases its lease when pausing for confirmation or exiting;
the operating system also releases flock after process exit. The frontend reads
`manifest.json` and `resource_status.json` to distinguish queuing from execution.
Do not count resource waiting as molecular computation time.

## Publication and Recovery

Normal failure, timeout, interruption, and resource unavailability all update
`published/`. The entrypoint returns `execution_report` (Markdown) and
`published_dir`; HTML and JSON reports are also available in that directory.
Scientific data remain in the working round, and the publication index links to
those files. Preserve both the working round and publication to keep the links
usable. A failed stage registers only files added/changed by that execution or
explicitly declared in its result, not leftover files from before a rewind as
new results.

`publication_mode=pipeline_partial` means partial evidence can be read. It does
not establish success of the complete pipeline or permit nomination of molecules
whose scientific evaluation is incomplete. After stopping, report existing results
and the failure reason. Do not rerun scoring just to generate a report or invent
missing metrics. After resolving the problem, resume with the same `--project-dir`
under the original confirmation contract.

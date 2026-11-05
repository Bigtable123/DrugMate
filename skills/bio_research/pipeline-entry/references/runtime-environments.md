# Runtime Environment Troubleshooting

Read this reference after `runtime_preflight_failed`, when a Conda path is stale,
or when invoking an individual stage tool. Routine execution rules are in the
[entry skill](../SKILL.md#runtime-environments-separate-the-main-and-tool-environments).
For path migration, read the [deployment configuration](deployment.md) first:
`SBDD_*` variables select environment and model directories for the current host.

## Environment Declarations and Launch

Use `TOOL_ENVIRONMENTS`, step runtime, and tool runtime declarations in
`../../shared/pipeline_registry.py` as the authority;
`../../shared/tool_runtime.py` resolves and launches them. Choosing a stage
script's interpreter and choosing the environments of tools within that stage
are separate decisions.

The following current declarations illustrate the mapping; they are not a
complete environment inventory to maintain manually:

| Component | Declared environment |
|---|---|
| Stage scripts with explicit runtime declarations | `drugmate` |
| coordination / reporting stage scripts | No stage runtime declared; inherit the entrypoint interpreter |
| UniDock | Environment name `unidock` |
| Pocket2Mol | Environment name `Pocket2Mol`, distinct from registry ID `pocket2mol` |
| RTMScore | Configured prefix `${SBDD_ENV_RTMSCORE_PREFIX}` |

The runtime layer uses `conda run -n <name>` or `conda run -p <prefix>` when needed.
It can launch directly when already in a matching environment. For fixed-prefix
environments, it can also use a checked executable inside the prefix when a
Conda frontend is unavailable. Named environments are not replaced by guessed
filesystem paths.

## Conda Migration and Stale Activation State

After Conda is moved, IDEs, terminals, or agent processes started before the move
may retain `CONDA_*` variables pointing to a deleted prefix. The runtime layer
clears stale activation state before launching the Conda frontend. If the
entrypoint itself still uses the old interpreter, reload/restart first and launch
it with the configured drugmate scientific environment's interpreter.

## Diagnose the Specific Tool

1. Read `missing`, `errors`, and the corresponding `tools` records in the current
   round's `tool_runtime_report.json` to identify the failed step/tool and its
   probe results.
2. Check the environment name or prefix, commands, resources, and probe conditions
   against the registry. An existing environment directory or a successful import
   in the main environment does not establish that the target tool's execution
   requirements are met.
3. The runner preflights only enabled tools that declare a runtime in stages
   pending execution. During normal pause/resume, completed steps retain their
   accumulated evidence and are not probed again merely because confirmation
   caused a pause.
4. `project.python_executable` (or a top-level key with the same name) overrides
   stage interpreters ahead of runtime selection, bypassing that selection.
   Normal experiment configs must not set it; it also does not uniformly change
   the individual tools' runtimes.

Preserve failure evidence. Do not silently switch environments, change tool
selection, or mix dependencies. Files and logs generated while troubleshooting
must remain under the deployed `SBDD_DATA_ROOT`, following the entry skill's
output-directory requirements.

## Reduce2 Chemical Data Dependencies

The `reduce` environment requires the official CCTBX `chem_data` package. It
provides geometry restraints in `chem_data/mon_lib` and `chem_data/geostd`, and
CCD chemical components in `chem_data/chemical_components`; CCTBX locates them
automatically. See the [official CCTBX instructions](https://github.com/cctbx/cctbx_project#monomer-library)
for installation. Use a release matching the installed CCTBX version and verify
its SHA-256.

Do not rename geostd to `chemical_components` or point `MMTBX_CCP4_MONOMER_LIB`
at CCD. When a custom restraint library is needed, that variable must point to
actual mon_lib/geostd data; CCD still requires a separate installation. Preflight
checks both restraint-library initialization and actual CCD → RDKit molecule
construction for ALA. Reading a nonempty CIF dictionary alone is insufficient
evidence of usability.

For `NoneType.iterrows`, missing restraints for standard residues, or
`CCTBX_ASSERT(j_sym >= 0)`, first check the chemical-data source and format. Do not respond by directly translating
the receptor, changing the unit cell, ignoring missing restraints, or disabling
flips. Acceptance after repair must include hydrogen addition to the actual
receptor, heavy-atom retention, and coordinate checks.

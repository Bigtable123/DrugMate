# Deployment Paths and Credentials

See the [README](../../../../README.md) for basic installation. This reference is
for deployments whose scientific tools are already installed and that need to
configure server paths or multiple Conda environments. It does not change the
scientific protocol, gates, or authorization scope. Set these variables before
starting the Python CLI so child processes inherit the same configuration;
restart the CLI after changes.

## Local deployment file

Save the following as a local `deployment.sh` outside the source checkout,
replacing example paths with actual installation locations. This file only sets
paths; it does not install tools. Add only the variables your existing deployment
needs.

```bash
# Resolve the checkout directory passed when sourcing this file.
drugmate_checkout_root="$(cd "${1:?Pass the checkout directory}" && pwd)" || return 1
export BIO_SKILLS_ROOT="$drugmate_checkout_root/skills/bio_research"

# Small indexes distributed with the source and separately downloaded bundles.
export SBDD_MEMORY_ROOT="$drugmate_checkout_root/memory"
export SBDD_MEMORY_BUNDLES_ROOT="$drugmate_checkout_root/memory/bundles"

# Keep all experiment outputs outside the source checkout and frozen data.
export SBDD_DATA_ROOT="$HOME/drugmate-runs"

# The installed scientific main environment, separate from the agent application.
# Do not also set SBDD_ENV_DRUGMATE_NAME when using PREFIX.
export SBDD_ENV_DRUGMATE_PREFIX="$HOME/.conda/envs/drugmate-science"
export SBDD_PYTHON="$SBDD_ENV_DRUGMATE_PREFIX/bin/python"

# Set these shared roots only when tools are already deployed at these locations.
# export SBDD_SOURCE_ROOT="$HOME/drugmate-tools"
# export SBDD_EXTERNAL_ROOT="$HOME/drugmate-tools"
# export SBDD_CONDA_ENVS_ROOT="$HOME/.conda/envs"

# Per-tool overrides take precedence over shared roots; enable as needed.
# export SBDD_POCKET2MOL_SOURCE_DIR="$HOME/drugmate-tools/Pocket2Mol"
# export SBDD_ENV_POCKET2MOL_PREFIX="$HOME/.conda/envs/Pocket2Mol"
# export SBDD_ENV_TARGETDIFF_NAME="targetdiff-gpu"
# export SBDD_AIZYNTHFINDER_CONFIG="$HOME/drugmate-tools/models/aizynthfinder/config.yml"
# export SBDD_BENCHMARK_ROOT="$HOME/drugmate-benchmarks"
```

Load it from the source root, then start the CLI in the activated application
environment:

```bash
export DEPLOYMENT_FILE="$HOME/.config/drugmate/deployment.sh"
source "$DEPLOYMENT_FILE" "$PWD"
python drugmate --cli
```

The script sets `BIO_SKILLS_ROOT` from its argument. If variables are supplied
through `$DRUGMATE_HOME/.env` or a service configuration, specify the directory's
actual absolute path there; those mechanisms do not execute `deployment.sh`.
These are host-deployment variables, not experiment JSON fields.

## Deployment variables

| Variable | Meaning and default behavior |
| --- | --- |
| `BIO_SKILLS_ROOT` | Absolute path to this checkout's `skills/bio_research`, used by command examples |
| `SBDD_DATA_ROOT` | Experiment output root for configs, normalized inputs, caches, locks, logs, and reports; defaults to `$HOME/drugmate-runs` |
| `SBDD_MEMORY_ROOT` | Root of the supplied selector indexes; defaults to the active source or installed package's `memory/`. An explicit `--selector-index` takes precedence |
| `SBDD_MEMORY_BUNDLES_ROOT` | Parent of the downloaded `drugmate-generation-memory/` and `drugmate-docking-memory/` directories; defaults to `bundles/` below the resolved memory root |
| `SBDD_SOURCE_ROOT` | Root for tool sources and model inputs; defaults to `$HOME/drugmate-tools` and does not move experiment outputs |
| `SBDD_EXTERNAL_ROOT` | Additional tool source root; defaults to `$HOME/drugmate-tools` |
| `SBDD_CONDA_ENVS_ROOT` | Shared parent directory of registered scientific environments, retaining their actual registry environment names |
| `SBDD_ENV_<ID>_PREFIX` | Absolute directory of one registered environment; takes precedence over the shared environment root |
| `SBDD_ENV_<ID>_NAME` | Conda name of one registered environment; takes precedence over the shared root and cannot be set together with that environment's `_PREFIX` |
| `SBDD_<TOOL>_SOURCE_DIR` | Source directory of one tool; relocates its default entry point and checkpoint together |
| `SBDD_AIZYNTHFINDER_CONFIG` | Absolute path to the deployed AiZynthFinder configuration |
| `SBDD_BENCHMARK_ROOT` | External benchmark-resource root; defaults to `$HOME/drugmate-benchmarks` |
| `CONDA_EXE` | Installed Conda executable; if unset, discovered through `PATH` and common installation locations |
| `SBDD_PYTHON` | Shell shorthand for the scientific orchestration interpreter, which must be chosen explicitly; does not override individual tool environments |

Path variables must contain nonempty absolute paths; `~` expansion is supported.
Unset roots use the portable defaults above. `<ID>` is the uppercase
registry environment ID, for example `DRUGMATE`, `VINA`, `REDUCE`, `POCKET2MOL`,
`RTMSCORE`, `ADFRSUITE`, or `MMGBSA`. The actual environment directory may be
named `targetdiff-gpu` or `rtmscore_clean`; do not infer it from a tool name.

`drugmate` is the logical ID of the scientific main environment. Select the actual installation with
`SBDD_ENV_DRUGMATE_PREFIX` or `_NAME`; this does not rename an existing environment.
`_NAME` does not automatically derive `SBDD_PYTHON`, which must still point to that
environment's actual interpreter. Do not use `project.python_executable` to merge
or bypass tool environments.

Supported source-override `<TOOL>` values are `POCKETMINER`, `POCKET2MOL`,
`DIFFLINKER`, `DIFFSBDD`, `REINVENT4`, `RESGEN`, `TARGETDIFF`, `DEEPPOCKET`,
`DIFFDOCK`, `DOCK6`, `LINGO3DMOL`, `TAMGEN`, `MOLCRAFT`, `RTMSCORE`, and `CHEMFH`.
The MolCRAFT override must point to the inner source directory containing
`sample_for_pocket_local.py`. Explicit tool parameters are not rewritten by
default-path remapping.

## Memory and experiment paths

For bundle extraction, see the [unified memory protocol](benchmark-memory.md#downloads-and-deployment-location).
`memory/bundles/` is excluded from Git and Python packages. Keep downloaded bundle
bytes, version identities and content hashes unchanged. Small indexes, large bundles, and
experiment outputs occupy three different locations and are not interchangeable.
For encoder settings such as `SAPROT_MODEL` and `SAPROT_FOLDSEEK`, see the
[encoder entry-point guide](../../../../memory/saprot_embedding/README.md).
Setting the main environment path does not install these resources.

To inspect resolved directories and environments without running tool probes:

```bash
python3 "${BIO_SKILLS_ROOT}/shared/deployment_config.py"
```

Explicitly set the actual `project.data_root` in every experiment config; JSON
does not expand shell variables. Allocate directories through the pipeline entry's
`--init-project-dir` and request-ID contract before writing all experiment outputs
under that root. Resolve deployment variables into actual paths before writing
JSON configuration. Read-only inputs may remain external. Path changes do not alter
Gate A/Gate B, checkpoint, campaign, or artifact-lineage requirements, and do not
bypass tool preflight, model-hash, or source-version checks. Inspecting
configuration does not establish tool availability; execution failures must still
be diagnosed through `tool_runtime_report.json`. The committed tool catalog
uses deployment-variable templates; do not replace them with local machine
paths in committed files.

## Credentials

Keep model/provider API keys in the operator's local Drugmate configuration or
process environment, never in experiment JSON, this skill, or the release
repository. Scientific data routes such as ChEMBL, RCSB, and ProteinsPlus may
require network access; they are separate from the removed general-purpose
network plugins. Scientific environments, weights, frozen memory data, and
compound libraries are provisioned separately by the deployer.

"""Portable deployment paths, separate from scientific run configurations.

Set deployment variables before starting the CLI or a pipeline process.
Resolving paths never creates directories, probes tools or reads credentials.
"""
from __future__ import annotations

import json
import os
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
import re


SOURCE_LOCATIONS = {
    "POCKETMINER": ("SBDD_SOURCE_ROOT", "gvp"),
    "POCKET2MOL": ("SBDD_SOURCE_ROOT", "Pocket2Mol"),
    "DIFFLINKER": ("SBDD_SOURCE_ROOT", "DiffLinker"),
    "DIFFSBDD": ("SBDD_SOURCE_ROOT", "DiffSBDD"),
    "REINVENT4": ("SBDD_SOURCE_ROOT", "REINVENT4"),
    "RESGEN": ("SBDD_SOURCE_ROOT", "ResGen"),
    "TARGETDIFF": ("SBDD_SOURCE_ROOT", "targetdiff"),
    "DEEPPOCKET": ("SBDD_SOURCE_ROOT", "DeepPocket"),
    "DIFFDOCK": ("SBDD_EXTERNAL_ROOT", "DiffDock"),
    "DOCK6": ("SBDD_EXTERNAL_ROOT", "dock6"),
    "LINGO3DMOL": ("SBDD_EXTERNAL_ROOT", "Lingo3DMol"),
    "TAMGEN": ("SBDD_EXTERNAL_ROOT", "TamGen"),
    "MOLCRAFT": ("SBDD_SOURCE_ROOT", "MolCRAFT/MolCRAFT"),
    "RTMSCORE": ("SBDD_EXTERNAL_ROOT", "RTMScore"),
    "CHEMFH": ("SBDD_EXTERNAL_ROOT", "ChemFH"),
}
MEMORY_BUNDLES = ("drugmate-generation-memory", "drugmate-docking-memory")
_TEMPLATE = re.compile(r"^\$\{([A-Z_][A-Z_0-9]*)\}(?:/(.*))?$")
_SYMBOLIC_PATHS = ContextVar("deployment_symbolic_paths", default=False)


@contextmanager
def symbolic_deployment_paths():
    """Build public registry metadata without reading process placement settings."""
    token = _SYMBOLIC_PATHS.set(True)
    try:
        yield
    finally:
        _SYMBOLIC_PATHS.reset(token)


def _absolute_env(name: str) -> Path | None:
    value = os.environ.get(name)
    if value is None:
        return None
    path = Path(value).expanduser()
    if not value.strip() or not path.is_absolute():
        raise ValueError(f"{name} must be a nonempty absolute path")
    return path


def deployment_root(variable: str) -> Path:
    """Resolve a supported deployment variable with a portable default."""
    source = next((location for tool, location in SOURCE_LOCATIONS.items()
                   if variable == f"SBDD_{tool}_SOURCE_DIR"), None)
    base_variables = {
        "BIO_SKILLS_ROOT", "SBDD_DATA_ROOT", "SBDD_SOURCE_ROOT", "SBDD_EXTERNAL_ROOT",
        "SBDD_MEMORY_ROOT", "SBDD_MEMORY_BUNDLES_ROOT", "SBDD_BENCHMARK_ROOT",
        "SBDD_CONDA_ENVS_ROOT", "SBDD_AIZYNTHFINDER_CONFIG",
    }
    if source is None and variable not in base_variables:
        raise ValueError(f"unsupported deployment variable: {variable}")
    override = _absolute_env(variable)
    if override is not None:
        return override
    if source is not None:
        root, relative = source
        return deployment_root(root) / relative
    if variable == "BIO_SKILLS_ROOT":
        return Path(__file__).resolve().parent.parent
    if variable == "SBDD_MEMORY_ROOT":
        return Path(__file__).resolve().parents[3] / "memory"
    if variable == "SBDD_MEMORY_BUNDLES_ROOT":
        return deployment_root("SBDD_MEMORY_ROOT") / "bundles"
    if variable == "SBDD_AIZYNTHFINDER_CONFIG":
        return deployment_root("SBDD_SOURCE_ROOT") / "models" / "aizynthfinder" / "config.yml"
    relative = {
        "SBDD_DATA_ROOT": "drugmate-runs",
        "SBDD_SOURCE_ROOT": "drugmate-tools",
        "SBDD_EXTERNAL_ROOT": "drugmate-tools",
        "SBDD_BENCHMARK_ROOT": "drugmate-benchmarks",
        "SBDD_CONDA_ENVS_ROOT": "drugmate-envs",
    }
    return Path.home() / relative[variable]


def deployment_path(default: str | Path) -> Path:
    """Resolve a deployment template; explicit caller paths remain unchanged."""
    raw = str(default)
    match = _TEMPLATE.fullmatch(raw)
    if match:
        variable, suffix = match.groups()
        relative = Path(suffix or "")
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("deployment template suffix must remain within its root")
        if _SYMBOLIC_PATHS.get():
            return Path(raw)
        return deployment_root(variable) / relative
    if "${" in raw:
        raise ValueError("deployment paths require one leading ${VARIABLE} template")
    return Path(default)


def data_root() -> Path:
    return deployment_root("SBDD_DATA_ROOT")


def configure_environments(environments: dict[str, dict]) -> dict[str, dict]:
    """Apply placement overrides without replacing environment/tool identities."""
    if _SYMBOLIC_PATHS.get():
        return {env_id: dict(spec) for env_id, spec in environments.items()}
    configured = {}
    common_root = _absolute_env("SBDD_CONDA_ENVS_ROOT")
    for env_id, original in environments.items():
        spec = dict(original)
        if spec.get("kind") == "conda":
            variable = f"SBDD_ENV_{env_id.upper()}"
            prefix = _absolute_env(variable + "_PREFIX")
            name = os.environ.get(variable + "_NAME")
            if prefix is not None and name is not None:
                raise ValueError(f"set only one of {variable}_PREFIX and {variable}_NAME")
            if name is not None:
                if not name.strip() or name != name.strip() or any(c in name for c in "/\\"):
                    raise ValueError(f"{variable}_NAME must be a Conda environment name")
                spec.pop("prefix", None)
                spec["name"] = name
            elif prefix is not None or common_root is not None:
                if prefix is None:
                    directory = original.get("name") or Path(original["prefix"]).name
                    prefix = common_root / directory
                spec.pop("name", None)
                spec["prefix"] = str(prefix)
        configured[env_id] = spec
    return configured


def describe() -> dict:
    """Report placement data without process credentials or availability claims."""
    import pipeline_registry as registry

    return {
        "data_root": str(data_root()),
        "sources": {tool.lower(): str(deployment_root(f"SBDD_{tool}_SOURCE_DIR"))
                    for tool in SOURCE_LOCATIONS},
        "aizynthfinder_config": str(deployment_root("SBDD_AIZYNTHFINDER_CONFIG")),
        "benchmark_root": str(deployment_root("SBDD_BENCHMARK_ROOT")),
        "memory_root": str(deployment_root("SBDD_MEMORY_ROOT")),
        "memory_bundles": {name: str(deployment_root("SBDD_MEMORY_BUNDLES_ROOT") / name)
                           for name in MEMORY_BUNDLES},
        "tool_environments": registry.TOOL_ENVIRONMENTS,
        "availability_checked": False,
    }


if __name__ == "__main__":
    print(json.dumps(describe(), indent=2, ensure_ascii=False))

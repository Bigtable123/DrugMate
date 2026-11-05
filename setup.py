from __future__ import annotations

from collections import defaultdict
import os
from pathlib import Path

from setuptools import setup


REPO_ROOT = Path(__file__).parent.resolve()


def _data_file_tree(root_name: str) -> list[tuple[str, list[str]]]:
    root = REPO_ROOT / root_name
    grouped: defaultdict[str, list[str]] = defaultdict(list)
    excluded_dirs = {"__pycache__", "node_modules", ".venv", ".pytest_cache", ".local"}
    for directory, dirs, filenames in os.walk(root):
        current = Path(directory)
        # Prune local downloads before walking them: memory bundles can contain
        # gigabytes of data and are distributed separately from Python packages.
        dirs[:] = sorted(name for name in dirs
                         if name not in excluded_dirs
                         and not (root_name == "memory" and (
                             (current == root and name == "bundles")
                             or name in {"drugmate-docking-memory", "drugmate-generation-memory"}
                         )))
        for filename in sorted(filenames):
            path = current / filename
            if not path.is_file():
                continue
            if (path.name == ".DS_Store"
                    or path.suffix in {".pyc", ".pyo", ".log", ".db", ".sqlite", ".sqlite3"}
                    or path.name.startswith((
                        "drugmate-generation-memory.tar.gz", "drugmate-docking-memory.tar.gz",
                        "drugmate-generation-memory.zip", "drugmate-docking-memory.zip",
                    ))
                    or path.name.startswith("deployment.local.")
                    or path.name == ".env"
                    or (path.name.startswith(".env.") and path.name != ".env.example")):
                continue
            rel_path = path.relative_to(REPO_ROOT)
            grouped[str(rel_path.parent)].append(str(rel_path))
    return sorted(grouped.items())


setup(
    data_files=[
        ("share/drugmate", ["SOUL.md"]),
        ("share/drugmate/licenses/viewers", [
            "licenses/viewers/3Dmol-min.js.LICENSE.txt",
            "licenses/viewers/molstar.LICENSE.txt",
        ]),
        # Keep all data_files together: a separate pyproject.toml data-files
        # table can replace these dynamically collected scientific assets.
        *_data_file_tree("locales"),
        *_data_file_tree("skills"),
        *_data_file_tree("memory"),
    ]
)

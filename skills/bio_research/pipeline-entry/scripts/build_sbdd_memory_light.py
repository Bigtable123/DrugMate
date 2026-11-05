"""Build, verify and install a planning sidecar without rewriting a frozen release."""
from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import statistics
import sys
import tempfile

from sbdd_memory_light import (
    DEFAULT_BUNDLE, FORMAT_VERSION, LightweightMemory, close_reader,
    digest, load_source_reader, planning_view, project_card, serialize, source_binding,
)


LAUNCHER = '''"""Read lightweight planning cards; the frozen reader remains read_memory.py."""
from pathlib import Path
import runpy
if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parent / "planning/reader.py"), run_name="__main__")
'''
GUIDE = '''# SBDD Lightweight Planning Memory

Each source task context has one lightweight card in cards.jsonl; index.json
locates it by context_id. manifest.json binds the cards, index and reader to
the frozen source bundle and runtime supplement through SHA-256 hashes.
Cards omit pose coordinates, embedding arrays, full parameters, logs and
protein sequences.

Use these cards for pocket_prediction or docking_memory_mode=legacy recall,
with benchmark-view=planning by default. Docking modes ligand_only and
ligand_protein use select-docking; off reads no docking history.
read_memory.py retains its existing behavior; read_memory_light.py returns
lightweight cards by default. Run from the installed memory bundle root, with
SBDD_PYTHON pointing to an environment compatible with its reader:

```bash
"${SBDD_PYTHON}" read_memory_light.py show --task docking --dataset astex_diverse --target 1G9V_RQ3
"${SBDD_PYTHON}" read_memory_light.py neighbors --task pocket_prediction --dataset astex_diverse --target 1G9V_RQ3 -k 3
```

Stored card metrics align with metric_columns; each cell contains
[value, status, optional reason]. Command output also uses a column-oriented
tool table defined by tool_columns. Its values align with metric_columns.
In states, '*' sets the row's default status and numeric indexes identify
exceptional statuses/reasons; unspecified values use common.metric_status.
run contains [execution_status, parsing_status, return_code]. Null run, flags
or signals values inherit the corresponding common value; protocol and cost
merge the corresponding common defaults. Timing values align with time_columns;
scope defaults to common.scope and extra may override it.

false, 0 and null have distinct meanings. Exclusion of a quality cohort does
not mean every tool failed. RMSD, DCA and DCC use angstroms; timing uses seconds.
Other metrics are booleans, counts or ratios. signals contains recorded
observations, not recommendations or success rates for the current target.
Cases sharing a group_id are not independent replicates. Top-K RMSD does not
include Top-K PoseBusters validation. PocketMiner spatial metrics are
not_applicable; execution failure, zero pockets and missing metrics remain
distinct. cost_scope identifies the timing boundary: local durations from
different stages or tools must not be summed into an end-to-end budget.
additional_observation_count indicates further recorded attempts; inspect
their details when they matter to a decision.

Use the card's details to locate read_memory.py show --view execution --tool
TOOL_ID for exact parameters. Read source evidence through its evidence command
using an artifact ID or selector. Similarity retrieval uses
DecisionMemory.neighbors, preserving encoder compatibility, current-case
exclusion and bundle access boundaries. Retrieval has no minimum similarity
threshold; the Agent must check the current structure, site and experimental
protocol before adopting nearest-neighbor evidence. Lightweight cards are
read-only: they neither supply current-query features nor modify the fixed split.
'''


def build_planning(bundle: Path, output: Path) -> dict:
    module, reader = load_source_reader(bundle)
    output = output.expanduser().resolve()
    try:
        if output == bundle.resolve() or bundle.resolve() in output.parents:
            raise ValueError("Build into a staging directory; use install to add it to the source bundle")
        output.mkdir(parents=True, exist_ok=False)
        index, task_counts, tool_counts = {}, Counter(), Counter()
        light_lengths, display_lengths, source_lengths, compare_lengths = [], [], [], []
        with (output / "cards.jsonl").open("wb") as handle:
            for row in reader.index.execute("SELECT * FROM cards ORDER BY context_id"):
                card = reader.load(row)  # Includes checked runtime supplement.
                light = project_card(card)
                payload = serialize(light).encode("utf-8")
                cid = card["identity"]["context_id"]
                index[cid] = [handle.tell(), len(payload)]
                handle.write(payload)
                task_counts[card["task_context"]["task_type"]] += 1
                tool_counts.update(tool["tool"] for tool in light["tools"])
                light_lengths.append(len(serialize(light)))
                display_lengths.append(len(serialize(planning_view(light))))
                source_lengths.append(row["length"])
                compare_lengths.append(len(serialize(module.comparison_view(card))))
        (output / "index.json").write_text(serialize(index), encoding="utf-8")
        shutil.copy2(Path(__file__).with_name("sbdd_memory_light.py"), output / "reader.py")
        (output / "launcher.py").write_text(LAUNCHER, encoding="utf-8")
        (output / "README.md").write_text(GUIDE, encoding="utf-8")
        names = ["cards.jsonl", "index.json", "reader.py", "launcher.py", "README.md"]
        manifest = {
            "schema_version": FORMAT_VERSION, "status": "passed",
            "source": source_binding(reader),
            "counts": {"contexts": len(index), "tasks": dict(task_counts), "tools": dict(tool_counts)},
            "artifacts": {name: {"sha256": digest(output / name),
                                 "bytes": (output / name).stat().st_size} for name in names},
            "generator": {"script": Path(__file__).name, "sha256": digest(Path(__file__)),
                          "projection_sha256": digest(Path(__file__).with_name("sbdd_memory_light.py"))},
        }
        (output / "manifest.json").write_text(serialize(manifest), encoding="utf-8")
        report = verify_planning(bundle, output, reader=reader)
        report["sizes"] = {
            "source_jsonl_bytes": sum(source_lengths),
            "light_jsonl_bytes": (output / "cards.jsonl").stat().st_size,
            "light_card_chars": {"median": statistics.median(light_lengths),
                                 "min": min(light_lengths), "max": max(light_lengths)},
            "planning_view_chars": {"median": statistics.median(display_lengths),
                                    "min": min(display_lengths), "max": max(display_lengths)},
            "source_compare_compact_median_chars": statistics.median(compare_lengths),
            "source_raw_median_bytes": statistics.median(source_lengths),
            "jsonl_reduction_percent": round(100 * (1 - (output / "cards.jsonl").stat().st_size
                                                       / sum(source_lengths)), 3),
        }
        (output / "validation_report.json").write_text(serialize(report), encoding="utf-8")
        return report
    finally:
        close_reader(reader)


def verify_planning(bundle: Path, directory: Path, reader=None) -> dict:
    owns_reader = reader is None
    if owns_reader:
        _, reader = load_source_reader(bundle)
    try:
        store = LightweightMemory(bundle, reader, directory)
        for name, expected in store.manifest["artifacts"].items():
            # Only the fixed files produced by this builder may be verified.
            if name not in {"cards.jsonl", "index.json", "reader.py", "launcher.py", "README.md", "README_zh.md"}:
                raise RuntimeError("Unexpected sidecar artifact")
            if not (directory / name).is_file() or digest(directory / name) != expected["sha256"]:
                raise RuntimeError("Sidecar artifact hash mismatch: " + name)
        checked, counts = 0, Counter()
        for row in reader.index.execute("SELECT * FROM cards ORDER BY context_id"):
            expected = project_card(reader.load(row))
            actual = store.get(row["context_id"])
            # Location expansion is a reader concern; stored facts remain byte-derived.
            actual["details"]["script"] = "read_memory.py"
            actual["evidence_record"].pop("path", None)
            actual["evidence_record"].pop("sha256", None)
            if actual != expected:
                raise RuntimeError("Planning projection differs from source: " + row["context_id"])
            checked += 1
            counts[row["task_type"]] += 1
        return {"status": "passed", "checked_contexts": checked, "tasks": dict(counts),
                "checks": ["source release and supplement hashes", "exact source context coverage",
                           "all card projections match source", "all sidecar artifact hashes"],
                "source_bundle": str(bundle.resolve()), "planning_directory": str(directory.resolve())}
    finally:
        if owns_reader:
            close_reader(reader)


def install_planning(bundle: Path, staged: Path) -> dict:
    report = verify_planning(bundle, staged)
    destination = bundle / "planning"
    launcher = bundle / "read_memory_light.py"
    launcher_bytes = (staged / "launcher.py").read_bytes()
    if launcher.exists() and launcher.read_bytes() != launcher_bytes:
        raise FileExistsError("Refusing to overwrite an existing read_memory_light.py")
    if destination.exists():
        old = json.loads((destination / "manifest.json").read_text(encoding="utf-8"))
        new = json.loads((staged / "manifest.json").read_text(encoding="utf-8"))
        if old != new:
            raise FileExistsError("Refusing to overwrite a different planning release")
        verify_planning(bundle, destination)
    else:
        temp = Path(tempfile.mkdtemp(prefix=".planning-install-", dir=bundle))
        shutil.copytree(staged, temp, dirs_exist_ok=True)
        os.replace(temp, destination)
    if not launcher.exists():
        with launcher.open("xb") as handle:
            handle.write(launcher_bytes)
    report["installed"] = {"planning": str(destination), "reader": str(launcher)}
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("build").add_argument("--output", required=True, type=Path)
    sub.add_parser("verify").add_argument("--planning-dir", required=True, type=Path)
    sub.add_parser("install").add_argument("--staged", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            result = build_planning(args.bundle, args.output)
        elif args.command == "verify":
            result = verify_planning(args.bundle, args.planning_dir)
        else:
            result = install_planning(args.bundle, args.staged)
        print(serialize(result), end="")
        return 0
    except (OSError, ValueError, RuntimeError, KeyError, ImportError) as exc:
        print(serialize({"status": "error", "error": str(exc)}), end="")
        return 1


if __name__ == "__main__":
    sys.exit(main())

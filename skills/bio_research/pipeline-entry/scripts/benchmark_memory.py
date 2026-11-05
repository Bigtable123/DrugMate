"""Query frozen selectors and append/recall evaluated SBDD research memory."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3
import sys
import zipfile


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from deployment_config import deployment_path


DOCKING_MODES = ("off", "legacy", "ligand_only", "ligand_protein")
DEFAULT_SELECTOR_INDEX = deployment_path("${SBDD_MEMORY_ROOT}/docking_memory/ligand_protein_v2.json")
DEFAULT_GENERATION_INDEX = deployment_path("${SBDD_MEMORY_ROOT}/generation_memory/ligand_protein_v1.json")
DEFAULT_POCKET_INDEX = deployment_path("${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json")


class QueryError(ValueError):
    def __init__(self, message: str, code: str = "invalid_arguments"):
        super().__init__(message)
        self.code = code


class JsonArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise QueryError(message)


def load_views():
    return load_sibling("sbdd_memory_recall")


def load_sibling(name):
    path = Path(__file__).with_name(name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    old_bytecode = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = old_bytecode
    return module


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bind_query_metadata(encoding_dir: Path, receptor_file: Path, output: Path) -> dict:
    """Verify encoder completion/provenance and write a separate query binding.

    Supports the actual run_batch and run_target_assets output contracts, without
    importing either scientific encoder or reconstructing provenance from names.
    """
    directory = Path(encoding_dir).resolve()
    receptor, output = Path(receptor_file).resolve(), Path(output).resolve()
    metadata_path, prepared_path, complete_path = (
        directory / name for name in ("metadata.json", "prepared.json", "complete.json")
    )
    metadata, prepared, complete = (
        json.loads(path.read_text(encoding="utf-8")) for path in (metadata_path, prepared_path, complete_path)
    )
    if any(not isinstance(value, dict) for value in (metadata, prepared, complete)):
        raise ValueError("encoder metadata, preparation and completion must be objects")

    def fingerprint(value):
        return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                         allow_nan=False).encode()).hexdigest()

    recipe, encoder, source = prepared.get("preparation_recipe"), metadata.get("encoder"), metadata.get("source_record")
    if not isinstance(recipe, dict) or not isinstance(encoder, dict) or not encoder or not isinstance(source, dict):
        raise ValueError("missing original preparation recipe, encoder or source record")
    if source != prepared.get("source_record") or source != recipe.get("source_record"):
        raise ValueError("encoder metadata and preparation source records disagree")
    task = source.get("task_type")
    if task not in ("docking", "pocket_prediction", "molecule_generation") or recipe.get("task_type") != task:
        raise ValueError("unsupported or inconsistent original query task type")
    prepared_fingerprint = fingerprint(recipe)
    if prepared.get("preparation_fingerprint") != prepared_fingerprint:
        raise ValueError("preparation fingerprint does not match the original recipe")
    expected_files = {"metadata.json", "embeddings.npz"}
    if "protein_preparation_fingerprint" in metadata:
        if task != "molecule_generation" or metadata["protein_preparation_fingerprint"] != prepared_fingerprint:
            raise ValueError("generation protein preparation fingerprint mismatch")
        pocket = json.loads((directory / "pocket_prepared.json").read_text(encoding="utf-8"))
        pocket_recipe = pocket.get("recipe")
        if not isinstance(pocket_recipe, dict) or pocket_recipe.get("core_preparation_fingerprint") != prepared_fingerprint:
            raise ValueError("generation pocket preparation is not bound to the original protein")
        pocket_fingerprint = fingerprint(pocket_recipe)
        if pocket.get("preparation_fingerprint") != pocket_fingerprint or metadata.get("pocket_preparation_fingerprint") != pocket_fingerprint:
            raise ValueError("generation pocket preparation fingerprint mismatch")
        pocket_status = pocket.get("status")
        if pocket_status not in ("available", "missing_cavity_asset") or metadata.get("pocket_status") != pocket_status:
            raise ValueError("inconsistent generation pocket status")
        expected_files |= {"prepared.json", "pocket_prepared.json"}
        if pocket_status == "available":
            expected_files |= {"cavity_points.npz", "pocket_embeddings.npz"}
        encoding_fingerprint = fingerprint({"protein": prepared_fingerprint, "pocket": pocket_fingerprint, "encoder": encoder})
        encoder_format = "run_target_assets"
    else:
        if metadata.get("preparation_fingerprint") != prepared_fingerprint:
            raise ValueError("encoder metadata preparation fingerprint mismatch")
        encoding_fingerprint = fingerprint({"preparation_fingerprint": prepared_fingerprint, "encoder": encoder})
        encoder_format = "run_batch"
    if metadata.get("encoding_fingerprint") != encoding_fingerprint or complete.get("encoding_fingerprint") != encoding_fingerprint:
        raise ValueError("encoding fingerprint does not match preparation and encoder")
    completed_hashes = complete.get("output_sha256")
    if not isinstance(completed_hashes, dict) or set(completed_hashes) != expected_files:
        raise ValueError("encoder completion has missing or unexpected output hashes")
    for name in sorted(expected_files):
        if file_digest(directory / name) != completed_hashes[name]:
            raise ValueError(f"encoder completion SHA-256 mismatch for {name}")
    receptor_hash = file_digest(receptor)
    if receptor_hash != recipe.get("protein_sha256") or (source.get("protein_sha256") and source["protein_sha256"] != receptor_hash):
        raise ValueError("current receptor does not match the original encoder preparation input")
    if output in {receptor, *(path.resolve() for path in (directory / name for name in
                                                        expected_files | {"prepared.json", "complete.json"}))}:
        raise ValueError("write bound metadata to a separate output; never overwrite encoder inputs")
    binding = {"encoder_format": encoder_format, "encoding_directory": str(directory),
               "raw_metadata_sha256": file_digest(metadata_path), "prepared_sha256": file_digest(prepared_path),
               "complete_sha256": file_digest(complete_path), "preparation_fingerprint": prepared_fingerprint,
               "encoding_fingerprint": encoding_fingerprint, "receptor_file": str(receptor)}
    result = {**metadata, "embedding_sha256": completed_hashes["embeddings.npz"],
              "source_input_sha256": receptor_hash, "query_binding": binding}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    return {"success": True, "status": "bound", "output": str(output), "task_type": task,
            "embedding_sha256": result["embedding_sha256"], "source_input_sha256": receptor_hash,
            "metadata_sha256": file_digest(output), "provenance": binding}


def recall(args, views) -> dict:
    if not 1 <= args.k <= 10 or not 4000 <= args.max_chars <= 60000:
        raise QueryError("recall requires 1 <= -k <= 10 and 4000 <= --max-chars <= 60000")
    if bool(args.query_npz) != bool(args.query_metadata):
        raise QueryError("Supply both --query-npz and --query-metadata")
    if bool(args.dataset) != bool(args.target):
        raise QueryError("Supply both --dataset and --target")
    if args.query_npz and args.dataset:
        raise QueryError("Use query NPZ/metadata or a bundled dataset/target, not both")
    args.benchmark_bundle = args.benchmark_bundle or views.DEFAULT_BENCHMARK
    requested = bool(args.query_npz or args.dataset)
    result = {
        "success": True, "schema_version": "sbdd_recall/1", "cards": [],
        "sources": {"benchmark": {
            "status": "not_requested", "reason": "No embedding or bundled query supplied",
        }},
        "errors": [], "summary_incomplete": False,
        "usage": (
            "Historical transfer is tentative; cosine is not success probability. "
            "Check current structure, site and protocol before adopting neighbor evidence. "
            "Fetch details for exact parameters. Benchmark stays frozen."
        ),
    }
    if args.benchmark_view == "summary":
        result["metric_columns"] = ["value", "status", "unit", "optional reason"]
    if requested:
        try:
            cards, status = views.benchmark_recall(args, file_digest)
            result["cards"] = cards
            result["sources"]["benchmark"] = status
            if args.benchmark_view == "planning":
                result["usage"] = views.load_planning_module().PLANNING_USAGE
        except (OSError, ValueError, RuntimeError, ImportError, KeyError, TypeError,
                AttributeError, RecursionError, sqlite3.Error, zipfile.BadZipFile) as exc:
            result["success"] = False
            result["sources"]["benchmark"] = {
                "status": "error", "bundle": str(args.benchmark_bundle), "reason": str(exc),
            }
            result["errors"].append({"source": "benchmark", "error": str(exc)})
    return views.bound_output(result, args.max_chars)


def select_docking(args, views) -> dict:
    """Planning only: the Agent applies the result to its explicit tool matrix."""
    selector = load_sibling("docking_selector")
    task = getattr(args, "task", "redocking")
    policy = selector.task_policy(task, getattr(args, "metric", None))
    metric, baseline = policy["metric"], policy["default_tool"]
    if metric not in selector.METRICS:
        raise QueryError("Unknown docking selection metric")
    configured_mode = None
    if args.config:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        if not isinstance(config, dict) or not isinstance(config.get("context", {}), dict):
            raise QueryError("--config must contain a pipeline object with an object context")
        configured_mode = config.get("context", {}).get("docking_memory_mode", "ligand_protein")
        if configured_mode not in DOCKING_MODES:
            raise QueryError("context.docking_memory_mode must be one of " + ", ".join(DOCKING_MODES))
        if args.mode is not None and args.mode != configured_mode:
            raise QueryError("--mode conflicts with context.docking_memory_mode; use one experiment mode")
    mode = args.mode or configured_mode or "ligand_protein"
    # Registration is the default capability boundary, not a runtime probe.
    # The planner narrows this list to installed engines compatible with the inputs.
    registry_path = Path(__file__).resolve().parents[2] / "shared"
    sys.path.insert(0, str(registry_path))
    from pipeline_registry import DOCKING_PRIMARY_ENGINES
    available = list(DOCKING_PRIMARY_ENGINES) if args.available_tools is None else args.available_tools
    unknown = set(available) - set(DOCKING_PRIMARY_ENGINES)
    if unknown:
        raise QueryError("Unregistered docking engines in --available-tools: " + ", ".join(sorted(unknown)))
    if mode == "off":
        # No query, selector index, historical cards, or scientific dependency is opened.
        result = selector.fallback_decision("memory_disabled", available, mode, metric, task)
        result.update({
            "schema_version": "sbdd_docking_selection/1", "status": "disabled",
            "reason": f"Docking memory disabled; fixed {baseline} baseline" if baseline in available
                      else f"Docking memory disabled; {baseline} unavailable",
            "available_tools": available,
        })
        return result
    if mode == "legacy":
        if task != "redocking" or metric != "top1_rmsd_lt_2a":
            raise QueryError("legacy supports redocking evidence cards only, not pose-generation or metric-specific selection")
        if args.exclude_group or args.query_group:
            raise QueryError("legacy recall does not support group exclusion; supply --exclude-case for every group member")
        legacy_args = argparse.Namespace(**vars(args))
        legacy_args.task = "docking"
        if args.query_case:
            legacy_args.exclude_case = [*args.exclude_case, args.query_case]
        result = recall(legacy_args, views)
        result.update({"mode": mode, "selected_tool": None, "task": task, "task_policy": policy, "metric": metric})
        return result

    try:
        if args.dataset or args.target:
            raise ValueError("New docking memory requires current --ligand-sdf and query features; dataset/target is legacy-only")
        if not args.ligand_sdf:
            raise ValueError("Supply --ligand-sdf for one prepared query ligand (never a reference bound pose)")
        index = selector.load_index(args.selector_index)
        excluded_cases = set(args.exclude_case)
        excluded_groups = set(args.exclude_group)
        if args.query_case:
            excluded_cases.add(args.query_case)
        if args.query_group:
            excluded_groups.add(args.query_group)
        query = {"ligand_sdf": str(args.ligand_sdf.resolve()),
                 "ligand_sha256": file_digest(args.ligand_sdf)}
        if mode == "ligand_protein":
            if not args.query_npz or not args.query_metadata:
                raise ValueError("ligand_protein requires --query-npz and --query-metadata")
            binding = selector.validate_query_metadata(
                args.query_metadata, args.query_npz, index, receptor_file=args.receptor_file,
            )
        elif args.query_metadata:
            # Identity is shared across ablations; ligand-only never opens NPZ
            # or checks a protein encoder, but must exclude the same query group.
            binding = selector.query_identity(args.query_metadata)
        else:
            binding = {}
        query.update(binding)
        if binding.get("case_id"):
            excluded_cases.add(binding["case_id"])
        if binding.get("group_id"):
            excluded_groups.add(binding["group_id"])
        features = selector.extract_query_features(
            args.ligand_sdf, args.query_npz if mode == "ligand_protein" else None,
        )
        decision = selector.select_docking_tool(
            features, index, mode=mode, available_tools=available,
            exclude_cases=excluded_cases, exclude_groups=excluded_groups,
            metric=metric, task=task,
        )
        decision["query"] = query
        decision["index"] = {"path": str(args.selector_index.resolve()),
                             "sha256": file_digest(args.selector_index)}
        return decision
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, AttributeError,
            IndexError, OverflowError, zipfile.BadZipFile) as exc:
        return selector.fallback_decision(str(exc), available_tools=available, mode=mode,
                                          metric=metric, task=task)


def select_generation(args) -> dict:
    selector = load_sibling("generation_selector")
    if set(args.priority_metric) - set(selector.METRICS):
        raise QueryError("Unknown generation priority metric; choose from " + ", ".join(selector.METRICS))
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
    from pipeline_registry import STEPS
    registered = set(STEPS["molecule_generation"]["tools"])
    supported = registered & set(selector.TOOLS)
    available = sorted(supported) if args.available_tools is None else args.available_tools
    if set(available) - supported:
        raise QueryError("--available-tools must contain registered generators with memory labels")
    try:
        if not args.ligand_sdf or not args.query_npz or not args.query_metadata:
            raise ValueError("generation memory requires one --ligand-sdf, --query-npz and --query-metadata")
        index = selector.load_index(args.selector_index)
        binding = selector.validate_query_metadata(
            args.query_metadata, args.query_npz, index,
            receptor_file=args.receptor_file, ligand_sdf=args.ligand_sdf,
        )
        excluded_cases, excluded_groups = set(args.exclude_case), set(args.exclude_group)
        excluded_cases.update(value for value in (args.query_case, binding.get("case_id")) if value)
        excluded_groups.update(value for value in (args.query_group, binding.get("group_id")) if value)
        result = selector.select_generation_tools(
            selector.extract_query_features(args.ligand_sdf, args.query_npz), index,
            priority_metrics=args.priority_metric, available_tools=available,
            exclude_cases=excluded_cases, exclude_groups=excluded_groups,
        )
        result["query"] = {**binding, "ligand_sdf": str(args.ligand_sdf.resolve()),
                           "ligand_sha256": file_digest(args.ligand_sdf)}
        result["index"] = {"path": str(args.selector_index.resolve()),
                           "sha256": file_digest(args.selector_index)}
    except (OSError, ValueError, RuntimeError, ImportError, KeyError, TypeError,
            AttributeError, IndexError, OverflowError, zipfile.BadZipFile) as exc:
        result = selector.fallback_decision(str(exc), available, args.priority_metric)
    result["available_tools"] = available
    return result


def select_pocket(args) -> dict:
    selector = load_sibling("pocket_memory")
    return selector.select_pocket(
        index_path=args.selector_index,
        query_npz=args.query_npz, query_metadata=args.query_metadata,
        receptor_file=args.receptor_file, candidates_path=args.candidates,
        memory_mode=args.memory_mode,
        preset=args.preset, query_case=args.query_case, query_group=args.query_group,
        exclude_cases=args.exclude_case, exclude_groups=args.exclude_group,
        export_box=args.export_box, box_padding=args.box_padding,
        min_box_size=args.min_box_size,
    )


def main(argv: list[str] | None = None) -> int:
    views = load_views()
    parser = JsonArgumentParser(description=__doc__, allow_abbrev=False)
    commands = parser.add_subparsers(dest="command", required=True)
    append_research = commands.add_parser("append-research-run", help="Append hash-bound evaluated observations to separate research memory", allow_abbrev=False)
    append_research.add_argument("--record", type=Path, required=True)
    append_research.add_argument("--store", type=Path)
    append_research.add_argument("--mode", choices=("research", "frozen"), default="research")
    research = commands.add_parser("recall-research", help="Retrieve prior research observations for planner review without changing frozen selectors", allow_abbrev=False)
    research.add_argument("--task", required=True, choices=("docking", "pocket_prediction", "molecule_generation"))
    research.add_argument("--receptor-file", type=Path, required=True)
    research.add_argument("--query-npz", type=Path, required=True)
    research.add_argument("--query-metadata", type=Path, required=True)
    research.add_argument("--ligand-sdf", type=Path)
    research.add_argument("--store", type=Path)
    research.add_argument("--mode", choices=("research", "frozen"), default="research")
    research.add_argument("-k", type=int, help="Default: 3 ligand profiles or 5 protein groups; overrides are diagnostic")
    research.add_argument("--exclude-case", action="append", default=[])
    research.add_argument("--exclude-group", action="append", default=[])
    research.add_argument("--exclude-run", action="append", default=[])
    bind = commands.add_parser("bind-query-metadata", help="Verify raw encoder provenance and write bound query metadata",
                               allow_abbrev=False)
    bind.add_argument("--encoding-dir", type=Path, required=True)
    bind.add_argument("--receptor-file", type=Path, required=True)
    bind.add_argument("--output", type=Path, required=True)
    query = commands.add_parser("recall", help="Read lightweight benchmark planning cards",
                                allow_abbrev=False)
    query.add_argument("--task", required=True, choices=("docking", "pocket_prediction"))
    query.add_argument("--receptor-file", type=Path)
    query.add_argument("--query-npz", type=Path)
    query.add_argument("--query-metadata", type=Path)
    query.add_argument("--dataset")
    query.add_argument("--target")
    query.add_argument("--benchmark-bundle", type=Path, help="Default: the fixed split memory release")
    query.add_argument("--benchmark-view", choices=("planning", "summary"), default="planning")
    query.add_argument("--exclude-case", action="append", default=[])
    query.add_argument("-k", type=int, default=3)
    query.add_argument("--max-chars", type=int, default=48000)
    select = commands.add_parser("select-docking", help="Deterministic docking memory with ablation modes",
                                 allow_abbrev=False)
    select.add_argument("--config", type=Path, help="Read context.docking_memory_mode from this pipeline JSON")
    select.add_argument("--mode", choices=DOCKING_MODES)
    select.add_argument("--task", choices=("redocking", "pose_generation"), default="redocking",
                        help="Planning task: redocking defaults to GNINA; pose_generation defaults to rDock")
    select.add_argument("--metric", default=None, choices=(
        "top1_rmsd_lt_2a", "posebusters_valid", "joint_success", "top1_rmsd_a", "median_rmsd_a",
    ), help="Default: Top-1 RMSD <2 A for redocking; Top-1 PoseBusters validity for pose_generation")
    select.add_argument("--ligand-sdf", type=Path)
    select.add_argument("--query-npz", type=Path)
    select.add_argument("--query-metadata", type=Path)
    select.add_argument("--receptor-file", type=Path)
    select.add_argument("--selector-index", type=Path, default=DEFAULT_SELECTOR_INDEX)
    select.add_argument("--available-tools", nargs="*", help="Registered engines verified usable for this query; empty means none")
    select.add_argument("--query-case")
    select.add_argument("--query-group")
    select.add_argument("--exclude-case", action="append", default=[])
    select.add_argument("--exclude-group", action="append", default=[])
    select.add_argument("--dataset", help="Legacy recall only")
    select.add_argument("--target", help="Legacy recall only")
    select.add_argument("--benchmark-bundle", type=Path)
    select.add_argument("--benchmark-view", choices=("planning", "summary"), default="planning")
    select.add_argument("-k", type=int, default=3, help="Legacy recall only; new selector always uses k=3")
    select.add_argument("--max-chars", type=int, default=48000, help="Legacy recall output limit")
    build = commands.add_parser("build-docking-index", help="Build a derived index from memory only",
                                allow_abbrev=False)
    build.add_argument("--benchmark-bundle", type=Path, default=views.DEFAULT_BENCHMARK)
    build.add_argument("--output", type=Path, default=deployment_path("${SBDD_DATA_ROOT}/docking_memory/ligand_protein_v2.json"))
    generation = commands.add_parser("select-generation", help="Five independent generation metric recommendations",
                                     allow_abbrev=False)
    generation.add_argument("--ligand-sdf", type=Path)
    generation.add_argument("--query-npz", type=Path)
    generation.add_argument("--query-metadata", type=Path)
    generation.add_argument("--receptor-file", type=Path)
    generation.add_argument("--selector-index", type=Path, default=DEFAULT_GENERATION_INDEX)
    generation.add_argument("--priority-metric", action="append", default=[],
                            help="Repeat to require an objective's recommendation; all five retain prior strength 10 and switch margin 0.03")
    generation.add_argument("--available-tools", nargs="*", help="Verified usable generators; empty means none")
    generation.add_argument("--query-case")
    generation.add_argument("--query-group")
    generation.add_argument("--exclude-case", action="append", default=[])
    generation.add_argument("--exclude-group", action="append", default=[])
    generation_build = commands.add_parser("build-generation-index", help="Build generation index from memory-only rows",
                                           allow_abbrev=False)
    generation_build.add_argument("--benchmark-bundle", type=Path)
    generation_build.add_argument("--property-winners", type=Path)
    generation_build.add_argument("--groups", type=Path)
    generation_build.add_argument("--output", type=Path, default=deployment_path("${SBDD_DATA_ROOT}/generation_memory/ligand_protein_v1.json"))
    pocket_build = commands.add_parser("build-pocket-index", help="Build protein-only pocket reliability index from the memory split",
                                       allow_abbrev=False)
    pocket_build.add_argument("--benchmark-bundle", type=Path, default=views.DEFAULT_BENCHMARK)
    pocket_build.add_argument("--output", type=Path, default=DEFAULT_POCKET_INDEX)
    pocket_plan = commands.add_parser("plan-pocket-tools", help="Retrieve pocket detector recommendations before running detectors",
                                      allow_abbrev=False)
    pocket_plan.add_argument("--memory-mode", choices=("auto", "off", "required"), default="auto",
                             help="auto records unavailable-memory fallback; required forbids fallback; off skips memory")
    pocket_plan.add_argument("--selector-index", type=Path, default=DEFAULT_POCKET_INDEX,
                             help="Derived protein-only pocket reliability index")
    pocket_plan.add_argument("--query-npz", type=Path, help="Current receptor SaProt protein embedding")
    pocket_plan.add_argument("--query-metadata", type=Path, help="Hash-bound pocket_prediction encoder metadata")
    pocket_plan.add_argument("--receptor-file", type=Path, required=True)
    pocket_plan.add_argument("--preset", choices=("conservative", "exploratory_memory_guard"), default="conservative")
    pocket_plan.add_argument("--query-case")
    pocket_plan.add_argument("--query-group")
    pocket_plan.add_argument("--exclude-case", action="append", default=[])
    pocket_plan.add_argument("--exclude-group", action="append", default=[])
    pocket = commands.add_parser("select-pocket", help="Fuse current point5 pockets using protein memory or an explicit fixed-anchor fallback",
                                 allow_abbrev=False)
    pocket.add_argument("--memory-mode", choices=("auto", "off", "required"), default="auto",
                        help="auto falls back only when memory is unavailable; required forbids fallback; off skips memory")
    pocket.add_argument("--selector-index", type=Path, default=DEFAULT_POCKET_INDEX,
                        help="Derived protein-only pocket reliability index")
    pocket.add_argument("--query-npz", type=Path, help="Current receptor SaProt protein embedding")
    pocket.add_argument("--query-metadata", type=Path, help="Hash-bound pocket_prediction encoder metadata")
    pocket.add_argument("--receptor-file", type=Path, required=True)
    pocket.add_argument("--candidates", type=Path, required=True)
    pocket.add_argument("--preset", choices=("conservative", "exploratory_memory_guard"), default="conservative",
                        help="Historical reliability and conditional P2Rank anchor-switch policy")
    pocket.add_argument("--query-case")
    pocket.add_argument("--query-group")
    pocket.add_argument("--exclude-case", action="append", default=[])
    pocket.add_argument("--exclude-group", action="append", default=[])
    pocket.add_argument("--export-box", type=Path)
    pocket.add_argument("--box-padding", type=float, default=4.0)
    pocket.add_argument("--min-box-size", type=float, default=12.0)
    normalize = commands.add_parser("normalize-pocket-candidates", help="Map detector support points to common point5 residues",
                                    allow_abbrev=False)
    normalize.add_argument("--candidates", type=Path, required=True)
    normalize.add_argument("--receptor-file", type=Path, required=True)
    normalize.add_argument("--output", type=Path, required=True)
    try:
        args = parser.parse_args(argv)
        if args.command == "append-research-run":
            result = load_sibling("research_memory").append_run(args.record, args.store, mode=args.mode)
        elif args.command == "recall-research":
            inputs = {"receptor_file": args.receptor_file, "query_npz": args.query_npz,
                      "query_metadata": args.query_metadata}
            if args.ligand_sdf is not None:
                inputs["ligand_sdf"] = args.ligand_sdf
            result = load_sibling("research_memory").recall(args.task, inputs, args.store,
                mode=args.mode, k=args.k, exclude_cases=args.exclude_case,
                exclude_groups=args.exclude_group, exclude_runs=args.exclude_run)
        elif args.command == "bind-query-metadata":
            result = bind_query_metadata(args.encoding_dir, args.receptor_file, args.output)
        elif args.command == "select-docking":
            result = select_docking(args, views)
        elif args.command == "build-docking-index":
            result = load_sibling("docking_selector").build_index(args.benchmark_bundle, args.output)
        elif args.command == "select-generation":
            result = select_generation(args)
        elif args.command == "build-generation-index":
            selector = load_sibling("generation_selector")
            result = selector.build_index(
                args.benchmark_bundle or selector.DEFAULT_BUNDLE, args.output,
                property_winners_csv=args.property_winners or selector.DEFAULT_PROPERTY_WINNERS,
                groups_csv=args.groups or selector.DEFAULT_GROUPS,
            )
        elif args.command == "select-pocket":
            result = select_pocket(args)
        elif args.command == "plan-pocket-tools":
            result = load_sibling("pocket_memory").plan_tools(
                receptor_file=args.receptor_file, index_path=args.selector_index,
                query_npz=args.query_npz, query_metadata=args.query_metadata,
                memory_mode=args.memory_mode, preset=args.preset,
                query_case=args.query_case, query_group=args.query_group,
                exclude_cases=args.exclude_case, exclude_groups=args.exclude_group,
            )
        elif args.command == "build-pocket-index":
            result = load_sibling("pocket_index").build_index(args.benchmark_bundle, args.output)
        elif args.command == "normalize-pocket-candidates":
            normalized = load_sibling("pocket_memory").normalize_candidates(
                args.candidates, args.receptor_file, output=args.output,
            )
            result = {"success": True, "status": "normalized", "output": str(args.output.resolve()),
                      "target": normalized.get("target"),
                      "candidate_counts": {key: len(value) for key, value in normalized["detectors"].items()}}
        else:
            result = recall(args, views)
        print(views.serialize(result), end="")
        return 0 if result["success"] else 1
    except (OSError, ValueError, TypeError, RuntimeError, RecursionError, ImportError,
            KeyError, IndexError, AttributeError, OverflowError, zipfile.BadZipFile) as exc:
        code = getattr(exc, "code", "io_error" if isinstance(exc, OSError) else "invalid_input")
        print(views.serialize({"success": False, "error": {"code": code, "message": str(exc)}}), end="")
        return 1


if __name__ == "__main__":
    sys.exit(main())

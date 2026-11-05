"""Bounded planning views over the frozen SBDD benchmark.

The release's own reader verifies files and performs vector retrieval. This
adapter only selects decision-relevant fields; it never writes the release.
"""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from deployment_config import deployment_path


DEFAULT_BENCHMARK = deployment_path(
    "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-docking-memory"
)
METRICS = {
    "rmsd", "rmsd_le_2a", "posebusters_valid", "joint_success",
    "pb_failed_check_count", "pb_unavailable_check_count", "best_observed_rmsd",
    "top3_rmsd_le_2a", "top5_rmsd_le_2a", "pocket_count",
    "gold_mapping_coverage", "top1_dca_angstrom", "top1_dcc_angstrom",
    "top1_dca_le_4", "top3_dca_le_4", "top5_dca_le_4",
    "any_best_dca_angstrom", "first_hit_rank_dca_le_4", "mrr_dca_le_4",
}


def serialize(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n"


def metric_view(metric: dict) -> list:
    # Preserve false, zero, null and availability status. No missing value is
    # converted into a scientific failure or a success probability.
    values = [metric.get("value"), metric.get("status"), metric.get("unit")]
    if metric.get("reason"):
        values.append(metric["reason"])
    return values


def benchmark_summary(card: dict, similarity: float, bundle: Path, manifest: dict) -> dict:
    selected = {s["observation_id"] for s in card["published_selections"]}
    protocols = {p["observation_id"]: p["protocol_labels"]
                 for p in card["comparability"].get("per_tool_protocols", [])}
    lessons = []
    for obs in card["tool_observations"]:
        if obs["observation_id"] not in selected and obs["scope"] != "original_full_run":
            continue
        candidates = obs["outputs"]["candidates"]
        assessments = []
        for assessment in card["assessments"]:
            if assessment["observation_id"] != obs["observation_id"]:
                continue
            kind = assessment["kind"]
            is_top1 = (kind == "pose_rmsd" and candidates
                       and assessment["candidate_ids"] == [candidates[0]["candidate_id"]])
            if not is_top1 and kind not in {"top1_posebusters", "pose_set_summary", "pocket_target_summary"}:
                continue
            metrics = {name: metric_view(value) for name, value in assessment["metrics"].items()
                       if name in METRICS}
            if metrics:
                assessments.append({"kind": kind, "metrics": metrics})
        action = obs["action"]
        protocol = protocols.get(obs["observation_id"], {})
        scope = {"tool_id": action["tool_id"], "runtime_version": action["runtime_version"],
                 "version_evidence_status": action["version_evidence_status"],
                 "observation_scope": obs["scope"]}
        for name in ("search_region_class", "box_relation_to_declared_target_box", "attempt_class"):
            if name in protocol:
                scope[name] = protocol[name]
        lessons.append({
            "id": obs["observation_id"], "status": "supported",
            "statement": action["tool_id"] + ": " + obs["scope"],
            "scope": scope,
            "execution": {name: obs["execution"].get(name)
                          for name in ("execution_status", "parsing_status", "return_code", "error")},
            "returned_count": obs["outputs"]["returned_count"], "assessments": assessments,
            "phase_seconds": {name: metric_view(value) for name, value in obs["cost"]["phase_seconds"].items()},
            "quality_flags": obs["quality_flags"],
            "parameters_available": bool(action["parameters"]),
        })
    task = card["task_context"]
    identity = card["identity"]
    return {
        "source": "benchmark", "identity": identity,
        "cosine_similarity": similarity,
        "scope": {
            "task_type": task["task_type"], "scenario": task["scenario"],
            "receptor_condition": task["receptor_condition"],
            "site_provided": task["site_information"]["provided"],
            "ligand_provided": task["ligand_information"]["provided_for_prediction"],
            "source_input_sha256": card["embeddings"][0]["source_input_sha256"],
            "group_id": card["extensions"].get("dataset_split", {}).get("group_id"),
            "applicability": "neighbor evidence is a tentative proposal for the current protein",
        },
        "summary": card["comparability"]["supported_interpretation"],
        "protocol_differences": card["comparability"]["differences"],
        "knowledge": {"lessons": lessons, "parameter_sets": [],
                      "open_questions": card["comparability"]["unresolved_items"]},
        "evidence_record": {
            "path": str(bundle / "memory/agent_memory.jsonl"),
            "sha256": manifest["artifacts"]["agent_memory.jsonl"]["sha256"],
            "context_id": identity["context_id"], "record_revision": identity["record_revision"],
        },
        "details": {"script": str(bundle / "read_memory.py"),
                    "args": ["show", "--task", task["task_type"], "--dataset", identity["dataset_id"],
                             "--target", identity["target_id"], "--view", "execution"],
                    "note": "add --tool TOOL_ID for parameters; use evidence for original sources"},
    }


def load_planning_module():
    """Load the sibling helper for both direct CLI and file-based imports."""
    path = Path(__file__).with_name("sbdd_memory_light.py")
    spec = importlib.util.spec_from_file_location("sbdd_memory_light", path)
    module = importlib.util.module_from_spec(spec)
    old_bytecode = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = old_bytecode
    return module


def load_reader(bundle: Path):
    """Use the release reader, including its fixed-root and hash checks."""
    return load_planning_module().load_source_reader(bundle)


def benchmark_recall(args, file_digest) -> tuple[list[dict], dict]:
    bundle = args.benchmark_bundle.expanduser().resolve()
    module, reader = load_reader(bundle)
    excluded = set(args.exclude_case)
    try:
        benchmark_view = getattr(args, "benchmark_view", "planning")
        planning_memory = None
        planning_sources = set()
        if benchmark_view == "planning":
            planning = load_planning_module()

            try:
                planning_memory = planning.LightweightMemory(bundle, reader)
            except FileNotFoundError:
                # A not-yet-generated sidecar may be projected from the verified
                # original. Invalid or stale materialized content must fail.
                pass
        if args.query_npz:
            metadata = json.loads(args.query_metadata.read_text(encoding="utf-8"))
            encoder = module.encoder_id_from_metadata(metadata)
            source = metadata["source_record"]
            if source.get("task_type") and source["task_type"] != args.task:
                raise ValueError("Query metadata task_type differs from --task")
            declared_npz = metadata.get("embedding_sha256")
            if declared_npz and file_digest(args.query_npz) != declared_npz:
                raise ValueError("Query NPZ hash differs from metadata")
            declared_receptor = metadata.get("source_input_sha256") or source.get("protein_sha256")
            if args.receptor_file and declared_receptor and file_digest(args.receptor_file) != declared_receptor:
                raise ValueError("Current receptor hash differs from query embedding source")
            source_case = None
            if source.get("dataset") and source.get("target"):
                dataset = module.ALIASES.get(source["dataset"], source["dataset"])
                source_case = dataset + "/" + source["target"]
                excluded.add(source_case)
            with module.np.load(args.query_npz, allow_pickle=False) as arrays:
                vector = arrays["protein_embedding_normalized"]
            query = {"npz": str(args.query_npz.resolve()), "metadata": str(args.query_metadata.resolve()),
                     "npz_hash_checked": bool(declared_npz),
                     "receptor_hash_checked": bool(args.receptor_file and declared_receptor),
                     "case_id": source_case}
        else:
            dataset = module.ALIASES.get(args.dataset, args.dataset)
            card = reader.get(args.task, dataset, args.target)
            vector = reader.vector_for(card)
            encoder = card["embeddings"][0]["encoder"]["compatibility_id"]
            excluded.add(card["identity"]["case_id"])
            query = {"case_id": card["identity"]["case_id"]}
        compatible = reader.index.execute(
            "SELECT COUNT(*) FROM cards WHERE task_type=? AND encoder_id=?", (args.task, encoder)
        ).fetchone()[0]
        if not compatible:
            raise ValueError("No compatible benchmark encoder for this task; check model, dtype, windows and pLDDT masking")
        neighbors = reader.neighbors(args.task, vector, encoder, args.k, excluded)
        summaries = []
        for neighbor in neighbors:
            identity = neighbor["historical_memory"]["identity"]
            if benchmark_view == "summary":
                card = reader.get(args.task, identity["dataset_id"], identity["target_id"])
                summaries.append(benchmark_summary(card, neighbor["cosine_similarity"], bundle, reader.manifest))
                continue
            try:
                if planning_memory is None:
                    raise FileNotFoundError("Planning sidecar has not been generated")
                summary = copy.deepcopy(planning_memory.get(identity["context_id"]))
                planning_sources.add("materialized")
            except FileNotFoundError:
                card = reader.get(args.task, identity["dataset_id"], identity["target_id"])
                summary = planning.project_card(card)
                planning_sources.add("dynamic")
            summary["view"] = "planning"
            summary["schema_version"] = planning.FORMAT_VERSION
            summary["cosine_similarity"] = neighbor["cosine_similarity"]
            summary["details"]["script"] = str(bundle / "read_memory.py")
            summary["evidence_record"].update({
                "path": str(bundle / "memory/agent_memory.jsonl"),
                "sha256": reader.manifest["artifacts"]["agent_memory.jsonl"]["sha256"],
            })
            summaries.append(planning.planning_view(summary))
        status = {
            "status": "ok" if summaries else "no_matches", "bundle": str(bundle),
            "release": reader.manifest.get("memory_release"), "query": query,
            "task_type": args.task, "encoder_compatibility_id": encoder,
            "excluded_cases": sorted(excluded), "requested_k": args.k,
            "retrieved": len(summaries), "output_cards": len(summaries), "omitted_cards": 0,
            "view": benchmark_view,
            "selection_policy": "cosine neighbors; no success probability or best-tool label",
        }
        if benchmark_view == "planning":
            status["planning_directory"] = str(bundle / "planning")
            status["planning_source"] = (
                next(iter(planning_sources)) if len(planning_sources) == 1
                else "mixed" if planning_sources else "not_used"
            )
        return summaries, status
    finally:
        for name in ("index", "locations", "hash_cache", "runtime_index"):
            connection = getattr(reader, name, None)
            if connection is not None:
                connection.close()


def bound_output(result: dict, max_chars: int) -> dict:
    """Remove whole summary entries with explicit counts, keeping scope intact.

    Source facts are never changed. Planning cards keep every tool and are
    omitted as whole cards when necessary. Omitted content is accessible through
    detail locators and must not be treated as absent evidence.
    """
    for card in result["cards"]:
        if "knowledge" in card:
            card["omitted"] = {kind: 0 for kind in ("lessons", "parameter_sets", "open_questions")}
    while len(serialize(result)) > max_chars:
        candidates = [(len(serialize(card["knowledge"][kind][-1])), card, kind)
                      for card in result["cards"]
                      for kind in ("lessons", "parameter_sets", "open_questions")
                      if card.get("knowledge", {}).get(kind)]
        if candidates:
            _, card, kind = max(candidates, key=lambda item: item[0])
            card["knowledge"][kind].pop()
            card["omitted"][kind] += 1
            continue
        planning_cards = [card for card in result["cards"] if card.get("view") == "planning"]
        if not planning_cards:
            raise ValueError("Summary metadata exceeds --max-chars; use a smaller -k or narrower identity query")
        # Retrieval is already ordered by similarity. Drop the last whole card,
        # preserving complete tool tables for the strongest remaining neighbors.
        result["cards"].remove(planning_cards[-1])
        benchmark = result["sources"]["benchmark"]
        benchmark["omitted_cards"] += 1
        benchmark["output_cards"] -= 1
    result["summary_incomplete"] = (
        any(any(c.get("omitted", {}).values()) for c in result["cards"])
        or bool(result["sources"]["benchmark"].get("omitted_cards", 0))
    )
    # The final flag replaces the preallocated false value with a shorter true.
    return result

"""Small, traceable planning cards derived from the frozen SBDD memory release.

The source release reader owns validation, runtime recovery and retrieval.
This module projects historical facts without inferring a best tool or changing
missing metrics, and reads the separately hashed planning sidecar.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


FORMAT_VERSION = "sbdd_planning_card/1"
_SCRIPT_PATH = Path(__file__).resolve()
if _SCRIPT_PATH.name == "reader.py":
    # build_sbdd_memory_light copies this file to bundle/planning/reader.py.
    # That standalone copy must not depend on the source checkout's shared/.
    # A staged copy outside the bundle can still use the explicit --bundle flag.
    DEFAULT_BUNDLE = _SCRIPT_PATH.parent.parent
else:
    sys.path.insert(0, str(_SCRIPT_PATH.parents[2] / "shared"))
    from deployment_config import deployment_path

    DEFAULT_BUNDLE = deployment_path(
        "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-docking-memory"
    )
METRIC_COLUMNS = {
    "docking": [
        "rmsd", "posebusters_valid", "joint_success", "top3_rmsd_le_2a",
        "top5_rmsd_le_2a", "best_observed_rmsd", "pb_failed_check_count",
        "pb_unavailable_check_count",
    ],
    "pocket_prediction": [
        "top1_dca_angstrom", "top1_dcc_angstrom", "top1_dca_le_4",
        "top3_dca_le_4", "top5_dca_le_4", "any_best_dca_angstrom",
        "first_hit_rank_dca_le_4", "gold_mapping_coverage",
        "point5_top1_rro", "point5_top1_jaccard",
    ],
}
USAGE = (
    "Historical descriptive evidence; transfer to the current protein is tentative. "
    "Check current structure, site and protocol before adopting neighbor evidence. "
    "Cosine similarity is not success probability. Tools have different protocols. "
    "RMSD/DCA/DCC are angstrom, time is seconds; other metrics are booleans/counts/fractions. "
    "metrics follow metric_columns; cells are [value,status,optional reason]. "
    "run is [execution_status,parsing_status,return_code]. "
    "Top-K/best RMSD has no Top-K PoseBusters guarantee. "
    "Times are partial observed costs on uncontrolled hardware, not end-to-end budgets. "
    "Fetch execution/evidence for exact parameters, retries and omitted metrics."
)
PLANNING_USAGE = (
    "Historical transfer is tentative; cosine is not success probability. "
    "Check current structure, site and protocol before adopting neighbor evidence. "
    "RMSD/DCA/DCC: angstrom. Time: seconds, partial costs on uncontrolled hardware, not end-to-end budgets. "
    "Top-K RMSD has no Top-K PoseBusters guarantee. quality_cohorts are source comparison eligibility. "
    "tools rows follow tool_columns. Null run/flags/signals use common; protocol/cost overlay common; "
    "scope is extra.scope or common.scope. values follow metric_columns; states[index] gives "
    "[status,optional reason], default states['*'] or [common.metric_status]. "
    "seconds keys index time_columns; cells=[value,status,optional reason]. "
    "Fetch details for parameters, full timer exclusions and source refs."
)


def serialize(value) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n"


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def metric_cell(metric: dict | None, absent_status="not_recorded") -> list:
    if metric is None:
        return [None, absent_status]
    cell = [copy.deepcopy(metric.get("value")), metric.get("status")]
    if metric.get("reason"):
        cell.append(metric["reason"])
    return cell


def cost_view(cost: dict, task: str) -> tuple[dict, dict]:
    phases = cost.get("phase_seconds", {})
    preferred = (
        ["engine_elapsed_seconds", "inference_loop_seconds", "total_seconds"]
        if task == "docking" else
        ["prepare_plus_predict_seconds", "prediction_seconds", "predict_seconds"]
    )
    chosen = next((name for name in preferred if name in phases
                   and phases[name].get("value") is not None), None)
    if chosen is None:
        chosen = next((name for name in preferred if name in phases), None)
    seconds = {chosen: metric_cell(phases[chosen])} if chosen else {}
    # The combined preparation/prediction timer already contains both aliases.
    # Parsing is separately observed and is never silently added to another timer.
    if task == "pocket_prediction" and "parse_seconds" in phases:
        seconds["parse_seconds"] = metric_cell(phases["parse_seconds"])
    scope = {"cumulative": cost.get("cumulative_cost_status"),
             "timeout_censored": cost.get("timeout_censored")}
    if cost.get("timeout_censored"):
        scope["timeout_limit_seconds"] = cost.get("timeout_limit_seconds")
    recovery = cost.get("runtime_recovery", {})
    measurement = next((m for m in recovery.get("measurements", [])
                        if m["metric_name"] == chosen), None)
    if measurement:
        scope["definition_id"] = phases[chosen].get("definition_id")
        scope["timer_scope"] = measurement["timer_scope"]
        scope["clock_kind"] = measurement["clock_kind"]
        scope["source_ref"] = copy.deepcopy(measurement["source_ref"])
    elif chosen == "prepare_plus_predict_seconds":
        scope["timer_scope"] = "prepare+predict; excludes parse/queue/uncaptured history"
    elif chosen:
        scope["timer_scope"] = chosen + "; no end-to-end wall-time claim"
    return seconds, scope


def project_card(card: dict) -> dict:
    task = card["task_context"]
    task_type = task["task_type"]
    columns = METRIC_COLUMNS[task_type]
    selected = {s["observation_id"] for s in card["published_selections"]}
    protocols = {p["observation_id"]: p["protocol_labels"]
                 for p in card["comparability"].get("per_tool_protocols", [])}
    assessments = {}
    for assessment in card["assessments"]:
        assessments.setdefault(assessment["observation_id"], []).append(assessment)
    claims = {}
    for claim in card.get("derived_experience", []):
        for oid in claim.get("supporting_observation_ids", []):
            claims.setdefault(oid, []).append(claim["claim"])
    counts = {}
    for obs in card["tool_observations"]:
        tool = obs["action"]["tool_id"]
        counts[tool] = counts.get(tool, 0) + 1
    tools = []
    for obs in card["tool_observations"]:
        oid = obs["observation_id"]
        if oid not in selected and obs["scope"] != "original_full_run":
            continue
        action, execution = obs["action"], obs["execution"]
        candidates = obs["outputs"]["candidates"]
        metrics = {}
        for assessment in assessments.get(oid, []):
            kind = assessment["kind"]
            top1 = (kind == "pose_rmsd" and candidates
                    and assessment["candidate_ids"] == [candidates[0]["candidate_id"]])
            if top1 or kind in {"top1_posebusters", "pose_set_summary", "pocket_target_summary"}:
                for name, metric in assessment["metrics"].items():
                    if name in columns:
                        if name in metrics and metrics[name] != metric:
                            raise ValueError("Conflicting selected metrics for " + oid + ": " + name)
                        metrics[name] = metric
        protocol = protocols.get(oid, {})
        seconds, cost_scope = cost_view(obs["cost"], task_type)
        flags = copy.deepcopy(obs.get("quality_flags", []))
        spatial_inapplicable = "per_residue_scores_not_spatial_pocket_candidates" in flags
        tool = {
            "tool": action["tool_id"], "observation_id": oid, "scope": obs["scope"],
            "run": [execution.get("execution_status"), execution.get("parsing_status"),
                    execution.get("return_code")],
            "returned_count": obs["outputs"]["returned_count"],
            "metrics": [metric_cell(metrics.get(name),
                                    "not_applicable" if spatial_inapplicable else "not_recorded")
                        for name in columns],
            "seconds": seconds, "cost_scope": cost_scope,
            "protocol": {
                "region": protocol.get("search_region_class"),
                "box": protocol.get("box_relation_to_declared_target_box"),
                "attempt": protocol.get("attempt_class"),
                "version": action.get("runtime_version"),
                "version_status": action.get("version_evidence_status"),
                "reference_used": obs.get("effective_inputs", {}).get("reference_used_for_prediction"),
            },
            "flags": flags, "signals": list(dict.fromkeys(claims.get(oid, []))),
        }
        if spatial_inapplicable:
            tool["output_kind"] = "residue_scores"
        if counts[action["tool_id"]] > 1:
            tool["additional_observation_count"] = counts[action["tool_id"]] - 1
        if execution.get("error"):
            error = str(execution["error"])
            tool["error"] = error[:240]
            if len(error) > 240:
                tool["error_truncated"] = True
        tools.append(tool)
    features = {}
    feature_names = ["residue_count", "chain_count", "pocket_residue_count",
                     "modified_amino_acid_count", "incomplete_backbone_residue_count",
                     "unknown_amino_acid_count", "confidence_masked_residue_count",
                     "confidence_mask_enabled", "plddt_summary"]
    for name in feature_names:
        feature = card.get("features", {}).get(name)
        if feature and feature.get("availability", {}).get("usage_role") == "query_feature" \
                and feature["availability"].get("available_after") == "input_ready":
            features[name] = metric_cell(feature)
    ligand = card.get("features", {}).get("intended_ligand_descriptors")
    if ligand and ligand.get("availability", {}).get("usage_role") == "query_feature" \
            and ligand["availability"].get("available_after") == "input_ready":
        features["ligand"] = {
            "status": ligand.get("status"),
            **{key: ligand["value"].get(key) for key in
               ["heavy_atom_count", "rotatable_bonds", "formal_charge", "molecular_weight"]}
        } if isinstance(ligand.get("value"), dict) else {"status": ligand.get("status"), "value": None}
    identity = card["identity"]
    comparison = card["comparability"]
    return {
        "schema_version": FORMAT_VERSION, "source": "benchmark", "view": "planning",
        "identity": {name: identity[name] for name in
                     ["case_id", "context_id", "dataset_id", "target_id", "record_revision"]},
        "task": {"task_type": task_type, "scenario": task["scenario"],
                 "receptor_condition": task["receptor_condition"],
                 "site_provided": task["site_information"]["provided"],
                 "ligand_provided": task["ligand_information"]["provided_for_prediction"]},
        "group_id": card.get("extensions", {}).get("dataset_split", {}).get("group_id"),
        "quality_cohorts": {
            membership["cohort_id"]: [membership["included"], membership.get("reason")]
            for membership in card.get("cohort_membership", [])
            if membership["cohort_id"] in {
                "docking_697_final", "top1_pb_5576", "pocket_five_shared_success_575"
            }
        },
        "features": features, "metric_columns": list(columns), "tools": tools,
        "comparability": {"interpretation": comparison["supported_interpretation"],
                          "differences": copy.deepcopy(comparison["differences"]),
                          "unresolved": copy.deepcopy(comparison["unresolved_items"])},
        "evidence_record": {"context_id": identity["context_id"],
                            "record_revision": identity["record_revision"]},
        "details": {"script": "read_memory.py", "args": ["show", "--task", task_type,
                    "--dataset", identity["dataset_id"], "--target", identity["target_id"],
                    "--view", "execution"],
                    "note": "Add --tool TOOL_ID for parameters; evidence resolves artifact IDs."},
    }


def planning_view(card: dict) -> dict:
    """Show every tool with compact defaults, retaining exact metric values/states.

    Materialized cards keep their original dictionary schema. This presentation
    omits only the source-reference detail and timer exclusion prose already
    available through the card's execution/evidence locator.
    """
    result = copy.deepcopy(card)
    tools = result.get("tools", [])
    if not tools or isinstance(tools[0], list):
        return result

    def most_common(values):
        groups = {}
        for value in values:
            key = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            if key not in groups:
                groups[key] = [0, value]
            groups[key][0] += 1
        return max(groups.values(), key=lambda pair: pair[0])

    def defaults(items):
        keys = set(items[0]).intersection(*(set(item) for item in items[1:]))
        shared = {}
        for key in sorted(keys):
            count, value = most_common([item[key] for item in items])
            if count > 1:
                shared[key] = value
        return shared

    def differences(item, shared):
        return {key: value for key, value in item.items()
                if key not in shared or value != shared[key]}

    common = {"metric_status": most_common([
        cell[1] for tool in tools for cell in tool["metrics"]
    ])[1]}
    for name in ("scope", "run", "flags", "signals"):
        count, value = most_common([tool[name] for tool in tools])
        if count > 1:
            common[name] = value
    common["protocol"] = defaults([tool["protocol"] for tool in tools])
    costs = []
    for tool in tools:
        cost = {key: value for key, value in tool.get("cost_scope", {}).items()
                if key not in {"source_ref", "definition_id"}}
        if isinstance(cost.get("timer_scope"), str):
            cost["timer_scope"] = cost["timer_scope"].split(";", 1)[0]
        costs.append(cost)
    common["cost"] = defaults(costs)
    columns = ["tool", "observation_id", "run", "returned_count", "values", "states",
               "seconds", "protocol", "cost", "flags", "signals", "extra"]
    consumed = {"tool", "observation_id", "run", "returned_count", "metrics", "seconds",
                "protocol", "cost_scope", "flags", "signals", "scope"}
    rows = []
    phases = list(dict.fromkeys(phase for tool in tools for phase in tool["seconds"]))
    for tool, cost in zip(tools, costs):
        states = {str(index): cell[1:] for index, cell in enumerate(tool["metrics"])
                  if cell[1] != common["metric_status"] or len(cell) > 2}
        row_statuses = [cell[1:] for cell in tool["metrics"]]
        if row_statuses and all(state == row_statuses[0] for state in row_statuses) and states:
            states = {"*": row_statuses[0]}
        extra = {key: value for key, value in tool.items() if key not in consumed}
        if "scope" not in common or tool["scope"] != common["scope"]:
            extra["scope"] = tool["scope"]
        row = {
            "tool": tool["tool"], "observation_id": tool["observation_id"],
            "returned_count": tool["returned_count"],
            "values": [cell[0] for cell in tool["metrics"]], "states": states,
            "seconds": {str(phases.index(phase)): cell for phase, cell in tool["seconds"].items()},
            "protocol": differences(tool["protocol"], common["protocol"]),
            "cost": differences(cost, common["cost"]), "extra": extra,
        }
        for name in ("run", "flags", "signals"):
            row[name] = None if name in common and tool[name] == common[name] else tool[name]
        rows.append([row[name] for name in columns])
    result["schema_version"] = "sbdd_planning_view/1"
    result["tool_columns"] = columns
    result["time_columns"] = phases
    result["common"] = common
    result["tools"] = rows
    result["details_omitted"] = ["cost.source_ref", "cost.definition_id", "cost.timer_scope_exclusion_prose"]
    # Identity already binds the evidence locator; avoid repeating its digest.
    for name in ("context_id", "record_revision"):
        result.get("evidence_record", {}).pop(name, None)
    return result


def load_source_reader(bundle: Path):
    """Import the release's reader without leaving module-name collisions."""
    bundle = bundle.expanduser().resolve()
    path = bundle / "code/bundle_reader.py"
    spec = importlib.util.spec_from_file_location("sbdd_benchmark_reader", path)
    module = importlib.util.module_from_spec(spec)
    old_path, old_bytecode = list(sys.path), sys.dont_write_bytecode
    local_modules = {p.stem for p in path.parent.glob("*.py")}
    saved = {name: sys.modules.pop(name) for name in local_modules if name in sys.modules}
    try:
        sys.path.insert(0, str(path.parent))
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
        reader = module.DecisionMemory(bundle)
    finally:
        sys.path[:] = old_path
        sys.dont_write_bytecode = old_bytecode
        for name in local_modules:
            sys.modules.pop(name, None)
        sys.modules.update(saved)
    return module, reader


def close_reader(reader):
    for name in ("index", "locations", "hash_cache", "runtime_index"):
        connection = getattr(reader, name, None)
        if connection is not None:
            connection.close()


def source_binding(reader) -> dict:
    return {
        "memory_release": reader.manifest.get("memory_release"),
        "memory_artifacts": {name: reader.manifest["artifacts"][name]["sha256"]
                             for name in ("agent_memory.jsonl", "memory_index.sqlite")},
        "essential_files": copy.deepcopy(reader.bundle_manifest["essential_files"]),
    }


class LightweightMemory:
    def __init__(self, bundle: Path, reader, directory: Path | None = None):
        self.bundle = bundle.expanduser().resolve()
        self.directory = directory or self.bundle / "planning"
        manifest_path = self.directory / "manifest.json"
        if not manifest_path.is_file():
            if self.directory.exists():
                raise RuntimeError("Planning directory exists but its manifest is missing")
            raise FileNotFoundError("Planning cards are not materialized")
        try:
            self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if self.manifest["schema_version"] != FORMAT_VERSION or self.manifest["status"] != "passed":
                raise RuntimeError("Unsupported or unvalidated planning release")
            if self.manifest["source"] != source_binding(reader):
                raise RuntimeError("Planning cards do not bind to this source release")
            for name in ("cards.jsonl", "index.json"):
                path = self.directory / name
                if not path.is_file() or digest(path) != self.manifest["artifacts"][name]["sha256"]:
                    raise RuntimeError("Planning content hash mismatch: " + name)
            self.entries = json.loads((self.directory / "index.json").read_text(encoding="utf-8"))
            expected = {row[0] for row in reader.index.execute("SELECT context_id FROM cards")}
            if set(self.entries) != expected or len(self.entries) != self.manifest["counts"]["contexts"]:
                raise RuntimeError("Planning cards do not cover the exact source context set")
            size = (self.directory / "cards.jsonl").stat().st_size
            for offset, length in self.entries.values():
                if (type(offset) is not int or type(length) is not int
                        or offset < 0 or length < 1 or offset + length > size):
                    raise RuntimeError("Planning index has an invalid byte range")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise RuntimeError("Invalid planning release: " + str(exc)) from exc

    def get(self, context_id: str) -> dict:
        offset, length = self.entries[context_id]
        try:
            with (self.directory / "cards.jsonl").open("rb") as handle:
                handle.seek(offset)
                card = json.loads(handle.read(length))
            if card["identity"]["context_id"] != context_id or card["schema_version"] != FORMAT_VERSION:
                raise RuntimeError("Planning index points to a different context")
            card["details"]["script"] = str(self.bundle / "read_memory.py")
            card["evidence_record"].update({
                "path": str(self.bundle / "memory/agent_memory.jsonl"),
                "sha256": self.manifest["source"]["memory_artifacts"]["agent_memory.jsonl"],
            })
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise RuntimeError("Invalid planning card: " + str(exc)) from exc
        return card


def query_vector(args, module, reader):
    excluded = set(args.exclude_case)
    if args.query_npz:
        metadata = json.loads(args.query_metadata.read_text(encoding="utf-8"))
        source = metadata["source_record"]
        if source.get("task_type") and source["task_type"] != args.task:
            raise ValueError("Query metadata task_type differs from --task")
        declared = metadata.get("embedding_sha256")
        if declared and digest(args.query_npz) != declared:
            raise ValueError("Query NPZ hash differs from metadata")
        receptor_hash = metadata.get("source_input_sha256") or source.get("protein_sha256")
        if args.receptor_file and receptor_hash and digest(args.receptor_file) != receptor_hash:
            raise ValueError("Current receptor hash differs from query embedding source")
        if source.get("dataset") and source.get("target"):
            excluded.add(module.ALIASES.get(source["dataset"], source["dataset"]) + "/" + source["target"])
        with module.np.load(args.query_npz, allow_pickle=False) as arrays:
            vector = arrays["protein_embedding_normalized"]
        encoder = module.encoder_id_from_metadata(metadata)
        query = {"npz": str(args.query_npz.resolve()), "case_excluded": sorted(excluded),
                 "npz_hash_checked": bool(declared),
                 "receptor_hash_checked": bool(args.receptor_file and receptor_hash)}
    else:
        dataset = module.ALIASES.get(args.dataset, args.dataset)
        card = reader.get(args.task, dataset, args.target)
        vector = reader.vector_for(card)
        encoder = card["embeddings"][0]["encoder"]["compatibility_id"]
        excluded.add(card["identity"]["case_id"])
        query = {"case_id": card["identity"]["case_id"]}
    if not reader.index.execute("SELECT 1 FROM cards WHERE task_type=? AND encoder_id=?",
                                (args.task, encoder)).fetchone():
        raise ValueError("No compatible benchmark encoder for this task")
    return vector, encoder, excluded, query


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--planning-dir", type=Path)
    sub = parser.add_subparsers(dest="command", required=True)
    show = sub.add_parser("show")
    show.add_argument("--task", required=True, choices=METRIC_COLUMNS)
    show.add_argument("--dataset", required=True)
    show.add_argument("--target", required=True)
    show.add_argument("--tool")
    show.add_argument("--view", choices=("planning", "compare", "execution"), default="planning")
    near = sub.add_parser("neighbors")
    near.add_argument("--task", required=True, choices=METRIC_COLUMNS)
    near.add_argument("--dataset")
    near.add_argument("--target")
    near.add_argument("--query-npz", type=Path)
    near.add_argument("--query-metadata", type=Path)
    near.add_argument("--receptor-file", type=Path)
    near.add_argument("--exclude-case", action="append", default=[])
    near.add_argument("-k", type=int, default=3)
    args = parser.parse_args(argv)
    reader = None
    try:
        module, reader = load_source_reader(args.bundle)
        store = LightweightMemory(args.bundle, reader, args.planning_dir)
        if args.command == "show":
            dataset = module.ALIASES.get(args.dataset, args.dataset)
            entry = reader.index.execute("SELECT context_id FROM cards WHERE task_type=? AND dataset=? AND target=?",
                                         (args.task, dataset, args.target)).fetchone()
            if entry is None:
                raise ValueError("Unknown task/dataset/target")
            if args.view == "planning":
                card = store.get(entry[0])
                if args.tool:
                    all_tools = card["tools"]
                    card["tools"] = [row for row in all_tools if row["tool"] == args.tool]
                    if not card["tools"]:
                        raise ValueError("Tool is not present in this card")
                    card["omitted_tools"] = len(all_tools) - len(card["tools"])
                result = {"usage": PLANNING_USAGE, "card": planning_view(card)}
            else:
                result = reader.view(reader.get(args.task, dataset, args.target), args.view, args.tool)
        else:
            if not 1 <= args.k <= 10:
                raise ValueError("-k must be between 1 and 10")
            if bool(args.query_npz) != bool(args.query_metadata):
                raise ValueError("Supply both --query-npz and --query-metadata")
            if bool(args.dataset) != bool(args.target) or bool(args.query_npz) == bool(args.dataset):
                raise ValueError("Supply dataset/target or NPZ/metadata, exclusively")
            vector, encoder, excluded, query = query_vector(args, module, reader)
            neighbors = reader.neighbors(args.task, vector, encoder, args.k, excluded)
            cards = []
            for neighbor in neighbors:
                card = store.get(neighbor["historical_memory"]["identity"]["context_id"])
                card["cosine_similarity"] = neighbor["cosine_similarity"]
                cards.append(planning_view(card))
            result = {"usage": PLANNING_USAGE, "query": query, "encoder_compatibility_id": encoder,
                      "excluded_cases": sorted(excluded), "cards": cards}
        print(serialize(result), end="")
        return 0
    except (OSError, ValueError, RuntimeError, KeyError, ImportError) as exc:
        print(serialize({"success": False, "error": str(exc)}), end="")
        return 1
    finally:
        if reader is not None:
            close_reader(reader)


if __name__ == "__main__":
    sys.exit(main())

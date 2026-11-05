"""Protein-memory pocket fusion, with an explicit memory-unavailable fallback.

Protein similarity selects detectors and weights pocket reliability. No detector,
evaluator, parameter search, or pipeline runner is imported or run here.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

TOOLS = ("deeppocket", "fpocket", "p2rank", "pykvfinder")
PAPER_TOOLS = (*TOOLS, "autosite")
STRATEGY = "deeppocket_top1_local_rank_fusion_v1"
MEMORY_STRATEGY = "protein_memory_local_fusion_v1"
INDEX_TYPE = "pocket_protein_reliability"
BEST_DETECTOR_POLICY = "minimum_native_top1_dca_then_lexicographic_tool_id_v1"
INFERENCE_NOTES = [
    "Best detector is inferred from minimum raw native Top-1 DCA; exact ties use lexicographic tool ID.",
    "Only cases complete for all five memory detectors enter retrieval and priors; exclusions apply to both.",
    "Repeated recommendations are retained in the five-entry multiset; each unique detector executes once.",
    "DeepPocket always executes to supply the default anchor, even if no neighbor recommends it.",
    "Every matching native Top-3 pocket contributes to memory fusion, including lower-ranked anchor-tool pockets.",
]
POINT5 = "receptor_protein_heavy_atoms_within_5A_of_native_support_points"
PARAMETERS = {"top_k": 3, "center_threshold_A": 5.0,
              "residue_IoU_threshold": 0.3, "residue_policy": "anchor_point5",
              "weight_policy": "reciprocal_native_rank", "family_weight_policy": "max_member_weight"}


def _sibling(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def default_index_path():
    spec = importlib.util.spec_from_file_location(
        "deployment_config", Path(__file__).resolve().parents[2] / "shared" / "deployment_config.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.deployment_path("${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json")


def _active_tools(payload):
    detectors = payload.get("detectors")
    if not isinstance(detectors, dict) or set(detectors) - set(PAPER_TOOLS) or "deeppocket" not in detectors:
        raise ValueError("detectors must contain DeepPocket and only supported pocket tools")
    return tuple(tool for tool in PAPER_TOOLS if tool in detectors)


def _sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write_json(path, payload):
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _xyz(value):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError("expected a three-coordinate vector")
    result = [float(x) for x in value]
    if not all(math.isfinite(x) for x in result):
        raise ValueError("nonfinite pocket coordinates")
    return result


def _identity(record):
    if record.get("identity_ambiguous"):
        raise ValueError("ambiguous residue identity cannot establish point5 overlap")
    number = record["resnum"]
    if isinstance(number, bool) or int(number) != float(number):
        raise ValueError("invalid residue number")
    return (str(record.get("chain", "")).strip(), int(number),
            str(record.get("icode", "")).strip(), str(record["resname"]).strip().upper())


def _protein_atoms(path):
    """First model, amino acids, heavy atoms, highest-occupancy alternate atoms."""
    from Bio.Data.PDBData import protein_letters_3to1_extended

    residues, serials = {}, {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        xyz = _xyz([line[30:38], line[38:46], line[46:54]])
        serials[int(line[6:11])] = xyz
        name, atom = line[17:20].strip().upper(), line[12:16].strip()
        element = line[76:78].strip().upper() or atom.lstrip("0123456789")[:1]
        if name not in protein_letters_3to1_extended or element in {"H", "D", "T"}:
            continue
        key = (line[21].strip(), int(line[22:26]), line[26].strip(), name)
        priority = (float(line[54:60].strip() or 0), line[16] == " ", line[16] == "A")
        atoms = residues.setdefault(key, {})
        if atom not in atoms or priority > atoms[atom][0]:
            atoms[atom] = (priority, xyz)
    residues = {key: atoms for key, atoms in residues.items()
                if "CA" in atoms or {"N", "C"} <= atoms.keys()}
    if not residues:
        raise ValueError("receptor has no usable protein heavy atoms")
    return {key: [value[1] for value in atoms.values()] for key, atoms in residues.items()}, serials


def _point_file(path, resname=None):
    points = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")) and (resname is None or line[17:20].strip() == resname):
            points.append(_xyz([line[30:38], line[38:46], line[46:54]]))
    return points


def _report_candidates(payload, base):
    """Adapt the existing Step-1 report without using its conventional consensus."""
    report_path = (base / payload["step1_report"]).resolve()
    report = _read_json(report_path)
    sources = payload["sources"]
    if report.get("enabled_tool_failures"):
        raise ValueError("detector preparation report contains enabled tool failures")
    detectors = {}
    detection = report.get("pocket_detection")
    if not isinstance(detection, dict):
        raise ValueError("Step-1 report must contain pocket_detection results")
    tools = [tool for tool in PAPER_TOOLS if detection.get(tool) is not None]
    if any(not isinstance(detection[tool], dict) for tool in tools):
        raise ValueError("Step-1 detector results must be objects or null when disabled")
    plan = detection.get("memory_plan")
    if plan is not None:
        planned = plan.get("execution_tools") if isinstance(plan, dict) else None
        if (not isinstance(planned, list) or any(not isinstance(tool, str) for tool in planned)
                or len(set(planned)) != len(planned) or set(planned) - set(PAPER_TOOLS)
                or "deeppocket" not in planned or set(planned) != set(tools)):
            raise ValueError("Step-1 detector results must exactly match memory_plan.execution_tools")
    if "deeppocket" not in tools:
        raise ValueError("Step-1 report must include the DeepPocket default anchor")
    p2 = {}
    if "p2rank" in tools:
        with (base / sources["p2rank_predictions"]).open(newline="") as handle:
            p2 = [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(handle)]
        p2 = {int(row["rank"]): row for row in p2}
    for tool in tools:
        result = detection[tool]
        if result.get("success") is False:
            raise ValueError(f"{tool} detector failed")
        if tool in ("deeppocket", "autosite") and (result.get("input") or {}).get("protein_sha256") != payload["receptor_sha256"]:
            raise ValueError(f"{tool} report input SHA-256 does not match candidate receptor")
        candidates = []
        for source in sorted(result.get("all_pockets", []), key=lambda item: item["rank"]):
            item = dict(source)
            rank = int(item["rank"])
            if rank > 3:
                continue
            if tool in ("fpocket", "deeppocket"):
                seed = int(item["original_fpocket_rank"]) if tool == "deeppocket" else rank
                item["support_file"] = str(base / sources[f"{tool}_pockets"] / f"pocket{seed}_vert.pqr")
            elif tool == "p2rank":
                row = p2[rank]
                serials = [int(value) for value in row["surf_atom_ids"].split()]
                if len(serials) != int(row["surf_atoms"]):
                    raise ValueError("P2Rank surface atom count mismatch")
                item["support_atom_serials"] = serials
            elif tool == "pykvfinder":
                item["support_file"] = str(base / sources["pykvfinder_points"])
                item["support_resname"] = item["cavity_id"]
            elif tool == "autosite" and not any(key in item for key in ("support_points", "support_file")):
                if not item.get("points_file"):
                    raise ValueError("AutoSite candidates require exported cluster support points")
                item["support_file"] = item["points_file"]
            candidates.append(item)
        detectors[tool] = candidates
    return {**payload, "detectors": detectors, "source_report_sha256": _sha(report_path)}


def normalize_payload(payload, receptor_file, base=None):
    """Map native support points to point5 residues in the current receptor.

    Accept detector arrays or a Step-1 report/source-path manifest. Relative
    support paths resolve against ``base``; native residue lists are not point5.
    """
    from scipy.spatial import cKDTree

    receptor, base = Path(receptor_file), Path(base or ".").resolve()
    digest = _sha(receptor)
    if payload.get("receptor_sha256") != digest:
        raise ValueError("candidate manifest receptor SHA-256 does not match current receptor")
    if "step1_report" in payload:
        payload = _report_candidates(payload, base)
    atoms, serials = _protein_atoms(receptor)
    detectors = {}
    for tool in _active_tools(payload):
        candidates = payload.get("detectors", {}).get(tool)
        if not isinstance(candidates, list) or not candidates:
            raise ValueError(f"{tool}: current detector candidates are required")
        result, ranks = [], set()
        for source in candidates:
            rank = source.get("rank")
            if isinstance(rank, bool) or not isinstance(rank, int) or rank < 1 or rank in ranks:
                raise ValueError(f"{tool}: unique positive integer native ranks required")
            ranks.add(rank)
            if rank > 3:
                continue
            item = {"rank": rank, "center": _xyz(source.get("center"))}
            if tool in ("deeppocket", "p2rank"):
                probability = float(source["probability"])
                if not math.isfinite(probability) or not 0 <= probability <= 1:
                    raise ValueError(f"{tool}: invalid native probability")
                item["probability"] = probability
            provenance = {}
            if "support_points" in source:
                points = [_xyz(point) for point in source["support_points"]]
            elif "support_atom_serials" in source:
                points = [serials[int(serial)] for serial in source["support_atom_serials"]]
            elif "support_file" in source:
                support_file = (base / source["support_file"]).resolve()
                support_digest = _sha(support_file)
                if "support_sha256" in source and source["support_sha256"] != support_digest:
                    raise ValueError(f"{tool}#{rank}: support file SHA-256 does not match detector evidence")
                points = _point_file(support_file, source.get("support_resname"))
                provenance = {"support_file": str(support_file), "support_sha256": support_digest}
            else:
                raise ValueError(f"{tool}#{rank}: native support points required; native residue lists are not point5")
            if not points:
                raise ValueError(f"{tool}#{rank}: empty support point cloud")
            tree = cKDTree(points)
            records = [{"chain": key[0], "resnum": key[1], "icode": key[2], "resname": key[3]}
                       for key, coords in atoms.items() if min(tree.query(coords)[0]) <= 5.0]
            if not records:
                raise ValueError(f"{tool}#{rank}: support cloud has no point5 receptor residues")
            item.update(residue_records=records, support_point_count=len(points), **provenance)
            result.append(item)
        result.sort(key=lambda item: item["rank"])
        if not result or result[0]["rank"] != 1:
            raise ValueError(f"{tool}: native rank 1 is required; do not relabel lower ranks")
        detectors[tool] = result
    normalized = {"schema_version": 1, "target": payload.get("target"),
                  "receptor_sha256": digest, "point5_definition": POINT5,
                  "detectors": detectors}
    if payload.get("source_report_sha256"):
        normalized["source_report_sha256"] = payload["source_report_sha256"]
    return normalized


def normalize_candidates(candidates_path, receptor_file, output=None):
    """Normalize a manifest and retain its source binding without modifying it."""
    path, receptor = Path(candidates_path), Path(receptor_file)
    normalized = normalize_payload(_read_json(path), receptor, path.parent)
    normalized.update(source_candidates=str(path.resolve()), source_candidates_sha256=_sha(path))
    if output is not None:
        if Path(output).resolve() in {path.resolve(), receptor.resolve()}:
            raise ValueError("normalization output must not overwrite its inputs")
        _write_json(output, normalized)
    return normalized


def _validate_candidates(payload, receptor_file):
    if payload.get("receptor_sha256") != _sha(receptor_file):
        raise ValueError("normalized candidates are not bound to the current receptor")
    if payload.get("point5_definition") != POINT5:
        raise ValueError("normalized point5 candidates required; run normalize-pocket-candidates")
    atoms, _ = _protein_atoms(receptor_file)
    for tool in _active_tools(payload):
        candidates = payload.get("detectors", {}).get(tool)
        if not isinstance(candidates, list) or not 1 <= len(candidates) <= 3:
            raise ValueError(f"{tool}: one to three valid native candidates required")
        ranks = [p.get("rank") for p in candidates]
        if any(isinstance(r, bool) or not isinstance(r, int) or r not in (1, 2, 3) for r in ranks):
            raise ValueError(f"{tool}: invalid native rank")
        if ranks != sorted(set(ranks)) or ranks[0] != 1:
            raise ValueError(f"{tool}: unique ascending native ranks starting at 1 required")
        for candidate in candidates:
            _xyz(candidate.get("center"))
            if tool in ("deeppocket", "p2rank"):
                probability = float(candidate["probability"])
                if not math.isfinite(probability) or not 0 <= probability <= 1:
                    raise ValueError(f"{tool}: invalid native confidence")
            records = candidate.get("residue_records") or []
            if not records or any(_identity(record) not in atoms for record in records):
                raise ValueError(f"{tool}: point5 residues are absent from current receptor")
            if int(candidate.get("support_point_count", 0)) < 1:
                raise ValueError(f"{tool}: point5 support count missing")
    return atoms


def _agreement(first, second):
    a = {_identity(r) for r in first["residue_records"]}
    b = {_identity(r) for r in second["residue_records"]}
    return math.dist(first["center"], second["center"]), len(a & b) / len(a | b) if a | b else 0.0


def _family(tool):
    return "fpocket_family" if tool in ("deeppocket", "fpocket") else tool


def _protein_vector(value):
    if not isinstance(value, (list, tuple)) or len(value) != 1280 or any(
            isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
            for x in value):
        raise ValueError("protein embedding must contain 1280 finite numbers")
    norm = math.sqrt(sum(x * x for x in value))
    if not math.isfinite(norm) or abs(norm - 1.0) > 1e-5:
        raise ValueError("protein embedding must be L2-normalized")
    return value


def _name(value, label, *, optional=False):
    if optional and value is None:
        return None
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"invalid {label}")
    return value


def _case_id(value):
    value = _name(value, "case identity", optional=True)
    if value and "/" in value:
        dataset, target = value.split("/", 1)
        aliases = {"astex": "astex_diverse", "posebusters": "posebusters_benchmark"}
        value = f"{aliases.get(dataset, dataset)}/{target}"
    return value


def infer_best_detector(metrics):
    """Choose minimum raw DCA, with stable tool-ID ties.

    Incomplete five-tool profiles do not receive an inferred winner.
    """
    if any(metrics.get(tool) is None for tool in PAPER_TOOLS):
        return None
    return min(PAPER_TOOLS, key=lambda tool: (metrics[tool], tool))


def validate_index(index):
    """Validate a compact derived index; never open its source release online."""
    if not isinstance(index, dict) or type(index.get("schema_version")) is not int or index["schema_version"] != 1:
        raise ValueError("incompatible pocket index schema")
    if index.get("index_type") != INDEX_TYPE or index.get("policy_version") != MEMORY_STRATEGY:
        raise ValueError("incompatible pocket index type or policy")
    if index.get("best_detector_policy", BEST_DETECTOR_POLICY) != BEST_DETECTOR_POLICY:
        raise ValueError("incompatible pocket best-detector derivation")
    _name(index.get("encoder_compatibility_id"), "encoder compatibility identity")
    if not isinstance(index.get("source"), dict) or index["source"].get("split", "memory") != "memory":
        raise ValueError("pocket index must identify its memory source partition")
    tools = index.get("tools")
    if not isinstance(tools, list) or not tools or any(not isinstance(t, str) for t in tools) or (
            len(tools) != len(set(tools)) or not set(TOOLS) <= set(tools) <= set(PAPER_TOOLS)):
        raise ValueError("invalid pocket index tools")
    cases = index.get("cases")
    if not isinstance(cases, list):
        raise ValueError("pocket index cases must be an array")
    seen = set()
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("invalid pocket memory case")
        identity = _case_id(case.get("case_id"))
        if identity is None or identity in seen:
            raise ValueError("missing or duplicate pocket memory case identity")
        seen.add(identity)
        _name(case.get("group_id"), "protein group")
        _protein_vector(case.get("protein"))
        metrics = case.get("top1_dca_angstrom")
        if not isinstance(metrics, dict) or set(metrics) != set(tools) or any(
                value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))
                                       or not math.isfinite(value) or value < 0)
                for value in metrics.values()):
            raise ValueError("invalid native Top-1 DCA measurements; missing values must be null")
        if "best_detector" in case and case["best_detector"] != infer_best_detector(metrics):
            raise ValueError("stored best detector disagrees with its declared raw-DCA derivation")
    return cases


def _unavailable(mode, reason):
    if mode == "required":
        raise ValueError(f"pocket memory is required but unavailable: {reason}")
    return {"mode": mode, "status": "unavailable", "reason": reason}


def _memory_evidence(index, embedding, mode, query_case, query_group, exclude_cases, exclude_groups):
    if mode == "off":
        return {"mode": mode, "status": "disabled", "reason": "memory_disabled"}
    # Validate each supplied input even when its companion is missing: invalid
    # provenance or numbers must never silently become an ordinary fallback.
    cases = validate_index(index) if index is not None else None
    if embedding is not None:
        _protein_vector(embedding)
    if cases is None or embedding is None:
        return _unavailable(mode, "missing_index" if cases is None else "missing_query_embedding")
    query_case = _case_id(query_case)
    _name(query_group, "query group", optional=True)
    excluded_cases = {_case_id(value) for value in exclude_cases}
    excluded_groups = {_name(value, "excluded group") for value in exclude_groups}
    if None in excluded_cases:
        raise ValueError("excluded case identities cannot be null")
    for case in cases:
        if query_case == _case_id(case["case_id"]):
            if query_group is not None and query_group != case["group_id"]:
                raise ValueError("query group conflicts with its indexed case")
            query_group = case["group_id"]
    if query_case:
        excluded_cases.add(query_case)
    if query_group:
        excluded_groups.add(query_group)
    eligible = [case for case in cases
                if _case_id(case["case_id"]) not in excluded_cases and case["group_id"] not in excluded_groups
                and all(case["top1_dca_angstrom"].get(tool) is not None for tool in PAPER_TOOLS)]
    ranked = sorted(((sum(a * b for a, b in zip(embedding, case["protein"])), case)
                     for case in eligible), key=lambda item: (-item[0], _case_id(item[1]["case_id"])))
    groups = {}
    for similarity, case in ranked:
        groups.setdefault(case["group_id"], (similarity, case))
    details = {"eligible_cases": len(eligible), "eligible_groups": len(groups),
               "excluded_cases": sorted(excluded_cases), "excluded_groups": sorted(excluded_groups),
               "population_policy": "complete_five_tool_cases_after_case_and_group_exclusions",
               "neighbor_policy": "closest_case_per_group_then_cosine_descending_case_id_ascending",
               "best_detector_policy": BEST_DETECTOR_POLICY}
    if len(groups) < 5:
        return {**_unavailable(mode, "fewer_than_five_complete_protein_groups"), **details}
    neighbors = list(groups.values())[:5]
    raw = [math.exp((similarity - neighbors[0][0]) / 0.02) / math.log2(j + 2)
           for j, (similarity, _) in enumerate(neighbors)]
    weights = [5.0 * value / sum(raw) for value in raw]
    priors = {tool: sum(case["top1_dca_angstrom"][tool] < 4.0 for case in eligible) / len(eligible)
              for tool in PAPER_TOOLS}
    reliability = {tool: (5.0 * priors[tool] + sum(
        weight * (case["top1_dca_angstrom"][tool] < 4.0)
        for weight, (_, case) in zip(weights, neighbors))) / (5.0 + sum(weights)) for tool in PAPER_TOOLS}
    return {"mode": mode, "status": "used", "reason": "compatible_history_available", **details,
            "priors": priors, "reliability": reliability,
            "neighbors": [{"case_id": case["case_id"], "group_id": case["group_id"],
                           "cosine_similarity": similarity, "weight": weight,
                           "best_detector": infer_best_detector(case["top1_dca_angstrom"]),
                           "best_detector_origin": "stored_derived_label" if "best_detector" in case else "derived_from_raw_dca",
                           "top1_dca_angstrom": dict(case["top1_dca_angstrom"]),
                           "top1_labels": {tool: int(case["top1_dca_angstrom"][tool] < 4.0) for tool in PAPER_TOOLS}}
                          for weight, (similarity, case) in zip(weights, neighbors)]}


def _matching(detectors, anchor_tool, *, all_matches=False):
    anchor = detectors[anchor_tool][0]
    members = [(anchor_tool, anchor)]
    for tool, candidates in detectors.items():
        if tool == anchor_tool and not all_matches:
            continue
        options = []
        for pocket in candidates:
            if tool == anchor_tool and pocket["rank"] == 1:
                continue
            distance, iou = _agreement(anchor, pocket)
            if distance <= 5.0 and iou >= 0.30:
                if all_matches:
                    members.append((tool, pocket))
                else:
                    options.append((iou, -distance, -pocket["rank"], pocket))
        if options:
            members.append((tool, max(options, key=lambda option: option[:3])[3]))
    return members


def _tool_plan(memory):
    used = memory["status"] == "used"
    multiset = [neighbor["best_detector"] for neighbor in memory["neighbors"]] if used else []
    selected = {"deeppocket", *multiset} if used else set(TOOLS)
    execution = [tool for tool in PAPER_TOOLS if tool in selected]
    return {"strategy": MEMORY_STRATEGY if used else STRATEGY,
            "implementation_status": "deterministic_memory_policy",
            "detector_selection_policy": "neighbor_best_detector_multiset_plus_deeppocket_anchor" if used else "fixed_four_detector_fallback",
            "recommended_multiset": multiset, "execution_tools": execution,
            "anchor_forced": used and "deeppocket" not in multiset,
            "inference_notes": list(INFERENCE_NOTES) if used else [],
            "memory_used": used, "memory_mode": memory["mode"], "memory_reason": memory["reason"], "memory": memory}


def predict(payload, embedding=None, index=None, *, preset="conservative", query_case=None,
            query_group=None, exclude_cases=(), exclude_groups=(), memory_mode="auto"):
    """Pure selection from normalized candidates and optional in-memory history.

    Both presets use the same memory policy. Missing memory falls back in
    auto mode; corrupt supplied inputs raise. Candidate detectors must match
    the historical recommendation plan, including the mandatory default anchor.
    """
    _validate_options(memory_mode, preset)
    memory = _memory_evidence(index, embedding, memory_mode, query_case, query_group, exclude_cases, exclude_groups)
    return _predict_with_memory(payload, memory, preset)


def _validate_options(memory_mode, preset):
    if preset not in ("conservative", "exploratory_memory_guard"):
        raise ValueError("unknown pocket preset")
    if memory_mode not in ("auto", "off", "required"):
        raise ValueError("pocket memory mode must be auto, off or required")


def _predict_with_memory(payload, memory, preset):
    plan = _tool_plan(memory)
    active_tools = _active_tools(payload)
    if set(active_tools) != set(plan["execution_tools"]):
        raise ValueError("candidate detectors differ from the pocket execution plan; expected "
                         + ", ".join(plan["execution_tools"]))
    detectors = {}
    for tool in active_tools:
        candidates = payload["detectors"].get(tool)
        if not isinstance(candidates, list) or not 1 <= len(candidates) <= 3:
            raise ValueError(f"{tool}: one to three normalized candidates required")
        ranks = [p.get("rank") for p in candidates]
        if any(type(rank) is not int or rank not in (1, 2, 3) for rank in ranks) or len(set(ranks)) != len(ranks):
            raise ValueError(f"{tool}: unique native ranks 1--3 required")
        detectors[tool] = sorted(candidates, key=lambda pocket: pocket["rank"])
        if detectors[tool][0]["rank"] != 1:
            raise ValueError(f"{tool}: native rank 1 is required; do not relabel lower ranks")
        for candidate in candidates:
            _xyz(candidate.get("center"))
            if not candidate.get("residue_records"):
                raise ValueError(f"{tool}: point5 residues required")
            for record in candidate["residue_records"]:
                _identity(record)
            if tool in ("deeppocket", "p2rank"):
                probability = candidate.get("probability")
                if isinstance(probability, bool) or not isinstance(probability, (int, float)) or not math.isfinite(probability) or not 0 <= probability <= 1:
                    raise ValueError(f"{tool}: invalid native confidence")
    memory_used = memory["status"] == "used"
    anchor_tool, guard = "deeppocket", {}
    if memory_used and "p2rank" in detectors:
        dp, p2 = detectors["deeppocket"][0], detectors["p2rank"][0]
        p2_matches = _matching(detectors, "p2rank", all_matches=True)
        values = {"deeppocket_native_probability": dp["probability"], "p2rank_native_probability": p2["probability"],
                  "anchor_distance_A": math.dist(dp["center"], p2["center"]),
                  "reliability_advantage": memory["reliability"]["p2rank"] - memory["reliability"]["deeppocket"],
                  "nearest_cosine_similarity": memory["neighbors"][0]["cosine_similarity"],
                  "p2rank_matching_families": len({_family(tool) for tool, _ in p2_matches})}
        guard = {"deeppocket_confidence_lt_0_95": values["deeppocket_native_probability"] < 0.95,
                 "p2rank_confidence_ge_0_80": values["p2rank_native_probability"] >= 0.80,
                 "anchor_distance_gt_5A": values["anchor_distance_A"] > 5.0,
                 "reliability_advantage_gt_0_15": values["reliability_advantage"] > 0.15,
                 "nearest_similarity_ge_0_88": values["nearest_cosine_similarity"] >= 0.88,
                 "p2rank_family_support_ge_2": values["p2rank_matching_families"] >= 2}
        memory["guard_values"] = values
        if all(guard.values()):
            anchor_tool = "p2rank"
    elif memory_used:
        memory["guard_values"] = {"p2rank_available": False, "reason": "p2rank_not_recommended"}
    anchor = detectors[anchor_tool][0]
    members = _matching(detectors, anchor_tool, all_matches=memory_used)
    byfamily = defaultdict(list)
    for tool, pocket in members:
        value = (memory["reliability"][tool] + 0.05) ** 2 if memory_used else 1.0
        byfamily[_family(tool)].append((tool, pocket, value / pocket["rank"]))
    numerator, denominator, audit = [0.0, 0.0, 0.0], 0.0, []
    for family, entries in byfamily.items():
        maximum, total = max(v for _, _, v in entries), sum(v for _, _, v in entries)
        for tool, pocket, v in entries:
            effective = v * maximum / total
            for axis in range(3):
                numerator[axis] += effective * pocket["center"][axis]
            denominator += effective
            distance, iou = _agreement(anchor, pocket)
            member = {"tool": tool, "rank": pocket["rank"], "family": family,
                      "effective_weight": effective, "raw_weight": v, "native_rank_weight": 1.0 / pocket["rank"],
                      "center_distance_A": distance, "residue_IoU": iou}
            if "probability" in pocket:
                member["native_probability"] = pocket["probability"]
            audit.append(member)
    parameters = dict(PARAMETERS)
    parameters["matching_policy"] = "one_per_detector_max_IoU_then_min_distance_then_native_rank"
    if memory_used:
        parameters.update(weight_policy="squared_reliability_plus_0_05_over_native_rank", neighbor_groups=5,
                          neighbor_weight_sum=5.0, temperature=0.02, prior_strength=5.0, success_rule="native_top1_DCA_lt_4A",
                          matching_policy="all_matching_native_top3_candidates",
                          best_detector_policy=BEST_DETECTOR_POLICY)
    return {**plan, "success": True, "status": "selected", "preset": preset,
            "target": payload.get("target"), "center": [value / denominator for value in numerator],
            "residue_records": anchor["residue_records"], "anchor_tool": anchor_tool,
            "anchor_rank": 1, "guard_switched": anchor_tool != "deeppocket", "guard_conditions": guard, "members": audit,
            "independent_support": len(byfamily), "parameters": parameters,
            "active_tools": list(active_tools), "unselected_tools": sorted(set(PAPER_TOOLS) - set(active_tools))}


def _query_identity(metadata):
    if not isinstance(metadata, dict) or not isinstance(metadata.get("source_record"), dict):
        raise ValueError("query metadata must include its original source_record")
    source = metadata["source_record"]
    if source.get("task_type") != "pocket_prediction":
        raise ValueError("query embedding metadata must be for pocket_prediction")
    case = _case_id(source.get("case_id"))
    if source.get("dataset") and source.get("target"):
        derived = _case_id(f"{source['dataset']}/{source['target']}")
        if case and case != derived:
            raise ValueError("query case identity disagrees with dataset/target metadata")
        case = derived
    return case, _name(source.get("group_id"), "query group", optional=True)


def _load_memory(*, receptor_file, index_path=None, query_npz=None, query_metadata=None,
                 memory_mode="auto", preset="conservative", query_case=None, query_group=None,
                 exclude_cases=(), exclude_groups=()):
    """Shared pre-detection and post-detection history/provenance loader."""
    _validate_options(memory_mode, preset)
    options = dict(query_case=query_case, query_group=query_group,
                   exclude_cases=exclude_cases, exclude_groups=exclude_groups)
    if memory_mode == "off":
        return {"mode": memory_mode, "status": "disabled", "reason": "memory_disabled"}
    paths = {"index": Path(index_path).expanduser() if index_path is not None else default_index_path(),
             "query_npz": Path(query_npz).expanduser() if query_npz is not None else None,
             "query_metadata": Path(query_metadata).expanduser() if query_metadata is not None else None}
    missing = []
    for name, path in paths.items():
        if path is None or not path.exists():
            missing.append(name)
        elif not path.is_file():
            raise ValueError(f"pocket {name} must be a file")
    index = _read_json(paths["index"]) if "index" not in missing else None
    if index is not None:
        validate_index(index)
    selector = _sibling("docking_selector")
    embedding, metadata, query_provenance = None, None, {}
    if "query_npz" not in missing:
        shape, vector = selector.read_npz_float32(paths["query_npz"], "protein_embedding_normalized")
        if shape != (1280,):
            raise ValueError("query protein embedding must have shape (1280,)")
        embedding = list(vector)
        _protein_vector(embedding)
        query_provenance.update(npz=str(paths["query_npz"].resolve()), embedding_sha256=_sha(paths["query_npz"]))
    if "query_metadata" not in missing:
        metadata = _read_json(paths["query_metadata"])
        case, group = _query_identity(metadata)
        encoder = selector.encoder_id_from_metadata(metadata)
        if index is not None and encoder != index["encoder_compatibility_id"]:
            raise ValueError("query encoder is incompatible with pocket memory")
        receptor_sha = _sha(receptor_file)
        if metadata.get("source_input_sha256") != receptor_sha:
            raise ValueError("query receptor SHA-256 does not match embedding metadata")
        source_sha = metadata["source_record"].get("protein_sha256")
        if source_sha is not None and source_sha != receptor_sha:
            raise ValueError("query source_record protein SHA-256 conflicts with its receptor binding")
        claimed_npz = metadata.get("embedding_sha256")
        if not isinstance(claimed_npz, str) or len(claimed_npz) != 64 or any(c not in "0123456789abcdef" for c in claimed_npz):
            raise ValueError("query metadata is missing a valid embedding SHA-256")
        if embedding is not None and claimed_npz != query_provenance["embedding_sha256"]:
            raise ValueError("query embedding SHA-256 does not match its metadata")
        explicit_case = _case_id(query_case)
        _name(query_group, "explicit query group", optional=True)
        if explicit_case and case and explicit_case != case:
            raise ValueError("explicit query case conflicts with encoder metadata")
        if query_group and group and query_group != group:
            raise ValueError("explicit query group conflicts with encoder metadata")
        options.update(query_case=explicit_case or case, query_group=query_group or group)
        query_provenance.update(metadata=str(paths["query_metadata"].resolve()),
                                metadata_sha256=_sha(paths["query_metadata"]),
                                encoder_compatibility_id=encoder, receptor_sha256=receptor_sha,
                                case_id=options["query_case"], group_id=options["query_group"])
    if missing:
        memory = _unavailable(memory_mode, "missing_assets:" + ",".join(missing))
    else:
        memory = _memory_evidence(index, embedding, memory_mode, **options)
    memory["requested_assets"] = {name: str(path.resolve()) if path else None for name, path in paths.items()}
    if index is not None:
        memory["index"] = {"path": str(paths["index"].resolve()), "sha256": _sha(paths["index"]),
                           "source": index["source"], "encoder_compatibility_id": index["encoder_compatibility_id"],
                           "best_detector_policy": BEST_DETECTOR_POLICY}
    if query_provenance:
        memory["query"] = query_provenance
    return memory


def plan_tools(*, receptor_file, index_path=None, query_npz=None, query_metadata=None,
               memory_mode="auto", preset="conservative", query_case=None, query_group=None,
               exclude_cases=(), exclude_groups=()):
    """Retrieve five profile recommendations before running any detector.

    The five-entry multiset preserves repeated historical recommendations;
    unique tools execute once, with DeepPocket additionally supplying the
    default anchor. This function never reads current detector outputs.
    """
    memory = _load_memory(receptor_file=receptor_file, index_path=index_path, query_npz=query_npz,
                          query_metadata=query_metadata, memory_mode=memory_mode, preset=preset,
                          query_case=query_case, query_group=query_group,
                          exclude_cases=exclude_cases, exclude_groups=exclude_groups)
    return {**_tool_plan(memory), "success": True, "status": "planned", "preset": preset,
            "receptor_sha256": None if memory_mode == "off" else _sha(receptor_file)}


def decide(payload, *, receptor_file, index_path=None, query_npz=None, query_metadata=None,
           memory_mode="auto", preset="conservative", query_case=None, query_group=None,
           exclude_cases=(), exclude_groups=()):
    """Verify current candidates, then apply the same history plan as routing.

    Invalid provenance fails closed; only absent assets/insufficient complete
    history permit auto fallback. Off mode never opens historical inputs.
    """
    _validate_candidates(payload, receptor_file)
    memory = _load_memory(receptor_file=receptor_file, index_path=index_path, query_npz=query_npz,
                          query_metadata=query_metadata, memory_mode=memory_mode, preset=preset,
                          query_case=query_case, query_group=query_group,
                          exclude_cases=exclude_cases, exclude_groups=exclude_groups)
    result = _predict_with_memory(payload, memory, preset)
    result["receptor_sha256"] = _sha(receptor_file)
    return result


def make_box(decision, receptor_file, *, padding=4.0, min_size=12.0):
    """Explicit downstream policy, separate from the benchmark's pocket metrics."""
    if not math.isfinite(padding) or padding < 0 or not math.isfinite(min_size) or min_size <= 0:
        raise ValueError("box padding must be finite/nonnegative and minimum size positive")
    atoms, _ = _protein_atoms(receptor_file)
    center = _xyz(decision["center"])
    points = [point for record in decision["residue_records"] for point in atoms[_identity(record)]]
    if not points:
        raise ValueError("anchor residues have no current receptor heavy atoms")
    size = [max(min_size, 2 * (max(abs(point[d] - center[d]) for point in points) + padding)) for d in range(3)]
    return {"schema_version": 1, "center": center, "size": size, "unit": "angstrom",
            "source": decision["strategy"], "geometry_source": "anchor_point5_residue_extent_about_fused_center",
            "padding": padding, "min_size": min_size, "receptor_sha256": _sha(receptor_file),
            "anchor_tool": decision["anchor_tool"],
            "note": "Explicit downstream box policy; not a box produced or validated by the original pocket benchmark."}


def select_pocket(*, receptor_file, candidates_path, index_path=None, query_npz=None,
                  query_metadata=None, preset="conservative", query_case=None, query_group=None,
                  exclude_cases=(), exclude_groups=(), export_box=None, box_padding=4.0,
                  min_box_size=12.0, memory_mode="auto"):
    """File-based selection and optional box export through the shared policy."""
    protected_index = Path(index_path).expanduser() if index_path is not None else default_index_path()
    if export_box is not None and Path(export_box).resolve() in {
            Path(p).resolve() for p in (protected_index, query_npz, query_metadata, receptor_file, candidates_path)
            if p is not None}:
        raise ValueError("exported box must not overwrite an input")
    payload = _read_json(candidates_path)
    result = decide(payload, receptor_file=receptor_file, index_path=index_path, query_npz=query_npz,
                    query_metadata=query_metadata, preset=preset, query_case=query_case, query_group=query_group,
                    exclude_cases=exclude_cases, exclude_groups=exclude_groups, memory_mode=memory_mode)
    result.update(receptor_sha256=_sha(receptor_file),
                  candidates={"path": str(Path(candidates_path).resolve()), "sha256": _sha(candidates_path)})
    if export_box is not None:
        box = make_box(result, receptor_file, padding=box_padding, min_size=min_box_size)
        _write_json(export_box, box)
        result["exported_box"] = {"path": str(Path(export_box).resolve()), "sha256": _sha(export_box), **box}
    return result

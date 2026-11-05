"""Fixed five-metric generation memory; no query labels or parameter search."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
from statistics import pstdev
import tempfile
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from deployment_config import deployment_path


_spec = importlib.util.spec_from_file_location(
    "generation_ligand_features", Path(__file__).with_name("docking_selector.py"))
_features = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_features)

DEFAULT_BUNDLE = deployment_path("${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-generation-memory/memory")
DEFAULT_PROPERTY_WINNERS = DEFAULT_BUNDLE.parent / "evaluation_v2/target_best_tools.csv"
DEFAULT_GROUPS = DEFAULT_BUNDLE.parent / "audit/similarity_groups.csv"
FORMAT_VERSION = "generation_selector/1"
POLICY_VERSION = "ligand_protein_generation_reweight_v1"
FEATURE_VERSION = _features.FEATURE_VERSION
TOOLS = ("difflinker", "diffsbdd", "lingo3dmol", "molcraft", "pocket2mol", "resgen", "tamgen", "targetdiff")
METRICS = ("median_qed", "median_normalized_sa", "molecular_diversity_n50", "scaffold_entropy_n50", "gate_pass_rate_unique")
PARAMETERS = {"profile": "property_focused", "k": 3, "protein_beta": 0.10,
              "prior_strength": 10.0, "switch_margin": 0.03, "default_tool": "pocket2mol"}


def _encoder_id(metadata: dict) -> str:
    # The generation release uses a narrower hash than the docking release.
    encoder = metadata.get("encoder", {})
    keys = ("model_files_sha256", "window_size", "window_overlap", "window_merge", "pooling", "context")
    if any(key not in encoder for key in keys):
        raise ValueError("incomplete generation encoder metadata")
    payload = {key: encoder[key] for key in keys}
    return "encoder_compatibility:" + hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def extract_query_features(ligand_sdf: Path, embeddings_npz: Path) -> dict:
    """Use the frozen graph schema, including the release's V3000 Fe complex."""
    path = Path(ligand_sdf)
    return _extract_features_text(path.read_text(encoding="utf-8"), embeddings_npz, str(path))


def _extract_features_text(text: str, embeddings_npz: Path, source_label: str) -> dict:
    records = text.split("$$$$")
    if len(records) > 2 or (len(records) == 2 and records[1].strip()):
        raise ValueError("generation memory requires one legal reference ligand, not a library")
    lines = text.splitlines()
    if len(lines) < 4 or "V3000" not in lines[3]:
        return _features.extract_query_features_text(text, embeddings_npz, source_label)
    # RDKit is only needed for V3000. Preserve the research adapter's order-9
    # dative -> single graph convention; this view is never a pipeline ligand.
    from rdkit import Chem
    mol = Chem.MolFromMolBlock(records[0], sanitize=False, removeHs=False)
    if mol is None or not mol.GetNumConformers():
        raise ValueError("invalid V3000 reference ligand")
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DATIVE:
            bond.SetBondType(Chem.BondType.SINGLE)
    view = Chem.MolToMolBlock(mol, kekulize=False, forceV3000=False) + "\n$$$$\n"
    result = _features.extract_query_features_text(view, embeddings_npz, source_label)
    # V2000 serialization rounds coordinates; compute shape from original data.
    conformer = mol.GetConformer()
    xyz = [tuple(conformer.GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetSymbol() not in ("H", "D")]
    if not xyz:
        raise ValueError("reference ligand has no heavy atoms")
    center = tuple(sum(point[d] for point in xyz) / len(xyz) for d in range(3))
    result["ligand_desc"][-2:] = [
        math.sqrt(sum(math.dist(point, center) ** 2 for point in xyz) / len(xyz)),
        max((math.dist(a, b) for i, a in enumerate(xyz) for b in xyz[i + 1:]), default=0.0),
    ]
    _features._validate_features(result)
    return result


def validate_query_metadata(metadata_path: Path, embeddings_npz: Path, index: dict,
                            receptor_file: Path | None = None, ligand_sdf: Path | None = None) -> dict:
    """Bind current inputs using generation query.json or hashed encoder metadata."""
    path = Path(metadata_path)
    query = json.loads(path.read_text(encoding="utf-8"))
    if query.get("task_type") == "molecule_generation":
        identity = query.get("identity", {})
        case_id = identity.get("case_id")
        group_id = identity.get("group_id")
        files, hashes = query["files"], query["file_hashes"]
        metadata_file = (path.parent / files["embedding_metadata"]).resolve()
        if path.parent.resolve() not in metadata_file.parents:
            raise ValueError("query encoder metadata must remain inside its input directory")
        if _features._digest(metadata_file) != hashes.get("embedding_metadata"):
            raise ValueError("query encoder metadata SHA-256 mismatch")
        metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
        expected_embedding, expected_receptor = hashes.get("embedding"), hashes.get("protein")
        if ligand_sdf is not None and _features._digest(Path(ligand_sdf)) != hashes.get("reference_ligand"):
            raise ValueError("query reference ligand SHA-256 mismatch")
    else:
        metadata = query
        source = metadata.get("source_record", {})
        if source.get("task_type") not in (None, "molecule_generation"):
            raise ValueError("query metadata is not for molecule_generation")
        case_id = source.get("case_id")
        if not case_id and source.get("dataset") and source.get("target"):
            case_id = f"{source['dataset']}/{source['target']}"
        group_id = source.get("group_id")
        expected_embedding = metadata.get("embedding_sha256")
        expected_receptor = metadata.get("source_input_sha256") or source.get("protein_sha256")
    source = metadata.get("source_record", {})
    if source.get("task_type") not in (None, "molecule_generation"):
        raise ValueError("query encoder metadata is not for molecule_generation")
    source_hash = metadata.get("source_input_sha256") or source.get("protein_sha256")
    if source_hash is not None and source_hash != expected_receptor:
        raise ValueError("query receptor provenance differs from encoder source")
    source_case = source.get("case_id")
    if not source_case and source.get("dataset") and source.get("target"):
        source_case = f"{source['dataset']}/{source['target']}"
    if source_case and case_id and source_case != case_id:
        raise ValueError("query case differs from encoder source identity")
    for value in (case_id, group_id):
        if value is not None and (not isinstance(value, str) or not value):
            raise ValueError("invalid query case or group identity")
    compatibility = _encoder_id(metadata)
    if compatibility != index.get("encoder_compatibility_id"):
        raise ValueError("query encoder is incompatible with generation memory")
    contract_id = _features.encoder_id_from_metadata(metadata)
    if index.get("encoder_contract_id", contract_id) != contract_id:
        raise ValueError("query encoder dtype, pLDDT mask or preprocessing differs from memory")
    embedding_hash = _features._digest(Path(embeddings_npz))
    if embedding_hash != expected_embedding:
        raise ValueError("query embedding SHA-256 does not match metadata")
    if receptor_file is not None and _features._digest(Path(receptor_file)) != expected_receptor:
        raise ValueError("query receptor SHA-256 does not match metadata")
    return {"case_id": case_id, "group_id": group_id, "encoder_compatibility_id": compatibility,
            "encoder_contract_id": contract_id, "embedding_sha256": embedding_hash, "source_input_sha256": expected_receptor,
            "receptor_sha256_checked": receptor_file is not None}


def _memory_rows(path: Path, allowed: set[str], column: str) -> dict:
    result = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            case_id = f"{row['dataset']}/{row['target']}"
            # Check membership and split BEFORE accessing labels or group values.
            if case_id not in allowed:
                continue
            if row.get("split") != "memory" or case_id in result:
                raise ValueError(f"invalid memory-only source row for {case_id}")
            result[case_id] = row[column]
    if set(result) != allowed:
        raise ValueError(f"incomplete memory source {path}: missing {len(allowed - set(result))} cases")
    return result


def _winner_labels(winners) -> dict:
    if not isinstance(winners, list) or not winners or any(tool not in TOOLS for tool in winners):
        raise ValueError("missing or invalid published metric winners")
    return {tool: tool in winners for tool in TOOLS}


def _memory_features(directory: Path) -> dict:
    # The frozen research extractor used record 0 for these historical assets.
    # This compatibility rule never applies to a new user library/query.
    text = (directory / "ligand.sdf").read_text(encoding="utf-8")
    if text.count("$$$$") <= 1:
        return extract_query_features(directory / "ligand.sdf", directory / "embeddings.npz")
    return _extract_features_text(text.split("$$$$", 1)[0] + "$$$$\n", directory / "embeddings.npz",
                                  str(directory / "ligand.sdf") + " record 0")


def build_index(memory_bundle: Path, output: Path,
                property_winners_csv: Path = DEFAULT_PROPERTY_WINNERS,
                groups_csv: Path = DEFAULT_GROUPS) -> dict:
    """Offline only: allowlist frozen memory IDs before joining mixed source CSVs."""
    bundle, output = Path(memory_bundle).resolve(), Path(output).resolve()
    if output == bundle.parent or bundle.parent in output.parents:
        raise ValueError("write the generation index outside the frozen release")
    records_path = bundle / "memory/agent_memory.jsonl"
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records_hash = _features._digest(records_path)
    if manifest.get("hashes", {}).get("memory/agent_memory.jsonl") != records_hash:
        raise ValueError("generation memory records SHA-256 mismatch")
    records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    allowed = {record["identity"]["case_id"] for record in records}
    if len(allowed) != len(records) or len(records) < 3:
        raise ValueError("generation memory requires at least three unique cases")
    properties = _memory_rows(property_winners_csv, allowed, "best_gate_pass_rate_tools")
    groups = _memory_rows(groups_csv, allowed, "group_id")
    inputs_root = (bundle / "inputs/molecule_generation").resolve()
    cases, encoder_ids, encoder_contracts = [], set(), set()
    input_digest = hashlib.sha256()
    for record in records:
        case_id = record["identity"]["case_id"]
        directory = (inputs_root / case_id).resolve()
        if inputs_root not in directory.parents:
            raise ValueError("generation case escapes memory input directory")
        query_path = directory / "query.json"
        query = json.loads(query_path.read_text(encoding="utf-8"))
        if query.get("split") != "memory" or query.get("identity", {}).get("case_id") != case_id:
            raise ValueError("generation index accepts matching memory queries only")
        encoder = record["encoder_compatibility_id"]
        binding = validate_query_metadata(query_path, directory / "embeddings.npz",
                    {"encoder_compatibility_id": encoder}, directory / "protein.pdb", directory / "ligand.sdf")
        if record["inputs"]["embedding_sha256"] != binding["embedding_sha256"]:
            raise ValueError("memory card embedding hash mismatch")
        if record["inputs"].get("protein_sha256", binding["source_input_sha256"]) != binding["source_input_sha256"]:
            raise ValueError("memory card receptor hash mismatch")
        if record["inputs"]["reference_ligand_sha256"] != _features._digest(directory / "ligand.sdf"):
            raise ValueError("memory card ligand hash mismatch")
        encoder_ids.add(encoder)
        encoder_contracts.add(binding["encoder_contract_id"])
        labels = {metric: _winner_labels(record["metric_winners"][metric]["best_tools"]) for metric in METRICS[:-1]}
        labels[METRICS[-1]] = _winner_labels(properties[case_id].split(";"))
        cases.append({"case_id": case_id, "group_id": groups[case_id], "labels": labels,
                      **_memory_features(directory)})
        for name in ("query.json", "embedding_metadata.json", "embeddings.npz", "ligand.sdf", "protein.pdb"):
            input_digest.update(f"{case_id}/{name}\0{_features._digest(directory / name)}\n".encode())
    if len(encoder_ids) != 1 or len(encoder_contracts) != 1:
        raise ValueError("generation memory contains incompatible encoders")
    # Preserve the research selector's target-id order, including deterministic ties.
    cases.sort(key=lambda case: (case["case_id"].split("/", 1)[-1], case["case_id"]))
    source = {"split": "memory", "bundle": str(bundle), "records_sha256": records_hash,
              "feature_inputs_sha256": input_digest.hexdigest(), "manifest_sha256": _features._digest(manifest_path),
              "property_winners_csv": str(property_winners_csv), "property_winners_sha256": _features._digest(Path(property_winners_csv)),
              "groups_csv": str(groups_csv), "groups_sha256": _features._digest(Path(groups_csv)),
              "historical_ligand_record": 0,
              "join_policy": "memory record allowlist AND split=memory; only selected IDs and fields retained"}
    index = {"schema_version": FORMAT_VERSION, "feature_version": FEATURE_VERSION,
             "policy_version": POLICY_VERSION, "parameters": PARAMETERS, "tools": list(TOOLS),
             "metrics": list(METRICS), "source": source, "cases": cases,
             "encoder_compatibility_id": next(iter(encoder_ids)),
             "encoder_contract_id": next(iter(encoder_contracts)),
             "memory_version": hashlib.sha256(json.dumps(source, sort_keys=True).encode()).hexdigest()}
    _validate_index(index)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(index, handle, separators=(",", ":"), allow_nan=False)
            handle.write("\n")
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return {"success": True, "status": "built", "path": str(output), "case_count": len(cases),
            "memory_version": index["memory_version"], "policy_version": POLICY_VERSION}


def _validate_index(index: dict) -> None:
    if not isinstance(index, dict) or index.get("schema_version") != FORMAT_VERSION:
        raise ValueError("incompatible generation index schema")
    if (index.get("feature_version") != FEATURE_VERSION or index.get("policy_version") != POLICY_VERSION
            or index.get("parameters") != PARAMETERS or index.get("tools") != list(TOOLS)
            or index.get("metrics") != list(METRICS) or index.get("source", {}).get("split") != "memory"):
        raise ValueError("incompatible generation index features, policy or split")
    if not index.get("memory_version") or not index.get("encoder_compatibility_id") or not index.get("encoder_contract_id"):
        raise ValueError("generation index requires memory and encoder versions")
    seen = set()
    if not isinstance(index.get("cases"), list) or not index["cases"]:
        raise ValueError("empty generation index")
    for case in index["cases"]:
        _features._validate_features(case)
        case_id, group_id = case.get("case_id"), case.get("group_id")
        if (not isinstance(case_id, str) or not case_id or case_id in seen
                or not isinstance(group_id, str) or not group_id):
            raise ValueError("invalid memory case or protein group")
        seen.add(case_id)
        labels = case.get("labels", {})
        if set(labels) != set(METRICS):
            raise ValueError("missing generation metric labels")
        for metric in METRICS:
            if (set(labels[metric]) != set(TOOLS) or not any(labels[metric].values())
                    or any(type(label) is not bool for label in labels[metric].values())):
                raise ValueError("generation labels must be published per-metric winner booleans")


def load_index(path: Path) -> dict:
    index = json.loads(Path(path).read_text(encoding="utf-8"))
    _validate_index(index)
    return index


def fallback_decision(reason: str, available_tools=None, priority_metrics=()) -> dict:
    available = set(TOOLS if available_tools is None else available_tools)
    selected = "pocket2mol" if "pocket2mol" in available else None
    return {"success": selected is not None, "status": "fallback" if selected else "unavailable",
            "policy_version": POLICY_VERSION, "feature_version": FEATURE_VERSION,
            "memory_version": None, "fallback_reason": reason, "neighbors": [],
            "priority_metrics": list(priority_metrics),
            "required_metrics": [metric for metric in METRICS if metric in priority_metrics] or ["gate_pass_rate_unique"],
            "recommendations": {metric: {"selected_tool": selected, "proposed_tool": None,
                "status": "fallback" if selected else "unavailable",
                "scores": {}, "fallback_reason": reason} for metric in METRICS},
            "default_recommendation": selected}


def select_generation_tools(features: dict, index: dict, priority_metrics=(),
                            available_tools=None, exclude_cases=(), exclude_groups=()) -> dict:
    """Paper Eq. 4–5 per metric; priority selects objectives, never hyperparameters."""
    priorities = set(priority_metrics)
    if priorities - set(METRICS):
        raise ValueError("unknown generation priority metric: " + ", ".join(sorted(priorities - set(METRICS))))
    _validate_index(index)
    _features._validate_features(features)
    available = set(TOOLS if available_tools is None else available_tools)
    excluded_cases, excluded_groups = set(exclude_cases), set(exclude_groups)
    excluded_groups.update(case["group_id"] for case in index["cases"] if case["case_id"] in excluded_cases)
    cases = [case for case in index["cases"] if case["case_id"] not in excluded_cases and case["group_id"] not in excluded_groups]
    if len(cases) < 3:
        return fallback_decision("fewer than three eligible memory neighbors", available, priority_metrics)
    scales = [pstdev(case["ligand_desc"][j] for case in cases) for j in range(18)]
    scales = [scale if scale > 0.0 else 1.0 for scale in scales]
    query_fp = set(features["ligand_fp"])
    protein = [sum(a * b for a, b in zip(features["protein"], case["protein"])) for case in cases]
    protein_order = sorted(range(len(cases)), key=lambda i: (protein[i], -i))
    percentile = {i: rank / max(len(cases) - 1, 1) for rank, i in enumerate(protein_order)}
    similarities = []
    for case in cases:
        fp = set(case["ligand_fp"])
        graph = len(query_fp & fp) / len(query_fp | fp)
        squared = sum(((a - b) / scale) ** 2 for a, b, scale in zip(features["ligand_desc"], case["ligand_desc"], scales))
        descriptor = math.exp(-0.5 * squared / 18)
        similarities.append((0.8 * graph + 0.2 * descriptor, graph, descriptor))
    indices = sorted(range(len(cases)), key=lambda i: (-similarities[i][0], i))[:3]
    weights = [math.exp(0.10 * (percentile[i] - 0.5)) / math.log2(rank + 2) for rank, i in enumerate(indices)]
    neighbors = [{"case_id": cases[i]["case_id"], "group_id": cases[i]["group_id"],
                  "ligand_similarity": similarities[i][0], "ligand_fingerprint_tanimoto": similarities[i][1],
                  "ligand_descriptor_rbf": similarities[i][2], "protein_cosine": protein[i],
                  "protein_percentile": percentile[i], "weight": weight}
                 for i, weight in zip(indices, weights)]
    recommendations = {}
    for metric in METRICS:
        alpha, margin = PARAMETERS["prior_strength"], PARAMETERS["switch_margin"]
        priors = {tool: sum(case["labels"][metric][tool] for case in cases) / len(cases) for tool in TOOLS}
        scores = {tool: (alpha * priors[tool] + sum(weight * cases[i]["labels"][metric][tool]
                    for i, weight in zip(indices, weights))) / (alpha + sum(weights)) for tool in TOOLS}
        best = max(TOOLS, key=lambda tool: (scores[tool], tool == "pocket2mol", tool))
        proposed = best if scores[best] - scores["pocket2mol"] > margin else "pocket2mol"
        selected = proposed if proposed in available else ("pocket2mol" if "pocket2mol" in available else None)
        recommendations[metric] = {"selected_tool": selected, "proposed_tool": proposed,
            "status": "unavailable" if selected is None else "selected" if selected == proposed else "fallback",
            "best_scoring_tool": best, "scores": scores, "tool_priors": priors,
            "parameters": {"prior_strength": alpha, "switch_margin": margin},
            "fallback_reason": None if selected == proposed else f"recommended tool {proposed} unavailable"}
    required_metrics = [metric for metric in METRICS if metric in priorities] or ["gate_pass_rate_unique"]
    required = [recommendations[metric] for metric in required_metrics]
    success = all(row["selected_tool"] is not None for row in required)
    return {"success": success,
            "status": ("unavailable" if not success else "selected"
                       if all(row["selected_tool"] == row["proposed_tool"] for row in required) else "fallback"),
            "policy_version": POLICY_VERSION, "feature_version": FEATURE_VERSION,
            "memory_version": index["memory_version"], "parameters": PARAMETERS,
            "priority_metrics": [metric for metric in METRICS if metric in priorities],
            "required_metrics": required_metrics,
            "candidate_count": len(cases), "neighbors": neighbors, "recommendations": recommendations,
            "default_recommendation": recommendations["gate_pass_rate_unique"]["selected_tool"],
            "excluded_cases": sorted(excluded_cases), "excluded_groups": sorted(excluded_groups)}

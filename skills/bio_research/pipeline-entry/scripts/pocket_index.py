"""Build protein-only pocket reliability indexes from frozen memory cards.

This offline adapter reads published native Top-1 DCA measurements and stores
an explicitly inferred best detector (minimum raw DCA, stable tool-ID ties).
It never treats inclusive ``top1_dca_le_4`` labels as strict success labels,
modifies the release, or promotes new experimental observations.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tempfile


SCHEMA_VERSION = 1
INDEX_TYPE = "pocket_protein_reliability"
POLICY_VERSION = "protein_memory_local_fusion_v1"
TOOLS = ("deeppocket", "fpocket", "p2rank", "pykvfinder")
OPTIONAL_TOOLS = ("autosite",)


def _features_module():
    path = Path(__file__).with_name("docking_selector.py")
    spec = importlib.util.spec_from_file_location("_pocket_index_features", path)
    module = importlib.util.module_from_spec(spec)
    old = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = old
    return module


def _policy_module():
    path = Path(__file__).with_name("pocket_memory.py")
    spec = importlib.util.spec_from_file_location("_pocket_index_policy", path)
    module = importlib.util.module_from_spec(spec)
    old = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = old
    return module


def _sha256(value, name: str) -> str:
    if (not isinstance(value, str) or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)):
        raise ValueError(f"missing or invalid {name} SHA-256")
    return value


def _case_identity(source: dict) -> str | None:
    case_id = source.get("case_id")
    if not case_id and source.get("dataset") and source.get("target"):
        case_id = f"{source['dataset']}/{source['target']}"
    if not isinstance(case_id, str) or not case_id:
        return None
    if "/" in case_id:
        dataset, target = case_id.split("/", 1)
        aliases = {"astex": "astex_diverse", "posebusters": "posebusters_benchmark"}
        case_id = f"{aliases.get(dataset, dataset)}/{target}"
    return case_id


def _published_dca(record: dict) -> tuple[dict, set[str], set[str]]:
    """Read only selected observations; missing/invalid assessments stay null."""
    known = set(TOOLS + OPTIONAL_TOOLS)
    selections, tool_observations, releases = {}, {}, set()
    for item in record.get("published_selections", []):
        tool, observation_id = item.get("tool_id"), item.get("observation_id")
        if tool not in known:
            continue
        if not isinstance(observation_id, str) or not observation_id:
            raise ValueError("published pocket selection has no observation identity")
        if observation_id in selections or tool in tool_observations:
            raise ValueError(f"duplicate or ambiguous published pocket observation for {tool}")
        selections[observation_id] = tool
        tool_observations[tool] = observation_id
        release = item.get("release_id")
        if release is not None:
            if not isinstance(release, str) or not release:
                raise ValueError("invalid pocket publication release identity")
            releases.add(release)
    observations = {}
    for item in record.get("tool_observations", []):
        observation_id = item.get("observation_id")
        if observation_id not in selections:
            continue
        if observation_id in observations:
            raise ValueError(f"duplicate selected pocket observation: {observation_id}")
        if item.get("action", {}).get("tool_id") != selections[observation_id]:
            raise ValueError("published selection and observed pocket tool disagree")
        observations[observation_id] = item
    if set(observations) != set(selections):
        raise ValueError("published pocket selection refers to a missing observation")

    values, seen = {tool: None for tool in tool_observations}, set()
    for assessment in record.get("assessments", []):
        observation_id = assessment.get("observation_id")
        if (observation_id not in selections
                or assessment.get("kind") != "pocket_target_summary"):
            continue
        metric = assessment.get("metrics", {}).get("top1_dca_angstrom")
        if metric is None:
            continue
        if observation_id in seen:
            raise ValueError(f"duplicate native Top-1 DCA assessment: {observation_id}")
        seen.add(observation_id)
        if not isinstance(metric, dict):
            raise ValueError("native Top-1 DCA must preserve its value and assessment status")
        if metric.get("status") != "valid":
            continue
        value = metric.get("value")
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or value < 0):
            raise ValueError("valid native Top-1 DCA must be finite and nonnegative")
        unit = metric.get("unit")
        if unit is not None and unit not in ("angstrom", "angstroms", "A", "Å", "Angstrom"):
            raise ValueError(f"unsupported native Top-1 DCA unit: {unit!r}")
        values[selections[observation_id]] = float(value)
    return values, set(tool_observations), releases


def _protein_features(record: dict, case_dir: Path, case_id: str,
                      group_id: str, helpers) -> tuple[list[float], str, dict]:
    paths = {name: case_dir / name for name in ("embeddings.npz", "metadata.json")}
    for name, path in paths.items():
        if not path.is_file():
            raise ValueError(
                f"missing pocket input {name} for {case_id}: expected {path}; "
                "restore the bundle's inputs/pocket_prediction/<case_id> layout"
            )
        if case_dir not in path.resolve().parents:
            raise ValueError(f"pocket input escapes its case directory: {path}")
    metadata = json.loads(paths["metadata.json"].read_text(encoding="utf-8"))
    source = metadata.get("source_record")
    if (not isinstance(source, dict) or source.get("task_type") != "pocket_prediction"
            or _case_identity(source) != case_id):
        raise ValueError("pocket embedding metadata has a different task or case identity")
    if source.get("group_id") is not None and source["group_id"] != group_id:
        raise ValueError("pocket embedding metadata and memory protein group disagree")
    encoder_id = helpers.encoder_id_from_metadata(metadata)
    embedding_hash = helpers._digest(paths["embeddings.npz"])
    if metadata.get("embedding_sha256") != embedding_hash:
        raise ValueError("pocket embedding SHA-256 does not match its metadata")
    source_hash = _sha256(metadata.get("source_input_sha256"), "pocket source receptor")
    if source.get("protein_sha256") is not None and source["protein_sha256"] != source_hash:
        raise ValueError("pocket source receptor hashes disagree")
    matched = [item for item in record.get("embeddings", [])
               if item.get("artifact_ref", {}).get("artifact_id") == f"sha256:{embedding_hash}"
               and item.get("source_input_sha256") == source_hash
               and item.get("encoder", {}).get("compatibility_id") == encoder_id]
    if not matched:
        raise ValueError("pocket embedding, source receptor or encoder differs from memory card")
    shape, vector = helpers.read_npz_float32(paths["embeddings.npz"], "protein_embedding_normalized")
    if shape not in ((1280,), (1, 1280)):
        raise ValueError(f"unexpected pocket protein embedding shape: {shape}")
    protein = helpers._finite_vector(vector, 1280, "pocket protein embedding")
    if abs(math.sqrt(sum(value * value for value in protein)) - 1.0) > 1e-5:
        raise ValueError("pocket protein embedding is not L2 normalized")
    return protein, encoder_id, {name: helpers._digest(path) for name, path in paths.items()}


def build_index(memory_bundle: Path, output: Path) -> dict:
    """Derive a separate index from declared pocket_prediction memory records."""
    bundle, output = Path(memory_bundle).expanduser().resolve(), Path(output).expanduser().resolve()
    if output == bundle or bundle in output.parents:
        raise ValueError("write the pocket index outside the frozen memory bundle")
    records_path = bundle / "memory/agent_memory.jsonl"
    if not records_path.is_file():
        raise ValueError(f"missing frozen memory records: {records_path}")
    if bundle not in records_path.resolve().parents:
        raise ValueError("memory records escape the frozen bundle")
    helpers = _features_module()
    records_hash = helpers._digest(records_path)
    manifest_path = bundle / "memory/memory_manifest.json"
    if manifest_path.is_file():
        if bundle not in manifest_path.resolve().parents:
            raise ValueError("memory manifest escapes the frozen bundle")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("artifacts", {}).get("agent_memory.jsonl", {}).get("sha256") != records_hash:
            raise ValueError("frozen memory records do not match the release manifest SHA-256")
    inputs_root = (bundle / "inputs/pocket_prediction").resolve()
    if bundle not in inputs_root.parents:
        raise ValueError("pocket inputs root escapes the frozen memory bundle")
    cases, seen, encoders, observed_tools, releases = [], set(), set(), set(), set()
    inputs_digest = hashlib.sha256()
    with records_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                if not isinstance(record, dict):
                    raise ValueError("memory record must be an object")
                if record.get("task_context", {}).get("task_type") != "pocket_prediction":
                    continue
                split = record.get("extensions", {}).get("dataset_split", {})
                if split.get("split") != "memory":
                    raise ValueError("pocket index accepts only declared memory split records")
                case_id, group_id = record.get("identity", {}).get("case_id"), split.get("group_id")
                if not isinstance(case_id, str) or not case_id or case_id in seen:
                    raise ValueError(f"invalid or duplicate pocket case identity: {case_id!r}")
                if not isinstance(group_id, str) or not group_id:
                    raise ValueError(f"missing protein group: {case_id}")
                case_dir = (inputs_root / case_id).resolve()
                if inputs_root not in case_dir.parents:
                    raise ValueError(f"pocket case escapes memory inputs: {case_id}")
                protein, encoder_id, hashes = _protein_features(record, case_dir, case_id, group_id, helpers)
                values, tools, publication_releases = _published_dca(record)
                cases.append({"case_id": case_id, "group_id": group_id, "protein": protein,
                              "top1_dca_angstrom": values})
                seen.add(case_id)
                encoders.add(encoder_id)
                observed_tools.update(tools)
                releases.update(publication_releases)
                for name, digest in sorted(hashes.items()):
                    inputs_digest.update(f"{case_id}/{name}\0{digest}\n".encode())
            except (ValueError, KeyError, TypeError, AttributeError, OSError) as exc:
                raise ValueError(f"invalid pocket memory record at line {line_number}: {exc}") from exc
    if not cases:
        raise ValueError("no declared pocket_prediction memory cases found")
    if len(encoders) != 1:
        raise ValueError("pocket memory contains incompatible protein encoders")
    policy = _policy_module()
    tools = list(policy.PAPER_TOOLS)
    for case in cases:
        case["top1_dca_angstrom"] = {tool: case["top1_dca_angstrom"].get(tool) for tool in tools}
        case["best_detector"] = policy.infer_best_detector(case["top1_dca_angstrom"])
    cases.sort(key=lambda case: case["case_id"])
    source = {"bundle": str(bundle), "records": str(records_path), "records_sha256": records_hash,
              "feature_inputs_sha256": inputs_digest.hexdigest(), "split": "memory",
              "release_manifest_sha256": helpers._digest(manifest_path) if manifest_path.is_file() else None,
              "release_manifest_verified": manifest_path.is_file(), "published_releases": sorted(releases),
              "observed_tools": sorted(observed_tools),
              "assessment_policy": "published pocket_target_summary raw native top1_dca_angstrom only",
              "best_detector_derivation": policy.BEST_DETECTOR_POLICY,
              "inference_note": "Best detector is the minimum raw DCA winner; exact ties use lexicographic tool ID."}
    index = {"schema_version": SCHEMA_VERSION, "index_type": INDEX_TYPE,
             "policy_version": POLICY_VERSION, "encoder_compatibility_id": next(iter(encoders)),
             "best_detector_policy": policy.BEST_DETECTOR_POLICY,
             "memory_version": hashlib.sha256(json.dumps(source, sort_keys=True).encode()).hexdigest(),
             "source": source, "tools": tools, "cases": cases}
    policy.validate_index(index)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                         prefix=".pocket-index-", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(index, handle, separators=(",", ":"), allow_nan=False)
            handle.write("\n")
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return {"success": True, "status": "built", "path": str(output), "case_count": len(cases),
            "group_count": len({case["group_id"] for case in cases}), "tools": tools,
            "memory_version": index["memory_version"], "policy_version": POLICY_VERSION,
            "best_detector_policy": policy.BEST_DETECTOR_POLICY,
            "best_detector_counts": {tool: sum(case["best_detector"] == tool for case in cases) for tool in tools},
            "missing_labels": {tool: sum(case["top1_dca_angstrom"][tool] is None for case in cases)
                               for tool in tools}}

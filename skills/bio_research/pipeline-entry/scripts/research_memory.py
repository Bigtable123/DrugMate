"""Append and retrieve evaluated research observations without changing snapshots.

This store supplies measured evidence to the planner. It does not manufacture
cross-tool winners, average different protocols, or update frozen selectors.
Evidence references use {path, sha256, pointer}, with RFC 6901 JSON pointers.
"""
from __future__ import annotations

from contextlib import contextmanager
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
from statistics import pstdev
import sys
import tempfile

TASKS = ("docking", "pocket_prediction", "molecule_generation")
FINISHED = {"complete", "completed", "success", "failed", "partial", "interrupted", "awaiting_final_confirmation", "publication_failed"}
SCHEMA_VERSION = "research_observation/1"
MAX_VALUE_BYTES = 512 * 1024


def _module(name, path=None):
    spec = importlib.util.spec_from_file_location(name, path or Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    previous = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


def _json_bytes(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False,
                      separators=(",", ":")).encode("utf-8")


def _identity(value):
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _event_identity(payload):
    def bound_value(item):
        return {key: item[key] for key in ("sha256", "pointer", "value")}
    def unordered(items):
        unique = {_identity(bound_value(item)): bound_value(item) for item in items}
        return [unique[key] for key in sorted(unique)]
    return _identity({"profile_id": payload["profile_id"], "run_id": payload["run_id"],
                      "execution": unordered(payload["execution"]),
                      "evidence": {name: bound_value(item) for name, item in payload["evidence"].items()},
                      "measurements": unordered(payload["measurements"])})


def _semantic_event_identity(payload):
    """A cached execution is not a new experiment when its manifest is rewritten.

    File hashes still seal every imported reference. Deduplication additionally
    uses the execution identities and selected values, so unrelated timestamps
    or JSON formatting cannot turn the same observations into additional runs.
    """
    def selected(item):
        return {"pointer": item["pointer"], "value": item["value"]}
    execution = {receipt["value"]["execution_id"] for receipt in payload["execution"]}
    measurements = {_identity(selected(item)): selected(item) for item in payload["measurements"]}
    return _identity({"profile_id": payload["profile_id"], "run_id": payload["run_id"],
                      "execution_ids": sorted(execution),
                      "evidence": {name: selected(item) for name, item in payload["evidence"].items()},
                      "measurements": [measurements[key] for key in sorted(measurements)]})


def _terminal_status(status):
    # The runner records a terminal stage failure as failed_<step_id>, not
    # simply failed. Intermediate <step_id>_failed/best-effort states do not
    # match this prefix and must wait until the round finishes.
    return isinstance(status, str) and (status in FINISHED or
                                       (status.startswith("failed_") and bool(status[7:])))


def _tools(value):
    if isinstance(value, str):
        return [_text(value, "observed tool")]
    if not isinstance(value, list) or not value:
        raise ValueError("observed tools must be a tool ID or a nonempty list of tool IDs")
    names = [_text(item, "observed tool") for item in value]
    if len(set(names)) != len(names):
        raise ValueError("observed tool list must not contain duplicates")
    return names


def _text(value, name):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be nonempty text without outer whitespace")
    return value


def _digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _disabled(mode):
    if mode not in ("research", "frozen"):
        raise ValueError("research memory mode must be research or frozen")
    return {"success": True, "status": "disabled", "mode": "frozen", "profiles": []}


def _locations(store):
    deployment = _module("deployment_config", Path(__file__).resolve().parents[2] / "shared/deployment_config.py")
    root = deployment.data_root().resolve()
    checkout = Path(__file__).resolve().parents[4]
    if root == checkout or checkout in root.parents:
        raise ValueError("SBDD_DATA_ROOT must remain outside the source checkout")
    directory = Path(store).expanduser().resolve() if store is not None else root / "research_memory"
    _runtime_path(directory, root, checkout, is_file=False)
    if directory == root:
        raise ValueError("research store must be a dedicated directory under SBDD_DATA_ROOT")
    return root, checkout, directory


def _runtime_path(value, root, checkout, *, is_file=True):
    path = Path(value).expanduser().resolve()
    if root not in path.parents or path == checkout or checkout in path.parents:
        raise ValueError("research inputs and store must be under SBDD_DATA_ROOT, outside the checkout")
    for parent in (path, *path.parents):
        if parent == root:
            break
        if (parent / "bundle_manifest.json").is_file() or (parent / "memory_manifest.json").is_file():
            raise ValueError("frozen memory bundles cannot be research inputs or stores")
    if is_file and not path.is_file():
        raise ValueError(f"research input is not a file: {path}")
    return path


def _pointer(document, pointer):
    if not isinstance(pointer, str) or (pointer and not pointer.startswith("/")):
        raise ValueError("evidence pointer must be an RFC 6901 JSON pointer")
    value = document
    for encoded in pointer.split("/")[1:] if pointer else ():
        # Reject invalid escape sequences instead of interpreting a different key.
        rest = encoded.replace("~0", "").replace("~1", "")
        if "~" in rest:
            raise ValueError("invalid JSON pointer escape")
        token = encoded.replace("~1", "/").replace("~0", "~")
        if isinstance(value, list):
            if not token.isdigit() or (len(token) > 1 and token.startswith("0")):
                raise ValueError("invalid JSON pointer array index")
            value = value[int(token)]
        elif isinstance(value, dict):
            value = value[token]
        else:
            raise ValueError("JSON pointer traverses a scalar")
    if len(_json_bytes(value)) > MAX_VALUE_BYTES:
        raise ValueError("selected evidence is too large; point to a bounded measurement block")
    return value


def _evidence(reference, root, checkout):
    if not isinstance(reference, dict) or set(reference) != {"path", "sha256", "pointer"}:
        raise ValueError("evidence requires exactly path, sha256 and pointer")
    path = _runtime_path(reference["path"], root, checkout)
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if reference["sha256"] != digest:
        raise ValueError(f"research evidence SHA-256 mismatch: {path}")
    value = _pointer(json.loads(data), reference["pointer"])
    return {"path": str(path), "sha256": digest, "pointer": reference["pointer"], "value": value}


def _library_ligand(record, inputs, receipts, root, checkout):
    """Bind a single original library record through the runner's preparation map."""
    source = record.get("ligand_source")
    if source is None:
        return None
    if (record["task"] != "docking" or not isinstance(source, dict)
            or set(source) != {"path", "sha256", "record_index", "identity"}):
        raise ValueError("docking ligand_source requires path, sha256, record_index and identity")
    index = source["record_index"]
    if type(index) is not int or index < 0:
        raise ValueError("ligand_source record_index must be a nonnegative integer")
    path = _runtime_path(source["path"], root, checkout)
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    primary = [receipt["value"] for receipt in receipts if receipt["value"]["step_id"] == "docking"]
    identity = _evidence(source["identity"], root, checkout)
    # The input library and its preparation map must belong to one execution,
    # not merely to two unrelated stages in the same project directory.
    if source["sha256"] != digest or not any(
        receipt["input_file_signatures"].get(str(path), {}).get("sha256") == digest
        and receipt["evidence_file_signatures"].get(identity["path"], {}).get("sha256") == identity["sha256"]
        for receipt in primary
    ):
        raise ValueError("ligand library and preparation identity are not bound to the same docking execution")
    records, current = [], []
    for line in data.splitlines(keepends=True):
        current.append(line)
        if line.strip() == b"$$$$":
            records.append(b"".join(current))
            current = []
    if b"".join(current).strip() or index >= len(records):
        raise ValueError("ligand source requires complete SDF records and an existing record_index")
    if records[index] != Path(inputs["ligand_sdf"]["path"]).read_bytes():
        raise ValueError("query ligand is not the byte-preserving source library record")
    document = _read(identity["path"])
    prepared = document.get("records")
    if (Path(identity["path"]).name != "microstate_identity_map.json"
            or identity["pointer"] != f"/records/{index}"
            or document.get("policy") != "preserve_declared"
            or document.get("source", {}).get("path") != str(path)
            or document.get("source", {}).get("sha256") != digest
            or not isinstance(prepared, list) or len(prepared) != len(records)
            or any(not isinstance(row, dict) or type(row.get("source_molecule_index")) is not int
                   or row["source_molecule_index"] != position for position, row in enumerate(prepared))):
        raise ValueError("preparation map does not preserve this source library's record identities")
    ids = [_text(row.get("compound_id"), "prepared compound ID") for row in prepared]
    if len(set(ids)) != len(ids):
        raise ValueError("preparation map compound IDs are not unique")
    return {"path": str(path), "sha256": digest, "record_index": index,
            "record_sha256": inputs["ligand_sdf"]["sha256"], "identity": identity}


def _molecule_measurements(measurements, ligand_source):
    if ligand_source is None:
        return
    compound = ligand_source["identity"]["value"]["compound_id"]
    index = ligand_source["record_index"]
    for item in measurements:
        rows = item["value"] if isinstance(item["value"], list) else [item["value"]]
        for row in rows:
            if not isinstance(row, dict) or row.get("compound_id") != compound:
                raise ValueError("per-ligand measurements must retain the matching prepared compound_id")
            if "source_molecule_index" in row and (type(row["source_molecule_index"]) is not int
                                                   or row["source_molecule_index"] != index):
                raise ValueError("measurement source_molecule_index differs from the prepared ligand identity")


def _execution(record, inputs, root, checkout):
    references = record.get("execution")
    if not isinstance(references, list):
        references = [references]
    if not references:
        raise ValueError("research import requires actual pipeline execution provenance")
    receipts = [_evidence(reference, root, checkout) for reference in references]
    run = _runtime_path(record["run_id"], root, checkout, is_file=False)
    if record["run_id"] != str(run):
        raise ValueError("run_id must be the canonical absolute round project directory")
    input_signatures, output_signatures, stages = {}, {}, set()
    for receipt in receipts:
        value = receipt["value"]
        if (receipt["pointer"] != "/research_provenance" or Path(receipt["path"]).name != "result.json"
                or run not in Path(receipt["path"]).parents or not isinstance(value, dict)
                or type(value.get("schema_version")) is not int or value["schema_version"] != 1
                or value.get("project_dir") != str(run)):
            raise ValueError("execution receipt is not bound to this run's original result.json")
        _text(value.get("execution_id"), "execution identity")
        stages.add(value.get("step_id"))
        if not isinstance(value.get("input_file_signatures"), dict) or not isinstance(value.get("evidence_file_signatures"), dict):
            raise ValueError("execution receipt has no original input/output signatures")
        input_signatures.update(value["input_file_signatures"])
        output_signatures.update(value["evidence_file_signatures"])
        output_signatures[receipt["path"]] = {"sha256": receipt["sha256"]}
    primary = {"docking": "docking", "pocket_prediction": "pocket_analysis",
               "molecule_generation": "molecule_generation"}[record["task"]]
    if primary not in stages:
        raise ValueError(f"research import requires the actual {primary} execution receipt")
    ligand_source = _library_ligand(record, inputs, receipts, root, checkout)
    for name in ("receptor_file", "ligand_sdf"):
        if name == "ligand_sdf" and ligand_source is not None:
            continue
        if name in inputs and not any(signature.get("sha256") == inputs[name]["sha256"]
                                      for signature in input_signatures.values()):
            raise ValueError(f"{name} is not a captured input of the evaluated run; do not assume receptor or ligand equivalence")
    return receipts, run, output_signatures, ligand_source


def _check_run_evidence(evidence, measurements, run, signatures):
    status = evidence["status"]
    if Path(status["path"]) != run / "manifest.json" or status["pointer"] != "/status":
        raise ValueError("run status must come from this round's manifest.json /status")
    if _read(status["path"]).get("project_dir") != str(run):
        raise ValueError("status manifest identifies a different run")
    for item in [evidence["tool"], evidence["parameters"], *measurements]:
        path = Path(item["path"])
        if run not in path.parents or signatures.get(str(path), {}).get("sha256") != item["sha256"]:
            raise ValueError("research evidence was not sealed as an output of the supplied execution receipts")


def _features(task, inputs, root, checkout):
    if task not in TASKS or not isinstance(inputs, dict):
        raise ValueError("research task/inputs are invalid")
    names = {"receptor_file", "query_npz", "query_metadata"}
    if task != "pocket_prediction":
        names.add("ligand_sdf")
    if set(inputs) != names:
        raise ValueError("research inputs require " + ", ".join(sorted(names)))
    paths = {name: _runtime_path(inputs[name], root, checkout) for name in names}
    selector = _module("docking_selector")
    metadata = _read(paths["query_metadata"])
    source = metadata.get("source_record", {})
    if source.get("task_type") != task:
        raise ValueError("research task differs from original encoder metadata")
    encoder = selector.encoder_id_from_metadata(metadata)
    hashes = {name: _digest(path) for name, path in paths.items()}
    if metadata.get("embedding_sha256") != hashes["query_npz"]:
        raise ValueError("research embedding SHA-256 does not match metadata")
    if metadata.get("source_input_sha256") != hashes["receptor_file"]:
        raise ValueError("research receptor SHA-256 does not match metadata")
    if source.get("protein_sha256") not in (None, hashes["receptor_file"]):
        raise ValueError("research receptor differs from encoder source_record")
    shape, vector = selector.read_npz_float32(paths["query_npz"], "protein_embedding_normalized")
    protein = selector._finite_vector(vector, 1280, "research protein embedding")
    if shape != (1280,) or abs(math.sqrt(sum(x * x for x in protein)) - 1.0) > 1e-5:
        raise ValueError("research protein vector must have shape (1280,) and unit L2 norm")
    features = {"feature_version": "protein1280_v1", "protein": protein}
    if task != "pocket_prediction":
        extractor = _module("generation_selector") if task == "molecule_generation" else selector
        features = extractor.extract_query_features(paths["ligand_sdf"], paths["query_npz"])
    case = source.get("case_id")
    derived = f"{source['dataset']}/{source['target']}" if source.get("dataset") and source.get("target") else None
    if case and derived and case != derived:
        raise ValueError("research encoder case differs from dataset/target")
    identity = {"task": task, "receptor_sha256": hashes["receptor_file"],
                "ligand_input_sha256": hashes.get("ligand_sdf"),
                "feature_version": features["feature_version"], "encoder_compatibility_id": encoder}
    return features, identity, {name: {"path": str(paths[name]), "sha256": hashes[name]} for name in names}, {
        "case_id": case or derived, "group_id": source.get("group_id")}


@contextmanager
def _lock(directory):
    directory.mkdir(parents=True, exist_ok=True)
    lock_path = directory / ".append.lock"
    if lock_path.is_symlink():
        raise ValueError("research store lock must not be a symbolic link")
    with lock_path.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _observations(directory):
    entries = directory / "observations"
    if entries.is_symlink() or (entries.exists() and not entries.is_dir()):
        raise ValueError("research observations must be a local store directory")
    for path in sorted(entries.glob("*.json")):
        if path.is_symlink():
            raise ValueError("research observation must not be a symbolic link")
        payload = _read(path)
        if not isinstance(payload, dict):
            raise ValueError(f"research observation must be an object: {path}")
        integrity = payload.pop("content_sha256", None)
        if (payload.get("schema_version") != SCHEMA_VERSION
                or _identity(payload) != integrity or _event_identity(payload) != path.stem
                or _identity(payload.get("profile_identity")) != payload.get("profile_id")):
            raise ValueError(f"research observation integrity check failed: {path}")
        yield path.stem, payload


def append_run(record_path, store=None, *, mode="research"):
    disabled = _disabled(mode)
    if mode == "frozen":
        return disabled
    root, checkout, directory = _locations(store)
    record_path = _runtime_path(record_path, root, checkout)
    record = _read(record_path)
    if not isinstance(record, dict) or type(record.get("schema_version")) is not int or record["schema_version"] != 1:
        raise ValueError("research import record requires schema_version=1")
    run_id, case_id, group_id = (_text(record.get(name), name) for name in ("run_id", "case_id", "group_id"))
    features, identity, inputs, binding = _features(record.get("task"), record.get("inputs"), root, checkout)
    execution, run, signatures, ligand_source = _execution(record, inputs, root, checkout)
    if binding["case_id"] not in (None, case_id) or binding["group_id"] not in (None, group_id):
        raise ValueError("research case/group differs from bound encoder metadata")
    references = record.get("evidence")
    if not isinstance(references, dict) or set(references) != {"status", "tool", "parameters"}:
        raise ValueError("research evidence requires status, tool and parameters references")
    evidence = {name: _evidence(reference, root, checkout) for name, reference in references.items()}
    status = _text(evidence["status"]["value"], "observed run status")
    if not _terminal_status(status):
        raise ValueError("only terminal research run statuses may be appended: " + ", ".join(sorted(FINISHED)) + ", failed_<step_id>")
    _tools(evidence["tool"]["value"])
    if not isinstance(evidence["parameters"]["value"], dict):
        raise ValueError("observed protocol parameters must be an object")
    if not isinstance(record.get("measurements"), list) or not record["measurements"]:
        raise ValueError("supply at least one actual measurement or explicit unavailable/failed assessment")
    measurements = [_evidence(reference, root, checkout) for reference in record["measurements"]]
    if any(not isinstance(item["value"], (dict, list)) or not item["value"] for item in measurements):
        raise ValueError("measurement references must retain named fields and any status/unit information")
    _check_run_evidence(evidence, measurements, run, signatures)
    _molecule_measurements(measurements, ligand_source)
    payload = {"schema_version": SCHEMA_VERSION, "profile_id": _identity(identity),
               "profile_identity": identity, "features": features, "run_id": run_id,
               "case_id": case_id, "group_id": group_id, "run_status": status,
               "inputs": inputs, "evidence": evidence, "measurements": measurements,
               "execution": execution,
               "import_record": {"path": str(record_path), "sha256": _digest(record_path)}}
    if ligand_source is not None:
        payload["ligand_source"] = ligand_source
    observation_id = _event_identity(payload)
    with _lock(directory):
        previous = list(_observations(directory))
        # Different aliases may identify the same exact prepared inputs, but a
        # protein group cannot change silently for an existing profile.
        if any(item["profile_id"] == payload["profile_id"] and item["group_id"] != group_id for _, item in previous):
            raise ValueError("research profile already has a different protein group")
        if any(item["profile_id"] == payload["profile_id"] and item["features"] != features for _, item in previous):
            raise ValueError("research profile has different features for the same bound inputs and encoder")
        semantic_identity = _semantic_event_identity(payload)
        for previous_id, item in previous:
            if _semantic_event_identity(item) == semantic_identity:
                return {"success": True, "status": "already_recorded", "profile_id": payload["profile_id"],
                        "observation_id": previous_id, "store": str(directory)}
        destination = directory / "observations" / f"{observation_id}.json"
        if destination.exists():
            return {"success": True, "status": "already_recorded", "profile_id": payload["profile_id"],
                    "observation_id": observation_id, "store": str(directory)}
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".observation-", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(_json_bytes({**payload, "content_sha256": _identity(payload)}) + b"\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return {"success": True, "status": "appended", "profile_id": payload["profile_id"],
            "observation_id": observation_id, "store": str(directory), "run_status": status,
            "profile_created": not any(item["profile_id"] == payload["profile_id"] for _, item in previous)}


def recall(task, inputs, store=None, *, mode="research", k=None,
           exclude_cases=(), exclude_groups=(), exclude_runs=()):
    disabled = _disabled(mode)
    if mode == "frozen":
        return disabled
    expected_k = 5 if task == "pocket_prediction" else 3
    k = expected_k if k is None else k
    if type(k) is not int or not 1 <= k <= 20:
        raise ValueError("research recall requires 1 <= k <= 20")
    root, checkout, directory = _locations(store)
    query, identity, query_inputs, _ = _features(task, inputs, root, checkout)
    cases, groups, runs = set(exclude_cases), set(exclude_groups), set(exclude_runs)
    profiles = {}
    for observation_id, item in _observations(directory):
        context = item["profile_identity"]
        if (context["task"] != task or context["feature_version"] != identity["feature_version"]
                or context["encoder_compatibility_id"] != identity["encoder_compatibility_id"]
                or item["case_id"] in cases or item["group_id"] in groups
                or item["run_id"] in runs or not _terminal_status(item["run_status"])):
            continue
        profile = profiles.setdefault(item["profile_id"], {"profile_id": item["profile_id"],
            "profile_identity": context, "features": item["features"], "group_id": item["group_id"],
            "case_ids": set(), "observations": []})
        profile["case_ids"].add(item["case_id"])
        profile["observations"].append({"observation_id": observation_id,
            "execution": [{"path": receipt["path"], "sha256": receipt["sha256"],
                           "execution_id": receipt["value"]["execution_id"],
                           "step_id": receipt["value"]["step_id"]} for receipt in item["execution"]],
            **{name: item[name] for name in ("run_id", "run_status", "case_id", "evidence", "measurements", "import_record")},
            **({"ligand_source": item["ligand_source"]} if "ligand_source" in item else {})})
    entries = sorted(profiles.values(), key=lambda item: item["profile_id"])
    scales = [pstdev(item["features"]["ligand_desc"][j] for item in entries) or 1.0 for j in range(18)] if entries and task != "pocket_prediction" else []
    for item in entries:
        features = item["features"]
        protein = sum(a * b for a, b in zip(query["protein"], features["protein"]))
        similarity = {"protein_cosine": protein}
        if scales:
            a, b = set(query["ligand_fp"]), set(features["ligand_fp"])
            tanimoto = len(a & b) / len(a | b) if a | b else 0.0
            descriptor = math.exp(-0.5 * sum(((x-y)/s)**2 for x, y, s in zip(query["ligand_desc"], features["ligand_desc"], scales)) / 18)
            similarity.update(ligand_tanimoto=tanimoto, ligand_descriptor_rbf=descriptor,
                              ligand_similarity=0.8*tanimoto + 0.2*descriptor)
        item["similarity"] = similarity
    protein_order = sorted(entries, key=lambda item: (item["similarity"]["protein_cosine"], item["profile_id"]))
    percentiles = {item["profile_id"]: rank / max(len(entries)-1, 1) for rank, item in enumerate(protein_order)}
    entries.sort(key=lambda item: (-item["similarity"].get("ligand_similarity", item["similarity"]["protein_cosine"]), item["profile_id"]))
    if task == "pocket_prediction":
        unique_groups = {}
        for item in entries:
            unique_groups.setdefault(item["group_id"], item)
        selected = list(unique_groups.values())[:k]
    else:
        selected = entries[:k]
    pocket_weights = [math.exp((item["similarity"]["protein_cosine"] - selected[0]["similarity"]["protein_cosine"]) / 0.02) / math.log2(rank+2)
                      for rank, item in enumerate(selected)] if task == "pocket_prediction" else []
    for rank, item in enumerate(selected):
        item.pop("features")
        item["case_ids"] = sorted(item["case_ids"])
        item["similarity"]["protein_percentile"] = percentiles[item["profile_id"]]
        if task != "pocket_prediction":
            item["similarity"]["evidence_weight"] = math.exp(0.10*(percentiles[item["profile_id"]]-0.5)) / math.log2(rank+2)
        else:
            item["similarity"]["evidence_weight"] = len(selected) * pocket_weights[rank] / sum(pocket_weights)
    return {"success": True, "status": "recalled", "mode": "research", "store": str(directory),
            "query_inputs": query_inputs, "eligible_profiles": len(entries), "profiles": selected,
            "requested_k": k, "diagnostic_neighborhood": k != expected_k or len(selected) < expected_k,
            "policy": "ligand similarity then protein reweighting" if scales else "protein cosine",
            "note": "Research evidence for planner review; no cross-run averaging, inferred winners, or frozen-selector updates."}

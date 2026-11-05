"""Frozen ligand-first docking memory policy; no query labels or evaluation inputs.

The graph fingerprint uses a custom circular encoding. Changing either feature
algorithm requires a new feature version.
"""
from __future__ import annotations

import array
import ast
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
from statistics import median, pstdev
import struct
import sys
import tempfile
from typing import Any, Sequence
import zipfile

FORMAT_VERSION = "sbdd_docking_selector_index_v2"
LEGACY_FORMAT_VERSION = "sbdd_docking_selector_index_v1"
FEATURE_VERSION = "custom_circular_r0_3_2048_desc18_protein1280_v1"
POLICY_VERSION = "ligand_protein_reweight_v1"
PUBLISHED_RELEASE = "known_pocket_v2_final_20261009"
TOOLS = ("autodock4", "autodock_gpu", "autodock_vina", "diffdock", "dock6", "gnina", "rdock", "unidock")
# AutoDock-GPU has no current Step 4 adapter; do not relabel it as CPU AutoDock4.
AVAILABLE_TOOLS = tuple(tool for tool in TOOLS if tool != "autodock_gpu")
PARAMETERS = {"k": 3, "protein_beta": 0.10, "prior_strength": 2.0,
              "switch_margin": 0.10, "default_tool": "gnina"}
DEFAULT_METRIC = "top1_rmsd_lt_2a"
TASK_DEFAULTS = {
    "redocking": {"default_tool": "gnina", "default_metric": DEFAULT_METRIC},
    "pose_generation": {"default_tool": "rdock", "default_metric": "posebusters_valid"},
}
METRICS = {DEFAULT_METRIC: "binary", "posebusters_valid": "binary", "joint_success": "binary",
           "top1_rmsd_a": "rmsd", "median_rmsd_a": "rmsd"}
RMSD_LOSS_CAP = 10.0


def _read_npy_stream(stream: Any) -> tuple[tuple[int, ...], array.array]:
    """Read a C-contiguous little-endian float32 NPY array from a stream."""
    magic = stream.read(6)
    if magic != b"\x93NUMPY":
        raise ValueError("not an NPY stream")
    major, minor = struct.unpack("BB", stream.read(2))
    if major == 1:
        header_len = struct.unpack("<H", stream.read(2))[0]
    elif major in (2, 3):
        header_len = struct.unpack("<I", stream.read(4))[0]
    else:
        raise ValueError(f"unsupported NPY version {major}.{minor}")
    header = ast.literal_eval(stream.read(header_len).decode("latin1"))
    if not isinstance(header, dict):
        raise ValueError("invalid NPY header")
    if header.get("fortran_order"):
        raise ValueError("Fortran-order NPY arrays are not supported")
    if header.get("descr") not in ("<f4", "=f4", "|f4"):
        raise ValueError(f"expected float32 NPY array, got {header.get('descr')!r}")
    shape = tuple(int(x) for x in header["shape"])
    if any(dimension < 0 for dimension in shape):
        raise ValueError("invalid NPY dimensions")
    count = math.prod(shape)
    values = array.array("f")
    values.frombytes(stream.read(count * 4))
    if len(values) != count:
        raise ValueError(f"truncated NPY data: expected {count}, got {len(values)}")
    if sys.byteorder == "big" and header.get("descr") == "<f4":
        values.byteswap()
    return shape, values


def read_npz_float32(path: Path, array_name: str) -> tuple[tuple[int, ...], array.array]:
    member = f"{array_name}.npy"
    try:
        with zipfile.ZipFile(path) as archive:
            with archive.open(member) as handle:
                return _read_npy_stream(handle)
    except (struct.error, SyntaxError, EOFError, KeyError, TypeError, OverflowError,
            zipfile.BadZipFile) as exc:
        raise ValueError(f"invalid protein embedding NPZ: {exc}") from exc


def stable_bit(text: str, nbits: int = 2048) -> int:
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % nbits


def _sdf_graph(path: Path) -> tuple[frozenset[int], tuple[float, ...]]:
    """Read the first V2000 molecule and return a circular fingerprint/descriptors."""
    return _sdf_graph_lines(path.read_text(encoding="utf-8", errors="replace").splitlines(), str(path))


def _sdf_graph_lines(lines: list[str], path: str = "<ligand>") -> tuple[frozenset[int], tuple[float, ...]]:
    """Share the frozen graph parser with in-memory generation compatibility views."""
    if len(lines) < 4 or "V2000" not in lines[3]:
        raise ValueError(f"expected V2000 SDF: {path}")
    atom_count = int(lines[3][0:3])
    bond_count = int(lines[3][3:6])
    if atom_count <= 0 or bond_count < 0 or len(lines) < 4 + atom_count + bond_count:
        raise ValueError(f"invalid atom/bond counts: {path}")
    atoms: list[tuple[str, float, float, float, int]] = []
    for line in lines[4 : 4 + atom_count]:
        atoms.append(
            (
                line[31:34].strip(),
                float(line[0:10]),
                float(line[10:20]),
                float(line[20:30]),
                int(line[36:39] or 0),
            )
        )
    if any(not atom[0] or any(not math.isfinite(value) for value in atom[1:4]) for atom in atoms):
        raise ValueError(f"invalid atoms or coordinates: {path}")
    bonds: list[tuple[int, int, int]] = []
    for line in lines[4 + atom_count : 4 + atom_count + bond_count]:
        bonds.append((int(line[0:3]) - 1, int(line[3:6]) - 1, int(line[6:9])))
    if any(a < 0 or b < 0 or a >= atom_count or b >= atom_count or order not in (1, 2, 3, 4)
           for a, b, order in bonds):
        raise ValueError(f"invalid bond indices or unsupported bond order: {path}")

    charges = [0] * atom_count
    # V2000 atom charge codes; M  CHG below overrides these values.
    charge_code = {1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3}
    for index, atom in enumerate(atoms):
        charges[index] = charge_code.get(atom[4], 0)
    for line in lines[4 + atom_count + bond_count :]:
        if line.startswith("M  CHG"):
            fields = line.split()
            for offset in range(int(fields[2])):
                charges[int(fields[3 + 2 * offset]) - 1] = int(fields[4 + 2 * offset])
        if line.startswith("M  END"):
            break

    heavy_old = [index for index, atom in enumerate(atoms) if atom[0].upper() not in ("H", "D")]
    old_to_new = {old: new for new, old in enumerate(heavy_old)}
    heavy_atoms = [atoms[old] for old in heavy_old]
    heavy_charges = [charges[old] for old in heavy_old]
    adjacency: list[list[tuple[int, int]]] = [[] for _ in heavy_atoms]
    heavy_bonds: list[tuple[int, int, int]] = []
    for old_a, old_b, order in bonds:
        if old_a not in old_to_new or old_b not in old_to_new:
            continue
        a, b = old_to_new[old_a], old_to_new[old_b]
        adjacency[a].append((b, order))
        adjacency[b].append((a, order))
        heavy_bonds.append((a, b, order))

    labels: list[str] = []
    bits: set[int] = set()
    for index, atom in enumerate(heavy_atoms):
        element = atom[0]
        degree = len(adjacency[index])
        valence = sum(min(order, 3) for _, order in adjacency[index])
        aromatic = int(any(order == 4 for _, order in adjacency[index]))
        label = f"{element}|{degree}|{valence}|{heavy_charges[index]}|{aromatic}"
        labels.append(label)
        bits.add(stable_bit("r0|" + label))
    for radius in range(1, 4):
        updated = []
        for index, label in enumerate(labels):
            neighborhood = ";".join(
                sorted(f"{order}:{labels[neighbor]}" for neighbor, order in adjacency[index])
            )
            new_label = hashlib.blake2b(
                f"{label}[{neighborhood}]".encode("utf-8"), digest_size=12
            ).hexdigest()
            updated.append(new_label)
            bits.add(stable_bit(f"r{radius}|{new_label}"))
        labels = updated

    elements = Counter(atom[0] for atom in heavy_atoms)
    n = max(len(heavy_atoms), 1)
    edge_count = len(heavy_bonds)
    # Cyclomatic number for the heavy-atom graph.
    components = 0
    seen: set[int] = set()
    for start in range(len(heavy_atoms)):
        if start in seen:
            continue
        components += 1
        stack = [start]
        seen.add(start)
        while stack:
            current = stack.pop()
            for neighbor, _ in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
    ring_rank = max(0, edge_count - len(heavy_atoms) + components)

    coords = [(a[1], a[2], a[3]) for a in heavy_atoms]
    if coords:
        center = tuple(sum(point[d] for point in coords) / len(coords) for d in range(3))
        radii = [math.dist(point, center) for point in coords]
        radius_gyration = math.sqrt(sum(r * r for r in radii) / len(radii))
        max_span = max((math.dist(a, b) for i, a in enumerate(coords) for b in coords[i + 1 :]), default=0.0)
    else:
        radius_gyration = 0.0
        max_span = 0.0

    descriptor = (
        float(len(heavy_atoms)),
        float(edge_count),
        float(ring_rank),
        sum(heavy_charges),
        sum(abs(value) for value in heavy_charges),
        sum(1 for _, _, order in heavy_bonds if order == 2),
        sum(1 for _, _, order in heavy_bonds if order == 3),
        sum(1 for _, _, order in heavy_bonds if order == 4),
        sum(1 for neighbors in adjacency if len(neighbors) == 1) / n,
        sum(1 for neighbors in adjacency if len(neighbors) >= 3) / n,
        elements["C"] / n,
        elements["N"] / n,
        elements["O"] / n,
        elements["S"] / n,
        elements["P"] / n,
        sum(elements[element] for element in ("F", "Cl", "Br", "I")) / n,
        radius_gyration,
        max_span,
    )
    return frozenset(bits), tuple(float(x) for x in descriptor)


def _finite_vector(values, length: int, name: str) -> list[float]:
    if not isinstance(values, (list, tuple, array.array)) or len(values) != length:
        raise ValueError(f"{name} must have {length} values")
    if any(isinstance(value, bool) or not isinstance(value, (int, float))
           or not math.isfinite(value) for value in values):
        raise ValueError(f"{name} contains a non-finite or non-numeric value")
    return list(values)


def _validate_features(features: dict, require_protein: bool = True) -> None:
    if not isinstance(features, dict):
        raise ValueError("query features must be an object")
    if features.get("feature_version") != FEATURE_VERSION:
        raise ValueError("incompatible feature version")
    fp = features.get("ligand_fp")
    if not isinstance(fp, (list, tuple)) or not fp or any(
        type(bit) is not int or bit < 0 or bit >= 2048 for bit in fp
    ) or len(set(fp)) != len(fp):
        raise ValueError("invalid ligand fingerprint")
    _finite_vector(features.get("ligand_desc"), 18, "ligand descriptors")
    if not require_protein:
        return
    protein = _finite_vector(features.get("protein"), 1280, "protein embedding")
    norm = math.sqrt(sum(value * value for value in protein))
    if abs(norm - 1.0) > 1e-5:
        raise ValueError("protein embedding is not L2 normalized")


def extract_query_features(ligand_sdf: Path, embeddings_npz: Path | None = None) -> dict:
    """Extract only query-time inputs; never read a reference pose or labels."""
    return extract_query_features_text(Path(ligand_sdf).read_text(encoding="utf-8", errors="replace"),
                                       embeddings_npz, str(ligand_sdf))


def extract_query_features_text(ligand_text: str, embeddings_npz: Path | None = None,
                                source_label: str = "<ligand>") -> dict:
    """Extract a single molecular feature view without writing temporary files."""
    molecules = ligand_text.split("$$$$")
    if len(molecules) > 2 or (len(molecules) == 2 and molecules[1].strip()):
        raise ValueError("docking memory selection requires a single query ligand, not a ligand library")
    fp, desc = _sdf_graph_lines(ligand_text.splitlines(), source_label)
    features = {"feature_version": FEATURE_VERSION, "ligand_fp": sorted(fp),
                "ligand_desc": list(desc), "protein": None}
    if embeddings_npz is None:
        _validate_features(features, require_protein=False)
        return features
    shape, protein = read_npz_float32(Path(embeddings_npz), "protein_embedding_normalized")
    if shape != (1280,):
        raise ValueError(f"unexpected protein embedding shape {shape}")
    norm = math.sqrt(sum(value * value for value in protein))
    if not math.isfinite(norm) or norm == 0:
        raise ValueError("invalid protein embedding norm")
    # Preserve the float32 normalization and tolerance of the frozen experiment.
    if abs(norm - 1.0) > 1e-5:
        protein = array.array("f", (value / norm for value in protein))
    features["protein"] = list(protein)
    _validate_features(features)
    return features


def _published_metrics(record: dict) -> dict:
    """Read the frozen observation's measured outcomes, retaining missing values.

    RMSDs remain raw in the index; clipping is only a selection loss transform.
    The median is over valid candidate RMSDs within one case/tool observation.
    """
    selections = {}
    for item in record.get("published_selections", []):
        if item.get("release_id") != PUBLISHED_RELEASE or item.get("tool_id") not in TOOLS:
            continue
        observation_id, tool = item["observation_id"], item["tool_id"]
        if (observation_id in selections and selections[observation_id] != tool) or any(
            previous_tool == tool and previous_id != observation_id
            for previous_id, previous_tool in selections.items()
        ):
            raise ValueError(f"conflicting published observations for {tool}")
        selections[observation_id] = tool
    candidates = {observation_id: {} for observation_id in selections}
    for observation in record.get("tool_observations", []):
        observation_id = observation.get("observation_id")
        if observation_id in selections:
            for candidate in observation.get("outputs", {}).get("candidates", []):
                rank, candidate_id = candidate.get("comparison_rank"), candidate.get("candidate_id")
                if type(rank) is int and rank > 0 and isinstance(candidate_id, str) and candidate_id:
                    previous = candidates[observation_id].get(candidate_id)
                    if previous is not None and previous != rank:
                        raise ValueError(f"conflicting published candidate ranks for {candidate_id}")
                    candidates[observation_id][candidate_id] = rank
    rmsds = {observation_id: {} for observation_id in selections}
    labels = {metric: {tool: None for tool in TOOLS} for metric in METRICS}
    for assessment in record.get("assessments", []):
        observation_id = assessment.get("observation_id")
        if observation_id not in selections:
            continue
        tool = selections[observation_id]
        if assessment.get("kind") == "pose_rmsd":
            metric = assessment.get("metrics", {}).get("rmsd", {})
            value = metric.get("value")
            if metric.get("status") != "valid" or isinstance(value, bool) or not isinstance(
                value, (int, float)
            ) or not math.isfinite(value) or value < 0:
                continue
            for candidate_id in assessment.get("candidate_ids", []):
                if candidate_id not in candidates[observation_id]:
                    continue
                previous = rmsds[observation_id].get(candidate_id)
                if previous is not None and not math.isclose(previous, value, rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f"conflicting published RMSDs for {candidate_id}")
                rmsds[observation_id][candidate_id] = value
        elif assessment.get("kind") == "top1_posebusters":
            for name in ("posebusters_valid", "joint_success"):
                metric = assessment.get("metrics", {}).get(name, {})
                if metric.get("status") != "valid":
                    continue
                value = metric.get("value")
                if type(value) is not bool:
                    raise ValueError(f"invalid published {name} for {tool}: expected boolean")
                previous = labels[name][tool]
                if previous is not None and previous != value:
                    raise ValueError(f"conflicting published {name} for {tool}")
                labels[name][tool] = value
    for observation_id, tool in selections.items():
        ranked = candidates[observation_id]
        rank1 = [candidate_id for candidate_id, rank in ranked.items() if rank == 1]
        if len(rank1) > 1:
            raise ValueError(f"multiple published rank-one poses for {tool}")
        top1 = rmsds[observation_id].get(rank1[0]) if rank1 else None
        if top1 is not None:
            labels["top1_rmsd_a"][tool] = top1
            labels[DEFAULT_METRIC][tool] = top1 < 2.0
        if rmsds[observation_id]:
            labels["median_rmsd_a"][tool] = median(rmsds[observation_id].values())
        pb, joint = labels["posebusters_valid"][tool], labels["joint_success"][tool]
        if top1 is not None and pb is not None and joint is not None and joint != (top1 < 2.0 and pb):
            raise ValueError(f"published joint success disagrees with Top-1 RMSD and PoseBusters for {tool}")
    return labels


def _published_labels(record: dict) -> dict:
    """Compatibility view of the original strict Top-1 success labels."""
    return _published_metrics(record)[DEFAULT_METRIC]


def _case_labels(case: dict, metric: str) -> dict:
    return case["labels"] if metric == DEFAULT_METRIC else case["metric_labels"][metric]


def _loss(value, metric: str):
    return min(value, RMSD_LOSS_CAP) if METRICS[metric] == "rmsd" else value


def _statistics(cases: Sequence[dict], metric: str = DEFAULT_METRIC) -> tuple[list[float], dict]:
    if not cases:
        return [1.0] * 18, {tool: None for tool in TOOLS}
    scales = [pstdev(case["ligand_desc"][j] for case in cases) for j in range(18)]
    # Eq. 1 replaces only zero-variance dimensions; small nonzero scales still
    # carry information and must not be silently discarded.
    scales = [scale if scale > 0.0 else 1.0 for scale in scales]
    priors = {}
    for tool in TOOLS:
        observed = [_loss(_case_labels(case, metric)[tool], metric) for case in cases
                    if _case_labels(case, metric)[tool] is not None]
        priors[tool] = sum(observed) / len(observed) if observed else None
    return scales, priors


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _scalar(value):
    if not isinstance(value, str):
        return value
    if value in ("", "NA", "NaN", "nan", "None", "null"):
        return None
    if value.lower() in ("true", "false"):
        return value.lower() == "true"
    try:
        number = float(value)
        if not math.isfinite(number):
            return None
        return int(number) if number.is_integer() and "." not in value and "e" not in value.lower() else number
    except ValueError:
        return value


def encoder_id_from_metadata(metadata: dict) -> str:
    """Preserve the existing memory reader's exact encoder compatibility key."""
    if not isinstance(metadata, dict) or not isinstance(metadata.get("source_record"), dict) or not isinstance(
        metadata.get("encoder"), dict
    ) or not metadata["encoder"]:
        raise ValueError("invalid encoder metadata")
    source = metadata["source_record"]
    value = {"encoder": {k: v for k, v in metadata["encoder"].items() if k not in ("model_path", "device_type")},
             "plddt_mask": _scalar(source["plddt_mask"]), "plddt_threshold": _scalar(source["plddt_threshold"])}
    canonical = json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return "encoder_compatibility:" + hashlib.sha256(canonical.encode()).hexdigest()


def _query_identity(metadata: dict) -> dict:
    if not isinstance(metadata, dict) or not isinstance(metadata.get("source_record"), dict):
        raise ValueError("query metadata must contain a source_record object")
    source = metadata["source_record"]
    if source.get("task_type") not in (None, "docking"):
        raise ValueError("query embedding metadata is not for docking")
    case_id = source.get("case_id")
    if case_id is not None and (not isinstance(case_id, str) or not case_id):
        raise ValueError("invalid query case identity")
    if source.get("group_id") is not None and (not isinstance(source["group_id"], str) or not source["group_id"]):
        raise ValueError("invalid query protein group")
    if not case_id and source.get("dataset") and source.get("target"):
        case_id = f"{source['dataset']}/{source['target']}"
    if case_id and "/" in case_id:
        dataset, target = case_id.split("/", 1)
        aliases = {"astex": "astex_diverse", "posebusters": "posebusters_benchmark"}
        case_id = f"{aliases.get(dataset, dataset)}/{target}"
    return {"case_id": case_id, "group_id": source.get("group_id")}


def query_identity(metadata_path: Path) -> dict:
    """Read identity alone for ligand-only exclusion, without touching embeddings."""
    return _query_identity(json.loads(Path(metadata_path).read_text(encoding="utf-8")))


def validate_query_metadata(metadata_path: Path, embeddings_npz: Path, index: dict,
                            receptor_file: Path | None = None) -> dict:
    metadata = json.loads(Path(metadata_path).read_text(encoding="utf-8"))
    identity = _query_identity(metadata)
    compatibility_id = encoder_id_from_metadata(metadata)
    if compatibility_id != index.get("encoder_compatibility_id"):
        raise ValueError("query encoder is incompatible with docking memory")
    embedding_hash = _digest(Path(embeddings_npz))
    if metadata.get("embedding_sha256") != embedding_hash:
        raise ValueError("query embedding SHA-256 does not match its metadata")
    source_hash = metadata.get("source_input_sha256")
    if receptor_file is not None and source_hash != _digest(Path(receptor_file)):
        raise ValueError("query receptor SHA-256 does not match embedding metadata")
    return {**identity,
            "encoder_compatibility_id": compatibility_id, "embedding_sha256": embedding_hash,
            "source_input_sha256": source_hash, "receptor_sha256_checked": receptor_file is not None}


def build_index(memory_bundle: Path, output: Path) -> dict:
    """Build outside the frozen bundle, using its memory partition exclusively."""
    bundle = Path(memory_bundle).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    if output == bundle or bundle in output.parents:
        raise ValueError("write the docking index outside the frozen memory bundle")
    records_path = bundle / "memory/agent_memory.jsonl"
    source_hash = _digest(records_path)
    manifest_path = bundle / "memory/memory_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else None
    if manifest is not None and manifest.get("artifacts", {}).get("agent_memory.jsonl", {}).get("sha256") != source_hash:
        raise ValueError("frozen memory records do not match the release manifest SHA-256")
    inputs_root = (bundle / "inputs/docking").resolve()
    cases, seen, inputs_digest, encoder_ids = [], set(), hashlib.sha256(), set()
    with records_path.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("task_context", {}).get("task_type") != "docking":
                continue
            split = record.get("extensions", {}).get("dataset_split", {})
            if split.get("split") != "memory":
                raise ValueError("docking index accepts only declared memory split records")
            case_id = record["identity"]["case_id"]
            group_id = split.get("group_id")
            if not isinstance(case_id, str) or not case_id or case_id in seen:
                raise ValueError(f"invalid or duplicate docking case: {case_id!r}")
            if not isinstance(group_id, str) or not group_id:
                raise ValueError(f"missing protein group: {case_id}")
            seen.add(case_id)
            case_dir = (inputs_root / case_id).resolve()
            if inputs_root not in case_dir.parents:
                raise ValueError(f"case outside memory inputs: {case_id}")
            try:
                features = extract_query_features(case_dir / "ligand_00.sdf", case_dir / "embeddings.npz")
                metadata_path = case_dir / "metadata.json"
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                encoder_id = encoder_id_from_metadata(metadata)
                encoder_ids.add(encoder_id)
                binding = validate_query_metadata(metadata_path, case_dir / "embeddings.npz",
                                                  {"encoder_compatibility_id": encoder_id})
                if binding["case_id"] != case_id:
                    raise ValueError("query metadata case identity does not match memory record")
                record_encoder_ids = {item.get("encoder", {}).get("compatibility_id")
                                      for item in record.get("embeddings", [])}
                if encoder_id not in record_encoder_ids:
                    raise ValueError("metadata encoder does not match frozen memory record")
                ligand_hash = _digest(case_dir / "ligand_00.sdf")
                ligand_refs = {item.get("artifact_id") for item in record.get("inputs", {}).get("ligand_refs", [])}
                if f"sha256:{ligand_hash}" not in ligand_refs:
                    raise ValueError("ligand input does not match the frozen memory record")
                if not any(item.get("artifact_ref", {}).get("artifact_id") == f"sha256:{binding['embedding_sha256']}"
                           and item.get("source_input_sha256") == binding["source_input_sha256"]
                           for item in record.get("embeddings", [])):
                    raise ValueError("embedding or source receptor does not match the frozen memory record")
            except (OSError, ValueError, KeyError, IndexError, struct.error, zipfile.BadZipFile) as exc:
                raise ValueError(f"invalid memory features for {case_id}: {exc}") from exc
            metric_labels = _published_metrics(record)
            cases.append({"case_id": case_id, "group_id": group_id, **features,
                          "labels": metric_labels[DEFAULT_METRIC], "metric_labels": metric_labels})
            for filename in ("ligand_00.sdf", "embeddings.npz", "metadata.json"):
                inputs_digest.update(f"{case_id}/{filename}\0{_digest(case_dir / filename)}\n".encode())
    cases.sort(key=lambda item: item["case_id"])
    if len(cases) < 3:
        raise ValueError("at least three memory cases are required")
    if len(encoder_ids) != 1:
        raise ValueError("memory contains incompatible protein encoders")
    scales, priors = _statistics(cases)
    metric_priors = {metric: _statistics(cases, metric)[1] for metric in METRICS}
    index = {"schema_version": FORMAT_VERSION, "feature_version": FEATURE_VERSION,
             "policy_version": POLICY_VERSION, "memory_version": source_hash,
             "published_release": PUBLISHED_RELEASE, "parameters": PARAMETERS,
             "encoder_compatibility_id": next(iter(encoder_ids)),
             "source": {"bundle": str(bundle), "records": str(records_path),
                        "records_sha256": source_hash, "feature_inputs_sha256": inputs_digest.hexdigest(),
                        "release_manifest_sha256": _digest(manifest_path) if manifest is not None else None,
                        "release_manifest_verified": manifest is not None,
                        "split": "memory"},
             "tools": list(TOOLS), "ligand_scales": scales, "tool_priors": priors,
             "metrics": METRICS, "rmsd_loss_cap": RMSD_LOSS_CAP,
             "metric_priors": metric_priors, "cases": cases}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                         prefix=".docking-index-", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(index, handle, separators=(",", ":"), allow_nan=False)
            handle.write("\n")
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return {"success": True, "status": "built", "path": str(output), "case_count": len(cases),
            "feature_version": FEATURE_VERSION, "memory_version": source_hash,
            "published_release": PUBLISHED_RELEASE, "tool_priors": priors,
            "metrics": list(METRICS), "metric_priors": metric_priors,
            "missing_labels": {tool: sum(case["labels"][tool] is None for case in cases) for tool in TOOLS}}


def _validate_index(index: dict) -> None:
    if not isinstance(index, dict) or index.get("schema_version") not in (FORMAT_VERSION, LEGACY_FORMAT_VERSION):
        raise ValueError("incompatible docking index schema")
    if index.get("feature_version") != FEATURE_VERSION or index.get("policy_version") != POLICY_VERSION:
        raise ValueError("incompatible docking index feature or policy version")
    if not isinstance(index.get("encoder_compatibility_id"), str) or not index["encoder_compatibility_id"]:
        raise ValueError("missing encoder compatibility identity")
    if not isinstance(index.get("memory_version"), str) or not index["memory_version"]:
        raise ValueError("missing memory version")
    if index.get("published_release") != PUBLISHED_RELEASE or not isinstance(index.get("source"), dict) or index["source"].get("split") != "memory":
        raise ValueError("index is not bound to the frozen memory release")
    if index.get("tools") != list(TOOLS) or index.get("parameters") != PARAMETERS:
        raise ValueError("incompatible docking index policy parameters or tools")
    multi_metric = index["schema_version"] == FORMAT_VERSION
    if multi_metric and (index.get("metrics") != METRICS or index.get("rmsd_loss_cap") != RMSD_LOSS_CAP):
        raise ValueError("incompatible docking index metrics or RMSD loss transform")
    scales = _finite_vector(index.get("ligand_scales"), 18, "descriptor scales")
    if any(value <= 0 for value in scales):
        raise ValueError("descriptor scales must be positive")
    priors = index.get("tool_priors")
    if not isinstance(priors, dict) or set(priors) != set(TOOLS) or any(
        value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))
                               or not math.isfinite(value) or not 0 <= value <= 1)
        for value in priors.values()
    ):
        raise ValueError("invalid tool priors")
    if multi_metric:
        metric_priors = index.get("metric_priors")
        if not isinstance(metric_priors, dict) or set(metric_priors) != set(METRICS):
            raise ValueError("invalid metric priors")
        for metric, kind in METRICS.items():
            values = metric_priors[metric]
            upper = RMSD_LOSS_CAP if kind == "rmsd" else 1
            if not isinstance(values, dict) or set(values) != set(TOOLS) or any(
                value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))
                                       or not math.isfinite(value) or not 0 <= value <= upper)
                for value in values.values()
            ):
                raise ValueError(f"invalid metric priors for {metric}")
        if metric_priors[DEFAULT_METRIC] != priors:
            raise ValueError("default metric priors disagree with compatibility priors")
    cases = index.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("empty or invalid docking index")
    seen = set()
    for case in cases:
        _validate_features(case)
        case_id, group_id = case.get("case_id"), case.get("group_id")
        if not isinstance(case_id, str) or not case_id or case_id in seen or not isinstance(group_id, str) or not group_id:
            raise ValueError("invalid case identity or protein group")
        seen.add(case_id)
        labels = case.get("labels")
        if not isinstance(labels, dict) or set(labels) != set(TOOLS) or any(
            value is not None and type(value) is not bool for value in labels.values()
        ):
            raise ValueError("labels must be boolean success or null for missing")
        if multi_metric:
            metrics = case.get("metric_labels")
            if not isinstance(metrics, dict) or set(metrics) != set(METRICS):
                raise ValueError("invalid case metric labels")
            for metric, kind in METRICS.items():
                values = metrics[metric]
                if not isinstance(values, dict) or set(values) != set(TOOLS):
                    raise ValueError(f"invalid case labels for {metric}")
                for value in values.values():
                    if value is None:
                        continue
                    if kind == "binary" and type(value) is not bool:
                        raise ValueError(f"{metric} labels must be boolean or null")
                    if kind == "rmsd" and (isinstance(value, bool) or not isinstance(value, (int, float))
                                           or not math.isfinite(value) or value < 0):
                        raise ValueError(f"{metric} labels must be finite nonnegative RMSD or null")
            if metrics[DEFAULT_METRIC] != labels:
                raise ValueError("default metric labels disagree with compatibility labels")
            for tool in TOOLS:
                top1, success = metrics["top1_rmsd_a"][tool], labels[tool]
                pb, joint = metrics["posebusters_valid"][tool], metrics["joint_success"][tool]
                if top1 is not None and success != (top1 < 2.0):
                    raise ValueError("Top-1 success label disagrees with raw RMSD")
                if success is not None and pb is not None and joint is not None and joint != (success and pb):
                    raise ValueError("joint success label disagrees with Top-1 success and PoseBusters")


def load_index(path: Path) -> dict:
    index = json.loads(Path(path).read_text(encoding="utf-8"))
    _validate_index(index)
    return index


def task_policy(task: str = "redocking", metric: str | None = None) -> dict:
    """Return the task's baseline, objective and guarded memory settings.

    Pose generation uses Top-1 PoseBusters labels from the frozen observations
    and shares the redocking prior strength and switching margin.
    """
    if task not in TASK_DEFAULTS:
        raise ValueError("selection task must be redocking or pose_generation")
    defaults = TASK_DEFAULTS[task]
    metric = defaults["default_metric"] if metric is None else metric
    inferred = task == "pose_generation"
    return {
        "task": task, "metric": metric, **defaults,
        "encoder_task": "docking",
        "policy_basis": "ligand_retrieval_protein_reweight_guarded_switch",
        "objective_scope": (
            "explicit_extension" if metric != defaults["default_metric"]
            else "inferred_pose_generation" if inferred else "redocking_success"
        ),
        "assumptions": ([
            "The default pose-generation objective uses Top-1 PoseBusters-valid labels from the frozen docking observations.",
            "Pose generation uses the rDock baseline and shares redocking prior_strength=2.0 and switch_margin=0.10.",
        ] if inferred else []),
        "reproduces_published_pose_quality_benchmark": False,
    }


def _available_tools(available_tools) -> set[str]:
    available = set(AVAILABLE_TOOLS if available_tools is None else available_tools)
    # ADFR is executable but has no frozen historical labels; it cannot be a
    # memory recommendation. AutoDock-GPU has no current executable adapter.
    unknown = available - (set(AVAILABLE_TOOLS) | {"adfr"})
    if unknown:
        raise ValueError("unsupported available docking engines: " + ", ".join(sorted(unknown)))
    return available


def fallback_decision(reason: str, available_tools=None, mode: str = "ligand_protein",
                      metric: str | None = None, task: str = "redocking") -> dict:
    policy = task_policy(task, metric)
    metric, baseline = policy["metric"], policy["default_tool"]
    available = _available_tools(available_tools)
    return {"success": baseline in available, "status": "fallback" if baseline in available else "unavailable",
            "mode": mode, "metric": metric, "metric_kind": METRICS.get(metric),
            "task": task, "task_policy": policy, "default_tool": baseline,
            "policy_version": POLICY_VERSION if task == "redocking" else "pose_generation_inferred_v1",
            "index_policy_version": POLICY_VERSION, "feature_version": FEATURE_VERSION,
            "parameters": {**PARAMETERS, "default_tool": baseline,
                           "protein_beta": 0.10 if mode == "ligand_protein" else 0.0},
            "memory_version": None, "selected_tool": baseline if baseline in available else None,
            "proposed_tool": None, "fallback_reason": reason, "scores": {}, "neighbors": []}


def select_docking_tool(features: dict, index: dict, mode: str = "ligand_protein",
                        available_tools=None, exclude_cases=(), exclude_groups=(),
                        metric: str | None = None, task: str = "redocking") -> dict:
    """Rank all historical tools, then gate the selected tool by availability.

    Exclusions are applied before ranks, priors and descriptor scales are fitted.
    Missing labels are masked independently per tool, never counted as failures.
    """
    available = _available_tools(available_tools)
    decision = fallback_decision("invalid_input", available, mode, metric, task)
    metric, baseline = decision["metric"], decision["default_tool"]
    if metric not in METRICS:
        decision["fallback_reason"] = "unsupported_metric"
        return decision
    if mode not in ("ligand_protein", "ligand_only"):
        decision["fallback_reason"] = "unsupported_mode"
        return decision
    try:
        _validate_index(index)
        _validate_features(features, require_protein=mode == "ligand_protein")
    except (ValueError, TypeError, KeyError, OverflowError) as exc:
        decision["detail"] = str(exc)
        return decision
    decision["memory_version"] = index["memory_version"]
    decision["published_release"] = index["published_release"]
    if metric != DEFAULT_METRIC and index["schema_version"] == LEGACY_FORMAT_VERSION:
        decision["fallback_reason"] = "unsupported_metric_in_legacy_index"
        decision["detail"] = "Rebuild the index from the frozen memory partition for this objective"
        return decision
    binary = METRICS[metric] == "binary"
    decision["score_direction"] = "higher_is_better" if binary else "lower_is_better"
    decision["score_units"] = "success_fraction" if binary else "angstrom"
    decision["rmsd_loss_cap"] = None if binary else RMSD_LOSS_CAP
    beta = 0.10 if mode == "ligand_protein" else 0.0
    decision["parameters"] = {**PARAMETERS, "protein_beta": beta, "default_tool": baseline}
    excluded_cases, excluded_groups = set(exclude_cases), set(exclude_groups)
    excluded_groups.update(case["group_id"] for case in index["cases"] if case["case_id"] in excluded_cases)
    decision["excluded_cases"] = sorted(excluded_cases)
    decision["excluded_groups"] = sorted(excluded_groups)
    cases = [case for case in index["cases"] if case["case_id"] not in excluded_cases
             and case["group_id"] not in excluded_groups]
    cases.sort(key=lambda case: case["case_id"])
    decision["candidate_count"] = len(cases)
    if len(cases) < 3:
        decision["fallback_reason"] = "insufficient_neighbors"
        return decision
    if excluded_cases or excluded_groups:
        scales, priors = _statistics(cases, metric)
    else:
        scales = index["ligand_scales"]
        priors = index["tool_priors"] if metric == DEFAULT_METRIC else index["metric_priors"][metric]
    fingerprint = set(features["ligand_fp"])
    similarities = []
    for case in cases:
        candidate_fp = set(case["ligand_fp"])
        graph = len(fingerprint & candidate_fp) / len(fingerprint | candidate_fp)
        distance = sum(((x - y) / scale) ** 2 for x, y, scale in
                       zip(features["ligand_desc"], case["ligand_desc"], scales)) / 18
        ligand = 0.8 * graph + 0.2 * math.exp(-0.5 * distance)
        protein = sum(x * y for x, y in zip(features["protein"], case["protein"])) if beta else None
        similarities.append((protein, ligand))
    protein_order = sorted(range(len(cases)), key=lambda i: (similarities[i][0], -i)) if beta else []
    protein_rank = {i: rank / max(len(cases) - 1, 1) for rank, i in enumerate(protein_order)}
    neighbors = sorted(range(len(cases)), key=lambda i: (-similarities[i][1], i))[:3]
    weights = [math.exp(beta * (protein_rank.get(i, 0.5) - 0.5)) / math.log2(rank + 2)
               for rank, i in enumerate(neighbors)]
    decision["neighbors"] = [{"case_id": cases[i]["case_id"], "group_id": cases[i]["group_id"],
                              "ligand_similarity": similarities[i][1], "protein_similarity": similarities[i][0],
                              "protein_percentile": protein_rank.get(i), "ligand_rank": rank,
                              "weight": weight, "labels": _case_labels(cases[i], metric)}
                             for rank, (i, weight) in enumerate(zip(neighbors, weights))]
    scores = {}
    labeled_neighbor_count = {}
    for tool in TOOLS:
        observed = [(_loss(_case_labels(cases[i], metric)[tool], metric), weight)
                    for i, weight in zip(neighbors, weights) if _case_labels(cases[i], metric)[tool] is not None]
        labeled_neighbor_count[tool] = len(observed)
        alpha = decision["parameters"]["prior_strength"]
        scores[tool] = ((alpha * priors[tool] + sum(label * weight for label, weight in observed))
                        / (alpha + sum(weight for _, weight in observed))) if priors[tool] is not None else None
    decision["scores"] = scores
    decision["labeled_neighbor_count"] = labeled_neighbor_count
    decision["tool_priors"] = priors
    if scores[baseline] is None:
        decision["fallback_reason"] = f"missing_{baseline}_labels"
        return decision
    scored_tools = [tool for tool in TOOLS if scores[tool] is not None]
    if binary:
        proposed = max(scored_tools, key=lambda tool: (scores[tool], tool == baseline, tool))
        advantage = scores[proposed] - scores[baseline]
    else:
        proposed = min(scored_tools, key=lambda tool: (scores[tool], tool != baseline, tool))
        advantage = scores[baseline] - scores[proposed]
    selected = proposed if advantage > decision["parameters"]["switch_margin"] else baseline
    decision.update({"proposed_tool": proposed, "advantage_over_default": advantage})
    if baseline == "gnina":
        decision["advantage_over_gnina"] = advantage
    if selected not in available:
        decision["fallback_reason"] = "selected_tool_unavailable"
        decision["unavailable_tool"] = selected
        return decision
    decision.update({"success": True, "status": "selected", "selected_tool": selected, "fallback_reason": None,
                     "selection_reason": "margin_exceeded" if selected != baseline else f"{baseline}_margin_guard"})
    return decision



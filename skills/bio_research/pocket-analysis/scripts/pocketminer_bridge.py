#!/usr/bin/env python3
"""Strict PocketMiner worker executed inside the dedicated Conda environment.

PocketMiner predicts a cryptic-opening propensity for each protein residue.  It
does not predict a cavity volume or a docking box.  This worker deliberately
does not use upstream ``util.load_checkpoint`` because that helper discards the
TensorFlow restore status and can therefore report success with random weights.

The module keeps TensorFlow imports inside the inference path so its PDB parser
and JSON contract can be tested from the Drugmate environment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Any


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from deployment_config import deployment_path


EXPECTED_REPOSITORY_COMMIT = "187062df3c94127e991669768009141a08fd5d8b"
EXPECTED_CHECKPOINT_SHA256 = {
    "pocketminer.index": "09e36a62a987b14bc49cf7a1fa53fa734f85d85d088777f98c6aca253b47bf46",
    "pocketminer.data-00000-of-00001": (
        "6f5ab62b9fb38b54040053ac321f3362d616bae9b31569249a032dbb50890b70"
    ),
}
BACKBONE_ATOMS = ("N", "CA", "C", "O")
STANDARD_RESIDUES = {
    "ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4,
    "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8, "ILE": 9,
    "LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14,
    "SER": 15, "THR": 16, "TRP": 17, "TYR": 18, "VAL": 19,
}
PARENT_RESIDUE = {
    "ASH": "ASP", "CYM": "CYS", "CYX": "CYS", "GLH": "GLU",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS",
    "HSE": "HIS", "HSP": "HIS", "LYN": "LYS", "MSE": "MET",
    "SEC": "CYS",
}
HETATM_POLYMER_RESIDUES = {"MSE", "SEC"}


class PocketMinerBridgeError(RuntimeError):
    """An input, checkpoint, or inference invariant was violated."""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def _canonical_record(chain: str, resname: str, resnum: int, icode: str,
                      inference_resname: str) -> dict[str, Any]:
    legacy = f"{resname}{resnum}{icode}"
    return {
        "chain": chain,
        "resname": resname,
        "resnum": resnum,
        "icode": icode,
        "label": f"{chain}:{legacy}" if chain else legacy,
        "legacy_label": legacy,
        "inference_resname": inference_resname,
    }


def _atom_priority(altloc: str, occupancy: float) -> tuple[int, float]:
    """Prefer blank then A altloc, using occupancy only within that class."""
    altloc = altloc.strip()
    rank = 0 if not altloc else 1 if altloc == "A" else 2
    return rank, -occupancy


def parse_chain_backbone(protein_file: Path, chain: str) -> list[dict[str, Any]]:
    """Parse one PDB chain into an ordered, chain/iCode-safe residue manifest.

    Coordinates are retained in Angstrom here and converted to nm only when
    constructing the model tensor, matching MDTraj's unit conversion.
    """
    if not protein_file.is_file():
        raise PocketMinerBridgeError(
            "input_missing", f"protein PDB does not exist: {protein_file}"
        )

    residues: OrderedDict[tuple[str, int, str], dict[str, Any]] = OrderedDict()
    first_model_started = False
    with protein_file.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            record_type = line[0:6].strip().upper()
            if record_type == "MODEL":
                if first_model_started:
                    break
                first_model_started = True
                continue
            if record_type == "ENDMDL":
                break
            if record_type not in {"ATOM", "HETATM"}:
                continue
            if len(line) < 54:
                raise PocketMinerBridgeError(
                    "invalid_pdb_record", f"short coordinate record at line {line_number}"
                )

            pdb_chain = line[21:22].strip()
            if pdb_chain != chain:
                continue
            resname = line[17:20].strip().upper()
            if record_type == "HETATM" and resname not in HETATM_POLYMER_RESIDUES:
                continue
            inference_resname = PARENT_RESIDUE.get(resname, resname)
            if inference_resname not in STANDARD_RESIDUES:
                raise PocketMinerBridgeError(
                    "unsupported_residue",
                    f"PocketMiner cannot encode {resname!r} at PDB line {line_number}",
                )
            try:
                resnum = int(line[22:26].strip())
                xyz = tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
                occupancy_text = line[54:60].strip()
                occupancy = float(occupancy_text) if occupancy_text else 0.0
            except ValueError as exc:
                raise PocketMinerBridgeError(
                    "invalid_pdb_record", f"invalid numeric field at line {line_number}: {exc}"
                ) from exc
            if not all(math.isfinite(value) for value in xyz):
                raise PocketMinerBridgeError(
                    "invalid_coordinates", f"non-finite coordinate at line {line_number}"
                )

            icode = line[26:27].strip()
            atom_name = line[12:16].strip().upper()
            if atom_name not in BACKBONE_ATOMS:
                continue
            key = (pdb_chain, resnum, icode)
            residue = residues.setdefault(
                key,
                {
                    "record": _canonical_record(
                        pdb_chain, resname, resnum, icode, inference_resname
                    ),
                    "atoms": {},
                },
            )
            prior_record = residue["record"]
            if prior_record["resname"] != resname:
                raise PocketMinerBridgeError(
                    "ambiguous_residue_identity",
                    f"multiple residue names at {pdb_chain or '<blank>'}:{resnum}{icode}",
                )
            candidate = {
                "coord": xyz,
                "priority": _atom_priority(line[16:17], occupancy),
            }
            current = residue["atoms"].get(atom_name)
            if current is None or candidate["priority"] < current["priority"]:
                residue["atoms"][atom_name] = candidate

    if not residues:
        raise PocketMinerBridgeError(
            "chain_not_found",
            f"no supported protein residues found for chain {chain!r} in {protein_file}",
        )

    incomplete = []
    manifest = []
    for residue in residues.values():
        missing = [name for name in BACKBONE_ATOMS if name not in residue["atoms"]]
        if missing:
            incomplete.append(f"{residue['record']['label']}({','.join(missing)})")
            continue
        manifest.append({
            **residue["record"],
            "backbone_coordinates_angstrom": [
                list(residue["atoms"][name]["coord"]) for name in BACKBONE_ATOMS
            ],
        })
    if incomplete:
        preview = ", ".join(incomplete[:20])
        suffix = " ..." if len(incomplete) > 20 else ""
        raise PocketMinerBridgeError(
            "incomplete_backbone",
            f"target chain has residues missing N/CA/C/O: {preview}{suffix}",
        )
    return manifest


def _repository_commit(source_dir: Path) -> str:
    try:
        completed = subprocess.run(
            ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PocketMinerBridgeError(
            "repository_unverifiable", f"could not verify PocketMiner commit: {exc}"
        ) from exc
    commit = (completed.stdout or "").strip()
    if completed.returncode != 0 or len(commit) != 40:
        detail = (completed.stderr or completed.stdout or "git rev-parse failed").strip()
        raise PocketMinerBridgeError("repository_unverifiable", detail)
    if commit != EXPECTED_REPOSITORY_COMMIT:
        raise PocketMinerBridgeError(
            "repository_commit_mismatch",
            f"expected {EXPECTED_REPOSITORY_COMMIT}, found {commit}",
        )
    dirty = subprocess.run(
        [
            "git", "-C", str(source_dir), "status", "--porcelain",
            "--untracked-files=all", "--", "src", "models/pocketminer.index",
            "models/pocketminer.data-00000-of-00001",
        ],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    if dirty.returncode != 0:
        detail = (dirty.stderr or dirty.stdout or "git status failed").strip()
        raise PocketMinerBridgeError("repository_unverifiable", detail)
    dirty_entries = [
        line for line in dirty.stdout.splitlines()
        if not _is_generated_python_cache(line)
    ]
    if dirty_entries:
        raise PocketMinerBridgeError(
            "repository_dirty",
            "PocketMiner source/checkpoint files differ from the pinned commit",
        )
    return commit


def _is_generated_python_cache(status_line: str) -> bool:
    """Allow only import-generated, untracked bytecode in the pinned source tree."""
    if not status_line.startswith("?? "):
        return False
    path = status_line[3:].strip()
    return "/__pycache__/" in f"/{path}" and path.endswith((".pyc", ".pyo"))


def _validate_installation(source_dir: Path, checkpoint_prefix: Path) -> dict[str, Any]:
    source_file = source_dir / "src" / "models.py"
    index_file = Path(f"{checkpoint_prefix}.index")
    data_files = sorted(checkpoint_prefix.parent.glob(f"{checkpoint_prefix.name}.data-*"))
    missing = [str(path) for path in (source_file, index_file) if not path.is_file()]
    if not data_files:
        missing.append(f"{checkpoint_prefix}.data-*")
    if missing:
        raise PocketMinerBridgeError(
            "installation_incomplete", "missing PocketMiner resource(s): " + ", ".join(missing)
        )
    checkpoint_files = [index_file, *data_files]
    observed_hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in checkpoint_files
    }
    if observed_hashes != EXPECTED_CHECKPOINT_SHA256:
        raise PocketMinerBridgeError(
            "checkpoint_hash_mismatch",
            "checkpoint files do not match the pinned PocketMiner assets "
            f"(observed={observed_hashes})",
        )
    return {
        "repository_commit": _repository_commit(source_dir),
        "checkpoint_prefix": str(checkpoint_prefix.resolve()),
        "checkpoint_index": str(index_file.resolve()),
        "checkpoint_data_files": [str(path.resolve()) for path in data_files],
        "checkpoint_sha256": observed_hashes,
    }


def _synthetic_inputs(np_module):
    length = 32
    x = np_module.zeros((1, length, 4, 3), dtype=np_module.float32)
    for index in range(length):
        theta = index * 1.1
        base = np_module.asarray(
            [0.15 * index, 0.25 * math.sin(theta), 0.25 * math.cos(theta)],
            dtype=np_module.float32,
        )
        x[0, index] = np_module.asarray([
            base + [-0.05, 0.02, 0.00],
            base,
            base + [0.05, 0.02, 0.01],
            base + [0.08, 0.05, 0.01],
        ])
    sequence = np_module.asarray(
        [[index % 20 for index in range(length)]], dtype=np_module.int32
    )
    mask = np_module.ones((1, length), dtype=np_module.float32)
    return x, sequence, mask


def _manifest_inputs(manifest: list[dict[str, Any]], np_module):
    coordinates = np_module.asarray(
        [record["backbone_coordinates_angstrom"] for record in manifest],
        dtype=np_module.float32,
    ) / np_module.float32(10.0)
    sequence = np_module.asarray(
        [STANDARD_RESIDUES[record["inference_resname"]] for record in manifest],
        dtype=np_module.int32,
    )
    return (
        coordinates[np_module.newaxis, ...],
        sequence[np_module.newaxis, ...],
        np_module.ones((1, len(manifest)), dtype=np_module.float32),
    )


def _validate_repeated_predictions(first, second, np_module) -> float:
    """Validate both inference passes and return their maximum absolute delta."""
    if first.shape != second.shape or first.ndim != 2 or first.shape[0] != 1:
        raise PocketMinerBridgeError(
            "prediction_shape_mismatch", f"unexpected prediction shape {first.shape}"
        )
    if not np_module.isfinite(first).all() or not np_module.isfinite(second).all():
        raise PocketMinerBridgeError(
            "nonfinite_predictions", "one or both repeated predictions contain NaN or infinity"
        )
    if (
        ((first < 0.0) | (first > 1.0)).any()
        or ((second < 0.0) | (second > 1.0)).any()
    ):
        raise PocketMinerBridgeError(
            "prediction_range_error", "one or both repeated sigmoid predictions are outside [0, 1]"
        )
    repeat_delta = float(np_module.max(np_module.abs(first - second))) if first.size else 0.0
    if not math.isfinite(repeat_delta) or repeat_delta > 1e-7:
        raise PocketMinerBridgeError(
            "nondeterministic_inference",
            f"two train=False forwards differ (max_abs_delta={repeat_delta:.9g})",
        )
    return repeat_delta


def _load_model_and_predict(source_dir: Path, checkpoint_prefix: Path, inputs):
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    import numpy as np
    import tensorflow as tf

    tf.random.set_seed(1729)
    np.random.seed(1729)
    try:
        tf.config.experimental.enable_op_determinism()
    except (AttributeError, RuntimeError):
        pass

    source_path = str((source_dir / "src").resolve())
    if source_path not in sys.path:
        sys.path.insert(0, source_path)
    try:
        from models import MQAModel
    except Exception as exc:
        raise PocketMinerBridgeError(
            "source_import_failed", f"could not import PocketMiner model: {exc}"
        ) from exc

    x, sequence, mask = inputs(np) if callable(inputs) else inputs
    model = MQAModel(
        node_features=(8, 50), edge_features=(1, 32),
        hidden_dim=(16, 100), num_layers=4, dropout=0.1,
    )

    # Materialize every model variable before reading the object checkpoint.
    # The checkpoint was created with Checkpoint.write(), so read() is the
    # matching API; optimizer-only entries may remain unused, model entries may not.
    try:
        model(x, sequence, mask, train=False, res_level=True)
        status = tf.train.Checkpoint(model=model).read(str(checkpoint_prefix))
        if hasattr(status, "assert_nontrivial_match"):
            status.assert_nontrivial_match()
        status.assert_existing_objects_matched()
    except (AssertionError, ValueError, tf.errors.OpError) as exc:
        detail = " ".join(str(exc).split())[:4000]
        raise PocketMinerBridgeError(
            "checkpoint_existing_objects_mismatch",
            "PocketMiner checkpoint did not populate the complete model graph: " + detail,
        ) from exc

    # assert_existing_objects_matched() proves every currently built object was
    # restored, but it permits checkpoint-only entries. Explicitly account for
    # every checkpoint model variable before allowing optimizer slots to remain
    # partial. Together, object matching + an exact model key-to-shape map from
    # a freshly written current checkpoint establish bidirectional coverage.
    checkpoint_variables = tf.train.list_variables(str(checkpoint_prefix))
    checkpoint_model_variables = {
        name: tuple(shape)
        for name, shape in checkpoint_variables
        if name.startswith("model/")
        and name.endswith("/.ATTRIBUTES/VARIABLE_VALUE")
        and "/.OPTIMIZER_SLOT/" not in name
    }
    with tempfile.TemporaryDirectory(prefix="pocketminer-model-graph-") as temp_dir:
        current_prefix = str(Path(temp_dir) / "current_model")
        tf.train.Checkpoint(model=model).write(current_prefix)
        current_model_variables = {
            name: tuple(shape)
            for name, shape in tf.train.list_variables(current_prefix)
            if name.startswith("model/")
            and name.endswith("/.ATTRIBUTES/VARIABLE_VALUE")
            and "/.OPTIMIZER_SLOT/" not in name
        }
    if checkpoint_model_variables != current_model_variables:
        missing = sorted(checkpoint_model_variables.keys() - current_model_variables.keys())
        extra = sorted(current_model_variables.keys() - checkpoint_model_variables.keys())
        raise PocketMinerBridgeError(
            "checkpoint_model_coverage_mismatch",
            "checkpoint/current model key-to-shape maps differ "
            f"(checkpoint={len(checkpoint_model_variables)}, current={len(current_model_variables)}, "
            f"missing_current_keys={missing[:10]}, extra_current_keys={extra[:10]})",
        )
    unexpected_checkpoint_keys = [
        name for name, _ in checkpoint_variables
        if name != "_CHECKPOINTABLE_OBJECT_GRAPH"
        and not name.startswith("model/")
        and not name.startswith("optimizer/")
        and "/.OPTIMIZER_SLOT/" not in name
    ]
    if unexpected_checkpoint_keys:
        raise PocketMinerBridgeError(
            "checkpoint_unclassified_entries",
            "checkpoint contains unclassified entries: "
            + ", ".join(unexpected_checkpoint_keys[:20]),
        )
    status.expect_partial()  # only optimizer state is permitted to remain unused

    first = np.asarray(
        model(x, sequence, mask, train=False, res_level=True).numpy(),
        dtype=np.float64,
    )
    second = np.asarray(
        model(x, sequence, mask, train=False, res_level=True).numpy(),
        dtype=np.float64,
    )
    repeat_delta = _validate_repeated_predictions(first, second, np)
    return first[0], {
        "tensorflow_version": tf.__version__,
        "numpy_version": np.__version__,
        "model_weight_objects": len(model.weights),
        "model_parameter_count": int(model.count_params()),
        "checkpoint_existing_objects_matched": True,
        "checkpoint_model_variable_count": len(checkpoint_model_variables),
        "checkpoint_model_key_shape_map_match": True,
        "unused_checkpoint_entries_policy": "optimizer_state_only",
        "all_scores_finite": True,
        "scores_in_unit_interval": True,
        "repeat_max_abs_delta": repeat_delta,
    }


def _percentile(values: list[float], quantile: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = quantile * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _cluster_hypotheses(scored_records: list[dict[str, Any]], quantile: float,
                        cutoff: float, min_residues: int,
                        max_hypotheses: int) -> tuple[list[dict], float]:
    scores = [float(record["score"]) for record in scored_records]
    threshold = _percentile(scores, quantile)
    selected = [record for record in scored_records if float(record["score"]) >= threshold]
    remaining = set(range(len(selected)))
    components = []
    cutoff_squared = cutoff * cutoff
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        component = [seed]
        frontier = [seed]
        while frontier:
            left = frontier.pop()
            left_ca = selected[left]["ca_coordinate_angstrom"]
            neighbors = []
            for right in remaining:
                right_ca = selected[right]["ca_coordinate_angstrom"]
                distance_squared = sum(
                    (float(left_ca[axis]) - float(right_ca[axis])) ** 2 for axis in range(3)
                )
                if distance_squared <= cutoff_squared:
                    neighbors.append(right)
            for right in neighbors:
                remaining.remove(right)
                frontier.append(right)
                component.append(right)
        if len(component) >= min_residues:
            components.append([selected[index] for index in sorted(component)])

    hypotheses = []
    for component in components:
        center = [
            round(sum(float(record["ca_coordinate_angstrom"][axis]) for record in component)
                  / len(component), 4)
            for axis in range(3)
        ]
        hypotheses.append({
            "residue_records": [
                {key: value for key, value in record.items() if key != "ca_coordinate_angstrom"}
                for record in component
            ],
            "residue_count": len(component),
            "mean_score": round(sum(float(record["score"]) for record in component)
                                / len(component), 8),
            "max_score": round(max(float(record["score"]) for record in component), 8),
            "visualization_centroid": center,
            "geometry_status": "provisional_residue_cluster_only",
            "cavity_volume": None,
            "docking_box": None,
        })
    hypotheses.sort(key=lambda item: (-item["mean_score"], -item["residue_count"]))
    hypotheses = hypotheses[:max_hypotheses]
    for index, hypothesis in enumerate(hypotheses, start=1):
        hypothesis["hypothesis_id"] = f"PMH{index:03d}"
    return hypotheses, threshold


def run_probe(source_dir: Path, checkpoint_prefix: Path) -> dict[str, Any]:
    installation = _validate_installation(source_dir, checkpoint_prefix)
    _, validation = _load_model_and_predict(
        source_dir, checkpoint_prefix, _synthetic_inputs,
    )
    return {
        "schema_version": 1,
        "success": True,
        "mode": "probe",
        "producer": {"name": "PocketMiner", **installation},
        "validation": validation,
    }


def run_prediction(args) -> dict[str, Any]:
    source_dir = Path(args.source_dir).resolve()
    checkpoint_prefix = Path(args.checkpoint).resolve()
    protein_file = Path(args.protein_file).resolve()
    installation = _validate_installation(source_dir, checkpoint_prefix)
    manifest = parse_chain_backbone(protein_file, args.chain)

    import numpy as np
    predictions, validation = _load_model_and_predict(
        source_dir,
        checkpoint_prefix,
        _manifest_inputs(manifest, np),
    )
    if len(predictions) != len(manifest):
        raise PocketMinerBridgeError(
            "prediction_count_mismatch",
            f"received {len(predictions)} scores for {len(manifest)} residues",
        )

    scored_records = []
    for index, (record, score) in enumerate(zip(manifest, predictions)):
        scored_records.append({
            "index": index,
            **{key: value for key, value in record.items()
               if key != "backbone_coordinates_angstrom"},
            "score": round(float(score), 10),
            "ca_coordinate_angstrom": list(record["backbone_coordinates_angstrom"][1]),
        })
    hypotheses, heuristic_threshold = _cluster_hypotheses(
        scored_records,
        quantile=args.provisional_quantile,
        cutoff=args.cluster_ca_cutoff,
        min_residues=args.min_cluster_residues,
        max_hypotheses=args.max_hypotheses,
    )
    public_scores = [
        {key: value for key, value in record.items() if key != "ca_coordinate_angstrom"}
        for record in scored_records
    ]
    score_values = [record["score"] for record in public_scores]
    validation["prediction_count_matches"] = True
    return {
        "schema_version": 1,
        "success": True,
        "mode": "predict",
        "producer": {"name": "PocketMiner", **installation},
        "input": {
            "protein_file": str(protein_file),
            "protein_sha256": hashlib.sha256(protein_file.read_bytes()).hexdigest(),
            "chain": args.chain,
            "residue_count": len(manifest),
        },
        "validation": validation,
        "score_semantics": {
            "meaning": "per-residue cryptic-pocket opening propensity",
            "official_fixed_inference_threshold": None,
            "note": "PocketMiner does not provide cavity geometry, volume, or a docking box.",
        },
        "score_summary": {
            "minimum": round(min(score_values), 10),
            "maximum": round(max(score_values), 10),
            "mean": round(sum(score_values) / len(score_values), 10),
        },
        "residue_scores": public_scores,
        "provisional_cryptic_hypotheses": hypotheses,
        "provisional_selection_policy": {
            "kind": "within_structure_score_quantile_then_ca_connected_components",
            "quantile": args.provisional_quantile,
            "computed_score_cutoff": round(heuristic_threshold, 10),
            "cluster_ca_cutoff_angstrom": args.cluster_ca_cutoff,
            "min_cluster_residues": args.min_cluster_residues,
            "max_hypotheses": args.max_hypotheses,
            "status": "pipeline_heuristic_not_an_upstream_calibrated_threshold",
            "promotion_policy": "never_publish_as_canonical_pocket_or_docking_box_without_cavity_geometry",
        },
    }


def _write_payload(payload: dict[str, Any], output_json: str | None) -> None:
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, default=str)
    if output_json:
        destination = Path(output_json)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name + ".tmp")
        temporary.write_text(rendered + "\n", encoding="utf-8")
        os.replace(temporary, destination)
    print(rendered)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Strict PocketMiner pipeline bridge")
    subparsers = parser.add_subparsers(dest="mode", required=True)
    for mode in ("probe", "predict"):
        command = subparsers.add_parser(mode)
        command.add_argument("--source-dir", default=str(deployment_path("${SBDD_POCKETMINER_SOURCE_DIR}")))
        command.add_argument("--checkpoint", default=str(deployment_path("${SBDD_POCKETMINER_SOURCE_DIR}/models/pocketminer")))
        command.add_argument("--output-json")
        if mode == "predict":
            command.add_argument("--protein-file", required=True)
            command.add_argument("--chain", default="A")
            command.add_argument("--provisional-quantile", type=float, default=0.9)
            command.add_argument("--cluster-ca-cutoff", type=float, default=8.0)
            command.add_argument("--min-cluster-residues", type=int, default=3)
            command.add_argument("--max-hypotheses", type=int, default=5)
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.mode == "predict":
        if not 0.0 < args.provisional_quantile < 1.0:
            parser.error("--provisional-quantile must be between 0 and 1")
        if args.cluster_ca_cutoff <= 0 or args.min_cluster_residues < 1 or args.max_hypotheses < 1:
            parser.error("clustering parameters must be positive")
    try:
        source_dir = Path(args.source_dir).resolve()
        checkpoint = Path(args.checkpoint).resolve()
        payload = (
            run_probe(source_dir, checkpoint)
            if args.mode == "probe" else run_prediction(args)
        )
        _write_payload(payload, args.output_json)
        return 0
    except PocketMinerBridgeError as exc:
        payload = {
            "schema_version": 1,
            "success": False,
            "mode": args.mode,
            "error_code": exc.code,
            "error": str(exc),
        }
        _write_payload(payload, args.output_json)
        return 2
    except Exception as exc:
        payload = {
            "schema_version": 1,
            "success": False,
            "mode": args.mode,
            "error_code": "unexpected_bridge_error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        _write_payload(payload, args.output_json)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())

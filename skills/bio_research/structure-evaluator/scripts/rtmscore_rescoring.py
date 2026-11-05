"""Pose-identity-safe RTMScore ensemble adapter for Step 5."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from statistics import fmean, pstdev
from typing import Any, Mapping

from tool_runtime import ToolRuntimeError, run_tool


SCHEMA_VERSION = "rtmscore-rescoring/v1"
WORKER = Path(__file__).resolve().with_name("rtmscore_worker.py")


class RTMScoreError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _load_and_relabel_poses(source: Path, destination: Path) -> list[dict]:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RTMScoreError("RDKit is required to bind RTMScore inputs") from exc

    supplier = Chem.SDMolSupplier(str(source), removeHs=False, sanitize=True)
    writer = Chem.SDWriter(str(destination))
    mapping: list[dict] = []
    seen_pose_ids: set[str] = set()
    try:
        for index, mol in enumerate(supplier):
            if mol is None:
                raise RTMScoreError(f"selected pose {index} cannot be parsed")
            if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
                raise RTMScoreError(f"selected pose {index} has no 3D conformer")
            pose_id = mol.GetProp("pose_id").strip() if mol.HasProp("pose_id") else ""
            if not pose_id:
                raise RTMScoreError(
                    f"selected pose {index} is missing canonical pose_id"
                )
            if pose_id in seen_pose_ids:
                raise RTMScoreError(
                    f"selected pose cohort contains duplicate pose_id {pose_id!r}"
                )
            seen_pose_ids.add(pose_id)
            token = f"RTMPOSE_{index:06d}"
            upstream_id = f"{token}-{index}"
            mol.SetProp("_Name", token)
            writer.write(mol)
            mapping.append({
                "record_index": index,
                "pose_id": pose_id,
                "compound_id": (
                    mol.GetProp("compound_id") if mol.HasProp("compound_id") else token
                ),
                "token": token,
                "upstream_id": upstream_id,
                "coordinates": [
                    [float(value) for value in position]
                    for position in mol.GetConformer().GetPositions()
                ],
            })
    finally:
        writer.close()
    if not mapping or not destination.is_file() or destination.stat().st_size == 0:
        raise RTMScoreError("no selected poses were available for RTMScore")
    return mapping


def _crop_receptor(receptor: Path, mapping: list[dict], destination: Path, cutoff: float) -> dict:
    import numpy as np

    ligand_coordinates = np.asarray(
        [point for record in mapping for point in record["coordinates"]], dtype=float
    )
    if ligand_coordinates.ndim != 2 or ligand_coordinates.shape[1] != 3:
        raise RTMScoreError("invalid ligand coordinates for pocket cropping")

    residues: dict[tuple[str, str, str, str], list[tuple[str, tuple[float, float, float]]]] = {}
    order: list[tuple[str, str, str, str]] = []
    for line in receptor.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            coordinate = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        key = (line[21:22], line[22:26], line[26:27], line[17:20])
        if key not in residues:
            residues[key] = []
            order.append(key)
        residues[key].append((line, coordinate))

    cutoff_sq = float(cutoff) ** 2
    selected: list[tuple[str, str, str, str]] = []
    for key in order:
        residue_coordinates = np.asarray([item[1] for item in residues[key]], dtype=float)
        keep = False
        for start in range(0, len(ligand_coordinates), 2048):
            ligand_chunk = ligand_coordinates[start:start + 2048]
            delta = residue_coordinates[:, None, :] - ligand_chunk[None, :, :]
            if bool(np.any(np.sum(delta * delta, axis=2) <= cutoff_sq)):
                keep = True
                break
        if keep:
            selected.append(key)
    if not selected:
        raise RTMScoreError("pocket crop selected no protein residues")

    lines = [line for key in selected for line, _coordinate in residues[key]]
    destination.write_text("\n".join([*lines, "TER", "END", ""]), encoding="utf-8")
    return {
        "cutoff_angstrom": float(cutoff),
        "residue_count": len(selected),
        "atom_count": len(lines),
        "residues": [
            {"chain": key[0].strip(), "resid": key[1].strip(), "icode": key[2].strip(), "resname": key[3].strip()}
            for key in selected
        ],
    }


def _run_rtmscore_rescoring(
    receptor_pdb: str | Path,
    selected_sdf: str | Path,
    output_dir: str | Path,
    runtime: Mapping[str, Any],
    *,
    model_ensemble: list[int] | tuple[int, ...] = (1, 2, 3),
    pocket_cutoff_angstrom: float = 10.0,
    batch_size: int = 128,
    num_workers: int = 0,
    parallel_graph_build: bool = False,
    device: str = "auto",
    timeout_seconds: int = 86400,
) -> dict:
    receptor = Path(receptor_pdb).resolve()
    poses = Path(selected_sdf).resolve()
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "rtmscore_manifest.json"
    if not receptor.is_file() or not poses.is_file():
        raise RTMScoreError("RTMScore receptor or selected-pose input is missing")
    if not WORKER.is_file():
        raise RTMScoreError(f"RTMScore worker is missing: {WORKER}")

    entrypoint = Path(str(runtime.get("entrypoint") or "")).resolve()
    checkpoints = [Path(str(path)).resolve() for path in runtime.get("checkpoints") or []]
    checksums = [str(value) for value in runtime.get("asset_sha256") or []]
    if not entrypoint.is_file() or len(checkpoints) != 3 or len(checksums) != 3:
        raise RTMScoreError("registry must provide the RTMScore entrypoint and three models")
    for checkpoint, expected in zip(checkpoints, checksums):
        if not checkpoint.is_file() or _sha256(checkpoint) != expected:
            raise RTMScoreError(f"RTMScore checkpoint checksum mismatch: {checkpoint}")

    requested_models = [int(value) for value in model_ensemble]
    if not requested_models or len(set(requested_models)) != len(requested_models):
        raise RTMScoreError("model_ensemble must contain unique model ids")
    if any(value not in {1, 2, 3} for value in requested_models):
        raise RTMScoreError("RTMScore model ids must be 1, 2, or 3")

    input_sdf = destination / "rtmscore_input.sdf"
    pocket_pdb = destination / "rtmscore_pocket.pdb"
    mapping = _load_and_relabel_poses(poses, input_sdf)
    pocket = _crop_receptor(
        receptor, mapping, pocket_pdb, float(pocket_cutoff_angstrom)
    )
    for record in mapping:
        record.pop("coordinates", None)

    expected_ids = {record["upstream_id"] for record in mapping}
    model_results: list[dict] = []
    scores_by_id: dict[str, list[float]] = {identifier: [] for identifier in expected_ids}
    progress_path = destination / "model_progress.json"
    _atomic_json(progress_path, {
        "schema_version": SCHEMA_VERSION,
        "status": "running",
        "mapping": mapping,
        "pocket": pocket,
        "models": model_results,
    })
    for model_id in requested_models:
        request_path = destination / f"model{model_id}_request.json"
        response_path = destination / f"model{model_id}_response.json"
        request = {
            "entrypoint": str(entrypoint),
            "pocket_pdb": str(pocket_pdb),
            "ligand_sdf": str(input_sdf),
            "checkpoint": str(checkpoints[model_id - 1]),
            "pocket_cutoff_angstrom": float(pocket_cutoff_angstrom),
            "batch_size": int(batch_size),
            "num_workers": int(num_workers),
            "parallel_graph_build": bool(parallel_graph_build),
            "device": device,
        }
        _atomic_json(request_path, request)
        try:
            execution = run_tool(
                dict(runtime),
                [WORKER, "--input-json", request_path, "--output-json", response_path],
                command=str(runtime.get("command") or "python"),
                cwd=entrypoint.parent,
                timeout=int(timeout_seconds),
            )
        except ToolRuntimeError as exc:
            raise RTMScoreError(f"RTMScore model {model_id} execution failed: {exc}") from exc
        if not response_path.is_file():
            raise RTMScoreError(f"RTMScore model {model_id} wrote no response")
        response = json.loads(response_path.read_text(encoding="utf-8"))
        if response.get("success") is not True:
            raise RTMScoreError(
                f"RTMScore model {model_id} returned success=false"
            )
        ids = [str(value) for value in response.get("ids") or []]
        values = [float(value) for value in response.get("scores") or []]
        if len(ids) != len(values):
            raise RTMScoreError(
                f"RTMScore model {model_id} returned {len(ids)} ids but "
                f"{len(values)} scores"
            )
        if set(ids) != expected_ids or len(ids) != len(expected_ids):
            raise RTMScoreError(
                f"RTMScore model {model_id} output does not cover the exact pose cohort"
            )
        for identifier, value in zip(ids, values):
            if not math.isfinite(value):
                raise RTMScoreError(f"RTMScore model {model_id} emitted a non-finite score")
            scores_by_id[identifier].append(value)
        model_results.append({
            "model_id": model_id,
            "checkpoint": str(checkpoints[model_id - 1]),
            "checkpoint_sha256": checksums[model_id - 1],
            "device": response.get("device"),
            "torch_version": response.get("torch_version"),
            "runtime_versions": response.get("runtime_versions"),
            "execution": {
                "command": execution.get("command"),
                "environment": execution.get("environment"),
                "returncode": execution.get("returncode"),
                "elapsed_seconds": execution.get("elapsed_seconds"),
                "stdout_tail": str(execution.get("stdout") or "")[-2000:],
                "stderr_tail": str(execution.get("stderr") or "")[-2000:],
            },
            "response_file": str(response_path),
        })
        _atomic_json(progress_path, {
            "schema_version": SCHEMA_VERSION,
            "status": "running",
            "mapping": mapping,
            "pocket": pocket,
            "models": model_results,
        })

    protocol_payload = {
        "method": "RTMScore",
        "receptor_sha256": _sha256(receptor),
        "pose_cohort_sha256": _sha256(poses),
        "pocket_sha256": _sha256(pocket_pdb),
        "pocket_cutoff_angstrom": float(pocket_cutoff_angstrom),
        "model_sha256": [checksums[index - 1] for index in requested_models],
        "entrypoint_sha256": _sha256(entrypoint),
        "worker_sha256": _sha256(WORKER),
        "batch_size": int(batch_size),
        "num_workers": int(num_workers),
        "parallel_graph_build": bool(parallel_graph_build),
        "requested_device": str(device),
        "resolved_devices": [entry.get("device") for entry in model_results],
        "runtime_versions": [
            entry.get("runtime_versions") for entry in model_results
        ],
    }
    protocol_id = "rtmscore:" + hashlib.sha256(
        json.dumps(protocol_payload, sort_keys=True).encode("utf-8")
    ).hexdigest()
    by_upstream = {record["upstream_id"]: record for record in mapping}
    pose_scores = []
    records = []
    for upstream_id in sorted(expected_ids, key=lambda item: by_upstream[item]["record_index"]):
        record = by_upstream[upstream_id]
        values = scores_by_id[upstream_id]
        mean_score = fmean(values)
        score = {
            "schema_version": 1,
            "metric_id": "rtmscore.ensemble_score",
            "value": mean_score,
            "unit": "RTMScore arbitrary unit",
            "direction": "higher",
            "method": "rtmscore",
            "method_version": "three-model-local" if len(values) == 3 else "subset-ensemble-local",
            "stage": "pose_evaluation_ml_rescore",
            "applicability": "applicable",
            "reliability": "medium",
            "pose_id": record["pose_id"],
            "pose_index": record["record_index"],
            "protocol_id": protocol_id,
            # The exact entrypoint, bridge, models, inputs, parameters, device,
            # and package versions are bound above.  RTMScore imports further
            # source-tree modules, so do not overclaim a fully captured runtime.
            "protocol_identity_complete": False,
            "protocol_identity_note": (
                "Pinned local assets and runtime versions are recorded, but the "
                "complete transitive source/environment closure is not hashed."
            ),
            "uncertainty": {
                "kind": "ensemble_standard_deviation",
                "value": pstdev(values) if len(values) > 1 else 0.0,
                "n_models": len(values),
            },
        }
        pose_scores.append(score)
        records.append({
            **record,
            "model_scores": values,
            "ensemble_mean": mean_score,
            "ensemble_stddev": score["uncertainty"]["value"],
        })

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "success": True,
        "status": "passed",
        "result_semantics": (
            "Machine-learning pose rescoring; higher is better within this exact "
            "receptor/pocket/model cohort. It is not an experimental affinity."
        ),
        "protocol": {**protocol_payload, "protocol_id": protocol_id},
        "pocket": pocket,
        "counts": {"requested": len(mapping), "scored": len(pose_scores)},
        "models": model_results,
        "records": records,
        "pose_scores": pose_scores,
        "files": {
            "input_sdf": str(input_sdf),
            "pocket_pdb": str(pocket_pdb),
            "model_progress": str(progress_path),
            "manifest": str(manifest_path),
        },
        "manifest_file": str(manifest_path),
    }
    _atomic_json(manifest_path, manifest)
    return manifest


def run_rtmscore_rescoring(
    receptor_pdb: str | Path,
    selected_sdf: str | Path,
    output_dir: str | Path,
    runtime: Mapping[str, Any],
    *,
    model_ensemble: list[int] | tuple[int, ...] = (1, 2, 3),
    pocket_cutoff_angstrom: float = 10.0,
    batch_size: int = 128,
    num_workers: int = 0,
    parallel_graph_build: bool = False,
    device: str = "auto",
    timeout_seconds: int = 86400,
) -> dict:
    """Run RTMScore and always publish a success or failure manifest."""
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "rtmscore_manifest.json"
    try:
        return _run_rtmscore_rescoring(
            receptor_pdb,
            selected_sdf,
            destination,
            runtime,
            model_ensemble=model_ensemble,
            pocket_cutoff_angstrom=pocket_cutoff_angstrom,
            batch_size=batch_size,
            num_workers=num_workers,
            parallel_graph_build=parallel_graph_build,
            device=device,
            timeout_seconds=timeout_seconds,
        )
    except Exception as exc:
        progress_path = destination / "model_progress.json"
        progress: dict[str, Any] = {}
        if progress_path.is_file():
            try:
                progress = json.loads(progress_path.read_text(encoding="utf-8"))
            except (OSError, ValueError, TypeError):
                progress = {}

        input_files = {
            "receptor_pdb": Path(receptor_pdb).resolve(),
            "selected_sdf": Path(selected_sdf).resolve(),
        }
        input_sha256 = {
            name: _sha256(path)
            for name, path in input_files.items()
            if path.is_file()
        }
        failure = {
            "schema_version": SCHEMA_VERSION,
            "success": False,
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "result_semantics": (
                "RTMScore did not produce a complete exact-cohort ensemble; "
                "no partial score is rankable."
            ),
            "counts": {
                "requested": len(progress.get("mapping") or []),
                "scored": 0,
            },
            "models": list(progress.get("models") or []),
            "pocket": progress.get("pocket"),
            "records": [],
            "pose_scores": [],
            "input_sha256": input_sha256,
            "files": {
                **{name: str(path) for name, path in input_files.items()},
                "model_progress": (
                    str(progress_path) if progress_path.is_file() else None
                ),
                "manifest": str(manifest_path),
            },
            "manifest_file": str(manifest_path),
        }
        _atomic_json(manifest_path, failure)
        return failure


__all__ = ["RTMScoreError", "SCHEMA_VERSION", "run_rtmscore_rescoring"]

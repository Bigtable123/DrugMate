#!/usr/bin/env python3
"""Strict DeepPocket worker executed inside the dedicated Conda environment.

The upstream ``predict.py`` script produces several loosely related files and
does not validate the return status of its internal fpocket command.  This
bridge runs the real upstream workflow, parses the ranked classifier output and
segmented-residue PDB files, and emits one normalized JSON contract.  A zero
exit status alone is never treated as a successful prediction.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from deployment_config import deployment_path


EXPECTED_REPOSITORY_COMMIT = "04ba9f564a81dbf35c6afbb1db5fb9ca255cedc0"
EXPECTED_CHECKPOINT_SHA256 = {
    "first_model_fold1_best_test_auc_85001.pth.tar": (
        "d166c8a36b2475297193dbe31e5a88ff9b853e343281ef21df537029ba053f9a"
    ),
    "seg0_best_test_IOU_91.pth.tar": (
        "aa6d29beda70e79b39504c3bdaf5b42380e9c1e5bc6b0557cb81e9075243faee"
    ),
}


class DeepPocketBridgeError(RuntimeError):
    """An installation, input, execution, or output invariant was violated."""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
        raise DeepPocketBridgeError(
            "repository_unverifiable", f"could not inspect DeepPocket repository: {exc}"
        ) from exc
    commit = (completed.stdout or "").strip()
    if completed.returncode != 0 or len(commit) != 40:
        detail = (completed.stderr or completed.stdout or "git rev-parse failed").strip()
        raise DeepPocketBridgeError("repository_unverifiable", detail)
    if commit != EXPECTED_REPOSITORY_COMMIT:
        raise DeepPocketBridgeError(
            "repository_commit_mismatch",
            f"expected DeepPocket commit {EXPECTED_REPOSITORY_COMMIT}, found {commit}",
        )
    tracked_inputs = [
        "predict.py", "clean_pdb.py", "get_centers.py", "model.py",
        "rank_pockets.py", "segment_pockets.py", "types_and_gninatyper.py",
        "unet.py", "gninamap",
    ]
    dirty = subprocess.run(
        [
            "git", "-C", str(source_dir), "status", "--porcelain",
            "--untracked-files=no", "--", *tracked_inputs,
        ],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    if dirty.returncode != 0:
        detail = (dirty.stderr or dirty.stdout or "git status failed").strip()
        raise DeepPocketBridgeError("repository_unverifiable", detail)
    if dirty.stdout.strip():
        raise DeepPocketBridgeError(
            "repository_dirty",
            "tracked DeepPocket inference source differs from the pinned commit",
        )
    return commit


def _validate_installation(
    source_dir: Path,
    class_checkpoint: Path,
    segmentation_checkpoint: Path,
) -> dict[str, Any]:
    required_source = [
        source_dir / "predict.py",
        source_dir / "model.py",
        source_dir / "unet.py",
        source_dir / "rank_pockets.py",
        source_dir / "segment_pockets.py",
        source_dir / "types_and_gninatyper.py",
        source_dir / "gninamap",
    ]
    required = [*required_source, class_checkpoint, segmentation_checkpoint]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise DeepPocketBridgeError(
            "installation_incomplete",
            "missing DeepPocket resource(s): " + ", ".join(missing),
        )
    if shutil.which("fpocket") is None:
        raise DeepPocketBridgeError(
            "fpocket_missing",
            "DeepPocket requires the fpocket executable in its isolated environment",
        )

    observed_hashes = {
        class_checkpoint.name: _sha256(class_checkpoint),
        segmentation_checkpoint.name: _sha256(segmentation_checkpoint),
    }
    if observed_hashes != EXPECTED_CHECKPOINT_SHA256:
        raise DeepPocketBridgeError(
            "checkpoint_hash_mismatch",
            "DeepPocket checkpoints do not match the validated installation "
            f"(observed={observed_hashes})",
        )
    return {
        "name": "DeepPocket",
        "repository": str(source_dir),
        "repository_commit": _repository_commit(source_dir),
        "class_checkpoint": str(class_checkpoint),
        "segmentation_checkpoint": str(segmentation_checkpoint),
        "checkpoint_sha256": observed_hashes,
    }


def _validate_input(protein_file: Path) -> dict[str, Any]:
    if not protein_file.is_file():
        raise DeepPocketBridgeError(
            "input_missing", f"prepared receptor does not exist: {protein_file}"
        )
    if protein_file.suffix.lower() != ".pdb":
        raise DeepPocketBridgeError(
            "unsupported_input_format", "DeepPocket requires a prepared receptor PDB"
        )
    atom_count = 0
    with protein_file.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_count += 1
    if atom_count == 0:
        raise DeepPocketBridgeError(
            "input_has_no_atoms", f"prepared receptor contains no coordinate records: {protein_file}"
        )
    return {
        "protein_file": str(protein_file),
        "protein_sha256": _sha256(protein_file),
        "coordinate_record_count": atom_count,
    }


def _parse_confidences(path: Path) -> list[float]:
    if not path.is_file():
        raise DeepPocketBridgeError(
            "confidence_output_missing", f"DeepPocket confidence output is missing: {path}"
        )
    try:
        raw = ast.literal_eval(path.read_text(encoding="utf-8").strip())
    except (OSError, SyntaxError, ValueError) as exc:
        raise DeepPocketBridgeError(
            "confidence_output_invalid", f"could not parse {path.name}: {exc}"
        ) from exc
    if not isinstance(raw, (list, tuple)) or not raw:
        raise DeepPocketBridgeError(
            "confidence_output_empty", "DeepPocket returned no classifier probabilities"
        )
    values = []
    for index, item in enumerate(raw, start=1):
        if isinstance(item, bool):
            raise DeepPocketBridgeError(
                "confidence_output_invalid", f"probability {index} is boolean"
            )
        try:
            value = float(item)
        except (TypeError, ValueError) as exc:
            raise DeepPocketBridgeError(
                "confidence_output_invalid", f"probability {index} is not numeric"
            ) from exc
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise DeepPocketBridgeError(
                "confidence_output_invalid",
                f"probability {index} must be finite and within [0, 1], found {value}",
            )
        values.append(value)
    if any(values[index] > values[index - 1] + 1e-8 for index in range(1, len(values))):
        raise DeepPocketBridgeError(
            "confidence_order_invalid", "DeepPocket probabilities are not rank-sorted"
        )
    return values


def _parse_ranked_types(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise DeepPocketBridgeError(
            "ranked_output_missing", f"DeepPocket ranked types output is missing: {path}"
        )
    records = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line:
            continue
        fields = line.split(maxsplit=4)
        if len(fields) != 5:
            raise DeepPocketBridgeError(
                "ranked_output_invalid",
                f"{path.name} line {line_number} has {len(fields)} fields, expected 5",
            )
        try:
            source_rank = int(fields[0])
            center = [float(value) for value in fields[1:4]]
        except ValueError as exc:
            raise DeepPocketBridgeError(
                "ranked_output_invalid", f"invalid rank/center at line {line_number}: {exc}"
            ) from exc
        if source_rank < 1 or not all(math.isfinite(value) for value in center):
            raise DeepPocketBridgeError(
                "ranked_output_invalid", f"invalid rank/center at line {line_number}"
            )
        records.append({
            "original_fpocket_rank": source_rank,
            "center": center,
            "receptor_feature_file": fields[4],
        })
    if not records:
        raise DeepPocketBridgeError(
            "ranked_output_empty", "DeepPocket returned no ranked fpocket candidates"
        )
    source_ranks = [record["original_fpocket_rank"] for record in records]
    if len(source_ranks) != len(set(source_ranks)):
        raise DeepPocketBridgeError(
            "ranked_output_invalid", "DeepPocket ranked output repeats an fpocket candidate"
        )
    return records


def _parse_fpocket_info(path: Path) -> dict[int, dict[str, float]]:
    """Read only upstream fields with stable labels; absence remains explicit."""
    if not path.is_file():
        return {}
    pockets: dict[int, dict[str, float]] = {}
    current: dict[str, float] | None = None
    current_rank: int | None = None
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        match = re.match(r"Pocket\s+(\d+)\s*:", line, flags=re.IGNORECASE)
        if match:
            current_rank = int(match.group(1))
            current = pockets.setdefault(current_rank, {})
            continue
        if current is None:
            continue
        field_match = re.match(
            r"(Score|Druggability Score|Volume)\s*:\s*([-+0-9.eE]+)",
            line,
            flags=re.IGNORECASE,
        )
        if not field_match:
            continue
        key = {
            "score": "fpocket_score",
            "druggability score": "fpocket_druggability_score",
            "volume": "volume",
        }[field_match.group(1).lower()]
        try:
            value = float(field_match.group(2))
        except ValueError:
            continue
        if math.isfinite(value):
            current[key] = value
    return pockets


def _parse_residue_records(path: Path) -> list[dict[str, Any]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    records: list[dict[str, Any]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip()
            icode = line[26:27].strip()
            if not resname:
                continue
            identity = (chain, resnum, icode, resname)
            if identity in seen:
                continue
            seen.add(identity)
            legacy = f"{resname}{resnum}{icode}"
            records.append({
                "resname": resname,
                "resnum": resnum,
                "chain": chain,
                "icode": icode,
                "label": f"{chain}:{legacy}" if chain else legacy,
                "legacy_label": legacy,
            })
    return records


def _copy_raw_run(run_dir: Path, raw_output_dir: Path) -> None:
    raw_output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = raw_output_dir.with_name(f".{raw_output_dir.name}.tmp-{os.getpid()}")
    if staging.exists():
        shutil.rmtree(staging)
    shutil.copytree(run_dir, staging)
    if raw_output_dir.exists():
        shutil.rmtree(raw_output_dir)
    os.replace(staging, raw_output_dir)


def _relative_output(path: Path, run_dir: Path, raw_output_dir: Path) -> str:
    return str((raw_output_dir / path.relative_to(run_dir)).resolve())


def _tail(text: str, limit: int = 4000) -> str:
    return " ".join((text or "").split())[-limit:]


def _torch_load(torch_module, path: Path):
    """Load trusted, hash-pinned local checkpoints across PyTorch versions."""
    try:
        return torch_module.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch_module.load(path, map_location="cpu")


def run_probe(
    source_dir: Path,
    class_checkpoint: Path,
    segmentation_checkpoint: Path,
) -> dict[str, Any]:
    """Validate imports, checkpoint/model compatibility, and required CUDA."""
    producer = _validate_installation(
        source_dir, class_checkpoint, segmentation_checkpoint
    )
    source_path = str(source_dir)
    if source_path not in sys.path:
        sys.path.insert(0, source_path)
    try:
        import torch
        import torch.nn as nn
        import molgrid
        import Bio
        import prody
        from model import Model
        from unet import Unet
    except Exception as exc:
        raise DeepPocketBridgeError(
            "runtime_import_failed", f"DeepPocket dependency import failed: {exc}"
        ) from exc
    if not torch.cuda.is_available():
        raise DeepPocketBridgeError(
            "cuda_unavailable",
            "the installed DeepPocket predict.py requires a functional CUDA runtime",
        )
    try:
        class_payload = _torch_load(torch, class_checkpoint)
        segmentation_payload = _torch_load(torch, segmentation_checkpoint)
        class_state = class_payload["model_state_dict"]
        segmentation_state = segmentation_payload["model_state_dict"]
        classifier = Model()
        classifier.load_state_dict(class_state, strict=True)
        segmenter = nn.DataParallel(Unet(1, None))
        segmenter.load_state_dict(segmentation_state, strict=True)
        # Exercise allocation on the same backend the upstream code hardcodes.
        cuda_probe = torch.zeros((1,), dtype=torch.float32, device="cuda")
        cuda_probe.add_(1.0)
        torch.cuda.synchronize()
        cuda_device = torch.cuda.get_device_name(torch.cuda.current_device())
        del cuda_probe, classifier, segmenter, class_payload, segmentation_payload
        torch.cuda.empty_cache()
    except Exception as exc:
        raise DeepPocketBridgeError(
            "checkpoint_model_probe_failed",
            f"DeepPocket checkpoint/model/CUDA probe failed: {type(exc).__name__}: {exc}",
        ) from exc
    return {
        "schema_version": 1,
        "success": True,
        "mode": "probe",
        "producer": producer,
        "validation": {
            "dependency_imports": True,
            "classifier_state_dict_strict_match": True,
            "segmenter_state_dict_strict_match": True,
            "cuda_available": True,
            "cuda_allocation_succeeded": True,
        },
        "runtime": {
            "python": sys.version.split()[0],
            "torch": str(torch.__version__),
            "molgrid": getattr(molgrid, "__version__", None),
            "biopython": getattr(Bio, "__version__", None),
            "prody": getattr(prody, "__version__", None),
            "cuda_device": cuda_device,
        },
    }


def run_prediction(args: argparse.Namespace) -> dict[str, Any]:
    source_dir = Path(args.source_dir).resolve()
    class_checkpoint = Path(args.class_checkpoint).resolve()
    segmentation_checkpoint = Path(args.segmentation_checkpoint).resolve()
    protein_file = Path(args.protein_file).resolve()
    raw_output_dir = Path(args.raw_output_dir).resolve()
    producer = _validate_installation(
        source_dir, class_checkpoint, segmentation_checkpoint
    )
    input_record = _validate_input(protein_file)

    with tempfile.TemporaryDirectory(prefix="drugmate-deeppocket-") as temporary:
        run_dir = Path(temporary)
        local_protein = run_dir / "receptor.pdb"
        shutil.copy2(protein_file, local_protein)
        command = [
            sys.executable,
            str(source_dir / "predict.py"),
            "--protein", str(local_protein),
            "--class_checkpoint", str(class_checkpoint),
            "--seg_checkpoint", str(segmentation_checkpoint),
            "--rank", str(args.max_candidates),
            "--threshold", str(args.threshold),
            "--mask_dist", str(args.mask_distance),
        ]
        try:
            completed = subprocess.run(
                command,
                cwd=str(source_dir),
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise DeepPocketBridgeError(
                "upstream_start_failed", f"could not start DeepPocket: {exc}"
            ) from exc

        if completed.returncode != 0:
            _copy_raw_run(run_dir, raw_output_dir)
            detail = _tail(completed.stderr) or _tail(completed.stdout)
            raise DeepPocketBridgeError(
                "upstream_failed",
                f"DeepPocket predict.py exited with code {completed.returncode}"
                + (f": {detail}" if detail else ""),
            )

        nowat_stem = "receptor_nowat"
        pockets_dir = run_dir / f"{nowat_stem}_out" / "pockets"
        ranked_path = pockets_dir / "bary_centers_ranked.types"
        confidence_path = pockets_dir / "bary_centers_confidence.txt"
        info_path = run_dir / f"{nowat_stem}_out" / f"{nowat_stem}_info.txt"
        ranked = _parse_ranked_types(ranked_path)
        confidences = _parse_confidences(confidence_path)
        if len(ranked) != len(confidences):
            raise DeepPocketBridgeError(
                "rank_probability_count_mismatch",
                f"ranked candidates ({len(ranked)}) and probabilities "
                f"({len(confidences)}) differ",
            )

        fpocket_metrics = _parse_fpocket_info(info_path)
        candidates = []
        segmentation_records = []
        requested = min(args.max_candidates, len(ranked))
        for position, (ranked_record, confidence) in enumerate(
            zip(ranked[:requested], confidences[:requested]), start=1
        ):
            mask_path = run_dir / f"{nowat_stem}_{position}.dx"
            pocket_path = run_dir / f"{nowat_stem}_pocket{position}.pdb"
            if not mask_path.is_file() or mask_path.stat().st_size == 0:
                raise DeepPocketBridgeError(
                    "segmentation_mask_missing",
                    f"segmentation mask for DeepPocket rank {position} is missing",
                )
            residue_records = _parse_residue_records(pocket_path)
            source_rank = ranked_record["original_fpocket_rank"]
            support_path = pockets_dir / f"pocket{source_rank}_vert.pqr"
            if not support_path.is_file() or support_path.stat().st_size == 0:
                raise DeepPocketBridgeError(
                    "proposal_cloud_missing",
                    f"native proposal support for DeepPocket rank {position} is missing",
                )
            segmentation_records.append({
                "rank": position,
                "original_fpocket_rank": source_rank,
                "classifier_probability": confidence,
                "mask_file": _relative_output(mask_path, run_dir, raw_output_dir),
                "residue_file": (
                    _relative_output(pocket_path, run_dir, raw_output_dir)
                    if pocket_path.is_file() else None
                ),
                "residue_count": len(residue_records),
                "status": "geometry_available" if residue_records else "no_predicted_residues",
            })
            # Keep the classifier's native ranks even when CNN segmentation
            # is empty. Memory fusion derives point5 from the proposal cloud;
            # explicit-site consensus separately skips empty segmented sets.
            metrics = fpocket_metrics.get(source_rank, {})
            candidates.append({
                "rank": position,
                "name": f"DeepPocket_rank_{position}",
                "score": confidence,
                "probability": confidence,
                "center": list(ranked_record["center"]),
                "residue_records": residue_records,
                "residues": [record["legacy_label"] for record in residue_records],
                "original_fpocket_rank": source_rank,
                "volume": metrics.get("volume"),
                "support_file": _relative_output(
                    support_path, run_dir, raw_output_dir,
                ),
                "support_sha256": _sha256(support_path),
                "fpocket_seed_metrics": metrics,
                "candidate_dependency": "fpocket_proposal_then_deeppocket_cnn_ranking_and_segmentation",
                "segmentation_threshold": args.threshold,
                "mask_distance_angstrom": args.mask_distance,
                "segmentation_mask_file": _relative_output(
                    mask_path, run_dir, raw_output_dir
                ),
                "segmented_residue_file": (
                    _relative_output(pocket_path, run_dir, raw_output_dir)
                    if pocket_path.is_file() else None
                ),
                "segmentation_residue_count": len(residue_records),
            })

        _copy_raw_run(run_dir, raw_output_dir)
        return {
            "schema_version": 1,
            "success": True,
            "producer": producer,
            "input": input_record,
            "parameters": {
                "max_candidates": args.max_candidates,
                "segmentation_threshold": args.threshold,
                "mask_distance_angstrom": args.mask_distance,
            },
            "method": {
                "candidate_generation": "fpocket",
                "candidate_ranking": "DeepPocket classifier",
                "pocket_geometry": "CNN-segmented residues when available; native proposal clouds retained for point5 fusion",
                "independence_group": "fpocket_family",
                "note": (
                    "DeepPocket and fpocket are not independent detector votes because "
                    "DeepPocket starts from fpocket proposals."
                ),
            },
            "num_ranked_candidates": len(ranked),
            "num_segmented_candidates_requested": requested,
            "num_pockets": len(candidates),
            "top_pocket": candidates[0],
            "all_pockets": candidates,
            "segmentation_records": segmentation_records,
            "raw_output_dir": str(raw_output_dir),
            "validation": {
                "upstream_exit_zero": True,
                "rank_probability_count_matches": True,
                "probabilities_finite_and_in_unit_interval": True,
                "probabilities_rank_sorted": True,
                "segmentation_masks_present": True,
                "proposal_support_files_present": True,
                "at_least_one_segmented_residue_set": any(
                    candidate["residue_records"] for candidate in candidates
                ),
                "normalized_candidate_count": len(candidates),
            },
            "execution": {
                "upstream_command": command,
                "upstream_returncode": completed.returncode,
                "stdout_tail": _tail(completed.stdout),
                "stderr_tail": _tail(completed.stderr),
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
    parser = argparse.ArgumentParser(description="Strict DeepPocket pipeline bridge")
    parser.add_argument("mode", choices=["probe", "predict"])
    parser.add_argument("--source-dir", default=str(deployment_path("${SBDD_DEEPPOCKET_SOURCE_DIR}")))
    parser.add_argument(
        "--class-checkpoint",
        default=(
            str(deployment_path("${SBDD_DEEPPOCKET_SOURCE_DIR}/models/"
            "first_model_fold1_best_test_auc_85001.pth.tar"))
        ),
    )
    parser.add_argument(
        "--segmentation-checkpoint",
        default=str(deployment_path("${SBDD_DEEPPOCKET_SOURCE_DIR}/models/seg0_best_test_IOU_91.pth.tar")),
    )
    parser.add_argument("--protein-file")
    parser.add_argument("--raw-output-dir")
    parser.add_argument("--output-json")
    parser.add_argument("--max-candidates", type=int, default=3)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--mask-distance", type=float, default=3.5)
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.mode == "predict":
        if not args.protein_file or not args.raw_output_dir:
            parser.error("predict requires --protein-file and --raw-output-dir")
        if not 1 <= args.max_candidates <= 10:
            parser.error("--max-candidates must be between 1 and 10")
        if not 0.0 < args.threshold < 1.0:
            parser.error("--threshold must be between 0 and 1")
        if not math.isfinite(args.mask_distance) or args.mask_distance <= 0:
            parser.error("--mask-distance must be a positive finite number")
    try:
        payload = (
            run_probe(
                Path(args.source_dir).resolve(),
                Path(args.class_checkpoint).resolve(),
                Path(args.segmentation_checkpoint).resolve(),
            )
            if args.mode == "probe" else run_prediction(args)
        )
        _write_payload(payload, args.output_json)
        return 0
    except DeepPocketBridgeError as exc:
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

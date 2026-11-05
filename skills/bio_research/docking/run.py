#!/usr/bin/env python3
"""Docking Skill - registry-selected molecular docking engine.

Usage:
    python run.py --input-file explicit_docking_config.json --output-dir ./output/

The input must declare ``search_space_source`` and a complete explicit
primary-engine plus ``redocking_validation`` enable matrix; see SKILL.md.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import shutil
import sys
import uuid
from collections import Counter, defaultdict
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from deployment_config import deployment_path
from gnina_utils import (
    check_gnina_available,
    load_docking_box,
    normalize_docking_box,
    run_gnina_docking,
)
from docking_common import compound_name
from autodock_suite_utils import run_adfr_docking, run_autodock4_docking
from diffdock_utils import run_diffdock_docking
from dock6_utils import run_dock6_docking
from io_utils import tool_config
from constants import pose_quota_for_compound_count
from receptor_preparation import run_docking_prep
from pipeline_registry import DOCKING_PRIMARY_ENGINES, METRIC_RULES_VERSION
from pose_scores import LEGACY_SCORE_SPECS
from pose_identity import (
    LIGAND_LINEAGE_PROPERTIES, matches_declared_molecule, sdf_molblock_sha256,
)
from docking_artifact_contract import (
    REDOCKING_SEALED_FIELDS,
    redocking_measurement_error,
)
from rdock_utils import run_rdock_docking
from unidock_utils import run_unidock_docking
from vina_utils import run_autodock_vina_docking


PRIMARY_ENGINES = tuple(DOCKING_PRIMARY_ENGINES)
SEARCH_SPACE_SOURCES = {
    "reference_ligand",
    "docking_box",
    "whole_receptor",
    "precomputed_grid",
}
SCIENTIFIC_COHORT_COMPONENT_FIELDS = (
    "protein_state_id",
    "pocket_id",
    "round_id",
    "protocol_profile_id",
    "metric_catalog_version",
)
SCIENTIFIC_COHORT_FIELDS = (
    "cohort_id", *SCIENTIFIC_COHORT_COMPONENT_FIELDS,
)
SCIENTIFIC_COHORT_PROVENANCE_FIELDS = (
    "cohort_identity_complete",
    "cohort_identity_missing_fields",
    "cohort_identity_derived_fields",
    "cohort_provenance_sources",
)
SCIENTIFIC_COHORT_OUTPUT_FIELDS = (
    *SCIENTIFIC_COHORT_FIELDS, *SCIENTIFIC_COHORT_PROVENANCE_FIELDS,
)
DOCKING_CANONICAL_ARTIFACTS = (
    "docking_report.json",
    "docking_scores.json",
    "docked_poses.sdf",
    "redocked_reference.sdf",
)


class DockedPoseContractError(ValueError):
    """A docked SDF and its score table cannot be joined safely."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        details: dict[str, object] | None = None,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.details = dict(details or {})


def _invalidate_previous_docking_artifacts(output_dir: Path) -> None:
    """Remove prior flat outputs before a same-directory rerun begins."""
    for name in DOCKING_CANONICAL_ARTIFACTS:
        (output_dir / name).unlink(missing_ok=True)


def _stage_reported_docked_sdf(result: dict, output_dir: Path) -> None:
    """Move an engine's flat SDF away until its score contract is committed."""
    reported = result.get("docked_sdf")
    if reported in (None, ""):
        return
    source = Path(str(reported))
    canonical = output_dir / "docked_poses.sdf"
    if not source.is_file() or source.resolve() != canonical.resolve():
        return
    staging_dir = output_dir / ".docking_uncommitted"
    staging_dir.mkdir(parents=True, exist_ok=True)
    staged = staging_dir / f"{uuid.uuid4().hex}-docked_poses.sdf"
    source.replace(staged)
    result["docked_sdf"] = str(staged)


def _preserve_failed_redocking(result: dict, output_dir: Path) -> None:
    """Keep failed-attempt coordinates without exposing a canonical artifact."""
    source = output_dir / "redocked_reference.sdf"
    if not source.is_file():
        return
    directory = output_dir / ".docking_uncommitted"
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f"{uuid.uuid4().hex}-redocked_reference.sdf"
    source.replace(target)
    result["unvalidated_redocked_sdf"] = str(target)
    validation = result.get("redocking_validation")
    if isinstance(validation, dict):
        validation["redocked_sdf"] = str(target)
        validation["artifact_available"] = False
        validation.pop("redocked_sdf_artifact", None)


def _mark_redocking_measurement_unavailable(
    validation: dict,
    reason: str,
) -> None:
    """Remove an unprovable RMSD verdict while preserving attempt diagnostics."""
    reported_reason = validation.get("reason")
    if reported_reason not in (None, "") and reported_reason != reason:
        validation["reported_reason"] = reported_reason
    validation.update({
        "measured": False,
        "applicable": False,
        "not_applicable": True,
        "passed": None,
        "gate": "not_applicable",
        "best_rmsd_angstrom": None,
        "best_rmsd_angstrom_display": None,
        "best_rmsd_pose_index": None,
        "artifact_available": False,
        "reason": reason,
    })
    validation.pop("redocked_sdf_artifact", None)


def _read_physical_sdf_record_bytes(path: Path) -> list[bytes]:
    """Read exact delimiter-framed SDF records without shifting bad rows."""
    blocks: list[bytes] = []
    current: list[bytes] = []
    saw_delimiter = False
    with path.open("rb") as handle:
        for line in handle:
            if line.rstrip(b"\r\n") == b"$$$$":
                saw_delimiter = True
                blocks.append(b"".join(current))
                current = []
            else:
                current.append(line)
    trailing = b"".join(current)
    if trailing.strip():
        raise DockedPoseContractError(
            "unterminated_redocked_sdf_record",
            "redocked SDF has data after its last $$$$ delimiter",
        )
    if not saw_delimiter or not blocks:
        raise DockedPoseContractError(
            "empty_or_unframed_redocked_sdf",
            "redocked SDF contains no framed molecule records",
        )
    return blocks


def _validated_redocking_record_identity(
    path: Path,
    pose_index: object,
) -> dict[str, object]:
    """Prove that a best-RMSD index names a valid physical SDF record."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise DockedPoseContractError(
            "rdkit_unavailable_for_redocking_validation",
            "RDKit is required to validate the sealed redocking pose",
        ) from exc

    try:
        index = _strict_integer(pose_index)
    except ValueError as exc:
        raise DockedPoseContractError(
            "invalid_redocking_pose_index",
            "best_rmsd_pose_index is not a canonical integer",
        ) from exc
    blocks = _read_physical_sdf_record_bytes(path)
    if index < 0 or index >= len(blocks):
        raise DockedPoseContractError(
            "redocking_pose_index_out_of_range",
            f"best_rmsd_pose_index {index} is outside 0..{len(blocks) - 1}",
        )

    block = blocks[index]
    framed = block if block.endswith((b"\n", b"\r")) else block + b"\n"
    supplier = Chem.ForwardSDMolSupplier(
        io.BytesIO(framed + b"$$$$\n"),
        removeHs=False,
        strictParsing=False,
    )
    molecules = list(supplier)
    molecule = molecules[0] if len(molecules) == 1 else None
    if (
        molecule is None
        or molecule.GetNumAtoms() == 0
        or molecule.GetNumConformers() != 1
    ):
        raise DockedPoseContractError(
            "invalid_redocking_pose_record",
            f"physical redocked SDF record {index} is not a valid pose",
        )
    conformer = molecule.GetConformer()
    for atom_index in range(molecule.GetNumAtoms()):
        position = conformer.GetAtomPosition(atom_index)
        if not all(
            math.isfinite(value)
            for value in (position.x, position.y, position.z)
        ):
            raise DockedPoseContractError(
                "invalid_redocking_pose_geometry",
                "best redocked pose contains non-finite coordinates",
            )
    return {
        "physical_sdf_records": len(blocks),
        "best_rmsd_pose_index": index,
        "best_pose_sdf_record_sha256": hashlib.sha256(block).hexdigest(),
    }


def _seal_redocking_artifact(
    result: dict,
    output_dir: Path,
) -> dict | None:
    """Seal a measured, physically valid current-run redocking pose."""
    validation = result.get("redocking_validation")
    if not isinstance(validation, dict):
        return None
    canonical = Path(output_dir) / "redocked_reference.sdf"
    if (
        validation.get("performed") is not True
        or validation.get("measured") is not True
    ):
        _mark_redocking_measurement_unavailable(
            validation,
            str(validation.get("reason") or "redocking RMSD was not measured"),
        )
        _preserve_failed_redocking(result, output_dir)
        return None
    path_value = validation.get("redocked_sdf")
    if path_value in (None, ""):
        _mark_redocking_measurement_unavailable(
            validation,
            "measured redocking did not declare its SDF artifact",
        )
        _preserve_failed_redocking(result, output_dir)
        return None
    path = Path(str(path_value))
    try:
        same_path = path.resolve() == canonical.resolve()
    except OSError:
        same_path = False
    if not same_path or not path.is_file():
        _mark_redocking_measurement_unavailable(
            validation,
            "measured redocking SDF is missing or is not the current canonical artifact",
        )
        _preserve_failed_redocking(result, output_dir)
        return None
    try:
        rmsd = float(validation.get("best_rmsd_angstrom"))
        if not math.isfinite(rmsd):
            raise ValueError("non-finite RMSD")
        record_identity = _validated_redocking_record_identity(
            path,
            validation.get("best_rmsd_pose_index"),
        )
        identity = _file_identity(path)
    except (DockedPoseContractError, OSError, TypeError, ValueError) as exc:
        _mark_redocking_measurement_unavailable(
            validation,
            f"redocking measurement artifact could not be verified: {exc}",
        )
        _preserve_failed_redocking(result, output_dir)
        return None
    if identity is None:
        _mark_redocking_measurement_unavailable(
            validation,
            "redocking measurement artifact could not be sealed",
        )
        _preserve_failed_redocking(result, output_dir)
        return None
    identity.update(record_identity)
    measurement = {
        field: validation.get(field) for field in REDOCKING_SEALED_FIELDS
    }
    semantic_error = redocking_measurement_error(
        measurement,
        redocked_reference_present=True,
    )
    if semantic_error is not None:
        _mark_redocking_measurement_unavailable(
            validation,
            f"redocking measurement is internally inconsistent: {semantic_error}",
        )
        _preserve_failed_redocking(result, output_dir)
        return None
    validation["redocked_sdf_artifact"] = identity
    validation["artifact_available"] = True
    return identity


def _publish_validated_docking_artifacts(
    result: dict,
    output_dir: Path,
    *,
    redocked_identity: dict | None = None,
    evaluation_receptor_identity: dict | None = None,
) -> tuple[Path, Path]:
    """Publish SDF+scores, then seal their hashes in the report commit marker."""
    source = Path(str(result.get("docked_sdf") or ""))
    if not source.is_file():
        raise DockedPoseContractError(
            "validated_docked_sdf_missing",
            f"validated docked SDF is missing: {source}",
        )
    # Validate every non-output member of the commit before exposing either
    # canonical output path.  Otherwise a late contract error could leave a
    # seemingly consumable SDF/JSON pair from a failed publication attempt.
    protocol_fingerprint = str(
        (result.get("score_protocol") or {}).get("fingerprint_sha256")
        or ""
    )
    try:
        valid_protocol_fingerprint = bool(
            len(protocol_fingerprint) == 64
            and int(protocol_fingerprint, 16) >= 0
        )
    except ValueError:
        valid_protocol_fingerprint = False
    if not valid_protocol_fingerprint:
        raise DockedPoseContractError(
            "docking_score_protocol_unavailable",
            "cannot commit docking artifacts without a score-protocol fingerprint",
        )
    if not isinstance(evaluation_receptor_identity, dict):
        raise DockedPoseContractError(
            "docking_evaluation_receptor_unavailable",
            "cannot commit docking artifacts without the evaluation receptor identity",
        )
    receptor_identity_path = str(
        evaluation_receptor_identity.get("path") or ""
    ).strip()
    try:
        current_receptor_identity = _file_identity(receptor_identity_path)
    except OSError as exc:
        raise DockedPoseContractError(
            "docking_evaluation_receptor_unavailable",
            "cannot revalidate the evaluation receptor before artifact commit",
        ) from exc
    if (
        not isinstance(current_receptor_identity, dict)
        or current_receptor_identity.get("sha256")
        != evaluation_receptor_identity.get("sha256")
        or current_receptor_identity.get("size_bytes")
        != evaluation_receptor_identity.get("size_bytes")
    ):
        raise DockedPoseContractError(
            "docking_evaluation_receptor_changed",
            "the evaluation receptor changed before docking artifact commit",
        )
    evaluation_receptor_identity = current_receptor_identity

    target_sdf = output_dir / "docked_poses.sdf"
    target_scores = output_dir / "docking_scores.json"
    publish_id = uuid.uuid4().hex
    temporary_sdf = output_dir / f".{publish_id}.docked_poses.sdf.publish.tmp"
    temporary_scores = (
        output_dir / f".{publish_id}.docking_scores.json.publish.tmp"
    )
    try:
        shutil.copy2(source, temporary_sdf)
        temporary_scores.write_text(
            json.dumps(
                result.get("scores") or [],
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            ),
            encoding="utf-8",
        )
        temporary_sdf.replace(target_sdf)
        temporary_scores.replace(target_scores)
    except (OSError, TypeError, ValueError) as exc:
        target_sdf.unlink(missing_ok=True)
        target_scores.unlink(missing_ok=True)
        raise DockedPoseContractError(
            "docking_artifact_publish_failed",
            f"could not publish validated docking artifact set: {exc}",
        ) from exc
    finally:
        temporary_sdf.unlink(missing_ok=True)
        temporary_scores.unlink(missing_ok=True)

    sdf_identity = _file_identity(target_sdf)
    score_identity = _file_identity(target_scores)
    contract = result.setdefault("docked_pose_artifact_contract", {})
    redocked_present = isinstance(redocked_identity, dict)
    redocking = result.get("redocking_validation")
    redocking_measurement = {
        field: redocking.get(field)
        if isinstance(redocking, dict) else None
        for field in REDOCKING_SEALED_FIELDS
    }
    seal_payload = {
        "docked_sdf_sha256": sdf_identity["sha256"],
        "docking_scores_sha256": score_identity["sha256"],
        "physical_sdf_records": contract.get("physical_sdf_records"),
        "score_records": contract.get("score_records"),
        "redocked_reference_present": redocked_present,
        "redocked_reference_sdf_sha256": (
            redocked_identity.get("sha256") if redocked_present else None
        ),
        "redocked_reference_sdf_size_bytes": (
            redocked_identity.get("size_bytes") if redocked_present else None
        ),
        "redocked_reference_physical_sdf_records": (
            redocked_identity.get("physical_sdf_records")
            if redocked_present else None
        ),
        "redocked_reference_best_rmsd_pose_index": (
            redocked_identity.get("best_rmsd_pose_index")
            if redocked_present else None
        ),
        "redocked_reference_best_pose_sdf_record_sha256": (
            redocked_identity.get("best_pose_sdf_record_sha256")
            if redocked_present else None
        ),
        # The redocked SDF hash proves which coordinates were measured; this
        # snapshot also commits the reported RMSD, thresholds, and verdict.
        # Without it, editing docking_report.json alone could turn the same
        # sealed pose into an arbitrary pass/fail gate.
        "redocking_measurement": redocking_measurement,
        # Step 4 evaluates receptor/pose geometry, not ligand coordinates in
        # isolation.  Seal the exact prepared receptor and score protocol so a
        # same-ligand artifact cannot be evaluated against another target/run.
        "evaluation_receptor_sha256": evaluation_receptor_identity.get(
            "sha256"
        ),
        "evaluation_receptor_size_bytes": evaluation_receptor_identity.get(
            "size_bytes"
        ),
        "score_protocol_fingerprint_sha256": protocol_fingerprint,
    }
    contract.update({
        "artifact_set_schema_version": 4,
        "status": "passed",
        "committed": True,
        "docked_sdf": sdf_identity,
        "docking_scores_json": score_identity,
        "redocked_reference_present": redocked_present,
        "redocking_measurement": redocking_measurement,
        "evaluation_receptor": dict(evaluation_receptor_identity),
        "score_protocol_fingerprint_sha256": protocol_fingerprint,
        "artifact_set_id": "sha256:" + hashlib.sha256(
            json.dumps(
                seal_payload, sort_keys=True, separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest(),
    })
    if redocked_present:
        contract["redocked_reference_sdf"] = dict(redocked_identity)
    else:
        contract.pop("redocked_reference_sdf", None)
    result["docked_sdf"] = str(target_sdf)
    if source.parent == output_dir / ".docking_uncommitted":
        source.unlink(missing_ok=True)
    return target_sdf, target_scores


def run(input_data: dict, output_dir: Path) -> dict:
    """Run molecular docking pipeline."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Calling this function for a same-directory rerun invalidates the previous
    # commit marker immediately.  Configuration/input validation can return
    # early, so delaying this until engine dispatch would leave a stale success
    # report consumable after the failed attempt.
    _invalidate_previous_docking_artifacts(output_dir)
    try:
        engine = _selected_engine(input_data)
        validate_redocking = _selected_redocking(input_data)
    except ValueError as exc:
        return {"success": False, "error": str(exc)}

    raw_search_space_source = input_data.get("search_space_source")
    if raw_search_space_source is None:
        return {
            "success": False,
            "error": "search_space_source is required",
        }
    search_space_source = str(raw_search_space_source).strip()
    if search_space_source not in SEARCH_SPACE_SOURCES:
        return {
            "success": False,
            "error": (
                "search_space_source must be one of "
                + ", ".join(sorted(SEARCH_SPACE_SOURCES))
            ),
        }
    if engine == "diffdock" and search_space_source != "whole_receptor":
        return {
            "success": False,
            "error": "DiffDock requires search_space_source='whole_receptor'",
        }
    if engine == "dock6" and search_space_source != "precomputed_grid":
        return {
            "success": False,
            "error": "DOCK 6 requires search_space_source='precomputed_grid'",
        }
    if engine not in {"diffdock", "dock6"} and search_space_source not in {
        "reference_ligand", "docking_box",
    }:
        return {
            "success": False,
            "error": (
                f"{engine} requires search_space_source='reference_ligand' "
                "or 'docking_box'"
            ),
        }

    if input_data.get("is_covalent") is True:
        return {
            "success": False,
            "error": (
                "is_covalent=true cannot run the current Step 4 protocol: "
                "protonation preparation is wired, but covalent atom selection, "
                "bond geometry, and covalent docking are not yet implemented"
            ),
        }

    receptor_file = input_data.get("receptor_file")
    receptor_pdbqt_file = input_data.get("receptor_pdbqt_file")
    ligand_file = input_data.get("ligand_file")
    reference_ligand_file = input_data.get("reference_ligand_file")
    ligand_preparation_manifest_file = input_data.get(
        "ligand_preparation_manifest_file"
    )
    dock6_receptor_site_file = input_data.get("dock6_receptor_site_file")
    dock6_grid_nrg_file = input_data.get("dock6_grid_nrg_file")
    dock6_grid_bmp_file = input_data.get("dock6_grid_bmp_file")
    docking_box_file = input_data.get("docking_box_file")
    docking_box_inline = input_data.get("docking_box")

    if not receptor_file or not ligand_file:
        return {"success": False, "error": "receptor_file and ligand_file are required"}

    receptor_path = Path(receptor_file)
    receptor_pdbqt_path = (
        Path(receptor_pdbqt_file) if receptor_pdbqt_file else None
    )
    ligand_path = Path(ligand_file)
    ref_path = Path(reference_ligand_file) if reference_ligand_file else None
    ligand_manifest_path = (
        Path(ligand_preparation_manifest_file)
        if ligand_preparation_manifest_file else None
    )
    dock6_receptor_site_path = (
        Path(dock6_receptor_site_file) if dock6_receptor_site_file else None
    )
    dock6_grid_nrg_path = (
        Path(dock6_grid_nrg_file) if dock6_grid_nrg_file else None
    )
    dock6_grid_bmp_path = (
        Path(dock6_grid_bmp_file) if dock6_grid_bmp_file else None
    )

    required_files = [(receptor_path, "receptor_file"), (ligand_path, "ligand_file")]
    if engine == "dock6":
        dock6_assets = (
            (dock6_receptor_site_path, "dock6_receptor_site_file"),
            (dock6_grid_nrg_path, "dock6_grid_nrg_file"),
            (dock6_grid_bmp_path, "dock6_grid_bmp_file"),
        )
        missing_asset_names = [
            name for path, name in dock6_assets if path is None
        ]
        if missing_asset_names:
            return {
                "success": False,
                "error": (
                    "DOCK 6 requires injected receptor-site/grid artifacts: "
                    + ", ".join(missing_asset_names)
                ),
            }
        required_files.extend(
            (path, name) for path, name in dock6_assets if path is not None
        )
    if ref_path is not None:
        required_files.append((ref_path, "reference_ligand_file"))
    if ligand_manifest_path is not None:
        required_files.append((
            ligand_manifest_path,
            "ligand_preparation_manifest_file",
        ))
    for f, name in required_files:
        if not f.exists():
            return {"success": False, "error": f"{name} not found: {f}"}

    docking_box = None
    engine_reference_path = None
    if search_space_source == "reference_ligand":
        if ref_path is None:
            return {
                "success": False,
                "error": (
                    "search_space_source='reference_ligand' requires "
                    "reference_ligand_file"
                ),
            }
        engine_reference_path = ref_path
    elif search_space_source == "docking_box":
        try:
            if docking_box_file:
                docking_box = load_docking_box(Path(docking_box_file))
            elif docking_box_inline is not None:
                docking_box = normalize_docking_box(docking_box_inline)
            else:
                return {
                    "success": False,
                    "error": (
                        "search_space_source='docking_box' requires "
                        "docking_box_file or docking_box"
                    ),
                }
        except ValueError as exc:
            return {"success": False, "error": str(exc)}

    if engine != "gnina" and validate_redocking:
        return {
            "success": False,
            "error": (
                "redocking_validation currently depends on the GNINA engine; "
                "disable docking.redocking_validation for the selected "
                f"{engine} engine"
            ),
        }
    if validate_redocking and search_space_source != "reference_ligand":
        return {
            "success": False,
            "error": (
                "redocking_validation currently requires "
                "search_space_source='reference_ligand'"
            ),
        }

    preparation = run_docking_prep(input_data, output_dir / "docking_prep")
    if not preparation["success"]:
        report = {"success": False, "error": "docking_prep failed: " + preparation["error"],
                  "docking_prep": preparation, "files": preparation["files"]}
        _write_report(output_dir, report)
        return report
    receptor_path = Path(preparation["files"]["prepared_receptor"])
    receptor_pdbqt_path = (Path(preparation["files"]["receptor_pdbqt"])
                           if preparation["files"].get("receptor_pdbqt") else None)
    ligand_path = Path(preparation["files"]["docking_species_library"])
    input_data = dict(input_data)
    input_data["_pose_quota"] = pose_quota_for_compound_count(
        _count_sdf_records_up_to(ligand_path, 1001)
    )
    input_data["protein_state_id"] = "receptor_sha256:" + hashlib.sha256(receptor_path.read_bytes()).hexdigest()
    input_data["cohort_provenance_sources"] = {
        **(input_data.get("cohort_provenance_sources") or {}),
        "protein_state_id": "docking_prep:prepared_receptor_sha256",
    }

    parameters = {
        "engine": engine,
        "pose_quota_per_compound": input_data["_pose_quota"],
        "receptor_file": str(receptor_path),
        "receptor_pdbqt_file": (
            str(receptor_pdbqt_path)
            if receptor_pdbqt_path is not None else None
        ),
        "ligand_file": str(ligand_path),
        "ligand_preparation_manifest_file": (
            str(ligand_manifest_path)
            if ligand_manifest_path is not None else None
        ),
        "reference_ligand_file": str(ref_path) if ref_path is not None else None,
        "dock6_receptor_site_file": (
            str(dock6_receptor_site_path)
            if dock6_receptor_site_path is not None else None
        ),
        "dock6_grid_nrg_file": (
            str(dock6_grid_nrg_path)
            if dock6_grid_nrg_path is not None else None
        ),
        "dock6_grid_bmp_file": (
            str(dock6_grid_bmp_path)
            if dock6_grid_bmp_path is not None else None
        ),
        "docking_box_file": str(docking_box_file) if docking_box_file else None,
        "docking_box": docking_box,
        "search_space_source": search_space_source,
        "search_space_mode": search_space_source,
        "validate_redocking": validate_redocking,
    }
    parameters.update({
        field: input_data[field]
        for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS
        if input_data.get(field) not in (None, "")
    })

    result, engine_parameters = ENGINE_RUNNERS[engine](
        input_data=input_data,
        receptor_file=receptor_path,
        receptor_pdbqt_file=receptor_pdbqt_path,
        ligand_file=ligand_path,
        reference_ligand_file=engine_reference_path,
        output_dir=output_dir,
        docking_box=docking_box,
        validate_redocking=validate_redocking,
    )
    parameters["engine_params"] = engine_parameters
    try:
        _stage_reported_docked_sdf(result, output_dir)
    except OSError as exc:
        result["success"] = False
        result["error"] = f"could not stage docking output for validation: {exc}"
        result.setdefault("operational_errors", []).append({
            "tool": engine,
            "step": "docked_pose_artifact_contract",
            "code": "docked_pose_artifact_staging_failed",
            "scope": "run",
            "retryable": True,
            "reason": result["error"],
        })

    if result.get("success"):
        try:
            _annotate_score_provenance(
                result,
                engine,
                engine_parameters,
                parameters,
            )
            try:
                _attach_scientific_cohort(result, input_data, parameters)
            except DockedPoseContractError as exc:
                if exc.code not in {
                    "docked_sdf_record_parse_failed", "invalid_docked_pose_geometry",
                    "docked_molecule_identity_changed",
                }:
                    raise
                _isolate_invalid_docked_molecules(result, ligand_path, output_dir, engine)
                _attach_scientific_cohort(result, input_data, parameters)
        except (DockedPoseContractError, OSError) as exc:
            code = getattr(exc, "code", "docked_pose_artifact_write_failed")
            result["success"] = False
            result["error"] = str(exc)
            operational_error = {
                "tool": engine,
                "step": "docked_pose_artifact_contract",
                "code": code,
                "scope": "run",
                "retryable": True,
                "reason": str(exc),
            }
            details = getattr(exc, "details", None)
            if details:
                operational_error["details"] = details
            result.setdefault("operational_errors", []).append(
                operational_error
            )

    # The ledger remains useful diagnostics on a failed run. Consumable score
    # and SDF paths are published below only after the artifact contract passes.
    _publish_molecule_ledger(
        result,
        ligand_file=ligand_path,
        output_dir=output_dir,
        engine=engine,
    )
    if result.get("success") and result.get(
        "molecule_ledger_unmatched_score_records"
    ):
        unmatched = result["molecule_ledger_unmatched_score_records"]
        result["success"] = False
        result["error"] = (
            f"{len(unmatched)} docking score record(s) cannot be linked to "
            "an input molecule"
        )
        result.setdefault("operational_errors", []).append({
            "tool": engine,
            "step": "docked_pose_artifact_contract",
            "code": "unmatched_docking_score_identity",
            "scope": "run",
            "retryable": True,
            "reason": result["error"],
        })

    published_sdf: Path | None = None
    published_scores: Path | None = None
    redocked_identity: dict | None = None
    if result.get("success"):
        redocked_identity = _seal_redocking_artifact(result, output_dir)
        try:
            published_sdf, published_scores = (
                _publish_validated_docking_artifacts(
                    result,
                    output_dir,
                    redocked_identity=redocked_identity,
                    evaluation_receptor_identity=_file_identity(
                        receptor_path
                    ),
                )
            )
        except DockedPoseContractError as exc:
            result["success"] = False
            result["error"] = str(exc)
            result.setdefault("operational_errors", []).append({
                "tool": engine,
                "step": "docked_pose_artifact_publish",
                "code": exc.code,
                "scope": "run",
                "retryable": True,
                "reason": str(exc),
            })

    if not result.get("success"):
        _preserve_failed_redocking(result, output_dir)
        unvalidated_sdf = result.pop("docked_sdf", None)
        unvalidated_scores = result.pop("scores", None)
        result["unvalidated_artifacts"] = {
            "docked_sdf": unvalidated_sdf,
            "score_record_count": (
                len(unvalidated_scores)
                if isinstance(unvalidated_scores, list) else 0
            ),
        }

    report = {**result, "engine": engine, "parameters": parameters}
    ledger_files = {
        key: value
        for key, value in (
            (
                "docking_molecule_ledger_json",
                result.get("molecule_ledger_json"),
            ),
            (
                "docking_molecule_ledger_tsv",
                result.get("molecule_ledger_tsv"),
            ),
        )
        if value
    }
    if ledger_files:
        report["files"] = ledger_files
    for key in ("docking_recovery_manifest", "raw_docked_sdf", "unvalidated_redocked_sdf"):
        if result.get(key):
            report.setdefault("files", {})[key] = result[key]

    if result["success"]:
        report_file = output_dir / "docking_report.json"
        report.setdefault("files", {}).update({
            "docked_poses_sdf": str(published_sdf),
            "docking_scores_json": str(published_scores),
            "docking_report_json": str(report_file),
            "engine_log": result.get("log_file", ""),
            "gpu_diagnostics_json": result.get("gpu_diagnostics_file", ""),
            "gpu_usage_timeseries_json": result.get("gpu_usage_file", ""),
        })
        if redocked_identity is not None:
            report["files"]["redocked_reference_sdf"] = redocked_identity[
                "path"
            ]
        if engine == "gnina":
            report["files"]["gnina_log"] = result.get("log_file", "")
        elif engine == "unidock":
            report["files"]["unidock_log"] = result.get("log_file", "")
            report["files"]["unidock_native_dir"] = result.get(
                "native_output_dir", "",
            )
        elif engine == "autodock_vina":
            report["files"]["autodock_vina_log"] = result.get(
                "log_file", ""
            )
            report["files"]["vina_native_dir"] = result.get(
                "native_output_dir", ""
            )
        elif engine == "diffdock":
            report["files"]["diffdock_log"] = result.get("log_file", "")
            report["files"]["diffdock_native_dir"] = result.get(
                "native_output_dir", ""
            )
        elif engine == "dock6":
            report["files"]["dock6_log"] = result.get("log_file", "")
            report["files"]["dock6_native_dir"] = result.get(
                "native_output_dir", ""
            )
        elif engine == "rdock":
            report["files"]["rdock_log"] = result.get("log_file", "")
            report["files"]["rdock_native_dir"] = result.get(
                "native_output_dir", ""
            )
        elif engine == "autodock4":
            report["files"]["autodock4_log"] = result.get("log_file", "")
            report["files"]["autodock4_native_dir"] = result.get(
                "native_output_dir", ""
            )
        elif engine == "adfr":
            report["files"]["adfr_log"] = result.get("log_file", "")
            report["files"]["adfr_native_dir"] = result.get(
                "native_output_dir", ""
            )

    report["docking_prep"] = preparation
    report.setdefault("files", {}).update(preparation["files"])
    _write_report(output_dir, report)
    return report


def _input_molecule_rows(ligand_file: Path) -> list[dict]:
    """Read every SDF record into a stable, engine-neutral input ledger."""
    try:
        from rdkit import Chem
    except ImportError:
        return _input_molecule_rows_text(ligand_file)

    rows: list[dict] = []
    for source_index, block in enumerate(_read_sdf_record_blocks(ligand_file)):
        framed = block if block.endswith(("\n", "\r")) else block + "\n"
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO((framed + "$$$$\n").encode("utf-8")),
            removeHs=False,
            strictParsing=False,
        )
        molecules = list(supplier)
        molecule = molecules[0] if len(molecules) == 1 else None
        if molecule is None:
            title = block.splitlines()[0].strip() if block.splitlines() else ""
            rows.append({
                "source_molecule_index": source_index,
                "compound_id": title or f"unparseable_record_{source_index:06d}",
                "mol_name": title or f"unparseable_record_{source_index:06d}",
                "chemical_identity_id": None,
                "chemical_id": None,
                "canonical_smiles": None,
                "input_status": "invalid",
                "_input_aliases": [],
            })
            continue

        aliases = []
        for prop in (
            "compound_id", "CompoundID", "compound_name", "_Name", "ID", "id",
        ):
            if molecule.HasProp(prop):
                value = molecule.GetProp(prop).strip()
                if value and value not in aliases:
                    aliases.append(value)
        identifier = compound_name(molecule, source_index)
        if identifier not in aliases:
            aliases.append(identifier)
        try:
            smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True)
        except Exception:
            smiles = None
        row = {
            "source_molecule_index": source_index,
            "compound_id": identifier,
            "chemical_identity_id": _molecule_property(
                molecule, "chemical_identity_id"
            ),
            "chemical_id": _molecule_property(molecule, "chemical_id"),
            "canonical_smiles": smiles or None,
            "input_status": "parseable",
            "_input_aliases": aliases,
        }
        row.update({
            field: molecule.GetProp(field).strip()
            for field in LIGAND_LINEAGE_PROPERTIES
            if molecule.HasProp(field) and molecule.GetProp(field).strip()
        })
        row["mol_name"] = _molecule_property(molecule, "_Name") or identifier
        rows.append(row)
    return rows


def _input_molecule_rows_text(ligand_file: Path) -> list[dict]:
    """Best-effort identity ledger for the GNINA text-parser fallback."""
    rows: list[dict] = []
    for record_index, block in enumerate(_read_sdf_record_blocks(ligand_file)):
        lines = block.splitlines()
        properties: dict[str, str] = {}
        for line_index, line in enumerate(lines):
            if not line.startswith(">") or "<" not in line or ">" not in line:
                continue
            start = line.find("<") + 1
            end = line.find(">", start)
            if end <= start:
                continue
            value = lines[line_index + 1].strip() if line_index + 1 < len(lines) else ""
            if value:
                properties[line[start:end].strip()] = value
        aliases = []
        for key in (
            "compound_id", "CompoundID", "compound_name", "_Name", "ID", "id",
        ):
            value = properties.get(key)
            if value and value not in aliases:
                aliases.append(value)
        title = lines[0].strip() if lines else ""
        if title and title not in aliases:
            aliases.append(title)
        identifier = aliases[0] if aliases else f"compound_{record_index:06d}"
        row = {
            "source_molecule_index": record_index,
            "compound_id": identifier,
            "mol_name": title or identifier,
            "chemical_identity_id": properties.get("chemical_identity_id"),
            "chemical_id": properties.get("chemical_id"),
            "canonical_smiles": properties.get("canonical_smiles")
            or properties.get("SMILES"),
            "input_status": "text_parsed",
            "_input_aliases": aliases or [identifier],
        }
        row.update({
            field: properties[field]
            for field in LIGAND_LINEAGE_PROPERTIES
            if properties.get(field)
        })
        rows.append(row)
    return rows


def _molecule_property(molecule, name: str) -> str | None:
    if not molecule.HasProp(name):
        return None
    value = molecule.GetProp(name).strip()
    return value or None


def _score_source_index(
    score: dict,
    *,
    rows_by_index: dict[int, dict],
    rows_by_alias: dict[str, list[int]],
) -> int | None:
    """Resolve one score to an input without trusting a shifted engine index."""
    score_name = str(
        score.get("compound_id") or score.get("mol_name") or ""
    ).strip()
    raw_index = score.get("source_molecule_index")
    try:
        source_index = int(raw_index)
    except (TypeError, ValueError):
        source_index = None

    indexed = rows_by_index.get(source_index) if source_index is not None else None
    if indexed is not None and indexed.get("input_status") != "invalid":
        aliases = indexed.get("_input_aliases") or []
        if not score_name or score_name in aliases:
            return source_index

    alias_matches = rows_by_alias.get(score_name, []) if score_name else []
    if len(alias_matches) == 1:
        return alias_matches[0]

    # Engines that preserve only a numeric source index still have an
    # unambiguous mapping. A conflicting named identity is deliberately not
    # accepted here because that can attach a pose to the wrong molecule.
    if indexed is not None and not score_name:
        return source_index
    return None


def _diagnostic_reasons(result: dict, source_index: int, compound_id: str) -> list[str]:
    """Collect per-input reasons already emitted by an engine adapter."""
    reasons: list[str] = []

    def visit(value) -> None:
        if isinstance(value, dict):
            matches_index = value.get("source_molecule_index") == source_index
            matches_indices = source_index in (
                value.get("source_molecule_indices") or []
            ) if isinstance(value.get("source_molecule_indices"), list) else False
            matches_name = any(
                str(value.get(key) or "").strip() == compound_id
                for key in ("compound_id", "mol_name")
            )
            if matches_index or matches_indices or matches_name:
                reason = value.get("reason") or value.get("error")
                if reason:
                    reasons.append(str(reason))
                else:
                    details = {
                        key: item
                        for key, item in value.items()
                        if key not in {
                            "source_molecule_index", "source_molecule_indices",
                            "compound_id", "mol_name", "file",
                        }
                    }
                    if details:
                        reasons.append(json.dumps(
                            details, sort_keys=True, ensure_ascii=False,
                        ))
            for item in value.values():
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(result.get("diagnostics") or result.get("merge_diagnostics") or {})
    visit(result.get("operational_errors") or [])
    return list(dict.fromkeys(reasons))[:3]


def _publish_molecule_ledger(
    result: dict,
    *,
    ligand_file: Path,
    output_dir: Path,
    engine: str,
) -> None:
    """Publish a lossless input ledger for success, failure, and skipped rows.

    Poses remain the only rows eligible for downstream pose ranking. Molecules
    without a pose stay visible here as explicit docking failures instead of
    silently disappearing from the campaign.
    """
    rows = _input_molecule_rows(ligand_file)
    if not rows:
        return

    rows_by_index = {
        int(row["source_molecule_index"]): row for row in rows
    }
    rows_by_alias: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        if row.get("input_status") == "invalid":
            continue
        source_index = int(row["source_molecule_index"])
        for alias in row.get("_input_aliases") or []:
            rows_by_alias[alias].append(source_index)

    pose_counts: Counter[int] = Counter()
    unmatched_scores: list[dict] = []
    for score in result.get("scores") or []:
        if not isinstance(score, dict):
            continue
        source_index = _score_source_index(
            score,
            rows_by_index=rows_by_index,
            rows_by_alias=rows_by_alias,
        )
        if source_index is None:
            unmatched_scores.append({
                "pose_index": score.get("pose_index"),
                "source_molecule_index": score.get("source_molecule_index"),
                "compound_id": score.get("compound_id") or score.get("mol_name"),
                "reason": "score could not be joined unambiguously to an input molecule",
            })
            continue
        pose_counts[source_index] += 1

    published_rows: list[dict] = []
    failed_rows: list[dict] = []
    for row in rows:
        source_index = int(row["source_molecule_index"])
        pose_count = int(pose_counts.get(source_index, 0))
        published = {
            key: value for key, value in row.items() if not key.startswith("_")
        }
        published["engine"] = engine
        published["pose_count"] = pose_count
        if row.get("input_status") == "invalid":
            published["docking_status"] = "input_invalid"
            published["failure_reason"] = "input SDF record could not be parsed"
            published["retryable"] = True
            failed_rows.append(published)
        elif pose_count > 0:
            published["docking_status"] = "scored"
            published["failure_reason"] = None
            published["retryable"] = False
        else:
            reasons = _diagnostic_reasons(
                result, source_index, str(row["compound_id"]),
            )
            if not reasons:
                reasons = [
                    str(result.get("error") or "No finite scored pose was published")
                ]
            published["docking_status"] = "failed"
            published["failure_reason"] = "; ".join(reasons)
            published["retryable"] = True
            failed_rows.append(published)
        published_rows.append(published)

    status_counts = Counter(row["docking_status"] for row in published_rows)
    total = len(published_rows)
    accounted = sum(status_counts.values())
    summary = {
        "input_molecules": total,
        "scored_molecules": status_counts.get("scored", 0),
        "failed_molecules": status_counts.get("failed", 0),
        "input_invalid_molecules": status_counts.get("input_invalid", 0),
        "accounted_molecules": accounted,
        "accounting_complete": accounted == total,
        "unmatched_score_records": len(unmatched_scores),
    }

    json_path = output_dir / "docking_molecule_ledger.json"
    tsv_path = output_dir / "docking_molecule_ledger.tsv"
    json_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "engine": engine,
                "summary": summary,
                "molecules": published_rows,
                "unmatched_score_records": unmatched_scores,
            },
            indent=2,
            ensure_ascii=False,
            default=str,
        ) + "\n",
        encoding="utf-8",
    )
    columns = list(dict.fromkeys([
        "source_molecule_index", "compound_id", "mol_name", "chemical_identity_id",
        "chemical_id", "canonical_smiles", *LIGAND_LINEAGE_PROPERTIES,
        "engine", "input_status",
        "docking_status", "pose_count", "failure_reason", "retryable",
    ]))
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {key: row.get(key) for key in columns}
            for row in published_rows
        )

    result["molecule_ledger"] = published_rows
    result["molecule_ledger_summary"] = summary
    result["molecule_ledger_json"] = str(json_path)
    result["molecule_ledger_tsv"] = str(tsv_path)
    if unmatched_scores:
        result["molecule_ledger_unmatched_score_records"] = unmatched_scores

    if failed_rows:
        affected = [str(row["compound_id"]) for row in failed_rows]
        affected_indices = [
            int(row["source_molecule_index"]) for row in failed_rows
        ]
        errors = result.setdefault("operational_errors", [])
        aggregate = next(
            (
                item for item in errors
                if isinstance(item, dict) and item.get("tool") == engine
                and item.get("code") == "docking_input_unscored"
            ),
            None,
        )
        if aggregate is None:
            aggregate = {
                "tool": engine,
                "code": "docking_input_unscored",
                "scope": "molecule",
                "reason": (
                    f"{len(failed_rows)}/{total} input molecule(s) have no "
                    "finite scored pose"
                ),
            }
            errors.append(aggregate)
        aggregate["affected_molecules"] = affected
        aggregate["source_molecule_indices"] = affected_indices
        aggregate["row_ids"] = affected
        aggregate["retryable"] = True


def _selected_engine(input_data: dict) -> str:
    """Resolve exactly one primary engine from a complete explicit matrix."""
    tools = input_data.get("tools")
    if not isinstance(tools, dict):
        raise ValueError("tools is required for explicit docking tool selection")

    missing = [engine for engine in PRIMARY_ENGINES if engine not in tools]
    if missing:
        raise ValueError(
            "tools must explicitly configure every primary docking engine; "
            "missing: " + ", ".join(missing)
        )

    enabled = []
    for engine in PRIMARY_ENGINES:
        config = tools.get(engine)
        if not isinstance(config, dict):
            raise ValueError(f"tools.{engine} must be an object")
        if not isinstance(config.get("enabled"), bool):
            raise ValueError(f"tools.{engine}.enabled must be an explicit boolean")
        if config["enabled"]:
            enabled.append(engine)
    if len(enabled) != 1:
        raise ValueError(
            "docking requires exactly one primary engine enabled "
            f"from {{{', '.join(PRIMARY_ENGINES)}}}; got {enabled or 'none'}"
        )
    return enabled[0]


def _selected_redocking(input_data: dict) -> bool:
    tools = input_data.get("tools")
    config = tools.get("redocking_validation") if isinstance(tools, dict) else None
    if not isinstance(config, dict):
        raise ValueError("tools.redocking_validation must be explicitly configured")
    if not isinstance(config.get("enabled"), bool):
        raise ValueError(
            "tools.redocking_validation.enabled must be an explicit boolean"
        )
    return config["enabled"]


def _count_sdf_records_up_to(path: Path, limit: int) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.rstrip("\r\n") == "$$$$":
                count += 1
                if count >= limit:
                    break
    return count


def _engine_param(input_data: dict, engine: str, name: str, default):
    params = (tool_config(input_data, engine).get("params") or {})
    value = (
        params[name] if name in params and params[name] is not None
        else input_data.get(name, default)
    )
    if name in {"num_modes", "samples_per_complex"}:
        return min(int(value), int(input_data.get("_pose_quota", 20)))
    return value


def _redocking_threshold(input_data: dict) -> float:
    params = (tool_config(input_data, "redocking_validation").get("params") or {})
    value = params.get(
        "redocking_rmsd_threshold",
        input_data.get("redocking_rmsd_threshold", 2.0),
    )
    return float(value)


def _redocking_hard_threshold(input_data: dict) -> float:
    params = (tool_config(input_data, "redocking_validation").get("params") or {})
    value = params.get(
        "redocking_rmsd_hard_threshold",
        input_data.get("redocking_rmsd_hard_threshold", 2.5),
    )
    return float(value)


def _redocking_timeout(input_data: dict) -> int:
    params = (
        tool_config(input_data, "redocking_validation").get("params") or {}
    )
    value = params.get("timeout_seconds", 7200)
    return int(value)


def _run_gnina_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    params = {
        "autobox_add": float(_engine_param(input_data, "gnina", "autobox_add", 4.0)),
        "exhaustiveness": int(_engine_param(input_data, "gnina", "exhaustiveness", 16)),
        "num_modes": int(_engine_param(input_data, "gnina", "num_modes", 20)),
        "device": int(_engine_param(input_data, "gnina", "device", 0)),
        "seed": int(_engine_param(input_data, "gnina", "seed", 42)),
        "pose_sort_order": str(
            _engine_param(input_data, "gnina", "pose_sort_order", "CNNscore")
        ),
        "cnn_scoring": str(
            _engine_param(input_data, "gnina", "cnn_scoring", "rescore")
        ),
        "gpu_monitor_interval_seconds": int(
            _engine_param(
                input_data, "gnina", "gpu_monitor_interval_seconds", 5,
            )
        ),
        "timeout_seconds": int(
            _engine_param(input_data, "gnina", "timeout_seconds", 86400)
        ),
        "redocking_rmsd_threshold": _redocking_threshold(input_data),
        "redocking_rmsd_hard_threshold": _redocking_hard_threshold(input_data),
        "redocking_timeout_seconds": _redocking_timeout(input_data),
    }
    if not check_gnina_available():
        return ({
            "success": False,
            "error": (
                "GNINA is not installed or not in PATH. Install from "
                "https://github.com/gnina/gnina/releases"
            ),
        }, params)

    result = run_gnina_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        validate_redocking=validate_redocking,
        **params,
    )
    return result, params


def _run_unidock_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "unidock")
    params = {
        "autobox_add": float(
            _engine_param(input_data, "unidock", "autobox_add", 4.0)
        ),
        "batch_size": int(
            _engine_param(input_data, "unidock", "batch_size", 18000)
        ),
        "scoring_function": str(
            _engine_param(input_data, "unidock", "scoring_function", "vina")
        ),
        "search_mode": str(
            _engine_param(input_data, "unidock", "search_mode", "")
        ),
        "exhaustiveness": int(
            _engine_param(input_data, "unidock", "exhaustiveness", 16)
        ),
        "max_step": int(_engine_param(input_data, "unidock", "max_step", 20)),
        "num_modes": int(_engine_param(input_data, "unidock", "num_modes", 20)),
        "refine_step": int(
            _engine_param(input_data, "unidock", "refine_step", 3)
        ),
        "energy_range": float(
            _engine_param(input_data, "unidock", "energy_range", 3.0)
        ),
        "device": int(_engine_param(input_data, "unidock", "device", 0)),
        "seed": int(_engine_param(input_data, "unidock", "seed", 42)),
        "timeout_seconds": int(
            _engine_param(input_data, "unidock", "timeout_seconds", 86400)
        ),
    }
    result = run_unidock_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        runtime=config.get("runtime") or None,
        **params,
    )
    return result, params


def _run_autodock_vina_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "autodock_vina")
    params = {
        "autobox_add": float(
            _engine_param(input_data, "autodock_vina", "autobox_add", 4.0)
        ),
        "scoring_function": str(
            _engine_param(
                input_data, "autodock_vina", "scoring_function", "vina"
            )
        ),
        "exhaustiveness": int(
            _engine_param(
                input_data, "autodock_vina", "exhaustiveness", 16
            )
        ),
        "num_modes": int(
            _engine_param(input_data, "autodock_vina", "num_modes", 20)
        ),
        "energy_range": float(
            _engine_param(input_data, "autodock_vina", "energy_range", 3.0)
        ),
        "min_rmsd": float(
            _engine_param(input_data, "autodock_vina", "min_rmsd", 1.0)
        ),
        "max_evals": int(
            _engine_param(input_data, "autodock_vina", "max_evals", 0)
        ),
        "cpu": int(_engine_param(input_data, "autodock_vina", "cpu", 1)),
        "seed": int(_engine_param(input_data, "autodock_vina", "seed", 42)),
        "spacing": float(
            _engine_param(input_data, "autodock_vina", "spacing", 0.375)
        ),
        "force_even_voxels": bool(
            _engine_param(
                input_data,
                "autodock_vina",
                "force_even_voxels",
                False,
            )
        ),
        "no_refine": bool(
            _engine_param(input_data, "autodock_vina", "no_refine", False)
        ),
        "timeout_seconds": int(
            _engine_param(
                input_data, "autodock_vina", "timeout_seconds", 86400
            )
        ),
    }
    if receptor_pdbqt_file is None:
        return ({
            "success": False,
            "error": "receptor_pdbqt_file is required for AutoDock Vina",
        }, params)
    result = run_autodock_vina_docking(
        receptor_pdbqt_file=receptor_pdbqt_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        runtime=config.get("runtime") or None,
        **params,
    )
    return result, params


def _run_diffdock_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "diffdock")
    params = {
        "source_root": str(
            _engine_param(
                input_data, "diffdock", "source_root", str(deployment_path("${SBDD_DIFFDOCK_SOURCE_DIR}")),
            )
        ),
        "model_dir": _engine_param(
            input_data, "diffdock", "model_dir", None,
        ),
        "confidence_model_dir": _engine_param(
            input_data, "diffdock", "confidence_model_dir", None,
        ),
        "checkpoint": str(
            _engine_param(
                input_data,
                "diffdock",
                "checkpoint",
                "best_ema_inference_epoch_model.pt",
            )
        ),
        "confidence_checkpoint": str(
            _engine_param(
                input_data,
                "diffdock",
                "confidence_checkpoint",
                "best_model_epoch75.pt",
            )
        ),
        "samples_per_complex": int(
            _engine_param(input_data, "diffdock", "samples_per_complex", 20)
        ),
        "batch_size": int(
            _engine_param(input_data, "diffdock", "batch_size", 10)
        ),
        "inference_steps": int(
            _engine_param(input_data, "diffdock", "inference_steps", 20)
        ),
        "actual_steps": _engine_param(
            input_data, "diffdock", "actual_steps", None,
        ),
        "old_confidence_model": bool(
            _engine_param(
                input_data, "diffdock", "old_confidence_model", True,
            )
        ),
        "device": int(_engine_param(input_data, "diffdock", "device", 0)),
        "seed": int(_engine_param(input_data, "diffdock", "seed", 42)),
        "timeout_seconds": int(
            _engine_param(input_data, "diffdock", "timeout_seconds", 86400)
        ),
        "scoring_function": "diffdock_confidence",
    }
    call_params = {
        key: value for key, value in params.items()
        if key != "scoring_function"
    }
    result = run_diffdock_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        runtime=config.get("runtime") or None,
        **call_params,
    )
    return result, params


def _run_dock6_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "dock6")
    params = {
        "dock_root": str(
            _engine_param(input_data, "dock6", "dock_root", str(deployment_path("${SBDD_DOCK6_SOURCE_DIR}")))
        ),
        "antechamber_path": _engine_param(
            input_data, "dock6", "antechamber_path", None,
        ),
        "charge_method": str(
            _engine_param(input_data, "dock6", "charge_method", "gas")
        ),
        "num_modes": int(
            _engine_param(input_data, "dock6", "num_modes", 20)
        ),
        "max_orientations": int(
            _engine_param(input_data, "dock6", "max_orientations", 500)
        ),
        "min_anchor_size": int(
            _engine_param(input_data, "dock6", "min_anchor_size", 5)
        ),
        "pruning_max_orients": int(
            _engine_param(input_data, "dock6", "pruning_max_orients", 100)
        ),
        "score_threshold": float(
            _engine_param(input_data, "dock6", "score_threshold", 1000.0)
        ),
        "seed": int(_engine_param(input_data, "dock6", "seed", 42)),
        "timeout_seconds": int(
            _engine_param(input_data, "dock6", "timeout_seconds", 86400)
        ),
        "ligand_prep_timeout_seconds": int(
            _engine_param(
                input_data,
                "dock6",
                "ligand_prep_timeout_seconds",
                600,
            )
        ),
        "scoring_function": "dock6_grid_score",
    }
    call_params = {
        key: value for key, value in params.items()
        if key != "scoring_function"
    }
    result = run_dock6_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        receptor_site_file=input_data.get("dock6_receptor_site_file"),
        grid_nrg_file=input_data.get("dock6_grid_nrg_file"),
        grid_bmp_file=input_data.get("dock6_grid_bmp_file"),
        runtime=config.get("runtime") or None,
        **call_params,
    )
    return result, params


def _run_rdock_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "rdock")
    params = {
        "cavity_radius": float(
            _engine_param(input_data, "rdock", "cavity_radius", 6.0)
        ),
        "num_runs": int(
            _engine_param(input_data, "rdock", "num_runs", 25)
        ),
        "num_modes": int(
            _engine_param(input_data, "rdock", "num_modes", 20)
        ),
        "receptor_flex": float(
            _engine_param(input_data, "rdock", "receptor_flex", 3.0)
        ),
        "seed": int(_engine_param(input_data, "rdock", "seed", 42)),
        "timeout_seconds": int(
            _engine_param(input_data, "rdock", "timeout_seconds", 86400)
        ),
        "scoring_function": "rdock_score",
    }
    result = run_rdock_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        runtime=config.get("runtime") or None,
        **{
            key: value for key, value in params.items()
            if key != "scoring_function"
        },
    )
    return result, params


def _run_autodock4_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "autodock4")
    params = {
        "autobox_add": float(
            _engine_param(input_data, "autodock4", "autobox_add", 4.0)
        ),
        "spacing": float(
            _engine_param(input_data, "autodock4", "spacing", 0.375)
        ),
        "num_runs": int(
            _engine_param(input_data, "autodock4", "num_runs", 20)
        ),
        "num_modes": int(
            _engine_param(input_data, "autodock4", "num_modes", 20)
        ),
        "max_evals": int(
            _engine_param(input_data, "autodock4", "max_evals", 2500000)
        ),
        "population_size": int(
            _engine_param(input_data, "autodock4", "population_size", 150)
        ),
        "cluster_rmsd": float(
            _engine_param(input_data, "autodock4", "cluster_rmsd", 2.0)
        ),
        "seed": int(_engine_param(input_data, "autodock4", "seed", 42)),
        "timeout_seconds": int(
            _engine_param(input_data, "autodock4", "timeout_seconds", 86400)
        ),
        "scoring_function": "autodock4_binding_energy",
    }
    result = run_autodock4_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        runtime=config.get("runtime") or None,
        **{
            key: value for key, value in params.items()
            if key != "scoring_function"
        },
    )
    return result, params


def _run_adfr_engine(
    *,
    input_data: dict,
    receptor_file: Path,
    receptor_pdbqt_file: Path | None,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    validate_redocking: bool,
) -> tuple[dict, dict]:
    del receptor_pdbqt_file
    del validate_redocking  # rejected before dispatch; kept in runner protocol
    config = tool_config(input_data, "adfr")
    params = {
        "autobox_add": float(
            _engine_param(input_data, "adfr", "autobox_add", 4.0)
        ),
        "spacing": float(
            _engine_param(input_data, "adfr", "spacing", 0.375)
        ),
        "num_runs": int(
            _engine_param(input_data, "adfr", "num_runs", 20)
        ),
        "num_modes": int(
            _engine_param(input_data, "adfr", "num_modes", 20)
        ),
        "max_evals": int(
            _engine_param(input_data, "adfr", "max_evals", 2500000)
        ),
        "population_size": int(
            _engine_param(input_data, "adfr", "population_size", 150)
        ),
        "max_cores": int(
            _engine_param(input_data, "adfr", "max_cores", 1)
        ),
        "seed": int(_engine_param(input_data, "adfr", "seed", 42)),
        "timeout_seconds": int(
            _engine_param(input_data, "adfr", "timeout_seconds", 86400)
        ),
        "scoring_function": "adfr_feb",
    }
    result = run_adfr_docking(
        receptor_file=receptor_file,
        ligand_file=ligand_file,
        reference_ligand_file=reference_ligand_file,
        output_dir=output_dir,
        docking_box=docking_box,
        runtime=config.get("runtime") or None,
        **{
            key: value for key, value in params.items()
            if key != "scoring_function"
        },
    )
    return result, params


ENGINE_RUNNERS = {
    "gnina": _run_gnina_engine,
    "unidock": _run_unidock_engine,
    "autodock_vina": _run_autodock_vina_engine,
    "diffdock": _run_diffdock_engine,
    "dock6": _run_dock6_engine,
    "rdock": _run_rdock_engine,
    "autodock4": _run_autodock4_engine,
    "adfr": _run_adfr_engine,
}


def _annotate_score_provenance(
    result: dict,
    engine: str,
    engine_parameters: dict,
    run_parameters: dict,
) -> None:
    """Attach the protocol identity needed for safe cohort normalization."""
    scoring_function = engine_parameters.get("scoring_function")
    if scoring_function in (None, ""):
        scoring_function = {
            "gnina": "gnina",
            "diffdock": "diffdock_confidence",
            "dock6": "dock6_grid_score",
            "rdock": "rdock_score",
            "autodock4": "autodock4_binding_energy",
            "adfr": "adfr_feb",
        }.get(engine, "unknown")
    software = result.get("software") or {}
    version = (
        software.get("gnina_version")
        or software.get("engine_version")
        or software.get("version")
        or "unknown"
    )
    score_parameters = {
        key: value
        for key, value in engine_parameters.items()
        if key not in {
            "timeout_seconds",
            "redocking_timeout_seconds",
            "redocking_rmsd_threshold",
            "redocking_rmsd_hard_threshold",
            "gpu_monitor_interval_seconds",
            # Resource content identities live in protocol.engine_assets;
            # installation locations must not change the fingerprint.
            "source_root",
            "model_dir",
            "confidence_model_dir",
            "dock_root",
            "antechamber_path",
        }
    }
    receptor_identity = (
        None
        if engine in {"autodock_vina", "dock6"}
        else _file_identity(run_parameters.get("receptor_file"))
    )
    receptor_pdbqt_identity = (
        _file_identity(run_parameters.get("receptor_pdbqt_file"))
        if engine == "autodock_vina"
        else None
    )
    ligand_batch_identity = _file_identity(run_parameters.get("ligand_file"))
    engine_assets = result.get("protocol_assets")
    identity_warnings = []
    if str(version).strip().casefold() in {"", "unknown", "none"}:
        identity_warnings.append("missing_engine_version")
    if str(scoring_function).strip().casefold() in {"", "unknown", "none"}:
        identity_warnings.append("missing_scoring_function")
    if (
        receptor_identity is None
        and receptor_pdbqt_identity is None
        and not engine_assets
    ):
        identity_warnings.append("missing_scoring_receptor_identity")
    if ligand_batch_identity is None:
        identity_warnings.append("missing_ligand_batch_identity")
    protocol = {
        "engine": engine,
        "engine_version": str(version),
        "software_identity": software,
        "scoring_function": str(scoring_function),
        "receptor": receptor_identity,
        "receptor_pdbqt": receptor_pdbqt_identity,
        "ligand_batch": ligand_batch_identity,
        "ligand_preparation_manifest": _file_identity(
            run_parameters.get("ligand_preparation_manifest_file")
        ),
        "reference_ligand": _file_identity(
            run_parameters.get("reference_ligand_file")
        ),
        "engine_assets": engine_assets,
        "search_space": result.get("search_space"),
        "score_parameters": score_parameters,
        "protocol_identity_complete": not identity_warnings,
        "protocol_identity_warnings": identity_warnings,
    }
    fingerprint_payload = _without_audit_paths(protocol)
    fingerprint = hashlib.sha256(
        json.dumps(
            fingerprint_payload, sort_keys=True, default=str,
        ).encode("utf-8")
    ).hexdigest()
    protocol_id = f"{engine}:{fingerprint[:16]}"
    result["score_protocol"] = {
        **protocol,
        "protocol_id": protocol_id,
        "fingerprint_sha256": fingerprint,
    }
    for score in result.get("scores") or []:
        if not isinstance(score, dict):
            continue
        score.setdefault("engine", engine)
        score.setdefault("method", engine)
        score.setdefault("method_version", str(version))
        score.setdefault("stage", "docking")
        score.setdefault("score_reliability", "supported")
        score.setdefault("scoring_function", str(scoring_function))
        score.setdefault("protocol_id", protocol_id)
        score.setdefault(
            "protocol_identity_complete", not identity_warnings,
        )


def _attach_scientific_cohort(
    result: dict,
    input_data: dict,
    run_parameters: dict,
) -> dict[str, object]:
    """Finalize the five-field cohort after the real protocol ID is known."""
    protocol = result.get("score_protocol") or {}
    fingerprint = str(protocol.get("fingerprint_sha256") or "")
    if not fingerprint:
        fingerprint = hashlib.sha256(
            json.dumps(protocol, sort_keys=True, default=str).encode("utf-8")
        ).hexdigest()
    source_token = fingerprint

    receptor = protocol.get("receptor") or protocol.get("receptor_pdbqt") or {}
    receptor_digest = receptor.get("sha256") if isinstance(receptor, dict) else None
    search_space = protocol.get("search_space") or result.get("search_space") or {}
    search_space_digest = hashlib.sha256(
        json.dumps(search_space, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()

    input_protein_state = input_data.get("protein_state_id")
    if (
        str(input_protein_state or "").startswith("unknown:")
        or "|state_source:" in str(input_protein_state or "")
    ):
        input_protein_state = None
    input_pocket = input_data.get("pocket_id")
    if str(input_pocket or "").startswith("unknown:"):
        input_pocket = None
    components: dict[str, object] = {
        "protein_state_id": input_protein_state
        or (
            f"receptor_sha256:{receptor_digest}"
            if receptor_digest
            else f"unknown:protein_state_id:{source_token}"
        ),
        "pocket_id": input_pocket
        or f"search_space_sha256:{search_space_digest}",
        "round_id": input_data.get("round_id")
        or f"unknown:round_id:{source_token}",
        # The adapter may carry a pre-execution placeholder.  Once docking has
        # run, the content-derived score protocol is the authoritative profile.
        "protocol_profile_id": protocol.get("protocol_id")
        or input_data.get("protocol_profile_id")
        or f"unknown:protocol_profile_id:{source_token}",
        "metric_catalog_version": input_data.get("metric_catalog_version")
        or METRIC_RULES_VERSION,
    }
    cohort_id = "|".join(
        f"{field}={components[field]}"
        for field in SCIENTIFIC_COHORT_COMPONENT_FIELDS
    )
    upstream_derived = input_data.get("cohort_identity_derived_fields")
    derived_fields = (
        list(upstream_derived) if isinstance(upstream_derived, list) else []
    )
    derived_here = ["protocol_profile_id"]
    if input_protein_state is None:
        derived_here.append("protein_state_id")
    if input_pocket is None:
        derived_here.append("pocket_id")
    if input_data.get("round_id") in (None, ""):
        derived_here.append("round_id")
    if input_data.get("metric_catalog_version") in (None, ""):
        derived_here.append("metric_catalog_version")
    derived_fields.extend(derived_here)
    missing_fields = [
        field
        for field, value in components.items()
        if str(value).startswith(("unknown:", "pending:"))
        or "|state_source:" in str(value)
    ]
    upstream_sources = input_data.get("cohort_provenance_sources")
    provenance_sources = (
        dict(upstream_sources) if isinstance(upstream_sources, dict) else {}
    )
    provenance_sources.update({
        "protocol_profile_id": "docking:score_protocol_fingerprint",
        "protein_state_id": (
            provenance_sources.get("protein_state_id")
            or "docking:input_protein_state_id"
            if input_protein_state is not None
            else "docking:scoring_receptor_sha256"
        ),
        "pocket_id": (
            provenance_sources.get("pocket_id")
            or "docking:input_pocket_id"
            if input_pocket is not None
            else "docking:search_space_sha256"
        ),
        "round_id": provenance_sources.get("round_id")
        or "docking:input_round_id",
        "metric_catalog_version": provenance_sources.get(
            "metric_catalog_version"
        ) or "registry:metric_catalog",
    })
    cohort = {
        **components,
        "cohort_id": cohort_id,
        "cohort_identity_complete": not missing_fields,
        "cohort_identity_missing_fields": missing_fields,
        "cohort_identity_derived_fields": sorted(set(derived_fields)),
        "cohort_provenance_sources": provenance_sources,
    }
    result.update(cohort)
    result["scientific_cohort"] = dict(cohort)
    run_parameters.update(cohort)
    for score in result.get("scores") or []:
        if isinstance(score, dict):
            score.update(cohort)
    docked_sdf = result.get("docked_sdf")
    if docked_sdf:
        result["docked_pose_artifact_contract"] = _set_sdf_cohort_properties(
            Path(str(docked_sdf)),
            cohort,
            scores=result.get("scores") or [],
            ligand_file=Path(run_parameters["ligand_file"]) if run_parameters.get("ligand_file") else None,
        )
    return cohort


def _read_sdf_record_blocks(path: Path) -> list[str]:
    """Read physical SDF records using delimiter lines, never substrings."""
    blocks: list[str] = []
    current: list[str] = []
    saw_delimiter = False
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.rstrip("\r\n") == "$$$$":
                saw_delimiter = True
                # Every delimiter closes one physical record.  Keep even an
                # empty record so validation fails at its real index instead
                # of silently shifting every later score one position left.
                blocks.append("".join(current))
                current = []
            else:
                current.append(line)
    if current and any(part.strip() for part in current):
        raise DockedPoseContractError(
            "unterminated_sdf_record",
            f"docked SDF has data after its last $$$$ delimiter: {path}",
        )
    if not saw_delimiter or not blocks:
        raise DockedPoseContractError(
            "empty_or_unframed_docked_sdf",
            f"docked SDF contains no framed molecule records: {path}",
        )
    return blocks


def _score_index(scores: list[dict], record_count: int) -> dict[int, dict]:
    """Build a strict index; duplicate/missing/out-of-range keys are fatal."""
    indexed: dict[int, dict] = {}
    for row_number, score in enumerate(scores):
        if not isinstance(score, dict):
            raise DockedPoseContractError(
                "invalid_docking_score_record",
                f"docking score row {row_number} is not an object",
            )
        if score.get("pose_index") in (None, ""):
            raise DockedPoseContractError(
                "missing_docking_pose_index",
                f"docking score row {row_number} has no pose_index",
            )
        try:
            pose_index = _strict_integer(score["pose_index"])
        except ValueError as exc:
            raise DockedPoseContractError(
                "invalid_docking_pose_index",
                f"docking score row {row_number} has invalid pose_index",
            ) from exc
        if pose_index < 0 or pose_index >= record_count:
            raise DockedPoseContractError(
                "docking_pose_index_out_of_range",
                f"docking score pose_index {pose_index} is outside "
                f"0..{record_count - 1}",
                details={"row_number": row_number, "pose_index": pose_index},
            )
        if pose_index in indexed:
            raise DockedPoseContractError(
                "duplicate_docking_pose_index",
                f"duplicate docking score pose_index: {pose_index}",
            )
        indexed[pose_index] = score
    return indexed


def _strict_integer(value: object) -> int:
    """Parse a JSON integer without truncating floats or accepting booleans."""
    if isinstance(value, bool):
        raise ValueError("boolean is not an integer index")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped and stripped.lstrip("-").isdigit():
            return int(stripped)
    raise ValueError("value is not a canonical integer")


def _validate_docked_pose_contract(
    path: Path,
    scores: list[dict],
    *,
    expected_record_count: int | None = None,
    ligand_file: Path | None = None,
    molecule_failures: dict[int, str] | None = None,
) -> dict[str, object]:
    """Prove that physical SDF records and score rows share one identity map."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise DockedPoseContractError(
            "rdkit_unavailable_for_docked_pose_validation",
            "RDKit is required to validate docked pose/score identity",
        ) from exc

    blocks = _read_sdf_record_blocks(path)
    if expected_record_count is not None and len(blocks) != expected_record_count:
        raise DockedPoseContractError(
            "docked_sdf_record_count_changed",
            "cohort annotation changed the number of physical SDF records",
            details={"before": expected_record_count, "after": len(blocks)},
        )
    scores_by_index = _score_index(scores, len(blocks))
    parsed_molecules = []
    for pose_index, block in enumerate(blocks):
        framed = block if block.endswith(("\n", "\r")) else block + "\n"
        payload = (framed + "$$$$\n").encode("utf-8")
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(payload), removeHs=False, strictParsing=False,
        )
        molecules = list(supplier)
        if len(molecules) != 1 or molecules[0] is None:
            error = DockedPoseContractError(
                "docked_sdf_record_parse_failed",
                f"RDKit cannot parse physical docked SDF record {pose_index}",
                details={"pose_index": pose_index},
            )
            # Only sanitization failures retain enough structure/metadata to
            # prove a local failure. Broken framing or unreadable atoms remain
            # a batch error, even when their title looks familiar.
            if molecule_failures is None:
                raise error
            molecules = list(Chem.ForwardSDMolSupplier(
                io.BytesIO(payload), removeHs=False, sanitize=False, strictParsing=False,
            ))
            if len(molecules) != 1 or molecules[0] is None:
                raise error
            molecule_failures[pose_index] = str(error)
        molecule = molecules[0]
        for field in ("pose_index", "source_pose_index"):
            if molecule.HasProp(field):
                try:
                    valid_index = _strict_integer(molecule.GetProp(field)) == pose_index
                except ValueError:
                    valid_index = False
                if not valid_index:
                    raise DockedPoseContractError(
                        "pose_score_identity_mismatch",
                        f"SDF {field} disagrees with physical record {pose_index}",
                    )
        if molecule.GetNumAtoms() == 0 or molecule.GetNumConformers() != 1:
            error = DockedPoseContractError(
                "invalid_docked_pose_geometry",
                "docked SDF record must contain a non-empty molecule with "
                f"exactly one conformer at pose_index {pose_index}",
                details={
                    "pose_index": pose_index,
                    "atom_count": molecule.GetNumAtoms(),
                    "conformer_count": molecule.GetNumConformers(),
                },
            )
            if molecule_failures is None:
                raise error
            molecule_failures[pose_index] = str(error)
            parsed_molecules.append(molecule)
            continue
        conformer = molecule.GetConformer()
        for atom_index in range(molecule.GetNumAtoms()):
            position = conformer.GetAtomPosition(atom_index)
            if not all(
                math.isfinite(value)
                for value in (position.x, position.y, position.z)
            ):
                error = DockedPoseContractError(
                    "invalid_docked_pose_geometry",
                    "docked SDF record contains non-finite coordinates at "
                    f"pose_index {pose_index}",
                    details={
                        "pose_index": pose_index,
                        "atom_index": atom_index,
                    },
                )
                if molecule_failures is None:
                    raise error
                molecule_failures[pose_index] = str(error)
        parsed_molecules.append(molecule)

    properties_by_index = [
        {
            str(field): molecule.GetProp(str(field)).strip()
            for field in molecule.GetPropNames(
                includePrivate=False, includeComputed=False,
            )
        }
        for molecule in parsed_molecules
    ]
    molblock_digests = [sdf_molblock_sha256(block) for block in blocks]
    names_by_index = []
    for pose_index, (molecule, properties) in enumerate(
        zip(parsed_molecules, properties_by_index, strict=True)
    ):
        names = {
            value
            for value in (
                molecule.GetProp("_Name").strip()
                if molecule.HasProp("_Name") else "",
                properties.get("mol_name", ""),
            )
            if value
        }
        if not names:
            names.add(f"pose_{pose_index}")
        names_by_index.append(names)

    for pose_index, molecule in enumerate(parsed_molecules):
        score = scores_by_index.get(pose_index)
        if score is None:
            # A genuinely absent score remains NA; it must never borrow a later
            # row. Coverage policy decides whether the run may continue.
            continue

        current_properties = properties_by_index[pose_index]
        candidate_indices = set(range(len(parsed_molecules)))
        evidence_fields: list[str] = []
        strong_identity_matched = False

        score_molblock_digest = score.get("sdf_molblock_sha256")
        if score_molblock_digest not in (None, ""):
            score_molblock_digest = str(score_molblock_digest).strip()
            if score_molblock_digest != molblock_digests[pose_index]:
                raise DockedPoseContractError(
                    "pose_score_identity_mismatch",
                    "docking score mol-block witness disagrees with physical "
                    f"record at pose_index {pose_index}",
                    details={"pose_index": pose_index},
                )
            candidate_indices &= {
                index
                for index, digest in enumerate(molblock_digests)
                if digest == score_molblock_digest
            }
            evidence_fields.append("sdf_molblock_sha256")
            strong_identity_matched = True

        def apply_exact_field(field: str, *, integer: bool = False) -> None:
            nonlocal candidate_indices, strong_identity_matched
            score_value = score.get(field)
            sdf_value = current_properties.get(field)
            if score_value in (None, "") or sdf_value in (None, ""):
                return
            try:
                normalized_score = (
                    _strict_integer(score_value)
                    if integer else str(score_value).strip()
                )
                normalized_sdf = (
                    _strict_integer(sdf_value)
                    if integer else str(sdf_value).strip()
                )
            except ValueError as exc:
                raise DockedPoseContractError(
                    "invalid_docked_pose_metadata",
                    f"invalid {field} at pose_index {pose_index}",
                    details={"pose_index": pose_index, "field": field},
                ) from exc
            if integer:
                minimum = 1 if field == "engine_pose_rank" else 0
                if normalized_score < minimum or normalized_sdf < minimum:
                    raise DockedPoseContractError(
                        "invalid_docked_pose_metadata",
                        f"{field} is outside its valid range at pose_index "
                        f"{pose_index}",
                        details={"pose_index": pose_index, "field": field},
                    )
            if normalized_sdf != normalized_score:
                raise DockedPoseContractError(
                    "pose_score_identity_mismatch",
                    f"SDF and score {field} disagree at pose_index {pose_index}",
                    details={
                        "pose_index": pose_index,
                        "field": field,
                        "sdf_value": sdf_value,
                        "score_value": score_value,
                    },
                )
            matching = set()
            for index, properties in enumerate(properties_by_index):
                candidate_value = properties.get(field)
                if candidate_value in (None, ""):
                    continue
                try:
                    candidate_value = (
                        _strict_integer(candidate_value)
                        if integer else str(candidate_value).strip()
                    )
                except ValueError:
                    continue
                if candidate_value == normalized_score:
                    matching.add(index)
            candidate_indices &= matching
            evidence_fields.append(field)
            if field in {
                "pose_id", "chemical_identity_id", "chemical_id",
                "source_identity_id",
            }:
                strong_identity_matched = True

        for field in (
            "pose_id", "chemical_identity_id", "chemical_id",
            "source_identity_id", "compound_id",
        ):
            apply_exact_field(field)
        for field in (
            "source_molecule_index", "preparation_variant_index",
            "engine_pose_index",
            "engine_pose_rank",
        ):
            apply_exact_field(field, integer=True)

        score_display_names = {
            str(score[field]).strip()
            for field in ("mol_name", "name")
            if score.get(field) not in (None, "")
            and str(score[field]).strip()
        }
        if (
            score.get("compound_id") not in (None, "")
            and current_properties.get("compound_id") in (None, "")
        ):
            score_display_names.add(str(score["compound_id"]).strip())
        if score_display_names:
            current_names = names_by_index[pose_index]
            if current_names & score_display_names:
                candidate_indices &= {
                    index for index, names in enumerate(names_by_index)
                    if names & score_display_names
                }
                evidence_fields.append("display_name")
            elif not strong_identity_matched:
                raise DockedPoseContractError(
                    "pose_score_identity_mismatch",
                    "docked SDF and score molecule labels disagree at "
                    f"pose_index {pose_index}",
                    details={
                        "pose_index": pose_index,
                        "sdf_names": sorted(current_names),
                        "score_names": sorted(score_display_names),
                    },
                )

        for field in LEGACY_SCORE_SPECS:
            if field not in score or score[field] is None:
                continue
            value = score[field]
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
            ):
                raise DockedPoseContractError(
                    "invalid_docking_score_value",
                    f"docking score {field} must be a finite JSON number at "
                    f"pose_index {pose_index}",
                    details={"pose_index": pose_index, "field": field},
                )

        # Numeric score properties catch swaps between two poses of the same
        # compound, where molecule names alone cannot distinguish them.
        for field, raw_score_value in score.items():
            if isinstance(raw_score_value, bool) or not isinstance(
                raw_score_value, (int, float)
            ):
                continue
            if not math.isfinite(float(raw_score_value)):
                raise DockedPoseContractError(
                    "non_finite_docking_score",
                    f"docking score {field} is non-finite at pose_index "
                    f"{pose_index}",
                    details={"pose_index": pose_index, "field": field},
                )
            if field in {
                "pose_index", "source_pose_index", "source_molecule_index",
                "engine_pose_index", "engine_pose_rank", "round_id",
            } or not molecule.HasProp(field):
                continue
            try:
                sdf_value = float(molecule.GetProp(field))
                score_value = float(raw_score_value)
            except (TypeError, ValueError):
                continue
            if not (
                math.isfinite(sdf_value)
                and math.isfinite(score_value)
                and math.isclose(
                    sdf_value, score_value, rel_tol=1e-7, abs_tol=1e-7,
                )
            ):
                raise DockedPoseContractError(
                    "pose_score_value_mismatch",
                    f"SDF and JSON {field} disagree at pose_index {pose_index}",
                    details={
                        "pose_index": pose_index,
                        "field": field,
                        "sdf_value": sdf_value,
                        "score_value": score_value,
                    },
                )
            matching = set()
            for index, properties in enumerate(properties_by_index):
                candidate_value = properties.get(field)
                if candidate_value in (None, ""):
                    continue
                try:
                    candidate_number = float(candidate_value)
                except (TypeError, ValueError):
                    continue
                if (
                    math.isfinite(candidate_number)
                    and math.isclose(
                        candidate_number,
                        score_value,
                        rel_tol=1e-7,
                        abs_tol=1e-7,
                    )
                ):
                    matching.add(index)
            candidate_indices &= matching
            evidence_fields.append(field)

        for structured_index, structured in enumerate(
            score.get("pose_scores") or []
        ):
            if not isinstance(structured, dict):
                continue
            value = structured.get("value")
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
            ):
                raise DockedPoseContractError(
                    "invalid_docking_score_value",
                    "structured docking score must be a finite JSON number at "
                    f"{pose_index}",
                    details={
                        "pose_index": pose_index,
                        "structured_score_index": structured_index,
                    },
                )
            nested_pose_index = structured.get("pose_index")
            if nested_pose_index not in (None, ""):
                try:
                    nested_pose_index = _strict_integer(nested_pose_index)
                except ValueError as exc:
                    raise DockedPoseContractError(
                        "invalid_docking_pose_index",
                        "structured docking score has invalid pose_index",
                        details={"pose_index": pose_index},
                    ) from exc
                if nested_pose_index != pose_index:
                    raise DockedPoseContractError(
                        "pose_score_identity_mismatch",
                        "structured docking score pose_index disagrees with "
                        f"its outer row at pose_index {pose_index}",
                        details={
                            "pose_index": pose_index,
                            "structured_pose_index": nested_pose_index,
                        },
                    )
            nested_pose_id = structured.get("pose_id")
            outer_pose_id = score.get("pose_id")
            if (
                nested_pose_id not in (None, "")
                and outer_pose_id not in (None, "")
                and str(nested_pose_id) != str(outer_pose_id)
            ):
                raise DockedPoseContractError(
                    "pose_score_identity_mismatch",
                    "structured docking score pose_id disagrees with its "
                    f"outer row at pose_index {pose_index}",
                )
            current_pose_id = current_properties.get("pose_id")
            if nested_pose_id not in (None, "") and current_pose_id:
                if str(nested_pose_id) != str(current_pose_id):
                    raise DockedPoseContractError(
                        "pose_score_identity_mismatch",
                        "structured docking score pose_id disagrees with its "
                        f"SDF record at pose_index {pose_index}",
                    )
                candidate_indices &= {
                    index
                    for index, properties in enumerate(properties_by_index)
                    if str(properties.get("pose_id") or "")
                    == str(nested_pose_id)
                }
                evidence_fields.append("pose_scores.pose_id")
            metric_id = str(structured.get("metric_id") or "")
            flat_field = next(
                (
                    field for field, spec in LEGACY_SCORE_SPECS.items()
                    if spec.get("metric_id") == metric_id
                ),
                None,
            )
            if flat_field and value is not None:
                flat_value = score.get(flat_field)
                if flat_value is not None and not math.isclose(
                    float(flat_value), float(value),
                    rel_tol=1e-7, abs_tol=1e-7,
                ):
                    raise DockedPoseContractError(
                        "pose_score_value_mismatch",
                        "structured and flat docking scores disagree at "
                        f"pose_index {pose_index}",
                        details={
                            "pose_index": pose_index,
                            "metric_id": metric_id,
                            "flat_field": flat_field,
                        },
                    )
                sdf_metric_value = current_properties.get(flat_field)
                if sdf_metric_value not in (None, ""):
                    try:
                        sdf_metric_number = float(sdf_metric_value)
                    except (TypeError, ValueError) as exc:
                        raise DockedPoseContractError(
                            "invalid_docked_pose_metadata",
                            f"SDF {flat_field} is not numeric at pose_index "
                            f"{pose_index}",
                        ) from exc
                    if not (
                        math.isfinite(sdf_metric_number)
                        and math.isclose(
                            sdf_metric_number,
                            float(value),
                            rel_tol=1e-7,
                            abs_tol=1e-7,
                        )
                    ):
                        raise DockedPoseContractError(
                            "pose_score_value_mismatch",
                            "structured docking score and SDF property "
                            f"disagree at pose_index {pose_index}",
                        )
                    matching = set()
                    for index, properties in enumerate(properties_by_index):
                        candidate_value = properties.get(flat_field)
                        if candidate_value in (None, ""):
                            continue
                        try:
                            candidate_number = float(candidate_value)
                        except (TypeError, ValueError):
                            continue
                        if (
                            math.isfinite(candidate_number)
                            and math.isclose(
                                candidate_number,
                                float(value),
                                rel_tol=1e-7,
                                abs_tol=1e-7,
                            )
                        ):
                            matching.add(index)
                    candidate_indices &= matching
                    evidence_fields.append(f"pose_scores.{metric_id}")

        if not evidence_fields or candidate_indices != {pose_index}:
            raise DockedPoseContractError(
                "insufficient_pose_identity_evidence",
                "score evidence does not uniquely identify its declared "
                f"physical record at pose_index {pose_index}",
                details={
                    "pose_index": pose_index,
                    "evidence_fields": evidence_fields,
                    "candidate_indices": sorted(candidate_indices),
                },
            )
    for pose_index, score in scores_by_index.items():
        molecule = parsed_molecules[pose_index]
        for field in LIGAND_LINEAGE_PROPERTIES:
            if molecule.HasProp(field):
                score.setdefault(field, molecule.GetProp(field))

    if ligand_file is not None:
        source_rows = _input_molecule_rows(ligand_file)
        by_index = {row["source_molecule_index"]: row for row in source_rows}
        declared_molecules = {
            index: Chem.MolFromSmiles(row["canonical_smiles"])
            if row["canonical_smiles"] else None
            for index, row in by_index.items()
        }
        by_alias: dict[str, list[int]] = defaultdict(list)
        for row in source_rows:
            for alias in row.get("_input_aliases") or []:
                by_alias[alias].append(row["source_molecule_index"])
        source_indices = {}
        indices = range(len(parsed_molecules)) if molecule_failures is not None else scores_by_index
        for pose_index in indices:
            actual = parsed_molecules[pose_index]
            score = scores_by_index.get(pose_index) or {
                "mol_name": compound_name(actual, pose_index),
            }
            source_index = _score_source_index(
                score, rows_by_index=by_index, rows_by_alias=by_alias,
            )
            if molecule_failures is not None:
                # Recovery requires agreement of raw record identity and score
                # identity; never infer a failed molecule from its batch slot.
                aliases = [actual.GetProp(key).strip() for key in ("compound_id", "_Name")
                           if actual.HasProp(key) and actual.GetProp(key).strip()]
                if not aliases or any(by_alias.get(alias) != [source_index] for alias in aliases):
                    source_index = None
                for identity in (score, properties_by_index[pose_index]):
                    if identity.get("source_molecule_index") not in (None, ""):
                        try:
                            if _strict_integer(identity["source_molecule_index"]) != source_index:
                                source_index = None
                        except ValueError:
                            source_index = None
            if source_index is None:
                raise DockedPoseContractError(
                    "unmatched_docking_score_identity",
                    f"docking score cannot be linked to one input: {pose_index}",
                )
            source_indices[pose_index] = source_index
            if molecule_failures is not None and pose_index in molecule_failures:
                continue
            if not matches_declared_molecule(actual, declared_molecules[source_index]):
                error = DockedPoseContractError(
                    "docked_molecule_identity_changed",
                    "docking changed the input molecular graph or specified "
                    f"stereochemistry at pose_index {pose_index}",
                    details={
                        "pose_index": pose_index,
                        "compound_id": by_index[source_index]["compound_id"],
                        "source_molecule_index": source_index,
                        "input_smiles": by_index[source_index]["canonical_smiles"],
                        "docked_smiles": Chem.MolToSmiles(Chem.RemoveHs(actual)),
                    },
                )
                if molecule_failures is None:
                    raise error
                molecule_failures[pose_index] = str(error)

    return {
        "status": "failed" if molecule_failures else "passed",
        **({"source_molecule_indices": source_indices} if molecule_failures is not None and ligand_file else {}),
        "physical_sdf_records": len(blocks),
        "rdkit_parsed_records": len(parsed_molecules),
        "score_records": len(scores_by_index),
        "unscored_sdf_records": len(blocks) - len(scores_by_index),
        "join_key": "pose_index_with_identity_validation",
    }


def _isolate_invalid_docked_molecules(
    result: dict, ligand_file: Path, output_dir: Path, engine: str,
) -> None:
    """Isolate proven local chemistry failures, preserving the raw batch.

    The ordinary validator still checks every score/record join, including
    quarantined rows. Only chemical/coordinate failures may be collected;
    identity, framing and score errors must fail the entire run.
    """
    source = Path(result["docked_sdf"])
    scores = result.get("scores") or []
    failures: dict[int, str] = {}
    contract = _validate_docked_pose_contract(
        source, scores, ligand_file=ligand_file, molecule_failures=failures,
    )
    if not failures:
        raise DockedPoseContractError("docking_recovery_unavailable", "No local molecule failure was established")
    sources = contract["source_molecule_indices"]
    excluded = {sources[index] for index in failures}
    blocks = _read_sdf_record_blocks(source)
    rows = {row["source_molecule_index"]: row for row in _input_molecule_rows(ligand_file)}
    retained = [index for index in range(len(blocks)) if sources[index] not in excluded]
    mapping = {old: new for new, old in enumerate(retained)}
    recovered_scores = []
    for original in scores:
        old = _strict_integer(original["pose_index"])
        for row in [original, *(original.get("pose_scores") or [])]:
            try:
                valid_source = "source_pose_index" not in row or _strict_integer(row["source_pose_index"]) == old
            except ValueError:
                valid_source = False
            if not valid_source:
                raise DockedPoseContractError(
                    "pose_score_identity_mismatch", "Score source_pose_index disagrees with its raw record",
                )
        if old not in mapping:
            continue
        # JSON cloning also keeps nested score identity tied to its new physical
        # record. The raw coordinate witness is unchanged by SD properties.
        score = json.loads(json.dumps(original, allow_nan=False))
        for row in [score, *(score.get("pose_scores") or [])]:
            if "pose_index" in row:
                row["pose_index"] = mapping[old]
            if "source_pose_index" in row:
                row["source_pose_index"] = mapping[old]
        score["raw_pose_index"] = old
        recovered_scores.append(score)
    for source_index in sorted(excluded):
        raw_indices = [index for index in sources if sources[index] == source_index]
        reason = next(failures[index] for index in raw_indices if index in failures)
        compound_id = str(rows[source_index]["compound_id"])
        result.setdefault("operational_errors", []).append({
            "tool": engine, "step": "docked_pose_artifact_contract",
            "code": "docked_molecule_quarantined", "scope": "molecule",
            "retryable": True, "reason": reason,
            "affected_molecules": [compound_id], "row_ids": [compound_id],
            "source_molecule_indices": [source_index], "raw_pose_indices": raw_indices,
        })
    result["raw_docked_sdf"] = str(source)
    manifest = {
        "raw_docked_sdf_artifact": _file_identity(source),
        "physical_record_count": len(blocks),
        "quarantined_molecule_indices": sorted(excluded),
        "quarantined_pose_indices": [index for index in sources if sources[index] in excluded],
        "invalid_pose_reasons": failures,
        "pose_index_to_raw_pose_index": retained,
    }
    manifest_path = output_dir / "docking_recovery.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    result["docking_recovery_manifest"] = str(manifest_path)
    result["docking_recovery"] = manifest
    result.setdefault("score_coverage", {}).update({
        "num_quarantined_poses": len(blocks) - len(retained),
        "num_published_score_records": len(recovered_scores),
    })
    if not recovered_scores:
        result["scores"] = []
        result["num_poses"] = 0
        raise DockedPoseContractError("no_usable_docked_poses", "No scored pose remains after molecule quarantine")

    rewritten = []
    for old in retained:
        lines = blocks[old].rstrip("\r\n").splitlines()
        kept = []
        had_source_index = False
        index = 0
        while index < len(lines):
            line = lines[index]
            field = line.partition("<")[2].partition(">")[0].strip() if line.startswith(">") else ""
            if field in {"pose_index", "source_pose_index", "raw_pose_index"}:
                had_source_index |= field == "source_pose_index"
                index += 1
                while index < len(lines) and lines[index].strip():
                    index += 1
            else:
                kept.append(line)
            index += 1
        kept.extend(["", "> <pose_index>", str(mapping[old]), "", "> <raw_pose_index>", str(old), ""])
        if had_source_index:
            kept.extend(["> <source_pose_index>", str(mapping[old]), ""])
        rewritten.append("\n".join(kept).rstrip("\n") + "\n\n$$$$\n")
    filtered = source.with_name(f"{uuid.uuid4().hex}-recovered.sdf")
    filtered.write_text("".join(rewritten), encoding="utf-8")
    _validate_docked_pose_contract(filtered, recovered_scores, ligand_file=ligand_file)
    result["docked_sdf"] = str(filtered)
    result["scores"] = recovered_scores
    result["num_poses"] = len(recovered_scores)


def _set_sdf_cohort_properties(
    path: Path,
    cohort: dict[str, object],
    *,
    scores: list[dict] | None = None,
    ligand_file: Path | None = None,
) -> dict[str, object]:
    """Replace cohort properties atomically and validate pose/score identity."""
    records = _read_sdf_record_blocks(path)
    rewritten: list[str] = []
    names = {*SCIENTIFIC_COHORT_OUTPUT_FIELDS, "pose_index"}
    for record_index, record in enumerate(records):
        lines = record.rstrip("\r\n").splitlines()
        kept: list[str] = []
        index = 0
        while index < len(lines):
            line = lines[index]
            property_name = ""
            if line.startswith(">") and "<" in line and ">" in line:
                start = line.find("<") + 1
                end = line.find(">", start)
                property_name = line[start:end].strip() if end > start else ""
            if property_name in names:
                if property_name == "pose_index":
                    existing = lines[index + 1] if index + 1 < len(lines) else ""
                    try:
                        existing_index = _strict_integer(existing)
                    except ValueError as exc:
                        raise DockedPoseContractError(
                            "invalid_docking_pose_index",
                            "SDF pose_index property is not a canonical integer "
                            f"at physical record {record_index}",
                        ) from exc
                    if existing_index != record_index:
                        raise DockedPoseContractError(
                            "pose_score_identity_mismatch",
                            "SDF pose_index property disagrees with physical "
                            f"record {record_index}",
                            details={
                                "pose_index": record_index,
                                "sdf_pose_index": existing_index,
                            },
                        )
                index += 1
                while index < len(lines) and lines[index].strip():
                    index += 1
                if index < len(lines):
                    index += 1
                continue
            kept.append(line)
            index += 1
        while kept and not kept[-1].strip():
            kept.pop()
        for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS:
            value = cohort[field]
            if isinstance(value, (dict, list, bool)):
                value = json.dumps(
                    value, sort_keys=True, separators=(",", ":")
                )
            kept.extend(["", f"> <{field}>", str(value), ""])
        kept.extend(["", "> <pose_index>", str(record_index), ""])
        # SD properties are terminated by a blank line before the record
        # delimiter.  Without it RDKit can silently drop the final molecule.
        rewritten.append("\n".join(kept).rstrip("\n") + "\n\n$$$$\n")
    temporary = path.with_name(path.name + ".cohort.tmp")
    try:
        temporary.write_text("".join(rewritten), encoding="utf-8")
        contract = _validate_docked_pose_contract(
            temporary,
            list(scores or []),
            expected_record_count=len(records),
            ligand_file=ligand_file,
        )
        # Validation above first proves that every engine score row belongs to
        # the declared physical SDF record.  Only after that proof succeeds do
        # we mint a coordinate-bearing witness for every score row.  Requiring
        # this witness downstream prevents a Step 4 report from one docking run
        # being joined to a same-name/same-index score from another run.
        rewritten_records = _read_sdf_record_blocks(temporary)
        for pose_index, score in _score_index(
            list(scores or []), len(rewritten_records),
        ).items():
            score["sdf_molblock_sha256"] = sdf_molblock_sha256(
                rewritten_records[pose_index]
            )
        temporary.replace(path)
        return contract
    finally:
        if temporary.exists():
            temporary.unlink()


def _file_identity(path_value: object) -> dict | None:
    """Return an auditable content identity for protocol-defining files."""
    if path_value in (None, ""):
        return None
    path = Path(str(path_value))
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": path.stat().st_size,
    }


def _without_audit_paths(value: object) -> object:
    """Remove filesystem locations while retaining content identities."""
    if isinstance(value, dict):
        return {
            key: _without_audit_paths(item)
            for key, item in value.items()
            if key != "path"
        }
    if isinstance(value, list):
        return [_without_audit_paths(item) for item in value]
    return value


def _write_report(output_dir: Path, report: dict):
    report_file = output_dir / "docking_report.json"
    temporary = output_dir / (
        f".{uuid.uuid4().hex}.docking_report.json.publish.tmp"
    )
    try:
        temporary.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
                default=str,
                allow_nan=False,
            ),
            encoding="utf-8",
        )
        temporary.replace(report_file)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    from io_utils import parse_skill_args, print_skill_result
    input_data, output_dir = parse_skill_args("Docking Skill")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

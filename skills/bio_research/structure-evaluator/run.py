#!/usr/bin/env python3
"""Structure Evaluator Skill - Multi-tool pose quality assessment.

Usage:
    python run.py --input '{"receptor_file": "protein.pdb", "docked_sdf": "docked.sdf", "interaction_blueprint_file": "blueprint.json"}' --output-dir ./output/
    python run.py --input-file params.json --output-dir ./output/
"""
from __future__ import annotations

import argparse
import csv
import importlib.metadata
import io
import json
import math
import sys
from time import monotonic
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from posebusters_utils import check_posebusters_available, run_posebusters
from posecheck_utils import check_posecheck_available, run_posecheck
from plip_eval_utils import run_plip_per_pose
from prolif_eval_utils import run_prolif_per_pose
from blueprint_match import (
    load_blueprint,
    match_pose_to_blueprint,
    unresolved_match_result,
)
from pocket_occupancy import compute_pocket_occupancy, load_subpockets, prepare_subpocket_coords
from evaluation_runtime import cached_file, profiler_worker_count, run_profiler_tasks, write_json
from diagnosis_engine import generate_diagnosis
from constants import (
    DEFAULT_POSECHECK_THRESHOLDS,
    DEFAULT_POSECHECK_NATIVE_BASELINE_LIMITS,
    POSEBUSTERS_VIOLATION_THRESHOLD,
    pose_quota_for_compound_count,
)
from metric_scoring import (
    DEFAULT_SCREENING_MODE,
    LABEL_TO_LEGACY_GATE as _LABEL_TO_LEGACY_GATE,
    NON_MOLECULAR_EFFECTS,
    POSE_EVALUATION_METRIC_RULES,
    RANKING_BLOCKING_EFFECTS,
    aggregate_pose_gate,
    assess_record,
    canonical_posebusters_check_id,
    metric_label,
    molecule_gate_report,
    resolve_screening_mode,
)
from io_utils import number as _number, results_by_pose_index as _results_by_pose_index, load_smiles_map as _load_smiles_map
from pose_scores import (
    LEGACY_SCORE_SPECS,
    POSE_SELECTION_PREFERENCE,
    collect_pose_scores,
    find_score,
    is_rankable_score,
    resolve_primary_metric,
    score_sort_key,
)
from pose_identity import (
    PoseIdentityError,
    artifact_signature,
    build_batch_fingerprint,
    build_batch_output_digest,
    build_pose_identity,
    sdf_molblock_sha256,
    validate_batch_metadata,
    validate_batch_output_digest,
)
from pose_identity import LIGAND_LINEAGE_PROPERTIES
from docking_artifact_contract import load_committed_docking_artifacts
from candidate_topology import build_candidate_topologies
from amber_mmgbsa import run_amber_mmgbsa
from openfe_rbfe import run_openfe_rbfe
from ties_rbfe import plan_ties_rbfe
from rtmscore_rescoring import run_rtmscore_rescoring
from pose_metric_extensions import attach_pose_geometry, merge_md_metrics
from supplemental_pose_interactions import attach_supplemental_interactions
from openmm_short_md import run_openmm_short_md


DEFAULT_OVERALL_TIMEOUT_SECONDS = 86400
DEFAULT_BATCH_SIZE = 50
DEFAULT_CHECKPOINT_INTERVAL_SECONDS = 300
EVALUATION_CANONICAL_ARTIFACTS = (
    "evaluation_report.json",
    "evaluation_report.partial.json",
    "evaluation_summary.json",
    "evaluation_input_poses.sdf",
    "evaluation_input_manifest.csv",
    "selected_poses.sdf",
    "selection_manifest.csv",
)
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


def _scientific_cohort_for(*sources: dict | None) -> dict[str, object]:
    """Resolve row provenance, preferring score/SDF evidence over defaults."""
    components: dict[str, object] = {}
    provenance: dict[str, object] = {}
    explicit_cohort = None
    for source in sources:
        if not isinstance(source, dict):
            continue
        if explicit_cohort in (None, "") and source.get("cohort_id") not in (
            None, "",
        ):
            explicit_cohort = source["cohort_id"]
        for field in SCIENTIFIC_COHORT_COMPONENT_FIELDS:
            if field not in components and source.get(field) not in (None, ""):
                components[field] = source[field]
        for field in SCIENTIFIC_COHORT_PROVENANCE_FIELDS:
            if field not in provenance and source.get(field) not in (None, ""):
                provenance[field] = source[field]
        if (
            "protocol_profile_id" not in components
            and source.get("protocol_id") not in (None, "")
        ):
            components["protocol_profile_id"] = source["protocol_id"]
    if len(components) == len(SCIENTIFIC_COHORT_COMPONENT_FIELDS):
        missing = provenance.get("cohort_identity_missing_fields")
        if not isinstance(missing, list):
            missing = [
                field for field, value in components.items()
                if str(value).startswith(("unknown:", "pending:"))
                or "|state_source:" in str(value)
            ]
        return {
            **components,
            "cohort_id": "|".join(
                f"{field}={components[field]}"
                for field in SCIENTIFIC_COHORT_COMPONENT_FIELDS
            ),
            **provenance,
            "cohort_identity_complete": not missing,
            "cohort_identity_missing_fields": missing,
        }
    # Legacy standalone inputs may only carry their historical opaque label.
    # Keep accepting it, but do not invent comparable scientific components.
    return {"cohort_id": explicit_cohort} if explicit_cohort not in (None, "") else {}

# Evaluator tools that can be toggled via interface-v2 `tools` block.
EVALUATOR_TOOL_DEFAULTS = {
    "identity_qc": True,
    "posebusters": True,
    "posecheck": True,
    "plip": True,
    # PLIP is the primary static-pose Blueprint profiler. ProLIF remains an
    # explicit secondary review route for detector conflicts, MD/multi-pose
    # fingerprints, or Blueprint interactions that PLIP cannot evaluate.
    "prolif": False,
    "blueprint_match": True,
    "pocket_occupancy": True,
    "arpeggio": False,
    "luna": False,
    "oddt": False,
    "rdkit_pose": False,
    "freesasa": False,
    "rtmscore": False,
    "rescore_rtmscore": False,
    # Expensive candidate-specific parameterization is an explicit escalation,
    # never part of the default pose-QC protocol.
    "candidate_topology": False,
    "openmm_md": False,
    "mdanalysis": False,
    # Endpoint rescoring consumes only the auditable candidate topology bundle.
    # It is never enabled by routine pose QC.
    "mmgbsa": False,
    # RBFE backends are contextual planning escalations. They never run as part
    # of routine pose QC and never manufacture a ΔΔG from a plan.
    "openfe_rbfe": False,
    "ties_rbfe": False,
}
EVALUATOR_TOOLS = tuple(EVALUATOR_TOOL_DEFAULTS)


def _scoring_plan_tools(
    switches: dict,
    *,
    blueprint_expected: bool,
    pocket_expected: bool,
) -> list[str]:
    """Translate evaluator switches into scorer-facing planned tool IDs."""
    enabled = {
        tool_id for tool_id in (
            "posebusters", "posecheck", "plip", "prolif",
            "rdkit_pose", "freesasa", "arpeggio", "mdanalysis", "rtmscore",
        )
        if switches.get(tool_id)
    }
    if blueprint_expected:
        enabled.add("blueprint_match")
    if pocket_expected:
        enabled.add("pocket_occupancy")
    if switches.get("mmgbsa"):
        enabled.add("amber_mmpbsa")
    return sorted(enabled)

_POSE_METRIC_RULE_BY_ID = {
    rule["metric_id"]: rule for rule in POSE_EVALUATION_METRIC_RULES
}
_DOCKING_METRIC_RULES = tuple(
    rule for rule in POSE_EVALUATION_METRIC_RULES
    if rule.get("stage") == "docking"
)
_DOCKING_GATE_ORDER = {
    "keep": 0,
    "modify": 1,
    # An unmeasured gate is not a worse gate than a rejected one.  This order
    # drives the pre-QC ``best_per_compound`` choice and the ``max_poses``
    # truncation, so ranking "not_evaluated" last used to discard a pose for a
    # missing score before any evaluator had looked at it -- the "missing as the
    # worst value" pattern the executable policy rules out.
    "not_evaluated": 2,
    "reject": 3,
}
# A pose whose evidence never resolved is not a rejected pose. The executable policy asks
# for ``Score=NA`` plus a coverage/confidence note, so "not_evaluated" stays in
# the survivor cohort -- ordered behind every measured pose -- and only a real
# hard-gate failure leaves.
_POSE_SURVIVOR_VERDICTS = frozenset({"pass", "warn", "not_evaluated"})
_EVALUATED_VERDICT_ORDER = {"pass": 0, "warn": 1, "not_evaluated": 2}
_POSEBUSTERS_CHECK_RULES = tuple(
    rule
    for rule in POSE_EVALUATION_METRIC_RULES
    if (rule.get("extractor") or {}).get("kind") == "posebusters_check"
)
_POSEBUSTERS_CHECK_RULE_BY_CHECK_ID = {
    canonical_posebusters_check_id((rule.get("extractor") or {}).get("check")): rule
    for rule in _POSEBUSTERS_CHECK_RULES
}


def _tool_switches(input_data: dict) -> dict:
    """Resolve {tool: enabled} for evaluator tools.

    `tools` (interface v2) maps tool -> {"enabled": bool}. When the block is
    absent, core pose-QC tools default on while expensive topology, MM-GBSA,
    and RBFE escalations default off. A disabled tool is skipped without
    penalty; it can never lift a pose to `pass` because the evidence invariant
    in _compute_verdict still holds.
    """
    from io_utils import tool_enabled
    return {
        tool_id: tool_enabled(input_data, tool_id, default=default)
        for tool_id, default in EVALUATOR_TOOL_DEFAULTS.items()
    }


def _build_evaluator_batch_fingerprint(
    *,
    input_data: dict,
    receptor_path: Path,
    batch_sdf: Path,
    batch_poses: list[dict],
    posecheck_thresholds: dict,
    posecheck_calibration: dict,
    screening_mode: dict,
    tool_switches: dict,
    run_cache: dict | None = None,
) -> dict:
    """Fingerprint every input that can change a resumable batch result."""
    raw_tools = input_data.get("tools") or {}
    tool_contracts = {}
    for tool_id in EVALUATOR_TOOLS:
        raw = raw_tools.get(tool_id) if isinstance(raw_tools, dict) else None
        raw = raw if isinstance(raw, dict) else {}
        contract = {
            "enabled": bool(tool_switches.get(tool_id, False)),
            "params": dict(raw.get("params") or {}),
        }
        if contract["enabled"]:
            contract["runtime"] = dict(raw.get("runtime") or {})
            if tool_id == "posecheck":
                # Old batches rebuilt receptor hydrogens with Hydride.
                contract["receptor_loading"] = "prepared_pdb_geometry_v1"
        tool_contracts[tool_id] = contract

    artifact_inputs = {}
    for field in (
        "interaction_blueprint_file",
        "pocket_features_file",
        "reference_ligand_file",
        "docking_scores_file",
        "docking_report_file",
        "library_manifest_file",
    ):
        value = input_data.get(field)
        if not value:
            continue
        path = Path(str(value))
        artifact_inputs[field] = (
            _cached_signature(path, label=field, run_cache=run_cache)
            if path.is_file()
            else {"status": "missing", "path": str(path)}
        )

    runtime_metadata = {
        tool_id: contract.get("runtime") or {}
        for tool_id, contract in tool_contracts.items()
        if contract["enabled"]
    }
    model_metadata = _collect_named_metadata(
        tool_contracts,
        ("model", "model_id", "cnn_model", "scorer"),
    )
    checkpoint_metadata = _collect_named_metadata(
        tool_contracts,
        ("checkpoint", "weights", "model_path", "checkpoint_path"),
        hash_files=True,
        run_cache=run_cache,
    )
    return build_batch_fingerprint(
        receptor_artifact=_cached_signature(receptor_path, label="receptor_artifact", run_cache=run_cache),
        pose_artifacts=batch_sdf,
        enabled_tools=tool_contracts,
        params={
            "chain": input_data.get("chain", "A"),
            "docking_attempt": _docking_attempt_number(input_data),
            "posecheck_thresholds": posecheck_thresholds,
            "posecheck_calibration": posecheck_calibration,
            "posecheck_num_conformers": _posecheck_num_conformers(input_data),
            # The mode decides which per-category fail count/ratio becomes a
            # rejection, so it changes every ``metric_gate``/``verdict`` in the
            # batch.  Without it here, resuming a V2 batch under V3 replayed the
            # cached rows verbatim and silently kept the old thresholds.  The
            # full spec is hashed, not just ``mode_id``, so a retuned threshold
            # invalidates too.
            "screening_mode": {
                key: screening_mode.get(key)
                for key in (
                    "mode_id",
                    "eliminates",
                    "fail_count_threshold",
                    "fail_ratio_threshold",
                )
            },
        },
        runtime_metadata=runtime_metadata,
        model_metadata=model_metadata,
        checkpoint_metadata=checkpoint_metadata,
        version_metadata=_enabled_tool_versions(tool_switches),
        extra_metadata={
            "input_artifacts": artifact_inputs,
            "batch_pose_contract": [
                {
                    "source_pose_index": pose.get("source_pose_index"),
                    "pose_id": pose.get("pose_id"),
                    "chemical_identity_id": (
                        (pose.get("chemical_identity") or {}).get("identity_id")
                    ),
                    "docking_scores": pose.get("docking_scores") or {},
                }
                for pose in batch_poses
            ],
            "implementation_files": _implementation_signatures(tool_switches, run_cache=run_cache),
        },
    )


_BATCH_POSE_BINDING_FIELDS = (
    "mol_name",
    "parent_pose_id",
    "coordinate_revision",
    "lineage_metadata",
    "chemical_identity",
    "identity_status",
    "sdf_molblock_sha256",
    "compound_id",
    "source_molecule_index",
    "preparation_variant_index",
    "engine_pose_index",
    "engine_pose_rank",
    "selection_rank",
    "compound_pose_rank",
    "selection_strategy",
    "pre_qc_docking_gate",
    "primary_score_metric_id",
    "docking_scores",
    "pose_scores",
)


def _strict_json_equal(left: object, right: object) -> bool:
    """Compare JSON data without Python's bool/int equality ambiguity."""
    try:
        return json.dumps(
            left,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ) == json.dumps(
            right,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError):
        return False


def _resume_batch_failure(status: str, reason: str, **details) -> dict:
    return {
        "reusable": False,
        "status": status,
        "reason": reason,
        **details,
    }


def _validate_resumable_batch_evaluation(
    batch_report: dict,
    batch_poses: list[dict],
) -> dict:
    """Bind every cached evaluation row to exactly one current batch pose."""
    entries = batch_report.get("evaluation")
    if not isinstance(entries, list):
        return _resume_batch_failure(
            "invalid_batch_evaluation",
            "cached batch evaluation must be a list",
        )

    output_validation = validate_batch_output_digest(batch_report, entries)
    if not output_validation.get("reusable"):
        return output_validation
    if len(entries) != len(batch_poses):
        return _resume_batch_failure(
            "batch_evaluation_count_mismatch",
            "cached evaluation row count differs from the current batch",
            expected_count=len(batch_poses),
            actual_count=len(entries),
        )

    expected_source_indices: set[int] = set()
    expected_pose_ids: set[str] = set()
    observed_source_indices: set[int] = set()
    observed_pose_ids: set[str] = set()
    for local_index, (entry, pose) in enumerate(zip(entries, batch_poses)):
        if not isinstance(entry, dict):
            return _resume_batch_failure(
                "invalid_batch_evaluation_row",
                f"cached evaluation row {local_index} is not an object",
            )
        try:
            expected_source = _strict_integer(pose.get("source_pose_index"))
            observed_source = _strict_integer(entry.get("source_pose_index"))
            observed_pose_index = _strict_integer(entry.get("pose_index"))
            observed_batch_index = _strict_integer(entry.get("batch_pose_index"))
        except ValueError:
            return _resume_batch_failure(
                "invalid_batch_evaluation_index",
                f"cached evaluation row {local_index} has a missing or invalid index",
            )
        expected_pose_id = pose.get("pose_id")
        observed_pose_id = entry.get("pose_id")
        if not isinstance(expected_pose_id, str) or not expected_pose_id.strip():
            return _resume_batch_failure(
                "invalid_current_batch_pose_identity",
                f"current batch pose {local_index} has no canonical pose_id",
            )
        if not isinstance(observed_pose_id, str) or not observed_pose_id.strip():
            return _resume_batch_failure(
                "missing_cached_pose_identity",
                f"cached evaluation row {local_index} has no canonical pose_id",
            )
        if expected_source in expected_source_indices:
            return _resume_batch_failure(
                "duplicate_current_batch_source_pose_index",
                f"current batch repeats source_pose_index {expected_source}",
            )
        if expected_pose_id in expected_pose_ids:
            return _resume_batch_failure(
                "duplicate_current_batch_pose_id",
                f"current batch repeats pose_id {expected_pose_id}",
            )
        if observed_source in observed_source_indices:
            return _resume_batch_failure(
                "duplicate_cached_source_pose_index",
                f"cached evaluation repeats source_pose_index {observed_source}",
            )
        if observed_pose_id in observed_pose_ids:
            return _resume_batch_failure(
                "duplicate_cached_pose_id",
                f"cached evaluation repeats pose_id {observed_pose_id}",
            )
        expected_source_indices.add(expected_source)
        expected_pose_ids.add(expected_pose_id)
        observed_source_indices.add(observed_source)
        observed_pose_ids.add(observed_pose_id)

        if (
            observed_batch_index != local_index
            or observed_source != expected_source
            or observed_pose_index != expected_source
            or observed_pose_id != expected_pose_id
        ):
            return _resume_batch_failure(
                "batch_evaluation_pose_mismatch",
                f"cached evaluation row {local_index} is bound to another pose",
                expected_source_pose_index=expected_source,
                observed_source_pose_index=observed_source,
                expected_pose_id=expected_pose_id,
                observed_pose_id=observed_pose_id,
            )

        # These values are copied from the current pose, not produced by an
        # evaluator.  Compare their exact JSON representations so a cache row
        # cannot borrow chemical/source identity, scores, or selection context.
        for field in _BATCH_POSE_BINDING_FIELDS:
            if not _strict_json_equal(entry.get(field), pose.get(field)):
                return _resume_batch_failure(
                    "batch_evaluation_lineage_mismatch",
                    f"cached evaluation row {local_index} disagrees on {field}",
                    row_index=local_index,
                    field=field,
                )

        # State and cohort fields use symmetric presence checks: no lineage
        # field may be silently dropped, added, or replaced during resume.
        for field in (
            *LIGAND_LINEAGE_PROPERTIES,
            *SCIENTIFIC_COHORT_OUTPUT_FIELDS,
        ):
            expected_present = pose.get(field) not in (None, "")
            observed_present = entry.get(field) not in (None, "")
            if expected_present != observed_present or (
                expected_present
                and not _strict_json_equal(entry[field], pose[field])
            ):
                return _resume_batch_failure(
                    "batch_evaluation_lineage_mismatch",
                    f"cached evaluation row {local_index} disagrees on {field}",
                    row_index=local_index,
                    field=field,
                )

    if observed_source_indices != expected_source_indices:
        return _resume_batch_failure(
            "batch_evaluation_source_set_mismatch",
            "cached evaluation does not cover the current source pose set",
        )
    if observed_pose_ids != expected_pose_ids:
        return _resume_batch_failure(
            "batch_evaluation_pose_set_mismatch",
            "cached evaluation does not cover the current pose identity set",
        )
    return {
        "reusable": True,
        "status": "match",
        "input_binding": "exact",
        "row_count": len(entries),
        "output_digest": output_validation,
    }


def _cached_signature(path: Path, *, label: str, run_cache: dict | None = None) -> dict:
    return cached_file(run_cache, Path(path), "sha256", lambda item: artifact_signature(item, label=label))


def _collect_named_metadata(
    value,
    key_fragments: tuple[str, ...],
    *,
    hash_files: bool = False,
    run_cache: dict | None = None,
) -> dict:
    """Collect model/checkpoint declarations while preserving their path."""
    collected = {}

    def visit(item, path: str) -> None:
        if isinstance(item, dict):
            for key, child in item.items():
                child_path = f"{path}.{key}" if path else str(key)
                if any(fragment in str(key).lower() for fragment in key_fragments):
                    payload = child
                    if hash_files and isinstance(child, str) and Path(child).is_file():
                        payload = {
                            "path": child,
                            "content": _cached_signature(Path(child), label=child_path, run_cache=run_cache),
                        }
                    collected[child_path] = payload
                visit(child, child_path)
        elif isinstance(item, list):
            for index, child in enumerate(item):
                visit(child, f"{path}[{index}]")

    visit(value, "")
    return dict(sorted(collected.items()))


def _enabled_tool_versions(tool_switches: dict) -> dict:
    package_names = {
        "identity_qc": "rdkit",
        "posebusters": "posebusters",
        "posecheck": "posecheck",
        "plip": "plip",
        "prolif": "prolif",
    }
    versions = {}
    for tool_id, package in package_names.items():
        if not tool_switches.get(tool_id, False):
            continue
        try:
            versions[tool_id] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[tool_id] = "unavailable"
    return versions


def _implementation_signatures(tool_switches: dict, *, run_cache: dict | None = None) -> dict:
    paths = {
        "evaluation_runtime": SKILL_DIR / "scripts" / "evaluation_runtime.py",
        "structure_evaluator": Path(__file__),
        "pose_identity": SHARED_DIR / "pose_identity.py",
        "pose_scores": SHARED_DIR / "pose_scores.py",
        "metric_scoring": SHARED_DIR / "metric_scoring.py",
        "constants": SHARED_DIR / "constants.py",
        "posebusters": SKILL_DIR / "scripts" / "posebusters_utils.py",
        "posecheck": SKILL_DIR / "scripts" / "posecheck_utils.py",
        "plip": SKILL_DIR / "scripts" / "plip_eval_utils.py",
        "prolif": SKILL_DIR / "scripts" / "prolif_eval_utils.py",
        "prolif_fingerprint": SKILL_DIR.parent / "pocket-analysis/scripts/prolif_utils.py",
        "blueprint_match": SKILL_DIR / "scripts" / "blueprint_match.py",
        "pocket_occupancy": SKILL_DIR / "scripts" / "pocket_occupancy.py",
        "candidate_topology": SKILL_DIR / "scripts" / "candidate_topology.py",
        "pose_metric_extensions": SKILL_DIR / "scripts" / "pose_metric_extensions.py",
        "supplemental_pose_interactions": SKILL_DIR / "scripts" / "supplemental_pose_interactions.py",
        "arpeggio": SKILL_DIR.parent / "pocket-analysis/scripts/arpeggio_utils.py",
        "luna": SKILL_DIR.parent / "pocket-analysis/scripts/luna_utils.py",
        "oddt": SKILL_DIR.parent / "pocket-analysis/scripts/oddt_utils.py",
        "supplemental_native_worker": SKILL_DIR.parent / "pocket-analysis/scripts/native_worker.py",
        "rtmscore": SKILL_DIR / "scripts" / "rtmscore_rescoring.py",
        "rtmscore_worker": SKILL_DIR / "scripts" / "rtmscore_worker.py",
        "openmm_md": SKILL_DIR / "scripts" / "openmm_short_md.py",
        "openmm_md_worker": SKILL_DIR / "scripts" / "openmm_md_worker.py",
        "mmgbsa": SKILL_DIR / "scripts" / "amber_mmgbsa.py",
        "openfe_rbfe": SKILL_DIR / "scripts" / "openfe_rbfe.py",
        "ties_rbfe": SKILL_DIR / "scripts" / "ties_rbfe.py",
    }
    signatures = {}
    for name, path in paths.items():
        if name in tool_switches and not tool_switches[name]:
            continue
        if name == "rtmscore_worker" and not tool_switches.get(
            "rtmscore", False
        ):
            continue
        if name == "openmm_md_worker" and not tool_switches.get(
            "openmm_md", False
        ):
            continue
        if path.is_file():
            signatures[name] = _cached_signature(path, label=f"implementation.{name}", run_cache=run_cache)
    return signatures


def _run_candidate_topology(
    *,
    input_data: dict,
    receptor_path: Path,
    selected_sdf: Path,
    output_dir: Path,
) -> dict:
    """Build the optional Amber topology bundle for post-QC survivors only."""
    config = ((input_data.get("tools") or {}).get("candidate_topology") or {})
    params = dict(config.get("params") or {})
    runtime = dict(config.get("runtime") or {})
    openmm_config = ((input_data.get("tools") or {}).get("openmm_md") or {})
    openmm_enabled = bool(openmm_config.get("enabled", False))
    openmm_params = dict(openmm_config.get("params") or {})
    return build_candidate_topologies(
        receptor_path,
        selected_sdf,
        output_dir / "candidate_topologies",
        runtime,
        execution_mode=params.get("execution_mode", "build"),
        protein_force_field=params.get("protein_force_field", "ff14SB"),
        ligand_force_field=params.get("ligand_force_field", "gaff2"),
        ligand_charge_method=params.get("ligand_charge_method", "bcc"),
        max_candidates=params.get("max_candidates", 20),
        timeout_per_candidate_seconds=params.get(
            "timeout_per_candidate_seconds", 2400
        ),
        is_covalent=bool(input_data.get("is_covalent", False)),
        has_metal=bool(input_data.get("has_metal", False)),
        allow_receptor_hetero=bool(params.get("allow_receptor_hetero", False)),
        fail_fast=bool(params.get("fail_fast", True)),
        solvate_for_md=openmm_enabled,
        water_model=openmm_params.get("water_model", "tip3p"),
        solvent_padding_angstrom=openmm_params.get(
            "solvent_padding_angstrom", 10.0
        ),
    )


def _run_rtmscore(
    *,
    input_data: dict,
    receptor_path: Path,
    selected_sdf: Path,
    output_dir: Path,
) -> dict:
    config = ((input_data.get("tools") or {}).get("rtmscore") or {})
    params = dict(config.get("params") or {})
    runtime = dict(config.get("runtime") or {})
    return run_rtmscore_rescoring(
        receptor_path,
        selected_sdf,
        output_dir / "rtmscore",
        runtime,
        model_ensemble=params.get("model_ensemble", [1, 2, 3]),
        pocket_cutoff_angstrom=params.get("pocket_cutoff_angstrom", 10.0),
        batch_size=params.get("batch_size", 128),
        num_workers=params.get("num_workers", 0),
        parallel_graph_build=bool(params.get("parallel_graph_build", False)),
        device=params.get("device", "auto"),
        timeout_seconds=params.get("timeout_seconds", 86400),
    )


def _run_openmm_md(
    *,
    input_data: dict,
    candidate_topologies_manifest: str | Path,
    output_dir: Path,
) -> dict:
    openmm_config = ((input_data.get("tools") or {}).get("openmm_md") or {})
    mda_config = ((input_data.get("tools") or {}).get("mdanalysis") or {})
    params = dict(openmm_config.get("params") or {})
    mda_params = dict(mda_config.get("params") or {})
    return run_openmm_short_md(
        candidate_topologies_manifest,
        output_dir / "short_md",
        dict(openmm_config.get("runtime") or {}),
        mdanalysis_runtime=dict(mda_config.get("runtime") or {}),
        run_trajectory_analysis=bool(mda_config.get("enabled", False)),
        top_k=params.get("top_k", 3),
        platform=params.get("platform", "CUDA"),
        precision=params.get("precision", "mixed"),
        device_index=params.get("device_index", 0),
        temperature_kelvin=params.get("temperature_kelvin", 300.0),
        pressure_bar=params.get("pressure_bar", 1.0),
        timestep_femtoseconds=params.get("timestep_femtoseconds", 2.0),
        friction_per_picosecond=params.get("friction_per_picosecond", 1.0),
        minimization_max_iterations=params.get(
            "minimization_max_iterations", 500
        ),
        equilibration_picoseconds=params.get("equilibration_picoseconds", 20.0),
        production_picoseconds=params.get("production_picoseconds", 100.0),
        report_interval_picoseconds=params.get(
            "report_interval_picoseconds", 2.0
        ),
        seed=params.get("seed", 42),
        timeout_per_candidate_seconds=params.get(
            "timeout_per_candidate_seconds", 86400
        ),
        fail_fast=bool(params.get("fail_fast", False)),
        alignment_selection=mda_params.get(
            "alignment_selection", "protein and backbone"
        ),
        ligand_selection=mda_params.get(
            "ligand_selection", "resname LIG and not name H*"
        ),
        contact_cutoff_angstrom=mda_params.get(
            "contact_cutoff_angstrom", 4.5
        ),
        analysis_stride=mda_params.get("analysis_stride", 1),
        analysis_timeout_per_candidate_seconds=mda_params.get(
            "timeout_per_candidate_seconds", 3600
        ),
    )


def _run_mmgbsa(
    *,
    input_data: dict,
    candidate_topologies_manifest: str | Path,
    output_dir: Path,
) -> dict:
    """Run optional single-frame Amber MM-GBSA on built post-QC topologies."""
    config = ((input_data.get("tools") or {}).get("mmgbsa") or {})
    params = dict(config.get("params") or {})
    runtime = dict(config.get("runtime") or {})
    return run_amber_mmgbsa(
        candidate_topologies_manifest,
        output_dir / "mmgbsa",
        runtime,
        execution_mode=params.get("execution_mode", "preflight"),
        top_k=params.get("top_k", 20),
        gb_model=params.get("gb_model", 5),
        salt_concentration_molar=params.get(
            "salt_concentration_molar", 0.15
        ),
        timeout_per_candidate_seconds=params.get(
            "timeout_per_candidate_seconds", 900
        ),
        fail_fast=bool(params.get("fail_fast", False)),
        keep_intermediate_files=bool(
            params.get("keep_intermediate_files", False)
        ),
        use_as_primary_metric=bool(
            params.get("use_as_primary_metric", True)
        ),
    )


def _run_ties_rbfe(
    *,
    input_data: dict,
    receptor_path: Path,
    selected_sdf: Path,
    output_dir: Path,
) -> dict:
    """Write the optional planning-only TIES campaign manifest."""
    config = ((input_data.get("tools") or {}).get("ties_rbfe") or {})
    params = dict(config.get("params") or {})
    runtime = dict(config.get("runtime") or {})
    return plan_ties_rbfe(
        receptor_path,
        selected_sdf,
        output_dir / "ties_rbfe",
        runtime,
        execution_mode=params.get("execution_mode", "preflight"),
        network_strategy=params.get("network_strategy", "star"),
        reference_pose_id=params.get("reference_pose_id"),
        max_candidates=params.get("max_candidates", 8),
        is_covalent=bool(input_data.get("is_covalent", False)),
        has_metal=bool(input_data.get("has_metal", False)),
        require_post_qc_lineage=bool(
            params.get("require_post_qc_lineage", True)
        ),
        runtime_probe_timeout_seconds=params.get(
            "runtime_probe_timeout_seconds", 30
        ),
        simulation_platform=params.get("simulation_platform", "auto"),
    )


def _run_openfe_rbfe(
    *,
    input_data: dict,
    receptor_path: Path,
    selected_sdf: Path,
    output_dir: Path,
) -> dict:
    """Run OpenFE applicability checks or create a mapping-only network."""
    config = ((input_data.get("tools") or {}).get("openfe_rbfe") or {})
    params = dict(config.get("params") or {})
    runtime = dict(config.get("runtime") or {})
    return run_openfe_rbfe(
        receptor_path,
        selected_sdf,
        output_dir / "openfe_rbfe",
        runtime,
        execution_mode=params.get("execution_mode", "preflight"),
        network_method=params.get("network_method", "minimal_spanning"),
        mapper=params.get("mapper", "KartografAtomMapper"),
        mapper_settings=params.get("mapper_settings"),
        max_candidates=params.get("max_candidates", 8),
        min_ligands=params.get("min_ligands", 2),
        allow_charge_changes=bool(params.get("allow_charge_changes", False)),
        has_metal=bool(input_data.get("has_metal", False)),
        is_covalent=bool(input_data.get("is_covalent", False)),
        timeout_seconds=params.get("timeout_seconds", 600),
        n_processes=params.get("n_processes", 1),
        mst_num=params.get("mst_num", 2),
        central_pose_id=params.get("central_pose_id"),
        min_mapping_score=params.get("min_mapping_score", 0.3),
        min_mapped_heavy_fraction=params.get(
            "min_mapped_heavy_fraction", 0.5
        ),
        require_explicit_hydrogens=bool(
            params.get("require_explicit_hydrogens", True)
        ),
        allow_undefined_stereochemistry=bool(
            params.get("allow_undefined_stereochemistry", False)
        ),
    )


def _receptor_identity(input_data: dict, receptor_path: Path) -> dict:
    """Name the protein this cohort was docked against.

    ``source_accession`` is the first column of the published final report, so it
    needs a producer.  An explicit accession wins; otherwise the receptor file's
    own stem is the only honest answer, and ``source_type`` says which of the two
    the reader is looking at.
    """
    accession = str(
        input_data.get("source_accession")
        or input_data.get("pdb_id")
        or ""
    ).strip()
    return {
        "receptor_file": str(receptor_path),
        "source_accession": accession or receptor_path.stem,
        "source_type": "accession" if accession else "file",
    }


def _invalidate_previous_evaluation_artifacts(output_dir: Path) -> None:
    """Remove prior commit markers/flat outputs before a same-dir rerun."""
    for name in EVALUATION_CANONICAL_ARTIFACTS:
        (output_dir / name).unlink(missing_ok=True)


def run(input_data: dict, output_dir: Path) -> dict:
    """Run structure evaluation pipeline on docked poses."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Invalidate the previous flat result before *any* validation return.  A
    # missing input on a rerun must not leave the earlier evaluation report
    # looking like the outcome of the new attempt.
    _invalidate_previous_evaluation_artifacts(output_dir)
    receptor_file = input_data.get("receptor_file")
    docked_sdf = input_data.get("docked_sdf")
    blueprint_file = input_data.get("interaction_blueprint_file")
    pocket_features_file = input_data.get("pocket_features_file")
    reference_ligand_file = input_data.get("reference_ligand_file")
    docking_scores_file = input_data.get("docking_scores_file")
    docking_report_file = input_data.get("docking_report_file")
    library_manifest_file = input_data.get("library_manifest_file")
    top_n = int(input_data.get("top_n", 10))

    if not receptor_file or not docked_sdf:
        return {"success": False, "error": "receptor_file and docked_sdf are required"}

    receptor_path = Path(receptor_file)
    docked_path = Path(docked_sdf)

    if not receptor_path.exists():
        return {"success": False, "error": f"receptor_file not found: {receptor_path}"}
    if not docked_path.exists():
        return {"success": False, "error": f"docked_sdf not found: {docked_path}"}

    # Resolved once: every exit path -- timeout, per-batch, and final -- publishes
    # the same receptor, so the partial reports name it too.
    scientific_cohort = _scientific_cohort_for(input_data)
    receptor_identity = {
        **_receptor_identity(input_data, receptor_path),
        **scientific_cohort,
    }
    posecheck_thresholds = _load_posecheck_thresholds(input_data)
    # Level 2 of the gate contract: which per-category fail count/ratio turns a
    # set of level-1 marks into a rejection.  ``candidate_source``
    # only supplies the default when no mode is named explicitly.
    screening_mode = resolve_screening_mode(
        input_data.get("screening_mode"),
        context=input_data.get("candidate_source"),
    )
    posecheck_baseline_limits = _load_posecheck_baseline_limits(input_data)
    calibrate_posecheck = bool(
        input_data.get("calibrate_posecheck_against_reference", False)
    )
    posecheck_calibration = {
        "mode": "absolute_thresholds",
        "requested": calibrate_posecheck,
        "reference_ligand_file": reference_ligand_file,
        "native_baseline": None,
    }

    docking_artifact_lineage = None
    try:
        if docking_report_file not in (None, ""):
            (
                docking_artifact_lineage,
                docking_scores,
                docking_scores_file,
            ) = _load_verified_docking_inputs(
                docking_report_file,
                docked_path=docked_path,
                docking_scores_file=docking_scores_file,
                receptor_path=receptor_path,
            )
        else:
            docking_scores = _load_docking_scores_contract(
                docking_scores_file
            )
    except PoseIdentityError as exc:
        return {
            "success": False,
            "error": str(exc),
            "operational_errors": [{
                "tool": "structure-evaluator",
                "step": "docking_scores_artifact_contract",
                "code": exc.code,
                "scope": "run",
                "retryable": True,
                "reason": str(exc),
                **({"details": exc.details} if exc.details else {}),
            }],
        }
    smiles_map = _load_smiles_map(library_manifest_file)
    selection_strategy = str(input_data.get("selection_strategy") or "ranked_per_compound").lower()
    sequential_selection = selection_strategy == "first_passing_per_compound"
    budgeted_selection = selection_strategy in {"ranked_per_compound", "first_passing_per_compound"}
    poses_per_compound = max(1, int(input_data.get("poses_per_compound", 20)))
    # Per-tool switches (interface v2). Core QC defaults on; expensive topology,
    # MM-GBSA, and RBFE escalations default off. A disabled tool is
    # "disabled_by_config": skipped without penalty, but it cannot lift a pose
    # to pass (see _compute_verdict evidence invariant).
    tool_switches = _tool_switches(input_data)
    if sequential_selection and any(tool_switches[name] for name in ("mdanalysis", "mmgbsa")):
        return {
            "success": False,
            "error": "first_passing_per_compound supports initial pose QC and RTMScore; "
                     "use selection_strategy=all when enabling mdanalysis or mmgbsa post-QC gates",
        }
    if not tool_switches["identity_qc"]:
        return {
            "success": False,
            "error": (
                "identity_qc cannot be disabled: pose evidence, scores, and "
                "resume metadata require canonical pose identity"
            ),
        }
    # Every pose identity is built by RDKit and ``identity_qc`` cannot be switched
    # off, so an absent RDKit leaves the run with no way to name a pose. That is a
    # fault of this environment, not a finding about any molecule: the previous
    # behaviour stamped ``rdkit_unavailable`` on every record, which reached
    # ``evidence_status["identity_qc"] == "failed"`` and rejected the entire
    # cohort on chemistry grounds. The executable policy keeps a tool fault out of the
    # molecule label, so declare it once against the run and stop.
    if not _rdkit_available():
        return {
            "success": False,
            "error": (
                "RDKit is required for stable pose identity and is not "
                "importable in this environment"
            ),
            "operational_errors": [{
                "tool": "rdkit",
                "step": "pose_identity",
                "code": "rdkit_unavailable",
                "scope": "run",
                "retryable": True,
                "reason": (
                    "RDKit is required for stable pose identity; install it "
                    "and rerun structure-evaluator"
                ),
            }],
        }
    max_poses = input_data.get("max_poses")
    max_poses = int(max_poses) if max_poses is not None else None
    batch_size = int(input_data.get("batch_size", DEFAULT_BATCH_SIZE))
    batch_size = max(1, batch_size)
    # Input caches belong to this invocation; never reuse across experiments.
    run_cache: dict = {}
    checkpoint_interval = max(0.0, float(input_data.get(
        "checkpoint_interval_seconds", DEFAULT_CHECKPOINT_INTERVAL_SECONDS,
    )))
    profiler_worker_count(input_data.get("profiler_workers", 2))
    last_checkpoint_at = None
    resume = bool(input_data.get("resume", True))
    overall_timeout_seconds = int(
        input_data.get("overall_timeout_seconds", DEFAULT_OVERALL_TIMEOUT_SECONDS)
    )
    started_at = monotonic()

    try:
        pose_records = _read_pose_records(docked_path, run_cache=run_cache)
        for record in pose_records:
            record.update(_scientific_cohort_for(record, scientific_cohort))
            record["docking_attempt"] = _docking_attempt_number(input_data, record)
        evaluation_input_poses = _select_pose_records_from_metadata(
            pose_records,
            docking_scores,
            selection_strategy=selection_strategy,
            poses_per_compound=poses_per_compound,
            max_poses=None if budgeted_selection else max_poses,
            screening_mode=screening_mode,
        )
    except PoseIdentityError as exc:
        operational_error = {
            "tool": "structure-evaluator",
            "step": "docked_pose_artifact_contract",
            "code": exc.code,
            "scope": "run",
            "retryable": True,
            "reason": str(exc),
        }
        if exc.details:
            operational_error["details"] = exc.details
        return {
            "success": False,
            "error": str(exc),
            "operational_errors": [operational_error],
        }
    if not evaluation_input_poses:
        return {"success": False, "error": "no poses selected for evaluation"}

    evaluation_input_sdf = output_dir / "evaluation_input_poses.sdf"
    try:
        _write_selected_sdf(
            docked_path, evaluation_input_poses, evaluation_input_sdf,
            run_cache=run_cache,
        )
    except PoseIdentityError as exc:
        operational_error = {
            "tool": "structure-evaluator",
            "step": "pose_export_identity",
            "code": exc.code,
            "scope": "run",
            "retryable": True,
            "reason": str(exc),
        }
        if exc.details:
            operational_error["details"] = exc.details
        return {
            "success": False,
            "error": str(exc),
            "operational_errors": [operational_error],
        }
    evaluation_input_manifest = output_dir / "evaluation_input_manifest.csv"
    _write_selection_manifest(evaluation_input_manifest, evaluation_input_poses)

    # ``selected_poses.sdf`` is a post-QC artifact. It is refreshed with each
    # cumulative checkpoint and contains at most one eligible pose per compound.
    selected_sdf = output_dir / "selected_poses.sdf"
    selection_manifest = output_dir / "selection_manifest.csv"
    _write_selected_sdf(docked_path, [], selected_sdf, run_cache=run_cache)
    _write_selection_manifest(selection_manifest, [])

    # Calibrate once for the whole Step 4 run. Batch workers reuse the same
    # thresholds, which keeps resumed runs consistent and avoids repeated native
    # ligand PoseCheck calls.
    global_completed = []
    global_skipped = []
    if (
        calibrate_posecheck
        and tool_switches["posecheck"]
        and check_posecheck_available()
        and reference_ligand_file
    ):
        posecheck_thresholds, posecheck_calibration = _calibrate_posecheck_thresholds(
            receptor_path,
            reference_ligand_file,
            posecheck_thresholds,
            posecheck_baseline_limits,
            global_completed,
            global_skipped,
            num_conformers=_posecheck_num_conformers(input_data),
        )
    elif calibrate_posecheck and tool_switches["posecheck"] and reference_ligand_file:
        global_skipped.append({
            "step": "posecheck_native_calibration",
            "reason": "PoseCheck unavailable",
        })
    elif calibrate_posecheck and tool_switches["posecheck"]:
        global_skipped.append({
            "step": "posecheck_native_calibration",
            "reason": "reference_ligand_file is required for calibration",
        })

    all_evaluation = []
    batch_reports = []
    profiler_operational_errors: list[dict] = []
    batches_dir = output_dir / "batches"
    batches_dir.mkdir(parents=True, exist_ok=True)

    for batch_number, batch_poses in enumerate(_evaluation_batches(
        evaluation_input_poses, all_evaluation, batch_size,
        sequential=sequential_selection,
        dynamic_pending=selection_strategy == "ranked_per_compound",
        resume_batch_sizes=(
            _cached_batch_sizes(batches_dir)
            if resume and selection_strategy == "ranked_per_compound" else ()
        ),
    ), start=1):
        elapsed = monotonic() - started_at
        stop_reason = None
        if elapsed >= overall_timeout_seconds:
            stop_reason = f"structure evaluation timed out after {overall_timeout_seconds}s"
        if budgeted_selection and max_poses is not None:
            remaining = max_poses - len(all_evaluation)
            if remaining <= 0:
                stop_reason = f"structure evaluation reached max_poses={max_poses} before resolving all compounds"
            else:
                batch_poses = batch_poses[:remaining]
        if stop_reason:
            post_qc_selected = _refresh_post_qc_selection(
                docked_path,
                evaluation_input_poses,
                all_evaluation,
                selected_sdf,
                selection_manifest,
                run_cache=run_cache,
            )
            partial = _write_final_outputs(
                output_dir=output_dir,
                evaluation=all_evaluation,
                evaluation_input_poses=evaluation_input_poses,
                selected_poses=post_qc_selected,
                posecheck_thresholds=posecheck_thresholds,
                posecheck_calibration=posecheck_calibration,
                posecheck_num_conformers=_posecheck_num_conformers(input_data),
                reference_ligand_file=reference_ligand_file,
                docking_scores_file=docking_scores_file,
                docking_artifact_lineage=docking_artifact_lineage,
                receptor_identity=receptor_identity,
                steps_completed=global_completed,
                steps_skipped=global_skipped + [{
                    "step": "structure_evaluation_batches",
                    "reason": stop_reason,
                }],
                files_extra={
                    "selected_poses_sdf": str(selected_sdf),
                    "selection_manifest_csv": str(selection_manifest),
                    "evaluation_input_poses_sdf": str(evaluation_input_sdf),
                    "evaluation_input_manifest_csv": str(evaluation_input_manifest),
                    "batches_dir": str(batches_dir),
                },
                complete=False,
                batch_reports=batch_reports,
                selection_strategy=selection_strategy,
                poses_per_compound=poses_per_compound,
                batch_size=batch_size,
                max_poses=max_poses,
                overall_timeout_seconds=overall_timeout_seconds,
                top_n=top_n,
                operational_errors=profiler_operational_errors,
                screening_mode=screening_mode,
            )
            partial["success"] = False
            partial["error"] = f"{stop_reason}; increase the budget and rerun with resume=true"
            write_json(output_dir / "evaluation_report.json", partial,
                       additional_paths=[output_dir / "evaluation_report.partial.json"])
            return partial

        batch_dir = batches_dir / f"batch_{batch_number:04d}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        batch_report_file = batch_dir / "evaluation_report.json"
        batch_sdf = batch_dir / "poses.sdf"
        _write_selected_sdf(docked_path, batch_poses, batch_sdf, run_cache=run_cache)
        batch_fingerprint = _build_evaluator_batch_fingerprint(
            input_data=input_data,
            receptor_path=receptor_path,
            batch_sdf=batch_sdf,
            batch_poses=batch_poses,
            posecheck_thresholds=posecheck_thresholds,
            posecheck_calibration=posecheck_calibration,
            screening_mode=screening_mode,
            tool_switches=tool_switches,
            run_cache=run_cache,
        )
        resume_validation = None

        if resume and batch_report_file.exists():
            try:
                batch_report = json.loads(batch_report_file.read_text())
                resume_validation = validate_batch_metadata(
                    batch_report,
                    expected_fingerprint=batch_fingerprint,
                )
                if resume_validation.get("reusable"):
                    resume_validation = _validate_resumable_batch_evaluation(
                        batch_report, batch_poses,
                    )
                if resume_validation.get("reusable"):
                    batch_eval = batch_report.get("evaluation", [])
                    for entry in batch_eval:
                        entry.update(_scientific_cohort_for(
                            entry,
                            entry.get("docking_scores") or {},
                            scientific_cohort,
                        ))
                    all_evaluation.extend(batch_eval)
                    batch_reports.append({
                        "batch_number": batch_number,
                        "status": "resumed",
                        "report": str(batch_report_file),
                        "poses": len(batch_eval),
                        "timings_seconds": batch_report.get("timings_seconds", {}),
                        "resume_validation": resume_validation,
                    })
                    global_completed.extend(batch_report.get("steps_completed", []))
                    global_skipped.extend(batch_report.get("steps_skipped", []))
                    profiler_operational_errors.extend(
                        _batch_operational_errors(batch_report, batch_number)
                    )
                    # Cached shards are already durable. Refresh the cumulative
                    # selection at the next checkpoint/finalization, not once
                    # per resumed shard (another quadratic export path).
                    continue
            except (OSError, json.JSONDecodeError) as exc:
                resume_validation = {
                    "reusable": False,
                    "status": "invalid_batch_report",
                    "reason": str(exc),
                }

        batch_scores = _scores_for_selected_batch(batch_poses)
        batch_report = _evaluate_single_sdf(
            input_data=input_data,
            receptor_path=receptor_path,
            docked_path=batch_sdf,
            output_dir=batch_dir,
            posecheck_thresholds=posecheck_thresholds,
            posecheck_calibration=posecheck_calibration,
            reference_ligand_file=reference_ligand_file,
            blueprint_file=blueprint_file,
            pocket_features_file=pocket_features_file,
            docking_scores=batch_scores,
            run_cache=run_cache,
            tool_switches=tool_switches,
            screening_mode=screening_mode,
        )

        batch_eval = _remap_batch_entries(
            batch_report.get("evaluation", []),
            batch_poses,
            posecheck_thresholds=posecheck_thresholds,
            identity_qc_enabled=tool_switches["identity_qc"],
            screening_mode=screening_mode,
        )
        for entry in batch_eval:
            entry.update(_scientific_cohort_for(
                entry,
                entry.get("docking_scores") or {},
                scientific_cohort,
            ))
        _inject_smiles(batch_eval, smiles_map)
        batch_report["evaluation"] = batch_eval
        batch_report["summary"] = _summary_for_evaluation(
            batch_eval,
            batch_report.get("steps_completed", []),
            batch_report.get("steps_skipped", []),
        )
        batch_report["batch"] = {
            "batch_number": batch_number,
            "batch_size": len(batch_poses),
            "source_pose_indices": [p["source_pose_index"] for p in batch_poses],
            "resume_validation": resume_validation,
        }
        batch_report["batch_fingerprint"] = batch_fingerprint
        batch_report["batch_fingerprint_schema"] = batch_fingerprint.get(
            "schema_version"
        )
        batch_report["batch_output_digest"] = build_batch_output_digest(
            batch_eval
        )
        profiler_operational_errors.extend(
            _batch_operational_errors(batch_report, batch_number)
        )
        write_json(batch_report_file, batch_report)
        write_json(batch_dir / "evaluation_summary.json", batch_report["summary"])

        all_evaluation.extend(batch_eval)
        batch_reports.append({
            "batch_number": batch_number,
            "status": "recomputed" if resume_validation else "completed",
            "report": str(batch_report_file),
            "poses": len(batch_eval),
            "timings_seconds": batch_report.get("timings_seconds", {}),
            "resume_validation": resume_validation,
        })
        global_completed.extend(batch_report.get("steps_completed", []))
        global_skipped.extend(batch_report.get("steps_skipped", []))

        # Every shard is durable immediately. Avoid quadratic cumulative JSON /
        # SDF rewrites; publish full checkpoints on the first/last batch or timer.
        now = monotonic()
        if (
            last_checkpoint_at is not None
            and now - last_checkpoint_at < checkpoint_interval
            and len(all_evaluation) < len(evaluation_input_poses)
        ):
            continue
        last_checkpoint_at = now
        post_qc_selected = _refresh_post_qc_selection(
            docked_path,
            evaluation_input_poses,
            all_evaluation,
            selected_sdf,
            selection_manifest,
            run_cache=run_cache,
        )
        _write_final_outputs(
            output_dir=output_dir,
            evaluation=all_evaluation,
            evaluation_input_poses=evaluation_input_poses,
            selected_poses=post_qc_selected,
            posecheck_thresholds=posecheck_thresholds,
            posecheck_calibration=posecheck_calibration,
            posecheck_num_conformers=_posecheck_num_conformers(input_data),
            reference_ligand_file=reference_ligand_file,
            docking_scores_file=docking_scores_file,
            docking_artifact_lineage=docking_artifact_lineage,
            receptor_identity=receptor_identity,
            steps_completed=global_completed,
            steps_skipped=global_skipped,
            files_extra={
                "selected_poses_sdf": str(selected_sdf),
                "selection_manifest_csv": str(selection_manifest),
                "evaluation_input_poses_sdf": str(evaluation_input_sdf),
                "evaluation_input_manifest_csv": str(evaluation_input_manifest),
                "batches_dir": str(batches_dir),
            },
            complete=(
                len(all_evaluation) == len(evaluation_input_poses)
                and not profiler_operational_errors
            ),
            batch_reports=batch_reports,
            selection_strategy=selection_strategy,
            poses_per_compound=poses_per_compound,
            batch_size=batch_size,
            max_poses=max_poses,
            overall_timeout_seconds=overall_timeout_seconds,
            partial=True,
            top_n=top_n,
            operational_errors=profiler_operational_errors,
            screening_mode=screening_mode,
        )

    post_qc_selected = _refresh_post_qc_selection(
        docked_path,
        evaluation_input_poses,
        all_evaluation,
        selected_sdf,
        selection_manifest,
        run_cache=run_cache,
    )
    final_files = {
        "selected_poses_sdf": str(selected_sdf),
        "selection_manifest_csv": str(selection_manifest),
        "evaluation_input_poses_sdf": str(evaluation_input_sdf),
        "evaluation_input_manifest_csv": str(evaluation_input_manifest),
        "batches_dir": str(batches_dir),
    }
    auxiliary_results = {}
    operational_errors: list[dict] = list(profiler_operational_errors)
    topology_result = None
    if tool_switches["rtmscore"]:
        try:
            rtmscore_cohort = _post_qc_survivor_pose_records(
                all_evaluation, evaluation_input_poses
            )
            rtmscore_cohort_sdf = output_dir / "rtmscore_survivor_poses.sdf"
            _write_selected_sdf(
                docked_path, rtmscore_cohort, rtmscore_cohort_sdf
            )
            final_files["rtmscore_survivor_poses_sdf"] = str(
                rtmscore_cohort_sdf
            )
            if not rtmscore_cohort:
                empty_manifest = output_dir / "rtmscore" / "rtmscore_manifest.json"
                empty_manifest.parent.mkdir(parents=True, exist_ok=True)
                rtmscore_result = {
                    "schema_version": "rtmscore-rescoring/v1",
                    "success": True,
                    "status": "not_evaluable",
                    "counts": {"requested": 0, "scored": 0},
                    "result_semantics": "No post-QC survivor was available to rescore.",
                    "pose_scores": [],
                    "manifest_file": str(empty_manifest),
                }
                empty_manifest.write_text(
                    json.dumps(rtmscore_result, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            else:
                rtmscore_result = _run_rtmscore(
                    input_data=input_data,
                    receptor_path=receptor_path,
                    selected_sdf=rtmscore_cohort_sdf,
                    output_dir=output_dir,
                )
            manifest_file = rtmscore_result.get("manifest_file")
            if manifest_file:
                final_files["rtmscore_manifest"] = manifest_file
            auxiliary_results["rtmscore"] = {
                key: rtmscore_result.get(key)
                for key in (
                    "schema_version",
                    "success",
                    "status",
                    "simulation_status",
                    "analysis_status",
                    "counts",
                    "analysis_counts",
                    "protocol",
                    "result_semantics",
                    "manifest_file",
                )
            }
            if (
                rtmscore_result.get("success") is False
                or rtmscore_result.get("status") == "failed"
            ):
                failure_reason = (
                    rtmscore_result.get("error")
                    or "RTMScore returned failed status"
                )
                global_skipped.append({
                    "step": "rtmscore", "reason": failure_reason,
                })
                operational_errors.append({
                    "tool": "rtmscore",
                    "reason": failure_reason,
                })
                rescore_enabled = bool(tool_switches["rescore_rtmscore"])
                global_skipped.append({
                    "step": "rescore_rtmscore",
                    "reason": "RTMScore dependency failed",
                })
                if rescore_enabled:
                    auxiliary_results["rescore_rtmscore"] = {
                        "success": True,
                        "status": "not_evaluable",
                        "reason": "RTMScore dependency failed",
                    }
            elif rtmscore_result.get("status") == "passed":
                global_completed.append("rtmscore")
                rescore_config = (
                    ((input_data.get("tools") or {}).get("rescore_rtmscore") or {})
                )
                rescore_params = dict(rescore_config.get("params") or {})
                use_as_primary = bool(
                    rescore_config.get("enabled", False)
                    and rescore_params.get("use_as_primary_metric", True)
                )
                merged_count = _merge_rtmscore_pose_scores(
                    all_evaluation,
                    rtmscore_result,
                    use_as_primary_metric=use_as_primary,
                )
                auxiliary_results["rtmscore"]["pose_scores_merged"] = merged_count
                if rescore_config.get("enabled", False):
                    global_completed.append("rescore_rtmscore")
                    auxiliary_results["rescore_rtmscore"] = {
                        "success": True,
                        "status": "passed",
                        "use_as_primary_metric": use_as_primary,
                        "pose_scores_merged": merged_count,
                    }
                else:
                    global_skipped.append({
                        "step": "rescore_rtmscore",
                        "reason": "disabled_by_config",
                    })
                # Re-select from every evaluated QC survivor, not merely from
                # the already collapsed one-pose-per-compound set.
                post_qc_selected = _refresh_post_qc_selection(
                    docked_path,
                    evaluation_input_poses,
                    all_evaluation,
                    selected_sdf,
                    selection_manifest,
                    run_cache=run_cache,
                )
            elif rtmscore_result.get("status") == "not_evaluable":
                global_completed.append("rtmscore")
                rescore_enabled = bool(
                    ((input_data.get("tools") or {}).get(
                        "rescore_rtmscore", {}
                    )).get("enabled", False)
                )
                global_skipped.append({
                    "step": "rescore_rtmscore",
                    "reason": "no post-QC survivor was available",
                })
                if rescore_enabled:
                    auxiliary_results["rescore_rtmscore"] = {
                        "success": True,
                        "status": "not_evaluable",
                        "reason": "no post-QC survivor was available",
                    }
        except Exception as exc:
            global_skipped.append({"step": "rtmscore", "reason": str(exc)})
            auxiliary_results["rtmscore"] = {
                "success": False,
                "status": "failed",
                "error": str(exc),
            }
            operational_errors.append({
                "tool": "rtmscore",
                "reason": str(exc),
            })
            global_skipped.append({
                "step": "rescore_rtmscore",
                "reason": "RTMScore dependency failed",
            })
            if tool_switches["rescore_rtmscore"]:
                auxiliary_results["rescore_rtmscore"] = {
                    "success": True,
                    "status": "not_evaluable",
                    "reason": "RTMScore dependency failed",
                }
    else:
        global_skipped.append({"step": "rtmscore", "reason": "disabled_by_config"})
        global_skipped.append({
            "step": "rescore_rtmscore", "reason": "disabled_by_config"
        })

    # Candidate-specific topologies must be built after any RTMScore pose
    # reselection so downstream OpenMM/MM-GBSA top-K cohorts reflect the final
    # post-QC ordering.
    if tool_switches["candidate_topology"]:
        try:
            topology_result = _run_candidate_topology(
                input_data=input_data,
                receptor_path=receptor_path,
                selected_sdf=selected_sdf,
                output_dir=output_dir,
            )
            global_completed.append("candidate_topology")
            manifest_file = topology_result.get("manifest_file")
            if manifest_file:
                final_files["candidate_topologies_manifest"] = manifest_file
            auxiliary_results["candidate_topology"] = {
                key: topology_result.get(key)
                for key in (
                    "schema_version",
                    "success",
                    "status",
                    "counts",
                    "manifest_file",
                )
            }
        except Exception as exc:
            global_skipped.append({
                "step": "candidate_topology",
                "reason": str(exc),
            })
            auxiliary_results["candidate_topology"] = {
                "success": False,
                "status": "not_evaluable",
                "error": str(exc),
            }
    else:
        global_skipped.append({
            "step": "candidate_topology",
            "reason": "disabled_by_config",
        })

    if tool_switches["openmm_md"]:
        try:
            topology_manifest = (
                topology_result.get("manifest_file")
                if isinstance(topology_result, dict)
                else None
            )
            if not topology_manifest:
                raise RuntimeError(
                    "OpenMM short MD requires candidate_topologies.json from "
                    "the enabled candidate_topology tool"
                )
            short_md_result = _run_openmm_md(
                input_data=input_data,
                candidate_topologies_manifest=topology_manifest,
                output_dir=output_dir,
            )
            if tool_switches["mdanalysis"]:
                merge_md_metrics(all_evaluation, short_md_result)
            global_completed.append("openmm_md")
            manifest_file = short_md_result.get("manifest_file")
            if manifest_file:
                final_files["short_md_manifest"] = manifest_file
            stability_file = short_md_result.get("trajectory_stability_report")
            if stability_file:
                final_files["trajectory_stability_report"] = stability_file
            auxiliary_results["openmm_md"] = {
                key: short_md_result.get(key)
                for key in (
                    "schema_version",
                    "success",
                    "status",
                    "counts",
                    "protocol",
                    "result_semantics",
                    "manifest_file",
                    "trajectory_stability_report",
                )
            }
            if tool_switches["mdanalysis"]:
                analysis_status = short_md_result.get("analysis_status")
                auxiliary_results["mdanalysis"] = {
                    "success": analysis_status != "failed",
                    "status": analysis_status,
                    "counts": short_md_result.get("analysis_counts"),
                    "trajectory_stability_report": stability_file,
                }
                if analysis_status == "passed":
                    global_completed.append("mdanalysis")
                else:
                    global_skipped.append({
                        "step": "mdanalysis",
                        "reason": f"trajectory analysis status={analysis_status}",
                    })
                if analysis_status == "failed":
                    reasons = [
                        str(entry.get("analysis_reason"))
                        for entry in short_md_result.get("candidates") or []
                        if entry.get("analysis_status") == "failed"
                    ]
                    operational_errors.append({
                        "tool": "mdanalysis",
                        "reason": "; ".join(reasons)
                        or "one or more trajectory analyses failed",
                    })
            else:
                global_skipped.append({
                    "step": "mdanalysis", "reason": "disabled_by_config"
                })
            if (
                short_md_result.get("success") is False
                or short_md_result.get("status") == "failed"
            ):
                operational_errors.append({
                    "tool": "openmm_md",
                    "reason": short_md_result.get("error")
                    or "one or more short-MD candidates failed",
                })
        except Exception as exc:
            global_skipped.append({"step": "openmm_md", "reason": str(exc)})
            auxiliary_results["openmm_md"] = {
                "success": False,
                "status": "failed",
                "error": str(exc),
            }
            operational_errors.append({
                "tool": "openmm_md",
                "reason": str(exc),
            })
            if tool_switches["mdanalysis"]:
                global_skipped.append({
                    "step": "mdanalysis",
                    "reason": "OpenMM dependency failed: " + str(exc),
                })
                auxiliary_results["mdanalysis"] = {
                    "success": True,
                    "status": "not_evaluable",
                    "reason": "OpenMM dependency failed: " + str(exc),
                }
    else:
        global_skipped.append({"step": "openmm_md", "reason": "disabled_by_config"})
        global_skipped.append({"step": "mdanalysis", "reason": "disabled_by_config"})

    if tool_switches["mmgbsa"]:
        try:
            topology_manifest = (
                topology_result.get("manifest_file")
                if isinstance(topology_result, dict)
                else None
            )
            if not topology_manifest:
                raise RuntimeError(
                    "MM-GBSA requires a candidate_topologies.json artifact "
                    "from the enabled candidate_topology tool"
                )
            mmgbsa_result = _run_mmgbsa(
                input_data=input_data,
                candidate_topologies_manifest=topology_manifest,
                output_dir=output_dir,
            )
            global_completed.append("mmgbsa")
            manifest_file = mmgbsa_result.get("manifest_file")
            if manifest_file:
                final_files["mmgbsa_manifest"] = manifest_file
            auxiliary_results["mmgbsa"] = {
                key: mmgbsa_result.get(key)
                for key in (
                    "schema_version",
                    "success",
                    "status",
                    "reason_code",
                    "execution_mode",
                    "counts",
                    "result_semantics",
                    "ranking",
                    "manifest_file",
                )
            }
            if (
                mmgbsa_result.get("success") is False
                or mmgbsa_result.get("status") == "failed"
            ):
                operational_errors.append({
                    "tool": "mmgbsa",
                    "reason": (
                        mmgbsa_result.get("error")
                        or mmgbsa_result.get("reason_code")
                        or "Amber MM-GBSA returned failed status"
                    ),
                })
            elif mmgbsa_result.get("status") == "passed":
                mmgbsa_params = dict(
                    (
                        ((input_data.get("tools") or {}).get("mmgbsa") or {})
                    ).get("params")
                    or {}
                )
                merged_count = _merge_mmgbsa_pose_scores(
                    all_evaluation,
                    mmgbsa_result,
                    use_as_primary_metric=bool(
                        mmgbsa_params.get("use_as_primary_metric", True)
                    ),
                )
                auxiliary_results["mmgbsa"][
                    "pose_scores_merged"
                ] = merged_count
                # MM-GBSA is a registered DYN metric with its own weak sanity
                # gate.  It is attached after the initial QC pass, so the old
                # verdict/assessment cannot be reused: a newly failing pose
                # must be removed and the next eligible pose for that compound
                # considered.
                for entry in all_evaluation:
                    entry["verdict"] = _compute_verdict(
                        entry,
                        posecheck_thresholds,
                        screening_mode=screening_mode,
                    )
                    _attach_diagnosis(
                        entry,
                        subpockets=None,
                        posecheck_thresholds=posecheck_thresholds,
                    )
                post_qc_selected = _refresh_post_qc_selection(
                    docked_path,
                    evaluation_input_poses,
                    all_evaluation,
                    selected_sdf,
                    selection_manifest,
                    run_cache=run_cache,
                )
        except Exception as exc:
            global_skipped.append({
                "step": "mmgbsa",
                "reason": str(exc),
            })
            auxiliary_results["mmgbsa"] = {
                "success": False,
                "status": "failed",
                "error": str(exc),
            }
            operational_errors.append({
                "tool": "mmgbsa",
                "reason": str(exc),
            })
    else:
        global_skipped.append({
            "step": "mmgbsa",
            "reason": "disabled_by_config",
        })

    if tool_switches["openfe_rbfe"]:
        try:
            openfe_result = _run_openfe_rbfe(
                input_data=input_data,
                receptor_path=receptor_path,
                selected_sdf=selected_sdf,
                output_dir=output_dir,
            )
            global_completed.append("openfe_rbfe")
            manifest_file = openfe_result.get("manifest_file")
            if manifest_file:
                final_files["openfe_rbfe_manifest"] = manifest_file
            auxiliary_results["openfe_rbfe"] = {
                key: openfe_result.get(key)
                for key in (
                    "schema_version",
                    "success",
                    "status",
                    "execution_mode",
                    "counts",
                    "applicability",
                    "result_semantics",
                    "network",
                    "manifest_file",
                )
            }
            if (
                openfe_result.get("success") is False
                or openfe_result.get("status") == "failed"
            ):
                operational_errors.append({
                    "tool": "openfe_rbfe",
                    "reason": (
                        openfe_result.get("error")
                        or openfe_result.get("reason_code")
                        or "OpenFE planning returned failed status"
                    ),
                })
        except Exception as exc:
            global_skipped.append({
                "step": "openfe_rbfe",
                "reason": str(exc),
            })
            auxiliary_results["openfe_rbfe"] = {
                "success": False,
                "status": "failed",
                "error": str(exc),
            }
            operational_errors.append({
                "tool": "openfe_rbfe",
                "reason": str(exc),
            })
    else:
        global_skipped.append({
            "step": "openfe_rbfe",
            "reason": "disabled_by_config",
        })

    if tool_switches["ties_rbfe"]:
        try:
            ties_result = _run_ties_rbfe(
                input_data=input_data,
                receptor_path=receptor_path,
                selected_sdf=selected_sdf,
                output_dir=output_dir,
            )
            global_completed.append("ties_rbfe")
            manifest_file = ties_result.get("manifest_file")
            if manifest_file:
                final_files["ties_rbfe_manifest"] = manifest_file
            auxiliary_results["ties_rbfe"] = {
                key: ties_result.get(key)
                for key in (
                    "schema_version",
                    "success",
                    "status",
                    "reason_code",
                    "candidate_counts",
                    "transformation_counts",
                    "receptor_validation",
                    "free_energy_result",
                    "manifest_file",
                )
            }
            if (
                ties_result.get("success") is False
                or ties_result.get("status") == "failed"
            ):
                operational_errors.append({
                    "tool": "ties_rbfe",
                    "reason": (
                        ties_result.get("error")
                        or ties_result.get("reason_code")
                        or "TIES planning returned failed status"
                    ),
                })
        except Exception as exc:
            global_skipped.append({
                "step": "ties_rbfe",
                "reason": str(exc),
            })
            auxiliary_results["ties_rbfe"] = {
                "success": False,
                "status": "failed",
                "error": str(exc),
            }
            operational_errors.append({
                "tool": "ties_rbfe",
                "reason": str(exc),
            })
    else:
        global_skipped.append({
            "step": "ties_rbfe",
            "reason": "disabled_by_config",
        })

    for entry in all_evaluation:
        entry["verdict"] = _compute_verdict(entry, posecheck_thresholds, screening_mode=screening_mode)
    post_qc_selected = _refresh_post_qc_selection(
        docked_path, evaluation_input_poses, all_evaluation, selected_sdf, selection_manifest,
        run_cache=run_cache,
    )

    final_report = _write_final_outputs(
        output_dir=output_dir,
        evaluation=all_evaluation,
        evaluation_input_poses=evaluation_input_poses,
        selected_poses=post_qc_selected,
        posecheck_thresholds=posecheck_thresholds,
        posecheck_calibration=posecheck_calibration,
        posecheck_num_conformers=_posecheck_num_conformers(input_data),
        reference_ligand_file=reference_ligand_file,
        docking_scores_file=docking_scores_file,
        docking_artifact_lineage=docking_artifact_lineage,
        receptor_identity=receptor_identity,
        steps_completed=global_completed,
        steps_skipped=global_skipped,
        files_extra=final_files,
        complete=not operational_errors,
        batch_reports=batch_reports,
        selection_strategy=selection_strategy,
        poses_per_compound=poses_per_compound,
        batch_size=batch_size,
        max_poses=max_poses,
        overall_timeout_seconds=overall_timeout_seconds,
        top_n=top_n,
        auxiliary_results=auxiliary_results,
        partial=bool(operational_errors),
        operational_errors=operational_errors,
        screening_mode=screening_mode,
    )
    if operational_errors:
        final_report["error"] = (
            "one or more explicitly enabled Step 5 tools failed "
            "operationally"
        )
        final_report["operational_errors"] = _dedupe_dicts(operational_errors)
        write_json(output_dir / "evaluation_report.json", final_report,
                   additional_paths=[output_dir / "evaluation_report.partial.json"])
    return final_report


def _batch_operational_errors(batch_report: dict, batch_number: int) -> list[dict]:
    """Attach batch provenance to operational failures from one shard."""
    errors = []
    for raw_error in batch_report.get("operational_errors") or []:
        if not isinstance(raw_error, dict):
            continue
        error = dict(raw_error)
        error.setdefault("batch_number", batch_number)
        errors.append(error)
    return errors


def _failed_profiler_results(
    profiler: str,
    pose_names: list[str],
    reason: str,
) -> tuple[list[dict], dict]:
    """Materialize a run-wide profiler fault on every applicable pose."""
    results = [
        {
            "pose_index": pose_index,
            "mol_name": mol_name,
            "success": False,
            "error": reason,
            "num_interactions": None,
            "interactions": [],
        }
        for pose_index, mol_name in enumerate(pose_names)
    ]
    return results, _profiler_operational_error(profiler, results, len(pose_names))


def _validated_profiler_results(
    profiler: str,
    results: list[dict] | None,
    pose_names: list[str],
) -> tuple[list[dict], dict | None]:
    """Require complete per-pose coverage from an explicitly enabled profiler.

    An explicit scientific ``not_evaluable``/``not_applicable`` row remains a
    valid result. Missing rows, exceptions, and rows that report failure are
    operational faults; they are kept as per-pose diagnostics rather than
    converted into zero interactions or a scientific rejection.
    """
    indexed = _results_by_pose_index(results)
    normalized: list[dict] = []
    failed: list[dict] = []
    for pose_index, mol_name in enumerate(pose_names):
        raw = indexed.get(pose_index)
        if raw is None:
            row = {
                "pose_index": pose_index,
                "mol_name": mol_name,
                "success": False,
                "error": f"{profiler} returned no result for pose {pose_index}",
                "num_interactions": None,
                "interactions": [],
            }
        else:
            row = {
                **raw,
                "pose_index": pose_index,
                "mol_name": raw.get("mol_name") or mol_name,
            }
        normalized.append(row)
        if _profiler_result_failed_operationally(row):
            failed.append(row)

    error = (
        _profiler_operational_error(profiler, failed, len(pose_names))
        if failed else None
    )
    return normalized, error


def _profiler_result_failed_operationally(result: dict) -> bool:
    """Distinguish a tool fault from an explicit scientific non-evaluation."""
    status = str(result.get("status") or "").strip().lower()
    if (
        status in {"not_evaluable", "not_applicable"}
        and result.get("success") is not False
        and not result.get("error")
    ):
        return False
    return result.get("success") is False or bool(result.get("error"))


def _profiler_has_usable_evidence(results: list[dict] | None) -> bool:
    """Whether at least one pose contains successful interaction evidence."""
    return any(
        result.get("success") is True
        and not result.get("error")
        and str(result.get("status") or "").strip().lower()
        not in {"not_evaluable", "not_applicable"}
        for result in results or []
        if isinstance(result, dict)
    )


def _profiler_operational_error(
    profiler: str,
    failed_results: list[dict],
    applicable_pose_count: int,
) -> dict:
    """Summarize profiler failure without erasing per-pose error records."""
    pose_indices = [
        result.get("pose_index") for result in failed_results
        if result.get("pose_index") is not None
    ]
    reasons = _dedupe([
        str(result.get("error") or f"{profiler} pose evaluation failed")
        for result in failed_results
    ])
    reason = (
        f"{profiler} failed for {len(failed_results)}/"
        f"{applicable_pose_count} applicable poses: "
        + "; ".join(reasons)
    )
    return {
        "tool": profiler,
        "step": f"{profiler}_interactions",
        "code": "profiler_operational_failure",
        "scope": "batch",
        "failed_pose_count": len(failed_results),
        "applicable_pose_count": applicable_pose_count,
        "pose_indices": pose_indices,
        "reason": reason,
    }


def _evaluate_single_sdf(
    input_data: dict,
    receptor_path: Path,
    docked_path: Path,
    output_dir: Path,
    posecheck_thresholds: dict,
    posecheck_calibration: dict,
    reference_ligand_file: str | None,
    blueprint_file: str | None,
    pocket_features_file: str | None,
    docking_scores: list[dict] | None,
    tool_switches: dict | None = None,
    screening_mode: object = DEFAULT_SCREENING_MODE,
    run_cache: dict | None = None,
) -> dict:
    """Evaluate one SDF shard. Pose indices are local to the shard."""

    batch_started = monotonic()
    timings = {}
    switches = tool_switches or dict(EVALUATOR_TOOL_DEFAULTS)
    steps_completed = []
    steps_skipped = []
    operational_errors: list[dict] = []
    pose_names = _read_pose_names(docked_path)

    # These profilers have independent inputs and output objects. Submit only
    # enabled/available work, then merge in fixed order so scheduling cannot
    # change reporting, identity joins or scientific gates. PLIP retains its
    # separate interpreter (OpenBabel must never be imported alongside RDKit).
    tasks = {}
    for name, checker, runner in (
        ("posebusters", check_posebusters_available, run_posebusters),
        ("posecheck", check_posecheck_available, run_posecheck),
    ):
        if not switches[name]:
            steps_skipped.append({"step": name, "reason": "disabled_by_config"})
        elif checker():
            if name == "posecheck":
                tasks[name] = lambda: run_posecheck(
                    docked_path, receptor_path,
                    num_conformers=_posecheck_num_conformers(input_data),
                )
            else:
                tasks[name] = lambda runner=runner: runner(docked_path, receptor_path)
        else:
            steps_skipped.append({"step": name, "reason": f"not installed (pip install {name})"})
    if switches["plip"]:
        tasks["plip"] = lambda: run_plip_per_pose(receptor_path, docked_path)
    else:
        steps_skipped.append({"step": "plip_interactions", "reason": "disabled_by_config"})
    if switches["prolif"]:
        prolif_params = ((input_data.get("tools") or {}).get("prolif") or {}).get("params") or {}
        tasks["prolif"] = lambda: run_prolif_per_pose(
            receptor_path, docked_path,
            chain=input_data.get("chain", "A"),
            vicinity_cutoff=float(prolif_params.get("vicinity_cutoff", 6.5)),
            interaction_types=prolif_params.get("interaction_types"),
            reference_ligand_file=Path(reference_ligand_file) if reference_ligand_file else None,
        )
    else:
        steps_skipped.append({"step": "prolif_interactions", "reason": "disabled_by_config"})

    profiler_results = {}
    for name, (results, error) in run_profiler_tasks(
        tasks, profiler_worker_count(input_data.get("profiler_workers", 2)),
        timings=timings,
    ).items():
        if name in {"plip", "prolif"}:
            try:
                if error is not None:
                    raise error
                results, profiler_error = _validated_profiler_results(name, results, pose_names)
            except Exception as exc:
                results, profiler_error = _failed_profiler_results(name, pose_names, str(exc))
            if profiler_error:
                if name != "prolif":
                    operational_errors.append(profiler_error)
                steps_skipped.append({"step": f"{name}_interactions", "reason": profiler_error["reason"]})
            else:
                steps_completed.append(f"{name}_interactions")
        elif error is not None:
            steps_skipped.append({"step": name, "reason": str(error)})
            operational_errors.append({"tool": name, "reason": str(error)})
            results = [
                {"pose_index": index, "mol_name": pose_name, "success": False, "error": str(error)}
                for index, pose_name in enumerate(pose_names)
            ]
        else:
            steps_completed.append(name)
            for row in results or []:
                if _profiler_result_failed_operationally(row):
                    operational_errors.append({
                        "tool": name, "pose_index": row.get("pose_index"),
                        "reason": row.get("error") or f"{name} pose evaluation failed",
                    })
        profiler_results[name] = results
    pb_results = profiler_results.get("posebusters")
    pc_results = profiler_results.get("posecheck")
    plip_results = profiler_results.get("plip")
    prolif_results = profiler_results.get("prolif")


    # === Blueprint matching (PLIP + ProLIF evidence) ===
    blueprint_results = None
    blueprint = None
    blueprint_applicable = True
    blueprint_has_constraints = False
    if not switches["blueprint_match"]:
        steps_skipped.append({"step": "blueprint_match", "reason": "disabled_by_config"})
    elif not blueprint_file:
        steps_skipped.append({
            "step": "blueprint_match",
            "status": "not_applicable",
            "reason": "optional interaction_blueprint_file was not provided",
        })
    else:
        bp_path = Path(blueprint_file)
        if bp_path.exists():
            try:
                blueprint = cached_file(run_cache, bp_path, "blueprint", load_blueprint)
                bp_validation = blueprint.get("validation") or {}
                blueprint_applicable = not bool(
                    bp_validation.get("not_applicable") is True
                    or bp_validation.get("applicable") is False
                    or blueprint.get("not_applicable") is True
                    or blueprint.get("applicable") is False
                )
                blueprint_has_constraints = bool(
                    any(blueprint.get(level) for level in ("required", "desired", "optional"))
                    or blueprint.get("required_groups")
                    or (blueprint.get("constraints") or {}).get("required_groups")
                )
                if not blueprint_applicable:
                    steps_skipped.append({
                        "step": "blueprint_match",
                        "reason": bp_validation.get("reason", "blueprint not applicable"),
                    })
                elif not blueprint_has_constraints:
                    steps_skipped.append({
                        "step": "blueprint_match",
                        "reason": "blueprint contains no interaction constraints",
                    })
                else:
                    if not switches["plip"]:
                        operational_errors.append({
                            "tool": "blueprint_match",
                            "step": "blueprint_match",
                            "code": "blueprint_profiler_required",
                            "reason": (
                                "an applicable interaction blueprint requires "
                                "PLIP to be enabled"
                            ),
                        })
                    blueprint_results = _match_blueprint_with_profilers(
                        blueprint,
                        pose_names,
                        {
                            "plip": plip_results,
                            "prolif": prolif_results,
                        },
                        {
                            "plip": switches["plip"],
                            "prolif": switches["prolif"],
                        },
                    )
                    if not _profiler_has_usable_evidence(plip_results):
                        steps_skipped.append({
                            "step": "blueprint_match",
                            "reason": (
                                "PLIP results unavailable; "
                                "hard constraints unresolved"
                            ),
                        })
                    if blueprint_results:
                        steps_completed.append("blueprint_match")
            except Exception as e:
                reason = str(e)
                steps_skipped.append({"step": "blueprint_match", "reason": reason})
                operational_errors.append({
                    "tool": "blueprint_match",
                    "step": "blueprint_match",
                    "code": "blueprint_match_failed",
                    "reason": reason,
                })
        else:
            reason = f"blueprint file not found: {bp_path}"
            steps_skipped.append({"step": "blueprint_match", "reason": reason})
            operational_errors.append({
                "tool": "blueprint_match",
                "step": "blueprint_match",
                "code": "blueprint_file_missing",
                "reason": reason,
            })

    # === Combine per-molecule verdicts ===
    # blueprint_expected drives the "blueprint_missing" soft-warn. A user who
    # turned blueprint_match off did not "expect" it, so opting out must not add
    # a warn penalty. It also cannot lift a pose to pass: the
    # evidence-count invariant in _compute_verdict is unchanged.
    blueprint_expected = (
        bool(blueprint_file)
        and switches["blueprint_match"]
        and blueprint_applicable
        and blueprint_has_constraints
    )
    pocket_expected = bool(pocket_features_file) and switches["pocket_occupancy"]
    evaluation = _combine_results(
        pb_results,
        pc_results,
        plip_results,
        blueprint_results,
        pose_names,
        posecheck_thresholds,
        docking_scores=docking_scores,
        blueprint_expected=blueprint_expected,
        pocket_expected=pocket_expected,
        posecheck_calibrated=posecheck_calibration.get("mode") == "delta_from_reference",
        tool_switches=switches,
        prolif_results=prolif_results,
        screening_mode=screening_mode,
    )

    # === Pocket occupancy + Diagnosis ===
    subpockets = None
    if not switches["pocket_occupancy"]:
        steps_skipped.append({"step": "pocket_occupancy", "reason": "disabled_by_config"})
    elif not pocket_features_file:
        steps_skipped.append({
            "step": "pocket_occupancy",
            "status": "not_applicable",
            "reason": "optional pocket_features_file was not provided",
        })
    else:
        pf_path = Path(pocket_features_file)
        if pf_path.exists():
            try:
                subpockets = cached_file(run_cache, pf_path, "subpockets", load_subpockets)
                if subpockets:
                    chain = input_data.get("chain", "A")
                    coords_key = ("subpocket_coords", chain, json.dumps(subpockets, sort_keys=True))
                    subpocket_coords = cached_file(
                        run_cache, receptor_path, coords_key,
                        lambda path: prepare_subpocket_coords(path, subpockets, chain),
                    )
                    from rdkit import Chem
                    supplier = Chem.SDMolSupplier(str(docked_path), removeHs=False)
                    for i, entry in enumerate(evaluation):
                        mol = supplier[i] if i < len(supplier) else None
                        if mol is None:
                            raise ValueError(
                                f"RDKit failed to parse pose {i} for pocket occupancy"
                            )
                        occ = compute_pocket_occupancy(
                            receptor_path, mol, subpockets,
                            chain=chain, subpocket_coords=subpocket_coords,
                        )
                        entry["pocket_contact_fraction"] = occ
                    steps_completed.append("pocket_occupancy")
                else:
                    steps_skipped.append({
                        "step": "pocket_occupancy",
                        "reason": "pocket features contain no key subpockets",
                    })
            except Exception as e:
                reason = str(e)
                steps_skipped.append({"step": "pocket_occupancy", "reason": reason})
                operational_errors.append({
                    "tool": "pocket_occupancy",
                    "step": "pocket_occupancy",
                    "code": "pocket_occupancy_failed",
                    "reason": reason,
                })
        else:
            reason = f"pocket features file not found: {pf_path}"
            steps_skipped.append({"step": "pocket_occupancy", "reason": reason})
            operational_errors.append({
                "tool": "pocket_occupancy",
                "step": "pocket_occupancy",
                "code": "pocket_features_file_missing",
                "reason": reason,
            })

    operational_errors.extend(attach_supplemental_interactions(
        evaluation, receptor_path, docked_path, output_dir,
        {key: {**(((input_data.get("tools") or {}).get(key)) or {}),
               "enabled": switches.get(key, False)} for key in ("arpeggio", "luna", "oddt")},
    ))
    for key in ("arpeggio", "luna", "oddt"):
        if switches.get(key, False) and not any(error.get("tool") == key for error in operational_errors):
            steps_completed.append(key)

    geometry_started = monotonic()
    operational_errors.extend(attach_pose_geometry(
        evaluation, receptor_path, docked_path,
        {key: {**(((input_data.get("tools") or {}).get(key)) or {}),
               "enabled": switches.get(key, False)} for key in ("rdkit_pose", "freesasa")},
    ))
    timings["pose_geometry"] = round(monotonic() - geometry_started, 3)
    for key in ("rdkit_pose", "freesasa"):
        if switches.get(key, False) and not any(error.get("tool") == key for error in operational_errors):
            steps_completed.append(key)

    # Pocket occupancy is attached after the initial result merge, so recompute
    # the verdict once every currently available registered scalar is present.
    for entry in evaluation:
        entry["verdict"] = _compute_verdict(
            entry, posecheck_thresholds, screening_mode=screening_mode
        )
        _attach_diagnosis(
            entry,
            subpockets=subpockets,
            posecheck_thresholds=posecheck_thresholds,
        )

    summary = _summary_for_evaluation(evaluation, steps_completed, steps_skipped)

    # === Write outputs ===
    report = {
        "success": not operational_errors,
        "timings_seconds": {**timings, "batch_wall": round(monotonic() - batch_started, 3)},
        "summary": summary,
        "evaluation": evaluation,
        "parameters": {
            "posecheck_thresholds": posecheck_thresholds,
            "posecheck_num_conformers": _posecheck_num_conformers(input_data),
            "posecheck_calibration": posecheck_calibration,
            "reference_ligand_file": reference_ligand_file,
            **_scientific_cohort_for(input_data),
        },
        "steps_completed": steps_completed,
        "steps_skipped": steps_skipped,
    }
    if operational_errors:
        report["operational_errors"] = operational_errors
    return report


def _match_blueprint_with_profilers(
    blueprint: dict,
    pose_names: list[str],
    profiler_results: dict[str, list[dict] | None],
    profiler_enabled: dict[str, bool],
) -> list[dict]:
    """Evaluate a blueprint per profiler, then combine without unioning hits."""
    indexed = {
        name: _results_by_pose_index(results)
        for name, results in profiler_results.items()
    }
    combined_results = []
    for pose_index, mol_name in enumerate(pose_names):
        matches: dict[str, dict] = {}
        for profiler in ("plip", "prolif"):
            if not profiler_enabled.get(profiler, False):
                continue
            pose_data = indexed.get(profiler, {}).get(pose_index)
            profiler_status = str(
                (pose_data or {}).get("status") or ""
            ).strip().lower()
            if (
                pose_data
                and pose_data.get("success", False)
                and profiler_status not in {"not_evaluable", "not_applicable"}
            ):
                matches[profiler] = match_pose_to_blueprint(
                    pose_data.get("interactions", []),
                    blueprint,
                    evaluator_capabilities=pose_data.get("evaluator_capabilities"),
                )
            else:
                reason = (
                    (pose_data or {}).get("error")
                    or (pose_data or {}).get("reason")
                    or f"{profiler} pose analysis unavailable"
                )
                matches[profiler] = unresolved_match_result(blueprint, reason)
                if profiler_status:
                    matches[profiler]["profiler_status"] = profiler_status
        combined_results.append({
            "pose_index": pose_index,
            "mol_name": mol_name,
            **_combine_profiler_blueprint_matches(matches, blueprint),
        })
    return combined_results


def _combine_profiler_blueprint_matches(
    matches: dict[str, dict],
    blueprint: dict,
) -> dict:
    """Use PLIP for blueprint gates and keep ProLIF as independent audit evidence."""
    plip = matches.get("plip")
    prolif = matches.get("prolif")
    if plip is None:
        base = unresolved_match_result(
            blueprint, "PLIP unavailable; ProLIF is record-only",
        )
    else:
        base = dict(plip)

    if plip is not None and prolif is not None:
        plip_verdict = plip.get("verdict")
        prolif_verdict = prolif.get("verdict")
        if plip_verdict in {"pass", "fail"} and prolif_verdict in {"pass", "fail"}:
            agreement = "consensus" if plip_verdict == prolif_verdict else "disagreed"
        else:
            agreement = "unresolved"
    elif plip is not None:
        agreement = "single_profiler"
    elif prolif is not None:
        agreement = "audit_only"
    else:
        agreement = "unresolved"

    base.update({
        "profiler_results": matches,
        "profiler_agreement": agreement,
        "reliability": {
            "consensus": "supported",
            "disagreed": "profiler_disagreement",
            "single_profiler": "single_profiler",
            "audit_only": "not_evaluable",
            "unresolved": "not_evaluable",
        }[agreement],
        "profiler_verdicts": {
            name: result.get("verdict")
            for name, result in matches.items()
        },
        "desired_profiler_scores": {
            name: result.get("desired_ranking_score")
            for name, result in matches.items()
        },
    })
    if agreement == "disagreed":
        base["reason"] = "PLIP blueprint verdict retained; ProLIF disagreement recorded for audit"
    return base

def _rdkit_available() -> bool:
    """Whether RDKit can be imported in this environment."""
    try:
        from rdkit import Chem  # noqa: F401
    except ImportError:
        return False
    return True


def _read_pose_records(docked_sdf: Path, *, run_cache: dict | None = None) -> list[dict]:
    """Read each physical SDF block without allowing parser index compression."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        # ``run`` refuses to start without RDKit and says so once, against the
        # run.  Raising here keeps a direct caller on the same footing instead of
        # handing back a record set whose every pose carries an environment fault
        # as its identity verdict.
        raise PoseIdentityError(
            "rdkit_unavailable",
            "RDKit is required for stable pose identity",
        ) from exc

    blocks = cached_file(run_cache, docked_sdf, "sdf_blocks", _read_sdf_blocks)
    if not blocks:
        raise PoseIdentityError(
            "empty_or_unframed_docked_sdf",
            "docked SDF contains no physical molecule records",
            details={"path": str(docked_sdf)},
        )

    def _pose_record_from_block(i: int, block: str) -> dict:
        """Bind one physical SDF block to a verified pose identity.

        Raises ``PoseIdentityError`` for that block alone.  The caller decides
        whether the defect is a per-pose IDQC reject or a systemic artifact
        failure.
        """
        framed = block if block.endswith(("\n", "\r")) else block + "\n"
        payload = (framed + "$$$$\n").encode("utf-8")
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(payload), removeHs=False, strictParsing=False,
        )
        molecules = list(supplier)
        if len(molecules) != 1 or molecules[0] is None:
            raise PoseIdentityError(
                "docked_sdf_record_parse_failed",
                f"RDKit cannot parse physical docked SDF record {i}",
                details={"path": str(docked_sdf), "source_pose_index": i},
            )
        mol = molecules[0]
        name = (
            mol.GetProp("_Name").strip()
            if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
            else f"pose_{i}"
        )
        record = {
            "source_pose_index": i,
            "mol_name": name,
            "_sdf_molblock_sha256": sdf_molblock_sha256(block),
            # Where this pose physically lives.  The published final report
            # has a ``pose_file`` column whose whole purpose is letting a
            # reviewer open the coordinates behind a row; paired with
            # ``source_pose_index`` it addresses exactly one pose.
            "pose_file": str(docked_sdf),
        }
        sdf_properties = {
            str(prop): mol.GetProp(str(prop)).strip()
            for prop in mol.GetPropNames(
                includePrivate=False, includeComputed=False,
            )
        }
        record["_sdf_join_properties"] = sdf_properties
        sdf_cohort = {
            field: sdf_properties[field]
            for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS
            if sdf_properties.get(field)
        }
        for field in SCIENTIFIC_COHORT_PROVENANCE_FIELDS:
            value = sdf_cohort.get(field)
            if not isinstance(value, str):
                continue
            try:
                sdf_cohort[field] = json.loads(value)
            except json.JSONDecodeError:
                pass
        if sdf_properties.get("protocol_id"):
            sdf_cohort.setdefault(
                "protocol_profile_id", sdf_properties["protocol_id"]
            )
        record.update(_scientific_cohort_for(sdf_cohort))

        compound_id = sdf_properties.get("compound_id")
        lineage_source_pose_index = i
        stored_source_pose_index = sdf_properties.get("source_pose_index")
        if stored_source_pose_index not in (None, ""):
            try:
                lineage_source_pose_index = _strict_integer(
                    stored_source_pose_index
                )
            except ValueError as exc:
                raise PoseIdentityError(
                    "invalid_docked_pose_metadata",
                    "source_pose_index is not an integer at physical record "
                    f"{i}",
                ) from exc
            if lineage_source_pose_index < 0:
                raise PoseIdentityError(
                    "invalid_docked_pose_metadata",
                    "source_pose_index is negative at physical record "
                    f"{i}",
                )
        stored_pose_index = sdf_properties.get("pose_index")
        if stored_pose_index not in (None, ""):
            try:
                stored_pose_index = _strict_integer(stored_pose_index)
            except ValueError as exc:
                raise PoseIdentityError(
                    "invalid_docked_pose_metadata",
                    f"pose_index is not an integer at physical record {i}",
                ) from exc
            if stored_pose_index != i:
                raise PoseIdentityError(
                    "docked_sdf_identity_mismatch",
                    "SDF pose_index property disagrees with physical record "
                    f"{i}",
                    details={
                        "source_pose_index": i,
                        "sdf_pose_index": stored_pose_index,
                    },
                )
            record["artifact_pose_index"] = stored_pose_index
        for field in (
            "compound_id", "source_molecule_index",
            "preparation_variant_index", "engine_pose_index",
            "engine_pose_rank",
        ):
            value = sdf_properties.get(field)
            if value in (None, ""):
                continue
            if field != "compound_id":
                try:
                    value = _strict_integer(value)
                except ValueError as exc:
                    raise PoseIdentityError(
                        "invalid_docked_pose_metadata",
                        f"{field} is not an integer at physical record {i}",
                        details={"source_pose_index": i, "field": field},
                    ) from exc
                minimum = 1 if field == "engine_pose_rank" else 0
                if value < minimum:
                    raise PoseIdentityError(
                        "invalid_docked_pose_metadata",
                        f"{field} is outside its valid range at physical "
                        f"record {i}",
                        details={"source_pose_index": i, "field": field},
                    )
            record[field] = value
        for field in LIGAND_LINEAGE_PROPERTIES:
            value = sdf_properties.get(field)
            if value not in (None, ""):
                record[field] = value

        parent_pose_id = sdf_properties.get("parent_pose_id") or None
        coordinate_revision = 0
        if parent_pose_id:
            coordinate_revision = 1
            if sdf_properties.get("coordinate_revision"):
                try:
                    coordinate_revision = max(
                        1, int(sdf_properties["coordinate_revision"]),
                    )
                except ValueError as exc:
                    raise PoseIdentityError(
                        "invalid_coordinate_revision",
                        "coordinate_revision SDF property is not an integer",
                    ) from exc
        source_identity = _source_identity(mol, name)
        if (
            compound_id
            and source_identity.get("compound_id") not in (None, "")
            and str(source_identity["compound_id"]) != str(compound_id)
        ):
            raise PoseIdentityError(
                "docked_sdf_identity_mismatch",
                "source_identity.compound_id disagrees with the explicit SDF "
                f"compound_id at record {i}",
                details={"source_pose_index": i},
            )
        identity = build_pose_identity(
            mol,
            source_identity=source_identity,
            parent_pose_id=parent_pose_id,
            coordinate_revision=coordinate_revision,
            lineage_metadata={
                "operation": "docking_output",
                "source_pose_index": lineage_source_pose_index,
            },
        )
        computed_fields = {
            "chemical_identity_id": identity["chemical_identity"]["identity_id"],
            "chemical_id": identity["chemical_identity"]["chemical_id"],
            "source_identity_id": identity["chemical_identity"][
                "source_identity_id"
            ],
        }
        for field, computed in computed_fields.items():
            stored = sdf_properties.get(field)
            if stored and stored != str(computed):
                raise PoseIdentityError(
                    "docked_sdf_identity_mismatch",
                    f"stored {field} does not match the docked pose at record {i}",
                    details={
                        "source_pose_index": i,
                        "field": field,
                        "stored": stored,
                        "computed": computed,
                    },
                )
        stored_pose_id = sdf_properties.get("pose_id")
        stored_coordinate_signature = sdf_properties.get(
            "coordinate_signature"
        )
        identity_basis = sdf_properties.get("pose_identity_basis")
        serialized_identity_conflicts = {
            field: (stored, computed)
            for field, stored, computed in (
                ("pose_id", stored_pose_id, identity["pose_id"]),
                (
                    "coordinate_signature",
                    stored_coordinate_signature,
                    identity["coordinates"]["signature"],
                ),
            )
            if stored and stored != computed
        }
        if serialized_identity_conflicts:
            raise PoseIdentityError(
                "docked_sdf_identity_mismatch",
                "serialized docked coordinates do not match their stored "
                f"pose identity at record {i}",
                details={
                    "source_pose_index": i,
                    "pose_identity_basis": identity_basis,
                    "conflicts": serialized_identity_conflicts,
                },
            )
        record.update({
            "pose_id": identity["pose_id"],
            "parent_pose_id": identity.get("parent_pose_id"),
            "coordinate_revision": identity["coordinate_revision"],
            "lineage_metadata": identity.get("lineage_metadata") or {},
            "chemical_identity": identity["chemical_identity"],
            "pose_identity_basis": identity_basis or "legacy_or_recomputed",
            "identity_status": {"status": "ok"},
        })
        return record

    records: list[dict] = []
    failures: list[dict] = []
    for i, block in enumerate(blocks):
        try:
            records.append(_pose_record_from_block(i, block))
        except PoseIdentityError as exc:
            # A small backtrack affects an individual candidate. Rejecting one
            # failing pose must not discard the batch, so one
            # unusable block may not discard every other pose in the file.  The
            # pose is kept with a failed ``identity_status`` -- the IDQC hard
            # reject the downstream gate already looks for -- so it is reported
            # and rejected instead of vanishing with the whole run.
            failure = {
                "source_pose_index": i,
                "code": exc.code,
                "reason": str(exc),
            }
            if exc.details:
                failure["details"] = exc.details
            failures.append(failure)
            records.append(_quarantined_pose_record(i, block, docked_sdf, exc))
    if failures and len(failures) == len(blocks):
        # Systemic failure: no record in this artifact is trustworthy.
        # Escalate to protein/ligand preparation rather than a per-pose rejection.
        raise PoseIdentityError(
            "docked_sdf_identity_systemic_failure",
            "every physical record in the docked SDF failed its identity "
            f"contract ({len(failures)} of {len(blocks)})",
            details={"path": str(docked_sdf), "failures": failures},
        )
    return records


def _quarantined_pose_record(
    i: int, block: str, docked_sdf: Path, exc: PoseIdentityError,
) -> dict:
    """Represent a pose whose identity contract failed, without inventing one.

    IDQC failure is one of the few genuine hard rejects, so this
    record deliberately carries no ``pose_id``/``chemical_identity``: fabricating
    either would let the pose join a cohort it may not belong to.  ``mol_name``
    is taken from the block textually, matching what ``_write_selected_sdf``
    reads back, because RDKit may have refused to parse the record at all.
    """
    error = {"code": exc.code, "message": str(exc)}
    if exc.details:
        error["details"] = exc.details
    lines = block.rstrip("\n\r").splitlines()
    physical_name = lines[0].strip() if lines and lines[0].strip() else ""
    try:
        witness = sdf_molblock_sha256(block)
    except PoseIdentityError:
        # An unframed block has no witness to compute.  That is exactly the
        # class of defect being quarantined, so it must not raise again out of
        # the handler and take the rest of the artifact down with it.
        witness = None
    record = {
        "source_pose_index": i,
        "mol_name": physical_name or f"pose_{i}",
        "pose_file": str(docked_sdf),
        "identity_status": {
            "status": "failed",
            "code": exc.code,
            "error": error,
        },
        "identity_qc_failed": True,
        "structural_reject": True,
        "structural_reject_reason": "pose_identity_contract_failed",
    }
    if witness is not None:
        record["_sdf_molblock_sha256"] = witness
    return record


def _source_identity(mol, mol_name: str) -> dict:
    """Bind chemical identity to compound/microstate provenance, not pose rank."""
    if mol.HasProp("source_identity") and mol.GetProp("source_identity").strip():
        try:
            stored = json.loads(mol.GetProp("source_identity"))
        except json.JSONDecodeError as exc:
            raise PoseIdentityError(
                "invalid_source_identity",
                "source_identity SDF property is not valid JSON",
            ) from exc
        if not isinstance(stored, dict) or not stored:
            raise PoseIdentityError(
                "invalid_source_identity",
                "source_identity SDF property must be a non-empty object",
            )
        return stored
    identity = {"compound_id": mol_name}
    aliases = {
        "source": ("source", "generator", "origin"),
        "microstate_id": ("microstate_id", "protonation_state_id", "tautomer_id"),
        "variant_id": ("variant_id", "stereoisomer_id"),
        "input_compound_id": ("input_compound_id", "parent_compound_id"),
    }
    for target, properties in aliases.items():
        for prop in properties:
            if mol.HasProp(prop) and mol.GetProp(prop).strip():
                identity[target] = mol.GetProp(prop).strip()
                break
    return identity


def _strict_integer(value: object) -> int:
    """Parse an index without bool acceptance or lossy float truncation."""
    if isinstance(value, bool):
        raise ValueError("boolean is not an integer index")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped and stripped.lstrip("-").isdigit():
            return int(stripped)
    raise ValueError("value is not a canonical integer")


def _load_docking_scores_contract(path_value: object) -> list[dict]:
    """Treat a declared-but-unreadable score artifact as an operational fault."""
    if path_value in (None, ""):
        return []
    path = Path(str(path_value))
    if not path.is_file():
        raise PoseIdentityError(
            "docking_scores_artifact_missing",
            f"declared docking_scores_file not found: {path}",
            details={"path": str(path)},
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise PoseIdentityError(
            "docking_scores_artifact_unreadable",
            f"could not read docking_scores_file: {path}",
            details={"path": str(path), "reason": str(exc)},
        ) from exc
    except json.JSONDecodeError as exc:
        raise PoseIdentityError(
            "docking_scores_artifact_invalid_json",
            f"docking_scores_file is invalid JSON: {path}",
            details={"path": str(path), "reason": str(exc)},
        ) from exc
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict) and isinstance(payload.get("scores"), list):
        return payload["scores"]
    raise PoseIdentityError(
        "docking_scores_artifact_invalid_shape",
        "docking_scores_file must contain a score list or an object with a "
        "score list",
        details={"path": str(path)},
    )


def _load_verified_docking_inputs(
    report_value: object,
    *,
    docked_path: Path,
    docking_scores_file: object,
    receptor_path: Path,
) -> tuple[dict, list[dict], str]:
    """Bind Step 4 inputs to one committed Step 3 artifact set.

    Matching names and artifact-local indices cannot establish run lineage.
    The Step 3 commit marker seals both files and supplies the set identity that
    the coordinator/report can later compare with this evaluation report.
    """
    report_path = Path(str(report_value)).expanduser().resolve()
    verified = load_committed_docking_artifacts(report_path)
    if not verified.get("available"):
        raise PoseIdentityError(
            "docking_artifact_set_unavailable",
            str(
                verified.get("reason")
                or "docking artifact set is unavailable"
            ),
            details={"docking_report_file": str(report_path)},
        )
    report = verified.get("docking_report") or {}
    contract = report.get("docked_pose_artifact_contract") or {}
    files = report.get("files") or {}
    try:
        schema = int(contract.get("artifact_set_schema_version"))
    except (TypeError, ValueError):
        schema = -1
    if schema < 4:
        raise PoseIdentityError(
            "docking_artifact_set_missing_receptor_binding",
            "docking artifact set does not seal the evaluation receptor; "
            "rerun docking",
        )
    evaluation_receptor = contract.get("evaluation_receptor") or {}
    try:
        receptor_digest = artifact_signature(
            receptor_path, label="evaluation receptor"
        )["digest"]
        receptor_size = Path(receptor_path).stat().st_size
        receptor_matches = bool(
            receptor_digest == evaluation_receptor.get("sha256")
            and not isinstance(evaluation_receptor.get("size_bytes"), bool)
            and receptor_size == int(evaluation_receptor.get("size_bytes"))
        )
    except (KeyError, OSError, TypeError, ValueError):
        receptor_matches = False
    if not receptor_matches:
        raise PoseIdentityError(
            "docking_artifact_set_input_mismatch",
            "receptor_file does not belong to docking_report_file",
            details={"supplied_receptor_file": str(receptor_path)},
        )
    verified_sdf = Path(
        str(verified.get("docked_poses_sdf") or "")
    ).expanduser().resolve()
    supplied_sdf = Path(docked_path).expanduser().resolve()
    if supplied_sdf != verified_sdf:
        raise PoseIdentityError(
            "docking_artifact_set_input_mismatch",
            "docked_sdf does not belong to docking_report_file",
            details={
                "supplied_docked_sdf": str(supplied_sdf),
                "sealed_docked_sdf": str(verified_sdf),
            },
        )
    sealed_scores_value = files.get("docking_scores_json")
    if sealed_scores_value in (None, ""):
        raise PoseIdentityError(
            "docking_artifact_set_input_mismatch",
            "docking report does not declare its sealed score artifact",
        )
    sealed_scores = Path(str(sealed_scores_value)).expanduser().resolve()
    if docking_scores_file not in (None, ""):
        supplied_scores = Path(
            str(docking_scores_file)
        ).expanduser().resolve()
        if supplied_scores != sealed_scores:
            raise PoseIdentityError(
                "docking_artifact_set_input_mismatch",
                "docking_scores_file does not belong to docking_report_file",
                details={
                    "supplied_docking_scores_file": str(supplied_scores),
                    "sealed_docking_scores_file": str(sealed_scores),
                },
            )
    lineage = {
        "status": "verified",
        "artifact_set_schema_version": contract.get(
            "artifact_set_schema_version"
        ),
        "artifact_set_id": contract.get("artifact_set_id"),
        "docking_report_file": str(report_path),
        "docked_sdf": dict(contract.get("docked_sdf") or {}),
        "docking_scores_json": dict(
            contract.get("docking_scores_json") or {}
        ),
        "evaluation_receptor": dict(evaluation_receptor),
        "score_protocol_fingerprint_sha256": contract.get(
            "score_protocol_fingerprint_sha256"
        ),
    }
    return (
        lineage,
        list(verified.get("docking_scores") or []),
        str(sealed_scores),
    )


def _record_identity_field(record: dict, field: str):
    if field in record and record.get(field) not in (None, ""):
        return record[field]
    chemical_identity = record.get("chemical_identity") or {}
    aliases = {
        "chemical_identity_id": "identity_id",
        "chemical_id": "chemical_id",
        "source_identity_id": "source_identity_id",
    }
    nested = aliases.get(field)
    return chemical_identity.get(nested) if nested else None


def _validate_score_record_join(
    record: dict,
    score: dict,
    pose_index: int,
    *,
    all_records: list[dict] | None = None,
) -> None:
    """Require score evidence to resolve uniquely to this physical record."""
    if not score:
        return
    records = list(all_records or [record])
    records_by_index = {
        int(item["source_pose_index"]): item for item in records
    }
    candidate_indices = set(records_by_index)
    evidence_fields: list[str] = []
    strong_identity_matched = False

    def apply_exact_field(field: str, *, integer: bool = False) -> None:
        nonlocal candidate_indices, strong_identity_matched
        score_value = score.get(field)
        record_value = _record_identity_field(record, field)
        if score_value in (None, "") or record_value in (None, ""):
            return
        try:
            normalized_score = (
                _strict_integer(score_value)
                if integer else str(score_value).strip()
            )
            normalized_record = (
                _strict_integer(record_value)
                if integer else str(record_value).strip()
            )
        except ValueError as exc:
            raise PoseIdentityError(
                "invalid_docked_pose_metadata",
                f"invalid {field} at pose_index {pose_index}",
            ) from exc
        if normalized_record != normalized_score:
            raise PoseIdentityError(
                "pose_score_identity_mismatch",
                f"pose record and docking score disagree on {field} at "
                f"pose_index {pose_index}",
                details={
                    "pose_index": pose_index,
                    "field": field,
                    "pose_value": record_value,
                    "score_value": score_value,
                },
            )
        matching = set()
        for candidate_index, candidate in records_by_index.items():
            candidate_value = _record_identity_field(candidate, field)
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
                matching.add(candidate_index)
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
        "engine_pose_index", "engine_pose_rank",
    ):
        apply_exact_field(field, integer=True)

    current_sdf_properties = record.get("_sdf_join_properties") or {}
    stored_local_index = current_sdf_properties.get("pose_index")
    if stored_local_index not in (None, ""):
        try:
            stored_local_index = _strict_integer(stored_local_index)
        except ValueError as exc:
            raise PoseIdentityError(
                "invalid_docked_pose_metadata",
                f"SDF pose_index is invalid at physical record {pose_index}",
            ) from exc
        if stored_local_index != pose_index:
            raise PoseIdentityError(
                "pose_score_identity_mismatch",
                "SDF artifact-local pose_index disagrees with its physical "
                f"record {pose_index}",
            )

    score_molblock_digest = score.get("sdf_molblock_sha256")
    if score_molblock_digest not in (None, ""):
        score_molblock_digest = str(score_molblock_digest).strip()
        record_digest = record.get("_sdf_molblock_sha256")
        if score_molblock_digest != record_digest:
            raise PoseIdentityError(
                "pose_score_identity_mismatch",
                "docking score mol-block witness disagrees with physical "
                f"record at pose_index {pose_index}",
                details={"pose_index": pose_index},
            )
        candidate_indices &= {
            candidate_index
            for candidate_index, candidate in records_by_index.items()
            if candidate.get("_sdf_molblock_sha256") == score_molblock_digest
        }
        evidence_fields.append("sdf_molblock_sha256")
        strong_identity_matched = True

    if score.get("source_pose_index") not in (None, ""):
        try:
            score_source_index = _strict_integer(score["source_pose_index"])
        except ValueError as exc:
            raise PoseIdentityError(
                "invalid_docking_pose_index",
                f"score source_pose_index is invalid at row {pose_index}",
            ) from exc
        stored_source_index = current_sdf_properties.get("source_pose_index")
        if stored_source_index not in (None, ""):
            try:
                stored_source_index = _strict_integer(stored_source_index)
            except ValueError as exc:
                raise PoseIdentityError(
                    "invalid_docked_pose_metadata",
                    "SDF source_pose_index is not a canonical integer at "
                    f"record {pose_index}",
                ) from exc
            if stored_source_index != score_source_index:
                raise PoseIdentityError(
                    "pose_score_identity_mismatch",
                    "SDF and score source_pose_index disagree at local "
                    f"pose_index {pose_index}",
                )
            matching = set()
            for candidate_index, candidate in records_by_index.items():
                candidate_source_index = (
                    candidate.get("_sdf_join_properties") or {}
                ).get("source_pose_index")
                if candidate_source_index in (None, ""):
                    continue
                try:
                    candidate_source_index = _strict_integer(
                        candidate_source_index
                    )
                except ValueError:
                    continue
                if candidate_source_index == score_source_index:
                    matching.add(candidate_index)
            candidate_indices &= matching
            evidence_fields.append("source_pose_index")

    record_names = {
        str(record[field]).strip()
        for field in ("mol_name", "name")
        if record.get(field) not in (None, "")
        and str(record[field]).strip()
    }
    score_names = {
        str(score[field]).strip()
        for field in ("mol_name", "name")
        if score.get(field) not in (None, "")
        and str(score[field]).strip()
    }
    if (
        score.get("compound_id") not in (None, "")
        and _record_identity_field(record, "compound_id") in (None, "")
    ):
        score_names.add(str(score["compound_id"]).strip())
    if score_names:
        if record_names & score_names:
            candidate_indices &= {
                candidate_index
                for candidate_index, candidate in records_by_index.items()
                if {
                    str(candidate[field]).strip()
                    for field in ("mol_name", "name")
                    if candidate.get(field) not in (None, "")
                    and str(candidate[field]).strip()
                } & score_names
            }
            evidence_fields.append("display_name")
        elif not strong_identity_matched:
            raise PoseIdentityError(
                "pose_score_identity_mismatch",
                "pose record and docking score molecule names disagree at "
                f"pose_index {pose_index}",
                details={
                    "pose_index": pose_index,
                    "pose_names": sorted(record_names),
                    "score_names": sorted(score_names),
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
            raise PoseIdentityError(
                "invalid_docking_score_value",
                f"docking score {field} must be a finite JSON number at "
                f"pose_index {pose_index}",
                details={"pose_index": pose_index, "field": field},
            )

    for field, raw_score_value in score.items():
        if isinstance(raw_score_value, bool) or not isinstance(
            raw_score_value, (int, float)
        ):
            continue
        score_value = float(raw_score_value)
        if not math.isfinite(score_value):
            raise PoseIdentityError(
                "non_finite_docking_score",
                f"docking score {field} is non-finite at pose_index "
                f"{pose_index}",
                details={"pose_index": pose_index, "field": field},
            )
        if field in {
            "pose_index", "source_pose_index", "source_molecule_index",
            "engine_pose_index", "engine_pose_rank", "round_id",
        } or field not in current_sdf_properties:
            continue
        try:
            sdf_value = float(current_sdf_properties[field])
        except (TypeError, ValueError):
            continue
        if not (
            math.isfinite(sdf_value)
            and math.isclose(
                sdf_value, score_value, rel_tol=1e-7, abs_tol=1e-7,
            )
        ):
            raise PoseIdentityError(
                "pose_score_value_mismatch",
                f"SDF and docking score {field} disagree at pose_index "
                f"{pose_index}",
                details={
                    "pose_index": pose_index,
                    "field": field,
                    "sdf_value": sdf_value,
                    "score_value": score_value,
                },
            )
        matching = set()
        for candidate_index, candidate in records_by_index.items():
            value = (candidate.get("_sdf_join_properties") or {}).get(field)
            if value in (None, ""):
                continue
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(value) and math.isclose(
                value, score_value, rel_tol=1e-7, abs_tol=1e-7,
            ):
                matching.add(candidate_index)
        candidate_indices &= matching
        evidence_fields.append(field)

    for structured_index, structured in enumerate(score.get("pose_scores") or []):
        if not isinstance(structured, dict):
            continue
        value = structured.get("value")
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
        ):
            raise PoseIdentityError(
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
                raise PoseIdentityError(
                    "invalid_docking_pose_index",
                    "structured score has invalid pose_index",
                ) from exc
            if nested_pose_index != pose_index:
                raise PoseIdentityError(
                    "pose_score_identity_mismatch",
                    "structured score pose_index disagrees with its outer row",
                )
        nested_pose_id = structured.get("pose_id")
        current_pose_id = _record_identity_field(record, "pose_id")
        if nested_pose_id not in (None, "") and current_pose_id not in (None, ""):
            if str(nested_pose_id) != str(current_pose_id):
                raise PoseIdentityError(
                    "pose_score_identity_mismatch",
                    "structured score pose_id disagrees with its pose record",
                )
            candidate_indices &= {
                candidate_index
                for candidate_index, candidate in records_by_index.items()
                if str(_record_identity_field(candidate, "pose_id") or "")
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
                raise PoseIdentityError(
                    "pose_score_value_mismatch",
                    "structured and flat docking scores disagree at "
                    f"pose_index {pose_index}",
                )
            sdf_value = current_sdf_properties.get(flat_field)
            if sdf_value not in (None, ""):
                try:
                    sdf_number = float(sdf_value)
                except (TypeError, ValueError) as exc:
                    raise PoseIdentityError(
                        "invalid_docked_pose_metadata",
                        f"SDF {flat_field} is not numeric at pose_index "
                        f"{pose_index}",
                    ) from exc
                if not (
                    math.isfinite(sdf_number)
                    and math.isclose(
                        sdf_number, float(value),
                        rel_tol=1e-7, abs_tol=1e-7,
                    )
                ):
                    raise PoseIdentityError(
                        "pose_score_value_mismatch",
                        "structured docking score and SDF property disagree "
                        f"at pose_index {pose_index}",
                    )

    if (
        "_sdf_join_properties" in record
        and (not evidence_fields or candidate_indices != {pose_index})
    ):
        raise PoseIdentityError(
            "insufficient_pose_identity_evidence",
            "docking score evidence does not uniquely resolve to physical "
            f"pose_index {pose_index}",
            details={
                "pose_index": pose_index,
                "evidence_fields": evidence_fields,
                "candidate_indices": sorted(candidate_indices),
            },
        )


def _strict_scores_by_pose_index(
    docking_scores: list[dict] | None,
    pose_records: list[dict],
) -> dict[int, dict]:
    record_indices = {int(record["source_pose_index"]) for record in pose_records}
    if len(record_indices) != len(pose_records):
        raise PoseIdentityError(
            "duplicate_source_pose_index",
            "pose records contain duplicate source_pose_index values",
        )
    indexed: dict[int, dict] = {}
    for row_number, score in enumerate(docking_scores or []):
        if not isinstance(score, dict):
            raise PoseIdentityError(
                "invalid_docking_score_record",
                f"docking score row {row_number} is not an object",
            )
        pose_value = score.get("pose_index")
        source_value = score.get("source_pose_index")
        if pose_value in (None, "") and source_value in (None, ""):
            raise PoseIdentityError(
                "missing_docking_pose_index",
                f"docking score row {row_number} has no pose index",
            )
        if pose_value in (None, ""):
            raise PoseIdentityError(
                "missing_docking_pose_index",
                f"docking score row {row_number} has no artifact-local "
                "pose_index",
            )
        try:
            pose_index = _strict_integer(pose_value)
            if source_value not in (None, ""):
                _strict_integer(source_value)
        except ValueError as exc:
            raise PoseIdentityError(
                "invalid_docking_pose_index",
                f"docking score row {row_number} has inconsistent pose indices",
            ) from exc
        if pose_index not in record_indices:
            raise PoseIdentityError(
                "docking_pose_index_out_of_range",
                f"docking score pose_index {pose_index} has no physical SDF record",
                details={"row_number": row_number, "pose_index": pose_index},
            )
        if pose_index in indexed:
            raise PoseIdentityError(
                "duplicate_docking_pose_index",
                f"duplicate docking score pose_index: {pose_index}",
            )
        indexed[pose_index] = score

    records_by_index = {
        int(record["source_pose_index"]): record for record in pose_records
    }
    for pose_index, score in indexed.items():
        record = records_by_index[pose_index]
        identity_status = record.get("identity_status")
        if (
            isinstance(identity_status, dict)
            and identity_status.get("status") != "ok"
        ):
            # The pose already failed its own identity contract and is carried
            # only so it can be reported and rejected.  There is
            # nothing trustworthy left to join against, so re-asserting the join
            # contract here would turn that per-pose reject back into the
            # whole-run abort this path exists to avoid.
            continue
        _validate_score_record_join(
            record, score, pose_index, all_records=pose_records,
        )
    return indexed


def _select_pose_records_from_metadata(
    pose_records: list[dict],
    docking_scores: list[dict] | None,
    selection_strategy: str = "all",
    poses_per_compound: int = 20,
    max_poses: int | None = None,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> list[dict]:
    """Bind scores to identities and choose the planned pose evaluation cohort.

    The pipeline default evaluates every docking-ranked pose within the
    cohort-size and configured per-compound caps. Legacy sequential QC retains
    its backups.
    """
    scores_by_index = _strict_scores_by_pose_index(
        docking_scores, pose_records,
    )
    enriched = []
    for record in pose_records:
        source_index = int(record["source_pose_index"])
        score = scores_by_index.get(source_index, {})
        raw_scores = {
            key: value
            for key, value in score.items()
            if key not in {"pose_index", "source_pose_index", "mol_name", "name"}
        }
        pose_scores = collect_pose_scores({
            **score,
            "pose_index": source_index,
            "pose_id": record.get("pose_id"),
            "stage": score.get("stage", "docking"),
        })
        if pose_scores:
            raw_scores["pose_scores"] = pose_scores
        enriched_record = {
            **record,
            # Public coordinate witness retained for the Step 4 -> Step 6
            # identity join.  The private copy is used while selecting input
            # rows, but dropping it entirely reopened cross-run score joins.
            "sdf_molblock_sha256": record.get("_sdf_molblock_sha256"),
            "docking_scores": raw_scores,
            "pose_scores": pose_scores,
            "_compound_selection_group": _compound_selection_group(
                record, score,
            ),
        }
        enriched_record.pop("_sdf_join_properties", None)
        enriched_record.pop("_sdf_molblock_sha256", None)
        # Score rows are emitted after the docking protocol fingerprint is
        # known, so they supersede any provisional SDF/input cohort label.
        enriched_record.update(_scientific_cohort_for(score, record))
        enriched_record["pre_qc_docking_gate"] = aggregate_pose_gate(
            assess_record(enriched_record, _DOCKING_METRIC_RULES),
            screening_mode=screening_mode,
        )
        enriched.append(enriched_record)

    strategy = (selection_strategy or "all").lower()
    primary_metric = resolve_primary_metric(
        (record.get("pose_scores") or [] for record in enriched),
        preferred_metric_ids=POSE_SELECTION_PREFERENCE,
    )
    for record in enriched:
        record["primary_score_metric_id"] = primary_metric
    if strategy in {"all", "input_order"}:
        selected = list(enriched)
    elif strategy in {"top_by_cnnscore", "top"}:
        selected = sorted(
            enriched,
            key=lambda record: _pose_score_sort_key(record, primary_metric),
        )
    elif strategy in {"ranked_per_compound", "first_passing_per_compound", "best_per_compound", "best_per_molecule", "unique_compounds"}:
        by_name: dict[str, list[dict]] = {}
        for record in enriched:
            name = str(
                record.get("_compound_selection_group")
                or record.get("mol_name")
                or f"pose_{record['source_pose_index']}"
            )
            by_name.setdefault(name, []).append(record)
        selected = []
        sequential = strategy == "first_passing_per_compound"
        adaptive = strategy == "ranked_per_compound"
        limit = None if sequential else max(1, int(poses_per_compound))
        if adaptive:
            limit = min(limit, pose_quota_for_compound_count(len(by_name)))
        for records in by_name.values():
            for compound_rank, record in enumerate(
                sorted(
                    records,
                    key=_docking_pose_rank_key if sequential or adaptive else (
                        lambda item: _pose_score_sort_key(item, primary_metric)
                    ),
                )[:limit],
                start=1,
            ):
                selected.append({**record, "compound_pose_rank": compound_rank})
        if adaptive:
            selected.sort(key=lambda record: (
                record["compound_pose_rank"],
                str(record["_compound_selection_group"]),
            ))
        elif not sequential:
            selected.sort(key=lambda record: _pose_score_sort_key(record, primary_metric))
    else:
        raise ValueError(
            "selection_strategy must be one of: ranked_per_compound, first_passing_per_compound, all, top_by_cnnscore, best_per_compound"
        )

    if max_poses is not None and strategy not in {"ranked_per_compound", "first_passing_per_compound"}:
        selected = selected[:max(0, max_poses)]
    for rank, record in enumerate(selected, start=1):
        record["selection_rank"] = rank
        record["selection_strategy"] = strategy
    return selected


def _docking_pose_rank_key(record: dict) -> tuple:
    rank = _number(record.get("engine_pose_rank"), None)
    if rank is None:
        rank = _number((record.get("docking_scores") or {}).get("engine_pose_rank"), None)
    return (
        rank if rank is not None and rank > 0 else float("inf"),
        _pose_score_sort_key(record),
    )


def _cached_batch_sizes(batches_dir: Path) -> list[int]:
    """Preserve prior batch boundaries when extending a partial evaluation."""
    sizes = []
    while True:
        report_file = batches_dir / f"batch_{len(sizes) + 1:04d}" / "evaluation_report.json"
        if not report_file.is_file():
            return sizes
        try:
            rows = json.loads(report_file.read_text()).get("evaluation")
        except (OSError, ValueError, TypeError):
            return sizes
        if not isinstance(rows, list) or not rows:
            return sizes
        sizes.append(len(rows))


def _evaluation_batches(
    poses, evaluation, batch_size, *, sequential=False, dynamic_pending=False,
    resume_batch_sizes=(),
):
    """Yield pose batches while preserving the complete planned inventory.

    The caller appends each completed/resumed batch to ``evaluation`` before
    advancing this iterator. Adaptive budgets recompute pending poses after a
    partially evaluated batch; legacy sequential QC advances only on a fail.
    """
    if dynamic_pending:
        batch_index = 0
        while True:
            observed = {
                int(row.get("source_pose_index", row.get("pose_index", -1)))
                for row in evaluation
            }
            pending = [
                pose for pose in poses
                if int(pose["source_pose_index"]) not in observed
            ]
            if not pending:
                return
            size = (
                resume_batch_sizes[batch_index]
                if batch_index < len(resume_batch_sizes) else batch_size
            )
            batch_index += 1
            yield pending[:size]
    if not sequential:
        yield from _chunks(poses, batch_size)
        return
    groups = {}
    for pose in poses:
        name = pose.get("_compound_selection_group") or pose.get("mol_name")
        groups.setdefault(name, []).append(pose)
    while True:
        observed = {
            int(row.get("source_pose_index", row.get("pose_index", -1))): row
            for row in evaluation
        }
        pending = []
        for records in groups.values():
            for pose in records:
                row = observed.get(int(pose["source_pose_index"]))
                if row is None:
                    pending.append(pose)
                    break
                # A tool fault or missing evidence is not a scientific rejection.
                if row.get("verdict") != "fail" or "tool_failed" in (
                    row.get("evidence_status") or {}
                ).values():
                    break
        if not pending:
            return
        yield from _chunks(pending, batch_size)


def _compound_selection_group(record: dict, score: dict) -> str:
    """Recover grouping identity from the index-matched score for bad SDF rows.

    A parse failure deliberately keeps its failed identity verdict, but its
    synthetic ``pose_N`` display name must not turn every pose into a separate
    compound.  Docking score rows are already joined by source pose index and
    retain the engine's compound name, which is sufficient for shortlist
    grouping without pretending the molecule itself parsed successfully.
    """
    identity_status = record.get("identity_status")
    if isinstance(identity_status, dict) and identity_status.get("status") != "ok":
        score_name = str(
            score.get("compound_id")
            or score.get("mol_name")
            or score.get("name")
            or ""
        ).strip()
        if score_name:
            return score_name
    return str(
        record.get("mol_name") or f"pose_{record['source_pose_index']}"
    )


def _pose_score_sort_key(
    record: dict,
    primary_metric_id: str | None = None,
) -> tuple:
    primary_metric = (
        primary_metric_id
        or record.get("primary_score_metric_id")
    )
    primary_score = find_score(
        collect_pose_scores(record, record.get("docking_scores") or {}),
        primary_metric,
    )
    if primary_metric:
        return (
            _docking_gate_order(record),
            score_sort_key(primary_score),
            int(record.get("source_pose_index", 0)),
        )
    scores = record.get("docking_scores") or {}
    return (
        _docking_gate_order(record),
        -_number(scores.get("CNNscore"), default=-1.0),
        -_number(scores.get("CNNaffinity"), default=-1.0),
        _number(scores.get("minimizedAffinity"), default=float("inf")),
        _number(scores.get("docking_score"), default=float("inf")),
        _number(scores.get("vina_affinity"), default=float("inf")),
        int(record.get("source_pose_index", 0)),
    )


def _docking_gate_order(record: dict) -> int:
    """Rank a pre-QC docking gate, treating an unknown gate as unmeasured."""
    gate = record.get("pre_qc_docking_gate") or "not_evaluated"
    return _DOCKING_GATE_ORDER.get(gate, _DOCKING_GATE_ORDER["not_evaluated"])


def _post_qc_survivor_pose_records(
    evaluation: list[dict],
    evaluation_input_poses: list[dict],
) -> list[dict]:
    """Return every surviving pose in deterministic evaluated-rank order.

    Survivors are pass, warn, and not_evaluated: a pose no tool could measure is
    carried forward with its verdict attached rather than deleted.
    """
    input_by_index = {
        int(record["source_pose_index"]): record
        for record in evaluation_input_poses
    }
    survivors: list[dict] = []
    for entry in sorted(evaluation, key=_evaluated_pose_sort_key):
        if entry.get("verdict") not in _POSE_SURVIVOR_VERDICTS:
            continue
        source_index = int(
            entry.get("source_pose_index", entry.get("pose_index", -1))
        )
        source_record = input_by_index.get(source_index)
        if source_record is None:
            continue
        survivors.append({
            **source_record,
            "post_qc_verdict": entry.get("verdict"),
            "pose_scores": collect_pose_scores(
                entry, entry.get("docking_scores") or {}
            ),
            "primary_score_metric_id": entry.get(
                "primary_score_metric_id"
            ),
        })
    for rank, record in enumerate(survivors, start=1):
        record["selection_rank"] = rank
        record["selection_strategy"] = "all_post_qc_survivors"
    return survivors


def _select_post_qc_pose_records(
    evaluation: list[dict],
    evaluation_input_poses: list[dict],
) -> list[dict]:
    """Select one QC-eligible pose per compound after all available evidence.

    Failed poses remain in ``evaluation_report.json`` for diagnosis but are not
    exported as downstream ``selected_poses``.  A pose that no enabled tool could
    measure is not a failed pose: it is carried with ``not_evaluated`` attached
    and ordered behind every measured pose.
    """
    input_by_index = {
        int(record["source_pose_index"]): record
        for record in evaluation_input_poses
    }
    grouped: dict[str, list[dict]] = {}
    for entry in evaluation:
        if entry.get("verdict") not in _POSE_SURVIVOR_VERDICTS:
            continue
        source_index = int(entry.get("source_pose_index", entry.get("pose_index", -1)))
        if source_index not in input_by_index:
            continue
        name = str(
            entry.get("mol_name")
            or input_by_index[source_index].get("mol_name")
            or f"pose_{source_index}"
        )
        grouped.setdefault(name, []).append(entry)

    selected: list[dict] = []
    for mol_name, entries in grouped.items():
        best = min(entries, key=_evaluated_pose_sort_key)
        source_index = int(best.get("source_pose_index", best.get("pose_index")))
        source_record = input_by_index[source_index]
        selected.append({
            **source_record,
            "mol_name": mol_name,
            "post_qc_verdict": best.get("verdict"),
            "post_qc_desired_score": _desired_ranking_score(
                best.get("blueprint_match") or {}
            ),
            "post_qc_selection": True,
            "pose_scores": collect_pose_scores(
                best,
                best.get("docking_scores") or {},
            ),
            "primary_score_metric_id": best.get(
                "primary_score_metric_id"
            ),
        })

    selected.sort(key=lambda record: (
        _EVALUATED_VERDICT_ORDER.get(record.get("post_qc_verdict"), 3),
        -_number(record.get("post_qc_desired_score"), 0.0),
        *_pose_score_sort_key(record),
    ))
    for rank, record in enumerate(selected, start=1):
        record["selection_rank"] = rank
        record["selection_strategy"] = "post_qc_best_per_compound"
    return selected


def _pose_exhausted_compounds(
    evaluation: list[dict],
    evaluation_input_poses: list[dict],
) -> list[dict]:
    """Report compounds only after every planned pose was evaluated and rejected.

    Exhaustion skips the molecule under the current protocol without redocking.
    Keep a separate diagnostic because it has no surviving pose to export.
    Missing verdicts, tool faults and unreadable structures are not evidence
    that every pose scientifically failed.
    """
    input_by_index = {
        int(record["source_pose_index"]): record
        for record in evaluation_input_poses
    }
    expected_indices_by_name: dict[str, set[int]] = {}
    for source_index, record in input_by_index.items():
        name = str(record.get("mol_name") or f"pose_{source_index}")
        expected_indices_by_name.setdefault(name, set()).add(source_index)
    grouped: dict[str, list[dict]] = {}
    for entry in evaluation:
        source_index = int(
            entry.get("source_pose_index", entry.get("pose_index", -1))
        )
        if source_index not in input_by_index:
            continue
        name = str(
            entry.get("mol_name")
            or input_by_index[source_index].get("mol_name")
            or f"pose_{source_index}"
        )
        grouped.setdefault(name, []).append(entry)

    exhausted: list[dict] = []
    for mol_name, entries in grouped.items():
        # Partial checkpoints contain the complete input pose inventory but only
        # the rows evaluated so far.  Do not mistake an all-fail prefix for an
        # exhausted molecule: every planned input pose must have a verdict first.
        observed_indices = {
            int(entry.get("source_pose_index", entry.get("pose_index", -1)))
            for entry in entries
        }
        expected_indices = expected_indices_by_name.get(mol_name)
        if not expected_indices or observed_indices != expected_indices:
            continue
        if any(
            entry.get("verdict") != "fail"
            or any(
                state == "tool_failed" and tool != "prolif"
                for tool, state in (entry.get("evidence_status") or {}).items()
            )
            or entry.get("operational_errors")
            or (
                isinstance(entry.get("identity_status"), dict)
                and entry["identity_status"].get("status") != "ok"
            )
            for entry in entries
        ):
            continue
        attempt = max(
            _docking_attempt_number(
                entry,
                input_by_index[
                    int(entry.get("source_pose_index", entry.get("pose_index", -1)))
                ],
            )
            for entry in entries
        )
        reasons: list[str] = []
        for entry in entries:
            verdict_reasons = entry.get("verdict_reasons") or {}
            reasons.extend(str(item) for item in verdict_reasons.get("hard") or [])
            reasons.extend(str(item) for item in verdict_reasons.get("fail") or [])
        first = entries[0]
        exhausted.append({
            "mol_name": mol_name,
            "compound_id": first.get("compound_id"),
            "chemical_identity_id": (
                (first.get("chemical_identity") or {}).get("identity_id")
            ),
            "cohort_id": first.get("cohort_id"),
            "evaluated_poses": len(entries),
            "pose_verdicts": sorted({
                str(entry.get("verdict")) for entry in entries
            }),
            "reject_reasons": list(dict.fromkeys(reasons)),
            "docking_attempt": attempt,
            "route": "reject_molecule",
            "reason_code": "all_poses_failed",
            "reason": "All selected poses failed QC under the current protocol; molecule skipped without redocking.",
            "route_basis": "pose ensemble exhaustion policy",
        })
    exhausted.sort(key=lambda item: str(item["mol_name"]))
    return exhausted


def _docking_attempt_number(entry: dict, source_record: dict | None = None) -> int:
    """How many docking passes this pose's molecule has had; 1 when undeclared.

    Retained as provenance for historical or explicitly requested docking runs;
    it does not control pose-exhaustion routing.
    """
    for source in (entry, source_record or {}):
        for key in ("docking_attempt", "docking_attempt_number"):
            value = source.get(key)
            if value in (None, ""):
                continue
            try:
                return max(1, int(value))
            except (TypeError, ValueError):
                continue
    return 1


def _merge_rtmscore_pose_scores(
    evaluation: list[dict],
    rtmscore_result: dict,
    *,
    use_as_primary_metric: bool,
) -> int:
    """Join RTMScore evidence by canonical pose_id without score arithmetic."""
    scores = rtmscore_result.get("pose_scores") or []
    if not isinstance(scores, list):
        raise ValueError("RTMScore pose_scores must be a list")
    evaluation_by_pose: dict[str, dict] = {}
    for entry in evaluation:
        pose_id = str(entry.get("pose_id") or "")
        if not pose_id:
            continue
        if pose_id in evaluation_by_pose:
            raise ValueError(f"duplicate pose_id in evaluation cohort: {pose_id}")
        evaluation_by_pose[pose_id] = entry

    expected_metric = "rtmscore.ensemble_score"
    seen: set[str] = set()
    staged: list[tuple[dict, dict]] = []
    for raw_score in scores:
        if not isinstance(raw_score, dict):
            raise ValueError("RTMScore emitted a non-object pose score")
        score = dict(raw_score)
        pose_id = str(score.get("pose_id") or "")
        if not pose_id or pose_id in seen:
            raise ValueError(f"missing or duplicate RTMScore pose_id: {pose_id!r}")
        seen.add(pose_id)
        entry = evaluation_by_pose.get(pose_id)
        if entry is None:
            raise ValueError(
                f"RTMScore score pose_id is absent from evaluation cohort: {pose_id}"
            )
        if score.get("metric_id") != expected_metric or not is_rankable_score(score):
            raise ValueError(
                f"RTMScore emitted an invalid {expected_metric!r} score for {pose_id}"
            )
        score["pose_index"] = entry.get(
            "source_pose_index", entry.get("pose_index")
        )
        staged.append((entry, score))

    # Commit only after the complete result has passed validation.  This keeps
    # a late malformed score from leaving a partially rescored cohort behind.
    for entry, score in staged:
        pose_id = str(score["pose_id"])
        existing = [
            item
            for item in collect_pose_scores(
                entry, entry.get("docking_scores") or {}
            )
            if not (
                item.get("metric_id") == expected_metric
                and item.get("pose_id") == pose_id
            )
        ]
        entry["pose_scores"] = [*existing, score]
        entry.setdefault("docking_scores", {})["pose_scores"] = list(
            entry["pose_scores"]
        )
    if use_as_primary_metric and staged:
        for entry, _score in staged:
            if entry.get("verdict") in {"pass", "warn"}:
                entry["primary_score_metric_id"] = expected_metric
    return len(scores)


def _merge_mmgbsa_pose_scores(
    evaluation: list[dict],
    mmgbsa_result: dict,
    *,
    use_as_primary_metric: bool,
) -> int:
    """Join MM-GBSA scores by exact canonical pose_id without score mixing."""
    scores = mmgbsa_result.get("pose_scores") or []
    if not isinstance(scores, list):
        raise ValueError("MM-GBSA pose_scores must be a list")

    evaluation_by_pose: dict[str, dict] = {}
    for entry in evaluation:
        pose_id = str(entry.get("pose_id") or "")
        if not pose_id:
            continue
        if pose_id in evaluation_by_pose:
            raise ValueError(
                f"duplicate pose_id in evaluation cohort: {pose_id}"
            )
        evaluation_by_pose[pose_id] = entry

    expected_metric = "amber_mmpbsa.delta_total_no_entropy"
    seen: set[str] = set()
    for raw_score in scores:
        if not isinstance(raw_score, dict):
            raise ValueError("MM-GBSA emitted a non-object pose score")
        score = dict(raw_score)
        pose_id = str(score.get("pose_id") or "")
        if not pose_id:
            raise ValueError("MM-GBSA score is missing canonical pose_id")
        if pose_id in seen:
            raise ValueError(f"duplicate MM-GBSA score for pose_id {pose_id}")
        seen.add(pose_id)
        if pose_id not in evaluation_by_pose:
            raise ValueError(
                "MM-GBSA score pose_id is absent from the evaluated cohort: "
                f"{pose_id}"
            )
        if score.get("metric_id") != expected_metric:
            raise ValueError(
                "MM-GBSA score metric_id must be "
                f"{expected_metric!r}, got {score.get('metric_id')!r}"
            )
        if not is_rankable_score(score):
            raise ValueError(
                f"MM-GBSA emitted a non-rankable score for {pose_id}"
            )

        entry = evaluation_by_pose[pose_id]
        source_index = entry.get(
            "source_pose_index", entry.get("pose_index")
        )
        score["pose_index"] = source_index
        existing = [
            item
            for item in collect_pose_scores(
                entry,
                entry.get("docking_scores") or {},
            )
            if not (
                item.get("metric_id") == expected_metric
                and item.get("pose_id") == pose_id
            )
        ]
        entry["pose_scores"] = [*existing, score]
        entry.setdefault("docking_scores", {})["pose_scores"] = list(
            entry["pose_scores"]
        )

    if use_as_primary_metric and scores:
        # One cohort has one primary metric. Poses outside top-K intentionally
        # have a missing value and therefore rank after scored survivors; no
        # arithmetic combination with docking scores is performed.
        for entry in evaluation:
            if entry.get("verdict") in {"pass", "warn"}:
                entry["primary_score_metric_id"] = expected_metric
    return len(scores)


def _refresh_existing_post_qc_ranking(
    *,
    docked_sdf: Path,
    selected_poses: list[dict],
    evaluation: list[dict],
    selected_sdf: Path,
    selection_manifest: Path,
) -> list[dict]:
    """Update/rerank the existing survivor set without changing RBFE nodes."""
    evaluation_by_pose = {
        str(entry.get("pose_id")): entry
        for entry in evaluation
        if entry.get("pose_id")
    }
    updated = []
    for pose in selected_poses:
        entry = evaluation_by_pose.get(str(pose.get("pose_id")))
        if entry is None:
            raise ValueError(
                "post-QC selected pose is absent from evaluation cohort: "
                f"{pose.get('pose_id')}"
            )
        updated.append({
            **pose,
            "pose_scores": collect_pose_scores(
                entry,
                entry.get("docking_scores") or {},
            ),
            "primary_score_metric_id": entry.get(
                "primary_score_metric_id"
            ),
        })

    updated.sort(key=lambda record: (
        0 if record.get("post_qc_verdict") == "pass" else 1,
        -_number(record.get("post_qc_desired_score"), -1.0),
        *_pose_score_sort_key(record),
    ))
    for rank, record in enumerate(updated, start=1):
        record["selection_rank"] = rank
    _write_selected_sdf(docked_sdf, updated, selected_sdf)
    _write_selection_manifest(selection_manifest, updated)
    return updated


def _evaluated_pose_sort_key(entry: dict) -> tuple:
    """Order eligible poses without allowing a score to override QC."""
    record = {
        "source_pose_index": entry.get("source_pose_index", entry.get("pose_index", 0)),
        "docking_scores": entry.get("docking_scores") or {},
        "pose_scores": entry.get("pose_scores") or [],
        "primary_score_metric_id": entry.get("primary_score_metric_id"),
    }
    return (
        _EVALUATED_VERDICT_ORDER.get(entry.get("verdict"), 3),
        -_desired_ranking_score(entry.get("blueprint_match") or {}),
        *_pose_score_sort_key(record),
    )


def _refresh_post_qc_selection(
    docked_sdf: Path,
    evaluation_input_poses: list[dict],
    evaluation: list[dict],
    selected_sdf: Path,
    selection_manifest: Path,
    *,
    run_cache: dict | None = None,
) -> list[dict]:
    selected = _select_post_qc_pose_records(evaluation, evaluation_input_poses)
    _write_selected_sdf(docked_sdf, selected, selected_sdf, run_cache=run_cache)
    _write_selection_manifest(selection_manifest, selected)
    return selected


def _read_sdf_blocks(sdf_path: Path) -> list[str]:
    blocks = []
    current = []
    saw_delimiter = False
    with open(sdf_path, errors="ignore") as handle:
        for line in handle:
            if line.rstrip("\n\r") == "$$$$":
                saw_delimiter = True
                # One delimiter means one physical record, including an empty
                # malformed record.  Retaining it lets the parser fail at that
                # index instead of assigning later scores to earlier poses.
                blocks.append("".join(current))
                current = []
            else:
                current.append(line)
    if current and any(item.strip() for item in current):
        raise PoseIdentityError(
            "unterminated_sdf_record",
            "docked SDF has data after its last $$$$ delimiter",
            details={"path": str(sdf_path)},
        )
    if not saw_delimiter:
        return []
    return blocks


def _without_sdf_properties(lines: list[str], names: set[str]) -> list[str]:
    """Remove named SDF property blocks before writing authoritative values."""
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
            index += 1
            while index < len(lines) and lines[index].strip():
                index += 1
            if index < len(lines):
                index += 1
            continue
        kept.append(line)
        index += 1
    return kept


def _write_selected_sdf(
    source_sdf: Path, selected_poses: list[dict], output_sdf: Path,
    *, run_cache: dict | None = None,
) -> None:
    blocks = cached_file(run_cache, source_sdf, "sdf_blocks", _read_sdf_blocks)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    selected_blocks: list[tuple[dict, int, list[str]]] = []
    for pose in selected_poses:
        source_index = int(pose["source_pose_index"])
        if source_index < 0 or source_index >= len(blocks):
            raise PoseIdentityError(
                "source_pose_index_out_of_range",
                f"source_pose_index out of range: {source_index}",
            )
        lines = blocks[source_index].rstrip("\n\r").splitlines()
        physical_name = (
            lines[0].strip() if lines and lines[0].strip()
            else f"pose_{source_index}"
        )
        selected_name = str(pose.get("mol_name") or "").strip()
        if selected_name and selected_name != physical_name:
            raise PoseIdentityError(
                "pose_export_identity_mismatch",
                "selected pose name does not match its physical source block",
                details={
                    "source_pose_index": source_index,
                    "selected_mol_name": selected_name,
                    "physical_mol_name": physical_name,
                },
            )
        selected_blocks.append((pose, source_index, lines))

    temporary = output_sdf.with_name(output_sdf.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            for local_index, (pose, source_index, lines) in enumerate(
                selected_blocks
            ):
                lines = _without_sdf_properties(
                    lines,
                    {
                        *SCIENTIFIC_COHORT_OUTPUT_FIELDS,
                        "pose_index",
                        "source_pose_index",
                        "selection_rank",
                        "pose_id",
                        "parent_pose_id",
                        "coordinate_revision",
                        "chemical_identity_id",
                        "chemical_id",
                        "source_identity_id",
                        "source_identity",
                        "lineage_metadata",
                    },
                )
                while lines and not lines[-1].strip():
                    lines.pop()
                # An SD property value must be terminated by a blank line.
                # Removing all trailing blanks and appending the next header
                # directly made RDKit consume ``> <pose_index>`` as part of the
                # preceding value, silently corrupting both fields.
                if lines:
                    lines.append("")
                lines.extend([
                    "> <pose_index>",
                    str(local_index),
                    "",
                    "> <source_pose_index>",
                    str(source_index),
                    "",
                    "> <selection_rank>",
                    str(pose.get("selection_rank", "")),
                    "",
                ])
                identity_properties = {
                    "pose_id": pose.get("pose_id"),
                    "parent_pose_id": pose.get("parent_pose_id"),
                    "coordinate_revision": pose.get("coordinate_revision"),
                    "chemical_identity_id": (
                        (pose.get("chemical_identity") or {}).get("identity_id")
                    ),
                    "chemical_id": (
                        (pose.get("chemical_identity") or {}).get("chemical_id")
                    ),
                    "source_identity_id": (
                        (pose.get("chemical_identity") or {}).get(
                            "source_identity_id"
                        )
                    ),
                    "source_identity": json.dumps(
                        (pose.get("chemical_identity") or {}).get(
                            "source_identity"
                        ),
                        sort_keys=True,
                        separators=(",", ":"),
                        ensure_ascii=False,
                    )
                    if (pose.get("chemical_identity") or {}).get(
                        "source_identity"
                    ) is not None
                    else None,
                    "lineage_metadata": json.dumps(
                        pose.get("lineage_metadata") or {},
                        sort_keys=True,
                        separators=(",", ":"),
                        ensure_ascii=False,
                    ),
                    **{
                        field: (
                            json.dumps(
                                pose.get(field),
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                            if isinstance(pose.get(field), (dict, list, bool))
                            else pose.get(field)
                        )
                        for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS
                    },
                }
                for key, value in identity_properties.items():
                    if value is None:
                        continue
                    if isinstance(value, (dict, list, bool)):
                        value = json.dumps(
                            value, sort_keys=True, separators=(",", ":")
                        )
                    lines.extend([f"> <{key}>", str(value), ""])
                handle.write("\n".join(lines))
                handle.write("\n$$$$\n")
        temporary.replace(output_sdf)
    finally:
        if temporary.exists():
            temporary.unlink()


def _write_selection_manifest(path: Path, selected_poses: list[dict]) -> None:
    fields = [
        "selection_rank",
        "compound_pose_rank",
        "source_pose_index",
        "pose_id",
        "parent_pose_id",
        "mol_name",
        "selection_strategy",
        "pre_qc_docking_gate",
        "post_qc_verdict",
        "post_qc_desired_score",
        "primary_score_metric_id",
        "primary_score_value",
        "primary_score_unit",
        "primary_score_direction",
        *SCIENTIFIC_COHORT_OUTPUT_FIELDS,
        "CNNscore",
        "CNNaffinity",
        "minimizedAffinity",
        "docking_score",
        "vina_affinity",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for pose in selected_poses:
            scores = pose.get("docking_scores") or {}
            pose_scores = collect_pose_scores(pose, scores)
            primary_metric = pose.get("primary_score_metric_id")
            primary_score = find_score(pose_scores, primary_metric)
            writer.writerow({
                "selection_rank": pose.get("selection_rank"),
                "compound_pose_rank": pose.get("compound_pose_rank"),
                "source_pose_index": pose.get("source_pose_index"),
                "pose_id": pose.get("pose_id"),
                "parent_pose_id": pose.get("parent_pose_id"),
                "mol_name": pose.get("mol_name"),
                "selection_strategy": pose.get("selection_strategy"),
                "pre_qc_docking_gate": pose.get("pre_qc_docking_gate"),
                "post_qc_verdict": pose.get("post_qc_verdict"),
                "post_qc_desired_score": pose.get("post_qc_desired_score"),
                "primary_score_metric_id": primary_metric,
                "primary_score_value": (
                    primary_score.get("value") if primary_score else None
                ),
                "primary_score_unit": (
                    primary_score.get("unit") if primary_score else None
                ),
                "primary_score_direction": (
                    primary_score.get("direction") if primary_score else None
                ),
                **{
                    field: (
                        json.dumps(
                            pose.get(field),
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        if isinstance(pose.get(field), (dict, list, bool))
                        else pose.get(field)
                    )
                    for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS
                },
                "CNNscore": scores.get("CNNscore"),
                "CNNaffinity": scores.get("CNNaffinity"),
                "minimizedAffinity": scores.get("minimizedAffinity"),
                "docking_score": scores.get("docking_score"),
                "vina_affinity": scores.get("vina_affinity"),
            })


def _scores_for_selected_batch(batch_poses: list[dict]) -> list[dict]:
    scores = []
    for local_index, pose in enumerate(batch_poses):
        entry = {
            "pose_index": local_index,
            "source_pose_index": pose.get("source_pose_index"),
            "pose_id": pose.get("pose_id"),
            "mol_name": pose.get("mol_name"),
        }
        score_payload = dict(pose.get("docking_scores") or {})
        structured = score_payload.pop("pose_scores", pose.get("pose_scores") or [])
        entry.update(score_payload)
        if structured:
            entry["pose_scores"] = [
                {
                    **score,
                    "pose_index": local_index,
                    "pose_id": pose.get("pose_id") or score.get("pose_id"),
                    "source_pose_index": pose.get("source_pose_index"),
                }
                for score in structured
            ]
        scores.append(entry)
    return scores


def _remap_batch_entries(
    entries: list[dict],
    batch_poses: list[dict],
    *,
    posecheck_thresholds: dict | None = None,
    identity_qc_enabled: bool = True,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> list[dict]:
    remapped = []
    for position, entry in enumerate(entries):
        local_index = int(entry.get("pose_index", position))
        pose = batch_poses[local_index] if 0 <= local_index < len(batch_poses) else {}
        source_index = int(pose.get("source_pose_index", entry.get("pose_index", position)))
        updated = dict(entry)
        updated["batch_pose_index"] = local_index
        updated["source_pose_index"] = source_index
        updated["pose_index"] = source_index
        updated["selection_rank"] = pose.get("selection_rank")
        updated["compound_pose_rank"] = pose.get("compound_pose_rank")
        updated["pre_qc_docking_gate"] = pose.get("pre_qc_docking_gate")
        updated["pose_id"] = pose.get("pose_id")
        updated["parent_pose_id"] = pose.get("parent_pose_id")
        updated["chemical_identity"] = pose.get("chemical_identity")
        updated["identity_status"] = pose.get("identity_status")
        for field in (
            "coordinate_revision",
            "lineage_metadata",
            "sdf_molblock_sha256",
            "compound_id",
            "source_molecule_index",
            "preparation_variant_index",
            "engine_pose_index",
            "engine_pose_rank",
            "selection_strategy",
            "primary_score_metric_id",
        ):
            if field in pose:
                updated[field] = pose[field]
            else:
                updated.pop(field, None)
        for field in LIGAND_LINEAGE_PROPERTIES:
            if pose.get(field) not in (None, ""):
                updated[field] = pose[field]
            else:
                updated.pop(field, None)
        updated.setdefault("evidence_status", {})["identity_qc"] = (
            "available"
            if (pose.get("identity_status") or {}).get("status") == "ok"
            else "failed"
            if identity_qc_enabled
            else "disabled_by_config"
        )
        if pose.get("mol_name"):
            updated["mol_name"] = pose["mol_name"]
        updated["docking_scores"] = dict(pose.get("docking_scores") or {})
        updated.update(_scientific_cohort_for(
            pose.get("docking_scores") or {},
            pose,
            updated.get("docking_scores") or {},
            updated,
        ))
        updated["pose_scores"] = list(pose.get("pose_scores") or [])
        if updated["pose_scores"]:
            updated["docking_scores"]["pose_scores"] = list(
                updated["pose_scores"]
            )
        if isinstance(updated.get("blueprint_match"), dict):
            updated["blueprint_match"] = {
                **updated["blueprint_match"],
                "batch_pose_index": local_index,
                "pose_index": source_index,
            }
        updated["verdict"] = _compute_verdict(
            updated,
            posecheck_thresholds,
            screening_mode=screening_mode,
        )
        if (
            identity_qc_enabled
            and (pose.get("identity_status") or {}).get("status") != "ok"
        ):
            identity_error = (
                ((pose.get("identity_status") or {}).get("error") or {}).get("message")
                or "stable pose identity unavailable"
            )
            updated.setdefault("diagnosis", []).append(
                f"IDENTITY_FAIL: {identity_error}"
            )
            updated["recommendation"] = "reject"
            updated["modification_hint"] = (
                "Repair molecule parsing, stereochemistry, coordinates, or "
                "source identity before rescoring"
            )
        remapped.append(updated)
    return remapped


def _inject_smiles(entries: list[dict], smiles_map: dict[str, str]) -> None:
    """Inject SMILES into evaluation entries by matching mol_name to compound_id."""
    for entry in entries:
        name = entry.get("mol_name", "")
        entry["smiles"] = smiles_map.get(name)


def _generate_top_candidates(evaluation: list[dict], top_n: int = 10) -> list[dict]:
    """Generate a Top-N candidate summary from evaluation results."""
    # "not_evaluated" joins pass/warn here: the sort below already places it last,
    # so an unmeasured pose is visible with its verdict rather than deleted from
    # the summary.
    passing = [e for e in evaluation if e.get("verdict") in _POSE_SURVIVOR_VERDICTS]
    primary_metric = resolve_primary_metric(
        (collect_pose_scores(entry, entry.get("docking_scores") or {}) for entry in passing),
        preferred_metric_ids=POSE_SELECTION_PREFERENCE,
    )
    # Hard-QC verdict first, then desired Blueprint satisfaction as a visible
    # tie-breaker before docking scores.
    passing.sort(key=lambda e: (
        _EVALUATED_VERDICT_ORDER.get(e.get("verdict"), 3),
        -_desired_ranking_score(e.get("blueprint_match") or {}),
        score_sort_key(find_score(
            collect_pose_scores(e, e.get("docking_scores") or {}),
            primary_metric,
        )),
    ))
    # Deduplicate by mol_name (keep best pose per compound)
    seen = set()
    candidates = []
    for entry in passing:
        name = entry.get("mol_name") or f"pose_{entry.get('pose_index', len(candidates))}"
        if name in seen:
            continue
        seen.add(name)
        scores = entry.get("docking_scores") or {}
        pose_scores = collect_pose_scores(entry, scores)
        blueprint_match = entry.get("blueprint_match") or {}
        candidates.append({
            "rank": len(candidates) + 1,
            "mol_name": name,
            "pose_index": entry.get("pose_index"),
            "source_pose_index": entry.get(
                "source_pose_index", entry.get("pose_index")
            ),
            "pose_id": entry.get("pose_id"),
            "chemical_identity_id": (
                (entry.get("chemical_identity") or {}).get("identity_id")
            ),
            "smiles": entry.get("smiles"),
            "CNNscore": scores.get("CNNscore"),
            "CNNaffinity": scores.get("CNNaffinity"),
            "minimizedAffinity": scores.get("minimizedAffinity"),
            "pose_scores": pose_scores,
            "primary_score_metric_id": primary_metric,
            "primary_score": find_score(pose_scores, primary_metric),
            "desired_blueprint_score": blueprint_match.get("desired_weighted_score"),
            "desired_blueprint_coverage": blueprint_match.get("desired_evaluable_coverage"),
            "desired_blueprint_ranking_score": blueprint_match.get("desired_ranking_score"),
            "verdict": entry.get("verdict"),
            **{
                field: entry.get(field)
                for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS
            },
        })
        if len(candidates) >= top_n:
            break
    return candidates


def _desired_ranking_score(blueprint_match: dict) -> float:
    value = blueprint_match.get("desired_ranking_score")
    if value is None:
        value = blueprint_match.get("desired_weighted_score")
    # blueprint_match.py computes this as ``weight_hit / all_weight``, so 0.0 is
    # the legal floor.  Falling back to -1.0 placed a pose with no blueprint
    # evidence below every measured pose, letting an absent profiler decide the
    # ``best_per_compound`` pick and the ``max_poses`` cut.  Landing on the floor
    # instead ties it with a measured-zero pose so the next key -- the docking
    # score -- decides, which is what the executable policy asks of a missing term.
    return _number(value, 0.0)


def _posecheck_num_conformers(input_data: dict) -> int:
    tool = (input_data.get("tools") or {}).get("posecheck") or {}
    params = tool.get("params") or {}
    return input_data.get("posecheck_num_conformers", params.get("num_conformers", 10))


def _calibrate_posecheck_thresholds(
    receptor_path: Path,
    reference_ligand_file: str,
    posecheck_thresholds: dict,
    baseline_limits: dict,
    steps_completed: list,
    steps_skipped: list,
    *,
    num_conformers: int = 10,
) -> tuple[dict, dict]:
    calibration = {
        "mode": "uncalibrated_absolute_thresholds",
        "requested": True,
        "reference_ligand_file": reference_ligand_file,
        "native_baseline": None,
    }
    ref_path = Path(reference_ligand_file)
    if not ref_path.exists():
        steps_skipped.append({
            "step": "posecheck_native_calibration",
            "reason": f"reference ligand not found: {ref_path}",
        })
        return posecheck_thresholds, calibration

    try:
        native_pc = run_posecheck(ref_path, receptor_path, num_conformers=num_conformers)
        native_entry = next(
            (x for x in native_pc if x.get("strain_energy_kcal") is not None),
            None,
        )
        if not native_entry:
            steps_skipped.append({
                "step": "posecheck_native_calibration",
                "reason": "native PoseCheck returned no numeric baseline",
            })
            return posecheck_thresholds, calibration

        baseline_strain = float(native_entry.get("strain_energy_kcal") or 0.0)
        baseline_clashes = int(native_entry.get("steric_clash_count") or 0)
        native_baseline = {
            "strain_energy_kcal": baseline_strain,
            "steric_clash_count": baseline_clashes,
            "steric_clash_count_all_atom": native_entry.get("steric_clash_count_all_atom"),
            "steric_clash_count_mode": native_entry.get("steric_clash_count_mode", "heavy_atom_only"),
        }
        if (
            baseline_strain > baseline_limits["strain_max_kcal"]
            or baseline_clashes > baseline_limits["clash_max_count"]
        ):
            calibration = {
                "mode": "native_baseline_rejected_absolute_thresholds",
                "requested": True,
                "reference_ligand_file": str(ref_path),
                "native_baseline": native_baseline,
                "baseline_limits": dict(baseline_limits),
                "reason": "native PoseCheck baseline exceeds configured sanity limits",
            }
            steps_skipped.append({
                "step": "posecheck_native_calibration",
                "reason": calibration["reason"],
                "native_baseline": native_baseline,
                "baseline_limits": dict(baseline_limits),
            })
            return posecheck_thresholds, calibration

        calibrated = {
            "strain_warn_kcal": baseline_strain + posecheck_thresholds["strain_warn_kcal"],
            "strain_fail_kcal": baseline_strain + posecheck_thresholds["strain_fail_kcal"],
            "clash_warn_count": baseline_clashes + posecheck_thresholds["clash_warn_count"],
            "clash_fail_count": baseline_clashes + posecheck_thresholds["clash_fail_count"],
        }
        calibration = {
            "mode": "delta_from_reference",
            "requested": True,
            "reference_ligand_file": str(ref_path),
            "native_baseline": native_baseline,
            "allowed_excess": dict(posecheck_thresholds),
            "baseline_limits": dict(baseline_limits),
        }
        steps_completed.append("posecheck_native_calibration")
        return calibrated, calibration
    except Exception as exc:
        steps_skipped.append({"step": "posecheck_native_calibration", "reason": str(exc)})
        return posecheck_thresholds, calibration


def _summary_for_evaluation(
    evaluation: list[dict],
    steps_completed: list,
    steps_skipped: list,
) -> dict:
    verdicts = [e.get("verdict", "unknown") for e in evaluation]
    return {
        "total_poses": len(evaluation),
        "pass": verdicts.count("pass"),
        "warn": verdicts.count("warn"),
        "fail": verdicts.count("fail"),
        "not_evaluated": verdicts.count("not_evaluated"),
        "steps_completed": _dedupe(steps_completed),
        "steps_skipped": _dedupe_dicts(steps_skipped),
    }


def _write_final_outputs(
    output_dir: Path,
    evaluation: list[dict],
    evaluation_input_poses: list[dict],
    selected_poses: list[dict],
    posecheck_thresholds: dict,
    posecheck_calibration: dict,
    reference_ligand_file: str | None,
    docking_scores_file: str | None,
    docking_artifact_lineage: dict | None,
    receptor_identity: dict,
    steps_completed: list,
    steps_skipped: list,
    files_extra: dict,
    complete: bool,
    batch_reports: list[dict],
    selection_strategy: str,
    poses_per_compound: int,
    batch_size: int,
    max_poses: int | None,
    overall_timeout_seconds: int,
    partial: bool = False,
    top_n: int = 10,
    auxiliary_results: dict | None = None,
    operational_errors: list[dict] | None = None,
    screening_mode: object = DEFAULT_SCREENING_MODE,
    posecheck_num_conformers: int = 10,
) -> dict:
    completed = _dedupe(steps_completed)
    skipped = _dedupe_dicts(steps_skipped)
    active_mode = resolve_screening_mode(screening_mode)
    summary = _summary_for_evaluation(evaluation, completed, skipped)
    pose_exhausted = _pose_exhausted_compounds(evaluation, evaluation_input_poses)
    observed = {int(row.get("source_pose_index", row.get("pose_index", -1))) for row in evaluation}
    unresolved = {
        row.get("mol_name") for row in evaluation
        if row.get("verdict") == "not_evaluated" or any(
            state == "tool_failed" and tool != "prolif"
            for tool, state in (row.get("evidence_status") or {}).items()
        )
    }
    accepted = {row.get("mol_name") for row in evaluation if row.get("verdict") in {"pass", "warn"}} - unresolved
    unevaluated = [
        {
            "mol_name": pose.get("mol_name"),
            "pose_id": pose.get("pose_id"),
            "source_pose_index": pose["source_pose_index"],
            "compound_pose_rank": pose.get("compound_pose_rank"),
            "status": "not_evaluated",
            "reason": (
                "earlier_pose_accepted" if pose.get("mol_name") in accepted
                else "earlier_pose_unresolved" if pose.get("mol_name") in unresolved
                else "pending_evaluation"
            ),
        }
        for pose in evaluation_input_poses
        if int(pose["source_pose_index"]) not in observed
    ]
    summary.update({
        "total_evaluation_input_poses": len(evaluation_input_poses),
        "total_selected_poses": len(selected_poses),
        "total_unevaluated_poses": len(unevaluated),
        "complete": complete,
        "batches_completed": len(batch_reports),
        # A compound with no surviving pose has no selection row, so without
        # these counters it would leave the run unaccounted for.
        "pose_exhausted_compounds": len(pose_exhausted),
        "pose_exhausted_rerun_docking": 0,  # Legacy summary field.
        "pose_exhausted_reject_molecule": sum(
            1 for item in pose_exhausted if item["route"] == "reject_molecule"
        ),
    })
    report = {
        "success": bool(complete),
        "partial": bool(not complete),
        "summary": summary,
        "selection": {
            "pre_qc_strategy": selection_strategy,
            "poses_per_compound": poses_per_compound,
            "post_qc_strategy": "best_eligible_per_compound",
            "eligible_verdicts": sorted(_POSE_SURVIVOR_VERDICTS),
            "authoritative_artifact": files_extra.get(
                "selection_manifest_csv"
            ),
            "authoritative_selected_poses": [
                {
                    "selection_rank": pose.get("selection_rank"),
                    "mol_name": pose.get("mol_name"),
                    "source_pose_index": pose.get("source_pose_index"),
                    "pose_id": pose.get("pose_id"),
                    "post_qc_verdict": pose.get("post_qc_verdict"),
                }
                for pose in selected_poses
            ],
            "max_poses": max_poses,
            "batch_size": batch_size,
            "evaluation_input_poses": len(evaluation_input_poses),
            "selected_poses": len(selected_poses),
            # Keep the rejected molecules' reasons even though they have no
            # row in ``authoritative_selected_poses``.
            "pose_exhausted_compounds": pose_exhausted,
            # Supplementary docking poses remain in the input SDF/manifest,
            # outside the evaluated cohort and all downstream candidate lists.
            "unevaluated_poses": unevaluated,
        },
        "batches": batch_reports,
        "timings_seconds": {
            name: round(sum(batch.get("timings_seconds", {}).get(name, 0.0) for batch in batch_reports), 3)
            for name in {key for batch in batch_reports for key in batch.get("timings_seconds", {})}
        },
        "evaluation": evaluation,
        "source_docking_artifact_set": docking_artifact_lineage,
        "parameters": {
            "posecheck_thresholds": posecheck_thresholds,
            "posecheck_calibration": posecheck_calibration,
            "posecheck_num_conformers": posecheck_num_conformers,
            "reference_ligand_file": reference_ligand_file,
            "docking_scores_file": docking_scores_file,
            "overall_timeout_seconds": overall_timeout_seconds,
            # Level 1 labels each metric; level 2 folds the labels per parent
            # category under this mode and only then decides the pose verdict.
            "gate_contract": "two_level_label_then_category_fold",
            "screening_mode": dict(active_mode),
            **receptor_identity,
        },
        "steps_completed": completed,
        "steps_skipped": skipped,
        "top_candidates": _generate_top_candidates(evaluation, top_n),
        "auxiliary_results": auxiliary_results or {},
    }
    if operational_errors:
        report["operational_errors"] = _dedupe_dicts(operational_errors)
    cohort_manifest: dict[str, dict] = {}
    for row in evaluation or [receptor_identity]:
        identity = {
            field: row.get(field)
            for field in SCIENTIFIC_COHORT_OUTPUT_FIELDS
        }
        if identity.get("cohort_id") not in (None, ""):
            cohort_manifest.setdefault(str(identity["cohort_id"]), identity)
    report["scientific_cohorts"] = list(cohort_manifest.values())
    report["files"] = {
        "evaluation_report_json": str(output_dir / "evaluation_report.json"),
        "evaluation_summary_json": str(output_dir / "evaluation_summary.json"),
        **files_extra,
    }

    eval_file = output_dir / "evaluation_report.json"
    summary_file = output_dir / "evaluation_summary.json"
    copies = [output_dir / "evaluation_report.partial.json"] if partial and not complete else []
    write_json(eval_file, report, additional_paths=copies)
    write_json(summary_file, summary)

    return report


def _chunks(items: list[dict], size: int):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def _dedupe(items: list) -> list:
    seen = set()
    result = []
    for item in items:
        key = json.dumps(item, sort_keys=True, default=str) if not isinstance(item, str) else item
        if key in seen:
            continue
        seen.add(key)
        result.append(item)
    return result


def _dedupe_dicts(items: list[dict]) -> list[dict]:
    return _dedupe(items)


def _canonical_posebusters_checks(raw_checks) -> dict[str, object]:
    """Return JSON-safe, stable PoseBusters check IDs.

    PoseBusters exposes DataFrame columns as tuples under some versions.  Such
    keys cannot be serialized by ``json.dumps`` and aliases can collapse onto
    the same registered check.  When that happens, a boolean failure wins so a
    duplicate label can never hide a hard gate.
    """
    if not isinstance(raw_checks, dict):
        return {}

    checks: dict[str, object] = {}
    for position, (raw_key, value) in enumerate(raw_checks.items()):
        check_id = canonical_posebusters_check_id(raw_key) or f"check_{position}"
        if check_id not in checks:
            checks[check_id] = value
            continue

        existing = checks[check_id]
        if isinstance(existing, bool) and isinstance(value, bool):
            checks[check_id] = existing and value
        elif isinstance(value, bool):
            checks[check_id] = value
    return checks


def _canonical_posebusters_violations(raw_violations) -> list[str]:
    """Normalize and de-duplicate PoseBusters violation labels."""
    if raw_violations is None:
        return []
    items = raw_violations if isinstance(raw_violations, list) else [raw_violations]
    violations: list[str] = []
    seen: set[str] = set()
    for raw_key in items:
        check_id = canonical_posebusters_check_id(raw_key)
        if check_id and check_id not in seen:
            seen.add(check_id)
            violations.append(check_id)
    return violations


def _combine_results(
    pb_results: list[dict] | None,
    pc_results: list[dict] | None,
    plip_results: list[dict] | None,
    blueprint_results: list[dict] | None,
    pose_names: list[str] | None = None,
    posecheck_thresholds: dict | None = None,
    docking_scores: list[dict] | None = None,
    blueprint_expected: bool = False,
    pocket_expected: bool = False,
    posecheck_calibrated: bool = False,
    tool_switches: dict | None = None,
    prolif_results: list[dict] | None = None,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> list[dict]:
    """Merge results from all tools into per-molecule evaluation entries."""
    switches = tool_switches or dict(EVALUATOR_TOOL_DEFAULTS)
    # Determine number of poses from whichever result is available
    n = len(pose_names or [])
    for results in [pb_results, pc_results, plip_results, prolif_results, blueprint_results]:
        if results:
            n = max(n, len(results))

    pb_by_index = _results_by_pose_index(pb_results)
    pc_by_index = _results_by_pose_index(pc_results)
    plip_by_index = _results_by_pose_index(plip_results)
    prolif_by_index = _results_by_pose_index(prolif_results)
    blueprint_by_index = _results_by_pose_index(blueprint_results)
    scores_by_index = _results_by_pose_index(docking_scores)

    evaluation = []
    for i in range(n):
        entry = {
            "pose_index": i,
            "mol_name": pose_names[i] if pose_names and i < len(pose_names) else f"pose_{i}",
            "metric_context": {
                "enabled_tools": _scoring_plan_tools(
                    switches,
                    blueprint_expected=blueprint_expected,
                    pocket_expected=pocket_expected,
                ),
            },
        }

        # PoseBusters
        if i in pb_by_index:
            pb = pb_by_index[i]
            checks = _canonical_posebusters_checks(pb.get("checks"))
            check_errors = {
                canonical_posebusters_check_id(key) or f"check_{position}": reason
                for position, (key, reason) in enumerate(
                    (pb.get("check_errors") or {}).items()
                )
            }
            violations = _canonical_posebusters_violations(pb.get("violations"))
            boolean_checks = [
                value for value in checks.values() if isinstance(value, bool)
            ]
            success_ratio = pb.get("success_ratio")
            if success_ratio is None and boolean_checks:
                success_ratio = sum(boolean_checks) / len(boolean_checks)
            try:
                success_ratio = (
                    float(success_ratio)
                    if success_ratio is not None
                    else None
                )
            except (TypeError, ValueError):
                success_ratio = None
            pb_success = pb.get("success")
            if not isinstance(pb_success, bool):
                # Older PoseBusters adapters may expose only an aggregate
                # ``all_passed``/violation summary.  That is still real tool
                # evidence even when per-check columns are unavailable.  A row
                # with neither check values nor an aggregate remains a tool
                # failure, so absence of evidence can never become a pass.
                pb_success = bool(boolean_checks) or any((
                    isinstance(pb.get("all_passed"), bool),
                    isinstance(pb.get("num_violations"), (int, float))
                    and not isinstance(pb.get("num_violations"), bool),
                    bool(pb.get("violations")),
                ))
            all_passed = pb.get("all_passed")
            if not isinstance(all_passed, bool):
                # Derive the aggregate from the checks that did resolve instead of
                # defaulting it to False.  A False here would fabricate an
                # aggregate failure from a missing summary key rather than from
                # any observed failing check.
                all_passed = bool(boolean_checks) and all(boolean_checks) and not violations
            if check_errors and not violations and all(boolean_checks):
                # An unavailable check is neither a pass nor an unattributed
                # physical failure. Keep partial measured evidence usable.
                all_passed = None
            # ``num_violations`` stays absent when the producer never reported it,
            # so a consumer can tell "no violations" from "never counted".
            num_violations = pb.get("num_violations")
            if num_violations is None and boolean_checks:
                num_violations = len(violations)
            entry["posebusters"] = {
                "success": pb_success,
                "num_checks": len(boolean_checks),
                "all_passed": None if all_passed is None else pb_success and all_passed,
                "violations": violations,
                "num_violations": num_violations,
                "checks": checks,
                "check_errors": check_errors,
                "success_ratio": success_ratio,
                "success_ratio_display": (
                    round(success_ratio, 6)
                    if success_ratio is not None else None
                ),
            }

        # PoseCheck
        if i in pc_by_index:
            pc = pc_by_index[i]
            pc_success = pc.get("success")
            if not isinstance(pc_success, bool):
                # Preserve compatibility with older successful records while
                # refusing to treat an error-only row as usable evidence.
                pc_success = not bool(pc.get("error"))
            entry["posecheck"] = {
                "success": pc_success,
                "error": pc.get("error"),
                "num_conformers": pc.get("num_conformers"),
                "steric_clash_count": pc.get("steric_clash_count"),
                "steric_clash_count_all_atom": pc.get("steric_clash_count_all_atom"),
                "steric_clash_count_mode": pc.get("steric_clash_count_mode", "heavy_atom_only"),
                "strain_energy_kcal": pc.get("strain_energy_kcal"),
                "strain_energy_kcal_display": pc.get(
                    "strain_energy_kcal_display"
                ),
            }

        # PLIP interactions
        if i in plip_by_index:
            plip = plip_by_index[i]
            # ``num_interactions`` is left absent rather than defaulted to 0: the
            # registered rule has hard_lower=2, so materializing a 0 for a worker
            # row that never reported a count hard-rejects the pose on the
            # strength of a missing key.  ProLIF below already reads it straight
            # through, and the extractor maps None to not_applicable.
            entry["plip"] = {
                "success": plip.get("success", True),
                "num_interactions": plip.get("num_interactions"),
                "interactions": plip.get("interactions", []),
                "error": plip.get("error"),
                "status": plip.get("status"),
                "reason": plip.get("reason"),
            }

        # ProLIF interactions
        if i in prolif_by_index:
            prolif = prolif_by_index[i]
            prolif_error = prolif.get("error")
            entry["prolif"] = {
                "success": bool(prolif.get("success", False)) and not prolif_error,
                "num_interactions": prolif.get("num_interactions"),
                "interactions": prolif.get("interactions", []),
                "evaluator_capabilities": prolif.get("evaluator_capabilities"),
                "diagnostics": prolif.get("diagnostics"),
                "residue_level_bits": prolif.get("residue_level_bits"),
                "tanimoto_to_reference": prolif.get("tanimoto_to_reference"),
                "reference_comparison": prolif.get("reference_comparison"),
                "error": prolif_error,
                "status": prolif.get("status"),
                "reason": prolif.get("reason"),
            }

        if ((entry.get("prolif") or {}).get("reference_comparison") or {}).get("status") == "not_applicable":
            entry["metric_context"].setdefault("not_applicable_metric_ids", []).append("prolif.tanimoto_to_reference")
        for profiler in ("plip", "prolif"):
            block = entry.get(profiler) or {}
            if block.get("success") is True and not block.get("error"):
                from collections import Counter
                block["interaction_type_counts"] = dict(Counter(
                    item.get("interaction_subtype") or item.get("interaction_type") or item.get("type")
                    for item in block.get("interactions") or [] if isinstance(item, dict)))
                if profiler == "plip":
                    block["interaction_fingerprint"] = block.get("interactions") or []

        # Blueprint match
        if i in blueprint_by_index:
            entry["blueprint_match"] = blueprint_by_index[i]

        if i in scores_by_index:
            score = scores_by_index[i]
            entry["docking_scores"] = {
                key: value
                for key, value in score.items()
                if key not in {"pose_index", "source_pose_index", "mol_name", "name"}
            }
            entry["pose_scores"] = collect_pose_scores(score)
            entry["primary_score_metric_id"] = resolve_primary_metric(
                [entry["pose_scores"]],
                preferred_metric_ids=POSE_SELECTION_PREFERENCE,
            )

        # Evidence status distinguishes "available" from "disabled_by_config"
        # (user opt-out, no penalty) and "missing" (expected but absent, soft-warn).
        entry["evidence_status"] = {
            "posebusters": _evidence_state(entry, "posebusters", switches["posebusters"]),
            "posecheck": _evidence_state(entry, "posecheck", switches["posecheck"]),
            "posecheck_calibration": (
                "delta_from_reference"
                if posecheck_calibrated
                else "absolute_thresholds"
            ),
            "plip": _evidence_state(entry, "plip", switches["plip"]),
            "prolif": _evidence_state(entry, "prolif", switches["prolif"]),
            "blueprint": _blueprint_evidence_state(
                entry,
                enabled=switches["blueprint_match"],
                expected=blueprint_expected,
            ),
        }

        # Final verdict
        entry["verdict"] = _compute_verdict(
            entry, posecheck_thresholds, screening_mode=screening_mode
        )
        evaluation.append(entry)

    return evaluation


def _evidence_state(entry: dict, key: str, enabled: bool) -> str:
    """Evidence state for a tool that produces its own entry block.

    available          — the tool ran and attached data
    not_evaluable      — tool ran, but the scientific question was not evaluable
    not_applicable     — tool ran, but the scientific question did not apply
    disabled_by_config — user turned it off; no penalty, no pass credit
    missing            — expected on but produced nothing (soft-warn)
    """
    result = entry.get(key) or {}
    if key in entry and result.get("success") is False:
        return "tool_failed"
    if key in entry:
        status = str(result.get("status") or "").strip().lower()
        if status in {"not_evaluable", "not_applicable"}:
            return status
        return "available"
    if not enabled:
        return "disabled_by_config"
    return "missing"


_AUXILIARY_EVIDENCE_METRIC_PREFIXES = {
    "pocket_occupancy": ("pocket_occupancy.",),
    "rtmscore": ("rtmscore.",),
    "mmgbsa": ("amber_mmpbsa.", "mmgbsa."),
}


def _auxiliary_evidence_present(entry: dict, name: str) -> bool:
    """True when an auxiliary evaluator attached a usable value to this pose.

    These three tools do not own an ``evidence_status`` key because they attach
    to the pose after the profiler merge -- pocket occupancy as
    ``pocket_contact_fraction``, RTMScore and MM-GBSA as rows in ``pose_scores``.
    """
    if name in {"rdkit_pose", "freesasa", "mdanalysis", "arpeggio"}:
        block = entry.get(name) or {}
        return block.get("success") is True and bool(block.get("metrics"))
    if name == "pocket_occupancy":
        occupancy = entry.get("pocket_contact_fraction")
        if isinstance(occupancy, dict):
            return any(value is not None for value in occupancy.values())
        return isinstance(occupancy, (int, float)) and not isinstance(occupancy, bool)
    prefixes = _AUXILIARY_EVIDENCE_METRIC_PREFIXES.get(name) or ()
    if not prefixes:
        return False
    for score in entry.get("pose_scores") or []:
        if not isinstance(score, dict):
            continue
        metric_id = str(score.get("metric_id") or "")
        if metric_id.startswith(prefixes) and score.get("value") is not None:
            return True
    return False


def _blueprint_evidence_state(entry: dict, *, enabled: bool, expected: bool) -> str:
    """Return a tri-state-aware status for the Blueprint comparison."""
    if "blueprint_match" in entry:
        verdict = (entry.get("blueprint_match") or {}).get("verdict")
        if verdict == "unresolved":
            return "unresolved"
        if verdict == "not_applicable":
            return "not_applicable"
        return "available"
    if not enabled:
        return "disabled_by_config"
    return "missing" if expected else "not_requested"


def _read_pose_names(docked_sdf: Path) -> list[str]:
    """Read molecule names from the docked SDF in pose order."""
    try:
        from rdkit import Chem

        names = []
        supplier = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)
        for i, mol in enumerate(supplier):
            if mol is not None and mol.HasProp("_Name") and mol.GetProp("_Name").strip():
                names.append(mol.GetProp("_Name").strip())
            else:
                names.append(f"pose_{i}")
        return names
    except ImportError:
        names = []
        for block in docked_sdf.read_text(errors="ignore").split("$$$$"):
            lines = [line.strip() for line in block.splitlines() if line.strip()]
            if not lines:
                continue
            names.append(lines[0] or f"pose_{len(names)}")
        return names


def _load_posecheck_thresholds(input_data: dict) -> dict:
    """Load calibrated PoseCheck thresholds, allowing per-run overrides."""
    thresholds = dict(DEFAULT_POSECHECK_THRESHOLDS)
    overrides = input_data.get("posecheck_thresholds") or {}
    for key in thresholds:
        if key in overrides:
            thresholds[key] = overrides[key]
        elif key in input_data:
            thresholds[key] = input_data[key]

    thresholds["strain_warn_kcal"] = float(thresholds["strain_warn_kcal"])
    thresholds["strain_fail_kcal"] = float(thresholds["strain_fail_kcal"])
    thresholds["clash_warn_count"] = int(thresholds["clash_warn_count"])
    thresholds["clash_fail_count"] = int(thresholds["clash_fail_count"])
    if thresholds["strain_warn_kcal"] > thresholds["strain_fail_kcal"]:
        raise ValueError(
            "PoseCheck strain thresholds require strain_warn_kcal <= "
            "strain_fail_kcal"
        )
    if thresholds["clash_warn_count"] > thresholds["clash_fail_count"]:
        raise ValueError(
            "PoseCheck clash thresholds require clash_warn_count <= "
            "clash_fail_count"
        )
    return thresholds


def _load_posecheck_baseline_limits(input_data: dict) -> dict:
    """Load native-ligand sanity limits for PoseCheck calibration."""
    limits = dict(DEFAULT_POSECHECK_NATIVE_BASELINE_LIMITS)
    overrides = input_data.get("posecheck_native_baseline_limits") or {}
    for key in limits:
        if key in overrides:
            limits[key] = overrides[key]
        elif key in input_data:
            limits[key] = input_data[key]

    limits["strain_max_kcal"] = float(limits["strain_max_kcal"])
    limits["clash_max_count"] = int(limits["clash_max_count"])
    return limits


def _registered_metric_assessments(
    entry: dict,
    posecheck_thresholds: dict,
    *,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> list[dict]:
    """Label emitted pose scalars, then fold the labels per parent category.

    Level 1 marks each metric ``pass``/``flag``/``fail``; level 2 requires
    a strict pass majority in every evaluated parent category for ``keep``
    unless a hard rejection applies.  PoseCheck rules use the same default
    thresholds as the shared table, but a valid native-ligand calibration may
    raise them for this run.  Reflect those effective bounds in both the label
    and the emitted assessment instead of applying the uncalibrated table a
    second time.
    """
    mode = resolve_screening_mode(screening_mode)
    # Stage 07 owns pose/docking evidence only.  Passing the combined
    # pose+property registry here used to append all Stage 08 property and
    # ADMET metrics to every pose as ``not_planned``.  Those rows did not gate
    # the verdict, but they obscured the stage boundary and greatly inflated
    # large evaluation reports.
    assessments = assess_record(entry, POSE_EVALUATION_METRIC_RULES)
    posecheck_overrides = {
        "posecheck.strain_energy_kcal": {
            "soft_upper": posecheck_thresholds["strain_warn_kcal"],
            "hard_upper": posecheck_thresholds["strain_fail_kcal"],
        },
        "posecheck.steric_clash_count": {
            "soft_upper": posecheck_thresholds["clash_warn_count"],
            "hard_upper": posecheck_thresholds["clash_fail_count"],
        },
    }
    for assessment in assessments:
        overrides = posecheck_overrides.get(assessment["metric_id"])
        if not overrides:
            continue
        effective_rule = dict(_POSE_METRIC_RULE_BY_ID[assessment["metric_id"]])
        effective_rule.update(overrides)
        label = metric_label(assessment.get("value"), effective_rule)
        molecule_label = label
        if (
            assessment.get("molecule_effect") in NON_MOLECULAR_EFFECTS
            and label == "fail"
        ):
            molecule_label = "flag"
        # Both spellings are rewritten together: the category fold reads
        # ``molecule_label`` first, so leaving the legacy keys alone here would
        # silently drop the calibration from the verdict.
        assessment.update({
            **overrides,
            "label": label,
            "molecule_label": molecule_label,
            "gate": _LABEL_TO_LEGACY_GATE[label],
            "molecule_gate": _LABEL_TO_LEGACY_GATE[molecule_label],
        })

    ranking_blocks = [
        assessment["metric_id"]
        for assessment in assessments
        if assessment.get("molecule_effect") in RANKING_BLOCKING_EFFECTS
        and assessment.get("label") == "fail"
    ]
    fold = molecule_gate_report(assessments, screening_mode=mode)
    entry["metric_assessments"] = assessments
    entry["metric_gate"] = str(fold["verdict"])
    entry["screening_mode"] = mode["mode_id"]
    entry["screening_mode_name"] = mode["name"]
    entry["screening_mode_eliminates"] = bool(mode["eliminates"])
    entry["screening_tripped_categories"] = list(fold["tripped_categories"])
    entry["metric_label_counts"] = {
        "fail": int(fold["count_fail"]),
        "flag": int(fold["count_flag"]),
        "labeled": int(fold["count_labeled"]),
    }
    entry["metric_category_folds"] = fold["categories"]
    entry["ranking_blocked"] = bool(ranking_blocks)
    entry["ranking_blocked_by"] = ranking_blocks
    return assessments


def _verdict_from_reasons(
    entry: dict,
    reasons_fail: list[str],
    reasons_warn: list[str],
    reasons_na: list[str] | None = None,
    *,
    metric_gate: str | None = None,
    hard_failures: list[str] | None = None,
) -> str:
    """Publish the audit trail, then derive the verdict from the two-level fold.

    ``reasons_fail`` holds level-1 ``fail`` marks.  Under the executable policy a single
    mark never decides the pose on its own, so it is recorded for audit and the
    verdict comes from ``metric_gate`` -- the per-category fold under the active
    screening mode.  ``hard_failures`` is the narrow set of non-metric defects
    (unreadable identity, a PoseBusters failure that cannot be attributed to any
    registered check) that are structural rather than threshold-based and so
    still fail outright.
    """
    fail = list(dict.fromkeys(reasons_fail))
    warn = list(dict.fromkeys(reasons_warn))
    hard = list(dict.fromkeys(hard_failures or []))
    # ``not_applicable`` records evidence that could not be judged at all.  It is
    # published for audit but drives no verdict, because the executable policy removes a
    # missing value from numerator and denominator rather than scoring it.
    entry["verdict_reasons"] = {
        "fail": fail,
        "warn": warn,
        "hard": hard,
        "not_applicable": list(dict.fromkeys(reasons_na or [])),
    }
    if hard:
        return "fail"
    if metric_gate == "reject":
        return "fail"
    if metric_gate == "modify":
        return "warn"
    if metric_gate == "not_evaluated":
        return "not_evaluated"
    return "pass"


def _attach_diagnosis(
    entry: dict,
    *,
    subpockets: list[dict] | None,
    posecheck_thresholds: dict,
) -> None:
    """Attach human-readable guidance without weakening a two-level verdict."""
    diagnosis, recommendation, modification_hint = generate_diagnosis(
        entry,
        subpockets,
        posecheck_thresholds,
    )
    reasons = entry.get("verdict_reasons") or {}
    hard_reasons = reasons.get("hard") or []
    fail_reasons = reasons.get("fail") or []
    if hard_reasons or entry.get("metric_gate") == "reject":
        # A structural defect or a category that tripped the screening mode
        # cannot be talked back into a modify by the legacy diagnosis heuristic.
        recommendation = "reject"
    if hard_reasons:
        diagnosis.append(
            "HARD_FAIL: " + ", ".join(str(item) for item in hard_reasons)
        )
    if fail_reasons:
        diagnosis.append(
            "METRIC_FAIL_MARKS: " + ", ".join(str(item) for item in fail_reasons)
        )
    tripped = entry.get("screening_tripped_categories") or []
    if tripped:
        diagnosis.append(
            "SCREENING_TRIPPED[%s]: %s"
            % (entry.get("screening_mode") or "", ", ".join(str(x) for x in tripped))
        )
    entry["diagnosis"] = diagnosis
    entry["recommendation"] = recommendation
    entry["modification_hint"] = modification_hint


def _attribute_legacy_posebusters_violations(entry: dict) -> bool:
    """Expose a legacy ``violations`` list to the registered per-check rules.

    Older PoseBusters adapters report only an aggregate (``all_passed`` /
    ``num_violations``) plus a ``violations`` list naming the checks that failed,
    with no per-check ``checks`` map.  ``posebusters_check`` extractors read
    ``checks`` only, so those named failures resolved to nothing and the aggregate
    could not be evaluated by the registered two-level fold in the executable policy.
    A *named* check is attributable, so it is republished
    here as a value the registered rules can resolve and the fold can count.

    The failure is recorded as the number ``0.0`` rather than ``False`` on purpose.
    ``posebusters_check`` accepts either, so every named check still earns its own
    level-1 ``fail`` in its own parent category; ``posebusters_success_ratio``
    averages only ``bool`` entries, so it keeps returning ``None`` and stays NA.
    That is the honest outcome: the list says *which* checks failed but not how
    many ran, and the executable policy keeps a value with an unknown denominator out of the
    fold instead of scoring it 0.

    Returns True when a derived map was installed.
    """
    pb = entry.get("posebusters")
    if not isinstance(pb, dict) or pb.get("success") is False:
        return False
    existing = pb.get("checks")
    if isinstance(existing, dict) and existing:
        return False
    named = {
        canonical_posebusters_check_id(name)
        for name in (pb.get("violations") or [])
    }
    named.discard("")
    derived = {
        check_id: 0.0
        for check_id in sorted(named)
        if check_id in _POSEBUSTERS_CHECK_RULE_BY_CHECK_ID
    }
    if not derived:
        return False
    pb["checks"] = derived
    # Marked so a reader never mistakes an attributed subset for the tool's own
    # per-check output; ``num_checks`` deliberately stays at the number of boolean
    # columns PoseBusters actually delivered.
    pb["checks_source"] = "derived_from_violations"
    return True


def _compute_verdict(
    entry: dict,
    posecheck_thresholds: dict | None = None,
    *,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> str:
    """Compute the overall verdict from evidence and the per-category fold."""
    thresholds = posecheck_thresholds or DEFAULT_POSECHECK_THRESHOLDS
    # Must run before the assessments are built: it is what lets a legacy
    # violation reach its registered rule instead of the hard-failure branch.
    _attribute_legacy_posebusters_violations(entry)
    reasons_fail: list[str] = []
    reasons_warn: list[str] = []
    reasons_na: list[str] = []
    hard_failures: list[str] = []
    assessments = _registered_metric_assessments(
        entry, thresholds, screening_mode=screening_mode
    )
    assessment_by_id = {
        assessment["metric_id"]: assessment for assessment in assessments
    }

    evidence = entry.get("evidence_status", {})
    if evidence.get("identity_qc") == "failed":
        # Not a threshold: the pose could not be identified at all.
        return _verdict_from_reasons(
            entry, [], [], hard_failures=["identity_qc"]
        )
    available = sum(
        evidence.get(name) == "available"
        for name in ("posebusters", "posecheck", "plip", "blueprint")
    )
    # Pose QC evidence is not limited to the core profilers. Pocket occupancy,
    # RTMScore and MM-GBSA each resolve registered pose metrics, so a run that
    # enabled one of them and nothing else does hold evidence about the pose.
    # Counting only the core profilers made every pose in such a run "not_evaluated" while
    # the survivor and top-candidate filters still dropped that verdict, so the
    # run silently produced no poses at all. The executable policy asks for Score=NA plus a
    # coverage/confidence note, never a deletion.
    available += sum(
        _auxiliary_evidence_present(entry, name)
        for name in ("pocket_occupancy", "rtmscore", "mmgbsa", "rdkit_pose", "freesasa", "mdanalysis", "arpeggio")
    )
    # Invariant: a pose with zero real evidence can never pass. Turning tools off
    # does not manufacture evidence — it just removes their penalties below.
    if available == 0:
        warnings = ["no_usable_evaluation_evidence"]
        for tool in ("posebusters", "posecheck"):
            if evidence.get(tool) in {"missing", "tool_failed"}:
                warnings.append(f"{tool}_{evidence.get(tool)}")
        entry["verdict_reasons"] = {
            "fail": [],
            "warn": warnings,
            "hard": [],
            "not_applicable": [],
        }
        return "not_evaluated"
    # A tool the user explicitly disabled must not add a warn penalty ("config
    # says don't use it, so its absence is not a defect"). Only a tool that was
    # expected (enabled) but produced nothing counts as a missing-evidence warn.
    if evidence.get("posebusters") in {"missing", "tool_failed"}:
        reasons_warn.append(f"posebusters_{evidence.get('posebusters')}")
    for check_id in sorted((entry.get("posebusters") or {}).get("check_errors") or {}):
        reasons_warn.append(f"posebusters_check_unavailable:{check_id}")
    if evidence.get("posecheck") in {"missing", "tool_failed"}:
        reasons_warn.append(f"posecheck_{evidence.get('posecheck')}")
    elif (
        evidence.get("posecheck") == "available"
        and evidence.get("posecheck_calibration") == "missing"
    ):
        reasons_warn.append("posecheck_uncalibrated")
    if evidence.get("blueprint") == "missing":
        reasons_warn.append("blueprint_missing")
    elif evidence.get("blueprint") == "unresolved":
        reasons_warn.append("blueprint_unresolved")
    if evidence.get("plip") in {"missing", "tool_failed"}:
        reasons_warn.append(f"plip_{evidence.get('plip')}")

    # PoseBusters checks with registered rules produce their own level-1 mark. A
    # check with no registered rule carries no threshold of its own, so
    # it is recorded NA rather than failed: PoseBusters reports several checks
    # (molecular_formula, mol_pred_loaded, "rmsd <= 2A") as False whenever no
    # reference ligand exists. A reference-dependent check without a reference
    # is NA and cannot be counted as fail. Legacy aggregate-only reports still remain strict because
    # they cannot distinguish a modify-only failure from an impossible one.
    pb = entry.get("posebusters", {})
    if pb and pb.get("success") is not False:
        checks = pb.get("checks")
        checks = checks if isinstance(checks, dict) else {}
        # A named violation is attributable evidence even when no per-check map
        # accompanies it, so it is tracked separately from ``checks``.
        named_violations = {
            canonical_posebusters_check_id(name)
            for name in (pb.get("violations") or [])
        }
        named_violations.discard("")
        aggregate_failed = (
            pb.get("all_passed") is False
            or _number(pb.get("num_violations"), 0) > POSEBUSTERS_VIOLATION_THRESHOLD
        )

        for assessment in assessments:
            if not assessment["metric_id"].startswith("posebusters."):
                continue
            if assessment.get("value") is None:
                continue
            if assessment.get("molecule_label") == "fail":
                reasons_fail.append(f"metric:{assessment['metric_id']}")
            elif assessment.get("molecule_label") == "flag":
                reasons_warn.append(f"metric:{assessment['metric_id']}")

        if checks:
            failed_check_ids = {
                canonical_posebusters_check_id(name)
                for name, value in checks.items()
                if isinstance(value, bool) and not value
            }
            failed_check_ids.update(
                canonical_posebusters_check_id(name)
                for name in (pb.get("violations") or [])
            )
            failed_check_ids.discard("")
            for check_id in sorted(failed_check_ids):
                rule = _POSEBUSTERS_CHECK_RULE_BY_CHECK_ID.get(check_id)
                if rule is None:
                    reasons_na.append(f"posebusters_unregistered:{check_id}")
                    continue
                assessment = assessment_by_id.get(rule["metric_id"])
                if assessment is None or assessment.get("value") is None:
                    reasons_na.append(f"posebusters_unresolved:{check_id}")
            if aggregate_failed and not failed_check_ids:
                # PoseBusters says the pose failed but names no registered
                # check. The evaluation contract defines no separate aggregate
                # hard gate, so preserve the
                # observation for audit without inventing a rejection rule.
                reasons_na.append("posebusters_unattributed_aggregate_failure")
        elif named_violations:
            # No per-check map, but the violation list names the failed checks, so
            # they were already attributed above (see
            # _attribute_legacy_posebusters_violations) or are reference-dependent
            # checks that must remain NA. Either way the
            # failure is accounted for by name; escalating the aggregate on top of
            # that would reject the pose on a single check.
            for check_id in sorted(named_violations):
                if check_id not in _POSEBUSTERS_CHECK_RULE_BY_CHECK_ID:
                    reasons_na.append(f"posebusters_unregistered:{check_id}")
        elif aggregate_failed:
            reasons_na.append("posebusters_unattributed_aggregate_failure")

    # PoseCheck uses the effective (possibly reference-calibrated) bounds placed
    # into its registered assessments above.
    for metric_id, fail_reason, warn_reason in (
        ("posecheck.strain_energy_kcal", "high_strain", "moderate_strain"),
        ("posecheck.steric_clash_count", "steric_clashes", "minor_steric_clashes"),
    ):
        assessment = assessment_by_id.get(metric_id)
        if assessment and assessment.get("molecule_label") == "fail":
            reasons_fail.append(fail_reason)
        elif assessment and assessment.get("molecule_label") == "flag":
            reasons_warn.append(warn_reason)

    # Blueprint verdict is a registered INTX metric, so a mismatch is a level-1
    # ``fail`` mark that only rejects once the INTX fold trips the active mode.
    bp = entry.get("blueprint_match", {})
    if bp:
        bp_verdict = bp.get("verdict", "")
        if bp_verdict == "fail":
            reasons_fail.append("blueprint_mismatch")
        elif bp_verdict == "warn":
            reasons_warn.append("blueprint_partial")
        elif bp_verdict == "unresolved":
            reasons_warn.append("blueprint_unresolved")
        if bp.get("profiler_agreement") == "disagreed":
            reasons_warn.append("blueprint_profiler_disagreement")

    # Remaining registered scalar metrics include PLIP interaction count,
    # required/desired Blueprint diagnostics, and available pocket occupancies.
    # Protocol-only and pose-export hard breaches are demoted to
    # molecule_label=flag and recorded in ranking_blocked_by.
    handled = {
        "blueprint_match.verdict",
        "posecheck.strain_energy_kcal",
        "posecheck.steric_clash_count",
    }
    for assessment in assessments:
        metric_id = assessment["metric_id"]
        if metric_id in handled or metric_id.startswith("posebusters."):
            continue
        if assessment.get("molecule_label") == "fail":
            reasons_fail.append(f"metric:{metric_id}")
        elif assessment.get("molecule_label") == "flag":
            reasons_warn.append(f"metric:{metric_id}")

    return _verdict_from_reasons(
        entry,
        reasons_fail,
        reasons_warn,
        reasons_na,
        metric_gate=entry.get("metric_gate"),
        hard_failures=hard_failures,
    )


def main():
    from io_utils import parse_skill_args, print_skill_result
    input_data, output_dir = parse_skill_args("Structure Evaluator Skill")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

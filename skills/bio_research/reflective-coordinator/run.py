#!/usr/bin/env python3
"""Reflective Coordinator Skill - Deterministic per-round assessment for the SBDD pipeline.

Usage:
    python run.py --input '{"evaluation_report_file": "eval.json", "property_filter_report_file": "prop.json", "round_number": 1}' --output-dir ./output/
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from decision_engine import (
    classify_molecule,
    evaluation_data_review,
    identify_bottleneck,
    is_scaffold_level,
    is_protocol_level,
    aggregate_modification_hints,
    property_data_review,
)
from ranking_engine import (
    prepare_pose_records,
    rank_candidates,
    screening_metrics,
    select_compound_poses,
)
from metric_scoring import (
    DOCKING_PROTOCOL_METRIC_RULES,
    DEFAULT_SCREENING_MODE,
    METRIC_RULES_VERSION,
    NORMALIZATION_METHOD,
    POSE_EVALUATION_METRIC_RULES,
    POSE_PROPERTY_METRIC_RULES,
    PROPERTY_FILTER_METRIC_RULES,
    SCORING_SCHEMA_VERSION,
    SCREENING_MODES,
    LABEL_TO_LEGACY_GATE as _LABEL_TO_LEGACY_GATE,
    metric_label,
    public_metric_rules,
    resolve_screening_mode,
    record_has_structural_reject,
    score_records,
    validate_metric_rules,
)
from pose_identity import LIGAND_LINEAGE_PROPERTIES, benchmark_label
from docking_artifact_contract import (
    evaluation_docking_lineage_error,
    load_committed_docking_artifacts,
)
from pipeline_registry import DOCKING_PRIMARY_ENGINES
from pose_scores import LEGACY_SCORE_SPECS


# ``not_evaluated`` sits between ``keep`` and ``modify``.  It is not a finding
# about the chemistry, so it must never dominate a real ``modify``/``reject``
#, but it must also stop a clean ``keep``: when one side of the
# combine says nothing was evaluated, the molecule has not earned nomination.
# This intentionally differs from ``metric_scoring._combine_gate_statuses``,
# where ``not_evaluated`` ranks lowest because there it merges partial views of
# one aggregate rather than two independent decisions about the same molecule.
_CLASS_SEVERITY = {"keep": 0, "not_evaluated": 1, "modify": 2, "reject": 3}
SCIENTIFIC_COHORT_COMPONENT_FIELDS = (
    "protein_state_id",
    "pocket_id",
    "round_id",
    "protocol_profile_id",
    "metric_catalog_version",
)
SCIENTIFIC_COHORT_FIELDS = (
    "cohort_id",
    *SCIENTIFIC_COHORT_COMPONENT_FIELDS,
    "cohort_identity_complete",
    "cohort_identity_missing_fields",
    "cohort_identity_derived_fields",
    "cohort_provenance_sources",
    "cohort_manifest_valid",
    "cohort_manifest_validation_errors",
)


def run(input_data: dict, output_dir: Path) -> dict:
    """Assess one round; campaign planning belongs to the Agent."""
    eval_file = input_data.get("evaluation_report_file")
    selection_manifest_file = input_data.get("selection_manifest_file")
    prop_file = input_data.get("property_filter_report_file")
    round_number = input_data.get("round_number", 1)
    max_rounds = input_data.get("max_rounds", 3)
    pass_rate_threshold = input_data.get("pass_rate_threshold", 0.3)
    min_candidates = input_data.get("min_candidates", 3)
    docking_scores_file = input_data.get("docking_scores_file")
    docking_report_file = input_data.get("docking_report_file")
    pocket_analysis_report_file = input_data.get("pocket_analysis_report_file")
    require_quality_gates = bool(input_data.get("require_quality_gates", True))
    min_validation_actives = int(input_data.get("min_validation_actives", 20))
    min_validation_decoys = int(input_data.get("min_validation_decoys", 50))
    min_validation_auroc = float(input_data.get("min_validation_auroc", 0.70))
    min_validation_auroc_hard = float(
        input_data.get("min_validation_auroc_hard", 0.50)
    )
    min_validation_ef1 = float(input_data.get("min_validation_ef1", 3.0))
    min_validation_ef1_hard = float(
        input_data.get("min_validation_ef1_hard", 1.0)
    )
    normalization_method = input_data.get(
        "normalization_method", NORMALIZATION_METHOD
    )
    metric_weight_overrides = input_data.get("metric_weight_overrides") or {}
    category_weight_overrides = input_data.get("category_weight_overrides") or {}
    # Level 2 of the gate contract.  ``candidate_source`` only
    # supplies the default when no mode is named explicitly: generated molecules
    # default to the strict V3, a sampled library to V2, and a user-supplied full
    # set to the marking-only V1.
    screening_mode = resolve_screening_mode(
        input_data.get("screening_mode"),
        context=input_data.get("candidate_source"),
    )
    round_history_file = input_data.get("round_history_file")

    # Interface-v2 tool switches. Absent `tools` ⇒ both on (legacy). The
    # require_quality_gates tool mirrors the legacy flat flag; screening_benchmark
    # gates the retrospective AUROC/EF computation.
    from io_utils import tool_enabled
    require_quality_gates = bool(
        tool_enabled(input_data, "require_quality_gates", default=require_quality_gates)
    )
    screening_benchmark_enabled = tool_enabled(input_data, "screening_benchmark", default=False)

    if not eval_file:
        return {"success": False, "error": "evaluation_report_file is required"}
    if normalization_method != NORMALIZATION_METHOD:
        return {
            "success": False,
            "error": (
                "normalization_method must be "
                f"{NORMALIZATION_METHOD!r}; got {normalization_method!r}"
            ),
        }
    if not isinstance(metric_weight_overrides, dict):
        return {
            "success": False,
            "error": "metric_weight_overrides must be an object mapping metric IDs to 0-5",
        }
    if not isinstance(category_weight_overrides, dict):
        return {
            "success": False,
            "error": "category_weight_overrides must be an object mapping category codes to 0-5",
        }
    for label, overrides in (
        ("metric_weight_overrides", metric_weight_overrides),
        ("category_weight_overrides", category_weight_overrides),
    ):
        for key, value in overrides.items():
            number = _number_or_none(value)
            if (
                not str(key).strip()
                or isinstance(value, bool)
                or number is None
                or not 0 <= number <= 5
            ):
                return {
                    "success": False,
                    "error": f"{label} entries must map non-empty IDs to numbers from 0 to 5",
                }
    if min_validation_auroc_hard > min_validation_auroc:
        return {
            "success": False,
            "error": "min_validation_auroc_hard must be <= min_validation_auroc",
        }
    if min_validation_ef1_hard > min_validation_ef1:
        return {
            "success": False,
            "error": "min_validation_ef1_hard must be <= min_validation_ef1",
        }

    eval_path = Path(eval_file).expanduser().resolve()
    eval_file = str(eval_path)
    if not eval_path.exists():
        return {"success": False, "error": f"evaluation report not found: {eval_path}"}

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load prior round history (optional)
    prior_history: list[dict] = []
    if round_history_file:
        hist_path = Path(round_history_file).expanduser().resolve()
        round_history_file = str(hist_path)
        if hist_path.exists():
            prior_history = json.loads(hist_path.read_text()).get("rounds", [])

    # Load evaluation results
    with open(eval_path) as f:
        eval_report = json.load(f)
    pose_exhaustion = {
        (row.get("cohort_id"), row.get("mol_name")): row
        for row in (eval_report.get("selection") or {}).get("pose_exhausted_compounds", [])
        if row.get("route") == "reject_molecule"
        and row.get("reason_code") == "all_poses_failed"
    }
    eval_results = eval_report.get("evaluation", [])
    # Attach using the evaluator's cohort identity before coordinator cohort
    # normalization; otherwise legacy/unspecified cohort IDs lose diagnostics.
    eval_results = [
        {
            **entry,
            "pose_exhaustion": deepcopy(pose_exhaustion.get(
                (entry.get("cohort_id"), entry.get("mol_name"))
            )),
        }
        for entry in eval_results
    ]
    try:
        (
            authoritative_pose_ids,
            authoritative_pose_records,
            pose_selection_metadata,
        ) = _load_authoritative_pose_selection(
            eval_report,
            eval_path=eval_path,
            selection_manifest_file=selection_manifest_file,
        )
    except (OSError, ValueError) as exc:
        return {
            "success": False,
            "error": f"invalid authoritative pose selection: {exc}",
        }
    docking_evidence_requested = bool(
        docking_scores_file
        or docking_report_file
        or any(
            _record_has_docking_evidence(row)
            for row in eval_results
            if isinstance(row, dict)
        )
    )
    if docking_scores_file:
        docking_scores_file = str(
            Path(docking_scores_file).expanduser().resolve(),
        )
    docking_report = {}
    docking_scores: list[dict] | None = (
        [] if docking_evidence_requested else None
    )
    docking_evidence_error: str | None = None
    docking_artifact_set_id: str | None = None
    if docking_report_file:
        docking_path = Path(docking_report_file).expanduser().resolve()
        docking_report_file = str(docking_path)
        verified_docking = load_committed_docking_artifacts(docking_path)
        docking_report = verified_docking.get("docking_report") or {}
        if not verified_docking.get("available"):
            docking_evidence_error = str(
                verified_docking.get("reason")
                or "docking artifact set is unavailable"
            )
        else:
            files = docking_report.get("files") or {}
            declared_value = files.get("docking_scores_json")
            declared_path = (
                Path(str(declared_value)).expanduser().resolve()
                if declared_value not in (None, "") else None
            )
            supplied_path = (
                Path(docking_scores_file).resolve()
                if docking_scores_file else None
            )
            if declared_path is None:
                docking_evidence_error = (
                    "committed docking report does not declare its score file"
                )
            elif supplied_path is not None and supplied_path != declared_path:
                docking_evidence_error = (
                    "docking_scores_file does not match the score artifact "
                    "sealed by docking_report_file"
                )
            else:
                lineage_error = evaluation_docking_lineage_error(
                    eval_report, docking_report,
                )
                if lineage_error is not None:
                    docking_evidence_error = lineage_error
                else:
                    docking_scores = list(
                        verified_docking.get("docking_scores") or []
                    )
                    docking_scores_file = str(declared_path)
                    contract = docking_report.get(
                        "docked_pose_artifact_contract"
                    ) or {}
                    docking_artifact_set_id = str(
                        contract.get("artifact_set_id") or ""
                    ) or None
    elif docking_evidence_requested:
        docking_evidence_error = (
            "docking evidence was supplied only through an uncommitted score "
            "file or downstream copy; docking_report_file is required"
        )

    # Load property filter results (optional)
    prop_results = []
    prop_report = {}
    if prop_file:
        prop_path = Path(prop_file).expanduser().resolve()
        prop_file = str(prop_path)
        if prop_path.exists():
            with open(prop_path) as f:
                prop_report = json.load(f)
            prop_results = prop_report.get("results", [])
    property_report_cohorts = _scientific_cohorts_from(prop_report)
    if (
        len(property_report_cohorts) == 1
        and _cohort_tuple_is_valid(property_report_cohorts[0])
    ):
        # Older property producers sometimes wrote the cohort once at report
        # level.  Project it only when the report is provably single-cohort;
        # doing this for a plural report would fabricate row membership.
        singular_property_cohort = property_report_cohorts[0]
        for row in prop_results:
            if not isinstance(row, dict):
                continue
            for field in SCIENTIFIC_COHORT_FIELDS:
                value = singular_property_cohort.get(field)
                if row.get(field) in (None, "") and value not in (None, ""):
                    row[field] = deepcopy(value)
    metric_rules = _effective_metric_rules(eval_report, prop_report)

    # Declare the steps this run is accountable for.  Coverage is measured
    # against these, so a step whose report was never supplied lowers coverage
    # instead of disappearing from both sides of the ratio.  All three are the
    # coordinator's premise: it ranks docked poses of property-filtered
    # molecules, so Step 4 is accountable even when its report was not wired --
    # deriving this from ``bool(prop_file)`` let an omitted argument narrow
    # accountability to 30 of 56 parent weight and still publish a `formal`
    # release.  Whether the report is *required* is a separate question, decided
    # by ``property_required`` below.
    active_steps: tuple[str, ...] = ("step4", "step5", "step6")

    # Attach canonical docking evidence to every evaluated pose before applying
    # the structure evaluator's authoritative post-QC selection.  The
    # coordinator must not choose a second, incompatible representative.
    prepared_eval_results = prepare_pose_records(
        eval_results,
        docking_scores,
        docking_evidence_unavailable_reason=docking_evidence_error,
    )
    if docking_artifact_set_id is not None:
        for entry in prepared_eval_results:
            if entry.get("docking_evidence_status") == "verified":
                entry["docking_artifact_set_id"] = docking_artifact_set_id
    evaluation_manifest = [
        item for item in _scientific_cohorts_from(eval_report, input_data)
        if _cohort_tuple_is_valid(item)
    ]
    if evaluation_manifest:
        candidate_cohort_ids = [
            str(item["cohort_id"]) for item in evaluation_manifest
        ]
        for entry in prepared_eval_results:
            compatible = any(
                _cohort_rows_are_compatible(entry, cohort)
                for cohort in evaluation_manifest
            )
            if compatible:
                continue
            has_scoped_identity = any(
                entry.get(field) not in (None, "")
                for field in (
                    "cohort_id", *SCIENTIFIC_COHORT_COMPONENT_FIELDS,
                )
            )
            if len(evaluation_manifest) == 1 and not has_scoped_identity:
                # A singular valid manifest can safely fill a legacy row.
                continue
            assignment_status = (
                "conflict" if has_scoped_identity else "ambiguous"
            )
            blocker = (
                "conflicting_scientific_cohort_assignment"
                if has_scoped_identity
                else "ambiguous_scientific_cohort_assignment"
            )
            entry["cohort_assignment_status"] = assignment_status
            entry["cohort_assignment_candidate_ids"] = candidate_cohort_ids
            context = entry.setdefault("metric_context", {})
            blockers = context.setdefault("protocol_blockers", [])
            if blocker not in blockers:
                blockers.append(blocker)
    for entry in prepared_eval_results:
        entry.update(_cohort_identity(
            entry, eval_report, round_number, input_data=input_data,
        ))
    step5_rules = tuple(
        rule for rule in metric_rules if rule.get("step_id") == "step5"
    )
    if authoritative_pose_ids is None and step5_rules and prepared_eval_results:
        # Legacy reports without a selection contract can choose by StepScore.
        # A modern manifest already binds the pose consumed by property_filter.
        try:
            prepared_eval_results, pose_selection_scoring = (
                _score_ranking_cohorts(
                    prepared_eval_results,
                    normalization_method=normalization_method,
                    metric_rules=step5_rules,
                    metric_weight_overrides={
                        metric_id: value
                        for metric_id, value in metric_weight_overrides.items()
                        if metric_id in {
                            str(rule.get("metric_id") or "")
                            for rule in step5_rules
                        }
                    },
                    category_weight_overrides=category_weight_overrides,
                    active_steps=("step5",),
                    screening_mode=screening_mode,
                )
            )
        except (TypeError, ValueError) as exc:
            return {"success": False, "error": f"invalid scoring policy: {exc}"}
        pose_selection_metadata["diagnostic_pose_scoring"] = (
            pose_selection_scoring
        )
    try:
        selected_eval_results = select_compound_poses(
            prepared_eval_results,
            authoritative_pose_ids=authoritative_pose_ids,
            eligible_verdicts=(
                set(pose_selection_metadata["eligible_verdicts"])
                if authoritative_pose_ids is not None
                else None
            ),
            authoritative_pose_records=authoritative_pose_records,
        )
    except ValueError as exc:
        return {
            "success": False,
            "error": f"invalid authoritative pose selection: {exc}",
        }

    evaluation_cohort_ids = {
        str(entry.get("cohort_id"))
        for entry in selected_eval_results
        if entry.get("cohort_id") not in (None, "")
    }
    property_cohort_ids = {
        str(entry.get("cohort_id"))
        for entry in prop_results
        if entry.get("cohort_id") not in (None, "")
    }
    # Legacy reports without row-level cohort provenance may be inferred only
    # inside a demonstrably single-cohort coordinator run.  A multi-cohort run,
    # or any report that already declares a different cohort, must fail closed.
    allow_legacy_property_cohort_inference = bool(
        len(evaluation_cohort_ids) == 1 and not property_cohort_ids
    )

    # Build cross-stage property lookups.  Stage-local row indices are
    # deliberately excluded: evaluator ``pose_index`` values refer to the
    # pre-QC docking-pose collection, while property ``record_index`` /
    # ``mol_index`` values refer to the post-QC selected-pose collection.  The
    # two integer domains can overlap for unrelated compounds.
    prop_by_name = {}
    prop_by_smiles = {}
    prop_by_chemical_identity = {}
    prop_by_chemical_id = {}
    prop_by_pose_id = {}
    for p in prop_results:
        chemical_identity_id = p.get("chemical_identity_id")
        if chemical_identity_id:
            _add_property_lookup(prop_by_chemical_identity, chemical_identity_id, p)
        chemical_id = p.get("chemical_id")
        if chemical_id:
            _add_property_lookup(prop_by_chemical_id, chemical_id, p)
        pose_id = p.get("pose_id")
        if pose_id:
            _add_property_lookup(prop_by_pose_id, pose_id, p)
        name = p.get("mol_name", "")
        if name:
            _add_property_lookup(prop_by_name, name, p)
        smiles = p.get("smiles")
        if smiles:
            _add_property_lookup(prop_by_smiles, smiles, p)

    # Classify each molecule
    # Intrinsic labels are populated from the current metric fold after scoring,
    # before coverage and experiment-quality blocks alter nomination status.
    intrinsic_classifications = []
    ranking_records = []
    for entry in selected_eval_results:
        mol_name = entry.get("mol_name", f"pose_{entry.get('pose_index', 0)}")
        chemical_identity_id = _chemical_identity_id(entry)
        prop_entry, property_join = _match_property_result(
            entry,
            by_chemical_identity=prop_by_chemical_identity,
            by_chemical_id=prop_by_chemical_id,
            by_pose_id=prop_by_pose_id,
            by_name=prop_by_name,
            by_smiles=prop_by_smiles,
            allow_legacy_cohort_inference=(
                allow_legacy_property_cohort_inference
            ),
        )
        property_review = property_data_review(
            prop_entry,
            required=bool(prop_file),
        )
        evaluation_review = evaluation_data_review(entry)
        # Pass the run's mode so V1 ("marking only") cannot inherit a
        # threshold-derived reject from a report written under V2/V3.
        cls = classify_molecule(
            entry,
            prop_entry,
            property_required=bool(prop_file),
            screening_mode=screening_mode,
        )
        blueprint_match = entry.get("blueprint_match") or {}
        capability_blocked = blueprint_match.get("verdict") == "unresolved"
        metric_context = deepcopy(entry.get("metric_context") or {})
        if evaluation_review is not None:
            protocol_blockers = list(
                metric_context.get("protocol_blockers") or []
            )
            blocker = "structure_evaluation_evidence_unavailable"
            if blocker not in protocol_blockers:
                protocol_blockers.append(blocker)
            metric_context["protocol_blockers"] = protocol_blockers
        record = {
            "mol_name": mol_name,
            "smiles": (
                entry.get("smiles")
                or (prop_entry.get("smiles", "") if prop_entry else "")
            ),
            "selected_pose_index": entry.get("pose_index"),
            "selection_rank": entry.get("selection_rank"),
            "source_pose_index": entry.get(
                "source_pose_index", entry.get("pose_index")
            ),
            "pose_id": entry.get("pose_id"),
            "chemical_identity_id": chemical_identity_id,
            "num_poses_considered": entry.get("num_poses_considered", 1),
            "pose_verdicts": entry.get("pose_verdicts", {}),
            "representative_pose_selection": entry.get(
                "representative_pose_selection"
            ),
            "selected_step5_score": entry.get("selected_step5_score"),
            "selected_step5_rank_status": entry.get(
                "selected_step5_rank_status"
            ),
            "structural_reject": record_has_structural_reject(entry) or record_has_structural_reject(prop_entry),
            "identity_status": entry.get("identity_status"),
            "evidence_status": entry.get("evidence_status"),
            "verdict_reasons": entry.get("verdict_reasons"),
            "legacy_classification": cls,
            "intrinsic_classification": cls,
            "classification": cls,
            "nomination_status": (
                "blocked_by_unresolved_blueprint"
                if capability_blocked
                else "blocked_by_structure_evaluation_protocol"
                if evaluation_review is not None
                else "blocked_by_property_protocol_evidence"
                if property_review is not None
                else "eligible"
            ),
            "eval_verdict": entry.get("verdict", "unknown"),
            "prop_verdict": prop_entry.get("verdict", "N/A") if prop_entry else "N/A",
            "modification_hint": (
                "Re-evaluate with a profiler capable of the required atom/role/geometry fields"
                if capability_blocked
                else evaluation_review["actions"][0]
                if evaluation_review is not None
                else entry.get("modification_hint", "") if cls == "modify" else ""
            ),
            "docking_scores": entry.get("docking_scores", {}),
            "docking_evidence_status": entry.get(
                "docking_evidence_status"
            ),
            "docking_score_join": entry.get("docking_score_join"),
            "docking_artifact_set_id": entry.get(
                "docking_artifact_set_id"
            ),
            "pose_scores": entry.get("pose_scores", []),
            "primary_score_metric_id": entry.get("primary_score_metric_id"),
            "primary_score": entry.get("primary_score"),
            "posecheck": entry.get("posecheck"),
            "posebusters": entry.get("posebusters"),
            "plip": entry.get("plip"),
            "prolif": entry.get("prolif"),
            "arpeggio": entry.get("arpeggio"),
            "rdkit_pose": entry.get("rdkit_pose"),
            "freesasa": entry.get("freesasa"),
            "mdanalysis": entry.get("mdanalysis"),
            "supplemental_interactions": entry.get("supplemental_interactions"),
            "blueprint_match": entry.get("blueprint_match"),
            "pocket_contact_fraction": entry.get("pocket_contact_fraction", entry.get("pocket_occupancy")),
            "metric_context": metric_context,
            "property": prop_entry,
            "property_join": property_join,
            "property_data_review": property_review,
            "evaluation_data_review": evaluation_review,
            "cohort_assignment_status": entry.get(
                "cohort_assignment_status", "verified"
            ),
            "cohort_assignment_candidate_ids": entry.get(
                "cohort_assignment_candidate_ids", []
            ),
        }
        record.update({
            field: deepcopy(entry.get(field))
            for field in LIGAND_LINEAGE_PROPERTIES
            if entry.get(field) not in (None, "")
        })
        record.update(_cohort_identity(
            entry, eval_report, round_number, input_data=input_data,
        ))
        exhaustion = entry.get("pose_exhaustion")
        if exhaustion is not None and evaluation_review is None:
            record["pose_exhaustion"] = deepcopy(exhaustion)
        ranking_records.append(record)

    # Fit one normalization profile for the nomination cohort.  DEC_* records
    # are scored in a separate benchmark cohort so adding retrospective decoys
    # cannot change candidate normalization or rank.  Raw docking scores retain
    # their single-metric contract; only normalized metric values are combined.
    try:
        ranking_records, scoring_metadata = _score_ranking_cohorts(
            ranking_records,
            normalization_method=normalization_method,
            metric_rules=metric_rules,
            metric_weight_overrides=metric_weight_overrides,
            category_weight_overrides=category_weight_overrides,
            active_steps=active_steps,
            screening_mode=screening_mode,
        )
    except (TypeError, ValueError) as exc:
        return {"success": False, "error": f"invalid scoring policy: {exc}"}
    for record in ranking_records:
        # Current category folds own the molecular decision. The shared gate
        # already includes explicit structural failures; legacy/coverage labels
        # remain diagnostics and cannot override new metric evidence.
        record["classification"] = record.get("metric_gate", "not_evaluated")
        record["intrinsic_classification"] = record["classification"]
        intrinsic_classifications.append(record["intrinsic_classification"])
        if record.get("pose_exhaustion") and record["classification"] == "reject":
            record["reason_codes"] = list(dict.fromkeys([
                *(record.get("ranking_blocked_by") or []), "all_poses_failed",
            ]))
            record["modification_hint"] = record["pose_exhaustion"]["reason"]
        if record["classification"] == "not_evaluated" and not str(
            record.get("nomination_status") or "eligible"
        ).startswith("blocked"):
            # A molecule no metric rule could label has not been shown to be
            # good OR bad.  It is not rejected, but it must not be
            # nominated on the strength of the intrinsic label alone.
            record["nomination_status"] = "blocked_by_missing_metric_evidence"
        if record.get("ranking_blocked"):
            # Keep the first, more specific block reason if one was already
            # recorded upstream; only annotate when nothing has claimed it.
            if not str(
                record.get("nomination_status") or "eligible"
            ).startswith("blocked"):
                record["nomination_status"] = "blocked_by_metric_protocol"
            record["modification_hint"] = (
                (record.get("modification_hint") or "")
                + "; Re-evaluate metrics that block trusted ranking: "
                + ", ".join(record.get("ranking_blocked_by") or [])
            ).strip("; ")
        else:
            _apply_coverage_consequences(record)
    # These molecules never reached pose evaluation. Keep them in the audit and
    # benchmark denominator, but outside normalization and candidate ranking.
    if not docking_evidence_error:
        ranking_records.extend(_docking_failure_records(
            docking_report, round_number=round_number, input_data=input_data,
        ))
    classifications = [record["classification"] for record in ranking_records]

    num_keep = classifications.count("keep")
    num_modify = classifications.count("modify")
    num_reject = classifications.count("reject")
    # Counted and reported separately so keep+modify+reject no longer has to
    # equal ``total``: a molecule no metric could label is neither advanced nor
    # rejected, and folding it into any of the three would misstate the round.
    num_not_evaluated = classifications.count("not_evaluated")
    total = len(classifications)
    pass_rate = num_keep / total if total > 0 else 0.0

    # Benchmark decoys remain available for AUROC/EF, but never drive molecule
    # optimization, bottleneck diagnosis, or next-round scaffold exclusions.
    candidate_pairs = [
        (entry, record)
        for entry, record in zip(selected_eval_results, ranking_records)
        if not _is_decoy_record(record)
    ]
    candidate_eval_results = [entry for entry, _record in candidate_pairs]
    candidate_classifications = [
        record["classification"] for _entry, record in candidate_pairs
    ]
    candidate_property_results = [
        record["property"]
        for _entry, record in candidate_pairs
        if record.get("property") is not None
    ]

    # Identify the optimization-cohort bottleneck.
    primary_bottleneck = identify_bottleneck(
        candidate_eval_results,
        candidate_property_results,
    )

    quality_gates = _quality_gates(docking_report_file, pocket_analysis_report_file)
    if docking_evidence_error:
        quality_gates["docking_artifact_identity"] = {
            "passed": None,
            "not_applicable": False,
            "gate": "modify",
            "reason": docking_evidence_error,
            "protocol_blocker": True,
        }
    benchmark = screening_metrics(ranking_records) if screening_benchmark_enabled else {
        "available": False, "reason": "screening_benchmark disabled_by_config",
    }
    if benchmark.get("available"):
        benchmark_gate = _screening_benchmark_gate(
            benchmark,
            min_actives=min_validation_actives,
            min_decoys=min_validation_decoys,
            auroc_soft=min_validation_auroc,
            auroc_hard=min_validation_auroc_hard,
            ef1_soft=min_validation_ef1,
            ef1_hard=min_validation_ef1_hard,
        )
        quality_gates["screening_benchmark"] = {
            "gate": benchmark_gate["gate"],
            "passed": benchmark_gate["gate"] == "keep",
            "reason": benchmark_gate["reason"],
            "details": benchmark,
            "criteria": {
                "min_actives": min_validation_actives,
                "min_decoys": min_validation_decoys,
                "AUROC_hard_lower": min_validation_auroc_hard,
                "AUROC_soft_lower": min_validation_auroc,
                "EF1pct_hard_lower": min_validation_ef1_hard,
                "EF1pct_soft_lower": min_validation_ef1,
            },
            "metric_assessments": benchmark_gate["metric_assessments"],
        }
    elif screening_benchmark_enabled:
        quality_gates["screening_benchmark"] = {
            "gate": "modify",
            "passed": False,
            "reason": benchmark.get(
                "reason",
                "labeled screening benchmark is incomplete",
            ),
            "details": benchmark,
            "criteria": {
                "requires_both_labeled_classes": True,
                "requires_complete_labeled_score_coverage": True,
                "requires_one_complete_shared_docking_protocol": True,
            },
        }
    quality_gate_status = _aggregate_quality_gate_status(quality_gates)
    quality_gate_passed = quality_gate_status == "keep"
    quality_gate_statuses = {
        name: _quality_gate_status(value)
        for name, value in quality_gates.items()
    }
    failed_quality_gates = [
        name
        for name, status in quality_gate_statuses.items()
        if status in {"modify", "reject"}
    ]
    not_applicable_quality_gates = [
        name
        for name, status in quality_gate_statuses.items()
        if status == "not_applicable"
    ]
    quality_gate_blockers = list(failed_quality_gates)
    if quality_gate_status == "not_applicable":
        quality_gate_blockers.extend(not_applicable_quality_gates)
    intrinsic_num_keep = intrinsic_classifications.count("keep")
    intrinsic_num_modify = intrinsic_classifications.count("modify")
    intrinsic_num_reject = intrinsic_classifications.count("reject")
    intrinsic_pass_rate = intrinsic_num_keep / total if total > 0 else 0.0

    # Protocol validation is orthogonal to the molecule gate.  Preserve the
    # intrinsic keep/modify/reject result and block nomination separately.
    if require_quality_gates and quality_gate_status != "keep":
        for record in ranking_records:
            record["ranking_blocked"] = True
            ranking_blocked_by = list(record.get("ranking_blocked_by") or [])
            ranking_blocked_by.extend(
                f"quality_gate:{name}" for name in quality_gate_blockers
            )
            record["ranking_blocked_by"] = list(
                dict.fromkeys(ranking_blocked_by)
            )
            record["weighted_score_rank"] = None
            record["overall_rank"] = None
            record["blocked_quality_gates"] = quality_gate_blockers
            if record["classification"] == "keep":
                if not str(
                    record.get("nomination_status") or "eligible"
                ).startswith("blocked"):
                    record["nomination_status"] = "blocked_by_quality_gate"
                record["protocol_block_reason"] = (
                    "Resolve experiment quality gates before hit nomination"
                )

    nominated_records = _candidate_records(ranking_records)
    num_candidates = len(nominated_records)
    num_candidate_keep = sum(
        1 for r in nominated_records if _eligible_keep(r)
    )
    num_candidate_modify = sum(1 for r in nominated_records if r.get("classification") == "modify")
    num_candidate_reject = sum(1 for r in nominated_records if r.get("classification") == "reject")
    num_candidate_not_evaluated = sum(
        1 for r in nominated_records
        if r.get("classification") == "not_evaluated"
    )
    candidate_pass_rate = num_candidate_keep / num_candidates if num_candidates > 0 else 0.0
    num_candidate_intrinsic_keep = sum(1 for r in nominated_records if r.get("intrinsic_classification") == "keep")
    num_candidate_intrinsic_modify = sum(1 for r in nominated_records if r.get("intrinsic_classification") == "modify")
    num_candidate_intrinsic_reject = sum(1 for r in nominated_records if r.get("intrinsic_classification") == "reject")
    candidate_intrinsic_pass_rate = num_candidate_intrinsic_keep / num_candidates if num_candidates > 0 else 0.0
    candidate_ranking_blocked = sum(
        1 for r in nominated_records if r.get("ranking_blocked") is True
    )
    candidate_metric_blocked = [
        r
        for r in nominated_records
        if r.get("classification") != "reject"
        and any(
            not str(reason).startswith("quality_gate:")
            for reason in (r.get("ranking_blocked_by") or [])
        )
    ]
    candidate_local_modify = [
        record
        for record in nominated_records
        if record.get("classification") == "modify"
        and record.get("intrinsic_classification") == "modify"
        and record.get("property_data_review") is None
        and record.get("evaluation_data_review") is None
    ]
    candidate_property_data_review = [
        r
        for r in nominated_records
        if r.get("classification") != "reject"
        and isinstance(r.get("property_data_review"), dict)
    ]
    blocked_keep = sum(
        1
        for r in ranking_records
        if r.get("classification") == "keep" and r.get("blocked_quality_gates")
    )
    candidate_blocked_keep = sum(
        1
        for r in nominated_records
        if r.get("classification") == "keep" and r.get("blocked_quality_gates")
    )

    # Round validity describes the evidence, never whether a campaign should stop.
    round_assessment = _assess_round(
        nominated_records,
        quality_gate_status=quality_gate_status,
        require_quality_gates=require_quality_gates,
        failed_quality_gates=failed_quality_gates,
        not_applicable_quality_gates=not_applicable_quality_gates,
    )
    valid_for_nomination = round_assessment["valid_for_nomination"]
    decision = "agent_review"
    decision_outcome = None
    goal_met = None
    rationale = "Round evidence is available for Agent campaign planning."
    rejected_quality_gates = [
        name for name, status in quality_gate_statuses.items() if status == "reject"
    ]
    rollback_instructions = (
        {
            "required": True,
            "scope": "protocol",
            "action": "rollback_invalid_protocol",
            "invalid_quality_gates": rejected_quality_gates,
            "preserve_molecule_evidence": True,
            "automatic": False,
            "instructions": [
                "Invalidate artifacts produced by the rejected protocol",
                "Return to the earliest affected protocol step before rerunning",
                "Do not reinterpret protocol failure as a molecule defect",
            ],
        }
        if round_assessment["status"] == "invalid" else None
    )

    # Diagnostic hints are evidence for the Agent, never an execution decision.
    diagnostics = {}
    if (
        require_quality_gates
        and quality_gate_status in {"modify", "reject", "not_applicable"}
    ):
        if quality_gate_status == "not_applicable":
            suggestions = [
                f"Acquire applicable validation evidence or explicitly revise "
                f"the protocol requirement: {name}"
                for name in not_applicable_quality_gates
            ]
        else:
            suggestions = [
                f"Rerun or repair protocol quality gate: {name}"
                for name in failed_quality_gates
            ]
        diagnostics = {
            "target_molecules": [],
            "suggestions": suggestions,
            "diagnosis": {
                "bottleneck_type": "protocol_level",
                "failed_scaffolds": [],
                "modification_hints": suggestions,
                "suggested_strategy": "re_evaluate",
            },
        }
    elif candidate_property_data_review:
        review_actions = list(dict.fromkeys(
            action
            for record in candidate_property_data_review
            for action in (record["property_data_review"].get("actions") or [])
        ))
        review_reasons = {
            record["mol_name"]: record["property_data_review"].get("reasons") or []
            for record in candidate_property_data_review
        }
        diagnostics = {
            "target_molecules": [],
            "affected_molecules": [
                record["mol_name"] for record in candidate_property_data_review
            ],
            "suggestions": review_actions,
            "property_data_review": {
                "reasons_by_molecule": review_reasons,
                "actions": review_actions,
            },
            "diagnosis": {
                "bottleneck_type": "protocol_level",
                "review_type": "property_data_review",
                "failed_scaffolds": [],
                "modification_hints": review_actions,
                "suggested_strategy": "re_evaluate_property_evidence",
            },
        }
    elif (
        candidate_metric_blocked
        and not candidate_local_modify
    ):
        metric_blockers = list(dict.fromkeys(
            reason
            for record in candidate_metric_blocked
            for reason in (record.get("ranking_blocked_by") or [])
            if not str(reason).startswith("quality_gate:")
        ))
        suggestions = [
            f"Re-evaluate metric/protocol blocker: {reason}"
            for reason in metric_blockers
        ]
        diagnostics = {
            "target_molecules": [],
            "affected_molecules": [
                record["mol_name"] for record in candidate_metric_blocked
            ],
            "suggestions": suggestions,
            "diagnosis": {
                "bottleneck_type": "protocol_level",
                "failed_scaffolds": [],
                "modification_hints": suggestions,
                "suggested_strategy": "re_evaluate",
            },
        }
    else:
        diagnostics = aggregate_modification_hints(
            candidate_eval_results,
            candidate_classifications,
        )
        scaffold_level = is_scaffold_level(primary_bottleneck)
        protocol_level = is_protocol_level(primary_bottleneck)
        # Extract Murcko scaffolds from rejected/modify molecules for exclusion
        failed_scaffolds = _extract_failed_scaffolds(
            candidate_eval_results,
            candidate_classifications,
        )
        diagnostics["diagnosis"] = {
            "bottleneck_type": (
                "protocol_level" if protocol_level
                else "scaffold_level" if scaffold_level
                else "local_level"
            ),
            "failed_scaffolds": failed_scaffolds,
            "modification_hints": diagnostics.get("suggestions", []),
            "suggested_strategy": (
                "re_evaluate" if protocol_level
                else "regenerate" if scaffold_level
                else "modify"
            ),
        }

        if candidate_metric_blocked:
            metric_blockers = list(dict.fromkeys(
                reason
                for record in candidate_metric_blocked
                for reason in (record.get("ranking_blocked_by") or [])
                if not str(reason).startswith("quality_gate:")
            ))
            diagnostics["parallel_protocol_review"] = {
                "affected_molecules": [
                    record["mol_name"] for record in candidate_metric_blocked
                ],
                "blockers": metric_blockers,
                "actions": [
                    f"Re-evaluate metric/protocol blocker: {reason}"
                    for reason in metric_blockers
                ],
                "preserve_molecule_modification_targets": True,
            }

    # Protocol/data-review work never asks a molecule generator or modifier to
    # change chemistry. Preserve the impacted IDs separately for audit/reruns.
    if (
        diagnostics
        and (diagnostics.get("diagnosis") or {}).get("bottleneck_type")
        == "protocol_level"
    ):
        affected = list(diagnostics.get("affected_molecules") or [])
        affected.extend(diagnostics.get("target_molecules") or [])
        diagnostics["target_molecules"] = []
        if affected:
            diagnostics["affected_molecules"] = list(dict.fromkeys(affected))

    diagnostics["advisory_only"] = True

    # Compute fix_rate from prior round history
    fix_rate = None
    if prior_history:
        prev_round = prior_history[-1]
        prev_molecules = {
            m["mol_name"]: m["classification"]
            for m in prev_round.get("molecules", [])
            if not _is_decoy_record(m)
        }
        fixed = sum(
            1 for rec in ranking_records
            if not _is_decoy_record(rec)
            and rec["mol_name"] in prev_molecules
            and prev_molecules[rec["mol_name"]] == "modify"
            # Fix-rate measures whether the molecule-level issue was repaired;
            # a separate protocol/coverage blocker may still prevent formal
            # nomination and must not turn a chemical fix back into zero.
            and rec.get("classification") == "keep"
        )
        prev_modify = sum(1 for c in prev_molecules.values() if c == "modify")
        fix_rate = round(fixed / prev_modify, 3) if prev_modify > 0 else 0.0

    # Build report
    round_summary = {
        "round": round_number,
        "assessment_schema_version": "round-assessment/v1",
        "round_assessment": round_assessment,
        "total_evaluated": total,
        "keep": num_keep,
        "modify": num_modify,
        "reject": num_reject,
        "not_evaluated": num_not_evaluated,
        "candidate_total": num_candidates,
        "candidate_keep": num_candidate_keep,
        "candidate_modify": num_candidate_modify,
        "candidate_reject": num_candidate_reject,
        "candidate_not_evaluated": num_candidate_not_evaluated,
        "pass_rate": round(pass_rate, 3),
        "candidate_pass_rate": round(candidate_pass_rate, 3),
        "intrinsic_keep": intrinsic_num_keep,
        "intrinsic_modify": intrinsic_num_modify,
        "intrinsic_reject": intrinsic_num_reject,
        "intrinsic_pass_rate": round(intrinsic_pass_rate, 3),
        "candidate_intrinsic_keep": num_candidate_intrinsic_keep,
        "candidate_intrinsic_modify": num_candidate_intrinsic_modify,
        "candidate_intrinsic_reject": num_candidate_intrinsic_reject,
        "candidate_intrinsic_pass_rate": round(candidate_intrinsic_pass_rate, 3),
        "quality_gate_blocked_keep": blocked_keep,
        "candidate_quality_gate_blocked_keep": candidate_blocked_keep,
        "candidate_ranking_blocked": candidate_ranking_blocked,
        "failed_quality_gates": failed_quality_gates,
        "not_applicable_quality_gates": not_applicable_quality_gates,
        "candidate_property_data_review": len(candidate_property_data_review),
        "fix_rate": fix_rate,
        "primary_bottleneck": primary_bottleneck,
        "quality_gate_passed": quality_gate_passed,
        "quality_gate_status": quality_gate_status,
        "goal_met": goal_met,
        "decision_outcome": decision_outcome,
        "valid_for_nomination": valid_for_nomination,
    }

    # Publish RANK.AVG after all nomination blocks have been applied, sorting
    # the full audit and non-decoy nomination sets independently.
    ranked_molecules = rank_candidates(deepcopy(ranking_records))
    # Gate first, then score: a rejected molecule does not advance, so it has no
    # place in the ranking.  It stays in ``molecules`` and in the published
    # tables, which are the audit trail for why it was dropped.
    ranked_candidates = rank_candidates(deepcopy([
        record
        for record in nominated_records
        if record.get("classification") != "reject" and not record.get("audit_only")
    ]))
    report = {
        "success": True,
        "assessment_schema_version": "round-assessment/v1",
        "round_assessment": round_assessment,
        "round_summary": round_summary,
        "molecules": ranked_molecules,
        "candidate_molecules": ranked_candidates,
        "ranking_contract": {
            "rank": "local_overall_rank_or_null",
            "overall_rank": (
                "rank_avg_with_ties_within_same_cohort_and_normalization_pool"
            ),
            "weighted_score_rank": "deprecated_alias_of_overall_rank",
            "gate_or_blocker_affects_rank": True,
            "candidate_ranking_scope": "non_decoy_nomination_candidates",
            "full_molecule_ranking_scope": "all_evaluated_records",
            "candidate_and_decoy_normalization_separate": True,
            "physical_order": (
                "cohort_then_normalization_pool_then_local_overall_rank"
            ),
            "cross_cohort_rank_comparison": False,
            "top_n_scope": "per_cohort_and_normalization_pool",
            "mixed_schema_policy": (
                "schema1_rows_are_unranked_audit_rows_when_schema2_is_present"
            ),
        },
        "representative_pose_selection": pose_selection_metadata,
        "metric_scoring": scoring_metadata,
        "metric_rules_version": METRIC_RULES_VERSION,
        "metric_rules": public_metric_rules(metric_rules),
        "protocol_metric_rules": public_metric_rules(
            DOCKING_PROTOCOL_METRIC_RULES
        ),
        "quality_gates": quality_gates,
        "quality_gate_status": quality_gate_status,
        "screening_validation": benchmark,
        "decision": decision,
        "decision_outcome": decision_outcome,
        "goal_met": goal_met,
        "valid_for_nomination": valid_for_nomination,
        "rollback_instructions": rollback_instructions,
        "decision_rationale": rationale,
        "next_round_instructions": None,
        "diagnostics": diagnostics,
        "parameters": {
            "round_number": round_number,
            "deprecated_campaign_parameters": {
                "ignored": True,
                "max_rounds": max_rounds,
                "pass_rate_threshold": pass_rate_threshold,
                "min_candidates": min_candidates,
            },
            "max_rounds": max_rounds,
            "pass_rate_threshold": pass_rate_threshold,
            "min_candidates": min_candidates,
            "require_quality_gates": require_quality_gates,
            "min_validation_actives": min_validation_actives,
            "min_validation_decoys": min_validation_decoys,
            "min_validation_auroc": min_validation_auroc,
            "min_validation_auroc_hard": min_validation_auroc_hard,
            "min_validation_ef1": min_validation_ef1,
            "min_validation_ef1_hard": min_validation_ef1_hard,
            "normalization_method": normalization_method,
            "gate_contract": "two_level_label_then_category_fold",
            "screening_mode": dict(screening_mode),
            "screening_modes": {
                mode_id: dict(spec) for mode_id, spec in SCREENING_MODES.items()
            },
            "metric_weight_overrides": dict(metric_weight_overrides),
            "category_weight_overrides": dict(category_weight_overrides),
            **{
                field: input_data.get(field)
                for field in (
                    "cohort_id", "protein_state_id", "pocket_id", "round_id",
                    "protocol_profile_id", "metric_catalog_version",
                    "cohort_identity_complete",
                    "cohort_identity_missing_fields",
                    "cohort_identity_derived_fields",
                    "cohort_provenance_sources",
                )
                if input_data.get(field) not in (None, "")
            },
        },
    }
    cohort_manifest: dict[str, dict] = {}
    upstream_cohorts = _scientific_cohorts_from(
        input_data, eval_report, prop_report, docking_report,
    )
    cohort_rows = [*upstream_cohorts, *ranking_records] or [
        _cohort_identity({}, eval_report, round_number, input_data=input_data)
    ]
    for row in cohort_rows:
        identity = {
            field: row.get(field)
            for field in SCIENTIFIC_COHORT_FIELDS
        }
        fingerprint = json.dumps(
            {
                field: identity.get(field)
                for field in ("cohort_id", *SCIENTIFIC_COHORT_COMPONENT_FIELDS)
            },
            sort_keys=True,
            default=str,
        )
        cohort_manifest.setdefault(fingerprint, identity)
    report["scientific_cohorts"] = list(cohort_manifest.values())

    # Build cross-round history snapshot for this round
    round_history_entry = {
        "round": round_number,
        "assessment_schema_version": "round-assessment/v1",
        "round_assessment": round_assessment,
        "decision": decision,
        "decision_outcome": decision_outcome,
        "goal_met": goal_met,
        "valid_for_nomination": valid_for_nomination,
        "primary_bottleneck": primary_bottleneck,
        "keep": num_keep,
        "modify": num_modify,
        "reject": num_reject,
        "not_evaluated": num_not_evaluated,
        "pass_rate": round(pass_rate, 3),
        "fix_rate": fix_rate,
        "molecules": [
            {
                "rank": r.get("rank"),
                "mol_name": r["mol_name"],
                "compound_id": r.get("compound_id"),
                "input_mol_name": r.get("input_mol_name"),
                "smiles": r.get("smiles", ""),
                "selected_pose_index": r.get("selected_pose_index"),
                "selection_rank": r.get("selection_rank"),
                "source_pose_index": r.get("source_pose_index"),
                "pose_id": r.get("pose_id"),
                "chemical_identity_id": r.get("chemical_identity_id"),
                "docking_artifact_set_id": r.get(
                    "docking_artifact_set_id"
                ),
                "cohort_id": r.get("cohort_id"),
                "protein_state_id": r.get("protein_state_id"),
                "pocket_id": r.get("pocket_id"),
                "round_id": r.get("round_id"),
                "protocol_profile_id": r.get("protocol_profile_id"),
                "metric_catalog_version": r.get("metric_catalog_version"),
                "cohort_identity_complete": r.get(
                    "cohort_identity_complete"
                ),
                "cohort_identity_missing_fields": r.get(
                    "cohort_identity_missing_fields"
                ),
                "cohort_identity_derived_fields": r.get(
                    "cohort_identity_derived_fields"
                ),
                "cohort_provenance_sources": r.get(
                    "cohort_provenance_sources"
                ),
                "normalization_pool": r.get("normalization_pool"),
                "rank_partition_id": r.get("rank_partition_id"),
                "rank_scope": r.get("rank_scope"),
                "ranking_schema": r.get("ranking_schema"),
                "scoring_schema_version": r.get("scoring_schema_version"),
                "intrinsic_classification": r.get("intrinsic_classification", r["classification"]),
                "legacy_classification": r.get("legacy_classification"),
                "classification": r["classification"],
                "nomination_status": r.get("nomination_status", "eligible"),
                "eval_verdict": r["eval_verdict"],
                "blueprint_match": r.get("blueprint_match"),
                "posecheck": r.get("posecheck"),
                "docking_scores": r.get("docking_scores"),
                "pose_scores": r.get("pose_scores"),
                "primary_score_metric_id": r.get("primary_score_metric_id"),
                "primary_score": r.get("primary_score"),
                "metric_gate": r.get("metric_gate"),
                "coverage_status": r.get("coverage_status"),
                "rank_status": r.get("rank_status"),
                "overall_score": r.get("overall_score"),
                "overall_rank": r.get("overall_rank"),
                "category_weighted_coverage": r.get(
                    "category_weighted_coverage"
                ),
                "step_scores": r.get("step_scores"),
                "category_scores": r.get("category_scores"),
                "weighted_total_score": r.get("weighted_total_score"),
                "weighted_score_fraction": r.get("weighted_score_fraction"),
                "weighted_score_rank": r.get("weighted_score_rank"),
                "score_coverage": r.get("score_coverage"),
                "ranking_blocked": r.get("ranking_blocked", False),
                "ranking_blocked_by": r.get("ranking_blocked_by", []),
                "modification_hint": r.get("modification_hint", ""),
                "property": {
                    "qed": _property_value(r["property"], "qed", "QED"),
                    "sa_score": _property_value(r["property"], "sa_score", "SA_score"),
                    "verdict": r["property"].get("verdict"),
                    "flags": r["property"].get("flags"),
                } if r.get("property") else None,
            }
            for r in report["molecules"]
        ],
    }
    full_history = {"rounds": prior_history + [round_history_entry]}
    report["round_history"] = full_history

    # Add fix_history to each molecule for the round report
    mol_history: dict[str, list[str]] = {}
    for rnd in full_history["rounds"]:
        for m in rnd.get("molecules", []):
            name = m.get("mol_name", "")
            cls = m.get("classification", "")
            hint = m.get("modification_hint", "")
            entry = f"R{rnd['round']}: {cls}"
            if hint:
                entry += f" ({hint})"
            mol_history.setdefault(name, []).append(entry)
    for rec in report["molecules"]:
        rec["fix_history"] = "; ".join(mol_history.get(rec.get("mol_name", ""), ["-"]))
    for rec in report["candidate_molecules"]:
        rec["fix_history"] = "; ".join(mol_history.get(rec.get("mol_name", ""), ["-"]))

    report_file = output_dir / "coordinator_report.json"

    summary_file = output_dir / "round_summary.json"
    summary_file.write_text(json.dumps(round_summary, indent=2, ensure_ascii=False))

    history_file = output_dir / "round_history.json"
    with history_file.open("w", encoding="utf-8") as handle:
        json.dump(full_history, handle, indent=2, ensure_ascii=False, default=str)

    cohort_manifest_file = output_dir / "scientific_cohorts.json"
    cohort_manifest_file.write_text(json.dumps(
        report.get("scientific_cohorts") or [],
        indent=2,
        ensure_ascii=False,
        default=str,
    ))

    ranking_json = output_dir / "candidate_ranking.json"
    with ranking_json.open("w", encoding="utf-8") as handle:
        json.dump(
            [_compact_ranking_record(row, prop_file)
             for row in report["candidate_molecules"]],
            handle, indent=2, ensure_ascii=False, default=str,
        )
    # Keep the compact CSV compatibility artifact alongside the richer
    # Reader-facing TSV tree.  Downstream consumers use this table to inspect a
    # provider-neutral primary docking score without relabeling it as a GNINA
    # CNN field.
    ranking_csv = output_dir / "candidate_ranking.csv"
    _write_ranking_csv(ranking_csv, report["candidate_molecules"])
    round_report = output_dir / "round_report.md"
    _write_round_report(round_report, report)

    report["files"] = {
        "coordinator_report_json": str(report_file),
        "round_summary_json": str(summary_file),
        "round_history_json": str(history_file),
        "scientific_cohorts_json": str(cohort_manifest_file),
        "candidate_ranking_json": str(ranking_json),
        "candidate_ranking_csv": str(ranking_csv),
        "round_report_md": str(round_report),
    }
    # The pipeline owns final publication.  It can only seal the reader-facing
    # view after this result and Step 7's reporting result have both been
    # written with their artifact fingerprints.
    report["publication_status"] = "deferred_to_pipeline_finalization"
    with report_file.open("w", encoding="utf-8") as handle:
        json.dump(
            _compact_coordinator_report(report, prop_file),
            handle, indent=2, ensure_ascii=False, default=str,
        )

    return report



def _compact_ranking_record(record: dict, property_report_file: str | None) -> dict:
    """Keep round scores and raw evidence; reference the upstream scoring audit.

    The property stage's metric/category/step scores are not the coordinator's
    scores: they use a different normalization pool. Preserve that source report
    rather than embedding its scoring matrices and aliases again in every row.
    Shallow copies deliberately leave the scientific in-memory result untouched.
    """
    compact = dict(record)
    properties = record.get("property")
    if isinstance(properties, dict) and property_report_file:
        compact["property"] = {
            key: value for key, value in properties.items()
            if key not in {
                "metric_assessments", "category_scores", "step_scores",
                "property_category_scores", "property_step_scores",
            }
        }
        compact["property"]["scoring_report_file"] = str(
            Path(property_report_file).expanduser().resolve()
        )
    return compact


def _compact_coordinator_report(report: dict, property_report_file: str | None) -> dict:
    """Store full molecule evidence once, with lightweight candidate summaries.

    Candidate order and rank remain those of the independently ranked candidate
    pool. Rich candidate rows and cross-round history already have their own
    artifacts, referenced by ``files`` in this report.
    """
    compact = {key: value for key, value in report.items() if key != "round_history"}
    compact["molecules"] = [
        _compact_ranking_record(row, property_report_file)
        for row in report.get("molecules") or []
    ]
    candidate_fields = (
        "mol_name", "compound_id", "input_mol_name", "smiles", "pose_id",
        "chemical_identity_id", "selected_pose_index", "source_pose_index",
        "selection_rank", "docking_artifact_set_id", *SCIENTIFIC_COHORT_FIELDS,
        "normalization_pool", "rank_partition_id", "rank_scope",
        "ranking_schema", "scoring_schema_version", "rank", "overall_rank",
        "overall_score", "weighted_score_rank", "weighted_total_score",
        "weighted_score_fraction", "category_weighted_coverage",
        "classification", "intrinsic_classification", "nomination_status",
        "metric_gate", "gate_status", "coverage_status", "rank_status",
        "ranking_blocked", "ranking_blocked_by", "audit_only",
        "eval_verdict", "prop_verdict", "modification_hint", "fix_history",
        "step4_score", "step4_rank", "step5_score", "step5_rank",
        "step6_score", "step6_rank",
    )
    compact["candidate_molecules"] = [
        {key: row[key] for key in candidate_fields if key in row}
        for row in report.get("candidate_molecules") or []
    ]
    compact["report_storage"] = {
        "schema_version": 1,
        "candidate_molecules": "summaries; full rows in files.candidate_ranking_json",
        "round_history": "files.round_history_json",
        "property_scoring": "molecules[].property.scoring_report_file",
    }
    return compact


def _docking_failure_records(docking_report: dict, *, round_number: int, input_data: dict) -> list[dict]:
    """Preserve molecule-local execution failures without inventing pose evidence."""
    if docking_report.get("success") is not True:
        return []
    records = []
    for row in docking_report.get("molecule_ledger") or []:
        if row.get("docking_status") not in {"failed", "input_invalid"}:
            continue
        record = {
            **{key: row[key] for key in LIGAND_LINEAGE_PROPERTIES if row.get(key) is not None},
            "mol_name": row.get("mol_name") or row.get("compound_id"),
            "compound_id": row.get("compound_id"),
            "source_molecule_index": row.get("source_molecule_index"),
            "chemical_identity_id": row.get("chemical_identity_id"),
            "chemical_id": row.get("chemical_id"),
            "smiles": row.get("canonical_smiles") or "",
            "docking_status": row["docking_status"],
            "failure_reason": row.get("failure_reason"),
            "retryable": row.get("retryable", True),
            "audit_only": True,
            "classification": "not_evaluated",
            "intrinsic_classification": "not_evaluated",
            "metric_gate": "not_evaluated",
            "nomination_status": "blocked_by_docking_failure",
            "ranking_blocked": True,
            "ranking_blocked_by": ["docking_molecule_failed"],
            "docking_evidence_status": "unavailable",
            "eval_verdict": "not_evaluated",
            "prop_verdict": "not_evaluated",
            "pose_scores": [], "docking_scores": {}, "property": None,
            "pose_id": None, "selected_pose_index": None, "source_pose_index": None,
            "primary_score": None, "overall_score": None, "overall_rank": None,
            "rank": None,
        }
        record.update(_cohort_identity(
            row, docking_report, round_number, input_data=input_data,
        ))
        records.append(record)
    return records


def _assess_round(
    records: list[dict],
    *,
    quality_gate_status: str,
    require_quality_gates: bool,
    failed_quality_gates: list[str],
    not_applicable_quality_gates: list[str],
) -> dict:
    """Keep deterministic evidence validity independent of campaign policy."""
    candidates = _candidate_records(records)
    eligible = [record for record in candidates if _eligible_keep(record)]
    protocol_valid = not require_quality_gates or quality_gate_status == "keep"
    unresolved = [
        record for record in candidates
        if record.get("classification") != "reject"
        and (record.get("ranking_blocked")
             or record.get("property_data_review")
             or record.get("evaluation_data_review")
             or record.get("classification") == "not_evaluated")
    ]
    evaluated = [
        record for record in candidates
        if record.get("classification") in {"keep", "modify", "reject"}
    ]
    if require_quality_gates and quality_gate_status == "reject":
        status = "invalid"
    elif not evaluated:
        status = "not_evaluated"
    elif not protocol_valid or unresolved:
        status = "review_required"
    else:
        status = "valid"
    return {
        "status": status,
        "valid_for_nomination": bool(eligible) and protocol_valid,
        "quality_gates_required": require_quality_gates,
        "protocol_valid": protocol_valid,
        "quality_gate_status": quality_gate_status,
        "failed_quality_gates": list(failed_quality_gates),
        "not_applicable_quality_gates": list(not_applicable_quality_gates),
        "candidate_total": len(candidates),
        "evaluated_candidates": len(evaluated),
        "eligible_keep_count": len(eligible),
        "review_required_molecules": [record.get("mol_name") for record in unresolved],
        "campaign_decision_owner": "agent",
    }


def _load_scores(path_value: str | None) -> list[dict]:
    from io_utils import load_scores
    return load_scores(path_value)


def _record_has_docking_evidence(record: dict) -> bool:
    """Detect score copies that require a committed Step 3 source."""
    if isinstance(record.get("docking_scores"), dict):
        return True
    if any(record.get(field) is not None for field in LEGACY_SCORE_SPECS):
        return True
    docking_metric_ids = {
        str(spec.get("metric_id"))
        for spec in LEGACY_SCORE_SPECS.values()
        if spec.get("metric_id") not in (None, "")
    }
    docking_methods = {
        str(method).casefold() for method in DOCKING_PRIMARY_ENGINES
    } | {"vina"}
    for score in record.get("pose_scores") or []:
        if not isinstance(score, dict):
            continue
        if (
            str(score.get("stage") or "").casefold() == "docking"
            or str(score.get("metric_id") or "") in docking_metric_ids
            or str(score.get("method") or "").casefold()
            in docking_methods
        ):
            return True
    return False


def _selection_manifest_integer(
    value: object,
    *,
    field: str,
    pose_id: str,
    minimum: int,
) -> int:
    if isinstance(value, bool):
        raise ValueError(
            f"selection manifest {field} is boolean for {pose_id!r}"
        )
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str) and value.strip().isdigit():
        parsed = int(value.strip())
    else:
        raise ValueError(
            f"selection manifest has invalid {field} for {pose_id!r}"
        )
    if parsed < minimum:
        raise ValueError(
            f"selection manifest {field} must be >= {minimum} for {pose_id!r}"
        )
    return parsed


def _load_authoritative_pose_selection(
    eval_report: dict,
    *,
    eval_path: Path,
    selection_manifest_file: object = None,
) -> tuple[set[str] | None, dict[str, dict], dict]:
    """Load and validate the evaluator-owned post-QC pose selection.

    ``None`` denotes a legacy report with no authoritative selection contract;
    an empty set is a valid modern manifest in which no pose survived QC.
    """
    selection = eval_report.get("selection") or {}
    files = eval_report.get("files") or {}
    declared = (
        selection_manifest_file
        or selection.get("authoritative_artifact")
        or files.get("selection_manifest_csv")
    )
    embedded = selection.get("authoritative_selected_poses")
    rows: list[dict] | None = None
    source = None
    if declared not in (None, ""):
        path = Path(str(declared)).expanduser()
        if not path.is_absolute():
            path = eval_path.parent / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f"selection manifest not found: {path}")
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            required_columns = {
                "selection_rank",
                "pose_id",
                "mol_name",
                "source_pose_index",
                "post_qc_verdict",
            }
            missing_columns = sorted(
                required_columns - set(reader.fieldnames or ())
            )
            if missing_columns:
                raise ValueError(
                    "selection manifest is missing required columns: "
                    + ", ".join(missing_columns)
                )
            rows = [dict(row) for row in reader]
        source = str(path)
    elif isinstance(embedded, list):
        if any(not isinstance(row, dict) for row in embedded):
            raise ValueError(
                "embedded authoritative selection must contain only objects"
            )
        rows = [dict(row) for row in embedded]
        source = "evaluation_report.selection.authoritative_selected_poses"
    elif selection.get("post_qc_strategy") not in (None, ""):
        raise ValueError(
            "structure-evaluator report declares post-QC selection but does "
            "not provide its authoritative pose manifest"
        )
    else:
        return None, {}, {
            "method": "structure_evaluator_qc_fallback",
            "source": "legacy_evaluation_report_without_selection_contract",
        }

    raw_eligible_verdicts = selection.get("eligible_verdicts")
    if not isinstance(raw_eligible_verdicts, list):
        raise ValueError(
            "authoritative selection requires selection.eligible_verdicts"
        )
    if any(not isinstance(value, str) for value in raw_eligible_verdicts):
        raise ValueError("selection.eligible_verdicts must contain only strings")
    eligible_verdict_list = [value.strip() for value in raw_eligible_verdicts]
    if any(not value for value in eligible_verdict_list):
        raise ValueError("selection.eligible_verdicts contains an empty value")
    if len(set(eligible_verdict_list)) != len(eligible_verdict_list):
        raise ValueError("selection.eligible_verdicts contains duplicates")
    eligible_verdicts = set(eligible_verdict_list)

    pose_ids = [str(row.get("pose_id") or "").strip() for row in rows]
    if any(not pose_id for pose_id in pose_ids):
        raise ValueError("selection manifest contains an empty pose_id")
    if any(not pose_id.startswith("pose:") for pose_id in pose_ids):
        raise ValueError("selection manifest contains a non-canonical pose_id")
    if len(set(pose_ids)) != len(pose_ids):
        raise ValueError("selection manifest contains duplicate pose_id values")

    normalized_rows: list[dict] = []
    selection_ranks: list[int] = []
    for row, pose_id in zip(rows, pose_ids):
        selection_rank = _selection_manifest_integer(
            row.get("selection_rank"),
            field="selection_rank",
            pose_id=pose_id,
            minimum=1,
        )
        selection_ranks.append(selection_rank)
        normalized_rows.append({
            "selection_rank": selection_rank,
            "pose_id": pose_id,
            "mol_name": str(row.get("mol_name") or "").strip(),
            "source_pose_index": _selection_manifest_integer(
                row.get("source_pose_index"),
                field="source_pose_index",
                pose_id=pose_id,
                minimum=0,
            ),
            "post_qc_verdict": str(
                row.get("post_qc_verdict") or ""
            ).strip(),
        })
    if sorted(selection_ranks) != list(range(1, len(rows) + 1)):
        raise ValueError(
            "selection manifest ranks must be unique and contiguous from 1"
        )

    evaluation = eval_report.get("evaluation") or []
    evaluation_by_pose: dict[str, list[dict]] = {}
    for entry in evaluation:
        if not isinstance(entry, dict) or entry.get("pose_id") in (None, ""):
            continue
        evaluation_by_pose.setdefault(str(entry["pose_id"]), []).append(entry)
    for row, pose_id in zip(normalized_rows, pose_ids):
        matches = evaluation_by_pose.get(pose_id) or []
        if len(matches) != 1:
            raise ValueError(
                "selected pose_id must resolve to exactly one evaluation row: "
                f"{pose_id!r} resolved to {len(matches)}"
            )
        entry = matches[0]
        verdict = str(entry.get("verdict") or "")
        if verdict not in eligible_verdicts:
            raise ValueError(
                f"selection manifest names ineligible pose {pose_id!r} "
                f"with verdict {verdict!r}"
            )
        declared_verdict = str(row.get("post_qc_verdict") or "").strip()
        if not declared_verdict:
            raise ValueError(
                f"selection manifest omits post_qc_verdict for {pose_id!r}"
            )
        if declared_verdict != verdict:
            raise ValueError(
                f"selection verdict disagrees with evaluation for {pose_id!r}"
            )
        declared_name = str(row.get("mol_name") or "").strip()
        observed_name = str(entry.get("mol_name") or "").strip()
        if not declared_name:
            raise ValueError(
                f"selection manifest omits mol_name for {pose_id!r}"
            )
        if not observed_name:
            raise ValueError(
                f"evaluation row omits mol_name for selected pose {pose_id!r}"
            )
        if declared_name != observed_name:
            raise ValueError(
                f"selection molecule name disagrees for {pose_id!r}"
            )
        declared_index = row.get("source_pose_index")
        observed_index = entry.get(
            "source_pose_index", entry.get("pose_index")
        )
        if declared_index in (None, ""):
            raise ValueError(
                f"selection manifest omits source_pose_index for {pose_id!r}"
            )
        if observed_index in (None, ""):
            raise ValueError(
                f"evaluation row omits source_pose_index for {pose_id!r}"
            )
        if declared_index not in (None, ""):
            declared_index = _selection_manifest_integer(
                declared_index,
                field="source_pose_index",
                pose_id=pose_id,
                minimum=0,
            )
            observed_index = _selection_manifest_integer(
                observed_index,
                field="evaluation source_pose_index",
                pose_id=pose_id,
                minimum=0,
            )
            index_matches = declared_index == observed_index
            if not index_matches:
                raise ValueError(
                    f"selection source_pose_index disagrees for {pose_id!r}"
                )

    if isinstance(embedded, list) and declared not in (None, ""):
        embedded_rows = []
        try:
            for row in embedded:
                if not isinstance(row, dict):
                    raise ValueError
                embedded_rows.append({
                    "selection_rank": _selection_manifest_integer(
                        row.get("selection_rank"),
                        field="selection_rank",
                        pose_id=str(row.get("pose_id") or ""),
                        minimum=1,
                    ),
                    "pose_id": str(row.get("pose_id") or "").strip(),
                    "mol_name": str(row.get("mol_name") or "").strip(),
                    "source_pose_index": _selection_manifest_integer(
                        row.get("source_pose_index"),
                        field="source_pose_index",
                        pose_id=str(row.get("pose_id") or ""),
                        minimum=0,
                    ),
                    "post_qc_verdict": str(
                        row.get("post_qc_verdict") or ""
                    ).strip(),
                })
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "embedded authoritative selection contains invalid lineage"
            ) from exc
        manifest_rows = [
            {**row, "source_pose_index": int(row["source_pose_index"])}
            for row in normalized_rows
        ]
        if embedded_rows != manifest_rows:
            raise ValueError(
                "selection manifest disagrees with the selection embedded in "
                "the evaluation report"
            )

    declared_selected_count = selection.get("selected_poses")
    if declared_selected_count not in (None, ""):
        selected_count = _selection_manifest_integer(
            declared_selected_count,
            field="selection.selected_poses",
            pose_id="<selection>",
            minimum=0,
        )
        count_matches = selected_count == len(rows)
        if not count_matches:
            raise ValueError(
                "selection.selected_poses disagrees with the authoritative "
                "selection manifest"
            )

    digest = hashlib.sha256(
        "\n".join(pose_ids).encode("utf-8")
    ).hexdigest()
    records_digest = hashlib.sha256(json.dumps(
        normalized_rows,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")).hexdigest()
    records_by_pose_id = {
        row["pose_id"]: {
            **row,
            "source_pose_index": int(row["source_pose_index"]),
        }
        for row in normalized_rows
    }
    return set(pose_ids), records_by_pose_id, {
        "method": "structure_evaluator_selection_manifest",
        "source": source,
        "selected_pose_count": len(pose_ids),
        "selected_pose_ids_sha256": digest,
        "selected_pose_records_sha256": records_digest,
        "eligible_verdicts": sorted(eligible_verdicts),
    }


def _scientific_cohorts_from(*sources: object) -> list[dict]:
    """Collect exact cohort manifests without collapsing conflicting tuples."""
    cohorts: list[dict] = []
    seen: set[str] = set()
    for source in sources:
        if not isinstance(source, dict):
            continue
        payload = source.get("scientific_cohorts")
        if payload is None and isinstance(source.get("scientific_cohort"), dict):
            payload = [source["scientific_cohort"]]
        if not isinstance(payload, list):
            continue
        for item in payload:
            if not isinstance(item, dict):
                continue
            normalized = {
                field: deepcopy(item.get(field))
                for field in SCIENTIFIC_COHORT_FIELDS
                if item.get(field) not in (None, "")
            }
            if not normalized:
                continue
            valid, errors = _cohort_tuple_validation(normalized)
            normalized["cohort_manifest_valid"] = valid
            normalized["cohort_manifest_validation_errors"] = errors
            fingerprint = json.dumps(
                {
                    field: normalized.get(field)
                    for field in (
                        "cohort_id", *SCIENTIFIC_COHORT_COMPONENT_FIELDS,
                    )
                },
                sort_keys=True,
                default=str,
            )
            if fingerprint in seen:
                continue
            seen.add(fingerprint)
            cohorts.append(normalized)
    return cohorts


def _cohort_tuple_validation(item: dict) -> tuple[bool, list[str]]:
    """Validate that a declared cohort ID contains the same five-tuple."""
    errors = [
        f"missing:{field}"
        for field in ("cohort_id", *SCIENTIFIC_COHORT_COMPONENT_FIELDS)
        if item.get(field) in (None, "")
    ]
    if not errors:
        canonical = "|".join(
            f"{field}={item[field]}"
            for field in SCIENTIFIC_COHORT_COMPONENT_FIELDS
        )
        cohort_id = str(item["cohort_id"])
        declared_wrapper = (
            cohort_id.startswith("declared_cohort_id=")
            and cohort_id.endswith(f"|{canonical}")
        )
        if cohort_id != canonical and not declared_wrapper:
            errors.append("cohort_id_component_mismatch")
    return not errors, errors


def _cohort_tuple_is_valid(item: dict) -> bool:
    return _cohort_tuple_validation(item)[0]


def _cohort_identity(
    entry: dict,
    eval_report: dict,
    round_number: int,
    input_data: dict | None = None,
) -> dict[str, object]:
    """Build the sbdd scientific-cohort identity without inventing evidence."""
    parameters = eval_report.get("parameters") or {}
    docking = entry.get("docking_scores") or {}
    adapter = input_data or {}
    report_cohorts = _scientific_cohorts_from(eval_report)
    report_cohort = (
        report_cohorts[0]
        if len(report_cohorts) == 1 and _cohort_tuple_is_valid(report_cohorts[0])
        else {}
    )
    adapter_cohorts = _scientific_cohorts_from(adapter)
    adapter_cohort = (
        adapter_cohorts[0]
        if len(adapter_cohorts) == 1
        and _cohort_tuple_is_valid(adapter_cohorts[0])
        else {}
    )
    source_token = hashlib.sha256(
        json.dumps({
            "evaluation_report_file": adapter.get("evaluation_report_file"),
            "docking_scores_file": adapter.get("docking_scores_file"),
            "receptor_file": parameters.get("receptor_file"),
            "source_accession": parameters.get("source_accession"),
            "round_number": round_number,
            "ambiguous_entry_identity": (
                {
                    "pose_index": entry.get("pose_index"),
                    "pose_id": entry.get("pose_id"),
                    "chemical_identity_id": entry.get(
                        "chemical_identity_id"
                    ),
                    "mol_name": entry.get("mol_name"),
                }
                if entry.get("cohort_assignment_status") in {
                    "ambiguous", "conflict",
                }
                else None
            ),
        }, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()

    def first(*values: object, default: object) -> object:
        return next(
            (value for value in values if value not in (None, "")),
            default,
        )

    cohort_fields = (
        "protein_state_id", "pocket_id", "round_id",
        "protocol_profile_id", "metric_catalog_version",
    )
    components = {
        "protein_state_id": first(
            entry.get("protein_state_id"),
            docking.get("protein_state_id"),
            parameters.get("protein_state_id"),
            report_cohort.get("protein_state_id"),
            adapter.get("protein_state_id"),
            adapter_cohort.get("protein_state_id"),
            default=f"unknown:protein_state_id:{source_token}",
        ),
        "pocket_id": first(
            entry.get("pocket_id"),
            docking.get("pocket_id"),
            parameters.get("pocket_id"),
            report_cohort.get("pocket_id"),
            adapter.get("pocket_id"),
            adapter_cohort.get("pocket_id"),
            default=f"unknown:pocket_id:{source_token}",
        ),
        "round_id": first(
            entry.get("round_id"),
            parameters.get("round_id"),
            report_cohort.get("round_id"),
            adapter.get("round_id"),
            adapter_cohort.get("round_id"),
            default=round_number,
        ),
        "protocol_profile_id": first(
            entry.get("protocol_profile_id"),
            docking.get("protocol_profile_id"),
            docking.get("protocol_id"),
            parameters.get("protocol_profile_id"),
            report_cohort.get("protocol_profile_id"),
            adapter.get("protocol_profile_id"),
            adapter_cohort.get("protocol_profile_id"),
            default=f"unknown:protocol_profile_id:{source_token}",
        ),
        "metric_catalog_version": first(
            entry.get("metric_catalog_version"),
            parameters.get("metric_catalog_version"),
            report_cohort.get("metric_catalog_version"),
            adapter.get("metric_catalog_version"),
            adapter_cohort.get("metric_catalog_version"),
            default=METRIC_RULES_VERSION,
        ),
    }
    explicit = first(
        entry.get("declared_cohort_id"),
        docking.get("declared_cohort_id"),
        entry.get("cohort_id"),
        docking.get("cohort_id"),
        parameters.get("cohort_id"),
        report_cohort.get("cohort_id"),
        adapter.get("cohort_id"),
        adapter_cohort.get("cohort_id"),
        default=None,
    )
    canonical_id = "|".join(
        f"{key}={components[key]}" for key in cohort_fields
    )
    provenance_sources = first(
        entry.get("cohort_provenance_sources"),
        docking.get("cohort_provenance_sources"),
        parameters.get("cohort_provenance_sources"),
        report_cohort.get("cohort_provenance_sources"),
        adapter.get("cohort_provenance_sources"),
        adapter_cohort.get("cohort_provenance_sources"),
        default={},
    )
    provenance_sources = (
        dict(provenance_sources)
        if isinstance(provenance_sources, dict) else {}
    )
    missing_fields = [
        field for field in cohort_fields
        if str(components[field]).startswith(("unknown:", "pending:"))
        or "|state_source:" in str(components[field])
    ]
    derived_fields = first(
        entry.get("cohort_identity_derived_fields"),
        docking.get("cohort_identity_derived_fields"),
        parameters.get("cohort_identity_derived_fields"),
        report_cohort.get("cohort_identity_derived_fields"),
        adapter.get("cohort_identity_derived_fields"),
        adapter_cohort.get("cohort_identity_derived_fields"),
        default=[],
    )
    components.update({
        "cohort_identity_complete": not missing_fields,
        "cohort_identity_missing_fields": missing_fields,
        "cohort_identity_derived_fields": (
            list(derived_fields) if isinstance(derived_fields, list) else []
        ),
        "cohort_provenance_sources": provenance_sources,
    })
    if explicit is None or str(explicit) == canonical_id:
        components["cohort_id"] = canonical_id
        return components

    declared_id = str(explicit)
    prefix = "declared_cohort_id="
    suffix = f"|{canonical_id}"
    if declared_id.startswith(prefix) and declared_id.endswith(suffix):
        declared_id = declared_id[len(prefix):-len(suffix)]
    components["declared_cohort_id"] = declared_id
    # A human-facing declared ID is provenance, not proof of comparability.
    # Include the five scientific components in the actual percentile key so
    # accidentally reusing one ID across proteins/pockets/profiles cannot mix
    # normalization pools.
    components["cohort_id"] = (
        f"declared_cohort_id={declared_id}|{canonical_id}"
    )
    return components


def _score_ranking_cohorts(
    records: list[dict],
    *,
    normalization_method: str,
    metric_rules: tuple[dict, ...] = POSE_PROPERTY_METRIC_RULES,
    metric_weight_overrides: dict | None = None,
    category_weight_overrides: dict | None = None,
    active_steps: tuple[str, ...] | None = None,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> tuple[list[dict], dict]:
    """Score each declared scientific cohort; keep decoys in separate pools."""
    mode = resolve_screening_mode(screening_mode)
    scored = [dict(record) for record in records]
    grouped: dict[str, list[int]] = {}
    for index, record in enumerate(records):
        cohort_id = str(record.get("cohort_id") or "cohort:unspecified")
        grouped.setdefault(cohort_id, []).append(index)

    cohorts: dict[str, dict] = {}
    for cohort_id, cohort_indices in grouped.items():
        candidate_indices = [
            index for index in cohort_indices
            if not _is_decoy_record(records[index])
        ]
        decoy_indices = [
            index for index in cohort_indices
            if _is_decoy_record(records[index])
        ]
        cohort_metadata: dict[str, dict] = {}
        if candidate_indices:
            candidate_scored, candidate_metadata = score_records(
                (records[index] for index in candidate_indices),
                rules=metric_rules,
                normalization_method=normalization_method,
                metric_weight_overrides=metric_weight_overrides,
                category_weight_overrides=category_weight_overrides,
                active_steps=active_steps,
                screening_mode=mode,
            )
            candidate_metadata["normalization_scope"] = (
                "nomination_candidates_in_scientific_cohort"
            )
            candidate_metadata["cohort_id"] = cohort_id
            for index, record in zip(candidate_indices, candidate_scored):
                record["normalization_pool"] = "nomination_candidates"
                scored[index] = record
            cohort_metadata["nomination_candidates"] = candidate_metadata
        if decoy_indices:
            decoy_scored, decoy_metadata = score_records(
                (records[index] for index in decoy_indices),
                rules=metric_rules,
                normalization_method=normalization_method,
                metric_weight_overrides=metric_weight_overrides,
                category_weight_overrides=category_weight_overrides,
                active_steps=active_steps,
                screening_mode=mode,
            )
            decoy_metadata["normalization_scope"] = (
                "benchmark_decoys_in_scientific_cohort"
            )
            decoy_metadata["cohort_id"] = cohort_id
            for index, record in zip(decoy_indices, decoy_scored):
                record["normalization_pool"] = "benchmark_decoys"
                scored[index] = record
            cohort_metadata["benchmark_decoys"] = decoy_metadata
        cohorts[cohort_id] = cohort_metadata

    return scored, {
        "schema_version": SCORING_SCHEMA_VERSION,
        "metric_rules_version": METRIC_RULES_VERSION,
        "normalization_method": normalization_method,
        # Level 1 marks each metric pass/flag/fail; level 2 folds those marks per
        # parent category under this mode and only then decides the molecule.
        "gate_contract": "two_level_label_then_category_fold",
        "screening_mode": dict(mode),
        "metric_weight_overrides": dict(metric_weight_overrides or {}),
        "category_weight_overrides": dict(category_weight_overrides or {}),
        "cohorts": cohorts,
        "cohort_count": len(cohorts),
        "active_steps": list(active_steps or ()),
        "cohort_fields": [
            "protein_state_id", "pocket_id", "round_id",
            "protocol_profile_id", "metric_catalog_version",
        ],
        "decoys_isolated_from_nomination_normalization": True,
    }


def _effective_metric_rules(
    eval_report: dict,
    property_report: dict | None = None,
) -> tuple[dict, ...]:
    """Combine evaluator-owned pose policy with Step 5's public property policy.

    Property reports intentionally publish rules without extractors.  Only known
    property-filter metric IDs may override policy fields on the local executable
    templates; pose rules and extraction paths always remain coordinator-owned.
    """
    pose_rules = list(deepcopy(POSE_EVALUATION_METRIC_RULES))
    thresholds = (
        ((eval_report.get("parameters") or {}).get("posecheck_thresholds"))
        or {}
    )
    overrides = {
        "posecheck.strain_energy_kcal": (
            thresholds.get("strain_warn_kcal"),
            thresholds.get("strain_fail_kcal"),
        ),
        "posecheck.steric_clash_count": (
            thresholds.get("clash_warn_count"),
            thresholds.get("clash_fail_count"),
        ),
    }
    for rule in pose_rules:
        soft_upper, hard_upper = overrides.get(
            rule.get("metric_id"),
            (None, None),
        )
        if soft_upper is not None:
            rule["soft_upper"] = soft_upper
        if hard_upper is not None:
            rule["hard_upper"] = hard_upper

    property_rule_by_id = {
        str(rule.get("metric_id")): rule
        for rule in deepcopy(PROPERTY_FILTER_METRIC_RULES)
    }
    property_rules: list[dict] = []
    public_rules = (
        property_report.get("metric_rules")
        if isinstance(property_report, dict)
        else None
    )
    seen_public_ids: set[str] = set()
    policy_fields = (
        "hard_lower",
        "soft_lower",
        "soft_upper",
        "hard_upper",
        "weight",
        "signed_weight",
        "direction",
        "rank_transform",
        "category_code",
        "step_id",
        "category_weight",
        "catalog_abs_weight",
        "plan_tool",
        "molecule_effect",
        "family_id",
        "family_weight_cap",
        "requires_property_tool",
        "covalent_review",
        "description",
    )
    if isinstance(public_rules, list):
        for public_rule in public_rules:
            if not isinstance(public_rule, dict):
                continue
            metric_id = str(public_rule.get("metric_id") or "")
            if not metric_id or metric_id in seen_public_ids:
                continue
            seen_public_ids.add(metric_id)
            template = property_rule_by_id.get(metric_id)
            if template is None:
                continue
            stage = public_rule.get("stage")
            if stage not in (None, "", "property_filter"):
                continue
            candidate = deepcopy(template)
            for field in policy_fields:
                if field in public_rule:
                    candidate[field] = deepcopy(public_rule[field])
            try:
                validate_metric_rules((candidate,))
            except (TypeError, ValueError):
                # A stale or malformed report cannot corrupt the executable
                # policy. Its metric remains represented by Step 5's already
                # aggregated classification, but is not re-scored here.
                continue
            property_rules.append(candidate)

    rules = tuple((*pose_rules, *property_rules))
    validate_metric_rules(rules)
    return rules


def _combine_classifications(*values: object) -> str:
    """Reduce applicable gate labels with reject dominance.

    Returns ``not_evaluated`` rather than ``reject`` when no value is a known
    label: an unreadable record is missing evidence, and the executable policy forbids
    turning missing evidence into a molecular rejection.
    """
    applicable = [
        str(value)
        for value in values
        if str(value) in _CLASS_SEVERITY
    ]
    if not applicable:
        return "not_evaluated"
    return max(applicable, key=_CLASS_SEVERITY.__getitem__)


def _apply_coverage_consequences(record: dict) -> None:
    """Route a coverage shortfall to the rank, never to the molecule label.

    The executable policy names two bands and one consequence each.  Between 0.10 and 0.50
    the rank is still published but only as ``provisional``, which ``rank_status``
    already carries; all this adds is a visible nomination annotation so a thin
    round is not silently nominated as if it were fully evidenced.  Below 0.10 --
    and likewise when an active step produced no evidence at all -- the spec asks
    to backtrack and refill the missing categories, so the record joins the same
    ``ranking_blocked``/``ranking_blocked_by`` path protocol failures use.

    Either way the keep/modify/reject label is untouched: a category this round
    never planned, or a tool that was switched off, is not a finding about the
    chemistry.
    """
    if str(record.get("coverage_gate") or "") != "modify":
        return
    rank_status = str(record.get("rank_status") or "")
    missing_steps = list(record.get("missing_active_steps") or [])
    blockers = [f"coverage:{step}" for step in missing_steps]
    missing_gate_ids = list(
        record.get("missing_planned_zero_weight_gate_ids") or []
    )
    blockers.extend(
        f"planned_gate_evidence:{metric_id}" for metric_id in missing_gate_ids
    )
    if rank_status == "insufficient_coverage":
        blockers.append("coverage:below_provisional_threshold")
    if blockers:
        record["ranking_blocked"] = True
        record["rank"] = record["overall_rank"] = record["weighted_score_rank"] = None
        record["ranking_blocked_by"] = list(
            dict.fromkeys(list(record.get("ranking_blocked_by") or []) + blockers)
        )
        if not str(
            record.get("nomination_status") or "eligible"
        ).startswith("blocked"):
            record["nomination_status"] = (
                "blocked_by_planned_gate_evidence"
                if missing_gate_ids
                else "blocked_by_metric_coverage"
            )
        record["modification_hint"] = (
            (record.get("modification_hint") or "")
            + "; Backfill the missing metric categories before nomination: "
            + ", ".join(blockers)
        ).strip("; ")
        return
    if str(record.get("nomination_status") or "eligible") == "eligible":
        # Deliberately not a ``blocked_`` status: a provisional rank is still a
        # publishable rank, so this annotates rather than vetoes.
        record["nomination_status"] = "coverage_review"


def _chemical_identity_id(entry: dict) -> object:
    """Read either the structured or flat chemical-identity contract safely."""
    identity = entry.get("chemical_identity")
    nested_id = identity.get("identity_id") if isinstance(identity, dict) else None
    return nested_id or entry.get("chemical_identity_id")


def _chemical_id(entry: dict) -> object:
    """Read the source-independent chemical ID from either identity shape."""
    identity = entry.get("chemical_identity")
    nested_id = identity.get("chemical_id") if isinstance(identity, dict) else None
    return nested_id or entry.get("chemical_id")


def _add_property_lookup(index: dict, key: object, result: dict) -> None:
    """Index a property result without silently resolving duplicate keys."""
    if key in (None, ""):
        return
    bucket = index.setdefault(key, [])
    if not any(existing is result for existing in bucket):
        bucket.append(result)


def _unique_property_result(
    index: dict,
    key: object,
    *,
    cohort: dict,
    allow_legacy_cohort_inference: bool,
    allow_exact_pose_without_cohort: bool = False,
) -> tuple[dict | None, str]:
    """Return a cohort-compatible unique row and its match state."""
    if key in (None, ""):
        return None, "not_found"
    matches = index.get(key) or []
    cohort_id = cohort.get("cohort_id")
    wanted = str(cohort_id) if cohort_id not in (None, "") else None
    exact = [
        item for item in matches
        if (
            str(item.get("cohort_id"))
            if item.get("cohort_id") not in (None, "") else None
        ) == wanted and wanted is not None
    ]
    compatible = [
        item for item in exact
        if _cohort_rows_are_compatible(cohort, item)
    ]
    if len(compatible) == 1 and len(exact) == 1:
        return compatible[0], "verified"
    if exact and not compatible:
        return None, "cohort_mismatch"
    if len(exact) > 1:
        return None, "ambiguous"
    missing_cohort = [
        item for item in matches if item.get("cohort_id") in (None, "")
    ]
    if (
        wanted is None
        and allow_exact_pose_without_cohort
        and len(missing_cohort) == 1
        and len(matches) == 1
    ):
        # A canonical pose_id already binds both chemical identity and exact
        # coordinates.  This keeps legacy standalone reports interoperable
        # without weakening any molecule/name based join.
        return missing_cohort[0], "exact_pose_without_cohort"
    if (
        wanted is not None
        and allow_legacy_cohort_inference
        and len(missing_cohort) == 1
        and len(matches) == 1
    ):
        return missing_cohort[0], "inferred"
    if missing_cohort:
        return None, "cohort_unscoped"
    if matches:
        return None, "cohort_mismatch"
    return None, "not_found"


def _cohort_rows_are_compatible(left: dict, right: dict) -> bool:
    """Require both modern rows to publish one internally consistent tuple."""
    if not _cohort_tuple_is_valid(left) or not _cohort_tuple_is_valid(right):
        return False
    return all(
        str(left.get(field)) == str(right.get(field))
        for field in SCIENTIFIC_COHORT_COMPONENT_FIELDS
    )


def _observed_property_cohort_ids(index: dict, key: object) -> list[str]:
    if key in (None, ""):
        return []
    values = {
        str(item.get("cohort_id"))
        if item.get("cohort_id") not in (None, "")
        else "<missing>"
        for item in (index.get(key) or [])
        if isinstance(item, dict)
    }
    return sorted(values)


def _match_property_result(
    entry: dict,
    *,
    by_chemical_identity: dict,
    by_chemical_id: dict,
    by_pose_id: dict,
    by_name: dict,
    by_smiles: dict,
    allow_legacy_cohort_inference: bool = False,
) -> tuple[dict | None, dict]:
    """Join property evidence only when it names this exact canonical pose."""
    cohort_id = entry.get("cohort_id")
    pose_id = entry.get("pose_id")
    if pose_id in (None, ""):
        # Molecule-level IDs cannot authorize a pose join, but they can still
        # explain why the available row was rejected.  Preserve that cohort
        # diagnosis instead of flattening every case into "pose ID missing".
        for field, key, index in (
            (
                "chemical_identity_id",
                _chemical_identity_id(entry),
                by_chemical_identity,
            ),
            ("chemical_id", _chemical_id(entry), by_chemical_id),
        ):
            _candidate, diagnostic_state = _unique_property_result(
                index,
                key,
                cohort=entry,
                allow_legacy_cohort_inference=(
                    allow_legacy_cohort_inference
                ),
            )
            if diagnostic_state == "cohort_mismatch":
                return None, {
                    "status": "cohort_mismatch",
                    "matched_by": field,
                    "reason": "property_cohort_mismatch",
                    "expected_cohort_id": cohort_id,
                    "observed_cohort_ids": _observed_property_cohort_ids(
                        index, key,
                    ),
                }
        return None, {
            "status": "not_found",
            "matched_by": None,
            "reason": "selected_pose_id_missing",
            "checked": ["pose_id"],
        }
    result, match_state = _unique_property_result(
        by_pose_id,
        pose_id,
        cohort=entry,
        allow_legacy_cohort_inference=allow_legacy_cohort_inference,
        allow_exact_pose_without_cohort=True,
    )
    if result is not None:
        return result, {
            "status": "matched",
            "matched_by": "pose_id",
            "exact_pose_id": str(pose_id),
            "cohort_verified": match_state == "verified",
            "cohort_inferred": match_state == "inferred",
        }
    if match_state == "cohort_mismatch":
        return None, {
            "status": "cohort_mismatch",
            "matched_by": "pose_id",
            "reason": "property_cohort_mismatch",
            "expected_cohort_id": cohort_id,
            "observed_cohort_ids": _observed_property_cohort_ids(
                by_pose_id, pose_id,
            ),
        }
    if match_state == "cohort_unscoped":
        return None, {
            "status": "cohort_unscoped",
            "matched_by": "pose_id",
            "reason": "property_cohort_unscoped",
            "expected_cohort_id": cohort_id,
        }
    if match_state == "ambiguous":
        return None, {
            "status": "ambiguous",
            "matched_by": "pose_id",
            "reason": "duplicate_property_pose_id",
        }

    # Molecule-level identifiers are diagnostic only: they may reveal that a
    # property row exists for another pose, but can never authorize the join.
    related_rows: list[dict] = []
    for key, index in (
        (_chemical_identity_id(entry), by_chemical_identity),
        (_chemical_id(entry), by_chemical_id),
    ):
        for candidate in index.get(key, []) if key not in (None, "") else []:
            if not any(existing is candidate for existing in related_rows):
                related_rows.append(candidate)
    if related_rows:
        observed_pose_ids = sorted({
            str(candidate.get("pose_id"))
            for candidate in related_rows
            if candidate.get("pose_id") not in (None, "")
        })
        return None, {
            "status": "pose_mismatch",
            "matched_by": None,
            "reason": "property_pose_id_mismatch",
            "expected_pose_id": str(pose_id),
            "observed_pose_ids": observed_pose_ids,
            "checked": ["pose_id"],
        }
    return None, {
        "status": "not_found",
        "matched_by": None,
        "reason": "no_matching_pose_id",
        "checked": ["pose_id"],
    }


def _extract_failed_scaffolds(eval_results: list[dict], classifications: list[str]) -> list[str]:
    """Extract Murcko scaffolds from rejected/failed molecules for exclusion in next round."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:
        return []

    scaffolds = set()
    for entry, cls in zip(eval_results, classifications):
        if cls != "reject":
            continue
        smiles = entry.get("smiles", "")
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            scaffolds.add(Chem.MolToSmiles(core))
        except Exception:
            continue
    return sorted(scaffolds)


def _is_decoy_record(record: dict) -> bool:
    return benchmark_label(record) == "DEC"


def _candidate_records(records: list[dict]) -> list[dict]:
    """Records eligible for nomination in this round.

    DEC_* molecules are retrospective benchmark decoys. They remain in the
    coordinator report and screening metrics, but they are not hit candidates.
    """
    return [dict(record) for record in records if not _is_decoy_record(record)]


def _eligible_keep(record: dict) -> bool:
    """Count only trusted, unblocked keeps toward nomination convergence."""
    if record.get("classification") != "keep":
        return False
    if record.get("ranking_blocked") is True:
        return False
    status = str(record.get("nomination_status") or "eligible").casefold()
    return not status.startswith("blocked")


def _gate_ok(gate: dict) -> bool:
    """Return whether a protocol gate permits nomination."""
    return _quality_gate_status(gate) in {"keep", "not_applicable"}


def _tri_state_from_passed(passed) -> str:
    """Map a binary/absent ``passed`` flag onto the three-tier gate vocabulary.

    ``True`` is a pass and ``False`` is a real protocol failure, but anything
    else (``None``, absent, non-boolean) means the check produced no verdict.
    Absence of evidence must not invalidate a run whose only defect is a report
    the coordinator was never given. It becomes ``modify``, a review finding
    that the Agent can resolve by obtaining the missing evidence.
    """
    if passed is True:
        return "keep"
    if passed is False:
        return "reject"
    return "modify"


def _quality_gate_status(gate: dict) -> str:
    """Normalize legacy binary and current three-tier quality-gate records."""
    if gate.get("not_applicable") is True:
        return "not_applicable"
    explicit = str(gate.get("gate") or "").casefold()
    if explicit == "not_applicable":
        return "not_applicable"
    if explicit in {"keep", "modify", "reject"}:
        return explicit
    return _tri_state_from_passed(gate.get("passed"))


def _aggregate_quality_gate_status(gates: dict[str, dict]) -> str:
    """Aggregate protocol gates with reject dominance and NA neutrality."""
    statuses = [_quality_gate_status(gate) for gate in gates.values()]
    applicable = [status for status in statuses if status != "not_applicable"]
    if "reject" in applicable:
        return "reject"
    if "modify" in applicable:
        return "modify"
    if applicable:
        return "keep"
    return "not_applicable"


def _screening_benchmark_gate(
    benchmark: dict,
    *,
    min_actives: int,
    min_decoys: int,
    auroc_soft: float,
    auroc_hard: float,
    ef1_soft: float,
    ef1_hard: float,
) -> dict:
    """Apply protocol-level AUROC/EF1 gates to a reliable labeled benchmark.

    This is a protocol gate, not a molecular one: it judges whether the docking
    protocol can discriminate at all. The executable policy keeps protocol conclusions on the
    single-metric level -- one hard breach invalidates the run and triggers a
    rollback -- so the two-level category fold, which exists to stop one bad
    metric from discarding a molecule, deliberately does not apply here.
    """
    rules = {
        rule["metric_id"]: dict(rule)
        for rule in DOCKING_PROTOCOL_METRIC_RULES
    }
    auroc_rule = rules["screening_benchmark.AUROC"]
    auroc_rule.update({"soft_lower": auroc_soft, "hard_lower": auroc_hard})
    ef1_rule = rules["screening_benchmark.EF1pct"]
    ef1_rule.update({"soft_lower": ef1_soft, "hard_lower": ef1_hard})
    assessments = []
    for rule, value, hard_lower, soft_lower in (
        (auroc_rule, benchmark.get("AUROC"), auroc_hard, auroc_soft),
        (ef1_rule, benchmark.get("EF1%"), ef1_hard, ef1_soft),
    ):
        label = metric_label(value, rule)
        assessments.append({
            "metric_id": rule["metric_id"],
            "value": value,
            "label": label,
            # Legacy spelling kept for report consumers that still read `gate`.
            "gate": _LABEL_TO_LEGACY_GATE.get(label, label),
            "hard_lower": hard_lower,
            "soft_lower": soft_lower,
        })
    enough_samples = (
        int(benchmark.get("num_actives", 0)) >= min_actives
        and int(benchmark.get("num_decoys", 0)) >= min_decoys
    )
    if not enough_samples:
        return {
            "gate": "modify",
            "reason": "insufficient labeled actives/decoys; extend benchmark",
            "metric_assessments": assessments,
        }
    labels = [assessment["label"] for assessment in assessments]
    aggregate = (
        "reject" if "fail" in labels
        else "modify" if "flag" in labels
        else "modify" if "not_applicable" in labels
        else "keep"
    )
    return {
        "gate": aggregate,
        "reason": "three-tier AUROC and EF1% protocol gate",
        "metric_assessments": assessments,
    }


def _is_not_applicable(details: dict) -> bool:
    """Recognize the common tri-state gate representation."""
    return bool(
        details.get("not_applicable") is True
        or details.get("applicable") is False
    )


def _quality_gates(docking_report_file: str | None, pocket_report_file: str | None) -> dict:
    gates = {}
    if docking_report_file:
        path = Path(docking_report_file)
        if path.exists():
            verified = load_committed_docking_artifacts(path)
            redocking = verified.get("redocking_validation") or {}
            evidence_status = verified.get("redocking_evidence_status")
            if not verified.get("available") or evidence_status == "invalid":
                gates["native_redocking"] = {
                    "passed": None,
                    "not_applicable": False,
                    "gate": "modify",
                    "reason": (
                        verified.get("reason")
                        or redocking.get("reason")
                        or "redocking evidence is not committed and verified"
                    ),
                    "protocol_blocker": True,
                    "details": redocking,
                }
            elif evidence_status == "verified":
                gates["native_redocking"] = {
                    "passed": redocking.get("passed"),
                    "not_applicable": False,
                    "gate": (
                        redocking.get("gate")
                        or _tri_state_from_passed(redocking.get("passed"))
                    ),
                    "details": redocking,
                }
            else:
                # A valid Step 3 contract explicitly records that native
                # redocking was disabled or could not produce a measurement.
                # It is N/A rather than a molecular failure; callers that
                # require an applicable quality gate can still block on it.
                gates["native_redocking"] = {
                    "passed": None,
                    "not_applicable": True,
                    "gate": "not_applicable",
                    "reason": redocking.get("reason"),
                    "details": redocking,
                }
        else:
            gates["native_redocking"] = {
                "passed": None,
                "gate": "modify",
                "reason": f"missing report: {path}",
            }
    else:
        gates["native_redocking"] = {
            "passed": None,
            "gate": "modify",
            "reason": "docking_report_file not provided",
        }
    if pocket_report_file:
        path = Path(pocket_report_file)
        if path.exists():
            report = json.loads(path.read_text())
            validation = report.get("validation") or {}
            blueprint = report.get("interaction_blueprint") or {}
            blueprint_validation = blueprint.get("validation") or {}
            pocket_not_applicable = _is_not_applicable(validation)
            blueprint_not_applicable = (
                _is_not_applicable(blueprint_validation)
                or _is_not_applicable(blueprint)
            )
            blueprint_valid = blueprint_validation.get("blueprint_valid")
            if blueprint_valid is None:
                # Legacy files conflate "valid" and "has a hard constraint".
                # Absence of a required constraint is not a malformed Blueprint.
                blueprint_valid = not bool(blueprint_validation.get("invalid"))
            required_hard_gate_applicable = blueprint_validation.get(
                "required_hard_gate_applicable"
            )
            if required_hard_gate_applicable is None:
                required_hard_gate_applicable = blueprint_validation.get(
                    "hard_gate_applicable"
                )
            if required_hard_gate_applicable is None:
                required_hard_gate_applicable = blueprint_validation.get(
                    "valid_for_hard_gating"
                )
            if required_hard_gate_applicable is None:
                required_hard_gate_applicable = bool(blueprint.get("required"))
            desired_metric_applicable = blueprint_validation.get(
                "desired_metric_applicable"
            )
            if desired_metric_applicable is None:
                desired_metric_applicable = bool(
                    not blueprint_not_applicable
                    and blueprint_valid
                    and blueprint.get("desired")
                )
            molecule_hard_gate_applicable = blueprint_validation.get(
                "molecule_hard_gate_applicable"
            )
            if molecule_hard_gate_applicable is None:
                molecule_hard_gate_applicable = bool(
                    required_hard_gate_applicable or desired_metric_applicable
                )
            gates["pocket_validation"] = {
                "passed": None if pocket_not_applicable else validation.get("passed"),
                "not_applicable": pocket_not_applicable,
                "gate": (
                    "not_applicable" if pocket_not_applicable
                    else _tri_state_from_passed(validation.get("passed"))
                ),
                "details": validation,
            }
            gates["blueprint_validation"] = {
                "passed": (
                    None if blueprint_not_applicable
                    else bool(blueprint_valid)
                ),
                "not_applicable": blueprint_not_applicable,
                "gate": (
                    "not_applicable" if blueprint_not_applicable
                    else "keep" if bool(blueprint_valid)
                    else "reject"
                ),
                "blueprint_valid": None if blueprint_not_applicable else bool(blueprint_valid),
                # Legacy alias: this field is scoped only to the curated
                # Required exact-match verdict.
                "hard_gate_applicable": (
                    False if blueprint_not_applicable
                    else bool(required_hard_gate_applicable)
                ),
                "hard_gate_status": (
                    "not_applicable"
                    if blueprint_not_applicable or not required_hard_gate_applicable
                    else "declared"
                ),
                "required_hard_gate_applicable": (
                    False if blueprint_not_applicable
                    else bool(required_hard_gate_applicable)
                ),
                "required_hard_gate_status": (
                    "not_applicable"
                    if blueprint_not_applicable or not required_hard_gate_applicable
                    else "declared"
                ),
                "desired_metric_applicable": (
                    False if blueprint_not_applicable
                    else bool(desired_metric_applicable)
                ),
                "molecule_hard_gate_applicable": (
                    False if blueprint_not_applicable
                    else bool(molecule_hard_gate_applicable)
                ),
                "details": blueprint_validation,
            }
        else:
            gates["pocket_validation"] = {
                "passed": None, "gate": "modify",
                "reason": f"missing report: {path}",
            }
            gates["blueprint_validation"] = {
                "passed": None, "gate": "modify",
                "reason": f"missing report: {path}",
            }
    else:
        gates["pocket_validation"] = {
            "passed": None, "gate": "modify",
            "reason": "pocket_analysis_report_file not provided",
        }
        gates["blueprint_validation"] = {
            "passed": None, "gate": "modify",
            "reason": "pocket_analysis_report_file not provided",
        }
    return gates


def _resolved_report_reference(
    value: Any,
    report_file: str | None,
) -> Any:
    """Resolve one producer path before making it publication-relative.

    Older reports may hold paths relative either to their own directory or to
    the producer's launch directory.  Prefer the report-local interpretation
    when it exists, then the current working directory; if neither exists, keep
    the report-local absolute candidate so the ambiguity is visible and stable.
    """
    if not isinstance(value, (str, Path)) or not value:
        return value
    candidate = Path(value).expanduser()
    if candidate.is_absolute():
        return str(candidate.resolve())
    report_candidate = (
        (Path(report_file).resolve().parent / candidate).resolve()
        if report_file else None
    )
    cwd_candidate = candidate.resolve()
    if report_candidate is not None and report_candidate.exists():
        return str(report_candidate)
    if cwd_candidate.exists():
        return str(cwd_candidate)
    return str(report_candidate or cwd_candidate)


def _publish_output_layout(
    output_dir: Path,
    report: dict,
    *,
    eval_file: str | None,
    prop_file: str | None,
    docking_file: str | None,
    eval_report: dict | None = None,
    prop_report: dict | None = None,
    docking_report: dict | None = None,
    publication_root: Path | None = None,
    publication_mode: str = "standalone",
) -> dict[str, str]:
    """Publish the round's reader tree without mixing it into the working tree.

    The tree is the reader-facing view of this round: eight workflow directories
    from inputs through final report, with scored stages drilling down from step
    to category to the tool that produced the value.  The legacy step
    directories keep the raw output; every table in ``published/`` points back
    to it with a relative reference.

    The raw upstream reports come along because they are the only witnesses to a
    tool that failed for the whole cohort: with no row holding a value, the tables
    cannot tell that gap apart from a metric this build never planned.
    """
    from failure_registry import collect_failures
    from publish_layout import publish_run, publish_run_transactionally
    from output_layout import PUBLISHED_DIR

    eval_file = (
        str(Path(eval_file).expanduser().resolve()) if eval_file else None
    )
    prop_file = (
        str(Path(prop_file).expanduser().resolve()) if prop_file else None
    )
    docking_file = (
        str(Path(docking_file).expanduser().resolve()) if docking_file else None
    )
    round_root = output_dir.parent.resolve()
    publish_root = (
        Path(publication_root).resolve()
        if publication_root is not None else round_root / PUBLISHED_DIR
    )

    records = []
    for source_record in report.get("molecules") or []:
        if not isinstance(source_record, dict):
            continue
        record = dict(source_record)
        if record.get("pose_file"):
            record["pose_file"] = _resolved_report_reference(
                record["pose_file"], eval_file,
            )
        records.append(record)
    rules = report.get("metric_rules") or POSE_PROPERTY_METRIC_RULES
    first = records[0] if records else {}
    parameters = (eval_report or {}).get("parameters") or {}
    identity = {
        "round_id": first.get("round_id"),
        "cohort_id": first.get("cohort_id"),
        "profile_id": first.get("protocol_profile_id"),
        "protein_state_id": first.get("protein_state_id"),
        "pocket_id": first.get("pocket_id"),
        "engine": first.get("docking_engine"),
        # The receptor these poses were docked against.  Every row shares it, so
        # it is run-level identity rather than a per-record field.
        "source_accession": parameters.get("source_accession"),
        "source_type": parameters.get("source_type"),
    }
    # References, not copies: what every candidate in the cohort was scored
    # against.  Without this the tree gave no way to tell which receptor a
    # published ranking belongs to.
    shared_assets = {
        key: value
        for key, value in (
            (
                "receptor_file",
                _resolved_report_reference(
                    parameters.get("receptor_file"), eval_file,
                ),
            ),
            (
                "reference_ligand_file",
                _resolved_report_reference(
                    parameters.get("reference_ligand_file"), eval_file,
                ),
            ),
            (
                "docking_scores_file",
                _resolved_report_reference(
                    parameters.get("docking_scores_file"), eval_file,
                ),
            ),
        )
        if value
    }
    step_sources: dict[str, list[str]] = {}
    if prop_file:
        step_sources["step4"] = [prop_file]
    if eval_file:
        step_sources["step5"] = [eval_file]
        step_sources["step6"] = [eval_file]
    if docking_file:
        step_sources.setdefault("step5", []).insert(0, docking_file)
    publisher = (
        publish_run
        if publication_root is not None else publish_run_transactionally
    )
    written = publisher(
        publish_root,
        records=records,
        metadata=report.get("metric_scoring") or {},
        rules=rules,
        identity={key: value for key, value in identity.items() if value},
        profile={
            "metric_rules_version": report.get("metric_rules_version"),
            "row_scope": "all_evaluated_records",
        },
        step_sources=step_sources,
        shared_assets=shared_assets,
        inputs={
            "evaluation_report_file": eval_file,
            "property_filter_report_file": prop_file,
            "docking_report_file": docking_file,
        },
        failures=collect_failures(
            report,
            eval_report=eval_report,
            prop_report=prop_report,
            docking_report=docking_report,
            eval_file=eval_file,
            prop_file=prop_file,
            docking_file=docking_file,
            rules=rules,
        ),
        session_id=round_root.name,
        working_root=round_root,
        publication_mode=publication_mode,
    )
    return {f"layout_{role}": path for role, path in written.items()}


def _property_value(prop: dict, *names: str) -> object:
    properties = prop.get("properties") or {}
    for name in names:
        if prop.get(name) not in (None, ""):
            return prop.get(name)
        if properties.get(name) not in (None, ""):
            return properties.get(name)
    return ""


def _write_ranking_csv(path: Path, records: list[dict]) -> None:
    """Write the stable compact ranking view retained for API compatibility."""
    fields = [
        "rank", "mol_name", "smiles", "selected_pose_index", "selection_rank",
        "source_pose_index", "pose_id", "chemical_identity_id", "classification",
        "docking_artifact_set_id",
        "cohort_id", "protein_state_id", "pocket_id", "round_id",
        "protocol_profile_id", "metric_catalog_version", "normalization_pool",
        "rank_partition_id", "rank_scope", "ranking_schema",
        "scoring_schema_version", "intrinsic_classification",
        "nomination_status", "eval_verdict", "prop_verdict", "metric_gate",
        "coverage_status", "rank_status", "overall_score", "overall_rank",
        "step4_score", "step4_rank", "step5_score", "step5_rank",
        "step6_score", "step6_rank", "category_weighted_coverage",
        "weighted_total_score", "weighted_score_max", "weighted_score_fraction",
        "weighted_score_rank", "score_coverage", "ranking_blocked",
        "ranking_blocked_by", "primary_score_metric_id", "primary_score_value",
        "primary_score_unit", "primary_score_direction", "primary_score_method",
        "CNNaffinity", "CNNscore", "minimizedAffinity", "blueprint_match_str",
        "desired_blueprint_score", "steric_clash_count", "strain_energy_kcal",
        "qed", "sa_score", "risk_flags", "modification_hint",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in records:
            scores = item.get("docking_scores") or {}
            posecheck = item.get("posecheck") or {}
            blueprint = item.get("blueprint_match") or {}
            prop = item.get("property") or {}
            required_hits = blueprint.get(
                "required_hits", blueprint.get("required_met", "")
            )
            required_total = blueprint.get("required_total", "")
            blueprint_summary = (
                f"{required_hits}/{required_total}R" if required_total else ""
            )
            desired_score = blueprint.get("desired_weighted_score")
            if desired_score is not None:
                blueprint_summary = (
                    f"{blueprint_summary} {float(desired_score):.0%}D"
                ).strip()
            primary_score = item.get("primary_score") or {}
            flags = prop.get("flags") or []
            writer.writerow({
                "rank": item.get("rank"),
                "mol_name": item.get("mol_name"),
                "smiles": item.get("smiles", ""),
                "selected_pose_index": item.get("selected_pose_index"),
                "selection_rank": item.get("selection_rank"),
                "source_pose_index": item.get("source_pose_index"),
                "pose_id": item.get("pose_id"),
                "chemical_identity_id": item.get("chemical_identity_id"),
                "docking_artifact_set_id": item.get(
                    "docking_artifact_set_id"
                ),
                "classification": item.get("classification"),
                "cohort_id": item.get("cohort_id"),
                "protein_state_id": item.get("protein_state_id"),
                "pocket_id": item.get("pocket_id"),
                "round_id": item.get("round_id"),
                "protocol_profile_id": item.get("protocol_profile_id"),
                "metric_catalog_version": item.get("metric_catalog_version"),
                "normalization_pool": item.get("normalization_pool"),
                "rank_partition_id": item.get("rank_partition_id"),
                "rank_scope": item.get("rank_scope"),
                "ranking_schema": item.get("ranking_schema"),
                "scoring_schema_version": item.get("scoring_schema_version"),
                "intrinsic_classification": item.get(
                    "intrinsic_classification", item.get("classification")
                ),
                "nomination_status": item.get("nomination_status", "eligible"),
                "eval_verdict": item.get("eval_verdict"),
                "prop_verdict": item.get("prop_verdict"),
                "metric_gate": item.get("metric_gate"),
                "coverage_status": item.get("coverage_status"),
                "rank_status": item.get("rank_status"),
                "overall_score": item.get("overall_score"),
                "overall_rank": item.get("overall_rank"),
                "step4_score": item.get("step4_score"),
                "step4_rank": item.get("step4_rank"),
                "step5_score": item.get("step5_score"),
                "step5_rank": item.get("step5_rank"),
                "step6_score": item.get("step6_score"),
                "step6_rank": item.get("step6_rank"),
                "category_weighted_coverage": item.get(
                    "category_weighted_coverage"
                ),
                "weighted_total_score": item.get("weighted_total_score"),
                "weighted_score_max": item.get("weighted_score_max"),
                "weighted_score_fraction": item.get("weighted_score_fraction"),
                "weighted_score_rank": item.get("weighted_score_rank"),
                "score_coverage": item.get("score_coverage"),
                "ranking_blocked": item.get("ranking_blocked", False),
                "ranking_blocked_by": "; ".join(
                    item.get("ranking_blocked_by") or []
                ),
                "primary_score_metric_id": item.get("primary_score_metric_id"),
                "primary_score_value": primary_score.get("value"),
                "primary_score_unit": primary_score.get("unit"),
                "primary_score_direction": primary_score.get("direction"),
                "primary_score_method": primary_score.get("method"),
                "CNNaffinity": scores.get("CNNaffinity"),
                "CNNscore": scores.get("CNNscore"),
                "minimizedAffinity": scores.get("minimizedAffinity"),
                "blueprint_match_str": blueprint_summary,
                "desired_blueprint_score": desired_score,
                "steric_clash_count": posecheck.get("steric_clash_count"),
                "strain_energy_kcal": posecheck.get("strain_energy_kcal"),
                "qed": _property_value(prop, "qed", "QED"),
                "sa_score": _property_value(prop, "sa_score", "SA_score"),
                "risk_flags": "; ".join(flags) if flags else "None",
                "modification_hint": item.get("modification_hint", ""),
            })


def _write_round_report(path: Path, report: dict) -> None:
    from candidate_review import build_candidate_reviews

    summary = report.get("round_summary") or {}
    validation = report.get("screening_validation") or {}
    gates = report.get("quality_gates") or {}
    all_molecules = report.get("molecules", [])
    candidates = report.get("candidate_molecules", [])
    keep_candidates = [
        m for m in candidates
        if _eligible_keep(m)
    ]
    # Gate first, then score.  The ranking table holds every molecule that
    # advanced -- keep and modify -- because that is the pool the scores were
    # computed over; rejects never enter it.  Matching ``_eligible_keep`` on the
    # status test keeps this table and the round summary counting the same rows:
    # an annotation such as ``coverage_review`` is not a nomination block.
    ranked_candidates = [
        m for m in candidates
        if m.get("classification") in {"keep", "modify"}
        and not str(m.get("nomination_status") or "").startswith("blocked")
    ]
    blocked_keep_candidates = [
        m for m in candidates
        if m.get("intrinsic_classification") == "keep"
        and str(m.get("nomination_status", "")).startswith("blocked_by_")
    ]
    decoy_keep = sum(
        1 for m in all_molecules
        if _is_decoy_record(m) and m.get("classification") == "keep"
    )
    high_strain_keep = sum(
        1 for m in keep_candidates
        if _number_or_none(((m.get("posecheck") or {}).get("strain_energy_kcal"))) is not None
        and _number_or_none(((m.get("posecheck") or {}).get("strain_energy_kcal")))
        > _report_metric_bound(
            report, "posecheck.strain_energy_kcal", "soft_upper", 8.0
        )
    )
    high_clash_keep = sum(
        1 for m in keep_candidates
        if _number_or_none(((m.get("posecheck") or {}).get("steric_clash_count"))) is not None
        and _number_or_none(((m.get("posecheck") or {}).get("steric_clash_count")))
        > _report_metric_bound(
            report, "posecheck.steric_clash_count", "soft_upper", 0.0
        )
    )
    lines = [
        f"# SBDD Round {summary.get('round', '')} Candidate Report".strip(),
        "",
        "## Executive Summary",
        "",
        f"- Round assessment: `{(report.get('round_assessment') or {}).get('status')}`",
        f"- Eligible candidates available: `{report.get('valid_for_nomination')}`; campaign planning is owned by the Agent.",
        f"- Quality gate status: `{summary.get('quality_gate_status')}` "
        f"(passed: `{summary.get('quality_gate_passed')}`)",
        f"- Screened compounds: {summary.get('total_evaluated', 0)} "
        f"(keep/modify/reject: {summary.get('keep', 0)}/{summary.get('modify', 0)}/{summary.get('reject', 0)})",
        f"- Non-decoy candidates screened: "
        f"{summary.get('candidate_total', len(candidates))}; "
        f"{summary.get('candidate_keep', len(keep_candidates))} eligible keep "
        f"(candidate pass rate {summary.get('candidate_pass_rate', 0):.1%})",
        f"- Non-rejected candidate ranking entries: {len(candidates)}",
        f"- Intrinsic molecule-level pass before protocol gating: "
        f"{summary.get('candidate_intrinsic_keep', 0)} candidate keep "
        f"(blocked by quality gates: {summary.get('candidate_quality_gate_blocked_keep', 0)})",
        f"- Benchmark decoys classified as keep: {decoy_keep} (used for validation only, not nomination)",
        f"- Main bottleneck: `{summary.get('primary_bottleneck')}`",
    ]
    if summary.get("fix_rate") is not None:
        lines.append(f"- Fix rate (modify→keep): `{summary['fix_rate']}`")
    if validation.get("available"):
        ci = validation.get("AUROC_bootstrap_95ci", [])
        lines.extend([
            "",
            "## Retrospective Screening Validation",
            "",
            f"- {validation.get('score', 'primary score')} AUROC: "
            f"{validation.get('AUROC')} (bootstrap 95% CI: {ci})",
            f"- Benchmark size: {validation.get('num_actives')} actives + {validation.get('num_decoys')} decoys",
            f"- EF1%/EF5%/EF10%: {validation.get('EF1%')}/{validation.get('EF5%')}/{validation.get('EF10%')}",
            f"- PR-AUC: {validation.get('PR_AUC')} "
            f"(prevalence {validation.get('active_prevalence')}, "
            f"relative lift {validation.get('PR_AUC_relative_to_prevalence')})",
            f"- Top 1%/5%/10% decoy fractions: "
            f"{validation.get('top1%_decoy_fraction')}/"
            f"{validation.get('top5%_decoy_fraction')}/"
            f"{validation.get('top10%_decoy_fraction')}",
        ])
    lines.extend([
        "",
        "## Quality Gates",
        "",
        "| Gate | Status | Key detail |",
        "|---|---:|---|",
    ])
    for name, gate in gates.items():
        details = gate.get("details") or {}
        detail = ""
        if name == "native_redocking":
            detail = f"best RMSD {details.get('best_rmsd_angstrom')} Å"
        elif name == "pocket_validation":
            detail = (
                f"known-site recall {details.get('known_site_recall')}, "
                f"IoU {details.get('iou_vs_known_site')}"
            )
        elif name == "blueprint_validation":
            detail = (
                f"blueprint valid: {gate.get('blueprint_valid')}; "
                f"Required exact gate: {gate.get('required_hard_gate_status')}; "
                f"Desired metrics: "
                f"{'declared' if gate.get('desired_metric_applicable') else 'not applicable'}"
            )
        elif name == "screening_benchmark":
            detail = f"AUROC {details.get('AUROC')}"
        lines.append(
            f"| {name} | {_quality_gate_status(gate)} | {detail} |"
        )

    lines.extend([
        "",
        "## Candidate Pool Notes",
        "",
        "- Candidate ranking below excludes benchmark decoys (`DEC_*`) and "
        "rejected molecules; rejects stay in `molecules` and the published "
        "tables as the audit trail for why they were dropped.",
        "- Overall ranks are local to one scientific cohort and normalization "
        "pool; ranks from different cohorts are not comparable.",
        "- The report shows up to 25 rows independently for each cohort.",
        f"- Among eligible keep candidates, {high_clash_keep} exceed the PoseCheck "
        f"soft clash bound and {high_strain_keep} exceed the soft strain bound.",
        "- Treat high-clash or high-strain molecules as optimization leads, not project-level picks.",
        "",
        "## Top Ranked Candidates (keep + modify)",
        "",
    ])
    _append_candidate_report_tables(
        lines,
        ranked_candidates,
        limit_per_partition=25,
        empty_message="No keep or modify candidates advanced to ranking.",
    )
    _append_candidate_review_tables(
        lines, build_candidate_reviews(candidates, source_round=summary.get("round"))
    )

    exhausted = [
        record for record in all_molecules
        if record.get("pose_exhaustion") and record.get("classification") == "reject"
    ]
    if exhausted:
        lines.extend([
            "", "## Molecules Skipped: All Poses Failed", "",
            "Every input pose failed QC under the current protocol. These molecules "
            "were skipped without redocking (`all_poses_failed`).", "",
            "| Molecule | Evaluated poses | QC failure reasons |",
            "|---|---:|---|",
        ])
        for record in exhausted:
            exhaustion = record["pose_exhaustion"]
            reasons = "; ".join(exhaustion.get("reject_reasons") or [])
            lines.append(
                f"| {record['mol_name']} | {exhaustion['evaluated_poses']} | {reasons} |"
            )

    if blocked_keep_candidates:
        lines.extend([
            "",
            "## Protocol-Blocked Keep Candidates",
            "",
            "These molecules passed molecule-level keep criteria but are not nominated because experiment quality gates failed.",
            "",
        ])
        _append_candidate_report_tables(
            lines,
            blocked_keep_candidates,
            limit_per_partition=25,
        )
    lines.extend([
        "",
        "## Campaign Review",
        "",
        "The Agent reviews this round's evidence, diagnostics, and campaign budget to plan the next action.",
        "Eligible keep candidates can be inspected while the campaign continues; this report does not declare campaign convergence.",
        "",
        "Property evidence, when present, is joined to structure and docking evidence only by the exact canonical pose_id.",
        "",
    ])
    path.write_text("\n".join(lines))


def _append_candidate_review_tables(lines: list[str], reviews: dict) -> None:
    lines.extend([
        "", "## KEEP and MODIFY Top 1–5 for Review", "",
        f"Source round: {reviews.get('source_round')}. Review only: each class is "
        "ordered by its existing OverallScore within the same cohort and pool. "
        "Review # is a display position, not OverallRank. No scores are "
        "recalculated, blocked rows remain blocked, and MODIFY rows are not nominations.", "",
    ])
    for group in reviews.get("groups") or []:
        lines.extend([
            f"### {group['classification'].upper()} Top 1–5 · Cohort `{group['cohort_id']}` · `{group['normalization_pool']}`",
            "",
            f"{group['scored_count']} scored; {group['unscored_count']} without a score (excluded from this review order).",
            "",
            "| Review # | Molecule | Existing OverallScore | OverallRank | Molecular category findings | Evidence gaps / protocol blockers | Rank status / nomination status |",
            "|---:|---|---:|---:|---|---|---|",
        ])
        for row in group["rows"]:
            evidence = list(row["evidence_gaps"]) + list(row["protocol_blockers"])
            cells = [
                row["review_position"], row["mol_name"], row["overall_score"], row["overall_rank"],
                "; ".join(row["molecular_findings"]) or "No flagged category evidence recorded",
                "; ".join(evidence) or "No gaps recorded",
                f"{row['rank_status'] or 'unreported'} / {row['nomination_status'] or 'unreported'}",
            ]
            lines.append("| " + " | ".join(
                str(cell).replace("|", "\\|").replace("\n", " ") if cell is not None else "—"
                for cell in cells
            ) + " |")
        if not group["rows"]:
            lines.append("| — | No candidates with an existing OverallScore. | | | | | |")
        lines.append("")


def _append_candidate_report_tables(
    lines: list[str],
    records: list[dict],
    *,
    limit_per_partition: int,
    empty_message: str | None = None,
) -> None:
    """Append independently truncated tables for non-comparable rank pools."""
    partitions: dict[tuple[str, str], list[dict]] = {}
    for record in records:
        cohort_id = str(record.get("cohort_id") or "cohort:unspecified")
        normalization_pool = str(
            record.get("normalization_pool") or "nomination_candidates"
        )
        partitions.setdefault((cohort_id, normalization_pool), []).append(record)
    if not partitions:
        if empty_message:
            lines.append(empty_message)
        return

    for cohort_id, normalization_pool in sorted(partitions):
        lines.extend([
            f"### Cohort `{cohort_id}` · `{normalization_pool}`",
            "",
            "| Overall rank | ID | Gate | Overall score | Rank status | Category coverage | Primary metric | Primary score | Blueprint | Clash | Strain | QED | SA | Flags | SMILES |",
            "|---:|---|---|---:|---|---:|---|---:|---|---:|---:|---:|---:|---|---|",
        ])
        partition_records = sorted(
            partitions[(cohort_id, normalization_pool)],
            key=lambda item: (
                _number_or_none(item.get("overall_rank")) is None,
                (
                    _number_or_none(item.get("overall_rank"))
                    if _number_or_none(item.get("overall_rank")) is not None
                    else float("inf")
                ),
                str(item.get("mol_name") or ""),
            ),
        )
        for item in partition_records[:limit_per_partition]:
            lines.append(_candidate_table_row(item))
        lines.append("")


def _candidate_table_row(item: dict) -> str:
    posecheck = item.get("posecheck") or {}
    bp = item.get("blueprint_match") or {}
    prop = item.get("property") or {}
    smiles = item.get("smiles", "")
    primary_score = item.get("primary_score") or {}
    primary_metric = (
        item.get("primary_score_metric_id")
        or primary_score.get("metric_id")
        or ""
    )
    primary_value = primary_score.get("value", "")
    req_met = bp.get("required_hits", bp.get("required_met", ""))
    req_total = bp.get("required_total", "")
    bp_str = f"{req_met}/{req_total} R" if req_total else item.get("eval_verdict", "")
    desired_score = bp.get("desired_weighted_score")
    if desired_score is not None:
        bp_str = f"{bp_str}; {float(desired_score):.0%} D"
    clash = posecheck.get("steric_clash_count", "")
    strain = posecheck.get("strain_energy_kcal", "")
    qed = _property_value(prop, "qed", "QED")
    sa = _property_value(prop, "sa_score", "SA_score")
    flags = prop.get("flags", [])
    risk_str = "; ".join(flags) if flags else "None"
    return (
        f"| {item.get('overall_rank')} | {item.get('mol_name')} | "
        f"{item.get('classification')} | "
        f"{item.get('overall_score')} | {item.get('rank_status')} | "
        f"{item.get('category_weighted_coverage')} | {primary_metric} | "
        f"{primary_value} | {bp_str} | {clash} | {strain} | "
        f"{qed} | {sa} | {risk_str} | {smiles} |"
    )


def _number_or_none(value: object) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _metric_bound(metric_id: str, field: str, default: float) -> float:
    for rule in POSE_PROPERTY_METRIC_RULES:
        if rule.get("metric_id") == metric_id:
            value = _number_or_none(rule.get(field))
            return value if value is not None else default
    return default


def _report_metric_bound(
    report: dict,
    metric_id: str,
    field: str,
    default: float,
) -> float:
    for rule in report.get("metric_rules") or []:
        if rule.get("metric_id") == metric_id:
            value = _number_or_none(rule.get(field))
            return value if value is not None else default
    return _metric_bound(metric_id, field, default)


def main():
    from io_utils import parse_skill_args, print_skill_result
    input_data, output_dir = parse_skill_args("Reflective Coordinator Skill")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

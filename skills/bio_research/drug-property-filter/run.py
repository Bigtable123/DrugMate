#!/usr/bin/env python3
"""Step 5 property filtering with isolated ADMET, ChemFH, and synthesis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from admet_utils import run_admet_worker
from admet_contract import (
    ADMET_NATIVE_ENDPOINT_IDS,
    validate_native_predictions,
)
from chemfh_utils import run_chemfh_worker
from rascore_utils import attach_rascore
from property_utils import (
    DEFAULT_ALERT_CATALOGS,
    DEFAULT_HARD_FAIL_CATALOGS,
    DEFAULT_LIPINSKI_MAX_VIOLATIONS,
    check_lipinski,
    check_pains_filter_available,
    check_rdkit_available,
    check_sa_score_available,
    check_structural_alerts,
    compute_drug_properties,
    compute_verdict,
)
from synthesis_utils import run_aizynth_worker
from metric_scoring import (
    CATEGORY_PLANNED_COVERAGE_FLOOR,
    DEFAULT_SCREENING_MODE,
    METRIC_RULES_VERSION,
    NORMALIZATION_METHOD,
    PROPERTY_FILTER_METRIC_RULES,
    SCORING_SCHEMA_VERSION,
    aggregate_molecule_gate,
    assess_record,
    extract_metric_value,
    public_metric_rules,
    record_has_structural_reject,
    resolve_screening_mode,
    score_records,
)
from pose_identity import LIGAND_LINEAGE_PROPERTIES, benchmark_label, complete_ligand_ids


LINEAGE_PROPERTIES = (
    "pose_id",
    "parent_pose_id",
    "coordinate_revision",
    "chemical_identity_id",
    "chemical_id",
    "source_identity_id",
    "source_identity",
    "lineage_metadata",
    "source_pose_index",
    "selection_rank",
    "compound_pose_rank",
    "cohort_id",
    "protein_state_id",
    "pocket_id",
    "round_id",
    "protocol_profile_id",
    "metric_catalog_version",
    "cohort_identity_complete",
    "cohort_identity_missing_fields",
    "cohort_identity_derived_fields",
    "cohort_provenance_sources",
    *LIGAND_LINEAGE_PROPERTIES,
)
JSON_PROPERTIES = {
    "source_identity",
    "lineage_metadata",
    "cohort_identity_complete",
    "cohort_identity_missing_fields",
    "cohort_identity_derived_fields",
    "cohort_provenance_sources",
}
INTEGER_PROPERTIES = {
    "coordinate_revision",
    "source_pose_index",
    "selection_rank",
    "compound_pose_rank",
}
# A small backtrack affects a single candidate, so one malformed record cannot
# discard the batch. A common integrity problem affecting much of the input
# indicates a possible upstream writer failure and requires review instead of
# continuing with suspect lineage keys. Escalation needs both a real fraction of the
# batch and more than one record.
SYSTEMIC_INTEGRITY_ABORT_RATIO = 0.5
# A graph conflict always implicates at least two records by construction, so a
# ratio alone would let one bad pair halt a small batch. Escalation additionally
# needs a third affected record before the problem is treated as systemic.
SYSTEMIC_INTEGRITY_ABORT_MIN_RECORDS = 3
# ``not_evaluated`` ranks between ``keep`` and ``modify``: no metric rule could
# label the molecule, which is not a molecular defect (so it must never dominate
# a real ``modify``/``reject``) but is also not a clean pass.
CLASSIFICATION_SEVERITY = {
    "keep": 0,
    "not_evaluated": 1,
    "modify": 2,
    "reject": 3,
}
CHEMFH_ENDPOINT_IDS = (
    "other_assay_interference",
    "blue_fluorescence",
    "fluc_inhibitors",
    "promiscuous_compounds",
    "colloidal_aggregators",
    "reactive_compounds",
    "green_fluorescence",
)
CHEMFH_RULE_NAMES = (
    "Aggregators",
    "Fluc",
    "Blue_fluorescence",
    "Green_fluorescence",
    "Reactive",
    "Other_assay_interference",
    "Promiscuous",
    "ALARM_NMR",
    "BMS",
    "Chelator_Rule",
    "GST_FHs_Rule",
    "His_FHs_Rule",
    "Luciferase_Inhibitor_Rule",
    "NTD",
    "PAINS",
    "Potential_Electrophilic_Rule",
    "Lilly",
)
CHEMFH_SUBSTRUCTURE_RULE_NAMES = CHEMFH_RULE_NAMES[7:]
CHEMFH_SUBSTRUCTURE_RULE_MAPPING_VERSION = (
    "chemfh-native-substructure-10-v1"
)
CHEMFH_POSITIVE_THRESHOLD = 0.5
ADMET_PARTIAL_COVERAGE_CATEGORIES = ("ABS", "TOX")
PROPERTY_PLAN_FIELDS = (
    "planned_metric_ids",
    "required_metric_ids",
    "required_admet_metric_ids",
    "not_applicable_metric_ids",
    "planned_category_codes",
    "planned_tools",
    "required_tools",
    "optional_tools",
    "failed_tools",
    "planned_but_failed_tools",
    "failed_metric_ids",
    "planned_but_failed_metric_ids",
    "operational_errors",
    "equivalent_fallbacks",
    "active_steps",
)
SYNTHESIS_ATTEMPTED_STATUSES = frozenset({
    "solved", "not_found_within_budget", "error",
})
SYNTHESIS_NOT_ATTEMPTED_STATUSES = frozenset({
    "disabled", "not_eligible", "not_selected_top_k",
})
AIZYNTHFINDER_METRIC_IDS = tuple(
    str(rule["metric_id"])
    for rule in PROPERTY_FILTER_METRIC_RULES
    if rule.get("requires_property_tool") == "aizynthfinder"
)
SYNTHESIS_RETRY_ACTION = (
    "Rerun or repair AiZynthFinder and capture a complete route-review result"
)


def run(input_data: dict, output_dir: Path) -> dict:
    """Run RDKit gates, batch property backends, and top-K retrosynthesis."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if not check_rdkit_available():
        return _failure_report(
            output_dir, "RDKit is required in the drugmate environment"
        )

    molecules_sdf = input_data.get("molecules_sdf")
    molecules_smiles = input_data.get("molecules_smiles")
    if not molecules_sdf and molecules_smiles is None:
        return _failure_report(
            output_dir, "molecules_sdf or molecules_smiles is required"
        )

    try:
        records = _load_records(
            molecules_sdf=molecules_sdf,
            molecules_smiles=molecules_smiles,
            remove_hs=_remove_hs_setting(input_data),
        )
    except (OSError, ValueError) as exc:
        return _failure_report(output_dir, str(exc))

    if not records:
        if input_data.get("allow_empty_input") is True:
            return _write_final_report(
                output_dir,
                success=True,
                summary={
                    "input_records": 0, "valid_input_records": 0,
                    "invalid_input_records": 0, "valid_records": 0,
                    "invalid_records": 0, "unique_candidates": 0,
                    "deduplicated_records": 0, "duplicate_pose_records": 0,
                },
                results=[],
                steps_completed=[],
                steps_skipped=[{
                    "step": "property_filter", "reason_code": "no_selected_poses",
                    "reason": "No poses survived QC; no molecules to evaluate.",
                }],
                parameters={"allow_empty_input": True},
                backend_reports={},
            )
        return _failure_report(
            output_dir,
            "Input contains zero molecule records",
            summary={
                "input_records": 0,
                "valid_input_records": 0,
                "invalid_input_records": 0,
                "unique_candidates": 0,
            },
        )

    _apply_input_cohort_defaults(records, input_data)

    metadata_errors = [
        {
            "record_index": record["record_index"],
            "mol_name": record["mol_name"],
            **error,
        }
        for record in records
        for error in (record.get("metadata_errors") or [])
    ]
    metadata_integrity = _integrity_scope(
        records,
        affected=sum(1 for record in records if record.get("metadata_errors")),
    )
    if metadata_integrity["systemic"]:
        return _failure_report(
            output_dir,
            "Identity/ranking/lineage metadata could not be parsed safely for "
            f"{metadata_integrity['affected_records']} of "
            f"{metadata_integrity['total_records']} records, which is a writer "
            "defect rather than a per-molecule problem",
            summary={
                "input_records": len(records),
                "valid_input_records": sum(
                    record["mol"] is not None for record in records
                ),
                "invalid_input_records": sum(
                    record["mol"] is None for record in records
                ),
                "unique_candidates": 0,
                "metadata_integrity_errors": len(metadata_errors),
            },
            diagnostics={
                "metadata_integrity_errors": metadata_errors,
                "metadata_integrity_scope": metadata_integrity,
            },
        )

    (
        enabled_checks,
        admet_enabled,
        chemfh_enabled,
        dedup_enabled,
        synthesis_enabled,
    ) = (
        _resolve_property_tools(input_data)
    )
    scoring_plan = _input_property_scoring_plan(input_data)
    thresholds = {
        "qed_threshold": _number(input_data.get("qed_threshold"), 0.45),
        "sa_threshold": _number(input_data.get("sa_threshold"), 6.0),
        "lipinski_max_violations": _integer(
            input_data.get("lipinski_max_violations"),
            DEFAULT_LIPINSKI_MAX_VIOLATIONS,
        ),
        "bbb_profile": str(
            (
                (_tool_config(input_data, "admet").get("params") or {}).get(
                    "bbb_profile",
                    input_data.get("bbb_profile", "non_cns"),
                )
            )
        ).strip().lower(),
    }
    property_metric_rules = _effective_property_metric_rules(thresholds)
    # Level 2 of the gate contract: which per-category fail count/ratio turns a
    # set of level-1 marks into a rejection.  ``candidate_source``
    # only supplies the default when no mode is named explicitly.
    screening_mode = resolve_screening_mode(
        input_data.get("screening_mode"),
        context=input_data.get("candidate_source"),
    )
    alert_config = _tool_config(input_data, "pains")
    alert_params = dict(alert_config.get("params") or {})
    alert_catalogs = alert_params.get("catalogs") or list(DEFAULT_ALERT_CATALOGS)
    hard_fail_catalogs = (
        alert_params.get("hard_fail_catalogs")
        if "hard_fail_catalogs" in alert_params
        else list(DEFAULT_HARD_FAIL_CATALOGS)
    )
    if hard_fail_catalogs is None:
        hard_fail_catalogs = list(DEFAULT_HARD_FAIL_CATALOGS)
    is_covalent = bool(input_data.get("is_covalent", False))

    identity_conflicts = _identity_conflicts(records)
    identity_scope = _apply_identity_conflict_quarantine(records, identity_conflicts)
    if identity_scope["systemic"]:
        return _failure_report(
            output_dir,
            "Conflicting molecular graphs were declared under the same chemical "
            f"identity for {identity_scope['affected_records']} of "
            f"{identity_scope['total_records']} records; the identity namespace "
            "itself is unreliable",
            summary={
                "input_records": len(records),
                "valid_input_records": sum(record["mol"] is not None for record in records),
                "invalid_input_records": sum(record["mol"] is None for record in records),
                "identity_conflicts": len(identity_conflicts),
            },
            diagnostics={
                "identity_conflicts": identity_conflicts,
                "identity_conflict_scope": identity_scope,
            },
        )

    groups = _group_records(records, deduplicate=dedup_enabled)
    results: list[dict] = []
    representatives: dict[str, dict] = {}
    valid_candidates = 0
    alert_backend_available = check_pains_filter_available()
    sa_available = check_sa_score_available()

    for group in groups:
        representative = group["representative"]
        result = _result_identity(group)
        enabled_metric_tools = set(enabled_checks)
        if admet_enabled:
            enabled_metric_tools.add("admet")
        if chemfh_enabled:
            enabled_metric_tools.add("chemfh")
        if synthesis_enabled:
            enabled_metric_tools.add("aizynthfinder")
        if ((input_data.get("tools") or {}).get("rascore") or {}).get("enabled", False):
            enabled_metric_tools.add("rascore")
        result["metric_context"] = {
            "enabled_tools": sorted(enabled_metric_tools),
            "is_covalent": is_covalent,
            "hard_fail_catalogs": list(hard_fail_catalogs),
            **deepcopy(scoring_plan),
        }
        if scoring_plan:
            result["scoring_plan"] = deepcopy(scoring_plan)
        result["is_covalent"] = is_covalent
        representatives[result["candidate_id"]] = representative
        group_metadata_errors = [
            {
                "record_index": member["record_index"],
                "mol_name": member["mol_name"],
                **error,
            }
            for member in group["records"]
            for error in (member.get("metadata_errors") or [])
        ]
        if group_metadata_errors:
            # Malformed lineage/ranking properties are a writer defect, not a
            # molecular defect, so the molecule keeps its three-tier
            # classification and only loses automatic nomination.
            result["metadata_integrity_errors"] = group_metadata_errors
        identity_conflict = [
            match
            for member in group["records"]
            for match in (member.get("identity_conflict") or [])
        ]
        if identity_conflict:
            # The executable policy admits only identity/valid_structure as hard rejects,
            # and the executable policy rejects the failing item, not the batch.  We
            # cannot tell which graph the identity names, so the candidate is
            # rejected here rather than silently scored under one of them.
            result.update({
                "valid_structure": True,
                "identity_reject": True,
                "structural_reject": True,
                "structural_reject_reason": "identity_conflict",
                "identity_conflict": identity_conflict,
                "properties": None,
                "lipinski": None,
                "structural_alerts": None,
                "pains": None,
                "admet": None,
                "chemfh": None,
                "synthesis": None,
                "verdict": "fail",
                "flags": [
                    "Chemical identity is associated with conflicting molecular "
                    "graphs"
                ],
                "suggestions": [
                    "Repair the identity assignment upstream: one "
                    "chemical_identity_id must name exactly one molecular graph"
                ],
                "nomination_eligibility": "rejected",
            })
            results.append(result)
            continue
        mol = representative.get("mol")
        if mol is None:
            result.update({
                "valid_structure": False,
                "structural_reject": True,
                "structural_reject_reason": "invalid_structure",
                "properties": None,
                "lipinski": None,
                "structural_alerts": None,
                "pains": None,
                "admet": None,
                "chemfh": None,
                "synthesis": None,
                "verdict": "fail",
                "flags": ["Molecule parse failed"],
                "suggestions": ["Repair or remove the invalid molecule record"],
                "nomination_eligibility": "rejected",
            })
            results.append(result)
            continue

        try:
            properties = compute_drug_properties(mol)
            lipinski = check_lipinski(properties)
            alerts = (
                check_structural_alerts(mol, alert_catalogs)
                if "pains" in enabled_checks
                else {
                    "available": None,
                    "catalogs": [],
                    "has_alert": None,
                    "matches": [],
                    "disabled": True,
                }
            )
            pains = _legacy_pains_projection(alerts)
            verdict, flags, suggestions, rdkit_evidence = compute_verdict(
                properties,
                lipinski,
                pains,
                thresholds,
                enabled_checks=enabled_checks,
                structural_alerts=alerts,
                hard_fail_catalogs=hard_fail_catalogs,
                is_covalent=is_covalent,
            )
        except Exception as exc:
            # The molecule itself parsed (``mol is not None``); what failed is
            # our own descriptor/verdict computation. The executable policy forbids
            # rejecting a molecule for a metric or tool coverage problem, so
            # this must NOT reuse the structural branch above: claiming
            # ``valid_structure: False`` would be a false statement about the
            # chemistry and would remove the row from every percentile pool.
            # It becomes an evidence gap instead -- no label, no nomination.
            result.update({
                "valid_structure": True,
                "properties": None,
                "lipinski": None,
                "structural_alerts": None,
                "pains": None,
                "admet": None,
                "chemfh": None,
                "synthesis": None,
                "verdict": "not_evaluated",
                "property_tool_error": {
                    "stage": "rdkit_property_pass",
                    "error": f"{type(exc).__name__}: {exc}",
                },
                "flags": [f"RDKit property calculation failed: {type(exc).__name__}"],
                "suggestions": [
                    "Repair the property-filter tool configuration or inputs and "
                    "rerun; the molecule is unscored, not rejected",
                    str(exc),
                ],
                "nomination_eligibility": "review_only",
                "molecule_review_reasons": [],
                "protocol_review_reasons": ["rdkit_property_pass_failed"],
                "nomination_eligibility_reasons": ["rdkit_property_pass_failed"],
            })
            results.append(result)
            continue

        valid_candidates += 1
        result.update({
            "valid_structure": True,
            "properties": properties,
            "lipinski": lipinski,
            "structural_alerts": alerts,
            "pains": pains,
            "admet": None,
            "chemfh": None,
            "synthesis": None,
            "verdict": verdict,
            "flags": flags,
            "suggestions": suggestions,
            "rdkit_evidence": rdkit_evidence,
            "nomination_eligibility": "review_only",
        })
        results.append(result)

    steps_completed, steps_skipped = _rdkit_step_status(
        valid_candidates,
        enabled_checks,
        sa_available=sa_available,
        alert_backend_available=alert_backend_available,
    )
    backend_reports: dict[str, dict] = {}

    if valid_candidates == 0:
        metric_policy = _apply_property_metric_policy(
            results,
            rules=property_metric_rules,
            screening_mode=screening_mode,
        )
        summary = _build_summary(
            records,
            groups,
            results,
            admet_enabled,
            chemfh_enabled,
            synthesis_enabled,
        )
        return _write_final_report(
            output_dir,
            success=False,
            error="All molecule records are invalid or failed RDKit processing",
            summary=summary,
            results=results,
            steps_completed=steps_completed,
            steps_skipped=steps_skipped,
            parameters=_parameters(
                thresholds,
                alert_catalogs,
                hard_fail_catalogs,
                is_covalent,
                dedup_enabled,
                input_data=input_data,
            ),
            backend_reports=backend_reports,
            metric_policy=metric_policy,
        )

    if admet_enabled:
        admet_report = _run_admet(
            input_data, output_dir, results, steps_completed, steps_skipped
        )
        backend_reports["admet_ai"] = _backend_provenance(admet_report)
    else:
        for result in results:
            if result.get("valid_structure"):
                result["admet"] = {"status": "disabled"}
        steps_skipped.append({
            "step": "admet",
            "reason": (
                "ADMET-AI not selected for this run; the categories it owns "
                "leave the planned coverage denominator"
            ),
        })

    if chemfh_enabled:
        chemfh_report = _run_chemfh(
            input_data, output_dir, results, steps_completed, steps_skipped
        )
        backend_reports["chemfh"] = _backend_provenance(chemfh_report)
    else:
        for result in results:
            if result.get("valid_structure"):
                result["chemfh"] = {"status": "disabled"}
        steps_skipped.append({
            "step": "chemfh",
            "reason": (
                "ChemFH not selected for this run; its zero-weight gates are "
                "N/A unless the scoring plan declares ChemFH required"
            ),
        })

    backend_reports["rascore"] = attach_rascore(
        input_data, output_dir, results, steps_completed, steps_skipped
    )

    # Backends attach evidence first. A single conservative reducer owns
    # promotion so one backend can never overwrite another's cap.  The
    # pre-synthesis fold intentionally excludes SYN: AiZynthFinder values are
    # still pending, and the ratio branch is not monotonic while its denominator
    # is incomplete.  Completed non-SYN categories may still reject before the
    # expensive route search; the complete SYN category is folded afterwards.
    _update_property_backend_eligibility(
        results,
        rules=property_metric_rules,
    )
    _apply_pre_synthesis_metric_gates(
        results, property_metric_rules, screening_mode=screening_mode
    )

    synthesis_report = _run_synthesis(
        input_data,
        output_dir,
        results,
        enabled=synthesis_enabled,
        # Retrosynthesis is intentionally downstream of complete ADMET even in
        # standalone calls; registry validation enforces the same dependency.
        admet_required=True,
        steps_completed=steps_completed,
        steps_skipped=steps_skipped,
    )
    if synthesis_report.get("backend"):
        backend_reports["aizynthfinder"] = _backend_provenance(synthesis_report)

    metric_policy = _apply_property_metric_policy(
        results,
        rules=property_metric_rules,
        screening_mode=screening_mode,
    )
    summary = _build_summary(
        records,
        groups,
        results,
        admet_enabled,
        chemfh_enabled,
        synthesis_enabled,
    )
    filtered_files: dict[str, str] = {}
    success = True
    error = None
    if input_data.get("write_filtered_sdf", False):
        filtered = _write_filtered_outputs(
            input_data,
            output_dir,
            representatives,
            results,
        )
        summary.update(filtered["summary"])
        filtered_files = filtered["files"]
        if (
            input_data.get("fail_if_no_filtered_molecules", False)
            and summary.get("filtered_total", 0) == 0
        ):
            success = False
            error = "No molecules remained after property prefilter"

    return _write_final_report(
        output_dir,
        success=success,
        error=error,
        summary=summary,
        results=results,
        steps_completed=steps_completed,
        steps_skipped=steps_skipped,
        parameters=_parameters(
            thresholds,
            alert_catalogs,
            hard_fail_catalogs,
            is_covalent,
            dedup_enabled,
            input_data=input_data,
        ),
        backend_reports=backend_reports,
        metric_policy=metric_policy,
        synthesis_report=synthesis_report,
        extra_files=filtered_files,
    )


def _load_records(
    *,
    molecules_sdf: str | None,
    molecules_smiles: Any,
    remove_hs: bool,
) -> list[dict]:
    from rdkit import Chem

    records: list[dict] = []
    if molecules_sdf:
        sdf_path = Path(str(molecules_sdf)).expanduser()
        if not sdf_path.is_file():
            raise ValueError(f"SDF not found: {sdf_path}")
        if sdf_path.stat().st_size == 0:
            return []
        supplier = Chem.SDMolSupplier(
            str(sdf_path), removeHs=bool(remove_hs), sanitize=True
        )
        for index, mol in enumerate(supplier):
            metadata = _sdf_metadata(mol)
            name = (
                str(mol.GetProp("_Name")).strip() or f"mol_{index}"
                if mol is not None and mol.HasProp("_Name")
                else f"mol_{index}"
            )
            records.append(_record(
                mol=mol,
                record_index=index,
                mol_name=name,
                input_smiles=None,
                metadata=metadata,
            ))
        return _complete_record_ids(records, hashlib.sha256(sdf_path.read_bytes()).hexdigest())

    if not isinstance(molecules_smiles, (list, tuple)):
        raise ValueError("molecules_smiles must be a list")
    for index, item in enumerate(molecules_smiles):
        if isinstance(item, dict):
            smiles = str(item.get("smiles") or "")
            name = str(item.get("mol_name") or item.get("name") or f"mol_{index}")
            metadata = {
                str(key): value
                for key, value in item.items()
                if key not in {"smiles", "mol_name", "name"}
            }
        else:
            smiles = str(item or "")
            name = f"mol_{index}"
            metadata = {}
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is not None:
            mol.SetProp("_Name", name)
        records.append(_record(
            mol=mol,
            record_index=index,
            mol_name=name,
            input_smiles=smiles,
            metadata=metadata,
        ))
    source_token = hashlib.sha256(
        json.dumps(molecules_smiles, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()
    return _complete_record_ids(records, source_token)


def _complete_record_ids(records: list[dict], source_token: str) -> list[dict]:
    for record in records:
        if record["mol"] is not None:
            for field in ("compound_id", "CompoundID", "compound_name", "input_compound_id", "input_mol_name", "source_record_id"):
                value = record["metadata"].get(field)
                if value not in (None, ""):
                    record["mol"].SetProp(field, str(value))
    complete_ligand_ids([record["mol"] for record in records], source_token=source_token)
    for record in records:
        mol = record["mol"]
        if mol is not None:
            for field in ("compound_id", "input_compound_id", "input_mol_name", "source_record_id"):
                if mol.HasProp(field):
                    record[field] = record["metadata"][field] = mol.GetProp(field)
            record["mol_name"] = record["compound_id"]
    return records


def _apply_input_cohort_defaults(records: list[dict], input_data: dict) -> None:
    """Fill absent row provenance from the adapter without overwriting SDF data."""
    defaults = {
        field: input_data.get(field)
        for field in LINEAGE_PROPERTIES
        if field.startswith("cohort_")
        or field in {
            "protein_state_id", "pocket_id", "round_id",
            "protocol_profile_id", "metric_catalog_version",
        }
    }
    source_descriptor = {
        "molecules_sdf": str(input_data.get("molecules_sdf") or ""),
        "molecules_smiles": input_data.get("molecules_smiles"),
    }
    source_token = hashlib.sha256(
        json.dumps(
            source_descriptor, sort_keys=True, default=str,
        ).encode("utf-8")
    ).hexdigest()
    component_fields = (
        "protein_state_id", "pocket_id", "round_id",
        "protocol_profile_id", "metric_catalog_version",
    )
    for record in records:
        for field, value in defaults.items():
            if record.get(field) in (None, "") and value not in (None, ""):
                record[field] = value
        missing_fields = [
            field for field in component_fields
            if record.get(field) in (None, "")
        ]
        for field in missing_fields:
            if field == "metric_catalog_version":
                record[field] = METRIC_RULES_VERSION
            else:
                record[field] = f"unknown:{field}:{source_token}"
        unresolved_fields = [
            field for field in component_fields
            if str(record.get(field) or "").startswith(("unknown:", "pending:"))
            or "|state_source:" in str(record.get(field) or "")
        ]
        if record.get("cohort_identity_complete") in (None, ""):
            record["cohort_identity_complete"] = not unresolved_fields
        if record.get("cohort_identity_missing_fields") in (None, ""):
            record["cohort_identity_missing_fields"] = unresolved_fields
        derived = record.get("cohort_identity_derived_fields")
        if not isinstance(derived, list):
            derived = []
        derived.extend(missing_fields)
        record["cohort_identity_derived_fields"] = sorted(set(derived))
        sources = record.get("cohort_provenance_sources")
        if not isinstance(sources, dict):
            sources = {}
        for field in missing_fields:
            sources.setdefault(
                field,
                "registry:metric_catalog"
                if field == "metric_catalog_version"
                else "derived:property_input_source_token",
            )
        record["cohort_provenance_sources"] = sources


def _record(
    *,
    mol: Any,
    record_index: int,
    mol_name: str,
    input_smiles: str | None,
    metadata: dict,
) -> dict:
    canonical_smiles = _canonical_smiles(mol)
    inchikey = _inchikey(mol)
    normalized = _normalize_metadata(metadata)
    identity_value = normalized.get("chemical_identity_id")
    chemical_id = normalized.get("chemical_id")
    if identity_value:
        identity_key = f"chemical_identity_id:{identity_value}"
    elif chemical_id:
        identity_key = f"chemical_id:{chemical_id}"
    elif inchikey:
        identity_key = f"inchikey:{inchikey}"
    elif canonical_smiles:
        identity_key = f"smiles:{canonical_smiles}"
    else:
        identity_key = f"invalid_record:{record_index}"
    return {
        "mol": mol,
        "record_index": record_index,
        "mol_index": record_index,
        "mol_name": mol_name,
        "input_smiles": input_smiles,
        "smiles": canonical_smiles or input_smiles,
        "canonical_smiles": canonical_smiles,
        "inchikey": inchikey,
        "identity_key": identity_key,
        "source": normalized.get("source", ""),
        "metadata": normalized,
        "metadata_errors": list(normalized.get("_metadata_errors") or []),
        **{key: normalized.get(key) for key in LINEAGE_PROPERTIES},
    }


def _sdf_metadata(mol: Any) -> dict:
    if mol is None:
        return {}
    return {str(name): mol.GetProp(name) for name in mol.GetPropNames()}


def _normalize_metadata(metadata: dict) -> dict:
    normalized = dict(metadata or {})
    errors: list[dict] = []
    for key in JSON_PROPERTIES:
        value = normalized.get(key)
        if isinstance(value, str) and value.strip():
            try:
                normalized[key] = json.loads(value)
            except json.JSONDecodeError:
                errors.append({
                    "field": key,
                    "code": "invalid_json",
                    "value": value,
                })
                normalized[key] = None
    for key in INTEGER_PROPERTIES:
        value = normalized.get(key)
        if value in (None, ""):
            normalized[key] = None
            continue
        try:
            normalized[key] = int(value)
        except (TypeError, ValueError):
            errors.append({
                "field": key,
                "code": "invalid_integer",
                "value": value,
            })
            normalized[key] = None
    normalized["_metadata_errors"] = errors
    return normalized


def _canonical_smiles(mol: Any) -> str | None:
    if mol is None:
        return None
    from rdkit import Chem

    try:
        prepared = Chem.RemoveHs(Chem.Mol(mol))
        Chem.SanitizeMol(prepared)
        return Chem.MolToSmiles(
            prepared, canonical=True, isomericSmiles=True
        )
    except Exception:
        return None


def _inchikey(mol: Any) -> str | None:
    if mol is None:
        return None
    from rdkit import Chem

    try:
        return str(Chem.MolToInchiKey(Chem.RemoveHs(Chem.Mol(mol)))) or None
    except Exception:
        return None


def _identity_conflicts(records: list[dict]) -> list[dict]:
    seen: dict[tuple[str, str], dict] = {}
    conflicts: list[dict] = []
    for record in records:
        if record.get("mol") is None:
            continue
        for field in ("chemical_identity_id", "chemical_id"):
            identity_id = record.get(field)
            if not identity_id:
                continue
            graph = record.get("canonical_smiles")
            previous = seen.setdefault((field, str(identity_id)), {
                "canonical_smiles": graph,
                "record_index": record["record_index"],
            })
            if previous["canonical_smiles"] != graph:
                conflicts.append({
                    "identity_field": field,
                    field: identity_id,
                    "first_record_index": previous["record_index"],
                    "first_canonical_smiles": previous["canonical_smiles"],
                    "conflicting_record_index": record["record_index"],
                    "conflicting_canonical_smiles": graph,
                })
    return conflicts


def _integrity_scope(records: list[dict], *, affected: int) -> dict:
    """Decide whether an integrity problem is per-record or systemic.

    A small backtrack affects a single candidate. The same defect affecting a
    sufficiently large fraction and count of records stops the batch for review.
    """
    total = len(records)
    ratio = (affected / total) if total else 0.0
    return {
        "total_records": total,
        "affected_records": affected,
        "affected_ratio": round(ratio, 6),
        "systemic": (
            affected >= SYSTEMIC_INTEGRITY_ABORT_MIN_RECORDS
            and ratio >= SYSTEMIC_INTEGRITY_ABORT_RATIO
        ),
    }


def _apply_identity_conflict_quarantine(
    records: list[dict], conflicts: list[dict]
) -> dict:
    """Quarantine only the records whose declared identity is contradicted.

    A contradicted identity is a structural hard rejection, alongside an invalid
    structure. Only the failing items are rejected. Because we cannot say which
    graph the identity actually names, every record sharing the contradicted
    value is quarantined -- not merely the second one encountered.
    """
    conflicted: dict[str, set[str]] = {}
    for conflict in conflicts:
        field = str(conflict["identity_field"])
        conflicted.setdefault(field, set()).add(str(conflict[field]))
    affected = 0
    for record in records:
        matches = [
            {"identity_field": field, "value": str(record.get(field))}
            for field, values in conflicted.items()
            if record.get(field) and str(record.get(field)) in values
        ]
        if matches:
            record["identity_conflict"] = matches
            affected += 1
    scope = _integrity_scope(records, affected=affected)
    scope["conflicts"] = len(conflicts)
    return scope


def _group_records(records: list[dict], *, deduplicate: bool) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[dict]] = {}
    order: list[tuple[str, str, str]] = []
    for record in records:
        cohort_id = _property_cohort_id(record)
        normalization_pool = _property_normalization_pool(record)
        record["cohort_id"] = cohort_id
        record["normalization_pool"] = normalization_pool
        identity_key = record["identity_key"]
        if not deduplicate:
            identity_key = f"{identity_key}:record:{record['record_index']}"
        group_key = (cohort_id, normalization_pool, identity_key)
        if group_key not in grouped:
            grouped[group_key] = []
            order.append(group_key)
        grouped[group_key].append(record)

    groups = []
    for cohort_id, normalization_pool, identity_key in order:
        members = grouped[(cohort_id, normalization_pool, identity_key)]
        representative = min(members, key=_representative_sort_key)
        digest_input = "|".join((cohort_id, normalization_pool, identity_key))
        candidate_digest = hashlib.sha256(digest_input.encode("utf-8")).hexdigest()
        groups.append({
            "identity_key": identity_key,
            "candidate_id": f"step5:{candidate_digest}",
            "cohort_id": cohort_id,
            "normalization_pool": normalization_pool,
            "representative": representative,
            "records": members,
        })
    return groups


def _representative_sort_key(record: dict) -> tuple:
    rank = record.get("selection_rank")
    rank_key = int(rank) if isinstance(rank, int) else math.inf
    return rank_key, int(record["record_index"])


def _result_identity(group: dict) -> dict:
    representative = group["representative"]
    members = group["records"]
    lineage = [
        {
            "record_index": item["record_index"],
            "mol_name": item["mol_name"],
            **{
                key: item.get(key)
                for key in LINEAGE_PROPERTIES
                if item.get(key) is not None
            },
        }
        for item in members
    ]
    return {
        "candidate_id": group["candidate_id"],
        "identity_key": group["identity_key"],
        "record_index": representative["record_index"],
        "mol_index": representative["mol_index"],
        "mol_name": representative["mol_name"],
        "smiles": representative.get("canonical_smiles") or representative.get("smiles"),
        "input_smiles": representative.get("input_smiles"),
        "inchikey": representative.get("inchikey"),
        "source": representative.get("source", ""),
        "chemical_identity_id": representative.get("chemical_identity_id"),
        "chemical_id": representative.get("chemical_id"),
        "source_identity_id": representative.get("source_identity_id"),
        "source_identity": representative.get("source_identity"),
        "pose_id": representative.get("pose_id"),
        "parent_pose_id": representative.get("parent_pose_id"),
        "coordinate_revision": representative.get("coordinate_revision"),
        "source_pose_index": representative.get("source_pose_index"),
        "selection_rank": representative.get("selection_rank"),
        "lineage_metadata": representative.get("lineage_metadata"),
        "pose_indices": [item["record_index"] for item in members],
        "source_pose_indices": _unique_non_null(
            item.get("source_pose_index") for item in members
        ),
        "pose_ids": _unique_non_null(item.get("pose_id") for item in members),
        "selection_ranks": _unique_non_null(
            item.get("selection_rank") for item in members
        ),
        "lineage": lineage,
        **{
            key: representative.get(key)
            for key in LIGAND_LINEAGE_PROPERTIES
        },
        "cohort_id": group["cohort_id"],
        "normalization_pool": group["normalization_pool"],
        **{
            key: representative.get(key)
            for key in (
                "protein_state_id",
                "pocket_id",
                "round_id",
                "protocol_profile_id",
                "metric_catalog_version",
                "cohort_identity_complete",
                "cohort_identity_missing_fields",
                "cohort_identity_derived_fields",
                "cohort_provenance_sources",
            )
        },
    }


def _property_cohort_id(record: dict) -> str:
    """Return a stable scientific-cohort key for property percentile pools."""
    fields = (
        "protein_state_id",
        "pocket_id",
        "round_id",
        "protocol_profile_id",
        "metric_catalog_version",
    )
    explicit = record.get("declared_cohort_id") or record.get("cohort_id")
    has_components = any(
        record.get(field) not in (None, "") for field in fields
    )
    if not has_components:
        return (
            str(explicit)
            if explicit not in (None, "") else "cohort:unspecified"
        )
    defaults = {
        "protein_state_id": "unspecified",
        "pocket_id": "unspecified",
        "round_id": "unspecified",
        "protocol_profile_id": "unspecified",
        "metric_catalog_version": METRIC_RULES_VERSION,
    }
    canonical_id = "|".join(
        f"{field}={record.get(field) if record.get(field) not in (None, '') else defaults[field]}"
        for field in fields
    )
    if explicit in (None, "") or str(explicit) == canonical_id:
        return canonical_id
    explicit_id = str(explicit)
    prefix = "declared_cohort_id="
    suffix = f"|{canonical_id}"
    if explicit_id.startswith(prefix) and explicit_id.endswith(suffix):
        return explicit_id
    # Treat a declared cohort label as provenance.  The five scientific
    # components remain part of the actual key so a reused label cannot merge
    # incomparable property-normalization pools or dedup groups.
    return f"declared_cohort_id={explicit_id}|{canonical_id}"


def _property_normalization_pool(record: dict) -> str:
    """Keep retrospective DEC_* records out of nomination percentiles."""
    return (
        "benchmark_decoys"
        if benchmark_label(record) == "DEC"
        else "nomination_candidates"
    )


def _run_admet(
    input_data: dict,
    output_dir: Path,
    results: list[dict],
    steps_completed: list[str],
    steps_skipped: list[dict],
) -> dict:
    config = _tool_config(input_data, "admet")
    params = dict(config.get("params") or {})
    candidates = [
        {"candidate_id": result["candidate_id"], "smiles": result["smiles"]}
        for result in results
        if result.get("valid_structure") and result.get("smiles")
    ]
    worker_report = run_admet_worker(
        candidates,
        output_dir / "admet_ai",
        runtime=config.get("runtime") or None,
        timeout_seconds=_integer(params.get("timeout_seconds"), 7200),
    )
    if not worker_report.get("success"):
        error = str(worker_report.get("error") or "ADMET-AI worker failed")
        for result in results:
            if result.get("valid_structure"):
                result["admet"] = {"status": "error", "error": error}
        _declare_tool_failure(results, "admet", error)
        steps_skipped.append({
            "step": "admet",
            "reason": f"{error}; planned but failed, see failure registry",
        })
        return worker_report

    try:
        indexed = _strict_backend_index(candidates, worker_report, "ADMET-AI")
        for candidate_id, prediction in indexed.items():
            _validate_admet_result(prediction, candidate_id=candidate_id)
    except ValueError as exc:
        worker_report["success"] = False
        worker_report["error"] = str(exc)
        for result in results:
            if result.get("valid_structure"):
                result["admet"] = {"status": "error", "error": str(exc)}
        _declare_tool_failure(results, "admet", str(exc))
        steps_skipped.append({
            "step": "admet",
            "reason": f"{exc}; planned but failed, see failure registry",
        })
        return worker_report

    for result in results:
        if not result.get("valid_structure"):
            continue
        prediction = indexed[result["candidate_id"]]
        result["admet"] = prediction
    steps_completed.append("admet")
    return worker_report


def _validate_admet_result(result: dict, *, candidate_id: str) -> None:
    """Fail closed when a nominally complete ADMET-AI record is inconsistent."""
    if str(result.get("status") or "") != "complete":
        return

    def invalid(detail: str) -> ValueError:
        return ValueError(
            "ADMET-AI complete-result contract violation for "
            f"{candidate_id}: {detail}"
        )

    predictions = result.get("predictions")
    issues = validate_native_predictions(predictions)
    if issues:
        raise invalid("; ".join(issues))
    if result.get("endpoint_count") != len(ADMET_NATIVE_ENDPOINT_IDS):
        raise invalid(
            "endpoint_count does not match the registered endpoint payload"
        )
    if result.get("missing_endpoints") not in ([], ()):
        raise invalid("complete result declares missing endpoints")
    if result.get("invalid_endpoints") not in ([], ()):
        raise invalid("complete result declares invalid endpoints")
    if result.get("endpoint_contract_issues") not in ([], ()):
        raise invalid("complete result declares endpoint contract issues")


def _run_chemfh(
    input_data: dict,
    output_dir: Path,
    results: list[dict],
    steps_completed: list[str],
    steps_skipped: list[dict],
) -> dict:
    config = _tool_config(input_data, "chemfh")
    params = dict(config.get("params") or {})
    candidates = [
        {"candidate_id": result["candidate_id"], "smiles": result["smiles"]}
        for result in results
        if result.get("valid_structure") and result.get("smiles")
    ]
    worker_report = run_chemfh_worker(
        candidates,
        output_dir / "chemfh",
        runtime=config.get("runtime") or None,
        timeout_seconds=_integer(params.get("timeout_seconds"), 7200),
        seed=_integer(params.get("seed"), 42),
    )
    if not worker_report.get("success"):
        error = str(worker_report.get("error") or "ChemFH worker failed")
        for result in results:
            if result.get("valid_structure"):
                result["chemfh"] = {"status": "error", "error": error}
        _declare_tool_failure(results, "chemfh", error)
        steps_skipped.append({
            "step": "chemfh",
            "reason": f"{error}; planned but failed, see failure registry",
        })
        return worker_report

    try:
        indexed = _strict_backend_index(candidates, worker_report, "ChemFH")
        for candidate_id, prediction in indexed.items():
            _validate_chemfh_result(prediction, candidate_id=candidate_id)
    except ValueError as exc:
        worker_report["success"] = False
        worker_report["error"] = str(exc)
        for result in results:
            if result.get("valid_structure"):
                result["chemfh"] = {"status": "error", "error": str(exc)}
        _declare_tool_failure(results, "chemfh", str(exc))
        steps_skipped.append({
            "step": "chemfh",
            "reason": f"{exc}; planned but failed, see failure registry",
        })
        return worker_report

    for result in results:
        if result.get("valid_structure"):
            result["chemfh"] = indexed[result["candidate_id"]]
    steps_completed.append("chemfh")
    return worker_report


def _declare_tool_failure(results: list[dict], tool: str, reason: str) -> None:
    """Record that a tool ran and failed, on every candidate it should have covered.

    The executable policy separates the two ways a metric can be blank: a tool nobody
    selected leaves Σ_planned, while a tool that was planned and failed stays in
    it and goes to the failure registry.  Only the second one can be told from
    the outside, and only if it is declared -- with no row anywhere holding a
    value, the scorer cannot distinguish a total outage from a metric this build
    never planned.  Declaring it here is what keeps the gap visible as a
    coverage shortfall instead of silently reading as full coverage.
    """
    for result in results:
        if not result.get("valid_structure"):
            continue
        context = result.setdefault("metric_context", {})
        errors = context.setdefault("operational_errors", [])
        if not any(
            isinstance(entry, dict) and entry.get("tool") == tool
            for entry in errors
        ):
            errors.append({"tool": tool, "reason": reason})


def _mark_synthesis_operational_failure(result: dict, reason: str) -> None:
    """Block nomination for an AiZynthFinder outage without judging chemistry."""
    reason = str(reason or "AiZynthFinder failed for this candidate")
    _declare_tool_failure([result], "aizynthfinder", reason)

    # This is a protocol/data gap, not evidence that the molecule needs chemical
    # modification. Keep it out of ``molecule_review_reasons`` so the legacy
    # classification reducer preserves the pre-synthesis molecular gate.
    result["nomination_eligibility"] = "review_only"
    molecule_reasons = result.setdefault("molecule_review_reasons", [])
    result["molecule_review_reasons"] = [
        item for item in molecule_reasons if item != "synthesis_route_unsolved"
    ]
    protocol_reasons = result.setdefault("protocol_review_reasons", [])
    if "synthesis_error" not in protocol_reasons:
        protocol_reasons.append("synthesis_error")
    eligibility_reasons = result.setdefault(
        "nomination_eligibility_reasons", [],
    )
    eligibility_reasons[:] = [
        item for item in eligibility_reasons
        if item != "synthesis_route_unsolved"
    ]
    if "synthesis_error" not in eligibility_reasons:
        eligibility_reasons.append("synthesis_error")

    # ``property_protocol_evidence`` is built before retrosynthesis. Amend the
    # existing object instead of leaving it claiming that automatic nomination
    # is safe after a selected AiZynthFinder call failed.
    evidence = result.get("property_protocol_evidence")
    if isinstance(evidence, dict):
        evidence["status"] = "review_required"
        evidence["blocks_automatic_nomination"] = True
        evidence_reasons = evidence.setdefault("reasons", [])
        if "synthesis_error" not in evidence_reasons:
            evidence_reasons.append("synthesis_error")
        evidence_actions = evidence.setdefault("actions", [])
        if SYNTHESIS_RETRY_ACTION not in evidence_actions:
            evidence_actions.append(SYNTHESIS_RETRY_ACTION)
        evidence["synthesis"] = {
            "planned": True,
            "attempted": True,
            "backend_status": "error",
            "evidence_status": "review_required",
            "review_required": True,
        }


def _mark_synthesis_not_attempted(
    result: dict,
    status: str,
    *,
    reason: str | None = None,
) -> None:
    """Publish an explicit N/A state for rows outside this round's search."""
    synthesis = {"status": status, "attempted": False}
    if reason:
        synthesis["reason"] = reason
    result["synthesis"] = synthesis

    # AiZynthFinder can be enabled for the stage while an individual row is
    # intentionally outside the top-K/eligibility slice.  Mark those scalar
    # metrics N/A per row so score_records emits ``not_planned`` rather than
    # ``planned_missing`` (the latter means an attempted tool failed).
    context = result.setdefault("metric_context", {})
    declared = context.get("not_applicable_metric_ids")
    not_applicable = (
        list(declared)
        if isinstance(declared, (list, tuple, set))
        else []
    )
    context["not_applicable_metric_ids"] = list(dict.fromkeys([
        *not_applicable,
        *AIZYNTHFINDER_METRIC_IDS,
    ]))


def _mark_synthesis_selected(result: dict, position: int) -> None:
    """Select a row for execution and clear any stale per-row N/A marker."""
    result["synthesis"] = {
        "status": "selected",
        "attempted": False,
        "top_k_rank": position,
    }
    context = result.get("metric_context")
    if not isinstance(context, dict):
        return
    declared = context.get("not_applicable_metric_ids")
    if not isinstance(declared, (list, tuple, set)):
        return
    synthesis_ids = set(AIZYNTHFINDER_METRIC_IDS)
    context["not_applicable_metric_ids"] = [
        metric_id for metric_id in declared if metric_id not in synthesis_ids
    ]


def _failed_synthesis_result(result: dict, error: str) -> dict:
    """Return the synthesis-report row for a selected batch-level failure."""
    return {
        "candidate_id": result["candidate_id"],
        "smiles": result["smiles"],
        "status": "error",
        "attempted": True,
        "is_solved": False,
        "route_count": 0,
        "statistics": None,
        "metrics": _empty_synthesis_metrics(),
        "error": error,
        "routes": [],
    }


def _report_operational_errors(results: list[dict]) -> list[dict]:
    """Aggregate row-local tool failures for the stage-level failure registry."""
    grouped: dict[tuple[str, str], list[str]] = {}
    for result in results:
        context = result.get("metric_context")
        if not isinstance(context, dict):
            continue
        row_id = str(
            result.get("chemical_identity_id")
            or result.get("chemical_id")
            or result.get("pose_id")
            or result.get("candidate_id")
            or result.get("mol_name")
            or ""
        )
        for entry in context.get("operational_errors") or []:
            if isinstance(entry, dict):
                tool = str(entry.get("tool") or "unknown")
                reason = str(entry.get("reason") or "tool execution failed")
            else:
                tool = str(entry or "unknown")
                reason = "tool execution failed"
            bucket = grouped.setdefault((tool, reason), [])
            if row_id and row_id not in bucket:
                bucket.append(row_id)
    return [
        {"tool": tool, "reason": reason, "row_ids": row_ids}
        for (tool, reason), row_ids in grouped.items()
    ]


def _update_property_backend_eligibility(
    results: list[dict],
    *,
    rules: tuple[dict, ...] = PROPERTY_FILTER_METRIC_RULES,
) -> None:
    """Reduce molecular findings and protocol evidence without conflating them.

    A required backend outage is not evidence that a molecule is bad, but it is
    still an evidence gap that must prevent automatic nomination.  The two
    reason lists below deliberately remain separate so downstream consumers can
    keep the molecule's three-tier classification while routing a protocol
    repair.  Tools intentionally omitted from ``enabled_tools`` are not gaps.
    """
    for result in results:
        if not result.get("valid_structure"):
            # Unreadable chemistry is structural, not a threshold breach, so it
            # is the only condition here that rejects outright.  A ``fail``
            # verdict from the RDKit screen is a level-1 mark and falls through
            # to ``molecule_reasons`` below, where it caps nomination without
            # deciding the molecule.
            result["structural_reject"] = True
            result.setdefault("structural_reject_reason", "invalid_structure")
            result["nomination_eligibility"] = "rejected"
            continue

        if result.get("identity_reject") or result.get("structural_reject"):
            # The other half of the executable policy: readable chemistry whose *identity*
            # is contradicted.  Already rejected upstream; recomputing from the
            # (absent) metric evidence would promote it back to review_only.
            result["nomination_eligibility"] = "rejected"
            continue

        molecule_reasons: list[str] = []
        tool_error = result.get("property_tool_error")
        _, rdkit_gaps = _rdkit_evidence_split(result.get("rdkit_evidence"))
        if (
            result.get("verdict") != "pass"
            # A non-pass RDKit verdict is only a molecular finding when the
            # screen actually marked something about *this* molecule.  Without
            # this test a missing sascorer, an unavailable FilterCatalog or a
            # fully disabled RDKit pass all appended the same
            # ``rdkit_property_review``, which caps nomination at review_only and
            # drives the legacy gate to MODIFY -- penalising the molecule for a
            # tool outage, which the executable policy forbids.  The gaps still block
            # automatic nomination, but through ``protocol_reasons`` below.
            and _rdkit_verdict_is_molecular(result)
        ):
            molecule_reasons.append("rdkit_property_review")

        chemfh = result.get("chemfh") or {}
        if str(chemfh.get("status") or "missing") == "complete":
            positive_state = _chemfh_positive_state(chemfh)
            if positive_state is True:
                # A real frequent-hitter prediction for this molecule.
                molecule_reasons.append("chemfh_positive_prediction")
            elif positive_state is None:
                # The run completed for this candidate yet its own verdict is
                # unreadable -- a per-molecule gap, not a tool outage.
                molecule_reasons.append("chemfh_positive_state_missing")

        protocol_evidence = _property_protocol_evidence(result, rules=rules)
        protocol_reasons = list(protocol_evidence["reasons"])
        protocol_reasons.extend(f"rdkit_{gap}" for gap in rdkit_gaps)
        if isinstance(tool_error, dict):
            # Our own RDKit pass raised.  That is a protocol repair, never a
            # molecular finding, so it blocks nomination without
            # touching the molecule's three-tier classification.
            protocol_reasons.append("rdkit_property_pass_failed")
        if result.get("metadata_integrity_errors"):
            # Same shape of problem one layer up: the chemistry is fine but its
            # lineage/ranking provenance could not be parsed, so cohort keys and
            # percentile pools for this row are suspect.
            protocol_reasons.append("metadata_integrity_error")
        protocol_reasons = list(dict.fromkeys(protocol_reasons))
        reasons = list(dict.fromkeys((*molecule_reasons, *protocol_reasons)))

        result["nomination_eligibility"] = (
            "review_only" if reasons else "auto_candidate"
        )
        result["molecule_review_reasons"] = molecule_reasons
        result["protocol_review_reasons"] = protocol_reasons
        result["nomination_eligibility_reasons"] = reasons
        result["property_protocol_evidence"] = protocol_evidence


def _rdkit_evidence_split(evidence: object) -> tuple[list[str] | None, list[str]]:
    """Split an RDKit screen result into molecular marks and coverage gaps.

    Returns ``(None, [])`` for rows written before ``compute_verdict`` reported
    the split -- reloaded reports, hand-built fixtures, the pre-RDKit rejection
    branches.  ``None`` means "unknown", and the caller keeps the older, more
    conservative reading of a non-pass verdict for those rows.
    """
    if not isinstance(evidence, dict):
        return None, []
    findings = evidence.get("molecular_findings")
    gaps = evidence.get("coverage_gaps")
    return (
        [str(item) for item in findings] if isinstance(findings, list) else None,
        [str(item) for item in gaps] if isinstance(gaps, list) else [],
    )


def _rdkit_verdict_is_molecular(result: dict) -> bool:
    """Did the RDKit screen mark this molecule, or only fail to look at it?

    ``True`` also when the answer is unknown: rows written before
    ``compute_verdict`` reported the split keep their previous, more
    conservative reading.
    """
    if isinstance(result.get("property_tool_error"), dict):
        return False
    findings, _ = _rdkit_evidence_split(result.get("rdkit_evidence"))
    return True if findings is None else bool(findings)


def _property_protocol_evidence(
    result: dict,
    *,
    rules: tuple[dict, ...] = PROPERTY_FILTER_METRIC_RULES,
) -> dict:
    """Describe whether planned ADMET/alert evidence permits nomination.

    ADMET ``partial`` is endpoint-local: valid endpoints continue to be marked
    and scored.  It is sufficient only when the explicit ABS/TOX planned-weighted
    floors are met and any metric explicitly declared required by the scoring
    plan is present.  The named toxicity hard risks remain non-compensable when
    measured, but the spec does not make all four mandatory by default.  Most
    ChemFH endpoints carry little or no score weight (only reactive compounds at
    ``1`` and promiscuous/colloidal at ``0.5``), so its outage would barely move
    a weighted total and needs an explicit protocol check instead; an ordinary
    RDKit FilterCatalog output is not silently treated as equivalent to ChemFH's
    fluorescence/aggregation models, and a replacement must be explicitly
    versioned and carry coverage and confidence.
    """
    reasons: list[str] = []
    actions: list[str] = []

    admet = result.get("admet") if isinstance(result.get("admet"), dict) else {}
    admet_status = str(admet.get("status") or "missing").strip().casefold()
    admet_planned = _property_tool_planned(result, "admet", admet_status)
    admet_required = _property_tool_required(result, "admet", admet_status)
    category_coverage = _admet_planned_weighted_coverage(result, rules=rules)
    required_admet_metrics = _required_admet_metric_ids(result, rules)
    not_applicable_metrics = _declared_property_metric_ids(
        _property_plan_contexts(result), "not_applicable_metric_ids",
    )
    admet_plan_conflicts = sorted(
        set(required_admet_metrics).intersection(not_applicable_metrics)
    )
    missing_required = [
        metric_id
        for metric_id in required_admet_metrics
        if _metric_value_by_id(result, rules, metric_id) is None
    ]
    if admet_plan_conflicts:
        reasons.append("admet_required_not_applicable_plan_conflict")
        actions.append(
            "Resolve the scoring plan conflict: the same ADMET endpoint "
            "cannot be both required and not_applicable"
        )
    accepted_partial = False
    admet_review_required = False
    if not admet_planned:
        admet_evidence_status = "not_planned"
        category_coverage = {
            category: None for category in ADMET_PARTIAL_COVERAGE_CATEGORIES
        }
        missing_required = []
    elif not admet_required:
        # Optional execution still contributes its observed values and missing
        # coverage.  Its failure does not create an independent nomination veto.
        accepted_partial = admet_status == "partial" and all(
            category_coverage.get(category) is None
            or category_coverage[category] >= CATEGORY_PLANNED_COVERAGE_FLOOR[category]
            for category in ADMET_PARTIAL_COVERAGE_CATEGORIES
        )
        admet_evidence_status = (
            "complete" if admet_status == "complete"
            else "sufficient_partial" if accepted_partial
            else "optional_incomplete"
        )
    elif admet_status == "complete":
        admet_review_required = bool(missing_required)
        admet_evidence_status = (
            "review_required" if missing_required else "complete"
        )
        if missing_required:
            reasons.append("admet_missing_required_metrics")
            actions.append(
                "Repair the declared complete ADMET payload so every "
                "explicitly required endpoint has valid evidence"
            )
    elif admet_status == "partial":
        low_categories = [
            category
            for category in ADMET_PARTIAL_COVERAGE_CATEGORIES
            if category_coverage.get(category) is not None
            and category_coverage[category]
            < CATEGORY_PLANNED_COVERAGE_FLOOR[category]
        ]
        for category in low_categories:
            reasons.append(f"admet_partial_low_{category.lower()}_coverage")
        if missing_required:
            reasons.append("admet_missing_required_metrics")
        admet_review_required = bool(low_categories or missing_required)
        accepted_partial = not admet_review_required
        admet_evidence_status = (
            "review_required" if admet_review_required else "sufficient_partial"
        )
        if admet_review_required:
            actions.append(
                "Rerun or replace ADMET-AI to restore the ABS/TOX planned "
                "weighted coverage floors and explicitly required endpoint "
                "evidence before automatic nomination"
            )
    else:
        admet_review_required = True
        admet_evidence_status = "review_required"
        reasons.append(f"admet_{admet_status}")
        actions.append(
            "Rerun or replace the planned ADMET-AI backend and capture usable "
            "endpoint evidence before automatic nomination"
        )
    if admet_plan_conflicts:
        admet_review_required = True
        accepted_partial = False
        admet_evidence_status = "review_required"

    chemfh = (
        result.get("chemfh") if isinstance(result.get("chemfh"), dict) else {}
    )
    chemfh_status = str(chemfh.get("status") or "missing").strip().casefold()
    chemfh_planned = _property_tool_planned(result, "chemfh", chemfh_status)
    chemfh_required = _property_tool_required(result, "chemfh", chemfh_status)
    required_chemfh_metrics = _required_chemfh_metric_ids(result)
    if chemfh_required and not required_chemfh_metrics:
        required_chemfh_metrics = sorted(_property_tool_metric_ids("chemfh"))
    chemfh_plan_conflicts = sorted(
        set(required_chemfh_metrics).intersection(not_applicable_metrics)
    )
    rdkit_fallback = _rdkit_filtercatalog_fallback(
        result, required_metric_ids=required_chemfh_metrics,
    )
    fallback_satisfied = bool(
        chemfh_required
        and chemfh_status != "complete"
        and not chemfh_plan_conflicts
        and rdkit_fallback["status"] == "complete"
    )
    chemfh_review_required = bool(
        chemfh_plan_conflicts
        or (
            chemfh_required
            and chemfh_status != "complete"
            and not fallback_satisfied
        )
    )
    if chemfh_plan_conflicts:
        reasons.append("chemfh_required_not_applicable_plan_conflict")
        actions.append(
            "Resolve the scoring plan conflict: the same ChemFH metric "
            "cannot be both required and not_applicable"
        )
    if not chemfh_planned:
        chemfh_evidence_status = "not_planned"
    elif chemfh_status == "complete":
        chemfh_evidence_status = "complete"
    elif fallback_satisfied:
        chemfh_evidence_status = "fallback_complete"
    else:
        chemfh_evidence_status = (
            "review_required" if chemfh_required else "optional_incomplete"
        )
        if chemfh_required:
            reasons.append(f"chemfh_{chemfh_status}")
            actions.append(
                "Rerun or replace the required ChemFH backend, or capture a "
                "versioned full-coverage equivalent fallback before automatic "
                "nomination"
            )

    contexts = _property_plan_contexts(result)
    rascore_required = (
        "rascore" in _declared_property_tools(contexts, "required_tools")
        or "rascore.score" in _declared_property_metric_ids(
            contexts, "required_metric_ids",
        )
    )
    rascore = result.get("rascore") or {}
    rascore_review_required = rascore_required and (
        rascore.get("status") != "available"
        or _finite_number_or_none(rascore.get("score")) is None
    )
    if rascore_review_required:
        reasons.append("rascore_prediction_unavailable")
        actions.append("Repair or rerun the explicitly required RAscore worker")
    blocks_nomination = (
        admet_review_required or chemfh_review_required or rascore_review_required
    )
    return {
        "status": "review_required" if blocks_nomination else "sufficient",
        "blocks_automatic_nomination": blocks_nomination,
        "reasons": list(dict.fromkeys(reasons)),
        "actions": list(dict.fromkeys(actions)),
        "missing_required_metrics": missing_required,
        "admet": {
            "planned": admet_planned,
            "required": admet_required,
            "backend_status": admet_status,
            "evidence_status": admet_evidence_status,
            "accepted_partial": accepted_partial,
            "category_planned_weighted_coverage": category_coverage,
            "coverage_floors": {
                category: CATEGORY_PLANNED_COVERAGE_FLOOR[category]
                for category in ADMET_PARTIAL_COVERAGE_CATEGORIES
            },
            "required_metric_ids": required_admet_metrics,
            "missing_required_metrics": missing_required,
            "plan_conflicts": admet_plan_conflicts,
            "review_required": admet_review_required,
        },
        "chemfh": {
            "planned": chemfh_planned,
            "required": chemfh_required,
            "backend_status": chemfh_status,
            "evidence_status": chemfh_evidence_status,
            "required_metric_ids": required_chemfh_metrics,
            "plan_conflicts": chemfh_plan_conflicts,
            "fallback_satisfied": fallback_satisfied,
            "fallback": rdkit_fallback,
            "review_required": chemfh_review_required,
        },
    }


def _property_tool_planned(result: dict, tool: str, status: str) -> bool:
    """Use explicit planned/required/N/A declarations before legacy status."""
    tool = str(tool).strip().casefold()
    contexts = _property_plan_contexts(result)
    known_metric_ids = _property_tool_metric_ids(tool)
    not_applicable = _declared_property_metric_ids(
        contexts, "not_applicable_metric_ids",
    )
    required_ids = _declared_property_metric_ids(
        contexts, "required_metric_ids", "required_admet_metric_ids",
    )
    planned_ids = _declared_property_metric_ids(
        contexts,
        "planned_metric_ids",
        "failed_metric_ids",
        "planned_but_failed_metric_ids",
    ) | required_ids

    if required_ids.intersection(known_metric_ids):
        return True
    if (planned_ids - not_applicable).intersection(known_metric_ids):
        return True
    for field in (
        "planned_tools", "required_tools", "optional_tools",
        "failed_tools", "planned_but_failed_tools",
    ):
        if tool in _declared_property_tools(contexts, field):
            return True
    if any(
        isinstance(context.get("operational_errors"), (list, tuple))
        and any(
            (
                str(item.get("tool") or "").strip().casefold()
                if isinstance(item, dict) else str(item).strip().casefold()
            ) == tool
            for item in context["operational_errors"]
        )
        for context in contexts
    ):
        return True

    explicit_metric_plan = any(
        "planned_metric_ids" in context for context in contexts
    )
    if explicit_metric_plan:
        # An explicitly empty/subset plan is authoritative.  In particular,
        # N/A metrics do not re-enter Σ_planned merely because a backend switch
        # exists in the same process configuration.
        return False

    if known_metric_ids and known_metric_ids.issubset(not_applicable):
        return False

    enabled_tools = _declared_property_tools(contexts, "enabled_tools")
    if enabled_tools:
        return tool in enabled_tools
    return status not in {"", "disabled", "missing", "not_applicable"}


def _property_tool_required(result: dict, tool: str, status: str) -> bool:
    """Separate required evidence from an optional backend execution."""
    tool = str(tool).strip().casefold()
    contexts = _property_plan_contexts(result)
    known_metric_ids = _property_tool_metric_ids(tool)
    required_ids = _declared_property_metric_ids(
        contexts, "required_metric_ids", "required_admet_metric_ids",
    )
    if required_ids.intersection(known_metric_ids):
        return True
    if tool in _declared_property_tools(contexts, "required_tools"):
        return True
    explicit_required_contract = any(
        "required_tools" in context
        or "optional_tools" in context
        or "required_metric_ids" in context
        or "required_admet_metric_ids" in context
        for context in contexts
    )
    if explicit_required_contract:
        return False
    # Backward compatibility is deliberately fail-closed: legacy runs exposed
    # only enabled/planned, not the required/optional distinction.
    return _property_tool_planned(result, tool, status)


def _property_plan_contexts(result: dict) -> list[dict]:
    return [
        context
        for context in (result.get("scoring_plan"), result.get("metric_context"))
        if isinstance(context, dict)
    ]


def _normalize_property_metric_id(value: object) -> str:
    metric_id = str(value).strip()
    if not metric_id:
        return ""
    if "." in metric_id:
        return metric_id
    if metric_id in ADMET_NATIVE_ENDPOINT_IDS:
        return f"admet_ai.{metric_id}"
    if metric_id in CHEMFH_ENDPOINT_IDS:
        return f"chemfh.{metric_id}"
    return metric_id


def _declared_property_metric_ids(
    contexts: list[dict], *fields: str,
) -> set[str]:
    declared: set[str] = set()
    for context in contexts:
        for field in fields:
            values = context.get(field)
            if not isinstance(values, (list, tuple, set)):
                continue
            declared.update(
                metric_id
                for metric_id in (
                    _normalize_property_metric_id(value) for value in values
                )
                if metric_id
            )
    return declared


def _declared_property_tools(contexts: list[dict], field: str) -> set[str]:
    declared: set[str] = set()
    for context in contexts:
        values = context.get(field)
        if isinstance(values, (list, tuple, set)):
            declared.update(str(item).strip().casefold() for item in values)
    declared.discard("")
    return declared


def _property_tool_metric_ids(tool: str) -> set[str]:
    prefix = "admet_ai." if tool == "admet" else f"{tool}."
    return {
        str(rule.get("metric_id"))
        for rule in PROPERTY_FILTER_METRIC_RULES
        if str(rule.get("metric_id") or "").startswith(prefix)
    }


def _admet_planned_weighted_coverage(
    result: dict,
    *,
    rules: tuple[dict, ...],
) -> dict[str, float | None]:
    """Compute the sbdd2 planned |weight| coverage for ABS and TOX."""
    contexts = _property_plan_contexts(result)
    planned_ids = _declared_property_metric_ids(
        contexts,
        "planned_metric_ids",
        "required_metric_ids",
        "required_admet_metric_ids",
        "failed_metric_ids",
        "planned_but_failed_metric_ids",
    )
    not_applicable = _declared_property_metric_ids(
        contexts, "not_applicable_metric_ids",
    )
    explicit_metric_plan = any(
        "planned_metric_ids" in context for context in contexts
    )
    coverage: dict[str, float | None] = {}
    for category in ADMET_PARTIAL_COVERAGE_CATEGORIES:
        category_rules = [
            rule
            for rule in rules
            if rule.get("category_code") == category
            and rule.get("requires_property_tool") == "admet"
            and abs(float(rule.get("weight", 0.0))) > 0
            and str(rule.get("metric_id")) not in not_applicable
            and (
                not explicit_metric_plan
                or str(rule.get("metric_id")) in planned_ids
            )
        ]
        planned_weight = sum(
            abs(float(rule.get("weight", 0.0))) for rule in category_rules
        )
        available_weight = sum(
            abs(float(rule.get("weight", 0.0)))
            for rule in category_rules
            if extract_metric_value(result, rule) is not None
        )
        coverage[category] = (
            round(available_weight / planned_weight, 6)
            if planned_weight > 0
            else None
        )
    return coverage


def _metric_value_by_id(
    result: dict,
    rules: tuple[dict, ...],
    metric_id: str,
) -> float | None:
    for rule in rules:
        if rule.get("metric_id") == metric_id:
            return extract_metric_value(result, rule)
    return None


def _required_admet_metric_ids(
    result: dict,
    rules: tuple[dict, ...],
) -> list[str]:
    """Read only explicitly declared required ADMET metrics from the plan."""
    known = {
        str(rule.get("metric_id"))
        for rule in rules
        if rule.get("requires_property_tool") == "admet"
    }
    declared: list[str] = []
    for context in (result.get("scoring_plan"), result.get("metric_context")):
        if not isinstance(context, dict):
            continue
        for field in ("required_metric_ids", "required_admet_metric_ids"):
            values = context.get(field)
            if not isinstance(values, (list, tuple, set)):
                continue
            for value in values:
                metric_id = str(value).strip()
                if metric_id and not metric_id.startswith("admet_ai."):
                    metric_id = f"admet_ai.{metric_id}"
                if metric_id in known and metric_id not in declared:
                    declared.append(metric_id)
    return declared


def _required_chemfh_metric_ids(result: dict) -> list[str]:
    known = _property_tool_metric_ids("chemfh")
    declared = _declared_property_metric_ids(
        _property_plan_contexts(result), "required_metric_ids",
    )
    return sorted(declared.intersection(known))


def _rdkit_filtercatalog_fallback(
    result: dict,
    *,
    required_metric_ids: list[str] | None = None,
) -> dict:
    """Return an auditable success state for an explicitly equivalent fallback.

    A normal PAINS/Brenk/NIH pass is not evidence for fluorescence, FLuc, or
    colloidal aggregation.  Automatic nomination is unblocked only when the
    replacement event pins its artifact/catalog mapping and explicitly covers
    every required ChemFH metric with reported coverage and confidence.
    """
    context = result.get("metric_context")
    enabled_tools = context.get("enabled_tools") if isinstance(context, dict) else None
    pains_planned = (
        "pains" in {str(item).strip().casefold() for item in enabled_tools}
        if isinstance(enabled_tools, (list, tuple, set))
        else False
    )
    alerts = (
        result.get("structural_alerts")
        if isinstance(result.get("structural_alerts"), dict)
        else {}
    )
    catalogs = {
        str(item).strip().upper()
        for item in (alerts.get("catalogs") or [])
        if str(item).strip()
    }
    expected_catalogs = {str(item).upper() for item in DEFAULT_ALERT_CATALOGS}
    catalog_errors = alerts.get("catalog_errors") or []
    fallback_event = _explicit_chemfh_fallback_event(result)
    required = {
        _normalize_property_metric_id(item)
        for item in (required_metric_ids or [])
        if _normalize_property_metric_id(item)
    }
    covered = {
        _normalize_property_metric_id(item)
        for item in (
            (fallback_event or {}).get("covered_metric_ids") or []
        )
        if _normalize_property_metric_id(item)
    }
    coverage = _finite_number_or_none(
        (fallback_event or {}).get("coverage")
    )
    confidence = _finite_number_or_none(
        (fallback_event or {}).get("confidence")
    )
    artifact_sha256 = str(
        (fallback_event or {}).get("artifact_sha256")
        or (fallback_event or {}).get("catalog_sha256")
        or ""
    ).strip().casefold()
    mapping_version = str(
        (fallback_event or {}).get("mapping_version") or ""
    ).strip()
    tool_version = str(
        (fallback_event or {}).get("tool_version")
        or (fallback_event or {}).get("catalog_version")
        or ""
    ).strip()
    attestation_errors: list[str] = []
    if fallback_event is not None:
        if not required:
            attestation_errors.append("required_metric_ids_empty")
        if not required.issubset(covered):
            attestation_errors.append("required_metric_coverage_incomplete")
        if coverage is None or not 0 <= coverage <= 1 or coverage < 1:
            attestation_errors.append("coverage_must_be_1.0")
        if confidence is None or not 0 < confidence <= 1:
            attestation_errors.append("confidence_must_be_in_(0,1]")
        if not tool_version:
            attestation_errors.append("tool_or_catalog_version_missing")
        if not mapping_version:
            attestation_errors.append("mapping_version_missing")
        if len(artifact_sha256) != 64 or any(
            character not in "0123456789abcdef"
            for character in artifact_sha256
        ):
            attestation_errors.append("artifact_sha256_invalid")
    complete = bool(
        fallback_event is not None
        and not attestation_errors
        and pains_planned
        and alerts.get("available") is True
        and isinstance(alerts.get("has_alert"), bool)
        and isinstance(alerts.get("matches"), list)
        and not catalog_errors
        and expected_catalogs.issubset(catalogs)
    )
    return {
        "tool": "rdkit_filtercatalog",
        "status": (
            "complete"
            if complete
            else "not_declared"
            if fallback_event is None
            else "incomplete"
        ),
        "explicit_fallback_event": fallback_event,
        "catalogs": sorted(catalogs),
        "required_catalogs": sorted(expected_catalogs),
        "catalog_errors": catalog_errors,
        "required_metric_ids": sorted(required),
        "covered_metric_ids": sorted(covered),
        "coverage": coverage,
        "confidence": confidence,
        "tool_version": tool_version or None,
        "mapping_version": mapping_version or None,
        "artifact_sha256": artifact_sha256 or None,
        "attestation_errors": attestation_errors,
    }


def _explicit_chemfh_fallback_event(result: dict) -> dict | None:
    """Read an explicitly approved ChemFH replacement event.

    RDKit alert screening normally runs in parallel with ChemFH.  Its mere
    success cannot silently rescue a ChemFH identity/schema failure: the run
    must record that this exact successful check was intentionally approved as
    an equivalent replacement.
    """
    for context in _property_plan_contexts(result):
        events = context.get("equivalent_fallbacks")
        if not isinstance(events, (list, tuple)):
            continue
        for event in events:
            if not isinstance(event, dict):
                continue
            replaces = str(
                event.get("replaces") or event.get("for_tool") or ""
            ).strip().casefold()
            tool = "".join(
                character
                for character in str(event.get("tool") or "").casefold()
                if character.isalnum()
            )
            status = str(event.get("status") or "").strip().casefold()
            if (
                replaces == "chemfh"
                and tool == "rdkitfiltercatalog"
                and status in {"complete", "success"}
                and event.get("equivalent") is True
                and event.get("approved") is True
            ):
                return dict(event)
    return None


def _validate_chemfh_result(result: dict, *, candidate_id: str) -> None:
    """Validate complete *and partial* ChemFH evidence endpoint by endpoint.

    ``partial`` is a payload completeness state, not permission to skip schema
    and range checks.  Every present endpoint remains usable only when its
    probability, uncertainty, and boolean call agree; listed missing endpoints
    remain NA and therefore cannot enter a numerical gate.
    """
    status = str(result.get("status") or "").strip().casefold()
    if status not in {"complete", "partial"}:
        raise ValueError(
            f"ChemFH result contract violation for {candidate_id}: "
            f"unsupported status {status or '<missing>'!r}"
        )

    def invalid(detail: str) -> ValueError:
        return ValueError(
            f"ChemFH {status}-result contract violation for "
            f"{candidate_id}: {detail}"
        )

    predictions = result.get("predictions")
    if not isinstance(predictions, dict) or set(predictions) != set(
        CHEMFH_ENDPOINT_IDS
    ):
        raise invalid(
            "expected exactly the seven registered model endpoints"
        )
    if result.get("endpoint_count") != len(CHEMFH_ENDPOINT_IDS):
        raise invalid("endpoint_count does not match the endpoint payload")
    declared_missing_endpoints = result.get("missing_endpoints")
    if not isinstance(declared_missing_endpoints, (list, tuple)):
        raise invalid("missing_endpoints must be a list")
    missing_endpoints = {str(item) for item in declared_missing_endpoints}
    if not missing_endpoints.issubset(CHEMFH_ENDPOINT_IDS):
        raise invalid("missing_endpoints contains an unregistered endpoint")

    positive_endpoints: list[str] = []
    observed_missing_endpoints: set[str] = set()
    for endpoint_id in CHEMFH_ENDPOINT_IDS:
        record = predictions.get(endpoint_id)
        if not isinstance(record, dict):
            raise invalid(f"endpoint {endpoint_id!r} is not an object")
        probability = _finite_number_or_none(record.get("probability"))
        uncertainty = _finite_number_or_none(record.get("uncertainty"))
        if probability is not None and not 0.0 <= probability <= 1.0:
            raise invalid(
                f"endpoint {endpoint_id!r} probability is outside [0, 1]"
            )
        if uncertainty is not None and uncertainty < 0.0:
            raise invalid(f"endpoint {endpoint_id!r} uncertainty is negative")
        endpoint_missing = probability is None or uncertainty is None
        if endpoint_missing:
            observed_missing_endpoints.add(endpoint_id)
        positive = record.get("positive")
        if endpoint_missing:
            if positive is not None:
                raise invalid(
                    f"missing endpoint {endpoint_id!r} must not publish a "
                    "positive call"
                )
            if record.get("high_confidence") is not None:
                raise invalid(
                    f"missing endpoint {endpoint_id!r} must not publish a "
                    "confidence call"
                )
            continue
        if not isinstance(positive, bool):
            raise invalid(f"endpoint {endpoint_id!r} has no boolean positive call")
        if positive != (probability > CHEMFH_POSITIVE_THRESHOLD):
            raise invalid(
                f"endpoint {endpoint_id!r} positive call disagrees with probability"
            )
        if positive:
            positive_endpoints.append(endpoint_id)
    if missing_endpoints != observed_missing_endpoints:
        raise invalid(
            "missing_endpoints does not match endpoints lacking complete "
            "probability/uncertainty evidence"
        )
    if status == "complete" and missing_endpoints:
        raise invalid("complete result declares missing endpoints")

    rules = result.get("rules")
    if not isinstance(rules, dict) or set(rules) != set(CHEMFH_RULE_NAMES):
        raise invalid("expected exactly the 17 registered native rules")
    if result.get("rule_count") != len(CHEMFH_RULE_NAMES):
        raise invalid("rule_count does not match the rule payload")
    declared_missing_rules = result.get("missing_rules")
    if not isinstance(declared_missing_rules, (list, tuple)):
        raise invalid("missing_rules must be a list")
    missing_rules = {str(item) for item in declared_missing_rules}
    if not missing_rules.issubset(CHEMFH_RULE_NAMES):
        raise invalid("missing_rules contains an unregistered native rule")
    if status == "complete" and missing_rules:
        raise invalid("complete result declares missing rules")
    for rule_name in CHEMFH_RULE_NAMES:
        rule = rules.get(rule_name)
        if not isinstance(rule, dict):
            raise invalid(f"rule {rule_name!r} is not an object")
        match_count = rule.get("match_count")
        atom_indices = rule.get("atom_indices")
        if (
            isinstance(match_count, bool)
            or not isinstance(match_count, int)
            or match_count < 0
            or not isinstance(atom_indices, list)
            or len(atom_indices) != match_count
        ):
            raise invalid(f"rule {rule_name!r} has inconsistent match evidence")
        if rule.get("decision_contribution") is not False:
            raise invalid(f"rule {rule_name!r} must remain evidence-only")

    pains = rules["PAINS"]
    if (
        pains.get("overlaps_with") != "pains"
        or pains.get("duplicate_decision_suppressed") is not True
    ):
        raise invalid("PAINS duplicate-decision suppression is missing")

    summary = result.get("summary")
    if not isinstance(summary, dict):
        raise invalid("summary is missing")
    if summary.get("has_positive_prediction") != bool(positive_endpoints):
        raise invalid("summary positive state disagrees with endpoint calls")
    if summary.get("positive_endpoint_count") != len(positive_endpoints):
        raise invalid("summary positive endpoint count is inconsistent")
    if summary.get("positive_endpoints") != positive_endpoints:
        raise invalid("summary positive endpoint list is inconsistent")
    substructure_hit_names = [
        rule_name
        for rule_name in CHEMFH_SUBSTRUCTURE_RULE_NAMES
        if rule_name not in missing_rules
        and rules[rule_name]["match_count"] > 0
    ]
    # New workers publish the explicit 7-model/10-substructure mapping.  Keep
    # old complete reports readable, but validate every field when present.
    optional_summary_contract = {
        "substructure_alert_count": len(substructure_hit_names),
        "substructure_alert_names": substructure_hit_names,
        "substructure_rule_names": list(CHEMFH_SUBSTRUCTURE_RULE_NAMES),
        "substructure_rule_mapping_version": (
            CHEMFH_SUBSTRUCTURE_RULE_MAPPING_VERSION
        ),
    }
    for field, expected in optional_summary_contract.items():
        if field in summary and summary.get(field) != expected:
            raise invalid(f"summary {field} is inconsistent")
    if summary.get("rules_contribute_to_decision") is not False:
        raise invalid("summary must keep native rules as evidence-only")
    if summary.get("pains_overlap_decision_suppressed") is not True:
        raise invalid("summary PAINS duplicate-decision suppression is missing")
    if status == "partial" and not (missing_endpoints or missing_rules):
        raise invalid("partial result declares no missing evidence")


def _chemfh_positive_state(chemfh: dict) -> bool | None:
    """Return ChemFH's model-positive state without using structural rules."""
    predictions = chemfh.get("predictions")
    if not isinstance(predictions, dict) or set(predictions) != set(
        CHEMFH_ENDPOINT_IDS
    ):
        return None
    calls = []
    for endpoint_id in CHEMFH_ENDPOINT_IDS:
        record = predictions.get(endpoint_id)
        if not isinstance(record, dict) or not isinstance(
            record.get("positive"), bool
        ):
            return None
        calls.append(record["positive"])
    return any(calls)


def _run_synthesis(
    input_data: dict,
    output_dir: Path,
    results: list[dict],
    *,
    enabled: bool,
    admet_required: bool,
    steps_completed: list[str],
    steps_skipped: list[dict],
) -> dict:
    config = _tool_config(input_data, "aizynthfinder")
    params = dict(config.get("params") or {})
    top_k = max(1, _integer(params.get("top_k"), 25))
    report: dict = {
        "success": True,
        "enabled": enabled,
        "selection_policy": {
            "top_k": top_k,
            "requires_complete_admet": admet_required,
            "accepts_sufficient_partial_admet": True,
            "requires_property_metric_gate_not_rejected": True,
            "excludes_benchmark_decoys": True,
            "sort": ["selection_rank", "record_index", "candidate_id"],
        },
        "selected_candidates": [],
        "results": [],
    }
    if not enabled:
        for result in results:
            _mark_synthesis_not_attempted(result, "disabled")
        steps_skipped.append({
            "step": "aizynthfinder",
            "reason": "AiZynthFinder explicitly disabled",
        })
        return report

    eligible: list[dict] = []
    excluded_decoys = 0
    for result in results:
        reason = None
        # Retrospective DEC_* decoys exist only to score the pipeline
        # (screening_metrics ranks them on docking scores alone) and can never
        # be nominated. Letting them compete for the top_k retrosynthesis budget
        # wastes AiZynthFinder time and pushes real candidates past top_k, where
        # they carry no route evidence at all.
        if _property_normalization_pool(result) != "nomination_candidates":
            reason = "benchmark_decoy"
            excluded_decoys += 1
        elif not result.get("valid_structure"):
            reason = "invalid_structure"
        elif result.get("pre_synthesis_classification") == "reject":
            # The category fold under the active screening mode has already
            # eliminated this molecule.  A bare ``verdict == "fail"`` mark is
            # deliberately not checked here: one threshold breach must not spend
            # the molecule's retrosynthesis slot before the fold has run.
            reason = "property_metric_gate_failed"
        elif admet_required and not _admet_evidence_allows_downstream(result):
            admet_protocol = (
                (result.get("property_protocol_evidence") or {}).get("admet")
                or {}
            )
            reason = (
                "admet_dependency_not_planned"
                if admet_protocol.get("planned") is False
                else "admet_protocol_evidence_incomplete"
            )
        if reason:
            _mark_synthesis_not_attempted(
                result, "not_eligible", reason=reason,
            )
        else:
            eligible.append(result)
    report["excluded_benchmark_decoys"] = excluded_decoys

    eligible.sort(key=_synthesis_sort_key)
    selected = eligible[:top_k]
    for position, result in enumerate(selected, start=1):
        _mark_synthesis_selected(result, position)
    for result in eligible[top_k:]:
        # top_k is a compute budget, not a finding about the molecule.  These
        # candidates were not rejected by any category fold; they simply sit
        # outside the slice this run could afford to search. The executable policy keeps
        # top-K as a selection convenience rather than a verdict, and the executable policy
        # only backtracks on evidence about the molecule, so the record keeps its
        # label and only ``synthesis.status`` says why no route evidence exists.
        _mark_synthesis_not_attempted(
            result,
            "not_selected_top_k",
            reason=f"outside top {top_k}",
        )

    candidates = [
        {"candidate_id": result["candidate_id"], "smiles": result["smiles"]}
        for result in selected
    ]
    report["selected_candidates"] = [
        {
            "top_k_rank": index,
            "candidate_id": result["candidate_id"],
            "mol_name": result["mol_name"],
            "selection_rank": result.get("selection_rank"),
            "smiles": result["smiles"],
        }
        for index, result in enumerate(selected, start=1)
    ]
    if not candidates:
        report["status"] = "no_eligible_candidates"
        steps_skipped.append({
            "step": "aizynthfinder",
            "reason": "No RDKit/ADMET-eligible candidate for retrosynthesis",
        })
        return report

    try:
        worker_report = run_aizynth_worker(
            candidates,
            output_dir / "aizynthfinder",
            runtime=config.get("runtime") or None,
            config_path=params.get("config_path"),
            timeout_seconds=_integer(params.get("timeout_seconds"), 21600),
            search_time_limit_seconds=_integer(
                params.get("search_time_limit_seconds"), 60
            ),
            iteration_limit=_integer(params.get("iteration_limit"), 100),
            max_routes=_integer(params.get("max_routes"), 5),
            return_first=bool(params.get("return_first", True)),
            expansion_policies=params.get("expansion_policies"),
            filter_policies=params.get("filter_policies"),
            stocks=params.get("stocks"),
        )
    except Exception as exc:
        # Runtime discovery, command construction, and process launch can fail
        # before the worker has a chance to serialize its normal error result.
        # Convert that exception into the same selected-row-scoped operational
        # failure contract; candidates outside the top-k remain unattempted.
        worker_report = {
            "success": False,
            "batch_status": "failed",
            "error_count": len(selected),
            "error": f"AiZynthFinder worker raised {type(exc).__name__}: {exc}",
            "files": {},
            "results": [],
        }
    report.update({
        key: value
        for key, value in worker_report.items()
        if key not in {"files"}
    })
    report["worker_files"] = worker_report.get("files") or {}
    returned_results = worker_report.get("results")
    structured_all_error = bool(
        worker_report.get("batch_status") == "failed"
        and worker_report.get("error_count") == len(selected)
        and isinstance(returned_results, list)
        and len(returned_results) == len(selected)
        and all(
            isinstance(entry, dict)
            and str(entry.get("status") or "").strip().casefold() == "error"
            for entry in returned_results
        )
    )
    if not worker_report.get("success") and not structured_all_error:
        error = str(worker_report.get("error") or "AiZynthFinder worker failed")
        report["batch_status"] = "failed"
        report["error_count"] = len(selected)
        for result in selected:
            result["synthesis"] = {
                **(result.get("synthesis") or {}),
                "status": "error",
                "attempted": True,
                "metrics": _empty_synthesis_metrics(),
                "error": error,
            }
            _mark_synthesis_operational_failure(result, error)
        report["results"] = [
            _failed_synthesis_result(result, error) for result in selected
        ]
        report["operational_errors"] = _report_operational_errors(selected)
        # The worker died for the whole batch, so the scope of the problem is
        # the protocol, not any one candidate: the executable policy scopes a single tool
        # backtrack to one candidate and the executable policy sends a batch-wide
        # systematic failure back to preparation instead.  Declaring it as
        # planned-but-failed keeps the synthesis metrics inside Sigma_planned so
        # coverage reports the gap and the failure registry drives the retry or
        # tool substitution.
        steps_skipped.append({
            "step": "aizynthfinder",
            "reason": f"{error}; planned but failed, see failure registry",
        })
        return report

    try:
        indexed = _strict_backend_index(
            candidates, worker_report, "AiZynthFinder"
        )
    except ValueError as exc:
        report["success"] = False
        report["batch_status"] = "failed"
        report["error_count"] = len(selected)
        report["error"] = str(exc)
        for result in selected:
            result["synthesis"] = {
                **(result.get("synthesis") or {}),
                "status": "error",
                "attempted": True,
                "metrics": _empty_synthesis_metrics(),
                "error": str(exc),
            }
            _mark_synthesis_operational_failure(result, str(exc))
        report["results"] = [
            _failed_synthesis_result(result, str(exc)) for result in selected
        ]
        report["operational_errors"] = _report_operational_errors(selected)
        # A malformed batch reply is a backend fault covering every candidate at
        # once, so it goes to the failure registry rather than onto each
        # molecule's label.
        steps_skipped.append({
            "step": "aizynthfinder",
            "reason": f"{exc}; planned but failed, see failure registry",
        })
        return report

    operational_error_count = 0
    for result in selected:
        route_result = indexed[result["candidate_id"]]
        route_status = str(route_result.get("status") or "").strip().casefold()
        route_error = str(route_result.get("error") or "").strip()
        if route_status not in SYNTHESIS_ATTEMPTED_STATUSES:
            route_error = (
                f"AiZynthFinder returned invalid status {route_status!r} for "
                f"{result['candidate_id']}"
            )
            route_status = "error"
        elif (
            route_status == "solved"
            and route_result.get("is_solved") is not True
        ):
            route_error = (
                "AiZynthFinder returned status='solved' without is_solved=true "
                f"for {result['candidate_id']}"
            )
            route_status = "error"
        elif route_status == "error" and not route_error:
            route_error = (
                "AiZynthFinder failed without an error message for "
                f"{result['candidate_id']}"
            )
        elif (
            route_status == "not_found_within_budget"
            and route_result.get("is_solved") is not False
        ):
            route_error = (
                "AiZynthFinder returned status='not_found_within_budget' "
                f"without is_solved=false for {result['candidate_id']}"
            )
            route_status = "error"

        raw_routes = route_result.get("routes")
        if route_status != "error" and not (
            raw_routes is None
            or isinstance(raw_routes, list)
            and all(isinstance(route, dict) for route in raw_routes)
        ):
            route_error = (
                "AiZynthFinder returned malformed routes for "
                f"{result['candidate_id']}"
            )
            route_status = "error"
        route_records = (
            (raw_routes or []) if route_status != "error" else []
        )
        supplied_metrics = route_result.get("metrics")
        if supplied_metrics is None:
            # Legacy workers predate the explicit scalar block; the fallback
            # below derives only values with stable meanings.
            supplied_metrics = {}
        elif not isinstance(supplied_metrics, dict):
            if route_status != "error":
                route_error = (
                    "AiZynthFinder returned malformed metrics for "
                    f"{result['candidate_id']}"
                )
                route_status = "error"
                route_records = []
            supplied_metrics = {}
        selected_route_index = _integer(
            supplied_metrics.get("selected_route_index"), 0
        )
        selected_route_index = (
            selected_route_index
            if 0 <= selected_route_index < len(route_records)
            else 0
        )
        starting_materials = (
            route_records[selected_route_index].get("starting_materials", [])
            if route_records
            else []
        )
        synthesis_metrics = (
            _empty_synthesis_metrics()
            if route_status == "error"
            else _synthesis_metrics(
                route_result,
                starting_materials=starting_materials,
            )
        )
        result["synthesis"] = {
            **(result.get("synthesis") or {}),
            "status": route_status,
            "attempted": True,
            "is_solved": route_status == "solved",
            "route_count": (
                0 if route_status == "error"
                else _integer(route_result.get("route_count"), 0)
            ),
            "statistics": (
                None if route_status == "error"
                else route_result.get("statistics")
            ),
            "metrics": synthesis_metrics,
            "error": route_error or None,
            "starting_materials": starting_materials,
        }
        # ``synthesis_report.json`` and the property row must agree: a failed
        # call carries explicit null scalar evidence and can never leak a stale
        # route_found=0 into downstream interpretation.
        route_result.update({
            "status": route_status,
            "attempted": True,
            "is_solved": route_status == "solved",
            "route_count": result["synthesis"]["route_count"],
            "statistics": result["synthesis"]["statistics"],
            "metrics": synthesis_metrics,
            "error": route_error or None,
            "routes": route_records,
        })
        if route_status == "error":
            operational_error_count += 1
            _mark_synthesis_operational_failure(result, route_error)
        elif route_status == "not_found_within_budget":
            # The search ran for this molecule and came back without a route.
            # That is evidence about the chemistry, so unlike a tool outage it
            # does belong on the record as MODIFY. An unsuccessful search is
            # a flag, not a hard rejection.
            result["nomination_eligibility"] = "review_only"
            molecule_reasons = result.setdefault("molecule_review_reasons", [])
            if "synthesis_route_unsolved" not in molecule_reasons:
                molecule_reasons.append("synthesis_route_unsolved")
            reasons = result.setdefault("nomination_eligibility_reasons", [])
            if "synthesis_route_unsolved" not in reasons:
                reasons.append("synthesis_route_unsolved")
    report["error_count"] = operational_error_count
    report["batch_status"] = (
        "complete" if operational_error_count == 0
        else "failed" if operational_error_count == len(selected)
        else "partial"
    )
    report["success"] = operational_error_count < len(selected)
    if operational_error_count:
        report["operational_errors"] = _report_operational_errors(selected)
    if operational_error_count == len(selected):
        report["error"] = (
            f"AiZynthFinder failed for all {len(selected)} candidates"
        )
        steps_skipped.append({
            "step": "aizynthfinder",
            "reason": f"{report['error']}; planned but failed, see failure registry",
        })
    else:
        steps_completed.append("aizynthfinder")
    return report


def _admet_evidence_allows_downstream(result: dict) -> bool:
    """Accept complete or coverage-sufficient planned ADMET evidence.

    An intentionally unplanned ADMET tool is outside the protocol denominator,
    so it is not converted into a backend failure; it nevertheless cannot
    satisfy this downstream dependency.  Current reports carry the explicit
    evidence object; ambiguous legacy partial/disabled reports remain held.
    """
    evidence = result.get("property_protocol_evidence")
    admet_evidence = evidence.get("admet") if isinstance(evidence, dict) else None
    if isinstance(admet_evidence, dict):
        if admet_evidence.get("planned") is False:
            return False
        return bool(
            admet_evidence.get("backend_status") == "complete"
            or admet_evidence.get("accepted_partial") is True
        )
    status = str((result.get("admet") or {}).get("status") or "missing").casefold()
    return status == "complete"


def _strict_backend_index(
    requested: list[dict],
    report: dict,
    backend_name: str,
) -> dict[str, dict]:
    returned = report.get("results")
    if not isinstance(returned, list):
        raise ValueError(f"{backend_name} returned no results list")
    requested_ids = [entry["candidate_id"] for entry in requested]
    returned_ids = [
        str(entry.get("candidate_id") or "")
        for entry in returned
        if isinstance(entry, dict)
    ]
    if len(returned_ids) != len(returned):
        raise ValueError(f"{backend_name} returned a non-object result")
    if len(set(returned_ids)) != len(returned_ids):
        raise ValueError(f"{backend_name} returned duplicate candidate_id values")
    if set(returned_ids) != set(requested_ids):
        missing = sorted(set(requested_ids) - set(returned_ids))
        extra = sorted(set(returned_ids) - set(requested_ids))
        raise ValueError(
            f"{backend_name} candidate_id mismatch: missing={missing}, extra={extra}"
        )
    indexed = {str(entry["candidate_id"]): entry for entry in returned}
    for candidate in requested:
        returned_smiles = str(
            indexed[candidate["candidate_id"]].get("smiles") or ""
        )
        if returned_smiles != candidate["smiles"]:
            raise ValueError(
                f"{backend_name} SMILES mismatch for {candidate['candidate_id']}"
            )
    return indexed


def _synthesis_metrics(
    route_result: dict,
    *,
    starting_materials: list[dict],
) -> dict:
    """Normalize the current AiZynthFinder worker's stable scalar contract."""
    supplied = route_result.get("metrics")
    if isinstance(supplied, dict):
        scorer_names = [
            str(name) for name in (supplied.get("state_score_scorers") or [])
        ]
        state_score_is_applicable = (
            bool(supplied.get("state_score_applicable", False))
            and [name.strip().casefold() for name in scorer_names]
            == ["state score"]
        )
        state_score = (
            _finite_number_or_none(supplied.get("state_score"))
            if state_score_is_applicable
            else None
        )
        return {
            "route_found": _finite_number_or_none(supplied.get("route_found")),
            "selected_route_index": _finite_number_or_none(
                supplied.get("selected_route_index")
            ),
            "selected_route_rank": _finite_number_or_none(
                supplied.get("selected_route_rank")
            ),
            "state_score": state_score,
            "state_score_applicable": state_score is not None,
            "state_score_scorers": scorer_names,
            "reaction_count": _finite_number_or_none(
                supplied.get("reaction_count")
            ),
            "precursor_count": _finite_number_or_none(
                supplied.get("precursor_count")
            ),
            "precursors_in_stock_count": _finite_number_or_none(
                supplied.get("precursors_in_stock_count")
            ),
            "precursors_in_stock_fraction": _finite_number_or_none(
                supplied.get("precursors_in_stock_fraction")
            ),
        }

    # Backward-compatible worker fixtures have no explicit scorer identity, so
    # do not infer a State score from an opaque top_score. The remaining values
    # have stable meanings in the compact route result.
    statistics = route_result.get("statistics") or {}
    material_count = len(starting_materials)
    in_stock_count = sum(
        bool(item.get("in_stock"))
        for item in starting_materials
        if isinstance(item, dict)
    )
    status = str(route_result.get("status") or "")
    return {
        "route_found": (
            1.0 if route_result.get("is_solved") is True
            else 0.0 if status == "not_found_within_budget"
            else None
        ),
        "selected_route_index": 0.0 if starting_materials else None,
        "selected_route_rank": 1.0 if starting_materials else None,
        "state_score": None,
        "state_score_applicable": False,
        "state_score_scorers": [],
        "reaction_count": _finite_number_or_none(
            statistics.get("number_of_steps")
        ),
        "precursor_count": float(material_count) if starting_materials else None,
        "precursors_in_stock_count": (
            float(in_stock_count) if starting_materials else None
        ),
        "precursors_in_stock_fraction": (
            in_stock_count / material_count if material_count else None
        ),
    }


def _empty_synthesis_metrics() -> dict:
    """Return the scalar schema for an attempted synthesis call with no result."""
    return {
        "route_found": None,
        "selected_route_index": None,
        "selected_route_rank": None,
        "state_score": None,
        "state_score_applicable": False,
        "state_score_scorers": [],
        "reaction_count": None,
        "precursor_count": None,
        "precursors_in_stock_count": None,
        "precursors_in_stock_fraction": None,
    }


def _synthesis_sort_key(result: dict) -> tuple:
    rank = result.get("selection_rank")
    return (
        int(rank) if isinstance(rank, int) else math.inf,
        int(result.get("record_index", 0)),
        str(result.get("candidate_id") or ""),
    )


def _input_property_scoring_plan(input_data: dict) -> dict:
    """Copy the resolved plan contract into every candidate result.

    The execution switches in ``tools`` answer which backend process runs.
    These fields answer the separate scientific questions "was this evidence
    planned?", "was it required?", and "was it explicitly N/A?".  Keeping the
    declaration on each row lets scoring, coordinator joins, and reports make
    the same decision instead of reconstructing it from a backend status.
    """
    contract: dict = {}
    nested = input_data.get("scoring_plan")
    for source in (input_data, nested if isinstance(nested, dict) else {}):
        for field in PROPERTY_PLAN_FIELDS:
            if field in source and source[field] is not None:
                contract[field] = deepcopy(source[field])
    return contract


def _resolve_property_tools(
    input_data: dict,
) -> tuple[set, bool, bool, bool, bool]:
    tools = input_data.get("tools")
    if isinstance(tools, dict) and tools:
        def enabled(name: str, default: bool) -> bool:
            value = tools.get(name)
            if not isinstance(value, dict):
                return default
            return bool(value.get("enabled", default))

        checks = {
            name
            for name in ("lipinski", "qed", "sa_score", "pains")
            if enabled(name, True)
        }
        return (
            checks,
            enabled("admet", True),
            # Older standalone tool blocks predate ChemFH. Keep those explicit
            # blocks opt-in; registry-resolved final plans always carry this key.
            enabled("chemfh", False),
            enabled("dedup", True),
            enabled("aizynthfinder", True),
        )

    checks = {"lipinski", "qed", "sa_score"}
    if input_data.get("check_pains", True):
        checks.add("pains")
    return (
        checks,
        bool(input_data.get("run_admet", True)),
        bool(input_data.get("run_chemfh", True)),
        bool(input_data.get("deduplicate_compounds", True)),
        bool(input_data.get("run_aizynthfinder", True)),
    )


def _tool_config(input_data: dict, tool_id: str) -> dict:
    tools = input_data.get("tools")
    value = tools.get(tool_id) if isinstance(tools, dict) else None
    return dict(value) if isinstance(value, dict) else {}


def _rdkit_step_status(
    valid_candidates: int,
    enabled_checks: set,
    *,
    sa_available: bool,
    alert_backend_available: bool,
) -> tuple[list[str], list[dict]]:
    completed: list[str] = []
    skipped: list[dict] = []
    if valid_candidates:
        for step in ("lipinski", "qed"):
            if step in enabled_checks:
                completed.append(step)
        if "sa_score" in enabled_checks:
            if sa_available:
                completed.append("sa_score")
            else:
                skipped.append({
                    "step": "sa_score",
                    "reason": "RDKit contrib sascorer is unavailable",
                })
        if "pains" in enabled_checks:
            if alert_backend_available:
                completed.append("structural_alerts")
            else:
                skipped.append({
                    "step": "structural_alerts",
                    "reason": "RDKit FilterCatalog is unavailable",
                })
    return completed, skipped


def _legacy_pains_projection(alerts: dict) -> dict:
    matches = [
        match
        for match in (alerts.get("matches") or [])
        if str(match.get("catalog") or "").startswith("PAINS_")
    ]
    is_pains = None if alerts.get("available") is False else bool(matches)
    return {
        "is_pains": is_pains,
        "matched_filters": [match.get("description") for match in matches],
        "matches": matches,
    }


def _build_summary(
    records: list[dict],
    groups: list[dict],
    results: list[dict],
    admet_enabled: bool,
    chemfh_enabled: bool,
    synthesis_enabled: bool,
) -> dict:
    verdicts = [str(result.get("verdict") or "unknown") for result in results]
    admet_statuses = [
        str((result.get("admet") or {}).get("status") or "missing")
        for result in results
        if result.get("valid_structure")
    ]
    chemfh_statuses = [
        str((result.get("chemfh") or {}).get("status") or "missing")
        for result in results
        if result.get("valid_structure")
    ]
    chemfh_positive = sum(
        _chemfh_positive_state(result.get("chemfh") or {}) is True
        for result in results
        if result.get("valid_structure")
    )
    synthesis_records = [
        result.get("synthesis") or {}
        for result in results
        if result.get("valid_structure")
    ]
    synthesis_statuses = [
        str(record.get("status") or "missing") for record in synthesis_records
    ]
    synthesis_attempted = [
        _synthesis_attempted(record) for record in synthesis_records
    ]
    valid_input_records = sum(record.get("mol") is not None for record in records)
    duplicate_pose_records = len(records) - len(groups)
    return {
        "input_records": len(records),
        "valid_input_records": valid_input_records,
        "invalid_input_records": len(records) - valid_input_records,
        "valid_records": valid_input_records,
        "invalid_records": len(records) - valid_input_records,
        "unique_candidates": len(groups),
        "deduplicated_records": duplicate_pose_records,
        "duplicate_pose_records": duplicate_pose_records,
        # Backward-compatible report keys.
        "total_compounds": len(results),
        "total_molecules": len(results),
        "pass": verdicts.count("pass"),
        "flag": verdicts.count("flag"),
        "fail": verdicts.count("fail"),
        "not_evaluated": verdicts.count("not_evaluated"),
        "review_only": sum(
            result.get("nomination_eligibility") == "review_only"
            for result in results
        ),
        "admet": {
            "enabled": admet_enabled,
            "eligible_candidates": len(admet_statuses),
            "complete": admet_statuses.count("complete"),
            "partial": admet_statuses.count("partial"),
            "error": admet_statuses.count("error"),
            "coverage_fraction": (
                round(admet_statuses.count("complete") / len(admet_statuses), 4)
                if admet_statuses
                else 0.0
            ),
        },
        "admet_requested": len(admet_statuses),
        "admet_complete": admet_statuses.count("complete"),
        "admet_partial": admet_statuses.count("partial"),
        "admet_error": admet_statuses.count("error"),
        "chemfh": {
            "enabled": chemfh_enabled,
            "eligible_candidates": len(chemfh_statuses),
            "complete": chemfh_statuses.count("complete"),
            "partial": chemfh_statuses.count("partial"),
            "error": chemfh_statuses.count("error"),
            "disabled": chemfh_statuses.count("disabled"),
            "positive_predictions": chemfh_positive,
            "coverage_fraction": (
                round(
                    chemfh_statuses.count("complete") / len(chemfh_statuses),
                    4,
                )
                if chemfh_statuses
                else 0.0
            ),
            "rule_decision_role": "provider_evidence_only",
            "pains_overlap_decision_suppressed": True,
        },
        "chemfh_requested": len(chemfh_statuses),
        "chemfh_complete": chemfh_statuses.count("complete"),
        "chemfh_partial": chemfh_statuses.count("partial"),
        "chemfh_error": chemfh_statuses.count("error"),
        "chemfh_positive": chemfh_positive,
        "synthesis": {
            "enabled": synthesis_enabled,
            "selected": synthesis_attempted.count(True),
            "attempted": synthesis_attempted.count(True),
            "not_attempted": synthesis_attempted.count(False),
            "attempt_status_unknown": synthesis_attempted.count(None),
            "solved": synthesis_statuses.count("solved"),
            "not_found_within_budget": synthesis_statuses.count(
                "not_found_within_budget"
            ),
            "error": synthesis_statuses.count("error"),
            "not_eligible": synthesis_statuses.count("not_eligible"),
            "not_selected_top_k": synthesis_statuses.count(
                "not_selected_top_k"
            ),
            "disabled": synthesis_statuses.count("disabled"),
        },
        "synthesis_selected": synthesis_attempted.count(True),
        "synthesis_attempted": synthesis_attempted.count(True),
        "synthesis_not_attempted": synthesis_attempted.count(False),
        "synthesis_solved": synthesis_statuses.count("solved"),
        "synthesis_not_found": synthesis_statuses.count(
            "not_found_within_budget"
        ),
        "synthesis_error": synthesis_statuses.count("error"),
        "synthesis_not_eligible": synthesis_statuses.count("not_eligible"),
        "synthesis_not_selected_top_k": synthesis_statuses.count(
            "not_selected_top_k"
        ),
        "synthesis_disabled": synthesis_statuses.count("disabled"),
    }


def _synthesis_attempted(synthesis: dict) -> bool | None:
    """Read the explicit attempt bit, with a status fallback for older reports."""
    attempted = synthesis.get("attempted")
    if isinstance(attempted, bool):
        return attempted
    status = str(synthesis.get("status") or "").strip().casefold()
    if status in SYNTHESIS_ATTEMPTED_STATUSES:
        return True
    if status in SYNTHESIS_NOT_ATTEMPTED_STATUSES:
        return False
    return None


def _effective_property_metric_rules(
    thresholds: dict,
) -> tuple[dict, ...]:
    """Apply explicit Step 5 soft-threshold overrides without weakening hard gates."""
    rules = list(deepcopy(PROPERTY_FILTER_METRIC_RULES))
    qed_soft = max(
        0.25,
        min(1.0, _number(thresholds.get("qed_threshold"), 0.45)),
    )
    sa_soft = max(
        1.0,
        min(8.0, _number(thresholds.get("sa_threshold"), 6.0)),
    )
    lipinski_max = max(
        0,
        min(
            4,
            _integer(
                thresholds.get("lipinski_max_violations"),
                DEFAULT_LIPINSKI_MAX_VIOLATIONS,
            ),
        ),
    )
    bbb_profile = str(thresholds.get("bbb_profile") or "non_cns").lower()
    if bbb_profile not in {"non_cns", "cns"}:
        raise ValueError("bbb_profile must be 'non_cns' or 'cns'")
    for rule in rules:
        metric_id = rule.get("metric_id")
        if metric_id == "rdkit.qed":
            rule["soft_lower"] = qed_soft
        elif metric_id == "rdkit.sa_score":
            rule["soft_upper"] = sa_soft
        elif metric_id == "lipinski.ro5_violations":
            # RO5 marks: 0-1 pass, 2-3 flag, >3 fail. The soft bound stays
            # at 1 so a single violation is not marked at all; only the hard
            # bound is operator-configurable, and it can never drop below the
            # soft bound.
            rule["soft_upper"] = min(1, lipinski_max)
            rule["hard_upper"] = max(lipinski_max, rule["soft_upper"])
        elif metric_id == "rdkit.tpsa" and bbb_profile == "cns":
            # The CNS profile fixes TPSA soft/hard upper bounds at 140/180,
            # currently the same window as the general oral profile. BBB_Martins
            # below supplies the barrier requirement. Restating the CNS bounds
            # keeps this profile stable if the general row is later retuned.
            rule.update({
                "soft_upper": 140.0,
                "hard_upper": 180.0,
                "description": (
                    "CNS profile: <=140 passes; 140-180 flags; "
                    ">180 fails -- the same upper window as the general oral "
                    "profile. BBB_Martins carries the barrier requirement."
                ),
            })
        elif metric_id == "admet_ai.BBB_Martins" and bbb_profile == "cns":
            # The executable policy -- the CNS profile mirrors the default row rather than
            # inventing a new window: the non-CNS band marks `<=0.30` pass and
            # `>0.90` fail, so reversing it means `>=0.90` pass and `<0.30`
            # fail with the same two numbers swapped onto the lower side.
            rule.update({
                "hard_lower": 0.3,
                "soft_lower": 0.9,
                "soft_upper": None,
                "hard_upper": None,
                "direction": "higher",
                "signed_weight": abs(float(rule.get("weight", 0.0))),
                "description": (
                    "CNS profile: >=0.90 passes; 0.30-0.90 flags; <0.30 "
                    "fails. High BBB penetration is rewarded."
                ),
            })
    return tuple(rules)


def _apply_property_metric_policy(
    results: list[dict],
    *,
    rules: tuple[dict, ...],
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> dict:
    """Gate and score independent scientific-cohort candidate/decoy pools."""
    mode = resolve_screening_mode(screening_mode)
    grouped: dict[tuple[str, str], list[int]] = {}
    for index, result in enumerate(results):
        cohort_id = _property_cohort_id(result)
        normalization_pool = _property_normalization_pool(result)
        result["cohort_id"] = cohort_id
        result["normalization_pool"] = normalization_pool
        grouped.setdefault((cohort_id, normalization_pool), []).append(index)

    cohorts: dict[str, dict] = {}
    for (cohort_id, normalization_pool), indices in grouped.items():
        scoring_inputs = []
        for index in indices:
            scoring_input = dict(results[index])
            # Invalid structures stay in the audit set.  Do not pass a prior
            # molecule classification into score_records: the current metric
            # fold must be recomputed for this round/mode.  The explicit
            # structural marker is the only upstream gate that remains hard.
            structural_reject = record_has_structural_reject(scoring_input)
            if structural_reject:
                scoring_input["structural_reject"] = True
            scoring_input.pop("intrinsic_classification", None)
            scoring_inputs.append(scoring_input)
        scored, pool_metadata = score_records(
            scoring_inputs,
            rules=rules,
            normalization_method=NORMALIZATION_METHOD,
            # This skill only ever produces step 4 evidence.  Declaring it keeps
            # coverage from being measured against steps this run was never
            # asked to execute.
            active_steps=("step4",),
            screening_mode=mode,
        )
        pool_metadata["normalization_scope"] = (
            "benchmark_decoys_in_scientific_cohort"
            if normalization_pool == "benchmark_decoys"
            else "nomination_candidates_in_scientific_cohort"
        )
        pool_metadata["cohort_id"] = cohort_id
        for index, scored_result in zip(indices, scored):
            structural_reject = record_has_structural_reject(scored_result)
            metric_gate = str(
                scored_result.get("metric_gate") or "not_evaluated"
            )
            # Coverage is deliberately absent here. The executable policy uses it
            # only for rank release/backtracking; ``classification`` is the
            # current total-indicator fold (plus an explicit structural hard
            # failure), never a stale legacy field or a coverage downgrade.
            classification = (
                "reject" if structural_reject else metric_gate
            )
            scored_result["legacy_gate"] = _legacy_gate_classification(
                scored_result
            )
            scored_result["property_metric_gate"] = metric_gate
            scored_result["classification"] = classification
            scored_result["gate"] = classification
            scored_result["structural_reject"] = bool(structural_reject)
            scored_result["property_weighted_total_score"] = scored_result.get(
                "weighted_total_score"
            )
            scored_result["property_weighted_score_fraction"] = scored_result.get(
                "weighted_score_fraction"
            )
            scored_result["property_weighted_score_rank"] = scored_result.get(
                "weighted_score_rank"
            )
            scored_result["property_overall_score"] = scored_result.get(
                "overall_score"
            )
            scored_result["property_overall_rank"] = scored_result.get(
                "overall_rank"
            )
            scored_result["property_rank_status"] = scored_result.get(
                "rank_status"
            )
            scored_result["property_category_scores"] = scored_result.get(
                "category_scores"
            )
            scored_result["property_step_scores"] = scored_result.get(
                "step_scores"
            )
            _apply_fold_nomination_eligibility(scored_result, classification)
            scored_result["metric_fold_reject"] = (
                classification == "reject" and not structural_reject
            )
            results[index] = scored_result
        cohorts.setdefault(cohort_id, {})[normalization_pool] = pool_metadata

    # Keep the physical report/SDF/manifest order aligned with local ranks,
    # never by comparing score magnitudes across scientific cohorts or against
    # the isolated retrospective decoy pool.
    results.sort(key=_property_weighted_total_sort_key)

    scoring_metadata = {
        "schema_version": SCORING_SCHEMA_VERSION,
        "metric_rules_version": METRIC_RULES_VERSION,
        "normalization_method": NORMALIZATION_METHOD,
        # Level 1 marks each metric pass/flag/fail; level 2 folds those marks
        # per parent category and only then decides keep/modify/reject.
        "gate_contract": "two_level_label_then_category_fold",
        "screening_mode": dict(mode),
        "normalization_scope": "scientific_cohort_candidate_and_decoy_pools",
        "cohorts": cohorts,
        "cohort_count": len(cohorts),
        "cohort_fields": [
            "protein_state_id",
            "pocket_id",
            "round_id",
            "protocol_profile_id",
            "metric_catalog_version",
        ],
        "decoys_isolated_from_nomination_normalization": True,
    }
    return {
        "schema_version": SCORING_SCHEMA_VERSION,
        "metric_rules_version": METRIC_RULES_VERSION,
        "normalization_method": NORMALIZATION_METHOD,
        "rules": public_metric_rules(rules),
        "scoring": scoring_metadata,
    }


def _apply_fold_nomination_eligibility(result: dict, classification: str) -> None:
    """Let the current category fold own molecular nomination eligibility.

    First-pass RDKit/ChemFH marks remain audit reasons; protocol evidence can
    still require review even when the molecule passes the category fold.
    """
    if classification == "reject":
        result["nomination_eligibility"] = "rejected"
    elif classification == "keep":
        protocol_reasons = list(result.get("protocol_review_reasons") or [])
        result["nomination_eligibility"] = (
            "review_only" if protocol_reasons else "auto_candidate"
        )
        result["nomination_eligibility_reasons"] = protocol_reasons
    else:
        result["nomination_eligibility"] = "review_only"


def _apply_pre_synthesis_metric_gates(
    results: list[dict],
    rules: tuple[dict, ...],
    *,
    screening_mode: object = DEFAULT_SCREENING_MODE,
) -> None:
    """Apply completed non-SYN raw gates before retrosynthesis.

    Percentile/category scoring is finalized after the synthesis evidence is
    attached. This pass exists solely to enforce the invariant that a molecule
    whose *completed* category fold already trips the active screening mode never
    advances to the next backend.  SYN is deliberately deferred because a
    RAScore/SA fail ratio computed before AiZynthFinder can fall when the pending
    route labels arrive.  This is distinct from a final NA: once a planned tool
    has completed or failed, NA is excluded from the denominator as usual.
    """
    mode = resolve_screening_mode(screening_mode)
    completed_rules = tuple(
        rule for rule in rules if str(rule.get("category_code")) != "SYN"
    )
    for result in results:
        metric_gate = aggregate_molecule_gate(
            assess_record(result, completed_rules), screening_mode=mode
        )
        structural_reject = record_has_structural_reject(result)
        classification = "reject" if structural_reject else metric_gate
        result["pre_synthesis_metric_gate"] = metric_gate
        result["pre_synthesis_classification"] = classification
        result["structural_reject"] = bool(structural_reject)
        result["screening_mode"] = mode["mode_id"]
        _apply_fold_nomination_eligibility(result, classification)
        result["metric_fold_reject"] = (
            classification == "reject" and not structural_reject
        )


def _parameters(
    thresholds: dict,
    alert_catalogs: list,
    hard_fail_catalogs: list,
    is_covalent: bool,
    dedup_enabled: bool,
    input_data: dict | None = None,
) -> dict:
    parameters = {
        **thresholds,
        "alert_catalogs": list(alert_catalogs),
        "hard_fail_catalogs": list(hard_fail_catalogs),
        "is_covalent": is_covalent,
        "deduplicate": dedup_enabled,
    }
    parameters.update({
        field: (input_data or {}).get(field)
        for field in (
            "cohort_id", "protein_state_id", "pocket_id", "round_id",
            "protocol_profile_id", "metric_catalog_version",
            "cohort_identity_complete", "cohort_identity_missing_fields",
            "cohort_identity_derived_fields", "cohort_provenance_sources",
        )
        if (input_data or {}).get(field) not in (None, "")
    })
    return parameters


def _write_filtered_outputs(
    input_data: dict,
    output_dir: Path,
    representatives: dict[str, dict],
    results: list[dict],
) -> dict:
    from rdkit import Chem

    # "not_evaluated" means the RDKit checks were all switched off, not that the
    # molecule failed. Keeping it in the default carry-forward set stops a
    # checks-off configuration from silently emptying the filtered library; the
    # classification is still "modify", so it cannot be nominated.
    # "fail" is likewise a level-1 mark on the RDKit screen, not a decision: the
    # ``classification != "reject"`` test below is the level-2 category fold and
    # is what actually removes a molecule.  Excluding "fail" here
    # would re-impose the old one-level gate on top of the fold.
    allowed_verdicts = set(
        input_data.get("filtered_verdicts")
        or ["pass", "flag", "fail", "not_evaluated"]
    )
    filtered_sdf = output_dir / (
        input_data.get("filtered_sdf_name") or "filtered_molecules.sdf"
    )
    filtered_manifest = output_dir / (
        input_data.get("filtered_manifest_name") or "filter_manifest.csv"
    )
    kept_rows: list[dict] = []
    removed_rows: list[dict] = []
    manifest_rows: list[dict] = []
    writer = Chem.SDWriter(str(filtered_sdf))
    try:
        for result in results:
            representative = representatives.get(result["candidate_id"]) or {}
            mol = representative.get("mol")
            row = _filter_manifest_row(result)
            if (
                result.get("verdict") in allowed_verdicts
                and result.get("classification") != "reject"
                and mol is not None
            ):
                output_mol = Chem.Mol(mol)
                output_mol.SetProp(
                    "property_prefilter_verdict", str(result["verdict"])
                )
                output_mol.SetProp(
                    "property_prefilter_flags",
                    "; ".join(result.get("flags") or []),
                )
                output_mol.SetProp(
                    "property_prefilter_gate",
                    str(result.get("classification") or "not_evaluated"),
                )
                for field in LINEAGE_PROPERTIES:
                    value = result.get(field)
                    if value in (None, ""):
                        continue
                    if isinstance(value, (dict, list, bool)):
                        value = json.dumps(
                            value,
                            sort_keys=True,
                            separators=(",", ":"),
                            ensure_ascii=False,
                        )
                    output_mol.SetProp(field, str(value))
                writer.write(output_mol)
                kept_rows.append(row)
                manifest_rows.append({**row, "kept": "true"})
            else:
                removed_rows.append(row)
                manifest_rows.append({**row, "kept": "false"})
    finally:
        writer.close()
    _write_filter_manifest(filtered_manifest, manifest_rows)
    return {
        "summary": {
            "filtered_total": len(kept_rows),
            "filtered_removed": len(removed_rows),
            "filtered_allowed_verdicts": sorted(allowed_verdicts),
            "filtered_pass": sum(row["verdict"] == "pass" for row in kept_rows),
            "filtered_flag": sum(row["verdict"] == "flag" for row in kept_rows),
            "filtered_fail": sum(row["verdict"] == "fail" for row in kept_rows),
        },
        "files": {
            "prefiltered_library": str(filtered_sdf),
            "prefiltered_library_manifest": str(filtered_manifest),
        },
    }


def _filter_manifest_row(result: dict) -> dict:
    properties = result.get("properties") or {}
    return {
        "record_index": result.get("record_index", ""),
        "mol_index": result.get("mol_index", ""),
        "candidate_id": result.get("candidate_id", ""),
        "compound_id": result.get("compound_id", ""),
        "chemical_identity_id": result.get("chemical_identity_id", ""),
        "pose_id": result.get("pose_id", ""),
        "selection_rank": result.get("selection_rank", ""),
        "cohort_id": result.get("cohort_id", ""),
        "normalization_pool": result.get("normalization_pool", ""),
        "protein_state_id": result.get("protein_state_id", ""),
        "pocket_id": result.get("pocket_id", ""),
        "round_id": result.get("round_id", ""),
        "protocol_profile_id": result.get("protocol_profile_id", ""),
        "metric_catalog_version": result.get("metric_catalog_version", ""),
        "cohort_identity_complete": result.get(
            "cohort_identity_complete", ""
        ),
        "cohort_identity_missing_fields": json.dumps(
            result.get("cohort_identity_missing_fields") or [],
            sort_keys=True,
            separators=(",", ":"),
        ),
        "cohort_identity_derived_fields": json.dumps(
            result.get("cohort_identity_derived_fields") or [],
            sort_keys=True,
            separators=(",", ":"),
        ),
        "cohort_provenance_sources": json.dumps(
            result.get("cohort_provenance_sources") or {},
            sort_keys=True,
            separators=(",", ":"),
        ),
        "mol_name": result.get("mol_name", ""),
        "source": result.get("source", ""),
        "smiles": result.get("smiles", ""),
        "verdict": result.get("verdict", ""),
        "classification": result.get("classification", ""),
        "flags": "; ".join(result.get("flags") or []),
        "MW": properties.get("MW", ""),
        "LogP": properties.get("LogP", ""),
        "QED": properties.get("QED", ""),
        "SA_score": properties.get("SA_score", ""),
        "Lipinski_violations": (result.get("lipinski") or {}).get(
            "violations", ""
        ),
        "property_weighted_total_score": result.get(
            "property_weighted_total_score", ""
        ),
        "property_weighted_score_rank": result.get(
            "property_weighted_score_rank", ""
        ),
        "property_overall_score": result.get("property_overall_score", ""),
        "property_overall_rank": result.get("property_overall_rank", ""),
        "property_rank_status": result.get("property_rank_status", ""),
        "kept": "",
    }


def _write_filter_manifest(
    path: Path,
    rows: list[dict],
) -> None:
    fieldnames = list(_filter_manifest_row({}).keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _property_weighted_total_sort_key(result: dict) -> tuple:
    """Sort by cohort/pool-local OverallRank; never compare cohort scores."""
    overall_rank = _finite_number_or_none(result.get("property_overall_rank"))
    rank_key = (0, overall_rank) if overall_rank is not None else (1, 0.0)
    return (
        str(result.get("cohort_id") or "cohort:unspecified"),
        1 if result.get("normalization_pool") == "benchmark_decoys" else 0,
        *rank_key,
        _integer(result.get("record_index"), 0),
        str(result.get("mol_name") or ""),
    )


def _write_final_report(
    output_dir: Path,
    *,
    success: bool,
    summary: dict,
    results: list[dict],
    steps_completed: list[str],
    steps_skipped: list[dict],
    parameters: dict,
    backend_reports: dict,
    error: str | None = None,
    diagnostics: dict | None = None,
    synthesis_report: dict | None = None,
    extra_files: dict | None = None,
    metric_policy: dict | None = None,
) -> dict:
    report_path = output_dir / "property_filter_report.json"
    summary_path = output_dir / "property_filter_summary.json"
    synthesis_path = output_dir / "synthesis_report.json"
    if metric_policy is None:
        metric_policy = _apply_property_metric_policy(
            results,
            rules=_effective_property_metric_rules(parameters),
        )
    _update_classification_summary(results, summary)
    synthesis_payload = synthesis_report or {
        "success": True,
        "enabled": False,
        "status": "not_run",
        "results": [],
    }
    _write_json(summary_path, summary)
    _write_json(synthesis_path, synthesis_payload)
    report = {
        "success": bool(success),
        "summary": summary,
        "results": results,
        "steps_completed": list(dict.fromkeys(steps_completed)),
        "steps_skipped": steps_skipped,
        "parameters": parameters,
        "backends": backend_reports,
        "metric_rules_version": metric_policy["metric_rules_version"],
        "metric_rules": metric_policy["rules"],
        "metric_scoring": metric_policy["scoring"],
        "files": {
            "property_filter_report_json": str(report_path),
            "property_filter_summary_json": str(summary_path),
            "synthesis_report_json": str(synthesis_path),
            **(extra_files or {}),
        },
    }
    operational_errors = _report_operational_errors(results)
    if operational_errors:
        report["operational_errors"] = operational_errors
        report["partial"] = True
    cohort_fields = (
        "cohort_id", "protein_state_id", "pocket_id", "round_id",
        "protocol_profile_id", "metric_catalog_version",
        "cohort_identity_complete", "cohort_identity_missing_fields",
        "cohort_identity_derived_fields", "cohort_provenance_sources",
    )
    cohort_manifest: dict[str, dict] = {}
    for row in results or [parameters]:
        identity = {field: row.get(field) for field in cohort_fields}
        if identity.get("cohort_id") not in (None, ""):
            cohort_manifest.setdefault(str(identity["cohort_id"]), identity)
    report["scientific_cohorts"] = list(cohort_manifest.values())
    if error:
        report["error"] = error
    if diagnostics:
        report["diagnostics"] = diagnostics
    _write_json(report_path, report)
    return report


def _legacy_gate_classification(result: dict) -> str:
    """Map the Step 4 legacy verdict/status fields to the stable three-tier gate.

    ``verdict == "fail"`` is a level-1 mark from the RDKit first-pass screen, not
    a molecule decision: under the executable policy only the per-category fail count and
    ratio may reject, so it is capped at ``modify`` here and the rejection has to
    come from ``metric_gate``.  Explicit structural markers (including
    ``valid_structure is False``) stay hard rejects; a bare
    ``nomination_eligibility="rejected"`` is not enough because older reports
    also used it for threshold-based metric folds.
    """
    verdict = str(result.get("verdict") or "").strip().lower()
    nomination = str(
        result.get("nomination_eligibility") or ""
    ).strip().lower()
    # ``nomination_eligibility=rejected`` is also used for a metric category
    # fold.  Only an explicit structural marker (or an identity/parse failure
    # recognized by the shared helper) is a hard rejection here; the metric
    # fold is recomputed by score_records under the active V1/V2/V3 mode.
    if record_has_structural_reject(result):
        return "reject"
    if verdict in ("fail", "flag", "not_evaluated"):
        # "not_evaluated" is absence of evidence, so it can never fall through to
        # "keep"; it blocks nomination and asks for review like any other flag.
        #
        # ...unless the whole verdict rests on a coverage gap.  A missing
        # sascorer, an unavailable FilterCatalog, a fully disabled RDKit pass or
        # a raised RDKit pass all produce a non-pass verdict without saying
        # anything about the chemistry, and the executable policy forbids letting a tool
        # problem rewrite the molecule.  Those rows fall through to the
        # review_only branch, which blocks nomination through
        # ``protocol_review_reasons`` while leaving the chemistry at "keep".
        if _rdkit_verdict_is_molecular(result):
            return "modify"
    if nomination == "review_only":
        molecule_reasons = result.get("molecule_review_reasons")
        protocol_reasons = result.get("protocol_review_reasons")
        if isinstance(molecule_reasons, list) or isinstance(protocol_reasons, list):
            # New records distinguish a molecular finding from a protocol-only
            # evidence gap.  The latter blocks nomination through its explicit
            # evidence object but must not rewrite the chemistry as MODIFY.
            return "modify" if molecule_reasons else "keep"
        # Reports written before the split used review_only only as a molecular
        # cap, so retain their historical interpretation.
        return "modify"
    return "keep"


def _update_classification_summary(results: list[dict], summary: dict) -> None:
    """Summarize final combined legacy and executable metric decisions.

    Counting only -- this must not rewrite ``classification``/``gate``.  It runs
    after ``_write_filtered_outputs``, so coercing an unrecognized label to
    ``reject`` here used to publish a molecule into ``filtered_molecules.sdf``
    while the report called it rejected.  A label this function does not know is
    missing evidence, not a molecular failure, so it is counted as
    ``not_evaluated`` and left alone.
    """
    counts = {"keep": 0, "modify": 0, "reject": 0, "not_evaluated": 0}
    for result in results:
        classification = str(result.get("classification") or "not_evaluated")
        if classification not in counts:
            classification = "not_evaluated"
        counts[classification] += 1
    summary.update(counts)


def _failure_report(
    output_dir: Path,
    error: str,
    *,
    summary: dict | None = None,
    diagnostics: dict | None = None,
) -> dict:
    normalized_summary = {
        "input_records": 0,
        "valid_input_records": 0,
        "invalid_input_records": 0,
        "valid_records": 0,
        "invalid_records": 0,
        "unique_candidates": 0,
        "deduplicated_records": 0,
        "duplicate_pose_records": 0,
    }
    normalized_summary.update(summary or {})
    if not summary or "valid_records" not in summary:
        normalized_summary["valid_records"] = normalized_summary.get(
            "valid_input_records", 0
        )
    if not summary or "invalid_records" not in summary:
        normalized_summary["invalid_records"] = normalized_summary.get(
            "invalid_input_records", 0
        )
    return _write_final_report(
        output_dir,
        success=False,
        error=error,
        summary=normalized_summary,
        results=[],
        steps_completed=[],
        steps_skipped=[],
        parameters={},
        backend_reports={},
        diagnostics=diagnostics,
    )


def _backend_provenance(report: dict) -> dict:
    return {
        "success": bool(report.get("success")),
        "error": report.get("error"),
        "backend": report.get("backend"),
        "execution": report.get("execution"),
        "requested": report.get("requested"),
        "returned": report.get("returned"),
        "files": report.get("files") or report.get("worker_files"),
    }


def _remove_hs_setting(input_data: dict) -> bool:
    configured = input_data.get("sdf_remove_hs")
    if configured is not None:
        return bool(configured)
    return not bool(input_data.get("write_filtered_sdf", False))


def _unique_non_null(values) -> list:
    output = []
    seen = set()
    for value in values:
        if value in (None, ""):
            continue
        marker = json.dumps(
            value, sort_keys=True, ensure_ascii=False, allow_nan=False
        )
        if marker not in seen:
            output.append(value)
            seen.add(marker)
    return output


def _number(value: Any, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _finite_number_or_none(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _integer(value: Any, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            value,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )


def main() -> None:
    from io_utils import parse_skill_args, print_skill_result

    input_data, output_dir = parse_skill_args("Drug Property Filter Skill")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

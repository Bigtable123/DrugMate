#!/usr/bin/env python3
"""Molecule Generator Skill - Standalone CLI.

Builds a compound library (actives + decoys) for virtual screening by querying
ChEMBL for bioactive compounds, generating or downloading decoys, standardizing
SMILES, and producing 3D conformers in a single SDF library.

Usage:
    python run.py --input '{"target_name": "CDK2", "activity_threshold_nm": 1000}' --output-dir ./output/
    python run.py --input '{"chembl_target_id": "CHEMBL301", "max_actives": 50}' --output-dir ./output/
    python run.py --input-file params.json --output-dir ./output/
"""
from __future__ import annotations

import json
import sys
import hashlib
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from deployment_config import deployment_path
from chembl_utils import (
    check_chembl_api_available,
    resolve_target_to_chembl_id,
    fetch_active_compounds,
)
from dude_utils import check_dude_available, download_dude_decoys
from decoy_gen_utils import check_rdkit_available, generate_property_matched_decoys
from pocket2mol_utils import (
    check_pocket2mol_status,
    prepare_pocket_pdb,
    run_pocket2mol,
    filter_by_blueprint,
)
from standardize_utils import (
    check_standardization_available,
    standardize_smiles_list,
    deduplicate_compounds,
)
from conformer_utils import check_conformer_gen_available, generate_conformers_sdf
from external_generator_utils import GENERATOR_IDS, run_registered_generator
from library_utils import assemble_library, validate_prepared_library
from io_utils import tool_enabled
from molecule_failures import MoleculeFailureTracker

CANDIDATE_SOURCE_TOOLS = (
    "chembl", "pocket2mol", "resgen", "diffsbdd", "targetdiff",
    "difflinker", "reinvent4", "lingo3dmol", "tamgen", "molcraft",
)
MOLGEN_TOOLS = (
    *CANDIDATE_SOURCE_TOOLS,
    "dude_decoys", "rdkit_decoys", "standardization", "conformer",
)
TOOL_DEFAULTS = {
    "chembl": True,
    "pocket2mol": True,
    "resgen": False,
    "diffsbdd": False,
    "targetdiff": False,
    "difflinker": False,
    "reinvent4": False,
    "lingo3dmol": False,
    "tamgen": False,
    "molcraft": False,
    "dude_decoys": True,
    "rdkit_decoys": False,
    "standardization": True,
    "conformer": True,
}


def run(input_data: dict, output_dir: Path) -> dict:
    """Persist isolated molecule failures and report shared execution faults."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    failure_file = output_dir / "molecule_failure_report.json"
    cohort_file = output_dir / "source_cohort_report.json"
    _write_json(cohort_file, {"status": "not_started", "cohorts": [], "sources": {}})

    def persist(tracker, **extra):
        _write_json(failure_file, {**tracker.report(), "status": "running", **extra})

    tracker = MoleculeFailureTracker(on_failure=persist)
    persist(tracker)
    try:
        result = _run(input_data, output_dir, tracker)
    except Exception as exc:
        result = {
            "success": False, "error": str(exc), "error_code": "molecule_tool_failed",
            "tool_failure": {"tool": "molecule_generation", "exception_type": type(exc).__name__,
                             "reason": str(exc)},
            "operational_errors": [{"tool": "molecule_generation", "reason": str(exc)}],
        }
    status = ("completed_with_molecule_errors" if tracker.failures else "complete")
    if not result.get("success"):
        status = "tool_failed" if result.get("tool_failure") else "failed"
    persist(tracker, status=status, tool_failure=result.get("tool_failure"),
            error=result.get("error"), released_summary=result.get("summary"))
    result.setdefault("files", {})["molecule_failure_report"] = str(failure_file)
    if cohort_file.exists():
        result["files"]["source_cohort_report"] = str(cohort_file)
    result["molecule_failure_summary"] = {
        "status": status, "failure_events": len(tracker.failures),
    }
    if result.get("success"):
        parameters = result.get("parameters", {})
        warnings = bool(tracker.failures or parameters.get("candidate_source_failures")
                        or parameters.get("preparation_warning_records")
                        or parameters.get("optimization_warnings"))
        result["outcome"] = "success_with_warnings" if warnings else "success"
    else:
        result.setdefault("outcome", "failed")
    result.setdefault("reason_codes", [result["error_code"]] if result.get("error_code") else [])
    result.setdefault("retryable", False)
    result.setdefault("next_action", "property_prefilter" if result.get("success") else "stop")
    return result


def _write_json(path: Path, data: dict):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str))
    temporary.replace(path)


def _run(input_data: dict, output_dir: Path, tracker: MoleculeFailureTracker) -> dict:
    """Main entry point for molecule library generation."""
    target_name = input_data.get("target_name")
    chembl_target_id = input_data.get("chembl_target_id")
    uniprot_id = input_data.get("uniprot_id")
    activity_type = input_data.get("activity_type", "IC50")
    activity_threshold_nm = input_data.get("activity_threshold_nm", 1000.0)
    max_actives = input_data.get("max_actives", 50)
    num_decoys = input_data.get("num_decoys", 200)
    seed = input_data.get("seed", 42)
    neutralize = bool(input_data.get("neutralize", False))
    # Exclusion constraints (for iterative rounds)
    exclude_scaffolds = input_data.get("exclude_scaffolds", [])
    max_heavy_atoms = input_data.get("max_heavy_atoms")
    max_mw = input_data.get("max_mw")

    # Interface-v2 tool switches. With no explicit override, the reviewed
    # defaults enable ChEMBL + Pocket2Mol and keep every optional generator off.
    switches = {
        tool_id: tool_enabled(input_data, tool_id, TOOL_DEFAULTS[tool_id])
        for tool_id in MOLGEN_TOOLS
    }
    if not any(switches[tool] for tool in CANDIDATE_SOURCE_TOOLS):
        return {"success": False, "outcome": "invalid_input", "error": "At least one candidate source must be enabled"}

    # library.sdf is the canonical Step 2 output and is written only by the
    # coordinate-preparation stage. A disabled conformer tool therefore makes
    # this step non-executable; fail before running candidate sources or RDKit.
    if not switches["conformer"]:
        return {
            "success": False,
            "outcome": "invalid_input",
            "error": (
                "molecule_generation.conformer must be enabled to produce the "
                "screening-ready library.sdf artifact"
            ),
        }

    if not switches["standardization"]:
        return {"success": False, "outcome": "invalid_input",
                "error": "molecule_generation.standardization must be enabled; bypass Stage 04 only for an already prepared injected library"}
    if not check_rdkit_available():
        return {"success": False, "outcome": "capability_gap", "error": "RDKit is required but not installed"}
    if not check_standardization_available():
        return {"success": False, "outcome": "capability_gap", "error": "RDKit MolStandardize is required for library preparation"}

    output_dir.mkdir(parents=True, exist_ok=True)

    steps_completed = []
    steps_skipped = []
    parameters = {
        "activity_type": activity_type,
        "activity_threshold_nm": activity_threshold_nm,
        "max_actives": max_actives,
        "num_decoys": num_decoys,
        "seed": seed,
        "neutralize": neutralize,
        "exclude_scaffolds": exclude_scaffolds,
        "max_heavy_atoms": max_heavy_atoms,
        "max_mw": max_mw,
        "microstate_preparation": "not_enumerated_in_stage04; docking_prep_may_expand",

    }

    actives = []
    target_info = {}
    source_failures = {}

    # ===== Step 1: Target Resolution (only for the ChEMBL route) =====
    if switches["chembl"]:
        raw_chembl_file = output_dir / "chembl_activities_raw.json"
        try:
            if not any([target_name, chembl_target_id, uniprot_id]):
                raise ValueError(
                    "Must provide target_name, chembl_target_id, or uniprot_id "
                    "for the selected ChEMBL source"
                )
            if chembl_target_id:
                target_info = {
                    "chembl_id": chembl_target_id,
                    "pref_name": target_name or "",
                    "organism": "",
                    "target_type": "",
                }
                if check_chembl_api_available():
                    try:
                        target_info = resolve_target_to_chembl_id(
                            chembl_target_id=chembl_target_id
                        )
                        steps_completed.append("target_resolution")
                    except Exception as exc:
                        steps_skipped.append({
                            "step": "target_resolution_validation",
                            "reason": str(exc),
                        })
                        steps_completed.append("target_resolution_from_input")
                else:
                    steps_completed.append("target_resolution_from_input")
            else:
                # The status endpoint can fail while target/activity reads
                # still work. Let the actual lookup use its bounded retries
                # instead of treating the advisory probe as a hard gate.
                target_info = resolve_target_to_chembl_id(
                    target_name=target_name, uniprot_id=uniprot_id
                )
                chembl_target_id = target_info["chembl_id"]
                steps_completed.append("target_resolution")

            # ===== Step 2: ChEMBL Activity Retrieval =====
            actives = fetch_active_compounds(
                chembl_target_id=chembl_target_id,
                activity_type=activity_type,
                threshold_nm=activity_threshold_nm,
                max_compounds=max_actives,
                raw_output_file=raw_chembl_file,
            )
            parameters["chembl_provenance"] = _chembl_provenance(
                raw_chembl_file
            )
            if not actives:
                raise RuntimeError(
                    f"No active compounds found for {chembl_target_id} with "
                    f"{activity_type} <= {activity_threshold_nm} nM"
                )
            steps_completed.append("chembl_retrieval")
            parameters["actives_retrieved"] = len(actives)
        except Exception as exc:
            actives = []
            if raw_chembl_file.is_file():
                try:
                    parameters["chembl_provenance"] = _chembl_provenance(
                        raw_chembl_file
                    )
                except Exception as provenance_exc:
                    parameters["chembl_provenance_error"] = str(provenance_exc)
            source_failures["chembl"] = str(exc)
            steps_skipped.append({
                "step": "chembl_retrieval",
                "reason": str(exc),
            })
    else:
        steps_skipped.append({"step": "chembl_retrieval", "reason": "disabled_by_config"})

    parameters["chembl_target_id"] = chembl_target_id
    parameters["target_info"] = target_info

    # ===== Step 2.5: Independently selectable molecule generators =====
    generated_compounds = []
    generation_reports = {}
    # pocket2mol on/off comes from the interface-v2 switch; the legacy
    # enable_pocket2mol flag is honored as a fallback for standalone calls only
    # when no tools block was supplied.
    enable_pocket2mol = switches["pocket2mol"]
    if "tools" not in input_data:
        enable_pocket2mol = enable_pocket2mol and bool(input_data.get("enable_pocket2mol", True))
    if enable_pocket2mol and input_data.get("pocket_pdb_file"):
        pocket2mol_compounds = _try_generate_pocket2mol(
            input_data=input_data,
            output_dir=output_dir,
            parameters=parameters,
            steps_completed=steps_completed,
            steps_skipped=steps_skipped,
            failure_tracker=tracker,
        )
        generated_compounds.extend(pocket2mol_compounds)
        generation_reports["pocket2mol"] = {
            "success": bool(pocket2mol_compounds),
            "num_valid": len(pocket2mol_compounds),
            **dict(parameters.get("pocket2mol") or {}),
        }
        if not pocket2mol_compounds:
            source_failures["pocket2mol"] = (
                generation_reports["pocket2mol"].get("error")
                or "produced no valid molecules"
            )
    elif not switches["pocket2mol"]:
        steps_skipped.append({"step": "pocket2mol", "reason": "disabled_by_config"})
    elif enable_pocket2mol:
        reason = "pocket_pdb_file not provided"
        source_failures["pocket2mol"] = reason
        steps_skipped.append({"step": "pocket2mol", "reason": reason})

    for generator_id in GENERATOR_IDS:
        if not switches[generator_id]:
            steps_skipped.append({
                "step": generator_id, "reason": "disabled_by_config",
            })
            continue
        try:
            compounds, report = run_registered_generator(
                generator_id, input_data, output_dir
            )
            generation_reports[generator_id] = report
            if compounds:
                generated_compounds.extend(compounds)
                steps_completed.append(f"{generator_id}_generation")
            else:
                reason = "produced no valid molecules"
                source_failures[generator_id] = reason
                steps_skipped.append({"step": generator_id, "reason": reason})
        except Exception as exc:
            source_failures[generator_id] = str(exc)
            generation_reports[generator_id] = {
                "tool": generator_id, "success": False, "error": str(exc),
            }
            steps_skipped.append({"step": generator_id, "reason": str(exc)})

    parameters["generation_by_source"] = generation_reports
    parameters["candidate_source_failures"] = source_failures
    if not actives and not generated_compounds:
        selected = [tool for tool in CANDIDATE_SOURCE_TOOLS if switches[tool]]
        _write_json(output_dir / "source_cohort_report.json", _source_cohort_report(
            [], [], switches, parameters, tracker, input_data,
        ))
        return {
            "success": False,
            "outcome": "terminal_no_candidates",
            "error": (
                "no candidate molecules were produced by selected sources "
                f"{selected}; failures={source_failures}"
            ),
        }

    # ===== Step 3: Decoy Generation =====
    decoys = []
    decoy_target_name = target_name or target_info.get("pref_name", "")

    if switches["dude_decoys"] and decoy_target_name and check_dude_available(decoy_target_name):
        dude_result = download_dude_decoys(
            target_name=decoy_target_name,
            output_dir=output_dir,
            max_decoys=num_decoys,
            seed=seed,
        )
        if dude_result["success"] and dude_result["decoys"]:
            decoys = dude_result["decoys"]
            steps_completed.append("dude_decoys")
            parameters["decoy_source"] = "dude"
            parameters["dude_target"] = dude_result.get("target_dude_name", "")
            parameters["dude_provenance"] = {
                "from_cache": dude_result.get("from_cache"),
                "cache_file": dude_result.get("cache_file"),
                "sampled_file": dude_result.get("file"),
                "cache_sha256": _sha256(Path(dude_result["cache_file"])) if dude_result.get("cache_file") else None,
                "sampled_sha256": _sha256(Path(dude_result["file"])) if dude_result.get("file") else None,
            }
        else:
            steps_skipped.append({"step": "dude_decoys", "reason": dude_result.get("error", "download failed")})
    elif not switches["dude_decoys"]:
        steps_skipped.append({"step": "dude_decoys", "reason": "disabled_by_config"})
    else:
        reason = "DUD-E not available for target" if decoy_target_name else "no target name for DUD-E lookup"
        steps_skipped.append({"step": "dude_decoys", "reason": reason})

    if not decoys and switches["rdkit_decoys"] and actives:
        active_smiles = [a["smiles"] for a in actives]
        decoy_result = generate_property_matched_decoys(
            active_smiles_list=active_smiles,
            num_decoys=num_decoys,
            seed=seed,
            failure_tracker=tracker,
        )
        if decoy_result["success"] and decoy_result["decoys"]:
            decoys = decoy_result["decoys"]
            steps_completed.append("rdkit_decoy_generation")
            parameters["decoy_source"] = "rdkit_property_matched"
            parameters["decoy_generation_params"] = decoy_result.get("generation_params", {})
        else:
            steps_skipped.append({"step": "rdkit_decoy_generation", "reason": decoy_result.get("error", "generation failed")})
    elif not decoys and not switches["rdkit_decoys"]:
        steps_skipped.append({"step": "rdkit_decoy_generation", "reason": "disabled_by_config"})
    elif not decoys and not actives:
        steps_skipped.append({
            "step": "rdkit_decoy_generation",
            "reason": "requires ChEMBL actives; generated candidates are not decoy templates",
        })

    if not decoys:
        steps_skipped.append({
            "step": "decoy_generation",
            "reason": "all decoy methods failed or disabled; library will contain candidates only",
        })

    # ===== Step 4: SMILES Standardization =====
    _mark_compound_sources(actives, decoys, generated_compounds, input_data=input_data)
    source_records = [dict(c) for c in actives + decoys + generated_compounds]
    parameters["input_source_records_by_role"] = {
        "active": len(actives), "validation_decoy": len(decoys),
        "generated": len(generated_compounds),
    }
    all_compounds = actives + decoys + generated_compounds

    standardized, failed = standardize_smiles_list(
        all_compounds, smiles_key="smiles", neutralize=neutralize, failure_tracker=tracker,
    )
    parameters["standardization_failed_count"] = len(failed)
    actives = [c for c in standardized if _compound_role(c) == "active"]
    decoys = [c for c in standardized if _compound_role(c) == "decoy"]
    generated_compounds = [c for c in standardized if _compound_role(c) == "generated"]
    parameters["standardized_by_cohort"] = _counts_by_cohort(standardized)
    steps_completed.append("standardization")

    actives, decoys, generated_compounds, deduplication = _deduplicate_inputs(
        actives, decoys, generated_compounds, failure_tracker=tracker,
    )
    parameters["deduplication"] = deduplication
    parameters["unique_by_cohort"] = _counts_by_cohort(actives + decoys + generated_compounds)
    cohort_file = output_dir / "source_cohort_report.json"
    _write_json(cohort_file, _source_cohort_report(source_records, [], switches, parameters, tracker, input_data))

    # Iterative-round exclusions apply to every source that can be nominated:
    # actives and generated candidates.  Validation decoys are deliberately
    # exempt -- they are labeled negatives for the retrospective screening
    # benchmark (AUROC/EF1%), never proposals, and dropping the ones that happen
    # to share an already-tried scaffold would shrink the negative set and
    # distort the very gate it feeds (``min_validation_decoys``).  The exemption
    # is published so a reviewer sees the decoy count was left untouched on
    # purpose rather than by omission.
    if exclude_scaffolds or max_heavy_atoms or max_mw:
        before_by_source = _counts_by_source(actives + generated_compounds)
        actives, excluded_actives = _apply_exclusion_filters(
            actives, exclude_scaffolds, max_heavy_atoms, max_mw
        )
        generated_compounds, excluded_generated = _apply_exclusion_filters(
            generated_compounds, exclude_scaffolds, max_heavy_atoms, max_mw
        )
        excluded_count = excluded_actives + excluded_generated
        if excluded_count:
            steps_completed.append("exclusion_filtering")
        parameters["exclusion_filtering"] = {
            "excluded_total": excluded_count,
            "before_by_source": before_by_source,
            "after_by_source": _counts_by_source(actives + generated_compounds),
            "exempt_validation_decoys": len(decoys),
            "exempt_reason": (
                "validation decoys are labeled negatives for the screening "
                "benchmark, not nomination candidates"
            ),
        }

    if not actives and not generated_compounds:
        return {
            "success": False,
            "outcome": "terminal_no_candidates",
            "error": "no candidate compounds remained after standardization/exclusion filtering",
        }

    # Assign stable IDs before conformer generation. Metadata files are written
    # only after conformer generation so their counts match library.sdf.
    _assign_compound_ids(actives, decoys, generated_compounds)
    for compound in actives + decoys + generated_compounds:
        compound["record_status"] = "ready"
        compound["molecule_failure_report"] = "molecule_failure_report.json"
        compound["source_cohort_report"] = "source_cohort_report.json"

    parameters["unique_records_by_role"] = {
        "active": len(actives), "validation_decoy": len(decoys),
        "generated": len(generated_compounds),
    }

    # ===== Step 5: 3D Conformer Generation =====
    all_for_conformers = actives + decoys + generated_compounds
    library_sdf = output_dir / "library.sdf"
    coordinate_sdf = output_dir / "conformer_candidates.sdf"

    if check_conformer_gen_available():
        conf_result = generate_conformers_sdf(
            compounds=all_for_conformers,
            output_sdf=coordinate_sdf,
            smiles_key="smiles_standardized",
            optimize=True,
            random_seed=seed,
            failure_tracker=tracker,
        )
        conf_result["output_sdf"] = str(library_sdf) if conf_result["success"] else None
        _write_json(output_dir / "conformer_report.json", conf_result)
        _write_json(cohort_file, _source_cohort_report(
            source_records, [c for c in all_for_conformers if c.get("conformer_written")],
            switches, parameters, tracker, input_data,
        ))
        if conf_result["success"]:
            steps_completed.append("conformer_generation")
            parameters["coordinate_preparation"] = (
                "preserve_compatible_generator_3d_else_ETKDGv3_MMFF_or_UFF"
            )
            parameters["conformers_written"] = conf_result["total_written"]
            parameters["conformers_failed"] = conf_result["failed_embedding"]
            parameters["raw_poses_preserved"] = conf_result.get("raw_poses_preserved", 0)
            parameters["raw_poses_rejected"] = conf_result.get("raw_poses_rejected", [])
            parameters["embedded_conformers"] = conf_result.get("embedded_conformers", 0)
            parameters["optimization_warnings"] = conf_result.get("optimization_warnings", [])
            actives = [compound for compound in actives if compound.get("conformer_written")]
            decoys = [compound for compound in decoys if compound.get("conformer_written")]
            generated_compounds = [
                compound for compound in generated_compounds if compound.get("conformer_written")
            ]
        else:
            return {"success": False, "outcome": "terminal_no_candidates", "error": "all conformer embeddings/optimizations failed"}
    else:
        return {"success": False, "error": "RDKit AllChem is required for screening-ready 3D conformers"}

    if not actives and not generated_compounds:
        return {
            "success": False,
            "outcome": "terminal_no_candidates",
            "error": "no candidate compounds remained after 3D coordinate preparation",
        }

    # Validate the materialized SDF once, after isolated molecule failures.
    # No pH-state enumeration, whole-parent quarantine, or ensemble replay.
    retained = actives + decoys + generated_compounds
    validate_prepared_library(coordinate_sdf, retained)
    parameters["rdkit_version"] = _rdkit_version()
    parameters["molecule_failure_report"] = str(output_dir / "molecule_failure_report.json")
    parameters["molecule_failures"] = tracker.report()
    source_cohorts = _source_cohort_report(
        source_records, retained, switches, parameters, tracker, input_data,
    )
    _write_json(cohort_file, source_cohorts)
    parameters["source_cohort_report"] = str(cohort_file)
    parameters["conformers_written"] = len(retained)
    parameters["raw_poses_preserved"] = sum(c.get("raw_pose_preserved", False) for c in retained)
    parameters["embedded_conformers"] = len(retained) - parameters["raw_poses_preserved"]
    parameters["preparation_warning_records"] = sum(bool(c.get("preparation_warnings_json")) for c in retained)
    # Stage the complete release before touching any previous SDF or metadata.
    # A failed replacement restores files already published by this attempt.
    with tempfile.TemporaryDirectory(prefix=".library-release-", dir=output_dir) as staging:
        staging = Path(staging)
        assembly_result = assemble_library(
            actives=actives,
            generated=generated_compounds,
            decoys=decoys,
            output_dir=staging,
            target_info=target_info,
            steps_completed=steps_completed,
            steps_skipped=steps_skipped,
            parameters=parameters,
        )
        coordinate_sdf.replace(staging / "library.sdf")
        _publish_library(staging, output_dir)
    assembly_result["manifest_csv"] = str(output_dir / "library_manifest.csv")
    assembly_result["report_json"] = str(output_dir / "library_report.json")

    # ===== Final Result =====
    result = {
        "success": True,
        "files": {
            "library_sdf": str(library_sdf) if library_sdf.exists() else None,
            # Canonical artifact key consumed downstream (structure-evaluator,
            # molecule-modifier) and by the orchestrator's artifact wiring.
            "library_manifest": assembly_result.get("manifest_csv"),
            "source_cohort_report": str(cohort_file),
            "report": assembly_result.get("report_json"),
            "conformer_report": str(output_dir / "conformer_report.json"),
        },
        "summary": {
            "num_actives": len(actives),
            "num_decoys": len(decoys),
            "num_generated": len(generated_compounds),
            "num_pocket2mol": sum(
                compound.get("source") == "pocket2mol"
                for compound in generated_compounds
            ),
            # Primary-source counts describe the retained representative.
            "generated_by_source": _counts_by_source(generated_compounds),
            # Contribution counts also credit generators collapsed into the
            # representative during cross-source SMILES deduplication.
            "generated_contribution_by_source": _contribution_counts_by_source(
                generated_compounds
            ),
            "total_compounds": len(actives) + len(decoys) + len(generated_compounds),
            "target_chembl_id": chembl_target_id,
            "target_name": target_info.get("pref_name", target_name or ""),
            "activity_type": activity_type,
            "activity_threshold_nm": activity_threshold_nm,
        },
        "steps_completed": steps_completed,
        "steps_skipped": steps_skipped,
        "parameters": parameters,
    }
    return result


def _publish_library(staging: Path, output_dir: Path) -> None:
    names = ("library.sdf", "library_manifest.csv", "library_report.json")
    for name in names:
        previous = output_dir / name
        if previous.exists():
            shutil.copy2(previous, staging / (name + ".previous"))
    replaced = []
    try:
        for name in names:
            (staging / name).replace(output_dir / name)
            replaced.append(name)
    except OSError:
        for name in reversed(replaced):
            backup = staging / (name + ".previous")
            if backup.exists():
                backup.replace(output_dir / name)
            else:
                (output_dir / name).unlink()
        raise


def _assign_compound_ids(actives: list[dict], decoys: list[dict], generated: list[dict] | None = None) -> None:
    for i, compound in enumerate(actives):
        compound["compound_id"] = f"ACT_{i + 1:04d}"
        compound["source"] = "active"
        compound["candidate_role"] = "active"
    for i, compound in enumerate(decoys):
        compound["compound_id"] = f"DEC_{i + 1:04d}"
        compound["candidate_role"] = "validation_decoy"
    for i, compound in enumerate(generated or []):
        compound["compound_id"] = f"GEN_{i + 1:04d}"
        compound["candidate_role"] = "generated"
    for compound in actives + decoys + (generated or []):
        compound["microstate_id"] = None
        compound["conformer_id"] = f"{compound['compound_id']}:conf:0001"


def _mark_compound_sources(
    actives: list[dict],
    decoys: list[dict],
    generated: list[dict] | None = None,
    *, input_data: dict | None = None,
) -> None:
    for compound in actives:
        compound["source"] = "active"
        compound["candidate_role"] = "active"
    for compound in decoys:
        compound.setdefault("source", "decoy")
        compound["candidate_role"] = "validation_decoy"
    for compound in generated or []:
        compound.setdefault("source", "generated")
        compound["candidate_role"] = "generated"
    for role, compounds in (("active", actives), ("decoy", decoys),
                             ("generated", generated or [])):
        for index, compound in enumerate(compounds, 1):
            compound.setdefault("source_record_id", f"{role}:{index:06d}")
            compound.setdefault("source_tool", "chembl" if role == "active" else {
                "dude": "dude_decoys", "rdkit_property_matched": "rdkit_decoys",
            }.get(compound["source"], compound["source"]))
            retrospective = role == "decoy" or (role == "active" and bool(decoys))
            compound.setdefault("campaign_mode", "retrospective" if retrospective else "prospective")
            if retrospective:
                dataset = ((input_data or {}).get("benchmark_dataset")
                           or ("DUD-E" if any(c.get("source") == "dude" for c in decoys) else "RDKit-property-matched"))
                compound.setdefault("benchmark_dataset", dataset)
                compound.setdefault("benchmark_split", (input_data or {}).get("benchmark_split") or "unspecified")
                compound.setdefault("cohort_id", f"retrospective:{dataset}")
            else:
                compound.setdefault("cohort_id", (input_data or {}).get("cohort_id") or "prospective")
            # A campaign may combine complementary generators while retaining
            # each generator's provenance. Only an explicit scientific cohort
            # boundary prevents exact-structure deduplication across sources.
            compound.setdefault("source_cohort_id", f"{compound['cohort_id']}:{compound['source_tool']}")


def _counts_by_cohort(compounds):
    counts = {}
    for compound in compounds:
        cohort = compound["cohort_id"]
        counts[cohort] = counts.get(cohort, 0) + 1
    return counts


def _source_cohort_report(inputs, retained, switches, parameters, tracker, config):
    """Index source records, failed molecules and zero-yield source cohorts."""
    cohorts = {}
    for record in inputs:
        cohort = cohorts.setdefault(record["cohort_id"], {
            key: record.get(key) for key in ("cohort_id", "campaign_mode", "benchmark_dataset", "benchmark_split")
        })
        cohort.setdefault("source_tools", set()).add(record["source_tool"])
        cohort.setdefault("source_record_ids", []).append(record["source_record_id"])
    ready = _counts_by_cohort(retained)
    failed_ids = {row.get("source_record_id") for row in tracker.failures}
    for cohort_id, cohort in cohorts.items():
        cohort["source_tools"] = sorted(cohort["source_tools"])
        cohort["counts"] = {
            "source_records": len(cohort["source_record_ids"]),
            "standardized": parameters.get("standardized_by_cohort", {}).get(cohort_id, 0),
            "unique": parameters.get("unique_by_cohort", {}).get(cohort_id, 0),
            "3d_ready": ready.get(cohort_id, 0),
        }
        cohort["failed_source_record_ids"] = [sid for sid in cohort["source_record_ids"] if sid in failed_ids]
        cohort["terminal_state"] = ("terminal_no_candidates" if not ready.get(cohort_id)
                                    else "success_with_warnings" if cohort["failed_source_record_ids"] else "success")
    realized_tools = {c["source_tool"] for c in inputs}
    sources = {}
    for tool in (*CANDIDATE_SOURCE_TOOLS, "dude_decoys", "rdkit_decoys"):
        if not switches[tool]:
            sources[tool] = {"terminal_state": "disabled_by_config"}
            continue
        records = [c for c in inputs if c["source_tool"] == tool]
        valid = [c for c in retained if tool in {
            parent.get("source_tool") for parent in c.get("parent_input_records", [c])
        }]
        error = parameters.get("candidate_source_failures", {}).get(tool)
        sources[tool] = {
            "requested": (config.get("max_actives", 50) if tool == "chembl"
                          else config.get("num_decoys", 200) if tool.endswith("decoys")
                          else parameters.get("generation_by_source", {}).get(tool, {}).get("num_samples_requested")
                          or config.get(f"{tool}_num_samples")
                          or config.get("tools", {}).get(tool, {}).get("params", {}).get("num_samples")),
            "source_records": len(records), "3d_ready": len(valid),
            "representative_3d_ready": sum(c["source_tool"] == tool for c in retained),
            "source_cohort_ids": sorted({c["source_cohort_id"] for c in records}) or [
                f"{config.get('cohort_id') or 'prospective'}:{tool}"
            ],
            "terminal_state": "failed" if error else "success" if valid else "terminal_no_candidates",
            "error": error,
            "provenance": parameters.get("generation_by_source", {}).get(tool, {}),
        }
        if tool not in realized_tools and tool in CANDIDATE_SOURCE_TOOLS:
            cohort_id = config.get("cohort_id") or "prospective"
            cohort = cohorts.setdefault(cohort_id, {"cohort_id": cohort_id, "campaign_mode": "prospective",
                "source_tools": [], "source_record_ids": [], "counts": {"source_records": 0, "standardized": 0, "unique": 0, "3d_ready": 0},
                "terminal_state": sources[tool]["terminal_state"]})
            cohort["source_tools"] = sorted(set(cohort["source_tools"]) | {tool})
            cohort.setdefault("failed_sources", {})[tool] = error or "No source records were produced"
            if cohort["counts"]["3d_ready"]:
                cohort["terminal_state"] = "success_with_warnings"
    return {"cohorts": list(cohorts.values()), "sources": sources,
            "molecule_failure_report": "molecule_failure_report.json",
            "microstate_preparation": "Stage 04 does not enumerate pH-dependent microstates; docking_prep may expand justified states"}


def _compound_role(compound: dict) -> str:
    role = compound.get("candidate_role")
    if role in ("active", "generated", "validation_decoy"):
        return "decoy" if role == "validation_decoy" else role
    source = compound.get("source")
    if source == "active" or compound.get("activity_nm") is not None:
        return "active"
    if source in CANDIDATE_SOURCE_TOOLS and source != "chembl":
        return "generated"
    if source in ("dude", "rdkit_property_matched", "decoy"):
        return "decoy"
    return "unknown"


def _deduplicate_inputs(actives, decoys, generated, *, failure_tracker=None):
    """Deduplicate exact standardized structures within their scientific cohort.

    Charge, tautomer and specified stereochemistry remain part of identity.
    A control that duplicates a candidate is excluded from the benchmark;
    candidate records remain unchanged and the control is never relabeled.
    """
    all_records = actives + decoys + generated
    candidates_by_smiles = {}
    for record in actives + generated:
        candidates_by_smiles.setdefault(record["smiles_standardized"], []).append(record)
    outputs, merged_groups = [], []
    for compounds, generated_role in ((actives, False), (decoys, False), (generated, True)):
        groups = {}
        for compound in compounds:
            smiles = compound["smiles_standardized"]
            if compound["candidate_role"] == "validation_decoy" and smiles in candidates_by_smiles:
                if failure_tracker:
                    failure_tracker.record(
                        compound, "standardization.role_identity", "candidate_control_identity_conflict",
                        "Benchmark control duplicates an existing candidate; excluded from the benchmark without changing candidate records",
                        details={"matched_candidate_source_record_ids": sorted({
                            c["source_record_id"] for c in candidates_by_smiles[smiles] if c.get("source_record_id")
                        })},
                    )
                continue
            key = (compound.get("cohort_id", ""), smiles)
            groups.setdefault(key, []).append(compound)
        result = []
        for entries in groups.values():
            representative = (_deduplicate_generated(entries)[0] if generated_role
                              else deduplicate_compounds(entries)[0])
            representative["parent_input_records"] = [
                {key: entry[key] for key in (
                    "source_record_id", "cohort_id", "source_cohort_id", "source_tool", "smiles", "smiles_standardized",
                    "source", "candidate_role", "chembl_id", "chembl_ids", "name",
                    "source_file", "source_mol_index",
                ) if key in entry} for entry in entries
            ]
            if len(entries) > 1:
                merged_groups.append({"source_record_ids": [c.get("source_record_id") for c in entries]})
            result.append(representative)
        outputs.append(result)
    retained = sum(map(len, outputs))
    conflicting = sum(c["smiles_standardized"] in candidates_by_smiles for c in decoys)
    return (*outputs, {
        "identity_method": "exact_standardized_isomeric_smiles_within_cohort",
        "input_records": len(all_records), "retained_records": retained,
        "merged_records": len(all_records) - retained - conflicting,
        "role_conflict_records": conflicting, "merged_groups": merged_groups,
    })


def _deduplicate_generated(
    compounds: list[dict], *, identity_key: str = "smiles_standardized",
) -> list[dict]:
    """Deduplicate across generators while retaining multi-source provenance."""
    grouped = {}
    for compound in compounds:
        smiles = compound.get(identity_key) or compound.get("smiles")
        if smiles:
            grouped.setdefault(smiles, []).append(compound)

    deduplicated = []
    for entries in grouped.values():
        representative = next(
            (entry for entry in entries if entry.get("raw_pose_molblock")),
            entries[0],
        ).copy()
        sources = sorted({
            str(entry.get("source"))
            for entry in entries
            if entry.get("source")
        })
        representative["source_tools"] = sources
        representative["also_generated_by"] = [
            source for source in sources if source != representative.get("source")
        ]
        representative["generator_entry_count"] = len(entries)
        deduplicated.append(representative)
    return deduplicated


def _counts_by_source(compounds: list[dict]) -> dict[str, int]:
    counts = {}
    for compound in compounds:
        source = str(compound.get("source") or "unknown")
        counts[source] = counts.get(source, 0) + 1
    return dict(sorted(counts.items()))


def _contribution_counts_by_source(compounds: list[dict]) -> dict[str, int]:
    counts = {}
    for compound in compounds:
        sources = compound.get("source_tools") or [compound.get("source") or "unknown"]
        for source in set(str(value) for value in sources if value):
            counts[source] = counts.get(source, 0) + 1
    return dict(sorted(counts.items()))


def _try_generate_pocket2mol(
    input_data: dict,
    output_dir: Path,
    parameters: dict,
    steps_completed: list,
    steps_skipped: list,
    failure_tracker=None,
) -> list[dict]:
    pocket_pdb_file = Path(input_data["pocket_pdb_file"])
    pocket_center = input_data.get("pocket_center")
    blueprint_file = input_data.get("interaction_blueprint_file")
    p2m_config = _pocket2mol_config(input_data)
    pocket2mol_dir = p2m_config["pocket2mol_dir"]
    conda_env = p2m_config["conda_env"]
    checkpoint_path = p2m_config["checkpoint"]
    gpu_device = int(p2m_config["gpu_device"])
    taskset_cpu = p2m_config["taskset_cpu"]
    conda_executable = p2m_config.get("conda_executable")
    checkpoint_sha256 = p2m_config.get("asset_sha256")

    if not pocket_center:
        reason = "pocket_center not provided"
        parameters["pocket2mol"] = {"success": False, "error": reason}
        steps_skipped.append({"step": "pocket2mol", "reason": reason})
        return []
    if not pocket_pdb_file.exists():
        reason = f"pocket_pdb_file not found: {pocket_pdb_file}"
        parameters["pocket2mol"] = {"success": False, "error": reason}
        steps_skipped.append({"step": "pocket2mol", "reason": reason})
        return []

    p2m_status = check_pocket2mol_status(
        pocket2mol_dir=pocket2mol_dir,
        conda_env=conda_env,
        checkpoint_path=checkpoint_path,
        conda_executable=conda_executable,
        checkpoint_sha256=checkpoint_sha256,
    )
    parameters["pocket2mol_status"] = p2m_status
    if not p2m_status.get("available"):
        reason = p2m_status.get("reason", "Pocket2Mol not ready")
        parameters["pocket2mol"] = {
            "success": False,
            "error": reason,
            "checkpoint_sha256": p2m_status.get("checkpoint_sha256"),
            "checkpoint_checksum_verified": p2m_status.get(
                "checkpoint_checksum_verified"
            ),
        }
        steps_skipped.append({"step": "pocket2mol", "reason": reason})
        return []

    try:
        blueprint_residues = _extract_blueprint_residues(blueprint_file)
        trimmed_pdb = prepare_pocket_pdb(
            pocket_pdb_file,
            blueprint_residues,
            chain=input_data.get("chain", "A"),
            output_dir=output_dir / "pocket2mol_work",
        )
        num_samples = int(input_data.get("pocket2mol_num_samples", 100))
        raw_generated, generation_report = run_pocket2mol(
            pocket_pdb=trimmed_pdb,
            center=pocket_center,
            num_samples=num_samples,
            pocket2mol_dir=pocket2mol_dir,
            conda_env=conda_env,
            checkpoint_path=checkpoint_path,
            gpu_device=gpu_device,
            taskset_cpu=taskset_cpu,
            timeout_seconds=int(input_data.get("pocket2mol_timeout_seconds", 86400)),
            output_dir=output_dir / "pocket2mol_work",
            conda_executable=conda_executable,
            seed=int(input_data.get("seed", 42)),
            checkpoint_sha256=checkpoint_sha256,
            return_report=True,
            failure_tracker=failure_tracker,
        )
        pocket2mol_compounds, filter_diagnostics = filter_by_blueprint(
            raw_generated,
            blueprint_file=blueprint_file,
            pocket_center=pocket_center,
            return_diagnostics=True,
            failure_tracker=failure_tracker,
        )
        if pocket2mol_compounds:
            steps_completed.append("pocket2mol_generation")
        else:
            steps_skipped.append({
                "step": "pocket2mol",
                "reason": "Pocket2Mol produced no candidates that passed output QC/filtering",
            })
        parameters["pocket2mol"] = {
            "success": bool(pocket2mol_compounds),
            "num_samples_requested": num_samples,
            "raw_generated": len(raw_generated),
            "after_filter": len(pocket2mol_compounds),
            "pocket_trimmed": str(trimmed_pdb),
            "blueprint_residues_used": len(blueprint_residues),
            "pocket2mol_dir": pocket2mol_dir,
            "conda_env": conda_env,
            "checkpoint": checkpoint_path,
            "checkpoint_sha256": generation_report.get("asset_sha256"),
            "gpu_device": gpu_device,
            "taskset_cpu": taskset_cpu,
            "seed": generation_report.get("seed"),
            "work_dir": generation_report.get("work_dir"),
            "output_dir": generation_report.get("output_dir"),
            "execution": generation_report.get("execution"),
            "parse_diagnostics": generation_report.get("diagnostics"),
            "filter_diagnostics": filter_diagnostics,
        }
        return pocket2mol_compounds
    except Exception as e:
        reason = str(e)
        parameters["pocket2mol"] = {
            "success": False,
            "error": reason,
            "checkpoint": checkpoint_path,
            "checkpoint_sha256": p2m_status.get("checkpoint_sha256"),
            "checkpoint_checksum_verified": p2m_status.get(
                "checkpoint_checksum_verified"
            ),
        }
        steps_skipped.append({"step": "pocket2mol", "reason": reason})
        return []


def _extract_blueprint_residues(blueprint_file: str | None) -> list[str]:
    """Extract Required + Desired residue names from a blueprint JSON file."""
    if not blueprint_file:
        return []
    try:
        bp = json.loads(Path(blueprint_file).read_text())
        residues = []
        for item in bp.get("required", []) + bp.get("desired", []):
            res = item.get("residue", "")
            if res:
                residues.append(res)
        return residues
    except Exception:
        return []


def _pocket2mol_config(input_data: dict) -> dict:
    """Read local Pocket2Mol execution settings from top-level input fields."""
    tool_config = dict(
        ((input_data.get("tools") or {}).get("pocket2mol") or {})
    )
    runtime = dict(tool_config.get("runtime") or {})
    return {
        "conda_env": _config_value(input_data, "pocket2mol_conda_env", "Pocket2Mol"),
        "pocket2mol_dir": _config_value(input_data, "pocket2mol_dir", str(deployment_path("${SBDD_POCKET2MOL_SOURCE_DIR}"))),
        "checkpoint": _config_value(input_data, "pocket2mol_checkpoint", str(deployment_path("${SBDD_POCKET2MOL_SOURCE_DIR}/pretrained_Pocket2Mol.pt"))),
        "asset_sha256": runtime.get("asset_sha256"),
        "gpu_device": _config_value(input_data, "pocket2mol_gpu_device", 0),
        "taskset_cpu": _config_value(input_data, "pocket2mol_taskset_cpu", 0),
        "conda_executable": _config_value(input_data, "pocket2mol_conda_executable", None),
    }


def _config_value(input_data: dict, key: str, default):
    value = input_data.get(key)
    if value is not None:
        return value
    return default


def _sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _chembl_provenance(raw_capture: Path) -> dict:
    """Read the explicit ChEMBL retrieval-completeness audit fields."""
    payload = json.loads(raw_capture.read_text(encoding="utf-8"))
    retrieval = payload.get("retrieval") or {}
    return {
        "raw_response_file": str(raw_capture),
        "sha256": _sha256(raw_capture),
        "retrieved_at": payload.get("retrieved_at"),
        "expected_records": retrieval.get("expected_records"),
        "retrieved_records": retrieval.get("retrieved_records"),
        "accepted_records": retrieval.get("accepted_records"),
        "flagged_data_validity_records": retrieval.get(
            "flagged_data_validity_records"
        ),
        "flagged_potential_duplicate_records": retrieval.get(
            "flagged_potential_duplicate_records"
        ),
        "quarantined_quality_records": retrieval.get(
            "quarantined_quality_records"
        ),
        "eligible_compound_groups": retrieval.get("eligible_compound_groups"),
        "quarantined_compound_groups": retrieval.get(
            "quarantined_compound_groups"
        ),
        "excluded_data_validity_records": retrieval.get(
            "excluded_data_validity_records"
        ),
        "excluded_potential_duplicate_records": retrieval.get(
            "excluded_potential_duplicate_records"
        ),
        "excluded_query_mismatch_records": retrieval.get(
            "excluded_query_mismatch_records"
        ),
        "pages_retrieved": retrieval.get("pages_retrieved"),
        "complete": retrieval.get("complete"),
        "truncated": retrieval.get("truncated"),
        "next": retrieval.get("next"),
        "error": payload.get("error"),
    }


def _apply_exclusion_filters(
    compounds: list[dict],
    exclude_scaffolds: list[str],
    max_heavy_atoms: int | None,
    max_mw: float | None,
) -> tuple[list[dict], int]:
    """Filter compounds by scaffold exclusion and size constraints.

    Returns (filtered_compounds, num_excluded).
    """
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem.Descriptors import ExactMolWt

    exclude_set = set()
    for smi in exclude_scaffolds:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            exclude_set.add(Chem.MolToSmiles(mol))

    original_count = len(compounds)
    filtered = []
    for comp in compounds:
        smiles = comp.get("smiles_standardized") or comp.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            # A compound we cannot parse cannot be verified against the
            # constraints. Dropping it is the safe choice for exclusion filters
            # (e.g. scaffold-hop rounds must not silently re-admit excluded or
            # oversized matter just because RDKit failed to parse it).
            continue
        # Scaffold exclusion
        if exclude_set:
            try:
                core = MurckoScaffold.GetScaffoldForMol(mol)
                core_smi = Chem.MolToSmiles(core)
                if core_smi in exclude_set:
                    continue
            except Exception:
                pass
        # Heavy atom limit (guard against a 0/None sentinel meaning "no limit")
        if max_heavy_atoms and mol.GetNumHeavyAtoms() > max_heavy_atoms:
            continue
        # Molecular weight limit
        if max_mw and ExactMolWt(mol) > max_mw:
            continue
        filtered.append(comp)

    return filtered, original_count - len(filtered)


def _rdkit_version() -> str | None:
    try:
        from rdkit import rdBase
        return rdBase.rdkitVersion
    except Exception:
        return None


def main():
    from io_utils import parse_skill_args, print_skill_result
    input_data, output_dir = parse_skill_args("Molecule Generator Skill")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

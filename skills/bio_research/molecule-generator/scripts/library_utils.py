"""Library assembly utilities.

Merges actives and decoys into final output files:
  - library.sdf (3D structures with tags)
  - library_manifest.csv (compound table with descriptors)
  - library_report.json (run summary and provenance)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

try:
    from pose_identity import LIGAND_LINEAGE_PROPERTIES
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
    from pose_identity import LIGAND_LINEAGE_PROPERTIES


PREPARATION_PROPERTIES = (
    *LIGAND_LINEAGE_PROPERTIES, "record_status",
)


def assemble_library(
    actives: list[dict],
    generated: list[dict],
    decoys: list[dict],
    output_dir: Path,
    target_info: dict,
    steps_completed: list[str],
    steps_skipped: list[dict],
    parameters: dict,
) -> dict:
    """Assign IDs and write manifest CSV + report JSON.

    The SDF is written separately by conformer_utils. This function handles
    the metadata files and ID assignment.

    Returns paths to generated files and summary stats.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, comp in enumerate(actives):
        comp.setdefault("compound_id", f"ACT_{i + 1:04d}")
        comp["source"] = "active"
        comp["candidate_role"] = "active"

    for i, comp in enumerate(generated):
        comp.setdefault("compound_id", f"GEN_{i + 1:04d}")
        comp.setdefault("source", "generated")
        comp["candidate_role"] = "generated"

    for i, comp in enumerate(decoys):
        comp.setdefault("compound_id", f"DEC_{i + 1:04d}")
        if "source" not in comp:
            comp["source"] = "decoy"
        comp["candidate_role"] = "validation_decoy"

    all_compounds = actives + generated + decoys

    manifest_path = output_dir / "library_manifest.csv"
    temporary_manifest = manifest_path.with_suffix(".csv.tmp")
    _write_manifest_csv(all_compounds, temporary_manifest)
    temporary_manifest.replace(manifest_path)

    report_path = output_dir / "library_report.json"
    report = _build_report(
        target_info=target_info,
        actives=actives,
        generated=generated,
        decoys=decoys,
        steps_completed=steps_completed,
        steps_skipped=steps_skipped,
        parameters=parameters,
    )
    temporary_report = report_path.with_suffix(".json.tmp")
    with open(temporary_report, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)
    temporary_report.replace(report_path)

    return {
        "manifest_csv": str(manifest_path),
        "report_json": str(report_path),
        "num_actives": len(actives),
        "num_generated": len(generated),
        "num_decoys": len(decoys),
        "total": len(all_compounds),
    }


def _write_manifest_csv(compounds: list[dict], output_path: Path):
    """Write CSV manifest with compound metadata and descriptors."""
    fieldnames = list(dict.fromkeys([
        "compound_id", "smiles", "smiles_original", "source", "candidate_role",
        *PREPARATION_PROPERTIES,
        "source_tools", "also_generated_by", "parent_input_records_json",
        "activity_nM", "standard_relation", "activity_lower_bound_nM",
        "activity_upper_bound_nM", "activity_lower_bound_inclusive",
        "activity_upper_bound_inclusive", "activity_is_censored",
        "activity_aggregation_method", "activity_measurement_count",
        "activity_aggregation_measurement_count",
        "activity_quarantined_measurement_count",
        "activity_data_validity_quarantined_count",
        "activity_potential_duplicate_quarantined_count",
        "activity_exact_measurement_count", "activity_censored_measurement_count",
        "activity_assay_heterogeneity", "activity_assay_types",
        "activity_bao_formats", "activity_assay_strata_json",
        "activity_document_chembl_ids", "activity_measurements_json",
        "assay_chembl_ids", "activity_eligible_assay_chembl_ids",
        "activity_quarantined_assay_chembl_ids", "chembl_ids",
        "chembl_id", "MW", "LogP", "HBD", "HBA",
        "RotBonds", "TPSA", "coordinate_origin", "raw_pose_preserved",
        "source_pose_file", "source_pose_index",
    ]))

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for comp in compounds:
            descriptors = _get_descriptors(comp)
            row = {
                "compound_id": comp.get("compound_id", ""),
                "smiles": comp.get("smiles_standardized", comp.get("smiles", "")),
                "smiles_original": comp.get("smiles", ""),
                "source": comp.get("source", ""),
                "candidate_role": comp.get("candidate_role", ""),
                **{
                    field: _manifest_value(comp.get(field))
                    for field in PREPARATION_PROPERTIES
                },
                "source_tools": ",".join(comp.get("source_tools") or []),
                "also_generated_by": ",".join(comp.get("also_generated_by") or []),
                "parent_input_records_json": json.dumps(
                    comp.get("parent_input_records", []),
                    sort_keys=True, separators=(",", ":"),
                ),
                "activity_nM": comp.get("activity_nm", ""),
                "standard_relation": comp.get("standard_relation", ""),
                "activity_lower_bound_nM": comp.get(
                    "activity_lower_bound_nm", ""
                ),
                "activity_upper_bound_nM": comp.get(
                    "activity_upper_bound_nm", ""
                ),
                "activity_lower_bound_inclusive": comp.get(
                    "activity_lower_bound_inclusive", ""
                ),
                "activity_upper_bound_inclusive": comp.get(
                    "activity_upper_bound_inclusive", ""
                ),
                "activity_is_censored": comp.get("activity_is_censored", ""),
                "activity_aggregation_method": comp.get(
                    "activity_aggregation_method", ""
                ),
                "activity_measurement_count": comp.get(
                    "activity_measurement_count", ""
                ),
                "activity_aggregation_measurement_count": comp.get(
                    "activity_aggregation_measurement_count", ""
                ),
                "activity_quarantined_measurement_count": comp.get(
                    "activity_quarantined_measurement_count", ""
                ),
                "activity_data_validity_quarantined_count": comp.get(
                    "activity_data_validity_quarantined_count", ""
                ),
                "activity_potential_duplicate_quarantined_count": comp.get(
                    "activity_potential_duplicate_quarantined_count", ""
                ),
                "activity_exact_measurement_count": comp.get(
                    "activity_exact_measurement_count", ""
                ),
                "activity_censored_measurement_count": comp.get(
                    "activity_censored_measurement_count", ""
                ),
                "activity_assay_heterogeneity": comp.get(
                    "activity_assay_heterogeneity", ""
                ),
                "activity_assay_types": ",".join(
                    str(value) for value in comp.get("activity_assay_types") or []
                ),
                "activity_bao_formats": ",".join(
                    str(value) for value in comp.get("activity_bao_formats") or []
                ),
                "activity_assay_strata_json": (
                    json.dumps(
                        comp.get("activity_assay_strata"),
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    if comp.get("activity_assay_strata") is not None else ""
                ),
                "activity_document_chembl_ids": ",".join(
                    str(value)
                    for value in comp.get("activity_document_chembl_ids") or []
                ),
                "activity_measurements_json": (
                    json.dumps(
                        comp.get("activity_measurements"),
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    if comp.get("activity_measurements") is not None else ""
                ),
                "assay_chembl_ids": ",".join(
                    str(value) for value in comp.get("assay_chembl_ids") or []
                ),
                "activity_eligible_assay_chembl_ids": ",".join(
                    str(value)
                    for value in comp.get("activity_eligible_assay_chembl_ids") or []
                ),
                "activity_quarantined_assay_chembl_ids": ",".join(
                    str(value)
                    for value in comp.get("activity_quarantined_assay_chembl_ids") or []
                ),
                "chembl_ids": ",".join(
                    str(value) for value in comp.get("chembl_ids") or []
                ),
                "chembl_id": comp.get("chembl_id", ""),
                "MW": descriptors.get("MW", ""),
                "LogP": descriptors.get("LogP", ""),
                "HBD": descriptors.get("HBD", ""),
                "HBA": descriptors.get("HBA", ""),
                "RotBonds": descriptors.get("RotBonds", ""),
                "TPSA": descriptors.get("TPSA", ""),
                "coordinate_origin": comp.get("coordinate_origin", ""),
                "raw_pose_preserved": comp.get("raw_pose_preserved", ""),
                "source_pose_file": comp.get("source_file", ""),
                "source_pose_index": comp.get("source_mol_index", ""),
            }
            writer.writerow(row)


def _manifest_value(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def validate_prepared_library(sdf_path: Path, compounds: list[dict]) -> dict:
    """Check the handoff once: IDs/roles, source lineage, exact chemistry and 3D.

    This does not enumerate chemical states or require a unique pH species.
    Undefined input stereochemistry may acquire a 3D orientation, while any
    explicitly specified stereochemistry must be retained.
    """
    import math
    from rdkit import Chem

    expected = {str(c["compound_id"]): c for c in compounds}
    if not expected or len(expected) != len(compounds):
        raise ValueError("Prepared library has empty or duplicate compound IDs")
    observed = set()
    for mol in Chem.SDMolSupplier(str(sdf_path), removeHs=False):
        if mol is None or not mol.HasProp("compound_id"):
            raise ValueError("Prepared library contains an unreadable or unidentified record")
        cid = mol.GetProp("compound_id")
        if cid not in expected or cid in observed:
            raise ValueError(f"Prepared library has an unexpected or duplicate compound: {cid}")
        record = expected[cid]
        for field in ("candidate_role", "source_record_id", "source_tool", "cohort_id", "campaign_mode", "conformer_id"):
            if not record.get(field) or not mol.HasProp(field) or mol.GetProp(field) != str(record[field]):
                raise ValueError(f"Prepared library {cid} lost {field}")
        if record["candidate_role"] not in {"active", "known", "generated", "validation_decoy"}:
            raise ValueError(f"Prepared library {cid} has an unknown role")
        if cid.startswith("DEC_") != (record["candidate_role"] == "validation_decoy"):
            raise ValueError(f"Prepared library {cid} changed its benchmark-control role")
        if record["candidate_role"] == "validation_decoy" and record["campaign_mode"] != "retrospective":
            raise ValueError(f"Prepared library {cid} incorrectly marks a control as prospective")
        if record["campaign_mode"] == "retrospective":
            for field in ("benchmark_dataset", "benchmark_split"):
                if not record.get(field) or not mol.HasProp(field) or mol.GetProp(field) != str(record[field]):
                    raise ValueError(f"Prepared library {cid} lost {field}")
        graph = Chem.RemoveHs(mol)
        intended = Chem.MolFromSmiles(record.get("smiles_standardized") or record.get("smiles", ""))
        if intended is None:
            raise ValueError(f"Prepared library {cid} has no standardized molecular graph")
        graph_without_stereo, intended_without_stereo = Chem.Mol(graph), Chem.Mol(intended)
        Chem.RemoveStereochemistry(graph_without_stereo)
        Chem.RemoveStereochemistry(intended_without_stereo)
        if (Chem.MolToSmiles(graph_without_stereo) != Chem.MolToSmiles(intended_without_stereo)
                or not graph.HasSubstructMatch(intended, useChirality=True)):
            raise ValueError(f"Prepared library {cid} changed connectivity, charge or specified stereochemistry")
        if record.get("coordinate_smiles") != Chem.MolToSmiles(graph):
            raise ValueError(f"Prepared library {cid} lost its prepared coordinate identity")
        if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
            raise ValueError(f"Prepared library {cid} is not 3D-ready")
        coordinates = [tuple(point) for point in mol.GetConformer().GetPositions()]
        if not all(math.isfinite(v) for point in coordinates for v in point):
            raise ValueError(f"Prepared library {cid} contains non-finite coordinates")
        if mol.GetNumAtoms() > 1 and len(set(coordinates)) < 2:
            raise ValueError(f"Prepared library {cid} contains collapsed coordinates")
        observed.add(cid)
    if observed != set(expected):
        raise ValueError("Prepared library lost expected coordinate records")
    return {"status": "passed", "records": len(observed)}


def _get_descriptors(comp: dict) -> dict:
    """Return descriptor values from metadata, computing from SMILES if needed."""
    descriptor_keys = ["MW", "LogP", "HBD", "HBA", "RotBonds", "TPSA"]
    descriptors = {key: comp.get(key, "") for key in descriptor_keys}
    if all(descriptors[key] not in (None, "") for key in descriptor_keys):
        return descriptors

    smiles = comp.get("smiles_standardized") or comp.get("smiles")
    if not smiles:
        return descriptors

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return descriptors
        computed = {
            "MW": f"{Descriptors.ExactMolWt(mol):.2f}",
            "LogP": f"{Descriptors.MolLogP(mol):.2f}",
            "HBD": str(rdMolDescriptors.CalcNumHBD(mol)),
            "HBA": str(rdMolDescriptors.CalcNumHBA(mol)),
            "RotBonds": str(rdMolDescriptors.CalcNumRotatableBonds(mol)),
            "TPSA": f"{Descriptors.TPSA(mol):.2f}",
        }
        for key, value in computed.items():
            if descriptors.get(key) in (None, ""):
                descriptors[key] = value
                comp[key] = value
    except Exception:
        pass

    return descriptors


def _build_report(
    target_info: dict,
    actives: list[dict],
    generated: list[dict],
    decoys: list[dict],
    steps_completed: list[str],
    steps_skipped: list[dict],
    parameters: dict,
) -> dict:
    """Build the library_report.json content."""
    exact_activity_values = [
        float(active["activity_nm"])
        for active in actives
        if active.get("activity_nm") is not None
        and str(active.get("standard_relation") or "=") == "="
        and not active.get("activity_is_censored", False)
    ]
    censored_activity_summaries = sum(
        1
        for active in actives
        if active.get("activity_nm") is not None
        and (
            str(active.get("standard_relation") or "=") != "="
            or active.get("activity_is_censored", False)
        )
    )
    activity_measurement_count = sum(
        int(active.get("activity_measurement_count") or 0) for active in actives
    )
    activity_aggregation_measurement_count = sum(
        int(active.get("activity_aggregation_measurement_count") or 0)
        for active in actives
    )
    activity_quarantined_measurement_count = sum(
        int(active.get("activity_quarantined_measurement_count") or 0)
        for active in actives
    )
    activity_data_validity_quarantined_count = sum(
        int(active.get("activity_data_validity_quarantined_count") or 0)
        for active in actives
    )
    activity_potential_duplicate_quarantined_count = sum(
        int(active.get("activity_potential_duplicate_quarantined_count") or 0)
        for active in actives
    )
    assay_heterogeneous_actives = [
        {
            "compound_id": active.get("compound_id"),
            "activity_assay_types": active.get("activity_assay_types") or [],
            "activity_bao_formats": active.get("activity_bao_formats") or [],
            "activity_assay_strata": active.get("activity_assay_strata") or [],
        }
        for active in actives
        if active.get("activity_assay_heterogeneity")
    ]

    decoy_sources = set(d.get("source", "unknown") for d in decoys)
    generated_by_source = {}
    generated_contribution_by_source = {}
    for compound in generated:
        source = str(compound.get("source") or "unknown")
        generated_by_source[source] = generated_by_source.get(source, 0) + 1
        sources = compound.get("source_tools") or [source]
        for contributing_source in set(str(value) for value in sources if value):
            generated_contribution_by_source[contributing_source] = (
                generated_contribution_by_source.get(contributing_source, 0) + 1
            )

    report = {
        "success": True,
        "target": target_info,
        "summary": {
            "num_actives": len(actives),
            "num_generated": len(generated),
            "generated_by_source": dict(sorted(generated_by_source.items())),
            "generated_contribution_by_source": dict(
                sorted(generated_contribution_by_source.items())
            ),
            "num_decoys": len(decoys),
            "num_candidates": len(actives) + len(generated),
            "total_compounds": len(actives) + len(generated) + len(decoys),
            "activity_type": parameters.get("activity_type", "IC50"),
            "activity_threshold_nm": parameters.get("activity_threshold_nm", 1000),
        },
        "library_stats": {
            # Backward-compatible names now summarize exact point estimates
            # only. Censored bounds are counted explicitly instead of being
            # silently inserted into an ordinary median.
            "actives_activity_range_nM": (
                [min(exact_activity_values), max(exact_activity_values)]
                if exact_activity_values else []
            ),
            "actives_median_nM": (
                _median(exact_activity_values) if exact_activity_values else None
            ),
            "actives_activity_summary_basis": "exact_point_estimates_only",
            "actives_censored_summary_count": censored_activity_summaries,
            "actives_activity_measurement_count": activity_measurement_count,
            "actives_activity_aggregation_measurement_count": (
                activity_aggregation_measurement_count
            ),
            "actives_activity_quarantined_measurement_count": (
                activity_quarantined_measurement_count
            ),
            "actives_activity_data_validity_quarantined_count": (
                activity_data_validity_quarantined_count
            ),
            "actives_activity_potential_duplicate_quarantined_count": (
                activity_potential_duplicate_quarantined_count
            ),
            "actives_assay_heterogeneity_count": len(
                assay_heterogeneous_actives
            ),
            "actives_assay_heterogeneity": assay_heterogeneous_actives,
            "decoys_source": list(decoy_sources),
        },
        "steps_completed": steps_completed,
        "steps_skipped": steps_skipped,
        "parameters": parameters,
    }

    return report


def _median(values: list[float]) -> float:
    """Return the conventional median for a non-empty numeric list."""
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[midpoint]
    return (ordered[midpoint - 1] + ordered[midpoint]) / 2

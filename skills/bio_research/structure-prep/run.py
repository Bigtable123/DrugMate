#!/usr/bin/env python3
"""Structure Preparation Skill - Standalone CLI.

Downloads PDB structure from RCSB, extracts receptor + ligand, then prepares
the protein for docking (PDBFixer, PROPKA, RDKit ligand conversion).

Supports three modes:
  1. Full pipeline: pdb_id/protein_name → download + prep
  2. Intake only:   pdb_id/protein_name + skip_prep=true → download only
  3. Prep only:     receptor_file → prep existing file (no download)

Usage:
    python run.py --input '{"pdb_id": "4KD1"}' --output-dir ./output/
    python run.py --input '{"pdb_id": "4KD1", "skip_prep": true}' --output-dir ./output/
    python run.py --input '{"receptor_file": "receptor.pdb"}' --output-dir ./output/
    python run.py --input-file params.json --output-dir ./output/
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from pdb_utils import (
    download_pdb_file,
    fetch_all_nonpolymer_entities,
    fetch_pdb_metadata,
    fetch_polymer_entities,
    parse_pdb_extract_components,
    search_pdb_by_protein_name,
)
from rcsb_http import RCSBRequestError


def _download_source_mmcif_best_effort(
    pdb_id: str,
    output_dir: Path,
) -> tuple[Path | None, dict]:
    """Fetch the untouched source complex for mmCIF-only downstream tools.

    The PDB file remains the required Step 0 intake.  Arpeggio's source mmCIF
    is supplemental evidence, so a transient RCSB/mmCIF failure must be
    explicit and auditable without turning an otherwise usable preparation
    into a failed run.
    """
    try:
        path = download_pdb_file(pdb_id, output_dir, fmt="cif")
    except Exception as exc:
        return None, {
            "status": "unavailable",
            "reason": f"best-effort mmCIF download failed: {exc}",
        }
    return path, {
        "status": "available",
        "path": str(path),
        "format": "mmCIF",
        "provenance": "RCSB source complex",
    }


def _existing_optional_path(value: object) -> Path | None:
    """Resolve an optional local evidence file without publishing stale paths."""
    if not isinstance(value, (str, Path)) or not str(value).strip():
        return None
    path = Path(value).expanduser()
    try:
        return path.resolve() if path.is_file() else None
    except OSError:
        return None


def _rcsb_intake_failure(
    *,
    stage: str,
    pdb_id: str | None,
    error: RCSBRequestError,
) -> dict:
    """Build the structured failure contract for exhausted RCSB requests."""
    details = {
        "service": "RCSB",
        "stage": stage,
        "pdb_id": pdb_id,
        **error.to_dict(),
    }
    return {
        "success": False,
        "error_code": "rcsb_request_failed",
        "error": f"RCSB {stage} failed for {pdb_id or 'protein search'}: {error}",
        "error_details": details,
    }


def _run_intake(input_data: dict, output_dir: Path) -> dict:
    """Phase 1: Download and parse PDB structure from RCSB."""
    from io_utils import tool_enabled

    pdb_id = input_data.get("pdb_id", "").upper()
    protein_name = input_data.get("protein_name", "")
    water_threshold = input_data.get("water_distance_threshold", 3.5)
    keep_auxiliary_hetero = bool(input_data.get("keep_auxiliary_hetero", False))
    preserve_hetero_resnames = input_data.get("preserve_hetero_resnames")
    water_protein_distance = input_data.get("water_protein_distance", 3.5)
    water_min_occupancy = input_data.get("water_min_occupancy", 0.5)
    water_max_b_factor = input_data.get("water_max_b_factor", 60.0)
    water_include_residues = input_data.get("water_include_residues", [])
    water_exclude_residues = input_data.get("water_exclude_residues", [])
    keep_selected_waters = False  # Inventory now; select waters only in docking_prep.

    if not pdb_id and not protein_name:
        return {"success": False, "error": "Must provide pdb_id or protein_name"}

    if not pdb_id and protein_name:
        try:
            results = search_pdb_by_protein_name(protein_name)
        except RCSBRequestError as exc:
            return _rcsb_intake_failure(
                stage="protein_search",
                pdb_id=None,
                error=exc,
            )
        if not results:
            return {"success": False, "error": f"No PDB entries found for '{protein_name}'"}
        pdb_id = results[0]["pdb_id"]

    if len(pdb_id) != 4 or not pdb_id.isalnum():
        return {"success": False, "error": f"Invalid PDB ID: {pdb_id}"}

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        metadata_raw = fetch_pdb_metadata(pdb_id)
    except RCSBRequestError as exc:
        return _rcsb_intake_failure(
            stage="entry_metadata",
            pdb_id=pdb_id,
            error=exc,
        )
    except Exception as e:
        return {"success": False, "error": f"Failed to fetch metadata: {e}"}

    entry_info = metadata_raw.get("rcsb_entry_info", {})
    struct = metadata_raw.get("struct", {})
    metadata = {
        "pdb_id": pdb_id,
        "title": struct.get("title", ""),
        "resolution": entry_info.get("resolution_combined", [None])[0],
        "method": entry_info.get("experimental_method", ""),
        "deposition_date": metadata_raw.get("rcsb_accession_info", {}).get("deposit_date", ""),
    }

    try:
        polymer_entities = fetch_polymer_entities(pdb_id)
    except RCSBRequestError as exc:
        return _rcsb_intake_failure(
            stage="polymer_entities",
            pdb_id=pdb_id,
            error=exc,
        )
    chains_metadata = []
    for ent in polymer_entities:
        entity_info = ent.get("rcsb_polymer_entity", {})
        container = ent.get("rcsb_polymer_entity_container_identifiers", {})
        gene_names = [g.get("value", "") for g in ent.get("rcsb_gene_name", [])] if ent.get("rcsb_gene_name") else []
        chains_metadata.append({
            "entity_id": container.get("entity_id", ""),
            "chain_ids": container.get("asym_ids", []),
            "author_chain_ids": container.get("auth_asym_ids", []),
            "description": entity_info.get("pdbx_description", ""),
            "organism": ent.get("rcsb_entity_source_organism", [{}])[0].get("ncbi_scientific_name", "") if ent.get("rcsb_entity_source_organism") else "",
            "gene_names": gene_names,
        })
    metadata["chains"] = chains_metadata
    metadata["selected_chain"] = input_data.get("chain")

    try:
        ligand_entities = fetch_all_nonpolymer_entities(pdb_id)
    except RCSBRequestError as exc:
        return _rcsb_intake_failure(
            stage="nonpolymer_entities",
            pdb_id=pdb_id,
            error=exc,
        )
    ligands_info = []
    for ent in ligand_entities:
        nonpoly = ent.get("rcsb_nonpolymer_entity", {})
        comp_id = ent.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("comp_id", "")
        if comp_id:
            ligands_info.append({"comp_id": comp_id, "name": nonpoly.get("pdbx_description", "")})
    metadata["ligands"] = ligands_info

    try:
        pdb_file = download_pdb_file(pdb_id, output_dir, fmt="pdb")
    except RCSBRequestError as exc:
        return _rcsb_intake_failure(
            stage="pdb_download",
            pdb_id=pdb_id,
            error=exc,
        )
    except Exception as e:
        return {"success": False, "error": f"Failed to download PDB: {e}"}
    complex_mmcif, complex_mmcif_status = _download_source_mmcif_best_effort(
        pdb_id, output_dir,
    )

    known_ligand_ids = [lig["comp_id"] for lig in ligands_info]
    extraction = parse_pdb_extract_components(
        pdb_file=pdb_file, output_dir=output_dir,
        ligand_ids=known_ligand_ids if known_ligand_ids else None,
        keep_conserved_waters=keep_selected_waters, water_distance_threshold=water_threshold,
        keep_auxiliary_hetero=keep_auxiliary_hetero,
        preserve_hetero_resnames=preserve_hetero_resnames,
        water_protein_distance=water_protein_distance,
        water_min_occupancy=water_min_occupancy,
        water_max_b_factor=water_max_b_factor,
        water_include_residues=water_include_residues,
        water_exclude_residues=water_exclude_residues,
    )

    if not ligands_info and extraction["selected_ligand"]:
        ligands_info = [{"comp_id": extraction["selected_ligand"], "name": ""}]
        metadata["ligands"] = ligands_info

    metadata["reference_ligand"] = extraction["selected_ligand"]
    metadata["conserved_waters_count"] = extraction["conserved_waters_count"]
    metadata["removed_auxiliary_hetero"] = extraction["removed_auxiliary_hetero"]
    metadata["preserved_auxiliary_hetero"] = extraction["preserved_auxiliary_hetero"]
    metadata["preserved_modified_polymer_residues"] = extraction[
        "preserved_modified_polymer_residues"
    ]
    metadata["chain_info"] = extraction["chain_info"]
    metadata["source_complex_mmcif"] = complex_mmcif_status

    water_report_file = output_dir / "water_report.json"
    water_report_file.write_text(json.dumps(
        extraction["water_selection"], indent=2, ensure_ascii=False,
    ))

    metadata_file = output_dir / "metadata.json"
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    return {
        "success": True,
        "pdb_id": pdb_id,
        "metadata": metadata,
        "files": {
            "pdb_file": str(pdb_file),
            "complex_mmcif": str(complex_mmcif) if complex_mmcif else None,
            "receptor_file": str(extraction["receptor_file"]),
            "ligand_file": str(extraction["ligand_file"]) if extraction["ligand_file"] else None,
            "metadata_file": str(metadata_file),
            "water_report": str(water_report_file),
        },
        "extraction": {
            "keep_auxiliary_hetero": extraction["keep_auxiliary_hetero"],
            "preserve_hetero_resnames": extraction["preserve_hetero_resnames"],
            "removed_auxiliary_hetero": extraction["removed_auxiliary_hetero"],
            "preserved_auxiliary_hetero": extraction["preserved_auxiliary_hetero"],
            "preserved_modified_polymer_residues": extraction[
                "preserved_modified_polymer_residues"
            ],
        },
    }


def _run_prep(receptor_file: Path, ligand_file: Path | None, input_data: dict, output_dir: Path) -> dict:
    """Phase 2: Prepare protein in a deterministic chemistry-aware order."""
    from prep_utils import (
        audit_structure,
        check_dependencies,
        fix_protein_with_pdbfixer,
        protonate_with_propka,
        validate_metals,
        validate_residues,
    )
    from rdkit_utils import check_rdkit_available, pdb_to_sdf, generate_2d_sdf
    from external_tools import run_reduce2
    from io_utils import tool_config, tool_enabled, tool_switches

    ph = input_data.get("ph")
    add_missing_residues = bool(input_data.get("add_missing_residues", False))
    add_missing_atoms = bool(input_data.get("add_missing_atoms", True))
    add_hydrogens = bool(input_data.get("add_hydrogens", False))
    required_residues = input_data.get("required_residues", None)
    ligand_comp_id = input_data.get("ligand_comp_id")
    if not ligand_comp_id and ligand_file:
        ligand_comp_id = _infer_ligand_comp_id(ligand_file)

    # Standalone geometry defaults match the registry. Full experiments still
    # require the complete explicit selection matrix.
    switches = tool_switches(
        input_data, ("pdbfixer", "ligand_sdf", "residue_validation"),
        default=True,
    )
    switches["propka"] = tool_enabled(input_data, "propka", default=False)
    switches["reduce"] = tool_enabled(input_data, "reduce", default=False)
    switches["metal_validation"] = tool_enabled(input_data, "metal_validation", default=True)
    switches["structure_qc"] = tool_enabled(input_data, "structure_qc", default=True)

    def config(tool_id: str) -> dict:
        return tool_config(input_data, tool_id)

    def params(tool_id: str) -> dict:
        return dict(config(tool_id).get("params") or {})

    def runtime(tool_id: str) -> dict:
        return dict(config(tool_id).get("runtime") or {})

    if not receptor_file.exists():
        return {"success": False, "error": f"Receptor file not found: {receptor_file}"}

    output_dir.mkdir(parents=True, exist_ok=True)

    deps = check_dependencies()
    rdkit_available = check_rdkit_available()

    report = {
        "ph": ph,
        "dependencies": {**deps, "rdkit": rdkit_available},
        "steps_completed": [],
        "steps_skipped": [],
        "parameters": {
            "ligand_comp_id": ligand_comp_id,
            "add_missing_residues": add_missing_residues,
            "add_missing_atoms": add_missing_atoms,
            "add_hydrogens": add_hydrogens,
        },
    }

    # Preserve a geometry-only source even if deposited coordinates contain H.
    from receptor_preparation import write_noh, file_identity, write_json
    prep_input_pdb = output_dir / "input_noH.pdb"
    write_noh(receptor_file, prep_input_pdb, remove_waters=True)
    protonation_manifest = None
    report["input_qc"] = audit_structure(prep_input_pdb, require_hydrogens=False)
    if not report["input_qc"]["passed"]:
        return _prep_failure(report, output_dir, "Input basic structure QC failed")

    # Step 1: PDBFixer
    prepared_pdb = output_dir / f"{receptor_file.stem}_prepared.pdb"

    if not switches["pdbfixer"]:
        report["steps_skipped"].append({"step": "pdbfixer", "reason": "disabled_by_config"})
        prepared_pdb = prep_input_pdb
    elif deps.get("pdbfixer") and deps.get("openmm"):
        try:
            fix_report = fix_protein_with_pdbfixer(
                input_pdb=prep_input_pdb, output_pdb=prepared_pdb,
                ph=ph, add_hydrogens=add_hydrogens,
                add_missing_residues=add_missing_residues,
                add_missing_atoms=add_missing_atoms,
                max_missing_residues=input_data.get("max_missing_residues", 3),
                seed=input_data.get("seed", 0),
            )
            report["pdbfixer"] = fix_report
            report["steps_completed"].append("pdbfixer")
            report["parameters"]["pdbfixer_ph"] = ph
            try:
                import openmm
                report["parameters"]["openmm_version"] = openmm.__version__
            except Exception:
                pass
        except Exception as e:
            return _prep_failure(
                report, output_dir,
                f"Enabled tool pdbfixer failed: {e}",
            )
    else:
        missing = [name for name in ("pdbfixer", "openmm") if not deps.get(name)]
        return _prep_failure(report, output_dir, "Enabled pdbfixer is unavailable: " + ", ".join(missing))

    report["post_repair_qc"] = audit_structure(prepared_pdb, require_hydrogens=False)
    if not report["post_repair_qc"]["passed"]:
        return _prep_failure(report, output_dir, "Post-repair basic structure QC failed")

    # Step 1b: Generate noH version (strip hydrogens for tools that need it)
    noH_pdb = output_dir / f"{receptor_file.stem}_noH.pdb"
    if prepared_pdb.exists():
        try:
            _write_noh_pdb(prepared_pdb, noH_pdb)
            report["steps_completed"].append("noH_generation")
        except Exception as e:
            report["steps_skipped"].append({"step": "noH_generation", "reason": str(e)})
            noH_pdb = None
    else:
        noH_pdb = None

    # Step 2: PROPKA
    if not switches["propka"]:
        report["steps_skipped"].append({"step": "propka", "reason": "disabled_by_config"})
    elif deps.get("propka"):
        try:
            # Run PROPKA on the PDBFixer-prepared structure (complete side chains
            # and neighbours) rather than the raw input. pKa predictions depend on
            # a complete environment; a structure with missing atoms/residues
            # yields incorrect protonation states. Falls back to the raw receptor
            # only when preparation was skipped/failed (prepared_pdb == receptor_file).
            propka_input = prepared_pdb if prepared_pdb.exists() else receptor_file
            propka_report = protonate_with_propka(propka_input, ph=ph)
            propka_report["input_structure"] = str(propka_input)
            propka_report["ran_on_prepared_structure"] = propka_input != receptor_file
            report["propka"] = propka_report
            report["steps_completed"].append("propka")
            report["parameters"]["propka_ph"] = ph
            try:
                import propka
                report["parameters"]["propka_version"] = getattr(propka, "__version__", "unknown")
            except Exception:
                pass
        except Exception as e:
            return _prep_failure(
                report, output_dir,
                f"Enabled tool propka failed: {e}",
            )
    else:
        # PROPKA is consumed only as ``propka_diagnostic`` in the protonation
        # manifest below; the backends that actually set protonation states are
        # reduce and pdb2pqr.  A missing diagnostic package therefore has no
        # bearing on the prepared structure, and aborting on it ended the round
        # for evidence nothing downstream gates on.
        report["steps_skipped"].append({
            "step": "propka",
            "reason": "Enabled tool propka is unavailable (missing: propka)",
        })

    # Only a separately recorded provisional Reduce derivative is permitted here.
    provisional_pdb = None
    if switches["reduce"]:
        provisional_pdb = output_dir / "provisional_receptor.pdb"
        reduce_params = params("reduce")
        reason = reduce_params.pop("reason", None)
        if not reason:
            return _prep_failure(report, output_dir, "Provisional reduce requires a hydrogen-dependent pocket-analysis reason")
        try:
            backend_result = run_reduce2(noH_pdb, provisional_pdb, runtime("reduce"), **reduce_params)
            protonation_manifest = output_dir / "provisional_protonation_manifest.json"
            write_json(protonation_manifest, {
                "backend": "reduce", "state": "provisional", "reason": reason,
                "source": file_identity(noH_pdb), "output": file_identity(provisional_pdb),
                "runtime_result": backend_result,
            })
            report["steps_completed"].append("reduce")
        except Exception as exc:
            return _prep_failure(report, output_dir, f"Provisional reduce failed: {exc}")

    # Step 3: Ligand → SDF
    ligand_sdf = None
    ligand_2d_sdf = None
    ligand_properties = None

    if not switches["ligand_sdf"]:
        report["steps_skipped"].append({"step": "ligand_sdf_conversion", "reason": "disabled_by_config"})
    elif ligand_file and ligand_file.exists() and rdkit_available:
        sdf_path = output_dir / f"{ligand_file.stem}_pose.sdf"
        try:
            sdf_result = pdb_to_sdf(ligand_file, sdf_path, comp_id=ligand_comp_id)
            if sdf_result["success"]:
                ligand_sdf = sdf_path
                ligand_properties = sdf_result["properties"]
                report["ligand_sdf_conversion"] = {
                    "bond_order_source": ligand_properties.get("bond_order_source"),
                    "template_comp_id": ligand_properties.get("template_comp_id"),
                    "template_file": ligand_properties.get("template_file"),
                    "template_graph_verified": ligand_properties.get("template_graph_verified", False),
                    "crystal_coordinates_preserved": ligand_properties.get("crystal_coordinates_preserved", False),
                }
                report["steps_completed"].append("ligand_sdf_conversion")
                if not ligand_properties.get("is_3d"):
                    report.setdefault("warnings", []).append(
                        "Ligand SDF has flat/2D coordinates — interaction analysis will fail. "
                        "Check that ligand PDB has original crystal coordinates."
                    )
            else:
                return _prep_failure(
                    report, output_dir,
                    "Enabled tool ligand_sdf failed: " + sdf_result.get("error", "unknown error"),
                )
        except Exception as e:
            return _prep_failure(
                report, output_dir,
                f"Enabled tool ligand_sdf failed: {e}",
            )

        sdf_2d_path = output_dir / f"{ligand_file.stem}_2d.sdf"
        try:
            result_2d = generate_2d_sdf(ligand_file, sdf_2d_path, comp_id=ligand_comp_id)
            if result_2d["success"]:
                ligand_2d_sdf = sdf_2d_path
        except Exception:
            pass
    elif ligand_file and not rdkit_available:
        return _prep_failure(report, output_dir, "Enabled tool ligand_sdf is unavailable (missing: rdkit)")
    elif not ligand_file:
        return _prep_failure(
            report, output_dir,
            "Enabled tool ligand_sdf requires a ligand file; disable it when no reference ligand is supplied",
        )
    else:
        return _prep_failure(report, output_dir, f"Ligand file not found: {ligand_file}")

    ligand_microstates = None
    openbabel_output = None

    metal_report_file = None
    if switches["metal_validation"]:
        try:
            metal_report = validate_metals(prepared_pdb, **params("metal_validation"))
            metal_report_file = output_dir / "metal_report.json"
            metal_report_file.write_text(json.dumps(
                metal_report, indent=2, ensure_ascii=False, default=str,
            ))
            report["metal_validation"] = metal_report
            report["steps_completed"].append("metal_validation")
        except Exception as exc:
            return _prep_failure(report, output_dir, f"Enabled tool metal_validation failed: {exc}")
    else:
        report["steps_skipped"].append({"step": "metal_validation", "reason": "disabled_by_config"})

    receptor_pdbqt = None
    amber_outputs = None

    # Step 4: Residue validation
    if not switches["residue_validation"]:
        report["steps_skipped"].append({"step": "residue_validation", "reason": "disabled_by_config"})
        validation = {"validation_passed": None}
    else:
        try:
            validation = validate_residues(prepared_pdb, required_residues=required_residues)
            report["validation"] = validation
            report["steps_completed"].append("residue_validation")
        except Exception as e:
            return _prep_failure(
                report, output_dir,
                f"Enabled tool residue_validation failed: {e}",
            )

        if validation.get("validation_passed") is False:
            return _prep_failure(
                report, output_dir,
                "Active-site residue validation failed",
            )

    structure_qc_file = None
    if switches["structure_qc"]:
        expected_outputs = {
            "prepared_noH_pdb": str(noH_pdb) if noH_pdb else None,
            "ligand_sdf": str(ligand_sdf) if ligand_sdf else None,
            "receptor_pdbqt": str(receptor_pdbqt) if receptor_pdbqt else None,
            "amber_prmtop": (amber_outputs or {}).get("prmtop"),
            "amber_inpcrd": (amber_outputs or {}).get("inpcrd"),
        }
        try:
            # ``add_hydrogens`` is the run's own instruction, so pass it through:
            # a preparation told not to add hydrogens must not then fail QC for
            # not having any.  An explicit ``structure_qc`` param
            # still wins, so a caller can demand hydrogens regardless -- but the
            # registry spells "unset" as ``None`` (shared/pipeline_registry.py
            # ``_p``), and ``None`` is falsy, so test for it rather than using
            # ``setdefault``: a materialized default must not silently switch the
            # hydrogen check off.
            qc_params = dict(params("structure_qc"))
            if qc_params.get("require_hydrogens") is None:
                qc_params["require_hydrogens"] = bool(add_hydrogens)
            structure_qc = audit_structure(
                prepared_pdb, expected_outputs=expected_outputs, **qc_params,
            )
            structure_qc_file = output_dir / "structure_qc_report.json"
            structure_qc_file.write_text(json.dumps(
                structure_qc, indent=2, ensure_ascii=False, default=str,
            ))
            report["structure_qc"] = structure_qc
            report["steps_completed"].append("structure_qc")
            if not structure_qc["passed"]:
                return _prep_failure(
                    report, output_dir,
                    "Final structure QC failed: " + "; ".join(structure_qc["issues"]),
                )
        except Exception as exc:
            return _prep_failure(report, output_dir, f"Enabled tool structure_qc failed: {exc}")
    else:
        report["steps_skipped"].append({"step": "structure_qc", "reason": "disabled_by_config"})

    prep_manifest = output_dir / "prep_manifest.json"
    write_json(prep_manifest, {
        "schema_version": 1, "stage": "geometry_only", "source": file_identity(receptor_file),
        "prepared_receptor_noH": file_identity(noH_pdb), "chain": input_data.get("chain"),
        "parameters": report["parameters"], "repair_reason": input_data.get("repair_reason"),
        "input_qc": report["input_qc"], "post_repair_qc": report["post_repair_qc"],
        "qc_scope": "basic structure/file checks; geometric/valence QC pending",
        "final_protonation": "deferred_to_docking_prep", "water_selection": "deferred_to_docking_prep",
    })
    report_file = output_dir / "preparation_report.json"
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    return {
        "success": True,
        "prepared_pdb": str(noH_pdb),
        "provisional_pdb": str(provisional_pdb) if provisional_pdb else None,
        "prep_manifest": str(prep_manifest),
        "prepared_noH_pdb": str(noH_pdb) if noH_pdb else None,
        "ligand_sdf": str(ligand_sdf) if ligand_sdf else None,
        "ligand_2d_sdf": str(ligand_2d_sdf) if ligand_2d_sdf else None,
        "ligand_properties": ligand_properties,
        "ligand_microstates": str(ligand_microstates) if ligand_microstates else None,
        "openbabel_output": str(openbabel_output) if openbabel_output else None,
        "receptor_pdbqt": str(receptor_pdbqt) if receptor_pdbqt else None,
        "protonation_manifest": str(protonation_manifest) if protonation_manifest else None,
        "metal_report": str(metal_report_file) if metal_report_file else None,
        "structure_qc_report": str(structure_qc_file) if structure_qc_file else None,
        "amber_outputs": amber_outputs,
        "preparation_report": str(report_file),
        "validation": validation,
        "steps_completed": report["steps_completed"],
        "steps_skipped": report["steps_skipped"],
        "parameters": report["parameters"],
    }


def _prep_failure(report: dict, output_dir: Path, message: str) -> dict:
    """Persist a failed preparation report and return the skill failure payload."""
    report["success"] = False
    report["error"] = message
    report_file = output_dir / "preparation_report.json"
    with open(report_file, "w") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False, default=str)
    return {
        "success": False,
        "error": message,
        "preparation_report": str(report_file),
        "steps_completed": report.get("steps_completed", []),
        "steps_skipped": report.get("steps_skipped", []),
        "parameters": report.get("parameters", {}),
    }


def run(input_data: dict, output_dir: Path) -> dict:
    """Main entry: determine mode and execute accordingly."""
    from io_utils import tool_enabled, tool_switches

    from receptor_preparation import DEFERRED_TOOLS, water_inventory
    for name in sorted(DEFERRED_TOOLS | {"amber_topology"}):
        if tool_enabled(input_data, name, default=False):
            destination = "pose_evaluation candidate MM preparation" if name == "amber_topology" else "docking"
            return {"success": False, "error": f"structure_prep.{name} moved to {destination}; regenerate the explicit tool matrix"}
    if input_data.get("add_hydrogens") is True:
        return {"success": False, "error": "structure_prep.add_hydrogens must be false; final hydrogens belong to docking_prep"}
    if input_data.get("add_missing_residues") and not input_data.get("repair_reason"):
        return {"success": False, "error": "add_missing_residues requires repair_reason for supported short internal gaps"}
    if tool_enabled(input_data, "propka", default=False) and input_data.get("ph") is None:
        return {"success": False, "error": "Diagnostic PROPKA requires an explicit justified pH"}

    has_pdb_id = bool(input_data.get("pdb_id") or input_data.get("protein_name"))
    has_receptor_file = bool(input_data.get("receptor_file"))
    skip_prep = input_data.get("skip_prep", False)
    intake_enabled = tool_enabled(input_data, "intake", default=has_pdb_id)

    if intake_enabled and not has_pdb_id:
        return {
            "success": False,
            "error": "intake is enabled but pdb_id/protein_name was not provided",
        }
    if not intake_enabled and not has_receptor_file:
        return {
            "success": False,
            "error": "intake is disabled but receptor_file was not provided",
        }

    prep_switches = tool_switches(
        input_data, (
            "pdbfixer", "propka", "reduce", "pdb2pqr", "protoss",
            "ligand_sdf", "ligand_microstates", "openbabel",
            "metal_validation", "meeko", "amber_topology", "structure_qc",
            "residue_validation",
        ),
        default=not skip_prep,
    )
    if skip_prep and any(prep_switches.values()):
        enabled = [name for name, value in prep_switches.items() if value]
        return {
            "success": False,
            "error": "skip_prep=true conflicts with enabled preparation tools: " + ", ".join(enabled),
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    result = {}

    # --- Phase 1: Intake (download + parse) ---
    if intake_enabled:
        intake_result = _run_intake(input_data, output_dir)
        if not intake_result["success"]:
            return intake_result

        result.update(intake_result)

        if skip_prep:
            return result

        receptor_file = Path(intake_result["files"]["receptor_file"])
        ligand_file_str = intake_result["files"].get("ligand_file")
        ligand_file = Path(ligand_file_str) if ligand_file_str else None
    else:
        receptor_file = Path(input_data["receptor_file"])
        ligand_file_str = input_data.get("ligand_file")
        ligand_file = Path(ligand_file_str) if ligand_file_str else None
        complex_mmcif = _existing_optional_path(input_data.get("complex_mmcif_file"))

    # --- Phase 2: Preparation ---
    source_lines = receptor_file.read_text().splitlines()
    if sum(line.startswith("MODEL ") for line in source_lines) > 1:
        return {"success": False, "error": "Multiple coordinate models require an explicitly selected local model before structure_prep"}
    selected_chain = input_data.get("chain")
    if selected_chain and not any(line.startswith(("ATOM  ", "HETATM")) and line[21].strip() == selected_chain
                                  for line in source_lines):
        return {"success": False, "error": f"Selected target chain {selected_chain!r} is absent from the receptor"}
    prep_result = _run_prep(receptor_file, ligand_file, input_data, output_dir)
    if not prep_result["success"]:
        return prep_result

    if "files" in result:
        result["files"]["prepared_pdb"] = prep_result["prepared_pdb"]
        result["files"]["prepared_noH_pdb"] = prep_result["prepared_noH_pdb"]
        result["files"]["ligand_sdf"] = prep_result["ligand_sdf"]
        result["files"]["ligand_2d_sdf"] = prep_result["ligand_2d_sdf"]
        result["files"]["preparation_report"] = prep_result["preparation_report"]
        result["files"]["ligand_microstates"] = prep_result["ligand_microstates"]
        result["files"]["openbabel_output"] = prep_result["openbabel_output"]
        result["files"]["receptor_pdbqt"] = prep_result["receptor_pdbqt"]
        result["files"]["protonation_manifest"] = prep_result["protonation_manifest"]
        result["files"]["metal_report"] = prep_result["metal_report"]
        result["files"]["structure_qc_report"] = prep_result["structure_qc_report"]
        result["files"]["amber_prmtop"] = (prep_result["amber_outputs"] or {}).get("prmtop")
        result["files"]["amber_inpcrd"] = (prep_result["amber_outputs"] or {}).get("inpcrd")
        result["files"]["amber_ligand_mol2"] = (prep_result["amber_outputs"] or {}).get("ligand_mol2")
        result["files"]["amber_ligand_frcmod"] = (prep_result["amber_outputs"] or {}).get("ligand_frcmod")
    else:
        result["success"] = True
        result["files"] = {
            "receptor_file": str(receptor_file),
            "ligand_file": str(ligand_file) if ligand_file else None,
            "complex_mmcif": str(complex_mmcif) if complex_mmcif else None,
            "prepared_pdb": prep_result["prepared_pdb"],
            "prepared_noH_pdb": prep_result["prepared_noH_pdb"],
            "ligand_sdf": prep_result["ligand_sdf"],
            "ligand_2d_sdf": prep_result["ligand_2d_sdf"],
            "preparation_report": prep_result["preparation_report"],
            "ligand_microstates": prep_result["ligand_microstates"],
            "openbabel_output": prep_result["openbabel_output"],
            "receptor_pdbqt": prep_result["receptor_pdbqt"],
            "protonation_manifest": prep_result["protonation_manifest"],
            "metal_report": prep_result["metal_report"],
            "structure_qc_report": prep_result["structure_qc_report"],
            "amber_prmtop": (prep_result["amber_outputs"] or {}).get("prmtop"),
            "amber_inpcrd": (prep_result["amber_outputs"] or {}).get("inpcrd"),
            "amber_ligand_mol2": (prep_result["amber_outputs"] or {}).get("ligand_mol2"),
            "amber_ligand_frcmod": (prep_result["amber_outputs"] or {}).get("ligand_frcmod"),
        }

    result["files"]["provisional_pdb"] = prep_result["provisional_pdb"]
    result["files"]["prep_manifest"] = prep_result["prep_manifest"]
    source_pdb = Path(result["files"].get("pdb_file") or receptor_file)
    result["files"]["source_complex_pdb"] = str(source_pdb)
    result["files"]["water_inventory"] = water_inventory(source_pdb, output_dir / "water_inventory.json")
    result["files"].pop("water_report", None)
    result["ligand_properties"] = prep_result["ligand_properties"]
    result["validation"] = prep_result["validation"]
    result["steps_completed"] = prep_result["steps_completed"]
    result["steps_skipped"] = prep_result["steps_skipped"]
    result["parameters"] = prep_result["parameters"]

    return result


def _infer_ligand_comp_id(ligand_file: Path) -> str | None:
    """Infer ligand component ID from files named like 4KD1_ligand_1QK.pdb."""
    marker = "_ligand_"
    stem = ligand_file.stem
    if marker not in stem:
        return None
    comp_id = stem.rsplit(marker, 1)[-1].strip()
    return comp_id.upper() if comp_id else None


def _pdb_line_is_hydrogen(line: str) -> bool:
    """Return True if a PDB ATOM/HETATM line describes a hydrogen atom.

    Prefers the element column (77-78). When absent, infers from the atom name
    (columns 13-16) using PDB naming rules: hydrogen names are frequently
    prefixed with a digit (e.g. "1HB", "2HG1", "HH11"), so a plain
    startswith("H") test misses them and leaves hydrogens in the "noH" file.
    """
    element = line[76:78].strip().upper() if len(line) > 76 else ""
    if element:
        return element == "H" or element == "D"  # include deuterium
    atom_name = line[12:16].strip()
    if not atom_name:
        return False
    # Strip a leading digit or altloc digit, then check the first alpha char.
    stripped = atom_name.lstrip("0123456789")
    return stripped[:1].upper() in ("H", "D")


def _write_noh_pdb(input_pdb: Path, output_pdb: Path) -> None:
    """Write a hydrogen-stripped copy while preserving all other PDB records."""
    with input_pdb.open() as source, output_pdb.open("w") as destination:
        for line in source:
            if line[:6] in ("ATOM  ", "HETATM") and _pdb_line_is_hydrogen(line):
                continue
            destination.write(line)


def main():
    from io_utils import parse_skill_args, print_skill_result
    try:
        input_data, output_dir = parse_skill_args("Structure Preparation Skill")
        result = run(input_data, output_dir)
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        result = {
            "success": False,
            "error_code": "structure_prep_unhandled_error",
            "error": f"Unhandled structure-prep error: {type(exc).__name__}: {exc}",
            "error_details": {
                "stage": "structure_prep",
                "error_type": type(exc).__name__,
                "message": str(exc),
            },
        }
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

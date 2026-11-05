"""3D coordinate preparation for generated and library compounds.

Generator-native 3D poses are preserved without force-field optimization.
Compounds without a compatible input pose use RDKit ETKDGv3 + MMFF/UFF.
"""
from __future__ import annotations

import json
import hashlib
import math
import sys
from pathlib import Path

try:
    from pose_identity import LIGAND_LINEAGE_PROPERTIES, matches_declared_molecule
except ModuleNotFoundError:
    # Utility modules are also imported directly by contract tests and external
    # callers, outside molecule-generator/run.py's path bootstrap.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
    from pose_identity import LIGAND_LINEAGE_PROPERTIES, matches_declared_molecule

from molecule_failures import MoleculeFailureTracker, is_rdkit_parse_error


def check_conformer_gen_available() -> bool:
    """Check if RDKit conformer generation is available."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        return True
    except ImportError:
        return False


def generate_conformers_sdf(
    compounds: list[dict],
    output_sdf: Path,
    smiles_key: str = "smiles_standardized",
    optimize: bool = True,
    random_seed: int = 42,
    failure_tracker: MoleculeFailureTracker | None = None,
) -> dict:
    """Generate 3D conformers for all compounds and write to SDF.

    Each compound dict should contain at minimum a SMILES field.
    SDF tags include source, activity value/relation/bounds, aggregation audit,
    SMILES, ChEMBL provenance, and standard molecular descriptors.

    Returns dict with success status and counts.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    output_sdf = Path(output_sdf)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    writer = Chem.SDWriter(str(output_sdf))
    total_written = 0
    failed_compounds = []
    written_compound_ids = []
    optimization_warnings = []
    raw_poses_preserved = 0
    raw_poses_rejected = []
    embedded_conformers = 0
    id_counters = {"active": 0, "decoy": 0, "generated": 0}
    tracker = failure_tracker or MoleculeFailureTracker()
    for comp in compounds:
        comp["conformer_written"] = False

    try:
        for comp in compounds:
            smiles = comp.get(smiles_key) or comp.get("smiles", "")
            try:
                mol, raw_pose_preserved, coordinate_origin, raw_pose_error, warning = (
                    _prepare_coordinates(comp, smiles, optimize, random_seed)
                )
                _record_coordinate_identity(comp, mol, smiles)
            except _MoleculeCoordinateError as exc:
                failed_compounds.append({**comp, "conformer_error": str(exc)})
                tracker.record(comp, "conformer.coordinates", exc.code, str(exc))
                continue
            except Exception as exc:
                failed_compounds.append({**comp, "conformer_error": str(exc)})
                tracker.record(comp, "conformer.coordinates", "conformer_generation_failed", str(exc),
                               exception=None if is_rdkit_parse_error(exc) else exc)
                continue
            if not raw_pose_preserved and comp.get("raw_pose_molblock"):
                raw_poses_rejected.append({"compound_id": comp.get("compound_id", ""),
                                           "reason": raw_pose_error})
            if warning:
                optimization_warnings.append({"compound_id": comp.get("compound_id", ""), **warning})

            comp_id = _next_compound_id(comp, id_counters)
            comp["compound_id"] = comp_id
            comp["raw_pose_preserved"] = raw_pose_preserved
            comp["coordinate_origin"] = coordinate_origin
            comp.setdefault("conformer_id", f"{comp_id}:conf:0001")
            mol.SetProp("_Name", str(comp_id))
            mol.SetProp("compound_id", str(comp_id))
            mol.SetProp("source", comp.get("source", "unknown"))
            mol.SetProp("candidate_role", comp.get("candidate_role", "unknown"))
            mol.SetProp("smiles", smiles)
            mol.SetProp("coordinate_origin", str(coordinate_origin))
            mol.SetProp("raw_pose_preserved", str(raw_pose_preserved).lower())
            if comp.get("parent_input_records"):
                mol.SetProp("parent_input_records_json", json.dumps(
                    comp["parent_input_records"], sort_keys=True, separators=(",", ":"),
                ))
            if comp.get("source_file"):
                mol.SetProp("source_pose_file", str(comp["source_file"]))
            if comp.get("source_mol_index") is not None:
                mol.SetProp("source_pose_index", str(comp["source_mol_index"]))
            if comp.get("source_tools"):
                mol.SetProp("source_tools", ",".join(str(v) for v in comp["source_tools"]))
            if comp.get("also_generated_by"):
                mol.SetProp(
                    "also_generated_by",
                    ",".join(str(v) for v in comp["also_generated_by"]),
                )
            for property_name in LIGAND_LINEAGE_PROPERTIES:
                if mol.HasProp(property_name):
                    continue
                value = comp.get(property_name)
                if value is None:
                    continue
                if isinstance(value, bool):
                    encoded = str(value).lower()
                elif isinstance(value, (dict, list, tuple)):
                    encoded = json.dumps(
                        value, sort_keys=True, separators=(",", ":")
                    )
                else:
                    encoded = str(value)
                mol.SetProp(property_name, encoded)

            activity = comp.get("activity_nm")
            mol.SetProp("activity_nM", str(activity) if activity is not None else "")
            activity_properties = {
                "standard_relation": comp.get("standard_relation"),
                "activity_lower_bound_nM": comp.get("activity_lower_bound_nm"),
                "activity_upper_bound_nM": comp.get("activity_upper_bound_nm"),
                "activity_lower_bound_inclusive": comp.get(
                    "activity_lower_bound_inclusive"
                ),
                "activity_upper_bound_inclusive": comp.get(
                    "activity_upper_bound_inclusive"
                ),
                "activity_is_censored": comp.get("activity_is_censored"),
                "activity_aggregation_method": comp.get("activity_aggregation_method"),
                "activity_measurement_count": comp.get("activity_measurement_count"),
                "activity_aggregation_measurement_count": comp.get(
                    "activity_aggregation_measurement_count"
                ),
                "activity_quarantined_measurement_count": comp.get(
                    "activity_quarantined_measurement_count"
                ),
                "activity_data_validity_quarantined_count": comp.get(
                    "activity_data_validity_quarantined_count"
                ),
                "activity_potential_duplicate_quarantined_count": comp.get(
                    "activity_potential_duplicate_quarantined_count"
                ),
                "activity_exact_measurement_count": comp.get(
                    "activity_exact_measurement_count"
                ),
                "activity_censored_measurement_count": comp.get(
                    "activity_censored_measurement_count"
                ),
            }
            for property_name, value in activity_properties.items():
                if value is not None:
                    encoded_value = (
                        str(value).lower() if isinstance(value, bool) else str(value)
                    )
                    mol.SetProp(property_name, encoded_value)
            list_properties = {
                "activity_assay_types": comp.get("activity_assay_types"),
                "activity_bao_formats": comp.get("activity_bao_formats"),
                "activity_document_chembl_ids": comp.get(
                    "activity_document_chembl_ids"
                ),
                "activity_eligible_assay_chembl_ids": comp.get(
                    "activity_eligible_assay_chembl_ids"
                ),
                "activity_quarantined_assay_chembl_ids": comp.get(
                    "activity_quarantined_assay_chembl_ids"
                ),
            }
            if comp.get("activity_assay_heterogeneity") is not None:
                mol.SetProp(
                    "activity_assay_heterogeneity",
                    str(bool(comp.get("activity_assay_heterogeneity"))).lower(),
                )
            for property_name, values in list_properties.items():
                if values:
                    mol.SetProp(
                        property_name,
                        ",".join(str(value) for value in values),
                    )
            if comp.get("activity_assay_strata") is not None:
                mol.SetProp(
                    "activity_assay_strata_json",
                    json.dumps(
                        comp["activity_assay_strata"],
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                )
            if comp.get("activity_measurements") is not None:
                mol.SetProp(
                    "activity_measurements_json",
                    json.dumps(
                        comp["activity_measurements"],
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                )
            if comp.get("assay_chembl_ids"):
                mol.SetProp(
                    "assay_chembl_ids",
                    ",".join(str(value) for value in comp["assay_chembl_ids"]),
                )
            if comp.get("chembl_ids"):
                mol.SetProp(
                    "chembl_ids",
                    ",".join(str(value) for value in comp["chembl_ids"]),
                )

            chembl_id = comp.get("chembl_id", "")
            mol.SetProp("chembl_id", str(chembl_id))

            try:
                mol_noH = Chem.RemoveHs(mol)
                mol.SetProp("MW", f"{Descriptors.ExactMolWt(mol_noH):.2f}")
                mol.SetProp("LogP", f"{Descriptors.MolLogP(mol_noH):.2f}")
                mol.SetProp("HBD", str(rdMolDescriptors.CalcNumHBD(mol_noH)))
                mol.SetProp("HBA", str(rdMolDescriptors.CalcNumHBA(mol_noH)))
                mol.SetProp("RotBonds", str(rdMolDescriptors.CalcNumRotatableBonds(mol_noH)))
                mol.SetProp("TPSA", f"{Descriptors.TPSA(mol_noH):.2f}")
            except Exception as exc:
                failed_compounds.append({**comp, "conformer_error": str(exc)})
                tracker.record(comp, "conformer.descriptors", "descriptor_calculation_failed", str(exc),
                               exception=None if is_rdkit_parse_error(exc) else exc)
                continue

            writer.write(mol)
            comp["conformer_written"] = True
            if raw_pose_preserved:
                raw_poses_preserved += 1
            else:
                embedded_conformers += 1
            total_written += 1
            written_compound_ids.append(str(comp_id))
    finally:
        writer.close()

    return {
        "success": total_written > 0,
        "output_sdf": str(output_sdf),
        "total_written": total_written,
        "failed_embedding": len(failed_compounds),
        "written_compound_ids": written_compound_ids,
        "raw_poses_preserved": raw_poses_preserved,
        "raw_poses_rejected": raw_poses_rejected,
        "embedded_conformers": embedded_conformers,
        "optimization_warnings": optimization_warnings,
        "failed_compounds": [
            {**{key: c[key] for key in ("compound_id", "input_compound_id", "parent_entity_id", "source_record_id", "candidate_role", "smiles", "smiles_standardized") if key in c},
             "error": c.get("conformer_error", "")}
            for c in failed_compounds
        ],
    }


class _MoleculeCoordinateError(ValueError):
    def __init__(self, code, reason):
        self.code = code
        super().__init__(reason)


def _record_coordinate_identity(comp: dict, mol, smiles: str) -> None:
    """Check the written graph and record a conformer, without enumerating states."""
    from rdkit import Chem

    # Match what downstream SDF readers will see, including stereochemistry
    # inferred from 3D. Unassigned input stereo remains a recorded uncertainty.
    restored = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), removeHs=False)
    expected = Chem.MolFromSmiles(smiles)
    if (restored is None or expected is None or not restored.GetNumConformers()
            or not restored.GetConformer().Is3D()):
        raise _MoleculeCoordinateError("conformer_invalid", "cannot read generated 3D conformer")
    positions = _usable_coordinates(restored)
    actual = Chem.RemoveHs(restored)
    expected = Chem.RemoveHs(expected)
    if not matches_declared_molecule(actual, expected):
        raise _MoleculeCoordinateError(
            "conformer_identity_mismatch", "conformer changed connectivity, charge, isotope or assigned stereochemistry",
        )
    unassigned = any(str(info.specified) == "Unspecified"
                     for info in Chem.FindPotentialStereo(expected))
    comp["coordinate_smiles"] = Chem.MolToSmiles(actual, isomericSmiles=True)
    comp["stereochemistry_status"] = "unassigned_input" if unassigned else "specified_or_not_applicable"
    if unassigned:
        comp["preparation_warnings_json"] = [
            "Input stereochemistry is unspecified; this conformer is not a stereochemical assignment."
        ]
    # This identifies the original ligand conformer, not later docked poses.
    payload = {"smiles": comp["coordinate_smiles"], "coordinates": positions}
    comp["conformer_coordinate_identity"] = "sha256:" + hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _usable_coordinates(mol):
    positions = [tuple(p) for p in mol.GetConformer().GetPositions()]
    if not all(math.isfinite(v) for xyz in positions for v in xyz):
        raise _MoleculeCoordinateError("conformer_invalid", "non-finite conformer coordinates")
    if len(positions) > 1 and all(math.dist(positions[0], p) <= 0.001 for p in positions[1:]):
        raise _MoleculeCoordinateError("conformer_invalid", "collapsed conformer coordinates")
    if any(math.dist(positions[b.GetBeginAtomIdx()], positions[b.GetEndAtomIdx()]) <= 0.001
           for b in mol.GetBonds()):
        raise _MoleculeCoordinateError("conformer_invalid", "zero-length conformer bond")
    return positions


def _prepare_coordinates(comp, smiles, optimize, random_seed):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    smiles_mol = Chem.MolFromSmiles(smiles) if smiles else None
    if smiles_mol is None:
        raise _MoleculeCoordinateError("molecule_parse_failed", "RDKit parse failed or empty SMILES")
    mol, raw_pose_error = _compatible_raw_pose(comp, smiles_mol)
    if mol is not None:
        return mol, True, comp.get("coordinate_origin", "generator_raw"), None, None
    mol = Chem.AddHs(smiles_mol)
    if not _embed_molecule(mol, random_seed):
        mol = Chem.AddHs(smiles_mol)
        if not _embed_molecule(mol, random_seed, use_random_coords=True):
            raise _MoleculeCoordinateError("conformer_embedding_failed", "embedding failed")
    warning = None
    if optimize:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            status, method = AllChem.MMFFOptimizeMolecule(mol, maxIters=500), "MMFF"
        else:
            status, method = AllChem.UFFOptimizeMolecule(mol, maxIters=500), "UFF"
        if status != 0:
            warning = {"method": method, "status": status}
    return mol, False, "rdkit_etkdg", raw_pose_error, warning


def _compatible_raw_pose(comp: dict, smiles_mol):
    """Return a generator pose only when it matches the standardized graph."""
    from rdkit import Chem

    molblock = comp.get("raw_pose_molblock")
    if not molblock:
        return None, "no raw pose"
    try:
        raw_mol = Chem.MolFromMolBlock(
            str(molblock), removeHs=False, sanitize=True, strictParsing=True
        )
    except Chem.MolSanitizeException as exc:
        return None, f"raw pose parse failed: {exc}"
    if raw_mol is None:
        return None, "raw pose parse failed"
    if raw_mol.GetNumConformers() < 1:
        return None, "raw pose has no conformer"
    if not raw_mol.GetConformer().Is3D():
        return None, "raw input contains 2D coordinates; rebuild a 3D conformer"

    try:
        _usable_coordinates(raw_mol)
        compatible = matches_declared_molecule(raw_mol, smiles_mol)
    except _MoleculeCoordinateError as exc:
        return None, str(exc)
    except Chem.MolSanitizeException as exc:
        return None, f"raw pose graph comparison failed: {exc}"
    if not compatible:
        return None, "raw pose graph differs after standardization"
    return raw_mol, None


def _next_compound_id(comp: dict, counters: dict[str, int]) -> str:
    # Keep IDs assigned before preparation, including gaps left by failed
    # records. Renumbering survivors would break their source lineage.
    if comp.get("compound_id"):
        return str(comp["compound_id"])
    source = comp.get("source", "")
    role = comp.get("candidate_role")
    if role == "active" or source == "active":
        counters["active"] += 1
        return f"ACT_{counters['active']:04d}"
    if role == "generated":
        counters["generated"] += 1
        return f"GEN_{counters['generated']:04d}"
    if role == "validation_decoy" or source in ("decoy", "dude", "rdkit_property_matched"):
        counters["decoy"] += 1
        return f"DEC_{counters['decoy']:04d}"
    return comp.get("compound_id", comp.get("name", f"MOL_{sum(counters.values()) + 1:04d}"))


def _embed_molecule(mol, seed: int = 42, use_random_coords: bool = False) -> bool:
    """Embed a single molecule with ETKDGv3 parameters."""
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if use_random_coords:
        params.useRandomCoords = True

    return AllChem.EmbedMolecule(mol, params) == 0

"""Preparation phase contracts shared by the validator and both skill runners.

Geometry QC remains the existing inexpensive audit. Passing this module's
checks does not certify bond valence, clashes, or repaired-loop geometry.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

from ligand_geometry import audit_sdf_geometry, geometry_failure_message

STRUCTURE_SCRIPTS = Path(__file__).resolve().parents[1] / "structure-prep" / "scripts"
DEFERRED_TOOLS = {
    "water_selection", "pdb2pqr", "protoss", "ligand_microstates", "openbabel", "meeko",
}


def preparation_module(name: str):
    """Load the existing wrappers without polluting generic module names."""
    spec = importlib.util.spec_from_file_location(
        f"sbdd_preparation_{name}", STRUCTURE_SCRIPTS / f"{name}.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(STRUCTURE_SCRIPTS))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(STRUCTURE_SCRIPTS))
    return module


def file_identity(path: Path) -> dict:
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def write_json(path: Path, value: dict) -> str:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str))
    return str(path)


def is_hydrogen(line: str) -> bool:
    element = line[76:78].strip().upper()
    if element:
        return element in {"H", "D"}
    return line[12:16].strip().lstrip("0123456789").startswith(("H", "D"))


def write_noh(source: Path, output: Path, *, remove_waters: bool = False) -> None:
    lines = []
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            if is_hydrogen(line) or (remove_waters and line[17:20].strip() in {"HOH", "WAT", "H2O"}):
                continue
        # Connectivity/serial references from the source may refer to removed atoms.
        if not line.startswith(("CONECT", "MASTER")):
            lines.append(line)
    output.write_text("".join(lines))


def water_inventory(source: Path, output: Path) -> str:
    waters = [line for line in source.read_text().splitlines()
              if line.startswith(("ATOM  ", "HETATM"))
              and line[17:20].strip() in {"HOH", "WAT", "H2O"}]
    return write_json(output, {
        "schema_version": 1, "stage": "structure_prep", "selection": "deferred_to_docking_prep",
        "source": file_identity(source), "water_atom_records": waters,
        "water_residue_count": len({line[21:27] for line in waters}),
    })


def docking_prep_errors(tools: dict) -> list[str]:
    errors = []
    prep = tools.get("docking_prep") or {}
    params = prep.get("params") or {}
    if prep.get("enabled") is not True:
        errors.append("docking.docking_prep must be explicitly enabled before the engine")
    for name in ("reason", "water_reason", "ligand_state_reason"):
        if not str(params.get(name) or "").strip():
            errors.append(f"docking.docking_prep.params.{name} must contain a justification")
    backends = [name for name in ("reduce", "pdb2pqr", "protoss")
                if (tools.get(name) or {}).get("enabled") is True]
    if len(backends) != 1:
        errors.append("docking_prep requires exactly one of reduce/pdb2pqr/protoss; no automatic fallback")
    if "protoss" in backends:
        errors.append("docking.protoss is unavailable: the current PDB-ID API cannot preserve the prepared complex; prepared-complex upload is required")
    ph = params.get("ph")
    if "pdb2pqr" in backends or (tools.get("ligand_microstates") or {}).get("enabled"):
        if not isinstance(ph, (int, float)) or isinstance(ph, bool) or not 0 <= ph <= 14:
            errors.append("docking_prep requires an explicit pH in [0,14] for pdb2pqr/ligand_microstates")
    if "reduce" in backends and (tools["reduce"].get("params") or {}).get("build") is False:
        errors.append("docking.reduce.build must be true when starting from prepared_receptor_noH")
    water_policy = params.get("water_policy")
    selected = (tools.get("water_selection") or {}).get("enabled") is True
    if water_policy not in {"dry", "selected"} or selected != (water_policy == "selected"):
        errors.append("docking_prep.water_policy must be dry or selected, matching docking.water_selection.enabled")
    if params.get("ligand_state_policy", "preserve_declared") != "preserve_declared":
        errors.append("docking_prep supports preserve_declared ligand states; supply any site-specific preparations as an explicit ligand library")
    vina = (tools.get("autodock_vina") or {}).get("enabled") is True
    meeko = (tools.get("meeko") or {}).get("enabled") is True
    if vina and not meeko:
        errors.append("AutoDock Vina requires docking.meeko after final receptor preparation")
    if meeko and not vina:
        errors.append("docking.meeko is only consumed by Vina; UniDock owns its native format preparation")
    if meeko and (tools["meeko"].get("params") or {}).get("allow_bad_residues") is True:
        errors.append("docking.meeko.allow_bad_residues must be false; do not hide preparation failures")
    if (tools.get("openbabel") or {}).get("enabled"):
        errors.append("docking.openbabel has no selected format consumer in this route; keep it disabled")
    return errors


def run_docking_prep(input_data: dict, output_dir: Path) -> dict:
    """Run final preparation once; the caller invokes its engine only on success."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tools = input_data.get("tools") or {}
    report = {"schema_version": 1, "stage": "docking_prep", "success": False,
              "qc_scope": "basic structure/file checks; geometric/valence QC pending"}
    files = {}
    try:
        errors = docking_prep_errors(tools)
        if errors:
            raise ValueError("; ".join(errors))
        params = tools["docking_prep"].get("params") or {}
        report["plan"] = params
        # Reject unusable coordinates before protonation or any docking engine.
        ligand = Path(input_data["ligand_file"])
        report["ligand_geometry_qc"] = audit_sdf_geometry(ligand)
        if not report["ligand_geometry_qc"]["success"]:
            raise ValueError(geometry_failure_message(report["ligand_geometry_qc"]))
        library = output_dir / "docking_species_library.sdf"
        ligand_records = _declared_ligand_records(ligand, output_path=library)
        source = Path(input_data.get("receptor_noH_file") or input_data["receptor_file"])
        if any(is_hydrogen(line) for line in source.read_text().splitlines()
               if line.startswith(("ATOM  ", "HETATM"))):
            raise ValueError("prepared_receptor_noH contains hydrogens; regenerate structure_prep with add_hydrogens=false")
        source_identity = file_identity(source)
        report["inputs"] = {"prepared_receptor_noH": source_identity}
        for key in ("prep_manifest_file", "pocket_definition_file", "docking_box_file",
                    "reference_ligand_file", "water_inventory_file", "interaction_blueprint_file"):
            if input_data.get(key):
                report["inputs"][key] = file_identity(Path(input_data[key]))
        report["search_space_source"] = input_data.get("search_space_source")
        report["docking_box_inline"] = input_data.get("docking_box")
        audit = preparation_module("prep_utils").audit_structure
        report["input_qc"] = audit(source, require_hydrogens=False)
        if not report["input_qc"]["passed"]:
            raise ValueError("docking_prep input basic QC failed")
        backend_input = output_dir / "receptor_noH_for_docking.pdb"
        write_noh(source, backend_input, remove_waters=True)
        water = {"policy": params["water_policy"], "reason": params["water_reason"], "selected_count": 0}
        if params["water_policy"] == "selected":
            inventory_path = input_data.get("water_inventory_file")
            if not inventory_path:
                raise ValueError("selected water policy requires the original water_inventory")
            inventory = json.loads(Path(inventory_path).read_text())
            pdb_utils = preparation_module("pdb_utils")
            water_source = output_dir / "source_waters.pdb"
            water_source.write_text("\n".join(inventory["water_atom_records"]) + "\n")
            atoms = pdb_utils.parse_pdb_atoms(backend_input) + pdb_utils.parse_pdb_atoms(water_source)
            ligand_key = None
            if input_data.get("reference_ligand_file"):
                from rdkit import Chem
                mol = next(iter(Chem.SDMolSupplier(input_data["reference_ligand_file"], removeHs=False)), None)
                if mol is None or not mol.GetNumConformers():
                    raise ValueError("Water selection requires a readable 3D reference ligand")
                ligand_key = ("LIG", "~", -1, "")
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() == 1:
                        continue
                    xyz = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                    atoms.append({"record": "HETATM", "res_name": "LIG", "chain_id": "~", "res_seq": -1,
                                  "i_code": "", "x": xyz.x, "y": xyz.y, "z": xyz.z})
            water_params = tools["water_selection"].get("params") or {}
            if ligand_key is None and not water_params.get("include_residues"):
                raise ValueError("Water selection without a reference ligand requires explicit include_residues")
            water.update(pdb_utils.select_binding_site_waters(atoms, ligand_key, **water_params))
            selected_atoms = water.pop("selected_atoms")
            base_lines = [line for line in backend_input.read_text().splitlines(keepends=True)
                          if not line.startswith("END")]
            backend_input.write_text("".join(base_lines) + "".join(a["raw_line"] for a in selected_atoms) + "END\n")
        files["water_report"] = write_json(output_dir / "water_report.json", water)
        backend = next(name for name in ("reduce", "pdb2pqr") if tools.get(name, {}).get("enabled"))
        backend_config = tools[backend]
        backend_params = dict(backend_config.get("params") or {})
        runtime = backend_config.get("runtime") or {}
        external = preparation_module("external_tools")
        final = output_dir / "prepared_receptor.pdb"
        if backend == "reduce":
            backend_result = external.run_reduce2(backend_input, final, runtime, **backend_params)
        else:
            backend_result = external.run_pdb2pqr(backend_input, final, output_dir / "receptor.pqr", runtime,
                                                 ph=params["ph"], **backend_params)
        report["backend_result"] = backend_result
        report["output_qc"] = audit(final, require_hydrogens=True)
        if not report["output_qc"]["passed"]:
            raise ValueError("Final receptor basic QC failed: " + "; ".join(report["output_qc"]["issues"]))
        final_waters = {line[21:27] for line in final.read_text().splitlines()
                        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() in {"HOH", "WAT", "H2O"}}
        if len(final_waters) != water["selected_count"]:
            raise ValueError("Protonation backend changed the selected water inventory; review the backend output")
        if file_identity(source) != source_identity:
            raise ValueError("Canonical noH source changed during docking_prep")
        files["prepared_receptor"] = str(final)
        blueprint_review = "requires_final_interaction_recalculation" if input_data.get("interaction_blueprint_file") else "not_supplied"
        files["receptor_state_manifest"] = write_json(output_dir / "receptor_state_manifest.json", {
            "schema_version": 1, "backend": backend, "ph": params.get("ph"),
            "ph_applied": backend == "pdb2pqr", "reason": params["reason"],
            "source": source_identity, "prepared_receptor": file_identity(final),
            "runtime": runtime, "parameters": backend_params, "water_report": file_identity(Path(files["water_report"])),
            "interaction_blueprint_status": blueprint_review,
        })
        if tools.get("meeko", {}).get("enabled"):
            pdbqt = output_dir / "receptor.pdbqt"
            meeko = tools["meeko"]
            report["meeko"] = external.run_meeko_receptor(final, pdbqt, meeko.get("runtime") or {}, **(meeko.get("params") or {}))
            files["receptor_pdbqt"] = str(pdbqt)
        # The annotated library preserves supplied structures and record order.
        files["docking_species_library"] = str(library)
        files["microstate_identity_map"] = write_json(output_dir / "microstate_identity_map.json", {
            "schema_version": 1, "policy": "preserve_declared", "reason": params["ligand_state_reason"],
            "source": file_identity(ligand), "output": file_identity(library),
            "records": ligand_records,
        })
        if tools.get("ligand_microstates", {}).get("enabled"):
            from rdkit import Chem
            reference = input_data.get("reference_ligand_file")
            if not reference:
                raise ValueError("Reference ligand microstates require reference_ligand_file")
            molecule = next(iter(Chem.SDMolSupplier(reference)), None)
            if molecule is None:
                raise ValueError("Cannot read reference ligand for microstate diagnostics")
            microstates = output_dir / "reference_ligand_microstates.smi"
            conf = tools["ligand_microstates"]
            report["reference_microstates"] = external.run_dimorphite(Chem.MolToSmiles(molecule), microstates,
                conf.get("runtime") or {}, ph=params["ph"], **(conf.get("params") or {}))
            files["reference_ligand_microstates"] = str(microstates)
        report["success"] = True
    except Exception as exc:
        report["error"] = str(exc)
    files["docking_prep_qc"] = str(output_dir / "docking_prep_qc.json")
    report["files"] = files
    write_json(Path(files["docking_prep_qc"]), report)
    return report


def _declared_ligand_records(path: Path, *, output_path: Path | None = None) -> list[dict]:
    """Preserve ordinary library identities without requiring state enumeration."""
    from rdkit import Chem
    from pose_identity import LIGAND_LINEAGE_PROPERTIES, complete_ligand_ids

    records = []
    seen = set()
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
    complete_ligand_ids(molecules, source_token=hashlib.sha256(path.read_bytes()).hexdigest())
    for index, molecule in enumerate(molecules):
        if molecule is None:
            raise ValueError(f"Cannot read ligand record {index}")
        compound_id = next((molecule.GetProp(key).strip() for key in ("compound_id", "_Name")
                            if molecule.HasProp(key) and molecule.GetProp(key).strip()), None)
        if not compound_id or compound_id in seen:
            raise ValueError(
                f"Missing or duplicate ligand compound ID {compound_id!r} at record {index}; "
                "supply a unique compound_id per docking species and preserve parent/state IDs separately"
            )
        seen.add(compound_id)
        # GNINA and other native tools join scores by the SDF title. Give every
        # engine the same unique key while retaining the supplied display name.
        molecule.SetProp("_Name", compound_id)
        record = {key: molecule.GetProp(key) for key in LIGAND_LINEAGE_PROPERTIES
                  if molecule.HasProp(key)}
        record.update(compound_id=compound_id, source_molecule_index=index,
                      microstate_id=record.get("microstate_id") or None,
                      conformer_id=record.get("conformer_id") or None,
                      coordinate_smiles=Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True))
        records.append(record)
    if output_path is not None:
        writer = Chem.SDWriter(str(output_path))
        try:
            for molecule in molecules:
                writer.write(molecule)
        finally:
            writer.close()
    return records

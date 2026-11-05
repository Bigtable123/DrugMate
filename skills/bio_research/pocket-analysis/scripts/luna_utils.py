"""LUNA 0.14 core-API wrapper for interaction and hashed IFP extraction."""
from __future__ import annotations

import math
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any


@contextmanager
def _heavy_atom_property_targets():
    """LUNA 0.14 passes explicit ligand H atoms to a heavy-atom-only mapper.

    Keep the molecular graph and hydrogens unchanged. Normalize only the target
    list used for property attribution, then restore the class method. The
    existing native worker isolates this compatibility patch from other tools.
    """
    from luna.mol.groups import AtomGroupPerceiver

    original = AtomGroupPerceiver._assign_properties

    def assign(self, mol_obj, target_atoms=None):
        if target_atoms is not None:
            target_atoms = [atom for atom in target_atoms if atom.element != "H"]
        return original(self, mol_obj, target_atoms)

    AtomGroupPerceiver._assign_properties = assign
    try:
        yield
    finally:
        AtomGroupPerceiver._assign_properties = original


_TYPE_MAP = {
    "Hydrogen bond": "HBond",
    "Weak hydrogen bond": "WeakHBond",
    "Hydrophobic": "Hydrophobic",
    "Van der Waals": "VdW",
    "Ionic": "SaltBridge",
    "Cation-pi": "PiCation",
    "Face-to-face pi-stacking": "PiStacking",
    "Displaced face-to-face pi-stacking": "PiStacking",
    "Edge-to-face pi-stacking": "EdgeToFace",
    "Halogen bond": "HalogenBond",
    "Metal coordination": "MetalComplex",
    "Water-bridged hydrogen bond": "WaterBridge",
}


def _version(distribution: str) -> str | None:
    try:
        from importlib.metadata import version

        return version(distribution)
    except Exception:
        return None


def _finite_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _atom_payload(atom: Any, *, ligand: bool) -> dict:
    data = atom.as_json()
    name_value = data.get("name") or ["", ""]
    name = name_value[0] if isinstance(name_value, (list, tuple)) else name_value
    payload = {
        "name": str(name or ""),
        "namespace": "luna_biopython_serial",
    }
    serial = getattr(getattr(atom, "atom", None), "serial_number", None)
    if serial is not None:
        try:
            serial = int(serial)
            payload["index"] = serial
            if ligand:
                payload["input_order"] = serial
                payload["residue_index"] = serial - 1
        except (TypeError, ValueError):
            pass
    return payload


def _residue_payload(group: Any) -> dict | None:
    atoms = list(group.atoms)
    if not atoms:
        return None
    data = atoms[0].as_json()
    resid = data.get("resid") or ["", None, ""]
    try:
        number = int(resid[1])
    except (IndexError, TypeError, ValueError):
        return None
    icode = ""
    if len(resid) > 2:
        icode = str(resid[2] or "").strip().replace("\x00", "")
    return {
        "resname": str(data.get("resname") or "").upper(),
        "resnum": number,
        "chain": str(data.get("chain") or ""),
        "icode": icode,
    }


def _group_features(group: Any) -> set[str]:
    return {str(getattr(feature, "name", feature)).lower() for feature in (group.features or [])}


def _roles(raw_type: str, ligand_is_source: bool, ligand_group: Any, protein_group: Any) -> tuple[str | None, str | None]:
    if raw_type in {"Hydrogen bond", "Weak hydrogen bond"}:
        return ("donor", "acceptor") if ligand_is_source else ("acceptor", "donor")
    if raw_type == "Water-bridged hydrogen bond":
        return ("donor", "acceptor") if ligand_is_source else ("acceptor", "donor")
    ligand_features = _group_features(ligand_group)
    protein_features = _group_features(protein_group)
    if raw_type == "Ionic":
        if any("pos" in feature or "cation" in feature for feature in ligand_features):
            return "cation", "anion"
        if any("neg" in feature or "anion" in feature for feature in ligand_features):
            return "anion", "cation"
    if raw_type == "Cation-pi":
        ligand_is_cation = any("pos" in feature or "cation" in feature for feature in ligand_features)
        protein_is_cation = any("pos" in feature or "cation" in feature for feature in protein_features)
        if ligand_is_cation and not protein_is_cation:
            return "cation", "aromatic"
        if protein_is_cation and not ligand_is_cation:
            return "aromatic", "cation"
    return None, None


def _geometry(params: dict) -> dict:
    geometry = {}
    preferred_distances = (
        ("da_dist_hb_inter", "heavy_atom_distance"),
        ("ha_dist_hb_inter", "acceptor_hydrogen_distance"),
        ("dist_hydrop_inter", "distance"),
        ("dist_atom_atom", "distance"),
        ("dist_attract_inter", "distance"),
        ("dist_cation_pi_inter", "distance"),
        ("cc_dist_pi_pi_inter", "distance"),
        ("xa_dist_xbond_inter", "distance"),
        ("ma_dist_metal_coord", "distance"),
    )
    for source, target in preferred_distances:
        value = _finite_number(params.get(source))
        if value is not None:
            geometry[target] = value
    angle_map = {
        "dha_ang_hb_inter": "donor_angle",
        "har_ang_hb_inter": "acceptor_angle",
        "dar_ang_hb_inter": "angle",
        "dihed_ang_pi_pi_inter": "dihedral_angle",
        "disp_ang_pi_pi_inter": "displacement_angle",
        "cxa_ang_xbond_inter": "donor_angle",
        "xar_ang_xbond_inter": "acceptor_angle",
    }
    for source, target in angle_map.items():
        value = _finite_number(params.get(source))
        if value is not None:
            geometry[target] = value
    if "distance" not in geometry and "heavy_atom_distance" not in geometry:
        for key, raw in params.items():
            if "dist" in str(key).lower() and (value := _finite_number(raw)) is not None:
                geometry["distance"] = value
                break
    return geometry


def _ligand_explicit_hydrogens(ligand_file: Path) -> int:
    from rdkit import Chem

    mol = next(
        (
            m
            for m in Chem.SDMolSupplier(
                str(ligand_file), removeHs=False, sanitize=False,
            )
            if m
        ),
        None,
    )
    if mol is None:
        raise ValueError(f"RDKit failed to parse LUNA ligand: {ligand_file}")
    return sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())


def parse_luna_interactions(interactions: list[Any]) -> list[dict]:
    details = []
    for interaction in interactions:
        source_ref = bool(interaction.src_grp.has_reference())
        target_ref = bool(interaction.trgt_grp.has_reference())
        if source_ref == target_ref:
            continue
        raw_type = str(interaction.type)
        mapped_type = _TYPE_MAP.get(raw_type)
        if mapped_type is None:
            continue
        ligand_is_source = source_ref
        ligand_group = interaction.src_grp if source_ref else interaction.trgt_grp
        protein_group = interaction.trgt_grp if source_ref else interaction.src_grp
        residue = _residue_payload(protein_group)
        if residue is None:
            continue
        ligand_role, protein_role = _roles(
            raw_type, ligand_is_source, ligand_group, protein_group,
        )
        weak = raw_type == "Weak hydrogen bond"
        details.append({
            "type": mapped_type,
            "interaction_subtype": raw_type,
            "residue_record": residue,
            "ligand_role": ligand_role,
            "protein_role": protein_role,
            "ligand_atoms": [_atom_payload(atom, ligand=True) for atom in ligand_group.atoms],
            "protein_atoms": [_atom_payload(atom, ligand=False) for atom in protein_group.atoms],
            "geometry": _geometry(dict(interaction.params or {})),
            "mediator": "water" if mapped_type == "WaterBridge" else "metal" if mapped_type == "MetalComplex" else "direct",
            "blueprint_eligible": not weak,
            "quality_flags": ["weak_hbond"] if weak else [],
            "raw_metadata": {
                "luna_params": dict(interaction.params or {}),
                "directional": bool(interaction.is_directional()),
                "ligand_is_source": ligand_is_source,
            },
        })
    return details


def run_luna_fingerprint(
    protein_file: Path,
    ligand_file: Path,
    ligand_resname: str,
    *,
    ifp_length: int = 4096,
    ifp_type: str = "EIFP",
    add_dependent_interactions: bool = False,
) -> dict:
    from luna.interaction.calc import InteractionCalculator
    from luna.interaction.filter import InteractionFilter
    from luna.interaction.fp.type import IFPType
    from luna.mol.entry import MolFileEntry
    from luna.projects import LocalProject

    if protein_file.suffix.lower() != ".pdb":
        raise ValueError("LUNA wrapper requires a PDB receptor")
    fingerprint_type = getattr(IFPType, str(ifp_type).upper(), IFPType.EIFP)
    entry = MolFileEntry.from_mol_file(
        protein_file.stem,
        ligand_resname or "ligand",
        str(ligand_file),
        is_multimol_file=False,
        mol_obj_type="rdkit",
        overwrite_mol_name=True,
    )
    calculator = InteractionCalculator(
        inter_filter=InteractionFilter.new_pli_filter(),
        add_cov=False,
        add_dependent_inter=bool(add_dependent_interactions),
    )
    with tempfile.TemporaryDirectory(prefix="luna_step1_") as tmp:
        project = LocalProject(
            entries=[entry],
            working_path=tmp,
            pdb_path=str(protein_file.parent),
            overwrite_path=True,
            add_h=False,
            inter_calc=calculator,
            calc_ifp=True,
            ifp_length=int(ifp_length),
            ifp_type=fingerprint_type,
            out_pse=False,
            nproc=None,
            verbosity=0,
            logging_enabled=False,
        )
        with _heavy_atom_property_targets():
            project.run()
        if project.errors:
            raise RuntimeError("; ".join(str(error) for error in project.errors))
        results = project.get_entry_results(entry)
        if results is None:
            raise RuntimeError("LUNA produced no EntryResults")
        details = parse_luna_interactions(list(results.interactions_mngr.interactions))
        fingerprint = results.ifp
        counts = {
            str(int(index)): int(count)
            for index, count in dict(fingerprint.counts).items()
        } if fingerprint is not None else {}

    explicit_hydrogens = _ligand_explicit_hydrogens(ligand_file)
    capability_warnings = []
    if explicit_hydrogens == 0:
        capability_warnings.append("ligand_donor_hbond_not_evaluable_without_explicit_hydrogens")
    return {
        "success": True,
        "tool": "LUNA",
        "tool_version": _version("luna"),
        "num_interactions_total": len(details),
        "interactions_detail": details,
        "interaction_fingerprint": {
            "type": str(ifp_type).upper(),
            "length": int(ifp_length),
            "count_fingerprint": True,
            "counts": counts,
            "nonzero_count": len(counts),
        },
        "quality_flags": capability_warnings,
        "chemical_representation": {
            "protein_file": str(protein_file),
            "ligand_file": str(ligand_file),
            "ligand_explicit_hydrogen_count": explicit_hydrogens,
        },
        "capabilities": {
            "atom_identity": True,
            "ligand_input_order": True,
            "directional_roles": True,
            "geometry": True,
            "ligand_donor_hbond": (
                "evaluable" if explicit_hydrogens else "not_evaluable"
            ),
            "pymol_visualization": False,
        },
    }

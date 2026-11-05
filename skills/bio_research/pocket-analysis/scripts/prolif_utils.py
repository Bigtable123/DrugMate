"""ProLIF utilities: compute protein-ligand interaction fingerprints."""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _protein_hydrogen_count(pdb_block: str) -> int:
    """Count explicit amino-acid hydrogens; waters/ligands cannot satisfy QC."""
    from Bio.PDB.Polypeptide import is_aa

    count = 0
    for line in pdb_block.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or not is_aa(line[17:20]):
            continue
        element = line[76:78].strip().upper()
        if not element:
            element = line[12:16].strip().lstrip("0123456789")[:1].upper()
        count += element in {"H", "D"}
    return count


def _validate_ligand_3d(coords, ligand_path: Path, conformer_is_3d: bool | None = None) -> dict:
    """Validate finite, non-degenerate coordinates without rejecting planar ligands.

    Planarity is rotation-dependent and chemically legitimate (aromatic ligands
    are often nearly planar), so it is diagnostic only and never a fatal 2D test.
    """
    import numpy as np
    coords = np.asarray(coords, dtype=float)
    finite = bool(coords.ndim == 2 and coords.shape[1:] == (3,) and np.isfinite(coords).all())
    if not finite or len(coords) == 0:
        return {
            "is_3d": False,
            "coordinates_valid": False,
            "reason": "coordinates missing, malformed, or non-finite",
        }
    z_std = float(np.std(coords[:, 2]))
    coord_range = coords.max(axis=0) - coords.min(axis=0)
    centered = coords - coords.mean(axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    coordinate_rank = int(sum(value > 1e-3 for value in singular_values))
    is_planar = coordinate_rank <= 2
    all_z_near_zero = bool(np.abs(coords[:, 2]).max() < 0.5)
    has_spatial_extent = bool(float(np.linalg.norm(coord_range)) > 0.1)

    diagnostics = {
        "is_3d": has_spatial_extent,
        "coordinates_valid": has_spatial_extent,
        "conformer_marked_3d": conformer_is_3d,
        "planar_geometry": is_planar,
        "coordinate_rank": coordinate_rank,
        "z_std": round(z_std, 4),
        "coord_range_xyz": [round(float(r), 2) for r in coord_range],
        "all_z_near_zero": all_z_near_zero,
    }

    if is_planar:
        logger.info(
            "Ligand %s is planar in the supplied pose (rank=%d); this is diagnostic, not an error",
            ligand_path.name,
            coordinate_rank,
        )
    return diagnostics


def _validate_ligand_proximity(
    lig_coords, protein_file: Path, max_contact_distance: float = 12.0
) -> dict:
    """Check nearest ligand-protein atoms, not distance to whole-protein center."""
    import numpy as np
    lig_coords = np.asarray(lig_coords, dtype=float)
    prot_coords = []
    with open(protein_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    prot_coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue

    if not prot_coords:
        return {"proximity_valid": None, "reason": "no protein atoms found"}

    protein = np.asarray(prot_coords, dtype=float)
    min_distance = float(np.sqrt(np.min(np.sum(
        (protein[:, None, :] - lig_coords[None, :, :]) ** 2,
        axis=2,
    ))))

    return {
        "proximity_valid": min_distance <= max_contact_distance,
        "minimum_ligand_protein_atom_distance": round(min_distance, 2),
        "max_contact_distance": max_contact_distance,
        "method": "minimum_interatomic_distance",
    }


def _nearby_protein_pdb_block(
    protein_file: Path,
    ligand_coordinates,
    cutoff: float,
    preserve_hetero_residues: list | None = None,
) -> tuple[str, int, bool]:
    """Return complete nearby residues using distance to any ligand atom."""
    import numpy as np

    ligand_coordinates = np.asarray(ligand_coordinates, dtype=float)
    cutoff_sq = float(cutoff) ** 2
    residues: dict[tuple[str, str, str, str], list[str]] = {}
    selected_keys = set()
    preserve_keys = _hetero_key_set(preserve_hetero_residues or [])

    with open(protein_file) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            key = (line[21:22], line[22:26], line[26:27], line[17:20])
            if line.startswith("HETATM") and _line_hetero_key(line) not in preserve_keys:
                continue
            residues.setdefault(key, []).append(line)
            try:
                coord = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            except ValueError:
                continue
            if float(np.min(np.sum((ligand_coordinates - coord) ** 2, axis=1))) <= cutoff_sq:
                selected_keys.add(key)

    if not residues:
        return "", 0, False
    if not selected_keys:
        selected_keys = set(residues)

    block_lines = []
    serial = 1
    for key in sorted(selected_keys, key=lambda item: (item[0], int(item[1]), item[2], item[3])):
        for line in residues[key]:
            block_lines.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        block_lines.append(f"TER   {serial:5d}\n")
        serial += 1
    block_lines.append("END\n")
    block = "".join(block_lines)
    return block, len(selected_keys), _protein_hydrogen_count(block) > 0


def _line_hetero_key(line: str) -> tuple[str, str, str, str]:
    return (
        line[21:22].strip(), line[22:26].strip(), line[26:27].strip(),
        line[17:20].strip().upper(),
    )


def _hetero_key_set(values: list) -> set[tuple[str, str, str, str]]:
    keys = set()
    for value in values:
        if isinstance(value, dict):
            chain = str(value.get("chain") or "").strip()
            resnum = str(value.get("resnum") or value.get("residue_number") or "").strip()
            icode = str(value.get("icode") or "").strip()
            resname = str(value.get("resname") or value.get("residue_name") or "").strip().upper()
            if resnum and resname:
                keys.add((chain, resnum, icode, resname))
        elif isinstance(value, str):
            # Canonical compact form: CHAIN:RESNAME123[ICODE]
            import re
            match = re.match(r"^(?:(?P<chain>[^:]+):)?(?P<name>[A-Za-z0-9]{1,3})(?P<num>-?\d+)(?P<icode>[A-Za-z]?)$", value.strip())
            if match:
                keys.add((
                    match.group("chain") or "", match.group("num"),
                    match.group("icode") or "", match.group("name").upper(),
                ))
    return keys


def run_prolif_fingerprint(
    protein_file: Path,
    ligand_file: Path,
    ligand_resname: str = "LIG",
    interaction_types: Optional[list[str]] = None,
    vicinity_cutoff: float = 6.5,
    preserve_hetero_residues: Optional[list] = None,
) -> dict:
    """Compute interaction fingerprint using ProLIF.

    Uses RDKit to load SDF ligands (MDAnalysis cannot parse SDF as topology).
    Limits protein to residues within vicinity_cutoff of ligand to avoid
    segfaults from large RDKit bond-guessing operations.

    Validates ligand 3D-ness and proximity before computing fingerprints.
    If validation fails, returns diagnostic info instead of silently returning empty.
    """
    import prolif
    from rdkit import Chem

    ligand_path = Path(ligand_file)

    # --- Load ligand via RDKit ---
    rdkit_mol = None
    if ligand_path.suffix.lower() == ".sdf":
        suppl = Chem.SDMolSupplier(str(ligand_path), removeHs=False)
        rdkit_mol = next((m for m in suppl if m is not None), None)
        if rdkit_mol is None:
            raise ValueError(f"RDKit failed to read any molecule from {ligand_path}")
        if not any(a.GetAtomicNum() == 1 for a in rdkit_mol.GetAtoms()):
            rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
        lig_mol = prolif.Molecule.from_rdkit(rdkit_mol, resname=ligand_resname)
    else:
        import MDAnalysis as mda

        u = mda.Universe(str(ligand_path))
        lig_sel = u.select_atoms(f"resname {ligand_resname}")
        if len(lig_sel) == 0:
            lig_sel = u.select_atoms("all")
        lig_mol = prolif.Molecule.from_mda(lig_sel, force=True)

    # --- Validate ligand 3D coordinates ---
    lig_coords = lig_mol.GetConformer().GetPositions()
    conformer_is_3d = bool(lig_mol.GetConformer().Is3D())
    coord_validation = _validate_ligand_3d(lig_coords, ligand_path, conformer_is_3d)
    if not coord_validation["is_3d"]:
        return {
            "residues": [],
            "interaction_types": [],
            "bitvector": [],
            "summary": {},
            "num_interacting_residues": 0,
            "num_interactions_total": 0,
            "validation_failed": True,
            "diagnostics": {
                "error": "ligand_coordinates_are_2d",
                "detail": (
                    f"Ligand '{ligand_path.name}' has invalid or degenerate coordinates. "
                    "Provide finite coordinates in the receptor coordinate frame."
                ),
                "coord_validation": coord_validation,
            },
        }

    # --- Validate ligand is near the protein ---
    centroid = lig_coords.mean(axis=0)
    proximity = _validate_ligand_proximity(lig_coords, Path(protein_file))
    if proximity.get("proximity_valid") is False:
        return {
            "residues": [],
            "interaction_types": [],
            "bitvector": [],
            "summary": {},
            "num_interacting_residues": 0,
            "num_interactions_total": 0,
            "validation_failed": True,
            "diagnostics": {
                "error": "ligand_not_near_protein",
                "detail": (
                    f"Nearest ligand/protein atoms are {proximity['minimum_ligand_protein_atom_distance']}Å apart "
                    f"(max allowed: {proximity['max_contact_distance']}Å). The ligand may have been re-embedded or "
                    f"its coordinates don't match the protein coordinate frame."
                ),
                "proximity": proximity,
            },
        }

    # --- Load protein near ligand without MDAnalysis->RDKit conversion ---
    # ProLIF's Molecule.from_mda path can segfault inside native RDKit splitting
    # for prepared/noH PDBs.  Build a small nearby-residue PDB block ourselves
    # and let RDKit read it directly, preserving PDB residue metadata.
    protein_block, nearby_residue_count, used_protein_h = _nearby_protein_pdb_block(
        Path(protein_file),
        lig_coords,
        vicinity_cutoff,
        preserve_hetero_residues,
    )
    if not protein_block:
        raise ValueError(f"No protein atoms found in {protein_file}")
    protein_input = {
        "protein_file_used": str(Path(protein_file).resolve()),
        "protein_h_used": used_protein_h,
        "protein_hydrogen_count": _protein_hydrogen_count(protein_block),
        "protein_hydrogen_count_scope": "selected_amino_acid_residues",
        "input_strategy": (
            "explicit_protein_hydrogens" if used_protein_h
            else "missing_protein_hydrogens"
        ),
    }
    if not used_protein_h:
        return {
            **protein_input,
            "status": "tool_failed",
            "tool": "ProLIF",
            "tool_version": getattr(prolif, "__version__", None),
            "validation_failed": True,
            "diagnostics": {
                "error": "protein_explicit_hydrogens_missing",
                "detail": (
                    "Selected protein residues contain no explicit hydrogens. "
                    "Use the separately prepared provisional receptor from "
                    "structure_prep.reduce; do not replace the noH geometry baseline."
                ),
                "nearby_residue_count": nearby_residue_count,
            },
        }

    prot_rdkit = Chem.MolFromPDBBlock(
        protein_block,
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    if prot_rdkit is None:
        return {
            **protein_input,
            "residues": [],
            "interaction_types": [],
            "bitvector": [],
            "summary": {},
            "num_interacting_residues": 0,
            "num_interactions_total": 0,
            "validation_failed": True,
            "diagnostics": {
                "error": "protein_rdkit_pdb_parse_failed",
                "detail": (
                    f"RDKit failed to parse nearby protein PDB block from {protein_file}. "
                    f"Nearby residue count: {nearby_residue_count}."
                ),
            },
        }

    try:
        prot_mol = prolif.Molecule(prot_rdkit)
    except Exception as e:
        return {
            **protein_input,
            "residues": [],
            "interaction_types": [],
            "bitvector": [],
            "summary": {},
            "num_interacting_residues": 0,
            "num_interactions_total": 0,
            "validation_failed": True,
            "diagnostics": {
                "error": "protein_prolif_molecule_failed",
                "detail": f"ProLIF failed to build protein Molecule from RDKit PDB block: {e}",
                "nearby_residue_count": nearby_residue_count,
            },
        }

    # --- Compute fingerprint ---
    fp = prolif.Fingerprint(interaction_types) if interaction_types else prolif.Fingerprint()
    fp.run_from_iterable([lig_mol], prot_mol, progress=False, n_jobs=1)

    df = fp.to_dataframe()

    # --- Parse results ---
    interaction_details = _serialize_ifp(fp, lig_mol, prot_mol)
    summary = {}
    residues = []
    all_interaction_types = set()

    if not df.empty:
        for col in df.columns:
            if isinstance(col, tuple) and len(col) >= 3:
                # ProLIF 2.x: (ligand_res, protein_res, interaction_type)
                res_id = str(col[1])
                itype = str(col[2])
            elif isinstance(col, tuple) and len(col) >= 2:
                res_id = str(col[0])
                itype = str(col[1])
            else:
                continue

            if res_id not in summary:
                summary[res_id] = []
                residues.append(res_id)

            if df[col].any():
                summary[res_id].append(itype)
                all_interaction_types.add(itype)

    bitvector = []
    itypes_sorted = sorted(all_interaction_types)
    for res in residues:
        row = [1 if it in summary.get(res, []) else 0 for it in itypes_sorted]
        bitvector.append(row)

    num_interacting = len([r for r in residues if summary.get(r)])
    num_total = len(interaction_details) or sum(len(v) for v in summary.values())

    result = {
        **protein_input,
        "status": "detected" if num_total else "not_detected",
        "tool": "ProLIF",
        "tool_version": getattr(prolif, "__version__", None),
        "residues": residues,
        "interaction_types": itypes_sorted,
        "bitvector": bitvector,
        "summary": summary,
        "interactions_detail": interaction_details,
        "num_interacting_residues": num_interacting,
        "num_interactions_total": num_total,
        "coord_validation": coord_validation,
        "preserved_hetero_residues": preserve_hetero_residues or [],
        "capabilities": {
            "enabled_pairwise_interaction_types": sorted(
                str(name) for name in fp.interactions.keys()
            ),
            "enabled_bridged_interaction_types": sorted(
                str(name) for name in (fp.bridged_interactions or {}).keys()
            ),
            "water_bridge_in_default_pairwise_run": False,
            "water_bridge_requirement": (
                "separate bridged-interaction analysis with explicit water route"
            ),
        },
    }

    if num_total == 0:
        has_h = None
        num_atoms = lig_mol.GetNumAtoms()
        if rdkit_mol is not None:
            has_h = any(a.GetAtomicNum() == 1 for a in rdkit_mol.GetAtoms())
            num_atoms = rdkit_mol.GetNumAtoms()
        logger.warning(
            f"ProLIF returned ZERO interactions for {ligand_path.name}. "
            f"This is unusual for a co-crystal ligand. Possible causes: "
            f"(1) ligand coordinates are 2D/flat, (2) wrong ligand_resname selection, "
            f"(3) missing hydrogens on protein/ligand, (4) ligand not in pocket."
        )
        result["diagnostics"] = {
            "warning": "zero_interactions_detected",
            "checks": {
                "ligand_is_3d": coord_validation["is_3d"],
                "ligand_has_hydrogens": has_h,
                "ligand_num_atoms": num_atoms,
                "protein_residues_near_ligand": nearby_residue_count,
                "vicinity_cutoff_used": vicinity_cutoff,
                "ligand_centroid": [round(float(c), 2) for c in centroid],
            },
            "suggestions": [
                "Verify ligand coordinates are finite, non-degenerate, and in the receptor coordinate frame",
                "Check that protein PDB has explicit hydrogens (use prepared PDB from protein-prep)",
                "Verify ligand_resname matches actual residue name in structure",
                "Try PLIP as alternative interaction profiler for cross-validation",
            ],
        }

    return result


def _json_safe(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        try:
            return _json_safe(value.item())
        except (TypeError, ValueError):
            pass
    return str(value)


def _rdkit_atom_record(molecule, atom_index, local_index=None) -> dict:
    try:
        index = int(atom_index)
        atom = molecule.GetAtomWithIdx(index)
    except (TypeError, ValueError, RuntimeError, IndexError):
        return {"index": _json_safe(atom_index)}
    record = {
        "index": index,
        "residue_index": _json_safe(local_index),
        "element": atom.GetSymbol(),
        "atomic_number": atom.GetAtomicNum(),
    }
    info = atom.GetPDBResidueInfo()
    if info is not None:
        name = info.GetName().strip()
        chain = info.GetChainId().strip()
        resname = info.GetResidueName().strip()
        resnum = info.GetResidueNumber()
        icode = info.GetInsertionCode().strip()
        if name:
            record["name"] = name
            record["canonical_id"] = f"{chain}:{resname}{resnum}{icode}:{name}"
    return record


def _prolif_roles(interaction_type: str) -> tuple[str | None, str | None]:
    return {
        "HBDonor": ("donor", "acceptor"),
        "HBAcceptor": ("acceptor", "donor"),
        "Cationic": ("cation", "anion"),
        "Anionic": ("anion", "cation"),
        "PiCation": ("aromatic", "cation"),
        "CationPi": ("cation", "aromatic"),
        "HalogenBond": ("halogen_donor", "halogen_acceptor"),
        "XBDonor": ("halogen_donor", "halogen_acceptor"),
        "XBAcceptor": ("halogen_acceptor", "halogen_donor"),
        # ProLIF naming is ligand-centric: MetalDonor means the ligand is the
        # metal center; MetalAcceptor means the protein residue is the metal.
        "MetalDonor": ("metal_center", "metal_coordinator"),
        "MetalAcceptor": ("metal_coordinator", "metal_center"),
    }.get(interaction_type, (None, None))


def _serialize_ifp(fp, ligand_molecule, protein_molecule) -> list[dict]:
    """Serialize ProLIF 2.x ``fp.ifp`` metadata instead of discarding it."""
    details = []
    for frame_index, frame in sorted((getattr(fp, "ifp", {}) or {}).items()):
        for residue_pair, interaction_map in frame.items():
            if not isinstance(residue_pair, tuple) or len(residue_pair) != 2:
                continue
            ligand_residue, protein_residue = residue_pair
            for interaction_type, metadata_items in interaction_map.items():
                for metadata in metadata_items or ():
                    indices = metadata.get("indices") or {}
                    parent = metadata.get("parent_indices") or indices
                    ligand_local = list(indices.get("ligand") or [])
                    protein_local = list(indices.get("protein") or [])
                    ligand_parent = list(parent.get("ligand") or ligand_local)
                    protein_parent = list(parent.get("protein") or protein_local)
                    ligand_role, protein_role = _prolif_roles(str(interaction_type))
                    geometry = {
                        key: _json_safe(value)
                        for key, value in metadata.items()
                        if key not in {"indices", "parent_indices"}
                        and isinstance(value, (int, float))
                    }
                    details.append({
                        "type": str(interaction_type),
                        "residue": str(protein_residue),
                        "ligand_residue": str(ligand_residue),
                        "ligand_role": ligand_role,
                        "protein_role": protein_role,
                        "ligand_atoms": [
                            _rdkit_atom_record(ligand_molecule, index, local)
                            for index, local in zip(ligand_parent, ligand_local)
                        ],
                        "protein_atoms": [
                            _rdkit_atom_record(protein_molecule, index, local)
                            for index, local in zip(protein_parent, protein_local)
                        ],
                        "geometry": geometry,
                        "raw_metadata": _json_safe(metadata),
                        "frame_index": _json_safe(frame_index),
                        "atom_index_namespace": "prolif_parent_rdkit_index",
                    })
    return details

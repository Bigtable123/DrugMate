"""PLIP utilities: protein-ligand interaction profiling directly from PDB.

PLIP works natively with PDB complex files — no MDAnalysis→RDKit bond-guessing
needed. This makes it far more robust than ProLIF for prepared protein structures
where PDBFixer-added hydrogens can cause sanitization failures.

Install: pip install plip
"""
from __future__ import annotations

import logging
import tempfile
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _plip_version() -> str | None:
    try:
        return version("plip")
    except PackageNotFoundError:
        return None


def run_plip_interactions(
    protein_file: Path,
    ligand_file: Path,
    ligand_resname: str = "LIG",
    chain: str = "A",
) -> dict:
    """Analyze protein-ligand interactions using PLIP.

    PLIP requires a single PDB file containing both protein and ligand.
    If they're in separate files, we merge them before analysis.
    """
    from plip.structure.preparation import PDBComplex

    protein_path = Path(protein_file)
    ligand_path = Path(ligand_file)

    # PLIP needs a combined PDB with both protein and ligand
    complex_pdb = _make_complex_pdb(protein_path, ligand_path, ligand_resname)

    mol = PDBComplex()
    mol.load_pdb(str(complex_pdb))
    mol.analyze()

    # Find the binding site for our ligand
    binding_site = None
    for bsid, site in mol.interaction_sets.items():
        # bsid format: "RESNAME:CHAIN:RESNUMBER"
        parts = bsid.split(":")
        if len(parts) >= 1 and parts[0].strip() == ligand_resname:
            binding_site = site
            break

    if binding_site is None:
        if len(mol.interaction_sets) == 1:
            binding_site = next(iter(mol.interaction_sets.values()))
        else:
            return {
                "success": False,
                "status": "tool_failed",
                "tool": "PLIP",
                "tool_version": _plip_version(),
                "error": f"PLIP found no unambiguous binding site for ligand '{ligand_resname}'",
                "num_binding_sites_found": len(mol.interaction_sets),
                "binding_site_ids": list(mol.interaction_sets.keys()),
            }

    if binding_site is None:
        return {
            "success": False,
            "status": "tool_failed",
            "tool": "PLIP",
            "tool_version": _plip_version(),
            "error": f"PLIP found no binding site for ligand '{ligand_resname}'",
            "num_binding_sites_found": len(mol.interaction_sets),
            "binding_site_ids": list(mol.interaction_sets.keys()),
        }

    # Extract all interactions
    interactions = _extract_interactions(binding_site)

    # Build residue summary
    summary = {}
    for interaction in interactions:
        res_id = interaction["residue"]
        itype = interaction["type"]
        if res_id not in summary:
            summary[res_id] = []
        if itype not in summary[res_id]:
            summary[res_id].append(itype)

    residues = sorted(summary.keys())
    all_types = sorted({t for types in summary.values() for t in types})

    bitvector = []
    for res in residues:
        row = [1 if t in summary[res] else 0 for t in all_types]
        bitvector.append(row)

    # Clean up temp file if created
    if complex_pdb != protein_path:
        try:
            complex_pdb.unlink()
        except Exception:
            pass

    return {
        "success": True,
        "status": "detected" if interactions else "not_detected",
        "tool": "PLIP",
        "tool_version": _plip_version(),
        "residues": residues,
        "interaction_types": all_types,
        "bitvector": bitvector,
        "summary": summary,
        "interactions_detail": interactions,
        "num_interacting_residues": len(residues),
        "num_interactions_total": len(interactions),
        "chemical_representation": {
            "protein_format": protein_path.suffix.lower().lstrip("."),
            "ligand_input_format": ligand_path.suffix.lower().lstrip("."),
            "plip_complex_format": "pdb",
            "ligand_conversion": (
                "internal_sdf_to_pdb_coordinates_elements_conect"
                if ligand_path.suffix.lower() == ".sdf" else "pdb_records_reused"
            ),
            "bond_order_and_formal_charge_policy": (
                "not_preserved_explicitly; OpenBabel/PLIP may re-perceive chemistry"
                if ligand_path.suffix.lower() == ".sdf" else "as represented/perceived from PDB"
            ),
        },
        "quality_flags": (
            ["ligand_sdf_to_pdb_chemical_reperception"]
            if ligand_path.suffix.lower() == ".sdf" else []
        ),
    }


def _make_complex_pdb(
    protein_file: Path, ligand_file: Path, ligand_resname: str
) -> Path:
    """Merge protein PDB and ligand PDB/SDF into a single complex PDB for PLIP."""
    ligand_path = Path(ligand_file)

    # If ligand is SDF, convert to PDB block first
    if ligand_path.suffix.lower() == ".sdf":
        ligand_lines = _sdf_to_pdb_lines(ligand_path, ligand_resname)
    else:
        ligand_lines = []
        with open(ligand_path) as f:
            for line in f:
                if line.startswith(("HETATM", "ATOM", "CONECT")):
                    ligand_lines.append(line)

    if not ligand_lines:
        raise ValueError(f"No atom lines extracted from ligand file: {ligand_path}")

    # Read protein lines (skip END/TER at the very end)
    protein_lines = []
    with open(protein_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TER")):
                protein_lines.append(line)
    ligand_lines = _renumber_ligand_lines(ligand_lines, _max_atom_serial(protein_lines) + 1)

    # Write combined complex
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    tmp.writelines(protein_lines)
    tmp.write("TER\n")
    tmp.writelines(ligand_lines)
    tmp.write("END\n")
    tmp.close()

    return Path(tmp.name)


def _sdf_to_pdb_lines(sdf_file: Path, resname: str = "LIG") -> list[str]:
    """Convert SDF ligand to HETATM PDB lines preserving 3D coordinates."""
    # Do not import RDKit here.  In the drugmate SBDD env, PLIP/OpenBabel and
    # RDKit can abort the interpreter when loaded in the same worker process.
    # A V2000 SDF has all data PLIP needs: atom coordinates/elements + bonds.
    try:
        text = Path(sdf_file).read_text(errors="ignore")
        first_record = text.split("$$$$", 1)[0].splitlines()
        if len(first_record) < 4:
            return []

        counts = first_record[3]
        if "V3000" in counts:
            atom_records, bond_records = _parse_v3000_sdf(first_record)
        else:
            num_atoms = int(counts[0:3])
            num_bonds = int(counts[3:6])
            atom_records = [
                {
                    "x": float(line[0:10]),
                    "y": float(line[10:20]),
                    "z": float(line[20:30]),
                    "element": line[31:34].strip() or "C",
                }
                for line in first_record[4:4 + num_atoms]
            ]
            bond_records = [
                (int(line[0:3]), int(line[3:6]))
                for line in first_record[4 + num_atoms:4 + num_atoms + num_bonds]
            ]
    except (OSError, ValueError):
        return []

    lines = []
    for idx, atom in enumerate(atom_records, start=1):
        x = atom["x"]
        y = atom["y"]
        z = atom["z"]
        element = atom["element"]
        atom_name = _format_pdb_atom_name(element, idx)
        # HETATM serial, atom name, altLoc, resname, chain, resseq, insertion,
        # coords, occupancy, temp factor, element.
        lines.append(
            f"HETATM{idx:5d} {atom_name:<4} {resname[:3]:>3} L   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element[:2]:>2}\n"
        )

    num_atoms = len(atom_records)
    for a1, a2 in bond_records:
        if 1 <= a1 <= num_atoms and 1 <= a2 <= num_atoms:
            lines.append(f"CONECT{a1:5d}{a2:5d}\n")

    return lines


def _parse_v3000_sdf(record_lines: list[str]) -> tuple[list[dict], list[tuple[int, int]]]:
    atoms = []
    bonds = []
    in_atoms = False
    in_bonds = False
    for line in record_lines:
        text = line.strip()
        if text == "M  V30 BEGIN ATOM":
            in_atoms = True
            continue
        if text == "M  V30 END ATOM":
            in_atoms = False
            continue
        if text == "M  V30 BEGIN BOND":
            in_bonds = True
            continue
        if text == "M  V30 END BOND":
            in_bonds = False
            continue
        if in_atoms and text.startswith("M  V30 "):
            parts = text.split()
            if len(parts) >= 7:
                atoms.append({
                    "element": parts[3],
                    "x": float(parts[4]),
                    "y": float(parts[5]),
                    "z": float(parts[6]),
                })
        elif in_bonds and text.startswith("M  V30 "):
            parts = text.split()
            if len(parts) >= 6:
                bonds.append((int(parts[4]), int(parts[5])))
    return atoms, bonds


def _format_pdb_atom_name(element: str, atom_index: int) -> str:
    element = (element or "C").strip()
    if len(element) == 1:
        return f"{element}{atom_index % 1000:>3d}"[:4]
    return f"{element[:2]}{atom_index % 100:>2d}"[:4]


def _max_atom_serial(lines: list[str]) -> int:
    max_serial = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "TER")):
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except ValueError:
                continue
    return max_serial


def _renumber_ligand_lines(lines: list[str], first_serial: int) -> list[str]:
    """Renumber ligand atoms and CONECT records to avoid protein serial clashes."""
    serial_map = {}
    next_serial = first_serial
    renumbered_atoms = []

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            old_serial = int(line[6:11])
        except ValueError:
            old_serial = len(serial_map) + 1
        serial_map[old_serial] = next_serial
        renumbered_atoms.append(f"{line[:6]}{next_serial:5d}{line[11:]}")
        next_serial += 1

    renumbered_conect = []
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        old_serials = []
        for pos in range(6, len(line), 5):
            chunk = line[pos:pos + 5].strip()
            if not chunk:
                continue
            try:
                old_serials.append(int(chunk))
            except ValueError:
                pass
        mapped = [serial_map[s] for s in old_serials if s in serial_map]
        if len(mapped) >= 2:
            renumbered_conect.append("CONECT" + "".join(f"{s:5d}" for s in mapped) + "\n")

    return renumbered_atoms + renumbered_conect


def _extract_interactions(site) -> list[dict]:
    """Extract all interaction details from a PLIP binding site."""
    interactions = []

    # Hydrogen bonds
    for hb in site.hbonds_pdon + site.hbonds_ldon:
        protein_atom = hb.d if hb.protisdon else hb.a
        ligand_atom = hb.a if hb.protisdon else hb.d
        protein_index = hb.d_orig_idx if hb.protisdon else hb.a_orig_idx
        ligand_index = hb.a_orig_idx if hb.protisdon else hb.d_orig_idx
        interactions.append({
            "type": "HBond",
            "residue": f"{hb.restype}{hb.resnr}",
            "chain": hb.reschain,
            "distance": round(hb.distance_ad, 2),
            "distance_ad": round(hb.distance_ad, 2),
            "distance_ah": round(hb.distance_ah, 2),
            "angle": round(hb.angle, 1),
            "donor_type": "protein" if hb.protisdon else "ligand",
            "ligand_role": "acceptor" if hb.protisdon else "donor",
            "protein_role": "donor" if hb.protisdon else "acceptor",
            "protein_atoms": [_atom_record(protein_atom, protein_index)],
            "ligand_atoms": [_atom_record(ligand_atom, ligand_index)],
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Hydrophobic contacts
    for hp in site.hydrophobic_contacts:
        interactions.append({
            "type": "Hydrophobic",
            "residue": f"{hp.restype}{hp.resnr}",
            "chain": hp.reschain,
            "distance": round(hp.distance, 2),
            "protein_atoms": [_atom_record(hp.bsatom, hp.bsatom_orig_idx)],
            "ligand_atoms": [_atom_record(hp.ligatom, hp.ligatom_orig_idx)],
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Pi-stacking
    for ps in site.pistacking:
        interactions.append({
            "type": "PiStacking",
            "residue": f"{ps.restype}{ps.resnr}",
            "chain": ps.reschain,
            "distance": round(ps.distance, 2),
            "angle": round(ps.angle, 1),
            "stacking_type": ps.type,
            "protein_role": "aromatic",
            "ligand_role": "aromatic",
            "protein_atoms": _atom_records(getattr(ps.proteinring, "atoms", [])),
            "ligand_atoms": _atom_records(
                getattr(ps.ligandring, "orig_atoms", None)
                or getattr(ps.ligandring, "atoms", [])
            ),
            "offset": round(ps.offset, 2),
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Pi-cation
    for pc in site.pication_paro + site.pication_laro:
        if pc.protcharged:
            protein_atoms = _atom_records(getattr(pc.charge, "atoms", []))
            ligand_atoms = _atom_records(
                getattr(pc.ring, "orig_atoms", None) or getattr(pc.ring, "atoms", [])
            )
            protein_role, ligand_role = "cation", "aromatic"
        else:
            protein_atoms = _atom_records(getattr(pc.ring, "atoms", []))
            ligand_atoms = _atom_records(
                getattr(pc.charge, "orig_atoms", None) or getattr(pc.charge, "atoms", [])
            )
            protein_role, ligand_role = "aromatic", "cation"
        interactions.append({
            "type": "PiCation",
            "residue": f"{pc.restype}{pc.resnr}",
            "chain": pc.reschain,
            "distance": round(pc.distance, 2),
            "offset": round(pc.offset, 2),
            "protein_role": protein_role,
            "ligand_role": ligand_role,
            "protein_atoms": protein_atoms,
            "ligand_atoms": ligand_atoms,
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Salt bridges
    for sb in site.saltbridge_pneg + site.saltbridge_lneg:
        if sb.protispos:
            protein_group, ligand_group = sb.positive, sb.negative
            protein_role, ligand_role = "cation", "anion"
        else:
            protein_group, ligand_group = sb.negative, sb.positive
            protein_role, ligand_role = "anion", "cation"
        interactions.append({
            "type": "SaltBridge",
            "residue": f"{sb.restype}{sb.resnr}",
            "chain": sb.reschain,
            "distance": round(sb.distance, 2),
            "protein_charge": "positive" if sb.protispos else "negative",
            "protein_role": protein_role,
            "ligand_role": ligand_role,
            "protein_atoms": _atom_records(getattr(protein_group, "atoms", [])),
            "ligand_atoms": _atom_records(
                getattr(ligand_group, "orig_atoms", None)
                or getattr(ligand_group, "atoms", [])
            ),
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Halogen bonds
    for hal in site.halogen_bonds:
        interactions.append({
            "type": "HalogenBond",
            "residue": f"{hal.restype}{hal.resnr}",
            "chain": hal.reschain,
            "distance": round(hal.distance, 2),
            "donor_angle": round(hal.don_angle, 1),
            "acceptor_angle": round(hal.acc_angle, 1),
            "protein_role": "halogen_acceptor",
            "ligand_role": "halogen_donor",
            "protein_atoms": [_atom_record(hal.acc.o, hal.acc_orig_idx)],
            "ligand_atoms": [_atom_record(hal.don.x, hal.don_orig_idx)],
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Water bridges
    for wb in site.water_bridges:
        protein_atom = wb.d if wb.protisdon else wb.a
        ligand_atom = wb.a if wb.protisdon else wb.d
        protein_index = wb.d_orig_idx if wb.protisdon else wb.a_orig_idx
        ligand_index = wb.a_orig_idx if wb.protisdon else wb.d_orig_idx
        interactions.append({
            "type": "WaterBridge",
            "residue": f"{wb.restype}{wb.resnr}",
            "chain": wb.reschain,
            "distance_aw": round(wb.distance_aw, 2),
            "distance_dw": round(wb.distance_dw, 2),
            "donor_angle": round(wb.d_angle, 1),
            "water_angle": round(wb.w_angle, 1),
            "donor_type": "protein" if wb.protisdon else "ligand",
            "protein_role": "donor" if wb.protisdon else "acceptor",
            "ligand_role": "acceptor" if wb.protisdon else "donor",
            "protein_atoms": [_atom_record(protein_atom, protein_index)],
            "ligand_atoms": [_atom_record(ligand_atom, ligand_index)],
            "mediator": {
                "type": "water",
                "atom": _atom_record(wb.water, wb.water_orig_idx),
            },
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    # Metal coordination
    for mc in site.metal_complexes:
        location = str(getattr(mc, "location", "")).lower()
        if location == "water":
            continue
        if location == "ligand":
            # The detected metal belongs to the receptor and the target atom
            # belongs to the ligand.
            residue = f"{mc.restype_l}{mc.resnr_l}"
            residue_chain = mc.reschain_l
            protein_atom = getattr(mc.metal, "m", mc.metal)
            protein_index = mc.metal_orig_idx
            ligand_atom = getattr(mc.target, "atom", mc.target)
            ligand_index = mc.target_orig_idx
            protein_role, ligand_role = "metal_center", "metal_coordinator"
        else:
            # Normal PLIP case: ligand metal coordinated by a protein target.
            residue = f"{mc.restype}{mc.resnr}"
            residue_chain = mc.reschain
            protein_atom = getattr(mc.target, "atom", mc.target)
            protein_index = mc.target_orig_idx
            ligand_atom = getattr(mc.metal, "m", mc.metal)
            ligand_index = mc.metal_orig_idx
            protein_role, ligand_role = "metal_coordinator", "metal_center"
        interactions.append({
            "type": "MetalComplex",
            "residue": residue,
            "chain": residue_chain,
            "distance": round(mc.distance, 2),
            "metal": mc.metal_type,
            "protein_role": protein_role,
            "ligand_role": ligand_role,
            "protein_atoms": [_atom_record(protein_atom, protein_index)],
            "ligand_atoms": [_atom_record(ligand_atom, ligand_index)],
            "mediator": {"type": "metal", "element": mc.metal_type},
            "metal_location": location or "unknown",
            "coordination_number": mc.coordination_num,
            "coordination_geometry": mc.geometry,
            "rms": round(float(mc.rms), 3) if mc.rms == mc.rms else None,
            "atom_index_namespace": "plip_complex_orig_idx",
        })

    return interactions


def _atom_record(atom, original_index=None) -> dict:
    """Extract stable atom identity from a PLIP/OpenBabel atom defensively."""
    record = {}
    if original_index is not None:
        try:
            record["index"] = int(original_index)
        except (TypeError, ValueError):
            record["index"] = str(original_index)
    if atom is None:
        return record
    try:
        record.setdefault("index", int(atom.idx))
    except (AttributeError, TypeError, ValueError):
        pass
    obatom = getattr(atom, "OBAtom", None)
    if obatom is not None:
        try:
            residue = obatom.GetResidue()
            name = residue.GetAtomID(obatom).strip() if residue is not None else ""
            if name:
                record["name"] = name
        except Exception:
            pass
        try:
            record["atomic_number"] = int(obatom.GetAtomicNum())
        except Exception:
            pass
    atom_type = getattr(atom, "type", None)
    if atom_type:
        record["atom_type"] = str(atom_type)
    return record


def _atom_records(atoms) -> list[dict]:
    records = []
    for atom in atoms or []:
        original = getattr(atom, "orig_idx", None)
        records.append(_atom_record(atom, original))
    return records

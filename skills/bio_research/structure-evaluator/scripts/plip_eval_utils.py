"""PLIP evaluation utilities: per-pose interaction analysis.

Runs PLIP on each docked pose individually to extract protein-ligand
interactions for comparison against the interaction blueprint.

PLIP/OpenBabel and RDKit can crash the interpreter when loaded together in some
Drugmate SBDD environments.  This module therefore runs PLIP in a subprocess and
parses docked SDF V2000 records as text inside that PLIP-only worker.

Install: pip install plip
"""
from __future__ import annotations

import json
import logging
import os
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)


def run_plip_per_pose(receptor_file: Path, docked_sdf: Path) -> list[dict]:
    """Run PLIP on each pose in the docked SDF in a PLIP-only subprocess."""
    timeout = int(os.environ.get("SBDD_PLIP_EVAL_TIMEOUT_SECONDS", "600"))
    output_fd, output_name = tempfile.mkstemp(
        prefix="structure_eval_plip_",
        suffix=".json",
    )
    os.close(output_fd)
    output_path = Path(output_name)

    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--run",
        str(receptor_file),
        str(docked_sdf),
        str(output_path),
    ]

    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(f"PLIP per-pose worker timed out after {timeout}s") from exc

    try:
        payload = json.loads(output_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = None
    finally:
        output_path.unlink(missing_ok=True)

    if completed.returncode != 0:
        if payload and payload.get("error_type") in {"ImportError", "ModuleNotFoundError"}:
            raise ImportError(payload.get("error", "plip not installed"))
        if payload and payload.get("error"):
            raise RuntimeError(payload["error"])
        raise RuntimeError(_format_worker_failure(completed.returncode, completed.stderr))

    if not payload or not payload.get("ok"):
        error = (payload or {}).get("error", "PLIP per-pose worker returned no result")
        raise RuntimeError(error)

    return payload.get("results", [])


def _run_plip_per_pose_in_worker(receptor_file: Path, docked_sdf: Path) -> list[dict]:
    """Worker implementation. Import PLIP here, never in the parent process."""
    from plip.structure.preparation import PDBComplex

    receptor_lines = _read_protein_lines(receptor_file)
    sdf_records = _read_sdf_records(docked_sdf)

    results = []
    for i, record in enumerate(sdf_records):
        mol_name = _sdf_record_name(record, i)

        try:
            ligand_lines = _sdf_record_to_pdb_lines(record, "LIG")
            if not ligand_lines:
                raise ValueError("Failed to parse SDF V2000 atom/bond block")
            complex_pdb = _write_complex(receptor_lines, ligand_lines)

            try:
                plip_mol = PDBComplex()
                plip_mol.load_pdb(str(complex_pdb))
                plip_mol.analyze()
                interactions = _extract_all_interactions(plip_mol)
            finally:
                complex_pdb.unlink(missing_ok=True)

            results.append({
                "pose_index": i,
                "mol_name": mol_name,
                "success": True,
                "interactions": interactions,
                "num_interactions": len(interactions),
            })
        except Exception as e:
            results.append({
                "pose_index": i,
                "mol_name": mol_name,
                "success": False,
                "error": str(e),
            })

    return results


def _format_worker_failure(returncode: int, stderr: str) -> str:
    if returncode < 0:
        signal_number = -returncode
        try:
            signal_name = signal.Signals(signal_number).name
        except ValueError:
            signal_name = f"signal {signal_number}"
        reason = f"PLIP per-pose worker terminated by {signal_name} ({signal_number})"
    else:
        reason = f"PLIP per-pose worker exited with code {returncode}"

    stderr_lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    if stderr_lines:
        reason += ": " + " | ".join(stderr_lines[-8:])
    return reason


def _read_protein_lines(pdb_file: Path) -> list[str]:
    lines = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TER")):
                lines.append(line)
    return lines


def _read_sdf_records(sdf_file: Path) -> list[list[str]]:
    records = []
    current = []
    with open(sdf_file, errors="ignore") as f:
        for line in f:
            if line.rstrip("\n\r") == "$$$$":
                if current:
                    records.append(current)
                current = []
            else:
                current.append(line.rstrip("\n\r"))
    if current and any(line.strip() for line in current):
        records.append(current)
    return records


def _sdf_record_name(record: list[str], pose_index: int) -> str:
    if record:
        name = record[0].strip()
        if name:
            return name
    return f"pose_{pose_index}"


def _sdf_record_to_pdb_lines(record: list[str], resname: str) -> list[str]:
    """Convert one SDF V2000 record to PDB ligand lines without RDKit."""
    if len(record) < 4:
        return []

    counts = record[3]
    if "V3000" in counts:
        raise ValueError("SDF V3000 is not supported by the PLIP text parser")

    try:
        num_atoms = int(counts[0:3])
        num_bonds = int(counts[3:6])
    except ValueError:
        return []

    atom_records = record[4:4 + num_atoms]
    bond_records = record[4 + num_atoms:4 + num_atoms + num_bonds]
    if len(atom_records) != num_atoms:
        return []

    lines = []
    for idx, atom_line in enumerate(atom_records, start=1):
        try:
            x = float(atom_line[0:10])
            y = float(atom_line[10:20])
            z = float(atom_line[20:30])
            element = atom_line[31:34].strip() or "C"
        except ValueError:
            return []

        atom_name = _format_pdb_atom_name(element, idx)
        lines.append(
            f"HETATM{idx:5d} {atom_name:<4} {resname[:3]:>3} L   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element[:2]:>2}\n"
        )

    for bond_line in bond_records:
        try:
            a1 = int(bond_line[0:3])
            a2 = int(bond_line[3:6])
        except ValueError:
            continue
        if 1 <= a1 <= num_atoms and 1 <= a2 <= num_atoms:
            lines.append(f"CONECT{a1:5d}{a2:5d}\n")

    return lines


def _format_pdb_atom_name(element: str, atom_index: int) -> str:
    element = (element or "C").strip()
    if len(element) == 1:
        return f"{element}{atom_index % 1000:>3d}"[:4]
    return f"{element[:2]}{atom_index % 100:>2d}"[:4]


def _write_complex(protein_lines: list[str], ligand_lines: list[str]) -> Path:
    ligand_lines = _renumber_ligand_lines(ligand_lines, _max_atom_serial(protein_lines) + 1)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    tmp.writelines(protein_lines)
    tmp.write("TER\n")
    tmp.writelines(ligand_lines)
    tmp.write("END\n")
    tmp.close()
    return Path(tmp.name)


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
    """Renumber ligand atom serials and CONECT records after receptor atoms."""
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


def _extract_all_interactions(plip_mol) -> list[dict]:
    """Serialize all interaction types for the candidate written as LIG:L:1."""
    ligand_site_id = "LIG:L:1"
    if ligand_site_id not in plip_mol.interaction_sets:
        raise ValueError(
            f"PLIP found no binding site for candidate {ligand_site_id}; "
            f"found {sorted(plip_mol.interaction_sets)}"
        )
    interactions = []
    for _bsid, site in plip_mol.interaction_sets.items():
        if _bsid != ligand_site_id:
            continue
        for hb in site.hbonds_pdon:
            interactions.append(_hbond_record(hb, protein_is_donor=True))
        for hb in site.hbonds_ldon:
            interactions.append(_hbond_record(hb, protein_is_donor=False))
        for hp in site.hydrophobic_contacts:
            distance = _round_or_none(getattr(hp, "distance", None))
            protein_atom = _atom_record(
                getattr(hp, "bsatom", None), getattr(hp, "bsatom_orig_idx", None),
            )
            interactions.append({
                "type": "Hydrophobic",
                "residue": f"{hp.restype}{hp.resnr}",
                "chain": hp.reschain,
                "icode": protein_atom.get("icode", ""),
                "distance": distance,
                "protein_atoms": [protein_atom],
                "ligand_atoms": [_atom_record(
                    getattr(hp, "ligatom", None), getattr(hp, "ligatom_orig_idx", None),
                )],
                "geometry": {
                    "distance": distance,
                    "distance_metric": "atom_atom",
                    "distance_unit": "angstrom",
                },
                "mediator": "direct",
            })
        for ps in site.pistacking:
            distance = _round_or_none(getattr(ps, "distance", None))
            angle = _round_or_none(getattr(ps, "angle", None))
            protein_atoms = _atom_records(
                getattr(getattr(ps, "proteinring", None), "atoms", []),
            )
            interactions.append({
                "type": "PiStacking",
                "residue": f"{ps.restype}{ps.resnr}",
                "chain": ps.reschain,
                "icode": _icode_from_atoms(protein_atoms),
                "protein_role": "aromatic",
                "ligand_role": "aromatic",
                "protein_atoms": protein_atoms,
                "ligand_atoms": _atom_records(
                    getattr(getattr(ps, "ligandring", None), "orig_atoms", None)
                    or getattr(getattr(ps, "ligandring", None), "atoms", []),
                ),
                "distance": distance,
                "angle": angle,
                "geometry": {
                    "distance": distance,
                    "angle": angle,
                    "distance_metric": "ring_centroid",
                    "distance_unit": "angstrom",
                    "angle_unit": "degree",
                },
                "mediator": "direct",
            })
        for pc in site.pication_paro + site.pication_laro:
            if pc.protcharged:
                protein_group, ligand_group = pc.charge, pc.ring
                protein_role, ligand_role = "cation", "aromatic"
            else:
                protein_group, ligand_group = pc.ring, pc.charge
                protein_role, ligand_role = "aromatic", "cation"
            protein_atoms = _atom_records(
                getattr(protein_group, "atoms", []),
            )
            distance = _round_or_none(getattr(pc, "distance", None))
            interactions.append({
                "type": "PiCation",
                "residue": f"{pc.restype}{pc.resnr}",
                "chain": pc.reschain,
                "icode": _icode_from_atoms(protein_atoms),
                "protein_role": protein_role,
                "ligand_role": ligand_role,
                "protein_atoms": protein_atoms,
                "ligand_atoms": _atom_records(
                    getattr(ligand_group, "orig_atoms", None)
                    or getattr(ligand_group, "atoms", []),
                ),
                "distance": distance,
                "geometry": {
                    "distance": distance,
                    "distance_metric": "charge_ring_center",
                    "distance_unit": "angstrom",
                },
                "mediator": "direct",
            })
        for sb in site.saltbridge_pneg:
            interactions.append(_salt_bridge_record(
                sb, protein_role="anion", ligand_role="cation",
                protein_group=sb.negative, ligand_group=sb.positive,
            ))
        for sb in site.saltbridge_lneg:
            interactions.append(_salt_bridge_record(
                sb, protein_role="cation", ligand_role="anion",
                protein_group=sb.positive, ligand_group=sb.negative,
            ))
        for wb in site.water_bridges:
            protein_atom = wb.d if wb.protisdon else wb.a
            ligand_atom = wb.a if wb.protisdon else wb.d
            protein_index = wb.d_orig_idx if wb.protisdon else wb.a_orig_idx
            ligand_index = wb.a_orig_idx if wb.protisdon else wb.d_orig_idx
            protein_atom_record = _atom_record(protein_atom, protein_index)
            distance = _round_or_none(getattr(wb, "distance_aw", None))
            interactions.append({
                "type": "WaterBridge",
                "residue": f"{wb.restype}{wb.resnr}",
                "chain": wb.reschain,
                "icode": protein_atom_record.get("icode", ""),
                "protein_role": "donor" if wb.protisdon else "acceptor",
                "ligand_role": "acceptor" if wb.protisdon else "donor",
                "protein_atoms": [protein_atom_record],
                "ligand_atoms": [_atom_record(ligand_atom, ligand_index)],
                "distance": distance,
                "geometry": {
                    "distance": distance,
                    "donor_angle": _round_or_none(getattr(wb, "d_angle", None)),
                    "water_angle": _round_or_none(getattr(wb, "w_angle", None)),
                    "distance_metric": "acceptor_water",
                    "distance_unit": "angstrom",
                    "angle_unit": "degree",
                },
                "mediator": "water",
            })
        for hal in getattr(site, "halogen_bonds", []):
            protein_atom = _atom_record(
                getattr(getattr(hal, "acc", None), "o", None),
                getattr(hal, "acc_orig_idx", None),
            )
            distance = _round_or_none(getattr(hal, "distance", None))
            interactions.append({
                "type": "HalogenBond",
                "residue": f"{hal.restype}{hal.resnr}",
                "chain": hal.reschain,
                "icode": protein_atom.get("icode", ""),
                "protein_role": "halogen_acceptor",
                "ligand_role": "halogen_donor",
                "protein_atoms": [protein_atom],
                "ligand_atoms": [_atom_record(
                    getattr(getattr(hal, "don", None), "x", None),
                    getattr(hal, "don_orig_idx", None),
                )],
                "distance": distance,
                "geometry": {
                    "distance": distance,
                    "donor_angle": _round_or_none(getattr(hal, "don_angle", None)),
                    "acceptor_angle": _round_or_none(getattr(hal, "acc_angle", None)),
                    "distance_metric": "atom_atom",
                    "distance_unit": "angstrom",
                    "angle_unit": "degree",
                },
                "mediator": "direct",
            })
        for mc in getattr(site, "metal_complexes", []):
            location = str(getattr(mc, "location", "")).lower()
            if location == "water":
                continue
            if location == "ligand":
                residue = f"{mc.restype_l}{mc.resnr_l}"
                residue_chain = mc.reschain_l
                protein_atom = _atom_record(
                    getattr(getattr(mc, "metal", None), "m", getattr(mc, "metal", None)),
                    getattr(mc, "metal_orig_idx", None),
                )
                ligand_atom = _atom_record(
                    getattr(getattr(mc, "target", None), "atom", getattr(mc, "target", None)),
                    getattr(mc, "target_orig_idx", None),
                )
                protein_role, ligand_role = "metal_center", "metal_coordinator"
            else:
                residue = f"{mc.restype}{mc.resnr}"
                residue_chain = mc.reschain
                protein_atom = _atom_record(
                    getattr(getattr(mc, "target", None), "atom", getattr(mc, "target", None)),
                    getattr(mc, "target_orig_idx", None),
                )
                ligand_atom = _atom_record(
                    getattr(getattr(mc, "metal", None), "m", getattr(mc, "metal", None)),
                    getattr(mc, "metal_orig_idx", None),
                )
                protein_role, ligand_role = "metal_coordinator", "metal_center"
            distance = _round_or_none(getattr(mc, "distance", None))
            interactions.append({
                "type": "MetalComplex",
                "residue": residue,
                "chain": residue_chain,
                "icode": protein_atom.get("icode", ""),
                "protein_role": protein_role,
                "ligand_role": ligand_role,
                "protein_atoms": [protein_atom],
                "ligand_atoms": [ligand_atom],
                "distance": distance,
                "geometry": {
                    "distance": distance,
                    "distance_metric": "atom_atom",
                    "distance_unit": "angstrom",
                },
                "mediator": "metal",
                "metal_location": location or "unknown",
            })
    return interactions


def _hbond_record(hbond, *, protein_is_donor: bool) -> dict:
    heavy_distance = _round_or_none(getattr(hbond, "distance_ad", None))
    hydrogen_distance = _round_or_none(getattr(hbond, "distance_ah", None))
    # PLIP 3.x exposes the D-H-A angle as ``angle`` (not ``don_angle``).
    angle = _round_or_none(getattr(hbond, "angle", None))
    protein_atom = hbond.d if protein_is_donor else hbond.a
    ligand_atom = hbond.a if protein_is_donor else hbond.d
    protein_index = hbond.d_orig_idx if protein_is_donor else hbond.a_orig_idx
    ligand_index = hbond.a_orig_idx if protein_is_donor else hbond.d_orig_idx
    protein_atom_record = _atom_record(protein_atom, protein_index)
    geometry = {
        "distance": heavy_distance,
        "heavy_atom_distance": heavy_distance,
        "acceptor_hydrogen_distance": hydrogen_distance,
        "distance_metric": "heavy_atom",
        "distance_unit": "angstrom",
        "angle_unit": "degree",
    }
    if angle is not None:
        geometry["angle"] = angle
        geometry["donor_angle"] = angle
    return {
        "type": "HBond",
        "residue": f"{hbond.restype}{hbond.resnr}",
        "chain": hbond.reschain,
        "icode": protein_atom_record.get("icode", ""),
        "distance": heavy_distance,
        "angle": angle,
        "geometry": geometry,
        "protein_role": "donor" if protein_is_donor else "acceptor",
        "ligand_role": "acceptor" if protein_is_donor else "donor",
        "protein_atoms": [protein_atom_record],
        "ligand_atoms": [_atom_record(ligand_atom, ligand_index)],
        "mediator": "direct",
    }


def _salt_bridge_record(
    sb, *, protein_role: str, ligand_role: str, protein_group, ligand_group,
) -> dict:
    protein_atoms = _atom_records(getattr(protein_group, "atoms", []))
    distance = _round_or_none(getattr(sb, "distance", None))
    return {
        "type": "SaltBridge",
        "residue": f"{sb.restype}{sb.resnr}",
        "chain": sb.reschain,
        "icode": _icode_from_atoms(protein_atoms),
        "protein_role": protein_role,
        "ligand_role": ligand_role,
        "protein_atoms": protein_atoms,
        "ligand_atoms": _atom_records(
            getattr(ligand_group, "orig_atoms", None)
            or getattr(ligand_group, "atoms", []),
        ),
        "distance": distance,
        "geometry": {
            "distance": distance,
            "distance_metric": "charge_center",
            "distance_unit": "angstrom",
        },
        "mediator": "direct",
    }


def _atom_record(atom, original_index=None) -> dict:
    """Serialize PLIP atom identity without crossing atom-index namespaces.

    PLIP's ``orig_idx`` refers to its corrected, hydrogen-filtered complex.  It
    must not be used as an ordinal or serial into the input receptor PDB.  The
    OpenBabel atom and residue attached to the interaction are the authoritative
    source for atom and residue identity.
    """
    record = {}
    try:
        index = int(original_index)
        record["index"] = index
    except (TypeError, ValueError):
        pass
    if atom is None:
        return record
    try:
        record.setdefault("index", int(atom.idx))
    except (AttributeError, TypeError, ValueError):
        pass
    obatom = getattr(atom, "OBAtom", None)
    if obatom is not None:
        record.update(_ob_residue_atom_metadata(obatom))
        try:
            record["atomic_number"] = int(obatom.GetAtomicNum())
        except Exception:
            pass
    atom_type = getattr(atom, "type", None)
    if atom_type:
        record["atom_type"] = str(atom_type)
    return record


def _ob_residue_atom_metadata(obatom) -> dict:
    """Return independently recoverable identity fields from an OBResidue."""
    try:
        residue = obatom.GetResidue()
    except Exception:
        residue = None
    if residue is None:
        return {}

    metadata = {}
    try:
        name = _clean_ob_text(residue.GetAtomID(obatom))
        if name:
            metadata["name"] = name
    except Exception:
        pass
    try:
        resname = _clean_ob_text(residue.GetName())
        if resname:
            metadata["resname"] = resname
    except Exception:
        pass
    try:
        chain = _clean_ob_text(residue.GetChain())
        metadata["chain"] = chain
    except Exception:
        pass
    try:
        metadata["resnum"] = int(residue.GetNum())
    except (AttributeError, TypeError, ValueError):
        pass
    try:
        metadata["icode"] = _clean_ob_text(residue.GetInsertionCode())
    except Exception:
        pass
    return metadata


def _clean_ob_text(value) -> str:
    return str(value or "").replace("\x00", "").strip()


def _atom_records(atoms) -> list[dict]:
    return [
        _atom_record(atom, getattr(atom, "orig_idx", None))
        for atom in (atoms or [])
    ]


def _icode_from_atoms(atoms: list[dict]) -> str:
    return next((str(atom.get("icode")) for atom in atoms if atom.get("icode")), "")


def _round_or_none(value, digits: int = 2):
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return None


def _main() -> int:
    if len(sys.argv) == 5 and sys.argv[1] == "--run":
        output_path = Path(sys.argv[4])
        try:
            results = _run_plip_per_pose_in_worker(Path(sys.argv[2]), Path(sys.argv[3]))
            payload = {"ok": True, "results": results}
            exit_code = 0
        except Exception as exc:
            payload = {
                "ok": False,
                "error_type": type(exc).__name__,
                "error": f"{type(exc).__name__}: {exc}",
            }
            exit_code = 1

        output_path.write_text(
            json.dumps(payload, ensure_ascii=False, default=str),
            encoding="utf-8",
        )
        return exit_code

    print(
        "Usage: plip_eval_utils.py --run RECEPTOR_PDB DOCKED_SDF OUTPUT_JSON",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    sys.exit(_main())

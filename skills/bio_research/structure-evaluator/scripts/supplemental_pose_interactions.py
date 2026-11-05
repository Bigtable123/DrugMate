"""Candidate-pose Arpeggio/LUNA/ODDT using the existing isolated workers."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from tool_runtime import run_tool

ROOT = Path(__file__).resolve().parents[2]
WORKER = ROOT / "pocket-analysis/scripts/native_worker.py"


def candidate_complex_mmcif(receptor_file, ligand, output_file):
    """Build an explicitly derived complex without pretending it is source data."""
    import gemmi

    structure = gemmi.read_structure(str(receptor_file))
    if len(structure) != 1:
        raise ValueError("candidate Arpeggio analysis requires one receptor coordinate model")
    used_chains = {chain.name for chain in structure[0]}
    chain_id = next((name for name in "ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789" if name not in used_chains), None)
    if chain_id is None:
        raise ValueError("no free chain identifier for the candidate ligand")
    used_resnames = {residue.name for chain in structure[0] for residue in chain}
    resname = next((name for name in ("LIG", "DRG", "MOL") if name not in used_resnames), None)
    if resname is None:
        raise ValueError("no unambiguous ligand residue name is available")
    chain = gemmi.Chain(chain_id)
    residue = gemmi.Residue()
    residue.name, residue.seqid, residue.het_flag = resname, gemmi.SeqId(1, " "), "H"
    conformer = ligand.GetConformer()
    names = []
    for index, atom in enumerate(ligand.GetAtoms()):
        target = gemmi.Atom()
        target.name = f"{atom.GetSymbol()}{index + 1}"
        names.append(target.name)
        target.element = gemmi.Element(atom.GetSymbol())
        target.charge = atom.GetFormalCharge()
        xyz = conformer.GetAtomPosition(index)
        target.pos = gemmi.Position(xyz.x, xyz.y, xyz.z)
        target.occ = 1.0
        residue.add_atom(target)
    chain.add_residue(residue)
    structure[0].add_chain(chain)
    structure.setup_entities()
    structure.assign_label_seq_id()
    document = structure.make_mmcif_document()
    block = document.sole_block()
    types = {res.name: "non-polymer" if res.het_flag == "H" else "L-peptide linking"
             for ch in structure[0] for res in ch}
    block.set_mmcif_category("_chem_comp.", {"id": list(types), "type": list(types.values()),
                                                     "name": ["WATER" if name in {"HOH", "WAT"} else name for name in types]})
    block.set_mmcif_category("_chem_comp_atom.", {
        "comp_id": [resname] * len(names), "atom_id": names,
        "type_symbol": [atom.GetSymbol() for atom in ligand.GetAtoms()],
        "charge": [str(atom.GetFormalCharge()) for atom in ligand.GetAtoms()],
        "pdbx_aromatic_flag": ["Y" if atom.GetIsAromatic() else "N" for atom in ligand.GetAtoms()],
    })
    orders = {"SINGLE": "SING", "DOUBLE": "DOUB", "TRIPLE": "TRIP", "AROMATIC": "AROM"}
    bonds = list(ligand.GetBonds())
    if bonds:
        block.set_mmcif_category("_chem_comp_bond.", {
            "comp_id": [resname] * len(bonds),
            "atom_id_1": [names[b.GetBeginAtomIdx()] for b in bonds],
            "atom_id_2": [names[b.GetEndAtomIdx()] for b in bonds],
            "value_order": [orders.get(str(b.GetBondType()), "SING") for b in bonds],
            "pdbx_aromatic_flag": ["Y" if b.GetIsAromatic() else "N" for b in bonds],
        })
    document.write_file(str(output_file))
    return resname, chain_id


def attach_supplemental_interactions(evaluation, receptor_file, docked_sdf, output_dir, configs):
    from rdkit import Chem

    enabled = {key: configs.get(key) or {} for key in ("arpeggio", "luna", "oddt")
               if (configs.get(key) or {}).get("enabled", False)}
    if not enabled:
        return []
    destination = Path(output_dir) / "supplemental_interactions"
    destination.mkdir(parents=True, exist_ok=True)
    receptor_copy = destination / "receptor.pdb"
    shutil.copyfile(receptor_file, receptor_copy)
    receptor_hash = hashlib.sha256(receptor_copy.read_bytes()).hexdigest()
    molecules = list(Chem.SDMolSupplier(str(docked_sdf), removeHs=False))
    errors = []
    for entry in evaluation:
        index = int(entry["pose_index"])
        pose_dir = destination / f"pose_{index:06d}"
        pose_dir.mkdir(parents=True, exist_ok=True)
        ligand_file = pose_dir / "ligand.sdf"
        molecule = molecules[index]
        if molecule is not None:
            with Chem.SDWriter(str(ligand_file)) as writer:
                writer.write(molecule)
        for tool, config in enabled.items():
            response = pose_dir / f"{tool}.json"
            response.unlink(missing_ok=True)
            params = dict(config.get("params") or {})
            timeout = int(params.pop("timeout_seconds", 300))
            try:
                if molecule is None:
                    raise ValueError("unreadable docked ligand")
                payload = {"protein_file": str(receptor_copy), "ligand_file": str(ligand_file),
                           "ligand_resname": "LIG", **params}
                if tool == "arpeggio":
                    complex_file = pose_dir / "candidate_complex.cif"
                    resname, chain = candidate_complex_mmcif(receptor_copy, molecule, complex_file)
                    payload.update(complex_mmcif_file=str(complex_file), ligand_resname=resname,
                                   chain=chain, tool_timeout_seconds=timeout)
                runtime = {"kind": "command", "command": "python",
                           "environment": "drugmate" if tool == "arpeggio" else tool,
                           **(config.get("runtime") or {})}
                execution = run_tool(
                    runtime, [WORKER, "--step", tool, "--input", json.dumps(payload), "--output", response],
                    command="python", cwd=pose_dir, timeout=timeout + 30, check=False,
                    env={"PYTHONDONTWRITEBYTECODE": "1", "TMPDIR": str(pose_dir),
                         "XDG_CACHE_HOME": str(pose_dir / "cache")})
                if not response.is_file():
                    raise RuntimeError(str(execution.get("stderr") or "native worker produced no response")[-2000:])
                wrapped = json.loads(response.read_text())
                if not wrapped.get("ok"):
                    raise RuntimeError(wrapped.get("error") or "native interaction worker failed")
                result = wrapped["result"]
                if result.get("success") is not True:
                    raise RuntimeError(result.get("error") or "interaction profiler returned failure")
                result["input_provenance"] = {
                    "receptor_sha256": receptor_hash,
                    "ligand_sdf_sha256": hashlib.sha256(ligand_file.read_bytes()).hexdigest(),
                    "pose_index": index, "coordinate_source": "current docked pose",
                    "complex_kind": "derived_candidate_complex" if tool == "arpeggio" else "separate_receptor_and_ligand",
                }
                result["execution"] = {key: execution.get(key) for key in ("command", "elapsed_seconds", "returncode")}
                if tool == "arpeggio":
                    qc = result.get("qc_counts") or {}
                    count = qc.get("atom_atom_records", 0)
                    result["metrics"] = {
                        "clash_contacts": qc.get("clash"), "vdw_clash_contacts": qc.get("vdw_clash"),
                        "proximal_contact_fraction": qc.get("proximal", 0) / count if count else None,
                    } if result.get("quality") == "valid" else {}
                    if not result["metrics"]:
                        result["metric_unavailability_reason"] = "Arpeggio chemistry quality is degraded"
                    entry["arpeggio"] = result
                entry.setdefault("supplemental_interactions", {})[tool] = result
            except Exception as exc:
                entry.setdefault("supplemental_interactions", {})[tool] = {
                    "success": False, "status": "failed", "error": str(exc)}
                errors.append({"tool": tool, "pose_index": index, "reason": str(exc)})
    return errors

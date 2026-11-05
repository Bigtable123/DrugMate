"""PoseCheck utilities: strain energy and steric clash analysis.

PoseCheck computes ligand strain energy and counts steric clashes
between docked poses and the receptor protein.

Install: pip install posecheck
"""
from __future__ import annotations

from pathlib import Path


def check_posecheck_available() -> bool:
    try:
        from posecheck import PoseCheck
        return True
    except ImportError:
        return False


def _count_heavy_atom_clashes(protein_mol, ligand_mol, tolerance: float) -> int:
    """Count protein-ligand clashes using PoseCheck's rule, excluding hydrogens."""
    from rdkit import Chem

    pt = Chem.GetPeriodicTable()
    protein_conf = protein_mol.GetConformer()
    ligand_conf = ligand_mol.GetConformer()
    protein_heavy = [
        atom.GetIdx() for atom in protein_mol.GetAtoms()
        if atom.GetAtomicNum() != 1
    ]
    ligand_heavy = [
        atom.GetIdx() for atom in ligand_mol.GetAtoms()
        if atom.GetAtomicNum() != 1
    ]

    clashes = 0
    for lig_idx in ligand_heavy:
        lig_atom = ligand_mol.GetAtomWithIdx(lig_idx)
        lig_pos = ligand_conf.GetAtomPosition(lig_idx)
        lig_vdw = pt.GetRvdw(lig_atom.GetAtomicNum())

        for prot_idx in protein_heavy:
            prot_atom = protein_mol.GetAtomWithIdx(prot_idx)
            cutoff = lig_vdw + pt.GetRvdw(prot_atom.GetAtomicNum()) - tolerance
            if cutoff <= 0:
                continue

            prot_pos = protein_conf.GetAtomPosition(prot_idx)
            dx = lig_pos.x - prot_pos.x
            dy = lig_pos.y - prot_pos.y
            dz = lig_pos.z - prot_pos.z
            if (dx * dx) + (dy * dy) + (dz * dz) < cutoff * cutoff:
                clashes += 1

    return clashes


def run_posecheck(
    docked_sdf: Path, receptor_file: Path, *, num_conformers: int = 10,
) -> list[dict]:
    """Run PoseCheck on docked poses to get strain energy and clash counts.

    ``steric_clash_count`` is the decision metric and counts only heavy-atom
    protein-ligand clashes. ``steric_clash_count_all_atom`` preserves
    PoseCheck's native all-atom clash count for diagnostics.

    Reuse the final prepared receptor, including its explicit hydrogens. Protein
    clashes need only elements and coordinates, so skip bond inference and
    sanitization rather than rebuilding the receptor with Hydride/ProLIF.

    The installed PoseCheck class hides its strain function's sampling budget.
    Call that same function with an explicit budget, keeping its local
    relaxation, force field, and starting-pose relaxation unchanged.
    """
    if (
        isinstance(num_conformers, bool)
        or not isinstance(num_conformers, int)
        or num_conformers < 1
    ):
        raise ValueError("PoseCheck num_conformers must be a positive integer")
    from rdkit import Chem

    protein = Chem.MolFromPDBFile(
        str(receptor_file), sanitize=False, removeHs=False, proximityBonding=False,
    )
    if protein is None or not protein.GetNumAtoms() or not protein.GetNumConformers():
        raise ValueError(f"Failed to load prepared receptor PDB: {receptor_file}")
    if not any(atom.GetAtomicNum() == 1 for atom in protein.GetAtoms()):
        raise ValueError(
            "PoseCheck requires a prepared receptor PDB with explicit hydrogens; "
            "use the final docking prepared_receptor."
        )

    from posecheck import PoseCheck
    from posecheck.utils.strain import calculate_strain_energy

    pc = PoseCheck()
    pc.protein = protein

    suppl = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)

    results = []
    for i, mol in enumerate(suppl):
        if mol is None:
            results.append({
                "pose_index": i,
                "mol_name": f"pose_{i}",
                "success": False,
                "error": "Failed to parse molecule",
            })
            continue

        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"pose_{i}"

        try:
            pc.load_ligands_from_mols([mol])
            clashes = pc.calculate_clashes()
            strain_energy = calculate_strain_energy(
                pc.ligands[0], num_confs=num_conformers,
            )

            all_atom_clash_count = clashes[0] if clashes else 0
            heavy_atom_clash_count = _count_heavy_atom_clashes(
                pc.protein,
                pc.ligands[0],
                tolerance=pc.clash_tolerance,
            )
            strain_energy_raw = (
                float(strain_energy) if strain_energy is not None else None
            )

            results.append({
                "pose_index": i,
                "mol_name": mol_name,
                "success": True,
                "steric_clash_count": heavy_atom_clash_count,
                "steric_clash_count_all_atom": all_atom_clash_count,
                "steric_clash_count_mode": "heavy_atom_only",
                "num_conformers": num_conformers,
                "strain_energy_kcal": strain_energy_raw,
                "strain_energy_kcal_display": (
                    round(strain_energy_raw, 2)
                    if strain_energy_raw is not None else None
                ),
            })
        except Exception as e:
            results.append({
                "pose_index": i,
                "mol_name": mol_name,
                "success": False,
                "error": str(e),
            })

    return results

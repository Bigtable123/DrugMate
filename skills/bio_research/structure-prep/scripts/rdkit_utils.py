"""RDKit utilities: PDB to SDF conversion (preserves 3D crystal coordinates)."""

from __future__ import annotations

import logging
import urllib.request
from pathlib import Path
from typing import Any

from rcsb_http import RCSBRequestError, request_bytes

logger = logging.getLogger(__name__)

RCSB_LIGAND_VIEW = "https://files.rcsb.org/ligands/view"


def check_rdkit_available() -> bool:
    try:
        from rdkit import Chem  # noqa: F401
        import numpy  # noqa: F401
        return True
    except ImportError:
        return False


def _is_3d(mol) -> bool:
    """Check that the molecule has non-degenerate 3D coordinates."""
    import numpy as np

    try:
        conf = mol.GetConformer()
    except ValueError:
        return False
    if mol.GetNumAtoms() == 0:
        return False
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    return bool(np.std(coords[:, 2]) > 0.01)


def pdb_to_sdf(
    input_pdb: Path,
    output_sdf: Path,
    sanitize: bool = True,
    comp_id: str | None = None,
    template_cache_dir: Path | None = None,
) -> dict[str, Any]:
    """Convert ligand PDB to SDF, preserving the original 3D crystal pose.

    If ``comp_id`` is supplied, first fetch the wwPDB/RCSB CCD ideal SDF and use
    it as a chemical template. The resulting molecule keeps the PDB crystal
    coordinates but takes bond orders, aromaticity, and formal charges from CCD.
    """
    from rdkit import Chem

    logger.info(f"Converting {input_pdb} to SDF format (preserving 3D pose)")

    if comp_id:
        try:
            template = _load_ccd_template(comp_id, template_cache_dir or output_sdf.parent / "ccd_templates")
            mol = _pdb_to_template_assigned_mol(input_pdb, template["mol"])
            _write_sdf(mol, output_sdf)
            props = _mol_properties(mol)
            props.update({
                "bond_order_source": "ccd_ideal_sdf",
                "template_comp_id": comp_id.upper(),
                "template_file": str(template["path"]),
                "template_graph_verified": True,
                "crystal_coordinates_preserved": True,
            })
            logger.info(
                "  SDF saved: %s (3D=%s, template=%s)",
                output_sdf,
                props["is_3d"],
                comp_id.upper(),
            )
            logger.info(f"  SMILES: {props['smiles']}")
            return {"success": True, "output_sdf": str(output_sdf), "properties": props}
        except Exception as exc:
            error = f"CCD template assignment failed for {comp_id}: {exc}"
            logger.error("%s", error)
            return {"success": False, "error": error}

    mol = Chem.MolFromPDBFile(str(input_pdb), sanitize=sanitize, removeHs=False)
    if mol is None:
        mol = Chem.MolFromPDBFile(str(input_pdb), sanitize=False, removeHs=False)
        if mol is None:
            return {"success": False, "error": "RDKit failed to parse PDB file"}

    is_3d = _is_3d(mol)
    if not is_3d:
        logger.warning("Ligand coordinates appear to be 2D — interaction analysis may fail")

    _write_sdf(mol, output_sdf)
    props = _mol_properties(mol)
    props["bond_order_source"] = "rdkit_pdb_perception"

    logger.info(f"  SDF saved: {output_sdf} (3D={is_3d})")
    logger.info(f"  SMILES: {props['smiles']}")
    return {"success": True, "output_sdf": str(output_sdf), "properties": props}


def generate_2d_sdf(
    input_pdb: Path,
    output_sdf: Path,
    sanitize: bool = True,
    comp_id: str | None = None,
    template_cache_dir: Path | None = None,
) -> dict[str, Any]:
    """Generate a 2D depiction SDF for visualization (NOT for interaction analysis)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = None
    if comp_id:
        try:
            mol = Chem.Mol(_load_ccd_template(comp_id, template_cache_dir or output_sdf.parent / "ccd_templates")["mol"])
        except Exception as exc:
            return {"success": False, "error": f"CCD 2D template load failed for {comp_id}: {exc}"}

    if mol is None:
        mol = Chem.MolFromPDBFile(str(input_pdb), sanitize=sanitize, removeHs=True)
    if mol is None:
        mol = Chem.MolFromPDBFile(str(input_pdb), sanitize=False, removeHs=True)
        if mol is None:
            return {"success": False, "error": "RDKit failed to parse PDB file"}

    AllChem.Compute2DCoords(mol)

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()

    return {"success": True, "output_sdf": str(output_sdf)}


def _load_ccd_template(comp_id: str, cache_dir: Path) -> dict[str, Any]:
    from rdkit import Chem

    ccd_id = comp_id.upper()
    cache_dir.mkdir(parents=True, exist_ok=True)
    template_path = cache_dir / f"{ccd_id}_ideal.sdf"

    if not template_path.exists() or template_path.stat().st_size == 0:
        url = f"{RCSB_LIGAND_VIEW}/{ccd_id}_ideal.sdf"
        logger.info("Fetching CCD template %s -> %s", url, template_path)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sbdd-structure-prep/1.0"})
            template_path.write_bytes(request_bytes(req))
        except RCSBRequestError as exc:
            raise RuntimeError(f"failed to fetch CCD ideal SDF: {exc}") from exc

    suppl = Chem.SDMolSupplier(str(template_path), removeHs=True, sanitize=True)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise RuntimeError(f"failed to parse CCD ideal SDF: {template_path}")
    return {"mol": mol, "path": template_path}


def _pdb_to_template_assigned_mol(input_pdb: Path, template_mol):
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

    pdb_mol = Chem.MolFromPDBFile(str(input_pdb), sanitize=False, removeHs=False)
    if pdb_mol is None:
        raise RuntimeError("RDKit failed to parse ligand PDB")

    # removeHs on the PDB reader is ineffective when sanitization is disabled.
    # Assigning a heavy-atom CCD template with deposited H still present can
    # double-count H and produce impossible valences (for example ERW in 6G9H).
    pdb_mol = Chem.RemoveHs(pdb_mol, sanitize=False)
    if pdb_mol.GetNumAtoms() != template_mol.GetNumAtoms():
        raise ValueError("Ligand heavy-atom count does not match CCD template")
    coordinates = pdb_mol.GetConformer().GetPositions()
    if not np.isfinite(coordinates).all():
        raise ValueError("Ligand crystal coordinates are not finite")

    mol = AssignBondOrdersFromTemplate(template_mol, pdb_mol)
    Chem.SanitizeMol(mol)
    # Check the full graph, not merely a substructure match. This comparison
    # covers bonds/charges/H; stereochemistry requires separate validation.
    if Chem.MolToSmiles(mol, isomericSmiles=False) != Chem.MolToSmiles(template_mol, isomericSmiles=False):
        raise ValueError("Assigned ligand chemical graph does not match CCD template")
    if ([atom.GetAtomicNum() for atom in mol.GetAtoms()] !=
            [atom.GetAtomicNum() for atom in pdb_mol.GetAtoms()] or
            not np.allclose(mol.GetConformer().GetPositions(), coordinates, rtol=0, atol=1e-6)):
        raise ValueError("CCD assignment changed crystal atom order or coordinates")
    if not _is_3d(mol):
        logger.warning("Template-assigned ligand coordinates appear to be 2D")
    mol.SetProp("_Name", input_pdb.stem)
    return mol


def _write_sdf(mol, output_sdf: Path) -> None:
    from rdkit import Chem

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()


def _mol_properties(mol) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    try:
        smiles = Chem.MolToSmiles(mol)
    except Exception:
        smiles = ""
    return {
        "smiles": smiles,
        "num_atoms": mol.GetNumAtoms(),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "molecular_weight": round(Descriptors.ExactMolWt(mol), 2),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "aromatic_bonds": sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic()),
        "is_3d": _is_3d(mol),
    }

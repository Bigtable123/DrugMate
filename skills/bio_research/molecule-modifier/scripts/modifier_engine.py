"""Molecule Modifier Engine — apply rule-based modifications to molecules."""
from __future__ import annotations

import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parents[1]
_SHARED_DIR = _SKILL_DIR.parent / "shared"
sys.path.insert(0, str(_SKILL_DIR / "scripts"))
sys.path.insert(0, str(_SHARED_DIR))

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors
from reaction_rules import REACTION_RULES, match_keywords, matched_any


def modify_molecule(smiles: str, fix_instruction: str, max_variants: int = 3) -> list[dict]:
    """Apply rule-based modifications to a molecule.

    Returns list of dicts: [{smiles, modification_type, rule_applied, instruction_recognized}, ...]
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    rule_keys = match_keywords(fix_instruction)
    recognized = matched_any(fix_instruction)
    variants = []
    seen = {Chem.MolToSmiles(mol)}

    for rule_key in rule_keys:
        reactions = REACTION_RULES.get(rule_key, [])
        for pattern_sma, product_sma in reactions:
            products = _apply_replacement(mol, pattern_sma, product_sma)
            for product_smi in products:
                if product_smi in seen:
                    continue
                seen.add(product_smi)
                if _validate(smiles, product_smi):
                    variants.append({
                        "smiles": product_smi,
                        "modification_type": rule_key,
                        "rule_applied": f"{pattern_sma} -> {product_sma}",
                        "instruction_recognized": recognized,
                    })
                    if len(variants) >= max_variants:
                        return variants
    return variants


def _apply_replacement(mol, pattern_sma: str, product_sma: str) -> list[str]:
    """Apply a single SMARTS replacement, returning up to 3 unique products."""
    try:
        reaction = AllChem.ReactionFromSmarts(f"{pattern_sma}>>{product_sma}")
    except Exception:
        return []
    if reaction is None:
        return []
    try:
        reaction.Initialize()
    except Exception:
        return []

    results = []
    try:
        product_sets = reaction.RunReactants((mol,))
    except Exception:
        return []

    for product_set in product_sets:
        if not product_set:
            continue
        new_mol = product_set[0]
        try:
            Chem.SanitizeMol(new_mol)
            smi = Chem.MolToSmiles(new_mol, canonical=True)
        except Exception:
            continue
        if smi and smi not in results:
            results.append(smi)
        if len(results) >= 3:
            break
    return results


def _validate(original_smi: str, modified_smi: str) -> bool:
    """Validate a modification is chemically reasonable."""
    orig = Chem.MolFromSmiles(original_smi)
    mod = Chem.MolFromSmiles(modified_smi)
    if orig is None or mod is None:
        return False
    # Atom count change < 30%
    orig_atoms = orig.GetNumHeavyAtoms()
    mod_atoms = mod.GetNumHeavyAtoms()
    if abs(mod_atoms - orig_atoms) > orig_atoms * 0.3:
        return False
    # Tanimoto similarity > 0.4
    fp1 = AllChem.GetMorganFingerprintAsBitVect(orig, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mod, 2, nBits=2048)
    if DataStructs.TanimotoSimilarity(fp1, fp2) < 0.4:
        return False
    return True


def generate_modified_library(
    instructions: dict,
    smiles_map: dict[str, str],
    max_variants_per_mol: int = 3,
) -> list[dict]:
    """Generate modified compounds from Agent-authored modification instructions.

    Args:
        instructions: target_molecules, per_molecule, and suggestions from the Agent
        smiles_map: compound_id -> SMILES mapping

    Returns:
        List of compound dicts ready for conformer generation.
    """
    per_molecule = instructions.get("per_molecule", {})
    target_molecules = instructions.get("target_molecules", [])
    compounds = []
    idx = 0

    for mol_name in target_molecules:
        smiles = smiles_map.get(mol_name)
        if not smiles:
            continue
        mol_instructions = per_molecule.get(mol_name, {})
        fix = mol_instructions.get("fix_instruction", "")
        if not fix:
            # Use global suggestions as fallback
            suggestions = instructions.get("suggestions", [])
            fix = suggestions[0] if suggestions else "H-bond donor"

        variants = modify_molecule(smiles, fix, max_variants=max_variants_per_mol)
        for v in variants:
            idx += 1
            compounds.append({
                "compound_id": f"MOD_{idx:04d}",
                "smiles": v["smiles"],
                "parent_id": mol_name,
                "parent_smiles": smiles,
                "modification_type": v["modification_type"],
                "rule_applied": v["rule_applied"],
                "source": "molecule_modifier",
            })

    return compounds

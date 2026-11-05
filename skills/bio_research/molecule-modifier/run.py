#!/usr/bin/env python3
"""Molecule Modifier Skill — rule-based molecular modification.

Applies RDKit SMARTS-based transformations to molecules based on
coordinator fix_instructions, generating modified variants for the next
docking round.

Usage:
    python run.py --input '{"library_manifest_file": "step2/library_manifest.csv", "modification_instructions_file": "step6/modification_instructions.json"}' --output-dir ./modified/
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))

from modifier_engine import generate_modified_library
from io_utils import parse_skill_args, load_smiles_map


def run(input_data: dict, output_dir: Path) -> dict:
    """Generate modified molecules from coordinator instructions."""
    manifest_file = input_data.get("library_manifest_file")
    instructions_file = input_data.get("modification_instructions_file")
    instructions_inline = input_data.get("modification_instructions")
    max_variants = int(input_data.get("max_variants_per_mol", 3))
    seed = int(input_data.get("seed", 42))

    # Load SMILES map
    smiles_map = load_smiles_map(manifest_file)
    if not smiles_map:
        return {"success": False, "error": f"No SMILES loaded from {manifest_file}"}

    # Load instructions
    if instructions_inline:
        instructions = instructions_inline
    elif instructions_file:
        path = Path(instructions_file)
        if not path.exists():
            return {"success": False, "error": f"Instructions file not found: {path}"}
        instructions = json.loads(path.read_text())
    else:
        return {"success": False, "error": "Must provide modification_instructions or modification_instructions_file"}

    if not instructions.get("target_molecules"):
        return {"success": False, "error": "No target_molecules in instructions"}

    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate modifications
    compounds = generate_modified_library(
        instructions=instructions,
        smiles_map=smiles_map,
        max_variants_per_mol=max_variants,
    )

    if not compounds:
        return {
            "success": False,
            "error": "No valid modifications could be generated",
            "target_molecules": instructions.get("target_molecules"),
        }

    # Surface molecules whose fix_instruction matched no known rule (the engine
    # falls back to a mild H-bond-donor edit). The Agent should treat these as
    # not genuinely directed, rather than assuming the intended chemistry ran.
    unrecognized = sorted({
        c.get("parent_id", "")
        for c in compounds
        if c.get("instruction_recognized") is False
    } - {""})

    # Generate 3D conformers and write SDF
    output_sdf = output_dir / "modified_library.sdf"
    manifest_csv = output_dir / "modified_manifest.csv"

    sdf_result = _write_sdf_with_conformers(compounds, output_sdf, seed)
    _write_manifest(compounds, manifest_csv)

    # Write report
    report = {
        "success": True,
        "num_target_molecules": len(instructions.get("target_molecules", [])),
        "num_modifications_generated": len(compounds),
        "num_conformers_written": sdf_result["written"],
        "num_conformers_failed": sdf_result["failed"],
        "unrecognized_instruction_parents": unrecognized,
        "files": {
            "modified_library_sdf": str(output_sdf),
            "modified_manifest_csv": str(manifest_csv),
        },
        "modifications": compounds,
    }
    (output_dir / "modification_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False)
    )
    return report


def _write_sdf_with_conformers(compounds: list[dict], output_sdf: Path, seed: int) -> dict:
    """Generate 3D conformers and write to SDF."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    writer = Chem.SDWriter(str(output_sdf))
    written = 0
    failed = 0

    for comp in compounds:
        mol = Chem.MolFromSmiles(comp["smiles"])
        if mol is None:
            failed += 1
            continue
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            failed += 1
            continue
        # Optimize
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        mol = Chem.RemoveHs(mol)
        # Set properties
        mol.SetProp("_Name", comp.get("compound_id", ""))
        mol.SetProp("smiles", comp.get("smiles", ""))
        mol.SetProp("source", "molecule_modifier")
        mol.SetProp("parent_id", comp.get("parent_id", ""))
        mol.SetProp("parent_smiles", comp.get("parent_smiles", ""))
        mol.SetProp("modification_type", comp.get("modification_type", ""))
        mw = f"{Descriptors.ExactMolWt(mol):.2f}"
        logp = f"{Descriptors.MolLogP(mol):.2f}"
        mol.SetProp("MW", mw)
        mol.SetProp("LogP", logp)
        # Persist descriptors back onto the compound dict so the manifest carries
        # them too (parity with the molecule-generator manifest schema).
        comp["MW"] = mw
        comp["LogP"] = logp
        writer.write(mol)
        written += 1

    writer.close()
    return {"written": written, "failed": failed}


def _write_manifest(compounds: list[dict], path: Path) -> None:
    """Write modification manifest CSV."""
    import csv
    fields = ["compound_id", "smiles", "parent_id", "parent_smiles", "modification_type", "rule_applied", "source", "MW", "LogP"]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for comp in compounds:
            writer.writerow({k: comp.get(k, "") for k in fields})


def main():
    input_data, output_dir = parse_skill_args("Molecule Modifier Skill")
    result = run(input_data, output_dir)
    print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

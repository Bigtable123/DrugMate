"""SMILES standardization utilities using RDKit.

Handles salt removal, neutralization, canonicalization, and deduplication.
"""
from __future__ import annotations

from typing import Tuple
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from molecule_failures import MoleculeFailureTracker, is_rdkit_parse_error


def check_standardization_available() -> bool:
    """Check if RDKit standardization modules are available."""
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
        return True
    except ImportError:
        return False


def standardize_smiles_list(
    compounds: list[dict],
    smiles_key: str = "smiles",
    neutralize: bool = False,
    failure_tracker: MoleculeFailureTracker | None = None,
) -> Tuple[list[dict], list[dict]]:
    """Standardize SMILES for a list of compound dicts.

    For each compound:
      1. Parse SMILES to mol
      2. Remove salts (keep largest fragment)
      3. Preserve formal charges by default (optional neutralization)
      4. Canonicalize SMILES

    Returns: (standardized_compounds, failed_compounds)
    Each compound in standardized list gets 'smiles_standardized' field added.
    Failed compounds get 'standardization_error' field.
    """
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    lfc = rdMolStandardize.LargestFragmentChooser()
    uc = rdMolStandardize.Uncharger()

    standardized = []
    failed = []
    tracker = failure_tracker or MoleculeFailureTracker()
    operation = "standardization.smiles"

    for comp in compounds:
        raw_smiles = comp.get(smiles_key, "")
        if not raw_smiles:
            failed.append({**comp, "standardization_error": "empty SMILES"})
            tracker.record(comp, operation, "molecule_parse_failed", "empty SMILES")
            continue

        try:
            mol = Chem.MolFromSmiles(raw_smiles)
            if mol is None:
                failed.append({**comp, "standardization_error": "RDKit parse failed"})
                tracker.record(comp, operation, "molecule_parse_failed", "RDKit parse failed")
                continue

            mol = lfc.choose(mol)
            if neutralize:
                mol = uc.uncharge(mol)
            canonical = Chem.MolToSmiles(mol, canonical=True)

            if not canonical:
                failed.append({**comp, "standardization_error": "canonicalization produced empty SMILES"})
                tracker.record(comp, operation, "molecule_parse_failed", "canonicalization produced empty SMILES")
                continue

            comp_out = {**comp, "smiles_standardized": canonical}
            standardized.append(comp_out)
        except Exception as e:
            failed.append({**comp, "standardization_error": str(e)})
            tracker.record(comp, operation, "standardization_failed", str(e),
                           exception=None if is_rdkit_parse_error(e) else e)

    return standardized, failed


def deduplicate_compounds(
    compounds: list[dict],
    smiles_key: str = "smiles_standardized",
    activity_key: str = "activity_nm",
) -> list[dict]:
    """Deduplicate compounds by standardized SMILES.

    For actives, preserve ChEMBL censoring semantics while aggregating repeated
    measurements. Exact values may contribute to a median; ``<``/``<=`` values
    remain bounds and are never treated as exact observations.
    For decoys (no activity), keeps the first occurrence.
    """
    grouped = {}
    for comp in compounds:
        smi = comp.get(smiles_key, comp.get("smiles", ""))
        if not smi:
            continue
        grouped.setdefault(smi, []).append(comp)

    deduplicated = []
    for entries in grouped.values():
        activity_entries = [
            entry for entry in entries if entry.get(activity_key) is not None
        ]
        if not activity_entries:
            deduplicated.append(entries[0])
            continue
        # Local import keeps the generic standardization helpers importable in
        # contexts that do not use ChEMBL while sharing one censoring contract
        # across pre- and post-standardization deduplication.
        from chembl_utils import aggregate_activity_measurements

        aggregation_entries = activity_entries
        if activity_key != "activity_nm":
            aggregation_entries = [
                {**entry, "activity_nm": entry[activity_key]}
                for entry in activity_entries
            ]
        representative = aggregate_activity_measurements(aggregation_entries)
        if activity_key != "activity_nm":
            representative[activity_key] = representative["activity_nm"]
        representative["standardized_entry_count"] = len(entries)
        deduplicated.append(representative)
    return deduplicated

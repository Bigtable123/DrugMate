"""RDKit-based property-matched decoy generation.

Generates decoys by BRICS fragmentation and recombination, filtered by
property similarity and topological dissimilarity to active compounds.
Used as fallback when DUD-E decoys are not available.
"""
from __future__ import annotations

import random
from itertools import islice
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from molecule_failures import MoleculeFailureTracker, is_rdkit_parse_error


def check_rdkit_available() -> bool:
    """Check if RDKit is available for decoy generation."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, BRICS
        from rdkit import DataStructs
        return True
    except ImportError:
        return False


def _compute_properties(molecule) -> dict | None:
    """Compute properties on the same sanitized object used for fingerprints."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(molecule) if isinstance(molecule, str) else molecule
    if mol is None:
        return None
    return {
        "MW": Descriptors.ExactMolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "TPSA": Descriptors.TPSA(mol),
    }


def _compute_fingerprint(mol):
    """Compute Morgan fingerprint (ECFP4, radius=2, 2048 bits)."""
    try:
        from rdkit.Chem import rdFingerprintGenerator
    except ImportError:
        rdFingerprintGenerator = None
    if rdFingerprintGenerator is None or not hasattr(rdFingerprintGenerator, "GetMorganGenerator"):
        from rdkit.Chem import AllChem
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return gen.GetFingerprint(mol)


def _max_tanimoto_to_actives(candidate_fp, active_fps) -> float:
    """Compute maximum Tanimoto similarity to any active compound."""
    from rdkit import DataStructs
    if not active_fps:
        return 0.0
    sims = DataStructs.BulkTanimotoSimilarity(candidate_fp, active_fps)
    return max(sims) if sims else 0.0


def _property_matches(props: dict, target_props: dict, tolerances: dict) -> bool:
    """Check if properties are within tolerance of target distribution."""
    for key, tol in tolerances.items():
        if key not in props or key not in target_props:
            continue
        target_val = target_props[key]
        if abs(props[key] - target_val) > tol:
            return False
    return True


def generate_property_matched_decoys(
    active_smiles_list: list[str],
    num_decoys: int = 200,
    max_tanimoto: float = 0.35,
    seed: int = 42,
    failure_tracker: MoleculeFailureTracker | None = None,
) -> dict:
    """Build sanitized BRICS decoys and audit each failed molecule/attempt."""
    from rdkit import Chem
    from rdkit.Chem import BRICS

    tracker = failure_tracker or MoleculeFailureTracker()
    first_failure = len(tracker.failures)
    rng = random.Random(seed)
    active_mols, active_props_list, active_fps = [], [], []
    tolerances = {"MW": 125.0, "LogP": 1.5, "HBD": 2, "HBA": 3,
                  "RotBonds": 3, "TPSA": 40.0}
    mean_props, all_fragments, decoys = {}, set(), []
    attempts = 0

    def attempt(record, operation, function):
        operation = "rdkit_decoys." + operation
        try:
            value = function()
        except Exception as exc:
            tracker.record(record, operation, "decoy_processing_failed", str(exc),
                           exception=None if is_rdkit_parse_error(exc) else exc)
            return None
        return value

    def result(error=None):
        return {
            "success": bool(decoys), "error": error, "decoys": decoys,
            "source": "rdkit_property_matched", "count": len(decoys),
            "molecule_failures": tracker.failures[first_failure:],
            "generation_params": {
                "max_tanimoto": max_tanimoto, "seed": seed,
                "num_fragments_used": len(all_fragments), "attempts": attempts,
                "target_num_decoys": num_decoys, "tolerances": tolerances,
                "mean_active_properties": mean_props,
                "reagent_order": "sorted_fragments_local_rng_no_global_scrambling",
            },
        }

    for index, smi in enumerate(active_smiles_list, 1):
        record = {"source_record_id": f"rdkit_decoys.active_{index}",
                  "smiles": smi, "candidate_role": "active"}
        before = len(tracker.failures)
        mol = attempt(record, "parse", lambda: Chem.MolFromSmiles(smi))
        if mol is None:
            if len(tracker.failures) == before:
                tracker.record(record, "rdkit_decoys.parse", "molecule_parse_failed",
                               "RDKit could not parse the active SMILES")
            continue
        props = attempt(record, "properties", lambda: _compute_properties(mol))
        fp = attempt(record, "fingerprint", lambda: _compute_fingerprint(mol))
        if props is None or fp is None:
            continue
        active_mols.append((mol, record))
        active_props_list.append(props)
        active_fps.append(fp)

    if not active_mols:
        return result("No active molecules with valid properties and fingerprints")
    for key in tolerances:
        mean_props[key] = sum(p[key] for p in active_props_list) / len(active_props_list)

    for mol, record in active_mols:
        frags = attempt(record, "fragmentation", lambda: BRICS.BRICSDecompose(mol))
        if frags:
            all_fragments.update(frags)
    if len(all_fragments) < 3:
        return result("Insufficient BRICS fragments for recombination")

    frag_list = sorted(all_fragments)
    fragment_mols = {s: Chem.MolFromSmiles(s) for s in frag_list}
    seen_smiles = {Chem.MolToSmiles(mol) for mol, _ in active_mols}
    while len(decoys) < num_decoys and attempts < num_decoys * 50:
        attempts += 1
        selected = rng.sample(frag_list, rng.randint(2, min(4, len(frag_list))))
        record = {"source_record_id": f"rdkit_decoys.attempt_{attempts}",
                  "candidate_role": "validation_decoy", "source": "rdkit_property_matched"}
        products = attempt(record, "recombination", lambda: list(islice(
            BRICS.BRICSBuild([fragment_mols[f] for f in selected],
                             scrambleReagents=False), 200)))
        if not products:
            continue
        candidate = Chem.Mol(rng.choice(products))
        # Reaction products have not initialized valence/ring caches. Sanitizing
        # this exact object is required before either descriptors or ECFP4.
        def sanitize():
            Chem.SanitizeMol(candidate)
            return candidate
        candidate = attempt(record, "sanitize", sanitize)
        if candidate is None:
            continue
        candidate_smi = attempt(record, "canonicalize", lambda: Chem.MolToSmiles(candidate))
        if not candidate_smi or candidate_smi in seen_smiles:
            continue
        record["smiles"] = candidate_smi
        props = attempt(record, "properties", lambda: _compute_properties(candidate))
        if props is None or not _property_matches(props, mean_props, tolerances):
            continue
        fp = attempt(record, "fingerprint", lambda: _compute_fingerprint(candidate))
        if fp is None:
            continue
        max_tc = attempt(record, "similarity", lambda: _max_tanimoto_to_actives(fp, active_fps))
        if max_tc is None or max_tc > max_tanimoto:
            continue
        seen_smiles.add(candidate_smi)
        decoys.append({"smiles": candidate_smi, "name": f"DECOY_{len(decoys):04d}",
                       "source": "rdkit_property_matched", **props})
    return result(None if decoys else "No decoys passed the property/similarity filters")

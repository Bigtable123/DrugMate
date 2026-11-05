"""CPU RAscore XGB batch worker; run in the registered RAscore environment."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def run(payload: dict) -> dict:
    from RAscore import RAscore_XGB
    from rdkit import Chem

    model = Path(RAscore_XGB.__file__).parent / "models/XGB_chembl_ecfp_counts/model.pkl"
    scorer = RAscore_XGB.RAScorerXGB(str(model))
    scorer.xgb_model.set_params(n_jobs=1)
    results = []
    seen = set()
    for candidate in payload["candidates"]:
        identity = str(candidate["candidate_id"])
        if not identity or identity in seen:
            raise ValueError("RAscore requires unique nonempty candidate_id values")
        seen.add(identity)
        row = {"candidate_id": identity, "smiles": candidate["smiles"]}
        if "model_representation" in candidate:
            row.update(model_representation=candidate["model_representation"],
                       model_smiles=candidate.get("model_smiles"))
        try:
            if candidate.get("representation_error"):
                raise ValueError(candidate["representation_error"])
            if "model_representation" in candidate:
                if candidate["model_representation"] != "stereo_independent_ecfp6_counts_v1":
                    raise ValueError("unsupported RAscore model representation")
                model_smiles = candidate.get("model_smiles")
            else:
                model_smiles = candidate["smiles"]  # Standalone legacy requests.
            mol = Chem.MolFromSmiles(model_smiles) if model_smiles else None
            if mol is None or not mol.GetNumAtoms():
                raise ValueError("unparseable molecule")
            value = float(scorer.predict(model_smiles))
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("RAscore returned a non-finite or out-of-range probability")
            row.update(status="available", score=value)
        except Exception as exc:
            row.update(status="failed", score=None, error=str(exc))
        results.append(row)
    return {
        "success": all(r["status"] == "available" for r in results),
        "backend": "RAscore_XGB_ChEMBL_ECFP6_counts",
        "model_sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
        "result_semantics": "Predicted retrosynthetic accessibility; no actual synthesis route is claimed.",
        "results": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    result = run(json.loads(args.input_json.read_text()))
    args.output_json.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

"""Run and join RAscore by candidate identity before the SYN category fold."""
from __future__ import annotations

import json
import math
from pathlib import Path

from tool_runtime import run_tool

MODEL_REPRESENTATION = "stereo_independent_ecfp6_counts_v1"


def model_candidate(result):
    """Keep chemical identity separate from the installed model's achiral input.

    RAscore_XGB uses Morgan counts without chirality. Only this model's copy
    loses stereo; isotope labels, charges and the original SMILES are retained.
    This also avoids sending newer non-tetrahedral SMILES to its older RDKit.
    """
    from rdkit import Chem

    row = {"candidate_id": result["candidate_id"], "smiles": result["smiles"],
           "model_representation": MODEL_REPRESENTATION, "model_smiles": None}
    try:
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None or not mol.GetNumAtoms():
            raise ValueError("unparseable molecule")
        model_mol = Chem.Mol(mol)
        Chem.RemoveStereochemistry(model_mol)
        row["model_smiles"] = Chem.MolToSmiles(
            model_mol, canonical=True, isomericSmiles=True,
        )
    except Exception as exc:
        row["representation_error"] = str(exc)
    return row


def attach_rascore(input_data, output_dir, results, completed, skipped):
    config = (input_data.get("tools") or {}).get("rascore") or {}
    if not config.get("enabled", False):
        skipped.append({"step": "rascore", "reason": "disabled_by_config"})
        return {"status": "disabled"}
    candidates = [model_candidate(r)
                  for r in results if r.get("valid_structure") and r.get("smiles")
                  and not r.get("structural_reject")]
    destination = Path(output_dir) / "rascore"
    destination.mkdir(parents=True, exist_ok=True)
    request, response = destination / "input.json", destination / "output.json"
    request.write_text(json.dumps({"candidates": candidates}, indent=2))
    response.unlink(missing_ok=True)
    runtime = {"kind": "command", "environment": "rascore", "command": "python",
               **(config.get("runtime") or {})}
    try:
        execution = run_tool(
            runtime, [Path(__file__).with_name("rascore_worker.py"),
                      "--input-json", request, "--output-json", response],
            command="python", cwd=destination,
            timeout=int((config.get("params") or {}).get("timeout_seconds", 600)),
            env={"PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1"},
        )
        report = json.loads(response.read_text())
        rows = report["results"]
        expected = {r["candidate_id"]: r["smiles"] for r in candidates}
        representations = {r["candidate_id"]: r for r in candidates}
        indexed = {r["candidate_id"]: r for r in rows}
        if len(indexed) != len(rows) or set(indexed) != set(expected):
            raise ValueError("RAscore result identities do not match the request")
        for identity, row in indexed.items():
            if row.get("smiles") != expected[identity]:
                raise ValueError("RAscore result SMILES does not match the request")
            for key in ("model_smiles", "model_representation"):
                if row.get(key) != representations[identity][key]:
                    raise ValueError(f"RAscore result {key} does not match the request")
            if row.get("status") == "available":
                value = float(row["score"])
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("invalid RAscore probability")
            elif row.get("score") is not None:
                raise ValueError("failed RAscore row must not contain a numeric score")
        report["execution"] = {key: execution.get(key) for key in
                               ("command", "returncode", "elapsed_seconds")}
        for result in results:
            if result.get("candidate_id") in indexed:
                result["rascore"] = {**indexed[result["candidate_id"]],
                                     "model_sha256": report["model_sha256"],
                                     "backend": report["backend"]}
        if report.get("success"):
            completed.append("rascore")
        else:
            skipped.append({"step": "rascore", "reason": "one or more predictions failed"})
    except Exception as exc:
        report = {"status": "failed", "success": False, "error": str(exc)}
        for result in results:
            if result.get("valid_structure"):
                result["rascore"] = {"status": "failed", "score": None, "error": str(exc)}
        skipped.append({"step": "rascore", "reason": str(exc)})
    return report

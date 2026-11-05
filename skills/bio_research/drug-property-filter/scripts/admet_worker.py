#!/usr/bin/env python
"""Batch ADMET-AI worker executed inside the dedicated ``admet_ai`` env."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from admet_contract import (  # noqa: E402
    ADMET_NATIVE_ENDPOINT_IDS,
    endpoint_value_error,
    finite_number,
    validate_native_predictions,
)

# The installed torch build targets CUDA 13 while this host currently exposes an
# older driver.  ADMET-AI is fast enough on CPU for the post-pose candidate set,
# so force CPU before importing torch/admet_ai.
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


PERCENTILE_SUFFIX = "_drugbank_approved_percentile"


def _json_value(value: Any) -> Any:
    """Convert pandas/numpy values to strict, finite JSON values."""
    if value is None:
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return str(value)


def _endpoint_resources() -> tuple[Path, dict[str, dict], str]:
    import admet_ai
    import pandas as pd

    path = Path(admet_ai.__file__).resolve().parent / "resources" / "data" / "admet.csv"
    table = pd.read_csv(path)
    metadata: dict[str, dict] = {}
    for row in table.to_dict(orient="records"):
        endpoint_id = str(row.get("id") or "")
        if not endpoint_id:
            continue
        metadata[endpoint_id] = {
            "category": _json_value(row.get("category")),
            "name": _json_value(row.get("name")),
            "task_type": _json_value(row.get("task_type")),
            "units": _json_value(row.get("units")),
            "species": _json_value(row.get("species")),
        }
    missing_registered = sorted(set(ADMET_NATIVE_ENDPOINT_IDS) - set(metadata))
    if missing_registered:
        raise ValueError(
            "ADMET-AI provider schema is missing registered endpoints: "
            f"{missing_registered}"
        )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return path, metadata, digest


def _summary(predictions: dict[str, dict]) -> dict:
    def value(endpoint_id: str) -> float | None:
        raw = (predictions.get(endpoint_id) or {}).get("value")
        return float(raw) if isinstance(raw, (int, float)) and math.isfinite(float(raw)) else None

    herg = value("hERG")
    cyp3a4 = value("CYP3A4_Veith")
    ames = value("AMES")
    bioavailability = value("Bioavailability_Ma")
    return {
        "hERG": {
            "probability": herg,
            "risk": (
                "unknown" if herg is None
                else "high" if herg > 0.5
                else "moderate" if herg > 0.2
                else "low"
            ),
        },
        "CYP3A4_Veith": {
            "probability": cyp3a4,
            "inhibitor": (
                "unknown" if cyp3a4 is None else "yes" if cyp3a4 > 0.5 else "no"
            ),
            "gate_band": (
                "unknown" if cyp3a4 is None
                else "reject" if cyp3a4 > 0.9
                else "modify" if cyp3a4 > 0.3
                else "keep"
            ),
        },
        "AMES": {
            "probability": ames,
            "call": "unknown" if ames is None else "positive" if ames > 0.5 else "negative",
            "risk": (
                "unknown" if ames is None
                else "high" if ames > 0.5
                else "moderate" if ames > 0.2
                else "low"
            ),
        },
        "Bioavailability_Ma": {
            "probability": bioavailability,
            "call": (
                "unknown" if bioavailability is None
                else "high" if bioavailability >= 0.7
                else "moderate" if bioavailability >= 0.3
                else "low"
            ),
        },
    }


def predict(payload: dict) -> dict:
    from admet_ai import ADMETModel

    entries = payload.get("candidates")
    if not isinstance(entries, list) or not entries:
        raise ValueError("candidates must be a non-empty list")

    candidate_ids: list[str] = []
    smiles: list[str] = []
    seen: set[str] = set()
    for position, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"candidate {position} must be an object")
        candidate_id = str(entry.get("candidate_id") or "")
        candidate_smiles = str(entry.get("smiles") or "")
        if not candidate_id or not candidate_smiles:
            raise ValueError(f"candidate {position} needs candidate_id and smiles")
        if candidate_id in seen:
            raise ValueError(f"duplicate candidate_id: {candidate_id}")
        seen.add(candidate_id)
        candidate_ids.append(candidate_id)
        smiles.append(candidate_smiles)

    metadata_path, endpoint_metadata, schema_sha256 = _endpoint_resources()
    model = ADMETModel()
    frame = model.predict(smiles)

    # Batch predictions must preserve a one-to-one, ordered association.  The
    # backend silently drops invalid SMILES, so refuse any cardinality/index drift.
    if not hasattr(frame, "iloc") or not hasattr(frame, "index"):
        raise TypeError("ADMET-AI batch prediction did not return a DataFrame")
    if len(frame) != len(smiles):
        raise ValueError(
            f"ADMET-AI cardinality mismatch: requested {len(smiles)}, returned {len(frame)}"
        )
    returned_smiles = [str(item) for item in frame.index.tolist()]
    if returned_smiles != smiles:
        raise ValueError(
            "ADMET-AI output index does not exactly match the ordered input SMILES"
        )

    results = []
    for candidate_id, candidate_smiles, (_, row) in zip(
        candidate_ids, smiles, frame.iterrows(), strict=True
    ):
        raw = {str(key): _json_value(value) for key, value in row.to_dict().items()}
        predictions: dict[str, dict] = {
            endpoint_id: {
                **endpoint_metadata[endpoint_id],
                "value": raw.get(endpoint_id),
            }
            for endpoint_id in ADMET_NATIVE_ENDPOINT_IDS
        }
        unregistered_predictions: dict[str, dict] = {
            endpoint_id: {
                **endpoint_record,
                "value": raw.get(endpoint_id),
            }
            for endpoint_id, endpoint_record in endpoint_metadata.items()
            if endpoint_id not in predictions
        }
        percentiles: dict[str, float | None] = {}
        unknown_columns: dict[str, Any] = {}
        for endpoint_id, value in raw.items():
            if endpoint_id.endswith(PERCENTILE_SUFFIX):
                percentiles[endpoint_id.removesuffix(PERCENTILE_SUFFIX)] = value
            elif endpoint_id not in endpoint_metadata:
                unknown_columns[endpoint_id] = value

        missing = sorted(
            endpoint_id
            for endpoint_id, record in predictions.items()
            if finite_number(record.get("value")) is None
        )
        invalid = sorted(
            endpoint_id
            for endpoint_id, record in predictions.items()
            if finite_number(record.get("value")) is not None
            and endpoint_value_error(endpoint_id, record.get("value")) is not None
        )
        validation_issues = validate_native_predictions(predictions)
        status = "complete" if not validation_issues else "partial"
        results.append({
            "candidate_id": candidate_id,
            "smiles": candidate_smiles,
            "status": status,
            "endpoint_count": len(predictions),
            "missing_endpoints": missing,
            "invalid_endpoints": invalid,
            "endpoint_contract_issues": validation_issues,
            "predictions": predictions,
            "unregistered_predictions": unregistered_predictions,
            "drugbank_approved_percentiles": percentiles,
            "unmapped_columns": unknown_columns,
            "summary": _summary(predictions),
        })

    return {
        "success": True,
        "backend": {
            "id": "admet_ai",
            "version": importlib.metadata.version("admet-ai"),
            "device": "cpu",
            "python": platform.python_version(),
            "endpoint_schema_file": str(metadata_path),
            "endpoint_schema_sha256": schema_sha256,
            "provider_endpoint_count": len(endpoint_metadata),
            "registered_endpoint_count": len(ADMET_NATIVE_ENDPOINT_IDS),
            "registered_endpoint_ids": list(ADMET_NATIVE_ENDPOINT_IDS),
        },
        "requested": len(entries),
        "returned": len(results),
        "results": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="ADMET-AI batch JSON worker")
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output-file", required=True)
    args = parser.parse_args()
    output_path = Path(args.output_file)
    try:
        payload = json.loads(Path(args.input_file).read_text(encoding="utf-8"))
        result = predict(payload)
        exit_code = 0
    except Exception as exc:
        result = {
            "success": False,
            "error": f"{type(exc).__name__}: {exc}",
            "backend": {
                "id": "admet_ai",
                "device": "cpu",
                "python": platform.python_version(),
            },
        }
        exit_code = 1
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    print(json.dumps({"success": result["success"], "output_file": str(output_path)}))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

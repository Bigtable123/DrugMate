#!/usr/bin/env python
"""Batch ChemFH worker executed inside the dedicated ``chemfh`` environment."""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import random
import subprocess
import sys
from pathlib import Path
from typing import Any


os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ChemFH.settings")
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True


POSITIVE_THRESHOLD = 0.5
ENDPOINTS: tuple[tuple[str, str, float], ...] = (
    (
        "other_assay_interference",
        "Other assay interference",
        9.781904685012926e-05,
    ),
    ("blue_fluorescence", "Blue fluorescence", 0.0009872382265281),
    ("fluc_inhibitors", "FLuc inhibitors", 0.004291940348845),
    ("promiscuous_compounds", "Promiscuous compounds", 0.0061948686233965),
    ("colloidal_aggregators", "Colloidal aggregators", 0.0010918042703688),
    ("reactive_compounds", "Reactive compounds", 0.0016870113071982),
    ("green_fluorescence", "Green fluorescence", 0.0023240092127062),
)
RULE_NAMES: tuple[str, ...] = (
    "Aggregators",
    "Fluc",
    "Blue_fluorescence",
    "Green_fluorescence",
    "Reactive",
    "Other_assay_interference",
    "Promiscuous",
    "ALARM_NMR",
    "BMS",
    "Chelator_Rule",
    "GST_FHs_Rule",
    "His_FHs_Rule",
    "Luciferase_Inhibitor_Rule",
    "NTD",
    "PAINS",
    "Potential_Electrophilic_Rule",
    "Lilly",
)
SUBSTRUCTURE_RULE_NAMES: tuple[str, ...] = RULE_NAMES[7:]
SUBSTRUCTURE_RULE_MAPPING_VERSION = "chemfh-native-substructure-10-v1"


def _json_value(value: Any) -> Any:
    """Convert pandas/numpy/RDKit values to strict finite JSON values."""
    if value is None:
        return None
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes, dict)):
        try:
            value = value.tolist()
        except (TypeError, ValueError):
            pass
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
    if isinstance(value, (list, tuple, set)):
        return [_json_value(item) for item in value]
    return str(value)


def _finite_number(value: Any) -> float | None:
    value = _json_value(value)
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rule_bundle_sha256(rule_data_dir: Path) -> str:
    """Match the verified digest of sorted sha256sum records for rule data."""
    files = sorted(path.resolve() for path in rule_data_dir.rglob("*") if path.is_file())
    digest = hashlib.sha256()
    for path in files:
        digest.update(f"{_sha256(path)}  {path}\n".encode("utf-8"))
    return digest.hexdigest()


def _git_commit(source_dir: Path) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    commit = completed.stdout.strip()
    return commit if completed.returncode == 0 and len(commit) == 40 else None


def _package_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _verify_assets(
    source_dir: Path,
    checkpoint: Path,
    *,
    expected_checkpoint_sha256: str | None,
    expected_rule_bundle_sha256: str | None,
    expected_source_commit: str | None,
) -> dict:
    source_dir = source_dir.resolve()
    checkpoint = checkpoint.resolve()
    expected_checkpoint = (
        source_dir
        / "static"
        / "media"
        / "chemprop"
        / "checkpoint"
        / "fold_0"
        / "model_0"
        / "model.pt"
    ).resolve()
    rule_data_dir = (source_dir / "static" / "media" / "rule" / "data").resolve()
    required = {
        "ChemFH source": source_dir / "ChemFH" / "settings.py",
        "ChemFH predictor": (
            source_dir / "static" / "media" / "chemprop" / "scripts" / "predict.py"
        ),
        "ChemFH checkpoint": checkpoint,
        "ChemFH rule data": rule_data_dir,
    }
    missing = [label for label, path in required.items() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing required ChemFH resource(s): " + ", ".join(missing))
    if checkpoint != expected_checkpoint:
        raise ValueError(
            "ChemFH native predictor is pinned to "
            f"{expected_checkpoint}, not {checkpoint}"
        )

    checkpoint_sha256 = _sha256(checkpoint)
    if (
        expected_checkpoint_sha256
        and checkpoint_sha256 != expected_checkpoint_sha256
    ):
        raise ValueError(
            "ChemFH checkpoint SHA256 mismatch: "
            f"expected {expected_checkpoint_sha256}, got {checkpoint_sha256}"
        )
    rule_bundle_sha256 = _rule_bundle_sha256(rule_data_dir)
    if (
        expected_rule_bundle_sha256
        and rule_bundle_sha256 != expected_rule_bundle_sha256
    ):
        raise ValueError(
            "ChemFH rule bundle SHA256 mismatch: "
            f"expected {expected_rule_bundle_sha256}, got {rule_bundle_sha256}"
        )
    source_commit = _git_commit(source_dir)
    if expected_source_commit and source_commit != expected_source_commit:
        raise ValueError(
            "ChemFH source commit mismatch: "
            f"expected {expected_source_commit}, got {source_commit or 'unknown'}"
        )
    return {
        "source_dir": str(source_dir),
        "source_commit": source_commit,
        "expected_source_commit": expected_source_commit,
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": checkpoint_sha256,
        "expected_checkpoint_sha256": expected_checkpoint_sha256,
        "rule_data_dir": str(rule_data_dir),
        "rule_bundle_sha256": rule_bundle_sha256,
        "expected_rule_bundle_sha256": expected_rule_bundle_sha256,
    }


def _prepare_source(source_dir: Path, seed: int) -> None:
    source_dir = source_dir.resolve()
    os.chdir(source_dir)
    if str(source_dir) not in sys.path:
        sys.path.insert(0, str(source_dir))
    os.environ["DJANGO_SETTINGS_MODULE"] = "ChemFH.settings"

    random.seed(seed)
    import numpy as np

    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass

    import django

    django.setup()


def _native_rule_outputs(smiles: list[str]) -> tuple[Any, Any]:
    """Keep native rule evidence without generating its unused SVG tables."""
    import pandas as pd
    from static.media.rule.filter_rule import FrequentHitterFilter, readSmiles

    mols = readSmiles(smiles)
    screener = FrequentHitterFilter()
    smarts_columns: dict[str, list] = {}
    index_columns: dict[str, list] = {}
    for rule_name in RULE_NAMES:
        matches = screener.screening(mols, endpoint=rule_name)
        if len(matches) != len(smiles):
            raise ValueError(
                f"ChemFH {rule_name} cardinality mismatch: "
                f"requested {len(smiles)}, returned {len(matches)}"
            )
        # Match ChemFH._disposedRes exactly: preserve SMARTS order and the
        # first atom match for each SMARTS, without its depiction side effects.
        smarts_columns[f"{rule_name}_SMARTS"] = [
            " | ".join(result.keys()) for result in matches
        ]
        index_columns[f"{rule_name}_index"] = [
            [atoms[0] for atoms in result.values()] for result in matches
        ]
    return pd.DataFrame(smarts_columns), pd.DataFrame(index_columns)


def _native_outputs(smiles: list[str]) -> tuple[Any, Any, Any, Any]:
    from static.media.chemprop.scripts.predict import predict as native_predict

    predictions, uncertainties = native_predict(smiles)
    matched_smarts, atom_indices = _native_rule_outputs(smiles)
    return predictions, uncertainties, matched_smarts, atom_indices


def _frame_rows(frame: Any, expected: int, label: str) -> list[dict]:
    if not hasattr(frame, "iloc") or len(frame) != expected:
        returned = len(frame) if hasattr(frame, "__len__") else "unknown"
        raise ValueError(
            f"ChemFH {label} cardinality mismatch: requested {expected}, returned {returned}"
        )
    rows = []
    for position in range(expected):
        row = frame.iloc[position]
        if not hasattr(row, "to_dict"):
            raise TypeError(f"ChemFH {label} row {position} is not tabular")
        rows.append({str(key): value for key, value in row.to_dict().items()})
    return rows


def _literal_collection(value: Any) -> Any:
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        try:
            return ast.literal_eval(text)
        except (SyntaxError, ValueError):
            return value
    return value


def _atom_matches(value: Any) -> list[list[int]]:
    value = _literal_collection(_json_value(value))
    if not isinstance(value, list):
        return []
    matches: list[list[int]] = []
    for match in value:
        if not isinstance(match, (list, tuple)):
            continue
        normalized: list[int] = []
        valid = True
        for atom_index in match:
            try:
                normalized.append(int(atom_index))
            except (TypeError, ValueError):
                valid = False
                break
        if valid:
            matches.append(normalized)
    return matches


def _matched_smarts(value: Any) -> list[str]:
    value = _json_value(value)
    if isinstance(value, str):
        return [item.strip() for item in value.split("|") if item.strip()]
    if isinstance(value, list):
        return [str(item) for item in value if str(item).strip()]
    return []


def _candidate_result(
    candidate_id: str,
    smiles: str,
    prediction_row: dict,
    uncertainty_row: dict,
    smarts_row: dict,
    atom_index_row: dict,
) -> dict:
    predictions: dict[str, dict] = {}
    missing_endpoints: list[str] = []
    positive_endpoints: list[str] = []
    high_confidence_positive_endpoints: list[str] = []
    for endpoint_id, native_name, uncertainty_threshold in ENDPOINTS:
        probability = _finite_number(prediction_row.get(native_name))
        uncertainty = _finite_number(
            uncertainty_row.get(f"{native_name} uncertainty")
        )
        endpoint_complete = probability is not None and uncertainty is not None
        positive = (
            probability > POSITIVE_THRESHOLD if endpoint_complete else None
        )
        high_confidence = (
            uncertainty < uncertainty_threshold if endpoint_complete else None
        )
        if not endpoint_complete:
            missing_endpoints.append(endpoint_id)
        elif positive:
            positive_endpoints.append(endpoint_id)
            if high_confidence:
                high_confidence_positive_endpoints.append(endpoint_id)
        predictions[endpoint_id] = {
            "name": native_name,
            "probability": probability,
            "uncertainty": uncertainty,
            "positive": positive,
            "high_confidence": high_confidence,
            "probability_threshold": POSITIVE_THRESHOLD,
            "probability_operator": ">",
            "uncertainty_threshold": uncertainty_threshold,
            "uncertainty_operator": "<",
        }

    rules: dict[str, dict] = {}
    missing_rules: list[str] = []
    rule_hit_names: list[str] = []
    total_rule_matches = 0
    for rule_name in RULE_NAMES:
        index_key = f"{rule_name}_index"
        smarts_key = f"{rule_name}_SMARTS"
        if index_key not in atom_index_row:
            missing_rules.append(rule_name)
        atom_indices = _atom_matches(atom_index_row.get(index_key))
        matched_smarts = _matched_smarts(smarts_row.get(smarts_key))
        if atom_indices:
            rule_hit_names.append(rule_name)
            total_rule_matches += len(atom_indices)
        rule_record = {
            "match_count": len(atom_indices),
            "atom_indices": atom_indices,
            "matched_smarts": matched_smarts,
            "decision_contribution": False,
            "decision_role": "provider_evidence_only",
        }
        if rule_name == "PAINS":
            rule_record.update({
                "overlaps_with": "pains",
                "duplicate_decision_suppressed": True,
            })
        rules[rule_name] = rule_record

    status = "complete" if not missing_endpoints and not missing_rules else "partial"
    substructure_hit_names = [
        rule_name
        for rule_name in SUBSTRUCTURE_RULE_NAMES
        if rules[rule_name]["match_count"] > 0
    ]
    return {
        "candidate_id": candidate_id,
        "smiles": smiles,
        "status": status,
        "endpoint_count": len(predictions),
        "missing_endpoints": missing_endpoints,
        "predictions": predictions,
        "rule_count": len(rules),
        "missing_rules": missing_rules,
        "rules": rules,
        "summary": {
            "has_positive_prediction": bool(positive_endpoints),
            "positive_endpoint_count": len(positive_endpoints),
            "positive_endpoints": positive_endpoints,
            "high_confidence_positive_endpoints": high_confidence_positive_endpoints,
            "rule_hit_count": len(rule_hit_names),
            "rule_match_count": total_rule_matches,
            "rule_hit_names": rule_hit_names,
            # The executable policy describes ten substructure-rule catalogs.  The
            # provider also emits seven mechanism-named rules that correspond
            # to its seven model endpoints; keep those out of this aggregate so
            # they are not counted twice.
            "substructure_alert_count": len(substructure_hit_names),
            "substructure_alert_names": substructure_hit_names,
            "substructure_rule_names": list(SUBSTRUCTURE_RULE_NAMES),
            "substructure_rule_mapping_version": (
                SUBSTRUCTURE_RULE_MAPPING_VERSION
            ),
            "rules_contribute_to_decision": False,
            "pains_overlap_decision_suppressed": True,
        },
        "decision_policy": {
            "model_positive_is_review_floor": True,
            # Level-1 marking bounds only.  A ``fail`` mark here
            # is one metric on one molecule; the ALERT category fold in the
            # shared scorer is the only place a molecule is rejected.
            "probability_gate": {
                "pass_max_inclusive": 0.3,
                "flag_max_inclusive": 0.9,
                "fail_above": 0.9,
            },
            "uncertainty_changes_positive_call": False,
            # No individual rule carries its own gate.  Ten of them are
            # aggregated into the single ``chemfh.substructure_alert_count``
            # metric; the rest stay pure audit rows.
            "native_rules_are_evidence_only": True,
        },
    }


def predict(
    payload: dict,
    *,
    source_dir: Path,
    checkpoint: Path,
    expected_checkpoint_sha256: str | None,
    expected_rule_bundle_sha256: str | None,
    expected_source_commit: str | None,
    seed: int,
) -> dict:
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

    assets = _verify_assets(
        source_dir,
        checkpoint,
        expected_checkpoint_sha256=expected_checkpoint_sha256,
        expected_rule_bundle_sha256=expected_rule_bundle_sha256,
        expected_source_commit=expected_source_commit,
    )
    _prepare_source(source_dir, seed)
    prediction_frame, uncertainty_frame, smarts_frame, atom_index_frame = (
        _native_outputs(smiles)
    )
    prediction_rows = _frame_rows(prediction_frame, len(entries), "prediction")
    uncertainty_rows = _frame_rows(uncertainty_frame, len(entries), "uncertainty")
    smarts_rows = _frame_rows(smarts_frame, len(entries), "rule SMARTS")
    atom_index_rows = _frame_rows(atom_index_frame, len(entries), "rule indices")

    results = [
        _candidate_result(
            candidate_id,
            candidate_smiles,
            prediction_row,
            uncertainty_row,
            smarts_row,
            atom_index_row,
        )
        for (
            candidate_id,
            candidate_smiles,
            prediction_row,
            uncertainty_row,
            smarts_row,
            atom_index_row,
        ) in zip(
            candidate_ids,
            smiles,
            prediction_rows,
            uncertainty_rows,
            smarts_rows,
            atom_index_rows,
        )
    ]
    # The native rule loader can regenerate incompatible pickles. Refuse to
    # publish a run if any pinned scientific asset changed during execution.
    final_assets = _verify_assets(
        source_dir,
        checkpoint,
        expected_checkpoint_sha256=expected_checkpoint_sha256,
        expected_rule_bundle_sha256=expected_rule_bundle_sha256,
        expected_source_commit=expected_source_commit,
    )
    if final_assets != assets:
        raise ValueError("ChemFH scientific asset identity changed during prediction")

    return {
        "success": True,
        "backend": {
            "id": "chemfh",
            "version": expected_source_commit or assets.get("source_commit"),
            "device": "cpu",
            "python": platform.python_version(),
            "chemprop": _package_version("chemprop"),
            "django": _package_version("Django"),
            "rdkit": _package_version("rdkit"),
            "torch": _package_version("torch"),
            "pandas": _package_version("pandas"),
            "seed": int(seed),
            **assets,
            "endpoint_schema": [
                {
                    "id": endpoint_id,
                    "name": native_name,
                    "probability_threshold": POSITIVE_THRESHOLD,
                    "probability_operator": ">",
                    "uncertainty_threshold": uncertainty_threshold,
                    "uncertainty_operator": "<",
                }
                for endpoint_id, native_name, uncertainty_threshold in ENDPOINTS
            ],
            "rules": list(RULE_NAMES),
            "rule_depiction_enabled": False,
            "rule_decision_role": "provider_evidence_only",
            "pains_overlap_policy": (
                "retain ChemFH PAINS evidence; do not duplicate the RDKit pains gate"
            ),
        },
        "requested": len(entries),
        "returned": len(results),
        "results": results,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ChemFH batch JSON worker")
    parser.add_argument("--input-file")
    parser.add_argument("--output-file")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--expected-checkpoint-sha256")
    parser.add_argument("--expected-rule-bundle-sha256")
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--probe", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    source_dir = Path(args.source_dir)
    checkpoint = Path(args.checkpoint)
    payload = {
        "schema_version": 1,
        "candidates": [{"candidate_id": "chemfh-probe", "smiles": "CCO"}],
    }
    output_path = Path(args.output_file) if args.output_file else None
    try:
        if not args.probe:
            if not args.input_file or output_path is None:
                raise ValueError("--input-file and --output-file are required")
            payload = json.loads(Path(args.input_file).read_text(encoding="utf-8"))
        result = predict(
            payload,
            source_dir=source_dir,
            checkpoint=checkpoint,
            expected_checkpoint_sha256=args.expected_checkpoint_sha256,
            expected_rule_bundle_sha256=args.expected_rule_bundle_sha256,
            expected_source_commit=args.expected_source_commit,
            seed=int(args.seed),
        )
        exit_code = 0
    except Exception as exc:
        result = {
            "success": False,
            "error": f"{type(exc).__name__}: {exc}",
            "backend": {
                "id": "chemfh",
                "device": "cpu",
                "python": platform.python_version(),
                "source_dir": str(source_dir),
                "checkpoint": str(checkpoint),
            },
        }
        exit_code = 1

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False),
            encoding="utf-8",
        )
    print(json.dumps({
        "success": result["success"],
        "backend": "chemfh",
        "output_file": str(output_path) if output_path else None,
    }))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

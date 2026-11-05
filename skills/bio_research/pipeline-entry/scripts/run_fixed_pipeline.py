#!/usr/bin/env python3
"""Run one frozen SBDD configuration without an LLM controller.

Default: normalize optional inputs and validate only. --execute also calls the
existing pipeline once. Scientific choices must already be explicit in config.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
from datetime import datetime, timezone

ENTRY = Path(__file__).resolve().parents[1]
BIO_ROOT = ENTRY.parent
sys.path.insert(0, str(BIO_ROOT / "shared"))

from deployment_config import deployment_path
import pipeline_config as pconfig

DATA_ROOT = deployment_path("${SBDD_DATA_ROOT}")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _write(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def _policy_errors(config: dict) -> list[str]:
    errors = list(pconfig.explicit_selection_errors(config))
    project = config.get("project") or {}
    if not isinstance(project, dict):
        return [*errors, "project must be an object"]
    if project.get("execution_mode") != "auto":
        errors.append("fixed baseline requires explicit project.execution_mode='auto'")
    if project.get("round_number", 1) != 1:
        errors.append("fixed baseline starts at round_number=1 and runs one round")
    if project.get("data_root") != str(DATA_ROOT):
        errors.append(f"project.data_root must explicitly be {DATA_ROOT}")
    for key in ("upstream_project_dir", "project_dir"):
        if project.get(key) or config.get(key):
            errors.append(f"fixed baseline requires a fresh run; remove {key}")
    return errors


def _merge_normalization(config: dict, patch: dict) -> dict:
    """Accept input wiring, but never silently change a preselected route."""
    for step, update in patch.get("steps", {}).items():
        original = config["steps"][step]
        if "enabled" in update and original["enabled"] != update["enabled"]:
            raise ValueError(f"normalization conflicts with fixed step selection: {step}")
        for tool, tool_update in update.get("tools", {}).items():
            if original["tools"][tool]["enabled"] != tool_update["enabled"]:
                raise ValueError(f"normalization conflicts with fixed tool: {step}.{tool}")
    if "docking_box" in patch.get("input_artifacts", {}):
        if config["steps"]["docking"].get("params", {}).get("search_space_source") != "docking_box":
            raise ValueError("provided box requires predeclared search_space_source='docking_box'")
    return pconfig.deep_merge(config, patch)


def _input_hashes(config: dict, normalization: dict | None) -> dict:
    paths = list((config.get("input_artifacts") or {}).values())
    paths.extend(value for key, value in (config.get("context") or {}).items()
                 if key.endswith("_file") and value)
    paths.extend(section["path"] for section in (normalization or {}).values()
                 if isinstance(section, dict) and section.get("path"))
    identities = {}
    for value in paths:
        if not isinstance(value, str):
            raise ValueError("declared file inputs must be absolute path strings")
        path = Path(value)
        if not path.is_absolute() or not path.is_file():
            raise ValueError(f"declared input is not an existing absolute file: {value}")
        identities[str(path.resolve())] = _digest(path)
    return identities


def _report(output: Path, audit: dict) -> None:
    _write(output / "baseline_audit.json", audit)
    lines = [
        "# Fixed-config, no-LLM baseline", "",
        f"- Status: `{audit['status']}`",
        "- Controller: frozen configuration; no LLM controller invocation.",
        "- Outer rounds: one. Existing deterministic retries are retained.",
        f"- Pipeline invoked: `{audit['pipeline_invoked']}`",
        f"- Elapsed seconds: {audit['elapsed_seconds']}",
        f"- Scientific outcome: `{audit.get('scientific_status', 'not_evaluated')}`",
        f"- Valid for nomination: `{audit.get('valid_for_nomination', False)}`", "",
    ]
    if audit.get("error"):
        lines.extend(["## Error", "", str(audit["error"]), ""])
    manifest = output / "round1" / "manifest.json"
    if manifest.is_file():
        lines.extend(["## Pipeline stages", "", "| Stage | Status |", "|---|---|"])
        for step, result in _read(manifest).get("steps", {}).items():
            lines.append(f"| {step} | {result.get('status', 'unknown')} |")
        lines.append("")
    lines.extend([
        "## Evidence", "",
        "See `config_validation_report.json`, `frozen_config.json`, and `baseline_audit.json`.",
        "When executed, canonical outputs, runtime probes, logs and checkpoint are in `round1/`.",
        "A completed execution does not imply that the scientific objective was achieved.",
        "This template report is not an LLM-written scientific interpretation.", "",
    ])
    (output / "baseline_report.md").write_text("\n".join(lines), encoding="utf-8")


def run_fixed(config: dict, output_dir: Path, *, normalization: dict | None = None,
              execute: bool = False) -> dict:
    output = Path(output_dir).resolve()
    if output == DATA_ROOT.resolve() or not output.is_relative_to(DATA_ROOT.resolve()):
        raise ValueError(f"use a new experiment subdirectory under {DATA_ROOT}")
    # Atomic reservation: never overwrite another arm's outputs or prior audit.
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    audit = {
        "schema_version": 1, "controller": "fixed_config",
        "controller_llm_calls": 0, "cross_round_iteration": "disabled",
        "retry_policy": "unchanged_shared_runner_rules",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "output_dir": str(output), "audit_file": str(output / "baseline_audit.json"),
        "execute_requested": execute,
        "pipeline_invoked": False, "success": False, "status": "preparing",
    }
    try:
        config = copy.deepcopy(config)
        _write(output / "requested_config.json", config)
        audit["requested_config_sha256"] = _digest(output / "requested_config.json")
        errors = _policy_errors(config)
        if errors:
            audit.update(status="config_invalid", error="\n".join(errors))
            _write(output / "config_validation_report.json", {"valid": False, "errors": errors})
            return audit
        audit["source_input_sha256"] = _input_hashes(config, normalization)
        # Includes dirty source content, not just a git commit that misses edits.
        audit["implementation_sha256"] = {
            str(path.relative_to(BIO_ROOT)): _digest(path)
            for path in sorted(BIO_ROOT.rglob("*.py"))
            if "tests" not in path.parts and "__pycache__" not in path.parts
        }
        if normalization is not None:
            _write(output / "normalization_input.json", normalization)
            normalizer = _load("_fixed_baseline_normalization", BIO_ROOT / "input-normalization" / "run.py")
            normalized = normalizer.run(normalization, output / "normalization")
            if normalized.get("success") is not True:
                audit.update(status="normalization_failed", error=normalized.get("error", "normalization failed"))
                return audit
            config = _merge_normalization(config, normalized["config_patch"])
            audit["normalization_patch"] = normalized["config_patch"]
        validator = _load("_fixed_baseline_validator", ENTRY / "validate_config.py")
        validation = validator.validate_config(config)
        _write(output / "config_validation_report.json", validation)
        if validation.get("valid") is not True:
            audit.update(status="config_invalid", error="\n".join(validation.get("errors", [])))
            return audit
        # Resolve every parameter default now, before the experiment runs.
        frozen = pconfig.merge_user_config(config)
        frozen_file = output / "frozen_config.json"
        _write(frozen_file, frozen)
        _write(output / "resolved_run_plan.json", validation["resolved_plan"])
        audit["frozen_config_sha256"] = _digest(frozen_file)
        audit["execution_order"] = validation["resolved_plan"]["execution_order"]
        audit["normalized_input_sha256"] = _input_hashes(frozen, normalization)
        audit["pose_evaluation"] = (
            "scheduled" if "pose_evaluation" in audit["execution_order"]
            else "not_executed_by_fixed_configuration"
        )
        if not execute:
            audit.update(status="prepared", success=True)
            return audit
        if not audit["execution_order"]:
            raise ValueError("cannot execute an empty experiment; enable the required scientific stages")
        runner = _load("_fixed_baseline_pipeline", ENTRY / "run_pipeline.py")
        invocation = copy.deepcopy(frozen)
        invocation["project"]["project_dir"] = str(output / "round1")
        audit.update(status="running", pipeline_invoked=True)
        _write(output / "baseline_audit.json", audit)
        result = runner.run(invocation)
        # Match the normal CLI's handling of a --project-dir override.
        runner._archive_cli_input_config(result, frozen)
        _write(output / "pipeline_result.json", result)
        complete = (result.get("success") is True and result.get("pipeline_complete") is True
                    and result.get("status") == "complete" and not result.get("requires_confirmation"))
        audit.update(status=result.get("status", "unknown"), success=complete)
        for key in ("scientific_status", "decision_outcome", "goal_met", "valid_for_nomination", "error"):
            if key in result:
                audit[key] = result[key]
        # Do not loop on 'iterate', repair configurations, or replace tools.
        return audit
    except Exception as exc:
        audit.update(status="failed", success=False, error=f"{type(exc).__name__}: {exc}")
        return audit
    finally:
        audit["elapsed_seconds"] = round(time.monotonic() - started, 3)
        _report(output, audit)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="complete, predeclared pipeline config")
    parser.add_argument("--output-dir", required=True, type=Path, help="new baseline directory")
    parser.add_argument("--normalization-input", type=Path, help="optional normalization request JSON")
    parser.add_argument("--execute", action="store_true", help="also execute the reviewed fixed plan; default prepares and validates only")
    args = parser.parse_args()
    try:
        result = run_fixed(_read(args.config), args.output_dir,
                           normalization=_read(args.normalization_input) if args.normalization_input else None,
                           execute=args.execute)
    except (OSError, ValueError) as exc:
        result = {"success": False, "status": "input_error", "error": str(exc)}
    summary = {key: value for key, value in result.items() if key not in {
        "implementation_sha256", "source_input_sha256", "normalized_input_sha256", "normalization_patch",
    }}
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    raise SystemExit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

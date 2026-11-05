#!/usr/bin/env python3
"""Dry-run validator for a pipeline config — validate + resolve, never execute.

The Agent uses this as the self-check step in query→config: build a config with
an explicit step/tool matrix, validate it here, read the errors, fix, repeat —
all before spending compute in run_pipeline.py. It fills registered parameter
defaults, validates execution contracts, then prints the materialized run plan.

Usage:
    python pipeline-entry/validate_config.py --input '{...}'
    python pipeline-entry/validate_config.py --input-file config.json
    # Treat a resume/round-N config as satisfiable from an upstream round:
    python pipeline-entry/validate_config.py --input '{...}' --assume-upstream

Exit code 0 = valid; 1 = invalid (errors printed as JSON on stdout).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BIO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIO_ROOT / "shared"))

import pipeline_config as pconfig
import pipeline_validation as pvalidation
import pipeline_registry as R
import campaign

# Every canonical artifact — used when --assume-upstream lets a resume config
# treat disabled upstream steps as already satisfied.
ALL_ARTIFACTS = set(R.ARTIFACTS)


def validate_config(config: dict, assume_upstream: bool = False) -> dict:
    """Merge, validate, and resolve a config without executing anything.

    Returns a report dict including valid/errors, execution mode, the mandatory
    plan-confirmation signal, plan_summary, and the fully resolved plan.
    """
    if not isinstance(config, dict):
        return {
            "valid": False,
            "errors": ["pipeline config must be an object"],
            "assumed_upstream_artifacts": [],
            "requires_plan_confirmation": False,
            "execution_mode": "interrupt",
            "plan_summary": [],
            "resolved_plan": None,
        }

    try:
        merged = pconfig.merge_user_config(config)
    except ValueError as exc:
        return {
            "valid": False,
            "errors": [str(exc)],
            "assumed_upstream_artifacts": [],
            "requires_plan_confirmation": False,
            "execution_mode": (config.get("project") or {}).get(
                "execution_mode", "interrupt"
            ),
            "plan_summary": [],
            "resolved_plan": None,
        }

    available = set()
    upstream_configured = (
        config.get("upstream_project_dir")
        or (config.get("project") or {}).get("upstream_project_dir")
    )
    if assume_upstream or upstream_configured:
        # A resume run reuses upstream artifacts; the real orchestrator preloads
        # them from disk. For a static dry-run we optimistically assume every
        # artifact a *disabled* step would have produced is available upstream,
        # then apply the same producer-bundle invalidation semantics as an
        # explicit input_artifacts replacement.
        available = _artifacts_from_disabled_steps(merged)

    errors: list[str] = campaign.validate_campaign_config(config)
    try:
        pvalidation.validate(merged, available_artifacts=available)
    except pvalidation.ConfigValidationError as exc:
        errors.extend(exc.errors)

    plan = pconfig.resolve(merged)
    return {
        "valid": not errors,
        "errors": errors,
        "assumed_upstream_artifacts": sorted(available),
        "requires_plan_confirmation": not errors,
        "execution_mode": (plan.get("project") or {}).get(
            "execution_mode", "interrupt"
        ),
        "plan_summary": _plan_summary(plan),
        "resolved_plan": plan,
    }


def _artifacts_from_disabled_steps(merged: dict) -> set[str]:
    """Artifacts a disabled step may contribute from an assumed upstream run.

    Match ``ArtifactStore.inject`` bundle semantics: replacing one output from
    an upstream producer invalidates every old sibling output from that
    producer. Explicit companion injections are added back by ``validate``.
    """
    steps_cfg = merged.get("steps") or {}
    out: set[str] = set()
    for step_id, spec in R.STEPS.items():
        # The merge step has already enforced the complete explicit matrix.
        # Registry recommendations must never fill a missing execution switch.
        enabled = steps_cfg.get(step_id, {}).get("enabled", False)
        if not enabled:
            out.update(spec["produces"])

    injected = set((merged.get("input_artifacts") or {}).keys())
    invalidated: set[str] = set()
    for spec in R.STEPS.values():
        produced = set(spec["produces"])
        if injected & produced:
            invalidated.update(produced)

    # A real cached result records every canonical artifact wired into the
    # step.  Propagate invalidation through the registry DAG so the optimistic
    # dry-run cannot retain a descendant of a replaced upstream generation
    # (for example old prefiltered_library after injecting a new library).
    changed = True
    while changed:
        changed = False
        for spec in R.STEPS.values():
            inputs = set(spec.get("consumes", []))
            inputs.update(spec.get("optional_consumes", []))
            if not (inputs & invalidated):
                continue
            newly_invalid = set(spec["produces"]) - invalidated
            if newly_invalid:
                invalidated.update(newly_invalid)
                changed = True

    out.difference_update(invalidated)
    return out


def _plan_summary(plan: dict) -> list[dict]:
    summary = []
    for step_id in plan["execution_order"]:
        sp = plan["steps"][step_id]
        summary.append({
            "step": step_id,
            "skill": sp["skill"],
            "tools_enabled": sp["tools_enabled"],
            "tools_disabled": sp["tools_disabled"],
        })
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate + resolve a pipeline config (no execution)")
    parser.add_argument("--input", type=str, help="JSON config string")
    parser.add_argument("--input-file", type=str, help="Path to JSON config file")
    parser.add_argument("--assume-upstream", action="store_true",
                        help="Treat disabled upstream steps as satisfied (resume/round-N dry-run)")
    parser.add_argument("--quiet", action="store_true", help="Print only valid/errors, not the full plan")
    args = parser.parse_args()

    if args.input:
        config = json.loads(args.input)
    elif args.input_file:
        config = json.loads(Path(args.input_file).read_text())
    else:
        parser.error("Must provide --input or --input-file")

    report = validate_config(config, assume_upstream=args.assume_upstream)
    if args.quiet:
        report = {"valid": report["valid"], "errors": report["errors"]}
    print(json.dumps(report, indent=2, ensure_ascii=False, default=str))
    sys.exit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()

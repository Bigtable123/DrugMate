#!/usr/bin/env python3
"""Collect and compare measured campaign evidence; assemble an R3 input library.

No strategy is inferred and no threshold is relaxed. Chemistry normalization and
3D preparation reuse the existing pipeline implementations. Reports are read-only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

BIO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BIO_ROOT / "shared"))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def file_ref(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(data, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".campaign-")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def chemistry(smiles):
    from rdkit import Chem
    from pose_identity import build_chemical_identity
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("candidate is missing SMILES; names cannot establish chemical identity")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid candidate SMILES: {smiles!r}")
    identity = build_chemical_identity(mol, source_identity={"source": "campaign_evidence"})
    return identity["chemical_id"], identity["canonical_isomeric_smiles"]


def is_decoy(row):
    return (str(row.get("mol_name", "")).upper().startswith("DEC_")
            or str(row.get("candidate_role", "")).lower() in {"decoy", "benchmark_decoy"}
            or str(row.get("source", "")).lower() in {"decoy", "decoys"})


def eligible(row):
    return (row.get("classification") == "keep"
            and row.get("ranking_blocked") is not True
            and not str(row.get("nomination_status", "")).startswith("blocked"))


def protocol_valid(report):
    assessment = report.get("round_assessment") or {}
    summary = report.get("round_summary") or {}
    if assessment.get("status") in {"invalid", "not_evaluated"}:
        return False
    if (report.get("decision_outcome", summary.get("decision_outcome")) == "invalid"
            or report.get("decision") == "terminate_invalid"):
        return False
    if report.get("assessment_schema_version") == "round-assessment/v1":
        return assessment.get("protocol_valid") is True
    if (report.get("parameters") or {}).get("require_quality_gates", True):
        return report.get("quality_gate_status", summary.get("quality_gate_status")) == "keep"
    return True


def load_round(directory):
    directory = Path(directory).resolve()
    report_path = directory / "step6/coordinator_report.json"
    report = read_json(report_path)
    if report.get("success") is not True:
        raise ValueError(f"round assessment is not complete: {report_path}")
    rows = report.get("molecules")
    if not isinstance(rows, list):
        raise ValueError(f"round is missing its complete molecule audit: {report_path}")
    observed = []
    unavailable = []
    for row in rows:
        if is_decoy(row):
            continue
        try:
            chemical_id, smiles = chemistry(row.get("smiles_standardized") or row.get("smiles"))
        except ValueError as exc:
            # Rejected/unreadable molecules still have an audit entry; they
            # cannot silently become usable campaign input.
            unavailable.append({"mol_name": row.get("mol_name"), "reason": str(exc),
                                "classification": row.get("classification")})
            continue
        observed.append({**row, "campaign_candidate_id": chemical_id, "canonical_smiles": smiles})
    config_path = directory / "pipeline_config.json"
    return {"directory": str(directory), "report": report, "rows": observed,
            "unavailable": unavailable, "source": file_ref(report_path),
            "config": read_json(config_path) if config_path.exists() else None}


def collect(round_dirs):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    candidates = {}
    sources = []
    unavailable = []
    for directory in dict.fromkeys(str(Path(d).resolve()) for d in round_dirs):
        run = load_round(directory)
        sources.append(run["source"])
        unavailable.extend({**r, "round_dir": directory} for r in run["unavailable"])
        for row in run["rows"]:
            key = row["campaign_candidate_id"]
            entry = candidates.setdefault(key, {"candidate_id": key,
                "smiles": row["canonical_smiles"], "observations": []})
            observation = {"round_dir": directory, "mol_name": row.get("mol_name"),
                "chemical_identity_id": row.get("chemical_identity_id"),
                "pose_id": row.get("pose_id"), "classification": row.get("classification"),
                "eligible_keep": eligible(row) and protocol_valid(run["report"]),
                "protocol_valid": protocol_valid(run["report"]),
                "ranking_blocked_by": row.get("ranking_blocked_by") or [],
                "cohort_id": row.get("cohort_id"), "source": run["source"]}
            entry["observations"].append(observation)
    for entry in candidates.values():
        observations = entry["observations"]
        entry["pool_role"] = (
            "keep" if any(o["eligible_keep"] for o in observations) else
            "modify" if any(o["classification"] == "modify" and o["protocol_valid"] for o in observations) else
            "reject" if all(o["classification"] == "reject" for o in observations) else "review"
        )
        entry["requires_unified_reassessment"] = True
        scaffold = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(entry["smiles"]))
        entry["scaffold_smiles"] = Chem.MolToSmiles(scaffold, isomericSmiles=True)
    return {"schema_version": "campaign-evidence/v1", "sources": sources,
            "candidates": sorted(candidates.values(), key=lambda c: c["candidate_id"]),
            "unavailable_candidates": unavailable,
            "note": "KEEP is a union of valid source evidence, not a pooled cross-profile nomination."}


def finite(value):
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)


def raw_metrics(row):
    values = row.get("metric_assessments") or []
    if isinstance(values, dict):
        values = list(values.values())
    return {r["metric_id"]: r["value"] for r in values
            if isinstance(r, dict) and r.get("metric_id") and finite(r.get("value"))}


PROFILE_FIELDS = ("protein_state_id", "pocket_id", "protocol_profile_id", "metric_catalog_version")


def profile(row):
    values = tuple(str(row[k]) if k == "metric_catalog_version" and type(row.get(k)) is int
                   else row.get(k) for k in PROFILE_FIELDS)
    return values if all(isinstance(v, str) and v and not v.startswith("unknown:") for v in values) else None


def unique_rows(run):
    groups = {}
    for row in run["rows"]:
        groups.setdefault(row["campaign_candidate_id"], []).append(row)
    # A multi-pocket/multi-profile molecule has no single cross-round score.
    return {key: rows[0] for key, rows in groups.items() if len(rows) == 1}


def config_changes(before, after, prefix=""):
    if isinstance(before, dict) and isinstance(after, dict):
        out = []
        for key in sorted(set(before) | set(after)):
            out.extend(config_changes(before.get(key), after.get(key), f"{prefix}.{key}".strip(".")))
        return out
    return [] if before == after else [{"path": prefix, "before": before, "after": after}]


def pareto(rows, objectives):
    if not objectives or any(v not in {"min", "max"} for v in objectives.values()):
        raise ValueError("objectives must map at least one raw metric_id to min/max")
    pools = {}
    excluded = []
    for row in rows:
        values = raw_metrics(row)
        cohort = row.get("cohort_id")
        if (row.get("classification") == "reject" or row.get("ranking_blocked")
                or not cohort or not profile(row) or any(k not in values for k in objectives)):
            excluded.append(row["campaign_candidate_id"])
            continue
        vector = tuple(values[k] * (-1 if direction == "min" else 1) for k, direction in objectives.items())
        pools.setdefault((cohort, profile(row)), []).append((row["campaign_candidate_id"], vector))
    fronts = {}
    profiles = {}
    for (cohort, identity), points in pools.items():
        pool_id = str(cohort) + "::" + hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:16]
        profiles[pool_id] = {"cohort_id": cohort, **dict(zip(PROFILE_FIELDS, identity))}
        fronts[pool_id] = sorted({key for key, vector in points if not any(
            all(b >= a for a, b in zip(vector, other)) and any(b > a for a, b in zip(vector, other))
            for _, other in points)})
    return {"objectives": objectives, "fronts_by_cohort": fronts, "profiles": profiles,
            "excluded": sorted(set(excluded))}


def score_pool(run, identity):
    """The released ranking pool and measured metric mask must both match."""
    return sorted((row["campaign_candidate_id"], tuple(sorted(raw_metrics(row))))
                  for row in run["rows"] if profile(row) == identity
                  and row.get("classification") != "reject" and not row.get("ranking_blocked")
                  and finite(row.get("overall_rank")))


def compare(before_dir, after_dir, top_k=20, objectives=None):
    if top_k < 1:
        raise ValueError("top_k must be positive")
    before, after = load_round(before_dir), load_round(after_dir)
    a, b = unique_rows(before), unique_rows(after)
    all_before = {r["campaign_candidate_id"] for r in before["rows"]}
    all_after = {r["campaign_candidate_id"] for r in after["rows"]}
    common = set(a) & set(b)
    same_pool = set(a) == set(b) and len(a) == len(before["rows"]) and len(b) == len(after["rows"])
    changes = config_changes(before["config"], after["config"])
    # Profile ids identify protocols, not the entire scoring policy. Reject
    # relative-score comparison when the actual scoring/tool policy changed.
    scoring_changed = any(c["path"].startswith("steps.") for c in changes)
    pairs = []
    for key in sorted(common):
        x, y = a[key], b[key]
        comparable = (same_pool and not scoring_changed and before["config"] is not None
                      and after["config"] is not None and profile(x) is not None and profile(x) == profile(y)
                      and protocol_valid(before["report"]) and protocol_valid(after["report"])
                      and finite(x.get("overall_rank")) and finite(y.get("overall_rank"))
                      and x.get("classification") != "reject" and y.get("classification") != "reject"
                      and not x.get("ranking_blocked") and not y.get("ranking_blocked")
                      and score_pool(before, profile(x)) == score_pool(after, profile(y)))
        mx, my = raw_metrics(x), raw_metrics(y)
        pair = {"candidate_id": key, "before_classification": x.get("classification"),
                "after_classification": y.get("classification"), "rank_comparable": comparable,
                "rank_delta": y["overall_rank"] - x["overall_rank"]
                    if comparable and finite(x.get("overall_rank")) and finite(y.get("overall_rank")) else None,
                "raw_metric_changes": {k: {"before": mx[k], "after": my[k], "delta": my[k] - mx[k]}
                                       for k in sorted(set(mx) & set(my))}}
        pairs.append(pair)
    parent_child_pairs = []
    for child_key, child in b.items():
        parents = set()
        raw_parents = child.get("campaign_parent_ids_json")
        if raw_parents:
            parsed = json.loads(raw_parents) if isinstance(raw_parents, str) else raw_parents
            if not isinstance(parsed, list) or any(not isinstance(p, str) for p in parsed):
                raise ValueError("campaign_parent_ids_json must be an array of chemical ids")
            parents.update(parsed)
        if child.get("parent_smiles"):
            parents.add(chemistry(child["parent_smiles"])[0])
        for parent_key in sorted(parents & set(a)):
            if parent_key == child_key:
                continue
            parent = a[parent_key]
            mx, my = raw_metrics(parent), raw_metrics(child)
            parent_child_pairs.append({"parent_candidate_id": parent_key, "child_candidate_id": child_key,
                "same_protocol": profile(parent) is not None and profile(parent) == profile(child),
                "parent_classification": parent.get("classification"), "child_classification": child.get("classification"),
                "raw_metric_changes": {k: {"before": mx[k], "after": my[k], "delta": my[k] - mx[k]}
                                       for k in sorted(set(mx) & set(my))}})
    def top_set(run):
        groups = {}
        for row in run["rows"]:
            if eligible(row) and protocol_valid(run["report"]) and finite(row.get("overall_rank")):
                groups.setdefault(row.get("cohort_id"), []).append(row)
        # K per cohort, never rank unrelated pockets against each other.
        return {r["campaign_candidate_id"] for rows in groups.values()
                for r in rows if r["overall_rank"] <= top_k}
    ta, tb = top_set(before), top_set(after)
    result = {"schema_version": "campaign-comparison/v1", "sources": [before["source"], after["source"]],
        "top_k": top_k, "top_k_scope": "per_cohort_eligible_keep_with_ties",
        "top_k_overlap": {"before_count": len(ta), "after_count": len(tb), "intersection": len(ta & tb),
                          "jaccard": len(ta & tb) / len(ta | tb) if ta | tb else None},
        "pairs": pairs, "parent_child_pairs": parent_child_pairs,
        "new_candidates": sorted(all_after - all_before), "absent_candidates": sorted(all_before - all_after),
        "multiple_observations": {"before": sorted(all_before - set(a)), "after": sorted(all_after - set(b))},
        "parameter_changes": changes, "before_summary": before["report"].get("round_summary"),
        "after_summary": after["report"].get("round_summary"),
        "note": "Raw changes across protocols are descriptive, not evidence of improvement. Different pools/profiles have no rank delta."}
    if objectives is not None:
        result["pareto_before"] = pareto([
            {**r, "ranking_blocked": r.get("ranking_blocked") or not protocol_valid(before["report"])}
            for r in before["rows"]], objectives)
        result["pareto_after"] = pareto([
            {**r, "ranking_blocked": r.get("ranking_blocked") or not protocol_valid(after["report"])}
            for r in after["rows"]], objectives)
    return result


def assemble(archive_path, modifications, output_dir, seed=42, max_candidates=100):
    archive = read_json(archive_path)
    if archive.get("schema_version") != "campaign-evidence/v1":
        raise ValueError("unsupported candidate archive")
    # Sources are immutable evidence, not a cache silently refreshed on read.
    for ref in archive.get("sources", []):
        if file_ref(ref["path"])["sha256"] != ref["sha256"]:
            raise ValueError(f"candidate archive source changed: {ref['path']}")
    authoritative = collect([str(Path(ref["path"]).parent.parent) for ref in archive.get("sources", [])])
    verified = {c["candidate_id"]: c for c in authoritative["candidates"]}
    pool = {}
    aliases = {}
    for candidate in archive["candidates"]:
        original = verified.get(candidate.get("candidate_id"))
        if (original is None or original["pool_role"] != candidate.get("pool_role")
                or original["smiles"] != candidate.get("smiles")):
            raise ValueError("candidate archive role/chemistry does not match source evidence; select an eligible subset without changing gates")
        candidate = original
        if candidate["pool_role"] not in {"keep", "modify"}:
            continue
        key, smiles = chemistry(candidate["smiles"])
        if key != candidate["candidate_id"]:
            raise ValueError("candidate archive chemistry differs from its identity")
        pool[key] = {"compound_id": "CAMP_" + key.split(":")[-1], "campaign_candidate_id": key,
            "smiles": smiles, "smiles_standardized": smiles, "source": "campaign_retained",
            "candidate_role": "candidate", "parent_input_records": candidate["observations"],
            "campaign_parent_ids_json": "[]"}
        for observation in candidate["observations"]:
            if observation.get("mol_name"):
                aliases.setdefault(observation["mol_name"], set()).add(key)
    baseline = set(pool)
    modification_sources = []
    for path in modifications:
        report = read_json(path)
        if report.get("success") is not True:
            raise ValueError(f"modification report is unsuccessful: {path}")
        modification_ref = file_ref(path)
        modification_sources.append(modification_ref)
        for child in report.get("modifications", []):
            matches = aliases.get(child.get("parent_id"), set())
            if child.get("parent_smiles"):
                parent_key, _ = chemistry(child["parent_smiles"])
                matches = matches & {parent_key} if matches else ({parent_key} if parent_key in baseline else set())
            if len(matches) != 1:
                raise ValueError(f"missing or ambiguous modification parent: {child.get('parent_id')}")
            parent = next(iter(matches))
            key, smiles = chemistry(child["smiles"])
            comp = pool.setdefault(key, {"compound_id": "CAMP_" + key.split(":")[-1],
                "campaign_candidate_id": key, "smiles": smiles, "smiles_standardized": smiles,
                "source": "campaign_optimized", "candidate_role": "candidate", "parent_input_records": [],
                "parent_id": pool[parent]["compound_id"], "parent_smiles": pool[parent]["smiles"],
                "campaign_parent_ids_json": "[]"})
            parents = set(json.loads(comp["campaign_parent_ids_json"]))
            parents.add(parent)
            comp["campaign_parent_ids_json"] = json.dumps(sorted(parents))
            comp["parent_input_records"].append({"modification_report": modification_ref,
                "source_compound_id": child.get("compound_id"), "parent_candidate_id": parent,
                "rule_applied": child.get("rule_applied"), "instruction_recognized": child.get("instruction_recognized")})
    if not pool or len(pool) > max_candidates:
        raise ValueError(f"candidate pool size {len(pool)} is outside 1..{max_candidates}; explicitly revise selection/budget")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = [output_dir / name for name in ("library.sdf", "library_manifest.csv", "candidate_pool.json")]
    if any(path.exists() for path in targets):
        raise ValueError("assembly output already exists; use a new input revision directory")
    sys.path.insert(0, str(BIO_ROOT / "molecule-generator/scripts"))
    from conformer_utils import generate_conformers_sdf
    # A failure must never publish a partial input library as the complete pool.
    with tempfile.TemporaryDirectory(dir=output_dir, prefix=".assemble-") as temporary:
        staged = Path(temporary)
        compounds = [pool[key] for key in sorted(pool)]
        result = generate_conformers_sdf(compounds, staged / "library.sdf", random_seed=seed)
        if any(not comp.get("conformer_written") for comp in compounds):
            write_json(output_dir / "assembly_failure.json", result)
            raise ValueError("one or more conformers failed; see assembly_failure.json; no library published")
        fields = ["compound_id", "smiles", "source", "candidate_role", "parent_id", "parent_smiles",
                  "campaign_candidate_id", "campaign_parent_ids_json"]
        with (staged / "library_manifest.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(compounds)
        payload = {"schema_version": "campaign-pool/v1", "archive": file_ref(archive_path),
            "modifications": modification_sources, "seed": seed, "candidate_count": len(compounds),
            "requires_property_prefilter": True, "compounds": compounds,
            "files": {"library": str(targets[0].resolve()), "library_manifest": str(targets[1].resolve())}}
        write_json(staged / "candidate_pool.json", payload)
        for path in targets:
            os.replace(staged / path.name, path)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    collect_parser = commands.add_parser("collect")
    collect_parser.add_argument("--round-dir", action="append", required=True)
    collect_parser.add_argument("--output", required=True)
    compare_parser = commands.add_parser("compare")
    compare_parser.add_argument("--before", required=True)
    compare_parser.add_argument("--after", required=True)
    compare_parser.add_argument("--top-k", type=int, default=20)
    compare_parser.add_argument("--objectives")
    compare_parser.add_argument("--output", required=True)
    assemble_parser = commands.add_parser("assemble")
    assemble_parser.add_argument("--archive", required=True)
    assemble_parser.add_argument("--modifications", action="append", default=[])
    assemble_parser.add_argument("--output-dir", required=True)
    assemble_parser.add_argument("--seed", type=int, default=42)
    assemble_parser.add_argument("--max-candidates", type=int, default=100)
    args = parser.parse_args()
    try:
        if args.command == "collect":
            result = collect(args.round_dir)
            write_json(args.output, result)
        elif args.command == "compare":
            result = compare(args.before, args.after, args.top_k, read_json(args.objectives) if args.objectives else None)
            write_json(args.output, result)
        else:
            result = assemble(args.archive, args.modifications, args.output_dir, args.seed, args.max_candidates)
        print(json.dumps({"success": True, "output": getattr(args, "output", None) or args.output_dir,
                          "candidate_count": result.get("candidate_count", len(result.get("candidates", [])))}, ensure_ascii=False))
    except (ValueError, OSError, ImportError) as exc:
        print(json.dumps({"success": False, "error": str(exc)}, ensure_ascii=False))
        sys.exit(1)


if __name__ == "__main__":
    main()

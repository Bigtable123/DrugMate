#!/usr/bin/env python3
"""Pocket Analysis Skill - Standalone CLI.

Routes and fuses detectors using protein memory when available, analyzes protein-ligand
interactions, and computes fingerprints for downstream SBDD workflows.

Usage:
    python run.py --input-file explicit_pocket_config.json --output-dir ./output/
    python run.py --input-file params.json --output-dir ./output/

Both input files must declare ``site_selection_strategy`` and the complete
``tools.*.enabled`` matrix documented in SKILL.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
SHARED_DIR = SKILL_DIR.parent / "shared"
sys.path.insert(0, str(SKILL_DIR / "scripts"))
sys.path.insert(0, str(SHARED_DIR))
sys.path.insert(0, str(SKILL_DIR.parent / "pipeline-entry" / "scripts"))

from pocket_memory import (
    STRATEGY as LOCAL_FUSION_STRATEGY,
    MEMORY_STRATEGY,
    PAPER_TOOLS,
    plan_tools as plan_pocket_tools,
    make_box as make_fusion_box,
    normalize_payload,
    decide as decide_pocket,
)

from fpocket_utils import run_fpocket, parse_fpocket_output
from p2rank_utils import run_p2rank, parse_p2rank_output
from pykvfinder_utils import run_pykvfinder
from deeppocket_utils import run_deeppocket
from autosite_utils import run_autosite
from pocket_metrics import (
    build_pocket_box,
    compute_pocket_volume,
    compute_burial_ratio,
    classify_subpockets,
    compute_flexibility_zones,
    compute_iou,
    compute_consensus_pocket,
)
from ligand_site_utils import (
    build_box_site,
    build_ligand_site,
    box_from_coordinates,
    load_ligand_heavy_atom_coords,
)
from pocketminer_utils import run_pocketminer
from blueprint_utils import generate_blueprint
from interaction_observations import normalize_interaction_observations
from residue_utils import (
    deduplicate_residue_records,
    legacy_residue_label,
    parse_residue_record,
    residue_identity,
)
from io_utils import tool_config
from tool_runtime import ToolRuntimeError, build_command, build_process_env


import re


POCKET_TOOLS = (
    "ligand_site", "fpocket", "pykvfinder", "p2rank", "deeppocket", "autosite", "pocketminer",
    "plip", "prolif", "arpeggio", "luna", "oddt", "blueprint",
)

SITE_SELECTION_STRATEGIES = {
    "ligand_site",
    "known_site",
    "detector_consensus",
    "provided_box",
}


def _explicit_tool_switches(input_data: dict) -> dict[str, bool]:
    """Read a complete explicit pocket-tool enable matrix."""
    tools = input_data.get("tools")
    if not isinstance(tools, dict):
        raise ValueError("tools is required for explicit pocket tool selection")
    missing = [tool_id for tool_id in POCKET_TOOLS if tool_id not in tools]
    if missing:
        raise ValueError(
            "tools must explicitly configure every pocket tool; missing: "
            + ", ".join(missing)
        )
    switches = {}
    for tool_id in POCKET_TOOLS:
        config = tools[tool_id]
        if not isinstance(config, dict):
            raise ValueError(f"tools.{tool_id} must be an object")
        if not isinstance(config.get("enabled"), bool):
            raise ValueError(
                f"tools.{tool_id}.enabled must be an explicit boolean"
            )
        switches[tool_id] = config["enabled"]
    return switches

DEFAULT_PROLIF_INTERACTION_TYPES = [
    "Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking",
    "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact",
    "XBDonor", "XBAcceptor", "MetalDonor", "MetalAcceptor",
]


def _progress(message: str) -> None:
    print(f"[pocket-analysis] {message}", file=sys.stderr, flush=True)


def _format_native_failure(step: str, returncode: int, stderr: str) -> str:
    if returncode < 0:
        signal_number = -returncode
        try:
            signal_name = signal.Signals(signal_number).name
        except ValueError:
            signal_name = f"signal {signal_number}"
        reason = f"native worker terminated by {signal_name} ({signal_number})"
    else:
        reason = f"native worker exited with code {returncode}"

    stderr_lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    if stderr_lines:
        stderr_tail = " | ".join(stderr_lines[-8:])
        reason += f": {stderr_tail}"
    return f"{step} {reason}"


def _run_isolated_native_step(
    step: str,
    input_data: dict,
    timeout_seconds: int,
    runtime: dict | None = None,
) -> tuple[dict | None, str | None]:
    """Run a native profiler out of process so abort/segfault cannot kill Step1."""
    worker = SKILL_DIR / "scripts" / "native_worker.py"
    output_fd, output_name = tempfile.mkstemp(
        prefix=f"pocket_{step}_",
        suffix=".json",
    )
    os.close(output_fd)
    output_path = Path(output_name)
    worker_args = [
        str(worker), "--step", step, "--input", json.dumps(input_data),
        "--output", str(output_path),
    ]
    try:
        process_env = None
        if runtime:
            python_runtime = {
                "environment": runtime.get("environment", "current"),
                "command": "python",
            }
            command = build_command(python_runtime, worker_args)
            process_env = build_process_env(
                python_runtime, command, overrides=None,
            )
        else:
            command = [sys.executable, *worker_args]
    except ToolRuntimeError as exc:
        output_path.unlink(missing_ok=True)
        return None, f"{step} runtime resolution failed: {exc}"

    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
            env=process_env,
        )
    except subprocess.TimeoutExpired:
        output_path.unlink(missing_ok=True)
        return None, f"{step} native worker timed out after {timeout_seconds}s"
    except OSError as exc:
        output_path.unlink(missing_ok=True)
        return None, f"{step} native worker could not start: {exc}"

    try:
        payload = json.loads(output_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = None
    finally:
        output_path.unlink(missing_ok=True)

    if completed.returncode != 0:
        if payload and payload.get("error"):
            return None, payload["error"]
        return None, _format_native_failure(step, completed.returncode, completed.stderr)

    if not payload or not payload.get("ok"):
        return None, (payload or {}).get("error", f"{step} worker returned no result")
    return payload.get("result"), None


def _infer_resname_from_filename(ligand_file: Path) -> str:
    """Try to extract ligand resname from filename pattern: *_ligand_{COMP_ID}*.sdf/pdb"""
    stem = Path(ligand_file).stem
    m = re.search(r"_ligand_([A-Z0-9]{1,5})", stem, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return ""


_EVIDENCE_INPUTS = {
    "protonation_manifest": "protonation_manifest_file",
    "water_inventory": "water_inventory_file",
    "metal": "metal_report_file",
    "structure_qc": "structure_qc_report_file",
}


def _load_input_evidence(input_data: dict) -> dict:
    """Load optional Step-0 evidence without making any one format mandatory."""
    evidence = {}
    for name, field in _EVIDENCE_INPUTS.items():
        raw = input_data.get(field)
        if raw in (None, ""):
            evidence[name] = {"status": "not_provided", "path": None}
            continue
        path = Path(raw)
        entry = {"path": str(path)}
        if not path.is_file():
            entry.update({"status": "missing", "content": None})
            evidence[name] = entry
            continue
        try:
            payload = path.read_bytes()
        except OSError as exc:
            entry.update({"status": "unreadable", "error": str(exc), "content": None})
            evidence[name] = entry
            continue
        entry.update({
            "size_bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        })
        try:
            entry.update({"status": "loaded_json", "content": json.loads(payload)})
        except (UnicodeDecodeError, json.JSONDecodeError):
            entry.update({
                "status": "available",
                "content": None,
                "format": path.suffix.lower().lstrip(".") or "unknown",
            })
        evidence[name] = entry
    return evidence


def _selected_hetero_residues(
    input_evidence: dict, protein_file: Path | None = None,
) -> list:
    """Extract explicitly retained water/metal residue identities from reports."""
    selected = []

    water_resnames = {}
    if protein_file and Path(protein_file).is_file():
        try:
            with open(protein_file) as handle:
                for line in handle:
                    if not line.startswith("HETATM"):
                        continue
                    try:
                        key = (
                            line[21:22].strip(), int(line[22:26]),
                            line[26:27].strip(),
                        )
                    except ValueError:
                        continue
                    water_resnames[key] = line[17:20].strip()
        except OSError:
            water_resnames = {}

    def normalized_record(value: dict, context_type: str | None = None) -> dict | None:
        resnum = value.get("resnum")
        if resnum is None:
            resnum = value.get("residue_number")
        if resnum is None:
            resnum = value.get("res_seq")
        resname = (
            value.get("resname") or value.get("residue_name")
            or value.get("res_name") or value.get("metal")
        )
        if resnum is None or not resname:
            return None
        return {
            "chain": str(value.get("chain") or value.get("chain_id") or ""),
            "resnum": resnum,
            "icode": str(value.get("icode") or value.get("i_code") or ""),
            "resname": str(resname),
            "retained": True,
            **({"context_type": context_type} if context_type else {}),
        }

    # Structure-prep has stable, purpose-specific report contracts.  Parse
    # these first: water entries use ``kept`` rather than ``retained``, while
    # metal ``sites`` are observations made on the already prepared receptor.
    water_content = (input_evidence.get("water") or {}).get("content") or {}
    if isinstance(water_content, dict):
        for water in water_content.get("waters") or []:
            if isinstance(water, dict) and water.get("kept") is True:
                key = (
                    str(water.get("chain_id") or ""), water.get("res_seq"),
                    str(water.get("i_code") or ""),
                )
                selected.append(normalized_record({
                    **water,
                    "res_name": water.get("res_name") or water_resnames.get(key) or "HOH",
                }, "water"))

    metal_content = (input_evidence.get("metal") or {}).get("content") or {}
    if isinstance(metal_content, dict):
        for site in metal_content.get("sites") or []:
            if isinstance(site, dict) and site.get("res_seq") is not None and site.get("metal"):
                selected.append(normalized_record(site, "metal"))

    def visit(value, *, selected_context: bool = False):
        if isinstance(value, dict):
            chosen = (
                selected_context
                or value.get("selected") is True
                or value.get("retained") is True
                or value.get("kept") is True
            )
            decision = str(value.get("decision") or value.get("status") or "").lower()
            chosen = chosen or decision in {"keep", "kept", "retain", "retained", "selected"}
            record = normalized_record(value)
            has_identity = record is not None
            if chosen and has_identity:
                selected.append(record)
            propagate_selection = chosen and not has_identity
            for key, child in value.items():
                key_selected = any(token in str(key).lower() for token in ("selected", "retained", "kept"))
                visit(child, selected_context=propagate_selection or key_selected)
        elif isinstance(value, list):
            for child in value:
                visit(child, selected_context=selected_context)
        elif isinstance(value, str) and selected_context:
            selected.append(value)

    for name in ("water", "metal"):
        visit((input_evidence.get(name) or {}).get("content"))

    # Preserve order while removing records found by both the stable parser and
    # the generic compatibility walk.
    unique = []
    seen = set()
    for item in selected:
        if item is None:
            continue
        if isinstance(item, dict):
            key = (
                str(item.get("chain") or item.get("chain_id") or ""),
                str(item.get("resnum") or item.get("residue_number") or item.get("res_seq") or ""),
                str(item.get("icode") or item.get("i_code") or ""),
                str(item.get("resname") or item.get("residue_name") or item.get("res_name") or item.get("metal") or "").upper(),
            )
        else:
            key = (str(item),)
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique


def _finite_triplet(value, *, positive: bool = False) -> bool:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        return False
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            return False
        number = float(item)
        if not math.isfinite(number) or (positive and number <= 0):
            return False
    return True


def _make_docking_box(
    consensus: dict, protein_file: Path, chain: str, padding: float,
    *, site_selection_strategy: str | None = None,
    ligand_file: Path | None = None, min_box_size: float = 12.0,
) -> dict:
    if consensus.get("strategy") in (LOCAL_FUSION_STRATEGY, MEMORY_STRATEGY):
        return make_fusion_box(
            consensus, protein_file, padding=padding, min_size=min_box_size,
        )
    ligand_route = site_selection_strategy == "ligand_site" or (
        site_selection_strategy is None
        and consensus.get("source") == "co_crystal_ligand"
    )
    if ligand_route:
        # Ligand-site geometry has a single authority: the reference ligand.
        # Rebuild even valid cached boxes so detector residues and stale boxes
        # cannot silently enlarge or move this route's search space.
        box = {
            "schema_version": 1,
            "center": None, "size": None, "bounds": None,
            "unit": "angstrom", "source": "co_crystal_ligand",
            "geometry_source": "reference_ligand_heavy_atoms",
            "padding": padding, "min_size": min_box_size,
        }
        try:
            if ligand_file is None:
                raise ValueError("ligand_site box requires a reference ligand file")
            if not math.isfinite(padding) or padding < 0:
                raise ValueError("ligand_site box padding must be finite and non-negative")
            if not math.isfinite(min_box_size) or min_box_size < 0:
                raise ValueError("ligand_site minimum box size must be finite and non-negative")
            coords = load_ligand_heavy_atom_coords(ligand_file)
            if not coords or not all(_finite_triplet(point) for point in coords):
                raise ValueError("reference ligand has no valid finite heavy-atom geometry")
            center, size, bounds = box_from_coordinates(coords, padding=padding)
            box.update({
                "center": center,
                "size": [round(max(axis, min_box_size), 4) for axis in size],
                "bounds": bounds,
            })
        except Exception as exc:
            # Retain a diagnostic proposal; invalid geometry fails the existing
            # publication gate. Never fall back to receptor-residue geometry.
            box["error"] = str(exc)
        return box

    center = consensus.get("box_center")
    size = consensus.get("size")
    box_meta = consensus.get("box") or {}
    geometry_is_valid = _finite_triplet(center) and _finite_triplet(
        size, positive=True
    )
    embedded_padding = box_meta.get("padding")
    padding_changed = False
    try:
        padding_changed = not math.isclose(
            float(embedded_padding), float(padding), abs_tol=1e-9
        )
    except (TypeError, ValueError):
        pass
    # Detector adapters often attach a default-padded box before consensus.
    # An explicit run-level padding must still be authoritative; recompute from
    # the selected residues instead of silently retaining that stale geometry.
    if not geometry_is_valid or (
        padding_changed and (consensus.get("residues") or [])
    ):
        generated = build_pocket_box(
            protein_file,
            consensus.get("residues") or [],
            chain,
            residue_records=consensus.get("residue_records") or None,
            padding=padding,
        )
        generated_center = generated.get("center")
        generated_size = generated.get("size")
        if (
            _finite_triplet(generated_center)
            and _finite_triplet(generated_size, positive=True)
        ) or not geometry_is_valid:
            center = generated_center
            size = generated_size
            box_meta = generated
    return {
        "schema_version": 1,
        "center": center,
        "size": size,
        "unit": "angstrom",
        "source": consensus.get("source", "none"),
        "padding": box_meta.get("padding", padding),
        "bounds": box_meta.get("bounds"),
    }


def _known_site_metrics(detected_records: list[dict], known_values: list, chain: str) -> dict:
    known_records = deduplicate_residue_records(
        [
            record
            for value in known_values or []
            if (record := parse_residue_record(value, default_chain=chain)) is not None
        ],
        default_chain=chain,
    )
    detected = {
        identity
        for record in detected_records
        if (identity := residue_identity(record)) is not None
    }
    known = {
        identity
        for record in known_records
        if (identity := residue_identity(record)) is not None
    }
    overlap = detected & known
    union = detected | known
    precision = len(overlap) / len(detected) if detected else 0.0
    recall = len(overlap) / len(known) if known else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "iou_vs_known_site": round(len(overlap) / len(union), 3) if union else 0.0,
        "known_site_precision": round(precision, 3),
        "known_site_recall": round(recall, 3),
        "known_site_f1": round(f1, 3),
        "detected_residues_count": len(detected),
        "known_residues_count": len(known),
    }


def _select_known_site_candidate(
    consensus: dict,
    known_site_residues: list,
    chain: str,
) -> dict:
    """Select a geometry-bearing detector candidate by known-site recovery.

    Known residues guide selection but never manufacture geometry. Recall and
    overlap lead the ordering; independent detector support and native rank are
    deterministic tie-breakers only.
    """
    candidates = list(consensus.get("candidate_pockets") or [])
    if not candidates:
        result = dict(consensus)
        result["selection_method"] = "known_site_no_detector_candidates"
        return result

    annotated = []
    for original_index, candidate in enumerate(candidates):
        entry = dict(candidate)
        entry["known_site_metrics"] = _known_site_metrics(
            entry.get("residue_records") or [],
            known_site_residues,
            chain,
        )
        entry["_pre_known_site_order"] = original_index
        annotated.append(entry)

    def ranking_key(candidate: dict) -> tuple:
        metrics = candidate["known_site_metrics"]
        return (
            float(metrics.get("known_site_recall") or 0.0),
            float(metrics.get("iou_vs_known_site") or 0.0),
            float(metrics.get("known_site_f1") or 0.0),
            int(candidate.get("support_count") or 0),
            float(candidate.get("rank_score") or 0.0),
            -int(candidate.get("_pre_known_site_order") or 0),
        )

    selected = dict(max(annotated, key=ranking_key))
    selected.pop("_pre_known_site_order", None)
    for candidate in annotated:
        candidate.pop("_pre_known_site_order", None)
    selected["candidate_pockets"] = annotated
    selected["selection_method"] = (
        "known_site_recall_then_iou_then_detector_support_and_native_rank"
    )
    selected["known_site_selection"] = {
        "candidate_id": selected.get("candidate_id"),
        **dict(selected.get("known_site_metrics") or {}),
    }
    return selected


def _select_pocket_by_strategy(
    site_selection_strategy: str,
    fpocket_result: dict | None,
    p2rank_result: dict | None,
    ligand_file: Path | None,
    protein_file: Path,
    chain: str,
    *,
    pykvfinder_result: dict | None = None,
    deeppocket_result: dict | None = None,
    autosite_result: dict | None = None,
    ligand_site_result: dict | None = None,
    known_site_residues: list | None = None,
    top_k: int = 5,
    center_match_threshold: float = 5.0,
    residue_iou_threshold: float = 0.3,
    pocket_memory: dict | None = None,
    expected_pocket_plan: dict | None = None,
) -> dict:
    """Apply the declared strategy without inspecting artifact presence."""
    if site_selection_strategy == "detector_consensus":
        results = {
            "deeppocket": deeppocket_result, "fpocket": fpocket_result,
            "p2rank": p2rank_result, "pykvfinder": pykvfinder_result,
            "autosite": autosite_result,
        }
        plan = expected_pocket_plan or plan_pocket_tools(receptor_file=protein_file, **(pocket_memory or {}))
        results = {tool: results[tool] for tool in plan["execution_tools"]}
        digest = hashlib.sha256(protein_file.read_bytes()).hexdigest()
        if plan.get("receptor_sha256") not in (None, digest):
            raise ValueError("current receptor changed after pocket tool planning")
        for tool, result in results.items():
            if not result or result.get("success") is False:
                raise ValueError(f"{tool}: current detector candidates are required")
            claimed_sha = (result.get("input") or {}).get("protein_sha256")
            if (tool in {"deeppocket", "autosite"} or claimed_sha is not None) and claimed_sha != digest:
                raise ValueError(f"{tool}: input SHA-256 does not match the current receptor")
        payload = normalize_payload({
            "receptor_sha256": digest,
            "detectors": {tool: result.get("all_pockets", []) for tool, result in results.items()},
        }, protein_file)
        decision = decide_pocket(payload, receptor_file=protein_file, **(pocket_memory or {}))
        if expected_pocket_plan is not None:
            # Detection can take minutes. Reject changed history/query assets
            # instead of applying a different policy to the originally planned run.
            policy_keys = ("strategy", "execution_tools", "recommended_multiset",
                           "memory_used", "memory_mode", "memory_reason", "anchor_forced")
            decided_memory = {key: value for key, value in decision["memory"].items()
                              if key != "guard_values"}
            if (any(decision.get(key) != plan.get(key) for key in policy_keys)
                    or decided_memory != plan["memory"]):
                raise ValueError("pocket memory/query changed after detector planning; review and replan")
        members = decision["members"]
        return {
            **decision,
            "source": decision["strategy"],
            "selection_method": decision["strategy"],
            "candidate_id": "P001",
            "residues": [legacy_residue_label(record) for record in decision["residue_records"]],
            "detector_support": sorted({member["tool"] for member in members}),
            "support_count": decision["independent_support"],
            "raw_tool_support_count": len({member["tool"] for member in members}),
            "fusion_member_count": len(members),
            "independent_support_groups": sorted({member["family"] for member in members}),
            "point5_candidates": payload,
            **{
                f"selected_{member['tool']}": next(
                    pocket for pocket in results[member["tool"]]["all_pockets"]
                    if pocket["rank"] == member["rank"]
                )
                # Preserve every fusion member above, but expose the best native
                # rank per tool for legacy volume/annotation consumers.
                for member in sorted(members, key=lambda item: item["rank"], reverse=True)
            },
        }
    consensus = compute_consensus_pocket(
        fpocket_result,
        p2rank_result,
        ligand_file if site_selection_strategy == "ligand_site" else None,
        protein_file,
        chain,
        pykvfinder_result=pykvfinder_result,
        deeppocket_result=deeppocket_result,
        ligand_site_result=(
            ligand_site_result
            if site_selection_strategy == "ligand_site"
            else None
        ),
        top_k=top_k,
        center_match_threshold=center_match_threshold,
        residue_iou_threshold=residue_iou_threshold,
    )
    if site_selection_strategy == "known_site":
        return _select_known_site_candidate(
            consensus,
            known_site_residues or [],
            chain,
        )
    return consensus


def _pocketminer_enrichment(
    candidate_records: list,
    residue_scores: list[dict],
    *,
    min_score_coverage: float,
) -> dict:
    """Threshold-free Mann-Whitney AUC enrichment for one static candidate."""
    score_by_position = {}
    for record in residue_scores or []:
        identity = residue_identity(record, include_resname=False)
        try:
            score = float(record.get("score"))
        except (AttributeError, TypeError, ValueError):
            continue
        if identity is not None and math.isfinite(score):
            score_by_position[identity] = score

    candidate_positions = {
        identity
        for record in candidate_records or []
        if (identity := residue_identity(record, include_resname=False)) is not None
    }
    inside = [score_by_position[key] for key in candidate_positions if key in score_by_position]
    outside = [
        score for key, score in score_by_position.items() if key not in candidate_positions
    ]
    coverage = len(inside) / len(candidate_positions) if candidate_positions else 0.0
    auc = None
    if coverage >= min_score_coverage and inside and outside:
        wins = 0.0
        for positive in inside:
            for negative in outside:
                wins += 1.0 if positive > negative else 0.5 if positive == negative else 0.0
        auc = wins / (len(inside) * len(outside))
    return {
        "method": "mann_whitney_probability_of_superiority",
        "enrichment_auc": round(auc, 6) if auc is not None else None,
        "candidate_residue_count": len(candidate_positions),
        "scored_candidate_residue_count": len(inside),
        "score_coverage": round(coverage, 6),
        "min_score_coverage": min_score_coverage,
        "mean_candidate_score": round(sum(inside) / len(inside), 8) if inside else None,
        "max_candidate_score": round(max(inside), 8) if inside else None,
        "eligible_for_rerank": auc is not None,
        "threshold_used_for_auc": None,
    }


def _apply_pocketminer_evidence(
    consensus: dict,
    pocketminer_result: dict,
    *,
    site_selection_strategy: str,
    rerank_mode: str,
    min_score_coverage: float,
) -> dict:
    """Annotate cryptic evidence without changing the selected pocket."""
    residue_scores = pocketminer_result.get("residue_scores") or []
    candidates = list(consensus.get("candidate_pockets") or [])
    if not candidates and consensus.get("residue_records"):
        candidates = [{key: value for key, value in consensus.items()
                       if key != "candidate_pockets"}]

    annotated = []
    for candidate in candidates:
        entry = dict(candidate)
        entry["pocketminer_evidence"] = _pocketminer_enrichment(
            entry.get("residue_records") or entry.get("residues") or [],
            residue_scores,
            min_score_coverage=min_score_coverage,
        )
        annotated.append(entry)

    result = dict(consensus)
    selected_evidence = _pocketminer_enrichment(
        consensus.get("residue_records") or consensus.get("residues") or [],
        residue_scores,
        min_score_coverage=min_score_coverage,
    )
    result["pocketminer_evidence"] = selected_evidence
    result["pocketminer_role"] = "auxiliary_residue_level_cryptic_evidence"
    result["pocketminer_counts_as_detector_support"] = False

    result["candidate_pockets"] = annotated
    result["pocketminer_rerank"] = {
        "applied": False,
        "site_selection_strategy": site_selection_strategy,
        "reason": (
            "ligand_site remains ligand-defined"
            if site_selection_strategy == "ligand_site"
            else "known_site retains known-site validation priority"
            if site_selection_strategy == "known_site"
            else "fixed site selection; PocketMiner is annotation-only"
        ),
    }
    return result


def _validate_pocket(
    site_selection_strategy: str,
    consensus: dict,
    docking_box: dict,
    ligand_site_result: dict | None,
    detector_results: dict[str, dict | None],
    known_site_residues: list | None,
    chain: str,
) -> dict:
    records = deduplicate_residue_records(
        consensus.get("residue_records") or consensus.get("residues") or [],
        default_chain=chain,
    )
    box_valid = (
        _finite_triplet(docking_box.get("center"))
        and _finite_triplet(docking_box.get("size"), positive=True)
    )
    geometry_valid = bool(records) and box_valid
    detector_count = sum(
        1 for result in detector_results.values()
        if result and (result.get("all_pockets") or result.get("top_pocket"))
    )
    geometry_detector_ids = {"fpocket", "p2rank", "pykvfinder", "deeppocket", "autosite"}
    selected_geometry_sources = sorted(
        geometry_detector_ids & set(consensus.get("detector_support") or [])
    )
    selected_has_geometry_source = bool(selected_geometry_sources)
    base = {
        "site_selection_strategy": site_selection_strategy,
        "applicable": True,
        "not_applicable": False,
        "pocket_residue_count": len(records),
        "box_valid": box_valid,
        "detectors_with_candidates": detector_count,
        "selected_geometry_sources": selected_geometry_sources,
        "selected_has_geometry_source": selected_has_geometry_source,
    }

    if site_selection_strategy == "provided_box":
        return {
            **base,
            "passed": geometry_valid,
            "heavy_atom_count_inside_box": consensus.get("heavy_atom_count_inside_box", 0),
            "evidence_level": "provided_geometry",
            "pass_criterion": "finite positive box in angstrom; selected-chain heavy-atom overlap",
            "reason": (
                "provided box overlaps the selected receptor chain"
                if geometry_valid else "provided box has no usable overlap with the selected receptor chain"
            ),
            "scope": "box/receptor geometry only; not cavity, binding-site, or pose validation",
        }

    if site_selection_strategy == "ligand_site":
        ligand_residues = (ligand_site_result or {}).get("residue_records") or []
        ligand_center = (ligand_site_result or {}).get("center")
        pocket_center = consensus.get("center")
        center_valid = _finite_triplet(ligand_center) and _finite_triplet(pocket_center)
        center_distance = (
            math.sqrt(sum((float(ligand_center[i]) - float(pocket_center[i])) ** 2 for i in range(3)))
            if center_valid else None
        )
        ligand_coverage = 1.0 if ligand_site_result and _finite_triplet(
            (ligand_site_result or {}).get("box_center")
        ) and _finite_triplet((ligand_site_result or {}).get("size"), positive=True) else 0.0
        passed = (
            bool(ligand_residues)
            and consensus.get("source") == "co_crystal_ligand"
            and center_valid
            and center_distance <= 2.0
            and geometry_valid
            and ligand_coverage >= 0.95
        )
        return {
            **base,
            "passed": passed,
            "ligand_heavy_atom_box_coverage": ligand_coverage,
            "ligand_center_valid": center_valid,
            "pocket_to_ligand_center_distance": round(center_distance, 4) if center_distance is not None else None,
            "pass_criterion": (
                "ligand shell non-empty; pocket-to-ligand center distance <= 2A; "
                "ligand center and docking box valid; "
                "ligand heavy-atom box coverage >= 0.95"
            ),
            "reason": "co-crystal ligand site is valid" if passed else "invalid co-crystal ligand site",
        }

    if site_selection_strategy == "known_site":
        metrics = _known_site_metrics(records, known_site_residues or [], chain)
        passed = (
            geometry_valid
            and selected_has_geometry_source
            and metrics["known_site_recall"] >= 0.8
        )
        return {
            **base,
            **metrics,
            "passed": passed,
            "pass_criterion": (
                "known_site_recall >= 0.8; selected candidate comes from a registered "
                "geometry detector; docking-box geometry valid"
            ),
            "reason": "known site recovered" if passed else "known site recall or geometry below criterion",
        }

    support_count = int(consensus.get("support_count") or 0)
    passed = geometry_valid and selected_has_geometry_source
    evidence_level = "high" if support_count >= 3 else "medium" if support_count >= 2 else "low"
    return {
        **base,
        "passed": passed,
        "detector_support_count": support_count,
        "evidence_level": evidence_level,
        "pass_criterion": (
            "non-empty pocket; valid explicit box; selected candidate is supported by "
            "fpocket, P2Rank, pyKVFinder, or DeepPocket"
        ),
        "reason": (
            "detector-consensus pocket geometry is usable"
            if passed else "no usable detector-consensus pocket geometry"
        ),
        "warnings": ["single-detector pocket; review before prospective use"] if passed and support_count < 2 else [],
    }


def _make_pocket_definition(
    site_selection_strategy: str,
    consensus: dict,
    docking_box: dict,
    validation: dict,
) -> dict:
    support_count = int(consensus.get("support_count") or 0)
    if validation.get("passed") is not True:
        confidence = "low"
    elif site_selection_strategy == "provided_box":
        confidence = "not_assessed"
    elif site_selection_strategy == "ligand_site":
        confidence = "high"
    elif support_count >= 2 or site_selection_strategy == "known_site":
        confidence = "medium"
    else:
        confidence = "low"
    return {
        "schema_version": 1,
        "pocket_id": consensus.get("candidate_id", "P001"),
        "site_selection_strategy": site_selection_strategy,
        "source": consensus.get("source", "none"),
        "confidence": confidence,
        "center": consensus.get("center") or docking_box.get("center"),
        "size": docking_box.get("size"),
        "residues": consensus.get("residue_records") or [],
        "legacy_residue_labels": consensus.get("residues") or [],
        "detector_support": consensus.get("detector_support") or [],
        "support_count": support_count,
        "raw_tool_support_count": consensus.get("raw_tool_support_count", support_count),
        "independent_support_groups": (
            consensus.get("independent_support_groups") or []
        ),
        "selection_method": consensus.get("selection_method"),
        "anchor_tool": consensus.get("anchor_tool"),
        "anchor_rank": consensus.get("anchor_rank"),
        "memory_used": consensus.get("memory_used"),
        "memory_mode": consensus.get("memory_mode"),
        "memory_reason": consensus.get("memory_reason"),
        "memory": consensus.get("memory"),
        "guard_switched": consensus.get("guard_switched"),
        "guard_conditions": consensus.get("guard_conditions"),
        "active_tools": consensus.get("active_tools"),
        "unselected_tools": consensus.get("unselected_tools"),
        "recommended_multiset": consensus.get("recommended_multiset"),
        "execution_tools": consensus.get("execution_tools"),
        "inference_notes": consensus.get("inference_notes"),
        "anchor_forced": consensus.get("anchor_forced"),
        "implementation_status": consensus.get("implementation_status"),
        "detector_selection_policy": consensus.get("detector_selection_policy"),
        "residue_policy": (consensus.get("parameters") or {}).get("residue_policy"),
        "cryptic_evidence": consensus.get("pocketminer_evidence"),
        "cryptic_evidence_counts_as_detector_support": (
            consensus.get("pocketminer_counts_as_detector_support")
        ),
        "docking_box": docking_box,
        "validation": validation,
    }


def run(input_data: dict, output_dir: Path) -> dict:
    protein_file = Path(input_data["protein_file"])
    complex_mmcif_file = input_data.get("complex_mmcif_file")
    if complex_mmcif_file:
        complex_mmcif_file = Path(complex_mmcif_file)
    ligand_file = input_data.get("ligand_file")
    if ligand_file:
        ligand_file = Path(ligand_file)
    ligand_resname = input_data.get("ligand_resname", "")
    if not ligand_resname and ligand_file:
        ligand_resname = _infer_resname_from_filename(ligand_file)
    if not ligand_resname:
        ligand_resname = "LIG"
    chain = input_data.get("chain")
    known_site_residues = input_data.get("known_site_residues", None)
    raw_site_selection_strategy = input_data.get("site_selection_strategy")
    if raw_site_selection_strategy is None:
        return {
            "success": False,
            "error": "site_selection_strategy is required",
        }
    site_selection_strategy = str(raw_site_selection_strategy).strip()
    if site_selection_strategy not in SITE_SELECTION_STRATEGIES:
        return {
            "success": False,
            "error": (
                "site_selection_strategy must be one of "
                + ", ".join(sorted(SITE_SELECTION_STRATEGIES))
            ),
        }
    if site_selection_strategy != "provided_box" and (
        not isinstance(chain, str) or not chain.strip()
    ):
        return {"success": False, "error": "chain must explicitly identify the target after entity/assembly review; never assume A"}
    try:
        switches = _explicit_tool_switches(input_data)
    except ValueError as exc:
        return {"success": False, "error": str(exc)}
    pocket_memory_options = {
        "memory_mode": input_data.get("pocket_memory_mode", "auto"),
        "preset": input_data.get("pocket_memory_preset", "conservative"),
        "index_path": input_data.get("pocket_memory_index"),
        "query_npz": input_data.get("pocket_memory_query_npz"),
        "query_metadata": input_data.get("pocket_memory_query_metadata"),
        "query_case": input_data.get("pocket_memory_query_case"),
        "query_group": input_data.get("pocket_memory_query_group"),
        "exclude_cases": input_data.get("pocket_memory_exclude_cases") or (),
        "exclude_groups": input_data.get("pocket_memory_exclude_groups") or (),
    }
    pocket_plan = None
    if switches["autosite"] and site_selection_strategy != "detector_consensus":
        return {"success": False, "error": "AutoSite is registered for detector_consensus; disable it for explicit-site strategies"}
    if site_selection_strategy == "detector_consensus":
        try:
            pocket_plan = plan_pocket_tools(
                receptor_file=Path(input_data.get("protein_noH_file") or protein_file),
                **pocket_memory_options,
            )
            configured = {tool for tool in PAPER_TOOLS if switches[tool]}
            if configured != set(pocket_plan["execution_tools"]):
                raise ValueError(
                    "detector_consensus tool switches must match plan-pocket-tools.execution_tools; "
                    f"configured={sorted(configured)}, required={pocket_plan['execution_tools']}"
                )
        except (ValueError, KeyError, TypeError, OSError) as exc:
            for stale_name in ("pocket_definition.json", "docking_box.json"):
                (output_dir / stale_name).unlink(missing_ok=True)
            return {"success": False, "error": str(exc), "pocket_memory_plan": pocket_plan}
    known_required_interactions = input_data.get("known_required_interactions", None)
    known_site_provenance = input_data.get("known_site_provenance")
    blueprint_policy = dict(input_data.get("blueprint_policy") or {})
    input_evidence = _load_input_evidence(input_data)
    provided_box = None
    provided_box_file = None
    provided_box_site = None
    if site_selection_strategy == "provided_box":
        try:
            provided_box_file = Path(input_data["docking_box_file"]).resolve()
            box_bytes = provided_box_file.read_bytes()
            provided_box = json.loads(box_bytes)
            receptor_baseline = Path(input_data.get("protein_noH_file") or protein_file)
            provided_box_site = build_box_site(receptor_baseline, provided_box, chain)
            chain = provided_box_site["chain"]
            input_evidence["docking_box"] = {
                "path": str(provided_box_file),
                "sha256": hashlib.sha256(box_bytes).hexdigest(),
                "center": provided_box["center"], "size": provided_box["size"],
                "receptor_path": str(Path(input_data.get("protein_noH_file") or protein_file).resolve()),
                "receptor_sha256": hashlib.sha256(
                    Path(input_data.get("protein_noH_file") or protein_file).read_bytes()
                ).hexdigest(),
                "chain": chain,
                "chain_selection": provided_box_site["chain_selection"],
            }
        except (KeyError, TypeError, ValueError, OSError) as exc:
            # An invalid chain/box rerun must not leave an earlier canonical
            # pocket available to consumers. The supplied box stays untouched.
            for stale_name in ("pocket_definition.json", "docking_box.json"):
                stale_file = output_dir / stale_name
                if stale_file.resolve() != provided_box_file:
                    stale_file.unlink(missing_ok=True)
            return {"success": False, "error": f"invalid provided docking box: {exc}"}
    preserve_hetero_residues = _selected_hetero_residues(input_evidence, protein_file)
    prolif_interaction_types = (
        input_data.get("prolif_interaction_types")
        or list(DEFAULT_PROLIF_INTERACTION_TYPES)
    )
    fpocket_params = input_data.get("fpocket_params", {})
    pykvfinder_params = input_data.get("pykvfinder_params", {})
    ligand_site_params = input_data.get("ligand_site_params", {})
    deeppocket_config = tool_config(input_data, "deeppocket")
    deeppocket_params = {
        **dict(input_data.get("deeppocket_params") or {}),
        **dict(deeppocket_config.get("params") or {}),
    }
    deeppocket_runtime = dict(deeppocket_config.get("runtime") or {})
    if site_selection_strategy == "detector_consensus":
        deeppocket_params["max_candidates"] = 3
    pocketminer_config = tool_config(input_data, "pocketminer")
    pocketminer_params = {
        **dict(input_data.get("pocketminer_params") or {}),
        **dict(pocketminer_config.get("params") or {}),
    }
    pocketminer_runtime = dict(pocketminer_config.get("runtime") or {})
    profiler_configs = {
        tool_id: tool_config(input_data, tool_id)
        for tool_id in ("plip", "prolif", "arpeggio", "luna", "oddt")
    }
    profiler_runtimes = {
        tool_id: dict(config.get("runtime") or {})
        for tool_id, config in profiler_configs.items()
    }
    arpeggio_params = {
        **dict(input_data.get("arpeggio_params") or {}),
        **dict(profiler_configs["arpeggio"].get("params") or {}),
    }
    luna_params = {
        **{
            key: value for key, value in {
                "ifp_length": input_data.get("luna_ifp_length"),
                "ifp_type": input_data.get("luna_ifp_type"),
                "add_dependent_interactions": input_data.get(
                    "luna_add_dependent_interactions"
                ),
            }.items() if value is not None
        },
        **dict(input_data.get("luna_params") or {}),
        **dict(profiler_configs["luna"].get("params") or {}),
    }
    oddt_params = {
        **{
            key: value for key, value in {
                "strict": input_data.get("oddt_strict"),
                "splif_depth": input_data.get("oddt_splif_depth"),
                "splif_size": input_data.get("oddt_splif_size"),
                "splif_distance_cutoff": input_data.get(
                    "oddt_splif_distance_cutoff"
                ),
            }.items() if value is not None
        },
        **dict(input_data.get("oddt_params") or {}),
        **dict(profiler_configs["oddt"].get("params") or {}),
    }
    consensus_top_k = int(input_data.get("consensus_top_k", 3))
    center_match_threshold = float(input_data.get("center_match_threshold", 5.0))
    residue_iou_threshold = float(input_data.get("residue_iou_threshold", 0.3))
    if site_selection_strategy == "detector_consensus":
        consensus_top_k, center_match_threshold, residue_iou_threshold = 3, 5.0, 0.3
    box_padding = float(input_data.get("box_padding", ligand_site_params.get("box_padding", 4.0)))
    if site_selection_strategy == "detector_consensus":
        box_padding = float(input_data.get("detector_box_padding", box_padding))
    min_box_size = float(input_data.get("min_box_size", ligand_site_params.get("min_box_size", 12.0)))
    ligand_residue_cutoff = float(
        input_data.get("ligand_residue_cutoff", ligand_site_params.get("residue_cutoff", 5.0))
    )
    native_step_timeout = int(input_data.get("native_step_timeout_seconds", 300))

    if not protein_file.exists():
        return {"success": False, "error": f"Protein file not found: {protein_file}"}

    if ligand_file and not Path(ligand_file).exists():
        return {"success": False, "error": f"Ligand file not found: {ligand_file}"}

    ligand_dependent = [
        tool_id
        for tool_id in ("ligand_site", "plip", "arpeggio", "luna", "oddt")
        if switches[tool_id]
    ]
    if ligand_file is None and ligand_dependent:
        return {
            "success": False,
            "error": (
                "enabled tools require ligand_file: "
                + ", ".join(ligand_dependent)
            ),
        }
    if (
        switches["arpeggio"]
        and (complex_mmcif_file is None or not complex_mmcif_file.is_file())
    ):
        return {
            "success": False,
            "error": "enabled tool arpeggio requires complex_mmcif_file",
        }

    if site_selection_strategy == "ligand_site" and ligand_file is None:
        return {
            "success": False,
            "error": (
                "site_selection_strategy='ligand_site' requires ligand_file"
            ),
        }
    if site_selection_strategy == "ligand_site" and not switches["ligand_site"]:
        return {
            "success": False,
            "error": (
                "site_selection_strategy='ligand_site' requires the ligand_site "
                "tool enabled"
            ),
        }
    if site_selection_strategy == "known_site" and not known_site_residues:
        return {
            "success": False,
            "error": (
                "site_selection_strategy='known_site' requires "
                "known_site_residues"
            ),
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    detector_protein = (
        Path(input_data.get("protein_noH_file") or protein_file)
        if site_selection_strategy == "detector_consensus" else protein_file
    )

    steps_completed = []
    steps_skipped = []
    enabled_tool_failures = []
    parameters = {"site_selection_strategy": site_selection_strategy}

    # ===== Step 1: optional ligand-defined site tool =====
    ligand_site_result = None
    if not switches["ligand_site"]:
        steps_skipped.append({"step": "ligand_site", "reason": "disabled_by_config"})
    else:
        _progress("building co-crystal ligand site")
        try:
            ligand_site_result = build_ligand_site(
                protein_file,
                ligand_file,
                chain=chain,
                residue_cutoff=ligand_residue_cutoff,
                box_padding=box_padding,
                min_box_size=min_box_size,
            )
            steps_completed.append("ligand_site")
            parameters["ligand_site"] = {
                "residue_cutoff": ligand_residue_cutoff,
                "box_padding": box_padding,
                "min_box_size": min_box_size,
            }
        except Exception as exc:
            reason = str(exc)
            steps_skipped.append({"step": "ligand_site", "reason": reason})
            enabled_tool_failures.append({"tool": "ligand_site", "reason": reason})

    # ===== Step 2: geometric detectors =====
    fpocket_result = None
    if not switches["fpocket"]:
        steps_skipped.append({"step": "fpocket", "reason": "disabled_by_config"})
    else:
        _progress("starting fpocket")
        try:
            min_r = fpocket_params.get("min_radius", 3.0)
            max_r = fpocket_params.get("max_radius", 6.0)
            fpocket_out_dir = output_dir / "fpocket_output"
            run_fpocket(detector_protein, fpocket_out_dir, min_radius=min_r, max_radius=max_r)
            fpocket_result = parse_fpocket_output(fpocket_out_dir, detector_protein.stem)
            steps_completed.append("fpocket")
            parameters["fpocket"] = {
                "min_alpha_sphere_radius": min_r,
                "max_alpha_sphere_radius": max_r,
            }
        except Exception as e:
            reason = str(e)
            steps_skipped.append({"step": "fpocket", "reason": reason})
            enabled_tool_failures.append({"tool": "fpocket", "reason": reason})

    pykvfinder_result = None
    if not switches["pykvfinder"]:
        steps_skipped.append({"step": "pykvfinder", "reason": "disabled_by_config"})
    else:
        _progress("starting pyKVFinder")
        try:
            probe_in = float(pykvfinder_params.get("probe_in", 1.4))
            probe_out = float(pykvfinder_params.get("probe_out", 4.0))
            removal_distance = float(pykvfinder_params.get("removal_distance", 2.4))
            volume_cutoff = float(pykvfinder_params.get("volume_cutoff", 5.0))
            kvf_out_dir = output_dir / "pykvfinder_output"
            pykvfinder_result = run_pykvfinder(
                detector_protein, kvf_out_dir,
                probe_in=probe_in,
                probe_out=probe_out,
                removal_distance=removal_distance,
                volume_cutoff=volume_cutoff,
            )
            steps_completed.append("pykvfinder")
            parameters["pykvfinder"] = {
                "probe_in": probe_in,
                "probe_out": probe_out,
                "removal_distance": removal_distance,
                "volume_cutoff": volume_cutoff,
            }
        except Exception as e:
            reason = str(e)
            steps_skipped.append({"step": "pykvfinder", "reason": reason})
            enabled_tool_failures.append({"tool": "pykvfinder", "reason": reason})

    # ===== Step 3: P2Rank =====
    p2rank_result = None
    if not switches["p2rank"]:
        steps_skipped.append({"step": "p2rank", "reason": "disabled_by_config"})
    else:
        _progress("starting P2Rank")
        try:
            p2rank_out_dir = output_dir / "p2rank_output"
            run_p2rank(detector_protein, p2rank_out_dir)
            p2rank_result = parse_p2rank_output(p2rank_out_dir, detector_protein.name)
            steps_completed.append("p2rank")
            parameters["p2rank"] = {"model": "default"}
        except Exception as e:
            reason = str(e)
            steps_skipped.append({"step": "p2rank", "reason": reason})
            enabled_tool_failures.append({"tool": "p2rank", "reason": reason})

    # ===== Step 3a: optional DeepPocket CNN ranking + segmentation =====
    deeppocket_result = None
    deeppocket_file = output_dir / "deeppocket_predictions.json"
    if not switches["deeppocket"]:
        steps_skipped.append({"step": "deeppocket", "reason": "disabled_by_config"})
    else:
        _progress("starting isolated DeepPocket worker")
        try:
            deeppocket_result = run_deeppocket(
                detector_protein,
                deeppocket_file,
                runtime=deeppocket_runtime,
                params=deeppocket_params,
            )
        except Exception as exc:
            deeppocket_result = {
                "schema_version": 1,
                "success": False,
                "error_code": "deeppocket_launcher_error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        if deeppocket_result.get("success") is True:
            steps_completed.append("deeppocket")
            parameters["deeppocket"] = {
                "role": "fpocket_proposal_cnn_ranking_and_segmentation",
                "max_candidates": int(deeppocket_params.get("max_candidates", 3)),
                "segmentation_threshold": float(
                    deeppocket_params.get("segmentation_threshold", 0.5)
                ),
                "mask_distance": float(
                    deeppocket_params.get("mask_distance", 3.5)
                ),
                "independence_group": "fpocket_family",
            }
        else:
            reason = (
                deeppocket_result.get("error")
                or "DeepPocket returned an unsuccessful normalized contract"
            )
            steps_skipped.append({
                "step": "deeppocket",
                "reason": reason,
                "error_code": deeppocket_result.get("error_code"),
            })
            enabled_tool_failures.append({"tool": "deeppocket", "reason": reason})

    autosite_result = None
    if not switches["autosite"]:
        steps_skipped.append({"step": "autosite", "reason": "disabled_by_config"})
    else:
        _progress("starting AutoSite native cluster detection")
        autosite_config = tool_config(input_data, "autosite")
        try:
            autosite_result = run_autosite(
                detector_protein, output_dir / "autosite_output",
                runtime=dict(autosite_config.get("runtime") or {}),
                timeout=int((autosite_config.get("params") or {}).get("timeout_seconds", 1200)),
            )
            steps_completed.append("autosite")
            parameters["autosite"] = autosite_result["parameters"]
        except Exception as exc:
            reason = str(exc)
            steps_skipped.append({"step": "autosite", "reason": reason})
            enabled_tool_failures.append({"tool": "autosite", "reason": reason})

    # ===== Step 3b: optional residue-level cryptic evidence =====
    pocketminer_result = None
    pocketminer_error = None
    pocketminer_file = output_dir / "pocketminer_predictions.json"
    if not switches["pocketminer"]:
        steps_skipped.append({"step": "pocketminer", "reason": "disabled_by_config"})
    else:
        _progress("starting strict PocketMiner worker")
        try:
            pocketminer_result = run_pocketminer(
                protein_file,
                chain,
                pocketminer_file,
                runtime=pocketminer_runtime,
                params=pocketminer_params,
            )
        except Exception as exc:
            pocketminer_result = {
                "schema_version": 1,
                "success": False,
                "error_code": "pocketminer_launcher_error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        if pocketminer_result.get("success") is True:
            steps_completed.append("pocketminer")
            parameters["pocketminer"] = {
                "role": "auxiliary_residue_level_cryptic_evidence",
                "rerank_mode": pocketminer_params.get("rerank_mode", "annotate"),
                "min_score_coverage": float(
                    pocketminer_params.get("min_score_coverage", 0.8)
                ),
                "provisional_quantile": float(
                    pocketminer_params.get("provisional_quantile", 0.9)
                ),
                "counts_as_detector_support": False,
            }
        else:
            pocketminer_error = (
                pocketminer_result.get("error")
                or "strict PocketMiner execution returned an unsuccessful contract"
            )
            steps_skipped.append({
                "step": "pocketminer",
                "reason": pocketminer_error,
                "error_code": pocketminer_result.get("error_code"),
            })
            enabled_tool_failures.append({
                "tool": "pocketminer", "reason": pocketminer_error,
            })

    # ===== Consensus pocket =====
    _progress(f"selecting pocket with {site_selection_strategy}")
    if provided_box is not None:
        # Keep detector outputs as supplemental evidence; chain/box overlap was
        # checked before running any optional analysis tools.
        consensus = provided_box_site
    else:
        try:
            consensus = _select_pocket_by_strategy(
                site_selection_strategy,
                fpocket_result,
                p2rank_result,
                ligand_file,
                detector_protein,
                chain,
                pykvfinder_result=pykvfinder_result,
                autosite_result=autosite_result,
                deeppocket_result=(
                    deeppocket_result
                    if deeppocket_result and deeppocket_result.get("success") is True
                    else None
                ),
                ligand_site_result=ligand_site_result,
                known_site_residues=known_site_residues,
                top_k=consensus_top_k,
                center_match_threshold=center_match_threshold,
                residue_iou_threshold=residue_iou_threshold,
                pocket_memory=pocket_memory_options,
                expected_pocket_plan=pocket_plan,
            )
        except (ValueError, KeyError, TypeError, OSError) as exc:
            if site_selection_strategy != "detector_consensus":
                raise
            consensus = {"source": "none", "error": str(exc)}
            enabled_tool_failures.append({"tool": "pocket_selection", "reason": str(exc)})
    if pocketminer_result and pocketminer_result.get("success") is True:
        consensus = _apply_pocketminer_evidence(
            consensus,
            pocketminer_result,
            site_selection_strategy=site_selection_strategy,
            rerank_mode=str(pocketminer_params.get("rerank_mode", "annotate")),
            min_score_coverage=float(pocketminer_params.get("min_score_coverage", 0.8)),
        )
    parameters["consensus"] = {
        "top_k": consensus_top_k,
        "match_operator": "and",
        "center_match_threshold": center_match_threshold,
        "residue_iou_threshold": residue_iou_threshold,
    }

    # ===== Step 3: Interaction Analysis (PLIP primary, ProLIF secondary) =====
    plip_result = None
    prolif_result = None
    interaction_tool_statuses = {}
    profiler_failures = []

    # 3a: PLIP. OpenBabel can abort the interpreter, so this optional profiler
    # must run outside the main Step1 process.
    if not switches["plip"]:
        steps_skipped.append({"step": "plip", "reason": "disabled_by_config"})
        interaction_tool_statuses["plip"] = {
            "status": "not_evaluable", "reason": "disabled_by_config", "enabled": False,
        }
    else:
        _progress("starting isolated PLIP worker")
        plip_result, plip_error = _run_isolated_native_step(
            "plip",
            {
                "protein_file": str(protein_file),
                "ligand_file": str(ligand_file),
                "ligand_resname": ligand_resname,
                "chain": chain,
            },
            native_step_timeout,
            runtime=profiler_runtimes["plip"],
        )
        if plip_result and plip_result.get("success"):
            steps_completed.append("plip")
            parameters["plip"] = {"role": "primary_interaction_profiler"}
            interaction_tool_statuses["plip"] = {
                "status": "detected" if plip_result.get("num_interactions_total") else "not_detected",
                "enabled": True,
                "tool_version": plip_result.get("tool_version"),
                "interaction_count": int(plip_result.get("num_interactions_total") or 0),
            }
        else:
            reason = plip_error or (plip_result or {}).get("error", "unknown")
            steps_skipped.append({"step": "plip", "reason": reason})
            interaction_tool_statuses["plip"] = {
                "status": "tool_failed", "reason": reason, "enabled": True,
            }
            profiler_failures.append({"tool": "plip", "reason": reason})
            plip_result = None
    # 3b: ProLIF. MDAnalysis-to-RDKit conversion can segfault, so isolate it
    # for the same reason as PLIP.
    if not switches["prolif"]:
        steps_skipped.append({"step": "prolif", "reason": "disabled_by_config"})
        interaction_tool_statuses["prolif"] = {
            "status": "not_evaluable", "reason": "disabled_by_config", "enabled": False,
        }
    elif ligand_file is None:
        steps_skipped.append({"step": "prolif", "reason": "reference_ligand_missing"})
        interaction_tool_statuses["prolif"] = {
            "status": "not_evaluable", "reason": "reference_ligand_missing", "enabled": True,
        }
    else:
        _progress("starting isolated ProLIF worker")
        # Chemical perception and directional H-bonds require the prepared,
        # hydrogenated receptor.  A no-H receptor is not an equivalent input.
        prolif_protein = protein_file
        prolif_result, prolif_error = _run_isolated_native_step(
            "prolif",
            {
                "protein_file": str(prolif_protein),
                "ligand_file": str(ligand_file),
                "ligand_resname": ligand_resname,
                "interaction_types": prolif_interaction_types,
                "vicinity_cutoff": float(input_data.get("prolif_vicinity_cutoff", 6.5)),
                "preserve_hetero_residues": preserve_hetero_residues,
            },
            native_step_timeout,
            runtime=profiler_runtimes["prolif"],
        )
        if prolif_result:
            protein_h_used = prolif_result.get("protein_h_used")
            parameters["prolif"] = {
                "interaction_types": (
                    prolif_result.get("capabilities", {}).get("enabled_pairwise_interaction_types")
                    or prolif_interaction_types
                    or []
                ),
                "capabilities": prolif_result.get("capabilities") or {},
                "protein_h_used": protein_h_used,
                "protein_hydrogen_count": prolif_result.get("protein_hydrogen_count"),
                "protein_hydrogen_count_scope": prolif_result.get("protein_hydrogen_count_scope"),
                "protein_file_used": prolif_result.get("protein_file_used"),
                "used_noH": None if protein_h_used is None else not protein_h_used,
                "input_strategy": prolif_result.get("input_strategy", "unknown"),
                "requires_final_state_review": True,
                "preserved_hetero_residues": preserve_hetero_residues,
            }
        if prolif_result and prolif_result.get("validation_failed"):
            reason = prolif_result.get("diagnostics", {}).get("error", "validation_failed")
            steps_skipped.append({"step": "prolif", "reason": reason})
            interaction_tool_statuses["prolif"] = {
                "status": "tool_failed", "reason": reason, "enabled": True,
                "tool_version": prolif_result.get("tool_version"),
            }
            profiler_failures.append({"tool": "prolif", "reason": reason})
        elif prolif_result:
            steps_completed.append("prolif")
            interaction_tool_statuses["prolif"] = {
                "status": "detected" if prolif_result.get("num_interactions_total") else "not_detected",
                "enabled": True,
                "tool_version": prolif_result.get("tool_version"),
                "interaction_count": int(prolif_result.get("num_interactions_total") or 0),
                "protein_input_strategy": parameters["prolif"]["input_strategy"],
            }
        else:
            reason = prolif_error or "unknown"
            steps_skipped.append({"step": "prolif", "reason": reason})
            interaction_tool_statuses["prolif"] = {
                "status": "tool_failed", "reason": reason, "enabled": True,
            }
            profiler_failures.append({"tool": "prolif", "reason": reason})
    # 3c-e: optional supplemental profilers.  They enrich the observation and
    # fingerprint evidence but never enter the strict PLIP/ProLIF consensus.
    # Once explicitly enabled, an operational failure invalidates the Step-1
    # execution contract; chemistry-quality labels remain scientific evidence.
    supplemental_results: dict[str, dict | None] = {
        "arpeggio": None,
        "luna": None,
        "oddt": None,
    }
    supplemental_profiler_failures = []

    if not switches["arpeggio"]:
        steps_skipped.append({"step": "arpeggio", "reason": "disabled_by_config"})
        interaction_tool_statuses["arpeggio"] = {
            "status": "not_evaluable", "reason": "disabled_by_config", "enabled": False,
        }
    else:
        _progress("starting isolated Arpeggio worker")
        result, error = _run_isolated_native_step(
            "arpeggio",
            {
                "complex_mmcif_file": str(complex_mmcif_file),
                "ligand_file": str(ligand_file),
                "ligand_resname": ligand_resname,
                "chain": chain,
                "tool_timeout_seconds": int(
                    arpeggio_params.get("timeout_seconds", native_step_timeout)
                ),
            },
            native_step_timeout,
            runtime=profiler_runtimes["arpeggio"],
        )
        supplemental_results["arpeggio"] = result
        if result and result.get("success"):
            steps_completed.append("arpeggio")
            parameters["arpeggio"] = {
                "role": "supplemental_atom_contact_evidence",
                "source_complex_mmcif": str(complex_mmcif_file),
                "quality": result.get("quality"),
            }
            interaction_tool_statuses["arpeggio"] = {
                "status": "detected" if result.get("num_interactions_total") else "not_detected",
                "enabled": True,
                "tool_version": result.get("tool_version"),
                "interaction_count": int(result.get("num_interactions_total") or 0),
                "quality": result.get("quality"),
                "quality_flags": result.get("quality_flags") or [],
            }
        elif result and result.get("not_evaluable"):
            reason = result.get("error") or "Arpeggio input is not evaluable"
            steps_skipped.append({"step": "arpeggio", "reason": reason})
            interaction_tool_statuses["arpeggio"] = {
                "status": "not_evaluable", "reason": reason, "enabled": True,
                "tool_version": result.get("tool_version"),
            }
        else:
            reason = error or (result or {}).get("error") or "unknown"
            steps_skipped.append({"step": "arpeggio", "reason": reason})
            interaction_tool_statuses["arpeggio"] = {
                "status": "tool_failed", "reason": reason, "enabled": True,
            }
            supplemental_profiler_failures.append({"tool": "arpeggio", "reason": reason})

    if not switches["luna"]:
        steps_skipped.append({"step": "luna", "reason": "disabled_by_config"})
        interaction_tool_statuses["luna"] = {
            "status": "not_evaluable", "reason": "disabled_by_config", "enabled": False,
        }
    else:
        _progress("starting isolated LUNA worker")
        result, error = _run_isolated_native_step(
            "luna",
            {
                "protein_file": str(protein_file),
                "ligand_file": str(ligand_file),
                "ligand_resname": ligand_resname,
                "ifp_length": int(luna_params.get("ifp_length", 4096)),
                "ifp_type": str(luna_params.get("ifp_type", "EIFP")),
                "add_dependent_interactions": bool(
                    luna_params.get("add_dependent_interactions", False)
                ),
            },
            native_step_timeout,
            runtime=profiler_runtimes["luna"],
        )
        supplemental_results["luna"] = result
        if result and result.get("success"):
            steps_completed.append("luna")
            parameters["luna"] = {
                "role": "supplemental_atom_aware_ifp",
                "ifp_length": int(luna_params.get("ifp_length", 4096)),
                "ifp_type": str(luna_params.get("ifp_type", "EIFP")),
                "add_dependent_interactions": bool(
                    luna_params.get("add_dependent_interactions", False)
                ),
            }
            interaction_tool_statuses["luna"] = {
                "status": "detected" if result.get("num_interactions_total") else "not_detected",
                "enabled": True,
                "tool_version": result.get("tool_version"),
                "interaction_count": int(result.get("num_interactions_total") or 0),
                "quality_flags": result.get("quality_flags") or [],
            }
        else:
            reason = error or (result or {}).get("error") or "unknown"
            steps_skipped.append({"step": "luna", "reason": reason})
            interaction_tool_statuses["luna"] = {
                "status": "tool_failed", "reason": reason, "enabled": True,
            }
            supplemental_profiler_failures.append({"tool": "luna", "reason": reason})

    if not switches["oddt"]:
        steps_skipped.append({"step": "oddt", "reason": "disabled_by_config"})
        interaction_tool_statuses["oddt"] = {
            "status": "not_evaluable", "reason": "disabled_by_config", "enabled": False,
        }
    else:
        _progress("starting isolated ODDT worker")
        result, error = _run_isolated_native_step(
            "oddt",
            {
                "protein_file": str(protein_file),
                "ligand_file": str(ligand_file),
                "ligand_resname": ligand_resname,
                "strict": bool(oddt_params.get("strict", True)),
                "splif_depth": int(oddt_params.get("splif_depth", 1)),
                "splif_size": int(oddt_params.get("splif_size", 4096)),
                "splif_distance_cutoff": float(
                    oddt_params.get("splif_distance_cutoff", 4.5)
                ),
            },
            native_step_timeout,
            runtime=profiler_runtimes["oddt"],
        )
        supplemental_results["oddt"] = result
        if result and result.get("success"):
            steps_completed.append("oddt")
            parameters["oddt"] = {
                "role": "supplemental_residue_ifp_and_relative_splif",
                "backend": result.get("backend"),
                "strict": bool(oddt_params.get("strict", True)),
                "splif_depth": int(oddt_params.get("splif_depth", 1)),
                "splif_size": int(oddt_params.get("splif_size", 4096)),
                "splif_distance_cutoff": float(
                    oddt_params.get("splif_distance_cutoff", 4.5)
                ),
            }
            interaction_tool_statuses["oddt"] = {
                "status": "detected" if result.get("num_interactions_total") else "not_detected",
                "enabled": True,
                "tool_version": result.get("tool_version"),
                "interaction_count": int(result.get("num_interactions_total") or 0),
                "quality_flags": result.get("quality_flags") or [],
            }
        else:
            reason = error or (result or {}).get("error") or "unknown"
            steps_skipped.append({"step": "oddt", "reason": reason})
            interaction_tool_statuses["oddt"] = {
                "status": "tool_failed", "reason": reason, "enabled": True,
            }
            supplemental_profiler_failures.append({"tool": "oddt", "reason": reason})

    for tool_id, status in interaction_tool_statuses.items():
        status["critical_for_hard_gate"] = tool_id == "plip"

    interaction_observations = normalize_interaction_observations(
        plip_result,
        prolif_result,
        default_chain=chain,
        tool_statuses=interaction_tool_statuses,
        input_evidence=input_evidence,
        supplemental_results=supplemental_results,
    )

    # ===== Step 4: Blueprint Generation =====
    if not switches["blueprint"]:
        _progress("blueprint disabled_by_config")
        interaction_blueprint = {
            "required": [],
            "desired": [],
            "optional": [],
            "generation_policy": {"note": "disabled_by_config"},
            "validation": {
                "applicable": False,
                "not_applicable": True,
                "valid_for_hard_gating": None,
                "reason": "disabled_by_config",
                "warnings": [],
            },
        }
        steps_skipped.append({"step": "blueprint", "reason": "disabled_by_config"})
    else:
        _progress("generating interaction blueprint")
        interaction_blueprint = generate_blueprint(
            plip_result=plip_result,
            prolif_result=prolif_result,
            # Never fall back to the candidate's residue union: that geometry
            # envelope includes residues seen by only one detector.
            consensus_residues=(
                (consensus or {}).get("interaction_consensus_residue_records") or []
            ),
            known_site_residues=known_site_residues,
            known_required_interactions=known_required_interactions,
            chain=chain,
            blueprint_applicable=bool(ligand_file or known_required_interactions),
            interaction_observations=interaction_observations,
            blueprint_policy=blueprint_policy,
            known_site_provenance=known_site_provenance,
            pocket_support_records=(consensus or {}).get("interaction_pocket_support") or [],
        )
        if interaction_blueprint.get("required") or interaction_blueprint.get("desired") or interaction_blueprint.get("optional"):
            steps_completed.append("blueprint")
        elif plip_result or prolif_result or any(supplemental_results.values()):
            steps_skipped.append({"step": "blueprint", "reason": "no classifiable interactions found"})
        else:
            steps_skipped.append({"step": "blueprint", "reason": "interaction profilers unavailable"})

    # ===== Pocket metrics =====
    _progress("computing pocket metrics")
    pocket_features = {}
    metric_protein = (
        Path(input_data.get("protein_noH_file") or protein_file)
        if provided_box is not None else detector_protein
    )
    if consensus and consensus.get("residues"):
        try:
            # Only cavity-aware tools provide a physical pocket volume. A convex
            # hull over residue atoms is retained separately as an envelope metric.
            selected_fpocket = consensus.get("selected_fpocket") or {}
            selected_pykvfinder = consensus.get("selected_pykvfinder") or {}
            selected_deeppocket = consensus.get("selected_deeppocket") or {}
            if selected_fpocket.get("volume") is not None:
                volume = selected_fpocket.get("volume")
                volume_method = "fpocket_selected_pocket"
            elif selected_pykvfinder.get("volume") is not None:
                volume = selected_pykvfinder.get("volume")
                volume_method = "pykvfinder_selected_pocket"
            elif selected_deeppocket.get("volume") is not None:
                # DeepPocket starts with an fpocket proposal; its reported
                # volume is therefore seed-cavity provenance, not CNN-mask
                # voxel volume.
                volume = selected_deeppocket.get("volume")
                volume_method = "deeppocket_fpocket_seed_volume"
            else:
                volume = None
                volume_method = "unavailable"
            metric_residues = consensus.get("residue_records") or consensus["residues"]
            residue_envelope_volume = compute_pocket_volume(metric_protein, metric_residues, chain)
            burial = compute_burial_ratio(metric_protein, metric_residues, chain)
            subpockets = classify_subpockets(metric_protein, consensus["residues"], chain)
            flex_zones = compute_flexibility_zones(metric_protein, metric_residues, chain)
            pocket_features = {
                "volume": volume,
                "volume_unit": "angstrom3",
                "volume_method": volume_method,
                "residue_envelope_volume": residue_envelope_volume,
                "residue_envelope_volume_unit": "angstrom3",
                "burial_fraction_proxy": burial,
                "burial_method": "heavy_atom_neighbor_fraction_within_8A",
                "key_subpockets": subpockets,
                "flexibility_zones": flex_zones,
            }
            steps_completed.append("pocket_metrics")
        except Exception as e:
            reason = str(e)
            steps_skipped.append({"step": "pocket_metrics", "reason": reason})
            enabled_tool_failures.append({"tool": "pocket_metrics", "reason": reason})

    # ===== Validation =====
    docking_box = provided_box if provided_box is not None else _make_docking_box(
        consensus, detector_protein, chain, box_padding,
        site_selection_strategy=site_selection_strategy,
        ligand_file=ligand_file, min_box_size=min_box_size,
    )
    detector_results = {
        "autosite": autosite_result,
        "fpocket": fpocket_result,
        "p2rank": p2rank_result,
        "pykvfinder": pykvfinder_result,
        "deeppocket": (
            deeppocket_result
            if deeppocket_result and deeppocket_result.get("success") is True
            else None
        ),
    }
    validation = _validate_pocket(
        site_selection_strategy,
        consensus,
        docking_box,
        ligand_site_result,
        detector_results,
        known_site_residues,
        chain,
    )
    if provided_box is not None:
        validation["input_evidence"] = input_evidence["docking_box"]
    pocket_definition = _make_pocket_definition(
        site_selection_strategy, consensus, docking_box, validation,
    )
    docking_box_valid = bool(validation.get("box_valid"))
    blueprint_validation = interaction_blueprint.get("validation") or {}
    blueprint_contract_valid = bool(
        not switches["blueprint"]
        or blueprint_validation.get("applicable") is False
        or blueprint_validation.get("blueprint_valid", True) is True
    )
    # Geometry is checked independently so a failed run can retain a diagnostic
    # proposal.  Canonical artifacts are published only when the full explicitly
    # selected tool contract succeeds.
    geometry_valid = validation.get("passed") is True and docking_box_valid
    blocking_profiler_failures = [
        item for item in profiler_failures if item.get("tool") != "prolif"
    ]
    contract_valid = (
        geometry_valid
        and blueprint_contract_valid
        and not enabled_tool_failures
        and not blocking_profiler_failures
        and not supplemental_profiler_failures
    )
    contract_error = None
    if not contract_valid:
        if enabled_tool_failures:
            failures = "; ".join(
                f"{item['tool']}: {item['reason']}"
                for item in enabled_tool_failures
            )
            contract_error = (
                "An explicitly enabled pocket-analysis tool failed: "
                f"{failures}"
            )
        elif blocking_profiler_failures:
            failures = "; ".join(
                f"{item['tool']}: {item['reason']}" for item in blocking_profiler_failures
            )
            contract_error = (
                "An enabled interaction profiler failed: "
                f"{failures}"
            )
        elif supplemental_profiler_failures:
            failures = "; ".join(
                f"{item['tool']}: {item['reason']}"
                for item in supplemental_profiler_failures
            )
            contract_error = (
                "An enabled supplemental interaction profiler failed: "
                f"{failures}"
            )
        elif not blueprint_contract_valid:
            detail = "; ".join(
                blueprint_validation.get("curated_constraint_errors") or []
            )
            contract_error = (
                "Interaction Blueprint failed its evidence contract"
                + (f": {detail}" if detail else "")
            )
        else:
            contract_error = (
                "pocket-definition contract failed: "
                f"{validation.get('reason') or 'strategy validation did not pass'}"
            )
    if not contract_valid:
        # A forced rerun may reuse an output directory from an earlier success.
        # Remove only the two canonical downstream contracts so stale geometry
        # cannot be mistaken for output from this failed attempt.
        for stale_name in ("pocket_definition.json", "docking_box.json"):
            if (output_dir / stale_name).resolve() == provided_box_file:
                continue
            (output_dir / stale_name).unlink(missing_ok=True)

    # ===== Write output files =====
    _progress("writing output files")
    report = {
        "success": contract_valid,
        "error": contract_error,
        "pocket_detection": {
            "ligand_site": ligand_site_result,
            "fpocket": fpocket_result,
            "pykvfinder": pykvfinder_result,
            "p2rank": p2rank_result,
            "deeppocket": deeppocket_result,
            "autosite": autosite_result,
            "memory_plan": pocket_plan,
            "pocketminer": pocketminer_result,
            "provisional_cryptic_hypotheses": (
                (pocketminer_result or {}).get("provisional_cryptic_hypotheses") or []
            ),
            "consensus_pocket": consensus,
        },
        "interaction_analysis": {
            "plip": plip_result,
            "prolif": prolif_result,
            "arpeggio": supplemental_results["arpeggio"],
            "luna": supplemental_results["luna"],
            "oddt": supplemental_results["oddt"],
            "tool_statuses": interaction_tool_statuses,
            "profiler_failures": profiler_failures,
            "supplemental_profiler_failures": supplemental_profiler_failures,
        },
        "enabled_tool_failures": enabled_tool_failures,
        "interaction_observations": interaction_observations,
        "interaction_blueprint": interaction_blueprint,
        "input_evidence": input_evidence,
        "pocket_features": pocket_features,
        # A failed run may retain the computed proposal for diagnosis, but never
        # publishes canonical downstream artifacts from a degraded tool set.
        "pocket_definition": pocket_definition if contract_valid else None,
        "docking_box": docking_box if contract_valid else None,
        "proposed_pocket_definition": (
            pocket_definition if not contract_valid else None
        ),
        "proposed_docking_box": docking_box if not contract_valid else None,
        "validation": validation,
        "steps_completed": steps_completed,
        "steps_skipped": steps_skipped,
        "parameters": parameters,
        "interaction_policy": {
            "enabled_tool_failure_policy": "hard_fail",
            "blueprint_policy": blueprint_policy,
            "known_site_provenance": known_site_provenance,
            "prolif_input_strategy": (
                (parameters.get("prolif") or {}).get("input_strategy")
            ),
        },
    }

    report_file = output_dir / "pocket_analysis_report.json"
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    bp_file = output_dir / "interaction_blueprint.json"
    with open(bp_file, "w") as f:
        json.dump(interaction_blueprint, f, indent=2, ensure_ascii=False)

    observations_file = output_dir / "interaction_observations.json"
    with open(observations_file, "w") as f:
        json.dump(interaction_observations, f, indent=2, ensure_ascii=False, default=str)

    pocket_definition_file = None
    if contract_valid:
        pocket_definition_file = output_dir / "pocket_definition.json"
        with open(pocket_definition_file, "w") as f:
            json.dump(pocket_definition, f, indent=2, ensure_ascii=False, default=str)

    docking_box_file = provided_box_file if contract_valid else None
    if contract_valid and provided_box is None:
        docking_box_file = output_dir / "docking_box.json"
        with open(docking_box_file, "w") as f:
            json.dump(docking_box, f, indent=2, ensure_ascii=False)

    if plip_result and plip_result.get("success"):
        plip_file = output_dir / "plip_interactions.json"
        with open(plip_file, "w") as f:
            json.dump(plip_result, f, indent=2, ensure_ascii=False)

    if prolif_result and not prolif_result.get("validation_failed"):
        ifp_file = output_dir / "interaction_fingerprint.json"
        with open(ifp_file, "w") as f:
            json.dump(prolif_result, f, indent=2, ensure_ascii=False)

    supplemental_output_files = {
        "arpeggio": output_dir / "arpeggio_interactions.json",
        "luna": output_dir / "luna_interaction_fingerprint.json",
        "oddt": output_dir / "oddt_interaction_fingerprint.json",
    }
    for tool_id, output_file in supplemental_output_files.items():
        tool_result = supplemental_results.get(tool_id)
        if tool_result and tool_result.get("success"):
            with open(output_file, "w") as f:
                json.dump(tool_result, f, indent=2, ensure_ascii=False, default=str)
        else:
            output_file.unlink(missing_ok=True)

    if pocket_features:
        pf_file = output_dir / "pocket_features.json"
        with open(pf_file, "w") as f:
            json.dump(pocket_features, f, indent=2, ensure_ascii=False)

    result = {
        "success": contract_valid,
        "error": contract_error,
        "pocket_detection": report["pocket_detection"],
        "interaction_analysis": report["interaction_analysis"],
        "enabled_tool_failures": enabled_tool_failures,
        "interaction_observations": interaction_observations,
        "interaction_blueprint": interaction_blueprint,
        "pocket_features": pocket_features,
        "pocket_definition": pocket_definition if contract_valid else None,
        "docking_box": docking_box if contract_valid else None,
        "proposed_pocket_definition": (
            pocket_definition if not contract_valid else None
        ),
        "proposed_docking_box": docking_box if not contract_valid else None,
        "validation": validation,
        "steps_completed": steps_completed,
        "steps_skipped": steps_skipped,
        "parameters": parameters,
        "files": {
            "report": str(report_file),
            "blueprint": str(bp_file),
            "interaction_observations": str(observations_file),
            "plip": str(output_dir / "plip_interactions.json") if (plip_result and plip_result.get("success")) else None,
            "fingerprint": str(output_dir / "interaction_fingerprint.json") if (prolif_result and not prolif_result.get("validation_failed")) else None,
            "arpeggio": (
                str(supplemental_output_files["arpeggio"])
                if supplemental_results.get("arpeggio")
                and supplemental_results["arpeggio"].get("success")
                else None
            ),
            "luna": (
                str(supplemental_output_files["luna"])
                if supplemental_results.get("luna")
                and supplemental_results["luna"].get("success")
                else None
            ),
            "oddt": (
                str(supplemental_output_files["oddt"])
                if supplemental_results.get("oddt")
                and supplemental_results["oddt"].get("success")
                else None
            ),
            "pocket_features": str(output_dir / "pocket_features.json") if pocket_features else None,
            "pocketminer": (
                str(pocketminer_file)
                if switches["pocketminer"] and pocketminer_file.is_file()
                else None
            ),
            "deeppocket": (
                str(deeppocket_file)
                if switches["deeppocket"] and deeppocket_file.is_file()
                else None
            ),
            "deeppocket_raw": (
                str(output_dir / "deeppocket_output")
                if switches["deeppocket"]
                and (output_dir / "deeppocket_output").is_dir()
                else None
            ),
            "pocket_definition": str(pocket_definition_file) if pocket_definition_file else None,
            "docking_box": str(docking_box_file) if docking_box_file else None,
        },
    }
    return result


def main():
    from io_utils import parse_skill_args, print_skill_result
    input_data, output_dir = parse_skill_args("Pocket Analysis Skill")
    result = run(input_data, output_dir)
    print_skill_result(result)
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()

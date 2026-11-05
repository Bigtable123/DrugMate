"""Diagnosis engine: generate human-readable diagnosis, recommendation, and modification hints."""
from __future__ import annotations

import sys
from pathlib import Path

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from constants import (
    DEFAULT_POSECHECK_THRESHOLDS,
    POSEBUSTERS_VIOLATION_THRESHOLD,
    POCKET_OCCUPANCY_LOW,
    POCKET_OCCUPANCY_HIGH,
)


def generate_diagnosis(
    entry: dict,
    subpockets: list[dict] | None = None,
    posecheck_thresholds: dict | None = None,
) -> tuple[list[str], str, str]:
    """Generate diagnosis list, recommendation, and modification hint for a pose.

    Args:
        entry: Combined evaluation entry (posebusters, posecheck, plip, blueprint_match, pocket_contact_fraction)
        subpockets: Subpocket definitions for occupancy-based hints

    Returns:
        (diagnosis_list, recommendation, modification_hint)
    """
    thresholds = posecheck_thresholds or DEFAULT_POSECHECK_THRESHOLDS
    diagnosis = []
    issues = []

    if entry.get("verdict") == "not_evaluated":
        diagnosis.append("EVIDENCE_MISSING: pose could not be evaluated by the required tools")
        issues.append("evidence_missing")

    # PoseBusters
    pb = entry.get("posebusters", {})
    n_violations = pb.get("num_violations")
    if n_violations is not None and n_violations > POSEBUSTERS_VIOLATION_THRESHOLD:
        diagnosis.append(
            f"PHYSICAL_VALIDITY_FAIL: {n_violations} PoseBusters "
            "violation(s) (bond lengths/angles/clashes)"
        )
        issues.append("geometric_violations")

    # PoseCheck
    pc = entry.get("posecheck", {})
    strain = pc.get("strain_energy_kcal")
    clashes = pc.get("steric_clash_count")
    all_atom_clashes = pc.get("steric_clash_count_all_atom")
    all_atom_note = (
        f"; all-atom PoseCheck: {all_atom_clashes}"
        if all_atom_clashes is not None
        else ""
    )
    if clashes is not None and clashes > thresholds["clash_fail_count"]:
        diagnosis.append(
            f"CLASH_FAIL: {clashes} heavy-atom steric clash(es) with receptor "
            f"(fail threshold: >{thresholds['clash_fail_count']}{all_atom_note})"
        )
        issues.append("steric_clashes")
    elif clashes is not None and clashes > thresholds["clash_warn_count"]:
        diagnosis.append(
            f"CLASH_WARN: {clashes} heavy-atom steric clash(es) with receptor "
            f"(native-calibrated fail threshold: >{thresholds['clash_fail_count']}{all_atom_note})"
        )
        issues.append("steric_clash_warning")
    elif clashes is not None:
        diagnosis.append(
            "GOOD: No heavy-atom PoseCheck steric clashes above warning threshold"
            f" ({all_atom_clashes} all-atom)" if all_atom_clashes is not None
            else "GOOD: No heavy-atom PoseCheck steric clashes above warning threshold"
        )

    if strain is not None:
        if strain > thresholds["strain_fail_kcal"]:
            diagnosis.append(
                f"HIGH_STRAIN: {strain:.1f} kcal/mol "
                f"(fail threshold: {thresholds['strain_fail_kcal']:.1f})"
            )
            issues.append("high_strain")
        elif strain > thresholds["strain_warn_kcal"]:
            diagnosis.append(
                f"MODERATE_STRAIN: {strain:.1f} kcal/mol "
                f"(warn threshold: {thresholds['strain_warn_kcal']:.1f})"
            )
            issues.append("strain_warning")
        else:
            diagnosis.append(f"GOOD: Strain energy {strain:.1f} kcal/mol within threshold")

    # Blueprint match
    bp = entry.get("blueprint_match", {})
    required_hits = bp.get("required_hits")
    required_total = bp.get("required_total")
    missed = bp.get("required_missed", [])
    unresolved = bp.get("required_unresolved", [])
    bp_verdict = bp.get("verdict")
    if bp_verdict == "not_applicable" and bp:
        diagnosis.append("UNCONSTRAINED_BLUEPRINT: no required interactions defined")
    elif bp_verdict == "unresolved":
        diagnosis.append(
            f"UNRESOLVED_BLUEPRINT: {len(unresolved) or required_total or 0} required "
            "constraint(s) cannot be evaluated with the available interaction evidence"
        )
        issues.append("blueprint_unresolved")
    elif bp_verdict == "fail":
        diagnosis.append(
            f"MISSING: {required_hits or 0}/{required_total or 0} required constraints met "
            "under exact group logic"
        )
        issues.append("blueprint_mismatch")
    elif bp_verdict == "pass":
        diagnosis.append(
            f"GOOD: all required Blueprint groups satisfied "
            f"({required_hits or 0}/{required_total or 0} individual constraints matched)"
        )

    if missed:
        for m in missed[:3]:
            res = m.get("residue", "?")
            itype = m.get("interaction_type", "?")
            diagnosis.append(f"MISSING: Required {itype} with {res}")

    desired_score = bp.get("desired_weighted_score")
    if desired_score is not None:
        diagnosis.append(
            f"DESIRED_SCORE: {float(desired_score):.1%} weighted desired interactions satisfied"
        )

    # Pocket occupancy
    occ = entry.get("pocket_contact_fraction", entry.get("pocket_occupancy", {}))
    low_occ_regions = []
    for sp_id, val in occ.items():
        if sp_id == "overall":
            continue
        if val is None:
            continue
        if val < POCKET_OCCUPANCY_LOW:
            diagnosis.append(f"LOW_OCCUPANCY: {sp_id} ({int(val*100)}%)")
            low_occ_regions.append(sp_id)
        elif val >= POCKET_OCCUPANCY_HIGH:
            diagnosis.append(f"GOOD: {sp_id} well occupied ({int(val*100)}%)")

    if not diagnosis:
        diagnosis.append("GOOD: No issues detected")

    # Recommendation
    verdict = entry.get("verdict", "pass")
    recommendation = _compute_recommendation(
        verdict, issues, clashes, strain, n_violations, thresholds
    )

    # Modification hint
    modification_hint = _build_modification_hint(
        issues, missed, low_occ_regions, subpockets, strain, thresholds
    )

    return diagnosis, recommendation, modification_hint


def _compute_recommendation(
    verdict: str,
    issues: list[str],
    clashes: int | None,
    strain: float | None,
    n_violations: int | None,
    posecheck_thresholds: dict,
) -> str:
    """Determine keep/modify/reject based on verdict and issue severity.

    ``n_violations`` is ``None`` when the producer never counted violations (see
    the note in :func:`generate_diagnosis` and structure-evaluator/run.py where
    the key is deliberately left absent).  "Never counted" is absence of
    evidence, not a severe violation, so it must not be compared as a number.
    """
    if verdict == "fail":
        severe = (
            (clashes is not None and clashes > posecheck_thresholds["clash_fail_count"])
            or (strain is not None and strain > posecheck_thresholds["strain_fail_kcal"])
            or (
                n_violations is not None
                and n_violations > POSEBUSTERS_VIOLATION_THRESHOLD
            )
        )
        if severe:
            return "reject"
        return "modify"
    elif verdict == "warn":
        return "modify"
    elif verdict == "not_evaluated":
        return "not_evaluated"
    return "keep"


def _build_modification_hint(
    issues: list[str],
    missed_required: list[dict],
    low_occ_regions: list[str],
    subpockets: list[dict] | None,
    strain: float | None,
    posecheck_thresholds: dict,
) -> str:
    """Build actionable modification hint from identified issues."""
    hints = []

    if missed_required:
        for m in missed_required[:2]:
            res = m.get("residue", "?")
            itype = m.get("interaction_type", "?")
            if "H-bond" in itype:
                hints.append(f"Add H-bond donor/acceptor group toward {res}")
            elif "Hydrophobic" in itype:
                hints.append(f"Extend hydrophobic moiety toward {res}")
            elif "Pi" in itype:
                hints.append(f"Add aromatic ring in proximity to {res}")
            elif "Salt" in itype:
                hints.append(f"Introduce charged group near {res}")
            else:
                hints.append(f"Add functional group for {itype} with {res}")

    if low_occ_regions and subpockets:
        sp_map = {sp["id"]: sp for sp in subpockets}
        for region in low_occ_regions[:2]:
            sp = sp_map.get(region, {})
            residues = sp.get("residues", [])[:3]
            prop = sp.get("property", "")
            if residues:
                hints.append(
                    f"Extend scaffold into {prop} region near {', '.join(residues)}"
                )

    if strain is not None and strain > posecheck_thresholds["strain_warn_kcal"]:
        hints.append("Reduce ring strain or simplify core scaffold")

    if "steric_clashes" in issues or "steric_clash_warning" in issues:
        hints.append("Trim substituent causing receptor clash or adjust binding pose orientation")

    return "; ".join(hints) if hints else ""

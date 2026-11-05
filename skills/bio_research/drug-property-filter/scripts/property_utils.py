"""Fast RDKit drug-property, synthetic-accessibility, and alert utilities."""
from __future__ import annotations

import os
import sys
from functools import lru_cache
from typing import Any, Iterable


DEFAULT_ALERT_CATALOGS = ("PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH")
DEFAULT_HARD_FAIL_CATALOGS = DEFAULT_ALERT_CATALOGS

# The executable policy turns each alert catalog into a count metric: 0 hits is a
# ``pass`` mark, 1-4 hits a ``flag`` mark, and only >4 hits a ``fail`` mark.
# A single hit therefore no longer ends the molecule here; it is one level-1
# mark that the category fold in shared/metric_scoring.py weighs against the
# rest of the ALERT category.
ALERT_COUNT_HARD_UPPER = 4
# The executable policy -- Lipinski RO5 violations: 0-1 pass, 2-3 flag, >3 fail.
DEFAULT_LIPINSKI_MAX_VIOLATIONS = 3


def check_rdkit_available() -> bool:
    try:
        from rdkit import Chem  # noqa: F401
        return True
    except ImportError:
        return False


def check_sa_score_available() -> bool:
    try:
        from rdkit.Chem import RDConfig

        sa_path = os.path.join(
            RDConfig.RDContribDir, "SA_Score", "sascorer.py"
        )
        return os.path.isfile(sa_path)
    except Exception:
        return False


def check_pains_filter_available() -> bool:
    """Backward-compatible availability probe for all structural catalogs."""
    try:
        from rdkit.Chem.FilterCatalog import (  # noqa: F401
            FilterCatalog,
            FilterCatalogParams,
        )
        return True
    except ImportError:
        return False


def _catalog_enum(name: str):
    from rdkit.Chem.FilterCatalog import FilterCatalogParams

    normalized = str(name).strip().upper()
    try:
        return getattr(FilterCatalogParams.FilterCatalogs, normalized)
    except AttributeError as exc:
        raise ValueError(f"unknown RDKit FilterCatalog: {name!r}") from exc


@lru_cache(maxsize=None)
def _get_catalog(name: str):
    """Create one catalog at a time so every hit retains its catalog source."""
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    normalized = str(name).strip().upper()
    params = FilterCatalogParams()
    params.AddCatalog(_catalog_enum(normalized))
    return FilterCatalog(params)


def _match_atom_indices(entry: Any, mol: Any) -> list[int]:
    """Extract alert atoms defensively across supported RDKit releases."""
    atom_indices: set[int] = set()
    try:
        filter_matches = entry.GetFilterMatches(mol)
    except Exception:
        filter_matches = []
    for filter_match in filter_matches or []:
        for pair in getattr(filter_match, "atomPairs", ()) or ():
            target = getattr(pair, "target", None)
            if target is None and isinstance(pair, (tuple, list)) and len(pair) > 1:
                target = pair[1]
            try:
                atom_indices.add(int(target))
            except (TypeError, ValueError):
                continue
    return sorted(atom_indices)


def check_structural_alerts(
    mol: Any,
    catalogs: Iterable[str] | None = None,
) -> dict:
    """Return every PAINS/Brenk/NIH hit with catalog and matched atom indices."""
    requested = tuple(
        dict.fromkeys(
            str(name).strip().upper()
            for name in (catalogs or DEFAULT_ALERT_CATALOGS)
            if str(name).strip()
        )
    )
    if not check_pains_filter_available():
        return {
            "available": False,
            "catalogs": list(requested),
            "has_alert": None,
            "matches": [],
            "error": "RDKit FilterCatalog unavailable",
        }

    matches: list[dict] = []
    catalog_errors: list[dict] = []
    for catalog_name in requested:
        try:
            entries = list(_get_catalog(catalog_name).GetMatches(mol))
        except Exception as exc:
            catalog_errors.append({
                "catalog": catalog_name,
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        for entry in entries:
            matches.append({
                "catalog": catalog_name,
                "description": str(entry.GetDescription()),
                "atom_indices": _match_atom_indices(entry, mol),
            })

    return {
        "available": not catalog_errors,
        "catalogs": list(requested),
        "has_alert": bool(matches),
        "match_count": len(matches),
        "matches": matches,
        "catalog_errors": catalog_errors,
    }


def compute_drug_properties(mol: Any) -> dict:
    """Compute deterministic first-pass descriptors and Ertl SA score."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, rdMolDescriptors

    mol_no_h = Chem.RemoveHs(Chem.Mol(mol))
    mw = float(Descriptors.ExactMolWt(mol_no_h))
    logp = float(Descriptors.MolLogP(mol_no_h))
    tpsa = float(Descriptors.TPSA(mol_no_h))
    fraction_csp3 = float(rdMolDescriptors.CalcFractionCSP3(mol_no_h))
    qed = float(QED.qed(mol_no_h))
    props = {
        # Machine-readable descriptor fields retain full precision for gating.
        # Explicit display copies preserve the compact legacy presentation.
        "MW": mw,
        "MW_display": round(mw, 2),
        "LogP": logp,
        "LogP_display": round(logp, 2),
        "HBD": int(rdMolDescriptors.CalcNumHBD(mol_no_h)),
        "HBA": int(rdMolDescriptors.CalcNumHBA(mol_no_h)),
        "TPSA": tpsa,
        "TPSA_display": round(tpsa, 2),
        "RotBonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol_no_h)),
        "RingCount": int(rdMolDescriptors.CalcNumRings(mol_no_h)),
        "AromaticRingCount": int(rdMolDescriptors.CalcNumAromaticRings(mol_no_h)),
        "FractionCSP3": fraction_csp3,
        "FractionCSP3_display": round(fraction_csp3, 3),
        "FormalCharge": int(Chem.GetFormalCharge(mol_no_h)),
        "HeavyAtomCount": int(rdMolDescriptors.CalcNumHeavyAtoms(mol_no_h)),
        "ChiralCenters": len(
            Chem.FindMolChiralCenters(mol_no_h, includeUnassigned=True)
        ),
        "QED": qed,
        "QED_display": round(qed, 3),
    }
    if check_sa_score_available():
        try:
            sa_score = float(_compute_sa_score(mol_no_h))
            props["SA_score"] = sa_score
            props["SA_score_display"] = round(sa_score, 2)
        except Exception:
            props["SA_score"] = None
            props["SA_score_display"] = None
    else:
        props["SA_score"] = None
        props["SA_score_display"] = None
    return props


def check_lipinski(props: dict) -> dict:
    """Check the four classic Rule-of-Five bounds."""
    details = {
        "MW_ok": float(props.get("MW", 0)) <= 500,
        "LogP_ok": float(props.get("LogP", 0)) <= 5.0,
        "HBD_ok": int(props.get("HBD", 0)) <= 5,
        "HBA_ok": int(props.get("HBA", 0)) <= 10,
    }
    violations = sum(1 for value in details.values() if not value)
    return {
        "violations": violations,
        # Compatibility alias used by older reports.
        "num_violations": violations,
        "details": details,
    }


def check_pains(mol: Any) -> dict:
    """Backward-compatible PAINS projection over the richer alert contract."""
    alerts = check_structural_alerts(
        mol, catalogs=("PAINS_A", "PAINS_B", "PAINS_C")
    )
    if alerts.get("has_alert") is None:
        return {"is_pains": None, "matched_filters": [], "matches": []}
    matches = alerts.get("matches") or []
    return {
        "is_pains": bool(matches),
        "matched_filters": [item["description"] for item in matches],
        "matches": matches,
    }


def compute_verdict(
    props: dict,
    lipinski: dict,
    pains: dict,
    thresholds: dict,
    enabled_checks: set | None = None,
    *,
    structural_alerts: dict | None = None,
    hard_fail_catalogs: Iterable[str] | None = None,
    is_covalent: bool = False,
) -> tuple[str, list[str], list[str], dict]:
    """Compute a transparent RDKit first-pass verdict.

    This is the level-1 marking pass only: it returns ``pass``/``flag``/``fail``
    for the RDKit screen as a whole and never decides keep/modify/reject. Under
    The executable policy that decision belongs to the per-category fold, which counts
    the ``fail`` marks inside each parent metric category once the whole step is
    computed.

    The fourth return value separates the two reasons a verdict can be non-pass,
    which the human-readable ``flags`` list conflates: ``molecular_findings``
    names the marks this molecule actually earned, while ``coverage_gaps`` names
    the checks that could not be run at all. Callers must not read a coverage gap
    as a molecular finding -- the executable policy forbids penalising a molecule for a tool
    or coverage problem, so a missing sascorer must not be reported the same way
    as a high SA score.

    Alert catalogs hard-fail only above ``ALERT_COUNT_HARD_UPPER`` hits, so a
    single PAINS/Brenk/NIH match is a ``flag`` mark that keeps the molecule in
    play with an audit note. Declaring a workflow covalent adds review guidance
    but does not rescue an alert: rescue must be an explicit, auditable policy
    decision (for example by removing a catalog from hard_fail_catalogs).
    """
    checks = (
        enabled_checks
        if enabled_checks is not None
        else {"lipinski", "qed", "sa_score", "pains"}
    )
    if not checks:
        # Every RDKit check is off, so nothing below can run. An empty flag list
        # is then absence of evidence, not a clean molecule: say so instead of
        # returning "pass" and letting downstream gates read it as a pass.
        return (
            "not_evaluated",
            ["No RDKit property check enabled"],
            [
                "Enable at least one of lipinski/qed/sa_score/pains, "
                "or assess drug-likeness with another tool"
            ],
            {"molecular_findings": [], "coverage_gaps": ["no_check_enabled"]},
        )
    flags: list[str] = []
    suggestions: list[str] = []
    molecular_findings: list[str] = []
    coverage_gaps: list[str] = []
    hard_failure = False

    qed_thresh = float(thresholds.get("qed_threshold", 0.45))
    sa_thresh = float(thresholds.get("sa_threshold", 6.0))
    max_lip = int(
        thresholds.get(
            "lipinski_max_violations", DEFAULT_LIPINSKI_MAX_VIOLATIONS
        )
    )

    if "pains" in checks:
        alerts = structural_alerts
        if alerts is None:
            legacy_matches = [
                {
                    "catalog": "PAINS",
                    "description": description,
                    "atom_indices": [],
                }
                for description in (pains.get("matched_filters") or [])
            ]
            alerts = {
                "available": pains.get("is_pains") is not None,
                "has_alert": bool(legacy_matches),
                "matches": legacy_matches,
            }
        if not alerts.get("available", True):
            coverage_gaps.append("alert_screening_incomplete")
            flags.append("Structural-alert screening incomplete")
            suggestions.append(
                "Structural-alert screening incomplete; review before nomination"
            )
        matches = alerts.get("matches") or []
        if matches:
            catalogs = sorted({str(item.get("catalog") or "unknown") for item in matches})
            molecular_findings.append("structural_alerts_present")
            flags.append(
                f"Structural alerts ({len(matches)}): {', '.join(catalogs)}"
            )
            suggestions.append(
                "Review alert substructures, assay-interference risk, and chemical stability"
            )
            hard_catalogs = {
                str(name).strip().upper()
                for name in (
                    hard_fail_catalogs
                    if hard_fail_catalogs is not None
                    else DEFAULT_HARD_FAIL_CATALOGS
                )
            }
            catalog_counts: dict[str, int] = {}
            for item in matches:
                catalog = str(item.get("catalog") or "").strip().upper()
                if catalog:
                    catalog_counts[catalog] = catalog_counts.get(catalog, 0) + 1
            if _alert_counts_hard_fail(catalog_counts, hard_catalogs):
                hard_failure = True
            if is_covalent:
                suggestions.append(
                    "Covalent route: confirm the exact alert/warhead manually; "
                    "covalent status alone does not rescue a hard alert"
                )

    if "lipinski" in checks and int(lipinski.get("violations", 0)) > max_lip:
        hard_failure = True
        molecular_findings.append("lipinski_violations_above_max")
        flags.append(f"Lipinski violations: {lipinski['violations']}")
        failing = [key for key, value in lipinski.get("details", {}).items() if not value]
        if "MW_ok" in failing:
            suggestions.append("Reduce molecular weight below 500")
        if "LogP_ok" in failing:
            suggestions.append("Reduce lipophilicity")
        if "HBD_ok" in failing:
            suggestions.append("Reduce H-bond donors")
        if "HBA_ok" in failing:
            suggestions.append("Reduce H-bond acceptors")

    if "qed" in checks:
        qed = props.get("QED")
        if qed is not None and float(qed) < qed_thresh:
            molecular_findings.append("qed_below_threshold")
            flags.append(f"Low QED ({float(qed):.2f})")
            suggestions.append("Improve the balance of size, polarity, and lipophilicity")

    if "sa_score" in checks:
        sa_score = props.get("SA_score")
        if sa_score is None:
            coverage_gaps.append("sa_score_unavailable")
            flags.append("SA score unavailable")
            suggestions.append("SA score unavailable; retain for review only")
        elif float(sa_score) > sa_thresh:
            molecular_findings.append("sa_score_above_threshold")
            flags.append(f"High SA score ({float(sa_score):.1f})")
            suggestions.append("Simplify stereochemistry, rings, or uncommon motifs")

    evidence = {
        "molecular_findings": molecular_findings,
        "coverage_gaps": coverage_gaps,
    }
    if hard_failure:
        return "fail", flags, suggestions, evidence
    if flags:
        return "flag", flags, suggestions, evidence
    return "pass", flags, suggestions, evidence


def _alert_counts_hard_fail(
    catalog_counts: dict[str, int],
    hard_catalogs: set[str],
    *,
    hard_upper: int = ALERT_COUNT_HARD_UPPER,
) -> bool:
    """Report a ``fail`` mark only once a catalog family exceeds its count bound.

    The executable policy scores every alert catalog as a count (0 pass, 1-4 flag,
    >4 fail), so hits below the bound are flags rather than eliminations.
    PAINS_A/B/C are three RDKit catalogs behind the single
    ``pains_alert_count`` metric, so their hits are summed into one family
    count before the bound is applied instead of being judged one catalog at a
    time. ``hard_catalogs`` still acts as the configurable opt-out: a catalog
    removed from it contributes flags only.
    """
    family_counts: dict[str, int] = {}
    for catalog, count in catalog_counts.items():
        if count <= 0:
            continue
        if not _catalog_is_hard_gated(catalog, hard_catalogs):
            continue
        family = (
            "PAINS"
            if catalog == "PAINS" or catalog.startswith("PAINS_")
            else catalog
        )
        family_counts[family] = family_counts.get(family, 0) + int(count)
    return any(count > hard_upper for count in family_counts.values())


def _catalog_is_hard_gated(catalog: str, hard_catalogs: set[str]) -> bool:
    """Honor the opt-out list, treating the PAINS family as one policy unit."""
    if catalog in hard_catalogs:
        return True
    if catalog == "PAINS":
        # Legacy callers report one aggregated "PAINS" catalog, so any
        # configured PAINS_* subset keeps it under the hard bound.
        return any(name.startswith("PAINS_") for name in hard_catalogs)
    if catalog.startswith("PAINS_"):
        return "PAINS" in hard_catalogs
    return False


def _compute_sa_score(mol: Any) -> float:
    from rdkit.Chem import RDConfig

    sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if sa_dir not in sys.path:
        sys.path.insert(0, sa_dir)
    import sascorer

    return float(sascorer.calculateScore(mol))

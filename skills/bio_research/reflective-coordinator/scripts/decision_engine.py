"""Round evidence classification and advisory diagnostic aggregation."""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from constants import POSEBUSTERS_VIOLATION_THRESHOLD
from metric_scoring import record_has_structural_reject, resolve_screening_mode


# ``not_evaluated`` ranks between ``keep`` and ``modify``.  Missing evidence is
# not a molecular finding, so it may never dominate a real ``modify``/``reject``
#, but it may not resolve to a clean ``keep`` either.
_CLASS_SEVERITY = {"keep": 0, "not_evaluated": 1, "modify": 2, "reject": 3}
_CLASS_ALIASES = {
    "pass": "keep",
    "keep": "keep",
    "flag": "modify",
    "warn": "modify",
    "review_only": "modify",
    "modify": "modify",
    "fail": "reject",
    "rejected": "reject",
    "reject": "reject",
    # Recognized rather than dropped.  An unrecognized value used to fall through
    # ``_first_classification`` as ``None``, and ``_marking_classification`` then
    # let ``verdict == "not_evaluated"`` reach the ``or "keep"`` fallback -- a
    # molecule nothing could measure was read as a clean pass.
    "not_evaluated": "not_evaluated",
    "not_evaluable": "not_evaluated",
}


def classify_molecule(
    eval_entry: dict,
    prop_entry: dict | None,
    property_required: bool = False,
    screening_mode: object = None,
) -> str:
    """Classify a molecule as keep/not_evaluated/modify/reject.

    ``screening_mode`` is the run's resolved screening mode (or anything
    :func:`resolve_screening_mode` accepts).  When the mode does not eliminate --
    V1, "marking only" -- a threshold-derived ``reject`` read off an upstream
    report is capped at ``modify``: the operator asked for marks, and a report
    written in an earlier round under V2/V3 must not smuggle an elimination back
    in.  Structural rejects (unreadable structure, unattributable PoseBusters or
    identity failure, an explicit ``rejected`` eligibility) sit outside the
    mode-driven fold and still reject in every mode.  Leaving this ``None``
    keeps the pre-mode behaviour of honouring every reject.
    """
    eliminates = True
    if screening_mode is not None:
        eliminates = bool(resolve_screening_mode(screening_mode).get("eliminates"))
    classifications = [
        _evaluation_classification(eval_entry, eliminates=eliminates),
        _property_classification(
            prop_entry, required=property_required, eliminates=eliminates
        ),
    ]
    return max(classifications, key=_CLASS_SEVERITY.__getitem__)


def _cap_below_reject(classification: str, *, eliminates: bool, structural: bool) -> str:
    """Hold a threshold-derived reject at ``modify`` when the mode only marks."""
    if eliminates or structural or classification != "reject":
        return classification
    return "modify"


def _evaluation_classification(entry: dict, *, eliminates: bool = True) -> str:
    """Map structure-evaluation evidence to the common three-tier gate.

    A structure-evaluator report written under the two-level contract already
    derived ``verdict`` from the per-category fold, so its ``fail`` is a molecule
    decision and may reject.  Records that only carry the old one-level verdict
    are treated as level-1 marks and capped at ``modify``: under the executable policy a
    single threshold breach is not grounds for rejection, and the explicit gate
    fields above are where a rejection has to come from.  The one exception is a
    ``hard`` reason -- unreadable identity or an unattributable PoseBusters
    failure -- which is structural rather than threshold-based.

    ``eliminates`` is ``False`` under the V1 marking mode; see
    :func:`classify_molecule`.  Only ``hard_reasons`` survive it here.
    """
    verdict = str(entry.get("verdict", "unknown")).strip().casefold()
    recommendation = str(entry.get("recommendation", "")).strip().casefold()

    explicit = _first_classification(
        entry,
        *_molecule_gate_fields(entry),
        "gate_classification",
        "classification",
    )
    reasons = entry.get("verdict_reasons")
    hard_reasons = (
        list(reasons.get("hard") or []) if isinstance(reasons, dict) else []
    )
    structural_failure = bool(hard_reasons) or record_has_structural_reject(entry)
    if structural_failure:
        return "reject"
    metric_gate = _first_classification(entry, *_molecule_gate_fields(entry))
    if metric_gate is not None:
        return _cap_below_reject(
            metric_gate, eliminates=eliminates, structural=False,
        )
    two_level = bool(hard_reasons) or isinstance(reasons, dict) or bool(
        entry.get("metric_category_folds") or entry.get("screening_mode")
    )
    if verdict == "fail" or recommendation == "reject":
        legacy = "reject" if (two_level or hard_reasons) else "modify"
    elif recommendation == "modify" or verdict == "warn":
        legacy = "modify"
    elif verdict == "pass":
        legacy = "keep"
    else:
        # Missing/failed evaluation evidence is a protocol uncertainty, not a
        # molecular hard failure.  Keep it in the review tier; the coordinator
        # separately blocks automatic nomination and requests a rerun.
        legacy = "modify"
    combined = max(
        (value for value in (explicit, legacy) if value is not None),
        key=_CLASS_SEVERITY.__getitem__,
    )
    return _cap_below_reject(
        combined, eliminates=eliminates, structural=structural_failure
    )


def _rdkit_verdict_is_molecular(entry: dict) -> bool:
    """Did the Step 4 RDKit screen mark the molecule, or only fail to look at it?

    ``True`` also when the answer is unknown: reports written before
    ``compute_verdict`` published the split keep their previous reading.
    """
    if isinstance(entry.get("property_tool_error"), dict):
        return False
    evidence = entry.get("rdkit_evidence")
    if not isinstance(evidence, dict):
        return True
    findings = evidence.get("molecular_findings")
    if not isinstance(findings, list):
        return True
    return bool(findings)


def evaluation_data_review(entry: dict) -> dict | None:
    """Describe non-molecular gaps in structure-evaluation evidence."""
    verdict = str(entry.get("verdict") or "unknown").strip().casefold()
    if verdict not in {"not_evaluated", "not_evaluable", "unknown", ""}:
        return None

    evidence = entry.get("evidence_status")
    failed_tools: list[str] = []
    missing_tools: list[str] = []
    if isinstance(evidence, dict):
        for tool, raw_status in evidence.items():
            status = str(raw_status or "").strip().casefold()
            if status in {"tool_failed", "error", "failed"}:
                failed_tools.append(str(tool))
            elif status == "missing":
                missing_tools.append(str(tool))

    reasons = [f"structure_evaluation_{verdict or 'unknown'}"]
    reasons.extend(f"{tool}_tool_failed" for tool in sorted(failed_tools))
    reasons.extend(f"{tool}_missing" for tool in sorted(missing_tools))
    return {
        "reasons": reasons,
        "actions": [
            "Rerun or repair structure evaluation and obtain usable pose-QC "
            "evidence before automatic nomination"
        ],
        "failed_tools": sorted(failed_tools),
        "missing_tools": sorted(missing_tools),
    }


def _property_classification(
    entry: dict | None, *, required: bool, eliminates: bool = True
) -> str:
    """Map property output to the common gate, preferring new gate fields.

    ``eliminates`` is ``False`` under the V1 marking mode; see
    :func:`classify_molecule`.  Explicit structural markers (unreadable
    structure, identity failure, or hard evidence failure) survive it.  A bare
    ``nomination_eligibility="rejected"`` is not structural because older
    property reports also used it for metric-fold rejection.
    """
    if entry is None:
        return "modify" if required else "keep"

    # ``nomination_eligibility=rejected`` is used by the property filter for
    # both structural failures and a metric category fold.  Only the explicit
    # structural/evidence signals are hard; a metric rejection is recomputed
    # under the active screening mode and must not survive V1 as an elimination.
    structural = record_has_structural_reject(entry)
    if structural:
        return "reject"
    metric_gate = _first_classification(entry, *_molecule_gate_fields(entry))
    if metric_gate is not None:
        return _cap_below_reject(
            metric_gate, eliminates=eliminates, structural=False,
        )
    explicit = _first_classification(
        entry,
        *_molecule_gate_fields(entry),
        "gate_classification",
        "property_classification",
        "classification",
        "gate",
    )
    # ``verdict`` on a Step 4 row is the level-1 mark of the RDKit first-pass
    # screen, not a molecule decision, so a ``fail`` there is capped at
    # ``modify``.  Rejection has to arrive through the explicit gate fields
    # above -- which hold the per-category fold -- or through an explicit
    # ``rejected`` eligibility below.
    legacy = _marking_classification(entry.get("verdict")) or "keep"
    if legacy == "modify" and not _rdkit_verdict_is_molecular(entry):
        # The screen's ``flag`` rests only on a coverage gap (sascorer missing,
        # FilterCatalog unavailable, the pass raised). The executable policy forbids letting
        # a tool problem rewrite the chemistry, so it is demoted to
        # ``not_evaluated``: still not a clean keep, but it can no longer dominate
        # or invent a molecular ``modify``.  The gap blocks automatic nomination
        # through ``protocol_review_reasons`` instead.
        legacy = "not_evaluated"
    explicit = _cap_below_reject(
        max(
            (value for value in (explicit, legacy) if value is not None),
            key=_CLASS_SEVERITY.__getitem__,
        ),
        eliminates=eliminates,
        structural=structural,
    )

    # New property reports split molecular review reasons from protocol repair
    # reasons.  The latter block automatic nomination elsewhere but must not
    # rewrite a chemically passing molecule as MODIFY.  Legacy reports did not
    # carry that split, so their historical review_only => modify mapping stays.
    eligibility = str(entry.get("nomination_eligibility", "")).strip().casefold()
    if eligibility == "review_only":
        molecule_reasons = entry.get("molecule_review_reasons")
        protocol_reasons = entry.get("protocol_review_reasons")
        if isinstance(molecule_reasons, list) or isinstance(protocol_reasons, list):
            if not molecule_reasons:
                return explicit
        return max((explicit, "modify"), key=_CLASS_SEVERITY.__getitem__)
    return explicit


def property_data_review(
    entry: dict | None,
    *,
    required: bool = False,
) -> dict | None:
    """Return a protocol/data-review plan for property evidence that failed.

    A numeric property finding remains a molecule-level modification.  This
    helper is deliberately limited to missing rows and to backends that were
    planned and did not deliver, so an ordinary QED/Lipinski/etc. modify is not
    rerouted as a protocol problem.

    The executable policy draws the line this helper follows: a backend nobody selected,
    or a compute budget that only reached the head of the queue, is not a tool
    failure and needs no repair. Listing those here would mark every run
    configured without ADMET-AI or with a ``top_k`` smaller than the cohort as
    requiring evidence repair, confusing a configuration choice with a defect.
    """
    if entry is None:
        if not required:
            return None
        return {
            "reasons": ["missing_property_result"],
            "actions": [
                "Rerun property_filter and repair the identity join so the "
                "compound has a property result before nomination",
            ],
        }

    protocol_evidence = entry.get("property_protocol_evidence")
    reasons: list[str] = []
    actions: list[str] = []
    if isinstance(protocol_evidence, dict):
        if protocol_evidence.get("blocks_automatic_nomination") is True:
            reasons.extend(
                str(reason)
                for reason in (protocol_evidence.get("reasons") or [])
                if str(reason)
            )
            actions.extend(
                str(action)
                for action in (protocol_evidence.get("actions") or [])
                if str(action)
            )
            if not reasons:
                reasons.append("property_evidence_gap")
            if not actions:
                actions.append(
                    "Repair the planned property evidence before automatic nomination"
                )
        # When false, this explicitly accepts endpoint-local partial ADMET or
        # a ChemFH outage satisfied by RDKit.  Do not re-infer those statuses.
    else:
        eligibility = str(
            entry.get("nomination_eligibility") or ""
        ).strip().casefold()
        context = entry.get("metric_context")
        enabled_tools = {
            str(item).strip().casefold()
            for item in (
                context.get("enabled_tools")
                if isinstance(context, dict)
                and isinstance(context.get("enabled_tools"), (list, tuple, set))
                else []
            )
        }
        legacy_review_marker = eligibility == "review_only"

        # A backend status that means "ran for this molecule and did not deliver".
        # ``disabled`` is a configuration choice and ``not_selected_top_k`` /
        # ``not_eligible`` are decisions this run already recorded elsewhere.
        unusable = {"error", "partial", "missing", "incomplete"}

        admet = entry.get("admet")
        admet_status = str(
            (admet.get("status") if isinstance(admet, dict) else None)
            or "missing"
        ).strip().casefold()
        if admet_status in unusable and (
            admet_status != "missing"
            or "admet" in enabled_tools
            or legacy_review_marker
        ):
            reasons.append(f"admet_{admet_status}")
            actions.append(
                "Rerun or repair ADMET-AI and require status=complete before "
                "automatic nomination"
            )

        chemfh = entry.get("chemfh")
        chemfh_status = str(
            (chemfh.get("status") if isinstance(chemfh, dict) else None)
            or "missing"
        ).strip().casefold()
        if (
            chemfh_status in unusable
            and (
                chemfh_status != "missing"
                or "chemfh" in enabled_tools
                or legacy_review_marker
            )
            and not _legacy_rdkit_filtercatalog_fallback(entry)
        ):
            reasons.append(f"chemfh_{chemfh_status}")
            actions.append(
                "Rerun or repair ChemFH and require complete seven-endpoint "
                "evidence before automatic nomination"
            )

    synthesis = entry.get("synthesis")
    synthesis_status = str(
        (
            synthesis.get("status")
            if isinstance(synthesis, dict)
            else None
        )
        or "missing"
    ).strip().casefold()
    synthesis_review_statuses = (
        {"error"}
        if isinstance(protocol_evidence, dict)
        else {"error", "missing"}
    )
    if synthesis_status in synthesis_review_statuses:
        reasons.append(f"synthesis_{synthesis_status}")
        actions.append(
            "Rerun or repair AiZynthFinder and capture a complete route-review "
            "result"
        )

    if not reasons:
        return None
    return {
        "reasons": list(dict.fromkeys(reasons)),
        "actions": list(dict.fromkeys(actions)),
    }


def _legacy_rdkit_filtercatalog_fallback(entry: dict) -> bool:
    """Recognize the same explicit fallback in pre-evidence-object reports."""
    context = entry.get("metric_context")
    enabled = context.get("enabled_tools") if isinstance(context, dict) else None
    if not isinstance(enabled, (list, tuple, set)) or "pains" not in {
        str(item).strip().casefold() for item in enabled
    }:
        return False
    required_metrics = {
        "chemfh.other_assay_interference",
        "chemfh.blue_fluorescence",
        "chemfh.fluc_inhibitors",
        "chemfh.promiscuous_compounds",
        "chemfh.colloidal_aggregators",
        "chemfh.reactive_compounds",
        "chemfh.green_fluorescence",
        "chemfh.substructure_alert_count",
    }
    fallback_declared = False
    events = context.get("equivalent_fallbacks")
    if isinstance(events, (list, tuple)):
        for event in events:
            if not isinstance(event, dict):
                continue
            replaces = str(
                event.get("replaces") or event.get("for_tool") or ""
            ).strip().casefold()
            tool = "".join(
                character
                for character in str(event.get("tool") or "").casefold()
                if character.isalnum()
            )
            status = str(event.get("status") or "").strip().casefold()
            covered = {
                str(item).strip()
                for item in (event.get("covered_metric_ids") or [])
                if str(item).strip()
            }
            try:
                coverage = float(event.get("coverage"))
                confidence = float(event.get("confidence"))
            except (TypeError, ValueError):
                coverage = confidence = -1.0
            digest = str(
                event.get("artifact_sha256")
                or event.get("catalog_sha256")
                or ""
            ).strip().casefold()
            if (
                replaces == "chemfh"
                and tool == "rdkitfiltercatalog"
                and status in {"complete", "success"}
                and event.get("equivalent") is True
                and event.get("approved") is True
                and required_metrics.issubset(covered)
                and coverage == 1.0
                and 0 < confidence <= 1
                and bool(
                    event.get("tool_version")
                    or event.get("catalog_version")
                )
                and bool(event.get("mapping_version"))
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
            ):
                fallback_declared = True
                break
    if not fallback_declared:
        return False
    alerts = entry.get("structural_alerts")
    if not isinstance(alerts, dict):
        return False
    required_catalogs = {"PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"}
    catalogs = {
        str(item).strip().upper()
        for item in (alerts.get("catalogs") or [])
        if str(item).strip()
    }
    return bool(
        alerts.get("available") is True
        and isinstance(alerts.get("has_alert"), bool)
        and isinstance(alerts.get("matches"), list)
        and not (alerts.get("catalog_errors") or [])
        and required_catalogs.issubset(catalogs)
    )


def _molecule_gate_fields(entry: dict) -> tuple[str, ...]:
    """Name the gate field that holds a molecular conclusion and nothing else.

    ``metric_gate`` is the current two-level pass/flag/fail fold.  Coverage is
    published separately as ``coverage_gate``/``rank_status`` and never enters
    this molecular conclusion.  Prefer the explicit
    ``metric_gate_precoverage`` alias when reading records from a transition
    build; otherwise use the current ``metric_gate``.
    """
    if "metric_gate_precoverage" in entry:
        return ("metric_gate_precoverage",)
    return ("metric_gate",)


def _first_classification(entry: dict, *fields: str) -> str | None:
    for field in fields:
        value = _classification_value(entry.get(field))
        if value is not None:
            return value
    return None


def _classification_value(value: object) -> str | None:
    if value is None:
        return None
    return _CLASS_ALIASES.get(str(value).strip().casefold())


def _marking_classification(value: object) -> str | None:
    """Read a per-metric mark, capping it below rejection.

    Level-1 ``pass``/``flag``/``fail`` marks describe one metric on one molecule.
    Only the per-category fold may reject, so a ``fail`` read from
    a marking field becomes ``modify`` here.
    """
    classification = _classification_value(value)
    if classification == "reject":
        return "modify"
    return classification


def identify_bottleneck(eval_results: list[dict], prop_results: list[dict]) -> str:
    """Find the most common failure reason across all molecules."""
    reasons = Counter()

    for entry in eval_results:
        v = entry.get("verdict", "pass")
        if v in ("fail", "warn"):
            bp = entry.get("blueprint_match", {})
            pb = entry.get("posebusters", {})
            diagnosis = entry.get("diagnosis", [])

            if any(
                str(item).startswith("CLASH_FAIL") or str(item).startswith("CLASH:")
                for item in diagnosis
            ):
                reasons["steric_clashes"] += 1
            elif any(str(item).startswith("CLASH_WARN") for item in diagnosis):
                reasons["minor_steric_clashes"] += 1
            if any(str(item).startswith("HIGH_STRAIN") for item in diagnosis):
                reasons["high_strain"] += 1
            elif any(str(item).startswith("MODERATE_STRAIN") for item in diagnosis):
                reasons["moderate_strain"] += 1
            violation_count = pb.get("num_violations")
            if (
                violation_count is not None
                and violation_count > POSEBUSTERS_VIOLATION_THRESHOLD
            ):
                reasons["posebusters_violations"] += 1
            if bp.get("verdict") == "unresolved":
                reasons["blueprint_unresolved"] += 1
            elif bp.get("verdict") in ("fail", "warn"):
                reasons["blueprint_mismatch"] += 1

    for entry in prop_results:
        v = entry.get("verdict", "pass")
        if v in ("fail", "flag"):
            flags = entry.get("flags", [])
            for f in flags:
                if "PAINS" in f:
                    reasons["pains"] += 1
                elif "Lipinski" in f:
                    reasons["lipinski"] += 1
                elif "QED" in f:
                    reasons["low_qed"] += 1
                elif "SA" in f:
                    reasons["high_sa"] += 1

    if not reasons:
        return "none"
    return reasons.most_common(1)[0][0]


def is_scaffold_level(bottleneck: str) -> bool:
    """Determine if the bottleneck requires scaffold-level changes."""
    return bottleneck in ("steric_clashes", "posebusters_violations")


def is_protocol_level(bottleneck: str) -> bool:
    """Return whether another evaluation protocol is needed, not a new molecule."""
    return bottleneck in ("blueprint_unresolved",)


def aggregate_modification_hints(
    eval_results: list[dict],
    classifications: list[str],
) -> dict:
    """Aggregate modification hints for molecules classified as 'modify'."""
    target_molecules = []
    all_hints = []
    issue_counts = Counter()
    per_molecule: dict[str, dict] = {}

    for entry, cls in zip(eval_results, classifications):
        if cls != "modify":
            continue
        mol_name = entry.get("mol_name", "")
        target_molecules.append(mol_name)
        hint = entry.get("modification_hint", "")
        if hint:
            all_hints.append(hint)

        blueprint_match = entry.get("blueprint_match") or {}
        diagnosis = entry.get("diagnosis", [])
        blueprint_gaps = [
            d[len("MISSING:"):] for d in diagnosis if d.startswith("MISSING:")
        ]
        required_met = blueprint_match.get(
            "required_hits", blueprint_match.get("required_met", 0)
        )
        required_total = blueprint_match.get("required_total", 0)

        for diag in diagnosis:
            if diag.startswith("MISSING:"):
                issue_counts[diag] += 1
            elif diag.startswith("LOW_OCCUPANCY:"):
                issue_counts[diag] += 1

        primary_issue = ""
        if blueprint_match.get("verdict") == "unresolved":
            primary_issue = "UNRESOLVED_BLUEPRINT: evaluator capability mismatch"
            if not hint:
                hint = "Re-evaluate with a profiler capable of the required atom/role/geometry fields"
                all_hints.append(hint)
        elif blueprint_gaps:
            primary_issue = f"MISSING:{blueprint_gaps[0]}"
        elif any(d.startswith("HIGH_STRAIN") for d in diagnosis):
            primary_issue = next(d for d in diagnosis if d.startswith("HIGH_STRAIN"))
        elif any(d.startswith("CLASH") for d in diagnosis):
            primary_issue = next(d for d in diagnosis if d.startswith("CLASH"))

        per_molecule[mol_name] = {
            "primary_issue": primary_issue,
            "fix_instruction": hint,
            "required_met": f"{required_met}/{required_total}" if required_total else "N/A",
            "desired_weighted_score": blueprint_match.get("desired_weighted_score"),
            "hard_gate_status": blueprint_match.get("hard_gate_status"),
            "blueprint_gaps": blueprint_gaps,
            "diagnosis": diagnosis,
        }

    common_issues = [issue for issue, _ in issue_counts.most_common(5)]

    target_molecules = list(dict.fromkeys(target_molecules))
    suggestions = sorted(set(all_hints))[:5]

    return {
        "target_molecules": target_molecules,
        "common_issues": common_issues,
        "suggestions": suggestions,
        "per_molecule": per_molecule,
    }

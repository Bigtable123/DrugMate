"""Blueprint generation from normalized observations and reviewed evidence.

PLIP/ProLIF agreement raises detector confidence, but does not by itself establish
biological necessity.  ``required`` constraints therefore remain provenance-gated.
"""
from __future__ import annotations

import sys
import hashlib
import json
from pathlib import Path
from typing import Any

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import (
    legacy_residue_label,
    make_residue_record,
    parse_residue_record,
    residue_identity,
)
from constants import normalize_interaction_type
from interaction_observations import (
    PLIP_TO_UNIFIED,
    PROLIF_TO_UNIFIED,
    normalize_interaction_observations,
)


# --- Interaction type mapping ---

DISTANCE_TOLERANCE = 0.5
ANGLE_TOLERANCE = 15.0
DISTANCE_METRIC_BY_TYPE = {
    "H-bond": "heavy_atom",
    "Hydrophobic": "atom_atom",
    "VdW": "atom_atom",
    "Pi-stacking": "ring_centroid",
    "Edge-to-face": "ring_centroid",
    "Pi-cation": "charge_ring_center",
    "Salt-bridge": "charge_center",
    "Halogen-bond": "atom_atom",
    "Water-bridge": "acceptor_water",
    "Metal-complex": "atom_atom",
}


# --- Chain/insertion-code-safe residue identities ---


def _residue_record(value: Any, default_chain: str = "") -> dict | None:
    """Parse detector, curated, or canonical residue values conservatively."""
    if isinstance(value, dict) and "residue" in value:
        record = parse_residue_record(
            value.get("residue"),
            default_chain=str(value.get("chain") or default_chain),
        )
        if record is not None and value.get("icode") not in (None, ""):
            record = make_residue_record(
                record.get("resname", ""), record.get("resnum", ""),
                record.get("chain", ""), value.get("icode", ""),
            )
        return record
    return parse_residue_record(value, default_chain=default_chain)


def _residue_key(value: Any, default_chain: str = "") -> tuple[Any, ...] | None:
    record = _residue_record(value, default_chain)
    return residue_identity(record) if record is not None else None


def _residue_key_set(values: list[Any] | None, default_chain: str) -> set[tuple[Any, ...]]:
    return {
        key for value in (values or [])
        if (key := _residue_key(value, default_chain)) is not None
    }


def _record_from_key(key: tuple[Any, ...]) -> dict:
    chain, resnum, icode, resname = key
    return make_residue_record(resname, resnum, chain, icode)


def _residue_sort_key(key: tuple[Any, ...]) -> tuple[str, int, str, str]:
    chain, resnum, icode, resname = key
    try:
        number = int(resnum)
    except (TypeError, ValueError):
        number = 10**9
    return str(chain), number, str(icode), str(resname)


# --- Main generation logic ---

def generate_blueprint(
    plip_result: dict | None,
    prolif_result: dict | None,
    consensus_residues: list[Any],
    known_site_residues: list[Any] | None,
    known_required_interactions: list[dict] | None = None,
    chain: str = "A",
    blueprint_applicable: bool = True,
    interaction_observations: dict | None = None,
    blueprint_policy: dict | None = None,
    known_site_provenance: dict | None = None,
    pocket_support_records: list[dict] | None = None,
) -> dict:
    """Turn normalized observations into explicitly authorized constraints."""
    policy = {
        "auto_required_from_tool_consensus": False,
        "require_curated_provenance": True,
        **dict(blueprint_policy or {}),
    }
    # Hard constraints are consumed as exact downstream gates.  Their evidence
    # contract is therefore invariant even for direct (non-pipeline) callers.
    policy["require_curated_provenance"] = True
    observation_contract = interaction_observations or normalize_interaction_observations(
        plip_result, prolif_result, default_chain=chain,
    )
    norm_consensus = _residue_key_set(consensus_residues, chain)
    norm_known = _residue_key_set(known_site_residues, chain)
    pocket_support_counts = {
        key: int(item.get("support_count") or 0)
        for item in (pocket_support_records or [])
        if (key := _residue_key(item.get("residue_record") or item.get("residue"), chain)) is not None
    }
    observations = list(observation_contract.get("observations") or [])
    groups = list(observation_contract.get("consensus_groups") or [])
    profiler_failures = [
        tool for tool, status in (observation_contract.get("tool_statuses") or {}).items()
        if status.get("status") == "tool_failed"
        and status.get("critical_for_hard_gate", tool == "plip")
    ]
    supplemental_profiler_failures = [
        tool for tool, status in (observation_contract.get("tool_statuses") or {}).items()
        if status.get("status") == "tool_failed"
        and not status.get("critical_for_hard_gate", tool == "plip")
    ]

    if not observations and not known_required_interactions:
        if not blueprint_applicable:
            validation = {
                "applicable": False,
                "not_applicable": True,
                "valid_for_hard_gating": None,
                "hard_gate_applicable": False,
                "required_hard_gate_applicable": False,
                "desired_metric_applicable": False,
                "molecule_hard_gate_applicable": False,
                "blueprint_valid": True,
                "reason": "no reference ligand and no curated constraints",
                "warnings": _profiler_failure_warnings(
                    profiler_failures, supplemental_profiler_failures,
                ),
                "curated_constraint_errors": [],
                "profiler_failures": profiler_failures,
                "supplemental_profiler_failures": supplemental_profiler_failures,
            }
            note = "No reference ligand or curated constraints; blueprint is not applicable"
        else:
            validation = {
                "applicable": True,
                "not_applicable": False,
                "valid_for_hard_gating": False,
                "hard_gate_applicable": False,
                "required_hard_gate_applicable": False,
                "desired_metric_applicable": False,
                "molecule_hard_gate_applicable": False,
                "blueprint_valid": True,
                "reason": "no classifiable interaction evidence",
                "warnings": [
                    "No approved curated constraints were generated; automated observations never authorize a hard gate by default",
                    *_profiler_failure_warnings(
                        profiler_failures, supplemental_profiler_failures,
                    ),
                ],
                "curated_constraint_errors": [],
                "profiler_failures": profiler_failures,
                "supplemental_profiler_failures": supplemental_profiler_failures,
            }
            note = "No interactions detected by either tool — blueprint is empty"
        return {
            "schema_version": 2,
            "required": [],
            "required_groups": [],
            "desired": [],
            "optional": [],
            "rejected_curated_constraints": [],
            "evidence_summary": _make_observation_summary(
                observation_contract, norm_consensus, norm_known, known_site_provenance,
            ),
            "generation_policy": {
                "rule": "evidence_overlap_classification",
                "note": note,
            },
            "validation": validation,
        }

    required, rejected, curated_errors = _curated_constraints(
        known_required_interactions or [], observations, chain, policy,
    )
    curated_errors.extend(_curated_group_validation_errors(required))
    required_groups = _required_groups(required)
    desired, optional, unsafe_required = _automatic_constraints(
        observations,
        groups,
        norm_consensus,
        norm_known,
        observation_contract.get("input_evidence") or {},
        pocket_support_counts,
        allow_auto_required=_unsafe_auto_required_allowed(policy, known_site_provenance),
    )
    # A curated constraint owns its residue/type/role contract. Keep the linked
    # raw observations in the evidence artifact, but do not emit a duplicate
    # automatic desired/optional constraint for the same chemical feature.
    desired = [item for item in desired if not _covered_by_curated(item, required, chain)]
    optional = [item for item in optional if not _covered_by_curated(item, required, chain)]
    unsafe_required = [
        item for item in unsafe_required if not _covered_by_curated(item, required, chain)
    ]
    if unsafe_required:
        required.extend(unsafe_required)
        required_groups.extend(_required_groups(unsafe_required))

    blueprint_valid = not curated_errors
    valid = bool(required) and blueprint_valid and not profiler_failures
    desired_metric_applicable = bool(desired) and blueprint_valid
    warnings = []
    if not required:
        warnings.append(
            "No approved curated hard constraints; eligible automatic observations remain Desired or Optional"
        )
    if curated_errors:
        warnings.append("One or more curated constraints failed provenance/scope/review validation")
    warnings.extend(_profiler_failure_warnings(
        profiler_failures, supplemental_profiler_failures,
    ))
    return {
        "schema_version": 2,
        "required": required,
        "required_groups": required_groups,
        "desired": desired,
        "optional": optional,
        "rejected_curated_constraints": rejected,
        "evidence_summary": _make_observation_summary(
            observation_contract, norm_consensus, norm_known, known_site_provenance,
        ),
        "generation_policy": {
            "rule": "observations_then_authorized_constraints",
            "required_criteria": "approved_curated_constraint_with_provenance_and_scope",
            "candidate_required_criteria": "none; ProLIF consensus is audit-only",
            "desired_criteria": "eligible_observation + multi_detector_pocket_support",
            "optional_criteria": "outside_multi_detector_pocket OR weak_VdW OR unsupported_water/metal",
            **policy,
        },
        "validation": {
            "applicable": True,
            "not_applicable": False,
            "valid_for_hard_gating": valid,
            "hard_gate_applicable": valid,
            # The legacy fields above authorize only the curated Required
            # verdict.  Desired scoring is a separate metric contract whose
            # numeric hard/soft cells remain active even when Required is empty.
            "required_hard_gate_applicable": valid,
            "desired_metric_applicable": desired_metric_applicable,
            "molecule_hard_gate_applicable": bool(
                valid or desired_metric_applicable
            ),
            "blueprint_valid": blueprint_valid,
            "reason": (
                "approved curated Required hard gate available" if valid
                else "Required exact-match hard gating requires approved curated constraints and gate-eligible profiler evidence"
            ),
            "warnings": warnings,
            "curated_constraint_errors": curated_errors,
            "profiler_failures": profiler_failures,
            "supplemental_profiler_failures": supplemental_profiler_failures,
        },
    }


def _profiler_failure_warnings(
    profiler_failures: list[str],
    supplemental_profiler_failures: list[str],
) -> list[str]:
    warnings = []
    if profiler_failures:
        warnings.append(
            f"Interaction profiler failures: {', '.join(sorted(profiler_failures))}"
        )
    if supplemental_profiler_failures:
        warnings.append(
            "Supplemental profiler failures (hard-gate validity unaffected): "
            f"{', '.join(sorted(supplemental_profiler_failures))}"
        )
    return warnings


def _stable_id(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, ensure_ascii=True, default=str)
    return f"{prefix}_{hashlib.sha1(encoded.encode()).hexdigest()[:12]}"


def _covered_by_curated(candidate: dict, required: list[dict], chain: str) -> bool:
    candidate_key = _residue_key(candidate, chain)
    candidate_type = candidate.get("interaction_type")
    for curated in required:
        if _residue_key(curated, chain) != candidate_key:
            continue
        if curated.get("interaction_type") != candidate_type:
            continue
        roles_match = True
        for field in ("ligand_role", "protein_role"):
            expected = curated.get(field)
            if expected not in (None, "") and candidate.get(field) != expected:
                roles_match = False
                break
        if roles_match:
            return True
    return False


def _curated_validation_errors(item: dict, policy: dict, index: int) -> list[str]:
    errors = []
    if _residue_record(item) is None:
        errors.append("invalid residue")
    if not _normalize_explicit_type(str(item.get("interaction_type") or "")):
        errors.append("missing interaction_type")
    if item.get("interaction_subtype") not in (None, ""):
        errors.append(
            "interaction_subtype is not a supported hard-gate field; encode a supported "
            "canonical interaction_type or leave the subtype as observation metadata"
        )
    provenance = item.get("provenance")
    scope = item.get("scope")
    if provenance is not None and not isinstance(provenance, dict):
        errors.append("provenance must be an object")
        provenance = {}
    if scope is not None and not isinstance(scope, dict):
        errors.append("scope must be an object")
        scope = {}
    provenance = provenance or {}
    scope = scope or {}

    nested_reference = (
        provenance.get("reference") or provenance.get("source")
        or provenance.get("citation") or provenance.get("record_id")
    )
    if not str(item.get("source") or nested_reference or "").strip():
        errors.append("requires source or provenance reference")
    if not str(item.get("rationale") or provenance.get("rationale") or "").strip():
        errors.append("requires rationale")
    scope_value = item.get("binding_mode_scope") or next((
        scope.get(field) for field in
        ("binding_mode", "binding_mode_scope", "mechanism", "description", "applies_to")
        if str(scope.get(field) or "").strip()
    ), None)
    if not str(scope_value or "").strip():
        errors.append("requires binding-mode scope")
    if provenance and not str(provenance.get("source_type") or "").strip():
        errors.append("provenance.source_type is required")

    review = str(
        item.get("review_status")
        or (provenance or {}).get("review_status")
        or ""
    ).lower()
    if review != "approved":
        errors.append("review_status must be approved")
    return [f"known_required_interactions[{index}]: {message}" for message in errors]


def _normalized_curated_metadata(item: dict) -> tuple[dict, dict]:
    provenance = item.get("provenance")
    if not isinstance(provenance, dict):
        provenance = {
            "source_type": "curated",
            "source": str(item.get("source") or ""),
            "rationale": str(item.get("rationale") or ""),
        }
        for field in ("reference", "citation", "record_id"):
            if item.get(field):
                provenance[field] = item[field]
    else:
        provenance = dict(provenance)
        if item.get("rationale") and not provenance.get("rationale"):
            provenance["rationale"] = item["rationale"]
    scope = item.get("scope")
    if not isinstance(scope, dict):
        scope = {"binding_mode": str(item.get("binding_mode_scope") or "")}
    else:
        scope = dict(scope)
    return provenance, scope


def _curated_logic(item: dict) -> tuple[str | None, str, int | None]:
    logic_group = item.get("logic_group")
    if isinstance(logic_group, dict):
        group_id = logic_group.get("group_id") or logic_group.get("id")
        operator = logic_group.get("logic") or logic_group.get("operator")
        at_least_k = logic_group.get("at_least_k")
    else:
        group_id = logic_group or item.get("group_id")
        operator = None
        at_least_k = None
    operator = (
        item.get("logic_operator") or item.get("group_logic") or item.get("logic")
        or operator or "all_of"
    )
    if at_least_k is None:
        at_least_k = item.get("at_least_k")
    return (
        str(group_id) if group_id not in (None, "") else None,
        str(operator),
        int(at_least_k) if at_least_k is not None else None,
    )


def _matching_observations(item: dict, observations: list[dict], chain: str) -> list[dict]:
    residue_key = _residue_key(item, chain)
    interaction_type = _normalize_explicit_type(str(item.get("interaction_type") or ""))
    matches = []
    for observation in observations:
        if residue_identity(observation.get("residue")) != residue_key:
            continue
        if observation.get("interaction_type") != interaction_type:
            continue
        if item.get("ligand_role") and item.get("ligand_role") != observation.get("ligand_role"):
            continue
        if item.get("protein_role") and item.get("protein_role") != observation.get("protein_role"):
            continue
        matches.append(observation)
    return matches


def _range_from_observations(observations: list[dict], field: str, tolerance: float) -> list[float] | None:
    values = []
    for observation in observations:
        geometry = observation.get("geometry") or {}
        aliases = (
            ("heavy_atom_distance", "distance") if field == "distance"
            else (field,)
        )
        value = next((geometry.get(key) for key in aliases if geometry.get(key) is not None), None)
        if value is not None:
            try:
                values.append(float(value))
            except (TypeError, ValueError):
                pass
    if not values:
        return None
    center = sum(values) / len(values)
    if field == "angle":
        return [round(max(0.0, center - tolerance), 1), round(min(180.0, center + tolerance), 1)]
    return [round(max(0.0, center - tolerance), 2), round(center + tolerance, 2)]


def _curated_constraints(
    items: list[dict], observations: list[dict], chain: str, policy: dict,
) -> tuple[list[dict], list[dict], list[str]]:
    required = []
    rejected = []
    errors = []
    for index, raw_item in enumerate(items):
        item = dict(raw_item or {})
        item_errors = _curated_validation_errors(item, policy, index)
        if item_errors:
            errors.extend(item_errors)
            rejected.append({"input": item, "errors": item_errors})
            continue
        record = _residue_record(item, chain)
        matched = _matching_observations(item, observations, chain)
        provenance, scope = _normalized_curated_metadata(item)
        group_id, group_logic, at_least_k = _curated_logic(item)
        constraint_id = str(item.get("constraint_id") or _stable_id("req", {
            "residue": record,
            "type": item.get("interaction_type"),
            "roles": [item.get("ligand_role"), item.get("protein_role")],
            "provenance": provenance,
        }))
        evidence_tools = sorted({obs.get("tool") for obs in matched if obs.get("tool")})
        protein_atoms = item.get("protein_atoms") or [
            {"name": str(name)} for name in (item.get("protein_atom_names") or [])
        ]
        ligand_atoms = item.get("ligand_atoms") or [
            {"name": str(name)} for name in (item.get("ligand_atom_names") or [])
        ]
        entry = {
            "constraint_id": constraint_id,
            "residue": legacy_residue_label(record),
            "chain": record.get("chain", ""),
            "icode": record.get("icode", ""),
            "interaction_type": _normalize_explicit_type(str(item.get("interaction_type"))),
            "ligand_role": item.get("ligand_role"),
            "protein_role": item.get("protein_role"),
            "protein_atoms": protein_atoms,
            "ligand_atoms": ligand_atoms,
            "evidence_source": ["curated", *evidence_tools],
            "observation_ids": [obs["observation_id"] for obs in matched],
            "confidence": "high",
            "provenance": provenance,
            "scope": scope,
            "review_status": "approved",
            "hard_gate_authorized": True,
            "group_id": group_id,
            "group_logic": group_logic,
            "at_least_k": at_least_k,
        }
        if item.get("mediator") not in (None, ""):
            entry["mediator"] = item["mediator"]
        # Curated biological authority does not automatically authorize the
        # exact geometry of one static input structure as a hard constraint.
        # Compile geometry only when the curator supplied it explicitly.
        distance_range = item.get("distance_range")
        angle_range = item.get("angle_range")
        if distance_range:
            entry["distance_range"] = distance_range
            entry.setdefault("geometry", {}).update({
                "distance_metric": DISTANCE_METRIC_BY_TYPE.get(
                    entry["interaction_type"], "atom_atom",
                ),
                "distance_unit": "angstrom",
            })
        if angle_range:
            entry["angle_range"] = angle_range
            entry.setdefault("geometry", {})["angle_unit"] = "degree"
        required.append(entry)
    return required, rejected, errors


def _required_groups(required: list[dict]) -> list[dict]:
    grouped: dict[str, dict] = {}
    for entry in required:
        group_id = str(entry.get("group_id") or f"group_{entry['constraint_id']}")
        group = grouped.setdefault(group_id, {
            "group_id": group_id,
            "logic": str(entry.get("group_logic") or "all_of"),
            "members": [],
        })
        group["members"].append(entry["constraint_id"])
        if entry.get("at_least_k") is not None:
            group["at_least_k"] = int(entry["at_least_k"])
    for group in grouped.values():
        if group["logic"] == "any_of":
            group.setdefault("at_least_k", 1)
        elif group["logic"] == "all_of":
            group.setdefault("at_least_k", len(group["members"]))
        elif group["logic"] == "at_least_k":
            group.setdefault("at_least_k", 1)
    return list(grouped.values())


def _curated_group_validation_errors(required: list[dict]) -> list[str]:
    """Reject ambiguous group contracts instead of silently taking first/last."""
    errors = []
    seen_ids = set()
    grouped: dict[str, list[dict]] = {}
    for entry in required:
        constraint_id = str(entry.get("constraint_id") or "")
        if constraint_id in seen_ids:
            errors.append(f"duplicate required constraint_id: {constraint_id}")
        seen_ids.add(constraint_id)
        if entry.get("group_id"):
            grouped.setdefault(str(entry["group_id"]), []).append(entry)
    for group_id, members in grouped.items():
        operators = {str(member.get("group_logic") or "all_of") for member in members}
        if len(operators) != 1:
            errors.append(
                f"required group {group_id!r} has conflicting operators: "
                f"{', '.join(sorted(operators))}"
            )
        k_values = {
            int(member["at_least_k"])
            for member in members if member.get("at_least_k") is not None
        }
        if len(k_values) > 1:
            errors.append(
                f"required group {group_id!r} has conflicting at_least_k values: "
                f"{', '.join(map(str, sorted(k_values)))}"
            )
        if k_values and next(iter(k_values)) > len(members):
            errors.append(
                f"required group {group_id!r} at_least_k exceeds member count"
            )
        operator = next(iter(operators)) if len(operators) == 1 else None
        if operator in {"all_of", "any_of"} and k_values:
            errors.append(
                f"required group {group_id!r} may set at_least_k only with at_least_k logic"
            )
    return [f"known_required_interactions: {message}" for message in errors]


def _unsafe_auto_required_allowed(policy: dict, known_site_provenance: dict | None) -> bool:
    provenance = known_site_provenance or {}
    return bool(
        policy.get("auto_required_from_tool_consensus") is True
        and provenance.get("independent_of_reference_structure") is True
        and str(provenance.get("review_status") or "").lower() == "approved"
    )


def _context_supported(input_evidence: dict, name: str) -> bool:
    evidence = input_evidence.get(name) or input_evidence.get(f"{name}_report") or {}
    if not (
        isinstance(evidence, dict)
        and evidence.get("status") in {"loaded", "loaded_json", "available"}
        and isinstance(evidence.get("content"), dict)
    ):
        return False
    content = evidence["content"]
    if name == "water":
        try:
            if int(content.get("selected_count") or 0) > 0:
                return True
        except (TypeError, ValueError):
            pass
        return any(
            isinstance(item, dict) and item.get("kept") is True
            for item in (content.get("waters") or [])
        )
    if name == "metal":
        try:
            count = int(content.get("metal_count") or 0)
        except (TypeError, ValueError):
            count = 0
        return count > 0 and bool(content.get("sites"))
    return bool(content)


def _constraint_from_observations(
    member_observations: list[dict], *, confidence: str, pocket_support_count: int,
) -> dict:
    representative = member_observations[0]
    record = representative.get("residue") or {}
    tools = sorted({obs.get("tool") for obs in member_observations if obs.get("tool")})
    payload = {
        "residue": record,
        "interaction_type": representative.get("interaction_type"),
        "roles": [representative.get("ligand_role"), representative.get("protein_role")],
        "atoms": representative.get("protein_atoms") or [],
        "observations": sorted(obs["observation_id"] for obs in member_observations),
    }
    entry = {
        "constraint_id": _stable_id("auto", payload),
        "residue": legacy_residue_label(record),
        "chain": record.get("chain", ""),
        "icode": record.get("icode", ""),
        "interaction_type": representative.get("interaction_type"),
        "ligand_role": representative.get("ligand_role"),
        "protein_role": representative.get("protein_role"),
        "protein_atoms": _transferable_protein_atoms(
            representative.get("protein_atoms") or []
        ),
        # Atom identity on the protein is transferable across candidate
        # ligands.  Reference-ligand atom indices/names are not, so retain
        # them only as provenance rather than as a pose-matching constraint.
        "reference_ligand_atoms": representative.get("ligand_atoms") or [],
        "mediator": _mediator_type(representative.get("mediator")),
        "evidence_source": tools,
        "observation_ids": payload["observations"],
        "confidence": confidence,
        "pocket_support_count": pocket_support_count,
        "hard_gate_authorized": False,
    }
    distance_range = _range_from_observations(member_observations, "distance", DISTANCE_TOLERANCE)
    angle_range = _range_from_observations(member_observations, "angle", ANGLE_TOLERANCE)
    if distance_range:
        entry["distance_range"] = distance_range
        entry.setdefault("geometry", {}).update({
            "distance_metric": DISTANCE_METRIC_BY_TYPE.get(entry["interaction_type"], "atom_atom"),
            "distance_unit": "angstrom",
        })
    if angle_range:
        entry["angle_range"] = angle_range
        entry.setdefault("geometry", {})["angle_unit"] = "degree"
    return entry


def _transferable_protein_atoms(atoms: list[dict]) -> list[dict]:
    """Keep stable receptor atom names, never tool-local atom indices."""
    names = []
    seen = set()
    for atom in atoms or []:
        if not isinstance(atom, dict):
            continue
        name = str(atom.get("name") or atom.get("atom_name") or "").strip()
        if name and name.upper() not in seen:
            seen.add(name.upper())
            names.append({"name": name})
    return names


def _mediator_type(value: Any) -> str:
    if isinstance(value, dict):
        value = value.get("type") or value.get("name")
    return str(value or "direct")


def _automatic_constraints(
    observations: list[dict],
    groups: list[dict],
    norm_consensus: set,
    norm_known: set,
    input_evidence: dict,
    pocket_support_counts: dict,
    *,
    allow_auto_required: bool,
) -> tuple[list[dict], list[dict], list[dict]]:
    eligible_observations = [
        obs for obs in observations
        if str(obs.get("tool") or "").lower() != "prolif"
        and obs.get("blueprint_eligible", True) is True
    ]
    by_id = {obs["observation_id"]: obs for obs in eligible_observations}
    consumed = set()
    desired, optional, unsafe_required = [], [], []

    for group in groups:
        member_ids = list(group.get("member_observation_ids", []) or [])
        members = [by_id[item] for item in member_ids if item in by_id]
        if len(members) < 2 or len(members) != len(member_ids):
            continue
        consumed.update(obs["observation_id"] for obs in members)
        residue_key = residue_identity(members[0].get("residue"))
        in_known = residue_key in norm_known
        in_pocket = residue_key in norm_consensus
        entry = _constraint_from_observations(
            members,
            confidence="high",
            pocket_support_count=pocket_support_counts.get(residue_key, 2 if in_pocket else 0),
        )
        entry["consensus_id"] = group.get("consensus_id")
        if in_known:
            entry["candidate_required"] = True
            entry["candidate_required_reason"] = (
                "strict atom/role/geometry tool concordance at a known-site residue; "
                "curated biological authorization is still required"
            )
        if allow_auto_required and in_known:
            entry.update({
                "hard_gate_authorized": True,
                "authorization": "unsafe_explicit_auto_required_policy",
                "provenance": {"source": "known_site_provenance"},
                "scope": {"description": "explicit unsafe auto-required override"},
                "review_status": "approved",
                "group_logic": "all_of",
            })
            unsafe_required.append(entry)
        elif entry["interaction_type"] == "VdW":
            optional.append(entry)
        elif entry["interaction_type"] == "Water-bridge" and not _context_supported(input_evidence, "water"):
            entry["confidence"] = "medium"
            entry["conditional"] = True
            optional.append(entry)
        elif entry["interaction_type"] == "Metal-complex" and not _context_supported(input_evidence, "metal"):
            entry["confidence"] = "medium"
            entry["conditional"] = True
            optional.append(entry)
        elif in_known or in_pocket:
            desired.append(entry)
        else:
            entry["confidence"] = "medium"
            optional.append(entry)

    for observation in eligible_observations:
        if observation["observation_id"] in consumed:
            continue
        residue_key = residue_identity(observation.get("residue"))
        in_pocket = residue_key in norm_consensus
        entry = _constraint_from_observations(
            [observation], confidence="medium" if in_pocket else "low",
            pocket_support_count=pocket_support_counts.get(residue_key, 2 if in_pocket else 0),
        )
        interaction_type = entry["interaction_type"]
        if interaction_type == "VdW":
            optional.append(entry)
        elif interaction_type == "Water-bridge":
            entry["conditional"] = True
            (desired if in_pocket and _context_supported(input_evidence, "water") else optional).append(entry)
        elif interaction_type == "Metal-complex":
            entry["conditional"] = True
            (desired if in_pocket and _context_supported(input_evidence, "metal") else optional).append(entry)
        elif in_pocket:
            desired.append(entry)
        else:
            optional.append(entry)
    return desired, optional, unsafe_required


def _make_observation_summary(
    contract: dict, norm_consensus: set, norm_known: set,
    known_site_provenance: dict | None,
) -> dict:
    observations = contract.get("observations") or []
    tools = {str(obs.get("tool")) for obs in observations if obs.get("tool")}
    return {
        "detected_observations_count": len(observations),
        "blueprint_eligible_observations_count": sum(
            obs.get("blueprint_eligible", True) is True for obs in observations
        ),
        "evidence_only_observations_count": sum(
            obs.get("blueprint_eligible", True) is not True for obs in observations
        ),
        "strict_cross_tool_consensus_count": len(contract.get("consensus_groups") or []),
        "tools_with_detected_observations": sorted(tools),
        "tool_statuses": contract.get("tool_statuses") or {},
        "interaction_consensus_pocket_residues_count": len(norm_consensus),
        "known_site_residues_count": len(norm_known),
        "known_site_provenance": known_site_provenance,
    }


def _normalize_explicit_type(value: str) -> str:
    return normalize_interaction_type(value)


# --- Internal parsers ---

def _parse_plip(
    plip_result: dict | None, *, default_chain: str = "",
) -> dict[tuple[Any, ...], dict[str, dict]]:
    """Parse PLIP evidence without collapsing chains or insertion codes."""
    if not plip_result or not plip_result.get("success"):
        return {}

    result: dict[tuple[Any, ...], dict[str, dict]] = {}
    for interaction in plip_result.get("interactions_detail", []):
        raw_type = interaction.get("type", "")
        unified = PLIP_TO_UNIFIED.get(raw_type)
        if not unified:
            continue

        norm = _residue_key(interaction, default_chain)
        if norm is None:
            continue

        detail = {}
        if interaction.get("distance") is not None:
            detail["distance"] = interaction["distance"]
        if interaction.get("angle") is not None:
            detail["angle"] = interaction["angle"]
        if interaction.get("distance_aw") is not None:
            detail["distance"] = interaction["distance_aw"]

        if norm not in result:
            result[norm] = {}
        if unified not in result[norm]:
            result[norm][unified] = detail
        else:
            existing = result[norm][unified]
            if detail.get("distance") and existing.get("distance"):
                if detail["distance"] < existing["distance"]:
                    result[norm][unified] = detail

    return result


def _parse_prolif(
    prolif_result: dict | None, *, default_chain: str = "",
) -> dict[tuple[Any, ...], dict[str, dict]]:
    """Parse ProLIF evidence without collapsing ``LEU83.A`` and ``LEU83.B``."""
    if not prolif_result or prolif_result.get("validation_failed"):
        return {}

    result: dict[tuple[Any, ...], dict[str, dict]] = {}
    summary = prolif_result.get("summary", {})
    for res_id, interaction_types in summary.items():
        norm = _residue_key(res_id, default_chain)
        if norm is None:
            continue
        if norm not in result:
            result[norm] = {}
        for itype in interaction_types:
            unified = PROLIF_TO_UNIFIED.get(itype)
            if unified and unified not in result[norm]:
                result[norm][unified] = {}

    return result


def _make_evidence_summary(
    plip_map: dict, prolif_map: dict, norm_consensus: set, norm_known: set
) -> dict:
    plip_residues = set(plip_map.keys())
    prolif_residues = set(prolif_map.keys())
    both = plip_residues & prolif_residues

    return {
        "plip_residues_count": len(plip_residues),
        "prolif_residues_count": len(prolif_residues),
        "consensus_pocket_residues_count": len(norm_consensus),
        "known_site_residues_count": len(norm_known),
        "both_tools_detected": len(both),
        "plip_only": len(plip_residues - prolif_residues),
        "prolif_only": len(prolif_residues - plip_residues),
    }

"""Canonical, atom-aware interaction observations for Step 1.

Interaction profilers report *observations*, not design requirements.  This
module keeps their raw chemical roles and geometry, records whether each tool
was evaluable, and only calls two observations concordant when their features
are chemically compatible.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any


_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import (  # noqa: E402
    legacy_residue_label,
    make_residue_record,
    parse_residue_record,
    residue_identity,
)


PLIP_TO_UNIFIED = {
    "HBond": "H-bond",
    "Hydrophobic": "Hydrophobic",
    "PiStacking": "Pi-stacking",
    "PiCation": "Pi-cation",
    "SaltBridge": "Salt-bridge",
    "HalogenBond": "Halogen-bond",
    "XBDonor": "Halogen-bond",
    "XBAcceptor": "Halogen-bond",
    "MetalDonor": "Metal-complex",
    "MetalAcceptor": "Metal-complex",
    "WaterBridge": "Water-bridge",
    "MetalComplex": "Metal-complex",
}

PROLIF_TO_UNIFIED = {
    "HBDonor": "H-bond",
    "HBAcceptor": "H-bond",
    "Hydrophobic": "Hydrophobic",
    "PiStacking": "Pi-stacking",
    "PiCation": "Pi-cation",
    "CationPi": "Pi-cation",
    "Cationic": "Salt-bridge",
    "Anionic": "Salt-bridge",
    "VdWContact": "VdW",
    "FaceToFace": "Pi-stacking",
    "EdgeToFace": "Edge-to-face",
    "FaceToEdge": "Edge-to-face",
    "HalogenBond": "Halogen-bond",
    "XBDonor": "Halogen-bond",
    "XBAcceptor": "Halogen-bond",
    "MetalDonor": "Metal-complex",
    "MetalAcceptor": "Metal-complex",
    "WaterBridge": "Water-bridge",
}

SUPPLEMENTAL_TO_UNIFIED = {
    "HBond": "H-bond",
    "WeakHBond": "H-bond",
    "Hydrophobic": "Hydrophobic",
    "PiStacking": "Pi-stacking",
    "PiCation": "Pi-cation",
    "SaltBridge": "Salt-bridge",
    "VdW": "VdW",
    "EdgeToFace": "Edge-to-face",
    "HalogenBond": "Halogen-bond",
    "MetalComplex": "Metal-complex",
    "WaterBridge": "Water-bridge",
}

TOOL_TYPE_MAPS = {
    "PLIP": PLIP_TO_UNIFIED,
    "ProLIF": PROLIF_TO_UNIFIED,
    "Arpeggio": SUPPLEMENTAL_TO_UNIFIED,
    "LUNA": SUPPLEMENTAL_TO_UNIFIED,
    "ODDT": SUPPLEMENTAL_TO_UNIFIED,
}

ROLE_SENSITIVE_TYPES = {
    "H-bond", "Salt-bridge", "Pi-cation", "Halogen-bond", "Water-bridge",
    "Metal-complex",
}
GEOMETRY_DISTANCE_TOLERANCE = {
    "H-bond": 0.8,
    "Hydrophobic": 1.0,
    "Pi-stacking": 1.2,
    "Pi-cation": 1.2,
    "Salt-bridge": 1.2,
    "Halogen-bond": 0.8,
    "Water-bridge": 1.0,
    "Metal-complex": 0.5,
}
ANGLE_TOLERANCE = 25.0


def _record(value: Any, default_chain: str) -> dict | None:
    if isinstance(value, dict):
        residue_value = value.get("residue_record") or value.get("residue")
        chain = str(value.get("chain") or default_chain)
        if residue_value is None and value.get("resnum") is not None:
            residue_value = {
                "resname": value.get("resname", ""),
                "resnum": value.get("resnum"),
                "chain": chain,
                "icode": value.get("icode", ""),
            }
        parsed = parse_residue_record(residue_value, default_chain=chain)
        if parsed is not None and value.get("icode") not in (None, ""):
            parsed = make_residue_record(
                parsed.get("resname", ""), parsed.get("resnum", ""),
                parsed.get("chain", ""), value.get("icode", ""),
            )
        return parsed
    return parse_residue_record(value, default_chain=default_chain)


def _json_value(value: Any) -> Any:
    """Convert numpy/RDKit-ish scalar containers without importing either."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(v) for v in value]
    if hasattr(value, "item"):
        try:
            return _json_value(value.item())
        except (TypeError, ValueError):
            pass
    return str(value)


def _atoms(value: Any, namespace: str) -> list[dict]:
    atoms = []
    for atom in value or []:
        if isinstance(atom, dict):
            entry = {str(key): _json_value(val) for key, val in atom.items() if val is not None}
        else:
            entry = {"index": _json_value(atom)}
        entry.setdefault("namespace", namespace)
        atoms.append(entry)
    return atoms


def _geometry(detail: dict) -> dict:
    raw = dict(detail.get("geometry") or {})
    aliases = {
        "distance": "distance",
        "distance_ad": "heavy_atom_distance",
        "distance_ah": "acceptor_hydrogen_distance",
        "angle": "angle",
        "distance_aw": "acceptor_water_distance",
        "distance_dw": "donor_water_distance",
        "d_angle": "donor_angle",
        "w_angle": "water_angle",
        "offset": "offset",
        "acc_angle": "acceptor_angle",
        "don_angle": "donor_angle",
        "rms": "coordination_rms",
    }
    for source, target in aliases.items():
        if detail.get(source) is not None and target not in raw:
            raw[target] = detail[source]
    clean = {}
    for key, value in raw.items():
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            clean[str(key)] = round(number, 4)
    if clean:
        clean["distance_unit"] = "angstrom"
        clean["angle_unit"] = "degree"
    return clean


def _roles(tool: str, raw_type: str, detail: dict) -> tuple[str | None, str | None]:
    ligand_role = detail.get("ligand_role")
    protein_role = detail.get("protein_role")
    if ligand_role or protein_role:
        return ligand_role, protein_role

    if tool == "PLIP":
        donor = str(detail.get("donor_type") or "").lower()
        if raw_type == "HBond" and donor in {"protein", "ligand"}:
            return (("acceptor", "donor") if donor == "protein" else ("donor", "acceptor"))
        if raw_type == "SaltBridge":
            protein_charge = str(detail.get("protein_charge") or "").lower()
            if protein_charge in {"positive", "cation"}:
                return "anion", "cation"
            if protein_charge in {"negative", "anion"}:
                return "cation", "anion"
        if raw_type == "HalogenBond":
            return "halogen_donor", "halogen_acceptor"
        if raw_type == "WaterBridge" and donor in {"protein", "ligand"}:
            return (("acceptor", "donor") if donor == "protein" else ("donor", "acceptor"))

    if tool == "ProLIF":
        role_map = {
            "HBDonor": ("donor", "acceptor"),
            "HBAcceptor": ("acceptor", "donor"),
            "Cationic": ("cation", "anion"),
            "Anionic": ("anion", "cation"),
            "PiCation": ("aromatic", "cation"),
            "CationPi": ("cation", "aromatic"),
            "HalogenBond": ("halogen_donor", "halogen_acceptor"),
            "XBDonor": ("halogen_donor", "halogen_acceptor"),
            "XBAcceptor": ("halogen_acceptor", "halogen_donor"),
            "MetalDonor": ("metal_center", "metal_coordinator"),
            "MetalAcceptor": ("metal_coordinator", "metal_center"),
        }
        if raw_type in role_map:
            return role_map[raw_type]
    return None, None


def _observation(tool: str, detail: dict, default_chain: str) -> dict | None:
    raw_type = str(detail.get("type") or detail.get("interaction_type") or "")
    interaction_type = (TOOL_TYPE_MAPS.get(tool) or {}).get(raw_type)
    if not interaction_type:
        return None
    record = _record(detail, default_chain)
    if record is None:
        return None
    ligand_role, protein_role = _roles(tool, raw_type, detail)
    namespace = str(detail.get("atom_index_namespace") or f"{tool.lower()}_native")
    payload = {
        "status": "detected",
        "tool": tool,
        "tool_version": detail.get("tool_version"),
        "interaction_type": interaction_type,
        "interaction_subtype": raw_type,
        "residue": record,
        "residue_label": legacy_residue_label(record),
        "ligand_role": ligand_role,
        "protein_role": protein_role,
        "ligand_atoms": _atoms(detail.get("ligand_atoms"), namespace),
        "protein_atoms": _atoms(detail.get("protein_atoms"), namespace),
        "mediator": _json_value(detail.get("mediator")) or (
            "water" if interaction_type == "Water-bridge" else
            "metal" if interaction_type == "Metal-complex" else "direct"
        ),
        "geometry": _geometry(detail),
        "blueprint_eligible": tool != "ProLIF" and bool(detail.get("blueprint_eligible", True)),
        "quality_flags": sorted(set(detail.get("quality_flags") or [])),
        "raw_metadata": _json_value(detail.get("raw_metadata") or {}),
    }
    stable = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    payload["observation_id"] = (
        f"obs_{tool.lower()}_{hashlib.sha1(stable.encode('utf-8')).hexdigest()[:12]}"
    )
    return payload


def _tool_observations(tool: str, result: dict | None, default_chain: str) -> list[dict]:
    if not result:
        return []
    details = result.get("interactions_detail") or []
    observations = []
    for detail in details:
        if not isinstance(detail, dict):
            continue
        enriched = dict(detail)
        enriched.setdefault("tool_version", result.get("tool_version"))
        enriched.setdefault("blueprint_eligible", result.get("blueprint_eligible", True))
        if observation := _observation(tool, enriched, default_chain):
            observation["quality_flags"] = sorted(set(
                list(observation.get("quality_flags") or [])
                + list(result.get("quality_flags") or [])
            ))
            if result.get("chemical_representation"):
                observation["chemical_representation"] = _json_value(
                    result.get("chemical_representation")
                )
            observations.append(observation)

    # Backward compatibility for older ProLIF files that only contain summary.
    if tool == "ProLIF" and not details:
        for residue, types in (result.get("summary") or {}).items():
            for raw_type in types or []:
                if observation := _observation(
                    tool,
                    {"residue": residue, "type": raw_type,
                     "tool_version": result.get("tool_version")},
                    default_chain,
                ):
                    observation["quality_flags"] = ["legacy_summary_without_atom_geometry"]
                    observations.append(observation)
    return observations


def _atom_tokens(atoms: list[dict], *, names_only: bool = False) -> set[str]:
    tokens = set()
    for atom in atoms or []:
        name = str(atom.get("name") or "").strip().upper()
        if name:
            tokens.add(f"name:{name}")
        if not names_only and atom.get("canonical_id") not in (None, ""):
            tokens.add(f"canonical:{atom['canonical_id']}")
    return tokens


def _ligand_atom_tokens(atoms: list[dict]) -> set[str]:
    """Map cross-tool ligand atom records to input-order/element features."""
    import re

    tokens = set()
    for atom in atoms or []:
        if not isinstance(atom, dict):
            continue
        name = str(atom.get("name") or "").strip().upper()
        compact_name = re.sub(r"\s+", "", name)
        element = str(atom.get("element") or "").strip().upper()
        ordinal = None
        if atom.get("residue_index") is not None:
            try:
                ordinal = int(atom["residue_index"]) + 1
            except (TypeError, ValueError):
                pass
        if ordinal is None and compact_name:
            match = re.match(r"^[A-Z]+(\d+)$", compact_name)
            if match:
                ordinal = int(match.group(1))
                if not element:
                    element = re.match(r"^[A-Z]+", compact_name).group(0)
        if not element and atom.get("atomic_number") is not None:
            element = f"Z{atom['atomic_number']}"
        if ordinal is not None:
            tokens.add(f"ordinal:{ordinal}:{element}")
        elif compact_name:
            tokens.add(f"name:{compact_name}")
    return tokens


def _mediator_feature(value: Any) -> tuple[str, str]:
    if isinstance(value, dict):
        return (
            str(value.get("type") or value.get("name") or "").lower(),
            str(value.get("element") or "").upper(),
        )
    return str(value or "").lower(), ""


def _primary_distance(geometry: dict) -> float | None:
    for key in ("heavy_atom_distance", "distance", "acceptor_water_distance"):
        if geometry.get(key) is not None:
            try:
                return float(geometry[key])
            except (TypeError, ValueError):
                pass
    return None


def observations_match(first: dict, second: dict) -> tuple[bool, dict]:
    """Return strict cross-tool feature concordance plus match diagnostics."""
    checks = {
        "same_residue": residue_identity(first.get("residue")) == residue_identity(second.get("residue")),
        "same_interaction_type": first.get("interaction_type") == second.get("interaction_type"),
        "role_compatible": True,
        "atom_compatible": None,
        "ligand_atom_compatible": None,
        "mediator_compatible": None,
        "geometry_compatible": None,
    }
    if not checks["same_residue"] or not checks["same_interaction_type"]:
        return False, checks

    interaction_type = str(first.get("interaction_type"))
    if interaction_type in ROLE_SENSITIVE_TYPES:
        roles1 = (first.get("ligand_role"), first.get("protein_role"))
        roles2 = (second.get("ligand_role"), second.get("protein_role"))
        checks["role_compatible"] = bool(all(roles1) and all(roles2) and roles1 == roles2)
        if not checks["role_compatible"]:
            return False, checks

    atoms1 = _atom_tokens(first.get("protein_atoms") or [], names_only=True)
    atoms2 = _atom_tokens(second.get("protein_atoms") or [], names_only=True)
    if atoms1 and atoms2:
        checks["atom_compatible"] = bool(atoms1 & atoms2)
        if not checks["atom_compatible"]:
            return False, checks

    ligand_atoms1 = _ligand_atom_tokens(first.get("ligand_atoms") or [])
    ligand_atoms2 = _ligand_atom_tokens(second.get("ligand_atoms") or [])
    if ligand_atoms1 and ligand_atoms2:
        checks["ligand_atom_compatible"] = bool(ligand_atoms1 & ligand_atoms2)
        if not checks["ligand_atom_compatible"]:
            return False, checks

    mediator1 = _mediator_feature(first.get("mediator"))
    mediator2 = _mediator_feature(second.get("mediator"))
    if mediator1[0] and mediator2[0]:
        type_match = mediator1[0] == mediator2[0]
        element_match = not (mediator1[1] and mediator2[1]) or mediator1[1] == mediator2[1]
        checks["mediator_compatible"] = type_match and element_match
        if not checks["mediator_compatible"]:
            return False, checks

    distance1 = _primary_distance(first.get("geometry") or {})
    distance2 = _primary_distance(second.get("geometry") or {})
    if distance1 is not None and distance2 is not None:
        tolerance = GEOMETRY_DISTANCE_TOLERANCE.get(interaction_type, 1.0)
        checks["geometry_compatible"] = abs(distance1 - distance2) <= tolerance
        if not checks["geometry_compatible"]:
            return False, checks

    for angle_key in ("angle", "donor_angle", "acceptor_angle"):
        angle1 = (first.get("geometry") or {}).get(angle_key)
        angle2 = (second.get("geometry") or {}).get(angle_key)
        if angle1 is not None and angle2 is not None:
            compatible = abs(float(angle1) - float(angle2)) <= ANGLE_TOLERANCE
            checks["geometry_compatible"] = (
                compatible if checks["geometry_compatible"] is None
                else checks["geometry_compatible"] and compatible
            )
            if not compatible:
                return False, checks

    # Atom identity and geometry must not both be absent.  A residue/type-only
    # collision is not independent chemical agreement.
    feature_evidence = (
        checks["atom_compatible"] is True
        or checks["ligand_atom_compatible"] is True
        or checks["geometry_compatible"] is True
    )
    checks["feature_evidence"] = feature_evidence
    return bool(feature_evidence), checks


def _consensus_groups(observations: list[dict]) -> list[dict]:
    plip = [obs for obs in observations if obs.get("tool") == "PLIP"]
    prolif = [obs for obs in observations if obs.get("tool") == "ProLIF"]
    candidates = []
    for left in plip:
        for right in prolif:
            matched, checks = observations_match(left, right)
            if matched:
                candidates.append((left, right, checks))

    # Prefer pairs with atom agreement, then smallest distance difference, and
    # consume each raw observation at most once.
    def rank(item: tuple[dict, dict, dict]) -> tuple[int, float]:
        left, right, checks = item
        d1 = _primary_distance(left.get("geometry") or {})
        d2 = _primary_distance(right.get("geometry") or {})
        return (
            0 if checks.get("atom_compatible") is True else 1,
            abs(d1 - d2) if d1 is not None and d2 is not None else 999.0,
        )

    groups = []
    consumed = set()
    for left, right, checks in sorted(candidates, key=rank):
        members = [left["observation_id"], right["observation_id"]]
        if consumed.intersection(members):
            continue
        consumed.update(members)
        stable = "|".join(sorted(members))
        groups.append({
            "consensus_id": f"cons_{hashlib.sha1(stable.encode()).hexdigest()[:12]}",
            "status": "detected",
            "interaction_type": left.get("interaction_type"),
            "residue": left.get("residue"),
            "ligand_role": left.get("ligand_role") or right.get("ligand_role"),
            "protein_role": left.get("protein_role") or right.get("protein_role"),
            "member_observation_ids": members,
            "tools": ["PLIP", "ProLIF"],
            "match_checks": checks,
            "independence_assessment": {
                "independent_experimental_evidence": False,
                "shared_dependencies": [
                    "same static coordinates",
                    "same protonation/microstate inputs",
                    "overlapping geometric interaction definitions",
                ],
                "interpretation": "detection concordance, not biological necessity",
            },
        })
    return groups


def _infer_tool_status(tool: str, result: dict | None) -> dict:
    critical_for_hard_gate = tool == "PLIP"
    if result is None:
        return {
            "status": "not_evaluable",
            "reason": "no tool result",
            "critical_for_hard_gate": critical_for_hard_gate,
        }
    if result.get("not_evaluable") is True:
        return {
            "status": "not_evaluable",
            "reason": result.get("error") or "tool input is not evaluable",
            "tool_version": result.get("tool_version"),
            "critical_for_hard_gate": critical_for_hard_gate,
        }
    if result.get("validation_failed") or result.get("success") is False:
        return {
            "status": "tool_failed",
            "reason": (result.get("diagnostics") or {}).get("error") or result.get("error"),
            "tool_version": result.get("tool_version"),
            "critical_for_hard_gate": critical_for_hard_gate,
        }
    count = int(result.get("num_interactions_total") or 0)
    if not count:
        count = len(result.get("interactions_detail") or [])
    if not count:
        count = sum(len(types or []) for types in (result.get("summary") or {}).values())
    return {
        "status": "detected" if count else "not_detected",
        "interaction_count": count,
        "tool_version": result.get("tool_version"),
        "critical_for_hard_gate": critical_for_hard_gate,
    }


def normalize_interaction_observations(
    plip_result: dict | None,
    prolif_result: dict | None,
    *,
    default_chain: str = "A",
    tool_statuses: dict[str, dict] | None = None,
    input_evidence: dict | None = None,
    supplemental_results: dict[str, dict | None] | None = None,
) -> dict:
    """Create the canonical Step-1 interaction observation contract."""
    supplemental_results = supplemental_results or {}
    supplemental_tools = {
        "arpeggio": "Arpeggio",
        "luna": "LUNA",
        "oddt": "ODDT",
    }
    observations = (
        _tool_observations("PLIP", plip_result, default_chain)
        + _tool_observations("ProLIF", prolif_result, default_chain)
    )
    for key, display_name in supplemental_tools.items():
        observations.extend(
            _tool_observations(
                display_name, supplemental_results.get(key), default_chain,
            )
        )
    statuses = {
        "plip": _infer_tool_status("PLIP", plip_result),
        "prolif": _infer_tool_status("ProLIF", prolif_result),
    }
    for key, display_name in supplemental_tools.items():
        if key in supplemental_results:
            statuses[key] = _infer_tool_status(
                display_name, supplemental_results.get(key),
            )
    for tool, status in (tool_statuses or {}).items():
        key = str(tool).lower()
        merged = dict(statuses.get(key) or {})
        merged.update(_json_value(status))
        merged["critical_for_hard_gate"] = key == "plip"
        statuses[key] = merged
    groups = _consensus_groups(observations)
    return {
        "schema_version": 2,
        "observation_semantics": "detector_evidence_not_design_requirement",
        "tool_statuses": statuses,
        "observations": observations,
        "consensus_groups": groups,
        "summary": {
            "detected_observation_count": len(observations),
            "strict_cross_tool_consensus_count": len(groups),
            "tool_failed_count": sum(
                status.get("status") == "tool_failed" for status in statuses.values()
            ),
            "not_evaluable_count": sum(
                status.get("status") == "not_evaluable" for status in statuses.values()
            ),
        },
        "input_evidence": _json_value(input_evidence or {}),
    }

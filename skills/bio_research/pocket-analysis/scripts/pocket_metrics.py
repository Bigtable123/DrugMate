"""Pocket metrics: volume, burial ratio, sub-pocket classification, flexibility, IoU."""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any, Iterable, Optional

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import (
    ONE_TO_THREE,
    canonical_residue_label,
    deduplicate_residue_records,
    extract_resnum as _shared_extract_resnum,
    legacy_residue_label,
    make_residue_record,
    normalize_residue as _shared_normalize_residue,
    parse_residue_record,
    residue_identity,
)
from ligand_site_utils import box_from_coordinates, load_ligand_heavy_atom_coords


def compute_pocket_volume(
    protein_file: Path,
    pocket_residues: list[str],
    chain: str = "A",
) -> float | None:
    """Estimate the residue-envelope volume, not the physical cavity volume.

    Args:
        protein_file: Path to protein PDB
        pocket_residues: List of residue identifiers (e.g., ["Leu83", "Glu81"])
        chain: Chain identifier

    Returns:
        Estimated volume in Angstrom^3, or None when too few atoms resolved to
        define an envelope at all.  A missing measurement is reported as missing
        rather than as a zero-volume pocket.
    """
    coords = _get_pocket_atom_coords(protein_file, pocket_residues, chain)
    if len(coords) < 4:
        return None

    from scipy.spatial import ConvexHull
    import numpy as np

    points = np.array(coords)
    hull = ConvexHull(points)
    return round(hull.volume, 1)


def compute_burial_ratio(
    protein_file: Path,
    pocket_residues: list[str],
    chain: str = "A",
) -> float | None:
    """Estimate burial ratio of pocket residues.

    Burial ratio = fraction of pocket atoms that are buried (not surface-exposed).
    Approximated by counting neighbors within 8A shell.

    Returns:
        Float between 0.0 (fully exposed) and 1.0 (fully buried), or None when no
        coordinates resolved.  0.0 means "measured as fully exposed", so absence
        gets its own value instead of borrowing that one.
    """
    pocket_coords = _get_pocket_atom_coords(protein_file, pocket_residues, chain)
    all_coords = _get_all_atom_coords(protein_file)

    if not pocket_coords or not all_coords:
        return None

    import numpy as np

    pocket_arr = np.array(pocket_coords)
    all_arr = np.array(all_coords)

    buried_count = 0
    total = len(pocket_arr)
    shell_radius = 8.0

    for patom in pocket_arr:
        dists = np.linalg.norm(all_arr - patom, axis=1)
        neighbors = np.sum((dists > 0.5) & (dists < shell_radius))
        if neighbors > 20:
            buried_count += 1

    return round(buried_count / total, 2) if total > 0 else None


def classify_subpockets(
    protein_file: Path,
    pocket_residues: list[str],
    chain: str = "A",
) -> list[dict]:
    """Classify pocket residues into sub-pockets by property.

    Groups residues into: polar, hydrophobic, charged, aromatic categories.

    Returns:
        List of sub-pocket dictionaries
    """
    HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
    POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
    POSITIVE = {"LYS", "ARG"}
    NEGATIVE = {"ASP", "GLU"}
    AROMATIC = {"PHE", "TYR", "TRP", "HIS"}

    hydrophobic_residues = []
    polar_residues = []
    positive_residues = []
    negative_residues = []
    aromatic_residues = []

    for res in pocket_residues:
        resname = _extract_resname(res)
        resname_upper = resname.upper() if resname else ""

        three_letter = _one_to_three(resname_upper) if len(resname_upper) <= 1 else resname_upper[:3].upper()

        if three_letter in HYDROPHOBIC:
            hydrophobic_residues.append(res)
        if three_letter in POLAR:
            polar_residues.append(res)
        if three_letter in POSITIVE:
            positive_residues.append(res)
        if three_letter in NEGATIVE:
            negative_residues.append(res)
        if three_letter in AROMATIC:
            aromatic_residues.append(res)

    subpockets = []
    if hydrophobic_residues:
        subpockets.append({
            "id": "hydrophobic_region",
            "residues": hydrophobic_residues,
            "property": "hydrophobic",
        })
    if polar_residues:
        subpockets.append({
            "id": "polar_region",
            "residues": polar_residues,
            "property": "polar",
        })
    if positive_residues or negative_residues:
        subpockets.append({
            "id": "charged_region",
            "residues": positive_residues + negative_residues,
            "property": "charged",
        })
    if aromatic_residues:
        subpockets.append({
            "id": "aromatic_region",
            "residues": aromatic_residues,
            "property": "aromatic",
        })

    return subpockets


def compute_flexibility_zones(
    protein_file: Path,
    pocket_residues: list[str],
    chain: str = "A",
    bfactor_threshold: float = 30.0,
) -> list[dict]:
    """Identify flexible zones in the pocket based on B-factors.

    Returns:
        List of flexibility zone dictionaries
    """
    residue_bfactors = _get_residue_bfactors(protein_file, pocket_residues, chain)

    flexible = []
    for res, bfactor in residue_bfactors.items():
        if bfactor > bfactor_threshold:
            flexible.append(res)

    if not flexible:
        return []

    zones = _group_consecutive_residues(flexible)

    return [
        {
            "residues": zone,
            "note": f"B-factor > {bfactor_threshold}",
        }
        for zone in zones
    ]


def compute_iou(detected_residues: Iterable[Any], known_residues: Iterable[Any]) -> float:
    """Compute chain/insertion-code-safe residue Intersection over Union.

    Legacy labels remain supported.  Labels without an explicit chain compare in
    the historical chain-agnostic namespace; structured records and canonical
    labels compare by ``(chain, resnum, icode, resname)``.  Residues explicitly
    marked ``identity_ambiguous`` are excluded rather than creating false
    agreement.
    """
    detected_normalized = _residue_identity_set(detected_residues)
    known_normalized = _residue_identity_set(known_residues)
    intersection = detected_normalized & known_normalized
    union = detected_normalized | known_normalized
    return len(intersection) / len(union) if union else 0.0


def compute_consensus_pocket(
    fpocket_result: Optional[dict],
    p2rank_result: Optional[dict],
    ligand_file: Optional[Path],
    protein_file: Path,
    chain: str = "A",
    *,
    pykvfinder_result: Optional[dict] = None,
    deeppocket_result: Optional[dict] = None,
    ligand_site_result: Optional[dict] = None,
    top_k: int = 5,
    center_match_threshold: float = 5.0,
    residue_iou_threshold: float = 0.3,
) -> dict:
    """Build a Top-K cross-tool consensus with optional ligand-site evidence.

    The first five positional arguments are retained for existing callers.
    Additional detector results and all tuning controls are keyword-only so
    extending the interface cannot silently reorder legacy calls.
    """
    try:
        top_k = max(1, int(top_k))
    except (TypeError, ValueError):
        top_k = 5
    center_match_threshold = max(0.0, float(center_match_threshold))
    residue_iou_threshold = max(0.0, min(1.0, float(residue_iou_threshold)))

    detector_results = {
        "fpocket": fpocket_result,
        "p2rank": p2rank_result,
        "pykvfinder": pykvfinder_result,
        "deeppocket": deeppocket_result,
    }
    nodes = _detector_nodes(detector_results, top_k=top_k, default_chain=chain)
    candidates = _build_detector_candidates(
        nodes,
        center_match_threshold=center_match_threshold,
        residue_iou_threshold=residue_iou_threshold,
        protein_file=protein_file,
        chain=chain,
    )

    ligand_site = ligand_site_result

    if ligand_site and ligand_site.get("center") is not None:
        selected = _ligand_first_consensus(
            ligand_site,
            nodes,
            candidates,
            center_match_threshold=center_match_threshold,
            residue_iou_threshold=residue_iou_threshold,
            default_chain=chain,
        )
        return selected

    if candidates:
        selected = dict(candidates[0])
        selected["candidate_pockets"] = candidates
        selected["selection_method"] = "top_k_cross_tool_rank_aggregation"
        return selected

    return {
        "residues": [],
        "residue_records": [],
        "center": None,
        "box_center": None,
        "size": None,
        "iou_fpocket_p2rank": None,
        "source": "none",
        "candidate_pockets": [],
        "selection_method": "no_detected_pockets",
    }


def build_pocket_box(
    protein_file: Path | str,
    pocket_residues: Iterable[Any] | None = None,
    chain: str = "A",
    *,
    residue_records: Iterable[Any] | None = None,
    padding: float = 4.0,
    min_size: float = 12.0,
) -> dict:
    """Build a docking box from the atoms of selected pocket residues.

    ``residue_records`` is preferred; ``pocket_residues`` preserves the legacy
    calling convention.  Each axis is clamped to ``min_size`` so a thin residue
    envelope cannot produce a degenerate detector-derived docking box.
    """
    values = list(residue_records or pocket_residues or [])
    records = deduplicate_residue_records(values, default_chain=chain)
    coords = _get_pocket_atom_coords(Path(protein_file), records, chain)
    if not coords:
        return {"center": None, "size": None, "bounds": None, "padding": float(padding)}
    center, size, bounds = box_from_coordinates(coords, padding=max(0.0, float(padding)))
    minimum = max(0.0, float(min_size))
    size = [round(max(axis_size, minimum), 4) for axis_size in size]
    return {
        "center": center,
        "size": size,
        "bounds": bounds,
        "padding": float(padding),
        "min_size": minimum,
    }


_TOOL_PRIORITY = {
    "p2rank": 0,
    "deeppocket": 1,
    "fpocket": 2,
    "pykvfinder": 3,
}


def _independence_group(tool: str) -> str:
    """Return the independent-evidence group for one detector.

    DeepPocket uses fpocket cavities as its proposal set before CNN ranking and
    segmentation.  Agreement between those two tools is useful, but it must not
    be counted as two independent cavity-detector votes.
    """
    return "fpocket_family" if tool in {"fpocket", "deeppocket"} else tool


def _independent_support(tools: Iterable[str]) -> list[str]:
    return sorted({_independence_group(tool) for tool in tools})


def _detector_nodes(
    detector_results: dict[str, Optional[dict]],
    *,
    top_k: int,
    default_chain: str,
) -> list[dict]:
    nodes: list[dict] = []
    for tool, result in detector_results.items():
        if not result:
            continue
        pockets = list(result.get("all_pockets") or [])
        if not pockets and result.get("top_pocket"):
            pockets = [result["top_pocket"]]
        pockets.sort(key=lambda item: _rank_value(item.get("rank")))
        for index, pocket in enumerate(pockets[:top_k], start=1):
            rank = _rank_value(pocket.get("rank"), fallback=index)
            raw_records = pocket.get("residue_records")
            if raw_records is None:
                raw_records = pocket.get("residues") or []
            records = deduplicate_residue_records(raw_records, default_chain=default_chain)
            if default_chain:
                records = [
                    record for record in records
                    if str(record.get("chain", "")).strip() == str(default_chain).strip()
                ]
            # A detector may report pockets from every chain.  Once a chain is
            # requested, a pocket with no residues on that chain is not a site
            # candidate for this route; retaining only its center would create
            # false support and can select the wrong chain.
            if not records:
                continue
            nodes.append({
                "tool": tool,
                "rank": rank,
                "pocket": pocket,
                "center": _valid_center(pocket.get("center")),
                "residue_records": records,
            })
    return nodes


def _build_detector_candidates(
    nodes: list[dict],
    *,
    center_match_threshold: float,
    residue_iou_threshold: float,
    protein_file: Path,
    chain: str,
) -> list[dict]:
    """Create anchor-based candidates with at most one pocket per tool."""
    candidates: list[dict] = []
    seen: set[tuple] = set()
    tools = tuple(_TOOL_PRIORITY)

    for anchor in nodes:
        members = [anchor]
        for tool in tools:
            if tool == anchor["tool"]:
                continue
            choices = [node for node in nodes if node["tool"] == tool]
            match = _best_matching_node(
                anchor,
                choices,
                center_match_threshold=center_match_threshold,
                residue_iou_threshold=residue_iou_threshold,
            )
            if match is not None:
                members.append(match)

        member_key = tuple(sorted(_node_key(member) for member in members))
        if member_key in seen:
            continue
        seen.add(member_key)
        candidates.append(
            _candidate_from_members(members, protein_file=protein_file, chain=chain)
        )

    candidates.sort(key=_candidate_sort_key)
    for index, candidate in enumerate(candidates, start=1):
        candidate["candidate_id"] = f"P{index:03d}"
    return candidates


def _candidate_from_members(members: list[dict], *, protein_file: Path, chain: str) -> dict:
    members = sorted(
        members,
        key=lambda item: (_TOOL_PRIORITY.get(item["tool"], 99), item["rank"]),
    )
    records = deduplicate_residue_records(
        record for member in members for record in member["residue_records"]
    )
    centers = [member["center"] for member in members if member["center"] is not None]
    center = _mean_center(centers)
    docking_box = build_pocket_box(
        protein_file,
        residue_records=records,
        chain=chain,
    )
    family_rank_scores: dict[str, float] = {}
    for member in members:
        family = _independence_group(member["tool"])
        family_rank_scores[family] = max(
            family_rank_scores.get(family, 0.0),
            1.0 / max(1, member["rank"]),
        )
    rank_score = sum(family_rank_scores.values())
    tools = [member["tool"] for member in members]
    support_groups = _independent_support(tools)
    support_count = len(support_groups)
    if support_count > 1:
        source = "consensus"
    elif len(set(tools)) > 1:
        source = f"{support_groups[0]}_dependent_consensus"
    else:
        source = f"{tools[0]}_only"

    candidate: dict[str, Any] = {
        "residues": [legacy_residue_label(record) for record in records],
        "residue_records": records,
        "center": center,
        "box_center": docking_box.get("center"),
        "size": docking_box.get("size"),
        "box": docking_box,
        "source": source,
        "support_count": support_count,
        "raw_tool_support_count": len(set(tools)),
        "independent_support_groups": support_groups,
        "detector_support": tools,
        "rank_score": round(rank_score, 6),
        "members": [_member_summary(member) for member in members],
        "pairwise_agreement": _pairwise_agreement(members),
        "iou_fpocket_p2rank": None,
    }
    for member in members:
        candidate[f"selected_{member['tool']}"] = member["pocket"]

    # Candidate residues are a geometry union used to build a sufficiently
    # large docking box.  Interaction promotion is stricter: a residue must be
    # reported by at least two independent geometry detectors in this matched
    # candidate, not merely belong to that union.
    detector_residue_sets = {
        member["tool"]: {
            identity
            for record in member.get("residue_records") or []
            if (identity := residue_identity(record)) is not None
        }
        for member in members
    }
    interaction_support = []
    interaction_consensus_records = []
    for record in records:
        identity = residue_identity(record)
        supporting_tools = sorted(
            tool for tool, identities in detector_residue_sets.items()
            if identity is not None and identity in identities
        )
        supporting_groups = _independent_support(supporting_tools)
        interaction_support.append({
            "residue_record": record,
            "residue": legacy_residue_label(record),
            "support_count": len(supporting_groups),
            "raw_tool_support_count": len(supporting_tools),
            "detector_support": supporting_tools,
            "independent_support_groups": supporting_groups,
        })
        if len(supporting_groups) >= 2:
            interaction_consensus_records.append(record)
    candidate["interaction_consensus_residue_records"] = interaction_consensus_records
    candidate["interaction_consensus_residues"] = [
        legacy_residue_label(record) for record in interaction_consensus_records
    ]
    candidate["interaction_pocket_support"] = interaction_support
    candidate["interaction_consensus_min_detector_support"] = 2

    by_tool = {member["tool"]: member for member in members}
    if "fpocket" in by_tool and "p2rank" in by_tool:
        candidate["iou_fpocket_p2rank"] = round(
            compute_iou(
                by_tool["fpocket"]["residue_records"],
                by_tool["p2rank"]["residue_records"],
            ),
            3,
        )
        intersection, union = _intersection_and_union_records(
            by_tool["fpocket"]["residue_records"],
            by_tool["p2rank"]["residue_records"],
        )
        candidate["intersection_residues"] = [
            legacy_residue_label(record) for record in intersection
        ]
        candidate["union_residues"] = [legacy_residue_label(record) for record in union]
    return candidate


def _ligand_first_consensus(
    ligand_site: dict,
    nodes: list[dict],
    detector_candidates: list[dict],
    *,
    center_match_threshold: float,
    residue_iou_threshold: float,
    default_chain: str,
) -> dict:
    raw_records = ligand_site.get("residue_records")
    if raw_records is None:
        raw_records = ligand_site.get("residues") or []
    records = deduplicate_residue_records(raw_records, default_chain=default_chain)
    ligand_node = {
        "tool": "ligand_site",
        "rank": 1,
        "pocket": ligand_site,
        "center": _valid_center(ligand_site.get("center")),
        "residue_records": records,
    }

    matched: list[dict] = []
    for tool in _TOOL_PRIORITY:
        choices = [node for node in nodes if node["tool"] == tool]
        node = _best_matching_node(
            ligand_node,
            choices,
            center_match_threshold=center_match_threshold,
            residue_iou_threshold=residue_iou_threshold,
        )
        if node is not None:
            matched.append(node)

    matched_tools = [node["tool"] for node in matched]
    matched_support_groups = _independent_support(matched_tools)
    primary: dict[str, Any] = {
        "candidate_id": "P000",
        "residues": [legacy_residue_label(record) for record in records],
        "residue_records": records,
        # center is deliberately the heavy-atom centroid.  box_center may differ.
        "center": ligand_site.get("center"),
        "box_center": ligand_site.get("box_center", ligand_site.get("center")),
        "size": ligand_site.get("size"),
        "box": {
            "center": ligand_site.get("box_center", ligand_site.get("center")),
            "size": ligand_site.get("size"),
            "bounds": ligand_site.get("bounds"),
            "padding": ligand_site.get("box_padding"),
        },
        "source": "co_crystal_ligand",
        "selection_method": "ligand_defined_primary",
        "support_count": 1 + len(matched_support_groups),
        "raw_tool_support_count": 1 + len(set(matched_tools)),
        "independent_support_groups": ["ligand_site", *matched_support_groups],
        "detector_support": matched_tools,
        "members": [
            {
                "tool": "ligand_site",
                "rank": 1,
                "center": ligand_site.get("center"),
                "residue_count": len(records),
            },
            *[_member_summary(node) for node in matched],
        ],
        "ligand_site": ligand_site,
        "iou_fpocket_p2rank": None,
    }
    by_tool = {node["tool"]: node for node in matched}
    for node in matched:
        primary[f"selected_{node['tool']}"] = node["pocket"]

    # The ligand-defined 5 Å shell is pocket geometry, not multi-detector
    # evidence.  Interaction promotion may only use residues supported by at
    # least two independent geometry detectors matched to this ligand site.
    detector_residue_sets = {
        node["tool"]: {
            identity
            for record in node.get("residue_records") or []
            if (identity := residue_identity(record)) is not None
        }
        for node in matched
    }
    interaction_support = []
    interaction_consensus_records = []
    for record in records:
        identity = residue_identity(record)
        supporting_tools = sorted(
            tool for tool, identities in detector_residue_sets.items()
            if identity is not None and identity in identities
        )
        supporting_groups = _independent_support(supporting_tools)
        interaction_support.append({
            "residue_record": record,
            "residue": legacy_residue_label(record),
            "support_count": len(supporting_groups),
            "raw_tool_support_count": len(supporting_tools),
            "detector_support": supporting_tools,
            "independent_support_groups": supporting_groups,
        })
        if len(supporting_groups) >= 2:
            interaction_consensus_records.append(record)
    primary["interaction_consensus_residue_records"] = interaction_consensus_records
    primary["interaction_consensus_residues"] = [
        legacy_residue_label(record) for record in interaction_consensus_records
    ]
    primary["interaction_pocket_support"] = interaction_support
    primary["interaction_consensus_min_detector_support"] = 2
    if "fpocket" in by_tool and "p2rank" in by_tool:
        primary["iou_fpocket_p2rank"] = round(
            compute_iou(
                by_tool["fpocket"]["residue_records"],
                by_tool["p2rank"]["residue_records"],
            ),
            3,
        )
        intersection, union = _intersection_and_union_records(
            by_tool["fpocket"]["residue_records"],
            by_tool["p2rank"]["residue_records"],
        )
        primary["intersection_residues"] = [
            legacy_residue_label(record) for record in intersection
        ]
        primary["union_residues"] = [legacy_residue_label(record) for record in union]

    primary_summary = {key: value for key, value in primary.items() if key != "ligand_site"}
    primary["candidate_pockets"] = [primary_summary, *detector_candidates]
    primary["alternative_candidates"] = detector_candidates
    return primary


def _best_matching_node(
    anchor: dict,
    choices: list[dict],
    *,
    center_match_threshold: float,
    residue_iou_threshold: float,
) -> dict | None:
    accepted: list[tuple[tuple[float, float, float], dict]] = []
    for choice in choices:
        iou, distance = _node_agreement(anchor, choice)
        center_match = distance is not None and distance <= center_match_threshold
        residue_match = iou >= residue_iou_threshold and iou > 0.0
        # Both spatial proximity and residue agreement are required. Neither
        # nearby disjoint sites nor distant overlapping predictions may match.
        if not (center_match and residue_match):
            continue
        # Prefer residue agreement, then spatial proximity, then native rank.
        distance_score = -distance if distance is not None else float("-inf")
        accepted.append(((iou, distance_score, -choice["rank"]), choice))
    if not accepted:
        return None
    accepted.sort(key=lambda item: item[0], reverse=True)
    return accepted[0][1]


def _node_agreement(left: dict, right: dict) -> tuple[float, float | None]:
    iou = compute_iou(left.get("residue_records") or [], right.get("residue_records") or [])
    distance = _center_distance(left.get("center"), right.get("center"))
    return iou, distance


def _pairwise_agreement(members: list[dict]) -> list[dict]:
    agreement = []
    for left_index, left in enumerate(members):
        for right in members[left_index + 1:]:
            iou, distance = _node_agreement(left, right)
            agreement.append({
                "tools": [left["tool"], right["tool"]],
                "residue_iou": round(iou, 4),
                "center_distance": round(distance, 4) if distance is not None else None,
            })
    return agreement


def _member_summary(member: dict) -> dict:
    pocket = member["pocket"]
    summary = {
        "tool": member["tool"],
        "rank": member["rank"],
        "center": member.get("center"),
        "residue_count": len(member.get("residue_records") or []),
    }
    for key in (
        "name", "cavity_id", "score", "probability", "volume",
        "original_fpocket_rank", "candidate_dependency",
        "segmentation_residue_count",
    ):
        if pocket.get(key) is not None:
            summary[key] = pocket[key]
    return summary


def _candidate_sort_key(candidate: dict) -> tuple:
    tools = candidate.get("detector_support") or []
    best_priority = min((_TOOL_PRIORITY.get(tool, 99) for tool in tools), default=99)
    best_rank = min(
        (member.get("rank", 10**9) for member in candidate.get("members") or []),
        default=10**9,
    )
    return (
        -int(candidate.get("support_count") or 0),
        -float(candidate.get("rank_score") or 0.0),
        best_priority,
        best_rank,
    )


def _node_key(node: dict) -> tuple:
    pocket = node.get("pocket") or {}
    discriminator = pocket.get("name", pocket.get("cavity_id", node.get("rank")))
    return node.get("tool"), node.get("rank"), str(discriminator)


def _rank_value(value: Any, fallback: int = 10**9) -> int:
    try:
        rank = int(value)
        return rank if rank > 0 else fallback
    except (TypeError, ValueError):
        return fallback


def _valid_center(value: Any) -> list[float] | None:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        return None
    try:
        center = [float(axis) for axis in value]
    except (TypeError, ValueError):
        return None
    if any(axis != axis or axis in (float("inf"), float("-inf")) for axis in center):
        return None
    return center


def _center_distance(left: Any, right: Any) -> float | None:
    left_center = _valid_center(left)
    right_center = _valid_center(right)
    if left_center is None or right_center is None:
        return None
    return sum((left_center[index] - right_center[index]) ** 2 for index in range(3)) ** 0.5


def _mean_center(centers: list[list[float]]) -> list[float] | None:
    if not centers:
        return None
    return [
        round(sum(center[axis] for center in centers) / len(centers), 4)
        for axis in range(3)
    ]


def _residue_identity_set(values: Iterable[Any]) -> set[tuple[Any, ...]]:
    identities: set[tuple[Any, ...]] = set()
    for value in values or []:
        identity = residue_identity(value)
        if identity is not None:
            identities.add(identity)
            continue
        # Preserve historical comparison for unparsable legacy labels, but not
        # for records explicitly marked ambiguous.
        if isinstance(value, dict) and value.get("identity_ambiguous"):
            continue
        if not isinstance(value, dict):
            normalized = _shared_normalize_residue(str(value))
            if normalized:
                identities.add(("legacy", normalized))
    return identities


def _intersection_and_union_records(
    left: Iterable[Any], right: Iterable[Any],
) -> tuple[list[dict], list[dict]]:
    left_records = deduplicate_residue_records(left)
    right_records = deduplicate_residue_records(right)
    left_map = {
        identity: record
        for record in left_records
        if (identity := residue_identity(record)) is not None
    }
    right_map = {
        identity: record
        for record in right_records
        if (identity := residue_identity(record)) is not None
    }
    intersection = deduplicate_residue_records(
        left_map[identity] for identity in left_map.keys() & right_map.keys()
    )
    union = deduplicate_residue_records([*left_records, *right_records])
    return intersection, union


# ===== Helper functions =====

def _get_pocket_atom_coords(
    protein_file: Path, pocket_residues: Iterable[Any], chain: str
) -> list[list[float]]:
    """Extract heavy-atom coordinates using exact chain/number/iCode identity."""
    records = deduplicate_residue_records(pocket_residues, default_chain=chain)
    exact_targets = {
        identity
        for record in records
        if (identity := residue_identity(record, allow_ambiguous=True)) is not None
    }
    position_targets = {
        residue_identity(record, include_resname=False, allow_ambiguous=True)
        for record in records
        if not record.get("resname")
    }

    coords = []
    with open(protein_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() or atom_name[:1]
            if element.upper() == "H":
                continue
            record = _pdb_residue_record(line)
            if record is None:
                continue
            exact = residue_identity(record, allow_ambiguous=True)
            position = residue_identity(record, include_resname=False, allow_ambiguous=True)
            if exact not in exact_targets and position not in position_targets:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue
    return coords


def _get_all_atom_coords(protein_file: Path) -> list[list[float]]:
    """Extract all heavy-atom coordinates from PDB (for burial proxy)."""
    coords = []
    with open(protein_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() or atom_name[:1]
            if element.upper() == "H":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue
    return coords


def _get_residue_bfactors(
    protein_file: Path, pocket_residues: Iterable[Any], chain: str
) -> dict:
    """Get average B-factor per residue for pocket residues."""
    records = deduplicate_residue_records(pocket_residues, default_chain=chain)
    target_identities = {
        identity
        for record in records
        if (identity := residue_identity(record, allow_ambiguous=True)) is not None
    }

    bfactor_sums = {}
    bfactor_counts = {}

    with open(protein_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            record = _pdb_residue_record(line)
            if record is None:
                continue
            try:
                bfactor = float(line[60:66].strip())
            except (ValueError, IndexError):
                continue
            if residue_identity(record, allow_ambiguous=True) in target_identities:
                key = canonical_residue_label(record)
                bfactor_sums[key] = bfactor_sums.get(key, 0.0) + bfactor
                bfactor_counts[key] = bfactor_counts.get(key, 0) + 1

    return {
        k: round(bfactor_sums[k] / bfactor_counts[k], 1)
        for k in bfactor_sums
        if bfactor_counts[k] > 0
    }


def _select_relevant_pocket(result: Optional[dict], ligand_center: list[float] | None) -> dict | None:
    """Select the pocket nearest the co-crystal ligand; otherwise use tool rank 1."""
    if not result:
        return None
    pockets = result.get("all_pockets") or []
    if not pockets and result.get("top_pocket"):
        pockets = [result["top_pocket"]]
    if not pockets:
        return None
    if ligand_center is None:
        return min(pockets, key=lambda p: p.get("rank") or 10**9)

    def distance_sq(pocket: dict) -> float:
        center = pocket.get("center")
        if not center or len(center) != 3:
            return float("inf")
        return sum((float(center[i]) - ligand_center[i]) ** 2 for i in range(3))

    selected = min(pockets, key=distance_sq)
    if distance_sq(selected) == float("inf"):
        return min(pockets, key=lambda p: p.get("rank") or 10**9)
    return selected


def _get_ligand_center(ligand_file: Path | None) -> list[float] | None:
    if not ligand_file or not Path(ligand_file).exists():
        return None
    try:
        coords = load_ligand_heavy_atom_coords(Path(ligand_file))
    except Exception:
        return None
    if not coords:
        return None
    return [sum(point[i] for point in coords) / len(coords) for i in range(3)]


def _extract_resnum(residue_str: Any) -> Optional[int]:
    """Extract residue number from string like 'Leu83', 'LEU83', 'A_LEU83'."""
    record = parse_residue_record(residue_str)
    if record is not None:
        try:
            return int(record["resnum"])
        except (TypeError, ValueError):
            return None
    return _shared_extract_resnum(residue_str)


def _extract_resname(residue_str: Any) -> str:
    """Extract residue name from string like 'Leu83'."""
    record = parse_residue_record(residue_str)
    if record is not None and record.get("resname"):
        return str(record["resname"])
    normalized = _shared_normalize_residue(residue_str)
    match = re.match(r"([A-Za-z]+)", normalized)
    return match.group(1) if match else ""


def _normalize_residue(res: str) -> str:
    """Normalize residue name for comparison (e.g., 'Leu83' -> 'LEU83')."""
    return _shared_normalize_residue(res)


def _one_to_three(one_letter: str) -> str:
    """Convert one-letter amino acid code to three-letter."""
    return ONE_TO_THREE.get(one_letter.upper(), one_letter)


def _group_consecutive_residues(residues: list[str]) -> list[list[str]]:
    """Group residues into consecutive stretches."""
    if not residues:
        return []

    numbered = []
    for r in residues:
        num = _extract_resnum(r)
        if num is not None:
            numbered.append((num, r))

    numbered.sort(key=lambda x: x[0])

    groups = []
    current_group = [numbered[0][1]] if numbered else []
    for i in range(1, len(numbered)):
        if numbered[i][0] - numbered[i-1][0] <= 2:
            current_group.append(numbered[i][1])
        else:
            groups.append(current_group)
            current_group = [numbered[i][1]]
    if current_group:
        groups.append(current_group)

    return groups


def _pdb_residue_record(line: str) -> dict | None:
    try:
        return make_residue_record(
            line[17:20].strip(), int(line[22:26].strip()),
            line[21:22].strip(), line[26:27].strip(),
        )
    except (ValueError, IndexError):
        return None

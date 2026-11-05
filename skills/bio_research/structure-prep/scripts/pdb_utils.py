"""PDB utilities: download, parse, search, extract receptor/ligand."""

from __future__ import annotations

import json
import logging
import urllib.request
from pathlib import Path
from typing import Any

from rcsb_http import RCSBNotFoundError, request_bytes, request_json

logger = logging.getLogger(__name__)

RCSB_DATA_API = "https://data.rcsb.org/rest/v1"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILES = "https://files.rcsb.org/download"

COMMON_IONS = {"NA", "CL", "MG", "ZN", "CA", "FE", "MN", "K", "NI", "CU", "CO", "CD", "SO4", "PO4", "GOL", "EDO", "ACE", "NH2"}
DEFAULT_PRESERVED_HETERO = {"MG", "ZN", "CA", "FE", "MN", "K", "NI", "CU", "CO", "CD"}
WATER_CODES = {"HOH", "WAT", "H2O"}


def search_pdb_by_protein_name(protein_name: str, max_results: int = 5) -> list[dict[str, Any]]:
    query_json = json.dumps({
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {"attribute": "struct.title", "operator": "contains_words", "value": protein_name}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.experimental_method", "operator": "exact_match", "value": "X-ray"}},
            ],
        },
        "request_options": {"paginate": {"start": 0, "rows": max_results}, "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]},
        "return_type": "entry",
    }).encode()
    req = urllib.request.Request(RCSB_SEARCH_API, data=query_json, headers={"Content-Type": "application/json", "Accept": "application/json"})
    data = request_json(req)
    return [
        {"pdb_id": hit["identifier"], "score": hit.get("score", 0)}
        for hit in data.get("result_set", [])
    ]


def fetch_pdb_metadata(pdb_id: str) -> dict[str, Any]:
    url = f"{RCSB_DATA_API}/core/entry/{pdb_id.upper()}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    return request_json(req)


def fetch_polymer_entities(pdb_id: str) -> list[dict[str, Any]]:
    entities = []
    for entity_id in range(1, 20):
        url = f"{RCSB_DATA_API}/core/polymer_entity/{pdb_id.upper()}/{entity_id}"
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            entities.append(request_json(req))
        except RCSBNotFoundError:
            break
    return entities


def fetch_all_nonpolymer_entities(pdb_id: str, max_entities: int = 20) -> list[dict[str, Any]]:
    entities = []
    for entity_id in range(1, max_entities + 1):
        url = f"{RCSB_DATA_API}/core/nonpolymer_entity/{pdb_id.upper()}/{entity_id}"
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            entities.append(request_json(req))
        except RCSBNotFoundError:
            break
    return entities


def download_pdb_file(pdb_id: str, output_dir: Path, fmt: str = "pdb") -> Path:
    suffix = "pdb" if fmt == "pdb" else "cif"
    url = f"{RCSB_FILES}/{pdb_id.upper()}.{suffix}"
    output_path = output_dir / f"{pdb_id.upper()}.{suffix}"
    logger.info(f"Downloading {url} -> {output_path}")
    req = urllib.request.Request(
        url,
        headers={"Accept": "chemical/x-pdb,chemical/x-mmcif,text/plain,*/*"},
    )
    output_path.write_bytes(request_bytes(req))
    return output_path


def parse_pdb_atoms(pdb_file: Path) -> list[dict[str, Any]]:
    atoms = []
    with open(pdb_file) as f:
        for order, line in enumerate(f):
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            atoms.append({
                "record": record,
                "serial": int(line[6:11].strip()),
                "name": line[12:16].strip(),
                "alt_loc": line[16].strip() if len(line) > 16 else "",
                "res_name": line[17:20].strip(),
                "chain_id": line[21].strip(),
                "res_seq": int(line[22:26].strip()),
                "i_code": line[26].strip() if len(line) > 26 else "",
                "x": float(line[30:38].strip()),
                "y": float(line[38:46].strip()),
                "z": float(line[46:54].strip()),
                "occupancy": _parse_float(line[54:60].strip(), 0.0),
                "b_factor": _parse_float(line[60:66].strip(), 0.0),
                "element": line[76:78].strip() if len(line) > 76 else "",
                "raw_line": line,
                "_order": order,
            })
    return _select_altloc_atoms(atoms)


def _parse_float(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _select_altloc_atoms(atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep one atom per alternate-location site, preserving PDB order."""
    grouped: dict[tuple, list[dict[str, Any]]] = {}
    for atom in atoms:
        key = (
            atom["record"],
            atom["chain_id"],
            atom["res_seq"],
            atom.get("i_code", ""),
            atom["res_name"],
            atom["name"],
        )
        grouped.setdefault(key, []).append(atom)

    selected = []
    for group in grouped.values():
        selected.append(max(group, key=_altloc_rank))
    return sorted(selected, key=lambda atom: atom["_order"])


def _altloc_rank(atom: dict[str, Any]) -> tuple[int, int, float]:
    alt_loc = atom.get("alt_loc", "")
    blank = 1 if not alt_loc else 0
    preferred = 1 if alt_loc == "A" else 0
    return (blank, preferred, float(atom.get("occupancy") or 0.0))


def get_chain_info(
    atoms: list[dict],
    polymer_hetero_keys: set[tuple[str, str, int, str]] | None = None,
) -> dict[str, dict]:
    chains: dict[str, set] = {}
    polymer_hetero_keys = polymer_hetero_keys or set()
    for a in atoms:
        key = (
            a["res_name"], a["chain_id"], a["res_seq"], a.get("i_code", ""),
        )
        if a["record"] == "ATOM" or key in polymer_hetero_keys:
            chains.setdefault(a["chain_id"], set()).add(a["res_seq"])
    return {cid: {"chain_id": cid, "residue_count": len(res), "first_residue": min(res), "last_residue": max(res)} for cid, res in chains.items()}


def compute_distance(a1: dict, a2: dict) -> float:
    return ((a1["x"] - a2["x"])**2 + (a1["y"] - a2["y"])**2 + (a1["z"] - a2["z"])**2) ** 0.5


def find_conserved_waters(
    atoms: list[dict],
    ligand_id: str,
    distance_threshold: float = 3.5,
    ligand_key: tuple[str, str, int, str] | None = None,
) -> list[dict]:
    if ligand_key:
        ligand_atoms = [a for a in atoms if _atom_instance_key(a) == ligand_key]
    else:
        ligand_atoms = [a for a in atoms if a["record"] == "HETATM" and a["res_name"] == ligand_id]
    water_atoms = [a for a in atoms if a["record"] == "HETATM" and a["res_name"] in WATER_CODES]
    if not ligand_atoms or not water_atoms:
        return []

    water_residues: dict[tuple, list[dict]] = {}
    for w in water_atoms:
        water_residues.setdefault((w["chain_id"], w["res_seq"]), []).append(w)

    conserved = []
    for key, w_atoms in water_residues.items():
        keep = False
        for w in w_atoms:
            for lig in ligand_atoms:
                if compute_distance(w, lig) <= distance_threshold:
                    keep = True
                    break
            if keep:
                break
        if keep:
            conserved.extend(w_atoms)
    return conserved


def select_binding_site_waters(
    atoms: list[dict],
    ligand_key: tuple[str, str, int, str] | None,
    *,
    ligand_distance: float = 3.5,
    protein_distance: float = 3.5,
    min_occupancy: float = 0.5,
    max_b_factor: float = 60.0,
    include_residues: list[str] | None = None,
    exclude_residues: list[str] | None = None,
) -> dict[str, Any]:
    """Select likely bridging waters with explicit, auditable heuristics.

    "Conserved" cannot be established from one structure alone.  This selector
    therefore requires proximity to both ligand and protein plus crystal-quality
    thresholds, and labels the result ``binding_site_water`` in its report.
    """
    ligand_atoms = [a for a in atoms if ligand_key and _atom_instance_key(a) == ligand_key]
    protein_atoms = [a for a in atoms if a["record"] == "ATOM"]
    waters: dict[tuple[str, int], list[dict]] = {}
    for atom in atoms:
        if atom["record"] == "HETATM" and atom["res_name"] in WATER_CODES:
            waters.setdefault((atom["chain_id"], atom["res_seq"]), []).append(atom)

    manual_include = {_parse_residue_selector(value) for value in (include_residues or [])}
    manual_exclude = {_parse_residue_selector(value) for value in (exclude_residues or [])}
    manual_include.discard(None)
    manual_exclude.discard(None)
    selected_atoms: list[dict] = []
    entries = []
    for key, water_atoms in sorted(waters.items()):
        representative = water_atoms[0]
        min_ligand = min(
            (compute_distance(w, lig) for w in water_atoms for lig in ligand_atoms),
            default=None,
        )
        min_protein = min(
            (compute_distance(w, atom) for w in water_atoms for atom in protein_atoms),
            default=None,
        )
        occupancy = max(float(w.get("occupancy") or 0.0) for w in water_atoms)
        b_factor = min(float(w.get("b_factor") or 0.0) for w in water_atoms)
        automatic = bool(
            ligand_atoms
            and min_ligand is not None and min_ligand <= ligand_distance
            and min_protein is not None and min_protein <= protein_distance
            and occupancy >= min_occupancy
            and b_factor <= max_b_factor
        )
        kept = (automatic or key in manual_include) and key not in manual_exclude
        reason = "manual_exclude" if key in manual_exclude else (
            "manual_include" if key in manual_include else (
                "ligand_protein_bridge" if automatic else "failed_heuristic"
            )
        )
        if kept:
            selected_atoms.extend(water_atoms)
        entries.append({
            "chain_id": key[0], "res_seq": key[1],
            "occupancy": occupancy, "b_factor": b_factor,
            "min_ligand_distance": round(min_ligand, 3) if min_ligand is not None else None,
            "min_protein_distance": round(min_protein, 3) if min_protein is not None else None,
            "kept": kept, "reason": reason,
            "atom_serial": representative.get("serial"),
        })
    return {
        "selection_label": "binding_site_water",
        "note": "Single-structure geometry/quality heuristic; conservation requires an ensemble or literature evidence.",
        "parameters": {
            "ligand_distance": ligand_distance,
            "protein_distance": protein_distance,
            "min_occupancy": min_occupancy,
            "max_b_factor": max_b_factor,
            "include_residues": include_residues or [],
            "exclude_residues": exclude_residues or [],
        },
        "selected_count": len({(a["chain_id"], a["res_seq"]) for a in selected_atoms}),
        "waters": entries,
        "selected_atoms": selected_atoms,
    }


def parse_pdb_extract_components(
    pdb_file: Path,
    output_dir: Path,
    ligand_ids: list[str] | None = None,
    keep_conserved_waters: bool = True,
    water_distance_threshold: float = 3.5,
    keep_auxiliary_hetero: bool = False,
    preserve_hetero_resnames: list[str] | None = None,
    water_protein_distance: float = 3.5,
    water_min_occupancy: float = 0.5,
    water_max_b_factor: float = 60.0,
    water_include_residues: list[str] | None = None,
    water_exclude_residues: list[str] | None = None,
) -> dict[str, Any]:
    atoms = parse_pdb_atoms(pdb_file)
    modified_polymer = _parse_modres_records(pdb_file)
    modified_polymer_keys = set(modified_polymer)
    chain_info = get_chain_info(atoms, modified_polymer_keys)
    preserved_hetero = {
        res.upper()
        for res in (preserve_hetero_resnames if preserve_hetero_resnames is not None else DEFAULT_PRESERVED_HETERO)
    }

    all_hetatm = set(a["res_name"] for a in atoms if a["record"] == "HETATM")
    candidate_resnames = all_hetatm - WATER_CODES - COMMON_IONS
    ligand_instances = {
        key: value
        for key, value in _ligand_instances(atoms, candidate_resnames).items()
        if key not in modified_polymer_keys
    }
    ligand_candidates = {key[0] for key in ligand_instances}

    selected_ligand = None
    selected_ligand_key = None
    if ligand_ids:
        for lid in ligand_ids:
            matches = [key for key in ligand_instances if key[0] == lid]
            if matches:
                selected_ligand_key = max(matches, key=lambda key: len(ligand_instances[key]))
                selected_ligand = selected_ligand_key[0]
                break
    if selected_ligand_key is None and ligand_instances:
        selected_ligand_key = max(ligand_instances, key=lambda key: len(ligand_instances[key]))
        selected_ligand = selected_ligand_key[0]

    water_report = {
        "selection_label": "binding_site_water", "selected_count": 0,
        "waters": [], "selected_atoms": [],
    }
    if keep_conserved_waters and selected_ligand_key:
        water_report = select_binding_site_waters(
            atoms, selected_ligand_key,
            ligand_distance=water_distance_threshold,
            protein_distance=water_protein_distance,
            min_occupancy=water_min_occupancy,
            max_b_factor=water_max_b_factor,
            include_residues=water_include_residues,
            exclude_residues=water_exclude_residues,
        )
    conserved_waters = water_report.pop("selected_atoms", [])
    conserved_water_keys = {(w["chain_id"], w["res_seq"]) for w in conserved_waters}

    pdb_id = pdb_file.stem
    receptor_lines = []
    selected_ligand_serials = {
        a["serial"]
        for a in atoms
        if selected_ligand_key and _atom_instance_key(a) == selected_ligand_key
    }
    receptor_serials = set()
    removed_auxiliary: dict[tuple[str, str, int, str], int] = {}
    preserved_auxiliary: dict[tuple[str, str, int, str], int] = {}
    preserved_modified: dict[tuple[str, str, int, str], int] = {}
    for a in atoms:
        if a["record"] == "ATOM":
            receptor_lines.append(a["raw_line"])
            receptor_serials.add(a["serial"])
        elif a["record"] == "HETATM" and a["res_name"] in WATER_CODES:
            if (a["chain_id"], a["res_seq"]) in conserved_water_keys:
                receptor_lines.append(a["raw_line"])
                receptor_serials.add(a["serial"])
        elif a["record"] == "HETATM" and a["serial"] not in selected_ligand_serials:
            key = _atom_instance_key(a)
            is_modified_polymer = bool(key and key in modified_polymer_keys)
            keep_hetero = (
                is_modified_polymer
                or keep_auxiliary_hetero
                or a["res_name"].upper() in preserved_hetero
            )
            if keep_hetero:
                receptor_lines.append(a["raw_line"])
                receptor_serials.add(a["serial"])
                if key and is_modified_polymer:
                    preserved_modified[key] = preserved_modified.get(key, 0) + 1
                elif key:
                    preserved_auxiliary[key] = preserved_auxiliary.get(key, 0) + 1
            elif key:
                removed_auxiliary[key] = removed_auxiliary.get(key, 0) + 1

    receptor_conect_lines = _conect_lines_for_serials(pdb_file, receptor_serials)
    receptor_residue_keys = {
        (a["res_name"], a["chain_id"], a["res_seq"], a.get("i_code", ""))
        for a in atoms
        if a["serial"] in receptor_serials
    }
    receptor_annotation_lines = _modified_residue_annotation_lines(
        pdb_file,
        receptor_residue_keys=receptor_residue_keys,
        modified_polymer_keys=modified_polymer_keys,
    )

    receptor_file = output_dir / f"{pdb_id}_receptor.pdb"
    with open(receptor_file, "w") as f:
        f.writelines(receptor_annotation_lines)
        f.writelines(receptor_lines)
        if receptor_conect_lines:
            f.writelines(receptor_conect_lines)
        f.write("END\n")

    ligand_file = None
    if selected_ligand_key:
        lig_atoms = [a for a in atoms if _atom_instance_key(a) == selected_ligand_key]
        lig_lines = [a["raw_line"] for a in lig_atoms]
        if lig_lines:
            # Include CONECT records for proper RDKit bond topology perception
            lig_serials = {a["serial"] for a in lig_atoms}
            conect_lines = _conect_lines_for_serials(pdb_file, lig_serials)

            ligand_file = output_dir / f"{pdb_id}_ligand_{selected_ligand}.pdb"
            with open(ligand_file, "w") as f:
                f.writelines(lig_lines)
                if conect_lines:
                    f.writelines(conect_lines)
                f.write("END\n")

    return {
        "receptor_file": receptor_file,
        "ligand_file": ligand_file,
        "selected_ligand": selected_ligand,
        "selected_ligand_instance": {
            "res_name": selected_ligand_key[0],
            "chain_id": selected_ligand_key[1],
            "res_seq": selected_ligand_key[2],
            "i_code": selected_ligand_key[3],
        } if selected_ligand_key else None,
        "all_hetatm_residues": sorted(all_hetatm),
        "chain_info": chain_info,
        "conserved_waters_count": len(conserved_water_keys),
        "water_selection": water_report,
        "ligand_candidates": sorted(ligand_candidates),
        "keep_auxiliary_hetero": keep_auxiliary_hetero,
        "preserve_hetero_resnames": sorted(preserved_hetero),
        "removed_auxiliary_hetero": _summarize_hetero_instances(removed_auxiliary),
        "preserved_auxiliary_hetero": _summarize_hetero_instances(preserved_auxiliary),
        "preserved_modified_polymer_residues": (
            _summarize_modified_polymer_residues(
                modified_polymer, preserved_modified,
            )
        ),
    }


def _parse_residue_selector(value: str) -> tuple[str, int] | None:
    try:
        chain, res_seq = str(value).split(":", 1)
        return chain.strip(), int(res_seq)
    except (TypeError, ValueError):
        return None


def _atom_instance_key(atom: dict[str, Any]) -> tuple[str, str, int, str] | None:
    if atom.get("record") != "HETATM":
        return None
    return (atom["res_name"], atom["chain_id"], atom["res_seq"], atom.get("i_code", ""))


def _ligand_instances(
    atoms: list[dict[str, Any]], ligand_candidates: set[str]
) -> dict[tuple[str, str, int, str], list[dict[str, Any]]]:
    instances: dict[tuple[str, str, int, str], list[dict[str, Any]]] = {}
    for atom in atoms:
        key = _atom_instance_key(atom)
        if key and key[0] in ligand_candidates:
            instances.setdefault(key, []).append(atom)
    return instances


def _conect_lines_for_serials(pdb_file: Path, serials: set[int]) -> list[str]:
    if not serials:
        return []
    conect_lines = []
    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("CONECT"):
                continue
            try:
                line_serials = [
                    int(line[index:index + 5].strip())
                    for index in range(6, len(line), 5)
                    if line[index:index + 5].strip()
                ]
            except ValueError:
                continue
            if not line_serials or line_serials[0] not in serials:
                continue
            kept_targets = [serial for serial in line_serials[1:] if serial in serials]
            if kept_targets:
                conect_lines.append(
                    f"CONECT{line_serials[0]:5d}" + "".join(f"{serial:5d}" for serial in kept_targets) + "\n"
                )
    return conect_lines


def _summarize_hetero_instances(instances: dict[tuple[str, str, int, str], int]) -> list[dict[str, Any]]:
    return [
        {
            "res_name": key[0],
            "chain_id": key[1],
            "res_seq": key[2],
            "i_code": key[3],
            "atom_count": atom_count,
        }
        for key, atom_count in sorted(instances.items())
    ]


def _parse_modres_records(
    pdb_file: Path,
) -> dict[tuple[str, str, int, str], dict[str, str]]:
    """Return PDB MODRES entries keyed like parsed HETATM residues."""
    records: dict[tuple[str, str, int, str], dict[str, str]] = {}
    with open(pdb_file) as handle:
        for line in handle:
            if not line.startswith("MODRES"):
                continue
            try:
                key = (
                    line[12:15].strip(),
                    line[16].strip(),
                    int(line[18:22].strip()),
                    line[22].strip(),
                )
            except (IndexError, ValueError):
                continue
            if not key[0]:
                continue
            records[key] = {
                "standard_res_name": line[24:27].strip(),
                "description": line[29:70].strip(),
                "raw_line": line,
            }
    return records


def _parse_link_residue_keys(
    line: str,
) -> tuple[tuple[str, str, int, str], tuple[str, str, int, str]] | None:
    try:
        first = (
            line[17:20].strip(), line[21].strip(),
            int(line[22:26].strip()), line[26].strip(),
        )
        second = (
            line[47:50].strip(), line[51].strip(),
            int(line[52:56].strip()), line[56].strip(),
        )
    except (IndexError, ValueError):
        return None
    return first, second


def _modified_residue_annotation_lines(
    pdb_file: Path,
    *,
    receptor_residue_keys: set[tuple[str, str, int, str]],
    modified_polymer_keys: set[tuple[str, str, int, str]],
) -> list[str]:
    """Keep MODRES/LINK evidence for modified residues retained in receptor."""
    lines: list[str] = []
    with open(pdb_file) as handle:
        for line in handle:
            if line.startswith("MODRES"):
                parsed = _parse_modres_line_key(line)
                if parsed in receptor_residue_keys:
                    lines.append(line)
            elif line.startswith("LINK"):
                endpoints = _parse_link_residue_keys(line)
                if (
                    endpoints
                    and all(key in receptor_residue_keys for key in endpoints)
                    and any(key in modified_polymer_keys for key in endpoints)
                ):
                    lines.append(line)
    return lines


def _parse_modres_line_key(line: str) -> tuple[str, str, int, str] | None:
    try:
        return (
            line[12:15].strip(), line[16].strip(),
            int(line[18:22].strip()), line[22].strip(),
        )
    except (IndexError, ValueError):
        return None


def _summarize_modified_polymer_residues(
    records: dict[tuple[str, str, int, str], dict[str, str]],
    atom_counts: dict[tuple[str, str, int, str], int],
) -> list[dict[str, Any]]:
    return [
        {
            "res_name": key[0],
            "chain_id": key[1],
            "res_seq": key[2],
            "i_code": key[3],
            "atom_count": atom_count,
            "standard_res_name": records[key].get("standard_res_name", ""),
            "description": records[key].get("description", ""),
        }
        for key, atom_count in sorted(atom_counts.items())
    ]

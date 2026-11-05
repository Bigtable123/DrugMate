"""pyKVFinder utilities for independent cavity-detection evidence."""
from __future__ import annotations

import sys
from pathlib import Path

_SHARED_DIR = Path(__file__).resolve().parents[2] / "shared"
sys.path.insert(0, str(_SHARED_DIR))

from residue_utils import (  # noqa: E402
    deduplicate_residue_records,
    legacy_residue_label,
    make_residue_record,
)


def run_pykvfinder(
    protein_file: Path,
    output_dir: Path,
    probe_in: float = 1.4,
    probe_out: float = 4.0,
    removal_distance: float = 2.4,
    volume_cutoff: float = 5.0,
) -> dict:
    """Run pyKVFinder cavity detection on a PDB file.

    Compatible with pyKVFinder >= 0.9.x API.
    """
    import pyKVFinder
    import numpy as np

    output_dir.mkdir(parents=True, exist_ok=True)

    atomic = pyKVFinder.read_pdb(str(protein_file))

    vertices = pyKVFinder.get_vertices(atomic, probe_out=probe_out, step=0.6)

    ncavities, cavities = pyKVFinder.detect(
        atomic,
        vertices,
        step=0.6,
        probe_in=probe_in,
        probe_out=probe_out,
        removal_distance=removal_distance,
        volume_cutoff=volume_cutoff,
    )

    if ncavities == 0:
        return {
            "num_pockets": 0,
            "top_pocket": None,
            "all_pockets": [],
        }

    # spatial returns (surface, volume_dict, area_dict) in 0.9.x
    surface, volumes, _ = pyKVFinder.spatial(cavities, step=0.6)

    # constitutional returns dict[cavity_id, list[list[str]]]
    # each inner list is [resnum, chain, resname] in pyKVFinder 0.9.x
    residue_info = pyKVFinder.constitutional(cavities, atomic, vertices, step=0.6)
    pdb_residue_lookup = _pdb_residue_lookup(protein_file)

    pockets = []
    for cav_id in sorted(volumes.keys()):
        pocket_records = []
        if cav_id in residue_info:
            for res_entry in residue_info[cav_id]:
                if isinstance(res_entry, (list, tuple)) and len(res_entry) >= 3:
                    resnum = res_entry[0]
                    chain = str(res_entry[1]).strip()
                    resname = res_entry[2]
                    matches = pdb_residue_lookup.get(
                        _lookup_key(chain, resnum, resname), []
                    )
                    if len(matches) == 1:
                        pocket_records.append(matches[0])
                    elif len(matches) > 1:
                        # pyKVFinder 0.9.x drops PDB insertion codes internally.
                        # Preserve every possible source residue but mark the
                        # identity ambiguous so consensus code will not count it
                        # as residue-IoU evidence.
                        candidate_icodes = sorted({item.get("icode", "") for item in matches})
                        for item in matches:
                            pocket_records.append({
                                **item,
                                "identity_ambiguous": True,
                                "candidate_icodes": candidate_icodes,
                                "ambiguity_reason": "pykvfinder_drops_insertion_code",
                            })
                    else:
                        pocket_records.append(
                            make_residue_record(resname, resnum, chain, "")
                        )

        pocket_records = deduplicate_residue_records(pocket_records)
        pocket_residues = [legacy_residue_label(record) for record in pocket_records]

        center = _compute_cavity_center(cavities, cav_id, vertices, step=0.6)

        pockets.append({
            "rank": len(pockets) + 1,
            "cavity_id": cav_id,
            # Native volume determines rank; display rounding can create ties
            # and change which cavity is retained as Top-1.
            "volume": float(volumes.get(cav_id, 0.0)),
            "residues": pocket_residues,
            "residue_records": pocket_records,
            "center": center,
        })

    pockets.sort(key=lambda p: p["volume"], reverse=True)
    for i, p in enumerate(pockets):
        p["rank"] = i + 1

    result_pdb = output_dir / f"{protein_file.stem}_cavities.pdb"
    # The point5 adapter consumes this cloud; never reuse an earlier run after
    # a failed export into the same output directory.
    result_pdb.unlink(missing_ok=True)
    pyKVFinder.export(str(result_pdb), cavities, surface, vertices, step=0.6)

    for pocket in pockets:
        pocket["support_file"] = str(result_pdb.resolve())
        pocket["support_resname"] = pocket["cavity_id"]

    return {
        "num_pockets": len(pockets),
        "top_pocket": pockets[0] if pockets else None,
        "all_pockets": pockets[:10],
    }


def _compute_cavity_center(
    cavities, cavity_id: str, vertices, step: float = 0.6
) -> list[float] | None:
    """Compute geometric center of a cavity from the grid."""
    import numpy as np

    try:
        # cavities is int32 array; cavity_id like "KAF" maps to integer index
        # KAA=2, KAB=3, ..., KAZ=27, KBA=28, KBB=29, ...
        # General formula: (second_char_offset * 26) + third_char_offset + 2
        int_id = (ord(cavity_id[1]) - ord('A')) * 26 + (ord(cavity_id[2]) - ord('A')) + 2
        indices = np.argwhere(cavities == int_id)
        if len(indices) == 0:
            return None

        origin = np.array([vertices[0][0], vertices[0][1], vertices[0][2]])
        coords = origin + indices * step
        center = coords.mean(axis=0)
        return [float(value) for value in center]
    except Exception:
        return None


def _pdb_residue_lookup(protein_file: Path) -> dict[tuple[str, str, str], list[dict]]:
    """Index original PDB residues so unique insertion codes can be restored."""
    lookup: dict[tuple[str, str, str], list[dict]] = {}
    seen = set()
    with open(protein_file) as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                record = make_residue_record(
                    line[17:20].strip(), line[22:26].strip(),
                    line[21:22].strip(), line[26:27].strip(),
                )
            except (ValueError, IndexError):
                continue
            exact = (
                record["chain"], record["resnum"], record["icode"], record["resname"]
            )
            if exact in seen:
                continue
            seen.add(exact)
            lookup.setdefault(
                _lookup_key(record["chain"], record["resnum"], record["resname"]), []
            ).append(record)
    return lookup


def _lookup_key(chain, resnum, resname) -> tuple[str, str, str]:
    try:
        normalized_number = str(int(str(resnum).strip()))
    except (TypeError, ValueError):
        normalized_number = str(resnum).strip()
    return str(chain or "").strip(), normalized_number, str(resname or "").strip().upper()

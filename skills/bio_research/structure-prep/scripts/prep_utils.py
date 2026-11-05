"""Protein preparation utilities: PDBFixer, PROPKA, residue validation."""

from __future__ import annotations

import importlib
import logging
import math
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# The metals ``validate_metals`` recognizes when the caller names none.  This
# mirrors ``shared/pipeline_registry.py``'s ``metal_validation.allowed_metals``
# default, which only reaches this module through ``structure-prep/run.py``; a
# direct caller used to get an empty set and therefore never found a metal, so
# the report said "no metals" for a structure full of them.  Passing an explicit
# empty list still means "recognize nothing".
DEFAULT_ALLOWED_METALS: tuple[str, ...] = (
    "MG", "ZN", "CA", "FE", "MN", "NI", "CU", "CO", "NA", "K", "CD",
)


def check_dependencies() -> dict[str, bool]:
    deps = {}
    for mod in ["pdbfixer", "openmm", "propka"]:
        try:
            importlib.import_module(mod)
            deps[mod] = True
        except ImportError:
            deps[mod] = False
    return deps


def fix_protein_with_pdbfixer(
    input_pdb: Path, output_pdb: Path,
    ph: float | None = None, add_hydrogens: bool = False,
    add_missing_residues: bool = False, add_missing_atoms: bool = True,
    max_missing_residues: int = 3, seed: int = 0,
) -> dict[str, Any]:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    logger.info(f"Running PDBFixer on {input_pdb}")
    fixer = PDBFixer(filename=str(input_pdb))

    fixer.findMissingResidues()
    missing_residues = dict(fixer.missingResidues)
    logger.info(f"  Missing residues: {len(missing_residues)} gaps")

    if not add_missing_residues:
        fixer.missingResidues = {}
    else:
        chains = list(fixer.topology.chains())
        unsupported = {
            key: names for key, names in missing_residues.items()
            if len(names) > max_missing_residues
            or key[1] in (0, len(list(chains[key[0]].residues())))
        }
        if unsupported:
            raise ValueError(
                "Only explicitly justified short internal gaps are supported; "
                f"terminal/long gaps require another reviewed structure: {unsupported}"
            )

    fixer.findMissingAtoms()
    missing_atoms = {str(k): len(v) for k, v in fixer.missingAtoms.items()} if fixer.missingAtoms else {}
    missing_terminals = {str(k): len(v) for k, v in fixer.missingTerminals.items()} if fixer.missingTerminals else {}
    logger.info(f"  Missing atoms in {len(missing_atoms)} residues")

    if add_missing_atoms:
        fixer.addMissingAtoms(seed=seed)

    if add_hydrogens:
        raise ValueError("Structure preparation is geometry-only; add hydrogens in docking_prep")

    with open(output_pdb, "w") as f:
        # Preserve source PDB residue identifiers so interaction reports remain
        # mappable to PDB/UniProt numbering (including valid residue 0 and
        # insertion codes).  OpenMM otherwise renumbers every chain from 1.
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    logger.info(f"  Prepared protein saved to {output_pdb}")
    return {
        "missing_residues_count": len(missing_residues),
        "missing_residues": {str(k): list(v) for k, v in missing_residues.items()},
        "missing_atoms_residues": len(missing_atoms),
        "missing_terminals": len(missing_terminals),
        "hydrogens_added": add_hydrogens,
        "ph": ph,
        "seed": seed,
        "max_missing_residues": max_missing_residues,
    }


def protonate_with_propka(input_pdb: Path, ph: float = 7.4) -> dict[str, Any]:
    import propka.run

    logger.info(f"Running PROPKA at pH {ph} on {input_pdb}")
    mol = propka.run.single(str(input_pdb), write_pka=False)

    pka_results = []
    for conformation in mol.conformations:
        conf = mol.conformations[conformation]
        for group in conf.groups:
            # PROPKA exposes backbone/amide bookkeeping groups alongside real
            # titratable groups. Counting every group with a numeric pKa inflates
            # the report by hundreds of residues; only chemically titratable
            # groups belong in this diagnostic.
            if (
                getattr(group, "titratable", False)
                and hasattr(group, "pka_value")
                and group.pka_value is not None
            ):
                pka_results.append({
                    "residue": f"{group.residue_type}{group.atom.res_num}",
                    "chain": group.atom.chain_id,
                    "pka": round(group.pka_value, 2),
                    "protonated_at_target_ph": group.pka_value > ph,
                })
        break

    return {"ph": ph, "titratable_residues": pka_results, "total_titratable": len(pka_results)}


def validate_residues(pdb_file: Path, required_residues: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    atoms = _parse_atoms(pdb_file)

    residue_map: dict[tuple[str, int], str] = {}
    residue_atoms: dict[tuple[str, int], int] = {}
    for a in atoms:
        if a["record"] == "ATOM":
            key = (a["chain_id"], a["res_seq"])
            residue_map[key] = a["res_name"]
            residue_atoms[key] = residue_atoms.get(key, 0) + 1

    total_residues = len(residue_map)

    EXPECTED_MIN = {"GLY": 4, "ALA": 5}
    incomplete = []
    for key, count in residue_atoms.items():
        res_name = residue_map[key]
        if count < EXPECTED_MIN.get(res_name, 4):
            incomplete.append({"chain": key[0], "res_seq": key[1], "res_name": res_name, "atom_count": count})

    required_check = []
    if required_residues:
        for req in required_residues:
            chain = req.get("chain", "A")
            seq = req["seq"]
            expected_name = req.get("name", "").upper()
            key = (chain, seq)
            found = key in residue_map
            actual_name = residue_map.get(key, None)
            matches = actual_name == expected_name if expected_name and actual_name else None
            required_check.append({"expected": f"{expected_name}{seq} (chain {chain})", "found": found, "actual_name": actual_name, "name_matches": matches})

    all_present = all(r["found"] for r in required_check) if required_check else True
    no_mutations = all(r["name_matches"] is not False for r in required_check) if required_check else True

    return {
        "total_residues": total_residues,
        "incomplete_residues": incomplete,
        "required_residues_check": required_check,
        "all_required_present": all_present,
        "no_mutations_detected": no_mutations,
        "validation_passed": all_present and no_mutations,
    }


def _parse_atoms(pdb_file: Path) -> list[dict]:
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            atoms.append({
                "record": record,
                "name": line[12:16].strip(),
                "res_name": line[17:20].strip(),
                "chain_id": line[21].strip(),
                "res_seq": int(line[22:26].strip()),
                "alt_loc": line[16].strip() if len(line) > 16 else "",
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": (line[76:78].strip() if len(line) > 76 else ""),
            })
    return atoms


def validate_metals(
    pdb_file: Path,
    *,
    coordination_cutoff: float = 3.0,
    allowed_metals: list[str] | None = None,
) -> dict[str, Any]:
    """Report local metal coordination geometry without assigning charges."""
    if allowed_metals is None:
        allowed_metals = list(DEFAULT_ALLOWED_METALS)
    allowed = {value.upper() for value in allowed_metals}
    atoms = _parse_atoms(pdb_file)
    metals = [
        atom for atom in atoms
        if atom["record"] == "HETATM"
        and (atom["res_name"].upper() in allowed or atom["element"].upper() in allowed)
    ]
    donors = [
        atom for atom in atoms
        if atom not in metals and atom["element"].upper() in {"N", "O", "S"}
    ]
    sites = []
    for metal in metals:
        neighbours = []
        for donor in donors:
            distance = math.sqrt(sum(
                (metal[axis] - donor[axis]) ** 2 for axis in ("x", "y", "z")
            ))
            if distance <= coordination_cutoff:
                neighbours.append({
                    "atom": donor["name"], "res_name": donor["res_name"],
                    "chain_id": donor["chain_id"], "res_seq": donor["res_seq"],
                    "element": donor["element"], "distance": round(distance, 3),
                })
        neighbours.sort(key=lambda entry: entry["distance"])
        sites.append({
            "metal": metal["res_name"], "chain_id": metal["chain_id"],
            "res_seq": metal["res_seq"], "atom": metal["name"],
            "coordination_count": len(neighbours), "neighbours": neighbours,
        })
    return {
        "metal_count": len(metals), "coordination_cutoff": coordination_cutoff,
        "allowed_metals": sorted(allowed), "sites": sites,
        "warning": (
            "Geometry only: formal charge and force-field parameters were not assigned."
            if metals else None
        ),
    }


def audit_structure(
    pdb_file: Path,
    *,
    fail_on_altloc: bool = False,
    require_hydrogens: bool = True,
    expected_outputs: dict[str, str | None] | None = None,
) -> dict[str, Any]:
    """Run inexpensive final structure/file sanity checks."""
    atoms = _parse_atoms(pdb_file)
    hydrogens = [
        atom for atom in atoms
        if atom["element"].upper() in {"H", "D"}
        or atom["name"].lstrip("0123456789")[:1].upper() in {"H", "D"}
    ]
    altloc_atoms = [atom for atom in atoms if atom.get("alt_loc")]
    files = {}
    for name, value in (expected_outputs or {}).items():
        path = Path(value) if value else None
        files[name] = {
            "path": str(path) if path else None,
            "exists": bool(path and path.exists() and path.stat().st_size > 0),
        }
    issues = []
    notes: list[str] = []
    if not atoms:
        issues.append("prepared receptor contains no ATOM/HETATM records")
    # Only a run that asked for hydrogens can be missing them.  With
    # ``add_hydrogens=False`` the absence is the caller's configuration choice,
    # and the executable policy keeps a choice out of the fault column -- reporting it as a
    # QC issue failed the whole preparation for doing exactly what it was told.
    if require_hydrogens and not hydrogens:
        issues.append("prepared receptor contains no explicit hydrogens")
    elif not hydrogens:
        notes.append(
            "prepared receptor contains no explicit hydrogens; "
            "add_hydrogens was not requested"
        )
    if fail_on_altloc and altloc_atoms:
        issues.append(f"prepared receptor retains {len(altloc_atoms)} alternate-location atoms")
    missing_expected = [name for name, info in files.items() if info["path"] and not info["exists"]]
    if missing_expected:
        issues.append("missing/empty expected outputs: " + ", ".join(missing_expected))
    return {
        "passed": not issues,
        "atom_count": len(atoms),
        "hydrogen_count": len(hydrogens),
        "hydrogens_requested": bool(require_hydrogens),
        "alternate_location_atom_count": len(altloc_atoms),
        "issues": issues,
        "notes": notes,
        "files": files,
    }

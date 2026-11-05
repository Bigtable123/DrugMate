"""DOCK 6 classic-docking adapter for the canonical Step 4 contract.

The adapter deliberately does not synthesize a receptor grid from a prepared
PDB.  A trustworthy DOCK 6 route requires receptor-site spheres and the paired
``.nrg``/``.bmp`` grid generated from a charged, atom-typed receptor.  Those
artifacts are explicit pipeline inputs; ligand SDF records are converted with
the DOCK installation's own Antechamber binary and are audited one-for-one.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import tempfile
from pathlib import Path
from typing import Any

from deployment_config import deployment_path
from docking_common import (
    annotate_canonical_pose,
    compound_name,
    copy_identity_properties,
)
from tool_runtime import ToolRuntimeError, run_tool


DEFAULT_DOCK_ROOT = deployment_path("${SBDD_DOCK6_SOURCE_DIR}")
DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "dock6",
    "command": str(DEFAULT_DOCK_ROOT / "bin/dock6"),
}
_HEADER = re.compile(r"^##########\s+(?P<name>[^:]+):\s*(?P<value>.*?)\s*$")
_RECORD_START = re.compile(r"(?m)(?=^##########\s+Name:)")


def run_dock6_docking(
    *,
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    receptor_site_file: str | Path | None,
    grid_nrg_file: str | Path | None,
    grid_bmp_file: str | Path | None,
    runtime: dict | None = None,
    dock_root: str | Path = DEFAULT_DOCK_ROOT,
    antechamber_path: str | Path | None = None,
    charge_method: str = "gas",
    num_modes: int = 5,
    max_orientations: int = 500,
    min_anchor_size: int = 5,
    pruning_max_orients: int = 100,
    score_threshold: float = 1000.0,
    seed: int = 42,
    timeout_seconds: int = 86400,
    ligand_prep_timeout_seconds: int = 600,
) -> dict:
    """Run classic DOCK 6 using injected receptor-site/grid artifacts."""
    dock_root = Path(dock_root).resolve()
    receptor_site = Path(receptor_site_file).resolve() if receptor_site_file else None
    grid_nrg = Path(grid_nrg_file).resolve() if grid_nrg_file else None
    grid_bmp = Path(grid_bmp_file).resolve() if grid_bmp_file else None
    if receptor_site is None or grid_nrg is None or grid_bmp is None:
        return {
            "success": False,
            "error": (
                "DOCK 6 requires dock6_receptor_site_file, "
                "dock6_grid_nrg_file, and dock6_grid_bmp_file; a prepared PDB "
                "alone is not a charged/atom-typed DOCK receptor grid"
            ),
        }
    if grid_nrg.suffix.lower() != ".nrg" or grid_bmp.suffix.lower() != ".bmp":
        return {
            "success": False,
            "error": "DOCK 6 grid inputs must use .nrg and .bmp suffixes",
        }
    if grid_nrg.with_suffix("") != grid_bmp.with_suffix(""):
        return {
            "success": False,
            "error": "DOCK 6 .nrg and .bmp files must share one grid prefix",
        }
    grid_prefix = grid_nrg.with_suffix("")
    antechamber = Path(
        antechamber_path or dock_root / "bin/antechamber"
    ).resolve()
    parameters_dir = dock_root / "parameters"
    required = {
        "DOCK 6 receptor-site spheres": receptor_site,
        "DOCK 6 energy grid": grid_nrg,
        "DOCK 6 bump grid": grid_bmp,
        "DOCK 6 Antechamber": antechamber,
        "DOCK 6 VDW definition": parameters_dir / "vdw_AMBER_parm99.defn",
        "DOCK 6 flex definition": parameters_dir / "flex.defn",
        "DOCK 6 flex-drive table": parameters_dir / "flex_drive.tbl",
    }
    missing = [
        label
        for label, path in required.items()
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        return {
            "success": False,
            "error": "missing or empty required DOCK 6 resource(s): " + ", ".join(missing),
        }
    if charge_method not in {"gas", "bcc"}:
        return {
            "success": False,
            "error": "DOCK 6 ligand charge_method must be gas or bcc",
        }
    try:
        num_modes = int(num_modes)
        max_orientations = int(max_orientations)
        min_anchor_size = int(min_anchor_size)
        pruning_max_orients = int(pruning_max_orients)
        seed = int(seed)
        timeout_seconds = int(timeout_seconds)
        ligand_prep_timeout_seconds = int(ligand_prep_timeout_seconds)
        score_threshold = float(score_threshold)
        if min(
            num_modes,
            max_orientations,
            min_anchor_size,
            pruning_max_orients,
            timeout_seconds,
            ligand_prep_timeout_seconds,
        ) < 1:
            raise ValueError("integer search/timeout parameters must be >= 1")
        if not math.isfinite(score_threshold):
            raise ValueError("score_threshold must be finite")
    except (TypeError, ValueError) as exc:
        return {"success": False, "error": f"invalid DOCK 6 parameters: {exc}"}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "dock6_native"
    native_root.mkdir(parents=True, exist_ok=True)
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    log_file = output_dir / "dock6.log"
    docked_sdf = output_dir / "docked_poses.sdf"
    docked_sdf.unlink(missing_ok=True)
    prepared_ligands = native_output_dir / "ligands.mol2"
    resolved_runtime = dict(runtime or DEFAULT_RUNTIME)
    dock_environment = {
        "DOCK_HOME": str(dock_root),
        "AMBERHOME": str(dock_root),
    }

    try:
        molecules, preparation = _prepare_ligands(
            ligand_file=Path(ligand_file),
            output_file=prepared_ligands,
            native_output_dir=native_output_dir,
            runtime=resolved_runtime,
            antechamber_path=antechamber,
            charge_method=charge_method,
            timeout_seconds=ligand_prep_timeout_seconds,
            environment=dock_environment,
        )
    except Exception as exc:
        return {
            "success": False,
            "error": f"could not prepare DOCK 6 ligands: {exc}",
            "native_output_dir": str(native_output_dir),
        }

    output_prefix = native_output_dir / "dock6"
    input_file = native_output_dir / "dock6.in"
    input_file.write_text(
        build_dock6_input(
            ligand_atom_file=prepared_ligands,
            receptor_site_file=receptor_site,
            grid_prefix=grid_prefix,
            output_prefix=output_prefix,
            parameters_dir=parameters_dir,
            num_modes=num_modes,
            max_orientations=max_orientations,
            min_anchor_size=min_anchor_size,
            pruning_max_orients=pruning_max_orients,
            score_threshold=score_threshold,
            seed=seed,
        ),
        encoding="utf-8",
    )
    try:
        process = run_tool(
            resolved_runtime,
            ["-i", str(input_file)],
            cwd=native_output_dir,
            timeout=timeout_seconds,
            check=False,
            env=dock_environment,
        )
    except (OSError, ToolRuntimeError) as exc:
        return {
            "success": False,
            "error": f"could not execute DOCK 6 runtime: {exc}",
            "native_output_dir": str(native_output_dir),
        }
    _write_log(log_file, process, preparation)
    command_record = _command_record(process)
    if not process.get("success"):
        return {
            "success": False,
            "error": f"DOCK 6 exited with code {process.get('returncode')}",
            "log_file": str(log_file),
            "native_output_dir": str(native_output_dir),
            "command": command_record,
        }
    scored_mol2 = Path(f"{output_prefix}_scored.mol2")
    if not scored_mol2.is_file() or scored_mol2.stat().st_size == 0:
        return {
            "success": False,
            "error": "DOCK 6 exited successfully but produced no scored Mol2 output",
            "log_file": str(log_file),
            "native_output_dir": str(native_output_dir),
            "command": command_record,
        }

    try:
        merge = _merge_dock6_outputs(
            scored_mol2=scored_mol2,
            input_molecules=molecules,
            max_poses_per_molecule=num_modes,
            output_sdf=docked_sdf,
        )
    except Exception as exc:
        return {
            "success": False,
            "error": f"could not validate DOCK 6 outputs: {exc}",
            "log_file": str(log_file),
            "native_output_dir": str(native_output_dir),
            "command": command_record,
        }

    protocol_assets = {
        "receptor_site_spheres": _file_identity(receptor_site),
        "grid_energy": _file_identity(grid_nrg),
        "grid_bump": _file_identity(grid_bmp),
        "vdw_definition": _file_identity(parameters_dir / "vdw_AMBER_parm99.defn"),
        "flex_definition": _file_identity(parameters_dir / "flex.defn"),
        "flex_drive_table": _file_identity(parameters_dir / "flex_drive.tbl"),
    }
    return {
        "success": True,
        "engine": "dock6",
        "docked_sdf": str(docked_sdf),
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "num_poses": len(merge["scores"]),
        "num_input_molecules": len(molecules),
        "scores": merge["scores"],
        "score_metadata": score_metadata(),
        "redocking_validation": {
            "performed": False,
            "passed": None,
            "applicable": False,
            "not_applicable": True,
            "reason": "DOCK 6 route has no registered engine-neutral redocking protocol",
            "gate": "not_applicable",
        },
        "quality_gate_passed": None,
        "search_space": {
            "mode": "dock6_precomputed_site_grid",
            "declared_pocket_context": docking_box,
            "reference_ligand_supplied": reference_ligand_file is not None,
        },
        "software": {
            "engine": "DOCK 6",
            "engine_version": _dock6_version(
                resolved_runtime, dock_environment,
            ),
            "runtime_environment": resolved_runtime.get("environment", "dock6"),
            "frontend": "dock6",
        },
        "protocol_assets": protocol_assets,
        "command": command_record,
        "ligand_preparation": preparation,
        "merge_diagnostics": merge["diagnostics"],
        # Carried through so the failure registry can name the poses and inputs
        # DOCK 6 could not score; empty for a complete batch.
        "operational_errors": merge.get("operational_errors", []),
    }


def build_dock6_input(
    *,
    ligand_atom_file: Path,
    receptor_site_file: Path,
    grid_prefix: Path,
    output_prefix: Path,
    parameters_dir: Path,
    num_modes: int,
    max_orientations: int,
    min_anchor_size: int,
    pruning_max_orients: int,
    score_threshold: float,
    seed: int,
) -> str:
    """Build an explicit DOCK 6.13 flexible/grid-scoring input file."""
    values = [
        ("conformer_search_type", "flex"),
        ("write_fragment_libraries", "no"),
        ("user_specified_anchor", "no"),
        ("limit_max_anchors", "no"),
        ("min_anchor_size", min_anchor_size),
        ("pruning_use_clustering", "yes"),
        ("pruning_max_orients", pruning_max_orients),
        ("pruning_clustering_cutoff", 100),
        ("pruning_orient_score_cutoff", 1000.0),
        ("pruning_conformer_score_cutoff", 25.0),
        ("pruning_conformer_score_scaling_factor", 1.0),
        ("use_clash_overlap", "no"),
        ("write_growth_tree", "no"),
        ("use_internal_energy", "yes"),
        ("internal_energy_rep_exp", 12),
        ("internal_energy_cutoff", 100.0),
        ("ligand_atom_file", ligand_atom_file),
        ("limit_max_ligands", "no"),
        ("skip_molecule", "no"),
        ("read_mol_solvation", "no"),
        ("calculate_rmsd", "no"),
        ("use_database_filter", "no"),
        ("orient_ligand", "yes"),
        ("automated_matching", "yes"),
        ("receptor_site_file", receptor_site_file),
        ("max_orientations", max_orientations),
        ("critical_points", "no"),
        ("chemical_matching", "no"),
        ("use_ligand_spheres", "no"),
        ("bump_filter", "no"),
        ("score_molecules", "yes"),
        ("contact_score_primary", "no"),
        ("grid_score_primary", "yes"),
        ("grid_score_rep_rad_scale", 1),
        ("grid_score_vdw_scale", 1),
        ("grid_score_es_scale", 1),
        ("grid_lig_efficiency", "no"),
        ("grid_score_grid_prefix", grid_prefix),
        ("minimize_ligand", "yes"),
        ("minimize_anchor", "yes"),
        ("minimize_flexible_growth", "yes"),
        ("use_advanced_simplex_parameters", "no"),
        ("minimize_flexible_growth_ramp", "no"),
        ("simplex_max_cycles", 1),
        ("simplex_score_converge", 0.1),
        ("simplex_cycle_converge", 1.0),
        ("simplex_trans_step", 1.0),
        ("simplex_rot_step", 0.1),
        ("simplex_tors_step", 10.0),
        ("simplex_anchor_max_iterations", 500),
        ("simplex_grow_max_iterations", 500),
        ("simplex_grow_tors_premin_iterations", 0),
        ("simplex_final_min", "no"),
        ("simplex_random_seed", seed),
        ("simplex_restraint_min", "no"),
        ("atom_model", "all"),
        ("vdw_defn_file", parameters_dir / "vdw_AMBER_parm99.defn"),
        ("flex_defn_file", parameters_dir / "flex.defn"),
        ("flex_drive_file", parameters_dir / "flex_drive.tbl"),
        ("ligand_outfile_prefix", output_prefix),
        ("write_mol_solvation", "no"),
        ("write_orientations", "no"),
        ("num_final_scored_poses", num_modes),
        ("num_preclustered_conformers", num_modes),
        ("score_threshold", score_threshold),
        ("rank_ligands", "no"),
    ]
    return "".join(f"{name:<64}{value}\n" for name, value in values)


def _prepare_ligands(
    *,
    ligand_file: Path,
    output_file: Path,
    native_output_dir: Path,
    runtime: dict,
    antechamber_path: Path,
    charge_method: str,
    timeout_seconds: int,
    environment: dict[str, str],
) -> tuple[list[Any], dict]:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(
        str(ligand_file), removeHs=False, strictParsing=False,
    )
    molecules: list[Any] = []
    blocks: list[str] = []
    command_records: list[dict] = []
    invalid: list[int] = []
    preparation_root = native_output_dir / "ligand_preparation"
    preparation_root.mkdir()
    for index, molecule in enumerate(supplier):
        if molecule is None:
            invalid.append(index)
            continue
        if molecule.GetNumConformers() != 1:
            raise ValueError(f"ligand record {index} must contain one conformer")
        token = f"drugmate_src_{index:06d}"
        record_dir = preparation_root / token
        record_dir.mkdir()
        mdl_file = record_dir / f"{token}.mol"
        mol2_file = record_dir / f"{token}.mol2"
        mdl_file.write_text(Chem.MolToMolBlock(molecule), encoding="utf-8")
        net_charge = int(Chem.GetFormalCharge(molecule))
        result = run_tool(
            runtime,
            [
                "-i", str(mdl_file),
                "-fi", "mdl",
                "-o", str(mol2_file),
                "-fo", "mol2",
                "-c", charge_method,
                "-nc", str(net_charge),
                "-at", "sybyl",
                "-pf", "y",
                "-s", "0",
            ],
            command=str(antechamber_path),
            cwd=record_dir,
            timeout=timeout_seconds,
            check=False,
            env=environment,
        )
        command_records.append(_command_record(result))
        if not result.get("success"):
            detail = (result.get("stderr") or result.get("stdout") or "").strip()
            raise ValueError(
                f"Antechamber failed for ligand {index} with code "
                f"{result.get('returncode')}: {detail[-500:]}"
            )
        if not mol2_file.is_file() or mol2_file.stat().st_size == 0:
            raise ValueError(f"Antechamber produced no Mol2 for ligand {index}")
        block = _replace_mol2_name(
            mol2_file.read_text(encoding="utf-8", errors="replace"), token,
        )
        if block.count("@<TRIPOS>MOLECULE") != 1:
            raise ValueError(f"prepared ligand {index} is not one Mol2 record")
        blocks.append(block.rstrip() + "\n")
        molecules.append(molecule)
    if invalid:
        raise ValueError(f"input ligand SDF has unparsable record(s): {invalid}")
    if not molecules:
        raise ValueError("input ligand SDF contains no parseable molecules")
    output_file.write_text("".join(blocks), encoding="utf-8")
    return molecules, {
        "method": "DOCK 6 Antechamber",
        "charge_method": charge_method,
        "atom_type": "sybyl",
        "num_input_molecules": len(molecules),
        "num_prepared_molecules": len(blocks),
        "commands": command_records,
    }


def _replace_mol2_name(text: str, name: str) -> str:
    lines = text.splitlines()
    try:
        marker = lines.index("@<TRIPOS>MOLECULE")
    except ValueError as exc:
        raise ValueError("Antechamber output has no @<TRIPOS>MOLECULE record") from exc
    if marker + 1 >= len(lines):
        raise ValueError("Antechamber Mol2 record has no molecule name")
    lines[marker + 1] = name
    return "\n".join(lines) + "\n"


def _merge_dock6_outputs(
    *,
    scored_mol2: Path,
    input_molecules: list[Any],
    max_poses_per_molecule: int,
    output_sdf: Path,
) -> dict:
    from rdkit import Chem

    text = scored_mol2.read_text(encoding="utf-8", errors="replace")
    chunks = [
        chunk for chunk in _RECORD_START.split(text)
        if "@<TRIPOS>MOLECULE" in chunk
    ]
    if not chunks:
        raise ValueError("scored Mol2 contains no pose records")
    partial_sdf = output_sdf.with_suffix(".partial.sdf")
    partial_sdf.unlink(missing_ok=True)
    writer = Chem.SDWriter(str(partial_sdf))
    scores: list[dict] = []
    counts = {index: 0 for index in range(len(input_molecules))}
    invalid_records: list[dict] = []
    try:
        for native_index, chunk in enumerate(chunks):
            headers = _parse_headers(chunk)
            token = str(headers.get("Name") or "").strip()
            match = re.fullmatch(r"drugmate_src_(\d{6})", token)
            if match is None:
                raise ValueError(
                    f"DOCK 6 pose {native_index} has unknown source name {token!r}"
                )
            source_index = int(match.group(1))
            if source_index not in counts:
                raise ValueError(
                    f"DOCK 6 pose references absent input molecule {source_index}"
                )
            grid_score = _finite_header(headers, "Grid_Score", required=True)
            grid_vdw = _finite_header(headers, "Grid_vdw_energy", required=False)
            grid_es = _finite_header(headers, "Grid_es_energy", required=False)
            mol2_block = chunk[chunk.index("@<TRIPOS>MOLECULE"):]
            molecule = Chem.MolFromMol2Block(
                mol2_block,
                sanitize=False,
                removeHs=False,
                cleanupSubstructures=False,
            )
            if molecule is None or molecule.GetNumConformers() != 1:
                invalid_records.append({
                    "native_pose_index": native_index,
                    "source_molecule_index": source_index,
                    "reason": "RDKit could not parse exactly one 3D Mol2 pose",
                })
                continue
            try:
                Chem.SanitizeMol(molecule)
            except Exception as exc:
                invalid_records.append({
                    "native_pose_index": native_index,
                    "source_molecule_index": source_index,
                    "reason": f"RDKit sanitization failed: {exc}",
                })
                continue
            source_molecule = input_molecules[source_index]
            copy_identity_properties(source_molecule, molecule)
            compound_id = compound_name(source_molecule, source_index)
            engine_pose_index = counts[source_index]
            engine_pose_rank = engine_pose_index + 1
            counts[source_index] += 1
            if counts[source_index] > max_poses_per_molecule:
                raise ValueError(
                    f"DOCK 6 emitted more than {max_poses_per_molecule} poses "
                    f"for input molecule {source_index}"
                )
            pose_index = len(scores)
            molecule.SetDoubleProp("dock6_grid_score", grid_score)
            molecule.SetProp("score_direction", "lower")
            if grid_vdw is not None:
                molecule.SetDoubleProp("dock6_grid_vdw_energy", grid_vdw)
            if grid_es is not None:
                molecule.SetDoubleProp("dock6_grid_es_energy", grid_es)
            identity = annotate_canonical_pose(
                molecule,
                compound_id=compound_id,
                pose_index=pose_index,
                source_molecule_index=source_index,
                engine="dock6",
                engine_pose_index=engine_pose_index,
                engine_pose_rank=engine_pose_rank,
            )
            writer.write(molecule)
            score = {
                "pose_index": pose_index,
                "source_molecule_index": source_index,
                "engine_pose_index": engine_pose_index,
                "engine_pose_rank": engine_pose_rank,
                "compound_id": compound_id,
                "mol_name": compound_id,
                "engine": "dock6",
                "dock6_grid_score": grid_score,
                "score_reliability": "native_engine_only",
                "native_score_fields": headers,
                **identity,
            }
            if grid_vdw is not None:
                score["dock6_grid_vdw_energy"] = grid_vdw
            if grid_es is not None:
                score["dock6_grid_es_energy"] = grid_es
            scores.append(score)
    finally:
        writer.close()
    if not scores:
        raise ValueError("DOCK 6 produced no finite scored poses")
    # An invalid pose record, or an input molecule DOCK 6 placed no finite pose
    # for, is that pose's / that molecule's missing value -- not the batch's
    # failure.  Raising here discarded every pose DOCK 6 did produce, so one
    # unparseable record ended docking for the whole cohort, while the executable policy
    # scopes a single tool fault to the single candidate and the executable policy drops a
    # missing value from numerator and denominator instead.  Both gaps stay named
    # in ``diagnostics`` and are declared in ``operational_errors`` for the
    # failure registry (reflective-coordinator/SKILL.md:260-266).
    incomplete = [index for index, count in counts.items() if count == 0]
    partial_sdf.replace(output_sdf)
    merged = {
        "success": True,
        "scores": scores,
        "diagnostics": {
            "native_pose_records": len(chunks),
            "published_pose_records": len(scores),
            "poses_per_input_molecule": counts,
            "invalid_records": invalid_records,
            "unscored_source_molecule_indices": incomplete,
            "complete": not invalid_records and not incomplete,
        },
    }
    operational_errors = []
    if invalid_records:
        operational_errors.append({
            "tool": "dock6",
            "reason": (
                "DOCK 6 emitted invalid pose record(s): "
                + json.dumps(invalid_records, ensure_ascii=False)
            ),
        })
    if incomplete:
        operational_errors.append({
            "tool": "dock6",
            "reason": (
                "DOCK 6 output is incomplete; no finite parseable pose for "
                f"input molecule(s) {incomplete}"
            ),
        })
    if operational_errors:
        merged["operational_errors"] = operational_errors
    return merged


def _parse_headers(chunk: str) -> dict[str, Any]:
    headers: dict[str, Any] = {}
    for line in chunk.splitlines():
        match = _HEADER.match(line)
        if match is None:
            continue
        name = match.group("name").strip()
        value = match.group("value").strip()
        try:
            number = float(value)
            parsed: Any = number if math.isfinite(number) else value
        except ValueError:
            parsed = value
        if name in headers and headers[name] != parsed:
            # Repeated descriptors such as Heavy_Atoms are harmless only when
            # their values agree.
            raise ValueError(f"conflicting repeated DOCK 6 header {name!r}")
        headers[name] = parsed
    return headers


def _finite_header(
    headers: dict[str, Any], name: str, *, required: bool,
) -> float | None:
    if name not in headers:
        if required:
            raise ValueError(f"DOCK 6 pose is missing native {name}")
        return None
    try:
        value = float(headers[name])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"DOCK 6 {name} is not numeric") from exc
    if not math.isfinite(value):
        raise ValueError(f"DOCK 6 {name} is not finite")
    return value


def score_metadata() -> dict:
    return {
        "dock6_grid_score": {
            "meaning": "DOCK 6 native grid score",
            "unit": "DOCK 6 grid-score units",
            "better": "lower",
            "cross_engine_comparable": False,
        },
        "dock6_grid_vdw_energy": {
            "meaning": "DOCK 6 native grid van der Waals component",
            "unit": "DOCK 6 grid-score units",
            "better": "context_only",
        },
        "dock6_grid_es_energy": {
            "meaning": "DOCK 6 native grid electrostatic component",
            "unit": "DOCK 6 grid-score units",
            "better": "context_only",
        },
        "engine_pose_rank": {
            "meaning": "One-based output order within one input ligand",
            "unit": "rank",
            "better": "context_only",
        },
    }


def _dock6_version(runtime: dict, environment: dict[str, str]) -> str:
    try:
        result = run_tool(
            runtime,
            ["--version"],
            timeout=30,
            check=False,
            env=environment,
        )
    except (OSError, ToolRuntimeError):
        return "unknown"
    output = "\n".join((result.get("stdout") or "", result.get("stderr") or ""))
    match = re.search(r"(?i)DOCK\s*(?:version)?\s*v?(\d+(?:\.\d+){1,3})", output)
    return match.group(1) if match else "unknown"


def _write_log(path: Path, process: dict, preparation: dict) -> None:
    path.write_text(
        "[ligand preparation]\n"
        + json.dumps(preparation, indent=2, ensure_ascii=False)
        + "\n[dock6 stdout]\n"
        + (process.get("stdout") or "")
        + ("\n[dock6 stderr]\n" if process.get("stderr") else "")
        + (process.get("stderr") or ""),
        encoding="utf-8",
    )


def _command_record(process: dict) -> dict:
    return {
        key: process.get(key)
        for key in ("command", "returncode", "elapsed_seconds", "environment")
    }


def _file_identity(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": path.stat().st_size,
    }

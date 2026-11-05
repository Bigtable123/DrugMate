"""Classic AutoDock4 and ADFR adapters for the Step 4 pose contract.

Both engines consume PDBQT, but their native preparation tools can reorder
atoms while building a torsion tree.  The shared parser therefore establishes
an audited mapping from prepared-PDBQT atom serials back to the source RDKit
molecule before any docking output is accepted.  Poses are published only
when every input molecule has at least one finite native score and every heavy
atom can be mapped one-to-one.
"""
from __future__ import annotations

import hashlib
import math
import re
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Callable

from docking_common import (
    annotate_canonical_pose,
    compound_name,
    copy_identity_properties,
    reference_ligand_box,
)
from gnina_utils import normalize_docking_box
from tool_runtime import ToolRuntimeError, run_tool


DEFAULT_AUTODOCK4_RUNTIME = {
    "kind": "composite",
    "components": {
        "io": {
            "kind": "command",
            "environment": "autodock4",
            "command": "python",
        },
        "grid": {
            "kind": "command",
            "environment": "autodock4",
            "command": "autogrid4",
        },
        "engine": {
            "kind": "command",
            "environment": "autodock4",
            "command": "autodock4",
        },
    },
}

DEFAULT_ADFR_RUNTIME = {
    "kind": "composite",
    "engine_version": "1.1",
    "components": {
        "io": {
            "kind": "command",
            "environment": "adfrsuite",
            "command": "prepare_receptor",
        },
        "ligand_io": {
            "kind": "command",
            "environment": "adfrsuite",
            "command": "prepare_ligand",
        },
        "grid": {
            "kind": "command",
            "environment": "adfrsuite",
            "command": "agfr",
        },
        "engine": {
            "kind": "command",
            "environment": "adfrsuite",
            "command": "adfr",
        },
    },
}

_AD4_SCORE = re.compile(
    r"Estimated\s+Free\s+Energy\s+of\s+Binding\s*=\s*(?P<value>\S+)",
    re.IGNORECASE,
)
_ADFR_SCORE = re.compile(r"\bFEB:\s*(?P<value>\S+)", re.IGNORECASE)
_AUTODOCK_RELEASE_VERSION = re.compile(
    r"\bAutoDock\s+v?\d+(?:\.\d+)+\s+Release\s+"
    r"v?(?P<version>\d+(?:\.\d+)+(?:[-\w.]*)?)",
    re.IGNORECASE,
)
_AUTODOCK_VERSION = re.compile(
    r"\bAutoDock(?:\s+Version)?\s+v?(?P<version>\d+(?:\.\d+)+(?:[-\w.]*)?)",
    re.IGNORECASE,
)
_ADFR_VERSION = re.compile(
    r"\bADFR\s+v?(?P<version>\d+(?:\.\d+)+(?:[-\w.]*)?)",
    re.IGNORECASE,
)

_AD_TYPE_TO_ELEMENT = {
    "A": "C",
    "C": "C",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "OA": "O",
    "OS": "O",
    "O": "O",
    "SA": "S",
    "S": "S",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "HD": "H",
    "HS": "H",
    "H": "H",
    "MG": "Mg",
    "MN": "Mn",
    "ZN": "Zn",
    "CA": "Ca",
    "FE": "Fe",
    "CU": "Cu",
}


def run_autodock4_docking(
    *,
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    runtime: dict | None = None,
    num_runs: int = 20,
    num_modes: int = 5,
    max_evals: int = 2_500_000,
    population_size: int = 150,
    cluster_rmsd: float = 2.0,
    spacing: float = 0.375,
    autobox_add: float = 4.0,
    seed: int = 42,
    timeout_seconds: int = 86400,
) -> dict:
    """Run ADT -> AutoGrid4 -> AutoDock4 for every input ligand."""
    try:
        params = _validated_parameters(
            engine="AutoDock4",
            num_runs=num_runs,
            num_modes=num_modes,
            max_evals=max_evals,
            population_size=population_size,
            spacing=spacing,
            autobox_add=autobox_add,
            seed=seed,
            timeout_seconds=timeout_seconds,
            cluster_rmsd=cluster_rmsd,
        )
        resolved_box = _resolve_box(
            reference_ligand_file, docking_box, params["autobox_add"],
        )
        grid_points = _autogrid_points(resolved_box["size"], params["spacing"])
        components = _runtime_components(
            runtime or DEFAULT_AUTODOCK4_RUNTIME, "AutoDock4",
        )
        molecules = _load_input_molecules(Path(ligand_file))
        _require_nonempty_file(Path(receptor_file), "prepared receptor")
    except Exception as exc:
        return {"success": False, "engine": "autodock4", "error": str(exc)}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "autodock4_native"
    native_root.mkdir(parents=True, exist_ok=True)
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    log_file = output_dir / "autodock4.log"
    docked_sdf = output_dir / "docked_poses.sdf"
    docked_sdf.unlink(missing_ok=True)
    stages: list[dict] = []
    execute = _stage_executor(
        deadline=time.monotonic() + params["timeout_seconds"],
        stage_records=stages,
    )

    prepared_receptor = native_output_dir / "receptor.pdbqt"
    entries: list[dict] = []
    ligand_protocols: list[dict] = []
    try:
        execute(
            "adt_receptor_preparation",
            components["io"],
            ["-r", str(receptor_file), "-o", str(prepared_receptor)],
            command="prepare_receptor4",
            cwd=native_output_dir,
        )
        _require_nonempty_file(prepared_receptor, "ADT receptor PDBQT")

        for source_index, molecule in enumerate(molecules):
            token = f"drugmate_src_{source_index:06d}"
            record_dir = native_output_dir / token
            record_dir.mkdir()
            ligand_pdb = record_dir / f"{token}.pdb"
            ligand_pdbqt = record_dir / f"{token}.pdbqt"
            local_receptor = record_dir / "receptor.pdbqt"
            _write_ligand_pdb(molecule, ligand_pdb)
            shutil.copy2(prepared_receptor, local_receptor)

            execute(
                f"adt_ligand_preparation:{source_index}",
                components["io"],
                ["-l", str(ligand_pdb), "-o", str(ligand_pdbqt)],
                command="prepare_ligand4",
                cwd=record_dir,
            )
            _require_nonempty_file(ligand_pdbqt, "ADT ligand PDBQT")
            atom_map = _prepared_atom_map(molecule, ligand_pdbqt)

            gpf = record_dir / "autodock4.gpf"
            dpf = record_dir / "autodock4.dpf"
            center = resolved_box["center"]
            execute(
                f"adt_grid_parameter_preparation:{source_index}",
                components["io"],
                [
                    "-l", str(ligand_pdbqt),
                    "-r", str(local_receptor),
                    "-o", str(gpf),
                    "-p", "gridcenter=" + ",".join(_format_number(v) for v in center),
                    "-p", "npts=" + ",".join(str(v) for v in grid_points),
                    "-p", f"spacing={_format_number(params['spacing'])}",
                ],
                command="prepare_gpf4",
                cwd=record_dir,
            )
            _require_nonempty_file(gpf, "AutoDock4 grid parameter file")

            execute(
                f"adt_docking_parameter_preparation:{source_index}",
                components["io"],
                [
                    "-l", str(ligand_pdbqt),
                    "-r", str(local_receptor),
                    "-o", str(dpf),
                    "-p", f"ga_run={params['num_runs']}",
                    "-p", f"ga_num_evals={params['max_evals']}",
                    "-p", f"ga_pop_size={params['population_size']}",
                    "-p", f"rmstol={_format_number(params['cluster_rmsd'])}",
                ],
                command="prepare_dpf42",
                cwd=record_dir,
            )
            _require_nonempty_file(dpf, "AutoDock4 docking parameter file")
            _replace_dpf_seed(dpf, params["seed"] + source_index * 2)

            grid_log = record_dir / "autogrid4.glg"
            execute(
                f"autogrid4:{source_index}",
                components["grid"],
                ["-p", str(gpf), "-l", str(grid_log)],
                cwd=record_dir,
            )
            _require_nonempty_file(grid_log, "AutoGrid4 log")
            map_files = sorted(record_dir.glob("*.map"))
            field_files = sorted(record_dir.glob("*.maps.fld"))
            if not map_files or not field_files:
                raise ValueError(
                    f"AutoGrid4 produced an incomplete map bundle for input {source_index}"
                )
            for path in [*map_files, *field_files]:
                _require_nonempty_file(path, "AutoGrid4 map asset")

            dlg = record_dir / "autodock4.dlg"
            execute(
                f"autodock4:{source_index}",
                components["engine"],
                ["-p", str(dpf), "-l", str(dlg)],
                cwd=record_dir,
            )
            _require_nonempty_file(dlg, "AutoDock4 docking log")
            entries.append({
                "source_molecule_index": source_index,
                "source_molecule": molecule,
                "atom_map": atom_map,
                "pose_file": dlg,
            })
            ligand_protocols.append({
                "source_molecule_index": source_index,
                "prepared_ligand": _file_identity(ligand_pdbqt),
                "grid_parameters": _file_identity(gpf),
                "docking_parameters": _file_identity(dpf),
                "grid_field": [_file_identity(path) for path in field_files],
                "affinity_maps": [_file_identity(path) for path in map_files],
            })

        version = _engine_version(
            "autodock4", stages, [entry["pose_file"] for entry in entries],
            runtime or DEFAULT_AUTODOCK4_RUNTIME,
        )
        merge = _merge_pdbqt_outputs(
            entries=entries,
            engine="autodock4",
            output_sdf=docked_sdf,
            max_poses_per_molecule=params["num_modes"],
            method_version=version,
        )
    except Exception as exc:
        return _failure(
            engine="autodock4",
            error=f"could not execute AutoDock4 adapter: {exc}",
            log_file=log_file,
            native_output_dir=native_output_dir,
            stages=stages,
            docked_sdf=docked_sdf,
            search_space=resolved_box,
        )

    _write_log(log_file, stages)
    return {
        "success": True,
        "engine": "autodock4",
        "docked_sdf": str(docked_sdf),
        "scores": merge["scores"],
        "score_metadata": autodock4_score_metadata(),
        "protocol_assets": {
            "prepared_receptor": _file_identity(prepared_receptor),
            "ligand_protocols": ligand_protocols,
        },
        "software": {
            "engine": "AutoDock4",
            "engine_version": version,
            "runtime_environment": components["engine"].get(
                "environment", "autodock4",
            ),
            "grid_engine": "AutoGrid4",
            "preparation": "AutoDockTools",
        },
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "num_poses": len(merge["scores"]),
        "num_input_molecules": len(molecules),
        "search_space": _search_space_record(resolved_box, reference_ligand_file),
        "redocking_validation": _not_applicable_redocking("AutoDock4"),
        "quality_gate_passed": None,
        "commands": [_command_record(stage) for stage in stages],
        "merge_diagnostics": merge["diagnostics"],
        # Carried through so the failure registry can name the inputs this
        # engine covered no pose for; empty for a complete batch.
        "operational_errors": merge.get("operational_errors", []),
    }


def run_adfr_docking(
    *,
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    runtime: dict | None = None,
    num_runs: int = 20,
    num_modes: int = 5,
    max_evals: int = 2_500_000,
    population_size: int = 150,
    max_cores: int = 1,
    spacing: float = 0.375,
    autobox_add: float = 4.0,
    seed: int = 42,
    timeout_seconds: int = 86400,
) -> dict:
    """Run ligand/receptor preparation -> AGFR -> rigid-receptor ADFR."""
    try:
        params = _validated_parameters(
            engine="ADFR",
            num_runs=num_runs,
            num_modes=num_modes,
            max_evals=max_evals,
            population_size=population_size,
            spacing=spacing,
            autobox_add=autobox_add,
            seed=seed,
            timeout_seconds=timeout_seconds,
            max_cores=max_cores,
        )
        resolved_box = _resolve_box(
            reference_ligand_file, docking_box, params["autobox_add"],
        )
        components = _runtime_components(runtime or DEFAULT_ADFR_RUNTIME, "ADFR")
        molecules = _load_input_molecules(Path(ligand_file))
        _require_nonempty_file(Path(receptor_file), "prepared receptor")
    except Exception as exc:
        return {"success": False, "engine": "adfr", "error": str(exc)}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "adfr_native"
    native_root.mkdir(parents=True, exist_ok=True)
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    log_file = output_dir / "adfr.log"
    docked_sdf = output_dir / "docked_poses.sdf"
    docked_sdf.unlink(missing_ok=True)
    stages: list[dict] = []
    execute = _stage_executor(
        deadline=time.monotonic() + params["timeout_seconds"],
        stage_records=stages,
    )

    prepared_receptor = native_output_dir / "receptor.pdbqt"
    prepared_entries: list[dict] = []
    try:
        execute(
            "adfr_receptor_preparation",
            components["io"],
            ["-r", str(receptor_file), "-o", str(prepared_receptor)],
            command="prepare_receptor",
            cwd=native_output_dir,
        )
        _require_nonempty_file(prepared_receptor, "ADFR receptor PDBQT")

        ligand_root = native_output_dir / "ligands"
        ligand_root.mkdir()
        for source_index, molecule in enumerate(molecules):
            token = f"drugmate_src_{source_index:06d}"
            record_dir = ligand_root / token
            record_dir.mkdir()
            ligand_pdb = record_dir / f"{token}.pdb"
            ligand_pdbqt = record_dir / f"{token}.pdbqt"
            _write_ligand_pdb(molecule, ligand_pdb)
            execute(
                f"adfr_ligand_preparation:{source_index}",
                components["ligand_io"],
                ["-l", str(ligand_pdb), "-o", str(ligand_pdbqt)],
                command="prepare_ligand",
                cwd=record_dir,
            )
            _require_nonempty_file(ligand_pdbqt, "ADFR ligand PDBQT")
            prepared_entries.append({
                "source_molecule_index": source_index,
                "source_molecule": molecule,
                "record_dir": record_dir,
                "ligand_pdbqt": ligand_pdbqt,
                "atom_map": _prepared_atom_map(molecule, ligand_pdbqt),
            })

        target_file = native_output_dir / "receptor_target.trg"
        center = resolved_box["center"]
        size = resolved_box["size"]
        execute(
            "agfr_target_generation",
            components["grid"],
            [
                "-r", str(prepared_receptor),
                "-b", "user",
                *(_format_number(value) for value in [*center, *size]),
                "-m", "all",
                "-s", _format_number(params["spacing"]),
                "-o", str(target_file),
            ],
            cwd=native_output_dir,
        )
        _require_nonempty_file(target_file, "AGFR target archive")

        entries: list[dict] = []
        for prepared in prepared_entries:
            source_index = prepared["source_molecule_index"]
            record_dir = prepared["record_dir"]
            job_name = f"drugmate_src_{source_index:06d}_adfr"
            output_base = record_dir / job_name
            execute(
                f"adfr_docking:{source_index}",
                components["engine"],
                [
                    "-l", str(prepared["ligand_pdbqt"]),
                    "-t", str(target_file),
                    "-J", job_name,
                    "-c", str(params["max_cores"]),
                    "-o", str(output_base),
                    "-n", str(params["num_runs"]),
                    "-p", str(params["population_size"]),
                    "-e", str(params["max_evals"]),
                    "-S", str(params["seed"] + source_index * params["num_runs"]),
                    "-T",
                    "-O",
                ],
                cwd=record_dir,
            )
            # ADFR treats -o as an output *base*, then appends these native
            # suffixes; it does not treat -o as the summary filename despite
            # the wording in the command-line help.
            pose_file = Path(f"{output_base}_out.pdbqt")
            summary_file = Path(f"{output_base}_summary.dlg")
            _require_nonempty_file(pose_file, "ADFR pose PDBQT")
            _require_nonempty_file(summary_file, "ADFR summary DLG")
            entries.append({
                "source_molecule_index": source_index,
                "source_molecule": prepared["source_molecule"],
                "atom_map": prepared["atom_map"],
                "pose_file": pose_file,
            })

        version = _engine_version(
            "adfr", stages, [entry["pose_file"] for entry in entries],
            runtime or DEFAULT_ADFR_RUNTIME,
        )
        merge = _merge_pdbqt_outputs(
            entries=entries,
            engine="adfr",
            output_sdf=docked_sdf,
            max_poses_per_molecule=params["num_modes"],
            method_version=version,
        )
    except Exception as exc:
        return _failure(
            engine="adfr",
            error=f"could not execute ADFR adapter: {exc}",
            log_file=log_file,
            native_output_dir=native_output_dir,
            stages=stages,
            docked_sdf=docked_sdf,
            search_space=resolved_box,
        )

    _write_log(log_file, stages)
    return {
        "success": True,
        "engine": "adfr",
        "docked_sdf": str(docked_sdf),
        "scores": merge["scores"],
        "score_metadata": adfr_score_metadata(),
        "protocol_assets": {
            "prepared_receptor": _file_identity(prepared_receptor),
            "target_archive": _file_identity(target_file),
            "prepared_ligands": [
                {
                    "source_molecule_index": entry["source_molecule_index"],
                    "asset": _file_identity(entry["ligand_pdbqt"]),
                }
                for entry in prepared_entries
            ],
        },
        "software": {
            "engine": "ADFR",
            "engine_version": version,
            "runtime_environment": components["engine"].get(
                "environment", "adfrsuite",
            ),
            "grid_engine": "AGFR",
            "preparation": "ADFRsuite prepare_receptor/prepare_ligand",
            "receptor_mode": "rigid",
        },
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "num_poses": len(merge["scores"]),
        "num_input_molecules": len(molecules),
        "search_space": _search_space_record(resolved_box, reference_ligand_file),
        "redocking_validation": _not_applicable_redocking("ADFR"),
        "quality_gate_passed": None,
        "commands": [_command_record(stage) for stage in stages],
        "merge_diagnostics": merge["diagnostics"],
        # Carried through so the failure registry can name the inputs this
        # engine covered no pose for; empty for a complete batch.
        "operational_errors": merge.get("operational_errors", []),
    }


def autodock4_score_metadata() -> dict:
    return {
        "autodock4_binding_energy": {
            "metric_id": "autodock4.binding_energy",
            "meaning": "AutoDock4 estimated free energy of binding",
            "unit": "kcal/mol",
            "better": "lower",
            "cross_engine_comparable": False,
        },
        "engine_pose_rank": {
            "meaning": "Rank by native AutoDock4 binding energy for one ligand",
            "unit": "rank",
            "better": "context_only",
        },
    }


def adfr_score_metadata() -> dict:
    return {
        "adfr_affinity": {
            "metric_id": "adfr.affinity",
            "meaning": "ADFR FEB affinity from the native solution record",
            "unit": "kcal/mol",
            "better": "lower",
            "cross_engine_comparable": False,
        },
        "engine_pose_rank": {
            "meaning": "Rank by native ADFR affinity for one ligand",
            "unit": "rank",
            "better": "context_only",
        },
    }


def _runtime_components(runtime: dict, label: str) -> dict[str, dict]:
    if runtime.get("kind") != "composite":
        raise ToolRuntimeError(
            f"{label} runtime must be a composite io + grid + engine declaration"
        )
    components = runtime.get("components")
    if not isinstance(components, dict):
        raise ToolRuntimeError(f"{label} composite runtime is missing components")
    resolved = {}
    required_components = (
        ("io", "ligand_io", "grid", "engine")
        if label == "ADFR"
        else ("io", "grid", "engine")
    )
    for name in required_components:
        component = components.get(name)
        if not isinstance(component, dict):
            raise ToolRuntimeError(f"{label} runtime requires component {name!r}")
        resolved[name] = dict(component)
    return resolved


def _validated_parameters(
    *,
    engine: str,
    num_runs: Any,
    num_modes: Any,
    max_evals: Any,
    population_size: Any,
    spacing: Any,
    autobox_add: Any,
    seed: Any,
    timeout_seconds: Any,
    cluster_rmsd: Any | None = None,
    max_cores: Any | None = None,
) -> dict:
    try:
        result = {
            "num_runs": int(num_runs),
            "num_modes": int(num_modes),
            "max_evals": int(max_evals),
            "population_size": int(population_size),
            "spacing": float(spacing),
            "autobox_add": float(autobox_add),
            "seed": int(seed),
            "timeout_seconds": int(timeout_seconds),
        }
        if cluster_rmsd is not None:
            result["cluster_rmsd"] = float(cluster_rmsd)
        if max_cores is not None:
            result["max_cores"] = int(max_cores)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {engine} parameter: {exc}") from exc

    positive = (
        "num_runs", "num_modes", "max_evals", "population_size",
        "timeout_seconds",
    )
    if any(result[name] < 1 for name in positive):
        raise ValueError(f"{engine} run/search/timeout parameters must be >= 1")
    if result["seed"] < 0:
        raise ValueError(f"{engine} seed must be >= 0")
    if not math.isfinite(result["spacing"]) or result["spacing"] <= 0:
        raise ValueError(f"{engine} spacing must be finite and > 0")
    if not math.isfinite(result["autobox_add"]) or result["autobox_add"] < 0:
        raise ValueError(f"{engine} autobox_add must be finite and >= 0")
    if "cluster_rmsd" in result and (
        not math.isfinite(result["cluster_rmsd"])
        or result["cluster_rmsd"] < 0
    ):
        raise ValueError("AutoDock4 cluster_rmsd must be finite and >= 0")
    if "max_cores" in result and result["max_cores"] < 1:
        raise ValueError("ADFR max_cores must be >= 1")
    return result


def _resolve_box(
    reference_ligand_file: Path | None,
    docking_box: dict | None,
    autobox_add: float,
) -> dict:
    if reference_ligand_file is not None:
        return reference_ligand_box(Path(reference_ligand_file), autobox_add)
    if docking_box is None:
        raise ValueError("reference_ligand_file or docking_box is required")
    return normalize_docking_box(docking_box)


def _autogrid_points(size: list[float], spacing: float) -> tuple[int, int, int]:
    points = []
    for dimension in size:
        value = int(math.ceil(float(dimension) / spacing))
        value = max(2, value + value % 2)
        if value > 126:
            raise ValueError(
                "AutoDock4 box exceeds AutoGrid4's 126-point limit at the "
                f"declared spacing ({dimension} A / {spacing} A)"
            )
        points.append(value)
    return tuple(points)  # type: ignore[return-value]


def _stage_executor(
    *, deadline: float, stage_records: list[dict],
) -> Callable[..., dict]:
    def execute(
        label: str,
        runtime: dict,
        args: list[str],
        *,
        command: str | None = None,
        cwd: Path | None = None,
    ) -> dict:
        remaining = int(math.ceil(deadline - time.monotonic()))
        if remaining < 1:
            raise TimeoutError(f"total timeout exhausted before {label}")
        process = run_tool(
            runtime,
            args,
            command=command,
            cwd=cwd,
            timeout=remaining,
            check=False,
        )
        record = {"stage": label, **process}
        stage_records.append(record)
        if not process.get("success"):
            detail = (
                process.get("error")
                or str(process.get("stderr") or "").strip()
                or str(process.get("stdout") or "").strip()
                or f"return code {process.get('returncode')}"
            )
            raise ToolRuntimeError(f"{label} failed: {detail}", process)
        return process

    return execute


def _load_input_molecules(path: Path) -> list[Any]:
    from rdkit import Chem

    _require_nonempty_file(path, "ligand SDF")
    supplier = Chem.SDMolSupplier(
        str(path), removeHs=False, strictParsing=False,
    )
    molecules: list[Any] = []
    invalid: list[int] = []
    for index, molecule in enumerate(supplier):
        if molecule is None:
            invalid.append(index)
            continue
        _validate_3d_molecule(molecule, index)
        molecules.append(Chem.Mol(molecule))
    if invalid:
        raise ValueError(f"input ligand SDF has unparsable record(s): {invalid}")
    if not molecules:
        raise ValueError("input ligand SDF contains no parseable molecules")
    return molecules


def _validate_3d_molecule(molecule: Any, index: int) -> None:
    if molecule.GetNumHeavyAtoms() < 1:
        raise ValueError(f"ligand record {index} has no heavy atoms")
    if molecule.GetNumConformers() != 1:
        raise ValueError(f"ligand record {index} must contain exactly one conformer")
    conformer = molecule.GetConformer()
    if not conformer.Is3D():
        raise ValueError(f"ligand record {index} conformer is not marked as 3D")
    for atom_index in range(molecule.GetNumAtoms()):
        position = conformer.GetAtomPosition(atom_index)
        if not all(
            math.isfinite(value) for value in (position.x, position.y, position.z)
        ):
            raise ValueError(
                f"ligand record {index} contains non-finite coordinates"
            )


def _write_ligand_pdb(molecule: Any, path: Path) -> None:
    from rdkit import Chem

    Chem.MolToPDBFile(molecule, str(path))
    _require_nonempty_file(path, "serialized ligand PDB")


def _prepared_atom_map(source_molecule: Any, prepared_pdbqt: Path) -> dict[int, int]:
    """Map prepared PDBQT heavy-atom serials to original RDKit atom indices."""
    atoms = [atom for atom in _pdbqt_atoms(prepared_pdbqt) if atom["element"] != "H"]
    conformer = source_molecule.GetConformer()
    source_heavy = [
        atom.GetIdx()
        for atom in source_molecule.GetAtoms()
        if atom.GetAtomicNum() > 1
    ]
    if len(atoms) != len(source_heavy):
        raise ValueError(
            "PDBQT preparation changed heavy-atom coverage "
            f"({len(atoms)} prepared != {len(source_heavy)} source)"
        )

    periodic_table = _periodic_table()
    unused = set(source_heavy)
    mapping: dict[int, int] = {}
    for prepared in atoms:
        atomic_number = periodic_table.GetAtomicNumber(prepared["element"])
        candidates = []
        for source_index in unused:
            source_atom = source_molecule.GetAtomWithIdx(source_index)
            if source_atom.GetAtomicNum() != atomic_number:
                continue
            position = conformer.GetAtomPosition(source_index)
            distance = math.sqrt(sum(
                (left - right) ** 2
                for left, right in zip(
                    prepared["coordinates"],
                    (position.x, position.y, position.z),
                    strict=True,
                )
            ))
            if distance <= 0.05:
                candidates.append((distance, source_index))
        candidates.sort()
        if not candidates:
            raise ValueError(
                "could not map prepared PDBQT atom serial "
                f"{prepared['serial']} back to the source molecule"
            )
        if (
            len(candidates) > 1
            and math.isclose(candidates[0][0], candidates[1][0], abs_tol=1e-6)
        ):
            raise ValueError(
                f"ambiguous source mapping for PDBQT atom serial {prepared['serial']}"
            )
        source_index = candidates[0][1]
        if prepared["serial"] in mapping:
            raise ValueError(f"duplicate PDBQT atom serial {prepared['serial']}")
        mapping[prepared["serial"]] = source_index
        unused.remove(source_index)
    if unused:
        raise ValueError(f"prepared PDBQT omitted source heavy atoms {sorted(unused)}")
    return mapping


def _pdbqt_atoms(path: Path) -> list[dict]:
    _require_nonempty_file(path, "PDBQT file")
    atoms = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = _native_line(raw_line)
        if line.startswith(("ATOM  ", "HETATM")):
            atoms.append(_parse_atom_line(line))
    if not atoms:
        raise ValueError(f"PDBQT file contains no atom records: {path}")
    serials = [atom["serial"] for atom in atoms]
    if len(serials) != len(set(serials)):
        raise ValueError(f"PDBQT file contains duplicate atom serials: {path}")
    return atoms


def _parse_pdbqt_models(path: Path, engine: str) -> list[dict]:
    """Parse AD4 DLG or ADFR PDBQT models with their native affinity."""
    _require_nonempty_file(path, f"{engine} pose output")
    score_pattern = _AD4_SCORE if engine == "autodock4" else _ADFR_SCORE
    models: list[dict] = []
    current: dict | None = None
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if engine == "autodock4" and not raw_line.lstrip().startswith("DOCKED:"):
            continue
        line = _native_line(raw_line)
        if line.startswith("MODEL"):
            if current is not None:
                _finish_model(current, models, path)
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"invalid MODEL record in {path}: {line!r}")
            try:
                model_number = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"invalid MODEL number in {path}: {line!r}") from exc
            current = {"model_number": model_number, "atoms": [], "score": None}
            continue
        if current is None:
            continue
        score_match = score_pattern.search(line)
        if score_match is not None:
            value = _finite_float(score_match.group("value"), f"{engine} score")
            if current["score"] is not None and not math.isclose(
                current["score"], value, abs_tol=1e-8,
            ):
                raise ValueError(
                    f"conflicting native scores for MODEL {current['model_number']}"
                )
            current["score"] = value
        if line.startswith(("ATOM  ", "HETATM")):
            current["atoms"].append(_parse_atom_line(line))
        if line.startswith("ENDMDL"):
            _finish_model(current, models, path)
            current = None
    if current is not None:
        _finish_model(current, models, path)
    if not models:
        raise ValueError(f"{engine} output contains no scored PDBQT MODEL records")
    numbers = [model["model_number"] for model in models]
    if len(numbers) != len(set(numbers)):
        raise ValueError(f"{engine} output contains duplicate MODEL numbers")
    return models


def _finish_model(model: dict, models: list[dict], path: Path) -> None:
    if model["score"] is None:
        raise ValueError(
            f"MODEL {model['model_number']} in {path} has no native score"
        )
    if not model["atoms"]:
        raise ValueError(f"MODEL {model['model_number']} in {path} has no atoms")
    serials = [atom["serial"] for atom in model["atoms"]]
    if len(serials) != len(set(serials)):
        raise ValueError(
            f"MODEL {model['model_number']} in {path} has duplicate atom serials"
        )
    models.append(model)


def _parse_atom_line(line: str) -> dict:
    try:
        serial = int(line[6:11])
        coordinates = tuple(float(line[start:end]) for start, end in (
            (30, 38), (38, 46), (46, 54),
        ))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid PDBQT atom record: {line!r}") from exc
    if not all(math.isfinite(value) for value in coordinates):
        raise ValueError(f"non-finite PDBQT coordinates for atom serial {serial}")
    fields = line.split()
    atom_type = fields[-1] if fields else ""
    element = _element_from_autodock_type(atom_type, line[12:16])
    return {
        "serial": serial,
        "element": element,
        "coordinates": coordinates,
    }


def _element_from_autodock_type(atom_type: str, atom_name: str) -> str:
    normalized = re.sub(r"[^A-Za-z]", "", atom_type).upper()
    if normalized in _AD_TYPE_TO_ELEMENT:
        return _AD_TYPE_TO_ELEMENT[normalized]
    name = re.sub(r"[^A-Za-z]", "", atom_name)
    for length in (2, 1):
        candidate = name[:length].capitalize()
        if candidate and _periodic_table().GetAtomicNumber(candidate) > 0:
            return candidate
    raise ValueError(
        f"unsupported AutoDock atom type {atom_type!r} for atom {atom_name!r}"
    )


def _merge_pdbqt_outputs(
    *,
    entries: list[dict],
    engine: str,
    output_sdf: Path,
    max_poses_per_molecule: int,
    method_version: str,
) -> dict:
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    if engine not in {"autodock4", "adfr"}:
        raise ValueError(f"unsupported AutoDock-suite engine {engine!r}")
    partial = output_sdf.with_suffix(".partial.sdf")
    partial.unlink(missing_ok=True)
    writer = Chem.SDWriter(str(partial))
    scores: list[dict] = []
    diagnostics = {"molecules": []}
    covered: set[int] = set()
    try:
        for entry in sorted(entries, key=lambda value: value["source_molecule_index"]):
            source_index = int(entry["source_molecule_index"])
            try:
                source_molecule = entry["source_molecule"]
                compound_id = compound_name(source_molecule, source_index)
                models = _parse_pdbqt_models(Path(entry["pose_file"]), engine)
                selected = sorted(
                    models,
                    key=lambda model: (model["score"], model["model_number"]),
                )[:max_poses_per_molecule]
                if not selected:
                    raise ValueError(f"{engine} produced no pose for input {source_index}")
                atom_map = entry["atom_map"]
                expected_serials = set(atom_map)
                diagnostics["molecules"].append({
                    "source_molecule_index": source_index,
                    "compound_id": compound_id,
                    "native_models": len(models),
                    "published_models": len(selected),
                })
                source_to_output = {}
                next_output_index = 0
                for atom in source_molecule.GetAtoms():
                    if atom.GetAtomicNum() > 1:
                        source_to_output[atom.GetIdx()] = next_output_index
                        next_output_index += 1

                for engine_pose_rank, model in enumerate(selected, start=1):
                    pose_atoms = [
                        atom for atom in model["atoms"] if atom["element"] != "H"
                    ]
                    pose_serials = {atom["serial"] for atom in pose_atoms}
                    if pose_serials != expected_serials:
                        raise ValueError(
                            f"{engine} MODEL {model['model_number']} has incomplete or "
                            "unexpected heavy-atom serial coverage"
                        )
                    molecule = Chem.RemoveHs(Chem.Mol(source_molecule), sanitize=False)
                    if molecule.GetNumAtoms() != len(source_to_output):
                        raise ValueError("could not construct the source heavy-atom molecule")
                    conformer = molecule.GetConformer()
                    for atom in pose_atoms:
                        source_atom_index = atom_map[atom["serial"]]
                        output_atom_index = source_to_output[source_atom_index]
                        expected_element = molecule.GetAtomWithIdx(
                            output_atom_index,
                        ).GetSymbol()
                        if atom["element"] != expected_element:
                            raise ValueError(
                                f"{engine} MODEL {model['model_number']} atom serial "
                                f"{atom['serial']} changed element "
                                f"({atom['element']} != {expected_element})"
                            )
                        conformer.SetAtomPosition(
                            output_atom_index, Point3D(*atom["coordinates"]),
                        )
                    copy_identity_properties(source_molecule, molecule)
                    pose_index = len(scores)
                    native_score = _finite_float(model["score"], f"{engine} affinity")
                    flat_field = (
                        "autodock4_binding_energy"
                        if engine == "autodock4"
                        else "adfr_affinity"
                    )
                    metric_id = (
                        "autodock4.binding_energy"
                        if engine == "autodock4"
                        else "adfr.affinity"
                    )
                    molecule.SetDoubleProp(flat_field, native_score)
                    molecule.SetProp("score_direction", "lower")
                    molecule.SetIntProp("native_model_number", model["model_number"])
                    identity = annotate_canonical_pose(
                        molecule,
                        compound_id=compound_id,
                        pose_index=pose_index,
                        source_molecule_index=source_index,
                        engine=engine,
                        engine_pose_index=model["model_number"] - 1,
                        engine_pose_rank=engine_pose_rank,
                    )
                    structured_score = {
                        "metric_id": metric_id,
                        "value": native_score,
                        "unit": "kcal/mol",
                        "direction": "lower",
                        "method": engine,
                        "method_version": method_version,
                        "stage": "docking",
                        "applicability": "applicable",
                        "reliability": "supported",
                        "primary": True,
                        "pose_index": pose_index,
                        "pose_id": identity["pose_id"],
                    }
                    writer.write(molecule)
                    scores.append({
                        "pose_index": pose_index,
                        "source_molecule_index": source_index,
                        "engine_pose_index": model["model_number"] - 1,
                        "engine_pose_rank": engine_pose_rank,
                        "native_model_number": model["model_number"],
                        "compound_id": compound_id,
                        "mol_name": compound_id,
                        "engine": engine,
                        "method": engine,
                        "method_version": method_version,
                        "stage": "docking",
                        "score_reliability": "supported",
                        flat_field: native_score,
                        "pose_scores": [structured_score],
                        **identity,
                    })
                covered.add(source_index)
            except Exception as exc:
                # One input whose pose file the engine wrote unusably is that
                # candidate's missing value, not the batch's failure.  Letting it
                # reach the outer handler discarded every pose already merged, so a
                # single corrupt MODEL ended docking for the whole cohort -- while
                # The executable policy scopes a single tool fault to the single candidate.
                # The molecule lands in ``failed_molecules`` and in
                # ``uncovered_source_molecule_indices`` below, so it carries no pose
                # and coverage reads it as missing evidence.
                diagnostics.setdefault("failed_molecules", []).append({
                    "source_molecule_index": source_index,
                    "reason": str(exc),
                })
                continue
    except Exception:
        writer.close()
        partial.unlink(missing_ok=True)
        raise
    else:
        writer.close()

    expected = {int(entry["source_molecule_index"]) for entry in entries}
    if not scores:
        partial.unlink(missing_ok=True)
        raise ValueError(
            f"{engine} produced no parseable scored poses for any of the "
            f"{len(expected)} input molecules"
        )
    # An input this engine covered no pose for is that molecule's missing value,
    # not the batch's failure.  Raising on ``covered != expected`` discarded every
    # pose the engine did produce, so one ligand it could not place ended docking
    # for the whole cohort -- while the executable policy scopes a single tool fault to the
    # single candidate, and the executable policy drops a missing value from numerator and
    # denominator instead.  The gap stays named in ``diagnostics`` and is declared
    # in ``operational_errors`` for the failure registry.
    uncovered = sorted(expected - covered)
    diagnostics["covered_source_molecule_indices"] = sorted(covered)
    diagnostics["uncovered_source_molecule_indices"] = uncovered
    diagnostics["complete"] = not uncovered
    partial.replace(output_sdf)
    merged = {"success": True, "scores": scores, "diagnostics": diagnostics}
    if uncovered:
        merged["operational_errors"] = [{
            "tool": engine,
            "reason": (
                f"{engine} output is incomplete: covered inputs "
                f"{sorted(covered)}, expected {sorted(expected)}"
            ),
        }]
    return merged


def _replace_dpf_seed(path: Path, seed: int) -> None:
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    replacement = f"seed {int(seed)} {int(seed) + 1}"
    replaced = False
    output = []
    for line in lines:
        if line.strip().startswith("seed ") and not replaced:
            comment = ""
            if "#" in line:
                comment = " #" + line.split("#", 1)[1]
            output.append(replacement + comment)
            replaced = True
        else:
            output.append(line)
    if not replaced:
        output.insert(0, replacement)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def _engine_version(
    engine: str,
    stages: list[dict],
    native_outputs: list[Path],
    runtime: dict,
) -> str:
    declared = runtime.get("engine_version") or runtime.get("version")
    if declared not in (None, ""):
        return str(declared)
    text = "\n".join(
        str(stage.get(field) or "")
        for stage in stages
        for field in ("stdout", "stderr")
    )
    for path in native_outputs:
        try:
            text += "\n" + path.read_text(
                encoding="utf-8", errors="replace",
            )[:200_000]
        except OSError:
            continue
    if engine == "autodock4":
        # AutoDock4 DLGs advertise both the parameter-family version and the
        # executable release (for example ``4.2 Release 4.2.6``).  The latter
        # is the software identity; taking the first dotted token loses the
        # patch release and can merge scientifically distinct installations.
        match = _AUTODOCK_RELEASE_VERSION.search(text)
        if match is None:
            match = _AUTODOCK_VERSION.search(text)
    else:
        match = _ADFR_VERSION.search(text)
    return match.group("version") if match else "unknown"


def _failure(
    *,
    engine: str,
    error: str,
    log_file: Path,
    native_output_dir: Path,
    stages: list[dict],
    docked_sdf: Path,
    search_space: dict,
) -> dict:
    docked_sdf.unlink(missing_ok=True)
    _write_log(log_file, stages)
    return {
        "success": False,
        "engine": engine,
        "error": str(error),
        "log_file": str(log_file),
        "native_output_dir": str(native_output_dir),
        "search_space": search_space,
        "commands": [_command_record(stage) for stage in stages],
    }


def _write_log(path: Path, stages: list[dict]) -> None:
    chunks = []
    for stage in stages:
        chunks.append(f"[{stage.get('stage', 'unknown')}]\n")
        chunks.append("command: " + " ".join(
            str(value) for value in (stage.get("command") or [])
        ) + "\n")
        if stage.get("stdout"):
            chunks.append(str(stage["stdout"]))
            if not str(stage["stdout"]).endswith("\n"):
                chunks.append("\n")
        if stage.get("stderr"):
            chunks.append("[stderr]\n" + str(stage["stderr"]))
            if not str(stage["stderr"]).endswith("\n"):
                chunks.append("\n")
    path.write_text("".join(chunks), encoding="utf-8")


def _command_record(stage: dict) -> dict:
    return {
        key: stage.get(key)
        for key in (
            "stage", "command", "returncode", "elapsed_seconds", "environment",
        )
    }


def _file_identity(path: Path) -> dict:
    _require_nonempty_file(path, "protocol asset")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": path.stat().st_size,
    }


def _require_nonempty_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {label}: {path}")


def _finite_float(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {label}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"non-finite {label}: {value!r}")
    return number


def _native_line(raw_line: str) -> str:
    stripped = raw_line.lstrip()
    if stripped.startswith("DOCKED:"):
        return stripped[len("DOCKED:"):].lstrip()
    return raw_line


def _format_number(value: Any) -> str:
    number = _finite_float(value, "numeric parameter")
    return f"{number:.6f}".rstrip("0").rstrip(".")


def _periodic_table():
    from rdkit import Chem

    return Chem.GetPeriodicTable()


def _search_space_record(
    resolved_box: dict, reference_ligand_file: Path | None,
) -> dict:
    return {
        "mode": (
            "reference_autobox"
            if reference_ligand_file is not None
            else "explicit_box"
        ),
        **resolved_box,
    }


def _not_applicable_redocking(engine: str) -> dict:
    return {
        "performed": False,
        "passed": None,
        "applicable": False,
        "not_applicable": True,
        "reason": f"{engine} route has no registered engine-neutral redocking protocol",
        "gate": "not_applicable",
    }

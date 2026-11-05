"""rDock adapter for the canonical Step 4 docking contract."""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import tempfile
from pathlib import Path

from docking_common import (
    annotate_canonical_pose,
    compound_name,
    copy_identity_properties,
    reference_ligand_box,
)
from gnina_utils import normalize_docking_box
from tool_runtime import ToolRuntimeError, environment_spec, run_tool


DEFAULT_RUNTIME = {
    "kind": "command",
    "environment": "rdock",
    "command": "rbdock",
}


def _registered_rbt_root(runtime: dict) -> Path:
    """Derive rDock data files from the central environment declaration."""
    env = environment_spec(runtime)
    prefix = env.get("prefix")
    if env.get("kind") != "conda" or not prefix:
        raise ToolRuntimeError(
            "rDock runtime requires a prefix-registered Conda environment"
        )
    return (Path(str(prefix)) / "share" / "rDock").resolve()


def run_rdock_docking(
    *,
    receptor_file: Path,
    ligand_file: Path,
    reference_ligand_file: Path | None,
    output_dir: Path,
    docking_box: dict | None,
    runtime: dict | None = None,
    cavity_radius: float = 6.0,
    num_runs: int = 25,
    num_modes: int = 5,
    receptor_flex: float = 3.0,
    seed: int = 42,
    timeout_seconds: int = 86400,
) -> dict:
    """Run rDock with a ligand-defined or explicit sphere-mapped cavity."""
    try:
        num_runs = int(num_runs)
        num_modes = int(num_modes)
        seed = int(seed)
        timeout_seconds = int(timeout_seconds)
        receptor_flex = float(receptor_flex)
        cavity_radius = float(cavity_radius)
        if min(num_runs, num_modes, timeout_seconds) < 1:
            raise ValueError("num_runs, num_modes, and timeout_seconds must be >= 1")
        if seed < 0:
            raise ValueError("seed must be >= 0")
        if not math.isfinite(receptor_flex) or receptor_flex < 0:
            raise ValueError("receptor_flex must be a finite number >= 0")
        if not math.isfinite(cavity_radius) or cavity_radius <= 0:
            raise ValueError("cavity_radius must be a finite number > 0")
        resolved_box = (
            reference_ligand_box(reference_ligand_file, cavity_radius)
            if reference_ligand_file is not None
            else normalize_docking_box(docking_box)
        )
        if reference_ligand_file is not None:
            # This box is an audit envelope derived with the same numerical
            # radius; the native search contract is the ligand-site mapper
            # recorded below, not an AutoDock-style autobox.
            resolved_box.pop("autobox_add", None)
            resolved_box["source"] = "reference_ligand_cavity_envelope"
    except (TypeError, ValueError) as exc:
        return {"success": False, "error": f"invalid rDock parameters: {exc}"}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    native_root = output_dir / "rdock_native"
    native_root.mkdir(parents=True, exist_ok=True)
    native_output_dir = Path(tempfile.mkdtemp(prefix="run_", dir=native_root))
    log_file = output_dir / "rdock.log"
    docked_sdf = output_dir / "docked_poses.sdf"
    docked_sdf.unlink(missing_ok=True)
    resolved_runtime = dict(runtime or DEFAULT_RUNTIME)
    try:
        rbt_root = (
            Path(str(resolved_runtime["rbt_root"])).resolve()
            if resolved_runtime.get("rbt_root")
            else _registered_rbt_root(resolved_runtime)
        )
    except ToolRuntimeError as exc:
        return {
            "success": False,
            "error": f"rDock runtime resolution failed: {exc}",
            "native_output_dir": str(native_output_dir),
        }
    protocol_file = rbt_root / "data" / "scripts" / "dock.prm"
    if not protocol_file.is_file():
        return {
            "success": False,
            "error": f"rDock protocol file is missing: {protocol_file}",
            "native_output_dir": str(native_output_dir),
        }
    environment = {
        "RBT_ROOT": str(rbt_root),
        "RBT_HOME": str(native_output_dir),
    }
    stages: list[dict] = []

    try:
        staged_ligands = native_output_dir / "ligands.sdf"
        molecules = _stage_ligands(Path(ligand_file), staged_ligands)
        staged_receptor = native_output_dir / "receptor_input.pdb"
        shutil.copy2(Path(receptor_file), staged_receptor)
        receptor_mol2 = native_output_dir / "receptor.mol2"
        conversion = run_tool(
            resolved_runtime,
            [
                "-ipdb", str(staged_receptor),
                "-omol2", "-O", str(receptor_mol2),
                "--partialcharge", "gasteiger",
            ],
            command="obabel",
            cwd=native_output_dir,
            timeout=min(timeout_seconds, 1800),
            check=False,
            env=environment,
        )
        stages.append({"stage": "receptor_conversion", **conversion})
        if not conversion.get("success") or not receptor_mol2.is_file():
            raise ValueError("Open Babel could not create the rDock receptor MOL2")
        _canonicalize_mol2_title(receptor_mol2)

        receptor_prm = native_output_dir / "receptor.prm"
        staged_reference = None
        if reference_ligand_file is not None:
            staged_reference = native_output_dir / "reference_ligand.sdf"
            shutil.copy2(Path(reference_ligand_file), staged_reference)
        receptor_prm.write_text(
            _receptor_parameters(
                receptor_mol2=Path(receptor_mol2.name),
                reference_ligand_file=(
                    Path(staged_reference.name)
                    if staged_reference is not None
                    else None
                ),
                resolved_box=resolved_box,
                receptor_flex=receptor_flex,
                cavity_radius=cavity_radius,
            ),
            encoding="utf-8",
        )
        cavity = run_tool(
            resolved_runtime,
            [f"-r{receptor_prm}", "-was"],
            command="rbcavity",
            cwd=native_output_dir,
            timeout=min(timeout_seconds, 7200),
            check=False,
            env=environment,
        )
        stages.append({"stage": "cavity_generation", **cavity})
        if not cavity.get("success"):
            raise ValueError("rbcavity failed to generate the declared docking site")

        output_prefix = native_output_dir / "rdock"
        docking = run_tool(
            resolved_runtime,
            [
                "-i", str(staged_ligands),
                "-o", str(output_prefix),
                "-r", str(receptor_prm),
                "-p", str(protocol_file),
                "-n", str(num_runs),
                "-s", str(seed),
            ],
            cwd=native_output_dir,
            timeout=timeout_seconds,
            check=False,
            env=environment,
        )
        stages.append({"stage": "rdock", **docking})
        if not docking.get("success"):
            raise ValueError(f"rDock exited with code {docking.get('returncode')}")
        native_sdf = output_prefix.with_suffix(".sd")
        if not native_sdf.is_file() or native_sdf.stat().st_size == 0:
            raise ValueError("rDock completed without a non-empty output SD file")
        merge = _merge_outputs(
            native_sdf=native_sdf,
            input_molecules=molecules,
            output_sdf=docked_sdf,
            num_modes=num_modes,
        )
    except (OSError, ToolRuntimeError, ValueError) as exc:
        _write_log(log_file, stages)
        return {
            "success": False,
            "error": f"could not execute rDock adapter: {exc}",
            "log_file": str(log_file),
            "native_output_dir": str(native_output_dir),
            "commands": [_command_record(stage) for stage in stages],
        }

    _write_log(log_file, stages)
    return {
        "success": True,
        "engine": "rdock",
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
            "reason": "rDock route has no registered engine-neutral redocking protocol",
            "gate": "not_applicable",
        },
        "quality_gate_passed": None,
        "search_space": {
            "mode": (
                "reference_ligand_cavity"
                if reference_ligand_file is not None
                else "explicit_box_sphere_mapper"
            ),
            **resolved_box,
            "native_mapper": {
                "type": (
                    "RbtLigandSiteMapper"
                    if reference_ligand_file is not None
                    else "RbtSphereSiteMapper"
                ),
                "radius": (
                    cavity_radius
                    if reference_ligand_file is not None
                    else _containing_sphere_radius(resolved_box["size"])
                ),
            },
        },
        "software": {
            "engine": "rDock",
            "engine_version": _rdock_version(resolved_runtime, environment),
            "runtime_environment": resolved_runtime.get("environment", "rdock"),
        },
        "protocol_assets": {
            "receptor_mol2": _file_identity(receptor_mol2),
            "receptor_parameters": _file_identity(receptor_prm),
            "docking_protocol": _file_identity(protocol_file),
        },
        "commands": [_command_record(stage) for stage in stages],
        "merge_diagnostics": merge["diagnostics"],
    }


def _receptor_parameters(
    *,
    receptor_mol2: Path,
    reference_ligand_file: Path | None,
    resolved_box: dict,
    receptor_flex: float,
    cavity_radius: float,
) -> str:
    lines = [
        "RBT_PARAMETER_FILE_V1.00",
        "TITLE Drugmate rDock receptor",
        f"RECEPTOR_FILE {receptor_mol2}",
        f"RECEPTOR_FLEX {receptor_flex}",
        "",
        "SECTION MAPPER",
    ]
    if reference_ligand_file is not None:
        lines.extend([
            "    SITE_MAPPER RbtLigandSiteMapper",
            f"    REF_MOL {Path(reference_ligand_file)}",
            f"    RADIUS {cavity_radius}",
        ])
    else:
        center = [float(value) for value in resolved_box["center"]]
        radius = _containing_sphere_radius(resolved_box["size"])
        lines.extend([
            "    SITE_MAPPER RbtSphereSiteMapper",
            "    CENTER (%.6f,%.6f,%.6f)" % tuple(center),
            f"    RADIUS {radius:.6f}",
        ])
    lines.extend([
        "    SMALL_SPHERE 1.0",
        "    MIN_VOLUME 100",
        "    MAX_CAVITIES 1",
        "    VOL_INCR 0.0",
        "    GRIDSTEP 0.5",
        "END_SECTION",
        "",
        "SECTION CAVITY",
        "    SCORING_FUNCTION RbtCavityGridSF",
        "    WEIGHT 1.0",
        "END_SECTION",
        "",
    ])
    return "\n".join(lines)


def _containing_sphere_radius(size: list[float]) -> float:
    """Return the half-diagonal radius that contains an explicit box."""
    return math.sqrt(sum(float(value) ** 2 for value in size)) / 2.0


def _canonicalize_mol2_title(path: Path) -> None:
    """Remove Open Babel's input path from the protocol-defining MOL2 title."""
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    try:
        marker = lines.index("@<TRIPOS>MOLECULE")
    except ValueError as exc:
        raise ValueError("Open Babel receptor MOL2 has no molecule header") from exc
    if marker + 1 >= len(lines):
        raise ValueError("Open Babel receptor MOL2 has no molecule title")
    lines[marker + 1] = "Drugmate_receptor"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _stage_ligands(path: Path, output_path: Path) -> list:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(path), removeHs=False, strictParsing=False)
    molecules = []
    writer = Chem.SDWriter(str(output_path))
    invalid = []
    try:
        for index, molecule in enumerate(supplier):
            if molecule is None or molecule.GetNumConformers() != 1:
                invalid.append(index)
                continue
            staged = Chem.Mol(molecule)
            staged.SetProp("_Name", f"drugmate_src_{index:06d}")
            staged.SetIntProp("drugmate_source_molecule_index", index)
            writer.write(staged)
            molecules.append((index, molecule))
    finally:
        writer.close()
    if invalid:
        raise ValueError(f"rDock input contains invalid/non-3D records: {invalid}")
    if not molecules:
        raise ValueError("rDock input contains no valid 3D molecules")
    return molecules


def _merge_outputs(
    *,
    native_sdf: Path,
    input_molecules: list,
    output_sdf: Path,
    num_modes: int,
) -> dict:
    from rdkit import Chem

    source_by_token = {
        f"drugmate_src_{index:06d}": (index, molecule)
        for index, molecule in input_molecules
    }
    grouped: dict[str, list] = {token: [] for token in source_by_token}
    invalid_output = []
    supplier = Chem.SDMolSupplier(
        # rDock can emit chemically valid poses whose legacy MDL bond typing
        # fails modern RDKit sanitization (for example, a formally charged
        # nitrogen represented with four plain single bonds).  The native file
        # is coordinate/score evidence only; canonical chemistry comes from
        # the already-sanitized input molecule below.
        str(native_sdf), removeHs=False, sanitize=False, strictParsing=False,
    )
    for position, molecule in enumerate(supplier):
        if molecule is None:
            invalid_output.append(position)
            continue
        token = molecule.GetProp("_Name").strip() if molecule.HasProp("_Name") else ""
        if token not in grouped:
            raise ValueError(f"rDock returned unknown source molecule {token!r}")
        score = _finite_property(molecule, "SCORE")
        grouped[token].append((score, position, molecule))
    if invalid_output:
        raise ValueError(f"rDock returned invalid SD records: {invalid_output}")
    missing = sorted(token for token, poses in grouped.items() if not poses)
    if missing:
        raise ValueError(f"rDock omitted input molecules: {missing}")

    partial = output_sdf.with_suffix(".partial.sdf")
    partial.unlink(missing_ok=True)
    writer = Chem.SDWriter(str(partial))
    scores = []
    pose_index = 0
    try:
        for token, (source_index, source_molecule) in source_by_token.items():
            ranked = sorted(grouped[token], key=lambda item: (item[0], item[1]))
            for engine_pose_index, (score, _, molecule) in enumerate(
                ranked[:num_modes]
            ):
                output_molecule = _pose_on_source_topology(
                    source_molecule, molecule,
                )
                copy_identity_properties(source_molecule, output_molecule)
                compound_id = compound_name(source_molecule, source_index)
                identity = annotate_canonical_pose(
                    output_molecule,
                    compound_id=compound_id,
                    pose_index=pose_index,
                    source_molecule_index=source_index,
                    engine="rdock",
                    engine_pose_index=engine_pose_index,
                    engine_pose_rank=engine_pose_index + 1,
                )
                output_molecule.SetDoubleProp("rdock_score", score)
                score_record = {
                    "pose_index": pose_index,
                    "source_molecule_index": source_index,
                    "engine_pose_index": engine_pose_index,
                    "engine_pose_rank": engine_pose_index + 1,
                    "compound_id": compound_id,
                    "engine": "rdock",
                    "rdock_score": score,
                    "rdock_inter_score": _optional_finite_property(
                        molecule, "SCORE.INTER"
                    ),
                    "rdock_intra_score": _optional_finite_property(
                        molecule, "SCORE.INTRA"
                    ),
                    "rdock_restraint_score": _optional_finite_property(
                        molecule, "SCORE.RESTR"
                    ),
                    **identity,
                }
                score_record["pose_scores"] = [{
                    "metric_id": "rdock.score",
                    "value": score,
                    "unit": "rDock score",
                    "direction": "lower",
                    "method": "rdock",
                    "stage": "docking",
                    "pose_id": identity["pose_id"],
                    "pose_index": pose_index,
                    "applicability": "applicable",
                    "reliability": "supported",
                    "primary": True,
                }]
                for name, value in score_record.items():
                    if name.startswith("rdock_") and value is not None:
                        output_molecule.SetDoubleProp(name, float(value))
                writer.write(output_molecule)
                scores.append(score_record)
                pose_index += 1
    except Exception:
        writer.close()
        partial.unlink(missing_ok=True)
        raise
    else:
        writer.close()
    if not scores:
        partial.unlink(missing_ok=True)
        raise ValueError("rDock produced no finite scored poses")
    partial.replace(output_sdf)
    return {
        "scores": scores,
        "diagnostics": {
            "returned_records": sum(len(poses) for poses in grouped.values()),
            "exported_records": len(scores),
            "input_molecules": len(input_molecules),
        },
    }


def _pose_on_source_topology(source_molecule, native_pose):
    """Overlay native rDock coordinates on the canonical input topology."""
    from rdkit import Chem

    if source_molecule.GetNumAtoms() != native_pose.GetNumAtoms():
        raise ValueError(
            "rDock pose atom count does not match its source molecule: "
            f"{native_pose.GetNumAtoms()} != {source_molecule.GetNumAtoms()}"
        )
    source_elements = [atom.GetAtomicNum() for atom in source_molecule.GetAtoms()]
    native_elements = [atom.GetAtomicNum() for atom in native_pose.GetAtoms()]
    if source_elements != native_elements:
        raise ValueError("rDock pose atom order/elements do not match the input")
    if native_pose.GetNumConformers() != 1:
        raise ValueError("rDock pose is missing a single coordinate conformer")
    output_molecule = Chem.Mol(source_molecule)
    output_molecule.RemoveAllConformers()
    output_molecule.AddConformer(
        Chem.Conformer(native_pose.GetConformer()), assignId=True,
    )
    return output_molecule


def _finite_property(molecule, name: str) -> float:
    value = _optional_finite_property(molecule, name)
    if value is None:
        raise ValueError(f"rDock pose is missing finite {name}")
    return value


def _optional_finite_property(molecule, name: str) -> float | None:
    if not molecule.HasProp(name):
        return None
    try:
        value = float(molecule.GetProp(name))
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def score_metadata() -> dict:
    return {
        "rdock_score": {
            "meaning": "rDock total SCORE",
            "unit": "rDock score",
            "better": "lower",
        },
        "rdock_inter_score": {
            "meaning": "rDock intermolecular score component",
            "unit": "rDock score",
            "better": "lower",
        },
        "rdock_intra_score": {
            "meaning": "rDock intramolecular score component",
            "unit": "rDock score",
            "better": "context_only",
        },
        "rdock_restraint_score": {
            "meaning": "rDock restraint score component",
            "unit": "rDock score",
            "better": "context_only",
        },
    }


def _rdock_version(runtime: dict, environment: dict[str, str]) -> str:
    try:
        result = run_tool(
            runtime, ["-h"], timeout=20, check=False, env=environment,
        )
    except (OSError, ToolRuntimeError):
        return "unknown"
    text = f"{result.get('stdout', '')}\n{result.get('stderr', '')}"
    match = re.search(r"Executable:\s*rbdock\s*-\s*v([^\s]+)", text)
    return match.group(1) if match else "unknown"


def _file_identity(path: Path) -> dict:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"path": str(path), "sha256": digest, "size_bytes": path.stat().st_size}


def _command_record(stage: dict) -> dict:
    return {
        "stage": stage.get("stage"),
        "command": stage.get("command"),
        "returncode": stage.get("returncode"),
        "elapsed_seconds": stage.get("elapsed_seconds"),
        "environment": stage.get("environment"),
    }


def _write_log(path: Path, stages: list[dict]) -> None:
    payload = []
    for stage in stages:
        payload.append({
            **_command_record(stage),
            "stdout": stage.get("stdout", ""),
            "stderr": stage.get("stderr", ""),
        })
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

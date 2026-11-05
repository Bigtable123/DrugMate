"""Strict, planning-only TIES RBFE backend for Step 4.

This module deliberately does not execute molecular dynamics or free-energy
analysis.  It validates post-QC candidate identity/applicability, probes the
separate TIES runtime without modifying scientific inputs, and writes an
auditable transformation plan for a later TIES/TIES_MD execution stage.

``passed`` in this artifact means that a candidate or transformation is
eligible for planning.  It never means that a binding free energy was
calculated.  The manifest always records ``protocol_executed=false`` and does
not contain a numerical delta-delta-G field.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

from pose_identity import PoseIdentityError, build_pose_identity
from tool_runtime import ToolRuntimeError, run_tool
from candidate_topology import inspect_receptor_risks


SCHEMA_VERSION = "sbdd.ties-rbfe-plan.v1"
SUPPORTED_LIGAND_ELEMENTS = frozenset(
    {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
)
PROTEIN_RESNAMES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "ASH",
        "CYS",
        "CYM",
        "CYX",
        "GLN",
        "GLU",
        "GLH",
        "GLY",
        "HID",
        "HIE",
        "HIP",
        "HIS",
        "ILE",
        "LEU",
        "LYN",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)
REQUIRED_PYTHON_MODULES = ("ties", "TIES_MD", "openmm")
REQUIRED_API_CONTRACTS = (
    "ties.Pair",
    "ties.Config",
    "ties.Pair.superimpose",
    "TIES_MD.TIES.TIES",
    "TIES_MD.ties_analysis",
)
REQUIRED_COMMANDS = (
    "ties",
    "ties_md",
    "ties_ana",
    "antechamber",
    "parmchk2",
    "tleap",
    "sqm",
)
OPENMM_PLATFORMS = ("CUDA", "OpenCL", "CPU")
_CANONICAL_PREFIXES = {
    "pose_id": "pose:",
    "chemical_id": "chem:",
    "source_identity_id": "source:",
    "chemical_identity_id": "chemical-source:",
}


class TIESPlanningError(RuntimeError):
    """Raised for malformed calls, not per-candidate applicability failures."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _first_property(mol: Any, names: Iterable[str]) -> str | None:
    for name in names:
        try:
            if mol.HasProp(name):
                value = str(mol.GetProp(name)).strip()
                if value:
                    return value
        except Exception:
            continue
    return None


def _required_property(mol: Any, name: str) -> str:
    value = _first_property(mol, (name, name.upper()))
    prefix = _CANONICAL_PREFIXES.get(name)
    if value is None or (prefix and not value.startswith(prefix)):
        raise PoseIdentityError(
            "missing_canonical_pose_identity",
            f"post-QC SDF record requires canonical {name!r} ({prefix or ''}...)",
            details={"field": name},
        )
    return value


def _json_property(
    mol: Any,
    name: str,
    *,
    required: bool,
    default: Any,
) -> Any:
    raw = _first_property(mol, (name, name.upper()))
    if raw is None:
        if required:
            raise PoseIdentityError(
                "missing_canonical_pose_identity",
                f"post-QC SDF record requires {name!r}",
                details={"field": name},
            )
        return default
    try:
        return json.loads(raw)
    except (TypeError, ValueError) as exc:
        raise PoseIdentityError(
            "invalid_canonical_pose_identity",
            f"{name} must be valid JSON",
            details={"field": name},
        ) from exc


def _integer_property(
    mol: Any,
    name: str,
    *,
    required: bool,
    minimum: int,
    default: int,
) -> int:
    raw = _first_property(mol, (name, name.upper()))
    if raw is None:
        if required:
            raise PoseIdentityError(
                "missing_post_qc_lineage",
                f"post-QC SDF record requires {name!r}",
                details={"field": name},
            )
        return default
    try:
        value = int(raw)
    except (TypeError, ValueError) as exc:
        raise PoseIdentityError(
            "invalid_post_qc_lineage",
            f"{name} must be an integer",
            details={"field": name, "value": raw},
        ) from exc
    if value < minimum:
        raise PoseIdentityError(
            "invalid_post_qc_lineage",
            f"{name} must be >= {minimum}",
            details={"field": name, "value": value},
        )
    return value


def _canonical_record_identity(
    mol: Any,
    index: int,
    *,
    require_post_qc_lineage: bool,
) -> dict[str, Any]:
    declared = {
        name: _required_property(mol, name)
        for name in (
            "pose_id",
            "chemical_id",
            "source_identity_id",
            "chemical_identity_id",
        )
    }
    source_identity = _json_property(
        mol, "source_identity", required=True, default=None
    )
    lineage_metadata = _json_property(
        mol, "lineage_metadata", required=False, default={}
    )
    if not isinstance(lineage_metadata, dict):
        raise PoseIdentityError(
            "invalid_canonical_pose_identity",
            "lineage_metadata must be a JSON object",
            details={"field": "lineage_metadata"},
        )
    parent_pose_id = _first_property(
        mol, ("parent_pose_id", "PARENT_POSE_ID")
    )
    coordinate_revision = _integer_property(
        mol,
        "coordinate_revision",
        required=True,
        minimum=0,
        default=0,
    )
    selection_rank = _integer_property(
        mol,
        "selection_rank",
        required=require_post_qc_lineage,
        minimum=1,
        default=index + 1,
    )
    source_pose_index = _integer_property(
        mol,
        "source_pose_index",
        required=require_post_qc_lineage,
        minimum=0,
        default=index,
    )

    computed = build_pose_identity(
        mol,
        source_identity=source_identity,
        parent_pose_id=parent_pose_id,
        coordinate_revision=coordinate_revision,
        lineage_metadata=lineage_metadata,
    )
    chemical = computed["chemical_identity"]
    expected = {
        "pose_id": computed["pose_id"],
        "chemical_id": chemical["chemical_id"],
        "source_identity_id": chemical["source_identity_id"],
        "chemical_identity_id": chemical["identity_id"],
    }
    mismatches = {
        name: {"declared": declared[name], "computed": expected[name]}
        for name in declared
        if declared[name] != expected[name]
    }
    if mismatches:
        raise PoseIdentityError(
            "lineage_identity_mismatch",
            "declared canonical identity does not match the molecule/coordinates",
            details={"record_index": index, "mismatches": mismatches},
        )

    conformer = mol.GetConformer()
    if not conformer.Is3D():
        raise PoseIdentityError(
            "missing_3d_coordinates",
            "post-QC candidate is not marked as a 3D conformer",
            details={"record_index": index},
        )
    return {
        "record_index": index,
        "selection_rank": selection_rank,
        "source_pose_index": source_pose_index,
        "pose_id": computed["pose_id"],
        "parent_pose_id": computed["parent_pose_id"],
        "coordinate_revision": computed["coordinate_revision"],
        "coordinate_signature": computed["coordinates"]["signature"],
        "chemical_identity_id": chemical["identity_id"],
        "chemical_id": chemical["chemical_id"],
        "source_identity_id": chemical["source_identity_id"],
        "source_identity": chemical["source_identity"],
        "canonical_isomeric_smiles": chemical["canonical_isomeric_smiles"],
        "formal_charge": int(chemical["formal_charge"]),
        "heavy_atom_count": int(chemical["heavy_atom_count"]),
        "fragment_count": int(chemical["fragment_count"]),
    }


def _load_post_qc_candidates(
    candidate_sdf: Path,
    *,
    max_candidates: int,
    require_post_qc_lineage: bool,
) -> list[dict[str, Any]]:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise TIESPlanningError(
            "RDKit is required to preflight TIES candidates"
        ) from exc

    supplier = Chem.SDMolSupplier(
        str(candidate_sdf),
        removeHs=False,
        sanitize=False,
        strictParsing=True,
    )
    records: list[dict[str, Any]] = []
    for index, mol in enumerate(supplier):
        if index >= max_candidates:
            break
        if mol is None:
            records.append({
                "record_index": index,
                "status": "failed",
                "reason_code": "sdf_parse_failed",
                "reason": "RDKit could not parse this post-QC SDF record",
            })
            continue
        try:
            identity = _canonical_record_identity(
                mol,
                index,
                require_post_qc_lineage=require_post_qc_lineage,
            )
        except PoseIdentityError as exc:
            records.append({
                "record_index": index,
                "status": "failed",
                "reason_code": exc.code,
                "reason": str(exc),
                "identity_error": exc.as_result()["error"],
            })
            continue
        except Exception as exc:
            records.append({
                "record_index": index,
                "status": "failed",
                "reason_code": "candidate_identity_failed",
                "reason": str(exc),
            })
            continue

        elements = sorted({atom.GetSymbol() for atom in mol.GetAtoms()})
        unsupported = sorted(set(elements) - SUPPORTED_LIGAND_ELEMENTS)
        radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        record = {
            **identity,
            "status": "passed",
            "reason_code": None,
            "reason": None,
            "preflight": {
                "elements": elements,
                "radical_electron_count": int(radicals),
                "post_qc_lineage_present": bool(require_post_qc_lineage),
            },
        }
        if unsupported:
            record.update({
                "status": "not_evaluable",
                "reason_code": "unsupported_ligand_elements",
                "reason": (
                    "Generic TIES/GAFF2 planning route does not accept elements: "
                    + ", ".join(unsupported)
                ),
            })
        elif identity["fragment_count"] != 1:
            record.update({
                "status": "not_evaluable",
                "reason_code": "disconnected_ligand",
                "reason": "TIES planning requires one connected ligand entity",
            })
        elif radicals:
            record.update({
                "status": "not_evaluable",
                "reason_code": "radical_ligand",
                "reason": "Generic TIES/GAFF2 planning does not accept radicals",
            })
        records.append(record)

    seen_poses: set[str] = set()
    seen_chemicals: set[str] = set()
    for record in records:
        if record["status"] != "passed":
            continue
        pose_id = record["pose_id"]
        chemical_id = record["chemical_id"]
        if pose_id in seen_poses:
            record.update({
                "status": "failed",
                "reason_code": "duplicate_pose_identity",
                "reason": "post-QC candidate SDF contains a duplicate canonical pose_id",
            })
            continue
        seen_poses.add(pose_id)
        if chemical_id in seen_chemicals:
            record.update({
                "status": "not_evaluable",
                "reason_code": "duplicate_chemical_candidate",
                "reason": (
                    "RBFE network contains this chemical identity more than once; "
                    "keep one post-QC pose per chemical"
                ),
            })
            continue
        seen_chemicals.add(chemical_id)
    return records


def _command_record(result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "command": list(result.get("command") or []),
        "environment": result.get("environment"),
        "returncode": result.get("returncode"),
        "elapsed_seconds": result.get("elapsed_seconds"),
        "stdout_tail": str(result.get("stdout") or "")[-4000:],
        "stderr_tail": str(result.get("stderr") or "")[-4000:],
    }


def inspect_ties_receptor(receptor_pdb: str | Path) -> dict[str, Any]:
    """Validate the protein receptor and retain generic Amber risk evidence."""
    path = Path(receptor_pdb)
    atom_record_count = 0
    protein_atom_count = 0
    protein_residue_names: set[str] = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_record_count += 1
        resname = line[17:20].strip().upper()
        if resname in PROTEIN_RESNAMES:
            protein_atom_count += 1
            protein_residue_names.add(resname)

    risks = inspect_receptor_risks(path)
    blockers: list[dict[str, Any]] = []
    if protein_atom_count == 0:
        blockers.append({
            "code": "receptor_missing_protein_atoms",
            "reason": (
                "Receptor PDB must contain ATOM records for recognized protein "
                "residues"
            ),
            "evidence": {
                "atom_record_count": atom_record_count,
                "protein_atom_count": protein_atom_count,
            },
        })
    if risks["metal_atoms"]:
        blockers.append({
            "code": "metal_site_requires_preparameterized_protocol",
            "reason": (
                "Planning-only generic TIES route cannot infer a reviewed "
                "metal-center model"
            ),
            "evidence": {"metal_atoms": risks["metal_atoms"]},
        })
    if risks["unsupported_hetero_residues"]:
        blockers.append({
            "code": "unsupported_receptor_hetero_residues",
            "reason": (
                "Generic TIES planning does not automatically parameterize "
                "non-water receptor HETATM residues"
            ),
            "evidence": {
                "residue_names": risks["unsupported_hetero_residues"],
            },
        })
    return {
        "status": "passed" if not blockers else "not_evaluable",
        "protein_atom_records_present": protein_atom_count > 0,
        "atom_record_count": atom_record_count,
        "protein_atom_count": protein_atom_count,
        "protein_residue_names": sorted(protein_residue_names),
        "risks": risks,
        "blockers": blockers,
    }


def _runtime_probe_code(simulation_platform: str) -> str:
    requested = json.dumps(simulation_platform)
    return f"""
import importlib
import importlib.metadata as md
import json
from pathlib import Path
import shutil
import sys

conda_packages = {{}}
conda_meta_dir = Path(sys.prefix) / "conda-meta"
if conda_meta_dir.is_dir():
    for metadata_file in sorted(conda_meta_dir.glob("*.json")):
        try:
            metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
        except Exception:
            continue
        package_name = metadata.get("name")
        if package_name in {{"ties", "ties_md", "openmm", "ambertools"}}:
            conda_packages[package_name] = {{
                "version": metadata.get("version") or "unknown",
                "build": metadata.get("build") or "unknown",
            }}

modules = {{}}
loaded_modules = {{}}
missing = []
for label, name, distributions in (
    ("builder", "ties", ("ties", "ties-lib")),
    ("runner", "TIES_MD", ("TIES_MD", "ties-md")),
    ("simulation_engine", "openmm", ("OpenMM", "openmm")),
):
    try:
        module = importlib.import_module(name)
        loaded_modules[label] = module
        version = getattr(module, "__version__", None)
        if not version:
            for distribution in distributions:
                try:
                    version = md.version(distribution)
                    break
                except md.PackageNotFoundError:
                    pass
        modules[label] = {{"module": name, "version": version or "unknown"}}
    except Exception as exc:
        modules[label] = {{"module": name, "error": str(exc)}}
        missing.append(name)

api_contracts = {{}}
builder = loaded_modules.get("builder")
pair = getattr(builder, "Pair", None) if builder is not None else None
api_contracts["ties.Pair"] = callable(pair)
api_contracts["ties.Config"] = (
    callable(getattr(builder, "Config", None))
    if builder is not None
    else False
)
api_contracts["ties.Pair.superimpose"] = (
    callable(getattr(pair, "superimpose", None))
    if pair is not None
    else False
)

try:
    runner_api = importlib.import_module("TIES_MD.TIES")
    api_contracts["TIES_MD.TIES.TIES"] = callable(
        getattr(runner_api, "TIES", None)
    )
except Exception:
    api_contracts["TIES_MD.TIES.TIES"] = False

try:
    importlib.import_module("TIES_MD.ties_analysis")
    api_contracts["TIES_MD.ties_analysis"] = True
except Exception:
    api_contracts["TIES_MD.ties_analysis"] = False

missing.extend(
    "api:" + name
    for name, available in api_contracts.items()
    if not available
)

command_names = {list(REQUIRED_COMMANDS)!r}
commands = {{name: shutil.which(name) for name in command_names}}
missing.extend(name for name, path in commands.items() if not path)

requested_platform = {requested}
platforms = {{}}
selected_platform = None
try:
    import openmm
    from openmm import unit

    registered = [
        openmm.Platform.getPlatform(index).getName()
        for index in range(openmm.Platform.getNumPlatforms())
    ]
    candidates = (
        {list(OPENMM_PLATFORMS)!r}
        if requested_platform == "auto"
        else [requested_platform]
    )
    for name in candidates:
        if name not in registered:
            platforms[name] = {{
                "registered": False,
                "functional": False,
                "error": "platform is not registered",
            }}
            continue
        context = None
        integrator = None
        try:
            system = openmm.System()
            system.addParticle(39.9)
            force = openmm.CustomExternalForce("x*x+y*y+z*z")
            force.addParticle(0, [])
            system.addForce(force)
            integrator = openmm.VerletIntegrator(0.001)
            platform = openmm.Platform.getPlatformByName(name)
            context = openmm.Context(system, integrator, platform)
            context.setPositions([[0.0, 0.0, 0.0]])
            energy = context.getState(getEnergy=True).getPotentialEnergy()
            energy.value_in_unit(unit.kilojoule_per_mole)
            platforms[name] = {{
                "registered": True,
                "functional": True,
            }}
            if selected_platform is None:
                selected_platform = name
        except Exception as exc:
            platforms[name] = {{
                "registered": True,
                "functional": False,
                "error": str(exc),
            }}
        finally:
            del context
            del integrator
except Exception as exc:
    platforms["probe"] = {{
        "registered": False,
        "functional": False,
        "error": str(exc),
    }}

if selected_platform is None:
    missing.append("openmm_platform:" + requested_platform)

print(json.dumps({{
    "python": sys.version.split()[0],
    "executable": sys.executable,
    "conda_packages": conda_packages,
    "modules": modules,
    "api_contracts": api_contracts,
    "commands": commands,
    "openmm_platforms": platforms,
    "requested_simulation_platform": requested_platform,
    "selected_simulation_platform": selected_platform,
    "missing": sorted(set(missing)),
}}, sort_keys=True))
sys.exit(0 if not missing else 2)
"""


def probe_ties_runtime(
    runtime: Mapping[str, Any],
    *,
    timeout_seconds: int = 30,
    simulation_platform: str = "auto",
) -> dict[str, Any]:
    """Probe the complete TIES, OpenMM, and AmberTools planning runtime."""
    commands = dict(runtime.get("commands") or {})
    try:
        result = run_tool(
            dict(runtime),
            ["-c", _runtime_probe_code(simulation_platform)],
            command=commands.get("python", "python"),
            timeout=max(1, int(timeout_seconds)),
            check=False,
        )
    except (OSError, ToolRuntimeError) as exc:
        failed = getattr(exc, "result", None)
        return {
            "status": "not_evaluable",
            "available": False,
            "reason_code": "ties_runtime_unavailable",
            "reason": str(exc),
            "requirements": {
                "python_modules": list(REQUIRED_PYTHON_MODULES),
                "api_contracts": list(REQUIRED_API_CONTRACTS),
                "commands": list(REQUIRED_COMMANDS),
            },
            "discovered": {},
            "probe": _command_record(failed) if failed else {
                "command": ["python", "-c", "<ties-runtime-probe>"],
                "environment": runtime.get("environment"),
            },
        }

    payload: dict[str, Any] = {}
    for line in reversed(str(result.get("stdout") or "").splitlines()):
        try:
            value = json.loads(line)
        except (TypeError, ValueError):
            continue
        if isinstance(value, dict) and "modules" in value:
            payload = value
            break
    modules = payload.get("modules")
    api_contracts = payload.get("api_contracts")
    discovered_commands = payload.get("commands")
    module_names = {
        entry.get("module")
        for entry in modules.values()
        if isinstance(entry, dict) and not entry.get("error")
    } if isinstance(modules, dict) else set()
    commands_available = (
        isinstance(discovered_commands, dict)
        and all(discovered_commands.get(name) for name in REQUIRED_COMMANDS)
    )
    api_contracts_available = (
        isinstance(api_contracts, dict)
        and all(api_contracts.get(name) for name in REQUIRED_API_CONTRACTS)
    )
    complete_evidence = (
        set(REQUIRED_PYTHON_MODULES).issubset(module_names)
        and api_contracts_available
        and commands_available
    )
    available = (
        bool(result.get("success"))
        and complete_evidence
        and not payload.get("missing")
    )
    return {
        "status": "passed" if available else "not_evaluable",
        "available": available,
        "reason_code": None if available else "ties_runtime_unavailable",
        "reason": (
            None
            if available
            else (
                "TIES builder, TIES_MD/OpenMM runtime, or required "
                "TIES/AmberTools CLI is unavailable"
            )
        ),
        "requirements": {
            "python_modules": list(REQUIRED_PYTHON_MODULES),
            "api_contracts": list(REQUIRED_API_CONTRACTS),
            "commands": list(REQUIRED_COMMANDS),
        },
        "discovered": payload,
        "probe": _command_record(result),
    }


def _transformation_id(reference: Mapping[str, Any], target: Mapping[str, Any]) -> str:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "reference_pose_id": reference["pose_id"],
        "reference_chemical_id": reference["chemical_id"],
        "target_pose_id": target["pose_id"],
        "target_chemical_id": target["chemical_id"],
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return f"ties-transformation:{digest}"


def _plan_transformations(
    candidates: list[dict[str, Any]],
    *,
    reference_pose_id: str | None,
    runtime_available: bool,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    eligible = sorted(
        (item for item in candidates if item["status"] == "passed"),
        key=lambda item: (item["selection_rank"], item["record_index"]),
    )
    if len(eligible) < 2:
        return None, []
    if reference_pose_id is None:
        reference = eligible[0]
    else:
        reference = next(
            (item for item in eligible if item["pose_id"] == reference_pose_id),
            None,
        )
        if reference is None:
            raise TIESPlanningError(
                "reference_pose_id must identify a passed post-QC candidate"
            )

    transformations: list[dict[str, Any]] = []
    for target in eligible:
        if target["pose_id"] == reference["pose_id"]:
            continue
        status = "passed"
        reason_code = None
        reason = None
        charge_delta = target["formal_charge"] - reference["formal_charge"]
        if charge_delta:
            status = "not_evaluable"
            reason_code = "net_charge_change_unsupported"
            reason = (
                "This TIES backend does not support charge-changing "
                "transformations"
            )
        if status == "passed" and not runtime_available:
            status = "not_evaluable"
            reason_code = "ties_runtime_unavailable"
            reason = "TIES/TIES_MD runtime preflight did not pass"

        transformations.append({
            "transformation_id": _transformation_id(reference, target),
            "status": status,
            "reason_code": reason_code,
            "reason": reason,
            "reference": {
                "pose_id": reference["pose_id"],
                "chemical_id": reference["chemical_id"],
                "formal_charge": reference["formal_charge"],
            },
            "target": {
                "pose_id": target["pose_id"],
                "chemical_id": target["chemical_id"],
                "formal_charge": target["formal_charge"],
            },
            "charge_delta": charge_delta,
            "review_flags": [],
            "planned_stages": [
                {
                    "stage": "hybrid_mapping_and_parameterization",
                    "implementation": "ties.Pair.superimpose/prepare_inputs",
                    "status": "pending",
                },
                {
                    "stage": "complex_and_ligand_ensemble_simulations",
                    "implementation": "TIES_MD",
                    "status": "pending",
                },
                {
                    "stage": "thermodynamic_cycle_analysis",
                    "implementation": "ties_ana",
                    "status": "pending",
                },
            ],
        })
    return reference, transformations


def _status_counts(entries: list[Mapping[str, Any]]) -> dict[str, int]:
    counts = {
        status: sum(item.get("status") == status for item in entries)
        for status in ("passed", "failed", "not_evaluable")
    }
    counts["requested"] = len(entries)
    return counts


def plan_ties_rbfe(
    receptor_pdb: str | Path,
    candidate_sdf: str | Path,
    output_dir: str | Path,
    runtime: Mapping[str, Any],
    *,
    execution_mode: str = "preflight",
    network_strategy: str = "star",
    reference_pose_id: str | None = None,
    max_candidates: int = 10,
    allow_charge_changes: bool = False,
    is_covalent: bool = False,
    has_metal: bool = False,
    require_post_qc_lineage: bool = True,
    runtime_probe_timeout_seconds: int = 30,
    simulation_platform: str = "auto",
) -> dict[str, Any]:
    """Write a strict preflight or star-network plan for a later TIES run."""
    if execution_mode not in {"preflight", "plan"}:
        raise TIESPlanningError("execution_mode must be 'preflight' or 'plan'")
    if allow_charge_changes:
        raise TIESPlanningError(
            "charge-changing transformations are unsupported by this TIES backend"
        )
    if network_strategy != "star":
        raise TIESPlanningError("the planning-only backend supports star networks")
    if isinstance(max_candidates, bool) or int(max_candidates) < 2:
        raise TIESPlanningError("max_candidates must be an integer >= 2")
    if (
        isinstance(runtime_probe_timeout_seconds, bool)
        or int(runtime_probe_timeout_seconds) < 1
    ):
        raise TIESPlanningError("runtime_probe_timeout_seconds must be >= 1")
    if simulation_platform not in {"auto", *OPENMM_PLATFORMS}:
        raise TIESPlanningError(
            "simulation_platform must be auto, CUDA, OpenCL, or CPU"
        )

    receptor_path = Path(receptor_pdb).resolve()
    candidate_path = Path(candidate_sdf).resolve()
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "ties_rbfe_plan.json"
    missing = [
        str(path)
        for path in (receptor_path, candidate_path)
        if not path.is_file()
    ]
    if missing:
        raise TIESPlanningError(
            "TIES planning inputs not found: " + ", ".join(missing)
        )

    candidates = _load_post_qc_candidates(
        candidate_path,
        max_candidates=int(max_candidates),
        require_post_qc_lineage=bool(require_post_qc_lineage),
    )
    receptor_validation = inspect_ties_receptor(receptor_path)
    global_blockers: list[dict[str, Any]] = list(
        receptor_validation["blockers"]
    )
    if is_covalent:
        global_blockers.append({
            "code": "covalent_complex_requires_specialized_protocol",
            "reason": (
                "Planning-only generic TIES route has no reviewed covalent "
                "protein-ligand parameterization"
            ),
        })
    if has_metal and not receptor_validation["risks"]["metal_atoms"]:
        global_blockers.append({
            "code": "metal_site_requires_preparameterized_protocol",
            "reason": (
                "Planning-only generic TIES route cannot infer a reviewed "
                "metal-center model"
            ),
        })
    eligible_charges = {
        candidate["formal_charge"]
        for candidate in candidates
        if candidate["status"] == "passed"
    }
    if len(eligible_charges) > 1:
        global_blockers.append({
            "code": "mixed_formal_charge",
            "reason": (
                "This TIES backend does not support a cohort containing "
                "charge-changing transformations"
            ),
            "evidence": {"formal_charges": sorted(eligible_charges)},
        })
    if global_blockers:
        reason = "; ".join(item["reason"] for item in global_blockers)
        for candidate in candidates:
            if candidate["status"] == "passed":
                candidate.update({
                    "status": "not_evaluable",
                    "reason_code": "global_applicability_blocker",
                    "reason": reason,
                })

    runtime_probe = probe_ties_runtime(
        runtime,
        timeout_seconds=int(runtime_probe_timeout_seconds),
        simulation_platform=simulation_platform,
    )
    reference: dict[str, Any] | None = None
    transformations: list[dict[str, Any]] = []
    if execution_mode == "plan":
        reference, transformations = _plan_transformations(
            candidates,
            reference_pose_id=reference_pose_id,
            runtime_available=bool(runtime_probe["available"]),
        )

    candidate_counts = _status_counts(candidates)
    transformation_counts = _status_counts(transformations)
    if candidate_counts["failed"]:
        status = "failed"
        reason_code = "candidate_contract_failed"
    elif global_blockers:
        status = "not_evaluable"
        reason_code = global_blockers[0]["code"]
    elif candidate_counts["passed"] < 2:
        status = "not_evaluable"
        reason_code = "insufficient_distinct_eligible_ligands"
    elif not runtime_probe["available"]:
        status = "not_evaluable"
        reason_code = "ties_runtime_unavailable"
    elif execution_mode == "preflight":
        status = "passed"
        reason_code = None
    elif transformation_counts["passed"] == 0:
        status = "not_evaluable"
        reason_code = "no_evaluable_transformations"
    elif (
        transformation_counts["failed"]
        or transformation_counts["not_evaluable"]
    ):
        status = "failed"
        reason_code = "incomplete_transformation_plan"
    else:
        status = "passed"
        reason_code = None

    implementation_path = Path(__file__).resolve()
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "success": True,
        "status": status,
        "reason_code": reason_code,
        "execution_mode": execution_mode,
        "inputs": {
            "receptor_pdb": str(receptor_path),
            "receptor_sha256": _sha256(receptor_path),
            "post_qc_candidate_sdf": str(candidate_path),
            "post_qc_candidate_sdf_sha256": _sha256(candidate_path),
        },
        "protocol_plan": {
            "method": "TIES_RBFE",
            "network_strategy": network_strategy,
            "reference_pose_id": (
                reference["pose_id"] if reference is not None else reference_pose_id
            ),
            "charge_changes_supported": False,
            "thermodynamic_cycle_legs": ["complex", "ligand"],
            "requested_simulation_platform": simulation_platform,
            "selected_simulation_platform": (
                (runtime_probe.get("discovered") or {}).get(
                    "selected_simulation_platform"
                )
            ),
            "max_candidates": int(max_candidates),
            "global_blockers": global_blockers,
        },
        "receptor_validation": receptor_validation,
        "tool_provenance": {
            "environment": runtime.get("environment"),
            "registered_commands": dict(runtime.get("commands") or {}),
            "runtime_probe": runtime_probe,
            "wrapper": {
                "schema_version": SCHEMA_VERSION,
                "path": str(implementation_path),
                "sha256": _sha256(implementation_path),
            },
        },
        "candidates": candidates,
        "candidate_counts": candidate_counts,
        "transformations": transformations,
        "transformation_counts": transformation_counts,
        "counts": {
            "candidates": candidate_counts,
            "transformations": transformation_counts,
        },
        "free_energy_result": {
            "status": "not_computed",
            "protocol_executed": False,
            "reason_code": "planning_only_backend",
            "reason": (
                "TIES hybrid setup, both thermodynamic legs, ensemble "
                "simulations, and analysis have not been executed"
            ),
        },
        "manifest_file": str(manifest_path),
    }
    _atomic_json(manifest_path, manifest)
    return manifest


__all__ = [
    "SCHEMA_VERSION",
    "REQUIRED_API_CONTRACTS",
    "REQUIRED_COMMANDS",
    "REQUIRED_PYTHON_MODULES",
    "OPENMM_PLATFORMS",
    "SUPPORTED_LIGAND_ELEMENTS",
    "TIESPlanningError",
    "inspect_ties_receptor",
    "plan_ties_rbfe",
    "probe_ties_runtime",
]

"""Strict OpenFE RBFE preflight and ligand-network planning.

This backend deliberately stops before system parameterization, protocol
construction, and simulation.  ``preflight`` validates the OpenFE runtime,
receptor, canonical pose identities, and ligand applicability.  ``plan`` adds
an OpenFE ligand mapping network.  Neither mode produces or accepts a binding
free-energy estimate.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import re
import sys
from typing import Any, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
BIO_RESEARCH_DIR = SCRIPT_DIR.parents[1]
SHARED_DIR = BIO_RESEARCH_DIR / "shared"
if str(SHARED_DIR) not in sys.path:
    sys.path.insert(0, str(SHARED_DIR))

from pose_identity import (  # noqa: E402
    PoseIdentityError,
    build_chemical_identity,
)
from tool_runtime import run_tool  # noqa: E402

# Reuse the evaluator's canonical SDF-property parser.  OpenFE applies its own
# chemistry policy below; Amber-specific applicability decisions are not used.
from candidate_topology import (  # noqa: E402
    _first_property,
    _record_identity,
    inspect_receptor_risks,
)


SCHEMA_VERSION = "openfe-rbfe-preflight/v1"
WORKER_RESULT_PREFIX = "OPENFE_RBFE_RESULT="
SUPPORTED_ELEMENTS = frozenset(
    {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
)
EXECUTION_MODES = frozenset({"preflight", "plan"})
NETWORK_METHODS = frozenset(
    {"minimal_spanning", "minimal_redundant", "radial"}
)
MAPPERS = frozenset({"KartografAtomMapper", "LomapAtomMapper"})


class OpenFERBFEError(RuntimeError):
    """Raised for malformed backend invocations, not scientific inapplicability."""


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


def _bounded_invocation(result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "command": list(result.get("command") or []),
        "environment": result.get("environment"),
        "returncode": result.get("returncode"),
        "elapsed_seconds": result.get("elapsed_seconds"),
        "stdout_tail": str(result.get("stdout") or "")[-4000:],
        "stderr_tail": str(result.get("stderr") or "")[-4000:],
    }


def _require_file(value: str | Path, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise OpenFERBFEError(f"{label} does not exist: {path}")
    if path.stat().st_size == 0:
        raise OpenFERBFEError(f"{label} is empty: {path}")
    return path


def _not_evaluable_record(
    index: int,
    *,
    reason_code: str,
    reason: str,
    pose_id: str | None = None,
) -> dict[str, Any]:
    return {
        "record_index": index,
        "pose_id": pose_id,
        "chemical_id": None,
        "source_identity_id": None,
        "formal_charge": None,
        "status": "not_evaluable",
        "reason_code": reason_code,
        "reason": reason,
        "preflight": {},
        "_mol": None,
    }


def _load_openfe_candidates(
    candidate_sdf: Path,
    *,
    max_candidates: int,
    require_explicit_hydrogens: bool,
    allow_undefined_stereochemistry: bool,
) -> list[dict[str, Any]]:
    """Load candidates and require canonical post-QC identity properties."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise OpenFERBFEError(
            "RDKit is required to preflight OpenFE candidates"
        ) from exc

    supplier = Chem.SDMolSupplier(
        str(candidate_sdf),
        removeHs=False,
        sanitize=False,
        strictParsing=True,
    )
    source_hash = _sha256(candidate_sdf)
    records: list[dict[str, Any]] = []
    for index, mol in enumerate(supplier):
        if index >= max_candidates:
            break
        if mol is None:
            records.append(
                _not_evaluable_record(
                    index,
                    reason_code="sdf_parse_failed",
                    reason="RDKit could not parse this SDF record",
                )
            )
            continue

        declared_pose_id = _first_property(mol, ("pose_id", "POSE_ID"))
        try:
            Chem.SanitizeMol(mol)
            identity = _record_identity(
                mol,
                index,
                source_artifact_sha256=source_hash,
            )
            missing_canonical = [
                field
                for field in (
                    "input_declared_pose_id",
                    "input_declared_chemical_id",
                    "input_declared_source_identity_id",
                )
                if not str(identity.get(field) or "").strip()
            ]
            if missing_canonical:
                raise PoseIdentityError(
                    "missing_canonical_pose_identity",
                    "OpenFE requires post-QC canonical pose, chemical, and "
                    "source identity properties",
                    details={"missing_fields": missing_canonical},
                )
            prefixes = (
                ("input_declared_pose_id", "pose:"),
                ("input_declared_chemical_id", "chem:"),
                ("input_declared_source_identity_id", "source:"),
            )
            invalid_prefixes = [
                field
                for field, prefix in prefixes
                if not str(identity[field]).startswith(prefix)
            ]
            if invalid_prefixes:
                raise PoseIdentityError(
                    "noncanonical_pose_identity",
                    "OpenFE input identity properties are not canonical",
                    details={"invalid_fields": invalid_prefixes},
                )
            chemical_identity = build_chemical_identity(
                mol,
                source_identity=identity["source_identity"],
            )
        except PoseIdentityError as exc:
            records.append(
                _not_evaluable_record(
                    index,
                    reason_code=exc.code,
                    reason=str(exc),
                    pose_id=declared_pose_id,
                )
            )
            continue
        except Exception as exc:
            records.append(
                _not_evaluable_record(
                    index,
                    reason_code="chemical_sanitization_failed",
                    reason=str(exc),
                    pose_id=declared_pose_id,
                )
            )
            continue

        elements = sorted({atom.GetSymbol() for atom in mol.GetAtoms()})
        unsupported = sorted(set(elements) - SUPPORTED_ELEMENTS)
        fragment_count = len(Chem.GetMolFrags(mol))
        radical_count = sum(
            atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()
        )
        explicit_hydrogen_count = sum(
            atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()
        )
        has_3d = bool(mol.GetNumConformers()) and bool(
            mol.GetConformer().Is3D()
        )
        stereo_unassigned = bool(
            chemical_identity["stereochemistry"]["has_unassigned"]
        )
        entry: dict[str, Any] = {
            **identity,
            "record_index": index,
            "status": "passed",
            "reason_code": None,
            "reason": None,
            "preflight": {
                "elements": elements,
                "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
                "fragment_count": fragment_count,
                "radical_electron_count": int(radical_count),
                "explicit_hydrogen_count": explicit_hydrogen_count,
                "has_3d_coordinates": has_3d,
                "has_unassigned_stereochemistry": stereo_unassigned,
            },
            "_mol": mol,
        }
        if unsupported:
            entry.update(
                status="not_evaluable",
                reason_code="unsupported_ligand_elements",
                reason=(
                    "Default OpenFE/OpenFF route does not cover elements: "
                    + ", ".join(unsupported)
                ),
            )
        elif fragment_count != 1:
            entry.update(
                status="not_evaluable",
                reason_code="disconnected_ligand",
                reason="Candidate must contain one connected chemical entity",
            )
        elif radical_count:
            entry.update(
                status="not_evaluable",
                reason_code="radical_ligand",
                reason="Default OpenFE RBFE route does not accept radicals",
            )
        elif not has_3d:
            entry.update(
                status="not_evaluable",
                reason_code="missing_3d_coordinates",
                reason="Candidate has no declared 3D conformer",
            )
        elif require_explicit_hydrogens and not explicit_hydrogen_count:
            entry.update(
                status="not_evaluable",
                reason_code="missing_explicit_hydrogens",
                reason=(
                    "Default OpenFE RBFE input requires explicit hydrogens; "
                    "prepare them upstream without overwriting the pose identity"
                ),
            )
        elif stereo_unassigned and not allow_undefined_stereochemistry:
            entry.update(
                status="not_evaluable",
                reason_code="undefined_stereochemistry",
                reason=(
                    "Candidate contains undefined stereochemistry; enumerate "
                    "a defined microstate/stereoisomer before RBFE"
                ),
            )
        records.append(entry)
    return records


def _public_candidate(entry: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in entry.items()
        if key != "_mol" and not key.startswith("input_declared_")
    }


def _counts(candidates: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    result = {
        status: sum(entry.get("status") == status for entry in candidates)
        for status in ("passed", "failed", "not_evaluable")
    }
    result["requested"] = len(candidates)
    return result


def _overall_status(candidates: Sequence[Mapping[str, Any]]) -> str:
    statuses = [entry.get("status") for entry in candidates]
    if not statuses or all(status == "not_evaluable" for status in statuses):
        return "not_evaluable"
    if any(status == "failed" for status in statuses):
        return "failed"
    if any(status == "not_evaluable" for status in statuses):
        return "failed"
    return "passed"


def _global_block(
    candidates: list[dict[str, Any]],
    blockers: Sequence[Mapping[str, str]],
) -> None:
    if not blockers:
        return
    codes = ",".join(str(blocker["code"]) for blocker in blockers)
    reason = "; ".join(str(blocker["reason"]) for blocker in blockers)
    for candidate in candidates:
        if candidate.get("status") == "passed":
            candidate.update(
                status="not_evaluable",
                reason_code=f"global_blocker:{codes}",
                reason=reason,
            )


def _parse_worker_result(execution: Mapping[str, Any]) -> dict[str, Any] | None:
    for line in reversed(str(execution.get("stdout") or "").splitlines()):
        if line.startswith(WORKER_RESULT_PREFIX):
            try:
                result = json.loads(line[len(WORKER_RESULT_PREFIX) :])
            except json.JSONDecodeError:
                return None
            return result if isinstance(result, dict) else None
    return None


def _invoke_worker(
    runtime: Mapping[str, Any],
    payload: Mapping[str, Any],
    payload_file: Path,
    *,
    timeout_seconds: int,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    _atomic_json(payload_file, payload)
    execution = run_tool(
        dict(runtime),
        [str(Path(__file__).resolve()), "--worker", str(payload_file)],
        command="python",
        cwd=payload_file.parent,
        timeout=timeout_seconds,
        check=False,
    )
    return _parse_worker_result(execution), execution


def _worker_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for distribution in ("openfe", "gufe", "rdkit", "openff-toolkit"):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def _worker_load_components(
    payload: Mapping[str, Any],
) -> tuple[list[Any], list[dict[str, Any]]]:
    import openfe
    from rdkit import Chem

    wanted = {
        int(entry["record_index"]): dict(entry)
        for entry in payload["eligible_candidates"]
    }
    supplier = Chem.SDMolSupplier(
        str(payload["candidate_sdf"]),
        removeHs=False,
        sanitize=True,
        strictParsing=True,
    )
    components: list[Any] = []
    outcomes: list[dict[str, Any]] = []
    for index, mol in enumerate(supplier):
        if index not in wanted:
            continue
        expected = wanted[index]
        if mol is None:
            outcomes.append(
                {
                    "pose_id": expected["pose_id"],
                    "status": "failed",
                    "reason_code": "openfe_sdf_parse_failed",
                    "reason": "OpenFE-environment RDKit could not parse record",
                }
            )
            continue
        actual_pose_id = (
            mol.GetProp("pose_id").strip() if mol.HasProp("pose_id") else ""
        )
        actual_chemical_id = (
            mol.GetProp("chemical_id").strip()
            if mol.HasProp("chemical_id")
            else ""
        )
        if (
            actual_pose_id != expected["pose_id"]
            or actual_chemical_id != expected["chemical_id"]
        ):
            outcomes.append(
                {
                    "pose_id": expected["pose_id"],
                    "status": "failed",
                    "reason_code": "worker_identity_mismatch",
                    "reason": (
                        "OpenFE worker observed different pose/chemical "
                        "identity properties"
                    ),
                }
            )
            continue
        component_name = f"RBFE_{len(components):04d}"
        try:
            component = openfe.SmallMoleculeComponent.from_rdkit(
                mol,
                name=component_name,
            )
        except Exception as exc:
            outcomes.append(
                {
                    "pose_id": expected["pose_id"],
                    "status": "failed",
                    "reason_code": "openfe_component_failed",
                    "reason": str(exc),
                }
            )
            continue
        components.append(component)
        outcomes.append(
            {
                "pose_id": expected["pose_id"],
                "chemical_id": expected["chemical_id"],
                "formal_charge": expected["formal_charge"],
                "component_name": component_name,
                "status": "passed",
                "reason_code": None,
                "reason": None,
            }
        )
    missing = sorted(set(wanted) - {index for index in range(len(supplier))})
    for index in missing:
        expected = wanted[index]
        outcomes.append(
            {
                "pose_id": expected["pose_id"],
                "status": "failed",
                "reason_code": "openfe_record_missing",
                "reason": f"SDF record {index} was not present",
            }
        )
    return components, outcomes


def _worker_validate_receptor(receptor_file: str) -> dict[str, Any]:
    import openfe

    path = Path(receptor_file)
    if path.suffix.casefold() in {".cif", ".mmcif", ".pdbx"}:
        protein = openfe.ProteinComponent.from_pdbx_file(
            str(path),
            name="receptor",
        )
    else:
        protein = openfe.ProteinComponent.from_pdb_file(
            path,
            name="receptor",
        )
    molecule = protein.to_rdkit()
    return {
        "status": "passed",
        "atom_count": int(molecule.GetNumAtoms()),
        "component_name": protein.name,
    }


def _worker_mapper(name: str, settings: Mapping[str, Any]) -> Any:
    import openfe

    if name == "KartografAtomMapper":
        return openfe.KartografAtomMapper(**dict(settings))
    if name == "LomapAtomMapper":
        return openfe.LomapAtomMapper(**dict(settings))
    raise ValueError(f"unsupported mapper: {name}")


def _worker_plan_network(
    payload: Mapping[str, Any],
    components: Sequence[Any],
    outcomes: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    import openfe

    planner = openfe.ligand_network_planning
    mapper = _worker_mapper(payload["mapper"], payload["mapper_settings"])
    common = {
        "ligands": components,
        "mappers": mapper,
        "scorer": openfe.lomap_scorers.default_lomap_score,
        "progress": False,
        "n_processes": int(payload["n_processes"]),
    }
    method = payload["network_method"]
    if method == "minimal_spanning":
        network = planner.generate_minimal_spanning_network(**common)
    elif method == "minimal_redundant":
        network = planner.generate_minimal_redundant_network(
            **common,
            mst_num=int(payload["mst_num"]),
        )
    elif method == "radial":
        central_pose_id = payload.get("central_pose_id")
        passed = [entry for entry in outcomes if entry["status"] == "passed"]
        component_by_pose_id = {
            entry["pose_id"]: component
            for component, entry in zip(components, passed)
        }
        if central_pose_id is None:
            central_ligand = components[0]
        else:
            central_ligand = component_by_pose_id.get(central_pose_id)
            if central_ligand is None:
                raise ValueError(
                    "central_pose_id does not match an OpenFE-eligible "
                    f"component: {central_pose_id}"
                )
        network = planner.generate_radial_network(
            **common,
            central_ligand=central_ligand,
        )
    else:
        raise ValueError(f"unsupported network method: {method}")

    component_to_outcome = {
        component.name: entry
        for component, entry in zip(
            components,
            (item for item in outcomes if item["status"] == "passed"),
        )
    }
    nodes = [
        {
            "component_name": component.name,
            "pose_id": component_to_outcome[component.name]["pose_id"],
            "chemical_id": component_to_outcome[component.name]["chemical_id"],
            "formal_charge": component_to_outcome[component.name][
                "formal_charge"
            ],
        }
        for component in components
    ]
    edges: list[dict[str, Any]] = []
    for mapping in network.edges:
        molecule_a = mapping.componentA.to_rdkit()
        molecule_b = mapping.componentB.to_rdkit()
        mapped_heavy = sum(
            molecule_a.GetAtomWithIdx(index_a).GetAtomicNum() > 1
            and molecule_b.GetAtomWithIdx(index_b).GetAtomicNum() > 1
            for index_a, index_b in mapping.componentA_to_componentB.items()
        )
        heavy_denominator = min(
            molecule_a.GetNumHeavyAtoms(),
            molecule_b.GetNumHeavyAtoms(),
        )
        score = float(openfe.lomap_scorers.default_lomap_score(mapping))
        source = component_to_outcome[mapping.componentA.name]
        target = component_to_outcome[mapping.componentB.name]
        edges.append(
            {
                "source_pose_id": source["pose_id"],
                "target_pose_id": target["pose_id"],
                "source_chemical_id": source["chemical_id"],
                "target_chemical_id": target["chemical_id"],
                "source_formal_charge": source["formal_charge"],
                "target_formal_charge": target["formal_charge"],
                "mapped_atom_count": len(
                    mapping.componentA_to_componentB
                ),
                "mapped_heavy_atom_count": mapped_heavy,
                "mapped_heavy_fraction": (
                    mapped_heavy / heavy_denominator
                    if heavy_denominator
                    else 0.0
                ),
                "lomap_score": score,
            }
        )
    edges.sort(key=lambda entry: (entry["source_pose_id"], entry["target_pose_id"]))
    graphml_path = Path(payload["graphml_file"])
    graphml_path.write_text(network.to_graphml(), encoding="utf-8")
    return {
        "planning_level": "ligand_mapping_network",
        "simulation_ready": False,
        "node_count": len(nodes),
        "edge_count": len(edges),
        "nodes": nodes,
        "edges": edges,
        "graphml_file": str(graphml_path),
        "graphml_sha256": _sha256(graphml_path),
    }


def _worker(payload: Mapping[str, Any]) -> dict[str, Any]:
    import openfe

    versions = _worker_versions()
    receptor = _worker_validate_receptor(str(payload["receptor_file"]))
    components, candidates = _worker_load_components(payload)
    result: dict[str, Any] = {
        "success": True,
        "versions": versions,
        "python": {
            "version": platform.python_version(),
            "executable": sys.executable,
        },
        "openfe_version": getattr(openfe, "__version__", versions.get("openfe")),
        "receptor": receptor,
        "candidates": candidates,
        "network": None,
    }
    if payload["action"] == "plan":
        if len(components) < int(payload["min_ligands"]):
            result.update(
                success=False,
                error={
                    "code": "insufficient_openfe_components",
                    "message": (
                        "OpenFE could construct only "
                        f"{len(components)} usable ligand components"
                    ),
                },
            )
        else:
            result["network"] = _worker_plan_network(
                payload,
                components,
                candidates,
            )
    return result


def _worker_main(payload_file: Path) -> int:
    try:
        payload = json.loads(payload_file.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("worker payload must be a JSON object")
        result = _worker(payload)
        returncode = 0 if result.get("success") else 2
    except Exception as exc:
        result = {
            "success": False,
            "error": {
                "code": "openfe_worker_exception",
                "type": type(exc).__name__,
                "message": str(exc),
            },
        }
        returncode = 2
    print(
        WORKER_RESULT_PREFIX
        + json.dumps(result, sort_keys=True, separators=(",", ":"))
    )
    return returncode


def _merge_worker_candidates(
    candidates: list[dict[str, Any]],
    worker_candidates: Sequence[Mapping[str, Any]],
) -> None:
    by_pose = {
        str(entry.get("pose_id")): entry
        for entry in worker_candidates
        if entry.get("pose_id")
    }
    for candidate in candidates:
        if candidate.get("status") != "passed":
            continue
        outcome = by_pose.get(str(candidate.get("pose_id")))
        if outcome is None:
            candidate.update(
                status="failed",
                reason_code="openfe_worker_result_missing",
                reason="OpenFE worker did not return this canonical pose",
            )
        elif outcome.get("status") != "passed":
            candidate.update(
                status="failed",
                reason_code=outcome.get("reason_code") or "openfe_worker_failed",
                reason=outcome.get("reason") or "OpenFE component conversion failed",
            )


def _network_blockers(
    network: Mapping[str, Any],
    *,
    allow_charge_changes: bool,
    min_mapping_score: float,
    min_mapped_heavy_fraction: float,
) -> list[dict[str, Any]]:
    blockers: list[dict[str, Any]] = []
    for edge in network.get("edges") or []:
        edge_id = (
            f"{edge.get('source_pose_id')}->{edge.get('target_pose_id')}"
        )
        if (
            not allow_charge_changes
            and edge.get("source_formal_charge")
            != edge.get("target_formal_charge")
        ):
            blockers.append(
                {
                    "code": "charge_changing_edge",
                    "edge": edge_id,
                    "reason": "Default RBFE route forbids net-charge-changing edges",
                }
            )
        score = edge.get("lomap_score")
        if not isinstance(score, (int, float)) or not math.isfinite(score):
            blockers.append(
                {
                    "code": "invalid_mapping_score",
                    "edge": edge_id,
                    "reason": "Mapping edge has no finite LOMAP score",
                }
            )
        elif score < min_mapping_score:
            blockers.append(
                {
                    "code": "low_mapping_score",
                    "edge": edge_id,
                    "reason": (
                        f"LOMAP score {score:.3f} is below "
                        f"{min_mapping_score:.3f}"
                    ),
                }
            )
        fraction = edge.get("mapped_heavy_fraction")
        if not isinstance(fraction, (int, float)) or not math.isfinite(
            fraction
        ):
            blockers.append(
                {
                    "code": "invalid_mapped_heavy_fraction",
                    "edge": edge_id,
                    "reason": "Mapping edge has no finite heavy-atom fraction",
                }
            )
        elif fraction < min_mapped_heavy_fraction:
            blockers.append(
                {
                    "code": "low_mapped_heavy_fraction",
                    "edge": edge_id,
                    "reason": (
                        f"Mapped heavy-atom fraction {fraction:.3f} is below "
                        f"{min_mapped_heavy_fraction:.3f}"
                    ),
                }
            )
    return blockers


def run_openfe_rbfe(
    receptor_file: str | Path,
    candidates_sdf: str | Path,
    output_dir: str | Path,
    runtime: Mapping[str, Any],
    *,
    execution_mode: str = "preflight",
    network_method: str = "minimal_spanning",
    mapper: str = "KartografAtomMapper",
    mapper_settings: Mapping[str, Any] | None = None,
    max_candidates: int = 24,
    min_ligands: int = 2,
    allow_charge_changes: bool = False,
    has_metal: bool = False,
    is_covalent: bool = False,
    timeout_seconds: int = 600,
    n_processes: int = 1,
    mst_num: int = 2,
    central_pose_id: str | None = None,
    min_mapping_score: float = 0.3,
    min_mapped_heavy_fraction: float = 0.5,
    require_explicit_hydrogens: bool = True,
    allow_undefined_stereochemistry: bool = False,
) -> dict[str, Any]:
    """Preflight an RBFE cohort or create a mapping-only OpenFE network.

    The returned manifest always states that simulation was not executed.
    Scientific blockers use ``not_evaluable``; runtime/conversion failures use
    ``failed``.  Malformed paths or parameter types raise
    :class:`OpenFERBFEError`.
    """
    mode = str(execution_mode).strip().casefold()
    if mode not in EXECUTION_MODES:
        raise OpenFERBFEError(
            f"execution_mode must be one of {sorted(EXECUTION_MODES)}"
        )
    method = str(network_method).strip()
    if method not in NETWORK_METHODS:
        raise OpenFERBFEError(
            f"network_method must be one of {sorted(NETWORK_METHODS)}"
        )
    mapper_name = str(mapper).strip()
    if mapper_name not in MAPPERS:
        raise OpenFERBFEError(f"mapper must be one of {sorted(MAPPERS)}")
    try:
        maximum = int(max_candidates)
        minimum = int(min_ligands)
        processes = int(n_processes)
        redundant_msts = int(mst_num)
        timeout = int(timeout_seconds)
        minimum_score = float(min_mapping_score)
        minimum_fraction = float(min_mapped_heavy_fraction)
    except (TypeError, ValueError) as exc:
        raise OpenFERBFEError("numeric OpenFE parameters are invalid") from exc
    if not 2 <= minimum <= maximum <= 200:
        raise OpenFERBFEError(
            "require 2 <= min_ligands <= max_candidates <= 200"
        )
    if not 1 <= processes <= 32 or timeout < 1 or redundant_msts < 1:
        raise OpenFERBFEError(
            "n_processes, timeout_seconds, and mst_num are out of bounds"
        )
    if not 0.0 <= minimum_score <= 1.0:
        raise OpenFERBFEError("min_mapping_score must be between 0 and 1")
    if not 0.0 <= minimum_fraction <= 1.0:
        raise OpenFERBFEError(
            "min_mapped_heavy_fraction must be between 0 and 1"
        )
    if not isinstance(runtime, Mapping):
        raise OpenFERBFEError("runtime must be a mapping")
    if mapper_settings is not None and not isinstance(
        mapper_settings, Mapping
    ):
        raise OpenFERBFEError("mapper_settings must be a mapping or null")
    if central_pose_id is not None and not isinstance(central_pose_id, str):
        raise OpenFERBFEError("central_pose_id must be a string or null")
    central_id = (
        central_pose_id.strip() if central_pose_id is not None else None
    )
    if central_pose_id is not None and not central_id:
        raise OpenFERBFEError("central_pose_id must not be empty")

    receptor_path = _require_file(receptor_file, "receptor_file")
    candidates_path = _require_file(candidates_sdf, "candidates_sdf")
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "openfe_rbfe_manifest.json"
    payload_path = destination / "openfe_worker_input.json"
    graphml_path = destination / "ligand_network.graphml"

    candidates = _load_openfe_candidates(
        candidates_path,
        max_candidates=maximum,
        require_explicit_hydrogens=bool(require_explicit_hydrogens),
        allow_undefined_stereochemistry=bool(
            allow_undefined_stereochemistry
        ),
    )
    usable = [entry for entry in candidates if entry["status"] == "passed"]

    duplicate_ids: dict[str, list[dict[str, Any]]] = {}
    for entry in usable:
        duplicate_ids.setdefault(str(entry["chemical_id"]), []).append(entry)
    for duplicates in duplicate_ids.values():
        for duplicate in duplicates[1:]:
            duplicate.update(
                status="not_evaluable",
                reason_code="duplicate_chemical_identity",
                reason=(
                    "RBFE network already contains this canonical chemical "
                    "identity; only one pose per chemical entity is allowed"
                ),
            )
    usable = [entry for entry in candidates if entry["status"] == "passed"]
    if (
        method == "radial"
        and central_id is not None
        and central_id
        not in {
            str(entry["pose_id"])
            for entry in candidates
            if entry.get("pose_id")
        }
    ):
        raise OpenFERBFEError(
            "central_pose_id does not match any candidate pose: "
            f"{central_id}"
        )

    receptor_risks = inspect_receptor_risks(receptor_path)
    blockers: list[dict[str, str]] = []
    if is_covalent:
        blockers.append(
            {
                "code": "covalent_system",
                "reason": "Default OpenFE RBFE route does not parameterize covalent ligands",
            }
        )
    if has_metal or receptor_risks["metal_atoms"]:
        blockers.append(
            {
                "code": "metal_site",
                "reason": "Default OpenFE RBFE route requires explicit metal-site parameterization",
            }
        )
    unsupported_hetero_residues = receptor_risks[
        "unsupported_hetero_residues"
    ]
    if unsupported_hetero_residues:
        blockers.append(
            {
                "code": "unsupported_receptor_hetero_residues",
                "reason": (
                    "Default OpenFE RBFE route does not silently discard or "
                    "parameterize non-water receptor HETATM residues: "
                    + ", ".join(unsupported_hetero_residues)
                ),
            }
        )
    charges = sorted({int(entry["formal_charge"]) for entry in usable})
    if len(charges) > 1 and not allow_charge_changes:
        blockers.append(
            {
                "code": "mixed_formal_charge",
                "reason": (
                    "Default RBFE route forbids a cohort spanning formal "
                    f"charges {charges}"
                ),
            }
        )
    if len(usable) < minimum:
        blockers.append(
            {
                "code": "insufficient_usable_ligands",
                "reason": (
                    f"RBFE requires at least {minimum} unique usable ligands; "
                    f"found {len(usable)}"
                ),
            }
        )
    _global_block(candidates, blockers)

    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "success": True,
        "status": "not_evaluable",
        "execution_mode": mode,
        "result_semantics": {
            "simulation_executed": False,
            "free_energy_estimates_produced": False,
            "planning_level": (
                "applicability_preflight"
                if mode == "preflight"
                else "ligand_mapping_network"
            ),
            "simulation_ready": False,
        },
        "protocol": {
            "backend": "OpenFE",
            "network_method": method,
            "mapper": mapper_name,
            "mapper_settings": dict(mapper_settings or {}),
            "max_candidates": maximum,
            "min_ligands": minimum,
            "allow_charge_changes": bool(allow_charge_changes),
            "n_processes": processes,
            "mst_num": redundant_msts,
            "central_pose_id": central_id,
            "quality_gates": {
                "min_mapping_score": minimum_score,
                "min_mapped_heavy_fraction": minimum_fraction,
            },
            "input_requirements": {
                "explicit_hydrogens": bool(require_explicit_hydrogens),
                "defined_stereochemistry": not bool(
                    allow_undefined_stereochemistry
                ),
            },
        },
        "inputs": {
            "receptor_file": str(receptor_path),
            "receptor_sha256": _sha256(receptor_path),
            "candidates_sdf": str(candidates_path),
            "candidates_sdf_sha256": _sha256(candidates_path),
        },
        "applicability": {
            "has_metal_declared": bool(has_metal),
            "is_covalent": bool(is_covalent),
            "receptor_risks": receptor_risks,
            "global_blockers": blockers,
        },
        "tool_provenance": {
            "environment": runtime.get("environment", "current"),
            "runtime": dict(runtime),
            "versions": {},
            "worker": None,
        },
        "candidates": [],
        "counts": {},
        "network": None,
        "manifest_file": str(manifest_path),
    }
    if blockers:
        manifest["candidates"] = [
            _public_candidate(entry) for entry in candidates
        ]
        manifest["counts"] = _counts(manifest["candidates"])
        manifest["status"] = _overall_status(manifest["candidates"])
        _atomic_json(manifest_path, manifest)
        return manifest

    eligible = sorted(
        (
            {
                "record_index": int(entry["record_index"]),
                "pose_id": entry["pose_id"],
                "chemical_id": entry["chemical_id"],
                "formal_charge": int(entry["formal_charge"]),
            }
            for entry in usable
        ),
        key=lambda entry: (entry["chemical_id"], entry["pose_id"]),
    )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "action": mode,
        "receptor_file": str(receptor_path),
        "candidate_sdf": str(candidates_path),
        "eligible_candidates": eligible,
        "min_ligands": minimum,
        "network_method": method,
        "mapper": mapper_name,
        "mapper_settings": dict(mapper_settings or {}),
        "n_processes": processes,
        "mst_num": redundant_msts,
        "central_pose_id": central_id,
        "graphml_file": str(graphml_path),
    }
    worker_result, execution = _invoke_worker(
        runtime,
        payload,
        payload_path,
        timeout_seconds=timeout,
    )
    manifest["tool_provenance"]["worker"] = _bounded_invocation(execution)
    if worker_result is None:
        for candidate in candidates:
            if candidate["status"] == "passed":
                candidate.update(
                    status="failed",
                    reason_code="openfe_worker_output_invalid",
                    reason="OpenFE worker returned no structured result",
                )
    else:
        manifest["tool_provenance"]["versions"] = dict(
            worker_result.get("versions") or {}
        )
        manifest["tool_provenance"]["python"] = dict(
            worker_result.get("python") or {}
        )
        manifest["tool_provenance"]["receptor_validation"] = dict(
            worker_result.get("receptor") or {}
        )
        _merge_worker_candidates(
            candidates,
            worker_result.get("candidates") or [],
        )
        if not worker_result.get("success"):
            error = worker_result.get("error") or {}
            for candidate in candidates:
                if candidate["status"] == "passed":
                    candidate.update(
                        status="failed",
                        reason_code=error.get("code") or "openfe_worker_failed",
                        reason=error.get("message") or "OpenFE worker failed",
                    )
        network = worker_result.get("network")
        if isinstance(network, dict):
            quality_blockers = _network_blockers(
                network,
                allow_charge_changes=bool(allow_charge_changes),
                min_mapping_score=minimum_score,
                min_mapped_heavy_fraction=minimum_fraction,
            )
            network["quality_gate"] = {
                "status": "passed" if not quality_blockers else "not_evaluable",
                "blockers": quality_blockers,
            }
            manifest["network"] = network
            if quality_blockers:
                manifest["applicability"]["network_blockers"] = (
                    quality_blockers
                )

    manifest["candidates"] = [_public_candidate(entry) for entry in candidates]
    manifest["counts"] = _counts(manifest["candidates"])
    manifest["status"] = _overall_status(manifest["candidates"])
    if (
        manifest["status"] == "passed"
        and manifest["network"]
        and manifest["network"]["quality_gate"]["blockers"]
    ):
        manifest["status"] = "not_evaluable"
    manifest["result_semantics"]["simulation_ready"] = False
    _atomic_json(manifest_path, manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args(argv)
    if args.worker is None:
        parser.error("--worker is required when executing this module")
    return _worker_main(args.worker)


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "EXECUTION_MODES",
    "MAPPERS",
    "NETWORK_METHODS",
    "OpenFERBFEError",
    "SCHEMA_VERSION",
    "run_openfe_rbfe",
]

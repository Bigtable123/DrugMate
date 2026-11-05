#!/usr/bin/env python3
"""Normalize user-supplied SBDD inputs before pipeline config assembly.

This is an Agent-facing pre-config adapter, not a numbered scientific pipeline
step.  External receptor formats are normalized to PDB and fixed screening
SMI libraries are converted to SDF with Open Babel.  The original
inputs are never modified.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
import threading
import time
import traceback
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import closing
from dataclasses import replace
from pathlib import Path
from typing import Any


SKILL_DIR = Path(__file__).resolve().parent
BIO_ROOT = SKILL_DIR.parent
SHARED_DIR = BIO_ROOT / "shared"
sys.path.insert(0, str(SHARED_DIR))

import pipeline_registry as registry
from ligand_geometry import audit_sdf_geometry, geometry_failure_message
from process_lifecycle import PipelineCancelled, cancellation_scope
from tool_runtime import ToolRuntimeError, prepare_tool, run_prepared_tool, run_tool


SCHEMA_VERSION = "sbdd-input-normalization/v1"
DEFAULT_CONVERSION_TIMEOUT_SECONDS = 60
GEN3D_SPEED = "fast"
OPENBABEL_RUNTIME = {
    "kind": "command",
    "environment": "drugmate",
    "command": "obabel",
}
RECEPTOR_FORMATS = {"auto", "pdb", "mol2"}
RECEPTOR_STATES = {"raw", "prepared"}


class InputNormalizationError(ValueError):
    """Raised when an input cannot be normalized without losing its contract."""


def _json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )


def _write_progress(path: Path, progress: dict) -> None:
    """Publish a complete snapshot for readers polling a running conversion."""
    temporary = path.with_name(path.name + ".tmp")
    _json_dump(temporary, {**progress, "updated_at_unix": time.time()})
    temporary.replace(path)


def _finish_progress(
    path: Path, status: str, error: str | None = None, library: dict | None = None,
) -> None:
    progress = json.loads(path.read_text(encoding="utf-8"))
    progress.update(status=status, running_records=[])
    if library is not None:
        progress.update(
            total_records=library["input_record_count"],
            completed_records=library["input_record_count"],
            retained_records=library["retained_record_count"],
            skipped_records=library["skipped_record_count"],
        )
    if error:
        progress["error"] = error
    _write_progress(path, progress)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _input_file(value: Any, label: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise InputNormalizationError(f"{label}.path must be a non-empty file path")
    path = Path(str(value)).expanduser().resolve()
    if not path.is_file():
        raise InputNormalizationError(f"{label} file does not exist: {path}")
    return path


def _declared_choice(config: dict, key: str, choices: set[str], default: str) -> str:
    value = str(config.get(key, default)).strip().lower()
    if value not in choices:
        raise InputNormalizationError(
            f"{key} must be one of {sorted(choices)}, got {value!r}"
        )
    return value


def _reject_source_output_overlap(
    sources: list[Path | None], outputs: list[Path], label: str,
) -> None:
    source_paths = {path.resolve() for path in sources if path is not None}
    overlaps = sorted(
        str(path.resolve()) for path in outputs if path.resolve() in source_paths
    )
    if overlaps:
        raise InputNormalizationError(
            f"{label} output would overwrite a source file: {overlaps}"
        )


def detect_receptor_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".pdb", ".ent"}:
        return "pdb"
    if suffix == ".mol2":
        return "mol2"
    head = path.read_text(encoding="utf-8", errors="replace")[:65536]
    if "@<TRIPOS>MOLECULE" in head and "@<TRIPOS>ATOM" in head:
        return "mol2"
    if any(line.startswith(("ATOM  ", "HETATM")) for line in head.splitlines()):
        return "pdb"
    raise InputNormalizationError(
        f"cannot detect receptor format for {path}; supported formats are PDB and MOL2"
    )


def _pdb_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            coordinate = tuple(float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
        except (TypeError, ValueError) as exc:
            raise InputNormalizationError(
                f"invalid PDB coordinates at {path}:{line_number}"
            ) from exc
        if not all(math.isfinite(value) for value in coordinate):
            raise InputNormalizationError(
                f"non-finite PDB coordinates at {path}:{line_number}"
            )
        atoms.append(
            {
                "record": line[:6].strip(),
                "atom_name": line[12:16].strip(),
                "residue_name": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "residue_number": line[22:26].strip(),
                "coordinate": coordinate,
            }
        )
    if not atoms:
        raise InputNormalizationError(f"receptor contains no PDB ATOM/HETATM records: {path}")
    return atoms


def _mol2_residue_identity(label: str, subtype: str | None) -> tuple[str, str, str]:
    match = re.fullmatch(r"([A-Za-z][A-Za-z0-9]{0,2})(-?\d+)([A-Za-z]?)", label)
    if not match:
        raise InputNormalizationError(
            f"MOL2 substructure name does not encode a PDB residue identity: {label!r}"
        )
    label_name, residue_number, insertion_code = match.groups()
    residue_name = (
        str(subtype).strip()
        if subtype and str(subtype).strip() not in {"", "****"}
        else label_name
    )
    if len(residue_name) > 3:
        raise InputNormalizationError(
            f"MOL2 residue name cannot be represented in PDB: {residue_name!r}"
        )
    return residue_name.upper(), residue_number, insertion_code.upper()


def _mol2_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    substructures: dict[str, dict[str, str]] = {}
    in_substructure_section = False
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if line.upper() == "@<TRIPOS>SUBSTRUCTURE":
            in_substructure_section = True
            continue
        if in_substructure_section and line.upper().startswith("@<TRIPOS>"):
            break
        if not in_substructure_section or not line:
            continue
        fields = line.split()
        if len(fields) < 4:
            raise InputNormalizationError(
                f"invalid MOL2 substructure record at {path}:{line_number}"
            )
        residue_name, residue_number, insertion_code = _mol2_residue_identity(
            fields[1], fields[6] if len(fields) > 6 else None
        )
        substructures[fields[0]] = {
            "label": fields[1],
            "type": fields[3],
            "chain": fields[5] if len(fields) > 5 and fields[5] != "****" else "",
            "residue_name": residue_name,
            "residue_number": residue_number,
            "insertion_code": insertion_code,
        }

    in_atom_section = False
    for line_number, raw_line in enumerate(
        lines, start=1
    ):
        line = raw_line.strip()
        if line.upper() == "@<TRIPOS>ATOM":
            in_atom_section = True
            continue
        if in_atom_section and line.upper().startswith("@<TRIPOS>"):
            break
        if not in_atom_section or not line:
            continue
        fields = line.split()
        if len(fields) < 6:
            raise InputNormalizationError(f"invalid MOL2 atom record at {path}:{line_number}")
        try:
            coordinate = tuple(float(fields[index]) for index in (2, 3, 4))
        except ValueError as exc:
            raise InputNormalizationError(
                f"invalid MOL2 coordinates at {path}:{line_number}"
            ) from exc
        if not all(math.isfinite(value) for value in coordinate):
            raise InputNormalizationError(
                f"non-finite MOL2 coordinates at {path}:{line_number}"
            )
        substructure_id = fields[6] if len(fields) > 6 else None
        substructure_name = fields[7] if len(fields) > 7 else None
        substructure = substructures.get(str(substructure_id), {})
        if not substructure and substructure_name:
            residue_name, residue_number, insertion_code = _mol2_residue_identity(
                substructure_name, None
            )
            substructure = {
                "label": substructure_name,
                "type": "unknown",
                "chain": "",
                "residue_name": residue_name,
                "residue_number": residue_number,
                "insertion_code": insertion_code,
            }
        if not substructure:
            raise InputNormalizationError(
                f"MOL2 atom {fields[0]} has no resolvable substructure identity"
            )
        atoms.append(
            {
                "atom_id": fields[0],
                "atom_name": fields[1],
                "coordinate": coordinate,
                "atom_type": fields[5],
                "substructure_id": substructure_id,
                "substructure_name": substructure_name,
                "substructure_type": substructure.get("type"),
                "residue_name": substructure["residue_name"],
                "chain": substructure["chain"],
                "residue_number": substructure["residue_number"],
                "insertion_code": substructure["insertion_code"],
            }
        )
    if not atoms:
        raise InputNormalizationError(f"receptor contains no MOL2 atom records: {path}")
    return atoms


def _restore_mol2_residue_identities(path: Path, source_atoms: list[dict]) -> None:
    """Restore MOL2 residue identity after Open Babel's sequential renumbering."""
    output_lines = []
    atom_index = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            output_lines.append(line)
            continue
        if atom_index >= len(source_atoms):
            raise InputNormalizationError("PDB contains more atoms than the source MOL2")
        source_atom = source_atoms[atom_index]
        output_atom_name = line[12:16].strip()
        if output_atom_name != source_atom["atom_name"]:
            raise InputNormalizationError(
                "Open Babel changed MOL2 atom identity/order at atom "
                f"{atom_index + 1}: {source_atom['atom_name']!r} -> {output_atom_name!r}"
            )
        try:
            residue_number = int(source_atom["residue_number"])
        except (TypeError, ValueError) as exc:
            raise InputNormalizationError(
                f"invalid MOL2 residue number: {source_atom['residue_number']!r}"
            ) from exc
        if not -999 <= residue_number <= 9999:
            raise InputNormalizationError(
                f"MOL2 residue number is outside the PDB field range: {residue_number}"
            )
        identity = (
            f"{source_atom['residue_name']:>3} "
            f"{str(source_atom['chain'])[:1]:1}"
            f"{residue_number:4d}"
            f"{str(source_atom['insertion_code'])[:1]:1}"
        )
        output_lines.append(line[:17] + identity + line[27:])
        atom_index += 1
    if atom_index != len(source_atoms):
        raise InputNormalizationError(
            f"PDB contains fewer atoms than the source MOL2: {atom_index} vs {len(source_atoms)}"
        )
    path.write_text("".join(output_lines), encoding="utf-8")


def _runtime_record(result: dict) -> dict:
    return {
        key: result.get(key)
        for key in (
            "success",
            "command",
            "returncode",
            "elapsed_seconds",
            "environment",
            "error",
        )
        if result.get(key) is not None
    }


def _numeric_triplet(value: Any, label: str, *, positive: bool = False) -> list[float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise InputNormalizationError(f"{label} must contain exactly three numbers")
    try:
        triplet = [float(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise InputNormalizationError(f"{label} must contain exactly three numbers") from exc
    if not all(math.isfinite(item) for item in triplet):
        raise InputNormalizationError(f"{label} contains non-finite values")
    if positive and not all(item > 0 for item in triplet):
        raise InputNormalizationError(f"{label} values must be positive")
    return triplet


def _normalize_docking_box(config: dict, output_dir: Path) -> dict:
    allowed = {"path"}
    unknown = sorted(set(config) - allowed)
    if unknown:
        raise InputNormalizationError(f"unknown docking_box normalization keys: {unknown}")
    source = _input_file(config.get("path"), "docking_box")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise InputNormalizationError(f"docking box is not valid JSON: {source}") from exc
    if not isinstance(payload, dict):
        raise InputNormalizationError("docking box JSON root must be an object")
    if payload.get("unit", "angstrom") != "angstrom":
        raise InputNormalizationError("docking box unit must be angstrom; no implicit unit conversion")
    center_key = "center" if "center" in payload else "center_angstrom"
    size_key = "size" if "size" in payload else "size_angstrom"
    center = _numeric_triplet(payload.get(center_key), f"docking_box.{center_key}")
    size = _numeric_triplet(
        payload.get(size_key), f"docking_box.{size_key}", positive=True
    )
    padding = payload.get("padding", payload.get("padding_angstrom"))
    if padding is not None:
        try:
            padding = float(padding)
        except (TypeError, ValueError) as exc:
            raise InputNormalizationError("docking box padding must be numeric") from exc
        if not math.isfinite(padding) or padding < 0:
            raise InputNormalizationError("docking box padding must be finite and non-negative")

    box_dir = output_dir / "docking_box"
    box_dir.mkdir(parents=True, exist_ok=True)
    normalized = box_dir / "docking_box.json"
    normalized_payload = {
        "schema_version": 1,
        "center": center,
        "size": size,
        "unit": "angstrom",
        "source": payload.get("source", "user_provided_box"),
        "padding": padding,
        "source_file": str(source),
    }
    _json_dump(normalized, normalized_payload)
    pocket_definition = box_dir / "pocket_definition.json"
    _json_dump(
        pocket_definition,
        {
            "schema_version": "provided-docking-box-pocket/v1",
            "method": "provided_docking_box",
            "center": center,
            "size": size,
            "unit": "angstrom",
            "source": payload.get("source", "user_provided_box"),
            "source_file": str(source),
        },
    )
    return {
        "source_file": str(source),
        "source_sha256": _sha256(source),
        "normalized_file": str(normalized.resolve()),
        "normalized_sha256": _sha256(normalized),
        "pocket_definition_file": str(pocket_definition.resolve()),
        "center": center,
        "size": size,
        "unit": "angstrom",
        "source": normalized_payload["source"],
    }


def _normalize_receptor(config: dict, output_dir: Path) -> dict:
    allowed = {"path", "format", "preparation_state", "coordinate_tolerance"}
    unknown = sorted(set(config) - allowed)
    if unknown:
        raise InputNormalizationError(f"unknown receptor normalization keys: {unknown}")

    source = _input_file(config.get("path"), "receptor")
    declared_format = _declared_choice(config, "format", RECEPTOR_FORMATS, "auto")
    detected_format = detect_receptor_format(source)
    if declared_format != "auto" and declared_format != detected_format:
        raise InputNormalizationError(
            f"declared receptor format {declared_format!r} does not match detected "
            f"format {detected_format!r} for {source}"
        )
    receptor_format = detected_format if declared_format == "auto" else declared_format
    state = _declared_choice(
        config, "preparation_state", RECEPTOR_STATES, "raw"
    )
    try:
        coordinate_tolerance = float(config.get("coordinate_tolerance", 0.002))
    except (TypeError, ValueError) as exc:
        raise InputNormalizationError("receptor.coordinate_tolerance must be numeric") from exc
    if coordinate_tolerance < 0 or coordinate_tolerance > 0.1:
        raise InputNormalizationError(
            "receptor.coordinate_tolerance must be between 0 and 0.1 angstrom"
        )

    receptor_dir = output_dir / "receptor"
    receptor_dir.mkdir(parents=True, exist_ok=True)
    normalized = receptor_dir / "normalized_receptor.pdb"
    _reject_source_output_overlap([source], [normalized], "receptor normalization")
    runtime_result = None
    if receptor_format == "pdb":
        source_atoms = _pdb_atoms(source)
        shutil.copyfile(source, normalized)
    else:
        source_atoms = _mol2_atoms(source)
        runtime_result = run_tool(
            OPENBABEL_RUNTIME,
            ["-imol2", source, "-opdb", "-O", normalized],
            cwd=receptor_dir,
            timeout=300,
        )
        if not normalized.is_file() or normalized.stat().st_size == 0:
            raise InputNormalizationError("Open Babel produced no normalized receptor PDB")
        _restore_mol2_residue_identities(normalized, source_atoms)

    normalized_atoms = _pdb_atoms(normalized)
    if len(source_atoms) != len(normalized_atoms):
        raise InputNormalizationError(
            "receptor atom count changed during format normalization: "
            f"{len(source_atoms)} -> {len(normalized_atoms)}"
        )
    if receptor_format == "mol2":
        identity_mismatches = [
            {
                "atom_index": index,
                "source": {
                    key: source_atom.get(key)
                    for key in (
                        "atom_name",
                        "residue_name",
                        "chain",
                        "residue_number",
                    )
                },
                "normalized": {
                    key: normalized_atom.get(key)
                    for key in (
                        "atom_name",
                        "residue_name",
                        "chain",
                        "residue_number",
                    )
                },
            }
            for index, (source_atom, normalized_atom) in enumerate(
                zip(source_atoms, normalized_atoms), start=1
            )
            if any(
                str(source_atom.get(key, "")) != str(normalized_atom.get(key, ""))
                for key in ("atom_name", "residue_name", "chain", "residue_number")
            )
        ]
        if identity_mismatches:
            raise InputNormalizationError(
                "receptor atom/residue identity changed during MOL2 normalization; "
                f"first mismatch: {identity_mismatches[0]}"
            )
    maximum_coordinate_delta = max(
        abs(source_atom["coordinate"][axis] - normalized_atom["coordinate"][axis])
        for source_atom, normalized_atom in zip(source_atoms, normalized_atoms)
        for axis in range(3)
    )
    if maximum_coordinate_delta > coordinate_tolerance:
        raise InputNormalizationError(
            "receptor coordinates changed beyond the allowed normalization tolerance: "
            f"max_delta={maximum_coordinate_delta:.6f} A, "
            f"tolerance={coordinate_tolerance:.6f} A"
        )

    chains = sorted({atom["chain"] for atom in normalized_atoms if atom["chain"]})
    residue_keys = {
        (atom["chain"], atom["residue_name"], atom["residue_number"])
        for atom in normalized_atoms
    }
    warnings = []
    if not chains:
        warnings.append(
            "normalized receptor has no chain identifiers; context.chain must be reviewed"
        )
    if receptor_format == "mol2":
        warnings.append(
            "PDB does not preserve every MOL2 bond-type or partial-charge field; "
            "the original MOL2 remains the provenance source"
        )

    prepared_artifacts = {}
    if state == "prepared":
        from receptor_preparation import write_noh, water_inventory, file_identity, write_json
        no_h = receptor_dir / "prepared_receptor_noH.pdb"
        write_noh(normalized, no_h, remove_waters=True)
        prepared_artifacts = {
            "prepared_receptor_noH": str(no_h.resolve()),
            "source_complex_pdb": str(normalized.resolve()),
            "water_inventory": water_inventory(normalized, receptor_dir / "water_inventory.json"),
            "prep_manifest": write_json(receptor_dir / "prep_manifest.json", {
                "schema_version": 1, "stage": "external_prepared_geometry",
                "source": file_identity(normalized), "prepared_receptor_noH": file_identity(no_h),
                "preparation_evidence": "user declared prepared; no geometry repair performed",
                "final_protonation": "deferred_to_docking_prep",
                "qc_scope": "format/coordinate preservation only; no valence or clash certification",
            }),
        }
    return {
        "source_file": str(source),
        "source_format": receptor_format,
        "prepared_artifacts": prepared_artifacts,
        "preparation_state": state,
        "normalized_file": str(normalized.resolve()),
        "source_sha256": _sha256(source),
        "normalized_sha256": _sha256(normalized),
        "source_atom_count": len(source_atoms),
        "normalized_atom_count": len(normalized_atoms),
        "atom_and_residue_identity_preserved": True,
        "residue_count": len(residue_keys),
        "chains": chains,
        "maximum_coordinate_delta_angstrom": round(maximum_coordinate_delta, 6),
        "coordinate_tolerance_angstrom": coordinate_tolerance,
        "warnings": warnings,
        "runtime": _runtime_record(runtime_result) if runtime_result else None,
    }


def _sdf_records(path: Path):
    """Split at record boundaries without rewriting the original SDF bytes."""
    record = []
    with path.open("rb") as handle:
        for line in handle:
            record.append(line)
            if line.strip() == b"$$$$":
                yield b"".join(record)
                record = []
    if record and b"".join(record).strip():
        yield b"".join(record)


def _record_coordinate_error(record: bytes, scratch: Path) -> str | None:
    if not record.rstrip().endswith(b"$$$$"):
        return "incomplete_sdf_record"
    scratch.write_bytes(record)
    qc = audit_sdf_geometry(scratch)
    if qc["invalid_records"]:
        return qc["invalid_records"][0]["reason"]
    if qc["record_count"] != 1:
        return "unparseable_record"
    return None


def _smiles_rows(source: Path) -> list[tuple[dict, str]]:
    rows = []
    with source.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.strip().split(maxsplit=1)
            rows.append(({
                "source_molecule_index": len(rows),
                "record_number": len(rows) + 1,
                "source_line_number": line_number,
                "name": fields[1] if len(fields) > 1 else None,
            }, line.rstrip("\r\n") + "\n"))
    return rows


def _conversion_timeout(config: dict) -> int:
    timeout = config.get("conversion_timeout_seconds")
    # Older requests may explicitly contain null. They now use the same finite
    # default as omitted values rather than silently disabling the deadline.
    if timeout is None:
        timeout = DEFAULT_CONVERSION_TIMEOUT_SECONDS
    if type(timeout) is not int or timeout <= 0:
        raise InputNormalizationError(
            "library.conversion_timeout_seconds must be a positive integer; "
            "omitted or null uses the 60-second default"
        )
    return timeout


def _conversion_options(config: dict, count: int) -> dict:
    requested = config.get("workers", 8)
    if type(requested) is not int or not 1 <= requested <= 64:
        raise InputNormalizationError("library.workers must be an integer from 1 to 64")
    timeout = _conversion_timeout(config)
    try:
        cpus = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        cpus = os.cpu_count() or 1
    return {
        "requested_workers": requested,
        "workers": min(requested, max(1, cpus), count),
        "threads_per_worker": 1,
        "gen3d_speed": GEN3D_SPEED,
        "conversion_timeout_seconds": timeout,
    }


def _prepare_fast_converter(prepared, runtime_dir: Path, options: dict):
    """Detect the Open Babel 3.1 fast writeback defect before screening a library.

    A local driver uses the same native implementation with its missing copy
    back restored. It never replaces the installed Open Babel or changes speed.
    """
    runtime_dir.mkdir(parents=True, exist_ok=True)
    probe_input = runtime_dir / "fast_probe.smi"
    probe_input.write_text("CC sbdd_fast_probe\n", encoding="utf-8")

    def probe(tool, name):
        output = runtime_dir / f"{name}.sdf"
        output.unlink(missing_ok=True)
        result = run_prepared_tool(
            tool, ["-ismi", probe_input, "-osdf", "-O", output,
                   "-h", "--gen3d", GEN3D_SPEED],
            cwd=runtime_dir, timeout=10, check=False,
        )
        records = list(_sdf_records(output)) if output.exists() else []
        reason = (
            _record_coordinate_error(records[0], runtime_dir / "probe_qc.sdf")
            if len(records) == 1 else "missing_probe_record"
        )
        return result, reason

    native_probe, reason = probe(prepared, "native_fast_probe")
    options["fast_runtime_probe"] = {**_runtime_record(native_probe), "geometry_error": reason}
    if native_probe["success"] and reason is None:
        options["command_backend"] = "openbabel_native_fast"
        return prepared
    if not native_probe["success"] or reason != "collapsed_coordinates":
        raise ToolRuntimeError("Open Babel fast-mode functional probe failed", native_probe)

    # This known defect returns success and writes the untouched, all-zero
    # input molecule. Do not mistake a broken runtime for 200 bad ligands.
    source_dir = SKILL_DIR / "scripts" / "openbabel_fast_compat"
    sources = sorted(source_dir.rglob("*.cpp"))
    prefix = Path(prepared.command).parent.parent
    include = prefix / "include"
    library = prefix / "lib" / "libopenbabel.so"
    eigen_candidates = [include / "eigen3", Path("/usr/include/eigen3"),
                        Path("/usr/local/include/eigen3")]
    # Conda runtime packages need not ship Eigen's development headers. A
    # sibling environment may already provide them (including TensorFlow's
    # bundled Eigen). This is header discovery only; no environment is changed.
    eigen_candidates.extend(sorted(prefix.parent.glob("*/include/eigen3")))
    eigen_candidates.extend(sorted(prefix.parent.glob(
        "*/lib/python*/site-packages/tensorflow/include",
    )))
    eigen_include = None
    for candidate in eigen_candidates:
        macros = candidate / "Eigen/src/Core/util/Macros.h"
        if (candidate / "Eigen/Core").is_file() and macros.is_file() and re.search(
            r"#\s*define\s+EIGEN_WORLD_VERSION\s+3\b",
            macros.read_text(encoding="utf-8"),
        ):
            eigen_include = candidate
            break
    compiler = shutil.which("g++") or shutil.which("c++")
    if not compiler or not sources or eigen_include is None or not all(path.is_file() for path in (
        include / "openbabel3/openbabel/mol.h", library,
    )):
        raise InputNormalizationError(
            "Installed Open Babel fast mode loses generated coordinates; the local "
            "writeback fix needs a C++17 compiler plus matching Open Babel/Eigen 3 headers "
            "and libopenbabel.so. No input molecules were submitted for conversion."
        )
    executable = runtime_dir / "openbabel_fast_compat"
    manifest = runtime_dir / "fast_compat_build.json"
    fingerprint = {
        "compiler": compiler, "openbabel_library_sha256": _sha256(library),
        "compile_options": ["-std=c++17", "-O1", "-DHAVE_EIGEN"],
        "eigen_include": str(eigen_include),
        "eigen_macros_sha256": _sha256(eigen_include / "Eigen/src/Core/util/Macros.h"),
        "sources": {str(path.relative_to(source_dir)): _sha256(path)
                    for path in sorted(source_dir.rglob("*"))
                    if path.is_file() and path.suffix in {".cpp", ".h"}},
    }
    cached = {}
    if executable.is_file() and manifest.is_file():
        try:
            cached = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass
    if cached.get("fingerprint") != fingerprint:
        build = run_tool(
            {"kind": "command", "environment": "current", "command": compiler},
            [*fingerprint["compile_options"],
             f"-I{include / 'openbabel3'}", f"-I{eigen_include}", f"-I{source_dir}",
             *sources, library, f"-Wl,-rpath,{prefix / 'lib'}", "-o", executable],
            cwd=runtime_dir, timeout=60,
        )
        _json_dump(manifest, {"fingerprint": fingerprint, "build": build})
    fixed = replace(prepared, command=str(executable))
    fixed_probe, fixed_reason = probe(fixed, "compat_fast_probe")
    if not fixed_probe["success"] or fixed_reason is not None:
        raise ToolRuntimeError("Open Babel fast compatibility probe failed", fixed_probe)
    options.update(
        command_backend="openbabel_fast_writeback_compat",
        compatibility_manifest=str(manifest),
        compatibility_probe=_runtime_record(fixed_probe),
    )
    return fixed


def _prepare_library_converter(options: dict, runtime_dir: Path):
    started = time.monotonic()
    prepared = prepare_tool(OPENBABEL_RUNTIME, env={
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    })
    prepared = _prepare_fast_converter(prepared, runtime_dir, options)
    options["environment_setup_seconds"] = round(time.monotonic() - started, 3)
    return prepared


def _convert_smiles_row(item, *, scratch_dir, source_format, prepared, timeout,
                        cancel_event):
    row, line = item
    # Each process owns its files; duplicate compound names cannot collide.
    with tempfile.TemporaryDirectory(
        prefix=f"record-{row['record_number']}-", dir=scratch_dir,
    ) as directory:
        work = Path(directory)
        smi_file, sdf_file = work / "source.smi", work / "converted.sdf"
        smi_file.write_text(line, encoding="utf-8")
        result = run_prepared_tool(
            prepared,
            [f"-i{source_format}", smi_file, "-osdf", "-O", sdf_file,
             "-h", "--gen3d", GEN3D_SPEED],
            cwd=work, timeout=timeout, check=False, cancel_event=cancel_event,
        )
        data = sdf_file.read_bytes() if sdf_file.exists() else b""
        records = list(_sdf_records(sdf_file)) if data else []
        if not result["success"]:
            reason = "conversion_timeout" if result.get("returncode") == 124 else "conversion_failed"
        elif not records:
            reason = "conversion_produced_no_record"
        elif len(records) != 1:
            reason = "unexpected_conversion_record_count"
        else:
            reason = _record_coordinate_error(records[0], work / "qc.sdf")
        if reason is None and "*** Open Babel Error" in result.get("stderr", ""):
            reason = "conversion_failed"
        return row, result, data, records[0] if reason is None else None, reason


def _convert_smiles(rows, *, options, scratch_dir, source_format, progress_path):
    """Observe completions immediately while yielding SDFs in source order."""
    started = time.monotonic()
    progress = {
        "status": "preparing_environment", "total_records": len(rows),
        "completed_records": 0, "retained_records": 0, "skipped_records": 0,
        "running_records": [], "workers": options["workers"],
        "conversion_timeout_seconds": options["conversion_timeout_seconds"],
        "gen3d_speed": GEN3D_SPEED,
    }
    _write_progress(progress_path, progress)
    if not rows:
        return
    prepared = _prepare_library_converter(options, progress_path.parent / "runtime")
    cancel_event = threading.Event()
    active = {}
    active_lock = threading.Lock()

    def convert(item):
        row, _line = item
        index = row["source_molecule_index"]
        with active_lock:
            active[index] = (row, time.monotonic())
        try:
            return _convert_smiles_row(
                item, scratch_dir=scratch_dir, source_format=source_format,
                prepared=prepared, timeout=options["conversion_timeout_seconds"],
                cancel_event=cancel_event,
            )
        finally:
            with active_lock:
                active.pop(index, None)

    def publish_progress():
        now = time.monotonic()
        with active_lock:
            progress["running_records"] = [
                {**row, "elapsed_seconds": round(now - record_started, 3)}
                for row, record_started in sorted(
                    active.values(), key=lambda item: item[0]["source_molecule_index"],
                )
            ]
        progress["elapsed_seconds"] = round(now - started, 3)
        _write_progress(progress_path, progress)

    # The main thread handles signals; workers clean their own process trees.
    # On setup/I/O failures or cancellation, do not wait for all queued molecules.
    with cancellation_scope():
        pool = ThreadPoolExecutor(max_workers=options["workers"])
        try:
            pending = {pool.submit(convert, item) for item in rows}
            ready = {}
            next_index = 0
            progress["status"] = "converting"
            publish_progress()
            while pending:
                done, pending = wait(pending, timeout=1, return_when=FIRST_COMPLETED)
                for future in done:
                    result = future.result()
                    row, _runtime, _data, _record, reason = result
                    ready[row["source_molecule_index"]] = result
                    progress["completed_records"] += 1
                    progress["skipped_records" if reason else "retained_records"] += 1
                # A slow early molecule must not hide later completions from
                # the Agent. Heartbeats also show which records are still busy.
                publish_progress()
                while next_index in ready:
                    yield ready.pop(next_index)
                    next_index += 1
        finally:
            cancel_event.set()
            pool.shutdown(wait=True, cancel_futures=True)
    options["elapsed_seconds"] = round(time.monotonic() - started, 3)


def _normalize_library(config: dict, output_dir: Path) -> dict:
    source = _input_file(config.get("path"), "library")
    source_format = str(config.get("format", "auto")).lower()
    if source_format == "auto":
        source_format = "sdf" if source.suffix.lower() in {".sdf", ".sd"} else "smi"
    normalized = output_dir / "library" / "library_3d.sdf"
    unfiltered = normalized.with_name("library_unfiltered.sdf")
    normalized.parent.mkdir(parents=True, exist_ok=True)
    _reject_source_output_overlap(
        [source], [normalized, unfiltered], "library normalization"
    )
    runtime = None
    conversion_execution = None
    conversion_runs = []
    skipped = []
    retained_indices = []
    input_count = 0
    # Clear prior outputs even when this attempt fails before writing a record.
    normalized.write_bytes(b"")
    unfiltered.write_bytes(b"")
    if source_format in {"sdf", "sd"}:
        shutil.copy2(source, unfiltered)
    elif source_format not in {"smi", "smiles", "can"}:
        conversion_execution = {
            "gen3d_speed": GEN3D_SPEED,
            "conversion_timeout_seconds": _conversion_timeout(config),
            "timeout_scope": "file",
        }
        prepared = _prepare_library_converter(conversion_execution, output_dir / "runtime")
        runtime = run_prepared_tool(
            prepared,
            [f"-i{source_format}", source, "-osdf", "-O", unfiltered,
             "-h", "--gen3d", GEN3D_SPEED],
            cwd=normalized.parent,
            timeout=conversion_execution["conversion_timeout_seconds"],
        )

    with tempfile.TemporaryDirectory(prefix="normalize-", dir=normalized.parent) as directory:
        scratch_dir = Path(directory)
        qc_file = scratch_dir / "record.sdf"
        with normalized.open("wb") as accepted:
            if source_format in {"smi", "smiles", "can"}:
                rows = _smiles_rows(source)
                input_count = len(rows)
                conversion_execution = _conversion_options(config, input_count)
                with unfiltered.open("wb") as raw, closing(_convert_smiles(
                    rows, options=conversion_execution, scratch_dir=scratch_dir,
                    source_format=source_format,
                    progress_path=output_dir / "normalization_progress.json",
                )) as conversions:
                    for row, result, data, record, reason in conversions:
                        conversion_runs.append({**row, **result})
                        raw.write(data)
                        if reason:
                            skipped.append({**row, "reason": reason})
                        else:
                            accepted.write(record)
                            retained_indices.append(row["source_molecule_index"])
            else:
                for index, record in enumerate(_sdf_records(unfiltered)):
                    input_count += 1
                    reason = _record_coordinate_error(record, qc_file)
                    if reason:
                        skipped.append({
                            "source_molecule_index": index,
                            "record_number": index + 1,
                            "name": record.splitlines()[0].decode("utf-8", errors="replace") or None,
                            "reason": reason,
                        })
                    else:
                        accepted.write(record)
                        retained_indices.append(index)

    warnings = []
    if skipped:
        warnings.append(
            f"Skipped {len(skipped)}/{input_count} input molecule(s) during normalization; "
            "see library.skipped_records. Continue with the retained molecules without resampling."
        )
    return {
        "source_file": str(source),
        "source_format": source_format,
        "normalized_file": str(normalized),
        "unfiltered_file": str(unfiltered),
        "runtime": runtime,
        "conversion_execution": conversion_execution,
        "conversion_runs": conversion_runs,
        "input_record_count": input_count,
        "retained_record_count": len(retained_indices),
        "retained_source_indices": retained_indices,
        "skipped_record_count": len(skipped),
        "skipped_records": skipped,
        "warnings": warnings,
        "geometry_qc": audit_sdf_geometry(normalized),
    }


def _disabled_step_patch(step_id: str) -> dict:
    return {
        "enabled": False,
        "tools": {
            tool_id: {"enabled": False}
            for tool_id in registry.STEPS[step_id]["tools"]
        },
    }


def build_config_patch(
    receptor_result: dict | None,
    library_result: dict | None,
    docking_box_result: dict | None = None,
) -> tuple[dict, list[str]]:
    """Build a merge patch plus human-readable plan requirements.

    The returned object is intentionally not a complete runnable config.  It is
    merged into the Agent's complete explicit matrix before normal validation.
    """
    patch: dict[str, Any] = {}
    requirements: list[str] = []
    if receptor_result:
        receptor_file = receptor_result["normalized_file"]
        if receptor_result["preparation_state"] == "prepared":
            artifacts = receptor_result.get("prepared_artifacts") or {}
            if not artifacts.get("prepared_receptor_noH"):
                raise InputNormalizationError("Prepared receptor normalization must publish geometry-only artifacts; rerun normalization")
            patch.setdefault("input_artifacts", {}).update(artifacts)
            patch.setdefault("steps", {})["structure_prep"] = _disabled_step_patch(
                "structure_prep"
            )
            requirements.append(
                "Keep structure_prep and all of its tools disabled because the user "
                "declared the normalized receptor prepared. Use its separate noH derivative; "
                "final protonation and water selection still run in docking_prep."
            )
        else:
            patch.setdefault("context", {})["receptor_file"] = receptor_file
            patch.setdefault("steps", {}).setdefault(
                "structure_prep", {"tools": {}}
            ).setdefault("tools", {})["intake"] = {"enabled": False}
            requirements.append(
                "Run structure_prep in prep-only mode with intake disabled; the Agent "
                "must still select every other Step 0 tool explicitly."
            )
    if library_result:
        patch.setdefault("input_artifacts", {})["library"] = library_result[
            "normalized_file"
        ]
        patch.setdefault("steps", {})["molecule_generation"] = _disabled_step_patch(
            "molecule_generation"
        )
        requirements.append(
            "Keep molecule_generation and all of its tools disabled because a fixed "
            "SDF library is injected."
        )
        if library_result.get("skipped_record_count"):
            requirements.append(
                "Carry normalization_report.json library.skipped_records into the final "
                "report as preprocessing failures with null scores. Report input, retained "
                "and skipped counts separately; do not resample to replace failed molecules."
            )
    if docking_box_result:
        patch.setdefault("input_artifacts", {})["docking_box"] = docking_box_result[
            "normalized_file"
        ]
        pocket = _disabled_step_patch("pocket_analysis")
        pocket.update(enabled=True, params={"site_selection_strategy": "provided_box"})
        patch.setdefault("steps", {})["pocket_analysis"] = pocket
        patch["steps"]["docking"] = {"params": {"search_space_source": "docking_box"}}
        requirements.append(
            "Run pocket_analysis with site_selection_strategy='provided_box' in the same round; "
            "it characterizes the supplied box without replacing it. Detectors and interaction "
            "profilers are optional explicit selections. Keep docking.search_space_source='docking_box'."
        )
        prepared_receptor = (patch.get("input_artifacts") or {}).get("prepared_receptor_noH")
        if prepared_receptor:
            sys.path.insert(0, str(BIO_ROOT / "pocket-analysis" / "scripts"))
            from ligand_site_utils import build_box_site
            try:
                site = build_box_site(Path(prepared_receptor), docking_box_result)
            except ValueError as exc:
                message = str(exc)
                docking_box_result["chain_selection"] = {"source": "review_required", "error": message}
                docking_box_result.setdefault("warnings", []).append(message)
                requirements.append(message)
            else:
                patch.setdefault("context", {})["chain"] = site["chain"]
                docking_box_result["chain_selection"] = site["chain_selection"]
                requirements.append(
                    "context.chain is inferred from the prepared receptor's unique box overlap. "
                    "Preserve any user-specified chain when merging this patch and validate its overlap; "
                    "the inferred chain must not override an explicit target selection."
                )
    return patch, requirements


def run(input_data: dict, output_dir: Path) -> dict:
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "normalization_report.json"
    progress_path = output_dir / "normalization_progress.json"
    config_patch_path = output_dir / "config_patch.json"
    # A failed retry must not leave a previous success patch consumable.
    config_patch_path.unlink(missing_ok=True)
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "success": False,
        "receptor": None,
        "library": None,
        "docking_box": None,
        "warnings": [],
        "normalization_progress": str(progress_path),
    }
    _json_dump(report_path, report)
    _write_progress(progress_path, {
        "status": "starting", "total_records": None, "completed_records": 0,
        "retained_records": 0, "skipped_records": 0, "running_records": [],
    })
    try:
        if not isinstance(input_data, dict):
            raise InputNormalizationError("input root must be a JSON object")
        unknown = sorted(set(input_data) - {"receptor", "library", "docking_box"})
        if unknown:
            raise InputNormalizationError(f"unknown input normalization keys: {unknown}")
        receptor_config = input_data.get("receptor")
        library_config = input_data.get("library")
        docking_box_config = input_data.get("docking_box")
        if receptor_config is None and library_config is None and docking_box_config is None:
            raise InputNormalizationError(
                "provide at least one of receptor, library, or docking_box"
            )
        if receptor_config is not None and not isinstance(receptor_config, dict):
            raise InputNormalizationError("receptor must be an object")
        if library_config is not None and not isinstance(library_config, dict):
            raise InputNormalizationError("library must be an object")
        if docking_box_config is not None and not isinstance(docking_box_config, dict):
            raise InputNormalizationError("docking_box must be an object")

        if receptor_config is not None:
            report["receptor"] = _normalize_receptor(receptor_config, output_dir)
            report["warnings"].extend(report["receptor"].get("warnings") or [])
        if library_config is not None:
            report["library"] = _normalize_library(library_config, output_dir)
            report["warnings"].extend(report["library"].get("warnings") or [])
            if not report["library"]["retained_record_count"]:
                raise InputNormalizationError(
                    "No usable molecules remain after input normalization; "
                    "see library.skipped_records for individual failures."
                )
            if not report["library"]["geometry_qc"]["success"]:
                raise InputNormalizationError(
                    geometry_failure_message(report["library"]["geometry_qc"])
                )
        if docking_box_config is not None:
            report["docking_box"] = _normalize_docking_box(
                docking_box_config, output_dir
            )

        config_patch, requirements = build_config_patch(
            report["receptor"], report["library"], report["docking_box"]
        )
        report["warnings"].extend((report["docking_box"] or {}).get("warnings") or [])
        _json_dump(config_patch_path, config_patch)
        report.update(
            {
                "success": True,
                "config_patch": config_patch,
                "config_patch_file": str(config_patch_path),
                "plan_requirements": requirements,
            }
        )
    except (PipelineCancelled, KeyboardInterrupt) as exc:
        report.update(error_code="input_normalization_cancelled", error=str(exc))
        _json_dump(report_path, report)
        _finish_progress(progress_path, "cancelled", str(exc))
        raise
    except (InputNormalizationError, ToolRuntimeError, OSError, ValueError, RuntimeError) as exc:
        report.update(
            {
                "success": False,
                "error_code": "input_normalization_failed",
                "error": f"{type(exc).__name__}: {exc}",
            }
        )
        if isinstance(exc, ToolRuntimeError) and exc.result:
            report["failed_runtime"] = _runtime_record(exc.result)
    _json_dump(report_path, report)
    _finish_progress(
        progress_path, "completed" if report["success"] else "failed",
        report.get("error"), report["library"],
    )
    report["normalization_report"] = str(report_path)
    return report


def main() -> None:
    from io_utils import parse_skill_args

    try:
        input_data, output_dir = parse_skill_args("SBDD Input Normalization Skill")
        result = run(input_data, output_dir)
    except Exception as exc:  # pragma: no cover - last-resort CLI boundary
        traceback.print_exc(file=sys.stderr)
        result = {
            "schema_version": SCHEMA_VERSION,
            "success": False,
            "error_code": "input_normalization_unhandled_error",
            "error": f"{type(exc).__name__}: {exc}",
        }
    print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
    raise SystemExit(0 if result.get("success") else 1)


if __name__ == "__main__":
    main()

"""AmberTools MMPBSA.py single-frame endpoint reranking backend.

The backend consumes only the auditable ``candidate_topologies.json`` artifact
written by :mod:`candidate_topology`.  ``preflight`` checks the registered
runtime and topology bundle without calculating energies.  ``score`` runs one
independent MMPBSA.py job per eligible pose using that pose's
``complex.inpcrd`` as a single frame.

These values are deliberately labelled low-reliability, no-entropy endpoint
reranking signals.  They are neither ensemble averages nor absolute binding
free energies, and they must not be averaged with docking or RBFE metrics.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

from tool_runtime import ToolRuntimeError, probe_tool, run_tool


SCHEMA_VERSION = "amber-mmgbsa-single-frame/v1"
EXECUTION_MODES = frozenset({"preflight", "score"})
SUPPORTED_GB_MODELS = frozenset({1, 2, 5, 7, 8})
REQUIRED_TOPOLOGY_FILES = (
    "candidate_sdf",
    "complex_prmtop",
    "complex_inpcrd",
    "receptor_prmtop",
    "receptor_inpcrd",
    "ligand_prmtop",
    "ligand_inpcrd",
)
_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_DELTA_LINE = re.compile(
    rf"^\s*(DELTA(?:\s+[A-Za-z][A-Za-z0-9_.()+/\-]*)+?)\s+({_FLOAT})(?:\s|$)",
    re.IGNORECASE,
)
_DIFFERENCES = re.compile(
    r"^\s*Differences\s*\(\s*Complex\s*-\s*Receptor\s*-\s*Ligand\s*\)\s*:\s*$",
    re.IGNORECASE | re.MULTILINE,
)
_MMPBSA_VERSION = re.compile(
    r"^\|?\s*MMPBSA\.py\s+Version\s*=\s*([^\s|]+)",
    re.IGNORECASE | re.MULTILINE,
)
_CANONICAL_POSE_ID = re.compile(r"^pose:[0-9a-f]{64}$")


class AmberMMGBSAError(RuntimeError):
    """Raised for malformed invocations rather than per-pose outcomes."""


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


def _safe_component(value: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    return cleaned[:80] or fallback


def _bounded_invocation(result: Mapping[str, Any] | None) -> dict[str, Any]:
    result = result or {}
    return {
        "command": [str(value) for value in (result.get("command") or [])],
        "environment": result.get("environment"),
        "returncode": result.get("returncode"),
        "elapsed_seconds": result.get("elapsed_seconds"),
        "stdout_tail": str(result.get("stdout") or "")[-4000:],
        "stderr_tail": str(result.get("stderr") or "")[-4000:],
    }


def _load_manifest(path: Path) -> dict[str, Any]:
    if path.name != "candidate_topologies.json":
        raise AmberMMGBSAError(
            "input artifact must be named candidate_topologies.json"
        )
    if not path.is_file():
        raise AmberMMGBSAError(
            f"candidate_topologies.json does not exist: {path}"
        )
    if path.stat().st_size == 0:
        raise AmberMMGBSAError(
            f"candidate_topologies.json is empty: {path}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise AmberMMGBSAError(
            f"candidate_topologies.json is not valid JSON: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise AmberMMGBSAError(
            "candidate_topologies.json root must be a JSON object"
        )
    if not isinstance(payload.get("candidates"), list):
        raise AmberMMGBSAError(
            "candidate_topologies.json must contain a candidates array"
        )
    if payload.get("schema_version") != "candidate-topologies/v1":
        raise AmberMMGBSAError(
            "unsupported candidate_topologies.json schema_version; "
            "expected candidate-topologies/v1"
        )
    return payload


def _resolve_artifact_path(value: Any, manifest_dir: Path) -> Path | None:
    if not isinstance(value, str) or not value.strip():
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_dir / path
    return path.resolve()


def _file_evidence(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
    }


def _prmtop_atom_count(path: Path) -> int:
    """Read NATOM from the first POINTERS value in an Amber prmtop."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip().upper() != "%FLAG POINTERS":
            continue
        for data_line in lines[index + 1:]:
            stripped = data_line.strip()
            if not stripped or stripped.startswith("%FORMAT"):
                continue
            try:
                count = int(stripped.split()[0])
            except (IndexError, ValueError) as exc:
                raise AmberMMGBSAError(
                    f"cannot parse NATOM from Amber topology {path}"
                ) from exc
            if count < 1:
                raise AmberMMGBSAError(
                    f"Amber topology has invalid NATOM={count}: {path}"
                )
            return count
    raise AmberMMGBSAError(
        f"Amber topology has no parseable %FLAG POINTERS section: {path}"
    )


def _inpcrd_atom_count(path: Path) -> int:
    """Read the atom count from an Amber ASCII inpcrd title/count header."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        try:
            count = int(stripped.split()[0])
        except (IndexError, ValueError) as exc:
            raise AmberMMGBSAError(
                f"cannot parse atom count from Amber coordinates {path}"
            ) from exc
        if count < 1:
            raise AmberMMGBSAError(
                f"Amber coordinates have invalid atom count={count}: {path}"
            )
        return count
    raise AmberMMGBSAError(
        f"Amber coordinates have no parseable atom-count line: {path}"
    )


def _validate_topology_bundle(paths: Mapping[str, Path]) -> dict[str, Any]:
    topology_counts = {
        component: _prmtop_atom_count(paths[f"{component}_prmtop"])
        for component in ("complex", "receptor", "ligand")
    }
    coordinate_counts = {
        component: _inpcrd_atom_count(paths[f"{component}_inpcrd"])
        for component in ("complex", "receptor", "ligand")
    }
    mismatched = [
        component
        for component in ("complex", "receptor", "ligand")
        if topology_counts[component] != coordinate_counts[component]
    ]
    if mismatched:
        details = ", ".join(
            f"{component} prmtop={topology_counts[component]} "
            f"inpcrd={coordinate_counts[component]}"
            for component in mismatched
        )
        raise AmberMMGBSAError(
            f"Amber topology/coordinate atom counts disagree: {details}"
        )
    expected_complex = topology_counts["receptor"] + topology_counts["ligand"]
    if topology_counts["complex"] != expected_complex:
        raise AmberMMGBSAError(
            "complex atom count does not equal receptor + ligand "
            f"({topology_counts['complex']} != "
            f"{topology_counts['receptor']} + {topology_counts['ligand']})"
        )
    return {
        "status": "passed",
        "topology_atom_counts": topology_counts,
        "coordinate_atom_counts": coordinate_counts,
        "complex_equals_receptor_plus_ligand": True,
    }


def _base_candidate(
    candidate: Any,
    *,
    index: int,
    manifest_dir: Path,
) -> dict[str, Any]:
    if not isinstance(candidate, dict):
        return {
            "record_index": index,
            "selection_rank": index + 1,
            "pose_id": None,
            "chemical_id": None,
            "status": "failed",
            "reason_code": "invalid_candidate_record",
            "reason": "candidate topology record must be a JSON object",
            "selected_for_mmgbsa": False,
            "topology_files": {},
            "topology_validation": None,
            "execution": None,
            "delta_terms": {},
            "pose_score": None,
        }

    pose_id = str(candidate.get("pose_id") or "").strip() or None
    record = {
        "record_index": candidate.get("record_index", index),
        "selection_rank": candidate.get("selection_rank", index + 1),
        "pose_id": pose_id,
        "chemical_id": candidate.get("chemical_id"),
        "source_identity_id": candidate.get("source_identity_id"),
        "parent_pose_id": candidate.get("parent_pose_id"),
        "coordinate_revision": candidate.get("coordinate_revision"),
        "status": "passed",
        "reason_code": None,
        "reason": None,
        "selected_for_mmgbsa": False,
        "topology_files": {},
        "topology_validation": None,
        "execution": None,
        "delta_terms": {},
        "pose_score": None,
    }
    upstream_status = str(candidate.get("status") or "").strip().casefold()
    if upstream_status == "not_evaluable":
        record.update(
            status="not_evaluable",
            reason_code=(
                candidate.get("reason_code")
                or "upstream_topology_not_evaluable"
            ),
            reason=(
                candidate.get("reason")
                or "candidate topology was scientifically not evaluable"
            ),
        )
        return record
    if upstream_status != "passed":
        record.update(
            status="failed",
            reason_code="upstream_candidate_topology_failed",
            reason=(
                candidate.get("reason")
                or f"candidate topology status is {upstream_status or 'missing'}"
            ),
        )
        return record
    if pose_id is None or _CANONICAL_POSE_ID.fullmatch(pose_id) is None:
        record.update(
            status="failed",
            reason_code="missing_canonical_pose_id",
            reason=(
                "passed topology record requires a canonical pose:<64 lowercase "
                "hex digits> pose_id"
            ),
        )
        return record

    output_files = candidate.get("output_files")
    if not isinstance(output_files, dict):
        record.update(
            status="failed",
            reason_code="missing_topology_bundle",
            reason="passed topology record has no output_files object",
        )
        return record

    resolved: dict[str, Path] = {}
    missing_fields: list[str] = []
    for name in REQUIRED_TOPOLOGY_FILES:
        path = _resolve_artifact_path(output_files.get(name), manifest_dir)
        if path is None:
            missing_fields.append(name)
        else:
            resolved[name] = path
    if missing_fields:
        record.update(
            status="failed",
            reason_code="missing_topology_bundle_fields",
            reason=(
                "passed topology record is missing output files: "
                + ", ".join(missing_fields)
            ),
        )
        return record

    coordinate = resolved["complex_inpcrd"]
    if coordinate.suffix.casefold() != ".inpcrd":
        record.update(
            status="not_evaluable",
            reason_code="unsupported_trajectory_input",
            reason=(
                "current MM-GBSA backend accepts only a single-frame "
                "complex.inpcrd, not a trajectory or ensemble"
            ),
        )
        return record

    missing_files = [
        name
        for name, path in resolved.items()
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing_files:
        record.update(
            status="failed",
            reason_code="topology_bundle_files_missing",
            reason=(
                "candidate topology manifest points to missing or empty files: "
                + ", ".join(missing_files)
            ),
        )
        return record

    try:
        record["topology_validation"] = _validate_topology_bundle(resolved)
    except (OSError, AmberMMGBSAError) as exc:
        record.update(
            status="failed",
            reason_code="topology_coordinate_validation_failed",
            reason=str(exc),
        )
        return record

    record["topology_files"] = {
        name: _file_evidence(path) for name, path in resolved.items()
    }
    return record


def _mark_duplicate_pose_ids(candidates: list[dict[str, Any]]) -> None:
    by_pose: dict[str, list[dict[str, Any]]] = {}
    for candidate in candidates:
        if candidate["status"] == "passed" and candidate.get("pose_id"):
            by_pose.setdefault(str(candidate["pose_id"]), []).append(candidate)
    for pose_id, records in by_pose.items():
        if len(records) < 2:
            continue
        for record in records:
            record.update(
                status="failed",
                reason_code="duplicate_canonical_pose_id",
                reason=(
                    "candidate_topologies.json contains duplicate canonical "
                    f"pose_id {pose_id}"
                ),
            )


def _select_top_k(candidates: list[dict[str, Any]], top_k: int) -> None:
    selected = 0
    for candidate in candidates:
        if candidate["status"] != "passed":
            continue
        if selected < top_k:
            candidate["selected_for_mmgbsa"] = True
            selected += 1
        else:
            candidate.update(
                status="not_evaluable",
                reason_code="outside_requested_top_k",
                reason=f"candidate is outside the requested top-{top_k} cohort",
            )


def _write_mmgbsa_input(
    path: Path,
    *,
    gb_model: int,
    salt_concentration_molar: float,
    keep_intermediate_files: bool,
) -> None:
    path.write_text(
        "# Single-frame, no-entropy endpoint reranking only.\n"
        "&general\n"
        "  startframe=1, endframe=1, interval=1,\n"
        f"  verbose=1, keep_files={1 if keep_intermediate_files else 0},\n"
        "/\n"
        "&gb\n"
        f"  igb={gb_model}, saltcon={salt_concentration_molar:.6f},\n"
        "/\n",
        encoding="utf-8",
    )


def parse_mmpbsa_differences(text: str) -> dict[str, float]:
    """Parse DELTA averages from exactly one MMPBSA.py Differences section."""
    markers = list(_DIFFERENCES.finditer(text))
    sections: list[dict[str, float]] = []
    for index, marker in enumerate(markers):
        end = markers[index + 1].start() if index + 1 < len(markers) else len(text)
        terms: dict[str, float] = {}
        for line in text[marker.end():end].splitlines():
            match = _DELTA_LINE.match(line)
            if not match:
                continue
            key = re.sub(
                r"[^a-z0-9]+",
                "_",
                match.group(1).strip().casefold(),
            ).strip("_")
            value = float(match.group(2))
            if not math.isfinite(value):
                raise AmberMMGBSAError(
                    f"MMPBSA.py emitted a non-finite value for {key}"
                )
            if key in terms:
                raise AmberMMGBSAError(
                    f"MMPBSA.py emitted duplicate {key} terms"
                )
            terms[key] = value
        if terms:
            sections.append(terms)
    complete = [section for section in sections if "delta_total" in section]
    if not complete:
        raise AmberMMGBSAError(
            "MMPBSA.py output has no parseable Differences DELTA TOTAL term"
        )
    if len(complete) != 1:
        raise AmberMMGBSAError(
            "MMPBSA.py output contains multiple DELTA TOTAL Differences sections"
        )
    return complete[0]


def _method_version(text: str, fallback: str) -> str:
    match = _MMPBSA_VERSION.search(text)
    return match.group(1) if match else fallback


def _pose_score(
    candidate: Mapping[str, Any],
    delta_total: float,
    *,
    method_version: str,
    primary: bool,
) -> dict[str, Any]:
    score = {
        "schema_version": 1,
        "metric_id": "amber_mmpbsa.delta_total_no_entropy",
        "value": float(delta_total),
        "unit": "kcal/mol",
        "direction": "lower",
        "method": "amber_mmpbsa",
        "method_version": method_version,
        "stage": "pose_rescoring",
        "applicability": "applicable",
        "reliability": "low",
        "pose_id": candidate["pose_id"],
        "interpretation": (
            "single-frame no-entropy endpoint reranking signal; not an "
            "absolute binding free energy"
        ),
    }
    if primary:
        score["primary"] = True
    return score


def _run_candidate(
    candidate: dict[str, Any],
    candidate_dir: Path,
    runtime: Mapping[str, Any],
    *,
    gb_model: int,
    salt_concentration_molar: float,
    timeout_seconds: int,
    keep_intermediate_files: bool,
    use_as_primary_metric: bool,
    method_version: str,
) -> None:
    candidate_dir.mkdir(parents=True, exist_ok=True)
    input_file = candidate_dir / "mmpbsa.in"
    result_file = candidate_dir / "FINAL_RESULTS_MMPBSA.dat"
    _write_mmgbsa_input(
        input_file,
        gb_model=gb_model,
        salt_concentration_molar=salt_concentration_molar,
        keep_intermediate_files=keep_intermediate_files,
    )
    paths = {
        name: evidence["path"]
        for name, evidence in candidate["topology_files"].items()
    }
    commands = dict(runtime.get("commands") or {})
    executable = commands.get("mmpbsa") or runtime.get("command") or "MMPBSA.py"
    arguments = [
        "-O",
        "-i",
        input_file,
        "-o",
        result_file,
        "-sp",
        paths["complex_prmtop"],
        "-cp",
        paths["complex_prmtop"],
        "-rp",
        paths["receptor_prmtop"],
        "-lp",
        paths["ligand_prmtop"],
        "-y",
        paths["complex_inpcrd"],
    ]
    execution: dict[str, Any] | None = None
    try:
        execution = run_tool(
            dict(runtime),
            arguments,
            command=str(executable),
            cwd=candidate_dir,
            timeout=timeout_seconds,
        )
        candidate["execution"] = _bounded_invocation(execution)
        if not result_file.is_file() or result_file.stat().st_size == 0:
            raise AmberMMGBSAError(
                "MMPBSA.py completed without a non-empty FINAL_RESULTS_MMPBSA.dat"
            )
        result_text = result_file.read_text(
            encoding="utf-8",
            errors="replace",
        )
        terms = parse_mmpbsa_differences(result_text)
    except (OSError, ToolRuntimeError, AmberMMGBSAError) as exc:
        failed_result = getattr(exc, "result", None)
        if failed_result:
            candidate["execution"] = _bounded_invocation(failed_result)
        elif candidate["execution"] is None and execution is not None:
            candidate["execution"] = _bounded_invocation(execution)
        candidate.update(
            status="failed",
            reason_code="mmgbsa_execution_failed",
            reason=str(exc),
        )
        candidate["input_file"] = (
            _file_evidence(input_file) if input_file.is_file() else None
        )
        if result_file.is_file() and result_file.stat().st_size:
            candidate["result_file"] = _file_evidence(result_file)
        return

    candidate["input_file"] = _file_evidence(input_file)
    candidate["result_file"] = _file_evidence(result_file)
    candidate["delta_terms"] = terms
    candidate["method_version"] = _method_version(result_text, method_version)
    candidate["pose_score"] = _pose_score(
        candidate,
        terms["delta_total"],
        method_version=candidate["method_version"],
        primary=use_as_primary_metric,
    )


def _counts(candidates: list[dict[str, Any]]) -> dict[str, int]:
    counts = {
        status: sum(candidate["status"] == status for candidate in candidates)
        for status in ("passed", "failed", "not_evaluable")
    }
    counts["requested"] = len(candidates)
    counts["selected"] = sum(
        bool(candidate.get("selected_for_mmgbsa")) for candidate in candidates
    )
    counts["scored"] = sum(candidate.get("pose_score") is not None for candidate in candidates)
    return counts


def _overall_status(candidates: list[dict[str, Any]]) -> tuple[str, str | None]:
    if any(candidate["status"] == "failed" for candidate in candidates):
        return "failed", "operational_failure"
    if not any(candidate["status"] == "passed" for candidate in candidates):
        return "not_evaluable", "no_evaluable_candidates"
    return "passed", None


def _ranking(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    ranked = sorted(
        (
            {
                "pose_id": candidate["pose_id"],
                "value": candidate["pose_score"]["value"],
                "unit": candidate["pose_score"]["unit"],
            }
            for candidate in candidates
            if candidate.get("pose_score") is not None
        ),
        key=lambda entry: (float(entry["value"]), str(entry["pose_id"])),
    )
    for rank, entry in enumerate(ranked, start=1):
        entry["rank"] = rank
    return {
        "metric_id": "amber_mmpbsa.delta_total_no_entropy",
        "direction": "lower",
        "scope": "independent_mmgbsa_rerank",
        "entries": ranked,
    }


def run_amber_mmgbsa(
    candidate_topologies: str | Path,
    output_dir: str | Path,
    runtime: Mapping[str, Any],
    *,
    execution_mode: str = "preflight",
    top_k: int = 20,
    gb_model: int = 5,
    salt_concentration_molar: float = 0.15,
    timeout_per_candidate_seconds: int = 900,
    fail_fast: bool = False,
    keep_intermediate_files: bool = False,
    use_as_primary_metric: bool = True,
) -> dict[str, Any]:
    """Preflight or score candidate-specific Amber topology bundles."""
    mode = str(execution_mode).strip().casefold()
    if mode not in EXECUTION_MODES:
        raise AmberMMGBSAError("execution_mode must be 'preflight' or 'score'")
    if isinstance(top_k, bool) or int(top_k) < 1:
        raise AmberMMGBSAError("top_k must be >= 1")
    if isinstance(gb_model, bool) or int(gb_model) not in SUPPORTED_GB_MODELS:
        raise AmberMMGBSAError(
            f"gb_model must be one of {sorted(SUPPORTED_GB_MODELS)}"
        )
    try:
        salt = float(salt_concentration_molar)
    except (TypeError, ValueError) as exc:
        raise AmberMMGBSAError(
            "salt_concentration_molar must be a finite non-negative number"
        ) from exc
    if not math.isfinite(salt) or salt < 0:
        raise AmberMMGBSAError(
            "salt_concentration_molar must be a finite non-negative number"
        )
    if (
        isinstance(timeout_per_candidate_seconds, bool)
        or int(timeout_per_candidate_seconds) < 1
    ):
        raise AmberMMGBSAError(
            "timeout_per_candidate_seconds must be >= 1"
        )

    source_path = Path(candidate_topologies).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "mmgbsa_manifest.json"
    source = _load_manifest(source_path)

    candidates = [
        _base_candidate(
            candidate,
            index=index,
            manifest_dir=source_path.parent,
        )
        for index, candidate in enumerate(source["candidates"])
    ]
    _mark_duplicate_pose_ids(candidates)
    _select_top_k(candidates, int(top_k))

    runtime_probe = probe_tool(dict(runtime), timeout=30)
    runtime_available = bool(runtime_probe.get("available"))
    if not runtime_available:
        for candidate in candidates:
            if candidate["status"] == "passed":
                candidate.update(
                    status="failed",
                    reason_code="mmgbsa_runtime_unavailable",
                    reason=str(
                        runtime_probe.get("detail")
                        or "registered MM-GBSA runtime is unavailable"
                    ),
                )

    method_version = str(
        runtime.get("version")
        or runtime_probe.get("version")
        or "unknown"
    )
    protocol_executed = False
    if mode == "score" and runtime_available:
        failure_seen = False
        for index, candidate in enumerate(candidates):
            if candidate["status"] != "passed" or not candidate["selected_for_mmgbsa"]:
                continue
            if failure_seen and fail_fast:
                candidate.update(
                    status="not_evaluable",
                    reason_code="aborted_after_operational_failure",
                    reason="batch stopped after an earlier MMPBSA.py failure",
                )
                continue
            protocol_executed = True
            candidate_dir = destination / (
                f"candidate_{index:05d}_"
                f"{_safe_component(str(candidate['pose_id']), f'pose-{index}')}"
            )
            _run_candidate(
                candidate,
                candidate_dir,
                runtime,
                gb_model=int(gb_model),
                salt_concentration_molar=salt,
                timeout_seconds=int(timeout_per_candidate_seconds),
                keep_intermediate_files=bool(keep_intermediate_files),
                use_as_primary_metric=bool(use_as_primary_metric),
                method_version=method_version,
            )
            failure_seen = candidate["status"] == "failed"

    status, reason_code = _overall_status(candidates)
    pose_scores = [
        candidate["pose_score"]
        for candidate in candidates
        if candidate.get("pose_score") is not None
    ]
    implementation_path = Path(__file__).resolve()
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "success": status != "failed",
        "status": status,
        "reason_code": reason_code,
        "execution_mode": mode,
        "result_semantics": {
            "protocol_executed": protocol_executed,
            "scores_produced": bool(pose_scores),
            "frames_per_pose": 1,
            "ensemble_averaging": False,
            "entropy_included": False,
            "absolute_binding_free_energy": False,
            "reliability": "low",
            "intended_use": "within-cohort_reranking_only",
        },
        "protocol": {
            "method": "AmberTools MMPBSA.py MM-GBSA",
            "gb_model": int(gb_model),
            "salt_concentration_molar": salt,
            "top_k": int(top_k),
            "fail_fast": bool(fail_fast),
            "keep_intermediate_files": bool(keep_intermediate_files),
            "use_as_primary_metric": bool(use_as_primary_metric),
            "trajectory_policy": "single_frame_complex_inpcrd_only",
            "entropy_method": None,
        },
        "inputs": {
            "candidate_topologies_manifest": str(source_path),
            "candidate_topologies_manifest_sha256": _sha256(source_path),
            "candidate_topologies_schema_version": source.get("schema_version"),
        },
        "tool_provenance": {
            "environment": runtime.get("environment", "current"),
            "registered_runtime": dict(runtime),
            "registered_commands": dict(runtime.get("commands") or {}),
            "runtime_probe": runtime_probe,
            "wrapper": {
                "path": str(implementation_path),
                "sha256": _sha256(implementation_path),
                "schema_version": SCHEMA_VERSION,
            },
        },
        "candidates": candidates,
        "pose_scores": pose_scores,
        "ranking": _ranking(candidates),
        "counts": _counts(candidates),
        "manifest_file": str(manifest_path),
    }
    _atomic_json(manifest_path, manifest)
    return manifest


__all__ = [
    "AmberMMGBSAError",
    "EXECUTION_MODES",
    "REQUIRED_TOPOLOGY_FILES",
    "SCHEMA_VERSION",
    "SUPPORTED_GB_MODELS",
    "parse_mmpbsa_differences",
    "run_amber_mmgbsa",
]

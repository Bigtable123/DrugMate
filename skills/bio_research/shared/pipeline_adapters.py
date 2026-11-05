"""Skill-input adapters (registry plan → each skill's run.py input).

The orchestrator is registry/artifact-driven; each skill's run.py takes a flat
input dict. These adapters translate the resolved run plan (wired artifact inputs
+ tools + params + context) into the dict each skill expects, and attach a
`tools` block ({tool_id: {"enabled": bool}}) which each skill reads via
io_utils.tool_enabled() to gate its per-tool switches.

Tool enablement is carried by the `tools` block. Standalone skill invocations
also accept flat flags as a fallback; these adapters emit the `tools` block.
"""
from __future__ import annotations

import hashlib
import json
import mmap
from pathlib import Path
from typing import Any

import pipeline_registry as R


COHORT_COMPONENT_FIELDS = (
    "protein_state_id",
    "pocket_id",
    "round_id",
    "protocol_profile_id",
    "metric_catalog_version",
)
COHORT_FIELDS = ("cohort_id", *COHORT_COMPONENT_FIELDS)
COHORT_PROVENANCE_FIELDS = (
    "cohort_identity_complete",
    "cohort_identity_missing_fields",
    "cohort_identity_derived_fields",
    "cohort_provenance_sources",
)
COHORT_INPUT_FIELDS = (*COHORT_FIELDS, *COHORT_PROVENANCE_FIELDS)
PROPERTY_PLAN_FIELDS = (
    "planned_metric_ids",
    "required_metric_ids",
    "required_admet_metric_ids",
    "not_applicable_metric_ids",
    "planned_category_codes",
    "planned_tools",
    "required_tools",
    "optional_tools",
    "failed_tools",
    "planned_but_failed_tools",
    "failed_metric_ids",
    "planned_but_failed_metric_ids",
    "operational_errors",
    "equivalent_fallbacks",
    "active_steps",
)


def _present(value: Any) -> bool:
    return value not in (None, "")


def _file_sha256(path_value: Any) -> str | None:
    if not _present(path_value):
        return None
    path = Path(str(path_value)).expanduser()
    try:
        if not path.is_file():
            return None
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def _read_small_json(path_value: Any) -> Any:
    """Read provenance JSON without loading multi-megabyte result reports."""
    if not _present(path_value):
        return None
    path = Path(str(path_value)).expanduser()
    try:
        if not path.is_file() or path.stat().st_size > 2 * 1024 * 1024:
            return None
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None


def _first_sdf_provenance(path_value: Any) -> dict[str, Any]:
    """Read only the first SDF record's cohort properties (no RDKit needed)."""
    if not _present(path_value):
        return {}
    path = Path(str(path_value)).expanduser()
    wanted = set(COHORT_INPUT_FIELDS) | {"protocol_id"}
    values: dict[str, Any] = {}
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            lines = iter(handle)
            for line in lines:
                if line.rstrip("\r\n") == "$$$$":
                    break
                if not line.startswith(">") or "<" not in line or ">" not in line:
                    continue
                start = line.find("<") + 1
                end = line.find(">", start)
                name = line[start:end].strip() if end > start else ""
                try:
                    value = next(lines).strip()
                except StopIteration:
                    break
                if name in wanted and value:
                    if name in COHORT_PROVENANCE_FIELDS:
                        try:
                            values.setdefault(name, json.loads(value))
                        except json.JSONDecodeError:
                            values.setdefault(name, value)
                    else:
                        values.setdefault(name, value)
    except OSError:
        return {}
    if not _present(values.get("protocol_profile_id")) and _present(
        values.get("protocol_id")
    ):
        values["protocol_profile_id"] = values["protocol_id"]
    return values


def _mapping_provenance(value: Any) -> dict[str, Any]:
    """Extract cohort fields from the small, explicit report locations."""
    values: dict[str, Any] = {}

    def visit(item: Any, depth: int = 0) -> None:
        if depth > 3 or not isinstance(item, dict):
            return
        for field in COHORT_INPUT_FIELDS:
            if _present(item.get(field)):
                values.setdefault(field, item[field])
        if not _present(values.get("protocol_profile_id")) and _present(
            item.get("protocol_id")
        ):
            values["protocol_profile_id"] = item["protocol_id"]
        for auxiliary in (
            "receptor_file", "pdb_id", "source_accession", "chain",
            "docking_box_file",
        ):
            if _present(item.get(auxiliary)):
                values.setdefault(auxiliary, item[auxiliary])
        for key in (
            "scientific_cohort",
            "parameters",
            "score_protocol",
            "pocket_definition",
        ):
            visit(item.get(key), depth + 1)
        for key in ("scores", "evaluation", "results", "molecules"):
            rows = item.get(key)
            if isinstance(rows, list) and rows and isinstance(rows[0], dict):
                visit(rows[0], depth + 1)
        cohorts = item.get("scientific_cohorts")
        # A singular run manifest is a safe compatibility projection.  Multiple
        # scientific cohorts must remain partitioned; choosing the first would
        # falsely label every downstream row with that cohort.
        if (
            isinstance(cohorts, list)
            and len(cohorts) == 1
            and isinstance(cohorts[0], dict)
        ):
            visit(cohorts[0], depth + 1)

    if isinstance(value, list) and value and isinstance(value[0], dict):
        visit(value[0])
    else:
        visit(value)
    return values


def _project_provenance_paths(
    rt: dict, *, include_docking: bool,
) -> list[Path]:
    project_value = rt.get("project_dir")
    if not _present(project_value) or "<" in str(project_value):
        return []
    project_dir = Path(str(project_value)).expanduser()
    paths = [
        project_dir / "step1" / "pocket_definition.json",
        project_dir / "step1" / "pocket_analysis_report.json",
    ]
    if include_docking:
        paths.extend([
            project_dir / "step3" / "docking_report.json",
            project_dir / "step3" / "docking_scores.json",
        ])
    return paths


def _scientific_cohort_manifest(path_value: Any) -> list[dict[str, Any]]:
    """Read a coordinator's small cohort manifest without a report-size cap.

    New coordinators write ``scientific_cohorts.json`` next to the large report.
    For older reports, extract that one top-level JSON array through ``mmap`` so
    a 20+ MB molecule audit does not force the provenance adapter to discard the
    cohort merely because ``_read_small_json`` intentionally caps its input.
    """
    if not _present(path_value):
        return []
    report_path = Path(str(path_value)).expanduser()
    candidates = [report_path.with_name("scientific_cohorts.json"), report_path]
    payload: Any = None
    for candidate in candidates:
        if not candidate.is_file():
            continue
        try:
            if candidate.name == "scientific_cohorts.json":
                payload = json.loads(candidate.read_text(encoding="utf-8"))
            else:
                payload = _mmap_top_level_json_value(
                    candidate, "scientific_cohorts",
                )
        except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError):
            payload = None
        if payload is not None:
            break
    if isinstance(payload, dict):
        payload = payload.get("scientific_cohorts")
    if not isinstance(payload, list):
        return []
    unique: list[dict[str, Any]] = []
    positions: dict[str, int] = {}
    for item in payload:
        if not isinstance(item, dict):
            continue
        item = _validated_cohort_manifest_item(item)
        fingerprint = json.dumps({
            field: item.get(field)
            for field in ("cohort_id", *COHORT_COMPONENT_FIELDS)
        }, sort_keys=True, default=str)
        position = positions.get(fingerprint)
        if position is None:
            positions[fingerprint] = len(unique)
            unique.append(dict(item))
            continue
        # Same scientific identity, richer metadata from another upstream
        # report.  Merge only fields the first source did not provide.
        for field, value in item.items():
            if not _present(unique[position].get(field)) and _present(value):
                unique[position][field] = value
    return unique


def _merge_scientific_cohort_manifests(
    *path_values: Any,
) -> list[dict[str, Any]]:
    """Merge upstream manifests without choosing a run-global cohort.

    Evaluation and property reports can legitimately contain more than one
    scientific cohort.  The adapter must carry every exact five-tuple into the
    coordinator; projecting the first row would let later identity joins and
    normalization cross protein/pocket/protocol boundaries.
    """
    merged: list[dict[str, Any]] = []
    positions: dict[str, int] = {}
    for path_value in path_values:
        for item in _scientific_cohort_manifest(path_value):
            fingerprint = json.dumps({
                field: item.get(field)
                for field in ("cohort_id", *COHORT_COMPONENT_FIELDS)
            }, sort_keys=True, default=str)
            position = positions.get(fingerprint)
            if position is None:
                positions[fingerprint] = len(merged)
                merged.append(dict(item))
                continue
            for field, value in item.items():
                if not _present(merged[position].get(field)) and _present(value):
                    merged[position][field] = value
    return merged


def _validated_cohort_manifest_item(item: dict[str, Any]) -> dict[str, Any]:
    """Annotate whether a manifest row's cohort_id agrees with its five tuple."""
    normalized = dict(item)
    missing = [
        field for field in COHORT_COMPONENT_FIELDS
        if not _present(normalized.get(field))
    ]
    errors = [f"missing:{field}" for field in missing]
    if not missing:
        canonical = "|".join(
            f"{field}={normalized[field]}" for field in COHORT_COMPONENT_FIELDS
        )
        declared = str(normalized.get("cohort_id") or "")
        declared_wrapper = (
            declared.startswith("declared_cohort_id=")
            and declared.endswith(f"|{canonical}")
        )
        if declared != canonical and not declared_wrapper:
            errors.append("cohort_id_component_mismatch")
    elif not _present(normalized.get("cohort_id")):
        errors.append("missing:cohort_id")
    normalized["cohort_manifest_valid"] = not errors
    normalized["cohort_manifest_validation_errors"] = errors
    return normalized


def _mmap_top_level_json_value(path: Path, key: str) -> Any:
    """Extract one JSON array/object value from a potentially large report."""
    marker = json.dumps(str(key)).encode("utf-8")
    with path.open("rb") as handle, mmap.mmap(
        handle.fileno(), 0, access=mmap.ACCESS_READ,
    ) as mapped:
        marker_at = mapped.find(marker)
        if marker_at < 0:
            return None
        colon = mapped.find(b":", marker_at + len(marker))
        if colon < 0:
            raise ValueError(f"JSON key {key!r} has no value separator")
        start = colon + 1
        while start < len(mapped) and chr(mapped[start]).isspace():
            start += 1
        if start >= len(mapped) or mapped[start] not in (ord("["), ord("{")):
            raise ValueError(f"JSON key {key!r} is not an array/object")
        opening = mapped[start]
        closing = ord("]") if opening == ord("[") else ord("}")
        depth = 0
        in_string = False
        escaped = False
        for offset in range(start, len(mapped)):
            byte = mapped[offset]
            if in_string:
                if escaped:
                    escaped = False
                elif byte == ord("\\"):
                    escaped = True
                elif byte == ord('"'):
                    in_string = False
                continue
            if byte == ord('"'):
                in_string = True
            elif byte == opening:
                depth += 1
            elif byte == closing:
                depth -= 1
                if depth == 0:
                    return json.loads(mapped[start:offset + 1].decode("utf-8"))
        raise ValueError(f"JSON key {key!r} has an unterminated value")


def _cohort_source_token(ctx: dict, rt: dict) -> str:
    source = {
        "project_dir": str(rt.get("project_dir") or ""),
        "upstream_project_dir": str(rt.get("upstream_project_dir") or ""),
        "round_number": rt.get("round_number"),
        "pdb_id": ctx.get("pdb_id"),
        "receptor_file": ctx.get("receptor_file"),
        "chain": ctx.get("chain", "A"),
    }
    digest = hashlib.sha256(
        json.dumps(source, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()
    return digest[:16]


def _scientific_cohort(
    wired: dict[str, Any],
    plan: dict,
    ctx: dict,
    rt: dict,
) -> dict[str, Any]:
    """Resolve the five-field scientific cohort without global unknown pools.

    New pipeline runs get explicit component fields on every downstream skill
    input.  Older/external artifacts remain accepted: known values are recovered
    from their JSON/SDF metadata, while a genuinely missing component receives
    a run-source-scoped value so unrelated unknown sources cannot normalize
    together.
    """
    recovered: dict[str, Any] = {}
    sources: dict[str, str] = {}

    def recover(field: str, value: Any, source: str) -> None:
        if not _present(value):
            return
        current = recovered.get(field)
        current_placeholder = str(current or "").startswith(
            ("unknown:", "pending:")
        ) or "|state_source:" in str(current or "")
        value_placeholder = str(value or "").startswith(
            ("unknown:", "pending:")
        ) or "|state_source:" in str(value or "")
        if not _present(current) or (current_placeholder and not value_placeholder):
            recovered[field] = value
            sources[field] = source

    for field in COHORT_FIELDS:
        if _present(ctx.get(field)):
            recover(field, ctx[field], "context")

    for key in ("selected_poses", "docked_poses", "prefiltered_library", "library"):
        for field, value in _first_sdf_provenance(wired.get(key)).items():
            recover(field, value, f"canonical_artifact:{key}:sdf")

    for key in (
        "pocket_definition", "pocket_analysis_report", "docking_report",
        "docking_scores", "evaluation_report", "property_report",
        "coordinator_report",
    ):
        candidate = wired.get(key)
        if isinstance(candidate, (dict, list)):
            loaded = candidate
        else:
            loaded = _read_small_json(candidate)
        if loaded is not None:
            for field, item in _mapping_provenance(loaded).items():
                recover(field, item, f"canonical_artifact:{key}")
    include_project_docking = any(
        _present(wired.get(key))
        for key in (
            "docked_poses", "docking_scores", "docking_report",
            "selected_poses", "evaluation_report", "property_report",
            "coordinator_report",
        )
    )
    for path in _project_provenance_paths(
        rt, include_docking=include_project_docking,
    ):
        loaded = _read_small_json(path)
        if loaded is not None:
            for field, item in _mapping_provenance(loaded).items():
                recover(
                    field, item, f"legacy_project_fallback:{path.name}"
                )

    token = _cohort_source_token(ctx, rt)
    receptor_file = (wired.get("prepared_receptor") or wired.get("prepared_receptor_noH")) or recovered.get(
        "receptor_file"
    ) or ctx.get("receptor_file")
    if not _present(receptor_file):
        project_value = rt.get("project_dir")
        if _present(project_value) and "<" not in str(project_value):
            candidates = sorted(
                (Path(str(project_value)).expanduser() / "step0").glob(
                    "*_receptor_prepared.pdb"
                )
            )
            receptor_file = str(candidates[0]) if candidates else None
            if receptor_file:
                sources.setdefault(
                    "receptor_file", "legacy_project_fallback:step0"
                )
    receptor_digest = _file_sha256(receptor_file)
    upstream_derived = recovered.get("cohort_identity_derived_fields")
    derived_fields: list[str] = (
        list(upstream_derived) if isinstance(upstream_derived, list) else []
    )
    missing_fields: list[str] = []
    if not _present(recovered.get("protein_state_id")):
        accession = ctx.get("pdb_id") or recovered.get("pdb_id")
        chain = ctx.get("chain") or recovered.get("chain") or "unspecified"
        if _present(accession):
            protein_state = f"pdb:{str(accession).upper()}|chain:{chain}"
            protein_state += (
                f"|receptor_sha256:{receptor_digest}"
                if receptor_digest
                else f"|state_source:{token}"
            )
        elif receptor_digest:
            protein_state = f"receptor_sha256:{receptor_digest}"
        else:
            protein_state = f"unknown:protein_state_id:{token}"
        recovered["protein_state_id"] = protein_state
        derived_fields.append("protein_state_id")
        sources["protein_state_id"] = (
            "derived:prepared_receptor_sha256"
            if receptor_digest
            else "derived:run_source_token"
        )
        if not receptor_digest:
            missing_fields.append("protein_state_id")

    if not _present(recovered.get("pocket_id")):
        box_value = wired.get("docking_box") or recovered.get("docking_box_file")
        box_digest = _file_sha256(box_value)
        recovered["pocket_id"] = (
            f"docking_box_sha256:{box_digest}"
            if box_digest
            else f"unknown:pocket_id:{token}"
        )
        derived_fields.append("pocket_id")
        sources["pocket_id"] = (
            "derived:docking_box_sha256"
            if box_digest
            else "derived:run_source_token"
        )
        if not box_digest:
            missing_fields.append("pocket_id")

    if not _present(recovered.get("round_id")):
        recovered["round_id"] = (
            rt.get("round_number")
            if _present(rt.get("round_number"))
            else f"unknown:round_id:{token}"
        )
        derived_fields.append("round_id")
        sources["round_id"] = (
            "runtime:round_number"
            if _present(rt.get("round_number"))
            else "derived:run_source_token"
        )
        if not _present(rt.get("round_number")):
            missing_fields.append("round_id")

    if not _present(recovered.get("protocol_profile_id")):
        selected_engines = [
            engine for engine in R.DOCKING_PRIMARY_ENGINES
            if (plan.get("tools", {}).get(engine) or {}).get("enabled")
        ]
        engine = selected_engines[0] if len(selected_engines) == 1 else "docking"
        recovered["protocol_profile_id"] = (
            f"pending:{engine}:{token}"
        )
        derived_fields.append("protocol_profile_id")
        missing_fields.append("protocol_profile_id")
        sources["protocol_profile_id"] = "pending:docking_execution"

    if not _present(recovered.get("metric_catalog_version")):
        recovered["metric_catalog_version"] = R.METRIC_RULES_VERSION
        derived_fields.append("metric_catalog_version")
        sources["metric_catalog_version"] = "registry:metric_catalog"

    components = {
        field: recovered[field] for field in COHORT_COMPONENT_FIELDS
    }
    missing_fields.extend(
        field for field, value in components.items()
        if str(value).startswith(("unknown:", "pending:"))
        or "|state_source:" in str(value)
    )
    components["cohort_id"] = "|".join(
        f"{field}={components[field]}" for field in COHORT_COMPONENT_FIELDS
    )
    components.update({
        "cohort_identity_complete": not missing_fields,
        "cohort_identity_missing_fields": sorted(set(missing_fields)),
        "cohort_identity_derived_fields": sorted(set(derived_fields)),
        "cohort_provenance_sources": {
            field: sources.get(field, "unknown")
            for field in COHORT_COMPONENT_FIELDS
        },
    })
    return components


def _screening_policy(step_plan: dict, ctx: dict) -> dict:
    """Pass the level-2 screening policy to any step that folds metric marks.

    ``screening_mode`` is the explicit choice (V1/V2/V3) and ``candidate_source``
    only decides the default when it is absent, so a named mode always wins
.  ``_drop_none`` in each adapter removes both keys when
    neither is set, leaving the skill on its own V2 default.
    """
    return {
        "screening_mode": _sparams(step_plan).get("screening_mode"),
        "candidate_source": ctx.get("candidate_source"),
    }


def _drop_none(payload: dict) -> dict:
    return {k: v for k, v in payload.items() if v is not None}


def _enabled(step_plan: dict, tool_id: str) -> bool:
    return bool((step_plan["tools"].get(tool_id) or {}).get("enabled", False))


def _tparams(step_plan: dict, tool_id: str) -> dict:
    return dict((step_plan["tools"].get(tool_id) or {}).get("params") or {})


def _truntime(step_plan: dict, tool_id: str) -> dict:
    return dict((step_plan["tools"].get(tool_id) or {}).get("runtime") or {})


def _sparams(step_plan: dict) -> dict:
    return dict(step_plan.get("params") or {})


def _required_sparam(step_plan: dict, name: str) -> Any:
    """Return a required resolved step parameter.

    Pipeline adapters receive a resolved plan, not an Agent-authored partial
    config.  Selection modes therefore must already be materialized by config
    resolution; silently inferring one from the artifacts would reintroduce
    routing inside the execution path.
    """
    params = _sparams(step_plan)
    if name not in params or params[name] is None:
        raise ValueError(f"resolved step plan is missing required param {name!r}")
    return params[name]


def _tool_configs(step_plan: dict) -> dict:
    """Emit enabled/params/runtime for interface-v2 skills.

    Every adapter attaches this; each skill reads the switches it owns and
    ignores the rest. Runtime metadata is registry-owned and therefore appears
    only after config resolution, never as an Agent-editable config field.
    """
    configs = {}
    for tid, tool in (step_plan.get("tools") or {}).items():
        enabled = bool(tool.get("enabled", False))
        config = {
            "enabled": enabled,
            "params": dict(tool.get("params") or {}),
        }
        # Disabled tools must not make resume fingerprints depend on large
        # checkpoints that were neither probed nor executed.
        if enabled:
            config["runtime"] = dict(tool.get("runtime") or {})
        configs[tid] = config
    return configs


def to_legacy_input(
    step_id: str,
    wired: dict[str, Any],
    step_plan: dict,
    context: dict,
    runtime: dict,
) -> dict:
    """Build the skill input for one step (flat legacy fields + `tools` block).

    Args:
        wired: {canonical artifact: path/value} resolved from upstream producers.
        step_plan: resolved step dict (tools{enabled,params}, params, ...).
        context: identity/shared values (pdb_id, chain, seed, target ids, known_*).
        runtime: {project_dir, upstream_project_dir, round_number, inferred_target_name}.
    """
    fn = _ADAPTERS.get(step_id)
    if fn is None:
        raise KeyError(f"no adapter for step {step_id!r}")
    payload = fn(wired, step_plan, context, runtime)
    # Attach interface-v2 switches to every step uniformly.
    payload["tools"] = _tool_configs(step_plan)
    return payload


# --- structure_prep -----------------------------------------------------------
def _structure_prep(wired, plan, ctx, rt) -> dict:
    intake = _tparams(plan, "intake")
    water = _tparams(plan, "water_selection")
    pdbfixer = _tparams(plan, "pdbfixer")
    propka = _tparams(plan, "propka")
    ph = pdbfixer.get("ph")
    if ph is None:
        ph = propka.get("ph")
    if ph is None:
        ph = _sparams(plan).get("ph")
    return _drop_none({
        "pdb_id": ctx.get("pdb_id"),
        "protein_name": ctx.get("protein_name"),
        "chain": ctx.get("chain"),
        "receptor_file": ctx.get("receptor_file"),
        "ligand_file": ctx.get("ligand_file"),
        "complex_mmcif_file": ctx.get("complex_mmcif_file"),
        "skip_prep": intake.get("skip_prep"),
        "water_distance_threshold": water.get("ligand_distance", intake.get("water_distance_threshold")),
        "water_protein_distance": water.get("protein_distance"),
        "water_min_occupancy": water.get("min_occupancy"),
        "water_max_b_factor": water.get("max_b_factor"),
        "water_include_residues": water.get("include_residues"),
        "water_exclude_residues": water.get("exclude_residues"),
        "keep_auxiliary_hetero": intake.get("keep_auxiliary_hetero"),
        "preserve_hetero_resnames": intake.get("preserve_hetero_resnames"),
        "ph": ph,
        "add_missing_residues": pdbfixer.get("add_missing_residues"),
        "add_missing_atoms": pdbfixer.get("add_missing_atoms"),
        "max_missing_residues": pdbfixer.get("max_missing_residues"),
        "repair_reason": pdbfixer.get("repair_reason"),
        "seed": pdbfixer.get("seed"),
        "add_hydrogens": pdbfixer.get("add_hydrogens"),
        "ligand_comp_id": _tparams(plan, "ligand_sdf").get("ligand_comp_id"),
        "required_residues": _tparams(plan, "residue_validation").get("required_residues"),
    })


# --- pocket_analysis ----------------------------------------------------------
def _pocket_analysis(wired, plan, ctx, rt) -> dict:
    fpocket = _tparams(plan, "fpocket")
    pykvfinder = _tparams(plan, "pykvfinder")
    ligand_site = _tparams(plan, "ligand_site")
    step_params = _sparams(plan)
    memory_active = (
        step_params.get("site_selection_strategy") == "detector_consensus"
        and step_params.get("pocket_memory_mode", "auto") != "off"
    )
    fpocket_params = _drop_none({
        "min_radius": fpocket.get("min_radius"),
        "max_radius": fpocket.get("max_radius"),
    })
    pykvfinder_params = _drop_none({
        "probe_in": pykvfinder.get("probe_in"),
        "probe_out": pykvfinder.get("probe_out"),
        "removal_distance": pykvfinder.get("removal_distance"),
        "volume_cutoff": pykvfinder.get("volume_cutoff"),
    })
    ligand_site_params = _drop_none({
        "residue_cutoff": ligand_site.get("residue_cutoff"),
        "box_padding": ligand_site.get("box_padding"),
        "min_box_size": ligand_site.get("min_box_size"),
    })
    return _drop_none({
        "protein_file": wired.get("provisional_receptor") or wired.get("prepared_receptor_noH"),
        "protein_noH_file": wired.get("prepared_receptor_noH"),
        "docking_box_file": (
            wired.get("docking_box")
            if step_params.get("site_selection_strategy") == "provided_box" else None
        ),
        "ligand_file": wired.get("reference_ligand"),
        "complex_mmcif_file": wired.get("source_complex_mmcif"),
        "protonation_manifest_file": wired.get("protonation_manifest"),
        "water_inventory_file": wired.get("water_inventory"),
        "metal_report_file": wired.get("metal_report"),
        "structure_qc_report_file": wired.get("structure_qc_report"),
        "chain": ctx.get("chain"),
        "known_site_residues": ctx.get("known_site_residues"),
        "known_site_provenance": ctx.get("known_site_provenance"),
        "known_required_interactions": ctx.get("known_required_interactions"),
        "site_selection_strategy": _required_sparam(
            plan, "site_selection_strategy",
        ),
        "pocket_memory_mode": step_params.get("pocket_memory_mode", "auto"),
        "pocket_memory_preset": step_params.get("pocket_memory_preset", "conservative"),
        **{
            key: step_params.get(key) if memory_active else None
            for key in (
                "pocket_memory_index", "pocket_memory_query_npz",
                "pocket_memory_query_metadata", "pocket_memory_query_case",
                "pocket_memory_query_group", "pocket_memory_exclude_cases",
                "pocket_memory_exclude_groups",
            )
        },
        "fpocket_params": fpocket_params or None,
        "pykvfinder_params": pykvfinder_params or None,
        "ligand_site_params": ligand_site_params or None,
        "consensus_top_k": step_params.get("consensus_top_k"),
        "center_match_threshold": step_params.get("center_match_threshold"),
        "residue_iou_threshold": step_params.get("residue_iou_threshold"),
        "detector_box_padding": step_params.get("detector_box_padding"),
        "native_step_timeout_seconds": step_params.get("native_step_timeout_seconds"),
        "fail_on_profiler_error": step_params.get("fail_on_profiler_error"),
        "fail_on_supplemental_profiler_error": step_params.get(
            "fail_on_supplemental_profiler_error"
        ),
        "prolif_vicinity_cutoff": _tparams(plan, "prolif").get("vicinity_cutoff"),
        "prolif_interaction_types": _tparams(plan, "prolif").get("interaction_types"),
        "luna_ifp_length": _tparams(plan, "luna").get("ifp_length"),
        "luna_ifp_type": _tparams(plan, "luna").get("ifp_type"),
        "luna_add_dependent_interactions": _tparams(plan, "luna").get(
            "add_dependent_interactions"
        ),
        "oddt_strict": _tparams(plan, "oddt").get("strict"),
        "oddt_splif_depth": _tparams(plan, "oddt").get("splif_depth"),
        "oddt_splif_size": _tparams(plan, "oddt").get("splif_size"),
        "oddt_splif_distance_cutoff": _tparams(plan, "oddt").get(
            "splif_distance_cutoff"
        ),
        "blueprint_policy": _tparams(plan, "blueprint"),
    })

# --- molecule_generation ------------------------------------------------------
def _molecule_generation(wired, plan, ctx, rt) -> dict:
    chembl = _tparams(plan, "chembl")
    dude = _tparams(plan, "dude_decoys")
    rdkit_decoys = _tparams(plan, "rdkit_decoys")
    p2m = _tparams(plan, "pocket2mol")
    p2m_runtime = _truntime(plan, "pocket2mol")
    sparams = _sparams(plan)
    # num_decoys: dude and rdkit share the concept; prefer dude's value.
    num_decoys = dude.get("num_decoys")
    if num_decoys is None:
        num_decoys = rdkit_decoys.get("num_decoys")

    payload = {
        "target_name": (
            ctx.get("target_name")
            or wired.get("inferred_target_name")
            or rt.get("inferred_target_name")
        ),
        "chembl_target_id": ctx.get("chembl_target_id"),
        "uniprot_id": ctx.get("uniprot_id"),
        "fragment_file": ctx.get("fragment_file"),
        "pocket_pdb_file": wired.get("prepared_receptor_noH"),
        "reference_ligand_file": wired.get("reference_ligand"),
        "pocket_center": wired.get("pocket_center"),
        "pocket_definition_file": wired.get("pocket_definition"),
        "interaction_blueprint_file": wired.get("interaction_blueprint"),
        "chain": ctx.get("chain"),
        "activity_type": chembl.get("activity_type"),
        "activity_threshold_nm": chembl.get("activity_threshold_nm"),
        "max_actives": chembl.get("max_actives"),
        "num_decoys": num_decoys,
        "seed": ctx.get("seed", 42),
        "neutralize": _tparams(plan, "standardization").get("neutralize"),
        "exclude_scaffolds": sparams.get("exclude_scaffolds"),
        "max_heavy_atoms": sparams.get("max_heavy_atoms"),
        "max_mw": sparams.get("max_mw"),
    }
    # pocket2mol on/off travels in the `tools` block. Still forward the pocket
    # inputs + params so the skill can run it when the switch is on.
    if _enabled(plan, "pocket2mol"):
        env_spec = R.TOOL_ENVIRONMENTS.get(
            p2m_runtime.get("environment", ""), {}
        )
        payload.update({
            "pocket2mol_num_samples": p2m.get("num_samples"),
            "pocket2mol_timeout_seconds": p2m.get("timeout_seconds"),
            "pocket2mol_gpu_device": p2m.get("gpu_device"),
            "pocket2mol_taskset_cpu": p2m.get("taskset_cpu"),
            "pocket2mol_conda_env": env_spec.get("name"),
            "pocket2mol_dir": p2m_runtime.get("source_dir"),
            "pocket2mol_checkpoint": p2m_runtime.get("checkpoint"),
        })
    return _drop_none(payload)


# --- property_prefilter / property_filter (same skill) ------------------------
def _property_common(wired, plan, ctx, molecules_sdf, rt=None) -> dict:
    # Tool on/off (pains/admet/dedup) travels in the `tools` block, not flat flags.
    scientific_facts = ctx.get("scientific_facts") or {}
    payload = {
        "molecules_sdf": molecules_sdf,
        "qed_threshold": _tparams(plan, "qed").get("qed_threshold"),
        "sa_threshold": _tparams(plan, "sa_score").get("sa_threshold"),
        "lipinski_max_violations": _tparams(plan, "lipinski").get("lipinski_max_violations"),
        "is_covalent": scientific_facts.get("is_covalent"),
    }
    # Level 2 of the gate contract: `screening_mode` names the fail count/ratio
    # that turns level-1 marks into a rejection; `candidate_source` supplies the
    # default when it is unset.
    payload.update(_screening_policy(plan, ctx))
    scoring_plan: dict[str, Any] = {}
    sources = (
        ctx,
        ctx.get("scoring_plan") if isinstance(ctx.get("scoring_plan"), dict) else {},
        plan,
        plan.get("scoring_plan")
        if isinstance(plan.get("scoring_plan"), dict) else {},
        _sparams(plan),
    )
    for source in sources:
        for field in PROPERTY_PLAN_FIELDS:
            if field in source and source[field] is not None:
                scoring_plan[field] = source[field]
    if scoring_plan:
        payload["scoring_plan"] = scoring_plan
    payload.update(_scientific_cohort(wired, plan, ctx, rt or {}))
    return payload


def _property_prefilter(wired, plan, ctx, rt) -> dict:
    sparams = _sparams(plan)
    payload = _property_common(wired, plan, ctx, wired.get("library"), rt)
    payload.update({
        "write_filtered_sdf": True,
        "filtered_verdicts": sparams.get("filtered_verdicts"),
        "filtered_sdf_name": "prefiltered_library.sdf",
        "filtered_manifest_name": "property_prefilter_manifest.csv",
        "fail_if_no_filtered_molecules": sparams.get("fail_if_no_filtered_molecules"),
    })
    return _drop_none(payload)


def _property_filter(wired, plan, ctx, rt) -> dict:
    # Prefer selected_poses; fall back to docked_poses (matches legacy behavior).
    molecules_sdf = wired.get("selected_poses")
    if not molecules_sdf or not Path(str(molecules_sdf)).exists():
        molecules_sdf = wired.get("docked_poses")
    payload = _property_common(wired, plan, ctx, molecules_sdf, rt)
    # A valid post-QC selection can contain no survivors. Never fall back to
    # rejected docking poses just because its existing SDF is empty.
    if molecules_sdf and molecules_sdf == wired.get("selected_poses"):
        payload["allow_empty_input"] = True
    return _drop_none(payload)


# --- docking ------------------------------------------------------------------
def _docking(wired, plan, ctx, rt) -> dict:
    redock = _tparams(plan, "redocking_validation")
    ligand = wired.get("prefiltered_library") or wired.get("library")
    scientific_facts = ctx.get("scientific_facts") or {}
    # Engine switches, parameters, and runtimes travel only in the `tools`
    # block appended by to_legacy_input(). Keep this adapter engine-neutral.
    payload = {
        "receptor_file": wired.get("prepared_receptor_noH"),
        "receptor_noH_file": wired.get("prepared_receptor_noH"),
        "prep_manifest_file": wired.get("prep_manifest"),
        "water_inventory_file": wired.get("water_inventory"),
        "pocket_definition_file": wired.get("pocket_definition"),
        "interaction_blueprint_file": wired.get("interaction_blueprint"),
        "ligand_file": ligand,
        "ligand_preparation_manifest_file": wired.get("library_manifest"),
        "reference_ligand_file": wired.get("reference_ligand"),
        "dock6_receptor_site_file": wired.get("dock6_receptor_site"),
        "dock6_grid_nrg_file": wired.get("dock6_grid_nrg"),
        "dock6_grid_bmp_file": wired.get("dock6_grid_bmp"),
        "is_covalent": scientific_facts.get("is_covalent"),
        "search_space_source": _required_sparam(plan, "search_space_source"),
        "redocking_rmsd_threshold": redock.get("redocking_rmsd_threshold"),
        "redocking_rmsd_hard_threshold": redock.get(
            "redocking_rmsd_hard_threshold"
        ),
    }
    payload.update(_scientific_cohort(wired, plan, ctx, rt))
    docking_box = wired.get("docking_box")
    if isinstance(docking_box, dict):
        # Inline fallback is accepted for old/external producers.  The canonical
        # Step 1 artifact is a JSON file path.
        payload["docking_box"] = docking_box
    else:
        payload["docking_box_file"] = docking_box
    return _drop_none(payload)


# --- pose_evaluation ----------------------------------------------------------
def _pose_evaluation(wired, plan, ctx, rt) -> dict:
    sparams = _sparams(plan)
    pc = _tparams(plan, "posecheck")
    scientific_facts = ctx.get("scientific_facts") or {}
    posecheck_thresholds = _drop_none({
        "strain_warn_kcal": pc.get("strain_warn_kcal"),
        "strain_fail_kcal": pc.get("strain_fail_kcal"),
        "clash_warn_count": pc.get("clash_warn_count"),
        "clash_fail_count": pc.get("clash_fail_count"),
    })
    payload = {
        "receptor_file": wired.get("prepared_receptor"),
        "docked_sdf": wired.get("docked_poses"),
        "docking_scores_file": wired.get("docking_scores"),
        "docking_report_file": wired.get("docking_report"),
        "reference_ligand_file": wired.get("reference_ligand"),
        "interaction_blueprint_file": wired.get("interaction_blueprint"),
        "pocket_features_file": wired.get("pocket_features"),
        "library_manifest_file": wired.get("library_manifest"),
        "has_metal": scientific_facts.get("has_metal"),
        "is_covalent": scientific_facts.get("is_covalent"),
        "chain": ctx.get("chain", "A"),
        "selection_strategy": sparams.get("selection_strategy"),
        "poses_per_compound": sparams.get("poses_per_compound"),
        "max_poses": sparams.get("max_poses"),
        "batch_size": sparams.get("batch_size"),
        "profiler_workers": sparams.get("profiler_workers"),
        "checkpoint_interval_seconds": sparams.get("checkpoint_interval_seconds"),
        "overall_timeout_seconds": sparams.get("overall_timeout_seconds"),
        "posecheck_num_conformers": pc.get("num_conformers"),
        "posecheck_thresholds": posecheck_thresholds or None,
        "calibrate_posecheck_against_reference": pc.get(
            "calibrate_against_reference"
        ),
        "resume": True,
    }
    payload.update(_screening_policy(plan, ctx))
    payload.update(_scientific_cohort(wired, plan, ctx, rt))
    return _drop_none(payload)


# --- coordination -------------------------------------------------------------
def _coordination(wired, plan, ctx, rt) -> dict:
    sparams = _sparams(plan)
    round_history_file = None
    upstream = rt.get("upstream_project_dir")
    if upstream:
        candidate = Path(upstream) / "step6" / "round_history.json"
        if candidate.exists():
            round_history_file = str(candidate)
    payload = {
        "evaluation_report_file": wired.get("evaluation_report"),
        "selection_manifest_file": wired.get("selection_manifest"),
        "property_filter_report_file": wired.get("property_report"),
        "docking_scores_file": wired.get("docking_scores"),
        "docking_report_file": wired.get("docking_report"),
        "pocket_analysis_report_file": wired.get("pocket_analysis_report"),
        "round_history_file": round_history_file,
        # require_quality_gates / screening_benchmark on/off travel in `tools`.
        "round_number": rt.get("round_number", 1),
        "max_rounds": sparams.get("max_rounds"),
        "pass_rate_threshold": sparams.get("pass_rate_threshold"),
        "min_candidates": sparams.get("min_candidates"),
        "min_validation_actives": sparams.get("min_validation_actives"),
        "min_validation_decoys": sparams.get("min_validation_decoys"),
        "min_validation_auroc": sparams.get("min_validation_auroc"),
        "min_validation_auroc_hard": sparams.get(
            "min_validation_auroc_hard"
        ),
        "min_validation_ef1": sparams.get("min_validation_ef1"),
        "min_validation_ef1_hard": sparams.get("min_validation_ef1_hard"),
        "normalization_method": sparams.get("normalization_method"),
        "metric_weight_overrides": sparams.get("metric_weight_overrides"),
        "category_weight_overrides": sparams.get("category_weight_overrides"),
        **_screening_policy(plan, ctx),
    }
    cohort_manifest = _merge_scientific_cohort_manifests(
        wired.get("evaluation_report"),
        wired.get("property_report"),
        wired.get("docking_report"),
    )
    if cohort_manifest:
        payload["scientific_cohorts"] = cohort_manifest
    if (
        len(cohort_manifest) == 1
        and cohort_manifest[0].get("cohort_manifest_valid") is True
    ):
        # Compatibility projection is safe only for a genuinely singular,
        # internally consistent manifest.
        payload.update({
            field: cohort_manifest[0][field]
            for field in COHORT_INPUT_FIELDS
            if _present(cohort_manifest[0].get(field))
        })
    elif not cohort_manifest:
        payload.update(_scientific_cohort(wired, plan, ctx, rt))
    return _drop_none(payload)


# --- reporting ----------------------------------------------------------------
def _reporting(wired, plan, ctx, rt) -> dict:
    sparams = _sparams(plan)
    project_dir = Path(rt["project_dir"])
    # Renderer on/off (rdkit_2d/mol_3d/molstar_main/plots/pymol_cover) travels in
    # the `tools` block, read by the skill via tool_enabled().
    payload = {
        "mode": "round",
        "project_dir": str(project_dir),
        "project_root": str(project_dir.parent),
        "round_number": rt.get("round_number", 1),
        "top_n": sparams.get("top_n"),
        "title": sparams.get("title"),
        "assets_source_dir": sparams.get("assets_source_dir"),
    }
    cohort_manifest = _scientific_cohort_manifest(
        wired.get("coordinator_report")
    )
    if cohort_manifest:
        payload["scientific_cohorts"] = cohort_manifest
    if (
        len(cohort_manifest) == 1
        and cohort_manifest[0].get("cohort_manifest_valid") is True
    ):
        # Preserve the exact completed identity emitted by the coordinator;
        # never replace it with a newly derived pending protocol token.
        singular = cohort_manifest[0]
        payload.update({
            field: singular[field]
            for field in COHORT_INPUT_FIELDS
            if _present(singular.get(field))
        })
    elif cohort_manifest:
        # A multi-cohort report has no honest run-global five-tuple.  Candidate
        # rows and the manifest retain their own partitions.
        pass
    else:
        payload.update(_scientific_cohort(wired, plan, ctx, rt))
    return _drop_none(payload)


_ADAPTERS = {
    "structure_prep": _structure_prep,
    "pocket_analysis": _pocket_analysis,
    "molecule_generation": _molecule_generation,
    "property_prefilter": _property_prefilter,
    "docking": _docking,
    "pose_evaluation": _pose_evaluation,
    "property_filter": _property_filter,
    "coordination": _coordination,
    "reporting": _reporting,
}

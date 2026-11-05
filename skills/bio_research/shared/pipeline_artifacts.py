"""Artifact wiring for the SBDD orchestrator.

Steps communicate only through *canonical artifacts* (see ARTIFACTS in
pipeline_registry). This module maps between a skill's raw result dict (its
`files` block + inline values) and those canonical names, so the orchestrator
can wire step outputs to downstream step inputs without any step knowing another
step's file naming.

Two directions:
  extract_artifacts(step_id, result)  -> {canonical: value}   (after a step runs)
  ArtifactStore                        holds {canonical: value} across the run

Physical directory layout is kept on the legacy `step{N}` scheme via
LEGACY_STEP_DIR so the report aggregator and coordinator round-history wiring keep
working unchanged during the migration.
"""
from __future__ import annotations

import hashlib
import json
import os
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable


# Clean step id -> on-disk directory name (legacy numeric scheme preserved).
LEGACY_STEP_DIR: dict[str, str] = {
    "structure_prep": "step0",
    "pocket_analysis": "step1",
    "molecule_generation": "step2",
    "property_prefilter": "step2_5",
    "docking": "step3",
    "pose_evaluation": "step4",
    "property_filter": "step5",
    "coordination": "step6",
    "reporting": "step7",
}


def _files(result: dict) -> dict:
    return dict((result or {}).get("files") or {})


def _report_member_path(value: Any, *, base: Path) -> Path | None:
    if not isinstance(value, (str, os.PathLike)) or not str(value).strip():
        return None
    try:
        path = Path(os.fspath(value)).expanduser()
        return (path if path.is_absolute() else base / path).resolve()
    except (OSError, TypeError, ValueError):
        return None


def _declared_report_members(
    result: dict,
    *,
    manifest_path: Path,
) -> dict[str, Path] | None:
    """Resolve the complete set of files the report says its HTML consumes."""
    result_files = _files(result)
    base = manifest_path.parent
    declarations = {
        "report_html": result_files.get("report_html"),
        "report_data": result_files.get("report_data_json"),
        "assets_manifest": result_files.get("assets_manifest_json"),
    }
    if result_files.get("pymol_cover_png"):
        declarations["pymol_cover"] = result_files["pymol_cover_png"]
    expected: dict[str, Path] = {}
    for role, raw_path in declarations.items():
        path = _report_member_path(raw_path, base=base)
        if path is None:
            return None
        expected[role] = path

    assets_manifest = expected["assets_manifest"]
    try:
        assets_payload = json.loads(assets_manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    assets = assets_payload.get("assets") if isinstance(assets_payload, dict) else None
    if not isinstance(assets, dict):
        return None
    licenses = assets_payload.get("licenses", {})
    if not isinstance(licenses, dict):
        return None
    for members, prefix in ((assets, "asset"), (licenses, "license")):
        for filename, record in members.items():
            if not isinstance(record, dict):
                return None
            if record.get("present") is not True:
                continue
            path = _report_member_path(record.get("path"), base=assets_manifest.parent)
            role = f"{prefix}:{filename}"
            if path is None or role in expected:
                return None
            expected[role] = path

    plots_status = (result or {}).get("plots_status") or {}
    figures = plots_status.get("figures") if isinstance(plots_status, dict) else None
    if figures is not None and not isinstance(figures, dict):
        return None
    for name, raw_path in sorted((figures or {}).items()):
        if not raw_path:
            continue
        path = _report_member_path(raw_path, base=base)
        role = f"plot:{name}"
        if path is None or role in expected:
            return None
        expected[role] = path
    return expected


def _verified_report_artifact_manifest(result: dict) -> Any:
    """Return the report bundle seal only while every member still matches."""
    result_files = _files(result)
    manifest_value = result_files.get("report_artifact_manifest_json")
    if not manifest_value:
        return None
    try:
        manifest_path = Path(os.fspath(manifest_value)).expanduser().resolve()
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        return None
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        return None

    by_role: dict[str, Path] = {}
    seen_paths: set[Path] = set()
    for record in artifacts:
        if not isinstance(record, dict):
            return None
        try:
            path = Path(os.fspath(record.get("path"))).expanduser().resolve()
            expected_size = int(record.get("size_bytes"))
            expected_sha256 = str(record.get("sha256") or "")
        except (TypeError, ValueError):
            return None
        role = str(record.get("role") or "")
        if (
            not role
            or role in by_role
            or path in seen_paths
            or not path.is_file()
            or len(expected_sha256) != 64
        ):
            return None
        try:
            int(expected_sha256, 16)
        except ValueError:
            return None
        try:
            if path.stat().st_size != expected_size:
                return None
            actual_sha256 = _sha256_file(path).removeprefix("sha256:")
        except OSError:
            return None
        if actual_sha256 != expected_sha256:
            return None
        seen_paths.add(path)
        by_role[role] = path

    expected_members = _declared_report_members(
        result,
        manifest_path=manifest_path,
    )
    if expected_members is None or set(by_role) != set(expected_members):
        return None
    for role, declared in expected_members.items():
        if by_role.get(role) != declared:
            return None
    return str(manifest_path)


def _consensus_center(result: dict) -> Any:
    docking_box = (result or {}).get("docking_box") or {}
    if isinstance(docking_box, dict) and docking_box.get("center") is not None:
        return docking_box.get("center")
    pocket_definition = (result or {}).get("pocket_definition") or {}
    if isinstance(pocket_definition, dict) and pocket_definition.get("center") is not None:
        return pocket_definition.get("center")
    detection = (result or {}).get("pocket_detection") or {}
    consensus = detection.get("consensus_pocket") or {}
    return consensus.get("center")


def _file_or_inline(result: dict, name: str) -> Any:
    """Prefer the canonical JSON artifact path, with inline producer fallback.

    Step 1 publishes both forms: files are the durable inter-step contract, while
    inline values keep result inspection and older/external producers convenient.
    """
    return _files(result).get(name) or (result or {}).get(name)


# Per-step extractors: canonical artifact -> function(result_dict) -> value|None.
# Maps each skill's native `files` key to the canonical artifact name. This
# mapping layer is intentional (skills keep descriptive local names; the pipeline
# speaks the artifact vocabulary). Emission is deterministic — one source key per
# artifact — so there are no alias `or` fallbacks except where a skill genuinely
# has two modes (e.g. structure_prep intake-only has no prepared_pdb).
ARTIFACT_EXTRACTORS: dict[str, dict[str, Callable[[dict], Any]]] = {
    "structure_prep": {
        "provisional_receptor": lambda r: _files(r).get("provisional_pdb"),
        "prep_manifest": lambda r: _files(r).get("prep_manifest"),
        "source_complex_pdb": lambda r: _files(r).get("source_complex_pdb"),
        "water_inventory": lambda r: _files(r).get("water_inventory"),
        # prep phase emits prepared_pdb; intake-only (skip_prep) emits only receptor_file.
        "inferred_target_name": lambda r: next(
            (
                chain.get("description")
                for chain in ((r.get("metadata") or {}).get("chains") or [])
                if isinstance(chain, dict) and chain.get("description")
                and (r.get("metadata") or {}).get("selected_chain") in
                    (chain.get("author_chain_ids") or chain.get("chain_ids") or [])
            ),
            None,
        ),
        "prepared_receptor_noH": lambda r: _files(r).get("prepared_noH_pdb"),
        "source_complex_mmcif": lambda r: _files(r).get("complex_mmcif"),
        "reference_ligand": lambda r: _files(r).get("ligand_sdf"),
        "protonation_manifest": lambda r: _files(r).get("protonation_manifest"),
        "metal_report": lambda r: _files(r).get("metal_report"),
        "structure_qc_report": lambda r: _files(r).get("structure_qc_report"),
    },
    "pocket_analysis": {
        "pocket_center": _consensus_center,
        "pocket_definition": lambda r: _file_or_inline(r, "pocket_definition"),
        "docking_box": lambda r: (
            None if (r.get("parameters") or {}).get("site_selection_strategy") == "provided_box"
            else _file_or_inline(r, "docking_box")
        ),
        "pocket_features": lambda r: _files(r).get("pocket_features"),
        "interaction_observations": lambda r: _files(r).get("interaction_observations"),
        "interaction_blueprint": lambda r: _files(r).get("blueprint"),
        "pocket_analysis_report": lambda r: _files(r).get("report"),
        "deeppocket_predictions": lambda r: _files(r).get("deeppocket"),
    },
    "molecule_generation": {
        "library": lambda r: _files(r).get("library_sdf"),
        "library_manifest": lambda r: _files(r).get("library_manifest"),
    },
    "property_prefilter": {
        "prefiltered_library": lambda r: _files(r).get("prefiltered_library"),
    },
    "docking": {
        "prepared_receptor": lambda r: _files(r).get("prepared_receptor"),
        "receptor_pdbqt": lambda r: _files(r).get("receptor_pdbqt"),
        "receptor_state_manifest": lambda r: _files(r).get("receptor_state_manifest"),
        "water_report": lambda r: _files(r).get("water_report"),
        "docking_species_library": lambda r: _files(r).get("docking_species_library"),
        "microstate_identity_map": lambda r: _files(r).get("microstate_identity_map"),
        "docking_prep_qc": lambda r: _files(r).get("docking_prep_qc"),
        "reference_ligand_microstates": lambda r: _files(r).get("reference_ligand_microstates"),
        "docked_poses": lambda r: _files(r).get("docked_poses_sdf"),
        "docking_scores": lambda r: _files(r).get("docking_scores_json"),
        "docking_report": lambda r: _files(r).get("docking_report_json"),
    },
    "pose_evaluation": {
        "selected_poses": lambda r: _files(r).get("selected_poses_sdf"),
        "selection_manifest": lambda r: _files(r).get(
            "selection_manifest_csv"
        ),
        "evaluation_report": lambda r: _files(r).get("evaluation_report_json"),
        "candidate_topologies": lambda r: _files(r).get(
            "candidate_topologies_manifest"
        ),
        "rtmscore_manifest": lambda r: _files(r).get(
            "rtmscore_manifest"
        ),
        "short_md_manifest": lambda r: _files(r).get(
            "short_md_manifest"
        ),
        "trajectory_stability_report": lambda r: _files(r).get(
            "trajectory_stability_report"
        ),
        "mmgbsa_manifest": lambda r: _files(r).get(
            "mmgbsa_manifest"
        ),
        "openfe_rbfe_manifest": lambda r: _files(r).get(
            "openfe_rbfe_manifest"
        ),
        "ties_rbfe_manifest": lambda r: _files(r).get(
            "ties_rbfe_manifest"
        ),
    },
    "property_filter": {
        "property_report": lambda r: _files(r).get("property_filter_report_json"),
        "synthesis_report": lambda r: _files(r).get("synthesis_report_json"),
    },
    "coordination": {
        "coordinator_report": lambda r: _files(r).get("coordinator_report_json"),
    },
    "reporting": {
        "report_html": lambda r: _files(r).get("report_html"),
        "report_data": lambda r: _files(r).get("report_data_json"),
        "report_assets_manifest": lambda r: _files(r).get(
            "assets_manifest_json"
        ),
        "report_artifact_manifest": _verified_report_artifact_manifest,
        "report_pymol_cover": lambda r: _files(r).get("pymol_cover_png"),
    },
}


def extract_artifacts(step_id: str, result: dict) -> dict[str, Any]:
    """Pull all canonical artifacts a step produced from its raw result dict.

    Only non-empty values are returned. Values may be file paths (most artifacts)
    or inline values (pocket_center).
    """
    # A failed step may retain proposed/diagnostic values and stale file paths,
    # but none of them are valid canonical inputs for downstream steps.
    if (result or {}).get("success") is False:
        return {}
    out: dict[str, Any] = {}
    for name, fn in ARTIFACT_EXTRACTORS.get(step_id, {}).items():
        try:
            value = fn(result or {})
        except Exception:
            value = None
        if value not in (None, ""):
            out[name] = value
    return out


_INLINE_ARTIFACTS = {"pocket_center", "inferred_target_name"}
_HASH_CHUNK_SIZE = 1024 * 1024


def _is_empty(value: Any) -> bool:
    return value is None or (isinstance(value, str) and value == "")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK_SIZE), b""):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def _inline_fingerprint(value: Any) -> str:
    try:
        payload = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            default=str,
        ).encode("utf-8")
    except (TypeError, ValueError):
        payload = repr(value).encode("utf-8")
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def _step_specs() -> dict[str, dict]:
    """Load the registry without constraining how this module is imported."""
    try:
        from .pipeline_registry import STEPS
    except ImportError:
        from pipeline_registry import STEPS
    return STEPS


def _artifact_dependency_graph() -> dict[str, set[str]]:
    """Return upstream artifact -> directly-derived artifact edges.

    A step output conservatively depends on every artifact the step may consume.
    This is deliberately stricter than the runtime's ``requires_any`` selection:
    retaining an output produced with an older optional input is unsafe.
    """
    graph: dict[str, set[str]] = {}
    for spec in _step_specs().values():
        for consumed in spec.get("consumes", []):
            graph.setdefault(consumed, set()).update(spec.get("produces", []))
    return graph


def _declared_producer(artifact: str) -> str:
    for step_id, extractors in ARTIFACT_EXTRACTORS.items():
        if artifact in extractors:
            return step_id
    return f"initial:{artifact}"


class ArtifactStore:
    """Accumulate canonical artifacts while enforcing cross-round lineage.

    Replacing an upstream artifact invalidates every artifact derived from it.
    Re-running a producer also removes that producer's prior outputs when they
    are absent from the new result. Active provenance and an append-only event
    history make those decisions auditable.
    """

    def __init__(self, initial: dict[str, Any] | None = None):
        self._store: dict[str, Any] = {}
        self._provenance: dict[str, dict[str, Any]] = {}
        self._history: list[dict[str, Any]] = []
        self._event_sequence = 0
        self._static_dependencies = _artifact_dependency_graph()
        self._dynamic_dependencies: dict[str, set[str]] = {}
        self._step_producers = set(ARTIFACT_EXTRACTORS)

        for name, value in (initial or {}).items():
            if _is_empty(value):
                continue
            metadata = self._inspect(name, value)
            if metadata is None:
                self._record_event(
                    "reject",
                    name,
                    reason="initial path is missing or unreadable",
                )
                continue
            self._set(
                name,
                value,
                producer=_declared_producer(name),
                source="initial",
                metadata=metadata,
            )

    def update_from(
        self,
        step_id: str,
        result: dict,
        *,
        source: str | None = None,
        fingerprints: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """Publish one producer result and invalidate incompatible old lineage.

        ``source`` and ``fingerprints`` are optional so existing callers remain
        valid. Fingerprints are normally computed from file bytes/inline JSON.
        """
        self._refresh_all()
        extracted = extract_artifacts(step_id, result)
        source_name = source or f"step_result:{step_id}"
        accepted: dict[str, tuple[Any, dict[str, Any]]] = {}
        rejected: set[str] = set()

        for name, value in extracted.items():
            metadata = self._inspect(
                name,
                value,
                fingerprint=(fingerprints or {}).get(name),
            )
            if metadata is None:
                rejected.add(name)
                self._record_event(
                    "reject",
                    name,
                    reason="producer path is missing or unreadable",
                    producer=step_id,
                    source=source_name,
                )
            else:
                accepted[name] = (value, metadata)

        incoming = set(accepted)
        previous = {
            name
            for name, provenance in self._provenance.items()
            if provenance.get("producer") == step_id
        }

        # A new successful result is the complete current output set for this
        # producer. Optional outputs omitted now must not survive from a prior
        # round. Failed results extract no outputs and therefore clear them all.
        stale = previous - incoming
        if stale:
            self._invalidate_many(
                stale,
                reason=f"{step_id} no longer produced this output",
            )

        changed: set[str] = set(rejected)
        for name, (_, metadata) in accepted.items():
            old = self._provenance.get(name)
            if (
                old is None
                or old.get("fingerprint") != metadata["fingerprint"]
                or old.get("producer") != step_id
            ):
                changed.add(name)
        if changed:
            self._invalidate_many(
                changed,
                reason=f"upstream output replaced by {step_id}",
                exclude=incoming,
            )

        recorded_inputs = (result or {}).get("artifact_inputs")
        input_names = (
            set(recorded_inputs)
            if isinstance(recorded_inputs, dict)
            else set()
        )
        self._replace_dynamic_dependencies(previous | incoming, set())
        self._replace_dynamic_dependencies(incoming, input_names)
        for name, (value, metadata) in accepted.items():
            self._set(
                name,
                value,
                producer=step_id,
                source=source_name,
                metadata=metadata,
            )
        return {name: self._store[name] for name in accepted}

    def preload_from(
        self,
        step_id: str,
        result: dict,
        *,
        source: str,
        expected_producer_fingerprint: str | None = None,
        expected_optional_inputs: set[str] | None = None,
    ) -> dict[str, Any]:
        """Publish a cached result only when its recorded inputs still match.

        Producer steps without artifact inputs (for example
        ``structure_prep``) may record an empty input set, but every cached
        result must carry output fingerprints. When an expected producer-config
        fingerprint is supplied, it must match as well.
        """
        reasons: list[str] = []
        if expected_producer_fingerprint is not None:
            recorded = (result or {}).get("producer_config_fingerprint")
            if recorded != expected_producer_fingerprint:
                reasons.append(
                    "producer config fingerprint is missing or changed "
                    f"({recorded!r} != {expected_producer_fingerprint!r})"
                )
        inputs_match, input_reasons = self.verify_result_inputs(
            step_id,
            result,
            expected_optional_inputs=expected_optional_inputs,
        )
        outputs_match, output_reasons = self.verify_result_outputs(step_id, result)
        reasons.extend(input_reasons)
        reasons.extend(output_reasons)
        if reasons or not inputs_match or not outputs_match:
            self.reject_preload(
                step_id,
                source=source,
                reason="; ".join(reasons),
            )
            return {}
        return self.update_from(step_id, result, source=source)

    def verify_result_inputs(
        self,
        step_id: str,
        result: dict,
        *,
        expected_optional_inputs: set[str] | None = None,
    ) -> tuple[bool, list[str]]:
        """Compare cached input lineage with the active artifact generation."""
        self._refresh_all()
        spec = _step_specs().get(step_id, {})
        consumed = set(spec.get("consumes", []))
        optional = set(spec.get("optional_consumes", []))
        if not consumed and not optional:
            return True, []

        recorded = (result or {}).get("artifact_inputs")
        if not isinstance(recorded, dict):
            return False, [
                "cached derived result has no artifact_inputs lineage metadata"
            ]

        recorded_names = set(recorded)
        allowed_names = consumed | optional
        relevant_optional = (
            optional & recorded_names
            if expected_optional_inputs is None
            else optional & set(expected_optional_inputs)
        )
        active_names = (consumed & set(self._store)) | (
            relevant_optional & set(self._store)
        )
        reasons: list[str] = []
        if expected_optional_inputs is not None:
            undeclared_expected = sorted(
                set(expected_optional_inputs) - optional
            )
            if undeclared_expected:
                reasons.append(
                    "orchestrator expected undeclared optional inputs: "
                    f"{undeclared_expected}"
                )
        unexpected = sorted(recorded_names - allowed_names)
        if unexpected:
            reasons.append(
                f"cached result records undeclared artifact inputs: {unexpected}"
            )
        if recorded_names != active_names:
            missing = sorted(recorded_names - active_names)
            added = sorted(active_names - recorded_names)
            if missing:
                reasons.append(
                    f"recorded inputs are no longer active: {missing}"
                )
            if added:
                reasons.append(
                    f"active inputs were absent from the cached result: {added}"
                )

        for name in sorted(recorded_names & active_names):
            expected = recorded.get(name)
            expected = expected if isinstance(expected, dict) else {}
            expected_fingerprint = expected.get("fingerprint")
            active = self._provenance.get(name) or {}
            active_fingerprint = active.get("fingerprint")
            if (
                active.get("kind") == "path"
                and (
                    not active.get("fingerprint_supplied")
                    or str(active_fingerprint).startswith("sha256:")
                )
            ):
                try:
                    current_fingerprint = _sha256_file(Path(active["path"]))
                except (OSError, ValueError, KeyError):
                    self._invalidate_path_generation(
                        name,
                        reason="canonical artifact path became unreadable "
                        "during input revalidation",
                    )
                    reasons.append(
                        f"input {name!r} became unreadable during revalidation"
                    )
                    continue
                if current_fingerprint != active_fingerprint:
                    self._invalidate_path_generation(
                        name,
                        reason="canonical artifact file changed during step "
                        "execution or cache validation",
                    )
                    reasons.append(
                        f"input {name!r} fingerprint changed "
                        f"({expected_fingerprint} != {current_fingerprint})"
                    )
                    continue
            if not expected_fingerprint:
                reasons.append(
                    f"recorded input {name!r} has no content fingerprint"
                )
            elif expected_fingerprint != active_fingerprint:
                reasons.append(
                    f"input {name!r} fingerprint changed "
                    f"({expected_fingerprint} != {active_fingerprint})"
                )
            expected_kind = expected.get("kind")
            active_kind = active.get("kind")
            if expected_kind and expected_kind != active_kind:
                reasons.append(
                    f"input {name!r} kind changed "
                    f"({expected_kind} != {active_kind})"
                )
        return not reasons, reasons

    def verify_result_outputs(
        self,
        step_id: str,
        result: dict,
    ) -> tuple[bool, list[str]]:
        """Check cached outputs against fingerprints saved at publication."""
        extracted = extract_artifacts(step_id, result)
        recorded = (result or {}).get("artifact_outputs")
        if not isinstance(recorded, dict):
            return False, [
                "cached result has no artifact_outputs lineage metadata"
            ]

        extracted_names = set(extracted)
        recorded_names = set(recorded)
        reasons: list[str] = []
        if extracted_names != recorded_names:
            missing = sorted(recorded_names - extracted_names)
            added = sorted(extracted_names - recorded_names)
            if missing:
                reasons.append(
                    f"recorded outputs are no longer present: {missing}"
                )
            if added:
                reasons.append(
                    f"unrecorded outputs appeared in cached result: {added}"
                )

        for name in sorted(extracted_names & recorded_names):
            expected = recorded.get(name)
            expected = expected if isinstance(expected, dict) else {}
            current = self._inspect(name, extracted[name])
            if current is None:
                reasons.append(
                    f"cached output {name!r} is missing or unreadable"
                )
                continue
            expected_fingerprint = expected.get("fingerprint")
            if not expected_fingerprint:
                reasons.append(
                    f"recorded output {name!r} has no content fingerprint"
                )
            elif expected_fingerprint != current.get("fingerprint"):
                reasons.append(
                    f"output {name!r} fingerprint changed "
                    f"({expected_fingerprint} != {current.get('fingerprint')})"
                )
            expected_kind = expected.get("kind")
            if expected_kind and expected_kind != current.get("kind"):
                reasons.append(
                    f"output {name!r} kind changed "
                    f"({expected_kind} != {current.get('kind')})"
                )
        return not reasons, reasons

    def reject_preload(
        self,
        step_id: str,
        *,
        source: str,
        reason: str,
    ) -> set[str]:
        """Reject a newer cached result and tombstone older producer output."""
        self._record_event(
            "reject_result",
            f"step:{step_id}",
            reason=reason,
            producer=step_id,
            source=source,
        )
        return self.invalidate_producer(
            step_id,
            reason=f"newer cached result rejected: {reason}",
        )

    def invalidate_producer(
        self,
        step_id: str,
        *,
        reason: str,
    ) -> set[str]:
        """Invalidate every active output from one producer and derivatives."""
        roots = {
            name
            for name, provenance in self._provenance.items()
            if provenance.get("producer") == step_id
        }
        if not roots:
            return set()
        return self._invalidate_many(roots, reason=reason)

    def inject(
        self,
        artifacts: dict[str, Any],
        *,
        producer: str = "input_artifacts",
        source: str = "config.input_artifacts",
        fingerprints: dict[str, str] | None = None,
    ) -> None:
        """Seed externally-supplied artifacts (config `input_artifacts`).

        Used to feed a step without running its producer — e.g. a modified
        library SDF for an iteration round. Injection is atomic: explicitly
        supplied companion artifacts are retained, while all other outputs from
        the replaced producer bundle and every downstream derivative are
        invalidated.
        """
        self._refresh_all()
        requested = {
            name: value
            for name, value in (artifacts or {}).items()
            if not _is_empty(value)
        }
        accepted: dict[str, tuple[Any, dict[str, Any]]] = {}
        for name, value in requested.items():
            metadata = self._inspect(
                name,
                value,
                fingerprint=(fingerprints or {}).get(name),
            )
            if metadata is None:
                self._record_event(
                    "reject",
                    name,
                    reason="injected path is missing or unreadable",
                    producer=producer,
                    source=source,
                )
            else:
                accepted[name] = (value, metadata)

        incoming = set(accepted)
        invalidation_roots = set(requested)

        # Replacing one member of a step-produced bundle (for example the
        # receptor) invalidates its old optional siblings (PDBQT, ligand, box
        # inputs derived later, etc.). Callers can preserve a compatible sibling
        # by injecting it in the same atomic call.
        for name in requested:
            old_producer = (self._provenance.get(name) or {}).get("producer")
            if old_producer not in self._step_producers:
                continue
            invalidation_roots.update(
                artifact
                for artifact, provenance in self._provenance.items()
                if provenance.get("producer") == old_producer
            )

        if invalidation_roots:
            self._invalidate_many(
                invalidation_roots,
                reason="external artifact injection replaced prior lineage",
                exclude=incoming,
            )

        self._replace_dynamic_dependencies(incoming, set())
        for name, (value, metadata) in accepted.items():
            self._set(
                name,
                value,
                producer=producer,
                source=source,
                metadata=metadata,
            )

    def get(self, name: str, default: Any = None) -> Any:
        self._refresh(name)
        value = self._store.get(name, default)
        if (self._provenance.get(name) or {}).get("kind") == "inline":
            return deepcopy(value)
        return value

    def has(self, name: str) -> bool:
        self._refresh(name)
        return name in self._store

    def names(self) -> set[str]:
        self._refresh_all()
        return set(self._store)

    def as_dict(self) -> dict[str, Any]:
        self._refresh_all()
        return {
            name: (
                deepcopy(value)
                if (self._provenance.get(name) or {}).get("kind") == "inline"
                else value
            )
            for name, value in self._store.items()
        }

    def provenance(self, name: str | None = None) -> dict[str, Any]:
        """Return active provenance for one artifact or the whole store."""
        if name is not None:
            self._refresh(name)
            return deepcopy(self._provenance.get(name, {}))
        self._refresh_all()
        return deepcopy(self._provenance)

    def history(self) -> list[dict[str, Any]]:
        """Return append-only publication/invalidation audit events."""
        self._refresh_all()
        return deepcopy(self._history)

    def invalidate(self, name: str, *, reason: str = "manual invalidation") -> set[str]:
        """Explicitly invalidate an artifact and all of its derivatives."""
        return self._invalidate_many({name}, reason=reason)

    def downstream_of(self, name: str) -> set[str]:
        """Return the transitive declared derivative set for inspection/tests."""
        return self._closure({name}) - {name}

    def _inspect(
        self,
        name: str,
        value: Any,
        *,
        fingerprint: str | None = None,
    ) -> dict[str, Any] | None:
        if name in _INLINE_ARTIFACTS or not isinstance(value, (str, os.PathLike)):
            return {
                "kind": "inline",
                "fingerprint": fingerprint or _inline_fingerprint(value),
                "fingerprint_supplied": fingerprint is not None,
            }

        try:
            path = Path(os.fspath(value)).expanduser().resolve()
            if not path.is_file():
                return None
            stat = path.stat()
            computed = fingerprint or _sha256_file(path)
        except (OSError, ValueError):
            return None
        return {
            "kind": "path",
            "path": str(path),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "fingerprint": computed,
            "fingerprint_supplied": fingerprint is not None,
        }

    def _set(
        self,
        name: str,
        value: Any,
        *,
        producer: str,
        source: str,
        metadata: dict[str, Any],
    ) -> None:
        previous = deepcopy(self._provenance.get(name))
        if metadata.get("kind") == "inline":
            stored_value = deepcopy(value)
        else:
            stored_value = metadata.get("path", value)
        self._store[name] = stored_value
        self._event_sequence += 1
        provenance = {
            "producer": producer,
            "source": source,
            **metadata,
            "sequence": self._event_sequence,
        }
        self._provenance[name] = provenance
        self._record_event(
            "replace" if previous else "publish",
            name,
            producer=producer,
            source=source,
            provenance=provenance,
            previous=previous,
        )

    def _refresh(self, name: str) -> None:
        provenance = self._provenance.get(name)
        if not provenance or provenance.get("kind") != "path":
            return
        try:
            path = Path(provenance["path"])
            if not path.is_file():
                raise FileNotFoundError(path)
            stat = path.stat()
        except (OSError, ValueError):
            self._invalidate_path_generation(
                name,
                reason="canonical artifact path disappeared",
            )
            return

        if (
            stat.st_size == provenance.get("size")
            and stat.st_mtime_ns == provenance.get("mtime_ns")
        ):
            return

        # A caller-supplied digest may use an algorithm unavailable here. In
        # that case any in-place file mutation conservatively breaks lineage.
        if provenance.get("fingerprint_supplied"):
            self._invalidate_path_generation(
                name,
                reason="canonical artifact file changed in place",
            )
            return
        try:
            current_fingerprint = _sha256_file(path)
        except OSError:
            self._invalidate_path_generation(
                name,
                reason="canonical artifact path became unreadable",
            )
            return
        if current_fingerprint != provenance.get("fingerprint"):
            self._invalidate_path_generation(
                name,
                reason="canonical artifact file changed in place",
            )
            return
        provenance["size"] = stat.st_size
        provenance["mtime_ns"] = stat.st_mtime_ns

    def _invalidate_path_generation(self, name: str, *, reason: str) -> None:
        """Invalidate a mutated file, its producer bundle, and derivatives."""
        roots = {name}
        producer = (self._provenance.get(name) or {}).get("producer")
        if producer in self._step_producers:
            roots.update(
                artifact
                for artifact, provenance in self._provenance.items()
                if provenance.get("producer") == producer
            )
        self._invalidate_many(roots, reason=reason)

    def _refresh_all(self) -> None:
        for name in list(self._store):
            self._refresh(name)

    def _closure(self, roots: set[str]) -> set[str]:
        closure = set(roots)
        pending = list(roots)
        while pending:
            current = pending.pop()
            derived_artifacts = (
                self._static_dependencies.get(current, set())
                | self._dynamic_dependencies.get(current, set())
            )
            for derived in derived_artifacts:
                if derived not in closure:
                    closure.add(derived)
                    pending.append(derived)
        return closure

    def _replace_dynamic_dependencies(
        self,
        outputs: set[str],
        inputs: set[str],
    ) -> None:
        """Bind this output generation to the inputs it actually consumed."""
        if not outputs:
            return
        for name in list(self._dynamic_dependencies):
            self._dynamic_dependencies[name].difference_update(outputs)
            if not self._dynamic_dependencies[name]:
                self._dynamic_dependencies.pop(name, None)
        for name in inputs:
            self._dynamic_dependencies.setdefault(name, set()).update(outputs)

    def _invalidate_many(
        self,
        roots: set[str],
        *,
        reason: str,
        exclude: set[str] | None = None,
    ) -> set[str]:
        root_names = set(roots)
        closure = self._closure(root_names)
        # Explicit external seeds are lineage boundaries.  The static registry
        # graph describes what a step-produced artifact *can* depend on, so a
        # fresh upstream producer may reach an injected artifact that is
        # independent in this generation (for example, Step 1 reaches the
        # fixed input library through the optional Step 2 route).  Preserve the
        # seed itself unless it is a direct invalidation root; derivatives past
        # that boundary remain in ``closure`` and are still invalidated.
        protected_external_seeds = {
            name
            for name in closure - root_names
            if (self._provenance.get(name) or {}).get("producer")
            in {"input_artifacts", "context"}
        }
        invalidated = (
            closure
            - protected_external_seeds
            - set(exclude or ())
        ) & set(self._store)
        for name in sorted(invalidated):
            previous = self._provenance.pop(name, {})
            self._store.pop(name, None)
            self._record_event(
                "invalidate",
                name,
                reason=reason,
                previous=previous,
            )
        return invalidated

    def _record_event(
        self,
        action: str,
        artifact: str,
        *,
        reason: str | None = None,
        producer: str | None = None,
        source: str | None = None,
        provenance: dict[str, Any] | None = None,
        previous: dict[str, Any] | None = None,
    ) -> None:
        event: dict[str, Any] = {
            "event": len(self._history) + 1,
            "action": action,
            "artifact": artifact,
        }
        if reason:
            event["reason"] = reason
        if producer:
            event["producer"] = producer
        if source:
            event["source"] = source
        if provenance:
            event["provenance"] = deepcopy(provenance)
        if previous:
            event["previous"] = deepcopy(previous)
        self._history.append(event)

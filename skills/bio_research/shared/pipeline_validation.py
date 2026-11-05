"""SBDD pipeline config validation — declarative execution hard stops.

All checks run once, before execution, after the explicit selection matrix has
been verified and registered parameter defaults have been materialized. Any
violation raises ConfigValidationError with a collected list of messages.

Rules:
  VR1  unknown tool id in a step
  VR2  unknown param key (step-level, tool-level, context, or project)
  VR3  param type / range violation (params, context, project)
  VR4  step/tool DAG: an enabled step or tool's required artifact has no enabled
       producer (input_artifacts and upstream preloads count as satisfied)
  VR5  same-step tool dependency (e.g. blueprint_match needs PLIP)
  VR6  class "at least one" (e.g. >=1 pocket detector, >=1 actives source)
  VR7  tool context prerequisite (e.g. chembl needs a target id)
  VR8  exactly-one tool groups (e.g. one primary docking engine)
  VR9  an enabled tool cannot be hidden under a disabled parent step
  VR10 enabled MM-GBSA requires a built, sufficiently large candidate-topology
       cohort and cannot be silently removed with Step 4

  Also: unknown top-level keys and unknown/invalid input_artifacts hard-fail.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import pipeline_registry as R
import pipeline_config as C
import pipeline_timeouts as T
from constants import INTERACTION_TYPE_CANONICAL, normalize_interaction_type
from pose_identity import benchmark_label


# Agent-facing summaries of the hard checks implemented below.  The catalog
# generator imports this mapping when it renders TOOL_SELECTION_POLICY.md, so
# the policy points at the validator's contracts instead of maintaining a
# second, hand-written list of safety rules.
VALIDATION_RULES = {
    "selection_matrix": (
        "before validation, every registered step and tool must have an explicit "
        "enabled boolean; registry defaults never complete a missing selection"
    ),
    "preparation_phases": (
        "structure_prep is geometry-only (PDBFixer add_hydrogens=false; missing residues opt-in); "
        "PROPKA is diagnostic and Reduce provisional. Final protonation, waters and PDBQT "
        "belong to mandatory docking.docking_prep after prefilter; explicit chain, backend, pH "
        "when applicable, water policy and ligand-state justification are required."
    ),
    "VR1": "every tool id under a step must exist in that step's registry",
    "VR2": "every param/context/project key must exist in its schema",
    "VR3": (
        "param type / min-max / choices (including list item choices) must hold; "
        "configurable metric soft bounds must not exceed their hard bounds"
    ),
    "VR4": (
        "an enabled step/tool's required artifact must be produced by an enabled "
        "upstream step, preloaded from upstream_project_dir, or supplied via "
        "input_artifacts"
    ),
    "VR5": "same-step tool dependencies and tool-one-of requirements must hold",
    "VR6": "every at-least-one tool class must retain an enabled member",
    "VR7": "tool context prerequisites must be present",
    "VR8": "exactly one tool must be enabled in each mutually-exclusive class",
    "VR9": (
        "an explicitly enabled tool may not remain under a parent step disabled "
        "by explicit config; range-excluded upstream settings may remain for lineage"
    ),
    "VR10": (
        "enabled MM-GBSA requires pose_evaluation plus a built, sufficiently "
        "large candidate-topology cohort"
    ),
    "top_level": (
        "unknown top-level keys hard-fail; input_artifacts must be injectable file "
        "paths and may not collide with an enabled producer"
    ),
    "explicit_modes": (
        "pocket site-selection and docking search-space modes must be explicit, "
        "compatible with the selected tools, and backed by required artifacts"
    ),
    "capability": (
        "unsupported covalent docking, generic Amber topology for declared metal "
        "sites, and tools missing explicitly supplied artifacts hard-fail"
    ),
    "evaluation": (
        "pose identity, generated-3D QC, short-MD, MM-GBSA, and RBFE escalation "
        "contracts must remain complete; first_passing_per_compound cannot be "
        "combined with mdanalysis or mmgbsa because they revise verdicts after selection"
    ),
    "metric_gates": (
        "configurable soft/hard metric thresholds must preserve "
        "pass/flag/fail marking order"
    ),
}


class ConfigValidationError(ValueError):
    """Raised when a merged config violates one or more VR rules."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        body = "\n".join(f"  - {e}" for e in errors)
        super().__init__(f"pipeline config invalid ({len(errors)} error(s)):\n{body}")


def validate(merged_config: dict, available_artifacts: set[str] | None = None) -> None:
    """Validate a merged config. Raises ConfigValidationError if any rule fails.

    Args:
        merged_config: explicit selection with registered parameter defaults filled.
        available_artifacts: artifacts already present at pipeline entry (e.g. loaded
            from upstream_project_dir on a resume/round-N run). Consumed artifacts in
            this set are treated as satisfied even if their producing step is disabled,
            matching the orchestrator's _load_existing_result behavior. Empty/None ⇒
            strict full-pipeline validation.
    """
    errors: list[str] = []
    errors.extend(C.explicit_selection_errors(merged_config))
    steps_cfg = merged_config.get("steps") or {}
    context = merged_config.get("context") or {}
    available = set(available_artifacts or set())

    _check_top_level_keys(merged_config, errors)                     # unknown top-level
    _check_step_timeouts(merged_config.get("step_timeouts"), errors)
    _check_runtime_resources(merged_config.get("runtime_resources"), errors)
    _check_mapping_schema("context", context, R.CONTEXT_SCHEMA, errors)   # VR2/VR3
    _check_context_artifact_inputs(
        context, merged_config.get("input_artifacts") or {}, errors,
    )
    _check_known_site_provenance(context, errors)
    _check_known_required_interactions(context, errors)
    _check_mapping_schema("project", merged_config.get("project") or {}, R.PROJECT_SCHEMA, errors)
    _check_scientific_facts(context, errors)
    _check_input_artifacts(
        merged_config.get("input_artifacts") or {}, steps_cfg, errors,
    )
    _check_benchmark_labels(
        merged_config.get("input_artifacts") or {}, steps_cfg, errors,
    )

    # Unknown step ids (defensive; merge normally only introduces known steps).
    for step_id in steps_cfg:
        if step_id not in R.STEPS:
            errors.append(f"unknown step {step_id!r}")

    for step_id, step_spec in R.STEPS.items():
        cfg_step = steps_cfg.get(step_id, {})
        if not isinstance(cfg_step, dict):
            errors.append(f"step {step_id!r} config must be an object")
            continue

        _check_step_params(step_id, step_spec, cfg_step, errors)      # VR2/VR3
        _check_tools(step_id, step_spec, cfg_step, errors, context)   # VR1/VR2/VR3/VR7

    _check_timeout_budgets(merged_config, steps_cfg, errors)

    # Artifacts supplied via input_artifacts are available to the DAG.
    available |= set((merged_config.get("input_artifacts") or {}).keys())
    available |= _context_artifacts(context)
    _check_dag(steps_cfg, available, errors)                          # VR4
    _check_tool_dependencies(steps_cfg, errors)                       # VR5
    _check_one_of(steps_cfg, errors)                                  # VR6
    _check_exactly_one(steps_cfg, errors)                              # VR8
    _check_disabled_step_opt_ins(merged_config, steps_cfg, errors)
    _check_structure_prep_contract(steps_cfg, context, errors)
    if _step_enabled(steps_cfg, "docking"):
        from receptor_preparation import docking_prep_errors
        errors.extend(docking_prep_errors(steps_cfg["docking"].get("tools") or {}))
    if _step_enabled(steps_cfg, "pocket_analysis") and not str(context.get("chain") or "").strip():
        errors.append("context.chain must explicitly identify the target chain after entity/assembly review; never assume A")
    _check_blueprint_evidence_policy(steps_cfg, context, errors)
    _check_scientific_capabilities(context, steps_cfg, errors)
    _check_explicit_modes(steps_cfg, context, available, errors)
    _check_provided_box_chain(merged_config, errors)
    _check_molecule_generation_contract(steps_cfg, errors)
    _check_mmgbsa_contract(steps_cfg, errors)
    _check_short_md_contract(steps_cfg, errors)
    _check_pose_identity_contract(steps_cfg, errors)
    _check_pose_selection_contract(steps_cfg, errors)
    _check_generated_pose_qc(steps_cfg, errors)
    _check_metric_gate_contract(steps_cfg, errors)

    if errors:
        raise ConfigValidationError(errors)


def _context_artifacts(context: dict) -> set[str]:
    """Return canonical artifacts explicitly supplied by context file fields."""
    return {
        artifact
        for artifact, context_key in R.CONTEXT_ARTIFACTS.items()
        if _readable_file(context.get(context_key))
    }


def _check_provided_box_chain(config: dict, errors: list[str]) -> None:
    """Reject a supplied receptor/box/chain mismatch before executing tools."""
    steps = config.get("steps") or {}
    pocket = steps.get("pocket_analysis") or {}
    if not _step_enabled(steps, "pocket_analysis") or (pocket.get("params") or {}).get("site_selection_strategy") != "provided_box":
        return
    artifacts = config.get("input_artifacts") or {}
    protein, box = artifacts.get("prepared_receptor_noH"), artifacts.get("docking_box")
    if not protein or not box:
        return
    sys.path.insert(0, str(_HERE.parent / "pocket-analysis" / "scripts"))
    from ligand_site_utils import build_box_site
    try:
        payload = json.loads(Path(box).expanduser().read_text(encoding="utf-8"))
        build_box_site(Path(protein).expanduser(), payload, (config.get("context") or {}).get("chain"))
    except (OSError, ValueError, TypeError, KeyError) as exc:
        errors.append(f"provided_box receptor/chain validation failed: {exc}")


def _readable_file(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    try:
        with Path(value).expanduser().open("rb") as handle:
            handle.read(1)
    except (OSError, ValueError):
        return False
    return True


def _check_context_artifact_inputs(
    context: dict,
    injected: dict,
    errors: list[str],
) -> None:
    """Validate explicit context artifacts and reject ambiguous dual sources."""
    injected = injected if isinstance(injected, dict) else {}
    for artifact, context_key in R.CONTEXT_ARTIFACTS.items():
        value = context.get(context_key)
        if value in (None, ""):
            continue
        if not _readable_file(value):
            errors.append(
                f"context.{context_key} file does not exist or is not readable: "
                f"{value}"
            )
        if artifact in injected:
            errors.append(
                f"context.{context_key} and input_artifacts.{artifact} both "
                "supply the same canonical artifact; choose exactly one source"
            )


# --- top-level / context / project / input_artifacts -------------------------
def _check_top_level_keys(merged_config: dict, errors: list[str]) -> None:
    for key in merged_config:
        if key not in R.TOP_LEVEL_KEYS:
            errors.append(
                f"unknown top-level config key {key!r} "
                f"(allowed: {', '.join(sorted(R.TOP_LEVEL_KEYS))})"
            )


def _check_runtime_resources(value: Any, errors: list[str]) -> None:
    if value is None:
        return
    if not isinstance(value, dict):
        errors.append("runtime_resources must be an object")
        return
    _check_mapping_schema("runtime_resources", value, R.RUNTIME_RESOURCE_SCHEMA, errors)
    for key in ("slots", "wait_timeout_seconds"):
        if key in value and value[key] is None:
            errors.append(f"runtime_resources.{key} must not be null")
    if "lock_directory" in value and (
        not isinstance(value["lock_directory"], str)
        or not value["lock_directory"].strip()
        or not Path(value["lock_directory"]).expanduser().is_absolute()
    ):
        errors.append("runtime_resources.lock_directory must be an absolute path")


def _check_step_timeouts(value: Any, errors: list[str]) -> None:
    """Validate optional outer subprocess budgets keyed by canonical step id."""
    if value is None:
        return
    if not isinstance(value, dict):
        errors.append("step_timeouts must be an object")
        return
    for step_id, timeout in value.items():
        path = f"step_timeouts.{step_id}"
        if step_id not in R.STEPS:
            errors.append(f"unknown step timeout key {step_id!r}")
            continue
        if type(timeout) is not int or timeout < 1:
            errors.append(f"{path} must be an integer >= 1")
        elif timeout > T.MAX_STEP_TIMEOUT_SECONDS:
            errors.append(
                f"{path} must be <= {T.MAX_STEP_TIMEOUT_SECONDS} seconds (30 days)"
            )


def _check_timeout_budgets(
    config: dict,
    steps_cfg: dict,
    errors: list[str],
) -> None:
    """Ensure an outer hard cap cannot terminate declared serial inner work."""
    overrides = config.get("step_timeouts")
    overrides = overrides if isinstance(overrides, dict) else {}
    for step_id in R.STEPS:
        if not _step_enabled(steps_cfg, step_id):
            continue
        try:
            automatic = T.automatic_step_timeout(step_id, config)
        except T.TimeoutBudgetError as exc:
            errors.append(str(exc))
            continue
        explicit = overrides.get(step_id)
        if type(explicit) is int and 1 <= explicit <= T.MAX_STEP_TIMEOUT_SECONDS:
            if explicit < automatic:
                errors.append(
                    f"step_timeouts.{step_id}={explicit} is below the required "
                    f"{automatic}s budget for enabled serial tools"
                )


def _check_mapping_schema(scope: str, mapping: dict, schema: dict, errors: list[str]) -> None:
    """Validate a flat key→value mapping (context/project) against a param schema."""
    if not isinstance(mapping, dict):
        errors.append(f"{scope} must be an object")
        return
    for key, value in mapping.items():
        if key not in schema:
            errors.append(f"unknown {scope} key {key!r}")
            continue
        _check_param_value(f"{scope}.{key}", value, schema[key], errors)


_CURATED_INTERACTION_FIELDS = {
    "residue", "chain", "icode", "interaction_type",
    "constraint_id",
    "protein_atom_names", "ligand_atom_names", "protein_role", "ligand_role",
    "protein_atoms", "ligand_atoms", "distance_range", "angle_range",
    "mediator", "rationale", "source", "binding_mode_scope", "review_status",
    "provenance", "scope", "logic_group", "logic_operator", "at_least_k",
    "group_id", "group_logic", "logic",
}
_CURATED_REQUIRED_FIELDS = {
    "residue", "interaction_type",
}
_CANONICAL_INTERACTION_TYPES = set(INTERACTION_TYPE_CANONICAL.values())


def _check_known_site_provenance(context: dict, errors: list[str]) -> None:
    value = context.get("known_site_provenance")
    if value is None or not isinstance(value, dict):
        return
    path = "context.known_site_provenance"
    allowed = {"source", "method", "independent_of_reference_structure", "review_status"}
    for key in value:
        if key not in allowed:
            errors.append(f"{path}: unknown field {key!r}")
    for key in ("source", "method", "review_status"):
        if key in value and (not isinstance(value[key], str) or not value[key].strip()):
            errors.append(f"{path}.{key} must be a non-empty string")
    independent = value.get("independent_of_reference_structure")
    if independent is not None and not isinstance(independent, bool):
        errors.append(f"{path}.independent_of_reference_structure must be bool")


def _check_known_required_interactions(context: dict, errors: list[str]) -> None:
    """Validate the evidence contract for user-curated hard constraints.

    A detector observation is not a biological requirement.  Anything entering
    ``known_required_interactions`` must therefore carry provenance, scope and
    an explicit approval state before Step 1 may compile it into a hard gate.
    """
    values = context.get("known_required_interactions")
    if values is None or not isinstance(values, list):
        return
    group_specs: dict[str, list[tuple[str, str, int | None]]] = {}
    constraint_ids: dict[str, str] = {}
    for index, item in enumerate(values):
        path = f"context.known_required_interactions[{index}]"
        if not isinstance(item, dict):
            # The generic list-item validator already reports this.
            continue
        for key in item:
            if key not in _CURATED_INTERACTION_FIELDS:
                errors.append(f"{path}: unknown field {key!r}")
        constraint_id = item.get("constraint_id")
        if isinstance(constraint_id, str) and constraint_id:
            if constraint_id in constraint_ids:
                errors.append(
                    f"{path}.constraint_id duplicates {constraint_ids[constraint_id]}"
                )
            else:
                constraint_ids[constraint_id] = path
        for key in sorted(_CURATED_REQUIRED_FIELDS):
            value = item.get(key)
            if not isinstance(value, str) or not value.strip():
                errors.append(f"{path}.{key} must be a non-empty string")

        provenance = item.get("provenance")
        source = item.get("source")
        nested_source = None
        nested_rationale = None
        if provenance is not None:
            if not isinstance(provenance, dict):
                errors.append(f"{path}.provenance must be an object")
            else:
                nested_source = (
                    provenance.get("reference") or provenance.get("source")
                    or provenance.get("citation") or provenance.get("record_id")
                )
                nested_rationale = provenance.get("rationale")
                if not isinstance(provenance.get("source_type"), str) or not provenance["source_type"].strip():
                    errors.append(f"{path}.provenance.source_type must be a non-empty string")
        if not (
            isinstance(source, str) and source.strip()
            or isinstance(nested_source, str) and nested_source.strip()
        ):
            errors.append(
                f"{path} requires a non-empty source or provenance.reference"
            )
        rationale = item.get("rationale") or nested_rationale
        if not isinstance(rationale, str) or not rationale.strip():
            errors.append(f"{path} requires a non-empty rationale")

        scope = item.get("scope")
        binding_mode_scope = item.get("binding_mode_scope")
        nested_scope = None
        if scope is not None:
            if not isinstance(scope, dict):
                errors.append(f"{path}.scope must be an object")
            else:
                nested_scope = scope.get("binding_mode") or scope.get("binding_mode_scope")
        if not (
            isinstance(binding_mode_scope, str) and binding_mode_scope.strip()
            or isinstance(nested_scope, str) and nested_scope.strip()
        ):
            errors.append(
                f"{path} requires a non-empty binding_mode_scope or scope.binding_mode"
            )

        interaction_type = item.get("interaction_type")
        if isinstance(interaction_type, str) and interaction_type.strip():
            canonical = normalize_interaction_type(interaction_type)
            if canonical not in _CANONICAL_INTERACTION_TYPES:
                errors.append(
                    f"{path}.interaction_type {interaction_type!r} is not a supported canonical type"
                )
        review_status = item.get("review_status")
        if review_status is None and isinstance(provenance, dict):
            review_status = provenance.get("review_status")
        if review_status != "approved":
            errors.append(f"{path}.review_status must be 'approved' for a required constraint")
        for key in (
            "chain", "icode", "constraint_id",
            "protein_role", "ligand_role", "mediator", "group_id", "group_logic", "logic",
        ):
            if key in item and not isinstance(item[key], str):
                errors.append(f"{path}.{key} must be a string")
        for key in ("protein_atom_names", "ligand_atom_names"):
            value = item.get(key)
            if value is not None and (
                not isinstance(value, list)
                or any(not isinstance(atom, str) or not atom for atom in value)
            ):
                errors.append(f"{path}.{key} must be a list of non-empty strings")
        for key in ("protein_atoms", "ligand_atoms"):
            value = item.get(key)
            if value is not None and (
                not isinstance(value, list) or any(not isinstance(atom, dict) for atom in value)
            ):
                errors.append(f"{path}.{key} must be a list of atom objects")
        for key, maximum in (("distance_range", None), ("angle_range", 180.0)):
            value = item.get(key)
            if value is None:
                continue
            valid_range = (
                isinstance(value, list) and len(value) == 2
                and all(isinstance(number, (int, float)) and not isinstance(number, bool) for number in value)
                and 0.0 <= float(value[0]) <= float(value[1])
                and (maximum is None or float(value[1]) <= maximum)
            )
            if not valid_range:
                suffix = " between 0 and 180" if maximum is not None else ""
                errors.append(f"{path}.{key} must be an ascending two-number range{suffix}")

        logic_group = item.get("logic_group")
        group_id = item.get("group_id")
        nested_operator = None
        nested_at_least_k = None
        if logic_group is not None and not isinstance(logic_group, (str, dict)):
            errors.append(f"{path}.logic_group must be a string or object")
        elif isinstance(logic_group, dict):
            unknown = set(logic_group) - {"group_id", "id", "logic", "operator", "at_least_k"}
            for key in sorted(unknown):
                errors.append(f"{path}.logic_group: unknown field {key!r}")
            group_id = logic_group.get("group_id") or logic_group.get("id") or group_id
            if group_id is not None and not isinstance(group_id, str):
                errors.append(f"{path}.logic_group.group_id must be a string")
            nested_operator = logic_group.get("logic") or logic_group.get("operator")
            nested_at_least_k = logic_group.get("at_least_k")
        elif isinstance(logic_group, str):
            group_id = logic_group

        operator = (
            item.get("logic_operator") or item.get("group_logic")
            or item.get("logic") or nested_operator
        )
        if operator is not None and operator not in {"all_of", "any_of", "at_least_k"}:
            errors.append(
                f"{path}.logic_operator must be one of all_of, any_of, at_least_k"
            )
        at_least_k = item.get("at_least_k", nested_at_least_k)
        if at_least_k is not None and (
            not isinstance(at_least_k, int) or isinstance(at_least_k, bool) or at_least_k < 1
        ):
            errors.append(f"{path}.at_least_k must be a positive integer")
        if operator == "at_least_k" and at_least_k is None:
            errors.append(f"{path}.at_least_k is required for logic_operator='at_least_k'")
        if operator in {"all_of", "any_of"} and at_least_k is not None:
            errors.append(
                f"{path}.at_least_k is only valid with logic_operator='at_least_k'"
            )
        if isinstance(group_id, str) and group_id:
            group_specs.setdefault(group_id, []).append((
                path, str(operator or "all_of"),
                at_least_k if isinstance(at_least_k, int) and not isinstance(at_least_k, bool) else None,
            ))

    for group_id, specs in group_specs.items():
        operators = {operator for _path, operator, _k in specs}
        if len(operators) > 1:
            errors.append(
                f"known_required_interactions logic group {group_id!r} has conflicting "
                f"operators: {', '.join(sorted(operators))}"
            )
        k_values = {k for _path, _operator, k in specs if k is not None}
        if len(k_values) > 1:
            errors.append(
                f"known_required_interactions logic group {group_id!r} has conflicting "
                f"at_least_k values: {', '.join(map(str, sorted(k_values)))}"
            )
        if k_values and next(iter(k_values)) > len(specs):
            errors.append(
                f"known_required_interactions logic group {group_id!r} at_least_k "
                f"cannot exceed its {len(specs)} members"
            )


def _check_scientific_facts(context: dict, errors: list[str]) -> None:
    facts = context.get("scientific_facts") or {}
    if not isinstance(facts, dict):
        errors.append("context.scientific_facts must be an object")
        return
    _check_mapping_schema(
        "context.scientific_facts", facts, R.SCIENTIFIC_FACTS_SCHEMA, errors,
    )


def _check_input_artifacts(
    injected: dict,
    steps_cfg: dict,
    errors: list[str],
) -> None:
    if not isinstance(injected, dict):
        errors.append("input_artifacts must be an object")
        return
    for name, value in injected.items():
        if name not in R.INJECTABLE_ARTIFACTS:
            errors.append(
                f"input_artifacts.{name}: not an injectable artifact "
                f"(allowed: {', '.join(sorted(R.INJECTABLE_ARTIFACTS))})"
            )
        elif not isinstance(value, str) or not value:
            errors.append(f"input_artifacts.{name} must be a non-empty file path string")
        else:
            try:
                path = Path(value).expanduser()
                exists = path.is_file()
            except (OSError, ValueError):
                exists = False
            if not exists:
                errors.append(
                    f"input_artifacts.{name} file does not exist or is not readable: {value}"
                )
            enabled_producers = sorted(
                step_id
                for step_id in _all_producers(name)
                if _producer_enabled(steps_cfg, step_id, name)
            )
            if enabled_producers:
                errors.append(
                    f"input_artifacts.{name} conflicts with enabled producer(s) "
                    f"{enabled_producers}; disable/bound those producer steps so "
                    "the injected generation cannot be silently overwritten"
                )


def _check_benchmark_labels(injected: dict, steps_cfg: dict, errors: list[str]) -> None:
    """Check existing external SDF titles before spending compute on a benchmark."""
    if not _tool_enabled(steps_cfg, "coordination", "screening_benchmark"):
        return
    if not isinstance(injected, dict):
        return
    value = injected.get("prefiltered_library") or injected.get("library")
    if not isinstance(value, str) or not value:
        return  # Generated/reused cohorts are checked by coordination.
    records = []
    try:
        with Path(value).expanduser().open(encoding="utf-8") as handle:
            title = True
            for line in handle:
                if title:
                    record = {"mol_name": line.strip()}
                    records.append(record)
                    title = False
                elif line.strip() == "$$$$":
                    title = True
                elif line.startswith(">"):
                    for field in ("input_mol_name", "compound_id"):
                        if f"<{field}>" in line:
                            record[field] = next(handle, "").strip()
                            break
    except (OSError, UnicodeError, ValueError):
        return  # The artifact/file checks own unreadable inputs.
    counts = {prefix: sum(benchmark_label(record) == prefix[:3] for record in records)
              for prefix in ("ACT_", "DEC_")}
    if not all(counts.values()):
        errors.append(
            "coordination.screening_benchmark requires verified ACT_ and DEC_ "
            f"labels in the supplied SDF (found {counts['ACT_']} ACT_, "
            f"{counts['DEC_']} DEC_). Add prefixes from a traceable label source "
            "before normalization, or disable screening_benchmark for unlabeled screening."
        )


# --- VR2/VR3: params ----------------------------------------------------------
def _check_step_params(step_id, step_spec, cfg_step, errors) -> None:
    for key in cfg_step:
        if key not in {"enabled", "params", "tools"}:
            errors.append(f"unknown config key {key!r} for step {step_id!r}")
    if "enabled" in cfg_step and not isinstance(cfg_step["enabled"], bool):
        errors.append(f"{step_id}.enabled must be bool")
    if not isinstance(cfg_step.get("params") or {}, dict):
        errors.append(f"{step_id}.params must be an object")
        return
    schema = step_spec.get("params", {})
    params = cfg_step.get("params") or {}
    for key, spec in schema.items():
        if spec.get("required") and params.get(key) is None:
            errors.append(f"{step_id}.{key} is required")
    for key, value in params.items():
        if key not in schema:
            errors.append(f"unknown param {key!r} for step {step_id!r}")
            continue
        _check_param_value(f"{step_id}.{key}", value, schema[key], errors)


def _check_tools(step_id, step_spec, cfg_step, errors, context) -> None:
    tools_schema = step_spec["tools"]
    step_enabled = bool(cfg_step.get("enabled", False))
    for tool_id, cfg_tool in (cfg_step.get("tools") or {}).items():
        if tool_id not in tools_schema:  # VR1
            from receptor_preparation import DEFERRED_TOOLS
            if step_id == "structure_prep" and tool_id in DEFERRED_TOOLS | {"amber_topology"}:
                destination = "pose_evaluation.candidate_topology" if tool_id == "amber_topology" else f"docking.{tool_id}"
                errors.append(f"structure_prep.{tool_id} moved to {destination}; regenerate the explicit tool matrix")
            else:
                errors.append(f"unknown tool {tool_id!r} in step {step_id!r}")
            continue
        if not isinstance(cfg_tool, dict):
            errors.append(f"tool {step_id}.{tool_id} config must be an object")
            continue
        for key in cfg_tool:
            if key not in {"enabled", "params"}:
                errors.append(f"unknown config key {key!r} for tool {step_id}.{tool_id}")
        if "enabled" in cfg_tool and not isinstance(cfg_tool["enabled"], bool):
            errors.append(f"{step_id}.{tool_id}.enabled must be bool")
        if not isinstance(cfg_tool.get("params") or {}, dict):
            errors.append(f"{step_id}.{tool_id}.params must be an object")
            continue
        pschema = tools_schema[tool_id].get("params", {})
        tool_params = cfg_tool.get("params") or {}
        tool_enabled_now = bool(
            cfg_tool.get("enabled", False)
        )
        if step_enabled and tool_enabled_now:
            for key, spec in pschema.items():
                if spec.get("required") and tool_params.get(key) is None:
                    errors.append(f"{step_id}.{tool_id}.{key} is required")
        for key, value in tool_params.items():
            if key not in pschema:  # VR2
                errors.append(f"unknown param {key!r} for {step_id}.{tool_id}")
                continue
            _check_param_value(f"{step_id}.{tool_id}.{key}", value, pschema[key], errors)

        # VR7: context prerequisite when this tool is enabled
        if step_enabled and tool_enabled_now:
            need = tools_schema[tool_id].get("context_one_of")
            if need and not any(context.get(k) for k in need):
                errors.append(
                    f"{step_id}.{tool_id} needs one of context.{{{', '.join(need)}}}"
                )


def _check_param_value(path: str, value: Any, spec: dict, errors: list[str]) -> None:
    ptype = spec["type"]
    if value is None:
        if spec.get("required"):
            errors.append(f"{path} is required")
        return
    ok = _type_ok(value, ptype)
    if not ok:
        errors.append(f"{path} must be {ptype}, got {type(value).__name__}")
        return
    if ptype in ("int", "float"):
        if spec.get("min") is not None and value < spec["min"]:
            errors.append(f"{path} must be >= {spec['min']} (got {value})")
        if spec.get("max") is not None and value > spec["max"]:
            errors.append(f"{path} must be <= {spec['max']} (got {value})")
    if spec.get("choices") and value not in spec["choices"]:
        errors.append(f"{path} must be one of {spec['choices']} (got {value!r})")
    if ptype == "list" and spec.get("item"):
        for index, item in enumerate(value):
            if not _type_ok(item, spec["item"]):
                errors.append(
                    f"{path}[{index}] must be {spec['item']}, got {type(item).__name__}"
                )
            elif spec.get("item_choices") and item not in spec["item_choices"]:
                errors.append(
                    f"{path}[{index}] must be one of {spec['item_choices']} "
                    f"(got {item!r})"
                )


def _type_ok(value: Any, ptype: str) -> bool:
    if ptype == "bool":
        return isinstance(value, bool)
    if ptype == "int":
        return isinstance(value, int) and not isinstance(value, bool)
    if ptype == "float":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if ptype == "str":
        return isinstance(value, str)
    if ptype == "list":
        return isinstance(value, list)
    if ptype == "dict":
        return isinstance(value, dict)
    return True


# --- VR4: step DAG ------------------------------------------------------------
def _step_enabled(steps_cfg: dict, step_id: str) -> bool:
    step = steps_cfg.get(step_id, {})
    if not isinstance(step, dict):
        return False
    return bool(step.get("enabled", False))


def _check_dag(steps_cfg: dict, available: set[str], errors: list[str]) -> None:
    """Each enabled step's required artifacts must be satisfiable.

    A step declares its needs one of two ways:
      - requires_all / requires_any (preferred): hard list + at-least-one groups;
      - otherwise `consumes` is treated as an all-hard list (legacy default).

    An artifact is satisfied when it is pre-supplied at entry (`available` — from
    upstream_project_dir preload or input_artifacts injection) OR produced by an
    enabled step of lower order. An artifact that NO step produces at all is an
    external/context input and is never flagged.
    """
    for step_id, spec in R.STEPS.items():
        if not _step_enabled(steps_cfg, step_id):
            continue
        order = spec["order"]

        def satisfied(artifact: str) -> bool:
            if artifact in available:
                return True
            producers = _all_producers(artifact)
            if not producers:
                return True  # external/context artifact
            return any(
                _producer_enabled(steps_cfg, p, artifact)
                and R.STEPS[p]["order"] < order
                for p in producers
            )

        hard = spec.get("requires_all")
        any_groups = spec.get("requires_any")
        if hard is None and any_groups is None:
            hard = spec["consumes"]  # legacy: all consumed artifacts are hard
            any_groups = []

        for artifact in (hard or []):
            if not satisfied(artifact):
                errors.append(
                    f"{step_id} requires {artifact!r} but no enabled upstream step "
                    f"produces it and it was not supplied (producers: {sorted(_all_producers(artifact))})"
                )
        for group in (any_groups or []):
            if not any(satisfied(a) for a in group):
                errors.append(
                    f"{step_id} requires at least one of "
                    f"{{{', '.join(group)}}} but none is produced upstream or supplied"
                )

        for tool_id, tool_spec in spec["tools"].items():
            enabled = _tool_enabled(steps_cfg, step_id, tool_id)
            if not enabled:
                continue
            for artifact in tool_spec.get("requires", []):
                if not satisfied(artifact):
                    errors.append(
                        f"{step_id}.{tool_id} requires {artifact!r} but no enabled "
                        "upstream step produces it and it was not supplied "
                        f"(producers: {sorted(_all_producers(artifact))})"
                    )
            for artifact in tool_spec.get("requires_supplied", []):
                if artifact not in available:
                    errors.append(
                        f"{step_id}.{tool_id} requires {artifact!r} to be "
                        "supplied through input_artifacts or an upstream-project "
                        "preload; no reviewed in-pipeline producer exists"
                    )


def _all_producers(artifact: str) -> list[str]:
    return [name for name, spec in R.STEPS.items() if artifact in spec["produces"]]


def _producer_enabled(
    steps_cfg: dict,
    step_id: str,
    artifact: str,
) -> bool:
    """Honor tool-conditional outputs such as Meeko's receptor PDBQT."""
    if not _step_enabled(steps_cfg, step_id):
        return False
    if (
        step_id == "pocket_analysis" and artifact == "docking_box"
        and (steps_cfg[step_id].get("params") or {}).get("site_selection_strategy") == "provided_box"
    ):
        return False
    step_spec = R.STEPS[step_id]
    conditional_tools = [
        tool_id
        for tool_id, tool_spec in step_spec["tools"].items()
        if artifact in tool_spec.get("produces", [])
    ]
    if not conditional_tools:
        return True
    return any(
        _tool_enabled(steps_cfg, step_id, tool_id)
        for tool_id in conditional_tools
    )


# --- VR5: same-step tool dependency ------------------------------------------
def _check_tool_dependencies(steps_cfg: dict, errors: list[str]) -> None:
    for step_id, spec in R.STEPS.items():
        if not _step_enabled(steps_cfg, step_id):
            continue
        for tool_id, tool_spec in spec["tools"].items():
            enabled = _tool_enabled(steps_cfg, step_id, tool_id)
            if not enabled:
                continue
            for dep in tool_spec.get("depends_on", []):
                dep_enabled = _tool_enabled(steps_cfg, step_id, dep)
                if not dep_enabled:
                    errors.append(
                        f"{step_id}.{tool_id} requires {step_id}.{dep} to be enabled"
                    )
            # tool_one_of: this tool needs at least one of a group enabled
            for group in spec.get("tool_one_of", {}).get(tool_id, []):
                if not any(
                    _tool_enabled(steps_cfg, step_id, g)
                    for g in group
                ):
                    errors.append(
                        f"{step_id}.{tool_id} requires at least one of "
                        f"{step_id}.{{{', '.join(group)}}} enabled"
                    )


# --- VR6: class "at least one" -----------------------------------------------
def _check_one_of(steps_cfg: dict, errors: list[str]) -> None:
    for step_id, spec in R.STEPS.items():
        if not _step_enabled(steps_cfg, step_id):
            continue
        if (
            step_id == "pocket_analysis"
            and (steps_cfg[step_id].get("params") or {}).get("site_selection_strategy") == "provided_box"
        ):
            continue
        for group in spec.get("one_of", []):
            if not any(
                _tool_enabled(steps_cfg, step_id, g)
                for g in group
            ):
                errors.append(
                    f"{step_id} needs at least one of {{{', '.join(group)}}} enabled"
                )


# --- VR8: exactly one tool in a mutually-exclusive class ---------------------
def _check_exactly_one(steps_cfg: dict, errors: list[str]) -> None:
    for step_id, spec in R.STEPS.items():
        if not _step_enabled(steps_cfg, step_id):
            continue
        for group in spec.get("exactly_one", []):
            enabled = [
                tool_id for tool_id in group
                if _tool_enabled(steps_cfg, step_id, tool_id)
            ]
            if len(enabled) != 1:
                errors.append(
                    f"{step_id} needs exactly one of "
                    f"{{{', '.join(group)}}} enabled (got {enabled or 'none'})"
                )


def _check_disabled_step_opt_ins(
    config: dict,
    steps_cfg: dict,
    errors: list[str],
) -> None:
    """Reject explicit tool selections that a disabled parent step would ignore.

    Tools inside a disabled step do not execute. A tool explicitly enabled by
    the Agent, however, must not disappear silently under a disabled parent.
    Steps before start_step are exempt: round-N configs retain upstream producer
    tool/parameter choices for lineage and cache validation. Steps after
    stop_step are future work, so an enabled tool there remains an invalid no-op.
    """
    start = config.get("start_step")
    try:
        lower = float(start) if start is not None else float("-inf")
    except (TypeError, ValueError):
        # merge_user_config/apply_step_range normally reports this first. Keep
        # direct validate() calls fail-collecting rather than raising here.
        lower = float("-inf")
    for step_id, step_spec in R.STEPS.items():
        if _step_enabled(steps_cfg, step_id):
            continue
        if step_spec["order"] < lower:
            continue
        enabled_tools = [
            tool_id
            for tool_id in step_spec["tools"]
            if _tool_enabled(steps_cfg, step_id, tool_id)
        ]
        if enabled_tools:
            errors.append(
                f"disabled step {step_id!r} contains enabled tools "
                f"{enabled_tools}; enable the parent step or disable those tools"
            )


# --- Step 0 execution contract + capability gates ----------------------------
def _tool_enabled(steps_cfg: dict, step_id: str, tool_id: str) -> bool:
    step = steps_cfg.get(step_id, {})
    step = step if isinstance(step, dict) else {}
    tools = step.get("tools", {})
    tools = tools if isinstance(tools, dict) else {}
    tool = tools.get(tool_id, {})
    tool = tool if isinstance(tool, dict) else {}
    return bool(tool.get("enabled", False))


def _check_structure_prep_contract(steps_cfg: dict, context: dict, errors: list[str]) -> None:
    if not _step_enabled(steps_cfg, "structure_prep"):
        return
    intake_enabled = _tool_enabled(steps_cfg, "structure_prep", "intake")
    has_remote_source = bool(context.get("pdb_id") or context.get("protein_name"))
    has_local_source = bool(context.get("receptor_file"))
    if intake_enabled and not has_remote_source:
        errors.append(
            "structure_prep.intake is enabled but context.pdb_id/protein_name is missing"
        )
    if not intake_enabled and not has_local_source:
        errors.append(
            "structure_prep.intake is disabled but context.receptor_file is missing"
        )
    if (
        not intake_enabled
        and _tool_enabled(steps_cfg, "structure_prep", "ligand_sdf")
        and not context.get("ligand_file")
    ):
        errors.append(
            "structure_prep.ligand_sdf is enabled in prep-only mode but "
            "context.ligand_file is missing; provide it or explicitly disable "
            "structure_prep.ligand_sdf"
        )

    step_cfg = steps_cfg.get("structure_prep", {})
    intake_params = step_cfg.get("tools", {}).get("intake", {}).get("params", {})
    if intake_params.get("skip_prep") is True:
        if not intake_enabled:
            errors.append("structure_prep.intake.skip_prep=true requires intake enabled")
        for tool_id in (
            "pdbfixer", "propka", "reduce", "pdb2pqr", "protoss",
            "ligand_sdf", "ligand_microstates", "openbabel",
            "metal_validation", "meeko", "amber_topology", "structure_qc",
            "residue_validation",
        ):
            if _tool_enabled(steps_cfg, "structure_prep", tool_id):
                errors.append(
                    "structure_prep.intake.skip_prep=true conflicts with enabled "
                    f"structure_prep.{tool_id}"
                )

    pdbfixer_ph = (
        step_cfg.get("tools", {}).get("pdbfixer", {}).get("params", {}).get("ph")
    )
    propka_ph = (
        step_cfg.get("tools", {}).get("propka", {}).get("params", {}).get("ph")
    )
    if pdbfixer_ph is not None and propka_ph is not None and pdbfixer_ph != propka_ph:
        errors.append(
            "structure_prep pdbfixer.ph and propka.ph must match; prefer the shared "
            "structure_prep.params.ph"
        )

    pdbfixer_params = step_cfg.get("tools", {}).get("pdbfixer", {}).get("params", {})
    if pdbfixer_params.get("add_hydrogens") is True:
        errors.append("structure_prep.pdbfixer.params.add_hydrogens must be false; final protonation belongs to docking_prep")
    if (step_cfg.get("tools", {}).get("structure_qc", {}).get("params") or {}).get("require_hydrogens") is True:
        errors.append("structure_prep.structure_qc cannot require final hydrogens on the canonical noH geometry")
    if pdbfixer_params.get("add_missing_residues") and not str(pdbfixer_params.get("repair_reason") or "").strip():
        errors.append("structure_prep.pdbfixer.add_missing_residues requires repair_reason for supported short internal gaps")
    if _tool_enabled(steps_cfg, "structure_prep", "propka"):
        ph = propka_ph if propka_ph is not None else step_cfg.get("params", {}).get("ph")
        if ph is None:
            errors.append("Diagnostic structure_prep.propka requires an explicit justified pH")
    if intake_params.get("skip_prep") and any(_step_enabled(steps_cfg, name) for name in ("pocket_analysis", "molecule_generation", "docking")):
        errors.append("intake.skip_prep cannot supply the canonical noH receptor; run geometry preparation before downstream steps")


def _check_blueprint_evidence_policy(
    steps_cfg: dict, context: dict, errors: list[str],
) -> None:
    """Keep hard-gate authorization reviewed and non-circular."""
    if not _step_enabled(steps_cfg, "pocket_analysis"):
        return
    params = (
        steps_cfg.get("pocket_analysis", {})
        .get("tools", {})
        .get("blueprint", {})
        .get("params", {})
    )
    auto_required = params.get(
        "auto_required_from_tool_consensus",
        R.STEPS["pocket_analysis"]["tools"]["blueprint"]["params"]
        ["auto_required_from_tool_consensus"]["default"],
    )
    require_curated = params.get(
        "require_curated_provenance",
        R.STEPS["pocket_analysis"]["tools"]["blueprint"]["params"]
        ["require_curated_provenance"]["default"],
    )
    if require_curated is not True:
        errors.append(
            "pocket_analysis.blueprint.params.require_curated_provenance must remain true; "
            "unreviewed observations cannot authorize a hard gate"
        )
    if not auto_required:
        return
    provenance = context.get("known_site_provenance") or {}
    if provenance.get("independent_of_reference_structure") is not True:
        errors.append(
            "pocket_analysis.blueprint.params.auto_required_from_tool_consensus=true "
            "requires context.known_site_provenance.independent_of_reference_structure=true"
        )
    if provenance.get("review_status") != "approved":
        errors.append(
            "pocket_analysis.blueprint.params.auto_required_from_tool_consensus=true "
            "requires context.known_site_provenance.review_status='approved'"
        )


def _check_scientific_capabilities(
    context: dict,
    steps_cfg: dict,
    errors: list[str],
) -> None:
    """Reject unsupported chemistry without selecting a replacement tool."""
    facts = context.get("scientific_facts") or {}
    if not isinstance(facts, dict):
        return
    if (
        facts.get("has_metal") is True
        and _step_enabled(steps_cfg, "structure_prep")
        and _tool_enabled(steps_cfg, "structure_prep", "amber_topology")
    ):
        errors.append(
            "generic structure_prep.amber_topology does not parameterize declared "
            "metal sites; use a reviewed MCPB.py-specific workflow"
        )
    if facts.get("is_covalent") is True and _step_enabled(steps_cfg, "docking"):
        errors.append(
            "context.scientific_facts.is_covalent=true is incompatible with the "
            "registered non-covalent docking engines; reactive-atom selection, bond "
            "geometry, and covalent docking are not implemented"
        )


def _artifact_satisfied_for_step(
    artifact: str,
    consumer_step: str,
    steps_cfg: dict,
    available: set[str],
) -> bool:
    if artifact in available:
        return True
    order = R.STEPS[consumer_step]["order"]
    return any(
        R.STEPS[producer]["order"] < order
        and _producer_enabled(steps_cfg, producer, artifact)
        for producer in _all_producers(artifact)
    )


def _check_explicit_modes(
    steps_cfg: dict,
    context: dict,
    available: set[str],
    errors: list[str],
) -> None:
    """Validate explicit algorithm/input modes; never infer a mode from artifacts."""
    if _step_enabled(steps_cfg, "pocket_analysis"):
        pocket_cfg = steps_cfg.get("pocket_analysis") or {}
        strategy = (pocket_cfg.get("params") or {}).get("site_selection_strategy")
        geometry_tools = ("fpocket", "p2rank", "pykvfinder", "deeppocket")
        if strategy == "provided_box":
            if "docking_box" not in available:
                errors.append("pocket_analysis.site_selection_strategy='provided_box' requires a supplied docking_box")
            if (
                _step_enabled(steps_cfg, "docking")
                and (steps_cfg["docking"].get("params") or {}).get("search_space_source") != "docking_box"
            ):
                errors.append("provided_box pocket analysis requires docking.search_space_source='docking_box'")
        elif strategy == "ligand_site":
            if not _tool_enabled(steps_cfg, "pocket_analysis", "ligand_site"):
                errors.append(
                    "pocket_analysis.site_selection_strategy='ligand_site' requires "
                    "pocket_analysis.ligand_site enabled"
                )
        elif strategy == "known_site":
            if not context.get("known_site_residues"):
                errors.append(
                    "pocket_analysis.site_selection_strategy='known_site' requires "
                    "context.known_site_residues"
                )
            if not any(
                _tool_enabled(steps_cfg, "pocket_analysis", tool_id)
                for tool_id in geometry_tools
            ):
                errors.append(
                    "pocket_analysis.site_selection_strategy='known_site' requires at "
                    "least one geometry detector"
                )
        elif strategy == "detector_consensus":
            if not _tool_enabled(steps_cfg, "pocket_analysis", "deeppocket"):
                errors.append("detector_consensus requires the default DeepPocket anchor enabled")
            if (pocket_cfg.get("params") or {}).get("pocket_memory_mode", "auto") == "off":
                if not all(_tool_enabled(steps_cfg, "pocket_analysis", tool_id) for tool_id in geometry_tools):
                    errors.append("memory-off detector_consensus requires DeepPocket, fpocket, P2Rank and pyKVFinder")
                if _tool_enabled(steps_cfg, "pocket_analysis", "autosite"):
                    errors.append("memory-off detector_consensus does not execute AutoSite")
            # Exact memory-guided routing depends on the bound receptor/query.
            # Step 1 verifies every switch against plan-pocket-tools before execution.
        if strategy != "detector_consensus" and _tool_enabled(steps_cfg, "pocket_analysis", "autosite"):
            errors.append("AutoSite is currently registered only for detector_consensus")

    if not _step_enabled(steps_cfg, "docking"):
        return
    docking_cfg = steps_cfg.get("docking") or {}
    source = (docking_cfg.get("params") or {}).get("search_space_source")
    enabled_engines = [
        engine
        for engine in R.DOCKING_PRIMARY_ENGINES
        if _tool_enabled(steps_cfg, "docking", engine)
    ]
    if len(enabled_engines) != 1:
        return  # the generic exactly-one check reports this first
    engine = enabled_engines[0]
    allowed_sources = {
        "diffdock": {"whole_receptor"},
        "dock6": {"precomputed_grid"},
    }.get(engine, {"reference_ligand", "docking_box"})
    if source not in allowed_sources:
        errors.append(
            f"docking.search_space_source={source!r} is incompatible with "
            f"docking.{engine}; allowed: {sorted(allowed_sources)}"
        )
        return
    if (
        _tool_enabled(steps_cfg, "docking", "redocking_validation")
        and source != "reference_ligand"
    ):
        errors.append(
            "docking.redocking_validation requires "
            "docking.search_space_source='reference_ligand'"
        )

    required_artifacts = {
        "reference_ligand": ("reference_ligand",),
        "docking_box": ("docking_box",),
        "precomputed_grid": (
            "dock6_receptor_site",
            "dock6_grid_nrg",
            "dock6_grid_bmp",
        ),
        "whole_receptor": (),
    }.get(source, ())
    for artifact in required_artifacts:
        if not _artifact_satisfied_for_step(
            artifact, "docking", steps_cfg, available,
        ):
            errors.append(
                f"docking.search_space_source={source!r} requires {artifact!r}, "
                "but it is not produced upstream or supplied"
            )


def _check_pose_identity_contract(
    steps_cfg: dict,
    errors: list[str],
) -> None:
    """Do not merge pose evidence or scores without canonical identity QC."""
    if (
        _step_enabled(steps_cfg, "pose_evaluation")
        and not _tool_enabled(steps_cfg, "pose_evaluation", "identity_qc")
    ):
        errors.append(
            "pose_evaluation requires pose_evaluation.identity_qc enabled; "
            "pose evidence, score joins, and resumable batches must remain "
            "bound to canonical pose identity"
        )


def _check_pose_selection_contract(steps_cfg: dict, errors: list[str]) -> None:
    """Sequential fallback requires the final verdict within each QC batch."""
    if not _step_enabled(steps_cfg, "pose_evaluation"):
        return
    params = (steps_cfg.get("pose_evaluation") or {}).get("params") or {}
    strategy = params.get(
        "selection_strategy",
        R.STEPS["pose_evaluation"]["params"]["selection_strategy"]["default"],
    )
    if strategy != "first_passing_per_compound":
        return
    for tool_id in ("mdanalysis", "mmgbsa"):
        if _tool_enabled(steps_cfg, "pose_evaluation", tool_id):
            errors.append(
                "pose_evaluation.selection_strategy='first_passing_per_compound' "
                f"cannot be combined with pose_evaluation.{tool_id}: this tool "
                "can revise verdicts after pose selection; explicitly choose "
                "selection_strategy='all' for this evaluation plan"
            )


def _check_molecule_generation_contract(
    steps_cfg: dict,
    errors: list[str],
) -> None:
    """Require the writer of Step 2's canonical screening-library artifact."""
    if (
        _step_enabled(steps_cfg, "molecule_generation")
        and not _tool_enabled(steps_cfg, "molecule_generation", "conformer")
    ):
        errors.append(
            "molecule_generation requires molecule_generation.conformer enabled; "
            "it writes the canonical screening-ready library.sdf artifact"
        )


def _check_mmgbsa_contract(
    steps_cfg: dict,
    errors: list[str],
) -> None:
    """Keep an explicitly enabled MM-GBSA escalation executable and bounded."""
    if not _tool_enabled(steps_cfg, "pose_evaluation", "mmgbsa"):
        return
    if not _step_enabled(steps_cfg, "pose_evaluation"):
        errors.append(
            "pose_evaluation.mmgbsa requires pose_evaluation enabled; "
            "an explicit MM-GBSA request cannot be silently removed by a "
            "step range"
        )
        return
    if not _tool_enabled(
        steps_cfg,
        "pose_evaluation",
        "candidate_topology",
    ):
        errors.append(
            "pose_evaluation.mmgbsa requires "
            "pose_evaluation.candidate_topology enabled"
        )
        return

    tools = (
        (steps_cfg.get("pose_evaluation") or {}).get("tools") or {}
    )
    topology_params = (
        (tools.get("candidate_topology") or {}).get("params") or {}
    )
    mmgbsa_params = (tools.get("mmgbsa") or {}).get("params") or {}
    topology_spec = R.STEPS["pose_evaluation"]["tools"]["candidate_topology"]
    mmgbsa_spec = R.STEPS["pose_evaluation"]["tools"]["mmgbsa"]
    topology_mode = topology_params.get(
        "execution_mode",
        topology_spec["params"]["execution_mode"]["default"],
    )
    mmgbsa_mode = mmgbsa_params.get(
        "execution_mode",
        mmgbsa_spec["params"]["execution_mode"]["default"],
    )
    if topology_mode != "build":
        errors.append(
            "pose_evaluation.mmgbsa requires "
            "candidate_topology.params.execution_mode='build'"
        )
    topology_max = topology_params.get(
        "max_candidates",
        topology_spec["params"]["max_candidates"]["default"],
    )
    top_k = mmgbsa_params.get(
        "top_k",
        mmgbsa_spec["params"]["top_k"]["default"],
    )
    if (
        isinstance(topology_max, int)
        and not isinstance(topology_max, bool)
        and isinstance(top_k, int)
        and not isinstance(top_k, bool)
        and topology_max < top_k
    ):
        errors.append(
            "candidate_topology.params.max_candidates must be >= "
            "mmgbsa.params.top_k"
        )

    fail_fast = topology_params.get(
        "fail_fast",
        topology_spec["params"]["fail_fast"]["default"],
    )
    if mmgbsa_mode == "score" and fail_fast is not False:
        errors.append(
            "MM-GBSA score mode requires "
            "candidate_topology.params.fail_fast=false so one topology "
            "failure does not silently discard the remaining top-K cohort"
        )


def _check_short_md_contract(steps_cfg: dict, errors: list[str]) -> None:
    """Require a complete explicit-solvent topology route for short MD."""
    if not _tool_enabled(steps_cfg, "pose_evaluation", "openmm_md"):
        return
    if not _step_enabled(steps_cfg, "pose_evaluation"):
        errors.append(
            "pose_evaluation.openmm_md requires pose_evaluation enabled"
        )
        return
    if not _tool_enabled(
        steps_cfg, "pose_evaluation", "candidate_topology"
    ):
        errors.append(
            "pose_evaluation.openmm_md requires "
            "pose_evaluation.candidate_topology enabled"
        )
        return
    tools = ((steps_cfg.get("pose_evaluation") or {}).get("tools") or {})
    topology_params = (
        (tools.get("candidate_topology") or {}).get("params") or {}
    )
    md_params = (tools.get("openmm_md") or {}).get("params") or {}
    topology_spec = R.STEPS["pose_evaluation"]["tools"]["candidate_topology"]
    md_spec = R.STEPS["pose_evaluation"]["tools"]["openmm_md"]
    topology_mode = topology_params.get(
        "execution_mode",
        topology_spec["params"]["execution_mode"]["default"],
    )
    if topology_mode != "build":
        errors.append(
            "pose_evaluation.openmm_md requires "
            "candidate_topology.params.execution_mode='build'"
        )
    protein_force_field = topology_params.get(
        "protein_force_field",
        topology_spec["params"]["protein_force_field"]["default"],
    )
    if protein_force_field != "ff14SB":
        errors.append(
            "pose_evaluation.openmm_md currently validates only the "
            "ff14SB/TIP3P explicit-solvent protocol"
        )
    topology_max = topology_params.get(
        "max_candidates",
        topology_spec["params"]["max_candidates"]["default"],
    )
    top_k = md_params.get("top_k", md_spec["params"]["top_k"]["default"])
    if (
        isinstance(topology_max, int)
        and not isinstance(topology_max, bool)
        and isinstance(top_k, int)
        and not isinstance(top_k, bool)
        and topology_max < top_k
    ):
        errors.append(
            "candidate_topology.params.max_candidates must be >= "
            "openmm_md.params.top_k"
        )
    fail_fast = topology_params.get(
        "fail_fast", topology_spec["params"]["fail_fast"]["default"]
    )
    if fail_fast is not False:
        errors.append(
            "OpenMM short MD requires candidate_topology.params.fail_fast=false "
            "so one unsupported ligand does not hide the remaining bounded cohort"
        )
    if _tool_enabled(steps_cfg, "pose_evaluation", "mdanalysis"):
        production = md_params.get(
            "production_picoseconds",
            md_spec["params"]["production_picoseconds"]["default"],
        )
        interval = md_params.get(
            "report_interval_picoseconds",
            md_spec["params"]["report_interval_picoseconds"]["default"],
        )
        mda_params = (tools.get("mdanalysis") or {}).get("params") or {}
        mda_spec = R.STEPS["pose_evaluation"]["tools"]["mdanalysis"]
        stride = mda_params.get(
            "analysis_stride",
            mda_spec["params"]["analysis_stride"]["default"],
        )
        if (
            isinstance(production, (int, float))
            and not isinstance(production, bool)
            and isinstance(interval, (int, float))
            and not isinstance(interval, bool)
            and isinstance(stride, int)
            and not isinstance(stride, bool)
            and production < (stride + 1) * interval
        ):
            errors.append(
                "MDAnalysis stability analysis requires at least two production "
                "frames after applying analysis_stride: "
                "openmm_md.production_picoseconds must be >= "
                "(mdanalysis.analysis_stride + 1) * "
                "openmm_md.report_interval_picoseconds"
            )


def _check_generated_pose_qc(steps_cfg: dict, errors: list[str]) -> None:
    """Do not advance native 3D-generator candidates without Step 5 QC."""
    if not _step_enabled(steps_cfg, "molecule_generation"):
        return
    direct_3d_generators = (
        "pocket2mol", "resgen", "diffsbdd", "targetdiff", "difflinker",
        "lingo3dmol", "molcraft",
    )
    enabled = [
        tool_id for tool_id in direct_3d_generators
        if _tool_enabled(steps_cfg, "molecule_generation", tool_id)
    ]
    # A generation-only run may intentionally stop at Step 2. Once docking is
    # enabled, native generated geometries must reach both validators.
    if not enabled or not _step_enabled(steps_cfg, "docking"):
        return
    if not _step_enabled(steps_cfg, "pose_evaluation"):
        errors.append(
            "3D molecule generators "
            f"{enabled} require pose_evaluation enabled when docking is enabled"
        )
        return
    for tool_id in ("posebusters", "posecheck"):
        if not _tool_enabled(steps_cfg, "pose_evaluation", tool_id):
            errors.append(
                "3D molecule generators "
                f"{enabled} require pose_evaluation.{tool_id} enabled"
            )


def _check_metric_gate_contract(steps_cfg: dict, errors: list[str]) -> None:
    """Validate configurable soft/hard bounds preserve the marking order."""
    pose_spec = R.STEPS["pose_evaluation"]["tools"]["posecheck"]
    configured = (
        ((steps_cfg.get("pose_evaluation") or {}).get("tools") or {})
        .get("posecheck", {})
        .get("params", {})
    )
    params = {
        name: configured.get(name, definition.get("default"))
        for name, definition in pose_spec.get("params", {}).items()
    }
    for soft_name, hard_name in (
        ("strain_warn_kcal", "strain_fail_kcal"),
        ("clash_warn_count", "clash_fail_count"),
    ):
        soft = params.get(soft_name)
        hard = params.get(hard_name)
        if (
            isinstance(soft, (int, float))
            and not isinstance(soft, bool)
            and isinstance(hard, (int, float))
            and not isinstance(hard, bool)
            and soft > hard
        ):
            errors.append(
                "pose_evaluation.posecheck requires "
                f"{soft_name} <= {hard_name} (got {soft} > {hard})"
            )

    redocking_spec = R.STEPS["docking"]["tools"]["redocking_validation"]
    redocking_configured = (
        ((steps_cfg.get("docking") or {}).get("tools") or {})
        .get("redocking_validation", {})
        .get("params", {})
    )
    redocking_params = {
        name: redocking_configured.get(name, definition.get("default"))
        for name, definition in redocking_spec.get("params", {}).items()
    }
    redocking_soft = redocking_params.get("redocking_rmsd_threshold")
    redocking_hard = redocking_params.get("redocking_rmsd_hard_threshold")
    if (
        isinstance(redocking_soft, (int, float))
        and not isinstance(redocking_soft, bool)
        and isinstance(redocking_hard, (int, float))
        and not isinstance(redocking_hard, bool)
        and redocking_soft > redocking_hard
    ):
        errors.append(
            "docking.redocking_validation requires "
            "redocking_rmsd_threshold <= redocking_rmsd_hard_threshold "
            f"(got {redocking_soft} > {redocking_hard})"
        )

    coordination_spec = R.STEPS["coordination"]["params"]
    coordination_configured = (
        (steps_cfg.get("coordination") or {}).get("params") or {}
    )
    coordination_params = {
        name: coordination_configured.get(name, definition.get("default"))
        for name, definition in coordination_spec.items()
    }
    for hard_name, soft_name in (
        ("min_validation_auroc_hard", "min_validation_auroc"),
        ("min_validation_ef1_hard", "min_validation_ef1"),
    ):
        hard = coordination_params.get(hard_name)
        soft = coordination_params.get(soft_name)
        if (
            isinstance(hard, (int, float))
            and not isinstance(hard, bool)
            and isinstance(soft, (int, float))
            and not isinstance(soft, bool)
            and hard > soft
        ):
            errors.append(
                "coordination requires "
                f"{hard_name} <= {soft_name} (got {hard} > {soft})"
            )

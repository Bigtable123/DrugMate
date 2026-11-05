#!/usr/bin/env python3
"""Generate Agent-facing tool inventory and selection policy.

The registry is the source of truth for tool inventory, parameter schemas, and
declarative execution constraints.  ``pipeline_validation`` supplies the hard
contract summaries.  Tool choice itself is an Agent decision described by the
generated policy; this generator never derives a tool switch from scientific
facts.  It derives three synchronized artifacts:

  pipeline-entry/tool_catalog.json  — machine-readable: steps, tools, params,
                                       defaults, deps, artifacts, VR rules.
  pipeline-entry/TOOL_CATALOG.md    — human/Agent-readable summary tables.
  pipeline-entry/TOOL_SELECTION_POLICY.md — human/Agent-readable decision
                                             policy generated from those rules.

Run this whenever the registry or validator contract changes:

    python pipeline-entry/build_tool_catalog.py

A --check mode verifies all three committed artifacts are in sync with those
machine rules (for CI / pre-flight), exiting non-zero on drift.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

BIO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIO_ROOT / "shared"))

import pipeline_validation as V
from deployment_config import symbolic_deployment_paths


def _catalog_registry():
    """Load an independent registry view with symbolic deployment defaults.

    A separate module preserves the process's runtime registry and environment.
    Symbolic defaults remain unambiguous even when tools share a source path.
    """
    spec = importlib.util.spec_from_file_location(
        "_drugmate_catalog_registry", BIO_ROOT / "shared" / "pipeline_registry.py")
    registry = importlib.util.module_from_spec(spec)
    with symbolic_deployment_paths():
        spec.loader.exec_module(registry)
    return registry


R = _catalog_registry()

HERE = Path(__file__).resolve().parent
JSON_PATH = HERE / "tool_catalog.json"
MD_PATH = HERE / "TOOL_CATALOG.md"
POLICY_PATH = HERE / "TOOL_SELECTION_POLICY.md"


def build_catalog() -> dict:
    """Derive the machine-readable catalog dict from the registry."""
    steps = []
    for step_id in R.steps_in_order():
        spec = R.STEPS[step_id]
        tools = []
        for tool_id, tool in spec["tools"].items():
            tools.append({
                "id": tool_id,
                "default_enabled": tool["default_enabled"],
                "purpose": tool.get("purpose", ""),
                "when_to_use": tool.get("when_to_use", ""),
                "depends_on": tool.get("depends_on", []),
                "requires_artifacts": tool.get("requires", []),
                "requires_supplied_artifacts": tool.get(
                    "requires_supplied", []
                ),
                "produces_artifacts": tool.get("produces", []),
                "context_one_of": tool.get("context_one_of"),
                "runtime": tool.get("runtime"),
                "params": _params_catalog(tool.get("params", {})),
            })
        steps.append({
            "id": step_id,
            "order": spec["order"],
            "skill": spec["skill"],
            "runtime": spec.get("runtime"),
            "default_enabled": spec["default_enabled"],
            "best_effort": bool(spec.get("best_effort", False)),
            "consumes": spec["consumes"],
            "optional_consumes": spec.get("optional_consumes", []),
            "requires_all": spec.get("requires_all"),
            "requires_any": spec.get("requires_any", []),
            "produces": spec["produces"],
            "params": _params_catalog(spec.get("params", {})),
            "one_of": spec.get("one_of", []),
            "exactly_one": spec.get("exactly_one", []),
            "tool_one_of": spec.get("tool_one_of", {}),
            "tools": tools,
        })
    catalog = {
        "schema_version": 15,
        "generated_from": [
            "shared/pipeline_registry.py",
            "shared/pipeline_validation.py:VALIDATION_RULES",
        ],
        "artifacts": sorted(R.ARTIFACTS),
        "tool_environments": R.DEFAULT_TOOL_ENVIRONMENTS,
        "injectable_artifacts": sorted(R.INJECTABLE_ARTIFACTS),
        "top_level_keys": sorted(R.TOP_LEVEL_KEYS),
        "context_schema": _params_catalog(R.CONTEXT_SCHEMA),
        "scientific_facts_schema": _params_catalog(
            getattr(R, "SCIENTIFIC_FACTS_SCHEMA", {})
        ),
        "project_schema": _params_catalog(R.PROJECT_SCHEMA),
        "runtime_resources_schema": _params_catalog(R.RUNTIME_RESOURCE_SCHEMA),
        "selection_contract": {
            "tool_switches": "complete_explicit_matrix",
            "strategy_params": _selection_strategy_params(steps),
        },
        "metric_scoring": {
            "rules_version": R.METRIC_RULES_VERSION,
            "normalization_method": R.NORMALIZATION_METHOD,
            "rules": R.METRIC_RULES,
        },
        "validation_rules": V.VALIDATION_RULES,
        "config_shape": _CONFIG_SHAPE,
        "steps": steps,
    }
    return _portable_catalog_paths(catalog)


def _portable_catalog_paths(value):
    """Keep deployment paths symbolic in public Agent-facing catalogs."""
    if isinstance(value, dict):
        return {
            key: _portable_catalog_paths(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_portable_catalog_paths(item) for item in value]
    if isinstance(value, tuple):
        return [_portable_catalog_paths(item) for item in value]
    if isinstance(value, str):
        # Runtime placement is already symbolic; only bundled script paths
        # originate from __file__. Probe code may contain repr-escaped paths.
        # Respect the end of the directory component for either representation.
        roots = {str(BIO_ROOT), repr(str(BIO_ROOT))[1:-1]}
        alternatives = "|".join(re.escape(root) for root in sorted(roots, key=len, reverse=True))
        pattern = "(?:" + alternatives + r")(?=[/\\'\"]|$)"
        return re.sub(pattern, lambda _: "${BIO_SKILLS_ROOT}", value)
    return value


def _selection_strategy_params(steps: list[dict]) -> list[dict]:
    """Expose explicitly registered selection-strategy parameters, if present.

    Schema-driven discovery includes only fields declared in the registry.
    """
    strategy_names = {"site_selection_strategy", "search_space_source"}
    discovered = []
    for step in steps:
        for param in step.get("params", []):
            if param["name"] in strategy_names:
                discovered.append({"step": step["id"], **param})
    return discovered


def _params_catalog(param_schema: dict) -> list[dict]:
    out = []
    for name, spec in param_schema.items():
        entry = {
            "name": name,
            "type": spec["type"],
            "default": spec.get("default"),
            "required": bool(spec.get("required", False)),
        }
        for key in ("min", "max", "choices", "item", "item_choices", "desc"):
            if spec.get(key) is not None:
                entry[key] = spec[key]
        out.append(entry)
    return out


_CONFIG_SHAPE = {
    "project": {"project_slug": "str", "data_root": "str", "round_number": "int",
                "campaign_state_file": "absolute str|null",
                "campaign_stage": "R1|R2|R3|null",
                "campaign_attempt": "int >= 1|null",
                "round_decision_file": "absolute str|null",
                "upstream_project_dir": "str|null",
                "execution_mode": "interrupt|auto (default: interrupt)"},
    "context": {"pdb_id": "str", "target_name": "str", "chembl_target_id": "str",
                "uniprot_id": "str", "chain": "str", "seed": "int",
                "candidate_source": (
                    "molecule_generation|sampled_library|user_supplied_full|null"
                ),
                "scientific_facts": {
                    "has_metal": "bool|null",
                    "is_covalent": "bool|null",
                    "protonation_sensitive": "bool|null",
                },
                "known_site_residues": "list[str]",
                "known_site_provenance": "dict|null",
                "known_required_interactions": (
                    "list[{residue,interaction_type,rationale,source,"
                    "binding_mode_scope,review_status=approved,...}]"
                )},
    "step_timeouts": {
        "<step_id>": (
            "int >= automatic serial-tool budget and <= 2592000 "
            "(optional outer process cap)"
        ),
    },
    "runtime_resources": {"wait_timeout_seconds": "int >= 0 (default 0: fail fast)",
                          "slots": f"int 1-64 (shared host capacity; default {R.RUNTIME_RESOURCE_SCHEMA['slots']['default']})",
                          "lock_directory": "absolute shared path (optional)"},
    "steps": {"<step_id>": {
        "enabled": "bool",
        "tools": {"<tool_id>": {"enabled": "bool", "params": {"<name>": "<value>"}}},
        "params": {"<name>": "<value>"},
    }},
    "start_step": "int (0-8, optional convenience bound)",
    "stop_step": "int (0-8, optional convenience bound)",
}
# __CATALOG_MD__


def render_markdown(catalog: dict) -> str:
    """Render the human/Agent-readable catalog. Generated — do not edit by hand."""
    lines = [
        "# SBDD Tool Catalog (generated)",
        "",
        "> Generated from `shared/pipeline_registry.py` by "
        "`pipeline-entry/build_tool_catalog.py`. **Do not edit by hand** — "
        "re-run the generator after changing the registry.",
        "",
        "This is the inventory the Agent uses while building a pipeline config. "
        "For tool choice see the generated `TOOL_SELECTION_POLICY.md`; for assembling "
        "confirmed facts and choices into config see `QUERY_TO_CONFIG.md`; for the execution "
        "contract see `SKILL.md`.",
        "",
        "`${SBDD_*}` and `${BIO_SKILLS_ROOT}` are deployment templates in this catalog. "
        "Runtime defaults resolve them from the deployment settings. Explicit JSON "
        "path overrides must contain actual paths; JSON does not expand shell variables.",
        "",
        "## Top-level config keys",
        "",
        "`" + "`, `".join(catalog["top_level_keys"]) + "` — anything else hard-fails.",
        "",
        "## Tool runtime environments",
        "",
        "Runtime placement is registry-owned; Agent configs select tools but do not repeat Conda paths.",
        "",
        "| id | kind | location |",
        "|---|---|---|",
    ]
    for env_id, spec in catalog["tool_environments"].items():
        location = (
            spec.get("name")
            or spec.get("prefix")
            or spec.get("base_url")
            or "current process"
        )
        lines.append(f"| `{env_id}` | {spec.get('kind')} | `{location}` |")
    lines += [
        "",
        "Enabled runtimes are checked before any scientific step and recorded in `tool_runtime_report.json`.",
        "",
        "## `context` (shared identity / knowledge)",
        "",
        "| key | type | default | range/choices |",
        "|---|---|---|---|",
    ]
    for p in catalog["context_schema"]:
        lines.append(_schema_row(p))
    if catalog["scientific_facts_schema"]:
        lines += [
            "",
            "### `context.scientific_facts`",
            "",
            "These are declared scientific inputs for tools and capability checks. "
            "They never select or enable a tool.",
            "",
            "| key | type | default | range/choices |",
            "|---|---|---|---|",
        ]
        for p in catalog["scientific_facts_schema"]:
            lines.append(_schema_row(p))
    lines += [
        "",
        "## `project` (run / directory control)",
        "",
        "| key | type | default | range/choices |",
        "|---|---|---|---|",
    ]
    for p in catalog["project_schema"]:
        lines.append(_schema_row(p))
    lines += ["", "## `runtime_resources` (queue budget, separate from execution)", "",
              "| key | type | default | range/choices |", "|---|---|---|---|"]
    for p in catalog["runtime_resources_schema"]:
        lines.append(_schema_row(p))
    lines += [
        "",
        "## Explicit selection contract",
        "",
        "The catalog inventories choices; it does not select them. Agent-authored "
        "configs provide every step switch and a complete tool on/off matrix. "
        "Registry defaults are recommendation metadata for policy comparison; they do not "
        "participate in execution and cannot substitute for that explicit decision.",
        "Ordinarily, every tool under a disabled step is false. The only exception is "
        "a producer before `start_step` reused from an audited upstream project: its prior "
        "raw step/tool/parameter matrix remains unchanged for lineage, while materialization "
        "sets the resolved step to disabled and prevents execution.",
        "",
        "Route-changing site/search behavior must be an explicit registered step "
        "parameter. It must not be inferred from a holo/apo label, artifact presence, "
        "PDB identifier, or another scientific fact inside the runner.",
        "",
    ]
    strategy_params = catalog["selection_contract"]["strategy_params"]
    if strategy_params:
        lines += [
            "Registered selection strategies:",
            "",
            "| step | param | type | default | choices |",
            "|---|---|---|---|---|",
        ]
        for param in strategy_params:
            choices = ", ".join(map(str, param.get("choices", []))) or "—"
            default = "**required**" if param.get("required") else f"`{param.get('default')}`"
            lines.append(
                f"| `{param['step']}` | `{param['name']}` | {param['type']} | "
                f"{default} | {_md_cell(choices)} |"
            )
        lines.append("")
    lines += [
        "",
        "## `input_artifacts` (feed a step without running its producer)",
        "",
        "For iteration: inject a modified/new file so a downstream step consumes it "
        "directly (e.g. a modified library SDF for a re-dock round). Injectable "
        "artifacts:",
        "",
        "`" + "`, `".join(catalog["injectable_artifacts"]) + "`",
        "",
        "The producer of an injected artifact must be disabled (often via "
        "`start_step`); otherwise validation fails instead of silently "
        "overwriting the injected generation.",
        "",
        "## Artifacts (inter-step wiring vocabulary)",
        "",
        "Inter-step outputs are wired by these canonical names rather than copied into "
        "tool params. External paths enter only through registered context or "
        "`input_artifacts` fields.",
        "",
        "`" + "`, `".join(catalog["artifacts"]) + "`",
        "",
        "## Two-level metric gate and weighted scoring",
        "",
        "The gate has two levels. Level 1 marks one metric on one "
        "molecule: a value beyond a hard bound is `fail`, a value outside the "
        "soft interval but inside the hard interval is `flag`, and a value "
        "inside the inclusive soft interval is `pass`. Missing sides are "
        "unbounded and a missing value is `not_applicable`. A mark is never a "
        "molecule decision.",
        "",
        "Level 2 folds the marks of each parent metric category and is the only "
        "place a molecule is decided: `V1` (marking) never rejects, `V2` "
        "(standard, the default) rejects on `count(fail) >= 3` or "
        "`ratio(fail) >= 50%` inside one category, and `V3` (strict) rejects on "
        "`count(fail) >= 2` or `ratio(fail) >= 30%`. The count and ratio "
        "branches are independent alternatives and neither carries a minimum "
        "fail count, so a lone `fail` among few labeled metrics is ratio 100% "
        "and rejects. Missing values leave the ratio; no observed-coverage "
        "guard modifies that formula. Any tripped category gives `reject`; "
        "otherwise every evaluated category needs strictly more than half of "
        "its labeled metrics to be `pass` for `keep`. Other labeled cases give "
        "`modify`, and nothing marked gives `not_evaluated`. "
        "Rules whose `effect` is `protocol` or `pose_export` are deliberate "
        "single-level exceptions: they invalidate a run or a pose rather than "
        "folding into a molecule verdict.",
        "",
        "Records the fold rejected leave every percentile pool. Each applicable "
        "metric among the survivors forms its own "
        "pool inside one scientific cohort and uses "
        f"`{catalog['metric_scoring']['normalization_method']}`; category "
        "scores are percentile-ranked again before StepScore/OverallScore. "
        "Bounded zero-weight metrics still mark but contribute exactly zero; "
        "finite metrics with no declared bound are score-only and add no pass "
        "to the gate denominator.",
        "",
        f"Metric rules version: `{catalog['metric_scoring']['rules_version']}`.",
        "",
        "| metric_id | tool | direction | hard lower | soft lower | soft upper | hard upper | weight magnitude | signed weight | gate? | effect |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|:---:|---|",
    ]
    for rule in catalog["metric_scoring"]["rules"]:
        lines.append(
            f"| `{rule['metric_id']}` | {rule['tool']} | "
            f"{rule.get('direction') or '—'} | "
            f"{_display_bound(rule.get('hard_lower'))} | "
            f"{_display_bound(rule.get('soft_lower'))} | "
            f"{_display_bound(rule.get('soft_upper'))} | "
            f"{_display_bound(rule.get('hard_upper'))} | "
            f"{rule.get('weight', 0)} | {rule.get('signed_weight', 0)} | "
            f"{'yes' if rule.get('gate_participates') else 'no'} | "
            f"{rule.get('molecule_effect', 'gate')} |"
        )
    lines += [
        "",
        "## Steps",
        "",
        "| order | step_id | skill | default | consumes | optional consumes | produces |",
        "|---:|---|---|:---:|---|---|---|",
    ]
    for s in catalog["steps"]:
        on = "on" if s["default_enabled"] else "off"
        if s.get("best_effort"):
            on += " (best-effort)"
        lines.append(
            f"| {s['order']} | `{s['id']}` | {s['skill']} | {on} | "
            f"{_join(s['consumes'])} | {_join(s['optional_consumes'])} | "
            f"{_join(s['produces'])} |"
        )
    lines.append("")

    for s in catalog["steps"]:
        lines.append(f"### `{s['id']}` — {s['skill']}")
        lines.append("")
        constraints = []
        for group in s.get("one_of", []):
            constraints.append("at least one of " + _join(group) + " must be enabled")
        for group in s.get("exactly_one", []):
            constraints.append("exactly one of " + _join(group) + " must be enabled")
        if s.get("requires_all"):
            constraints.append("requires all artifacts " + _join(s["requires_all"]))
        for group in s.get("requires_any", []):
            constraints.append("requires at least one artifact from " + _join(group))
        for tool_id, groups in (s.get("tool_one_of") or {}).items():
            for group in groups:
                constraints.append(f"`{tool_id}` needs one of " + _join(group))
        if constraints:
            lines.append("Constraints: " + "; ".join(constraints) + ".")
            lines.append("")
        lines.append(
            "| tool | default | runtime | purpose | when to use | needs | produces |"
        )
        lines.append("|---|:---:|---|---|---|---|---|")
        for t in s["tools"]:
            needs = list(t.get("depends_on", []))
            needs += [
                f"artifact:{a}" for a in t.get("requires_artifacts", [])
            ]
            needs += [
                f"supplied-artifact:{a}"
                for a in t.get("requires_supplied_artifacts", [])
            ]
            if t.get("context_one_of"):
                needs.append("context:" + "/".join(t["context_one_of"]))
            lines.append(
                f"| `{t['id']}` | {'on' if t['default_enabled'] else 'off'} | "
                f"{_runtime_label(t.get('runtime'))} | {t['purpose']} | "
                f"{t['when_to_use']} | {_join(needs)} | "
                f"{_join(t.get('produces_artifacts', []))} |"
            )
        lines.append("")
        param_rows = _param_rows(s)
        if param_rows:
            lines.append("Parameters:")
            lines.append("")
            lines.append("| param | scope | type | default | range/choices |")
            lines.append("|---|---|---|---|---|")
            lines.extend(param_rows)
            lines.append("")
    return "\n".join(lines) + "\n"


def render_selection_policy(catalog: dict) -> str:
    """Render the generated decision policy without duplicating config syntax."""
    lines = [
        "# SBDD Tool Selection Policy (generated)",
        "",
        "> Generated by `pipeline-entry/build_tool_catalog.py` from the current "
        "registry and validator contract metadata. It contains Agent selection guidance. "
        "**Do not edit by hand.** Change the machine rules and rerun the generator.",
        "",
        "This document answers one question: once the facts, user objectives, and resource constraints have been extracted, which registered tools should this round enable or disable? `QUERY_TO_CONFIG.md` only assembles confirmed facts and selection decisions into a config; `TOOL_CATALOG.md` provides the complete inventory of parameters and capabilities; `shared/pipeline_validation.py` enforces the final hard constraints.",
        "",
        "## Selection order",
        "",
        "1. First satisfy the validator's hard constraints on safety, capability, dependencies, and artifacts.",
        "2. Prioritize registered tools explicitly requested or excluded by the user. If that conflicts with a hard constraint, stop and explain; do not silently rewrite the request.",
        "3. Exclude routes whose input, context, or upstream-artifact prerequisites are unmet. Scientific facts inform selection; they must not make the runner rewrite the tool matrix automatically.",
        "4. Match `when_to_use` in the tables below to the task, structural state, evidence requirements, hardware, and budget.",
        "5. Select the smallest sufficient tool set that satisfies the requested outputs and hard constraints. Every additional evidence route needs an explicit purpose or applicability rationale.",
        "6. Without a decisive signal, leave the choice unresolved and request clarification, or select the smallest route with a clear applicability rationale. Do not treat a registry default as a choice already made by the Agent.",
        "7. If multiple candidates remain in an exclusive group, resolve them by explicit user choice first, then by more specific applicability. Compare multiple exclusive routes in separate configs/campaigns.",
        "8. Keep registry parameter defaults unless an explicit need justifies an override. Override only for user constraints, observable resource budgets, or hard cross-tool contracts.",
        "",
        "Both selection decisions and the submitted config must contain a complete explicit matrix: `enabled` for every registered step and `enabled=true/false` for every registered tool under each step. Explicitly set every tool under a disabled step to false. The sole exception is an upstream producer reused before `start_step`: preserve its original step/tool/parameter matrix from the previous round, while `start_step` disables the resolved step during materialization. Those switches record lineage only and do not execute. Parameters may still specify only overrides relative to schema defaults.",
        "",
        "## Unified memory decisions before the main experiment",
        "",
        "Before finalizing the tool matrix, follow the [unified memory protocol](references/benchmark-memory.md) to check pocket, generation and docking memory applicability, prepare the required inputs, or record explicit inapplicability/disablement/fallback reasons. Unknown-site pocket selection uses protein memory within the main Step 1 when available; fixed DeepPocket#1 local-rank fusion is the unavailable-memory or explicitly disabled fallback. Record policies and actual decisions together in memory_consultation.md; do not mix supervised labels.",
        "",
        "Docking retains context.docking_memory_mode choices ligand_protein (default), ligand_only, legacy, and off. select-docking --task redocking uses Top-1 RMSD <2 Å with GNINA; --task pose_generation uses rDock with Top-1 PoseBusters labels and shared prior/margin settings. Off returns the selected task's baseline without memory IO. Generation select-generation uses separate best-tool labels for five metrics and a Pocket2Mol baseline: all objectives use prior strength 10 and switch margin 0.03; priority chooses the required recommendation without changing scores or parameters. For unknown sites, use plan-pocket-tools before Gate A: retrieve five protein groups, retain their five-entry best-detector multiset, then execute unique recommended tools plus the DeepPocket anchor. The registry supports DeepPocket, fpocket, P2Rank, AutoSite and pyKVFinder. Fuse all anchor-matched native Top-3 pockets by DCA reliability, rank and detector family; P2Rank#1 replaces DeepPocket#1 only when recommended and all six guards pass. Preserve anchor point5 residues. auto records unavailable-memory fallback; off uses the fixed four-tool DeepPocket#1 rank/family policy; required rejects unavailable memory. Invalid candidates or provenance fail. Explicit boxes/reviewed sites take precedence. Consult branch references for formulas, inference notes and review commands.",
        "",
        "Explicit user choices and validator hard constraints take precedence. The Agent first confirms runtimes, input protocols, and resources, then uses --available-tools to restrict generation/docking selection to the genuinely available set; registration alone does not mean a tool can execute. Save recommendations, actual selections, fallback/override reasons, and versions in the consultation record. Pocket memory inputs and presets use the registered pocket_memory_* step parameters; do not invent config fields for other selector options. Validate the final complete tool matrix and site/search strategy and present them at Gate A. Neither selectors nor the runner rewrite configs automatically; recommending multiple generators does not give each a full candidate budget.",
        "",
        "Ligand models handle one valid query: do not use the first molecule of a screening library to represent the entire library, or extract features from its true bound pose. Do not force a generation model when a complete candidate library already exists, a valid conditioning ligand is unavailable, or the task only modifies molecules through reactions. Keep evaluation snapshots frozen. In ordinary memory-enabled research, append evaluated observations with append-research-run, then use recall-research before the next plan and record how those observations inform tools and parameters. A report alone is not an ingested memory update. Historical scores are uncalibrated and do not replace current quality gates. See the unified protocol for ablation session/directory isolation and data boundaries.",
        "",
        "## Candidate-source fallback when ChEMBL is unavailable",
        "",
        "For candidate-discovery tasks, record a source failure when ChEMBL retrieval cannot finish because of a read timeout, connection error, or HTTP service error. Continue with existing usable candidates if available; otherwise transfer the original ChEMBL candidate budget to Pocket2Mol, explicitly setting `chembl.enabled=false` and `pocket2mol.enabled=true`. A ChEMBL-only plan then becomes entirely Pocket2Mol-generated. Respect the total molecule limit rather than allocating a full budget to each source.",
        "",
        "The Agent updates, validates, and executes this fallback plan; the runner does not enable tools implicitly. Account for one fallback within the existing authorization and attempt budgets, retain failure evidence, and reuse valid upstream artifacts. Stop explicitly if Pocket2Mol lacks the required structure/runtime or also fails. If the user requests only known ACT compounds, a complete ACT/DEC benchmark, or excludes generation, do not present GEN as ACT to satisfy the original objective; do not replace an existing local library either. See the [molecule-generator ChEMBL fallback procedure](../molecule-generator/SKILL.md#chembl-unavailable-pocket2mol-recovery) for switch, count, recovery, and reporting requirements.",
        "",
        "## Explicit strategies without implicit holo/apo routing",
        "",
        "The policy must explicitly select `pocket_analysis.params.site_selection_strategy`: `ligand_site`, `known_site`, `detector_consensus`, or `provided_box`. It must also explicitly select `docking.params.search_space_source`: `reference_ligand`, `docking_box`, `whole_receptor`, or `precomputed_grid`. Use only choices currently registered in the catalog. If a required artifact is missing, change the plan or stop; the runner does not change values automatically based on artifact availability.",
        "",
        "`context.scientific_facts.has_metal`, `is_covalent`, and `protonation_sensitive` are scientific inputs that inform Agent selection and validator capability/safety checks; they do not automatically enable or disable tools.",
        "",
        "### Choosing `site_selection_strategy`",
        "",
        "| choice | policy condition | required alignment |",
        "|---|---|---|",
        "| `provided_box` | The user supplies a fixed box whose region must be analyzed and validated in the same round | Retain `input_artifacts.docking_box`; Step 1 does not generate a replacement box; all detectors may be disabled; docking uses `search_space_source=docking_box` |",
        "| `ligand_site` | A usable reference ligand has been confirmed and is explicitly selected to define the primary site for this round | Enable `pocket_analysis.ligand_site` and ensure the `reference_ligand` artifact is available |",
        "| `known_site` | Independently reviewed `context.known_site_residues` are available and the objective is to define the search region from that site | Supply reviewed residues/provenance; geometric tools may provide evidence but must not replace this primary strategy |",
        "| `detector_consensus` | No reference ligand/reviewed site can serve as the primary definition, or the user explicitly requests cavity-detector consensus | Set the complete tool matrix from plan-pocket-tools.execution_tools; fuse all matching native Top-3 pockets under memory, using distance ≤5 Å and point5 IoU ≥0.30. Only unavailable/off memory uses the four-tool fixed DeepPocket#1 fallback with at most one match per other tool |",
        "",
        "### Choosing `search_space_source`",
        "",
        "| choice | policy condition | required alignment |",
        "|---|---|---|",
        "| `reference_ligand` | The user requests reference-guided/autobox search and the reference ligand has been reviewed as usable | Supply `reference_ligand` and select a docking engine supporting that input |",
        "| `docking_box` | Use an explicit center/size from pocket analysis or external review | Supply the `docking_box` artifact; the runner must not change this to autobox |",
        "| `whole_receptor` | The user explicitly requests blind/whole-receptor search and the selected engine and budget support it | Record the high-cost/low-focus limitations; a missing box must not automatically activate this mode |",
        "| `precomputed_grid` | The user supplies and requests reuse of an engine-specific grid/map bundle | Explicitly inject every required grid artifact for that engine and pass validation |",
        "",
        "If several choices apply, prioritize the explicit user request. Otherwise choose the option most directly aligned with the primary site evidence and selected engine input contract. If any prerequisite is unmet, stop or submit a new plan; the runner must not switch to another choice automatically.",
        "",
        "## Pose evaluation scope and cost",
        "",
        "The default `pose_evaluation.params.selection_strategy=ranked_per_compound` selects poses by docking `engine_pose_rank` (or existing scores when rank is absent). The per-molecule evaluation limit depends on the number of molecules in this round: 20 for at most 50 molecules, 5 for 51–1000, and 1 for more than 1000. Evaluate every selected pose; do not use the top 1 verdict to decide whether to evaluate the next pose.",
        "",
        "`poses_per_compound` can further tighten the per-molecule limit; `max_poses` is the actual evaluation budget. Return a resumable partial result if not all selected poses are evaluated. Tool errors and NA remain unresolved and do not count as scientific failures. The legacy `first_passing_per_compound` is retained only for explicit compatibility use; `all` supports full-ensemble comparison.",
        "",
        "`pose_evaluation.tools.posecheck.params.num_conformers` defaults to `10`, and reference-ligand calibration uses the same count. This parameter controls reference-conformer sampling for strain, independently of docking `num_modes`. Reducing sampling reduces coverage of low-energy reference conformers; keep it fixed within a cohort and retain the parameter record.",
        "",
        "`mdanalysis` and `mmgbsa` can change pose verdicts during postprocessing. The default strategy evaluates all selected poses and supports both tools. The legacy `first_passing_per_compound` remains incompatible with them and is rejected by the validator.",
        "",
        "## No silent fallback",
        "",
        "Disabled tools receive no runtime preflight and are not called. A missing or broken runtime for an enabled tool fails the current plan, except for record-only ProLIF, whose runtime failure is recorded in `optional_missing` and allows execution to continue. Propose alternative tools only in a new plan version, with renewed validation and confirmation. In particular, never automatically replace a tool explicitly named by the user.",
        "",
        "",
        "## Tool selection tables by step",
        "",
        "`default` records registry recommendation metadata for policy comparison only; it neither participates in execution nor represents the Agent's selection. `when_to_use`, dependencies, artifacts, and exclusivity relationships all come from the current registry.",
        "",
    ]
    for step in catalog["steps"]:
        lines.append(f"### `{step['id']}`")
        lines.append("")
        constraints = _step_selection_constraints(step)
        if constraints:
            lines.append("Declarative hard constraints for this step: " + "; ".join(constraints) + ".")
            lines.append("")
        lines.append("| canonical tool | default | when to use | hard prerequisites / group |")
        lines.append("|---|:---:|---|---|")
        for tool in step["tools"]:
            lines.append(
                f"| `{step['id']}.{tool['id']}` | "
                f"{'on' if tool['default_enabled'] else 'off'} | "
                f"{_md_cell(tool.get('when_to_use') or tool.get('purpose') or '—')} | "
                f"{_md_cell(_tool_selection_requirements(step, tool))} |"
            )
        lines.append("")

    lines += [
        "## Validator hard stops",
        "",
        "The following summaries are derived from `shared/pipeline_validation.py:VALIDATION_RULES`. The policy cannot relax them; validator code and dry-run results define exact behavior and errors.",
        "",
        "| contract | enforced rule |",
        "|---|---|",
    ]
    for rule_id, description in catalog["validation_rules"].items():
        lines.append(f"| `{rule_id}` | {_md_cell(description)} |")
    lines += [
        "",
        "## Agent output contract",
        "",
        "After selection, provide the config-assembly stage with:",
        "",
        "- An explicit enabled decision for every registered step;",
        "- The complete explicit on/off matrix for every registered tool under every step;",
        "- Explicit `site_selection_strategy` and `search_space_source` choices with reasons;",
        "- Parameter overrides and the user facts or resource evidence supporting them;",
        "- Hard-constraint or applicability reasons for excluded routes;",
        "- Unresolved facts that could change the route.",
        "",
        "Do not invent unregistered tools, environment paths, artifact names, or config fields in the selection result. Only the canonical `step.tool` entries in this document's tables are currently selectable switches; distribution names, internal helpers, and labels from other delivery packages are not registry tools.",
        "",
    ]
    return "\n".join(lines)


def _step_selection_constraints(step: dict) -> list[str]:
    constraints = []
    for group in step.get("one_of", []):
        constraints.append("at least one of " + _join(group))
    for group in step.get("exactly_one", []):
        constraints.append("exactly one of " + _join(group))
    if step.get("requires_all"):
        constraints.append("requires artifacts " + _join(step["requires_all"]))
    for group in step.get("requires_any", []):
        constraints.append("requires one artifact from " + _join(group))
    for tool_id, groups in (step.get("tool_one_of") or {}).items():
        for group in groups:
            constraints.append(f"`{tool_id}` needs one of " + _join(group))
    return constraints


def _tool_selection_requirements(step: dict, tool: dict) -> str:
    requirements = []
    if tool.get("depends_on"):
        requirements.append(
            "depends on "
            + ", ".join(f"`{step['id']}.{item}`" for item in tool["depends_on"])
        )
    if tool.get("requires_artifacts"):
        requirements.append("artifacts " + _join(tool["requires_artifacts"]))
    if tool.get("requires_supplied_artifacts"):
        requirements.append(
            "explicitly supplied " + _join(tool["requires_supplied_artifacts"])
        )
    if tool.get("context_one_of"):
        requirements.append("context one of " + _join(tool["context_one_of"]))
    for group in step.get("one_of", []):
        if tool["id"] in group:
            requirements.append("member of at-least-one group")
    for group in step.get("exactly_one", []):
        if tool["id"] in group:
            requirements.append("member of exactly-one group")
    return "; ".join(requirements) if requirements else "—"


def _md_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _display_bound(value: object) -> str:
    return "—" if value is None else str(value)


def _param_rows(step: dict) -> list[str]:
    rows = []
    for p in step.get("params", []):
        rows.append(_param_row(p, "step"))
    for t in step["tools"]:
        for p in t.get("params", []):
            rows.append(_param_row(p, f"{t['id']}"))
    return rows


def _param_row(p: dict, scope: str) -> str:
    rng = ""
    if "choices" in p:
        rng = "one of " + ", ".join(str(c) for c in p["choices"])
    elif "item_choices" in p:
        rng = "items from " + ", ".join(str(c) for c in p["item_choices"])
    elif "min" in p or "max" in p:
        rng = f"[{p.get('min', '')}, {p.get('max', '')}]"
    default = "**required**" if p.get("required") else f"`{p.get('default')}`"
    return f"| `{p['name']}` | {scope} | {p['type']} | {default} | {rng} |"


def _schema_row(p: dict) -> str:
    rng = ""
    if "choices" in p:
        rng = "one of " + ", ".join(str(c) for c in p["choices"])
    elif "item_choices" in p:
        rng = "items from " + ", ".join(str(c) for c in p["item_choices"])
    elif "min" in p or "max" in p:
        rng = f"[{p.get('min', '')}, {p.get('max', '')}]"
    default = "—" if p.get("required") else f"`{p.get('default')}`"
    return f"| `{p['name']}` | {p['type']} | {default} | {rng} |"


def _join(items) -> str:
    return ", ".join(f"`{i}`" for i in items) if items else "—"


def _runtime_label(runtime: dict | None) -> str:
    if not runtime:
        return "—"
    if runtime.get("kind") == "composite":
        labels = []
        for name, component in (runtime.get("components") or {}).items():
            env = component.get("environment", "current")
            target = (
                component.get("command")
                or component.get("module")
                or component.get("kind")
            )
            labels.append(f"{name}={env}:{target}")
        return "`" + "; ".join(labels) + "`" if labels else "`composite`"
    env = runtime.get("environment", "current")
    target = runtime.get("command") or runtime.get("module") or runtime.get("kind")
    return f"`{env}:{target}`"



def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Agent tool catalog and selection policy from machine rules"
    )
    parser.add_argument("--check", action="store_true",
                        help="Verify committed catalog matches the registry; exit 1 on drift")
    args = parser.parse_args()

    catalog = build_catalog()
    json_text = json.dumps(catalog, indent=2, ensure_ascii=False) + "\n"
    md_text = render_markdown(catalog)
    policy_text = render_selection_policy(catalog)

    if args.check:
        drift = []
        if not JSON_PATH.exists() or JSON_PATH.read_text() != json_text:
            drift.append(str(JSON_PATH))
        if not MD_PATH.exists() or MD_PATH.read_text() != md_text:
            drift.append(str(MD_PATH))
        if not POLICY_PATH.exists() or POLICY_PATH.read_text() != policy_text:
            drift.append(str(POLICY_PATH))
        if drift:
            print("catalog OUT OF SYNC with registry: " + ", ".join(drift), file=sys.stderr)
            print("run: python pipeline-entry/build_tool_catalog.py", file=sys.stderr)
            sys.exit(1)
        print("tool catalog and selection policy in sync with machine rules")
        return

    JSON_PATH.write_text(json_text)
    MD_PATH.write_text(md_text)
    POLICY_PATH.write_text(policy_text)
    print(f"wrote {JSON_PATH}")
    print(f"wrote {MD_PATH}")
    print(f"wrote {POLICY_PATH}")


if __name__ == "__main__":
    main()

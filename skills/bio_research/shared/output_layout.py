"""Canonical publication layout for one scored run.

The pipeline's working directories stay as they are: every skill keeps writing
its own ``step{N}/`` tree, because that is what the orchestrator wires together.
This module defines the isolated *published* view beside that tree -- eight
reader-facing workflow directories, the flat ``n x m`` scoring tables, the
numbered category directories and the mandated filenames that a reader is
allowed to depend on.

Two conventions supplement the tracked output-layout contract:

* **Reference, never copy.**  ``source_map.json`` and ``input_refs.json`` carry
  paths relative to the publication directory to the artifacts a table was
  derived from.  Poses,
  trajectories and topologies are large and already written once; duplicating
  them into ``02_tool_runs/`` would multiply a run's footprint for no new
  information.  The contract defines ``source_map.json`` as
  "tool_run_id, metric columns and relative paths".
* **An empty cell means "not measured".**  Never ``0``, never ``NA`` as a
  string, never a sentinel.  A missing metric leaves both sides of every ratio
  it participates in, so writing a number here would be a silent lie.

Round and backtracking artifacts (``checkpoints/round_snapshot_manifest.json``)
are deliberately absent: that logic does not exist in this pipeline yet.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from metric_scoring import CATEGORY_CONFIG, STEP_IDS


#: Written into every cell that has no measured value.
NA = ""

#: Reader-facing output is kept separate from the resumable working tree.
PUBLISHED_DIR = "published"

#: Tracked, skill-local description of the output contract implemented here.
LAYOUT_CONTRACT = "reflective-coordinator/references/output-layout.md"
LAYOUT_VERSION = 3

# ---------------------------------------------------------------------------
# Mandated filenames.  The output contract fixes these patterns by name; the
# rest are fixed here so a reader never has to guess.
# ---------------------------------------------------------------------------

ALL_METRICS_TABLE = "all_metrics_table.tsv"
FINAL_REPORT = "final_report.tsv"
FINAL_REPORT_FULL = "final_report_full.tsv"
METRIC_CATALOG = "metric_catalog_and_scoring.md"
RUN_MANIFEST = "run_manifest.json"
FAILURE_REGISTRY = "failure_registry.jsonl"
LINEAGE_GRAPH = "lineage_graph.json"

#: The publication is ordered by scientific workflow, independently of the
#: legacy runner's ``step0`` ... ``step7`` working-directory names.
INPUTS_DIR = "00_inputs"
STRUCTURE_PREP_DIR = "01_structure_prep"
POCKET_ANALYSIS_DIR = "02_pocket_analysis"
MOLECULE_GENERATION_DIR = "03_molecule_generation"
SM_PROPERTY_DIR = "04_sm_property"
DOCKING_DIR = "05_docking"
POST_DOCKING_DIR = "06_post_docking"
FINAL_REPORT_DIR = "07_final_report"

PUBLICATION_STAGE_DIRS: tuple[str, ...] = (
    INPUTS_DIR,
    STRUCTURE_PREP_DIR,
    POCKET_ANALYSIS_DIR,
    MOLECULE_GENERATION_DIR,
    SM_PROPERTY_DIR,
    DOCKING_DIR,
    POST_DOCKING_DIR,
    FINAL_REPORT_DIR,
)


def final_artifact(name: str) -> str:
    """Publication-relative path for a global/final artifact."""
    return f"{FINAL_REPORT_DIR}/{name}"


#: Files that make a publication usable even though the manifest cannot hash
#: itself (or the lineage file written immediately after it).
MANDATORY_PUBLICATION_FILES: tuple[str, ...] = (
    final_artifact(ALL_METRICS_TABLE),
    final_artifact(FINAL_REPORT),
    final_artifact(FINAL_REPORT_FULL),
    final_artifact(METRIC_CATALOG),
    final_artifact(RUN_MANIFEST),
    final_artifact(LINEAGE_GRAPH),
)

# Partial views have evidence indexes and an execution snapshot, not rankings.
PARTIAL_PUBLICATION_FILES: tuple[str, ...] = (
    final_artifact(RUN_MANIFEST),
    final_artifact(LINEAGE_GRAPH),
    final_artifact("publication_status.json"),
    final_artifact("README.md"),
    final_artifact("final_report.md"),
    final_artifact("final_report.html"),
    final_artifact("final_report_data.json"),
    *(f"{stage}/artifact_index.tsv" for stage in PUBLICATION_STAGE_DIRS),
    "00_inputs/input_manifest.json",
    *(f"{stage}/step_manifest.json" for stage in PUBLICATION_STAGE_DIRS[1:]),
)

STEP_TABLE = "step_table.tsv"
STEP_REPORT = "step_report.tsv"
STEP_MANIFEST = "step_manifest.json"
INPUT_REFS = "input_refs.json"

CATEGORY_TABLE = "category_table.tsv"
CATEGORY_REPORT = "category_report.tsv"
SOURCE_MAP = "source_map.json"

TOOL_TABLE = "tool_table.tsv"
ASSET_MANIFEST = "asset_manifest.json"
ARTIFACT_INDEX = "artifact_index.tsv"
INPUT_MANIFEST = "input_manifest.json"

SHARED_ASSETS_DIR = "00_shared_assets"
CATEGORIES_DIR = "01_categories"
TOOL_RUNS_DIR = "02_tool_runs"
INTERNAL_DIR = "_internal"

#: ``_internal/`` is standardised across every tool run.
TOOL_SUMMARY = "tool_summary.json"
INTERNAL_INDEX = "internal_index.tsv"
RAW_ARCHIVE = "raw_response_or_native.tar.zst"
RUN_LOG = "run.log"
COMMAND_AND_CONFIG = "command_and_config.json"
PROVENANCE = "provenance.json"

#: Scoring-policy step IDs map onto reader-facing workflow directories.  These
#: IDs and their table-column semantics remain stable even though the physical
#: names no longer expose them.
STEP_DIR: dict[str, str] = {
    "step4": SM_PROPERTY_DIR,
    "step5": DOCKING_DIR,
    "step6": POST_DOCKING_DIR,
}

#: Category code -> numbered directory, in the contract's report order.  The
#: number is part of the name so a directory listing sorts into report order.
CATEGORY_DIR: dict[str, str] = {
    # 04_sm_property
    "PHY": "01_molecular_properties",
    "ALERT": "02_structural_alerts",
    "SYN": "03_synthesizability",
    "ABS": "04_admet_absorption",
    "DIST": "05_admet_distribution",
    "MET": "06_admet_metabolism",
    "EXCR": "07_admet_excretion",
    "TOX": "08_admet_toxicity",
    # 05_docking
    "AFF": "01_docking_affinity",
    "DRAW": "02_docking_raw",
    "FLEX": "03_flexibility",
    "BENCH": "04_benchmark_controls",
    # 06_post_docking
    "IDQC": "01_identity_qc",
    "TOPO": "02_candidate_topology",
    "PBQ": "03_post_docking_summary",
    "CHEM": "04_chemical_validity",
    "INTRA": "05_pose_intramolecular_validity",
    "INTER": "06_pose_intermolecular_validity",
    "DQC": "07_docking_qc",
    "FP": "08_pose_fingerprint",
    "INTX": "09_interaction_information",
    "POCKET": "10_pose_pocket_validity",
    "DYN": "11_dynamics_energy",
}

#: Report order within a step, derived from the directory numbering above.
STEP_CATEGORY_ORDER: dict[str, tuple[str, ...]] = {
    step_id: tuple(
        code for code in sorted(
            (c for c in CATEGORY_DIR if CATEGORY_CONFIG[c]["step_id"] == step_id),
            key=lambda c: CATEGORY_DIR[c],
        )
    )
    for step_id in STEP_DIR
}

#: Every category in report order, across every step.
CATEGORY_ORDER: tuple[str, ...] = tuple(
    code for step_id in STEP_IDS for code in STEP_CATEGORY_ORDER[step_id]
)

#: Pipeline tool switch (the ``plan_tool`` of a metric rule, or a tool that
#: produces artifacts but no metric) -> ``02_tool_runs/`` directory name.
#: The contract's spelling wins wherever it names the tool; the rest follow the same
#: lowercase-with-underscores shape.
TOOL_DIR: dict[str, str] = {
    # Step 4 -- the RDKit facilities are split the way the contract splits them:
    # descriptor-style numbers in one archive, alert catalogs in the other.
    "lipinski": "rdkit_descriptors",
    "qed": "rdkit_descriptors",
    "sa_score": "rdkit_descriptors",
    "pains": "rdkit",
    "chemfh": "chemfh",
    "aizynthfinder": "aizynthfinder",
    "admet": "admet-ai",
    # Step 5
    "autodock_vina": "autodock_vina",
    "unidock": "uni-dock",
    "gnina": "gnina",
    "autodock4": "autodock_4",
    "rdock": "rdock",
    "adfr": "autodockfr",
    "dock6": "dock6",
    "diffdock": "diffdock",
    "redocking": "vina_redocking_validation",
    # Step 6
    "posebusters": "posebusters",
    "posecheck": "posecheck",
    "plip": "plip",
    "prolif": "prolif",
    "arpeggio": "arpeggio",
    "luna": "luna",
    "oddt": "oddt",
    "blueprint_match": "blueprint_match",
    "pocket_occupancy": "pocket_occupancy",
    "candidate_topology": "candidate_topology",
    "rtmscore": "rtmscore",
    "short_md": "openmm_gromacs",
    "amber_mmpbsa": "mmpbsa",
    "openfe_rbfe": "openfe_rbfe",
    "ties_rbfe": "ties",
}


def tool_dir_name(tool_id: str) -> str:
    """Directory name for a tool, falling back to its own id."""
    return TOOL_DIR.get(tool_id, str(tool_id))


def slug_tool_id(name: Any) -> str:
    """Fold a display name into an id: ``RDKit descriptors`` -> ``rdkit_descriptors``."""
    slug = "".join(
        char if char.isalnum() else "_" for char in str(name or "").strip().lower()
    )
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def rule_tool_id(rule: Mapping[str, Any]) -> str:
    """The tool id whose directory owns a metric rule.

    Most rules name their planner tool directly.  The RDKit descriptor block has
    no planner entry of its own -- it is computed inline -- so its display name
    is folded into an id instead of inventing a placeholder directory.
    """
    return str(
        rule.get("plan_tool") or slug_tool_id(rule.get("tool")) or "unknown"
    )

# ---------------------------------------------------------------------------
# Table headers.
# ---------------------------------------------------------------------------

#: The 19 columns displayed in candidate tables.
FINAL_REPORT_COLUMNS: tuple[str, ...] = (
    "source_accession", "pocket_id", "ligand_id", "canonical_SMILES", "pose_id",
    "overall_score", "overall_rank", "decision",
    "step4_score", "step4_rank", "step5_score", "step5_rank",
    "step6_score", "step6_rank",
    "pose_file", "session_id", "round_id", "cohort_id", "generated_at",
)

_FULL_IDENTITY_COLUMNS: tuple[str, ...] = (
    # run and comparability
    "session_id", "round_id", "cohort_id", "profile_id",
    "metric_catalog_version", "scoring_schema_version", "generated_at",
    # ranking partition: ``overall_rank`` is local to one cohort *and* one
    # normalization pool, so the pool has to travel with the rank for a reader
    # to know which "rank 1" they are looking at. ``benchmark_decoys`` rows are
    # never nomination candidates.
    "normalization_pool", "rank_partition_id", "rank_scope",
    # protein identity
    "protein_id", "source_type", "source_accession", "chain_ids",
    "sequence_hash", "structure_hash", "protein_state_id", "preparation_version",
    # pocket identity
    "pocket_id", "pocket_method", "residue_ids", "center_xyz_A", "box_xyz_A",
    "pocket_volume_A3", "pocket_rank", "pocket_version",
    # ligand identity
    "ligand_id", "parent_ligand_id", "molecule_state_id", "canonical_SMILES",
    "InChIKey", "tautomer_id", "protomer_id", "stereo_id",
    # pose identity
    "pose_id", "engine", "engine_version", "score_function", "engine_pose_rank",
    "seed", "receptor_conformer_id", "pose_hash", "pose_file",
    "selected_pose_index", "source_pose_index",
)

_FULL_GATE_COLUMNS: tuple[str, ...] = (
    "gate_status", "hard_fail_codes", "soft_alert_codes",
    "planned_weighted_coverage_overall", "category_weighted_coverage",
    # The catalog-referenced ratio is a pipeline-local audit column kept visible
    # beside the planned one; only the planned ratio the spec defines at
    # The executable policy decides a release band.
    "catalog_category_weighted_coverage",
    "confidence_grade", "missing_required_metrics", "missing_active_steps",
    "rank_status", "ranking_blocked", "ranking_blocked_by",
)

_FULL_STEP_COLUMNS: tuple[str, ...] = tuple(
    f"{step_id}_{suffix}"
    for step_id in STEP_IDS for suffix in ("score", "rank", "coverage")
)

#: The executable policy: a production table carries four columns per category.
#: The 8.3 sample merges score and rank into one cell purely to fit the page.
CATEGORY_COLUMN_SUFFIXES: tuple[str, ...] = (
    "score", "rank", "status", "coverage",
)

_FULL_CATEGORY_COLUMNS: tuple[str, ...] = tuple(
    f"{code}_{suffix}"
    for code in CATEGORY_ORDER for suffix in CATEGORY_COLUMN_SUFFIXES
)

_FULL_RANKING_COLUMNS: tuple[str, ...] = (
    "overall_score", "overall_rank", "decision", "reason_codes",
    "pareto_cluster", "diversity_cluster",
)

_FULL_SUPPLEMENT_COLUMNS: tuple[str, ...] = (
    "uncertainty_grade", "user_priority_snapshot", "fallback_used",
    "artifact_path", "notes",
)

FINAL_REPORT_FULL_COLUMNS: tuple[str, ...] = (
    _FULL_IDENTITY_COLUMNS
    + _FULL_GATE_COLUMNS
    + _FULL_STEP_COLUMNS
    + _FULL_CATEGORY_COLUMNS
    + _FULL_RANKING_COLUMNS
    + _FULL_SUPPLEMENT_COLUMNS
)

#: One row per step per candidate: the score, its rank and the evidence behind
#: it.  The contract marks this the user-facing table of a step.
STEP_REPORT_COLUMNS: tuple[str, ...] = (
    "row_id", "step_id", "step_score", "step_rank",
    "category_weighted_coverage", "planned_category_weighted_coverage",
    "catalog_category_weighted_coverage",
    "available_parent_weight", "reference_parent_weight",
    "available_categories", "missing_reference_categories",
    "rank_status", "gate_status",
)

#: One row per candidate per category.
CATEGORY_REPORT_COLUMNS: tuple[str, ...] = (
    "row_id", "category_code", "category_name", "step_id", "parent_weight",
    "score", "rank", "status", "coverage",
    "available_abs_weight", "planned_abs_weight", "catalog_abs_weight",
    "catalog_implementation_fraction", "evidence_strength",
    "available_metric_ids", "consequential_missing_metric_ids",
    "audit_only_missing_metric_ids",
)

# ---------------------------------------------------------------------------
# Stable row identity and the four-way decision.
# ---------------------------------------------------------------------------

#: Step 4 rows are molecule states; Step 5-6 rows are poses.  A table must join
#: on these, never on line order.
ROW_ID_FIELD: dict[str, str] = {
    "step4": "molecule_state_id",
    "step5": "pose_id",
    "step6": "pose_id",
}


def row_id(record: Mapping[str, Any], step_id: str | None = None) -> str:
    """The stable id for one row of a table at ``step_id``.

    Step 4 is keyed by molecule state, so two poses of the same molecule share a
    row; Steps 5-6 are keyed by pose.  The final tables are per representative
    pose and therefore use the pose id as well.
    """
    if step_id == "step4":
        for key in ("molecule_state_id", "chemical_identity_id", "mol_name"):
            value = record.get(key)
            if value:
                return str(value)
        return NA
    for key in ("pose_id", "chemical_identity_id", "mol_name"):
        value = record.get(key)
        if value:
            return str(value)
    return NA


def decision(metric_gate: Any, rank_status: Any) -> str:
    """Map the three-tier gate and release status onto four reporting verdicts.

    The report uses KEEP/MODIFY/RESERVE/REJECT. RESERVE has no gate of its own:
    the candidate cleared every gate, but its evidence is too thin for a formal
    release, so it is held rather than advanced.
    """
    gate = str(metric_gate or "").lower()
    if gate == "reject":
        return "REJECT"
    if gate == "modify":
        return "MODIFY"
    if gate == "keep":
        return "KEEP" if str(rank_status or "") == "formal" else "RESERVE"
    return NA

# ---------------------------------------------------------------------------
# Atomic writers.  Anything a reader may already be holding open is replaced,
# never truncated in place.
# ---------------------------------------------------------------------------


def _atomic_write(path: Path, payload: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(
        dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp",
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8", newline="") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    return path


def format_cell(value: Any) -> str:
    """Render one table cell.  ``None`` becomes an empty cell, never ``0``."""
    if value is None:
        return NA
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return NA
        # Six places keeps every rank-percentile distinction without pretending
        # to a precision the upstream tools do not have.
        return f"{value:.6f}".rstrip("0").rstrip(".") or "0"
    if isinstance(value, (list, tuple, set)):
        return ";".join(format_cell(item) for item in sorted(map(str, value)))
    if isinstance(value, Mapping):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    text = str(value)
    # A tab or newline inside a cell would break the flat-table contract.
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")

def write_tsv(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> Path:
    """Write one ``n x m`` flat table.  Columns are fixed; rows may be sparse."""
    lines = ["\t".join(columns)]
    for record in rows:
        lines.append("\t".join(format_cell(record.get(name)) for name in columns))
    return _atomic_write(path, "\n".join(lines) + "\n")


def write_json(path: Path, payload: Any) -> Path:
    return _atomic_write(
        path, json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
    )


def write_text(path: Path, payload: str) -> Path:
    return _atomic_write(path, payload if payload.endswith("\n") else payload + "\n")


def append_jsonl(path: Path, records: Iterable[Mapping[str, Any]]) -> Path:
    """Append to a registry.  Existing lines are never rewritten."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        for record in records:
            stream.write(
                json.dumps(record, sort_keys=True, default=str) + "\n",
            )
    return path


def file_digest(path: Path) -> dict[str, Any]:
    """Size and SHA-256 of one artifact, for the manifests and index tables."""
    if not path.is_file():
        return {"exists": False}
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return {
        "exists": True,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def directory_digest(path: Path) -> dict[str, Any]:
    """Content-address an exact directory tree, including member names.

    Directory-valued tool outputs (for example per-batch pose evidence) cannot
    be represented by hashing only the files that happened to exist at publish
    time: an added or removed member must also invalidate the publication.
    Symbolic links and special files are rejected because their targets or byte
    semantics are not stable publication evidence.
    """
    path = Path(path)
    if path.is_symlink() or not path.is_dir():
        return {"exists": False, "kind": "directory"}
    digest = hashlib.sha256()
    total_bytes = 0
    file_count = 0
    for member in sorted(path.rglob("*"), key=lambda item: item.relative_to(path).as_posix()):
        relative = member.relative_to(path).as_posix()
        if member.is_symlink():
            raise ValueError(f"directory evidence contains a symbolic link: {member}")
        if member.is_dir():
            digest.update(b"D\0")
            digest.update(relative.encode("utf-8"))
            digest.update(b"\0")
            continue
        if not member.is_file():
            raise ValueError(f"directory evidence contains a special entry: {member}")
        identity = file_digest(member)
        size = int(identity["bytes"])
        digest.update(b"F\0")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(size).encode("ascii"))
        digest.update(b"\0")
        digest.update(str(identity["sha256"]).encode("ascii"))
        digest.update(b"\0")
        total_bytes += size
        file_count += 1
    return {
        "exists": True,
        "kind": "directory",
        "files": file_count,
        "bytes": total_bytes,
        "sha256": digest.hexdigest(),
    }


def _verify_pipeline_source_snapshot(
    publication: Path,
    manifest: Mapping[str, Any],
) -> int:
    """Verify the live step/result generation sealed by a final publication."""
    finalization = manifest.get("finalization")
    if not isinstance(finalization, Mapping) or finalization.get("status") != "complete":
        raise ValueError("pipeline publication is not marked fully finalized")

    source_root_value = manifest.get("source_root")
    if not isinstance(source_root_value, str) or not source_root_value:
        raise ValueError("pipeline publication has no source_root")
    source_root_path = Path(source_root_value)
    if source_root_path.is_absolute():
        raise ValueError("pipeline publication source_root must be relative")
    source_root = (publication / source_root_path).resolve()
    if source_root != publication.parent.resolve():
        raise ValueError("pipeline publication source_root is not its round directory")
    if source_root.name != str(manifest.get("session_id") or ""):
        raise ValueError("pipeline publication source_root does not match session_id")

    snapshot = manifest.get("source_snapshot")
    if not isinstance(snapshot, Mapping) or not snapshot:
        raise ValueError("pipeline publication has no source snapshot")
    checked = 0
    snapshot_paths: dict[Path, str] = {}
    for raw_relative, raw_expected in snapshot.items():
        relative = Path(str(raw_relative))
        if relative.is_absolute():
            raise ValueError(
                f"pipeline source snapshot path must be relative: {raw_relative}"
            )
        path = (publication / relative).resolve()
        # Explicit input artifacts may live outside the session directory and
        # upstream rounds need not share a direct parent.  Their references are
        # still publication-relative and content-addressed; requiring project
        # containment here would make legitimate start-step runs impossible.
        try:
            path.relative_to(publication)
        except ValueError:
            pass
        else:
            raise ValueError(
                f"pipeline source snapshot points into publication: {raw_relative}"
            )
        if path in snapshot_paths:
            raise ValueError(
                "pipeline source snapshot aliases one file more than once: "
                f"{snapshot_paths[path]}, {raw_relative}"
            )
        snapshot_paths[path] = str(raw_relative)
        if not isinstance(raw_expected, Mapping):
            raise ValueError(
                f"pipeline publication source metadata is invalid: {raw_relative}"
            )
        expected_kind = str(raw_expected.get("kind") or "path")
        if expected_kind == "directory":
            if path.is_symlink() or not path.is_dir():
                raise ValueError(
                    f"pipeline publication source directory is missing: {raw_relative}"
                )
            actual = directory_digest(path)
        elif expected_kind == "path":
            if path.is_symlink() or not path.is_file():
                raise ValueError(
                    f"pipeline publication source is missing: {raw_relative}"
                )
            actual = file_digest(path)
        else:
            raise ValueError(
                f"pipeline publication source kind is invalid: {raw_relative}"
            )
        if raw_expected.get("bytes") != actual.get("bytes"):
            raise ValueError(f"pipeline publication source size changed: {raw_relative}")
        if raw_expected.get("sha256") != actual.get("sha256"):
            raise ValueError(f"pipeline publication source SHA-256 changed: {raw_relative}")
        checked += 1

    final_stage_path = publication / FINAL_REPORT_DIR / STEP_MANIFEST
    try:
        final_stage = json.loads(final_stage_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("pipeline final-stage manifest is invalid") from exc
    if not isinstance(final_stage, Mapping) or final_stage.get("status") != "complete":
        raise ValueError("pipeline final-stage evidence is not complete")
    source_results = final_stage.get("source_results")
    if not isinstance(source_results, list):
        raise ValueError("pipeline final-stage source results are missing")
    statuses = {
        str(item.get("pipeline_step_id")): str(item.get("status"))
        for item in source_results if isinstance(item, Mapping)
    }
    if statuses.get("coordination") != "complete" or statuses.get("reporting") != "complete":
        raise ValueError(
            "pipeline final-stage coordination/reporting evidence is incomplete"
        )
    expected_results = {
        "coordination": source_root / "step6" / "result.json",
        "reporting": source_root / "step7" / "result.json",
    }
    declared_results = {
        str(item.get("pipeline_step_id")): item.get("result_json")
        for item in source_results
        if isinstance(item, Mapping)
    }
    for step_id, expected_path in expected_results.items():
        raw_result = declared_results.get(step_id)
        if not isinstance(raw_result, str) or not raw_result:
            raise ValueError(
                f"pipeline final-stage {step_id} result reference is missing"
            )
        actual_path = (publication / raw_result).resolve()
        if actual_path != expected_path.resolve():
            raise ValueError(
                f"pipeline final-stage {step_id} result reference is invalid"
            )
        if actual_path not in snapshot_paths:
            raise ValueError(
                f"pipeline final-stage {step_id} result is not source-sealed"
            )
    return checked


def verify_publication(
    root: Path | str,
    *,
    expected_session_id: str | None = None,
    require_pipeline_final: bool = False,
) -> dict[str, Any]:
    """Validate a complete publication and reject stale, untracked files.

    Exact membership matters here.  Merely validating the files named by the
    manifest would let an incremental re-publish retain a stale tool directory
    or failure registry that readers could mistake for current evidence.
    """
    requested_root = Path(root).expanduser()
    if requested_root.is_symlink():
        raise ValueError(f"publication root must not be a symlink: {requested_root}")
    publication = requested_root.resolve()
    try:
        manifest = json.loads(
            (publication / final_artifact(RUN_MANIFEST)).read_text(
                encoding="utf-8",
            ),
        )
        lineage = json.loads(
            (publication / final_artifact(LINEAGE_GRAPH)).read_text(
                encoding="utf-8",
            ),
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("publication manifest or lineage is invalid JSON") from exc
    if not isinstance(manifest, Mapping) or not isinstance(lineage, Mapping):
        raise ValueError("publication manifest and lineage must be JSON objects")
    if manifest.get("layout_contract") != LAYOUT_CONTRACT:
        raise ValueError("publication layout contract is missing or stale")
    if manifest.get("layout_version") != LAYOUT_VERSION:
        raise ValueError("publication layout version is missing or stale")
    publication_mode = manifest.get("publication_mode")
    if publication_mode not in {"pipeline_final", "pipeline_partial", "standalone"}:
        raise ValueError("publication mode is missing or invalid")
    if require_pipeline_final and publication_mode != "pipeline_final":
        raise ValueError("publication is not a finalized pipeline deliverable")
    mandatory_files = (
        PARTIAL_PUBLICATION_FILES if publication_mode == "pipeline_partial"
        else MANDATORY_PUBLICATION_FILES
    )
    missing_mandatory = [
        name for name in mandatory_files if not (publication / name).is_file()
    ]
    if missing_mandatory:
        raise ValueError(
            "publication is missing mandatory files: " + ", ".join(missing_mandatory)
        )
    if (
        expected_session_id is not None
        and manifest.get("session_id") != expected_session_id
    ):
        raise ValueError("publication session_id does not match the round")
    if lineage.get("session_id") != manifest.get("session_id"):
        raise ValueError("publication lineage session_id does not match manifest")
    nodes = lineage.get("nodes")
    edges = lineage.get("edges")
    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise ValueError("publication lineage nodes and edges must be arrays")
    node_ids = [str(node) for node in nodes]
    if len(node_ids) != len(set(node_ids)):
        raise ValueError("publication lineage contains duplicate nodes")
    edge_pairs = []
    for edge in edges:
        if not isinstance(edge, Mapping) or "from" not in edge or "to" not in edge:
            raise ValueError("publication lineage contains an invalid edge")
        edge_pairs.append((str(edge["from"]), str(edge["to"])))
    if len(edge_pairs) != len(set(edge_pairs)):
        raise ValueError("publication lineage contains duplicate edges")

    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, Mapping):
        raise ValueError("publication manifest artifacts must be an object")
    expected_files = {
        final_artifact(RUN_MANIFEST),
        final_artifact(LINEAGE_GRAPH),
    }
    checked = 0
    for raw_relative, raw_expected in artifacts.items():
        relative = Path(str(raw_relative))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"manifest artifact escapes publication: {raw_relative}")
        path = (publication / relative).resolve()
        try:
            path.relative_to(publication)
        except ValueError as exc:
            raise ValueError(
                f"manifest artifact escapes publication: {raw_relative}",
            ) from exc
        if not path.is_file():
            raise ValueError(f"manifest artifact is missing: {raw_relative}")
        if not isinstance(raw_expected, Mapping):
            raise ValueError(f"manifest artifact metadata is invalid: {raw_relative}")
        actual = file_digest(path)
        if raw_expected.get("bytes") != actual["bytes"]:
            raise ValueError(f"manifest size mismatch: {raw_relative}")
        if raw_expected.get("sha256") != actual["sha256"]:
            raise ValueError(f"manifest SHA-256 mismatch: {raw_relative}")
        expected_files.add(relative.as_posix())
        checked += 1

    artifact_ids = {str(relative) for relative in artifacts}
    if not artifact_ids.issubset(set(node_ids)):
        raise ValueError("publication lineage is missing manifest artifact nodes")

    observed_entries = list(publication.rglob("*"))
    unexpected_links = sorted(
        path.relative_to(publication).as_posix()
        for path in observed_entries if path.is_symlink()
    )
    if unexpected_links:
        raise ValueError(
            "publication contains symbolic links: " + ",".join(unexpected_links),
        )
    unexpected_special = sorted(
        path.relative_to(publication).as_posix()
        for path in observed_entries
        if not path.is_file() and not path.is_dir()
    )
    if unexpected_special:
        raise ValueError(
            "publication contains special entries: " + ",".join(unexpected_special),
        )
    observed_files = {
        path.relative_to(publication).as_posix()
        for path in observed_entries if path.is_file()
    }
    untracked = sorted(observed_files - expected_files)
    missing = sorted(expected_files - observed_files)
    if untracked or missing:
        detail = []
        if untracked:
            detail.append("untracked=" + ",".join(untracked))
        if missing:
            detail.append("missing=" + ",".join(missing))
        raise ValueError(
            "publication file set does not match manifest: " + "; ".join(detail),
        )
    source_files_checked = (
        _verify_pipeline_source_snapshot(publication, manifest)
        if publication_mode == "pipeline_final" else 0
    )
    expected_dirs: set[str] = set()
    for relative in expected_files:
        parent = Path(relative).parent
        while parent != Path("."):
            expected_dirs.add(parent.as_posix())
            parent = parent.parent
    observed_dirs = {
        path.relative_to(publication).as_posix()
        for path in observed_entries if path.is_dir()
    }
    observed_top_level = {
        path.name for path in publication.iterdir() if path.is_dir()
    }
    required_top_level = set(PUBLICATION_STAGE_DIRS)
    if observed_top_level != required_top_level:
        detail = []
        extra_top_level = sorted(observed_top_level - required_top_level)
        missing_top_level = sorted(required_top_level - observed_top_level)
        if extra_top_level:
            detail.append("unexpected=" + ",".join(extra_top_level))
        if missing_top_level:
            detail.append("missing=" + ",".join(missing_top_level))
        raise ValueError(
            "publication top-level directories do not match the contract: "
            + "; ".join(detail),
        )
    extra_dirs = sorted(observed_dirs - expected_dirs)
    missing_dirs = sorted(expected_dirs - observed_dirs)
    if extra_dirs or missing_dirs:
        detail = []
        if extra_dirs:
            detail.append("untracked=" + ",".join(extra_dirs))
        if missing_dirs:
            detail.append("missing=" + ",".join(missing_dirs))
        raise ValueError(
            "publication directory set does not match manifest: "
            + "; ".join(detail),
        )
    return {
        "manifest_artifacts_checked": checked,
        "mandatory_files": len(mandatory_files),
        "files": len(observed_files),
        "directories": len(observed_dirs),
        "publication_mode": publication_mode,
        "source_files_checked": source_files_checked,
    }


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

# ---------------------------------------------------------------------------
# The layout object.
# ---------------------------------------------------------------------------


@dataclass
class RunLayout:
    """Path resolver and recorder for one run's published tree.

    Nothing is created until it is asked for, which gives the contract's
    "a tool that did not run this round gets no directory" for free: a caller
    that never publishes a tool never names its directory.

    Every write goes through this object so ``run_manifest.json`` and
    ``lineage_graph.json`` can describe the tree without a second filesystem
    walk that might disagree with it.
    """

    root: Path
    session_id: str = ""
    generated_at: str = field(default_factory=utc_timestamp)
    _written: dict[str, dict[str, Any]] = field(default_factory=dict, repr=False)
    _edges: list[dict[str, Any]] = field(default_factory=list, repr=False)

    def __post_init__(self) -> None:
        self.root = Path(self.root).resolve()
        if not self.session_id:
            # ``root`` *is* the run directory -- the one holding the step
            # directories -- so its own name is the session identity: it already
            # carries the timestamp, the run id and the project slug.  Resolve
            # first, or a relative root such as ``./`` or ``step6/..`` yields an
            # empty name and the column silently ships blank.
            resolved = self.root.resolve()
            self.session_id = resolved.name or resolved.parent.name

    # -- paths --------------------------------------------------------------

    def path(self, *parts: str) -> Path:
        return self.root.joinpath(*parts)

    def step_dir(self, step_id: str) -> Path:
        return self.path(STEP_DIR[step_id])

    def final_dir(self) -> Path:
        return self.path(FINAL_REPORT_DIR)

    def category_dir(self, code: str) -> Path:
        step_id = CATEGORY_CONFIG[code]["step_id"]
        return self.step_dir(step_id) / CATEGORIES_DIR / CATEGORY_DIR[code]

    def tool_run_dir(self, step_id: str, tool_id: str) -> Path:
        return self.step_dir(step_id) / TOOL_RUNS_DIR / tool_dir_name(tool_id)

    def internal_dir(self, step_id: str, tool_id: str) -> Path:
        return self.tool_run_dir(step_id, tool_id) / INTERNAL_DIR

    def shared_assets_dir(self, step_id: str = "step5") -> Path:
        return self.step_dir(step_id) / SHARED_ASSETS_DIR

    def relative(self, path: Path | str) -> str:
        """Publication-relative POSIX path, including ``..`` working refs."""
        candidate = Path(path)
        if not candidate.is_absolute():
            return candidate.as_posix()
        try:
            return candidate.resolve().relative_to(self.root.resolve()).as_posix()
        except ValueError:
            try:
                return Path(
                    os.path.relpath(candidate.resolve(), start=self.root.resolve())
                ).as_posix()
            except ValueError:
                # Different Windows drives cannot be expressed relatively.
                return candidate.resolve().as_posix()

    # -- recording ----------------------------------------------------------

    def record(
        self,
        path: Path,
        *,
        role: str,
        derived_from: Sequence[str] = (),
        digest: bool = True,
    ) -> Path:
        """Note one published artifact and where it came from."""
        key = self.relative(path)
        entry: dict[str, Any] = {"role": role}
        if digest:
            entry.update(file_digest(Path(path)))
        self._written[key] = entry
        for source in derived_from:
            self._edges.append({"from": source, "to": key})
        return path

    def publish_tsv(
        self,
        path: Path,
        columns: Sequence[str],
        rows: Iterable[Mapping[str, Any]],
        *,
        role: str,
        derived_from: Sequence[str] = (),
    ) -> Path:
        write_tsv(Path(path), columns, rows)
        return self.record(path, role=role, derived_from=derived_from)

    def publish_json(
        self,
        path: Path,
        payload: Any,
        *,
        role: str,
        derived_from: Sequence[str] = (),
    ) -> Path:
        write_json(Path(path), payload)
        return self.record(path, role=role, derived_from=derived_from)

    def publish_text(
        self,
        path: Path,
        payload: str,
        *,
        role: str,
        derived_from: Sequence[str] = (),
    ) -> Path:
        write_text(Path(path), payload)
        return self.record(path, role=role, derived_from=derived_from)

    # -- the standardised sidecars -----------------------------------------

    def write_source_map(
        self,
        code: str,
        *,
        metric_columns: Mapping[str, Sequence[str]],
        artifacts: Mapping[str, Sequence[Path | str]] | None = None,
    ) -> Path:
        """``source_map.json``: which tool produced which column, and from where.

        ``metric_columns`` maps a tool id to the category columns it filled;
        ``artifacts`` maps the same tool id to the files those values were read
        out of.  Paths are run-relative references, not copies.
        """
        step_id = CATEGORY_CONFIG[code]["step_id"]
        sources = []
        for tool_id in sorted(metric_columns):
            sources.append({
                "tool_run_id": tool_dir_name(tool_id),
                "tool_id": tool_id,
                "tool_run_path": self.relative(
                    self.tool_run_dir(step_id, tool_id),
                ),
                "metric_columns": list(metric_columns[tool_id]),
                "artifacts": [
                    self.relative(item)
                    for item in (artifacts or {}).get(tool_id, ())
                ],
            })
        return self.publish_json(
            self.category_dir(code) / SOURCE_MAP,
            {
                "category_code": code,
                "category_dir": CATEGORY_DIR[code],
                "step_id": step_id,
                "row_id_field": ROW_ID_FIELD[step_id],
                "generated_at": self.generated_at,
                "sources": sources,
            },
            role="source_map",
        )

    def write_internal_index(self, step_id: str, tool_id: str) -> Path:
        """``_internal/internal_index.tsv``: every file in the tool directory."""
        tool_root = self.tool_run_dir(step_id, tool_id)
        rows = []
        for item in sorted(tool_root.rglob("*")):
            if item.is_file():
                rows.append({
                    "path": self.relative(item),
                    **file_digest(item),
                })
        return self.publish_tsv(
            self.internal_dir(step_id, tool_id) / INTERNAL_INDEX,
            ("path", "bytes", "sha256"),
            rows,
            role="internal_index",
        )

    def register_failure(
        self,
        *,
        stage: str,
        tool_id: str | None,
        reason: str,
        detail: Any = None,
        row_ids: Sequence[str] = (),
    ) -> Path:
        """Append one entry to ``failure_registry.jsonl``.

        A failure belongs here rather than folded into a score: the whole point
        of the registry is that a missing value stays visible as a missing value.
        """
        return append_jsonl(self.final_dir() / FAILURE_REGISTRY, [{
            "recorded_at": utc_timestamp(),
            "stage": stage,
            "tool_id": tool_id,
            "tool_run_id": tool_dir_name(tool_id) if tool_id else None,
            "reason": reason,
            "detail": detail,
            "row_ids": list(row_ids),
        }])

    def finalize(
        self,
        *,
        profile: Mapping[str, Any] | None = None,
        inputs: Mapping[str, Any] | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        """Write ``run_manifest.json`` and ``lineage_graph.json``.

        Call once, after every table is on disk, so the manifest describes the
        tree that actually exists.
        """
        manifest = {
            "session_id": self.session_id,
            "generated_at": self.generated_at,
            "layout_contract": LAYOUT_CONTRACT,
            "layout_version": LAYOUT_VERSION,
            "profile": dict(profile or {}),
            "inputs": {
                key: self.relative(value)
                if isinstance(value, (str, Path)) and value else value
                for key, value in dict(inputs or {}).items()
            },
            "artifacts": {
                key: self._written[key] for key in sorted(self._written)
            },
        }
        manifest.update(extra or {})
        manifest_path = self.final_dir() / RUN_MANIFEST
        write_json(manifest_path, manifest)

        lineage_path = self.final_dir() / LINEAGE_GRAPH
        edges = [
            {"from": source, "to": target}
            for source, target in sorted({
                (edge["from"], edge["to"]) for edge in self._edges
            }, key=lambda pair: (pair[1], pair[0]))
        ]
        write_json(lineage_path, {
            "session_id": self.session_id,
            "generated_at": self.generated_at,
            "nodes": sorted(
                {key for key in self._written}
                | {edge["from"] for edge in edges}
                | {edge["to"] for edge in edges}
            ),
            "edges": edges,
        })
        return manifest_path, lineage_path

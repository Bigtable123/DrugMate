"""SBDD pipeline tool registry — single source of truth for steps/tools/params.

The orchestrator, default-template builder, config validator, and the Agent-facing
tool catalog all read from this module. Nothing here executes science; it only
*declares* what each step can run, what artifacts it consumes/produces, and which
parameters (with defaults + ranges) each tool accepts.

"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

from deployment_config import configure_environments, deployment_path
from runtime_resources import DEFAULT_RESOURCE_SLOTS
from metric_scoring import (
    ALL_METRIC_RULES,
    METRIC_RULES_VERSION,
    NORMALIZATION_METHOD,
    SCREENING_MODES,
    public_metric_rules,
)


# Executable scalar metric policy exposed through the generated catalog.  Only
# metrics currently emitted by the pipeline are included; table rows for
# unavailable tools remain documentation-only until their backends exist.
METRIC_RULES = public_metric_rules(ALL_METRIC_RULES)


# Runtime locations are infrastructure owned by the pipeline, not scientific
# parameters owned by an Agent-authored run config.  A tool references one of
# these environments through its ``runtime`` declaration below.  The shared
# runner resolves ``conda`` dynamically (CONDA_EXE/PATH/common install roots),
# so the registry remains portable across machines with the same env names.
TOOL_ENVIRONMENTS: dict[str, dict] = {
    "current": {"kind": "current"},
    # SBDD_ENV_DRUGMATE_PREFIX / SBDD_ENV_DRUGMATE_NAME select the scientific runtime.
    "drugmate": {"kind": "conda", "name": "hermes"},
    "admet_ai": {"kind": "conda", "name": "admet_ai"},
    "chemfh": {"kind": "conda", "name": "chemfh"},
    "aizynth": {"kind": "conda", "name": "aizynth"},
    "vina": {"kind": "conda", "name": "vina"},
    "unidock": {"kind": "conda", "name": "unidock"},
    "luna": {"kind": "conda", "name": "luna"},
    "oddt": {"kind": "conda", "name": "oddt"},
    "reduce": {
        "kind": "conda",
        "name": "reduce",
        # CCTBX discovers mon_lib/geostd and CCD independently from its official
        # chem_data package. Never redirect the monomer library to CCD: geostd
        # files placed there pass a shallow CIF check but break Reduce2.
    },
    "amber": {"kind": "conda", "name": "amber"},
    # This adapter uses AmberTools MMPBSA.py, not the gmx_MMPBSA CLI.
    "mmgbsa": {"kind": "conda", "name": "gmxMMPBSA"},
    "openfe": {"kind": "conda", "name": "openfe"},
    "ties": {"kind": "conda", "name": "ties"},
    "pocketminer": {"kind": "conda", "name": "pocketminer"},
    "pocket2mol": {"kind": "conda", "name": "Pocket2Mol"},
    "difflinker": {"kind": "conda", "name": "difflinker"},
    "diffsbdd": {"kind": "conda", "name": "diffsbdd"},
    "reinvent4": {"kind": "conda", "name": "reinvent4"},
    "resgen": {"kind": "conda", "name": "resgen"},
    # TargetDiff requires a compatible CUDA environment.
    "targetdiff": {"kind": "conda", "name": "targetdiff-gpu"},
    "deeppocket": {"kind": "conda", "name": "deeppocket"},
    "diffdock": {"kind": "conda", "name": "diffdock"},
    "dock6": {"kind": "conda", "name": "dock6"},
    "autodock4": {"kind": "conda", "name": "autodock4"},
    "rdock": {"kind": "conda", "name": "rdock"},
    "adfrsuite": {"kind": "conda", "name": "adfrsuite"},
    "lingo3dmol": {"kind": "conda", "name": "lingo3dmol"},
    "rascore": {"kind": "conda", "name": "rascore"},
    "rtmscore": {"kind": "conda", "name": "rtmscore_clean"},
    "tamgen": {"kind": "conda", "name": "tamgen"},
    "molcraft": {"kind": "conda", "name": "molcraft"},
    # Short MD requires an OpenMM environment with CUDA support.
    "openmm": {"kind": "conda", "name": "openmm"},
    "proteins_plus": {
        "kind": "http",
        "base_url": "https://proteins.plus/api",
    },
}


DEFAULT_TOOL_ENVIRONMENTS = deepcopy(TOOL_ENVIRONMENTS)
TOOL_ENVIRONMENTS = configure_environments(DEFAULT_TOOL_ENVIRONMENTS)

# Canonical primary engines for Step 4.  AutoGrid4 and AutoDockTools are
# components of the AutoDock4 route, while AGFR is a component of ADFR; none
# is a second scoring engine or an independent enable switch.
DOCKING_PRIMARY_ENGINES = (
    "gnina",
    "unidock",
    "autodock_vina",
    "diffdock",
    "dock6",
    "autodock4",
    "rdock",
    "adfr",
)


_BIO_RESEARCH_DIR = Path(__file__).resolve().parent.parent
_POCKETMINER_BRIDGE = (
    _BIO_RESEARCH_DIR / "pocket-analysis" / "scripts" / "pocketminer_bridge.py"
)
_POCKETMINER_SOURCE_DIR = deployment_path("${SBDD_POCKETMINER_SOURCE_DIR}")
_POCKETMINER_CHECKPOINT = _POCKETMINER_SOURCE_DIR / "models" / "pocketminer"
_POCKET2MOL_SOURCE_DIR = deployment_path("${SBDD_POCKET2MOL_SOURCE_DIR}")
_POCKET2MOL_CHECKPOINT = _POCKET2MOL_SOURCE_DIR / "pretrained_Pocket2Mol.pt"
_POCKET2MOL_CHECKPOINT_SHA256 = (
    "e9d8152dfdf034dabf48e566904e9d7c6cc1b96843bdc4d0026d891678913911"
)
_DIFFLINKER_SOURCE_DIR = deployment_path("${SBDD_DIFFLINKER_SOURCE_DIR}")
_DIFFLINKER_CHECKPOINT = _DIFFLINKER_SOURCE_DIR / "models" / "pockets_difflinker_full.ckpt"
_DIFFSBDD_SOURCE_DIR = deployment_path("${SBDD_DIFFSBDD_SOURCE_DIR}")
_DIFFSBDD_CHECKPOINT = _DIFFSBDD_SOURCE_DIR / "checkpoints" / "crossdocked_fullatom_cond.ckpt"
_REINVENT4_SOURCE_DIR = deployment_path("${SBDD_REINVENT4_SOURCE_DIR}")
_REINVENT4_PRIOR = _REINVENT4_SOURCE_DIR / "priors" / "reinvent_pubchem.prior"
_RESGEN_SOURCE_DIR = deployment_path("${SBDD_RESGEN_SOURCE_DIR}")
_RESGEN_CHECKPOINT = _RESGEN_SOURCE_DIR / "ckpt" / "val_172.pt"
_TARGETDIFF_SOURCE_DIR = deployment_path("${SBDD_TARGETDIFF_SOURCE_DIR}")
_TARGETDIFF_CHECKPOINT = _TARGETDIFF_SOURCE_DIR / "pretrained_models" / "pretrained_diffusion.pt"
_DEEPPOCKET_SOURCE_DIR = deployment_path("${SBDD_DEEPPOCKET_SOURCE_DIR}")
_DEEPPOCKET_RANK_CHECKPOINT = (
    _DEEPPOCKET_SOURCE_DIR / "models" / "first_model_fold1_best_test_auc_85001.pth.tar"
)
_DEEPPOCKET_SEGMENT_CHECKPOINT = (
    _DEEPPOCKET_SOURCE_DIR / "models" / "seg0_best_test_IOU_91.pth.tar"
)
_DIFFDOCK_SOURCE_DIR = deployment_path("${SBDD_DIFFDOCK_SOURCE_DIR}")
_DIFFDOCK_SCORE_CHECKPOINT = (
    _DIFFDOCK_SOURCE_DIR / "workdir" / "v1.1" / "score_model"
    / "best_ema_inference_epoch_model.pt"
)
_DIFFDOCK_CONFIDENCE_CHECKPOINT = (
    _DIFFDOCK_SOURCE_DIR / "workdir" / "v1.1" / "confidence_model"
    / "best_model_epoch75.pt"
)
_DOCK6_SOURCE_DIR = deployment_path("${SBDD_DOCK6_SOURCE_DIR}")
_DOCK6_BINARY = _DOCK6_SOURCE_DIR / "bin" / "dock6"
_LINGO3DMOL_SOURCE_DIR = deployment_path("${SBDD_LINGO3DMOL_SOURCE_DIR}")
_LINGO3DMOL_CONTACT_CHECKPOINT = _LINGO3DMOL_SOURCE_DIR / "checkpoint" / "contact.pkl"
_LINGO3DMOL_GENERATOR_CHECKPOINT = _LINGO3DMOL_SOURCE_DIR / "checkpoint" / "gen_mol.pkl"
_TAMGEN_SOURCE_DIR = deployment_path("${SBDD_TAMGEN_SOURCE_DIR}")
_TAMGEN_CHECKPOINT = (
    _TAMGEN_SOURCE_DIR / "checkpoints" / "crossdock_pdb_A10" / "checkpoint_best.pt"
)
_MOLCRAFT_SOURCE_DIR = deployment_path("${SBDD_MOLCRAFT_SOURCE_DIR}")
_MOLCRAFT_CHECKPOINT = _MOLCRAFT_SOURCE_DIR / "checkpoints" / "last.ckpt"
_MOLCRAFT_CHECKPOINT_CONFIG = (
    _MOLCRAFT_SOURCE_DIR / "checkpoints" / "last.ckpt.config.yaml"
)
_RTMSCORE_SOURCE_DIR = deployment_path("${SBDD_RTMSCORE_SOURCE_DIR}")
_RTMSCORE_ENTRYPOINT = _RTMSCORE_SOURCE_DIR / "example" / "rtmscore.py"
_RTMSCORE_CHECKPOINTS = tuple(
    _RTMSCORE_SOURCE_DIR / "trained_models" / f"rtmscore_model{index}.pth"
    for index in (1, 2, 3)
)
_AIZYNTHFINDER_CONFIG = deployment_path("${SBDD_AIZYNTHFINDER_CONFIG}")
_AIZYNTHFINDER_PROBE = (
    _BIO_RESEARCH_DIR
    / "drug-property-filter"
    / "scripts"
    / "aizynth_probe.py"
)
_CHEMFH_SOURCE_DIR = deployment_path("${SBDD_CHEMFH_SOURCE_DIR}")
_CHEMFH_CHECKPOINT = (
    _CHEMFH_SOURCE_DIR
    / "static"
    / "media"
    / "chemprop"
    / "checkpoint"
    / "fold_0"
    / "model_0"
    / "model.pt"
)
_CHEMFH_CHECKPOINT_SHA256 = (
    "9f5c82630bed2d787eaf16d99eeac500763748d54474818757207659433a8b6b"
)
_CHEMFH_RULE_BUNDLE_SHA256 = (
    "f67d0cf228af15cd3ccf200034f0363a1495aaba0f25388e0393379cd1e55519"
)
_CHEMFH_SOURCE_COMMIT = "f8c046fac9a1d331979025a104382d3d2b2614ca"
_CHEMFH_WORKER = (
    _BIO_RESEARCH_DIR
    / "drug-property-filter"
    / "scripts"
    / "chemfh_worker.py"
)


# --- Artifact vocabulary: the only medium for inter-step wiring ---------------
ARTIFACTS = {
    "provisional_receptor", "prep_manifest",
    "source_complex_pdb", "water_inventory", "receptor_state_manifest",
    "docking_species_library", "microstate_identity_map", "docking_prep_qc",
    "prepared_receptor",
    "inferred_target_name",
    "prepared_receptor_noH",
    "source_complex_mmcif",
    "receptor_pdbqt",
    "reference_ligand",
    "reference_ligand_microstates",
    "protonation_manifest",
    "water_report",
    "metal_report",
    "structure_qc_report",
    "amber_prmtop",
    "amber_inpcrd",
    "amber_ligand_mol2",
    "amber_ligand_frcmod",
    "pocket_center",
    "pocket_definition",
    "docking_box",
    "dock6_receptor_site",
    "dock6_grid_nrg",
    "dock6_grid_bmp",
    "pocket_features",
    "interaction_observations",
    "interaction_blueprint",
    "pocket_analysis_report",
    "deeppocket_predictions",
    "library",
    "library_manifest",
    "prefiltered_library",
    "docked_poses",
    "docking_scores",
    "docking_report",
    "selected_poses",
    "selection_manifest",
    "evaluation_report",
    "candidate_topologies",
    "mmgbsa_manifest",
    "openfe_rbfe_manifest",
    "ties_rbfe_manifest",
    "rtmscore_manifest",
    "short_md_manifest",
    "trajectory_stability_report",
    "property_report",
    "synthesis_report",
    "coordinator_report",
    "report_html",
    "report_data",
    "report_assets_manifest",
    "report_artifact_manifest",
    "report_pymol_cover",
}


def _p(default: Any, ptype: str, *, min=None, max=None, choices=None,
       item: str | None = None, item_choices=None,
       required: bool = False, desc: str = "") -> dict:
    """Declare one parameter.

    ``None`` means "no default", not "required". Requiredness is explicit so
    optional selectors such as ``ligand_comp_id`` are not mislabeled in the
    generated catalog.
    """
    return {
        "default": default, "type": ptype, "min": min, "max": max,
        "choices": choices, "item": item, "item_choices": item_choices,
        "required": required, "desc": desc,
    }


# --- context: shared identity / knowledge inputs (validated like params) ------
SCIENTIFIC_FACTS_SCHEMA: dict[str, dict] = {
    "has_metal": _p(
        None,
        "bool",
        desc="Metal or metallocofactor participates in the binding site.",
    ),
    "is_covalent": _p(
        None,
        "bool",
        desc="The requested chemistry uses a covalent mechanism.",
    ),
    "protonation_sensitive": _p(
        None,
        "bool",
        desc="Binding depends on uncertain ionization or tautomer states.",
    ),
}


# Every key a user may put under `context`. Unknown keys hard-fail (VR2-style),
# matching the "all config keys must be known" contract.
CONTEXT_SCHEMA: dict[str, dict] = {
    "pdb_id": _p(None, "str", desc="4-char RCSB PDB id"),
    "protein_name": _p(None, "str", desc="Protein name for RCSB search"),
    "receptor_file": _p(None, "str", desc="Existing receptor PDB (prep-only mode)"),
    "ligand_file": _p(None, "str", desc="Existing ligand PDB (prep-only mode)"),
    "fragment_file": _p(
        None, "str",
        desc=(
            "3D molecule file containing at least two disconnected fragments; "
            "required only for the DiffLinker route."
        ),
    ),
    "complex_mmcif_file": _p(
        None, "str",
        desc=(
            "Optional untouched source complex mmCIF for prep-only mode and "
            "Arpeggio supplemental interaction analysis."
        ),
    ),
    "target_name": _p(None, "str"),
    "chembl_target_id": _p(None, "str"),
    "uniprot_id": _p(None, "str"),
    "chain": _p(None, "str", desc="Explicit target author chain; resolve entity/assembly first, never assume A."),
    "seed": _p(42, "int"),
    "device": _p(None, "int", min=0, max=16, desc="Default CUDA device"),
    "docking_memory_mode": _p(
        "ligand_protein", "str", required=True,
        choices=["off", "legacy", "ligand_only", "ligand_protein"],
        desc=(
            "Agent-side docking planning mode; default ligand_protein. "
            "off uses no docking memory, legacy recalls protein neighbors, "
            "ligand_only disables protein reweighting. The runner never "
            "queries memory or changes explicit tool selections."
        ),
    ),
    "candidate_source": _p(
        None,
        "str",
        choices=["molecule_generation", "sampled_library", "user_supplied_full"],
        desc=(
            "Explicit molecule provenance for automatic V1/V2/V3 selection. "
            "Use sampled_library when an externally supplied large library was "
            "sampled before injection; omit only when artifact provenance is "
            "sufficient to derive the source."
        ),
    ),
    "known_site_residues": _p(None, "list", item="str"),
    "known_site_provenance": _p(
        None, "dict",
        desc=(
            "Optional provenance for known_site_residues: source, method, "
            "independent_of_reference_structure and review_status."
        ),
    ),
    "known_required_interactions": _p(
        None, "list", item="dict",
        desc=(
            "Approved curated hard constraints. Each item requires residue, "
            "interaction_type, rationale, source, binding_mode_scope and "
            "review_status='approved'; atom/features and logic groups are optional."
        ),
    ),
    "scientific_facts": _p(
        None,
        "dict",
        desc=(
            "Declared chemistry facts passed to capability checks and tools. "
            "They never enable, disable, or select tools."
        ),
    ),
}

# Explicit context file inputs that also satisfy canonical artifact contracts.
# This is a wiring declaration only: it must never be inferred from holo/apo
# state or other scientific facts.
CONTEXT_ARTIFACTS = {
    "source_complex_mmcif": "complex_mmcif_file",
}

# --- project: run/directory control -------------------------------------------
PROJECT_SCHEMA: dict[str, dict] = {
    "project_slug": _p(None, "str"),
    "data_root": _p(None, "str"),
    "project_dir": _p(None, "str", desc="Explicit override; skips auto-naming"),
    "round_number": _p(1, "int", min=1, max=100),
    "campaign_state_file": _p(None, "str", desc="Durable R1-R3 campaign authorization and execution budget."),
    "campaign_stage": _p(None, "str", choices=["R1", "R2", "R3"], desc="Scientific stage; independent of round_number and retry attempts."),
    "campaign_attempt": _p(None, "int", min=1, desc="Attempt within the scientific stage; checkpoint resumes reuse the same attempt."),
    "round_decision_file": _p(None, "str", desc="Agent-authored plan bound to this config; required for campaign execution."),
    "upstream_project_dir": _p(None, "str", desc="Prior round dir to reuse artifacts from"),
    "execution_mode": _p(
        "interrupt",
        "str",
        choices=["interrupt", "auto"],
        desc=(
            "Agent orchestration mode. 'interrupt' checkpoints after every "
            "newly processed step so the Agent can request explicit user "
            "confirmation; 'auto' runs all remaining steps. The initial "
            "steps/tools plan still requires confirmation in both modes."
        ),
    ),
    "resume_existing": _p(
        True,
        "bool",
        desc=(
            "Reuse a successful step result only when its resolved input/tool/runtime and "
            "local-file fingerprint matches; PocketMiner-enabled Step 1 always reruns. "
            "This controls ordinary cache reuse for a pending step, not the separately "
            "validated processed-step prefix in an interactive checkpoint."
        ),
    ),
    "python_executable": _p(None, "str"),
}

RUNTIME_RESOURCE_SCHEMA = {
    "wait_timeout_seconds": _p(0, "int", min=0, max=2592000),
    "slots": _p(DEFAULT_RESOURCE_SLOTS, "int", min=1, max=64),
    "lock_directory": _p(None, "str", desc="Shared absolute lock directory; defaults to data_root/_resources"),
}

# --- top-level convenience keys (not project/context/steps) -------------------
# start_step/stop_step bound the enabled range; input_artifacts injects
# externally-supplied artifacts (e.g. a modified SDF) into the wiring store.
TOP_LEVEL_KEYS = {
    "project", "context", "steps", "start_step", "stop_step",
    "input_artifacts", "step_timeouts", "runtime_resources",
}


# Artifacts a user may inject via top-level `input_artifacts` to feed a step
# without running its producer — the key mechanism for iterating with a new/
# modified library. Restricted to file-path artifacts that are safe to override.
INJECTABLE_ARTIFACTS = {
    "provisional_receptor", "prepared_receptor_noH", "prep_manifest",
    "source_complex_pdb", "water_inventory", "receptor_state_manifest",
    "docking_species_library", "microstate_identity_map", "docking_prep_qc",
    "library", "library_manifest", "prefiltered_library",
    "reference_ligand", "prepared_receptor",
    "receptor_pdbqt", "source_complex_mmcif", "docking_box", "docked_poses",
    "selected_poses", "dock6_receptor_site", "dock6_grid_nrg",
    "dock6_grid_bmp",
}


# STEPS is assembled below. Each step:
#   order, skill, default_enabled, consumes[], produces[], params{}, tools{}
# Each tool: default_enabled, purpose, when_to_use, depends_on[], requires[],
#            params{}, optional runtime{}
STEPS: dict[str, dict] = {}
# __REGISTRY_APPEND__

STEPS["structure_prep"] = {
    "order": 0,
    "skill": "structure-prep",
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    "consumes": [],
    "produces": [
        "inferred_target_name", "prepared_receptor_noH", "provisional_receptor",
        "source_complex_mmcif", "source_complex_pdb", "water_inventory",
        "reference_ligand", "protonation_manifest", "metal_report",
        "structure_qc_report", "prep_manifest",
    ],
    "params": {
        "ph": _p(None, "float", min=0.0, max=14.0,
                 desc="Optional diagnostic pH, with experimental justification; final state belongs to docking_prep."),
    },
    "tools": {
        "intake": {
            "default_enabled": True,
            "purpose": "Download PDB from RCSB and split receptor/ligand.",
            "when_to_use": "Starting from a pdb_id or protein_name.",
            "depends_on": [], "requires": [],
            "produces": ["inferred_target_name", "source_complex_mmcif"],
            "runtime": {"kind": "internal"},
            "params": {
                "water_distance_threshold": _p(3.5, "float", min=0.0, max=10.0),
                "keep_auxiliary_hetero": _p(False, "bool"),
                "preserve_hetero_resnames": _p(None, "list", item="str"),
                "skip_prep": _p(False, "bool", desc="Intake only, no preparation."),
            },
        },
        "pdbfixer": {
            "default_enabled": True,
            "purpose": "Repair supported standard heavy atoms; canonical receptor remains hydrogen-free.",
            "when_to_use": "Use for defensible heavy-atom repairs; short internal gaps require explicit justification.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module", "environment": "drugmate",
                "module": "pdbfixer",
            },
            "params": {
                # Optional legacy override. Prefer structure_prep.params.ph.
                "ph": _p(None, "float", min=0.0, max=14.0),
                "add_missing_residues": _p(False, "bool"),
                "max_missing_residues": _p(3, "int", min=1, max=3),
                "repair_reason": _p(None, "str"),
                "seed": _p(0, "int", min=0),
                "add_missing_atoms": _p(True, "bool"),
                "add_hydrogens": _p(False, "bool"),
            },
        },
        "propka": {
            "default_enabled": False,
            "purpose": "Report pKa / protonation states (diagnostic).",
            "when_to_use": "When protonation review matters; not a state rewrite.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module", "environment": "drugmate",
                "module": "propka",
            },
            "params": {
                # Optional legacy override. Prefer structure_prep.params.ph.
                "ph": _p(None, "float", min=0.0, max=14.0),
            },
        },
        "reduce": {
            "default_enabled": False,
            "purpose": "Provisional hydrogens/flips for an explicitly hydrogen-dependent pocket analysis.",
            "when_to_use": "Local flip backend when the Reduce2 monomer library is configured.",
            "depends_on": [], "requires": [],
            "produces": ["protonation_manifest", "provisional_receptor"],
            "runtime": {
                "kind": "command", "environment": "reduce",
                "command": "mmtbx.reduce2", "probe_args": ["--version"],
                "probe_code": (
                    "import shutil,sys;"
                    "from mmtbx.monomer_library import server;"
                    "from mmtbx import chemical_components;"
                    "from mmtbx.ligands.rdkit_utils import get_molecule_from_resname;"
                    "server.server();"
                    "assert chemical_components.get_cif_dictionary('ALA') is not None;"
                    "mol=get_molecule_from_resname('ALA');"
                    "assert mol is not None and mol.GetNumAtoms()>0, 'ALA CCD molecule could not be constructed';"
                    "sys.exit(0 if shutil.which('mmtbx.reduce2') else 1)"
                ),
                "probe_failure_hint": (
                    "mmtbx.reduce2, its CCP4/geostd monomer library, or the CCTBX "
                    "CCD chemical-components database is unavailable or incompatible. "
                    "Configure MMTBX_CCP4_MONOMER_LIB or CLIBD_MON for CCP4/geostd "
                    "restraints and a separate genuine CCD at the libtbx "
                    "chemical_components path; geostd data_*.cif files are not CCD files."
                ),
            },
            "params": {"build": _p(True, "bool"), "flip": _p(True, "bool"),
                       "reason": _p(None, "str", required=True, desc="Pocket detector/profile requiring provisional hydrogens.")},
        },
        "ligand_sdf": {
            "default_enabled": True,
            "purpose": "Reference ligand bound-pose SDF with CCD template bond orders.",
            "when_to_use": "Needed for autobox and PoseCheck calibration downstream.",
            "depends_on": [], "requires": [],
            "produces": ["reference_ligand"],
            "runtime": {
                "kind": "python_module", "environment": "drugmate",
                "module": "rdkit",
            },
            "params": {"ligand_comp_id": _p(None, "str")},
        },
        "metal_validation": {
            "default_enabled": True,
            "purpose": "Report retained metals, nearby donors, distances, and coordination counts.",
            "when_to_use": "Always safe; mandatory for declared metal sites.",
            "depends_on": [], "requires": [],
            "produces": ["metal_report"],
            "runtime": {"kind": "internal"},
            "params": {
                "coordination_cutoff": _p(3.0, "float", min=1.5, max=4.0),
                "allowed_metals": _p(
                    ["MG", "ZN", "CA", "FE", "MN", "NI", "CU", "CO", "NA", "K", "CD"],
                    "list", item="str",
                ),
            },
        },
        "structure_qc": {
            "default_enabled": True,
            "purpose": "Final audit of atoms, hydrogens, alternates, waters, metals, and generated files.",
            "when_to_use": "Always; provides the Step 0 acceptance report.",
            "depends_on": [], "requires": [],
            "produces": ["structure_qc_report"],
            "runtime": {"kind": "internal"},
            "params": {
                "fail_on_altloc": _p(False, "bool"),
                # Defaults to the run's own ``add_hydrogens`` when unset; set it
                # explicitly to demand hydrogens regardless of that choice.
                "require_hydrogens": _p(None, "bool"),
            },
        },
        "residue_validation": {
            "default_enabled": True,
            "purpose": "Verify active-site residues are present and unmutated.",
            "when_to_use": "When known_site residues must be confirmed.",
            "depends_on": [], "requires": [],
            "runtime": {"kind": "internal"},
            "params": {"required_residues": _p(None, "list", item="dict")},
        },
    },
}

STEPS["pocket_analysis"] = {
    "order": 1,
    "skill": "pocket-analysis",
    # Keep pocket detection and interaction analysis in one Step 1 while
    # guaranteeing that its native profilers execute in the validated Drugmate
    # environment (their isolated workers inherit this interpreter).
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    "consumes": [
        "provisional_receptor", "prepared_receptor_noH", "reference_ligand",
        "source_complex_mmcif",
        "protonation_manifest", "water_inventory",
        "metal_report", "structure_qc_report",
    ],
    "requires_all": ["prepared_receptor_noH"],
    "optional_consumes": ["docking_box"],
    "produces": ["pocket_center", "pocket_definition", "docking_box",
                 "pocket_features", "interaction_observations",
                 "interaction_blueprint", "pocket_analysis_report",
                 "deeppocket_predictions"],
    "params": {
        "site_selection_strategy": _p(
            None,
            "str",
            choices=["ligand_site", "known_site", "detector_consensus", "provided_box"],
            required=True,
            desc=(
                "Explicit pocket-selection algorithm chosen by the Agent. "
                "Execution never infers this value from ligand availability. "
                "provided_box characterizes the supplied docking_box without replacing it. "
                "detector_consensus uses plan-pocket-tools to select native Top-3 detectors, using protein "
                "memory when available and fixed DeepPocket#1 rank fusion only when "
                "memory is unavailable or explicitly disabled."
            ),
        ),
        "pocket_memory_mode": _p(
            "auto", "str", choices=["auto", "off", "required"],
            desc=(
                "auto uses compatible protein memory or records an unavailable-memory "
                "fallback; off explicitly skips memory; required rejects unavailable "
                "memory. Invalid current candidates or provenance always fail."
            ),
        ),
        "pocket_memory_index": _p(
            str(deployment_path("${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json")),
            "str", desc="Derived protein-only pocket reliability index; build from the memory split.",
        ),
        "pocket_memory_query_npz": _p(
            None, "str", desc="Current receptor SaProt embeddings.npz for pocket_prediction.",
        ),
        "pocket_memory_query_metadata": _p(
            None, "str", desc="Hash-bound pocket_prediction query metadata for the current receptor.",
        ),
        "pocket_memory_preset": _p(
            "conservative", "str", choices=["conservative", "exploratory_memory_guard"],
            desc="Pocket anchor-switch and reliability policy preset.",
        ),
        "pocket_memory_query_case": _p(None, "str", desc="Optional query case identity; must agree with metadata."),
        "pocket_memory_query_group": _p(None, "str", desc="Optional query protein-group identity; must agree with metadata."),
        "pocket_memory_exclude_cases": _p([], "list", item="str", desc="Additional historical cases excluded before retrieval."),
        "pocket_memory_exclude_groups": _p([], "list", item="str", desc="Additional protein groups excluded before retrieval."),
        "native_step_timeout_seconds": _p(300, "int", min=30, max=7200),
        "fail_on_profiler_error": _p(
            False, "bool",
            desc=(
                "Fail Step 1 when an enabled and applicable interaction profiler "
                "fails.  Default false: a profiler outage is a fault against the "
                "protocol, not evidence about the pocket, so it is recorded in "
                "the failure registry while sound geometry still publishes. "
                "Set true only for a benchmark that must abort "
                "unless every profiler ran."
            ),
        ),
        "fail_on_supplemental_profiler_error": _p(
            False, "bool",
            desc=(
                "Fail Step 1 when an enabled supplemental profiler fails. "
                "Default false keeps Arpeggio/LUNA/ODDT diagnostic and portable."
            ),
        ),
        "consensus_top_k": _p(3, "int", min=1, max=20),
        "center_match_threshold": _p(5.0, "float", min=0.0, max=30.0),
        "residue_iou_threshold": _p(0.3, "float", min=0.0, max=1.0),
        "detector_box_padding": _p(
            4.0,
            "float",
            min=0.0,
            max=30.0,
            desc=(
                "Padding added around detector-consensus pocket residues when "
                "publishing the docking box. Reduce only under an explicit "
                "engine grid-size constraint."
            ),
        ),
    },
    # provided_box uses the supplied geometry; other modes require a site tool.
    "one_of": [["ligand_site", "fpocket", "pykvfinder", "p2rank", "deeppocket", "autosite"]],
    "tools": {
        "ligand_site": {
            "default_enabled": True,
            "purpose": "Define a pocket from the observed reference-ligand shell.",
            "when_to_use": "Primary source whenever a usable reference ligand is present.",
            "depends_on": [], "requires": ["reference_ligand"],
            "runtime": {"kind": "internal"},
            "params": {
                "residue_cutoff": _p(5.0, "float", min=0.1, max=15.0),
                "box_padding": _p(4.0, "float", min=0.0, max=30.0),
                "min_box_size": _p(12.0, "float", min=1.0, max=100.0),
            },
        },
        "fpocket": {
            "default_enabled": True,
            "purpose": "Geometric (Voronoi) pocket detection.",
            "when_to_use": "When recommended by plan-pocket-tools, or in the four-tool no-memory fallback; shares the fpocket family with DeepPocket.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command", "environment": "drugmate", "command": "fpocket",
                "probe_args": ["--help"],
                "probe_failure_hint": "fpocket executable probe failed in the drugmate environment",
            },
            "params": {
                "min_radius": _p(3.0, "float", min=0.0, max=10.0),
                "max_radius": _p(6.0, "float", min=0.0, max=15.0),
            },
        },
        "pykvfinder": {
            "default_enabled": False,
            "purpose": "Pure-Python cavity detection.",
            "when_to_use": "When recommended by plan-pocket-tools, or in the four-tool no-memory fallback; contributes anchor-matched native Top-3 candidates.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module", "environment": "drugmate", "module": "pyKVFinder",
            },
            "params": {
                "probe_in": _p(1.4, "float", min=0.0, max=10.0),
                "probe_out": _p(4.0, "float", min=0.0, max=15.0),
                "removal_distance": _p(2.4, "float", min=0.0, max=15.0),
                "volume_cutoff": _p(5.0, "float", min=0.0, max=100.0),
            },
        },
        "p2rank": {
            "default_enabled": True,
            "purpose": "ML ligand-binding-site scoring and ranking.",
            "when_to_use": "When recommended by plan-pocket-tools, or in the four-tool no-memory fallback; may replace the DeepPocket anchor only when every memory guard passes.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command", "environment": "drugmate", "command": "prank",
                # P2Rank 2.6 prints help/version but exits non-zero; probe the
                # executable path instead so preflight does not reject a healthy install.
                "probe_code": (
                    "import shutil,sys; "
                    "sys.exit(0 if shutil.which('prank') else 1)"
                ),
                "probe_failure_hint": "P2Rank prank executable probe failed in the drugmate environment",
            },
            "params": {},
        },
        "deeppocket": {
            "default_enabled": False,
            "purpose": (
                "DeepPocket segmentation plus learned pocket ranking, normalized "
                "as geometry-bearing Step 1 detector evidence."
            ),
            "when_to_use": (
                "Required default native Top-1 anchor for detector_consensus; "
                "memory may authorize a P2Rank anchor, while the no-memory fallback "
                "fixes DeepPocket#1. Shares one detector family with fpocket."
            ),
            "depends_on": [],
            "requires": [],
            "produces": ["deeppocket_predictions"],
            "runtime": {
                "kind": "command",
                "environment": "deeppocket",
                "command": "python",
                "bridge_path": str(
                    _BIO_RESEARCH_DIR / "pocket-analysis" / "scripts"
                    / "deeppocket_bridge.py"
                ),
                "source_dir": str(_DEEPPOCKET_SOURCE_DIR),
                "class_checkpoint": str(_DEEPPOCKET_RANK_CHECKPOINT),
                "segmentation_checkpoint": str(_DEEPPOCKET_SEGMENT_CHECKPOINT),
                "asset_sha256": {
                    "class_checkpoint": "d166c8a36b2475297193dbe31e5a88ff9b853e343281ef21df537029ba053f9a",
                    "segmentation_checkpoint": "aa6d29beda70e79b39504c3bdaf5b42380e9c1e5bc6b0557cb81e9075243faee",
                },
                "probe_args": [
                    str(
                        _BIO_RESEARCH_DIR / "pocket-analysis" / "scripts"
                        / "deeppocket_bridge.py"
                    ),
                    "probe",
                    "--source-dir", str(_DEEPPOCKET_SOURCE_DIR),
                    "--class-checkpoint", str(_DEEPPOCKET_RANK_CHECKPOINT),
                    "--segmentation-checkpoint", str(_DEEPPOCKET_SEGMENT_CHECKPOINT),
                ],
                "probe_timeout_seconds": 180,
                "probe_failure_hint": (
                    "DeepPocket imports, CUDA allocation, fpocket executable, "
                    "or strict checkpoint restoration failed."
                ),
            },
            "params": {
                "max_candidates": _p(3, "int", min=1, max=10),
                # The upstream sigmoid mask contract is open-interval (0, 1).
                "segmentation_threshold": _p(0.5, "float", min=0.001, max=0.999),
                "mask_distance": _p(3.5, "float", min=0.1, max=20.0),
                "timeout_seconds": _p(1200, "int", min=60, max=86400),
            },
        },
        "autosite": {
            "default_enabled": False,
            "purpose": "AutoSite cluster detection with native ranks and exported point centroids.",
            "when_to_use": "Only for detector_consensus when included in plan-pocket-tools.execution_tools.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command", "environment": "adfrsuite", "command": "autosite",
                "commands": {"prepare_receptor": "prepare_receptor", "engine": "autosite"},
                "probe_code": (
                    "import shutil; assert shutil.which('autosite'); "
                    "assert shutil.which('prepare_receptor')"
                ),
                "probe_failure_hint": "AutoSite and prepare_receptor must be available in the configured ADFRsuite environment",
            },
            "params": {"timeout_seconds": _p(1200, "int", min=60, max=86400)},
        },
        "pocketminer": {
            "default_enabled": False,
            "purpose": (
                "Residue-level cryptic-opening propensity annotation; never reranks "
                "geometry-bearing candidates or adds detector support."
            ),
            "when_to_use": (
                "Explicit cryptic-pocket analysis on one prepared protein chain. "
                "Requires at least one separate cavity/site detector for a dockable definition."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "pocketminer",
                "command": "python",
                "bridge_path": str(_POCKETMINER_BRIDGE),
                "source_dir": str(_POCKETMINER_SOURCE_DIR),
                "checkpoint": str(_POCKETMINER_CHECKPOINT),
                "probe_args": [
                    str(_POCKETMINER_BRIDGE), "probe",
                    "--source-dir", str(_POCKETMINER_SOURCE_DIR),
                    "--checkpoint", str(_POCKETMINER_CHECKPOINT),
                ],
                "probe_timeout_seconds": 120,
                "probe_failure_hint": (
                    "PocketMiner strict checkpoint/model-graph validation failed; "
                    "do not use its residue scores"
                ),
            },
            "params": {
                "timeout_seconds": _p(600, "int", min=30, max=7200),
                "rerank_mode": _p(
                    "annotate", "str", choices=["annotate", "tie_break"],
                    desc=(
                        "Annotation only; legacy tie_break is accepted but cannot change "
                        "the selected pocket, fusion center, or anchor residues."
                    ),
                ),
                "min_score_coverage": _p(0.8, "float", min=0.0, max=1.0),
                "provisional_quantile": _p(
                    0.9, "float", min=0.5, max=0.999,
                    desc="Pipeline heuristic for provisional residue clusters; not an upstream cutoff.",
                ),
                "cluster_ca_cutoff": _p(8.0, "float", min=1.0, max=30.0),
                "min_cluster_residues": _p(3, "int", min=1, max=100),
                "max_hypotheses": _p(5, "int", min=1, max=50),
            },
        },
        "plip": {
            "default_enabled": True,
            "purpose": "Primary protein-ligand interaction profiler (isolated process).",
            "when_to_use": "Whenever a reference ligand is present.",
            "depends_on": [], "requires": ["reference_ligand"],
            "runtime": {
                "kind": "command",
                "environment": "drugmate",
                "command": "python",
                "probe_code": (
                    "import plip; "
                    "from plip.structure.preparation import PDBComplex"
                ),
                "probe_failure_hint": "PLIP import/API probe failed in the drugmate environment",
            },
            "params": {},
        },
        "prolif": {
            "default_enabled": True,
            "purpose": "Record-only interaction fingerprint audit (secondary, isolated process).",
            "when_to_use": "For a record-only fingerprint audit when a reference ligand and explicit-H provisional receptor are available; otherwise record not_evaluable and continue.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "drugmate",
                "command": "python",
                "probe_code": "import prolif, MDAnalysis, rdkit",
                "probe_failure_hint": "ProLIF/RDKit/MDAnalysis probe failed in the drugmate environment",
                "required_for_preflight": False,
            },
            "params": {
                "vicinity_cutoff": _p(6.5, "float", min=1.0, max=15.0),
                "interaction_types": _p(
                    [
                        "Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking",
                        "Anionic", "Cationic", "CationPi", "PiCation",
                        "VdWContact", "XBDonor", "XBAcceptor",
                        "MetalDonor", "MetalAcceptor",
                    ],
                    "list", item="str",
                    desc=(
                        "Explicit ProLIF 2.x pairwise interaction capability set. "
                        "WaterBridge uses a separate selected-water route."
                    ),
                ),
            },
        },
        "arpeggio": {
            "default_enabled": False,
            "purpose": (
                "Supplement PLIP/ProLIF with atom-level PDBe Arpeggio contacts "
                "from the untouched source complex mmCIF."
            ),
            "when_to_use": (
                "Holo structures with an RCSB/source mmCIF; degraded chemistry "
                "warnings remain observations and cannot promote hard constraints."
            ),
            "depends_on": [],
            "requires": ["reference_ligand", "source_complex_mmcif"],
            "runtime": {
                "kind": "command",
                "environment": "drugmate",
                "command": "python",
                "probe_code": (
                    "import shutil,arpeggio,gemmi; "
                    "from openbabel import openbabel; "
                    "assert shutil.which('pdbe-arpeggio')"
                ),
                "probe_failure_hint": (
                    "PDBe Arpeggio, Gemmi, OpenBabel, or pdbe-arpeggio CLI probe "
                    "failed in the drugmate environment"
                ),
            },
            "params": {},
        },
        "luna": {
            "default_enabled": False,
            "purpose": (
                "Supplement core consensus with atom-aware LUNA interactions "
                "and a hashed count fingerprint."
            ),
            "when_to_use": (
                "Independent interaction cross-check in the dedicated luna "
                "environment; PyMOL session generation is intentionally disabled."
            ),
            "depends_on": [],
            "requires": ["reference_ligand"],
            "runtime": {
                "kind": "command",
                "environment": "luna",
                "command": "python",
                "probe_code": (
                    "from luna.interaction.calc import InteractionCalculator; "
                    "from luna.interaction.filter import InteractionFilter; "
                    "from luna.mol.entry import MolFileEntry; "
                    "from luna.projects import LocalProject"
                ),
                "probe_failure_hint": (
                    "LUNA core API probe failed in the luna environment; the "
                    "run_luna/PyMOL CLI is not required"
                ),
            },
            "params": {
                "ifp_length": _p(4096, "int", min=128, max=65536),
                "ifp_type": _p(
                    "EIFP", "str", choices=["EIFP", "HIFP", "FIFP"],
                ),
                "add_dependent_interactions": _p(False, "bool"),
            },
        },
        "oddt": {
            "default_enabled": False,
            "purpose": (
                "Generate ODDT residue interaction counts plus SPLIF hashes as "
                "supporting fingerprint evidence."
            ),
            "when_to_use": (
                "Legacy independent cross-check in the dedicated oddt environment; "
                "never use residue-only observations or SPLIF as a hard gate."
            ),
            "depends_on": [],
            "requires": ["reference_ligand"],
            "runtime": {
                "kind": "command",
                "environment": "oddt",
                "command": "python",
                "probe_code": (
                    "import os; os.environ['ODDT_TOOLKIT']='ob'; "
                    "import numpy as np,oddt; from oddt import fingerprints; "
                    "assert hasattr(np,'in1d')"
                ),
                "probe_failure_hint": (
                    "ODDT/OpenBabel/NumPy compatibility probe failed in the oddt environment"
                ),
            },
            "params": {
                "strict": _p(True, "bool"),
                "splif_depth": _p(1, "int", min=0, max=4),
                "splif_size": _p(4096, "int", min=128, max=65536),
                "splif_distance_cutoff": _p(4.5, "float", min=1.0, max=10.0),
            },
        },
        "blueprint": {
            "default_enabled": True,
            "purpose": "Build 3-tier interaction constraints from eligible PLIP and curated knowledge; record ProLIF separately.",
            "when_to_use": "Needed for downstream blueprint_match hard gating.",
            # Curated known_required_interactions can define a valid blueprint
            # without either co-crystal interaction profiler.
            "depends_on": [],
            "requires": [],
            "runtime": {"kind": "internal"},
            "params": {
                "auto_required_from_tool_consensus": _p(
                    False, "bool",
                    desc=(
                        "Unsafe legacy promotion: permit tool consensus at a known-site "
                        "residue to become required without curated causal evidence."
                    ),
                ),
                "require_curated_provenance": _p(
                    True, "bool",
                    desc=(
                        "Mandatory safety invariant: require source, rationale, binding-mode "
                        "scope and approved review status before a curated interaction becomes "
                        "a hard constraint; pipeline validation rejects false."
                    ),
                ),
            },
        },
    },
}

STEPS["molecule_generation"] = {
    "order": 2,
    "skill": "molecule-generator",
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    "consumes": [
        "prepared_receptor_noH", "reference_ligand", "pocket_center",
        "pocket_definition", "interaction_blueprint",
    ],
    # Used only when context.target_name is absent. Keeping this separate from
    # the unconditional wiring set avoids invalidating an expensive library
    # when an unused chain description changes.
    "optional_consumes": ["inferred_target_name"],
    # Candidate sources have different prerequisites. Tool-level ``requires``
    # below enforces pocket artifacts only when Pocket2Mol is enabled.
    "requires_all": [],
    "produces": ["library", "library_manifest"],
    "params": {
        "exclude_scaffolds": _p([], "list", item="str"),
        "max_heavy_atoms": _p(None, "int", min=1, max=200),
        "max_mw": _p(None, "float", min=0.0, max=2000.0),
    },
    # Every member is a candidate source. Defaults intentionally remain
    # ChEMBL + Pocket2Mol; additional generators are opt-in query routes.
    "one_of": [[
        "chembl", "pocket2mol", "resgen", "diffsbdd", "targetdiff",
        "difflinker", "reinvent4", "lingo3dmol", "tamgen", "molcraft",
    ]],
    "tools": {
        "chembl": {
            "default_enabled": True,
            "purpose": "Retrieve known active compounds as screening candidates from ChEMBL bioactivity data.",
            "when_to_use": "Default evidence-backed candidate source; needs a resolvable target id.",
            "depends_on": [], "requires": [],
            "context_one_of": ["chembl_target_id", "target_name", "uniprot_id"],
            "params": {
                "activity_type": _p("IC50", "str"),
                "activity_threshold_nm": _p(1000.0, "float", min=0.0, max=1e6),
                "max_actives": _p(50, "int", min=1, max=5000),
            },
        },
        "pocket2mol": {
            "default_enabled": True,
            "purpose": "Structure-based de novo generation conditioned on the pocket.",
            "when_to_use": "When a prepared pocket PDB + pocket_center are available.",
            "depends_on": [], "requires": ["prepared_receptor_noH", "pocket_center"],
            "runtime": {
                "kind": "command",
                "environment": "pocket2mol",
                "command": "python",
                "source_dir": str(_POCKET2MOL_SOURCE_DIR),
                "checkpoint": str(_POCKET2MOL_CHECKPOINT),
                "asset_sha256": _POCKET2MOL_CHECKPOINT_SHA256,
                "probe_code": (
                    "from pathlib import Path;"
                    "import hashlib, torch, torch_geometric;"
                    f"assert Path({str(_POCKET2MOL_SOURCE_DIR / 'sample_for_pdb.py')!r}).is_file();"
                    f"p=Path({str(_POCKET2MOL_CHECKPOINT)!r});"
                    "assert p.is_file();"
                    f"assert hashlib.sha256(p.read_bytes()).hexdigest()=={_POCKET2MOL_CHECKPOINT_SHA256!r};"
                    "assert torch.cuda.is_available()"
                ),
                "probe_failure_hint": (
                    "Pocket2Mol environment, source/checkpoint, PyTorch-Geometric, "
                    "or CUDA is unavailable"
                ),
                "probe_timeout_seconds": 60,
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=10000),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "gpu_device": _p(0, "int", min=0, max=16),
                "taskset_cpu": _p(0, "int", min=-1, max=256),
            },
        },
        "resgen": {
            "default_enabled": False,
            "purpose": "Fast autoregressive 3D generation around a receptor pocket.",
            "when_to_use": (
                "Speed-oriented generation when a prepared receptor and reference "
                "ligand are available; complements or replaces Pocket2Mol."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH", "reference_ligand", "pocket_center"],
            "runtime": {
                "kind": "command",
                "environment": "resgen",
                "command": "python",
                "source_dir": str(_RESGEN_SOURCE_DIR),
                "entrypoint": str(_RESGEN_SOURCE_DIR / "gen.py"),
                "checkpoint": str(_RESGEN_CHECKPOINT),
                "asset_sha256": "45cac2f5186564e71b43fe43ce709ae2515a3a32847ef592d6f99810d32916df",
                "probe_code": (
                    "from pathlib import Path;"
                    "import torch,torch_geometric;"
                    "from rdkit import Chem;"
                    f"assert Path({str(_RESGEN_SOURCE_DIR / 'gen.py')!r}).is_file();"
                    f"assert Path({str(_RESGEN_CHECKPOINT)!r}).is_file();"
                    "assert torch.cuda.is_available()"
                ),
                "probe_failure_hint": (
                    "ResGen environment, source/checkpoint, PyTorch-Geometric, "
                    "RDKit, or CUDA is unavailable"
                ),
                "probe_timeout_seconds": 60,
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=10000),
                "beam_size": _p(500, "int", min=1, max=5000),
                "max_steps": _p(50, "int", min=1, max=200),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "gpu_device": _p(0, "int", min=0, max=16),
            },
        },
        "diffsbdd": {
            "default_enabled": False,
            "purpose": "Pocket-conditioned full-atom 3D diffusion generation.",
            "when_to_use": (
                "Explore a diverse pocket-conditioned 3D candidate set; slower "
                "than the speed route and subject to downstream pose QC."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH", "pocket_center"],
            "runtime": {
                "kind": "command",
                "environment": "diffsbdd",
                "command": "python",
                "source_dir": str(_DIFFSBDD_SOURCE_DIR),
                "entrypoint": str(_DIFFSBDD_SOURCE_DIR / "generate_ligands.py"),
                "checkpoint": str(_DIFFSBDD_CHECKPOINT),
                "asset_sha256": "07f86764bf569aafbc40a9c15fc02de8e2550437dd0f17f657eab3abe66c372c",
                "probe_code": (
                    "from pathlib import Path;"
                    "import torch,pytorch_lightning,torch_scatter;"
                    "from rdkit import Chem;"
                    f"assert Path({str(_DIFFSBDD_SOURCE_DIR / 'generate_ligands.py')!r}).is_file();"
                    f"assert Path({str(_DIFFSBDD_CHECKPOINT)!r}).is_file();"
                    "assert torch.cuda.is_available()"
                ),
                "probe_failure_hint": (
                    "DiffSBDD environment, source/checkpoint, torch-scatter, "
                    "RDKit, or CUDA is unavailable"
                ),
                "probe_timeout_seconds": 60,
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=10000),
                "batch_size": _p(20, "int", min=1, max=1000),
                "timesteps": _p(500, "int", min=1, max=500),
                "pocket_radius": _p(10.0, "float", min=3.0, max=20.0),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "gpu_device": _p(0, "int", min=0, max=16),
            },
        },
        "targetdiff": {
            "default_enabled": False,
            "purpose": "Pocket-conditioned 3D diffusion generation with TargetDiff.",
            "when_to_use": (
                "Diversity-focused 3D exploration when full 1000-step diffusion "
                "cost is acceptable; downstream pose QC is mandatory."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH", "pocket_center"],
            "runtime": {
                "kind": "command",
                "environment": "targetdiff",
                "command": "python",
                "source_dir": str(_TARGETDIFF_SOURCE_DIR),
                "module": "scripts.sample_for_pocket",
                "checkpoint": str(_TARGETDIFF_CHECKPOINT),
                "asset_sha256": "8977821339ba31b7103b46835ef4ddbf7f15f3f2efc5cb1f2a405e7905b09db6",
                "probe_code": (
                    "from pathlib import Path;"
                    "import torch,torch_geometric;"
                    "from rdkit import Chem;"
                    f"assert Path({str(_TARGETDIFF_SOURCE_DIR / 'scripts' / 'sample_for_pocket.py')!r}).is_file();"
                    f"assert Path({str(_TARGETDIFF_CHECKPOINT)!r}).is_file();"
                    "assert torch.cuda.is_available()"
                ),
                "probe_failure_hint": (
                    "TargetDiff GPU environment, source/checkpoint, "
                    "PyTorch-Geometric, RDKit, or CUDA is unavailable"
                ),
                "probe_timeout_seconds": 60,
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=10000),
                "batch_size": _p(20, "int", min=1, max=1000),
                "num_steps": _p(1000, "int", min=1, max=1000),
                "pocket_radius": _p(10.0, "float", min=3.0, max=20.0),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "gpu_device": _p(0, "int", min=0, max=16),
            },
        },
        "difflinker": {
            "default_enabled": False,
            "purpose": "Pocket-conditioned 3D linking of disconnected input fragments.",
            "when_to_use": (
                "Only for fragment-linking queries with a trusted 3D fragment "
                "file containing at least two disconnected components."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH"],
            "context_one_of": ["fragment_file"],
            "runtime": {
                "kind": "command",
                "environment": "difflinker",
                "command": "python",
                "source_dir": str(_DIFFLINKER_SOURCE_DIR),
                "entrypoint": str(_DIFFLINKER_SOURCE_DIR / "generate_with_protein.py"),
                "checkpoint": str(_DIFFLINKER_CHECKPOINT),
                "asset_sha256": "695ef7d634c09db4b9635c4e5efebad23385b8432f0687b6eae5e726ba7883f9",
                "probe_code": (
                    "from pathlib import Path;"
                    "import shutil,torch,pytorch_lightning,torch_scatter;"
                    "from rdkit import Chem;"
                    f"assert Path({str(_DIFFLINKER_SOURCE_DIR / 'generate_with_protein.py')!r}).is_file();"
                    f"assert Path({str(_DIFFLINKER_CHECKPOINT)!r}).is_file();"
                    "assert shutil.which('obabel');"
                    "assert torch.cuda.is_available()"
                ),
                "probe_failure_hint": (
                    "DiffLinker environment, source/checkpoint, torch-scatter, "
                    "Open Babel, RDKit, or CUDA is unavailable"
                ),
                "probe_timeout_seconds": 60,
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=10000),
                "linker_size": _p(5, "int", min=1, max=50),
                "n_steps": _p(None, "int", min=1, max=1000),
                "max_batch_size": _p(16, "int", min=1, max=1000),
                "anchor_indices": _p(None, "str"),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "gpu_device": _p(0, "int", min=0, max=16),
            },
        },
        "reinvent4": {
            "default_enabled": False,
            "purpose": "Ligand-only 2D de novo sampling from the reviewed PubChem prior.",
            "when_to_use": (
                "Add a fast pocket-agnostic diversity baseline or generate "
                "candidates when no structural conditioning is desired."
            ),
            "depends_on": [],
            "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "reinvent4",
                "command": "reinvent",
                "source_dir": str(_REINVENT4_SOURCE_DIR),
                "checkpoint": str(_REINVENT4_PRIOR),
                "asset_sha256": "fe8cd1678452ad292a8f93e97cb19a85959b729e17113e157180d5e69ae89ef3",
                "official_md5": "75052d5e93b2470e9c6684d369704f5f",
                "probe_code": (
                    "from pathlib import Path;"
                    "import importlib.util,torch;"
                    "from rdkit import Chem;"
                    f"assert Path({str(_REINVENT4_PRIOR)!r}).is_file();"
                    "assert importlib.util.find_spec('reinvent');"
                    "assert torch.cuda.is_available()"
                ),
                "probe_failure_hint": (
                    "REINVENT4 environment, reviewed PubChem prior, RDKit, "
                    "or CUDA is unavailable"
                ),
                "probe_timeout_seconds": 60,
            },
            "params": {
                "num_samples": _p(200, "int", min=1, max=100000),
                "unique_molecules": _p(True, "bool"),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "gpu_device": _p(0, "int", min=0, max=16),
            },
        },
        "lingo3dmol": {
            "default_enabled": False,
            "purpose": (
                "Generate pocket-conditioned molecules with Lingo3DMol and "
                "retain its native 3D coordinates in the canonical library."
            ),
            "when_to_use": (
                "Opt-in Step 2 source when an additional pocket-aware 3D "
                "generator is wanted; its outputs still pass the common "
                "identity, standardization, and downstream pose-QC contracts."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH", "pocket_center"],
            "runtime": {
                "kind": "command",
                "environment": "lingo3dmol",
                "command": "python",
                "source_dir": str(_LINGO3DMOL_SOURCE_DIR),
                "entrypoint": str(
                    _LINGO3DMOL_SOURCE_DIR / "inference"
                    / "inference_avoid_clash.py"
                ),
                "checkpoint": str(_LINGO3DMOL_GENERATOR_CHECKPOINT),
                "asset_sha256": (
                    "edc665e1b46ce1ff39d03a15c56b9c019ec44a8573acb51d"
                    "0698719ef87247ab"
                ),
                "auxiliary_checkpoint": str(_LINGO3DMOL_CONTACT_CHECKPOINT),
                "auxiliary_asset_sha256": (
                    "14fd26a256c8979cd4e4ae87c3d1575444fb88913475b002a"
                    "415a1544dad7990"
                ),
                "probe_code": (
                    "import sys,subprocess; from pathlib import Path; "
                    f"sys.path.insert(0,{str(_LINGO3DMOL_SOURCE_DIR)!r}); "
                    "import torch; "
                    "from rdkit import Chem; "
                    "from model.transformer_v1_res_mp1 import TransformerModel as ContactModel; "
                    "from model.transformer_v1_res_fac2 import TransformerModel as GeneratorModel; "
                    "from dataloader.dataloader_case_nci_res_merge import mydataset; "
                    f"entry=Path({str(_LINGO3DMOL_SOURCE_DIR / 'inference' / 'inference_avoid_clash.py')!r}); "
                    f"gen=Path({str(_LINGO3DMOL_GENERATOR_CHECKPOINT)!r}); "
                    f"contact=Path({str(_LINGO3DMOL_CONTACT_CHECKPOINT)!r}); "
                    "assert entry.is_file(); "
                    "assert subprocess.check_output(['sha256sum',str(gen)],text=True).split()[0]=="
                    "'edc665e1b46ce1ff39d03a15c56b9c019ec44a8573acb51d0698719ef87247ab'; "
                    "assert subprocess.check_output(['sha256sum',str(contact)],text=True).split()[0]=="
                    "'14fd26a256c8979cd4e4ae87c3d1575444fb88913475b002a415a1544dad7990'; "
                    "assert torch.cuda.is_available()"
                ),
                "probe_timeout_seconds": 120,
                "probe_failure_hint": (
                    "Lingo3DMol source, pinned contact/generation checkpoints, "
                    "RDKit, PyTorch, or CUDA is unavailable."
                ),
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=1000),
                "gpu_device": _p(0, "int", min=0, max=16),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "max_run_hours": _p(1, "int", min=1, max=24),
                "pocket_radius": _p(10.0, "float", min=3.0, max=12.0),
                "gen_frag_set": _p(40, "int", min=5, max=500),
                "prod_time": _p(1, "int", min=1, max=100),
                "topk": _p(5, "int", min=1, max=50),
                "nci_threshold": _p(0.7, "float", min=0.0, max=1.0),
                "collision_distance": _p(2.5, "float", min=0.5, max=5.0),
                "frag_len_add": _p(15, "int", min=0, max=50),
            },
        },
        "tamgen": {
            "default_enabled": False,
            "purpose": (
                "Generate pocket-conditioned SMILES with TamGen from a "
                "prepared receptor and the Step 1 pocket center."
            ),
            "when_to_use": (
                "Opt-in sequence/coordinate generation route. The default "
                "unconditional decoder avoids treating TamGen's placeholder "
                "target-side molecule as a real chemical seed."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH", "pocket_center"],
            "runtime": {
                "kind": "command",
                "environment": "tamgen",
                "command": "python",
                "source_dir": str(_TAMGEN_SOURCE_DIR),
                "entrypoint": str(_TAMGEN_SOURCE_DIR / "generate.py"),
                "checkpoint": str(_TAMGEN_CHECKPOINT),
                "asset_sha256": (
                    "526c6337132b93050eb2b23cd8177175ad3d7d2310177536"
                    "a1dca1949ea45b4b"
                ),
                "probe_code": (
                    "import sys,subprocess; from pathlib import Path; "
                    f"sys.path.insert(0,{str(_TAMGEN_SOURCE_DIR)!r}); "
                    "import torch,numpy; from Bio.PDB import PDBParser; "
                    "from rdkit import Chem; from fairseq import checkpoint_utils,tasks; "
                    "from fairseq.molecule_utils.external.fairseq_dataset_build_utils import dump_center_data; "
                    f"entry=Path({str(_TAMGEN_SOURCE_DIR / 'generate.py')!r}); "
                    f"model=Path({str(_TAMGEN_CHECKPOINT)!r}); "
                    f"root=Path({str(_TAMGEN_SOURCE_DIR)!r}); "
                    "assert entry.is_file() and (root/'dict/dict.tg.txt').is_file() "
                    "and (root/'dict/dict.m1.txt').is_file(); "
                    "assert subprocess.check_output(['sha256sum',str(model)],text=True).split()[0]=="
                    "'526c6337132b93050eb2b23cd8177175ad3d7d2310177536a1dca1949ea45b4b'; "
                    "assert torch.cuda.is_available()"
                ),
                "probe_timeout_seconds": 180,
                "probe_failure_hint": (
                    "TamGen source/checkpoint, Fairseq/BioPython/RDKit, or CUDA "
                    "is unavailable."
                ),
            },
            "params": {
                "num_samples": _p(100, "int", min=1, max=1000),
                "beam": _p(20, "int", min=1, max=100),
                "max_seeds": _p(5, "int", min=1, max=100),
                "pocket_radius": _p(10.0, "float", min=3.0, max=20.0),
                "sample_beta": _p(1.0, "float", min=0.01, max=10.0),
                "use_conditional": _p(
                    False, "bool",
                    desc=(
                        "Expert-only target-side seed conditioning. Keep false "
                        "for center-built inputs whose target token is a stub."
                    ),
                ),
                "gpu_device": _p(0, "int", min=0, max=16),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "molcraft": {
            "default_enabled": False,
            "purpose": (
                "Sample pocket-conditioned 3D molecules with MolCRAFT and "
                "preserve the native coordinates in the canonical library."
            ),
            "when_to_use": (
                "Opt-in 3D generation when a prepared receptor and reference "
                "ligand define the installed local-sampling route."
            ),
            "depends_on": [],
            "requires": ["prepared_receptor_noH", "reference_ligand"],
            "runtime": {
                "kind": "command",
                "environment": "molcraft",
                "command": "python",
                "source_dir": str(_MOLCRAFT_SOURCE_DIR),
                "entrypoint": str(_MOLCRAFT_SOURCE_DIR / "sample_for_pocket_local.py"),
                "checkpoint": str(_MOLCRAFT_CHECKPOINT),
                "asset_sha256": (
                    "0a54f9bfb847ea87e7aec538df90a009c8c49d12a7b5d099"
                    "47d04b7f9444e024"
                ),
                "checkpoint_config": str(_MOLCRAFT_CHECKPOINT_CONFIG),
                "checkpoint_config_sha256": (
                    "efd77b3dbc49a76ebdaf1ab00b0e5f043c4c0790bec6d046"
                    "bec604c95dab8805"
                ),
                "probe_code": (
                    "import sys,subprocess; from pathlib import Path; "
                    f"sys.path.insert(0,{str(_MOLCRAFT_SOURCE_DIR)!r}); "
                    "import torch,torch_geometric,pytorch_lightning,yaml; "
                    "from rdkit import Chem; import posecheck; "
                    "from core.models.sbdd_train_loop import SBDDTrainLoop; "
                    "from core.callbacks.validation_callback_for_sample import DockingTestCallback; "
                    f"entry=Path({str(_MOLCRAFT_SOURCE_DIR / 'sample_for_pocket_local.py')!r}); "
                    f"model=Path({str(_MOLCRAFT_CHECKPOINT)!r}); "
                    f"config=Path({str(_MOLCRAFT_CHECKPOINT_CONFIG)!r}); "
                    "assert entry.is_file(); "
                    "assert subprocess.check_output(['sha256sum',str(model)],text=True).split()[0]=="
                    "'0a54f9bfb847ea87e7aec538df90a009c8c49d12a7b5d09947d04b7f9444e024'; "
                    "assert subprocess.check_output(['sha256sum',str(config)],text=True).split()[0]=="
                    "'efd77b3dbc49a76ebdaf1ab00b0e5f043c4c0790bec6d046bec604c95dab8805'; "
                    "assert torch.cuda.is_available()"
                ),
                "probe_timeout_seconds": 120,
                "probe_failure_hint": (
                    "MolCRAFT local entrypoint, pinned checkpoint/config, "
                    "RDKit, PyTorch, YAML, or CUDA is unavailable."
                ),
            },
            "params": {
                "num_samples": _p(10, "int", min=1, max=1000),
                "sample_steps": _p(100, "int", min=1, max=1000),
                "sample_num_atoms": _p(
                    "prior", "str", choices=["prior", "ref"],
                ),
                "beta1": _p(1.5, "float", min=0.1, max=10.0),
                "sigma1_coord": _p(0.03, "float", min=0.001, max=1.0),
                "sampling_strategy": _p(
                    "end_back", "str",
                    choices=["vanilla", "end_back", "end_back_pmf", "end_back_mode"],
                ),
                "gpu_device": _p(0, "int", min=0, max=16),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "dude_decoys": {
            "default_enabled": True,
            "purpose": "Download DUD-E property-matched decoys.",
            "when_to_use": "Prefer cached/downloaded decoys for retrospective benchmarking. If unavailable or download fails, record the reason and continue without decoys; do not enable local generation automatically.",
            "depends_on": [], "requires": [],
            "params": {"num_decoys": _p(200, "int", min=0, max=10000)},
        },
        "rdkit_decoys": {
            "default_enabled": False,
            "purpose": "Optional local BRICS generation of property-matched decoys.",
            "when_to_use": "Only when local decoy generation is explicitly requested; unavailable DUD-E data alone does not enable this tool.",
            "depends_on": [], "requires": [],
            "params": {"num_decoys": _p(200, "int", min=0, max=10000)},
        },
        "standardization": {
            "default_enabled": True,
            "purpose": "Salt strip, charge handling, canonicalization.",
            "when_to_use": "Almost always before conformer generation.",
            "depends_on": [], "requires": [],
            "params": {"neutralize": _p(False, "bool")},
        },
        "conformer": {
            "default_enabled": True,
            "purpose": "Preserve generator-native 3D poses; otherwise generate ETKDGv3 + MMFF/UFF conformers.",
            "when_to_use": "Always; docking needs 3D coordinates and generator poses need provenance.",
            "depends_on": [], "requires": [],
            "params": {},
        },
    },
}


def _screening_mode_param() -> dict:
    """Declare the level-2 screening mode for a step that folds metric marks.

    Level 1 marks each metric ``pass``/``flag``/``fail``; level 2 counts the
    ``fail`` marks inside every parent metric category and only there decides
    keep/modify/reject. This parameter selects the level-2
    thresholds.  ``None`` means "auto": the step derives the mode from
    ``context.candidate_source`` -- generated molecules get the strict V3, a
    sampled library V2, and a user-supplied full set the marking-only V1 --
    falling back to V2 when the source is unknown.
    """
    return _p(
        None,
        "str",
        choices=sorted(SCREENING_MODES),
        desc="; ".join(
            f"{mode_id}={spec['name']} "
            f"(count(fail)>={spec['fail_count_threshold']} or "
            f"ratio(fail)>={spec['fail_ratio_threshold']:.0%} rejects)"
            if spec["eliminates"]
            else f"{mode_id}={spec['name']} (marks only, never eliminates)"
            for mode_id, spec in SCREENING_MODES.items()
        ),
    )


def _property_tools(
    admet_default: bool,
    chemfh_default: bool,
    dedup_default: bool,
    aizynth_default: bool,
) -> dict:
    """Shared tool block for property_prefilter and property_filter steps."""
    return {
        "lipinski": {
            "default_enabled": True,
            "purpose": "Lipinski Rule of Five violation check.",
            "when_to_use": "Default drug-likeness gate.",
            "depends_on": [], "requires": [],
            # RO5 bounds: 0-1 pass, 2-3 flag,
            # >3 fail.
            "params": {"lipinski_max_violations": _p(3, "int", min=0, max=4)},
        },
        "qed": {
            "default_enabled": True,
            "purpose": "Quantitative estimate of drug-likeness.",
            "when_to_use": "Default soft-flag on low QED.",
            "depends_on": [], "requires": [],
            "params": {"qed_threshold": _p(0.45, "float", min=0.0, max=1.0)},
        },
        "rascore": {
            "default_enabled": False,
            "purpose": "CPU RAscore XGB accessibility probability in the SYN category.",
            "when_to_use": "Low-cost learned accessibility check before retrosynthesis.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command", "environment": "rascore", "command": "python",
                "probe_code": "from RAscore.RAscore_XGB import RAScorerXGB; RAScorerXGB()",
            },
            "params": {"timeout_seconds": _p(600, "int", min=30, max=7200)},
        },
        "sa_score": {
            "default_enabled": True,
            "purpose": "Synthetic accessibility score.",
            "when_to_use": "Default soft-flag on hard-to-make molecules.",
            "depends_on": [], "requires": [],
            "params": {"sa_threshold": _p(6.0, "float", min=1.0, max=10.0)},
        },
        "pains": {
            "default_enabled": True,
            "purpose": "All-match PAINS A/B/C, Brenk, and NIH structural-alert screen.",
            "when_to_use": (
                "Default first pass. Each PAINS-family, Brenk, and NIH count uses "
                "the sbdd numeric band: 0 pass, 1-4 flag, >4 fail; the resulting "
                "marks are folded in ALERT before any molecule is rejected. "
                "Covalent status is not an automatic rescue; an explicit audited "
                "catalog policy is required."
            ),
            "depends_on": [], "requires": [],
            "params": {
                "catalogs": _p(
                    ["PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"],
                    "list",
                    item="str",
                    item_choices=["PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"],
                ),
                "hard_fail_catalogs": _p(
                    ["PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"],
                    "list",
                    item="str",
                    item_choices=["PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH"],
                ),
            },
        },
        "admet": {
            "default_enabled": admet_default,
            "purpose": "Strict batch ADMET-AI predictions with endpoint provenance.",
            "when_to_use": (
                "On for final filtering and off for pre-docking prefilter. "
                "Missing/incomplete predictions cap a molecule at review_only."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "admet_ai",
                "command": "python",
                "probe_timeout_seconds": 90,
                "probe_code": (
                    "import os;"
                    "os.environ['CUDA_VISIBLE_DEVICES']='-1';"
                    "from admet_ai import ADMETModel;"
                    "ADMETModel()"
                ),
                "probe_failure_hint": (
                    "ADMET-AI model initialization failed in the admet_ai "
                    "environment (CPU mode)."
                ),
            },
            "params": {
                "timeout_seconds": _p(7200, "int", min=30, max=7200),
                "bbb_profile": _p(
                    "non_cns",
                    "str",
                    choices=["non_cns", "cns"],
                    desc=(
                        "BBB_Martins direction profile. non_cns penalizes high "
                        "penetration; cns rewards high penetration."
                    ),
                ),
            },
        },
        "chemfh": {
            "default_enabled": chemfh_default,
            "purpose": (
                "Strict batch ChemFH frequent-hitter and assay-interference "
                "predictions with all 17 native rule results."
            ),
            "when_to_use": (
                "On for final filtering and off for the pre-docking prefilter. "
                "Seven table-matched risk probabilities use 0.30/0.90 level-1 "
                "ALERT bounds; their signed weights follow the metric catalog. "
                "Incomplete evidence is a protocol review condition, not a "
                "molecule fail. Native structural-rule hits remain evidence only; "
                "ChemFH PAINS does not duplicate the RDKit gate."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "chemfh",
                "command": "python",
                "source_dir": str(_CHEMFH_SOURCE_DIR),
                "source_commit": _CHEMFH_SOURCE_COMMIT,
                "checkpoint": str(_CHEMFH_CHECKPOINT),
                "checkpoint_sha256": _CHEMFH_CHECKPOINT_SHA256,
                "rule_bundle_sha256": _CHEMFH_RULE_BUNDLE_SHA256,
                "probe_timeout_seconds": 180,
                "probe_args": [
                    "-B", str(_CHEMFH_WORKER),
                    "--probe",
                    "--source-dir", str(_CHEMFH_SOURCE_DIR),
                    "--checkpoint", str(_CHEMFH_CHECKPOINT),
                    "--expected-checkpoint-sha256", _CHEMFH_CHECKPOINT_SHA256,
                    "--expected-rule-bundle-sha256", _CHEMFH_RULE_BUNDLE_SHA256,
                    "--expected-source-commit", _CHEMFH_SOURCE_COMMIT,
                ],
                "probe_failure_hint": (
                    "ChemFH source, checkpoint, rule bundle, or functional "
                    "CPU prediction probe failed in the configured chemfh environment."
                ),
            },
            "params": {
                "timeout_seconds": _p(7200, "int", min=30, max=7200),
                "seed": _p(42, "int", min=0, max=2147483647),
            },
        },
        "aizynthfinder": {
            "default_enabled": aizynth_default,
            "purpose": (
                "Top-K retrosynthesis planning with AiZynthFinder after RDKit "
                "and complete ADMET assessment."
            ),
            "when_to_use": (
                "Default for final candidates only; never run across the full "
                "pre-docking library."
            ),
            "depends_on": ["admet"], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "aizynth",
                "command": "python",
                "config_path": str(_AIZYNTHFINDER_CONFIG),
                "probe_timeout_seconds": 120,
                "probe_code": (
                    "import sys;"
                    f"sys.path.insert(0,{str(_AIZYNTHFINDER_PROBE.parent)!r});"
                    "from aizynth_probe import probe;"
                    f"probe({str(_AIZYNTHFINDER_CONFIG)!r})"
                ),
                "probe_failure_hint": (
                    "AiZynthFinder config, expansion/filter policy, or stock "
                    "failed to load in the aizynth environment."
                ),
            },
            "params": {
                "top_k": _p(25, "int", min=1, max=500),
                "timeout_seconds": _p(21600, "int", min=30, max=21600),
                "search_time_limit_seconds": _p(60, "int", min=1, max=600),
                "iteration_limit": _p(100, "int", min=1, max=1000),
                "max_routes": _p(5, "int", min=1, max=50),
                "return_first": _p(True, "bool"),
                "expansion_policies": _p(
                    ["uspto"], "list", item="str",
                    item_choices=["uspto", "ringbreaker"],
                ),
                "filter_policies": _p(
                    ["uspto"], "list", item="str", item_choices=["uspto"],
                ),
                "stocks": _p(
                    ["zinc"], "list", item="str", item_choices=["zinc"],
                ),
            },
        },
        "dedup": {
            "default_enabled": dedup_default,
            "purpose": "Collapse repeated docking poses to one record per compound.",
            "when_to_use": "On for post-docking filter; off for raw library prefilter.",
            "depends_on": [], "requires": [],
            "params": {},
        },
    }


STEPS["property_prefilter"] = {
    "order": 3,
    "skill": "drug-property-filter",
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    "consumes": ["library"],
    "produces": ["prefiltered_library"],
    "params": {
        # "not_evaluated" is carried forward too: it means the RDKit checks were
        # switched off, not that the molecule failed, so dropping it here would
        # silently empty the prefiltered library. Its classification stays
        # "modify", so it still cannot be nominated.  "fail" is likewise only a
        # level-1 mark on one metric; the skill's own two-level classification
        # is what removes a molecule.
        "filtered_verdicts": _p(
            ["pass", "flag", "fail", "not_evaluated"], "list", item="str"
        ),
        "fail_if_no_filtered_molecules": _p(True, "bool"),
        "screening_mode": _screening_mode_param(),
    },
    "tools": _property_tools(
        admet_default=False,
        chemfh_default=False,
        dedup_default=False,
        aizynth_default=False,
    ),
}

STEPS["docking"] = {
    "order": 4,
    "skill": "docking",
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    # `consumes` is the wiring superset (orchestrator wires whatever is present).
    # Validation uses requires_all (hard) + requires_any (at least one per group).
    "consumes": ["prepared_receptor_noH", "prep_manifest", "source_complex_pdb", "water_inventory",
                 "interaction_blueprint",
                 "prefiltered_library", "library", "reference_ligand",
                 "pocket_definition", "docking_box", "library_manifest",
                 "dock6_receptor_site",
                 "dock6_grid_nrg", "dock6_grid_bmp"],
    # External/prepared libraries may carry the same state lineage directly in
    # their SDF properties.  Consume and verify a manifest when supplied, but
    # do not make an optional molecule-generation sidecar a hard DAG input.
    "requires_all": ["prepared_receptor_noH"],
    "requires_any": [
        ["prefiltered_library", "library"],
    ],
    "produces": ["docked_poses", "docking_scores", "docking_report",
                 "prepared_receptor", "receptor_pdbqt", "receptor_state_manifest",
                 "water_report", "docking_species_library", "microstate_identity_map",
                 "docking_prep_qc", "reference_ligand_microstates"],
    "params": {
        "search_space_source": _p(
            None,
            "str",
            choices=[
                "reference_ligand",
                "docking_box",
                "whole_receptor",
                "precomputed_grid",
            ],
            required=True,
            desc=(
                "Explicit docking search-space source. The runner never chooses "
                "between a reference-ligand autobox and an explicit box."
            ),
        ),
    },
    # Until staged screening/refinement and cross-engine score fusion have an
    # explicit contract, Step 4 executes exactly one primary pose engine per
    # run. ``exactly_one`` already includes the old ``one_of`` constraint.
    "exactly_one": [list(DOCKING_PRIMARY_ENGINES), ["reduce", "pdb2pqr", "protoss"]],
    "tools": {
        "docking_prep": {
            "default_enabled": True,
            "purpose": "Mandatory final receptor/water/ligand-state preparation before the engine.",
            "when_to_use": "Every docking run, after property_prefilter and pocket selection.",
            "depends_on": [], "requires": [],
            "runtime": {"kind": "internal"},
            "params": {
                "ph": _p(None, "float", min=0.0, max=14.0),
                "reason": _p(None, "str", required=True, desc="Evidence for the selected protonation/flip backend and pH, when applicable."),
                "water_policy": _p(None, "str", required=True, choices=["dry", "selected"]),
                "water_reason": _p(None, "str", required=True),
                "ligand_state_policy": _p("preserve_declared", "str", choices=["preserve_declared"]),
                "ligand_state_reason": _p(None, "str", required=True, desc="Why the prefiltered exact ligand states are appropriate for this site."),
            },
        },
        "reduce": {
            "default_enabled": False,
            "purpose": "Add hydrogens and optimize His/Asn/Gln flips with mmtbx.reduce2.",
            "when_to_use": (
                "Preferred final preparation backend for routine non-covalent protein receptors "
                "when target-pH protonation assignment is not required, unless the user specifies "
                "another backend. Explicitly enable build=true and flip=true with a passing runtime "
                "preflight. Does not apply a target pH; never switch backends silently after failure."
            ),
            "depends_on": ["docking_prep"], "requires": [],
            "produces": ["receptor_state_manifest"],
            "runtime": {
                "kind": "command", "environment": "reduce",
                "command": "mmtbx.reduce2", "probe_args": ["--version"],
                "probe_code": (
                    "import shutil,sys;"
                    "from mmtbx.monomer_library import server;"
                    "from mmtbx import chemical_components;"
                    "from mmtbx.ligands.rdkit_utils import get_molecule_from_resname;"
                    "server.server();"
                    "assert chemical_components.get_cif_dictionary('ALA') is not None;"
                    "mol=get_molecule_from_resname('ALA');"
                    "assert mol is not None and mol.GetNumAtoms()>0, 'ALA CCD molecule could not be constructed';"
                    "sys.exit(0 if shutil.which('mmtbx.reduce2') else 1)"
                ),
                "probe_failure_hint": (
                    "mmtbx.reduce2, its CCP4/geostd monomer library, or the CCTBX "
                    "CCD chemical-components database is unavailable or incompatible. "
                    "Configure MMTBX_CCP4_MONOMER_LIB or CLIBD_MON for CCP4/geostd "
                    "restraints and a separate genuine CCD at the libtbx "
                    "chemical_components path; geostd data_*.cif files are not CCD files."
                ),
            },
            "params": {"build": _p(True, "bool"), "flip": _p(True, "bool")},
        },
        "pdb2pqr": {
            "default_enabled": False,
            "purpose": "Apply pH-aware receptor protonation with PDB2PQR+PROPKA and optimize hydrogen bonds and His/Asn/Gln flips.",
            "when_to_use": (
                "Choose for target-pH protonation assignment, protonation-sensitive protocols, "
                "or an explicit user request for PDB2PQR+PROPKA. "
                "Explicitly enable with justified docking_prep.params.ph and optimize_hbond_network=true. "
                "Check support for retained metals, cofactors and nonstandard residues before selection; "
                "no universal pH default or automatic runtime fallback."
            ),
            "depends_on": ["docking_prep"], "requires": [],
            "produces": ["receptor_state_manifest"],
            "runtime": {
                "kind": "command", "environment": "drugmate",
                "command": "pdb2pqr30", "probe_args": ["--version"],
            },
            "params": {
                "force_field": _p("AMBER", "str", choices=["AMBER", "CHARMM", "PARSE"]),
                "optimize_hbond_network": _p(True, "bool"),
            },
        },
        "protoss": {
            "default_enabled": False,
            "purpose": "Use ProteinsPlus Protoss for complex-aware protonation and flips.",
            "when_to_use": "Unavailable until prepared-complex upload preserves input identity; do not use PDB-ID substitution.",
            "depends_on": ["docking_prep"], "requires": [],
            "produces": ["receptor_state_manifest"],
            "runtime": {"kind": "http", "environment": "proteins_plus", "endpoint": "/protoss_rest"},
            "params": {
                "poll_interval_seconds": _p(5, "int", min=1, max=60),
                "timeout_seconds": _p(600, "int", min=30, max=3600),
            },
        },
        "water_selection": {
            "default_enabled": False,
            "purpose": "Select ligand-bridging crystallographic waters and emit an auditable report.",
            "when_to_use": "When a retained reference ligand can define bridging waters, or with manual water residue ids.",
            "depends_on": ["docking_prep"], "requires": [],
            "runtime": {"kind": "internal"},
            "params": {
                "ligand_distance": _p(3.5, "float", min=0.0, max=10.0),
                "protein_distance": _p(3.5, "float", min=0.0, max=10.0),
                "min_occupancy": _p(0.5, "float", min=0.0, max=1.0),
                "max_b_factor": _p(60.0, "float", min=0.0, max=500.0),
                "include_residues": _p([], "list", item="str", desc="Manual chain:resSeq water keeps."),
                "exclude_residues": _p([], "list", item="str", desc="Manual chain:resSeq water removals."),
            },
        },
        "meeko": {
            "default_enabled": False,
            "purpose": "Generate the receptor PDBQT consumed by the AutoDock Vina adapter.",
            "when_to_use": (
                "Enable for the autodock_vina route after final receptor protonation. AutoDock4 and ADFR use adapter-owned native preparation; "
                "GNINA PDB+SDF does not need Meeko."
            ),
            "depends_on": ["docking_prep"], "requires": [],
            "produces": ["receptor_pdbqt"],
            "runtime": {
                "kind": "command", "environment": "reduce",
                "command": "mk_prepare_receptor.py", "probe_args": ["--help"],
                "probe_code": (
                    "import shutil,sys;"
                    "from meeko import MoleculePreparation;"
                    "MoleculePreparation();"
                    "sys.exit(0 if shutil.which('mk_prepare_receptor.py') else 1)"
                ),
                "probe_failure_hint": (
                    "Meeko receptor preparation is not functional in the reduce environment; "
                    "check the command and bundled atom-parameter resources."
                ),
            },
            "params": {"allow_bad_residues": _p(False, "bool")},
        },
        "ligand_microstates": {
            "default_enabled": False,
            "purpose": "Enumerate reference-ligand protonation microstates with Dimorphite-DL.",
            "when_to_use": "Ionizable reference ligands or protonation-sensitive docking studies.",
            "depends_on": ["docking_prep"], "requires": [],
            "produces": ["reference_ligand_microstates"],
            "runtime": {
                "kind": "command", "environment": "drugmate",
                "command": "dimorphite_dl", "probe_args": ["--help"],
            },
            "params": {
                "ph_tolerance": _p(1.0, "float", min=0.0, max=4.0),
                "max_variants": _p(16, "int", min=1, max=128),
            },
        },
        "openbabel": {
            "default_enabled": False,
            "purpose": "Contextual format conversion fallback through Open Babel.",
            "when_to_use": "Only when a downstream format is unsupported by RDKit/Meeko.",
            "depends_on": ["docking_prep"], "requires": [],
            "runtime": {
                "kind": "command", "environment": "drugmate",
                "command": "obabel", "probe_args": ["-V"],
            },
            "params": {"output_format": _p("mol2", "str", choices=["mol2", "pdb", "sdf"])},
        },
        "gnina": {
            "default_enabled": True,
            "purpose": "GNINA docking with CNN scoring + Vina energy.",
            "when_to_use": "Default engine; use for routine docking and CNN rescoring.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "drugmate",
                "command": "gnina",
                "probe_args": ["--version"],
                "probe_failure_hint": "GNINA binary is unavailable from the drugmate runtime.",
            },
            "params": {
                "autobox_add": _p(4.0, "float", min=0.0, max=30.0),
                "exhaustiveness": _p(16, "int", min=1, max=64),
                "num_modes": _p(20, "int", min=1, max=100),
                "cnn_scoring": _p("rescore", "str",
                                  choices=["none", "rescore", "refinement", "all"]),
                "pose_sort_order": _p("CNNscore", "str",
                                      choices=["CNNscore", "CNNaffinity", "Energy"]),
                "device": _p(0, "int", min=0, max=16),
                "seed": _p(42, "int"),
                "gpu_monitor_interval_seconds": _p(5, "int", min=1, max=600),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "unidock": {
            "default_enabled": False,
            "purpose": (
                "GPU-batched Uni-Dock through its topology-aware "
                "unidocktools pipeline."
            ),
            "when_to_use": (
                "Large non-covalent libraries when throughput is more important "
                "than GNINA CNN scoring; current adapter is validated for "
                "shard-sized runs, not a monolithic million-molecule job."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "unidock",
                "command": "unidocktools",
                "probe_code": (
                    "import shutil,sys; "
                    "sys.exit(0 if all(shutil.which(name) for name in "
                    "('unidocktools','unidock')) else 1)"
                ),
                "probe_failure_hint": (
                    "Uni-Dock Tools or the Uni-Dock executable is unavailable "
                    "from the unidock runtime."
                ),
            },
            "params": {
                "autobox_add": _p(4.0, "float", min=0.0, max=30.0),
                "batch_size": _p(18000, "int", min=1, max=1000000),
                "scoring_function": _p(
                    "vina", "str", choices=["vina", "vinardo"],
                ),
                "search_mode": _p(
                    "", "str", choices=["", "fast", "balance", "detail"],
                ),
                "exhaustiveness": _p(16, "int", min=1, max=1024),
                "max_step": _p(20, "int", min=1, max=1000),
                "num_modes": _p(20, "int", min=1, max=100),
                "refine_step": _p(3, "int", min=1, max=100),
                "energy_range": _p(3.0, "float", min=0.0, max=100.0),
                "device": _p(0, "int", min=0, max=16),
                "seed": _p(42, "int"),
                "timeout_seconds": _p(86400, "int", min=1, max=86400),
            },
        },
        "autodock_vina": {
            "default_enabled": False,
            "purpose": (
                "Classic AutoDock Vina CPU docking through the Python API, "
                "with validated Meeko ligand preparation and SDF export."
            ),
            "when_to_use": (
                "Small-to-medium non-covalent libraries, a reproducible classic "
                "Vina baseline, or an explicit user request. It is not a strongly "
                "orthogonal validator for GNINA/Uni-Dock because all are in the "
                "Vina method family."
            ),
            "depends_on": [],
            "requires": [],
            "runtime": {
                "kind": "composite",
                "components": {
                    "engine": {
                        "kind": "command",
                        "environment": "vina",
                        "command": "python",
                        "probe_code": (
                            "from vina import Vina;"
                            "Vina(sf_name='vina',cpu=1,seed=42,verbosity=0)"
                        ),
                        "probe_failure_hint": (
                            "The AutoDock Vina Python API is not functional in "
                            "the vina environment."
                        ),
                    },
                    "io": {
                        "kind": "command",
                        "environment": "reduce",
                        "command": "python",
                        "probe_code": (
                            "from rdkit import Chem;"
                            "from meeko import MoleculePreparation,"
                            "PDBQTMolecule,RDKitMolCreate;"
                            "MoleculePreparation()"
                        ),
                        "probe_failure_hint": (
                            "RDKit/Meeko ligand preparation or export is not "
                            "functional in the reduce environment."
                        ),
                    },
                },
            },
            "params": {
                "autobox_add": _p(4.0, "float", min=0.0, max=30.0),
                "scoring_function": _p(
                    "vina", "str", choices=["vina", "vinardo"],
                ),
                "exhaustiveness": _p(16, "int", min=1, max=1024),
                "num_modes": _p(20, "int", min=1, max=100),
                "energy_range": _p(3.0, "float", min=0.01, max=100.0),
                "min_rmsd": _p(1.0, "float", min=0.0, max=20.0),
                "max_evals": _p(0, "int", min=0),
                "cpu": _p(1, "int", min=1, max=256),
                "seed": _p(42, "int"),
                "spacing": _p(0.375, "float", min=0.01, max=2.0),
                "force_even_voxels": _p(False, "bool"),
                "no_refine": _p(False, "bool"),
                "timeout_seconds": _p(86400, "int", min=1, max=86400),
            },
        },
        "diffdock": {
            "default_enabled": False,
            "purpose": (
                "Run generative DiffDock pose prediction and normalize every "
                "ranked pose plus confidence into the canonical docking contract."
            ),
            "when_to_use": (
                "Select as the one primary Step 4 engine for a generative "
                "docking route. Its confidence is method-specific and is not "
                "treated as a Vina-family binding energy."
            ),
            "depends_on": [],
            "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "diffdock",
                "command": "python",
                "source_dir": str(_DIFFDOCK_SOURCE_DIR),
                "score_checkpoint": str(_DIFFDOCK_SCORE_CHECKPOINT),
                "confidence_checkpoint": str(_DIFFDOCK_CONFIDENCE_CHECKPOINT),
                "asset_sha256": {
                    "score_checkpoint": (
                        "64afbdc97f6c700df0309162eb533ee6f2754c71412aca7a"
                        "3a534a7239667260"
                    ),
                    "confidence_checkpoint": (
                        "18d0fcc97fc55436786db02305ac61c734ad682aa5f8dee6"
                        "01555369c8461a57"
                    ),
                },
                "probe_code": (
                    "from pathlib import Path; import subprocess,torch,torch_geometric,e3nn; "
                    "from rdkit import Chem; "
                    f"root=Path({str(_DIFFDOCK_SOURCE_DIR)!r}); "
                    f"score=Path({str(_DIFFDOCK_SCORE_CHECKPOINT)!r}); "
                    f"confidence=Path({str(_DIFFDOCK_CONFIDENCE_CHECKPOINT)!r}); "
                    "required=[root/'inference.py',"
                    "root/'workdir/v1.1/score_model/model_parameters.yml',"
                    "root/'workdir/v1.1/confidence_model/model_parameters.yml']; "
                    "assert all(path.is_file() for path in required); "
                    "assert subprocess.check_output(['sha256sum',str(score)],text=True).split()[0]=="
                    "'64afbdc97f6c700df0309162eb533ee6f2754c71412aca7a3a534a7239667260'; "
                    "assert subprocess.check_output(['sha256sum',str(confidence)],text=True).split()[0]=="
                    "'18d0fcc97fc55436786db02305ac61c734ad682aa5f8dee601555369c8461a57'"
                ),
                "probe_timeout_seconds": 120,
                "probe_failure_hint": (
                    "DiffDock source, pinned score/confidence models, geometric "
                    "dependencies, or RDKit is unavailable. The adapter may "
                    "use DiffDock's CPU fallback when CUDA is not exposed."
                ),
            },
            "params": {
                "source_root": _p(str(_DIFFDOCK_SOURCE_DIR), "str"),
                "model_dir": _p(
                    str(_DIFFDOCK_SCORE_CHECKPOINT.parent), "str",
                ),
                "confidence_model_dir": _p(
                    str(_DIFFDOCK_CONFIDENCE_CHECKPOINT.parent), "str",
                ),
                "checkpoint": _p(
                    _DIFFDOCK_SCORE_CHECKPOINT.name, "str",
                ),
                "confidence_checkpoint": _p(
                    _DIFFDOCK_CONFIDENCE_CHECKPOINT.name, "str",
                ),
                "samples_per_complex": _p(20, "int", min=1, max=40),
                "batch_size": _p(10, "int", min=1, max=64),
                "inference_steps": _p(20, "int", min=1, max=100),
                "actual_steps": _p(None, "int", min=1, max=100),
                "old_confidence_model": _p(True, "bool"),
                "device": _p(0, "int", min=0, max=16),
                "seed": _p(42, "int"),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "dock6": {
            "default_enabled": False,
            "purpose": (
                "Run the classic DOCK 6 anchor-and-grow/grid-scoring route "
                "against explicit precomputed receptor-site and grid artifacts."
            ),
            "when_to_use": (
                "Select as the one primary Step 4 engine when an independently "
                "prepared and audited DOCK 6 sphere/grid bundle is available."
            ),
            "depends_on": [],
            "requires": [
                "dock6_receptor_site", "dock6_grid_nrg", "dock6_grid_bmp",
            ],
            # This installation has no reviewed receptor-surface/grid builder.
            # These assets must therefore be injected or preloaded explicitly;
            # treating unproduced artifacts as implicit context would be unsafe.
            "requires_supplied": [
                "dock6_receptor_site", "dock6_grid_nrg", "dock6_grid_bmp",
            ],
            "runtime": {
                "kind": "command",
                "environment": "dock6",
                "command": str(_DOCK6_BINARY),
                "source_dir": str(_DOCK6_SOURCE_DIR),
                "probe_code": (
                    "from pathlib import Path; import subprocess; "
                    f"root=Path({str(_DOCK6_SOURCE_DIR)!r}); "
                    "required=[root/'bin/dock6',root/'bin/antechamber',"
                    "root/'parameters/vdw_AMBER_parm99.defn',"
                    "root/'parameters/flex.defn',root/'parameters/flex_drive.tbl']; "
                    "assert all(path.is_file() and path.stat().st_size>0 for path in required); "
                    "assert subprocess.run([str(required[0]),'--version'],"
                    "capture_output=True,timeout=20).returncode==0"
                ),
                "probe_timeout_seconds": 30,
                "probe_failure_hint": (
                    "DOCK 6 executable, ligand preparation binary, or required "
                    "parameter tables are unavailable."
                ),
            },
            "params": {
                "dock_root": _p(str(_DOCK6_SOURCE_DIR), "str"),
                "antechamber_path": _p(
                    str(_DOCK6_SOURCE_DIR / "bin" / "antechamber"), "str",
                ),
                "charge_method": _p("gas", "str", choices=["gas", "bcc"]),
                "num_modes": _p(20, "int", min=1, max=100),
                "max_orientations": _p(500, "int", min=1, max=100000),
                "min_anchor_size": _p(5, "int", min=1, max=100),
                "pruning_max_orients": _p(100, "int", min=1, max=10000),
                "score_threshold": _p(
                    1000.0, "float", min=-10000.0, max=1000000.0,
                ),
                "seed": _p(42, "int"),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
                "ligand_prep_timeout_seconds": _p(
                    600, "int", min=30, max=7200,
                ),
            },
        },
        "autodock4": {
            "default_enabled": False,
            "purpose": (
                "Run the classic AutoDock4 Lamarckian genetic-algorithm "
                "workflow with native AutoDockTools preparation and AutoGrid4 "
                "maps."
            ),
            "when_to_use": (
                "Select as the one primary Step 4 engine for a reproducible "
                "AutoDock4 baseline. AutoGrid4 and AutoDockTools are helper "
                "components of this route, not independent scoring engines."
            ),
            "depends_on": [],
            "requires": [],
            "runtime": {
                "kind": "composite",
                "components": {
                    "io": {
                        "kind": "command",
                        "environment": "autodock4",
                        "command": "python",
                        "probe_code": (
                            "import shutil; import AutoDockTools; "
                            "required=('prepare_receptor4','prepare_ligand4',"
                            "'prepare_gpf4','prepare_dpf42'); "
                            "assert all(shutil.which(name) for name in required)"
                        ),
                        "probe_failure_hint": (
                            "AutoDockTools or its AutoDock4 preparation scripts "
                            "are unavailable."
                        ),
                    },
                    "grid": {
                        "kind": "command",
                        "environment": "autodock4",
                        "command": "autogrid4",
                        "probe_args": ["--version"],
                    },
                    "engine": {
                        "kind": "command",
                        "environment": "autodock4",
                        "command": "autodock4",
                        "probe_args": ["--version"],
                    },
                },
            },
            "params": {
                "autobox_add": _p(4.0, "float", min=0.0, max=30.0),
                "spacing": _p(0.375, "float", min=0.01, max=2.0),
                "num_runs": _p(20, "int", min=1, max=256),
                "num_modes": _p(20, "int", min=1, max=100),
                "max_evals": _p(2500000, "int", min=1, max=100000000),
                "population_size": _p(150, "int", min=2, max=10000),
                "cluster_rmsd": _p(2.0, "float", min=0.0, max=20.0),
                "seed": _p(42, "int", min=0),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "rdock": {
            "default_enabled": False,
            "purpose": (
                "Run rDock cavity generation and native docking, preserving "
                "its total and component scores without relabeling their units."
            ),
            "when_to_use": (
                "Select as the one primary Step 4 engine for an independent "
                "classic scoring/search route with a reference-ligand cavity "
                "or an explicit sphere-mapped search region."
            ),
            "depends_on": [],
            "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "rdock",
                "command": "rbdock",
                "commands": {
                    "engine": "rbdock",
                    "cavity": "rbcavity",
                    "converter": "obabel",
                },
                "probe_failure_hint": (
                    "rDock, rbcavity, or its Open Babel converter is unavailable."
                ),
            },
            "params": {
                "cavity_radius": _p(6.0, "float", min=0.1, max=50.0),
                "num_runs": _p(25, "int", min=1, max=1000),
                "num_modes": _p(20, "int", min=1, max=100),
                "receptor_flex": _p(3.0, "float", min=0.0, max=20.0),
                "seed": _p(42, "int", min=0),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "adfr": {
            "default_enabled": False,
            "purpose": (
                "Run the ADFRsuite non-covalent workflow with native ligand/"
                "receptor preparation, AGFR target maps, and ADFR search."
            ),
            "when_to_use": (
                "Select as the one primary Step 4 engine when an ADFR/AGFR "
                "maps-based route is requested. Flexible-residue and covalent "
                "contracts remain unsupported until explicitly implemented."
            ),
            "depends_on": [],
            "requires": [],
            "runtime": {
                "kind": "composite",
                "engine_version": "1.1",
                "components": {
                    "io": {
                        "kind": "command",
                        "environment": "adfrsuite",
                        "command": "prepare_receptor",
                        "probe_args": ["-h"],
                    },
                    "ligand_io": {
                        "kind": "command",
                        "environment": "adfrsuite",
                        "command": "prepare_ligand",
                        "probe_args": ["-h"],
                    },
                    "grid": {
                        "kind": "command",
                        "environment": "adfrsuite",
                        "command": "agfr",
                        "probe_args": ["-h"],
                    },
                    "engine": {
                        "kind": "command",
                        "environment": "adfrsuite",
                        "command": "adfr",
                        "probe_args": ["-h"],
                    },
                },
            },
            "params": {
                "autobox_add": _p(4.0, "float", min=0.0, max=30.0),
                "spacing": _p(0.375, "float", min=0.01, max=2.0),
                "num_runs": _p(20, "int", min=1, max=256),
                "num_modes": _p(20, "int", min=1, max=100),
                "max_evals": _p(2500000, "int", min=1, max=100000000),
                "population_size": _p(150, "int", min=2, max=10000),
                "max_cores": _p(1, "int", min=1, max=256),
                "seed": _p(42, "int", min=0),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "redocking_validation": {
            "default_enabled": True,
            "purpose": "Redock native ligand and compute symmetry-aware RMSD gate.",
            "when_to_use": (
                "Default GNINA protocol validation; disable for every other "
                "primary engine until engine-neutral redocking is implemented."
            ),
            "depends_on": ["gnina"], "requires": ["reference_ligand"],
            "params": {
                "redocking_rmsd_threshold": _p(
                    2.0, "float", min=0.0, max=10.0
                ),
                "redocking_rmsd_hard_threshold": _p(
                    2.5, "float", min=0.0, max=10.0
                ),
                "timeout_seconds": _p(
                    7200, "int", min=60, max=86400
                ),
            },
        },
    },
}

STEPS["pose_evaluation"] = {
    "order": 5,
    "skill": "structure-evaluator",
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    "consumes": ["prepared_receptor", "pocket_definition", "docked_poses", "docking_scores",
                 "docking_report",
                 "reference_ligand", "interaction_blueprint", "pocket_features",
                 "library_manifest"],
    "requires_all": ["prepared_receptor", "docked_poses", "docking_report"],
    "produces": [
        "selected_poses",
        "selection_manifest",
        "evaluation_report",
        "candidate_topologies",
        "rtmscore_manifest",
        "short_md_manifest",
        "trajectory_stability_report",
        "mmgbsa_manifest",
        "openfe_rbfe_manifest",
        "ties_rbfe_manifest",
    ],
    "tool_one_of": {
        "blueprint_match": [["plip"]],
    },
    "params": {
        "screening_mode": _screening_mode_param(),
        "selection_strategy": _p(
            "ranked_per_compound", "str",
            choices=["ranked_per_compound", "first_passing_per_compound", "all", "top_by_cnnscore", "best_per_compound"],
            desc=(
                "Evaluate each compound's docking-ranked poses up to the cohort-size "
                "quota (20 at <=50 compounds, 5 at <=1000, otherwise 1)."
            ),
        ),
        "poses_per_compound": _p(
            20, "int", min=1, max=100,
            desc=(
                "Per-compound cap for ranked_per_compound and shortlist size for "
                "best_per_compound. Ignored by first_passing_per_compound."
            ),
        ),
        "batch_size": _p(50, "int", min=1, max=5000),
        "profiler_workers": _p(2, "int", min=1, max=4),
        "checkpoint_interval_seconds": _p(300, "int", min=0, max=86400),
        "max_poses": _p(
            None, "int", min=1, max=100000,
            desc=(
                "Evaluation budget. With ranked_per_compound or legacy "
                "first_passing_per_compound, caps actual evaluations without "
                "discarding selected poses; exhaustion writes a partial report. "
                "Other strategies cap the selected input poses."
            ),
        ),
        "overall_timeout_seconds": _p(86400, "int", min=60, max=86400),
    },
    "tools": {
        "identity_qc": {
            "default_enabled": True,
            "purpose": (
                "Bind chemical identity, pose lineage, coordinate revision, and "
                "content hashes before any pose evidence is merged."
            ),
            "when_to_use": (
                "Default integrity gate for every evaluated pose and resumable batch."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module",
                "environment": "drugmate",
                "module": "rdkit",
            },
            "params": {},
        },
        "posebusters": {
            "default_enabled": True,
            "purpose": "Geometric validity (bonds, angles, clashes, flatness).",
            "when_to_use": "Default geometric sanity gate.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module", "environment": "drugmate", "module": "posebusters",
            },
            "params": {},
        },
        "posecheck": {
            "default_enabled": True,
            "purpose": "Strain energy + heavy-atom steric clash count.",
            "when_to_use": (
                "Default absolute energetic/clash gate. Optional native-ligand "
                "delta calibration must be enabled explicitly."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module", "environment": "drugmate", "module": "posecheck",
            },
            "params": {
                "num_conformers": _p(
                    10, "int", min=1, max=1000,
                    desc=(
                        "Reference conformers sampled for PoseCheck strain; "
                        "also used for native-ligand calibration."
                    ),
                ),
                "strain_warn_kcal": _p(8.0, "float", min=0.0, max=1000.0),
                "strain_fail_kcal": _p(60.0, "float", min=0.0, max=10000.0),
                "clash_warn_count": _p(0, "int", min=0, max=1000),
                "clash_fail_count": _p(2, "int", min=0, max=1000),
                "calibrate_against_reference": _p(
                    False,
                    "bool",
                    desc=(
                        "Add a validated native-ligand PoseCheck baseline to "
                        "the absolute 8/60 kcal/mol and 0/2 clash thresholds."
                    ),
                ),
            },
        },
        "plip": {
            "default_enabled": True,
            "purpose": "Per-pose protein-ligand interaction detection.",
            "when_to_use": "Needed on its own and as input to blueprint_match.",
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "command",
                "environment": "drugmate",
                "command": "python",
                "probe_code": (
                    "import plip; "
                    "from plip.structure.preparation import PDBComplex"
                ),
                "probe_failure_hint": "PLIP import/API probe failed in the drugmate environment",
            },
            "params": {},
        },
        "prolif": {
            "default_enabled": False,
            "purpose": (
                "Record-only per-pose RDKit-native interaction fingerprints "
                "as a cross-check for PLIP evidence."
            ),
            "when_to_use": (
                "Enable for secondary review or MD/multi-pose fingerprints; "
                "detector disagreement is reported for audit and does not gate."
            ),
            "depends_on": [], "requires": [],
            "runtime": {
                "kind": "python_module",
                "environment": "drugmate",
                "module": "prolif",
                "required_for_preflight": False,
            },
            "params": {
                "vicinity_cutoff": _p(6.0, "float", min=2.0, max=12.0),
                "interaction_types": _p(
                    None, "list", item="str",
                    desc="Optional ProLIF interaction classes; null uses its default set.",
                ),
            },
        },
        "blueprint_match": {
            "default_enabled": True,
            "purpose": "Compare pose interactions to the required blueprint.",
            "when_to_use": (
                "When enforcing required interactions (e.g. hinge H-bond); "
                "requires PLIP; ProLIF is an audit-only cross-check."
            ),
            "depends_on": [], "requires": [],
            "runtime": {"kind": "internal"},
            "params": {},
        },
        "pocket_occupancy": {
            "default_enabled": True,
            "purpose": "Sub-pocket contact/occupancy fractions.",
            "when_to_use": "To diagnose under-filled sub-pockets.",
            "depends_on": [], "requires": [],
            "runtime": {"kind": "internal"},
            "params": {},
        },
        "rdkit_pose": {
            "default_enabled": False,
            "purpose": "Independent MMFF94s strain, TFD and reference-conformer RMSD in INTRA.",
            "when_to_use": "CPU geometry review using a reproducible sampled conformer minimum.",
            "depends_on": [], "requires": [],
            "runtime": {"kind": "python_module", "environment": "drugmate", "module": "rdkit"},
            "params": {
                "num_conformers": _p(8, "int", min=1, max=64),
                "seed": _p(42, "int", min=0, max=2147483647),
                "max_iterations": _p(500, "int", min=1, max=10000),
            },
        },
        "freesasa": {
            "default_enabled": False,
            "purpose": "Ligand-side buried SASA and surface fractions in POCKET.",
            "when_to_use": "CPU ligand burial assessment from the actual docked coordinates.",
            "depends_on": [], "requires": [],
            "runtime": {"kind": "python_module", "environment": "drugmate", "module": "rdkit.Chem.rdFreeSASA"},
            "params": {"probe_radius": _p(1.4, "float", min=0.1, max=3.0)},
        },
        "rtmscore": {
            "default_enabled": False,
            "purpose": (
                "Score post-QC protein-ligand poses with the three-model "
                "RTMScore graph-transformer ensemble while preserving canonical "
                "pose identity."
            ),
            "when_to_use": (
                "Optional machine-learning rescoring after geometric and "
                "interaction QC. Scores are an independent method and are never "
                "arithmetically mixed with docking energies."
            ),
            "depends_on": ["identity_qc", "posebusters", "posecheck"],
            "requires": [],
            "produces": ["rtmscore_manifest"],
            "runtime": {
                "kind": "command",
                "environment": "rtmscore",
                "command": "python",
                "source_dir": str(_RTMSCORE_SOURCE_DIR),
                "entrypoint": str(_RTMSCORE_ENTRYPOINT),
                "checkpoints": [str(path) for path in _RTMSCORE_CHECKPOINTS],
                "asset_sha256": [
                    "399727a85de96c011a75a4c242ebe947753724f01e8d15139fae7bcad106bfc2",
                    "931bdff69f43982d2e10f8914f28c25b35a4d4b90367b10147e456b179da4d29",
                    "ab33a7e0ca2b68081669f3d3a2877846d02d7bc9f7755d31b0e8e0fac2ba8075",
                ],
                "probe_code": (
                    "from pathlib import Path; import hashlib,torch,dgl,MDAnalysis; "
                    "from rdkit import Chem; import prody; "
                    f"entry=Path({str(_RTMSCORE_ENTRYPOINT)!r}); assert entry.is_file(); "
                    f"paths={[str(path) for path in _RTMSCORE_CHECKPOINTS]!r}; "
                    "expected=['399727a85de96c011a75a4c242ebe947753724f01e8d15139fae7bcad106bfc2',"
                    "'931bdff69f43982d2e10f8914f28c25b35a4d4b90367b10147e456b179da4d29',"
                    "'ab33a7e0ca2b68081669f3d3a2877846d02d7bc9f7755d31b0e8e0fac2ba8075']; "
                    "assert all(Path(p).is_file() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==h "
                    "for p,h in zip(paths,expected))"
                ),
                "probe_timeout_seconds": 120,
                "probe_failure_hint": (
                    "RTMScore clean runtime, graph dependencies, entrypoint, or "
                    "one of the pinned ensemble checkpoints failed validation."
                ),
            },
            "params": {
                "model_ensemble": _p(
                    [1, 2, 3], "list", item="int", item_choices=[1, 2, 3],
                ),
                "pocket_cutoff_angstrom": _p(10.0, "float", min=3.0, max=20.0),
                "batch_size": _p(128, "int", min=1, max=4096),
                "num_workers": _p(
                    0, "int", min=0, max=64,
                    desc=(
                        "Graph-loader worker processes. Zero is the portable "
                        "default; increase only on a runtime that permits forked workers."
                    ),
                ),
                "parallel_graph_build": _p(False, "bool"),
                "device": _p("auto", "str", choices=["auto", "cuda", "cpu"]),
                "timeout_seconds": _p(86400, "int", min=60, max=86400),
            },
        },
        "rescore_rtmscore": {
            "default_enabled": False,
            "purpose": (
                "Join RTMScore ensemble outputs to the evaluated cohort by exact "
                "canonical pose_id and optionally rerank survivors by that one signal."
            ),
            "when_to_use": (
                "Enable together with rtmscore when RTMScore should become the "
                "independent post-QC primary ranking metric."
            ),
            "depends_on": ["rtmscore"],
            "requires": [],
            "runtime": {"kind": "internal"},
            "params": {
                "use_as_primary_metric": _p(True, "bool"),
            },
        },
        "candidate_topology": {
            "default_enabled": False,
            "purpose": (
                "Build candidate-specific Amber complex/ligand/receptor topologies "
                "with ff14SB/GAFF2 and AM1-BCC, emitting a per-pose manifest."
            ),
            "when_to_use": (
                "Explicit top-K escalation before MM-GBSA or short MD; run only "
                "on pose-QC survivors. Generic metal and covalent systems are "
                "reported as not_evaluable."
            ),
            "depends_on": ["identity_qc", "posebusters", "posecheck"],
            "requires": [],
            "produces": ["candidate_topologies"],
            "runtime": {
                "kind": "command",
                "environment": "amber",
                "command": "antechamber",
                "commands": {
                    "pdb4amber": "pdb4amber",
                    "antechamber": "antechamber",
                    "parmchk2": "parmchk2",
                    "tleap": "tleap",
                    "sqm": "sqm",
                },
                "probe_failure_hint": (
                    "Candidate topology preparation requires functional "
                    "pdb4amber, antechamber, parmchk2 and tleap commands in "
                    "the amber environment."
                ),
            },
            "params": {
                "execution_mode": _p(
                    "build", "str", choices=["preflight", "build"],
                    desc=(
                        "preflight writes applicability/identity evidence only; "
                        "build additionally invokes AmberTools."
                    ),
                ),
                "protein_force_field": _p(
                    "ff14SB", "str", choices=["ff14SB", "ff19SB"],
                ),
                "ligand_force_field": _p(
                    "gaff2", "str", choices=["gaff2"],
                ),
                "ligand_charge_method": _p(
                    "bcc", "str", choices=["bcc"],
                    desc="AM1-BCC charge assignment through antechamber.",
                ),
                "max_candidates": _p(20, "int", min=1, max=1000),
                "timeout_per_candidate_seconds": _p(
                    2400, "int", min=60, max=86400,
                ),
                "fail_fast": _p(True, "bool"),
                "allow_receptor_hetero": _p(
                    False, "bool",
                    desc=(
                        "Expert override for pre-parameterized receptor HETATM "
                        "records; false keeps the generic route conservative."
                    ),
                ),
            },
        },
        "openmm_md": {
            "default_enabled": False,
            "purpose": (
                "Run short explicit-solvent periodic OpenMM MD for a bounded "
                "post-QC candidate subset and retain trajectories/state data."
            ),
            "when_to_use": (
                "Optional pose-stability escalation after candidate topology "
                "construction. This is a short diagnostic trajectory, not a "
                "binding free-energy calculation."
            ),
            "depends_on": [
                "identity_qc", "posebusters", "posecheck", "candidate_topology",
            ],
            "requires": [],
            "produces": ["short_md_manifest"],
            "runtime": {
                "kind": "command",
                "environment": "openmm",
                "command": "python",
                "probe_code": (
                    "import openmm,openmm.app as app; from openmm import unit; "
                    "system=openmm.System(); system.addParticle(1.0); "
                    "integrator=openmm.VerletIntegrator(0.001); "
                    "platform=openmm.Platform.getPlatformByName('CUDA'); "
                    "context=openmm.Context(system,integrator,platform); "
                    "context.setPositions([[0,0,0]]); context.getState(getEnergy=True); "
                    "del context,integrator"
                ),
                "probe_timeout_seconds": 90,
                "probe_failure_hint": (
                    "The dedicated openmm environment could not create a real "
                    "CUDA Context; do not launch short MD on this host."
                ),
            },
            "params": {
                "top_k": _p(3, "int", min=1, max=100),
                "platform": _p(
                    "CUDA", "str", choices=["CUDA"],
                    desc=(
                        "The registered runtime is preflighted with a real CUDA "
                        "Context; other platforms require a separately validated route."
                    ),
                ),
                "precision": _p(
                    "mixed", "str", choices=["mixed", "single", "double"],
                ),
                "device_index": _p(0, "int", min=0, max=16),
                "temperature_kelvin": _p(300.0, "float", min=250.0, max=400.0),
                "pressure_bar": _p(1.0, "float", min=0.1, max=100.0),
                "timestep_femtoseconds": _p(
                    2.0, "float", min=0.5, max=2.0,
                    desc=(
                        "No hydrogen-mass repartitioning is used; the audited "
                        "HBonds-constrained route is capped at 2 fs."
                    ),
                ),
                "friction_per_picosecond": _p(1.0, "float", min=0.01, max=100.0),
                "minimization_max_iterations": _p(500, "int", min=0, max=100000),
                "equilibration_picoseconds": _p(20.0, "float", min=0.0, max=10000.0),
                "production_picoseconds": _p(100.0, "float", min=1.0, max=100000.0),
                "report_interval_picoseconds": _p(2.0, "float", min=0.01, max=1000.0),
                "solvent_padding_angstrom": _p(10.0, "float", min=6.0, max=30.0),
                "water_model": _p("tip3p", "str", choices=["tip3p"]),
                "seed": _p(42, "int"),
                "timeout_per_candidate_seconds": _p(
                    86400, "int", min=60, max=604800,
                ),
                "fail_fast": _p(False, "bool"),
            },
        },
        "mdanalysis": {
            "default_enabled": False,
            "purpose": (
                "Analyze OpenMM trajectories after protein alignment and report "
                "ligand RMSD, center displacement, and native-contact retention."
            ),
            "when_to_use": (
                "Enable with openmm_md when trajectory-stability screening "
                "evidence is required. Shared sbdd numeric bounds are heuristic "
                "gates, not universal experimental cutoffs."
            ),
            "depends_on": ["openmm_md"],
            "requires": [],
            "produces": ["trajectory_stability_report"],
            "runtime": {
                "kind": "command",
                "environment": "openmm",
                "command": "python",
                "probe_code": (
                    "import MDAnalysis,networkx; from MDAnalysis.analysis import align,rms; "
                    "from MDAnalysis.lib.distances import distance_array; "
                    "from MDAnalysis.transformations.nojump import NoJump"
                ),
                "probe_failure_hint": (
                    "MDAnalysis trajectory/alignment APIs are unavailable in the "
                    "dedicated openmm environment."
                ),
            },
            "params": {
                "alignment_selection": _p("protein and backbone", "str"),
                "ligand_selection": _p("resname LIG and not name H*", "str"),
                "contact_cutoff_angstrom": _p(4.5, "float", min=2.0, max=10.0),
                "analysis_stride": _p(1, "int", min=1, max=10000),
                "timeout_per_candidate_seconds": _p(
                    3600, "int", min=30, max=86400,
                ),
            },
        },
        "mmgbsa": {
            "default_enabled": False,
            "purpose": (
                "Run AmberTools MMPBSA.py on post-QC candidate topology "
                "bundles and emit a pose-bound, single-frame endpoint "
                "reranking manifest."
            ),
            "when_to_use": (
                "Explicit top-K MM-GBSA escalation after candidate topology "
                "build. Disabled for routine pose QC; current scores are "
                "single-frame, no-entropy estimates rather than absolute "
                "binding free energies."
            ),
            "depends_on": [
                "identity_qc",
                "posebusters",
                "posecheck",
                "candidate_topology",
            ],
            "requires": [],
            "produces": ["mmgbsa_manifest"],
            "runtime": {
                "kind": "command",
                "environment": "mmgbsa",
                "command": "MMPBSA.py",
                "commands": {
                    "python": "python",
                    "mmpbsa": "MMPBSA.py",
                    "cpptraj": "cpptraj",
                    "energy": "mmpbsa_py_energy",
                },
                "probe_code": (
                    "import shutil; import MMPBSA_mods, parmed; "
                    "assert all(shutil.which(name) for name in "
                    "('MMPBSA.py','cpptraj','mmpbsa_py_energy'))"
                ),
                "probe_failure_hint": (
                    "Amber MM-GBSA requires importable MMPBSA_mods/ParmEd and "
                    "functional MMPBSA.py, cpptraj and mmpbsa_py_energy "
                    "commands in the registered gmxMMPBSA environment."
                ),
            },
            "params": {
                "execution_mode": _p(
                    "preflight",
                    "str",
                    choices=["preflight", "score"],
                    desc=(
                        "preflight validates the runtime and built topology "
                        "bundle; score additionally runs single-frame "
                        "MMPBSA.py for each eligible top-K pose."
                    ),
                ),
                "top_k": _p(20, "int", min=1, max=1000),
                "gb_model": _p(
                    5,
                    "int",
                    choices=[1, 2, 5, 7, 8],
                    desc="Amber igb model used consistently across the cohort.",
                ),
                "salt_concentration_molar": _p(
                    0.15,
                    "float",
                    min=0.0,
                    max=5.0,
                ),
                "timeout_per_candidate_seconds": _p(
                    900,
                    "int",
                    min=30,
                    max=86400,
                ),
                "fail_fast": _p(False, "bool"),
                "keep_intermediate_files": _p(False, "bool"),
                "use_as_primary_metric": _p(
                    True,
                    "bool",
                    desc=(
                        "When score mode succeeds, rerank the post-QC top "
                        "candidate summary by MM-GBSA only; never average it "
                        "with docking or RBFE metrics."
                    ),
                ),
            },
        },
        "openfe_rbfe": {
            "default_enabled": False,
            "purpose": (
                "Validate post-QC ligand identities/chemistry and optionally "
                "build an OpenFE ligand-mapping network for later RBFE work."
            ),
            "when_to_use": (
                "Contextual RBFE escalation. Generic FEP/RBFE requests select "
                "this backend, but routine pose_evaluation leaves it disabled."
            ),
            "depends_on": ["identity_qc", "posebusters", "posecheck"],
            "requires": [],
            "produces": ["openfe_rbfe_manifest"],
            "runtime": {
                "kind": "command",
                "environment": "openfe",
                "command": "python",
                "commands": {"python": "python", "openfe": "openfe"},
                "probe_timeout_seconds": 60,
                "probe_code": (
                    "import shutil; import openfe, gufe, openmm, "
                    "openmmtools, rdkit; import openff.toolkit; "
                    "assert shutil.which('openfe')"
                ),
                "probe_failure_hint": (
                    "OpenFE planning requires functional OpenFE/GUFE/OpenMM/"
                    "OpenFF/OpenMMTools/RDKit imports and the openfe CLI in the "
                    "registered openfe environment."
                ),
            },
            "params": {
                "execution_mode": _p(
                    "preflight",
                    "str",
                    choices=["preflight", "plan"],
                    desc=(
                        "preflight checks runtime, identity and applicability; "
                        "plan additionally creates a mapping-only ligand network. "
                        "Neither mode runs quickrun or reports ΔΔG."
                    ),
                ),
                "network_method": _p(
                    "minimal_spanning",
                    "str",
                    choices=[
                        "minimal_spanning",
                        "minimal_redundant",
                        "radial",
                    ],
                ),
                "mapper": _p(
                    "KartografAtomMapper",
                    "str",
                    choices=["KartografAtomMapper", "LomapAtomMapper"],
                ),
                "mapper_settings": _p(None, "dict"),
                "max_candidates": _p(8, "int", min=2, max=200),
                "min_ligands": _p(2, "int", min=2, max=200),
                "allow_charge_changes": _p(
                    False,
                    "bool",
                    desc=(
                        "Charge-changing mappings are excluded unless the user "
                        "explicitly accepts the need for a reviewed correction protocol."
                    ),
                ),
                "timeout_seconds": _p(600, "int", min=1, max=86400),
                "n_processes": _p(1, "int", min=1, max=32),
                "mst_num": _p(2, "int", min=1, max=10),
                "central_pose_id": _p(None, "str"),
                "min_mapping_score": _p(
                    0.3, "float", min=0.0, max=1.0,
                ),
                "min_mapped_heavy_fraction": _p(
                    0.5, "float", min=0.0, max=1.0,
                ),
                "require_explicit_hydrogens": _p(True, "bool"),
                "allow_undefined_stereochemistry": _p(False, "bool"),
            },
        },
        "ties_rbfe": {
            "default_enabled": False,
            "purpose": (
                "Validate post-QC ligand identities and write an auditable TIES "
                "pair/campaign plan without claiming a computed free energy."
            ),
            "when_to_use": (
                "Only when the user explicitly requests TIES RBFE planning. The "
                "generic RBFE default is OpenFE; TIES is never a silent fallback."
            ),
            "depends_on": ["identity_qc", "posebusters", "posecheck"],
            "requires": [],
            "produces": ["ties_rbfe_manifest"],
            "runtime": {
                "kind": "command",
                "environment": "ties",
                "command": "python",
                "commands": {
                    "python": "python",
                    "ties": "ties",
                    "ties_md": "ties_md",
                    "ties_ana": "ties_ana",
                    "antechamber": "antechamber",
                    "parmchk2": "parmchk2",
                    "tleap": "tleap",
                    "sqm": "sqm",
                },
                "probe_code": (
                    "import importlib, shutil; import ties, TIES_MD, openmm; "
                    "assert callable(ties.Pair); "
                    "assert callable(ties.Config); "
                    "assert callable(ties.Pair.superimpose); "
                    "assert callable(importlib.import_module("
                    "'TIES_MD.TIES').TIES); "
                    "importlib.import_module('TIES_MD.ties_analysis'); "
                    "assert all(shutil.which(name) for name in "
                    "('ties','ties_md','ties_ana','antechamber','parmchk2',"
                    "'tleap','sqm'))"
                ),
                "probe_failure_hint": (
                    "TIES planning requires the builder API/CLI, "
                    "TIES_MD runner/analysis APIs and CLIs, its ABI-compatible "
                    "OpenMM, and AmberTools "
                    "antechamber/parmchk2/tleap/sqm in the registered ties "
                    "environment. Do not inject OpenMM from the Python-3.11 "
                    "drugmate environment into the Python-3.10 ties runtime."
                ),
            },
            "params": {
                "execution_mode": _p(
                    "preflight",
                    "str",
                    choices=["preflight", "plan"],
                    desc=(
                        "preflight checks environment/identity/applicability; "
                        "plan additionally writes pending pairwise stages. "
                        "Neither mode runs MD or reports ΔΔG."
                    ),
                ),
                "network_strategy": _p("star", "str", choices=["star"]),
                "reference_pose_id": _p(
                    None,
                    "str",
                    desc=(
                        "Optional canonical post-QC pose_id used as the star "
                        "network reference; null selects the highest-ranked survivor."
                    ),
                ),
                "max_candidates": _p(8, "int", min=2, max=100),
                "require_post_qc_lineage": _p(True, "bool"),
                "runtime_probe_timeout_seconds": _p(
                    30, "int", min=1, max=600,
                ),
                "simulation_platform": _p(
                    "auto",
                    "str",
                    choices=["auto", "CUDA", "OpenCL", "CPU"],
                    desc=(
                        "Requested future TIES_MD OpenMM platform. auto probes "
                        "CUDA→OpenCL→CPU and records the first "
                        "functional option; planning never launches MD."
                    ),
                ),
            },
        },
    },
}

# Reuse the same isolated profiler APIs on each candidate pose. These switches
# are independent of reference-complex analysis in pocket_analysis.
for _candidate_profiler in ("arpeggio", "luna", "oddt"):
    _candidate_spec = deepcopy(STEPS["pocket_analysis"]["tools"][_candidate_profiler])
    _candidate_spec["requires"] = []
    _candidate_spec["purpose"] = "Candidate pose " + _candidate_profiler + " interactions/fingerprints."
    _candidate_spec["when_to_use"] = "Optional per-pose assessment; Arpeggio numeric QC, LUNA/ODDT audit fingerprints."
    _candidate_spec["params"]["timeout_seconds"] = _p(300, "int", min=10, max=7200)
    STEPS["pose_evaluation"]["tools"][_candidate_profiler] = _candidate_spec


STEPS["property_filter"] = {
    "order": 6,
    "skill": "drug-property-filter",
    "runtime": {"environment": "drugmate"},
    "default_enabled": True,
    "consumes": ["selected_poses", "docked_poses"],
    "requires_any": [["selected_poses", "docked_poses"]],
    "produces": ["property_report", "synthesis_report"],
    "params": {
        "screening_mode": _screening_mode_param(),
    },
    "tools": _property_tools(
        admet_default=True,
        chemfh_default=True,
        dedup_default=True,
        aizynth_default=True,
    ),
}

STEPS["coordination"] = {
    "order": 7,
    "skill": "reflective-coordinator",
    "default_enabled": True,
    "consumes": ["evaluation_report", "selection_manifest",
                 "property_report", "docking_scores", "docking_report",
                 "pocket_analysis_report"],
    "produces": ["coordinator_report"],
    "params": {
        "pass_rate_threshold": _p(0.3, "float", min=0.0, max=1.0, desc="Legacy input, ignored by round assessment; the Agent owns campaign goals."),
        "min_candidates": _p(3, "int", min=1, max=1000, desc="Legacy input, ignored; nomination eligibility is assessed per candidate."),
        "max_rounds": _p(3, "int", min=1, max=20, desc="Legacy input, ignored; campaign max_runs bounds execution attempts."),
        "min_validation_actives": _p(20, "int", min=1, max=100000),
        "min_validation_decoys": _p(50, "int", min=1, max=1000000),
        "min_validation_auroc": _p(0.70, "float", min=0.0, max=1.0),
        "min_validation_auroc_hard": _p(
            0.50, "float", min=0.0, max=1.0,
        ),
        "min_validation_ef1": _p(3.0, "float", min=0.0),
        "min_validation_ef1_hard": _p(1.0, "float", min=0.0),
        "normalization_method": _p(
            NORMALIZATION_METHOD,
            "str",
            choices=[NORMALIZATION_METHOD],
            desc=(
                "Use descending RANK.AVG percentiles for metrics and category "
                "scores within one comparable round/profile cohort."
            ),
        ),
        "metric_weight_overrides": _p(
            None,
            "dict",
            desc=(
                "Session-local metric_id -> 0..5 magnitude overrides. The "
                "catalog direction/sign and all hard gates are preserved."
            ),
        ),
        "category_weight_overrides": _p(
            None,
            "dict",
            desc=(
                "Session-local category code -> 0..5 parent-weight overrides; "
                "hard gates are unchanged."
            ),
        ),
        "screening_mode": _screening_mode_param(),
    },
    "tools": {
        "require_quality_gates": {
            "default_enabled": True,
            "purpose": (
                "Block nomination on protocol gates; modify requests a rerun "
                "and reject invalidates the protocol."
            ),
            "when_to_use": "Default; keep on for trustworthy nominations.",
            "depends_on": [], "requires": [],
            "params": {},
        },
        "screening_benchmark": {
            "default_enabled": False,
            "purpose": "Compute AUROC/EF from ACT_/DEC_ labels when present.",
            "when_to_use": "Only with verified ACT_ and DEC_ labels in the actual screening cohort. Numeric or unknown IDs alone provide no activity labels; disable for unlabeled screening.",
            "depends_on": [], "requires": [],
            "params": {},
        },
    },
}

STEPS["reporting"] = {
    "order": 8,
    "skill": "report-generator",
    "default_enabled": True,
    "best_effort": True,
    "consumes": ["coordinator_report"],
    "produces": [
        "report_html",
        "report_data",
        "report_assets_manifest",
        "report_artifact_manifest",
        "report_pymol_cover",
    ],
    "params": {
        "top_n": _p(10, "int", min=1, max=1000),
        "title": _p(None, "str"),
        "assets_source_dir": _p(None, "str"),
    },
    "tools": {
        "rdkit_2d": {
            "default_enabled": True,
            "purpose": "RDKit 2D structure depictions on candidate cards.",
            "when_to_use": "Default.",
            "depends_on": [], "requires": [], "params": {},
        },
        "mol_3d": {
            "default_enabled": True,
            "purpose": "3Dmol.js per-pose viewer cards.",
            "when_to_use": "Default.",
            "depends_on": [], "requires": [], "params": {},
        },
        "molstar_main": {
            "default_enabled": True,
            "purpose": "Mol* main viewer (pocket + reference ligand).",
            "when_to_use": "Default top viewer.",
            "depends_on": [], "requires": [], "params": {},
        },
        "plots": {
            "default_enabled": True,
            "purpose": "matplotlib QC / druglikeness / ROC-EF panels.",
            "when_to_use": "Default.",
            "depends_on": [], "requires": [], "params": {},
        },
        "pymol_cover": {
            "default_enabled": False,
            "purpose": "Static PyMOL cover PNG.",
            "when_to_use": "Optional; off by default.",
            "depends_on": [], "requires": [],
            "produces": ["report_pymol_cover"],
            "params": {},
        },
    },
}


# --- Convenience accessors ----------------------------------------------------
def steps_in_order() -> list[str]:
    return [name for name, _ in sorted(STEPS.items(), key=lambda kv: kv[1]["order"])]


def step_by_order() -> dict[float, str]:
    return {spec["order"]: name for name, spec in STEPS.items()}


def get_step(step_id: str) -> dict:
    if step_id not in STEPS:
        raise KeyError(f"unknown step: {step_id!r}")
    return STEPS[step_id]


def get_tool(step_id: str, tool_id: str) -> dict:
    tools = get_step(step_id)["tools"]
    if tool_id not in tools:
        raise KeyError(f"unknown tool {tool_id!r} in step {step_id!r}")
    return tools[tool_id]

"""Three-tier gates and hierarchical rank-percentile scoring.

The rules in this module are the auditable, executable subset of the metric
tables currently produced by docking, ``pose_evaluation``, and the property
filter.  A candidate rule is included only when the pipeline emits a stable
molecule/pose-level scalar for it.  Protocol rules are kept separate so a
failed calibration can roll back a workflow without being mislabeled as a
failed molecule.

Two-level gate contract
-------------------------------------

Level 1 -- one bounded metric, one molecule.  A bound comparison produces a
*label* and nothing else:

* below ``hard_lower`` or above ``hard_upper`` -> ``fail``
* outside the soft interval but inside the hard interval -> ``flag``
* inside the soft interval (including the boundary) -> ``pass``
* no finite value, or no declared gate bound -> ``not_applicable`` (leaves both
  numerator and denominator)

A single ``fail`` never performs direct removal.  It can nevertheless be the
only observed label in the later category fold; in that case its fail ratio is
100% and the *category fold* rejects in V2/V3.

Level 2 -- one parent category, one molecule.  Once every planned metric
in a category has reached a terminal value-or-NA state, its observed labels are
folded by count and proportion::

    count(fail)
    ratio(fail) = count(fail) / (count(fail) + count(flag) + count(pass))

and the active screening mode turns those two numbers into the molecule
verdict (``V1``/``V2``/``V3``, see ``SCREENING_MODES``):

* any category trips its mode's ``count``/``ratio`` threshold -> ``reject``
* otherwise, every evaluated category has more than half ``pass`` labels
  -> ``keep``
* otherwise                                                -> ``modify``

``keep``/``modify``/``reject`` therefore exist only at molecule level and only
after the category fold: ``classification == metric_gate``.  A metric that is
still pending is not a missing value and must not trigger an interim fold of its
unfinished category.

The count and ratio thresholds are two independent alternatives and nothing
else conditions them.  In particular the ratio branch carries
no minimum fail count, so a lone ``fail`` in an otherwise sparsely labeled
category is ratio 100% and does eliminate. No additional observed-coverage
guard modifies the numeric formula.

An omitted lower/upper bound behaves as negative/positive infinity.  A zero
weight metric still carries a label, but never contributes to a category score.
Rejected molecules are removed before every percentile pool.  Missing evidence
is removed from that candidate's numerator and denominator rather than zeroing
a metric for the whole cohort, and a coverage or tool gap must never reject a
molecule or a pose.

For backward compatibility every per-metric assessment still publishes the
legacy ``gate``/``molecule_gate`` spelling of its label (``pass``->``keep``,
``flag``->``modify``, ``fail``->``reject``) alongside the new ``label``/
``molecule_label`` fields.  Those legacy values are *marks*, not verdicts.
"""
from __future__ import annotations

import json
import math
import re
from copy import deepcopy
from typing import Any, Iterable

from pose_scores import collect_pose_scores


METRIC_RULES_VERSION = 12
SCORING_SCHEMA_VERSION = 4
NORMALIZATION_METHOD = "rank_avg_percentile_v1"

# Coverage >=0.50 permits formal ranking; 0.10-0.50 permits provisional
# ranking; below 0.10 the overall score/rank is NA and evidence needs repair.
# Coverage and tool gaps cannot reject molecules or poses. These thresholds
# only choose a rank release band.
FORMAL_COVERAGE_THRESHOLD = 0.50
PROVISIONAL_COVERAGE_THRESHOLD = 0.10

# The pipeline's scored stages, in report order.  Step 4 and Step 7 share one
# small-molecule property stage in the report numbering.
STEP_IDS: tuple[str, ...] = ("step4", "step5", "step6")


# ---------------------------------------------------------------------------
# Level 1: per-metric labels
# ---------------------------------------------------------------------------

# The label a single metric puts on a single molecule.
METRIC_LABELS: tuple[str, ...] = ("pass", "flag", "fail")
NOT_APPLICABLE_LABEL = "not_applicable"

# Severity order used when one metric is described by more than one source.
_LABEL_SEVERITY: dict[str, int] = {"pass": 0, "flag": 1, "fail": 2}

# Legacy spelling of the same three labels.  Kept because published reports,
# category tables and several skills still key on it; it names a *mark*, never a
# molecule verdict.
LABEL_TO_LEGACY_GATE: dict[str, str] = {
    "pass": "keep",
    "flag": "modify",
    "fail": "reject",
    NOT_APPLICABLE_LABEL: NOT_APPLICABLE_LABEL,
}
LEGACY_GATE_TO_LABEL: dict[str, str] = {
    "keep": "pass",
    "modify": "flag",
    "reject": "fail",
    NOT_APPLICABLE_LABEL: NOT_APPLICABLE_LABEL,
}


# ---------------------------------------------------------------------------
# Level 2: category fold -> molecule verdict
# ---------------------------------------------------------------------------

# ``fail_count_threshold`` and ``fail_ratio_threshold`` are both
# inclusive (">=").  ``V1`` performs no elimination at all: labels stay as
# marks, so the worst verdict it can produce is ``modify``.
#
# The two branches are independent alternatives and nothing else conditions
# them. The ratio branch carries no minimum fail count: one fail out of one
# labeled metric is ratio 100% and rejects under both V2 and V3.

# Numeric V2/V3 count OR ratio thresholds are applied without a coverage guard.
# Missing values leave the fail ratio; coverage separately controls rank release.

# ``molecule_effect`` says what a level-1 ``fail`` on this row is *about*.
# ``gate`` is a defect of the molecule.  The rest describe the workflow or are
# excluded from independent molecular scoring, so their ``fail`` cannot enter a level-2
# category count -- it is demoted to ``flag`` and handled by rollback instead.
# ``protocol``/``pose_export`` additionally block publishing a rank, because a
# broken protocol or an unexportable pose makes the *number* untrustworthy;
# ``record``/``route`` do not, since the measurement itself is fine.
MOLECULE_EFFECTS = frozenset({
    "gate", "protocol", "pose_export", "route", "record",
})
NON_MOLECULAR_EFFECTS = frozenset(MOLECULE_EFFECTS - {"gate"})
RANKING_BLOCKING_EFFECTS = frozenset({"protocol", "pose_export"})
SCREENING_MODES: dict[str, dict[str, Any]] = {
    "V1": {
        "mode_id": "V1",
        "name": "marking",
        "eliminates": False,
        "fail_count_threshold": None,
        "fail_ratio_threshold": None,
        "description": (
            "No molecule is eliminated; pass/flag/fail stay per-metric marks."
        ),
    },
    "V2": {
        "mode_id": "V2",
        "name": "standard",
        "eliminates": True,
        "fail_count_threshold": 3,
        "fail_ratio_threshold": 0.50,
        "description": (
            "Default. Reject when any parent category reaches 3 fails or a "
            "50% fail ratio."
        ),
    },
    "V3": {
        "mode_id": "V3",
        "name": "strict",
        "eliminates": True,
        "fail_count_threshold": 2,
        "fail_ratio_threshold": 0.30,
        "description": (
            "Strict. Reject when any parent category reaches 2 fails or a "
            "30% fail ratio."
        ),
    },
}
DEFAULT_SCREENING_MODE = "V2"

# Selection rule from the executable policy: a run that generates its own molecules can
# afford the strict mode, a run that randomly samples a large user library uses
# the standard mode, and a run that consumes a small user-supplied set in full
# only marks.  Callers may always name the mode explicitly instead.
SCREENING_MODE_SELECTION_HINTS: dict[str, str] = {
    "molecule_generation": "V3",
    "sampled_library": "V2",
    "user_supplied_full": "V1",
}



_ADMET_AI_ENDPOINTS_BY_CATEGORY: dict[str, tuple[str, ...]] = {
    "ABS": (
        "HIA_Hou", "Caco2_Wang", "PAMPA_NCATS", "Pgp_Broccatelli",
        "Bioavailability_Ma", "Solubility_AqSolDB",
        "Lipophilicity_AstraZeneca", "HydrationFreeEnergy_FreeSolv",
    ),
    "DIST": ("BBB_Martins", "PPBR_AZ", "VDss_Lombardo"),
    "MET": (
        "CYP2C9_Substrate_CarbonMangels",
        "CYP2D6_Substrate_CarbonMangels",
        "CYP3A4_Substrate_CarbonMangels",
        "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith",
        "CYP2D6_Veith", "CYP3A4_Veith",
    ),
    "EXCR": (
        "Half_Life_Obach", "Clearance_Hepatocyte_AZ",
        "Clearance_Microsome_AZ",
    ),
    "TOX": (
        "hERG", "ClinTox", "AMES", "DILI", "Carcinogens_Lagunin",
        "LD50_Zhu", "Skin_Reaction", "NR-AR", "NR-AR-LBD", "NR-AhR",
        "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
        "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
    ),
}

# ADMET-AI returns a numeric payload even when an endpoint violates its
# scientific domain.  The worker publishes those violations in
# ``invalid_endpoints``, but old reports may predate that field, so the scorer
# also enforces the small set of provider-independent domains here.  Evidence
# validity precedes numeric policy: an invalid value is NA for both gating and
# ranking.  The table's numeric hard/soft bounds remain authoritative for values
# that pass this endpoint-domain check.
_ADMET_AI_NON_PROBABILITY_ENDPOINT_IDS = frozenset({
    "Caco2_Wang",
    "Solubility_AqSolDB",
    "Lipophilicity_AstraZeneca",
    "HydrationFreeEnergy_FreeSolv",
    "PPBR_AZ",
    "VDss_Lombardo",
    "Half_Life_Obach",
    "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ",
    "LD50_Zhu",
})
_ADMET_AI_PROBABILITY_ENDPOINT_IDS = frozenset(
    endpoint_id
    for endpoint_ids in _ADMET_AI_ENDPOINTS_BY_CATEGORY.values()
    for endpoint_id in endpoint_ids
    if endpoint_id not in _ADMET_AI_NON_PROBABILITY_ENDPOINT_IDS
)
_ADMET_AI_NONNEGATIVE_ENDPOINT_IDS = frozenset({
    "VDss_Lombardo",
    "Half_Life_Obach",
    "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ",
})

_ADMET_AI_DERIVED_METRIC_IDS_BY_CATEGORY = {
    "MET": ("admet_ai.cyp_substrate_count", "admet_ai.cyp_inhibition_count"),
}

_CHEMFH_ENDPOINTS: tuple[tuple[str, str], ...] = (
    ("other_assay_interference", "Other assay interference"),
    ("blue_fluorescence", "Blue fluorescence"),
    ("fluc_inhibitors", "FLuc inhibitors"),
    ("promiscuous_compounds", "Promiscuous compounds"),
    ("colloidal_aggregators", "Colloidal aggregators"),
    ("reactive_compounds", "Reactive compounds"),
    ("green_fluorescence", "Green fluorescence"),
)
# The executable policy gives per-endpoint weights inside ALERT; the four unlisted
# endpoints stay at 0 and remain pure gate/audit signals.
_CHEMFH_ENDPOINT_WEIGHTS: dict[str, float] = {
    "promiscuous_compounds": 0.5,
    "colloidal_aggregators": 0.5,
    "reactive_compounds": 1.0,
}
_CHEMFH_SUBSTRUCTURE_RULE_NAMES: tuple[str, ...] = (
    "ALARM_NMR",
    "BMS",
    "Chelator_Rule",
    "GST_FHs_Rule",
    "His_FHs_Rule",
    "Luciferase_Inhibitor_Rule",
    "NTD",
    "PAINS",
    "Potential_Electrophilic_Rule",
    "Lilly",
)


# Category weights and complete-catalog absolute child-weight totals are
# configured below. Catalog coverage is informational; planned/category
# coverage controls whether an overall rank may be released.
CATEGORY_CONFIG: dict[str, dict[str, Any]] = {
    "PHY": {"step_id": "step4", "parent_weight": 3.0, "catalog_abs_weight": 9.0},
    # ALERT has parent weight 1 and seven nonzero children with sum(|w|)=7.
    # It contributes to scoring with low weight and also supplies gate marks.
    "ALERT": {"step_id": "step4", "parent_weight": 1.0, "catalog_abs_weight": 7.0},
    "SYN": {"step_id": "step4", "parent_weight": 2.0, "catalog_abs_weight": 11.0},
    "ABS": {"step_id": "step4", "parent_weight": 4.0, "catalog_abs_weight": 29.0},
    "DIST": {"step_id": "step4", "parent_weight": 3.0, "catalog_abs_weight": 9.0},
    "MET": {"step_id": "step4", "parent_weight": 3.0, "catalog_abs_weight": 20.2},
    "EXCR": {"step_id": "step4", "parent_weight": 2.0, "catalog_abs_weight": 3.0},
    "TOX": {"step_id": "step4", "parent_weight": 5.0, "catalog_abs_weight": 77.5},
    "AFF": {"step_id": "step5", "parent_weight": 5.0, "catalog_abs_weight": 31.0},
    "DRAW": {"step_id": "step5", "parent_weight": 3.0, "catalog_abs_weight": 13.0},
    "FLEX": {"step_id": "step5", "parent_weight": 3.0, "catalog_abs_weight": 11.0},
    "BENCH": {"step_id": "step5", "parent_weight": 0.0, "catalog_abs_weight": 0.0},
    "IDQC": {"step_id": "step6", "parent_weight": 0.0, "catalog_abs_weight": 0.0},
    "TOPO": {"step_id": "step6", "parent_weight": 0.0, "catalog_abs_weight": 0.0},
    "PBQ": {"step_id": "step6", "parent_weight": 0.0, "catalog_abs_weight": 3.0},
    "CHEM": {"step_id": "step6", "parent_weight": 0.0, "catalog_abs_weight": 0.0},
    "INTRA": {"step_id": "step6", "parent_weight": 4.0, "catalog_abs_weight": 18.0},
    "INTER": {"step_id": "step6", "parent_weight": 5.0, "catalog_abs_weight": 20.0},
    "DQC": {"step_id": "step6", "parent_weight": 0.0, "catalog_abs_weight": 0.0},
    "FP": {"step_id": "step6", "parent_weight": 0.0, "catalog_abs_weight": 0.0},
    "INTX": {"step_id": "step6", "parent_weight": 5.0, "catalog_abs_weight": 17.0},
    "POCKET": {"step_id": "step6", "parent_weight": 5.0, "catalog_abs_weight": 23.0},
    "DYN": {"step_id": "step6", "parent_weight": 3.0, "catalog_abs_weight": 56.0},
}


# Per-category planned-coverage floors: PHY/ABS 0.60 and TOX 0.70.
# These implementation thresholds ask for evidence repair or hold a formal
# release; coverage can never reject a molecule or pose. Every other
# category relies on the global band alone rather than an invented default.
CATEGORY_PLANNED_COVERAGE_FLOOR: dict[str, float] = {
    "PHY": 0.60,
    "ABS": 0.60,
    "TOX": 0.70,
}


_CATEGORY_METRIC_IDS: dict[str, set[str]] = {
    "AFF": {
        "autodock_vina.vina_affinity", "unidock.docking_score",
        "gnina.minimizedAffinity", "diffdock.confidence", "dock6.grid_score",
        "autodock4.binding_energy", "rdock.score", "adfr.affinity",
    },
    "DRAW": {
        "autodock_vina.ligand_efficiency",
        "autodock_vina.engine_pose_rank", "autodock_vina.score_gap_1_2",
        "autodock_vina.vina_inter", "autodock_vina.vina_intra_best_pose",
        "gnina.intramol", "gnina.CNNscore", "gnina.CNNaffinity",
        "gnina.CNNvariance", "gnina.engine_pose_rank",
    },
    "IDQC": {"posebusters.sanitization", "posebusters.inchi_convertible"},
    "PBQ": {"posebusters.overall_success_ratio"},
    "CHEM": {"posebusters.all_atoms_connected"},
    "INTRA": {
        "posebusters.bond_lengths", "posebusters.bond_angles",
        "posebusters.internal_steric_clash",
        "posebusters.aromatic_ring_flatness",
        "posebusters.non_aromatic_ring_non_flatness",
        "posebusters.double_bond_flatness", "posebusters.internal_energy",
        "posecheck.strain_energy_kcal",
    },
    "INTER": {
        "posebusters.receptor_ligand_maximum_distance",
        "posebusters.minimum_distance_to_receptor",
        "posebusters.minimum_distance_to_organic_cofactors",
        "posebusters.minimum_distance_to_inorganic_cofactors",
        "posebusters.minimum_distance_to_waters",
        "posebusters.volume_overlap_with_receptor",
        "posebusters.volume_overlap_with_organic_cofactors",
        "posebusters.volume_overlap_with_inorganic_cofactors",
        "posebusters.volume_overlap_with_waters",
        "posecheck.steric_clash_count",
    },
    "INTX": {
        "plip.num_interactions", "prolif.num_interactions",
        "blueprint_match.verdict", "blueprint_match.required_hit_rate",
        "blueprint_match.desired_hit_rate",
        "blueprint_match.desired_ranking_score",
        "blueprint_match.desired_evaluable_coverage",
        "blueprint_match.optional_hits",
    },
    "POCKET": {
        "pocket_occupancy.hydrophobic_region", "pocket_occupancy.polar_region",
        "pocket_occupancy.charged_region", "pocket_occupancy.aromatic_region",
        "pocket_occupancy.overall",
    },
    "DYN": {"mmgbsa.delta_g_bind"},
    "PHY": {
        "rdkit.exact_mol_wt", "rdkit.clogp", "rdkit.tpsa", "rdkit.hbd",
        "rdkit.hba", "rdkit.rotatable_bonds", "rdkit.ring_count",
        "rdkit.aromatic_ring_count", "rdkit.fraction_csp3",
        "rdkit.formal_charge", "rdkit.heavy_atom_count",
        "rdkit.chiral_centers", "lipinski.ro5_violations", "rdkit.qed",
    },
    "ALERT": {
        "rdkit.pains_alert_count", "rdkit.brenk_alert_count",
        "rdkit.nih_alert_count",
        *(f"chemfh.{endpoint_id}" for endpoint_id, _name in _CHEMFH_ENDPOINTS),
        "chemfh.substructure_alert_count",
    },
    "SYN": {
        "rdkit.sa_score", "aizynthfinder.route_found",
        "aizynthfinder.state_score", "aizynthfinder.reaction_count",
        "aizynthfinder.precursor_count",
        "aizynthfinder.precursors_in_stock_complete",
    },
    "ABS": {
        f"admet_ai.{endpoint}"
        for endpoint in _ADMET_AI_ENDPOINTS_BY_CATEGORY["ABS"]
    },
    "DIST": {
        f"admet_ai.{endpoint}"
        for endpoint in _ADMET_AI_ENDPOINTS_BY_CATEGORY["DIST"]
    },
    "MET": {
        *(f"admet_ai.{endpoint}"
          for endpoint in _ADMET_AI_ENDPOINTS_BY_CATEGORY["MET"]),
        *_ADMET_AI_DERIVED_METRIC_IDS_BY_CATEGORY["MET"],
    },
    "EXCR": {
        f"admet_ai.{endpoint}"
        for endpoint in _ADMET_AI_ENDPOINTS_BY_CATEGORY["EXCR"]
    },
    "TOX": {
        f"admet_ai.{endpoint}"
        for endpoint in _ADMET_AI_ENDPOINTS_BY_CATEGORY["TOX"]
    },
    "BENCH": {
        "redocking.native_ligand_rmsd", "screening_benchmark.EF1pct",
        "screening_benchmark.AUROC",
    },
}

# Additional endpoints are assigned before constructing their shared rules.
_CATEGORY_METRIC_IDS["INTER"].update({"arpeggio.clash_contacts", "arpeggio.vdw_clash_contacts", "arpeggio.proximal_contact_fraction"})
_CATEGORY_METRIC_IDS.setdefault('SYN', set()).update(['rascore.score'])
_CATEGORY_METRIC_IDS.setdefault('INTRA', set()).update(['rdkit_pose.strain_energy', 'rdkit_pose.tfd', 'rdkit_pose.delta_rmsd'])
_CATEGORY_METRIC_IDS.setdefault('POCKET', set()).update(['freesasa.buriedness_fraction', 'freesasa.contact_surface_fraction'])
_CATEGORY_METRIC_IDS.setdefault('FP', set()).update(['prolif.tanimoto_to_reference'])
_CATEGORY_METRIC_IDS.setdefault('DYN', set()).update(['mdanalysis.ligand_rmsd', 'mdanalysis.pocket_com_drift', 'mdanalysis.contact_retention', 'rtmscore.ensemble_score'])

_TRAPEZOID_METRIC_IDS = {
    "rdkit.exact_mol_wt", "rdkit.clogp", "rdkit.tpsa",
    "rdkit.formal_charge", "rdkit.heavy_atom_count",
    "admet_ai.Lipophilicity_AstraZeneca", "admet_ai.VDss_Lombardo",
}

_PLAN_TOOL_BY_PREFIX = {
    "autodock_vina.": "autodock_vina",
    "unidock.": "unidock",
    "gnina.": "gnina",
    "diffdock.": "diffdock",
    "dock6.": "dock6",
    "autodock4.": "autodock4",
    "rdock.": "rdock",
    "adfr.": "adfr",
    "posebusters.": "posebusters",
    "posecheck.": "posecheck",
    "plip.": "plip",
    "arpeggio.": "arpeggio",
    "prolif.": "prolif",
    "blueprint_match.": "blueprint_match",
    "pocket_occupancy.": "pocket_occupancy",
    "mmgbsa.": "amber_mmpbsa",
    "aizynthfinder.": "aizynthfinder",
    "chemfh.": "chemfh",
}

_PLAN_TOOL_BY_PREFIX.update({'rascore.': 'rascore', 'rdkit_pose.': 'rdkit_pose', 'freesasa.': 'freesasa', 'prolif.': 'prolif', 'mdanalysis.': 'mdanalysis', 'rtmscore.': 'rtmscore'})


def _metric_category(metric_id: str) -> str:
    matches = [
        category
        for category, metric_ids in _CATEGORY_METRIC_IDS.items()
        if metric_id in metric_ids
    ]
    if len(matches) != 1:
        raise ValueError(
            f"{metric_id}: expected exactly one scoring category, got {matches}"
        )
    return matches[0]


def _plan_tool(metric_id: str, required_property_tool: str | None) -> str | None:
    if required_property_tool:
        return str(required_property_tool)
    for prefix, tool_id in _PLAN_TOOL_BY_PREFIX.items():
        if metric_id.startswith(prefix):
            return tool_id
    return None


def _rule(
    metric_id: str,
    tool: str,
    output_metric: str,
    *,
    extractor: dict[str, Any],
    hard_lower: float | None = None,
    soft_lower: float | None = None,
    soft_upper: float | None = None,
    hard_upper: float | None = None,
    weight: float = 0.0,
    description: str,
    direction: str | None = None,
    molecule_effect: str = "gate",
    stage: str = "pose_evaluation",
    family_id: str | None = None,
    family_weight_cap: float | None = None,
    requires_property_tool: str | None = None,
    covalent_review: bool = False,
    activates_category: bool = True,
) -> dict[str, Any]:
    if direction is None:
        has_lower = hard_lower is not None or soft_lower is not None
        has_upper = hard_upper is not None or soft_upper is not None
        direction = (
            "target_range"
            if has_lower and has_upper
            else "higher"
            if has_lower
            else "lower"
            if has_upper
            else None
        )
    category_code = _metric_category(metric_id)
    category = CATEGORY_CONFIG[category_code]
    weight_magnitude = abs(float(weight))
    signed_weight = (
        0.0
        if weight_magnitude == 0
        else -weight_magnitude if direction == "lower" else weight_magnitude
    )
    gate_participates = any(
        bound is not None
        for bound in (hard_lower, soft_lower, soft_upper, hard_upper)
    )
    return {
        "metric_id": metric_id,
        "stage": stage,
        "tool": tool,
        "output_metric": output_metric,
        "extractor": extractor,
        "hard_lower": hard_lower,
        "soft_lower": soft_lower,
        "soft_upper": soft_upper,
        "hard_upper": hard_upper,
        "weight": float(weight),
        "signed_weight": signed_weight,
        "description": description,
        "direction": direction,
        "rank_transform": (
            "trapezoid" if metric_id in _TRAPEZOID_METRIC_IDS else "identity"
        ),
        "category_code": category_code,
        "step_id": category["step_id"],
        "category_weight": category["parent_weight"],
        "catalog_abs_weight": category["catalog_abs_weight"],
        "plan_tool": _plan_tool(metric_id, requires_property_tool),
        "molecule_effect": molecule_effect,
        "gate_participates": gate_participates,
        "family_id": family_id,
        "family_weight_cap": family_weight_cap,
        "requires_property_tool": requires_property_tool,
        "covalent_review": bool(covalent_review),
        "activates_category": bool(activates_category),
    }


# These numeric bounds define the current executable metric configuration.
POSE_EVALUATION_METRIC_RULES: tuple[dict[str, Any], ...] = (
    _rule(
        "autodock_vina.vina_affinity", "AutoDock Vina",
        "docking score / affinity [R]",
        extractor={"kind": "pose_score", "metric_id": "autodock_vina.vina_affinity"},
        soft_upper=-7.0, hard_upper=-1.0, weight=5,
        description=(
            "Fallback within one target/protocol: <=-7 pass; -7 to -1 flag; "
            ">-1 fail. Cohort normalization remains the primary ranking "
            "treatment."
        ),
        stage="docking",
        family_id="vina_family_energy",
        family_weight_cap=5,
    ),
    _rule(
        "unidock.docking_score", "Uni-Dock", "docking score [R]",
        extractor={"kind": "pose_score", "metric_id": "unidock.docking_score"},
        soft_upper=-7.0, hard_upper=-1.0, weight=3,
        description=(
            "Vina-like fallback within one scoring function/protocol: <=-7 "
            "pass; -7 to -1 flag; >-1 fail."
        ),
        stage="docking",
        family_id="vina_family_energy",
        family_weight_cap=5,
    ),
    _rule(
        "gnina.minimizedAffinity", "GNINA", "affinity [R]",
        extractor={"kind": "pose_score", "metric_id": "gnina.minimizedAffinity"},
        soft_upper=-7.0, hard_upper=-1.0, weight=5,
        description=(
            "GNINA Vina-like minimized affinity fallback: <=-7 pass; -7 to -1 "
            "flag; >-1 fail; it is not CNNaffinity."
        ),
        stage="docking",
        family_id="vina_family_energy",
        family_weight_cap=5,
    ),
    _rule(
        "diffdock.confidence", "DiffDock", "confidence-model score [R]",
        extractor={"kind": "pose_score", "metric_id": "diffdock.confidence"},
        weight=0, direction="higher",
        description=(
            "No universal absolute cutoff: finite confidence values are ranked "
            "only within one receptor/model/checkpoint/protocol cohort. This "
            "runtime-only signal is audit-only because the canonical catalog "
            "does not assign it a child weight."
        ),
        stage="docking",
    ),
    _rule(
        "dock6.grid_score", "DOCK 6", "Grid_Score [R]",
        extractor={"kind": "pose_score", "metric_id": "dock6.grid_score"},
        weight=0, direction="lower",
        description=(
            "No universal absolute cutoff: finite native grid scores are ranked "
            "lower-is-better only within one sphere/grid/scoring-protocol cohort. "
            "It is audit-only because the canonical catalog has no DOCK 6 row."
        ),
        stage="docking",
    ),
    _rule(
        "autodock4.binding_energy", "AutoDock 4",
        "estimated free energy of binding [R]",
        extractor={
            "kind": "pose_score", "metric_id": "autodock4.binding_energy",
        },
        soft_upper=-7.0, hard_upper=-1.0, weight=4,
        description=(
            "Fallback within one grid/preparation/protocol cohort: <=-7 pass; "
            "-7 to -1 flag; >-1 fail. It is an estimated docking energy, not "
            "experimental binding free energy."
        ),
        stage="docking",
    ),
    _rule(
        "rdock.score", "rDock", "SCORE / total [R]",
        extractor={"kind": "pose_score", "metric_id": "rdock.score"},
        weight=5, direction="lower",
        description=(
            "No Vina-like absolute cutoff: rank lower-is-better only within one "
            "rDock scoring function, cavity, restraint, and parameter protocol."
        ),
        stage="docking",
    ),
    _rule(
        "adfr.affinity", "ADFR", "native FEB affinity [R]",
        extractor={"kind": "pose_score", "metric_id": "adfr.affinity"},
        soft_upper=-6.0, hard_upper=0.0, weight=4, direction="lower",
        description=(
            "AutoDockFR bounds: <=-6 pass; -6 to 0 flag; >0 fail, ranked "
            "lower-is-better within one protocol cohort. A flexible-receptor "
            "protocol requires a fixed flexible-residue set, while the registered runtime "
            "ADFR is rigid-receptor, so treat the value as a same-protocol rank "
            "and a coarse sanity bound rather than a flexible-receptor result."
        ),
        stage="docking",
    ),
    _rule(
        "autodock_vina.ligand_efficiency", "AutoDock Vina",
        "ligand efficiency [R]",
        extractor={
            "kind": "pose_score",
            "metric_id": "autodock_vina.ligand_efficiency",
            "scale": -1.0,
        },
        soft_upper=-0.25, hard_upper=-0.05, weight=2,
        description=(
            "Native score is -affinity/HAC; multiply by -1 for the numeric "
            "table's affinity/HAC convention: <=-0.25 pass, >-0.05 fail."
        ),
        stage="docking",
    ),
    _rule(
        "autodock_vina.vina_inter", "AutoDock Vina", "vina_inter",
        extractor={
            "kind": "pose_score",
            "metric_id": "autodock_vina.vina_inter",
        },
        hard_upper=0, weight=0, direction="lower",
        description=(
            "Intermolecular term of vina_affinity in kcal/mol. It is a "
            "sub-term of an already-weighted total, so it never adds score; "
            "a positive value is physically implausible for a bound pose and "
            "fails this metric so the pose can be screened out."
        ),
        stage="docking",
    ),
    _rule(
        "autodock_vina.vina_intra_best_pose", "AutoDock Vina",
        "vina_intra_best_pose",
        extractor={
            "kind": "pose_score",
            "metric_id": "autodock_vina.vina_intra_best_pose",
        },
        hard_upper=0, weight=0, direction="lower",
        description=(
            "Intramolecular energy of the best pose in kcal/mol, i.e. the "
            "quality ceiling this ligand can reach. Score-excluded; a positive "
            "value fails so a strained best pose is caught instead of being "
            "silently accepted."
        ),
        stage="docking",
    ),
    _rule(
        "autodock_vina.engine_pose_rank", "AutoDock Vina", "pose rank",
        extractor={
            "kind": "pose_score",
            "metric_id": "autodock_vina.engine_pose_rank",
        },
        soft_upper=20, hard_upper=20, weight=0,
        description=(
            "One-based rank inside one ligand-preparation run; ranks above 20 "
            "hard-fail the pose export set, are demoted to a flag for the "
            "parent molecule, and never add score."
        ),
        stage="docking",
        molecule_effect="pose_export",
    ),
    _rule(
        "autodock_vina.score_gap_1_2", "AutoDock Vina", "score gap",
        extractor={
            "kind": "pose_score",
            "metric_id": "autodock_vina.score_gap_1_2",
        },
        soft_lower=0.5, weight=0,
        description=(
            "score_2-score_1 within one preparation: >=0.5 pass; <0.5 flags "
            "as uncertain only, with no hard molecular failure."
        ),
        stage="docking",
    ),
    _rule(
        "gnina.intramol", "GNINA", "intramol",
        extractor={"kind": "pose_score", "metric_id": "gnina.intramol"},
        weight=0, direction="lower",
        description="Optional GNINA intramolecular term; audit-only strain context.",
        stage="docking",
    ),
    _rule(
        "gnina.CNNscore", "GNINA", "CNNpose score [R]",
        extractor={"kind": "pose_score", "metric_id": "gnina.CNNscore"},
        hard_lower=0.20, soft_lower=0.50, weight=4,
        description=(
            ">=0.50 pass; 0.20-0.50 flag; <0.20 fail. The value is a "
            "pose-quality score, not a calibrated probability."
        ),
        stage="docking",
    ),
    _rule(
        "gnina.CNNaffinity", "GNINA", "CNNaffinity [R]",
        extractor={"kind": "pose_score", "metric_id": "gnina.CNNaffinity"},
        hard_lower=2.0, soft_lower=6.0, weight=3,
        description=(
            "Fallback pK-scale gate: >=6 pass; 2-6 flag; <2 fail; never "
            "combine its raw value with kcal/mol affinity."
        ),
        stage="docking",
    ),
    _rule(
        "gnina.CNNvariance", "GNINA", "CNNvariance [R]",
        extractor={"kind": "pose_score", "metric_id": "gnina.CNNvariance"},
        weight=0, direction="lower",
        description=(
            "Optional ensemble uncertainty; audit/rank only within the same "
            "GNINA ensemble and model version."
        ),
        stage="docking",
    ),
    _rule(
        "gnina.engine_pose_rank", "GNINA", "pose rank by configured sort",
        extractor={
            "kind": "pose_score", "metric_id": "gnina.engine_pose_rank",
        },
        soft_upper=20, hard_upper=20, weight=0,
        description=(
            "One-based GNINA output rank under the configured pose_sort_order; "
            "ranks above 20 hard-fail the pose export set, are demoted to a "
            "flag for the parent molecule, and never add score."
        ),
        stage="docking",
        molecule_effect="pose_export",
    ),
    _rule(
        "posebusters.sanitization", "PoseBusters", "sanitization",
        extractor={"kind": "posebusters_check", "check": "sanitization"},
        hard_lower=1, soft_lower=1, weight=0,
        description="pass=1; fail=0 marks the pose as chemically unreadable.",
    ),
    _rule(
        "posebusters.inchi_convertible", "PoseBusters", "inchi_convertible",
        extractor={"kind": "posebusters_check", "check": "inchi_convertible"},
        hard_lower=1, soft_lower=1, weight=0,
        description="pass=1; fail=0 marks the pose and rolls back molecule preparation.",
    ),
    _rule(
        "posebusters.all_atoms_connected", "PoseBusters", "all_atoms_connected",
        extractor={"kind": "posebusters_check", "check": "all_atoms_connected"},
        hard_lower=1, soft_lower=1, weight=0,
        description="pass=1; fail=0 unless an allowed multicomponent route is explicit.",
    ),
    _rule(
        "posebusters.overall_success_ratio", "PoseBusters", "overall success ratio",
        extractor={"kind": "posebusters_success_ratio"},
        hard_lower=0.75, soft_lower=0.90, weight=3,
        description=">=0.90 pass; 0.75-0.90 flag; <0.75 fail; critical checks still label independently.",
    ),
    _rule(
        "posebusters.bond_lengths", "PoseBusters", "bond_lengths",
        extractor={"kind": "posebusters_check", "check": "bond_lengths"},
        hard_lower=1, soft_lower=1, weight=1,
        description="pass=1; fail=0.",
    ),
    _rule(
        "posebusters.bond_angles", "PoseBusters", "bond_angles",
        extractor={"kind": "posebusters_check", "check": "bond_angles"},
        hard_lower=1, soft_lower=1, weight=1,
        description="pass=1; fail=0.",
    ),
    _rule(
        "posebusters.internal_steric_clash", "PoseBusters", "internal_steric_clash",
        extractor={"kind": "posebusters_check", "check": "internal_steric_clash"},
        hard_lower=1, soft_lower=1, weight=2,
        description="pass=1; fail=0.",
    ),
    _rule(
        "posebusters.aromatic_ring_flatness", "PoseBusters", "aromatic_ring_flatness",
        extractor={"kind": "posebusters_check", "check": "aromatic_ring_flatness"},
        hard_lower=1, soft_lower=1, weight=1,
        description="pass=1; fail=0.",
    ),
    _rule(
        "posebusters.non_aromatic_ring_non_flatness", "PoseBusters", "non-aromatic_ring_non-flatness",
        extractor={"kind": "posebusters_check", "check": "non_aromatic_ring_non_flatness"},
        soft_lower=1, weight=1,
        description="pass=1; fail=0 only flags for local minimization, never a hard fail.",
    ),
    _rule(
        "posebusters.double_bond_flatness", "PoseBusters", "double_bond_flatness",
        extractor={"kind": "posebusters_check", "check": "double_bond_flatness"},
        hard_lower=1, soft_lower=1, weight=1,
        description="pass=1; fail=0.",
    ),
    _rule(
        "posebusters.internal_energy", "PoseBusters", "internal_energy",
        extractor={"kind": "posebusters_check", "check": "internal_energy"},
        soft_lower=1, weight=1,
        description="pass=1; fail=0 flags only and is confirmed with an independent strain calculation.",
    ),
    _rule(
        "posebusters.receptor_ligand_maximum_distance", "PoseBusters", "receptor-ligand_maximum_distance",
        extractor={"kind": "posebusters_check", "check": "receptor_ligand_maximum_distance"},
        hard_lower=1, soft_lower=1, weight=1,
        description="pass=1; fail=0 marks a pose far outside the receptor/pocket.",
    ),
    _rule(
        "posebusters.minimum_distance_to_receptor", "PoseBusters", "minimum_distance_to_receptor",
        extractor={"kind": "posebusters_check", "check": "minimum_distance_to_receptor"},
        hard_lower=1, soft_lower=1, weight=2,
        description="pass=1; fail=0 marks a severe receptor clash.",
    ),
    _rule(
        "posebusters.minimum_distance_to_organic_cofactors", "PoseBusters", "minimum_distance_to_organic_cofactors",
        extractor={"kind": "posebusters_check", "check": "minimum_distance_to_organic_cofactors"},
        soft_lower=1, weight=1,
        description="pass=1; fail=0 flags for review under the supplied numeric bounds.",
    ),
    _rule(
        "posebusters.minimum_distance_to_inorganic_cofactors", "PoseBusters", "minimum_distance_to_inorganic_cofactors",
        extractor={"kind": "posebusters_check", "check": "minimum_distance_to_inorganic_cofactors"},
        soft_lower=1, weight=1,
        description="pass=1; metal systems require a dedicated coordination profile.",
    ),
    _rule(
        "posebusters.minimum_distance_to_waters", "PoseBusters", "minimum_distance_to_waters",
        extractor={"kind": "posebusters_check", "check": "minimum_distance_to_waters"},
        soft_lower=1, weight=1,
        description="fail=0 flags only; a confirmed conserved-water clash is escalated by review, not by this bound.",
    ),
    _rule(
        "posebusters.volume_overlap_with_receptor", "PoseBusters", "volume_overlap_with_receptor",
        extractor={"kind": "posebusters_check", "check": "volume_overlap_with_receptor"},
        hard_lower=1, soft_lower=1, weight=2,
        description="pass=1; fail=0.",
    ),
    _rule(
        "posebusters.volume_overlap_with_organic_cofactors", "PoseBusters", "volume_overlap_with_organic_cofactors",
        extractor={"kind": "posebusters_check", "check": "volume_overlap_with_organic_cofactors"},
        soft_lower=1, weight=1,
        description="pass=1; fail=0 flags for review under the supplied numeric bounds.",
    ),
    _rule(
        "posebusters.volume_overlap_with_inorganic_cofactors", "PoseBusters", "volume_overlap_with_inorganic_cofactors",
        extractor={"kind": "posebusters_check", "check": "volume_overlap_with_inorganic_cofactors"},
        soft_lower=1, weight=1,
        description="pass=1; ordinary overlap failure only flags for metal systems.",
    ),
    _rule(
        "posebusters.volume_overlap_with_waters", "PoseBusters", "volume_overlap_with_waters",
        extractor={"kind": "posebusters_check", "check": "volume_overlap_with_waters"},
        soft_lower=1, weight=1,
        description="fail=0 flags only; hard treatment requires a non-displaceable water.",
    ),
    _rule(
        "posecheck.strain_energy_kcal", "PoseCheck", "strain energy [R]",
        extractor={"kind": "path", "path": "posecheck.strain_energy_kcal"},
        soft_upper=8, hard_upper=60, weight=4,
        description="<=8 kcal/mol pass; 8-60 flag; >60 fail.",
    ),
    _rule(
        "posecheck.steric_clash_count", "PoseCheck", "steric clash count",
        extractor={"kind": "path", "path": "posecheck.steric_clash_count"},
        soft_upper=0, hard_upper=2, weight=5,
        description="0 pass; 1-2 flag; >2 fail.",
    ),
    _rule(
        "plip.num_interactions", "PLIP", "num_interactions",
        extractor={"kind": "plip_num_interactions"},
        hard_lower=2, weight=0,
        description=(
            "Two or more interactions pass; fewer than two fail. "
            "The count remains weight zero and is not a binding-quality score."
        ),
    ),
    _rule(
        "prolif.num_interactions", "ProLIF", "interaction fingerprint observations",
        extractor={"kind": "prolif_num_interactions"},
        weight=0, direction="higher", molecule_effect="record",
        activates_category=False,
        description="Audit only; fingerprint bit count is not added to the score.",
    ),
    _rule(
        "blueprint_match.verdict", "blueprint_match", "verdict",
        extractor={"kind": "blueprint_verdict"},
        hard_lower=1, soft_lower=2, weight=5,
        description=(
            "pass=2; warn=1 flags for review; fail=0. "
            "unresolved carries no value: an absent or crashed profiler is a gap "
            "in evidence, so the metric leaves both numerator and denominator "
            "."
        ),
    ),
    _rule(
        "blueprint_match.required_hit_rate", "blueprint_match", "required_hit_rate [R]",
        extractor={"kind": "path", "path": "blueprint_match.required_hit_rate"},
        soft_lower=0.8, weight=0,
        description=(
            ">=0.80 passes and lower coverage flags; exact required-group logic "
            "remains the independent hard bound."
        ),
    ),
    _rule(
        "blueprint_match.desired_hit_rate", "blueprint_match", "desired_hit_rate",
        extractor={"kind": "path", "path": "blueprint_match.desired_hit_rate"},
        soft_lower=0.5, weight=0,
        description="Auxiliary explanation metric; below 0.5 requests modification.",
    ),
    _rule(
        "blueprint_match.desired_ranking_score", "blueprint_match", "desired_ranking_score [R]",
        extractor={"kind": "path", "path": "blueprint_match.desired_ranking_score"},
        hard_lower=0.1, soft_lower=0.5, weight=4,
        description=">=0.50 pass; 0.10-0.50 flag; <0.10 fail.",
    ),
    _rule(
        "blueprint_match.desired_evaluable_coverage", "blueprint_match", "desired_evaluable_coverage [R]",
        extractor={"kind": "path", "path": "blueprint_match.desired_evaluable_coverage"},
        soft_lower=0.8, weight=5,
        molecule_effect="protocol",
        description="Numeric table supplies only soft_lower=0.80; coverage flags request review.",
    ),
    _rule(
        "blueprint_match.optional_hits", "blueprint_match", "optional_hits [R]",
        extractor={"kind": "blueprint_optional_hits"},
        weight=1, direction="higher",
        description="Final tie-break evidence within the same blueprint.",
    ),
    _rule(
        "pocket_occupancy.hydrophobic_region", "pocket_occupancy", "hydrophobic_region",
        extractor={"kind": "path", "path": "pocket_contact_fraction.hydrophobic_region"},
        hard_lower=0.01, soft_lower=0.20, weight=0,
        description="Numeric bounds label this 0-1 fraction; weight zero excludes scoring.",
    ),
    _rule(
        "pocket_occupancy.polar_region", "pocket_occupancy", "polar_region",
        extractor={"kind": "path", "path": "pocket_contact_fraction.polar_region"},
        hard_lower=0.01, soft_lower=0.20, weight=0,
        description="Numeric bounds label this 0-1 fraction; weight zero excludes scoring.",
    ),
    _rule(
        "pocket_occupancy.charged_region", "pocket_occupancy", "charged_region",
        extractor={"kind": "path", "path": "pocket_contact_fraction.charged_region"},
        hard_lower=0.01, soft_lower=0.20, weight=0,
        description="Numeric bounds label this 0-1 fraction; weight zero excludes scoring.",
    ),
    _rule(
        "pocket_occupancy.aromatic_region", "pocket_occupancy", "aromatic_region",
        extractor={"kind": "path", "path": "pocket_contact_fraction.aromatic_region"},
        hard_lower=0.01, soft_lower=0.20, weight=0,
        description="Numeric bounds label this 0-1 fraction; weight zero excludes scoring.",
    ),
    _rule(
        "pocket_occupancy.overall", "pocket_occupancy", "overall [R]",
        extractor={"kind": "path", "path": "pocket_contact_fraction.overall"},
        hard_lower=0.05, soft_lower=0.20, weight=3,
        description=">=0.20 pass; 0.05-0.20 flag; <0.05 fail; no reference profile is NA.",
    ),
    _rule(
        "mmgbsa.delta_g_bind", "MM-GBSA", "deltaG_bind [R]",
        extractor={"kind": "pose_score", "metric_id": "amber_mmpbsa.delta_total_no_entropy"},
        soft_upper=0, hard_upper=10, weight=4,
        description="Same-protocol rank; <=0 pass, 0-10 flag, >10 fail as a weak sanity bound.",
    ),
)


# ADMET-AI endpoint IDs below match the installed provider schema exactly.  Do
# not map similarly named ADMETlab/SwissADME/pkCSM rows onto this backend.  The
# two count metrics are deterministic derivatives of the named ADMET-AI
# probabilities; every other ADMET-AI metric reads one native endpoint.
def _admet_ai_rule(
    endpoint: str,
    *,
    hard_lower: float | None = None,
    soft_lower: float | None = None,
    soft_upper: float | None = None,
    hard_upper: float | None = None,
    weight: float = 0.0,
    direction: str | None = None,
    description: str,
    extractor: dict[str, Any] | None = None,
    activates_category: bool = True,
) -> dict[str, Any]:
    metric_id = (
        endpoint
        if endpoint.startswith("admet_ai.")
        else f"admet_ai.{endpoint}"
    )
    output_metric = endpoint.removeprefix("admet_ai.")
    return _rule(
        metric_id,
        "ADMET-AI",
        output_metric,
        extractor=extractor or {
            "kind": "admet_ai_endpoint",
            "endpoint_id": output_metric,
        },
        hard_lower=hard_lower,
        soft_lower=soft_lower,
        soft_upper=soft_upper,
        hard_upper=hard_upper,
        weight=weight,
        direction=direction,
        stage="property_filter",
        requires_property_tool="admet",
        activates_category=activates_category,
        description=description,
    )


ADMET_AI_METRIC_RULES: tuple[dict[str, Any], ...] = (
    # Absorption
    _admet_ai_rule(
        "HIA_Hou", hard_lower=0.1, soft_lower=0.7, weight=4,
        description=">=0.70 pass; 0.10-0.70 flag; <0.10 fail.",
    ),
    _admet_ai_rule(
        "Caco2_Wang", hard_lower=-7, soft_lower=-5.15, weight=3.5,
        description=(
            "log10(Papp [cm/s]): >=-5.15 pass; -7 to -5.15 flag; <-7 fail."
        ),
    ),
    _admet_ai_rule(
        "PAMPA_NCATS", hard_lower=0.1, soft_lower=0.7, weight=2.5,
        description=">=0.70 pass; 0.10-0.70 flag; <0.10 fail.",
    ),
    _admet_ai_rule(
        "Pgp_Broccatelli", soft_upper=0.3, hard_upper=0.9, weight=3.5,
        description="<=0.30 pass; 0.30-0.90 flag; >0.90 fail.",
    ),
    _admet_ai_rule(
        "Bioavailability_Ma", hard_lower=0.1, soft_lower=0.7, weight=4,
        description=">=0.70 pass; 0.10-0.70 flag; <0.10 fail.",
    ),
    _admet_ai_rule(
        "Solubility_AqSolDB", hard_lower=-9, soft_lower=-4, weight=4,
        description=(
            "logS (mol/L): >=-4 pass; -9 to -4 flag; <-9 fail."
        ),
    ),
    _admet_ai_rule(
        "Lipophilicity_AstraZeneca",
        hard_lower=-2, soft_lower=-1, soft_upper=3.5, hard_upper=5,
        weight=2.5,
        description=(
            "Target-range transform: -1 to 3.5 pass; -2 to -1 or 3.5 to 5 "
            "flag; outside -2 to 5 fail."
        ),
    ),
    _admet_ai_rule(
        "HydrationFreeEnergy_FreeSolv", weight=0, direction="higher",
        description="Audit-only hydration free energy; no universal gate or weight.",
    ),

    # Distribution (default oral non-CNS profile).
    _admet_ai_rule(
        "BBB_Martins", soft_upper=0.3, hard_upper=0.9, weight=2,
        description=(
            "Default non-CNS profile: <=0.30 pass; 0.30-0.90 flag; >0.90 "
            "fail. A CNS profile must explicitly replace this direction."
        ),
    ),
    _admet_ai_rule(
        "PPBR_AZ", soft_upper=90, hard_upper=98, weight=2,
        description="Percent bound: <=90 pass; 90-98 flag; >98 fail.",
    ),
    _admet_ai_rule(
        "VDss_Lombardo", hard_lower=0, soft_lower=0.04, soft_upper=20,
        weight=1,
        description=(
            "Target-range transform in L/kg for domain-valid evidence: "
            "0.04-20 pass; 0-0.04 or values above 20 flag. Negative values "
            "are invalid evidence (NA), so they do not gate or rank."
        ),
    ),

    # Metabolism.  The three substrate metrics use the explicitly selected
    # negative policy: higher substrate probability is a penalty.  This is a
    # consequence of the signed scoring convention implemented here.
    _admet_ai_rule(
        "admet_ai.cyp_substrate_count", soft_lower=1, weight=0,
        extractor={
            "kind": "admet_ai_count",
            "endpoint_ids": [
                "CYP2C9_Substrate_CarbonMangels",
                "CYP2D6_Substrate_CarbonMangels",
                "CYP3A4_Substrate_CarbonMangels",
            ],
            "comparison": "lt",
            "threshold": 0.5,
        },
        description=(
            "Derived n_sum(probability < 0.5): >=1 pass; <1 flag; weight 0."
        ),
    ),
    _admet_ai_rule(
        "CYP2C9_Substrate_CarbonMangels", weight=0.1, direction="lower",
        description=(
            "No absolute gate. Explicit policy override: signed weight -0.1, "
            "so higher substrate probability is penalized."
        ),
    ),
    _admet_ai_rule(
        "CYP2D6_Substrate_CarbonMangels", weight=0.1, direction="lower",
        description=(
            "No absolute gate. Explicit policy override: signed weight -0.1, "
            "so higher substrate probability is penalized."
        ),
    ),
    _admet_ai_rule(
        "CYP3A4_Substrate_CarbonMangels", weight=0.5, direction="lower",
        description=(
            "No absolute gate. Explicit policy override: signed weight -0.5, "
            "so higher substrate probability is penalized."
        ),
    ),
    _admet_ai_rule(
        "admet_ai.cyp_inhibition_count",
        soft_upper=1, hard_upper=3, weight=4,
        extractor={
            "kind": "admet_ai_count",
            "endpoint_ids": [
                "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith",
                "CYP2D6_Veith", "CYP3A4_Veith",
            ],
            "comparison": "gt",
            "threshold": 0.7,
        },
        description=(
            "Derived n_sum(probability > 0.7): <=1 pass; 2-3 flag; >3 fail."
        ),
    ),
    *tuple(
        _admet_ai_rule(
            endpoint,
            soft_upper=0.3,
            hard_upper=0.9,
            weight=weight,
            description="<=0.30 pass; 0.30-0.90 flag; >0.90 fail.",
        )
        for endpoint, weight in (
            ("CYP1A2_Veith", 2.5), ("CYP2C19_Veith", 2.5),
            ("CYP2C9_Veith", 3), ("CYP2D6_Veith", 3.5),
            ("CYP3A4_Veith", 4),
        )
    ),

    # Excretion outputs are provider-native audit evidence. Their configured
    # weights and bounds are all zero/blank, so they do not by themselves
    # activate the weighted EXCR parent category.
    *tuple(
        _admet_ai_rule(
            endpoint,
            weight=0,
            direction="higher",
            activates_category=False,
            description="Audit-only provider output; no universal gate or weight.",
        )
        for endpoint in _ADMET_AI_ENDPOINTS_BY_CATEGORY["EXCR"]
    ),

    # Toxicity risks are lower-is-better and retain independent hard gates.
    *tuple(
        _admet_ai_rule(
            endpoint,
            soft_upper=soft_upper,
            hard_upper=hard_upper,
            weight=weight,
            description=(
                f"<={soft_upper:g} pass; {soft_upper:g}-{hard_upper:g} "
                f"flag; >{hard_upper:g} fail."
            ),
        )
        for endpoint, soft_upper, hard_upper, weight in (
            ("hERG", 0.2, 0.9, 5),
            ("ClinTox", 0.2, 0.9, 4),
            ("AMES", 0.2, 0.9, 5),
            ("DILI", 0.2, 0.9, 5),
            ("Carcinogens_Lagunin", 0.2, 0.9, 3.5),
            ("LD50_Zhu", 2.5, 6.5, 3),
            ("Skin_Reaction", 0.3, 0.9, 3),
            ("NR-AR", 0.2, 0.9, 1.5),
            ("NR-AR-LBD", 0.2, 0.9, 1.5),
            ("NR-AhR", 0.2, 0.9, 2),
            ("NR-Aromatase", 0.2, 0.9, 2),
            ("NR-ER", 0.2, 0.9, 2),
            ("NR-ER-LBD", 0.2, 0.9, 2),
            ("NR-PPAR-gamma", 0.2, 0.9, 1.5),
            ("SR-ARE", 0.2, 0.9, 2),
            ("SR-ATAD5", 0.2, 0.9, 2.5),
            ("SR-HSE", 0.2, 0.9, 1.5),
            ("SR-MMP", 0.2, 0.9, 3),
            ("SR-p53", 0.2, 0.9, 2.5),
        )
    ),
)


CHEMFH_METRIC_RULES: tuple[dict[str, Any], ...] = tuple(
    _rule(
        f"chemfh.{endpoint_id}",
        "ChemFH",
        endpoint_name,
        extractor={
            "kind": "chemfh_endpoint_probability",
            "endpoint_id": endpoint_id,
        },
        soft_upper=0.3,
        hard_upper=0.9,
        weight=_CHEMFH_ENDPOINT_WEIGHTS.get(endpoint_id, 0.0),
        direction="lower",
        stage="property_filter",
        requires_property_tool="chemfh",
        description=(
            "Native ChemFH risk probability: <=0.30 pass; 0.30-0.90 flag; "
            ">0.90 fail. Promiscuous/colloidal (|w|=0.5) and reactive (|w|=1) "
            "carry a low negative weight inside ALERT; the other four "
            "endpoints stay weight zero as gate and audit signals."
        ),
    )
    for endpoint_id, endpoint_name in _CHEMFH_ENDPOINTS
) + (
    _rule(
        "chemfh.substructure_alert_count",
        "ChemFH",
        "frequent hitter substructure alert",
        extractor={"kind": "chemfh_substructure_alert_count"},
        soft_upper=2,
        hard_upper=10,
        weight=2,
        direction="lower",
        stage="property_filter",
        requires_property_tool="chemfh",
        description=(
            "Severity of the ten native substructure-rule catalogs, one hit "
            "per rule regardless of SMARTS multiplicity, mapped none=0 / "
            "mild=2 / severe=10. Numeric sbdd contract: 0-2 pass; 3-10 flag; "
            ">10 fail. Since exactly ten rules are mapped, the tabulated hard "
            "bound is only reachable after a catalog-version change; the "
            "numeric columns take precedence over contradictory prose."
        ),
    ),
)


# Executable subset of the supplied property/ADMET/synthesis table. ChemFH
# emits provenance-bearing frequent-hitter evidence. Four of its seven model
# calls remain zero-weight gates; promiscuous/colloidal (|w|=0.5), reactive
# (|w|=1) and the mapped ten-rule substructure severity (|w|=2) carry the low
# negative ALERT weights of the executable policy, which together with the three RDKit
# catalog counts make up the category's Sigma|w|=7. Individual provider rule rows
# remain evidence-only (including the PAINS overlap). ADMETlab, SwissADME,
# RAScore, SCScore, SYBA, ASKCOS, and pkCSM rows are likewise not mapped onto a
# different backend merely because an endpoint name is similar.
PROPERTY_FILTER_METRIC_RULES: tuple[dict[str, Any], ...] = (
    _rule(
        "rdkit.exact_mol_wt", "RDKit descriptors", "MW / ExactMolWt",
        extractor={"kind": "property_path", "path": "properties.MW"},
        hard_lower=50, soft_lower=150, soft_upper=500, hard_upper=1050,
        weight=0, stage="property_filter",
        description=(
            "General oral profile: 150-500 pass; 50-150 or 500-1050 flag; "
            "outside 50-1050 fail. Fragment libraries need another profile."
        ),
    ),
    _rule(
        "rdkit.clogp", "RDKit descriptors", "cLogP",
        extractor={"kind": "property_path", "path": "properties.LogP"},
        hard_lower=-3, soft_lower=0, soft_upper=5, hard_upper=10,
        weight=3, stage="property_filter",
        description="0-5 pass; -3-0 or 5-10 flag; outside -3 to 10 fail.",
    ),
    _rule(
        "rdkit.tpsa", "RDKit descriptors", "TPSA",
        extractor={"kind": "property_path", "path": "properties.TPSA"},
        hard_lower=0, soft_lower=20, soft_upper=140, hard_upper=180,
        weight=3, stage="property_filter",
        description=(
            "General oral profile: 20-140 pass; 0-20 or 140-180 flag; "
            "outside 0-180 fail. CNS projects require a separate profile."
        ),
    ),
    _rule(
        "rdkit.hbd", "RDKit descriptors", "HBD",
        extractor={"kind": "property_path", "path": "properties.HBD"},
        soft_upper=5, hard_upper=10, weight=0, stage="property_filter",
        description="<=5 pass; 6-10 flag; >10 fail.",
    ),
    _rule(
        "rdkit.hba", "RDKit descriptors", "HBA",
        extractor={"kind": "property_path", "path": "properties.HBA"},
        soft_upper=10, hard_upper=20, weight=0, stage="property_filter",
        description="<=10 pass; 11-20 flag; >20 fail.",
    ),
    _rule(
        "rdkit.rotatable_bonds", "RDKit descriptors", "RotB",
        extractor={"kind": "property_path", "path": "properties.RotBonds"},
        soft_upper=10, hard_upper=15, weight=0, stage="property_filter",
        description="<=10 pass; 11-15 flag; >15 fail.",
    ),
    _rule(
        "rdkit.ring_count", "RDKit descriptors", "RingCount",
        extractor={"kind": "property_path", "path": "properties.RingCount"},
        soft_lower=2, soft_upper=4, hard_upper=20, weight=0,
        stage="property_filter",
        description=(
            "Oral drugs cluster at 2-4 rings with 3 the most balanced choice, "
            "so 2-4 passes, 0-1 or 5-20 flags, and >20 fails."
        ),
    ),
    _rule(
        "rdkit.aromatic_ring_count", "RDKit descriptors",
        "AromaticRingCount",
        extractor={
            "kind": "property_path", "path": "properties.AromaticRingCount",
        },
        soft_upper=4, hard_upper=6, weight=0, stage="property_filter",
        description="0-4 pass; 5-6 flag; >6 fail.",
    ),
    _rule(
        "rdkit.fraction_csp3", "RDKit descriptors", "FractionCSP3",
        extractor={
            "kind": "property_path", "path": "properties.FractionCSP3",
        },
        soft_lower=0.20, weight=0, stage="property_filter",
        description=(
            ">=0.20 pass; lower values flag for audit only and never hard-fail "
            "an intentionally planar ligand."
        ),
    ),
    _rule(
        "rdkit.formal_charge", "RDKit descriptors", "Formal charge",
        extractor={"kind": "property_path", "path": "properties.FormalCharge"},
        hard_lower=-3, soft_lower=-1, soft_upper=1, hard_upper=3,
        weight=0, stage="property_filter",
        description="Net charge -1 to +1 pass; +/-2 or +/-3 flag; outside fail.",
    ),
    _rule(
        "rdkit.heavy_atom_count", "RDKit descriptors", "HeavyAtomCount",
        extractor={
            "kind": "property_path", "path": "properties.HeavyAtomCount",
        },
        hard_lower=5, soft_lower=12, soft_upper=50, hard_upper=80,
        weight=0, stage="property_filter",
        description=(
            "12-50 pass; 5-12 or 50-80 flag; outside fail as library-range QC. "
            "It remains weight zero and also supports ligand efficiency."
        ),
    ),
    _rule(
        "rdkit.chiral_centers", "RDKit descriptors", "Chiral centers",
        extractor={"kind": "property_path", "path": "properties.ChiralCenters"},
        soft_upper=2, weight=0, stage="property_filter",
        description=(
            "Counts <=2 pass and higher counts flag for review; there is no "
            "hard failure from this count alone."
        ),
    ),
    _rule(
        "lipinski.ro5_violations", "Lipinski", "RO5 violation count",
        extractor={
            "kind": "property_path", "path": "lipinski.violations",
        },
        soft_upper=1, hard_upper=3, weight=0, stage="property_filter",
        requires_property_tool="lipinski",
        description=(
            "0-1 passes; 2-3 flags; >3 fails, following the authoritative "
            "numeric cells."
        ),
    ),
    _rule(
        "rdkit.qed", "RDKit.Chem.QED", "QED score",
        extractor={"kind": "property_path", "path": "properties.QED"},
        hard_lower=0.1, soft_lower=0.45, weight=3,
        stage="property_filter", requires_property_tool="qed",
        description="<0.1 fail; 0.1-0.45 flag; >=0.45 pass.",
    ),
    _rule(
        "rdkit.pains_alert_count", "RDKit PAINS", "PAINS alert count",
        extractor={
            "kind": "property_alert_count",
            "catalogs": ["PAINS_A", "PAINS_B", "PAINS_C", "PAINS"],
        },
        soft_upper=0, hard_upper=4, weight=1, direction="lower",
        stage="property_filter",
        requires_property_tool="pains", covalent_review=True,
        description=(
            "0 pass; 1-4 flag; >4 fail. If the project allows a rescue it must "
            "be confirmed by alert name rather than by count alone; covalent "
            "status alone is not a rescue."
        ),
    ),
    _rule(
        "rdkit.brenk_alert_count", "RDKit FilterCatalog", "Brenk alert count",
        extractor={
            "kind": "property_alert_count", "catalogs": ["BRENK"],
        },
        soft_upper=0, hard_upper=4, weight=1, direction="lower",
        stage="property_filter",
        requires_property_tool="pains", covalent_review=True,
        description=(
            "0 pass; 1-4 flag; >4 fail. A single highly reactive alert may be "
            "failed manually, and only an explicit catalog rescue policy can "
            "downgrade the hard result."
        ),
    ),
    _rule(
        "rdkit.nih_alert_count", "RDKit FilterCatalog", "NIH alert count",
        extractor={
            "kind": "property_alert_count", "catalogs": ["NIH"],
        },
        soft_upper=0, hard_upper=4, weight=1, direction="lower",
        stage="property_filter",
        requires_property_tool="pains", covalent_review=True,
        description=(
            "0 pass; 1-4 flag; >4 fail. Preserve the rule name and matched "
            "atoms; only an explicit catalog rescue policy can map a hit to "
            "review, and covalent status alone cannot."
        ),
    ),
    _rule(
        "rdkit.sa_score", "RDKit SA_score", "SA score [R]",
        extractor={"kind": "property_path", "path": "properties.SA_score"},
        soft_upper=6, hard_upper=9, weight=2, stage="property_filter",
        requires_property_tool="sa_score",
        description="<=6 pass; 6-9 flag; >9 fail as a synthesis heuristic.",
    ),
    _rule(
        "aizynthfinder.route_found", "AiZynthFinder", "route found",
        extractor={
            "kind": "property_path", "path": "synthesis.metrics.route_found",
        },
        soft_lower=1, weight=0, stage="property_filter",
        requires_property_tool="aizynthfinder",
        description=(
            "found=1 pass; 0 flags and extends search/stock/templates. One "
            "failed budget never hard-fails the molecule."
        ),
    ),
    _rule(
        "aizynthfinder.state_score", "AiZynthFinder", "State score",
        extractor={
            "kind": "property_path", "path": "synthesis.metrics.state_score",
        },
        weight=2, direction="higher", stage="property_filter",
        requires_property_tool="aizynthfinder",
        description=(
            "Current backend top-state score, normalized within one fixed "
            "AiZynthFinder version/config/search-budget cohort."
        ),
    ),
    _rule(
        "aizynthfinder.reaction_count", "AiZynthFinder",
        "Number of reactions [R]",
        extractor={
            "kind": "property_path", "path": "synthesis.metrics.reaction_count",
        },
        soft_upper=5, weight=0, stage="property_filter",
        requires_property_tool="aizynthfinder",
        description=(
            "<=5 pass; >5 flag; no hard-fail bound."
        ),
    ),
    _rule(
        "aizynthfinder.precursor_count", "AiZynthFinder",
        "Number of pre-cursors",
        extractor={
            "kind": "property_path", "path": "synthesis.metrics.precursor_count",
        },
        weight=0, direction="lower", stage="property_filter",
        requires_property_tool="aizynthfinder",
        description="Audit-only count from the exported best route.",
    ),
    _rule(
        "aizynthfinder.precursors_in_stock_complete", "AiZynthFinder",
        "All pre-cursors in stock",
        extractor={
            "kind": "property_precursors_in_stock_complete",
        },
        soft_lower=1, weight=0, stage="property_filter",
        requires_property_tool="aizynthfinder",
        description=(
            "1 passes only when precursor_count equals the in-stock count; a "
            "partial or zero-stock route only flags."
        ),
    ),
    *CHEMFH_METRIC_RULES,
    *ADMET_AI_METRIC_RULES,
)


_ADDITIONAL_METRIC_RULES = (
    _rule("arpeggio.clash_contacts", "Arpeggio", "clash_contacts",
          extractor={"kind": "path", "path": "arpeggio.metrics.clash_contacts"},
          soft_upper=0, hard_upper=2, weight=3, description="sbdd numeric row 820; validated candidate contacts."),
    _rule("arpeggio.vdw_clash_contacts", "Arpeggio", "vdw_clash_contacts",
          extractor={"kind": "path", "path": "arpeggio.metrics.vdw_clash_contacts"},
          soft_upper=8, hard_upper=15, weight=1, description="sbdd numeric row 821; validated candidate contacts."),
    _rule("arpeggio.proximal_contact_fraction", "Arpeggio", "proximal_contact_fraction",
          extractor={"kind": "path", "path": "arpeggio.metrics.proximal_contact_fraction"},
          hard_lower=0.2, soft_lower=0.8, weight=0, description="sbdd numeric row 822; proximal atom contacts / direct atom contacts."),
    _rule(
        'rascore.score', 'RAScore', 'rascore [R]',
        extractor={'kind': 'property_path', 'path': 'rascore.score'},
        hard_lower=0.2,
        soft_lower=0.5,
        weight=3.0, direction='higher',
        stage="property_filter", requires_property_tool="rascore",
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'rdkit_pose.strain_energy', 'RDKit UFF/MMFF', 'rdkit_strain_energy_kcal_mol [R]',
        extractor={'kind': 'path', 'path': 'rdkit_pose.metrics.strain_energy_kcal_mol'},
        soft_upper=15.0,
        hard_upper=80.0,
        weight=3.0, direction='lower',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'rdkit_pose.tfd', 'RDKit TFD', 'tfd [R]',
        extractor={'kind': 'path', 'path': 'rdkit_pose.metrics.tfd'},
        soft_upper=0.2,
        hard_upper=0.7,
        weight=2.0, direction='lower',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'rdkit_pose.delta_rmsd', 'RDKit conformer RMSD', 'delta_rmsd_A [R]',
        extractor={'kind': 'path', 'path': 'rdkit_pose.metrics.delta_rmsd_angstrom'},
        soft_upper=2.0,
        weight=1.0, direction='lower',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'freesasa.buriedness_fraction', 'RDKit FreeSASA', 'buriedness_fraction [R]',
        extractor={'kind': 'path', 'path': 'freesasa.metrics.buriedness_fraction'},
        hard_lower=0.1,
        soft_lower=0.3,
        weight=3.0, direction='higher',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'freesasa.contact_surface_fraction', 'RDKit FreeSASA', 'contact_surface_area_fraction [R]',
        extractor={'kind': 'path', 'path': 'freesasa.metrics.contact_surface_area_fraction'},
        hard_lower=0.1,
        soft_lower=0.3,
        weight=1.0, direction='higher',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'prolif.tanimoto_to_reference', 'ProLIF', 'prolif_tanimoto_to_reference [R]',
        extractor={'kind': 'path', 'path': 'prolif.tanimoto_to_reference'},
        weight=0.0, direction='higher', molecule_effect='record',
        activates_category=False,
        description='Audit-only reference fingerprint similarity; no gate or score weight.',
    ),
    _rule(
        'mdanalysis.ligand_rmsd', 'MDAnalysis', 'ligand RMSD [R]',
        extractor={'kind': 'path', 'path': 'mdanalysis.metrics.ligand_heavy_atom_rmsd_angstrom.median'},
        soft_upper=2.0,
        hard_upper=5.0,
        weight=4.0, direction='lower',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'mdanalysis.pocket_com_drift', 'MDAnalysis (supplemental)', 'ligand–pocket COM distance drift [R]',
        extractor={'kind': 'path', 'path': 'mdanalysis.metrics.ligand_pocket_com_distance_drift_angstrom.p95'},
        soft_upper=2.0,
        hard_upper=5.0,
        weight=3.0, direction='lower',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'mdanalysis.contact_retention', 'MDAnalysis', 'protein-ligand contacts [R]',
        extractor={'kind': 'path', 'path': 'mdanalysis.metrics.initial_contact_retention_fraction.mean'},
        hard_lower=0.2,
        soft_lower=0.5,
        weight=4.0, direction='higher',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
    _rule(
        'rtmscore.ensemble_score', 'RTMScore', 'score [R]',
        extractor={'kind': 'pose_score', 'metric_id': 'rtmscore.ensemble_score'},
        weight=3.0, direction='higher',
        description='See METRIC_POLICY.md for the producer definition.',
    ),
)
POSE_EVALUATION_METRIC_RULES += tuple(r for r in _ADDITIONAL_METRIC_RULES if r["stage"] == "pose_evaluation")
PROPERTY_FILTER_METRIC_RULES += tuple(r for r in _ADDITIONAL_METRIC_RULES if r["stage"] == "property_filter")

POSE_PROPERTY_METRIC_RULES: tuple[dict[str, Any], ...] = (
    *POSE_EVALUATION_METRIC_RULES,
    *PROPERTY_FILTER_METRIC_RULES,
)


# These values gate the target/pocket/protocol rather than a molecule.  They
# are published in the metric catalog and consumed by the docking/coordinator
# quality-gate paths, but are intentionally excluded from score_records().
DOCKING_PROTOCOL_METRIC_RULES: tuple[dict[str, Any], ...] = (
    _rule(
        "redocking.native_ligand_rmsd", "GNINA",
        "native ligand RMSD",
        extractor={"kind": "path", "path": "redocking_validation.best_rmsd_angstrom"},
        soft_upper=2.0, hard_upper=2.5, weight=0,
        molecule_effect="protocol",
        description=(
            "Runtime implementation is GNINA best-of-generated-modes "
            "symmetry-aware heavy-atom RMSD: <=2.0 A pass; 2.0-2.5 A "
            "flag; >2.5 A fails the protocol. The canonical table names "
            "Vina redocking, so matching numeric bounds do not imply adapter "
            "equivalence."
        ),
        stage="docking",
    ),
    _rule(
        "screening_benchmark.EF1pct", "screening_benchmark",
        "enrichment factor at 1%",
        extractor={"kind": "path", "path": "screening_validation.EF1%"},
        hard_lower=1.0, soft_lower=3.0, weight=0,
        molecule_effect="protocol",
        description=(
            "With labeled actives and matched decoys: EF1% >=3 pass; 1-3 "
            "flag; <1 fails the screening protocol."
        ),
        stage="coordination",
    ),
    _rule(
        "screening_benchmark.AUROC", "screening_benchmark", "ROC-AUC",
        extractor={"kind": "path", "path": "screening_validation.AUROC"},
        hard_lower=0.50, soft_lower=0.70, weight=0,
        molecule_effect="protocol",
        description=(
            "With labeled actives and decoys: AUROC >=0.70 pass; 0.50-0.70 "
            "flag; <0.50 fails the screening protocol."
        ),
        stage="coordination",
    ),
)


ALL_METRIC_RULES: tuple[dict[str, Any], ...] = (
    *POSE_PROPERTY_METRIC_RULES,
    *DOCKING_PROTOCOL_METRIC_RULES,
)


_POSEBUSTERS_ALIASES = {
    "passes_rdkit_sanity_checks": "sanitization",
    "bond_lengths_within_bounds": "bond_lengths",
    "bond_angles_within_bounds": "bond_angles",
    "no_internal_clash": "internal_steric_clash",
    "flatness_passes": "aromatic_ring_flatness",
    "energy_ratio_passes": "internal_energy",
    "protein_ligand_maximum_distance": "receptor_ligand_maximum_distance",
    "minimum_distance_to_protein": "minimum_distance_to_receptor",
    "volume_overlap_with_protein": "volume_overlap_with_receptor",
}


def canonical_posebusters_check_id(value: Any) -> str:
    """Normalize PoseBusters display/raw column names to stable metric IDs."""
    if isinstance(value, tuple) and value:
        value = value[-1]
    text = re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")
    return _POSEBUSTERS_ALIASES.get(text, text)


def validate_metric_rules(rules: Iterable[dict[str, Any]]) -> None:
    """Validate ordering, IDs, weights, directions, and effect semantics."""
    seen: set[str] = set()
    for rule in rules:
        metric_id = str(rule.get("metric_id") or "")
        if not metric_id or metric_id in seen:
            raise ValueError(f"metric_id must be non-empty and unique: {metric_id!r}")
        seen.add(metric_id)
        weight = _finite_number(rule.get("weight"))
        if weight is None or weight < 0 or weight > 5:
            raise ValueError(
                f"{metric_id}: weight must be a finite magnitude from 0 to 5"
            )
        bounds = [
            rule.get("hard_lower"), rule.get("soft_lower"),
            rule.get("soft_upper"), rule.get("hard_upper"),
        ]
        numeric = [_finite_number(value) for value in bounds]
        if any(value is None for value in numeric):
            invalid = [
                value for value, number in zip(bounds, numeric)
                if value is not None and number is None
            ]
            if invalid:
                raise ValueError(
                    f"{metric_id}: bounds must be finite numbers or null"
                )
        hard_lower, soft_lower, soft_upper, hard_upper = numeric
        ordered = (
            -math.inf if hard_lower is None else hard_lower,
            hard_lower if soft_lower is None and hard_lower is not None
            else -math.inf if soft_lower is None else soft_lower,
            hard_upper if soft_upper is None and hard_upper is not None
            else math.inf if soft_upper is None else soft_upper,
            math.inf if hard_upper is None else hard_upper,
        )
        if not (ordered[0] <= ordered[1] <= ordered[2] <= ordered[3]):
            raise ValueError(
                f"{metric_id}: require hard_lower <= soft_lower <= soft_upper <= hard_upper"
            )
        if rule.get("direction") not in (None, "higher", "lower", "target_range"):
            raise ValueError(
                f"{metric_id}: direction must be higher/lower/target_range/null"
            )
        if rule.get("molecule_effect", "gate") not in MOLECULE_EFFECTS:
            raise ValueError(f"{metric_id}: invalid molecule_effect")
        required_tool = rule.get("requires_property_tool")
        if required_tool is not None and not str(required_tool).strip():
            raise ValueError(
                f"{metric_id}: requires_property_tool must be non-empty or null"
            )
        family_id = rule.get("family_id")
        family_cap = _finite_number(rule.get("family_weight_cap"))
        if family_id and (family_cap is None or family_cap <= 0):
            raise ValueError(
                f"{metric_id}: family_weight_cap must be finite and > 0"
            )
        if family_cap is not None and not family_id:
            raise ValueError(
                f"{metric_id}: family_weight_cap requires family_id"
            )


def metric_label(value: Any, rule: dict[str, Any]) -> str:
    """Label one scalar against one rule's inclusive soft/hard bounds.

    This is the whole of level 1: the return value is a *mark* on
    one metric of one molecule -- ``pass``, ``flag``, ``fail``, or
    ``not_applicable`` when there is no finite value to judge.  It never decides
    whether the molecule survives; only the category fold in
    ``aggregate_molecule_gate`` does that.  A finite score with no declared
    hard/soft bound is score-only evidence, not an automatic ``pass`` mark.
    """
    number = _finite_number(value)
    if number is None:
        return NOT_APPLICABLE_LABEL
    if rule.get("gate_participates") is False:
        return NOT_APPLICABLE_LABEL
    hard_lower = _finite_number(rule.get("hard_lower"))
    soft_lower = _finite_number(rule.get("soft_lower"))
    soft_upper = _finite_number(rule.get("soft_upper"))
    hard_upper = _finite_number(rule.get("hard_upper"))
    if all(
        bound is None
        for bound in (hard_lower, soft_lower, soft_upper, hard_upper)
    ):
        return NOT_APPLICABLE_LABEL
    if (hard_lower is not None and number < hard_lower) or (
        hard_upper is not None and number > hard_upper
    ):
        return "fail"
    if (soft_lower is not None and number < soft_lower) or (
        soft_upper is not None and number > soft_upper
    ):
        return "flag"
    return "pass"


def gate_value(value: Any, rule: dict[str, Any]) -> str:
    """Legacy spelling of :func:`metric_label` (keep/modify/reject wording).

    Retained because published category tables, several skills and older reports
    key on these three strings.  The value is still only a per-metric mark: a
    ``"reject"`` here means "this one metric hard-failed", not "drop this
    molecule".
    """
    return LABEL_TO_LEGACY_GATE[metric_label(value, rule)]



def extract_metric_value(record: dict[str, Any], rule: dict[str, Any]) -> float | None:
    """Extract a finite numeric value for one declared metric rule."""
    if _metric_explicitly_not_applicable(record, rule):
        return None
    if not _required_property_tool_enabled(record, rule):
        return None
    if not _blueprint_metric_applicable(record, rule):
        return None
    extractor = rule.get("extractor") or {}
    kind = extractor.get("kind")
    if kind == "path":
        return _finite_number(_get_path(record, str(extractor.get("path") or "")))
    if kind == "property_path":
        return _finite_number(
            _get_path(
                _property_record(record),
                str(extractor.get("path") or ""),
            )
        )
    if kind == "property_alert_count":
        return _property_alert_count(record, extractor)
    if kind == "property_precursors_in_stock_complete":
        synthesis = _get_path(_property_record(record), "synthesis.metrics")
        if not isinstance(synthesis, dict):
            return None
        precursor_count = _finite_number(synthesis.get("precursor_count"))
        in_stock_count = _finite_number(
            synthesis.get("precursors_in_stock_count")
        )
        if precursor_count is None or in_stock_count is None:
            return None
        if precursor_count < 0 or in_stock_count < 0:
            return None
        return 1.0 if precursor_count == in_stock_count else 0.0
    if kind == "admet_ai_endpoint":
        endpoint_id = str(extractor.get("endpoint_id") or "")
        value, status = _admet_ai_endpoint_evidence(record, endpoint_id)
        return value if status == "valid" else None
    if kind == "chemfh_endpoint_probability":
        endpoint_id = str(extractor.get("endpoint_id") or "")
        value, status = _chemfh_endpoint_evidence(record, endpoint_id)
        return value if status == "valid" else None
    if kind == "chemfh_substructure_alert_count":
        value, status = _chemfh_substructure_alert_evidence(record)
        return value if status == "valid" else None
    if kind == "admet_ai_count":
        endpoint_ids = [
            str(item) for item in (extractor.get("endpoint_ids") or [])
        ]
        if not endpoint_ids:
            return None
        values: list[float] = []
        for endpoint_id in endpoint_ids:
            value, status = _admet_ai_endpoint_evidence(record, endpoint_id)
            if status != "valid" or value is None:
                return None
            values.append(value)
        threshold = _finite_number(extractor.get("threshold"))
        if threshold is None:
            return None
        comparison = str(extractor.get("comparison") or "")
        if comparison == "lt":
            return float(sum(value < threshold for value in values))
        if comparison == "gt":
            return float(sum(value > threshold for value in values))
        raise ValueError(
            f"unsupported ADMET-AI count comparison: {comparison!r}"
        )
    if kind == "posebusters_check":
        checks = ((record.get("posebusters") or {}).get("checks") or {})
        wanted = canonical_posebusters_check_id(extractor.get("check"))
        for key, value in checks.items():
            if canonical_posebusters_check_id(key) == wanted:
                if isinstance(value, bool):
                    return 1.0 if value else 0.0
                return _finite_number(value)
        return None
    if kind == "posebusters_success_ratio":
        posebusters = record.get("posebusters") or {}
        explicit = _finite_number(posebusters.get("success_ratio"))
        if explicit is not None:
            return explicit
        boolean_values = [
            value for value in (posebusters.get("checks") or {}).values()
            if isinstance(value, bool)
        ]
        if not boolean_values:
            return None
        return sum(boolean_values) / len(boolean_values)
    if kind == "blueprint_verdict":
        verdict = str((record.get("blueprint_match") or {}).get("verdict") or "")
        # ``unresolved`` is the absence of evidence, not evidence of a weak pose:
        # blueprint_match.py returns it when no enabled profiler produced
        # interactions (PLIP/ProLIF missing or crashed) or when a constraint is
        # outside the detector's capabilities.  Scoring it as 1.0 landed it
        # between hard_lower=1 and soft_lower=2, so a profiler outage soft-failed
        # a weight-5 metric and turned every pose into MODIFY. The executable policy removes
        # a missing value from numerator and denominator instead, which ``None``
        # does here -- ``gate_value`` maps it to ``not_applicable``.
        return {"fail": 0.0, "warn": 1.0, "pass": 2.0}.get(verdict)
    if kind == "plip_close_contact_count":
        plip = record.get("plip") or {}
        if plip.get("success") is False:
            return None
        interactions = plip.get("interactions")
        if not isinstance(interactions, list):
            return None
        threshold = float(extractor.get("distance", 1.36))
        count = 0
        for interaction in interactions:
            if not isinstance(interaction, dict):
                continue
            distance = _finite_number(interaction.get("distance"))
            if distance is None:
                distance = _finite_number((interaction.get("geometry") or {}).get("distance"))
            if distance is not None and distance < threshold:
                count += 1
        return float(count)
    if kind == "plip_num_interactions":
        plip = record.get("plip")
        if not isinstance(plip, dict) or plip.get("success") is False:
            return None
        return _finite_number(plip.get("num_interactions"))
    if kind == "prolif_num_interactions":
        prolif = record.get("prolif")
        if (
            not isinstance(prolif, dict)
            or prolif.get("success") is False
            or bool(prolif.get("error"))
        ):
            return None
        return _finite_number(prolif.get("num_interactions"))
    if kind == "blueprint_optional_hits":
        blueprint = record.get("blueprint_match")
        if not isinstance(blueprint, dict):
            return None
        optional_total = _finite_number(blueprint.get("optional_total"))
        if optional_total is None or optional_total <= 0:
            return None
        return _finite_number(blueprint.get("optional_hits"))
    if kind == "pose_score":
        wanted = str(extractor.get("metric_id") or "")
        for score in _pose_scores(record):
            if (
                str(score.get("metric_id") or "") == wanted
                and _score_is_usable_for_metric_policy(score)
            ):
                value = _finite_number(score.get("value"))
                return value * float(extractor.get("scale", 1.0)) if value is not None else None
        return None
    raise ValueError(f"unsupported metric extractor kind: {kind!r}")


def assess_record(
    record: dict[str, Any],
    rules: Iterable[dict[str, Any]] = POSE_PROPERTY_METRIC_RULES,
) -> list[dict[str, Any]]:
    """Return per-metric values and level-1 labels.

    Every entry carries both spellings of the same mark: ``label``/
    ``molecule_label`` in the ``pass``/``flag``/``fail`` vocabulary of the executable policy and
    the legacy ``gate``/``molecule_gate`` synonyms.  ``category_code``/``step_id``
    are attached here so that callers which fold raw ``assess_record`` output
    (without going through ``score_records``) can still group by parent category,
    which is what level 2 requires.

    ``availability_status`` is stamped here so reporting can distinguish
    unplanned metrics from planned metrics that produced no result, and a
    caller that folds raw output never reaches ``score_records``' cohort-aware
    planned matrix.  Single-record inference via ``_metric_plan_hint`` is weaker
    than that matrix -- it answers ``None`` whenever the caller declared no plan
    -- and ``None`` is read as planned, which is the pre-existing behaviour.
    """
    rules = tuple(rules)
    declared_failed = _declared_failed_metric_ids(record, rules)
    assessments: list[dict[str, Any]] = []
    for rule in rules:
        value = extract_metric_value(record, rule)
        label = metric_label(value, rule)
        invalid_raw_value = _invalid_metric_raw_value(record, rule)
        molecule_label = label
        # A protocol/pose-export/route metric describes the *workflow*, not the
        # molecule, so its hard failure must never enter a molecular category
        # count as a ``fail``; it is demoted to a mark and handled by rollback.
        #
        if (
            rule.get("molecule_effect") in NON_MOLECULAR_EFFECTS
            and label == "fail"
        ):
            molecule_label = "flag"
        if (
            label == "fail"
            and rule.get("covalent_review")
            and not _alert_catalog_hard_fail_enabled(record, rule)
        ):
            molecule_label = "flag"
        category_code = _metric_category(str(rule["metric_id"]))
        assessments.append({
            "metric_id": rule["metric_id"],
            "tool": rule["tool"],
            "output_metric": rule["output_metric"],
            "category_code": category_code,
            "step_id": _category_step_id(category_code),
            "value": value,
            "invalid_raw_value": invalid_raw_value,
            "evidence_status": (
                "invalid" if invalid_raw_value is not None
                else "available" if value is not None
                else "not_applicable"
            ),
            "availability_status": (
                "available" if value is not None
                else "not_planned"
                if _metric_plan_hint(
                    record, rule, value, declared_failed=declared_failed
                ) is False
                else "planned_missing"
            ),
            "label": label,
            "molecule_label": molecule_label,
            "gate": LABEL_TO_LEGACY_GATE[label],
            "molecule_gate": LABEL_TO_LEGACY_GATE[molecule_label],
            "molecule_effect": rule.get("molecule_effect", "gate"),
            "gate_participates": bool(rule.get(
                "gate_participates",
                any(
                    rule.get(bound) is not None
                    for bound in (
                        "hard_lower", "soft_lower", "soft_upper", "hard_upper",
                    )
                ),
            )),
            "hard_lower": rule.get("hard_lower"),
            "soft_lower": rule.get("soft_lower"),
            "soft_upper": rule.get("soft_upper"),
            "hard_upper": rule.get("hard_upper"),
            "weight": float(rule.get("weight", 0.0)),
            "effective_weight": float(rule.get("weight", 0.0)),
            "family_id": rule.get("family_id"),
            "family_weight_cap": rule.get("family_weight_cap"),
            "score_included": False,
            "normalized_value": None,
            "weighted_score": 0.0,
            "description": rule.get("description", ""),
        })
    return assessments


def resolve_screening_mode(
    mode: Any = None,
    *,
    context: Any = None,
    default: str = DEFAULT_SCREENING_MODE,
) -> dict[str, Any]:
    """Resolve a screening-mode spec into its ``SCREENING_MODES`` entry.

    ``mode`` accepts ``"V1"``/``"V2"``/``"V3"`` in any case, the English mode
    name, or an already-resolved mode dict.  When it is ``None`` the
    ``context`` hint is consulted (``molecule_generation`` -> V3,
    ``sampled_library`` -> V2, ``user_supplied_full`` -> V1) and
    finally ``default`` applies.
    """
    if isinstance(mode, dict):
        candidate = mode.get("mode_id") or mode.get("mode") or mode.get("name")
    else:
        candidate = mode
    resolved: str | None = None
    if candidate is not None:
        token = str(candidate).strip()
        if token:
            upper = token.upper()
            if upper in SCREENING_MODES:
                resolved = upper
            else:
                lowered = token.lower()
                for mode_id, spec in SCREENING_MODES.items():
                    if lowered == str(spec["name"]).lower():
                        resolved = mode_id
                        break
            if resolved is None:
                raise ValueError(
                    f"unknown screening mode: {candidate!r} "
                    f"(expected one of {sorted(SCREENING_MODES)})"
                )
    if resolved is None and context is not None:
        hint = str(context).strip().lower()
        resolved = SCREENING_MODE_SELECTION_HINTS.get(hint)
    if resolved is None:
        resolved = str(default).strip().upper()
    if resolved not in SCREENING_MODES:
        raise ValueError(f"unknown screening mode: {resolved!r}")
    return dict(SCREENING_MODES[resolved])


_LEGACY_KEY_FOR_LABEL: dict[str, str] = {
    "molecule_label": "molecule_gate",
    "label": "gate",
}


def category_label_counts(
    assessments: Iterable[dict[str, Any]],
    *,
    label_key: str = "molecule_label",
) -> dict[str, dict[str, Any]]:
    """Fold level-1 labels into per-category counts and fail ratios.

    ``ratio_fail = count(fail) / (count(fail) + count(flag) + count(pass))``
    exactly as the executable policy defines it, so ``not_applicable`` metrics leave both the
    numerator and the denominator untouched.

    ``count_assessed`` and ``observed_coverage`` are reporting fields for the
    metrics planned this round. They do not alter the numeric fail-ratio test.
    Unplanned and NA metrics do not manufacture labels.
    """
    legacy_key = _LEGACY_KEY_FOR_LABEL.get(label_key, "molecule_gate")
    folds: dict[str, dict[str, Any]] = {}
    assessed_by_category: dict[str, int] = {}
    unplanned_by_category: dict[str, int] = {}
    for item in assessments:
        # Score-only rows have no gate label by design.  They belong in metric
        # scoring coverage, not in this gate-label coverage denominator; counting
        # them here would make a fully observed category look spuriously missing.
        if item.get("gate_participates") is False:
            continue
        label = item.get(label_key)
        if label is None:
            # Tolerate assessment dicts produced by older callers that only
            # carry the legacy keep/modify/reject spelling.
            label = LEGACY_GATE_TO_LABEL.get(str(item.get(legacy_key)))
        category_code = str(
            item.get("category_code")
            or _metric_category(str(item.get("metric_id") or ""))
        )
        if (
            item.get("availability_status") == "not_planned"
            and label not in _LABEL_SEVERITY
        ):
            # Never planned *and* carrying no mark: outside this run's scope.
            # A ``not_planned`` row that still produced a label is a caller
            # inconsistency; keep the label and count it, since dropping an
            # observed fail would be the more dangerous of the two readings.
            unplanned_by_category[category_code] = (
                unplanned_by_category.get(category_code, 0) + 1
            )
            continue
        assessed_by_category[category_code] = (
            assessed_by_category.get(category_code, 0) + 1
        )
        if label not in _LABEL_SEVERITY:
            continue
        fold = folds.get(category_code)
        if fold is None:
            fold = folds[category_code] = {
                "category_code": category_code,
                "step_id": item.get("step_id") or _category_step_id(category_code),
                "count_pass": 0,
                "count_flag": 0,
                "count_fail": 0,
                "count_labeled": 0,
                "count_assessed": 0,
                "count_not_planned": 0,
                "ratio_fail": None,
                "observed_coverage": None,
                "fail_metric_ids": [],
                "flag_metric_ids": [],
            }
        fold[f"count_{label}"] += 1
        fold["count_labeled"] += 1
        if label in {"fail", "flag"}:
            metric_id = item.get("metric_id")
            if metric_id is not None:
                fold[f"{label}_metric_ids"].append(str(metric_id))
    for category_code, fold in folds.items():
        labeled = int(fold["count_labeled"])
        assessed = int(assessed_by_category.get(category_code, labeled))
        fold["count_assessed"] = assessed
        fold["count_not_planned"] = int(unplanned_by_category.get(category_code, 0))
        fold["ratio_fail"] = (
            float(fold["count_fail"]) / labeled if labeled > 0 else None
        )
        fold["observed_coverage"] = (
            round(labeled / assessed, 6) if assessed > 0 else None
        )
    return folds


def category_fold_trips_mode(fold: dict[str, Any], mode: dict[str, Any]) -> bool:
    """Apply the numeric count OR ratio formula; NA never counts as a fail."""
    if not mode.get("eliminates"):
        return False
    count_fail = int(fold.get("count_fail") or 0)
    if count_fail <= 0:
        return False
    count_threshold = mode.get("fail_count_threshold")
    if count_threshold is not None and count_fail >= int(count_threshold):
        return True
    ratio_threshold = mode.get("fail_ratio_threshold")
    ratio_fail = fold.get("ratio_fail")
    ratio_min_fail = mode.get("ratio_min_fail_count") or 1
    if (
        ratio_threshold is not None
        and ratio_fail is not None
        and count_fail >= int(ratio_min_fail)
        and float(ratio_fail) >= float(ratio_threshold)
    ):
        return True
    return False


def category_fold_passes_majority(fold: dict[str, Any]) -> bool:
    """Whether pass labels strictly outnumber all other observed labels."""
    labeled = int(fold.get("count_labeled") or 0)
    return labeled > 0 and 2 * int(fold.get("count_pass") or 0) > labeled


def fold_molecule_verdict(
    folds: dict[str, dict[str, Any]],
    mode: dict[str, Any],
) -> dict[str, Any]:
    """Turn per-category folds into one molecule verdict.

    * any category trips the mode's count/ratio threshold -> ``reject``
    * otherwise every evaluated category has >50% ``pass`` -> ``keep``
    * otherwise                                            -> ``modify``
    * nothing was labeled at all                           -> ``not_evaluated``
    """
    tripped: list[str] = []
    count_fail = 0
    count_flag = 0
    count_labeled = 0
    all_categories_pass_majority = True
    for category_code in sorted(folds):
        fold = folds[category_code]
        count_fail += int(fold.get("count_fail") or 0)
        count_flag += int(fold.get("count_flag") or 0)
        count_labeled += int(fold.get("count_labeled") or 0)
        if int(fold.get("count_labeled") or 0) and not category_fold_passes_majority(fold):
            all_categories_pass_majority = False
        if category_fold_trips_mode(fold, mode):
            tripped.append(category_code)
    if tripped:
        verdict = "reject"
    elif count_labeled and all_categories_pass_majority:
        verdict = "keep"
    elif count_labeled:
        verdict = "modify"
    else:
        verdict = "not_evaluated"
    return {
        "verdict": verdict,
        "screening_mode": mode.get("mode_id"),
        "tripped_categories": tripped,
        "count_fail": count_fail,
        "count_flag": count_flag,
        "count_labeled": count_labeled,
        "categories": {code: dict(fold) for code, fold in folds.items()},
    }


def molecule_gate_report(
    assessments: Iterable[dict[str, Any]],
    *,
    screening_mode: Any = DEFAULT_SCREENING_MODE,
    label_key: str = "molecule_label",
) -> dict[str, Any]:
    """Run the full level-2 fold and return verdict plus per-category evidence."""
    mode = resolve_screening_mode(screening_mode)
    folds = category_label_counts(assessments, label_key=label_key)
    return fold_molecule_verdict(folds, mode)


def aggregate_molecule_gate(
    assessments: Iterable[dict[str, Any]],
    *,
    screening_mode: Any = DEFAULT_SCREENING_MODE,
) -> str:
    """Molecule verdict from the two-level contract: ``classification == metric_gate``.

    Level 1 labels are folded per parent category and only the category
    count/ratio thresholds of the active screening mode can produce ``reject``.
    A single hard-failed metric never eliminates directly.  The later category
    ratio may still reject on that one fail (for example 1/1 observed labels in
    V2), which is an aggregate decision rather than a level-1 side effect.
    """
    return str(
        molecule_gate_report(
            assessments, screening_mode=screening_mode
        )["verdict"]
    )


def aggregate_pose_gate(
    assessments: Iterable[dict[str, Any]],
    *,
    screening_mode: Any = DEFAULT_SCREENING_MODE,
) -> str:
    """Same fold over raw ``label`` marks, before parent-molecule demotions.

    Used for pose-level roll-ups where protocol/pose-export metrics must still be
    visible as failures rather than being demoted to marks.
    """
    return str(
        molecule_gate_report(
            assessments, screening_mode=screening_mode, label_key="label"
        )["verdict"]
    )



def score_records(
    records: Iterable[dict[str, Any]],
    rules: Iterable[dict[str, Any]] = POSE_PROPERTY_METRIC_RULES,
    *,
    normalization_method: str = NORMALIZATION_METHOD,
    metric_weight_overrides: dict[str, Any] | None = None,
    category_weight_overrides: dict[str, Any] | None = None,
    active_steps: Iterable[str] | None = None,
    screening_mode: Any = DEFAULT_SCREENING_MODE,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Gate then score one comparable cohort using the sbdd hierarchy.

    The returned list preserves input order.  Gating is the two-level contract of
    The executable policy: each metric gets a ``pass``/``flag``/``fail`` label, each parent
    category is folded by ``count(fail)``/``ratio(fail)``, and only a category
    that trips the active ``screening_mode`` produces ``reject``.  Only records
    that are not rejected enter metric percentile pools.  Each metric uses
    descending ``RANK.AVG``; target-range metrics rank a hard/soft trapezoid
    desirability.  Signed child weights produce 0-100 category scores, category
    percentiles feed step and overall scores, and category coverage controls rank
    release -- never elimination.

    ``screening_mode`` selects ``V1`` (mark only), ``V2`` (default: 3 fails or a
    50% fail ratio in one category) or ``V3`` (2 fails or a 30% ratio).  The
    count and ratio alternatives are evaluated independently and neither carries
    a minimum fail count or observed-coverage guard.

    ``active_steps`` names the steps this run is expected to have executed.  A
    step named here with no evidence at all is reported as a structural hole in
    ``missing_active_steps`` and cannot be released as ``formal``.  Within an
    active step, coverage follows the executable policy and divides by the parent weight
    this round *planned*, so a category nobody planned leaves both sides of the
    ratio instead of charging every molecule for it.  The
    catalog-referenced ratio is a pipeline-local audit column published alongside
    as ``catalog_category_weighted_coverage`` for information only.
    When ``active_steps`` is omitted, the caller's per-record
    ``scoring_plan``/``metric_context`` declaration is used; failing that, a step
    counts as active only if it planned at least one metric, which cannot detect
    a step whose rules were never loaded.  Callers that know the intended
    pipeline shape should always pass it explicitly.
    """
    if normalization_method != NORMALIZATION_METHOD:
        raise ValueError(
            f"normalization_method must be {NORMALIZATION_METHOD!r}"
        )
    active_screening_mode = resolve_screening_mode(screening_mode)
    rules = _prepare_scoring_rules(
        tuple(deepcopy(tuple(rules))),
        metric_weight_overrides=metric_weight_overrides,
    )
    validate_metric_rules(rules)
    effective_category_weights = _effective_category_weights(
        rules,
        category_weight_overrides=category_weight_overrides,
    )
    scored = [dict(record) for record in records]
    assessments_by_record = [assess_record(record, rules) for record in scored]
    rule_by_id = {rule["metric_id"]: rule for rule in rules}

    # Resolve the evidence plan before making any molecular decision. The plan
    # describes this experiment, not its survivors: inferring it from an earlier
    # gate creates a circular dependency and can change the gate after ranking.
    planned_matrix, declared_failed_metric_ids = _planned_metric_matrix(
        scored, assessments_by_record, rules,
        [not record_has_structural_reject(record) for record in scored],
    )
    for index, assessments in enumerate(assessments_by_record):
        for assessment in assessments:
            assessment["availability_status"] = (
                "available" if assessment.get("value") is not None
                else "planned_missing"
                if planned_matrix[index][assessment["metric_id"]]
                else "not_planned"
            )

    for assessments in assessments_by_record:
        for assessment in assessments:
            rule = rule_by_id[assessment["metric_id"]]
            assessment.update({
                "category_code": rule["category_code"],
                "step_id": rule["step_id"],
                "base_weight": round(float(rule["base_weight"]), 6),
                "base_signed_weight": round(
                    float(rule["base_signed_weight"]), 6,
                ),
                "signed_weight": round(float(rule["signed_weight"]), 6),
                "effective_signed_weight": round(
                    float(rule["signed_weight"]), 6,
                ),
                "rank_transform": rule["rank_transform"],
                "rank_input": None,
                "desirability": None,
                "rank_avg": None,
                "rank_percentage": None,
            })

    # The two-level fold defines the ranking population before protocol
    # compatibility is evaluated.  A rejected record's protocol must not
    # contaminate the percentile pool that it never enters.
    initial_gate_statuses: list[str] = []
    for record, assessments in zip(scored, assessments_by_record):
        metric_gate = aggregate_molecule_gate(
            assessments, screening_mode=active_screening_mode
        )
        # Only a structural defect may survive a re-score as an existing gate.
        # Historical ``classification``/``gate``/``verdict`` values are metric
        # results from another round (and possibly another V1/V2/V3 mode); they
        # must not remove this record from the current percentile pool.
        initial_gate_statuses.append(
            _combine_gate_statuses(
                _preexisting_record_gate(record),
                metric_gate,
            )
        )

    gate_passing_indices = [
        index
        for index, initial_gate in enumerate(initial_gate_statuses)
        if initial_gate != "reject"
    ]
    protocol_audit = _cohort_docking_protocol_check([
        scored[index] for index in gate_passing_indices
    ])
    local_protocol_checks: dict[int, dict[str, Any]] = {}
    local_record_blocks: list[list[str]] = [
        [] for _record in scored
    ]
    protocol_comparable_indices: list[int] = []
    for record_index in gate_passing_indices:
        local_check = _cohort_docking_protocol_check([scored[record_index]])
        local_protocol_checks[record_index] = local_check
        metric_rank_blocks = [
            assessment["metric_id"]
            for assessment in assessments_by_record[record_index]
            if assessment.get("molecule_effect") in RANKING_BLOCKING_EFFECTS
            and assessment.get("label") == "fail"
        ]
        blockers = [
            *metric_rank_blocks,
            *_record_docking_score_blockers(scored[record_index]),
            *_declared_protocol_blockers(scored[record_index]),
            *local_check["blockers"],
        ]
        local_record_blocks[record_index] = list(dict.fromkeys(blockers))
        if not blockers:
            protocol_comparable_indices.append(record_index)

    protocol_comparison = _cohort_docking_protocol_check([
        scored[index] for index in protocol_comparable_indices
    ])
    mixed_protocol_sets: dict[str, set[str]] = {}
    for local_check in local_protocol_checks.values():
        for metric_id, protocol_ids in local_check["mixed_metric_protocols"].items():
            mixed_protocol_sets.setdefault(metric_id, set()).update(protocol_ids)
    for metric_id, protocol_ids in protocol_comparison[
        "mixed_metric_protocols"
    ].items():
        mixed_protocol_sets.setdefault(metric_id, set()).update(protocol_ids)
    mixed_metric_protocols = {
        metric_id: sorted(protocol_ids)
        for metric_id, protocol_ids in sorted(mixed_protocol_sets.items())
    }
    protocol_check = dict(protocol_audit)
    protocol_check["mixed_metric_protocols"] = mixed_metric_protocols
    protocol_check["blockers"] = [
        f"mixed_docking_protocols:{metric_id}"
        for metric_id in mixed_metric_protocols
    ]
    if protocol_check["missing_protocol_records"]:
        protocol_check["blockers"].append(
            "missing_docking_protocol_identity"
        )
    if protocol_check["incomplete_protocol_records"]:
        protocol_check["blockers"].append(
            "incomplete_docking_protocol_identity"
        )

    record_protocol_blocks: list[list[str]] = []
    comparable_index_set = set(protocol_comparable_indices)
    for record_index, record in enumerate(scored):
        if initial_gate_statuses[record_index] == "reject":
            record_protocol_blocks.append([])
            continue
        blockers = list(local_record_blocks[record_index])
        if (
            protocol_audit["engines"]
            and not _record_docking_engines(record)
        ):
            blockers.append("missing_primary_docking_score")
        if record_index in comparable_index_set:
            local_metric_ids = set(
                local_protocol_checks[record_index]["metric_protocol_ids"]
            )
            blockers.extend(
                f"mixed_docking_protocols:{metric_id}"
                for metric_id in protocol_comparison[
                    "mixed_metric_protocols"
                ]
                if metric_id in local_metric_ids
            )
        record_protocol_blocks.append(list(dict.fromkeys(blockers)))

    scoring_eligible = [
        initial_gate != "reject" and not blockers
        for initial_gate, blockers in zip(
            initial_gate_statuses,
            record_protocol_blocks,
        )
    ]

    planned_any = {
        rule["metric_id"]: any(
            scoring_eligible[index] and row[rule["metric_id"]]
            for index, row in enumerate(planned_matrix)
        )
        for rule in rules
    }
    effective_signed_weights, family_caps = _effective_signed_rule_weights(
        rules,
        planned_any,
    )

    # Which steps this run is accountable for.  An explicit declaration is
    # authoritative because a step whose rules were never loaded is invisible to
    # any data-driven inference.
    declared_active_steps = _declared_active_steps(scored, active_steps)
    if declared_active_steps is None:
        active_step_source = "inferred_from_planned_metrics"
        active_step_ids = {
            str(rule["step_id"])
            for rule in rules
            if planned_any[rule["metric_id"]]
            and str(rule["step_id"]) in STEP_IDS
        }
    else:
        active_step_source = "declared_by_caller"
        active_step_ids = set(declared_active_steps)

    reference_category_weights = _reference_category_weights(
        rules,
        effective_category_weights,
        category_weight_overrides=category_weight_overrides,
    )
    reference_step_categories: dict[str, list[str]] = {
        step_id: [
            code
            for code, weight in reference_category_weights.items()
            if weight > 0 and _step_of_category(code, rules) == step_id
        ]
        for step_id in STEP_IDS
    }
    reference_step_parent_weight: dict[str, float] = {
        step_id: sum(
            reference_category_weights[code]
            for code in reference_step_categories[step_id]
        )
        for step_id in STEP_IDS
    }

    normalization: dict[str, dict[str, Any]] = {}
    for rule_index, rule in enumerate(rules):
        metric_id = rule["metric_id"]
        pool = [
            (index, float(assessments_by_record[index][rule_index]["value"]))
            for index in range(len(scored))
            if scoring_eligible[index]
            and assessments_by_record[index][rule_index]["value"] is not None
        ]
        raw_values = [value for _index, value in pool]
        rank_inputs = [
            (
                index,
                _rank_input(
                    value,
                    rule,
                    cohort_min=min(raw_values) if raw_values else None,
                    cohort_max=max(raw_values) if raw_values else None,
                ),
            )
            for index, value in pool
        ]
        metric_percentiles = _rank_avg_percentiles(rank_inputs)
        for index, rank_input in rank_inputs:
            assessment = assessments_by_record[index][rule_index]
            rank_avg, percentile = metric_percentiles[index]
            assessment["rank_input"] = round(rank_input, 6)
            if rule["rank_transform"] == "trapezoid":
                assessment["desirability"] = round(rank_input, 6)
            assessment["rank_avg"] = round(rank_avg, 6)
            # Rank percentages are computational inputs to every higher level
            # of the hierarchy.  Keep their full precision here; presentation
            # layers may round them without changing category membership or
            # manufacturing a tie in a later RANK.AVG pool.
            assessment["rank_percentage"] = percentile
            assessment["normalized_value"] = percentile

        observed_count = len(pool)
        eligible_count = sum(scoring_eligible)
        planned_count = sum(row[metric_id] for row in planned_matrix)
        effective_signed_weight = effective_signed_weights[metric_id]
        normalization[metric_id] = {
            "active": bool(pool),
            "method": normalization_method,
            "rank_transform": rule["rank_transform"],
            "rank_input_min": min((value for _i, value in rank_inputs), default=None),
            "rank_input_max": max((value for _i, value in rank_inputs), default=None),
            # Deprecated aliases retained for one schema transition.
            "oriented_utility_min": min(
                (value for _i, value in rank_inputs), default=None,
            ),
            "oriented_utility_max": max(
                (value for _i, value in rank_inputs), default=None,
            ),
            "weight": float(rule.get("weight", 0.0)),
            "base_weight": round(float(rule["base_weight"]), 6),
            "base_signed_weight": round(
                float(rule["base_signed_weight"]), 6,
            ),
            "signed_weight": round(float(rule["signed_weight"]), 6),
            "effective_weight": round(abs(effective_signed_weight), 6),
            "effective_signed_weight": round(effective_signed_weight, 6),
            "observed_count": observed_count,
            "planned_count": planned_count,
            "cohort_size": len(scored),
            "eligible_count": eligible_count,
            "coverage_fraction": (
                round(observed_count / planned_count, 6)
                if planned_count > 0 else None
            ),
            "complete_cohort": observed_count == eligible_count,
            "score_comparable": bool(pool),
        }

    categories = list(dict.fromkeys(rule["category_code"] for rule in rules))
    for record_index, (record, assessments) in enumerate(
        zip(scored, assessments_by_record)
    ):
        category_scores: dict[str, dict[str, Any]] = {}
        for assessment in assessments:
            metric_id = assessment["metric_id"]
            planned = planned_matrix[record_index][metric_id]
            value_available = assessment.get("value") is not None
            assessment["availability_status"] = (
                "available" if value_available
                else "planned_missing" if planned
                else "not_planned"
            )
            signed_weight = effective_signed_weights[metric_id]
            assessment["effective_weight"] = round(abs(signed_weight), 6)
            assessment["effective_signed_weight"] = round(signed_weight, 6)
            assessment["score_included"] = bool(
                scoring_eligible[record_index]
                and assessment.get("rank_percentage") is not None
                and signed_weight != 0.0
            )
            assessment["weighted_score"] = (
                round(
                    signed_weight * float(assessment["rank_percentage"]),
                    6,
                )
                if assessment["score_included"] else 0.0
            )

        for category_code in categories:
            members = [
                assessment
                for assessment in assessments
                if assessment["category_code"] == category_code
            ]
            planned_members = [
                assessment
                for assessment in members
                if planned_matrix[record_index][assessment["metric_id"]]
            ]
            raw_available_members = [
                assessment
                for assessment in members
                if assessment.get("value") is not None
            ]
            available = [
                assessment
                for assessment in members
                if assessment.get("rank_percentage") is not None
                and effective_signed_weights[assessment["metric_id"]] != 0.0
            ]
            catalog_available = [
                assessment
                for assessment in members
                if assessment.get("value") is not None
                and abs(
                    float(
                        rule_by_id[assessment["metric_id"]][
                            "base_signed_weight"
                        ]
                    )
                ) > 0.0
            ]
            planned = [
                assessment
                for assessment in planned_members
                if effective_signed_weights[assessment["metric_id"]] != 0.0
            ]
            planned_available = [
                assessment
                for assessment in planned
                if assessment.get("value") is not None
            ]
            raw_score = sum(
                effective_signed_weights[assessment["metric_id"]]
                * float(assessment["rank_percentage"])
                for assessment in available
            )
            available_abs_weight = sum(
                abs(effective_signed_weights[assessment["metric_id"]])
                for assessment in available
            )
            planned_abs_weight = sum(
                abs(effective_signed_weights[assessment["metric_id"]])
                for assessment in planned
            )
            planned_available_abs_weight = sum(
                abs(effective_signed_weights[assessment["metric_id"]])
                for assessment in planned_available
            )
            negative_offset = sum(
                max(-effective_signed_weights[assessment["metric_id"]], 0.0)
                for assessment in available
            )
            category_score = (
                100.0 * (raw_score + negative_offset) / available_abs_weight
                if scoring_eligible[record_index] and available_abs_weight > 0
                else None
            )
            catalog_abs_weight = _category_catalog_abs_weight(
                category_code,
                rules,
            )
            implemented_abs_weight = _category_implemented_abs_weight(
                category_code,
                rules,
            )
            implementation_fraction = (
                min(1.0, implemented_abs_weight / catalog_abs_weight)
                if catalog_abs_weight > 0 else None
            )
            catalog_available_abs_weight = sum(
                abs(
                    float(
                        rule_by_id[assessment["metric_id"]][
                            "base_signed_weight"
                        ]
                    )
                )
                for assessment in catalog_available
            )
            category_planned = any(
                rule_by_id[item["metric_id"]].get(
                    "activates_category", True,
                )
                for item in planned_members
            )
            # A metric with zero effective weight and no bound can neither move
            # a score nor gate a molecule.  Missing such a metric is an audit
            # note, not evidence that changes the decision, so it is tracked
            # separately from the missing metrics that do carry consequences.
            missing_members = [
                item for item in planned_members
                if item.get("availability_status") == "planned_missing"
            ]
            consequential_missing = [
                item for item in missing_members
                if _rule_is_consequential(
                    rule_by_id[item["metric_id"]],
                    effective_signed_weights[item["metric_id"]],
                )
            ]
            audit_only_missing = [
                item for item in missing_members
                if item not in consequential_missing
            ]
            # Level 2 of the gate contract.  ``molecule_label`` is
            # the demoted mark, so protocol/pose-export failures are visible in
            # ``label_fold_raw`` but cannot make a molecular category trip.
            molecule_fold = category_label_counts(members).get(
                category_code,
                {
                    "category_code": category_code,
                    "step_id": _category_step_id(category_code, rules),
                    "count_pass": 0,
                    "count_flag": 0,
                    "count_fail": 0,
                    "count_labeled": 0,
                    "ratio_fail": None,
                    "fail_metric_ids": [],
                    "flag_metric_ids": [],
                },
            )
            raw_fold = category_label_counts(members, label_key="label").get(
                category_code, dict(molecule_fold),
            )
            category_verdict = (
                "reject"
                if category_fold_trips_mode(molecule_fold, active_screening_mode)
                else "keep"
                if category_fold_passes_majority(molecule_fold)
                else "modify"
                if molecule_fold["count_labeled"]
                else "not_evaluated"
            )
            category_scores[category_code] = {
                "step_id": _category_step_id(category_code, rules),
                # Worst raw label any member metric tripped, over every member
                # that produced a value.  Seven categories (ALERT/BENCH/IDQC/
                # TOPO/CHEM/DQC/FP) carry no child weight at all: they exist
                # purely to gate, so their score is legitimately NA and this
                # verdict is the only thing they have to report.  Without it a
                # published table cannot tell "this QC ran and passed" from
                # "this QC never ran".
                "gate_status": aggregate_pose_gate(
                    raw_available_members,
                    screening_mode=active_screening_mode,
                ),
                # The category fold itself, published so a reader can audit why
                # the molecule verdict came out the way it did.
                "screening_mode": active_screening_mode["mode_id"],
                "label_counts": {
                    "pass": int(molecule_fold["count_pass"]),
                    "flag": int(molecule_fold["count_flag"]),
                    "fail": int(molecule_fold["count_fail"]),
                    "labeled": int(molecule_fold["count_labeled"]),
                },
                "count_fail": int(molecule_fold["count_fail"]),
                "ratio_fail": (
                    round(float(molecule_fold["ratio_fail"]), 6)
                    if molecule_fold["ratio_fail"] is not None else None
                ),
                "fail_metric_ids": list(molecule_fold["fail_metric_ids"]),
                "flag_metric_ids": list(molecule_fold["flag_metric_ids"]),
                "raw_count_fail": int(raw_fold.get("count_fail") or 0),
                "raw_ratio_fail": (
                    round(float(raw_fold["ratio_fail"]), 6)
                    if raw_fold.get("ratio_fail") is not None else None
                ),
                "trips_screening_mode": category_fold_trips_mode(
                    molecule_fold, active_screening_mode,
                ),
                "category_verdict": category_verdict,
                # This value is reused by StepScore and OverallScore, so it is
                # machine precision rather than a display-rounded weight.
                "parent_weight": effective_category_weights[category_code],
                "planned": category_planned,
                "raw_score": round(raw_score, 6) if available else None,
                "negative_offset": round(negative_offset, 6),
                "available_abs_weight": round(available_abs_weight, 6),
                "planned_abs_weight": round(planned_abs_weight, 6),
                "planned_available_abs_weight": round(
                    planned_available_abs_weight, 6,
                ),
                "catalog_abs_weight": round(catalog_abs_weight, 6),
                "catalog_available_abs_weight": round(
                    catalog_available_abs_weight, 6,
                ),
                # How much of the category's declared child weight this build
                # can compute at all.  Parent weights follow the spec, so a
                # thinly implemented category still votes at full strength; that
                # has to be legible in the report.
                "catalog_implemented_abs_weight": round(
                    implemented_abs_weight, 6,
                ),
                "catalog_implementation_fraction": (
                    round(implementation_fraction, 6)
                    if implementation_fraction is not None else None
                ),
                "evidence_strength": _evidence_strength(implementation_fraction),
                # CategoryScore is ranked below.  Rounding it here can turn
                # genuinely distinct candidates into a false RANK.AVG tie.
                "score": category_score,
                "rank_avg": None,
                "rank_percentage": None,
                "planned_weighted_coverage": (
                    planned_available_abs_weight / planned_abs_weight
                    if planned_abs_weight > 0 else None
                ),
                "catalog_weighted_coverage": (
                    round(
                        min(
                            1.0,
                            catalog_available_abs_weight / catalog_abs_weight,
                        ),
                        6,
                    )
                    if catalog_abs_weight > 0 else None
                ),
                "planned_metric_ids": [
                    item["metric_id"] for item in planned_members
                ],
                "planned_weighted_metric_ids": [
                    item["metric_id"] for item in planned
                ],
                "available_metric_ids": [item["metric_id"] for item in available],
                "raw_available_metric_ids": [
                    item["metric_id"] for item in raw_available_members
                ],
                "catalog_available_metric_ids": [
                    item["metric_id"] for item in catalog_available
                ],
                "missing_metric_ids": [
                    item["metric_id"] for item in missing_members
                ],
                "consequential_missing_metric_ids": [
                    item["metric_id"] for item in consequential_missing
                ],
                "audit_only_missing_metric_ids": [
                    item["metric_id"] for item in audit_only_missing
                ],
                "unusable_for_ranking_metric_ids": [
                    item["metric_id"]
                    for item in planned_members
                    if item.get("value") is not None
                    and item.get("rank_percentage") is None
                ],
            }
        record["metric_assessments"] = assessments
        record["category_scores"] = category_scores

    for category_code in categories:
        pool = [
            (index, float(record["category_scores"][category_code]["score"]))
            for index, record in enumerate(scored)
            if record["category_scores"][category_code]["score"] is not None
        ]
        percentiles = _rank_avg_percentiles(pool)
        for index, _score in pool:
            rank_avg, percentile = percentiles[index]
            category = scored[index]["category_scores"][category_code]
            category["rank_avg"] = round(rank_avg, 6)
            # RP_ic feeds both StepScore and OverallScore and must not be
            # rounded before either calculation.
            category["rank_percentage"] = percentile

    for record_index, record in enumerate(scored):
        step_scores: dict[str, dict[str, Any]] = {}
        for step_id in STEP_IDS:
            step_categories = [
                (code, data)
                for code, data in record["category_scores"].items()
                if data["step_id"] == step_id and data["parent_weight"] > 0
            ]
            step_active = step_id in active_step_ids
            planned_parent_weight = sum(
                data["parent_weight"]
                for _code, data in step_categories
                if data["planned"]
            )
            # The executable policy admits only categories this run planned
            # and whose Score_ic is not NA. A category
            # nobody planned this round -- e.g. a stray PAINS value present while
            # only ``rdkit.qed`` was planned -- must leave numerator and
            # denominator together; counting it in the numerator
            # alone let coverage exceed 1.0 and let an unplanned category vote on
            # StepScore/OverallScore.
            available_parent_weight = sum(
                data["parent_weight"]
                for _code, data in step_categories
                if data["rank_percentage"] is not None and data["planned"]
            )
            # The catalog denominator is informational only: the spec defines
            # coverage against the *planned* weight, so a category
            # nobody planned may not penalize across profiles.  It is still
            # published so a reader
            # can see how much of the full catalog this profile covers, but it
            # never decides a release band or a gate.
            reference_parent_weight = (
                reference_step_parent_weight[step_id] if step_active else 0.0
            )
            missing_reference_categories = (
                [
                    code
                    for code in reference_step_categories[step_id]
                    if record["category_scores"].get(code, {}).get(
                        "rank_percentage"
                    ) is None
                ]
                if step_active else []
            )
            step_score = (
                100.0 * sum(
                    data["parent_weight"] * float(data["rank_percentage"])
                    for _code, data in step_categories
                    if data["rank_percentage"] is not None and data["planned"]
                ) / available_parent_weight
                if scoring_eligible[record_index] and available_parent_weight > 0
                else None
            )
            step_scores[step_id] = {
                # The coordinator uses Step5Score to choose a representative
                # pose, so retain the same full-precision ordering contract.
                "score": step_score,
                "rank_avg": None,
                "active": step_active,
                "planned_parent_weight": round(planned_parent_weight, 6),
                "available_parent_weight": round(available_parent_weight, 6),
                "reference_parent_weight": round(reference_parent_weight, 6),
                # The executable policy defines category_weighted_coverage
                # against the planned parent weight and includes only planned
                # categories whose Score_ic is not NA. An unplanned category leaves
                # numerator and denominator together.
                "category_weighted_coverage": (
                    round(available_parent_weight / planned_parent_weight, 6)
                    if planned_parent_weight > 0 else None
                ),
                # Same ratio under its previous name, kept because
                # publish_layout maps the spec's published column onto it.
                "planned_category_weighted_coverage": (
                    round(available_parent_weight / planned_parent_weight, 6)
                    if planned_parent_weight > 0 else None
                ),
                # Informational pipeline-local catalog coverage.  Says how much
                # of the published catalog this profile reaches; never gates.
                "catalog_category_weighted_coverage": (
                    round(available_parent_weight / reference_parent_weight, 6)
                    if reference_parent_weight > 0 else None
                ),
                "planned_categories": [
                    code for code, data in step_categories
                    if data["planned"]
                ],
                "available_categories": [
                    code for code, data in step_categories
                    if data["rank_percentage"] is not None and data["planned"]
                ],
                "reference_categories": list(
                    reference_step_categories[step_id]
                ) if step_active else [],
                "missing_reference_categories": missing_reference_categories,
            }

        parent_categories = [
            data
            for data in record["category_scores"].values()
            if data["parent_weight"] > 0
        ]
        planned_parent_weight = sum(
            data["parent_weight"]
            for data in parent_categories
            if data["planned"] and data["step_id"] in active_step_ids
        )
        available_parent_weight = sum(
            data["parent_weight"]
            for data in parent_categories
            if data["rank_percentage"] is not None
            and data["planned"]
            and data["step_id"] in active_step_ids
        )
        reference_parent_weight = sum(
            reference_step_parent_weight[step_id] for step_id in active_step_ids
        )
        # Release bands use coverage relative to the current round's plan.
        # The catalog ratio is
        # published beside it but must not penalize, because a
        # category nobody planned this round is identical for every molecule and
        # says nothing about this candidate's evidence.
        category_coverage = (
            available_parent_weight / planned_parent_weight
            if planned_parent_weight > 0 else None
        )
        catalog_category_coverage = (
            available_parent_weight / reference_parent_weight
            if reference_parent_weight > 0 else None
        )
        # An active step that produced no usable category at all is a structural
        # hole, not a partial one: coverage alone can stay above the release
        # threshold while an entire stage's evidence is absent.
        missing_active_steps = [
            step_id for step_id in STEP_IDS
            if step_scores[step_id]["active"]
            and step_scores[step_id]["score"] is None
        ] if scoring_eligible[record_index] else []
        raw_overall_score = (
            100.0 * sum(
                data["parent_weight"] * float(data["rank_percentage"])
                for data in parent_categories
                if data["rank_percentage"] is not None
                and data["planned"]
                and data["step_id"] in active_step_ids
            ) / available_parent_weight
            if scoring_eligible[record_index] and available_parent_weight > 0
            else None
        )

        if initial_gate_statuses[record_index] == "reject":
            ranking_status = "ineligible_gate"
            overall_score = None
        elif record_protocol_blocks[record_index]:
            ranking_status = "blocked_protocol"
            overall_score = None
        elif category_coverage is None:
            ranking_status = "not_evaluated"
            overall_score = None
        elif category_coverage >= FORMAL_COVERAGE_THRESHOLD:
            ranking_status = "formal"
            overall_score = raw_overall_score
        elif category_coverage >= PROVISIONAL_COVERAGE_THRESHOLD:
            ranking_status = "provisional"
            overall_score = raw_overall_score
        else:
            ranking_status = "insufficient_coverage"
            overall_score = None

        # A wholly missing active step can leave overall coverage above 0.50,
        # so the release status is capped independently.  The score is still
        # published, but never as a formal result.
        if missing_active_steps and ranking_status == "formal":
            ranking_status = "provisional"

        # Pipeline-local planned-coverage floors for PHY/ABS/TOX (see
        # CATEGORY_PLANNED_COVERAGE_FLOOR). Elsewhere, zero available weight
        # makes the category NA and triggers the coverage rule. The release
        # band already covers this, so do not invent a per-category default.
        planned_metric_coverage_incomplete = any(
            data["planned_abs_weight"] > 0
            and data["planned_weighted_coverage"] is not None
            and code in CATEGORY_PLANNED_COVERAGE_FLOOR
            and data["planned_weighted_coverage"]
            < CATEGORY_PLANNED_COVERAGE_FLOOR[code]
            for code, data in record["category_scores"].items()
        )
        planned_metric_evidence_incomplete = any(
            data["consequential_missing_metric_ids"]
            for data in record["category_scores"].values()
        )
        # Weight-zero metrics do not contribute to Score_ic, but the executable policy
        # keeps bounded binary/identity/protocol fields as gates or statuses.
        # If such a gate was planned and produced no value, treating coverage as
        # KEEP would erase the very evidence needed to evaluate the gate.  Keep
        # this separate from ordinary non-zero-weight partial coverage: those
        # metrics remain governed by the explicit PHY/ABS/TOX and overall floors.
        missing_planned_zero_weight_gate_ids = sorted({
            assessment["metric_id"]
            for assessment in assessments_by_record[record_index]
            if assessment.get("availability_status") == "planned_missing"
            and float(assessment.get("effective_signed_weight") or 0.0) == 0.0
            and any(
                assessment.get(bound) is not None
                for bound in (
                    "hard_lower", "soft_lower", "soft_upper", "hard_upper",
                )
            )
            and _planned_zero_weight_gate_requires_evidence(
                scored[record_index],
                rule_by_id[assessment["metric_id"]],
            )
        })
        planned_zero_weight_gate_evidence_incomplete = bool(
            missing_planned_zero_weight_gate_ids
        )
        if ranking_status in {
            "ineligible_gate", "blocked_protocol", "not_evaluated",
        }:
            # Gate and protocol outcomes are orthogonal.  A protocol block can
            # suppress every percentile even when all raw evidence is present;
            # do not relabel that as missing molecular evidence.
            coverage_gate = "not_evaluated"
        elif (
            # Apply the release band to the planned ratio for
            # to the categories this round actually planned.
            ranking_status in {"provisional", "insufficient_coverage"}
            or planned_metric_coverage_incomplete
            or planned_zero_weight_gate_evidence_incomplete
            or missing_active_steps
        ):
            # Coverage repair uses the three per-category floors above and
            # the global release band. A single missing planned
            # metric inside an otherwise well-covered category is reported
            # (planned_metric_evidence_incomplete, planned_but_failed_metric_ids,
            # the failure registry) but is not a threshold — treating it as one
            # meant, say, 31 of 32 toxicity endpoints demoting every candidate.
            # Coverage cannot cause rejection, so ``modify`` is the
            # hardest verdict this branch may ever produce.
            coverage_gate = "modify"
        else:
            coverage_gate = "keep"
        molecule_fold_report = molecule_gate_report(
            assessments_by_record[record_index],
            screening_mode=active_screening_mode,
        )
        metric_gate_precoverage = str(molecule_fold_report["verdict"])
        # Incomplete planned coverage is a soft MODIFY protocol condition,
        # while the release thresholds below independently decide whether an
        # overall score/rank may be published.
        # Coverage controls rank release and backtracking only;
        # it can never change the molecule's keep/modify/reject conclusion.
        # ``metric_gate`` is therefore the current two-level metric fold, while
        # ``coverage_gate`` remains an independent protocol/release status.
        structural_reject = record_has_structural_reject(record)
        metric_gate = (
            "reject" if structural_reject else metric_gate_precoverage
        )
        record.update({
            "scoring_schema_version": SCORING_SCHEMA_VERSION,
            "screening_mode": active_screening_mode["mode_id"],
            "screening_mode_name": active_screening_mode["name"],
            "screening_mode_eliminates": bool(
                active_screening_mode["eliminates"]
            ),
            "metric_label_counts": {
                "fail": int(molecule_fold_report["count_fail"]),
                "flag": int(molecule_fold_report["count_flag"]),
                "labeled": int(molecule_fold_report["count_labeled"]),
            },
            "screening_tripped_categories": list(
                molecule_fold_report["tripped_categories"]
            ),
            "metric_gate_precoverage": metric_gate_precoverage,
            "coverage_gate": coverage_gate,
            "coverage_status": ranking_status,
            "metric_gate": metric_gate,
            # The executable policy -- classification IS the two-level metric gate.
            "classification": metric_gate,
            "structural_reject": bool(structural_reject),
            "gate_status": metric_gate,
            "step_scores": step_scores,
            # OverallRank must be derived from the unrounded OverallScore.
            "overall_score": overall_score,
            "overall_rank": None,
            "rank_status": ranking_status,
            "category_weighted_coverage": (
                round(category_coverage, 6)
                if category_coverage is not None else None
            ),
            "catalog_category_weighted_coverage": (
                round(catalog_category_coverage, 6)
                if catalog_category_coverage is not None else None
            ),
            "reference_parent_weight": round(reference_parent_weight, 6),
            "available_parent_weight": round(available_parent_weight, 6),
            "planned_parent_weight": round(planned_parent_weight, 6),
            "planned_category_weighted_coverage": (
                round(category_coverage, 6)
                if category_coverage is not None else None
            ),
            "active_steps": sorted(active_step_ids),
            "active_steps_source": active_step_source,
            "missing_active_steps": missing_active_steps,
            "ranking_blocked": bool(record_protocol_blocks[record_index]),
            "ranking_blocked_by": record_protocol_blocks[record_index],
            "planned_metric_coverage_incomplete": (
                planned_metric_coverage_incomplete
            ),
            "planned_metric_evidence_incomplete": (
                planned_metric_evidence_incomplete
            ),
            "planned_zero_weight_gate_evidence_incomplete": (
                planned_zero_weight_gate_evidence_incomplete
            ),
            "missing_planned_zero_weight_gate_ids": (
                missing_planned_zero_weight_gate_ids
            ),
            # Metrics the caller declared as planned-but-failed.  They stay in
            # Σ_planned, so this list is the audit trail for why
            # coverage dropped when no row anywhere holds a value.
            "planned_but_failed_metric_ids": (
                declared_failed_metric_ids[record_index]
            ),
        })
        # Deprecated flat-score aliases retained for report consumers during
        # the schema transition.  Their meaning is now OverallScore/Rank.
        record["weighted_total_score"] = record["overall_score"]
        record["weighted_score_max"] = (
            100.0 if record["overall_score"] is not None else None
        )
        record["weighted_score_fraction"] = (
            round(record["overall_score"] / 100.0, 6)
            if record["overall_score"] is not None else None
        )
        record["weighted_score_rank"] = None
        record["score_coverage"] = record["category_weighted_coverage"]

    # Step ranks are local descriptive ranks, not another percentile layer in
    # OverallScore.  Rank the full-precision StepScore values only after every
    # record in this already-partitioned cohort/pool has been scored.  Rejected
    # or protocol-blocked records cannot enter a step pool even if a caller
    # supplied a stale score field.
    for step_id in STEP_IDS:
        step_pool: list[tuple[int, float]] = []
        for index, record in enumerate(scored):
            step_score = _finite_number(
                record["step_scores"][step_id]["score"],
            )
            if scoring_eligible[index] and step_score is not None:
                step_pool.append((index, step_score))
        step_ranks = _rank_avg_percentiles(step_pool)
        for index, _score in step_pool:
            rank_avg, _percentile = step_ranks[index]
            scored[index]["step_scores"][step_id]["rank_avg"] = rank_avg

    # Stable flat names are part of the final report contract in sbdd.  Keep
    # the nested structure as the machine-rich representation and expose these
    # aliases so CSV/report consumers do not have to decode JSON objects.
    for record in scored:
        for step_id in STEP_IDS:
            step = record["step_scores"][step_id]
            record[f"{step_id}_score"] = step.get("score")
            record[f"{step_id}_rank"] = step.get("rank_avg")

    overall_pool = [
        (index, float(record["overall_score"]))
        for index, record in enumerate(scored)
        if record["overall_score"] is not None
        and record["rank_status"] in {"formal", "provisional"}
        and record["gate_status"] != "reject"
        and not record["ranking_blocked"]
    ]
    overall_ranks = _rank_avg_percentiles(overall_pool)
    for index, _score in overall_pool:
        rank_avg, _percentile = overall_ranks[index]
        scored[index]["overall_rank"] = round(rank_avg, 6)
        scored[index]["weighted_score_rank"] = round(rank_avg, 6)

    active_weight = sum(
        abs(effective_signed_weights[rule["metric_id"]])
        for rule in rules
        if planned_any[rule["metric_id"]]
    )
    suggested_active_weight = sum(
        abs(float(rule["base_signed_weight"]))
        for rule in rules
        if planned_any[rule["metric_id"]]
    )
    eligible_count = sum(scoring_eligible)
    incomplete_weighted_metrics = [
        rule["metric_id"]
        for rule in rules
        if planned_any[rule["metric_id"]]
        and abs(effective_signed_weights[rule["metric_id"]]) > 0
        and normalization[rule["metric_id"]]["observed_count"] < eligible_count
    ]

    return scored, {
        "schema_version": SCORING_SCHEMA_VERSION,
        "metric_rules_version": METRIC_RULES_VERSION,
        "normalization_method": normalization_method,
        "normalization_scope": "all_evaluated_molecules_in_round",
        "active_weight": round(active_weight, 6),
        "suggested_active_weight": round(suggested_active_weight, 6),
        "incomplete_weighted_metrics": incomplete_weighted_metrics,
        "family_weight_caps": family_caps,
        "category_config": {
            code: {
                "step_id": _category_step_id(code, rules),
                "base_parent_weight": _category_base_parent_weight(code, rules),
                "effective_parent_weight": effective_category_weights[code],
                "catalog_abs_weight": _category_catalog_abs_weight(code, rules),
                "catalog_implemented_abs_weight": (
                    _category_implemented_abs_weight(code, rules)
                ),
                "evidence_strength": _evidence_strength(
                    min(
                        1.0,
                        _category_implemented_abs_weight(code, rules)
                        / _category_catalog_abs_weight(code, rules),
                    )
                    if _category_catalog_abs_weight(code, rules) > 0 else None
                ),
            }
            for code in categories
        },
        # Categories the spec assigns to an active step that this build has no
        # executable rule for.  They stay in the coverage denominator, which is
        # why the release status can never silently ignore them.
        "unimplemented_reference_categories": sorted(
            code
            for code, weight in reference_category_weights.items()
            if weight > 0
            and code not in categories
            and _step_of_category(code, rules) in active_step_ids
        ),
        "active_steps": sorted(active_step_ids),
        "active_steps_source": active_step_source,
        "reference_step_parent_weight": {
            step_id: round(reference_step_parent_weight[step_id], 6)
            for step_id in STEP_IDS
        },
        "coverage_release_thresholds": {
            "formal": FORMAL_COVERAGE_THRESHOLD,
            "provisional": PROVISIONAL_COVERAGE_THRESHOLD,
        },
        # The executable policy -- coverage only chooses a release band; it can never reject
        # a molecule or a pose.
        "coverage_can_reject": False,
        "screening_mode": dict(active_screening_mode),
        "screening_modes": {
            mode_id: dict(spec) for mode_id, spec in SCREENING_MODES.items()
        },
        "gate_contract": "two_level_label_then_category_fold",
        "metric_labels": list(METRIC_LABELS),
        "metric_weight_overrides": dict(metric_weight_overrides or {}),
        "category_weight_overrides": dict(category_weight_overrides or {}),
        "protocol_compatible": not protocol_check["blockers"],
        "protocol_blockers": protocol_check["blockers"],
        "docking_engines": protocol_check["engines"],
        "docking_protocol_ids": protocol_check["protocol_ids"],
        "docking_metric_protocol_ids": protocol_check["metric_protocol_ids"],
        "mixed_docking_metric_protocols": protocol_check[
            "mixed_metric_protocols"
        ],
        "docking_protocol_missing_records": protocol_check[
            "missing_protocol_records"
        ],
        "docking_protocol_incomplete_records": protocol_check[
            "incomplete_protocol_records"
        ],
        "normalization": normalization,
    }


def _prepare_scoring_rules(
    rules: tuple[dict[str, Any], ...],
    *,
    metric_weight_overrides: dict[str, Any] | None,
) -> tuple[dict[str, Any], ...]:
    overrides = dict(metric_weight_overrides or {})
    known_metric_ids = {str(rule.get("metric_id") or "") for rule in rules}
    unknown = sorted(set(overrides) - known_metric_ids)
    if unknown:
        raise ValueError(f"unknown metric weight override(s): {', '.join(unknown)}")

    prepared: list[dict[str, Any]] = []
    for rule in rules:
        metric_id = str(rule.get("metric_id") or "")
        direction = _rule_direction(rule)
        category_code = str(rule.get("category_code") or "UNCLASSIFIED")
        step_id = str(
            rule.get("step_id")
            or {
                "property_filter": "step4",
                "docking": "step5",
            }.get(str(rule.get("stage") or ""), "step6")
        )
        base_weight = abs(float(rule.get("weight", 0.0)))
        magnitude = (
            _runtime_weight(overrides[metric_id], f"metric {metric_id}")
            if metric_id in overrides else base_weight
        )
        explicit_signed = _finite_number(rule.get("signed_weight"))
        if explicit_signed is not None and explicit_signed != 0:
            sign = -1.0 if explicit_signed < 0 else 1.0
        else:
            sign = -1.0 if direction == "lower" else 1.0
        candidate = dict(rule)
        candidate.update({
            "direction": direction,
            "base_weight": base_weight,
            "base_signed_weight": (
                0.0 if base_weight == 0 else sign * base_weight
            ),
            "weight": magnitude,
            "signed_weight": 0.0 if magnitude == 0 else sign * magnitude,
            "category_code": category_code,
            "step_id": step_id,
            "category_weight": float(
                rule.get(
                    "category_weight",
                    CATEGORY_CONFIG.get(category_code, {}).get(
                        "parent_weight", 1.0,
                    ),
                )
            ),
            "catalog_abs_weight": float(
                rule.get(
                    "catalog_abs_weight",
                    CATEGORY_CONFIG.get(category_code, {}).get(
                        "catalog_abs_weight", 0.0,
                    ),
                )
            ),
            "rank_transform": str(
                rule.get("rank_transform")
                or ("trapezoid" if direction == "target_range" else "identity")
            ),
            "plan_tool": rule.get("plan_tool") or _plan_tool(
                metric_id,
                rule.get("requires_property_tool"),
            ),
        })
        if candidate["rank_transform"] not in {"identity", "trapezoid"}:
            raise ValueError(
                f"{metric_id}: rank_transform must be identity or trapezoid"
            )
        prepared.append(candidate)
    return tuple(prepared)


def _rule_direction(rule: dict[str, Any]) -> str | None:
    direction = rule.get("direction")
    if direction in {"higher", "lower", "target_range"}:
        return str(direction)
    has_lower = any(
        _finite_number(rule.get(field)) is not None
        for field in ("hard_lower", "soft_lower")
    )
    has_upper = any(
        _finite_number(rule.get(field)) is not None
        for field in ("soft_upper", "hard_upper")
    )
    if has_lower and has_upper:
        return "target_range"
    if has_lower:
        return "higher"
    if has_upper:
        return "lower"
    return "higher"


def _runtime_weight(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} weight override must be between 0 and 5")
    number = _finite_number(value)
    if number is None or number < 0 or number > 5:
        raise ValueError(f"{label} weight override must be between 0 and 5")
    return float(number)


def _reference_category_weights(
    rules: tuple[dict[str, Any], ...],
    effective_category_weights: dict[str, float],
    *,
    category_weight_overrides: dict[str, Any] | None,
) -> dict[str, float]:
    """Parent weights for every category the spec expects, not just the coded ones.

    Categories with no executable rule keep their published §3 parent weight so
    they stay in the coverage denominator.  A profile that genuinely drops one
    still can, by overriding its weight to zero.
    """
    overrides = dict(category_weight_overrides or {})
    result: dict[str, float] = {}
    for category_code, config in CATEGORY_CONFIG.items():
        if category_code in effective_category_weights:
            result[category_code] = float(
                effective_category_weights[category_code]
            )
        elif category_code in overrides:
            result[category_code] = _runtime_weight(
                overrides[category_code],
                f"category {category_code}",
            )
        else:
            result[category_code] = float(config["parent_weight"])
    for category_code, weight in effective_category_weights.items():
        result.setdefault(category_code, float(weight))
    return result


def _step_of_category(
    category_code: str,
    rules: tuple[dict[str, Any], ...],
) -> str:
    return _category_step_id(category_code, rules)


def _effective_category_weights(
    rules: tuple[dict[str, Any], ...],
    *,
    category_weight_overrides: dict[str, Any] | None,
) -> dict[str, float]:
    categories = list(dict.fromkeys(rule["category_code"] for rule in rules))
    overrides = dict(category_weight_overrides or {})
    unknown = sorted(set(overrides) - set(CATEGORY_CONFIG) - set(categories))
    if unknown:
        raise ValueError(
            f"unknown category weight override(s): {', '.join(unknown)}"
        )
    result: dict[str, float] = {}
    for category_code in categories:
        base = _category_base_parent_weight(category_code, rules)
        result[category_code] = (
            _runtime_weight(
                overrides[category_code],
                f"category {category_code}",
            )
            if category_code in overrides else base
        )
    return result


def record_has_structural_reject(record: dict[str, Any] | None) -> bool:
    """Return whether ``record`` has an sbdd structural hard failure.

    The field ``nomination_eligibility=rejected`` is intentionally *not* enough
    to establish this.  The property filter historically used that value for a
    category-fold metric rejection, which must be recomputed under the current
    screening mode (and must disappear entirely under V1).  New producers set
    ``structural_reject`` explicitly; the evidence/reason fallbacks keep legacy
    structure-evaluator records safe.
    """
    if not isinstance(record, dict):
        return False
    if record.get("structural_reject") is True:
        return True
    if record.get("valid_structure") is False:
        return True
    for field in (
        "identity_reject",
        "identity_failed",
        "identity_qc_failed",
        "structural_failure",
    ):
        if record.get(field) is True:
            return True

    # Legacy/flat pose reports exposed the identity checks as scalar fields
    # rather than an ``identity_status`` object.  These checks are explicitly
    # structural hard failures in the executable policy (identity/valence cannot be repaired
    # by a V1 marking mode), so preserve that meaning when re-scoring them.
    for field in (
        "canonical_smiles_consistent",
        "inchikey_consistent",
        "graph_isomorphic",
        "stereochemistry_consistent",
        "valence_check",
        "formal_charge_check",
        "tautomer_consistent",
        "protomer_charge_consistent",
    ):
        if record.get(field) is False:
            return True

    identity_qc = record.get("identity_qc")
    if identity_qc is False:
        return True
    if isinstance(identity_qc, dict):
        status = str(identity_qc.get("status") or "").strip().casefold()
        if status in {
            "failed", "fail", "invalid", "error", "mismatch", "rejected",
        }:
            return True

    identity_status = record.get("identity_status")
    if isinstance(identity_status, dict):
        status = str(identity_status.get("status") or "").strip().casefold()
        if status in {
            "failed", "fail", "invalid", "error", "mismatch", "rejected",
        }:
            return True

    evidence_status = record.get("evidence_status")
    if isinstance(evidence_status, dict):
        status = str(evidence_status.get("identity_qc") or "").strip().casefold()
        if status in {"failed", "fail", "error", "mismatch", "rejected"}:
            return True

    verdict_reasons = record.get("verdict_reasons")
    if isinstance(verdict_reasons, dict) and verdict_reasons.get("hard"):
        return True

    # A producer may provide only a reason string for a hard structural stop.
    reason = str(record.get("structural_reject_reason") or "").strip()
    if reason:
        return True

    # Do not infer a structural failure from ``nomination_eligibility`` alone.
    # Older property-filter reports used ``rejected`` for a metric category
    # fold, and the executable policy requires that threshold result to be recomputed
    # (and capped under V1).  Producers that mean a real parse/identity failure
    # must emit one of the explicit signals above.
    return False


def _preexisting_record_gate(record: dict[str, Any]) -> str:
    """Return only a structural gate carried in from an upstream record.

    The executable policy requires the current round to fold its own level-1 labels.
    Reusing old ``classification``, ``intrinsic_classification``, ``gate`` or
    ``verdict`` here made a V2/V3 rejection survive a V1 re-score and made a
    historical reject suppress a newly passing record.  Coverage and protocol
    review statuses are likewise orthogonal to the molecule gate.
    """
    return "reject" if record_has_structural_reject(record) else "keep"


def _normalize_gate_status(value: Any) -> str | None:
    text = str(value or "").strip().casefold()
    return {
        "keep": "keep", "pass": "keep", "passed": "keep",
        "modify": "modify", "warn": "modify", "warning": "modify",
        "flag": "modify", "review": "modify",
        "reject": "reject", "fail": "reject", "failed": "reject",
        "not_evaluated": "not_evaluated", "not_applicable": "not_evaluated",
        "unknown": "not_evaluated", "n/a": "not_evaluated",
    }.get(text)


def _combine_gate_statuses(*statuses: str) -> str:
    severity = {"not_evaluated": -1, "keep": 0, "modify": 1, "reject": 2}
    normalized = [
        status for status in statuses if status in severity
    ]
    if not normalized:
        return "not_evaluated"
    return max(normalized, key=severity.__getitem__)


def _planned_metric_matrix(
    records: list[dict[str, Any]],
    assessments_by_record: list[list[dict[str, Any]]],
    rules: tuple[dict[str, Any], ...],
    scoring_eligible: list[bool],
) -> tuple[list[dict[str, bool]], list[list[str]]]:
    """Decide, per record and metric, whether this run planned to produce it.

    Returns the planned matrix plus, per record, the metric ids the caller
    declared as *planned but failed*.  That second list is the answer to the
    hole the planned denominator would otherwise open: a tool that failed for
    the entire cohort leaves no row holding a value, so cohort inference alone
    would read it as "never planned", drop it from Σ_planned, and report full
    coverage for a total outage. The executable policy uses declarations:
    unselected or inapplicable metrics leave Σ_planned, while planned but failed
    metrics remain in it and enter the failure registry.
    """
    hints: list[dict[str, bool | None]] = []
    declared_failed: list[list[str]] = []
    for record, assessments in zip(records, assessments_by_record):
        failed = _declared_failed_metric_ids(record, rules)
        declared_failed.append(sorted(failed))
        hints.append({
            rule["metric_id"]: _metric_plan_hint(
                record,
                rule,
                assessment.get("value"),
                declared_failed=failed,
            )
            for rule, assessment in zip(rules, assessments)
        })

    planned_default = {
        rule["metric_id"]: any(
            scoring_eligible[index]
            and (
                assessments_by_record[index][rule_index].get("value") is not None
                or hints[index][rule["metric_id"]] is True
            )
            for index in range(len(records))
        )
        for rule_index, rule in enumerate(rules)
    }
    matrix = [
        {
            rule["metric_id"]: (
                False
                if hints[record_index][rule["metric_id"]] is False
                else True
                if assessments_by_record[record_index][rule_index].get("value")
                is not None
                else True
                if hints[record_index][rule["metric_id"]] is True
                else planned_default[rule["metric_id"]]
            )
            for rule_index, rule in enumerate(rules)
        }
        for record_index in range(len(records))
    ]
    return matrix, declared_failed


def _metric_plan_hint(
    record: dict[str, Any],
    rule: dict[str, Any],
    value: Any,
    *,
    declared_failed: set[str] | None = None,
) -> bool | None:
    metric_id = str(rule["metric_id"])
    category_code = str(rule.get("category_code") or _metric_category(metric_id))
    contexts = _metric_plan_contexts(record, rule)

    if _metric_explicitly_not_applicable(record, rule):
        return False
    if declared_failed and metric_id in declared_failed:
        # A tool the caller reports as having failed was planned, whatever the
        # rest of the cohort looks like.  Saying so here is what
        # keeps a whole-cohort outage inside Σ_planned instead of erasing it.
        return True
    for context in contexts:
        required_ids: set[str] = set()
        for field in ("required_metric_ids", "required_admet_metric_ids"):
            values = context.get(field)
            if isinstance(values, (list, tuple, set)):
                required_ids.update(str(item).strip() for item in values)
        bare_metric_id = metric_id.split(".", 1)[1] if "." in metric_id else metric_id
        if metric_id in required_ids or bare_metric_id in required_ids:
            return True
        planned_ids = context.get("planned_metric_ids")
        if isinstance(planned_ids, (list, tuple, set)):
            declared = {str(item).strip() for item in planned_ids}
            return metric_id in declared or bare_metric_id in declared
    if value is not None:
        return True

    for context in contexts:
        planned_categories = context.get("planned_category_codes")
        if isinstance(planned_categories, (list, tuple, set)) and category_code in {
            str(item) for item in planned_categories
        }:
            return True

    plan_tool = rule.get("plan_tool")
    if plan_tool:
        for context in contexts:
            enabled = context.get("enabled_tools")
            if isinstance(enabled, (list, tuple, set)):
                return str(plan_tool) in {str(item) for item in enabled}

    if rule.get("stage") == "docking":
        engines = _record_docking_engines(record)
        if engines and plan_tool:
            return str(plan_tool) in engines
    return None


def _metric_plan_contexts(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> list[dict[str, Any]]:
    # The coordinator's outer metric_context comes from pose evaluation.
    # Property metrics must keep their producer's plan when that row is nested;
    # an absent tool in another stage does not mean it was disabled here.
    # Explicit whole-record scoring_plan overrides retain first priority.
    stage_context = (
        _property_metric_context(record)
        if rule.get("stage") == "property_filter"
        else record.get("metric_context")
    )
    return [
        context for context in (
            record.get("scoring_plan"),
            stage_context,
        )
        if isinstance(context, dict)
    ]


def _planned_zero_weight_gate_requires_evidence(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> bool:
    """Keep optional ChemFH outages from masquerading as required gates.

    The executable policy makes an empty alert backend a backtracking trigger when
    that tool is *required*.  Older rows lacked a required/optional distinction, so
    they retain the conservative behavior; an explicit modern contract is
    authoritative.
    """
    if rule.get("requires_property_tool") != "chemfh":
        return True
    contexts = _metric_plan_contexts(record, rule)
    metric_id = str(rule.get("metric_id") or "")
    required_ids: set[str] = set()
    required_tools: set[str] = set()
    explicit = False
    for context in contexts:
        if any(
            field in context
            for field in (
                "required_metric_ids", "required_admet_metric_ids",
                "required_tools", "optional_tools",
            )
        ):
            explicit = True
        values = context.get("required_metric_ids")
        if isinstance(values, (list, tuple, set)):
            for item in values:
                declared = str(item).strip()
                known_bare = {
                    endpoint_id for endpoint_id, _name in _CHEMFH_ENDPOINTS
                } | {"substructure_alert_count"}
                if declared in known_bare:
                    declared = f"chemfh.{declared}"
                if declared:
                    required_ids.add(declared)
        values = context.get("required_tools")
        if isinstance(values, (list, tuple, set)):
            required_tools.update(
                str(item).strip().casefold() for item in values
            )
    if not explicit:
        return True
    return metric_id in required_ids or "chemfh" in required_tools


def _metric_explicitly_not_applicable(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> bool:
    metric_id = str(rule["metric_id"])
    bare_metric_id = metric_id.split(".", 1)[1] if "." in metric_id else metric_id
    context_declared = any(
        isinstance(context.get("not_applicable_metric_ids"), (list, tuple, set))
        and bool(
            {metric_id, bare_metric_id}.intersection({
                str(item).strip()
                for item in context["not_applicable_metric_ids"]
            })
        )
        for context in _metric_plan_contexts(record, rule)
    )
    if context_declared:
        return True
    if not metric_id.startswith("blueprint_match."):
        return False
    blueprint = record.get("blueprint_match")
    if not isinstance(blueprint, dict):
        return False
    if blueprint.get("blueprint_applicable") is False:
        return True
    if blueprint.get("blueprint_valid") is False:
        return False
    # A selected matcher does not plan every kind of interaction constraint.
    # Use declared scope/counts, not a missing value or an unresolved gate:
    # profiler failures must still remain planned-but-missing evidence.
    if metric_id in {"blueprint_match.verdict", "blueprint_match.required_hit_rate"}:
        if _finite_number(blueprint.get("required_total")) == 0:
            return True
        if metric_id == "blueprint_match.verdict":
            return (
                blueprint.get("required_hard_gate_declared") is False
                or blueprint.get("hard_gate_declared") is False
            )
    if metric_id.startswith("blueprint_match.desired_"):
        return (
            _finite_number(blueprint.get("desired_total")) == 0
            or blueprint.get("desired_metric_applicable") is False
        )
    if metric_id == "blueprint_match.optional_hits":
        return _finite_number(blueprint.get("optional_total")) == 0
    return False


def _blueprint_metric_applicable(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> bool:
    """Apply Blueprint validity without conflating Required and Desired gates."""
    metric_id = str(rule.get("metric_id") or "")
    if not metric_id.startswith("blueprint_match."):
        return True
    blueprint = record.get("blueprint_match")
    if not isinstance(blueprint, dict):
        return True
    if (
        blueprint.get("blueprint_applicable") is False
        or blueprint.get("blueprint_valid") is False
    ):
        return False
    if metric_id == "blueprint_match.verdict":
        return True
    if metric_id.startswith("blueprint_match.desired_"):
        if metric_id == "blueprint_match.desired_ranking_score":
            explicit_ranking = blueprint.get(
                "desired_ranking_gate_applicable"
            )
            if explicit_ranking is not None:
                return bool(explicit_ranking)
        if metric_id == "blueprint_match.desired_evaluable_coverage":
            explicit_coverage = blueprint.get(
                "desired_coverage_protocol_gate_applicable"
            )
            if explicit_coverage is not None:
                return bool(explicit_coverage)
        explicit = blueprint.get("desired_metric_applicable")
        if explicit is not None:
            return bool(explicit)
        total = _finite_number(blueprint.get("desired_total"))
        return total is None or total > 0
    if metric_id == "blueprint_match.optional_hits":
        explicit = blueprint.get("optional_metric_applicable")
        if explicit is not None:
            return bool(explicit)
    return True


def _tool_alias_key(value: Any) -> str:
    """Fold the several spellings one tool has across reports into one key.

    A rule carries both a plan name (``admet``) and a display name
    (``ADMET-AI``); upstream ``operational_errors`` entries use either, and
    sometimes a step-qualified form.  Comparing on letters and digits alone is
    enough to join them without a hand-maintained alias table.
    """
    return "".join(
        character
        for character in str(value or "").casefold()
        if character.isalnum()
    )


def _declared_failure_tool_keys(record: dict[str, Any]) -> set[str]:
    """Tool names the caller reported as having run and failed.

    Only failures are collected.  A tool the caller chose not to run is a plan
    decision and belongs in ``enabled_tools``/``not_applicable_metric_ids``,
    which already remove it from Σ_planned.
    """
    keys: set[str] = set()
    for context in (
        record.get("scoring_plan"),
        record.get("metric_context"),
        _property_metric_context(record),
    ):
        if not isinstance(context, dict):
            continue
        for field in ("failed_tools", "planned_but_failed_tools"):
            value = context.get(field)
            if isinstance(value, (list, tuple, set)):
                keys.update(_tool_alias_key(item) for item in value)
        errors = context.get("operational_errors")
        if isinstance(errors, (list, tuple)):
            for entry in errors:
                if isinstance(entry, dict):
                    keys.add(_tool_alias_key(entry.get("tool")))
                else:
                    keys.add(_tool_alias_key(entry))
    keys.discard("")
    return keys


def _declared_failed_metric_ids(
    record: dict[str, Any],
    rules: tuple[dict[str, Any], ...],
) -> set[str]:
    """Metric ids this record's caller declared as planned but not delivered."""
    declared: set[str] = set()
    for context in (
        record.get("scoring_plan"),
        record.get("metric_context"),
        _property_metric_context(record),
    ):
        if not isinstance(context, dict):
            continue
        for field in ("failed_metric_ids", "planned_but_failed_metric_ids"):
            value = context.get(field)
            if isinstance(value, (list, tuple, set)):
                declared.update(str(item) for item in value)

    tool_keys = _declared_failure_tool_keys(record)
    if tool_keys:
        for rule in rules:
            candidates = {
                _tool_alias_key(rule.get("plan_tool")),
                _tool_alias_key(rule.get("tool")),
            }
            candidates.discard("")
            if candidates & tool_keys:
                declared.add(str(rule["metric_id"]))

    known = {str(rule["metric_id"]) for rule in rules}
    return {metric_id for metric_id in declared if metric_id in known}


def _effective_signed_rule_weights(
    rules: tuple[dict[str, Any], ...],
    planned_any: dict[str, bool],
) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
    effective = {
        rule["metric_id"]: float(rule["signed_weight"])
        for rule in rules
    }
    grouped: dict[str, list[dict[str, Any]]] = {}
    for rule in rules:
        family_id = rule.get("family_id")
        if (
            family_id
            and planned_any[rule["metric_id"]]
            and abs(float(rule["signed_weight"])) > 0
        ):
            grouped.setdefault(str(family_id), []).append(rule)

    metadata: dict[str, dict[str, Any]] = {}
    for family_id, members in grouped.items():
        caps = {float(rule["family_weight_cap"]) for rule in members}
        if len(caps) != 1:
            raise ValueError(
                f"{family_id}: active rules must declare one family_weight_cap"
            )
        cap = caps.pop()
        declared_total = sum(abs(float(rule["signed_weight"])) for rule in members)
        scale = min(1.0, cap / declared_total) if declared_total > 0 else 1.0
        for rule in members:
            metric_id = rule["metric_id"]
            effective[metric_id] = float(rule["signed_weight"]) * scale
        metadata[family_id] = {
            "cap": cap,
            "active_declared_abs_weight": round(declared_total, 6),
            "active_effective_abs_weight": round(
                sum(abs(effective[rule["metric_id"]]) for rule in members), 6,
            ),
            "scale": round(scale, 6),
            "active_metrics": [rule["metric_id"] for rule in members],
        }
    return effective, metadata


def _rank_input(
    value: float,
    rule: dict[str, Any],
    *,
    cohort_min: float | None,
    cohort_max: float | None,
) -> float:
    if rule.get("rank_transform") != "trapezoid":
        return float(value)
    hard_lower = _finite_number(rule.get("hard_lower"))
    soft_lower = _finite_number(rule.get("soft_lower"))
    soft_upper = _finite_number(rule.get("soft_upper"))
    hard_upper = _finite_number(rule.get("hard_upper"))
    lower = hard_lower if hard_lower is not None else cohort_min
    upper = hard_upper if hard_upper is not None else cohort_max
    if lower is None or upper is None:
        return 1.0
    soft_lower = lower if soft_lower is None else soft_lower
    soft_upper = upper if soft_upper is None else soft_upper
    if value < lower or value > upper:
        return 0.0
    if value < soft_lower:
        span = soft_lower - lower
        return 1.0 if span <= 0 else _clamp01((value - lower) / span)
    if value <= soft_upper:
        return 1.0
    span = upper - soft_upper
    return 1.0 if span <= 0 else _clamp01((upper - value) / span)


def _rank_avg_percentiles(
    indexed_values: Iterable[tuple[int, float]],
) -> dict[int, tuple[float, float]]:
    values = sorted(
        ((int(index), float(value)) for index, value in indexed_values),
        key=lambda item: item[1],
        reverse=True,
    )
    count = len(values)
    result: dict[int, tuple[float, float]] = {}
    cursor = 0
    while cursor < count:
        end = cursor + 1
        while end < count and values[end][1] == values[cursor][1]:
            end += 1
        first_rank = cursor + 1
        last_rank = end
        average_rank = (first_rank + last_rank) / 2.0
        percentile = (
            0.5 if count == 1
            else (count - average_rank) / (count - 1)
        )
        for index, _value in values[cursor:end]:
            result[index] = (average_rank, _clamp01(percentile))
        cursor = end
    return result


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def _category_base_parent_weight(
    category_code: str,
    rules: tuple[dict[str, Any], ...],
) -> float:
    if category_code in CATEGORY_CONFIG:
        return float(CATEGORY_CONFIG[category_code]["parent_weight"])
    for rule in rules:
        if rule["category_code"] == category_code:
            return float(rule.get("category_weight", 1.0))
    return 1.0


def _category_step_id(
    category_code: str,
    rules: tuple[dict[str, Any], ...] = (),
) -> str:
    if category_code in CATEGORY_CONFIG:
        return str(CATEGORY_CONFIG[category_code]["step_id"])
    for rule in rules:
        if rule["category_code"] == category_code:
            return str(rule.get("step_id") or "step6")
    return "step6"


def _category_catalog_abs_weight(
    category_code: str,
    rules: tuple[dict[str, Any], ...],
) -> float:
    configured = CATEGORY_CONFIG.get(category_code, {}).get("catalog_abs_weight")
    if configured is not None:
        return float(configured)
    return sum(
        abs(float(rule.get("base_signed_weight", rule["signed_weight"])))
        for rule in rules
        if rule["category_code"] == category_code
    )


def _category_implemented_abs_weight(
    category_code: str,
    rules: tuple[dict[str, Any], ...],
) -> float:
    """Sum of catalog child weights this build can actually compute.

    Compared with ``_category_catalog_abs_weight`` this exposes how much of the
    spec's evidence for a category is backed by an executable extractor.  A
    category may keep its full parent weight while resting on a small slice of
    its declared child weight; that asymmetry must be visible in the report
    rather than inferred from the rule table.
    """
    return sum(
        abs(float(rule.get("base_signed_weight", rule["signed_weight"])))
        for rule in rules
        if rule["category_code"] == category_code
    )


def _rule_is_consequential(
    rule: dict[str, Any],
    effective_signed_weight: float,
) -> bool:
    """Return whether a missing value for this rule can change any outcome.

    A rule with zero effective weight and no bound is pure provenance: it cannot
    shift a category score and it cannot gate.  Treating its absence as missing
    evidence would demote KEEP to MODIFY for a field nothing depends on.
    """
    if float(effective_signed_weight) != 0.0:
        return True
    return any(
        rule.get(bound) is not None
        for bound in ("hard_lower", "soft_lower", "soft_upper", "hard_upper")
    )


def _evidence_strength(fraction: float | None) -> str:
    """Grade how much of a category's declared child weight is implemented."""
    if fraction is None:
        return "not_applicable"
    if fraction <= 0:
        return "absent"
    if fraction < 0.5:
        return "sparse"
    if fraction < 0.9:
        return "partial"
    return "full"


def catalog_step_parent_weight(step_id: str) -> float:
    """Static ``Σ W'_c`` for one step from configured parent-category weights.

    Coverage must not be measured against whatever the current build happens to
    implement, otherwise a category that is missing from the rule catalog leaves
    both the numerator and the denominator and the ratio reads 1.0.
    """
    return sum(
        float(config["parent_weight"])
        for config in CATEGORY_CONFIG.values()
        if config["step_id"] == step_id and float(config["parent_weight"]) > 0
    )


def _declared_active_steps(
    records: list[dict[str, Any]],
    declared: Iterable[str] | None,
) -> set[str] | None:
    """Read the caller's explicit statement of which steps this run executes."""
    if declared is not None:
        selected = {
            str(step_id) for step_id in declared if str(step_id) in STEP_IDS
        }
        return selected
    from_records: set[str] = set()
    found = False
    for record in records:
        for container_key in ("scoring_plan", "metric_context"):
            container = record.get(container_key)
            if not isinstance(container, dict):
                continue
            values = container.get("active_steps")
            if isinstance(values, (list, tuple, set)):
                found = True
                from_records.update(
                    str(step_id) for step_id in values
                    if str(step_id) in STEP_IDS
                )
    return from_records if found else None


def public_metric_rules(
    rules: Iterable[dict[str, Any]] = POSE_PROPERTY_METRIC_RULES,
) -> list[dict[str, Any]]:
    """Return catalog/report-safe rule metadata without extraction internals."""
    fields = (
        "metric_id", "stage", "tool", "output_metric", "hard_lower",
        "soft_lower", "soft_upper", "hard_upper", "weight", "direction",
        "signed_weight", "rank_transform", "category_code", "step_id",
        "category_weight", "catalog_abs_weight", "plan_tool",
        "molecule_effect", "gate_participates", "family_id",
        "family_weight_cap", "description",
        "requires_property_tool", "covalent_review", "activates_category",
    )
    return [{key: rule.get(key) for key in fields} for rule in rules]


def _property_record(record: dict[str, Any]) -> dict[str, Any]:
    """Return a Step 5 row whether it is direct or nested by coordinator."""
    nested = record.get("property")
    return nested if isinstance(nested, dict) else record


def _admet_ai_endpoint_evidence(
    record: dict[str, Any],
    endpoint_id: str,
) -> tuple[float | None, str]:
    """Return the raw finite endpoint value and its evidence status.

    Status is one of ``valid``, ``invalid``, or ``missing``.  A ``partial``
    backend result does not invalidate its other endpoints: validity remains
    endpoint-local, as required by the missing-aware scoring contract.
    """
    admet = _property_record(record).get("admet")
    if not isinstance(admet, dict):
        return None, "missing"
    predictions = admet.get("predictions")
    if not isinstance(predictions, dict) or endpoint_id not in predictions:
        return None, "missing"

    endpoint = predictions.get(endpoint_id)
    raw_value = endpoint.get("value") if isinstance(endpoint, dict) else endpoint
    # ``bool`` is a subclass of ``int`` in Python.  Older reports without an
    # endpoint status would otherwise turn True/False into 1.0/0.0 and let a
    # malformed payload participate in gates and ranking.
    if isinstance(raw_value, bool):
        return raw_value, "invalid"
    value = _finite_number(raw_value)
    missing = {
        str(item) for item in (admet.get("missing_endpoints") or [])
    }
    invalid = {
        str(item) for item in (admet.get("invalid_endpoints") or [])
    }
    if endpoint_id in missing or value is None:
        return None, "missing"

    if isinstance(endpoint, dict):
        explicit_status = str(endpoint.get("status") or "").strip().lower()
        if endpoint.get("valid") is False or explicit_status in {
            "invalid", "error", "failed", "missing", "unavailable",
            "not_applicable",
        }:
            invalid.add(endpoint_id)

    if endpoint_id in _ADMET_AI_PROBABILITY_ENDPOINT_IDS and not 0 <= value <= 1:
        invalid.add(endpoint_id)
    if endpoint_id in _ADMET_AI_NONNEGATIVE_ENDPOINT_IDS and value < 0:
        invalid.add(endpoint_id)
    if endpoint_id == "PPBR_AZ" and not 0 <= value <= 100:
        invalid.add(endpoint_id)
    return value, "invalid" if endpoint_id in invalid else "valid"


def _invalid_admet_ai_raw_value(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> float | None:
    """Return finite invalid raw evidence for audit, never for gate or rank."""
    extractor = rule.get("extractor") or {}
    if extractor.get("kind") != "admet_ai_endpoint":
        return None
    value, status = _admet_ai_endpoint_evidence(
        record,
        str(extractor.get("endpoint_id") or ""),
    )
    return value if status == "invalid" else None


def _chemfh_endpoint_evidence(
    record: dict[str, Any],
    endpoint_id: str,
) -> tuple[float | None, str]:
    """Return one ChemFH probability only when its local evidence is valid."""
    chemfh = _property_record(record).get("chemfh")
    if not isinstance(chemfh, dict):
        return None, "missing"
    status = str(chemfh.get("status") or "").strip().casefold()
    if status not in {"complete", "partial"}:
        return None, "missing"
    predictions = chemfh.get("predictions")
    if not isinstance(predictions, dict):
        return None, "missing"
    endpoint = predictions.get(endpoint_id)
    if not isinstance(endpoint, dict):
        return None, "missing"

    raw_probability = endpoint.get("probability")
    raw_uncertainty = endpoint.get("uncertainty")
    if isinstance(raw_probability, bool) or isinstance(raw_uncertainty, bool):
        return None, "invalid"
    probability = _finite_number(raw_probability)
    uncertainty = _finite_number(raw_uncertainty)
    if probability is not None and not 0 <= probability <= 1:
        return probability, "invalid"
    if uncertainty is not None and uncertainty < 0:
        return probability, "invalid"
    missing = {str(item) for item in (chemfh.get("missing_endpoints") or [])}
    if endpoint_id in missing or probability is None or uncertainty is None:
        return None, "missing"
    positive = endpoint.get("positive")
    if not isinstance(positive, bool) or positive != (probability > 0.5):
        return probability, "invalid"
    return probability, "valid"


def _chemfh_substructure_alert_evidence(
    record: dict[str, Any],
) -> tuple[float | None, str]:
    """Count hit catalogs across the ten-rule sbdd2 mapping, never matches."""
    chemfh = _property_record(record).get("chemfh")
    if not isinstance(chemfh, dict):
        return None, "missing"
    status = str(chemfh.get("status") or "").strip().casefold()
    if status not in {"complete", "partial"}:
        return None, "missing"
    missing_rules = {
        str(item) for item in (chemfh.get("missing_rules") or [])
    }
    if missing_rules.intersection(_CHEMFH_SUBSTRUCTURE_RULE_NAMES):
        return None, "missing"
    rules = chemfh.get("rules")
    if not isinstance(rules, dict):
        return None, "missing"
    hit_count = 0
    for rule_name in _CHEMFH_SUBSTRUCTURE_RULE_NAMES:
        rule = rules.get(rule_name)
        if not isinstance(rule, dict):
            return None, "missing"
        count = rule.get("match_count")
        atoms = rule.get("atom_indices")
        if (
            isinstance(count, bool)
            or not isinstance(count, int)
            or count < 0
            or not isinstance(atoms, list)
            or len(atoms) != count
        ):
            return None, "invalid"
        hit_count += int(count > 0)
    return float(hit_count), "valid"


def _invalid_metric_raw_value(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> float | None:
    """Expose malformed finite evidence for audit without gating on it."""
    admet_value = _invalid_admet_ai_raw_value(record, rule)
    if admet_value is not None:
        return admet_value
    extractor = rule.get("extractor") or {}
    kind = extractor.get("kind")
    if kind == "chemfh_endpoint_probability":
        value, status = _chemfh_endpoint_evidence(
            record, str(extractor.get("endpoint_id") or ""),
        )
        return value if status == "invalid" else None
    if kind == "chemfh_substructure_alert_count":
        value, status = _chemfh_substructure_alert_evidence(record)
        if status != "invalid":
            return None
        summary = (_property_record(record).get("chemfh") or {}).get("summary")
        if isinstance(summary, dict):
            return _finite_number(summary.get("substructure_alert_count"))
    return None


def _property_metric_context(record: dict[str, Any]) -> dict[str, Any]:
    context = _property_record(record).get("metric_context")
    return context if isinstance(context, dict) else {}


def _required_property_tool_enabled(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> bool:
    required = rule.get("requires_property_tool")
    if required in (None, ""):
        return True
    enabled = _property_metric_context(record).get("enabled_tools")
    if not isinstance(enabled, (list, tuple, set)):
        # Older reports did not publish tool context. A present scalar remains
        # usable for backward compatibility; current reports are explicit.
        return True
    return str(required) in {str(item) for item in enabled}


def _alert_catalog_hard_fail_enabled(
    record: dict[str, Any],
    rule: dict[str, Any],
) -> bool:
    configured = _property_metric_context(record).get("hard_fail_catalogs")
    if not isinstance(configured, (list, tuple, set)):
        return True
    hard_catalogs = {str(item).strip().upper() for item in configured}
    requested = {
        str(item).strip().upper()
        for item in ((rule.get("extractor") or {}).get("catalogs") or [])
    }
    base = _property_record(record)
    alerts = base.get("structural_alerts")
    if isinstance(alerts, dict) and isinstance(alerts.get("matches"), list):
        matched_catalogs = {
            str(item.get("catalog") or "").strip().upper()
            for item in alerts["matches"]
            if isinstance(item, dict)
            and str(item.get("catalog") or "").strip().upper() in requested
        }
        if matched_catalogs:
            return any(
                (
                    catalog in hard_catalogs
                    if catalog != "PAINS"
                    else "PAINS" in hard_catalogs
                    or any(name.startswith("PAINS_") for name in hard_catalogs)
                )
                for catalog in matched_catalogs
            )

    # Older PAINS projections may report only a boolean/count without the
    # originating A/B/C catalog. Fail closed unless the entire requested
    # catalog family was explicitly removed from the hard-fail policy.
    return bool(hard_catalogs & requested)


def _property_alert_count(
    record: dict[str, Any],
    extractor: dict[str, Any],
) -> float | None:
    base = _property_record(record)
    alerts = base.get("structural_alerts")
    requested = {
        str(item).strip().upper()
        for item in (extractor.get("catalogs") or [])
    }
    if isinstance(alerts, dict):
        if alerts.get("disabled") is True or alerts.get("available") is False:
            return None
        matches = alerts.get("matches")
        if isinstance(matches, list):
            return float(sum(
                1
                for item in matches
                if isinstance(item, dict)
                and str(item.get("catalog") or "").strip().upper() in requested
            ))

    if requested & {"PAINS", "PAINS_A", "PAINS_B", "PAINS_C"}:
        pains = base.get("pains")
        if isinstance(pains, dict) and pains.get("is_pains") is not None:
            matches = pains.get("matches")
            if isinstance(matches, list):
                return float(len(matches))
            filters = pains.get("matched_filters")
            if isinstance(filters, list):
                return float(len(filters))
            return 1.0 if pains.get("is_pains") else 0.0
    return None


def _get_path(value: Any, path: str) -> Any:
    current = value
    for part in path.split(".") if path else ():
        if not isinstance(current, dict):
            return None
        current = current.get(part)
    return current


def _pose_scores(record: dict[str, Any]) -> list[dict[str, Any]]:
    result = collect_pose_scores(record)
    primary = record.get("primary_score")
    if isinstance(primary, dict):
        result.extend(collect_pose_scores({"pose_scores": [primary]}))
    return result


def _score_is_usable_for_metric_policy(score: dict[str, Any]) -> bool:
    """Accept finite applicable evidence, including audit-only directions."""
    if score.get("applicability") != "applicable":
        return False
    if str(score.get("reliability") or "").casefold() in {
        "failed", "invalid", "not_evaluated", "unreliable", "unsupported",
    }:
        return False
    return _finite_number(score.get("value")) is not None


def _cohort_docking_protocol_check(
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    """Reject only like-for-like metrics mixed across protocol identities.

    Different engines expose different metric IDs and may coexist in one
    protocol profile. A single metric ID may not be percentile-ranked across
    multiple receptor/grid/version/checkpoint protocols.
    """
    engines: set[str] = set()
    protocol_ids: set[str] = set()
    metric_protocol_ids: dict[str, set[str]] = {}
    missing_protocol_records: list[str] = []
    incomplete_protocol_records: list[str] = []
    for record_index, record in enumerate(records):
        # The same scalar must use the same extension protocol/model across
        # a ranking cohort. Per-ligand coordinates are evidence, not protocol IDs.
        for block_name, prefix in (("rdkit_pose", "rdkit_pose."), ("freesasa", "freesasa."),
                                   ("mdanalysis", "mdanalysis.")):
            block = record.get(block_name) or {}
            protocol = block.get("protocol") or {}
            if protocol and block.get("success") is True:
                signature = json.dumps(protocol, sort_keys=True, default=str)
                for extension_rule in _ADDITIONAL_METRIC_RULES:
                    extension_id = extension_rule["metric_id"]
                    if extension_id.startswith(prefix) and extract_metric_value(record, extension_rule) is not None:
                        metric_protocol_ids.setdefault(extension_id, set()).add(signature)
        rascore = _property_record(record).get("rascore") or {}
        if rascore.get("status") == "available" and rascore.get("model_sha256"):
            metric_protocol_ids.setdefault("rascore.score", set()).add(str(rascore["model_sha256"]))
        for extension_score in _pose_scores(record):
            if (extension_score.get("metric_id") == "rtmscore.ensemble_score"
                    and _score_is_usable_for_metric_policy(extension_score)
                    and extension_score.get("protocol_id")):
                metric_protocol_ids.setdefault("rtmscore.ensemble_score", set()).add(str(extension_score["protocol_id"]))
        docking_scores = [
            score
            for score in _pose_scores(record)
            if _docking_engine_for_score(score) is not None
            and _score_is_usable_for_metric_policy(score)
        ]
        record_engines = {
            engine
            for score in docking_scores
            if (engine := _docking_engine_for_score(score)) is not None
        }
        declared_metadata = _record_docking_metadata(record)
        record_engines.update(declared_metadata["engines"])
        engines.update(record_engines)
        record_has_explicit_protocol = False
        record_has_incomplete_protocol = False
        for score in docking_scores:
            metric_id = str(score.get("metric_id") or "")
            explicit = score.get("protocol_id")
            if explicit not in (None, ""):
                protocol_id = str(explicit)
                protocol_ids.add(protocol_id)
                if metric_id:
                    metric_protocol_ids.setdefault(metric_id, set()).add(
                        protocol_id
                    )
                record_has_explicit_protocol = True
                if score.get("protocol_identity_complete") is False:
                    record_has_incomplete_protocol = True
                continue
            scoring_function = score.get("scoring_function")
            method = score.get("method")
            method_version = score.get("method_version")
            if method not in (None, "") and scoring_function not in (None, ""):
                protocol_id = (
                    f"{method}:{scoring_function}:{method_version or 'unknown'}"
                )
                protocol_ids.add(protocol_id)
                if metric_id:
                    metric_protocol_ids.setdefault(metric_id, set()).add(
                        protocol_id
                    )
        protocol_ids.update(declared_metadata["protocol_ids"])
        for engine in declared_metadata["engines"]:
            for protocol_id in declared_metadata["protocol_ids"]:
                metric_protocol_ids.setdefault(
                    f"engine:{engine}", set()
                ).add(protocol_id)
        if declared_metadata["protocol_ids"]:
            record_has_explicit_protocol = True
        if declared_metadata["incomplete_protocol_identity"]:
            record_has_incomplete_protocol = True
        record_name = str(
            record.get("mol_name")
            or record.get("compound_id")
            or f"record_{record_index}"
        )
        if record_engines and not record_has_explicit_protocol:
            missing_protocol_records.append(record_name)
        if record_has_incomplete_protocol:
            incomplete_protocol_records.append(record_name)
    blockers = []
    mixed_metric_protocols = {
        metric_id: sorted(ids)
        for metric_id, ids in metric_protocol_ids.items()
        if len(ids) > 1
    }
    blockers.extend(
        f"mixed_docking_protocols:{metric_id}"
        for metric_id in sorted(mixed_metric_protocols)
    )
    if missing_protocol_records:
        blockers.append("missing_docking_protocol_identity")
    if incomplete_protocol_records:
        blockers.append("incomplete_docking_protocol_identity")
    return {
        "engines": sorted(engines),
        "protocol_ids": sorted(protocol_ids),
        "metric_protocol_ids": {
            metric_id: sorted(ids)
            for metric_id, ids in sorted(metric_protocol_ids.items())
        },
        "mixed_metric_protocols": mixed_metric_protocols,
        "missing_protocol_records": missing_protocol_records,
        "incomplete_protocol_records": incomplete_protocol_records,
        "blockers": blockers,
    }


def _record_docking_engines(record: dict[str, Any]) -> set[str]:
    engines = {
        engine
        for score in _pose_scores(record)
        if _score_is_usable_for_metric_policy(score)
        and (engine := _docking_engine_for_score(score)) is not None
    }
    engines.update(_record_docking_metadata(record)["engines"])
    return engines


def _record_docking_metadata(record: dict[str, Any]) -> dict[str, Any]:
    """Collect declared docking identity even when the primary score is absent."""
    containers = [record]
    for key in ("docking_scores", "score_protocol"):
        nested = record.get(key)
        if isinstance(nested, dict):
            containers.append(nested)

    engines: set[str] = set()
    protocol_ids: set[str] = set()
    incomplete_protocol_identity = False
    for container in containers:
        engine = _normalize_docking_engine(
            container.get("engine")
            or container.get("docking_engine")
            or container.get("method")
        )
        if engine is None:
            continue
        engines.add(engine)
        protocol_id = container.get("protocol_id")
        if protocol_id not in (None, ""):
            protocol_ids.add(str(protocol_id))
        if container.get("protocol_identity_complete") is False:
            incomplete_protocol_identity = True
    return {
        "engines": engines,
        "protocol_ids": protocol_ids,
        "incomplete_protocol_identity": incomplete_protocol_identity,
    }


def _normalize_docking_engine(value: Any) -> str | None:
    text = str(value or "").strip().casefold()
    aliases = {
        "vina": "autodock_vina",
        "autodock-vina": "autodock_vina",
        "uni-dock": "unidock",
    }
    normalized = aliases.get(text, text)
    return normalized if normalized in {
        "gnina", "unidock", "autodock_vina", "diffdock", "dock6",
        "autodock4", "rdock", "adfr",
    } else None


def _docking_engine_for_score(score: dict[str, Any]) -> str | None:
    """Return the supported native docking engine for one score, if any."""
    metric_id = str(score.get("metric_id") or "").casefold()
    method = _normalize_docking_engine(score.get("method"))
    if metric_id.startswith("gnina.") or method == "gnina":
        return "gnina"
    if metric_id.startswith("unidock.") or method == "unidock":
        return "unidock"
    if metric_id.startswith("autodock_vina.") or method in {
        "autodock_vina", "vina",
    }:
        return "autodock_vina"
    if metric_id.startswith("diffdock.") or method == "diffdock":
        return "diffdock"
    if metric_id.startswith("dock6.") or method == "dock6":
        return "dock6"
    if metric_id.startswith("autodock4.") or method == "autodock4":
        return "autodock4"
    if metric_id.startswith("rdock.") or method == "rdock":
        return "rdock"
    if metric_id.startswith("adfr.") or method == "adfr":
        return "adfr"
    return None


def _record_docking_score_blockers(record: dict[str, Any]) -> list[str]:
    engines = _record_docking_engines(record)
    if not engines:
        return []
    available = {
        str(score.get("metric_id") or "")
        for score in _pose_scores(record)
        if _score_is_usable_for_metric_policy(score)
    }
    required_any = {
        "gnina": {
            "gnina.minimizedAffinity", "gnina.CNNscore", "gnina.CNNaffinity",
        },
        "unidock": {"unidock.docking_score"},
        "autodock_vina": {"autodock_vina.vina_affinity"},
        "diffdock": {"diffdock.confidence"},
        "dock6": {"dock6.grid_score"},
        "autodock4": {"autodock4.binding_energy"},
        "rdock": {"rdock.score"},
        "adfr": {"adfr.affinity"},
    }
    return [
        f"missing_primary_docking_score:{engine}"
        for engine in sorted(engines)
        if not available & required_any[engine]
    ]


def _declared_protocol_blockers(record: dict[str, Any]) -> list[str]:
    """Carry fail-closed provenance/protocol blockers into rank eligibility."""
    blockers: list[str] = []
    for context in (
        record,
        record.get("scoring_plan"),
        record.get("metric_context"),
        _property_metric_context(record),
    ):
        if not isinstance(context, dict):
            continue
        values = context.get("protocol_blockers")
        if isinstance(values, (list, tuple, set)):
            blockers.extend(str(item) for item in values if str(item))
    return list(dict.fromkeys(blockers))


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


validate_metric_rules(ALL_METRIC_RULES)

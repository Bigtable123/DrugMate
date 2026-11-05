"""Shared constants for SBDD skills — single source of truth for thresholds."""
from __future__ import annotations


# --- PoseCheck thresholds ---
# Clash thresholds apply to the evaluator's heavy-atom-only clash count.
# The raw PoseCheck all-atom count is still reported as diagnostics.
DEFAULT_POSECHECK_THRESHOLDS = {
    "strain_warn_kcal": 8.0,
    "strain_fail_kcal": 60.0,
    "clash_warn_count": 0,
    "clash_fail_count": 2,
}

DEFAULT_POSECHECK_NATIVE_BASELINE_LIMITS = {
    "strain_max_kcal": 100.0,
    "clash_max_count": 3,
}

# --- PoseBusters ---
# Legacy PoseBusters summaries without per-check evidence remain a strict Step
# 4 physical-validity gate.  When checks are available, registered metric rules
# distinguish hard failures from modify/review outcomes.
POSEBUSTERS_VIOLATION_THRESHOLD = 0

# --- Blueprint matching ---
# BLUEPRINT_PASS_RATE/BLUEPRINT_WARN_RATE retain their historical values only
# for readers of legacy reports. Schema-v2 required constraints do not use
# percentage gates; all-of/any-of/at-least-k logic is evaluated exactly.
BLUEPRINT_PASS_RATE = 0.6
BLUEPRINT_WARN_RATE = 0.4
BLUEPRINT_DISTANCE_TOLERANCE_ANGSTROM = 0.75
BLUEPRINT_ANGLE_TOLERANCE_DEGREES = 20.0
BLUEPRINT_DESIRED_DEFAULT_WEIGHT = 1.0
BLUEPRINT_CONFIDENCE_WEIGHTS = {
    "high": 1.0,
    "medium": 0.6,
    "low": 0.3,
}

# --- Pocket occupancy ---
POCKET_OCCUPANCY_LOW = 0.2
POCKET_OCCUPANCY_HIGH = 0.7

# --- Docking defaults (GNINA) ---
DEFAULT_DOCKING_PARAMS = {
    "autobox_add": 4.0,
    "exhaustiveness": 16,
    "num_modes": 20,
    "pose_sort_order": "CNNscore",
    "cnn_scoring": "rescore",
    "seed": 42,
}


def pose_quota_for_compound_count(compound_count: int) -> int:
    """Maximum docking poses carried into QC for a cohort of this size."""
    if compound_count <= 50:
        return 20
    if compound_count <= 1000:
        return 5
    return 1


# --- Pipeline step timeout floors (seconds) ---
# These are conservative outer-process floors for large libraries. The actual
# timeout is raised further from the sum of enabled serial tool budgets.
DEFAULT_STEP_TIMEOUTS = {
    "structure_prep": 7200,
    "pocket_analysis": 14400,
    "molecule_generation": 172800,
    "property_prefilter": 43200,
    "docking": 172800,
    "pose_evaluation": 172800,
    "property_filter": 43200,
    "coordination": 7200,
    "reporting": 7200,
}

# Legacy numeric aliases retained for standalone callers.
DEFAULT_TIMEOUTS = {
    0: DEFAULT_STEP_TIMEOUTS["structure_prep"],
    1: DEFAULT_STEP_TIMEOUTS["pocket_analysis"],
    2: DEFAULT_STEP_TIMEOUTS["molecule_generation"],
    "2.5": DEFAULT_STEP_TIMEOUTS["property_prefilter"],
    3: DEFAULT_STEP_TIMEOUTS["docking"],
    4: DEFAULT_STEP_TIMEOUTS["pose_evaluation"],
    5: DEFAULT_STEP_TIMEOUTS["property_filter"],
    6: DEFAULT_STEP_TIMEOUTS["coordination"],
    7: DEFAULT_STEP_TIMEOUTS["reporting"],
}

# --- Interaction type canonical mapping ---
# Canonical form: hyphenated (matches blueprint generation output)
INTERACTION_TYPE_CANONICAL = {
    "hbond": "H-bond",
    "h-bond": "H-bond",
    "hbdonor": "H-bond",
    "hbacceptor": "H-bond",
    "hydrophobic": "Hydrophobic",
    "vdw": "VdW",
    "vdwcontact": "VdW",
    "pistacking": "Pi-stacking",
    "pi-stacking": "Pi-stacking",
    "facetoface": "Pi-stacking",
    "pication": "Pi-cation",
    "cationpi": "Pi-cation",
    "pi-cation": "Pi-cation",
    "saltbridge": "Salt-bridge",
    "salt-bridge": "Salt-bridge",
    "salt_bridge": "Salt-bridge",
    "waterbridge": "Water-bridge",
    "water-bridge": "Water-bridge",
    "water_bridge": "Water-bridge",
    "halogen": "Halogen-bond",
    "halogenbond": "Halogen-bond",
    "xbdonor": "Halogen-bond",
    "xbacceptor": "Halogen-bond",
    "halogen-bond": "Halogen-bond",
    "metalcomplex": "Metal-complex",
    "metaldonor": "Metal-complex",
    "metalacceptor": "Metal-complex",
    "edgetoface": "Edge-to-face",
    "facetoedge": "Edge-to-face",
}


def normalize_interaction_type(itype: str) -> str:
    """Normalize interaction type to canonical hyphenated form."""
    key = itype.lower().replace(" ", "").replace("_", "").replace("-", "")
    # Direct lookup on compacted key
    for k, v in INTERACTION_TYPE_CANONICAL.items():
        if key == k.replace(" ", "").replace("_", "").replace("-", ""):
            return v
    return itype

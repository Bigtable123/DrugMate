"""Strict ADMET-AI endpoint contract shared by worker and Drugmate adapter.

Only these provider-native endpoints are consumed by the registered scorer.
Provider releases may expose additional columns, but those columns are audit
evidence and cannot silently expand or replace the scoring payload.
"""
from __future__ import annotations

import math
from typing import Any


ADMET_NATIVE_ENDPOINTS_BY_CATEGORY: dict[str, tuple[str, ...]] = {
    "ABS": (
        "HIA_Hou",
        "Caco2_Wang",
        "PAMPA_NCATS",
        "Pgp_Broccatelli",
        "Bioavailability_Ma",
        "Solubility_AqSolDB",
        "Lipophilicity_AstraZeneca",
        "HydrationFreeEnergy_FreeSolv",
    ),
    "DIST": ("BBB_Martins", "PPBR_AZ", "VDss_Lombardo"),
    "MET": (
        "CYP2C9_Substrate_CarbonMangels",
        "CYP2D6_Substrate_CarbonMangels",
        "CYP3A4_Substrate_CarbonMangels",
        "CYP1A2_Veith",
        "CYP2C19_Veith",
        "CYP2C9_Veith",
        "CYP2D6_Veith",
        "CYP3A4_Veith",
    ),
    "EXCR": (
        "Half_Life_Obach",
        "Clearance_Hepatocyte_AZ",
        "Clearance_Microsome_AZ",
    ),
    "TOX": (
        "hERG",
        "ClinTox",
        "AMES",
        "DILI",
        "Carcinogens_Lagunin",
        "LD50_Zhu",
        "Skin_Reaction",
        "NR-AR",
        "NR-AR-LBD",
        "NR-AhR",
        "NR-Aromatase",
        "NR-ER",
        "NR-ER-LBD",
        "NR-PPAR-gamma",
        "SR-ARE",
        "SR-ATAD5",
        "SR-HSE",
        "SR-MMP",
        "SR-p53",
    ),
}

ADMET_NATIVE_ENDPOINT_IDS: tuple[str, ...] = tuple(
    endpoint_id
    for category_endpoint_ids in ADMET_NATIVE_ENDPOINTS_BY_CATEGORY.values()
    for endpoint_id in category_endpoint_ids
)

# The scorer derives its two CYP count metrics solely from these native
# probabilities.  Keeping the dependencies explicit makes contract drift
# independently testable even if the rule ordering changes.
ADMET_DERIVED_COUNT_INPUTS: dict[str, tuple[str, ...]] = {
    "cyp_substrate_count": (
        "CYP2C9_Substrate_CarbonMangels",
        "CYP2D6_Substrate_CarbonMangels",
        "CYP3A4_Substrate_CarbonMangels",
    ),
    "cyp_inhibition_count": (
        "CYP1A2_Veith",
        "CYP2C19_Veith",
        "CYP2C9_Veith",
        "CYP2D6_Veith",
        "CYP3A4_Veith",
    ),
}

DOMAIN_FINITE = "finite"
DOMAIN_PROBABILITY = "probability_0_1"
DOMAIN_NONNEGATIVE = "nonnegative"
DOMAIN_PERCENT_0_100 = "percent_0_100"

ADMET_PROBABILITY_ENDPOINT_IDS = frozenset({
    "HIA_Hou",
    "PAMPA_NCATS",
    "Pgp_Broccatelli",
    "Bioavailability_Ma",
    "BBB_Martins",
    "CYP2C9_Substrate_CarbonMangels",
    "CYP2D6_Substrate_CarbonMangels",
    "CYP3A4_Substrate_CarbonMangels",
    "CYP1A2_Veith",
    "CYP2C19_Veith",
    "CYP2C9_Veith",
    "CYP2D6_Veith",
    "CYP3A4_Veith",
    "hERG",
    "ClinTox",
    "AMES",
    "DILI",
    "Carcinogens_Lagunin",
    "Skin_Reaction",
    "NR-AR",
    "NR-AR-LBD",
    "NR-AhR",
    "NR-Aromatase",
    "NR-ER",
    "NR-ER-LBD",
    "NR-PPAR-gamma",
    "SR-ARE",
    "SR-ATAD5",
    "SR-HSE",
    "SR-MMP",
    "SR-p53",
})

# sbdd explicitly marks each EXCR value below zero as nonphysical.  VDss is
# likewise a nonnegative physical quantity and PPBR is a percentage.
ADMET_NONNEGATIVE_ENDPOINT_IDS = frozenset({
    "VDss_Lombardo",
    "Half_Life_Obach",
    "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ",
})

ADMET_ENDPOINT_DOMAINS: dict[str, str] = {
    endpoint_id: DOMAIN_FINITE for endpoint_id in ADMET_NATIVE_ENDPOINT_IDS
}
ADMET_ENDPOINT_DOMAINS.update({
    endpoint_id: DOMAIN_PROBABILITY
    for endpoint_id in ADMET_PROBABILITY_ENDPOINT_IDS
})
ADMET_ENDPOINT_DOMAINS.update({
    endpoint_id: DOMAIN_NONNEGATIVE
    for endpoint_id in ADMET_NONNEGATIVE_ENDPOINT_IDS
})
ADMET_ENDPOINT_DOMAINS["PPBR_AZ"] = DOMAIN_PERCENT_0_100


def finite_number(value: Any) -> float | None:
    """Return a finite numeric value, rejecting booleans and non-numbers."""
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def endpoint_value_error(endpoint_id: str, value: Any) -> str | None:
    """Return the endpoint-domain violation, or ``None`` for valid evidence."""
    if endpoint_id not in ADMET_ENDPOINT_DOMAINS:
        return "endpoint is not registered"
    number = finite_number(value)
    if number is None:
        return "value is not finite numeric evidence"

    domain = ADMET_ENDPOINT_DOMAINS[endpoint_id]
    if domain == DOMAIN_PROBABILITY and not 0.0 <= number <= 1.0:
        return "probability is outside [0, 1]"
    if domain == DOMAIN_NONNEGATIVE and number < 0.0:
        return "value is negative"
    if domain == DOMAIN_PERCENT_0_100 and not 0.0 <= number <= 100.0:
        return "percentage is outside [0, 100]"
    return None


def validate_native_predictions(predictions: Any) -> list[str]:
    """Validate the exact 41-endpoint scorer payload, returning all issues."""
    if not isinstance(predictions, dict):
        return ["predictions is not an object"]

    expected = set(ADMET_NATIVE_ENDPOINT_IDS)
    actual = set(predictions)
    issues: list[str] = []
    missing = sorted(expected - actual)
    unexpected = sorted(actual - expected)
    if missing:
        issues.append(f"missing endpoints: {missing}")
    if unexpected:
        issues.append(f"unexpected endpoints: {unexpected}")

    for endpoint_id in ADMET_NATIVE_ENDPOINT_IDS:
        if endpoint_id not in predictions:
            continue
        record = predictions[endpoint_id]
        if not isinstance(record, dict):
            issues.append(f"endpoint {endpoint_id!r} is not an object")
            continue
        error = endpoint_value_error(endpoint_id, record.get("value"))
        if error:
            issues.append(f"endpoint {endpoint_id!r} {error}")
    return issues


assert len(ADMET_NATIVE_ENDPOINT_IDS) == 41
assert len(set(ADMET_NATIVE_ENDPOINT_IDS)) == len(ADMET_NATIVE_ENDPOINT_IDS)
assert set(ADMET_ENDPOINT_DOMAINS) == set(ADMET_NATIVE_ENDPOINT_IDS)
assert all(
    set(endpoint_ids) <= ADMET_PROBABILITY_ENDPOINT_IDS
    for endpoint_ids in ADMET_DERIVED_COUNT_INPUTS.values()
)

"""Predefined Reaction SMARTS rules for molecular modification.

Each rule maps a keyword from fix_instruction to a list of (SMARTS_pattern, SMARTS_product) pairs.
The engine tries each in order until one produces a valid product.
"""
from __future__ import annotations

REACTION_RULES: dict[str, list[tuple[str, str]]] = {
    "H-bond donor": [
        # -CH3 → -CH2OH
        ("[CH3:1]", "[CH2:1]O"),
        # aromatic C-H → C-OH
        ("[cH1:1]", "[c:1]O"),
        # -CH3 → -CH2NH2
        ("[CH3:1]", "[CH2:1]N"),
    ],
    "H-bond acceptor": [
        # -CH3 → -CH2OCH3
        ("[CH3:1]", "[CH2:1]OC"),
        # aromatic C-H → C-OMe
        ("[cH1:1]", "[c:1]OC"),
    ],
    "hydrophobic": [
        # aromatic C-H → C-Me
        ("[cH1:1]", "[c:1]C"),
        # -OH → -OC(C)C
        ("[OH1:1]", "[O:1]C(C)C"),
    ],
    "aromatic": [
        # -CH3 → -CH2-phenyl
        ("[CH3:1]", "[CH2:1]c1ccccc1"),
        # aromatic C-H → C-phenyl
        ("[cH1:1]", "[c:1]c1ccccc1"),
    ],
    "charged": [
        # -CH3 → -CH2NH3+ (protonated amine)
        ("[CH3:1]", "[CH2:1][NH3+]"),
        # aromatic C-H → C-COOH
        ("[cH1:1]", "[c:1]C(=O)O"),
    ],
    "reduce strain": [
        # remove one methyl from a gem-dimethyl carbon (demethylation).
        # [CX4] pins the central atom to sp3 carbon so ring/aromatic carbons and
        # heteroatoms are not matched; the deleted methyl is the unmapped [CH3].
        ("[CX4:1]([CH3:2])[CH3]", "[C:1][CH3:2]"),
    ],
    "trim clash": [
        # trim a TERMINAL ethyl to methyl: the [CH3] must be a true chain end,
        # and [CH2] must connect to exactly one heavy atom besides that methyl.
        # Degree/H specifiers ([CH2X4], [CH3X4]) keep this from matching internal
        # -CH2-CH3-like fragments and cutting a molecule in half.
        ("[CH2X4:1][CH3X4]", "[CH3:1]"),
        # remove a terminal methoxy -OCH3 → -OH ([OX2] = two-connection ether O).
        ("[OX2:1][CH3X4]", "[OH:1]"),
    ],
}


# Rules that only make sense as explicit, keyword-matched requests. When a
# fix_instruction matches nothing, we fall back to a mild H-bond-donor edit but
# the caller is told it was an unrecognized instruction (see match_keywords).
DEFAULT_FALLBACK_RULE = "H-bond donor"


def match_keywords(instruction: str) -> list[str]:
    """Extract matching rule keys from a fix_instruction string.

    Returns the matched rule keys. When nothing matches, returns the fallback
    rule; use `matched_any()` (or check against DEFAULT_FALLBACK_RULE) to detect
    that the instruction was unrecognized rather than silently assuming a donor
    edit was intended.
    """
    return _match(instruction) or [DEFAULT_FALLBACK_RULE]


def matched_any(instruction: str) -> bool:
    """Return True if the instruction matched at least one known rule keyword."""
    return bool(_match(instruction))


def _match(instruction: str) -> list[str]:
    instruction_lower = (instruction or "").lower()
    matched: list[str] = []
    keyword_map = {
        "h-bond donor": "H-bond donor",
        "hbond donor": "H-bond donor",
        "donor": "H-bond donor",
        "h-bond acceptor": "H-bond acceptor",
        "hbond acceptor": "H-bond acceptor",
        "acceptor": "H-bond acceptor",
        "hydrophobic": "hydrophobic",
        "lipophilic": "hydrophobic",
        "aromatic": "aromatic",
        "charged": "charged",
        "salt": "charged",
        "strain": "reduce strain",
        "simplify": "reduce strain",
        "clash": "trim clash",
        "trim": "trim clash",
    }
    for trigger, rule_key in keyword_map.items():
        if trigger in instruction_lower and rule_key not in matched:
            matched.append(rule_key)
    return matched

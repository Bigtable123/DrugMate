"""Shared residue utilities for SBDD skills.

``normalize_residue`` and ``extract_resnum`` are legacy helpers used by several
skills and intentionally keep their historical, chain-agnostic behaviour.  New
pocket code should use the record helpers below instead: a residue is only
unambiguous when chain and insertion code are part of its identity.
"""
from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from typing import Any, Optional


ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def normalize_residue(res: str) -> str:
    """Normalize residue ID to canonical form: 'LEU83.A' / 'Leu83' / 'L83' → 'LEU83'.

    Handles:
    - Chain suffix: "LEU83.A" → "LEU83"
    - Underscore prefix: "A_LEU83" → "LEU83"
    - Single-letter codes: "L83" → "LEU83"
    - Mixed case: "leu83" → "LEU83"
    """
    cleaned = res.split(".")[0].strip()
    if "_" in cleaned:
        cleaned = cleaned.split("_")[-1]
    match = re.match(r"([A-Za-z]+)(\d+)", cleaned)
    if not match:
        return cleaned.upper()
    name = match.group(1).upper()
    num = match.group(2)
    # Convert single-letter to three-letter
    if len(name) == 1:
        name = ONE_TO_THREE.get(name, name)
    else:
        name = name[:3]
    return f"{name}{num}"


def extract_resnum(residue_str: str) -> Optional[int]:
    """Extract residue number from a residue string."""
    match = re.search(r"(\d+)", residue_str)
    return int(match.group(1)) if match else None


# ---------------------------------------------------------------------------
# Chain/insertion-code-safe residue records

_NAMED_RESIDUE_RE = re.compile(
    # PDB component names are one to three alphanumeric characters, so a digit
    # can appear anywhere in them (``1PE``, ``ZN``).  ``resname`` is therefore
    # matched lazily: a greedy class swallowed the head of the residue number
    # instead, turning ``ZN301`` into resname ``ZN3`` at residue 1.  Lazy
    # matching hands the whole trailing digit run to ``resnum``, which is the
    # only reading in which the number survives.
    r"^(?P<resname>[A-Za-z0-9]{1,3}?)(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)$"
)
_NUMBERED_RESIDUE_RE = re.compile(r"^(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)$")


def make_residue_record(
    resname: str,
    resnum: int | str,
    chain: str = "",
    icode: str = "",
    **metadata: Any,
) -> dict[str, Any]:
    """Create a JSON-serializable, chain-safe residue record.

    The stable identity fields are ``chain``, ``resnum``, ``icode`` and
    ``resname``.  ``label`` is a human-readable canonical label while
    ``legacy_label`` preserves the old ``LEU83``-style representation.
    Additional metadata (for example ``identity_ambiguous``) is retained.
    """
    name = str(resname or "").strip().upper()
    if len(name) == 1:
        name = ONE_TO_THREE.get(name, name)
    elif name:
        name = name[:3]

    try:
        number: int | str = int(str(resnum).strip())
    except (TypeError, ValueError):
        number = str(resnum or "").strip()

    chain_id = str(chain or "").strip()
    insertion = str(icode or "").strip()
    legacy = f"{name}{number}{insertion}"
    label = f"{chain_id}:{legacy}" if chain_id else legacy
    record: dict[str, Any] = {
        "chain": chain_id,
        "resname": name,
        "resnum": number,
        "icode": insertion,
        "label": label,
        "legacy_label": legacy,
    }
    record.update(metadata)
    return record


def parse_residue_record(
    value: Any,
    default_chain: str = "",
    default_resname: str = "",
) -> dict[str, Any] | None:
    """Coerce common residue representations to a residue record.

    Accepted strings include ``LEU83``, ``LEU83A``, ``A:LEU83A``,
    ``A_LEU83``, ``LEU83.A`` (ProLIF), and numeric P2Rank labels such as
    ``A_83`` when ``default_resname`` is supplied.  Mapping inputs may use the
    canonical field names or the common PDB aliases ``chain_id``, ``res_seq``
    and ``i_code``.
    """
    if value is None:
        return None

    metadata: dict[str, Any] = {}
    if isinstance(value, Mapping):
        mapping = dict(value)
        chain = mapping.get("chain", mapping.get("chain_id", default_chain))
        resname = mapping.get(
            "resname", mapping.get("res_name", mapping.get("name", default_resname))
        )
        resnum = mapping.get(
            "resnum", mapping.get("res_seq", mapping.get("seq", mapping.get("number")))
        )
        icode = mapping.get(
            "icode", mapping.get("i_code", mapping.get("insertion_code", ""))
        )
        if resnum is None and mapping.get("label"):
            parsed = parse_residue_record(
                mapping["label"], default_chain=str(chain or default_chain),
                default_resname=str(resname or default_resname),
            )
            if parsed is None:
                return None
            for key, item in mapping.items():
                if key not in parsed:
                    parsed[key] = item
            return parsed
        if resnum is None:
            return None
        canonical_fields = {
            "chain", "chain_id", "resname", "res_name", "name", "resnum",
            "res_seq", "seq", "number", "icode", "i_code", "insertion_code",
            "label", "legacy_label",
        }
        metadata = {k: v for k, v in mapping.items() if k not in canonical_fields}
        return make_residue_record(resname, resnum, chain, icode, **metadata)

    raw = str(value).strip()
    if not raw:
        return None

    chain = str(default_chain or "").strip()
    token = raw

    # Canonical output form: A:LEU83A
    if ":" in token:
        prefix, remainder = token.split(":", 1)
        if prefix and remainder:
            chain, token = prefix.strip(), remainder.strip()
    # ProLIF form: LEU83.A
    elif "." in token:
        residue_part, suffix = token.rsplit(".", 1)
        if suffix:
            token, chain = residue_part.strip(), suffix.strip()
    # P2Rank forms: A_LEU83 / A_83
    elif "_" in token:
        prefix, remainder = token.split("_", 1)
        if prefix and remainder:
            chain, token = prefix.strip(), remainder.strip()

    token = token.replace(" ", "")
    # The bare-number form is tried first because it is the only unambiguous one.
    # A named pattern that also accepts digits matches ``45`` as component ``4``
    # at residue 5 and ``301A`` as component ``30`` at residue 1 -- a P2Rank or
    # PDB numeric label silently became a different residue.
    match = _NUMBERED_RESIDUE_RE.fullmatch(token)
    if match:
        return make_residue_record(
            default_resname, match.group("resnum"), chain, match.group("icode")
        )

    match = _NAMED_RESIDUE_RE.fullmatch(token)
    if match:
        return make_residue_record(
            match.group("resname"), match.group("resnum"), chain, match.group("icode")
        )
    return None


def residue_identity(
    value: Any,
    default_chain: str = "",
    *,
    include_resname: bool = True,
    allow_ambiguous: bool = False,
) -> tuple[Any, ...] | None:
    """Return a hashable residue identity, or ``None`` when it is ambiguous."""
    record = parse_residue_record(value, default_chain=default_chain)
    if record is None:
        return None
    if record.get("identity_ambiguous") and not allow_ambiguous:
        return None
    number = record.get("resnum")
    if number in (None, ""):
        return None
    base: tuple[Any, ...] = (
        record.get("chain", ""), number, record.get("icode", ""),
    )
    if include_resname:
        return base + (record.get("resname", ""),)
    return base


def legacy_residue_label(value: Any, default_chain: str = "") -> str:
    """Return the historical ``LEU83``-style label for a residue value."""
    record = parse_residue_record(value, default_chain=default_chain)
    if record is None:
        return str(value)
    return str(record.get("legacy_label") or record.get("label") or value)


def canonical_residue_label(value: Any, default_chain: str = "") -> str:
    """Return a chain-aware human-readable label such as ``A:LEU83A``."""
    record = parse_residue_record(value, default_chain=default_chain)
    if record is None:
        return str(value)
    return str(record.get("label") or value)


def deduplicate_residue_records(
    values: Iterable[Any], default_chain: str = "",
) -> list[dict[str, Any]]:
    """Normalize and deterministically de-duplicate residue records.

    Ambiguous records use their canonical label plus ambiguity metadata as a
    conservative key; they are never silently merged with an exact identity.
    """
    records: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    for value in values or []:
        record = parse_residue_record(value, default_chain=default_chain)
        if record is None:
            continue
        identity = residue_identity(record)
        if identity is None:
            key = (
                "ambiguous", record.get("label", ""),
                tuple(record.get("candidate_icodes") or []),
            )
        else:
            key = ("exact",) + identity
        if key in seen:
            continue
        seen.add(key)
        records.append(record)
    return sorted(
        records,
        key=lambda item: (
            str(item.get("chain", "")),
            _sortable_resnum(item.get("resnum")),
            str(item.get("icode", "")),
            str(item.get("resname", "")),
        ),
    )


def residue_records_from_values(
    values: Iterable[Any], default_chain: str = "",
) -> list[dict[str, Any]]:
    """Public alias for normalizing a list of legacy labels or records."""
    return deduplicate_residue_records(values, default_chain=default_chain)


def _sortable_resnum(value: Any) -> tuple[int, Any]:
    try:
        return 0, int(value)
    except (TypeError, ValueError):
        return 1, str(value)

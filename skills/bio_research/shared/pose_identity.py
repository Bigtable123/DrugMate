"""Deterministic chemical, pose-lineage, and batch-cache identities.

The chemical and pose hashing helpers are deliberately strict. A missing source
identity, an unparsable molecule, non-finite coordinates, or unsupported
metadata raises a typed error instead of producing a fallback identifier that
could collide with a valid record.

The JSON-shaped return values are suitable for embedding directly in SBDD
artifacts.  Hashes include an explicit schema namespace so future contract
changes invalidate old cache entries safely.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import hmac
import json
import math
import os
from pathlib import Path
import re
from typing import Any, Mapping, Sequence


CHEMICAL_IDENTITY_SCHEMA = "sbdd.chemical-identity.v1"
POSE_IDENTITY_SCHEMA = "sbdd.pose-identity.v1"
BATCH_FINGERPRINT_SCHEMA = "sbdd.structure-evaluator.batch-fingerprint.v1"
BATCH_OUTPUT_DIGEST_SCHEMA = "sbdd.structure-evaluator.batch-output.v1"
SDF_MOLBLOCK_DIGEST_SCHEMA = "sbdd.sdf-molblock.v1"
_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

# Ordinary prepared libraries carry these fields when available. They do not
# assert a pH population model or require an enumerated microstate ensemble.
LIGAND_LINEAGE_PROPERTIES = (
    "source_record_id", "source_tool", "input_compound_id", "input_mol_name", "compound_id",
    "microstate_id", "conformer_id", "conformer_coordinate_identity",
    "coordinate_smiles", "smiles_standardized", "candidate_role", "source",
    "source_tools", "also_generated_by", "cohort_id", "source_cohort_id", "campaign_mode",
    "benchmark_dataset", "benchmark_split", "coordinate_origin",
    "raw_pose_preserved", "stereochemistry_status", "preparation_warnings_json",
    "molecule_failure_report", "source_cohort_report", "parent_input_records_json",
    "parent_id", "parent_smiles", "campaign_candidate_id", "campaign_parent_ids_json",
)


def complete_ligand_ids(molecules: Sequence[Any], *, source_token: str) -> None:
    """Fill absent record IDs before selection, retaining declared identities.

    Titles are display names and need not be unique. Only an unambiguous title
    can serve as an ID; otherwise the source record supplies a stable fallback.
    Duplicate explicit IDs remain visible for the docking boundary to reject.
    """
    def prop(mol: Any, *keys: str) -> str:
        return next((mol.GetProp(key).strip() for key in keys
                     if mol is not None and mol.HasProp(key)
                     and mol.GetProp(key).strip()), "")

    declared = [prop(mol, "compound_id", "CompoundID", "compound_name")
                for mol in molecules]
    titles = [prop(mol, "_Name") for mol in molecules]
    counts = Counter(titles)
    declared_ids = set(declared)
    reserved = declared_ids | {name for name in titles if counts[name] == 1}
    for index, (mol, compound_id, title) in enumerate(zip(molecules, declared, titles)):
        if mol is None:
            continue
        if not mol.HasProp("input_mol_name"):
            mol.SetProp("input_mol_name", mol.GetProp("_Name") if mol.HasProp("_Name") else "")
        if not compound_id:
            if title and counts[title] == 1 and title not in declared_ids:
                compound_id = title
            else:
                source_record = prop(mol, "source_record_id") or f"{source_token}:record:{index}"
                digest = hashlib.sha256(source_record.encode("utf-8")).hexdigest()
                compound_id = f"ligand_{digest}"
                # Repeated source_record_id values or user IDs must never cause
                # two anonymous physical records to acquire the same ID.
                while compound_id in reserved:
                    compound_id += f"_{index}"
                mol.SetProp("source_record_id", source_record)
                if title and not prop(mol, "input_compound_id"):
                    mol.SetProp("input_compound_id", title)
            reserved.add(compound_id)
        mol.SetProp("compound_id", compound_id)


def benchmark_label(record: Mapping[str, Any]) -> str | None:
    """Read original ACT_/DEC_ labels without changing the engine join key."""
    for key in ("input_mol_name", "mol_name", "compound_id"):
        name = str(record.get(key) or "").strip()
        for label in ("ACT", "DEC"):
            if name.startswith(label + "_"):
                return label
    return None


def sdf_molblock_sha256(record: str) -> str:
    """Hash the coordinate-bearing mol block of one physical SDF record.

    SD properties are deliberately excluded: cohort annotations may be added
    after an engine emits its result, whereas the title, atom/bond tables, and
    coordinates must remain the same.  Newlines and trailing blank separator
    lines are normalized so a text-only property rewrite cannot invalidate the
    witness.
    """
    if not isinstance(record, str):
        raise TypeError("SDF record must be text")
    normalized = record.replace("\r\n", "\n").replace("\r", "\n")
    molblock_lines: list[str] = []
    found_terminator = False
    for line_index, line in enumerate(normalized.split("\n")):
        molblock_lines.append(line)
        # Both V2000 and V3000 mol blocks terminate with this CTAB line.  Do
        # not detect the property section from a line beginning with ``>``:
        # that text is legal in an SDF title and would truncate a valid block.
        # A zero-atom mol block still has the three header lines plus a counts
        # line before its terminator.  Requiring that prefix prevents a legal
        # title/comment literally equal to ``M  END`` from truncating the
        # witness into a one-line hash shared by unrelated poses.
        if line_index >= 4 and line.strip() == "M  END":
            found_terminator = True
            break
    if not found_terminator:
        raise PoseIdentityError(
            "missing_sdf_molblock_terminator",
            "cannot build a pose-record witness without an M  END terminator",
        )
    payload = (
        SDF_MOLBLOCK_DIGEST_SCHEMA
        + "\n"
        + "\n".join(molblock_lines)
        + "\n"
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


class PoseIdentityError(ValueError):
    """A molecule or coordinate set cannot be assigned a safe identity."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.details = dict(details or {})

    def as_result(self) -> dict[str, Any]:
        """Return a JSON-safe error result without any identity fields."""
        result: dict[str, Any] = {
            "status": "error",
            "error": {
                "code": self.code,
                "message": str(self),
            },
        }
        if self.details:
            result["error"]["details"] = _normalize_json(
                self.details, path="$.error.details"
            )
        return result


class BatchFingerprintError(ValueError):
    """A batch execution contract cannot be fingerprinted unambiguously."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.details = dict(details or {})


def _normalize_json(value: Any, *, path: str = "$") -> Any:
    """Convert supported values to a deterministic, strict JSON shape.

    Unknown Python objects and non-finite floats are rejected.  In particular,
    this function never falls back to ``str(object)`` because unrelated
    objects can share the same display string and silently collide.
    """
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise BatchFingerprintError(
                "non_finite_number",
                f"{path} contains a non-finite floating-point value",
            )
        return 0.0 if value == 0.0 else value
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if isinstance(value, (bytes, bytearray, memoryview)):
        raw = bytes(value)
        return {
            "$bytes_sha256": hashlib.sha256(raw).hexdigest(),
            "size_bytes": len(raw),
        }
    if isinstance(value, Mapping):
        normalized: dict[str, Any] = {}
        for key, child in value.items():
            if not isinstance(key, str):
                raise BatchFingerprintError(
                    "non_string_mapping_key",
                    f"{path} contains a non-string mapping key: {key!r}",
                )
            normalized[key] = _normalize_json(child, path=f"{path}.{key}")
        return {key: normalized[key] for key in sorted(normalized)}
    if isinstance(value, (list, tuple)):
        return [
            _normalize_json(child, path=f"{path}[{index}]")
            for index, child in enumerate(value)
        ]
    if isinstance(value, (set, frozenset)):
        raise BatchFingerprintError(
            "unordered_collection",
            f"{path} contains an unordered collection; pass a sorted list",
        )
    raise BatchFingerprintError(
        "unsupported_metadata_type",
        f"{path} contains unsupported type {type(value).__name__}",
    )


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize a supported value to stable UTF-8 JSON bytes."""
    normalized = _normalize_json(value)
    return json.dumps(
        normalized,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def stable_json_hash(value: Any, *, namespace: str) -> str:
    """Return a namespaced SHA-256 identifier for a JSON-shaped value."""
    if not isinstance(namespace, str) or not namespace.strip():
        raise BatchFingerprintError(
            "invalid_namespace", "hash namespace must be a non-empty string"
        )
    digest = hashlib.sha256(
        canonical_json_bytes({"namespace": namespace, "value": value})
    ).hexdigest()
    return f"sha256:{digest}"


def _rdkit_modules():
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import rdMolDescriptors
    except ImportError as exc:
        raise PoseIdentityError(
            "rdkit_unavailable",
            "RDKit is required to build chemical and pose identities",
        ) from exc
    return Chem, rdBase, rdMolDescriptors


def _prepare_molecule(mol: Any, *, remove_hydrogens: bool) -> Any:
    Chem, _, _ = _rdkit_modules()
    if mol is None:
        raise PoseIdentityError("missing_molecule", "molecule is None")
    if not hasattr(mol, "GetAtoms") or not hasattr(mol, "GetNumAtoms"):
        raise PoseIdentityError(
            "invalid_molecule_type",
            f"expected an RDKit molecule, got {type(mol).__name__}",
        )
    try:
        prepared = Chem.Mol(mol)
    except Exception as exc:
        raise PoseIdentityError(
            "molecule_copy_failed", f"could not copy RDKit molecule: {exc}"
        ) from exc
    if prepared.GetNumAtoms() == 0:
        raise PoseIdentityError("empty_molecule", "molecule has no atoms")

    for atom in prepared.GetAtoms():
        if atom.GetAtomicNum() <= 0 or atom.HasQuery():
            raise PoseIdentityError(
                "unsupported_query_atom",
                "query, wildcard, and zero-atomic-number atoms are not "
                "supported for stable chemical identity",
                details={"atom_index": atom.GetIdx()},
            )
        # Atom-map numbers describe a representation/provenance convention,
        # not chemical identity.  Clear them before canonicalization.
        atom.SetAtomMapNum(0)

    try:
        Chem.SanitizeMol(prepared)
        if remove_hydrogens:
            # Normal explicit hydrogens are representation details.  RDKit's
            # default RemoveHs parameters preserve isotopic hydrogens.
            prepared = Chem.RemoveHs(prepared, sanitize=True)
            Chem.SanitizeMol(prepared)
        Chem.AssignStereochemistry(prepared, cleanIt=True, force=True)
    except Exception as exc:
        raise PoseIdentityError(
            "molecule_sanitization_failed",
            f"RDKit could not sanitize/canonicalize molecule: {exc}",
        ) from exc

    if not any(atom.GetAtomicNum() > 1 for atom in prepared.GetAtoms()):
        raise PoseIdentityError(
            "no_heavy_atoms", "molecule contains no heavy atoms"
        )
    return prepared


def matches_declared_molecule(actual: Any, declared: Any) -> bool:
    """Match the full graph and only the stereochemistry declared by the input.

    Output coordinates may resolve unspecified input stereo. They must not
    invert or erase specified stereo, or change bonds, charge, or isotopes.
    This directional compatibility check does not replace exact pose hashes.
    """
    if actual is None or declared is None:
        return False
    Chem, _, _ = _rdkit_modules()
    actual = _prepare_molecule(actual, remove_hydrogens=True)
    declared = _prepare_molecule(declared, remove_hydrogens=True)
    if Chem.MolToSmiles(actual) == Chem.MolToSmiles(declared):
        return True
    graphs = []
    for molecule in (actual, declared):
        graph = Chem.Mol(molecule)
        Chem.RemoveStereochemistry(graph)
        # Keep isomericSmiles=True so isotope changes remain detectable.
        graphs.append(Chem.MolToSmiles(graph, isomericSmiles=True))
    return graphs[0] == graphs[1] and actual.HasSubstructMatch(
        declared, useChirality=True,
    )


def _require_source_identity(source_identity: Any) -> Any:
    if source_identity is None:
        raise PoseIdentityError(
            "missing_source_identity",
            "source_identity is required to prevent cross-source collisions",
        )
    if isinstance(source_identity, str) and not source_identity.strip():
        raise PoseIdentityError(
            "empty_source_identity", "source_identity cannot be an empty string"
        )
    if isinstance(source_identity, (Mapping, list, tuple)) and not source_identity:
        raise PoseIdentityError(
            "empty_source_identity",
            "source_identity cannot be an empty mapping or sequence",
        )
    if isinstance(source_identity, bool):
        raise PoseIdentityError(
            "invalid_source_identity",
            "source_identity cannot be a boolean",
        )
    try:
        return _normalize_json(source_identity, path="$.source_identity")
    except BatchFingerprintError as exc:
        raise PoseIdentityError(
            "invalid_source_identity",
            f"source_identity is not deterministic JSON metadata: {exc}",
        ) from exc


def _normalize_optional_mapping(
    value: Mapping[str, Any] | None,
    *,
    path: str,
) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise BatchFingerprintError(
            "invalid_metadata_mapping",
            f"{path} must be a mapping or null",
        )
    return _normalize_json(value, path=path)


def _stereochemistry_summary(Chem: Any, mol: Any) -> dict[str, Any]:
    try:
        centers = Chem.FindMolChiralCenters(
            mol,
            includeUnassigned=True,
            includeCIP=True,
            useLegacyImplementation=False,
        )
    except TypeError:
        centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, includeCIP=True
        )
    tetrahedral = Counter(str(label) for _, label in centers)

    double_bonds: Counter[str] = Counter()
    for bond in mol.GetBonds():
        stereo = str(bond.GetStereo())
        if stereo not in {"STEREONONE", "STEREOANY"}:
            double_bonds[stereo] += 1
        elif stereo == "STEREOANY":
            double_bonds["unassigned"] += 1

    potential: Counter[str] = Counter()
    try:
        for info in Chem.FindPotentialStereo(mol):
            key = ":".join(
                (
                    str(getattr(info, "type", "unknown")),
                    str(getattr(info, "specified", "unknown")),
                    str(getattr(info, "descriptor", "unknown")),
                )
            )
            potential[key] += 1
    except (AttributeError, RuntimeError):
        # Canonical isomeric SMILES remains the authoritative stereo identity.
        # This summary is diagnostic and older supported RDKit builds may not
        # expose FindPotentialStereo.
        pass

    return {
        "tetrahedral_center_count": len(centers),
        "tetrahedral_assignments": dict(sorted(tetrahedral.items())),
        "double_bond_assignments": dict(sorted(double_bonds.items())),
        "potential_stereo": dict(sorted(potential.items())),
        "has_unassigned": (
            tetrahedral.get("?", 0) > 0
            or double_bonds.get("unassigned", 0) > 0
            or any("Unspecified" in key for key in potential)
        ),
    }


def build_chemical_identity(
    mol: Any,
    *,
    source_identity: Any,
) -> dict[str, Any]:
    """Build a deterministic chemical identity from an RDKit molecule.

    ``chemical_id`` identifies the normalized molecular graph and remains the
    same across source records. ``source_identity_id`` identifies the caller's
    compound/microstate/variant provenance. ``identity_id`` binds both.
    Coordinates and conformer identifiers are intentionally excluded.
    """
    Chem, rdBase, rdMolDescriptors = _rdkit_modules()
    normalized_source = _require_source_identity(source_identity)
    prepared = _prepare_molecule(mol, remove_hydrogens=True)

    try:
        canonical_smiles = Chem.MolToSmiles(
            prepared,
            canonical=True,
            isomericSmiles=True,
            allHsExplicit=False,
            allBondsExplicit=False,
        )
        formula = rdMolDescriptors.CalcMolFormula(
            prepared, separateIsotopes=True, abbreviateHIsotopes=False
        )
    except Exception as exc:
        raise PoseIdentityError(
            "canonicalization_failed",
            f"RDKit could not create canonical chemical descriptors: {exc}",
        ) from exc
    if not canonical_smiles:
        raise PoseIdentityError(
            "canonicalization_failed", "canonical isomeric SMILES is empty"
        )

    isotopes = Counter(
        (atom.GetAtomicNum(), atom.GetIsotope())
        for atom in prepared.GetAtoms()
        if atom.GetIsotope()
    )
    isotope_summary = [
        {
            "atomic_number": atomic_number,
            "mass_number": mass_number,
            "count": count,
        }
        for (atomic_number, mass_number), count in sorted(isotopes.items())
    ]
    structural_payload = {
        "schema_version": CHEMICAL_IDENTITY_SCHEMA,
        "canonical_isomeric_smiles": canonical_smiles,
        "formal_charge": sum(
            atom.GetFormalCharge() for atom in prepared.GetAtoms()
        ),
        "heavy_atom_count": sum(
            atom.GetAtomicNum() > 1 for atom in prepared.GetAtoms()
        ),
        "fragment_count": len(Chem.GetMolFrags(prepared)),
        "isotopes": isotope_summary,
    }
    chemical_id = stable_json_hash(
        structural_payload, namespace=CHEMICAL_IDENTITY_SCHEMA
    ).replace("sha256:", "chem:", 1)
    source_identity_id = stable_json_hash(
        normalized_source, namespace=f"{CHEMICAL_IDENTITY_SCHEMA}.source"
    ).replace("sha256:", "source:", 1)
    identity_id = stable_json_hash(
        {
            "chemical_id": chemical_id,
            "source_identity_id": source_identity_id,
        },
        namespace=f"{CHEMICAL_IDENTITY_SCHEMA}.bound",
    ).replace("sha256:", "chemical-source:", 1)

    return {
        "status": "ok",
        "schema_version": CHEMICAL_IDENTITY_SCHEMA,
        "chemical_id": chemical_id,
        "source_identity_id": source_identity_id,
        "identity_id": identity_id,
        **structural_payload,
        "molecular_formula": formula,
        "stereochemistry": _stereochemistry_summary(Chem, prepared),
        "source_identity": normalized_source,
        "producer": {
            "library": "rdkit",
            "version": rdBase.rdkitVersion,
        },
    }


def chemical_identity_result(
    mol: Any,
    *,
    source_identity: Any,
) -> dict[str, Any]:
    """Return a chemical identity or a structured error with no fallback ID."""
    try:
        return build_chemical_identity(mol, source_identity=source_identity)
    except PoseIdentityError as exc:
        result = exc.as_result()
        result["schema_version"] = CHEMICAL_IDENTITY_SCHEMA
        return result


def _coordinate_signature(
    mol: Any,
    *,
    conformer_id: int,
    coordinate_precision: int,
) -> tuple[dict[str, Any], int]:
    Chem, _, _ = _rdkit_modules()
    prepared = _prepare_molecule(mol, remove_hydrogens=False)
    if isinstance(conformer_id, bool) or not isinstance(conformer_id, int):
        raise PoseIdentityError(
            "invalid_conformer_id", "conformer_id must be an integer"
        )
    if (
        isinstance(coordinate_precision, bool)
        or not isinstance(coordinate_precision, int)
        or not 0 <= coordinate_precision <= 12
    ):
        raise PoseIdentityError(
            "invalid_coordinate_precision",
            "coordinate_precision must be an integer between 0 and 12",
        )
    try:
        conformer = prepared.GetConformer(conformer_id)
    except Exception as exc:
        raise PoseIdentityError(
            "missing_conformer",
            f"molecule has no conformer with id {conformer_id}",
        ) from exc

    try:
        symmetry_classes = list(
            Chem.CanonicalRankAtoms(
                prepared,
                breakTies=False,
                includeChirality=True,
                includeIsotopes=True,
            )
        )
    except Exception as exc:
        raise PoseIdentityError(
            "coordinate_canonicalization_failed",
            f"could not compute atom symmetry classes: {exc}",
        ) from exc

    atom_records: list[list[Any]] = []
    for atom in prepared.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        coordinates = [float(position.x), float(position.y), float(position.z)]
        if not all(math.isfinite(value) for value in coordinates):
            raise PoseIdentityError(
                "non_finite_coordinate",
                "pose contains a non-finite coordinate",
                details={"atom_index": atom.GetIdx()},
            )
        rounded = [
            0.0 if round(value, coordinate_precision) == 0.0
            else round(value, coordinate_precision)
            for value in coordinates
        ]
        atom_records.append(
            [
                symmetry_classes[atom.GetIdx()],
                atom.GetAtomicNum(),
                atom.GetIsotope(),
                atom.GetFormalCharge(),
                str(atom.GetChiralTag()),
                *rounded,
            ]
        )
    atom_records.sort()
    coordinate_payload = {
        "atom_count": len(atom_records),
        "is_3d": bool(conformer.Is3D()),
        "precision_decimals": coordinate_precision,
        "atom_records": atom_records,
    }
    signature = stable_json_hash(
        coordinate_payload, namespace=f"{POSE_IDENTITY_SCHEMA}.coordinates"
    ).replace("sha256:", "coordinates:", 1)
    return (
        {
            "signature": signature,
            "atom_count": len(atom_records),
            "is_3d": bool(conformer.Is3D()),
            "precision_decimals": coordinate_precision,
        },
        conformer.GetId(),
    )


def build_pose_identity(
    mol: Any,
    *,
    source_identity: Any,
    parent_pose_id: str | None = None,
    coordinate_revision: int = 0,
    conformer_id: int = -1,
    coordinate_precision: int = 6,
    lineage_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a stable pose ID and explicit coordinate lineage.

    Root poses use ``coordinate_revision=0`` and no parent.  A
    coordinate-changing operation must provide its parent and a positive
    revision.  Conformer IDs are reported but excluded from the hash because
    RDKit may renumber them during serialization.
    """
    if (
        isinstance(coordinate_revision, bool)
        or not isinstance(coordinate_revision, int)
        or coordinate_revision < 0
    ):
        raise PoseIdentityError(
            "invalid_coordinate_revision",
            "coordinate_revision must be a non-negative integer",
        )
    if parent_pose_id is None and coordinate_revision != 0:
        raise PoseIdentityError(
            "missing_parent_pose_id",
            "a non-zero coordinate_revision requires parent_pose_id",
        )
    if parent_pose_id is not None:
        if not isinstance(parent_pose_id, str) or not parent_pose_id.strip():
            raise PoseIdentityError(
                "invalid_parent_pose_id",
                "parent_pose_id must be a non-empty string",
            )
        if coordinate_revision == 0:
            raise PoseIdentityError(
                "invalid_coordinate_revision",
                "a child pose must use coordinate_revision >= 1",
            )

    try:
        normalized_lineage = _normalize_optional_mapping(
            lineage_metadata, path="$.lineage_metadata"
        )
    except BatchFingerprintError as exc:
        raise PoseIdentityError(
            "invalid_lineage_metadata",
            f"lineage_metadata is not deterministic JSON metadata: {exc}",
        ) from exc
    chemical_identity = build_chemical_identity(
        mol, source_identity=source_identity
    )
    coordinates, resolved_conformer_id = _coordinate_signature(
        mol,
        conformer_id=conformer_id,
        coordinate_precision=coordinate_precision,
    )
    lineage = {
        "parent_pose_id": parent_pose_id,
        "coordinate_revision": coordinate_revision,
    }
    pose_payload = {
        "schema_version": POSE_IDENTITY_SCHEMA,
        "chemical_identity_id": chemical_identity["identity_id"],
        "coordinate_signature": coordinates["signature"],
        "lineage": lineage,
    }
    pose_id = stable_json_hash(
        pose_payload, namespace=POSE_IDENTITY_SCHEMA
    ).replace("sha256:", "pose:", 1)
    return {
        "status": "ok",
        "schema_version": POSE_IDENTITY_SCHEMA,
        "pose_id": pose_id,
        "parent_pose_id": parent_pose_id,
        "coordinate_revision": coordinate_revision,
        "conformer_id": resolved_conformer_id,
        "coordinates": coordinates,
        "lineage_metadata": normalized_lineage,
        "chemical_identity": chemical_identity,
    }


def _validate_sha256(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise BatchFingerprintError(
            "invalid_content_hash", f"{label} SHA-256 must be a string"
        )
    digest = value.removeprefix("sha256:")
    if not _SHA256_RE.fullmatch(digest):
        raise BatchFingerprintError(
            "invalid_content_hash",
            f"{label} is not a 64-character SHA-256 digest",
        )
    return digest.lower()


def artifact_signature(source: Any, *, label: str) -> dict[str, str]:
    """Return a SHA-256 content signature from bytes, a path, or a digest.

    Digest inputs must be either ``"sha256:<hex>"``/bare hex or a mapping with
    ``sha256`` (or ``algorithm="sha256", digest=...``).  Other strings are
    treated as file paths and must exist.
    """
    digest: str
    if isinstance(source, (bytes, bytearray, memoryview)):
        digest = hashlib.sha256(bytes(source)).hexdigest()
    elif isinstance(source, Mapping):
        allowed_keys = {"algorithm", "digest", "sha256", "size_bytes"}
        unknown_keys = set(source) - allowed_keys
        if unknown_keys:
            raise BatchFingerprintError(
                "invalid_content_signature",
                f"{label} signature has unsupported keys: "
                f"{sorted(map(repr, unknown_keys))!r}",
            )
        if "algorithm" in source and source["algorithm"] != "sha256":
            raise BatchFingerprintError(
                "invalid_content_signature",
                f"{label} signature algorithm must be 'sha256'",
            )
        if "size_bytes" in source and (
            isinstance(source["size_bytes"], bool)
            or not isinstance(source["size_bytes"], int)
            or source["size_bytes"] < 0
        ):
            raise BatchFingerprintError(
                "invalid_content_signature",
                f"{label} size_bytes must be a non-negative integer",
            )
        if "sha256" in source:
            digest = _validate_sha256(source["sha256"], label=label)
            if "digest" in source:
                second_digest = _validate_sha256(
                    source["digest"], label=label
                )
                if digest != second_digest:
                    raise BatchFingerprintError(
                        "invalid_content_signature",
                        f"{label} sha256 and digest fields disagree",
                    )
        elif source.get("algorithm") == "sha256" and "digest" in source:
            digest = _validate_sha256(source["digest"], label=label)
        else:
            raise BatchFingerprintError(
                "invalid_content_signature",
                f"{label} signature must contain sha256 or "
                "algorithm='sha256' plus digest",
            )
    elif isinstance(source, (str, os.PathLike)):
        text = os.fspath(source)
        if isinstance(source, str) and (
            text.startswith("sha256:") or _SHA256_RE.fullmatch(text)
        ):
            digest = _validate_sha256(text, label=label)
        else:
            path = Path(text).expanduser()
            try:
                if not path.is_file():
                    raise BatchFingerprintError(
                        "missing_content_file",
                        f"{label} file does not exist or is not a regular file: "
                        f"{path}",
                    )
                hasher = hashlib.sha256()
                with path.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        hasher.update(chunk)
                digest = hasher.hexdigest()
            except BatchFingerprintError:
                raise
            except OSError as exc:
                raise BatchFingerprintError(
                    "content_read_failed",
                    f"could not read {label} file {path}: {exc}",
                ) from exc
    else:
        raise BatchFingerprintError(
            "invalid_content_source",
            f"{label} must be bytes, a path, or a SHA-256 signature",
        )
    return {"algorithm": "sha256", "digest": digest}


def _normalize_pose_artifacts(pose_artifacts: Any) -> list[dict[str, Any]]:
    if isinstance(pose_artifacts, Mapping) and not (
        "sha256" in pose_artifacts
        or (
            pose_artifacts.get("algorithm") == "sha256"
            and "digest" in pose_artifacts
        )
    ):
        if not pose_artifacts:
            raise BatchFingerprintError(
                "missing_pose_artifacts", "pose_artifacts cannot be empty"
            )
        if any(
            not isinstance(name, str) or not name.strip()
            for name in pose_artifacts
        ):
            raise BatchFingerprintError(
                "invalid_pose_artifact_name",
                "pose artifact names must be non-empty strings",
            )
        records = []
        for name in sorted(pose_artifacts):
            records.append(
                {
                    "name": name,
                    "content": artifact_signature(
                        pose_artifacts[name], label=f"pose_artifacts.{name}"
                    ),
                }
            )
        return records
    if isinstance(pose_artifacts, Sequence) and not isinstance(
        pose_artifacts, (str, bytes, bytearray, memoryview)
    ):
        if not pose_artifacts:
            raise BatchFingerprintError(
                "missing_pose_artifacts", "pose_artifacts cannot be empty"
            )
        return [
            {
                "name": f"pose_{index:06d}",
                "content": artifact_signature(
                    source, label=f"pose_artifacts[{index}]"
                ),
            }
            for index, source in enumerate(pose_artifacts)
        ]
    return [
        {
            "name": "pose_batch",
            "content": artifact_signature(
                pose_artifacts, label="pose_artifacts"
            ),
        }
    ]


def _normalize_tool_matrix(enabled_tools: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(enabled_tools, Mapping):
        raise BatchFingerprintError(
            "invalid_tool_matrix", "enabled_tools must be a mapping"
        )
    normalized: dict[str, Any] = {}
    for tool_id, config in enabled_tools.items():
        if not isinstance(tool_id, str) or not tool_id.strip():
            raise BatchFingerprintError(
                "invalid_tool_id", "tool ids must be non-empty strings"
            )
        if isinstance(config, bool):
            normalized[tool_id] = {"enabled": config}
            continue
        if not isinstance(config, Mapping):
            raise BatchFingerprintError(
                "invalid_tool_config",
                f"enabled_tools.{tool_id} must be a bool or mapping",
            )
        if not isinstance(config.get("enabled"), bool):
            raise BatchFingerprintError(
                "missing_tool_enabled",
                f"enabled_tools.{tool_id}.enabled must be a boolean",
            )
        normalized[tool_id] = _normalize_json(
            config, path=f"$.enabled_tools.{tool_id}"
        )
    return {tool_id: normalized[tool_id] for tool_id in sorted(normalized)}


def build_batch_fingerprint(
    *,
    receptor_artifact: Any,
    pose_artifacts: Any,
    enabled_tools: Mapping[str, Any],
    params: Mapping[str, Any] | None = None,
    runtime_metadata: Mapping[str, Any] | None = None,
    model_metadata: Mapping[str, Any] | None = None,
    checkpoint_metadata: Mapping[str, Any] | None = None,
    version_metadata: Mapping[str, Any] | None = None,
    extra_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Fingerprint the complete scientific execution contract for one batch.

    Receptor and pose inputs are content-addressed, never path-addressed.
    Callers should populate the metadata groups with every enabled tool's
    resolved runtime, model/checkpoint identity, and version.  Empty groups are
    represented explicitly so omissions are visible in stored metadata.
    """
    payload = {
        "schema_version": BATCH_FINGERPRINT_SCHEMA,
        "receptor": artifact_signature(
            receptor_artifact, label="receptor_artifact"
        ),
        "poses": _normalize_pose_artifacts(pose_artifacts),
        "enabled_tools": _normalize_tool_matrix(enabled_tools),
        "params": _normalize_optional_mapping(params, path="$.params"),
        "runtime_metadata": _normalize_optional_mapping(
            runtime_metadata, path="$.runtime_metadata"
        ),
        "model_metadata": _normalize_optional_mapping(
            model_metadata, path="$.model_metadata"
        ),
        "checkpoint_metadata": _normalize_optional_mapping(
            checkpoint_metadata, path="$.checkpoint_metadata"
        ),
        "version_metadata": _normalize_optional_mapping(
            version_metadata, path="$.version_metadata"
        ),
        "extra_metadata": _normalize_optional_mapping(
            extra_metadata, path="$.extra_metadata"
        ),
    }
    return {
        "status": "ok",
        "schema_version": BATCH_FINGERPRINT_SCHEMA,
        "algorithm": "sha256",
        "fingerprint": stable_json_hash(
            payload, namespace=BATCH_FINGERPRINT_SCHEMA
        ),
        "payload": payload,
    }


def build_batch_output_digest(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Seal the complete, ordered JSON result rows emitted by one batch.

    The digest detects accidental or post-run modification of evaluator output.
    It is deliberately separate from the input fingerprint: callers still have
    to bind every cached row to the current input cohort before reuse.
    """
    if (
        not isinstance(records, Sequence)
        or isinstance(records, (str, bytes, bytearray))
    ):
        raise BatchFingerprintError(
            "invalid_batch_output",
            "batch output records must be a sequence of mappings",
        )
    normalized_records = []
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise BatchFingerprintError(
                "invalid_batch_output",
                f"batch output record {index} must be a mapping",
            )
        normalized_records.append(
            _normalize_json(record, path=f"$.records[{index}]")
        )
    payload = {
        "schema_version": BATCH_OUTPUT_DIGEST_SCHEMA,
        "records": normalized_records,
    }
    return {
        "schema_version": BATCH_OUTPUT_DIGEST_SCHEMA,
        "algorithm": "sha256",
        "record_count": len(normalized_records),
        "digest": stable_json_hash(
            payload, namespace=BATCH_OUTPUT_DIGEST_SCHEMA
        ),
    }


def validate_batch_output_digest(
    metadata: Any,
    records: Any,
) -> dict[str, Any]:
    """Verify the stored seal over the complete ordered batch evaluation."""
    if not isinstance(metadata, Mapping):
        return {
            "reusable": False,
            "status": "invalid_metadata",
            "reason": "batch metadata must be a mapping",
        }
    stored = metadata.get("batch_output_digest")
    if not isinstance(stored, Mapping):
        return {
            "reusable": False,
            "status": "missing_batch_output_digest",
            "reason": "batch metadata has no supported batch_output_digest",
        }
    if stored.get("schema_version") != BATCH_OUTPUT_DIGEST_SCHEMA:
        return {
            "reusable": False,
            "status": "unsupported_batch_output_digest_schema",
            "reason": "batch output digest uses an unsupported schema",
        }
    if stored.get("algorithm") != "sha256":
        return {
            "reusable": False,
            "status": "unsupported_batch_output_digest_algorithm",
            "reason": "batch output digest must use sha256",
        }
    try:
        stored_digest = _fingerprint_value(
            stored.get("digest"), label="batch_output_digest.digest"
        )
        expected = build_batch_output_digest(records)
        expected_digest = _fingerprint_value(
            expected["digest"], label="expected_batch_output_digest.digest"
        )
        stored_count = stored.get("record_count")
        if (
            isinstance(stored_count, bool)
            or not isinstance(stored_count, int)
            or stored_count < 0
        ):
            raise BatchFingerprintError(
                "invalid_batch_output_record_count",
                "batch output digest record_count must be a non-negative integer",
            )
    except BatchFingerprintError as exc:
        return {
            "reusable": False,
            "status": exc.code,
            "reason": str(exc),
        }
    if stored_count != expected["record_count"]:
        return {
            "reusable": False,
            "status": "batch_output_record_count_mismatch",
            "reason": "stored batch output record count differs from evaluation",
            "stored_record_count": stored_count,
            "actual_record_count": expected["record_count"],
        }
    if not hmac.compare_digest(stored_digest, expected_digest):
        return {
            "reusable": False,
            "status": "batch_output_digest_mismatch",
            "reason": "stored batch evaluation has changed since it was written",
            "stored_digest": stored_digest,
            "actual_digest": expected_digest,
        }
    return {
        "reusable": True,
        "status": "match",
        "record_count": stored_count,
        "actual_digest": expected_digest,
        "stored_digest": stored_digest,
    }


def _fingerprint_value(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise BatchFingerprintError(
            "invalid_batch_fingerprint",
            f"{label} fingerprint must be a string",
        )
    digest = _validate_sha256(value, label=label)
    return f"sha256:{digest}"


def _validate_fingerprint_record(
    record: Any,
    *,
    label: str,
    allow_string: bool,
    outer_schema: Any = None,
) -> tuple[str, str | None]:
    if isinstance(record, str):
        if not allow_string and outer_schema != BATCH_FINGERPRINT_SCHEMA:
            raise BatchFingerprintError(
                "missing_fingerprint_schema",
                f"{label} string fingerprint has no supported schema",
            )
        return _fingerprint_value(record, label=label), outer_schema
    if not isinstance(record, Mapping):
        raise BatchFingerprintError(
            "invalid_batch_fingerprint",
            f"{label} must be a fingerprint record",
        )
    schema = record.get("schema_version")
    if schema != BATCH_FINGERPRINT_SCHEMA:
        raise BatchFingerprintError(
            "unsupported_fingerprint_schema",
            f"{label} uses unsupported schema {schema!r}",
        )
    fingerprint = _fingerprint_value(
        record.get("fingerprint"), label=f"{label}.fingerprint"
    )
    if "payload" in record:
        payload = record["payload"]
        if not isinstance(payload, Mapping):
            raise BatchFingerprintError(
                "invalid_fingerprint_payload",
                f"{label}.payload must be a mapping",
            )
        if payload.get("schema_version") != BATCH_FINGERPRINT_SCHEMA:
            raise BatchFingerprintError(
                "unsupported_fingerprint_schema",
                f"{label}.payload has an unsupported schema",
            )
        recomputed = stable_json_hash(
            payload, namespace=BATCH_FINGERPRINT_SCHEMA
        )
        if not hmac.compare_digest(fingerprint, recomputed):
            raise BatchFingerprintError(
                "corrupt_fingerprint_payload",
                f"{label}.payload does not match its stored fingerprint",
            )
    return fingerprint, schema


def validate_batch_metadata(
    metadata: Any,
    *,
    expected_fingerprint: Mapping[str, Any] | str,
) -> dict[str, Any]:
    """Decide whether a completed batch report is safe to resume.

    Old reports without a fingerprint, failed reports, unsupported schemas,
    corrupted embedded payloads, and mismatches are all explicit non-reusable
    states.  The function never upgrades a legacy success-only report.
    """
    if not isinstance(metadata, Mapping):
        return {
            "reusable": False,
            "status": "invalid_metadata",
            "reason": "batch metadata must be a mapping",
        }
    if metadata.get("success") is not True:
        return {
            "reusable": False,
            "status": "batch_not_successful",
            "reason": "only a batch with success=true can be resumed",
        }
    if "batch_fingerprint" not in metadata:
        return {
            "reusable": False,
            "status": "missing_fingerprint",
            "reason": "legacy batch metadata has no batch_fingerprint",
        }

    try:
        expected, _ = _validate_fingerprint_record(
            expected_fingerprint,
            label="expected_fingerprint",
            allow_string=True,
        )
        actual, _ = _validate_fingerprint_record(
            metadata["batch_fingerprint"],
            label="metadata.batch_fingerprint",
            allow_string=False,
            outer_schema=metadata.get("batch_fingerprint_schema"),
        )
    except BatchFingerprintError as exc:
        return {
            "reusable": False,
            "status": exc.code,
            "reason": str(exc),
        }

    if not hmac.compare_digest(actual, expected):
        return {
            "reusable": False,
            "status": "fingerprint_mismatch",
            "reason": "stored batch inputs or execution metadata differ",
            "actual_fingerprint": actual,
            "expected_fingerprint": expected,
        }
    return {
        "reusable": True,
        "status": "match",
        "actual_fingerprint": actual,
        "expected_fingerprint": expected,
    }

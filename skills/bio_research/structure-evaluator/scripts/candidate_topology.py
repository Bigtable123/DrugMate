"""Candidate-specific Amber topology preparation for expensive Step 4 routes.

This module deliberately sits behind an opt-in evaluator tool.  The ordinary
pose-QC path must not pay the cost of charge assignment or topology generation.
When enabled, the caller should pass only QC-surviving, already-ranked
candidates (normally a top-K SDF).

The public :func:`build_candidate_topologies` entry point has two modes:

``preflight``
    Validate chemical identity, force-field applicability and pose lineage,
    then write a manifest without invoking AmberTools.

``build``
    Perform the same all-candidate preflight first, then generate a separate
    ff14SB/GAFF2 + AM1-BCC complex topology for every evaluable SDF record.

Known unsupported chemistry is represented as ``not_evaluable`` rather than a
false pass.  Tool execution failures are represented as ``failed``.  A written
manifest therefore survives partial batches and provides the evidence needed by
future MM-GBSA, short-MD and free-energy stages.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

from pose_identity import PoseIdentityError, build_pose_identity
from tool_runtime import ToolRuntimeError, run_tool


SCHEMA_VERSION = "candidate-topologies/v1"
SUPPORTED_LIGAND_ELEMENTS = frozenset(
    {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
)
METAL_ELEMENTS = frozenset(
    {
        "Li", "Na", "K", "Rb", "Cs",
        "Be", "Mg", "Ca", "Sr", "Ba",
        "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu",
    }
)
WATER_RESNAMES = frozenset({"HOH", "WAT", "H2O", "TIP3", "TIP3P"})


class CandidateTopologyError(RuntimeError):
    """Raised for malformed API calls, not for per-candidate chemistry failures."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _safe_component(value: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    return (cleaned[:80] or fallback)


def _first_property(mol: Any, names: Iterable[str]) -> str | None:
    for name in names:
        try:
            if mol.HasProp(name):
                value = str(mol.GetProp(name)).strip()
                if value:
                    return value
        except Exception:
            continue
    return None


def _integer_property(mol: Any, names: Iterable[str], default: int) -> int:
    value = _first_property(mol, names)
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        raise PoseIdentityError(
            "invalid_coordinate_revision",
            f"coordinate revision is not an integer: {value!r}",
        )


def _parse_source_identity(mol: Any) -> Any:
    raw = _first_property(mol, ("source_identity", "SOURCE_IDENTITY"))
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except (TypeError, ValueError):
        return raw


def _parse_lineage_metadata(mol: Any) -> dict:
    raw = _first_property(mol, ("lineage_metadata", "LINEAGE_METADATA"))
    if raw is None:
        return {}
    try:
        value = json.loads(raw)
    except (TypeError, ValueError) as exc:
        raise PoseIdentityError(
            "invalid_lineage_metadata",
            "lineage_metadata SDF property must be a JSON object",
        ) from exc
    if not isinstance(value, dict):
        raise PoseIdentityError(
            "invalid_lineage_metadata",
            "lineage_metadata SDF property must be a JSON object",
        )
    return value


def _record_identity(
    mol: Any,
    index: int,
    *,
    source_artifact_sha256: str,
) -> dict:
    compound_id = _first_property(
        mol, ("compound_id", "CompoundID", "compound_name", "_Name")
    )
    microstate_id = _first_property(
        mol, ("microstate_id", "MICROSTATE_ID", "protonation_state_id")
    )
    variant_id = _first_property(mol, ("variant_id", "VARIANT_ID"))
    declared_chemical_id = _first_property(
        mol, ("chemical_id", "CHEMICAL_ID")
    )
    declared_source_identity_id = _first_property(
        mol, ("source_identity_id", "SOURCE_IDENTITY_ID")
    )
    declared_pose_id = _first_property(mol, ("pose_id", "POSE_ID"))
    parent_pose_id = _first_property(
        mol, ("parent_pose_id", "PARENT_POSE_ID")
    )
    coordinate_revision = _integer_property(
        mol, ("coordinate_revision", "COORDINATE_REVISION"), 0
    )
    source_identity = _parse_source_identity(mol)
    if source_identity is None:
        source_identity = {
            "input_artifact_sha256": source_artifact_sha256,
            "record_index": index,
            "compound_id": compound_id,
            "microstate_id": microstate_id,
            "variant_id": variant_id,
        }
    pose_identity = build_pose_identity(
        mol,
        source_identity=source_identity,
        parent_pose_id=parent_pose_id,
        coordinate_revision=coordinate_revision,
        lineage_metadata=_parse_lineage_metadata(mol),
    )
    chemical_identity = pose_identity["chemical_identity"]

    # A canonical upstream identifier must describe exactly this graph/source/
    # coordinate lineage. Arbitrary legacy IDs remain audit-only metadata.
    comparisons = (
        ("chemical_id", declared_chemical_id, chemical_identity["chemical_id"], "chem:"),
        (
            "source_identity_id",
            declared_source_identity_id,
            chemical_identity["source_identity_id"],
            "source:",
        ),
        ("pose_id", declared_pose_id, pose_identity["pose_id"], "pose:"),
    )
    for label, declared, computed, prefix in comparisons:
        if declared and declared.startswith(prefix) and declared != computed:
            raise PoseIdentityError(
                "lineage_identity_mismatch",
                f"declared canonical {label} does not match recomputed identity",
                details={
                    "field": label,
                    "declared": declared,
                    "computed": computed,
                    "record_index": index,
                },
            )
    return {
        "identity_schema_version": pose_identity["schema_version"],
        "chemical_identity_id": chemical_identity["identity_id"],
        "chemical_id": chemical_identity["chemical_id"],
        "source_identity_id": chemical_identity["source_identity_id"],
        "source_identity": chemical_identity["source_identity"],
        "pose_id": pose_identity["pose_id"],
        "parent_pose_id": pose_identity["parent_pose_id"],
        "coordinate_revision": pose_identity["coordinate_revision"],
        "coordinate_signature": pose_identity["coordinates"]["signature"],
        "input_declared_chemical_id": declared_chemical_id,
        "input_declared_source_identity_id": declared_source_identity_id,
        "input_declared_pose_id": declared_pose_id,
        # Preserve legacy generator/docking identifiers while the canonical
        # identity fields above remain the downstream contract.
        "compound_id": compound_id,
        "microstate_id": microstate_id,
        "variant_id": variant_id,
        "record_index": index,
        "canonical_isomeric_smiles": chemical_identity[
            "canonical_isomeric_smiles"
        ],
        "formal_charge": chemical_identity["formal_charge"],
    }


def _load_candidates(candidate_sdf: Path, max_candidates: int | None) -> list[dict]:
    """Load and sanitize SDF records without silently dropping parse failures."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise CandidateTopologyError(
            "RDKit is required to preflight candidate topologies"
        ) from exc

    supplier = Chem.SDMolSupplier(
        str(candidate_sdf), removeHs=False, sanitize=False, strictParsing=True
    )
    source_artifact_sha256 = _sha256(candidate_sdf)
    records: list[dict] = []
    for index, mol in enumerate(supplier):
        if max_candidates is not None and index >= max_candidates:
            break
        if mol is None:
            records.append({
                "index": index,
                "mol": None,
                "status": "not_evaluable",
                "reason_code": "sdf_parse_failed",
                "reason": "RDKit could not parse this SDF record",
                "identity": {
                    "record_index": index,
                    "pose_id": f"unparsed:pose:{index}",
                },
            })
            continue

        try:
            Chem.SanitizeMol(mol)
            identity = _record_identity(
                mol,
                index,
                source_artifact_sha256=source_artifact_sha256,
            )
        except PoseIdentityError as exc:
            records.append({
                "index": index,
                "mol": mol,
                "status": "not_evaluable",
                "reason_code": exc.code,
                "reason": str(exc),
                "identity": {
                    "record_index": index,
                    "pose_id": _first_property(
                        mol, ("pose_id", "POSE_ID", "_Name")
                    ) or f"identity-error:pose:{index}",
                    "identity_error": exc.as_result()["error"],
                },
            })
            continue
        except Exception as exc:
            records.append({
                "index": index,
                "mol": mol,
                "status": "not_evaluable",
                "reason_code": "chemical_sanitization_failed",
                "reason": str(exc),
                "identity": {
                    "record_index": index,
                    "pose_id": _first_property(
                        mol, ("pose_id", "POSE_ID", "_Name")
                    ) or f"unsanitized:pose:{index}",
                },
            })
            continue

        elements = sorted({atom.GetSymbol() for atom in mol.GetAtoms()})
        unsupported = sorted(set(elements) - SUPPORTED_LIGAND_ELEMENTS)
        fragments = tuple(Chem.GetMolFrags(mol))
        radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        has_3d = bool(mol.GetNumConformers()) and bool(mol.GetConformer().Is3D())

        record = {
            "index": index,
            "mol": mol,
            "status": "passed",
            "reason_code": None,
            "reason": None,
            "identity": identity,
            "preflight": {
                "elements": elements,
                "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
                "fragment_count": len(fragments),
                "radical_electron_count": int(radicals),
                "has_3d_coordinates": has_3d,
            },
        }
        if unsupported:
            record.update({
                "status": "not_evaluable",
                "reason_code": "unsupported_ligand_elements",
                "reason": (
                    "Generic GAFF2 route does not support elements: "
                    + ", ".join(unsupported)
                ),
            })
        elif len(fragments) != 1:
            record.update({
                "status": "not_evaluable",
                "reason_code": "disconnected_ligand",
                "reason": (
                    "Candidate contains multiple disconnected fragments; "
                    "desalt/select a single chemical entity first"
                ),
            })
        elif radicals:
            record.update({
                "status": "not_evaluable",
                "reason_code": "radical_ligand",
                "reason": "Generic GAFF2/AM1-BCC route does not accept radicals",
            })
        elif not has_3d:
            record.update({
                "status": "not_evaluable",
                "reason_code": "missing_3d_coordinates",
                "reason": "Candidate has no declared 3D conformer",
            })
        records.append(record)
    return records


def inspect_receptor_risks(receptor_pdb: str | Path) -> dict:
    """Return conservative receptor blockers for the generic Amber route."""
    path = Path(receptor_pdb)
    metal_atoms: list[dict] = []
    unsupported_hetero_residues: set[str] = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        element = line[76:78].strip().title()
        if not element:
            atom_name = re.sub(r"[^A-Za-z]", "", line[12:16])
            element = atom_name[:2].title() if len(atom_name) >= 2 else atom_name.title()
        if element in METAL_ELEMENTS or resname.title() in METAL_ELEMENTS:
            metal_atoms.append({
                "serial": line[6:11].strip(),
                "atom_name": line[12:16].strip(),
                "resname": resname,
                "chain": line[21:22].strip(),
                "residue_number": line[22:26].strip(),
                "element": element,
            })
        elif resname and resname not in WATER_RESNAMES:
            unsupported_hetero_residues.add(resname)
    return {
        "metal_atoms": metal_atoms,
        "unsupported_hetero_residues": sorted(unsupported_hetero_residues),
    }


def _command_record(result: Mapping[str, Any]) -> dict:
    """Keep reproducibility fields while bounding potentially huge tool logs."""
    return {
        "command": list(result.get("command") or []),
        "environment": result.get("environment"),
        "returncode": result.get("returncode"),
        "elapsed_seconds": result.get("elapsed_seconds"),
        "stdout_tail": str(result.get("stdout") or "")[-4000:],
        "stderr_tail": str(result.get("stderr") or "")[-4000:],
    }


def _extract_version(text: str) -> str | None:
    patterns = (
        r"(?i)\bversion\s*[:=]?\s*(v?\d+(?:\.\d+){0,3}[A-Za-z0-9_.+-]*)",
        r"(?i)\b(?:ambertools|antechamber|parmchk2|tleap|pdb4amber)\s+"
        r"(v?\d+(?:\.\d+){0,3}[A-Za-z0-9_.+-]*)",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    return None


def collect_tool_versions(runtime: Mapping[str, Any]) -> dict[str, dict]:
    """Best-effort AmberTools version banners with exact probe commands."""
    commands = dict(runtime.get("commands") or {})
    versions: dict[str, dict] = {}

    # Amber executables do not expose a consistent --version interface. The
    # conda package record is the authoritative suite version for this
    # registered environment, while per-command banners below remain useful
    # diagnostics.
    package_code = (
        "import glob,json,os,sys;"
        "paths=glob.glob(os.path.join(sys.prefix,'conda-meta','ambertools-*.json'));"
        "data=json.load(open(sorted(paths)[-1])) if paths else {};"
        "print(json.dumps({'name':data.get('name'),'version':data.get('version'),"
        "'build':data.get('build'),'channel':data.get('channel')}))"
    )
    try:
        package_result = run_tool(
            dict(runtime),
            ["-c", package_code],
            command="python",
            timeout=15,
            check=False,
        )
        package = json.loads(str(package_result.get("stdout") or "{}").strip() or "{}")
        versions["ambertools"] = {
            "version": package.get("version") or "unknown",
            "build": package.get("build"),
            "channel": package.get("channel"),
            "probe": _command_record(package_result),
        }
    except Exception as exc:
        versions["ambertools"] = {
            "version": "unknown",
            "error": str(exc),
            "probe": {
                "command": ["python", "-c", "<conda ambertools package probe>"],
                "environment": runtime.get("environment"),
            },
        }

    probes = {
        "pdb4amber": ["--version"],
        "antechamber": ["-h"],
        "parmchk2": ["-h"],
        "tleap": ["-h"],
    }
    for name, args in probes.items():
        executable = commands.get(name, name)
        try:
            result = run_tool(
                dict(runtime),
                args,
                command=executable,
                timeout=15,
                check=False,
            )
            text = f"{result.get('stdout', '')}\n{result.get('stderr', '')}"
            nonempty = [line.strip() for line in text.splitlines() if line.strip()]
            versions[name] = {
                "version": _extract_version(text) or "unknown",
                "banner": "\n".join(nonempty[:8])[:1200],
                "probe": _command_record(result),
            }
        except Exception as exc:
            versions[name] = {
                "version": "unknown",
                "error": str(exc),
                "probe": {
                    "command": [str(executable), *args],
                    "environment": runtime.get("environment"),
                },
            }
    return versions


def _write_single_sdf(mol: Any, path: Path) -> None:
    from rdkit import Chem

    writer = Chem.SDWriter(str(path))
    try:
        writer.write(mol)
    finally:
        writer.close()
    if not path.is_file() or path.stat().st_size == 0:
        raise CandidateTopologyError(f"failed to write candidate SDF: {path}")


def _validate_mol2_charge(path: Path, expected_charge: int) -> dict:
    """Validate that written AM1-BCC atom charges sum to the requested net charge."""
    if not path.is_file() or path.stat().st_size == 0:
        raise CandidateTopologyError(
            "antechamber did not create a non-empty ligand MOL2"
        )
    in_atom_section = False
    charges: list[float] = []
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.upper() == "@<TRIPOS>ATOM":
            in_atom_section = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_atom_section:
            break
        if not in_atom_section or not stripped:
            continue
        fields = stripped.split()
        if len(fields) < 9:
            raise CandidateTopologyError(
                "antechamber MOL2 atom record has no partial-charge field"
            )
        try:
            charges.append(float(fields[-1]))
        except ValueError as exc:
            raise CandidateTopologyError(
                "antechamber MOL2 contains a non-numeric partial charge"
            ) from exc
    if not charges:
        raise CandidateTopologyError(
            "antechamber MOL2 contains no atom partial charges"
        )
    observed = sum(charges)
    delta = observed - float(expected_charge)
    result = {
        "atom_count": len(charges),
        "expected_formal_charge": int(expected_charge),
        "mol2_partial_charge_sum": round(observed, 6),
        "delta": round(delta, 6),
        "tolerance": 0.01,
    }
    if abs(delta) > result["tolerance"]:
        raise CandidateTopologyError(
            "AM1-BCC partial charges do not sum to the requested formal charge "
            f"({observed:.6f} versus {expected_charge:+d})"
        )
    return result


def _leap_path(path: Path) -> str:
    # LEaP accepts double-quoted paths. JSON string escaping is compatible for
    # ordinary filesystem paths and avoids unquoted whitespace/metacharacters.
    return json.dumps(str(path.resolve()))


def _run_candidate_build(
    receptor_cleaned: Path,
    record: dict,
    output_dir: Path,
    runtime: Mapping[str, Any],
    *,
    protein_force_field: str,
    timeout_seconds: int,
    solvate_for_md: bool = False,
    water_model: str = "tip3p",
    solvent_padding_angstrom: float = 10.0,
) -> dict:
    commands = dict(runtime.get("commands") or {})
    index = int(record["index"])
    identity = record["identity"]
    pose_id = str(identity.get("pose_id") or f"pose-{index}")
    candidate_dir = output_dir / (
        f"candidate_{index:05d}_{_safe_component(pose_id, f'pose-{index}')}"
    )
    candidate_dir.mkdir(parents=True, exist_ok=True)

    candidate_sdf = candidate_dir / "candidate.sdf"
    ligand_mol2 = candidate_dir / "ligand.gaff2.mol2"
    ligand_frcmod = candidate_dir / "ligand.gaff2.frcmod"
    complex_prmtop = candidate_dir / "complex.prmtop"
    complex_inpcrd = candidate_dir / "complex.inpcrd"
    receptor_prmtop = candidate_dir / "receptor.prmtop"
    receptor_inpcrd = candidate_dir / "receptor.inpcrd"
    ligand_prmtop = candidate_dir / "ligand.prmtop"
    ligand_inpcrd = candidate_dir / "ligand.inpcrd"
    solvated_complex_prmtop = candidate_dir / "solvated_complex.prmtop"
    solvated_complex_inpcrd = candidate_dir / "solvated_complex.inpcrd"
    leap_input = candidate_dir / "tleap.in"

    _write_single_sdf(record["mol"], candidate_sdf)
    formal_charge = int(identity["formal_charge"])
    command_records: list[dict] = []
    charge_validation: dict = {}

    try:
        result = run_tool(
            dict(runtime),
            [
                "-i", candidate_sdf, "-fi", "sdf",
                "-o", ligand_mol2, "-fo", "mol2",
                "-at", "gaff2", "-c", "bcc", "-nc", formal_charge,
                "-rn", "LIG", "-s", "2",
            ],
            command=commands.get("antechamber", "antechamber"),
            cwd=candidate_dir,
            timeout=timeout_seconds,
        )
        command_records.append(_command_record(result))
        charge_validation = _validate_mol2_charge(
            ligand_mol2, formal_charge
        )
        result = run_tool(
            dict(runtime),
            [
                "-i", ligand_mol2, "-f", "mol2",
                "-o", ligand_frcmod, "-s", "gaff2",
            ],
            command=commands.get("parmchk2", "parmchk2"),
            cwd=candidate_dir,
            timeout=min(timeout_seconds, 900),
        )
        command_records.append(_command_record(result))

        leap_lines = [
            f"source leaprc.protein.{protein_force_field}",
            "source leaprc.gaff2",
        ]
        if solvate_for_md:
            if water_model != "tip3p":
                raise CandidateTopologyError(
                    "short-MD topology currently supports only TIP3P water"
                )
            leap_lines.append("source leaprc.water.tip3p")
        leap_lines.extend([
            f"REC = loadpdb {_leap_path(receptor_cleaned)}",
            f"LIG = loadmol2 {_leap_path(ligand_mol2)}",
            f"loadamberparams {_leap_path(ligand_frcmod)}",
            "COM = combine { REC LIG }",
            f"saveamberparm COM {_leap_path(complex_prmtop)} {_leap_path(complex_inpcrd)}",
            f"saveamberparm REC {_leap_path(receptor_prmtop)} {_leap_path(receptor_inpcrd)}",
            f"saveamberparm LIG {_leap_path(ligand_prmtop)} {_leap_path(ligand_inpcrd)}",
        ])
        if solvate_for_md:
            leap_lines.extend([
                f"solvateOct COM TIP3PBOX {float(solvent_padding_angstrom):.3f}",
                "addIonsRand COM Na+ 0",
                "addIonsRand COM Cl- 0",
                f"saveamberparm COM {_leap_path(solvated_complex_prmtop)} {_leap_path(solvated_complex_inpcrd)}",
            ])
        leap_lines.append("quit")
        leap_input.write_text("\n".join(leap_lines) + "\n", encoding="utf-8")
        result = run_tool(
            dict(runtime),
            ["-f", leap_input],
            command=commands.get("tleap", "tleap"),
            cwd=candidate_dir,
            timeout=min(timeout_seconds, 900),
        )
        command_records.append(_command_record(result))
        leap_text = f"{result.get('stdout', '')}\n{result.get('stderr', '')}"
        leap_errors = re.findall(r"(?i)\bErrors\s*=\s*(\d+)", leap_text)
        if leap_errors and any(int(count) != 0 for count in leap_errors):
            raise CandidateTopologyError(
                "tleap reported one or more topology errors"
            )
    except (OSError, ToolRuntimeError, CandidateTopologyError) as exc:
        failed_result = getattr(exc, "result", None)
        if failed_result:
            command_records.append(_command_record(failed_result))
        return {
            "status": "failed",
            "reason_code": "ambertools_execution_failed",
            "reason": str(exc),
            "commands": command_records,
            "charge_validation": charge_validation,
            "output_files": {},
        }

    outputs = {
        "candidate_sdf": candidate_sdf,
        "ligand_mol2": ligand_mol2,
        "ligand_frcmod": ligand_frcmod,
        "complex_prmtop": complex_prmtop,
        "complex_inpcrd": complex_inpcrd,
        "receptor_prmtop": receptor_prmtop,
        "receptor_inpcrd": receptor_inpcrd,
        "ligand_prmtop": ligand_prmtop,
        "ligand_inpcrd": ligand_inpcrd,
        "tleap_input": leap_input,
    }
    if solvate_for_md:
        outputs.update({
            "solvated_complex_prmtop": solvated_complex_prmtop,
            "solvated_complex_inpcrd": solvated_complex_inpcrd,
        })
    missing = [
        name
        for name, path in outputs.items()
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        return {
            "status": "failed",
            "reason_code": "ambertools_outputs_missing",
            "reason": "AmberTools did not create: " + ", ".join(missing),
            "commands": command_records,
            "charge_validation": charge_validation,
            "output_files": {
                name: str(path) for name, path in outputs.items() if path.is_file()
            },
        }
    return {
        "status": "passed",
        "reason_code": None,
        "reason": None,
        "commands": command_records,
        "charge_validation": charge_validation,
        "output_files": {name: str(path) for name, path in outputs.items()},
        "output_sha256": {name: _sha256(path) for name, path in outputs.items()},
    }


def _overall_status(candidates: list[dict]) -> str:
    statuses = [entry.get("status") for entry in candidates]
    if not statuses or all(status == "not_evaluable" for status in statuses):
        return "not_evaluable"
    if any(status == "failed" for status in statuses):
        return "failed"
    if any(status == "not_evaluable" for status in statuses):
        # A requested batch is complete only when every requested candidate has
        # an evaluable topology. The per-candidate records still permit callers
        # to consume the successful subset explicitly.
        return "failed"
    return "passed"


def build_candidate_topologies(
    receptor_pdb: str | Path,
    candidate_sdf: str | Path,
    output_dir: str | Path,
    runtime: Mapping[str, Any],
    *,
    execution_mode: str = "build",
    protein_force_field: str = "ff14SB",
    ligand_force_field: str = "gaff2",
    ligand_charge_method: str = "bcc",
    max_candidates: int | None = 20,
    timeout_per_candidate_seconds: int = 2400,
    is_covalent: bool = False,
    has_metal: bool = False,
    allow_receptor_hetero: bool = False,
    fail_fast: bool = True,
    solvate_for_md: bool = False,
    water_model: str = "tip3p",
    solvent_padding_angstrom: float = 10.0,
) -> dict:
    """Preflight or build candidate-specific Amber topologies.

    ``success`` means the auditable manifest was written. Scientific
    applicability is intentionally separate in ``status`` and is always one of
    ``passed``, ``failed`` or ``not_evaluable``.
    """
    if execution_mode not in {"preflight", "build"}:
        raise CandidateTopologyError(
            "execution_mode must be 'preflight' or 'build'"
        )
    if protein_force_field not in {"ff14SB", "ff19SB"}:
        raise CandidateTopologyError("unsupported protein force field")
    if ligand_force_field != "gaff2" or ligand_charge_method != "bcc":
        raise CandidateTopologyError(
            "candidate topology contract currently supports only GAFF2 + AM1-BCC"
        )
    if max_candidates is not None and int(max_candidates) < 1:
        raise CandidateTopologyError("max_candidates must be >= 1 or None")
    if int(timeout_per_candidate_seconds) < 1:
        raise CandidateTopologyError("timeout_per_candidate_seconds must be >= 1")
    if water_model != "tip3p":
        raise CandidateTopologyError("short-MD topology currently supports only TIP3P")
    if not 6.0 <= float(solvent_padding_angstrom) <= 30.0:
        raise CandidateTopologyError(
            "solvent_padding_angstrom must be between 6 and 30"
        )

    receptor_path = Path(receptor_pdb).resolve()
    candidate_path = Path(candidate_sdf).resolve()
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "candidate_topologies.json"

    missing = [
        str(path) for path in (receptor_path, candidate_path) if not path.is_file()
    ]
    if missing:
        raise CandidateTopologyError(
            "candidate topology inputs not found: " + ", ".join(missing)
        )

    candidates = _load_candidates(
        candidate_path,
        int(max_candidates) if max_candidates is not None else None,
    )
    receptor_risks = inspect_receptor_risks(receptor_path)
    global_blockers: list[dict] = []
    if is_covalent:
        global_blockers.append({
            "code": "covalent_complex",
            "reason": (
                "Generic non-covalent GAFF2 topology generation cannot parameterize "
                "a protein-ligand covalent bond"
            ),
        })
    if has_metal or receptor_risks["metal_atoms"]:
        global_blockers.append({
            "code": "metal_site",
            "reason": (
                "Generic ff14SB/GAFF2 topology generation does not assign a reviewed "
                "metal-center model; use an MCPB.py-specific protocol"
            ),
        })
    if (
        receptor_risks["unsupported_hetero_residues"]
        and not allow_receptor_hetero
    ):
        global_blockers.append({
            "code": "unsupported_receptor_hetero",
            "reason": (
                "Generic route cannot parameterize retained receptor hetero residues: "
                + ", ".join(receptor_risks["unsupported_hetero_residues"])
            ),
        })

    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "success": True,
        "status": "not_evaluable",
        "execution_mode": execution_mode,
        "protocol": {
            "protein_force_field": protein_force_field,
            "ligand_force_field": ligand_force_field,
            "ligand_charge_method": "AM1-BCC",
            "max_candidates": max_candidates,
            "fail_fast": bool(fail_fast),
            "explicit_solvent_for_md": bool(solvate_for_md),
            "water_model": water_model if solvate_for_md else None,
            "solvent_padding_angstrom": (
                float(solvent_padding_angstrom) if solvate_for_md else None
            ),
        },
        "inputs": {
            "receptor_pdb": str(receptor_path),
            "receptor_sha256": _sha256(receptor_path),
            "candidate_sdf": str(candidate_path),
            "candidate_sdf_sha256": _sha256(candidate_path),
        },
        "applicability": {
            "is_covalent": bool(is_covalent),
            "has_metal_declared": bool(has_metal),
            "receptor_risks": receptor_risks,
            "global_blockers": global_blockers,
        },
        "tool_provenance": {
            "environment": runtime.get("environment"),
            "registered_commands": dict(runtime.get("commands") or {}),
            "versions": {},
        },
        "candidates": [],
        "counts": {},
        "manifest_file": str(manifest_path),
    }

    if global_blockers:
        reason = "; ".join(blocker["reason"] for blocker in global_blockers)
        for record in candidates:
            manifest["candidates"].append({
                **record["identity"],
                "status": "not_evaluable",
                "reason_code": "global_applicability_blocker",
                "reason": reason,
                "preflight": record.get("preflight", {}),
                "commands": [],
                "output_files": {},
            })
        manifest["counts"] = {
            "requested": len(candidates),
            "passed": 0,
            "failed": 0,
            "not_evaluable": len(candidates),
        }
        _atomic_json(manifest_path, manifest)
        return manifest

    for record in candidates:
        candidate_entry = {
            **record["identity"],
            "status": record["status"],
            "reason_code": record.get("reason_code"),
            "reason": record.get("reason"),
            "preflight": record.get("preflight", {}),
            "commands": [],
            "output_files": {},
        }
        manifest["candidates"].append(candidate_entry)

    if execution_mode == "preflight":
        manifest["status"] = _overall_status(manifest["candidates"])
        manifest["counts"] = {
            status: sum(
                entry["status"] == status for entry in manifest["candidates"]
            )
            for status in ("passed", "failed", "not_evaluable")
        }
        manifest["counts"]["requested"] = len(manifest["candidates"])
        _atomic_json(manifest_path, manifest)
        return manifest

    if not candidates:
        manifest["status"] = "not_evaluable"
        manifest["counts"] = {
            "requested": 0, "passed": 0, "failed": 0, "not_evaluable": 0,
        }
        _atomic_json(manifest_path, manifest)
        return manifest

    manifest["tool_provenance"]["versions"] = collect_tool_versions(runtime)
    commands = dict(runtime.get("commands") or {})
    receptor_cleaned = destination / "amber_receptor.pdb"
    try:
        receptor_result = run_tool(
            dict(runtime),
            [
                "-i", receptor_path, "-o", receptor_cleaned,
                "--dry", "--nohyd",
            ],
            command=commands.get("pdb4amber", "pdb4amber"),
            cwd=destination,
            timeout=900,
        )
        manifest["tool_provenance"]["receptor_preparation"] = _command_record(
            receptor_result
        )
        if not receptor_cleaned.is_file():
            raise CandidateTopologyError(
                "pdb4amber completed without creating amber_receptor.pdb"
            )
    except (OSError, ToolRuntimeError, CandidateTopologyError) as exc:
        failed_result = getattr(exc, "result", None)
        manifest["tool_provenance"]["receptor_preparation"] = (
            _command_record(failed_result) if failed_result else {
                "error": str(exc),
                "command": [
                    commands.get("pdb4amber", "pdb4amber"),
                    "-i", str(receptor_path), "-o", str(receptor_cleaned),
                    "--dry", "--nohyd",
                ],
                "environment": runtime.get("environment"),
            }
        )
        for entry in manifest["candidates"]:
            if entry["status"] == "passed":
                entry.update({
                    "status": "failed",
                    "reason_code": "receptor_preparation_failed",
                    "reason": str(exc),
                })
        manifest["status"] = _overall_status(manifest["candidates"])
        manifest["counts"] = {
            status: sum(
                entry["status"] == status for entry in manifest["candidates"]
            )
            for status in ("passed", "failed", "not_evaluable")
        }
        manifest["counts"]["requested"] = len(manifest["candidates"])
        _atomic_json(manifest_path, manifest)
        return manifest

    abort_reason: str | None = None
    for record, entry in zip(candidates, manifest["candidates"]):
        if entry["status"] != "passed":
            continue
        if abort_reason is not None:
            entry.update({
                "status": "not_evaluable",
                "reason_code": "aborted_after_candidate_failure",
                "reason": abort_reason,
            })
            continue
        outcome = _run_candidate_build(
            receptor_cleaned,
            record,
            destination,
            runtime,
            protein_force_field=protein_force_field,
            timeout_seconds=int(timeout_per_candidate_seconds),
            solvate_for_md=bool(solvate_for_md),
            water_model=water_model,
            solvent_padding_angstrom=float(solvent_padding_angstrom),
        )
        entry.update(outcome)
        _atomic_json(manifest_path, {
            **manifest,
            "status": _overall_status(manifest["candidates"]),
        })
        if outcome["status"] == "failed" and fail_fast:
            abort_reason = (
                f"Batch stopped after failure of pose {entry.get('pose_id')}: "
                f"{outcome.get('reason')}"
            )

    manifest["status"] = _overall_status(manifest["candidates"])
    manifest["counts"] = {
        status: sum(
            entry["status"] == status for entry in manifest["candidates"]
        )
        for status in ("passed", "failed", "not_evaluable")
    }
    manifest["counts"]["requested"] = len(manifest["candidates"])
    _atomic_json(manifest_path, manifest)
    return manifest


__all__ = [
    "CandidateTopologyError",
    "SCHEMA_VERSION",
    "SUPPORTED_LIGAND_ELEMENTS",
    "build_candidate_topologies",
    "collect_tool_versions",
    "inspect_receptor_risks",
]

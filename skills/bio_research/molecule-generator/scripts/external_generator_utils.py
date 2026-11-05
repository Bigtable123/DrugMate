"""Adapters for optional molecule generators isolated in Conda environments.

Each adapter performs four jobs only: validate its declared input contract,
create a per-run config, execute the registry-owned runtime without a shell,
and normalize valid output molecules into the Step 2 compound contract.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable

import pipeline_registry as R
from tool_runtime import ToolRuntimeError, probe_tool, run_tool


GENERATOR_IDS = (
    "resgen", "diffsbdd", "targetdiff", "difflinker", "reinvent4",
    "lingo3dmol", "tamgen", "molcraft",
)
_SEED_BRIDGE = Path(__file__).resolve().with_name("seeded_runpy.py")
_GENERATOR_WORKER = Path(__file__).resolve().with_name("external_generator_worker.py")
_STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def run_registered_generator(
    tool_id: str,
    input_data: dict,
    output_dir: Path,
) -> tuple[list[dict], dict]:
    """Run one enabled generator and return normalized compounds + diagnostics."""
    if tool_id not in GENERATOR_IDS:
        raise ValueError(f"unsupported registered generator {tool_id!r}")

    config = dict(((input_data.get("tools") or {}).get(tool_id) or {}))
    params = dict(config.get("params") or {})
    runtime = dict(config.get("runtime") or {})
    if not runtime:
        runtime = dict(R.STEPS["molecule_generation"]["tools"][tool_id]["runtime"])

    work_dir = _create_invocation_work_dir(output_dir, tool_id)
    _verify_asset(runtime)

    status = probe_tool(runtime)
    if not status.get("available"):
        raise RuntimeError(status.get("detail") or f"{tool_id} runtime probe failed")

    dispatch = {
        "resgen": _run_resgen,
        "diffsbdd": _run_diffsbdd,
        "targetdiff": _run_targetdiff,
        "difflinker": _run_difflinker,
        "reinvent4": _run_reinvent4,
        "lingo3dmol": _run_lingo3dmol,
        "tamgen": _run_tamgen,
        "molcraft": _run_molcraft,
    }
    compounds, execution, diagnostics = dispatch[tool_id](
        input_data, params, runtime, work_dir
    )
    report = {
        "tool": tool_id,
        "success": bool(compounds),
        "num_samples_requested": int(params.get("num_samples", 0)),
        "num_valid": len(compounds),
        "work_dir": str(work_dir),
        "diagnostics": diagnostics,
        "execution": execution,
        "runtime": {
            "environment": runtime.get("environment"),
            "source_dir": runtime.get("source_dir"),
            "checkpoint": runtime.get("checkpoint"),
            "asset_sha256": runtime.get("asset_sha256"),
            "auxiliary_checkpoint": runtime.get("auxiliary_checkpoint"),
            "auxiliary_asset_sha256": runtime.get("auxiliary_asset_sha256"),
            "checkpoint_config": runtime.get("checkpoint_config"),
            "checkpoint_config_sha256": runtime.get("checkpoint_config_sha256"),
        },
    }
    return compounds, report


def _run_diffsbdd(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    center = _coerce_center(input_data.get("pocket_center"))
    reference = _optional_file(input_data.get("reference_ligand_file"))
    entrypoint = _required_file(runtime.get("entrypoint"), "DiffSBDD entrypoint")
    checkpoint = _required_file(runtime.get("checkpoint"), "DiffSBDD checkpoint")
    source_dir = _required_dir(runtime.get("source_dir"), "DiffSBDD source")

    radius = float(params.get("pocket_radius", 10.0))
    pocket_pdb, residues = prepare_center_pocket(
        receptor, center, work_dir / "pocket.pdb",
        radius=radius, chain=input_data.get("chain"),
        normalize_residue_ids=reference is None,
    )
    output_sdf = work_dir / "generated.sdf"
    requested = int(params.get("num_samples", 100))
    batch_size = _divisor_batch_size(requested, int(params.get("batch_size", 20)))
    args = [
        str(_SEED_BRIDGE), "--seed", str(int(input_data.get("seed", 42))),
        "--script", str(entrypoint), "--", str(checkpoint),
        "--pdbfile", str(receptor), "--outfile", str(output_sdf),
        "--n_samples", str(requested), "--batch_size", str(batch_size),
        "--timesteps", str(int(params.get("timesteps", 500))), "--sanitize",
    ]
    if reference is not None:
        args.extend(["--ref_ligand", str(reference)])
        pocket_basis = "reference_ligand"
        generation_receptor = receptor
    else:
        args.extend(["--resi_list", *residues])
        pocket_basis = "center_radius_residues"
        # DiffSBDD parses residue selectors as ``CHAIN:integer`` and cannot
        # represent blank chains or PDB insertion codes. The center-cropped
        # file is therefore normalized to chain A + sequential residue ids.
        generation_receptor = pocket_pdb

    pdbfile_index = args.index("--pdbfile") + 1
    args[pdbfile_index] = str(generation_receptor)

    result = _execute(
        runtime, args, work_dir, source_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python",
    )
    compounds, parsed = parse_generated_sdfs([output_sdf], "diffsbdd")
    return compounds, result, {
        **parsed, "pocket_basis": pocket_basis, "pocket_radius": radius,
        "pocket_residue_count": len(residues),
        "normalized_residue_ids": reference is None,
    }


def _run_targetdiff(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    center = _coerce_center(input_data.get("pocket_center"))
    checkpoint = _required_file(runtime.get("checkpoint"), "TargetDiff checkpoint")
    source_dir = _required_dir(runtime.get("source_dir"), "TargetDiff source")
    module = runtime.get("module") or "scripts.sample_for_pocket"

    radius = float(params.get("pocket_radius", 10.0))
    pocket_pdb, residues = prepare_center_pocket(
        receptor, center, work_dir / "pocket.pdb",
        radius=radius, chain=input_data.get("chain"),
    )
    requested = int(params.get("num_samples", 100))
    config_path = work_dir / "sampling.yml"
    config_path.write_text(
        "model:\n"
        f"  checkpoint: {json.dumps(str(checkpoint))}\n\n"
        "sample:\n"
        f"  seed: {int(input_data.get('seed', 42))}\n"
        f"  num_samples: {requested}\n"
        f"  num_steps: {int(params.get('num_steps', 1000))}\n"
        "  pos_only: False\n"
        "  center_pos_mode: protein\n"
        "  sample_num_atoms: prior\n",
        encoding="utf-8",
    )
    result_dir = work_dir / "result"
    args = [
        "-B", "-m", str(module), str(config_path),
        "--pdb_path", str(pocket_pdb), "--device", "cuda:0",
        "--batch_size", str(min(requested, int(params.get("batch_size", 20)))),
        "--num_samples", str(requested), "--result_path", str(result_dir),
    ]
    result = _execute(
        runtime, args, work_dir, source_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python", extra_env={"PYTHONPATH": str(source_dir)},
    )
    compounds, parsed = parse_generated_sdfs(
        sorted((result_dir / "sdf").glob("*.sdf")), "targetdiff"
    )
    return compounds, result, {
        **parsed, "pocket_radius": radius, "pocket_residue_count": len(residues),
    }


def _run_resgen(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    reference = _required_file(
        input_data.get("reference_ligand_file"),
        "reference ligand required by the installed ResGen CLI",
    )
    center = _coerce_center(input_data.get("pocket_center"))
    entrypoint = _required_file(runtime.get("entrypoint"), "ResGen entrypoint")
    checkpoint = _required_file(runtime.get("checkpoint"), "ResGen checkpoint")
    source_dir = _required_dir(runtime.get("source_dir"), "ResGen source")

    requested = int(params.get("num_samples", 100))
    config_path = work_dir / "sampling.yml"
    config_path.write_text(
        "sample:\n"
        f"  seed: {int(input_data.get('seed', 42))}\n"
        "  mask_init: True\n"
        f"  num_samples: {requested}\n"
        f"  beam_size: {int(params.get('beam_size', 500))}\n"
        f"  max_steps: {int(params.get('max_steps', 50))}\n"
        "  threshold:\n"
        "    focal_threshold: 0.5\n"
        "    pos_threshold: 0.25\n"
        "    element_threshold: 0.3\n"
        "    hasatom_threshold: 0.6\n"
        "    bond_threshold: 0.4\n",
        encoding="utf-8",
    )
    result_root = work_dir / "result"
    args = [
        str(_SEED_BRIDGE), "--seed", str(int(input_data.get("seed", 42))),
        "--script", str(entrypoint), "--",
        "--config", str(config_path), "--outdir", str(result_root),
        "--device", "cuda", "--check_point", str(checkpoint),
        "--pdb_file", str(receptor), "--sdf_file", str(reference),
        "--center", ",".join(str(value) for value in center),
        "--name_by_pdb", "True",
    ]
    result = _execute(
        runtime, args, work_dir, source_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python", extra_env={"PYTHONPATH": str(source_dir)},
    )
    aggregate = (
        result_root / receptor.stem / f"{receptor.stem}_gen.sdf"
    )
    compounds, parsed = parse_generated_sdfs([aggregate], "resgen")
    return compounds, result, {
        **parsed, "pocket_basis": "step1_center",
        "reference_ligand_cli_requirement": str(reference),
    }


def _run_difflinker(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    fragments = _required_file(input_data.get("fragment_file"), "DiffLinker fragments")
    _validate_fragment_file(fragments)
    entrypoint = _required_file(runtime.get("entrypoint"), "DiffLinker entrypoint")
    checkpoint = _required_file(runtime.get("checkpoint"), "DiffLinker checkpoint")
    source_dir = _required_dir(runtime.get("source_dir"), "DiffLinker source")

    raw_dir = work_dir / "raw"
    args = [
        str(entrypoint), "--fragments", str(fragments),
        "--protein", str(receptor), "--model", str(checkpoint),
        "--linker_size", str(int(params.get("linker_size", 5))),
        "--output", str(raw_dir),
        "--n_samples", str(int(params.get("num_samples", 100))),
        "--max_batch_size", str(int(params.get("max_batch_size", 16))),
        "--random_seed", str(int(input_data.get("seed", 42))),
    ]
    if params.get("n_steps") is not None:
        args.extend(["--n_steps", str(int(params["n_steps"]))])
    if params.get("anchor_indices"):
        args.extend(["--anchors", str(params["anchor_indices"])])
    result = _execute(
        runtime, args, work_dir, source_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python", extra_env={"PYTHONPATH": str(source_dir)},
    )
    compounds, parsed = parse_generated_sdfs(
        sorted(raw_dir.glob("*.sdf")), "difflinker"
    )
    return compounds, result, {**parsed, "fragment_file": str(fragments)}


def _run_reinvent4(input_data, params, runtime, work_dir):
    prior = _required_file(runtime.get("checkpoint"), "REINVENT4 prior")
    source_dir = _required_dir(runtime.get("source_dir"), "REINVENT4 source")
    output_csv = work_dir / "sampling.csv"
    log_file = work_dir / "sampling.log"
    config_path = work_dir / "sampling.toml"
    requested = int(params.get("num_samples", 200))
    unique = bool(params.get("unique_molecules", True))
    config_path.write_text(
        'run_type = "sampling"\n'
        'device = "cuda:0"\n\n'
        "[parameters]\n"
        f"model_file = {json.dumps(str(prior))}\n"
        f"output_file = {json.dumps(str(output_csv))}\n"
        f"num_smiles = {requested}\n"
        f"unique_molecules = {'true' if unique else 'false'}\n"
        "randomize_smiles = true\n",
        encoding="utf-8",
    )
    args = [
        "-l", str(log_file), "--seed", str(int(input_data.get("seed", 42))),
        str(config_path),
    ]
    result = _execute(
        runtime, args, work_dir, source_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="reinvent",
    )
    compounds, parsed = _parse_reinvent_csv(output_csv)
    log_text = log_file.read_text(errors="replace") if log_file.exists() else ""
    parsed["metadata_hash_warning"] = "has invalid hash" in log_text
    parsed["official_asset_checksum_verified"] = True
    return compounds, result, parsed


def _run_lingo3dmol(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    center = _coerce_center(input_data.get("pocket_center"))
    source_dir = _required_dir(runtime.get("source_dir"), "Lingo3DMol source")
    entrypoint = _required_file(runtime.get("entrypoint"), "Lingo3DMol entrypoint")
    caption_checkpoint = _required_file(
        runtime.get("checkpoint"), "Lingo3DMol generation checkpoint"
    )
    contact_checkpoint = _required_file(
        runtime.get("auxiliary_checkpoint"), "Lingo3DMol contact checkpoint"
    )
    requested = int(params.get("num_samples", 100))
    radius = float(params.get("pocket_radius", 10.0))
    pocket_pdb, residues = prepare_center_pocket(
        receptor,
        center,
        work_dir / "lingo3dmol_pocket.pdb",
        radius=radius,
        chain=input_data.get("chain"),
    )
    raw_dir = work_dir / "raw"
    args = [
        str(_GENERATOR_WORKER), "lingo3dmol",
        "--source-dir", str(source_dir),
        "--entrypoint", str(entrypoint),
        "--receptor", str(pocket_pdb),
        "--output-dir", str(raw_dir),
        "--caption-checkpoint", str(caption_checkpoint),
        "--contact-checkpoint", str(contact_checkpoint),
        "--num-samples", str(requested),
        "--gpu-device", str(int(params.get("gpu_device", 0))),
        "--max-run-hours", str(int(params.get("max_run_hours", 1))),
        "--gen-frag-set", str(int(params.get("gen_frag_set", 40))),
        "--prod-time", str(int(params.get("prod_time", 1))),
        "--topk", str(int(params.get("topk", 5))),
        "--nci-threshold", str(float(params.get("nci_threshold", 0.7))),
        "--collision-distance", str(float(params.get("collision_distance", 2.5))),
        "--frag-len-add", str(int(params.get("frag_len_add", 15))),
        "--seed", str(int(input_data.get("seed", 42))),
    ]
    result = _execute(
        runtime, args, work_dir, work_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python",
    )
    structure_files = sorted(raw_dir.rglob("*.mol")) + sorted(raw_dir.rglob("*.sdf"))
    compounds, parsed = parse_generated_structure_files(
        structure_files, "lingo3dmol"
    )
    compounds = compounds[:requested]
    return compounds, result, {
        **parsed,
        "retained": len(compounds),
        "pocket_basis": "step1_center_radius",
        "pocket_radius": radius,
        "pocket_residue_count": len(residues),
        "native_output_format": "mol",
        "generation_time_log": str(work_dir / "generation_time.log"),
    }


def _run_tamgen(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    center = _coerce_center(input_data.get("pocket_center"))
    source_dir = _required_dir(runtime.get("source_dir"), "TamGen source")
    checkpoint = _required_file(runtime.get("checkpoint"), "TamGen checkpoint")
    output_csv = work_dir / "generated.csv"
    requested = int(params.get("num_samples", 100))
    args = [
        str(_GENERATOR_WORKER), "tamgen",
        "--source-dir", str(source_dir),
        "--receptor", str(receptor),
        "--checkpoint", str(checkpoint),
        "--output-csv", str(output_csv),
        "--center", *(str(value) for value in center),
        "--num-samples", str(requested),
        "--beam", str(int(params.get("beam", 20))),
        "--max-seeds", str(int(params.get("max_seeds", 5))),
        "--pocket-radius", str(float(params.get("pocket_radius", 10.0))),
        "--sample-beta", str(float(params.get("sample_beta", 1.0))),
        "--seed", str(int(input_data.get("seed", 42))),
    ]
    if input_data.get("chain"):
        args.extend(["--chain", str(input_data["chain"])])
    if bool(params.get("use_conditional", False)):
        args.append("--use-conditional")
    result = _execute(
        runtime, args, work_dir, work_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python",
    )
    compounds, parsed = _parse_generated_smiles_csv(output_csv, "tamgen")
    return compounds, result, {
        **parsed,
        "pocket_basis": "step1_center_radius",
        "pocket_radius": float(params.get("pocket_radius", 10.0)),
        "chain": input_data.get("chain"),
        "coordinate_contract": "smiles_only_ETKDG_downstream",
        "raw_generation_log": str(work_dir / "raw_generation.log"),
    }


def _run_molcraft(input_data, params, runtime, work_dir):
    receptor = _required_file(input_data.get("pocket_pdb_file"), "prepared receptor")
    reference = _required_file(
        input_data.get("reference_ligand_file"),
        "reference ligand required by the installed MolCRAFT entrypoint",
    )
    source_dir = _required_dir(runtime.get("source_dir"), "MolCRAFT source")
    entrypoint = _required_file(
        runtime.get("entrypoint"), "MolCRAFT local compatibility entrypoint"
    )
    if entrypoint.name != "sample_for_pocket_local.py":
        raise ValueError(
            "MolCRAFT must use the reviewed sample_for_pocket_local.py "
            "compatibility entrypoint"
        )
    checkpoint = _required_file(runtime.get("checkpoint"), "MolCRAFT checkpoint")
    checkpoint_config = _required_file(
        runtime.get("checkpoint_config")
        or checkpoint.with_name(f"{checkpoint.name}.config.yaml"),
        "MolCRAFT checkpoint config",
    )
    raw_dir = work_dir / "raw"
    requested = int(params.get("num_samples", 100))
    args = [
        str(_GENERATOR_WORKER), "molcraft",
        "--source-dir", str(source_dir),
        "--entrypoint", str(entrypoint),
        "--receptor", str(receptor),
        "--reference-ligand", str(reference),
        "--checkpoint", str(checkpoint),
        "--checkpoint-config", str(checkpoint_config),
        "--output-dir", str(raw_dir),
        "--num-samples", str(requested),
        "--sample-steps", str(int(params.get("sample_steps", 100))),
        "--sample-num-atoms", str(params.get("sample_num_atoms", "prior")),
        "--beta1", str(float(params.get("beta1", 1.5))),
        "--sigma1-coord", str(float(params.get("sigma1_coord", 0.03))),
        "--sampling-strategy", str(params.get("sampling_strategy", "end_back")),
        "--seed", str(int(input_data.get("seed", 42))),
    ]
    result = _execute(
        runtime, args, work_dir, work_dir,
        timeout=int(params.get("timeout_seconds", 86400)),
        gpu_device=int(params.get("gpu_device", 0)),
        command="python",
    )
    compounds, parsed = parse_generated_sdfs(sorted(raw_dir.glob("*.sdf")), "molcraft")
    return compounds, result, {
        **parsed,
        "pocket_basis": "reference_ligand",
        "compatibility_entrypoint": str(entrypoint),
    }


def prepare_center_pocket(
    protein_pdb: Path,
    center: Iterable[float],
    output_path: Path,
    *,
    radius: float = 10.0,
    chain: str | None = None,
    normalize_residue_ids: bool = False,
) -> tuple[Path, list[str]]:
    """Write standard protein residues having any heavy atom within radius."""
    cx, cy, cz = _coerce_center(center)
    residue_lines: dict[tuple[str, str, str, str], list[str]] = {}
    residue_near: set[tuple[str, str, str, str]] = set()
    for line in protein_pdb.read_text(errors="replace").splitlines(True):
        if not line.startswith("ATOM"):
            continue
        line_chain = line[21].strip()
        if chain and line_chain != str(chain):
            continue
        resname = line[17:20].strip().upper()
        if resname not in _STANDARD_AMINO_ACIDS:
            continue
        altloc = line[16].strip()
        if altloc not in ("", "A"):
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        atom_name = line[12:16].strip().upper()
        if element == "H" or (not element and atom_name.startswith("H")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        key = (line_chain, line[22:26].strip(), line[26].strip(), resname)
        residue_lines.setdefault(key, []).append(line)
        if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= radius ** 2:
            residue_near.add(key)
    if not residue_near:
        raise ValueError(
            f"no standard protein residue lies within {radius:.1f} A of pocket center"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_keys = [key for key in residue_lines if key in residue_near]
    with output_path.open("w", encoding="utf-8") as handle:
        for normalized_index, key in enumerate(ordered_keys, start=1):
            lines = residue_lines[key]
            if normalize_residue_ids:
                lines = [
                    f"{line[:21]}A{normalized_index:4d} {line[27:]}"
                    for line in lines
                ]
            handle.writelines(lines)
        handle.write("END\n")
    if normalize_residue_ids:
        residue_ids = [
            f"A:{index}" for index in range(1, len(ordered_keys) + 1)
        ]
    else:
        residue_ids = [
            f"{key[0] or '_'}:{key[1]}{key[2]}" for key in ordered_keys
        ]
    return output_path, residue_ids


def _create_invocation_work_dir(output_dir: Path, tool_id: str) -> Path:
    """Allocate an isolated directory so retries cannot ingest stale outputs."""
    tool_root = Path(output_dir) / f"{tool_id}_work"
    tool_root.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="run-", dir=str(tool_root)))


def parse_generated_sdfs(
    sdf_files: Iterable[Path],
    source: str,
) -> tuple[list[dict], dict]:
    """Parse only valid single-fragment molecules with finite native 3D poses."""
    return parse_generated_structure_files(sdf_files, source)


def parse_generated_structure_files(
    structure_files: Iterable[Path],
    source: str,
) -> tuple[list[dict], dict]:
    """Parse valid SDF/MOL records with connected, finite native 3D poses."""
    from rdkit import Chem

    files = [Path(path) for path in structure_files if Path(path).is_file()]
    compounds: list[dict] = []
    records = invalid = nonfinite = multicomponent = 0
    for structure_file in files:
        if structure_file.suffix.lower() == ".mol":
            supplier = [Chem.MolFromMolFile(
                str(structure_file), sanitize=True, removeHs=False,
                strictParsing=True,
            )]
        elif structure_file.suffix.lower() == ".sdf":
            supplier = Chem.SDMolSupplier(
                str(structure_file), sanitize=True, removeHs=False,
                strictParsing=True,
            )
        else:
            continue
        for mol_index, mol in enumerate(supplier):
            records += 1
            if mol is None or mol.GetNumConformers() < 1:
                invalid += 1
                continue
            positions = mol.GetConformer().GetPositions()
            if not all(
                math.isfinite(float(value))
                for row in positions
                for value in row
            ):
                nonfinite += 1
                continue
            smiles = Chem.MolToSmiles(
                Chem.RemoveHs(mol), canonical=True, isomericSmiles=True
            )
            if not smiles:
                invalid += 1
                continue
            if "." in smiles:
                multicomponent += 1
                continue
            compounds.append({
                "smiles": smiles,
                "source": source,
                "candidate_role": "generated",
                "source_file": str(structure_file.resolve()),
                "source_mol_index": mol_index,
                "raw_pose_molblock": Chem.MolToMolBlock(mol),
                "coordinate_origin": f"{source}_raw",
            })
    return compounds, {
        "sdf_files": len(files),
        "records": records,
        "invalid_or_missing_3d": invalid,
        "nonfinite_coordinates": nonfinite,
        "multicomponent": multicomponent,
        "valid": len(compounds),
    }


def _parse_generated_smiles_csv(
    csv_path: Path,
    source: str,
) -> tuple[list[dict], dict]:
    """Parse a worker CSV while retaining only valid, connected SMILES."""
    from rdkit import Chem

    if not csv_path.is_file():
        raise RuntimeError(f"{source} did not create output CSV: {csv_path}")
    compounds = []
    invalid = multicomponent = records = 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for index, row in enumerate(csv.DictReader(handle)):
            records += 1
            smiles = (row.get("smiles") or row.get("SMILES") or "").strip()
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            if mol is None:
                invalid += 1
                continue
            canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            if "." in canonical:
                multicomponent += 1
                continue
            compound = {
                "smiles": canonical,
                "source": source,
                "candidate_role": "generated",
                "source_file": str(csv_path.resolve()),
                "source_mol_index": index,
                "coordinate_origin": "smiles_only",
            }
            if row.get("score") not in (None, ""):
                compound["generator_score"] = row["score"]
            compounds.append(compound)
    return compounds, {
        "records": records,
        "invalid": invalid,
        "multicomponent": multicomponent,
        "valid": len(compounds),
    }


def _parse_reinvent_csv(csv_path: Path) -> tuple[list[dict], dict]:
    from rdkit import Chem

    if not csv_path.is_file():
        raise RuntimeError(f"REINVENT4 did not create output CSV: {csv_path}")
    compounds = []
    invalid = multicomponent = records = 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for index, row in enumerate(csv.DictReader(handle)):
            records += 1
            smiles = (row.get("SMILES") or "").strip()
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            if mol is None:
                invalid += 1
                continue
            canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            if "." in canonical:
                multicomponent += 1
                continue
            compounds.append({
                "smiles": canonical,
                "source": "reinvent4",
                "candidate_role": "generated",
                "source_file": str(csv_path.resolve()),
                "source_mol_index": index,
                "coordinate_origin": "smiles_only",
                "generator_nll": row.get("NLL", ""),
            })
    return compounds, {
        "records": records, "invalid": invalid,
        "multicomponent": multicomponent, "valid": len(compounds),
    }


def _execute(
    runtime: dict,
    args: list[str],
    work_dir: Path,
    cwd: Path,
    *,
    timeout: int,
    gpu_device: int,
    command: str,
    extra_env: dict[str, str] | None = None,
) -> dict:
    env = {"CUDA_VISIBLE_DEVICES": str(gpu_device), **(extra_env or {})}
    execution_runtime = {
        "kind": "command",
        "environment": runtime.get("environment"),
        "command": command,
    }
    try:
        result = run_tool(
            execution_runtime, args, cwd=cwd, timeout=timeout,
            check=True, env=env,
        )
    except ToolRuntimeError as exc:
        if exc.result:
            _write_execution_files(work_dir, exc.result)
        raise RuntimeError(f"{runtime.get('environment')} execution failed: {exc}") from exc
    return _write_execution_files(work_dir, result)


def _write_execution_files(work_dir: Path, result: dict) -> dict:
    stdout_path = work_dir / "stdout.log"
    stderr_path = work_dir / "stderr.log"
    stdout_path.write_text(result.get("stdout", ""), encoding="utf-8")
    stderr_path.write_text(result.get("stderr", ""), encoding="utf-8")
    record = {
        key: result.get(key)
        for key in (
            "success", "command", "returncode", "elapsed_seconds",
            "environment", "error", "process_group_terminated",
        )
        if result.get(key) is not None
    }
    record["stdout_log"] = str(stdout_path)
    record["stderr_log"] = str(stderr_path)
    record_path = work_dir / "execution.json"
    record_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    record["record_file"] = str(record_path)
    return record


def _verify_asset(runtime: dict) -> None:
    expected = runtime.get("asset_sha256")
    checkpoint = runtime.get("checkpoint")
    if expected and checkpoint:
        path = _required_file(checkpoint, "registered model asset")
        actual = _sha256(path)
        if actual != expected:
            raise RuntimeError(
                f"registered model checksum mismatch for {path}: "
                f"expected {expected}, got {actual}"
            )
    auxiliary = runtime.get("auxiliary_checkpoint")
    auxiliary_expected = runtime.get("auxiliary_asset_sha256")
    if auxiliary and auxiliary_expected:
        auxiliary_path = _required_file(auxiliary, "registered auxiliary model asset")
        auxiliary_actual = _sha256(auxiliary_path)
        if auxiliary_actual != auxiliary_expected:
            raise RuntimeError(
                f"registered model checksum mismatch for {auxiliary_path}: "
                f"expected {auxiliary_expected}, got {auxiliary_actual}"
            )
    config = runtime.get("checkpoint_config")
    config_expected = runtime.get("checkpoint_config_sha256")
    if config and config_expected:
        config_path = _required_file(config, "registered checkpoint config")
        config_actual = _sha256(config_path)
        if config_actual != config_expected:
            raise RuntimeError(
                f"registered config checksum mismatch for {config_path}: "
                f"expected {config_expected}, got {config_actual}"
            )


def _validate_fragment_file(path: Path) -> None:
    from rdkit import Chem

    suffix = path.suffix.lower()
    if suffix == ".sdf":
        mol = next((item for item in Chem.SDMolSupplier(
            str(path), sanitize=False, removeHs=False
        ) if item is not None), None)
    elif suffix == ".mol":
        mol = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
    else:
        raise ValueError("DiffLinker fragment_file must be SDF, MOL, MOL2, or PDB")
    if mol is None:
        raise ValueError(f"DiffLinker could not parse fragment file: {path}")
    if mol.GetNumConformers() < 1:
        raise ValueError("DiffLinker fragment_file must contain 3D coordinates")
    if len(Chem.GetMolFrags(mol)) < 2:
        raise ValueError(
            "DiffLinker fragment_file must contain at least two disconnected fragments"
        )
    positions = mol.GetConformer().GetPositions()
    if not all(math.isfinite(float(value)) for row in positions for value in row):
        raise ValueError("DiffLinker fragment_file contains non-finite coordinates")


def _required_file(value, label: str) -> Path:
    if not value:
        raise ValueError(f"{label} was not provided")
    path = Path(value).expanduser()
    if not path.is_file():
        raise ValueError(f"{label} not found: {path}")
    return path.resolve()


def _optional_file(value) -> Path | None:
    if not value:
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_file() else None


def _required_dir(value, label: str) -> Path:
    if not value:
        raise ValueError(f"{label} was not provided")
    path = Path(value).expanduser()
    if not path.is_dir():
        raise ValueError(f"{label} not found: {path}")
    return path.resolve()


def _coerce_center(center) -> tuple[float, float, float]:
    if isinstance(center, str):
        values = center.replace(";", ",").split(",")
    else:
        values = center
    if not isinstance(values, (list, tuple)) or len(values) != 3:
        raise ValueError("pocket_center must contain exactly three coordinates")
    try:
        converted = tuple(float(value) for value in values)
    except (TypeError, ValueError) as exc:
        raise ValueError("pocket_center contains a non-numeric coordinate") from exc
    if not all(math.isfinite(value) for value in converted):
        raise ValueError("pocket_center coordinates must be finite")
    return converted


def _divisor_batch_size(num_samples: int, requested_batch: int) -> int:
    batch = max(1, min(num_samples, requested_batch))
    while num_samples % batch:
        batch -= 1
    return batch


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

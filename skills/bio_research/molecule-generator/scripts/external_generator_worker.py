#!/usr/bin/env python3
"""Environment-local workers for optional Step 2 generators.

The parent adapter launches this file through the registry-owned Conda runtime.
Keeping upstream imports and process-global state here prevents incompatible
PyTorch/RDKit stacks from leaking into the Drugmate process.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
import random
import re
import runpy
import shutil
import subprocess
import sys
from pathlib import Path


def _lingo3dmol(args: argparse.Namespace) -> None:
    source_dir = Path(args.source_dir).resolve()
    work_dir = Path.cwd().resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # The upstream script incorporates the input basename into an ``rm -rf``
    # command.  Give it a fixed safe basename wholly inside this invocation.
    receptor = work_dir / "pocket.pdb"
    shutil.copyfile(Path(args.receptor).resolve(), receptor)
    input_list = work_dir / "input_list.csv"
    input_list.write_text(f",,{receptor}\n", encoding="utf-8")

    random.seed(args.seed)
    try:
        import numpy as np
        np.random.seed(args.seed)
    except ImportError:
        pass
    try:
        import torch
        torch.manual_seed(args.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(args.seed)
    except ImportError:
        pass

    savedir = str(output_dir) + os.sep
    argv = [
        str(args.entrypoint),
        "--cuda", str(args.gpu_device),
        "--cuda_list", str(args.gpu_device),
        "--input_list", str(input_list),
        "--savedir", savedir,
        "--contact_path", str(Path(args.contact_checkpoint).resolve()),
        "--caption_path", str(Path(args.caption_checkpoint).resolve()),
        "--max_run_hours", str(args.max_run_hours),
        "--gennums", str(args.num_samples),
        "--gen_frag_set", str(args.gen_frag_set),
        "--prod_time", str(args.prod_time),
        "--topk", str(args.topk),
        "--nci_thrs", str(args.nci_threshold),
        "--coc_dis", str(args.collision_distance),
        "--frag_len_add", str(args.frag_len_add),
        "--save_time_path", str(work_dir / "generation_time.log"),
    ]
    previous_cwd = Path.cwd()
    previous_argv = sys.argv
    try:
        os.chdir(source_dir)
        sys.path.insert(0, str(source_dir))
        sys.argv = argv
        runpy.run_path(str(Path(args.entrypoint).resolve()), run_name="__main__")
    finally:
        sys.argv = previous_argv
        os.chdir(previous_cwd)


def _tamgen(args: argparse.Namespace) -> None:
    source_dir = Path(args.source_dir).resolve()
    receptor = Path(args.receptor).resolve()
    data_dir = Path.cwd().resolve() / "tamgen_data"
    raw_log = Path.cwd().resolve() / "raw_generation.log"
    center = tuple(float(value) for value in args.center)
    _build_tamgen_dataset(
        source_dir=source_dir,
        receptor=receptor,
        data_dir=data_dir,
        center=center,
        radius=float(args.pocket_radius),
        chain=args.chain,
    )

    beam = max(1, min(int(args.beam), int(args.num_samples)))
    unique: dict[str, tuple[float, int]] = {}
    log_chunks: list[str] = []
    for offset in range(max(1, int(args.max_seeds))):
        if len(unique) >= int(args.num_samples):
            break
        seed = int(args.seed) + offset
        command = [
            sys.executable, str(source_dir / "generate.py"), str(data_dir),
            "-s", "tg", "-t", "m1",
            "--task", "translation_coord",
            "--path", str(Path(args.checkpoint).resolve()),
            "--gen-subset", "test",
            "--beam", str(beam), "--nbest", str(beam),
            "--max-tokens", "1024",
            "--seed", str(seed),
            "--sample-beta", str(args.sample_beta),
            "--use-src-coord",
        ]
        if args.use_conditional:
            command.append("--gen-vae")
        completed = subprocess.run(
            command,
            cwd=source_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        log_chunks.append(f"# seed={seed}\n{completed.stdout}")
        raw_log.write_text("\n".join(log_chunks), encoding="utf-8")
        if completed.returncode != 0:
            sys.stderr.write(completed.stdout)
            raise RuntimeError(
                f"TamGen generation failed for seed {seed} "
                f"with exit code {completed.returncode}"
            )
        for smiles, score in _parse_tamgen_stdout(completed.stdout):
            unique.setdefault(smiles, (score, seed))

    output_csv = Path(args.output_csv).resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["smiles", "score", "seed"])
        writer.writeheader()
        for smiles, (score, seed) in list(unique.items())[: int(args.num_samples)]:
            writer.writerow({"smiles": smiles, "score": score, "seed": seed})


def _build_tamgen_dataset(
    *,
    source_dir: Path,
    receptor: Path,
    data_dir: Path,
    center: tuple[float, float, float],
    radius: float,
    chain: str | None,
) -> None:
    """Build TamGen's binary pocket input directly from the prepared PDB."""
    import numpy as np
    from Bio.PDB import PDBParser

    # Upstream ``dump_center_data`` invokes a literal ``python`` subprocess.
    # Anchor that lookup to this worker's environment even when the prefix was
    # launched directly instead of through an activated shell.
    environment_bin = str(Path(sys.executable).resolve().parent)
    os.environ["PATH"] = environment_bin + os.pathsep + os.environ.get("PATH", "")

    sys.path.insert(0, str(source_dir))
    from fairseq.molecule_utils.coordinate_utils.atom_positions import (
        get_residue_atom_coordinates,
        get_residue_average_position,
    )
    from fairseq.molecule_utils.database.common_utils import aa_3to1
    from fairseq.molecule_utils.external.fairseq_dataset_build_utils import (
        dump_center_data,
    )

    structure = PDBParser(QUIET=True).get_structure("pipeline", str(receptor))
    model = next(structure.get_models(), None)
    if model is None:
        raise RuntimeError(f"TamGen could not parse receptor PDB: {receptor}")

    sequences: dict[str, str] = {}
    coordinates: dict[str, object] = {}
    near_center_masks: dict[str, object] = {}
    res_ids: dict[str, object] = {}
    center_array = np.asarray(center, dtype=np.float32)
    for pdb_chain in model:
        chain_id = str(pdb_chain.id).strip() or "_"
        if chain and chain_id != str(chain):
            continue
        residues = []
        for residue in pdb_chain:
            if residue.id[0] not in (" ", "W") and residue.get_resname() != "MSE":
                continue
            code = aa_3to1(residue.get_resname().upper())
            if code == "X":
                continue
            atom_coordinates = get_residue_atom_coordinates(residue)
            if atom_coordinates.size == 0:
                continue
            residues.append((residue, code, atom_coordinates))
        if not residues:
            continue
        sequences[chain_id] = "".join(item[1] for item in residues)
        coordinates[chain_id] = np.asarray([
            get_residue_average_position(
                item[0], dtype=np.float32, center_of_gravity=True
            )
            for item in residues
        ], dtype=np.float32)
        mask = np.asarray([
            bool(np.min(np.linalg.norm(item[2] - center_array, axis=1)) <= radius)
            for item in residues
        ], dtype=np.int32)
        near_center_masks[chain_id] = mask
        res_ids[chain_id] = np.asarray([
            int(item[0].id[1]) for item, selected in zip(residues, mask)
            if selected
        ], dtype=np.int32)

    if not sequences:
        suffix = f" for chain {chain!r}" if chain else ""
        raise RuntimeError(f"TamGen found no standard protein residues{suffix}")
    if not any(int(mask.sum()) for mask in near_center_masks.values()):
        raise RuntimeError(
            f"TamGen found no protein residue within {radius:.1f} A of pocket center"
        )

    data = {
        "index": 0,
        "pdb_id": receptor.stem.lower(),
        "sequences": sequences,
        "uniprot_refs": {},
        "coordinates": coordinates,
        "zero_masks": {
            key: np.any(value != 0, axis=1).astype(np.int32)
            for key, value in coordinates.items()
        },
        "near_center_masks": near_center_masks,
        "res_ids": res_ids,
        "center": center_array,
        "all_chain_ids": sorted(sequences),
    }
    dump_center_data(
        [data], "test", output_dir=data_dir,
        fairseq_root=source_dir, pre_dicts_root=source_dir / "dict",
        max_len=1023,
    )


def _parse_tamgen_stdout(text: str) -> list[tuple[str, float]]:
    from rdkit import Chem

    parsed: list[tuple[str, float]] = []
    for line in text.splitlines():
        if not line.startswith("H-"):
            continue
        fields = line.split("\t", 2)
        if len(fields) != 3:
            continue
        raw_smiles = fields[2].replace(" ", "").replace("[generation]", "")
        mol = Chem.MolFromSmiles(raw_smiles)
        if mol is None:
            continue
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if not smiles or "." in smiles:
            continue
        try:
            score = float(fields[1])
        except ValueError:
            score = math.nan
        parsed.append((smiles, score))
    return parsed


def _molcraft(args: argparse.Namespace) -> None:
    source_dir = Path(args.source_dir).resolve()
    entrypoint = Path(args.entrypoint).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # The reviewed local entrypoint resolves this exact relative config path.
    # Copy it into the invocation and set the requested deterministic seed.
    checkpoint_dir = Path.cwd().resolve() / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    config_text = Path(args.checkpoint_config).read_text(encoding="utf-8")
    config_text, replacements = re.subn(
        r"(?m)^seed:\s*.*$", f"seed: {int(args.seed)}", config_text, count=1
    )
    if replacements != 1:
        raise RuntimeError("MolCRAFT checkpoint config has no top-level seed field")
    (checkpoint_dir / "last.ckpt.config.yaml").write_text(
        config_text, encoding="utf-8"
    )

    sys.path.insert(0, str(source_dir))
    specification = importlib.util.spec_from_file_location(
        "pipeline_molcraft_entrypoint", entrypoint
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load MolCRAFT entrypoint: {entrypoint}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)

    # Redirect the upstream process-global output directory into this fresh
    # invocation instead of deleting/reusing the source checkout's output/.
    import core.callbacks.validation_callback_for_sample as callback_module
    callback_module.OUT_DIR = str(output_dir)
    module.call(
        str(Path(args.receptor).resolve()),
        str(Path(args.reference_ligand).resolve()),
        ckpt_path=str(Path(args.checkpoint).resolve()),
        num_samples=int(args.num_samples),
        sample_steps=int(args.sample_steps),
        sample_num_atoms=args.sample_num_atoms,
        beta1=float(args.beta1),
        sigma1_coord=float(args.sigma1_coord),
        sampling_strategy=args.sampling_strategy,
        seed=int(args.seed),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="worker", required=True)

    lingo = subparsers.add_parser("lingo3dmol")
    lingo.add_argument("--source-dir", required=True)
    lingo.add_argument("--entrypoint", required=True)
    lingo.add_argument("--receptor", required=True)
    lingo.add_argument("--output-dir", required=True)
    lingo.add_argument("--caption-checkpoint", required=True)
    lingo.add_argument("--contact-checkpoint", required=True)
    lingo.add_argument("--num-samples", type=int, required=True)
    lingo.add_argument("--gpu-device", type=int, default=0)
    lingo.add_argument("--max-run-hours", type=int, default=1)
    lingo.add_argument("--gen-frag-set", type=int, default=40)
    lingo.add_argument("--prod-time", type=int, default=1)
    lingo.add_argument("--topk", type=int, default=5)
    lingo.add_argument("--nci-threshold", type=float, default=0.7)
    lingo.add_argument("--collision-distance", type=float, default=2.5)
    lingo.add_argument("--frag-len-add", type=int, default=15)
    lingo.add_argument("--seed", type=int, default=42)
    lingo.set_defaults(function=_lingo3dmol)

    tamgen = subparsers.add_parser("tamgen")
    tamgen.add_argument("--source-dir", required=True)
    tamgen.add_argument("--receptor", required=True)
    tamgen.add_argument("--checkpoint", required=True)
    tamgen.add_argument("--output-csv", required=True)
    tamgen.add_argument("--center", type=float, nargs=3, required=True)
    tamgen.add_argument("--num-samples", type=int, required=True)
    tamgen.add_argument("--beam", type=int, default=20)
    tamgen.add_argument("--max-seeds", type=int, default=5)
    tamgen.add_argument("--pocket-radius", type=float, default=10.0)
    tamgen.add_argument("--sample-beta", type=float, default=1.0)
    tamgen.add_argument("--seed", type=int, default=42)
    tamgen.add_argument("--chain")
    tamgen.add_argument("--use-conditional", action="store_true")
    tamgen.set_defaults(function=_tamgen)

    molcraft = subparsers.add_parser("molcraft")
    molcraft.add_argument("--source-dir", required=True)
    molcraft.add_argument("--entrypoint", required=True)
    molcraft.add_argument("--receptor", required=True)
    molcraft.add_argument("--reference-ligand", required=True)
    molcraft.add_argument("--checkpoint", required=True)
    molcraft.add_argument("--checkpoint-config", required=True)
    molcraft.add_argument("--output-dir", required=True)
    molcraft.add_argument("--num-samples", type=int, required=True)
    molcraft.add_argument("--sample-steps", type=int, default=100)
    molcraft.add_argument("--sample-num-atoms", default="prior")
    molcraft.add_argument("--beta1", type=float, default=1.5)
    molcraft.add_argument("--sigma1-coord", type=float, default=0.03)
    molcraft.add_argument("--sampling-strategy", default="end_back")
    molcraft.add_argument("--seed", type=int, default=42)
    molcraft.set_defaults(function=_molcraft)
    return parser


def main() -> None:
    args = _parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()

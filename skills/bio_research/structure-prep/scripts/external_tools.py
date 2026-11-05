"""Thin, auditable wrappers for external Step 0 tools."""
from __future__ import annotations

import json
import time
import urllib.request
from pathlib import Path
from typing import Any

from tool_runtime import ToolRuntimeError, environment_spec, run_tool


def run_pdb2pqr(
    input_pdb: Path,
    output_pdb: Path,
    output_pqr: Path,
    runtime: dict,
    *,
    ph: float,
    force_field: str = "AMBER",
    optimize_hbond_network: bool = True,
) -> dict[str, Any]:
    args = [
        "--ff", force_field, "--keep-chain",
        "--titration-state-method", "propka", "--with-ph", ph,
        "--pdb-output", output_pdb,
    ]
    if not optimize_hbond_network:
        args.append("--noopt")
    args.extend([input_pdb, output_pqr])
    result = run_tool(runtime, args, cwd=output_pdb.parent, timeout=600)
    if not output_pdb.exists() or output_pdb.stat().st_size == 0:
        raise ToolRuntimeError("pdb2pqr completed without creating the requested PDB", result)
    return {**result, "output_pdb": str(output_pdb), "output_pqr": str(output_pqr)}


def run_reduce2(
    input_pdb: Path,
    output_pdb: Path,
    runtime: dict,
    *,
    build: bool = True,
    flip: bool = True,
) -> dict[str, Any]:
    args = [
        input_pdb,
        f"approach={'add' if build else 'optimize'}",
        f"add_flip_movers={str(flip)}",
        f"output.filename={output_pdb}",
        "output.overwrite=True",
        "--quiet",
    ]
    result = run_tool(runtime, args, cwd=output_pdb.parent, timeout=600)
    if not output_pdb.exists() or output_pdb.stat().st_size == 0:
        raise ToolRuntimeError("Reduce2 completed without creating a non-empty PDB", result)
    return {**result, "output_pdb": str(output_pdb)}


def run_protoss(
    pdb_id: str,
    output_pdb: Path,
    runtime: dict,
    *,
    poll_interval_seconds: int = 5,
    timeout_seconds: int = 600,
) -> dict[str, Any]:
    """Submit and retrieve a ProteinsPlus Protoss REST job by PDB code."""
    env = environment_spec(runtime)
    endpoint = env["base_url"].rstrip("/") + runtime.get("endpoint", "/protoss_rest")
    payload = json.dumps({"protoss": {"pdbCode": pdb_id.upper()}}).encode()
    request = urllib.request.Request(
        endpoint,
        data=payload,
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        submitted = json.loads(response.read().decode())
    location = submitted.get("location")
    if not location:
        raise RuntimeError(f"Protoss response did not include a job location: {submitted}")
    if not str(location).startswith("http"):
        location = env["base_url"].rstrip("/") + "/protoss_rest/" + str(location).strip("/").split("/")[-1]

    started = time.monotonic()
    job = submitted
    while time.monotonic() - started < timeout_seconds:
        with urllib.request.urlopen(str(location), timeout=30) as response:
            job = json.loads(response.read().decode())
        if job.get("status_code") == 200 and job.get("protein"):
            break
        time.sleep(poll_interval_seconds)
    else:
        raise RuntimeError(f"Protoss job timed out after {timeout_seconds}s: {location}")

    protein_url = job["protein"]
    if not str(protein_url).startswith("http"):
        protein_url = "https://proteins.plus" + "/" + str(protein_url).lstrip("/")
    urllib.request.urlretrieve(str(protein_url), str(output_pdb))
    return {
        "success": True,
        "pdb_id": pdb_id.upper(),
        "job_location": str(location),
        "protein_url": str(protein_url),
        "output_pdb": str(output_pdb),
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }


def run_dimorphite(
    smiles: str,
    output_smi: Path,
    runtime: dict,
    *,
    ph: float,
    ph_tolerance: float = 1.0,
    max_variants: int = 16,
) -> dict[str, Any]:
    args = [
        smiles,
        "--ph_min", max(0.0, ph - ph_tolerance),
        "--ph_max", min(14.0, ph + ph_tolerance),
        "--max_variants", max_variants,
        "--output_file", output_smi,
    ]
    result = run_tool(runtime, args, cwd=output_smi.parent, timeout=300)
    if not output_smi.exists():
        raise ToolRuntimeError("Dimorphite-DL did not produce its SMILES output", result)
    variants = [line.strip() for line in output_smi.read_text().splitlines() if line.strip()]
    return {**result, "output_smi": str(output_smi), "variant_count": len(variants)}


def run_openbabel(
    input_file: Path,
    output_file: Path,
    runtime: dict,
) -> dict[str, Any]:
    result = run_tool(runtime, [input_file, "-O", output_file], cwd=output_file.parent, timeout=300)
    if not output_file.exists():
        raise ToolRuntimeError("Open Babel did not produce the requested output", result)
    return {**result, "output_file": str(output_file)}


def run_meeko_receptor(
    input_pdb: Path,
    output_pdbqt: Path,
    runtime: dict,
    *,
    allow_bad_residues: bool = False,
) -> dict[str, Any]:
    args = ["--read_pdb", input_pdb, "--write_pdbqt", output_pdbqt]
    if allow_bad_residues:
        args.append("--allow_bad_res")
    result = run_tool(runtime, args, cwd=output_pdbqt.parent, timeout=600)
    if not output_pdbqt.exists():
        raise ToolRuntimeError("Meeko did not produce receptor PDBQT", result)
    return {**result, "output_pdbqt": str(output_pdbqt)}


def run_amber_topology(
    receptor_pdb: Path,
    ligand_sdf: Path,
    output_dir: Path,
    runtime: dict,
    *,
    protein_force_field: str = "ff14SB",
    ligand_charge_method: str = "bcc",
) -> dict[str, Any]:
    commands = runtime.get("commands") or {}
    cleaned = output_dir / "amber_receptor.pdb"
    ligand_mol2 = output_dir / "amber_ligand.mol2"
    ligand_frcmod = output_dir / "amber_ligand.frcmod"
    prmtop = output_dir / "complex.prmtop"
    inpcrd = output_dir / "complex.inpcrd"
    leap_input = output_dir / "tleap.in"

    records = []
    records.append(run_tool(
        runtime, ["-i", receptor_pdb, "-o", cleaned, "--dry", "--nohyd"],
        command=commands.get("pdb4amber", "pdb4amber"), cwd=output_dir, timeout=600,
    ))
    records.append(run_tool(
        runtime,
        ["-i", ligand_sdf, "-fi", "sdf", "-o", ligand_mol2, "-fo", "mol2",
         "-at", "gaff2", "-c", ligand_charge_method, "-rn", "LIG", "-s", "2"],
        command=commands.get("antechamber", "antechamber"), cwd=output_dir, timeout=1800,
    ))
    records.append(run_tool(
        runtime, ["-i", ligand_mol2, "-f", "mol2", "-o", ligand_frcmod, "-s", "gaff2"],
        command=commands.get("parmchk2", "parmchk2"), cwd=output_dir, timeout=600,
    ))
    leap_input.write_text(
        f"source leaprc.protein.{protein_force_field}\n"
        "source leaprc.gaff2\n"
        f"REC = loadpdb {cleaned}\n"
        f"LIG = loadmol2 {ligand_mol2}\n"
        f"loadamberparams {ligand_frcmod}\n"
        "COM = combine { REC LIG }\n"
        f"saveamberparm COM {prmtop} {inpcrd}\n"
        "quit\n"
    )
    records.append(run_tool(
        runtime, ["-f", leap_input], command=commands.get("tleap", "tleap"),
        cwd=output_dir, timeout=600,
    ))
    missing = [path for path in (prmtop, inpcrd, ligand_mol2, ligand_frcmod) if not path.exists()]
    if missing:
        raise ToolRuntimeError("Amber composite did not create: " + ", ".join(map(str, missing)))
    return {
        "success": True,
        "commands": records,
        "prmtop": str(prmtop),
        "inpcrd": str(inpcrd),
        "ligand_mol2": str(ligand_mol2),
        "ligand_frcmod": str(ligand_frcmod),
        "tleap_input": str(leap_input),
    }

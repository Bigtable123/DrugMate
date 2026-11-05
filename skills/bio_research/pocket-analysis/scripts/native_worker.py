#!/usr/bin/env python3
"""Run crash-prone native interaction profilers in an isolated process."""
from __future__ import annotations

import argparse
import faulthandler
import json
import sys
from pathlib import Path


faulthandler.enable(all_threads=True)


def _run_step(step: str, input_data: dict) -> dict:
    if step == "plip":
        from plip_utils import run_plip_interactions

        return run_plip_interactions(
            Path(input_data["protein_file"]),
            Path(input_data["ligand_file"]),
            input_data["ligand_resname"],
            input_data.get("chain", "A"),
        )

    if step == "prolif":
        from prolif_utils import run_prolif_fingerprint

        return run_prolif_fingerprint(
            Path(input_data["protein_file"]),
            Path(input_data["ligand_file"]),
            input_data["ligand_resname"],
            interaction_types=input_data.get("interaction_types") or None,
            vicinity_cutoff=float(input_data.get("vicinity_cutoff", 6.5)),
            preserve_hetero_residues=input_data.get("preserve_hetero_residues") or [],
        )

    if step == "arpeggio":
        from arpeggio_utils import run_arpeggio_interactions

        return run_arpeggio_interactions(
            Path(input_data["complex_mmcif_file"]),
            Path(input_data["ligand_file"]),
            input_data["ligand_resname"],
            input_data.get("chain", "A"),
            timeout_seconds=int(input_data.get("tool_timeout_seconds", 300)),
        )

    if step == "luna":
        from luna_utils import run_luna_fingerprint

        return run_luna_fingerprint(
            Path(input_data["protein_file"]),
            Path(input_data["ligand_file"]),
            input_data["ligand_resname"],
            ifp_length=int(input_data.get("ifp_length", 4096)),
            ifp_type=str(input_data.get("ifp_type", "EIFP")),
            add_dependent_interactions=bool(
                input_data.get("add_dependent_interactions", False)
            ),
        )

    if step == "oddt":
        import os

        # ODDT chooses a backend at import time.  Pin the tested OpenBabel
        # implementation before the wrapper import.
        os.environ["ODDT_TOOLKIT"] = "ob"
        from oddt_utils import run_oddt_fingerprints

        return run_oddt_fingerprints(
            Path(input_data["protein_file"]),
            Path(input_data["ligand_file"]),
            input_data["ligand_resname"],
            strict=bool(input_data.get("strict", True)),
            splif_depth=int(input_data.get("splif_depth", 1)),
            splif_size=int(input_data.get("splif_size", 4096)),
            splif_distance_cutoff=float(
                input_data.get("splif_distance_cutoff", 4.5)
            ),
        )

    raise ValueError(f"Unsupported native step: {step}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--step",
        choices=("plip", "prolif", "arpeggio", "luna", "oddt"),
        required=True,
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output_path = Path(args.output)
    try:
        result = _run_step(args.step, json.loads(args.input))
        payload = {"ok": True, "result": result}
        exit_code = 0
    except Exception as exc:
        payload = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        exit_code = 1

    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

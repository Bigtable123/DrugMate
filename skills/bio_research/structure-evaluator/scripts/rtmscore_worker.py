#!/usr/bin/env python3
"""Isolated RTMScore inference worker.

This file is executed inside the pinned ``rtmscore_clean`` Conda prefix.  The
parent adapter owns pose identity and output validation; this worker only loads
one pinned model and returns the upstream scores in input order.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path


def _load_upstream(entrypoint: Path):
    source_dir = entrypoint.parent.parent
    sys.path.insert(0, str(source_dir))
    spec = importlib.util.spec_from_file_location("pipeline_rtmscore_upstream", entrypoint)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load RTMScore entrypoint: {entrypoint}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(payload: dict) -> dict:
    entrypoint = Path(payload["entrypoint"]).resolve()
    receptor = Path(payload["pocket_pdb"]).resolve()
    ligands = Path(payload["ligand_sdf"]).resolve()
    checkpoint = Path(payload["checkpoint"]).resolve()
    for path in (entrypoint, receptor, ligands, checkpoint):
        if not path.is_file():
            raise FileNotFoundError(path)

    upstream = _load_upstream(entrypoint)
    import MDAnalysis
    import dgl
    import numpy
    import prody
    import rdkit
    import torch

    requested_device = str(payload.get("device", "auto"))
    if requested_device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    elif requested_device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("RTMScore CUDA was requested but torch.cuda is unavailable")
        device = "cuda"
    elif requested_device == "cpu":
        device = "cpu"
    else:
        raise ValueError("device must be auto, cuda, or cpu")

    kwargs = {
        "batch_size": int(payload.get("batch_size", 128)),
        "dist_threhold": 5,
        "device": device,
        "num_workers": int(payload.get("num_workers", 0)),
        "num_node_featsp": 41,
        "num_node_featsl": 41,
        "num_edge_featsp": 5,
        "num_edge_featsl": 10,
        "hidden_dim0": 128,
        "hidden_dim": 128,
        "n_gaussians": 10,
        "dropout_rate": 0.10,
    }
    identifiers, predictions = upstream.scoring(
        prot=str(receptor),
        lig=str(ligands),
        modpath=str(checkpoint),
        cut=float(payload.get("pocket_cutoff_angstrom", 10.0)),
        gen_pocket=False,
        reflig=None,
        explicit_H=False,
        use_chirality=True,
        parallel=bool(payload.get("parallel_graph_build", False)),
        **kwargs,
    )
    ids = [str(value) for value in identifiers]
    values = [float(value) for value in predictions]
    if len(ids) != len(values):
        raise RuntimeError("RTMScore returned mismatched identifiers and scores")
    if len(set(ids)) != len(ids):
        raise RuntimeError("RTMScore returned duplicate identifiers")
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError("RTMScore returned a non-finite score")
    return {
        "success": True,
        "device": device,
        "torch_version": str(torch.__version__),
        "runtime_versions": {
            "torch": str(torch.__version__),
            "dgl": str(dgl.__version__),
            "rdkit": str(rdkit.__version__),
            "numpy": str(numpy.__version__),
            "mdanalysis": str(MDAnalysis.__version__),
            "prody": str(prody.__version__),
        },
        "cuda_available": bool(torch.cuda.is_available()),
        "ids": ids,
        "scores": values,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    payload = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
    result = run(payload)
    output = Path(args.output_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

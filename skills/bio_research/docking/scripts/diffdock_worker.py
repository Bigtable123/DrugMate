#!/usr/bin/env python3
"""Deterministic launcher for the locally installed DiffDock inference CLI."""
from __future__ import annotations

import argparse
import os
import random
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("diffdock_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    forwarded = list(args.diffdock_args)
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    source_root = Path(args.source_root).resolve()
    if not (source_root / "inference.py").is_file():
        parser.error(f"DiffDock inference.py not found under {source_root}")

    # DiffDock does not expose a seed in its public parser.  Seed every RNG
    # before importing/invoking the model while retaining its native CLI.
    os.chdir(source_root)
    sys.path.insert(0, str(source_root))
    random.seed(args.seed)
    import numpy as np
    import torch

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)

    import inference

    sys.argv = ["inference", *forwarded]
    inference_args = inference.get_parser().parse_args()
    inference.main(inference_args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run a legacy Python entrypoint after seeding its available RNG backends."""
from __future__ import annotations

import argparse
import random
import runpy
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("script_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()

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

    forwarded = list(args.script_args)
    if forwarded[:1] == ["--"]:
        forwarded = forwarded[1:]
    sys.path.insert(0, str(args.script.resolve().parent))
    sys.argv = [str(args.script), *forwarded]
    runpy.run_path(str(args.script), run_name="__main__")


if __name__ == "__main__":
    main()

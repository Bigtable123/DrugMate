"""Lightweight functional resource probe for the isolated AiZynthFinder env."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any


def _resource_path(value: Any, config_path: Path) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    path = path.resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"missing or empty AiZynthFinder resource: {path}")
    return path


def probe(config_file: str | Path) -> dict:
    """Validate imports, config keys, ONNX graphs, templates, and HDF metadata.

    Loading the complete 1.3 GB stock takes roughly 100 seconds on this host.
    The actual worker necessarily pays that cost once; preflight validates the
    HDF table metadata instead of loading it a second time.
    """
    import importlib.metadata

    import onnxruntime as ort
    import pandas as pd
    import yaml
    from aizynthfinder.aizynthfinder import AiZynthFinder  # noqa: F401

    config_path = Path(config_file).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"AiZynthFinder config not found: {config_path}")
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}

    expansion = config.get("expansion") or {}
    filters = config.get("filter") or {}
    stocks = config.get("stock") or {}
    if "uspto" not in expansion:
        raise ValueError("AiZynthFinder config has no uspto expansion policy")
    if "uspto" not in filters:
        raise ValueError("AiZynthFinder config has no uspto filter policy")
    if "zinc" not in stocks:
        raise ValueError("AiZynthFinder config has no zinc stock")

    policy_resources: dict[str, dict] = {}
    for policy_name, values in expansion.items():
        if not isinstance(values, list) or len(values) != 2:
            raise ValueError(
                f"expansion.{policy_name} must contain model and template paths"
            )
        model_path = _resource_path(values[0], config_path)
        template_path = _resource_path(values[1], config_path)
        session = ort.InferenceSession(
            str(model_path), providers=["CPUExecutionProvider"]
        )
        with gzip.open(template_path, "rt", encoding="utf-8") as handle:
            header = handle.readline().strip()
        if not header:
            raise ValueError(f"empty reaction-template table: {template_path}")
        policy_resources[str(policy_name)] = {
            "model": str(model_path),
            "model_inputs": [item.name for item in session.get_inputs()],
            "templates": str(template_path),
        }

    filter_path = _resource_path(filters["uspto"], config_path)
    filter_session = ort.InferenceSession(
        str(filter_path), providers=["CPUExecutionProvider"]
    )
    stock_path = _resource_path(stocks["zinc"], config_path)
    with pd.HDFStore(str(stock_path), mode="r") as store:
        keys = store.keys()
        if "/table" not in keys:
            raise ValueError(f"ZINC stock has no /table dataset: {stock_path}")
        # The official stock is a fixed-format DataFrame. ``nrows`` on its
        # FrameFixed storer is None, while axis1 carries the 35M-row dimension
        # as HDF metadata and can be read without materializing the SMILES.
        axis = store._handle.get_node("/table/axis1")
        row_count = int(axis.nrows)
    if row_count <= 0:
        raise ValueError(f"ZINC stock table is empty: {stock_path}")

    return {
        "version": importlib.metadata.version("aizynthfinder"),
        "config_path": str(config_path),
        "expansion": policy_resources,
        "filter_model": str(filter_path),
        "filter_inputs": [item.name for item in filter_session.get_inputs()],
        "stock": str(stock_path),
        "stock_rows": row_count,
    }

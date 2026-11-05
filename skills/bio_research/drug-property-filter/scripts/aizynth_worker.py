#!/usr/bin/env python
"""AiZynthFinder batch worker executed inside the dedicated ``aizynth`` env."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any


def _json_value(value: Any) -> Any:
    if value is None:
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return str(value)


def _resource_fingerprints(config_path: Path) -> dict[str, dict]:
    """Fingerprint resources without rereading a multi-GB stock every batch."""
    import yaml

    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    paths: set[Path] = {config_path.resolve()}

    def collect(value: Any) -> None:
        if isinstance(value, dict):
            for item in value.values():
                collect(item)
        elif isinstance(value, list):
            for item in value:
                collect(item)
        elif isinstance(value, str):
            candidate = Path(value).expanduser()
            if not candidate.is_absolute():
                candidate = config_path.parent / candidate
            if candidate.is_file():
                paths.add(candidate.resolve())

    collect(config)
    fingerprints: dict[str, dict] = {}
    for path in sorted(paths, key=str):
        stat = path.stat()
        record: dict[str, Any] = {
            "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
        if stat.st_size <= 256 * 1024 * 1024:
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            record["sha256"] = digest.hexdigest()
        else:
            # The official stock is ~1.3 GB. Its checksum is verified when it is
            # installed; per-run stat metadata avoids a full extra disk scan.
            checksum_file = path.with_suffix(path.suffix + ".md5")
            if checksum_file.is_file():
                record["verified_install_md5"] = checksum_file.read_text(
                    encoding="utf-8"
                ).strip()
            record["hash_policy"] = "install-time checksum; per-run size+mtime"
        fingerprints[str(path)] = record
    return fingerprints


def _starting_material_occurrences(route: dict) -> list[dict]:
    materials: list[dict] = []

    def walk(node: Any) -> None:
        if isinstance(node, list):
            for item in node:
                walk(item)
            return
        if not isinstance(node, dict):
            return
        children = node.get("children") or []
        node_type = str(node.get("type") or "").lower()
        smiles = node.get("smiles")
        if node_type in {"mol", "molecule"} and smiles and not children:
            materials.append({
                "smiles": str(smiles),
                "in_stock": bool(node.get("in_stock")),
            })
        walk(children)

    walk(route)
    return materials


def _starting_materials(route: dict) -> list[dict]:
    """Return unique starting-material identities for the human route report."""
    unique: dict[str, dict] = {}
    for material in _starting_material_occurrences(route):
        unique.setdefault(material["smiles"], material)
    return list(unique.values())


def _route_is_solved(route: dict) -> bool:
    """A serialized route is solved only when every leaf molecule is in stock."""
    materials = _starting_material_occurrences(route)
    return bool(materials) and all(item["in_stock"] for item in materials)


def _reaction_count(route: dict) -> int:
    """Count reaction nodes in one exported route tree."""
    count = 0

    def walk(node: Any) -> None:
        nonlocal count
        if isinstance(node, list):
            for item in node:
                walk(item)
            return
        if not isinstance(node, dict):
            return
        if str(node.get("type") or "").casefold() == "reaction":
            count += 1
        walk(node.get("children") or [])

    walk(route)
    return count


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _analysis_scorer_names(analysis: Any) -> list[str]:
    """Return the scorers actually selected for route analysis."""
    scorers = getattr(analysis, "scorers", None)
    return [repr(scorer) for scorer in scorers] if isinstance(scorers, list) else []


def _route_metrics(
    *,
    routes: list[dict],
    route_scores: list[dict],
    scorer_names: list[str],
) -> dict:
    """Publish only scalar route evidence with an explicit backend meaning."""
    solved_indices = [
        index for index, route in enumerate(routes) if _route_is_solved(route)
    ]
    selected_index = solved_indices[0] if solved_indices else (0 if routes else None)
    selected_route = routes[selected_index] if selected_index is not None else None
    materials = (
        _starting_material_occurrences(selected_route) if selected_route else []
    )
    normalized_scorers = [name.strip().casefold() for name in scorer_names]
    selected_scores = (
        route_scores[selected_index]
        if selected_index is not None and selected_index < len(route_scores)
        else {}
    )
    state_score = (
        _finite_number(selected_scores.get(scorer_names[0]))
        if normalized_scorers == ["state score"]
        else None
    )
    return {
        "route_found": 1.0 if solved_indices else 0.0,
        "selected_route_index": selected_index,
        "selected_route_rank": selected_index + 1 if selected_index is not None else None,
        "state_score": state_score,
        "state_score_applicable": state_score is not None,
        "state_score_scorers": scorer_names,
        "reaction_count": _reaction_count(selected_route) if selected_route else None,
        "precursor_count": len(materials) if selected_route else None,
        "precursors_in_stock_count": (
            sum(bool(item.get("in_stock")) for item in materials)
            if selected_route
            else None
        ),
        "precursors_in_stock_fraction": (
            sum(bool(item.get("in_stock")) for item in materials) / len(materials)
            if materials
            else None
        ),
    }


def _select_component(component, requested: list[str], label: str) -> None:
    available = [str(item) for item in component.items]
    missing = [item for item in requested if item not in available]
    if missing:
        raise RuntimeError(
            f"AiZynthFinder config lacks requested {label}: {missing}; "
            f"available={available}"
        )
    component.select(requested)


def _select_loaded_backends(finder, payload: dict) -> None:
    if not finder.expansion_policy.items:
        raise RuntimeError("AiZynthFinder config loaded no expansion policy")
    if not finder.stock.items:
        raise RuntimeError("AiZynthFinder config loaded no stock")
    expansion = [
        str(item) for item in (payload.get("expansion_policies") or ["uspto"])
    ]
    filters = [
        str(item) for item in (payload.get("filter_policies") or ["uspto"])
    ]
    stocks = [str(item) for item in (payload.get("stocks") or ["zinc"])]
    _select_component(finder.expansion_policy, expansion, "expansion policies")
    if finder.filter_policy.items:
        _select_component(finder.filter_policy, filters, "filter policies")
    elif filters:
        raise RuntimeError("AiZynthFinder config loaded no filter policy")
    _select_component(finder.stock, stocks, "stocks")


def plan(payload: dict) -> dict:
    from aizynthfinder.aizynthfinder import AiZynthFinder
    from aizynthfinder.analysis import RouteSelectionArguments

    entries = payload.get("candidates")
    if not isinstance(entries, list) or not entries:
        raise ValueError("candidates must be a non-empty list")

    config_path = Path(str(payload.get("config_path") or "")).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"AiZynthFinder config not found: {config_path}")

    candidate_ids: set[str] = set()
    for position, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"candidate {position} must be an object")
        candidate_id = str(entry.get("candidate_id") or "")
        smiles = str(entry.get("smiles") or "")
        if not candidate_id or not smiles:
            raise ValueError(f"candidate {position} needs candidate_id and smiles")
        if candidate_id in candidate_ids:
            raise ValueError(f"duplicate candidate_id: {candidate_id}")
        candidate_ids.add(candidate_id)

    search_time_limit = max(1, int(payload.get("search_time_limit_seconds", 60)))
    iteration_limit = max(1, int(payload.get("iteration_limit", 100)))
    max_routes = max(1, int(payload.get("max_routes", 5)))
    return_first = bool(payload.get("return_first", True))

    finder = AiZynthFinder(configfile=str(config_path))
    _select_loaded_backends(finder, payload)
    scorer_names: list[str] = []
    finder.config.search.time_limit = search_time_limit
    finder.config.search.iteration_limit = iteration_limit
    finder.config.search.return_first = return_first
    finder.config.post_processing.min_routes = 1
    finder.config.post_processing.max_routes = max_routes

    results = []
    for entry in entries:
        candidate_id = str(entry["candidate_id"])
        smiles = str(entry["smiles"])
        try:
            finder.target_smiles = smiles
            finder.tree_search(show_progress=False)
            finder.build_routes(
                selection=RouteSelectionArguments(
                    nmin=1,
                    nmax=max_routes,
                    return_all=False,
                )
            )
            scorer_names = _analysis_scorer_names(finder.analysis)
            statistics = _json_value(finder.extract_statistics())
            routes = _json_value(list(finder.routes.dicts))
            route_scores = _json_value(list(finder.routes.scores))
            metrics = _route_metrics(
                routes=routes,
                route_scores=route_scores,
                scorer_names=scorer_names,
            )
            is_solved = metrics["route_found"] == 1.0
            status = "solved" if is_solved else "not_found_within_budget"
            route_records = [
                {
                    "rank": index + 1,
                    "tree": route,
                    "is_solved": _route_is_solved(route),
                    "scores": route_scores[index] if index < len(route_scores) else {},
                    "starting_materials": _starting_materials(route),
                }
                for index, route in enumerate(routes)
            ]
            results.append({
                "candidate_id": candidate_id,
                "smiles": smiles,
                "status": status,
                "is_solved": is_solved,
                "route_count": len(routes),
                "statistics": statistics,
                "metrics": metrics,
                "routes": route_records,
            })
        except Exception as exc:
            results.append({
                "candidate_id": candidate_id,
                "smiles": smiles,
                "status": "error",
                "is_solved": False,
                "route_count": 0,
                "metrics": {},
                "error": f"{type(exc).__name__}: {exc}",
                "routes": [],
            })

    error_count = sum(result.get("status") == "error" for result in results)
    success = error_count < len(results)
    return {
        "success": success,
        "error": (
            None
            if success
            else f"AiZynthFinder failed for all {len(results)} candidates"
        ),
        "batch_status": (
            "complete" if error_count == 0
            else "failed" if error_count == len(results)
            else "partial"
        ),
        "error_count": error_count,
        "backend": {
            "id": "aizynthfinder",
            "version": importlib.metadata.version("aizynthfinder"),
            "python": platform.python_version(),
            "config_path": str(config_path),
            "expansion_policies": list(finder.expansion_policy.selection),
            "filter_policies": list(finder.filter_policy.selection),
            "stocks": list(finder.stock.selection),
            "scorers": scorer_names,
            "resource_fingerprints": _resource_fingerprints(config_path),
        },
        "search_budget": {
            "time_limit_seconds_per_candidate": search_time_limit,
            "iteration_limit_per_candidate": iteration_limit,
            "max_routes_per_candidate": max_routes,
            "return_first": return_first,
        },
        "requested": len(entries),
        "returned": len(results),
        "results": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="AiZynthFinder batch JSON worker")
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output-file", required=True)
    args = parser.parse_args()
    output_path = Path(args.output_file)
    try:
        payload = json.loads(Path(args.input_file).read_text(encoding="utf-8"))
        result = plan(payload)
        exit_code = 0 if result.get("success") else 1
    except Exception as exc:
        result = {
            "success": False,
            "error": f"{type(exc).__name__}: {exc}",
            "backend": {
                "id": "aizynthfinder",
                "python": platform.python_version(),
            },
        }
        exit_code = 1
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    print(json.dumps({"success": result["success"], "output_file": str(output_path)}))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())

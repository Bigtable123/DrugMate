"""Manage one report's ``assets/`` — 3Dmol.js + Mol* offline bundles.

Each generated report gets its own snapshot so later rounds cannot mutate the
assets or manifest sealed by an earlier report.

Pure stdlib; graceful-skip when source files are missing.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any


REQUIRED_ASSETS = {
    "3Dmol-min.js": {
        "purpose": "per-molecule 3D card viewer",
        "cdn_hint": "https://3Dmol.org/build/3Dmol-min.js",
        "license_file": "3Dmol-min.js.LICENSE.txt",
    },
    "molstar.js": {
        "purpose": "main pocket viewer (top of report)",
        "cdn_hint": "https://cdn.jsdelivr.net/npm/molstar/build/viewer/molstar.js",
        "license_file": "molstar.LICENSE.txt",
    },
    "molstar.css": {
        "purpose": "Mol* stylesheet",
        "cdn_hint": "https://cdn.jsdelivr.net/npm/molstar/build/viewer/molstar.css",
        "license_file": "molstar.LICENSE.txt",
    },
}


def ensure_assets(
    project_root: Path,
    source_dir: Path | None = None,
    bundled_asset_dir: Path | None = None,
) -> dict[str, Any]:
    """Ensure this report root's ``assets/`` contains offline viewer bundles.

    Search order for each asset:
      1. `{source_dir}/<file>` when source_dir is provided
      2. `{bundled_asset_dir}/<file>` (report-generator/assets_src/, if present)
      3. Existing `{project_root}/assets/<file>` (from a previous run)

    Missing assets are recorded but do not fail the report — the HTML uses
    graceful-skip placeholders. Copy notices from each selected asset's source;
    another source's notice is usable only when its asset bytes are identical.
    Missing custom notices are reported without changing viewer selection.
    """
    project_root = Path(project_root).expanduser().resolve()
    assets_dir = project_root / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = assets_dir / "assets_manifest.json"
    try:
        previous_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(previous_manifest, dict):
            previous_manifest = {}
    except (OSError, ValueError):
        previous_manifest = {}

    manifest = {
        "assets_dir": str(assets_dir),
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "source_dir": str(source_dir) if source_dir else None,
        "bundled_asset_dir": str(bundled_asset_dir) if bundled_asset_dir else None,
        "assets": {},
        "licenses": {},
        "missing": [],
    }

    for filename, meta in REQUIRED_ASSETS.items():
        target = assets_dir / filename
        source = _resolve_source(filename, source_dir, bundled_asset_dir, target)

        if source is None:
            manifest["assets"][filename] = {
                "present": False,
                "purpose": meta["purpose"],
                "cdn_hint": meta["cdn_hint"],
            }
            manifest["missing"].append({
                "file": filename,
                "reason": "no source found in source_dir, bundled dir, or existing assets/",
                "cdn_hint": meta["cdn_hint"],
            })
            continue

        # Only copy when target absent or content differs
        if not target.exists() or _sha256(target) != _sha256(source):
            if source.resolve() != target.resolve():
                shutil.copyfile(source, target)

        manifest["assets"][filename] = {
            "present": True,
            "purpose": meta["purpose"],
            "path": str(target),
            "size_bytes": target.stat().st_size,
            "sha256": _sha256(target),
        }
        license_name = meta["license_file"]
        license_source = _resolve_license_source(
            filename, license_name, source, source_dir, bundled_asset_dir,
            assets_dir, previous_manifest,
        )
        if license_source is None:
            manifest["assets"][filename]["license_file"] = None
            manifest["missing"].append({
                "file": license_name,
                "asset": filename,
                "reason": "no accompanying notice for the selected viewer build",
            })
            continue
        digest = _sha256(license_source)
        existing = manifest["licenses"].get(license_name)
        if existing and existing["sha256"] != digest:
            # A custom Mol* JS and bundled CSS can carry different notices.
            # Keep both instead of overwriting the notice already assigned.
            license_name = f"{Path(license_name).stem}-{digest[:12]}.txt"
        license_target = assets_dir / license_name
        if not license_target.exists() or _sha256(license_target) != digest:
            shutil.copyfile(license_source, license_target)
        manifest["licenses"][license_name] = {
            "present": True,
            "path": str(license_target),
            "size_bytes": license_target.stat().st_size,
            "sha256": digest,
        }
        manifest["assets"][filename]["license_file"] = license_name

    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def relative_asset_path(html_path: Path, project_root: Path, asset_name: str) -> str:
    """Compute relative href from an HTML file to `{project_root}/assets/<asset_name>`."""
    asset = (Path(project_root) / "assets" / asset_name).resolve()
    try:
        rel = _relpath(asset, start=Path(html_path).parent.resolve())
    except ValueError:
        rel = str(asset)
    return rel


def asset_available(manifest: dict, filename: str) -> bool:
    return bool(((manifest or {}).get("assets") or {}).get(filename, {}).get("present"))


# ---------- helpers ----------

def _resolve_license_source(
    filename: str,
    license_name: str,
    source: Path,
    source_dir: Path | None,
    bundled_asset_dir: Path | None,
    assets_dir: Path,
    previous_manifest: dict,
) -> Path | None:
    previous = {}
    if source.parent.resolve() == assets_dir.resolve():
        previous_assets = previous_manifest.get("assets")
        if isinstance(previous_assets, dict) and isinstance(previous_assets.get(filename), dict):
            previous = previous_assets[filename]
        previous_name = previous.get("license_file")
        if isinstance(previous_name, str) and Path(previous_name).name == previous_name:
            notice = assets_dir / previous_name
            previous_licenses = previous_manifest.get("licenses")
            recorded = previous_licenses.get(previous_name) if isinstance(previous_licenses, dict) else None
            if (previous.get("sha256") == _sha256(source) and notice.is_file()
                    and isinstance(recorded, dict) and recorded.get("sha256") == _sha256(notice)):
                return notice
    accompanying = source.parent / license_name
    # A leftover notice from a replaced viewer is not its new build's notice.
    if not previous and accompanying.is_file() and accompanying.stat().st_size > 0:
        return accompanying
    for directory in (source_dir, bundled_asset_dir):
        if directory is None:
            continue
        asset, notice = Path(directory) / filename, Path(directory) / license_name
        if (asset.is_file() and notice.is_file() and notice.stat().st_size > 0
                and _sha256(asset) == _sha256(source)):
            return notice
    return None


def _resolve_source(
    filename: str,
    source_dir: Path | None,
    bundled_asset_dir: Path | None,
    existing_target: Path,
) -> Path | None:
    if source_dir:
        candidate = Path(source_dir) / filename
        if candidate.exists():
            return candidate
    if bundled_asset_dir:
        candidate = Path(bundled_asset_dir) / filename
        if candidate.exists():
            return candidate
    if existing_target.exists() and existing_target.stat().st_size > 0:
        return existing_target
    return None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _relpath(target: Path, start: Path) -> str:
    import os
    return os.path.relpath(str(target), start=str(start))

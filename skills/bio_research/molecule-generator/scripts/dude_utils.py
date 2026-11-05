"""DUD-E decoy download utilities.

Downloads pre-computed decoys from the DUD-E database (dude.docking.org).
Uses a shared cache directory under SBDD_DATA_ROOT/dude_cache/ so that
each target's decoys are only downloaded once across all experiments.
Falls back gracefully if the target is not available.
"""
from __future__ import annotations

import os
import random
import urllib.error
import urllib.request
from pathlib import Path


DUDE_BASE_URL = "https://dude.docking.org/targets"
_REQUEST_TIMEOUT = 30
_TARGET_ALIASES = {
    "cyclin-dependent kinase 2": "cdk2",
}

REPO_ROOT = Path(__file__).resolve().parents[4]
DATA_ROOT = Path(os.environ.get("SBDD_DATA_ROOT", REPO_ROOT.parent / "data")).expanduser()
DUDE_CACHE_DIR = DATA_ROOT / "dude_cache"


def _resolve_target_name(target_name: str) -> str:
    """Resolve known full target names while preserving DUD-E slug inputs."""
    normalized = target_name.strip().lower()
    return _TARGET_ALIASES.get(normalized, normalized)


def _get_cache_path(target_name: str) -> Path:
    """Return the cache file path for a given target."""
    return DUDE_CACHE_DIR / _resolve_target_name(target_name) / "decoys_final.ism"


def _is_cached(target_name: str) -> bool:
    """Check if decoys for this target are already cached locally."""
    cache_file = _get_cache_path(target_name)
    return cache_file.exists() and cache_file.stat().st_size > 0


def check_dude_available(target_name: str) -> bool:
    """Check if DUD-E decoys are available from cache or the DUD-E website.

    Accepts DUD-E slugs (e.g., "cdk2", "egfr", "aa2ar") and known full names.
    """
    target_lower = _resolve_target_name(target_name)
    if _is_cached(target_lower):
        return True

    url = f"{DUDE_BASE_URL}/{target_lower}/decoys_final.ism"
    try:
        req = urllib.request.Request(url, method="HEAD")
        resp = urllib.request.urlopen(req, timeout=10)
        return resp.status == 200
    except Exception:
        return False


def _download_to_cache(target_name: str) -> Path | None:
    """Download DUD-E .ism file and save to cache. Returns cache path or None."""
    target_lower = _resolve_target_name(target_name)
    cache_file = _get_cache_path(target_lower)
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    ism_url = f"{DUDE_BASE_URL}/{target_lower}/decoys_final.ism"
    try:
        req = urllib.request.Request(ism_url, headers={"Accept": "text/plain"})
        resp = urllib.request.urlopen(req, timeout=_REQUEST_TIMEOUT)
        content = resp.read()
        with open(cache_file, "wb") as f:
            f.write(content)
        return cache_file
    except Exception:
        if cache_file.exists():
            cache_file.unlink()
        return None


def download_dude_decoys(
    target_name: str,
    output_dir: Path,
    max_decoys: int = 200,
    seed: int = 42,
) -> dict:
    """Get DUD-E decoys for a target (from cache or download).

    Cache location: SBDD_DATA_ROOT/dude_cache/{target}/decoys_final.ism
    Once downloaded, subsequent runs read from cache without network access.
    """
    target_lower = _resolve_target_name(target_name)

    # Read from cache or download
    cache_file = _get_cache_path(target_lower)
    from_cache = False

    if _is_cached(target_lower):
        from_cache = True
    else:
        result_path = _download_to_cache(target_lower)
        if result_path is None:
            return {"success": False, "error": f"DUD-E download failed for target '{target_lower}'", "decoys": []}

    try:
        content = cache_file.read_text(encoding="utf-8")
    except UnicodeDecodeError as e:
        return {"success": False, "error": f"DUD-E file is not valid UTF-8: {e}", "decoys": []}
    except Exception as e:
        return {"success": False, "error": f"Failed to read cache file: {e}", "decoys": []}

    lines = [line.strip() for line in content.splitlines() if line.strip()]
    if not lines:
        return {"success": False, "error": "DUD-E file is empty", "decoys": []}

    all_decoys = []
    for line in lines:
        parts = line.split()
        if not parts:
            continue
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else f"DUDE_{len(all_decoys):04d}"
        all_decoys.append({"smiles": smiles, "name": name, "source": "dude"})

    if len(all_decoys) > max_decoys:
        rng = random.Random(seed)
        all_decoys = rng.sample(all_decoys, max_decoys)

    # Also write the sampled subset to output_dir for provenance
    output_dir.mkdir(parents=True, exist_ok=True)
    ism_file = output_dir / "dude_decoys_sampled.ism"
    with open(ism_file, "w") as f:
        for d in all_decoys:
            f.write(f"{d['smiles']} {d['name']}\n")

    return {
        "success": True,
        "decoys": all_decoys,
        "source": "dude",
        "count": len(all_decoys),
        "total_available": len(lines),
        "target_dude_name": target_lower,
        "from_cache": from_cache,
        "cache_file": str(cache_file),
        "file": str(ism_file),
    }

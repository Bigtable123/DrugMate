"""Optional PyMOL static cover renderer.

PyMOL open-source is a heavy dependency (needs EGL/OSMesa or X server) and
frequently unavailable in production containers. This module tries the two
common entry points and graceful-skips on any failure — a report will still
render without the PyMOL cover image.

Only invoked when `enable_pymol_cover=true`.
"""
from __future__ import annotations

import subprocess
import shutil
import sys
import tempfile
from pathlib import Path


def check_pymol_available() -> tuple[bool, str]:
    """Return (available, reason). Prefer `pymol` on PATH; fall back to `import pymol`."""
    if shutil.which("pymol"):
        return True, "pymol CLI on PATH"
    try:
        import pymol  # noqa: F401
        return True, "pymol python module"
    except ImportError:
        return False, "pymol not installed"


def render_cover(
    receptor_pdb_path: str | None,
    reference_ligand_sdf_path: str | None,
    pocket_residues: list[str] | None,
    output_dir: Path,
    width: int = 1024,
    height: int = 640,
) -> tuple[str | None, str | None]:
    """Render `figures/pymol_cover.png`.

    Returns (relative_path, error). One is None.
    """
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_png = figures_dir / "pymol_cover.png"

    if not receptor_pdb_path or not Path(receptor_pdb_path).exists():
        return None, "receptor pdb missing"

    available, why = check_pymol_available()
    if not available:
        return None, why

    resi_selection = _resi_selection(pocket_residues or [])
    ligand_load = ""
    if reference_ligand_sdf_path and Path(reference_ligand_sdf_path).exists():
        ligand_load = (
            f"load {reference_ligand_sdf_path}, refligand\n"
            "hide everything, refligand\n"
            "show sticks, refligand\n"
            "color yellow, refligand and elem C\n"
        )

    pocket_show = ""
    if resi_selection:
        pocket_show = (
            f"select pocket_res, receptor and resi {resi_selection}\n"
            "show sticks, pocket_res\n"
            "color orange, pocket_res and elem C\n"
        )

    pml = f"""
bg_color white
load {receptor_pdb_path}, receptor
hide everything, receptor
show cartoon, receptor
color grey70, receptor
{pocket_show}
{ligand_load}
orient
zoom refligand or pocket_res, buffer=6.0
ray {width}, {height}
png {out_png}, dpi=150
"""

    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as f:
        f.write(pml)
        pml_path = f.name

    try:
        if shutil.which("pymol"):
            proc = subprocess.run(
                ["pymol", "-cq", pml_path],
                capture_output=True,
                text=True,
                timeout=300,
            )
            if proc.returncode != 0:
                return None, f"pymol exited {proc.returncode}: {(proc.stderr or '').strip()[:400]}"
        else:
            # Fallback: run in-process via pymol module (headless -cq behavior)
            try:
                from pymol import cmd  # noqa: F401
                subprocess.run(
                    [sys.executable, "-c", f"from pymol import cmd; cmd.run('{pml_path}')"],
                    capture_output=True,
                    text=True,
                    timeout=300,
                    check=False,
                )
            except Exception as e:  # noqa: BLE001
                return None, f"pymol module invocation failed: {e}"
    except subprocess.TimeoutExpired:
        return None, "pymol render timed out (300s)"
    except Exception as e:  # noqa: BLE001
        return None, f"pymol invocation failed: {e}"
    finally:
        Path(pml_path).unlink(missing_ok=True)

    if not out_png.exists():
        return None, "pymol produced no output PNG"
    return "figures/pymol_cover.png", None


def _resi_selection(residues: list[str]) -> str:
    import re
    resi_pattern = re.compile(r"(\d+)")
    ids: list[int] = []
    for r in residues:
        if r is None:
            continue
        m = resi_pattern.search(str(r))
        if m:
            try:
                ids.append(int(m.group(1)))
            except ValueError:
                continue
    return "+".join(str(i) for i in sorted(set(ids)))

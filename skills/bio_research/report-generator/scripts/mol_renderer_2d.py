"""RDKit 2D structure renderer — graceful-skip when RDKit is unavailable.

Produces base64-encoded PNGs (embedded in the HTML) so a report is
self-contained without additional file hops for 2D pictures.
"""
from __future__ import annotations

import base64
import io
from typing import Any


def check_rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
        return True
    except ImportError:
        return False


def render_smiles_png_base64(
    smiles: str,
    size: tuple[int, int] = (280, 220),
) -> tuple[str | None, str | None]:
    """Render a SMILES to a base64 PNG data URI.

    Returns (data_uri, error). One of the two is None.
    """
    if not smiles:
        return None, "empty smiles"
    if not check_rdkit_available():
        return None, "rdkit not installed"

    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, f"rdkit failed to parse smiles: {smiles[:80]}"

        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.bondLineWidth = 1.5
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png_bytes = drawer.GetDrawingText()

        b64 = base64.b64encode(png_bytes).decode("ascii")
        return f"data:image/png;base64,{b64}", None
    except Exception as e:  # noqa: BLE001 — renderer must never crash the report
        # Cairo backend not always available — fall back to a bitmap via Pillow-free path
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, f"rdkit failed to parse smiles after cairo error: {e}"
            img = Draw.MolToImage(mol, size=size)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            b64 = base64.b64encode(buf.getvalue()).decode("ascii")
            return f"data:image/png;base64,{b64}", None
        except Exception as e2:  # noqa: BLE001
            return None, f"rdkit draw failed: {e} / fallback: {e2}"


def render_candidates(candidates: list[dict[str, Any]], size: tuple[int, int] = (280, 220)) -> dict[str, Any]:
    """Attach base64 2D PNGs to each candidate. Returns a status dict."""
    if not check_rdkit_available():
        for c in candidates:
            c["structure_2d_png"] = None
            c["structure_2d_error"] = "rdkit not installed"
        return {"available": False, "reason": "rdkit not installed", "rendered": 0, "failed": len(candidates)}

    rendered = 0
    failed = 0
    for c in candidates:
        smiles = c.get("smiles")
        data_uri, err = render_smiles_png_base64(smiles, size=size) if smiles else (None, "no smiles")
        c["structure_2d_png"] = data_uri
        c["structure_2d_error"] = err
        if data_uri:
            rendered += 1
        else:
            failed += 1
    return {"available": True, "rendered": rendered, "failed": failed}

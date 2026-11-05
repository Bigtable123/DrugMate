"""3D structure renderers.

Two flavors:
  1. per-candidate 3Dmol.js viewer HTML fragments (small, one per card)
  2. Mol* main viewer HTML (pocket + reference ligand + optional top-N overlay)

Pure stdlib. RDKit is used only if available for pose extraction from an SDF
containing multiple candidates; otherwise a text-splitter fallback is used.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any


# ---------- 3Dmol.js per-molecule cards ----------

CARD_TEMPLATE = """<div id="viewer-{viewer_id}" class="viewer-3dmol" style="position: relative; width: {width}px; height: {height}px;"></div>
<script>
(function() {{
  function boot() {{
    if (typeof $3Dmol === 'undefined') {{
      var el = document.getElementById('viewer-{viewer_id}');
      if (el) el.innerHTML = '<div class="viewer-error">3Dmol.js asset missing (see assets/assets_manifest.json)</div>';
      return;
    }}
    var receptorPdb = {receptor_json};
    var poseSdf = {pose_json};
    var v = $3Dmol.createViewer(document.getElementById('viewer-{viewer_id}'), {{ backgroundColor: 'white' }});
    if (receptorPdb) {{
      v.addModel(receptorPdb, 'pdb');
      v.setStyle({{ model: 0 }}, {{ cartoon: {{ color: 'spectrum', opacity: 0.85 }} }});
      {pocket_highlight_js}
    }}
    if (poseSdf) {{
      var poseModel = v.addModel(poseSdf, 'sdf');
      v.setStyle({{ model: poseModel }}, {{ stick: {{ colorscheme: 'greenCarbon', radius: 0.18 }} }});
      v.zoomTo({{ model: poseModel }});
    }} else {{
      v.zoomTo();
    }}
    v.render();
  }}
  if (document.readyState === 'loading') {{
    document.addEventListener('DOMContentLoaded', boot);
  }} else {{
    boot();
  }}
}})();
</script>"""


def render_candidate_card_viewer(
    viewer_id: str,
    receptor_pdb_path: str | None,
    pose_sdf_text: str | None,
    pocket_residues: list[str] | None = None,
    width: int = 380,
    height: int = 300,
) -> tuple[str, str | None]:
    """Return (html_fragment, error).

    Fragment references a global `$3Dmol` (loaded once at the top of the report).
    Data is embedded inline as JSON-escaped strings so no file fetch is needed.
    """
    import json as _json

    if not receptor_pdb_path or not Path(receptor_pdb_path).exists():
        return _placeholder(viewer_id, "receptor PDB missing", width, height), "receptor pdb missing"

    try:
        receptor_pdb_text = Path(receptor_pdb_path).read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        return _placeholder(viewer_id, f"could not read receptor pdb: {e}", width, height), str(e)

    pocket_highlight_js = ""
    if pocket_residues:
        resi_list = _residues_to_resi_list(pocket_residues)
        if resi_list:
            pocket_highlight_js = (
                f"v.setStyle({{ model: 0, resi: {resi_list} }}, {{ cartoon: {{ color: 'orange' }}, "
                "stick: { colorscheme: 'yellowCarbon', radius: 0.15 } });"
            )

    html = CARD_TEMPLATE.format(
        viewer_id=viewer_id,
        width=width,
        height=height,
        receptor_json=_json.dumps(receptor_pdb_text),
        pose_json=_json.dumps(pose_sdf_text) if pose_sdf_text else "null",
        pocket_highlight_js=pocket_highlight_js,
    )
    return html, None


def extract_pose_sdf(sdf_path: str | None, mol_name: str, pose_index: int | None = None) -> str | None:
    """Extract a single pose block from a multi-molecule SDF.

    Preference:
      1. Match by mol_name (SDF `_Name` header line, first line of each block)
      2. Match by pose_index using RDKit's SDMolSupplier ordering
    """
    if not sdf_path:
        return None
    path = Path(sdf_path)
    if not path.exists():
        return None

    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None

    # SDF blocks are separated by $$$$ ; the first line is the mol name
    blocks = [b for b in text.split("$$$$") if b.strip()]
    for i, block in enumerate(blocks):
        first_line = (block.lstrip("\n\r").splitlines() or [""])[0].strip()
        if first_line == mol_name:
            return block.strip() + "\n$$$$\n"

    # Fallback by pose_index if provided
    if pose_index is not None and 0 <= pose_index < len(blocks):
        return blocks[pose_index].strip() + "\n$$$$\n"

    return None


def render_candidate_cards(
    candidates: list[dict[str, Any]],
    receptor_pdb_path: str | None,
    selected_sdf_path: str | None,
    docked_sdf_path: str | None,
    pocket_residues: list[str] | None = None,
    width: int = 380,
    height: int = 300,
) -> dict[str, Any]:
    """Attach a `viewer_3d_html` string to every candidate.

    Prefers `selected_poses.sdf` (post-selection); falls back to `docked_poses.sdf`.
    """
    rendered = 0
    failed = 0
    for c in candidates:
        mol_name = c.get("mol_name")
        pose_index = c.get("pose_index")
        viewer_id = _safe_viewer_id(mol_name or f"pose_{pose_index}")

        pose_sdf = None
        for sdf_path in (selected_sdf_path, docked_sdf_path):
            pose_sdf = extract_pose_sdf(sdf_path, mol_name, pose_index) if sdf_path else None
            if pose_sdf:
                break

        html, err = render_candidate_card_viewer(
            viewer_id=viewer_id,
            receptor_pdb_path=receptor_pdb_path,
            pose_sdf_text=pose_sdf,
            pocket_residues=pocket_residues,
            width=width,
            height=height,
        )
        c["viewer_3d_html"] = html
        c["viewer_3d_error"] = err
        c["viewer_3d_has_pose"] = bool(pose_sdf)
        if err is None and pose_sdf:
            rendered += 1
        else:
            failed += 1

    return {"rendered": rendered, "failed": failed, "total": len(candidates)}


# ---------- Mol* main viewer ----------

MOLSTAR_MAIN_TEMPLATE = """<div id="molstar-main" class="molstar-main-viewer" style="position: relative; width: 100%; height: {height}px;"></div>
<script>
(function() {{
  function boot() {{
    var container = document.getElementById('molstar-main');
    if (typeof molstar === 'undefined') {{
      container.innerHTML = '<div class="viewer-error">Mol* asset missing (see assets/assets_manifest.json)</div>';
      return;
    }}
    molstar.Viewer.create(container, {{
      layoutIsExpanded: false,
      layoutShowControls: {show_controls},
      layoutShowRemoteState: false,
      layoutShowSequence: true,
      layoutShowLog: false,
      layoutShowLeftPanel: {show_left},
      viewportShowExpand: true,
      viewportShowSelectionMode: false,
      viewportShowAnimation: false,
    }}).then(function(viewer) {{
      var receptorUrl = {receptor_url_json};
      var ligandUrl = {ligand_url_json};
      var status = document.createElement('div');
      status.className = 'viewer-error';
      status.style.position = 'absolute';
      status.style.left = '8px';
      status.style.bottom = '8px';
      status.style.background = 'rgba(255,255,255,0.88)';
      status.style.display = 'none';
      container.appendChild(status);
      function showStatus(message) {{
        status.textContent = message;
        status.style.display = 'block';
      }}
      var chain = receptorUrl
        ? viewer.loadStructureFromUrl(receptorUrl, 'pdb')
        : Promise.resolve();
      chain.then(function() {{
        if (!ligandUrl) return null;
        return viewer.loadStructureFromUrl(ligandUrl, 'sdf').catch(function(err) {{
          showStatus('Reference ligand load warning: ' + err);
        }});
      }}).catch(function(err) {{
        container.innerHTML = '<div class="viewer-error">Mol* receptor load error: ' + err + '</div>';
      }});
    }}).catch(function(err) {{
      container.innerHTML = '<div class="viewer-error">Mol* init error: ' + err + '</div>';
    }});
  }}
  if (document.readyState === 'loading') {{
    document.addEventListener('DOMContentLoaded', boot);
  }} else {{
    boot();
  }}
}})();
</script>"""


def render_molstar_main(
    receptor_pdb_path: str | None,
    reference_ligand_sdf_path: str | None,
    html_output_path: Path,
    height: int = 520,
    show_controls: bool = True,
    show_left_panel: bool = False,
) -> tuple[str, str | None]:
    """Mol* main viewer. Loads receptor + reference ligand via relative URLs.

    Mol* prefers URL loading over inline text for large PDBs, so we compute
    relative paths from the HTML file to the receptor and reference SDF.
    """
    import json as _json

    def _rel_or_none(path_str: str | None) -> str | None:
        if not path_str:
            return None
        p = Path(path_str)
        if not p.exists():
            return None
        try:
            import os
            return os.path.relpath(str(p.resolve()), start=str(html_output_path.parent.resolve()))
        except ValueError:
            return str(p.resolve())

    receptor_rel = _rel_or_none(receptor_pdb_path)
    ligand_rel = _rel_or_none(reference_ligand_sdf_path)

    if not receptor_rel:
        return _placeholder("molstar-main", "receptor PDB missing", 800, height), "receptor pdb missing"

    html = MOLSTAR_MAIN_TEMPLATE.format(
        height=height,
        show_controls=_js_bool(show_controls),
        show_left=_js_bool(show_left_panel),
        receptor_url_json=_json.dumps(receptor_rel),
        ligand_url_json=_json.dumps(ligand_rel) if ligand_rel else "null",
    )
    return html, None


# ---------- helpers ----------

_RESI_RE = re.compile(r"([A-Za-z]+)\s*[-_]?\s*(\d+)")


def _residues_to_resi_list(residues: list[str]) -> list[int]:
    """Convert ['Leu83', 'Val18'] → [83, 18] for 3Dmol resi selection."""
    out: list[int] = []
    for r in residues or []:
        if r is None:
            continue
        m = _RESI_RE.search(str(r))
        if m:
            try:
                out.append(int(m.group(2)))
            except ValueError:
                continue
    return sorted(set(out))


def _safe_viewer_id(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", str(name))[:64] or "viewer"


def _placeholder(viewer_id: str, message: str, width: int, height: int) -> str:
    return (
        f'<div id="viewer-{viewer_id}" class="viewer-3dmol viewer-placeholder" '
        f'style="width: {width}px; height: {height}px; border: 1px dashed #999; '
        f'display: flex; align-items: center; justify-content: center; color: #666;">'
        f'<div class="viewer-error">{message}</div>'
        f"</div>"
    )


def _js_bool(b: bool) -> str:
    return "true" if b else "false"

# Third-party notices and source provenance

Drugmate includes code from Hermes Agent and other open-source projects.
The original copyright notices and licenses are retained. The upstream
[MIT license](LICENSE) does not replace the component-specific terms below.

## Hermes Agent

The agent runtime contains code from [Hermes Agent](https://github.com/NousResearch/hermes-agent),
licensed under MIT, copyright (c) 2025 Nous Research. The full license text is
retained in [LICENSE](LICENSE).

## Bundled components

| Component | Location and origin | License and retained notices |
|---|---|---|
| Mol* viewer | `skills/bio_research/report-generator/assets_src/molstar.js` and `.css`; <https://github.com/molstar/molstar> | [Mol* and incorporated component notices](viewers/molstar.LICENSE.txt) |
| 3Dmol.js viewer | `skills/bio_research/report-generator/assets_src/3Dmol-min.js`; <https://github.com/3dmol/3Dmol.js> | [BSD-3-Clause and incorporated component notices](viewers/3Dmol-min.js.LICENSE.txt) |
| Open Babel fast-coordinate compatibility program | `skills/bio_research/input-normalization/scripts/openbabel_fast_compat/`; derived from Open Babel `openbabel-3-1-1` | [GPL version 2](../skills/bio_research/input-normalization/scripts/openbabel_fast_compat/COPYING), preserved source headers and [source/change provenance](../skills/bio_research/input-normalization/scripts/openbabel_fast_compat/README.md) |

The Mol* viewer's pre-existing vendor record identifies version 5.10.1. The
viewer bundle includes React 18.3.1, Immutable.js and attributed shader code;
the stylesheet identifies normalize.css 3.0.3. Their notices accompany the
retained report viewer. Existing embedded legal comments remain intact.

The exact version of the copied 3Dmol.js build was not recorded. Its SHA-256
is documented in the [viewer asset inventory](../skills/bio_research/report-generator/assets_src/README.md).
The companion license file contains the upstream project's full BSD/MIT notices
and the shader attribution retained in the copied build; it does not assign a
version based on the upstream website's current release.

`licenses/viewers/` archives exact copies of the two notices accompanying the
bundled viewer builds. The original companion files remain in `assets_src/` so
source and installed report generators can copy them with the viewer assets.
These archived copies do not replace or update the viewer code or its notices.

The Open Babel compatibility sources carry different per-file GPL declarations:
`gen3d_fast.cpp` and `driver.cpp` use GPL-2.0-only, while the stereo helpers
allow GPL-2.0-or-later. The combined compatibility program is distributed under
GPL version 2. Preserve its source, headers and `COPYING` when redistributing
that component. It is not relicensed by the upstream MIT file.

## Installed dependencies, tools and data

Python dependencies are described in `pyproject.toml`. Installed dependency
trees, environments and downloaded binaries are excluded from this source
release; dependencies retain their own licenses when installed.

External scientific tools, model weights, protein structures, compound libraries
and benchmark data are governed by their respective source terms. The repository
license does not grant rights to those external materials. Cite the tools,
models and datasets actually used in an experiment, using their published
references. This inventory records known bundled components; it is not a license
grant for separately obtained software or data.

When distributing generated HTML reports that include copied viewer JavaScript
or CSS, include the applicable viewer license and notice files with the reports.
The report asset manager copies the notices for the selected viewer sources and
records their hashes in `assets_manifest.json`. Custom builds without matching
notices are identified in that manifest instead of receiving another build's
notice automatically.


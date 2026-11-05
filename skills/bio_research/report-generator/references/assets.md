# Assets — 3Dmol.js and Mol* offline bundles

The report generator does not download assets at runtime. It places them under
`{project_root}/assets/` once per project so every round's `report.html` (and
the final report) references them via a relative path.

## Files required

| File | Purpose | Approximate size |
|---|---|---|
| `3Dmol-min.js` | Per-molecule 3D card viewer | ~600 KB |
| `molstar.js` | Main pocket viewer at the top of the report | ~3-4 MB |
| `molstar.css` | Mol* stylesheet | ~50 KB |

## Where they end up

```
{data_root}/{timestamp}-{slug}/
├── assets/                    ← single copy, referenced by every round
│   ├── 3Dmol-min.js
│   ├── molstar.js
│   ├── molstar.css
│   └── assets_manifest.json   ← generated on every run
├── round1/step7/report.html   ← <script src="../../assets/3Dmol-min.js">
├── round2/step7/report.html
└── final_report.html          ← <script src="assets/3Dmol-min.js">
```

## Providing the assets

`report-generator` searches in this order for each file:

1. `input_data.assets_source_dir` (if provided when invoking the skill)
2. `report-generator/assets_src/` (a directory you can pre-populate inside the skill folder)
3. An existing `{project_root}/assets/<file>` (from a previous run)

If a file is not found in any of the three, `assets_manifest.json` records it
as missing and the HTML renders a graceful-skip placeholder for that viewer.

## Downloading offline copies

Any of the following approaches works — pick one and keep the file under
`report-generator/assets_src/` (the "bundled" location) so future runs pick
it up automatically without configuration.

### 3Dmol.js

```bash
# Latest release from 3Dmol.org
curl -L https://3Dmol.org/build/3Dmol-min.js \
  -o report-generator/assets_src/3Dmol-min.js

# Or npm mirror
curl -L https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js \
  -o report-generator/assets_src/3Dmol-min.js
```

### Mol\*

```bash
# JS bundle
curl -L https://cdn.jsdelivr.net/npm/molstar/build/viewer/molstar.js \
  -o report-generator/assets_src/molstar.js

# Stylesheet
curl -L https://cdn.jsdelivr.net/npm/molstar/build/viewer/molstar.css \
  -o report-generator/assets_src/molstar.css
```

Version-pinning is recommended for reproducible reports:
`https://cdn.jsdelivr.net/npm/molstar@<version>/build/viewer/molstar.js`.

## Invocation with an explicit source directory

Replace angle-bracket path placeholders with the configured directories before
saving this JSON; JSON does not expand shell variables.

```json
{
  "mode": "round",
  "project_dir": "<SBDD_DATA_ROOT>/20260626-120000-cdk2-4kd1/round1",
  "assets_source_dir": "<REPORT_ASSETS_SOURCE_DIR>"
}
```

Assets from the configured `assets_source_dir` are copied to `{project_root}/assets/` and
recorded in `assets_manifest.json` with SHA-256 and size for provenance.

## Content Security Policy note

The report is designed to work when opened as a local file or served behind a
tight CSP. It does not fetch anything cross-origin — the two 3D viewers live
in the same directory as the HTML and access to receptor PDB / reference
ligand SDF is via relative URL from the HTML file. If you serve the reports
behind a proxy, the proxy must serve `assets/` and each `round*/step7/` (or
`final_report.html`) under the same origin.

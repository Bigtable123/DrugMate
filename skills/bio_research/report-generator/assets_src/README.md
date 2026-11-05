# Offline molecular viewers

These files support the SBDD report generator's offline molecular viewers.
Keep them in the source release; they are runtime report resources.

| File | Provenance | SHA-256 |
|---|---|---|
| `3Dmol-min.js` | 3Dmol.js; copied build version was not recorded | `e5b7daf4264660c6514ce0474e2d7b8dc582fb7186fa5def60b6be553b87ea4a` |
| `molstar.js` | Mol* 5.10.1, as identified in the existing web vendor record | `5567eb19fa8e7a7b3b161d4b96807c4db244cecd2f3e7c87f99c052b8b5b5b30` |
| `molstar.css` | Stylesheet paired with the copied Mol* viewer | `5b68ceb6d3642549b4e9b2c071e58e41b98a5350ae269180587b39da86925d55` |

`3Dmol-min.js.LICENSE.txt` contains 3Dmol.js's upstream BSD-3-Clause license,
the GLmol/Three.js/jQuery notices supplied by that project, and the MIT shader
attribution preserved in this build. The upstream project license text was
retrieved from <https://raw.githubusercontent.com/3dmol/3Dmol.js/master/LICENSE>
on 2026-11-03; no viewer code was replaced during release preparation.

`molstar.LICENSE.txt` contains the existing Mol* MIT license and notices for
incorporated React, Immutable.js, normalize.css and attributed shader code.
The normalize.css license was retrieved from its recorded version:
<https://raw.githubusercontent.com/necolas/normalize.css/3.0.3/LICENSE.md>.
These report assets and their notices are part of the scientific workflow.

Exact copies of the two companion notices are archived under
[`licenses/viewers/`](../../../../licenses/viewers/) for repository-level review.
The copies here remain runtime resources for source and installed reports.

The report generator copies each selected viewer and its accompanying notice
into the report's asset snapshot. `assets_manifest.json` records the notice's
size and SHA-256 and associates it with the viewer. Custom viewer sources should
include the same companion filenames. A bundled notice is reused for a custom
asset only when the asset bytes match; otherwise a missing custom notice is
reported in the manifest. Mixed sources with different notices retain both.
Keep these notice files when publishing a generated report. See the repository's
[third-party notices](../../../../licenses/THIRD_PARTY_NOTICES.md) for the other bundled
components and external-dependency scope.

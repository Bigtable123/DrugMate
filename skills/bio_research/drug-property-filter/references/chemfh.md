# ChemFH backend

## Placement and runtime

ChemFH runs in the final `property_filter` step alongside ADMET-AI. It is off
in the pre-docking RDKit prefilter and on by default in the final filter. The
source project is `${SBDD_CHEMFH_SOURCE_DIR}`; the isolated runtime prefix is
`${SBDD_ENV_CHEMFH_PREFIX}`. The Drugmate process communicates with one CPU batch
worker through strict JSON and never imports the Django/Chemprop stack.

The runtime pins and verifies:

- source commit `f8c046fac9a1d331979025a104382d3d2b2614ca`;
- native Chemprop checkpoint SHA-256
  `9f5c82630bed2d787eaf16d99eeac500763748d54474818757207659433a8b6b`;
- combined 17-rule asset SHA-256
  `f67d0cf228af15cd3ccf200034f0363a1495aaba0f25388e0393379cd1e55519`.

The worker checks the assets again after prediction so a native rule loader
cannot silently mutate the scientific inputs during a run.

## Model endpoints

The seven native endpoints are other assay interference, blue fluorescence,
FLuc inhibitors, promiscuous compounds, colloidal aggregators, reactive
compounds, and green fluorescence. A model call is positive only when its
probability is strictly greater than `0.5`. Each endpoint retains the native
probability, uncertainty, its endpoint-specific uncertainty threshold, and the
resulting confidence annotation.

The scoring policy applies the exact table bounds to each raw probability as a
level-1 ALERT mark: `<=0.30` is `pass`, `0.30 < p <= 0.90` is `flag`, and
`p > 0.90` is `fail`. Endpoint weights are not uniform: reactive compounds weigh
`1`, promiscuous compounds and colloidal aggregators `0.5` each, and the four
fluorescence/interference endpoints `0`. A `fail` here is one mark on one metric
and never rejects the molecule by itself — only the ALERT category fold does
. The provider's native `>0.5` positive call remains a
conservative `review_only` floor that cannot weaken the `>0.90` mark.
Uncertainty is reported and identifies high-confidence positives, but it does
not reverse either raw probability decision. Missing endpoint/uncertainty data,
incomplete batch coverage, backend failure, duplicate IDs, or a
candidate-ID/SMILES mismatch also caps automatic nomination and is never
converted to favorable evidence.

## Native structural rules and overlap

All 17 rule families are retained with matched SMARTS and atom indices:
Aggregators, Fluc, Blue/Green fluorescence, Reactive,
Other-assay-interference, Promiscuous, ALARM_NMR, BMS, Chelator, GST/His FH,
Luciferase inhibitor, NTD, PAINS, Potential Electrophilic, and Lilly.

The worker calls the native `FrequentHitterFilter.screening` directly and
retains its SMARTS order and first atom match per SMARTS, as in ChemFH's native
tables. SVG rendering is disabled: the pipeline does not consume ChemFH's
depiction tables, and 2D layout failures must not interrupt rule screening.
This adapter leaves the input SMILES, model, rule assets, and asset checks
unchanged; provenance records `rule_depiction_enabled: false`.

Individual native rule hits are provider evidence only: no rule creates its own
gate. Ten of them — ALARM_NMR, BMS, Chelator, GST FH, His FH, Luciferase
inhibitor, NTD, PAINS, Potential Electrophilic, and Lilly — are aggregated into
the single `chemfh.substructure_alert_count` metric (`soft_upper=2`,
`hard_upper=10`, weight `2` inside ALERT), counting hit rules rather than SMARTS
multiplicity. The other seven rule families stay pure audit rows. In particular,
the ChemFH PAINS result adds no second PAINS gate: it is marked as overlapping
the existing RDKit `pains` tool, so its duplicate decision contribution is
suppressed and the fast RDKit catalogs remain the sole per-catalog
structural-alert layer.

## Nomination reducer

ADMET-AI and ChemFH attach evidence first. One central reducer then grants
`auto_candidate` only when the RDKit verdict passes, ADMET is complete, ChemFH
is complete, and no ChemFH model endpoint is positive. This prevents one
backend from overwriting a review cap set by another. AiZynthFinder remains a
separate optional downstream assessment and does not become a dependency of
the `chemfh` tool. The shared raw metric reducer then marks the ChemFH
probabilities at `0.30/0.90` and folds the ALERT category; only a fold rejection
keeps a candidate out of retrosynthesis, so an isolated `>0.90` endpoint caps the
candidate at `modify` rather than dropping it.

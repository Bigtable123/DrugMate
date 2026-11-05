# Retrospective benchmark labels

The existing benchmark reads molecule names, using `ACT_` for known actives and
`DEC_` for known negative/decoy controls. Numeric IDs alone are unlabeled.
This convention needs no additional pipeline config fields.

Before input normalization:

1. Locate an authoritative label source for the same assay/target and exact
   compound identities. Do not infer labels from docking scores, file order,
   numeric IDs, the target name, or a parent molecule's activity.
2. Create a derived input under the new round's `inputs/`; keep the original
   library untouched. Use `ACT_4256661` only if 4256661 is a verified active,
   or `DEC_4256661` only if it is a verified negative/decoy control. Unknown
   compounds retain their original names. Existing verified prefixes stay intact.
3. Preserve an ID map with `original_id`, `screening_id`, `label`, and the
   label source/assay reference. Record whether negatives are experimentally
   inactive or constructed decoys; the prefix is an implementation convention,
   not a claim that these evidence types are identical. For SDF, the title and
   any explicit compound ID fields must identify the same screening record;
   preserve `original_id` separately. Reject conflicting labels or ambiguous
   many-to-one ID joins instead of guessing.
4. Normalize this derived input and use the resulting library consistently in
   all later stages. Never rename only the scores or an already committed SDF.
5. Check both class counts in the actual sampled/retained cohort. Only then
   enable `coordination.screening_benchmark`; retain scientifically justified
   minimum counts. Do not lower them solely to obtain a passing gate.

For unlabeled candidate screening or a smoke test, set the benchmark tool to
`enabled=false`, keep the native-redocking and pose-quality checks, and report
AUROC/EF/PR-AUC as not evaluated. No fabricated label or metric is needed.

`DEC_*` records are excluded from candidate nomination. Use this prefix only
for compounds deliberately assigned to the benchmark control pool. A missing
score for a labeled compound, including a quarantined docking failure, prevents
validation on a selectively retained subset. Explicitly requested benchmarking
without both classes produces a protocol `modify`, not an implicit pass.

This accounting covers the docking input cohort. If normalization or sampling
removed compounds earlier, disclose those exclusions and the retained class
counts; do not describe the metrics as coverage of the original submitted library.

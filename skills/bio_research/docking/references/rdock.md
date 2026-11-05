# rDock adapter

## Runtime and cavity definition

The canonical primary-engine ID is `rdock`, registered at
`${SBDD_ENV_RDOCK_PREFIX}`. The adapter invokes `obabel`, `rbcavity`, and
`rbdock` from that prefix. Every invocation sets:

- `RBT_ROOT` points to `${SBDD_ENV_RDOCK_PREFIX}/share/rDock`
- `RBT_HOME=<isolated native run directory>`

`search_space_source=reference_ligand` uses `RbtLigandSiteMapper` and the
explicit reference-ligand artifact. `search_space_source=docking_box` uses the
installed `RbtSphereSiteMapper`, centered on the box
with a half-diagonal radius that contains all eight corners. The adapter stages
the receptor and reference ligand under canonical relative names so equivalent
runs do not acquire different protocol identities from temporary paths.
`rbcavity` writes the declared active site before `rbdock` runs the installed
`data/scripts/dock.prm` protocol.

## Scores and identity

The primary metric is `rdock.score`, lower is better. Its unit is explicitly
`rDock score`; it is not labeled kcal/mol. `SCORE.INTER`, `SCORE.INTRA`, and
`SCORE.RESTR` are preserved as diagnostic components. Native SD chemistry is
used only for coordinates and scores; poses are placed back on the sanitized
source topology after strict atom-count, element-order, molecule-coverage, and
finite-score checks.

`rdock.score` is a level-1 `AFF` mark with no universal absolute cutoff:
retain finite values for same-protocol ranking with suggested signed weight
`-5`; missing/non-finite values are `not_applicable` and an availability
blocker. The category fold, not a single finite score, determines the
molecule-level `keep`/`modify`/`reject` result.

The default request is `num_modes=20`, capped to 20/5/1 by cohort size; users may request fewer in the
experiment config.

## Boundaries

- The adapter is non-covalent and rigid-receptor only.
- rDock help commands may return nonzero even when installed; availability is
  preflighted from the complete executable set rather than `rbdock -h` alone.
  Installation acceptance should additionally run the documented minimal
  cavity+docking smoke.
- Native rDock values form their own score family and cannot be pooled with
  Vina, GNINA, AutoDock4, ADFR, DiffDock, or DOCK 6 raw values.

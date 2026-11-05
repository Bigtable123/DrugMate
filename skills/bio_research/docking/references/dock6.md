# DOCK 6 adapter

## Receptor artifact boundary

DOCK 6 classic grid docking requires a receptor-site sphere file and a matched
energy/bump grid pair. The current structure-preparation step publishes a PDB,
not a charged and atom-typed DOCK receptor Mol2. Although the local DOCK 6
installation contains `sphgen`, `sphere_selector`, `showbox`, and `grid`, those
programs do not replace the missing trusted receptor-preparation/molecular-
surface stage. The adapter therefore never fabricates a receptor grid.

An explicitly selected DOCK 6 plan must set
`docking.params.search_space_source=precomputed_grid` and inject all three
artifacts before execution:

- `dock6_receptor_site` -> `dock6_receptor_site_file` (`.sph`);
- `dock6_grid_nrg` -> `dock6_grid_nrg_file` (`.nrg`);
- `dock6_grid_bmp` -> `dock6_grid_bmp_file` (`.bmp`).

The `.nrg` and `.bmp` paths must be non-empty regular files with one shared
prefix. The validator blocks the run before worker execution if the source is
different or any artifact is absent; the adapter never fabricates or selects a
grid.

The default request is `num_modes=20`, capped to 20/5/1 by cohort size; users may request fewer in the
experiment config. Report actual output counts against the requested budget.

## Ligand preparation and execution

The adapter splits the canonical Step 2 library SDF one record at a time and
uses `${SBDD_DOCK6_SOURCE_DIR}/bin/antechamber` to write charged, SYBYL-typed Mol2
records. `charge_method=gas` is the fast default; `bcc` is available as an
explicit, slower choice. Every ligand must prepare successfully before the
combined Mol2 is passed to the DOCK 6.13 flexible anchor-and-grow/grid-scoring
protocol.

Native runs live under `dock6_native/run_*`, so stale scored files cannot enter
a new result. The adapter requires at least one finite, parseable `Grid_Score`
pose for every input molecule and rejects unknown source names, invalid Mol2,
non-finite scores, excess ranks, or incomplete molecule coverage.

## Score semantics

The canonical score fields are:

- `dock6_grid_score`: native DOCK 6 grid score, lower is better;
- `dock6_grid_vdw_energy`: native van der Waals component, diagnostic;
- `dock6_grid_es_energy`: native electrostatic component, diagnostic.

DOCK 6 grid-score units and calibration are protocol-specific. Raw values are
not interchangeable with Vina-family kcal/mol energies, DiffDock confidence,
or RTMScore outputs. Compare poses only within one protocol ID whose sphere,
grid, definition tables, ligand preparation, and search parameters match.

`dock6_grid_score` has no universal absolute cutoff in the executable policy and carries
suggested signed weight `0`; a finite value is retained as level-1 audit
evidence, while missing/non-finite output is `not_applicable` and an
availability blocker. It must not independently produce a molecule-level
`keep`, `modify`, or `reject` decision.

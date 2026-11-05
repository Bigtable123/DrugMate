# ADFR / AGFR composite adapter

## Runtime and components

The canonical primary-engine ID is `adfr`. ADFRsuite is the distribution name,
not a second tool ID. Its prefix is `${SBDD_ENV_ADFRSUITE_PREFIX}`, and the one
adapter owns:

1. `prepare_receptor` and `prepare_ligand` for native PDBQT preparation;
2. `agfr` for the maps-based `.trg` target;
3. `adfr` for genetic-algorithm docking.

The installed suite carries its own compatibility Python and AutoGrid build.
The adapter therefore invokes the suite wrappers directly and never imports
their Python 2 packages into the Drugmate process.

## Contract

The adapter consumes prepared receptor PDB, a 3D multi-record ligand SDF, and
the search-space artifact named by the explicit
`docking.params.search_space_source`: `reference_ligand` or `docking_box`.
AGFR receives that exact definition and emits one auditable target. Each ligand
is prepared and docked separately, then native PDBQT poses are mapped back onto
canonical source topology with complete-coverage and finite-score checks. The
adapter does not infer the source from which artifact happens to exist.

The primary metric is `adfr.affinity` in kcal/mol, lower is better. ADFR
components remain provider-namespaced context and are not interchangeable with
AutoDock4 or Vina energies.

## SBDD score marks

The registered `AutoDockFR` row uses these level-1 `AFF` marks and
signed weight:

| Metric | Pass mark | Flag mark | Fail mark | Suggested signed weight |
|---|---|---|---|---:|
| `adfr.affinity` (AutoDockFR flexible profile) | `<= -6` kcal/mol | `(-6, 0]` | `> 0` | -4 |

The category fold applies V1/V2/V3 after all planned metrics are marked; an
individual `fail` does not itself reject a molecule. Missing/non-finite values
are `not_applicable` and an availability blocker.

The default request is `num_modes=20`, capped to 20/5/1 by cohort size; users may request fewer in the
experiment config.

## Boundaries

- The registered adapter is rigid-receptor, non-covalent docking only.
- The installed CLI exposes flexible-residue and covalent options, but the
  pipeline does not yet have the required residue/atom/bond-geometry contract;
  `is_covalent=true` still hard-fails before docking.
- The canonical numeric row above describes a flexible-receptor profile. The
  current metric registry aliases the rigid adapter's `adfr.affinity` to that
  row and applies its `-6/0` bounds and `-4` signed weight. This preserves the
  source numbers but is a known applicability gap; a rigid result must not be
  described as validated flexible-receptor evidence.
- The installed `adfr --version` path is unreliable under its legacy argparse
  stack. Availability uses functional wrapper/help/import probes, while run
  provenance records the suite identity conservatively.

# AutoDock4 composite adapter

## Runtime and components

The canonical primary-engine ID is `autodock4`. Its registered prefix is
`${SBDD_ENV_AUTODOCK4_PREFIX}`. The adapter owns three installed capabilities:

1. AutoDockTools prepares receptor/ligand PDBQT plus GPF/DPF inputs.
2. AutoGrid4 creates the atom-type affinity and electrostatic maps.
3. AutoDock4 performs the Lamarckian genetic-algorithm search.

`autogrid4` and `autodocktools` are not separate `enabled` tools. Registering
them separately would create duplicate engine choices, time budgets, and score
evidence for one inseparable protocol. This adapter uses its native preparation;
Step 0 Meeko remains specific to the `autodock_vina` adapter.

## Contract

The adapter consumes prepared receptor PDB, a 3D multi-record ligand SDF, and
the artifact named by the explicit `docking.params.search_space_source`:
`reference_ligand` or `docking_box`. It stages every ligand separately, builds
maps covering the union of required AutoDock atom types, runs one DPF/DLG per
ligand, and exports the requested finite scored poses to the canonical
`docked_poses.sdf`. It never chooses between those sources itself.

The primary score is `autodock4.binding_energy` in kcal/mol, lower is better.
DLG/native component values remain provider-namespaced audit evidence. They
must not be copied into GNINA/Vina fields or compared across engines without a
calibrated protocol.

## SBDD score marks

The following are level-1 per-metric marks for the canonical `AFF` category;
the shared scorer folds them with the other category marks before assigning a
molecule-level classification.

| Metric | Pass mark | Flag mark | Fail mark | Suggested signed weight |
|---|---|---|---|---:|
| `autodock4.binding_energy` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -4 |

These fallback bounds apply only within one target/pocket cohort with the same
AutoDock4 version, maps, scoring parameters, receptor/ligand preparation, and
search settings. Missing or non-finite values are `not_applicable` and an
availability blocker, not a synthetic zero.

The default request is `num_modes=20`, capped to 20/5/1 by cohort size; users may request fewer in the
experiment config. AutoDock4Zn and flexible/covalent variants have separate
profile rows in the executable policy and are not implied by this rigid adapter.

## Boundaries

- Only ordinary non-covalent rigid-receptor docking is registered.
- AutoDock4Zn, flexible-receptor, covalent, smina, and QuickVina variants are
  not implied by this adapter.
- Missing maps, incomplete ligand coverage, atom-order mismatch, or a
  non-finite primary energy fails the whole cohort.
- Grid/map inputs and their content identities are retained in protocol
  provenance; installation paths do not define score comparability.

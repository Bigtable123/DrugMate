# PoseBusters Reference

## Overview

PoseBusters validates the geometric plausibility of docked poses by checking bond geometry, internal clashes, protein-ligand clashes, and stereochemistry.

- Paper: Buttenschoen et al., Chemical Science 2024
- GitHub: https://github.com/maabuu/posebusters
- Install: `pip install posebusters`

## Configuration Modes

| Config | Use Case |
|--------|----------|
| `"dock"` | Validate docked poses against a receptor |
| `"redock"` | Validate redocked poses (stricter RMSD checks) |
| `"mol"` | Validate molecule geometry only (no receptor) |

## Python API

```python
from posebusters import PoseBusters

pb = PoseBusters(config="dock")
df = pb.bust("docked.sdf", None, "receptor.pdb")
# df: DataFrame with bool columns per check, one row per molecule
```

## Checks Performed (dock config)

| Check | Description |
|-------|-------------|
| `mol_pred_loaded` | Molecule can be parsed |
| `sanitization` | RDKit sanitization passes |
| `all_atoms_connected` | No disconnected fragments |
| `bond_lengths` | All bonds within expected length range |
| `bond_angles` | All angles within expected range |
| `internal_steric_clash` | No intra-molecular atom overlaps |
| `aromatic_ring_flatness` | Aromatic rings are planar |
| `double_bond_flatness` | Double bonds have correct planarity |
| `internal_energy` | Strain energy within threshold |
| `protein_ligand_maximum_distance` | Ligand is near protein |
| `minimum_distance_to_protein` | No impossibly close contacts |
| `volume_overlap_with_protein` | Ligand doesn't penetrate protein |

## Output Format

DataFrame columns are boolean check results (`True` = check passed, `False` =
check failed). `df.all(axis=1)` may be retained only as a raw diagnostic; it is
not the final pose gate. Map each registered applicable check to its level-1
`pass`/`flag`/`fail` mark, then fold all marks within their parent metric
categories under the active V1/V2/V3 mode. Unregistered checks remain
`not_applicable` until a reviewed rule supplies their threshold and category.

## Known Issues

| Issue | Workaround |
|-------|-----------|
| Requires RDKit-parseable SDF | Ensure GNINA output has proper SDF format |
| Slow on large libraries (>500 mols) | Process in batches |
| Hydrogens can affect clash detection | GNINA output includes H by default |

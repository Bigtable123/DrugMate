# PoseCheck Reference

## Overview

PoseCheck evaluates docked poses by computing ligand strain energy and protein-ligand steric clashes.

- GitHub: https://github.com/cch1999/posecheck
- Install: `pip install posecheck`

## Python API

The local wrapper reuses the final explicit-H receptor produced by docking
preparation. Reduce is the usual backend; retain PDB2PQR's receptor when that is
the selected protocol. The receptor loading path is:

```python
from rdkit import Chem
from posecheck import PoseCheck
from posecheck.utils.strain import calculate_strain_energy

pc = PoseCheck()
protein = Chem.MolFromPDBFile(
    "docking_prep/prepared_receptor.pdb",
    sanitize=False, removeHs=False, proximityBonding=False,
)
if protein is None or protein.GetNumAtoms() == 0 or protein.GetNumConformers() == 0:
    raise ValueError("Cannot read prepared receptor atoms and coordinates")
if not any(atom.GetAtomicNum() == 1 for atom in protein.GetAtoms()):
    raise ValueError("Prepare the receptor with explicit hydrogens before PoseCheck")
pc.protein = protein

mol = Chem.MolFromMolFile("ligand.sdf")
pc.load_ligands_from_mols([mol])

clashes = pc.calculate_clashes()        # list[int]: clash atom count per mol
strain = [calculate_strain_energy(mol, num_confs=10) for mol in pc.ligands]
```

Use the evaluator's existing `run_posecheck()` wrapper in production. It validates
the receptor and preserves the native ligand preparation and strain calculation.
RDKit reads the existing coordinates without sanitizing or inferring proximity
bonds; collision distances and atom radii do not require receptor bond inference.
It preserves explicit H and PDB residue identifiers, including negative numbers
such as `ACE A -1`. This avoids both Hydride's nonconverging hydrogen relaxation
and the upstream loader's ProLIF conversion failure on negative residue numbers.
Do not call `pc.load_protein_from_pdb()`, which regenerates receptor hydrogens.

## Metrics

| Metric | Unit | Good | Warn | Bad |
|--------|------|------|------|-----|
| Strain energy | kcal/mol | <=8 | (8, 60] | >60 |
| Steric clashes (heavy atoms) | count | 0 | 1-2 | >2 |

These bands produce level-1 `pass`/`flag`/`fail` marks. A `fail` in either
PoseCheck metric does not by itself reject a pose; rejection occurs only after
all available metrics in its parent category are folded under the active
V1/V2/V3 screening mode.

## Strain Energy Calculation

The installed Drugmate PoseCheck runtime uses UFF: it first performs a locally
constrained relaxation (maximum displacement 0.1 Å), then compares that energy
with the lowest energy from independently embedded/relaxed conformers and an
unconstrained relaxation of the starting pose. High strain indicates an
energetically unfavorable conformation. This is separate from the optional
`rdkit_pose` tool's MMFF94s calculation.

The upstream `PoseCheck.calculate_strain_energy()` method samples 50 conformers
per pose and does not expose that count. The local wrapper calls the same
underlying strain function with **10 conformers by default**. Override it with
`tools.posecheck.params.num_conformers` (a positive integer); standalone evaluator
input also accepts `posecheck_num_conformers`. Explicit native-ligand calibration
uses the same count as candidate poses. Results record `num_conformers`, and
changing the count invalidates evaluation reuse. Set 50 to restore the previous
sampling budget. Reducing the sample count changes the estimated reference
minimum and can change strain values and downstream verdicts; it is a sampling
cost/coverage choice, not an equivalent-score optimization.
Reusing the prepared receptor does not change this conformer sampling budget;
`num_conformers` controls ligand strain sampling, not receptor hydrogen optimization.

## Steric Clash Detection

Counts atom pairs between ligand and protein where the distance is less than the sum of van der Waals radii minus a tolerance (typically 0.4 Å).

The local structure-evaluator wrapper uses PoseCheck's all-atom count as a
diagnostic field (`steric_clash_count_all_atom`) and reports `steric_clash_count`
as the heavy-atom-only count used for verdicts and native calibration. Both use
the same final receptor coordinates as docking. The all-atom count includes that
receptor's prepared hydrogens and can differ from older runs that regenerated H
with Hydride; the heavy-atom count is unchanged when heavy-atom coordinates are
the same.

## Known Issues

| Issue | Workaround |
|-------|-----------|
| Requires RDKit for mol loading | Always have RDKit installed |
| Receptor lacks explicit H or usable atom coordinates | Run final docking preparation and supply its `prepared_receptor`; the wrapper reports failed evidence instead of regenerating H |
| MMFF parameterization can fail for unusual atoms | Returns None for strain |
| Memory usage scales with mol count | Process one at a time for large sets |

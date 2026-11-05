# ProLIF Reference

## What it does

ProLIF (Protein-Ligand Interaction Fingerprints) computes interaction fingerprints (IFP) that encode which protein residues interact with a ligand, and through which interaction types. ProLIF 2.x also stores geometric metadata for detected interactions, including distances, angles, and atom indices.

## Installation

```bash
pip install prolif mdanalysis rdkit
```

Requires: RDKit (ligand chemistry and SDF loading) and MDAnalysis (protein structure loading).

## Basic Usage

### From separate protein PDB and ligand SDF

MDAnalysis does not parse SDF as a topology format. Load SDF ligands with RDKit, then convert to ProLIF.

```python
import prolif
import MDAnalysis as mda
from rdkit import Chem

# ligand SDF via RDKit
suppl = Chem.SDMolSupplier("ligand.sdf", removeHs=False)
rdkit_mol = next(m for m in suppl if m is not None)
ligand = prolif.Molecule.from_rdkit(rdkit_mol, resname="LIG")

# protein PDB via MDAnalysis
u = mda.Universe("protein_prepared.pdb")
protein_sel = u.select_atoms("protein")
protein = prolif.Molecule.from_mda(protein_sel, force=True)

fp = prolif.Fingerprint()
fp.run_from_iterable([ligand], protein, progress=False)
df = fp.to_dataframe()
```

### Recommended: select nearby protein residues

For large proteins, convert only residues near the ligand to avoid slow or unstable RDKit bond inference.

```python
lig_coords = ligand.GetConformer().GetPositions()
centroid = lig_coords.mean(axis=0)
nearby = u.select_atoms(
    f"byres (protein and point {centroid[0]} {centroid[1]} {centroid[2]} 6.5)"
)
protein = prolif.Molecule.from_mda(nearby, force=True)
```

### From a complex PDB

```python
import prolif
import MDAnalysis as mda

u = mda.Universe("complex.pdb")
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname LIG")

fp = prolif.Fingerprint()
fp.run(u.trajectory, ligand, protein, progress=False)
df = fp.to_dataframe()
```

## Interaction Types

| ProLIF Name | Description | Parameters |
|-------------|-------------|------------|
| HBDonor | H-bond with residue/ligand as donor | dist < 3.5Å, angle > 130° |
| HBAcceptor | H-bond with residue/ligand as acceptor | dist < 3.5Å, angle > 130° |
| Hydrophobic | Carbon-carbon hydrophobic contact | dist < 4.5Å |
| PiStacking | Face-to-face aromatic stacking | dist < 5.5Å, angle < 35° |
| EdgeToFace | T-shaped aromatic stacking | dist < 6.5Å, angle 50°-90° |
| Cationic | Positive charge interaction | dist < 4.5Å |
| Anionic | Negative charge interaction | dist < 4.5Å |
| CationPi | Cation-pi interaction | dist < 4.5Å |
| PiCation | Pi-cation interaction | dist < 4.5Å |
| MetalAcceptor | Metal coordination | dist < 2.8Å |
| VdWContact | Van der Waals contact | dist < sum(vdW) + 0.5Å |

## Specifying Interaction Types

```python
fp = prolif.Fingerprint(["HBDonor", "HBAcceptor", "Hydrophobic", "VdWContact"])
```

## Output DataFrame Structure

ProLIF 2.x commonly uses MultiIndex columns:

```
(ligand_residue, protein_residue, interaction_type)
("LIG1", "LEU83.A", "HBDonor")
```

Rows are frames; a single structure has one row. Values are booleans.

## Extracting Summary

```python
df = fp.to_dataframe()
for col in df.columns:
    if df[col].any():
        lig_res, prot_res, interaction_type = col
        print(prot_res, interaction_type)
```

## Extracting Geometry Metadata

ProLIF stores interaction metadata in `fp.ifp`:

```python
for frame, interactions in fp.ifp.items():
    for (lig_res, prot_res), by_type in interactions.items():
        for interaction_type, entries in by_type.items():
            for entry in entries:
                print(prot_res, interaction_type, entry)
                # entry may include distance, angle fields, indices, parent_indices
```

Use this metadata when building strict interaction constraints or an Interaction Blueprint.

## Pitfalls

- **SDF loading**: Do not use `mda.Universe("ligand.sdf")`; load SDF with RDKit first.
- **No hydrogens**: The wrapper requires explicit H/D atoms on selected amino-acid
  residues and reports missing protein hydrogens as a profiler failure. Water or
  ligand H atoms cannot satisfy the check. In pipeline runs, explicitly select
  `structure_prep.reduce` with `build=true` and a reason to generate the separate
  `provisional_receptor` analysis copy; retain the noH baseline. This is provisional
  interaction evidence requiring review against the final docking receptor.
- **Large proteins**: Convert only nearby residues to ProLIF molecule objects to avoid RDKit bond-inference instability.
- **Empty fingerprint**: Usually means ligand and protein are in different coordinate systems, the ligand selection is wrong, or the cutoff is too strict.
- **Column shape**: ProLIF 2.x columns are usually `(ligand_residue, protein_residue, interaction_type)`.

## References

- Bouysset & Fiorucci, *J. Cheminformatics* 13:72 (2021) — ProLIF
- Desaphy et al., *J. Chem. Inf. Model.* 53:623 (2013) — IFP concept

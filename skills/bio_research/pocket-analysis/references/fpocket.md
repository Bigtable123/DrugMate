# fpocket Reference

## What it does

fpocket detects protein binding pockets using Voronoi tessellation and alpha sphere clustering. It identifies cavities on the protein surface that could accommodate small-molecule ligands.

## Installation

```bash
# conda (recommended)
conda install -c conda-forge fpocket

# from source
git clone https://github.com/Discngine/fpocket.git
cd fpocket && make && sudo make install
```

## Key Parameters

| Parameter | Flag | Default | Recommended | Description |
|-----------|------|---------|-------------|-------------|
| Min alpha sphere radius | `-m` | 3.4 | 3.0 | Minimum radius to keep (smaller = more internal) |
| Max alpha sphere radius | `-M` | 6.2 | 6.0 | Maximum radius to keep (larger = surface) |
| Min spheres per pocket | `-i` | 30 | 30 | Minimum cluster size for a valid pocket |
| Clustering distance | `-D` | 1.73 | 1.73 | Distance for sphere clustering |
| Druggability | `-d` | on | on | Compute druggability score |

## Usage

```bash
fpocket -f protein.pdb -m 3.0 -M 6.0
```

## Output Structure

```
protein_out/
├── protein_info.txt          # Summary: all pockets with scores
├── protein_out.pdb           # PDB with pocket annotations (HETATM STP)
├── protein_pockets.pqr       # All pocket alpha spheres
└── pockets/
    ├── pocket1_atm.pdb       # Atoms forming pocket 1
    ├── pocket1_vert.pqr      # Alpha sphere centers for pocket 1
    ├── pocket2_atm.pdb
    └── ...
```

## Parsing _info.txt

Format:
```
Pocket 1 :
        Score :                 0.4532
        Druggability Score :    0.872
        Number of Alpha Spheres :       98
        Total SASA :            234.5
        Polar SASA :            45.2
        Apolar SASA :           189.3
        Volume (MC) :           452.3
        ...
```

## Druggability Score Interpretation

| Score | Interpretation |
|-------|---------------|
| > 0.8 | Highly druggable |
| 0.5 - 0.8 | Moderately druggable |
| < 0.5 | Difficult target |

## Pitfalls

- **No pockets found**: Check that PDB has proper ATOM records (not just CA trace). Also try relaxing `-m` to 2.8.
- **Too many pockets**: Increase `-i` to 40-50 to filter small cavities.
- **Volume underestimated**: fpocket volume (Monte Carlo) is approximate; for precise volume use POVME or mdpocket.
- **Membrane proteins**: fpocket may detect the membrane-facing cavity; filter by residue solvent accessibility.

## References

- Le Guilloux et al., *BMC Bioinformatics* 10:168 (2009) — fpocket algorithm
- Schmidtke et al., *J. Med. Chem.* 53:5858 (2010) — fpocket druggability

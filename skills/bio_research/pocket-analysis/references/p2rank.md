# P2Rank Reference

## What it does

P2Rank predicts ligand-binding sites on protein surfaces using a machine learning model (random forest) trained on surface point features (chemical properties, geometry, conservation). It outperforms geometric methods (fpocket, SiteMap) on most benchmarks.

## Installation

```bash
# Download release (Java required)
wget https://github.com/rdk/p2rank/releases/download/2.4/p2rank_2.4.tar.gz
tar xzf p2rank_2.4.tar.gz
export PATH=$PATH:$(pwd)/p2rank_2.4

# Verify
prank --version
```

Requires Java 11+.

## Usage

```bash
# Single protein
prank predict -f protein.pdb -o output_dir/

# Multiple proteins (batch)
prank predict -f proteins.ds -o output_dir/
```

## Output Structure

```
output_dir/
├── protein.pdb_predictions.csv     # Ranked pocket list
├── protein.pdb_residues.csv        # Per-residue pocket assignment
└── visualizations/
    └── protein.pdb.pml             # PyMOL visualization script
```

## Predictions CSV Format

```csv
name,rank,score,prob,sas_points,surf_atoms,center_x,center_y,center_z,residue_ids
pocket1,1,24.31,0.924,142,89,12.5,34.2,8.7,"A_LEU83 A_GLU81 A_PHE80 A_ILE10 A_VAL18"
pocket2,2,8.45,0.612,67,42,28.1,15.3,22.4,"A_ASP210 A_ARG213"
```

| Column | Description |
|--------|-------------|
| name | Pocket identifier |
| rank | Ranking by score (1 = best) |
| score | P2Rank confidence score (higher = better, typically 5-40) |
| prob | Probability of being a true binding site (0-1) |
| sas_points | Number of solvent-accessible surface points |
| surf_atoms | Number of surface atoms in pocket |
| center_x/y/z | Geometric center coordinates |
| residue_ids | Space-separated residue IDs (format: Chain_ResnameResnum) |

## Score Interpretation

| Score | Probability | Confidence |
|-------|-------------|------------|
| > 20 | > 0.85 | High confidence binding site |
| 10-20 | 0.6-0.85 | Likely binding site |
| 5-10 | 0.3-0.6 | Possible, needs validation |
| < 5 | < 0.3 | Unlikely |

## Residues CSV Format

Per-residue data:
```csv
chain,residue_label,residue_name,residue_number,pocket,score
A,A_LEU83,LEU,83,pocket1,0.95
A,A_GLU81,GLU,81,pocket1,0.92
```

## Pitfalls

- **Wrong top pocket on multi-domain proteins**: P2Rank picks the most "drugable" surface — if targeting an allosteric site, manually compare pocket centers with target residues.
- **Low score for shallow sites**: P2Rank is trained on deep pockets; flat binding surfaces (PPI interfaces) score lower.
- **Java memory for large complexes**: Use `prank predict -f file.pdb -Xmx4g` for proteins > 1000 residues.

## Performance vs fpocket

On COACH420 benchmark:
- P2Rank: 73% success rate (DCA < 4Å)
- fpocket: 58% success rate
- SiteMap: 65% success rate

P2Rank is generally preferred; fpocket provides complementary geometric information.

## References

- Krivak & Hoksza, *J. Cheminformatics* 10:39 (2018) — P2Rank algorithm
- Krivak & Hoksza, *J. Chem. Inf. Model.* 55:2423 (2015) — Ligand binding site prediction

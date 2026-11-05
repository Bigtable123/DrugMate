# DUD-E (Database of Useful Decoys: Enhanced) Reference

## Overview

DUD-E provides property-matched decoys for benchmarking virtual screening methods.
Each target has ~50 active ligands paired with ~39 decoys per active (property-matched
but topologically dissimilar).

**Website:** https://dude.docking.org/

---

## URL Patterns

### Decoys (SMILES format, smaller download)

```
https://dude.docking.org/targets/{target_name}/decoys_final.ism
```

### Decoys (SDF format, 3D coordinates)

```
https://dude.docking.org/targets/{target_name}/decoys_final.sdf.gz
```

### Actives (SMILES format)

```
https://dude.docking.org/targets/{target_name}/actives_final.ism
```

### Actives (SDF format)

```
https://dude.docking.org/targets/{target_name}/actives_final.sdf.gz
```

**Note:** Target names are lowercase (e.g., "cdk2", "egfr", "aa2ar").

---

## .ism File Format

Plain text, one compound per line:
```
SMILES_STRING ZINC_ID
```

Example:
```
CC(=O)Nc1ccc(O)cc1 ZINC00000001
COc1ccc(NC(=O)c2ccccc2)cc1 ZINC00000002
```

---

## Available Targets (subset relevant to SBDD)

| DUD-E Name | Full Name | # Actives | # Decoys |
|------------|-----------|-----------|----------|
| cdk2 | Cyclin-dependent kinase 2 | 72 | 2810 |
| egfr | Epidermal growth factor receptor | 542 | 35050 |
| braf | B-Raf kinase | 152 | 9950 |
| abl1 | Tyrosine kinase ABL1 | 182 | 10750 |
| jak2 | Janus kinase 2 | 107 | 6350 |
| src | Proto-oncogene tyrosine kinase Src | 524 | 34500 |
| vegfr2 | VEGF receptor 2 | 409 | 24950 |
| akt1 | RAC-alpha serine/threonine kinase | 293 | 16450 |
| erk2 | MAP kinase ERK2 | 76 | 4550 |
| p38 | MAP kinase p38 alpha | 454 | 25350 |

Full list: https://dude.docking.org/targets

---

## Property Matching Methodology

DUD-E decoys are selected to match actives on:
- Molecular weight (MW)
- Calculated LogP (cLogP)
- Number of rotatable bonds
- Number of hydrogen bond donors (HBD)
- Number of hydrogen bond acceptors (HBA)
- Net molecular charge

While being topologically dissimilar (Tanimoto < 0.35 using ECFP4 fingerprints).

---

## Usage Notes

1. **Not all targets available**: Only ~102 targets in DUD-E. If your target isn't there, use RDKit property-matched generation as fallback.

2. **Decoy count**: Typically ~39 decoys per active. For CDK2: 72 actives × 39 ≈ 2810 decoys available. We sample randomly if we need fewer.

3. **3D coordinates**: The .sdf.gz files contain pre-generated 3D conformers. The .ism files are smaller but require conformer generation.

4. **Bias warning**: DUD-E decoys have known biases (analogue bias, artificial enrichment). For rigorous benchmarking, consider also using DEKOIS 2.0 or MUV datasets.

5. **Network dependency**: DUD-E downloads require internet access. The fallback RDKit decoy generator works offline.

---

## Fallback: RDKit Property-Matched Generation

When DUD-E is unavailable, our `decoy_gen_utils.py` generates decoys by:

1. BRICS fragmentation of active compounds
2. Random recombination of fragments
3. Property filtering (MW, LogP, HBD, HBA, RotBonds within tolerance)
4. Tanimoto filtering (ECFP4 similarity < 0.35 to any active)

This follows the same principles as DUD-E but generates de novo compounds.

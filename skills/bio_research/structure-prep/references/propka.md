# PROPKA Reference

## What it does

PROPKA predicts pKa values for ionizable residues in proteins based on the 3D structure. This tells you which residues are protonated/deprotonated at a given pH — critical for correct molecular docking and electrostatic calculations.

## Installation

```bash
pip install propka
```

## Python API

```python
import propka.run

# Run pKa prediction on a PDB file
mol = propka.run.single("receptor.pdb", write_pka=False)

# Access results
for conformation in mol.conformations:
    conf = mol.conformations[conformation]
    for group in conf.groups:
        if group.pka_value is not None:
            print(f"{group.residue_type}{group.atom.res_num} "
                  f"chain {group.atom.chain_id}: pKa = {group.pka_value:.2f}")
    break  # Usually only need first conformation
```

## Key Titratable Residues

| Residue | Model pKa | Typical shifted range | Effect at pH 7.4 |
|---------|-----------|----------------------|-------------------|
| ASP | 3.80 | 2.0 – 7.0 | Usually deprotonated (negative) |
| GLU | 4.50 | 3.0 – 7.5 | Usually deprotonated (negative) |
| HIS | 6.50 | 4.0 – 8.5 | Often partially protonated |
| CYS | 9.00 | 6.0 – 12.0 | Usually protonated (neutral) |
| TYR | 10.00 | 8.0 – 13.0 | Usually protonated (neutral) |
| LYS | 10.50 | 8.0 – 12.5 | Usually protonated (positive) |
| ARG | 12.50 | 11.0 – 14.0 | Always protonated (positive) |
| N-term | 8.00 | 6.0 – 9.0 | Usually protonated (positive) |
| C-term | 3.20 | 2.0 – 5.5 | Usually deprotonated (negative) |

## Interpreting Results

- **pKa > target pH** → residue is protonated at that pH
- **pKa < target pH** → residue is deprotonated at that pH
- **pKa near target pH** → residue in equilibrium, may need special treatment in docking

## Key Parameters

| Parameter | Description |
|-----------|-------------|
| `write_pka` | If True, writes .pka output file (set False for API use) |
| Target pH | Compare pKa against this to determine protonation state |

## Notes

- PROPKA requires standard PDB ATOM format (not OpenMM output format)
- Run PROPKA after missing-atom completion. The current Step 0 implementation uses the PDBFixer-completed receptor when available, so its report reflects the structure sent to downstream preparation.
- Treat the result as diagnostic evidence. PROPKA in this pipeline reports likely states but does not itself rewrite every residue protonation/tautomer state.
- Backbone amide groups (BB) show pKa = 0.0 — these are non-titratable
- For docking: pay attention to HIS residues near the binding site (HID vs HIE vs HIP)
- PROPKA accuracy: ±1.0 pKa unit for well-resolved structures

## Common Issues

1. **"no molecular conformations" error**: The PDB file is missing proper ATOM records or is in non-standard format. Ensure the file has standard PDB column formatting.
2. **Many BB groups at pKa 0.0**: Normal — these are backbone groups that PROPKA reports but are not titratable.
3. **Stream parameter**: Do NOT pass `stream=io.StringIO()` — this makes PROPKA read from the empty stream instead of the file.

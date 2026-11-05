# PDBFixer in this pipeline

PDBFixer repairs supported standard heavy atoms in structure_prep. It is not the
final protonation backend. The executable wrapper is `scripts/prep_utils.py`;
use the registry-driven runner rather than copying an independent API recipe.

| Parameter | Default | Contract |
|---|---|---|
| add_missing_atoms | true | Review whether reconstruction is defensible |
| add_missing_residues | false | Opt in with repair_reason; internal gaps at most 3 residues |
| add_hydrogens | false | Required; final hydrogens belong to docking_prep |
| seed | 0 | PDBFixer addMissingAtoms seed |
| pH | unset | Optional diagnostic context; no protonation here |

The writer uses `keepIds=true`. Multi-model input requires explicit model
selection before running. Unsupported terminal or long missing segments fail
instead of being silently rebuilt. Missing dependencies for a selected repair
fail the step. PDBFixer reconstruction does not establish refinement-quality
geometry; existing basic QC does not certify bond valence or clashes.

See [preparation contract](../../pipeline-entry/references/preparation-contract.md)
for noH preservation, provisional Reduce, final docking preparation and migration.

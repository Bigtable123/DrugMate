---
name: structure-prep
description: Establish a hydrogen-free receptor geometry baseline, preserve source structure and water inventory, extract the bound reference ligand, and optionally supply provisional hydrogens for pocket analysis.
---

# Structure Preparation

Read [the preparation contract](../pipeline-entry/references/preparation-contract.md)
before choosing tools. Enter full experiments through `../pipeline-entry/SKILL.md`;
the complete explicit matrix is validated and runtimes come from the registry.

Step 0 performs geometry preparation. Resolve the target entity, assembly, author
chain, coordinate model and reference ligand instance first. No chain A default
exists. Provide a reviewed local structure for unsupported assemblies or model
selection; the runner does not infer these choices.

PDBFixer defaults to `add_missing_atoms=true`, `add_missing_residues=false`,
`add_hydrogens=false`, `seed=0`. Adding residues requires `repair_reason` and is
limited to supported internal gaps of at most three residues. Terminal or longer
gaps are rejected. Selected but unavailable PDBFixer fails preparation.

PROPKA is off by default and diagnostic only; an explicit pH is required when
enabled. Reduce is off by default and needs a pocket-analysis reason. Its output
is a separate provisional derivative; it never overwrites the canonical noH
structure. Final protonation/flip, water selection, reference microstate diagnostics
and Vina Meeko conversion belong to mandatory `docking.docking_prep`. Candidate
Amber topology preparation belongs to the explicit downstream MM workflow.

Outputs are `prepared_receptor_noH`, `prep_manifest`, the original source complex,
`water_inventory`, bound `reference_ligand`, `metal_report`, basic
`structure_qc_report`, and optional `provisional_receptor` / provisional
`protonation_manifest`. Retain original water coordinates; do not preselect final
waters. Pocket/generation tools consume the noH geometry; hydrogen-dependent
pocket profiles may use the explicitly selected provisional derivative.

The current QC covers basic structure/file checks before and after repair.
Valence, clashes and reconstructed geometry have not been upgraded. Inspect
suspicious repairs and do not equate basic QC success with chemical validity.

Reference ligands with a CCD component ID must match the CCD chemical graph
after explicit hydrogen removal and bond-order assignment, while preserving
crystal heavy-atom order and coordinates. The report records
`template_graph_verified` and `crystal_coordinates_preserved`. Template loading,
assignment or validation failure stops preparation; it never falls back to
unverified PDB bond perception. This graph check does not certify stereochemistry
or conformer quality. Local ligands without a component ID retain the explicitly
reported `rdkit_pdb_perception` route and require chemistry review.

The current exact tool/parameter/runtime schema is generated in
`../pipeline-entry/TOOL_CATALOG.md`. Old Step 0 protonation, Meeko, water and
reference-Amber tool keys produce migration errors; rebuild the matrix and
validate before starting a new round. Do not overwrite old experiments.

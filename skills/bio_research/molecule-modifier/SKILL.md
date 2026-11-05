---
name: molecule-modifier
description: Generate local molecular variants from Agent-authored optimization instructions using RDKit SMARTS rules, preserving parent identity for an SBDD R3 campaign.
metadata:
  skill-author: SBDD Agent
  version: "1.0"
  dependencies: [rdkit]
  optional-dependencies: []
  external-tools: []
  python: ">=3.10"
---

# Molecule Modifier

Rule-based molecular modification engine. Reads Agent-authored instructions,
applies SMARTS transformations, validates local changes, and generates a 3D SDF.
New variants must pass property prefiltering and downstream evaluation before nomination.

## When to Use

- When the Agent selects MODIFY parents for local optimization based on R1/R2 evidence.
- To generate variants with explicit parent links for R3 comparison. Coordinator
  diagnostics are evidence, not a command to continue or an optimization target.

## Workflow

```
Agent modification_instructions.json (targets + per-molecule goals)
  + library_manifest.csv (compound_id → SMILES)
         │
         ▼
┌─────────────────────────────┐
│ 1. Parse fix_instruction    │  Keyword matching → rule selection
│    keywords per molecule    │
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│ 2. Apply SMARTS rules       │  1-3 variants per molecule
│    (reaction_rules.py)      │
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│ 3. Validate modifications   │  Validity + Tanimoto > 0.4
│    (RDKit)                  │  + atom-count change < 30%
└─────────────┬───────────────┘
              ▼
┌─────────────────────────────┐
│ 4. 3D Conformer + SDF       │  ETKDGv3 + MMFF/UFF
└─────────────────────────────┘
         │
         ▼
   modified_library.sdf  →  property_prefilter (order 3) → docking/evaluation
```

## Input Schema

```json
{
  "library_manifest_file": "step2/library_manifest.csv",
  "modification_instructions_file": "inputs/R3-a1/modification_instructions.json",
  "max_variants_per_mol": 3,
  "seed": 42
}
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `library_manifest_file` | string | Yes | — | Path to library_manifest.csv (compound_id + smiles columns) |
| `modification_instructions_file` | string | Yes* | — | Path to Agent-authored instructions JSON with top-level target_molecules |
| `modification_instructions` | object | Yes* | — | Inline instructions (alternative to file) |
| `max_variants_per_mol` | int | No | 3 | Max modified variants per input molecule |
| `seed` | int | No | 42 | Random seed for conformer generation |

*One of `modification_instructions_file` or `modification_instructions` required.

The instruction object uses the existing worker schema:

```json
{
  "target_molecules": ["ACT_0003"],
  "per_molecule": {
    "ACT_0003": {"fix_instruction": "Trim bulky substituent to reduce clash"}
  },
  "suggestions": ["Preserve the parent scaffold"]
}
```

Use real parent IDs present in the supplied manifest. Do not pass an entire
coordinator report: new reports have no authoritative `next_round_instructions`.
The engine matches keywords and does not enforce arbitrary natural-language
constraints or place a donor toward a named residue in 3D; verify every variant
against the Agent's objectives after generation and re-evaluation.

## Supported Modification Rules

| Keyword in fix_instruction | Transformation |
|---|---|
| "H-bond donor" | -CH3 → -CH2OH / -NH2; aromatic C-H → C-OH |
| "H-bond acceptor" | -CH3 → -CH2OCH3; aromatic C-H → C-OMe |
| "hydrophobic" | aromatic C-H → C-Me; -OH → -OC(C)C |
| "aromatic" | -CH3 → -CH2-phenyl; C-H → C-phenyl |
| "charged" | -CH3 → -CH2NH2; C-H → C-COOH |
| "strain" / "simplify" | Remove gem-dimethyl → mono-methyl |
| "clash" / "trim" | -CH2CH3 → -CH3; -OCH3 → -OH |

## Output Files

| File | Description |
|------|-------------|
| `modified_library.sdf` | 3D SDF for property prefiltering and subsequent docking |
| `modified_manifest.csv` | compound_id, smiles, parent_id, parent_smiles, modification_type |
| `modification_report.json` | Full report with counts and per-molecule details |

## Usage

```bash
# Replace the example paths with this campaign's actual artifacts.
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/molecule-modifier/run.py" \
  --input "$(python3 -c 'import json, os; print(json.dumps({"library_manifest_file": os.path.join(os.environ["SBDD_DATA_ROOT"], "project-actual/round1/step2/library_manifest.csv"), "modification_instructions_file": os.path.join(os.environ["SBDD_DATA_ROOT"], "project-actual/inputs/R3-a1/modification_instructions.json")}))')" \
  --output-dir "${SBDD_DATA_ROOT}/project-actual/modified/R3-a1"
```

Preserve `parent_id` and `parent_smiles` in the report, manifest and downstream
library; compound names alone can collide across runs. When combining R1/R2
parents, retain their source runs and reject ambiguous parent mappings.

## Integration with Pipeline

The Agent chooses targets, modifications and the next scientific plan:

1. Read R1/R2 `round_assessment`, per-candidate diagnoses and protocol validity.
   Choose actual optimization parents; do not modify molecules to compensate
   for an invalid docking protocol or a tool outage.
2. Write instructions and run this skill with the matching parent manifest.
3. Use the [campaign evidence helper](../pipeline-entry/references/campaign.md#collect-compare-and-assemble-r3-inputs)
   to combine historical KEEP, MODIFY parent baselines and successful variants.
4. Inject the new `library` and run `property_prefilter` (order 3) or earlier,
   followed by docking and full downstream evaluation. Do not use `start_step: 4`
   to skip prefiltering new chemistry. Freeze the R3 plan and check existing
   authorization as described in the [iteration guide](../pipeline-entry/ITERATIVE_GUIDE.md).

## Validation Rules

Each modified molecule must pass:
- `Chem.MolFromSmiles()` returns valid mol (not None)
- Heavy atom count change ≤ 30% of original
- Tanimoto similarity to parent > 0.4 (ensures it's a modification, not a new molecule)
- Successful 3D embedding (ETKDGv3)

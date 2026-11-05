# Candidate-specific topology preparation

`candidate_topology` is an opt-in canonical Step 5 escalation for MM-GBSA and
short-MD workflows. Whether it is enabled is decided in
`../../pipeline-entry/TOOL_SELECTION_POLICY.md`, written explicitly during
config assembly, and checked by the validator. This reference only defines its
execution and evidence contract.

## Contract

The evaluator must supply an SDF containing only QC-surviving candidates,
normally a ranked top-K subset. Each record is handled separately and retains:

- `chemical_id`, `source_identity_id`, and optional `source_identity`
- `pose_id`, `parent_pose_id`, and `coordinate_revision`
- legacy `compound_id`, `microstate_id`, and `variant_id` when present
- canonical isomeric SMILES and the RDKit formal charge actually passed to
  `antechamber -nc`

Canonical IDs are recomputed with the shared `pose_identity` contract. Legacy
identifiers are audit fields. An input value using a canonical `chem:`,
`source:`, or `pose:` prefix must exactly match the recomputed identity;
otherwise the record is `not_evaluable` with a lineage mismatch instead of
silently accepting stale metadata.

The implementation is
`scripts/candidate_topology.py::build_candidate_topologies`. It supports:

- `execution_mode=preflight`: identity/applicability validation and manifest
  only
- `execution_mode=build`: ff14SB or ff19SB receptor, GAFF2 ligand, and AM1-BCC
  charges using the registered Amber environment

Build mode writes separate dry complex, receptor, and ligand `prmtop`/`inpcrd`
files for every candidate. `candidate_topologies.json` records input hashes,
exact commands, bounded command logs, version probes, lineage, output paths,
output hashes, and the status of every requested record.

When `openmm_md` is enabled, the same build additionally writes a neutralized
TIP3P periodic octahedral `solvated_complex.prmtop/inpcrd`. The validated short
MD route fixes the protein force field to ff14SB and consumes those files in
the dedicated `openmm` environment. The dry artifacts remain the MM-GBSA
inputs. Topology generation or a bare OpenMM Context smoke test is still not a
pose-stability MD result; only the explicit simulation manifest is.

## Evidence semantics

`success` means the manifest was written. Scientific/tool outcome is always
reported separately as one of:

- `passed`: preflight passed, or all requested topology files were generated
- `failed`: AmberTools execution/output verification failed, or a requested
  batch contains only a partial evaluable subset
- `not_evaluable`: the generic protocol does not apply

Errors and unsupported chemistry must never be converted to a pass.

## Applicability boundaries

The generic route stops before Amber execution for:

- a declared covalent protein-ligand complex
- a declared or receptor-detected metal site
- retained, unparameterized receptor hetero residues unless an expert
  explicitly asserts they have already been handled

Individual ligand records are `not_evaluable` for unsupported elements,
multiple disconnected fragments, radicals, invalid sanitization, or missing
3D coordinates. Metal/covalent systems require a separately reviewed protocol
(for example, an MCPB.py-specific metal model or explicit covalent
parameterization); enabling this tool does not silently invent one.

## Adapter parameter shape

The fragment below documents this tool's parameters; it is not a complete
pipeline config or a tool-selection recommendation:

```yaml
steps:
  pose_evaluation:
    tools:
      candidate_topology:
        enabled: true
        params:
          execution_mode: build
          max_candidates: 20
          protein_force_field: ff14SB
          ligand_force_field: gaff2
          ligand_charge_method: bcc
```

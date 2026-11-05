# Metric extension controls

Use the numeric bounds/child weights in the [metric policy](../../METRIC_POLICY.md)
and its executable tables; retain the existing parent
weights. Exact computation definitions are in [METRIC_POLICY.md](../../METRIC_POLICY.md).
All tools below are registered. New switches are off until explicitly selected.
Old complete pipeline configurations need explicit entries for the new tools;
generate the selection matrix from the current registry/tool catalog rather
than reusing an incomplete older matrix.

| Step / tool | Produced evidence | Cost and prerequisites |
| --- | --- | --- |
| property_prefilter or property_filter / `rascore` | SYN probability; hard/soft lower 0.2/0.5, +3 | CPU batch, installed XGB model; no route search |
| pose_evaluation / `rdkit_pose` | INTRA strain, TFD, low-energy-reference RMSD | CPU conformers; default 8 ETKDG samples, seed 42, 500 iterations; unsupported chemistry NA |
| pose_evaluation / `freesasa` | POCKET buriedness and ligand-side contact-surface fraction, absolute SASA audit fields | CPU, current receptor and pose coordinates; probe 1.4 Å |
| pose_evaluation / `prolif` | FP Tanimoto to reference, subtype counts and residue bits | Existing profiler; requires a usable same-receptor reference for Tanimoto only |
| pose_evaluation / `arpeggio` | INTER clash/vdW-clash/proximal fraction; raw contact records | Isolated CPU worker per pose, derived complex mmCIF; quality-degraded outputs leave numeric metrics missing |
| pose_evaluation / `luna` | Interaction records and hashed fingerprint | Isolated CPU worker per pose; audit only, no extra numeric score |
| pose_evaluation / `oddt` | Residue IFP, simple IFP and SPLIF | Isolated CPU worker per pose; audit only, no extra numeric score |
| pose_evaluation / `rtmscore` | DYN same-protocol percentile, +3 | Existing three-model inference; CPU or GPU; no absolute cutoff |
| pose_evaluation / `candidate_topology` + `openmm_md` + `mdanalysis` | DYN ligand RMSD median, COM-distance drift P95, mean initial-contact retention | Existing parameterization/simulation; defaults top-K 3 and 100 ps; no new simulation is launched just by loading rules |

The following is a **configuration fragment**, not a complete runnable plan.
Merge these choices into a complete explicit step/tool matrix and supply the
usual scientific inputs. Existing `prolif`, `rtmscore`, and MD switches remain
independent; the expensive MD/RBFE options should reflect the intended experiment.

```json
{
  "steps": {
    "property_filter": {
      "tools": {"rascore": {"enabled": true, "params": {"timeout_seconds": 600}}}
    },
    "pose_evaluation": {
      "tools": {
        "rdkit_pose": {"enabled": true, "params": {"num_conformers": 8}},
        "freesasa": {"enabled": true, "params": {"probe_radius": 1.4}},
        "arpeggio": {"enabled": false},
        "luna": {"enabled": false},
        "oddt": {"enabled": false}
      }
    }
  }
}
```

Runtime scales with the evaluated pose count. Arpeggio/LUNA/ODDT invoke isolated
workers per pose and are more expensive than simply folding an existing field.
Run them on an intentionally selected pose set rather than assuming a large
library will have the same cost as one smoke-test complex. Tool failures are
missing evidence and operational diagnostics, never a fabricated zero/pass.

Not implemented here: actual OpenFE/TIES/TIES_MD RBFE execution/analysis,
H-bond or predeclared interaction occupancy across trajectories, ΔRMSF against
apo/reference, restraint occupancy, a composite short-MD stability score,
multi-frame MMGBSA/SEM/convergence, and the remaining providers/geometric
endpoints listed in the integration audit. Successful OpenFE/TIES preflight or
planning is not a measured free energy.

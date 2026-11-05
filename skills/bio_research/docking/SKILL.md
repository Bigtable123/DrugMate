---
name: docking
description: Dock prepared molecules with one selected engine.
metadata:
  skill-author: SBDD Agent
  version: "1.0"
  dependencies: []
  optional-dependencies: [rdkit]
  external-tools: [gnina, unidock, unidocktools, vina, meeko, diffdock, dock6, antechamber, autodock4, autogrid4, autodocktools, rdock, adfr, agfr]
  python: ">=3.10"
---

# Docking

Read [the preparation contract](../pipeline-entry/references/preparation-contract.md)
first. Step 4 runs mandatory `docking_prep` on the canonical noH receptor and
prefiltered 3D ligand library before invoking an engine. Explicitly
select one final protonation backend, pH when applicable, water policy and
justifications. Final receptor and preparation manifests are separate artifacts;
pose evaluation consumes this final receptor. Before receptor tools run,
`docking_prep` checks every ligand record for usable coordinates, including
ordinary SDF inputs without state metadata. Empty/unreadable records, missing
or non-finite coordinates, collapsed coordinates and zero-length bonds are
reported in `docking_prep_qc.json` under `ligand_geometry_qc` and prevent engine
dispatch. Regenerate affected molecules from original declared SMILES and verify
specified stereochemistry; re-embedding a collapsed SDF can change the isomer.
A preparation failure prevents engine dispatch. Basic QC does not certify
valence, stereochemistry or full conformer quality.

For routine non-covalent protein receptors without a target-pH assignment
requirement, prefer Reduce with `build=true` and `flip=true` after runtime
preflight. Choose PDB2PQR+PROPKA for target-pH or protonation-sensitive protocols,
with an explicit, justified pH and hydrogen-bond optimization enabled. The
preparation contract covers retained-component compatibility; preserve an
explicit user choice and do not switch backends automatically after failure.

Step 4 executes exactly one of eight primary docking engines selected in config.
GNINA keeps its CNN-based scores. Uni-Dock supports GPU-batched
non-covalent screening and keeps its native lower-is-better `docking_score`.
AutoDock Vina provides a reproducible classic CPU baseline on small or
medium libraries and keeps its native `vina_*` energy fields. DiffDock is the
blind generative-pose mode and keeps its confidence-model output. DOCK 6 is
the classic grid-scoring mode and requires injected, precomputed receptor-site
spheres plus the matching `.nrg`/`.bmp` grid pair.
AutoDock4 is the classic Lamarckian-GA mode and owns its AutoGrid4 and
AutoDockTools helpers. rDock supplies an independent cavity/grid-scoring mode.
ADFR supplies the ADFRsuite maps-based search mode and owns AGFR plus its
native preparation helpers. The helpers are deliberately not separate primary
engines, enable switches, or score sources.

`search_space_source` is also required in config. Runtime never chooses it from
whichever artifact happens to be present and never falls back to another source.

Before finalizing the engine selection, follow the
[unified memory protocol](../pipeline-entry/references/benchmark-memory.md) and
[docking branch](../pipeline-entry/references/docking-memory.md)
for `context.docking_memory_mode` (default `ligand_protein`). The Agent runs the
selector after valid query inputs are ready and before docking, then records
the metric, proposed and actual engine in the plan. The default objective is
Top-1 RMSD <2 Å; the five supported objectives are separate selections, not
five simultaneous guarantees. Explicit user choices take priority.
The selector handles one query ligand; the first record of a screening library
must not stand in for the whole library. This runner executes the approved
config and never queries memory or changes engine switches.

## When to Use

- After molecule library generation (Step 2) to dock compounds into the target pocket
- When you have a prepared receptor PDB and a compound library SDF
- Use GNINA when you need CNN pose/affinity rescoring
- Use Uni-Dock when throughput is the priority and native Vina/Vinardo scores are acceptable
- Use AutoDock Vina for a classic CPU baseline or an explicit Vina request on a
  small-to-medium library; it is not a strongly orthogonal validator for GNINA
  or Uni-Dock because all three belong to the Vina method family
- Use DiffDock when a blind generative pose proposal is explicitly requested;
  its confidence is not binding affinity
- Use DOCK 6 for a classic non-Vina grid-docking mode when trusted receptor
  spheres and energy/bump grids are available
- Use AutoDock4 for the native ADT → AutoGrid4 → AutoDock4 classic baseline
- Use rDock for a reference-ligand or explicit sphere-mapped cavity mode; its
  score is an rDock-native unit, not kcal/mol
- Use ADFR for a native preparation → AGFR target → ADFR maps-based mode;
  flexible-residue and covalent modes are not yet part of this pipeline contract

## Input Schema

```json
{
  "receptor_noH_file": "path/to/prepared_receptor_noH.pdb",
  "receptor_file": "path/to/prepared_receptor_noH.pdb",
  "ligand_file": "path/to/library.sdf",
  "reference_ligand_file": "path/to/reference_ligand.sdf",
  "docking_box_file": "path/to/docking_box.json",
  "search_space_source": "reference_ligand",
  "is_covalent": false,
  "tools": {
    "docking_prep": {"enabled": true, "params": {
      "reason": "Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment",
      "water_policy": "dry", "water_reason": "Reviewed dry receptor protocol",
      "ligand_state_policy": "preserve_declared",
      "ligand_state_reason": "Library states were prepared for the selected site"
    }},
    "reduce": {"enabled": true, "params": {"build": true, "flip": true}},
    "pdb2pqr": {"enabled": false},
    "protoss": {"enabled": false},
    "water_selection": {"enabled": false}, "meeko": {"enabled": false},
    "ligand_microstates": {"enabled": false}, "openbabel": {"enabled": false},
    "gnina": {
      "enabled": true,
      "params": {"exhaustiveness": 16, "num_modes": 20}
    },
    "unidock": {"enabled": false},
    "autodock_vina": {"enabled": false},
    "diffdock": {"enabled": false},
    "dock6": {"enabled": false},
    "autodock4": {"enabled": false},
    "rdock": {"enabled": false},
    "adfr": {"enabled": false},
    "redocking_validation": {
      "enabled": true,
      "params": {
        "redocking_rmsd_threshold": 2.0,
        "redocking_rmsd_hard_threshold": 2.5
      }
    }
  }
}
```

This Reduce example does not assign protonation states at a target pH. For that
requirement, explicitly select PDB2PQR+PROPKA and justify `docking_prep.params.ph`.
PDB2PQR invokes PROPKA internally; the separate `structure_prep.propka` diagnostic
is not required for that route.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `receptor_file` | string | Yes | - | Canonical hydrogen-free protein geometry |
| `receptor_pdbqt_file` | string | Conditional | - | Produced internally after final receptor preparation by `docking.meeko` for Vina; old injected PDBQT cannot bypass preparation |
| `ligand_file` | string | Yes | - | Compound library SDF (multi-mol) |
| `reference_ligand_file` | string | Conditional | - | Required only when `search_space_source=reference_ligand`; used to define the box for box-based engines. An incidental reference artifact never overrides another source. |
| `docking_box_file` | string | Conditional | - | Required when `search_space_source=docking_box`; pocket-analysis JSON containing finite `center[3]` and positive `size[3]`. An incidental box never overrides another source. |
| `search_space_source` | string | Yes | - | Explicit source: `reference_ligand`, `docking_box`, `whole_receptor`, or `precomputed_grid`. There is no inferred or standalone default. |
| `tools.gnina.enabled` | bool | Yes | - | Explicit GNINA switch; exactly one primary engine must be enabled in the resolved config |
| `tools.unidock.enabled` | bool | Yes | - | Explicit Uni-Dock switch; GNINA and redocking validation must be disabled when true |
| `tools.autodock_vina.enabled` | bool | Yes | - | Explicit classic CPU AutoDock Vina switch; GNINA, Uni-Dock, and redocking validation must be disabled when true |
| `tools.diffdock.enabled` | bool | Yes | - | Explicit blind generative DiffDock switch; every other primary engine and redocking validation must be disabled when true |
| `tools.dock6.enabled` | bool | Yes | - | Explicit classic DOCK 6 grid-docking switch; every other primary engine and redocking validation must be disabled when true |
| `tools.autodock4.enabled` | bool | Yes | - | Explicit AutoDockTools + AutoGrid4 + AutoDock4 switch; helpers remain inside this one engine mode |
| `tools.rdock.enabled` | bool | Yes | - | Explicit rDock cavity-generation and docking switch |
| `tools.adfr.enabled` | bool | Yes | - | Explicit ADFRsuite preparation + AGFR + ADFR switch |
| `dock6_receptor_site_file` | string | DOCK 6 | - | Injected DOCK 6 receptor-site sphere file (`.sph`) |
| `dock6_grid_nrg_file` | string | DOCK 6 | - | Injected DOCK 6 energy grid (`.nrg`) |
| `dock6_grid_bmp_file` | string | DOCK 6 | - | Injected DOCK 6 bump grid (`.bmp`) sharing the `.nrg` prefix |
| `tools.*.params` | object | No | registry | Engine-native parameters; see the generated tool catalog |
| `tools.redocking_validation.enabled` | bool | Yes | - | Explicit GNINA native-ligand RMSD protocol switch; currently GNINA-specific |
| `is_covalent` | bool | No | false | Safety flag forwarded by the orchestrator; true hard-fails until a covalent protocol is wired |

The `tools` block must contain every registered docking tool with an explicit
boolean `enabled`: `docking_prep`, `reduce`, `pdb2pqr`, `protoss`,
`water_selection`, `meeko`, `ligand_microstates`, `openbabel`, all eight primary
engines, and `redocking_validation`. Enable `docking_prep` and exactly one
supported final receptor backend, with the preparation parameters shown above.
Vina additionally requires `meeko`. Exactly one of
`gnina`, `unidock`, `autodock_vina`, `diffdock`, `dock6`, `autodock4`, `rdock`,
and `adfr` must be enabled in the resolved config. Source
compatibility is explicit: DiffDock requires `whole_receptor`; DOCK 6 requires
`precomputed_grid`; the six box-based engines accept only `reference_ligand` or
`docking_box`. Missing source-specific artifacts fail before engine execution.

By default, request up to **20 poses per molecule**. After ligand preparation,
the request is capped by the number of input molecules: 20 for at most 50,
5 for 51–1000, and 1 above 1000. The same cap applies to `num_modes` for
GNINA, Uni-Dock, AutoDock Vina, DOCK 6, AutoDock4, rDock, and ADFR, and to
`samples_per_complex` for DiffDock. Users may request fewer poses in the
experiment config. Record requested and returned counts; downstream evaluation
uses the same cohort-based cap by default.

## Workflow

```
receptor.pdb + library.sdf + explicit search_space_source
                              │
          ┌───────────────────┼──────────────────────┐
 reference_ligand / docking_box  whole_receptor  precomputed_grid
          │                         │                │
 six box-based engines          DiffDock          DOCK 6
          └─────────────────────────┴────────────────┘
                              │
                   exactly one primary engine
                              ▼
             docked_poses.sdf + scores/report JSON
```

GNINA, Uni-Dock, AutoDock Vina, AutoDock4, rDock, and ADFR use exactly one
declared box source. `reference_ligand` requires the reference ligand and
`docking_box` requires an explicit box; the runtime does not substitute one for
the other. DiffDock performs blind whole-receptor inference only. DOCK 6
consumes its injected site/grid as the actual search-space definition; it never
fabricates a receptor grid from the prepared PDB.

## Output Files

| File | Description |
|------|-------------|
| `docked_poses.sdf` | All docked poses with native engine score properties |
| `docking_scores.json` | Per-pose native scores plus `engine` and pose lineage |
| `docking_report.json` | Run parameters, status, file paths |
| `docking_molecule_ledger.json` | Lossless per-input accounting: `scored`, `failed`, or `input_invalid`, with pose counts and failure reasons |
| `docking_molecule_ledger.tsv` | Tabular view of the same per-input accounting for Agent/report inspection |
| `gpu_diagnostics.json` | GNINA-only startup GPU snapshot |
| `gpu_usage_timeseries.json` | GNINA-only runtime GPU samples |
| `unidock_native/` | Uni-Dock-only audited per-ligand native outputs |
| `vina_native/` | AutoDock Vina-only prepared ligand and native PDBQT outputs |
| `diffdock_native/` | DiffDock-only isolated inputs and ranked native predictions |
| `dock6_native/` | DOCK 6-only ligand preparation, input deck, and scored Mol2 |
| `autodock4_native/` | AutoDock4-only PDBQT, GPF/DPF, maps, and DLG files |
| `rdock_native/` | rDock-only receptor PRM/cavity and native scored SD files |
| `adfr_native/` | ADFR-only prepared PDBQT, AGFR target, DLG/DRO, and native poses |
| `redocked_reference.sdf` | GNINA native-ligand redocking poses when applicable |

### Canonical pose score fields

Native and derived fields appear in `docked_poses.sdf` and/or
`docking_scores.json` according to the selected adapter:

| Property | Description |
|----------|-------------|
| `minimizedAffinity` | Vina-like binding affinity (kcal/mol, lower = better) |
| `CNNscore` | CNN pose quality score (0-1, higher = better) |
| `CNNaffinity` | CNN-predicted binding affinity on a pK scale (higher = stronger predicted affinity) |
| `intramol`, `CNNvariance` | Optional GNINA audit fields; weight 0 |
| `docking_score` | Uni-Dock Vina/Vinardo energy (kcal/mol, lower = better) |
| `vina_affinity` | AutoDock Vina/Vinardo total energy (kcal/mol, lower = better) |
| `heavy_atom_count`, `ligand_efficiency` | AutoDock Vina derived size normalization inputs/output |
| `engine_pose_rank`, `score_gap_1_2` | Audit-only mode ordering fields when available; weight 0 |
| `vina_inter`, `vina_intra`, `vina_torsions`, `vina_intra_best_pose` | AutoDock Vina energy components; diagnostic context, not interchangeable ranking fields |
| `diffdock_confidence` | DiffDock confidence-model output (higher is better); not affinity and not cross-engine comparable |
| `dock6_grid_score` | DOCK 6 native grid score (lower is better); compare only within one DOCK 6 protocol |
| `dock6_grid_vdw_energy`, `dock6_grid_es_energy` | DOCK 6 diagnostic grid components; context only |
| `autodock4_binding_energy` | AutoDock4 estimated free energy of binding (kcal/mol, lower is better) |
| `rdock_score` | rDock native total score (lower is better; not labeled kcal/mol) |
| `rdock_inter_score`, `rdock_intra_score`, `rdock_restraint_score` | rDock native components; diagnostic context |
| `adfr_affinity` | ADFR native affinity/energy (kcal/mol, lower is better) |

Do not compare score directions interchangeably: `minimizedAffinity` and
Uni-Dock `docking_score` and AutoDock Vina `vina_affinity` are
lower-is-better; `CNNscore` and `CNNaffinity` are higher-is-better. Apply the
protocol-level native-redocking gate below only when GNINA,
`search_space_source=reference_ligand`, and `redocking_validation` are all
explicitly selected. Uni-Dock and AutoDock Vina do not currently have an
engine-neutral redocking protocol, so their redocking gate must be disabled and
is `not_applicable`. For every other search-space source, native redocking is
also `not_applicable`, not a pass; judge the run using pocket geometry and the
remaining pose-quality evidence.

DiffDock confidence, DOCK 6 grid score, AutoDock4 energy, rDock score, and ADFR
affinity remain in separate score families. No raw value from different engine
families may be ranked against another without an independently calibrated,
declared protocol.

### Executable score gates and weights

Apply these fallback gates only within one comparable target/pocket protocol.
They are early-screening sanity bounds, not experimental binding-affinity
claims.

These are **level-1 per-metric marks**: each raw value is labeled `pass`,
`flag`, or `fail`. They do not independently keep, modify, or reject a
molecule. A finite score with no declared bound remains `not_applicable` for the
gate while still being rankable. After all planned metrics in a parent category have reached a terminal value/NA state, the
shared scorer counts `fail` marks and applies the selected V1/V2/V3 category
fold to produce the molecule-level `keep`, `modify`, or `reject` decision.

| Metric | Pass mark | Flag mark | Fail mark | Suggested signed weight |
|---|---|---|---|---:|
| GNINA `minimizedAffinity` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -5 |
| GNINA `CNNscore` | `>= 0.50` | `[0.20, 0.50)` | `< 0.20` | +4 |
| GNINA `CNNaffinity` | `>= 6` pK | `[2, 6)` | `< 2` | +3 |
| Uni-Dock `docking_score` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -3 |
| AutoDock Vina `vina_affinity` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -5 |
| AutoDock Vina ligand efficiency, `-vina_affinity / HeavyAtomCount` | `>= 0.25` | `[0.05, 0.25)` | `< 0.05` | +2 |
| DiffDock `confidence` | N/A (finite audit value retained) | no absolute tier | no molecular gate | 0 |
| DOCK 6 `grid_score` | N/A (finite audit value retained) | no absolute tier | no molecular gate | 0 |
| AutoDock4 `binding_energy` | `<= -7` kcal/mol | `(-7, -1]` | `> -1` | -4 |
| rDock `score` | N/A (finite same-protocol rank retained) | no absolute tier | no molecular gate | -5 |
| ADFR `affinity` (current registry alias of the AutoDockFR flexible row) | `<= -6` kcal/mol | `(-6, 0]` | `> 0` | -4 |

The supplied table uses negative suggested weights for lower-is-better. The
public executable policy may expose both magnitude and signed weight; scoring
uses the signed value shown above. The positive ligand-efficiency runtime field
`-vina_affinity/HAC` with weight `+2` is algebraically equivalent to the source
table's negative `vina_affinity/HAC` with weight `-2`. DiffDock and
DOCK 6 are absent from that table, so their runtime confidence/grid values are
weight-zero audit evidence unless a session-local `0..5` override explicitly
activates them. The canonical `AutoDockFR` row is a flexible-receptor score
(`soft_upper=-6`, `hard_upper=0`, signed weight `-4`, magnitude `4`). The
registered adapter is rigid-receptor ADFR, but the current metric registry maps
`adfr.affinity` onto that row and applies those bounds and weight. This matches
the source numbers but not the row's flexible-receptor applicability. Treat it
as a documented implementation conformance gap; do not describe a rigid ADFR
result as validated flexible-receptor evidence.

The Step 5 parent-category weights are `AFF=5`, `DRAW=3`, `FLEX=3`, and
`BENCH=0`. Child weights above are independent of those parent weights.

Ligand efficiency is derived from the exact pose score and the ligand heavy
atom count; it is not a native Vina field, although the adapter publishes the
derived `heavy_atom_count` and `ligand_efficiency`. Pose/mode rank, score gap,
GNINA `intramol`, and `CNNvariance` have weight 0. Rank `<=20` remains inside
the declared export window; rank `>20` is a hard pose-export rejection. It is
remapped to `modify` for the parent compound, allowing a rank-`<=20` alternate
pose or a rerun instead of rejecting the chemistry. Preserve the other audit
fields, or use a weak score-gap warning as
described in the engine reference, but never add their raw values to the
weighted total.

For every scored cohort:

1. require the same target, receptor/pocket, engine and version, scoring
   function, preparation/protonation policy, search box, and material search
   settings;
2. mark every raw value against its hard/soft bounds first, then fold the marks
   per parent category; fold-rejected poses do not enter any percentile pool and
   receive no score/rank;
3. calculate descending RANK.AVG percentile independently for each finite
   metric among the surviving poses (`n=1 → 0.5`; ties use average rank); missing
   values leave only that pose's numerator and denominator;
4. combine signed child weights into 0–100 AFF/DRAW category scores, percentile-
   rank the category scores, then apply parent weights to Step 5 StepScore;
   weight-zero values remain mark/audit evidence and contribute exactly zero;
5. do not stack homologous Vina-family energies. GNINA
   `minimizedAffinity`, Uni-Dock `docking_score`, and AutoDock Vina
   `vina_affinity` occupy one family slot; use the predeclared primary engine's
   value rather than summing results from multiple Vina-family engines; the
   effective family child-weight total remains capped at 5.

The adapter publishes a `protocol_id` built from engine/software identity,
scoring/search parameters, search space, and content SHA-256 identities of the
receptor representation actually consumed by that engine, the ligand batch,
the optional ligand-preparation manifest, and the reference ligand when used.
Paths are retained for audit but do not change the fingerprint. Unknown engine
version, missing required content identity, or mixed protocol IDs block trusted
cohort ranking.

GNINA `CNNscore` and `CNNaffinity` are distinct CNN evidence and may accompany
the selected GNINA energy within a GNINA cohort. They must not be renamed,
copied into another engine's records, or compared as if they shared units.

### Native-redocking protocol gate

The current validation mode is GNINA-specific and requires
`search_space_source=reference_ligand`. It computes the minimum
symmetry-aware heavy-atom RMSD across all generated redocked reference-ligand
modes in the receptor coordinate frame, without independently aligning the
ligand (`CalcRMS`, not `GetBestRMS`): it is a **best-of-generated-modes** diagnostic, not score-top-1 accuracy
and not a top-N active/decoy success rate. Redocking reuses the main run's
`num_modes`, `pose_sort_order`, `cnn_scoring`, exhaustiveness, device, and seed;
changing any of them defines a different protocol.

| Best redocking RMSD | Protocol action |
|---|---|
| `<= 2.0 Å` | keep protocol |
| `(2.0, 2.5] Å` | modify/recalibrate protocol |
| `> 2.5 Å` | reject and roll back the target×pocket×protocol |

This gate applies to the docking protocol, never as evidence to reject an
individual candidate molecule. A disabled or inapplicable redocking run is NA,
not a pass and not part of the weighted molecular score. A redocking attempt
that crashed, timed out, exited non-zero, or produced no atom-compatible pose
measured no RMSD, so it publishes the same NA shape
(`gate="not_applicable"`, `passed=null`, `measured=false`) with `performed=true`
plus a `reason` — never `reject`, which would be indistinguishable from a
genuine RMSD past the hard threshold and would roll the protocol back on a
benchmark that never produced a number.

Downstream pose-score collection recognizes the engine-native fields above.
Native scores remain namespaced and must never be copied into fake CNN fields.

## Dependencies

| Tool | Install | Notes |
|------|---------|-------|
| GNINA | Download binary from GitHub releases | Requires CUDA for GPU acceleration |
| Uni-Dock + Uni-Dock Tools | Dedicated `unidock` Conda environment | Tools wrapper is mandatory so SDF ligands receive flexible-topology properties |
| AutoDock Vina Python API | Dedicated `vina` Conda environment | CPU engine; current adapter exposes Vina/Vinardo scoring |
| Meeko + RDKit | Registered `reduce` environment | Step 4 `docking.meeko` creates PDBQT from the current final receptor; the Vina adapter prepares SDF ligands and exports poses |
| RDKit | `pip install rdkit` | Required for symmetry-aware native-redocking RMSD; score parsing alone has a text fallback |
| DiffDock | Prefix environment `${SBDD_ENV_DIFFDOCK_PREFIX}` and source `${SBDD_DIFFDOCK_SOURCE_DIR}` | Requires local score/confidence checkpoints; adapter performs blind whole-receptor inference |
| DOCK 6 + Antechamber | Prefix environment `${SBDD_ENV_DOCK6_PREFIX}`, binaries under `${SBDD_DOCK6_SOURCE_DIR}/bin` | Requires injected trusted `.sph`, `.nrg`, and `.bmp` receptor assets; the adapter never fabricates a grid from PDB |
| AutoDock4 + AutoGrid4 + AutoDockTools | Prefix `${SBDD_ENV_AUTODOCK4_PREFIX}` | One composite mode; native ADT preparation does not use `docking.meeko` |
| rDock + Open Babel | Prefix `${SBDD_ENV_RDOCK_PREFIX}` | Adapter sets `RBT_ROOT` and an isolated per-run `RBT_HOME` |
| ADFRsuite (ADFR + AGFR + preparation tools) | Prefix `${SBDD_ENV_ADFRSUITE_PREFIX}` | One composite mode using the suite's native wrappers |

## Usage

Choose one example for the reviewed task. These standalone inputs assume a
routine non-covalent receptor without target-pH assignment, a reviewed dry
protocol and an already prefiltered 3D library. Replace the example paths and
justifications with the current inputs and approved preparation choices. Each
example declares the complete docking tool matrix; Vina enables `meeko` within
Step 4. For pH-sensitive tasks, use the PDB2PQR preparation contract above.

Set `SBDD_PYTHON` and `BIO_SKILLS_ROOT` through the deployment configuration.
Run from the allocated request directory under `SBDD_DATA_ROOT`, outside the
checkout, so these relative input and output paths resolve within that request.
Use a separate output directory for each execution.

```bash
cd "${SBDD_REQUEST_DIR:?Set SBDD_REQUEST_DIR to the allocated request directory}"
# Standard docking
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","reference_ligand_file":"inputs/reference_ligand.sdf","search_space_source":"reference_ligand","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":false},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":true,"params":{"num_modes":20}},"unidock":{"enabled":false},"autodock_vina":{"enabled":false},"diffdock":{"enabled":false},"dock6":{"enabled":false},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":true}}}' \
  --output-dir ./step4/

# Higher exhaustiveness for final poses
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","reference_ligand_file":"inputs/reference_ligand.sdf","search_space_source":"reference_ligand","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":false},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":true,"params":{"exhaustiveness":32,"num_modes":20}},"unidock":{"enabled":false},"autodock_vina":{"enabled":false},"diffdock":{"enabled":false},"dock6":{"enabled":false},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":true}}}' \
  --output-dir ./step4_final/

# Explicit pocket-analysis box
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","docking_box_file":"inputs/docking_box.json","search_space_source":"docking_box","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":false},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":true,"params":{"num_modes":20}},"unidock":{"enabled":false},"autodock_vina":{"enabled":false},"diffdock":{"enabled":false},"dock6":{"enabled":false},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":false}}}' \
  --output-dir ./step4_box/

# Uni-Dock: explicitly turn off GNINA and its GNINA-specific redocking protocol
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","docking_box_file":"inputs/docking_box.json","search_space_source":"docking_box","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":false},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":false},"unidock":{"enabled":true,"params":{"search_mode":"fast","num_modes":20}},"autodock_vina":{"enabled":false},"diffdock":{"enabled":false},"dock6":{"enabled":false},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":false}}}' \
  --output-dir ./step4_unidock/

# AutoDock Vina: Step 4 Meeko prepares the current final receptor PDBQT
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","docking_box_file":"inputs/docking_box.json","search_space_source":"docking_box","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":true},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":false},"unidock":{"enabled":false},"autodock_vina":{"enabled":true,"params":{"cpu":8,"exhaustiveness":16,"num_modes":20}},"diffdock":{"enabled":false},"dock6":{"enabled":false},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":false}}}' \
  --output-dir ./step4_vina/

# DiffDock: blind generative poses (declared pocket is provenance context only)
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","search_space_source":"whole_receptor","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":false},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":false},"unidock":{"enabled":false},"autodock_vina":{"enabled":false},"diffdock":{"enabled":true,"params":{"samples_per_complex":20,"device":0}},"dock6":{"enabled":false},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":false}}}' \
  --output-dir ./step4_diffdock/

# DOCK 6: inject a matched receptor-site/energy-grid/bump-grid set
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/docking/run.py" \
  --input '{"receptor_file":"inputs/prepared_receptor_noH.pdb","ligand_file":"inputs/prefiltered_library.sdf","search_space_source":"precomputed_grid","dock6_receptor_site_file":"inputs/site.sph","dock6_grid_nrg_file":"inputs/grid.nrg","dock6_grid_bmp_file":"inputs/grid.bmp","tools":{"docking_prep":{"enabled":true,"params":{"reason":"Reduce hydrogen/flip preparation for a routine non-covalent site without target-pH assignment","water_policy":"dry","water_reason":"Reviewed dry receptor protocol","ligand_state_policy":"preserve_declared","ligand_state_reason":"Prefiltered library states were prepared for the selected site"}},"reduce":{"enabled":true,"params":{"build":true,"flip":true}},"pdb2pqr":{"enabled":false},"protoss":{"enabled":false},"water_selection":{"enabled":false},"meeko":{"enabled":false},"ligand_microstates":{"enabled":false},"openbabel":{"enabled":false},"gnina":{"enabled":false},"unidock":{"enabled":false},"autodock_vina":{"enabled":false},"diffdock":{"enabled":false},"dock6":{"enabled":true,"params":{"num_modes":20}},"autodock4":{"enabled":false},"rdock":{"enabled":false},"adfr":{"enabled":false},"redocking_validation":{"enabled":false}}}' \
  --output-dir ./step4_dock6/
```

## Error Handling

### Filtered inputs and retries

Stage 04 libraries and ordinary external 3D SDFs are accepted without
candidate-state tags or a seal. Stage 04 standardizes molecules and prepares
conformers; it does not claim exhaustive pH microstates. Preserve compound IDs,
source records, conformer identity, and cohorts when filtering or preparing
retry batches. `microstate_id` may be null when no state expansion was performed.

Legacy candidate-state tags are ordinary input metadata and do not trigger
enumeration, ensemble-completeness checks, or a manifest requirement.

Docked molecules must preserve the input SDF's molecular graph and specified
stereochemistry. Output coordinates may resolve unspecified input E/Z or
tetrahedral stereo; canonical isomeric SMILES need not be identical. Changed
connectivity, bond order, charge, isotopes, or lost/inverted specified stereo
still fail validation. Do not rewrite input stereochemistry to force string
equality. Exact pose/score and coordinate-file identity checks remain required.

A sanitization, molecular-identity or coordinate failure is isolated only when
the raw record and every associated score map unambiguously to the same input.
All poses of that molecule are quarantined; surviving poses pass the complete
artifact contract again before publication. Unknown identities, broken framing,
score mismatches and zero remaining scored poses still fail the batch.
`docking_recovery.json` retains the original file hash, excluded physical record
indices and the published-to-raw index map. `raw_pose_index` identifies the raw
record; `pose_index` and `source_pose_index` address the published SDF.
Failed molecules remain operational failures in the ledger, without molecular
scores or chemical rejection. Raw docking and failed redocking coordinates
remain in `.docking_uncommitted` for audit; they are not downstream inputs.
Score coverage distinguishes physical, invalid, quarantined and published
records; `complete` describes score fields on readable poses only.

Choose stable `compound_id` values for a retry from the audited docking ledger,
then prepare a new input directory:

```bash
python "${BIO_SKILLS_ROOT}/docking/scripts/select_states.py" \
  --input-file "${SBDD_DATA_ROOT}/example/round1/logs/docking_docking_input.json" \
  --state-ids-file "${SBDD_DATA_ROOT}/example/retry_state_ids.json" \
  --output-dir "${SBDD_DATA_ROOT}/example/retry_inputs"
```

The IDs file is a JSON array. The helper checks unambiguous IDs, preserves
source record order and all SDF metadata, and writes `docking_input.json`,
`selected_ligands.sdf` and `selection_manifest.json`. It does not run an engine
or infer completion from an incomplete pose file. Use the resulting input with
`docking/run.py` and a new output directory, retaining previous attempts. A
changed ligand batch has its own protocol fingerprint; keep provenance when
combining results and do not treat a successful retry as full-run completion.

Retain the original source alongside filtered molecules.

| Scenario | Behavior |
|----------|----------|
| `tools` matrix or any `enabled` boolean omitted | Worker fails before inspecting artifacts or selecting an engine |
| Selected engine unavailable | Runtime preflight or worker returns `success: false`; disabled engines are not probed or run |
| Zero or multiple primary engines selected | Validation and worker both fail before docking |
| `search_space_source` omitted or incompatible with the selected engine | Validation and worker fail before docking; no source is inferred |
| Input files missing | `success: false` with file path |
| AutoDock Vina selected without `docking.meeko` | Validation or preparation fails; enable Meeko in Step 4 to prepare the current final receptor PDBQT |
| GNINA timeout (>2h) | `success: false` |
| No poses generated | `success: false` |
| Any input molecule has no published finite pose | It remains in `docking_molecule_ledger.json/.tsv` with stable input index/identity, explicit failure reason, and `retryable`; the accounting invariant is `scored + failed + input_invalid == input` |
| GNINA emits any valid SDF pose without a finite primary GNINA score | The complete poses publish; the shortfall is itemized in `score_coverage` and declared in `operational_errors` for the failure registry. Only a batch with no parseable scored pose at all fails |
| `reference_ligand` selected without a valid reference ligand | `success: false` before the selected engine starts; an available box is not a fallback |
| `docking_box` selected without a valid explicit box | `success: false` before the selected engine starts; an available reference ligand is not a fallback |
| RDKit not available to GNINA score parsing | Fallback to text-based SDF parsing |
| AutoDock Vina Meeko/RDKit runtime unavailable | Runtime preflight hard-fails; there is no lossy format-conversion fallback |
| `nvidia-smi` unavailable during run | Docking continues; GPU usage samples record the error |
| Uni-Dock silently skips a molecule | The molecules it did score publish; the skipped ones carry no pose, are counted in `diagnostics.unscored_molecules`, and are declared in `operational_errors` |
| AutoDock Vina ligand preparation or engine stage fails for an input record/variant | The batch fails before canonical pose export; native artifacts and diagnostics are retained |
| AutoDock Vina cannot export an input record/variant | If usable scored poses remain, export retains them and records failed records in `export_diagnostics.failures` plus incomplete coverage in `operational_errors`; zero usable exported poses fails. The common docking artifact checks still apply |
| DiffDock skips a molecule/rank or emits a non-finite confidence | The molecule lands in `diagnostics.incomplete_molecules`/`invalid_records` and is declared in `operational_errors`; the remaining cohort still publishes. Only zero usable poses fails the run |
| DOCK 6 site/grid artifacts absent, empty, or have mismatched prefixes | Validation/worker hard-fails before docking; no receptor grid is synthesized |
| DOCK 6 skips a molecule or emits an invalid/non-finite scored Mol2 pose | The valid poses publish; `diagnostics.invalid_records` / `unscored_source_molecule_indices` name the gap and `operational_errors` declares it. Only zero finite scored poses fails the run |
| AutoDock4, rDock, or ADFR omits a molecule or emits a missing/non-finite primary score | The covered molecules publish; `diagnostics.uncovered_source_molecule_indices` / `failed_molecules` name the rest and `operational_errors` declares it. Only zero parseable scored poses fails the run |
| `is_covalent=true` | Fails before ordinary non-covalent docking; protonation preparation alone is not treated as a covalent protocol |

## Upstream / Downstream

### Upstream
| From | File | Usage |
|------|------|-------|
| structure-prep | `<receptor_stem>_noH.pdb` | Canonical hydrogen-free receptor input |
| structure-prep | `reference_ligand.sdf` | Required only for the explicitly selected `reference_ligand` source |
| pocket-analysis | `docking_box.json` | Required only for the explicitly selected `docking_box` source |
| molecule-generator | `library.sdf` | Ligand library |
| injected DOCK 6 assets | `*.sph`, paired `*.nrg`/`*.bmp` | Trusted classic-docking site and receptor grid; `precomputed_grid` with DOCK 6 only |

For Vina, Step 4 produces its own `receptor_pdbqt` through `docking.meeko` after
final receptor preparation; it is not a Step 0 input artifact.

### Downstream
| To | File | Usage |
|----|------|-------|
| structure-evaluator | `docked_poses.sdf` | Pose quality assessment |
| structure-evaluator | `docking_scores.json` | Score integration |
| Agent / final report | `docking_molecule_ledger.json` | Audit every input molecule, including candidates that produced no pose and therefore cannot enter pose ranking |

## References

- `references/gnina.md` — GNINA CLI options, scoring details, and known issues
- `references/unidock.md` — topology-aware invocation, output contract, and scale boundary
- `references/autodock_vina.md` — Step 4 receptor preparation/PDBQT contract, CPU adapter, native scores, and scientific limitations
- `references/diffdock.md` — blind generative inference, confidence semantics, and completeness checks
- `references/dock6.md` — injected receptor-grid contract, ligand preparation, and native grid-score semantics
- `references/autodock4.md` — ADT/AutoGrid4/AutoDock4 workflow and DLG score semantics
- `references/rdock.md` — cavity mapping, environment variables, and native score semantics
- `references/adfr.md` — ADFRsuite preparation, AGFR targets, and ADFR limitations

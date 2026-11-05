---
name: pocket-analysis
description: Analyze binding pockets and fixed docking boxes.
---

# Pocket Analysis

Read [the preparation contract](../pipeline-entry/references/preparation-contract.md).
Use canonical noH geometry or an explicitly justified provisional Reduce derivative.
Target chain must be selected from entity/assembly evidence; example A is not a default.
For an already prepared receptor and fixed box, unique chain overlap can resolve
an omitted chain as described below; it does not establish biological identity.
Interaction evidence precedes final docking protonation and requires final-state review.
To obtain a ProLIF fingerprint, supply a reference ligand and an explicit-H
analysis receptor: enable `structure_prep.reduce` with `build=true` and a reason,
or inject an audited hydrogenated `provisional_receptor`. When either input is
unavailable, record ProLIF as not evaluable or failed and continue with the
other evidence. The analysis copy remains separate from the noH baseline.

Define protein binding pockets, analyze reference-ligand interactions when available, and emit a normalized docking-box contract for downstream SBDD steps.

This remains one pipeline **Step 1**: it covers site identification/geometric characterization and, when a reference ligand exists, normalized interaction observations plus an Interaction Blueprint. It does not implement a druggability-scoring branch or call paid, commercial, or Web-only services.

For an unknown site, use `site_selection_strategy=detector_consensus` in the main
Step 1. Run `benchmark_memory.py plan-pocket-tools` before detection. Compatible
protein memory retrieves five distinct protein-group neighbors, whose best
detectors form a five-entry multiset from DeepPocket, fpocket, P2Rank, pyKVFinder,
and AutoSite. Enable exactly the returned `execution_tools`: each distinct
recommendation executes once, with DeepPocket always included as the default
anchor. Step 1 validates this explicit matrix before execution. It normalizes
native Top-3 candidates to point5 and applies `protein_memory_local_fusion_v1`.
A six-condition guard may switch to P2Rank#1 when P2Rank is among the selected
tools. Every matching Top-3 pocket contributes reliability/rank weights with
family de-correlation; output residues remain the anchor's exact point5 set.
With unavailable memory or explicit `pocket_memory_mode=off`, the plan selects
the four-tool fixed DeepPocket#1 fallback, with AutoSite disabled, and reports
the reason. `required` rejects unavailable memory. AutoSite is integrated through
its native CSV/cluster-point adapter in `adfrsuite` and is currently supported
only for `detector_consensus`.

The selection policy records its implementation choices: minimum
native-top-1 DCA defines each historical best detector, ties use lexicographic
tool IDs, all five historical measurements must be present, repeated
recommendations execute once, DeepPocket is forced as anchor, and all matching
Top-3 candidates contribute. Method checks do not establish experimentally
measured performance for a new query.
See the [memory-fusion contract](../pipeline-entry/references/pocket-memory.md)
for support-point provenance and independent replay. User-provided boxes and
explicitly reviewed sites retain priority. Do not preselect a fused pocket and
route it back through `provided_box`, which would recompute its residue set.

## When to Use

- After protein preparation (typically after `structure-prep` skill)
- Before molecular docking to define the binding box and constraints
- When you need to locate and rank candidate binding pockets by detector consensus
- When you need a cross-checked interaction profile or fingerprint of a reference ligand
- When generating interaction constraints or fingerprint-based filters for virtual screening

## Input Schema

```json
{
  "protein_file": "path/to/protein_prepared.pdb",
  "protein_noH_file": "path/to/protein_noH.pdb",
  "complex_mmcif_file": "path/to/source_complex.cif",
  "ligand_file": "path/to/reference_ligand.sdf",
  "site_selection_strategy": "ligand_site",
  "ligand_resname": "LIG",
  "chain": "A",
  "known_site_residues": ["Ile10", "Val18", "Leu83"],
  "detector_box_padding": 4.0,
  "known_required_interactions": [
    {
      "residue": "Leu83", "interaction_type": "H-bond", "chain": "A",
      "protein_atom_names": ["N"], "ligand_role": "acceptor", "protein_role": "donor",
      "source": "PMID:12345678", "rationale": "Conserved hinge anchor supported by SAR",
      "binding_mode_scope": "ATP-competitive inhibitors", "review_status": "approved"
    }
  ],
  "fail_on_profiler_error": false,
  "fail_on_supplemental_profiler_error": false,
  "blueprint_policy": {"auto_required_from_tool_consensus": false, "require_curated_provenance": true},
  "fpocket_params": {
    "min_radius": 3.0,
    "max_radius": 6.0
  },
  "tools": {
    "ligand_site": {"enabled": true},
    "fpocket": {"enabled": true},
    "pykvfinder": {"enabled": true},
    "p2rank": {"enabled": true},
    "autosite": {"enabled": false},
    "deeppocket": {
      "enabled": false,
      "params": {
        "max_candidates": 3,
        "segmentation_threshold": 0.5,
        "mask_distance": 3.5,
        "timeout_seconds": 1200
      }
    },
    "pocketminer": {
      "enabled": false,
      "params": {
        "rerank_mode": "annotate",
        "min_score_coverage": 0.8,
        "provisional_quantile": 0.9,
        "cluster_ca_cutoff": 8.0
      }
    },
    "plip": {"enabled": true},
    "prolif": {"enabled": true},
    "arpeggio": {"enabled": false},
    "luna": {"enabled": false},
    "oddt": {"enabled": false},
    "blueprint": {"enabled": true}
  }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `protein_file` | string | Yes | Receptor PDB; a successful ProLIF audit uses an explicitly hydrogenated analysis copy (`provisional_receptor`) |
| `protein_noH_file` | string | No | Canonical geometry baseline retained for provenance; never a ProLIF fallback |
| `complex_mmcif_file` | string | No | Source complex mmCIF from Step 0; required only by enabled/applicable Arpeggio and must retain `_chem_comp` ligand chemistry |
| `ligand_file` | string | Conditional | Reference ligand SDF/PDB. Required by `site_selection_strategy=ligand_site` and enabled gate-bearing ligand profilers; ProLIF records `not_evaluable` if absent. Its presence never changes the selected strategy. |
| `site_selection_strategy` | string | Yes | Explicit mode: `ligand_site`, `known_site`, `detector_consensus`, or `provided_box`. There is no inferred or standalone default. |
| `docking_box_file` | string | Conditional | Required by `provided_box`: canonical JSON with finite `center[3]`, positive `size[3]`, and angstrom units. Pipeline supplies `input_artifacts.docking_box`. |
| `ligand_resname` | string | No | Residue name of ligand (auto-inferred from filename if omitted) |
| `chain` | string | Conditional; no default | Reviewed target author chain. `provided_box` can resolve an omitted value only when exactly one named receptor chain overlaps the box. |
| `known_site_residues` | array | Conditional | Required by `site_selection_strategy=known_site`; used for recall validation. Its presence never changes the selected strategy. |
| `center_match_threshold` | float | No (default: 5.0) | Other site strategies may configure this; `detector_consensus` fixes the center-distance limit at 5 Å. |
| `residue_iou_threshold` | float | No (default: 0.3) | Other site strategies may configure this; `detector_consensus` fixes point5 residue IoU at 0.30. |
| `consensus_top_k` | integer | No (default: 3) | `detector_consensus` always uses native Top-3; other site strategies may configure their candidate count. |
| `detector_box_padding` | float | No (default: 4.0) | Padding in Å around detector-consensus residues when publishing the docking box. Reduce only for an explicit engine grid-size constraint; the resolved plan records the choice. |
| `known_required_interactions` | array | No | Curated hard constraints; each requires source, rationale, binding-mode scope, and approved review status |
| `fail_on_profiler_error` | boolean | No (default: false) | Legacy compatibility field; ProLIF failures are audit-only. An enabled PLIP worker failure currently invalidates the Step 1 contract regardless of this field. |
| `fail_on_supplemental_profiler_error` | boolean | No (default: false) | Legacy compatibility field; enabled Arpeggio/LUNA/ODDT worker failures currently invalidate the Step 1 contract regardless of this field. |
| `blueprint_policy` | object | No | Keep automatic tool-consensus hard requirements disabled by default |
| `fpocket_params` | object | No | fpocket parameters (`min_radius`, `max_radius`) |
| `pykvfinder_params` | object | No | pyKVFinder parameters; configured independently from fpocket |
| `pocket_memory_mode` | string | No | `auto` (default), `off`, or `required`; unavailable memory falls back only in `auto` |
| `pocket_memory_index` | string | No | Separate derived pocket DCA/protein index; not the docking ligand/protein index |
| `pocket_memory_query_npz` / `pocket_memory_query_metadata` | string | For memory | Compatible current receptor encoding and bound `pocket_prediction` metadata |
| `pocket_memory_query_case` / `pocket_memory_query_group` | string | No | Known query identities; exclude self and protein group |
| `pocket_memory_exclude_cases` / `pocket_memory_exclude_groups` | array | No | Additional history exclusions applied before retrieval and priors |
| `tools` | object | Yes | Complete matrix containing every registered pocket tool; every entry must have an explicit boolean `enabled`. Missing entries or values fail before execution. |
| `tools.ligand_site/fpocket/pykvfinder/p2rank` | object | Yes | Explicit built-in/geometric detector switches; artifact presence never changes them. |
| `tools.deeppocket` | object | Yes | Required enabled for `detector_consensus`; isolated DeepPocket classifier/segmentation. Its proposals originate from fpocket, so the two share one independence group. |
| `tools.autosite` | object | Yes | Explicit switch for AutoSite native cluster detection in the registered `adfrsuite` runtime. Enable only when the `detector_consensus` plan selects it; disable for all explicit-site strategies. |
| `tools.pocketminer` | object | Yes | Explicit isolated PocketMiner switch. Its scores are auxiliary cryptic evidence, not cavity geometry. |
| `tools.plip/prolif` | object | Yes | Explicit profiler switches. PLIP requires `ligand_file`; ProLIF is record-only and skips when it is absent. |
| `tools.arpeggio` | object | Yes | Supplemental atom-contact adjudication in `drugmate`; enabling it requires a chemically complete source mmCIF and reference ligand. |
| `tools.luna` | object | Yes | Supplemental atom-aware interactions and hashed IFP in the dedicated `luna` Conda environment. |
| `tools.oddt` | object | Yes | Legacy residue IFP/SPLIF evidence in the dedicated `oddt` Conda environment; never authorizes Blueprint constraints. |
| `tools.blueprint` | object | Yes | Explicit Interaction Blueprint compilation switch. |
| `native_step_timeout_seconds` | integer | No (default: 300) | Timeout for each isolated interaction-profiler worker |

## Workflow

```
prepared receptor + explicit site_selection_strategy
  ├─ ligand_site       → observed ligand shell, center, and box
  ├─ known_site        → detector geometry checked against known residues
  ├─ provided_box      → characterize the supplied box without resizing it
  └─ detector_consensus→ validate memory + plan detector tools before execution
                              │
                              ├─ memory: five recommendations → distinct tools + DeepPocket
                              │    repertoire: DeepPocket / fpocket / P2Rank / pyKVFinder / AutoSite
                              │    off/unavailable: DeepPocket / fpocket / P2Rank / pyKVFinder
                              │    explicit tool matrix must match the plan
                              ├─ run planned detectors → native Top-3 → point5
                              │    memory: DCA reliability + guarded anchor switch
                              │    fallback: DeepPocket#1 + reciprocal native rank
                              │    distance ≤5 Å AND residue IoU ≥0.30
                              │    center fusion with family de-correlation
                              │    retain the selected anchor's point5 residues
                              │
                              ├─ optional PocketMiner (separate Conda env)
                              │    residue enrichment annotation only
                              │    provisional clusters only when no cavity is found
                              │
                              ▼
            pocket_definition.json + docking_box.json
                              │
                              └─ explicitly enabled ligand profilers
                                   ├─ PLIP: Blueprint evidence; ProLIF: audit-only
                                   ├─ optional supplemental: Arpeggio / LUNA / ODDT
                                   └─ observations + fingerprints + blueprint
```

See the corresponding files under `references/` for each tool's parameters, API usage, and caveats. Read `references/interaction-profilers.md` before selecting, installing, or routing interaction-analysis tools.

---

## Explicit Site Selection and Consensus

`site_selection_strategy` is a required decision made in pipeline config. The runtime validates and executes that decision; it never derives the strategy from `ligand_file`, `known_site_residues`, or another artifact.

- **`provided_box`**: use when the user supplies a fixed docking box. Keep Step 1
  enabled in the same experiment. Select residues with at least one heavy atom
  inside the box (inclusive boundaries) on the reviewed target chain, using the
  noH baseline. Preserve chain/insertion codes and the exact box center/size.
  Check overlap before running detectors or profilers. An explicit chain with
  no overlap is a configuration error; inspect the reported per-chain counts
  and correct the choice after reviewing the target/assembly. Never silently
  replace it with another chain. If chain is omitted, select only the unique
  named chain with overlap and record the evidence; zero/multiple matching
  chains require a reviewed selection. Do not default to A or select whichever
  chain has most atoms. In a pipeline, normalization can supply this geometric
  selection in `config_patch.context.chain`; merge the patch and review it
  before validation. Preserve necessary assembly partners in the receptor.
  No detector is required; explicitly enabled detectors remain supplemental and
  cannot move the site or box. Interaction profilers still require a reference
  ligand and their ordinary prerequisites. Geometry validation checks finite,
  positive box dimensions and nonempty receptor overlap; it does not establish
  cavity quality or binding, so confidence is `not_assessed` and cavity volume
  remains unavailable. The report binds the box and receptor paths/SHA-256 and
  selected chain. `files.docking_box` references the input file; Step 1 does not
  publish a replacement canonical box. Downstream docking must use
  `search_space_source=docking_box`. Do not create an auxiliary pocket project.
- **`ligand_site`**: requires `ligand_file` and the `ligand_site` tool enabled. It uses the observed ligand to define the primary residue shell, center, and box. Enabled detectors may provide supporting evidence and alternatives, but no predictor may replace the selected ligand-defined site.
  The final docking box is rebuilt from reference-ligand heavy atoms: its center is the midpoint of the coordinate minima/maxima, and each side is `max(ligand_extent + 2 * box_padding, min_box_size)` (defaults: 4 Å per side and 12 Å minimum). The ligand atom centroid remains the pocket center; it need not equal the box center. Residue shells, detector unions, and cached candidate boxes do not set the ligand-site box extent, including when padding changes or cached geometry is invalid. Missing/invalid reference-ligand geometry prevents canonical publication; it never triggers a residue-box fallback. The output records `geometry_source=reference_ligand_heavy_atoms`; `bounds` stores the unpadded ligand bounds.
- **`known_site`**: requires non-empty `known_site_residues`. Run the enabled geometry detectors, choose by their declared consensus procedure, and validate that the selected geometry recovers the known residues. The residue list does not itself manufacture a box.
- **`detector_consensus`**: run the [pocket planner](../pipeline-entry/references/pocket-memory.md#integration-with-main-step-1), then explicitly enable exactly its `execution_tools` in the complete matrix. Retrieval and priors require complete historical native-top-1 DCA for all five detectors, even if the current plan uses fewer tools. Normalize the planned tools' native Top-3 support points to point5 on the current noH receptor. With memory, include the anchor and every candidate within 5 Å **AND** point5 residue IoU ≥0.30, including matching lower ranks from the anchor's detector. With unavailable/off memory, run the four fallback detectors and retain at most one match from each non-anchor detector, ordered by greatest IoU, shortest distance, then native rank. Members need not match one another. Top-K/threshold overrides cannot alter the policy.
- **Center and residue contract**: memory uses `v_t=(p_hat_detector+0.05)^2/rank_t` for each member pocket; fallback uses `v_t=1/rank_t`. Families are `{DeepPocket, fpocket}`, `{P2Rank}`, `{pyKVFinder}`, and `{AutoSite}`; each family's total weight is `max(v_t)` and is divided among all its members in proportion to `v_t`. Normalize effective weights for the center average. No additional match means the selected anchor's original center. Publish exactly that anchor's point5 residues, preserving chain/insertion codes; never union members' residues or recompute them around the fused center. Native probabilities enter the guarded anchor decision, not the fusion weights.
- **DeepPocket**: required for `detector_consensus`, optional for explicit-site strategies. Run its real fpocket-proposal → CNN-ranking → CNN-segmentation workflow in the dedicated `deeppocket` environment. Its point5 support is the original fpocket proposal PQR mapped by classifier rank; segmentation outputs remain validated provenance and do not replace those support points. DeepPocket and fpocket share `fpocket_family` and cannot alone create two independent votes.
- **AutoSite**: use the registered `adfrsuite` runtime and `scripts/autosite_utils.py`. Preserve native ranks from `receptor_summary.csv`; map exported `receptor_cl_<rank>.pdb` points to point5 residues. Retain cluster-point hashes and input-receptor provenance. AutoSite is a separate detector family; its native score is not a cross-tool fusion weight. Keep `tools.autosite.enabled=false` for `ligand_site`, `known_site`, and `provided_box`.
- **Availability**: report failures and preserve successful raw outputs. Every planned detector must provide valid native rank 1; naturally fewer than three pockets are allowed. Missing runtimes, failed detectors, invalid candidates, or corrupted memory fail the procedure instead of silently selecting the fallback or renumbering candidates. Only off/unavailable memory produces the fixed four-tool plan.
- **Cryptic evidence (opt-in)**: PocketMiner emits per-residue cryptic-opening propensities, not pockets. It never joins the geometry-detector `one_of`, increments detector `support_count`, supplies a box, or changes the selected pocket. Candidate enrichment uses a threshold-free Mann–Whitney probability-of-superiority AUC for annotation only, including when a legacy config specifies `rerank_mode=tie_break`.
- **No static cavity**: PocketMiner may report high-score CA-connected residue clusters as `provisional_cryptic_hypotheses`. These are visualization/sampling hypotheses only: they cannot produce canonical `pocket_definition.json` or `docking_box.json` until a geometry-bearing method (for example a detector on an opened conformation or an MD pocket analysis) confirms a cavity.
- **Strict runtime gate**: an explicitly enabled PocketMiner run hard-fails Step 1 if the repository commit, complete checkpoint-object match, residue mapping, finite/range checks, prediction count, or deterministic repeated inference fails. An upstream process exiting zero is not sufficient evidence of a valid checkpoint restore.
- **No usable site**: fail the pocket-definition contract. Do not invent a docking center.

Residues retain chain ID and insertion code in the normalized output so similarly numbered residues from different chains cannot collapse into one site.

**Cross-tool comparison rule (all site-selection strategies):** Never average unmatched centers, add heterogeneous scores, or average volumes defined by different algorithms. DeepPocket's score is not directly comparable with P2Rank or fpocket scores.

- Unmatched predictions remain separate candidates; never resolve their disagreement by publishing a midpoint. Detector-center averaging is restricted to members admitted by the declared consensus matching rule. In `ligand_site`, detector centers never replace the reference-ligand geometry.
- Preserve native scores with tool provenance. For unknown sites, use the planned memory policy and its six-condition P2Rank anchor guard, or the fixed DeepPocket#1 fallback when memory is off/unavailable. Independent support alone cannot authorize an anchor switch. Other explicit-site strategies retain their declared selection rules.
- Preserve native cavity volumes separately. A published scalar cavity volume must identify its single source in `volume_method`; do not average or sum volumes across tools. The residue-envelope volume is a separately labeled descriptor, not a replacement cavity volume.

---

## Interaction Fingerprint

Step 1 uses PLIP for Blueprint evidence and ProLIF for records only; Arpeggio,
LUNA, and ODDT are supplemental profilers disabled by default. Save each tool's
raw results separately, then aggregate them into `interaction_observations.json`.
These tools remain part of the same `pocket_analysis` Step 1; do not create a
separately resumable Step 1b.

Canonical observations retain the tool/version, chain/iCode, protein/ligand
atoms or features, chemical roles, mediators, geometry, `blueprint_eligible`,
and quality flags. Tool statuses distinguish `detected`, `not_detected`,
`not_evaluable`, and `tool_failed`. Cross-tool consensus requires compatible
atoms/features, roles, and geometry; matching only the same residue and a
coarse interaction type is insufficient.

| Profiler | Runtime | Step 1 role | Blueprint policy |
|---|---|---|---|
| PLIP | `drugmate` isolated worker | Core atom-aware static profiler | Eligible for automatic Desired/Optional constraints |
| ProLIF | `drugmate` isolated worker | IFP and geometry audit | `blueprint_eligible=false`; observations and PLIP agreement are recorded only |
| Arpeggio | `drugmate` CLI worker | Supplemental atom-contact adjudication from source mmCIF | Eligible only when chemistry QC passes; OpenBabel kekulization/aromaticity warnings mark output degraded and observation-only |
| LUNA | `luna` core-API worker | Supplemental atom-aware interactions and hashed/count IFP | Eligible for Desired evidence when atom roles/geometry are evaluable; weak H-bonds remain observation-only |
| ODDT | `oddt` Python worker | Supplemental residue IFP and SPLIF comparison | Observation/fingerprint only; residue-level channels and SPLIF never compile into Blueprint constraints |

Arpeggio must read the original complex mmCIF retained by Step 0; a renamed PDB
or a converted file lacking `_chem_comp` chemical definitions is not a substitute.
Without an mmCIF, the status is `not_evaluable`, not a tool failure. LUNA runs
headlessly through its core API in the dedicated `luna` environment, without
invoking the PyMOL-dependent CLI/PSE route. If the SDF contains only heavy atoms,
the ligand's H-bond donor capability must be marked not evaluable; missing
information must not be treated as a true negative. ODDT uses the validated
OpenBabel backend. Its 8-channel IFP is a count vector, not a bit vector; exclude
the metal channel from interpretation. Use SPLIF only for relative ranking or
similarity, never as a hard gate.

ProLIF also produces a sparse bit vector for each protein residue against the
reference ligand. If explicitly enabled without `ligand_file`, record
`not_evaluable` and skip ProLIF. Do not treat an empty fingerprint as negative
evidence or switch the site-selection strategy.

```json
{
  "interaction_fingerprint": {
    "residues": ["LEU83.A", "LYS89.A", "ALA31.A"],
    "interaction_types": ["HBDonor", "HBAcceptor", "Hydrophobic", "VdWContact"],
    "bitvector": [[1,0,0,1], [0,0,1,1], [0,0,0,1]],
    "summary": {
      "LEU83.A": ["HBDonor", "VdWContact"],
      "LYS89.A": ["Hydrophobic", "VdWContact"],
      "ALA31.A": ["VdWContact"]
    }
  }
}
```

ProLIF 2.x also retains geometric metadata for each interaction (distance, angle, and atom indices). Use these data for audits and review; they do not automatically compile into interaction constraints.

Failure handling has separate levels. Except for record-only ProLIF, explicitly
enabled tools must pass pipeline runtime preflight; missing Conda environments,
imports, or CLIs block entry into Step 1. A missing ProLIF runtime is recorded in
`optional_missing` and execution continues. ProLIF analysis failures are recorded
in `profiler_failures` and the failure registry without vetoing Blueprint,
pocket geometry, or canonical artifacts. PLIP analysis failures continue to
follow the existing Step 1 contract. Structure-specific worker or analysis
failures in Arpeggio/LUNA/ODDT are recorded in `supplemental_profiler_failures`
and likewise follow the existing Step 1 contract. `fail_on_profiler_error` and
`fail_on_supplemental_profiler_error` are retained for configuration compatibility
and do not change these runtime decisions. PocketMiner declares
`counts_as_detector_support=false`, so its failure does not veto the contract.

Publish `pocket_definition.json` and `docking_box.json` when the Step 1 contract
passes. If geometry is invalid, retain `proposed_pocket_definition` /
`proposed_docking_box` for diagnostics. A ProLIF `tool_failed` status does not
invalidate a passing geometry or Blueprint contract.

---

## Interaction Blueprint Generation

Blueprint defines design constraints, rather than raw detection results. By
default, treat computational results as observations first, then compile three
levels of constraints according to the strength of the evidence:

| Level | Classification rule | Confidence |
|------|---------|-----------|
| **Required** | `known_required_interactions` provides provenance, a mechanistic/SAR rationale, binding-mode scope, and `review_status=approved` | high |
| **Desired** | Blueprint-eligible PLIP, LUNA, or chemistry-QC-passing Arpeggio observations fall within the actual multi-detector interaction-consensus pocket | medium |
| **Optional** | Blueprint-eligible observations outside the pocket, weak VdW contacts, or water bridges/metal coordination without Step 0 evidence of retained mediators | low/conditional |

**Special rules**:

- Retain PLIP∩ProLIF agreement in the observation artifact for audit; it does not generate `candidate_required` or other Blueprint constraints.
- Arpeggio/LUNA/ODDT do not add votes for `candidate_required`. LUNA and chemistry-QC-passing Arpeggio can strengthen at most Desired evidence; degraded Arpeggio, LUNA weak H-bonds, and all ODDT records remain only in the evidence artifact.
- Manually authorize mechanistically required interactions as Required through `known_required_interactions`; `known_site_residues` alone cannot authorize a hard gate.
- The legacy-compatible `auto_required_from_tool_consensus=true` is permitted only when known-site evidence is independent of the reference structure and has human approval; it is off by default.
- WaterBridge/Metal-complex can be promoted to Desired only when the Step 0 report explicitly retains the corresponding class of water/metal; otherwise it remains conditional Optional.
- Without approved Required constraints, mark Blueprint as unusable for hard gating; evaluation must not pass by treating 0/0 as 100%.
- With only one available tool, evidence can reach at most Desired; distinguish a tool failure from genuine `not_detected` evidence.
- Without a reference ligand or known interaction evidence, explicitly mark the blueprint and downstream `blueprint_match` as `not_applicable`; N/A counts neither as a pass nor as a quality-gate failure.
- Blueprint internals and the matcher retain chain/iCode, but PLIP's interaction API often exposes only `restype`/`resnr`/`reschain`. If the original PLIP record lacks an insertion code, iCode-specific hard constraints must rely on independently curated evidence or manual review. ProLIF may retain audit records but is not a gating basis. The matcher does not treat a missing iCode as a wildcard.

**Downstream scoring boundaries (see the [metric policy](../METRIC_POLICY.md) and its machine-readable source `shared/metric_scoring.py`)**:

- `desired_ranking_score >= 0.50` is `pass`; `0.10 <= score < 0.50` is `flag`; `score < 0.10` is `fail`.
- `desired_evaluable_coverage` has only `soft_lower=0.80`, with no hard lower bound. Values below `0.80` are `flag`; at `<0.30`, stop trusted ranking and request reevaluation. Coverage alone must not reject a molecule or pose.
- Each downstream Required/Desired/Optional result first produces a per-metric mark. If any aggregate metric meets the selected V1/V2/V3 rejection condition, return `reject`. Otherwise, return `keep` only when strictly more than half of the marks within every evaluated aggregate metric are `pass`; return `modify` for the remainder. ProLIF records do not enter these counts.

**Output format**:

This abbreviated JSON-with-comments illustration omits optional details;
actual output contains the complete records and uses plain JSON.

```jsonc
{
  "schema_version": 2,
  "required": [
    {
      "constraint_id": "hinge_anchor_1",
      "residue": "LEU83",
      "chain": "A",
      "interaction_type": "H-bond",
      "protein_role": "donor",
      "ligand_role": "acceptor",
      "protein_atoms": [{"name": "N"}],
      "evidence_source": ["curated", "PLIP", "ProLIF"],
      "confidence": "high",
      "provenance": {
        "source_type": "literature",
        "reference": "PMID:12345678",
        "rationale": "Conserved hinge anchor supported by SAR"
      },
      "scope": {"binding_mode": "ATP-competitive inhibitors"},
      "review_status": "approved",
      "hard_gate_authorized": true
    }
  ],
  "required_groups": [
    {"group_id": "hinge_anchor_1", "logic": "all_of", "members": ["hinge_anchor_1"], "at_least_k": 1}
  ],
  "desired": [/* desired constraints omitted */],
  "optional": [/* optional constraints omitted */],
  "evidence_summary": { /* evidence details omitted */ },
  "generation_policy": { /* generation-policy details omitted */ }
}
```

**Implementation**: `scripts/blueprint_utils.py` → `generate_blueprint()`

---

## Pocket Features Output

```json
{
  "pocket_features": {
    "volume": 450.0,
    "volume_unit": "angstrom3",
    "volume_method": "fpocket_selected_pocket",
    "residue_envelope_volume": 1250.0,
    "burial_fraction_proxy": 0.72,
    "key_subpockets": [
      {
        "id": "hinge_cleft",
        "residues": ["Glu81", "Phe82", "Leu83", "His84"],
        "property": "H-bond rich"
      },
      {
        "id": "hydrophobic_back",
        "residues": ["Ile10", "Val18", "Ala31", "Val64", "Leu134"],
        "property": "hydrophobic"
      }
    ],
    "flexibility_zones": [
      {
        "residues": "Gly-rich loop (11-18)",
        "note": "moderate flexibility, B-factor > 30"
      }
    ]
  }
}
```

Feature calculation methods:

| Feature | Method |
|------|------|
| Volume | Use only the selected fpocket/pyKVFinder cavity's volume; unavailable with P2Rank alone |
| Residue envelope volume | Convex hull of pocket-residue heavy atoms; a spatial-envelope reference, not a cavity-volume estimate |
| Burial fraction proxy | Heavy-atom neighbor-shell approximation, explicitly labeled as a proxy |
| Sub-pockets | Cluster by residue properties (polar / hydrophobic / charged / aromatic) |
| Flexibility zones | Regions with B-factor > 30 |

---

## Validation

### 1. Strategy-aware pocket validation

| `site_selection_strategy` | Validation evidence |
|---|---|
| `ligand_site` | Ligand coverage, ligand-center distance from the pocket/box center, and whether the box encloses the ligand |
| `known_site` | `known_site_residues` recall (default target `>= 0.8`), independent geometric-detector provenance, and valid box geometry |
| `detector_consensus` | Five-entry recommendation multiset and exact execution-tool plan (or explicit four-tool fallback reason), point5 provenance for selected tools' Top-3 candidates, memory/guard/inference evidence, local matching and reliability/rank/family fusion, retained anchor residues, and valid box geometry; the anchor center is allowed when no other candidate matches |
| `provided_box` | Valid user box values in Å and target-chain heavy atoms inside the box; validates input geometry only, without requiring detector consensus or establishing that the binding site is correct |

Every strategy requires a three-dimensional `center` with finite values, a strictly positive three-dimensional `size`, and nonempty pocket residues before passing the definition to docking. If strategy validation fails, Step 1 returns top-level `success: false`, retains diagnostics only in the report's `proposed_*` fields, and neither writes nor publishes canonical `pocket_definition` / `docking_box` artifacts.

### 2. Interaction validation

- Key active-site residues should appear in the ProLIF fingerprint (for example, Leu83 in a kinase hinge region).
- If H-bonds are missing, first check whether `protein_prepared.pdb` contains explicit hydrogens.
- ProLIF checks the actual selected amino-acid residues for explicit H/D atoms;
  water or ligand hydrogens do not satisfy this check. Missing protein hydrogens
  produce `validation_failed`, never a successful empty/partial fingerprint.
  `protein_hydrogen_count` and `input_strategy` describe the actual input; they
  do not certify complete protonation or the final docking state.
- If the fingerprint is empty, check that ligand and protein coordinates use the same coordinate frame.

### 3. Blueprint completeness for manual or rule-based construction

- Required: ≥ 1 interaction confirmed by literature, SAR, or a crystal structure; otherwise hard gating is prohibited.
- Desired: ≥ 3 interactions.
- Report a warning when a known hinge-region hydrogen bond is missing.

---

## Upstream / Downstream

### Upstream inputs from structure-prep

```json
{
  "prepared_pdb": "prep_out/4KD1_receptor_prepared.pdb",
  "source_complex_mmcif": "prep_out/4KD1.cif",
  "ligand_sdf": "prep_out/4KD1_ligand_1QK.sdf",
  "reference_ligand_microstates": "prep_out/reference_ligand_microstates.smi",
  "protonation_manifest": "prep_out/protonation_manifest.json",
  "water_report": "prep_out/water_report.json",
  "metal_report": "prep_out/metal_report.json",
  "structure_qc_report": "prep_out/structure_qc_report.json"
}
```

### Downstream use by docking/scoring

| Field | Downstream use |
|------|---------|
| `pocket_definition.json` | Normalized site source, confidence, center, size, residues, and detector evidence |
| `docking_box.json` | Explicit `center[3]` / `size[3]`; consumed downstream only when `search_space_source=docking_box` is declared |
| `interaction_fingerprint` | Pose comparison, interaction screening, and ML features |
| `interaction_observations.json` | Atom/role/geometry-aware profiler evidence and tool-status audit trail |
| `interaction_blueprint.json` | Reviewed hard gates plus desired/optional optimization objectives |
| `arpeggio_interactions.json` | Optional Arpeggio atom contacts plus mmCIF/chemistry quality status |
| `luna_interaction_fingerprint.json` | Optional LUNA atom-aware interactions and hashed/count fingerprint |
| `oddt_interaction_fingerprint.json` | Optional ODDT residue-count IFP and SPLIF comparison evidence |
| `pocket_features` | Descriptive features such as volume, burial proxy, and sub-pockets |
| `pocketminer_predictions.json` | Opt-in residue scores, strict validation provenance, and provisional cryptic hypotheses; retained inside Step 1 and not used as a docking box |
| `deeppocket_predictions.json` | Validated DeepPocket ranks, probabilities, segmented-residue geometry, fpocket seed mapping, hashes, and worker provenance |
| `pocket_analysis_report.json` → `pocket_detection.autosite` | AutoSite native ranks, cluster-point centroids, support-file hashes, input provenance, and command records when selected |

---

## Dependencies

| Package | Step | Install |
|---------|------|---------|
| fpocket (CLI) | Pocket detection | `conda install -c conda-forge fpocket` |
| pyKVFinder | Independent geometric pocket detector | `pip install pyKVFinder` |
| P2Rank (CLI) | ML pocket scoring | Download from GitHub releases |
| PocketMiner (isolated) | Residue-level cryptic-opening evidence | Dedicated `pocketminer` Conda env; see `references/pocketminer.md` |
| DeepPocket (isolated) | fpocket proposal ranking and 3D CNN pocket segmentation | Dedicated `deeppocket` Conda env; see `references/deeppocket.md` |
| AutoSite (ADFRsuite) | Native cluster detection and exported support points | Registered `adfrsuite` runtime with `prepare_receptor` and `autosite`; configure existing installation paths through the tool registry |
| plip | Static interaction profiler | Pin in `drugmate`; invoked in an isolated subprocess |
| prolif + MDAnalysis | Interaction fingerprint | `pip install prolif mdanalysis` |
| pdbe-arpeggio | Supplemental atom contacts | Installed in `drugmate`; CLI worker requires the source complex mmCIF |
| LUNA 0.14 | Supplemental atom-aware interactions / hashed IFP | Dedicated `luna` Conda env; invoke the core Python API headlessly |
| ODDT 0.7 | Legacy residue IFP / SPLIF | Dedicated `oddt` Conda env with Python 3.9, NumPy 1.26, and OpenBabel backend |
| rdkit | Ligand SDF loading for ProLIF | `pip install rdkit` |
| biopython | PDB parsing | `pip install biopython` |
| scipy | Volume calculation | `pip install scipy` |
| numpy | Geometry | `pip install numpy` |

`ligand_site` is built-in coordinate analysis and requires no additional software. For unknown sites, `detector_consensus` runs the exact subset selected by `plan-pocket-tools` from the five-detector repertoire, including the mandatory DeepPocket anchor. The unavailable/off-memory plan uses DeepPocket, fpocket, P2Rank, and pyKVFinder. AutoSite is disabled in that fallback.

PocketMiner is integrated as supplemental evidence, disabled by default and installed in an isolated environment; see `references/pocketminer.md`. Enable it explicitly only after the strict runtime probe passes; successful installation or an upstream script exit code of 0 does not establish that the checkpoint is usable.

DeepPocket is disabled by default in the registry; `detector_consensus` must
explicitly enable it and invoke its configured installation through the strict worker
in the `deeppocket` environment. Classification probabilities, ranked types,
segmentation DX, and native proposal support must agree before candidates are
published; exit code 0 alone is not success. Preserve native classifier ranks
when segmentation returns no residues: record the empty segmented set, retain
its proposal cloud for point5 fusion, and skip that empty set in explicit-site
consensus. DeepPocket depends on fpocket
proposals, so it does not add independent detector support on top of fpocket.

Arpeggio can coexist with the core tools in `drugmate`. Registry-owned workers
run LUNA and ODDT across Conda environments in `luna` and `oddt`, respectively,
without exposing absolute environment paths in Agent configs or mixing their
dependencies into `drugmate`. LigPlot+ is not installed, registered, or integrated
into the pipeline; do not use that tool ID in a config.

---

## Usage

For memory-guided unknown-site analysis, first run `plan-pocket-tools` with the
current noH receptor and bound `pocket_prediction` encoding as shown in the
[planning example](../pipeline-entry/references/pocket-memory.md#integration-with-main-step-1).
Copy its memory input paths/mode into Step 1 and enable exactly `execution_tools`;
every other registered tool still needs an explicit switch. Do not use a fixed
four-tool matrix for a memory plan that selects a different subset. The examples
below use chain `A` only as a reviewed example value; select the real target chain.

```bash
# Full analysis with reference ligand
python skills/bio_research/pocket-analysis/run.py \
  --input '{"protein_file":"4KD1_prepared.pdb","ligand_file":"dinaciclib.sdf","site_selection_strategy":"ligand_site","chain":"A","tools":{"ligand_site":{"enabled":true},"fpocket":{"enabled":true},"pykvfinder":{"enabled":true},"p2rank":{"enabled":true},"deeppocket":{"enabled":false},"autosite":{"enabled":false},"pocketminer":{"enabled":false},"plip":{"enabled":true},"prolif":{"enabled":true},"arpeggio":{"enabled":false},"luna":{"enabled":false},"oddt":{"enabled":false},"blueprint":{"enabled":true}}}' \
  --output-dir ./pocket_out/

# Unknown-site fusion with an explicit memory-off ablation
python skills/bio_research/pocket-analysis/run.py \
  --input '{"protein_file":"receptor.pdb","site_selection_strategy":"detector_consensus","chain":"A","pocket_memory_mode":"off","tools":{"ligand_site":{"enabled":false},"fpocket":{"enabled":true},"pykvfinder":{"enabled":true},"p2rank":{"enabled":true},"deeppocket":{"enabled":true,"params":{"max_candidates":3}},"autosite":{"enabled":false},"pocketminer":{"enabled":false},"plip":{"enabled":false},"prolif":{"enabled":false},"arpeggio":{"enabled":false},"luna":{"enabled":false},"oddt":{"enabled":false},"blueprint":{"enabled":false}}}' \
  --output-dir ./pocket_out/

# With validation against known site
python skills/bio_research/pocket-analysis/run.py \
  --input '{"protein_file":"4KD1_prepared.pdb","site_selection_strategy":"known_site","chain":"A","known_site_residues":["Ile10","Val18","Ala31","Lys33","Glu81","Leu83"],"tools":{"ligand_site":{"enabled":false},"fpocket":{"enabled":true},"pykvfinder":{"enabled":true},"p2rank":{"enabled":true},"deeppocket":{"enabled":false},"autosite":{"enabled":false},"pocketminer":{"enabled":false},"plip":{"enabled":false},"prolif":{"enabled":false},"arpeggio":{"enabled":false},"luna":{"enabled":false},"oddt":{"enabled":false},"blueprint":{"enabled":false}}}' \
  --output-dir ./pocket_out/
```

---

## Output Files

| File | Description |
|------|-------------|
| `pocket_analysis_report.json` | Full analysis report with all results |
| `pocket_definition.json` | Normalized selected pocket, source/confidence, chain-safe residues, center/size, and detector evidence |
| `docking_box.json` | Validated `center[3]` and positive `size[3]` contract consumed when docking explicitly selects `search_space_source=docking_box` |
| `interaction_fingerprint.json` | Per-residue interaction bit vector |
| `interaction_observations.json` | Canonical atom/role/geometry observations, profiler statuses, and strict consensus groups |
| `interaction_blueprint.json` | Schema-v2 required/desired/optional constraints and hard-gate validity |
| `plip_interactions.json` | Rich PLIP result when PLIP succeeds |
| `arpeggio_interactions.json` | Optional Arpeggio result, source-mmCIF provenance, and chemistry-quality flags |
| `luna_interaction_fingerprint.json` | Optional LUNA interactions and hashed/count fingerprint |
| `oddt_interaction_fingerprint.json` | Optional ODDT residue-count IFP and SPLIF records |
| `pocket_features.json` | Volume, burial ratio, sub-pockets |
| `fpocket_output/` | Raw fpocket output directory |
| `p2rank_output/` | Raw P2Rank output directory |
| `pocketminer_predictions.json` | Opt-in strict provenance, residue scores, and provisional cryptic hypotheses; never a standalone docking box |
| `deeppocket_predictions.json` | Normalized DeepPocket detector result with strict provenance |
| `deeppocket_output/` | Raw DeepPocket cleaned receptor, fpocket seeds, classifier ranks/probabilities, masks, and segmented residue PDBs |
| `autosite_output/` | Selected AutoSite run's private receptor/preparation, native summary CSV, cluster-point PDBs, and stdout/stderr logs |

---

## References

Detailed parameters, code examples, and caveats for each tool:

- `references/fpocket.md` — fpocket parameters, output files, and pitfalls
- `references/p2rank.md` — P2Rank usage and CSV parsing
- `references/prolif.md` — ProLIF fingerprint calculation and interaction types
- `references/pykvfinder.md` — pyKVFinder Python cavity detection
- `references/pocketminer.md` — PocketMiner's isolated environment, strict checkpoint probe, and pipeline integration semantics
- `references/deeppocket.md` — DeepPocket's isolated environment, actual output parsing, dependency-based evidence grouping, and failure semantics
- [Memory and AutoSite integration](../pipeline-entry/references/pocket-memory.md) — detector planning, native support-point provenance, fusion equations, and selection rules
- `references/interaction-profilers.md` — environment, capability, and routing matrix for PLIP/ProLIF/Arpeggio/LUNA/ODDT, and the boundary excluding LigPlot+

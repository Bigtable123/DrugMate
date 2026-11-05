# Pocket Selection: Protein Memory and Local Fusion

For unknown sites, main pipeline Step 1 uses
`site_selection_strategy=detector_consensus`. Its default
`pocket_memory_mode=auto` uses `protein_memory_local_fusion_v1` when compatible
query encoding and historical pocket evidence are available. Only unavailable
memory, or an explicit memory-off ablation, selects
`deeppocket_top1_local_rank_fusion_v1`. Record the actual strategy, memory status,
and fallback reason; a missing query is not a successful retrieval with no hits.
`required` makes unavailable memory an error. Corrupt evidence, inconsistent
hashes and invalid current detector results fail instead of silently falling back.

The detector repertoire is DeepPocket, fpocket, P2Rank, pyKVFinder, and AutoSite. Retrieval
happens **before detection**: each of five neighbors recommends its best
detector, producing a five-entry multiset. Run each distinct recommended tool
once and always run DeepPocket to supply the default anchor. AutoSite uses
`pocket-analysis/scripts/autosite_utils.py` in the registered `adfrsuite` runtime;
it preserves native CSV cluster ranks and exported cluster-point geometry.
Selection and normalization rules are recorded in `inference_notes`.
Method checks do not establish experimentally measured performance for a new query.

User-specified fixed boxes, valid ligand-defined sites, and reviewed known sites
still use `provided_box`, `ligand_site`, and `known_site`, respectively, without
being overridden by predictions. Memory fusion applies only to unknown sites.

## Memory-guided algorithm

1. Use the current prepared receptor's compatible SaProt protein embedding.
   Retrieve one nearest case per protein group, then the five most similar groups,
   using cosine similarity. Exclude the query case and known group before
   retrieval and prior estimation. No ligand embedding or current-query ground
   truth is used for this selection. Eligible cases must have raw native-top-1
   DCA for **all five** detectors. Break equal cosine similarities by case
   ID. Each neighbor's `best_detector` is the detector with the smallest raw
   native-top-1 DCA, with equal DCA broken by lexicographic tool ID. Preserve all
   five recommendations, including repeats; execute their distinct tool IDs
   plus DeepPocket. P2Rank is not forced when it is absent from the multiset.
2. For neighbors ranked `j=1..5`, compute
   `raw_j=exp((s_j-s_1)/0.02)/log2(j+1)` and
   `w_j=5*raw_j/sum(raw)`. Weights sum to **5**, not 1.
3. Use raw historical native-top-1 DCA to define
   `y_jt=1[DCA_jt < 4 Å]`. The legacy `top1_dca_le_4` boolean is insufficient
   at the boundary. Let `pi_t` be the eligible historical mean label and compute
   `p_hat_t=(5*pi_t + sum(w_j*y_jt))/(5+sum(w_j))`.
4. Start with DeepPocket#1. If P2Rank was selected and executed, switch to
   P2Rank#1 only if **all** conditions hold:
   DeepPocket native probability `<0.95`; P2Rank probability `>=0.80`;
   their center distance `>5 Å`; `p_hat_P2Rank-p_hat_DeepPocket>0.15`;
   highest retrieved cosine similarity `>=0.88`; and the P2Rank-centered matching
   set has at least two independent detector families, including P2Rank.
5. Include the anchor and **every** selected detector's native Top-3 candidate
   satisfying center distance `<=5 Å` **and** point5 residue IoU `>=0.30` relative
   to the anchor. This includes matching rank-2/3 pockets from the anchor's own
   detector. Count the anchor only once; members need not match one another.
6. For each member pocket, set `v_t=(p_hat_detector+0.05)^2/native_rank_t`.
   Families are `{DeepPocket, fpocket}`, `{P2Rank}`, `{pyKVFinder}`, and `{AutoSite}`.
   Set `V_f=max(v_t in family)` and
   `effective_t=v_t*V_f/sum(v in family)`. Normalize effective weights to compute
   the center; preserve exactly the selected anchor's point5 residue set.

Apply these selection rules consistently: derive
`best_detector` from minimum raw DCA with the stable tie rule above; use complete
five-detector cases for both retrieval and priors; execute repeated recommended
tools once; force DeepPocket when not recommended; and retain all matching Top-3
pockets instead of choosing one per tool. Priors use that same complete population
after exclusions. Missing measurements are never converted into successes or
failures; fewer than five eligible groups makes memory unavailable. An older
derived index without `best_detector` uses the same documented derivation during
reading; a stored winner inconsistent with that rule is rejected. The compatibility
presets `conservative` and `exploratory_memory_guard` both use this policy.

## Fallback when memory is unavailable or disabled

1. Run DeepPocket, fpocket, P2Rank, and pyKVFinder, taking each tool's native
   Top-3 and normalizing them into point5 residue sets on the current
   hydrogen-free receptor. Preserve native ranks; do not relabel rank 2 as rank 1.
2. In this fallback, DeepPocket#1 is always the anchor. A candidate from another tool must satisfy
   both center distance `≤5 Å` and point5 residue IoU `≥0.30` relative to the
   anchor. Choose at most one candidate per tool, ordered by highest IoU, then
   smallest distance, then best native rank. Candidates need to match the anchor,
   not one another.
3. Each member's weight is `v_t=1/rank_t`. Families are `{DeepPocket, fpocket}`,
   `{P2Rank}`, and `{pyKVFinder}`. A family's total weight is its members'
   `max(v_t)`, allocated among those members in proportion to `v_t`. Normalize
   the effective weights and compute the fused center. Do not weight by historical
   reliability or native probability.
4. If no other pocket matches, use the DeepPocket#1 center directly. Output
   residues must always equal exactly the DeepPocket#1 point5 set, preserving
   `(chain, resnum, icode, resname)`. Do not take a union or recompute residues
   around the fused center. PocketMiner supplies annotations only and must not
   rerank candidates.

All four tools must be enabled and each must have a valid native rank-1 candidate.
If a tool naturally returns fewer than three candidates, retain the actual count.
Missing tools, failures, or invalid candidates must be reported and stop this
workflow, without silently falling back to a different selection algorithm.
Top-K and matching-threshold overrides cannot change the registered policy.

## Prepare the pocket index and current query

The supplied docking `ligand_protein_v2.json` and generation index do not contain
the pocket DCA labels needed here. Build a separate index from the downloaded
docking/pocket bundle; preserve the frozen bundle and write the result outside
the checkout:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" build-pocket-index \
  --benchmark-bundle "${SBDD_MEMORY_BUNDLES_ROOT}/drugmate-docking-memory" \
  --output "${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json"
```

The builder reads only declared `pocket_prediction` memory cases, raw published
native-top-1 DCA and compatible bound protein features. It stores all five DCA
measurements and the derived `best_detector`, including its derivation policy.
It does not run detectors or evaluate the current query. Real-bundle rebuilding
has been checked: 407 cases and 267 groups; 404 cases / 265 groups have complete
measurements. The three incomplete cases retain missing labels rather than
fabricated failures. This checks the offline data adapter, not experimental results.

Prepare a current-input-only CSV and use the [SaProt encoder](../../../../memory/saprot_embedding/README.md)
with `--task-type pocket_prediction`, matching the index's encoder recipe.
Bind the finished encoding with `benchmark_memory.py bind-query-metadata`;
retain the original encoder metadata and completion evidence. Do not relabel a
docking query as a pocket query or change hashes to force compatibility.

## Integration with main Step 1

Use the current hydrogen-free receptor, assembly, author chain, and coordinate
frame reviewed under the [preparation contract](preparation-contract.md). Prepare
memory assets first, then inspect the plan before enabling any detectors:

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" plan-pocket-tools \
  --memory-mode required \
  --selector-index "${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json" \
  --query-npz "${SBDD_DATA_ROOT}/REQUEST/inputs/embeddings/pocket_prediction/QUERY/TARGET/embeddings.npz" \
  --query-metadata "${SBDD_DATA_ROOT}/REQUEST/inputs/query_metadata.json" \
  --receptor-file "${SBDD_DATA_ROOT}/REQUEST/inputs/receptor.pdb" \
  > "${SBDD_DATA_ROOT}/REQUEST/inputs/pocket_tool_plan.json"
```

Use `auto` instead of `required` only when an explicitly reported unavailable-memory
fallback is acceptable. Planning runs no detectors and launches no encoder. In the
main config:

- Set `steps.pocket_analysis.enabled=true` and
  `steps.pocket_analysis.params.site_selection_strategy="detector_consensus"`.
- Include all five detector IDs in the complete tool matrix: `deeppocket`,
  `fpocket`, `p2rank`, `pykvfinder`, and `autosite`. Enable **exactly** the IDs in
  the plan's `execution_tools`; explicitly disable the other detectors. For
  example, recommendations `[p2rank, p2rank, autosite, autosite, p2rank]` require
  DeepPocket, P2Rank, and AutoSite; fpocket and pyKVFinder stay disabled.
  DeepPocket's internal fpocket proposal generation still runs. Its
  `max_candidates` is fixed at 3. Decide every other registered tool explicitly.
- Set `pocket_memory_mode="auto"` (or `"required"` for an experiment that must
  use memory). Supply `pocket_memory_index`, `pocket_memory_query_npz`, and
  `pocket_memory_query_metadata` under `steps.pocket_analysis.params`, using
  actual absolute paths. Supply known case/group identities with
  `pocket_memory_query_case` / `pocket_memory_query_group`; additional exclusions
  use `pocket_memory_exclude_cases` / `pocket_memory_exclude_groups`.
- Downstream stages use Step 1's published `pocket_definition`, `pocket_center`,
  and `docking_box`. Docking sets `search_space_source="docking_box"`.

Step 1 repeats planning against the same receptor and memory inputs, checks that
the enabled detector set equals `execution_tools`, and performs detection,
point5 normalization, and fusion in the same run. A mismatch fails before
detector execution; review and update the explicit matrix if the plan changes.
Before publishing the fused result, Step 1 also checks that the receptor,
memory index, query provenance and recommendation policy still match its
pre-detection plan. Changes during detector execution fail rather than mixing
two memory versions in one result.
Fallback planning returns DeepPocket, fpocket, P2Rank, and pyKVFinder, with
AutoSite disabled. A selected tool failing or being absent from its registered
runtime is an execution error, not a reason to ignore usable memory. Do not feed a preselected result
back into the main workflow through `provided_box`: that strategy regenerates
residues from heavy atoms inside the box, violating the anchor point5 contract.

DeepPocket still follows the real classifier/segmentation, checkpoint, and
runtime-environment contracts. Include every planned tool's computations in the main
experiment's budget, authorization, and lineage. Temporary files obey this
request's output-root constraint as well. Existing predictions may be independently
inspected through the entry points below without rerunning scientific tools.

## point5 support points and independent inspection

point5 contains residues whose **current receptor protein heavy atoms are within
5 Å of the tool's original support points**. This 5 Å definition differs from the
center-matching distance. A native tool residue list cannot simply be renamed
point5.

| Tool | Support points and native rank |
|---|---|
| fpocket | Alpha spheres in `pocket<rank>_vert.pqr` for the native rank |
| DeepPocket | Proposal PQR for the classifier rank's `original_fpocket_rank`; do not substitute segmentation PDB/DX for support points |
| P2Rank | `surf_atom_ids` in the original predictions CSV, mapped to atom serials in the same input PDB; do not substitute visualization labels |
| pyKVFinder | Grid points in the output cavity PDB whose resname equals `cavity_id`; retain the existing rank in descending cavity-volume order |
| AutoSite | Exported points in `receptor_cl_<rank>.pdb`, with cluster rank from `receptor_summary.csv`; center is the point centroid, not a cavity-volume estimate |

If original predictions for the current receptor already exist, create
`candidate_sources.json` with the original detector-input SHA-256 and actual
paths. Relative paths are resolved against this manifest. Verify `<receptor_stem>`
and DeepPocket's cleaned `<nowat_stem>`; matching directory names do not establish
matching provenance:

```json
{
  "target": "QUERY",
  "receptor_sha256": "SHA256_RECORDED_FOR_DETECTOR_INPUT",
  "step1_report": "pocket_evidence/pocket_analysis_report.json",
  "sources": {
    "fpocket_pockets": "pocket_evidence/fpocket_output/<receptor_stem>_out/pockets",
    "deeppocket_pockets": "pocket_evidence/deeppocket_output/<nowat_stem>_out/pockets",
    "p2rank_predictions": "pocket_evidence/p2rank_output/<receptor_filename>_predictions.csv",
    "pykvfinder_points": "pocket_evidence/pykvfinder_output/<receptor_stem>_cavities.pdb"
  }
}
```

The example shows all four fallback detector sources. For a memory-guided run,
include source entries only for tools actually selected. AutoSite pockets in the
Step 1 report carry their own `support_file` paths and native ranks. A native
`detectors.autosite` array may supply those same fields for independent inspection.

```bash
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" normalize-pocket-candidates \
  --candidates "${SBDD_DATA_ROOT}/REQUEST/inputs/candidate_sources.json" \
  --receptor-file "${SBDD_DATA_ROOT}/REQUEST/inputs/receptor.pdb" \
  --output "${SBDD_DATA_ROOT}/REQUEST/inputs/pocket_candidates.json"
"${SBDD_PYTHON}" "${BIO_SKILLS_ROOT}/pipeline-entry/scripts/benchmark_memory.py" select-pocket \
  --memory-mode required \
  --selector-index "${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json" \
  --query-npz "${SBDD_DATA_ROOT}/REQUEST/inputs/embeddings/pocket_prediction/QUERY/TARGET/embeddings.npz" \
  --query-metadata "${SBDD_DATA_ROOT}/REQUEST/inputs/query_metadata.json" \
  --receptor-file "${SBDD_DATA_ROOT}/REQUEST/inputs/receptor.pdb" \
  --candidates "${SBDD_DATA_ROOT}/REQUEST/inputs/pocket_candidates.json" \
  --export-box "${SBDD_DATA_ROOT}/REQUEST/inputs/fused_docking_box.json" \
  > "${SBDD_DATA_ROOT}/REQUEST/inputs/pocket_fusion_decision.json"
```

For a deliberate no-memory ablation, use `--memory-mode off` without memory or
encoding paths. The adapter
also accepts native `detectors.<tool>` arrays. Each item supplies `rank`, `center`,
and one of `support_points` (an xyz array), `support_file` (PDB/PQR), or
`support_atom_serials` (an array of serials from the current input PDB). A shared
cavity file also supplies `support_resname`. DeepPocket/P2Rank must still include
native `probability` to verify detector artifacts and evaluate the anchor guard.
Probabilities are not fusion weights. The adapter emits the point5
definition, input/support-file hashes, and support-point counts without modifying
original detector files. When the receptor changes, regenerate matching
candidates; editing recorded hashes does not make old predictions reusable.

Candidate tool IDs must equal the plan's `execution_tools`; extra or missing
detectors do not silently alter the selected policy. AutoSite is currently
supported only by `detector_consensus`; disable it for explicit-site strategies.

Outputs retain the policy name, implementation/inference notes, five-entry
`recommended_multiset`, `execution_tools`, forced-anchor status, anchor, matching
members, native ranks, family weights, final center, anchor residues, memory provenance, retrieved groups,
reliabilities, guard conditions and any fallback reason for report and artifact-lineage
inspection. `--export-box` applies an explicit downstream geometry rule: retain
the fused center `c`, and use **current receptor heavy atoms of the anchor point5
residues** to calculate each axis as
`size_d=max(min_size, 2 × (max_atom |x_d-c_d| + padding))`, with defaults of 4 Å
padding on each side and a minimum edge length of 12 Å. This reads neither an
experimental ligand nor a reference box and does not establish experimental
validation of the site. Standalone commands only inspect/export the specified
artifacts; they do not edit the config or launch the pipeline.

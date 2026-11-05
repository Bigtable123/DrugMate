# DeepPocket pipeline integration

DeepPocket is required for Step 1 `detector_consensus` and optional for explicit
site strategies. It runs in the dedicated
`deeppocket` Conda environment against `${SBDD_DEEPPOCKET_SOURCE_DIR}`; its PyTorch,
libmolgrid, ProDy, and fpocket dependencies are never imported into the Drugmate
process.

## Scientific role

The upstream workflow first generates candidates with fpocket, ranks those
candidates with the DeepPocket classifier, and segments the requested Top-K
pockets with its 3D CNN. Consequently, DeepPocket and a separately enabled
fpocket run are assigned to the same `fpocket_family` independence group. They
remain separate provenance-bearing results, but agreement between them does not
count as two independent detector votes or satisfy the two-detector interaction
consensus by itself.

The normalized candidate center and classifier probability come from
`bary_centers_ranked.types` and `bary_centers_confidence.txt`. Pocket residues
in the native detector report come from the segmented `*_pocketN.pdb` file.
For `detector_consensus`, point5 residues are computed from the original proposal
`pocket<original_fpocket_rank>_vert.pqr` support points; the final residues are
the selected anchor's point5 set. The anchor is DeepPocket#1 by default, or
P2Rank#1 when all six memory-guided switch conditions are satisfied. The bridge
retains native classifier ranks even if a segmentation residue set is empty,
records that empty result, and never substitutes proposal residues for CNN
segmentation. Proposal-based point5 fusion can still use that native candidate;
explicit-site consensus skips its empty segmented-residue set. A volume,
when present, is labeled `deeppocket_fpocket_seed_volume`; it is not represented
as a CNN-mask volume.

## Adapter parameters

Whether DeepPocket is enabled is decided in
`../../pipeline-entry/TOOL_SELECTION_POLICY.md` and written explicitly in the
complete config matrix. When it is selected, its adapter parameter shape is:

```json
{
  "tools": {
    "deeppocket": {
      "enabled": true,
      "params": {
        "max_candidates": 3,
        "segmentation_threshold": 0.5,
        "mask_distance": 3.5,
        "timeout_seconds": 1200
      }
    }
  }
}
```

`max_candidates` must be 1–10, `segmentation_threshold` must be strictly
between 0 and 1, and `mask_distance` is in Å. Runtime paths and the Conda
environment belong to the central registry, not scientific run configuration.

## Strict output contract

The isolated bridge validates all of the following before returning success:

- repository commit and both checkpoint SHA-256 digests;
- importable torch/libmolgrid/Biopython/ProDy, strict state-dict/model matches,
  an available fpocket CLI, and an actual CUDA allocation during preflight;
- a prepared PDB containing coordinate records;
- successful real `predict.py` execution;
- non-empty ranked candidates and classifier probabilities with matching
  counts, finite `[0, 1]` values, and rank-sorted order;
- a non-empty DX mask and original proposal support file for every requested rank.

Segmented-residue availability is reported separately and can be false after a
successful native execution. Empty segmentation does not remove a native Top-3
classifier candidate or renumber lower ranks. The point5 adapter validates the
proposal cloud against the current receptor before fusion; explicit-site
consensus requires actual segmented residues and reports no detected pocket if
no enabled detector provides usable geometry.

The normalized artifact is `deeppocket_predictions.json`; raw cleaned receptor,
fpocket proposals, ranked types, confidences, DX masks, and segmented pocket PDBs
are retained under `deeppocket_output/`. The result records input/checkpoint
hashes, repository commit, scientific parameters, subprocess provenance, and
the mapping from DeepPocket ranks back to original fpocket proposal ranks.

If an enabled DeepPocket run fails on a particular structure, Step 1 records the
error code and diagnostic while preserving successful independent geometry
detectors. If no enabled detector supplies usable geometry, the existing pocket
contract fails and canonical `pocket_definition.json` / `docking_box.json` are
not published. Missing environment/import/CLI prerequisites are expected to be
rejected earlier by pipeline runtime preflight.

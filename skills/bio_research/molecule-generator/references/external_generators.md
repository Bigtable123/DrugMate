# External Generator Runtime Reference

This file specifies Step 2 runtime and model requirements. Configure source
paths and environments through the [deployment variables](../../pipeline-entry/references/deployment.md).
Model checksums and tool contracts are declared in `shared/pipeline_registry.py`;
experiment configs do not replace those requirements.

## Runtime requirements

| Tool | Source directory setting | Environment setting or ID | Fixed model asset | SHA-256 | Step 2 role |
|---|---|---|---|---|---|
| ResGen | `${SBDD_RESGEN_SOURCE_DIR}` | `resgen` | `ckpt/val_172.pt` | `45cac2f5186564e71b43fe43ce709ae2515a3a32847ef592d6f99810d32916df` | Fast pocket-conditioned 3D generation; requires a reference ligand |
| DiffSBDD | `${SBDD_DIFFSBDD_SOURCE_DIR}` | `diffsbdd` | `checkpoints/crossdocked_fullatom_cond.ckpt` | `07f86764bf569aafbc40a9c15fc02de8e2550437dd0f17f657eab3abe66c372c` | Pocket-conditioned full-atom 3D diffusion |
| TargetDiff | `${SBDD_TARGETDIFF_SOURCE_DIR}` | `targetdiff-gpu` | `pretrained_models/pretrained_diffusion.pt` | `8977821339ba31b7103b46835ef4ddbf7f15f3f2efc5cb1f2a405e7905b09db6` | Pocket-conditioned 3D diffusion |
| DiffLinker | `${SBDD_DIFFLINKER_SOURCE_DIR}` | `difflinker` | `models/pockets_difflinker_full.ckpt` | `695ef7d634c09db4b9635c4e5efebad23385b8432f0687b6eae5e726ba7883f9` | Pocket-conditioned linking of supplied 3D fragments |
| REINVENT4 | `${SBDD_REINVENT4_SOURCE_DIR}` | `reinvent4` | `priors/reinvent_pubchem.prior` | `fe8cd1678452ad292a8f93e97cb19a85959b729e17113e157180d5e69ae89ef3` | Unconditional ligand-only 2D sampling |
| Lingo3DMol | `${SBDD_LINGO3DMOL_SOURCE_DIR}` | `${SBDD_ENV_LINGO3DMOL_PREFIX}` | `checkpoint/gen_mol.pkl` + `checkpoint/contact.pkl` | `edc665e1b46ce1ff39d03a15c56b9c019ec44a8573acb51d0698719ef87247ab` + `14fd26a256c8979cd4e4ae87c3d1575444fb88913475b002a415a1544dad7990` | Pocket-conditioned native 3D language-model generation |
| TamGen | `${SBDD_TAMGEN_SOURCE_DIR}` | `${SBDD_ENV_TAMGEN_PREFIX}` | `checkpoints/crossdock_pdb_A10/checkpoint_best.pt` | `526c6337132b93050eb2b23cd8177175ad3d7d2310177536a1dca1949ea45b4b` | Target-aware conditional SMILES generation |
| MolCRAFT | `${SBDD_MOLCRAFT_SOURCE_DIR}` | `molcraft` | `checkpoints/last.ckpt` + checkpoint-specific config | `0a54f9bfb847ea87e7aec538df90a009c8c49d12a7b5d09947d04b7f9444e024` | Continuous-space pocket-conditioned native 3D generation |

TargetDiff requires a compatible GPU environment selected through
`SBDD_ENV_TARGETDIFF_PREFIX` or `SBDD_ENV_TARGETDIFF_NAME`.

## Tool-specific contracts

### ResGen

- Inputs: prepared receptor, Step 1 pocket center, and reference-ligand SDF.
- Output: native 3D molecules parsed from ResGen's aggregate SDF.
- The installed upstream CLI requires `--sdf_file` even when a center is
  supplied. Do not invent a placeholder ligand to bypass that contract.
- Preserve NumPy/SciPy compatibility warnings in the run log.

### DiffSBDD

- Inputs: prepared receptor and pocket center; a reference ligand is used when
  available.
- Output: full-atom native 3D SDF records.
- Without a reference ligand, the adapter writes a center-cropped pocket with
  chain `A` and sequential residue numbers because the upstream selector cannot
  parse blank chains or PDB insertion codes.
- `timesteps` may be reduced for a smoke test, but scientific generation should
  use the configured reviewed value.

### TargetDiff

- Inputs: prepared receptor and pocket center.
- The adapter creates a standard-residue pocket crop around the center and
  launches `python -m scripts.sample_for_pocket` with the repository on
  `PYTHONPATH`. The upstream README-style direct script launch is not used
  because it cannot resolve the repository's `utils` package in this install.
- Output: reconstructed native 3D SDF records. A short diffusion run can exit
  successfully with zero reconstructable molecules; zero valid output is still
  a source failure.

### DiffLinker

- Inputs: prepared receptor and `context.fragment_file`.
- The fragment input must contain a single 3D record with at least two
  disconnected components.
- `linker_size` is explicit in this adapter. Optional anchor indices must refer
  to atoms in the supplied fragment record.
- Open Babel/RDKit bond perception can reject part of a generated batch. Filter
  records individually and keep only sanitized, connected, finite 3D molecules.

The registered checkpoint corresponds to the official pocket-conditioned
DiffLinker model asset (official MD5:
`b8edbaf204933bb1d452e5cdd09b3435`).

### REINVENT4

- Inputs: no receptor or pocket artifact.
- Output: validated 2D SMILES; Step 2 creates downstream conformers.
- The installed reviewed asset is the unconditional PubChem REINVENT prior.
  This route is therefore a diversity baseline, not pocket-conditioned
  generation, lead optimization, or reinforcement learning.
- Do not claim Mol2Mol, LibInvent, or LinkInvent behavior unless the
  corresponding reviewed prior and adapter are separately registered.
- The official prior MD5 is
  `75052d5e93b2470e9c6684d369704f5f`. An upstream internal metadata-hash warning
  under the installed Python version is retained in provenance; the registry
  SHA-256 remains the hard local asset check.

### Lingo3DMol

- Inputs: prepared receptor PDB and Step 1 pocket center. The adapter writes a
  center/radius-cropped standard-residue pocket, and the worker copies it to a
  fixed safe basename before calling the upstream inference script.
- Both the contact and generation checkpoints are registry assets.
- Output: upstream MOL files recursively parsed as connected, finite native 3D
  molecules. The requested count is enforced after parsing.

### TamGen

- Inputs: prepared receptor PDB and Step 1 pocket center.
- The worker builds TamGen's coordinate-aware fairseq test dataset directly
  from the prepared PDB, avoiding an RCSB download and preserving the exact
  Step 0 receptor used by the pipeline.
- Output: unique valid SMILES from real `H-*` model hypotheses. These records
  intentionally have `coordinate_origin=smiles_only`; Step 2 supplies ETKDG
  conformers afterward.
- The default does not enable TamGen's ligand-latent `--gen-vae` mode because
  Step 1 supplies no seed molecule for that latent. Enable `use_conditional`
  only for an explicitly reviewed conditional-ligand route.

### MolCRAFT

- Inputs: prepared receptor and reference ligand, which the reviewed local
  entrypoint uses to define the pocket.
- Use only `sample_for_pocket_local.py` with its checkpoint-specific config.
- The worker redirects MolCRAFT's process-global output path to the fresh
  invocation so the upstream source checkout is never deleted or reused.
- Output: connected, finite native 3D SDF molecules.

## Generation-only boundary

None of these adapters invokes Vina, Meeko, PDBQT preparation, rescoring, or
optimization. Upstream repositories may document those optional workflows, but
this Step 2 integration only creates candidate molecules. Docking is owned by
the pipeline docking step.

Every call receives a fresh `<tool>_work/run-*` directory. Parsers inspect only
that invocation's files, so a retry in the same Step 2 output directory cannot
mix previous SDF, MOL, or CSV records into the new candidate set.

All native 3D outputs remain hypotheses. If the pipeline continues through
docking, both PoseBusters and PoseCheck are mandatory before nomination.

## Upstream projects and model records

- DiffLinker: <https://github.com/igashov/DiffLinker>
- DiffLinker pocket model record: <https://zenodo.org/records/10988017>
- DiffSBDD: <https://github.com/arneschneuing/DiffSBDD>
- TargetDiff: <https://github.com/guanjq/targetdiff>
- ResGen: <https://github.com/OdinZhang/ResGen>
- REINVENT4: <https://github.com/MolecularAI/REINVENT4>
- REINVENT4 priors: <https://doi.org/10.5281/zenodo.15641296>
- Lingo3DMol: <https://github.com/stonewiseAIDrugDesign/Lingo3DMol>
- TamGen: <https://github.com/microsoft/TamGen>
- MolCRAFT: <https://github.com/GenSI-THUAIR/MolCRAFT>

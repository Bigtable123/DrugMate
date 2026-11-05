# Scientific tools

Install only the tools selected for your workflow. Start with the two installation
steps in the [README](../README.md); a complete SBDD run needs the dependencies
of every enabled stage. This list describes the adapters in this source release,
not software already installed by DrugMATE.

Environment IDs below are registry IDs, not package names. `drugmate` means the
**scientific runtime**, separate from the Agent environment. Point each ID at its
actual Conda environment using `SBDD_ENV_<ID>_PREFIX` or `_NAME`; see the
[deployment reference](../skills/bio_research/pipeline-entry/references/deployment.md).
Model-based tools also need their compatible source trees and weights. Use the
[machine-readable catalog](../skills/bio_research/pipeline-entry/tool_catalog.json)
for exact entrypoints, checkpoints, dependencies and parameters, and the
[selection policy](../skills/bio_research/pipeline-entry/TOOL_SELECTION_POLICY.md)
for supported uses.

## Structure and ligand preparation

| Software | Tool ID / use | Environment and additional requirements |
| --- | --- | --- |
| RDKit | `ligand_sdf`; molecular preparation and analysis | `drugmate`; also required by several dedicated tool environments |
| Open Babel | `openbabel`; input format and 3D conversion | `drugmate`; executable and Python bindings. The input-normalization compatibility fallback additionally needs C++17, matching Open Babel development files and Eigen headers |
| PDBFixer | `pdbfixer` | `drugmate`; supported heavy-atom repair |
| PROPKA | `propka` | `drugmate`; diagnostic pKa/protonation reports |
| Reduce2 / CCTBX | `reduce` | `reduce`; `mmtbx.reduce2`, compatible monomer library and CCD chemical-components data |
| PDB2PQR + PROPKA | `pdb2pqr` | `drugmate`; `pdb2pqr30` final receptor preparation route |
| Meeko | `meeko` | `reduce`; receptor preparation plus RDKit/Meeko ligand preparation and SDF export for Vina |
| Dimorphite-DL | `ligand_microstates` | `drugmate`; reference-ligand protonation states |

The two supported local final-receptor backends are Reduce2 and PDB2PQR; select
one explicitly. Protoss is listed under network services below because its
current adapter is unavailable for this role.

## Pockets and interaction analysis

| Software | Tool ID | Environment and additional requirements |
| --- | --- | --- |
| fpocket | `fpocket` | `drugmate`; `fpocket` executable |
| pyKVFinder | `pykvfinder` | `drugmate`; `pyKVFinder` Python module |
| P2Rank | `p2rank` | `drugmate`; working `prank` executable and its runtime |
| DeepPocket | `deeppocket` | `deeppocket`; source, classification/segmentation checkpoints, fpocket and CUDA |
| AutoSite | `autosite` | `adfrsuite`; `prepare_receptor` and `autosite` commands; native cluster CSV/PDB outputs |
| PocketMiner | `pocketminer` | `pocketminer`; source and checkpoint; cryptic-opening annotation, not a standalone pocket definition |
| PLIP | `plip` | `drugmate`; primary interaction profiling |
| ProLIF | `prolif` | `drugmate`; RDKit/MDAnalysis dependencies; supporting fingerprints |
| PDBe Arpeggio | `arpeggio` | `drugmate`; `pdbe-arpeggio`, Gemmi and Open Babel |
| LUNA | `luna` | `luna`; dedicated compatible Python environment |
| ODDT | `oddt` | `oddt`; compatible Open Babel and NumPy versions |

PLIP, ProLIF, Arpeggio, LUNA and ODDT are used in both pocket and pose analysis.
For unknown-site `detector_consensus`, run `plan-pocket-tools` first and install
every returned `execution_tools` dependency. Memory chooses from the five cavity
detectors and retains DeepPocket as the default anchor. Only the unavailable/off
memory fallback requires DeepPocket, fpocket, P2Rank and pyKVFinder together.

## Molecule generation

| Software | Tool ID / environment ID | Additional requirements |
| --- | --- | --- |
| Pocket2Mol | `pocket2mol` | Source, checkpoint, PyTorch-Geometric and CUDA |
| ResGen | `resgen` | Source, checkpoint and compatible graph-learning/CUDA environment |
| DiffSBDD | `diffsbdd` | Source, checkpoint and compatible PyTorch/CUDA environment |
| TargetDiff | `targetdiff` | Source, checkpoint and compatible graph-learning/CUDA environment |
| DiffLinker | `difflinker` | Source/checkpoint and actual disconnected fragments with linking anchors |
| REINVENT4 | `reinvent4` | `reinvent` entrypoint and reviewed sampling prior; ligand-only generation |
| Lingo3DMol | `lingo3dmol` | Source and model assets; retains generated 3D coordinates |
| TamGen | `tamgen` | Source/model assets; pocket-conditioned SMILES generation |
| MolCRAFT | `molcraft` | Source/checkpoint and compatible model environment; retains generated 3D coordinates |

ChEMBL sourcing and DUD-E decoys use data services rather than a local generation
model. Local standardization, conformers and BRICS decoys use RDKit.

## Docking engines

| Software | Tool ID | Environment and additional requirements |
| --- | --- | --- |
| GNINA | `gnina` | `drugmate`; working GNINA executable with its scoring assets/runtime |
| Uni-Dock | `unidock` | `unidock`; `unidocktools` topology pipeline and GPU runtime |
| AutoDock Vina | `autodock_vina` | `vina` for the Python `vina.Vina` API; `reduce` for RDKit/Meeko preparation and export |
| DiffDock | `diffdock` | `diffdock`; source, models and compatible inference environment |
| DOCK6 | `dock6` | `dock6`; executable plus prepared receptor-site and grid artifacts |
| AutoDock4 | `autodock4` | `autodock4`; AutoDockTools, AutoGrid4 and AutoDock4 |
| rDock | `rdock` | `rdock`; `rbcavity`, `rbdock` and cavity configuration |
| ADFRsuite | `adfr` | `adfrsuite`; preparation commands, AGFR and ADFR; current route is non-covalent |

Choose one primary engine per docking configuration. Installing an engine alone
does not supply receptor preparation, ligand states or pose assessment.

## Pose assessment, scoring and simulation

| Software | Tool ID | Environment and scope |
| --- | --- | --- |
| PoseBusters | `posebusters` | `drugmate`; geometric validity checks |
| PoseCheck | `posecheck` | `drugmate`; strain and clash assessment |
| RTMScore | `rtmscore` | `rtmscore`; source, graph dependencies and three ensemble checkpoints |
| AmberTools | `candidate_topology` | `amber`; pdb4amber, antechamber, parmchk2, tleap and sqm |
| OpenMM | `openmm_md` | `openmm`; current short-MD adapter requires CUDA and prepared candidate topologies |
| MDAnalysis | `mdanalysis` | `openmm`; trajectory analysis and its supporting libraries |
| AmberTools MM-GBSA | `mmgbsa` | `mmgbsa`; MMPBSA.py, cpptraj, mmpbsa_py_energy and ParmEd; current adapter gives single-frame estimates |
| OpenFE | `openfe_rbfe` | `openfe`; OpenFE/GUFE, OpenFF, OpenMM/OpenMMTools and RDKit; current adapter performs identity checks and mapping/planning |
| TIES | `ties_rbfe` | `ties`; TIES builder, TIES_MD APIs/commands, compatible OpenMM and AmberTools; current adapter produces a pair/campaign plan |

OpenFE/TIES registration does not mean these adapters compute a completed RBFE
result. Molecular simulation and rescoring are optional escalations, not minimum
Agent or input-normalization requirements.

## Property filtering and synthesis

| Software | Tool ID | Environment and additional requirements |
| --- | --- | --- |
| RDKit | `lipinski`, `qed`, `sa_score`, `pains`, `dedup` | `drugmate`; built-in chemistry analyses |
| RAscore | `rascore` | `rascore`; XGBoost model assets |
| ADMET-AI | `admet` | `admet_ai`; compatible models; adapter uses CPU inference |
| ChemFH | `chemfh` | `chemfh`; matching source revision, checkpoint and native rule bundle |
| AiZynthFinder | `aizynthfinder` | `aizynth`; configuration, expansion/filter policies and stock database |

These tool IDs appear in both property stages; their allowed settings differ
between pre-docking screening and final candidate assessment.

## Data services, reporting and memory support

| Component | Tool ID / use | Installation or availability |
| --- | --- | --- |
| RCSB PDB / CCD | `intake` and reference structures | Network data service; client code is included |
| ChEMBL | `chembl` | Network bioactivity service; client code is included |
| DUD-E | `dude_decoys` | Downloaded/cached benchmark decoys when available |
| ProteinsPlus Protoss | `protoss` | HTTP service; **currently unavailable** as a final receptor backend because prepared-complex identity preservation is not implemented |
| RDKit depictions | `rdkit_2d` | RDKit in the reporting runtime |
| matplotlib | `plots` | Plotting library in the reporting runtime |
| PyMOL | `pymol_cover` | Optional static-cover renderer |
| 3Dmol.js / Mol* | `mol_3d`, `molstar_main` | Viewer assets included for HTML reports; no web server installation |
| SaProt / Foldseek | Query encoding, not pipeline tool switches | Separate compatible environment, local model weights and Foldseek; see the [encoder guide](../memory/saprot_embedding/README.md) |
| Benchmark memories | Historical selector/reader inputs | Supplied small indexes and optional large downloads; see [memory deployment](memory-and-paths.md) |

## Included analyses that do not need another tool installation

These registered switches use project code and the scientific libraries above:

| Stage | Included tool IDs |
| --- | --- |
| Structure preparation | `metal_validation`, `structure_qc`, `residue_validation` |
| Pocket analysis | `ligand_site`, `blueprint` |
| Molecule generation | `rdkit_decoys`, `standardization`, `conformer` |
| Docking preparation / checks | `docking_prep`, `water_selection`, `redocking_validation` (currently GNINA-specific) |
| Pose assessment | `identity_qc`, `blueprint_match`, `pocket_occupancy`, `rdkit_pose`, `freesasa` (RDKit's `rdFreeSASA`), `rescore_rtmscore` |
| Coordination | `require_quality_gates`, `screening_benchmark` |


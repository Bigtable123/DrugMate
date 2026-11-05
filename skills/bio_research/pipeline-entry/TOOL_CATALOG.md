# SBDD Tool Catalog (generated)

> Generated from `shared/pipeline_registry.py` by `pipeline-entry/build_tool_catalog.py`. **Do not edit by hand** — re-run the generator after changing the registry.

This is the inventory the Agent uses while building a pipeline config. For tool choice see the generated `TOOL_SELECTION_POLICY.md`; for assembling confirmed facts and choices into config see `QUERY_TO_CONFIG.md`; for the execution contract see `SKILL.md`.

`${SBDD_*}` and `${BIO_SKILLS_ROOT}` are deployment templates in this catalog. Runtime defaults resolve them from the deployment settings. Explicit JSON path overrides must contain actual paths; JSON does not expand shell variables.

## Top-level config keys

`context`, `input_artifacts`, `project`, `runtime_resources`, `start_step`, `step_timeouts`, `steps`, `stop_step` — anything else hard-fails.

## Tool runtime environments

Runtime placement is registry-owned; Agent configs select tools but do not repeat Conda paths.

| id | kind | location |
|---|---|---|
| `current` | current | `current process` |
| `drugmate` | conda | `hermes` |
| `admet_ai` | conda | `admet_ai` |
| `chemfh` | conda | `chemfh` |
| `aizynth` | conda | `aizynth` |
| `vina` | conda | `vina` |
| `unidock` | conda | `unidock` |
| `luna` | conda | `luna` |
| `oddt` | conda | `oddt` |
| `reduce` | conda | `reduce` |
| `amber` | conda | `amber` |
| `mmgbsa` | conda | `gmxMMPBSA` |
| `openfe` | conda | `openfe` |
| `ties` | conda | `ties` |
| `pocketminer` | conda | `pocketminer` |
| `pocket2mol` | conda | `Pocket2Mol` |
| `difflinker` | conda | `difflinker` |
| `diffsbdd` | conda | `diffsbdd` |
| `reinvent4` | conda | `reinvent4` |
| `resgen` | conda | `resgen` |
| `targetdiff` | conda | `targetdiff-gpu` |
| `deeppocket` | conda | `deeppocket` |
| `diffdock` | conda | `diffdock` |
| `dock6` | conda | `dock6` |
| `autodock4` | conda | `autodock4` |
| `rdock` | conda | `rdock` |
| `adfrsuite` | conda | `adfrsuite` |
| `lingo3dmol` | conda | `lingo3dmol` |
| `rascore` | conda | `rascore` |
| `rtmscore` | conda | `rtmscore_clean` |
| `tamgen` | conda | `tamgen` |
| `molcraft` | conda | `molcraft` |
| `openmm` | conda | `openmm` |
| `proteins_plus` | http | `https://proteins.plus/api` |

Enabled runtimes are checked before any scientific step and recorded in `tool_runtime_report.json`.

## `context` (shared identity / knowledge)

| key | type | default | range/choices |
|---|---|---|---|
| `pdb_id` | str | `None` |  |
| `protein_name` | str | `None` |  |
| `receptor_file` | str | `None` |  |
| `ligand_file` | str | `None` |  |
| `fragment_file` | str | `None` |  |
| `complex_mmcif_file` | str | `None` |  |
| `target_name` | str | `None` |  |
| `chembl_target_id` | str | `None` |  |
| `uniprot_id` | str | `None` |  |
| `chain` | str | `None` |  |
| `seed` | int | `42` |  |
| `device` | int | `None` | [0, 16] |
| `docking_memory_mode` | str | — | one of off, legacy, ligand_only, ligand_protein |
| `candidate_source` | str | `None` | one of molecule_generation, sampled_library, user_supplied_full |
| `known_site_residues` | list | `None` |  |
| `known_site_provenance` | dict | `None` |  |
| `known_required_interactions` | list | `None` |  |
| `scientific_facts` | dict | `None` |  |

### `context.scientific_facts`

These are declared scientific inputs for tools and capability checks. They never select or enable a tool.

| key | type | default | range/choices |
|---|---|---|---|
| `has_metal` | bool | `None` |  |
| `is_covalent` | bool | `None` |  |
| `protonation_sensitive` | bool | `None` |  |

## `project` (run / directory control)

| key | type | default | range/choices |
|---|---|---|---|
| `project_slug` | str | `None` |  |
| `data_root` | str | `None` |  |
| `project_dir` | str | `None` |  |
| `round_number` | int | `1` | [1, 100] |
| `campaign_state_file` | str | `None` |  |
| `campaign_stage` | str | `None` | one of R1, R2, R3 |
| `campaign_attempt` | int | `None` | [1, ] |
| `round_decision_file` | str | `None` |  |
| `upstream_project_dir` | str | `None` |  |
| `execution_mode` | str | `interrupt` | one of interrupt, auto |
| `resume_existing` | bool | `True` |  |
| `python_executable` | str | `None` |  |

## `runtime_resources` (queue budget, separate from execution)

| key | type | default | range/choices |
|---|---|---|---|
| `wait_timeout_seconds` | int | `0` | [0, 2592000] |
| `slots` | int | `3` | [1, 64] |
| `lock_directory` | str | `None` |  |

## Explicit selection contract

The catalog inventories choices; it does not select them. Agent-authored configs provide every step switch and a complete tool on/off matrix. Registry defaults are recommendation metadata for policy comparison; they do not participate in execution and cannot substitute for that explicit decision.
Ordinarily, every tool under a disabled step is false. The only exception is a producer before `start_step` reused from an audited upstream project: its prior raw step/tool/parameter matrix remains unchanged for lineage, while materialization sets the resolved step to disabled and prevents execution.

Route-changing site/search behavior must be an explicit registered step parameter. It must not be inferred from a holo/apo label, artifact presence, PDB identifier, or another scientific fact inside the runner.

Registered selection strategies:

| step | param | type | default | choices |
|---|---|---|---|---|
| `pocket_analysis` | `site_selection_strategy` | str | **required** | ligand_site, known_site, detector_consensus, provided_box |
| `docking` | `search_space_source` | str | **required** | reference_ligand, docking_box, whole_receptor, precomputed_grid |


## `input_artifacts` (feed a step without running its producer)

For iteration: inject a modified/new file so a downstream step consumes it directly (e.g. a modified library SDF for a re-dock round). Injectable artifacts:

`dock6_grid_bmp`, `dock6_grid_nrg`, `dock6_receptor_site`, `docked_poses`, `docking_box`, `docking_prep_qc`, `docking_species_library`, `library`, `library_manifest`, `microstate_identity_map`, `prefiltered_library`, `prep_manifest`, `prepared_receptor`, `prepared_receptor_noH`, `provisional_receptor`, `receptor_pdbqt`, `receptor_state_manifest`, `reference_ligand`, `selected_poses`, `source_complex_mmcif`, `source_complex_pdb`, `water_inventory`

The producer of an injected artifact must be disabled (often via `start_step`); otherwise validation fails instead of silently overwriting the injected generation.

## Artifacts (inter-step wiring vocabulary)

Inter-step outputs are wired by these canonical names rather than copied into tool params. External paths enter only through registered context or `input_artifacts` fields.

`amber_inpcrd`, `amber_ligand_frcmod`, `amber_ligand_mol2`, `amber_prmtop`, `candidate_topologies`, `coordinator_report`, `deeppocket_predictions`, `dock6_grid_bmp`, `dock6_grid_nrg`, `dock6_receptor_site`, `docked_poses`, `docking_box`, `docking_prep_qc`, `docking_report`, `docking_scores`, `docking_species_library`, `evaluation_report`, `inferred_target_name`, `interaction_blueprint`, `interaction_observations`, `library`, `library_manifest`, `metal_report`, `microstate_identity_map`, `mmgbsa_manifest`, `openfe_rbfe_manifest`, `pocket_analysis_report`, `pocket_center`, `pocket_definition`, `pocket_features`, `prefiltered_library`, `prep_manifest`, `prepared_receptor`, `prepared_receptor_noH`, `property_report`, `protonation_manifest`, `provisional_receptor`, `receptor_pdbqt`, `receptor_state_manifest`, `reference_ligand`, `reference_ligand_microstates`, `report_artifact_manifest`, `report_assets_manifest`, `report_data`, `report_html`, `report_pymol_cover`, `rtmscore_manifest`, `selected_poses`, `selection_manifest`, `short_md_manifest`, `source_complex_mmcif`, `source_complex_pdb`, `structure_qc_report`, `synthesis_report`, `ties_rbfe_manifest`, `trajectory_stability_report`, `water_inventory`, `water_report`

## Two-level metric gate and weighted scoring

The gate has two levels. Level 1 marks one metric on one molecule: a value beyond a hard bound is `fail`, a value outside the soft interval but inside the hard interval is `flag`, and a value inside the inclusive soft interval is `pass`. Missing sides are unbounded and a missing value is `not_applicable`. A mark is never a molecule decision.

Level 2 folds the marks of each parent metric category and is the only place a molecule is decided: `V1` (marking) never rejects, `V2` (standard, the default) rejects on `count(fail) >= 3` or `ratio(fail) >= 50%` inside one category, and `V3` (strict) rejects on `count(fail) >= 2` or `ratio(fail) >= 30%`. The count and ratio branches are independent alternatives and neither carries a minimum fail count, so a lone `fail` among few labeled metrics is ratio 100% and rejects. Missing values leave the ratio; no observed-coverage guard modifies that formula. Any tripped category gives `reject`; otherwise every evaluated category needs strictly more than half of its labeled metrics to be `pass` for `keep`. Other labeled cases give `modify`, and nothing marked gives `not_evaluated`. Rules whose `effect` is `protocol` or `pose_export` are deliberate single-level exceptions: they invalidate a run or a pose rather than folding into a molecule verdict.

Records the fold rejected leave every percentile pool. Each applicable metric among the survivors forms its own pool inside one scientific cohort and uses `rank_avg_percentile_v1`; category scores are percentile-ranked again before StepScore/OverallScore. Bounded zero-weight metrics still mark but contribute exactly zero; finite metrics with no declared bound are score-only and add no pass to the gate denominator.

Metric rules version: `12`.

| metric_id | tool | direction | hard lower | soft lower | soft upper | hard upper | weight magnitude | signed weight | gate? | effect |
|---|---|---|---:|---:|---:|---:|---:|---:|:---:|---|
| `autodock_vina.vina_affinity` | AutoDock Vina | lower | — | — | -7.0 | -1.0 | 5.0 | -5.0 | yes | gate |
| `unidock.docking_score` | Uni-Dock | lower | — | — | -7.0 | -1.0 | 3.0 | -3.0 | yes | gate |
| `gnina.minimizedAffinity` | GNINA | lower | — | — | -7.0 | -1.0 | 5.0 | -5.0 | yes | gate |
| `diffdock.confidence` | DiffDock | higher | — | — | — | — | 0.0 | 0.0 | no | gate |
| `dock6.grid_score` | DOCK 6 | lower | — | — | — | — | 0.0 | 0.0 | no | gate |
| `autodock4.binding_energy` | AutoDock 4 | lower | — | — | -7.0 | -1.0 | 4.0 | -4.0 | yes | gate |
| `rdock.score` | rDock | lower | — | — | — | — | 5.0 | -5.0 | no | gate |
| `adfr.affinity` | ADFR | lower | — | — | -6.0 | 0.0 | 4.0 | -4.0 | yes | gate |
| `autodock_vina.ligand_efficiency` | AutoDock Vina | lower | — | — | -0.25 | -0.05 | 2.0 | -2.0 | yes | gate |
| `autodock_vina.vina_inter` | AutoDock Vina | lower | — | — | — | 0 | 0.0 | 0.0 | yes | gate |
| `autodock_vina.vina_intra_best_pose` | AutoDock Vina | lower | — | — | — | 0 | 0.0 | 0.0 | yes | gate |
| `autodock_vina.engine_pose_rank` | AutoDock Vina | lower | — | — | 20 | 20 | 0.0 | 0.0 | yes | pose_export |
| `autodock_vina.score_gap_1_2` | AutoDock Vina | higher | — | 0.5 | — | — | 0.0 | 0.0 | yes | gate |
| `gnina.intramol` | GNINA | lower | — | — | — | — | 0.0 | 0.0 | no | gate |
| `gnina.CNNscore` | GNINA | higher | 0.2 | 0.5 | — | — | 4.0 | 4.0 | yes | gate |
| `gnina.CNNaffinity` | GNINA | higher | 2.0 | 6.0 | — | — | 3.0 | 3.0 | yes | gate |
| `gnina.CNNvariance` | GNINA | lower | — | — | — | — | 0.0 | 0.0 | no | gate |
| `gnina.engine_pose_rank` | GNINA | lower | — | — | 20 | 20 | 0.0 | 0.0 | yes | pose_export |
| `posebusters.sanitization` | PoseBusters | higher | 1 | 1 | — | — | 0.0 | 0.0 | yes | gate |
| `posebusters.inchi_convertible` | PoseBusters | higher | 1 | 1 | — | — | 0.0 | 0.0 | yes | gate |
| `posebusters.all_atoms_connected` | PoseBusters | higher | 1 | 1 | — | — | 0.0 | 0.0 | yes | gate |
| `posebusters.overall_success_ratio` | PoseBusters | higher | 0.75 | 0.9 | — | — | 3.0 | 3.0 | yes | gate |
| `posebusters.bond_lengths` | PoseBusters | higher | 1 | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.bond_angles` | PoseBusters | higher | 1 | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.internal_steric_clash` | PoseBusters | higher | 1 | 1 | — | — | 2.0 | 2.0 | yes | gate |
| `posebusters.aromatic_ring_flatness` | PoseBusters | higher | 1 | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.non_aromatic_ring_non_flatness` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.double_bond_flatness` | PoseBusters | higher | 1 | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.internal_energy` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.receptor_ligand_maximum_distance` | PoseBusters | higher | 1 | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.minimum_distance_to_receptor` | PoseBusters | higher | 1 | 1 | — | — | 2.0 | 2.0 | yes | gate |
| `posebusters.minimum_distance_to_organic_cofactors` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.minimum_distance_to_inorganic_cofactors` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.minimum_distance_to_waters` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.volume_overlap_with_receptor` | PoseBusters | higher | 1 | 1 | — | — | 2.0 | 2.0 | yes | gate |
| `posebusters.volume_overlap_with_organic_cofactors` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.volume_overlap_with_inorganic_cofactors` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posebusters.volume_overlap_with_waters` | PoseBusters | higher | — | 1 | — | — | 1.0 | 1.0 | yes | gate |
| `posecheck.strain_energy_kcal` | PoseCheck | lower | — | — | 8 | 60 | 4.0 | -4.0 | yes | gate |
| `posecheck.steric_clash_count` | PoseCheck | lower | — | — | 0 | 2 | 5.0 | -5.0 | yes | gate |
| `plip.num_interactions` | PLIP | higher | 2 | — | — | — | 0.0 | 0.0 | yes | gate |
| `prolif.num_interactions` | ProLIF | higher | — | — | — | — | 0.0 | 0.0 | no | record |
| `blueprint_match.verdict` | blueprint_match | higher | 1 | 2 | — | — | 5.0 | 5.0 | yes | gate |
| `blueprint_match.required_hit_rate` | blueprint_match | higher | — | 0.8 | — | — | 0.0 | 0.0 | yes | gate |
| `blueprint_match.desired_hit_rate` | blueprint_match | higher | — | 0.5 | — | — | 0.0 | 0.0 | yes | gate |
| `blueprint_match.desired_ranking_score` | blueprint_match | higher | 0.1 | 0.5 | — | — | 4.0 | 4.0 | yes | gate |
| `blueprint_match.desired_evaluable_coverage` | blueprint_match | higher | — | 0.8 | — | — | 5.0 | 5.0 | yes | protocol |
| `blueprint_match.optional_hits` | blueprint_match | higher | — | — | — | — | 1.0 | 1.0 | no | gate |
| `pocket_occupancy.hydrophobic_region` | pocket_occupancy | higher | 0.01 | 0.2 | — | — | 0.0 | 0.0 | yes | gate |
| `pocket_occupancy.polar_region` | pocket_occupancy | higher | 0.01 | 0.2 | — | — | 0.0 | 0.0 | yes | gate |
| `pocket_occupancy.charged_region` | pocket_occupancy | higher | 0.01 | 0.2 | — | — | 0.0 | 0.0 | yes | gate |
| `pocket_occupancy.aromatic_region` | pocket_occupancy | higher | 0.01 | 0.2 | — | — | 0.0 | 0.0 | yes | gate |
| `pocket_occupancy.overall` | pocket_occupancy | higher | 0.05 | 0.2 | — | — | 3.0 | 3.0 | yes | gate |
| `mmgbsa.delta_g_bind` | MM-GBSA | lower | — | — | 0 | 10 | 4.0 | -4.0 | yes | gate |
| `arpeggio.clash_contacts` | Arpeggio | lower | — | — | 0 | 2 | 3.0 | -3.0 | yes | gate |
| `arpeggio.vdw_clash_contacts` | Arpeggio | lower | — | — | 8 | 15 | 1.0 | -1.0 | yes | gate |
| `arpeggio.proximal_contact_fraction` | Arpeggio | higher | 0.2 | 0.8 | — | — | 0.0 | 0.0 | yes | gate |
| `rdkit_pose.strain_energy` | RDKit UFF/MMFF | lower | — | — | 15.0 | 80.0 | 3.0 | -3.0 | yes | gate |
| `rdkit_pose.tfd` | RDKit TFD | lower | — | — | 0.2 | 0.7 | 2.0 | -2.0 | yes | gate |
| `rdkit_pose.delta_rmsd` | RDKit conformer RMSD | lower | — | — | 2.0 | — | 1.0 | -1.0 | yes | gate |
| `freesasa.buriedness_fraction` | RDKit FreeSASA | higher | 0.1 | 0.3 | — | — | 3.0 | 3.0 | yes | gate |
| `freesasa.contact_surface_fraction` | RDKit FreeSASA | higher | 0.1 | 0.3 | — | — | 1.0 | 1.0 | yes | gate |
| `prolif.tanimoto_to_reference` | ProLIF | higher | — | — | — | — | 0.0 | 0.0 | no | record |
| `mdanalysis.ligand_rmsd` | MDAnalysis | lower | — | — | 2.0 | 5.0 | 4.0 | -4.0 | yes | gate |
| `mdanalysis.pocket_com_drift` | MDAnalysis (supplemental) | lower | — | — | 2.0 | 5.0 | 3.0 | -3.0 | yes | gate |
| `mdanalysis.contact_retention` | MDAnalysis | higher | 0.2 | 0.5 | — | — | 4.0 | 4.0 | yes | gate |
| `rtmscore.ensemble_score` | RTMScore | higher | — | — | — | — | 3.0 | 3.0 | no | gate |
| `rdkit.exact_mol_wt` | RDKit descriptors | target_range | 50 | 150 | 500 | 1050 | 0.0 | 0.0 | yes | gate |
| `rdkit.clogp` | RDKit descriptors | target_range | -3 | 0 | 5 | 10 | 3.0 | 3.0 | yes | gate |
| `rdkit.tpsa` | RDKit descriptors | target_range | 0 | 20 | 140 | 180 | 3.0 | 3.0 | yes | gate |
| `rdkit.hbd` | RDKit descriptors | lower | — | — | 5 | 10 | 0.0 | 0.0 | yes | gate |
| `rdkit.hba` | RDKit descriptors | lower | — | — | 10 | 20 | 0.0 | 0.0 | yes | gate |
| `rdkit.rotatable_bonds` | RDKit descriptors | lower | — | — | 10 | 15 | 0.0 | 0.0 | yes | gate |
| `rdkit.ring_count` | RDKit descriptors | target_range | — | 2 | 4 | 20 | 0.0 | 0.0 | yes | gate |
| `rdkit.aromatic_ring_count` | RDKit descriptors | lower | — | — | 4 | 6 | 0.0 | 0.0 | yes | gate |
| `rdkit.fraction_csp3` | RDKit descriptors | higher | — | 0.2 | — | — | 0.0 | 0.0 | yes | gate |
| `rdkit.formal_charge` | RDKit descriptors | target_range | -3 | -1 | 1 | 3 | 0.0 | 0.0 | yes | gate |
| `rdkit.heavy_atom_count` | RDKit descriptors | target_range | 5 | 12 | 50 | 80 | 0.0 | 0.0 | yes | gate |
| `rdkit.chiral_centers` | RDKit descriptors | lower | — | — | 2 | — | 0.0 | 0.0 | yes | gate |
| `lipinski.ro5_violations` | Lipinski | lower | — | — | 1 | 3 | 0.0 | 0.0 | yes | gate |
| `rdkit.qed` | RDKit.Chem.QED | higher | 0.1 | 0.45 | — | — | 3.0 | 3.0 | yes | gate |
| `rdkit.pains_alert_count` | RDKit PAINS | lower | — | — | 0 | 4 | 1.0 | -1.0 | yes | gate |
| `rdkit.brenk_alert_count` | RDKit FilterCatalog | lower | — | — | 0 | 4 | 1.0 | -1.0 | yes | gate |
| `rdkit.nih_alert_count` | RDKit FilterCatalog | lower | — | — | 0 | 4 | 1.0 | -1.0 | yes | gate |
| `rdkit.sa_score` | RDKit SA_score | lower | — | — | 6 | 9 | 2.0 | -2.0 | yes | gate |
| `aizynthfinder.route_found` | AiZynthFinder | higher | — | 1 | — | — | 0.0 | 0.0 | yes | gate |
| `aizynthfinder.state_score` | AiZynthFinder | higher | — | — | — | — | 2.0 | 2.0 | no | gate |
| `aizynthfinder.reaction_count` | AiZynthFinder | lower | — | — | 5 | — | 0.0 | 0.0 | yes | gate |
| `aizynthfinder.precursor_count` | AiZynthFinder | lower | — | — | — | — | 0.0 | 0.0 | no | gate |
| `aizynthfinder.precursors_in_stock_complete` | AiZynthFinder | higher | — | 1 | — | — | 0.0 | 0.0 | yes | gate |
| `chemfh.other_assay_interference` | ChemFH | lower | — | — | 0.3 | 0.9 | 0.0 | 0.0 | yes | gate |
| `chemfh.blue_fluorescence` | ChemFH | lower | — | — | 0.3 | 0.9 | 0.0 | 0.0 | yes | gate |
| `chemfh.fluc_inhibitors` | ChemFH | lower | — | — | 0.3 | 0.9 | 0.0 | 0.0 | yes | gate |
| `chemfh.promiscuous_compounds` | ChemFH | lower | — | — | 0.3 | 0.9 | 0.5 | -0.5 | yes | gate |
| `chemfh.colloidal_aggregators` | ChemFH | lower | — | — | 0.3 | 0.9 | 0.5 | -0.5 | yes | gate |
| `chemfh.reactive_compounds` | ChemFH | lower | — | — | 0.3 | 0.9 | 1.0 | -1.0 | yes | gate |
| `chemfh.green_fluorescence` | ChemFH | lower | — | — | 0.3 | 0.9 | 0.0 | 0.0 | yes | gate |
| `chemfh.substructure_alert_count` | ChemFH | lower | — | — | 2 | 10 | 2.0 | -2.0 | yes | gate |
| `admet_ai.HIA_Hou` | ADMET-AI | higher | 0.1 | 0.7 | — | — | 4.0 | 4.0 | yes | gate |
| `admet_ai.Caco2_Wang` | ADMET-AI | higher | -7 | -5.15 | — | — | 3.5 | 3.5 | yes | gate |
| `admet_ai.PAMPA_NCATS` | ADMET-AI | higher | 0.1 | 0.7 | — | — | 2.5 | 2.5 | yes | gate |
| `admet_ai.Pgp_Broccatelli` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 3.5 | -3.5 | yes | gate |
| `admet_ai.Bioavailability_Ma` | ADMET-AI | higher | 0.1 | 0.7 | — | — | 4.0 | 4.0 | yes | gate |
| `admet_ai.Solubility_AqSolDB` | ADMET-AI | higher | -9 | -4 | — | — | 4.0 | 4.0 | yes | gate |
| `admet_ai.Lipophilicity_AstraZeneca` | ADMET-AI | target_range | -2 | -1 | 3.5 | 5 | 2.5 | 2.5 | yes | gate |
| `admet_ai.HydrationFreeEnergy_FreeSolv` | ADMET-AI | higher | — | — | — | — | 0.0 | 0.0 | no | gate |
| `admet_ai.BBB_Martins` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 2.0 | -2.0 | yes | gate |
| `admet_ai.PPBR_AZ` | ADMET-AI | lower | — | — | 90 | 98 | 2.0 | -2.0 | yes | gate |
| `admet_ai.VDss_Lombardo` | ADMET-AI | target_range | 0 | 0.04 | 20 | — | 1.0 | 1.0 | yes | gate |
| `admet_ai.cyp_substrate_count` | ADMET-AI | higher | — | 1 | — | — | 0.0 | 0.0 | yes | gate |
| `admet_ai.CYP2C9_Substrate_CarbonMangels` | ADMET-AI | lower | — | — | — | — | 0.1 | -0.1 | no | gate |
| `admet_ai.CYP2D6_Substrate_CarbonMangels` | ADMET-AI | lower | — | — | — | — | 0.1 | -0.1 | no | gate |
| `admet_ai.CYP3A4_Substrate_CarbonMangels` | ADMET-AI | lower | — | — | — | — | 0.5 | -0.5 | no | gate |
| `admet_ai.cyp_inhibition_count` | ADMET-AI | lower | — | — | 1 | 3 | 4.0 | -4.0 | yes | gate |
| `admet_ai.CYP1A2_Veith` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 2.5 | -2.5 | yes | gate |
| `admet_ai.CYP2C19_Veith` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 2.5 | -2.5 | yes | gate |
| `admet_ai.CYP2C9_Veith` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 3.0 | -3.0 | yes | gate |
| `admet_ai.CYP2D6_Veith` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 3.5 | -3.5 | yes | gate |
| `admet_ai.CYP3A4_Veith` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 4.0 | -4.0 | yes | gate |
| `admet_ai.Half_Life_Obach` | ADMET-AI | higher | — | — | — | — | 0.0 | 0.0 | no | gate |
| `admet_ai.Clearance_Hepatocyte_AZ` | ADMET-AI | higher | — | — | — | — | 0.0 | 0.0 | no | gate |
| `admet_ai.Clearance_Microsome_AZ` | ADMET-AI | higher | — | — | — | — | 0.0 | 0.0 | no | gate |
| `admet_ai.hERG` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 5.0 | -5.0 | yes | gate |
| `admet_ai.ClinTox` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 4.0 | -4.0 | yes | gate |
| `admet_ai.AMES` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 5.0 | -5.0 | yes | gate |
| `admet_ai.DILI` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 5.0 | -5.0 | yes | gate |
| `admet_ai.Carcinogens_Lagunin` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 3.5 | -3.5 | yes | gate |
| `admet_ai.LD50_Zhu` | ADMET-AI | lower | — | — | 2.5 | 6.5 | 3.0 | -3.0 | yes | gate |
| `admet_ai.Skin_Reaction` | ADMET-AI | lower | — | — | 0.3 | 0.9 | 3.0 | -3.0 | yes | gate |
| `admet_ai.NR-AR` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 1.5 | -1.5 | yes | gate |
| `admet_ai.NR-AR-LBD` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 1.5 | -1.5 | yes | gate |
| `admet_ai.NR-AhR` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.0 | -2.0 | yes | gate |
| `admet_ai.NR-Aromatase` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.0 | -2.0 | yes | gate |
| `admet_ai.NR-ER` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.0 | -2.0 | yes | gate |
| `admet_ai.NR-ER-LBD` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.0 | -2.0 | yes | gate |
| `admet_ai.NR-PPAR-gamma` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 1.5 | -1.5 | yes | gate |
| `admet_ai.SR-ARE` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.0 | -2.0 | yes | gate |
| `admet_ai.SR-ATAD5` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.5 | -2.5 | yes | gate |
| `admet_ai.SR-HSE` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 1.5 | -1.5 | yes | gate |
| `admet_ai.SR-MMP` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 3.0 | -3.0 | yes | gate |
| `admet_ai.SR-p53` | ADMET-AI | lower | — | — | 0.2 | 0.9 | 2.5 | -2.5 | yes | gate |
| `rascore.score` | RAScore | higher | 0.2 | 0.5 | — | — | 3.0 | 3.0 | yes | gate |
| `redocking.native_ligand_rmsd` | GNINA | lower | — | — | 2.0 | 2.5 | 0.0 | 0.0 | yes | protocol |
| `screening_benchmark.EF1pct` | screening_benchmark | higher | 1.0 | 3.0 | — | — | 0.0 | 0.0 | yes | protocol |
| `screening_benchmark.AUROC` | screening_benchmark | higher | 0.5 | 0.7 | — | — | 0.0 | 0.0 | yes | protocol |

## Steps

| order | step_id | skill | default | consumes | optional consumes | produces |
|---:|---|---|:---:|---|---|---|
| 0 | `structure_prep` | structure-prep | on | — | — | `inferred_target_name`, `prepared_receptor_noH`, `provisional_receptor`, `source_complex_mmcif`, `source_complex_pdb`, `water_inventory`, `reference_ligand`, `protonation_manifest`, `metal_report`, `structure_qc_report`, `prep_manifest` |
| 1 | `pocket_analysis` | pocket-analysis | on | `provisional_receptor`, `prepared_receptor_noH`, `reference_ligand`, `source_complex_mmcif`, `protonation_manifest`, `water_inventory`, `metal_report`, `structure_qc_report` | `docking_box` | `pocket_center`, `pocket_definition`, `docking_box`, `pocket_features`, `interaction_observations`, `interaction_blueprint`, `pocket_analysis_report`, `deeppocket_predictions` |
| 2 | `molecule_generation` | molecule-generator | on | `prepared_receptor_noH`, `reference_ligand`, `pocket_center`, `pocket_definition`, `interaction_blueprint` | `inferred_target_name` | `library`, `library_manifest` |
| 3 | `property_prefilter` | drug-property-filter | on | `library` | — | `prefiltered_library` |
| 4 | `docking` | docking | on | `prepared_receptor_noH`, `prep_manifest`, `source_complex_pdb`, `water_inventory`, `interaction_blueprint`, `prefiltered_library`, `library`, `reference_ligand`, `pocket_definition`, `docking_box`, `library_manifest`, `dock6_receptor_site`, `dock6_grid_nrg`, `dock6_grid_bmp` | — | `docked_poses`, `docking_scores`, `docking_report`, `prepared_receptor`, `receptor_pdbqt`, `receptor_state_manifest`, `water_report`, `docking_species_library`, `microstate_identity_map`, `docking_prep_qc`, `reference_ligand_microstates` |
| 5 | `pose_evaluation` | structure-evaluator | on | `prepared_receptor`, `pocket_definition`, `docked_poses`, `docking_scores`, `docking_report`, `reference_ligand`, `interaction_blueprint`, `pocket_features`, `library_manifest` | — | `selected_poses`, `selection_manifest`, `evaluation_report`, `candidate_topologies`, `rtmscore_manifest`, `short_md_manifest`, `trajectory_stability_report`, `mmgbsa_manifest`, `openfe_rbfe_manifest`, `ties_rbfe_manifest` |
| 6 | `property_filter` | drug-property-filter | on | `selected_poses`, `docked_poses` | — | `property_report`, `synthesis_report` |
| 7 | `coordination` | reflective-coordinator | on | `evaluation_report`, `selection_manifest`, `property_report`, `docking_scores`, `docking_report`, `pocket_analysis_report` | — | `coordinator_report` |
| 8 | `reporting` | report-generator | on (best-effort) | `coordinator_report` | — | `report_html`, `report_data`, `report_assets_manifest`, `report_artifact_manifest`, `report_pymol_cover` |

### `structure_prep` — structure-prep

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `intake` | on | `current:internal` | Download PDB from RCSB and split receptor/ligand. | Starting from a pdb_id or protein_name. | — | `inferred_target_name`, `source_complex_mmcif` |
| `pdbfixer` | on | `drugmate:pdbfixer` | Repair supported standard heavy atoms; canonical receptor remains hydrogen-free. | Use for defensible heavy-atom repairs; short internal gaps require explicit justification. | — | — |
| `propka` | off | `drugmate:propka` | Report pKa / protonation states (diagnostic). | When protonation review matters; not a state rewrite. | — | — |
| `reduce` | off | `reduce:mmtbx.reduce2` | Provisional hydrogens/flips for an explicitly hydrogen-dependent pocket analysis. | Local flip backend when the Reduce2 monomer library is configured. | — | `protonation_manifest`, `provisional_receptor` |
| `ligand_sdf` | on | `drugmate:rdkit` | Reference ligand bound-pose SDF with CCD template bond orders. | Needed for autobox and PoseCheck calibration downstream. | — | `reference_ligand` |
| `metal_validation` | on | `current:internal` | Report retained metals, nearby donors, distances, and coordination counts. | Always safe; mandatory for declared metal sites. | — | `metal_report` |
| `structure_qc` | on | `current:internal` | Final audit of atoms, hydrogens, alternates, waters, metals, and generated files. | Always; provides the Step 0 acceptance report. | — | `structure_qc_report` |
| `residue_validation` | on | `current:internal` | Verify active-site residues are present and unmutated. | When known_site residues must be confirmed. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `ph` | step | float | `None` | [0.0, 14.0] |
| `water_distance_threshold` | intake | float | `3.5` | [0.0, 10.0] |
| `keep_auxiliary_hetero` | intake | bool | `False` |  |
| `preserve_hetero_resnames` | intake | list | `None` |  |
| `skip_prep` | intake | bool | `False` |  |
| `ph` | pdbfixer | float | `None` | [0.0, 14.0] |
| `add_missing_residues` | pdbfixer | bool | `False` |  |
| `max_missing_residues` | pdbfixer | int | `3` | [1, 3] |
| `repair_reason` | pdbfixer | str | `None` |  |
| `seed` | pdbfixer | int | `0` | [0, ] |
| `add_missing_atoms` | pdbfixer | bool | `True` |  |
| `add_hydrogens` | pdbfixer | bool | `False` |  |
| `ph` | propka | float | `None` | [0.0, 14.0] |
| `build` | reduce | bool | `True` |  |
| `flip` | reduce | bool | `True` |  |
| `reason` | reduce | str | **required** |  |
| `ligand_comp_id` | ligand_sdf | str | `None` |  |
| `coordination_cutoff` | metal_validation | float | `3.0` | [1.5, 4.0] |
| `allowed_metals` | metal_validation | list | `['MG', 'ZN', 'CA', 'FE', 'MN', 'NI', 'CU', 'CO', 'NA', 'K', 'CD']` |  |
| `fail_on_altloc` | structure_qc | bool | `False` |  |
| `require_hydrogens` | structure_qc | bool | `None` |  |
| `required_residues` | residue_validation | list | `None` |  |

### `pocket_analysis` — pocket-analysis

Constraints: at least one of `ligand_site`, `fpocket`, `pykvfinder`, `p2rank`, `deeppocket`, `autosite` must be enabled; requires all artifacts `prepared_receptor_noH`.

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `ligand_site` | on | `current:internal` | Define a pocket from the observed reference-ligand shell. | Primary source whenever a usable reference ligand is present. | `artifact:reference_ligand` | — |
| `fpocket` | on | `drugmate:fpocket` | Geometric (Voronoi) pocket detection. | When recommended by plan-pocket-tools, or in the four-tool no-memory fallback; shares the fpocket family with DeepPocket. | — | — |
| `pykvfinder` | off | `drugmate:pyKVFinder` | Pure-Python cavity detection. | When recommended by plan-pocket-tools, or in the four-tool no-memory fallback; contributes anchor-matched native Top-3 candidates. | — | — |
| `p2rank` | on | `drugmate:prank` | ML ligand-binding-site scoring and ranking. | When recommended by plan-pocket-tools, or in the four-tool no-memory fallback; may replace the DeepPocket anchor only when every memory guard passes. | — | — |
| `deeppocket` | off | `deeppocket:python` | DeepPocket segmentation plus learned pocket ranking, normalized as geometry-bearing Step 1 detector evidence. | Required default native Top-1 anchor for detector_consensus; memory may authorize a P2Rank anchor, while the no-memory fallback fixes DeepPocket#1. Shares one detector family with fpocket. | — | `deeppocket_predictions` |
| `autosite` | off | `adfrsuite:autosite` | AutoSite cluster detection with native ranks and exported point centroids. | Only for detector_consensus when included in plan-pocket-tools.execution_tools. | — | — |
| `pocketminer` | off | `pocketminer:python` | Residue-level cryptic-opening propensity annotation; never reranks geometry-bearing candidates or adds detector support. | Explicit cryptic-pocket analysis on one prepared protein chain. Requires at least one separate cavity/site detector for a dockable definition. | — | — |
| `plip` | on | `drugmate:python` | Primary protein-ligand interaction profiler (isolated process). | Whenever a reference ligand is present. | `artifact:reference_ligand` | — |
| `prolif` | on | `drugmate:python` | Record-only interaction fingerprint audit (secondary, isolated process). | For a record-only fingerprint audit when a reference ligand and explicit-H provisional receptor are available; otherwise record not_evaluable and continue. | — | — |
| `arpeggio` | off | `drugmate:python` | Supplement PLIP/ProLIF with atom-level PDBe Arpeggio contacts from the untouched source complex mmCIF. | Holo structures with an RCSB/source mmCIF; degraded chemistry warnings remain observations and cannot promote hard constraints. | `artifact:reference_ligand`, `artifact:source_complex_mmcif` | — |
| `luna` | off | `luna:python` | Supplement core consensus with atom-aware LUNA interactions and a hashed count fingerprint. | Independent interaction cross-check in the dedicated luna environment; PyMOL session generation is intentionally disabled. | `artifact:reference_ligand` | — |
| `oddt` | off | `oddt:python` | Generate ODDT residue interaction counts plus SPLIF hashes as supporting fingerprint evidence. | Legacy independent cross-check in the dedicated oddt environment; never use residue-only observations or SPLIF as a hard gate. | `artifact:reference_ligand` | — |
| `blueprint` | on | `current:internal` | Build 3-tier interaction constraints from eligible PLIP and curated knowledge; record ProLIF separately. | Needed for downstream blueprint_match hard gating. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `site_selection_strategy` | step | str | **required** | one of ligand_site, known_site, detector_consensus, provided_box |
| `pocket_memory_mode` | step | str | `auto` | one of auto, off, required |
| `pocket_memory_index` | step | str | `${SBDD_DATA_ROOT}/docking_memory/pocket_protein_v1.json` |  |
| `pocket_memory_query_npz` | step | str | `None` |  |
| `pocket_memory_query_metadata` | step | str | `None` |  |
| `pocket_memory_preset` | step | str | `conservative` | one of conservative, exploratory_memory_guard |
| `pocket_memory_query_case` | step | str | `None` |  |
| `pocket_memory_query_group` | step | str | `None` |  |
| `pocket_memory_exclude_cases` | step | list | `[]` |  |
| `pocket_memory_exclude_groups` | step | list | `[]` |  |
| `native_step_timeout_seconds` | step | int | `300` | [30, 7200] |
| `fail_on_profiler_error` | step | bool | `False` |  |
| `fail_on_supplemental_profiler_error` | step | bool | `False` |  |
| `consensus_top_k` | step | int | `3` | [1, 20] |
| `center_match_threshold` | step | float | `5.0` | [0.0, 30.0] |
| `residue_iou_threshold` | step | float | `0.3` | [0.0, 1.0] |
| `detector_box_padding` | step | float | `4.0` | [0.0, 30.0] |
| `residue_cutoff` | ligand_site | float | `5.0` | [0.1, 15.0] |
| `box_padding` | ligand_site | float | `4.0` | [0.0, 30.0] |
| `min_box_size` | ligand_site | float | `12.0` | [1.0, 100.0] |
| `min_radius` | fpocket | float | `3.0` | [0.0, 10.0] |
| `max_radius` | fpocket | float | `6.0` | [0.0, 15.0] |
| `probe_in` | pykvfinder | float | `1.4` | [0.0, 10.0] |
| `probe_out` | pykvfinder | float | `4.0` | [0.0, 15.0] |
| `removal_distance` | pykvfinder | float | `2.4` | [0.0, 15.0] |
| `volume_cutoff` | pykvfinder | float | `5.0` | [0.0, 100.0] |
| `max_candidates` | deeppocket | int | `3` | [1, 10] |
| `segmentation_threshold` | deeppocket | float | `0.5` | [0.001, 0.999] |
| `mask_distance` | deeppocket | float | `3.5` | [0.1, 20.0] |
| `timeout_seconds` | deeppocket | int | `1200` | [60, 86400] |
| `timeout_seconds` | autosite | int | `1200` | [60, 86400] |
| `timeout_seconds` | pocketminer | int | `600` | [30, 7200] |
| `rerank_mode` | pocketminer | str | `annotate` | one of annotate, tie_break |
| `min_score_coverage` | pocketminer | float | `0.8` | [0.0, 1.0] |
| `provisional_quantile` | pocketminer | float | `0.9` | [0.5, 0.999] |
| `cluster_ca_cutoff` | pocketminer | float | `8.0` | [1.0, 30.0] |
| `min_cluster_residues` | pocketminer | int | `3` | [1, 100] |
| `max_hypotheses` | pocketminer | int | `5` | [1, 50] |
| `vicinity_cutoff` | prolif | float | `6.5` | [1.0, 15.0] |
| `interaction_types` | prolif | list | `['Hydrophobic', 'HBDonor', 'HBAcceptor', 'PiStacking', 'Anionic', 'Cationic', 'CationPi', 'PiCation', 'VdWContact', 'XBDonor', 'XBAcceptor', 'MetalDonor', 'MetalAcceptor']` |  |
| `ifp_length` | luna | int | `4096` | [128, 65536] |
| `ifp_type` | luna | str | `EIFP` | one of EIFP, HIFP, FIFP |
| `add_dependent_interactions` | luna | bool | `False` |  |
| `strict` | oddt | bool | `True` |  |
| `splif_depth` | oddt | int | `1` | [0, 4] |
| `splif_size` | oddt | int | `4096` | [128, 65536] |
| `splif_distance_cutoff` | oddt | float | `4.5` | [1.0, 10.0] |
| `auto_required_from_tool_consensus` | blueprint | bool | `False` |  |
| `require_curated_provenance` | blueprint | bool | `True` |  |

### `molecule_generation` — molecule-generator

Constraints: at least one of `chembl`, `pocket2mol`, `resgen`, `diffsbdd`, `targetdiff`, `difflinker`, `reinvent4`, `lingo3dmol`, `tamgen`, `molcraft` must be enabled.

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `chembl` | on | — | Retrieve known active compounds as screening candidates from ChEMBL bioactivity data. | Default evidence-backed candidate source; needs a resolvable target id. | `context:chembl_target_id/target_name/uniprot_id` | — |
| `pocket2mol` | on | `pocket2mol:python` | Structure-based de novo generation conditioned on the pocket. | When a prepared pocket PDB + pocket_center are available. | `artifact:prepared_receptor_noH`, `artifact:pocket_center` | — |
| `resgen` | off | `resgen:python` | Fast autoregressive 3D generation around a receptor pocket. | Speed-oriented generation when a prepared receptor and reference ligand are available; complements or replaces Pocket2Mol. | `artifact:prepared_receptor_noH`, `artifact:reference_ligand`, `artifact:pocket_center` | — |
| `diffsbdd` | off | `diffsbdd:python` | Pocket-conditioned full-atom 3D diffusion generation. | Explore a diverse pocket-conditioned 3D candidate set; slower than the speed route and subject to downstream pose QC. | `artifact:prepared_receptor_noH`, `artifact:pocket_center` | — |
| `targetdiff` | off | `targetdiff:python` | Pocket-conditioned 3D diffusion generation with TargetDiff. | Diversity-focused 3D exploration when full 1000-step diffusion cost is acceptable; downstream pose QC is mandatory. | `artifact:prepared_receptor_noH`, `artifact:pocket_center` | — |
| `difflinker` | off | `difflinker:python` | Pocket-conditioned 3D linking of disconnected input fragments. | Only for fragment-linking queries with a trusted 3D fragment file containing at least two disconnected components. | `artifact:prepared_receptor_noH`, `context:fragment_file` | — |
| `reinvent4` | off | `reinvent4:reinvent` | Ligand-only 2D de novo sampling from the reviewed PubChem prior. | Add a fast pocket-agnostic diversity baseline or generate candidates when no structural conditioning is desired. | — | — |
| `lingo3dmol` | off | `lingo3dmol:python` | Generate pocket-conditioned molecules with Lingo3DMol and retain its native 3D coordinates in the canonical library. | Opt-in Step 2 source when an additional pocket-aware 3D generator is wanted; its outputs still pass the common identity, standardization, and downstream pose-QC contracts. | `artifact:prepared_receptor_noH`, `artifact:pocket_center` | — |
| `tamgen` | off | `tamgen:python` | Generate pocket-conditioned SMILES with TamGen from a prepared receptor and the Step 1 pocket center. | Opt-in sequence/coordinate generation route. The default unconditional decoder avoids treating TamGen's placeholder target-side molecule as a real chemical seed. | `artifact:prepared_receptor_noH`, `artifact:pocket_center` | — |
| `molcraft` | off | `molcraft:python` | Sample pocket-conditioned 3D molecules with MolCRAFT and preserve the native coordinates in the canonical library. | Opt-in 3D generation when a prepared receptor and reference ligand define the installed local-sampling route. | `artifact:prepared_receptor_noH`, `artifact:reference_ligand` | — |
| `dude_decoys` | on | — | Download DUD-E property-matched decoys. | Prefer cached/downloaded decoys for retrospective benchmarking. If unavailable or download fails, record the reason and continue without decoys; do not enable local generation automatically. | — | — |
| `rdkit_decoys` | off | — | Optional local BRICS generation of property-matched decoys. | Only when local decoy generation is explicitly requested; unavailable DUD-E data alone does not enable this tool. | — | — |
| `standardization` | on | — | Salt strip, charge handling, canonicalization. | Almost always before conformer generation. | — | — |
| `conformer` | on | — | Preserve generator-native 3D poses; otherwise generate ETKDGv3 + MMFF/UFF conformers. | Always; docking needs 3D coordinates and generator poses need provenance. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `exclude_scaffolds` | step | list | `[]` |  |
| `max_heavy_atoms` | step | int | `None` | [1, 200] |
| `max_mw` | step | float | `None` | [0.0, 2000.0] |
| `activity_type` | chembl | str | `IC50` |  |
| `activity_threshold_nm` | chembl | float | `1000.0` | [0.0, 1000000.0] |
| `max_actives` | chembl | int | `50` | [1, 5000] |
| `num_samples` | pocket2mol | int | `100` | [1, 10000] |
| `timeout_seconds` | pocket2mol | int | `86400` | [60, 86400] |
| `gpu_device` | pocket2mol | int | `0` | [0, 16] |
| `taskset_cpu` | pocket2mol | int | `0` | [-1, 256] |
| `num_samples` | resgen | int | `100` | [1, 10000] |
| `beam_size` | resgen | int | `500` | [1, 5000] |
| `max_steps` | resgen | int | `50` | [1, 200] |
| `timeout_seconds` | resgen | int | `86400` | [60, 86400] |
| `gpu_device` | resgen | int | `0` | [0, 16] |
| `num_samples` | diffsbdd | int | `100` | [1, 10000] |
| `batch_size` | diffsbdd | int | `20` | [1, 1000] |
| `timesteps` | diffsbdd | int | `500` | [1, 500] |
| `pocket_radius` | diffsbdd | float | `10.0` | [3.0, 20.0] |
| `timeout_seconds` | diffsbdd | int | `86400` | [60, 86400] |
| `gpu_device` | diffsbdd | int | `0` | [0, 16] |
| `num_samples` | targetdiff | int | `100` | [1, 10000] |
| `batch_size` | targetdiff | int | `20` | [1, 1000] |
| `num_steps` | targetdiff | int | `1000` | [1, 1000] |
| `pocket_radius` | targetdiff | float | `10.0` | [3.0, 20.0] |
| `timeout_seconds` | targetdiff | int | `86400` | [60, 86400] |
| `gpu_device` | targetdiff | int | `0` | [0, 16] |
| `num_samples` | difflinker | int | `100` | [1, 10000] |
| `linker_size` | difflinker | int | `5` | [1, 50] |
| `n_steps` | difflinker | int | `None` | [1, 1000] |
| `max_batch_size` | difflinker | int | `16` | [1, 1000] |
| `anchor_indices` | difflinker | str | `None` |  |
| `timeout_seconds` | difflinker | int | `86400` | [60, 86400] |
| `gpu_device` | difflinker | int | `0` | [0, 16] |
| `num_samples` | reinvent4 | int | `200` | [1, 100000] |
| `unique_molecules` | reinvent4 | bool | `True` |  |
| `timeout_seconds` | reinvent4 | int | `86400` | [60, 86400] |
| `gpu_device` | reinvent4 | int | `0` | [0, 16] |
| `num_samples` | lingo3dmol | int | `100` | [1, 1000] |
| `gpu_device` | lingo3dmol | int | `0` | [0, 16] |
| `timeout_seconds` | lingo3dmol | int | `86400` | [60, 86400] |
| `max_run_hours` | lingo3dmol | int | `1` | [1, 24] |
| `pocket_radius` | lingo3dmol | float | `10.0` | [3.0, 12.0] |
| `gen_frag_set` | lingo3dmol | int | `40` | [5, 500] |
| `prod_time` | lingo3dmol | int | `1` | [1, 100] |
| `topk` | lingo3dmol | int | `5` | [1, 50] |
| `nci_threshold` | lingo3dmol | float | `0.7` | [0.0, 1.0] |
| `collision_distance` | lingo3dmol | float | `2.5` | [0.5, 5.0] |
| `frag_len_add` | lingo3dmol | int | `15` | [0, 50] |
| `num_samples` | tamgen | int | `100` | [1, 1000] |
| `beam` | tamgen | int | `20` | [1, 100] |
| `max_seeds` | tamgen | int | `5` | [1, 100] |
| `pocket_radius` | tamgen | float | `10.0` | [3.0, 20.0] |
| `sample_beta` | tamgen | float | `1.0` | [0.01, 10.0] |
| `use_conditional` | tamgen | bool | `False` |  |
| `gpu_device` | tamgen | int | `0` | [0, 16] |
| `timeout_seconds` | tamgen | int | `86400` | [60, 86400] |
| `num_samples` | molcraft | int | `10` | [1, 1000] |
| `sample_steps` | molcraft | int | `100` | [1, 1000] |
| `sample_num_atoms` | molcraft | str | `prior` | one of prior, ref |
| `beta1` | molcraft | float | `1.5` | [0.1, 10.0] |
| `sigma1_coord` | molcraft | float | `0.03` | [0.001, 1.0] |
| `sampling_strategy` | molcraft | str | `end_back` | one of vanilla, end_back, end_back_pmf, end_back_mode |
| `gpu_device` | molcraft | int | `0` | [0, 16] |
| `timeout_seconds` | molcraft | int | `86400` | [60, 86400] |
| `num_decoys` | dude_decoys | int | `200` | [0, 10000] |
| `num_decoys` | rdkit_decoys | int | `200` | [0, 10000] |
| `neutralize` | standardization | bool | `False` |  |

### `property_prefilter` — drug-property-filter

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `lipinski` | on | — | Lipinski Rule of Five violation check. | Default drug-likeness gate. | — | — |
| `qed` | on | — | Quantitative estimate of drug-likeness. | Default soft-flag on low QED. | — | — |
| `rascore` | off | `rascore:python` | CPU RAscore XGB accessibility probability in the SYN category. | Low-cost learned accessibility check before retrosynthesis. | — | — |
| `sa_score` | on | — | Synthetic accessibility score. | Default soft-flag on hard-to-make molecules. | — | — |
| `pains` | on | — | All-match PAINS A/B/C, Brenk, and NIH structural-alert screen. | Default first pass. Each PAINS-family, Brenk, and NIH count uses the sbdd numeric band: 0 pass, 1-4 flag, >4 fail; the resulting marks are folded in ALERT before any molecule is rejected. Covalent status is not an automatic rescue; an explicit audited catalog policy is required. | — | — |
| `admet` | off | `admet_ai:python` | Strict batch ADMET-AI predictions with endpoint provenance. | On for final filtering and off for pre-docking prefilter. Missing/incomplete predictions cap a molecule at review_only. | — | — |
| `chemfh` | off | `chemfh:python` | Strict batch ChemFH frequent-hitter and assay-interference predictions with all 17 native rule results. | On for final filtering and off for the pre-docking prefilter. Seven table-matched risk probabilities use 0.30/0.90 level-1 ALERT bounds; their signed weights follow the metric catalog. Incomplete evidence is a protocol review condition, not a molecule fail. Native structural-rule hits remain evidence only; ChemFH PAINS does not duplicate the RDKit gate. | — | — |
| `aizynthfinder` | off | `aizynth:python` | Top-K retrosynthesis planning with AiZynthFinder after RDKit and complete ADMET assessment. | Default for final candidates only; never run across the full pre-docking library. | `admet` | — |
| `dedup` | off | — | Collapse repeated docking poses to one record per compound. | On for post-docking filter; off for raw library prefilter. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `filtered_verdicts` | step | list | `['pass', 'flag', 'fail', 'not_evaluated']` |  |
| `fail_if_no_filtered_molecules` | step | bool | `True` |  |
| `screening_mode` | step | str | `None` | one of V1, V2, V3 |
| `lipinski_max_violations` | lipinski | int | `3` | [0, 4] |
| `qed_threshold` | qed | float | `0.45` | [0.0, 1.0] |
| `timeout_seconds` | rascore | int | `600` | [30, 7200] |
| `sa_threshold` | sa_score | float | `6.0` | [1.0, 10.0] |
| `catalogs` | pains | list | `['PAINS_A', 'PAINS_B', 'PAINS_C', 'BRENK', 'NIH']` | items from PAINS_A, PAINS_B, PAINS_C, BRENK, NIH |
| `hard_fail_catalogs` | pains | list | `['PAINS_A', 'PAINS_B', 'PAINS_C', 'BRENK', 'NIH']` | items from PAINS_A, PAINS_B, PAINS_C, BRENK, NIH |
| `timeout_seconds` | admet | int | `7200` | [30, 7200] |
| `bbb_profile` | admet | str | `non_cns` | one of non_cns, cns |
| `timeout_seconds` | chemfh | int | `7200` | [30, 7200] |
| `seed` | chemfh | int | `42` | [0, 2147483647] |
| `top_k` | aizynthfinder | int | `25` | [1, 500] |
| `timeout_seconds` | aizynthfinder | int | `21600` | [30, 21600] |
| `search_time_limit_seconds` | aizynthfinder | int | `60` | [1, 600] |
| `iteration_limit` | aizynthfinder | int | `100` | [1, 1000] |
| `max_routes` | aizynthfinder | int | `5` | [1, 50] |
| `return_first` | aizynthfinder | bool | `True` |  |
| `expansion_policies` | aizynthfinder | list | `['uspto']` | items from uspto, ringbreaker |
| `filter_policies` | aizynthfinder | list | `['uspto']` | items from uspto |
| `stocks` | aizynthfinder | list | `['zinc']` | items from zinc |

### `docking` — docking

Constraints: exactly one of `gnina`, `unidock`, `autodock_vina`, `diffdock`, `dock6`, `autodock4`, `rdock`, `adfr` must be enabled; exactly one of `reduce`, `pdb2pqr`, `protoss` must be enabled; requires all artifacts `prepared_receptor_noH`; requires at least one artifact from `prefiltered_library`, `library`.

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `docking_prep` | on | `current:internal` | Mandatory final receptor/water/ligand-state preparation before the engine. | Every docking run, after property_prefilter and pocket selection. | — | — |
| `reduce` | off | `reduce:mmtbx.reduce2` | Add hydrogens and optimize His/Asn/Gln flips with mmtbx.reduce2. | Preferred final preparation backend for routine non-covalent protein receptors when target-pH protonation assignment is not required, unless the user specifies another backend. Explicitly enable build=true and flip=true with a passing runtime preflight. Does not apply a target pH; never switch backends silently after failure. | `docking_prep` | `receptor_state_manifest` |
| `pdb2pqr` | off | `drugmate:pdb2pqr30` | Apply pH-aware receptor protonation with PDB2PQR+PROPKA and optimize hydrogen bonds and His/Asn/Gln flips. | Choose for target-pH protonation assignment, protonation-sensitive protocols, or an explicit user request for PDB2PQR+PROPKA. Explicitly enable with justified docking_prep.params.ph and optimize_hbond_network=true. Check support for retained metals, cofactors and nonstandard residues before selection; no universal pH default or automatic runtime fallback. | `docking_prep` | `receptor_state_manifest` |
| `protoss` | off | `proteins_plus:http` | Use ProteinsPlus Protoss for complex-aware protonation and flips. | Unavailable until prepared-complex upload preserves input identity; do not use PDB-ID substitution. | `docking_prep` | `receptor_state_manifest` |
| `water_selection` | off | `current:internal` | Select ligand-bridging crystallographic waters and emit an auditable report. | When a retained reference ligand can define bridging waters, or with manual water residue ids. | `docking_prep` | — |
| `meeko` | off | `reduce:mk_prepare_receptor.py` | Generate the receptor PDBQT consumed by the AutoDock Vina adapter. | Enable for the autodock_vina route after final receptor protonation. AutoDock4 and ADFR use adapter-owned native preparation; GNINA PDB+SDF does not need Meeko. | `docking_prep` | `receptor_pdbqt` |
| `ligand_microstates` | off | `drugmate:dimorphite_dl` | Enumerate reference-ligand protonation microstates with Dimorphite-DL. | Ionizable reference ligands or protonation-sensitive docking studies. | `docking_prep` | `reference_ligand_microstates` |
| `openbabel` | off | `drugmate:obabel` | Contextual format conversion fallback through Open Babel. | Only when a downstream format is unsupported by RDKit/Meeko. | `docking_prep` | — |
| `gnina` | on | `drugmate:gnina` | GNINA docking with CNN scoring + Vina energy. | Default engine; use for routine docking and CNN rescoring. | — | — |
| `unidock` | off | `unidock:unidocktools` | GPU-batched Uni-Dock through its topology-aware unidocktools pipeline. | Large non-covalent libraries when throughput is more important than GNINA CNN scoring; current adapter is validated for shard-sized runs, not a monolithic million-molecule job. | — | — |
| `autodock_vina` | off | `engine=vina:python; io=reduce:python` | Classic AutoDock Vina CPU docking through the Python API, with validated Meeko ligand preparation and SDF export. | Small-to-medium non-covalent libraries, a reproducible classic Vina baseline, or an explicit user request. It is not a strongly orthogonal validator for GNINA/Uni-Dock because all are in the Vina method family. | — | — |
| `diffdock` | off | `diffdock:python` | Run generative DiffDock pose prediction and normalize every ranked pose plus confidence into the canonical docking contract. | Select as the one primary Step 4 engine for a generative docking route. Its confidence is method-specific and is not treated as a Vina-family binding energy. | — | — |
| `dock6` | off | `dock6:${SBDD_DOCK6_SOURCE_DIR}/bin/dock6` | Run the classic DOCK 6 anchor-and-grow/grid-scoring route against explicit precomputed receptor-site and grid artifacts. | Select as the one primary Step 4 engine when an independently prepared and audited DOCK 6 sphere/grid bundle is available. | `artifact:dock6_receptor_site`, `artifact:dock6_grid_nrg`, `artifact:dock6_grid_bmp`, `supplied-artifact:dock6_receptor_site`, `supplied-artifact:dock6_grid_nrg`, `supplied-artifact:dock6_grid_bmp` | — |
| `autodock4` | off | `io=autodock4:python; grid=autodock4:autogrid4; engine=autodock4:autodock4` | Run the classic AutoDock4 Lamarckian genetic-algorithm workflow with native AutoDockTools preparation and AutoGrid4 maps. | Select as the one primary Step 4 engine for a reproducible AutoDock4 baseline. AutoGrid4 and AutoDockTools are helper components of this route, not independent scoring engines. | — | — |
| `rdock` | off | `rdock:rbdock` | Run rDock cavity generation and native docking, preserving its total and component scores without relabeling their units. | Select as the one primary Step 4 engine for an independent classic scoring/search route with a reference-ligand cavity or an explicit sphere-mapped search region. | — | — |
| `adfr` | off | `io=adfrsuite:prepare_receptor; ligand_io=adfrsuite:prepare_ligand; grid=adfrsuite:agfr; engine=adfrsuite:adfr` | Run the ADFRsuite non-covalent workflow with native ligand/receptor preparation, AGFR target maps, and ADFR search. | Select as the one primary Step 4 engine when an ADFR/AGFR maps-based route is requested. Flexible-residue and covalent contracts remain unsupported until explicitly implemented. | — | — |
| `redocking_validation` | on | — | Redock native ligand and compute symmetry-aware RMSD gate. | Default GNINA protocol validation; disable for every other primary engine until engine-neutral redocking is implemented. | `gnina`, `artifact:reference_ligand` | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `search_space_source` | step | str | **required** | one of reference_ligand, docking_box, whole_receptor, precomputed_grid |
| `ph` | docking_prep | float | `None` | [0.0, 14.0] |
| `reason` | docking_prep | str | **required** |  |
| `water_policy` | docking_prep | str | **required** | one of dry, selected |
| `water_reason` | docking_prep | str | **required** |  |
| `ligand_state_policy` | docking_prep | str | `preserve_declared` | one of preserve_declared |
| `ligand_state_reason` | docking_prep | str | **required** |  |
| `build` | reduce | bool | `True` |  |
| `flip` | reduce | bool | `True` |  |
| `force_field` | pdb2pqr | str | `AMBER` | one of AMBER, CHARMM, PARSE |
| `optimize_hbond_network` | pdb2pqr | bool | `True` |  |
| `poll_interval_seconds` | protoss | int | `5` | [1, 60] |
| `timeout_seconds` | protoss | int | `600` | [30, 3600] |
| `ligand_distance` | water_selection | float | `3.5` | [0.0, 10.0] |
| `protein_distance` | water_selection | float | `3.5` | [0.0, 10.0] |
| `min_occupancy` | water_selection | float | `0.5` | [0.0, 1.0] |
| `max_b_factor` | water_selection | float | `60.0` | [0.0, 500.0] |
| `include_residues` | water_selection | list | `[]` |  |
| `exclude_residues` | water_selection | list | `[]` |  |
| `allow_bad_residues` | meeko | bool | `False` |  |
| `ph_tolerance` | ligand_microstates | float | `1.0` | [0.0, 4.0] |
| `max_variants` | ligand_microstates | int | `16` | [1, 128] |
| `output_format` | openbabel | str | `mol2` | one of mol2, pdb, sdf |
| `autobox_add` | gnina | float | `4.0` | [0.0, 30.0] |
| `exhaustiveness` | gnina | int | `16` | [1, 64] |
| `num_modes` | gnina | int | `20` | [1, 100] |
| `cnn_scoring` | gnina | str | `rescore` | one of none, rescore, refinement, all |
| `pose_sort_order` | gnina | str | `CNNscore` | one of CNNscore, CNNaffinity, Energy |
| `device` | gnina | int | `0` | [0, 16] |
| `seed` | gnina | int | `42` |  |
| `gpu_monitor_interval_seconds` | gnina | int | `5` | [1, 600] |
| `timeout_seconds` | gnina | int | `86400` | [60, 86400] |
| `autobox_add` | unidock | float | `4.0` | [0.0, 30.0] |
| `batch_size` | unidock | int | `18000` | [1, 1000000] |
| `scoring_function` | unidock | str | `vina` | one of vina, vinardo |
| `search_mode` | unidock | str | `` | one of , fast, balance, detail |
| `exhaustiveness` | unidock | int | `16` | [1, 1024] |
| `max_step` | unidock | int | `20` | [1, 1000] |
| `num_modes` | unidock | int | `20` | [1, 100] |
| `refine_step` | unidock | int | `3` | [1, 100] |
| `energy_range` | unidock | float | `3.0` | [0.0, 100.0] |
| `device` | unidock | int | `0` | [0, 16] |
| `seed` | unidock | int | `42` |  |
| `timeout_seconds` | unidock | int | `86400` | [1, 86400] |
| `autobox_add` | autodock_vina | float | `4.0` | [0.0, 30.0] |
| `scoring_function` | autodock_vina | str | `vina` | one of vina, vinardo |
| `exhaustiveness` | autodock_vina | int | `16` | [1, 1024] |
| `num_modes` | autodock_vina | int | `20` | [1, 100] |
| `energy_range` | autodock_vina | float | `3.0` | [0.01, 100.0] |
| `min_rmsd` | autodock_vina | float | `1.0` | [0.0, 20.0] |
| `max_evals` | autodock_vina | int | `0` | [0, ] |
| `cpu` | autodock_vina | int | `1` | [1, 256] |
| `seed` | autodock_vina | int | `42` |  |
| `spacing` | autodock_vina | float | `0.375` | [0.01, 2.0] |
| `force_even_voxels` | autodock_vina | bool | `False` |  |
| `no_refine` | autodock_vina | bool | `False` |  |
| `timeout_seconds` | autodock_vina | int | `86400` | [1, 86400] |
| `source_root` | diffdock | str | `${SBDD_DIFFDOCK_SOURCE_DIR}` |  |
| `model_dir` | diffdock | str | `${SBDD_DIFFDOCK_SOURCE_DIR}/workdir/v1.1/score_model` |  |
| `confidence_model_dir` | diffdock | str | `${SBDD_DIFFDOCK_SOURCE_DIR}/workdir/v1.1/confidence_model` |  |
| `checkpoint` | diffdock | str | `best_ema_inference_epoch_model.pt` |  |
| `confidence_checkpoint` | diffdock | str | `best_model_epoch75.pt` |  |
| `samples_per_complex` | diffdock | int | `20` | [1, 40] |
| `batch_size` | diffdock | int | `10` | [1, 64] |
| `inference_steps` | diffdock | int | `20` | [1, 100] |
| `actual_steps` | diffdock | int | `None` | [1, 100] |
| `old_confidence_model` | diffdock | bool | `True` |  |
| `device` | diffdock | int | `0` | [0, 16] |
| `seed` | diffdock | int | `42` |  |
| `timeout_seconds` | diffdock | int | `86400` | [60, 86400] |
| `dock_root` | dock6 | str | `${SBDD_DOCK6_SOURCE_DIR}` |  |
| `antechamber_path` | dock6 | str | `${SBDD_DOCK6_SOURCE_DIR}/bin/antechamber` |  |
| `charge_method` | dock6 | str | `gas` | one of gas, bcc |
| `num_modes` | dock6 | int | `20` | [1, 100] |
| `max_orientations` | dock6 | int | `500` | [1, 100000] |
| `min_anchor_size` | dock6 | int | `5` | [1, 100] |
| `pruning_max_orients` | dock6 | int | `100` | [1, 10000] |
| `score_threshold` | dock6 | float | `1000.0` | [-10000.0, 1000000.0] |
| `seed` | dock6 | int | `42` |  |
| `timeout_seconds` | dock6 | int | `86400` | [60, 86400] |
| `ligand_prep_timeout_seconds` | dock6 | int | `600` | [30, 7200] |
| `autobox_add` | autodock4 | float | `4.0` | [0.0, 30.0] |
| `spacing` | autodock4 | float | `0.375` | [0.01, 2.0] |
| `num_runs` | autodock4 | int | `20` | [1, 256] |
| `num_modes` | autodock4 | int | `20` | [1, 100] |
| `max_evals` | autodock4 | int | `2500000` | [1, 100000000] |
| `population_size` | autodock4 | int | `150` | [2, 10000] |
| `cluster_rmsd` | autodock4 | float | `2.0` | [0.0, 20.0] |
| `seed` | autodock4 | int | `42` | [0, ] |
| `timeout_seconds` | autodock4 | int | `86400` | [60, 86400] |
| `cavity_radius` | rdock | float | `6.0` | [0.1, 50.0] |
| `num_runs` | rdock | int | `25` | [1, 1000] |
| `num_modes` | rdock | int | `20` | [1, 100] |
| `receptor_flex` | rdock | float | `3.0` | [0.0, 20.0] |
| `seed` | rdock | int | `42` | [0, ] |
| `timeout_seconds` | rdock | int | `86400` | [60, 86400] |
| `autobox_add` | adfr | float | `4.0` | [0.0, 30.0] |
| `spacing` | adfr | float | `0.375` | [0.01, 2.0] |
| `num_runs` | adfr | int | `20` | [1, 256] |
| `num_modes` | adfr | int | `20` | [1, 100] |
| `max_evals` | adfr | int | `2500000` | [1, 100000000] |
| `population_size` | adfr | int | `150` | [2, 10000] |
| `max_cores` | adfr | int | `1` | [1, 256] |
| `seed` | adfr | int | `42` | [0, ] |
| `timeout_seconds` | adfr | int | `86400` | [60, 86400] |
| `redocking_rmsd_threshold` | redocking_validation | float | `2.0` | [0.0, 10.0] |
| `redocking_rmsd_hard_threshold` | redocking_validation | float | `2.5` | [0.0, 10.0] |
| `timeout_seconds` | redocking_validation | int | `7200` | [60, 86400] |

### `pose_evaluation` — structure-evaluator

Constraints: requires all artifacts `prepared_receptor`, `docked_poses`, `docking_report`; `blueprint_match` needs one of `plip`.

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `identity_qc` | on | `drugmate:rdkit` | Bind chemical identity, pose lineage, coordinate revision, and content hashes before any pose evidence is merged. | Default integrity gate for every evaluated pose and resumable batch. | — | — |
| `posebusters` | on | `drugmate:posebusters` | Geometric validity (bonds, angles, clashes, flatness). | Default geometric sanity gate. | — | — |
| `posecheck` | on | `drugmate:posecheck` | Strain energy + heavy-atom steric clash count. | Default absolute energetic/clash gate. Optional native-ligand delta calibration must be enabled explicitly. | — | — |
| `plip` | on | `drugmate:python` | Per-pose protein-ligand interaction detection. | Needed on its own and as input to blueprint_match. | — | — |
| `prolif` | off | `drugmate:prolif` | Record-only per-pose RDKit-native interaction fingerprints as a cross-check for PLIP evidence. | Enable for secondary review or MD/multi-pose fingerprints; detector disagreement is reported for audit and does not gate. | — | — |
| `blueprint_match` | on | `current:internal` | Compare pose interactions to the required blueprint. | When enforcing required interactions (e.g. hinge H-bond); requires PLIP; ProLIF is an audit-only cross-check. | — | — |
| `pocket_occupancy` | on | `current:internal` | Sub-pocket contact/occupancy fractions. | To diagnose under-filled sub-pockets. | — | — |
| `rdkit_pose` | off | `drugmate:rdkit` | Independent MMFF94s strain, TFD and reference-conformer RMSD in INTRA. | CPU geometry review using a reproducible sampled conformer minimum. | — | — |
| `freesasa` | off | `drugmate:rdkit.Chem.rdFreeSASA` | Ligand-side buried SASA and surface fractions in POCKET. | CPU ligand burial assessment from the actual docked coordinates. | — | — |
| `rtmscore` | off | `rtmscore:python` | Score post-QC protein-ligand poses with the three-model RTMScore graph-transformer ensemble while preserving canonical pose identity. | Optional machine-learning rescoring after geometric and interaction QC. Scores are an independent method and are never arithmetically mixed with docking energies. | `identity_qc`, `posebusters`, `posecheck` | `rtmscore_manifest` |
| `rescore_rtmscore` | off | `current:internal` | Join RTMScore ensemble outputs to the evaluated cohort by exact canonical pose_id and optionally rerank survivors by that one signal. | Enable together with rtmscore when RTMScore should become the independent post-QC primary ranking metric. | `rtmscore` | — |
| `candidate_topology` | off | `amber:antechamber` | Build candidate-specific Amber complex/ligand/receptor topologies with ff14SB/GAFF2 and AM1-BCC, emitting a per-pose manifest. | Explicit top-K escalation before MM-GBSA or short MD; run only on pose-QC survivors. Generic metal and covalent systems are reported as not_evaluable. | `identity_qc`, `posebusters`, `posecheck` | `candidate_topologies` |
| `openmm_md` | off | `openmm:python` | Run short explicit-solvent periodic OpenMM MD for a bounded post-QC candidate subset and retain trajectories/state data. | Optional pose-stability escalation after candidate topology construction. This is a short diagnostic trajectory, not a binding free-energy calculation. | `identity_qc`, `posebusters`, `posecheck`, `candidate_topology` | `short_md_manifest` |
| `mdanalysis` | off | `openmm:python` | Analyze OpenMM trajectories after protein alignment and report ligand RMSD, center displacement, and native-contact retention. | Enable with openmm_md when trajectory-stability screening evidence is required. Shared sbdd numeric bounds are heuristic gates, not universal experimental cutoffs. | `openmm_md` | `trajectory_stability_report` |
| `mmgbsa` | off | `mmgbsa:MMPBSA.py` | Run AmberTools MMPBSA.py on post-QC candidate topology bundles and emit a pose-bound, single-frame endpoint reranking manifest. | Explicit top-K MM-GBSA escalation after candidate topology build. Disabled for routine pose QC; current scores are single-frame, no-entropy estimates rather than absolute binding free energies. | `identity_qc`, `posebusters`, `posecheck`, `candidate_topology` | `mmgbsa_manifest` |
| `openfe_rbfe` | off | `openfe:python` | Validate post-QC ligand identities/chemistry and optionally build an OpenFE ligand-mapping network for later RBFE work. | Contextual RBFE escalation. Generic FEP/RBFE requests select this backend, but routine pose_evaluation leaves it disabled. | `identity_qc`, `posebusters`, `posecheck` | `openfe_rbfe_manifest` |
| `ties_rbfe` | off | `ties:python` | Validate post-QC ligand identities and write an auditable TIES pair/campaign plan without claiming a computed free energy. | Only when the user explicitly requests TIES RBFE planning. The generic RBFE default is OpenFE; TIES is never a silent fallback. | `identity_qc`, `posebusters`, `posecheck` | `ties_rbfe_manifest` |
| `arpeggio` | off | `drugmate:python` | Candidate pose arpeggio interactions/fingerprints. | Optional per-pose assessment; Arpeggio numeric QC, LUNA/ODDT audit fingerprints. | — | — |
| `luna` | off | `luna:python` | Candidate pose luna interactions/fingerprints. | Optional per-pose assessment; Arpeggio numeric QC, LUNA/ODDT audit fingerprints. | — | — |
| `oddt` | off | `oddt:python` | Candidate pose oddt interactions/fingerprints. | Optional per-pose assessment; Arpeggio numeric QC, LUNA/ODDT audit fingerprints. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `screening_mode` | step | str | `None` | one of V1, V2, V3 |
| `selection_strategy` | step | str | `ranked_per_compound` | one of ranked_per_compound, first_passing_per_compound, all, top_by_cnnscore, best_per_compound |
| `poses_per_compound` | step | int | `20` | [1, 100] |
| `batch_size` | step | int | `50` | [1, 5000] |
| `profiler_workers` | step | int | `2` | [1, 4] |
| `checkpoint_interval_seconds` | step | int | `300` | [0, 86400] |
| `max_poses` | step | int | `None` | [1, 100000] |
| `overall_timeout_seconds` | step | int | `86400` | [60, 86400] |
| `num_conformers` | posecheck | int | `10` | [1, 1000] |
| `strain_warn_kcal` | posecheck | float | `8.0` | [0.0, 1000.0] |
| `strain_fail_kcal` | posecheck | float | `60.0` | [0.0, 10000.0] |
| `clash_warn_count` | posecheck | int | `0` | [0, 1000] |
| `clash_fail_count` | posecheck | int | `2` | [0, 1000] |
| `calibrate_against_reference` | posecheck | bool | `False` |  |
| `vicinity_cutoff` | prolif | float | `6.0` | [2.0, 12.0] |
| `interaction_types` | prolif | list | `None` |  |
| `num_conformers` | rdkit_pose | int | `8` | [1, 64] |
| `seed` | rdkit_pose | int | `42` | [0, 2147483647] |
| `max_iterations` | rdkit_pose | int | `500` | [1, 10000] |
| `probe_radius` | freesasa | float | `1.4` | [0.1, 3.0] |
| `model_ensemble` | rtmscore | list | `[1, 2, 3]` | items from 1, 2, 3 |
| `pocket_cutoff_angstrom` | rtmscore | float | `10.0` | [3.0, 20.0] |
| `batch_size` | rtmscore | int | `128` | [1, 4096] |
| `num_workers` | rtmscore | int | `0` | [0, 64] |
| `parallel_graph_build` | rtmscore | bool | `False` |  |
| `device` | rtmscore | str | `auto` | one of auto, cuda, cpu |
| `timeout_seconds` | rtmscore | int | `86400` | [60, 86400] |
| `use_as_primary_metric` | rescore_rtmscore | bool | `True` |  |
| `execution_mode` | candidate_topology | str | `build` | one of preflight, build |
| `protein_force_field` | candidate_topology | str | `ff14SB` | one of ff14SB, ff19SB |
| `ligand_force_field` | candidate_topology | str | `gaff2` | one of gaff2 |
| `ligand_charge_method` | candidate_topology | str | `bcc` | one of bcc |
| `max_candidates` | candidate_topology | int | `20` | [1, 1000] |
| `timeout_per_candidate_seconds` | candidate_topology | int | `2400` | [60, 86400] |
| `fail_fast` | candidate_topology | bool | `True` |  |
| `allow_receptor_hetero` | candidate_topology | bool | `False` |  |
| `top_k` | openmm_md | int | `3` | [1, 100] |
| `platform` | openmm_md | str | `CUDA` | one of CUDA |
| `precision` | openmm_md | str | `mixed` | one of mixed, single, double |
| `device_index` | openmm_md | int | `0` | [0, 16] |
| `temperature_kelvin` | openmm_md | float | `300.0` | [250.0, 400.0] |
| `pressure_bar` | openmm_md | float | `1.0` | [0.1, 100.0] |
| `timestep_femtoseconds` | openmm_md | float | `2.0` | [0.5, 2.0] |
| `friction_per_picosecond` | openmm_md | float | `1.0` | [0.01, 100.0] |
| `minimization_max_iterations` | openmm_md | int | `500` | [0, 100000] |
| `equilibration_picoseconds` | openmm_md | float | `20.0` | [0.0, 10000.0] |
| `production_picoseconds` | openmm_md | float | `100.0` | [1.0, 100000.0] |
| `report_interval_picoseconds` | openmm_md | float | `2.0` | [0.01, 1000.0] |
| `solvent_padding_angstrom` | openmm_md | float | `10.0` | [6.0, 30.0] |
| `water_model` | openmm_md | str | `tip3p` | one of tip3p |
| `seed` | openmm_md | int | `42` |  |
| `timeout_per_candidate_seconds` | openmm_md | int | `86400` | [60, 604800] |
| `fail_fast` | openmm_md | bool | `False` |  |
| `alignment_selection` | mdanalysis | str | `protein and backbone` |  |
| `ligand_selection` | mdanalysis | str | `resname LIG and not name H*` |  |
| `contact_cutoff_angstrom` | mdanalysis | float | `4.5` | [2.0, 10.0] |
| `analysis_stride` | mdanalysis | int | `1` | [1, 10000] |
| `timeout_per_candidate_seconds` | mdanalysis | int | `3600` | [30, 86400] |
| `execution_mode` | mmgbsa | str | `preflight` | one of preflight, score |
| `top_k` | mmgbsa | int | `20` | [1, 1000] |
| `gb_model` | mmgbsa | int | `5` | one of 1, 2, 5, 7, 8 |
| `salt_concentration_molar` | mmgbsa | float | `0.15` | [0.0, 5.0] |
| `timeout_per_candidate_seconds` | mmgbsa | int | `900` | [30, 86400] |
| `fail_fast` | mmgbsa | bool | `False` |  |
| `keep_intermediate_files` | mmgbsa | bool | `False` |  |
| `use_as_primary_metric` | mmgbsa | bool | `True` |  |
| `execution_mode` | openfe_rbfe | str | `preflight` | one of preflight, plan |
| `network_method` | openfe_rbfe | str | `minimal_spanning` | one of minimal_spanning, minimal_redundant, radial |
| `mapper` | openfe_rbfe | str | `KartografAtomMapper` | one of KartografAtomMapper, LomapAtomMapper |
| `mapper_settings` | openfe_rbfe | dict | `None` |  |
| `max_candidates` | openfe_rbfe | int | `8` | [2, 200] |
| `min_ligands` | openfe_rbfe | int | `2` | [2, 200] |
| `allow_charge_changes` | openfe_rbfe | bool | `False` |  |
| `timeout_seconds` | openfe_rbfe | int | `600` | [1, 86400] |
| `n_processes` | openfe_rbfe | int | `1` | [1, 32] |
| `mst_num` | openfe_rbfe | int | `2` | [1, 10] |
| `central_pose_id` | openfe_rbfe | str | `None` |  |
| `min_mapping_score` | openfe_rbfe | float | `0.3` | [0.0, 1.0] |
| `min_mapped_heavy_fraction` | openfe_rbfe | float | `0.5` | [0.0, 1.0] |
| `require_explicit_hydrogens` | openfe_rbfe | bool | `True` |  |
| `allow_undefined_stereochemistry` | openfe_rbfe | bool | `False` |  |
| `execution_mode` | ties_rbfe | str | `preflight` | one of preflight, plan |
| `network_strategy` | ties_rbfe | str | `star` | one of star |
| `reference_pose_id` | ties_rbfe | str | `None` |  |
| `max_candidates` | ties_rbfe | int | `8` | [2, 100] |
| `require_post_qc_lineage` | ties_rbfe | bool | `True` |  |
| `runtime_probe_timeout_seconds` | ties_rbfe | int | `30` | [1, 600] |
| `simulation_platform` | ties_rbfe | str | `auto` | one of auto, CUDA, OpenCL, CPU |
| `timeout_seconds` | arpeggio | int | `300` | [10, 7200] |
| `ifp_length` | luna | int | `4096` | [128, 65536] |
| `ifp_type` | luna | str | `EIFP` | one of EIFP, HIFP, FIFP |
| `add_dependent_interactions` | luna | bool | `False` |  |
| `timeout_seconds` | luna | int | `300` | [10, 7200] |
| `strict` | oddt | bool | `True` |  |
| `splif_depth` | oddt | int | `1` | [0, 4] |
| `splif_size` | oddt | int | `4096` | [128, 65536] |
| `splif_distance_cutoff` | oddt | float | `4.5` | [1.0, 10.0] |
| `timeout_seconds` | oddt | int | `300` | [10, 7200] |

### `property_filter` — drug-property-filter

Constraints: requires at least one artifact from `selected_poses`, `docked_poses`.

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `lipinski` | on | — | Lipinski Rule of Five violation check. | Default drug-likeness gate. | — | — |
| `qed` | on | — | Quantitative estimate of drug-likeness. | Default soft-flag on low QED. | — | — |
| `rascore` | off | `rascore:python` | CPU RAscore XGB accessibility probability in the SYN category. | Low-cost learned accessibility check before retrosynthesis. | — | — |
| `sa_score` | on | — | Synthetic accessibility score. | Default soft-flag on hard-to-make molecules. | — | — |
| `pains` | on | — | All-match PAINS A/B/C, Brenk, and NIH structural-alert screen. | Default first pass. Each PAINS-family, Brenk, and NIH count uses the sbdd numeric band: 0 pass, 1-4 flag, >4 fail; the resulting marks are folded in ALERT before any molecule is rejected. Covalent status is not an automatic rescue; an explicit audited catalog policy is required. | — | — |
| `admet` | on | `admet_ai:python` | Strict batch ADMET-AI predictions with endpoint provenance. | On for final filtering and off for pre-docking prefilter. Missing/incomplete predictions cap a molecule at review_only. | — | — |
| `chemfh` | on | `chemfh:python` | Strict batch ChemFH frequent-hitter and assay-interference predictions with all 17 native rule results. | On for final filtering and off for the pre-docking prefilter. Seven table-matched risk probabilities use 0.30/0.90 level-1 ALERT bounds; their signed weights follow the metric catalog. Incomplete evidence is a protocol review condition, not a molecule fail. Native structural-rule hits remain evidence only; ChemFH PAINS does not duplicate the RDKit gate. | — | — |
| `aizynthfinder` | on | `aizynth:python` | Top-K retrosynthesis planning with AiZynthFinder after RDKit and complete ADMET assessment. | Default for final candidates only; never run across the full pre-docking library. | `admet` | — |
| `dedup` | on | — | Collapse repeated docking poses to one record per compound. | On for post-docking filter; off for raw library prefilter. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `screening_mode` | step | str | `None` | one of V1, V2, V3 |
| `lipinski_max_violations` | lipinski | int | `3` | [0, 4] |
| `qed_threshold` | qed | float | `0.45` | [0.0, 1.0] |
| `timeout_seconds` | rascore | int | `600` | [30, 7200] |
| `sa_threshold` | sa_score | float | `6.0` | [1.0, 10.0] |
| `catalogs` | pains | list | `['PAINS_A', 'PAINS_B', 'PAINS_C', 'BRENK', 'NIH']` | items from PAINS_A, PAINS_B, PAINS_C, BRENK, NIH |
| `hard_fail_catalogs` | pains | list | `['PAINS_A', 'PAINS_B', 'PAINS_C', 'BRENK', 'NIH']` | items from PAINS_A, PAINS_B, PAINS_C, BRENK, NIH |
| `timeout_seconds` | admet | int | `7200` | [30, 7200] |
| `bbb_profile` | admet | str | `non_cns` | one of non_cns, cns |
| `timeout_seconds` | chemfh | int | `7200` | [30, 7200] |
| `seed` | chemfh | int | `42` | [0, 2147483647] |
| `top_k` | aizynthfinder | int | `25` | [1, 500] |
| `timeout_seconds` | aizynthfinder | int | `21600` | [30, 21600] |
| `search_time_limit_seconds` | aizynthfinder | int | `60` | [1, 600] |
| `iteration_limit` | aizynthfinder | int | `100` | [1, 1000] |
| `max_routes` | aizynthfinder | int | `5` | [1, 50] |
| `return_first` | aizynthfinder | bool | `True` |  |
| `expansion_policies` | aizynthfinder | list | `['uspto']` | items from uspto, ringbreaker |
| `filter_policies` | aizynthfinder | list | `['uspto']` | items from uspto |
| `stocks` | aizynthfinder | list | `['zinc']` | items from zinc |

### `coordination` — reflective-coordinator

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `require_quality_gates` | on | — | Block nomination on protocol gates; modify requests a rerun and reject invalidates the protocol. | Default; keep on for trustworthy nominations. | — | — |
| `screening_benchmark` | off | — | Compute AUROC/EF from ACT_/DEC_ labels when present. | Only with verified ACT_ and DEC_ labels in the actual screening cohort. Numeric or unknown IDs alone provide no activity labels; disable for unlabeled screening. | — | — |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `pass_rate_threshold` | step | float | `0.3` | [0.0, 1.0] |
| `min_candidates` | step | int | `3` | [1, 1000] |
| `max_rounds` | step | int | `3` | [1, 20] |
| `min_validation_actives` | step | int | `20` | [1, 100000] |
| `min_validation_decoys` | step | int | `50` | [1, 1000000] |
| `min_validation_auroc` | step | float | `0.7` | [0.0, 1.0] |
| `min_validation_auroc_hard` | step | float | `0.5` | [0.0, 1.0] |
| `min_validation_ef1` | step | float | `3.0` | [0.0, ] |
| `min_validation_ef1_hard` | step | float | `1.0` | [0.0, ] |
| `normalization_method` | step | str | `rank_avg_percentile_v1` | one of rank_avg_percentile_v1 |
| `metric_weight_overrides` | step | dict | `None` |  |
| `category_weight_overrides` | step | dict | `None` |  |
| `screening_mode` | step | str | `None` | one of V1, V2, V3 |

### `reporting` — report-generator

| tool | default | runtime | purpose | when to use | needs | produces |
|---|:---:|---|---|---|---|---|
| `rdkit_2d` | on | — | RDKit 2D structure depictions on candidate cards. | Default. | — | — |
| `mol_3d` | on | — | 3Dmol.js per-pose viewer cards. | Default. | — | — |
| `molstar_main` | on | — | Mol* main viewer (pocket + reference ligand). | Default top viewer. | — | — |
| `plots` | on | — | matplotlib QC / druglikeness / ROC-EF panels. | Default. | — | — |
| `pymol_cover` | off | — | Static PyMOL cover PNG. | Optional; off by default. | — | `report_pymol_cover` |

Parameters:

| param | scope | type | default | range/choices |
|---|---|---|---|---|
| `top_n` | step | int | `10` | [1, 1000] |
| `title` | step | str | `None` |  |
| `assets_source_dir` | step | str | `None` |  |


"""Additional pose scalars with explicit coordinates, references and NA reasons."""
from __future__ import annotations

import hashlib
from pathlib import Path


def rdkit_pose_metrics(mol, *, num_conformers=8, seed=42, max_iterations=500):
    """MMFF94s strain/TFD/RMSD relative to a converged sampled minimum.

    Relax hydrogen positions with heavy atoms fixed before measuring docked
    energy. The reference is an independently embedded ETKDG ensemble and does
    not include the docked pose; it is a sampled minimum, not a claim of the
    global minimum.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, TorsionFingerprints, rdMolAlign

    if mol is None or not mol.GetNumConformers():
        raise ValueError("a parsed 3D pose is required")
    if not 1 <= int(num_conformers) <= 64:
        raise ValueError("num_conformers must be between 1 and 64")
    pose = Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)
    organic_elements = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}
    if any(atom.GetAtomicNum() not in organic_elements for atom in pose.GetAtoms()) or not AllChem.MMFFHasAllMoleculeParams(pose):
        return {"success": True, "status": "not_applicable",
                "reason": "metal/nonstandard chemistry outside this organic MMFF94s protocol", "metrics": {}}
    props = AllChem.MMFFGetMoleculeProperties(pose, mmffVariant="MMFF94s")
    force = AllChem.MMFFGetMoleculeForceField(pose, props)
    for atom in pose.GetAtoms():
        if atom.GetAtomicNum() > 1:
            force.AddFixedPoint(atom.GetIdx())
    if force.Minimize(maxIts=int(max_iterations)) != 0:
        raise RuntimeError("fixed-heavy-atom hydrogen minimization did not converge")
    docked_energy = float(force.CalcEnergy())
    reference = Chem.Mol(pose)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = 1
    params.clearConfs = True
    embedded_ids = list(AllChem.EmbedMultipleConfs(
        reference, numConfs=int(num_conformers), params=params))
    if not embedded_ids or not reference.GetNumConformers():
        return {
            "success": True, "status": "not_evaluable", "metrics": {},
            "reason": "ETKDGv3 generated no reference conformers; strain, TFD and RMSD cannot be evaluated",
            "reference": {"method": "independent ETKDGv3 ensemble; MMFF94s minimum",
                          "seed": int(seed), "conformers_requested": int(num_conformers),
                          "conformers_generated": 0, "conformers_converged": 0},
        }
    energies = AllChem.MMFFOptimizeMoleculeConfs(
        reference, numThreads=1, maxIters=int(max_iterations), mmffVariant="MMFF94s")
    converged = [(energy, conformer.GetId())
                 for conformer, (status, energy) in zip(reference.GetConformers(), energies)
                 if status == 0]
    if not converged:
        raise RuntimeError("no sampled MMFF94s reference conformer converged")
    minimum, best = min(converged)
    light_pose, light_reference = Chem.RemoveHs(pose), Chem.RemoveHs(reference)
    rmsd = float(rdMolAlign.GetBestRMS(
        Chem.Mol(light_pose), light_reference, prbId=0, refId=best, maxMatches=10000))
    non_ring, ring = TorsionFingerprints.CalculateTorsionLists(light_pose)
    tfd = None
    if non_ring or ring:
        tfd = float(TorsionFingerprints.GetTFDBetweenMolecules(
            light_pose, light_reference, confId1=0, confId2=best))
    return {
        "success": True, "status": "available",
        "metrics": {"strain_energy_kcal_mol": max(0.0, docked_energy - minimum),
                    "tfd": tfd, "delta_rmsd_angstrom": rmsd},
        "protocol": {"force_field": "MMFF94s", "conformer_method": "ETKDGv3",
                     "num_conformers": int(num_conformers), "seed": int(seed),
                     "max_iterations": int(max_iterations)},
        "reference": {"method": "independent ETKDGv3 ensemble; MMFF94s minimum",
                      "seed": int(seed), "conformers_requested": int(num_conformers),
                      "conformers_converged": len(converged), "conformer_id": best,
                      "minimum_energy_kcal_mol": float(minimum),
                      "hydrogen_relaxed_docked_energy_kcal_mol": docked_energy},
        "not_applicable_metrics": ["tfd"] if tfd is None else [],
    }


class LigandSASA:
    """Heavy-atom SASA with one consistent RDKit van der Waals radius scheme."""

    def __init__(self, receptor_file, probe_radius=1.4):
        from rdkit import Chem
        from rdkit.Chem import rdFreeSASA

        self.receptor = Chem.MolFromPDBFile(str(receptor_file), sanitize=False, removeHs=False)
        if self.receptor is None or not self.receptor.GetNumConformers():
            raise ValueError("could not read receptor coordinates for FreeSASA")
        self.receptor = Chem.RemoveHs(self.receptor, sanitize=False)
        self.probe_radius = float(probe_radius)
        if not 0 < self.probe_radius <= 3:
            raise ValueError("SASA probe radius must be in (0, 3] angstrom")
        self.options = rdFreeSASA.SASAOpts()
        self.options.probeRadius = self.probe_radius
        table = Chem.GetPeriodicTable()
        self.receptor_radii = [
            table.GetRvdw(atom.GetAtomicNum()) for atom in self.receptor.GetAtoms()
        ]
        if not all(value > 0 for value in self.receptor_radii):
            raise ValueError("SASA radius unavailable for a receptor atom")
        self.free_receptor_sasa = float(
            rdFreeSASA.CalcSASA(
                self.receptor, self.receptor_radii, opts=self.options,
            )
        )
        if self.free_receptor_sasa <= 0:
            raise ValueError("free receptor SASA is zero")
        self.receptor_hash = hashlib.sha256(Path(receptor_file).read_bytes()).hexdigest()

    def calculate(self, mol):
        from rdkit import Chem
        from rdkit.Chem import rdFreeSASA

        ligand = Chem.RemoveHs(mol)
        table = Chem.GetPeriodicTable()

        def radii(molecule):
            values = [table.GetRvdw(atom.GetAtomicNum()) for atom in molecule.GetAtoms()]
            if not all(value > 0 for value in values):
                raise ValueError("SASA radius unavailable for an atom")
            return values

        ligand_radii = radii(ligand)
        free = float(
            rdFreeSASA.CalcSASA(ligand, ligand_radii, opts=self.options)
        )
        complex_mol = Chem.CombineMols(self.receptor, ligand)
        complex_sasa = float(
            rdFreeSASA.CalcSASA(
                complex_mol,
                [*self.receptor_radii, *ligand_radii],
                opts=self.options,
            )
        )
        receptor_atoms = range(self.receptor.GetNumAtoms())
        ligand_atoms = range(
            self.receptor.GetNumAtoms(), complex_mol.GetNumAtoms()
        )
        bound_receptor = sum(
            complex_mol.GetAtomWithIdx(i).GetDoubleProp("SASA")
            for i in receptor_atoms
        )
        bound = sum(
            complex_mol.GetAtomWithIdx(i).GetDoubleProp("SASA")
            for i in ligand_atoms
        )
        if free <= 0:
            raise ValueError("free ligand SASA is zero")
        # Numerical quadrature can introduce tiny negative differences.
        buried = min(free, max(0.0, free - bound))
        # The executable policy defines the interface contact area as one half of the SASA
        # lost by the two separated partners.  This is not generally identical
        # to the ligand-side burial used by ``buriedness_fraction``.
        contact_surface = 0.5 * (
            self.free_receptor_sasa + free - complex_sasa
        )
        contact_surface = min(free, max(0.0, contact_surface))
        return {
            "success": True, "status": "available",
            "metrics": {"free_receptor_sasa_A2": self.free_receptor_sasa,
                        "bound_receptor_sasa_A2": bound_receptor,
                        "free_ligand_sasa_A2": free, "bound_ligand_sasa_A2": bound,
                        "complex_sasa_A2": complex_sasa,
                        "buried_ligand_sasa_A2": buried, "buriedness_fraction": buried / free,
                        "contact_surface_area_A2": contact_surface,
                        "contact_surface_area_fraction": contact_surface / free},
            "protocol": {"radii": "RDKit periodic-table van der Waals; heavy atoms",
                         "probe_radius_angstrom": self.probe_radius,
                         "contact_surface_definition": "0.5 * (free receptor SASA + free ligand SASA - complex SASA)",
                         "receptor_sha256": self.receptor_hash},
        }


def attach_pose_geometry(evaluation, receptor_file, docked_sdf, tool_configs):
    """Attach before the final pose fold; failures retain missing-value status."""
    from rdkit import Chem

    configs = {key: tool_configs.get(key) or {} for key in ("rdkit_pose", "freesasa")}
    if not any(config.get("enabled", False) for config in configs.values()):
        return []
    sasa, sasa_error = None, None
    if configs["freesasa"].get("enabled", False):
        try:
            sasa = LigandSASA(receptor_file, **(configs["freesasa"].get("params") or {}))
        except Exception as exc:
            sasa_error = str(exc)
    supplier = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)
    molecules = list(supplier)
    errors = []
    for entry in evaluation:
        index = int(entry["pose_index"])
        for tool, config in configs.items():
            if not config.get("enabled", False):
                continue
            try:
                mol = molecules[index]
                if mol is None or not mol.GetNumConformers():
                    raise ValueError("unreadable ligand coordinates")
                if tool == "rdkit_pose":
                    result = rdkit_pose_metrics(mol, **(config.get("params") or {}))
                else:
                    if sasa is None:
                        raise RuntimeError(sasa_error)
                    result = sasa.calculate(mol)
                entry[tool] = result
                metric_ids = {
                    "rdkit_pose": ["rdkit_pose.strain_energy", "rdkit_pose.tfd", "rdkit_pose.delta_rmsd"],
                    "freesasa": ["freesasa.buriedness_fraction", "freesasa.contact_surface_fraction"],
                }[tool]
                missing = metric_ids if result.get("status") in {"not_applicable", "not_evaluable"} else (
                    ["rdkit_pose.tfd"] if "tfd" in result.get("not_applicable_metrics", []) else [])
                if missing:
                    entry.setdefault("metric_context", {}).setdefault("not_applicable_metric_ids", []).extend(missing)
            except Exception as exc:
                entry[tool] = {"success": False, "status": "failed", "metrics": {}, "error": str(exc)}
                errors.append({"tool": tool, "pose_index": index, "reason": str(exc)})
            entry.setdefault("evidence_status", {})[tool] = (
                "tool_failed" if entry[tool].get("success") is False else entry[tool].get("status", "available"))
    return errors


def merge_md_metrics(evaluation, result):
    """Join analyzed trajectories to their exact canonical input pose IDs."""
    indexed = {entry["pose_id"]: entry for entry in evaluation if entry.get("pose_id")}
    if len(indexed) != sum(bool(entry.get("pose_id")) for entry in evaluation):
        raise ValueError("duplicate evaluation pose_id")
    seen, staged = set(), []
    for candidate in result.get("candidates") or []:
        pose_id = candidate.get("pose_id")
        if pose_id not in indexed or pose_id in seen:
            raise ValueError("unknown or duplicate MD result pose_id")
        seen.add(pose_id)
        analysis = candidate.get("analysis") or {}
        if candidate.get("analysis_status") == "passed":
            if analysis.get("success") is not True or int(analysis.get("frame_count", 0)) < 2:
                raise ValueError("MD result claims success without a valid analyzed trajectory")
            analysis = {**analysis, "protocol": result.get("protocol") or {}}
            staged.append((indexed[pose_id], analysis))
    for entry, analysis in staged:
        entry["mdanalysis"] = analysis
        entry.setdefault("evidence_status", {})["mdanalysis"] = "available"
        for field, metric_id in (("ligand_pocket_com_distance_drift_angstrom", "mdanalysis.pocket_com_drift"),
                                 ("initial_contact_retention_fraction", "mdanalysis.contact_retention")):
            block = (analysis.get("metrics") or {}).get(field) or {}
            if block.get("applicability") in {"not_applicable", "not_evaluable"}:
                entry.setdefault("metric_context", {}).setdefault("not_applicable_metric_ids", []).append(metric_id)
    # Selection budgets are intentional exclusions, not failed calculations.
    for entry in evaluation:
        if entry.get("pose_id") not in seen:
            entry.setdefault("metric_context", {}).setdefault("not_applicable_metric_ids", []).extend(
                ["mdanalysis.ligand_rmsd", "mdanalysis.pocket_com_drift", "mdanalysis.contact_retention"])
    return len(staged)

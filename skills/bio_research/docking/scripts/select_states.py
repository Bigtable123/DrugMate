#!/usr/bin/env python3
"""Build an audited docking subset while preserving source and graph identity."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def select_states(input_file: Path, state_ids: list[str], output_dir: Path) -> dict:
    """Select stable compound IDs from an existing docking input.

    This does not infer completed work from a partially written pose file. The
    caller chooses IDs from the docking ledger and retains previous attempts.
    """
    from rdkit import Chem

    input_file = Path(input_file).resolve()
    config = json.loads(input_file.read_text())
    if (
        not isinstance(state_ids, list)
        or not state_ids
        or any(not isinstance(item, str) or not item.strip() for item in state_ids)
        or len(set(state_ids)) != len(state_ids)
    ):
        raise ValueError("state IDs must be a nonempty list of unique strings")
    source = Path(config["ligand_file"]).resolve()
    keep = set(state_ids)
    selected_molecules = []
    source_ids = []
    seen_ids = set()
    for index, molecule in enumerate(Chem.SDMolSupplier(str(source), removeHs=False)):
        if molecule is None:
            raise ValueError(f"unparseable docking input record: {index}")
        compound_id = next((
            molecule.GetProp(prop).strip()
            for prop in ("compound_id", "CompoundID", "compound_name", "_Name")
            if molecule.HasProp(prop) and molecule.GetProp(prop).strip()
        ), None)
        if compound_id is None:
            raise ValueError(f"docking input record has no stable compound ID: {index}")
        if compound_id in seen_ids:
            raise ValueError(f"ambiguous duplicate compound ID in docking input: {compound_id}")
        if compound_id in keep:
            selected_molecules.append(molecule)
        source_ids.append(compound_id)
        seen_ids.add(compound_id)
    unknown = set(state_ids) - set(source_ids)
    if unknown:
        raise ValueError(f"state IDs absent from docking input: {sorted(unknown)}")
    output_dir = Path(output_dir).resolve()
    # A new attempt directory prevents a retry from replacing any original
    # receptor, library, configuration or docking result.
    output_dir.mkdir(parents=True, exist_ok=False)
    selected = output_dir / "selected_ligands.sdf"
    writer = Chem.SDWriter(str(selected))
    try:
        for molecule in selected_molecules:
            writer.write(molecule)
    finally:
        writer.close()
    config["ligand_file"] = str(selected)
    config_path = output_dir / "docking_input.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    report = {
        "source_input": str(input_file),
        "source_library": str(source),
        "source_library_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "selected_state_ids": [compound_id for compound_id in source_ids if compound_id in keep],
        "omitted_from_input_state_ids": [
            compound_id for compound_id in source_ids if compound_id not in keep
        ],
        "selected_state_records": len(keep),
        "docking_input_file": str(config_path),
    }
    (output_dir / "selection_manifest.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-file", type=Path, required=True)
    parser.add_argument(
        "--state-ids-file",
        type=Path,
        required=True,
        help="JSON array of stable compound IDs to execute (including explicit state IDs)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="New directory for this attempt's inputs",
    )
    args = parser.parse_args()
    report = select_states(
        args.input_file, json.loads(args.state_ids_file.read_text()), args.output_dir
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

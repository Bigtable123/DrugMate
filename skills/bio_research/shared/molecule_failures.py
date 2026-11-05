"""Auditable per-molecule failures without batch-wide rejection thresholds."""
from __future__ import annotations

from collections import Counter


class MoleculeFailureTracker:
    """Retain each failed record and persist it through the optional callback."""

    def __init__(self, on_failure=None):
        self.failures: list[dict] = []
        self.on_failure = on_failure

    def record(self, compound: dict, operation: str, code: str, reason: str,
               *, exception: Exception | None = None, details: dict | None = None):
        entry = {
            "status": "error", "tool": operation.split(".")[0],
            "operation": operation, "code": code, "reason": reason,
            "kind": "unexpected_exception" if exception else "molecule_error",
            **{key: compound[key] for key in (
                "compound_id", "input_compound_id", "parent_entity_id",
                "source_record_id", "candidate_role", "source", "name",
                "chembl_id", "smiles", "smiles_standardized",
                "source_file", "source_mol_index",
            ) if compound.get(key) is not None},
        }
        if details:
            entry["details"] = details
        if exception is not None:
            root = exception
            while root.__cause__ is not None:
                root = root.__cause__
            entry.update(exception_type=type(root).__name__, exception_message=str(root))
        self.failures.append(entry)
        if self.on_failure:
            self.on_failure(self)

    def report(self):
        return {
            "schema_version": 1,
            "failure_events": len(self.failures),
            "counts_by_operation": dict(Counter(x["operation"] for x in self.failures)),
            "failures": self.failures,
        }


def is_rdkit_parse_error(exc: Exception) -> bool:
    """Recognized invalid molecular chemistry, not arbitrary runtime errors."""
    from rdkit import Chem
    return isinstance(exc, Chem.MolSanitizeException)

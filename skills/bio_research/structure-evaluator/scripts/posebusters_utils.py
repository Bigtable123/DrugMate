"""PoseBusters utilities: geometric validation of docked poses.

PoseBusters checks bond lengths, bond angles, internal clashes,
protein-ligand clashes, flatness, and double bond geometry.

Install: pip install posebusters
"""
from __future__ import annotations

import logging
import math
from numbers import Real
from pathlib import Path

logger = logging.getLogger(__name__)


def check_posebusters_available() -> bool:
    try:
        from posebusters import PoseBusters
        return True
    except ImportError:
        return False


def run_posebusters(docked_sdf: Path, receptor_file: Path) -> list[dict]:
    """Run PoseBusters on each pose in the docked SDF.

    Returns per-pose checks: True/False for measured results, None plus a
    check_errors reason when the check could not produce a boolean result.
    """
    from posebusters import PoseBusters

    pose_names = _read_sdf_names(docked_sdf)
    pb = PoseBusters(config="dock")
    df = pb.bust(str(docked_sdf), None, str(receptor_file))
    rows = list(df.iterrows())
    upcast_boolean_columns = {
        col
        for col in df.columns
        if _is_upcast_boolean_column(row[col] for _idx, row in rows)
    }

    results = []
    for position, (_idx, row) in enumerate(rows):
        # PoseBusters returns one row per input pose in SDF order. Recover the
        # true SDF record index from the DataFrame's own row label when it
        # encodes it (PoseBusters uses a (file, molecule_index) MultiIndex), so
        # a pose that PoseBusters could not load does not silently shift every
        # subsequent result onto the wrong molecule.
        pose_index = _row_pose_index(_idx, position)
        mol_name = pose_names[pose_index] if pose_index < len(pose_names) else f"pose_{pose_index}"
        checks = {}
        check_errors = {}
        violations = []
        for col in df.columns:
            val = row[col]
            # PoseBusters check columns are booleans, but pandas yields them as
            # numpy.bool_, for which isinstance(val, bool) is False. Coerce via
            # a real bool test so failing checks are actually counted.
            bool_val = _as_bool(
                val,
                allow_binary_numeric=col in upcast_boolean_columns,
            )
            if bool_val is not None:
                checks[col] = bool_val
                if not bool_val:
                    violations.append(col)
            else:
                # Energy/conformer calculation failures can return NaN even
                # when other checks succeeded. Keep unavailable evidence out
                # of both pass/fail counts and strict JSON batch fingerprints.
                checks[col] = None
                if isinstance(val, Real) and not math.isfinite(val):
                    check_errors[col] = "non_finite_result"
                elif val is None or type(val).__name__ in {"NAType", "NaTType"}:
                    check_errors[col] = "missing_result"
                else:
                    check_errors[col] = "non_boolean_result"

        # A pose PoseBusters could not load yields an all-NaN row, so no check
        # resolves to a bool. "Zero violations" is then absence of evidence, not
        # a pass: report it as a tool failure for this pose instead of a silent
        # all_passed=True.
        num_checks = sum(1 for value in checks.values() if isinstance(value, bool))
        all_passed = None
        if violations:
            all_passed = False
        elif num_checks and not check_errors:
            all_passed = True
        results.append({
            "pose_index": pose_index,
            "mol_name": mol_name,
            "success": num_checks > 0,
            "num_checks": num_checks,
            "all_passed": all_passed,
            "violations": violations,
            "num_violations": len(violations),
            "checks": checks,
            "check_errors": check_errors,
        })

    return results


def _as_bool(value, *, allow_binary_numeric: bool = False) -> bool | None:
    """Return a Python bool for genuine boolean values, else None.

    Handles numpy.bool_ (whose isinstance(bool) is False) and pandas NA. Numeric
    or string cells are treated as non-boolean by default.  ``allow_binary_numeric``
    is reserved for a check column proven to be pandas-upcast booleans because
    one or more rows are missing.
    """
    if isinstance(value, bool):
        return value
    # numpy.bool_ and pandas nullable booleans expose a real bool via bool(),
    # but numeric types (0/1) would also coerce, so restrict to boolean dtypes.
    type_name = type(value).__name__
    if type_name in ("bool_", "BooleanDtype") or type_name.startswith("bool"):
        try:
            return bool(value)
        except (TypeError, ValueError):
            return None
    if allow_binary_numeric and isinstance(value, Real):
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError):
            return None
        if math.isfinite(numeric) and numeric in (0.0, 1.0):
            return bool(int(numeric))
    return None


def _is_upcast_boolean_column(values) -> bool:
    """Detect a boolean check column promoted to floats by missing values.

    Pandas represents ``[True, NaN, False]`` as ``[1.0, NaN, 0.0]`` in some
    PoseBusters result frames.  Requiring both a missing value and exclusively
    binary finite numeric observations avoids accepting arbitrary numeric tool
    output as boolean evidence.
    """
    saw_missing = False
    saw_binary_numeric = False
    for value in values:
        if value is None or type(value).__name__ in {"NAType", "NaTType"}:
            saw_missing = True
            continue
        if isinstance(value, Real):
            try:
                numeric = float(value)
            except (TypeError, ValueError, OverflowError):
                return False
            if not math.isfinite(numeric):
                saw_missing = True
                continue
            if _as_bool(value) is not None:
                continue
            if numeric in (0.0, 1.0):
                saw_binary_numeric = True
                continue
        elif _as_bool(value) is not None:
            continue
        return False
    return saw_missing and saw_binary_numeric


def _row_pose_index(row_label, fallback: int) -> int:
    """Extract the SDF record index from a PoseBusters DataFrame row label."""
    candidate = row_label
    if isinstance(row_label, tuple):
        candidate = row_label[-1]
    try:
        return int(candidate)
    except (TypeError, ValueError):
        return fallback


def _read_sdf_names(sdf_path: Path) -> list[str]:
    from rdkit import Chem

    names = []
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    for i, mol in enumerate(supplier):
        if mol is not None and mol.HasProp("_Name") and mol.GetProp("_Name").strip():
            names.append(mol.GetProp("_Name").strip())
        else:
            names.append(f"pose_{i}")
    return names

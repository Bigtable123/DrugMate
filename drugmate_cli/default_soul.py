"""Locate the fixed project persona, shared by every Drugmate profile."""

from pathlib import Path
import sysconfig


def get_fixed_soul_path() -> Path:
    """Find this checkout's SOUL.md or the copy shipped with an installation.

    Do not search the working directory or DRUGMATE_HOME: neither can override
    the identity that belongs to the running Drugmate installation.
    """
    source_root = Path(__file__).resolve().parent.parent
    if (source_root / "drugmate").is_file() and (source_root / "pyproject.toml").is_file():
        # If a checkout has lost its SOUL.md, report that error instead of
        # silently selecting the persona from a separate installed version.
        return source_root / "SOUL.md"

    relative_path = Path("share") / "drugmate" / "SOUL.md"
    # pip --target places data files alongside the installed Python packages.
    candidates = [source_root / relative_path]
    user_scheme = sysconfig.get_preferred_scheme("user")
    user_library = sysconfig.get_path("purelib", scheme=user_scheme)
    user_data = sysconfig.get_path("data", scheme=user_scheme)
    if user_library and user_data and source_root == Path(user_library).resolve():
        candidates.append(Path(user_data) / relative_path)
    data_root = sysconfig.get_path("data")
    if data_root:
        candidates.append(Path(data_root) / relative_path)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Drugmate's fixed SOUL.md is missing. Restore the checkout's SOUL.md "
        "or reinstall Drugmate with its share/drugmate/SOUL.md data file."
    )


def read_fixed_soul_md() -> str:
    """Read the canonical persona; a broken installation must not change identity."""
    soul_path = get_fixed_soul_path()
    try:
        content = soul_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise RuntimeError(f"Cannot read Drugmate's fixed persona at {soul_path}: {exc}") from exc
    if not content.strip():
        raise RuntimeError(f"Drugmate's fixed persona is empty: {soul_path}")
    return content

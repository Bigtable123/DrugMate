"""Read legacy local pin metadata to preserve explicit deletion protection."""

import json
from typing import Any, Dict

from drugmate_runtime.drugmate_constants import get_drugmate_home


def get_record(skill_name: str) -> Dict[str, Any]:
    """Return a skill's existing metadata without recording usage activity."""
    path = get_drugmate_home() / "skills" / ".usage.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(data, dict):
        return {}
    record = data.get(skill_name)
    return record if isinstance(record, dict) else {}

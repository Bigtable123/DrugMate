"""Local skill listing, reading and bundled-source synchronization."""
from __future__ import annotations

import json
import shlex


def _result(action: str, name: str | None = None) -> dict:
    if action == "sync":
        from tools.skills_sync import sync_skills
        return sync_skills(quiet=True)
    from tools.skills_tool import skills_list, skill_view
    if action in {"view", "inspect"} and name:
        return json.loads(skill_view(name, preprocess=False))
    if action in {"list", ""}:
        return json.loads(skills_list())
    return {"error": "Use skills list, skills view <name>, or skills sync."}


def skills_command(args) -> None:
    action = getattr(args, "skills_action", None) or "list"
    result = _result(action, getattr(args, "name", None))
    if action in {"view", "inspect"} and "content" in result:
        print(result["content"])
    else:
        print(json.dumps(result, ensure_ascii=False, indent=2))


def handle_skills_slash(command: str, console) -> None:
    try:
        parts = shlex.split(command)
    except ValueError as exc:
        console.print(str(exc), markup=False)
        return
    args = parts[1:]
    result = _result(args[0] if args else "list", args[1] if len(args) > 1 else None)
    if "content" in result:
        console.print(result["content"], markup=False)
    else:
        console.print(json.dumps(result, ensure_ascii=False, indent=2), markup=False)

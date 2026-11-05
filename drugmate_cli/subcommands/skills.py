"""Local-only ``drugmate skills`` commands."""
from __future__ import annotations

from typing import Callable


def build_skills_parser(subparsers, *, cmd_skills: Callable) -> None:
    parser = subparsers.add_parser("skills", help="List, read, sync, or configure local skills")
    actions = parser.add_subparsers(dest="skills_action")
    actions.add_parser("list", help="List local skills")
    view = actions.add_parser("view", aliases=["inspect"], help="Read a local skill")
    view.add_argument("name", help="Skill name or path")
    actions.add_parser("sync", help="Sync bundled skills from this source checkout")
    actions.add_parser("config", help="Enable or disable local skills")
    parser.set_defaults(func=cmd_skills)

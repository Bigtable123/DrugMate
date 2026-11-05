"""``drugmate tools`` subcommand parser.

Extracted from ``drugmate_cli/main.py:main()`` (god-file Phase 2 follow-up).
Handler injected to avoid importing ``main``.
"""

from __future__ import annotations

from typing import Callable


def build_tools_parser(subparsers, *, cmd_tools: Callable) -> None:
    """Attach the ``tools`` subcommand to ``subparsers``."""
    tools_parser = subparsers.add_parser(
        "tools",
        help="Configure which tools are enabled in the CLI",
        description=(
            "Enable, disable, or list tools for the Python CLI.\n\n"
            "Built-in toolsets use plain names (e.g. terminal, file, skills).\n"
            "Run 'python drugmate tools' with no subcommand for interactive configuration."
        ),
    )
    tools_parser.add_argument(
        "--summary",
        action="store_true",
        help="Print a summary of enabled CLI tools and exit",
    )
    tools_sub = tools_parser.add_subparsers(dest="tools_action")

    # drugmate tools list [--platform cli]
    tools_list_p = tools_sub.add_parser(
        "list",
        help="Show all tools and their enabled/disabled status",
    )
    tools_list_p.add_argument(
        "--platform",
        default="cli",
        choices=["cli"],
        help="Platform to show (default: cli)",
    )

    # drugmate tools disable <name...> [--platform cli]
    tools_disable_p = tools_sub.add_parser(
        "disable",
        help="Disable built-in toolsets",
    )
    tools_disable_p.add_argument(
        "names",
        nargs="+",
        metavar="NAME",
        help="Toolset name (e.g. terminal)",
    )
    tools_disable_p.add_argument(
        "--platform",
        default="cli",
        choices=["cli"],
        help="Platform to apply to (default: cli)",
    )

    # drugmate tools enable <name...> [--platform cli]
    tools_enable_p = tools_sub.add_parser(
        "enable",
        help="Enable built-in toolsets",
    )
    tools_enable_p.add_argument(
        "names",
        nargs="+",
        metavar="NAME",
        help="Toolset name",
    )
    tools_enable_p.add_argument(
        "--platform",
        default="cli",
        choices=["cli"],
        help="Platform to apply to (default: cli)",
    )

    tools_parser.set_defaults(func=cmd_tools)

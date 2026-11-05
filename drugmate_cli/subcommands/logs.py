"""``drugmate logs`` subcommand parser.

Extracted verbatim from ``drugmate_cli/main.py:main()`` (god-file Phase 2).
Handler injected to avoid importing ``main``.
"""

from __future__ import annotations

import argparse
from typing import Callable


def build_logs_parser(subparsers, *, cmd_logs: Callable) -> None:
    """Attach the ``logs`` subcommand to ``subparsers``."""
    # =========================================================================
    # logs command
    # =========================================================================
    logs_parser = subparsers.add_parser(
        "logs",
        help="View and filter Drugmate log files",
        description="View, tail, and filter agent.log / errors.log",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
    drugmate logs                    Show last 50 lines of agent.log
    drugmate logs -f                 Follow agent.log in real time
    drugmate logs errors             Show last 50 lines of errors.log
    drugmate logs --level WARNING    Only show WARNING and above
    drugmate logs --session abc123   Filter by session ID
    drugmate logs --component tools  Only show tool-related lines
    drugmate logs --since 1h         Lines from the last hour
    drugmate logs --since 30m -f     Follow, starting from 30 min ago
    drugmate logs list               List available log files with sizes
""",
    )
    logs_parser.add_argument(
        "log_name",
        nargs="?",
        default="agent",
        help="Log to view: agent (default), errors, or 'list' to show available files",
    )
    logs_parser.add_argument(
        "-n",
        "--lines",
        type=int,
        default=50,
        help="Number of lines to show (default: 50)",
    )
    logs_parser.add_argument(
        "-f",
        "--follow",
        action="store_true",
        help="Follow the log in real time (like tail -f)",
    )
    logs_parser.add_argument(
        "--level",
        metavar="LEVEL",
        help="Minimum log level to show (DEBUG, INFO, WARNING, ERROR)",
    )
    logs_parser.add_argument(
        "--session",
        metavar="ID",
        help="Filter lines containing this session ID substring",
    )
    logs_parser.add_argument(
        "--since",
        metavar="TIME",
        help="Show lines since TIME ago (e.g. 1h, 30m, 2d)",
    )
    logs_parser.add_argument(
        "--component",
        metavar="NAME",
        help="Filter by component: agent, tools, cli",
    )
    logs_parser.set_defaults(func=cmd_logs)

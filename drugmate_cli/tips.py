"""Random tips shown at CLI session start to help users discover features."""

import random


# ---------------------------------------------------------------------------
# Tip corpus — one-liners covering slash commands, CLI flags, config,
# keybindings, tools, skills, profiles, and workflow tricks.
# ---------------------------------------------------------------------------

TIPS = [
    '/title "my project" names your session — resume it later with /resume or python drugmate --cli -c.',
    '/resume picks up where you left off in a previously named session.',
    '/branch forks the current session so you can explore a different direction without losing progress.',
    '/compress manually compresses conversation context when things get long.',
    '/queue <prompt> queues a message for the next turn without interrupting the current one.',
    '/usage shows token usage, cost breakdown, and session duration.',
    '/model lets you select a model and provider for the conversation.',
    '/tools lists the tools available to your session.',
    '/reload-skills reloads the local skill catalog after you edit a skill.',
    '/agents shows active delegated tasks.',
    '@file:path/to/file.py includes a local file in your message.',
    'Use a separate DRUGMATE_HOME for development so configuration and sessions stay separate.',
    'The terminal and process tools let the agent monitor and stop long-running commands.',
    'Review the requested inputs and execution plan before confirming a scientific run.',
    'The bio_research/pipeline-entry skill defines the scientific workflow and its confirmation gates.',
    'Keep scientific tool paths and credentials in local configuration, outside the published source.',
    'Start a research conversation with python drugmate --cli.',
]


def get_random_tip(exclude_recent: int = 0) -> str:
    """Return a random tip string.

    Args:
        exclude_recent: not used currently; reserved for future
            deduplication across sessions.
    """
    return random.choice(TIPS)

"""Registry of built-in inference-provider profiles.

Only modules in the fixed source-release list below are loaded. Profiles retain
provider-specific credentials, endpoints, request shaping and model catalogs;
local directories and installed entry points cannot override them.
"""

from __future__ import annotations

import importlib
import threading

from providers.base import OMIT_TEMPERATURE, ProviderProfile  # noqa: F401

_REGISTRY: dict[str, ProviderProfile] = {}
_ALIASES: dict[str, str] = {}
_discovered = False
_lock = threading.RLock()

_BUILTIN_MODULES = (
    "alibaba",
    "alibaba_coding_plan",
    "anthropic",
    "arcee",
    "azure_foundry",
    "bedrock",
    "copilot",
    "copilot_acp",
    "custom",
    "deepseek",
    "gemini",
    "gmi",
    "huggingface",
    "kilocode",
    "kimi_coding",
    "minimax",
    "nous",
    "novita",
    "nvidia",
    "ollama_cloud",
    "openai_codex",
    "opencode_zen",
    "openrouter",
    "qwen_oauth",
    "stepfun",
    "xai",
    "xiaomi",
    "zai",
)


def register_provider(profile: ProviderProfile) -> None:
    """Register a shipped profile by its canonical name and aliases."""
    with _lock:
        _REGISTRY[profile.name] = profile
        for alias in profile.aliases:
            _ALIASES[alias] = profile.name


def _load_builtin_providers() -> None:
    global _discovered
    with _lock:
        if _discovered:
            return
        for module_name in _BUILTIN_MODULES:
            importlib.import_module(f"providers.builtin.{module_name}")
        _discovered = True


def get_provider_profile(name: str) -> ProviderProfile | None:
    """Return a built-in profile by name or alias, or None for generic APIs."""
    _load_builtin_providers()
    return _REGISTRY.get(_ALIASES.get(name, name))


def list_providers() -> list[ProviderProfile]:
    """Return one profile for each built-in canonical provider name."""
    _load_builtin_providers()
    return list(_REGISTRY.values())

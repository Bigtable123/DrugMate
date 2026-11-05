"""Nous model-account helpers and optional managed terminal status."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Set
from drugmate_cli.config import load_config
from drugmate_cli.nous_account import (
    NousPortalAccountInfo, format_nous_portal_entitlement_message,
    get_nous_portal_account_info,
)
from tools.managed_tool_gateway import is_managed_tool_gateway_ready
from tools.tool_backend_helpers import has_direct_modal_credentials, normalize_modal_mode, resolve_modal_backend_state

_DEFAULT_PLATFORM_TOOLSETS = {"cli": "drugmate-cli"}
MANAGED_FEATURE_COVERAGE_CATEGORY = {"modal": "modal"}

@dataclass(frozen=True)
class NousFeatureState:
    key: str
    label: str
    included_by_default: bool
    available: bool
    active: bool
    managed_by_nous: bool
    direct_override: bool
    toolset_enabled: bool
    current_provider: str = ""
    explicit_configured: bool = False

@dataclass(frozen=True)
class NousSubscriptionFeatures:
    subscribed: bool
    nous_auth_present: bool
    provider_is_nous: bool
    features: Dict[str, NousFeatureState]
    account_info: Optional[NousPortalAccountInfo] = None

    @property
    def modal(self) -> NousFeatureState:
        return self.features["modal"]

    def items(self) -> Iterable[NousFeatureState]:
        return self.features.values()


def _model_config_dict(config: Dict[str, object]) -> Dict[str, object]:
    model_cfg = config.get("model")
    if isinstance(model_cfg, dict):
        return dict(model_cfg)
    if isinstance(model_cfg, str) and model_cfg.strip():
        return {"default": model_cfg.strip()}
    return {}



def _toolset_enabled(config: Dict[str, object], toolset_key: str) -> bool:
    from drugmate_runtime.toolsets import resolve_toolset

    platform_toolsets = config.get("platform_toolsets")
    if not isinstance(platform_toolsets, dict) or not platform_toolsets:
        platform_toolsets = {"cli": [_DEFAULT_PLATFORM_TOOLSETS["cli"]]}

    target_tools = set(resolve_toolset(toolset_key))
    if not target_tools:
        return False

    for platform, raw_toolsets in platform_toolsets.items():
        if isinstance(raw_toolsets, list):
            toolset_names = list(raw_toolsets)
        else:
            default_toolset = _DEFAULT_PLATFORM_TOOLSETS.get(platform)
            toolset_names = [default_toolset] if default_toolset else []
        if not toolset_names:
            default_toolset = _DEFAULT_PLATFORM_TOOLSETS.get(platform)
            if default_toolset:
                toolset_names = [default_toolset]

        available_tools: Set[str] = set()
        for toolset_name in toolset_names:
            if not isinstance(toolset_name, str) or not toolset_name:
                continue
            try:
                available_tools.update(resolve_toolset(toolset_name))
            except Exception:
                continue

        if target_tools and target_tools.issubset(available_tools):
            return True

    return False



def ensure_nous_portal_access(
    *,
    capability: str = "the Nous Tool Gateway",
    coverage_category: Optional[str] = None,
) -> bool:
    """Make sure the user is entitled to the Nous Tool Gateway, logging in if
    needed.

    Used by ``drugmate tools`` when a user selects a Nous-managed Tool Gateway
    backend (e.g. "Firecrawl (Nous Portal)").  Unlike ``drugmate model``'s Nous
    login, this:

    - does NOT change the inference provider (``model.provider`` is untouched),
    - does NOT run model selection, and
    - does NOT offer the bulk "enable for all tools" Tool Gateway prompt.

    It only performs the Nous Portal device-code OAuth (when the user isn't
    already logged in) and refreshes entitlement, so the caller can enable the
    single tool the user picked.

    Entitlement is satisfied by paid service access OR a live free tool pool.
    When ``coverage_category`` is given (e.g. ``"fal"`` for image gen), the pool
    must cover that category specifically — so a pool user selecting video
    (``"fal-video"``, not pool-funded) is correctly denied.

    Returns ``True`` when the account is entitled after the flow, ``False``
    otherwise (declined login, login failed, or no entitlement).
    """

    def _entitled(account) -> bool:
        if account is None:
            return False
        if coverage_category is not None:
            return account.tool_gateway_entitled_for(coverage_category)
        return account.tool_gateway_entitled

    # Fast path: already entitled.
    try:
        info = get_nous_portal_account_info(force_fresh=True)
    except Exception:
        info = None
    if _entitled(info):
        return True

    # If not logged in at all, run the device-code login (auth only).
    if info is None or not info.logged_in:
        if not _run_nous_portal_login_only(capability=capability):
            return False
        try:
            info = get_nous_portal_account_info(force_fresh=True)
        except Exception:
            info = None

    if _entitled(info):
        return True

    # Logged in but not entitled for this capability — surface neutral billing
    # guidance, do not enable. coverage_category keeps a pool user who lacks this
    # one category from being told their credits are exhausted.
    message = format_nous_portal_entitlement_message(
        info, capability=capability, coverage_category=coverage_category
    )
    if message:
        for line in message.splitlines():
            print(f"  {line}")
    return False



def _run_nous_portal_login_only(*, capability: str) -> bool:
    """Run the Nous Portal device-code OAuth and persist credentials only.

    No model selection, no provider switch, no Tool Gateway bulk prompt.
    Returns ``True`` on a successful login, ``False`` if the user declined or
    the flow failed.
    """
    try:
        from drugmate_cli.auth import (
            _auth_store_lock,
            _load_auth_store,
            _nous_device_code_login,
            _read_shared_nous_state,
            _save_auth_store,
            _save_provider_state,
            _sync_nous_pool_from_auth_store,
            _try_import_shared_nous_state,
            _write_shared_nous_state,
        )
    except Exception as exc:  # pragma: no cover - defensive
        print(f"  Could not start Nous Portal login: {exc}")
        return False

    print()
    print(f"  {capability} requires a Nous Portal login.")
    try:
        proceed = input("  Log in to Nous Portal now? [Y/n]: ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print()
        return False
    if proceed not in {"", "y", "yes"}:
        print("  Skipped Nous Portal login.")
        return False

    try:
        # Snapshot the active_provider so a tool-config login never silently
        # switches the user's inference provider to Nous.
        with _auth_store_lock():
            prior_active_provider = _load_auth_store().get("active_provider")

        auth_state = None
        shared = _read_shared_nous_state()
        if shared:
            try:
                do_import = input(
                    "  Found existing Nous OAuth credentials. Import them? [Y/n]: "
                ).strip().lower()
            except (EOFError, KeyboardInterrupt):
                do_import = "y"
            if do_import in {"", "y", "yes"}:
                auth_state = _try_import_shared_nous_state(timeout_seconds=15.0)

        if auth_state is None:
            auth_state = _nous_device_code_login()

        with _auth_store_lock():
            auth_store = _load_auth_store()
            _save_provider_state(auth_store, "nous", auth_state)
            # Preserve the user's existing inference provider — this login is
            # for tool entitlement only, not a provider switch.
            if prior_active_provider:
                auth_store["active_provider"] = prior_active_provider
            else:
                auth_store.pop("active_provider", None)
            _save_auth_store(auth_store)

        _write_shared_nous_state(auth_state)
        _sync_nous_pool_from_auth_store()
        print("  Nous Portal login successful.")
        return True
    except KeyboardInterrupt:
        print("\n  Login cancelled.")
        return False
    except SystemExit:
        # _nous_device_code_login raises SystemExit on subscription_required;
        # it already printed billing guidance.
        return False
    except Exception as exc:
        print(f"  Nous Portal login failed: {exc}")
        return False



def get_nous_subscription_features(config: Optional[Dict[str, object]] = None, *, force_fresh: bool = False) -> NousSubscriptionFeatures:
    config = dict(config if config is not None else load_config() or {})
    provider_is_nous = str(_model_config_dict(config).get("provider") or "").lower() == "nous"
    try:
        account = get_nous_portal_account_info(force_fresh=force_fresh)
    except Exception:
        account = None
    logged_in = bool(account and account.logged_in)
    entitled = bool(logged_in and account.tool_gateway_entitled_for("modal"))
    terminal = config.get("terminal")
    terminal = terminal if isinstance(terminal, dict) else {}
    backend = str(terminal.get("backend") or "local").lower()
    mode = normalize_modal_mode(terminal.get("modal_mode"))
    enabled = _toolset_enabled(config, "terminal")
    direct = has_direct_modal_credentials()
    managed = bool(entitled and is_managed_tool_gateway_ready("modal"))
    state = resolve_modal_backend_state(mode, has_direct=direct, managed_ready=managed, managed_enabled=entitled)
    selected = state["selected_backend"]
    available = backend != "modal" or selected in {"direct", "managed"}
    feature = NousFeatureState(
        key="modal", label="Terminal execution", included_by_default=False,
        available=available, active=bool(enabled and available),
        managed_by_nous=bool(enabled and backend == "modal" and selected == "managed"),
        direct_override=bool(enabled and backend == "modal" and selected == "direct"),
        toolset_enabled=enabled, current_provider=backend,
        explicit_configured=backend == "modal",
    )
    return NousSubscriptionFeatures(provider_is_nous or logged_in, logged_in, provider_is_nous, {"modal": feature}, account)


def apply_nous_managed_defaults(config: Dict[str, object], *, enabled_toolsets: Optional[Iterable[str]] = None, force_fresh: bool = False) -> set[str]:
    """Model setup does not select or install execution backends."""
    return set()


def prompt_enable_tool_gateway(config: Dict[str, object], *, force_fresh: bool = True) -> set[str]:
    """Retain the model setup hook without offering removed web/media tools."""
    return set()

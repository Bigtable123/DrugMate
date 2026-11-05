# Drugmate development guide

This source checkout contains a Hermes-derived SBDD research agent. Start with
README.md, docs/tools.md and docs/architecture.md.

This file governs repository maintenance. The deployed research assistant's
identity, task scope and experiment behavior are defined in `SOUL.md` and the
scientific skills; do not duplicate its persona here. The assistant must not
self-modify those protocols while conducting an experiment. Maintainer-requested
source edits and translations are separate repository work. The cleanup-only
restrictions below do not disable validation or execution required by an
explicitly requested scientific experiment under its skill protocol.

- Keep skills/bio_research/ at its established path. Its pipeline entry, stage
  policies, authorization gates and artifact contracts define scientific behavior.
- Keep root SOUL.md as the distributed research persona. Preserve its experiment
  boundaries when maintaining runtime prompts. Maintain the canonical SKILL.md
  and its supporting instructions in English without changing its
  tool IDs, thresholds, authorization labels or command contracts.
- Keep only the Python CLI, launched with python drugmate --cli. TUI, WebUI,
  dashboard, their gateways and frontend dependencies are excluded at the
  maintainer's request. Preserve scientific HTML report generation and its
  bundled molecular viewers; those are research artifacts, not interface servers.
- Keep the supporting agent runtime, model providers, session persistence,
  per-conversation prompt caching and message-role ordering. Shared runtime
  modules live in drugmate_runtime/, including constants, time, logging, state,
  utilities, model_tools and toolsets. They are runtime dependencies. Keep the
  drugmate launcher, cli.py and run_agent.py at the root for their entrypoint
  and source-path conventions.
- General web/browser tools, cron, social/multimedia integrations, external skill
  catalogs, MCP clients/servers, generic plugin loading, legacy installers and unrelated frontend templates are excluded.
- Tests and test-only dependencies/scripts are excluded at the maintainer's
  request. Use source, documentation, and isolated local smoke checks for release
  verification; scientific experiments require explicit authorization.
- Do not download scientific tools or execute experiments for source cleanup.
- Configure host-specific scientific paths through the deployment variables in
  skills/bio_research/pipeline-entry/references/deployment.md. Keep credentials
  in local Drugmate settings.
- Public documentation and catalogs use deployment variables or relative paths.
  Do not publish personal/research-host paths or internal design-source references.
  Preserve operating-system paths and runtime-resolved artifact paths needed for
  execution; these are not deployment examples.
- Never commit installed environments, node_modules, caches, session databases,
  private configuration or experiment outputs. Use a separate DRUGMATE_HOME.
- Large memory downloads may be placed under memory/bundles/, which is excluded
  from Git and Python packages. Set SBDD_MEMORY_BUNDLES_ROOT before launch;
  keep frozen input bytes unchanged and experiment outputs outside the checkout.
- Preserve licenses/LICENSE, licenses/THIRD_PARTY_NOTICES.md and component
  copyright/license files. The upstream MIT license does not override GPL,
  BSD, Apache or OFL components.

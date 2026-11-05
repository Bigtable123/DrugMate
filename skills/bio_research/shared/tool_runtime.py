"""Cross-environment tool discovery, probing, and execution.

Scientific run configs choose tools and parameters.  Runtime placement stays in
``pipeline_registry.TOOL_ENVIRONMENTS`` and each tool's ``runtime`` declaration.
Commands are always passed as argv lists and never through a shell.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import pipeline_registry as R
from process_lifecycle import run_command


class ToolRuntimeError(RuntimeError):
    """Raised when a registered tool cannot be located or executed."""

    def __init__(self, message: str, result: dict | None = None):
        super().__init__(message)
        self.result = result or {}


@dataclass(frozen=True)
class PreparedTool:
    """Reusable command placement; the captured environment is private, never report data."""

    command: str
    environment: str
    _process_env: dict[str, str] = field(repr=False)


def find_conda() -> str | None:
    """Find conda without assuming the pipeline itself runs in base."""
    candidates = [
        os.environ.get("CONDA_EXE"),
        shutil.which("conda"),
        Path.home() / "miniconda3" / "bin" / "conda",
        Path.home() / "anaconda3" / "bin" / "conda",
        Path.home() / "miniforge3" / "bin" / "conda",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        try:
            path = Path(candidate)
            if path.is_file():
                return str(path.resolve())
        except OSError:
            # An inaccessible common installation root must not make runtime
            # discovery itself crash; continue to the remaining candidates.
            continue
    return None


def environment_spec(runtime: dict) -> dict:
    """Return the referenced environment declaration, validating its id."""
    env_id = runtime.get("environment", "current")
    try:
        return R.TOOL_ENVIRONMENTS[env_id]
    except KeyError as exc:
        raise ToolRuntimeError(f"unknown registered tool environment {env_id!r}") from exc


def _resolved_existing_dir(value: Any) -> Path | None:
    """Resolve an existing directory without trusting a stale activation path."""
    if not value:
        return None
    try:
        path = Path(str(value)).resolve()
        return path if path.is_dir() else None
    except OSError:
        return None


def build_process_env(
    runtime: dict,
    argv: list[str],
    overrides: dict[str, str] | None,
) -> dict[str, str]:
    """Build a subprocess environment isolated from stale Conda activation state."""
    process_env = {**os.environ, **(overrides or {})}
    try:
        env_spec = environment_spec(runtime)
    except ToolRuntimeError:
        return process_env
    if env_spec.get("kind") != "conda":
        return process_env

    prefix = env_spec.get("prefix")
    resolved_prefix = Path(str(prefix)).resolve() if prefix else None
    uses_conda_frontend = len(argv) >= 2 and argv[1] == "run"
    if uses_conda_frontend:
        # A long-lived IDE/agent process can retain CONDA_* values for a base
        # installation that has since been moved or removed.  Passing those
        # values to the new frontend makes ``conda run`` generate an activation
        # wrapper that calls the stale CONDA_EXE.  Start the frontend from a
        # neutral activation state while preserving unrelated tool variables.
        for key in tuple(process_env):
            if key in {
                "CONDA_DEFAULT_ENV",
                "CONDA_EXE",
                "CONDA_PREFIX",
                "CONDA_PROMPT_MODIFIER",
                "CONDA_PYTHON_EXE",
                "CONDA_SHLVL",
                "_CE_CONDA",
                "_CE_M",
            } or key.startswith("CONDA_PREFIX_"):
                process_env.pop(key, None)

        conda_bin = Path(argv[0]).resolve().parent
        conda_root = conda_bin.parent
        path_entries = [str(conda_bin)]
        condabin = conda_root / "condabin"
        if condabin.is_dir():
            path_entries.append(str(condabin))
        existing_path = process_env.get("PATH", "")
        if existing_path:
            path_entries.append(existing_path)
        process_env["PATH"] = os.pathsep.join(path_entries)
    elif resolved_prefix is not None:
        # When a prefix-registered environment is already active, or when the
        # Conda frontend is unavailable and build_command selected a binary
        # inside that prefix directly, expose the target prefix consistently.
        process_env["CONDA_PREFIX"] = str(resolved_prefix)
        process_env["CONDA_DEFAULT_ENV"] = resolved_prefix.name
        process_env["PATH"] = os.pathsep.join(
            [str(resolved_prefix / "bin"), process_env.get("PATH", "")]
        ).rstrip(os.pathsep)

    declared_env = env_spec.get("env") or {}
    if not isinstance(declared_env, dict):
        raise ToolRuntimeError("registered environment 'env' must be an object")
    for key, value in declared_env.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ToolRuntimeError(
                "registered environment variables must use string keys and values"
            )
        if "{prefix}" in value:
            if resolved_prefix is None:
                raise ToolRuntimeError(
                    f"registered environment variable {key!r} requires a prefix"
                )
            value = value.replace("{prefix}", str(resolved_prefix))
        process_env[key] = value
    return process_env


def build_command(runtime: dict, args: Iterable[Any] = (), *, command: str | None = None) -> list[str]:
    """Build an argv list for a registered command in its declared environment."""
    executable = command or runtime.get("command")
    if not executable:
        raise ToolRuntimeError("command runtime is missing its executable name")
    argv = [str(executable), *(str(arg) for arg in args)]
    env = environment_spec(runtime)
    kind = env.get("kind")
    if kind == "current":
        return argv
    if kind == "conda":
        name = env.get("name")
        prefix = env.get("prefix")
        if bool(name) == bool(prefix):
            raise ToolRuntimeError(
                "conda environment must declare exactly one of 'name' or 'prefix'"
            )
        active_env = os.environ.get("CONDA_DEFAULT_ENV", "")
        active_prefix = os.environ.get("CONDA_PREFIX", "")
        active_prefix_path = _resolved_existing_dir(active_prefix)
        if (
            name
            and active_prefix_path is not None
            and active_prefix_path.name == str(name)
            and (active_env == str(name) or Path(active_env).name == str(name))
        ):
            return argv
        if prefix:
            resolved_prefix = Path(str(prefix)).resolve()
            if active_prefix_path == resolved_prefix:
                return argv
            if not resolved_prefix.is_dir():
                raise ToolRuntimeError(
                    f"registered Conda prefix does not exist: {resolved_prefix}"
                )
        conda = find_conda()
        if conda:
            selector = ["-n", str(name)] if name else ["-p", str(resolved_prefix)]
            return [conda, "run", "--no-capture-output", *selector, *argv]

        # Prefix-registered environments remain executable on stripped-down
        # worker hosts where the environment was copied into place but the
        # Conda frontend itself is not mounted.  Resolve only the declared
        # executable inside that exact prefix; named environments still need a
        # Conda frontend because their installation path is intentionally not
        # guessed.
        if prefix:
            executable_path = Path(str(executable))
            if not executable_path.is_absolute():
                executable_path = resolved_prefix / "bin" / executable_path
            try:
                executable_path = executable_path.resolve()
                executable_path.relative_to(resolved_prefix)
            except (OSError, ValueError) as exc:
                raise ToolRuntimeError(
                    "prefix runtime executable must resolve inside its "
                    f"registered environment: {executable_path}"
                ) from exc
            if executable_path.is_file() and os.access(executable_path, os.X_OK):
                return [str(executable_path), *(str(arg) for arg in args)]
            raise ToolRuntimeError(
                "conda executable was not found and the registered prefix "
                f"contains no executable {executable_path}"
            )
        raise ToolRuntimeError(
            "conda executable was not found (checked CONDA_EXE, PATH, common roots)"
        )
    raise ToolRuntimeError(
        f"environment kind {kind!r} cannot execute local commands"
    )


def run_tool(
    runtime: dict,
    args: Iterable[Any] = (),
    *,
    command: str | None = None,
    cwd: str | Path | None = None,
    timeout: int = 300,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> dict:
    """Run one registered command and return a reproducibility record."""
    argv = build_command(runtime, args, command=command)
    process_env = build_process_env(runtime, argv, env)
    return _run_tool_command(
        argv, process_env, runtime.get("environment", "current"),
        cwd=cwd, timeout=timeout, check=check,
    )


def prepare_tool(
    runtime: dict,
    *,
    env: dict[str, str] | None = None,
    timeout: int = 60,
) -> PreparedTool:
    """Activate a registered runtime once, retaining hooks and env for repeated calls.

    The activation probe contains the full child environment (possibly secrets),
    so neither its output nor its exception details may enter a report or log.
    """
    argv = build_command(runtime)
    process_env = build_process_env(runtime, argv, env)
    environment = runtime.get("environment", "current")
    if len(argv) >= 2 and argv[1] == "run":
        marker = "SBDD_PREPARED_TOOL="
        code = (
            "import json,os,shutil,sys;"
            f"print({marker!r}+json.dumps({{'command':shutil.which(sys.argv[1]),"
            "'env':dict(os.environ)}))"
        )
        probe_argv = [*argv[:-1], "python", "-c", code, str(runtime["command"])]
        try:
            probe = run_command(probe_argv, env=process_env, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise ToolRuntimeError(
                f"tool environment preparation timed out after {timeout}s"
            ) from None
        except OSError:
            raise ToolRuntimeError("could not start tool environment preparation") from None
        except BaseException as exc:
            # Cancellation may carry partial captured output, including secrets.
            for attribute in ("stdout", "stderr", "output"):
                if hasattr(exc, attribute):
                    setattr(exc, attribute, "")
            raise
        if probe.returncode != 0:
            raise ToolRuntimeError(
                f"tool environment preparation failed (exit {probe.returncode})"
            )
        try:
            payload = next(line[len(marker):] for line in reversed(probe.stdout.splitlines())
                           if line.startswith(marker))
            captured = json.loads(payload)
            executable = captured["command"]
            captured_env = captured["env"]
            if not isinstance(captured_env, dict) or not all(
                isinstance(key, str) and isinstance(value, str)
                for key, value in captured_env.items()
            ):
                raise ValueError
            if executable is not None and not isinstance(executable, str):
                raise ValueError
        except (StopIteration, KeyError, TypeError, ValueError):
            raise ToolRuntimeError("tool environment preparation returned invalid data") from None
        # The preparation subprocess owns a temporary scope. Repeated commands
        # inherit our caller's scope, never that already-finished probe's scope.
        scope_key = "SBDD_PROCESS_SCOPES"
        if scope_key in process_env:
            captured_env[scope_key] = process_env[scope_key]
        else:
            captured_env.pop(scope_key, None)
        process_env = captured_env
    else:
        executable = shutil.which(argv[0], path=process_env.get("PATH"))
    if not executable:
        raise ToolRuntimeError(f"registered command {runtime['command']!r} was not found")
    # Do not dereference symlinks: some launchers select behavior from argv[0].
    return PreparedTool(os.path.abspath(executable), environment, process_env)


def run_prepared_tool(
    prepared: PreparedTool,
    args: Iterable[Any] = (),
    *,
    cwd: str | Path | None = None,
    timeout: int = 300,
    check: bool = True,
    cancel_event=None,
) -> dict:
    """Execute directly in a previously prepared environment, without reactivation."""
    return _run_tool_command(
        [prepared.command, *(str(arg) for arg in args)],
        prepared._process_env, prepared.environment,
        cwd=cwd, timeout=timeout, check=check, cancel_event=cancel_event,
    )


def _run_tool_command(
    argv, process_env, environment, *, cwd, timeout, check, cancel_event=None,
) -> dict:
    started = time.monotonic()
    try:
        options = {"cancel_event": cancel_event} if cancel_event is not None else {}
        process = run_command(argv, cwd=cwd, env=process_env, timeout=timeout, **options)
        result = {
            "success": process.returncode == 0,
            "command": argv,
            "returncode": process.returncode,
            "stdout": process.stdout or "",
            "stderr": process.stderr or "",
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "environment": environment,
        }
    except subprocess.TimeoutExpired as exc:
        def text_output(value):
            return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")

        result = {
            "success": False,
            "command": argv,
            "returncode": 124,
            "stdout": text_output(exc.stdout),
            "stderr": text_output(exc.stderr),
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "environment": environment,
            "error": f"tool timed out after {timeout}s",
            "process_group_terminated": True,
        }
    if check and not result["success"]:
        detail = result.get("error") or result["stderr"].strip() or result["stdout"].strip()
        raise ToolRuntimeError(detail or "tool command failed", result)
    return result


def probe_tool(runtime: dict, *, timeout: int = 30) -> dict:
    """Probe one runtime without changing scientific data or calling HTTP APIs."""
    try:
        timeout = max(1, int(runtime.get("probe_timeout_seconds", timeout)))
    except (TypeError, ValueError):
        timeout = max(1, int(timeout))
    kind = runtime.get("kind", "internal")
    base = {
        "runtime": runtime,
        "available": False,
    }
    try:
        if kind == "internal":
            return {**base, "available": True, "detail": "pipeline-internal implementation"}
        if kind == "composite":
            components = runtime.get("components")
            if not isinstance(components, dict) or not components:
                raise ToolRuntimeError(
                    "composite runtime must declare one or more components"
                )
            component_probes = []
            for label, component in components.items():
                if not isinstance(component, dict):
                    raise ToolRuntimeError(
                        f"composite runtime component {label!r} must be an object"
                    )
                component_probes.append({
                    "label": label,
                    **probe_tool(component, timeout=timeout),
                })
            available = all(
                component["available"] for component in component_probes
            )
            return {
                **base,
                "available": available,
                "detail": (
                    "all composite runtime components are available"
                    if available
                    else "one or more composite runtime components are unavailable"
                ),
                "component_probes": component_probes,
            }
        if kind == "http":
            env = environment_spec(runtime)
            if env.get("kind") != "http" or not env.get("base_url"):
                raise ToolRuntimeError("HTTP runtime has no registered base_url")
            return {
                **base,
                "available": True,
                "detail": "HTTP endpoint is registered; network is checked only when invoked",
                "base_url": env["base_url"],
            }
        if kind == "python_module":
            module = runtime.get("module")
            if not module:
                raise ToolRuntimeError("python_module runtime is missing module")
            env = environment_spec(runtime)
            if env.get("kind") == "current":
                available = importlib.util.find_spec(str(module)) is not None
                return {
                    **base, "available": available,
                    "detail": f"Python module {module!r} " + ("found" if available else "not found"),
                }
            probe_runtime = {
                "kind": "command",
                "environment": runtime.get("environment"),
                "command": "python",
            }
            code = f"import importlib.util,sys;sys.exit(0 if importlib.util.find_spec({module!r}) else 1)"
            result = run_tool(probe_runtime, ["-c", code], timeout=timeout, check=False)
            return {
                **base, "available": result["success"],
                "detail": f"Python module {module!r} " + ("found" if result["success"] else "not found"),
                "probe": result,
            }
        if kind == "command":
            if runtime.get("probe_code"):
                check_runtime = {
                    "kind": "command",
                    "environment": runtime.get("environment"),
                    "command": "python",
                }
                result = run_tool(
                    check_runtime, ["-c", runtime["probe_code"]],
                    timeout=timeout, check=False,
                )
                return {
                    **base, "available": result["success"],
                    "detail": (
                        "functional resource probe succeeded" if result["success"]
                        else runtime.get("probe_failure_hint", "functional resource probe failed")
                    ),
                    "probe": result,
                }
            commands = runtime.get("commands") or {}
            if commands:
                probes = []
                for label, executable in commands.items():
                    check_runtime = {
                        "kind": "command",
                        "environment": runtime.get("environment"),
                        "command": "python",
                    }
                    code = (
                        "import shutil,sys;"
                        f"sys.exit(0 if shutil.which({str(executable)!r}) else 1)"
                    )
                    command_probe = run_tool(
                        check_runtime, ["-c", code], timeout=timeout, check=False,
                    )
                    probes.append({
                        "label": label, "command": executable,
                        "available": command_probe["success"], "probe": command_probe,
                    })
                available = all(entry["available"] for entry in probes)
                return {
                    **base, "available": available,
                    "detail": "all composite commands found" if available else "one or more composite commands missing",
                    "command_probes": probes,
                }
            result = run_tool(
                runtime,
                runtime.get("probe_args", ["--help"]),
                timeout=timeout,
                check=False,
            )
            return {
                **base, "available": result["success"],
                "detail": (
                    "command probe succeeded" if result["success"]
                    else runtime.get("probe_failure_hint", "command probe failed")
                ),
                "probe": result,
            }
        raise ToolRuntimeError(f"unsupported tool runtime kind {kind!r}")
    except (OSError, ToolRuntimeError) as exc:
        result = getattr(exc, "result", None)
        return {**base, "available": False, "detail": str(exc), "probe": result}


def probe_plan(plan: dict, *, timeout: int = 30) -> dict:
    """Probe enabled runtimes, reporting record-only failures without blocking."""
    tools: list[dict] = []
    for step_id in plan.get("execution_order", []):
        for tool_id, tool in (plan["steps"][step_id].get("tools") or {}).items():
            if not tool.get("enabled"):
                continue
            runtime = tool.get("runtime")
            if not runtime:
                continue
            probe = probe_tool(runtime, timeout=timeout)
            tools.append({
                "step": step_id, "tool": tool_id,
                "required_for_preflight": runtime.get("required_for_preflight") is not False,
                **probe,
            })
    missing = [
        f"{entry['step']}.{entry['tool']}" for entry in tools
        if not entry["available"] and entry["required_for_preflight"]
    ]
    optional_missing = [
        f"{entry['step']}.{entry['tool']}" for entry in tools
        if not entry["available"] and not entry["required_for_preflight"]
    ]
    errors = [
        {"tool": f"{entry['step']}.{entry['tool']}", "detail": entry.get("detail")}
        for entry in tools if not entry["available"]
    ]
    return {
        "success": not missing,
        "tools": tools,
        "missing": missing,
        "optional_missing": optional_missing,
        "errors": errors,
        "conda_executable": find_conda(),
    }

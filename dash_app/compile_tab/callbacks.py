"""Callbacks for the compile tab."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import time

from dash import ALL, Input, Output, State, callback_context, dcc, html, no_update

from .discovery import (
    available_modules_after_stack,
    command_in_env,
    compiler_from_env,
    discover_compile_state,
    module_is_compiler_like,
    netcdf_modules_after_stack,
    resolve_lmod_stack,
)
from .runtime import (
    append_log_tail,
    build_compile_spec,
    build_warnings,
    cancel_compile_job,
    detect_build_statuses,
    finish_compile_job,
    format_runtime,
    job_process_is_live,
    poll_compile_job,
    read_log_increment,
    rebuild_failed_path_from_log,
    start_compile_job,
    start_rebuild_job,
    update_active_job,
)
from .state import BUILD_DIR, INSTALL_DIR


def options_from_flags(flag_values):
    """Convert checklist values to the compile options dict."""
    flags = set(flag_values or [])
    return {
        "debug": "debug" in flags,
        "run_tests": "run_tests" in flags,
        "python": "python" in flags,
        "fresh": "fresh" in flags,
        "openmp": "openmp" in flags,
        "tuning": "tuning" in flags,
        "gptl": "gptl" in flags,
        "disable_netcdf": "disable_netcdf" in flags,
        "disable_silhs": "disable_silhs" in flags,
    }


def collect_options(precision, gpu, toolchain, flag_values, extra_args, module_stack=None):
    """Collect all UI values into compile.py options."""
    options = options_from_flags(flag_values)
    options.update(
        {
            "precision": precision or "double",
            "gpu": gpu or "none",
            "toolchain": toolchain or "auto",
            "extra_args": extra_args or "",
            "module_stack": module_stack or [],
        }
    )
    return options


def lmod_available(discovery):
    """Return whether the current discovery payload has Lmod support."""
    return bool((discovery or {}).get("lmod", {}).get("available"))


def selected_module_stack(discovery, compiler_module, netcdf_module, extra_modules):
    """Build the active module stack for Lmod mode."""
    if not lmod_available(discovery):
        return []
    stack = []
    for module in [compiler_module, netcdf_module, *(extra_modules or [])]:
        if module and module not in stack:
            stack.append(module)
    return stack


def lmod_compiler_options(discovery):
    """Build compiler-module dropdown options."""
    options = []
    for item in (discovery or {}).get("lmod", {}).get("compiler_modules", []):
        label = item.get("module", "")
        canonical = item.get("canonical")
        if canonical:
            label = f"{label} ({canonical})"
        elif item.get("unverified"):
            label = f"{label} (unverified)"
        if item.get("enabled") is False:
            reason = item.get("disabled_reason") or "does not expose a supported compiler"
            label = f"{label} - {reason}"
        options.append(
            {
                "label": label,
                "value": item.get("module"),
                "disabled": item.get("enabled") is False,
            }
        )
    return options


def module_options(module_names, selected=()):
    """Build generic module dropdown options."""
    selected_set = set(selected or [])
    return [
        {"label": module, "value": module, "disabled": module in selected_set}
        for module in module_names
    ]


def lmod_compiler_module_names(discovery):
    """Return module names that should only be selectable as compilers."""
    names = set()
    for item in (discovery or {}).get("lmod", {}).get("compiler_modules", []):
        module = item.get("module")
        if module and not item.get("unverified"):
            names.add(module)
    return names


def first_enabled_value(options):
    """Return the first enabled dropdown value."""
    for option in options:
        if not option.get("disabled"):
            return option.get("value")
    return None


def default_netcdf_module(options):
    """Pick the most useful NetCDF module by default."""
    enabled = [option["value"] for option in options if not option.get("disabled")]
    if not enabled:
        return None
    for value in enabled:
        if "fortran" in value.lower():
            return value
    return enabled[0]


def render_lmod_stack_summary(discovery, compiler_module, netcdf_module, extra_modules):
    """Render loaded-stack diagnostics for Lmod mode."""
    if not lmod_available(discovery):
        return ""
    stack = selected_module_stack(discovery, compiler_module, netcdf_module, extra_modules)
    if not stack:
        return html.Div("Select a compiler module to inspect the resulting environment.", className="compile-muted")
    result = resolve_lmod_stack(stack)
    if not result.get("ok"):
        return html.Div(
            [
                html.Div("Module stack failed to load.", className="compile-warning"),
                html.Div(result.get("stderr", "").strip(), className="compile-path"),
            ]
        )
    env = result.get("env", {})
    rows = [
        ("Loaded", " ".join(stack)),
        ("Compiler", compiler_from_env(env) or "not detected"),
        ("FC", env.get("FC") or "not set"),
        ("CC", env.get("CC") or "not set"),
        ("nf-config", command_in_env("nf-config", env) or "not found"),
        ("nc-config", command_in_env("nc-config", env) or "not found"),
    ]
    children = [
        html.Div([html.Span(key, className="compile-key"), html.Span(value)])
        for key, value in rows
    ]
    stderr = (result.get("stderr") or "").strip()
    if stderr:
        children.append(html.Div(stderr.splitlines()[0], className="compile-warning"))
    return html.Div(children, className="compile-lmod-diagnostics")


def warning_is_blocking(warning):
    """Return whether a warning should prevent launching a compile."""
    return any(
        text in warning
        for text in [
            "requires NetCDF",
            "not currently available",
            "cmake was not found",
            "reject GPU",
            "Cannot load",
            "Select a compiler module",
            "No matching CLUBB toolchain",
        ]
    )


def environment_options(discovery):
    """Build dropdown options for environments."""
    options = []
    for environment in (discovery or {}).get("environments", []):
        label = environment.get("label", environment.get("id"))
        detail = environment.get("detail")
        if detail:
            label = f"{label} - {detail}"
        options.append(
            {
                "label": label,
                "value": environment.get("id"),
                "disabled": environment.get("enabled") is False,
            }
        )
    return options


def toolchain_options(discovery):
    """Build dropdown options for CMake toolchains."""
    options = [{"label": "Auto from environment", "value": "auto"}]
    for toolchain in (discovery or {}).get("toolchains", []):
        label = toolchain["name"]
        if not toolchain.get("matches_host"):
            label = f"{label} (other platform)"
        options.append({"label": label, "value": toolchain["path"]})
    return options


def render_detection(discovery):
    """Render compiler/toolchain discovery summary."""
    env = (discovery or {}).get("env", {})
    tools = (discovery or {}).get("tools", {})
    lmod = (discovery or {}).get("lmod", {})
    rows = [
        html.Div([html.Span("Host", className="compile-key"), html.Span((discovery or {}).get("host_platform", ""))]),
        html.Div([html.Span("cmake", className="compile-key"), html.Span(tools.get("cmake") or "not found")]),
        html.Div([html.Span("ninja", className="compile-key"), html.Span(tools.get("ninja") or "not found")]),
        html.Div([html.Span("FC", className="compile-key"), html.Span(env.get("FC") or "not set")]),
        html.Div([html.Span("LMOD", className="compile-key"), html.Span(lmod.get("message") or "not detected")]),
    ]
    return rows


def clicked_trigger_id():
    """Return the triggered id only when a button was actually clicked."""
    triggered_id = callback_context.triggered_id
    if not triggered_id:
        return None
    for item in callback_context.triggered:
        if not item.get("prop_id", "").endswith(".n_clicks"):
            continue
        value = item.get("value")
        if isinstance(value, (int, float)) and value > 0:
            return triggered_id
    return None


REBUILDABLE_BUILD_STATUSES = {"needs_configure", "needs_rebuild", "toolchain_dirty"}


def build_status_for(statuses, build_path):
    """Return status metadata for one build path."""
    return ((statuses or {}).get("statuses") or {}).get(
        build_path,
        {"status": "checking", "label": "checking", "detail": "Checking build freshness."},
    )


def render_build_status_badge(status_info):
    """Render a compact build freshness badge."""
    status = (status_info or {}).get("status", "checking")
    label = (status_info or {}).get("label") or status.replace("_", " ")
    class_name = "compile-badge"
    if status == "current":
        class_name += " compile-badge-good"
    elif status in REBUILDABLE_BUILD_STATUSES:
        class_name += " compile-badge-warn"
    elif status == "unknown":
        class_name += " compile-badge-bad"
    return html.Span(label, className=class_name, title=(status_info or {}).get("detail", ""))


def build_card_status_class(build, status, failed_rebuild_paths=None):
    """Return the visual state class for one build card."""
    if build.get("path") in (failed_rebuild_paths or set()):
        return "compile-build-card-failed"
    if not build.get("install_exists") or status == "unknown":
        return "compile-build-card-failed"
    if status == "current":
        return "compile-build-card-current"
    if status in REBUILDABLE_BUILD_STATUSES:
        return "compile-build-card-stale"
    return "compile-build-card-checking"


def render_build_status_summary(discovery, statuses):
    """Render a one-line summary of current/stale build status."""
    builds = (discovery or {}).get("builds", [])
    status_map = (statuses or {}).get("statuses") or {}
    if not builds:
        return ""
    if not status_map:
        return "Checking build freshness..."
    counts = {"current": 0, "stale": 0, "unknown": 0, "checking": 0}
    for build in builds:
        status = status_map.get(build.get("path"), {}).get("status", "checking")
        if status == "current":
            counts["current"] += 1
        elif status in REBUILDABLE_BUILD_STATUSES:
            counts["stale"] += 1
        elif status == "unknown":
            counts["unknown"] += 1
        else:
            counts["checking"] += 1
    parts = []
    if counts["stale"]:
        parts.append(f"{counts['stale']} stale")
    if counts["current"]:
        parts.append(f"{counts['current']} current")
    if counts["unknown"]:
        parts.append(f"{counts['unknown']} unknown")
    if counts["checking"]:
        parts.append(f"{counts['checking']} checking")
    return " | ".join(parts)


def completed_rebuild_paths_for_ui(job, statuses=None):
    """Return paths from a rebuild job that finished successfully outside store state."""
    if (job or {}).get("kind") != "rebuild":
        return set()
    if job.get("status") == "running":
        status = poll_compile_job(job)
        if status == 0:
            return set(job.get("build_paths") or [])
        return set()
    if job.get("returncode") == 0:
        checked_at = float((statuses or {}).get("checked_at") or 0)
        start_time = float(job.get("start_time") or 0)
        if checked_at < start_time:
            return set(job.get("build_paths") or [])
    return set()


def render_delete_confirmation(build):
    """Render the guarded delete confirmation row for one build."""
    return html.Div(
        [
            html.Div(f"Delete build {build.get('name')}?", className="compile-delete-title"),
            html.Div(
                [
                    html.Div("This removes the CMake build directory and its install directory.", className="compile-warning"),
                    html.Div(build.get("path", ""), className="compile-path"),
                    html.Div(build.get("install_prefix") or "No install prefix", className="compile-path"),
                ],
                className="compile-delete-body",
            ),
            html.Div(
                [
                    html.Button(
                        "Cancel",
                        id={"type": "compile-build-delete-cancel", "index": build["path"]},
                        type="button",
                        n_clicks=0,
                        className="compile-build-confirm-button",
                    ),
                    html.Button(
                        "Delete",
                        id={"type": "compile-build-delete-confirm", "index": build["path"]},
                        type="button",
                        n_clicks=0,
                        className="compile-build-confirm-button compile-build-confirm-delete",
                    ),
                ],
                className="compile-build-confirm-actions",
            ),
        ],
        className="compile-build-delete-confirm",
    )


def render_build_list(discovery, statuses=None, failures=None, delete_target=None, job=None):
    """Render existing builds as selectable rows with maintenance actions."""
    builds = (discovery or {}).get("builds", [])
    if not builds:
        return html.Div("No CMake build directories found.", className="compile-muted")
    job_running = job_process_is_live(job)
    rebuilding_paths = set((job or {}).get("build_paths") or []) if job_running else set()
    completed_rebuild_paths = completed_rebuild_paths_for_ui(job, statuses)
    failed_rebuild_paths = set((failures or {}).keys())
    if (job or {}).get("kind") == "rebuild" and (job or {}).get("status") == "failed":
        if job.get("failed_build_path"):
            failed_rebuild_paths.add(job["failed_build_path"])
        else:
            failed_path = rebuild_failed_path_from_log(job.get("log"))
            if failed_path:
                failed_rebuild_paths.add(failed_path)
    items = []
    for build in builds:
        if delete_target == build.get("path"):
            items.append(render_delete_confirmation(build))
            continue
        status_info = build_status_for(statuses, build.get("path"))
        if build.get("path") in completed_rebuild_paths:
            status_info = {
                **status_info,
                "status": "current",
                "label": "current",
                "detail": "Rebuild finished; freshness check pending.",
            }
        status = status_info.get("status", "checking")
        badges = []
        if build.get("is_latest"):
            badges.append(html.Span("latest", className="compile-badge compile-badge-good"))
        if build.get("is_selected"):
            badges.append(html.Span("selected", className="compile-badge compile-badge-selected"))
        badges.append(render_build_status_badge(status_info))
        if build.get("path") in rebuilding_paths:
            badges.append(html.Span("building", className="compile-badge compile-badge-warn"))
        if build.get("install_exists"):
            badges.append(html.Span("installed", className="compile-badge"))
        meta = " | ".join(
            part
            for part in [
                build.get("build_type"),
                f"precision {build.get('precision')}" if build.get("precision") else "",
                f"gpu {build.get('gpu')}" if build.get("gpu") else "",
                "openmp" if build.get("openmp") == "ON" else "",
                "python" if build.get("python") == "ON" else "",
            ]
            if part
        )
        card_class = " ".join(
            [
                "compile-build-card",
                "compile-build-card-button",
                build_card_status_class(build, status, failed_rebuild_paths),
                "compile-build-card-selected" if build.get("is_selected") else "",
            ]
        ).strip()
        row_children = [
            html.Button(
                [
                    html.Div(
                        [
                            html.Div([html.Span(build["name"], className="compile-build-name"), html.Span(badges, className="compile-badges")]),
                        ],
                        className="compile-build-card-header",
                    ),
                    html.Div(meta or "No cache metadata", className="compile-muted"),
                    html.Div(build.get("install_prefix") or "No install prefix", className="compile-path"),
                ],
                id={"type": "compile-build-select", "index": build["path"]},
                type="button",
                n_clicks=0,
                disabled=bool(job_running or build.get("is_selected") or not build.get("install_exists")),
                className=card_class,
                title="Use this build for run_scm.py" if build.get("install_exists") else "Install directory not found",
            ),
            html.Div(
                [
                    html.Button(
                        "↻",
                        id={"type": "compile-build-rebuild", "index": build["path"]},
                        type="button",
                        n_clicks=0,
                        disabled=job_running,
                        className="compile-build-row-button compile-build-rebuild-button",
                        title=f"Rebuild {build.get('name')}",
                    ),
                    html.Button(
                        "×",
                        id={"type": "compile-build-delete-request", "index": build["path"]},
                        type="button",
                        n_clicks=0,
                        disabled=job_running,
                        className="compile-build-row-button compile-build-delete-button",
                        title=f"Delete {build.get('name')}",
                    ),
                ],
                className="compile-build-actions",
            ),
        ]
        items.append(
            html.Div(
                row_children,
                className="compile-build-row",
            )
        )
    return items


def render_warnings(warnings):
    """Render warning lines."""
    if not warnings:
        return ""
    return [html.Div(warning, className="compile-warning") for warning in warnings]


def update_build_failures(failures, job, returncode, failed_path=None):
    """Persist failed rebuild markers until that build later rebuilds successfully."""
    updated = dict(failures or {})
    if (job or {}).get("kind") != "rebuild":
        return updated
    build_paths = [path for path in (job or {}).get("build_paths", []) if path]
    if returncode == 0:
        for path in build_paths:
            updated.pop(path, None)
        return updated
    if failed_path and failed_path in build_paths:
        for path in build_paths[:build_paths.index(failed_path)]:
            updated.pop(path, None)
    if failed_path:
        updated[failed_path] = {"returncode": returncode, "failed_at": time.time()}
    return updated


def set_selected_install(install_prefix):
    """Point install/selected at install_prefix."""
    target = Path(install_prefix).resolve()
    if not target.is_dir():
        raise RuntimeError(f"Install directory not found: {target}")
    selected_path = Path(INSTALL_DIR) / "selected"
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    selected_str = str(selected_path)
    if os.path.lexists(selected_str):
        if not selected_path.is_symlink():
            raise RuntimeError(f"{selected_path} exists and is not a symlink.")
        selected_path.unlink()
    os.symlink(str(target), selected_str)
    return selected_path


def find_build(discovery, build_path):
    """Return one discovered build by path."""
    return next(
        (item for item in (discovery or {}).get("builds", []) if item.get("path") == build_path),
        None,
    )


def discovery_with_selected_install(discovery, install_prefix):
    """Return discovery data with only selected-build flags updated."""
    try:
        selected_target = str(Path(install_prefix).resolve())
    except (OSError, TypeError):
        selected_target = ""

    updated = dict(discovery or {})
    updated_builds = []
    for build in (discovery or {}).get("builds", []):
        item = dict(build)
        try:
            item_target = str(Path(item.get("install_prefix") or "").resolve())
        except (OSError, TypeError):
            item_target = ""
        item["is_selected"] = bool(selected_target and item_target == selected_target)
        updated_builds.append(item)
    updated["builds"] = updated_builds
    return updated


def safe_child_path(path_value, root_value, label):
    """Resolve path_value and require it to be inside root_value."""
    root = Path(root_value).resolve()
    path = Path(path_value or "").resolve()
    if path == root or not path.is_relative_to(root):
        raise RuntimeError(f"{label} is outside {root}: {path}")
    return path


def unlink_install_aliases(install_path):
    """Remove selected/latest symlinks if they point at install_path."""
    removed = []
    target = Path(install_path).resolve()
    for name in ["selected", "latest"]:
        alias = Path(INSTALL_DIR) / name
        try:
            if alias.is_symlink() and alias.resolve() == target:
                alias.unlink()
                removed.append(name)
        except OSError:
            continue
    return removed


def delete_existing_build(build):
    """Delete one discovered build and its install directory when safe."""
    build_path = safe_child_path(build.get("path"), BUILD_DIR, "Build directory")
    install_prefix = build.get("install_prefix") or ""
    install_path = safe_child_path(install_prefix, INSTALL_DIR, "Install directory") if install_prefix else None
    removed = []

    if install_path and install_path.exists():
        unlink_install_aliases(install_path)
        if install_path.is_symlink():
            install_path.unlink()
        elif install_path.is_dir():
            shutil.rmtree(install_path)
        else:
            install_path.unlink()
        removed.append(str(install_path))

    if build_path.is_symlink():
        build_path.unlink()
    elif build_path.exists():
        shutil.rmtree(build_path)
    removed.append(str(build_path))
    return removed


def register_compile_callbacks(app):
    """Register compile tab callbacks."""

    @app.callback(
        Output("compile-discovery", "data"),
        Input("compile-refresh", "n_clicks"),
    )
    def refresh_discovery(_n_clicks):
        return discover_compile_state()

    @app.callback(
        Output("compile-discovery", "data", allow_duplicate=True),
        Output("compile-build-list", "children", allow_duplicate=True),
        Output("compile-build-status-summary", "children", allow_duplicate=True),
        Output("compile-build-select-message", "children"),
        Input({"type": "compile-build-select", "index": ALL}, "n_clicks"),
        State("compile-discovery", "data"),
        State("compile-build-statuses", "data"),
        State("compile-build-failures", "data"),
        State("compile-build-delete-target", "data"),
        State("compile-job", "data"),
        prevent_initial_call=True,
    )
    def select_existing_build(_n_clicks, discovery, statuses, failures, delete_target, job):
        triggered = clicked_trigger_id()
        if not triggered:
            return no_update, no_update, no_update, no_update
        build_path = triggered.get("index")
        build = find_build(discovery, build_path)
        if not build:
            return no_update, no_update, no_update, html.Div("Build no longer exists. Refresh and try again.", className="compile-warning")
        if not build.get("install_exists"):
            return no_update, no_update, no_update, html.Div("This build does not have an install directory.", className="compile-warning")
        try:
            set_selected_install(build.get("install_prefix", ""))
        except RuntimeError as exc:
            return no_update, no_update, no_update, html.Div(str(exc), className="compile-warning")
        updated_discovery = discovery_with_selected_install(discovery, build.get("install_prefix", ""))
        return (
            updated_discovery,
            render_build_list(updated_discovery, statuses, failures, delete_target, job),
            render_build_status_summary(updated_discovery, statuses),
            "",
        )

    @app.callback(
        Output("compile-build-delete-target", "data"),
        Input({"type": "compile-build-delete-request", "index": ALL}, "n_clicks"),
        Input({"type": "compile-build-delete-cancel", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def set_delete_target(_request_clicks, _cancel_clicks):
        triggered = clicked_trigger_id()
        if not triggered:
            return no_update
        if triggered.get("type") == "compile-build-delete-request":
            return triggered.get("index")
        return None

    @app.callback(
        Output("compile-discovery", "data", allow_duplicate=True),
        Output("compile-build-delete-target", "data", allow_duplicate=True),
        Output("compile-build-select-message", "children", allow_duplicate=True),
        Input({"type": "compile-build-delete-confirm", "index": ALL}, "n_clicks"),
        State("compile-discovery", "data"),
        State("compile-job", "data"),
        prevent_initial_call=True,
    )
    def confirm_delete_build(_confirm_clicks, discovery, job):
        triggered = clicked_trigger_id()
        if not triggered:
            return no_update, no_update, no_update
        if job_process_is_live(job):
            return no_update, no_update, html.Div("Cannot delete a build while a compile job is running.", className="compile-warning")
        build = find_build(discovery, triggered.get("index"))
        if not build:
            return no_update, None, html.Div("Build no longer exists. Refresh and try again.", className="compile-warning")
        try:
            delete_existing_build(build)
        except RuntimeError as exc:
            return no_update, no_update, html.Div(str(exc), className="compile-warning")
        return discover_compile_state(), None, html.Div(f"Deleted build: {build.get('name')}", className="compile-muted")

    @app.callback(
        Output("compile-job", "data", allow_duplicate=True),
        Output("compile-log", "data", allow_duplicate=True),
        Output("compile-log-offset", "data", allow_duplicate=True),
        Output("compile-interval", "disabled", allow_duplicate=True),
        Output("compile-interval", "n_intervals", allow_duplicate=True),
        Output("compile-build-select-message", "children", allow_duplicate=True),
        Input("compile-rebuild-all", "n_clicks"),
        Input({"type": "compile-build-rebuild", "index": ALL}, "n_clicks"),
        State("compile-discovery", "data"),
        prevent_initial_call=True,
    )
    def start_rebuild(_all_clicks, _rebuild_clicks, discovery):
        triggered = clicked_trigger_id()
        if not triggered:
            return no_update, no_update, no_update, no_update, no_update, no_update
        builds = (discovery or {}).get("builds", [])
        if triggered == "compile-rebuild-all":
            build_paths = [build.get("path") for build in builds if build.get("path")]
            label = "all builds"
            name = "all builds"
        else:
            build = find_build(discovery, triggered.get("index"))
            if not build:
                message = html.Div("Build no longer exists. Refresh and try again.", className="compile-warning")
                return no_update, no_update, no_update, no_update, no_update, message
            build_paths = [build.get("path")]
            label = build.get("name") or "selected build"
            name = label
        try:
            job = start_rebuild_job(build_paths, label)
        except RuntimeError as exc:
            return no_update, append_log_tail("", f"{exc}\n"), no_update, no_update, no_update, html.Div(str(exc), className="compile-warning")
        header = f"--- Running rebuild job ---\n{job['command']}\n\n"
        message = html.Div(f"Started rebuild: {name}", className="compile-muted")
        return job, header, 0, False, 0, message

    @app.callback(
        Output("compile-env-select", "options"),
        Output("compile-env-select", "value"),
        Output("compile-toolchain-select", "options"),
        Output("compile-toolchain-select", "value"),
        Output("compile-detection-summary", "children"),
        Output("compile-native-env-field", "style"),
        Output("compile-lmod-panel", "style"),
        Output("compile-lmod-compiler", "options"),
        Output("compile-lmod-compiler", "value"),
        Input("compile-discovery", "data"),
        State("compile-env-select", "value"),
        State("compile-toolchain-select", "value"),
        State("compile-lmod-compiler", "value"),
    )
    def update_discovery_views(discovery, selected_env, selected_toolchain, selected_lmod_compiler):
        env_options = environment_options(discovery)
        env_values = [option["value"] for option in env_options if not option.get("disabled")]
        env_value = selected_env if selected_env in env_values else (env_values[0] if env_values else None)
        tc_options = toolchain_options(discovery)
        tc_values = [option["value"] for option in tc_options]
        tc_value = selected_toolchain if selected_toolchain in tc_values else "auto"
        compiler_options = lmod_compiler_options(discovery)
        compiler_values = [option["value"] for option in compiler_options if not option.get("disabled")]
        compiler_value = selected_lmod_compiler if selected_lmod_compiler in compiler_values else first_enabled_value(compiler_options)
        native_style = {"display": "none"} if lmod_available(discovery) else {}
        lmod_style = {} if lmod_available(discovery) else {"display": "none"}
        return (
            env_options,
            env_value,
            tc_options,
            tc_value,
            render_detection(discovery),
            native_style,
            lmod_style,
            compiler_options,
            compiler_value,
        )

    @app.callback(
        Output("compile-build-statuses", "data"),
        Input("compile-discovery", "data"),
        Input("compile-build-status-interval", "n_intervals"),
        Input("compile-job", "data"),
    )
    def update_build_statuses(discovery, _tick, job):
        if not discovery:
            return {}
        if job_process_is_live(job):
            return no_update
        return detect_build_statuses(discovery)

    @app.callback(
        Output("compile-build-list", "children"),
        Output("compile-build-status-summary", "children"),
        Input("compile-discovery", "data"),
        Input("compile-build-statuses", "data"),
        Input("compile-build-failures", "data"),
        Input("compile-build-delete-target", "data"),
        Input("compile-job", "data"),
        Input("compile-interval", "n_intervals"),
    )
    def update_build_list(discovery, statuses, failures, delete_target, job, _compile_tick):
        return (
            render_build_list(discovery, statuses, failures, delete_target, job),
            render_build_status_summary(discovery, statuses),
        )

    @app.callback(
        Output("compile-rebuild-all", "disabled"),
        Input("compile-job", "data"),
        Input("compile-discovery", "data"),
        Input("compile-interval", "n_intervals"),
    )
    def update_rebuild_all_disabled(job, discovery, _compile_tick):
        return job_process_is_live(job) or not bool((discovery or {}).get("builds"))

    @app.callback(
        Output("compile-lmod-netcdf", "options"),
        Output("compile-lmod-netcdf", "value"),
        Input("compile-discovery", "data"),
        Input("compile-lmod-compiler", "value"),
        State("compile-lmod-netcdf", "value"),
    )
    def update_lmod_netcdf_options(discovery, compiler_module, selected_netcdf):
        if not lmod_available(discovery) or not compiler_module:
            return [], None
        netcdf_names = netcdf_modules_after_stack([compiler_module])
        netcdf_options = module_options(netcdf_names)
        netcdf_values = [option["value"] for option in netcdf_options if not option.get("disabled")]
        netcdf_value = selected_netcdf if selected_netcdf in netcdf_values else default_netcdf_module(netcdf_options)
        return netcdf_options, netcdf_value

    @app.callback(
        Output("compile-lmod-extra", "options"),
        Output("compile-lmod-extra", "value"),
        Input("compile-discovery", "data"),
        Input("compile-lmod-compiler", "value"),
        Input("compile-lmod-netcdf", "value"),
        State("compile-lmod-extra", "value"),
    )
    def update_lmod_extra_options(discovery, compiler_module, netcdf_value, selected_extra):
        if not lmod_available(discovery) or not compiler_module:
            return [], []
        stack_for_extra = [compiler_module] + ([netcdf_value] if netcdf_value else [])
        compiler_modules = lmod_compiler_module_names(discovery)
        extra_names = [
            module for module in available_modules_after_stack(stack_for_extra)
            if module not in stack_for_extra
            and module not in compiler_modules
            and not module_is_compiler_like(module)
        ]
        extra_value = [module for module in (selected_extra or []) if module in extra_names]
        return module_options(extra_names, stack_for_extra), extra_value

    @app.callback(
        Output("compile-lmod-stack-summary", "children"),
        Input("compile-discovery", "data"),
        Input("compile-lmod-compiler", "value"),
        Input("compile-lmod-netcdf", "value"),
        Input("compile-lmod-extra", "value"),
    )
    def update_lmod_stack_summary(discovery, compiler_module, netcdf_module, extra_modules):
        return render_lmod_stack_summary(discovery, compiler_module, netcdf_module, extra_modules)

    @app.callback(
        Output("compile-command-preview", "children"),
        Output("compile-command-copy", "content"),
        Output("compile-warnings", "children"),
        Input("compile-discovery", "data"),
        Input("compile-env-select", "value"),
        Input("compile-toolchain-select", "value"),
        Input("compile-precision-select", "value"),
        Input("compile-gpu-select", "value"),
        Input("compile-feature-flags", "value"),
        Input("compile-extra-args", "value"),
        Input("compile-lmod-compiler", "value"),
        Input("compile-lmod-netcdf", "value"),
        Input("compile-lmod-extra", "value"),
    )
    def update_command_preview(discovery, env_id, toolchain, precision, gpu, flags, extra_args, compiler_module, netcdf_module, extra_modules):
        module_stack = selected_module_stack(discovery, compiler_module, netcdf_module, extra_modules)
        options = collect_options(precision, gpu, toolchain, flags, extra_args, module_stack)
        spec = build_compile_spec(discovery, env_id, options)
        warnings = build_warnings(discovery, env_id, options)
        return spec["display"], spec["display"], render_warnings(warnings)

    @app.callback(
        Output("compile-job", "data", allow_duplicate=True),
        Output("compile-log", "data", allow_duplicate=True),
        Output("compile-log-offset", "data", allow_duplicate=True),
        Output("compile-interval", "disabled", allow_duplicate=True),
        Output("compile-interval", "n_intervals"),
        Input("compile-start", "n_clicks"),
        State("compile-discovery", "data"),
        State("compile-env-select", "value"),
        State("compile-toolchain-select", "value"),
        State("compile-precision-select", "value"),
        State("compile-gpu-select", "value"),
        State("compile-feature-flags", "value"),
        State("compile-extra-args", "value"),
        State("compile-lmod-compiler", "value"),
        State("compile-lmod-netcdf", "value"),
        State("compile-lmod-extra", "value"),
        prevent_initial_call=True,
    )
    def start_compile(_n_clicks, discovery, env_id, toolchain, precision, gpu, flags, extra_args, compiler_module, netcdf_module, extra_modules):
        module_stack = selected_module_stack(discovery, compiler_module, netcdf_module, extra_modules)
        options = collect_options(precision, gpu, toolchain, flags, extra_args, module_stack)
        warnings = build_warnings(discovery, env_id, options)
        if any(warning_is_blocking(warning) for warning in warnings):
            log = "Cannot start compile with the current selection:\n" + "\n".join(f"- {warning}" for warning in warnings) + "\n"
            return {}, log, 0, True, no_update
        try:
            job = start_compile_job(discovery, env_id, options)
        except RuntimeError as exc:
            return no_update, append_log_tail("", f"{exc}\n"), no_update, no_update, no_update
        header = f"--- Running compile job ---\n{job['command']}\n\n"
        return job, header, 0, False, 0

    @app.callback(
        Output("compile-job", "data", allow_duplicate=True),
        Output("compile-log", "data", allow_duplicate=True),
        Output("compile-interval", "disabled", allow_duplicate=True),
        Input("compile-cancel", "n_clicks"),
        State("compile-job", "data"),
        State("compile-log", "data"),
        prevent_initial_call=True,
    )
    def cancel_compile(_n_clicks, job, log_text):
        if not job or job.get("status") != "running":
            return no_update, no_update, no_update
        if not job_process_is_live(job):
            updated = finish_compile_job(job, 1)
            if log_text and not log_text.endswith("\n"):
                log_text += "\n"
            return updated, append_log_tail(log_text or "", "--- Job process is no longer running; marking it failed. ---\n"), True
        cancelled = cancel_compile_job(job)
        if not cancelled:
            return no_update, no_update, no_update
        runtime_txt = format_runtime(time.time() - float(job.get("start_time", time.time())))
        updated = dict(job)
        updated["status"] = "cancelled"
        updated["returncode"] = 1
        if log_text and not log_text.endswith("\n"):
            log_text += "\n"
        return updated, append_log_tail(log_text or "", f"--- Cancelled; runtime: {runtime_txt} ---\n"), True

    @app.callback(
        Output("compile-log", "data", allow_duplicate=True),
        Output("compile-job", "data", allow_duplicate=True),
        Output("compile-log-offset", "data", allow_duplicate=True),
        Output("compile-interval", "disabled", allow_duplicate=True),
        Output("compile-build-failures", "data"),
        Input("compile-interval", "n_intervals"),
        State("compile-job", "data"),
        State("compile-log", "data"),
        State("compile-log-offset", "data"),
        State("compile-build-failures", "data"),
        prevent_initial_call=True,
    )
    def poll_compile(_tick, job, log_text, offset, build_failures):
        if not job:
            return no_update, no_update, no_update, True, no_update
        chunk, new_offset = read_log_increment(job.get("log"), int(offset or 0))
        updated_log = append_log_tail(log_text or "", chunk)
        if chunk:
            job = update_active_job(job, {"last_output_time": time.time()})
        status = poll_compile_job(job)
        if status is None:
            return updated_log if chunk else no_update, job if chunk else no_update, new_offset, False, no_update
        lost_job = status == "lost"
        returncode = 1 if lost_job else status
        chunk, new_offset = read_log_increment(job.get("log"), int(new_offset or 0))
        if chunk:
            updated_log = append_log_tail(updated_log, chunk)
            job = update_active_job(job, {"last_output_time": time.time()})
        runtime_txt = format_runtime(time.time() - float(job.get("start_time", time.time())))
        updated_job = finish_compile_job(job, returncode)
        failed_path = None
        if job.get("kind") == "rebuild" and returncode != 0:
            failed_path = rebuild_failed_path_from_log(job.get("log"))
            if failed_path:
                updated_job = dict(updated_job)
                updated_job["failed_build_path"] = failed_path
        updated_failures = update_build_failures(build_failures, job, returncode, failed_path)
        job_label = "Rebuild" if job.get("kind") == "rebuild" else "Compile"
        label = "completed" if returncode == 0 else f"failed (exit {returncode})"
        if updated_log and not updated_log.endswith("\n"):
            updated_log += "\n"
        if lost_job:
            updated_log = append_log_tail(updated_log, "--- Job process is no longer running; marking it failed. ---\n")
        updated_log = append_log_tail(updated_log, f"--- {job_label} {label}; runtime: {runtime_txt} ---\n")
        return updated_log, updated_job, new_offset, True, updated_failures

    @app.callback(
        Output("compile-discovery", "data", allow_duplicate=True),
        Input("compile-job", "data"),
        prevent_initial_call=True,
    )
    def refresh_discovery_after_job(job):
        if not job or job.get("status") == "running":
            return no_update
        if job.get("kind") not in {"compile", "rebuild"}:
            return no_update
        return discover_compile_state()

    @app.callback(
        Output("compile-job", "data", allow_duplicate=True),
        Output("compile-log", "data", allow_duplicate=True),
        Output("compile-log-offset", "data", allow_duplicate=True),
        Output("compile-interval", "disabled", allow_duplicate=True),
        Input("compile-clear", "n_clicks"),
        State("compile-job", "data"),
        prevent_initial_call=True,
    )
    def clear_compile(_n_clicks, job):
        if job_process_is_live(job):
            return no_update, no_update, no_update, no_update
        return {}, "", 0, True

    @app.callback(
        Output("compile-console", "children"),
        Output("compile-job-summary", "children"),
        Input("compile-log", "data"),
        Input("compile-job", "data"),
        Input("compile-interval", "n_intervals"),
    )
    def render_console(log_text, job, _tick):
        if not job:
            summary = html.Div("No active compile job.", className="compile-muted")
            return log_text or "No compile jobs yet.", summary
        status = job.get("status", "running")
        if status == "running":
            now = time.time()
            runtime_txt = format_runtime(now - float(job.get("start_time", now)))
            quiet_txt = format_runtime(now - float(job.get("last_output_time") or job.get("start_time", now)))
            summary_text = f"running | pid {job.get('pid')} | elapsed {runtime_txt} | last output {quiet_txt} ago"
            status_class = "compile-status-running"
        elif status == "completed":
            summary_text = f"completed | exit 0"
            status_class = "compile-status-good"
        elif status == "cancelled":
            summary_text = "cancelled"
            status_class = "compile-status-bad"
        else:
            summary_text = f"failed | exit {job.get('returncode')}"
            status_class = "compile-status-bad"
        summary = html.Div(summary_text, className=f"compile-job-state {status_class}")
        return log_text or "Waiting for output...", summary

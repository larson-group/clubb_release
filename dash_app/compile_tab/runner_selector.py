"""Provide a shared runner chooser independent of execution-tab bodies.

Call register_runner_selector alongside build_selector.build_selector_overlay
when constructing the dashboard. The component owns build/runtime metadata and
publishes the existing compile-run-* preferences plus runner-selected-build.
Execution tabs consume these stores; launches still resolve install/selected.
environment_change lets the existing handoff heartbeat detect managed setup
completion from file timestamps and publish a readiness revision.

Read-only discovery and JAX inspection use the bounded native task transport's
runner lane once at page initialization. Compile publishes subsequent catalog
updates; managed environment changes refresh readiness using the same devices.
Explicit refresh reacquires both inventories. Opening the menu and changing
preferences never request discovery.
assets/35_runner_selector.js owns native click/change handling and paints cached
panels and badges before publishing preferences. It flushes pending values
before another control can submit a job; no selection mutation waits on Dash's
callback scheduler.
"""
from __future__ import annotations

import time
from pathlib import Path
from dash import ALL, Input, Output, State, html, no_update
from dash_app.plot_tab.async_callbacks import task_callback
from .build_selector import (
    BUILD_SELECTOR_TRIGGER_IDS, inspect_jax_runtime_profiles, jax_device_info,
    normalize_run_implementation, normalize_jax_profile,
    build_implementation_capability, selected_build_info,
    register_build_selector_position_callback,
)
from .callbacks import (
    clicked_trigger_id, render_build_selector_help, jax_preallocation_available,
    build_visual_state, render_compact_build_selector, selected_build_from_discovery,
    displayed_build_status, build_card_status_class, find_build, set_selected_install,
    discovery_with_selected_install,
)
from .discovery import discover_compile_state
from .runtime import detect_build_statuses


def environment_change(runtime_info, previous):
    """Detect completed runtime setup from managed files, without querying hardware."""
    signatures = []
    for info in (runtime_info or {}).values():
        venv = (info.get("runtime") or {}).get("venv")
        if not venv:
            continue
        for path in (Path(venv) / ".clubb-jax-requirements.sha256", Path(venv) / "bin" / "python"):
            try:
                stat = path.stat()
                signatures.append([str(path), stat.st_mtime_ns, stat.st_size])
            except OSError:
                signatures.append([str(path), None, None])
    if not signatures or signatures == (previous or {}).get("signature"):
        return no_update
    return {"signature": signatures, "revision": int((previous or {}).get("revision", -1)) + 1}



def _selector_state(
    anchor,
    discovery,
    statuses,
    failures,
    job,
    implementation,
    jax_profile,
    jax_runtime_info,
    jax_gpu,
    jax_xla_prealloc,
    _tab,
):
    """Render one cached selection and return its badge metadata."""
    job = (job or {}).get("compile") or {}
    job = dict(job.get("job") or {}) if job.get("state") in {"running", "stopping"} else {}
    implementation = normalize_run_implementation(implementation)
    jax_profile = normalize_jax_profile(jax_profile)
    jax_runtime_info = dict(jax_runtime_info or {})
    if jax_profile == "gpu" and jax_gpu:
        jax_runtime_info["gpu"] = jax_device_info(jax_runtime_info.get("gpu") or {}, jax_gpu)
    selected_info = selected_build_info() if implementation != "jax" else {}
    if selected_info.get("target"):
        discovery = discovery_with_selected_install(discovery, selected_info["target"])
    visual = build_visual_state(statuses, failures, job)
    menu = render_compact_build_selector(
        discovery,
        statuses,
        failures,
        job,
        visual,
        implementation,
        (anchor or {}).get("trigger_id", ""),
        jax_profile,
        jax_runtime_info,
        jax_gpu,
        jax_xla_prealloc,
    )
    selected_build = selected_build_from_discovery(discovery) if implementation != "jax" else None
    if implementation == "jax":
        available, reason = build_implementation_capability(
            "", "jax", jax_profile=jax_profile
        )
        runtime_info = dict((jax_runtime_info or {}).get(jax_profile) or {})
        runtime_status = str(runtime_info.get("status") or "checking")
        if runtime_info.get("selectable") is False:
            available = False
            reason = str(runtime_info.get("reason") or "JAX profile unavailable")
        selected_name = "CPU" if jax_profile == "cpu" else "GPU"
        if jax_profile == "gpu" and jax_gpu:
            gpus = (runtime_info.get("hardware") or {}).get("gpus") or []
            selected_gpu = next((gpu for gpu in gpus if gpu.get("uuid") == jax_gpu), {})
            selected_name = f"GPU {selected_gpu.get('index', '?')} · {selected_gpu.get('name', 'unavailable')}"
        if not available:
            status_class = "compile-build-card-failed"
            status_label = reason
        elif runtime_status == "ready":
            status_class = "compile-build-card-current"
            status_label = "ready"
        else:
            status_class = "compile-build-card-checking"
            status_label = runtime_status.replace("_", " ")
    elif selected_build:
        selected_name = selected_info["name"]
        status_info = displayed_build_status(selected_build, statuses, visual)
        status = status_info.get("status", "checking")
        status_class = build_card_status_class(selected_build, status, visual["failed_paths"])
        selected_path = selected_build.get("path")
        if selected_path == visual["active_path"]:
            status_class += " compile-build-card-rebuild-active"
        elif selected_path in visual["queued_paths"] or selected_path in visual["rebuilding_paths"]:
            status_class += " compile-build-card-rebuild-queued"
        status_label = status_info.get("label") or status.replace("_", " ")
    else:
        selected_name = selected_info["name"]
        status_class = (
            "compile-build-card-failed"
            if selected_info["status"] in {"missing", "broken"}
            else "compile-build-card-checking"
        )
        status_label = selected_info["status"]
    return menu, {
        "name": selected_name, "implementation": implementation,
        "className": "selected-build-badge " + status_class,
        "title": f"{implementation.title()} using {selected_name} — {status_label}. Click to configure.",
    }


def register_runner_selector(app):
    """Register the selector state owner without requiring any execution tab."""
    register_build_selector_position_callback(app)
    app.clientside_callback(
        "function(jobs, previous) { return window.clubbRunnerSelector.compileState(jobs, previous); }",
        Output("runner-compile-state", "data"), Input("dashboard-broker-jobs", "data"),
        State("runner-compile-state", "data"),
    )

    @app.callback(
        Output("compile-build-selector-help", "children"),
        Input({"type": "compile-selector-help-open", "index": ALL}, "n_clicks"),
        Input({"type": "compile-selector-help-close", "index": ALL}, "n_clicks"),
        State("compile-jax-runtime-info", "data"),
        prevent_initial_call=True,
    )
    def toggle_selector_help(_open, _close, runtime_info):
        """Open or dismiss cached runtime help after a real button click."""
        trigger = clicked_trigger_id()
        if not isinstance(trigger, dict):
            return no_update
        if trigger.get("type") == "compile-selector-help-close":
            return ""
        return render_build_selector_help(runtime_info)

    app.clientside_callback(
        "function(refresh, environment, previous) { return window.clubbRunnerSelector.runtimeRequest(refresh, environment, previous); }",
        Output("runner-runtime-request", "data"),
        Input("runner-refresh", "n_clicks"), Input("runner-environment-change", "data"),
        State("runner-runtime-request", "data"),
    )
    app.clientside_callback(
        "function(refresh, jobs, mounted, previous) { return window.clubbRunnerSelector.catalogRequest(refresh, jobs, mounted, previous); }",
        Output("runner-catalog-request", "data"), Input("runner-refresh", "n_clicks"),
        Input("runner-compile-state", "data"), Input("compile-tab-layout", "id", allow_optional=True),
        State("runner-catalog-request", "data"),
    )
    app.clientside_callback(
        "function(discovery) { window.clubbPlotTasks.cancel('runner_catalog'); return dash_clientside.no_update; }",
        Output("runner-catalog-signal", "data"), Input("compile-discovery", "data"),
    )

    @task_callback(app, "runner_info",
        Output("compile-jax-runtime-info", "data"),
        Input("runner-runtime-request", "data"), State("compile-jax-runtime-info", "data"),
    )
    def refresh_jax_runtime_info(request, current):
        """Discover devices once, or refresh readiness after a managed environment changes."""
        if not request:
            return no_update
        if request.get("kind") != "environment" or not current:
            return inspect_jax_runtime_profiles()
        from clubb_jax.runtime_info import refresh_environment
        from clubb_jax.run_jax import runtime_paths
        updated = {}
        for profile, info in current.items():
            accelerator = (info.get("runtime") or {}).get("accelerator")
            if not accelerator:
                updated[profile] = info
                continue
            requirements, _, _ = runtime_paths(accelerator)
            updated[profile] = refresh_environment(info, requirements)
        return updated

    @app.callback(
        Output("compile-build-selector-menu", "children"), Output("runner-badge-state", "data"),
        Input("compile-discovery", "data"), Input("compile-build-statuses", "data"),
        Input("runner-failures", "data"), Input("runner-compile-state", "data"),
        Input("compile-jax-runtime-info", "data"),
        State("compile-build-selector-anchor", "data"), State("compile-run-implementation", "data"),
        State("compile-run-jax-profile", "data"), State("compile-run-jax-gpu", "data"),
        State("compile-run-jax-xla-prealloc", "data"), State("dashboard-tabs", "value"),
    )
    def update_build_selector(discovery, statuses, failures, job, info, anchor, implementation, profile, gpu, prealloc, tab):
        """Build both cached panels when metadata changes; preferences are browser-local."""
        def state(implementation, profile="cpu", gpu=""):
            return _selector_state(None, discovery, statuses, failures, job,
                                   implementation, profile, info, gpu, prealloc, tab)
        fortran, fortran_badge = state("fortran")
        jax, _ = state("jax", profile, gpu)
        _, cpu_badge = state("jax")
        views = {"fortran": fortran_badge, "jax": {"cpu": cpu_badge}}
        for device in dict.fromkeys(["", gpu, *[item.get("uuid") for item in ((info or {}).get("gpu", {}).get("hardware") or {}).get("gpus", []) if item.get("uuid")]]):
            _, badge = state("jax", "gpu", device)
            badge["can_preallocate"] = jax_preallocation_available(info, "gpu", device)
            views["jax"][device or "gpu"] = badge
        return [fortran[0],
                html.Div(fortran[1:], id="runner-fortran-options"),
                html.Div(jax[1:], id="runner-jax-options", style={"display": "none"}),
                html.Div(id="runner-context-note", className="compile-run-implementation-note")], views

    @task_callback(app, "runner_catalog",
        Output("compile-discovery", "data", allow_duplicate=True),
        Output("compile-build-statuses", "data", allow_duplicate=True),
        Input("runner-catalog-request", "data"), prevent_initial_call=True,
    )
    def discover_runner(request):
        """Initialize or explicitly invalidate the shared build catalog in its worker."""
        if not request:
            return no_update, no_update
        discovery = discover_compile_state()
        return discovery, detect_build_statuses(discovery)

    @app.callback(
        Output("runner-selected-build", "data"),
        Output("runner-message", "children"),
        Output("compile-discovery", "data", allow_duplicate=True),
        Input({"type": "compile-build-selector-select", "index": ALL}, "n_clicks"),
        State("compile-discovery", "data"), prevent_initial_call=True,
    )
    def choose_build(clicks, discovery):
        """Select the installed build and publish its metadata to consumers."""
        trigger = clicked_trigger_id()
        if not isinstance(trigger, dict):
            return no_update, no_update, no_update
        build = find_build(discovery, trigger.get("index"))
        if not build or not build.get("install_exists") or build.get("install_prefix_mismatch"):
            return no_update, "Build unavailable. Refresh the runner list.", no_update
        try:
            set_selected_install(build["install_prefix"])
        except RuntimeError as exc:
            return no_update, str(exc), no_update
        return ({"build": build, "revision": time.time_ns()}, "",
                discovery_with_selected_install(discovery, build["install_prefix"]))

    @app.callback(
        Output("runner-message", "children", allow_duplicate=True),
        Input({"type": "compile-build-selector-rebuild", "index": ALL}, "n_clicks"),
        State("compile-discovery", "data"), prevent_initial_call=True,
    )
    def rebuild(clicks, discovery):
        """Submit the existing broker rebuild action for one discovered build."""
        trigger = clicked_trigger_id()
        if not isinstance(trigger, dict):
            return no_update
        build = find_build(discovery, trigger.get("index"))
        if not build:
            return "Build unavailable. Refresh the runner list."
        from dash_app.shared.broker_client import perform_action
        try:
            perform_action("launch_rebuild_request", {
                "builds": [build], "discovery": discovery, "label": build.get("name"),
            }, internal=True)
        except (OSError, RuntimeError, ValueError) as exc:
            return str(exc)
        return "Rebuild submitted. See Compile for progress."

    app.clientside_callback(
        "function(views, anchor, tab, implementation, profile, gpu, prealloc, menu, a, b, c) { return window.clubbRunnerSelector.synchronize(views, anchor, implementation, profile, gpu, prealloc); }",
        Output("runner-paint-signal", "data"), Input("runner-badge-state", "data"),
        Input("compile-build-selector-anchor", "data"), Input("dashboard-tabs", "value"),
        Input("compile-run-implementation", "data"), Input("compile-run-jax-profile", "data"),
        Input("compile-run-jax-gpu", "data"), Input("compile-run-jax-xla-prealloc", "data"),
        Input("compile-build-selector-menu", "children"),
        *[Input(component, "id", allow_optional=True) for component in BUILD_SELECTOR_TRIGGER_IDS],
    )

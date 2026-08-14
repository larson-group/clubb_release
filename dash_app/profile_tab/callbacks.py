"""Dash callbacks for Profile launch, storage, comparison, and plotting."""

from __future__ import annotations

import time
from datetime import datetime, timezone
from typing import Any

from dash import ALL, Input, Output, State, callback_context, dcc, html, no_update

from dash_app.shared.broker_client import perform_action

from .figures import (
    METRICS,
    decomposition_figure,
    profile_figure,
    relative_figure,
    timer_options,
    variability_figure,
)
from .library import (
    comparison_warnings,
    discover_profile_library,
    enrich_summary_rows,
    export_profile_library,
    import_profile_upload,
    load_profile_library_data,
    profile_option,
    resolve_library_path,
)
from .layout import active_profile_items, available_profile_buttons
from .runtime import (
    overwrite_confirmation,
    profile_command_display,
    profile_save_target,
    read_log_tail,
    read_profile_data,
    read_profile_results,
)


SETTING_STATE_IDS = (
    "profile-case",
    "profile-processes",
    "profile-columns",
    "profile-warmups",
    "profile-repeats",
    "profile-output",
    "profile-name",
    "profile-timeout",
    "profile-continue-on-error",
    "profile-config",
    "profile-override",
    "profile-executable",
    "profile-extra-args",
)


def collect_profile_settings(values: list[Any]) -> dict[str, Any]:
    mapped = dict(zip(SETTING_STATE_IDS, values))
    return {
        "case_name": mapped.get("profile-case"),
        "processes": mapped.get("profile-processes"),
        "columns": mapped.get("profile-columns"),
        "warmups": mapped.get("profile-warmups"),
        "repetitions": mapped.get("profile-repeats"),
        "output": mapped.get("profile-output"),
        "name": mapped.get("profile-name"),
        "timeout": mapped.get("profile-timeout"),
        "continue_on_error": "continue" in (mapped.get("profile-continue-on-error") or []),
        "config": mapped.get("profile-config"),
        "override": mapped.get("profile-override"),
        "executable": mapped.get("profile-executable"),
        "extra_args": mapped.get("profile-extra-args"),
    }


def _setting_states() -> list[State]:
    return [State(component_id, "value") for component_id in SETTING_STATE_IDS]


def _status_view(job: dict[str, Any], row_count: int) -> tuple[str, str, bool, bool]:
    state = str(job.get("state") or "idle")
    if state == "running":
        return (
            f"Running · {row_count} summary row{'s' if row_count != 1 else ''}",
            "profile-status profile-status-running",
            True,
            False,
        )
    if state == "stopping":
        return "Stopping…", "profile-status profile-status-running", True, True
    if state == "finished":
        return "Complete", "profile-status profile-status-success", False, True
    if state == "stopped":
        return "Cancelled", "profile-status profile-status-idle", False, True
    if state == "error":
        code = job.get("returncode")
        suffix = "" if code is None else f" (exit {code})"
        return f"Failed{suffix}", "profile-status profile-status-error", False, True
    return "Idle", "profile-status profile-status-idle", False, True


def reconcile_profile_selection(
    catalog: list[dict[str, Any]],
    selected: list[str] | None,
    baseline: str | None,
    detail: str | None,
    preferred: list[str] | None = None,
    replace_existing: str | None = None,
    select_default: bool = False,
) -> tuple[list[str], str | None, str | None]:
    """Preserve valid user choices while selecting newly created/imported runs."""
    valid = [str(record.get("run_id") or "") for record in catalog if record.get("run_id")]
    valid_set = set(valid)
    chosen = [str(value) for value in (selected or []) if str(value) in valid_set]
    baseline_replaced = False
    detail_replaced = False
    if replace_existing in valid_set:
        records = {
            str(record.get("run_id") or ""): record
            for record in catalog
            if record.get("run_id")
        }
        replacement = records[str(replace_existing)]
        identity = (
            str(replacement.get("label") or ""),
            str(replacement.get("case_name") or ""),
        )
        superseded = {
            run_id
            for run_id in chosen
            if run_id != replace_existing
            and (
                str(records[run_id].get("label") or ""),
                str(records[run_id].get("case_name") or ""),
            )
            == identity
        }
        baseline_replaced = baseline in superseded
        detail_replaced = detail in superseded
        chosen = [
            run_id
            for run_id in chosen
            if run_id not in superseded
        ]
    for run_id in preferred or []:
        if run_id in valid_set and run_id not in chosen:
            chosen.append(run_id)
    # Older releases created timestamped directories when a label was reused.
    # Treat equal label/case pairs as one logical save in the active selection,
    # keeping the newest catalog entry without deleting the older data.
    chosen_set = set(chosen)
    newest_by_identity: dict[tuple[str, str], str] = {}
    identity_by_run: dict[str, tuple[str, str]] = {}
    for record in catalog:
        run_id = str(record.get("run_id") or "")
        if run_id not in chosen_set:
            continue
        identity = (str(record.get("label") or ""), str(record.get("case_name") or ""))
        identity_by_run[run_id] = identity
        newest_by_identity.setdefault(identity, run_id)
    collapsed: list[str] = []
    collapsed_to: dict[str, str] = {}
    for run_id in chosen:
        survivor = newest_by_identity.get(identity_by_run.get(run_id, ()), run_id)
        collapsed_to[run_id] = survivor
        if survivor not in collapsed:
            collapsed.append(survivor)
    chosen = collapsed
    if not chosen and valid and select_default:
        chosen = [valid[0]]
    chosen_set = set(chosen)
    collapsed_baseline = collapsed_to.get(str(baseline or ""), baseline)
    collapsed_detail = collapsed_to.get(str(detail or ""), detail)
    baseline_value = (
        replace_existing
        if baseline_replaced
        else collapsed_baseline
        if collapsed_baseline in chosen_set
        else chosen[0]
        if chosen
        else None
    )
    detail_value = (
        replace_existing
        if detail_replaced
        else collapsed_detail
        if collapsed_detail in chosen_set
        else chosen[-1]
        if chosen
        else None
    )
    return chosen, baseline_value, detail_value


def profile_selection_preferences(
    job: dict[str, Any] | None,
    library_action: dict[str, Any] | None,
    *,
    library_action_triggered: bool,
) -> tuple[list[str], str | None]:
    """Return transient auto-selections without pinning completed profiles."""
    action = dict(library_action or {})
    preferred = (
        [str(value) for value in action.get("run_ids") or []]
        if library_action_triggered
        else []
    )
    current_job = dict(job or {})
    active_run_id = (
        str(current_job.get("run_id") or "")
        if str(current_job.get("state") or "") in {"running", "stopping"}
        else ""
    )
    if active_run_id:
        preferred.append(active_run_id)
    return preferred, active_run_id or None


def detail_process_options(
    rows: list[dict[str, Any]],
    detail_id: str | None,
    current: Any = None,
) -> tuple[list[dict[str, int]], int | None]:
    values = sorted(
        {
            int(str(row.get("process_count")))
            for row in rows
            if str(row.get("profile_id") or row.get("run_id") or "") == str(detail_id or "")
            and str(row.get("process_count") or "").isdigit()
        }
    )
    options = [{"label": str(value), "value": value} for value in values]
    try:
        selected = int(current)
    except (TypeError, ValueError):
        selected = None
    return options, selected if selected in values else (values[0] if values else None)


def load_profile_plot_data(
    output_value: Any,
    selected: list[str] | None,
    active_job: dict[str, Any] | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Load stored comparisons and replace the active profile with its live rows."""
    job = dict(active_job or {})
    is_active = str(job.get("state") or "") in {"running", "stopping"}
    known_active_id = str(job.get("run_id") or "") if is_active else ""
    stored_ids = [
        str(run_id)
        for run_id in (selected or [])
        if str(run_id) and str(run_id) != known_active_id
    ]
    summary_rows, process_rows = load_profile_library_data(output_value, stored_ids)
    if not is_active:
        return enrich_summary_rows(summary_rows), process_rows

    run_id, live_summary, live_processes = read_profile_data(job)
    if not run_id:
        return enrich_summary_rows(summary_rows), process_rows

    def belongs_to_active(row: dict[str, Any]) -> bool:
        return str(row.get("profile_id") or row.get("run_id") or "") == run_id

    # The active profile may already be selected in the library.  Replace that
    # disk snapshot instead of appending a duplicate, and include a newly
    # launched profile immediately even before its chooser pill is reconciled.
    summary_rows = [row for row in summary_rows if not belongs_to_active(row)]
    process_rows = [row for row in process_rows if not belongs_to_active(row)]
    summary_rows.extend(live_summary)
    process_rows.extend(live_processes)
    return enrich_summary_rows(summary_rows), process_rows


def profile_control_rows(
    summary_rows: list[dict[str, Any]],
    process_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    """Return only the small row subsets needed by timer/process controls."""
    timers = sorted(
        {
            str(row.get("timer_name") or "")
            for row in summary_rows
            if str(row.get("timer_name") or "")
        }
    )
    processes = sorted(
        {
            (
                str(row.get("profile_id") or row.get("run_id") or ""),
                int(row.get("process_count")),
            )
            for row in process_rows
            if str(row.get("profile_id") or row.get("run_id") or "")
            and str(row.get("process_count") or "").isdigit()
        }
    )
    return (
        [{"timer_name": timer} for timer in timers],
        [
            {"profile_id": profile_id, "process_count": process_count}
            for profile_id, process_count in processes
        ],
    )


def build_profile_graph_outputs(
    rows: list[dict[str, Any]],
    process_rows: list[dict[str, Any]],
    timer_name: Any,
    metric: Any,
    baseline_id: Any,
    detail_id: Any,
    process_count: Any,
    top_n: Any,
    percentage_values: Any,
    x_axis: Any,
    x_scale: Any,
    theme_name: Any,
    catalog: list[dict[str, Any]] | None,
    selected: list[str] | None,
) -> tuple[Any, Any, Any, Any, str, list[Any]]:
    """Build every Profile figure and its shared summary from one row snapshot."""
    metric = metric if metric in METRICS else "timer_max_seconds"
    try:
        top_n = max(1, min(20, int(top_n)))
    except (TypeError, ValueError):
        top_n = 8
    try:
        process_count = int(process_count) if process_count is not None else None
    except (TypeError, ValueError):
        process_count = None
    figures = (
        profile_figure(rows, timer_name, metric, theme_name, x_axis=x_axis, x_scale=x_scale),
        relative_figure(rows, timer_name, metric, baseline_id, theme_name, x_axis=x_axis, x_scale=x_scale),
        decomposition_figure(
            process_rows,
            detail_id,
            timer_name,
            process_count,
            top_n,
            "percentage" in (percentage_values or []),
            theme_name,
        ),
        variability_figure(process_rows, detail_id, timer_name, process_count, theme_name),
    )
    matching = [row for row in rows if str(row.get("timer_name") or "") == str(timer_name or "")]
    successful = sum(str(row.get("status") or "") == "success" for row in matching)
    selected_count = len(
        {
            str(row.get("profile_id") or row.get("run_id") or "")
            for row in rows
            if str(row.get("profile_id") or row.get("run_id") or "")
        }
    )
    summary = (
        "No measured results selected."
        if not rows
        else f"{selected_count} profile{'s' if selected_count != 1 else ''} · {successful} successful measurement{'s' if successful != 1 else ''} · {METRICS[metric]}"
    )
    warnings = comparison_warnings(list(catalog or []), list(selected or []))
    warning_view = [html.Span(warning, className="profile-comparison-warning") for warning in warnings]
    return *figures, summary, warning_view


def register_profile_callbacks(app) -> None:
    @app.callback(
        Output("profile-picker-expanded", "data"),
        Input("profile-picker-toggle", "n_clicks_timestamp"),
        State("profile-picker-expanded", "data"),
        prevent_initial_call=True,
    )
    def toggle_profile_picker(timestamp, expanded):
        if not isinstance(timestamp, (int, float)) or timestamp <= 0:
            return no_update
        return not bool(expanded)

    @app.callback(
        Output("profile-available-list", "children"),
        Output("profile-active-list", "children"),
        Output("profile-picker-menu", "className"),
        Input("profile-catalog", "data"),
        Input("profile-selected-runs", "value"),
        Input("profile-picker-expanded", "data"),
    )
    def render_profile_picker(catalog, selected, expanded):
        menu_class = "profile-picker-menu profile-picker-menu-expanded" if expanded else "profile-picker-menu"
        return (
            available_profile_buttons(catalog or [], selected or [], expanded=bool(expanded)),
            active_profile_items(catalog or [], selected or []),
            menu_class,
        )

    @app.callback(
        Output("profile-selected-runs", "value", allow_duplicate=True),
        Output("profile-baseline-run", "options", allow_duplicate=True),
        Output("profile-baseline-run", "value", allow_duplicate=True),
        Output("profile-detail-run", "options", allow_duplicate=True),
        Output("profile-detail-run", "value", allow_duplicate=True),
        Input({"type": "profile-add-run", "run_id": ALL}, "n_clicks_timestamp"),
        Input({"type": "profile-remove-run", "run_id": ALL}, "n_clicks_timestamp"),
        State("profile-catalog", "data"),
        State("profile-selected-runs", "value"),
        State("profile-baseline-run", "value"),
        State("profile-detail-run", "value"),
        prevent_initial_call=True,
    )
    def update_profile_selection(
        _add_timestamps,
        _remove_timestamps,
        catalog,
        current_selected,
        current_baseline,
        current_detail,
    ):
        trigger = callback_context.triggered_id
        if not isinstance(trigger, dict):
            return no_update, no_update, no_update, no_update, no_update
        triggered_value = (callback_context.triggered[0].get("value") if callback_context.triggered else None)
        if not isinstance(triggered_value, (int, float)) or triggered_value <= 0:
            return no_update, no_update, no_update, no_update, no_update

        valid = {
            str(record.get("run_id") or "")
            for record in (catalog or [])
            if str(record.get("run_id") or "")
        }
        selected = [str(run_id) for run_id in (current_selected or []) if str(run_id) in valid]
        run_id = str(trigger.get("run_id") or "")
        if trigger.get("type") == "profile-add-run" and run_id in valid and run_id not in selected:
            selected.append(run_id)
        elif trigger.get("type") == "profile-remove-run" and run_id in selected:
            selected = [value for value in selected if value != run_id]
        else:
            return no_update, no_update, no_update, no_update, no_update

        selected_set = set(selected)
        options = [profile_option(record) for record in (catalog or []) if str(record.get("run_id") or "") in selected_set]
        baseline = current_baseline if current_baseline in selected_set else (selected[0] if selected else None)
        detail = current_detail if current_detail in selected_set else (selected[-1] if selected else None)
        return selected, options, baseline, options, detail

    @app.callback(
        Output("profile-action", "data"),
        Output("profile-pending-run", "data"),
        Output("profile-overwrite-confirm", "displayed"),
        Output("profile-overwrite-confirm", "message"),
        Input("profile-start", "n_clicks"),
        Input("profile-cancel", "n_clicks"),
        Input("profile-overwrite-confirm", "submit_n_clicks"),
        Input("profile-overwrite-confirm", "cancel_n_clicks"),
        *_setting_states(),
        State("profile-pending-run", "data"),
        prevent_initial_call=True,
    )
    def run_profile_action(start_clicks, cancel_clicks, _confirm_clicks, _decline_clicks, *values):
        trigger = callback_context.triggered_id
        triggered_property = (
            callback_context.triggered[0].get("prop_id", "").rsplit(".", 1)[-1]
            if callback_context.triggered
            else ""
        )
        setting_values = values[:-1]
        pending = dict(values[-1] or {})
        if trigger == "profile-start" and start_clicks:
            try:
                settings = collect_profile_settings(list(setting_values))
                profile_command_display(settings)
                confirmation = overwrite_confirmation(settings)
                if confirmation:
                    return no_update, settings, True, confirmation
                response = perform_action("launch_profile_request", {"settings": settings}, internal=True)
                return (
                    {"kind": "started", "at": time.time(), "job": dict(response.get("job") or {})},
                    {},
                    False,
                    no_update,
                )
            except (OSError, RuntimeError, ValueError) as exc:
                return {"kind": "error", "at": time.time(), "message": str(exc)}, {}, False, no_update
        if trigger == "profile-overwrite-confirm" and triggered_property == "submit_n_clicks":
            if not pending:
                return no_update, {}, False, no_update
            try:
                settings = {**pending, "overwrite": True}
                response = perform_action("launch_profile_request", {"settings": settings}, internal=True)
                return (
                    {"kind": "started", "at": time.time(), "job": dict(response.get("job") or {})},
                    {},
                    False,
                    no_update,
                )
            except (OSError, RuntimeError, ValueError) as exc:
                return {"kind": "error", "at": time.time(), "message": str(exc)}, {}, False, no_update
        if trigger == "profile-overwrite-confirm" and triggered_property == "cancel_n_clicks":
            return no_update, {}, False, no_update
        if trigger == "profile-cancel" and cancel_clicks:
            try:
                perform_action("stop_profile", {}, internal=True)
                return {"kind": "stopping", "at": time.time()}, no_update, False, no_update
            except (OSError, RuntimeError, ValueError) as exc:
                return {"kind": "error", "at": time.time(), "message": str(exc)}, no_update, False, no_update
        return no_update, no_update, no_update, no_update

    @app.callback(
        Output("profile-name-status", "children"),
        Output("profile-name-status", "className"),
        Output("profile-name-control", "className"),
        Input("profile-name", "value"),
        Input("profile-case", "value"),
        Input("profile-output", "value"),
        Input("profile-catalog", "data"),
    )
    def update_profile_name_status(name, case_name, output, _catalog):
        target = profile_save_target({"name": name, "case_name": case_name, "output": output})
        if target["exists"] and target["replaceable"]:
            return (
                "Saved profile exists — running will ask before replacing it.",
                "profile-name-status profile-name-status-exists",
                "profile-name-control profile-name-control-exists",
            )
        if target["exists"]:
            return (
                "Name unavailable — this directory is not a timing profile.",
                "profile-name-status profile-name-status-conflict",
                "profile-name-control profile-name-control-conflict",
            )
        return (
            "Name available",
            "profile-name-status profile-name-status-available",
            "profile-name-control",
        )

    @app.callback(
        Output("profile-command-preview", "children"),
        *[Input(component_id, "value") for component_id in SETTING_STATE_IDS],
    )
    def update_command_preview(*setting_values):
        try:
            return profile_command_display(collect_profile_settings(list(setting_values)))
        except ValueError as exc:
            return f"Configuration incomplete: {exc}"

    @app.callback(
        Output("profile-job", "data"),
        Output("profile-active-results", "data"),
        Output("profile-status", "children"),
        Output("profile-status", "className"),
        Output("profile-log", "children"),
        Output("profile-start", "disabled"),
        Output("profile-cancel", "disabled"),
        Output("profile-graph", "figure"),
        Output("profile-relative-graph", "figure"),
        Output("profile-decomposition-graph", "figure"),
        Output("profile-variability-graph", "figure"),
        Output("profile-result-summary", "children"),
        Output("profile-comparison-warnings", "children"),
        Input("profile-interval", "n_intervals"),
        Input("profile-action", "data"),
        Input("dashboard-broker-jobs", "data"),
        State("profile-job", "data"),
        State("profile-selected-runs", "value"),
        State("profile-output", "value"),
        State("profile-result-timer", "value"),
        State("profile-result-metric", "value"),
        State("profile-baseline-run", "value"),
        State("profile-detail-run", "value"),
        State("profile-detail-process", "value"),
        State("profile-top-timers", "value"),
        State("profile-cost-percentage", "value"),
        State("profile-x-axis", "value"),
        State("profile-x-scale", "value"),
        State("theme-store", "data"),
        State("profile-catalog", "data"),
    )
    def refresh_profile(
        _ticks,
        action,
        broker_snapshot,
        current_job,
        selected,
        output_value,
        timer_name,
        metric,
        baseline_id,
        detail_id,
        process_count,
        top_n,
        percentage_values,
        x_axis,
        x_scale,
        theme_name,
        catalog,
    ):
        action = dict(action or {})
        broker_job = dict((broker_snapshot or {}).get("profile") or {})
        launch_job = dict(action.get("job") or {})
        if launch_job and float(launch_job.get("start_time") or 0.0) > float(broker_job.get("start_time") or 0.0):
            job = launch_job
        else:
            job = broker_job or launch_job or dict(current_job or {})
        run_id = ""
        rows: list[dict[str, Any]] = []
        if action.get("kind") == "error" and not broker_job:
            message = str(action.get("message") or "Profile action failed")
            status = message
            class_name = "profile-status profile-status-error"
            log_text = message
            start_disabled, cancel_disabled = False, True
        elif not job:
            job = {}
            status = "Idle"
            class_name = "profile-status profile-status-idle"
            log_text = "No profile runs yet."
            start_disabled, cancel_disabled = False, True
        else:
            run_id, rows = read_profile_results(job)
            if run_id and run_id != job.get("run_id"):
                job = {**job, "run_id": run_id, "result_rows": len(rows)}
            status, class_name, start_disabled, cancel_disabled = _status_view(job, len(rows))
            log_text = str(job.get("log_tail") or "") or read_log_tail(job.get("log"))
        progress = {
            "run_id": run_id,
            "row_count": len(rows),
            "tick": int(_ticks or 0),
        }
        try:
            plot_rows, process_rows = load_profile_plot_data(output_value, selected, job)
        except OSError:
            plot_rows, process_rows = [], []
        graph_outputs = build_profile_graph_outputs(
            plot_rows,
            process_rows,
            timer_name,
            metric,
            baseline_id,
            detail_id,
            process_count,
            top_n,
            percentage_values,
            x_axis,
            x_scale,
            theme_name,
            catalog,
            selected,
        )
        return (
            job,
            progress,
            status,
            class_name,
            log_text or "Waiting for timing output…",
            start_disabled,
            cancel_disabled,
            *graph_outputs,
        )

    @app.callback(
        Output("profile-library-action", "data"),
        Input("profile-library-refresh", "n_clicks"),
        Input("profile-library-import", "contents"),
        State("profile-library-import", "filename"),
        State("profile-output", "value"),
        prevent_initial_call=True,
    )
    def profile_library_action(_refresh_clicks, upload_contents, upload_name, output_value):
        if callback_context.triggered_id == "profile-library-import":
            try:
                imported = import_profile_upload(output_value, upload_name, upload_contents)
            except (OSError, ValueError) as exc:
                return {"kind": "error", "at": time.time(), "message": f"Import failed: {exc}"}
            return {
                "kind": "imported",
                "at": time.time(),
                "run_ids": imported,
                "message": f"Imported {len(imported)} profile{'s' if len(imported) != 1 else ''}.",
            }
        return {"kind": "refreshed", "at": time.time(), "message": "Profile library refreshed."}

    @app.callback(
        Output("profile-library-download", "data"),
        Output("profile-export-status", "children"),
        Input("profile-library-export", "n_clicks"),
        State("profile-selected-runs", "value"),
        State("profile-output", "value"),
        prevent_initial_call=True,
    )
    def export_selected_profiles(clicks, selected, output_value):
        if not clicks:
            return no_update, no_update
        try:
            data = export_profile_library(output_value, selected or [])
        except (OSError, ValueError) as exc:
            return no_update, f"Export failed: {exc}"
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        filename = f"clubb-timing-profiles-{stamp}.zip"
        return dcc.send_bytes(data, filename, type="application/zip"), f"Exported {len(selected or [])} profile{'s' if len(selected or []) != 1 else ''}."

    @app.callback(
        Output("profile-catalog", "data"),
        Output("profile-selected-runs", "options"),
        Output("profile-selected-runs", "value"),
        Output("profile-baseline-run", "options"),
        Output("profile-baseline-run", "value"),
        Output("profile-detail-run", "options"),
        Output("profile-detail-run", "value"),
        Output("profile-library-status", "children"),
        Output("profile-selection-library", "data"),
        Input("profile-output", "value"),
        Input("profile-library-action", "data"),
        Input("profile-job", "data"),
        Input("profile-active-results", "data"),
        State("profile-selected-runs", "value"),
        State("profile-baseline-run", "value"),
        State("profile-detail-run", "value"),
        State("profile-selection-library", "data"),
    )
    def refresh_profile_library(
        output_value,
        library_action,
        job,
        _active_rows,
        current_selected,
        current_baseline,
        current_detail,
        current_library,
    ):
        library_key = str(resolve_library_path(output_value))
        try:
            catalog = discover_profile_library(output_value)
        except OSError as exc:
            return [], [], [], [], None, [], None, f"Unable to read profile library: {exc}", library_key
        options = [profile_option(record) for record in catalog]
        action = dict(library_action or {})
        preferred, replace_existing = profile_selection_preferences(
            job,
            action,
            library_action_triggered=callback_context.triggered_id == "profile-library-action",
        )
        library_changed = str(current_library or "") != library_key
        selected, baseline, detail = reconcile_profile_selection(
            catalog,
            [] if library_changed else list(current_selected or []),
            current_baseline,
            current_detail,
            preferred,
            replace_existing,
            select_default=library_changed and not preferred,
        )
        selected_set = set(selected)
        selected_options = [option for option in options if option["value"] in selected_set]
        message = str(action.get("message") or "")
        if action.get("kind") == "error":
            status = message
        else:
            count = len(catalog)
            status = f"{count} stored profile{'s' if count != 1 else ''}"
            if message:
                status += f" · {message}"
        return (
            catalog,
            options,
            selected,
            selected_options,
            baseline,
            selected_options,
            detail,
            status,
            library_key,
        )

    @app.callback(
        Output("profile-results", "data"),
        Output("profile-process-results", "data"),
        Input("profile-selected-runs", "value"),
        Input("profile-output", "value"),
        Input("profile-interval", "n_intervals"),
        State("profile-job", "data"),
    )
    def load_selected_profile_data(selected, output_value, _ticks, active_job):
        try:
            summary_rows, process_rows = load_profile_plot_data(output_value, selected, active_job)
        except OSError:
            return [], []
        return profile_control_rows(summary_rows, process_rows)

    @app.callback(
        Output("profile-result-timer", "options"),
        Output("profile-result-timer", "value"),
        Input("profile-results", "data"),
        State("profile-result-timer", "value"),
    )
    def update_timer_choices(rows, selected_timer):
        default_timer = "advance_clubb_to_end"
        options = timer_options(rows or [], default_timer)
        values = [option["value"] for option in options]
        if selected_timer in values:
            return options, selected_timer
        if default_timer in values:
            return options, default_timer
        return options, (values[0] if values else None)

    @app.callback(
        Output("profile-detail-process", "options"),
        Output("profile-detail-process", "value"),
        Input("profile-process-results", "data"),
        Input("profile-detail-run", "value"),
        State("profile-detail-process", "value"),
    )
    def update_detail_processes(rows, detail_id, current):
        return detail_process_options(list(rows or []), detail_id, current)

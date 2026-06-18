"""Callbacks for starting, stopping, and polling dashboard tuning runs."""

from __future__ import annotations

from dash import ALL, Input, Output, State, ctx, html, no_update

from .runtime import (
    active_tuning_job,
    active_job_exited,
    clear_finished_job,
    empty_status_payload,
    poll_loss_runs,
    read_tuning_results,
    read_tuning_status,
    start_loss_run,
    start_tuning_job,
    stop_tuning_job,
)

from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_LOSS_MODE,
    LOSS_MODE_NAMES,
)


def _case_list(case_names):
    if isinstance(case_names, str):
        case_names = [case_names]
    return [str(case_name).strip() for case_name in (case_names or []) if str(case_name).strip()]


def _integer_seconds(raw_value):
    value = float(raw_value)
    if int(value) != value:
        raise ValueError("time values must be integer seconds")
    return int(value)


def _case_configs_from_rows(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
):
    configs = []
    for raw_name, raw_t0, raw_t1, raw_average in zip(
        case_names or [],
        time_start_values or [],
        time_end_values or [],
        average_time_values or [],
    ):
        case_name = "" if raw_name is None else str(raw_name).strip()
        if not case_name:
            continue
        time_start = _integer_seconds(raw_t0)
        time_end = _integer_seconds(raw_t1)
        average_time_seconds = _integer_seconds(raw_average)
        if time_end <= time_start:
            raise ValueError("time end must be after time start")
        if average_time_seconds < 1:
            raise ValueError("average_time_seconds must be >= 1")
        if (time_end - time_start) % average_time_seconds != 0:
            raise ValueError("average_time_seconds must divide the time range evenly")
        num_time_windows = (time_end - time_start) // average_time_seconds
        configs.append(
            {
                "case_name": case_name,
                "time_average_range": [time_start, time_end],
                "average_time_seconds": average_time_seconds,
                "num_time_windows": int(num_time_windows),
            }
        )
    return configs


def format_tuning_error(status, active_job):
    """Render a tuning failure with enough context to diagnose without opening JSON."""
    error_message = str((status or {}).get("error_message") or "Unknown tuning error.")
    log_path = (active_job or {}).get("log_path") or (status or {}).get("log_path")
    visible_lines = error_message.splitlines()[:80]
    if len(error_message.splitlines()) > len(visible_lines):
        visible_lines.append("... truncated in UI; see worker.log for the full traceback ...")

    children = [
        html.Div("Tuning job failed.", style={"fontWeight": "700", "marginBottom": "6px"}),
        html.Pre(
            "\n".join(visible_lines),
            style={
                "whiteSpace": "pre-wrap",
                "maxHeight": "360px",
                "overflowY": "auto",
                "padding": "8px",
                "border": "1px solid #fca5a5",
                "borderRadius": "4px",
                "backgroundColor": "#fff1f2",
                "color": "#7f1d1d",
            },
        ),
    ]
    if log_path:
        children.append(
            html.Div(
                f"Worker log: {log_path}",
                style={"marginTop": "6px", "fontFamily": "monospace", "fontSize": "12px"},
            )
        )
    return html.Div(children)


def build_validation_message(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
    selected_fields,
    param_values,
    min_values,
    max_values,
    batch_size,
    max_workers,
    strategy_mode,
    loss_mode,
    aggregation_mode,
    random_max_samples,
    resolve_spacing,
    case_data,
):
    """Validate the requested tuning configuration and return an error string if needed."""
    cases = _case_list(case_names)
    if not cases:
        return "Select at least one case."
    if len(set(cases)) != len(cases):
        return "Each tuning case can only be added once."
    for case_name in cases:
        if not (case_data or {}).get(case_name, {}):
            return f"No tuning defaults found for case {case_name}."
    try:
        case_configs = _case_configs_from_rows(
            case_names,
            time_start_values,
            time_end_values,
            average_time_values,
        )
    except TypeError:
        return "Each case row needs numeric time and average-time values."
    except ValueError as exc:
        message = str(exc)
        if "divide" in message:
            return "Average time must divide the selected time range evenly."
        return "Each case row needs integer time and average-time values."
    if len(case_configs) != len(cases):
        return "Each case row needs a case, time range, and average time."
    for config in case_configs:
        case_name = config["case_name"]
        time_start, time_end = config["time_average_range"]
        average_time_seconds = config["average_time_seconds"]
        if time_end <= time_start:
            return f"{case_name} requires time end > time start."
        if average_time_seconds < 1:
            return f"{case_name} requires average time >= 1 second."
        if (time_end - time_start) % average_time_seconds != 0:
            return f"{case_name} requires average time to divide the selected time range evenly."
    if not selected_fields:
        return "Select at least one CLUBB field to tune."
    for case_name in cases:
        case_fields = set(((case_data or {}).get(case_name, {}) or {}).get("clubb_fields", []))
        missing_fields = [field_name for field_name in selected_fields if field_name not in case_fields]
        if missing_fields:
            return f"Field(s) not available for {case_name}: {', '.join(missing_fields)}"
    if batch_size in (None, ""):
        return "Enter a batch size."
    try:
        batch_size_value = float(batch_size)
    except (TypeError, ValueError):
        return "Batch size must be an integer."
    if int(batch_size_value) != batch_size_value:
        return "Batch size must be an integer."
    batch_size = int(batch_size_value)
    if batch_size < 1:
        return "Batch size must be >= 1."
    if max_workers in (None, ""):
        return "Enter a max worker count."
    try:
        max_workers_value = float(max_workers)
    except (TypeError, ValueError):
        return "Max workers must be an integer."
    if int(max_workers_value) != max_workers_value:
        return "Max workers must be an integer."
    if int(max_workers_value) < 1:
        return "Max workers must be >= 1."
    if strategy_mode not in {"random", "resolve"}:
        return "Select a tuning mode."
    loss_mode = loss_mode or DEFAULT_LOSS_MODE
    aggregation_mode = aggregation_mode or DEFAULT_AGGREGATION_MODE
    if loss_mode not in LOSS_MODE_NAMES:
        return "Select a loss mode."
    if aggregation_mode not in AGGREGATION_MODE_NAMES:
        return "Select an aggregation mode."
    if strategy_mode == "random":
        if random_max_samples in (None, ""):
            return "Enter a random max sample count."
        try:
            max_samples_value = float(random_max_samples)
        except (TypeError, ValueError):
            return "Random max samples must be an integer."
        if int(max_samples_value) != max_samples_value:
            return "Random max samples must be an integer."
        if int(max_samples_value) < 1:
            return "Random max samples must be >= 1."
    if strategy_mode == "resolve":
        if resolve_spacing in (None, ""):
            return "Enter a resolve spacing."
        try:
            spacing_value = float(resolve_spacing)
        except (TypeError, ValueError):
            return "Resolve spacing must be numeric."
        if spacing_value <= 0.0:
            return "Resolve spacing must be > 0."

    ranges = []
    seen_names = set()
    for param_name, min_text, max_text in zip(param_values or [], min_values or [], max_values or []):
        param_name = "" if param_name is None else str(param_name).strip()
        min_text = "" if min_text is None else str(min_text).strip()
        max_text = "" if max_text is None else str(max_text).strip()
        if not any((param_name, min_text, max_text)):
            continue
        if not all((param_name, min_text, max_text)):
            return "Each parameter row must include parameter, min, and max."
        if param_name in seen_names:
            return f"Duplicate tuning parameter: {param_name}"
        seen_names.add(param_name)
        try:
            min_value = float(min_text.replace("D", "E").replace("d", "e"))
            max_value = float(max_text.replace("D", "E").replace("d", "e"))
        except ValueError:
            return f"Invalid numeric range for {param_name}."
        if min_value > max_value:
            return f"{param_name} requires min <= max."
        ranges.append({"name": param_name, "min": min_value, "max": max_value})

    if not ranges:
        return "Add at least one parameter range to tune."
    return ""


def build_request_payload(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
    selected_fields,
    param_values,
    min_values,
    max_values,
    batch_size,
    max_workers,
    strategy_mode,
    loss_mode,
    aggregation_mode,
    random_max_samples,
    resolve_spacing,
):
    """Build the worker request payload from the live UI values."""
    batch_size = int(float(batch_size))
    max_workers = int(float(max_workers))
    ranges = []
    for param_name, min_text, max_text in zip(param_values or [], min_values or [], max_values or []):
        param_name = "" if param_name is None else str(param_name).strip()
        min_text = "" if min_text is None else str(min_text).strip()
        max_text = "" if max_text is None else str(max_text).strip()
        if not all((param_name, min_text, max_text)):
            continue
        ranges.append(
            {
                "name": param_name,
                "min": float(min_text.replace("D", "E").replace("d", "e")),
                "max": float(max_text.replace("D", "E").replace("d", "e")),
            }
        )
    strategy_options = {}
    if strategy_mode == "random":
        strategy_options["max_samples"] = int(float(random_max_samples))
    if strategy_mode == "resolve":
        strategy_options["spacing"] = float(resolve_spacing)

    case_configs = _case_configs_from_rows(
        case_names,
        time_start_values,
        time_end_values,
        average_time_values,
    )

    return {
        "cases": [config["case_name"] for config in case_configs],
        "case_configs": case_configs,
        "selected_fields": list(selected_fields or []),
        "batch_size": int(batch_size),
        "max_workers": int(max_workers),
        "loss_mode": loss_mode or DEFAULT_LOSS_MODE,
        "aggregation_mode": aggregation_mode or DEFAULT_AGGREGATION_MODE,
        "time_window_aggregation_mode": aggregation_mode or DEFAULT_AGGREGATION_MODE,
        "parameter_ranges": ranges,
        "strategy": {
            "name": strategy_mode,
            "options": strategy_options,
        },
    }


def register_run_callbacks(app):
    """Register tuning start/stop/poll callbacks."""

    @app.callback(
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-best-results", "data", allow_duplicate=True),
        Output("tune-best-results-by-case", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-interval", "n_intervals"),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-start-button", "n_clicks"),
        State({"type": "tune-case-name", "index": ALL}, "value"),
        State({"type": "tune-case-time-start", "index": ALL}, "value"),
        State({"type": "tune-case-time-end", "index": ALL}, "value"),
        State({"type": "tune-case-average-time", "index": ALL}, "value"),
        State("tune-field-selector", "value"),
        State({"type": "tune-range-param", "index": ALL}, "value"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State("tune-batch-size", "value"),
        State("tune-max-workers", "value"),
        State("tune-strategy-mode", "data"),
        State("tune-loss-mode", "data"),
        State("tune-aggregation-mode", "data"),
        State("tune-random-max-samples", "value"),
        State("tune-resolve-spacing", "value"),
        State("tune-case-data", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def start_tuning(
        _n_clicks,
        case_names,
        time_start_values,
        time_end_values,
        average_time_values,
        selected_fields,
        param_values,
        min_values,
        max_values,
        batch_size,
        max_workers,
        strategy_mode,
        loss_mode,
        aggregation_mode,
        random_max_samples,
        resolve_spacing,
        case_data,
        active_job,
    ):
        """Validate the tuning inputs, then launch the background worker."""
        if active_job:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

        validation_message = build_validation_message(
            case_names,
            time_start_values,
            time_end_values,
            average_time_values,
            selected_fields,
            param_values,
            min_values,
            max_values,
            batch_size,
            max_workers,
            strategy_mode,
            loss_mode,
            aggregation_mode,
            random_max_samples,
            resolve_spacing,
            case_data,
        )
        if validation_message:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, validation_message

        try:
            active_job = start_tuning_job(
                build_request_payload(
                    case_names,
                    time_start_values,
                    time_end_values,
                    average_time_values,
                    selected_fields,
                    param_values,
                    min_values,
                    max_values,
                    batch_size,
                    max_workers,
                    strategy_mode,
                    loss_mode,
                    aggregation_mode,
                    random_max_samples,
                    resolve_spacing,
                )
            )
        except Exception as exc:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, str(exc)

        return active_job, empty_status_payload(), [], [], {}, False, 0, ""

    @app.callback(
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-stop-button", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def stop_tuning(_n_clicks, active_job):
        """Request a graceful stop and keep polling until the tuner exits."""
        if not active_job:
            return no_update, no_update, no_update, no_update, no_update
        stop_tuning_job(active_job)
        return no_update, no_update, no_update, no_update, ""

    @app.callback(
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-best-results", "data", allow_duplicate=True),
        Output("tune-best-results-by-case", "data", allow_duplicate=True),
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-interval", "n_intervals"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def poll_tuning_status(_tick, active_job):
        """Poll the worker status file and retire finished jobs."""
        if not active_job:
            return no_update, no_update, no_update, no_update, no_update, True, no_update

        active_tuning_job(active_job).heartbeat()
        status = read_tuning_status(active_job.get("status_path"))
        top_results = status.get("top_results", [])
        retained_results = read_tuning_results(active_job.get("results_path"))
        best_results = retained_results.get("best_results", [])
        best_results_by_case = retained_results.get("best_results_by_case", {})
        state = status.get("state", "idle")
        validation_message = ""
        active_job_out = active_job
        interval_disabled = False

        if state in {"stopped", "error", "finished"} or active_job_exited(active_job):
            clear_finished_job(active_job)
            active_job_out = {}
            interval_disabled = True
            if state == "error":
                validation_message = format_tuning_error(status, active_job)

        return status, top_results, best_results, best_results_by_case, active_job_out, interval_disabled, validation_message

    @app.callback(
        Output("tune-loss-runs", "data", allow_duplicate=True),
        Output("tune-loss-run-interval", "disabled", allow_duplicate=True),
        Output("tune-loss-run-message", "children"),
        Input({"type": "tune-loss-run-button", "action": ALL}, "n_clicks_timestamp"),
        State("tune-top-results", "data"),
        State({"type": "tune-case-name", "index": ALL}, "value"),
        State({"type": "tune-case-time-start", "index": ALL}, "value"),
        State({"type": "tune-case-time-end", "index": ALL}, "value"),
        State({"type": "tune-case-average-time", "index": ALL}, "value"),
        State("tune-field-selector", "value"),
        State("tune-loss-runs", "data"),
        prevent_initial_call=True,
    )
    def run_result_loss(
        _button_timestamps,
        top_results,
        case_names,
        time_start_values,
        time_end_values,
        average_time_values,
        selected_fields,
        loss_runs,
    ):
        """Run listed results through the window loss runner or normal SCM runner."""
        triggered_id = ctx.triggered_id
        if not isinstance(triggered_id, dict):
            return no_update, no_update, no_update
        trigger_value = ctx.triggered[0].get("value") if ctx.triggered else None
        if isinstance(trigger_value, list):
            values = trigger_value
        else:
            values = [trigger_value]
        try:
            click_timestamp = max(float(value) for value in values if value is not None)
        except (TypeError, ValueError):
            click_timestamp = -1.0
        if click_timestamp <= 0:
            return no_update, no_update, no_update

        action = triggered_id.get("action")
        if action not in {"window", "complete"}:
            return no_update, no_update, no_update
        rank_key = str(action)
        loss_runs = dict(loss_runs or {})
        existing_run = loss_runs.get(rank_key) or {}
        if float(existing_run.get("click_timestamp", -1.0) or -1.0) >= click_timestamp:
            return no_update, no_update, no_update
        if existing_run.get("state") == "running":
            return no_update, no_update, no_update

        param_sets = [result.get("params", {}) for result in (top_results or [])[:16]]
        if not param_sets:
            return no_update, no_update, "No result rows are available to run."
        try:
            if action == "window":
                case_configs = _case_configs_from_rows(
                    case_names,
                    time_start_values,
                    time_end_values,
                    average_time_values,
                )
                selected_case_names = [config["case_name"] for config in case_configs]
                run_case_configs = case_configs
                run_fields = list(selected_fields or [])
                message_prefix = "Started windowed loss run"
            else:
                selected_case_names = _case_list(case_names)
                run_case_configs = None
                run_fields = []
                message_prefix = "Started complete loss run"
            run_data = start_loss_run(
                selected_case_names,
                run_fields,
                param_sets,
                rank=action,
                case_configs=run_case_configs,
                run_mode=action,
            )
        except Exception as exc:
            return no_update, no_update, str(exc)

        run_data["click_timestamp"] = click_timestamp
        loss_runs[rank_key] = run_data
        return loss_runs, False, (
            f"{message_prefix} for {len(param_sets)} result rows "
            f"(pid {run_data['pid']}). Log: {run_data['log_path']}"
        )

    @app.callback(
        Output("tune-loss-runs", "data", allow_duplicate=True),
        Output("tune-loss-run-interval", "disabled", allow_duplicate=True),
        Input("tune-loss-run-interval", "n_intervals"),
        State("tune-loss-runs", "data"),
        prevent_initial_call=True,
    )
    def poll_result_loss_runs(_tick, loss_runs):
        """Poll ad-hoc loss runs and update row button states."""
        updated_runs, any_running = poll_loss_runs(loss_runs or {})
        return updated_runs, not any_running

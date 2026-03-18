"""Callbacks for launching, polling, cancelling, and clearing SCM runs."""

import time

from dash import ALL, Input, Output, State, callback_context, no_update

from .namelist import build_override_updates, cleanup_temp_files
from .runtime import (
    append_log_tail,
    build_case_command,
    clean_cli_option,
    cleanup_log_file,
    clear_case_status,
    format_runtime,
    get_cached_status,
    get_proc,
    launch_from_queue,
    read_log_increment,
    record_case_finish,
    snapshot_active_cases,
    snapshot_status_lists,
    start_case_process,
)
from .state import DEFAULT_STATS_NAME, RUN_STREAM_LOCK


def normalize_multicol_text(value):
    """Normalize one live multicol field into a stripped string."""
    return "" if value is None else str(value).strip()


def build_multicol_spec(param_values, min_values, max_values, npoint_values):
    """Serialize valid live multicol inputs into the `-hr` CLI format."""
    specs = []
    for row in zip(
        param_values or [],
        min_values or [],
        max_values or [],
        npoint_values or [],
    ):
        param, min_value, max_value, npoints_value = row
        row_data = {
            "param": normalize_multicol_text(param),
            "min": normalize_multicol_text(min_value),
            "max": normalize_multicol_text(max_value),
            "npoints": normalize_multicol_text(npoints_value),
        }

        if not all((row_data["param"], row_data["min"], row_data["max"], row_data["npoints"])):
            continue
        try:
            float(row_data["min"].replace("D", "E").replace("d", "e"))
            float(row_data["max"].replace("D", "E").replace("d", "e"))
            npoints = int(row_data["npoints"])
        except ValueError:
            continue
        if npoints <= 0:
            continue
        specs.append(
            f"{row_data['param']}/{row_data['min']}:{row_data['max']}/{npoints}"
        )
    return ",".join(specs)


def register_run_callbacks(app):
    """Register callbacks that mutate running, queued, completed, and failed run state."""

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Output("run-case-order", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Input("run-clear", "n_clicks"),
        State("run-running-cases", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-cases", "data"),
        State("run-case-runtimes", "data"),
        State("run-log-offsets", "data"),
        State("run-case-order", "data"),
        State("run-case-commands", "data"),
        prevent_initial_call=True,
    )
    def clear_non_running_tabs(_n_clicks, running_cases, case_logs, completed_cases, failed_cases, selected_cases, case_runtimes, log_offsets, case_order, case_commands):
        """Drop finished and selected-only panels while preserving currently running consoles."""
        running = set((running_cases or {}).keys())
        if not running and not (case_logs or completed_cases or failed_cases or selected_cases):
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update

        logs_in = dict(case_logs or {})
        completed_in = list(completed_cases or [])
        failed_in = list(failed_cases or [])
        selected_in = list(selected_cases or [])
        runtimes_in = dict(case_runtimes or {})
        offsets_in = dict(log_offsets or {})
        order_in = list(case_order or [])
        commands_in = dict(case_commands or {})

        logs_out = {k: v for k, v in logs_in.items() if k in running}
        completed_out = [k for k in completed_in if k in running]
        failed_out = [k for k in failed_in if k in running]
        runtimes_out = {k: v for k, v in runtimes_in.items() if k in running}
        order_out = [k for k in order_in if k in running]
        commands_out = {k: v for k, v in commands_in.items() if k in running}

        running_logs = {v.get("log") for v in (running_cases or {}).values() if isinstance(v, dict) and v.get("log")}
        offsets_out = {k: v for k, v in offsets_in.items() if k in running or k in running_logs}

        known_cases = (
            set(logs_in.keys())
            | set(completed_in)
            | set(failed_in)
            | set(selected_in)
            | set(runtimes_in.keys())
            | set(order_in)
            | set(commands_in.keys())
        )
        for case_name in [name for name in known_cases if name not in running]:
            clear_case_status(case_name)

        return (
            logs_out if logs_out != logs_in else no_update,
            completed_out if completed_out != completed_in else no_update,
            failed_out if failed_out != failed_in else no_update,
            runtimes_out if runtimes_out != runtimes_in else no_update,
            offsets_out if offsets_out != offsets_in else no_update,
            order_out if order_out != order_in else no_update,
            commands_out if commands_out != commands_in else no_update,
        )

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-interval", "n_intervals"),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Input("run-button", "n_clicks"),
        State("run-selected-cases", "data"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-stats-file", "data"),
        State("run-opt-max-iters", "value"),
        State("run-opt-debug", "value"),
        State("run-opt-dt-main", "value"),
        State("run-opt-dt-rad", "value"),
        State("run-opt-tout", "value"),
        State("run-opt-out-dir", "value"),
        State("run-case-runtimes", "data"),
        State("run-case-commands", "data"),
        State({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-hr-min", "index": ALL}, "value"),
        State({"type": "run-hr-max", "index": ALL}, "value"),
        State({"type": "run-hr-npoints", "index": ALL}, "value"),
        State({"type": "run-flag", "name": ALL}, "value"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        State("run-defaults", "data"),
        State("run-flag-names", "data"),
        State("run-param-meta", "data"),
        prevent_initial_call=True,
    )
    def run_selected_cases(
        _selected_clicks,
        selected_cases,
        running_cases,
        queued_cases,
        case_logs,
        completed_cases,
        failed_cases,
        selected_stats,
        opt_max_iters,
        opt_debug,
        opt_dt_main,
        opt_dt_rad,
        opt_tout,
        opt_out_dir,
        case_runtimes,
        case_commands,
        multicol_param_values,
        multicol_min_values,
        multicol_max_values,
        multicol_npoint_values,
        flag_values,
        param_values,
        defaults_data,
        flag_names_data,
        param_meta,
    ):
        """Queue selected cases, start as many as allowed, and enable interval polling."""
        if callback_context.triggered_id != "run-button":
            return (no_update,) * 9
        cases_to_run = list(selected_cases or [])
        if not cases_to_run:
            return (no_update,) * 9

        running = dict(running_cases or {})
        queued = list(queued_cases or [])
        queued_names = {item.get("case") for item in queued}
        logs = dict(case_logs or {})
        completed = list(completed_cases or [])
        failed = list(failed_cases or [])
        runtimes = dict(case_runtimes or {})
        commands = dict(case_commands or {})
        stats_name = selected_stats or DEFAULT_STATS_NAME

        overrides = build_override_updates(flag_values, param_values, defaults_data, flag_names_data, param_meta)
        cli_options = {}
        hr_cleaned = clean_cli_option(
            build_multicol_spec(
                multicol_param_values,
                multicol_min_values,
                multicol_max_values,
                multicol_npoint_values,
            )
        )
        if hr_cleaned:
            cli_options["hr"] = hr_cleaned
        for key, raw_value in (("max_iters", opt_max_iters), ("debug", opt_debug), ("dt_main", opt_dt_main), ("dt_rad", opt_dt_rad), ("tout", opt_tout)):
            cleaned = clean_cli_option(raw_value)
            if cleaned:
                cli_options[key] = cleaned
        out_dir_cleaned = clean_cli_option(opt_out_dir)
        if out_dir_cleaned and out_dir_cleaned != "output":
            cli_options["out_dir"] = out_dir_cleaned

        for case_name in cases_to_run:
            if case_name in running or case_name in queued_names:
                continue
            if case_name in snapshot_active_cases():
                continue
            clear_case_status(case_name)
            if case_name in completed:
                completed.remove(case_name)
            if case_name in failed:
                failed.remove(case_name)
            runtimes.pop(case_name, None)
            overrides_copy = {key: dict(value) for key, value in overrides.items()}
            cli_options_copy = dict(cli_options)
            queued.append({"case": case_name, "stats": stats_name, "overrides": overrides_copy, "cli_options": cli_options_copy})
            queued_names.add(case_name)
            logs[case_name] = f"--- Queued {case_name} ({stats_name}) ---\n"
            commands[case_name] = build_case_command(case_name, stats_name, cli_options_copy)

        queued, started_any = launch_from_queue(running, queued, logs)
        interval_disabled = not bool(running or queued)
        n_intervals = 0 if (started_any or queued) else no_update
        return logs, running, queued, completed, failed, interval_disabled, n_intervals, runtimes, commands

    @app.callback(
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Input("run-cancel", "n_clicks"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-case-logs", "data"),
        State("run-case-runtimes", "data"),
        prevent_initial_call=True,
    )
    def cancel_runs(_n_clicks, running_cases, queued_cases, failed_cases, case_logs, case_runtimes):
        """Terminate active runs, mark queued runs cancelled, and stop polling."""
        with RUN_STREAM_LOCK:
            running = dict(running_cases or {})
            for case_name, proc_data in snapshot_active_cases().items():
                if case_name not in running:
                    running[case_name] = proc_data
            queued = list(queued_cases or [])
            if not running and not queued:
                return no_update, no_update, no_update, no_update, True, no_update
            failed = list(failed_cases or [])
            logs = dict(case_logs or {})
            runtimes = dict(case_runtimes or {})
            now = time.time()

            for case_name, proc_data in running.items():
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                cleanup_temp_files(proc_data.get("temp_files"))
                proc = get_proc(pid)
                if proc is not None:
                    proc.terminate()
                    if proc.poll() is None:
                        proc.kill()
                record_case_finish(case_name, pid, 1)
                cleanup_log_file(log_path)
                if case_name not in failed:
                    failed.append(case_name)
                runtime_secs = now - float(proc_data.get("start_time", now))
                runtimes[case_name] = runtime_secs
                runtime_txt = format_runtime(runtime_secs)
                existing = logs.get(case_name, "")
                if existing and not existing.endswith("\n"):
                    existing += "\n"
                logs[case_name] = append_log_tail(existing, f"--- Cancelled (runtime: {runtime_txt}) ---\n")

            for item in queued:
                case_name = item.get("case")
                if not case_name:
                    continue
                record_case_finish(case_name, None, 1)
                if case_name not in failed:
                    failed.append(case_name)
                existing = logs.get(case_name, "")
                if existing and not existing.endswith("\n"):
                    existing += "\n"
                logs[case_name] = append_log_tail(existing, "--- Cancelled (never started) ---\n")

            _, failed_global = snapshot_status_lists()
            return {}, [], failed_global, logs, True, runtimes

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-selected-cases", "data", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Input("run-interval", "n_intervals"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-cases", "data"),
        State("run-case-runtimes", "data"),
        State("run-log-offsets", "data"),
        prevent_initial_call=True,
    )
    def stream_output(_tick, running_cases, queued_cases, case_logs, completed_cases, failed_cases, selected_cases, case_runtimes, log_offsets):
        """Poll running processes, stream log output, and advance the queued run list."""
        if not RUN_STREAM_LOCK.acquire(blocking=False):
            return (no_update,) * 9

        try:
            running = dict(running_cases or {})
            queued = list(queued_cases or [])
            if not running and not queued:
                return no_update, True, no_update, no_update, no_update, no_update, no_update, no_update, no_update

            active_cases = snapshot_active_cases()
            active_names = set(active_cases.keys())
            for case_name, proc_data in active_cases.items():
                if case_name not in running:
                    running[case_name] = proc_data
            if queued:
                queued = [
                    item
                    for item in queued
                    if item.get("case")
                    and item.get("case") not in active_names
                    and get_cached_status(item.get("case")) is None
                ]

            running_in = dict(running_cases or {})
            queued_in = list(queued_cases or [])
            completed_in = list(completed_cases or [])
            failed_in = list(failed_cases or [])
            selected_in = list(selected_cases or [])
            runtimes_in = dict(case_runtimes or {})
            offsets_in = dict(log_offsets or {})

            logs = dict(case_logs or {})
            completed = list(completed_in)
            failed = list(failed_in)
            selected = list(selected_in)
            runtimes = dict(runtimes_in)
            offsets = dict(offsets_in)
            finished = []
            now = time.time()
            logs_changed = False

            for case_name, proc_data in running.items():
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                offset_key = log_path or case_name
                chunk, new_offset = read_log_increment(log_path, int(offsets.get(offset_key, 0)))
                offsets[offset_key] = new_offset
                if chunk:
                    logs[case_name] = append_log_tail(logs.get(case_name, ""), chunk)
                    logs_changed = True
                proc = get_proc(pid)
                if proc is None:
                    cached_status = get_cached_status(case_name)
                    if cached_status is None:
                        misses = int(proc_data.get("missing_polls", 0)) + 1
                        proc_data["missing_polls"] = misses
                        running[case_name] = proc_data
                        if misses >= 5:
                            runtime_secs = now - float(proc_data.get("start_time", now))
                            finished.append((case_name, 1, runtime_secs))
                        continue
                    runtime_secs = now - float(proc_data.get("start_time", now))
                    finished.append((case_name, cached_status, runtime_secs))
                    continue
                status = proc.poll()
                if status is not None:
                    runtime_secs = now - float(proc_data.get("start_time", now))
                    finished.append((case_name, status, runtime_secs))

            for case_name, status, runtime_secs in finished:
                proc_data = running.pop(case_name, None) or {}
                runtimes[case_name] = runtime_secs
                runtime_txt = format_runtime(runtime_secs)
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                cleanup_temp_files(proc_data.get("temp_files"))
                offset_key = log_path or case_name
                chunk, _ = read_log_increment(log_path, int(offsets.get(offset_key, 0)))
                if chunk:
                    logs[case_name] = append_log_tail(logs.get(case_name, ""), chunk)
                    logs_changed = True
                offsets.pop(offset_key, None)
                already_finalized = record_case_finish(case_name, pid, status)
                cleanup_log_file(log_path)
                if not already_finalized:
                    if status == 0:
                        if case_name not in completed:
                            completed.append(case_name)
                        if case_name in failed:
                            failed.remove(case_name)
                    else:
                        if case_name not in failed:
                            failed.append(case_name)
                        if case_name in completed:
                            completed.remove(case_name)
                    existing = logs.get(case_name, "")
                    if existing and not existing.endswith("\n"):
                        existing += "\n"
                    result_txt = "completed" if status == 0 else f"failed (exit {status})"
                    logs[case_name] = append_log_tail(existing, f"--- {result_txt}; runtime: {runtime_txt} ---\n")
                    logs_changed = True

            running_before = dict(running)
            queued_before = list(queued)
            queued, _ = launch_from_queue(running, queued, logs)
            if running != running_before or queued != queued_before:
                logs_changed = True

            interval_disabled = not bool(running or queued)
            completed_global, failed_global = snapshot_status_lists()
            return (
                logs if logs_changed else no_update,
                interval_disabled,
                running if running != running_in else no_update,
                queued if queued != queued_in else no_update,
                completed_global if completed_global != completed_in else no_update,
                failed_global if failed_global != failed_in else no_update,
                selected if selected != selected_in else no_update,
                runtimes if runtimes != runtimes_in else no_update,
                offsets if offsets != offsets_in else no_update,
            )
        finally:
            RUN_STREAM_LOCK.release()

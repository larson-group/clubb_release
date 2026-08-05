"""Callbacks for launching, polling, cancelling, and clearing SCM runs."""

import hashlib
import json
import secrets
import time

from dash import ALL, Input, Output, State, callback_context, no_update

from .namelist import cleanup_temp_files
from .runtime import (
    append_log_tail,
    build_case_command,
    clean_cli_option,
    cleanup_log_file,
    clear_case_status,
    format_runtime,
    get_cached_status,
    normalize_task_limit,
    get_proc,
    pid_is_alive,
    read_log_increment,
    record_case_finish,
    refresh_running_runtimes,
    snapshot_active_cases,
    snapshot_status_lists,
    split_extra_cli_args,
)
from .state import DEFAULT_STATS_NAME, RUN_STREAM_LOCK
from dash_app.shared.activity import broker_jobs, set_broker_run
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name
from utilities.clubb_settings_validation import (
    apply_linked_parameter_values,
    evaluate_settings,
    values_by_name,
    values_by_setting_key,
)


def normalize_multicol_text(value):
    """Normalize one live multicol field into a stripped string."""
    return "" if value is None else str(value).strip()


def fresh_batch_request_id(request_material):
    """Return a fresh broker idempotency key for one Run Selected invocation."""
    request_hash = hashlib.sha256(str(request_material).encode()).hexdigest()[:24]
    return f"dash-run-batch-{secrets.token_urlsafe(12)}-{request_hash}"


def build_multicol_spec(
    param_values,
    min_values,
    max_values,
    npoint_values,
    available_names=None,
    linked_groups=None,
):
    """Serialize valid live multicol inputs into the `-multicol` hypergrid format.

    A persisted Run workspace can contain a documented predecessor name after
    a parameter rename. Resolve that small compatibility set before writing
    the command; reject any other complete unknown row here, before an SCM
    child and its temporary namelists are created.
    """
    specs = []
    available = set(available_names or [])
    linked_by_member = {
        member: tuple(group)
        for group in (linked_groups or [])
        for member in group
    }
    claimed_targets = set()
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
        raw_names = [name.strip() for name in row_data["param"].split("=") if name.strip()]
        parameter_names = [canonical_tunable_parameter_name(name, available) for name in raw_names]
        if not parameter_names or (available and any(name not in available for name in parameter_names)):
            raise ValueError(
                f"multicol parameter '{row_data['param']}' is not available in the selected configuration"
            )
        # A normal dropdown selection names one physical parameter.  If it is
        # part of a required-equality group, promote it to one logical axis so
        # create_multi_col_params writes identical values for all members.
        if len(parameter_names) == 1 and parameter_names[0] in linked_by_member:
            parameter_names = list(linked_by_member[parameter_names[0]])
        duplicate_targets = claimed_targets.intersection(parameter_names)
        if duplicate_targets:
            names = ", ".join(sorted(duplicate_targets))
            raise ValueError(f"linked multicol parameter already selected: {names}")
        claimed_targets.update(parameter_names)
        try:
            float(row_data["min"].replace("D", "E").replace("d", "e"))
            float(row_data["max"].replace("D", "E").replace("d", "e"))
            npoints = int(row_data["npoints"])
        except ValueError:
            continue
        if npoints <= 0:
            continue
        specs.append(
            f"{'='.join(parameter_names)}/{row_data['min']}:{row_data['max']}/{npoints}"
        )
    return ",".join(specs)


def expand_linked_parameter_values(parameter_ids, parameter_values, linked_ids, linked_values):
    """Return physical values after expanding shared linked controls.

    The Run tab only supplies component identifiers here.  The equality
    interpretation itself lives in the shared settings validator.
    """
    by_group = {
        str((component_id or {}).get("group") or ""): value
        for component_id, value in zip(linked_ids or [], linked_values or [])
    }
    by_key = values_by_setting_key(parameter_ids, parameter_values)
    expanded_by_key = apply_linked_parameter_values(by_key, by_group)
    return [
        expanded_by_key.get(
            f"{str((component_id or {}).get('file') or '')}:{str((component_id or {}).get('name') or '')}",
            value,
        )
        for component_id, value in zip(parameter_ids or [], parameter_values or [])
    ]


def launch_broker_queue(running, queued, logs, max_run_procs):
    """Start queued native-UI cases through the same broker lifecycle as MCP.

    The queue remains a Run-tab presentation concern, but subprocess ownership
    and global admission live in the broker.  That keeps reload/cancel/retry
    semantics identical for dashboard and agent initiated work.
    """
    from dash_app.shared.broker_client import perform_action

    queue = list(queued or [])
    launched = False
    failures: list[str] = []
    limit = normalize_task_limit(max_run_procs)
    while queue and len(running) < limit:
        item = dict(queue.pop(0) or {})
        case_name = str(item.get("case") or "")
        if not case_name or case_name in running:
            continue
        try:
            result = perform_action(
                "launch_dashboard_scm_request",
                {
                    "case": case_name,
                    "stats_file": item.get("stats") or DEFAULT_STATS_NAME,
                    "config": item.get("config") or "default",
                    "overrides": item.get("overrides") or {},
                    "cli_options": item.get("cli_options") or {},
                },
                internal=True,
            )
            proc_data = dict(result.get("proc_data") or {})
            if not proc_data:
                raise RuntimeError("dashboard broker did not return SCM process metadata")
            proc_data["broker_managed"] = True
            running[case_name] = proc_data
            logs[case_name] = f"--- Running {case_name} ({result.get('stats_file') or DEFAULT_STATS_NAME}, config {result.get('config') or 'default'}) ---\n"
            launched = True
        except (OSError, RuntimeError, ValueError) as exc:
            logs[case_name] = append_log_tail(logs.get(case_name, ""), f"--- Could not start {case_name}: {exc}\n")
            failures.append(case_name)
    return queue, launched, failures


def launch_broker_batch(running, queued, logs, max_workers, request_id):
    """Submit one native multi-case selection through the shared batch service."""
    from dash_app.shared.broker_client import perform_action

    queue = list(queued or [])
    if not queue:
        return queue, False, []
    first = dict(queue[0] or {})
    cases = [str(item.get("case") or "") for item in queue if str(item.get("case") or "")]
    if not cases:
        return [], False, []
    cli_options = dict(first.get("cli_options") or {})
    common_key = json.dumps(
        {
            "stats": first.get("stats") or DEFAULT_STATS_NAME,
            "config": first.get("config") or "default",
            "overrides": first.get("overrides") or {},
            "cli_options": cli_options,
        },
        sort_keys=True,
        default=str,
    )
    if any(
        json.dumps(
            {
                "stats": item.get("stats") or DEFAULT_STATS_NAME,
                "config": item.get("config") or "default",
                "overrides": item.get("overrides") or {},
                "cli_options": item.get("cli_options") or {},
            },
            sort_keys=True,
            default=str,
        ) != common_key
        for item in queue
    ):
        # A queue containing different settings is represented as separate
        # native selections; do not silently merge incompatible requests.
        return queue, False, cases

    typed_options = {}
    for key in ("max_iters", "dt_main", "dt_rad", "tout"):
        value = cli_options.get(key)
        if value in (None, ""):
            continue
        typed_options[key] = int(value) if key == "max_iters" else float(value)
    batch_request_id = str(request_id or "").strip()
    try:
        result = perform_action(
            "domain_submit_scm_batch",
            {
                "request": {
                    "request_id": batch_request_id,
                    "cases": cases,
                    "stats_file": first.get("stats") or DEFAULT_STATS_NAME,
                    "config": first.get("config") or "default",
                    "overrides": {},
                    "run_options": typed_options,
                    "max_workers": normalize_task_limit(max_workers),
                },
                "native_overrides": first.get("overrides") or {},
                "native_cli_options": cli_options,
            },
            internal=True,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        message = f"--- Could not submit SCM batch: {exc}\n"
        for case_name in cases:
            logs[case_name] = append_log_tail(logs.get(case_name, ""), message)
        return [], False, cases
    if result.get("error"):
        message = f"--- Could not submit SCM batch: {result['error'].get('message') or 'broker rejected request'}\n"
        for case_name in cases:
            logs[case_name] = append_log_tail(logs.get(case_name, ""), message)
        return [], False, cases

    for child in result.get("children") or []:
        case_name = str(child.get("case") or "")
        if not case_name:
            continue
        state = str(child.get("state") or "queued")
        logs[case_name] = f"--- Batch {result.get('batch_id')}: {state} {case_name} ---\n"
        runtime = dict(child.get("runtime") or {})
        proc_data = dict(runtime.get("proc_data") or {})
        if not proc_data and runtime.get("pid"):
            proc_data = runtime
        if proc_data and state in {"starting", "running"}:
            proc_data["broker_managed"] = True
            running[case_name] = proc_data
    # The broker owns any not-yet-started children; keeping them in the Dash
    # queue would make the legacy callback submit them a second time.
    return [], bool(running), []


def discard_terminal_broker_runs(records, remove_record):
    """Forget completed broker SCM records after an explicit Run-tab clear.

    The broker deliberately retains terminal records so a Dash reload can
    reconnect a just-finished job.  Once the user presses Clear, that durable
    history must no longer rehydrate a console.  Running/stopping records are
    left alone because Clear promises to preserve live work.
    """
    removed = []
    for raw_case, raw_record in (records or {}).items():
        record = dict(raw_record or {})
        state = str(record.get("state") or "").lower()
        if state in {"running", "stopping"}:
            continue
        case_name = str(record.get("case") or raw_case or "").strip()
        if not case_name:
            continue
        remove_record(case_name, None)
        removed.append(case_name)
    return removed


def register_run_callbacks(app):
    @app.callback(
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-interval", "n_intervals", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Input("dashboard-request", "data"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-case-logs", "data"),
        State("run-case-runtimes", "data"),
        State("run-case-commands", "data"),
        State("run-log-offsets", "data"),
        prevent_initial_call=True,
    )
    def adopt_agent_run(
        request,
        running_cases,
        queued_cases,
        completed_cases,
        failed_cases,
        case_logs,
        case_runtimes,
        case_commands,
        log_offsets,
    ):
        """Display an agent-started normal SCM process in the existing Run tab."""
        if not request or request.get("tab") != "run" or request.get("operation") != "start":
            return (no_update,) * 10
        case_name = str(request.get("case") or "")
        if not case_name:
            return (no_update,) * 10
        running = dict(running_cases or {})
        queued = list(queued_cases or [])
        completed = list(completed_cases or [])
        failed = list(failed_cases or [])
        logs = dict(case_logs or {})
        runtimes = dict(case_runtimes or {})
        commands = dict(case_commands or {})
        offsets = dict(log_offsets or {})
        proc_data = dict(request.get("proc_data") or {})
        proc_data["broker_managed"] = True
        header = f"--- Agent-started {case_name} ({request.get('stats_file') or DEFAULT_STATS_NAME}, config {request.get('config') or 'default'}) ---\n"
        logs.setdefault(case_name, header)
        commands.setdefault(
            case_name,
            build_case_command(
                case_name,
                request.get("stats_file") or DEFAULT_STATS_NAME,
                request.get("cli_options") or {},
                request.get("config") or "default",
            ),
        )

        # SCM children now belong to the durable broker, not this Dash
        # process. ``snapshot_active_cases`` is process-local and therefore
        # normally empty here for a healthy broker-owned run. Adopt the
        # broker-supplied PID directly so the console immediately becomes
        # active and its interval can derive elapsed time and ETA.
        if pid_is_alive(proc_data.get("pid")):
            running[case_name] = proc_data
            if case_name in completed:
                completed.remove(case_name)
            if case_name in failed:
                failed.remove(case_name)
            if proc_data.get("start_time") is not None:
                try:
                    runtimes[case_name] = max(0.0, time.time() - float(proc_data["start_time"]))
                except (TypeError, ValueError):
                    pass
            return running, queued, completed, failed, logs, False, 0, runtimes, commands, offsets

        active = snapshot_active_cases().get(case_name)
        if active:
            running[case_name] = active
            if case_name in completed:
                completed.remove(case_name)
            if case_name in failed:
                failed.remove(case_name)
            return running, queued, completed, failed, logs, False, 0, runtimes, commands, offsets

        # A short SCM case can finish between the server action and the
        # browser handoff.  Preserve that terminal run instead of silently
        # dropping it just because it is no longer in RUN_ACTIVE_CASES.
        returncode = get_cached_status(case_name)
        if returncode is None:
            return (no_update,) * 10
        log_path = proc_data.get("log")
        chunk, offset = read_log_increment(log_path, 0)
        if chunk:
            logs[case_name] = append_log_tail(logs.get(case_name, header), chunk)
        if log_path:
            offsets[log_path] = offset
        if returncode == 0:
            if case_name not in completed:
                completed.append(case_name)
            if case_name in failed:
                failed.remove(case_name)
        else:
            if case_name not in failed:
                failed.append(case_name)
            if case_name in completed:
                completed.remove(case_name)
        if proc_data.get("start_time") is not None:
            runtimes[case_name] = max(0.0, time.time() - float(proc_data["start_time"]))
        return running, queued, completed, failed, logs, not bool(running or queued), 0, runtimes, commands, offsets

    @app.callback(
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-interval", "n_intervals", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Input("dashboard-broker-jobs", "data"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-case-logs", "data"),
        State("run-case-runtimes", "data"),
        State("run-case-commands", "data"),
        State("run-log-offsets", "data"),
        prevent_initial_call=True,
    )
    def restore_broker_runs(
        jobs,
        running_cases,
        queued_cases,
        completed_cases,
        failed_cases,
        case_logs,
        case_runtimes,
        case_commands,
        log_offsets,
    ):
        """Rehydrate Run consoles from broker-owned SCM jobs after Dash reloads."""
        records = (jobs or {}).get("runs") or {}
        if not isinstance(records, dict):
            return (no_update,) * 10
        running = dict(running_cases or {})
        queued = list(queued_cases or [])
        completed = list(completed_cases or [])
        failed = list(failed_cases or [])
        logs = dict(case_logs or {})
        runtimes = dict(case_runtimes or {})
        commands = dict(case_commands or {})
        offsets = dict(log_offsets or {})
        saw_record = False

        for raw_case, raw_record in records.items():
            record = dict(raw_record or {})
            proc_data = record.get("proc_data")
            if not isinstance(proc_data, dict):
                continue
            saw_record = True
            case_name = str(record.get("case") or raw_case)
            state = str(record.get("state") or "")
            proc = dict(proc_data)
            proc["broker_managed"] = True
            header = f"--- Agent-broker {case_name} ({record.get('stats_file') or DEFAULT_STATS_NAME}, config {record.get('config') or 'default'}) ---\n"
            tail = str(record.get("log_tail") or "")
            logs[case_name] = header + tail
            log_path = proc.get("log")
            if log_path:
                offsets[log_path] = int(record.get("log_offset") or len(tail))
            commands.setdefault(
                case_name,
                build_case_command(
                    case_name,
                    record.get("stats_file") or DEFAULT_STATS_NAME,
                    record.get("cli_options") or {},
                    record.get("config") or "default",
                ),
            )
            if proc.get("start_time") is not None:
                runtimes[case_name] = max(0.0, time.time() - float(proc["start_time"]))
            if state in {"running", "stopping"}:
                running[case_name] = proc
                if case_name in completed:
                    completed.remove(case_name)
                if case_name in failed:
                    failed.remove(case_name)
                continue
            running.pop(case_name, None)
            if state == "finished" or record.get("returncode") == 0:
                if case_name not in completed:
                    completed.append(case_name)
                if case_name in failed:
                    failed.remove(case_name)
            elif state:
                if case_name not in failed:
                    failed.append(case_name)
                if case_name in completed:
                    completed.remove(case_name)

        if not saw_record:
            return (no_update,) * 10
        return running, queued, completed, failed, logs, not bool(running or queued), 0, runtimes, commands, offsets

    """Register callbacks that mutate running, queued, completed, and failed run state."""

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Output("run-case-order", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Output("run-selected-cases", "data", allow_duplicate=True),
        Output("run-open-cases", "data", allow_duplicate=True),
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
        removed_broker_runs = discard_terminal_broker_runs(
            (broker_jobs().get("runs") or {}),
            set_broker_run,
        )
        if not running and not (case_logs or completed_cases or failed_cases or selected_cases or removed_broker_runs):
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, [], []

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
            [],
            [],
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
        Output("run-max-tasks-active", "data", allow_duplicate=True),
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
        State("run-opt-extra-args", "value"),
        State("run-max-tasks", "value"),
        State("run-batch-size", "value"),
        State("run-case-runtimes", "data"),
        State("run-case-commands", "data"),
        State({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-hr-min", "index": ALL}, "value"),
        State({"type": "run-hr-max", "index": ALL}, "value"),
        State({"type": "run-hr-npoints", "index": ALL}, "value"),
        State({"type": "run-flag", "name": ALL}, "value"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "id"),
        State({"type": "run-linked-param", "group": ALL}, "id"),
        State({"type": "run-linked-param", "group": ALL}, "value"),
        State("run-settings-schema", "data"),
        State({"type": "run-flag", "name": ALL}, "id"),
        State("run-selected-config", "data"),
        State("run-tunable-names", "data"),
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
        opt_extra_args,
        max_tasks_value,
        batch_size_value,
        case_runtimes,
        case_commands,
        multicol_param_values,
        multicol_min_values,
        multicol_max_values,
        multicol_npoint_values,
        flag_values,
        param_values,
        param_ids,
        linked_ids,
        linked_values,
        settings_schema,
        flag_ids,
        selected_config,
        tunable_names,
    ):
        """Queue selected cases, start as many as allowed, and enable interval polling."""
        if callback_context.triggered_id != "run-button":
            return (no_update,) * 10
        cases_to_run = list(selected_cases or [])
        if not cases_to_run:
            return (no_update,) * 10

        running = dict(running_cases or {})
        queued = list(queued_cases or [])
        queued_names = {item.get("case") for item in queued}
        logs = dict(case_logs or {})
        completed = list(completed_cases or [])
        failed = list(failed_cases or [])
        runtimes = dict(case_runtimes or {})
        commands = dict(case_commands or {})
        stats_name = selected_stats or DEFAULT_STATS_NAME
        config_name = clean_cli_option(selected_config) or "default"
        max_tasks = normalize_task_limit(max_tasks_value)

        param_values = expand_linked_parameter_values(param_ids, param_values, linked_ids, linked_values)
        evaluation = evaluate_settings(
            settings_schema or {},
            flag_values={name: bool(value) for name, value in values_by_name(flag_ids, flag_values).items()},
            parameter_values=values_by_setting_key(param_ids, param_values),
        )
        overrides = {
            file_name: dict(values)
            for file_name, values in dict(evaluation.get("overrides") or {}).items()
        }
        resolution_issues = list(evaluation.get("issues") or [])
        resolver_errors = [issue for issue in resolution_issues if issue.get("severity") == "error"]
        if resolver_errors:
            message = "; ".join(str(issue.get("message") or "Invalid CLUBB settings.") for issue in resolver_errors)
            for case_name in cases_to_run:
                logs[case_name] = f"--- Settings resolver blocked launch ---\n{message}\n"
                if case_name not in failed:
                    failed.append(case_name)
            return logs, running, queued, completed, failed, not bool(running or queued), no_update, runtimes, commands, max_tasks
        resolver_notes = [
            str(issue.get("message") or "")
            for issue in resolution_issues
            if issue.get("severity") == "warning"
        ]
        cli_options = {}
        try:
            multicol_cleaned = clean_cli_option(
                build_multicol_spec(
                    multicol_param_values,
                    multicol_min_values,
                    multicol_max_values,
                    multicol_npoint_values,
                    tunable_names,
                    evaluation.get("linked_parameter_groups") or [],
                )
            )
        except ValueError as exc:
            for case_name in cases_to_run:
                logs[case_name] = f"--- Multicol setup error ---\n{exc}\n"
                if case_name not in failed:
                    failed.append(case_name)
                if case_name in completed:
                    completed.remove(case_name)
            return (
                logs,
                running,
                queued,
                completed,
                failed,
                not bool(running or queued),
                no_update,
                runtimes,
                commands,
                max_tasks,
            )
        if multicol_cleaned:
            cli_options["multicol"] = multicol_cleaned
            batch_size_cleaned = clean_cli_option(batch_size_value)
            if batch_size_cleaned:
                try:
                    batch_size_float = float(batch_size_cleaned)
                except (TypeError, ValueError):
                    batch_size_float = None
                if (
                    batch_size_float is not None
                    and int(batch_size_float) == batch_size_float
                    and batch_size_float >= 1
                ):
                    batch_size_int = int(batch_size_float)
                    cli_options["batch_size"] = str(batch_size_int)
        for key, raw_value in (("max_iters", opt_max_iters), ("debug", opt_debug), ("dt_main", opt_dt_main), ("dt_rad", opt_dt_rad), ("tout", opt_tout)):
            cleaned = clean_cli_option(raw_value)
            if cleaned:
                cli_options[key] = cleaned
        out_dir_cleaned = clean_cli_option(opt_out_dir)
        if out_dir_cleaned and out_dir_cleaned != "output":
            cli_options["out_dir"] = out_dir_cleaned
        try:
            extra_args = split_extra_cli_args(opt_extra_args)
        except ValueError as exc:
            for case_name in cases_to_run:
                logs[case_name] = f"--- Optional run args parse error ---\n{exc}\n"
                if case_name not in failed:
                    failed.append(case_name)
                if case_name in completed:
                    completed.remove(case_name)
            return (
                logs,
                running,
                queued,
                completed,
                failed,
                not bool(running or queued),
                no_update,
                runtimes,
                commands,
                max_tasks,
            )
        if extra_args:
            cli_options["extra_args"] = extra_args

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
            queued.append({"case": case_name, "stats": stats_name, "config": config_name, "overrides": overrides_copy, "cli_options": cli_options_copy})
            queued_names.add(case_name)
            note_text = ("\nSettings resolver: " + "; ".join(resolver_notes)) if resolver_notes else ""
            logs[case_name] = f"--- Queued {case_name} ({stats_name}, config {config_name}) ---{note_text}\n"
            commands[case_name] = build_case_command(case_name, stats_name, cli_options_copy, config_name)

        request_material = json.dumps(
            {
                "cases": cases_to_run,
                "stats": stats_name,
                "config": config_name,
                "overrides": overrides,
                "cli_options": cli_options,
            },
            sort_keys=True,
            default=str,
        )
        batch_request_id = fresh_batch_request_id(request_material)
        queued, started_any, launch_failures = launch_broker_batch(
            running,
            queued,
            logs,
            max_tasks,
            batch_request_id,
        )
        for case_name in launch_failures:
            if case_name not in failed:
                failed.append(case_name)
            if case_name in completed:
                completed.remove(case_name)
        interval_disabled = not bool(running or queued)
        n_intervals = 0 if (started_any or queued) else no_update
        return logs, running, queued, completed, failed, interval_disabled, n_intervals, runtimes, commands, max_tasks

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
                if proc_data.get("broker_managed"):
                    # The broker process owns the Popen handle.  Ask it to
                    # stop the process, but retain this console until the
                    # broker snapshot reports the terminal state.
                    try:
                        from dash_app.shared.broker_client import perform_action

                        perform_action("stop_run", {"case": case_name}, internal=True)
                        existing = logs.get(case_name, "")
                        if existing and not existing.endswith("\n"):
                            existing += "\n"
                        logs[case_name] = append_log_tail(existing, "--- Stop requested from dashboard; waiting for broker ---\n")
                    except (OSError, RuntimeError, ValueError):
                        pass
                    continue
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
            broker_waiting = any(bool(proc_data.get("broker_managed")) for proc_data in running.values())
            return (running if broker_waiting else {}), [], failed_global, logs, not broker_waiting, runtimes

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
        State("run-max-tasks-active", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-cases", "data"),
        State("run-case-runtimes", "data"),
        State("run-log-offsets", "data"),
        prevent_initial_call=True,
    )
    def stream_output(_tick, running_cases, queued_cases, max_tasks_active, case_logs, completed_cases, failed_cases, selected_cases, case_runtimes, log_offsets):
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
            now = time.time()
            runtimes = refresh_running_runtimes(runtimes_in, running, now)
            offsets = dict(offsets_in)
            finished = []
            logs_changed = False
            broker_records = dict((broker_jobs().get("runs") or {}))

            for case_name, proc_data in running.items():
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                offset_key = log_path or case_name
                chunk, new_offset = read_log_increment(log_path, int(offsets.get(offset_key, 0)))
                offsets[offset_key] = new_offset
                if chunk:
                    logs[case_name] = append_log_tail(logs.get(case_name, ""), chunk)
                    logs_changed = True
                if proc_data.get("broker_managed"):
                    broker_record = dict(broker_records.get(case_name) or {})
                    broker_state = str(broker_record.get("state") or "")
                    if broker_state in {"running", "stopping"}:
                        continue
                    if broker_state:
                        status = broker_record.get("returncode")
                        if status is None:
                            status = 0 if broker_state == "finished" else 1
                        runtime_secs = now - float(proc_data.get("start_time", now))
                        finished.append((case_name, int(status), runtime_secs))
                        continue
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

            running_before = dict(running)
            queued_before = list(queued)
            queued, _, launch_failures = launch_broker_queue(running, queued, logs, max_tasks_active)
            for case_name in launch_failures:
                if case_name not in failed:
                    failed.append(case_name)
                if case_name in completed:
                    completed.remove(case_name)
            if running != running_before or queued != queued_before:
                logs_changed = True

            interval_disabled = not bool(running or queued)
            completed_global, failed_global = snapshot_status_lists()
            # Broker-owned runs finalize in the broker process, so retain the
            # UI's local terminal list as well as any legacy in-process state.
            completed_global = sorted(set(completed_global) | set(completed))
            failed_global = sorted(set(failed_global) | set(failed))
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

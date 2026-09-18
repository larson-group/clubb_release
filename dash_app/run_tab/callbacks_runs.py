"""Compact broker-backed lifecycle for native SCM runs."""

import hashlib
import json
import secrets
import time

from dash import ALL, Input, Output, State, callback_context, html, no_update

from .runtime import (
    clean_cli_option,
    normalize_task_limit,
    output_directory_details,
    split_extra_cli_args,
)
from .state import DEFAULT_STATS_NAME
from dash_app.compile_tab.build_selector import selected_launch_target
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name
from utilities.output_paths import resolve_output_dir
from utilities.clubb_settings_validation import (
    apply_linked_parameter_values,
    evaluate_settings,
    format_setting_value,
    values_by_name,
    values_by_setting_key,
)


RUN_OVERWRITE_OPEN = "run-overwrite-modal"
RUN_OVERWRITE_CLOSED = "run-overwrite-modal run-overwrite-modal-hidden"


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
    """Serialize valid live multicol inputs into the `-multicol` format."""
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
        if not all(row_data.values()):
            continue
        raw_names = [name.strip() for name in row_data["param"].split("=") if name.strip()]
        parameter_names = [
            canonical_tunable_parameter_name(name, available) for name in raw_names
        ]
        if not parameter_names or (
            available and any(name not in available for name in parameter_names)
        ):
            raise ValueError(
                f"multicol parameter '{row_data['param']}' is not available in the selected configuration"
            )
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
        if npoints > 0:
            specs.append(
                f"{'='.join(parameter_names)}/{row_data['min']}:{row_data['max']}/{npoints}"
            )
    return ",".join(specs)


def expand_linked_parameter_values(
    parameter_ids, parameter_values, linked_ids, linked_values
):
    """Return physical values after expanding shared linked controls."""
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


def complete_run_overrides(evaluation):
    """Freeze every effective native setting instead of mutable config deltas."""
    normalized_parameters = dict(evaluation.get("normalized_parameters") or {})
    return {
        "flags": {
            str(name): format_setting_value(value)
            for name, value in dict(evaluation.get("normalized_flags") or {}).items()
        },
        **{
            file_name: {
                str(name): format_setting_value(value)
                for name, value in dict(normalized_parameters.get(file_name) or {}).items()
            }
            for file_name in ("tunable", "silhs")
        },
    }


def run_output_rename_available(proposed_output, pending):
    """Return whether a pending Run can move to a different output target."""
    proposed = clean_cli_option(proposed_output)
    current = clean_cli_option((pending or {}).get("output_dir"))
    if not proposed or not current:
        return False
    try:
        return resolve_output_dir(proposed).resolve() != resolve_output_dir(current).resolve()
    except (OSError, TypeError, ValueError):
        return False


def run_output_detail_components(details):
    """Render the requested output-folder facts for the overwrite dialog."""
    case_count = int(details.get("case_count") or 0)
    cases = f"{case_count} case" if case_count == 1 else f"{case_count} cases"
    return [
        html.Div(
            [
                html.Span(label, className="run-overwrite-fact-label"),
                html.Span(value, className="run-overwrite-fact-value"),
            ],
            className="run-overwrite-fact",
        )
        for label, value in (
            ("Folder", str(details.get("path") or "")),
            ("Created", str(details.get("created") or "Unknown")),
            ("Last edited", str(details.get("last_edited") or "Unknown")),
            ("Cases", f"{cases} with *_stats.nc output"),
        )
    ]


def prepared_run_with_output(pending, output_dir):
    """Copy a frozen Run submission while changing only its output target."""
    updated = dict(pending or {})
    updated["output_dir"] = clean_cli_option(output_dir) or "output"
    updated["cli_options"] = dict(updated.get("cli_options") or {})
    if updated["output_dir"] == "output":
        updated["cli_options"].pop("out_dir", None)
    else:
        updated["cli_options"]["out_dir"] = updated["output_dir"]
    return updated


def submit_prepared_run(pending, perform_action):
    """Submit one previously validated and frozen Run request."""
    request_material = json.dumps(pending, sort_keys=True, default=str)
    result = perform_action(
        "domain_submit_scm_batch",
        {
            "request": {
                "request_id": fresh_batch_request_id(request_material),
                "cases": list(pending.get("cases") or []),
                "implementation": pending.get("implementation") or "fortran",
                "jax_profile": pending.get("jax_profile") or "cpu",
                "jax_gpu": pending.get("jax_gpu") or "",
                "jax_xla_prealloc": pending.get("jax_xla_prealloc"),
                "stats_file": pending.get("stats") or DEFAULT_STATS_NAME,
                "config": pending.get("config") or "default",
                "overrides": dict(pending.get("typed_overrides") or {}),
                "run_options": dict(pending.get("typed_options") or {}),
                "max_workers": int(pending.get("max_workers") or 1),
            },
            "native_overrides": dict(pending.get("overrides") or {}),
            "native_cli_options": dict(pending.get("cli_options") or {}),
            "submission_origin": "dash",
        },
        internal=True,
    )
    return {
        "action": "run",
        "at": time.time(),
        "job_id": result.get("job_id"),
    }


def register_run_callbacks(app):
    """Keep user commands responsive while one reducer owns lifecycle state."""

    @app.callback(
        Output("run-action-result", "data"),
        Output("run-pending-request", "data"),
        Output("run-overwrite-modal", "className"),
        Output("run-overwrite-name", "value"),
        Output("run-overwrite-message", "children"),
        Output("run-overwrite-details", "children"),
        Output("run-opt-out-dir", "value"),
        Input("run-button", "n_clicks"),
        Input("run-cancel", "n_clicks"),
        Input("run-clear", "n_clicks"),
        Input("run-overwrite-button", "n_clicks"),
        Input("run-rename-button", "n_clicks"),
        Input("run-overwrite-cancel-button", "n_clicks"),
        State("run-selected-cases", "data"),
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
        State("compile-run-implementation", "data"),
        State("compile-run-jax-profile", "data"),
        State("compile-run-jax-gpu", "data"),
        State("compile-run-jax-xla-prealloc", "data"),
        State("run-pending-request", "data"),
        State("run-overwrite-name", "value"),
        prevent_initial_call=True,
    )
    def execute_run_action(
        _run_clicks,
        _cancel_clicks,
        _clear_clicks,
        _overwrite_clicks,
        _rename_clicks,
        _overwrite_cancel_clicks,
        selected_cases,
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
        run_implementation,
        run_jax_profile,
        run_jax_gpu,
        run_jax_xla_prealloc,
        pending_request,
        proposed_output,
    ):
        from dash_app.shared.broker_client import perform_action

        trigger = callback_context.triggered_id
        if trigger == "run-clear":
            perform_action("clear_terminal_scm_session", {}, internal=True)
            return (
                {"action": "clear", "at": time.time()},
                {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update,
            )

        if trigger == "run-cancel":
            result = perform_action("domain_cancel_all_scm", {}, internal=True)
            return (
                {"action": "cancel", "at": time.time(), "result": result},
                {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update,
            )

        pending = dict(pending_request or {})
        if trigger == "run-overwrite-button":
            if not pending:
                return (no_update,) * 7
            try:
                action = submit_prepared_run(pending, perform_action)
                return (
                    action, {}, RUN_OVERWRITE_CLOSED, "", "", [],
                    pending["output_dir"],
                )
            except (OSError, RuntimeError, ValueError) as exc:
                return (
                    {"action": "error", "at": time.time(), "message": str(exc)},
                    pending,
                    RUN_OVERWRITE_OPEN,
                    proposed_output,
                    str(exc),
                    no_update,
                    no_update,
                )

        if trigger == "run-rename-button":
            if not pending or not run_output_rename_available(proposed_output, pending):
                return (no_update,) * 7
            try:
                renamed = prepared_run_with_output(pending, proposed_output)
                details = output_directory_details(renamed["output_dir"])
                if details["nonempty"]:
                    return (
                        no_update,
                        renamed,
                        RUN_OVERWRITE_OPEN,
                        proposed_output,
                        "That folder also contains files. Choose another folder or overwrite it.",
                        run_output_detail_components(details),
                        no_update,
                    )
                action = submit_prepared_run(renamed, perform_action)
                return (
                    action, {}, RUN_OVERWRITE_CLOSED, "", "", [], proposed_output,
                )
            except (OSError, RuntimeError, TypeError, ValueError) as exc:
                return (
                    {"action": "error", "at": time.time(), "message": str(exc)},
                    pending,
                    RUN_OVERWRITE_OPEN,
                    proposed_output,
                    str(exc),
                    no_update,
                    no_update,
                )

        if trigger == "run-overwrite-cancel-button":
            return no_update, {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update

        if trigger != "run-button":
            return (no_update,) * 7

        cases_to_run = list(selected_cases or [])
        if not cases_to_run:
            return (no_update,) * 7
        stats_name = selected_stats or DEFAULT_STATS_NAME
        config_name = clean_cli_option(selected_config) or "default"
        max_tasks = normalize_task_limit(max_tasks_value)
        param_values = expand_linked_parameter_values(
            param_ids, param_values, linked_ids, linked_values
        )
        evaluation = evaluate_settings(
            settings_schema or {},
            flag_values={
                name: bool(value)
                for name, value in values_by_name(flag_ids, flag_values).items()
            },
            parameter_values=values_by_setting_key(param_ids, param_values),
        )
        errors = [
            issue
            for issue in (evaluation.get("issues") or [])
            if issue.get("severity") == "error"
        ]
        if errors:
            message = "; ".join(
                str(issue.get("message") or "Invalid CLUBB settings.")
                for issue in errors
            )
            return (
                {
                    "action": "error",
                    "at": time.time(),
                    "cases": cases_to_run,
                    "message": message,
                },
                {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update,
            )
        overrides = complete_run_overrides(evaluation)
        typed_overrides = {
            name: value
            for values in overrides.values()
            for name, value in values.items()
        }
        cli_options = {}
        try:
            launch_target = selected_launch_target(
                run_implementation, jax_profile=run_jax_profile, jax_gpu=run_jax_gpu,
                jax_xla_prealloc=run_jax_xla_prealloc,
            )
            cli_options.update(
                implementation=launch_target["implementation"],
                install_dir=launch_target["install_dir"],
            )
            if launch_target["implementation"] == "jax":
                cli_options["jax_profile"] = launch_target["jax_profile"]
                cli_options["jax_gpu"] = launch_target.get("jax_gpu") or ""
                cli_options["jax_xla_prealloc"] = launch_target.get("jax_xla_prealloc")
            multicol = clean_cli_option(
                build_multicol_spec(
                    multicol_param_values,
                    multicol_min_values,
                    multicol_max_values,
                    multicol_npoint_values,
                    tunable_names,
                    evaluation.get("linked_parameter_groups") or [],
                )
            )
            extra_args = split_extra_cli_args(opt_extra_args)
        except ValueError as exc:
            return (
                {
                    "action": "error",
                    "at": time.time(),
                    "cases": cases_to_run,
                    "message": str(exc),
                },
                {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update,
            )
        if multicol:
            cli_options["multicol"] = multicol
            try:
                batch_size = int(float(clean_cli_option(batch_size_value)))
            except (TypeError, ValueError):
                batch_size = 0
            if batch_size > 0:
                cli_options["batch_size"] = str(batch_size)
        for key, raw_value in (
            ("max_iters", opt_max_iters),
            ("debug", opt_debug),
            ("dt_main", opt_dt_main),
            ("dt_rad", opt_dt_rad),
            ("tout", opt_tout),
        ):
            value = clean_cli_option(raw_value)
            if value:
                cli_options[key] = value
        output_dir = clean_cli_option(opt_out_dir)
        if output_dir and output_dir != "output":
            cli_options["out_dir"] = output_dir
        if extra_args:
            cli_options["extra_args"] = extra_args

        typed_options = {}
        for key in ("max_iters", "dt_main", "dt_rad", "tout"):
            value = cli_options.get(key)
            if value in (None, ""):
                continue
            typed_options[key] = int(value) if key == "max_iters" else float(value)
        prepared = {
            "cases": cases_to_run,
            "stats": stats_name,
            "config": config_name,
            "overrides": overrides,
            "typed_overrides": typed_overrides,
            "cli_options": cli_options,
            "typed_options": typed_options,
            "max_workers": max_tasks,
            "output_dir": output_dir or "output",
            "implementation": launch_target["implementation"],
            "jax_profile": launch_target.get("jax_profile") or "cpu",
            "jax_gpu": launch_target.get("jax_gpu") or "",
            "jax_xla_prealloc": launch_target.get("jax_xla_prealloc"),
        }
        try:
            details = output_directory_details(prepared["output_dir"])
            if details["nonempty"]:
                return (
                    no_update,
                    prepared,
                    RUN_OVERWRITE_OPEN,
                    output_dir or "output",
                    "Running here may replace matching case output files.",
                    run_output_detail_components(details),
                    no_update,
                )
            action = submit_prepared_run(prepared, perform_action)
            return action, {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            return (
                {
                    "action": "error",
                    "at": time.time(),
                    "cases": cases_to_run,
                    "message": str(exc),
                },
                {}, RUN_OVERWRITE_CLOSED, "", "", [], no_update,
            )

    @app.callback(
        Output("run-rename-button", "disabled"),
        Output("run-rename-button", "title"),
        Input("run-overwrite-name", "value"),
        State("run-pending-request", "data"),
    )
    def update_run_rename_action(proposed_output, pending):
        available = run_output_rename_available(proposed_output, pending)
        return (
            not available,
            "Rename the output folder and start the run."
            if available
            else "Enter a different output folder to rename and run.",
        )

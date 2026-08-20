"""Compact broker-backed lifecycle for native SCM runs."""

import hashlib
import json
import secrets
import time

from dash import ALL, Input, Output, State, callback_context, no_update

from .runtime import (
    clean_cli_option,
    normalize_task_limit,
    split_extra_cli_args,
)
from .state import DEFAULT_STATS_NAME
from dash_app.compile_tab.build_selector import selected_launch_target
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name
from utilities.clubb_settings_validation import (
    apply_linked_parameter_values,
    evaluate_settings,
    format_setting_value,
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


def register_run_callbacks(app):
    """Keep user commands responsive while one reducer owns lifecycle state."""

    @app.callback(
        Output("run-action-result", "data"),
        Input("run-button", "n_clicks"),
        Input("run-cancel", "n_clicks"),
        Input("run-clear", "n_clicks"),
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
        prevent_initial_call=True,
    )
    def execute_run_action(
        _run_clicks,
        _cancel_clicks,
        _clear_clicks,
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
    ):
        from dash_app.shared.broker_client import perform_action

        trigger = callback_context.triggered_id
        if trigger == "run-clear":
            perform_action("clear_terminal_scm_session", {}, internal=True)
            return {"action": "clear", "at": time.time()}

        if trigger == "run-cancel":
            result = perform_action("domain_cancel_all_scm", {}, internal=True)
            return {"action": "cancel", "at": time.time(), "result": result}

        if trigger != "run-button":
            return no_update

        cases_to_run = list(selected_cases or [])
        if not cases_to_run:
            return no_update
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
            return {
                "action": "error",
                "at": time.time(),
                "cases": cases_to_run,
                "message": message,
            }
        overrides = complete_run_overrides(evaluation)
        typed_overrides = {
            name: value
            for values in overrides.values()
            for name, value in values.items()
        }
        cli_options = {}
        try:
            launch_target = selected_launch_target(run_implementation)
            cli_options.update(
                implementation=launch_target["implementation"],
                install_dir=launch_target["install_dir"],
            )
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
            return {
                "action": "error",
                "at": time.time(),
                "cases": cases_to_run,
                "message": str(exc),
            }
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
        typed_options = {}
        for key in ("max_iters", "dt_main", "dt_rad", "tout"):
            value = cli_options.get(key)
            if value in (None, ""):
                continue
            typed_options[key] = int(value) if key == "max_iters" else float(value)
        try:
            result = perform_action(
                "domain_submit_scm_batch",
                {
                    "request": {
                        "request_id": fresh_batch_request_id(request_material),
                        "cases": cases_to_run,
                        "stats_file": stats_name,
                        "config": config_name,
                        "overrides": typed_overrides,
                        "run_options": typed_options,
                        "max_workers": max_tasks,
                    },
                    "native_overrides": overrides,
                    "native_cli_options": cli_options,
                    "submission_origin": "dash",
                },
                internal=True,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            return {
                "action": "error",
                "at": time.time(),
                "cases": cases_to_run,
                "message": str(exc),
            }
        return {
            "action": "run",
            "at": time.time(),
            "job_id": result.get("job_id"),
        }

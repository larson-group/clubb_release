"""Callbacks for run-tab flag/parameter editing and multicol UI state."""

from dash import (
    ALL,
    Input,
    Output,
    Patch,
    State,
    ClientsideFunction,
    callback_context,
    no_update,
)

from dash_app.shared.tunable_configs import tunable_config_names
from utilities.output_paths import resolve_output_dir
from utilities.clubb_settings_validation import (
    evaluate_settings,
    is_independently_tunable,
    linked_flag_updates,
    linked_parameter_members as utility_linked_parameter_members,
    values_by_name,
    values_by_setting_key,
)

from .config_state import build_tunable_config_state
from .layout import (
    build_multicol_row,
    build_right_pane,
    build_tunable_controls,
    compute_width_hints,
)
from .namelist import normalize_numeric_display
from .state import RUN_FINALIZED, RUN_LOCK, RUN_STATUS


def blank_multicol_row(row_id):
    """Return one empty multicol hypergrid row record."""
    return {"id": row_id, "param": "", "min": "", "max": "", "npoints": "4"}


def default_multicol_row(row_id, tunable_names, default_ranges=None, selected_params=None):
    """Return a new multicol row initialized to the first useful tunable parameter."""
    selected = {str(value or "").strip() for value in (selected_params or [])}
    selected.discard("")
    available = [name for name in (tunable_names or []) if name]
    param = next((name for name in available if name not in selected), available[0] if available else "")
    row = blank_multicol_row(row_id)
    row["param"] = param
    derived = (default_ranges or {}).get(param, {})
    row["min"] = derived.get("min", "")
    row["max"] = derived.get("max", "")
    return row


def unchanged_multicol_bounds(param_values):
    """Return wildcard-shaped no-op values for the multicol min/max controls."""
    count = len(param_values or [])
    return [no_update] * count, [no_update] * count


def synchronized_flag_values(flag_ids, flag_values, trigger_id=None):
    """Return live flags with only a declared checkbox-pair update applied."""
    resolved = {
        name: bool(value)
        for name, value in values_by_name(flag_ids, flag_values).items()
    }
    if isinstance(trigger_id, dict) and trigger_id.get("type") == "run-flag":
        name = str(trigger_id.get("name") or "")
        resolved.update(linked_flag_updates(name, bool(resolved.get(name))))
    return resolved


def linked_groups_from_resolution(resolution, param_meta):
    """Extract active equal-value tunable groups from the shared resolver."""
    available = {
        str(meta.get("name"))
        for meta in (param_meta or [])
        if isinstance(meta, dict) and meta.get("file") == "tunable" and meta.get("name")
    }
    groups = []
    for members in (resolution or {}).get("linked_parameter_groups") or []:
        members = [str(member) for member in members or []]
        if len(members) >= 2 and all(member in available for member in members):
            groups.append(members)
    return groups


def _fallback_rows_by_id(rows):
    """Return saved multicol rows indexed by row id."""
    return {row.get("id"): row for row in (rows or []) if isinstance(row, dict)}


def multicol_rows_from_values(
    row_order,
    param_values,
    min_values,
    max_values,
    npoint_values,
    tunable_names,
    default_ranges=None,
    fallback_rows=None,
):
    """Preserve live multicol rows across config changes."""
    available = set(tunable_names or [])
    fallback_by_id = _fallback_rows_by_id(fallback_rows)
    row_ids = list(row_order or [])
    if not row_ids and fallback_by_id:
        row_ids = [row.get("id") for row in (fallback_rows or []) if isinstance(row, dict)]
    rows = []
    for index, row_id in enumerate(row_ids):
        fallback = fallback_by_id.get(row_id, {})
        param = fallback.get("param", "") if index >= len(param_values or []) else param_values[index]
        param = str(param).strip() if param is not None else ""
        min_value = fallback.get("min", "") if index >= len(min_values or []) or min_values[index] is None else str(min_values[index])
        max_value = fallback.get("max", "") if index >= len(max_values or []) or max_values[index] is None else str(max_values[index])
        if param and param in available:
            derived = (default_ranges or {}).get(param, {})
            if not min_value.strip() and derived.get("min"):
                min_value = derived["min"]
            if not max_value.strip() and derived.get("max"):
                max_value = derived["max"]
        rows.append(
            {
                "id": row_id,
                "param": param,
                "min": min_value,
                "max": max_value,
                "npoints": fallback.get("npoints", "4") if index >= len(npoint_values or []) or npoint_values[index] is None else str(npoint_values[index]),
            }
        )
    if rows:
        return rows
    return [blank_multicol_row(0)]


def invalidate_cached_run_results(completed_cases, failed_cases):
    """Clear cached run status and completed/failed stores after a settings edit."""
    with RUN_LOCK:
        if RUN_STATUS or RUN_FINALIZED:
            RUN_STATUS.clear()
            RUN_FINALIZED.clear()

    completed_update = [] if completed_cases else no_update
    failed_update = [] if failed_cases else no_update
    return completed_update, failed_update


def output_directory_warning(value):
    """Return a concise overwrite warning for the selected output directory.

    The SCM runner intentionally reuses an output directory. This helper is
    advisory rather than a launch guard and never creates the path.
    """
    try:
        path = resolve_output_dir(value)
    except (TypeError, ValueError):
        return "⚠ Output path is invalid.", "run-output-dir-warning run-output-dir-warning--warning"

    try:
        if path.exists() and not path.is_dir():
            return "⚠ Output path is an existing file.", "run-output-dir-warning run-output-dir-warning--warning"
        if not path.is_dir():
            return "", "run-output-dir-warning run-output-dir-warning--hidden"
        count = sum(1 for _entry in path.iterdir())
    except OSError:
        return "⚠ Could not inspect this output directory.", "run-output-dir-warning run-output-dir-warning--warning"

    if count:
        noun = "item" if count == 1 else "items"
        return (
            f"⚠ {count} existing {noun}; matching output files will be overwritten.",
            "run-output-dir-warning run-output-dir-warning--warning",
        )
    return "", "run-output-dir-warning run-output-dir-warning--hidden"


def register_settings_callbacks(app):
    """Register callbacks that track dirty state and changed-field styling."""

    @app.callback(
        Output("run-output-dir-warning", "children"),
        Output("run-output-dir-warning", "className"),
        Input("run-opt-out-dir", "value"),
    )
    def show_output_directory_warning(output_directory):
        """Make intentional output reuse visible before a Run launch."""
        return output_directory_warning(output_directory)

    @app.callback(
        Output("run-settings-resolution", "data"),
        Input({"type": "run-flag", "name": ALL}, "value"),
        Input({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        State("run-settings-schema", "data"),
        State({"type": "run-flag", "name": ALL}, "id"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "id"),
    )
    def resolve_live_settings(flag_values, param_values, schema, flag_ids, param_ids):
        """Delegate all Run settings semantics to the shared evaluator.

        An explicit checkbox-pair click has a companion update in a separate
        Dash callback.  Apply that same declared update to this in-flight
        value map so the displayed validation state never lags one callback
        behind the checkbox state.
        """
        resolved_flags = synchronized_flag_values(
            flag_ids,
            flag_values,
            callback_context.triggered_id,
        )
        return evaluate_settings(
            schema or {},
            flag_values=resolved_flags,
            parameter_values=values_by_setting_key(param_ids, param_values),
        )

    @app.callback(
        Output("run-settings-resolution-note", "children"),
        Output("run-settings-resolution-note", "className"),
        Input("run-settings-resolution", "data"),
        prevent_initial_call=True,
    )
    def show_settings_resolution_note(resolution):
        """Render settings consequences server-side for asset-reload safety."""
        issues = [
            item for item in (resolution or {}).get("issues", [])
            if item.get("severity") == "error" or str(item.get("code") or "").startswith("forced_")
        ]
        if not issues:
            return "", "run-settings-resolution-note"
        message = " ".join(str(item.get("message") or "") for item in issues)
        has_error = any(item.get("severity") == "error" for item in issues)
        suffix = "run-settings-resolution-note--error" if has_error else "run-settings-resolution-note--derived"
        return message, f"run-settings-resolution-note {suffix}"

    @app.callback(
        Output({"type": "run-flag", "name": ALL}, "value", allow_duplicate=True),
        Output({"type": "run-param", "file": ALL, "name": ALL}, "value", allow_duplicate=True),
        Input("run-settings-resolution", "data"),
        Input({"type": "run-flag", "name": ALL}, "value"),
        State({"type": "run-flag", "name": ALL}, "id"),
        State({"type": "run-flag", "name": ALL}, "value"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "id"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        prevent_initial_call=True,
    )
    def apply_forced_settings(resolution, live_flag_values, flag_ids, flag_values, param_ids, param_values):
        """Apply forced values and make opposite checkbox groups switchable."""
        trigger_id = callback_context.triggered_id
        if isinstance(trigger_id, dict) and trigger_id.get("type") == "run-flag":
            try:
                source_index = list(flag_ids or []).index(trigger_id)
            except ValueError:
                source_index = -1
            source_is_enabled = source_index >= 0 and bool((live_flag_values or [])[source_index])
            related_updates = linked_flag_updates(trigger_id.get("name"), source_is_enabled)
            if related_updates:
                return (
                    [
                        (["on"] if desired else [])
                        if (desired := related_updates.get((component_id or {}).get("name"))) is not None
                        and (live_flag_values or [])[index] != (["on"] if desired else [])
                        else no_update
                        for index, component_id in enumerate(flag_ids or [])
                    ],
                    [no_update for _ in (param_ids or [])],
                )
            return [no_update for _ in (flag_ids or [])], [no_update for _ in (param_ids or [])]

        # All remaining triggers come from the shared resolver.
        forced_parameters = dict((resolution or {}).get("forced_parameters") or {})
        updated_flags = [no_update for _ in (flag_ids or [])]
        updated_params = [no_update for _ in (param_ids or [])]
        for index, component_id in enumerate(param_ids or []):
            if (component_id or {}).get("file") != "tunable":
                continue
            name = component_id.get("name")
            if name in forced_parameters:
                desired = str(forced_parameters[name])
                if index >= len(param_values or []) or str((param_values or [])[index]) != desired:
                    updated_params[index] = desired
        return updated_flags, updated_params

    @app.callback(
        Output({"type": "run-param", "file": "tunable", "name": ALL}, "value", allow_duplicate=True),
        Input({"type": "run-linked-param", "group": ALL}, "value"),
        State({"type": "run-linked-param", "group": ALL}, "id"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "id"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "value"),
        prevent_initial_call=True,
    )
    def mirror_linked_parameter_value(linked_values, linked_ids, param_ids, param_values):
        """Expand one visible linked control into its physical namelist members."""
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "run-linked-param":
            return [no_update for _ in (param_ids or [])]
        group = str(trigger_id.get("group") or "")
        try:
            position = list(linked_ids or []).index(trigger_id)
        except ValueError:
            return [no_update for _ in (param_ids or [])]
        if position >= len(linked_values or []):
            return [no_update for _ in (param_ids or [])]
        desired = "" if (linked_values or [])[position] is None else str((linked_values or [])[position])
        members = set(utility_linked_parameter_members(group))
        return [
            desired
            if (component_id or {}).get("name") in members
            and str((param_values or [])[index]) != desired
            else no_update
            for index, component_id in enumerate(param_ids or [])
        ]

    @app.callback(
        Output("run-tunable-controls", "children"),
        Output("run-linked-parameter-groups", "data"),
        Input("run-settings-resolution", "data"),
        State("run-linked-parameter-groups", "data"),
        State("run-param-meta", "data"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "id"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "value"),
        prevent_initial_call=True,
    )
    def rebuild_linked_tunable_controls(resolution, current_groups, param_meta, param_ids, param_values):
        """Rebuild only when a flag activates or releases a link rule.

        Ordinary value edits leave the component tree alone, preserving focus
        and browser persistence.  The resolver is the only source of group
        membership, so the conditional C1/C14 rule remains consistent with
        CLUBB's validation logic.
        """
        desired_groups = linked_groups_from_resolution(resolution, param_meta)
        if desired_groups == list(current_groups or []):
            return no_update, no_update
        values_by_name = {
            str((component_id or {}).get("name") or ""): value
            for component_id, value in zip(param_ids or [], param_values or [])
        }
        forced = dict((resolution or {}).get("forced_parameters") or {})
        for group in desired_groups:
            forced_value = next((forced[name] for name in group if name in forced), None)
            if forced_value is not None:
                values_by_name[group[0]] = forced_value
        entries = [
            {"name": str(meta["name"]), "value": values_by_name.get(str(meta["name"]), "")}
            for meta in (param_meta or [])
            if isinstance(meta, dict) and meta.get("file") == "tunable" and meta.get("name")
        ]
        label_width_px, value_width_px, _ = compute_width_hints([entry["name"] for entry in entries])
        return (
            build_tunable_controls(entries, desired_groups, label_width_px, value_width_px, normalize_numeric_display),
            desired_groups,
        )

    @app.callback(
        Output("run-selected-config", "data"),
        Output("run-settings-schema", "data"),
        Output("run-param-meta", "data"),
        Output("run-tunable-names", "data"),
        Output("run-tunable-default-ranges", "data"),
        Output("run-linked-parameter-groups", "data", allow_duplicate=True),
        Output("run-settings-resolution", "data", allow_duplicate=True),
        Output("run-right-pane", "children"),
        Output("run-multicol-next-id", "data", allow_duplicate=True),
        Output("run-multicol-row-order", "data", allow_duplicate=True),
        Output("run-multicol-rows-state", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input("run-config-reset-signal", "data"),
        State("run-tunable-configs", "data"),
        State("run-selected-config", "data"),
        State("run-multicol-next-id", "data"),
        State("run-multicol-row-order", "data"),
        State("run-multicol-rows-state", "data"),
        State({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-hr-min", "index": ALL}, "value"),
        State({"type": "run-hr-max", "index": ALL}, "value"),
        State({"type": "run-hr-npoints", "index": ALL}, "value"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        prevent_initial_call=True,
    )
    def select_tunable_config(
        reset_signal,
        configs,
        current_config,
        next_id,
        row_order,
        saved_multicol_rows,
        multicol_param_values,
        multicol_min_values,
        multicol_max_values,
        multicol_npoint_values,
        completed_cases,
        failed_cases,
    ):
        """Switch the run tab to another tunable config and rebuild config-dependent controls."""
        if not isinstance(reset_signal, dict):
            return (no_update,) * 13

        config_name = str(reset_signal.get("name", "")).strip()
        if not config_name or config_name not in tunable_config_names(configs):
            return (no_update,) * 13
        config_state = build_tunable_config_state(config_name, configs)
        # A config button is an explicit reset to that config's complete
        # checked-in state.  Do not preserve a previous workspace's parameter
        # grid or live right-pane values under a newly selected config.
        multicol_rows = [
            default_multicol_row(
                0,
                config_state["tunable_names"],
                config_state["tunable_default_ranges"],
            )
        ]
        multicol_row_order = [row["id"] for row in multicol_rows]
        multicol_next_id = 1
        right_pane_data = {"tunable_configs": configs or [], "multicol_rows": multicol_rows, **config_state}
        completed_update, failed_update = invalidate_cached_run_results(completed_cases, failed_cases)

        return (
            config_state["selected_config"],
            config_state["settings_schema"],
            config_state["param_meta"],
            config_state["tunable_names"],
            config_state["tunable_default_ranges"],
            config_state["linked_parameter_groups"],
            config_state["settings_resolution"],
            build_right_pane(right_pane_data),
            multicol_next_id,
            multicol_row_order,
            multicol_rows,
            completed_update,
            failed_update,
        )

    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace",
            function_name="resetRunConfigControls",
        ),
        Output("run-config-reset-signal", "data"),
        Input({"type": "run-config-button", "name": ALL}, "n_clicks"),
        State({"type": "run-config-button", "name": ALL}, "id"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="runTab",
            function_name="syncAllParamRowClasses",
        ),
        Output({"type": "run-param-container", "file": ALL, "name": ALL}, "className"),
        Input("run-settings-resolution", "data"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "id"),
    )

    app.clientside_callback(
        ClientsideFunction(namespace="runTab", function_name="syncAllFlagRowClasses"),
        Output({"type": "run-flag-container", "name": ALL}, "className"),
        Input("run-settings-resolution", "data"),
        State({"type": "run-flag", "name": ALL}, "id"),
    )

    app.clientside_callback(
        ClientsideFunction(namespace="runTab", function_name="syncLinkedParamRowClasses"),
        Output({"type": "run-linked-param-container", "group": ALL}, "className"),
        Input("run-settings-resolution", "data"),
        State({"type": "run-linked-param-container", "group": ALL}, "id"),
    )

    @app.callback(
        Output("run-multicol-rows", "children"),
        Output("run-multicol-next-id", "data"),
        Output("run-multicol-row-order", "data"),
        Output("run-multicol-rows-state", "data", allow_duplicate=True),
        Input("run-multicol-add", "n_clicks"),
        Input({"type": "run-hr-remove", "index": ALL}, "n_clicks"),
        State("run-multicol-next-id", "data"),
        State("run-multicol-row-order", "data"),
        State("run-multicol-rows-state", "data"),
        State("run-tunable-names", "data"),
        State("run-tunable-default-ranges", "data"),
        State({"type": "run-hr-param", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def sync_multicol_rows(add_clicks, remove_clicks, next_id, row_order, rows_state, tunable_names, default_ranges, param_values):
        """Add or remove multicol rows without mirroring row values through the server."""
        current_order = list(row_order or [])
        current_rows = list(rows_state or [])
        trigger_id = callback_context.triggered_id

        if trigger_id == "run-multicol-add":
            if not add_clicks:
                return no_update, no_update, no_update, no_update
            row_id = int(next_id or 0)
            patch = Patch()
            row = default_multicol_row(row_id, tunable_names or [], default_ranges or {}, param_values or [])
            patch.append(build_multicol_row(row, tunable_names or []))
            return patch, row_id + 1, current_order + [row_id], current_rows + [row]

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "run-hr-remove":
            remove_id = trigger_id.get("index")
            if remove_id not in current_order:
                return no_update, no_update, no_update, no_update
            remove_pos = current_order.index(remove_id)
            if remove_pos >= len(remove_clicks or []) or not (remove_clicks or [])[remove_pos]:
                return no_update, no_update, no_update, no_update
            patch = Patch()
            del patch[remove_pos]
            current_order.remove(remove_id)
            updated_rows = [row for row in current_rows if row.get("id") != remove_id]
            return patch, no_update, current_order, updated_rows

        return no_update, no_update, no_update, no_update

    @app.callback(
        Output("run-multicol-rows-state", "data", allow_duplicate=True),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        Input({"type": "run-hr-min", "index": ALL}, "value"),
        Input({"type": "run-hr-max", "index": ALL}, "value"),
        Input({"type": "run-hr-npoints", "index": ALL}, "value"),
        State("run-multicol-row-order", "data"),
        State("run-multicol-rows-state", "data"),
        prevent_initial_call=True,
    )
    def persist_multicol_rows(param_values, min_values, max_values, npoint_values, row_order, saved_multicol_rows):
        """Persist live multicol values outside the right-pane subtree."""
        if callback_context.triggered_id is None:
            return no_update
        return multicol_rows_from_values(
            row_order,
            param_values,
            min_values,
            max_values,
            npoint_values,
            [],
            {},
            saved_multicol_rows,
        )

    @app.callback(
        Output({"type": "run-hr-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "run-hr-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-hr-min", "index": ALL}, "value"),
        State({"type": "run-hr-max", "index": ALL}, "value"),
        State("run-multicol-row-order", "data"),
        State("run-tunable-default-ranges", "data"),
        prevent_initial_call=True,
    )
    def autofill_multicol_bounds(param_values, min_values, max_values, row_order, default_ranges):
        """Fill multicol min/max from the selected parameter's default value."""
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "run-hr-param":
            return unchanged_multicol_bounds(param_values)

        row_id = trigger_id.get("index")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return unchanged_multicol_bounds(param_values)
        row_pos = current_order.index(row_id)

        values = list(param_values or [])
        if row_pos >= len(values):
            return unchanged_multicol_bounds(param_values)
        param_name = values[row_pos]
        derived_range = (default_ranges or {}).get(param_name)
        if not derived_range:
            return unchanged_multicol_bounds(param_values)

        updated_min = list(min_values or [])
        updated_max = list(max_values or [])
        while len(updated_min) < len(values):
            updated_min.append("")
        while len(updated_max) < len(values):
            updated_max.append("")

        updated_min[row_pos] = derived_range.get("min", "")
        updated_max[row_pos] = derived_range.get("max", "")
        return updated_min, updated_max

    @app.callback(
        Output("run-batch-size", "disabled"),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
    )
    def sync_batch_size_disabled(multicol_param_values):
        """Enable batch size only when a multicol parameter is selected."""
        return not any(str(value or "").strip() for value in (multicol_param_values or []))

    @app.callback(
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input({"type": "run-flag", "name": ALL}, "value"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        prevent_initial_call=True,
    )
    def invalidate_after_flag_edit(_flag_values, completed_cases, failed_cases):
        """Invalidate results after a user changes a boolean model flag."""
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "run-flag":
            return no_update, no_update
        return invalidate_cached_run_results(completed_cases, failed_cases)

    @app.callback(
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        Input("run-multicol-add", "n_clicks"),
        Input({"type": "run-hr-remove", "index": ALL}, "n_clicks"),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        Input({"type": "run-hr-min", "index": ALL}, "value"),
        Input({"type": "run-hr-max", "index": ALL}, "value"),
        Input({"type": "run-hr-npoints", "index": ALL}, "value"),
        Input("run-opt-max-iters", "value"),
        Input("run-opt-debug", "value"),
        Input("run-opt-dt-main", "value"),
        Input("run-opt-dt-rad", "value"),
        Input("run-opt-tout", "value"),
        Input("run-opt-out-dir", "value"),
        Input("run-batch-size", "value"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        prevent_initial_call=True,
    )
    def mark_nonboolean_settings_dirty(
        _param_values,
        _multicol_add_clicks,
        _multicol_remove_clicks,
        _multicol_params,
        _multicol_mins,
        _multicol_maxes,
        _multicol_npoints,
        _opt_max_iters,
        _opt_debug,
        _opt_dt_main,
        _opt_dt_rad,
        _opt_tout,
        _opt_out_dir,
        _batch_size,
        completed_cases,
        failed_cases,
    ):
        """Invalidate stale run results on any non-boolean settings change."""
        if callback_context.triggered_id is None:
            return no_update, no_update
        return invalidate_cached_run_results(completed_cases, failed_cases)

    @app.callback(
        Output({"type": "run-param", "file": "tunable", "name": ALL}, "disabled"),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        Input("run-settings-resolution", "data"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "id"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "disabled"),
        prevent_initial_call=True,
    )
    def sync_tunable_disabled(multicol_param_values, resolution, tunable_ids, current_disabled):
        """Disable hypergrid-owned, inactive, and derived tunable controls."""
        claimed = {str(value).strip() for value in (multicol_param_values or []) if str(value).strip()}
        states = dict((resolution or {}).get("parameter_states") or {})
        output = []
        for index, component_id in enumerate(tunable_ids or []):
            name = (component_id or {}).get("name")
            state = dict(states.get(name) or {})
            disabled = name in claimed or not is_independently_tunable(state)
            current = bool((current_disabled or [])[index]) if index < len(current_disabled or []) else False
            output.append(no_update if disabled == current else disabled)
        return output

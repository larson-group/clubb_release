from dash import ALL, Input, Output, State, dcc, html, no_update
import itertools
import math

from .plot_types.shared import (
    closest_column,
    load_compare_flag_values,
    load_compare_param_values,
    load_flag_values,
    load_param_values,
)
from .state import format_column_values


def _column_slider(ncols, value):
    """Render the shared column-index slider used when free column selection is allowed."""
    return html.Div(
        [
            dcc.Slider(
                id="plots-column-index-slider",
                min=1,
                max=ncols,
                value=min(max(int(value) + 1, 1), ncols),
                step=1,
                marks={1: "1", ncols: str(ncols)},
                tooltip={"always_visible": True, "placement": "bottom"},
            ),
        ],
        className="plots-slider-block",
    )


def _constant_param_grid(constant_params):
    """Render constant parameter values in the responsive two-column name/value grid."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div(name, className="plots-constant-param-name"),
                    html.Div(f"{value:g}", className="plots-constant-param-value"),
                ],
                className="plots-constant-param-row",
            )
            for name, value in constant_params
        ],
        className="plots-constant-param-grid",
    )


def _read_only_value_grid(entries):
    """Render arbitrary read-only name/value strings using the parameter grid."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div(name, className="plots-constant-param-name"),
                    html.Div(str(value), className="plots-constant-param-value"),
                ],
                className="plots-constant-param-row",
            )
            for name, value in entries
        ],
        className="plots-constant-param-grid",
    )


def _read_only_disclosure(title, content, count=None):
    """Return one collapsed-by-default optional-information section."""
    label = f"{title} ({count})" if count is not None else title
    return html.Details(
        [
            html.Summary(label, className="plots-readonly-summary"),
            html.Div(content, className="plots-readonly-content"),
        ],
        open=False,
        className="plots-readonly-disclosure",
    )


def _append_read_only_sections(children, constant_params, param_data):
    """Append constant tunables and saved model flags below interactive controls."""
    if constant_params:
        children.append(
            _read_only_disclosure(
                "Constant parameters",
                _constant_param_grid(constant_params),
                len(constant_params),
            )
        )
    flags = sorted((param_data.get("flags") or {}).items())
    flag_mismatches = list(param_data.get("flag_mismatches") or [])
    flag_rows = flags + flag_mismatches
    if flag_rows:
        children.append(
            _read_only_disclosure(
                "Configurable flags",
                _read_only_value_grid(flag_rows),
                len(flag_rows),
            )
        )
    else:
        message = (
            "No saved configurable flag namelist was found beside this output."
            if not param_data.get("has_flags")
            else "No configurable flags were recorded."
        )
        children.append(
            _read_only_disclosure("Configurable flags", html.Div(message))
        )
    return children


def _column_param_values(params, col_idx, ncols):
    """Return parameter values for the currently selected column."""
    column_idx = min(max(int(col_idx), 0), max(int(ncols) - 1, 0))
    values = []
    for name in sorted(params or {}):
        param_values = params.get(name, [])
        if len(param_values) == ncols:
            values.append((name, float(param_values[column_idx])))
    return values


def _param_error_grid(mismatched_params):
    """Render parameter rows whose values do not agree across compared outputs."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div(name, className="plots-constant-param-name"),
                    html.Div(message, className="plots-constant-param-value"),
                ],
                className="plots-constant-param-row",
            )
            for name, message in mismatched_params
        ],
        className="plots-constant-param-grid",
    )


def _regular_unique_values(values, *, atol=1.0e-10):
    """Return sorted unique values when they are evenly spaced."""
    unique_vals = sorted(set(float(value) for value in values))
    if len(unique_vals) <= 2:
        return unique_vals
    step = unique_vals[1] - unique_vals[0]
    scale = max(1.0, max(abs(value) for value in unique_vals))
    tolerance = atol * scale
    for left, right in zip(unique_vals[1:-1], unique_vals[2:]):
        if abs((right - left) - step) > tolerance:
            return None
    return unique_vals


def _param_value_groups(values, ncols):
    """Return sorted unique parameter values with the columns that use each value."""
    groups = {}
    for col_idx, value in enumerate(values[:ncols]):
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            continue
        groups.setdefault(numeric_value, []).append(col_idx)
    return [(value, groups[value]) for value in sorted(groups)]


def _slider_marks_for_groups(groups, active_index=None):
    """Return readable marks for a discrete parameter-value slider."""
    if len(groups) <= 8:
        return {idx: f"{value:g}" for idx, (value, _columns) in enumerate(groups)}
    marks = {
        0: f"{groups[0][0]:g}",
        len(groups) - 1: f"{groups[-1][0]:g}",
    }
    if active_index is not None and 0 <= active_index < len(groups):
        marks[int(active_index)] = f"{groups[int(active_index)][0]:g}"
    return marks


def _column_indices_for_filters(param_data, fixed_filters):
    """Return the column indices that match all active parameter filters."""
    if not param_data or not fixed_filters:
        return None, {}
    params = param_data.get("params") or {}
    ncols = max(int(param_data.get("ngrdcol") or 1), 1)
    selected = set(range(ncols))
    normalized_filters = {}
    for name, requested_index in fixed_filters.items():
        groups = _param_value_groups(params.get(name, []), ncols)
        if not groups:
            continue
        group_index = max(0, min(int(requested_index or 0), len(groups) - 1))
        group_value, group_columns = groups[group_index]
        selected &= set(group_columns)
        normalized_filters[name] = {"index": group_index, "value": group_value}
    if not normalized_filters:
        return None, {}
    return sorted(selected), normalized_filters


def _overplot_filter_row(name, values, ncols, column_filters):
    """Render one overplot parameter row with an optional fixed-value slider."""
    groups = _param_value_groups(values, ncols)
    filter_state = ((column_filters or {}).get("filters") or {}).get(name) or {}
    fixed = bool(filter_state)
    active_index = max(0, min(int(filter_state.get("index", 0) or 0), max(len(groups) - 1, 0)))
    value_summary = format_column_values(values)
    header = html.Div(
        [
            html.Div(name, className="plots-constant-param-name"),
            dcc.Checklist(
                id={"type": "plots-column-filter-enabled", "name": name},
                options=[{"label": "Fix", "value": "fixed"}],
                value=["fixed"] if fixed else [],
                labelStyle={"display": "inline-flex", "alignItems": "center", "gap": "4px", "whiteSpace": "nowrap"},
                inputStyle={"marginRight": "4px"},
                className="plots-column-filter-check",
            ),
        ],
        className="plots-column-filter-header",
    )
    children = [
        header,
        html.Div(value_summary, className="plots-constant-param-value"),
    ]
    if fixed and groups:
        children.append(
            html.Div(
                dcc.Slider(
                    id={"type": "plots-column-filter-slider", "name": name},
                    min=0,
                    max=len(groups) - 1,
                    value=active_index,
                    step=1,
                    marks=_slider_marks_for_groups(groups, active_index),
                ),
                className="plots-slider-block plots-column-filter-slider",
            )
        )
    return html.Div(children, className="plots-column-filter-row")


def _is_regular_hypergrid(params, varied_names, ncols):
    """Return whether varied column parameters form a complete regular hypergrid."""
    if not params or not varied_names:
        return False

    unique_by_name = {}
    for name in varied_names:
        values = params.get(name, [])
        if len(values) != ncols:
            return False
        unique_vals = _regular_unique_values(values)
        if unique_vals is None or len(unique_vals) < 2:
            return False
        unique_by_name[name] = unique_vals

    expected_count = math.prod(len(unique_by_name[name]) for name in varied_names)
    if expected_count != ncols:
        return False

    actual_points = {
        tuple(float(params[name][idx]) for name in varied_names)
        for idx in range(ncols)
    }
    expected_points = set(itertools.product(*(unique_by_name[name] for name in varied_names)))
    return actual_points == expected_points


def _varying_param_slider(name, current_val, unique_vals):
    """Render one parameter slider row for a column-varying tunable parameter."""
    return html.Div(
        [
            html.Div(
                name,
                style={
                    "minWidth": "72px",
                    "marginRight": "10px",
                    "fontWeight": "600",
                    "whiteSpace": "nowrap",
                },
            ),
            html.Div(
                dcc.Slider(
                    id={"type": "plots-param", "name": name},
                    min=min(unique_vals),
                    max=max(unique_vals),
                    value=current_val,
                    step=None,
                    marks={val: f"{val:g}" for val in unique_vals},
                    tooltip={"always_visible": True, "placement": "bottom"},
                ),
                style={"flex": "1 1 auto", "minWidth": 0},
            ),
        ],
        className="plots-slider-block",
        style={"display": "flex", "alignItems": "center"},
    )


def register_param_callbacks(app):
    """Register callbacks for loading, presenting, and selecting column parameters."""
    @app.callback(
        Output("plots-param-data", "data"),
        Output("plots-param-names", "data"),
        Input("plots-case-data", "data"),
    )
    def load_params(case_data):
        """Load tunable parameter data and split it into varying and constant groups."""
        if not case_data:
            return None, None
        ncols = max(int(case_data.get("columns_len") or 1), 1)
        if case_data.get("compare_mode"):
            files = case_data.get("files") or []
            params, has_clubb_params, has_param_names, mismatched_params, _per_file = load_compare_param_values(files)
            flags, flag_mismatches, has_flags = load_compare_flag_values(files)
            varied = [name for name, values in params.items() if len(values) == ncols and len(set(values)) > 1]
            allow_param_selection = _is_regular_hypergrid(params, varied, ncols)
            constant_params = [
                (name, params[name][0])
                for name in params.keys()
                if name not in varied and len(params[name]) == ncols and len(set(params[name])) <= 1
            ]
            return {
                "params": params,
                "ngrdcol": ncols,
                "compare_mode": True,
                "has_clubb_params": has_clubb_params,
                "has_param_names": has_param_names,
                "constant_params": constant_params,
                "mismatched_params": mismatched_params,
                "flags": flags,
                "flag_mismatches": flag_mismatches,
                "has_flags": has_flags,
                "allow_column_param_selection": allow_param_selection,
            }, varied
        files = case_data.get("files") or []
        params, has_clubb_params, has_param_names = load_param_values(files)
        flags, has_flags = load_flag_values(files)
        varied = [name for name, values in params.items() if len(values) == ncols and len(set(values)) > 1]
        allow_param_selection = _is_regular_hypergrid(params, varied, ncols)
        constant_params = [
            (name, params[name][0])
            for name in params.keys()
            if name not in varied and len(params[name]) == ncols and len(set(params[name])) <= 1
        ]
        return {
            "params": params,
            "ngrdcol": ncols,
            "compare_mode": False,
            "has_clubb_params": has_clubb_params,
            "has_param_names": has_param_names,
            "constant_params": constant_params,
            "mismatched_params": [],
            "flags": flags,
            "flag_mismatches": [],
            "has_flags": has_flags,
            "allow_column_param_selection": allow_param_selection,
        }, varied

    @app.callback(
        Output("plots-column-heading", "children"),
        Input("plots-selected-column", "data"),
        Input("plots-param-data", "data"),
        Input("plots-column-mode", "value"),
        Input("plots-column-filters", "data"),
        Input("plots-case-data", "data"),
    )
    def update_column_label(col_idx, param_data, column_mode, column_filters, _case_data):
        """Show the active column or overplot range in the column section header."""
        ncols = 1 if not param_data else param_data.get("ngrdcol", 1)
        if column_mode == "all":
            active_count = (column_filters or {}).get("active_count")
            if active_count is not None:
                return f"Columns: {active_count} of {ncols}"
            return f"Columns: 1 - {ncols}"
        return f"Column: {int(col_idx) + 1}"

    @app.callback(
        Output("plots-param-panel", "children"),
        Input("plots-param-data", "data"),
        Input("plots-param-names", "data"),
        Input("plots-selected-column", "data"),
        Input("plots-column-mode", "value"),
        Input("plots-case-data", "data"),
        Input("plots-column-filters", "data"),
    )
    def render_param_panel(param_data, param_names, col_idx, column_mode, case_data, column_filters):
        """Render the column-selection UI for the current case and plot mode."""
        if not param_data:
            return [html.Div("Select a case to enable column controls.")]
        ncols = param_data.get("ngrdcol", 1)
        compare_mode = bool(case_data and case_data.get("compare_mode"))
        has_clubb_params = bool(param_data.get("has_clubb_params"))
        has_param_names = bool(param_data.get("has_param_names"))
        constant_params = param_data.get("constant_params") or []
        mismatched_params = param_data.get("mismatched_params") or []
        allow_param_selection = bool(param_data.get("allow_column_param_selection"))
        children = []
        if compare_mode:
            children.append(html.Div("Compare mode uses a shared column index across matching outputs."))
        params = param_data["params"]
        displayed_varying_count = 0
        show_column_slider = (not has_clubb_params) or (not param_names)
        if compare_mode:
            show_column_slider = column_mode != "all" and (show_column_slider or bool(params))
        if column_mode != "all" and param_names and not allow_param_selection:
            show_column_slider = True
        if show_column_slider and column_mode != "all" and ncols > 1:
            children.append(_column_slider(ncols, col_idx))
        for name in param_names or []:
            if column_mode != "all" and not allow_param_selection:
                break
            values = params.get(name, [])
            if len(values) != ncols:
                continue
            unique_vals = sorted(set(values))
            if len(unique_vals) < 2:
                continue
            if column_mode == "all":
                displayed_varying_count += 1
                children.append(_overplot_filter_row(name, values, ncols, column_filters))
                continue
            current_val = values[min(max(int(col_idx), 0), ncols - 1)]
            children.append(_varying_param_slider(name, current_val, unique_vals))
        if column_mode == "all" and displayed_varying_count:
            children.insert(0, html.Div("Displayed varying parameters", style={"fontWeight": "600", "marginTop": "6px", "marginBottom": "8px"}))
        if mismatched_params:
            children.append(html.Div("Mismatched parameters", style={"fontWeight": "600", "marginTop": "12px", "marginBottom": "8px"}))
            children.append(_param_error_grid(mismatched_params))
        if param_names:
            if column_mode != "all" and not allow_param_selection:
                column_params = _column_param_values(params, col_idx, ncols)
                return _append_read_only_sections(children, column_params, param_data)
            return _append_read_only_sections(children, constant_params, param_data)
        if not has_clubb_params:
            children.append(html.Div("No clubb_params field found in the selected stats file."))
        elif not has_param_names:
            children.append(html.Div("clubb_params found, but param_name is missing from the selected stats file."))
        else:
            children.append(html.Div("All clubb_params are constant across columns."))
        return _append_read_only_sections(children, constant_params, param_data)

    @app.callback(
        Output("plots-column-filters", "data"),
        Input("plots-param-data", "data"),
        Input("plots-column-mode", "value"),
        Input({"type": "plots-column-filter-enabled", "name": ALL}, "value"),
        Input({"type": "plots-column-filter-slider", "name": ALL}, "value"),
        State({"type": "plots-column-filter-enabled", "name": ALL}, "id"),
        State({"type": "plots-column-filter-slider", "name": ALL}, "id"),
    )
    def update_column_filters(param_data, column_mode, enabled_values, slider_values, enabled_ids, slider_ids):
        """Compute the active overplot column subset from fixed-parameter controls."""
        empty_filter = {"indices": None, "filters": {}, "active_count": None}
        if column_mode != "all" or not param_data:
            return empty_filter
        slider_by_name = {
            item.get("name"): value
            for item, value in zip(slider_ids or [], slider_values or [])
            if isinstance(item, dict) and item.get("name") is not None
        }
        requested_filters = {}
        for item, values in zip(enabled_ids or [], enabled_values or []):
            if not isinstance(item, dict):
                continue
            name = item.get("name")
            if not name or "fixed" not in (values or []):
                continue
            requested_filters[name] = slider_by_name.get(name, 0)
        indices, normalized_filters = _column_indices_for_filters(param_data, requested_filters)
        if not normalized_filters:
            return empty_filter
        return {
            "indices": indices,
            "filters": normalized_filters,
            "active_count": len(indices or []),
        }

    @app.callback(
        Output("plots-selected-column", "data", allow_duplicate=True),
        Input("plots-column-index-slider", "value"),
        State("plots-param-data", "data"),
        prevent_initial_call=True,
    )
    def update_column_from_slider(value, param_data):
        """Map the visible column slider value back to the zero-based selected column."""
        if value is None or not param_data:
            return no_update
        ncols = param_data.get("ngrdcol", 1)
        return max(0, min(int(value) - 1, ncols - 1))

    @app.callback(
        Output("plots-selected-column", "data", allow_duplicate=True),
        Input({"type": "plots-param", "name": ALL}, "value"),
        State("plots-param-names", "data"),
        State("plots-param-data", "data"),
        State("plots-column-mode", "value"),
        prevent_initial_call=True,
    )
    def update_column_from_params(values, names, param_data, column_mode):
        """Choose the nearest matching column from the current parameter slider values."""
        if column_mode != "single" or not param_data or not names or not param_data.get("allow_column_param_selection"):
            return no_update
        params = param_data["params"]
        selection = {}
        for name, value in zip(names, values or []):
            if name in params and value is not None:
                selection[name] = value
        return closest_column(params, selection)

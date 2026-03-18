from dash import ALL, Input, Output, State, dcc, html, no_update

from .plot_types.shared import closest_column, load_compare_param_values, load_param_values
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
            params, has_clubb_params, has_param_names, mismatched_params, _per_file = load_compare_param_values(case_data.get("files") or [])
            varied = [name for name, values in params.items() if len(values) == ncols and len(set(values)) > 1]
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
                "allow_column_param_selection": bool(params),
            }, varied
        params, has_clubb_params, has_param_names = load_param_values(case_data.get("files") or [])
        varied = [name for name, values in params.items() if len(values) == ncols and len(set(values)) > 1]
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
            "allow_column_param_selection": True,
        }, varied

    @app.callback(
        Output("plots-column-heading", "children"),
        Input("plots-selected-column", "data"),
        Input("plots-param-data", "data"),
        Input("plots-column-mode", "value"),
        Input("plots-case-data", "data"),
    )
    def update_column_label(col_idx, param_data, column_mode, _case_data):
        """Show the active column or overplot range in the column section header."""
        ncols = 1 if not param_data else param_data.get("ngrdcol", 1)
        if column_mode == "all":
            return f"Columns: 1 - {ncols}"
        return f"Column: {int(col_idx) + 1}"

    @app.callback(
        Output("plots-param-panel", "children"),
        Input("plots-param-data", "data"),
        Input("plots-param-names", "data"),
        Input("plots-selected-column", "data"),
        Input("plots-column-mode", "value"),
        Input("plots-case-data", "data"),
    )
    def render_param_panel(param_data, param_names, col_idx, column_mode, case_data):
        """Render the column-selection UI for the current case and plot mode."""
        if not param_data:
            return [html.Div("Select a case to enable column controls.")]
        ncols = param_data.get("ngrdcol", 1)
        compare_mode = bool(case_data and case_data.get("compare_mode"))
        has_clubb_params = bool(param_data.get("has_clubb_params"))
        has_param_names = bool(param_data.get("has_param_names"))
        constant_params = param_data.get("constant_params") or []
        mismatched_params = param_data.get("mismatched_params") or []
        children = []
        if compare_mode:
            children.append(html.Div("Compare mode uses a shared column index across matching outputs."))
        params = param_data["params"]
        show_column_slider = (not has_clubb_params) or (not param_names)
        if compare_mode:
            show_column_slider = column_mode != "all" and (show_column_slider or bool(params))
        if show_column_slider and column_mode != "all" and ncols > 1:
            children.append(_column_slider(ncols, col_idx))
        for name in param_names or []:
            values = params.get(name, [])
            if len(values) != ncols:
                continue
            unique_vals = sorted(set(values))
            if len(unique_vals) < 2:
                continue
            if column_mode == "all":
                children.append(
                    html.Div(
                        [
                            html.Div(name, className="plots-constant-param-name"),
                            html.Div(format_column_values(values), className="plots-constant-param-value"),
                        ],
                        className="plots-constant-param-row",
                    )
                )
                continue
            current_val = values[min(max(int(col_idx), 0), ncols - 1)]
            children.append(_varying_param_slider(name, current_val, unique_vals))
        if column_mode == "all" and children:
            children.insert(0, html.Div("Displayed varying parameters", style={"fontWeight": "600", "marginTop": "6px", "marginBottom": "8px"}))
        if mismatched_params:
            children.append(html.Div("Mismatched parameters", style={"fontWeight": "600", "marginTop": "12px", "marginBottom": "8px"}))
            children.append(_param_error_grid(mismatched_params))
        if param_names:
            if constant_params:
                children.append(html.Div("Constant parameters", style={"fontWeight": "600", "marginTop": "6px", "marginBottom": "8px"}))
                children.append(_constant_param_grid(constant_params))
            return children
        if not has_clubb_params:
            children.append(html.Div("No clubb_params field found in the selected stats file."))
        elif not has_param_names:
            children.append(html.Div("clubb_params found, but param_name is missing from the selected stats file."))
        else:
            children.append(html.Div("All clubb_params are constant across columns."))
            if constant_params:
                children.append(_constant_param_grid(constant_params))
        return children

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

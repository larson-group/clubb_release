import json
import os
import glob

import numpy as np
import plotly.graph_objects as go
from dash import dcc, html, Input, Output, State, MATCH, ALL, callback_context, no_update
from netCDF4 import Dataset

from case_definitions import load_case_definitions

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTPUT_DIR = os.path.join(REPO_ROOT, "output")
PARAM_FILE_CANDIDATES = [
    os.path.join(REPO_ROOT, "run_scripts", "clubb.in"),
    os.path.join(REPO_ROOT, "run_scripts", "clubb_params_multi_col.in"),
    os.path.join(REPO_ROOT, "run_scripts", "tmp_dup_params.in"),
]

Z_DIM_NAMES = {"z", "altitude", "height", "lev"}
T_DIM_NAMES = {"t", "time"}
X_DIM_NAMES = {"x", "longitude", "lon"}
Y_DIM_NAMES = {"y", "latitude", "lat"}
COL_DIM_NAMES = {"columns", "column", "col", "ngrdcol"}

DEFAULT_VARS = ["thlm", "um", "wp2"]


def _find_dim(dim_names, candidates):
    for dim_name in dim_names:
        if dim_name.lower() in candidates:
            return dim_name
    return None


def _collect_vars(ds, require_columns=True):
    vars_with_z = []
    for name, var in ds.variables.items():
        if name in ds.dimensions:
            continue
        z_dim = _find_dim(var.dimensions, Z_DIM_NAMES)
        if z_dim is None:
            continue
        if require_columns:
            col_dim = _find_dim(var.dimensions, COL_DIM_NAMES)
            if col_dim is None:
                continue
        vars_with_z.append(name)
    return sorted(vars_with_z)


class DatasetInfo:
    def __init__(self, path):
        self.path = path
        self.ds = Dataset(path, "r")
        self.vars = _collect_vars(self.ds)
        self.var_info = {}
        for name in self.vars:
            dims = list(self.ds.variables[name].dimensions)
            self.var_info[name] = {
                "z_dim": _find_dim(dims, Z_DIM_NAMES),
                "t_dim": _find_dim(dims, T_DIM_NAMES),
                "col_dim": _find_dim(dims, COL_DIM_NAMES),
            }
        self.time_dim = _find_dim(list(self.ds.dimensions.keys()), T_DIM_NAMES)
        self.time_len = 1
        self.time_values = None
        self.time_units = None
        if self.time_dim:
            self.time_len = len(self.ds.dimensions[self.time_dim])
            if self.time_dim in self.ds.variables:
                time_var = self.ds.variables[self.time_dim]
                self.time_values = np.asarray(time_var[:])
                self.time_units = getattr(time_var, "units", None)
        self.columns_dim = _find_dim(list(self.ds.dimensions.keys()), COL_DIM_NAMES)
        self.columns_len = 1
        if self.columns_dim:
            self.columns_len = len(self.ds.dimensions[self.columns_dim])

    def close(self):
        self.ds.close()


class DatasetCollection:
    def __init__(self, paths):
        self.paths = paths
        self.datasets = [DatasetInfo(path) for path in paths]
        self.var_to_ds = {}
        for dataset in self.datasets:
            for var_name in dataset.vars:
                if var_name not in self.var_to_ds:
                    self.var_to_ds[var_name] = dataset
        self.vars = sorted(self.var_to_ds.keys())
        self.time_len = max((ds.time_len for ds in self.datasets), default=1)
        self.columns_len = max((ds.columns_len for ds in self.datasets), default=1)

    def dataset_for_var(self, var_name):
        return self.var_to_ds.get(var_name)

    def primary_time_source(self):
        if not self.datasets:
            return None
        return max(self.datasets, key=lambda ds: ds.time_len)

    def close(self):
        for ds in self.datasets:
            ds.close()


def _time_units_factor(units):
    if not units:
        return 1.0
    unit = units.split()[0].lower()
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return 1.0
    if unit in {"min", "mins", "minute", "minutes"}:
        return 60.0
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return 3600.0
    return 1.0


def _time_label(ds_info, t_range):
    t0, t1 = t_range
    if ds_info is None or ds_info.time_values is None:
        return f"Time: {t0} - {t1} s"
    tmax = len(ds_info.time_values) - 1
    t0 = max(0, min(t0, tmax))
    t1 = max(0, min(t1, tmax))
    factor = _time_units_factor(ds_info.time_units)
    t0_sec = float(ds_info.time_values[t0]) * factor
    t1_sec = float(ds_info.time_values[t1]) * factor
    return f"Time: {t0_sec:.0f} - {t1_sec:.0f} s"


def _minutes_to_slider_range(ds_info, start_min, end_min, slider_max):
    if start_min is None or end_min is None:
        return [1, slider_max]
    if ds_info is None or ds_info.time_values is None or len(ds_info.time_values) == 0:
        start_idx = max(1, min(int(start_min), slider_max))
        end_idx = max(1, min(int(end_min), slider_max))
    else:
        factor = _time_units_factor(ds_info.time_units)
        time_sec = ds_info.time_values * factor
        start_sec = float(start_min) * 60.0
        end_sec = float(end_min) * 60.0
        start_idx = int(np.argmin(np.abs(time_sec - start_sec))) + 1
        end_idx = int(np.argmin(np.abs(time_sec - end_sec))) + 1
        start_idx = max(1, min(start_idx, slider_max))
        end_idx = max(1, min(end_idx, slider_max))
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    return [start_idx, end_idx]


def _time_avg_profile(ds_info, var_name, t_range, col_index):
    var = ds_info.ds.variables[var_name]
    dims = list(var.dimensions)
    z_dim = ds_info.var_info[var_name]["z_dim"]
    t_dim = ds_info.var_info[var_name]["t_dim"]
    col_dim = ds_info.var_info[var_name]["col_dim"]

    t0, t1 = t_range
    if t_dim is None:
        t0, t1 = 0, 0
    else:
        tmax = len(ds_info.ds.dimensions[t_dim]) - 1
        t0 = max(0, min(t0, tmax))
        t1 = max(0, min(t1, tmax))
        if t1 < t0:
            t0, t1 = t1, t0

    slices = []
    kept_dims = []
    for dim in dims:
        if dim == t_dim:
            slices.append(slice(t0, t1 + 1))
            kept_dims.append(dim)
        elif dim == z_dim:
            slices.append(slice(None))
            kept_dims.append(dim)
        elif dim == col_dim:
            slices.append(col_index)
        elif dim.lower() in X_DIM_NAMES or dim.lower() in Y_DIM_NAMES:
            slices.append(0)
        else:
            slices.append(0)

    data = var[tuple(slices)]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    data = np.asarray(data)

    if t_dim in kept_dims:
        t_axis = kept_dims.index(t_dim)
        mean = np.mean(data, axis=t_axis)
        kept_dims.pop(t_axis)
    else:
        mean = data

    if mean.ndim > 1:
        mean = np.reshape(mean, (mean.shape[0], -1))[:, 0]

    z_vals = ds_info.ds.variables[z_dim][:]
    if np.ma.isMaskedArray(z_vals):
        z_vals = z_vals.filled(np.nan)

    return np.asarray(z_vals), mean


def _make_figure(z_vals, profile, title, field_name, z_units):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=profile,
            y=z_vals,
            mode="lines",
            name=field_name,
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title=field_name,
        yaxis_title=f"{z_units}" if z_units else "z",
        margin=dict(l=60, r=20, t=40, b=50),
    )
    return fig


def _scan_output_cases():
    cases = {}
    for path in glob.glob(os.path.join(OUTPUT_DIR, "*_multi_col_*.nc")):
        base = os.path.basename(path)
        if "_multi_col_" not in base:
            continue
        case_name = base.split("_multi_col_")[0]
        cases.setdefault(case_name, []).append(path)
    for files in cases.values():
        files.sort()
    return cases


def _parse_param_file(path, ncols):
    params = {}
    in_params = False
    current = None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.split("!")[0].strip()
            if not line:
                continue
            if line.startswith("&"):
                in_params = line.lower().startswith("&clubb_params_nl")
                current = None
                continue
            if line.startswith("/"):
                in_params = False
                current = None
                continue
            if not in_params:
                continue
            if "=" in line:
                name, rest = line.split("=", 1)
                current = name.strip()
                params[current] = _parse_param_values(rest)
            elif current:
                params[current] += _parse_param_values(line)

    for name, values in params.items():
        if not values:
            continue
        if len(values) < ncols:
            values.extend([values[-1]] * (ncols - len(values)))
        if len(values) > ncols:
            params[name] = values[:ncols]
    return params


def _parse_param_values(text):
    tokens = []
    for part in text.replace(",", " ").split():
        try:
            tokens.append(float(part))
        except ValueError:
            continue
    return tokens


def _load_param_values(ncols):
    for path in PARAM_FILE_CANDIDATES:
        if os.path.exists(path):
            params = _parse_param_file(path, ncols)
            if params:
                return params
    return {}


def _closest_column(param_values, selection):
    if not param_values or not selection:
        return 0
    ncols = len(next(iter(param_values.values())))
    best_idx = 0
    best_score = None
    for idx in range(ncols):
        score = 0.0
        for name, target in selection.items():
            col_val = param_values[name][idx]
            score += (col_val - target) ** 2
        if best_score is None or score < best_score:
            best_score = score
            best_idx = idx
    return best_idx


def _plot_card(idx, default_var):
    id_prefix = "multicol"
    return html.Div(
        [
            html.Div(
                [
                    html.Label("Field"),
                    dcc.Dropdown(
                        id={"type": f"{id_prefix}-var-select", "index": idx},
                        options=[{"label": v, "value": v} for v in default_var["options"]],
                        value=default_var["value"],
                        clearable=False,
                    ),
                ],
                style={"width": "55%", "display": "inline-block"},
            ),
            html.Div(
                [
                    dcc.Graph(id={"type": f"{id_prefix}-graph", "index": idx}),
                ],
                style={
                    "border": "1px solid #bdbdbd",
                    "borderRadius": "6px",
                    "padding": "6px",
                    "backgroundColor": "#ffffff",
                },
            ),
        ],
        style={
            "padding": "10px",
            "border": "1px solid #c9c9c9",
            "borderRadius": "6px",
            "backgroundColor": "#f9f9f9",
        },
    )


def build_tab(app):
    case_defs = load_case_definitions()
    case_by_name = {case["name"]: case for case in case_defs if case.get("name")}
    output_cases = _scan_output_cases()

    case_buttons = []
    for case in case_defs:
        name = case.get("name")
        if not name:
            continue
        files = output_cases.get(name, [])
        available = bool(files)
        if available:
            style = {
                "backgroundColor": "#2563eb",
                "color": "#ffffff",
                "border": "none",
                "padding": "6px 10px",
                "margin": "4px",
                "borderRadius": "4px",
                "cursor": "pointer",
            }
        else:
            style = {
                "backgroundColor": "#c9c9c9",
                "color": "#5f5f5f",
                "border": "none",
                "padding": "6px 10px",
                "margin": "4px",
                "borderRadius": "4px",
                "cursor": "not-allowed",
            }
        case_buttons.append(
            html.Button(
                name,
                id={"type": "multicol-case-button", "name": name},
                disabled=not available,
                style=style,
            )
        )

    layout = html.Div(
        [
            html.Div(
                [
                    html.Div("Cases:"),
                    html.Div(case_buttons),
                ],
                style={"marginBottom": "12px"},
            ),
            dcc.Store(id="multicol-case-data", data=None),
            dcc.Store(id="multicol-param-data", data=None),
            dcc.Store(id="multicol-param-names", data=None),
            dcc.Store(id="multicol-selected-column", data=0),
            dcc.Store(id="multicol-plot-ids", data=[]),
            dcc.Interval(id="multicol-case-refresh", interval=3000, n_intervals=0),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(id="multicol-time-label"),
                            dcc.RangeSlider(
                                id="multicol-time-range",
                                min=1,
                                max=1,
                                value=[1, 1],
                                step=1,
                                allowCross=False,
                                marks={1: "1"},
                                tooltip={"always_visible": True, "placement": "bottom"},
                            ),
                            html.Button("Add plot", id="multicol-add-plot", n_clicks=0),
                            html.Div(
                                id="multicol-plot-container",
                                style={
                                    "display": "grid",
                                    "gridTemplateColumns": "repeat(auto-fill, minmax(400px, 1fr))",
                                    "gap": "16px",
                                    "marginTop": "16px",
                                },
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            html.H4("Parameters"),
                            html.Div(id="multicol-column-label"),
                            html.Div(
                                id="multicol-param-panel",
                                style={
                                    "overflowY": "auto",
                                    "maxHeight": "700px",
                                    "paddingRight": "8px",
                                },
                            ),
                        ],
                        style={"borderLeft": "1px solid #d0d0d0", "paddingLeft": "16px"},
                    ),
                ],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "1fr 360px",
                    "gap": "16px",
                },
            ),
        ],
        style={"padding": "10px"},
    )

    @app.callback(
        Output("multicol-case-data", "data"),
        Output("multicol-plot-ids", "data"),
        Output("multicol-selected-column", "data"),
        Input({"type": "multicol-case-button", "name": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _select_case(_clicks):
        if not callback_context.triggered:
            return no_update, no_update, no_update
        button_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        try:
            case_name = json.loads(button_id).get("name")
        except Exception:
            return no_update, no_update, no_update
        current_defs = load_case_definitions()
        current_by_name = {case["name"]: case for case in current_defs if case.get("name")}
        current_output = _scan_output_cases()
        case = current_by_name.get(case_name)
        files = current_output.get(case_name, [])
        if not case or not files:
            return no_update, no_update, no_update
        case_data = {
            "name": case_name,
            "files": files,
            "start_time": case.get("start_time"),
            "end_time": case.get("end_time"),
        }
        return case_data, [0, 1, 2], 0

    @app.callback(
        Output({"type": "multicol-case-button", "name": ALL}, "style"),
        Output({"type": "multicol-case-button", "name": ALL}, "disabled"),
        Input("multicol-case-refresh", "n_intervals"),
        State({"type": "multicol-case-button", "name": ALL}, "id"),
    )
    def _refresh_case_buttons(_tick, ids):
        if not ids:
            return [], []
        output_cases = _scan_output_cases()
        styles = []
        disabled = []
        for case_id in ids:
            name = case_id.get("name")
            available = bool(output_cases.get(name, []))
            if available:
                style = {
                    "backgroundColor": "#2563eb",
                    "color": "#ffffff",
                    "border": "none",
                    "padding": "6px 10px",
                    "margin": "4px",
                    "borderRadius": "4px",
                    "cursor": "pointer",
                }
            else:
                style = {
                    "backgroundColor": "#c9c9c9",
                    "color": "#5f5f5f",
                    "border": "none",
                    "padding": "6px 10px",
                    "margin": "4px",
                    "borderRadius": "4px",
                    "cursor": "not-allowed",
                }
            styles.append(style)
            disabled.append(not available)
        return styles, disabled

    @app.callback(
        Output("multicol-param-data", "data"),
        Output("multicol-param-names", "data"),
        Input("multicol-case-data", "data"),
        prevent_initial_call=True,
    )
    def _load_params(case_data):
        if not case_data:
            return no_update, no_update
        collection = DatasetCollection(case_data["files"])
        try:
            ncols = max(collection.columns_len, 1)
        finally:
            collection.close()
        params = _load_param_values(ncols)
        varied = [
            name
            for name, values in params.items()
            if len(set(values)) > 1 and len(values) == ncols
        ]
        param_data = {"params": params, "ngrdcol": ncols}
        return param_data, varied

    @app.callback(
        Output("multicol-time-range", "max"),
        Output("multicol-time-range", "value"),
        Output("multicol-time-range", "marks"),
        Input("multicol-case-data", "data"),
    )
    def _configure_time(case_data):
        if not case_data:
            return 1, [1, 1], {1: "1"}
        collection = DatasetCollection(case_data["files"])
        try:
            slider_max = max(collection.time_len, 1)
            time_source = collection.primary_time_source()
            slider_value = _minutes_to_slider_range(
                time_source,
                case_data.get("start_time"),
                case_data.get("end_time"),
                slider_max,
            )
        finally:
            collection.close()
        marks = {1: "1", slider_max: str(slider_max)}
        return slider_max, slider_value, marks

    @app.callback(
        Output("multicol-time-label", "children"),
        Input("multicol-time-range", "value"),
        Input("multicol-case-data", "data"),
    )
    def _update_time_label(t_range, case_data):
        if not case_data:
            return ""
        collection = DatasetCollection(case_data["files"])
        try:
            time_source = collection.primary_time_source()
            t_range_zero = [max(0, t_range[0] - 1), max(0, t_range[1] - 1)]
            return _time_label(time_source, t_range_zero)
        finally:
            collection.close()

    @app.callback(
        Output("multicol-column-label", "children"),
        Input("multicol-selected-column", "data"),
    )
    def _update_column_label(col_idx):
        return f"Selected column: {int(col_idx) + 1}"

    @app.callback(
        Output("multicol-param-panel", "children"),
        Input("multicol-param-data", "data"),
        Input("multicol-param-names", "data"),
        Input("multicol-selected-column", "data"),
    )
    def _render_params(param_data, param_names, col_idx):
        if not param_data or not param_names:
            return [html.Div("No varying parameters found.")]
        params = param_data["params"]
        ncols = param_data.get("ngrdcol", 1)
        sliders = []
        for name in param_names:
            values = params.get(name, [])
            if len(values) != ncols:
                continue
            unique_vals = sorted(set(values))
            if len(unique_vals) < 2:
                continue
            marks = {val: f"{val:g}" for val in unique_vals}
            current_val = values[min(max(int(col_idx), 0), ncols - 1)]
            sliders.append(
                html.Div(
                    [
                        html.Label(name),
                        dcc.Slider(
                            id={"type": "multicol-param", "name": name},
                            min=min(unique_vals),
                            max=max(unique_vals),
                            value=current_val,
                            step=None,
                            marks=marks,
                            tooltip={"always_visible": True, "placement": "bottom"},
                        ),
                    ],
                    style={"marginBottom": "16px"},
                )
            )
        return sliders

    @app.callback(
        Output("multicol-selected-column", "data", allow_duplicate=True),
        Input({"type": "multicol-param", "name": ALL}, "value"),
        State("multicol-param-names", "data"),
        State("multicol-param-data", "data"),
        State("multicol-selected-column", "data"),
        prevent_initial_call=True,
    )
    def _update_selected_column(values, names, param_data, current_col):
        if not param_data or not names:
            return current_col
        params = param_data["params"]
        selection = {}
        for name, value in zip(names, values):
            if name in params and value is not None:
                selection[name] = value
        return _closest_column(params, selection)

    @app.callback(
        Output("multicol-plot-ids", "data", allow_duplicate=True),
        Input("multicol-add-plot", "n_clicks"),
        State("multicol-plot-ids", "data"),
        prevent_initial_call=True,
    )
    def _add_plot(n_clicks, plot_ids):
        if plot_ids is None:
            plot_ids = []
        next_id = max(plot_ids) + 1 if plot_ids else 0
        return plot_ids + [next_id]

    @app.callback(
        Output("multicol-plot-container", "children"),
        Input("multicol-plot-ids", "data"),
        Input("multicol-case-data", "data"),
    )
    def _render_plots(plot_ids, case_data):
        if not case_data:
            return [html.Div("Select a case to view plots.")]
        collection = DatasetCollection(case_data["files"])
        try:
            if not collection.vars:
                return [html.Div("No variables with columns dimension found.")]
            defaults = []
            for idx in range(len(plot_ids or [])):
                if idx < len(DEFAULT_VARS) and DEFAULT_VARS[idx] in collection.vars:
                    defaults.append(DEFAULT_VARS[idx])
                else:
                    defaults.append(collection.vars[0])
            cards = []
            for idx, default_name in zip(plot_ids or [], defaults):
                cards.append(
                    _plot_card(
                        idx,
                        {"value": default_name, "options": collection.vars},
                    )
                )
            return cards
        finally:
            collection.close()

    @app.callback(
        Output({"type": "multicol-graph", "index": MATCH}, "figure"),
        Input({"type": "multicol-var-select", "index": MATCH}, "value"),
        Input("multicol-time-range", "value"),
        Input("multicol-case-data", "data"),
        Input("multicol-selected-column", "data"),
    )
    def _update_graph(var_name, t_range, case_data, col_idx):
        if not case_data or var_name is None:
            return go.Figure()
        collection = DatasetCollection(case_data["files"])
        try:
            ds = collection.dataset_for_var(var_name) or collection.primary_time_source()
            if ds is None:
                return go.Figure()
            t_range_zero = [max(0, t_range[0] - 1), max(0, t_range[1] - 1)]
            z_vals, profile = _time_avg_profile(ds, var_name, t_range_zero, int(col_idx))
            z_units = getattr(ds.ds.variables[ds.var_info[var_name]["z_dim"]], "units", "")
            fig = _make_figure(z_vals, profile, var_name, var_name, z_units)
            fig.update_layout(width=400, height=400)
            return fig
        finally:
            collection.close()

    return dcc.Tab(label="Multi-Column", value="multicol", children=layout)

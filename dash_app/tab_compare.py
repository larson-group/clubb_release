import base64
import json
import os
import tempfile

import numpy as np
import plotly.graph_objects as go
from dash import dcc, html, Input, Output, State, MATCH, ALL, callback_context, no_update
from netCDF4 import Dataset

from case_definitions import load_case_definitions

Z_DIM_NAMES = {"z", "altitude", "height", "lev"}
T_DIM_NAMES = {"t", "time"}
X_DIM_NAMES = {"x", "longitude", "lon"}
Y_DIM_NAMES = {"y", "latitude", "lat"}


def _find_dim(dim_names, candidates):
    for dim_name in dim_names:
        if dim_name.lower() in candidates:
            return dim_name
    return None


def _convert_units(values, units):
    if units is None:
        return values
    units = units.strip()
    if units == "g/kg":
        return values / 1000.0
    if units == "K/day":
        return values / 86400.0
    return values


def _collect_vars(ds):
    vars_with_z = []
    for name, var in ds.variables.items():
        if name in ds.dimensions:
            continue
        z_dim = _find_dim(var.dimensions, Z_DIM_NAMES)
        if z_dim is not None:
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

    def dataset_for_var(self, var_name):
        return self.var_to_ds.get(var_name)

    def primary_time_source(self):
        if not self.datasets:
            return None
        return max(self.datasets, key=lambda ds: ds.time_len)

    def close(self):
        for ds in self.datasets:
            ds.close()


def _save_upload(contents, filename, prefix):
    if contents is None:
        return None
    _, content_string = contents.split(",", 1)
    data = base64.b64decode(content_string)
    safe_name = os.path.basename(filename) if filename else "upload.nc"
    tmp = tempfile.NamedTemporaryFile(
        delete=False,
        prefix=f"clubb_dash_{prefix}_",
        suffix=f"_{safe_name}",
        dir="/tmp",
    )
    tmp.write(data)
    tmp.close()
    return tmp.name


def _save_uploads(contents_list, filename_list, prefix):
    if contents_list is None:
        return None
    paths = []
    for contents, filename in zip(contents_list, filename_list or []):
        path = _save_upload(contents, filename, prefix)
        if path:
            paths.append(path)
    return paths


def _format_selected_label(paths):
    if not paths:
        return ""
    if isinstance(paths, list):
        names = ", ".join(paths)
    else:
        names = paths
    return names


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


def _time_avg_profile(ds_info, var_name, t_range):
    var = ds_info.ds.variables[var_name]
    dims = list(var.dimensions)
    z_dim = ds_info.var_info[var_name]["z_dim"]
    t_dim = ds_info.var_info[var_name]["t_dim"]

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
        elif dim.lower() in X_DIM_NAMES or dim.lower() in Y_DIM_NAMES:
            slices.append(0)
        else:
            slices.append(0)

    data = var[tuple(slices)]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    data = np.asarray(data)
    data = _convert_units(data, getattr(var, "units", None))

    if t_dim in kept_dims:
        t_axis = kept_dims.index(t_dim)
        mean = np.mean(data, axis=t_axis)
        mean2 = np.mean(data**2, axis=t_axis)
        kept_dims.pop(t_axis)
    else:
        mean = data
        mean2 = data**2

    if mean.ndim > 1:
        mean = np.reshape(mean, (mean.shape[0], -1))[:, 0]
        mean2 = np.reshape(mean2, (mean2.shape[0], -1))[:, 0]

    z_vals = ds_info.ds.variables[z_dim][:]
    if np.ma.isMaskedArray(z_vals):
        z_vals = z_vals.filled(np.nan)

    return np.asarray(z_vals), mean, mean2


def _time_label(ds_info, t_range):
    t0, t1 = t_range
    if ds_info.time_values is None:
        return f"Time: {t0} - {t1} s"
    tmax = len(ds_info.time_values) - 1
    t0 = max(0, min(t0, tmax))
    t1 = max(0, min(t1, tmax))
    factor = _time_units_factor(ds_info.time_units)
    t0_sec = float(ds_info.time_values[t0]) * factor
    t1_sec = float(ds_info.time_values[t1]) * factor
    return f"Time: {t0_sec:.0f} - {t1_sec:.0f} s"


def _time_label_point(ds_info, t_index):
    if ds_info.time_values is None:
        return f"Time: {t_index} s"
    tmax = len(ds_info.time_values) - 1
    idx = max(0, min(t_index, tmax))
    factor = _time_units_factor(ds_info.time_units)
    t_sec = float(ds_info.time_values[idx]) * factor
    return f"Time: {t_sec:.0f} s"


def _minutes_to_slider_range(ds_info, start_min, end_min, slider_max):
    if start_min is None or end_min is None:
        return [1, slider_max]
    if ds_info.time_values is None or len(ds_info.time_values) == 0:
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


def _make_figure(z_les, les_profile, z_clubb, clubb_profile, title, field_name, z_units):
    fig = go.Figure()
    if z_les is not None and les_profile is not None:
        fig.add_trace(
            go.Scatter(
                x=les_profile,
                y=z_les,
                mode="lines",
                name="les",
                line={"dash": "dash"},
            )
        )
    fig.add_trace(
        go.Scatter(
            x=clubb_profile,
            y=z_clubb,
            mode="lines",
            name="clubb",
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title=field_name,
        yaxis_title=f"{z_units}" if z_units else "z",
        margin=dict(l=60, r=20, t=40, b=50),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def _compute_error(obs_profile, obs_z, model_profile, model_profile2, model_z, field_name):
    obs_profile = np.asarray(obs_profile)
    model_profile = np.asarray(model_profile)
    model_profile2 = np.asarray(model_profile2)
    obs_z = np.asarray(obs_z)
    model_z = np.asarray(model_z)

    if obs_profile.shape != model_profile.shape:
        valid = np.isfinite(obs_z) & np.isfinite(obs_profile)
        if np.count_nonzero(valid) >= 2:
            obs_profile = np.interp(model_z, obs_z[valid], obs_profile[valid])
        else:
            obs_profile = np.full_like(model_profile, np.nan)

    les_minmax = np.nanmax(obs_profile) - np.nanmin(obs_profile)
    if field_name == "cloud_frac":
        les_minmax = max(0.01, les_minmax)
    elif field_name == "rcm":
        les_minmax = max(1e-6, les_minmax)
    return np.nansum(
        (model_profile2 - 2.0 * model_profile * obs_profile + obs_profile**2)
        / (les_minmax**2)
    )


def _case_clubb_files(case):
    clubb_files = case.get("clubb_files") or []
    clubb_file = case.get("clubb_file")
    if clubb_file:
        clubb_files = clubb_files + [clubb_file]
    seen = set()
    results = []
    for path in clubb_files:
        if path and path not in seen:
            results.append(path)
            seen.add(path)
    return results


def _plot_card(idx, default_var, slider_max, slider_value, time_point, use_range):
    id_prefix = "compare"
    options = default_var["options"]
    selected = default_var["value"] if options else None
    return html.Div(
        [
            html.Div(
                [
                    html.Label("Field"),
                    dcc.Dropdown(
                        id={"type": f"{id_prefix}-var-select", "index": idx},
                        options=[{"label": v, "value": v} for v in options],
                        value=selected,
                        clearable=False,
                    ),
                ],
                style={"width": "40%", "display": "inline-block"},
            ),
            html.Div(
                [
                    html.Div(id={"type": f"{id_prefix}-time-label", "index": idx}),
                    dcc.Checklist(
                        id={"type": f"{id_prefix}-time-mode", "index": idx},
                        options=[{"label": "Average over range", "value": "range"}],
                        value=["range"] if use_range else [],
                        style={"marginBottom": "6px"},
                    ),
                    # RangeSlider/Slider don't accept "style"; use wrapper div + callback for display toggling.
                    html.Div(
                        dcc.RangeSlider(
                            id={"type": f"{id_prefix}-time-range", "index": idx},
                            min=1,
                            max=slider_max,
                            value=slider_value,
                            step=1,
                            allowCross=False,
                            marks={1: "1", slider_max: str(slider_max)},
                            tooltip={"always_visible": True, "placement": "bottom"},
                            className="compare-time-range",
                        ),
                        id={"type": f"{id_prefix}-time-range-wrapper", "index": idx},
                        style={"display": "block" if use_range else "none"},
                    ),
                    html.Div(
                        dcc.Slider(
                            id={"type": f"{id_prefix}-time-point", "index": idx},
                            min=1,
                            max=slider_max,
                            value=time_point,
                            step=1,
                            marks={1: "1", slider_max: str(slider_max)},
                            tooltip={"always_visible": True, "placement": "bottom"},
                            included=False,
                            className="compare-time-point",
                        ),
                        id={"type": f"{id_prefix}-time-point-wrapper", "index": idx},
                        style={"display": "block" if not use_range else "none"},
                    ),
                ],
                style={"width": "48%", "display": "inline-block"},
            ),
            html.Div(
                [
                    dcc.Graph(id={"type": f"{id_prefix}-graph", "index": idx}),
                    html.Div(id={"type": f"{id_prefix}-error-text", "index": idx}),
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

    case_buttons = []
    for case in case_defs:
        if not case.get("name"):
            continue
        has_clubb = bool(_case_clubb_files(case))
        if has_clubb:
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
                case["name"],
                id={"type": "compare-case-button", "name": case["name"]},
                disabled=not has_clubb,
                style=style,
            )
        )

    default_var = {"value": None, "options": []}
    initial_plots = []
    layout = html.Div(
        [
            html.Div(
                [
                    html.Div("Cases:"),
                    html.Div(case_buttons),
                ],
                style={"marginBottom": "12px"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("LES file"),
                            dcc.Upload(
                                id="compare-upload-les",
                                children=html.Button("Choose LES file"),
                                multiple=False,
                            ),
                            html.Div(id="compare-file-les-name"),
                        ],
                        style={"width": "48%", "display": "inline-block"},
                    ),
                    html.Div(
                        [
                            html.Label("CLUBB file"),
                            dcc.Upload(
                                id="compare-upload-clubb",
                                children=html.Button("Choose CLUBB file"),
                                multiple=True,
                            ),
                            html.Div(id="compare-file-clubb-name"),
                        ],
                        style={"width": "48%", "display": "inline-block"},
                    ),
                ],
                style={"marginBottom": "10px"},
            ),
            html.Button("Add plot", id="compare-add-plot", n_clicks=0),
            dcc.Store(id="compare-plot-ids", data=initial_plots),
            dcc.Store(id="compare-file-les", data=None),
            dcc.Store(id="compare-file-clubb", data=None),
            dcc.Store(id="compare-case-defaults", data=None),
            dcc.Store(id="compare-selected-case", data=None),
            dcc.Store(id="compare-plot-state", data={}),
            dcc.Interval(id="compare-case-refresh", interval=3000, n_intervals=0),
            html.Div(
                id="compare-plot-container",
                style={
                    "display": "grid",
                    "gridTemplateColumns": "repeat(auto-fill, minmax(400px, 1fr))",
                    "gap": "16px",
                    "marginTop": "16px",
                },
            ),
        ],
        style={"padding": "10px"},
    )

    @app.callback(
        Output("compare-file-les", "data"),
        Output("compare-file-clubb", "data"),
        Output("compare-file-les-name", "children"),
        Output("compare-file-clubb-name", "children"),
        Output("compare-case-defaults", "data"),
        Output("compare-selected-case", "data"),
        Input("compare-upload-les", "contents"),
        Input("compare-upload-clubb", "contents"),
        Input({"type": "compare-case-button", "name": ALL}, "n_clicks"),
        State("compare-upload-les", "filename"),
        State("compare-upload-clubb", "filename"),
        State("compare-file-les", "data"),
        State("compare-file-clubb", "data"),
        prevent_initial_call=True,
    )
    def _select_case(
        les_contents,
        clubb_contents,
        _clicks,
        les_filename,
        clubb_filename,
        current_les,
        current_clubb,
    ):
        if not callback_context.triggered:
            return no_update, no_update, no_update, no_update, no_update
        button_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        if button_id == "compare-upload-les":
            path = _save_upload(les_contents, les_filename, "les")
            if path is None:
                return no_update, no_update, no_update, no_update, no_update, no_update
            label = _format_selected_label(path)
            return path, current_clubb, label, no_update, None, None
        if button_id == "compare-upload-clubb":
            if isinstance(clubb_contents, str):
                clubb_contents = [clubb_contents]
                clubb_filename = [clubb_filename]
            paths = _save_uploads(clubb_contents, clubb_filename, "clubb")
            if not paths:
                return no_update, no_update, no_update, no_update, no_update, no_update
            label = _format_selected_label(paths)
            return current_les, paths, no_update, label, None, None
        try:
            case_name = json.loads(button_id).get("name")
        except Exception:
            return no_update, no_update, no_update, no_update, no_update, no_update
        current_defs = load_case_definitions()
        current_by_name = {case["name"]: case for case in current_defs if case.get("name")}
        case = current_by_name.get(case_name)
        if not case:
            return no_update, no_update, no_update, no_update, no_update, no_update
        clubb_files = _case_clubb_files(case)
        if not clubb_files:
            return no_update, no_update, no_update, no_update, no_update, no_update
        les_file = case.get("les_file")
        les_label = _format_selected_label(les_file) if les_file else ""
        clubb_label = _format_selected_label(clubb_files)
        defaults = {
            "case_name": case_name,
            "start_time": case.get("start_time"),
            "end_time": case.get("end_time"),
            "les_file": les_file,
            "clubb_files": clubb_files,
        }
        return les_file, clubb_files, les_label, clubb_label, defaults, case_name

    @app.callback(
        Output({"type": "compare-case-button", "name": ALL}, "style"),
        Output({"type": "compare-case-button", "name": ALL}, "disabled"),
        Input("compare-case-refresh", "n_intervals"),
        Input("compare-selected-case", "data"),
        State({"type": "compare-case-button", "name": ALL}, "id"),
    )
    def _refresh_case_buttons(_tick, selected_case, ids):
        if not ids:
            return [], []
        case_defs = load_case_definitions()
        availability = {
            case["name"]: bool(_case_clubb_files(case))
            for case in case_defs
            if case.get("name")
        }
        styles = []
        disabled = []
        for case_id in ids:
            name = case_id.get("name")
            available = availability.get(name, False)
            is_selected = selected_case == name
            if available:
                style = {
                    "backgroundColor": "#2563eb",
                    "color": "#ffffff",
                    "border": "2px solid #f59e0b" if is_selected else "1px solid transparent",
                    "padding": "6px 10px",
                    "margin": "4px",
                    "borderRadius": "4px",
                    "cursor": "pointer",
                }
            else:
                style = {
                    "backgroundColor": "#c9c9c9",
                    "color": "#5f5f5f",
                    "border": "2px solid #f59e0b" if is_selected else "1px solid transparent",
                    "padding": "6px 10px",
                    "margin": "4px",
                    "borderRadius": "4px",
                    "cursor": "not-allowed",
                }
            styles.append(style)
            disabled.append(not available)
        return styles, disabled

    @app.callback(
        Output("compare-plot-ids", "data"),
        Input("compare-add-plot", "n_clicks"),
        Input("compare-file-les", "data"),
        Input("compare-file-clubb", "data"),
        Input("compare-case-defaults", "data"),
        State("compare-plot-ids", "data"),
        prevent_initial_call=True,
    )
    def _update_plot_ids(n_clicks, file_a, file_b, case_defaults, plot_ids):
        trigger = callback_context.triggered[0]["prop_id"].split(".")[0]
        if case_defaults and file_b and not (plot_ids or []):
            return [0]
        if trigger == "compare-case-defaults":
            if case_defaults and file_b:
                return [0]
            return plot_ids or []
        if trigger == "compare-file-clubb":
            if file_b:
                return plot_ids or []
            return []
        if trigger == "compare-file-les":
            return plot_ids or []
        if trigger == "compare-add-plot":
            if not file_b:
                return plot_ids or []
            if plot_ids is None:
                plot_ids = []
            next_id = max(plot_ids) + 1 if plot_ids else 0
            return plot_ids + [next_id]
        return plot_ids

    @app.callback(
        Output("compare-plot-container", "children"),
        Input("compare-plot-ids", "data"),
        Input("compare-file-les", "data"),
        Input("compare-file-clubb", "data"),
        Input("compare-case-defaults", "data"),
        State("compare-plot-state", "data"),
    )
    def _render_plots(plot_ids, file_a, file_b, case_defaults, plot_state):
        if not file_b:
            return [html.Div("Select CLUBB files to start.")]
        clubb_paths = file_b if isinstance(file_b, list) else [file_b]
        info_a = DatasetInfo(file_a) if file_a else None
        info_b = DatasetCollection(clubb_paths)
        try:
            clubb_vars = info_b.vars
            if not clubb_vars:
                return [html.Div("No variables with vertical levels found in CLUBB files.")]
            slider_b_max = max(info_b.time_len, 1)
            slider_a_max = max(info_a.time_len, 1) if info_a else 1
            slider_max = max(slider_a_max, slider_b_max)
            slider_value = [1, slider_max]
            if case_defaults and case_defaults.get("clubb_files") == clubb_paths:
                time_source = info_b.primary_time_source() or info_a
                slider_value = _minutes_to_slider_range(
                    time_source,
                    case_defaults.get("start_time"),
                    case_defaults.get("end_time"),
                    slider_max,
                )
            plot_state = plot_state or {}
            cards = []
            for idx in (plot_ids or []):
                state = plot_state.get(str(idx), {})
                selected = state.get("var")
                if selected not in clubb_vars:
                    selected = "thlm" if "thlm" in clubb_vars else clubb_vars[0]
                mode = state.get("time_mode", "range")
                time_range = state.get("time_range")
                if time_range:
                    t0, t1 = time_range
                    t0 = max(1, min(int(t0), slider_max))
                    t1 = max(1, min(int(t1), slider_max))
                    if t1 < t0:
                        t0, t1 = t1, t0
                    time_range = [t0, t1]
                else:
                    time_range = slider_value
                time_point = state.get("time_point")
                if time_point is not None:
                    time_point = max(1, min(int(time_point), slider_max))
                else:
                    time_point = time_range[0]
                local_default = {"value": selected, "options": clubb_vars}
                cards.append(
                    _plot_card(
                        idx,
                        local_default,
                        slider_max,
                        time_range,
                        time_point,
                        mode != "point",
                    )
                )
            return cards
        finally:
            if info_a is not None:
                info_a.close()
            info_b.close()

    @app.callback(
        Output("compare-plot-state", "data", allow_duplicate=True),
        Input({"type": "compare-var-select", "index": ALL}, "value"),
        Input({"type": "compare-time-range", "index": ALL}, "value"),
        Input({"type": "compare-time-point", "index": ALL}, "value"),
        Input({"type": "compare-time-mode", "index": ALL}, "value"),
        State({"type": "compare-var-select", "index": ALL}, "id"),
        State({"type": "compare-time-range", "index": ALL}, "id"),
        State({"type": "compare-time-point", "index": ALL}, "id"),
        State({"type": "compare-time-mode", "index": ALL}, "id"),
        prevent_initial_call=True,
    )
    def _sync_plot_state(
        var_values,
        time_ranges,
        time_points,
        time_modes,
        var_ids,
        time_ids,
        point_ids,
        mode_ids,
    ):
        plot_state = {}
        for meta, value in zip(var_ids or [], var_values or []):
            idx = meta.get("index")
            entry = dict(plot_state.get(str(idx), {}))
            if value is not None:
                entry["var"] = value
            plot_state[str(idx)] = entry
        for meta, value in zip(time_ids or [], time_ranges or []):
            idx = meta.get("index")
            entry = dict(plot_state.get(str(idx), {}))
            if value is not None:
                entry["time_range"] = value
            plot_state[str(idx)] = entry
        for meta, value in zip(point_ids or [], time_points or []):
            idx = meta.get("index")
            entry = dict(plot_state.get(str(idx), {}))
            if value is not None:
                entry["time_point"] = value
            plot_state[str(idx)] = entry
        for meta, value in zip(mode_ids or [], time_modes or []):
            idx = meta.get("index")
            entry = dict(plot_state.get(str(idx), {}))
            if value is not None and "range" in value:
                entry["time_mode"] = "range"
            else:
                entry["time_mode"] = "point"
            plot_state[str(idx)] = entry
        return plot_state

    @app.callback(
        Output("compare-plot-state", "data", allow_duplicate=True),
        Input("compare-file-clubb", "data"),
        prevent_initial_call=True,
    )
    def _reset_plot_state(_clubb_file):
        return {}

    @app.callback(
        Output({"type": "compare-time-range-wrapper", "index": MATCH}, "style"),
        Output({"type": "compare-time-point-wrapper", "index": MATCH}, "style"),
        Input({"type": "compare-time-mode", "index": MATCH}, "value"),
    )
    def _toggle_time_mode(t_mode):
        use_range = bool(t_mode) and "range" in t_mode
        if use_range:
            return {"display": "block"}, {"display": "none"}
        return {"display": "none"}, {"display": "block"}

    @app.callback(
        Output({"type": "compare-graph", "index": MATCH}, "figure"),
        Output({"type": "compare-error-text", "index": MATCH}, "children"),
        Output({"type": "compare-time-label", "index": MATCH}, "children"),
        Input({"type": "compare-var-select", "index": MATCH}, "value"),
        Input({"type": "compare-time-range", "index": MATCH}, "value"),
        Input({"type": "compare-time-point", "index": MATCH}, "value"),
        Input({"type": "compare-time-mode", "index": MATCH}, "value"),
        Input("compare-file-les", "data"),
        Input("compare-file-clubb", "data"),
    )
    def _update_graphs(var_name, t_range, t_point, t_mode, file_a, file_b):
        if not file_b or var_name is None:
            return go.Figure(), "", ""
        clubb_paths = file_b if isinstance(file_b, list) else [file_b]
        info_a = DatasetInfo(file_a) if file_a else None
        info_b = DatasetCollection(clubb_paths)
        try:
            use_range = bool(t_mode) and "range" in t_mode
            if use_range:
                t_range_zero = [max(0, t_range[0] - 1), max(0, t_range[1] - 1)]
            else:
                t_idx = max(0, (t_point or 1) - 1)
                t_range_zero = [t_idx, t_idx]
            clubb_ds = info_b.dataset_for_var(var_name) or info_b.primary_time_source()
            if clubb_ds is None:
                return go.Figure(), "Variable not found in CLUBB files.", ""
            z_clubb, clubb_profile, clubb_profile2 = _time_avg_profile(
                clubb_ds, var_name, t_range_zero
            )
            z_les = None
            les_profile = None
            if info_a is not None and var_name in info_a.vars:
                z_les, les_profile, _ = _time_avg_profile(info_a, var_name, t_range_zero)

            z_units = getattr(
                clubb_ds.ds.variables[clubb_ds.var_info[var_name]["z_dim"]], "units", ""
            )
            fig = _make_figure(
                z_les,
                les_profile,
                z_clubb,
                clubb_profile,
                f"{var_name}",
                var_name,
                z_units,
            )
            fig.update_layout(width=400, height=400)
            error_text = ""
            if les_profile is not None and z_les is not None:
                error_val = _compute_error(
                    les_profile,
                    z_les,
                    clubb_profile,
                    clubb_profile2,
                    z_clubb,
                    var_name,
                )
                error_text = f"Error: {error_val}"
            if use_range:
                time_label = _time_label(clubb_ds, t_range_zero)
            else:
                time_label = _time_label_point(clubb_ds, t_range_zero[0])
            return fig, error_text, time_label
        finally:
            if info_a is not None:
                info_a.close()
            info_b.close()

    return dcc.Tab(label="Plots", value="plots", children=layout)

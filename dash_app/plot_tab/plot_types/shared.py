import glob
import os
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
from dash import dcc, html
from netCDF4 import Dataset, chartostring

from ..case_definitions import load_case_definitions
from .budget_groups import BUDGET_GROUPS

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUTPUT_DIR = os.path.join(REPO_ROOT, "output")
OUTPUT_FILE_SUFFIXES = ["_stats.nc"]

Z_DIM_NAMES = {"z", "zm", "zt", "lh_zt", "altitude", "height", "lev"}
T_DIM_NAMES = {"t", "time"}
X_DIM_NAMES = {"x", "longitude", "lon"}
Y_DIM_NAMES = {"y", "latitude", "lat"}
COL_DIM_NAMES = {"columns", "column", "col", "ngrdcol"}
SUBCOLUMN_DIM_NAMES = {
    "subcolumn",
    "subcolumns",
    "sample_point",
    "sample_points",
    "lh_sample_point",
    "silhs_sample_point",
}

DEFAULT_PROFILE_VARS = ["thlm", "rtm", "wp2"]
TIMEHEIGHT_DEFAULT_VARS = ["corr_rt_thl_1", "corr_rt_thl_2", "corr_w_chi_1"]
SUBCOLUMN_DEFAULT_VARS = ["w", "chi", "eta"]
DEFAULT_TIMESERIES_VARS = ["lwp", "T_in_K", "lh_rrm_auto", "lh_rrm_mc_nonadj"]

TIMEHEIGHT_CATALOG = [
    {"name": "corr_chi_eta_1", "title": "Correlation of chi and eta, 1st component", "axis_title": "corr_chi_eta_1 [$-$]"},
    {"name": "corr_chi_eta_2", "title": "Correlation of chi and eta, 2nd component", "axis_title": "corr_chi_eta_2 [$-$]"},
    {"name": "corr_chi_eta_1_ca", "title": "Correlation of chi and eta (ca), 1st component", "axis_title": "corr_chi_eta_1_ca [$-$]"},
    {"name": "corr_chi_eta_2_ca", "title": "Correlation of chi and eta (ca), 2nd component", "axis_title": "corr_chi_eta_2_ca [$-$]"},
    {"name": "corr_chi_Ncn_1", "title": "Correlation (in-precip) of w and Ncn, 1st component", "axis_title": "corr_chi_Ncn_1 [$-$]"},
    {"name": "corr_chi_Ncn_2", "title": "Correlation (in-precip) of w and Ncn, 2nd component", "axis_title": "corr_chi_Ncn_2 [$-$]"},
    {"name": "corr_eta_Ncn_1", "title": "Correlation of eta and N_cn, 1st component", "axis_title": "corr_eta_Ncn_1 [$-$]"},
    {"name": "corr_eta_Ncn_2", "title": "Correlation of eta and N_cn, 2nd component", "axis_title": "corr_eta_Ncn_2 [$-$]"},
    {"name": "corr_rt_thl_1", "title": "Correlation of rt and thl, 1st component", "axis_title": "corr_rt_thl_1 [$-$]"},
    {"name": "corr_rt_thl_2", "title": "Correlation of rt and thl, 2nd component", "axis_title": "corr_rt_thl_2 [$-$]"},
    {"name": "corr_w_chi_1", "title": "Correlation of w and chi, 1st component", "axis_title": "corr_w_chi_1 [$-$]"},
    {"name": "corr_w_chi_2", "title": "Correlation of w and chi, 2nd component", "axis_title": "corr_w_chi_2 [$-$]"},
    {"name": "corr_w_chi_1_ca", "title": "Correlation of w and chi (ca), 1st component", "axis_title": "corr_w_chi_1_ca [$-$]"},
    {"name": "corr_w_chi_2_ca", "title": "Correlation of w and chi (ca), 2nd component", "axis_title": "corr_w_chi_2_ca [$-$]"},
    {"name": "corr_w_Nr_1", "title": "Correlation (in-precip) of w and Nr, 1st component", "axis_title": "corr_w_Nr_1 [$-$]"},
    {"name": "corr_w_Nr_2", "title": "Correlation (in-precip) of w and Nr, 2nd component", "axis_title": "corr_w_Nr_2 [$-$]"},
    {"name": "corr_w_Ncn_1", "title": "Correlation (in-precip) of w and Ncn, 1st component", "axis_title": "corr_w_Ncn_1 [$-$]"},
    {"name": "corr_w_Ncn_2", "title": "Correlation (in-precip) of w and Ncn, 2nd component", "axis_title": "corr_w_Ncn_2 [$-$]"},
    {"name": "corr_w_eta_1", "title": "Correlation (in-precip) of w and eta, 1st component", "axis_title": "corr_w_eta_1 [$-$]"},
    {"name": "corr_w_eta_2", "title": "Correlation (in-precip) of w and eta, 2nd component", "axis_title": "corr_w_eta_2 [$-$]"},
    {"name": "corr_w_eta_1_ca", "title": "Correlation (in-precip) of w and eta (ca), 1st component", "axis_title": "corr_w_eta_1_ca [$-$]"},
    {"name": "corr_w_eta_2_ca", "title": "Correlation (in-precip) of w and eta (ca), 2nd component", "axis_title": "corr_w_eta_2_ca [$-$]"},
    {"name": "corr_w_rt_1", "title": "Correlation (in-precip) of w and rt, 1st component", "axis_title": "corr_w_rt_1 [$-$]"},
    {"name": "corr_w_rt_2", "title": "Correlation (in-precip) of w and rt, 2nd component", "axis_title": "corr_w_rt_2 [$-$]"},
    {"name": "corr_w_thl_1", "title": "Correlation (in-precip) of w and thl, 1st component", "axis_title": "corr_w_thl_1 [$-$]"},
    {"name": "corr_w_thl_2", "title": "Correlation (in-precip) of w and thl, 2nd component", "axis_title": "corr_w_thl_2 [$-$]"},
    {"name": "corr_w_rr_1", "title": "Correlation (in-precip) of w and rr, 1st component", "axis_title": "corr_w_rr_1 [-]"},
    {"name": "corr_w_rr_2", "title": "Correlation (in-precip) of w and rr, 2nd component", "axis_title": "corr_w_rr_2 [-]"},
    {"name": "corr_chi_rr_1", "title": "Correlation (in-precip) of chi and rr, 1st component", "axis_title": "corr_chi_rr_1 [-]"},
    {"name": "corr_chi_rr_2", "title": "Correlation (in-precip) of chi and rr, 2nd component", "axis_title": "corr_chi_rr_2 [-]"},
    {"name": "corr_chi_Nr_1", "title": "Correlation (in-precip) of chi and Nr, 1st component", "axis_title": "corr_chi_Nr_1 [-]"},
    {"name": "corr_chi_Nr_2", "title": "Correlation (in-precip) of chi and Nr, 2nd component", "axis_title": "corr_chi_Nr_2 [-]"},
    {"name": "corr_rr_Nr_1", "title": "Correlation (in-precip) of rr and Nr, 1st component", "axis_title": "corr_rr_Nr_1 [-]"},
    {"name": "corr_rr_Nr_2", "title": "Correlation (in-precip) of rr and Nr, 2nd component", "axis_title": "corr_rr_Nr_2 [-]"},
]

SUBCOLUMN_CATALOG = [
    {"name": "w", "label": "w"},
    {"name": "chi", "label": "chi"},
    {"name": "eta", "label": "eta"},
    {"name": "Ncn", "label": "Ncn"},
    {"name": "rr", "label": "rr"},
    {"name": "Nr", "label": "Nr"},
    {"name": "ri", "label": "ri"},
    {"name": "Ni", "label": "Ni"},
    {"name": "rg", "label": "rg"},
    {"name": "Ng", "label": "Ng"},
    {"name": "rs", "label": "rs"},
    {"name": "Ns", "label": "Ns"},
]

_CACHE_MAX_ENTRIES = 256
_EXTRACT_CACHE = OrderedDict()
_PRELOAD_CACHE = OrderedDict()
_PLOT_DATA_CACHE = OrderedDict()
_DATASET_META_CACHE = OrderedDict()
_COLLECTION_META_CACHE = OrderedDict()

PLOT_SIZE_HEIGHTS = {
    "normal": 450,
    "large": 700,
}

COLUMN_LINE_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]

SOURCE_LINE_DASHES = [
    "solid",
    "dash",
    "dot",
    "dashdot",
    "longdash",
    "longdashdot",
]


def _file_signature(path):
    try:
        stat = os.stat(path)
    except OSError:
        return (path, None, None)
    return (path, stat.st_mtime_ns, stat.st_size)


def _freeze_cached_value(value):
    if isinstance(value, np.ndarray):
        frozen = np.asarray(value)
        frozen.setflags(write=False)
        return frozen
    if isinstance(value, dict):
        return {key: _freeze_cached_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return tuple(_freeze_cached_value(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze_cached_value(item) for item in value)
    return value


def _cached_extract(key, builder):
    cached = _EXTRACT_CACHE.get(key)
    if cached is not None:
        _EXTRACT_CACHE.move_to_end(key)
        return cached
    value = _freeze_cached_value(builder())
    _EXTRACT_CACHE[key] = value
    _EXTRACT_CACHE.move_to_end(key)
    while len(_EXTRACT_CACHE) > _CACHE_MAX_ENTRIES:
        _EXTRACT_CACHE.popitem(last=False)
    return value


def _cached_preload(key, builder):
    cached = _PRELOAD_CACHE.get(key)
    if cached is not None:
        _PRELOAD_CACHE.move_to_end(key)
        return cached
    value = _freeze_cached_value(builder())
    _PRELOAD_CACHE[key] = value
    _PRELOAD_CACHE.move_to_end(key)
    return value


def _cached_plot_data(key, builder):
    cached = _PLOT_DATA_CACHE.get(key)
    if cached is not None:
        _PLOT_DATA_CACHE.move_to_end(key)
        return cached
    value = _freeze_cached_value(builder())
    _PLOT_DATA_CACHE[key] = value
    _PLOT_DATA_CACHE.move_to_end(key)
    return value


def _cached_lru(cache, max_entries, key, builder):
    cached = cache.get(key)
    if cached is not None:
        cache.move_to_end(key)
        return cached
    value = builder()
    cache[key] = value
    cache.move_to_end(key)
    while len(cache) > max_entries:
        cache.popitem(last=False)
    return value


def clear_all_caches():
    _PLOT_DATA_CACHE.clear()
    _PRELOAD_CACHE.clear()
    _EXTRACT_CACHE.clear()
    _DATASET_META_CACHE.clear()
    _COLLECTION_META_CACHE.clear()


def _dataset_metadata(ds, path):
    key = ("dataset_meta", _file_signature(path))

    def _build():
        var_info = {}
        vars_found = []
        dim_units = {}
        for dim_name in ds.dimensions:
            if dim_name in ds.variables:
                dim_units[dim_name] = getattr(ds.variables[dim_name], "units", "") or ""
        for name, var in ds.variables.items():
            if name in ds.dimensions:
                continue
            dims = tuple(var.dimensions)
            var_info[name] = {
                "dims": dims,
                "z_dim": find_dim(dims, Z_DIM_NAMES),
                "t_dim": find_dim(dims, T_DIM_NAMES),
                "col_dim": find_dim(dims, COL_DIM_NAMES),
                "subcol_dim": find_subcolumn_dim(dims),
                "units": getattr(var, "units", ""),
                "long_name": getattr(var, "long_name", "") or "",
            }
            vars_found.append(name)
        vars_found = tuple(sorted(vars_found))
        profile_vars = tuple(name for name in vars_found if var_info[name]["z_dim"] is not None)
        dim_names = tuple(ds.dimensions.keys())
        time_dim = find_dim(dim_names, T_DIM_NAMES)
        time_len = len(ds.dimensions[time_dim]) if time_dim else 1
        time_values = None
        time_units = None
        if time_dim and time_dim in ds.variables:
            time_var = ds.variables[time_dim]
            time_values = _freeze_cached_value(as_array(time_var[:]))
            time_units = getattr(time_var, "units", None)
        columns_dim = find_dim(dim_names, COL_DIM_NAMES)
        columns_len = len(ds.dimensions[columns_dim]) if columns_dim else 1
        return {
            "vars": vars_found,
            "var_info": var_info,
            "profile_vars": profile_vars,
            "time_dim": time_dim,
            "time_len": time_len,
            "time_values": time_values,
            "time_units": time_units,
            "columns_dim": columns_dim,
            "columns_len": columns_len,
            "dim_units": dim_units,
        }

    return _cached_lru(_DATASET_META_CACHE, _CACHE_MAX_ENTRIES, key, _build)


def _finite_bounds(values):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return None
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None
    return float(np.min(finite)), float(np.max(finite))


def _collection_metadata(datasets):
    key = ("collection_meta", tuple(_file_signature(dataset.path) for dataset in datasets))

    def _build():
        var_to_index = {}
        for dataset_index, dataset in enumerate(datasets):
            for var_name in dataset.vars:
                var_to_index.setdefault(var_name, dataset_index)
        return {
            "vars": tuple(sorted(var_to_index.keys())),
            "var_to_index": var_to_index,
            "time_len": max((dataset.time_len for dataset in datasets), default=1),
            "columns_len": max((dataset.columns_len for dataset in datasets), default=1),
        }

    return _cached_lru(_COLLECTION_META_CACHE, _CACHE_MAX_ENTRIES, key, _build)


def column_line_color(index):
    """Return a stable line color for a column index."""
    return COLUMN_LINE_COLORS[index % len(COLUMN_LINE_COLORS)]


def source_line_dash(index):
    """Return a stable dash style for one compare-mode source index."""
    return SOURCE_LINE_DASHES[index % len(SOURCE_LINE_DASHES)]


def source_line_color(index):
    """Return a stable line color for one compare-mode source index."""
    return COLUMN_LINE_COLORS[index % len(COLUMN_LINE_COLORS)]


def add_directory_legend_traces(fig, labels):
    """Add dummy legend entries so compare/overplot legends show only directory styles."""
    for source_idx, label in enumerate(labels):
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines",
                name=label,
                line={"width": 2.0, "dash": source_line_dash(source_idx), "color": source_line_color(source_idx)},
                showlegend=True,
                hoverinfo="skip",
            )
        )


@dataclass
class DatasetInfo:
    path: str

    def __post_init__(self):
        self.ds = Dataset(self.path, "r")
        # Some CLUBB stats files mark 0 as a fill/missing value. For plotting, 0 is
        # a valid value and should remain visible rather than being auto-masked away.
        self.ds.set_auto_mask(False)
        meta = _dataset_metadata(self.ds, self.path)
        self.vars = meta["vars"]
        self.var_info = meta["var_info"]
        self.profile_vars = meta["profile_vars"]
        self.time_dim = meta["time_dim"]
        self.time_len = meta["time_len"]
        self.time_values = meta["time_values"]
        self.time_units = meta["time_units"]
        self.columns_dim = meta["columns_dim"]
        self.columns_len = meta["columns_len"]
        self.dim_units = meta["dim_units"]

    def close(self):
        self.ds.close()


class DatasetCollection:
    def __init__(self, paths):
        self.paths = list(paths or [])
        self.datasets = [DatasetInfo(path) for path in self.paths]
        meta = _collection_metadata(self.datasets)
        self.vars = meta["vars"]
        self.var_to_ds = {
            var_name: self.datasets[dataset_index]
            for var_name, dataset_index in meta["var_to_index"].items()
        }
        self.time_len = meta["time_len"]
        self.columns_len = meta["columns_len"]

    def dataset_for_var(self, var_name):
        return self.var_to_ds.get(var_name)

    def primary_time_source(self):
        if not self.datasets:
            return None
        return max(self.datasets, key=lambda dataset: dataset.time_len)

    def close(self):
        for dataset in self.datasets:
            dataset.close()


def dataset_metadata_for_path(path):
    key = ("dataset_meta", _file_signature(path))
    cached = _DATASET_META_CACHE.get(key)
    if cached is not None:
        _DATASET_META_CACHE.move_to_end(key)
        return cached
    ds = Dataset(path, "r")
    try:
        ds.set_auto_mask(False)
        return _dataset_metadata(ds, path)
    finally:
        ds.close()


def dataset_info_for_var(paths, var_name):
    for path in paths or []:
        meta = dataset_metadata_for_path(path)
        if var_name in meta["vars"]:
            return path, meta
    return None, None


def dataset_infos_for_var(paths, var_name):
    matches = []
    for path in paths or []:
        meta = dataset_metadata_for_path(path)
        if var_name in meta["vars"]:
            matches.append((path, meta))
    return matches


def find_dim(dim_names, candidates):
    for dim_name in dim_names:
        if dim_name.lower() in candidates:
            return dim_name
    return None


def find_subcolumn_dim(dim_names):
    for dim_name in dim_names:
        lowered = dim_name.lower()
        if lowered in SUBCOLUMN_DIM_NAMES:
            return dim_name
    return None


def as_array(data):
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    return np.asarray(data)


def mean_over_axis(data, axis):
    """Average over one axis without warning on zero-length slices."""
    arr = np.asarray(data)
    axis = int(axis)
    if arr.ndim == 0:
        return arr
    if axis < 0:
        axis += arr.ndim
    if axis < 0 or axis >= arr.ndim:
        raise np.AxisError(axis, arr.ndim)
    if arr.shape[axis] == 0:
        out_shape = arr.shape[:axis] + arr.shape[axis + 1 :]
        return np.full(out_shape, np.nan, dtype=float)
    return np.mean(arr, axis=axis)


def time_units_factor(units):
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


def time_values_seconds(ds_info):
    if ds_info is None:
        return None
    if ds_info.time_values is None:
        return None
    key = ("time_seconds", _file_signature(ds_info.path), ds_info.time_dim, ds_info.time_units)
    return _cached_extract(key, lambda: np.asarray(ds_info.time_values * time_units_factor(ds_info.time_units), dtype=float))


def time_values_pyplotgen_minutes(ds_info):
    if ds_info is None or ds_info.time_values is None:
        return None
    key = ("time_pyplotgen_minutes", _file_signature(ds_info.path), ds_info.time_dim, ds_info.time_units)
    def _build():
        time_values = np.asarray(ds_info.time_values, dtype=float)
        units = (ds_info.time_units or "").lower()
        if "day" in units:
            time_values = time_values * 24.0 * 60.0
        elif "hour" in units:
            time_values = time_values * 60.0
        elif "sec" in units:
            time_values = time_values / 60.0
        if len(time_values) <= 1:
            return np.array([1.0], dtype=float) if len(time_values) == 1 else np.array([], dtype=float)
        delta_t = float(time_values[1] - time_values[0])
        return np.array([delta_t * (idx + 1) for idx in range(len(time_values))], dtype=float)
    return _cached_extract(key, _build)


def slider_value_to_index(value, max_len):
    if max_len <= 0:
        return 0
    return max(0, min(int(value or 1) - 1, max_len - 1))


def slider_range_to_indices(value, max_len):
    if max_len <= 0:
        return [0, 0]
    if not value or len(value) != 2:
        return [0, max_len - 1]
    start = slider_value_to_index(value[0], max_len)
    end = slider_value_to_index(value[1], max_len)
    if end < start:
        start, end = end, start
    return [start, end]


def serialize_time_seconds(ds_info):
    """Return the primary time axis as a JSON-safe list of seconds."""
    time_sec = time_values_seconds(ds_info)
    if time_sec is None:
        return []
    return [float(value) for value in np.asarray(time_sec, dtype=float).tolist()]


def minutes_to_slider_range(ds_info, start_min, end_min, slider_max):
    if slider_max <= 0:
        return [1, 1]
    if start_min is None or end_min is None:
        return [1, slider_max]
    time_minutes = time_values_pyplotgen_minutes(ds_info)
    if time_minutes is None or len(time_minutes) == 0:
        start_idx = max(1, min(int(start_min), slider_max))
        end_idx = max(1, min(int(end_min), slider_max))
    else:
        start_idx = int(np.argmin(np.abs(time_minutes - float(start_min)))) + 1
        end_idx = int(np.argmin(np.abs(time_minutes - float(end_min)))) + 1
        start_idx = max(1, min(start_idx, slider_max))
        end_idx = max(1, min(end_idx, slider_max))
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    return [start_idx, end_idx]


def time_label_for_range(time_seconds, slider_range):
    if time_seconds is None:
        return "Time range"
    time_sec = np.asarray(time_seconds, dtype=float)
    if len(time_sec) == 0:
        start_idx, end_idx = slider_range_to_indices(slider_range, 0)
        return f"Time: {start_idx + 1} - {end_idx + 1}"
    start_idx, end_idx = slider_range_to_indices(slider_range, len(time_sec))
    start_sec = float(time_sec[start_idx])
    end_sec = float(time_sec[end_idx])
    start_min = start_sec / 60.0
    end_min = end_sec / 60.0
    return f"Time: {start_sec:g}-{end_sec:g}s ({start_min:g}-{end_min:g}m)"


def time_label_for_point(time_seconds, slider_value):
    if time_seconds is None:
        return "Time point"
    time_sec = np.asarray(time_seconds, dtype=float)
    if len(time_sec) == 0:
        idx = slider_value_to_index(slider_value, 1)
        return f"Time: {idx + 1}"
    idx = slider_value_to_index(slider_value, len(time_sec))
    time_sec_value = float(time_sec[idx])
    time_min_value = time_sec_value / 60.0
    return f"Time: {time_sec_value:g}s ({time_min_value:g}m)"


def plot_theme_colors(theme_name):
    theme = (theme_name or "dark").lower()
    if theme == "light":
        return {
            "paper_bg": "#ffffff",
            "plot_bg": "#ffffff",
            "font": "#0f172a",
            "grid": "#cbd5e1",
            "line": "#94a3b8",
        }
    return {
        "paper_bg": "#0f172a",
        "plot_bg": "#0f172a",
        "font": "#e5e7eb",
        "grid": "#334155",
        "line": "#64748b",
    }


def apply_plot_theme(fig, theme_name):
    colors = plot_theme_colors(theme_name)
    fig.update_layout(
        paper_bgcolor=colors["paper_bg"],
        plot_bgcolor=colors["plot_bg"],
        font={"color": colors["font"]},
    )
    fig.update_xaxes(
        gridcolor=colors["grid"],
        linecolor=colors["line"],
        zerolinecolor=colors["grid"],
        color=colors["font"],
    )
    fig.update_yaxes(
        gridcolor=colors["grid"],
        linecolor=colors["line"],
        zerolinecolor=colors["grid"],
        color=colors["font"],
    )
    return fig


def make_empty_figure(message, theme_name="dark"):
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"color": plot_theme_colors(theme_name)["font"]},
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(margin=dict(l=40, r=20, t=40, b=40), height=430)
    return apply_plot_theme(fig, theme_name)


def variable_long_name(ds_info, var_name):
    if ds_info is None:
        return ""
    info = ds_info.var_info.get(var_name)
    if info is None:
        return ""
    return info.get("long_name", "") or ""


def variable_units(ds_info, var_name):
    if ds_info is None:
        return ""
    info = ds_info.var_info.get(var_name)
    if info is None:
        return ""
    return info.get("units", "") or ""


def dimension_units(ds_info, dim_name):
    if ds_info is None or not dim_name:
        return ""
    return (getattr(ds_info, "dim_units", {}) or {}).get(dim_name, "")


def format_plot_title(name, long_name=""):
    if long_name:
        return f"{long_name} ({name})"
    return name


def format_axis_title(name, units=""):
    return f"{name} [{units}]" if units else name


def format_height_axis_title(units=""):
    return format_axis_title("Height", units)


def figure_uirevision(plot_type, case_data, var_name, *extra_parts):
    case_name = (case_data or {}).get("name") or ""
    file_key = "|".join((case_data or {}).get("files") or [])
    parts = [plot_type, case_name, file_key, var_name or ""]
    parts.extend("" if part is None else str(part) for part in extra_parts)
    return "::".join(parts)


def apply_figure_chrome(fig, title, showlegend=False, height=430, uirevision=None):
    top_margin = 56
    bottom_margin = 88 if showlegend else 52
    fig.update_layout(
        title={"text": title, "x": 0.5, "xanchor": "center"},
        margin=dict(l=60, r=20, t=top_margin, b=bottom_margin),
        showlegend=showlegend,
        legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="center", x=0.5),
        height=height,
        uirevision=uirevision,
    )
    return fig


def padded_range(min_value, max_value, fraction=0.03):
    min_value = float(min_value)
    max_value = float(max_value)
    if max_value < min_value:
        min_value, max_value = max_value, min_value
    span = max_value - min_value
    if span == 0.0:
        pad = max(abs(min_value) * fraction, 1.0)
    else:
        pad = span * fraction
    return [min_value - pad, max_value + pad]


def padded_data_range(values, fraction=0.03):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return None
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None
    return padded_range(np.min(finite), np.max(finite), fraction=fraction)


def apply_axis_bounds(fig, axis_name, bounds):
    if not bounds or len(bounds) != 2:
        return fig
    low, high = float(bounds[0]), float(bounds[1])
    if high < low:
        low, high = high, low
    axis_update = {"range": [low, high], "minallowed": low, "maxallowed": high}
    if axis_name == "x":
        fig.update_xaxes(**axis_update)
    else:
        fig.update_yaxes(**axis_update)
    return fig


def normalize_plot_size(size_value):
    if size_value in PLOT_SIZE_HEIGHTS:
        return size_value
    return "normal"


def figure_height_for_size(size_value):
    return PLOT_SIZE_HEIGHTS[normalize_plot_size(size_value)]


def plot_card_class_name(size_value):
    size_name = normalize_plot_size(size_value)
    return f"plots-card plots-card-size-{size_name}"


def plot_size_button_props(size_value):
    size_name = normalize_plot_size(size_value)
    if size_name == "large":
        return "-", "plots-card-size-toggle plots-card-size-toggle-large"
    return "+", "plots-card-size-toggle plots-card-size-toggle-normal"


def plot_graph_shell_style(size_value):
    return {"height": f"{figure_height_for_size(size_value)}px"}


def make_plot_card(controls, graph_id, size_button=None, size_value="normal", subtitle=None, card_id=None, close_button_id=None, render_signal_id=None, graph_shell_id=None, size_store_id=None):
    title_children = []
    if subtitle:
        title_children.append(
            html.Span(
                "i",
                className="plots-card-info",
                title=subtitle,
            )
        )
    header = html.Div(
        [
            html.Div(controls, style={"flex": "1", "minWidth": 0}),
            html.Div(
                [
                    *title_children,
                    size_button if size_button is not None else html.Div(),
                ],
                style={"display": "flex", "alignItems": "center", "gap": "8px"},
            ),
            html.Button("x", id=close_button_id, className="plots-card-close", title="Close plot") if close_button_id is not None else html.Div(),
        ],
        style={"display": "flex", "justifyContent": "space-between", "alignItems": "center", "gap": "10px", "marginBottom": "10px"},
    )
    return html.Div(
        [
            header,
            dcc.Store(id=size_store_id, data=normalize_plot_size(size_value)) if size_store_id is not None else html.Div(style={"display": "none"}),
            html.Div(id=render_signal_id, style={"display": "none"}) if render_signal_id is not None else html.Div(style={"display": "none"}),
            html.Div(
                dcc.Graph(id=graph_id, config={"displayModeBar": True}, className="plots-graph", style={"height": "100%"}),
                id=graph_shell_id,
                className="plots-graph-shell",
                style=plot_graph_shell_style(size_value),
            ),
        ],
        id=card_id,
        key=str(card_id),
        className=plot_card_class_name(size_value),
        style={"padding": "10px"},
    )


def normalize_output_directory(directory):
    expanded = os.path.expanduser(directory)
    if os.path.isabs(expanded):
        return os.path.abspath(expanded)
    cwd_path = os.path.abspath(expanded)
    repo_path = os.path.abspath(os.path.join(REPO_ROOT, expanded))
    if os.path.isdir(cwd_path):
        return cwd_path
    if os.path.isdir(repo_path):
        return repo_path
    return repo_path


def source_labels_for_dirs(directories):
    normalized = [normalize_output_directory(directory) for directory in directories or []]
    basenames = [os.path.basename(path.rstrip(os.sep)) or path for path in normalized]
    labels = []
    for path, basename in zip(normalized, basenames):
        if basenames.count(basename) == 1:
            labels.append(basename)
        else:
            labels.append(path)
    return labels


def _scan_cases_in_directory(directory):
    cases = {}
    directory = normalize_output_directory(directory)
    sorted_suffixes = sorted(OUTPUT_FILE_SUFFIXES, key=len, reverse=True)
    for path in glob.glob(os.path.join(directory, "*.nc")):
        base = os.path.basename(path)
        case_name = None
        for suffix in sorted_suffixes:
            if base.endswith(suffix):
                case_name = base[: -len(suffix)]
                break
        if case_name:
            cases[case_name] = path
    return cases


def scan_output_cases(directories=None):
    if directories is None:
        raw_dirs = [OUTPUT_DIR]
    else:
        raw_dirs = list(directories)
    if not raw_dirs:
        return {}
    normalized_dirs = [normalize_output_directory(directory) for directory in raw_dirs]
    dir_case_maps = [_scan_cases_in_directory(directory) for directory in normalized_dirs]
    if not dir_case_maps:
        return {}
    if len(dir_case_maps) == 1:
        return {case_name: [path] for case_name, path in dir_case_maps[0].items()}
    common_cases = set(dir_case_maps[0].keys())
    for case_map in dir_case_maps[1:]:
        common_cases &= set(case_map.keys())
    return {
        case_name: [case_map[case_name] for case_map in dir_case_maps]
        for case_name in sorted(common_cases)
    }


def ordered_case_names(available_case_names=None):
    case_defs = load_case_definitions()
    if available_case_names is None:
        available = set(scan_output_cases().keys())
    else:
        available = set(available_case_names)
    ordered = []
    seen = set()
    for case in case_defs:
        name = case.get("name")
        if name and name in available and name not in seen:
            ordered.append(name)
            seen.add(name)
    for name in sorted(available):
        if name not in seen:
            ordered.append(name)
            seen.add(name)
    return ordered


def case_definition_by_name():
    return {case["name"]: case for case in load_case_definitions() if case.get("name")}


def list_profile_vars(collection, require_all=False):
    if not require_all:
        return [{"label": name, "value": name} for name in collection.vars if has_profile_shape(collection, name)]
    common_vars = None
    for dataset in collection.datasets:
        dataset_vars = {name for name in dataset.vars if dataset.var_info[name]["z_dim"] is not None}
        common_vars = dataset_vars if common_vars is None else common_vars & dataset_vars
    common_vars = sorted(common_vars or [])
    return [{"label": name, "value": name} for name in common_vars]


def list_budget_groups(collection):
    options = []
    for name, group in BUDGET_GROUPS.items():
        if any(collection.dataset_for_var(term) is not None for term in group["terms"]):
            options.append({"label": group["label"], "value": name})
    return options


def build_case_data(case_name, files, directories=None):
    case_def = case_definition_by_name().get(case_name, {})
    from ..benchmark_overlay import build_case_benchmark_info

    collection = DatasetCollection(files)
    try:
        compare_mode = len(files) > 1
        source_dirs = [normalize_output_directory(directory) for directory in (directories or [os.path.dirname(path) for path in files])]
        profile_vars = list_profile_vars(collection, require_all=compare_mode)
        budget_groups = [] if compare_mode else list_budget_groups(collection)
        timeheight_vars = [] if compare_mode else list_timeheight_vars(collection)
        subcolumn_vars = [] if compare_mode else list_subcolumn_vars(collection)
        timeseries_vars = list_timeseries_vars(collection, require_all=compare_mode)
        slider_max = max(collection.time_len, 1)
        data_height_min, data_height_max = collection_height_bounds(collection)
        case_height_min = case_def.get("height_min_value")
        case_height_max = case_def.get("height_max_value")
        height_slider_min = data_height_min
        height_slider_max = data_height_max
        if case_height_min is not None:
            height_slider_min = case_height_min if height_slider_min is None else min(height_slider_min, case_height_min)
        if case_height_max is not None:
            height_slider_max = case_height_max if height_slider_max is None else max(height_slider_max, case_height_max)
        if height_slider_min is None or height_slider_max is None:
            height_slider_min, height_slider_max = 0.0, 1.0
        if height_slider_max < height_slider_min:
            height_slider_min, height_slider_max = height_slider_max, height_slider_min
        default_height_range = [
            float(case_height_min if case_height_min is not None else height_slider_min),
            float(case_height_max if case_height_max is not None else height_slider_max),
        ]
        if default_height_range[1] < default_height_range[0]:
            default_height_range.reverse()
        height_step = max((float(height_slider_max) - float(height_slider_min)) / 500.0, 1.0)
        time_source = collection.primary_time_source()
        default_range = minutes_to_slider_range(
            time_source,
            case_def.get("start_time"),
            case_def.get("end_time"),
            slider_max,
        )
        time_seconds = serialize_time_seconds(time_source)
        benchmarks = build_case_benchmark_info(case_def)
        return {
            "name": case_name,
            "files": list(files),
            "compare_mode": compare_mode,
            "output_dirs": source_dirs,
            "source_labels": source_labels_for_dirs(source_dirs),
            "start_time": case_def.get("start_time"),
            "end_time": case_def.get("end_time"),
            "height_min_value": case_height_min,
            "height_max_value": case_height_max,
            "height_slider_min": float(height_slider_min),
            "height_slider_max": float(height_slider_max),
            "height_step": float(height_step),
            "default_height_range": default_height_range,
            "columns_len": max(collection.columns_len, 1),
            "time_len": slider_max,
            "time_seconds": time_seconds,
            "default_time_range": default_range,
            "benchmarks": benchmarks,
            "profile_vars": profile_vars,
            "budget_groups": budget_groups,
            "timeheight_vars": timeheight_vars,
            "subcolumn_vars": subcolumn_vars,
            "timeseries_vars": timeseries_vars,
        }
    finally:
        collection.close()


def has_profile_shape(collection, var_name):
    ds_info = collection.dataset_for_var(var_name)
    if ds_info is None:
        return False
    info = ds_info.var_info.get(var_name, {})
    return info.get("z_dim") is not None


def list_timeheight_vars(collection, curated_catalog=None):
    catalog = curated_catalog or TIMEHEIGHT_CATALOG
    options = []
    seen = set()
    for entry in catalog:
        name = entry["name"]
        ds_info = collection.dataset_for_var(name)
        if ds_info is None:
            continue
        info = ds_info.var_info.get(name, {})
        if info.get("t_dim") is None or info.get("z_dim") is None:
            continue
        options.append({"label": entry["title"], "value": name})
        seen.add(name)
    for name in collection.vars:
        if name in seen:
            continue
        ds_info = collection.dataset_for_var(name)
        if ds_info is None:
            continue
        info = ds_info.var_info.get(name, {})
        if info.get("t_dim") is None or info.get("z_dim") is None:
            continue
        meta = timeheight_meta(name)
        options.append({"label": meta["title"], "value": name})
    return options


def has_timeseries_shape(collection, var_name):
    ds_info = collection.dataset_for_var(var_name)
    if ds_info is None:
        return False
    info = ds_info.var_info.get(var_name, {})
    return info.get("t_dim") is not None and info.get("z_dim") is None


def list_timeseries_vars(collection, require_all=False):
    if not require_all:
        return [{"label": name, "value": name} for name in collection.vars if has_timeseries_shape(collection, name)]
    common_vars = None
    for dataset in collection.datasets:
        dataset_vars = {
            name
            for name in dataset.vars
            if dataset.var_info[name]["t_dim"] is not None and dataset.var_info[name]["z_dim"] is None
        }
        common_vars = dataset_vars if common_vars is None else common_vars & dataset_vars
    common_vars = sorted(common_vars or [])
    return [{"label": name, "value": name} for name in common_vars]


def resolve_subcolumn_var(base_name, ds_info):
    return resolve_subcolumn_var_from_names(base_name, ds_info.vars)


def resolve_subcolumn_var_from_names(base_name, vars_found):
    candidates = [
        f"lh_nl_{base_name}",
        f"lh_u_{base_name}",
        base_name,
        f"lh_{base_name}",
        f"lh_{base_name}m",
        f"lh_{base_name}m_mc",
        f"lh_{base_name}_mc",
    ]
    for candidate in candidates:
        if candidate in vars_found:
            return candidate
    return None


def subcolumn_source_for_paths(paths, base_name):
    for path in paths or []:
        meta = dataset_metadata_for_path(path)
        resolved = resolve_subcolumn_var_from_names(base_name, meta["vars"])
        if resolved:
            return path, resolved, meta
    return None, None, None


def list_subcolumn_vars(collection, curated_catalog=None):
    catalog = curated_catalog or SUBCOLUMN_CATALOG
    options = []
    for entry in catalog:
        resolved = None
        for dataset in collection.datasets:
            resolved = resolve_subcolumn_var(entry["name"], dataset)
            if resolved:
                break
        if resolved:
            options.append({"label": entry["label"], "value": entry["name"]})
    return options


def case_has_subcolumn_vars(collection, curated_catalog=None):
    return bool(list_subcolumn_vars(collection, curated_catalog=curated_catalog))


def timeheight_meta(var_name):
    for entry in TIMEHEIGHT_CATALOG:
        if entry["name"] == var_name:
            return entry
    return {"name": var_name, "title": var_name, "axis_title": var_name}


def pick_default_var(options, preferred_names, fallback_index=0):
    values = [option["value"] for option in options]
    for name in preferred_names:
        if name in values:
            return name
    if not values:
        return None
    index = min(max(int(fallback_index), 0), len(values) - 1)
    return values[index]


def get_z_values(ds_info, z_dim):
    key = ("z_values", _file_signature(ds_info.path), z_dim)
    def _build():
        if z_dim in ds_info.ds.variables:
            return as_array(ds_info.ds.variables[z_dim][:])
        return np.arange(len(ds_info.ds.dimensions[z_dim]))
    return _cached_extract(key, _build)


def collection_height_bounds(collection):
    height_min = None
    height_max = None
    for ds_info in collection.datasets:
        seen_z_dims = set()
        for info in ds_info.var_info.values():
            z_dim = info.get("z_dim")
            if z_dim is None or z_dim in seen_z_dims:
                continue
            seen_z_dims.add(z_dim)
            z_vals = np.asarray(get_z_values(ds_info, z_dim), dtype=float)
            if z_vals.size == 0:
                continue
            local_min = float(np.nanmin(z_vals))
            local_max = float(np.nanmax(z_vals))
            height_min = local_min if height_min is None else min(height_min, local_min)
            height_max = local_max if height_max is None else max(height_max, local_max)
    return height_min, height_max


def _tz_plot_data_key(path, var_name, meta):
    info = meta["var_info"][var_name]
    return (
        "tz_plot_data",
        _file_signature(path),
        var_name,
        info["dims"],
        info["t_dim"],
        info["z_dim"],
        info["col_dim"],
    )


def ensure_tz_plot_data(path, var_name):
    meta = dataset_metadata_for_path(path)
    info = meta["var_info"].get(var_name)
    if not info or info["t_dim"] is None or info["z_dim"] is None:
        return None
    key = _tz_plot_data_key(path, var_name, meta)

    def _build():
        ds_info = DatasetInfo(path)
        try:
            info = ds_info.var_info[var_name]
            dims = tuple(info["dims"])
            t_dim = info["t_dim"]
            z_dim = info["z_dim"]
            col_dim = info["col_dim"]
            var = ds_info.ds.variables[var_name]
            slices = []
            kept_dims = []
            for dim in dims:
                lowered = dim.lower()
                if dim == t_dim or dim == z_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                elif dim == col_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                elif lowered in X_DIM_NAMES or lowered in Y_DIM_NAMES:
                    slices.append(0)
                else:
                    slices.append(0)
            data = as_array(var[tuple(slices)])
            ordered_dims = list(kept_dims)
            time_axis = ordered_dims.index(t_dim)
            if time_axis != 0:
                data = np.moveaxis(data, time_axis, 0)
                ordered_dims.insert(0, ordered_dims.pop(time_axis))
            z_axis = ordered_dims.index(z_dim)
            if z_axis != 1:
                data = np.moveaxis(data, z_axis, 1)
                ordered_dims.insert(1, ordered_dims.pop(z_axis))
            remaining_dims = [dim for dim in ordered_dims if dim not in {t_dim, z_dim}]
            if not remaining_dims:
                cube = np.asarray(data, dtype=float)[:, :, np.newaxis]
                line_labels = (var_name,)
                line_columns = (0,)
            else:
                line_labels = tuple(_line_labels_from_dims(remaining_dims, data.shape[2:], col_dim=col_dim))
                cube = np.reshape(data, (data.shape[0], data.shape[1], -1))
                if col_dim in remaining_dims:
                    col_axis = remaining_dims.index(col_dim)
                    line_columns = tuple(int(index_tuple[col_axis]) for index_tuple in np.ndindex(*data.shape[2:]))
                else:
                    line_columns = tuple(0 for _ in range(cube.shape[2]))
            time_vals = time_values_seconds(ds_info)
            if time_vals is None:
                time_vals = np.arange(cube.shape[0], dtype=float)
            return {
                "time_values": np.asarray(time_vals, dtype=float),
                "z_values": np.asarray(get_z_values(ds_info, z_dim), dtype=float),
                "cube": np.asarray(cube, dtype=float),
                "line_labels": line_labels,
                "line_columns": line_columns,
                "columns_len": int(ds_info.columns_len or 1),
                "units": variable_units(ds_info, var_name),
                "long_name": variable_long_name(ds_info, var_name),
                "z_units": dimension_units(ds_info, z_dim),
                "bounds": _finite_bounds(cube),
            }
        finally:
            ds_info.close()

    return _cached_plot_data(key, _build)


def extract_time_avg_profile_for_path(path, var_name, time_range, col_index=0, column_mode="single"):
    data = ensure_tz_plot_data(path, var_name)
    if data is None:
        return None
    cube = np.asarray(data["cube"], dtype=float)
    if cube.ndim != 3 or cube.shape[0] == 0:
        return None
    start_idx, end_idx = int(time_range[0]), int(time_range[1])
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    tmax = cube.shape[0] - 1
    start = max(0, min(start_idx, tmax))
    end = max(0, min(end_idx, tmax))
    if end < start:
        start, end = end, start
    window = cube[start : end + 1, :, :]
    if column_mode == "all":
        profiles = np.mean(window, axis=0).T
        labels = tuple(
            f"col {idx + 1}"
            for idx in range(profiles.shape[0])
        )
        bounds = _finite_bounds(profiles)
    else:
        line_columns = tuple(data["line_columns"])
        matching = [idx for idx, line_col in enumerate(line_columns) if int(line_col) == int(col_index)]
        if not matching:
            matching = [max(0, min(int(col_index), cube.shape[2] - 1))]
        profile = np.mean(window[:, :, matching], axis=(0, 2))
        profiles = np.asarray(profile, dtype=float)[np.newaxis, :]
        labels = (f"col {max(0, min(int(col_index), int(data['columns_len']) - 1)) + 1}",)
        bounds = data["bounds"]
    return {
        "z_values": np.asarray(data["z_values"], dtype=float),
        "profiles": np.asarray(profiles, dtype=float),
        "labels": labels,
        "units": data["units"],
        "long_name": data["long_name"],
        "z_units": data["z_units"],
        "bounds": bounds,
    }


def extract_time_height_for_path(path, var_name, col_index=0, column_mode="single"):
    data = ensure_tz_plot_data(path, var_name)
    if data is None:
        return None
    cube = np.asarray(data["cube"], dtype=float)
    if cube.ndim != 3:
        return None
    if column_mode == "all":
        image = np.mean(cube, axis=2)
    else:
        line_columns = tuple(data["line_columns"])
        matching = [idx for idx, line_col in enumerate(line_columns) if int(line_col) == int(col_index)]
        if not matching:
            matching = [max(0, min(int(col_index), cube.shape[2] - 1))]
        image = np.mean(cube[:, :, matching], axis=2)
    return (
        np.asarray(data["time_values"], dtype=float),
        np.asarray(data["z_values"], dtype=float),
        np.asarray(image, dtype=float),
        data["units"],
        data["long_name"],
        data["z_units"],
    )


def _timeseries_plot_data_key(path, var_name, meta):
    info = meta["var_info"][var_name]
    return (
        "timeseries_plot_data",
        _file_signature(path),
        var_name,
        info["dims"],
        info["t_dim"],
        info["col_dim"],
    )


def ensure_timeseries_plot_data(path, var_name):
    meta = dataset_metadata_for_path(path)
    info = meta["var_info"].get(var_name)
    if not info or info["t_dim"] is None or info["z_dim"] is not None:
        return None
    key = _timeseries_plot_data_key(path, var_name, meta)

    def _build():
        ds_info = DatasetInfo(path)
        try:
            info = ds_info.var_info[var_name]
            dims = tuple(info["dims"])
            t_dim = info["t_dim"]
            col_dim = info["col_dim"]
            var = ds_info.ds.variables[var_name]
            slices = []
            kept_dims = []
            for dim in dims:
                lowered = dim.lower()
                if dim == t_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                elif dim == col_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                elif lowered in X_DIM_NAMES or lowered in Y_DIM_NAMES:
                    slices.append(0)
                else:
                    slices.append(0)
            data = as_array(var[tuple(slices)])
            if data.ndim == 0:
                lines = np.asarray([data], dtype=float)[:, np.newaxis]
                line_labels = (var_name,)
                line_columns = (0,)
            else:
                time_axis = kept_dims.index(t_dim) if t_dim in kept_dims else 0
                if time_axis != 0:
                    data = np.moveaxis(data, time_axis, 0)
                    kept_dims = [kept_dims[time_axis]] + kept_dims[:time_axis] + kept_dims[time_axis + 1 :]
                remaining_dims = [dim for dim in kept_dims if dim != t_dim]
                if data.ndim == 1:
                    lines = np.asarray(data, dtype=float)[:, np.newaxis]
                    line_labels = (var_name,)
                    line_columns = (0,)
                else:
                    line_labels = tuple(_line_labels_from_dims(remaining_dims, data.shape[1:], col_dim=col_dim))
                    lines = np.reshape(data, (data.shape[0], -1))
                    if col_dim in remaining_dims:
                        col_axis = remaining_dims.index(col_dim)
                        line_columns = tuple(int(index_tuple[col_axis]) for index_tuple in np.ndindex(*data.shape[1:]))
                    else:
                        line_columns = tuple(0 for _ in range(lines.shape[1]))
            time_vals = time_values_seconds(ds_info)
            if time_vals is None:
                time_vals = np.arange(lines.shape[0], dtype=float)
            return {
                "time_values": np.asarray(time_vals, dtype=float),
                "lines": np.asarray(lines, dtype=float),
                "line_labels": line_labels,
                "line_columns": line_columns,
                "columns_len": int(ds_info.columns_len or 1),
                "units": variable_units(ds_info, var_name),
                "long_name": variable_long_name(ds_info, var_name),
                "time_units": ds_info.time_units,
            }
        finally:
            ds_info.close()

    return _cached_plot_data(key, _build)


def extract_timeseries_for_path(path, var_name, col_index=0, column_mode="single"):
    data = ensure_timeseries_plot_data(path, var_name)
    if data is None:
        return None
    lines = np.asarray(data["lines"], dtype=float)
    if column_mode == "all":
        selected = lines
        labels = tuple(data["line_labels"])
    else:
        line_columns = tuple(data["line_columns"])
        matching = [idx for idx, line_col in enumerate(line_columns) if int(line_col) == int(col_index)]
        if not matching:
            matching = [max(0, min(int(col_index), lines.shape[1] - 1))]
        selected = lines[:, matching]
        labels = tuple(data["line_labels"][idx] for idx in matching)
    return (
        np.asarray(data["time_values"], dtype=float),
        np.asarray(selected, dtype=float),
        labels,
        data["units"],
        data["long_name"],
        data["time_units"],
    )


def _profile_preload_key(ds_info, var_name, col_index):
    info = ds_info.var_info[var_name]
    dims = tuple(ds_info.ds.variables[var_name].dimensions)
    return ("profile_preload", _file_signature(ds_info.path), var_name, dims, info["t_dim"], info["z_dim"], info["col_dim"], int(col_index))


def get_profile_preload(ds_info, var_name, col_index=0):
    """Return the cached full-time profile series for one variable/column, if present."""
    return _PRELOAD_CACHE.get(_profile_preload_key(ds_info, var_name, col_index))


def ensure_profile_preload(ds_info, var_name, col_index=0):
    """Load and cache the full-time profile series for one variable/column."""
    info = ds_info.var_info[var_name]
    dims = tuple(ds_info.ds.variables[var_name].dimensions)
    t_dim = info["t_dim"]
    z_dim = info["z_dim"]
    col_dim = info["col_dim"]
    key = _profile_preload_key(ds_info, var_name, col_index)

    def _build():
        var = ds_info.ds.variables[var_name]
        slices = []
        kept_dims = []
        for dim in dims:
            if dim == t_dim:
                slices.append(slice(None))
                kept_dims.append(dim)
            elif dim == z_dim:
                slices.append(slice(None))
                kept_dims.append(dim)
            elif dim == col_dim:
                cmax = len(ds_info.ds.dimensions[col_dim]) - 1
                slices.append(max(0, min(int(col_index), cmax)))
            elif dim.lower() in X_DIM_NAMES or dim.lower() in Y_DIM_NAMES:
                slices.append(0)
            else:
                slices.append(0)
        data = as_array(var[tuple(slices)])
        if data.ndim >= 2 and t_dim in kept_dims and z_dim in kept_dims:
            ordered_dims = list(kept_dims)
            time_axis = ordered_dims.index(t_dim)
            if time_axis != 0:
                data = np.moveaxis(data, time_axis, 0)
                ordered_dims.insert(0, ordered_dims.pop(time_axis))
            z_axis = ordered_dims.index(z_dim)
            if z_axis != 1:
                data = np.moveaxis(data, z_axis, 1)
                ordered_dims.insert(1, ordered_dims.pop(z_axis))
            if data.ndim > 2:
                data = np.reshape(data, (data.shape[0], data.shape[1], -1))[:, :, 0]
        elif data.ndim == 1:
            data = np.asarray(data, dtype=float)[:, np.newaxis]
        bounds = _finite_bounds(data)
        return np.asarray(data, dtype=float), bounds

    return _cached_preload(key, _build)


def profile_preload_x_range(ds_info, var_name, col_index=0):
    """Return a stable padded x-axis range derived from the full-time profile preload."""
    _series, bounds = ensure_profile_preload(ds_info, var_name, col_index)
    if bounds is None:
        return None
    return padded_range(bounds[0], bounds[1])


def profile_preload_bounds(ds_info, var_name, col_index=0):
    """Return the raw finite bounds for a cached full-time profile preload."""
    _series, bounds = ensure_profile_preload(ds_info, var_name, col_index)
    return bounds


def _subcolumn_plot_data_key(path, resolved_name, meta):
    info = meta["var_info"][resolved_name]
    return (
        "subcolumn_plot_data",
        _file_signature(path),
        resolved_name,
        info["dims"],
        info["t_dim"],
        info["z_dim"],
        info["col_dim"],
        info["subcol_dim"],
    )


def ensure_subcolumn_plot_data(path, base_name):
    meta = dataset_metadata_for_path(path)
    resolved_name = resolve_subcolumn_var_from_names(base_name, meta["vars"])
    if not resolved_name:
        return None
    info = meta["var_info"].get(resolved_name)
    if not info or info["t_dim"] is None or info["z_dim"] is None:
        return None
    key = _subcolumn_plot_data_key(path, resolved_name, meta)

    def _build():
        ds_info = DatasetInfo(path)
        try:
            info = ds_info.var_info[resolved_name]
            dims = tuple(info["dims"])
            t_dim = info["t_dim"]
            z_dim = info["z_dim"]
            col_dim = info["col_dim"]
            subcol_dim = info["subcol_dim"]
            var = ds_info.ds.variables[resolved_name]
            slices = []
            kept_dims = []
            for dim in dims:
                lowered = dim.lower()
                if dim == t_dim or dim == z_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                elif dim == col_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                elif lowered in X_DIM_NAMES or lowered in Y_DIM_NAMES:
                    slices.append(0)
                else:
                    slices.append(slice(None))
                    kept_dims.append(dim)
            data = as_array(var[tuple(slices)])
            ordered_dims = list(kept_dims)
            time_axis = ordered_dims.index(t_dim)
            if time_axis != 0:
                data = np.moveaxis(data, time_axis, 0)
                ordered_dims.insert(0, ordered_dims.pop(time_axis))
            z_axis = ordered_dims.index(z_dim)
            if z_axis != 1:
                data = np.moveaxis(data, z_axis, 1)
                ordered_dims.insert(1, ordered_dims.pop(z_axis))
            remaining_dims = [dim for dim in ordered_dims if dim not in {t_dim, z_dim}]
            if not remaining_dims:
                profiles = np.asarray(data, dtype=float)[:, :, np.newaxis]
                labels = (resolved_name,)
                line_columns = (0,)
                has_true_subcolumns = False
            else:
                labels = tuple(_line_labels_from_dims(remaining_dims, data.shape[2:], col_dim=col_dim, subcol_dim=subcol_dim))
                profiles = np.reshape(data, (data.shape[0], data.shape[1], -1))
                if col_dim in remaining_dims:
                    col_axis = remaining_dims.index(col_dim)
                    line_columns = tuple(int(index_tuple[col_axis]) for index_tuple in np.ndindex(*data.shape[2:]))
                else:
                    line_columns = tuple(0 for _ in range(profiles.shape[2]))
                has_true_subcolumns = bool(subcol_dim is not None and subcol_dim in remaining_dims)
            return {
                "resolved_name": resolved_name,
                "z_values": np.asarray(get_z_values(ds_info, z_dim), dtype=float),
                "profiles": np.asarray(profiles, dtype=float),
                "labels": labels,
                "line_columns": line_columns,
                "has_true_subcolumns": has_true_subcolumns,
                "units": variable_units(ds_info, resolved_name),
                "long_name": variable_long_name(ds_info, resolved_name),
                "z_units": dimension_units(ds_info, z_dim),
                "bounds": _finite_bounds(profiles),
            }
        finally:
            ds_info.close()

    return _cached_plot_data(key, _build)


def subcolumn_x_range_for_path(path, base_name, col_index=0, column_mode="single"):
    data = ensure_subcolumn_plot_data(path, base_name)
    if data is None:
        return None
    if column_mode == "all":
        bounds = data["bounds"]
    else:
        line_columns = tuple(data["line_columns"])
        matching = [idx for idx, line_col in enumerate(line_columns) if int(line_col) == int(col_index)]
        if not matching:
            matching = [max(0, min(int(col_index), len(line_columns) - 1))]
        bounds = _finite_bounds(np.asarray(data["profiles"], dtype=float)[:, :, matching])
    if bounds is None:
        return None
    return padded_range(bounds[0], bounds[1])


def extract_time_avg_profile(ds_info, var_name, time_range, col_index=0):
    info = ds_info.var_info[var_name]
    dims = tuple(ds_info.ds.variables[var_name].dimensions)
    t_dim = info["t_dim"]
    z_dim = info["z_dim"]
    col_dim = info["col_dim"]
    start_idx, end_idx = int(time_range[0]), int(time_range[1])
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    key = ("profile", _file_signature(ds_info.path), var_name, dims, t_dim, z_dim, col_dim, start_idx, end_idx, int(col_index))

    preload = get_profile_preload(ds_info, var_name, col_index)
    if preload is not None:
        full_data, _bounds = preload
        tmax = full_data.shape[0] - 1
        start = max(0, min(start_idx, tmax))
        end = max(0, min(end_idx, tmax))
        if end < start:
            start, end = end, start
        window = np.asarray(full_data[start : end + 1, :], dtype=float)
        if window.ndim == 1:
            window = window[np.newaxis, :]
        return get_z_values(ds_info, z_dim), np.asarray(np.mean(window, axis=0), dtype=float)

    def _build():
        var = ds_info.ds.variables[var_name]
        slices = []
        kept_dims = []
        for dim in dims:
            if dim == t_dim:
                tmax = len(ds_info.ds.dimensions[t_dim]) - 1
                start = max(0, min(start_idx, tmax))
                end = max(0, min(end_idx, tmax))
                if end < start:
                    start, end = end, start
                slices.append(slice(start, end + 1))
                kept_dims.append(dim)
            elif dim == z_dim:
                slices.append(slice(None))
                kept_dims.append(dim)
            elif dim == col_dim:
                cmax = len(ds_info.ds.dimensions[col_dim]) - 1
                slices.append(max(0, min(int(col_index), cmax)))
            elif dim.lower() in X_DIM_NAMES or dim.lower() in Y_DIM_NAMES:
                slices.append(0)
            else:
                slices.append(0)
        data = as_array(var[tuple(slices)])
        if t_dim in kept_dims:
            data = mean_over_axis(data, kept_dims.index(t_dim))
        if data.ndim > 1:
            z_axis = kept_dims.index(z_dim) if z_dim in kept_dims else 0
            data = np.moveaxis(data, z_axis, 0)
            data = np.reshape(data, (data.shape[0], -1))[:, 0]
        return get_z_values(ds_info, z_dim), np.asarray(data)

    return _cached_extract(key, _build)


def extract_time_height(ds_info, var_name, col_index=0, column_mode="single"):
    info = ds_info.var_info[var_name]
    t_dim = info["t_dim"]
    z_dim = info["z_dim"]
    col_dim = info["col_dim"]
    if t_dim is None or z_dim is None:
        return None
    dims = tuple(ds_info.ds.variables[var_name].dimensions)
    key = ("timeheight", _file_signature(ds_info.path), var_name, dims, t_dim, z_dim, col_dim, int(col_index), column_mode)

    def _build():
        var = ds_info.ds.variables[var_name]
        slices = []
        kept_dims = []
        for dim in dims:
            if dim == t_dim or dim == z_dim:
                slices.append(slice(None))
                kept_dims.append(dim)
            elif dim == col_dim:
                if column_mode == "all":
                    slices.append(slice(None))
                    kept_dims.append(dim)
                else:
                    cmax = len(ds_info.ds.dimensions[col_dim]) - 1
                    slices.append(max(0, min(int(col_index), cmax)))
            elif dim.lower() in X_DIM_NAMES or dim.lower() in Y_DIM_NAMES:
                slices.append(0)
            else:
                slices.append(0)
        data = as_array(var[tuple(slices)])
        if col_dim in kept_dims:
            col_axis = kept_dims.index(col_dim)
            data = mean_over_axis(data, col_axis)
            kept_dims = [dim for dim in kept_dims if dim != col_dim]
        if data.ndim != 2:
            return None
        time_axis = kept_dims.index(t_dim)
        z_axis = kept_dims.index(z_dim)
        if time_axis != 0:
            data = np.moveaxis(data, time_axis, 0)
            if z_axis == 0:
                z_axis = 1
            else:
                z_axis -= 1
        if z_axis != 1:
            data = np.moveaxis(data, z_axis, 1)
        time_vals = time_values_seconds(ds_info)
        if time_vals is None:
            time_vals = np.arange(data.shape[0])
        z_vals = get_z_values(ds_info, z_dim)
        return time_vals, z_vals, data

    return _cached_extract(key, _build)


def extract_timeseries(ds_info, var_name, col_index=0, column_mode="single"):
    info = ds_info.var_info[var_name]
    t_dim = info["t_dim"]
    z_dim = info["z_dim"]
    col_dim = info["col_dim"]
    if t_dim is None or z_dim is not None:
        return None, None, None
    dims = tuple(ds_info.ds.variables[var_name].dimensions)
    key = ("timeseries", _file_signature(ds_info.path), var_name, dims, t_dim, col_dim, int(col_index), column_mode)

    def _build():
        var = ds_info.ds.variables[var_name]
        slices = []
        kept_dims = []
        for dim in dims:
            lowered = dim.lower()
            if dim == t_dim:
                slices.append(slice(None))
                kept_dims.append(dim)
            elif dim == col_dim:
                if column_mode == "all":
                    slices.append(slice(None))
                    kept_dims.append(dim)
                else:
                    cmax = len(ds_info.ds.dimensions[col_dim]) - 1
                    slices.append(max(0, min(int(col_index), cmax)))
            elif lowered in X_DIM_NAMES or lowered in Y_DIM_NAMES:
                slices.append(0)
            else:
                slices.append(0)
        data = as_array(var[tuple(slices)])
        if data.ndim == 0:
            data = np.asarray([data], dtype=float)
            time_vals = np.asarray([0.0], dtype=float)
            return time_vals, data[:, np.newaxis], (var_name,)
        time_axis = kept_dims.index(t_dim) if t_dim in kept_dims else 0
        if time_axis != 0:
            data = np.moveaxis(data, time_axis, 0)
            kept_dims = [kept_dims[time_axis]] + kept_dims[:time_axis] + kept_dims[time_axis + 1 :]
        remaining_dims = [dim for dim in kept_dims if dim != t_dim]
        if data.ndim == 1:
            lines = np.asarray(data)[:, np.newaxis]
            labels = (var_name,)
        else:
            labels = tuple(_line_labels_from_dims(remaining_dims, data.shape[1:], col_dim=col_dim))
            lines = np.reshape(data, (data.shape[0], -1))
        time_vals = time_values_seconds(ds_info)
        if time_vals is None:
            time_vals = np.arange(lines.shape[0])
        return np.asarray(time_vals), np.asarray(lines), labels

    return _cached_extract(key, _build)


def _line_labels_from_dims(remaining_dims, shape, col_dim=None, subcol_dim=None):
    if not remaining_dims or not shape:
        return ["line 1"]
    labels = []
    for index_tuple in np.ndindex(*shape):
        parts = []
        for dim_name, dim_index in zip(remaining_dims, index_tuple):
            if dim_name == col_dim:
                parts.append(f"col {dim_index + 1}")
            elif dim_name == subcol_dim:
                parts.append(f"subcol {dim_index + 1}")
            else:
                parts.append(f"{dim_name} {dim_index + 1}")
        labels.append(" / ".join(parts))
    return labels


def extract_subcolumn_profiles(ds_info, base_name, time_range, col_index=0, column_mode="single"):
    return extract_subcolumn_profiles_for_path(ds_info.path, base_name, time_range, col_index=col_index, column_mode=column_mode)


def extract_subcolumn_profiles_for_path(path, base_name, time_range, col_index=0, column_mode="single"):
    data = ensure_subcolumn_plot_data(path, base_name)
    if data is None:
        return None, None, None, False
    start_idx, end_idx = int(time_range[0]), int(time_range[1])
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    full_profiles = np.asarray(data["profiles"], dtype=float)
    tmax = full_profiles.shape[0] - 1
    start = max(0, min(start_idx, tmax))
    end = max(0, min(end_idx, tmax))
    if end < start:
        start, end = end, start
    line_indices = list(range(full_profiles.shape[2]))
    if column_mode != "all":
        line_columns = tuple(data["line_columns"])
        line_indices = [idx for idx, line_col in enumerate(line_columns) if int(line_col) == int(col_index)]
        if not line_indices:
            line_indices = [max(0, min(int(col_index), full_profiles.shape[2] - 1))]
    window = full_profiles[start : end + 1, :, line_indices]
    if window.ndim == 2:
        window = window[np.newaxis, :, :]
    averaged = np.asarray(np.mean(window, axis=0), dtype=float)
    labels = tuple(data["labels"][idx] for idx in line_indices)
    return np.asarray(data["z_values"], dtype=float), averaged, labels, bool(data["has_true_subcolumns"])


def _decode_param_names(ds):
    if "param_name" not in ds.variables:
        return []
    var = ds.variables["param_name"]
    raw = as_array(var[:])
    if raw.ndim != 2:
        return []
    dims = list(var.dimensions)
    if "param" in dims:
        param_axis = dims.index("param")
        if param_axis != 0:
            raw = np.moveaxis(raw, param_axis, 0)
    names = chartostring(raw)
    if names.ndim != 1:
        names = np.reshape(names, (-1,))
    decoded = []
    for name in names.tolist():
        if isinstance(name, bytes):
            name = name.decode("utf-8", errors="replace")
        decoded.append(str(name).strip())
    return decoded


def _load_param_values_for_path(path):
    """Load clubb_params and param_name from one stats file."""
    try:
        ds_info = DatasetInfo(path)
    except OSError:
        return {}, False, False
    try:
        ds = ds_info.ds
        if "clubb_params" not in ds.variables:
            return {}, False, False
        param_names = _decode_param_names(ds)
        if not param_names:
            return {}, True, False
        var = ds.variables["clubb_params"]
        data = as_array(var[:])
        dims = list(var.dimensions)
        if data.ndim != 2:
            return {}, True, bool(param_names)
        param_axis = dims.index("param") if "param" in dims else 0
        col_axis = None
        for idx, dim in enumerate(dims):
            if dim.lower() in COL_DIM_NAMES:
                col_axis = idx
                break
        if col_axis is None:
            col_axis = 1 if param_axis == 0 else 0
        if param_axis != 0 or col_axis != 1:
            data = np.moveaxis(data, [param_axis, col_axis], [0, 1])
        nparams = min(data.shape[0], len(param_names))
        return {
            param_names[i]: np.asarray(data[i, :], dtype=float).tolist()
            for i in range(nparams)
            if param_names[i]
        }, True, True
    finally:
        ds_info.close()


def load_param_values(paths):
    """Load parameter values from the first valid stats file in a path list."""
    for path in paths or []:
        params, has_clubb_params, has_param_names = _load_param_values_for_path(path)
        if has_clubb_params or has_param_names or params:
            return params, has_clubb_params, has_param_names
    return {}, False, False


def load_compare_param_values(paths):
    """Load and reconcile parameter values across multiple compared outputs."""
    per_file = []
    any_clubb_params = False
    all_have_clubb_params = True
    all_have_param_names = True
    for path in paths or []:
        params, has_clubb_params, has_param_names = _load_param_values_for_path(path)
        per_file.append(params)
        any_clubb_params = any_clubb_params or has_clubb_params
        all_have_clubb_params = all_have_clubb_params and has_clubb_params
        all_have_param_names = all_have_param_names and has_param_names

    if not per_file:
        return {}, False, False, [], []

    if not all_have_clubb_params:
        return {}, any_clubb_params, all_have_param_names, [], []
    if not all_have_param_names:
        return {}, True, False, [], []

    name_sets = [set(params.keys()) for params in per_file]
    common_names = sorted(set.intersection(*name_sets)) if name_sets else []
    matched = {}
    mismatched = []
    for name in common_names:
        reference = np.asarray(per_file[0][name], dtype=float)
        same = True
        for params in per_file[1:]:
            candidate = np.asarray(params[name], dtype=float)
            if reference.shape != candidate.shape or not np.allclose(reference, candidate, rtol=0.0, atol=0.0):
                same = False
                break
        if same:
            matched[name] = reference.tolist()
        else:
            mismatched.append((name, "Mismatch across compared outputs"))
    missing_names = sorted(set.union(*name_sets) - set(common_names)) if name_sets else []
    mismatched.extend((name, "Missing in at least one compared output") for name in missing_names)
    return matched, True, True, mismatched, per_file


def closest_column(param_values, selection):
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

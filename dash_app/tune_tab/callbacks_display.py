"""Callbacks for rendering tuning status and result tables."""

from __future__ import annotations

import json
import math
from pathlib import Path

from dash import ALL, Input, Output, State, dcc, html, no_update
import numpy as np
import plotly.graph_objects as go

from .layout import (
    action_button_style,
    build_results_placeholder,
    config_button_style,
    mode_button_style,
    mode_options_block_style,
)
from tuner.sample_history import sample_history_paths
from tuner.taylor_metrics import (
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_LOSS_MODE,
    LOSS_METRIC_NAMES,
    aggregate_losses,
    compute_field_loss_diagnostics,
)


_CORRELATION_GRID = [-1.0, -0.75, -0.5, 0.0, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
_POSITIVE_CORRELATION_GRID = [0.0, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
_FIELD_MARKER_SYMBOLS = [
    "circle",
    "square",
    "diamond",
    "cross",
    "x",
    "triangle-up",
    "triangle-down",
    "triangle-right",
    "triangle-left",
    "pentagon",
    "hexagon",
    "hexagon2",
    "hourglass",
    "bowtie",
]
_DIAGNOSTIC_PENALTY_THRESHOLD = 1.0e20
_TAYLOR_PLOT_VALUE_LIMIT = 50.0
_PARAMETER_GROUP_COLORS = [
    "#2563eb",
    "#dc2626",
    "#16a34a",
    "#9333ea",
    "#ea580c",
    "#0891b2",
    "#be123c",
    "#4f46e5",
]
_LANDSCAPE_ALL = "__all__"
_LANDSCAPE_METADATA_NAMES = [
    "schema_version",
    "param_names",
    "case_names",
    "field_names",
    "metric_names",
    "case_window_counts",
    "obs_case",
    "obs_field",
    "obs_window",
    "obs_window_start_seconds",
    "obs_window_end_seconds",
]
_LANDSCAPE_AGGREGATION_NAMES = {"mean", "median", "min", "p90", "max"}
_LANDSCAPE_MODE_NAMES = {"samples", "binned"}
_LANDSCAPE_BINS = 28
_LANDSCAPE_LOSS_LIKE_METRICS = {"total_loss", "raw:scaled_rmse", "raw:centered_rmse_norm"}


def _finite_float(value):
    """Return a finite float, or None for absent/non-finite values."""
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result):
        return None
    return result


def _plot_diagnostic_float(value):
    """Return a finite, non-penalty diagnostic value for plotting."""
    result = _finite_float(value)
    if result is None or abs(result) >= _DIAGNOSTIC_PENALTY_THRESHOLD:
        return None
    return result


def _plot_metric_mask(values):
    """Return values finite enough to scale landscape plots."""
    array = np.asarray(values, dtype=float)
    return np.isfinite(array) & (np.abs(array) < _DIAGNOSTIC_PENALTY_THRESHOLD)


def _level_step(max_value):
    """Choose a readable grid interval for Taylor-space distances."""
    if max_value <= 2.5:
        return 0.5
    if max_value <= 5.0:
        return 1.0
    return 2.0


def _levels_up_to(max_value):
    """Return positive grid levels through max_value."""
    step = _level_step(max_value)
    count = int(math.ceil(max_value / step))
    count = min(count, 32)
    return [round(step * idx, 6) for idx in range(1, count + 1)]


def _arc_points(radius, theta_start, theta_end, count=100):
    """Return x/y points for an origin-centered Taylor std-ratio arc."""
    if count < 2:
        count = 2
    theta_values = [
        theta_start + (theta_end - theta_start) * idx / (count - 1)
        for idx in range(count)
    ]
    return (
        [radius * math.cos(theta) for theta in theta_values],
        [radius * math.sin(theta) for theta in theta_values],
    )


def _crmse_points(radius, x_min, x_max, y_max, count=160):
    """Return visible points for a centered-RMSE contour around the reference point."""
    x_values = []
    y_values = []
    for idx in range(count):
        theta = math.pi * idx / (count - 1)
        x_value = 1.0 + radius * math.cos(theta)
        y_value = radius * math.sin(theta)
        if x_min <= x_value <= x_max and 0.0 <= y_value <= y_max:
            x_values.append(x_value)
            y_values.append(y_value)
        elif x_values and x_values[-1] is not None:
            x_values.append(None)
            y_values.append(None)
    return x_values, y_values


def _collect_taylor_points(best_results):
    """Flatten retained best results into Taylor diagram point records."""
    points = []
    for result_index, result in enumerate(best_results or [], start=1):
        rank = int(result.get("rank", result_index))
        total_loss = _plot_diagnostic_float(result.get("total_loss"))
        if total_loss is None:
            continue
        field_metrics = result.get("field_metrics", {})
        if not isinstance(field_metrics, dict):
            continue
        for case_name, case_fields in sorted(field_metrics.items()):
            if not isinstance(case_fields, dict):
                continue
            for field_name, metrics in sorted(case_fields.items()):
                if not isinstance(metrics, dict):
                    continue
                correlation = _plot_diagnostic_float(metrics.get("correlation"))
                std_ratio = _plot_diagnostic_float(metrics.get("std_ratio"))
                centered_rmse_norm = _plot_diagnostic_float(metrics.get("centered_rmse_norm"))
                bias_norm = _plot_diagnostic_float(metrics.get("bias_norm"))
                field_loss = _plot_diagnostic_float(metrics.get("loss"))
                scaled_rmse = _plot_diagnostic_float(metrics.get("scaled_rmse", metrics.get("simple_rms")))
                if (
                    correlation is None
                    or std_ratio is None
                    or centered_rmse_norm is None
                    or bias_norm is None
                    or field_loss is None
                ):
                    continue
                correlation = max(-1.0, min(1.0, correlation))
                if (
                    std_ratio < 0.0
                    or std_ratio > _TAYLOR_PLOT_VALUE_LIMIT
                    or centered_rmse_norm > _TAYLOR_PLOT_VALUE_LIMIT
                    or abs(bias_norm) > _TAYLOR_PLOT_VALUE_LIMIT
                ):
                    continue
                points.append(
                    {
                        "rank": rank,
                        "case": str(case_name),
                        "field": str(field_name),
                        "correlation": correlation,
                        "std_ratio": std_ratio,
                        "centered_rmse_norm": centered_rmse_norm,
                        "bias_norm": bias_norm,
                        "field_loss": field_loss,
                        "scaled_rmse": scaled_rmse,
                        "total_loss": total_loss,
                        "x": std_ratio * correlation,
                        "y": std_ratio * math.sqrt(max(0.0, 1.0 - correlation * correlation)),
                    }
                )
    return points


def _fmt_metric(value, precision=4):
    """Format a float for hover text."""
    if value is None:
        return "--"
    return f"{value:.{precision}g}"


def _taylor_hover_text(point):
    """Return hover text for one Taylor diagnostic point."""
    return "<br>".join(
        [
            f"Rank {point['rank']}",
            f"Case {point['case']}",
            f"Field {point['field']}",
            f"Total smart loss {_fmt_metric(point['total_loss'], 6)}",
            f"Field smart loss {_fmt_metric(point['field_loss'], 6)}",
            f"scaled_rmse {_fmt_metric(point['scaled_rmse'], 6)}",
            f"correlation {_fmt_metric(point['correlation'])}",
            f"std_ratio {_fmt_metric(point['std_ratio'])}",
            f"centered_rmse_norm {_fmt_metric(point['centered_rmse_norm'])}",
            f"bias_norm {_fmt_metric(point['bias_norm'])}",
        ]
    )


def _empty_diagnostics_figure(title, message):
    """Return a stable empty diagnostics figure."""
    fig = go.Figure()
    fig.update_layout(
        title={"text": title, "x": 0.02, "xanchor": "left"},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 58, "r": 18, "t": 48, "b": 54},
        height=430,
        annotations=[
            {
                "text": message,
                "xref": "paper",
                "yref": "paper",
                "x": 0.5,
                "y": 0.5,
                "showarrow": False,
                "font": {"size": 13, "color": "#64748b"},
            }
        ],
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig


def build_taylor_figure(best_results):
    """Build a Taylor-diagram-ready diagnostics figure for retained top results."""
    points = _collect_taylor_points(best_results)
    if not points:
        return _empty_diagnostics_figure(
            "Taylor Diagnostics",
            "Taylor diagnostics will appear after tuning results are available.",
        )

    show_negative = any(point["correlation"] < -1.0e-8 for point in points)
    theta_end = math.pi if show_negative else math.pi / 2.0
    max_std_ratio = max([1.0] + [point["std_ratio"] for point in points])
    max_y = max([1.0] + [point["y"] for point in points])
    max_x_abs = max([1.0] + [abs(point["x"]) for point in points])
    radius_limit = max(1.5, 1.15 * max(max_std_ratio, max_y, max_x_abs))
    x_min = -radius_limit if show_negative else 0.0
    x_max = radius_limit
    y_max = radius_limit

    fig = go.Figure()
    grid_color = "rgba(71, 85, 105, 0.28)"
    crmse_color = "rgba(202, 138, 4, 0.34)"

    for radius in _levels_up_to(radius_limit):
        x_values, y_values = _arc_points(radius, 0.0, theta_end)
        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="lines",
                line={"color": grid_color, "width": 1},
                hoverinfo="skip",
                showlegend=False,
            )
        )
        if radius <= radius_limit:
            label_theta = theta_end * 0.55
            fig.add_annotation(
                x=radius * math.cos(label_theta),
                y=radius * math.sin(label_theta),
                text=f"{radius:g}",
                showarrow=False,
                font={"size": 10, "color": "#64748b"},
                xanchor="center",
                yanchor="middle",
            )

    correlation_grid = _CORRELATION_GRID if show_negative else _POSITIVE_CORRELATION_GRID
    for correlation in correlation_grid:
        theta = math.acos(max(-1.0, min(1.0, correlation)))
        if theta > theta_end + 1.0e-12:
            continue
        x_end = radius_limit * math.cos(theta)
        y_end = radius_limit * math.sin(theta)
        fig.add_trace(
            go.Scatter(
                x=[0.0, x_end],
                y=[0.0, y_end],
                mode="lines",
                line={"color": grid_color, "width": 1, "dash": "dot"},
                hoverinfo="skip",
                showlegend=False,
            )
        )
        fig.add_annotation(
            x=1.04 * x_end,
            y=1.04 * y_end,
            text=f"{correlation:g}",
            showarrow=False,
            font={"size": 10, "color": "#334155"},
            xanchor="center",
            yanchor="middle",
        )

    max_crmse = max([1.0] + [point["centered_rmse_norm"] or 0.0 for point in points])
    for radius in _levels_up_to(max(radius_limit, max_crmse)):
        x_values, y_values = _crmse_points(radius, x_min, x_max, y_max)
        if len([value for value in x_values if value is not None]) < 2:
            continue
        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="lines",
                line={"color": crmse_color, "width": 1, "dash": "dash"},
                hoverinfo="skip",
                showlegend=False,
            )
        )

    rank_max = max(point["rank"] for point in points)
    fig.add_trace(
        go.Scatter(
            x=[1.0],
            y=[0.0],
            mode="markers",
            marker={"symbol": "star", "size": 14, "color": "#dc2626", "line": {"color": "#7f1d1d", "width": 1}},
            name="Reference",
            hovertemplate="Reference<br>std ratio 1<br>correlation 1<extra></extra>",
        )
    )

    field_names = sorted({point["field"] for point in points})
    field_symbols = {
        field_name: _FIELD_MARKER_SYMBOLS[idx % len(_FIELD_MARKER_SYMBOLS)]
        for idx, field_name in enumerate(field_names)
    }
    for idx, field_name in enumerate(field_names):
        field_points = [point for point in points if point["field"] == field_name]
        fig.add_trace(
            go.Scatter(
                x=[point["x"] for point in field_points],
                y=[point["y"] for point in field_points],
                mode="markers",
                marker={
                    "symbol": field_symbols[field_name],
                    "size": 10,
                    "color": [point["rank"] for point in field_points],
                    "colorscale": "Viridis_r",
                    "cmin": 1,
                    "cmax": max(2, rank_max),
                    "colorbar": {"title": "Rank", "len": 0.8},
                    "showscale": idx == 0,
                    "line": {"color": "#0f172a", "width": 0.7},
                    "opacity": 0.88,
                },
                text=[_taylor_hover_text(point) for point in field_points],
                hovertemplate="%{text}<extra></extra>",
                name=field_name,
            )
        )

    fig.update_layout(
        title={"text": "Taylor Diagnostics", "x": 0.02, "xanchor": "left"},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 58, "r": 18, "t": 48, "b": 54},
        height=430,
        hovermode="closest",
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1.0},
        uirevision="tune-taylor-diagnostics",
    )
    fig.update_xaxes(
        range=[x_min, x_max],
        zeroline=True,
        zerolinecolor="#94a3b8",
        gridcolor="rgba(148, 163, 184, 0.24)",
        title="std ratio * correlation",
    )
    fig.update_yaxes(
        range=[0.0, y_max],
        scaleanchor="x",
        scaleratio=1,
        zeroline=True,
        zerolinecolor="#94a3b8",
        gridcolor="rgba(148, 163, 184, 0.24)",
        title="centered std ratio",
    )

    return fig


def build_taylor_diagram(best_results):
    """Render Taylor-diagram-ready diagnostics for the retained top results."""
    return dcc.Graph(
        id="tune-taylor-diagram",
        figure=build_taylor_figure(best_results),
        className="tune-taylor-graph",
        config={"responsive": True, "displaylogo": False},
        style={"width": "100%", "minWidth": 0, "height": "430px"},
    )


def _result_params(result):
    params = result.get("selected_params")
    if not isinstance(params, dict):
        params = result.get("params", {})
    return params if isinstance(params, dict) else {}


def _parameter_group_specs(best_results, best_results_by_case=None, selected_groups=None):
    """Return selected parameter-spread groups in plotting order."""
    by_case = best_results_by_case if isinstance(best_results_by_case, dict) else {}
    selected = list(selected_groups or ["aggregate"])
    specs = []
    if "aggregate" in selected:
        specs.append(("aggregate", "Aggregate", list(best_results or []), "total_loss"))
    for case_name in sorted(by_case):
        key = f"case:{case_name}"
        if key in selected:
            specs.append((key, str(case_name), list(by_case.get(case_name) or []), "case_loss"))
    return specs


def parameter_box_group_options(best_results_by_case, current_value=None):
    """Build toggle options and a valid selected group list for the parameter box plot."""
    by_case = best_results_by_case if isinstance(best_results_by_case, dict) else {}
    options = [{"label": "Aggregate", "value": "aggregate"}]
    options.extend({"label": str(case_name), "value": f"case:{case_name}"} for case_name in sorted(by_case))
    valid_values = {option["value"] for option in options}
    selected = [value for value in (current_value or ["aggregate"]) if value in valid_values]
    if not selected:
        selected = ["aggregate"]
    return options, selected


def _collect_top_parameter_values(best_results, best_results_by_case=None, selected_groups=None):
    """Return top-result parameter values grouped by selected parameter name and result group."""
    param_names = []
    grouped = {}
    for group_key, group_label, group_results, loss_key in _parameter_group_specs(best_results, best_results_by_case, selected_groups):
        for result_index, result in enumerate(group_results, start=1):
            rank = int(result.get("rank", result_index))
            group_loss = _finite_float(result.get(loss_key))
            total_loss = _finite_float(result.get("total_loss"))
            scaled_rmse_sum = _finite_float(result.get("scaled_rmse_sum", result.get("simple_rms_sum")))
            params = _result_params(result)
            for name, raw_value in params.items():
                value = _finite_float(raw_value)
                if value is None:
                    continue
                if name not in grouped:
                    param_names.append(name)
                    grouped[name] = {}
                if group_key not in grouped[name]:
                    grouped[name][group_key] = {
                        "label": group_label,
                        "values": [],
                        "hover": [],
                    }
                grouped[name][group_key]["values"].append(value)
                grouped[name][group_key]["hover"].append(
                    "<br>".join(
                        [
                            f"{group_label} rank {rank}",
                            f"{name} {_fmt_metric(value, 6)}",
                            f"{loss_key} {_fmt_metric(group_loss, 6)}",
                            f"Total smart loss {_fmt_metric(total_loss, 6)}",
                            f"Scaled RMSE sum {_fmt_metric(scaled_rmse_sum, 6)}",
                        ]
                    )
                )
    return param_names, grouped


def build_parameter_box_figure(best_results, best_results_by_case=None, selected_groups=None):
    """Build a box plot figure of selected parameter values in retained top results."""
    param_names, grouped_values = _collect_top_parameter_values(best_results, best_results_by_case, selected_groups)
    if not param_names:
        return _empty_diagnostics_figure(
            "Top-16 Parameter Spread",
            "Parameter distributions will appear after tuning results are available.",
        )

    fig = go.Figure()
    group_keys = [group_key for group_key, _label, _results, _loss_key in _parameter_group_specs(best_results, best_results_by_case, selected_groups)]
    legend_seen = set()
    for group_idx, group_key in enumerate(group_keys):
        color = _PARAMETER_GROUP_COLORS[group_idx % len(_PARAMETER_GROUP_COLORS)]
        for name in param_names:
            group_data = grouped_values.get(name, {}).get(group_key)
            if not group_data:
                continue
            fig.add_trace(
                go.Box(
                    x=[name] * len(group_data["values"]),
                    y=group_data["values"],
                    name=group_data["label"],
                    legendgroup=group_key,
                    showlegend=group_key not in legend_seen,
                    boxpoints="all",
                    jitter=0.35,
                    pointpos=0,
                    marker={"size": 7, "color": color, "opacity": 0.72},
                    line={"color": color, "width": 1.2},
                    fillcolor="rgba(148, 163, 184, 0.20)",
                    text=group_data["hover"],
                    hovertemplate="%{text}<extra></extra>",
                )
            )
            legend_seen.add(group_key)

    fig.update_layout(
        title={"text": "Top-16 Parameter Spread", "x": 0.02, "xanchor": "left"},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 58, "r": 18, "t": 48, "b": 68},
        height=430,
        showlegend=True,
        boxmode="group",
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1.0},
        uirevision="tune-parameter-box-plot",
    )
    fig.update_xaxes(
        tickangle=-25,
        gridcolor="rgba(148, 163, 184, 0.18)",
    )
    fig.update_yaxes(
        title="parameter value",
        zeroline=True,
        zerolinecolor="#94a3b8",
        gridcolor="rgba(148, 163, 184, 0.24)",
    )

    return fig


def build_parameter_box_plot(best_results, best_results_by_case=None, selected_groups=None):
    """Render a box plot of selected parameter values in the retained top results."""
    return dcc.Graph(
        id="tune-parameter-box-plot",
        figure=build_parameter_box_figure(best_results, best_results_by_case, selected_groups),
        className="tune-parameter-box-graph",
        config={"responsive": True, "displaylogo": False},
        style={"width": "100%", "minWidth": 0, "height": "430px"},
    )


def build_diagnostics_row(best_results):
    """Render the side-by-side Taylor and parameter-spread diagnostics."""
    return html.Div(
        [
            html.Div(build_taylor_diagram(best_results), className="tune-diagnostics-taylor"),
            html.Div(build_parameter_box_plot(best_results), className="tune-diagnostics-params"),
        ],
        className="tune-diagnostics-row",
    )


def _job_dir_from_status(status):
    """Return the current tuner job directory from a status payload."""
    job_dir = (status or {}).get("job_dir")
    if not job_dir:
        return None
    return Path(job_dir)


def _sample_history_signature(job_dir):
    """Return a cheap signature for all available sample-history chunks."""
    if job_dir is None:
        return ""
    try:
        paths = sample_history_paths(Path(job_dir))
    except Exception:
        return ""
    parts = []
    for path in paths:
        try:
            stat = path.stat()
        except OSError:
            continue
        parts.append(f"{path.name}:{stat.st_size}:{stat.st_mtime_ns}")
    return "|".join(parts)


def _load_sample_history(job_dir):
    """Load sample-history chunks from one job directory."""
    if job_dir is None:
        return None
    paths = sample_history_paths(Path(job_dir))
    if not paths:
        return None

    chunks = []
    try:
        chunks = [np.load(path) for path in paths]
        first = chunks[0]
        history = {
            name: first[name].copy()
            for name in _LANDSCAPE_METADATA_NAMES
        }
        history["all_params"] = np.concatenate([chunk["all_params"] for chunk in chunks], axis=0)
        history["loss_metrics"] = np.concatenate([chunk["loss_metrics"] for chunk in chunks], axis=0)
        history["sample_id"] = np.concatenate([chunk["sample_id"] for chunk in chunks], axis=0)
        history["batch_id"] = np.concatenate([chunk["batch_id"] for chunk in chunks], axis=0)
        history["chunk_paths"] = paths
        return history
    finally:
        for chunk in chunks:
            try:
                chunk.close()
            except Exception:
                pass


def _read_job_request(job_dir):
    """Return the request.json payload for a tuner job, if available."""
    if job_dir is None:
        return {}
    try:
        with (Path(job_dir) / "request.json").open(encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _string_array_values(history, name):
    """Return a metadata string array as plain Python strings."""
    return [str(value) for value in list((history or {}).get(name, []))]


def _clean_param_names(values):
    """Return non-empty parameter names preserving order."""
    names = []
    seen = set()
    for value in values or []:
        if not isinstance(value, str):
            continue
        name = value.strip()
        if not name or name in seen:
            continue
        names.append(name)
        seen.add(name)
    return names


def _request_tuned_param_names(request):
    """Return parameter names listed in the active tuning request."""
    ranges = (request or {}).get("parameter_ranges", [])
    names = []
    seen = set()
    for spec in ranges or []:
        if not isinstance(spec, dict):
            continue
        name = str(spec.get("name", "")).strip()
        if name and name not in seen:
            names.append(name)
            seen.add(name)
    return names


def _varying_param_names(history):
    """Return parameters that actually vary across loaded samples."""
    param_names = _string_array_values(history, "param_names")
    params = np.asarray((history or {}).get("all_params", []), dtype=float)
    if params.ndim != 2 or params.shape[0] == 0:
        return param_names
    varying = []
    for idx, name in enumerate(param_names):
        if idx >= params.shape[1]:
            continue
        column = params[:, idx]
        finite = column[np.isfinite(column)]
        if finite.size > 1 and float(np.max(finite) - np.min(finite)) > 0.0:
            varying.append(name)
    return varying or param_names[: min(2, len(param_names))]


def _dropdown_options(values, all_label=None):
    """Return Dash dropdown options for string values."""
    options = []
    if all_label is not None:
        options.append({"label": all_label, "value": _LANDSCAPE_ALL})
    options.extend({"label": value, "value": value} for value in values)
    return options


def _first_valid(current_value, valid_values, fallback=None):
    """Return current_value if valid, otherwise a valid fallback."""
    valid = list(valid_values or [])
    if current_value in valid:
        return current_value
    if fallback in valid:
        return fallback
    return valid[0] if valid else None


def _landscape_metric_options(history):
    """Return metric selector options for loaded sample history."""
    metric_names = _string_array_values(history, "metric_names")
    options = [{"label": "Total selected loss", "value": "total_loss"}]
    options.extend({"label": f"Raw {name}", "value": f"raw:{name}"} for name in metric_names)
    return options


def _landscape_window_options(history):
    """Return time-window selector options from the observation axis."""
    obs_window = np.asarray((history or {}).get("obs_window", []), dtype=int)
    if obs_window.size == 0:
        return [{"label": "All windows", "value": _LANDSCAPE_ALL}]
    options = [{"label": "All windows", "value": _LANDSCAPE_ALL}]
    for window_idx in range(int(np.max(obs_window)) + 1):
        options.append({"label": f"Window {window_idx + 1}", "value": str(window_idx)})
    return options


def _aggregate_finite(values, mode):
    """Aggregate finite values with the requested reduction."""
    finite = np.asarray(values, dtype=float)
    finite = finite[_plot_metric_mask(finite)]
    if finite.size == 0:
        return float("nan")
    mode = mode if mode in _LANDSCAPE_AGGREGATION_NAMES else "mean"
    if mode == "median":
        return float(np.median(finite))
    if mode == "min":
        return float(np.min(finite))
    if mode == "p90":
        return float(np.percentile(finite, 90.0))
    if mode == "max":
        return float(np.max(finite))
    return float(np.mean(finite))


def _aggregate_sample_matrix(matrix, mode):
    """Aggregate a sample x observation matrix to one value per sample."""
    values = np.asarray(matrix, dtype=float)
    if values.ndim == 1:
        return values
    if values.size == 0:
        return np.full(values.shape[0], np.nan)
    return np.asarray([_aggregate_finite(row, mode) for row in values], dtype=float)


def _observation_mask(history, case_value, field_value, window_value):
    """Return an observation-axis mask for the selected case/field/window."""
    case_names = _string_array_values(history, "case_names")
    field_names = _string_array_values(history, "field_names")
    obs_case = np.asarray((history or {}).get("obs_case", []), dtype=int)
    obs_field = np.asarray((history or {}).get("obs_field", []), dtype=int)
    obs_window = np.asarray((history or {}).get("obs_window", []), dtype=int)
    mask = np.ones(obs_case.shape, dtype=bool)

    if case_value not in (None, _LANDSCAPE_ALL):
        if case_value not in case_names:
            return np.zeros(obs_case.shape, dtype=bool)
        mask &= obs_case == case_names.index(case_value)

    if field_value not in (None, _LANDSCAPE_ALL):
        if field_value not in field_names:
            return np.zeros(obs_case.shape, dtype=bool)
        mask &= obs_field == field_names.index(field_value)

    if window_value not in (None, _LANDSCAPE_ALL):
        try:
            window_idx = int(window_value)
        except (TypeError, ValueError):
            return np.zeros(obs_case.shape, dtype=bool)
        mask &= obs_window == window_idx

    return mask


def _request_weight(mapping, name):
    """Return a finite weight from a request mapping, defaulting to one."""
    if not isinstance(mapping, dict):
        return 1.0
    value = _finite_float(mapping.get(name, 1.0))
    return 1.0 if value is None else float(value)


def _total_selected_loss_series(history, request):
    """Compute scheduler-equivalent total selected loss for every history sample."""
    case_names = _string_array_values(history, "case_names")
    field_names = _string_array_values(history, "field_names")
    metric_names = _string_array_values(history, "metric_names")
    metric_index = {name: idx for idx, name in enumerate(metric_names)}
    missing_metrics = [name for name in LOSS_METRIC_NAMES if name not in metric_index]
    loss_metrics = np.asarray((history or {}).get("loss_metrics", []), dtype=float)
    if missing_metrics or loss_metrics.ndim != 3:
        return np.full(loss_metrics.shape[0] if loss_metrics.ndim else 0, np.nan)

    request = request if isinstance(request, dict) else {}
    loss_mode = request.get("loss_mode") or DEFAULT_LOSS_MODE
    aggregation_mode = request.get("aggregation_mode") or DEFAULT_AGGREGATION_MODE
    time_window_aggregation_mode = request.get("time_window_aggregation_mode") or aggregation_mode
    case_weights = request.get("case_weights", {})
    field_weights = request.get("field_weights", {})

    obs_by_axis = {}
    for obs_idx, (case_idx, field_idx, window_idx) in enumerate(
        zip(history["obs_case"], history["obs_field"], history["obs_window"])
    ):
        obs_by_axis[(int(case_idx), int(field_idx), int(window_idx))] = obs_idx

    case_window_counts = np.asarray(history.get("case_window_counts", []), dtype=int)
    output = np.empty(loss_metrics.shape[0], dtype=float)
    for sample_idx in range(loss_metrics.shape[0]):
        smart_losses = []
        smart_weights = []
        for case_idx, case_name in enumerate(case_names):
            if case_idx >= len(case_window_counts):
                continue
            window_count = int(case_window_counts[case_idx])
            case_weight = _request_weight(case_weights, case_name)
            for field_idx, field_name in enumerate(field_names):
                subwindow_losses = []
                for window_idx in range(window_count):
                    obs_idx = obs_by_axis.get((case_idx, field_idx, window_idx))
                    if obs_idx is None:
                        continue
                    metrics = {
                        metric_name: float(loss_metrics[sample_idx, obs_idx, metric_index[metric_name]])
                        for metric_name in LOSS_METRIC_NAMES
                    }
                    try:
                        diagnostics = compute_field_loss_diagnostics(metrics)
                        subwindow_losses.append(float(diagnostics["per_field_losses"][loss_mode]))
                    except Exception:
                        subwindow_losses.append(float("nan"))
                if not subwindow_losses:
                    continue
                window_aggregation = aggregate_losses(
                    subwindow_losses,
                    [1.0] * len(subwindow_losses),
                    time_window_aggregation_mode,
                )
                smart_losses.append(float(window_aggregation["loss"]))
                smart_weights.append(case_weight * _request_weight(field_weights, field_name))
        try:
            output[sample_idx] = float(aggregate_losses(smart_losses, smart_weights, aggregation_mode)["loss"])
        except Exception:
            output[sample_idx] = float("nan")
    return output


def _landscape_metric_series(history, request, metric_key, case_value, field_value, window_value, aggregation):
    """Return one plotted metric value per sample."""
    metric_key = metric_key or "total_loss"
    if metric_key == "total_loss":
        return _total_selected_loss_series(history, request)
    if not metric_key.startswith("raw:"):
        return np.full(np.asarray((history or {}).get("all_params", [])).shape[0], np.nan)

    metric_name = metric_key.split(":", 1)[1]
    metric_names = _string_array_values(history, "metric_names")
    if metric_name not in metric_names:
        return np.full(np.asarray((history or {}).get("all_params", [])).shape[0], np.nan)
    mask = _observation_mask(history, case_value, field_value, window_value)
    if not np.any(mask):
        return np.full(np.asarray((history or {}).get("all_params", [])).shape[0], np.nan)
    values = np.asarray(history["loss_metrics"], dtype=float)[:, mask, metric_names.index(metric_name)]
    return _aggregate_sample_matrix(values, aggregation)


def _landscape_metric_label(metric_key, request):
    """Return a readable label for a landscape metric selector value."""
    if metric_key == "total_loss":
        loss_mode = (request or {}).get("loss_mode") or DEFAULT_LOSS_MODE
        return f"total selected loss ({loss_mode})"
    if isinstance(metric_key, str) and metric_key.startswith("raw:"):
        return metric_key.split(":", 1)[1]
    return "metric"


def _robust_percentile_limits(values, lower=2.0, upper=98.0):
    """Return robust finite color limits, or ``(None, None)`` when unavailable."""
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size < 2:
        return None, None
    cmin, cmax = np.percentile(finite_values, [lower, upper])
    if np.isfinite(cmin) and np.isfinite(cmax) and cmin < cmax:
        return float(cmin), float(cmax)
    return None, None


def _symmetric_robust_limit(values, percentile=98.0):
    """Return a robust symmetric absolute color limit around zero."""
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size < 2:
        return None
    limit = float(np.percentile(np.abs(finite_values), percentile))
    if np.isfinite(limit) and limit > 0.0:
        return limit
    return None


def _fill_invalid_color_values(values, color_settings):
    """Keep plotted points visible even when their color transform is invalid."""
    color_values = np.asarray(values, dtype=float)
    if np.all(np.isfinite(color_values)):
        return color_values
    fallback = color_settings.get("cmin", 0.0)
    if fallback is None or not np.isfinite(float(fallback)):
        fallback = 0.0
    color_values = color_values.copy()
    color_values[~np.isfinite(color_values)] = float(fallback)
    return color_values


def _landscape_color_settings(metric_key, values, metric_label, *, fill_invalid=True):
    """Return transformed Plotly color values and settings for one metric."""
    raw_values = np.asarray(values, dtype=float)
    if metric_key == "raw:correlation":
        return {
            "values": raw_values,
            "colorscale": "RdBu",
            "cmin": -1.0,
            "cmax": 1.0,
            "title": metric_label,
        }

    if metric_key == "raw:std_ratio":
        with np.errstate(divide="ignore", invalid="ignore"):
            color_values = np.where(raw_values > 0.0, np.log(raw_values), np.nan)
        limit = _symmetric_robust_limit(color_values)
        settings = {
            "values": color_values,
            "colorscale": "RdBu",
            "title": f"log({metric_label})",
        }
        if limit is not None:
            settings.update(
                {
                    "cmin": -limit,
                    "cmax": limit,
                    "title": f"log({metric_label}) (centered 0, clipped p2-p98)",
                }
            )
        if fill_invalid:
            settings["values"] = _fill_invalid_color_values(settings["values"], settings)
        return settings

    if metric_key == "raw:bias_norm":
        color_values = raw_values
        limit = _symmetric_robust_limit(color_values)
        settings = {
            "values": color_values,
            "colorscale": "RdBu",
            "title": metric_label,
        }
        if limit is not None:
            settings.update(
                {
                    "cmin": -limit,
                    "cmax": limit,
                    "title": f"{metric_label} (centered 0, clipped p2-p98)",
                }
            )
        if fill_invalid:
            settings["values"] = _fill_invalid_color_values(settings["values"], settings)
        return settings

    if metric_key in _LANDSCAPE_LOSS_LIKE_METRICS:
        with np.errstate(invalid="ignore"):
            color_values = np.where(raw_values >= 0.0, np.log1p(raw_values), np.nan)
        cmin, cmax = _robust_percentile_limits(color_values)
        settings = {
            "values": color_values,
            "colorscale": "Viridis",
            "title": f"log1p({metric_label})",
        }
        if cmin is not None and cmax is not None:
            settings.update(
                {
                    "cmin": cmin,
                    "cmax": cmax,
                    "title": f"log1p({metric_label}) (clipped p2-p98)",
                }
            )
        if fill_invalid:
            settings["values"] = _fill_invalid_color_values(settings["values"], settings)
        return settings

    cmin, cmax = _robust_percentile_limits(raw_values)
    settings = {"values": raw_values, "colorscale": "Viridis", "title": metric_label}
    if cmin is not None and cmax is not None:
        settings.update({"cmin": cmin, "cmax": cmax, "title": f"{metric_label} (clipped p2-p98)"})
    if fill_invalid:
        settings["values"] = _fill_invalid_color_values(settings["values"], settings)
    return settings


def _landscape_hover_text(sample_id, batch_id, x_name, x_value, y_name, y_value, metric_label, metric_value):
    """Return hover text for one landscape sample."""
    return "<br>".join(
        [
            f"sample {int(sample_id)}",
            f"batch {int(batch_id)}",
            f"{x_name} {_fmt_metric(x_value, 6)}",
            f"{y_name} {_fmt_metric(y_value, 6)}",
            f"{metric_label} {_fmt_metric(metric_value, 6)}",
        ]
    )


def _empty_landscape_figure(title, message):
    """Return a stable empty landscape figure."""
    return _empty_diagnostics_figure(title, message)


def _binned_landscape(x_values, y_values, z_values, aggregation):
    """Return binned landscape heatmap arrays."""
    if x_values.size < 1 or y_values.size < 1:
        return None
    x_min = float(np.min(x_values))
    x_max = float(np.max(x_values))
    y_min = float(np.min(y_values))
    y_max = float(np.max(y_values))
    if x_min == x_max or y_min == y_max:
        return None

    x_edges = np.linspace(x_min, x_max, _LANDSCAPE_BINS + 1)
    y_edges = np.linspace(y_min, y_max, _LANDSCAPE_BINS + 1)
    x_bins = np.searchsorted(x_edges, x_values, side="right") - 1
    y_bins = np.searchsorted(y_edges, y_values, side="right") - 1
    x_bins = np.clip(x_bins, 0, _LANDSCAPE_BINS - 1)
    y_bins = np.clip(y_bins, 0, _LANDSCAPE_BINS - 1)

    cells = [[[] for _ in range(_LANDSCAPE_BINS)] for _ in range(_LANDSCAPE_BINS)]
    for x_bin, y_bin, z_value in zip(x_bins, y_bins, z_values):
        cells[int(y_bin)][int(x_bin)].append(float(z_value))

    grid = np.full((_LANDSCAPE_BINS, _LANDSCAPE_BINS), np.nan, dtype=float)
    for y_bin in range(_LANDSCAPE_BINS):
        for x_bin in range(_LANDSCAPE_BINS):
            grid[y_bin, x_bin] = _aggregate_finite(cells[y_bin][x_bin], aggregation)

    return {
        "x": 0.5 * (x_edges[:-1] + x_edges[1:]),
        "y": 0.5 * (y_edges[:-1] + y_edges[1:]),
        "z": grid,
    }


def build_landscape_figure(
    history,
    request,
    x_name,
    y_name,
    metric_key,
    aggregation,
    plot_mode,
    case_value,
    field_value,
    window_value,
    z_values=None,
):
    """Build the sample landscape scatter or binned heatmap."""
    if history is None:
        return _empty_landscape_figure(
            "Parameter Landscape",
            "Parameter landscapes will appear after sample-history chunks are available.",
        )
    param_names = _string_array_values(history, "param_names")
    if x_name not in param_names or y_name not in param_names:
        return _empty_landscape_figure("Parameter Landscape", "Choose two parameters to plot.")
    if x_name == y_name:
        return _empty_landscape_figure("Parameter Landscape", "Choose different X and Y parameters.")

    metric_label = _landscape_metric_label(metric_key, request)
    params = np.asarray(history["all_params"], dtype=float)
    x_values = params[:, param_names.index(x_name)]
    y_values = params[:, param_names.index(y_name)]
    if z_values is None:
        z_values = _landscape_metric_series(
            history,
            request,
            metric_key,
            case_value,
            field_value,
            window_value,
            aggregation,
        )
    else:
        z_values = np.asarray(z_values, dtype=float)
    sample_id = np.asarray(history.get("sample_id", np.arange(len(z_values))), dtype=int)
    batch_id = np.asarray(history.get("batch_id", np.zeros(len(z_values))), dtype=int)
    finite = np.isfinite(x_values) & np.isfinite(y_values) & _plot_metric_mask(z_values)
    if not np.any(finite):
        return _empty_landscape_figure("Parameter Landscape", "No plottable samples match the selected metric filters.")

    x_plot = x_values[finite]
    y_plot = y_values[finite]
    z_plot = z_values[finite]
    sample_plot = sample_id[finite]
    batch_plot = batch_id[finite]
    aggregation = aggregation if aggregation in _LANDSCAPE_AGGREGATION_NAMES else "mean"
    plot_mode = plot_mode if plot_mode in _LANDSCAPE_MODE_NAMES else "samples"
    color_settings = _landscape_color_settings(metric_key, z_plot, metric_label)
    color_plot = color_settings["values"]

    fig = go.Figure()
    if plot_mode == "binned":
        binned = _binned_landscape(x_plot, y_plot, z_plot, aggregation)
        if binned is None:
            return _empty_landscape_figure("Parameter Landscape", "Binned heatmap needs variation in both parameters.")
        binned_color_settings = _landscape_color_settings(
            metric_key,
            binned["z"],
            metric_label,
            fill_invalid=False,
        )
        fig.add_trace(
            go.Heatmap(
                x=binned["x"],
                y=binned["y"],
                z=binned_color_settings["values"],
                customdata=binned["z"],
                colorscale=color_settings["colorscale"],
                zmin=color_settings.get("cmin"),
                zmax=color_settings.get("cmax"),
                colorbar={"title": color_settings["title"]},
                hovertemplate=(
                    f"{x_name} bin %{{x:.6g}}<br>"
                    f"{y_name} bin %{{y:.6g}}<br>"
                    f"{aggregation} {metric_label} %{{customdata:.6g}}<extra></extra>"
                ),
            )
        )
        fig.add_trace(
            go.Scattergl(
                x=x_plot,
                y=y_plot,
                mode="markers",
                marker={"size": 5, "color": "rgba(15, 23, 42, 0.52)", "line": {"width": 0}},
                name="samples",
                hoverinfo="skip",
                showlegend=False,
            )
        )
    else:
        fig.add_trace(
            go.Scattergl(
                x=x_plot,
                y=y_plot,
                mode="markers",
                marker={
                    "size": 9,
                    "color": color_plot,
                    "colorscale": color_settings["colorscale"],
                    "cmin": color_settings.get("cmin"),
                    "cmax": color_settings.get("cmax"),
                    "colorbar": {"title": color_settings["title"]},
                    "line": {"color": "#0f172a", "width": 0.6},
                    "opacity": 0.86,
                },
                text=[
                    _landscape_hover_text(sample, batch, x_name, x_value, y_name, y_value, metric_label, z_value)
                    for sample, batch, x_value, y_value, z_value in zip(
                        sample_plot,
                        batch_plot,
                        x_plot,
                        y_plot,
                        z_plot,
                    )
                ],
                hovertemplate="%{text}<extra></extra>",
                name=metric_label,
            )
        )

    fig.update_layout(
        title={"text": "Parameter Landscape", "x": 0.02, "xanchor": "left"},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 58, "r": 18, "t": 48, "b": 62},
        height=430,
        hovermode="closest",
        uirevision="tune-parameter-landscape",
    )
    fig.update_xaxes(title=x_name, gridcolor="rgba(148, 163, 184, 0.24)")
    fig.update_yaxes(title=y_name, gridcolor="rgba(148, 163, 184, 0.24)")
    return fig


def build_parameter_correlation_figure(
    history,
    request,
    metric_key,
    aggregation,
    case_value,
    field_value,
    window_value,
    z_values=None,
):
    """Build a parameter-to-selected-metric correlation bar chart."""
    if history is None:
        return _empty_landscape_figure(
            "Parameter Correlations",
            "Parameter correlations will appear after sample-history chunks are available.",
        )
    metric_label = _landscape_metric_label(metric_key, request)
    if z_values is None:
        z_values = _landscape_metric_series(
            history,
            request,
            metric_key,
            case_value,
            field_value,
            window_value,
            aggregation,
        )
    else:
        z_values = np.asarray(z_values, dtype=float)
    params = np.asarray(history["all_params"], dtype=float)
    finite_z = _plot_metric_mask(z_values)
    if np.count_nonzero(finite_z) < 2 or float(np.std(z_values[finite_z])) == 0.0:
        return _empty_landscape_figure("Parameter Correlations", "At least two non-constant metric values are required.")

    param_names = _string_array_values(history, "param_names")
    tuned_names = _request_tuned_param_names(request)
    candidate_names = [name for name in tuned_names if name in param_names]
    if tuned_names and not candidate_names:
        return _empty_landscape_figure("Parameter Correlations", "No tuned parameters are present in sample history.")
    if not candidate_names:
        candidate_names = _varying_param_names(history)

    correlations = []
    for name in candidate_names:
        idx = param_names.index(name)
        if idx >= params.shape[1]:
            continue
        column = params[:, idx]
        finite = finite_z & np.isfinite(column)
        if np.count_nonzero(finite) < 2:
            continue
        column_values = column[finite]
        if float(np.std(column_values)) == 0.0:
            continue
        z_subset = z_values[finite]
        if float(np.std(z_subset)) == 0.0:
            continue
        correlation = float(
            np.mean((column_values - np.mean(column_values)) * (z_subset - np.mean(z_subset)))
            / (float(np.std(column_values)) * float(np.std(z_subset)))
        )
        if math.isfinite(correlation):
            correlations.append((name, correlation))

    if not correlations:
        return _empty_landscape_figure("Parameter Correlations", "No tuned parameters vary enough for correlation.")

    correlations = sorted(correlations, key=lambda item: item[1])
    names = [name for name, _value in correlations]
    values = [value for _name, value in correlations]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=values,
            y=names,
            orientation="h",
            marker={
                "color": values,
                "colorscale": "RdBu",
                "cmin": -1.0,
                "cmax": 1.0,
                "line": {"color": "#0f172a", "width": 0.4},
            },
            hovertemplate="parameter %{y}<br>correlation %{x:.4f}<extra></extra>",
        )
    )
    fig.update_layout(
        title={"text": "Parameter Correlations", "x": 0.02, "xanchor": "left"},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 98, "r": 18, "t": 48, "b": 54},
        height=430,
        showlegend=False,
        uirevision="tune-parameter-correlations",
        annotations=[
            {
                "text": metric_label,
                "xref": "paper",
                "yref": "paper",
                "x": 1.0,
                "y": 1.08,
                "showarrow": False,
                "font": {"size": 11, "color": "#475569"},
                "xanchor": "right",
            }
        ],
    )
    fig.update_xaxes(
        title="Pearson correlation",
        range=[-1.0, 1.0],
        zeroline=True,
        zerolinecolor="#0f172a",
        gridcolor="rgba(148, 163, 184, 0.24)",
    )
    fig.update_yaxes(gridcolor="rgba(148, 163, 184, 0.12)")
    return fig


def _field_bias_target(history, aggregation, case_value, field_value, window_value):
    """Return selected bias_norm values for field-response diagnostics."""
    if field_value in (None, _LANDSCAPE_ALL):
        return None
    metric_names = _string_array_values(history, "metric_names")
    if "bias_norm" not in metric_names:
        return None
    return _landscape_metric_series(
        history,
        {},
        "raw:bias_norm",
        case_value,
        field_value,
        window_value,
        aggregation,
    )


def _field_bias_selection_label(case_value, field_value, window_value):
    """Return a compact label for the selected field response."""
    case_label = "all cases" if case_value in (None, _LANDSCAPE_ALL) else str(case_value)
    window_label = "all windows"
    if window_value not in (None, _LANDSCAPE_ALL):
        try:
            window_label = f"window {int(window_value) + 1}"
        except (TypeError, ValueError):
            window_label = str(window_value)
    return f"{field_value}, {case_label}, {window_label}"


def _sensitivity_param_matrix(history, request):
    """Return candidate parameter names and values for sensitivity fits."""
    param_names = _string_array_values(history, "param_names")
    params = np.asarray(history.get("all_params", []), dtype=float)
    if params.ndim != 2 or params.shape[0] == 0:
        return [], np.empty((0, 0), dtype=float)

    tuned_names = _request_tuned_param_names(request)
    candidate_names = [name for name in tuned_names if name in param_names]
    if not candidate_names:
        candidate_names = _varying_param_names(history)

    columns = []
    names = []
    for name in candidate_names:
        idx = param_names.index(name)
        if idx >= params.shape[1]:
            continue
        columns.append(params[:, idx])
        names.append(name)
    if not columns:
        return [], np.empty((params.shape[0], 0), dtype=float)
    return names, np.column_stack(columns)


def _standardized_sensitivity_data(history, request, target_values):
    """Return finite standardized parameter data and centered target values."""
    names, values = _sensitivity_param_matrix(history, request)
    target = np.asarray(target_values, dtype=float)
    if not names or values.size == 0 or target.size == 0:
        return [], np.empty((0, 0), dtype=float), np.empty(0, dtype=float), 0
    if values.shape[0] != target.shape[0]:
        return [], np.empty((0, 0), dtype=float), np.empty(0, dtype=float), 0

    finite = _plot_metric_mask(target) & np.all(np.isfinite(values), axis=1)
    sample_count = int(np.count_nonzero(finite))
    if sample_count < 3:
        return [], np.empty((0, 0), dtype=float), np.empty(0, dtype=float), sample_count

    values = values[finite]
    target = target[finite]
    means = np.mean(values, axis=0)
    stds = np.std(values, axis=0)
    keep = np.isfinite(stds) & (stds > 0.0)
    if not np.any(keep):
        return [], np.empty((values.shape[0], 0), dtype=float), target - np.mean(target), sample_count

    kept_names = [name for name, selected in zip(names, keep) if bool(selected)]
    standardized = (values[:, keep] - means[keep]) / stds[keep]
    centered_target = target - np.mean(target)
    return kept_names, standardized, centered_target, sample_count


def _ridge_fit(design, target):
    """Return small-ridge least-squares coefficients for a centered design."""
    x = np.asarray(design, dtype=float)
    y = np.asarray(target, dtype=float)
    if x.ndim != 2 or x.shape[0] == 0 or x.shape[1] == 0:
        return np.empty(0, dtype=float)
    xtx = x.T @ x
    scale = float(np.trace(xtx) / max(1, xtx.shape[0]))
    alpha = max(scale, 1.0) * 1.0e-8
    try:
        return np.linalg.solve(xtx + alpha * np.eye(xtx.shape[0]), x.T @ y)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(x, y, rcond=None)[0]


def _r_squared(target, prediction):
    """Return R^2 for centered target values."""
    y = np.asarray(target, dtype=float)
    pred = np.asarray(prediction, dtype=float)
    denom = float(np.sum(y**2))
    if denom <= 0.0 or y.shape != pred.shape:
        return None
    return max(0.0, min(1.0, 1.0 - float(np.sum((y - pred) ** 2)) / denom))


def build_field_sensitivity_figure(history, request, aggregation, case_value, field_value, window_value):
    """Build main-effect sensitivity of selected field bias to tuned parameters."""
    if history is None:
        return _empty_landscape_figure(
            "Field Sensitivity",
            "Field sensitivity will appear after sample-history chunks are available.",
        )
    if field_value in (None, _LANDSCAPE_ALL):
        return _empty_landscape_figure("Field Sensitivity", "Choose one field to estimate field sensitivity.")

    target_values = _field_bias_target(history, aggregation, case_value, field_value, window_value)
    if target_values is None:
        return _empty_landscape_figure("Field Sensitivity", "bias_norm is not available for this sample history.")

    names, design, target, sample_count = _standardized_sensitivity_data(history, request, target_values)
    if sample_count < 3:
        return _empty_landscape_figure("Field Sensitivity", "At least three plottable samples are required.")
    if not names or design.shape[1] == 0:
        return _empty_landscape_figure("Field Sensitivity", "No tuned parameters vary enough for sensitivity fitting.")
    if float(np.std(target)) == 0.0:
        return _empty_landscape_figure("Field Sensitivity", "Selected field bias is constant across plottable samples.")

    coefficients = _ridge_fit(design, target)
    prediction = design @ coefficients
    r2 = _r_squared(target, prediction)
    order = np.argsort(coefficients)
    ordered_names = [names[idx] for idx in order]
    ordered_coefficients = [float(coefficients[idx]) for idx in order]
    limit = max([abs(value) for value in ordered_coefficients] or [1.0])
    if limit <= 0.0:
        limit = 1.0

    title = f"Field Sensitivity: {_field_bias_selection_label(case_value, field_value, window_value)}"
    subtitle = f"n={sample_count}"
    if r2 is not None:
        subtitle += f", linear R2={r2:.2f}"

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=ordered_coefficients,
            y=ordered_names,
            orientation="h",
            marker={
                "color": ordered_coefficients,
                "colorscale": "RdBu",
                "cmin": -limit,
                "cmax": limit,
                "line": {"color": "#0f172a", "width": 0.4},
            },
            hovertemplate="parameter %{y}<br>effect %{x:.4g} bias_norm / parameter sigma<extra></extra>",
        )
    )
    fig.update_layout(
        title={"text": title, "x": 0.02, "xanchor": "left"},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 98, "r": 18, "t": 48, "b": 62},
        height=430,
        showlegend=False,
        uirevision="tune-field-sensitivity",
        annotations=[
            {
                "text": subtitle,
                "xref": "paper",
                "yref": "paper",
                "x": 1.0,
                "y": 1.08,
                "showarrow": False,
                "font": {"size": 11, "color": "#475569"},
                "xanchor": "right",
            }
        ],
    )
    fig.update_xaxes(
        title="bias_norm change per parameter sigma",
        zeroline=True,
        zerolinecolor="#0f172a",
        gridcolor="rgba(148, 163, 184, 0.24)",
    )
    fig.update_yaxes(gridcolor="rgba(148, 163, 184, 0.12)")
    return fig


def build_field_interaction_figure(history, request, aggregation, case_value, field_value, window_value):
    """Build pairwise parameter-interaction sensitivity for selected field bias."""
    if history is None:
        return _empty_landscape_figure(
            "Field Interactions",
            "Field interactions will appear after sample-history chunks are available.",
        )
    if field_value in (None, _LANDSCAPE_ALL):
        return _empty_landscape_figure("Field Interactions", "Choose one field to estimate parameter interactions.")

    target_values = _field_bias_target(history, aggregation, case_value, field_value, window_value)
    if target_values is None:
        return _empty_landscape_figure("Field Interactions", "bias_norm is not available for this sample history.")

    names, design, target, sample_count = _standardized_sensitivity_data(history, request, target_values)
    if sample_count < 5:
        return _empty_landscape_figure("Field Interactions", "At least five plottable samples are required.")
    if len(names) < 2 or design.shape[1] < 2:
        return _empty_landscape_figure("Field Interactions", "At least two tuned parameters must vary.")
    if float(np.std(target)) == 0.0:
        return _empty_landscape_figure("Field Interactions", "Selected field bias is constant across plottable samples.")

    interaction = np.full((len(names), len(names)), np.nan, dtype=float)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pair_design = np.column_stack(
                [
                    design[:, i],
                    design[:, j],
                    design[:, i] * design[:, j],
                ]
            )
            coefficients = _ridge_fit(pair_design, target)
            if coefficients.size == 3 and math.isfinite(float(coefficients[2])):
                interaction[i, j] = float(coefficients[2])
                interaction[j, i] = float(coefficients[2])

    finite = interaction[np.isfinite(interaction)]
    if finite.size == 0:
        return _empty_landscape_figure("Field Interactions", "No finite pairwise interactions could be estimated.")
    limit = float(np.max(np.abs(finite)))
    if limit <= 0.0:
        limit = 1.0

    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            x=names,
            y=names,
            z=interaction,
            colorscale="RdBu",
            zmin=-limit,
            zmax=limit,
            colorbar={"title": "interaction"},
            hovertemplate="%{y} x %{x}<br>interaction %{z:.4g}<extra></extra>",
        )
    )
    fig.update_layout(
        title={
            "text": f"Field Interactions: {_field_bias_selection_label(case_value, field_value, window_value)}",
            "x": 0.02,
            "xanchor": "left",
        },
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#f8fafc",
        font={"color": "#0f172a"},
        margin={"l": 88, "r": 18, "t": 48, "b": 82},
        height=430,
        uirevision="tune-field-interactions",
    )
    fig.update_xaxes(tickangle=-35, gridcolor="rgba(148, 163, 184, 0.12)")
    fig.update_yaxes(gridcolor="rgba(148, 163, 184, 0.12)")
    return fig


def _diagnostics_signature(best_results, best_results_by_case=None, selected_groups=None):
    """Return a stable signature for the plotted top-result diagnostics."""
    try:
        return json.dumps(
            {
                "best_results": best_results or [],
                "best_results_by_case": best_results_by_case or {},
                "selected_groups": selected_groups or [],
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    except TypeError:
        return repr((best_results or [], best_results_by_case or {}, selected_groups or []))


def format_status_text(status):
    """Build the compact tuning-status summary text."""
    state = (status or {}).get("state", "idle")
    samples = int((status or {}).get("samples_evaluated", 0))
    total_samples = (status or {}).get("total_samples")
    elapsed_seconds = float((status or {}).get("elapsed_seconds", 0.0) or 0.0)
    best_total_loss = (status or {}).get("best_total_loss")
    if best_total_loss is None:
        best_text = "best smart loss: --"
    else:
        best_text = f"best smart loss: {float(best_total_loss):.6E}"
    uptime_text = f"uptime: {elapsed_seconds:.1f}s"
    if elapsed_seconds > 0.0:
        rate_text = f"samples/s: {samples / elapsed_seconds:.2f}"
    else:
        rate_text = "samples/s: --"
    worker_text = (
        f"active: {int((status or {}).get('active_evaluations', 0))} | "
        f"idle: {int((status or {}).get('idle_workers', 0))} | "
        f"initialized: {int((status or {}).get('initialized_workers', 0))} | "
        f"queued: {int((status or {}).get('queued_case_jobs', 0))}"
    )
    if total_samples is None:
        sample_text = f"samples: {samples}"
    else:
        sample_text = f"samples: {samples}/{int(total_samples)}"
    loss_mode = (status or {}).get("loss_mode") or DEFAULT_LOSS_MODE
    aggregation_mode = (status or {}).get("aggregation_mode") or DEFAULT_AGGREGATION_MODE
    case_window_counts = (status or {}).get("case_window_counts") or {}
    if case_window_counts:
        window_text = ", ".join(
            f"{case_name}:{count}" for case_name, count in sorted(case_window_counts.items())
        )
    else:
        window_text = "1"
    policy_text = f"loss: {loss_mode} | aggregation: {aggregation_mode} | windows: {window_text}"
    return f"state: {state} | {sample_text} | {uptime_text} | {rate_text} | {worker_text} | {policy_text} | {best_text}"


def loss_run_button_state(loss_runs, key, default_label="Run"):
    """Return display state for an ad-hoc loss-run button."""
    run_state = ((loss_runs or {}).get(str(key)) or {}).get("state")
    if run_state == "running":
        return "Running", action_button_style("#ca8a04"), True
    return default_label, action_button_style("#2563eb"), False


def build_results_table(top_results, selected_param_names):
    """Render the current best tuning results as an HTML table."""
    if not top_results:
        return build_results_placeholder()

    header = [
        html.Th("Rank"),
        html.Th("Smart Loss"),
        html.Th("Scaled RMSE Sum"),
    ]
    header.extend(html.Th(name) for name in selected_param_names)

    body_rows = []
    for row_index, result in enumerate(top_results, start=1):
        rank = result.get("rank", row_index)
        row = [
            html.Td(rank),
            html.Td(f"{float(result.get('total_loss', 0.0)):.6E}"),
            html.Td(
                ""
                if result.get("scaled_rmse_sum", result.get("simple_rms_sum")) is None
                else f"{float(result.get('scaled_rmse_sum', result.get('simple_rms_sum', 0.0))):.6E}"
            ),
        ]
        params = result.get("params", {})
        for name in selected_param_names:
            value = params.get(name)
            row.append(html.Td("" if value is None else f"{float(value):.6g}"))
        body_rows.append(html.Tr(row))

    return html.Table(
        [
            html.Thead(html.Tr(header)),
            html.Tbody(body_rows),
        ],
        className="tune-results-table",
    )


def mode_options_ready(
    strategy_mode,
    random_max_samples,
    resolve_spacing,
    simann_max_iters,
    simann_initial_temp,
    simann_final_temp,
):
    """Return whether the selected mode has the required options."""
    if strategy_mode == "random":
        if random_max_samples in (None, ""):
            return False
        try:
            value = float(random_max_samples)
        except (TypeError, ValueError):
            return False
        return int(value) == value and int(value) >= 1
    if strategy_mode == "resolve":
        if resolve_spacing in (None, ""):
            return False
        try:
            return float(resolve_spacing) > 0.0
        except (TypeError, ValueError):
            return False
    if strategy_mode == "simann":
        if simann_max_iters in (None, ""):
            return False
        try:
            max_iters_value = float(simann_max_iters)
        except (TypeError, ValueError):
            return False
        if int(max_iters_value) != max_iters_value or int(max_iters_value) < 1:
            return False
        try:
            return float(simann_initial_temp) > 0.0 and float(simann_final_temp) > 0.0
        except (TypeError, ValueError):
            return False
    return False


def case_window_setup_ready(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
    altitude_min_values,
    altitude_max_values,
):
    """Return whether every selected case has valid time and altitude window setup."""
    selected_cases = [
        value.strip()
        for value in (case_names or [])
        if isinstance(value, str) and value.strip()
    ]
    if not selected_cases:
        return False
    valid_rows = 0
    for raw_name, raw_start, raw_end, raw_average, raw_altitude_min, raw_altitude_max in zip(
        case_names or [],
        time_start_values or [],
        time_end_values or [],
        average_time_values or [],
        altitude_min_values or [],
        altitude_max_values or [],
    ):
        case_name = raw_name.strip() if isinstance(raw_name, str) else ""
        if not case_name:
            continue
        try:
            start_value = float(raw_start)
            end_value = float(raw_end)
            average_value = float(raw_average)
            altitude_min = float(raw_altitude_min)
            altitude_max = float(raw_altitude_max)
        except (TypeError, ValueError):
            return False
        if int(start_value) != start_value or int(end_value) != end_value or int(average_value) != average_value:
            return False
        if not math.isfinite(altitude_min) or not math.isfinite(altitude_max) or altitude_max < altitude_min:
            return False
        start = int(start_value)
        end = int(end_value)
        average = int(average_value)
        if end <= start or average < 1 or (end - start) % average != 0:
            return False
        valid_rows += 1
    return valid_rows == len(selected_cases)


def resolve_total_samples_text(resolve_spacing, param_names, min_values, max_values):
    """Estimate resolve grid size from current ranges without materializing points."""
    if resolve_spacing in (None, ""):
        return ""
    try:
        spacing = float(resolve_spacing)
    except (TypeError, ValueError):
        return ""
    if spacing <= 0.0:
        return ""

    total = 1
    complete_ranges = 0
    for param_name, min_text, max_text in zip(param_names or [], min_values or [], max_values or []):
        param_name = "" if param_name is None else str(param_name).strip()
        min_text = "" if min_text is None else str(min_text).strip()
        max_text = "" if max_text is None else str(max_text).strip()
        if not all((param_name, min_text, max_text)):
            continue
        try:
            min_value = float(min_text.replace("D", "E").replace("d", "e"))
            max_value = float(max_text.replace("D", "E").replace("d", "e"))
        except ValueError:
            return ""
        if min_value > max_value:
            return ""
        span = max_value - min_value
        count = 1 if span == 0.0 else int(math.ceil(span / spacing)) + 1
        total *= count
        complete_ranges += 1

    if complete_ranges == 0:
        return ""
    return f"Total samples: {total:,}"


def register_display_callbacks(app):
    """Register result-table and status-display callbacks."""

    @app.callback(
        Output("tune-landscape-x-param", "options"),
        Output("tune-landscape-x-param", "value"),
        Output("tune-landscape-y-param", "options"),
        Output("tune-landscape-y-param", "value"),
        Output("tune-landscape-metric", "options"),
        Output("tune-landscape-metric", "value"),
        Output("tune-landscape-case", "options"),
        Output("tune-landscape-case", "value"),
        Output("tune-landscape-field", "options"),
        Output("tune-landscape-field", "value"),
        Output("tune-landscape-window", "options"),
        Output("tune-landscape-window", "value"),
        Input("tune-status", "data"),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        State("tune-landscape-x-param", "value"),
        State("tune-landscape-y-param", "value"),
        State("tune-landscape-metric", "value"),
        State("tune-landscape-case", "value"),
        State("tune-landscape-field", "value"),
        State("tune-landscape-window", "value"),
    )
    def sync_landscape_controls(
        status,
        selected_param_names,
        current_x,
        current_y,
        current_metric,
        current_case,
        current_field,
        current_window,
    ):
        """Keep landscape controls matched to available sample-history axes."""
        job_dir = _job_dir_from_status(status)
        history = None
        if _sample_history_signature(job_dir):
            try:
                history = _load_sample_history(job_dir)
            except Exception:
                history = None

        if history is None:
            param_names = _clean_param_names(selected_param_names)
            metric_options = _landscape_metric_options({"metric_names": list(LOSS_METRIC_NAMES)})
            case_names = sorted(str(name) for name in ((status or {}).get("case_window_counts") or {}))
            case_options = _dropdown_options(case_names, "All cases")
            field_options = _dropdown_options([], "All fields")
            max_window = max([int(value) for value in ((status or {}).get("case_window_counts") or {}).values()] or [1])
            window_options = [{"label": "All windows", "value": _LANDSCAPE_ALL}]
            window_options.extend({"label": f"Window {idx + 1}", "value": str(idx)} for idx in range(max_window))
        else:
            history_param_names = _string_array_values(history, "param_names")
            preferred = [
                name
                for name in _clean_param_names(selected_param_names)
                if name in history_param_names and name in _varying_param_names(history)
            ]
            preferred.extend(name for name in _varying_param_names(history) if name not in preferred)
            param_names = preferred
            metric_options = _landscape_metric_options(history)
            case_options = _dropdown_options(_string_array_values(history, "case_names"), "All cases")
            field_options = _dropdown_options(_string_array_values(history, "field_names"), "All fields")
            window_options = _landscape_window_options(history)

        param_options = _dropdown_options(param_names)
        param_values = [option["value"] for option in param_options]
        x_value = _first_valid(current_x, param_values, param_names[0] if param_names else None)
        y_fallback = param_names[1] if len(param_names) > 1 else (param_names[0] if param_names else None)
        y_value = _first_valid(current_y, param_values, y_fallback)
        if x_value == y_value and len(param_names) > 1:
            y_value = next((name for name in param_names if name != x_value), y_value)

        metric_values = [option["value"] for option in metric_options]
        case_values = [option["value"] for option in case_options]
        field_values = [option["value"] for option in field_options]
        window_values = [option["value"] for option in window_options]
        return (
            param_options,
            x_value,
            param_options,
            y_value,
            metric_options,
            _first_valid(current_metric, metric_values, "total_loss"),
            case_options,
            _first_valid(current_case, case_values, _LANDSCAPE_ALL),
            field_options,
            _first_valid(current_field, field_values, _LANDSCAPE_ALL),
            window_options,
            _first_valid(current_window, window_values, _LANDSCAPE_ALL),
        )

    @app.callback(
        Output("tune-landscape-plot", "figure"),
        Output("tune-parameter-correlation-plot", "figure"),
        Output("tune-field-sensitivity-plot", "figure"),
        Output("tune-field-interaction-plot", "figure"),
        Output("tune-landscape-signature", "data"),
        Input("tune-status", "data"),
        Input("tune-landscape-x-param", "value"),
        Input("tune-landscape-y-param", "value"),
        Input("tune-landscape-metric", "value"),
        Input("tune-landscape-aggregation", "value"),
        Input("tune-landscape-mode", "value"),
        Input("tune-landscape-case", "value"),
        Input("tune-landscape-field", "value"),
        Input("tune-landscape-window", "value"),
        State("tune-landscape-signature", "data"),
    )
    def render_landscape(
        status,
        x_name,
        y_name,
        metric_key,
        aggregation,
        plot_mode,
        case_value,
        field_value,
        window_value,
        previous_signature,
    ):
        """Render all-sample parameter landscape diagnostics."""
        job_dir = _job_dir_from_status(status)
        history_signature = _sample_history_signature(job_dir)
        signature = json.dumps(
            {
                "job_dir": str(job_dir) if job_dir is not None else None,
                "history": history_signature,
                "x": x_name,
                "y": y_name,
                "metric": metric_key,
                "aggregation": aggregation,
                "plot_mode": plot_mode,
                "case": case_value,
                "field": field_value,
                "window": window_value,
                "loss_mode": (status or {}).get("loss_mode"),
                "aggregation_mode": (status or {}).get("aggregation_mode"),
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        if signature == (previous_signature or ""):
            return no_update, no_update, no_update, no_update, no_update
        if not history_signature:
            return (
                _empty_landscape_figure(
                    "Parameter Landscape",
                    "Parameter landscapes will appear after sample-history chunks are available.",
                ),
                _empty_landscape_figure(
                    "Parameter Correlations",
                    "Parameter correlations will appear after sample-history chunks are available.",
                ),
                _empty_landscape_figure(
                    "Field Sensitivity",
                    "Field sensitivity will appear after sample-history chunks are available.",
                ),
                _empty_landscape_figure(
                    "Field Interactions",
                    "Field interactions will appear after sample-history chunks are available.",
                ),
                signature,
            )

        try:
            history = _load_sample_history(job_dir)
            request = _read_job_request(job_dir)
            if (status or {}).get("loss_mode"):
                request["loss_mode"] = (status or {}).get("loss_mode")
            if (status or {}).get("aggregation_mode"):
                request["aggregation_mode"] = (status or {}).get("aggregation_mode")
            if (status or {}).get("time_window_aggregation_mode"):
                request["time_window_aggregation_mode"] = (status or {}).get("time_window_aggregation_mode")
        except Exception as exc:
            return (
                _empty_landscape_figure("Parameter Landscape", f"Unable to load sample history: {exc}"),
                _empty_landscape_figure("Parameter Correlations", f"Unable to load sample history: {exc}"),
                _empty_landscape_figure("Field Sensitivity", f"Unable to load sample history: {exc}"),
                _empty_landscape_figure("Field Interactions", f"Unable to load sample history: {exc}"),
                signature,
            )

        z_values = _landscape_metric_series(
            history,
            request,
            metric_key,
            case_value,
            field_value,
            window_value,
            aggregation,
        )
        return (
            build_landscape_figure(
                history,
                request,
                x_name,
                y_name,
                metric_key,
                aggregation,
                plot_mode,
                case_value,
                field_value,
                window_value,
                z_values=z_values,
            ),
            build_parameter_correlation_figure(
                history,
                request,
                metric_key,
                aggregation,
                case_value,
                field_value,
                window_value,
                z_values=z_values,
            ),
            build_field_sensitivity_figure(
                history,
                request,
                aggregation,
                case_value,
                field_value,
                window_value,
            ),
            build_field_interaction_figure(
                history,
                request,
                aggregation,
                case_value,
                field_value,
                window_value,
            ),
            signature,
        )

    @app.callback(
        Output("tune-taylor-diagram", "figure"),
        Output("tune-parameter-box-plot", "figure"),
        Output("tune-diagnostics-signature", "data"),
        Input("tune-best-results", "data"),
        Input("tune-best-results-by-case", "data"),
        Input("tune-parameter-box-groups", "value"),
        State("tune-diagnostics-signature", "data"),
    )
    def render_diagnostics(best_results, best_results_by_case, selected_box_groups, previous_signature):
        """Update diagnostics figures without rebuilding the graph components."""
        best_results = best_results or []
        best_results_by_case = best_results_by_case or {}
        signature = _diagnostics_signature(best_results, best_results_by_case, selected_box_groups)
        if signature == (previous_signature or ""):
            return no_update, no_update, no_update
        return build_taylor_figure(best_results), build_parameter_box_figure(best_results, best_results_by_case, selected_box_groups), signature

    @app.callback(
        Output("tune-parameter-box-groups", "options"),
        Output("tune-parameter-box-groups", "value"),
        Input("tune-best-results-by-case", "data"),
        State("tune-parameter-box-groups", "value"),
    )
    def sync_parameter_box_group_options(best_results_by_case, current_value):
        """Expose aggregate and per-case top-result lists as plot toggles."""
        return parameter_box_group_options(best_results_by_case, current_value)

    @app.callback(
        Output("tune-results-summary", "children"),
        Output("tune-results-container", "children"),
        Output("tune-start-button", "disabled"),
        Output("tune-stop-button", "disabled"),
        Output("tune-start-button", "style"),
        Output("tune-stop-button", "style"),
        Output({"type": "tune-loss-run-button", "action": "window"}, "children"),
        Output({"type": "tune-loss-run-button", "action": "window"}, "style"),
        Output({"type": "tune-loss-run-button", "action": "window"}, "disabled"),
        Output({"type": "tune-loss-run-button", "action": "complete"}, "children"),
        Output({"type": "tune-loss-run-button", "action": "complete"}, "style"),
        Output({"type": "tune-loss-run-button", "action": "complete"}, "disabled"),
        Output("tune-mode-random", "disabled"),
        Output("tune-mode-resolve", "disabled"),
        Output("tune-mode-simann", "disabled"),
        Output("tune-mode-random", "style"),
        Output("tune-mode-resolve", "style"),
        Output("tune-mode-simann", "style"),
        Output("tune-loss-mode-scaled-rmse", "disabled"),
        Output("tune-loss-mode-centered-rmse-bias", "disabled"),
        Output("tune-loss-mode-taylor-components", "disabled"),
        Output("tune-loss-mode-taylor-components-squared", "disabled"),
        Output("tune-loss-mode-shape-first", "disabled"),
        Output("tune-loss-mode-bias-light-taylor", "disabled"),
        Output("tune-loss-mode-decomposed-taylor", "disabled"),
        Output("tune-aggregation-mean-max", "disabled"),
        Output("tune-aggregation-mean-worst-quantile", "disabled"),
        Output("tune-loss-mode-scaled-rmse", "style"),
        Output("tune-loss-mode-centered-rmse-bias", "style"),
        Output("tune-loss-mode-taylor-components", "style"),
        Output("tune-loss-mode-taylor-components-squared", "style"),
        Output("tune-loss-mode-shape-first", "style"),
        Output("tune-loss-mode-bias-light-taylor", "style"),
        Output("tune-loss-mode-decomposed-taylor", "style"),
        Output("tune-aggregation-mean-max", "style"),
        Output("tune-aggregation-mean-worst-quantile", "style"),
        Output("tune-random-options", "style"),
        Output("tune-resolve-options", "style"),
        Output("tune-simann-options", "style"),
        Output("tune-no-mode-options", "style"),
        Output("tune-resolve-total-samples", "children"),
        Output("tune-random-max-samples", "disabled"),
        Output("tune-resolve-spacing", "disabled"),
        Output("tune-simann-max-iters", "disabled"),
        Output("tune-simann-initial-temp", "disabled"),
        Output("tune-simann-final-temp", "disabled"),
        Output("tune-case-add", "disabled"),
        Output({"type": "tune-case-name", "index": ALL}, "disabled"),
        Output({"type": "tune-case-time-start", "index": ALL}, "disabled"),
        Output({"type": "tune-case-time-end", "index": ALL}, "disabled"),
        Output({"type": "tune-case-average-time", "index": ALL}, "disabled"),
        Output({"type": "tune-case-altitude-min", "index": ALL}, "disabled"),
        Output({"type": "tune-case-altitude-max", "index": ALL}, "disabled"),
        Output({"type": "tune-case-remove", "index": ALL}, "disabled"),
        Output("tune-field-selector", "disabled"),
        Output("tune-batch-size", "disabled"),
        Output("tune-max-workers", "disabled"),
        Output("tune-range-add", "disabled"),
        Output({"type": "tune-range-param", "index": ALL}, "disabled"),
        Output({"type": "tune-range-min", "index": ALL}, "disabled"),
        Output({"type": "tune-range-max", "index": ALL}, "disabled"),
        Output({"type": "tune-range-remove", "index": ALL}, "disabled"),
        Output({"type": "tune-config-button", "name": ALL}, "disabled"),
        Output({"type": "tune-config-button", "name": ALL}, "style"),
        Input("tune-status", "data"),
        Input("tune-top-results", "data"),
        Input("tune-loss-runs", "data"),
        Input("tune-active-job", "data"),
        Input("tune-strategy-mode", "data"),
        Input("tune-loss-mode", "data"),
        Input("tune-aggregation-mode", "data"),
        Input("tune-selected-config", "data"),
        Input("tune-tunable-configs", "data"),
        Input("tune-random-max-samples", "value"),
        Input("tune-resolve-spacing", "value"),
        Input("tune-simann-max-iters", "value"),
        Input("tune-simann-initial-temp", "value"),
        Input("tune-simann-final-temp", "value"),
        Input("tune-tunable-names", "data"),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        Input({"type": "tune-case-time-start", "index": ALL}, "value"),
        Input({"type": "tune-case-time-end", "index": ALL}, "value"),
        Input({"type": "tune-case-average-time", "index": ALL}, "value"),
        Input({"type": "tune-case-altitude-min", "index": ALL}, "value"),
        Input({"type": "tune-case-altitude-max", "index": ALL}, "value"),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        Input({"type": "tune-range-min", "index": ALL}, "value"),
        Input({"type": "tune-range-max", "index": ALL}, "value"),
    )
    def render_tuning_state(
        status,
        top_results,
        loss_runs,
        active_job,
        strategy_mode,
        loss_mode,
        aggregation_mode,
        selected_config,
        tunable_configs,
        random_max_samples,
        resolve_spacing,
        simann_max_iters,
        simann_initial_temp,
        simann_final_temp,
        tunable_names,
        selected_case_names,
        time_start_values,
        time_end_values,
        average_time_values,
        altitude_min_values,
        altitude_max_values,
        selected_param_names,
        min_values,
        max_values,
    ):
        """Render the current tuning status and best-results table."""
        param_names = [
            value.strip()
            for value in (selected_param_names or [])
            if isinstance(value, str) and value.strip()
        ]
        status_text = format_status_text(status)
        results_table = build_results_table(top_results or [], param_names)
        running = bool(active_job)
        case_ready = case_window_setup_ready(
            selected_case_names,
            time_start_values,
            time_end_values,
            average_time_values,
            altitude_min_values,
            altitude_max_values,
        )
        start_disabled = (
            running
            or not mode_options_ready(
                strategy_mode,
                random_max_samples,
                resolve_spacing,
                simann_max_iters,
                simann_initial_temp,
                simann_final_temp,
            )
            or not case_ready
        )
        has_results = bool(top_results)
        selected_cases = [
            value.strip()
            for value in (selected_case_names or [])
            if isinstance(value, str) and value.strip()
        ]
        window_label, window_style, window_disabled = loss_run_button_state(
            loss_runs,
            "window",
            default_label="Run loss window",
        )
        complete_label, complete_style, complete_disabled = loss_run_button_state(
            loss_runs,
            "complete",
            default_label="Run complete",
        )
        window_disabled = window_disabled or not has_results or not case_ready or not param_names
        complete_disabled = complete_disabled or not has_results or not selected_cases
        if window_disabled and window_label != "Running":
            window_style = action_button_style("#2563eb", disabled=True)
        if complete_disabled and complete_label != "Running":
            complete_style = action_button_style("#2563eb", disabled=True)
        selected_count = len(set(param_names))
        add_disabled = running or selected_count >= len(tunable_names or [])
        case_disabled = [running] * len(selected_case_names or [])
        range_disabled = [running] * len(selected_param_names or [])
        config_values = [
            str(config.get("value", "")).strip()
            for config in (tunable_configs or [])
            if str(config.get("value", "")).strip()
        ]
        config_disabled = [running] * len(config_values)
        config_styles = [
            config_button_style(selected=value == selected_config, disabled=running)
            for value in config_values
        ]
        resolve_total_text = resolve_total_samples_text(
            resolve_spacing,
            selected_param_names,
            min_values,
            max_values,
        )
        return (
            status_text,
            results_table,
            start_disabled,
            not running,
            action_button_style("#16a34a", disabled=start_disabled),
            action_button_style("#b91c1c", disabled=not running),
            window_label,
            window_style,
            window_disabled,
            complete_label,
            complete_style,
            complete_disabled,
            running,
            running,
            running,
            mode_button_style(selected=strategy_mode == "random", disabled=running),
            mode_button_style(selected=strategy_mode == "resolve", disabled=running),
            mode_button_style(selected=strategy_mode == "simann", disabled=running),
            running,
            running,
            running,
            running,
            running,
            running,
            running,
            running,
            running,
            mode_button_style(selected=loss_mode == "scaled_rmse", disabled=running),
            mode_button_style(selected=(loss_mode or DEFAULT_LOSS_MODE) == "centered_rmse_bias", disabled=running),
            mode_button_style(selected=loss_mode == "taylor_components", disabled=running),
            mode_button_style(selected=loss_mode == "taylor_components_squared", disabled=running),
            mode_button_style(selected=loss_mode == "shape_first", disabled=running),
            mode_button_style(selected=loss_mode == "bias_light_taylor", disabled=running),
            mode_button_style(selected=loss_mode == "decomposed_taylor", disabled=running),
            mode_button_style(selected=(aggregation_mode or DEFAULT_AGGREGATION_MODE) == "mean_max", disabled=running),
            mode_button_style(selected=aggregation_mode == "mean_worst_quantile", disabled=running),
            mode_options_block_style(strategy_mode == "random"),
            mode_options_block_style(strategy_mode == "resolve"),
            mode_options_block_style(strategy_mode == "simann"),
            mode_options_block_style(strategy_mode not in {"random", "resolve", "simann"}),
            resolve_total_text,
            running,
            running,
            running,
            running,
            running,
            running,
            case_disabled,
            case_disabled,
            case_disabled,
            case_disabled,
            case_disabled,
            case_disabled,
            case_disabled,
            running,
            running,
            running,
            add_disabled,
            range_disabled,
            range_disabled,
            range_disabled,
            range_disabled,
            config_disabled,
            config_styles,
        )

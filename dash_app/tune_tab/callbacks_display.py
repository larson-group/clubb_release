"""Callbacks for rendering tuning status and result tables."""

from __future__ import annotations

import json
import math

from dash import ALL, Input, Output, State, dcc, html, no_update
import plotly.graph_objects as go

from .layout import action_button_style, build_results_placeholder, mode_button_style, mode_options_block_style
from tuner.taylor_metrics import DEFAULT_AGGREGATION_MODE, DEFAULT_LOSS_MODE


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


def mode_options_ready(strategy_mode, random_max_samples, resolve_spacing):
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
    return False


def case_window_setup_ready(case_names, time_start_values, time_end_values, average_time_values):
    """Return whether every selected case has an integral average-time window setup."""
    selected_cases = [
        value.strip()
        for value in (case_names or [])
        if isinstance(value, str) and value.strip()
    ]
    if not selected_cases:
        return False
    valid_rows = 0
    for raw_name, raw_start, raw_end, raw_average in zip(
        case_names or [],
        time_start_values or [],
        time_end_values or [],
        average_time_values or [],
    ):
        case_name = raw_name.strip() if isinstance(raw_name, str) else ""
        if not case_name:
            continue
        try:
            start_value = float(raw_start)
            end_value = float(raw_end)
            average_value = float(raw_average)
        except (TypeError, ValueError):
            return False
        if int(start_value) != start_value or int(end_value) != end_value or int(average_value) != average_value:
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
        Output("tune-mode-random", "style"),
        Output("tune-mode-resolve", "style"),
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
        Output("tune-no-mode-options", "style"),
        Output("tune-resolve-total-samples", "children"),
        Output("tune-random-max-samples", "disabled"),
        Output("tune-resolve-spacing", "disabled"),
        Output("tune-case-add", "disabled"),
        Output({"type": "tune-case-name", "index": ALL}, "disabled"),
        Output({"type": "tune-case-time-start", "index": ALL}, "disabled"),
        Output({"type": "tune-case-time-end", "index": ALL}, "disabled"),
        Output({"type": "tune-case-average-time", "index": ALL}, "disabled"),
        Output({"type": "tune-case-remove", "index": ALL}, "disabled"),
        Output("tune-field-selector", "disabled"),
        Output("tune-batch-size", "disabled"),
        Output("tune-max-workers", "disabled"),
        Output("tune-range-add", "disabled"),
        Output({"type": "tune-range-param", "index": ALL}, "disabled"),
        Output({"type": "tune-range-min", "index": ALL}, "disabled"),
        Output({"type": "tune-range-max", "index": ALL}, "disabled"),
        Output({"type": "tune-range-remove", "index": ALL}, "disabled"),
        Input("tune-status", "data"),
        Input("tune-top-results", "data"),
        Input("tune-loss-runs", "data"),
        Input("tune-active-job", "data"),
        Input("tune-strategy-mode", "data"),
        Input("tune-loss-mode", "data"),
        Input("tune-aggregation-mode", "data"),
        Input("tune-random-max-samples", "value"),
        Input("tune-resolve-spacing", "value"),
        Input("tune-tunable-names", "data"),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        Input({"type": "tune-case-time-start", "index": ALL}, "value"),
        Input({"type": "tune-case-time-end", "index": ALL}, "value"),
        Input({"type": "tune-case-average-time", "index": ALL}, "value"),
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
        random_max_samples,
        resolve_spacing,
        tunable_names,
        selected_case_names,
        time_start_values,
        time_end_values,
        average_time_values,
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
        )
        start_disabled = (
            running
            or not mode_options_ready(strategy_mode, random_max_samples, resolve_spacing)
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
            mode_button_style(selected=strategy_mode == "random", disabled=running),
            mode_button_style(selected=strategy_mode == "resolve", disabled=running),
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
            mode_options_block_style(strategy_mode not in {"random", "resolve"}),
            resolve_total_text,
            running,
            running,
            running,
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
        )

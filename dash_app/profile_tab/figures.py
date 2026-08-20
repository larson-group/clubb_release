"""Plotly figures for stored and actively running timing profiles."""

from __future__ import annotations

import statistics
import zlib
from collections import defaultdict
from typing import Any, Iterable, Sequence

import plotly.graph_objects as go
from plotly.colors import qualitative

from utilities.timing_profiles import GROUP_WALL_TIMER


METRICS = {
    "throughput_columns_per_second": "Throughput (runs / second)",
    "throughput_column_steps_per_second": "Throughput (iterations / second)",
    "throughput_grid_box_iterations_per_second": "Throughput (grid box iterations / second)",
    "timer_max_seconds": "Runtime (seconds)",
    "process_imbalance_ratio": "Process imbalance (max / mean)",
}

X_AXES = {
    "total_columns": "Overall batch size",
    "columns_per_process": "Batch size per process",
}

_DASHES = ("solid", "dash", "dot", "dashdot", "longdash", "longdashdot")


def _series_color(identifier: str) -> str:
    index = zlib.crc32(str(identifier).encode("utf-8")) % len(qualitative.Plotly)
    return qualitative.Plotly[index]


def _theme(theme_name: str) -> tuple[str, str, str]:
    if str(theme_name or "dark").lower() == "light":
        return "plotly_white", "#0f172a", "#cbd5e1"
    return "plotly_dark", "#e5e7eb", "#334155"


def _number(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed


def _integer(value: Any) -> int | None:
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return None


def _profile_id(row: dict[str, Any]) -> str:
    return str(row.get("profile_id") or row.get("run_id") or "active")


def _profile_label(row: dict[str, Any]) -> str:
    return str(row.get("profile_label") or row.get("label") or _profile_id(row))


def _apply_layout(
    figure: go.Figure,
    *,
    theme_name: str,
    title: str,
    x_title: str,
    y_title: str,
    x_scale: str = "linear",
) -> go.Figure:
    template, color, grid = _theme(theme_name)
    figure.update_layout(
        template=template,
        title={"text": title, "x": 0.035, "y": 0.96, "font": {"size": 15}},
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"color": color, "size": 12},
        hovermode="closest",
        hoverlabel={"font": {"size": 12}},
        legend={
            "orientation": "h",
            "yanchor": "top",
            "y": -0.2,
            "xanchor": "left",
            "x": 0,
            "title": {"text": ""},
            "font": {"size": 11},
            "traceorder": "normal",
        },
        margin={"l": 68, "r": 24, "t": 64, "b": 108},
        xaxis_title=x_title,
        yaxis_title=y_title,
    )
    figure.update_xaxes(
        gridcolor=grid,
        gridwidth=1,
        linecolor=grid,
        tickfont={"size": 11},
        title_font={"size": 12},
        type="log" if x_scale == "log" else "linear",
        zerolinecolor=grid,
    )
    figure.update_yaxes(
        gridcolor=grid,
        gridwidth=1,
        linecolor=grid,
        tickfont={"size": 11},
        title_font={"size": 12},
        zerolinecolor=grid,
    )
    return figure


def empty_profile_figure(message: str, theme_name: str = "dark") -> go.Figure:
    template, color, grid = _theme(theme_name)
    figure = go.Figure()
    figure.update_layout(
        template=template,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"color": color},
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[
            {
                "text": message,
                "xref": "paper",
                "yref": "paper",
                "x": 0.5,
                "y": 0.5,
                "showarrow": False,
                "font": {"size": 14, "color": color},
                "bgcolor": "rgba(100,116,139,0.08)",
                "bordercolor": grid,
                "borderpad": 12,
            }
        ],
        margin={"l": 32, "r": 24, "t": 50, "b": 36},
    )
    figure.update_xaxes(gridcolor=grid)
    figure.update_yaxes(gridcolor=grid)
    return figure


def timer_options(rows: Iterable[dict[str, object]], preferred_timer: str = "") -> list[dict[str, str]]:
    names = sorted(
        {
            str(row.get("timer_name") or "").strip()
            for row in rows
            if str(row.get("timer_name") or "").strip()
        }
    )
    preferred = str(preferred_timer or "").strip()
    if preferred in names:
        names.remove(preferred)
        names.insert(0, preferred)
    elif preferred:
        names.insert(0, preferred)
    return [
        {
            "label": "Process-group wall time" if name == GROUP_WALL_TIMER else name,
            "value": name,
        }
        for name in names
    ]


def profile_figure(
    rows: Iterable[dict[str, object]],
    timer_name: str,
    metric: str,
    theme_name: str = "dark",
    *,
    x_axis: str = "total_columns",
    x_scale: str = "linear",
) -> go.Figure:
    """Plot repetition means for one timer across one or more profiles."""
    rows = list(rows or [])
    selected_timer = str(timer_name or "").strip()
    if not rows:
        return empty_profile_figure("Select or run a profile to see timing results.", theme_name)
    if metric not in METRICS:
        metric = "throughput_columns_per_second"
    if x_axis not in X_AXES:
        x_axis = "total_columns"

    grouped: dict[tuple[str, str, int, int, int], list[float]] = defaultdict(list)
    for row in rows:
        if str(row.get("timer_name") or "") != selected_timer:
            continue
        if str(row.get("status") or "") != "success":
            continue
        process_count = _integer(row.get("process_count"))
        total_columns = _integer(row.get("total_columns"))
        columns_per_process = _integer(row.get("columns_per_process"))
        value = _number(row.get(metric))
        if None in (process_count, total_columns, columns_per_process, value):
            continue
        x_value = total_columns if x_axis == "total_columns" else columns_per_process
        grouped[
            (
                _profile_id(row),
                _profile_label(row),
                int(process_count),
                int(x_value),
                int(columns_per_process),
            )
        ].append(float(value))

    if not grouped:
        return empty_profile_figure(f"No successful {selected_timer or 'timer'} results yet.", theme_name)

    figure = go.Figure()
    profiles = list(dict.fromkeys((key[0], key[1]) for key in grouped))
    colors = {profile_id: _series_color(profile_id) for profile_id, _label in profiles}
    process_counts = sorted({key[2] for key in grouped})
    dashes = {count: _DASHES[index % len(_DASHES)] for index, count in enumerate(process_counts)}
    multiple_profiles = len(profiles) > 1
    for profile_id, label in profiles:
        for process_count in process_counts:
            points = []
            for (row_profile, _row_label, count, x_value, columns_per_process), values in grouped.items():
                if row_profile != profile_id or count != process_count:
                    continue
                points.append(
                    {
                        "x": x_value,
                        "y": statistics.fmean(values),
                        "error": statistics.pstdev(values),
                        "columns_per_process": columns_per_process,
                        "repetitions": len(values),
                    }
                )
            if not points:
                continue
            points.sort(key=lambda point: point["x"])
            trace_name = f"{label} · {process_count} proc" if multiple_profiles else f"{process_count} process{'es' if process_count != 1 else ''}"
            figure.add_trace(
                go.Scatter(
                    x=[point["x"] for point in points],
                    y=[point["y"] for point in points],
                    error_y={
                        "type": "data",
                        "array": [point["error"] for point in points],
                        "visible": any(point["error"] > 0 for point in points),
                    },
                    customdata=[[point["columns_per_process"], point["repetitions"]] for point in points],
                    mode="lines+markers",
                    name=trace_name,
                    line={"color": colors[profile_id], "dash": dashes[process_count]},
                    marker={"color": colors[profile_id]},
                    hovertemplate=(
                        "%{fullData.name}<br>"
                        f"{X_AXES[x_axis].lower()}=%{{x}}"
                        "<br>batch size/process=%{customdata[0]}"
                        f"<br>{METRICS[metric]}=%{{y:.6g}}"
                        "<br>repetitions=%{customdata[1]}<extra></extra>"
                    ),
                )
            )

    figure = _apply_layout(
        figure,
        theme_name=theme_name,
        title=f"{'Process-group wall time' if selected_timer == GROUP_WALL_TIMER else selected_timer} · {METRICS[metric]}",
        x_title=X_AXES[x_axis],
        y_title=METRICS[metric],
        x_scale=x_scale,
    )
    if x_scale != "log":
        figure.update_xaxes(rangemode="tozero")
    figure.update_yaxes(rangemode="tozero")
    return figure


def relative_figure(
    rows: Iterable[dict[str, Any]],
    timer_name: str,
    metric: str,
    baseline_id: str,
    theme_name: str = "dark",
    *,
    x_axis: str = "total_columns",
    x_scale: str = "linear",
) -> go.Figure:
    """Plot percent changes at exactly matched workloads against a baseline."""
    if not baseline_id:
        return empty_profile_figure("Choose a baseline profile for relative comparisons.", theme_name)
    metric = metric if metric in METRICS else "throughput_columns_per_second"
    x_axis = x_axis if x_axis in X_AXES else "total_columns"
    grouped: dict[tuple[str, str, int, int, int], list[float]] = defaultdict(list)
    for row in rows:
        if str(row.get("timer_name") or "") != str(timer_name or "") or str(row.get("status") or "") != "success":
            continue
        process_count = _integer(row.get("process_count"))
        columns = _integer(row.get("columns_per_process"))
        total = _integer(row.get("total_columns"))
        value = _number(row.get(metric))
        if None in (process_count, columns, total, value):
            continue
        grouped[(_profile_id(row), _profile_label(row), int(process_count), int(columns), int(total))].append(float(value))
    means = {key: statistics.fmean(values) for key, values in grouped.items()}
    baseline = {
        (process_count, columns, total): value
        for (profile_id, _label, process_count, columns, total), value in means.items()
        if profile_id == baseline_id
    }
    candidates = sorted({(key[0], key[1]) for key in means if key[0] != baseline_id})
    if not baseline or not candidates:
        return empty_profile_figure("Select a baseline and at least one comparison profile.", theme_name)

    figure = go.Figure()
    for profile_id, label in candidates:
        by_process: dict[int, list[tuple[int, float, int]]] = defaultdict(list)
        for (row_profile, _row_label, process_count, columns, total), value in means.items():
            if row_profile != profile_id:
                continue
            base = baseline.get((process_count, columns, total))
            if base is None or base == 0:
                continue
            x_value = total if x_axis == "total_columns" else columns
            by_process[process_count].append((x_value, 100.0 * (value / base - 1.0), columns))
        for process_count, points in sorted(by_process.items()):
            points.sort()
            figure.add_trace(
                go.Scatter(
                    x=[point[0] for point in points],
                    y=[point[1] for point in points],
                    customdata=[[point[2]] for point in points],
                    mode="lines+markers",
                    name=f"{label} · {process_count} proc",
                    line={"color": _series_color(profile_id)},
                    hovertemplate=(
                        "%{fullData.name}<br>"
                        f"{X_AXES[x_axis].lower()}=%{{x}}<br>batch size/process=%{{customdata[0]}}"
                        "<br>change=%{y:+.2f}%<extra></extra>"
                    ),
                )
            )
    if not figure.data:
        return empty_profile_figure("Selected profiles have no exactly matched workloads.", theme_name)
    figure.add_hline(y=0, line_width=1, line_dash="dot")
    figure = _apply_layout(
        figure,
        theme_name=theme_name,
        title=f"Change from baseline · {METRICS[metric]}",
        x_title=X_AXES[x_axis],
        y_title="Change from baseline (%)",
        x_scale=x_scale,
    )
    return figure


def variability_figure(
    process_rows: Iterable[dict[str, Any]],
    profile_id: str,
    timer_name: str,
    process_count: int | None,
    theme_name: str = "dark",
) -> go.Figure:
    """Show the raw cross-process and repetition timing distribution."""
    critical_threads: dict[tuple[int, int, int, int], float] = {}
    label = profile_id
    for row in process_rows:
        if _profile_id(row) != profile_id or str(row.get("timer_name") or "") != str(timer_name or ""):
            continue
        if str(row.get("status") or "") != "success":
            continue
        count = _integer(row.get("process_count"))
        total = _integer(row.get("total_columns"))
        repetition = _integer(row.get("repetition"))
        process_index = _integer(row.get("process_index"))
        value = _number(row.get("inclusive_seconds"))
        if None in (count, total, repetition, process_index, value) or (process_count is not None and count != process_count):
            continue
        label = _profile_label(row)
        key = (int(count), int(total), int(repetition), int(process_index))
        critical_threads[key] = max(float(value), critical_threads.get(key, 0.0))
    grouped: dict[tuple[int, int], list[float]] = defaultdict(list)
    for (count, total, _repetition, _process_index), value in critical_threads.items():
        grouped[(count, total)].append(value)
    if not grouped:
        return empty_profile_figure("No process-level values are available for this selection.", theme_name)
    figure = go.Figure()
    for index, count in enumerate(sorted({key[0] for key in grouped})):
        x_values = []
        y_values = []
        for (row_count, total), values in sorted(grouped.items()):
            if row_count != count:
                continue
            x_values.extend([str(total)] * len(values))
            y_values.extend(values)
        figure.add_trace(
            go.Box(
                x=x_values,
                y=y_values,
                name=f"{count} proc",
                boxpoints="outliers",
                marker={"color": qualitative.Plotly[index % len(qualitative.Plotly)]},
                legendgroup=str(count),
            )
        )
    return _apply_layout(
        figure,
        theme_name=theme_name,
        title=f"Process spread · {label}",
        x_title="Overall batch size",
        y_title="Inclusive timer time (s)",
    )


def decomposition_figure(
    process_rows: Iterable[dict[str, Any]],
    profile_id: str,
    selected_timer: str,
    process_count: int | None,
    top_n: int = 8,
    percentage: bool = False,
    theme_name: str = "dark",
) -> go.Figure:
    """Stack exclusive costs from the critical process of each repetition."""
    if process_count is None:
        process_counts = [
            count
            for row in process_rows
            if _profile_id(row) == profile_id
            and (count := _integer(row.get("process_count"))) is not None
        ]
        process_count = min(process_counts) if process_counts else None
    if process_count is None:
        return empty_profile_figure("No process-level data is available for cost decomposition.", theme_name)
    relevant = [
        row
        for row in process_rows
        if _profile_id(row) == profile_id
        and _integer(row.get("process_count")) == process_count
        and str(row.get("status") or "") == "success"
    ]
    if not relevant:
        return empty_profile_figure("No process-level data is available for cost decomposition.", theme_name)

    # Keep all cost components from the same process/thread that supplied the
    # slowest primary timer.  Independent maxima would create a run that never
    # actually occurred.
    entities: dict[tuple[int, int, int, int], dict[str, float]] = defaultdict(dict)
    primary_values: dict[tuple[int, int, int, int], float] = {}
    labels = []
    for row in relevant:
        repetition = _integer(row.get("repetition"))
        columns = _integer(row.get("columns_per_process"))
        total = _integer(row.get("total_columns"))
        process_index = _integer(row.get("process_index"))
        thread = _integer(row.get("thread")) or 0
        inclusive = _number(row.get("inclusive_seconds"))
        exclusive = _number(row.get("exclusive_seconds"))
        timer = str(row.get("timer_name") or "")
        if None in (repetition, columns, total, process_index) or not timer:
            continue
        key = (int(repetition), int(columns), int(process_index), int(thread))
        if exclusive is not None:
            entities[key][timer] = float(exclusive)
        if timer == selected_timer and inclusive is not None:
            primary_values[key] = float(inclusive)
        labels.append(_profile_label(row))

    critical: dict[tuple[int, int], tuple[int, int, int, int]] = {}
    for key, value in primary_values.items():
        repetition, columns, _process_index, _thread = key
        config = (repetition, columns)
        if config not in critical or value > primary_values[critical[config]]:
            critical[config] = key
    if not critical or not any(entities.get(key) for key in critical.values()):
        return empty_profile_figure(
            "Exclusive costs are unavailable for this timer output.",
            theme_name,
        )

    by_point: dict[tuple[int, str], list[float]] = defaultdict(list)
    all_timers: set[str] = set()
    for (_repetition, columns), entity_key in critical.items():
        total = process_count * columns
        for timer, value in entities.get(entity_key, {}).items():
            by_point[(total, timer)].append(value)
            all_timers.add(timer)
    scores = {
        timer: sum(value for (total, name), values in by_point.items() if name == timer for value in values)
        for timer in all_timers
    }
    selected_timers = [timer for timer, _score in sorted(scores.items(), key=lambda item: item[1], reverse=True)[: max(1, int(top_n))]]
    totals = sorted({key[0] for key in by_point})
    averaged: dict[tuple[int, str], float] = {
        key: statistics.fmean(values) for key, values in by_point.items()
    }
    has_other = any(timer not in selected_timers for timer in all_timers)
    figure = go.Figure()
    series = list(selected_timers) + (["Other"] if has_other else [])
    for index, timer in enumerate(series):
        values = []
        for total in totals:
            if timer == "Other":
                value = sum(averaged.get((total, name), 0.0) for name in all_timers if name not in selected_timers)
            else:
                value = averaged.get((total, timer), 0.0)
            values.append(value)
        figure.add_trace(
            go.Scatter(
                x=totals,
                y=values,
                name=timer,
                mode="lines",
                stackgroup="cost",
                groupnorm="percent" if percentage else None,
                line={"width": 1.5, "color": qualitative.Plotly[index % len(qualitative.Plotly)]},
                hovertemplate=(
                    "%{fullData.name}<br>overall batch size=%{x}<br>cost=%{y:.6g}"
                    + ("%" if percentage else " s")
                    + "<extra></extra>"
                ),
            )
        )
    label = labels[0] if labels else profile_id
    return _apply_layout(
        figure,
        theme_name=theme_name,
        title=f"Exclusive cost decomposition · {label} · {process_count} proc",
        x_title="Overall batch size",
        y_title="Share of exclusive cost (%)" if percentage else "Exclusive timer cost (s)",
    )

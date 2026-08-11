"""Plotly figures for direct CLUBB hole-filling snapshots."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .analysis import MOMENTS, below_threshold_counts, load_moment


def _finish(figure, *, height: int, uirevision: str):
    figure.update_layout(
        template="plotly_dark",
        height=height,
        margin={"l": 72, "r": 34, "t": 62, "b": 92},
        paper_bgcolor="#0f172a",
        plot_bgcolor="#111827",
        font={"color": "#e5e7eb", "size": 12},
        hovermode="closest",
        uirevision=uirevision,
        legend={
            "orientation": "h",
            "y": -0.2,
            "yanchor": "top",
            "x": 0,
            "xanchor": "left",
        },
    )
    figure.update_xaxes(gridcolor="#334155", zerolinecolor="#94a3b8")
    figure.update_yaxes(gridcolor="#334155", zerolinecolor="#94a3b8")
    return figure


def empty_figure(message: str):
    figure = go.Figure()
    figure.add_annotation(text=message, showarrow=False, font={"color": "#cbd5e1"})
    return _finish(figure, height=360, uirevision="hf-empty")


def profile_figure(
    data: dict,
    moment: str,
    record: int,
    column: int,
    field_scale: str = "log",
):
    before = data["before"][record, :, column]
    after = data["after"][record, :, column]
    height = data["height"]
    threshold = float(data["threshold"][record, column])
    changed = np.isfinite(before) & np.isfinite(after) & ~np.isclose(
        before, after, rtol=1.0e-12, atol=0.0
    )

    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=before,
            y=height,
            mode="lines+markers",
            name="before",
            line={"color": "#f59e0b", "width": 2},
            marker={"size": 4},
            hovertemplate="before=%{x:.7g}<br>z=%{y:.2f} m<extra></extra>",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=after,
            y=height,
            mode="lines",
            name="after",
            line={"color": "#34d399", "width": 2},
            hovertemplate="after=%{x:.7g}<br>z=%{y:.2f} m<extra></extra>",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=after[changed],
            y=height[changed],
            mode="markers",
            name="changed levels",
            marker={"color": "#f43f5e", "size": 8, "symbol": "diamond"},
            customdata=before[changed],
            hovertemplate=(
                "before=%{customdata:.7g}<br>after=%{x:.7g}"
                "<br>z=%{y:.2f} m<extra></extra>"
            ),
        )
    )
    if np.isfinite(threshold):
        figure.add_vline(
            x=threshold,
            line={"color": "#fb7185", "dash": "dash", "width": 1.5},
            annotation_text="threshold",
            annotation_position="top",
        )
    figure.update_layout(
        title=f"{MOMENTS[moment]['label']} at one output record",
        xaxis_title=f"{MOMENTS[moment]['label']} ({data['units'] or 'file units'})",
        xaxis_type="linear" if field_scale == "linear" else "log",
        yaxis_title="Altitude (m)",
    )
    return _finish(
        figure,
        height=610,
        uirevision=f"hf-profile-{moment}-{field_scale}",
    )


def holes_over_time_figure(data: dict, moment: str, column: int):
    threshold = data["threshold"][:, column, np.newaxis]
    before = data["before"][:, :, column]
    after = data["after"][:, :, column]
    before_holes = np.where(np.isfinite(before), (before < threshold).astype(float), np.nan)
    after_holes = np.where(np.isfinite(after), (after < threshold).astype(float), np.nan)
    time_minutes = data["times"] / 60.0

    figure = make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True,
        subplot_titles=("Before fill", "Still below threshold after fill"),
        horizontal_spacing=0.08,
    )
    colorscale = [[0.0, "#111827"], [0.499, "#111827"], [0.5, "#f43f5e"], [1.0, "#f43f5e"]]
    for column_index, values in enumerate((before_holes, after_holes), start=1):
        figure.add_trace(
            go.Heatmap(
                x=time_minutes,
                y=data["height"],
                z=values.T,
                zmin=0,
                zmax=1,
                colorscale=colorscale,
                showscale=False,
                hovertemplate="t=%{x:.2f} min<br>z=%{y:.2f} m<extra></extra>",
            ),
            row=1,
            col=column_index,
        )
    figure.update_xaxes(title_text="Model time (min)")
    figure.update_yaxes(title_text="Altitude (m)", row=1, col=1)
    figure.update_layout(title=f"Where {MOMENTS[moment]['label']} falls below its fill threshold")
    return _finish(figure, height=470, uirevision=f"hf-map-{moment}-{column}")


def summary_figure(path: str, moments: list[str], column: int):
    figure = go.Figure()
    colors = ["#38bdf8", "#a78bfa", "#f59e0b", "#34d399", "#f472b6", "#fb7185", "#22d3ee"]
    for color, moment in zip(colors, moments):
        data = load_moment(path, moment)
        before_count, after_count = below_threshold_counts(data, column)
        figure.add_trace(
            go.Scatter(
                x=data["times"] / 60.0,
                y=before_count,
                mode="lines",
                name=MOMENTS[moment]["label"],
                line={"color": color, "width": 1.8},
                customdata=after_count,
                hovertemplate=(
                    "t=%{x:.2f} min<br>before=%{y} levels"
                    "<br>after=%{customdata} levels<extra>%{fullData.name}</extra>"
                ),
            )
        )
    figure.update_layout(
        title="Levels below threshold before each fill",
        xaxis_title="Model time (min)",
        yaxis_title="Level count",
    )
    return _finish(figure, height=420, uirevision=f"hf-summary-{column}")

"""Plotly figures for the interactive mixing-length trajectory explorer."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from plotly.subplots import make_subplots

from .analysis import GRAV

# Every stop is deliberately high-contrast against the dark plot background.
# Unlike the low end of Turbo/Viridis, launch heights near the surface remain
# clearly visible without changing trajectory opacity or line width.
BRIGHT_HEIGHT_SCALE = [
    "#60a5fa",
    "#22d3ee",
    "#2dd4bf",
    "#4ade80",
    "#a3e635",
    "#fde047",
    "#fb923c",
    "#f87171",
    "#e879f9",
]


def _finish(figure, *, height=760):
    figure.update_layout(
        template="plotly_dark",
        height=height,
        margin={"l": 72, "r": 35, "t": 70, "b": 62},
        paper_bgcolor="#0f172a",
        plot_bgcolor="#111827",
        font={"color": "#e5e7eb", "size": 13},
        hovermode="closest",
        uirevision="mixing-length-trajectories",
    )
    figure.update_xaxes(gridcolor="#334155", zerolinecolor="#94a3b8")
    figure.update_yaxes(gridcolor="#334155", zerolinecolor="#94a3b8")
    return figure


def trajectory_figure(result, direction):
    """Draw every launch's remaining-energy path for one direction."""
    if direction == "up":
        paths = result.upward_paths
        title = "Upward parcels"
    else:
        paths = result.downward_paths
        title = "Downward parcels"
    colorscale = BRIGHT_HEIGHT_SCALE

    figure = go.Figure()
    z_min = float(result.z[0])
    z_span = max(float(result.z[-1] - result.z[0]), 1.0)
    for path in paths:
        fraction = (float(result.z[path.launch_index]) - z_min) / z_span
        color = sample_colorscale(colorscale, [fraction])[0]
        launch = float(result.z[path.launch_index])
        figure.add_trace(
            go.Scattergl(
                x=path.energy,
                y=path.altitude,
                mode="lines",
                line={"color": color, "width": 1.25},
                opacity=0.78,
                showlegend=False,
                customdata=np.full(path.energy.shape, launch),
                hovertemplate=(
                    "launch z=%{customdata:.0f} m"
                    "<br>parcel z=%{y:.1f} m"
                    "<br>K=%{x:.6g} m² s⁻²<extra></extra>"
                ),
            )
        )

    # The launch markers provide a continuous color key without adding one
    # legend item for every spaghetti line.
    launch_indices = np.asarray([path.launch_index for path in paths])
    figure.add_trace(
        go.Scatter(
            x=result.tke[launch_indices],
            y=result.z[launch_indices],
            mode="markers",
            name="launch",
            marker={
                "size": 5,
                "color": result.z[launch_indices],
                "colorscale": colorscale,
                "showscale": True,
                "colorbar": {"title": "launch z (m)", "thickness": 13},
                "line": {"width": 0},
            },
            hovertemplate=(
                "launch z=%{y:.0f} m"
                "<br>initial K=%{x:.6g} m² s⁻²<extra></extra>"
            ),
        )
    )
    figure.add_vline(
        x=0.0, line={"color": "#f8fafc", "width": 1.2, "dash": "dot"}
    )
    figure.update_layout(
        title=f"{title}: remaining parcel energy K",
        showlegend=False,
    )
    figure.update_xaxes(title_text="Remaining K (m² s⁻²)")
    figure.update_yaxes(title_text="Altitude (m)")
    return _finish(figure)


def buoyancy_trajectory_figure(result, direction):
    """Draw the entraining parcels' diagnosed buoyancy along every path."""
    if direction == "up":
        paths = result.upward_paths
        title = "Upward parcels"
    else:
        paths = result.downward_paths
        title = "Downward parcels"

    figure = go.Figure()
    z_min = float(result.z[0])
    z_span = max(float(result.z[-1] - result.z[0]), 1.0)
    for path in paths:
        fraction = (float(result.z[path.launch_index]) - z_min) / z_span
        color = sample_colorscale(BRIGHT_HEIGHT_SCALE, [fraction])[0]
        launch = float(result.z[path.launch_index])
        figure.add_trace(
            go.Scattergl(
                x=path.buoyancy,
                y=path.altitude,
                mode="lines",
                line={"color": color, "width": 1.25},
                opacity=0.78,
                showlegend=False,
                customdata=np.full(path.buoyancy.shape, launch),
                hovertemplate=(
                    "launch z=%{customdata:.0f} m"
                    "<br>parcel z=%{y:.1f} m"
                    "<br>B=%{x:.6g} m s⁻²<extra></extra>"
                ),
            )
        )

    launch_indices = np.asarray([path.launch_index for path in paths])
    figure.add_trace(
        go.Scatter(
            x=np.zeros(launch_indices.size),
            y=result.z[launch_indices],
            mode="markers",
            name="launch",
            marker={
                "size": 5,
                "color": result.z[launch_indices],
                "colorscale": BRIGHT_HEIGHT_SCALE,
                "showscale": True,
                "colorbar": {"title": "launch z (m)", "thickness": 13},
                "line": {"width": 0},
            },
            hovertemplate="launch z=%{y:.0f} m<br>B=0 m s⁻²<extra></extra>",
        )
    )
    figure.add_vline(
        x=0.0, line={"color": "#f8fafc", "width": 1.2, "dash": "dot"}
    )
    figure.update_layout(
        title=f"{title}: parcel buoyancy along the diagnosed path",
        showlegend=False,
    )
    figure.update_xaxes(title_text="Parcel buoyancy B (m s⁻²)")
    figure.update_yaxes(title_text="Altitude (m)")
    return _finish(figure)


def parcel_state_figure(record, result, field):
    """Plot one upward-parcel state family alongside the grid mean."""
    field_spec = {
        "thv": {
            "title": "Virtual potential temperature",
            "axis": "θᵥ (K)",
            "mean_label": "Grid mean thvm",
            "mean": record.thvm,
            # This is algebraically identical to the buoyancy used by the
            # parcel-energy integral.  Interpolating thvm also supports the
            # diagnosed stopping point inside a model layer.
            "path": lambda path: np.interp(
                path.altitude, result.z, record.thvm
            )
            * (1.0 + path.buoyancy / GRAV),
            "hover": "θᵥ=%{x:.5g} K",
        },
        "rt": {
            "title": "Total water",
            "axis": "rₜ (g kg⁻¹)",
            "mean_label": "Grid mean rtm",
            "mean": 1000.0 * record.mean_rtm,
            "path": lambda path: 1000.0 * path.parcel_rt,
            "hover": "rₜ=%{x:.5g} g kg⁻¹",
        },
        "thl": {
            "title": "Liquid-water potential temperature",
            "axis": "θₗ (K)",
            "mean_label": "Grid mean thlm",
            "mean": record.mean_thlm,
            "path": lambda path: path.parcel_thl,
            "hover": "θₗ=%{x:.5g} K",
        },
    }
    if field not in field_spec:
        raise ValueError(f"Unknown parcel-state field: {field}")
    spec = field_spec[field]
    figure = go.Figure()
    z_min = float(result.z[0])
    z_span = max(float(result.z[-1] - result.z[0]), 1.0)

    for path in result.upward_paths:
        fraction = (float(result.z[path.launch_index]) - z_min) / z_span
        color = sample_colorscale(BRIGHT_HEIGHT_SCALE, [fraction])[0]
        launch = float(result.z[path.launch_index])
        values = np.asarray(spec["path"](path))
        figure.add_trace(
            go.Scattergl(
                x=values,
                y=path.altitude,
                mode="lines",
                line={"color": color, "width": 1.25},
                opacity=0.78,
                showlegend=False,
                customdata=np.full(values.shape, launch),
                hovertemplate=(
                    "launch z=%{customdata:.0f} m"
                    "<br>parcel z=%{y:.1f} m"
                    f"<br>{spec['hover']}<extra></extra>"
                ),
            )
        )

    figure.add_trace(
        go.Scatter(
            x=spec["mean"],
            y=result.z,
            mode="lines",
            name=spec["mean_label"],
            line={"color": "#f8fafc", "width": 3.0},
            hovertemplate=(
                "z=%{y:.1f} m"
                "<br>mean=%{x:.5g}<extra>"
                + spec["mean_label"]
                + "</extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Upward parcels: {spec['title']}",
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "left",
            "x": 0,
        },
    )
    figure.update_xaxes(title_text=spec["axis"])
    figure.update_yaxes(title_text="Altitude (m)")
    return _finish(figure)


def length_profile_figure(record, result):
    """Compare raw, finalized, and stored directional/final profiles."""
    figure = make_subplots(
        rows=1,
        cols=3,
        shared_yaxes=True,
        horizontal_spacing=0.065,
        subplot_titles=("Upward length", "Downward length", "Final Lscale"),
    )
    panels = (
        (
            result.raw_up,
            result.enveloped_up,
            result.lscale_up,
            record.fortran_up,
            "L↑ (m)",
        ),
        (
            result.raw_down,
            result.enveloped_down,
            result.lscale_down,
            record.fortran_down,
            "L↓ (m)",
        ),
        (
            None,
            None,
            result.lscale,
            record.fortran_lscale,
            "L (m)",
        ),
    )
    for column, (
        raw,
        enveloped,
        calculated,
        reference,
        axis_title,
    ) in enumerate(
        panels, start=1
    ):
        if raw is not None:
            figure.add_trace(
                go.Scatter(
                    x=raw,
                    y=result.z,
                    mode="lines",
                    name="Own parcel endpoint + 0.1 m",
                    legendgroup="own-parcel",
                    showlegend=column == 1,
                    line={"color": "#64748b", "width": 1.2, "dash": "dot"},
                    hovertemplate=(
                        "z=%{y:.0f} m<br>own reach=%{x:.4g} m"
                        "<extra>own parcel endpoint</extra>"
                    ),
                ),
                row=1,
                col=column,
            )
        if enveloped is not None:
            figure.add_trace(
                go.Scatter(
                    x=enveloped,
                    y=result.z,
                    mode="lines",
                    name="After nonlocal envelope",
                    legendgroup="envelope",
                    showlegend=column == 1,
                    line={"color": "#a78bfa", "width": 1.5, "dash": "dot"},
                    hovertemplate=(
                        "z=%{y:.0f} m<br>enveloped=%{x:.4g} m"
                        "<extra>nonlocal envelope</extra>"
                    ),
                ),
                row=1,
                col=column,
            )
        figure.add_trace(
            go.Scatter(
                x=calculated,
                y=result.z,
                mode="lines",
                name="Python replica",
                legendgroup="python",
                showlegend=column == 1,
                line={"color": "#38bdf8", "width": 2.5},
                hovertemplate=(
                    "z=%{y:.0f} m<br>length=%{x:.5g} m"
                    "<extra>Python replica</extra>"
                ),
            ),
            row=1,
            col=column,
        )
        figure.add_trace(
            go.Scatter(
                x=reference,
                y=result.z,
                mode="lines",
                name="Fortran NetCDF",
                legendgroup="fortran",
                showlegend=column == 1,
                line={"color": "#f97316", "width": 1.7, "dash": "dash"},
                hovertemplate=(
                    "z=%{y:.0f} m<br>length=%{x:.5g} m"
                    "<extra>Fortran NetCDF</extra>"
                ),
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text=axis_title, row=1, col=column)
    figure.update_yaxes(title_text="Altitude (m)", row=1, col=1)
    figure.update_layout(
        title=(
            "Directional parcel reach and final geometric-mean mixing length"
        ),
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.04,
            "xanchor": "left",
            "x": 0,
        },
    )
    return _finish(figure, height=720)

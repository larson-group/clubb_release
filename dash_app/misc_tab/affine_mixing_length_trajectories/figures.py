"""Plot parcel paths and coherent mixing lengths from saved-column diagnostics.

Usage:
  - ``parcel_matrix_figure(fields)`` plots parcel energy and buoyancy.
  - ``coherence_matrix_figure(result)`` plots energy barriers and survival weights.
  - ``coherence_sum_figure(result)`` sums arriving weights by target level.
  - ``coherence_length_profile_figure(record, result)`` compares calculated and
    saved upward, downward, and combined lengths.
  - ``coherence_temporal_figure(series)`` compares nearby records and step changes.
  Each function returns a Plotly figure; none writes files.

Understanding:
  Analysis matrices use (launch,target) axes and are transposed for display:
  launch altitude is horizontal and target altitude vertical, both in m.
  Energy/barrier units are m^2/s^2, buoyancy m/s^2, coherence dimensionless,
  length m, and time min. Log color transforms retain physical hover values
  and labeled ticks. Saved Fortran curves remain as recorded when controls
  change; the figures are diagnostics rather than numerical parity tests.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dash_app.misc_tab.mixing_length_trajectories.figures import _finish


def _balanced_signed_log(values):
    """Return signed-log colors, tick positions, and physical-value tick labels.

    Scale negative and positive values independently onto [-1,0] and [0,1],
    using finite-value percentiles and clipping extremes. Keep original values
    outside this transform for physical hover data.
    """

    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    nonzero = np.abs(finite[np.abs(finite) > np.finfo(float).eps])
    if nonzero.size == 0:
        return np.zeros_like(values), [0.0], ["0"]

    threshold = max(float(np.percentile(nonzero, 5.0)), np.finfo(float).eps)
    negative = -finite[finite < 0.0]
    positive = finite[finite > 0.0]
    negative_limit = max(
        float(np.percentile(negative, 99.0)) if negative.size else 0.0,
        threshold,
    )
    positive_limit = max(
        float(np.percentile(positive, 99.0)) if positive.size else 0.0,
        threshold,
    )
    negative_denominator = np.log1p(negative_limit / threshold)
    positive_denominator = np.log1p(positive_limit / threshold)

    transformed = np.zeros_like(values)
    negative_mask = values < 0.0
    positive_mask = values > 0.0
    transformed[negative_mask] = -np.log1p(
        -values[negative_mask] / threshold
    ) / negative_denominator
    transformed[positive_mask] = np.log1p(
        values[positive_mask] / threshold
    ) / positive_denominator
    transformed = np.clip(transformed, -1.0, 1.0)

    tick_values = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    physical_ticks = np.concatenate(
        (
            -threshold
            * np.expm1(-tick_values[:2] * negative_denominator),
            [0.0],
            threshold * np.expm1(tick_values[3:] * positive_denominator),
        )
    )
    return (
        transformed,
        tick_values.tolist(),
        [f"{value:.3g}" for value in physical_ticks],
    )


def _matrix_customdata(raw_values, direction):
    """Return (target,launch,2) hover data containing direction labels and original
    physical values.
    """
    custom = np.empty((*raw_values.shape, 2), dtype=object)
    custom[..., 0] = direction
    custom[..., 1] = raw_values
    return np.transpose(custom, (1, 0, 2))


def _positive_log(values):
    """Return nonnegative log colors, tick positions, and physical labels using
    finite-value percentile scaling.
    """

    values = np.asarray(values, dtype=float)
    positive = values[np.isfinite(values) & (values > np.finfo(float).eps)]
    if positive.size == 0:
        return np.zeros_like(values), [0.0], ["0"]
    scale = max(float(np.percentile(positive, 5.0)), np.finfo(float).eps)
    limit = max(float(np.percentile(positive, 99.0)), scale)
    denominator = np.log1p(limit / scale)
    transformed = np.clip(
        np.log1p(np.maximum(values, 0.0) / scale) / denominator, 0.0, 1.0
    )
    tick_values = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    physical_ticks = scale * np.expm1(tick_values * denominator)
    return (
        transformed,
        tick_values.tolist(),
        [f"{value:.3g}" for value in physical_ticks],
    )


def parcel_matrix_figure(fields):
    """Return energy/buoyancy heatmaps and deterministic stopping-reference markers.

    fields provides (launch,target) energy [m^2/s^2], upward-positive buoyancy
    [m/s^2], and stop heights [m]. Colors use separate signed logarithmic
    scales; hover data retains the untransformed physical quantities.
    """

    energy_color, energy_ticks, energy_tick_labels = _balanced_signed_log(
        fields.energy_matrix
    )
    buoyancy_color, buoyancy_ticks, buoyancy_tick_labels = (
        _balanced_signed_log(fields.buoyancy_matrix)
    )
    launch, target = np.indices(fields.energy_matrix.shape)
    direction = np.where(
        target > launch,
        "up",
        np.where(target < launch, "down", "launch"),
    )
    energy_custom = _matrix_customdata(fields.energy_matrix, direction)
    buoyancy_custom = _matrix_customdata(fields.buoyancy_matrix, direction)
    figure = make_subplots(
        rows=1,
        cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=0.11,
        subplot_titles=("Parcel energy", "Parcel buoyancy"),
    )
    figure.add_trace(
        go.Heatmap(
            x=fields.z,
            y=fields.z,
            z=energy_color.T,
            zmin=-1.0,
            zmax=1.0,
            zmid=0.0,
            colorscale="RdBu_r",
            colorbar={
                "title": "K<br>m² s⁻²<br>signed log",
                "x": 0.43,
                "len": 0.82,
                "tickvals": energy_ticks,
                "ticktext": energy_tick_labels,
            },
            customdata=energy_custom,
            hovertemplate=(
                "launch=%{x:.0f} m<br>target=%{y:.0f} m"
                "<br>K=%{customdata[1]:.6g} m² s⁻²"
                "<br>%{customdata[0]}<extra></extra>"
            ),
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Heatmap(
            x=fields.z,
            y=fields.z,
            z=buoyancy_color.T,
            zmin=-1.0,
            zmax=1.0,
            zmid=0.0,
            colorscale="BrBG",
            colorbar={
                "title": "B<br>m s⁻²<br>signed log",
                "x": 1.0,
                "len": 0.82,
                "tickvals": buoyancy_ticks,
                "ticktext": buoyancy_tick_labels,
            },
            customdata=buoyancy_custom,
            hovertemplate=(
                "launch=%{x:.0f} m<br>target=%{y:.0f} m"
                "<br>B=%{customdata[1]:.6g} m s⁻²"
                "<br>%{customdata[0]}<extra></extra>"
            ),
        ),
        row=1,
        col=2,
    )
    for endpoints, name, color in (
        (
            fields.legacy_upward_stops,
            "Legacy first upward stop",
            "#000000",
        ),
        (
            fields.legacy_downward_stops,
            "Legacy first downward stop",
            "#064e3b",
        ),
    ):
        active = np.isfinite(endpoints)
        for column in (1, 2):
            figure.add_trace(
                go.Scattergl(
                    x=fields.z[active],
                    y=endpoints[active],
                    mode="markers",
                    name=name,
                    legendgroup=name,
                    showlegend=column == 1,
                    marker={"size": 4, "color": color},
                    hovertemplate=(
                        "launch=%{x:.0f} m<br>stop=%{y:.1f} m"
                        f"<extra>{name}</extra>"
                    ),
                ),
                row=1,
                col=column,
            )
    for column in (1, 2):
        figure.update_xaxes(title_text="Launch altitude (m)", row=1, col=column)
    figure.update_yaxes(title_text="Target altitude (m)", row=1, col=1)
    figure.update_layout(
        title=(
            "Dense affine parcel fields"
            "<br><sup>Ascent lies above the diagonal and descent below it; "
            "each sign receives equal color range.</sup>"
        ),
        legend={
            "orientation": "h",
            "yanchor": "top",
            "y": -0.13,
            "xanchor": "center",
            "x": 0.5,
        },
    )
    figure = _finish(figure, height=700)
    figure.update_layout(margin={"l": 72, "r": 35, "t": 88, "b": 112})
    return figure


def coherence_matrix_figure(result):
    """Return barrier and retained-coherence heatmaps with a q=0.5 contour.

    result contains (launch,target) barriers [m^2/s^2] and endpoint survival
    weights [-]. The barrier color scale is logarithmic; coherence is linear.
    """

    barrier_color, barrier_ticks, barrier_tick_labels = _positive_log(
        result.barrier_matrix
    )
    launch, target = np.indices(result.barrier_matrix.shape)
    direction = np.where(
        target > launch,
        "up",
        np.where(target < launch, "down", "launch"),
    )
    barrier_custom = _matrix_customdata(result.barrier_matrix, direction)
    coherence_custom = _matrix_customdata(result.coherence_matrix, direction)
    figure = make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True,
        horizontal_spacing=0.08,
        subplot_titles=(
            "Required adverse-work barrier",
            "Inverse-square coherent fraction",
        ),
    )
    figure.add_trace(
        go.Heatmap(
            x=result.z,
            y=result.z,
            z=barrier_color.T,
            zmin=0.0,
            zmax=1.0,
            colorscale="Magma",
            colorbar={
                "title": "H<br>m² s⁻²<br>log",
                "x": 0.46,
                "len": 0.82,
                "tickvals": barrier_ticks,
                "ticktext": barrier_tick_labels,
            },
            customdata=barrier_custom,
            hovertemplate=(
                "launch=%{x:.0f} m<br>target=%{y:.0f} m"
                "<br>H=%{customdata[1]:.6g} m² s⁻²"
                "<br>%{customdata[0]}<extra></extra>"
            ),
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Heatmap(
            x=result.z,
            y=result.z,
            z=result.coherence_matrix.T,
            zmin=0.0,
            zmax=1.0,
            colorscale="Viridis",
            colorbar={"title": "coherent<br>fraction", "x": 1.01, "len": 0.82},
            customdata=coherence_custom,
            hovertemplate=(
                "launch=%{x:.0f} m<br>target=%{y:.0f} m"
                "<br>q=%{customdata[1]:.6f}"
                "<br>%{customdata[0]}<extra></extra>"
            ),
        ),
        row=1,
        col=2,
    )
    figure.add_trace(
        go.Contour(
            x=result.z,
            y=result.z,
            z=result.coherence_matrix.T,
            contours={
                "start": 0.5,
                "end": 0.5,
                "size": 0.1,
                "coloring": "none",
            },
            line={"color": "#000000", "width": 2.0},
            showscale=False,
            hoverinfo="skip",
        ),
        row=1,
        col=2,
    )
    for column in (1, 2):
        figure.update_xaxes(title_text="Launch altitude (m)", row=1, col=column)
    figure.update_yaxes(title_text="Target altitude (m)", row=1, col=1)
    figure.update_layout(
        title=(
            "Affine barrier and inverse-square coherence"
            f"<br><sup>TKE multiplier={result.tke_multiplier:.3g}; "
            "paths include the first layer at "
            "or below 1% coherence, then stop.</sup>"
        )
    )
    return _finish(figure, height=650)


def coherence_sum_figure(result):
    """Return upward, downward, and total survival-weight profiles by target.

    result supplies the retained (launch,target) coherence matrix. Each
    separate directional sum includes the unit local launch; the total counts
    it once. Displayed sums include the lid launch, unlike the active-domain
    pooled-length calculation.
    """

    def directional_sums(coherence):
        """Return upward/downward sums including the diagonal, and the total counting it once."""
        upward = np.array(
            [
                np.sum(coherence[: target + 1, target])
                for target in range(result.z.size)
            ]
        )
        downward = np.array(
            [
                np.sum(coherence[target:, target])
                for target in range(result.z.size)
            ]
        )
        return upward, downward, np.sum(coherence, axis=0)

    inverse_sums = directional_sums(result.coherence_matrix)
    figure = make_subplots(
        rows=1,
        cols=3,
        shared_yaxes=True,
        horizontal_spacing=0.065,
        subplot_titles=(
            "Upward coherence sum",
            "Downward coherence sum",
            "Total coherence sum",
        ),
    )
    for column, values in enumerate(inverse_sums, start=1):
        figure.add_trace(
            go.Scatter(
                x=values,
                y=result.z,
                mode="lines",
                name="Inverse-square coherence",
                legendgroup="inverse-square",
                showlegend=column == 1,
                line={"color": "#f97316", "width": 2.7},
                hovertemplate=(
                    "target=%{y:.0f} m<br>Σ coherence=%{x:.6g}"
                    "<extra>Inverse-square coherence</extra>"
                ),
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text="Σ coherence", row=1, col=column)
    figure.update_yaxes(title_text="Target altitude (m)", row=1, col=1)
    figure.update_layout(
        title="Coherent launch counts reaching each target level",
        legend={
            "orientation": "h",
            "yanchor": "top",
            "y": -0.13,
            "xanchor": "center",
            "x": 0.5,
        },
    )
    figure = _finish(figure, height=650)
    figure.update_layout(margin={"l": 72, "r": 35, "t": 82, "b": 108})
    return figure


def coherence_length_profile_figure(record, result):
    """Return three length panels [m] comparing result with the unchanged saved record profiles."""

    figure = make_subplots(
        rows=1,
        cols=3,
        shared_yaxes=True,
        horizontal_spacing=0.065,
        subplot_titles=("Upward length", "Downward length", "Final Lscale"),
    )
    panels = (
        (
            result.lscale_up,
            record.fortran_up,
        ),
        (
            result.lscale_down,
            record.fortran_down,
        ),
        (
            result.lscale,
            record.fortran_lscale,
        ),
    )
    styles = (
        (
            "Inverse-square coherence",
            "inverse-square",
            "#f97316",
            None,
            2.8,
        ),
        ("Saved Fortran", "fortran", "#22c55e", "longdash", 2.4),
    )
    for column, profiles in enumerate(panels, start=1):
        for profile, (name, group, color, dash, width) in zip(profiles, styles):
            line = {"color": color, "width": width}
            if dash:
                line["dash"] = dash
            figure.add_trace(
                go.Scatter(
                    x=profile,
                    y=result.z,
                    mode="lines",
                    name=name,
                    legendgroup=group,
                    showlegend=column == 1,
                    line=line,
                    hovertemplate=(
                        "z=%{y:.0f} m<br>length=%{x:.6g} m"
                        f"<extra>{name}</extra>"
                    ),
                ),
                row=1,
                col=column,
            )
        figure.update_xaxes(title_text="Length (m)", row=1, col=column)
    figure.update_yaxes(title_text="Altitude (m)", row=1, col=1)
    overwrite_label = {
        "no_overwrite": "no overwrite",
        "legacy_max": "legacy max overwrite",
        "coherence_averaged": "coherence-averaged overwrite",
    }.get(result.overwrite_mode, result.overwrite_mode)
    figure.update_layout(
        title=(
            "Inverse-square mixing lengths"
            f"<br><sup>Calculated profiles use {overwrite_label}; saved Fortran "
            "output is unaffected by this selection.</sup>"
        ),
        legend={
            "orientation": "h",
            "yanchor": "top",
            "y": -0.13,
            "xanchor": "center",
            "x": 0.5,
        },
    )
    figure = _finish(figure, height=760)
    figure.update_layout(margin={"l": 72, "r": 35, "t": 88, "b": 145})
    return figure


def coherence_temporal_figure(series):
    """Return selected-level lengths and profile RMS step changes [m] over time [min].

    series is the mapping from callbacks._temporal_series; mark its selected
    time and compare recomputed values with the saved Fortran series.
    """

    figure = make_subplots(
        rows=1,
        cols=2,
        horizontal_spacing=0.1,
        subplot_titles=(
            f"Lscale at z={series['sensitive_z']:.0f} m",
            "Whole-profile timestep roughness",
        ),
    )
    for values, name, color, dash in (
        (
            series["inverse_square_at_level"],
            "Inverse-square coherence",
            "#f97316",
            None,
        ),
        (series["fortran_at_level"], "Saved Fortran", "#38bdf8", "dash"),
    ):
        line = {"color": color, "width": 2.7}
        if dash:
            line["dash"] = dash
        figure.add_trace(
            go.Scatter(
                x=series["time_minutes"],
                y=values,
                mode="lines+markers",
                name=name,
                legendgroup=name,
                line=line,
            ),
            row=1,
            col=1,
        )
    step_times = 0.5 * (
        series["time_minutes"][:-1] + series["time_minutes"][1:]
    )
    for values, name, color, dash in (
        (
            series["inverse_square_roughness"],
            "Inverse-square coherence",
            "#f97316",
            None,
        ),
        (series["fortran_roughness"], "Saved Fortran", "#38bdf8", "dash"),
    ):
        line = {"color": color, "width": 2.7}
        if dash:
            line["dash"] = dash
        figure.add_trace(
            go.Scatter(
                x=step_times,
                y=values,
                mode="lines+markers",
                name=name,
                legendgroup=name,
                showlegend=False,
                line=line,
            ),
            row=1,
            col=2,
        )
    figure.add_vline(
        x=series["selected_time_minutes"],
        line={"color": "#64748b", "width": 1.2, "dash": "dot"},
        row=1,
        col=1,
    )
    figure.update_xaxes(title_text="Model time (min)", row=1, col=1)
    figure.update_xaxes(title_text="Model time (min)", row=1, col=2)
    figure.update_yaxes(title_text="Lscale (m)", row=1, col=1)
    figure.update_yaxes(title_text="RMS ΔLscale (m)", row=1, col=2)
    figure.update_layout(
        title=(
            f"Temporal window, records {series['record_indices'][0]}–"
            f"{series['record_indices'][-1]}"
        ),
        legend={
            "orientation": "h",
            "yanchor": "top",
            "y": -0.14,
            "xanchor": "center",
            "x": 0.5,
        },
    )
    figure = _finish(figure, height=620)
    figure.update_layout(margin={"l": 72, "r": 35, "t": 82, "b": 108})
    return figure


__all__ = [
    "coherence_length_profile_figure",
    "coherence_matrix_figure",
    "coherence_sum_figure",
    "coherence_temporal_figure",
    "parcel_matrix_figure",
]

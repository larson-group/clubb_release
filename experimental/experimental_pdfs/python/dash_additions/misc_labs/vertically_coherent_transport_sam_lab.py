"""Raw-SAM laboratory for iterative, Lscale-coupled PDF coherence."""

from __future__ import annotations

from dash import Input, Output, dcc, html
import numpy as np
import plotly.graph_objects as go

if __package__ and __package__.startswith("dash_app."):
    from dash_app.shared.components import styled_dropdown
    from dash_app.shared.reporting import empty_state, fact_grid, report_header
else:  # Script-style imports used by dash_app/app.py.
    from shared.components import styled_dropdown
    from shared.reporting import empty_state, fact_grid, report_header

from ..registry import ReportSpec, register_report
from ..transport_2g_prototype import (
    diagnose_transport_2g_from_moments,
    diagnose_transport_2g_with_mixed_center_alignment,
)
from ..vertical_coherence import (
    VerticalCoherenceSettings,
    apply_iterative_vertical_coherence_column,
    standardized_displacement,
)
from ..sam_lscale import (
    SamLscaleError,
    compute_sam_lscale_profile,
    compute_sam_plume_lscale_profile,
)
from . import trivariate_transport_sam_lab as base
from utilities.sam_3d_reference import load_profile_moments


CASE_ID = "notes-vertical-coherence-sam-case"
TIME_ID = "notes-vertical-coherence-sam-time"
HEIGHT_ID = "notes-vertical-coherence-sam-height"
TIME_LABEL_ID = "notes-vertical-coherence-sam-time-label"
HEIGHT_LABEL_ID = "notes-vertical-coherence-sam-height-label"
ENABLED_ID = "notes-vertical-coherence-enabled"
BLEND_ID = "notes-vertical-coherence-blend"
BLEND_VALUE_ID = "notes-vertical-coherence-blend-value"
ITERATIONS_ID = "notes-vertical-coherence-iterations"
ITERATIONS_VALUE_ID = "notes-vertical-coherence-iterations-value"
MIXED_ENABLED_ID = "notes-vertical-coherence-mixed-enabled"
MIXED_STRENGTH_ID = "notes-vertical-coherence-mixed-strength"
MIXED_STRENGTH_VALUE_ID = "notes-vertical-coherence-mixed-strength-value"
FACTS_ID = "notes-vertical-coherence-facts"
STATUS_ID = "notes-vertical-coherence-status"
FIGURE_IDS = {
    projection: f"notes-vertical-coherence-{projection}" for projection in base.PROJECTIONS
}
RT_THL_FIGURE_ID = "notes-vertical-coherence-rt-thl"
LSCALE_FIGURE_ID = "notes-vertical-coherence-lscale"


def _input_id(name: str) -> str:
    return f"notes-vertical-coherence-{name}"


def _value_id(name: str) -> str:
    return f"notes-vertical-coherence-{name}-value"


def _tuning_from_values(values: tuple[float, ...]):
    return base._tuning_from_values(values)


def _coherence_controls():
    return html.Section(
        [
            html.H2("Vertical coherence · Design A"),
            html.P(
                "Iterate between the diagnosed G1 center path, native CLUBB direct-Lscale reach, and a bounded source-to-target center prior. Every pass is anchored to each level's original moment-only PDF.",
                className="transport-sam-controls-intro",
            ),
            dcc.Checklist(
                id=ENABLED_ID,
                options=[
                    {
                        "label": " Enable iterative vertical coherence",
                        "value": "enabled",
                    }
                ],
                value=[],
                className="transport-sam-coherence-toggle",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Maximum local blend", htmlFor=BLEND_ID),
                            html.Output("15%", id=BLEND_VALUE_ID),
                        ],
                        className="transport-sam-control-heading",
                    ),
                    dcc.Slider(
                        id=BLEND_ID,
                        min=0.0,
                        max=0.25,
                        step=0.01,
                        value=0.15,
                        marks=None,
                        updatemode="drag",
                    ),
                    html.Small(
                        "Hard capped at 25%. The local selected-plane diagnosis remains dominant."
                    ),
                ],
                className="transport-sam-control",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("PDF ↔ Lscale passes", htmlFor=ITERATIONS_ID),
                            html.Output("2", id=ITERATIONS_VALUE_ID),
                        ],
                        className="transport-sam-control-heading",
                    ),
                    dcc.Slider(
                        id=ITERATIONS_ID,
                        min=1,
                        max=5,
                        step=1,
                        value=2,
                        marks={value: str(value) for value in range(1, 6)},
                    ),
                    html.Small(
                        "Each pass rebuilds the full G1 thermodynamic path and reruns native calc_Lscale_directly. The solve stops early when the center geometry converges."
                    ),
                ],
                className="transport-sam-control",
            ),
            html.P(
                "All directly reachable lower sources use their Lscale_up; upper sources use their Lscale_down. Distance, activity, and center-direction compatibility determine support.",
                className="transport-sam-complexity-note",
            ),
            html.P(
                "The native calculation follows the current G1 rₜ–θₗ center curve against the raw-SAM environmental θᵥ and TKE proxy. It is a laboratory fixed point, not a new prognostic variable or matrix solve.",
                className="transport-sam-complexity-note",
            ),
        ],
        className="transport-sam-rail-section",
    )


def _mixed_moment_controls():
    """Render the bounded mixed-moment center-direction experiment."""

    return html.Section(
        [
            html.H2("Mixed-moment center direction"),
            html.P(
                "Use raw-SAM conditional-tail information to make one bounded rₜ–θₗ center-direction correction. It preserves the local center-separation budget instead of trying to force more internal tilt into a PSD-limited PDF.",
                className="transport-sam-controls-intro",
            ),
            dcc.Checklist(
                id=MIXED_ENABLED_ID,
                options=[
                    {
                        "label": " Use SAM WP²RTP + WPRTPTHLP targets",
                        "value": "enabled",
                    }
                ],
                value=[],
                className="transport-sam-coherence-toggle",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Center-direction strength", htmlFor=MIXED_STRENGTH_ID),
                            html.Output("70%", id=MIXED_STRENGTH_VALUE_ID),
                        ],
                        className="transport-sam-control-heading",
                    ),
                    dcc.Slider(
                        id=MIXED_STRENGTH_ID,
                        min=0.0,
                        max=1.0,
                        step=0.05,
                        value=0.70,
                        marks=None,
                        updatemode="drag",
                    ),
                    html.Small(
                        "One local 2-by-2 compatibility solve rotates only the thermodynamic center direction toward WP²RTP and WPRTPTHLP. WPRTP² is an independent check."
                    ),
                ],
                className="transport-sam-control",
            ),
            html.P(
                "This is raw-SAM-only evidence. It preserves the selected mean, full covariance, and center-separation metric length. The ordinary PDF-10 width and internal-covariance diagnosis is then rebuilt at the rotated center direction.",
                className="transport-sam-complexity-note",
            ),
        ],
        className="transport-sam-rail-section",
    )


def _tuning_control(definition):
    name, label, minimum, maximum, step, default, note = definition
    return html.Div(
        [
            html.Div(
                [
                    html.Label(label, htmlFor=_input_id(name)),
                    html.Output(base._format_tuning(name, default), id=_value_id(name)),
                ],
                className="transport-sam-control-heading",
            ),
            dcc.Slider(
                id=_input_id(name),
                min=minimum,
                max=maximum,
                step=step,
                value=default,
                marks=None,
                updatemode="drag",
            ),
            html.Small(note),
        ],
        className="transport-sam-control",
    )


def _control_rail(configuration, source_error: str):
    """Render the original lab's controls with independent Dash IDs."""

    return html.Aside(
        [
            html.Section(
                [
                    html.H2("SAM snapshot"),
                    html.Label("Case", htmlFor=CASE_ID),
                    styled_dropdown(
                        id=CASE_ID,
                        options=[
                            {"label": value["label"], "value": key}
                            for key, value in base.CASE_DEFINITIONS.items()
                        ],
                        value="arm",
                        clearable=False,
                    ),
                    html.Div(
                        [html.Span("LES time"), html.Strong(id=TIME_LABEL_ID)],
                        className="transport-sam-control-heading",
                    ),
                    dcc.Slider(
                        id=TIME_ID,
                        min=0,
                        max=configuration["time_max"],
                        step=1,
                        value=configuration["time_value"],
                        marks=configuration["time_marks"],
                        updatemode="drag",
                    ),
                    html.Div(
                        [html.Span("Altitude"), html.Strong(id=HEIGHT_LABEL_ID)],
                        className="transport-sam-control-heading",
                    ),
                    dcc.Slider(
                        id=HEIGHT_ID,
                        min=0,
                        max=configuration["height_max"],
                        step=1,
                        value=configuration["height_value"],
                        marks=configuration["height_marks"],
                        updatemode="drag",
                    ),
                    html.P(
                        source_error,
                        className="transport-sam-source-error",
                        style={} if source_error else {"display": "none"},
                    ),
                ],
                className="transport-sam-rail-section",
            ),
            html.Section(
                [
                    html.H2("Closure tuning"),
                    html.P(
                        "These remain calibration choices in the local diagnosis. The vertical experiment adjusts only the diagnosed center geometry before reconstructing this selected plane.",
                        className="transport-sam-controls-intro",
                    ),
                    *[_tuning_control(control) for control in base.TUNING_CONTROLS],
                    html.P(
                        "The component-width contrasts are still obtained algebraically from this plane's marginal third moments; the same analytic PSD cap remains active.",
                        className="transport-sam-complexity-note",
                    ),
                ],
                className="transport-sam-rail-section",
            ),
            _coherence_controls(),
            _mixed_moment_controls(),
        ],
        className="transport-sam-control-rail",
    )


def _rt_thl_ellipse(mean: np.ndarray, covariance: np.ndarray, radius: float):
    """Return one r_t--theta_l component footprint in display units."""

    transform = np.diag((1000.0, 1.0))
    projected = transform @ np.asarray(covariance[np.ix_((1, 2), (1, 2))], float) @ transform
    values, vectors = np.linalg.eigh(0.5 * (projected + projected.T))
    angles = np.linspace(0.0, 2.0 * np.pi, 181)
    points = np.asarray((float(mean[1]) * 1000.0, float(mean[2])))[:, None] + radius * vectors @ np.diag(np.sqrt(np.maximum(values, 0.0))) @ np.vstack(
        (np.cos(angles), np.sin(angles))
    )
    return points[0], points[1]


def _rt_thl_figure(snapshot, result):
    """Plot the raw thermodynamic shape that direction-space priors must respect."""

    counts, xedges, yedges = np.histogram2d(
        np.asarray(snapshot.samples[:, 1], float) * 1000.0,
        np.asarray(snapshot.samples[:, 2], float),
        bins=64,
    )
    figure = go.Figure(
        go.Heatmap(
            x=0.5 * (xedges[:-1] + xedges[1:]),
            y=0.5 * (yedges[:-1] + yedges[1:]),
            z=np.log1p(counts.T),
            colorscale="Cividis",
            showscale=False,
            hovertemplate="rₜ=%{x:.3f} g/kg<br>θₗ=%{y:.3f} K<br>raw SAM count=%{customdata:.0f}<extra></extra>",
            customdata=counts.T,
        )
    )
    for number, (mean, covariance, color) in enumerate(
        zip(result.mixture.means, result.mixture.covariances, base.COMPONENT_COLORS),
        start=1,
    ):
        for radius, dash, width in ((1.0, "solid", 1.8), (2.0, "dash", 2.5)):
            x, y = _rt_thl_ellipse(mean, covariance, radius)
            figure.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line={"color": color, "dash": dash, "width": width},
                    name=f"G{number} · {radius:.0f}σ",
                    showlegend=False,
                    hovertemplate=f"G{number} rₜ–θₗ footprint · {radius:.0f}σ<extra></extra>",
                )
            )
        figure.add_trace(
            go.Scatter(
                x=[float(mean[1]) * 1000.0],
                y=[float(mean[2])],
                mode="markers+text",
                marker={"color": color, "size": 8, "symbol": "x"},
                text=[f"G{number}"],
                textposition="top center",
                showlegend=False,
                hovertemplate=(
                    f"G{number}<br>weight={result.mixture.weights[number - 1]:.3f}"
                    "<extra></extra>"
                ),
            )
        )
    figure.update_layout(
        template="plotly_dark",
        height=480,
        margin={"l": 54, "r": 22, "t": 66, "b": 54},
        title=f"Raw-SAM probability with PDF-10 rₜ–θₗ footprints · {snapshot.height_m:.0f} m",
        uirevision="vertical-coherence-rt-thl",
    )
    figure.update_xaxes(title_text="Total water, rₜ [g/kg]")
    figure.update_yaxes(title_text="Liquid-water potential temperature, θₗ [K]")
    return figure


def _lscale_figure(profile, selected_height_m: float):
    """Use PDF-9's source/arrival arrow grammar for direct-Lscale reach."""

    z = np.asarray(profile.height_m, float)
    step = max(1, len(z) // 28)
    figure = go.Figure()
    for reach, color, sign, name, symbol in (
        (profile.lscale_up_m, "#38bdf8", 1.0, "Lscale_up", "arrow-up"),
        (profile.lscale_down_m, "#fb7185", -1.0, "Lscale_down", "arrow-down"),
    ):
        for index in range(0, len(z), step):
            arrival = z[index] + sign * max(float(reach[index]), 0.0)
            figure.add_trace(
                go.Scatter(
                    x=[z[index], z[index]],
                    y=[z[index], arrival],
                    mode="lines+markers",
                    line={"color": color, "width": 2},
                    marker={"size": [4, 8], "symbol": ["square-open", symbol]},
                    name=name,
                    legendgroup=name,
                    showlegend=index == 0,
                    hovertemplate=(
                        "launch=%{x:.0f} m<br>arrival=%{y:.0f} m"
                        f"<br>{name}={abs(arrival-z[index]):.0f} m<extra>{name}</extra>"
                    ),
                )
            )
    figure.add_shape(
        type="line", x0=float(z[0]), y0=float(z[0]), x1=float(z[-1]), y1=float(z[-1]),
        line={"color": "#64748b", "dash": "dot", "width": 1},
    )
    figure.add_hline(y=float(selected_height_m), line_color="#fbbf24", line_dash="dot")
    figure.update_layout(
        template="plotly_dark",
        height=480,
        margin={"l": 54, "r": 22, "t": 66, "b": 54},
        title=(
            "Iterated G1-path direct CLUBB Lscale_up/down · source-to-arrival reach"
            if profile.path_kind == "g1"
            else "Raw-SAM mean-state direct CLUBB Lscale_up/down · source-to-arrival reach"
        ),
        legend={"orientation": "h", "y": -0.18},
        uirevision="vertical-coherence-lscale",
    )
    figure.update_xaxes(title_text="Parcel launch altitude [m]")
    figure.update_yaxes(title_text="Estimated arrival altitude [m]")
    return figure


def _mixed_moment_facts(allocation):
    """Render target, baseline, and adjusted conditional-tail diagnostics."""

    before = np.asarray(allocation.before_normalized, dtype=float)
    after = np.asarray(allocation.after_normalized, dtype=float)
    target = np.asarray(allocation.target_normalized, dtype=float)
    g1, g2 = allocation.result.mixture.covariances

    def _rho(covariance):
        return float(
            covariance[1, 2]
            / np.sqrt(max(float(covariance[1, 1] * covariance[2, 2]), 1.0e-30))
        )

    return fact_grid(
        [
            {
                "label": "WP²RTP · energetic-motion moisture",
                "value": f"SAM {target[0]:+.2f} · 2G {before[0]:+.2f} → {after[0]:+.2f}",
                "detail": "Targeted only when the switch is enabled.",
                "tone": "good" if abs(after[0] - target[0]) <= abs(before[0] - target[0]) else "warning",
            },
            {
                "label": "WPRTPTHLP · transported thermodynamic partition",
                "value": f"SAM {target[2]:+.2f} · 2G {before[2]:+.2f} → {after[2]:+.2f}",
                "detail": "Targeted only when the switch is enabled.",
                "tone": "good" if abs(after[2] - target[2]) <= abs(before[2] - target[2]) else "warning",
            },
            {
                "label": "WPRTP² · independent consistency check",
                "value": f"SAM {target[1]:+.2f} · 2G {before[1]:+.2f} → {after[1]:+.2f}",
                "detail": "Not included in the allocation solve.",
                "tone": "neutral",
            },
            {
                "label": "Mixed center-direction step",
                "value": f"ρ₁={_rho(g1):+.2f} · ρ₂={_rho(g2):+.2f}",
                "detail": allocation.message,
                "tone": "info" if allocation.applied else "neutral",
            },
        ]
    )


def build_layout():
    configuration, source_error = base._initial_slider_configuration()
    return html.Article(
        [
            report_header(
                "Trivariate transport 2G · vertically coherent SAM laboratory",
                "Compare the moment-driven two-Gaussian diagnosis with a bounded fixed point between its G1 center path and CLUBB's native directional mixing-length calculation. Every selected plane retains its supplied mean and covariance.",
                eyebrow="Design A · vertical reach + conditional-tail allocation",
                badges=(
                    "raw LES moments",
                    "native Lscale passes",
                    "local PSD reconstruction",
                    "raw mixed-third targets",
                ),
            ),
            html.Div(
                [
                    html.Section(
                        [
                            html.Div(
                                [
                                    dcc.Graph(
                                        id=FIGURE_IDS[projection],
                                        config={"displaylogo": False, "responsive": True},
                                        className="transport-sam-graph",
                                    )
                                    for projection in base.PROJECTIONS
                                ]
                                + [
                                    dcc.Graph(
                                        id=RT_THL_FIGURE_ID,
                                        config={"displaylogo": False, "responsive": True},
                                        className="transport-sam-graph",
                                    ),
                                    dcc.Graph(
                                        id=LSCALE_FIGURE_ID,
                                        config={"displaylogo": False, "responsive": True},
                                        className="transport-sam-graph",
                                    ),
                                ],
                                className="transport-sam-figure-row",
                            ),
                            html.P(
                                "The heatmaps are raw SAM at the selected plane. Coherence changes only candidate center-separation geometry using directionally reachable source levels. The mixed-moment switch then applies one bounded raw-SAM thermodynamic direction correction before rebuilding the ordinary local PDF-10 allocation; every reconstruction retains its supplied mean and covariance. Gold/cyan mark G1/G2 and solid/dashed loops show 1σ/2σ footprints.",
                                className="transport-sam-figure-key",
                            ),
                            html.Div(id=FACTS_ID),
                            html.P(id=STATUS_ID, className="transport-sam-status"),
                        ],
                        className="transport-sam-visual-panel",
                    ),
                    _control_rail(configuration, source_error),
                ],
                className="transport-sam-workspace",
            ),
            html.Details(
                [
                    html.Summary("What the vertical and mixed-moment switches test"),
                    html.P(
                        "This extends Design A from doc/vertically_coherent_pdf_diagnosis.md without changing CLUBB's closure. Pass zero diagnoses independent PDF-10 mixtures at every raw-SAM level. Each enabled pass reruns native calc_Lscale_directly along the current G1 rₜ–θₗ center curve, gathers compatible source geometry only inside directional reach, blends back toward the frozen local diagnosis, and invokes the ordinary local reconstruction and PSD cap. The iteration count is an experimental numerical control; strong sensitivity after several passes is a warning, not a tuning success.",
                    ),
                    html.P(
                        "The optional mixed-moment experiment is intentionally aimed at the low-level center-placement failure. It reads raw-SAM WP²RTP and WPRTPTHLP, evaluates their two-Gaussian response to small rₜ and θₗ rotations of the center direction, and takes one bounded minimum-norm correction while retaining the same covariance-metric center budget. The ordinary PDF-10 width and internal-covariance allocation is rebuilt afterward. Therefore an improved footprint is evidence that the missing information belongs in center geometry; a weak or counterproductive response is evidence that the existing two-component allocation remains the limit.",
                    ),
                ],
                className="transport-sam-method-note",
            ),
        ],
        className="notes-report transport-sam-report",
    )


def register_callbacks(app):
    tuning_input_ids = tuple(_input_id(control[0]) for control in base.TUNING_CONTROLS)
    tuning_value_ids = tuple(_value_id(control[0]) for control in base.TUNING_CONTROLS)

    @app.callback(
        Output(TIME_ID, "max"),
        Output(TIME_ID, "marks"),
        Output(TIME_ID, "value"),
        Output(HEIGHT_ID, "max"),
        Output(HEIGHT_ID, "marks"),
        Output(HEIGHT_ID, "value"),
        Input(CASE_ID, "value"),
    )
    def change_case(case_name):
        try:
            configuration = base._slider_configuration(case_name or "arm")
        except Exception:
            return 0, {}, 0, 0, {}, 0
        return (
            configuration["time_max"],
            configuration["time_marks"],
            configuration["time_value"],
            configuration["height_max"],
            configuration["height_marks"],
            configuration["height_value"],
        )

    @app.callback(
        Output(TIME_LABEL_ID, "children"),
        Output(HEIGHT_LABEL_ID, "children"),
        *(Output(FIGURE_IDS[projection], "figure") for projection in base.PROJECTIONS),
        Output(RT_THL_FIGURE_ID, "figure"),
        Output(LSCALE_FIGURE_ID, "figure"),
        Output(FACTS_ID, "children"),
        Output(STATUS_ID, "children"),
        Output(BLEND_VALUE_ID, "children"),
        Output(ITERATIONS_VALUE_ID, "children"),
        Output(MIXED_STRENGTH_VALUE_ID, "children"),
        *(Output(identifier, "children") for identifier in tuning_value_ids),
        Input(CASE_ID, "value"),
        Input(TIME_ID, "value"),
        Input(HEIGHT_ID, "value"),
        Input(ENABLED_ID, "value"),
        Input(BLEND_ID, "value"),
        Input(ITERATIONS_ID, "value"),
        Input(MIXED_ENABLED_ID, "value"),
        Input(MIXED_STRENGTH_ID, "value"),
        *(Input(identifier, "value") for identifier in tuning_input_ids),
    )
    def render(
        case_name,
        time_index,
        height_index,
        enabled_values,
        max_blend,
        iterations,
        mixed_enabled_values,
        mixed_strength,
        *values,
    ):
        case_name = case_name if case_name in base.CASE_DEFINITIONS else "arm"
        settings = tuple(float(value) for value in values)
        max_blend = float(max_blend if max_blend is not None else 0.15)
        iterations = int(np.clip(iterations if iterations is not None else 2, 1, 5))
        mixed_strength = float(np.clip(mixed_strength if mixed_strength is not None else 0.70, 0.0, 1.0))
        try:
            configuration = base._slider_configuration(case_name)
            time_index = int(np.clip(time_index or 0, 0, configuration["time_max"]))
            height_index = int(
                np.clip(height_index or 0, 0, configuration["height_max"])
            )
            time_seconds = int(configuration["times"][time_index])
            heights = configuration["heights"]
            snapshot = base._raw_snapshot(
                case_name, time_seconds, float(heights[height_index])
            )
            tuning = _tuning_from_values(settings)
            definition = base.CASE_DEFINITIONS[case_name]
            profile_moments = load_profile_moments(
                str(definition["run_dir"]), time_seconds
            )
            profile_heights = np.asarray(profile_moments.height_m, dtype=float)
            profile_index = int(
                np.argmin(np.abs(profile_heights - float(snapshot.height_m)))
            )
            local_results = tuple(
                diagnose_transport_2g_from_moments(
                    profile_moments.mean[level],
                    profile_moments.covariance[level],
                    profile_moments.third[level],
                    tuning,
                )
                for level in range(profile_heights.size)
            )
            lscale_error = ""
            mean_lscale_profile = None
            try:
                mean_lscale_profile = compute_sam_lscale_profile(
                    str(definition["run_dir"]), time_seconds
                )
            except SamLscaleError as error:
                lscale_error = str(error)

            reach_state = {
                "profile": mean_lscale_profile,
                "fallback": "",
            }

            def reach_provider(current_results, _iteration):
                plume_rt = np.asarray(
                    [result.mixture.means[0, 1] for result in current_results],
                    dtype=float,
                )
                plume_thl = np.asarray(
                    [result.mixture.means[0, 2] for result in current_results],
                    dtype=float,
                )
                try:
                    profile = compute_sam_plume_lscale_profile(
                        str(definition["run_dir"]),
                        time_seconds,
                        plume_rt=plume_rt,
                        plume_thl=plume_thl,
                    )
                    reach_state["profile"] = profile
                    return profile.lscale_up_m, profile.lscale_down_m
                except SamLscaleError as error:
                    reach_state["fallback"] = str(error)
                    if mean_lscale_profile is not None:
                        reach_state["profile"] = mean_lscale_profile
                        return (
                            mean_lscale_profile.lscale_up_m,
                            mean_lscale_profile.lscale_down_m,
                        )
                    zeros = np.zeros(profile_heights.size, dtype=float)
                    return zeros, zeros

            enabled = "enabled" in (enabled_values or [])
            coherence = apply_iterative_vertical_coherence_column(
                profile_moments.mean,
                profile_moments.covariance,
                profile_moments.third,
                profile_heights,
                tuning,
                local_results,
                reach_provider=reach_provider,
                settings=VerticalCoherenceSettings(
                    enabled=enabled,
                    max_blend=max_blend,
                    iterations=iterations,
                ),
            )
            mixed_enabled = "enabled" in (mixed_enabled_values or [])
            allocation = diagnose_transport_2g_with_mixed_center_alignment(
                profile_moments.mean[profile_index],
                profile_moments.covariance[profile_index],
                profile_moments.third[profile_index],
                tuning,
                profile_moments.mixed_third[profile_index],
                mixed_strength if mixed_enabled else 0.0,
                standardized_displacement_override=standardized_displacement(
                    coherence.results[profile_index]
                ),
            )
            result = allocation.result
            state = "on" if enabled else "off"
            status = f"Vertical coherence {state}. {coherence.message} " + base._status(
                snapshot, result
            )
            status += " " + allocation.message
            lscale_profile = reach_state["profile"]
            if lscale_profile is not None:
                local_lscale = float(lscale_profile.lscale_m[profile_index])
                selected_blend = float(coherence.level_blends[profile_index])
                lower_support = float(coherence.lower_support[profile_index])
                upper_support = float(coherence.upper_support[profile_index])
                status += (
                    f" Selected-level direct Lscale={local_lscale:.0f} m; final "
                    f"blend={selected_blend:.0%}, lower/upper aggregate support "
                    f"{lower_support:.2f}/{upper_support:.2f}."
                )
                if reach_state["fallback"]:
                    status += (
                        " G1-path recalculation fell back to mean-state reach: "
                        f"{reach_state['fallback']}"
                    )
                lscale_figure = _lscale_figure(lscale_profile, snapshot.height_m)
            else:
                status += f" Direct Lscale unavailable: {lscale_error}"
                lscale_figure = base._empty_figure("Direct Lscale", lscale_error)
            return (
                f"{snapshot.elapsed_minutes:.0f} min",
                f"{snapshot.height_m:.0f} m",
                *(
                    base.make_projection_figure(snapshot, result, projection)
                    for projection in base.PROJECTIONS
                ),
                _rt_thl_figure(snapshot, result),
                lscale_figure,
                html.Div(
                    [
                        base._facts(snapshot, result),
                        _mixed_moment_facts(allocation),
                        html.P(
                            lscale_profile.approximation if lscale_profile is not None else lscale_error,
                            className="transport-sam-complexity-note",
                        ),
                    ]
                ),
                status,
                f"{max_blend:.0%}",
                str(iterations),
                f"{mixed_strength:.0%}",
                *(
                    base._format_tuning(name, value)
                    for (name, *_rest), value in zip(base.TUNING_CONTROLS, settings)
                ),
            )
        except Exception as error:
            figure = base._empty_figure("Vertically coherent SAM laboratory", str(error))
            return (
                "—",
                "—",
                *(figure for _ in base.PROJECTIONS),
                figure,
                figure,
                empty_state("SAM data unavailable", str(error)),
                "The local vertical-coherence experiment has no raw-SAM plane to diagnose.",
                f"{max_blend:.0%}",
                str(iterations),
                f"{mixed_strength:.0%}",
                *(
                    base._format_tuning(name, value)
                    for (name, *_rest), value in zip(base.TUNING_CONTROLS, settings)
                ),
            )


REPORT = register_report(
    ReportSpec(
        slug="vertically-coherent-transport-sam-lab",
        title="Vertically coherent transport 2G · SAM laboratory",
        summary="Test bounded fixed-point coupling between PDF-10 G1 placement and native directional Lscale reach.",
        category="PDF development",
        updated="Iterative Design A prototype",
        tags=("transport 2G", "SAM", "vertical coherence", "Lscale", "fixed point"),
        order=36,
        build_layout=build_layout,
        register_callbacks=register_callbacks,
    )
)

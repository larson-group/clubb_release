"""Raw-SAM placement laboratory for the moment-driven transport 2G prototype."""

from __future__ import annotations

from functools import lru_cache
from math import erf, exp, pi, sqrt
from pathlib import Path

from dash import Input, Output, dcc, html
import numpy as np
import plotly.graph_objects as go

if __package__ and __package__.startswith("dash_app."):
    from dash_app.shared.bivariate_heatmap import (
        BIVARIATE_LEVELS,
        signed_bivariate_codes,
        signed_bivariate_reference_colorscale,
    )
    from dash_app.shared.reporting import empty_state, fact_grid, report_header
else:  # Script-style imports used by dash_app/app.py.
    from shared.bivariate_heatmap import (
        BIVARIATE_LEVELS,
        signed_bivariate_codes,
        signed_bivariate_reference_colorscale,
    )
    from shared.reporting import empty_state, fact_grid, report_header

from utilities.sam_3d_reference import (
    DEFAULT_BOMEX_SAM_RUN,
    DEFAULT_SAM_RUN,
    diagnose_cloud,
    inventory_run,
    load_snapshot,
)

from ..registry import ReportSpec, register_report
from ..transport_2g_prototype import (
    Transport2GTuning,
    Transport2GResult,
    diagnose_transport_2g_from_moments,
)


CASE_ID = "notes-transport-sam-case"
TIME_ID = "notes-transport-sam-time"
HEIGHT_ID = "notes-transport-sam-height"
TIME_LABEL_ID = "notes-transport-sam-time-label"
HEIGHT_LABEL_ID = "notes-transport-sam-height-label"
FACTS_ID = "notes-transport-sam-facts"
STATUS_ID = "notes-transport-sam-status"
FIGURE_IDS = {
    "rt": "notes-transport-sam-w-rt",
    "chi": "notes-transport-sam-w-chi",
    "rc": "notes-transport-sam-w-rc",
}

TUNING_CONTROLS = (
    (
        "moisture-tail-gain",
        "Moisture-tail weight gain",
        0.0,
        2.5,
        0.05,
        1.00,
        "Maps |rₜ skewness| to G1's diagnosed mass. Higher values make the transport component rarer.",
    ),
    (
        "center-budget",
        "Center-separation budget",
        0.10,
        0.92,
        0.01,
        0.72,
        "The allowed covariance-metric room used by the diagnosed center vector. The PSD bound remains analytic.",
    ),
    (
        "w-direction-scale",
        "Mean-w direction scale",
        0.0,
        0.70,
        0.01,
        0.20,
        "Keeps the two mean w values close while w skewness can instead change their internal widths.",
    ),
    (
        "g1-wrt-capture",
        "G1 positive w–rₜ tilt capture",
        0.0,
        1.0,
        0.01,
        1.00,
        "For positive grid w–rₜ covariance, requests that residual fraction in G1. Negative grid covariance is kept out of G1.",
    ),
    (
        "plume-structure-strength",
        "Coherent plume-structure strength",
        0.0,
        1.0,
        0.01,
        0.50,
        "For a positive-w, positive-w–rₜ plume, blends placement toward a signed transport direction and asks G2 to remain comparatively well mixed.",
    ),
)

CASE_DEFINITIONS = {
    "arm": {
        "label": "ARM",
        "run_dir": Path(DEFAULT_SAM_RUN),
        "default_minutes": 453.0,
        "default_height": 1300.0,
    },
    "bomex": {
        "label": "BOMEX",
        "run_dir": Path(DEFAULT_BOMEX_SAM_RUN),
        "default_minutes": 180.0,
        "default_height": 800.0,
    },
}
PROJECTIONS = ("rt", "chi", "rc")
PROJECTION_DETAILS = {
    "rt": ("w–rₜ", "Total water, rₜ [g/kg]"),
    "chi": ("w–χ", "Extended cloud water, χ [g/kg]"),
    "rc": ("w–r_c", "Cloud water, r_c [g/kg]"),
}
COMPONENT_COLORS = ("#fbbf24", "#22d3ee")


def _input_id(name: str) -> str:
    return f"notes-transport-sam-{name}"


def _value_id(name: str) -> str:
    return f"notes-transport-sam-{name}-value"


def _format_tuning(name: str, value: float) -> str:
    if name in {"center-budget", "g1-wrt-capture", "plume-structure-strength"}:
        return f"{value:.0%}"
    return f"{value:.2f}"


def _nearest_index(values, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


def _slider_marks(values, formatter) -> dict[int, str]:
    values = np.asarray(values, dtype=float)
    indices = np.unique(
        np.rint(np.linspace(0, len(values) - 1, min(6, len(values)))).astype(int)
    )
    return {int(index): formatter(float(values[index])) for index in indices}


@lru_cache(maxsize=2)
def _inventory(case_name: str):
    definition = CASE_DEFINITIONS.get(case_name, CASE_DEFINITIONS["arm"])
    return inventory_run(definition["run_dir"])


@lru_cache(maxsize=32)
def _raw_snapshot(case_name: str, seconds: int, height_m: float):
    definition = CASE_DEFINITIONS.get(case_name, CASE_DEFINITIONS["arm"])
    return load_snapshot(definition["run_dir"], int(seconds), float(height_m))


def _slider_configuration(case_name: str):
    """Return the time/height controls for one raw-SAM run."""

    case_name = case_name if case_name in CASE_DEFINITIONS else "arm"
    inventory = _inventory(case_name)
    definition = CASE_DEFINITIONS[case_name]
    times = np.asarray(inventory.steps_seconds, dtype=float)
    heights = np.asarray(inventory.heights_m, dtype=float)
    time_value = _nearest_index(times, 60.0 * definition["default_minutes"])
    height_value = _nearest_index(heights, definition["default_height"])
    return {
        "time_max": len(times) - 1,
        "time_marks": _slider_marks(times, lambda value: f"{value / 60.0:.0f}m"),
        "time_value": time_value,
        "height_max": len(heights) - 1,
        "height_marks": _slider_marks(heights, lambda value: f"{value:.0f}m"),
        "height_value": height_value,
        "times": times,
        "heights": heights,
    }


def _initial_slider_configuration():
    try:
        return _slider_configuration("arm"), ""
    except Exception as error:  # The report still mounts off the SAM filesystem.
        return (
            {
                "time_max": 0,
                "time_marks": {},
                "time_value": 0,
                "height_max": 0,
                "height_marks": {},
                "height_value": 0,
                "times": np.array((0.0,)),
                "heights": np.array((0.0,)),
            },
            str(error),
        )


def _tuning_from_values(values: tuple[float, ...]) -> Transport2GTuning:
    settings = {
        name: float(value)
        for (name, *_rest), value in zip(TUNING_CONTROLS, values)
    }
    return Transport2GTuning(
        moisture_tail_gain=settings["moisture-tail-gain"],
        center_budget=settings["center-budget"],
        w_direction_scale=settings["w-direction-scale"],
        g1_wrt_capture=settings["g1-wrt-capture"],
        plume_structure_strength=settings["plume-structure-strength"],
    )


def _tuning_control(definition):
    name, label, minimum, maximum, step, default, note = definition
    return html.Div(
        [
            html.Div(
                [
                    html.Label(label, htmlFor=_input_id(name)),
                    html.Output(_format_tuning(name, default), id=_value_id(name)),
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


def _raw_heatmap(snapshot, projection: str) -> go.Heatmap:
    """Use the same raw-SAM probability/signed-transport background as PDF-9."""

    if projection == "chi":
        scalar = np.asarray(snapshot.chi_samples, dtype=float) * 1000.0
        symbol = "χ"
    elif projection == "rc":
        scalar = np.asarray(snapshot.rc_samples, dtype=float) * 1000.0
        symbol = "r_c"
    else:
        scalar = np.asarray(snapshot.samples[:, 1], dtype=float) * 1000.0
        symbol = "rₜ"
    w = np.asarray(snapshot.samples[:, 0], dtype=float)
    signed_cloud_transport = (w - float(snapshot.mean[0])) * snapshot.rc_samples * 1000.0
    counts, xedges, yedges = np.histogram2d(scalar, w, bins=56)
    transport, _, _ = np.histogram2d(
        scalar, w, bins=(xedges, yedges), weights=signed_cloud_transport
    )
    codes, _, _ = signed_bivariate_codes(counts.T, transport.T)
    return go.Heatmap(
        x=0.5 * (xedges[:-1] + xedges[1:]),
        y=0.5 * (yedges[:-1] + yedges[1:]),
        z=codes,
        colorscale=signed_bivariate_reference_colorscale(),
        zmin=0,
        zmax=BIVARIATE_LEVELS * (2 * BIVARIATE_LEVELS - 1) - 1,
        showscale=False,
        hovertemplate=(
            f"{symbol}=%{{x:.3f}} g/kg<br>w=%{{y:.3f}} m/s"
            "<extra>raw SAM bin</extra>"
        ),
    )


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + erf(value / sqrt(2.0)))


def _normal_pdf(value: float) -> float:
    return exp(-0.5 * value * value) / sqrt(2.0 * pi)


def _projection_component(snapshot, covariance: np.ndarray, mean: np.ndarray, projection: str):
    """Return one physical w-scalar component projection in display units."""

    if projection == "rt":
        return {
            "mean_x": float(mean[1]) * 1000.0,
            "mean_w": float(mean[0]),
            "var_x": float(covariance[1, 1]) * 1.0e6,
            "var_w": float(covariance[0, 0]),
            "cov_xw": float(covariance[0, 1]) * 1000.0,
            "cloudy_branch": False,
        }

    mean_chi = float(
        snapshot.chi_mean
        + snapshot.chi_coeff_rt * (mean[1] - snapshot.mean[1])
        - snapshot.chi_coeff_thl * (mean[2] - snapshot.mean[2])
    )
    variance_chi = max(
        float(
            snapshot.chi_coeff_rt**2 * covariance[1, 1]
            + snapshot.chi_coeff_thl**2 * covariance[2, 2]
            - 2.0
            * snapshot.chi_coeff_rt
            * snapshot.chi_coeff_thl
            * covariance[1, 2]
        ),
        1.0e-30,
    )
    covariance_w_chi = float(
        snapshot.chi_coeff_rt * covariance[0, 1]
        - snapshot.chi_coeff_thl * covariance[0, 2]
    )
    if projection == "chi":
        return {
            "mean_x": mean_chi * 1000.0,
            "mean_w": float(mean[0]),
            "var_x": variance_chi * 1.0e6,
            "var_w": float(covariance[0, 0]),
            "cov_xw": covariance_w_chi * 1000.0,
            "cloudy_branch": False,
        }

    # r_c is max(chi, 0), hence its continuous cloudy branch is not Gaussian.
    # The plot draws the positive part of the Gaussian chi footprint rather
    # than claiming that the r_c marginal is itself an ellipse.
    sigma_chi = sqrt(variance_chi)
    normalized = float(np.clip(mean_chi / sigma_chi, -12.0, 12.0))
    cloudy_probability = _normal_cdf(normalized)
    mean_rc = mean_chi * cloudy_probability + sigma_chi * _normal_pdf(normalized)
    return {
        "mean_x": mean_rc * 1000.0,
        "ellipse_mean_x": mean_chi * 1000.0,
        "mean_w": float(mean[0]),
        "var_x": variance_chi * 1.0e6,
        "var_w": float(covariance[0, 0]),
        "cov_xw": covariance_w_chi * 1000.0,
        "cloudy_branch": True,
        "chi_mean_x": mean_chi * 1000.0,
    }


def _ellipse(component: dict, radius: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    covariance = np.array(
        (
            (component["var_x"], component["cov_xw"]),
            (component["cov_xw"], component["var_w"]),
        ),
        dtype=float,
    )
    values, vectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    transform = vectors @ np.diag(np.sqrt(np.maximum(values, 0.0)))
    angles = np.linspace(0.0, 2.0 * pi, 181)
    points = np.array([[component.get("ellipse_mean_x", component["mean_x"])], [component["mean_w"]]]) + radius * transform @ np.vstack(
        (np.cos(angles), np.sin(angles))
    )
    return points[0], points[1]


def _component_correlation(covariance: np.ndarray) -> float:
    return float(
        covariance[0, 1]
        / sqrt(max(float(covariance[0, 0] * covariance[1, 1]), 1.0e-30))
    )


def make_projection_figure(snapshot, result: Transport2GResult, projection: str) -> go.Figure:
    """Overlay both diagnosed components on one selected raw-SAM plane."""

    title, x_title = PROJECTION_DETAILS[projection]
    figure = go.Figure(_raw_heatmap(snapshot, projection))
    for number, (mean, covariance, color) in enumerate(
        zip(
            result.mixture.means,
            result.mixture.covariances,
            COMPONENT_COLORS,
        ),
        start=1,
    ):
        component = _projection_component(snapshot, covariance, mean, projection)
        for radius, width, dash in ((1.0, 1.6, "solid"), (2.0, 2.8, "dash")):
            x, y = _ellipse(component, radius)
            if component["cloudy_branch"]:
                # Map only the physically cloudy chi portion into the r_c
                # panel.  A broken line is more honest than a fake r_c
                # Gaussian ellipse, while still exposing the component shape.
                x = np.where(x >= 0.0, x, np.nan)
                label = f"G{number} cloudy χ branch · {radius:.0f}σ"
            else:
                label = f"G{number} · {radius:.0f}σ"
            figure.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line={"color": color, "width": width, "dash": dash},
                    name=label,
                    showlegend=False,
                    hovertemplate=(
                        f"{label}<br>weight={result.mixture.weights[number - 1]:.3f}"
                        "<extra></extra>"
                    ),
                )
            )
        figure.add_trace(
            go.Scatter(
                x=[component["mean_x"]],
                y=[component["mean_w"]],
                mode="markers+text",
                text=[f"G{number}"],
                textposition="top center",
                marker={
                    "color": color,
                    "size": 11,
                    "symbol": "x",
                    "line": {"color": "#0f172a", "width": 1.2},
                },
                name=f"G{number} center",
                showlegend=False,
                hovertemplate=(
                    f"G{number} center<br>{x_title.split(' [')[0]}=%{{x:.3f}}"
                    "<br>w=%{y:.3f}<extra></extra>"
                ),
            )
        )

    figure.add_hline(
        y=float(snapshot.mean[0]),
        line={"color": "#94a3b8", "width": 1.0, "dash": "dot"},
    )
    if projection == "rt":
        figure.add_vline(
            x=float(snapshot.mean[1]) * 1000.0,
            line={"color": "#94a3b8", "width": 1.0, "dash": "dot"},
        )
    else:
        figure.add_vline(
            x=0.0,
            line={"color": "#f59e0b", "width": 1.3, "dash": "dot"},
        )
    suffix = "cloudy χ branches" if projection == "rc" else "Gaussian footprints"
    figure.update_xaxes(title_text=x_title, zeroline=False)
    figure.update_yaxes(title_text="Vertical velocity, w [m/s]", zeroline=False)
    figure.update_layout(
        template="plotly_dark",
        title={
            "text": f"{title} · raw SAM + {suffix}",
            "x": 0.5,
            "font": {"size": 15},
        },
        height=410,
        margin={"l": 58, "r": 14, "t": 48, "b": 54},
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        showlegend=False,
        uirevision=(
            f"notes-transport-sam-{projection}-{snapshot.run_dir}-"
            f"{snapshot.elapsed_seconds}-{snapshot.height_m:.0f}"
        ),
    )
    return figure


def _empty_figure(title: str, message: str) -> go.Figure:
    figure = go.Figure()
    figure.add_annotation(
        text=message,
        showarrow=False,
        font={"color": "#94a3b8", "size": 14},
    )
    figure.update_layout(
        template="plotly_dark",
        title=title,
        height=410,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
    )
    return figure


def _wrt_covariance_budget(result: Transport2GResult) -> tuple[float, float, float, float]:
    weights = result.mixture.weights
    center = weights[0] * weights[1] * result.displacement[0] * result.displacement[1]
    g1 = weights[0] * result.mixture.covariances[0, 0, 1]
    g2 = weights[1] * result.mixture.covariances[1, 0, 1]
    return center, g1, g2, float(result.target_covariance[0, 1])


def _facts(snapshot, result: Transport2GResult):
    g1, g2 = result.mixture.covariances
    center, g1_wrt, g2_wrt, total = _wrt_covariance_budget(result)
    standard_deviations = np.sqrt(np.maximum(np.diag(result.target_covariance), 1.0e-30))
    third_error = np.max(
        np.abs(result.reconstructed_third - result.target_third) / standard_deviations**3
    )
    cloud = diagnose_cloud(
        result.mixture.weights,
        result.mixture.means,
        result.mixture.covariances,
        snapshot,
    )
    g1_offset = result.mixture.means[0] - snapshot.mean
    return fact_grid(
        [
            {
                "label": "Diagnosed G1 mass",
                "value": f"{result.mixture.weights[0]:.1%}",
                "detail": (
                    f"G1 anomaly: Δw={g1_offset[0]:+.2f} m/s, "
                    f"Δrₜ={g1_offset[1] * 1000.0:+.2f} g/kg"
                ),
                "tone": "good" if result.mixture.weights[0] < 0.25 else "info",
            },
            {
                "label": "Center covariance budget",
                "value": f"{result.center_metric_fraction:.0%}",
                "detail": (
                    f"G1–G2: Δw={result.displacement[0]:+.2f} m/s, "
                    f"Δrₜ={result.displacement[1] * 1000.0:+.2f} g/kg"
                ),
                "tone": "info",
            },
            {
                "label": "Internal w–rₜ tilt",
                "value": f"ρ₁={_component_correlation(g1):+.2f} · ρ₂={_component_correlation(g2):+.2f}",
                "detail": (
                    "covariance [m/s·g/kg]: "
                    f"center {center * 1000.0:+.3f}; G1 {g1_wrt * 1000.0:+.3f}; "
                    f"G2 {g2_wrt * 1000.0:+.3f}; total {total * 1000.0:+.3f}"
                ),
                "tone": "good" if g1[0, 1] >= -1.0e-12 else "warning",
            },
            {
                "label": "Coherent plume branch",
                "value": f"active blend {result.plume_blend:.0%}",
                "detail": (
                    "Positive w skewness and w–rₜ correlation open the gate; "
                    "the slider supplies its maximum strength."
                ),
                "tone": "info" if result.plume_blend > 1.0e-6 else "neutral",
            },
            {
                "label": "Marginal-third fit",
                "value": f"max |ΔSk| = {third_error:.2f}",
                "detail": (
                    f"width/tilt contrast used {result.contrast_scale:.0%}; "
                    "mean and covariance are conserved to roundoff"
                ),
                "tone": "good" if third_error < 0.05 else "warning",
            },
            {
                "label": "Cloud fraction (evaluation)",
                "value": f"SAM {snapshot.cloud_fraction:.1%} · 2G {cloud['cloud_fraction']:.1%}",
                "detail": "Linearized χ cloud integral; not supplied to the closure diagnosis.",
                "tone": "neutral",
            },
        ]
    )


def _status(snapshot, result: Transport2GResult) -> str:
    covariance_error = float(
        np.max(np.abs(result.reconstructed_covariance - result.target_covariance))
    )
    messages = [
        "Inputs are instantaneous raw-SAM plane moments: the full w/rₜ/θₗ mean and covariance plus w′³, rₜ′³, and θₗ′³. ",
        "The backgrounds are the same plane: gold is probability, blue/red is positive/negative signed cloud transport. ",
        f"Mean/covariance reconstruction error is {covariance_error:.2e}. ",
        "Gold is G1; cyan is G2. Solid and dashed loops are 1σ and 2σ; the r_c panel shows only the positive χ branch because r_c is nonlinear. ",
    ]
    if result.flags:
        messages.append(" ".join(result.flags))
    else:
        messages.append("No prototype PSD cap was active at this selection.")
    return "".join(messages)


def _control_rail(configuration, source_error: str):
    return html.Aside(
        [
            html.Section(
                [
                    html.H2("SAM snapshot"),
                    html.Label("Case", htmlFor=CASE_ID),
                    dcc.Dropdown(
                        id=CASE_ID,
                        options=[
                            {"label": value["label"], "value": key}
                            for key, value in CASE_DEFINITIONS.items()
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
                        "These are calibration choices in the diagnosis, not direct component-center sliders. Every remaining center, width, and correlation is diagnosed from the selected LES moments.",
                        className="transport-sam-controls-intro",
                    ),
                    *[_tuning_control(control) for control in TUNING_CONTROLS],
                    html.P(
                        "The three component-width contrasts are then obtained algebraically from the three diagonal third moments. Any nonphysical request is softened by one analytic PSD cap, shown beside the plots.",
                        className="transport-sam-complexity-note",
                    ),
                ],
                className="transport-sam-rail-section",
            ),
        ],
        className="transport-sam-control-rail",
    )


def build_layout():
    configuration, source_error = _initial_slider_configuration()
    return html.Article(
        [
            report_header(
                "Trivariate transport 2G · SAM placement laboratory",
                "Diagnose the proposed two-Gaussian transport PDF directly from one raw LES plane, then inspect where its two components land over w–rₜ, w–χ, and w–r_c SAM heatmaps.",
                eyebrow="Moment-driven closure prototype",
                badges=(
                    "raw LES moments",
                    "analytic PSD caps",
                    "three projections",
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
                                    for projection in PROJECTIONS
                                ],
                                className="transport-sam-figure-row",
                            ),
                            html.P(
                                "The plotted heatmaps are raw SAM. Gold/cyan × marks are the diagnosed G1/G2 centers; solid/dashed loops are their 1σ/2σ projections. The r_c view intentionally clips to the cloudy χ branch rather than pretending r_c is Gaussian. No legend is embedded inside the plots, so it cannot consume their vertical layout.",
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
                    html.Summary("What this test does — and does not — establish"),
                    html.P(
                        "This is a Python prototype of the new trivariate transport closure, not a fit to the plotted heatmap. The same plane's mean, covariance, and three marginal third moments are the only inputs. Its two Gaussian centers are mean-preserving, its covariance is conserved exactly, and its width contrasts use the closed-form marginal-third relation. The mixed third moments and cloud quantities shown here remain held-out evaluation signals. A good overlay demonstrates whether the closure family has the necessary placement degrees of freedom; it does not yet validate a prognostic CLUBB implementation.",
                    ),
                ],
                className="transport-sam-method-note",
            ),
        ],
        className="notes-report transport-sam-report",
    )


def register_callbacks(app):
    tuning_input_ids = tuple(_input_id(control[0]) for control in TUNING_CONTROLS)
    tuning_value_ids = tuple(_value_id(control[0]) for control in TUNING_CONTROLS)

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
            configuration = _slider_configuration(case_name or "arm")
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
        *(Output(FIGURE_IDS[projection], "figure") for projection in PROJECTIONS),
        Output(FACTS_ID, "children"),
        Output(STATUS_ID, "children"),
        *(Output(identifier, "children") for identifier in tuning_value_ids),
        Input(CASE_ID, "value"),
        Input(TIME_ID, "value"),
        Input(HEIGHT_ID, "value"),
        *(Input(identifier, "value") for identifier in tuning_input_ids),
    )
    def render(case_name, time_index, height_index, *values):
        case_name = case_name if case_name in CASE_DEFINITIONS else "arm"
        settings = tuple(float(value) for value in values)
        try:
            configuration = _slider_configuration(case_name)
            time_index = int(
                np.clip(time_index or 0, 0, configuration["time_max"])
            )
            height_index = int(
                np.clip(height_index or 0, 0, configuration["height_max"])
            )
            snapshot = _raw_snapshot(
                case_name,
                int(configuration["times"][time_index]),
                float(configuration["heights"][height_index]),
            )
            result = diagnose_transport_2g_from_moments(
                snapshot.mean,
                snapshot.covariance,
                np.array((snapshot.wp3, snapshot.rtp3, snapshot.thlp3), dtype=float),
                _tuning_from_values(settings),
            )
            return (
                f"{snapshot.elapsed_minutes:.0f} min",
                f"{snapshot.height_m:.0f} m",
                *(make_projection_figure(snapshot, result, projection) for projection in PROJECTIONS),
                _facts(snapshot, result),
                _status(snapshot, result),
                *(
                    _format_tuning(name, value)
                    for (name, *_rest), value in zip(TUNING_CONTROLS, settings)
                ),
            )
        except Exception as error:
            figure = _empty_figure("SAM placement laboratory", str(error))
            return (
                "—",
                "—",
                *(figure for _ in PROJECTIONS),
                empty_state("SAM data unavailable", str(error)),
                "The closure prototype has no raw-SAM plane to diagnose.",
                *(
                    _format_tuning(name, value)
                    for (name, *_rest), value in zip(TUNING_CONTROLS, settings)
                ),
            )


REPORT = register_report(
    ReportSpec(
        slug="trivariate-transport-sam-lab",
        title="Trivariate transport 2G · SAM placement laboratory",
        summary="Diagnose the proposed two-Gaussian transport closure from raw LES moments and overlay it on three SAM heatmaps.",
        category="PDF development",
        updated="Interactive prototype",
        tags=("transport 2G", "SAM", "w–rₜ", "trivariate"),
        order=35,
        build_layout=build_layout,
        register_callbacks=register_callbacks,
    )
)

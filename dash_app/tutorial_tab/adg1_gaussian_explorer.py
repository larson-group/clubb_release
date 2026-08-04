"""Small interactive view of the two Gaussians diagnosed by ADG1.

This is intentionally a visualization of the normalized ADG1 teaching
diagnosis. It is not a replacement for the full CLUBB Fortran PDF diagnosis.
"""

from __future__ import annotations

from math import pi, sqrt

from dash import Input, Output, dcc, html
import numpy as np
import plotly.graph_objects as go

from dash_app.shared.reporting import report_header, report_section
from dash_app.shared.two_gaussian_model import (
    Mixture,
    diagnose_adg1,
    grid_covariance,
    projected_components,
)

FIGURE_ID = "notes-adg1-gaussian-figure"
FIGURE_IDS = {
    "w_rt": FIGURE_ID,
    "w_chi": "notes-adg1-w-chi-figure",
    "w_rc": "notes-adg1-w-rc-figure",
}
STATUS_ID = "notes-adg1-gaussian-status"
SUMMARY_ID = "notes-adg1-component-summary"
RESET_ID = "notes-adg1-reset"
ADG1_CONTROL_GROUP_ID = "notes-adg1-adg-controls"
CHI_MEAN = -0.25

DEFAULTS = {
    "var-w": 1.00,
    "var-rt": 1.00,
    "var-thl": 0.80,
    "corr-w-rt": 0.55,
    "corr-w-thl": -0.35,
    "corr-rt-thl": -0.15,
    "skew-w": 1.10,
    "gamma": 0.40,
    "beta": 1.50,
}

COMMON_CONTROLS = (
    ("var-w", "w variance", 0.20, 2.00, 0.05, "Sets the vertical spread."),
    ("var-rt", "rₜ variance", 0.20, 2.00, 0.05, "Sets the horizontal spread."),
    ("var-thl", "θₗ variance", 0.20, 2.00, 0.05, "Completes the shared trivariate covariance."),
    ("corr-w-rt", "corr(w, rₜ)", -0.95, 0.95, 0.01, "Places moist and rising air together or apart."),
    ("corr-w-thl", "corr(w, θₗ)", -0.95, 0.95, 0.01, "Sets the shared w–θₗ covariance target."),
    ("corr-rt-thl", "corr(rₜ, θₗ)", -0.95, 0.95, 0.01, "Sets the shared rₜ–θₗ covariance target."),
)
ADG1_CONTROLS = (
    ("skew-w", "w skewness", -2.40, 2.40, 0.05, "Sets the component weights and vertical separation."),
    ("gamma", "γ: within-Gaussian w-width ceiling", 0.05, 0.95, 0.01, "Actual within-Gaussian w variance is γ[1 − max(corr(w,rₜ)², corr(w,θₗ)²)]. Higher γ broadens both Gaussians and reduces their center separation."),
    ("beta", "β: scalar-width allocation", 0.10, 2.90, 0.05, "Redistributes the remaining rₜ and θₗ variance between G1 and G2 while preserving the grid variance. It changes G1's available room for internal tilt."),
)
CONTROLS = COMMON_CONTROLS + ADG1_CONTROLS
COMPONENT_COLORS = ("#fbbf24", "#22d3ee")
PROJECTIONS = ("w_rt", "w_chi", "w_rc")
PROJECTION_DETAILS = {
    "w_rt": ("w–rₜ", "Normalized total-water anomaly, rₜ", "rₜ", -5.0, 5.0),
    "w_chi": ("w–χ", "Normalized extended cloud water, χ", "χ", -5.0, 5.0),
    "w_rc": ("w–rᶜ", "Normalized cloud water, rᶜ", "rᶜ", 0.0, 5.0),
}


def _input_id(name: str) -> str:
    return f"notes-adg1-{name}"


def _value_id(name: str) -> str:
    return f"notes-adg1-{name}-value"


def _format_value(name: str, value: float) -> str:
    if name == "transport-weight" or name.endswith("variance-share"):
        return f"{value:.0%}"
    return f"{value:+.2f}" if "corr" in name or name.startswith("delta-") or name == "skew-w" else f"{value:.2f}"


def _control(name: str, label: str, minimum: float, maximum: float, step: float, note: str):
    value = DEFAULTS[name]
    return html.Div(
        [
            html.Div(
                [
                    html.Label(label, htmlFor=_input_id(name)),
                    html.Output(_format_value(name, value), id=_value_id(name)),
                ],
                className="adg1-control-heading",
            ),
            dcc.Slider(
                id=_input_id(name),
                min=minimum,
                max=maximum,
                step=step,
                value=value,
                marks=None,
                updatemode="drag",
            ),
            html.Small(note),
        ],
        className="adg1-control",
    )


def _settings(values: tuple[float, ...]) -> dict[str, float]:
    return {control[0]: float(value) for control, value in zip(CONTROLS, values)}


def _is_psd(matrix: np.ndarray, tolerance: float = 1.0e-8) -> bool:
    return bool(np.min(np.linalg.eigvalsh(0.5 * (matrix + matrix.T))) >= -tolerance)


def _adg_weight(skew_w: float, within_w: float) -> float:
    if abs(skew_w) <= 1.0e-7:
        return 0.5
    denominator = sqrt(max(4.0 * (1.0 - within_w) ** 3 + skew_w**2, 1.0e-12))
    return float(np.clip(0.5 * (1.0 - skew_w / denominator), 0.04, 0.96))


def _between_covariance(weight: float, displacement: np.ndarray) -> np.ndarray:
    return weight * (1.0 - weight) * np.outer(displacement, displacement)


def _scale_displacement_to_psd(
    covariance: np.ndarray, weight: float, displacement: np.ndarray
) -> tuple[np.ndarray, float]:
    scale = 1.0
    for _ in range(80):
        candidate = displacement * scale
        residual = covariance - _between_covariance(weight, candidate)
        floor = 1.0e-6 * max(float(np.max(np.diag(residual))), 1.0)
        if float(np.min(np.linalg.eigvalsh(residual))) >= floor:
            return candidate, scale
        scale *= 0.96
    return np.zeros_like(displacement), 0.0


def _scalar_variance_factors(weight: float, beta: float) -> tuple[float, float]:
    """Return the CLUBB beta allocation factors for G1 and G2 scalar widths."""

    first = (2.0 / 3.0) * beta + 2.0 * weight * (1.0 - (2.0 / 3.0) * beta)
    return float(first), float(2.0 - first)


def _component_covariances(
    residual: np.ndarray, weight: float, beta: float, g1_w_rt_covariance: float = 0.0
) -> np.ndarray | None:
    """Allocate residual scalar widths with CLUBB's beta parameterization."""

    complement = 1.0 - weight
    first_factor, second_factor = _scalar_variance_factors(weight, beta)
    if min(first_factor, second_factor, *np.diag(residual)) <= 1.0e-10:
        return None
    components = np.zeros((2, 3, 3), dtype=float)
    components[:, 0, 0] = residual[0, 0]
    for variable in (1, 2):
        components[0, variable, variable] = first_factor * residual[variable, variable] / (2.0 * weight)
        components[1, variable, variable] = second_factor * residual[variable, variable] / (2.0 * complement)
    denominator = (
        weight * sqrt(components[0, 1, 1] * components[0, 2, 2])
        + complement * sqrt(components[1, 1, 1] * components[1, 2, 2])
    )
    scalar_rho = float(residual[1, 2] / max(denominator, 1.0e-12))
    if abs(scalar_rho) >= 0.999:
        return None
    for component in range(2):
        scalar_covariance = scalar_rho * sqrt(
            components[component, 1, 1] * components[component, 2, 2]
        )
        components[component, 1, 2] = scalar_covariance
        components[component, 2, 1] = scalar_covariance
    components[0, 0, 1] = components[0, 1, 0] = g1_w_rt_covariance
    return components if all(_is_psd(component) for component in components) else None


def diagnose_parameterized_adg1(
    covariance: np.ndarray, skew_w: float, gamma: float, beta: float
) -> Mixture:
    """Normalized ADG1 reconstruction with active CLUBB gamma and beta controls."""

    variances = np.maximum(np.diag(covariance), 1.0e-12)
    correlations = covariance / np.sqrt(np.outer(variances, variances))
    strongest_w_scalar_correlation_sqd = min(
        max(float(correlations[0, 1] ** 2), float(correlations[0, 2] ** 2)), 1.0
    )
    within_w = float(np.clip(gamma, 0.0, 1.0)) * (1.0 - strongest_w_scalar_correlation_sqd)
    within_w = float(np.clip(within_w, 1.0e-4, 0.99))
    beta = float(np.clip(beta, 1.0e-4, 2.9999))
    weight = _adg_weight(skew_w, within_w)
    between_w = max((1.0 - within_w) * covariance[0, 0], 1.0e-10)
    displacement_w = sqrt(between_w / (weight * (1.0 - weight)))
    displacement = np.array(
        [
            displacement_w,
            covariance[0, 1] / max(weight * (1.0 - weight) * displacement_w, 1.0e-12),
            covariance[0, 2] / max(weight * (1.0 - weight) * displacement_w, 1.0e-12),
        ]
    )
    displacement, scale = _scale_displacement_to_psd(covariance, weight, displacement)
    residual = covariance - _between_covariance(weight, displacement)
    components = _component_covariances(residual, weight, beta)
    if components is None:
        # This is rare for the normal slider range.  The tutorial's fixed
        # ADG1 reconstruction remains a safe diagnostic fallback.
        fallback = diagnose_adg1(covariance, skew_w)
        return Mixture(
            "adg1",
            "ADG1",
            fallback.weights,
            fallback.means,
            fallback.covariances,
            fallback.note,
            "The requested γ/β allocation was not realizable; the explorer used its baseline ADG1 allocation.",
        )
    adjustment = ""
    if scale < 0.999:
        adjustment = "Component separation was reduced to keep the ADG1 matrices realizable."
    return Mixture(
        "adg1",
        "ADG1",
        np.array([weight, 1.0 - weight]),
        np.vstack(((1.0 - weight) * displacement, -weight * displacement)),
        components,
        "γ sets the resolved within-Gaussian w width; β allocates remaining scalar width between G1 and G2. Each component's w–scalar tilt is zero.",
        adjustment,
    )


def _diagnosis(values: tuple[float, ...]):
    settings = _settings(values)
    covariance, correlation, adjusted = grid_covariance(
        (settings["var-w"], settings["var-rt"], settings["var-thl"]),
        (settings["corr-w-rt"], settings["corr-w-thl"], settings["corr-rt-thl"]),
    )
    mixture = diagnose_parameterized_adg1(
        covariance, settings["skew-w"], settings["gamma"], settings["beta"]
    )
    return mixture, correlation, adjusted


def _density(x_grid: np.ndarray, w_grid: np.ndarray, component: dict) -> np.ndarray:
    """Evaluate one displayed projected Gaussian, with a numerical PSD margin."""

    covariance = np.array(
        [
            [component["var_x"], component["cov_wx"]],
            [component["cov_wx"], component["var_w"]],
        ]
    )
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    covariance = eigenvectors @ np.diag(np.maximum(eigenvalues, 1.0e-7)) @ eigenvectors.T
    inverse = np.linalg.inv(covariance)
    dx = x_grid - component["mean_x"]
    dw = w_grid - component["mean_w"]
    exponent = inverse[0, 0] * dx**2 + 2.0 * inverse[0, 1] * dx * dw + inverse[1, 1] * dw**2
    return np.exp(-0.5 * exponent) / (2.0 * pi * sqrt(max(float(np.linalg.det(covariance)), 1.0e-15)))


def _ellipse(component: dict, radius: float = 2.0) -> np.ndarray:
    covariance = np.array(
        [
            [component["var_x"], component["cov_wx"]],
            [component["cov_wx"], component["var_w"]],
        ]
    )
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    transform = eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))
    angles = np.linspace(0.0, 2.0 * pi, 160)
    return np.array([[component["mean_x"]], [component["mean_w"]]]) + radius * transform @ np.vstack((np.cos(angles), np.sin(angles)))


def make_figure(values: tuple[float, ...], projection: str = "w_rt") -> go.Figure:
    """Plot one projected ADG1 mixture and its two diagnosed Gaussians."""

    title, x_title, x_symbol, x_min, x_max = PROJECTION_DETAILS[projection]
    mixture, _, _ = _diagnosis(values)
    components = projected_components(mixture, projection, chi_mean=CHI_MEAN)
    x_values = np.linspace(x_min, x_max, 161)
    w_values = np.linspace(-5.0, 5.0, 161)
    x_grid, w_grid = np.meshgrid(x_values, w_values)
    component_densities = [_density(x_grid, w_grid, component) for component in components]
    mixture_density = sum(component["weight"] * density for component, density in zip(components, component_densities))

    figure = go.Figure()
    figure.add_trace(
        go.Contour(
            x=x_values,
            y=w_values,
            z=mixture_density,
            colorscale=[[0.0, "#172554"], [0.45, "#1d4ed8"], [1.0, "#bfdbfe"]],
            contours={"coloring": "heatmap", "showlabels": False},
            line={"width": 0},
            showscale=False,
            hovertemplate=f"{x_symbol}=%{{x:.3f}}<br>w=%{{y:.3f}}<br>density=%{{z:.3g}}<extra></extra>",
            name=f"{mixture.label} density",
        )
    )
    for index, (component, density, color) in enumerate(zip(components, component_densities, COMPONENT_COLORS), start=1):
        figure.add_trace(
            go.Contour(
                x=x_values,
                y=w_values,
                z=component["weight"] * density,
                contours={"coloring": "lines", "start": float(np.max(density)) * component["weight"] * 0.04, "end": float(np.max(density)) * component["weight"] * 0.75, "size": float(np.max(density)) * component["weight"] * 0.18},
                line={"color": color, "width": 2.2},
                showscale=False,
                hovertemplate=f"Component {index}<br>weighted density=%{{z:.3g}}<extra></extra>",
                name=f"Component {index} density",
                showlegend=False,
            )
        )
        ellipse = _ellipse(component)
        ellipse_x = np.where(ellipse[0] >= 0.0, ellipse[0], np.nan) if projection == "w_rc" else ellipse[0]
        figure.add_trace(
            go.Scatter(
                x=ellipse_x, y=ellipse[1], mode="lines",
                line={"color": color, "width": 2.8, "dash": "dash"},
                name=f"Component {index} · 2σ",
                legendgroup=f"component-{index}",
                showlegend=False,
                hovertemplate=f"Component {index}<br>weight={component['weight']:.3f}<extra></extra>",
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[max(component["mean_x"], 0.0) if projection == "w_rc" else component["mean_x"]], y=[component["mean_w"]], mode="markers+text",
                text=[f"G{index}"], textposition="top center",
                marker={"color": color, "size": 11, "symbol": "x", "line": {"width": 1, "color": "#0f172a"}},
                name=f"Component {index} center",
                legendgroup=f"component-{index}",
                showlegend=False,
                hovertemplate=f"Component {index} center<br>{x_symbol}=%{{x:.3f}}<br>w=%{{y:.3f}}<extra></extra>",
            )
        )
    figure.add_hline(y=0.0, line={"color": "#94a3b8", "width": 1, "dash": "dot"})
    if projection in {"w_chi", "w_rc"}:
        figure.add_vline(x=0.0, line={"color": "#f59e0b", "width": 1.5, "dash": "dot"})
    else:
        figure.add_vline(x=0.0, line={"color": "#94a3b8", "width": 1, "dash": "dot"})
    figure.update_xaxes(title_text=x_title, range=[x_min, x_max], zeroline=False)
    figure.update_yaxes(title_text="Normalized vertical-velocity anomaly, w", range=[-5.0, 5.0], zeroline=False)
    figure.update_layout(
        template="plotly_dark",
        title={"text": title, "x": 0.5, "font": {"size": 17}},
        height=420,
        margin={"l": 56, "r": 14, "t": 42, "b": 52},
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        showlegend=False,
        uirevision=f"notes-adg1-{projection}",
    )
    return figure


def _component_correlation(covariance: np.ndarray, first: int, second: int) -> float:
    return float(
        covariance[first, second]
        / sqrt(max(covariance[first, first] * covariance[second, second], 1.0e-12))
    )


def _component_summary(values: tuple[float, ...]):
    mixture, _, _ = _diagnosis(values)
    return html.Div(
        [
            html.Div(
                [
                    html.Span(f"Gaussian {index}", className="adg1-component-title"),
                    html.Strong(f"weight {weight:.1%}"),
                    html.Span(
                        f"center  w={mean[0]:+.2f}, rₜ={mean[1]:+.2f}, θₗ={mean[2]:+.2f}"
                    ),
                    html.Span(
                        f"widths  σw={sqrt(max(covariance[0, 0], 0.0)):.2f}, σrₜ={sqrt(max(covariance[1, 1], 0.0)):.2f}, σθₗ={sqrt(max(covariance[2, 2], 0.0)):.2f}",
                        className="adg1-component-detail",
                    ),
                    html.Span(
                        f"internal ρ(w,rₜ)={_component_correlation(covariance, 0, 1):+.2f} · "
                        f"ρ(w,θₗ)={_component_correlation(covariance, 0, 2):+.2f} · "
                        f"ρ(rₜ,θₗ)={_component_correlation(covariance, 1, 2):+.2f}",
                        className="adg1-component-detail",
                    ),
                ],
                className="adg1-component-card",
                style={"--adg1-component-color": color},
            )
            for index, (weight, mean, covariance, color) in enumerate(
                zip(mixture.weights, mixture.means, mixture.covariances, COMPONENT_COLORS),
                start=1,
            )
        ],
        className="adg1-component-summary",
    )


def _status(values: tuple[float, ...]):
    mixture, _, adjusted = _diagnosis(values)
    messages = []
    if adjusted:
        messages.append("The requested correlations were moved to the nearest physical correlation matrix.")
    if mixture.adjustment:
        messages.append(mixture.adjustment)
    messages.append(mixture.note)
    messages.append("ADG1 fixes every within-component w–scalar tilt to zero; the grid covariance comes from the separation between their centers.")
    messages.append("The χ and rᶜ panels use χ̄ = −0.25 only to display the cloud boundary.")
    return " ".join(messages)


def _control_group(title: str, controls: tuple, **properties):
    return html.Div(
        [
            html.H3(title, className="adg1-control-group-title"),
            html.Div([_control(*control) for control in controls], className="adg1-controls-grid"),
        ],
        className="adg1-control-group",
        **properties,
    )


def build_layout():
    return html.Div(
        [
            report_header(
                "ADG1 two-Gaussian explorer",
                "Explore ADG1 across w–rₜ, w–χ, and cloudy-only w–rᶜ. The dashed loops are 2σ component footprints; × marks their centers.",
                eyebrow="Interactive PDF geometry",
                badges=("ADG1", "three projections"),
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H2("Inputs ADG1 receives"),
                                    html.Button(
                                        "Reset",
                                        id=RESET_ID,
                                        n_clicks=0,
                                        className="adg1-reset-button",
                                    ),
                                ],
                                className="adg1-controls-heading",
                            ),
                            html.P(
                                "All three variances and correlations enter the common covariance matrix. w skewness controls the two weights and their vertical spacing. γ and β are ADG1's advanced width controls.",
                                className="adg1-controls-intro",
                            ),
                            _control_group("Shared target covariance", COMMON_CONTROLS),
                            _control_group("ADG1 diagnosis controls", ADG1_CONTROLS, id=ADG1_CONTROL_GROUP_ID),
                        ],
                        className="adg1-controls-panel",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    dcc.Graph(id=FIGURE_IDS[projection], config={"displaylogo": False, "responsive": True}, className="adg1-gaussian-graph")
                                    for projection in PROJECTIONS
                                ],
                                className="adg1-figure-row",
                            ),
                            html.P("Gold is Gaussian 1; cyan is Gaussian 2. Colored contours show each weighted component density. In the χ and rᶜ panels, the orange dotted line is the cloud boundary; rᶜ shows only its cloudy side.", className="adg1-figure-key"),
                            html.Div(id=SUMMARY_ID),
                            html.P(id=STATUS_ID, className="adg1-status"),
                        ],
                        className="adg1-visual-panel",
                    ),
                ],
                className="adg1-workspace",
            ),
            report_section(
                "What this view is showing",
                html.P("This is a compact, normalized geometry workspace. γ drives within-Gaussian w width after the resolved-correlation limiter, and β allocates remaining scalar width. CLUBB's Fortran PDF diagnosis remains authoritative for physical cases."),
                class_name="adg1-method-note",
            ),
        ],
        className="notes-report adg1-report",
    )


def register_callbacks(app):
    input_ids = tuple(_input_id(control[0]) for control in CONTROLS)
    value_ids = tuple(_value_id(control[0]) for control in CONTROLS)

    @app.callback(
        *(Output(FIGURE_IDS[projection], "figure") for projection in PROJECTIONS),
        Output(SUMMARY_ID, "children"),
        Output(STATUS_ID, "children"),
        *(Output(identifier, "children") for identifier in value_ids),
        *(Input(identifier, "value") for identifier in input_ids),
    )
    def _update(*values):
        value_tuple = tuple(float(value) for value in values)
        settings = _settings(value_tuple)
        return (
            *(make_figure(value_tuple, projection) for projection in PROJECTIONS),
            _component_summary(value_tuple),
            _status(value_tuple),
            *(_format_value(name, settings[name]) for name, *_ in CONTROLS),
        )

    @app.callback(
        *(Output(identifier, "value") for identifier in input_ids),
        Input(RESET_ID, "n_clicks"),
        prevent_initial_call=True,
    )
    def _reset(_clicks):
        return tuple(DEFAULTS[name] for name, *_ in CONTROLS)

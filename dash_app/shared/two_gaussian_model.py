"""Small, standardized double-Gaussian teaching model.

The routines in this module preserve the structural constraints that matter for
comparing ADG1, the new PDF, and a flexible binormal.  They intentionally work
in normalized coordinates; the full CLUBB Fortran remains the authoritative
implementation of each closure.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, exp, pi, sqrt

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dash_app.shared.bivariate_heatmap import (
    BIVARIATE_LEVELS,
    bivariate_codes,
    bivariate_reference_colorscale as _bivariate_colorscale,
)


VARIABLES = ("w", "rₜ", "θₗ")
MODEL_COLORS = {"adg1": "#f97316", "new": "#3b82f6", "free": "#a855f7"}

@dataclass(frozen=True)
class Mixture:
    key: str
    label: str
    weights: np.ndarray
    means: np.ndarray
    covariances: np.ndarray
    note: str
    adjustment: str = ""


def _is_psd(matrix: np.ndarray, tolerance: float = 1.0e-9) -> bool:
    return bool(np.min(np.linalg.eigvalsh(0.5 * (matrix + matrix.T))) >= -tolerance)


def nearest_correlation(values: np.ndarray) -> tuple[np.ndarray, bool]:
    """Return a positive-semidefinite correlation matrix with unit diagonal."""
    matrix = np.asarray(values, dtype=float)
    matrix = 0.5 * (matrix + matrix.T)
    np.fill_diagonal(matrix, 1.0)
    changed = not _is_psd(matrix)
    for _ in range(4):
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        eigenvalues = np.maximum(eigenvalues, 1.0e-7)
        matrix = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T
        scales = np.sqrt(np.maximum(np.diag(matrix), 1.0e-12))
        matrix = matrix / np.outer(scales, scales)
        np.fill_diagonal(matrix, 1.0)
    np.fill_diagonal(matrix, 1.0)
    return matrix, changed


def grid_covariance(
    variances: tuple[float, float, float],
    correlations: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray, bool]:
    corr = np.array(
        [
            [1.0, correlations[0], correlations[1]],
            [correlations[0], 1.0, correlations[2]],
            [correlations[1], correlations[2], 1.0],
        ],
        dtype=float,
    )
    corr, adjusted = nearest_correlation(corr)
    standard_deviations = np.sqrt(np.maximum(np.asarray(variances, dtype=float), 0.05))
    return corr * np.outer(standard_deviations, standard_deviations), corr, adjusted


def _adg_weight(skew_w: float, within_fraction: float) -> float:
    if abs(skew_w) <= 1.0e-7:
        return 0.5
    denominator = sqrt(max(4.0 * (1.0 - within_fraction) ** 3 + skew_w**2, 1.0e-12))
    return float(np.clip(0.5 * (1.0 - skew_w / denominator), 0.04, 0.96))


def _component_means(weight: float, displacement: np.ndarray) -> np.ndarray:
    return np.vstack(((1.0 - weight) * displacement, -weight * displacement))


def _between_covariance(weight: float, displacement: np.ndarray) -> np.ndarray:
    return weight * (1.0 - weight) * np.outer(displacement, displacement)


def _scale_displacement_to_psd(
    covariance: np.ndarray, weight: float, displacement: np.ndarray
) -> tuple[np.ndarray, float]:
    scale = 1.0
    for _ in range(80):
        candidate = displacement * scale
        residual = covariance - _between_covariance(weight, candidate)
        # Leave a small interior margin.  A merely semidefinite residual can
        # become slightly negative when the component contrast is applied.
        eigenvalue_floor = 1.0e-6 * max(float(np.max(np.diag(residual))), 1.0)
        if float(np.min(np.linalg.eigvalsh(residual))) >= eigenvalue_floor:
            return candidate, scale
        scale *= 0.96
    return np.zeros_like(displacement), 0.0


def _shared_correlation_covariances(
    residual: np.ndarray,
    weight: float,
    diagonal_contrast: np.ndarray,
) -> tuple[np.ndarray, str]:
    diagonal = np.maximum(np.diag(residual), 1.0e-8)
    contrast = np.clip(np.asarray(diagonal_contrast, dtype=float), -0.9, 0.9)
    variance_1 = diagonal * (1.0 + (1.0 - weight) * contrast)
    variance_2 = diagonal * (1.0 - weight * contrast)
    corr = np.eye(3)
    for row in range(3):
        for column in range(row + 1, 3):
            denominator = (
                weight * sqrt(variance_1[row] * variance_1[column])
                + (1.0 - weight) * sqrt(variance_2[row] * variance_2[column])
            )
            corr[row, column] = corr[column, row] = residual[row, column] / max(denominator, 1.0e-12)
    corr_requested = corr.copy()
    corr, adjusted = nearest_correlation(corr)
    result = np.empty((2, 3, 3), dtype=float)
    for component, variances in enumerate((variance_1, variance_2)):
        result[component] = corr * np.outer(np.sqrt(variances), np.sqrt(variances))
    note = ""
    if adjusted or np.max(np.abs(corr - corr_requested)) > 1.0e-5:
        note = "Shared correlations were moved to the nearest realizable matrix."
    return result, note


def diagnose_adg1(covariance: np.ndarray, skew_w: float) -> Mixture:
    within_w = 0.40
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

    component_covariances = np.zeros((2, 3, 3), dtype=float)
    component_covariances[:, 0, 0] = residual[0, 0]
    # ADG1's beta=1.5 form gives this particularly transparent allocation.
    for variable in (1, 2):
        component_covariances[0, variable, variable] = residual[variable, variable] / (2.0 * weight)
        component_covariances[1, variable, variable] = residual[variable, variable] / (2.0 * (1.0 - weight))
    denominator = (
        weight
        * sqrt(component_covariances[0, 1, 1] * component_covariances[0, 2, 2])
        + (1.0 - weight)
        * sqrt(component_covariances[1, 1, 1] * component_covariances[1, 2, 2])
    )
    scalar_rho = float(np.clip(residual[1, 2] / max(denominator, 1.0e-12), -0.98, 0.98))
    for component in range(2):
        scalar_covariance = scalar_rho * sqrt(
            component_covariances[component, 1, 1] * component_covariances[component, 2, 2]
        )
        component_covariances[component, 1, 2] = scalar_covariance
        component_covariances[component, 2, 1] = scalar_covariance

    adjustment = ""
    if scale < 0.999:
        adjustment = "Component separation was reduced to keep the ADG1 matrices realizable."
    return Mixture(
        "adg1",
        "ADG1",
        np.array([weight, 1.0 - weight]),
        _component_means(weight, displacement),
        component_covariances,
        "Grid w–scalar covariance must come from separation of the component centers; each component's w–scalar tilt is zero.",
        adjustment,
    )


def diagnose_new_pdf(
    covariance: np.ndarray,
    skewness: tuple[float, float, float],
) -> Mixture:
    within_w = 0.32
    weight = _adg_weight(skewness[0], within_w)
    displacement_w = sqrt(
        max((1.0 - within_w) * covariance[0, 0], 1.0e-10)
        / (weight * (1.0 - weight))
    )
    displacement = np.zeros(3)
    displacement[0] = displacement_w
    for variable in (1, 2):
        sign = 1.0 if covariance[0, variable] >= 0.0 else -1.0
        spread_fraction = np.clip(0.12 + 0.16 * abs(skewness[variable]), 0.08, 0.62)
        displacement[variable] = sign * sqrt(
            spread_fraction * covariance[variable, variable] / (weight * (1.0 - weight))
        )
    displacement, scale = _scale_displacement_to_psd(covariance, weight, displacement)
    residual = covariance - _between_covariance(weight, displacement)
    diagonal_contrast = np.tanh(np.asarray(skewness) * 0.35)
    # The weight and w separation above were solved from Skw assuming equal
    # within-component w widths.  Retain that equality so the reconstructed
    # mixture preserves the supplied third moment exactly.
    diagonal_contrast[0] = 0.0
    component_covariances, correlation_note = _shared_correlation_covariances(
        residual, weight, diagonal_contrast
    )
    notes = [correlation_note] if correlation_note else []
    if scale < 0.999:
        notes.append("Mean separation was reduced to remain realizable.")
    return Mixture(
        "new",
        "New PDF",
        np.array([weight, 1.0 - weight]),
        _component_means(weight, displacement),
        component_covariances,
        "Skewness sets component geometry; the unexplained grid covariance becomes the same correlation coefficient in both components.",
        " ".join(notes),
    )


def diagnose_flexible(
    covariance: np.ndarray,
    weight: float,
    separation: tuple[float, float, float],
    variance_contrast: float,
    tilt_contrast: float,
) -> Mixture:
    weight = float(np.clip(weight, 0.04, 0.96))
    standard_deviations = np.sqrt(np.maximum(np.diag(covariance), 1.0e-8))
    displacement = np.asarray(separation, dtype=float) * standard_deviations / sqrt(
        weight * (1.0 - weight)
    )
    displacement, separation_scale = _scale_displacement_to_psd(
        covariance, weight, displacement
    )
    residual = covariance - _between_covariance(weight, displacement)
    contrast_matrix = np.zeros((3, 3), dtype=float)
    contrast_matrix[np.diag_indices(3)] = variance_contrast * np.diag(residual)
    for row in range(3):
        for column in range(row + 1, 3):
            value = tilt_contrast * 0.65 * sqrt(
                max(residual[row, row] * residual[column, column], 0.0)
            )
            # Alternate the scalar-scalar direction so the control visibly
            # reallocates rather than uniformly inflating all covariance.
            if (row, column) == (1, 2):
                value *= -0.5
            contrast_matrix[row, column] = contrast_matrix[column, row] = value
    contrast_scale = 1.0
    while contrast_scale > 1.0e-5:
        component_1 = residual + (1.0 - weight) * contrast_scale * contrast_matrix
        component_2 = residual - weight * contrast_scale * contrast_matrix
        if _is_psd(component_1) and _is_psd(component_2):
            break
        contrast_scale *= 0.92
    if not (_is_psd(component_1) and _is_psd(component_2)):
        contrast_scale = 0.0
        component_1 = residual.copy()
        component_2 = residual.copy()
    component_covariances = np.stack((component_1, component_2))
    messages = []
    if separation_scale < 0.999:
        messages.append("Separation was reduced to preserve a physical residual covariance.")
    if contrast_scale < 0.999:
        messages.append("Component contrast was reduced to keep both matrices physical.")
    return Mixture(
        "free",
        "Flexible 2G",
        np.array([weight, 1.0 - weight]),
        _component_means(weight, displacement),
        component_covariances,
        "Component geometry is adjustable, but the mixture mean and covariance are still held equal to the common CLUBB inputs.",
        " ".join(messages),
    )


def rare_updraft_reference() -> Mixture:
    """Return the LES-diagnostic-fitted flexible 2G used as the background."""
    # Local import avoids a module cycle: oracle_presets constructs Mixture.
    from .oracle_presets import rare_updraft_family_oracles

    reference = rare_updraft_family_oracles()[-1]
    return Mixture(
        "reference",
        "LES-diagnostic-fitted flexible 2G",
        reference.weights,
        reference.means,
        reference.covariances,
        "A constrained inference from LES moments, not the raw LES joint PDF.",
        reference.adjustment,
    )


def mixture_moments(mixture: Mixture) -> tuple[np.ndarray, np.ndarray, float]:
    mean = np.sum(mixture.weights[:, None] * mixture.means, axis=0)
    covariance = np.zeros((3, 3))
    wp3 = 0.0
    for weight, component_mean, component_covariance in zip(
        mixture.weights, mixture.means, mixture.covariances
    ):
        offset = component_mean - mean
        covariance += weight * (component_covariance + np.outer(offset, offset))
        wp3 += weight * (offset[0] ** 3 + 3.0 * offset[0] * component_covariance[0, 0])
    return mean, covariance, float(wp3)


def projected_components(
    mixture: Mixture,
    projection: str,
    chi_mean: float,
    chi_coefficients: tuple[float, float] = (1.0, 1.0),
):
    rt_coefficient, thl_coefficient = chi_coefficients
    projected = []
    for weight, mean, covariance in zip(mixture.weights, mixture.means, mixture.covariances):
        if projection == "w_rt":
            x_mean = mean[1]
            x_variance = covariance[1, 1]
            wx_covariance = covariance[0, 1]
        else:
            x_mean = chi_mean + rt_coefficient * mean[1] - thl_coefficient * mean[2]
            x_variance = (
                rt_coefficient**2 * covariance[1, 1]
                + thl_coefficient**2 * covariance[2, 2]
                - 2.0 * rt_coefficient * thl_coefficient * covariance[1, 2]
            )
            wx_covariance = (
                rt_coefficient * covariance[0, 1]
                - thl_coefficient * covariance[0, 2]
            )
        projected.append(
            {
                "weight": float(weight),
                "mean_w": float(mean[0]),
                "mean_x": float(x_mean),
                "var_w": max(float(covariance[0, 0]), 1.0e-8),
                "var_x": max(float(x_variance), 1.0e-8),
                "cov_wx": float(wx_covariance),
            }
        )
    return projected


def _density(x_grid, w_grid, component):
    covariance = np.array(
        [
            [component["var_x"], component["cov_wx"]],
            [component["cov_wx"], component["var_w"]],
        ]
    )
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    covariance = eigenvectors @ np.diag(np.maximum(eigenvalues, 1.0e-7)) @ eigenvectors.T
    inverse = np.linalg.inv(covariance)
    delta_x = x_grid - component["mean_x"]
    delta_w = w_grid - component["mean_w"]
    exponent = (
        inverse[0, 0] * delta_x**2
        + 2.0 * inverse[0, 1] * delta_x * delta_w
        + inverse[1, 1] * delta_w**2
    )
    return np.exp(-0.5 * exponent) / (2.0 * pi * sqrt(max(np.linalg.det(covariance), 1.0e-15)))


def _positive_normal_mean(mean, variance):
    variance = max(float(variance), 1.0e-12)
    sigma = sqrt(variance)
    normalized = np.asarray(mean, dtype=float) / sigma
    probability = 0.5 * (
        1.0 + np.vectorize(erf, otypes=[float])(normalized / sqrt(2.0))
    )
    normal_density = np.exp(-0.5 * normalized**2) / sqrt(2.0 * pi)
    return np.asarray(mean, dtype=float) * probability + sigma * normal_density


def reference_transport_importance(
    mixture: Mixture,
    projection: str,
    chi_mean: float,
    x_grid: np.ndarray,
    w_grid: np.ndarray,
    chi_coefficients: tuple[float, float] = (1.0, 1.0),
) -> np.ndarray:
    """Density weighted by positive upward cloud-water transport.

    For w-rt, cloud water is integrated over the conditional thl distribution.
    For w-chi and w-rc, the cloud coordinate is explicit on the horizontal axis.
    The field is an explanatory importance map, not a probability density.
    """
    importance = np.zeros_like(x_grid, dtype=float)
    upward_velocity = np.maximum(w_grid, 0.0)
    if projection in {"w_chi", "w_rc"}:
        components = projected_components(
            mixture, projection, chi_mean, chi_coefficients
        )
        cloud_water = np.maximum(x_grid, 0.0)
        for component in components:
            importance += (
                component["weight"]
                * _density(x_grid, w_grid, component)
                * upward_velocity
                * cloud_water
            )
        return importance

    for weight, mean, covariance in zip(
        mixture.weights, mixture.means, mixture.covariances
    ):
        component = {
            "mean_w": float(mean[0]),
            "mean_x": float(mean[1]),
            "var_w": max(float(covariance[0, 0]), 1.0e-8),
            "var_x": max(float(covariance[1, 1]), 1.0e-8),
            "cov_wx": float(covariance[0, 1]),
        }
        joint_density = _density(x_grid, w_grid, component)
        observed_covariance = covariance[np.ix_((0, 1), (0, 1))]
        hidden_covariance = covariance[2, (0, 1)]
        regression = hidden_covariance @ np.linalg.inv(observed_covariance)
        conditional_thl_mean = (
            mean[2]
            + regression[0] * (w_grid - mean[0])
            + regression[1] * (x_grid - mean[1])
        )
        conditional_thl_variance = max(
            float(covariance[2, 2] - regression @ hidden_covariance), 1.0e-12
        )
        rt_coefficient, thl_coefficient = chi_coefficients
        conditional_chi_mean = (
            chi_mean
            + rt_coefficient * x_grid
            - thl_coefficient * conditional_thl_mean
        )
        conditional_cloud_water = _positive_normal_mean(
            conditional_chi_mean,
            thl_coefficient**2 * conditional_thl_variance,
        )
        importance += (
            float(weight)
            * joint_density
            * upward_velocity
            * conditional_cloud_water
        )
    return importance


def _ellipse(component, radius=2.0):
    covariance = np.array(
        [
            [component["var_x"], component["cov_wx"]],
            [component["cov_wx"], component["var_w"]],
        ]
    )
    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    transform = eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))
    angles = np.linspace(0.0, 2.0 * pi, 160)
    points = np.array([[component["mean_x"]], [component["mean_w"]]]) + radius * transform @ np.vstack(
        (np.cos(angles), np.sin(angles))
    )
    return points


def cloud_diagnostics(
    mixture: Mixture,
    chi_mean: float,
    chi_coefficients: tuple[float, float] = (1.0, 1.0),
) -> tuple[float, float]:
    cloud_fraction = 0.0
    wprcp = 0.0
    overall_w = float(np.sum(mixture.weights * mixture.means[:, 0]))
    for weight, component in zip(
        mixture.weights,
        projected_components(mixture, "w_chi", chi_mean, chi_coefficients),
    ):
        sigma_chi = sqrt(component["var_x"])
        normalized_mean = component["mean_x"] / sigma_chi
        probability = 0.5 * (1.0 + erf(normalized_mean / sqrt(2.0)))
        normal_density = exp(-0.5 * normalized_mean**2) / sqrt(2.0 * pi)
        mean_positive_chi = component["mean_x"] * probability + sigma_chi * normal_density
        cloud_fraction += weight * probability
        wprcp += weight * (
            (component["mean_w"] - overall_w) * mean_positive_chi
            + component["cov_wx"] * probability
        )
    return float(cloud_fraction), float(wprcp)


def covariance_matrix_rows(mixture: Mixture):
    rows = []
    for component in range(2):
        covariance = mixture.covariances[component]
        std = np.sqrt(np.maximum(np.diag(covariance), 1.0e-12))
        correlation = covariance / np.outer(std, std)
        rows.append(
            {
                "component": component + 1,
                "weight": mixture.weights[component],
                "means": mixture.means[component],
                "covariance": covariance,
                "correlation": correlation,
                "min_eigenvalue": float(np.min(np.linalg.eigvalsh(covariance))),
            }
        )
    return rows


def make_comparison_figure(
    mixtures: list[Mixture],
    projection: str,
    chi_mean: float,
    reference: Mixture | None = None,
    reference_chi_mean: float | None = None,
    chi_coefficients: tuple[float, float] = (1.0, 1.0),
    reference_chi_coefficients: tuple[float, float] | None = None,
) -> go.Figure:
    labels = [mixture.label for mixture in mixtures]
    figure = make_subplots(rows=1, cols=3, subplot_titles=labels, shared_yaxes=True)
    all_components = [
        projected_components(mixture, projection, chi_mean, chi_coefficients)
        for mixture in mixtures
    ]
    reference_components = None
    if reference is not None:
        reference_components = projected_components(
            reference,
            projection,
            chi_mean if reference_chi_mean is None else reference_chi_mean,
            (
                chi_coefficients
                if reference_chi_coefficients is None
                else reference_chi_coefficients
            ),
        )
    # Interactive ranges remain stationary while sliders move.  The fixed
    # snapshot is allowed one wider frame so its broad rare component is not
    # cropped; it then remains stable while projections are inspected.
    if reference_components is None:
        x_min, x_max = (0.0, 5.0) if projection == "w_rc" else (-5.0, 5.0)
        w_min, w_max = -5.0, 5.0
    else:
        visible_components = [
            component
            for components in (*all_components, reference_components)
            for component in components
        ]
        w_low = min(
            component["mean_w"] - 2.2 * sqrt(component["var_w"])
            for component in visible_components
        )
        w_high = max(
            component["mean_w"] + 2.2 * sqrt(component["var_w"])
            for component in visible_components
        )
        w_min, w_max = min(-5.0, w_low - 0.4), max(5.0, w_high + 0.4)
        x_low = min(
            component["mean_x"] - 2.2 * sqrt(component["var_x"])
            for component in visible_components
        )
        x_high = max(
            component["mean_x"] + 2.2 * sqrt(component["var_x"])
            for component in visible_components
        )
        if projection == "w_rc":
            x_min, x_max = 0.0, max(5.0, x_high + 0.4)
        else:
            x_min, x_max = min(-5.0, x_low - 0.4), max(5.0, x_high + 0.4)
    x_values = np.linspace(x_min, x_max, 145)
    w_values = np.linspace(w_min, w_max, 145)
    x_grid, w_grid = np.meshgrid(x_values, w_values)

    for column, (mixture, components) in enumerate(zip(mixtures, all_components), start=1):
        density = np.zeros_like(x_grid)
        for component in components:
            density += component["weight"] * _density(x_grid, w_grid, component)
        if reference_components is not None:
            reference_density = np.zeros_like(x_grid)
            for component in reference_components:
                reference_density += component["weight"] * _density(x_grid, w_grid, component)
            reference_importance = reference_transport_importance(
                reference,
                projection,
                chi_mean if reference_chi_mean is None else reference_chi_mean,
                x_grid,
                w_grid,
                (
                    chi_coefficients
                    if reference_chi_coefficients is None
                    else reference_chi_coefficients
                ),
            )
            bivariate_code, normalized_probability, normalized_importance = bivariate_codes(
                reference_density, reference_importance
            )
            figure.add_trace(
                go.Heatmap(
                    x=x_values,
                    y=w_values,
                    z=bivariate_code,
                    zmin=0,
                    zmax=BIVARIATE_LEVELS * BIVARIATE_LEVELS - 1,
                    colorscale=_bivariate_colorscale(),
                    customdata=np.stack(
                        (
                            normalized_probability,
                            normalized_importance,
                            reference_density,
                            reference_importance,
                        ),
                        axis=-1,
                    ),
                    showscale=False,
                    zsmooth=False,
                    showlegend=column == 1,
                    name="Target: probability + transport",
                    hovertemplate=(
                        "x=%{x:.3f}<br>w=%{y:.3f}"
                        "<br>normalized probability=%{customdata[0]:.2f}"
                        "<br>normalized transport importance=%{customdata[1]:.2f}"
                        "<br>density=%{customdata[2]:.3g}"
                        "<br>transport importance=%{customdata[3]:.3g}"
                        f"<extra>{reference.label}</extra>"
                    ),
                ),
                row=1,
                col=column,
            )
        figure.add_trace(
            go.Contour(
                x=x_values,
                y=w_values,
                z=density,
                colorscale=[[0.0, MODEL_COLORS[mixture.key]], [1.0, MODEL_COLORS[mixture.key]]],
                contours={"coloring": "lines"} if reference is not None else {"coloring": "heatmap"},
                line={"width": 2.0 if reference is not None else 0, "color": MODEL_COLORS[mixture.key]},
                showscale=False,
                showlegend=False,
                hovertemplate="x=%{x:.3f}<br>w=%{y:.3f}<br>density=%{z:.3g}<extra>analytic mixture density</extra>",
            ),
            row=1,
            col=column,
        )
        for index, component in enumerate(components):
            ellipse = _ellipse(component)
            if projection == "w_rc":
                mask = ellipse[0] >= 0.0
                ellipse_x = np.where(mask, ellipse[0], np.nan)
            else:
                ellipse_x = ellipse[0]
            figure.add_trace(
                go.Scatter(
                    x=ellipse_x,
                    y=ellipse[1],
                    mode="lines",
                    line={"color": "#fbbf24" if index == 0 else "#22d3ee", "width": 2.5, "dash": "dash"},
                    name=f"Component {index + 1}",
                    legendgroup=f"component-{index}",
                    showlegend=column == 1,
                    hovertemplate=f"Component {index + 1}<br>weight={component['weight']:.3f}<extra></extra>",
                ),
                row=1,
                col=column,
            )
            figure.add_trace(
                go.Scatter(
                    x=[max(component["mean_x"], 0.0) if projection == "w_rc" else component["mean_x"]],
                    y=[component["mean_w"]],
                    mode="markers",
                    marker={"color": "#fbbf24" if index == 0 else "#22d3ee", "size": 9, "symbol": "x"},
                    showlegend=False,
                    hovertemplate=f"Component {index + 1} center<extra></extra>",
                ),
                row=1,
                col=column,
            )
        if projection in {"w_chi", "w_rc"}:
            figure.add_vline(x=0.0, line={"color": "#f59e0b", "dash": "dot", "width": 2}, row=1, col=column)

    x_title = {"w_rt": "Normalized total water, rₜ", "w_chi": "Normalized extended cloud water, χ", "w_rc": "Normalized cloud water, rᶜ"}[projection]
    figure.update_xaxes(title_text=x_title, range=[x_min, x_max])
    figure.update_yaxes(title_text="Normalized vertical velocity, w", range=[w_min, w_max])
    figure.update_layout(
        height=445,
        margin={"l": 52, "r": 18, "t": 48, "b": 58},
        template="plotly_dark",
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        legend={"orientation": "h", "y": -0.16, "x": 0.02},
        uirevision=projection,
        autosize=True,
    )
    return figure

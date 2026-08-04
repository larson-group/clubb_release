from collections import OrderedDict
from functools import lru_cache
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from dash import ALL, Input, MATCH, Output, State, callback_context, dcc, html, no_update
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter
from scipy.special import ndtr

from dash_app.shared.bivariate_heatmap import (
    BIVARIATE_LEVELS,
    robust_field_upper,
    signed_bivariate_codes,
    signed_bivariate_reference_colorscale,
)
from dash_app.shared.components import styled_dropdown
from utilities.sam_3d_reference import (
    DEFAULT_BOMEX_SAM_RUN,
    DEFAULT_SAM_RUN,
    inventory_run,
    load_snapshot,
)

from .. import benchmark_overlay
from . import shared
from .base_plot import BasePlotType


COMMON_PDF_FIELDS = (
    "mixt_frac",
    "w_1",
    "w_2",
    "varnce_w_1",
    "varnce_w_2",
)

PROJECTIONS = {
    "w_chi": {
        "label": "w – chi (cloud diagnostic)",
        "x_title": "Extended cloud water, χ [g/kg]",
        "y_title": "Vertical velocity, w [m/s]",
        "axes": ("chi", "w"),
        "x_mean": ("chi_1", "chi_2"),
        "x_spread": ("stdev_chi_1", "stdev_chi_2"),
        "spread_is_stdev": True,
        "correlation": ("corr_w_chi_1", "corr_w_chi_2"),
        "cloud_threshold": True,
    },
    "w_rt": {
        "label": "w – rt",
        "x_title": "Total water, rₜ [g/kg]",
        "y_title": "Vertical velocity, w [m/s]",
        "axes": ("rt", "w"),
        "x_mean": ("rt_1", "rt_2"),
        "x_spread": ("varnce_rt_1", "varnce_rt_2"),
        "spread_is_stdev": False,
        "correlation": ("corr_w_rt_1", "corr_w_rt_2"),
        "cloud_threshold": False,
    },
    "w_rc": {
        "label": "w – rc (cloud-water covariance)",
        "x_title": "Cloud water, r꜀ [g/kg]",
        "y_title": "Vertical velocity, w [m/s]",
        "axes": ("rc", "w"),
        "x_mean": ("chi_1", "chi_2"),
        "x_spread": ("stdev_chi_1", "stdev_chi_2"),
        "spread_is_stdev": True,
        "correlation": ("corr_w_chi_1", "corr_w_chi_2"),
        "cloud_threshold": True,
    },
    "w_thl": {
        "label": "w – θₗ",
        "x_title": "Liquid-water potential temperature, θₗ [K]",
        "y_title": "Vertical velocity, w [m/s]",
        "axes": ("thl", "w"),
        "x_mean": ("thl_1", "thl_2"),
        "x_spread": ("varnce_thl_1", "varnce_thl_2"),
        "spread_is_stdev": False,
        "correlation": ("corr_w_thl_1", "corr_w_thl_2"),
        "cloud_threshold": False,
        "integrate_w_transport": True,
    },
    "rt_thl": {
        "label": "rₜ – θₗ",
        "x_title": "Total water, rₜ [g/kg]",
        "y_title": "Liquid-water potential temperature, θₗ [K]",
        "axes": ("rt", "thl"),
        "x_mean": ("rt_1", "rt_2"),
        "x_spread": ("varnce_rt_1", "varnce_rt_2"),
        "spread_is_stdev": False,
        "correlation": ("corr_rt_thl_1", "corr_rt_thl_2"),
        "cloud_threshold": False,
        "integrate_w_transport": True,
    },
    "chi_thl": {
        "label": "χ – θₗ",
        "x_title": "Extended cloud water, χ [g/kg]",
        "y_title": "Liquid-water potential temperature, θₗ [K]",
        "axes": ("chi", "thl"),
        "x_mean": ("chi_1", "chi_2"),
        "x_spread": ("stdev_chi_1", "stdev_chi_2"),
        "spread_is_stdev": True,
        # corr(χ, θₗ) is reconstructed from the trivariate component
        # covariance, using the fields below rather than a stored scalar.
        "correlation": ("corr_rt_thl_1", "corr_rt_thl_2"),
        "cloud_threshold": True,
        "integrate_w_transport": True,
    },
}
PROJECTION_OPTIONS = [
    {"label": config["label"], "value": name}
    for name, config in PROJECTIONS.items()
]
# This is visually equivalent to the former 161x161 dashboard raster while
# reducing each time-sample density calculation by about 61 percent.
PDF_GRID_POINTS = 101
ENCLOSED_MASSES = (0.50, 0.80, 0.995)
ENCLOSED_CONTOUR_OPACITY = {0.50: 1.0, 0.80: 0.72, 0.995: 0.48}
ENCLOSED_CONTOUR_WIDTH = {0.50: 2.7, 0.80: 2.1, 0.995: 1.5}
DEFAULT_RAW_CONTOUR_SMOOTHING_BINS = 1.25
COMPONENT_COLORS = ("#ef4444", "#06b6d4", "#22c55e")
MASS_CONTOUR_COLOR = "#7c3aed"
RAW_MASS_CONTOUR_COLOR = "#facc15"
RAW_ENCLOSED_PROBABILITY_LABEL = "Raw SAM enclosed probability (50/80/99.5%)"
TRANSPORT_VIEWS = {"signed": {"label": "Signed"}}
COLOR_SIGNALS = {
    "probability": {
        "control_label": "Probability only",
        "label": "probability",
        "units": "",
        "probability_only": True,
    },
    "wprcp": {
        "control_label": "w′r꜀′ cloud",
        "label": "w′r꜀′ cloud-water covariance",
        "units": "m/s·g/kg",
    },
    "wprtp": {
        "control_label": "w′rₜ′ total water",
        "label": "w′rₜ′ total-water covariance",
        "units": "m/s·g/kg",
    },
    "wpchi": {
        "control_label": "w′χ′ saturation",
        "label": "w′χ′ saturation-excess covariance",
        "units": "m/s·g/kg",
    },
    "wpthlp": {
        "control_label": "w′θₗ′ thermal",
        "label": "w′θₗ′ liquid-water-potential-temperature covariance",
        "units": "m K/s",
    },
    "wprtp2": {
        "control_label": "w′rₜ′² rtp2",
        "label": "w′rₜ′² total-water third-moment contribution",
        "units": "m/s·(g/kg)²",
    },
}
COLOR_SIGNAL_OPTIONS = [
    {"label": config["control_label"], "value": name}
    for name, config in COLOR_SIGNALS.items()
]
ENCLOSED_PROBABILITY_LABEL = "Enclosed probability (50/80/99.5%)"
CLUBB_ENCLOSED_WEIGHTED_LABEL = "CLUBB enclosed |{}| (50/80/99.5%)"
RAW_3D_RUNS = {
    "arm": Path(DEFAULT_SAM_RUN),
    "bomex": Path(DEFAULT_BOMEX_SAM_RUN),
}
_DENSITY_CACHE_MAX_ENTRIES = 48
_DENSITY_CACHE = OrderedDict()
_RAW_HISTOGRAM_CACHE_MAX_ENTRIES = 48
_RAW_HISTOGRAM_CACHE = OrderedDict()
CLUBB_CLOUD_CONDITIONING_FIELDS = (
    "chi_1",
    "chi_2",
    "stdev_chi_1",
    "stdev_chi_2",
    "corr_w_chi_1",
    "corr_w_chi_2",
    "varnce_thl_1",
    "varnce_thl_2",
    "corr_rt_thl_1",
    "corr_rt_thl_2",
    "crt_1",
    "crt_2",
    "cthl_1",
    "cthl_2",
)

THIRD_COMPONENT_FIELDS = (
    "mixt_frac_3",
    "w_3",
    "varnce_w_3",
)

_TRIVARIATE_COMPONENT_FIELDS = (
    "w",
    "rt",
    "thl",
    "chi",
    "stdev_chi",
    "varnce_w",
    "varnce_rt",
    "varnce_thl",
    "corr_w_rt",
    "corr_w_thl",
    "corr_w_chi",
    "corr_rt_thl",
    "crt",
    "cthl",
)


def _trivariate_component_field_names(component):
    return tuple(f"{field}_{int(component)}" for field in _TRIVARIATE_COMPONENT_FIELDS)


class PdfContourError(ValueError):
    pass


def constrained_pdf_heights(height_range, current_values):
    """Clamp local PDF-card heights to the active global height window."""
    values = list(current_values or [])
    if not values:
        return [], [], []
    try:
        low, high = (float(height_range[0]), float(height_range[1]))
    except (TypeError, ValueError, IndexError):
        return [no_update] * len(values), [no_update] * len(values), [no_update] * len(values)
    if not np.isfinite(low) or not np.isfinite(high):
        return [no_update] * len(values), [no_update] * len(values), [no_update] * len(values)
    low, high = min(low, high), max(low, high)
    clamped = []
    for value in values:
        try:
            clamped.append(min(high, max(low, float(value))))
        except (TypeError, ValueError):
            clamped.append(0.5 * (low + high))
    return [low] * len(values), [high] * len(values), clamped


def _projection_fields(projection):
    config = PROJECTIONS[projection]
    if config.get("integrate_w_transport"):
        return tuple(
            dict.fromkeys(
                COMMON_PDF_FIELDS
                + tuple(
                    field_name
                    for component in (1, 2)
                    for field_name in (
                        f"rt_{component}",
                        f"thl_{component}",
                        f"chi_{component}",
                        f"stdev_chi_{component}",
                        f"varnce_rt_{component}",
                        f"varnce_thl_{component}",
                        f"corr_w_rt_{component}",
                        f"corr_w_thl_{component}",
                        f"corr_w_chi_{component}",
                        f"corr_rt_thl_{component}",
                        f"crt_{component}",
                        f"cthl_{component}",
                    )
                )
            )
        )
    return tuple(
        dict.fromkeys(
            COMMON_PDF_FIELDS
            + config["x_mean"]
            + config["x_spread"]
            + config["correlation"]
        )
    )


def _component_field(field_names, component_number):
    """Return a component-suffixed field using component 1 as the template."""
    return f"{field_names[0].rsplit('_', 1)[0]}_{int(component_number)}"


def _third_component_projection_fields(projection):
    """Fields required to locate and size component 3 in a projection."""
    config = PROJECTIONS[projection]
    return THIRD_COMPONENT_FIELDS + (
        _component_field(config["x_mean"], 3),
        _component_field(config["x_spread"], 3),
    )


def _component_weights(mixture, third_weight=None):
    """Return absolute PDF weights for a 2G or nested 3G closure.

    In CLUBB's 3G closures, ``mixt_frac`` remains the conditional split of the
    outer ADG1 pair and ``mixt_frac_3`` is the absolute third-component mass.
    """
    mixture = np.clip(np.asarray(mixture, dtype=float), 0.0, 1.0)
    if third_weight is None:
        return (mixture, 1.0 - mixture)
    third_weight = np.clip(np.asarray(third_weight, dtype=float), 0.0, 1.0)
    outer_weight = 1.0 - third_weight
    return (
        outer_weight * mixture,
        outer_weight * (1.0 - mixture),
        third_weight,
    )


def _bivariate_normal_density(x_grid, y_grid, mean_x, mean_y, var_x, var_y, correlation):
    values = np.asarray([mean_x, mean_y, var_x, var_y, correlation], dtype=float)
    if not np.all(np.isfinite(values)):
        return np.zeros_like(x_grid, dtype=float)
    var_x = max(float(var_x), 0.0)
    var_y = max(float(var_y), 0.0)
    correlation = float(np.clip(correlation, -0.999999, 0.999999))
    covariance = correlation * np.sqrt(var_x * var_y)
    matrix = np.array([[var_x, covariance], [covariance, var_y]], dtype=float)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    scale = max(float(np.max(eigenvalues)), 1.0)
    eigenvalues = np.maximum(eigenvalues, scale * 1.0e-10)
    matrix = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T
    inverse = np.linalg.inv(matrix)
    determinant = max(float(np.linalg.det(matrix)), np.finfo(float).tiny)
    dx = x_grid - float(mean_x)
    dy = y_grid - float(mean_y)
    exponent = inverse[0, 0] * dx * dx + 2.0 * inverse[0, 1] * dx * dy + inverse[1, 1] * dy * dy
    return np.exp(-0.5 * exponent) / (2.0 * np.pi * np.sqrt(determinant))


def _positive_normal_mean(mean, variance):
    """Return E[max(X, 0)] for a normally distributed field on a grid."""
    mean = np.asarray(mean, dtype=float)
    variance = max(float(variance), 0.0)
    if variance <= np.finfo(float).tiny:
        return np.maximum(mean, 0.0)
    standard_deviation = np.sqrt(variance)
    standardized = mean / standard_deviation
    normal_density = np.exp(-0.5 * standardized**2) / np.sqrt(2.0 * np.pi)
    return mean * ndtr(standardized) + standard_deviation * normal_density


def _conditional_cloud_water(
    x_grid,
    w_grid,
    *,
    mean_x,
    mean_w,
    var_x,
    var_w,
    corr_w_x,
    mean_chi,
    var_chi,
    covar_x_chi,
    covar_w_chi,
):
    """Return E[max(chi,0) | x,w] for a trivariate Gaussian component."""
    covariance_w_x = float(corr_w_x) * np.sqrt(
        max(float(var_w), 0.0) * max(float(var_x), 0.0)
    )
    hidden_covariance = np.array(
        [[float(var_x), covariance_w_x], [covariance_w_x, float(var_w)]],
        dtype=float,
    )
    hidden_covariance = 0.5 * (hidden_covariance + hidden_covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(hidden_covariance)
    floor = max(float(np.max(np.abs(eigenvalues))), 1.0) * 1.0e-10
    hidden_covariance = eigenvectors @ np.diag(np.maximum(eigenvalues, floor)) @ eigenvectors.T
    cross_covariance = np.array([covar_x_chi, covar_w_chi], dtype=float)
    regression = cross_covariance @ np.linalg.inv(hidden_covariance)
    conditional_mean = (
        float(mean_chi)
        + regression[0] * (x_grid - float(mean_x))
        + regression[1] * (w_grid - float(mean_w))
    )
    conditional_variance = max(
        float(var_chi) - float(regression @ cross_covariance), 0.0
    )
    return _positive_normal_mean(conditional_mean, conditional_variance)


def _cloud_transport_fields(
    density,
    x_grid,
    w_grid,
    projection,
    mean_w,
    mean_rc,
    conditional_cloud_water=None,
):
    """Return upward/downward cloud importance and local w'rc' contribution."""
    if projection in {"w_chi", "w_rc"}:
        cloud_water = np.maximum(x_grid, 0.0)
    elif conditional_cloud_water is not None:
        cloud_water = np.maximum(np.asarray(conditional_cloud_water, dtype=float), 0.0)
    else:
        cloud_water = np.zeros_like(density)
    w_prime = w_grid - float(mean_w)
    upward_importance = density * np.maximum(w_prime, 0.0) * cloud_water
    downward_importance = density * np.maximum(-w_prime, 0.0) * cloud_water
    signed_contribution = density * w_prime * (cloud_water - float(mean_rc))
    return upward_importance, downward_importance, signed_contribution


def _integrated_transport_summary(bundle):
    """Return mapped transport integrals in the displayed m/s·g/kg units."""
    x_values = np.asarray(bundle["x"], dtype=float)
    y_values = np.asarray(bundle["y"], dtype=float)
    dx = abs(float(x_values[1] - x_values[0])) if x_values.size > 1 else 1.0
    dy = abs(float(y_values[1] - y_values[0])) if y_values.size > 1 else 1.0
    signed = np.asarray(bundle["signed_transport"], dtype=float)
    return {
        "upward": float(np.sum(bundle["upward_importance"]) * dx * dy),
        "downward": float(np.sum(bundle["downward_importance"]) * dx * dy),
        "positive": float(np.sum(np.maximum(signed, 0.0)) * dx * dy),
        "negative": float(np.sum(np.minimum(signed, 0.0)) * dx * dy),
        "net": float(np.sum(signed) * dx * dy),
    }


def _transport_view_summary(bundle, transport_view):
    """Format only the integral relevant to the active coloring mode."""
    integrals = bundle["transport_integrals"]
    if transport_view == "upward":
        return [f"upward cloud transport = {integrals['upward']:.3e} m/s·g/kg"]
    if transport_view == "downward":
        return [f"downward cloud transport = {integrals['downward']:.3e} m/s·g/kg"]
    return [
        f"signed positive = {integrals['positive']:.3e} m/s·g/kg",
        f"signed negative = {integrals['negative']:.3e} m/s·g/kg",
        f"signed net = {integrals['net']:.3e} m/s·g/kg",
    ]


def covariance_ellipse_points(mean, covariance, radius=1.0, point_count=181):
    """Return an RMS covariance ellipse without asserting Gaussian coverage."""
    mean = np.asarray(mean, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if mean.shape != (2,) or covariance.shape != (2, 2):
        raise PdfContourError("Covariance ellipse inputs have invalid shapes.")
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(covariance)):
        raise PdfContourError("Covariance ellipse inputs contain non-finite values.")
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = max(float(np.max(np.abs(eigenvalues))), 1.0) * 1.0e-10
    if np.min(eigenvalues) < -tolerance:
        raise PdfContourError("Covariance matrix is not positive semidefinite.")
    eigenvalues = np.maximum(eigenvalues, 0.0)
    angles = np.linspace(0.0, 2.0 * np.pi, int(point_count))
    circle = np.vstack((np.cos(angles), np.sin(angles)))
    transform = eigenvectors @ np.diag(np.sqrt(eigenvalues))
    points = mean[:, np.newaxis] + float(radius) * transform @ circle
    return points[0], points[1]


def enclosed_mass_thresholds(density, x_values, y_values, masses=ENCLOSED_MASSES):
    """Return density levels enclosing the requested highest-density masses."""
    density = np.asarray(density, dtype=float)
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    if density.size == 0 or not np.any(np.isfinite(density)):
        return []
    dx = abs(float(x_values[1] - x_values[0])) if x_values.size > 1 else 1.0
    dy = abs(float(y_values[1] - y_values[0])) if y_values.size > 1 else 1.0
    flattened = np.nan_to_num(density, nan=0.0, posinf=0.0, neginf=0.0).ravel()
    order = np.argsort(flattened)[::-1]
    cumulative = np.cumsum(flattened[order]) * dx * dy
    total = float(cumulative[-1]) if cumulative.size else 0.0
    if total <= 0.0:
        return []
    thresholds = []
    for mass in masses:
        position = min(int(np.searchsorted(cumulative, float(mass) * total)), order.size - 1)
        thresholds.append(float(flattened[order[position]]))
    return thresholds


def smoothed_raw_enclosure(raw_weight, smoothing_bins, masses=ENCLOSED_MASSES):
    """Return smooth contours whose cutoffs retain raw-weight coverage.

    The Gaussian field determines the geometry and ranking, while the original
    nonnegative weights determine when 50% and 99.5% have been accumulated.
    This regularizes sparse histogram islands without redefining the reported
    probability or absolute-contribution fraction.
    """
    raw = np.nan_to_num(
        np.asarray(raw_weight, dtype=float), nan=0.0, posinf=0.0, neginf=0.0
    )
    raw = np.maximum(raw, 0.0)
    total = float(np.sum(raw))
    if raw.size == 0 or total <= 0.0:
        return raw, []
    sigma = min(max(float(smoothing_bins), 0.0), 4.0)
    smooth = (
        gaussian_filter(raw, sigma=sigma, mode="nearest", truncate=3.0)
        if sigma > 1.0e-9
        else raw.copy()
    )
    order = np.argsort(smooth.ravel())[::-1]
    cumulative = np.cumsum(raw.ravel()[order])
    thresholds = []
    for mass in masses:
        position = min(
            int(np.searchsorted(cumulative, float(mass) * total)),
            order.size - 1,
        )
        thresholds.append(float(smooth.ravel()[order[position]]))
    return smooth, thresholds


def aggregate_gaussian_component(weights, mean_x, mean_y, var_x, var_y, correlation):
    """Moment-match a time sequence of one mixture component to one ellipse."""
    weights = np.asarray(weights, dtype=float)
    mean_x = np.asarray(mean_x, dtype=float)
    mean_y = np.asarray(mean_y, dtype=float)
    var_x = np.asarray(var_x, dtype=float)
    var_y = np.asarray(var_y, dtype=float)
    correlation = np.asarray(correlation, dtype=float)
    finite = (
        np.isfinite(weights)
        & np.isfinite(mean_x)
        & np.isfinite(mean_y)
        & np.isfinite(var_x)
        & np.isfinite(var_y)
        & np.isfinite(correlation)
        & (weights >= 0.0)
    )
    if not np.any(finite) or float(np.sum(weights[finite])) <= 0.0:
        raise PdfContourError("No valid samples are available for this PDF component.")
    normalized = weights[finite] / np.sum(weights[finite])
    means = np.column_stack((mean_x[finite], mean_y[finite]))
    mean = np.sum(normalized[:, np.newaxis] * means, axis=0)
    cov_xy = np.clip(correlation[finite], -0.999999, 0.999999) * np.sqrt(
        np.maximum(var_x[finite], 0.0) * np.maximum(var_y[finite], 0.0)
    )
    covariance = np.zeros((2, 2), dtype=float)
    for sample, sample_weight in enumerate(normalized):
        within = np.array(
            [
                [max(float(var_x[finite][sample]), 0.0), float(cov_xy[sample])],
                [float(cov_xy[sample]), max(float(var_y[finite][sample]), 0.0)],
            ]
        )
        offset = means[sample] - mean
        covariance += sample_weight * (within + np.outer(offset, offset))
    denominator = np.sqrt(max(covariance[0, 0], 0.0) * max(covariance[1, 1], 0.0))
    corr = covariance[0, 1] / denominator if denominator > 0.0 else np.nan
    ellipse_x, ellipse_y = covariance_ellipse_points(mean, covariance, radius=2.0)
    return {
        "weight": float(np.mean(weights[finite])),
        "mean": mean,
        "covariance": covariance,
        "correlation": float(corr),
        "ellipse_x": ellipse_x,
        "ellipse_y": ellipse_y,
    }


def _stats_series_at_height(path, field_names, time_indices, height, col_index):
    extracted = {}
    selected_heights = []
    for field_name in field_names:
        result = shared.extract_time_height_for_path(
            path,
            field_name,
            col_index=col_index,
            column_mode="single",
        )
        if result is None:
            raise PdfContourError(f"{field_name} is required for this PDF projection.")
        _time_values, z_values, image, _units, _long_name, _z_units = result
        z_values = np.asarray(z_values, dtype=float)
        image = np.asarray(image, dtype=float)
        if image.ndim != 2 or image.shape[0] == 0 or image.shape[1] == 0:
            raise PdfContourError(f"{field_name} does not contain time-height data.")
        z_index = int(np.nanargmin(np.abs(z_values - float(height))))
        start = max(0, min(int(time_indices[0]), image.shape[0] - 1))
        end = max(0, min(int(time_indices[1]), image.shape[0] - 1))
        if end < start:
            start, end = end, start
        extracted[field_name] = np.asarray(image[start : end + 1, z_index], dtype=float)
        selected_heights.append(float(z_values[z_index]))
    count = min(len(values) for values in extracted.values())
    if count <= 0:
        raise PdfContourError("The selected PDF time window is empty.")
    return {name: values[:count] for name, values in extracted.items()}, float(np.median(selected_heights))


def _finite_axis_range(values, fallback):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return list(fallback)
    low = float(np.min(finite))
    high = float(np.max(finite))
    if high <= low:
        padding = max(abs(low) * 0.1, 1.0)
        return [low - padding, high + padding]
    padding = 0.08 * (high - low)
    return [low - padding, high + padding]


def _finite_mean(values):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else np.nan


@lru_cache(maxsize=8)
def _raw_run_inventory(run_dir):
    """Cache the expensive raw-3D directory scan for the life of the app."""
    return inventory_run(Path(run_dir))


@lru_cache(maxsize=48)
def _raw_snapshot(run_dir, elapsed_seconds, height_m):
    """Cache planes because slider revisits are common and NetCDF I/O dominates."""
    return load_snapshot(Path(run_dir), int(elapsed_seconds), float(height_m))


def _normalize_color_signal(value):
    return value if value in COLOR_SIGNALS else "wprcp"


def _matching_raw_snapshot(case_data, global_context, enabled_sources, height):
    """Return a shared raw SAM plane when one is registered for the case."""
    case_name = str((case_data or {}).get("name") or "").lower()
    run_dir = RAW_3D_RUNS.get(case_name)
    if "sam" not in enabled_sources or run_dir is None or not run_dir.exists():
        return None
    try:
        inventory = _raw_run_inventory(str(run_dir))
        initial = float((case_data or {}).get("model_time_initial_seconds") or 0.0)
        start = float(global_context.get("time_point") or initial)
        duration = float(global_context.get("time_range") or 0.0) * 60.0
        elapsed = max(0.0, start + 0.5 * duration - initial)
        step = min(inventory.steps_seconds, key=lambda value: abs(value - elapsed))
        raw_height = min(inventory.heights_m, key=lambda value: abs(value - float(height)))
        return _raw_snapshot(str(run_dir), int(step), float(raw_height))
    except (FileNotFoundError, OSError, ValueError):
        return None


def _raw_coordinates(snapshot, projection):
    def axis_values(axis):
        if axis == "w":
            return snapshot.samples[:, 0]
        if axis == "rt":
            return snapshot.samples[:, 1] * 1000.0
        if axis == "thl":
            return snapshot.samples[:, 2]
        if axis == "chi":
            return snapshot.chi_samples * 1000.0
        return snapshot.rc_samples * 1000.0

    x_axis, y_axis = PROJECTIONS[projection]["axes"]
    return np.asarray(axis_values(x_axis), dtype=float), np.asarray(axis_values(y_axis), dtype=float)


def _raw_axis_extent(snapshot, projection):
    x_values, y_values = _raw_coordinates(snapshot, projection)
    if projection == "w_rc":
        x_low = 0.0
        x_high = float(np.percentile(x_values, 99.9))
    else:
        x_low, x_high = np.percentile(x_values, (0.2, 99.8))
    y_low, y_high = np.percentile(y_values, (0.2, 99.8))
    cloudy = np.asarray(snapshot.rc_samples) > 0.0
    if np.any(cloudy):
        x_low = min(x_low, float(np.min(x_values[cloudy])))
        x_high = max(x_high, float(np.max(x_values[cloudy])))
        y_low = min(y_low, float(np.min(y_values[cloudy])))
        y_high = max(y_high, float(np.max(y_values[cloudy])))
    return float(x_low), float(x_high), float(y_low), float(y_high)


def _centers_to_edges(values):
    values = np.asarray(values, dtype=float)
    midpoints = 0.5 * (values[:-1] + values[1:])
    return np.concatenate(
        ([values[0] - (midpoints[0] - values[0])], midpoints, [values[-1] + (values[-1] - midpoints[-1])])
    )


def _raw_histogram_geometry(snapshot, projection, x_values, y_values):
    """Cache bin membership; color switches then need only one weighted sum."""
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    key = (
        str(getattr(snapshot, "source_file", f"in-memory-{id(snapshot)}")),
        int(getattr(snapshot, "elapsed_seconds", round(float(snapshot.elapsed_minutes) * 60.0))),
        round(float(snapshot.height_m), 6),
        projection,
        x_values.tobytes(),
        y_values.tobytes(),
    )
    cached = _RAW_HISTOGRAM_CACHE.get(key)
    if cached is not None:
        _RAW_HISTOGRAM_CACHE.move_to_end(key)
        return cached
    raw_x, raw_w = _raw_coordinates(snapshot, projection)
    x_edges = _centers_to_edges(x_values)
    y_edges = _centers_to_edges(y_values)
    x_bin = np.searchsorted(x_edges, raw_x, side="right") - 1
    y_bin = np.searchsorted(y_edges, raw_w, side="right") - 1
    # numpy.histogram2d includes the final right edge in the final bin.
    x_bin = np.where(raw_x == x_edges[-1], x_values.size - 1, x_bin)
    y_bin = np.where(raw_w == y_edges[-1], y_values.size - 1, y_bin)
    valid = (
        np.isfinite(raw_x)
        & np.isfinite(raw_w)
        & (x_bin >= 0)
        & (x_bin < x_values.size)
        & (y_bin >= 0)
        & (y_bin < y_values.size)
    )
    sample_indices = np.flatnonzero(valid)
    flat_bins = y_bin[valid] * x_values.size + x_bin[valid]
    probability = np.bincount(
        flat_bins, minlength=x_values.size * y_values.size
    ).reshape(y_values.size, x_values.size) / raw_x.size
    cached = {
        "sample_indices": sample_indices,
        "flat_bins": flat_bins,
        "probability": probability,
        "sample_count": raw_x.size,
    }
    _RAW_HISTOGRAM_CACHE[key] = cached
    _RAW_HISTOGRAM_CACHE.move_to_end(key)
    while len(_RAW_HISTOGRAM_CACHE) > _RAW_HISTOGRAM_CACHE_MAX_ENTRIES:
        _RAW_HISTOGRAM_CACHE.popitem(last=False)
    return cached


def _raw_color_signal_samples(snapshot, color_signal):
    """Return exact raw-plane local contributions for one displayed signal."""
    color_signal = _normalize_color_signal(color_signal)
    w = np.asarray(snapshot.samples[:, 0], dtype=float)
    w_prime = w - float(snapshot.mean[0])
    rt = np.asarray(snapshot.samples[:, 1], dtype=float) * 1000.0
    rt_prime = rt - float(np.mean(rt))
    if color_signal == "wprcp":
        rc = np.asarray(snapshot.rc_samples, dtype=float) * 1000.0
        return w_prime * (rc - float(np.mean(rc)))
    if color_signal == "wprtp":
        return w_prime * rt_prime
    if color_signal == "wpchi":
        chi = np.asarray(snapshot.chi_samples, dtype=float) * 1000.0
        return w_prime * (chi - float(np.mean(chi)))
    if color_signal == "wpthlp":
        thl = np.asarray(snapshot.samples[:, 2], dtype=float)
        return w_prime * (thl - float(np.mean(thl)))
    return w_prime * np.square(rt_prime)


def _raw_histogram_bundle(
    snapshot,
    projection,
    x_values,
    y_values,
    color_signal="wprcp",
    contour_smoothing_bins=DEFAULT_RAW_CONTOUR_SMOOTHING_BINS,
):
    """Build a probability/signed-moment raster from one cached raw SAM plane."""
    color_signal = _normalize_color_signal(color_signal)
    geometry = _raw_histogram_geometry(snapshot, projection, x_values, y_values)
    probability_only = bool(COLOR_SIGNALS[color_signal].get("probability_only"))
    if probability_only:
        weighted = np.zeros_like(geometry["probability"])
    else:
        local_signal = _raw_color_signal_samples(snapshot, color_signal)
        weighted = np.bincount(
            geometry["flat_bins"],
            weights=local_signal[geometry["sample_indices"]],
            minlength=len(x_values) * len(y_values),
        ).reshape(len(y_values), len(x_values)) / geometry["sample_count"]
    enclosed_weight = (
        geometry["probability"] if probability_only else np.abs(weighted)
    )
    enclosed_label = (
        RAW_ENCLOSED_PROBABILITY_LABEL
        if probability_only
        else (
            "Raw SAM enclosed "
            f"|{COLOR_SIGNALS[color_signal]['label']}| (50/80/99.5%)"
        )
    )
    smooth_enclosure, enclosure_thresholds = smoothed_raw_enclosure(
        enclosed_weight, contour_smoothing_bins
    )
    return {
        "x": np.asarray(x_values, dtype=float),
        "y": np.asarray(y_values, dtype=float),
        "total": geometry["probability"],
        "signed_transport": weighted,
        "color_signal": color_signal,
        "color_signal_label": COLOR_SIGNALS[color_signal]["label"],
        "color_signal_units": COLOR_SIGNALS[color_signal]["units"],
        "color_available": not probability_only,
        "probability_only": probability_only,
        # A signed local contribution cannot define an enclosed fraction.  Its
        # magnitude does: each contour therefore contains the stated fraction
        # of total absolute contribution, not a signed cancellation.
        "enclosed_mass_field": smooth_enclosure,
        "enclosed_mass_thresholds": enclosure_thresholds,
        "enclosed_mass_label": enclosed_label,
        "enclosed_smoothing_bins": min(
            max(float(contour_smoothing_bins), 0.0), 4.0
        ),
        "raw_sampled": True,
        "raw_time_minutes": float(snapshot.elapsed_minutes),
        "raw_height": float(snapshot.height_m),
    }


def _analytic_background_bundle(bundle, color_signal):
    """Use only the exact analytic cloud-water color; leave other signals gold."""
    color_signal = _normalize_color_signal(color_signal)
    result = dict(bundle)
    result["color_signal"] = color_signal
    result["color_signal_label"] = COLOR_SIGNALS[color_signal]["label"]
    result["color_signal_units"] = COLOR_SIGNALS[color_signal]["units"]
    result["color_available"] = color_signal == "wprcp"
    result["probability_only"] = (
        bool(COLOR_SIGNALS[color_signal].get("probability_only"))
        or not result["color_available"]
    )
    if not result["color_available"]:
        result["signed_transport"] = np.zeros_like(bundle["total"])
    return result


def _component_trivariate_statistics(series, index, component):
    """Return one component's w/r_t/theta_l/chi moments in plotted units."""
    suffix = str(component)
    means = {
        "w": float(series[f"w_{suffix}"][index]),
        "rt": float(series[f"rt_{suffix}"][index] * 1000.0),
        "thl": float(series[f"thl_{suffix}"][index]),
        "chi": float(series[f"chi_{suffix}"][index] * 1000.0),
    }
    variances = {
        "w": float(series[f"varnce_w_{suffix}"][index]),
        "rt": float(series[f"varnce_rt_{suffix}"][index] * 1.0e6),
        "thl": float(series[f"varnce_thl_{suffix}"][index]),
        "chi": float((series[f"stdev_chi_{suffix}"][index] * 1000.0) ** 2),
    }
    covariances = {
        ("w", "rt"): float(series[f"corr_w_rt_{suffix}"][index]) * np.sqrt(max(variances["w"], 0.0) * max(variances["rt"], 0.0)),
        ("w", "thl"): float(series[f"corr_w_thl_{suffix}"][index]) * np.sqrt(max(variances["w"], 0.0) * max(variances["thl"], 0.0)),
        ("w", "chi"): float(series[f"corr_w_chi_{suffix}"][index]) * np.sqrt(max(variances["w"], 0.0) * max(variances["chi"], 0.0)),
        ("rt", "thl"): float(series[f"corr_rt_thl_{suffix}"][index]) * np.sqrt(max(variances["rt"], 0.0) * max(variances["thl"], 0.0)),
    }
    covariances[("rt", "chi")] = (
        float(series[f"crt_{suffix}"][index]) * variances["rt"]
        - float(series[f"cthl_{suffix}"][index]) * covariances[("rt", "thl")]
    )
    covariances[("thl", "chi")] = (
        float(series[f"crt_{suffix}"][index]) * covariances[("rt", "thl")]
        - float(series[f"cthl_{suffix}"][index]) * variances["thl"]
    )
    return means, variances, covariances


def _covariance(covariances, variances, left, right):
    """Fetch an unordered covariance, including a variance on the diagonal."""
    if left == right:
        return variances[left]
    return covariances.get((left, right), covariances.get((right, left), np.nan))


def _integrated_wprcp_field(density, x_grid, y_grid, axes, means, variances, covariances, mean_w, mean_rc):
    """Integrate local w'r_c' over w for a selected two-dimensional marginal.

    Conditioning the component's joint (w, chi) Gaussian on the plotted axes
    gives an exact one-dimensional w integral for the Gaussian cloud-water
    diagnostic used by this display.  This is the analytic counterpart to
    summing raw-SAM w'r_c' samples into the same r_t/theta_l or chi/theta_l bin.
    """
    first, second = axes
    covariance_z = np.array(
        [
            [_covariance(covariances, variances, first, first), _covariance(covariances, variances, first, second)],
            [_covariance(covariances, variances, second, first), _covariance(covariances, variances, second, second)],
        ],
        dtype=float,
    )
    values = np.asarray((means["w"], means["chi"], mean_w, mean_rc, *covariance_z.ravel()), dtype=float)
    if not np.all(np.isfinite(values)):
        return np.zeros_like(density)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_z)
    scale = max(float(np.max(eigenvalues)), 1.0)
    covariance_z = eigenvectors @ np.diag(np.maximum(eigenvalues, scale * 1.0e-10)) @ eigenvectors.T
    inverse_z = np.linalg.inv(covariance_z)
    delta = np.stack((x_grid - means[first], y_grid - means[second]), axis=0)
    covariance_w_z = np.asarray(
        [_covariance(covariances, variances, "w", first), _covariance(covariances, variances, "w", second)],
        dtype=float,
    )
    covariance_chi_z = np.asarray(
        [_covariance(covariances, variances, "chi", first), _covariance(covariances, variances, "chi", second)],
        dtype=float,
    )
    w_gain = covariance_w_z @ inverse_z
    chi_gain = covariance_chi_z @ inverse_z
    mean_w_given_z = means["w"] + np.einsum("i,ijk->jk", w_gain, delta)
    mean_chi_given_z = means["chi"] + np.einsum("i,ijk->jk", chi_gain, delta)
    variance_chi_given_z = max(
        variances["chi"] - covariance_chi_z @ inverse_z @ covariance_chi_z,
        0.0,
    )
    covariance_w_chi_given_z = (
        _covariance(covariances, variances, "w", "chi")
        - covariance_w_z @ inverse_z @ covariance_chi_z
    )
    cloud_water = _positive_normal_mean(mean_chi_given_z, variance_chi_given_z)
    if variance_chi_given_z > np.finfo(float).tiny:
        cloud_probability = ndtr(mean_chi_given_z / np.sqrt(variance_chi_given_z))
    else:
        cloud_probability = (mean_chi_given_z > 0.0).astype(float)
    contribution = (
        (mean_w_given_z - float(mean_w)) * (cloud_water - float(mean_rc))
        + covariance_w_chi_given_z * cloud_probability
    )
    return density * contribution


def _conditional_gaussian_signal_field(
    density,
    x_grid,
    y_grid,
    axes,
    means,
    variances,
    covariances,
    grid_means,
    color_signal,
):
    """Return a binned local moment contribution of one Gaussian component.

    The raw-SAM contour weights are absolute local contributions.  This is the
    analytic counterpart: condition the full Gaussian on the displayed pair,
    then take the appropriate conditional product moment before multiplying by
    that pair's density.
    """
    color_signal = _normalize_color_signal(color_signal)
    if color_signal == "probability":
        return np.asarray(density, dtype=float)
    if color_signal == "wprcp":
        conditioned_axes = tuple("chi" if axis == "rc" else axis for axis in axes)
        return _integrated_wprcp_field(
            density, x_grid, y_grid, conditioned_axes, means, variances, covariances,
            grid_means["w"], grid_means["rc"],
        )

    first, second = ("chi" if axis == "rc" else axis for axis in axes)
    covariance_z = np.asarray(
        (
            (_covariance(covariances, variances, first, first), _covariance(covariances, variances, first, second)),
            (_covariance(covariances, variances, second, first), _covariance(covariances, variances, second, second)),
        ),
        dtype=float,
    )
    if not np.all(np.isfinite(covariance_z)):
        return np.zeros_like(density)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_z)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    covariance_z = eigenvectors @ np.diag(np.maximum(eigenvalues, scale * 1.0e-10)) @ eigenvectors.T
    inverse_z = np.linalg.inv(covariance_z)
    delta = np.stack((x_grid - means[first], y_grid - means[second]), axis=0)

    def conditional(variable):
        covariance_variable_z = np.asarray(
            (
                _covariance(covariances, variances, variable, first),
                _covariance(covariances, variances, variable, second),
            ),
            dtype=float,
        )
        gain = covariance_variable_z @ inverse_z
        mean = means[variable] + np.einsum("i,ijk->jk", gain, delta)
        variance = max(
            variances[variable] - covariance_variable_z @ inverse_z @ covariance_variable_z,
            0.0,
        )
        return mean, variance, covariance_variable_z

    target = {"wprtp": "rt", "wpchi": "chi", "wpthlp": "thl", "wprtp2": "rt"}.get(color_signal)
    if target is None:
        return np.zeros_like(density)
    mean_w, _variance_w, covariance_w_z = conditional("w")
    mean_target, variance_target, covariance_target_z = conditional(target)
    covariance_w_target = (
        _covariance(covariances, variances, "w", target)
        - covariance_w_z @ inverse_z @ covariance_target_z
    )
    w_prime = mean_w - grid_means["w"]
    target_prime = mean_target - grid_means[target]
    if color_signal == "wprtp2":
        contribution = (
            w_prime * (target_prime**2 + variance_target)
            + 2.0 * target_prime * covariance_w_target
        )
    else:
        contribution = w_prime * target_prime + covariance_w_target
    return density * contribution


def _clubb_weighted_enclosure(
    x_values, y_values, axes, component_times, color_signal
):
    """Return a CLUBB local-signal field and matching enclosed-mass metadata."""
    color_signal = _normalize_color_signal(color_signal)
    if color_signal == "probability" or not component_times:
        return None, None, ENCLOSED_PROBABILITY_LABEL
    x_grid, y_grid = np.meshgrid(x_values, y_values)
    weighted_field = np.zeros_like(x_grid)
    valid_count = 0
    for components in component_times:
        weights = np.asarray([weight for weight, _means, _variances, _covariances in components], dtype=float)
        if not np.all(np.isfinite(weights)) or np.sum(weights) <= 0.0:
            continue
        grid_means = {
            variable: float(sum(weight * means[variable] for weight, means, _variances, _covariances in components))
            for variable in ("w", "rt", "thl", "chi")
        }
        grid_means["rc"] = float(sum(
            weight * _positive_normal_mean(means["chi"], variances["chi"])
            for weight, means, variances, _covariances in components
        ))
        for weight, means, variances, covariances in components:
            first, second = ("chi" if axis == "rc" else axis for axis in axes)
            covariance_xy = _covariance(covariances, variances, first, second)
            correlation_xy = covariance_xy / np.sqrt(
                max(variances[first], 0.0) * max(variances[second], 0.0)
            )
            density = _bivariate_normal_density(
                x_grid, y_grid, means[first], means[second], variances[first],
                variances[second], correlation_xy,
            )
            weighted_field += _conditional_gaussian_signal_field(
                weight * density, x_grid, y_grid, axes, means, variances,
                covariances, grid_means, color_signal,
            )
        valid_count += 1
    if valid_count == 0 or not np.any(np.isfinite(weighted_field)):
        return None, None, ENCLOSED_PROBABILITY_LABEL
    weighted_field /= valid_count
    enclosed_field = np.abs(weighted_field)
    thresholds = enclosed_mass_thresholds(enclosed_field, x_values, y_values)
    if not thresholds:
        return None, None, ENCLOSED_PROBABILITY_LABEL
    return enclosed_field, thresholds, CLUBB_ENCLOSED_WEIGHTED_LABEL.format(
        COLOR_SIGNALS[color_signal]["label"]
    )


class PdfContourPlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="pdf_contour",
            default_vars=["w_chi"],
            var_input_type="pdf-projection",
            graph_type="pdf-contour-graph",
            case_data_var_key="profile_vars",
            subtitle=(
                "Gold shows probability, blue/red show positive/negative local w′r꜀′ transport, "
                "and pale colors show where probability and transport overlap. When matching raw SAM 3-D "
                "data are available, that sampled plane replaces the analytic CLUBB heatmap. "
                "The local height applies to this card; global time and selected column still apply."
            ),
        )

    def case_data_options(self, case_data):
        case_data = case_data or {}
        available = {
            option.get("value")
            for option in case_data.get("profile_vars") or []
            if option.get("value")
        }
        return [
            option
            for option in PROJECTION_OPTIONS
            if set(_projection_fields(option["value"])).issubset(available)
        ]

    def get_variable_options(self, collection, case_data):
        return self.case_data_options(case_data)

    def clear_render_state(self):
        """Release raw-plane and analytic rasters when the Plot context changes."""
        _DENSITY_CACHE.clear()
        _RAW_HISTOGRAM_CACHE.clear()
        _raw_snapshot.cache_clear()
        _raw_run_inventory.cache_clear()

    def make_default_state(self, case_data, plot_id):
        case_data = case_data or {}
        height_range = case_data.get("default_height_range") or [0.0, 1.0]
        return {
            "plot_type": self.plot_type_id,
            "var": "w_chi",
            "height": 0.5 * (float(height_range[0]) + float(height_range[1])),
            "transport_view": "signed",
            "color_signal": "wprcp",
            "contour_smoothing_bins": DEFAULT_RAW_CONTOUR_SMOOTHING_BINS,
            "size": "large",
        }

    def height_input_id(self, plot_id):
        return {"type": "pdf-height", "index": plot_id}

    def transport_view_input_id(self, plot_id):
        return {"type": "pdf-transport-view", "index": plot_id}

    def color_signal_input_id(self, plot_id):
        return {"type": "pdf-color-signal", "index": plot_id}

    def contour_smoothing_input_id(self, plot_id):
        return {"type": "pdf-contour-smoothing", "index": plot_id}

    def render_card(self, plot_id, state, global_context):
        case_data = global_context.get("case_data") or {}
        options = self.case_data_options(case_data)
        projection = state.get("var") or "w_chi"
        transport_view = "signed"
        if transport_view not in TRANSPORT_VIEWS:
            transport_view = "upward"
        height_min = float(case_data.get("height_slider_min") or 0.0)
        height_max = float(case_data.get("height_slider_max") or max(height_min, 1.0))
        if height_max < height_min:
            height_min, height_max = height_max, height_min
        default_height = 0.5 * (height_min + height_max)
        height = min(height_max, max(height_min, float(state.get("height", default_height))))
        height_step = max(float(case_data.get("height_step") or 1.0), 1.0e-6)
        controls = html.Div(
            [
                styled_dropdown(
                    id=self.var_input_id(plot_id),
                    options=options,
                    value=projection,
                    clearable=False,
                ),
                dcc.Input(
                    id=self.transport_view_input_id(plot_id),
                    value="signed",
                    type="hidden",
                ),
                html.Div(
                    [
                        html.Label("SAM color signal", className="pdf-color-signal-label"),
                        dcc.RadioItems(
                            id=self.color_signal_input_id(plot_id),
                            options=COLOR_SIGNAL_OPTIONS,
                            value=_normalize_color_signal(state.get("color_signal")),
                            inline=True,
                            className="pdf-color-signal-switch",
                        ),
                    ],
                    className="pdf-inline-control-row pdf-color-signal-control",
                ),
                html.Div(
                    [
                        html.Label(
                            "SAM contour smoothing [bins]",
                            className="pdf-contour-smoothing-label",
                        ),
                        dcc.Slider(
                            id=self.contour_smoothing_input_id(plot_id),
                            min=0.0,
                            max=3.0,
                            step=0.25,
                            value=min(
                                max(
                                    float(
                                        state.get(
                                            "contour_smoothing_bins",
                                            DEFAULT_RAW_CONTOUR_SMOOTHING_BINS,
                                        )
                                    ),
                                    0.0,
                                ),
                                3.0,
                            ),
                            marks=None,
                            tooltip={"placement": "bottom", "always_visible": False},
                            className="pdf-compact-slider",
                        ),
                    ],
                    className="pdf-inline-control-row pdf-contour-smoothing-control",
                ),
                html.Div(
                    [
                        html.Label("Height [m]", className="pdf-compact-control-label"),
                        dcc.Slider(
                            id=self.height_input_id(plot_id),
                            min=height_min,
                            max=height_max,
                            step=height_step,
                            value=height,
                            marks=None,
                            tooltip={"placement": "bottom", "always_visible": False},
                            className="pdf-compact-slider",
                        ),
                    ],
                    className="pdf-inline-control-row pdf-height-control",
                ),
            ],
            className="pdf-contour-controls",
        )
        size_value = shared.normalize_plot_size(state.get("size") or "large")
        size_text, size_class = shared.plot_size_button_props(size_value)
        size_button = html.Button(
            size_text,
            id=self.size_toggle_id(plot_id),
            className=size_class,
            title="Toggle plot size",
        )
        return shared.make_plot_card(
            subtitle=self.card_subtitle(case_data),
            help_button_id=self.help_button_id(plot_id),
            controls=controls,
            size_button=size_button,
            size_value=size_value,
            size_store_id=self.size_store_id(plot_id),
            graph_id=self.graph_id(plot_id),
            graph_shell_id=self.graph_shell_id(plot_id),
            card_id=self.card_id(plot_id),
            close_button_id=self.close_button_id(plot_id),
            render_signal_id=self.render_signal_id(plot_id),
        )

    def _paths_for_projection(self, files, projection):
        required = _projection_fields(projection)
        paths = []
        unreadable_output = False
        for path in files or []:
            metadata = shared.readable_dataset_metadata_for_path(path)
            if metadata is None:
                unreadable_output = True
                continue
            if all(name in metadata["vars"] for name in required):
                paths.append(path)
        return paths, unreadable_output

    def _integrated_w_density_bundle(
        self,
        path,
        projection,
        height,
        case_data,
        global_context,
        raw_extent,
        color_signal,
    ):
        """Build a non-w marginal and analytically integrate w'r_c' in each bin."""
        config = PROJECTIONS[projection]
        time_indices = shared.selected_time_indices(
            case_data,
            global_context.get("time_range"),
            global_context.get("time_point"),
            global_context.get("time_mode") or "range",
        )
        available_fields = shared.dataset_metadata_for_path(path)["vars"]
        component_fields = tuple(
            field_name
            for component in (1, 2)
            for field_name in (
                f"w_{component}",
                f"rt_{component}",
                f"thl_{component}",
                f"chi_{component}",
                f"stdev_chi_{component}",
                f"varnce_w_{component}",
                f"varnce_rt_{component}",
                f"varnce_thl_{component}",
                f"corr_w_rt_{component}",
                f"corr_w_thl_{component}",
                f"corr_w_chi_{component}",
                f"corr_rt_thl_{component}",
                f"crt_{component}",
                f"cthl_{component}",
            )
        )
        third_component_fields = tuple(
            field_name
            for component in (3,)
            for field_name in (
                "mixt_frac_3",
                f"w_{component}",
                f"rt_{component}",
                f"thl_{component}",
                f"chi_{component}",
                f"stdev_chi_{component}",
                f"varnce_w_{component}",
                f"varnce_rt_{component}",
                f"varnce_thl_{component}",
                f"corr_w_rt_{component}",
                f"corr_w_thl_{component}",
                f"corr_w_chi_{component}",
                f"corr_rt_thl_{component}",
                f"crt_{component}",
                f"cthl_{component}",
            )
        )
        has_third_component = all(name in available_fields for name in third_component_fields)
        fields = list(dict.fromkeys(component_fields + third_component_fields if has_third_component else component_fields))
        fields.extend(name for name in ("mixt_frac", "wm_zt", "rcm") if name in available_fields)
        series, selected_height = _stats_series_at_height(
            path,
            fields,
            time_indices,
            height,
            int(global_context.get("selected_column") or 0),
        )
        mixture = np.clip(series["mixt_frac"], 0.0, 1.0)
        third_weight = np.clip(series["mixt_frac_3"], 0.0, 1.0) if has_third_component else None
        if third_weight is not None and not np.any(third_weight > 1.0e-12):
            third_weight = None
            has_third_component = False
        component_numbers = (1, 2, 3) if has_third_component else (1, 2)
        component_weights = _component_weights(mixture, third_weight)
        axes = config["axes"]
        statistics = [
            [_component_trivariate_statistics(series, index, component) for component in component_numbers]
            for index in range(len(mixture))
        ]
        extent_x = np.asarray(raw_extent[:2], dtype=float) if raw_extent is not None else np.empty(0)
        extent_y = np.asarray(raw_extent[2:], dtype=float) if raw_extent is not None else np.empty(0)
        x_candidates = np.asarray(
            [
                candidate
                for values in statistics
                for means, variances, _covariances in values
                for candidate in (
                    means[axes[0]] - 4.0 * np.sqrt(max(variances[axes[0]], 0.0)),
                    means[axes[0]] + 4.0 * np.sqrt(max(variances[axes[0]], 0.0)),
                )
            ]
            + extent_x.tolist(),
            dtype=float,
        )
        y_candidates = np.asarray(
            [
                candidate
                for values in statistics
                for means, variances, _covariances in values
                for candidate in (
                    means[axes[1]] - 4.0 * np.sqrt(max(variances[axes[1]], 0.0)),
                    means[axes[1]] + 4.0 * np.sqrt(max(variances[axes[1]], 0.0)),
                )
            ]
            + extent_y.tolist(),
            dtype=float,
        )
        x_range = _finite_axis_range(x_candidates, [-1.0, 1.0])
        y_range = _finite_axis_range(y_candidates, [-1.0, 1.0])
        x_values = np.linspace(x_range[0], x_range[1], PDF_GRID_POINTS)
        y_values = np.linspace(y_range[0], y_range[1], PDF_GRID_POINTS)
        x_grid, y_grid = np.meshgrid(x_values, y_values)
        total_density = np.zeros_like(x_grid)
        signed_transport = np.zeros_like(x_grid)
        component_densities = [np.zeros_like(x_grid) for _component in component_numbers]
        valid_count = 0
        for index, component_stats in enumerate(statistics):
            weights = [weight[index] for weight in component_weights]
            if not np.all(np.isfinite(weights)):
                continue
            mean_w = float(series["wm_zt"][index]) if "wm_zt" in series else float(sum(weight * values[0]["w"] for weight, values in zip(weights, component_stats)))
            mean_rc = float(series["rcm"][index] * 1000.0) if "rcm" in series else float(sum(weight * _positive_normal_mean(values[0]["chi"], values[1]["chi"]) for weight, values in zip(weights, component_stats)))
            for component_index, (weight, values) in enumerate(zip(weights, component_stats)):
                means, variances, covariances = values
                covariance_xy = _covariance(covariances, variances, axes[0], axes[1])
                correlation_xy = covariance_xy / np.sqrt(max(variances[axes[0]], 0.0) * max(variances[axes[1]], 0.0))
                density = _bivariate_normal_density(
                    x_grid,
                    y_grid,
                    means[axes[0]],
                    means[axes[1]],
                    variances[axes[0]],
                    variances[axes[1]],
                    correlation_xy,
                )
                weighted_density = weight * density
                total_density += weighted_density
                component_densities[component_index] += weighted_density
                signed_transport += weight * _integrated_wprcp_field(
                    density,
                    x_grid,
                    y_grid,
                    axes,
                    means,
                    variances,
                    covariances,
                    mean_w,
                    mean_rc,
                )
            valid_count += 1
        if valid_count == 0:
            raise PdfContourError("No valid PDF parameters were found in the selected window.")
        total_density /= valid_count
        signed_transport /= valid_count
        component_densities = [density / valid_count for density in component_densities]
        components = []
        component_means = []
        for component_index in range(len(component_numbers)):
            means = np.asarray([values[component_index][0][axes[0]] for values in statistics])
            means_y = np.asarray([values[component_index][0][axes[1]] for values in statistics])
            variances_x = np.asarray([values[component_index][1][axes[0]] for values in statistics])
            variances_y = np.asarray([values[component_index][1][axes[1]] for values in statistics])
            correlations = np.asarray([
                _covariance(values[component_index][2], values[component_index][1], axes[0], axes[1])
                / np.sqrt(max(values[component_index][1][axes[0]], 0.0) * max(values[component_index][1][axes[1]], 0.0))
                for values in statistics
            ])
            try:
                components.append(aggregate_gaussian_component(component_weights[component_index], means, means_y, variances_x, variances_y, correlations))
            except PdfContourError:
                components.append(None)
            component_means.append((_finite_mean(means), _finite_mean(means_y)))
        component_times = [
            [
                (component_weights[component_index][time_index], *component_stats[component_index])
                for component_index in range(len(component_numbers))
            ]
            for time_index, component_stats in enumerate(statistics)
        ]
        enclosed_field, enclosed_thresholds, enclosed_label = _clubb_weighted_enclosure(
            x_values, y_values, axes, component_times, color_signal
        )
        if enclosed_field is None:
            enclosed_field = total_density
            enclosed_thresholds = enclosed_mass_thresholds(
                total_density, x_values, y_values
            )
        bundle = {
            "x": x_values,
            "y": y_values,
            "total": total_density,
            "signed_transport": signed_transport,
            "upward_importance": np.zeros_like(total_density),
            "downward_importance": np.zeros_like(total_density),
            "components": components,
            "component_means": tuple(component_means),
            "mass_thresholds": enclosed_thresholds,
            "enclosed_mass_field": enclosed_field,
            "enclosed_mass_label": enclosed_label,
            "selected_height": selected_height,
            "summary": [f"z = {selected_height:.1f} m", "w′r꜀′ is integrated over w at each plotted thermodynamic bin."],
        }
        bundle.update({f"component_{index + 1}": density for index, density in enumerate(component_densities)})
        bundle["transport_integrals"] = _integrated_transport_summary(bundle)
        return bundle

    def _density_bundle(
        self, path, projection, height, case_data, global_context, raw_extent=None,
        color_signal="wprcp",
    ):
        config = PROJECTIONS[projection]
        time_indices = shared.selected_time_indices(
            case_data,
            global_context.get("time_range"),
            global_context.get("time_point"),
            global_context.get("time_mode") or "range",
        )
        try:
            stat = Path(path).stat()
            file_signature = (stat.st_mtime_ns, stat.st_size)
        except OSError:
            file_signature = (None, None)
        extent_key = tuple(round(float(value), 8) for value in raw_extent or ())
        cache_key = (
            str(path),
            file_signature,
            projection,
            round(float(height), 6),
            tuple(int(value) for value in time_indices),
            int(global_context.get("selected_column") or 0),
            extent_key,
            _normalize_color_signal(color_signal),
        )
        cached = _DENSITY_CACHE.get(cache_key)
        if cached is not None:
            _DENSITY_CACHE.move_to_end(cache_key)
            return cached
        if config.get("integrate_w_transport"):
            result = self._integrated_w_density_bundle(
                path, projection, height, case_data, global_context, raw_extent,
                color_signal,
            )
            _DENSITY_CACHE[cache_key] = result
            _DENSITY_CACHE.move_to_end(cache_key)
            while len(_DENSITY_CACHE) > _DENSITY_CACHE_MAX_ENTRIES:
                _DENSITY_CACHE.popitem(last=False)
            return result
        fields = list(_projection_fields(projection))
        available_fields = shared.dataset_metadata_for_path(path)["vars"]
        third_projection_fields = _third_component_projection_fields(projection)
        has_third_component = all(
            field_name in available_fields for field_name in third_projection_fields
        )
        if has_third_component:
            fields.extend(third_projection_fields)
            third_correlation_field = _component_field(config["correlation"], 3)
            if third_correlation_field in available_fields:
                fields.append(third_correlation_field)
            elif projection in {"w_chi", "w_rc"}:
                fields.extend(
                    field_name
                    for field_name in (
                        "corr_w_rt_3",
                        "corr_w_thl_3",
                        "varnce_rt_3",
                        "varnce_thl_3",
                        "crt_1",
                        "crt_2",
                        "cthl_1",
                        "cthl_2",
                    )
                    if field_name in available_fields
                )
        if _normalize_color_signal(color_signal) != "probability":
            contour_components = (1, 2, 3) if has_third_component else (1, 2)
            for component in contour_components:
                component_fields = _trivariate_component_field_names(component)
                if all(field_name in available_fields for field_name in component_fields):
                    fields.extend(component_fields)
        for field_name in ("wm_zt", "rcm"):
            if field_name in available_fields:
                fields.append(field_name)
        if projection in {"w_chi", "w_rc"}:
            cloud_components = (1, 2, 3) if has_third_component else (1, 2)
            for field_name in (
                "wm_zt",
                "rcm",
                *(f"rc_{component}" for component in cloud_components),
                *(f"cloud_frac_{component}" for component in cloud_components),
            ):
                if field_name in available_fields:
                    fields.append(field_name)
        elif projection == "w_rt":
            conditioning_fields = list(CLUBB_CLOUD_CONDITIONING_FIELDS)
            if has_third_component:
                conditioning_fields.extend(
                    field_name
                    for field_name in (
                        "chi_3",
                        "stdev_chi_3",
                        "corr_w_chi_3",
                        "corr_w_thl_3",
                        "varnce_thl_3",
                        "corr_rt_thl_3",
                        "crt_3",
                        "cthl_3",
                    )
                    if field_name in available_fields
                )
            for field_name in conditioning_fields:
                if field_name in available_fields:
                    fields.append(field_name)
        fields = list(dict.fromkeys(fields))
        series, selected_height = _stats_series_at_height(
            path,
            fields,
            time_indices,
            height,
            int(global_context.get("selected_column") or 0),
        )
        x_mean_1 = series[config["x_mean"][0]] * 1000.0
        x_mean_2 = series[config["x_mean"][1]] * 1000.0
        if config["spread_is_stdev"]:
            x_var_1 = np.square(series[config["x_spread"][0]] * 1000.0)
            x_var_2 = np.square(series[config["x_spread"][1]] * 1000.0)
        else:
            x_var_1 = series[config["x_spread"][0]] * 1.0e6
            x_var_2 = series[config["x_spread"][1]] * 1.0e6
        y_mean_1 = series["w_1"]
        y_mean_2 = series["w_2"]
        y_var_1 = series["varnce_w_1"]
        y_var_2 = series["varnce_w_2"]
        corr_1 = series[config["correlation"][0]]
        corr_2 = series[config["correlation"][1]]
        mixture = np.clip(series["mixt_frac"], 0.0, 1.0)
        component_x_means = [x_mean_1, x_mean_2]
        component_x_vars = [x_var_1, x_var_2]
        component_y_means = [y_mean_1, y_mean_2]
        component_y_vars = [y_var_1, y_var_2]
        component_correlations = [corr_1, corr_2]
        third_weight = None
        used_legacy_third_correlation = False
        if has_third_component:
            third_weight = np.clip(series["mixt_frac_3"], 0.0, 1.0)
            has_third_component = bool(
                np.any(np.isfinite(third_weight) & (third_weight > 1.0e-12))
            )
        if has_third_component:
            if "corr_w_chi_3" not in series and all(
                field_name in series
                for field_name in (
                    "corr_w_rt_3",
                    "corr_w_thl_3",
                    "varnce_rt_3",
                    "varnce_thl_3",
                    "stdev_chi_3",
                    "crt_1",
                    "crt_2",
                    "cthl_1",
                    "cthl_2",
                )
            ):
                series["crt_3"] = (
                    mixture * series["crt_1"]
                    + (1.0 - mixture) * series["crt_2"]
                )
                series["cthl_3"] = (
                    mixture * series["cthl_1"]
                    + (1.0 - mixture) * series["cthl_2"]
                )
                series["corr_w_chi_3"] = np.divide(
                    series["crt_3"]
                    * np.sqrt(np.maximum(series["varnce_rt_3"], 0.0))
                    * series["corr_w_rt_3"]
                    - series["cthl_3"]
                    * np.sqrt(np.maximum(series["varnce_thl_3"], 0.0))
                    * series["corr_w_thl_3"],
                    series["stdev_chi_3"],
                    out=np.zeros_like(series["stdev_chi_3"]),
                    where=series["stdev_chi_3"] > 1.0e-12,
                )
                series["corr_w_chi_3"] = np.clip(
                    series["corr_w_chi_3"], -0.999999, 0.999999
                )
                used_legacy_third_correlation = True
            x_mean_3 = series[_component_field(config["x_mean"], 3)] * 1000.0
            if config["spread_is_stdev"]:
                x_var_3 = np.square(
                    series[_component_field(config["x_spread"], 3)] * 1000.0
                )
            else:
                x_var_3 = (
                    series[_component_field(config["x_spread"], 3)] * 1.0e6
                )
            component_x_means.append(x_mean_3)
            component_x_vars.append(x_var_3)
            component_y_means.append(series["w_3"])
            component_y_vars.append(series["varnce_w_3"])
            third_correlation_field = _component_field(config["correlation"], 3)
            if third_correlation_field in series:
                third_correlation = series[third_correlation_field]
            elif projection in {"w_chi", "w_rc"} and all(
                field_name in series
                for field_name in (
                    "corr_w_rt_3",
                    "corr_w_thl_3",
                    "varnce_rt_3",
                    "varnce_thl_3",
                    "crt_1",
                    "crt_2",
                    "cthl_1",
                    "cthl_2",
                )
            ):
                # Older 3G stats files did not expose corr(w,chi) for
                # component 3.  Approximate its transform coefficients from
                # the outer pair so the component is not silently omitted.
                crt_3_approx = mixture * series["crt_1"] + (1.0 - mixture) * series["crt_2"]
                cthl_3_approx = mixture * series["cthl_1"] + (1.0 - mixture) * series["cthl_2"]
                sigma_chi_3 = np.sqrt(np.maximum(x_var_3, 0.0)) / 1000.0
                third_correlation = np.divide(
                    crt_3_approx
                    * np.sqrt(np.maximum(series["varnce_rt_3"], 0.0))
                    * series["corr_w_rt_3"]
                    - cthl_3_approx
                    * np.sqrt(np.maximum(series["varnce_thl_3"], 0.0))
                    * series["corr_w_thl_3"],
                    sigma_chi_3,
                    out=np.zeros_like(sigma_chi_3),
                    where=sigma_chi_3 > 1.0e-12,
                )
                third_correlation = np.clip(third_correlation, -0.999999, 0.999999)
                used_legacy_third_correlation = True
            else:
                third_correlation = np.zeros_like(third_weight)
                used_legacy_third_correlation = True
            component_correlations.append(third_correlation)
        component_weights = _component_weights(
            mixture, third_weight if has_third_component else None
        )

        extent_x = (
            np.asarray(raw_extent[:2], dtype=float)
            if raw_extent is not None
            else np.empty(0, dtype=float)
        )
        extent_y = (
            np.asarray(raw_extent[2:], dtype=float)
            if raw_extent is not None
            else np.empty(0, dtype=float)
        )
        x_candidates = np.concatenate(
            tuple(
                candidate
                for mean, variance in zip(component_x_means, component_x_vars)
                for candidate in (
                    mean - 4.0 * np.sqrt(np.maximum(variance, 0.0)),
                    mean + 4.0 * np.sqrt(np.maximum(variance, 0.0)),
                )
            )
            + (extent_x,)
        )
        y_candidates = np.concatenate(
            tuple(
                candidate
                for mean, variance in zip(component_y_means, component_y_vars)
                for candidate in (
                    mean - 4.0 * np.sqrt(np.maximum(variance, 0.0)),
                    mean + 4.0 * np.sqrt(np.maximum(variance, 0.0)),
                )
            )
            + (extent_y,)
        )
        x_range = _finite_axis_range(x_candidates, [-1.0, 1.0])
        y_range = _finite_axis_range(y_candidates, [-1.0, 1.0])
        if projection == "w_rc":
            positive_candidates = x_candidates[np.isfinite(x_candidates) & (x_candidates >= 0.0)]
            positive_max = float(np.max(positive_candidates)) if positive_candidates.size else 0.0
            if positive_max <= 0.0:
                positive_max = 0.1
            x_range = [0.0, 1.08 * positive_max]
        x_values = np.linspace(x_range[0], x_range[1], PDF_GRID_POINTS)
        y_values = np.linspace(y_range[0], y_range[1], PDF_GRID_POINTS)
        x_grid, y_grid = np.meshgrid(x_values, y_values)
        total_density = np.zeros_like(x_grid)
        component_densities = [
            np.zeros_like(x_grid) for _component in component_weights
        ]
        upward_importance = np.zeros_like(x_grid)
        downward_importance = np.zeros_like(x_grid)
        signed_transport = np.zeros_like(x_grid)
        valid_count = 0
        for index in range(len(mixture)):
            parameters = tuple(
                value[index]
                for component_values in zip(
                    component_weights,
                    component_x_means,
                    component_x_vars,
                    component_y_means,
                    component_y_vars,
                    component_correlations,
                )
                for value in component_values
            )
            if not np.all(np.isfinite(parameters)):
                continue
            densities = [
                _bivariate_normal_density(
                    x_grid,
                    y_grid,
                    component_x_means[component_index][index],
                    component_y_means[component_index][index],
                    component_x_vars[component_index][index],
                    component_y_vars[component_index][index],
                    component_correlations[component_index][index],
                )
                for component_index in range(len(component_weights))
            ]
            for component_index, (weights, density) in enumerate(
                zip(component_weights, densities)
            ):
                weighted_density = weights[index] * density
                component_densities[component_index] += weighted_density
                total_density += weighted_density
            mean_w = (
                float(series["wm_zt"][index])
                if "wm_zt" in series
                else float(
                    sum(
                        weights[index] * means[index]
                        for weights, means in zip(
                            component_weights, component_y_means
                        )
                    )
                )
            )
            mean_rc = (
                float(series["rcm"][index] * 1000.0)
                if "rcm" in series
                else 0.0
            )
            for component_index, (weights, component_density) in enumerate(
                zip(component_weights, densities), start=1
            ):
                weight = weights[index]
                conditional_cloud = None
                if projection == "w_rt":
                    suffix = str(component_index)
                    cloud_moment_fields = (
                        f"chi_{suffix}",
                        f"stdev_chi_{suffix}",
                    )
                    if all(name in series for name in cloud_moment_fields):
                        mean_chi = float(series[f"chi_{suffix}"][index] * 1000.0)
                        var_chi = float((series[f"stdev_chi_{suffix}"][index] * 1000.0) ** 2)
                        conditioning_fields = (
                            f"corr_w_chi_{suffix}",
                            f"varnce_thl_{suffix}",
                            f"corr_rt_thl_{suffix}",
                            f"crt_{suffix}",
                            f"cthl_{suffix}",
                        )
                        if all(name in series for name in conditioning_fields):
                            covar_w_chi = float(
                                series[f"corr_w_chi_{suffix}"][index]
                                * np.sqrt(
                                    max(var_chi, 0.0)
                                    * max(
                                        component_y_vars[component_index - 1][index],
                                        0.0,
                                    )
                                )
                            )
                            var_rt_native = float(
                                series[f"varnce_rt_{suffix}"][index]
                            )
                            var_thl = float(series[f"varnce_thl_{suffix}"][index])
                            covar_rt_thl = float(
                                series[f"corr_rt_thl_{suffix}"][index]
                                * np.sqrt(max(var_rt_native, 0.0) * max(var_thl, 0.0))
                            )
                            covar_x_chi = 1.0e6 * (
                                series[f"crt_{suffix}"][index] * var_rt_native
                                - series[f"cthl_{suffix}"][index] * covar_rt_thl
                            )
                            conditional_cloud = _conditional_cloud_water(
                                x_grid,
                                y_grid,
                                mean_x=component_x_means[component_index - 1][index],
                                mean_w=component_y_means[component_index - 1][index],
                                var_x=component_x_vars[component_index - 1][index],
                                var_w=component_y_vars[component_index - 1][index],
                                corr_w_x=component_correlations[component_index - 1][index],
                                mean_chi=mean_chi,
                                var_chi=var_chi,
                                covar_x_chi=covar_x_chi,
                                covar_w_chi=covar_w_chi,
                            )
                        else:
                            conditional_cloud = np.full_like(
                                x_grid, _positive_normal_mean(mean_chi, var_chi)
                            )
                component_upward, component_downward, component_signed = _cloud_transport_fields(
                    weight * component_density,
                    x_grid,
                    y_grid,
                    projection,
                    mean_w,
                    mean_rc,
                    conditional_cloud,
                )
                upward_importance += component_upward
                downward_importance += component_downward
                signed_transport += component_signed
            valid_count += 1
        if valid_count == 0:
            raise PdfContourError("No valid PDF parameters were found in the selected window.")
        total_density /= valid_count
        component_densities = [
            density / valid_count for density in component_densities
        ]
        upward_importance /= valid_count
        downward_importance /= valid_count
        signed_transport /= valid_count
        if projection == "w_rc":
            clear_mask = x_values < 0.0
            total_density[:, clear_mask] = 0.0
            for density in component_densities:
                density[:, clear_mask] = 0.0

        components = []
        for weights, mean_x, mean_y, var_x, var_y, correlation in zip(
            component_weights,
            component_x_means,
            component_y_means,
            component_x_vars,
            component_y_vars,
            component_correlations,
        ):
            try:
                components.append(
                    aggregate_gaussian_component(
                        weights, mean_x, mean_y, var_x, var_y, correlation
                    )
                )
            except PdfContourError:
                components.append(None)

        bundle = {
            "x": x_values,
            "y": y_values,
            "total": total_density,
            "upward_importance": upward_importance,
            "downward_importance": downward_importance,
            "signed_transport": signed_transport,
        }
        bundle.update(
            {
                f"component_{component_index}": density
                for component_index, density in enumerate(
                    component_densities, start=1
                )
            }
        )
        integrals = _integrated_transport_summary(bundle)
        bundle["transport_integrals"] = integrals
        summary = [
            f"z = {selected_height:.1f} m",
            "component weights = "
            + ", ".join(
                f"ξ{component_index}={np.nanmean(weights):.3f}"
                for component_index, weights in enumerate(
                    component_weights, start=1
                )
            ),
        ]
        if used_legacy_third_correlation:
            summary.append(
                "component-3 cloud transform uses a legacy-file approximation; rerun to output corr_w_chi_3 exactly"
            )
        component_numbers = tuple(range(1, len(component_weights) + 1))
        required_cloud_fields = tuple(
            field_name
            for component_number in component_numbers
            for field_name in (
                f"rc_{component_number}",
                f"cloud_frac_{component_number}",
                f"chi_{component_number}",
                f"stdev_chi_{component_number}",
            )
        )
        if (
            projection in {"w_chi", "w_rc"}
            and all(name in series for name in ("wm_zt", "rcm"))
            and all(name in series for name in required_cloud_fields)
        ):
            contributions = []
            cloud_fractions = []
            cloudy_w_terms = []
            for component_number, weights, mean_w_component, var_w_component, correlation in zip(
                component_numbers,
                component_weights,
                component_y_means,
                component_y_vars,
                component_correlations,
            ):
                cloud_fraction_component = np.clip(
                    series[f"cloud_frac_{component_number}"], 0.0, 1.0
                )
                contribution = weights * (
                    mean_w_component - series["wm_zt"]
                ) * (
                    series[f"rc_{component_number}"] - series["rcm"]
                )
                contribution += (
                    weights
                    * correlation
                    * np.sqrt(np.maximum(var_w_component, 0.0))
                    * series[f"stdev_chi_{component_number}"]
                    * cloud_fraction_component
                )
                standard = np.divide(
                    series[f"chi_{component_number}"],
                    series[f"stdev_chi_{component_number}"],
                    out=np.zeros_like(series[f"chi_{component_number}"]),
                    where=series[f"stdev_chi_{component_number}"] > 0.0,
                )
                normal_pdf = np.exp(-0.5 * standard**2) / np.sqrt(2.0 * np.pi)
                contributions.append(contribution)
                cloud_fractions.append(weights * cloud_fraction_component)
                cloudy_w_terms.append(
                    weights
                    * (
                        mean_w_component * cloud_fraction_component
                        + correlation
                        * np.sqrt(np.maximum(var_w_component, 0.0))
                        * normal_pdf
                    )
                )
            modeled_transport = np.sum(contributions, axis=0)
            cloud_fraction = np.sum(cloud_fractions, axis=0)
            cloudy_w_numerator = np.sum(cloudy_w_terms, axis=0)
            summary.append(
                f"CLUBB w′r꜀′ = {1000.0 * np.nanmean(modeled_transport):.3e} m/s·g/kg"
            )
            if np.any(np.isfinite(cloud_fraction)):
                rcm_in_cloud = np.divide(
                    series["rcm"],
                    cloud_fraction,
                    out=np.full_like(cloud_fraction, np.nan),
                    where=cloud_fraction > 1.0e-12,
                )
                w_in_cloud = np.divide(
                    cloudy_w_numerator,
                    cloud_fraction,
                    out=np.full_like(cloud_fraction, np.nan),
                    where=cloud_fraction > 1.0e-12,
                )
                cloud_mean_transport = cloud_fraction * (
                    w_in_cloud - series["wm_zt"]
                ) * rcm_in_cloud
                within_cloud_covariance = modeled_transport - cloud_mean_transport
                mean_cloud_fraction = _finite_mean(cloud_fraction)
                mean_cloud_transport = _finite_mean(cloud_mean_transport)
                mean_within_covariance = _finite_mean(within_cloud_covariance)
                if (
                    np.isfinite(mean_cloud_transport)
                    and np.isfinite(mean_within_covariance)
                ):
                    summary.extend(
                        (
                            f"cloudy-mean transport = {1000.0 * mean_cloud_transport:.3e} m/s·g/kg",
                            f"within-cloud covariance = {1000.0 * mean_within_covariance:.3e} m/s·g/kg",
                        )
                    )
                else:
                    summary.append("cloudy transport decomposition unavailable (no resolved cloud)")
                summary.append(f"cloud fraction = {mean_cloud_fraction:.3f}")
        contour_components = (1, 2, 3) if has_third_component else (1, 2)
        all_contour_fields_available = all(
            field_name in series
            for component in contour_components
            for field_name in _trivariate_component_field_names(component)
        )
        component_times = None
        if (
            _normalize_color_signal(color_signal) != "probability"
            and all_contour_fields_available
        ):
            component_times = [
                [
                    (
                        component_weights[component_index][time_index],
                        *_component_trivariate_statistics(
                            series, time_index, component_number
                        ),
                    )
                    for component_index, component_number in enumerate(
                        contour_components
                    )
                ]
                for time_index in range(len(mixture))
            ]
        enclosed_field, enclosed_thresholds, enclosed_label = _clubb_weighted_enclosure(
            x_values, y_values, config["axes"], component_times, color_signal
        )
        if enclosed_field is None:
            enclosed_field = total_density
            enclosed_thresholds = enclosed_mass_thresholds(
                total_density, x_values, y_values
            )
        result = {
            **bundle,
            "mass_thresholds": enclosed_thresholds,
            "enclosed_mass_field": enclosed_field,
            "enclosed_mass_label": enclosed_label,
            "components": components,
            "component_means": tuple(
                (
                    float(np.nanmean(series[f"rc_{component_number}"] * 1000.0))
                    if projection == "w_rc"
                    and f"rc_{component_number}" in series
                    else float(np.nanmean(mean_x)),
                    float(np.nanmean(mean_y)),
                )
                for component_number, mean_x, mean_y in zip(
                    component_numbers, component_x_means, component_y_means
                )
            ),
            "selected_height": selected_height,
            "summary": summary,
        }
        _DENSITY_CACHE[cache_key] = result
        _DENSITY_CACHE.move_to_end(cache_key)
        while len(_DENSITY_CACHE) > _DENSITY_CACHE_MAX_ENTRIES:
            _DENSITY_CACHE.popitem(last=False)
        return result

    @staticmethod
    def _axis_domain_reference(axis, column):
        suffix = "" if column == 1 else str(column)
        return f"{axis}{suffix} domain"

    def _add_transport_trace(
        self,
        fig,
        bundle,
        column,
        panel_label,
        probability_upper,
        transport_upper,
        transport_view,
        show_color_legend,
    ):
        color_signal = _normalize_color_signal(bundle.get("color_signal"))
        signal_label = bundle.get(
            "color_signal_label", COLOR_SIGNALS["wprcp"]["label"]
        )
        color_available = bool(bundle.get("color_available", True))
        codes, normalized_probability, normalized_signed = signed_bivariate_codes(
            bundle["total"],
            bundle["signed_transport"],
            probability_upper,
            transport_upper,
        )
        customdata = np.stack(
            (
                bundle["total"],
                bundle["signed_transport"],
                normalized_probability,
                normalized_signed,
            ),
            axis=-1,
        )
        probability_name = (
            "bin probability" if bundle.get("raw_sampled") else "probability density"
        )
        probability_only = bool(bundle.get("probability_only"))
        signed_label = (
            f"signed local {signal_label} contribution"
            if color_available
            else f"{signal_label} coloring unavailable for analytic CLUBB background"
        )
        transport_hover = (
            "transport contribution=probability-only view"
            if probability_only and color_signal == "probability"
            else (
                f"{signal_label} contribution=unavailable for analytic CLUBB background"
                if probability_only
                else f"{signed_label}=%{{customdata[1]:+.4g}}"
            )
        )
        fig.add_trace(
            go.Heatmap(
                x=bundle["x"],
                y=bundle["y"],
                z=codes,
                zmin=0,
                zmax=BIVARIATE_LEVELS * (2 * BIVARIATE_LEVELS - 1) - 1,
                colorscale=signed_bivariate_reference_colorscale(),
                customdata=customdata,
                showscale=False,
                showlegend=False,
                zsmooth=False,
                name=(
                    f"{panel_label} probability + signed {signal_label}"
                    if color_available
                    else f"{panel_label} probability"
                ),
                hovertemplate=(
                    "x=%{x:.4g}<br>w=%{y:.4g}"
                    f"<br>{probability_name}=%{{customdata[0]:.4g}}"
                    f"<br>{transport_hover}"
                    f"<extra>{panel_label}</extra>"
                ),
            ),
            row=1,
            col=column,
        )

    @staticmethod
    def _add_enclosed_mass_contours(
        fig,
        *,
        x,
        y,
        field,
        thresholds,
        label,
        color,
        legendgroup,
        column,
        showlegend,
    ):
        """Add highest-weight contours and one compact shared legend entry."""
        for mass, threshold in zip(ENCLOSED_MASSES, thresholds):
            opacity = ENCLOSED_CONTOUR_OPACITY.get(float(mass), 0.7)
            width = ENCLOSED_CONTOUR_WIDTH.get(float(mass), 2.0)
            fig.add_trace(
                go.Contour(
                    x=x,
                    y=y,
                    z=field,
                    contours={
                        "start": threshold,
                        "end": threshold,
                        "size": max(abs(threshold) * 1.0e-8, 1.0e-30),
                        "coloring": "lines",
                    },
                    colorscale=[[0.0, color], [1.0, color]],
                    line={"color": color, "width": width, "dash": "solid"},
                    opacity=opacity,
                    showscale=False,
                    name=label,
                    legendgroup=legendgroup,
                    showlegend=False,
                    hoverinfo="skip",
                ),
                row=1,
                col=column,
            )
        if thresholds and showlegend:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="lines",
                    name=label,
                    legendgroup=legendgroup,
                    showlegend=True,
                    line={"color": color, "width": ENCLOSED_CONTOUR_WIDTH[0.50]},
                    hoverinfo="skip",
                ),
                row=1,
                col=column,
            )

    def _add_panel_traces(
        self,
        fig,
        bundle,
        projection,
        column,
        panel_label,
        probability_upper,
        transport_upper,
        transport_view,
        show_color_legend,
        show_clubb_legend,
        background_bundle=None,
    ):
        self._add_transport_trace(
            fig,
            background_bundle or bundle,
            column,
            "Raw SAM 3-D"
            if background_bundle is not None and background_bundle.get("raw_sampled")
            else panel_label,
            probability_upper,
            transport_upper,
            transport_view,
            show_color_legend,
        )
        self._add_enclosed_mass_contours(
            fig,
            x=bundle["x"],
            y=bundle["y"],
            field=bundle.get("enclosed_mass_field", bundle["total"]),
            thresholds=bundle["mass_thresholds"],
            label=bundle.get("enclosed_mass_label", ENCLOSED_PROBABILITY_LABEL),
            color=MASS_CONTOUR_COLOR,
            legendgroup="clubb_enclosed_mass",
            column=column,
            showlegend=show_clubb_legend,
        )
        if background_bundle is not None and background_bundle.get("raw_sampled"):
            self._add_enclosed_mass_contours(
                fig,
                x=background_bundle["x"],
                y=background_bundle["y"],
                field=background_bundle["enclosed_mass_field"],
                thresholds=background_bundle["enclosed_mass_thresholds"],
                label=background_bundle["enclosed_mass_label"],
                color=RAW_MASS_CONTOUR_COLOR,
                legendgroup="raw_sam_enclosed_mass",
                column=column,
                showlegend=show_clubb_legend,
            )
        for component_index, component in enumerate(bundle["components"]):
            if component is None:
                continue
            color = COMPONENT_COLORS[component_index]
            ellipse_x = component["ellipse_x"]
            if projection == "w_rc":
                ellipse_x = np.maximum(ellipse_x, 0.0)
            covariance = component["covariance"]
            corr_text = (
                f"{component['correlation']:.3f}"
                if np.isfinite(component["correlation"])
                else "n/a"
            )
            hover_text = (
                f"component {component_index + 1}<br>"
                f"mean time weight={component['weight']:.3f}<br>"
                f"mean x={component['mean'][0]:.4g}<br>"
                f"mean y={component['mean'][1]:.4g}<br>"
                f"σx={np.sqrt(max(covariance[0, 0], 0.0)):.4g}<br>"
                f"σy={np.sqrt(max(covariance[1, 1], 0.0)):.4g}<br>"
                f"ρ(x,y)={corr_text}"
            )
            fig.add_trace(
                go.Scatter(
                    x=ellipse_x,
                    y=component["ellipse_y"],
                    mode="lines",
                    name=f"CLUBB Gaussian {component_index + 1} shape (2× RMS)",
                    legendgroup=f"component_{component_index}",
                    showlegend=show_clubb_legend,
                    line={"color": color, "width": 3.0, "dash": "dash"},
                    text=[hover_text] * len(ellipse_x),
                    hovertemplate="%{text}<extra></extra>",
                ),
                row=1,
                col=column,
            )
            mean_x, mean_y = bundle["component_means"][component_index]
            fig.add_trace(
                go.Scatter(
                    x=[mean_x],
                    y=[mean_y],
                    mode="markers",
                    name=f"CLUBB component {component_index + 1} mean",
                    marker={"color": color, "size": 10, "symbol": "x", "line": {"width": 2}},
                    showlegend=False,
                    text=[hover_text],
                    hovertemplate="%{text}<extra></extra>",
                ),
                row=1,
                col=column,
            )

        # The binned moment is intentionally hover-only.  A persistent
        # annotation obscured the data and repeated information already shown
        # more precisely at the cursor.

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        projection = state.get("var") or "w_chi"
        transport_view = "signed"
        color_signal = _normalize_color_signal(state.get("color_signal"))
        try:
            contour_smoothing_bins = float(
                state.get(
                    "contour_smoothing_bins", DEFAULT_RAW_CONTOUR_SMOOTHING_BINS
                )
            )
        except (TypeError, ValueError):
            contour_smoothing_bins = DEFAULT_RAW_CONTOUR_SMOOTHING_BINS
        contour_smoothing_bins = min(max(contour_smoothing_bins, 0.0), 3.0)
        if projection not in PROJECTIONS:
            return shared.make_empty_figure("Select a valid PDF projection.", theme_name)
        files = list(case_data.get("files") or [])
        paths, unreadable_output = self._paths_for_projection(files, projection)
        if unreadable_output:
            return shared.make_empty_figure(
                "The selected output is being replaced by a run. PDF contours will refresh when its stats file is readable.",
                theme_name,
            )
        if not paths or (case_data.get("compare_mode") and len(paths) != len(files)):
            return shared.make_empty_figure(
                "The selected output does not contain the required PDF parameters.",
                theme_name,
            )
        height = float(state.get("height") or 0.0)
        enabled_sources = benchmark_overlay.sanitize_enabled_sources(
            case_data, global_context.get("enabled_benchmarks")
        )
        raw_snapshot = _matching_raw_snapshot(
            case_data, global_context, enabled_sources, height
        )
        raw_extent = (
            _raw_axis_extent(raw_snapshot, projection)
            if raw_snapshot is not None
            else None
        )
        try:
            bundles = [
                self._density_bundle(
                    path,
                    projection,
                    height,
                    case_data,
                    global_context,
                    raw_extent,
                    color_signal,
                )
                for path in paths
            ]
        except (PdfContourError, OSError) as exc:
            message = str(exc)
            if isinstance(exc, OSError):
                message = "The selected output changed while its PDF contour was loading. It will refresh when the run finishes."
            return shared.make_empty_figure(message, theme_name)

        source_labels = case_data.get("source_labels") or [
            f"output {index + 1}" for index in range(len(paths))
        ]
        if len(paths) == 1:
            source_labels = [source_labels[0] if case_data.get("compare_mode") else "CLUBB"]
        x_range = [
            min(float(bundle["x"][0]) for bundle in bundles),
            max(float(bundle["x"][-1]) for bundle in bundles),
        ]
        y_range = [
            min(float(bundle["y"][0]) for bundle in bundles),
            max(float(bundle["y"][-1]) for bundle in bundles),
        ]
        common_x = np.linspace(x_range[0], x_range[1], PDF_GRID_POINTS)
        common_y = np.linspace(y_range[0], y_range[1], PDF_GRID_POINTS)
        raw_bundle = (
            _raw_histogram_bundle(
                raw_snapshot,
                projection,
                common_x,
                common_y,
                color_signal,
                contour_smoothing_bins,
            )
            if raw_snapshot is not None
            else None
        )
        background_bundles = (
            [raw_bundle]
            if raw_bundle is not None
            else [_analytic_background_bundle(bundle, color_signal) for bundle in bundles]
        )
        probability_upper = robust_field_upper(
            bundle["total"] for bundle in background_bundles
        )
        transport_upper = robust_field_upper(
            np.abs(bundle["signed_transport"]) for bundle in background_bundles
        )
        transport_upper = max(transport_upper, np.finfo(float).tiny)
        fig = make_subplots(
            rows=1,
            cols=len(bundles),
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.035 if len(bundles) > 1 else 0.0,
            subplot_titles=source_labels,
        )
        for offset, (bundle, label) in enumerate(zip(bundles, source_labels)):
            column = offset + 1
            color_background = (
                raw_bundle if raw_bundle is not None else background_bundles[offset]
            )
            self._add_panel_traces(
                fig,
                bundle,
                projection,
                column,
                label,
                probability_upper,
                transport_upper,
                transport_view,
                show_color_legend=False,
                show_clubb_legend=offset == 0,
                background_bundle=color_background,
            )
        for column in range(1, len(bundles) + 1):
            if PROJECTIONS[projection]["cloud_threshold"]:
                fig.add_vline(
                    x=0.0,
                    line={"color": "#f59e0b", "width": 2.0, "dash": "dash"},
                    row=1,
                    col=column,
                )
            fig.update_xaxes(
                title_text=PROJECTIONS[projection]["x_title"],
                range=x_range,
                row=1,
                col=column,
            )
            fig.update_yaxes(
                title_text=PROJECTIONS[projection]["y_title"] if column == 1 else "",
                range=y_range,
                row=1,
                col=column,
            )
        selected_height = float(np.median([bundle["selected_height"] for bundle in bundles]))
        plot_height = shared.figure_height_for_size(global_context.get("size") or "large")
        shared.apply_figure_chrome(
            fig,
            title=(
                f"Probability + {COLOR_SIGNALS[color_signal]['label']}: {PROJECTIONS[projection]['label']} "
                f"at {selected_height:.1f} m"
            ),
            showlegend=True,
            height=plot_height,
            uirevision=shared.figure_uirevision(
                self.plot_type_id,
                case_data,
                projection,
                transport_view,
                color_signal,
                contour_smoothing_bins,
                round(selected_height, 3),
            ),
        )
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def register_callbacks(self, app):
        @app.callback(
            Output(self.height_input_id(ALL), "min"),
            Output(self.height_input_id(ALL), "max"),
            Output(self.height_input_id(ALL), "value"),
            Input("plots-global-height-range", "value"),
            State(self.height_input_id(ALL), "value"),
            prevent_initial_call=True,
        )
        def constrain_pdf_height_controls(global_height_range, local_heights):
            """Keep local contour heights inside the global plot-height selection."""
            return constrained_pdf_heights(global_height_range, local_heights)

        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Output(self.render_signal_id(MATCH), "children"),
            Input(self.var_input_id(MATCH), "value"),
            Input(self.transport_view_input_id(MATCH), "value"),
            Input(self.color_signal_input_id(MATCH), "value"),
            Input(self.contour_smoothing_input_id(MATCH), "value"),
            Input(self.height_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-global-time-range", "value"),
            Input("plots-global-time-point", "value"),
            Input("plots-time-override", "data"),
            Input("plots-selected-column", "data"),
            Input("plots-enabled-benchmarks", "data"),
            Input("theme-store", "data"),
            Input(self.size_store_id(MATCH), "data"),
            State(self.graph_id(MATCH), "relayoutData"),
        )
        def _update_pdf_contour_graph(
            projection,
            transport_view,
            color_signal,
            contour_smoothing_bins,
            height,
            case_data,
            time_range,
            time_point,
            time_override,
            selected_column,
            enabled_benchmarks,
            theme_name,
            size_store_value,
            relayout_data,
        ):
            size_value = shared.normalize_plot_size(size_store_value or "large")
            active_time = shared.resolve_active_time_values(case_data, time_range, time_point, time_override)
            signal = int(active_time["start_seconds"]) if active_time["start_seconds"] is not None else ""
            fig = self.build_figure(
                {
                    "var": projection,
                    "transport_view": transport_view,
                    "color_signal": color_signal,
                    "contour_smoothing_bins": contour_smoothing_bins,
                    "height": height,
                    "size": size_value,
                },
                {
                    "case_data": case_data,
                    "time_range": active_time["duration_minutes"],
                    "time_point": active_time["start_seconds"],
                    "selected_column": selected_column,
                    "enabled_benchmarks": enabled_benchmarks,
                    "size": size_value,
                    "theme_name": theme_name,
                },
            )
            if callback_context.triggered_id == "plots-case-data" and (case_data or {}).get("preserve_plot_view"):
                shared.apply_relayout_ranges(fig, relayout_data)
            return fig, signal


PLOT = PdfContourPlotType()

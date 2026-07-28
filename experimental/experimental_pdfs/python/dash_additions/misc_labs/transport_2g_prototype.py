"""Analytic, moment-driven trivariate transport two-Gaussian prototype.

This module is deliberately independent of Dash.  It is a small Python
laboratory for the proposed transport closure: given the moments which CLUBB
already predicts (a three-variable mean/covariance and the three marginal
third moments), it constructs two Gaussian components without a nonlinear
solve.  The construction works in standardized coordinates, preserves the
supplied mean and covariance, and uses analytic PSD caps for every requested
component contrast.

It is a prototype, not an implementation of a production CLUBB PDF type.
The result exposes caps and marginal-third-moment mismatch explicitly so that
the SAM Misc lab can distinguish an unsupported moment combination from a
good fit.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import sqrt

import numpy as np


VARIABLE_NAMES = ("w", "rₜ", "θₗ")
MIXED_THIRD_NAMES = ("WP2RTP", "WPRTP2", "WPRTPTHLP")
_PSD_MARGIN = 1.0e-8
_DIRECTION_FLOOR = 1.0e-5


@dataclass(frozen=True)
class Transport2GTuning:
    """The calibration set shared with the current PDF-10 diagnosis.

    ``moisture_tail_gain`` diagnoses the rare-component mass from ``rₜ``
    skewness.  ``center_budget`` sets how much covariance-metric room can be
    consumed by the component-center vector.  ``w_direction_scale`` limits
    center separation in ``w`` relative to the scalar direction, leaving
    vertical-velocity asymmetry to the component widths.  ``g1_wrt_capture``
    assigns positive residual ``w-rₜ`` covariance to G1; a negative grid
    covariance is instead kept out of G1 by construction.
    ``plume_structure_strength`` gates a coherent positive-w, positive-w-rₜ
    plume branch.  That branch blends placement toward a signed transport
    direction and asks G2 to retain a comparatively well-mixed covariance.
    """

    moisture_tail_gain: float = 1.0
    center_budget: float = 0.72
    w_direction_scale: float = 0.20
    g1_wrt_capture: float = 1.0
    plume_structure_strength: float = 0.50


@dataclass(frozen=True)
class TwoGaussianMixture:
    """A physical-unit trivariate two-Gaussian mixture."""

    weights: np.ndarray
    means: np.ndarray
    covariances: np.ndarray
    label: str = "Trivariate transport 2G"


@dataclass(frozen=True)
class Transport2GResult:
    """Closure result plus the reconstruction diagnostics used by the lab."""

    mixture: TwoGaussianMixture
    target_mean: np.ndarray
    target_covariance: np.ndarray
    target_third: np.ndarray
    reconstructed_mean: np.ndarray
    reconstructed_covariance: np.ndarray
    reconstructed_third: np.ndarray
    standardized_skewness: np.ndarray
    weight: float
    displacement: np.ndarray
    center_metric_fraction: float
    negative_tilt_scale: float
    contrast_scale: float
    plume_blend: float
    flags: tuple[str, ...]
    fallback_used: bool


@dataclass(frozen=True)
class MixedMomentAllocationResult:
    """Result of a bounded mixed-third-moment covariance allocation request.

    The laboratory treats raw-SAM mixed moments as extra conditional-tail
    information, not as mandatory constraints.  Each laboratory mode retains
    the grid mean and covariance.  The covariance-allocation mode fixes
    centers and preserves the existing marginal-third reconstruction exactly;
    the center-alignment mode instead preserves center-metric length and then
    rebuilds the ordinary PDF-10 width allocation.
    """

    result: Transport2GResult
    target_normalized: np.ndarray
    before_normalized: np.ndarray
    after_normalized: np.ndarray
    requested_strength: float
    psd_scale: float
    applied: bool
    message: str


def _symmetric(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T)


def _mixture_moments(mixture: TwoGaussianMixture) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the mean, covariance, and all three marginal third moments."""

    weights = np.asarray(mixture.weights, dtype=float)
    means = np.asarray(mixture.means, dtype=float)
    covariances = np.asarray(mixture.covariances, dtype=float)
    mean = np.sum(weights[:, None] * means, axis=0)
    covariance = np.zeros((3, 3), dtype=float)
    third = np.zeros(3, dtype=float)
    for weight, component_mean, component_covariance in zip(
        weights, means, covariances
    ):
        offset = component_mean - mean
        covariance += weight * (component_covariance + np.outer(offset, offset))
        third += weight * (
            offset**3 + 3.0 * offset * np.diag(component_covariance)
        )
    return mean, _symmetric(covariance), third


def mixed_third_moments(mixture: TwoGaussianMixture) -> np.ndarray:
    """Return ``<w'^2 r_t'>``, ``<w' r_t'^2>``, and ``<w'r_t'theta_l'>``.

    Each Gaussian has zero intrinsic third central moment.  Its contribution
    about the mixture mean is therefore the cubic center offset plus the
    three offset--covariance products.  This exact identity makes the Misc
    experiment a controlled test of covariance allocation rather than a
    sampled or fitted approximation.
    """

    weights = np.asarray(mixture.weights, dtype=float)
    means = np.asarray(mixture.means, dtype=float)
    covariances = np.asarray(mixture.covariances, dtype=float)
    mean = np.sum(weights[:, None] * means, axis=0)
    moments = np.zeros(3, dtype=float)
    for weight, component_mean, covariance in zip(weights, means, covariances):
        offset = component_mean - mean
        w, rt, thl = offset
        cov_w_rt = covariance[0, 1]
        cov_w_thl = covariance[0, 2]
        cov_rt_thl = covariance[1, 2]
        moments += weight * np.asarray(
            (
                w * w * rt + 2.0 * w * cov_w_rt + rt * covariance[0, 0],
                w * rt * rt + 2.0 * rt * cov_w_rt + w * covariance[1, 1],
                w * rt * thl + w * cov_rt_thl + rt * cov_w_thl + thl * cov_w_rt,
            ),
            dtype=float,
        )
    return moments


def _normalize_mixed_third(
    moments: np.ndarray,
    covariance: np.ndarray,
) -> np.ndarray:
    """Return variance-aware nondimensional mixed third moments."""

    sigma = np.sqrt(np.maximum(np.diag(np.asarray(covariance, dtype=float)), 1.0e-30))
    scale = np.asarray(
        (
            sigma[0] ** 2 * sigma[1],
            sigma[0] * sigma[1] ** 2,
            sigma[0] * sigma[1] * sigma[2],
        ),
        dtype=float,
    )
    return np.asarray(moments, dtype=float) / scale


def _mixed_moment_sensitivity(
    mixture: TwoGaussianMixture,
    covariance: np.ndarray,
) -> tuple[np.ndarray, tuple[tuple[int, int], ...], np.ndarray]:
    """Return exact normalized target sensitivity to standardized tilt contrast.

    The returned 2-by-3 matrix maps a contrast in ``(w,r_t)``,
    ``(w,theta_l)``, and ``(r_t,theta_l)`` covariance to normalized
    ``(WP2RTP, WPRTPTHLP)``.  It is valid at fixed weights and centers,
    independent of the current internal covariance allocation.
    """

    weights = np.asarray(mixture.weights, dtype=float)
    means = np.asarray(mixture.means, dtype=float)
    grid_mean = np.sum(weights[:, None] * means, axis=0)
    offsets = means - grid_mean
    difference = offsets[0] - offsets[1]
    standard_deviations = np.sqrt(
        np.maximum(np.diag(np.asarray(covariance, dtype=float)), 1.0e-30)
    )
    physical_scale = np.outer(standard_deviations, standard_deviations)
    basis_pairs = ((0, 1), (0, 2), (1, 2))
    normalized_scales = np.asarray(
        (
            standard_deviations[0] ** 2 * standard_deviations[1],
            standard_deviations[0]
            * standard_deviations[1]
            * standard_deviations[2],
        ),
        dtype=float,
    )
    raw_sensitivity = np.asarray(
        (
            (2.0 * weights[0] * difference[0], 0.0, 0.0),
            (
                weights[0] * difference[2],
                weights[0] * difference[1],
                weights[0] * difference[0],
            ),
        ),
        dtype=float,
    )
    contrast_scales = np.asarray(
        [physical_scale[left, right] for left, right in basis_pairs], dtype=float
    )
    sensitivity = (
        raw_sensitivity * contrast_scales[None, :] / normalized_scales[:, None]
    )
    return sensitivity, basis_pairs, physical_scale


def apply_mixed_moment_covariance_allocation(
    result: Transport2GResult,
    target_moments: np.ndarray,
    strength: float,
) -> MixedMomentAllocationResult:
    """Apply a PSD-capped, center-preserving mixed-moment allocation request.

    ``WP2RTP`` and ``WPRTPTHLP`` are used as the two soft targets.  At fixed
    component centers and variances, their response to a three-entry
    off-diagonal covariance contrast is linear.  A tiny least-squares solve
    therefore finds the minimum-norm contrast request.  The existing
    algebraic PSD cap then limits that request simultaneously for both
    component covariances.  ``WPRTP2`` is intentionally not targeted: it is
    shown as an independent consistency check.

    This preserves the grid mean/covariance and all three *already diagnosed*
    marginal third moments exactly up to roundoff.  (A prior PDF-10 PSD cap
    may already have softened a requested marginal third moment.)  In
    particular, it cannot make a bad component-center placement look good by
    moving a center.
    """

    target_moments = np.asarray(target_moments, dtype=float).reshape(3)
    target_normalized = _normalize_mixed_third(target_moments, result.target_covariance)
    before = _normalize_mixed_third(
        mixed_third_moments(result.mixture), result.target_covariance
    )
    requested_strength = float(np.clip(strength, 0.0, 1.0))
    unchanged = MixedMomentAllocationResult(
        result=result,
        target_normalized=target_normalized,
        before_normalized=before,
        after_normalized=before,
        requested_strength=requested_strength,
        psd_scale=0.0,
        applied=False,
        message="Mixed-moment allocation is off; the footprints use only marginal moments.",
    )
    if (
        result.fallback_used
        or requested_strength <= 0.0
        or not np.all(np.isfinite(target_normalized))
    ):
        return unchanged

    weights = np.asarray(result.mixture.weights, dtype=float)
    if np.min(weights) <= 1.0e-12:
        return replace(
            unchanged,
            message="Mixed-moment allocation was skipped because one component has negligible weight.",
        )

    first, second = np.asarray(result.mixture.covariances, dtype=float)
    ratio = float(weights[0] / weights[1])

    # The three independent off-diagonal contrast directions are expressed
    # in standardized covariance units.  At fixed centers their response is
    # analytic: a C_wr contrast changes <w'^2 r_t'> by
    # 2 a (d_w1-d_w2) dC_wr.  The r_t--theta_l mixed moment has one analogous
    # center-offset coefficient for each covariance entry.  This avoids
    # finite-difference noise in the small least-squares solve.
    sensitivity, basis_pairs, physical_scale = _mixed_moment_sensitivity(
        result.mixture, result.target_covariance
    )

    requested_delta = requested_strength * np.clip(
        target_normalized[[0, 2]] - before[[0, 2]], -6.0, 6.0
    )
    if np.max(np.abs(sensitivity)) <= 1.0e-12:
        return replace(
            unchanged,
            message="The fixed centers provide no mixed-moment sensitivity to an internal covariance allocation.",
        )
    coefficients, *_ = np.linalg.lstsq(sensitivity, requested_delta, rcond=1.0e-10)
    standardized_contrast = np.zeros((3, 3), dtype=float)
    for coefficient, (left, right) in zip(coefficients, basis_pairs):
        standardized_contrast[left, right] = coefficient
        standardized_contrast[right, left] = coefficient
    physical_contrast = standardized_contrast * physical_scale
    psd_scale = _maximum_contrast_scale(first, second, physical_contrast, ratio)
    if psd_scale <= 1.0e-12:
        return replace(
            unchanged,
            psd_scale=psd_scale,
            message="The requested mixed-moment covariance contrast was blocked by component PSD limits.",
        )

    adjusted_mixture = TwoGaussianMixture(
        weights=weights,
        means=np.asarray(result.mixture.means, dtype=float),
        covariances=np.stack(
            (
                _symmetric(first + psd_scale * physical_contrast),
                _symmetric(second - ratio * psd_scale * physical_contrast),
            )
        ),
        label=result.mixture.label,
    )
    reconstructed_mean, reconstructed_covariance, reconstructed_third = _mixture_moments(
        adjusted_mixture
    )
    adjusted = replace(
        result,
        mixture=adjusted_mixture,
        reconstructed_mean=reconstructed_mean,
        reconstructed_covariance=reconstructed_covariance,
        reconstructed_third=reconstructed_third,
        flags=result.flags
        + (
            "Raw-SAM mixed moments requested a center-preserving internal covariance allocation.",
        ),
    )
    after = _normalize_mixed_third(
        mixed_third_moments(adjusted_mixture), adjusted.target_covariance
    )
    return MixedMomentAllocationResult(
        result=adjusted,
        target_normalized=target_normalized,
        before_normalized=before,
        after_normalized=after,
        requested_strength=requested_strength,
        psd_scale=psd_scale,
        applied=True,
        message=(
            "Mixed-moment allocation changed only internal off-diagonal covariance "
            f"(PSD-safe contrast retained {psd_scale:.0%})."
        ),
    )


def _regularized_correlation(covariance: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    """Standardize an input covariance, applying a tiny PSD repair if needed."""

    covariance = _symmetric(covariance)
    if covariance.shape != (3, 3) or not np.all(np.isfinite(covariance)):
        raise ValueError("The input covariance must be a finite 3-by-3 matrix.")
    diagonal = np.diag(covariance)
    if np.any(diagonal <= 0.0):
        raise ValueError("The input covariance must have positive marginal variances.")
    standard_deviations = np.sqrt(diagonal)
    correlation = covariance / np.outer(standard_deviations, standard_deviations)
    correlation = _symmetric(correlation)
    np.fill_diagonal(correlation, 1.0)
    values, vectors = np.linalg.eigh(correlation)
    repaired = bool(np.min(values) < _PSD_MARGIN)
    if repaired:
        values = np.maximum(values, _PSD_MARGIN)
        correlation = vectors @ np.diag(values) @ vectors.T
        scale = np.sqrt(np.maximum(np.diag(correlation), _PSD_MARGIN))
        correlation = correlation / np.outer(scale, scale)
        correlation = _symmetric(correlation)
        np.fill_diagonal(correlation, 1.0)
    # A positive-definite Cholesky factor is used below, so retain a tiny
    # interior margin even for an exactly semidefinite source covariance.
    np.linalg.cholesky(correlation)
    return correlation, standard_deviations, repaired


def _maximum_contrast_scale(
    first: np.ndarray,
    second: np.ndarray,
    contrast: np.ndarray,
    second_ratio: float,
) -> float:
    """Return an analytic PSD-safe scale for a two-component contrast path.

    The path is ``C1(t)=first+t*A`` and
    ``C2(t)=second-second_ratio*t*A``.  Whiten each base covariance and use
    the limiting eigenvalue; no iteration or root solve is required.
    """

    if float(np.max(np.abs(contrast))) <= 1.0e-15:
        return 1.0

    limits: list[float] = []
    for covariance, direction in (
        (_symmetric(first), _symmetric(contrast)),
        (_symmetric(second), _symmetric(-second_ratio * contrast)),
    ):
        lower = np.linalg.cholesky(covariance)
        whitened = np.linalg.solve(lower, direction)
        whitened = np.linalg.solve(lower, whitened.T).T
        eigenvalues = np.linalg.eigvalsh(_symmetric(whitened))
        negative = eigenvalues[eigenvalues < 0.0]
        if negative.size:
            limits.append(float((1.0 - _PSD_MARGIN) / -np.min(negative)))

    limit = min(limits, default=np.inf)
    if limit >= 1.0:
        return 1.0
    return max(0.0, 0.995 * limit)


def _collapsed_fallback(
    mean: np.ndarray,
    covariance: np.ndarray,
    third: np.ndarray,
    reason: str,
) -> Transport2GResult:
    """Return the robust single-Gaussian limit represented by two copies."""

    safe_mean = np.asarray(mean, dtype=float).reshape(3)
    safe_covariance = _symmetric(np.asarray(covariance, dtype=float))
    if (
        safe_covariance.shape != (3, 3)
        or not np.all(np.isfinite(safe_covariance))
        or np.any(np.diag(safe_covariance) <= 0.0)
    ):
        safe_covariance = np.eye(3, dtype=float)
    safe_third = np.asarray(third, dtype=float).reshape(3)
    safe_third = np.where(np.isfinite(safe_third), safe_third, 0.0)
    mixture = TwoGaussianMixture(
        weights=np.array((0.5, 0.5), dtype=float),
        means=np.vstack((safe_mean, safe_mean)),
        covariances=np.stack((safe_covariance, safe_covariance)),
        label="Collapsed transport 2G fallback",
    )
    reconstructed_mean, reconstructed_covariance, reconstructed_third = _mixture_moments(
        mixture
    )
    standard_deviations = np.sqrt(np.maximum(np.diag(safe_covariance), 1.0e-30))
    return Transport2GResult(
        mixture=mixture,
        target_mean=safe_mean,
        target_covariance=safe_covariance,
        target_third=safe_third,
        reconstructed_mean=reconstructed_mean,
        reconstructed_covariance=reconstructed_covariance,
        reconstructed_third=reconstructed_third,
        standardized_skewness=safe_third / standard_deviations**3,
        weight=0.5,
        displacement=np.zeros(3, dtype=float),
        center_metric_fraction=0.0,
        negative_tilt_scale=0.0,
        contrast_scale=0.0,
        plume_blend=0.0,
        flags=(reason,),
        fallback_used=True,
    )


def diagnose_transport_2g_from_moments(
    mean: np.ndarray,
    covariance: np.ndarray,
    third: np.ndarray,
    tuning: Transport2GTuning = Transport2GTuning(),
    *,
    standardized_displacement_override: np.ndarray | None = None,
) -> Transport2GResult:
    """Diagnose a two-Gaussian transport PDF from trivariate LES moments.

    The component centers are constrained to the mean-preserving form
    ``mu1 = mean + b*d`` and ``mu2 = mean - a*d``.  The direction is inferred
    from all three marginal skewnesses, while moisture skewness sets the
    transport tail mass.  Component width differences then use the exact
    marginal-third-moment identity when it is realizable:

    ``M3_j = a*b*((b-a)*d_j**3 + 3*d_j*(V1_j - V2_j))``.

    ``standardized_displacement_override`` is an optional pre-reconstruction
    geometry request used by the vertically coherent Misc experiment.  It is
    capped to the local displacement's covariance-metric size, so it cannot
    consume more local covariance than the unmodified diagnosis.  The local
    mean, covariance, and PSD checks still run normally.

    Only linear algebra is used: first an analytic cap enforces zero G1
    ``w-rₜ`` tilt for a negative grid covariance, then a second analytic
    cap keeps the requested width/positive-tilt contrast physical.  The
    mixture mean and covariance are preserved exactly apart from roundoff.
    """

    try:
        mean = np.asarray(mean, dtype=float).reshape(3)
        covariance = _symmetric(np.asarray(covariance, dtype=float))
        third = np.asarray(third, dtype=float).reshape(3)
        if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(third)):
            raise ValueError("The input mean and third moments must be finite.")

        correlation, standard_deviations, repaired = _regularized_correlation(
            covariance
        )
        target_covariance = correlation * np.outer(
            standard_deviations, standard_deviations
        )
        flags: list[str] = []
        if repaired:
            flags.append(
                "The supplied covariance was moved a tiny distance inside the PSD cone."
            )

        raw_skewness = third / standard_deviations**3
        skewness = np.clip(raw_skewness, -12.0, 12.0)
        if not np.allclose(skewness, raw_skewness, rtol=0.0, atol=1.0e-12):
            flags.append("An extreme marginal skewness was soft-capped for this prototype.")

        tail_gain = float(np.clip(tuning.moisture_tail_gain, 0.0, 3.0))
        center_budget = float(np.clip(tuning.center_budget, 0.02, 0.97))
        w_direction_scale = float(np.clip(tuning.w_direction_scale, 0.0, 1.2))
        wrt_capture = float(np.clip(tuning.g1_wrt_capture, 0.0, 1.0))
        plume_structure = float(
            np.clip(tuning.plume_structure_strength, 0.0, 1.0)
        )

        # The maps are bounded to keep rare weights and center excursions
        # smooth across raw LES samples.  Positive versus negative moisture
        # skewness changes the direction, not which component is named G1.
        softened_skewness = np.tanh(skewness / 1.25)
        moisture_signal = abs(float(softened_skewness[1]))
        w_signal = abs(float(softened_skewness[0]))
        tail_weight_signal = np.tanh(tail_gain * abs(float(skewness[1])) / 1.25)
        weight = float(np.clip(0.5 * (1.0 - tail_weight_signal), 0.035, 0.5))
        plume_gate = min(
            max(float(softened_skewness[0]), 0.0)
            * max(float(correlation[0, 1]), 0.0)
            / 0.25**2,
            1.0,
        )
        plume_blend = plume_structure * plume_gate
        w_tail_weight = float(
            np.clip(
                0.5
                * (
                    1.0
                    - abs(float(skewness[0]))
                    / sqrt(4.0 + float(skewness[0]) ** 2)
                ),
                0.035,
                0.5,
            )
        )
        weight = (1.0 - plume_blend) * weight + plume_blend * w_tail_weight
        complement = 1.0 - weight

        raw_direction = np.array(
            (
                w_direction_scale * softened_skewness[0],
                softened_skewness[1],
                softened_skewness[2],
            ),
            dtype=float,
        )
        plume_direction = np.array(
            (1.0, correlation[0, 1], correlation[0, 2]), dtype=float
        )
        raw_direction = (
            (1.0 - plume_blend) * raw_direction
            + plume_blend * plume_direction
        )
        if float(np.linalg.norm(raw_direction)) <= _DIRECTION_FLOOR:
            direction = np.zeros(3, dtype=float)
            flags.append("The diagnosed skewness direction is nearly zero; the component centers coincide.")
        else:
            metric_length = sqrt(
                max(
                    float(raw_direction @ np.linalg.solve(correlation, raw_direction)),
                    _PSD_MARGIN,
                )
            )
            direction = raw_direction / metric_length

        separation_signal = (1.0 - plume_blend) * moisture_signal + plume_blend * max(
            moisture_signal, w_signal
        )
        local_displacement = (
            center_budget * separation_signal / sqrt(weight * complement) * direction
        )
        displacement = local_displacement
        if standardized_displacement_override is not None:
            candidate = np.asarray(
                standardized_displacement_override, dtype=float
            ).reshape(3)
            if not np.all(np.isfinite(candidate)):
                raise ValueError("The standardized displacement override must be finite.")
            local_metric = sqrt(
                max(
                    float(
                        local_displacement
                        @ np.linalg.solve(correlation, local_displacement)
                    ),
                    0.0,
                )
            )
            candidate_metric = sqrt(
                max(
                    float(candidate @ np.linalg.solve(correlation, candidate)),
                    0.0,
                )
            )
            if candidate_metric > _DIRECTION_FLOOR and local_metric > 0.0:
                displacement = candidate * min(1.0, local_metric / candidate_metric)
                if not np.allclose(displacement, local_displacement, atol=1.0e-12):
                    flags.append("The center geometry includes the bounded vertical-coherence request.")
            elif local_metric <= _DIRECTION_FLOOR:
                flags.append("The local center separation is negligible; vertical coherence left it unchanged.")
        between = weight * complement * np.outer(displacement, displacement)
        residual = _symmetric(correlation - between)
        np.linalg.cholesky(residual)
        center_metric_fraction = float(
            weight
            * complement
            * displacement
            @ np.linalg.solve(correlation, displacement)
        )

        # When grid w-r_t covariance is negative, G1 starts with exactly zero
        # internal w-r_t covariance.  The complementary tilt is carried by G2.
        # An analytic cap exposes the rare case where that strict allocation
        # itself would not be PSD.
        base_contrast = np.zeros((3, 3), dtype=float)
        negative_tilt_scale = 1.0
        if correlation[0, 1] < -1.0e-12:
            base_contrast[0, 1] = base_contrast[1, 0] = -residual[0, 1]
            negative_tilt_scale = _maximum_contrast_scale(
                residual, residual, base_contrast, weight / complement
            )
            if negative_tilt_scale < 0.999:
                flags.append(
                    "The negative residual w-rₜ allocation was limited before G1 tilt could reach zero."
                )
        first_base = _symmetric(residual + negative_tilt_scale * base_contrast)
        second_base = _symmetric(
            residual - (weight / complement) * negative_tilt_scale * base_contrast
        )
        np.linalg.cholesky(first_base)
        np.linalg.cholesky(second_base)

        # ``contrast`` has units of a normalized covariance.  Its diagonal
        # entries create the exact marginal width difference requested by
        # each supplied third moment.  Its w-r_t entry is only used for a
        # positive residual covariance when the grid covariance is positive;
        # negative grid cases deliberately retain the zero-G1-tilt base
        # allocation above.
        contrast = np.zeros((3, 3), dtype=float)
        for variable, name in enumerate(VARIABLE_NAMES):
            residual_variance = float(residual[variable, variable])
            separation = float(displacement[variable])
            if abs(separation) <= _DIRECTION_FLOOR:
                if abs(skewness[variable]) > 1.0e-5:
                    flags.append(
                        f"{name} has no resolved center separation, so its third moment is not matched."
                    )
                continue
            difference = (
                skewness[variable] / (weight * complement)
                - (complement - weight) * separation**3
            ) / (3.0 * separation)
            margin = min(1.0e-6, 0.01 * residual_variance)
            lower = (-residual_variance + margin) / complement
            upper = (residual_variance - margin) / weight
            bounded_difference = float(np.clip(difference, lower, upper))
            if abs(bounded_difference - difference) > 1.0e-10:
                flags.append(
                    f"{name} width contrast was limited before matching its third moment."
                )
            contrast[variable, variable] = complement * bounded_difference

        if (
            correlation[0, 1] > 1.0e-12
            and residual[0, 1] > 1.0e-12
            and max(wrt_capture, plume_blend) > 0.0
        ):
            contrast[0, 1] = contrast[1, 0] = (
                max(wrt_capture, plume_blend) * complement / weight * residual[0, 1]
            )

        # Mirror PDF-10's coherent-plume division of roles.  Its signed G1
        # thermodynamic structure is balanced by an equal-and-opposite
        # weighted contrast, preserving the supplied grid covariance.  G2 is
        # asked to lose w-theta_l tilt of either sign and only positive
        # r_t-theta_l tilt; an already-negative thermodynamic background is
        # deliberately retained.
        contrast[0, 2] = contrast[2, 0] = (
            plume_blend * complement / weight * second_base[0, 2]
        )
        contrast[1, 2] = contrast[2, 1] = (
            plume_blend * complement / weight * max(second_base[1, 2], 0.0)
        )

        contrast_scale = _maximum_contrast_scale(
            first_base, second_base, contrast, weight / complement
        )
        if contrast_scale < 0.999:
            flags.append(
                f"The requested width/positive-tilt contrast was reduced to {contrast_scale:.0%} by the PSD cap."
            )
        first = _symmetric(first_base + contrast_scale * contrast)
        second = _symmetric(
            second_base - (weight / complement) * contrast_scale * contrast
        )
        np.linalg.cholesky(first)
        np.linalg.cholesky(second)

        normalized_means = np.vstack(
            (complement * displacement, -weight * displacement)
        )
        physical_means = mean[None, :] + normalized_means * standard_deviations
        physical_covariances = np.stack((first, second)) * np.outer(
            standard_deviations, standard_deviations
        )
        mixture = TwoGaussianMixture(
            weights=np.array((weight, complement), dtype=float),
            means=physical_means,
            covariances=physical_covariances,
        )
        reconstructed_mean, reconstructed_covariance, reconstructed_third = _mixture_moments(
            mixture
        )
        return Transport2GResult(
            mixture=mixture,
            target_mean=mean,
            target_covariance=target_covariance,
            target_third=third,
            reconstructed_mean=reconstructed_mean,
            reconstructed_covariance=reconstructed_covariance,
            reconstructed_third=reconstructed_third,
            standardized_skewness=skewness,
            weight=weight,
            displacement=displacement * standard_deviations,
            center_metric_fraction=center_metric_fraction,
            negative_tilt_scale=negative_tilt_scale,
            contrast_scale=contrast_scale,
            plume_blend=plume_blend,
            flags=tuple(flags),
            fallback_used=False,
        )
    except (ValueError, np.linalg.LinAlgError) as error:
        return _collapsed_fallback(
            np.asarray(mean, dtype=float),
            np.asarray(covariance, dtype=float),
            np.asarray(third, dtype=float),
            f"Collapsed to the robust single-Gaussian limit: {error}",
        )


def diagnose_transport_2g_with_mixed_center_alignment(
    mean: np.ndarray,
    covariance: np.ndarray,
    third: np.ndarray,
    tuning: Transport2GTuning,
    target_moments: np.ndarray,
    strength: float,
    *,
    standardized_displacement_override: np.ndarray | None = None,
) -> MixedMomentAllocationResult:
    """Apply one bounded mixed-moment center-direction compatibility step.

    The low-level V problem is chiefly a center-placement problem, while the
    existing width/tilt contrasts commonly sit near the PSD boundary.  This
    Misc-only experiment therefore holds the center metric length fixed and
    uses one finite-difference 2-by-2 solve to rotate only its r_t and
    theta_l components toward ``WP2RTP`` and ``WPRTPTHLP`` compatibility.
    It is a single local linearization, not an iterative fit or a production
    closure proposal.
    """

    requested_strength = float(np.clip(strength, 0.0, 1.0))
    baseline = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        tuning,
        standardized_displacement_override=standardized_displacement_override,
    )
    target_moments = np.asarray(target_moments, dtype=float).reshape(3)
    target_normalized = _normalize_mixed_third(
        target_moments, baseline.target_covariance
    )
    before = _normalize_mixed_third(
        mixed_third_moments(baseline.mixture), baseline.target_covariance
    )
    unchanged = MixedMomentAllocationResult(
        result=baseline,
        target_normalized=target_normalized,
        before_normalized=before,
        after_normalized=before,
        requested_strength=requested_strength,
        psd_scale=baseline.contrast_scale,
        applied=False,
        message="Mixed-moment targets are off; the ordinary coherent-plume prior is shown.",
    )
    if (
        baseline.fallback_used
        or requested_strength <= 0.0
        or not np.all(np.isfinite(target_normalized))
    ):
        return unchanged

    correlation, _standard_deviations, _repaired = _regularized_correlation(
        covariance
    )
    reference = np.asarray(baseline.displacement, dtype=float) / np.sqrt(
        np.maximum(np.diag(baseline.target_covariance), 1.0e-30)
    )
    reference_metric = sqrt(
        max(float(reference @ np.linalg.solve(correlation, reference)), 0.0)
    )
    if reference_metric <= _DIRECTION_FLOOR:
        return replace(
            unchanged,
            message="The center separation is too small for a mixed-moment direction adjustment.",
        )

    def same_metric(candidate: np.ndarray) -> np.ndarray:
        candidate_metric = sqrt(
            max(float(candidate @ np.linalg.solve(correlation, candidate)), 0.0)
        )
        if candidate_metric <= _DIRECTION_FLOOR:
            return reference
        return candidate * (reference_metric / candidate_metric)

    def diagnose_candidate(candidate: np.ndarray) -> Transport2GResult:
        return diagnose_transport_2g_from_moments(
            mean,
            covariance,
            third,
            tuning,
            standardized_displacement_override=same_metric(candidate),
        )

    perturbation = 0.03
    sensitivity_columns: list[np.ndarray] = []
    for component in (1, 2):
        offset = np.zeros(3, dtype=float)
        offset[component] = perturbation
        positive = _normalize_mixed_third(
            mixed_third_moments(diagnose_candidate(reference + offset).mixture),
            baseline.target_covariance,
        )
        negative = _normalize_mixed_third(
            mixed_third_moments(diagnose_candidate(reference - offset).mixture),
            baseline.target_covariance,
        )
        sensitivity_columns.append(
            (positive[[0, 2]] - negative[[0, 2]]) / (2.0 * perturbation)
        )
    sensitivity = np.column_stack(sensitivity_columns)
    if np.max(np.abs(sensitivity)) <= 1.0e-12:
        return replace(
            unchanged,
            message="The local center direction has no usable mixed-moment sensitivity.",
        )
    requested_delta = requested_strength * np.clip(
        target_normalized[[0, 2]] - before[[0, 2]], -6.0, 6.0
    )
    scalar_change, *_ = np.linalg.lstsq(sensitivity, requested_delta, rcond=1.0e-10)
    # A half standardized unit is intentionally a large but bounded rotation.
    # Keeping this one-shot cap makes the laboratory responsive near PSD/cap
    # branches and makes the result a falsifiable local correction.
    scalar_change *= min(1.0, 0.50 / max(float(np.linalg.norm(scalar_change)), 1.0e-30))

    # PDF-10's ordinary reconstruction has intentional PSD/cap branches.
    # A linear response predicted on one side of such a branch can be worse
    # after a full step.  This is not an iterative fit: evaluate a fixed,
    # bounded set of fractions of the *one* local correction and retain only
    # the fraction that lowers their combined normalized residual.
    target_pair = target_normalized[[0, 2]]
    candidates: list[tuple[float, Transport2GResult, np.ndarray]] = []
    for fraction in (0.0, 0.25, 0.50, 0.75, 1.0):
        candidate = reference.copy()
        candidate[1:] += fraction * scalar_change
        candidate_result = diagnose_candidate(candidate)
        candidate_after = _normalize_mixed_third(
            mixed_third_moments(candidate_result.mixture),
            candidate_result.target_covariance,
        )
        candidates.append((fraction, candidate_result, candidate_after))
    fraction, adjusted, after = min(
        candidates,
        key=lambda item: float(np.linalg.norm(item[2][[0, 2]] - target_pair)),
    )
    changed = fraction > 0.0 and not np.allclose(
        adjusted.displacement, baseline.displacement, atol=1.0e-12
    )
    if not changed:
        adjusted = baseline
        after = before
    return MixedMomentAllocationResult(
        result=adjusted,
        target_normalized=target_normalized,
        before_normalized=before,
        after_normalized=after,
        requested_strength=requested_strength,
        psd_scale=adjusted.contrast_scale,
        applied=changed,
        message=(
            "Raw-SAM mixed targets applied one bounded r_t/theta_l center-direction rotation "
            f"({fraction:.0%} of its local step); the ordinary diagnosis retained "
            f"{adjusted.contrast_scale:.0%} of its width/tilt contrast."
            if changed
            else "No bounded mixed-moment center-direction step lowered the combined target residual."
        ),
    )


__all__ = [
    "Transport2GTuning",
    "Transport2GResult",
    "MixedMomentAllocationResult",
    "TwoGaussianMixture",
    "MIXED_THIRD_NAMES",
    "VARIABLE_NAMES",
    "apply_mixed_moment_covariance_allocation",
    "diagnose_transport_2g_from_moments",
    "diagnose_transport_2g_with_mixed_center_alignment",
    "mixed_third_moments",
]

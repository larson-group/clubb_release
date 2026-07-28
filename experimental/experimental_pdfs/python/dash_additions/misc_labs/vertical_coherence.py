"""Bounded vertical-coherence experiments for the transport 2G Misc lab.

This implements Design A from ``doc/vertically_coherent_pdf_diagnosis.md``.
The original helper couples immediate neighboring diagnosed separation
vectors.  The column helper repeats a source-to-target version of that
operation while a caller recomputes directional reach from the current G1
thermodynamic path.  Every pass remains anchored to the original local
diagnosis and reconstructs through the ordinary transport-2G moment/PSD
constraints.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp
from typing import Callable, Sequence

import numpy as np

from .transport_2g_prototype import (
    Transport2GResult,
    Transport2GTuning,
    diagnose_transport_2g_from_moments,
)


@dataclass(frozen=True)
class VerticalCoherenceSettings:
    """Bounded controls for the Design-A experiments."""

    enabled: bool = False
    max_blend: float = 0.15
    iterations: int = 1


@dataclass(frozen=True)
class VerticalCoherenceResult:
    """Adjusted local diagnosis and enough detail for an honest UI explanation."""

    result: Transport2GResult
    applied: bool
    blend: float
    lower_weight: float
    upper_weight: float
    message: str


@dataclass(frozen=True)
class IterativeVerticalCoherenceResult:
    """Whole-column fixed-point result and convergence diagnostics."""

    results: tuple[Transport2GResult, ...]
    applied: bool
    iterations_requested: int
    iterations_completed: int
    convergence_history: tuple[float, ...]
    level_blends: np.ndarray
    lower_support: np.ndarray
    upper_support: np.ndarray
    lscale_up_m: np.ndarray
    lscale_down_m: np.ndarray
    message: str


ReachProvider = Callable[
    [tuple[Transport2GResult, ...], int],
    tuple[np.ndarray, np.ndarray],
]


def standardized_displacement(result: Transport2GResult) -> np.ndarray:
    """Return the component-center separation in the local standardized units."""

    standard_deviations = np.sqrt(
        np.maximum(np.diag(result.target_covariance), 1.0e-30)
    )
    return np.asarray(result.displacement, dtype=float) / standard_deviations


def _neighbor_weight(
    local: Transport2GResult,
    neighbor: Transport2GResult | None,
    distance_m: float | None,
    reach_m: float,
) -> float:
    """Return a bounded one-neighbor Design-A support weight.

    The gate needs both diagnoses to contain a meaningful center separation.
    It smoothly suppresses an opposite direction, while still allowing a
    gradual rotation to share some support.  This is intentionally a weak
    geometrical consistency check, not a plume-tracking algorithm.
    """

    if neighbor is None or local.fallback_used or neighbor.fallback_used:
        return 0.0
    if distance_m is None or distance_m <= 0.0 or reach_m <= 0.0:
        return 0.0
    if distance_m > reach_m:
        return 0.0
    local_q = standardized_displacement(local)
    neighbor_q = standardized_displacement(neighbor)
    local_norm = float(np.linalg.norm(local_q))
    neighbor_norm = float(np.linalg.norm(neighbor_q))
    if local_norm <= 1.0e-7 or neighbor_norm <= 1.0e-7:
        return 0.0
    similarity = max(
        0.0,
        0.5
        * (
            1.0
            + float(local_q @ neighbor_q) / max(local_norm * neighbor_norm, 1.0e-30)
        ),
    )
    activity = min(1.0, local_norm / 0.20) * min(1.0, neighbor_norm / 0.20)
    distance_decay = exp(-float(distance_m / reach_m) ** 2)
    return float(np.clip(similarity * activity * distance_decay, 0.0, 1.0))


def _source_supported_result(
    mean: np.ndarray,
    covariance: np.ndarray,
    third: np.ndarray,
    tuning: Transport2GTuning,
    anchor: Transport2GResult,
    current: Sequence[Transport2GResult],
    heights_m: np.ndarray,
    target_index: int,
    lscale_up_m: np.ndarray,
    lscale_down_m: np.ndarray,
    max_blend: float,
) -> tuple[Transport2GResult, float, float, float]:
    """Blend all directly reachable source geometries into one target level."""

    target_height = float(heights_m[target_index])
    requested = np.zeros(3, dtype=float)
    lower_support = 0.0
    upper_support = 0.0
    for source_index, source in enumerate(current):
        if source_index == target_index:
            continue
        source_height = float(heights_m[source_index])
        distance = abs(target_height - source_height)
        if source_index < target_index:
            reach = float(lscale_up_m[source_index])
            weight = _neighbor_weight(anchor, source, distance, reach)
            lower_support += weight
        else:
            reach = float(lscale_down_m[source_index])
            weight = _neighbor_weight(anchor, source, distance, reach)
            upper_support += weight
        if weight > 0.0:
            requested += weight * standardized_displacement(source)

    support = lower_support + upper_support
    if support <= 1.0e-12 or max_blend <= 0.0:
        return anchor, 0.0, lower_support, upper_support
    requested /= support
    blend = max_blend * min(1.0, support)
    # Always return toward the frozen local diagnosis.  Repeated passes can
    # update the source path and support, but cannot recursively diffuse the
    # target away from the moments-only PDF.
    adjusted_geometry = (
        (1.0 - blend) * standardized_displacement(anchor) + blend * requested
    )
    result = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        tuning,
        standardized_displacement_override=adjusted_geometry,
    )
    return result, blend, lower_support, upper_support


def apply_local_vertical_coherence(
    mean: np.ndarray,
    covariance: np.ndarray,
    third: np.ndarray,
    tuning: Transport2GTuning,
    local: Transport2GResult,
    *,
    lower: Transport2GResult | None = None,
    lower_distance_m: float | None = None,
    lower_reach_m: float | None = None,
    upper: Transport2GResult | None = None,
    upper_distance_m: float | None = None,
    upper_reach_m: float | None = None,
    settings: VerticalCoherenceSettings = VerticalCoherenceSettings(),
) -> VerticalCoherenceResult:
    """Apply Design A and reconstruct a locally conservative two-Gaussian PDF."""

    if not settings.enabled:
        return VerticalCoherenceResult(
            result=local,
            applied=False,
            blend=0.0,
            lower_weight=0.0,
            upper_weight=0.0,
            message="Vertical coherence is off; this is the independent local diagnosis.",
        )

    max_blend = float(np.clip(settings.max_blend, 0.0, 0.25))
    # A lower source can influence this plane only through its upward parcel
    # reach; an upper source only through its downward reach.  These are
    # supplied by the native direct-Lscale diagnostic, never reconstructed in
    # this geometrical helper.
    lower_weight = _neighbor_weight(
        local, lower, lower_distance_m, float(lower_reach_m or 0.0)
    )
    upper_weight = _neighbor_weight(
        local, upper, upper_distance_m, float(upper_reach_m or 0.0)
    )
    support = lower_weight + upper_weight
    if support <= 1.0e-12 or max_blend <= 0.0:
        return VerticalCoherenceResult(
            result=local,
            applied=False,
            blend=0.0,
            lower_weight=lower_weight,
            upper_weight=upper_weight,
            message="No adjacent plane passed the local coherence gate; the diagnosis remains local.",
        )

    requested = np.zeros(3, dtype=float)
    if lower_weight:
        requested += lower_weight * standardized_displacement(lower)
    if upper_weight:
        requested += upper_weight * standardized_displacement(upper)
    requested /= support
    blend = max_blend * min(1.0, support)
    adjusted_geometry = (1.0 - blend) * standardized_displacement(local) + blend * requested
    result = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        tuning,
        standardized_displacement_override=adjusted_geometry,
    )
    return VerticalCoherenceResult(
        result=result,
        applied=not np.allclose(
            standardized_displacement(result), standardized_displacement(local), atol=1.0e-12
        ),
        blend=blend,
        lower_weight=lower_weight,
        upper_weight=upper_weight,
        message=(
            "Design A blends only the standardized center geometry from immediate "
            f"neighbors (blend {blend:.0%}; lower/upper support "
            f"{lower_weight:.2f}/{upper_weight:.2f}). Local moments and PSD caps "
            "are reconstructed afterward."
        ),
    )


def apply_iterative_vertical_coherence_column(
    means: np.ndarray,
    covariances: np.ndarray,
    thirds: np.ndarray,
    heights_m: np.ndarray,
    tuning: Transport2GTuning,
    local_results: Sequence[Transport2GResult],
    *,
    reach_provider: ReachProvider,
    settings: VerticalCoherenceSettings = VerticalCoherenceSettings(),
    convergence_tolerance: float = 1.0e-4,
) -> IterativeVerticalCoherenceResult:
    """Couple a column through bounded, locally anchored fixed-point passes.

    ``reach_provider`` is intentionally injected.  The Misc page supplies a
    native CLUBB direct-Lscale calculation driven by the current G1 center
    path, while unit tests can use a deterministic light-weight provider.
    Each source can influence a target only when its directional reach and
    distance-decay gate permit it.
    """

    means = np.asarray(means, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    thirds = np.asarray(thirds, dtype=float)
    heights_m = np.asarray(heights_m, dtype=float)
    anchors = tuple(local_results)
    nlevels = len(anchors)
    if (
        means.shape != (nlevels, 3)
        or covariances.shape != (nlevels, 3, 3)
        or thirds.shape != (nlevels, 3)
        or heights_m.shape != (nlevels,)
    ):
        raise ValueError("Column moments, heights, and local diagnoses must share one level dimension.")
    if nlevels and not np.all(np.diff(heights_m) > 0.0):
        raise ValueError("Vertical-coherence heights must be strictly increasing.")

    requested_iterations = int(np.clip(settings.iterations, 0, 5))
    empty = np.zeros(nlevels, dtype=float)
    if not settings.enabled or requested_iterations == 0 or nlevels == 0:
        return IterativeVerticalCoherenceResult(
            results=anchors,
            applied=False,
            iterations_requested=requested_iterations,
            iterations_completed=0,
            convergence_history=(),
            level_blends=empty,
            lower_support=empty,
            upper_support=empty,
            lscale_up_m=empty,
            lscale_down_m=empty,
            message="Iterative vertical coherence is off; every level remains local.",
        )

    max_blend = float(np.clip(settings.max_blend, 0.0, 0.25))
    current = anchors
    history: list[float] = []
    level_blends = empty.copy()
    lower_support = empty.copy()
    upper_support = empty.copy()
    lscale_up_m = empty.copy()
    lscale_down_m = empty.copy()
    for iteration in range(requested_iterations):
        lscale_up_m, lscale_down_m = (
            np.asarray(value, dtype=float)
            for value in reach_provider(current, iteration)
        )
        if lscale_up_m.shape != (nlevels,) or lscale_down_m.shape != (nlevels,):
            raise ValueError("The reach provider must return one upward/downward value per level.")
        if not np.all(np.isfinite(lscale_up_m)) or not np.all(np.isfinite(lscale_down_m)):
            raise ValueError("The reach provider returned non-finite directional lengths.")
        lscale_up_m = np.maximum(lscale_up_m, 0.0)
        lscale_down_m = np.maximum(lscale_down_m, 0.0)

        updated: list[Transport2GResult] = []
        level_blends = np.zeros(nlevels, dtype=float)
        lower_support = np.zeros(nlevels, dtype=float)
        upper_support = np.zeros(nlevels, dtype=float)
        for level in range(nlevels):
            result, blend, lower, upper = _source_supported_result(
                means[level],
                covariances[level],
                thirds[level],
                tuning,
                anchors[level],
                current,
                heights_m,
                level,
                lscale_up_m,
                lscale_down_m,
                max_blend,
            )
            updated.append(result)
            level_blends[level] = blend
            lower_support[level] = lower
            upper_support[level] = upper

        updated_tuple = tuple(updated)
        max_change = max(
            (
                float(
                    np.linalg.norm(
                        standardized_displacement(after)
                        - standardized_displacement(before)
                    )
                )
                for before, after in zip(current, updated_tuple)
            ),
            default=0.0,
        )
        history.append(max_change)
        current = updated_tuple
        if max_change <= max(float(convergence_tolerance), 0.0):
            break

    applied = any(
        not np.allclose(
            standardized_displacement(result),
            standardized_displacement(anchor),
            atol=1.0e-12,
        )
        for anchor, result in zip(anchors, current)
    )
    supported_levels = int(np.count_nonzero(level_blends > 0.0))
    final_change = history[-1] if history else 0.0
    converged = final_change <= max(float(convergence_tolerance), 0.0)
    return IterativeVerticalCoherenceResult(
        results=current,
        applied=applied,
        iterations_requested=requested_iterations,
        iterations_completed=len(history),
        convergence_history=tuple(history),
        level_blends=level_blends,
        lower_support=lower_support,
        upper_support=upper_support,
        lscale_up_m=lscale_up_m,
        lscale_down_m=lscale_down_m,
        message=(
            f"Completed {len(history)}/{requested_iterations} bounded PDF↔Lscale "
            f"passes across {supported_levels}/{nlevels} supported levels; final "
            f"standardized-geometry change {final_change:.2e}"
            + (" (converged)." if converged else ".")
        ),
    )


__all__ = [
    "IterativeVerticalCoherenceResult",
    "VerticalCoherenceResult",
    "VerticalCoherenceSettings",
    "apply_iterative_vertical_coherence_column",
    "apply_local_vertical_coherence",
    "standardized_displacement",
]

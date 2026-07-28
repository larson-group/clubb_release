"""Focused contracts for the Misc-only vertical-coherence Design-A experiment."""

import numpy as np

from .transport_2g_prototype import Transport2GTuning, diagnose_transport_2g_from_moments
from .vertical_coherence import (
    VerticalCoherenceSettings,
    apply_iterative_vertical_coherence_column,
    apply_local_vertical_coherence,
    standardized_displacement,
)


def _moments(skewness=(0.40, 0.85, -0.35)):
    standard_deviations = np.array((1.2, 1.8e-3, 1.5))
    correlation = np.array(
        (
            (1.0, 0.45, -0.25),
            (0.45, 1.0, -0.30),
            (-0.25, -0.30, 1.0),
        )
    )
    covariance = correlation * np.outer(standard_deviations, standard_deviations)
    return (
        np.array((0.15, 0.011, 299.5)),
        covariance,
        np.asarray(skewness, dtype=float) * standard_deviations**3,
    )


def test_design_a_reconstructs_the_selected_plane_after_a_neighbor_geometry_blend():
    mean, covariance, third = _moments()
    _lower_mean, _lower_covariance, lower_third = _moments((0.70, 0.65, -0.20))
    _upper_mean, _upper_covariance, upper_third = _moments((0.05, 1.10, -0.70))
    local = diagnose_transport_2g_from_moments(mean, covariance, third)
    lower = diagnose_transport_2g_from_moments(mean, covariance, lower_third)
    upper = diagnose_transport_2g_from_moments(mean, covariance, upper_third)

    coherence = apply_local_vertical_coherence(
        mean,
        covariance,
        third,
        tuning=Transport2GTuning(),
        local=local,
        lower=lower,
        lower_distance_m=40.0,
        lower_reach_m=180.0,
        upper=upper,
        upper_distance_m=40.0,
        upper_reach_m=180.0,
        settings=VerticalCoherenceSettings(enabled=True, max_blend=0.20),
    )

    assert coherence.blend <= 0.20
    assert coherence.lower_weight > 0.0
    assert coherence.upper_weight > 0.0
    assert coherence.applied
    np.testing.assert_allclose(coherence.result.reconstructed_mean, mean, atol=1.0e-12)
    np.testing.assert_allclose(
        coherence.result.reconstructed_covariance,
        coherence.result.target_covariance,
        atol=1.0e-11,
    )
    assert not np.allclose(
        standardized_displacement(coherence.result),
        standardized_displacement(local),
    )
    for component in coherence.result.mixture.covariances:
        assert np.min(np.linalg.eigvalsh(component)) >= -1.0e-10


def test_design_a_disabled_returns_the_unchanged_local_diagnosis():
    mean, covariance, third = _moments()
    local = diagnose_transport_2g_from_moments(mean, covariance, third)
    coherence = apply_local_vertical_coherence(
        mean,
        covariance,
        third,
        tuning=Transport2GTuning(),
        local=local,
        settings=VerticalCoherenceSettings(enabled=False),
    )

    assert coherence.result is local
    assert not coherence.applied
    assert coherence.blend == 0.0


def test_iterative_column_recomputes_reach_and_preserves_every_local_moment():
    variants = (
        (0.70, 0.55, -0.10),
        (0.35, 0.95, 0.20),
        (0.05, 1.15, -0.60),
        (-0.20, 0.75, -0.90),
    )
    moments = [_moments(skewness) for skewness in variants]
    means = np.asarray([item[0] for item in moments])
    covariances = np.asarray([item[1] for item in moments])
    thirds = np.asarray([item[2] for item in moments])
    local = tuple(
        diagnose_transport_2g_from_moments(mean, covariance, third)
        for mean, covariance, third in moments
    )
    provider_calls = []

    def reach_provider(current, iteration):
        provider_calls.append(
            (
                iteration,
                np.asarray([result.mixture.means[0, 1] for result in current]),
            )
        )
        return np.full(4, 260.0), np.full(4, 220.0)

    coherence = apply_iterative_vertical_coherence_column(
        means,
        covariances,
        thirds,
        np.asarray((100.0, 160.0, 230.0, 310.0)),
        Transport2GTuning(),
        local,
        reach_provider=reach_provider,
        settings=VerticalCoherenceSettings(
            enabled=True,
            max_blend=0.20,
            iterations=3,
        ),
        convergence_tolerance=0.0,
    )

    assert coherence.iterations_completed == 3
    assert len(provider_calls) == 3
    assert coherence.applied
    assert np.count_nonzero(coherence.level_blends) >= 2
    assert not np.allclose(provider_calls[0][1], provider_calls[-1][1])
    for level, result in enumerate(coherence.results):
        np.testing.assert_allclose(result.reconstructed_mean, means[level], atol=1.0e-12)
        np.testing.assert_allclose(
            result.reconstructed_covariance,
            result.target_covariance,
            atol=1.0e-11,
        )
        for component in result.mixture.covariances:
            assert np.min(np.linalg.eigvalsh(component)) >= -1.0e-10


def test_iterative_column_does_not_call_reach_provider_when_disabled():
    mean, covariance, third = _moments()
    local = (diagnose_transport_2g_from_moments(mean, covariance, third),)

    def unexpected_provider(_current, _iteration):
        raise AssertionError("disabled coherence must not run a reach calculation")

    coherence = apply_iterative_vertical_coherence_column(
        mean[None, :],
        covariance[None, :, :],
        third[None, :],
        np.asarray((100.0,)),
        Transport2GTuning(),
        local,
        reach_provider=unexpected_provider,
        settings=VerticalCoherenceSettings(enabled=False, iterations=5),
    )

    assert coherence.results == local
    assert coherence.iterations_completed == 0
    assert not coherence.applied


def test_iterative_column_does_not_leak_beyond_directional_reach():
    moments = [_moments((0.20, 0.55, -0.10)), _moments((0.70, 1.10, -0.60))]
    means = np.asarray([item[0] for item in moments])
    covariances = np.asarray([item[1] for item in moments])
    thirds = np.asarray([item[2] for item in moments])
    local = tuple(
        diagnose_transport_2g_from_moments(mean, covariance, third)
        for mean, covariance, third in moments
    )

    coherence = apply_iterative_vertical_coherence_column(
        means,
        covariances,
        thirds,
        np.asarray((100.0, 300.0)),
        Transport2GTuning(),
        local,
        reach_provider=lambda _current, _iteration: (
            np.asarray((199.0, 199.0)),
            np.asarray((199.0, 199.0)),
        ),
        settings=VerticalCoherenceSettings(enabled=True, iterations=3),
    )

    assert not coherence.applied
    assert np.count_nonzero(coherence.level_blends) == 0
    assert all(result is anchor for result, anchor in zip(coherence.results, local))

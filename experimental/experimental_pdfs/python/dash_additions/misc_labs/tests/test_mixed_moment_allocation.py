"""Contracts for the Misc-only mixed-third-moment covariance experiment."""

import numpy as np

from .transport_2g_prototype import (
    Transport2GTuning,
    apply_mixed_moment_covariance_allocation,
    diagnose_transport_2g_from_moments,
    diagnose_transport_2g_with_mixed_center_alignment,
    mixed_third_moments,
)


def _moments():
    standard_deviations = np.array((1.3, 1.7e-3, 1.6))
    correlation = np.array(
        (
            (1.0, 0.46, -0.22),
            (0.46, 1.0, -0.34),
            (-0.22, -0.34, 1.0),
        )
    )
    covariance = correlation * np.outer(standard_deviations, standard_deviations)
    third = np.array((0.45, 0.90, -0.30)) * standard_deviations**3
    return np.array((0.10, 0.011, 300.0)), covariance, third


def _from_normalized(normalized, covariance):
    sigma = np.sqrt(np.diag(covariance))
    return np.asarray(normalized) * np.array(
        (
            sigma[0] ** 2 * sigma[1],
            sigma[0] * sigma[1] ** 2,
            sigma[0] * sigma[1] * sigma[2],
        )
    )


def test_mixed_allocation_preserves_local_moments_and_component_centers():
    mean, covariance, third = _moments()
    initial = diagnose_transport_2g_from_moments(mean, covariance, third)
    baseline = _from_normalized(
        np.array((0.30, -0.10, 0.22)), covariance
    )

    allocation = apply_mixed_moment_covariance_allocation(
        initial, baseline, strength=1.0
    )

    assert allocation.applied
    np.testing.assert_allclose(allocation.result.mixture.means, initial.mixture.means)
    np.testing.assert_allclose(allocation.result.mixture.weights, initial.mixture.weights)
    np.testing.assert_allclose(allocation.result.reconstructed_mean, mean, atol=1.0e-12)
    np.testing.assert_allclose(
        allocation.result.reconstructed_covariance, initial.target_covariance, atol=1.0e-11
    )
    np.testing.assert_allclose(
        allocation.result.reconstructed_third, initial.reconstructed_third, atol=1.0e-11
    )
    assert not np.allclose(
        allocation.result.mixture.covariances, initial.mixture.covariances
    )
    assert abs(allocation.after_normalized[0] - allocation.target_normalized[0]) < abs(
        allocation.before_normalized[0] - allocation.target_normalized[0]
    )
    assert abs(allocation.after_normalized[2] - allocation.target_normalized[2]) < abs(
        allocation.before_normalized[2] - allocation.target_normalized[2]
    )
    np.testing.assert_allclose(
        mixed_third_moments(allocation.result.mixture)
        / _from_normalized(np.ones(3), covariance),
        allocation.after_normalized,
    )
    for component in allocation.result.mixture.covariances:
        assert np.min(np.linalg.eigvalsh(component)) >= -1.0e-10


def test_zero_strength_is_an_exact_noop():
    mean, covariance, third = _moments()
    initial = diagnose_transport_2g_from_moments(mean, covariance, third)
    allocation = apply_mixed_moment_covariance_allocation(
        initial, np.ones(3), strength=0.0
    )

    assert allocation.result is initial
    assert not allocation.applied
    assert allocation.psd_scale == 0.0


def test_center_alignment_preserves_moments_and_improves_a_reachable_target():
    mean, covariance, third = _moments()
    initial = diagnose_transport_2g_from_moments(mean, covariance, third)
    target = mixed_third_moments(initial.mixture).copy()
    target[[0, 2]] *= 1.6

    allocation = diagnose_transport_2g_with_mixed_center_alignment(
        mean,
        covariance,
        third,
        Transport2GTuning(),
        target,
        strength=1.0,
    )

    assert allocation.applied
    np.testing.assert_allclose(
        allocation.result.reconstructed_mean, mean, atol=1.0e-12
    )
    np.testing.assert_allclose(
        allocation.result.reconstructed_covariance,
        allocation.result.target_covariance,
        atol=1.0e-11,
    )
    assert not np.allclose(allocation.result.displacement, initial.displacement)
    for index in (0, 2):
        assert abs(allocation.after_normalized[index] - allocation.target_normalized[index]) < abs(
            allocation.before_normalized[index] - allocation.target_normalized[index]
        )

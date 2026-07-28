"""Contracts for the analytic moment-driven transport two-Gaussian prototype."""

from types import SimpleNamespace

import numpy as np

from .reports.trivariate_transport_sam_lab import make_projection_figure
from .transport_2g_prototype import (
    Transport2GTuning,
    diagnose_transport_2g_from_moments,
)


def _moments(correlation_w_rt: float = 0.45):
    standard_deviations = np.array((1.2, 1.8e-3, 1.5))
    correlation = np.array(
        (
            (1.0, correlation_w_rt, -0.25),
            (correlation_w_rt, 1.0, -0.30),
            (-0.25, -0.30, 1.0),
        )
    )
    covariance = correlation * np.outer(standard_deviations, standard_deviations)
    skewness = np.array((0.35, 0.85, -0.45))
    return (
        np.array((0.15, 0.011, 299.5)),
        covariance,
        skewness * standard_deviations**3,
    )


def _snapshot_for_figures():
    rng = np.random.default_rng(19)
    mean, covariance, third = _moments()
    samples = rng.multivariate_normal(mean, covariance, size=1024)
    chi = 0.85 * (samples[:, 1] - mean[1]) - 0.00025 * (samples[:, 2] - mean[2])
    rc = np.maximum(chi, 0.0)
    return SimpleNamespace(
        samples=samples,
        chi_samples=chi,
        rc_samples=rc,
        mean=mean,
        covariance=covariance,
        wp3=third[0],
        rtp3=third[1],
        thlp3=third[2],
        chi_mean=0.0,
        chi_coeff_rt=0.85,
        chi_coeff_thl=0.00025,
        run_dir="test-sam",
        elapsed_seconds=600,
        height_m=800.0,
    )


def test_transport_prototype_preserves_supplied_mean_and_covariance_with_psd_components():
    mean, covariance, third = _moments()
    result = diagnose_transport_2g_from_moments(mean, covariance, third)

    assert not result.fallback_used
    np.testing.assert_allclose(result.reconstructed_mean, mean, atol=1.0e-12)
    np.testing.assert_allclose(
        result.reconstructed_covariance, result.target_covariance, atol=1.0e-11
    )
    assert 0.035 <= result.weight <= 0.5
    assert 0.0 <= result.center_metric_fraction < 1.0
    assert 0.0 <= result.contrast_scale <= 1.0
    for component in result.mixture.covariances:
        assert np.min(np.linalg.eigvalsh(component)) >= -1.0e-10


def test_more_moisture_skewness_diagnoses_a_rarer_transport_component():
    mean, covariance, third = _moments()
    stronger_tail = third.copy()
    stronger_tail[1] *= 2.0
    baseline = diagnose_transport_2g_from_moments(mean, covariance, third)
    stronger = diagnose_transport_2g_from_moments(mean, covariance, stronger_tail)

    assert stronger.weight < baseline.weight
    assert abs(stronger.displacement[1]) > abs(baseline.displacement[1])


def test_negative_residual_w_rt_covariance_is_kept_out_of_g1():
    mean, covariance, _third = _moments(correlation_w_rt=-0.25)
    standard_deviations = np.sqrt(np.diag(covariance))
    # A pure moisture skew keeps the diagnosed center direction horizontal in
    # r_t, so the residual w-r_t covariance remains negative.
    third = np.array((0.0, 0.75, 0.0)) * standard_deviations**3
    result = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        Transport2GTuning(w_direction_scale=0.2),
    )

    assert not result.fallback_used
    assert result.target_covariance[0, 1] < 0.0
    assert result.negative_tilt_scale > 0.999
    assert abs(result.mixture.covariances[0, 0, 1]) < 1.0e-10
    assert result.mixture.covariances[1, 0, 1] < 0.0


def test_negative_grid_w_rt_covariance_keeps_g1_untilted_despite_center_spread():
    mean, covariance, _third = _moments(correlation_w_rt=-0.10)
    standard_deviations = np.sqrt(np.diag(covariance))
    # Opposite w and r_t skewness directions make the center-spread contribution
    # negative too, so the residual w-r_t covariance is positive.  The global
    # negative-covariance rule must still keep G1 untilted.
    third = np.array((-0.85, 0.85, -0.45)) * standard_deviations**3
    result = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        Transport2GTuning(w_direction_scale=1.2),
    )

    assert not result.fallback_used
    assert result.target_covariance[0, 1] < 0.0
    assert result.negative_tilt_scale > 0.999
    assert abs(result.mixture.covariances[0, 0, 1]) < 1.0e-10
    assert result.mixture.covariances[1, 0, 1] < 0.0


def test_coherent_plume_branch_moves_geometry_and_relaxes_g2_w_thl_tilt():
    mean, covariance, third = _moments()
    baseline = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        Transport2GTuning(plume_structure_strength=0.0),
    )
    structured = diagnose_transport_2g_from_moments(
        mean,
        covariance,
        third,
        Transport2GTuning(plume_structure_strength=1.0),
    )

    assert structured.plume_blend > 0.99
    assert abs(structured.displacement[0]) > abs(baseline.displacement[0])
    assert abs(structured.mixture.covariances[1, 0, 2]) < abs(
        baseline.mixture.covariances[1, 0, 2]
    )
    np.testing.assert_allclose(structured.reconstructed_mean, mean, atol=1.0e-12)
    np.testing.assert_allclose(
        structured.reconstructed_covariance, structured.target_covariance, atol=1.0e-11
    )


def test_sam_projection_figures_overlay_two_components_on_all_three_raw_backgrounds():
    snapshot = _snapshot_for_figures()
    result = diagnose_transport_2g_from_moments(
        snapshot.mean,
        snapshot.covariance,
        np.array((snapshot.wp3, snapshot.rtp3, snapshot.thlp3)),
    )

    figures = [make_projection_figure(snapshot, result, projection) for projection in ("rt", "chi", "rc")]

    assert [figure.layout.title.text for figure in figures] == [
        "w–rₜ · raw SAM + Gaussian footprints",
        "w–χ · raw SAM + Gaussian footprints",
        "w–r_c · raw SAM + cloudy χ branches",
    ]
    assert [len(figure.data) for figure in figures] == [7, 7, 7]
    assert all(figure.data[0].type == "heatmap" for figure in figures)
    assert all(figure.layout.showlegend is False for figure in figures)

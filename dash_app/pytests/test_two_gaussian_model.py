import numpy as np

from dash_app.shared.bivariate_heatmap import bivariate_reference_rgb
from dash_app.shared.two_gaussian_model import (
    diagnose_adg1,
    diagnose_flexible,
    diagnose_new_pdf,
    grid_covariance,
    mixture_moments,
)


def _covariance():
    return grid_covariance((1.0, 1.0, 0.8), (0.5, -0.3, -0.1))[0]


def test_all_demo_families_are_realizable():
    covariance = _covariance()
    mixtures = (
        diagnose_adg1(covariance, 1.0),
        diagnose_new_pdf(covariance, (1.0, 0.6, -0.3)),
        diagnose_flexible(covariance, 0.3, (0.55, 0.3, -0.2), 0.2, 0.3),
    )
    for mixture in mixtures:
        assert np.all(mixture.weights > 0.0)
        assert np.isclose(np.sum(mixture.weights), 1.0)
        for component_covariance in mixture.covariances:
            assert np.min(np.linalg.eigvalsh(component_covariance)) >= -1.0e-8


def test_flexible_family_preserves_grid_mean_and_covariance():
    covariance = _covariance()
    mixture = diagnose_flexible(covariance, 0.31, (0.5, 0.25, -0.15), 0.4, -0.35)
    mean, diagnosed_covariance, _wp3 = mixture_moments(mixture)
    assert np.allclose(mean, 0.0, atol=1.0e-12)
    assert np.allclose(diagnosed_covariance, covariance, atol=1.0e-8)


def test_adg1_has_zero_within_component_w_scalar_covariance():
    mixture = diagnose_adg1(_covariance(), 1.0)
    assert np.allclose(mixture.covariances[:, 0, 1:], 0.0)
    assert np.allclose(mixture.covariances[:, 1:, 0], 0.0)


def test_new_pdf_uses_equal_component_correlations():
    mixture = diagnose_new_pdf(_covariance(), (1.0, 0.6, -0.3))
    correlations = []
    for covariance in mixture.covariances:
        standard_deviations = np.sqrt(np.diag(covariance))
        correlations.append(covariance / np.outer(standard_deviations, standard_deviations))
    assert np.allclose(correlations[0], correlations[1], atol=1.0e-8)


def test_adg1_and_new_pdf_preserve_supplied_wp3():
    skew_w = 1.0
    covariance = _covariance()
    target = skew_w * covariance[0, 0] ** 1.5
    for mixture in (
        diagnose_adg1(covariance, skew_w),
        diagnose_new_pdf(covariance, (skew_w, 0.6, -0.3)),
    ):
        assert np.isclose(mixture_moments(mixture)[2], target, atol=1.0e-8)


def test_bivariate_reference_colors_distinguish_probability_transport_and_overlap():
    probability = bivariate_reference_rgb(1.0, 0.0)
    transport = bivariate_reference_rgb(0.0, 1.0)
    overlap = bivariate_reference_rgb(1.0, 1.0)
    assert probability[0] > probability[2]
    assert probability[1] > probability[2]
    assert transport[2] > transport[0]
    assert min(overlap) > 175

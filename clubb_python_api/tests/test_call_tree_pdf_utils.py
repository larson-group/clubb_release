"""Test wrappers for pdf_utilities routines in the advance call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_smooth_corr_quotient_bounds():
    """smooth_corr_quotient should keep magnitudes near correlation bounds."""
    numerator = np.array([[1.0, -2.0, 0.0]], dtype=np.float64)
    denominator = np.array([[1.0e-12, 1.0e-12, 1.0]], dtype=np.float64)
    ngrdcol, nz = numerator.shape
    out = clubb_api.smooth_corr_quotient(
        ngrdcol=ngrdcol, nz=nz, numerator=numerator, denominator=denominator, denom_thresh=1.0e-8)

    assert out.shape == numerator.shape
    assert np.all(np.isfinite(out))
    # Correlations should remain bounded in magnitude.
    assert np.all(np.abs(out) <= 1.0 + 1.0e-12)


def test_calc_comp_corrs_binormal_basic():
    """calc_comp_corrs_binormal should return equal, finite component correlations."""
    shape = (1, 3)
    xpyp = np.full(shape, 0.2, dtype=np.float64)
    xm = np.zeros(shape, dtype=np.float64)
    ym = np.zeros(shape, dtype=np.float64)
    mu_x_1 = np.ones(shape, dtype=np.float64)
    mu_x_2 = -np.ones(shape, dtype=np.float64)
    mu_y_1 = np.ones(shape, dtype=np.float64)
    mu_y_2 = -np.ones(shape, dtype=np.float64)
    sigma_x_1_sqd = np.ones(shape, dtype=np.float64)
    sigma_x_2_sqd = np.ones(shape, dtype=np.float64)
    sigma_y_1_sqd = np.ones(shape, dtype=np.float64)
    sigma_y_2_sqd = np.ones(shape, dtype=np.float64)
    mixt_frac = np.full(shape, 0.5, dtype=np.float64)
    ngrdcol, nz = shape

    corr_1, corr_2 = clubb_api.calc_comp_corrs_binormal(
        nz=nz, ngrdcol=ngrdcol, xpyp=xpyp, xm=xm, ym=ym, mu_x_1=mu_x_1, mu_x_2=mu_x_2, mu_y_1=mu_y_1,
        mu_y_2=mu_y_2, sigma_x_1_sqd=sigma_x_1_sqd, sigma_x_2_sqd=sigma_x_2_sqd,
        sigma_y_1_sqd=sigma_y_1_sqd, sigma_y_2_sqd=sigma_y_2_sqd, mixt_frac=mixt_frac,
    )

    assert corr_1.shape == shape
    assert corr_2.shape == shape
    assert np.all(np.isfinite(corr_1))
    assert np.all(np.isfinite(corr_2))
    np.testing.assert_allclose(corr_1, corr_2)
    assert np.all(np.abs(corr_1) <= 1.0 + 1.0e-12)


def test_lognormal_conversion_helpers_match_formulas():
    """Scalar lognormal conversion wrappers should match the analytic formulas."""
    mu_x = 2.0
    sigma2_on_mu2 = 0.25

    mu_x_n = clubb_api.mean_l2n(mu_x, sigma2_on_mu2)
    sigma_x_n = clubb_api.stdev_l2n(sigma2_on_mu2)

    np.testing.assert_allclose(mu_x_n, np.log(mu_x / np.sqrt(1.0 + sigma2_on_mu2)))
    np.testing.assert_allclose(sigma_x_n, np.sqrt(np.log(1.0 + sigma2_on_mu2)))


def test_corr_nn2n_family_matches_expected_formulas():
    """Correlation-conversion wrappers should match the documented formulas."""
    corr_x_y_n = np.array([[0.2, -0.3]], dtype=np.float64)
    sigma_y_n = np.array([[0.5, 0.0]], dtype=np.float64)
    y_sigma2_on_mu2 = np.array([[0.1, 0.0]], dtype=np.float64)
    ngrdcol, nz = corr_x_y_n.shape

    corr_x_y = clubb_api.corr_nn2nl(
        nz=nz, ngrdcol=ngrdcol, corr_x_y_n=corr_x_y_n, sigma_y_n=sigma_y_n, y_sigma2_on_mu2=y_sigma2_on_mu2
    )
    expected_nn2nl = np.array(
        [[corr_x_y_n[0, 0] * sigma_y_n[0, 0] / np.sqrt(y_sigma2_on_mu2[0, 0]), corr_x_y_n[0, 1]]],
        dtype=np.float64,
    )
    np.testing.assert_allclose(corr_x_y, expected_nn2nl)

    corr_ll_n = np.array([[0.1, -0.4]], dtype=np.float64)
    sigma_x_n = np.array([[0.6, 0.0]], dtype=np.float64)
    sigma_y_n2 = np.array([[0.5, 0.2]], dtype=np.float64)
    x_sigma2_on_mu2 = np.array([[0.2, 0.3]], dtype=np.float64)
    y_sigma2_on_mu2_2 = np.array([[0.1, 0.4]], dtype=np.float64)

    corr_ll = clubb_api.corr_nn2ll(
        nz=nz,
        ngrdcol=ngrdcol,
        corr_x_y_n=corr_ll_n,
        sigma_x_n=sigma_x_n,
        sigma_y_n=sigma_y_n2,
        x_sigma2_on_mu2=x_sigma2_on_mu2,
        y_sigma2_on_mu2=y_sigma2_on_mu2_2,
    )
    expected_nn2ll = np.array(
        [[
            (np.exp(sigma_x_n[0, 0] * sigma_y_n2[0, 0] * corr_ll_n[0, 0]) - 1.0)
            / np.sqrt(x_sigma2_on_mu2[0, 0] * y_sigma2_on_mu2_2[0, 0]),
            corr_ll_n[0, 1],
        ]],
        dtype=np.float64,
    )
    np.testing.assert_allclose(corr_ll, expected_nn2ll)


def test_binormal_and_chi_eta_helpers_match_formulas():
    """Scalar binormal and chi/eta correlation helpers should match analytic forms."""
    xm = clubb_api.compute_mean_binormal(1.0, 3.0, 0.25)
    xp2 = clubb_api.compute_variance_binormal(xm, 1.0, 3.0, 0.5, 1.0, 0.25)
    np.testing.assert_allclose(xm, 0.25 * 1.0 + 0.75 * 3.0)
    np.testing.assert_allclose(
        xp2,
        0.25 * ((1.0 - xm) ** 2 + 0.5**2) + 0.75 * ((3.0 - xm) ** 2 + 1.0**2),
    )

    corr_chi_x = clubb_api.calc_corr_chi_x(0.8, 0.2, 0.3, 0.4, 0.5, 0.6, -0.25)
    corr_eta_x = clubb_api.calc_corr_eta_x(0.8, 0.2, 0.3, 0.4, 0.5, 0.6, -0.25)
    np.testing.assert_allclose(corr_chi_x, 0.8 * (0.3 / 0.5) * 0.6 - 0.2 * (0.4 / 0.5) * (-0.25))
    np.testing.assert_allclose(corr_eta_x, 0.8 * (0.3 / 0.5) * 0.6 + 0.2 * (0.4 / 0.5) * (-0.25))

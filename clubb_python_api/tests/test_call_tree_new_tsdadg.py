"""Test wrappers for new_tsdadg_pdf routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_calc_l_x_skx_fnc_sign_branching_matches_formula():
    """calc_l_x_skx_fnc should swap l_x_1/l_x_2 when Skx*sgn_wpxp < 0."""
    skx = 0.8
    small_l_x_1 = 0.9
    small_l_x_2 = 0.4

    factor = abs(skx) / np.sqrt(4.0 + skx**2)

    big_l_x_1_pos, big_l_x_2_pos = clubb_api.calc_l_x_skx_fnc(
        skx=skx, sgn_wpxp=1.0, small_l_x_1=small_l_x_1, small_l_x_2=small_l_x_2,
    )
    np.testing.assert_allclose(big_l_x_1_pos, small_l_x_1 * factor)
    np.testing.assert_allclose(big_l_x_2_pos, small_l_x_2 * factor)

    big_l_x_1_neg, big_l_x_2_neg = clubb_api.calc_l_x_skx_fnc(
        skx=skx, sgn_wpxp=-1.0, small_l_x_1=small_l_x_1, small_l_x_2=small_l_x_2,
    )
    np.testing.assert_allclose(big_l_x_1_neg, small_l_x_2 * factor)
    np.testing.assert_allclose(big_l_x_2_neg, small_l_x_1 * factor)


def test_calc_setter_parameters_tsdadg_matches_formula():
    """calc_setter_parameters_tsdadg should match the closed-form equations."""
    xm = 1.2
    xp2 = 0.7
    skx = 0.5
    sgn_wpxp = 1.0
    big_l_x_1 = 0.6
    big_l_x_2 = 0.3

    (
        mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, mixt_frac,
        coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    ) = clubb_api.calc_setter_parameters_tsdadg(
        xm=xm, xp2=xp2, skx=skx, sgn_wpxp=sgn_wpxp, big_l_x_1=big_l_x_1, big_l_x_2=big_l_x_2,
    )

    factor_plus = 1.0 + skx * sgn_wpxp / np.sqrt(4.0 + skx**2)
    factor_minus = 1.0 - skx * sgn_wpxp / np.sqrt(4.0 + skx**2)
    sqrt_factor_plus_ov_minus = np.sqrt(factor_plus / factor_minus)
    sqrt_factor_minus_ov_plus = np.sqrt(factor_minus / factor_plus)

    mu_x_1_nrmlized = big_l_x_1 * sqrt_factor_plus_ov_minus * sgn_wpxp
    mu_x_2_nrmlized = -big_l_x_2 * sqrt_factor_minus_ov_plus * sgn_wpxp

    expected_mu_x_1 = xm + mu_x_1_nrmlized * np.sqrt(xp2)
    expected_mu_x_2 = xm + mu_x_2_nrmlized * np.sqrt(xp2)
    expected_mixt_frac = 1.0 / (1.0 + abs(mu_x_1_nrmlized / mu_x_2_nrmlized))

    term = (
        skx / (3.0 * expected_mixt_frac * mu_x_1_nrmlized)
        - mu_x_1_nrmlized**2 / 3.0
        + mu_x_2_nrmlized**2 / 3.0
    )
    expected_coef_sigma_x_1_sqd = (
        1.0
        - expected_mixt_frac * mu_x_1_nrmlized**2
        - (1.0 - expected_mixt_frac) * mu_x_2_nrmlized**2
        + (1.0 - expected_mixt_frac) * term
    )
    expected_coef_sigma_x_2_sqd = (
        1.0
        - expected_mixt_frac * mu_x_1_nrmlized**2
        - (1.0 - expected_mixt_frac) * mu_x_2_nrmlized**2
        - expected_mixt_frac * term
    )
    expected_sigma_x_1_sqd = expected_coef_sigma_x_1_sqd * xp2
    expected_sigma_x_2_sqd = expected_coef_sigma_x_2_sqd * xp2

    np.testing.assert_allclose(mu_x_1, expected_mu_x_1)
    np.testing.assert_allclose(mu_x_2, expected_mu_x_2)
    np.testing.assert_allclose(mixt_frac, expected_mixt_frac)
    np.testing.assert_allclose(coef_sigma_x_1_sqd, expected_coef_sigma_x_1_sqd)
    np.testing.assert_allclose(coef_sigma_x_2_sqd, expected_coef_sigma_x_2_sqd)
    np.testing.assert_allclose(sigma_x_1_sqd, expected_sigma_x_1_sqd)
    np.testing.assert_allclose(sigma_x_2_sqd, expected_sigma_x_2_sqd)


def test_tsdadg_pdf_driver_column_consistency():
    """tsdadg_pdf_driver multi-column call should match per-column calls."""
    wm = np.array([[0.10, 0.20, 0.15], [-0.05, 0.00, 0.04]], dtype=np.float64)
    rtm = np.array([[0.012, 0.011, 0.010], [0.009, 0.010, 0.011]], dtype=np.float64)
    thlm = np.array([[300.0, 300.5, 301.0], [299.5, 300.0, 300.4]], dtype=np.float64)
    wp2 = np.array([[0.20, 0.22, 0.21], [0.18, 0.17, 0.19]], dtype=np.float64)
    rtp2 = np.array([[2e-6, 2.2e-6, 2.1e-6], [1.8e-6, 1.7e-6, 1.9e-6]], dtype=np.float64)
    thlp2 = np.array([[0.08, 0.09, 0.085], [0.07, 0.072, 0.074]], dtype=np.float64)
    skw = np.array([[0.2, -0.1, 0.15], [0.05, -0.2, 0.1]], dtype=np.float64)
    skrt = np.array([[0.1, -0.05, 0.08], [0.02, -0.1, 0.04]], dtype=np.float64)
    skthl = np.array([[0.12, -0.08, 0.05], [0.03, -0.06, 0.02]], dtype=np.float64)
    wprtp = np.array([[0.001, -0.0008, 0.0009], [0.0006, -0.0007, 0.0005]], dtype=np.float64)
    wpthlp = np.array([[0.02, -0.015, 0.018], [0.012, -0.010, 0.011]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    got = clubb_api.tsdadg_pdf_driver(
        nz=nz, ngrdcol=ngrdcol, wm=wm, rtm=rtm, thlm=thlm, wp2=wp2, rtp2=rtp2, thlp2=thlp2, skw=skw,
        skrt=skrt, skthl=skthl, wprtp=wprtp, wpthlp=wpthlp)

    ncols = wm.shape[0]
    for i in range(ncols):
        got_i = clubb_api.tsdadg_pdf_driver(
            nz=nz, ngrdcol=1, wm=wm[i:i+1, :], rtm=rtm[i:i+1, :], thlm=thlm[i:i+1, :],
            wp2=wp2[i:i+1, :], rtp2=rtp2[i:i+1, :], thlp2=thlp2[i:i+1, :], skw=skw[i:i+1, :],
            skrt=skrt[i:i+1, :], skthl=skthl[i:i+1, :], wprtp=wprtp[i:i+1, :], wpthlp=wpthlp[i:i+1, :],
        )
        for arr_multi, arr_single in zip(got, got_i):
            np.testing.assert_allclose(arr_multi[i:i+1, :], arr_single)

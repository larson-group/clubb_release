"""Test wrappers for LY93 PDF routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def _calc_params_formula(xm, xp2, skx, mixt_frac):
    """Closed-form LY93 parameter equations."""
    sgn_skx = np.where(skx >= 0.0, 1.0, -1.0)
    b_x = sgn_skx * np.sqrt(xp2) * (np.abs(skx) / (1.0 - mixt_frac)) ** (1.0 / 3.0)

    mu_x_1 = xm - b_x * (1.0 - mixt_frac)
    mu_x_2 = xm + b_x * mixt_frac
    sigma_x_1_sqd = xp2 - (
        b_x**2
        * (1.0 - mixt_frac)
        * (1.0 + mixt_frac + mixt_frac**2)
        / (3.0 * mixt_frac)
    )
    sigma_x_2_sqd = xp2 + b_x**2 * (1.0 - mixt_frac) ** 2 / 3.0
    return mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd


def test_calc_params_ly93_matches_formula():
    """calc_params_ly93 should match the closed-form Lewellen-Yoh equations."""
    xm = np.array([[300.0, 301.0, 302.0]], dtype=np.float64)
    xp2 = np.array([[0.8, 0.5, 1.2]], dtype=np.float64)
    skx = np.array([[0.6, -0.3, 0.1]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.7, 0.55]], dtype=np.float64)
    nz = xm.shape[1]
    ngrdcol = xm.shape[0]

    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd = clubb_api.calc_params_ly93(
        nz=nz, ngrdcol=ngrdcol, xm=xm, xp2=xp2, skx=skx, mixt_frac=mixt_frac,
    )

    expected_mu_x_1, expected_mu_x_2, expected_sigma_x_1_sqd, expected_sigma_x_2_sqd = (
        _calc_params_formula(xm, xp2, skx, mixt_frac)
    )

    np.testing.assert_allclose(mu_x_1, expected_mu_x_1)
    np.testing.assert_allclose(mu_x_2, expected_mu_x_2)
    np.testing.assert_allclose(sigma_x_1_sqd, expected_sigma_x_1_sqd)
    np.testing.assert_allclose(sigma_x_2_sqd, expected_sigma_x_2_sqd)


def test_ly93_driver_matches_formula_for_low_skewness():
    """For |Sk|max <= 0.84, LY93 driver should use fixed mixt_frac = 0.75."""
    wm = np.array([[0.2, -0.1, 0.4]], dtype=np.float64)
    rtm = np.array([[0.010, 0.011, 0.012]], dtype=np.float64)
    thlm = np.array([[299.0, 300.0, 301.0]], dtype=np.float64)
    wp2 = np.array([[0.4, 0.5, 0.6]], dtype=np.float64)
    rtp2 = np.array([[1.5e-6, 1.1e-6, 1.7e-6]], dtype=np.float64)
    thlp2 = np.array([[0.12, 0.15, 0.11]], dtype=np.float64)
    skw = np.array([[0.2, -0.3, 0.1]], dtype=np.float64)
    skrt = np.array([[0.1, 0.2, -0.2]], dtype=np.float64)
    skthl = np.array([[0.05, -0.1, 0.15]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    (
        mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2,
        sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd,
        sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac,
    ) = clubb_api.ly93_driver(
        nz=nz, ngrdcol=ngrdcol, wm=wm, rtm=rtm, thlm=thlm, wp2=wp2, rtp2=rtp2, thlp2=thlp2, skw=skw,
        skrt=skrt, skthl=skthl)

    expected_mixt_frac = np.full_like(wm, 0.75)
    np.testing.assert_allclose(mixt_frac, expected_mixt_frac)

    expected_w = _calc_params_formula(wm, wp2, skw, expected_mixt_frac)
    expected_rt = _calc_params_formula(rtm, rtp2, skrt, expected_mixt_frac)
    expected_thl = _calc_params_formula(thlm, thlp2, skthl, expected_mixt_frac)

    np.testing.assert_allclose(mu_w_1, expected_w[0])
    np.testing.assert_allclose(mu_w_2, expected_w[1])
    np.testing.assert_allclose(sigma_w_1_sqd, expected_w[2])
    np.testing.assert_allclose(sigma_w_2_sqd, expected_w[3])

    np.testing.assert_allclose(mu_rt_1, expected_rt[0])
    np.testing.assert_allclose(mu_rt_2, expected_rt[1])
    np.testing.assert_allclose(sigma_rt_1_sqd, expected_rt[2])
    np.testing.assert_allclose(sigma_rt_2_sqd, expected_rt[3])

    np.testing.assert_allclose(mu_thl_1, expected_thl[0])
    np.testing.assert_allclose(mu_thl_2, expected_thl[1])
    np.testing.assert_allclose(sigma_thl_1_sqd, expected_thl[2])
    np.testing.assert_allclose(sigma_thl_2_sqd, expected_thl[3])

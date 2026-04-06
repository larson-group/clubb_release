"""Test wrappers for new_hybrid_pdf helper routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_calculate_w_params_matches_formula_branches():
    """calculate_w_params should match the closed-form equations and branch behavior."""
    wm = np.array([[0.1, -0.2, 0.3, 0.0]], dtype=np.float64)
    wp2 = np.array([[0.5, 0.7, 0.4, 0.6]], dtype=np.float64)
    skw = np.array([[0.25, -0.15, 0.0, 0.2]], dtype=np.float64)
    f_w = np.array([[0.4, 0.3, 0.0, 0.0]], dtype=np.float64)
    zeta_w = np.array([[0.2, 0.5, 0.4, 0.1]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    (
        mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac,
        coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
    ) = clubb_api.calculate_w_params(
        nz=nz, ngrdcol=ngrdcol, wm=wm, wp2=wp2, skw=skw, f_w=f_w, zeta_w=zeta_w)

    expected_mixt_frac = np.empty_like(mixt_frac)
    positive_f_w = f_w > 0.0
    expected_mixt_frac[positive_f_w] = (
        4.0 * f_w[positive_f_w] ** 3
        + 18.0 * f_w[positive_f_w] * (zeta_w[positive_f_w] + 1.0) * (1.0 - f_w[positive_f_w]) / (zeta_w[positive_f_w] + 2.0)
        + 6.0 * f_w[positive_f_w] ** 2 * (1.0 - f_w[positive_f_w]) / (zeta_w[positive_f_w] + 2.0)
        + skw[positive_f_w] ** 2
        - skw[positive_f_w] * np.sqrt(
            4.0 * f_w[positive_f_w] ** 3
            + 12.0 * f_w[positive_f_w] ** 2 * (1.0 - f_w[positive_f_w])
            + 36.0 * f_w[positive_f_w] * (zeta_w[positive_f_w] + 1.0)
            * (1.0 - f_w[positive_f_w]) ** 2 / (zeta_w[positive_f_w] + 2.0) ** 2
            + skw[positive_f_w] ** 2
        )
    ) / (2.0 * f_w[positive_f_w] * (f_w[positive_f_w] - 3.0) ** 2 + 2.0 * skw[positive_f_w] ** 2)

    zero_f_w_zero_skw = (f_w == 0.0) & (np.abs(skw) == 0.0)
    expected_mixt_frac[zero_f_w_zero_skw] = (zeta_w[zero_f_w_zero_skw] + 1.0) / (zeta_w[zero_f_w_zero_skw] + 2.0)

    zero_f_w_nonzero_skw = (f_w == 0.0) & (np.abs(skw) > 0.0)
    expected_mixt_frac[zero_f_w_nonzero_skw] = -1.0

    expected_mu_w_1 = np.zeros_like(mu_w_1)
    expected_mu_w_2 = np.zeros_like(mu_w_2)
    expected_sigma_w_1 = np.zeros_like(sigma_w_1)
    expected_sigma_w_2 = np.zeros_like(sigma_w_2)
    expected_coef_sigma_w_1_sqd = np.zeros_like(coef_sigma_w_1_sqd)
    expected_coef_sigma_w_2_sqd = np.zeros_like(coef_sigma_w_2_sqd)

    valid = (expected_mixt_frac > 0.0) & (expected_mixt_frac < 1.0)
    expected_mu_w_1[valid] = wm[valid] + np.sqrt(f_w[valid] * ((1.0 - expected_mixt_frac[valid]) / expected_mixt_frac[valid]) * wp2[valid])
    expected_mu_w_2[valid] = wm[valid] - (expected_mixt_frac[valid] / (1.0 - expected_mixt_frac[valid])) * (expected_mu_w_1[valid] - wm[valid])
    expected_coef_sigma_w_1_sqd[valid] = ((zeta_w[valid] + 1.0) * (1.0 - f_w[valid])) / ((zeta_w[valid] + 2.0) * expected_mixt_frac[valid])
    expected_coef_sigma_w_2_sqd[valid] = (1.0 - f_w[valid]) / ((zeta_w[valid] + 2.0) * (1.0 - expected_mixt_frac[valid]))
    expected_sigma_w_1[valid] = np.sqrt(expected_coef_sigma_w_1_sqd[valid] * wp2[valid])
    expected_sigma_w_2[valid] = np.sqrt(expected_coef_sigma_w_2_sqd[valid] * wp2[valid])

    np.testing.assert_allclose(mixt_frac, expected_mixt_frac)
    np.testing.assert_allclose(mu_w_1, expected_mu_w_1)
    np.testing.assert_allclose(mu_w_2, expected_mu_w_2)
    np.testing.assert_allclose(sigma_w_1, expected_sigma_w_1)
    np.testing.assert_allclose(sigma_w_2, expected_sigma_w_2)
    np.testing.assert_allclose(coef_sigma_w_1_sqd, expected_coef_sigma_w_1_sqd)
    np.testing.assert_allclose(coef_sigma_w_2_sqd, expected_coef_sigma_w_2_sqd)


def test_calculate_responder_params_matches_formula_branches():
    """calculate_responder_params should match both the |wpxp|>0 and zero-covariance branches."""
    xm = np.array([[0.8, -0.3, 0.2, 0.1]], dtype=np.float64)
    xp2 = np.array([[0.6, 0.9, 0.4, 0.7]], dtype=np.float64)
    skx = np.array([[0.25, -0.15, 0.35, 0.0]], dtype=np.float64)
    wpxp = np.array([[0.12, -0.08, 0.0, 0.0]], dtype=np.float64)
    wp2 = np.array([[0.5, 0.7, 0.4, 0.6]], dtype=np.float64)
    f_w = np.array([[0.4, 0.3, 0.0, 0.0]], dtype=np.float64)
    mixt_frac = np.array([[0.3, 0.6, 0.45, 0.5]], dtype=np.float64)
    nz = xm.shape[1]
    ngrdcol = xm.shape[0]

    (
        mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd,
        coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    ) = clubb_api.calculate_responder_params(
        nz=nz, ngrdcol=ngrdcol, xm=xm, xp2=xp2, skx=skx, wpxp=wpxp, wp2=wp2, f_w=f_w, mixt_frac=mixt_frac)

    expected_mu_x_1 = np.array(xm, copy=True)
    expected_mu_x_2 = np.array(xm, copy=True)
    expected_sigma_x_1_sqd = np.array(xp2, copy=True)
    expected_sigma_x_2_sqd = np.array(xp2, copy=True)
    expected_coef_sigma_x_1_sqd = np.ones_like(xp2)
    expected_coef_sigma_x_2_sqd = np.ones_like(xp2)

    nonzero = np.abs(wpxp) > 0.0
    expected_mu_x_1[nonzero] = xm[nonzero] + np.sqrt((1.0 - mixt_frac[nonzero]) / mixt_frac[nonzero]) * wpxp[nonzero] / np.sqrt(f_w[nonzero] * wp2[nonzero])
    expected_mu_x_2[nonzero] = xm[nonzero] - np.sqrt(mixt_frac[nonzero] / (1.0 - mixt_frac[nonzero])) * wpxp[nonzero] / np.sqrt(f_w[nonzero] * wp2[nonzero])
    expected_coef_sigma_x_1_sqd[nonzero] = (
        1.0
        + np.sqrt((1.0 - mixt_frac[nonzero]) / mixt_frac[nonzero]) * skx[nonzero] * np.sqrt(f_w[nonzero] * wp2[nonzero] * xp2[nonzero]) / (3.0 * wpxp[nonzero])
        - ((1.0 + mixt_frac[nonzero]) / mixt_frac[nonzero]) * wpxp[nonzero] ** 2 / (3.0 * f_w[nonzero] * wp2[nonzero] * xp2[nonzero])
    )
    expected_coef_sigma_x_2_sqd[nonzero] = (
        1.0
        - np.sqrt(mixt_frac[nonzero] / (1.0 - mixt_frac[nonzero])) * skx[nonzero] * np.sqrt(f_w[nonzero] * wp2[nonzero] * xp2[nonzero]) / (3.0 * wpxp[nonzero])
        + ((mixt_frac[nonzero] - 2.0) / (1.0 - mixt_frac[nonzero])) * wpxp[nonzero] ** 2 / (3.0 * f_w[nonzero] * wp2[nonzero] * xp2[nonzero])
    )
    expected_coef_sigma_x_1_sqd[nonzero] = np.maximum(expected_coef_sigma_x_1_sqd[nonzero], 0.0)
    expected_coef_sigma_x_2_sqd[nonzero] = np.maximum(expected_coef_sigma_x_2_sqd[nonzero], 0.0)
    expected_sigma_x_1_sqd[nonzero] = expected_coef_sigma_x_1_sqd[nonzero] * xp2[nonzero]
    expected_sigma_x_2_sqd[nonzero] = expected_coef_sigma_x_2_sqd[nonzero] * xp2[nonzero]

    np.testing.assert_allclose(mu_x_1, expected_mu_x_1)
    np.testing.assert_allclose(mu_x_2, expected_mu_x_2)
    np.testing.assert_allclose(sigma_x_1_sqd, expected_sigma_x_1_sqd)
    np.testing.assert_allclose(sigma_x_2_sqd, expected_sigma_x_2_sqd)
    np.testing.assert_allclose(coef_sigma_x_1_sqd, expected_coef_sigma_x_1_sqd)
    np.testing.assert_allclose(coef_sigma_x_2_sqd, expected_coef_sigma_x_2_sqd)


def test_calculate_coef_wp4_implicit_matches_formula():
    """calculate_coef_wp4_implicit should match the closed-form coefficient expression."""
    mixt_frac = np.array([[0.2, 0.6, 0.8]], dtype=np.float64)
    f_w = np.array([[0.4, 0.3, 0.5]], dtype=np.float64)
    coef_sigma_w_1_sqd = np.array([[0.7, 0.4, 0.9]], dtype=np.float64)
    coef_sigma_w_2_sqd = np.array([[0.6, 0.8, 0.5]], dtype=np.float64)
    nz = mixt_frac.shape[1]
    ngrdcol = mixt_frac.shape[0]

    got = clubb_api.calculate_coef_wp4_implicit(
        nz=nz, ngrdcol=ngrdcol, mixt_frac=mixt_frac, f_w=f_w,
        coef_sigma_w_1_sqd=coef_sigma_w_1_sqd, coef_sigma_w_2_sqd=coef_sigma_w_2_sqd,
    )

    expected = (
        3.0 * mixt_frac * coef_sigma_w_1_sqd**2
        + 6.0 * f_w * (1.0 - mixt_frac) * coef_sigma_w_1_sqd
        + f_w**2 * (1.0 - mixt_frac) ** 2 / mixt_frac
        + 3.0 * (1.0 - mixt_frac) * coef_sigma_w_2_sqd**2
        + 6.0 * f_w * mixt_frac * coef_sigma_w_2_sqd
        + f_w**2 * mixt_frac**2 / (1.0 - mixt_frac)
    )

    np.testing.assert_allclose(got, expected)


def test_calc_coef_wp2xp_implicit_matches_formula_branches():
    """calc_coef_wp2xp_implicit should match both the F_w>0 and F_w=0 branches."""
    wp2 = np.array([[0.5, 0.7, 0.4, 0.6]], dtype=np.float64)
    mixt_frac = np.array([[0.3, 0.6, 0.45, 0.5]], dtype=np.float64)
    f_w = np.array([[0.4, 0.3, 0.0, 0.0]], dtype=np.float64)
    coef_sigma_w_1_sqd = np.array([[0.6, 0.4, 0.0, 0.0]], dtype=np.float64)
    coef_sigma_w_2_sqd = np.array([[0.5, 0.3, 0.0, 0.0]], dtype=np.float64)
    nz = wp2.shape[1]
    ngrdcol = wp2.shape[0]

    got = clubb_api.calc_coef_wp2xp_implicit(
        nz=nz, ngrdcol=ngrdcol, wp2=wp2, mixt_frac=mixt_frac, f_w=f_w,
        coef_sigma_w_1_sqd=coef_sigma_w_1_sqd, coef_sigma_w_2_sqd=coef_sigma_w_2_sqd,
    )

    expected = np.zeros_like(got)
    positive = f_w > 0.0
    expected[positive] = (
        np.sqrt(mixt_frac[positive] * (1.0 - mixt_frac[positive]))
        * (
            f_w[positive] * ((1.0 - mixt_frac[positive]) / mixt_frac[positive] - mixt_frac[positive] / (1.0 - mixt_frac[positive]))
            + coef_sigma_w_1_sqd[positive] - coef_sigma_w_2_sqd[positive]
        )
        * np.sqrt(wp2[positive] / f_w[positive])
    )

    np.testing.assert_allclose(got, expected)


def test_calc_coefs_wpxp2_semiimpl_matches_formula_branches():
    """calc_coefs_wpxp2_semiimpl should match both the active and zero branches."""
    wp2 = np.array([[0.5, 0.7, 0.0, 0.6]], dtype=np.float64)
    wpxp = np.array([[0.12, -0.08, 0.05, 0.0]], dtype=np.float64)
    mixt_frac = np.array([[0.3, 0.6, 0.45, 0.5]], dtype=np.float64)
    f_w = np.array([[0.4, 0.3, 0.2, 0.0]], dtype=np.float64)
    coef_sigma_x_1_sqd = np.array([[0.9, 0.7, 0.3, 1.0]], dtype=np.float64)
    coef_sigma_x_2_sqd = np.array([[0.4, 0.5, 0.2, 1.0]], dtype=np.float64)
    nz = wp2.shape[1]
    ngrdcol = wp2.shape[0]

    got_coef, got_term = clubb_api.calc_coefs_wpxp2_semiimpl(
        nz=nz, ngrdcol=ngrdcol, wp2=wp2, wpxp=wpxp, mixt_frac=mixt_frac, f_w=f_w,
        coef_sigma_x_1_sqd=coef_sigma_x_1_sqd, coef_sigma_x_2_sqd=coef_sigma_x_2_sqd,
    )

    expected_coef = np.zeros_like(got_coef)
    expected_term = np.zeros_like(got_term)
    active = (f_w > 0.0) & (wp2 > 0.0)
    expected_coef[active] = (
        np.sqrt(mixt_frac[active] * (1.0 - mixt_frac[active])) * np.sqrt(f_w[active] * wp2[active])
        * (coef_sigma_x_1_sqd[active] - coef_sigma_x_2_sqd[active])
    )
    expected_term[active] = (
        np.sqrt(mixt_frac[active] * (1.0 - mixt_frac[active])) * wpxp[active] ** 2 / np.sqrt(f_w[active] * wp2[active])
        * ((1.0 - mixt_frac[active]) / mixt_frac[active] - mixt_frac[active] / (1.0 - mixt_frac[active]))
    )

    np.testing.assert_allclose(got_coef, expected_coef)
    np.testing.assert_allclose(got_term, expected_term)

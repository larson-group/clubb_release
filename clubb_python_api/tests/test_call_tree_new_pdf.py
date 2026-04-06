"""Test wrappers for new_pdf routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def _cubic_solve_roots(a_coef, b_coef, c_coef, d_coef):
    """Replicate calc_roots:cubic_solve root ordering for one cubic."""
    q_coef = (3.0 * (c_coef / a_coef) - (b_coef / a_coef) ** 2) / 9.0
    r_coef = (
        9.0 * (b_coef / a_coef) * (c_coef / a_coef)
        - 27.0 * (d_coef / a_coef)
        - 2.0 * (b_coef / a_coef) ** 3
    ) / 54.0
    determinant = q_coef**3 + r_coef**2
    sqrt_det = complex(determinant) ** 0.5
    s_coef = (complex(r_coef) + sqrt_det) ** (1.0 / 3.0)
    t_coef = (complex(r_coef) - sqrt_det) ** (1.0 / 3.0)

    root1 = -(b_coef / a_coef) / 3.0 + (s_coef + t_coef)
    root2 = -(b_coef / a_coef) / 3.0 - 0.5 * (s_coef + t_coef) + 0.5j * np.sqrt(3.0) * (s_coef - t_coef)
    root3 = -(b_coef / a_coef) / 3.0 - 0.5 * (s_coef + t_coef) - 0.5j * np.sqrt(3.0) * (s_coef - t_coef)
    return np.array([root1, root2, root3], dtype=np.complex128)


def _expected_calc_limits_f_x_responder(
    mixt_frac, skx, sgn_wpxp, max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp,
):
    """Replicate the new_pdf.F90 root-selection logic for responder F_x limits."""
    min_f_x = np.zeros_like(mixt_frac)
    max_f_x = np.zeros_like(mixt_frac)

    for i in range(mixt_frac.shape[0]):
        for j in range(mixt_frac.shape[1]):
            mf = mixt_frac[i, j]
            sk = skx[i, j]
            sgn = sgn_wpxp[i, j]
            max_pos = max_skx2_pos_skx_sgn_wpxp[i, j]
            max_neg = max_skx2_neg_skx_sgn_wpxp[i, j]

            roots_1 = _cubic_solve_roots(
                -(1.0 + mf),
                0.0,
                3.0 * mf,
                np.sqrt(mf * (1.0 - mf)) * sk * sgn,
            )
            roots_2 = _cubic_solve_roots(
                -(2.0 - mf),
                0.0,
                3.0 * (1.0 - mf),
                -np.sqrt(mf * (1.0 - mf)) * sk * sgn,
            )
            roots_1_sorted = np.sort(roots_1.real)
            roots_2_sorted = np.sort(roots_2.real)

            if sk * sgn >= 0.0:
                min_sqrt_f_x = roots_2_sorted[1]
                if sk**2 > max_neg:
                    max_sqrt_f_x = min(roots_1[0].real, roots_2_sorted[2])
                else:
                    max_sqrt_f_x = min(roots_1_sorted[2], roots_2_sorted[2])
            else:
                min_sqrt_f_x = roots_1_sorted[1]
                if sk**2 > max_pos:
                    max_sqrt_f_x = min(roots_2[0].real, roots_1_sorted[2])
                else:
                    max_sqrt_f_x = min(roots_1_sorted[2], roots_2_sorted[2])

            min_sqrt_f_x = max(min_sqrt_f_x, 0.0)
            max_sqrt_f_x = min(max_sqrt_f_x, 1.0)
            min_f_x[i, j] = min_sqrt_f_x**2
            max_f_x[i, j] = max_sqrt_f_x**2

    return min_f_x, max_f_x


def test_calc_setter_var_params_matches_formula():
    """calc_setter_var_params should match the closed-form equations in new_pdf.F90."""
    xm = np.array([[0.8, -0.3, 0.2]], dtype=np.float64)
    xp2 = np.array([[0.6, 0.9, 0.4]], dtype=np.float64)
    skx = np.array([[0.25, -0.15, 0.35]], dtype=np.float64)
    sgn_wpxp = np.array([[1.0, -1.0, 1.0]], dtype=np.float64)
    f_x = np.array([[0.4, 0.55, 0.3]], dtype=np.float64)
    zeta_x = np.array([[0.2, 0.5, 0.1]], dtype=np.float64)
    nz = xm.shape[1]
    ngrdcol = xm.shape[0]

    (
        mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, mixt_frac,
        coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    ) = clubb_api.calc_setter_var_params(
        nz=nz, ngrdcol=ngrdcol, xm=xm, xp2=xp2, skx=skx, sgn_wpxp=sgn_wpxp, f_x=f_x, zeta_x=zeta_x)

    term_under_sqrt = (
        4.0 * f_x**3
        + 12.0 * f_x**2 * (1.0 - f_x)
        + 36.0 * f_x * (zeta_x + 1.0) * (1.0 - f_x) ** 2 / (zeta_x + 2.0) ** 2
        + skx**2
    )
    expected_mixt_frac = (
        4.0 * f_x**3
        + 18.0 * f_x * (zeta_x + 1.0) * (1.0 - f_x) / (zeta_x + 2.0)
        + 6.0 * f_x**2 * (1.0 - f_x) / (zeta_x + 2.0)
        + skx**2
        - skx * sgn_wpxp * np.sqrt(term_under_sqrt)
    ) / (2.0 * f_x * (f_x - 3.0) ** 2 + 2.0 * skx**2)

    expected_mu_x_1 = xm + np.sqrt(f_x * ((1.0 - expected_mixt_frac) / expected_mixt_frac) * xp2) * sgn_wpxp
    expected_mu_x_2 = xm - (expected_mixt_frac / (1.0 - expected_mixt_frac)) * (expected_mu_x_1 - xm)

    expected_coef_sigma_x_1_sqd = ((zeta_x + 1.0) * (1.0 - f_x)) / ((zeta_x + 2.0) * expected_mixt_frac)
    expected_sigma_x_1 = np.sqrt(expected_coef_sigma_x_1_sqd * xp2)

    expected_coef_sigma_x_2_sqd = (1.0 - f_x) / ((zeta_x + 2.0) * (1.0 - expected_mixt_frac))
    expected_sigma_x_2 = np.sqrt(expected_coef_sigma_x_2_sqd * xp2)

    np.testing.assert_allclose(mixt_frac, expected_mixt_frac)
    np.testing.assert_allclose(mu_x_1, expected_mu_x_1)
    np.testing.assert_allclose(mu_x_2, expected_mu_x_2)
    np.testing.assert_allclose(coef_sigma_x_1_sqd, expected_coef_sigma_x_1_sqd)
    np.testing.assert_allclose(coef_sigma_x_2_sqd, expected_coef_sigma_x_2_sqd)
    np.testing.assert_allclose(sigma_x_1, expected_sigma_x_1)
    np.testing.assert_allclose(sigma_x_2, expected_sigma_x_2)


def test_calc_responder_params_matches_formula_branches():
    """calc_responder_params should match both the F_x>0 and F_x=0 branches."""
    xm = np.array([[0.8, -0.3, 0.2, 0.1]], dtype=np.float64)
    xp2 = np.array([[0.6, 0.9, 0.4, 0.7]], dtype=np.float64)
    skx = np.array([[0.25, -0.15, 0.35, 0.0]], dtype=np.float64)
    sgn_wpxp = np.array([[1.0, -1.0, 1.0, 1.0]], dtype=np.float64)
    f_x = np.array([[0.4, 0.55, 0.3, 0.0]], dtype=np.float64)
    mixt_frac = np.array([[0.3, 0.6, 0.45, 0.5]], dtype=np.float64)
    nz = xm.shape[1]
    ngrdcol = xm.shape[0]

    (
        mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd,
        coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    ) = clubb_api.calc_responder_params(
        nz=nz, ngrdcol=ngrdcol, xm=xm, xp2=xp2, skx=skx, sgn_wpxp=sgn_wpxp, f_x=f_x, mixt_frac=mixt_frac)

    expected_mu_x_1 = np.array(mu_x_1, copy=True)
    expected_mu_x_2 = np.array(mu_x_2, copy=True)
    expected_sigma_x_1_sqd = np.array(sigma_x_1_sqd, copy=True)
    expected_sigma_x_2_sqd = np.array(sigma_x_2_sqd, copy=True)
    expected_coef_sigma_x_1_sqd = np.array(coef_sigma_x_1_sqd, copy=True)
    expected_coef_sigma_x_2_sqd = np.array(coef_sigma_x_2_sqd, copy=True)

    positive_mask = f_x > 0.0
    expected_mu_x_1[positive_mask] = (
        xm[positive_mask]
        + np.sqrt(f_x[positive_mask] * ((1.0 - mixt_frac[positive_mask]) / mixt_frac[positive_mask]) * xp2[positive_mask])
        * sgn_wpxp[positive_mask]
    )
    expected_mu_x_2[positive_mask] = (
        xm[positive_mask]
        - (mixt_frac[positive_mask] / (1.0 - mixt_frac[positive_mask]))
        * (expected_mu_x_1[positive_mask] - xm[positive_mask])
    )
    expected_coef_sigma_x_1_sqd[positive_mask] = (
        np.sqrt(mixt_frac[positive_mask] * (1.0 - mixt_frac[positive_mask])) * skx[positive_mask] * sgn_wpxp[positive_mask]
        / (3.0 * mixt_frac[positive_mask] * np.sqrt(f_x[positive_mask]))
        - (1.0 + mixt_frac[positive_mask]) * f_x[positive_mask] / (3.0 * mixt_frac[positive_mask])
        + 1.0
    )
    expected_sigma_x_1_sqd[positive_mask] = expected_coef_sigma_x_1_sqd[positive_mask] * xp2[positive_mask]
    expected_coef_sigma_x_2_sqd[positive_mask] = (
        ((1.0 - f_x[positive_mask]) - mixt_frac[positive_mask] * expected_coef_sigma_x_1_sqd[positive_mask])
        / (1.0 - mixt_frac[positive_mask])
    )
    expected_sigma_x_2_sqd[positive_mask] = expected_coef_sigma_x_2_sqd[positive_mask] * xp2[positive_mask]

    zero_mask = ~positive_mask
    expected_mu_x_1[zero_mask] = xm[zero_mask]
    expected_mu_x_2[zero_mask] = xm[zero_mask]
    expected_sigma_x_1_sqd[zero_mask] = xp2[zero_mask]
    expected_sigma_x_2_sqd[zero_mask] = xp2[zero_mask]
    expected_coef_sigma_x_1_sqd[zero_mask] = 1.0
    expected_coef_sigma_x_2_sqd[zero_mask] = 1.0

    np.testing.assert_allclose(mu_x_1, expected_mu_x_1)
    np.testing.assert_allclose(mu_x_2, expected_mu_x_2)
    np.testing.assert_allclose(sigma_x_1_sqd, expected_sigma_x_1_sqd)
    np.testing.assert_allclose(sigma_x_2_sqd, expected_sigma_x_2_sqd)
    np.testing.assert_allclose(coef_sigma_x_1_sqd, expected_coef_sigma_x_1_sqd)
    np.testing.assert_allclose(coef_sigma_x_2_sqd, expected_coef_sigma_x_2_sqd)


def test_calc_limits_f_x_responder_matches_root_selection_logic():
    """calc_limits_f_x_responder should match the cubic-root selection logic in new_pdf.F90."""
    mixt_frac = np.array([[0.3, 0.6, 0.45, 0.4]], dtype=np.float64)
    skx = np.array([[0.25, -0.15, 0.35, -0.22]], dtype=np.float64)
    sgn_wpxp = np.array([[1.0, -1.0, -1.0, 1.0]], dtype=np.float64)
    max_skx2_pos_skx_sgn_wpxp = np.array([[0.20, 0.18, 0.16, 0.12]], dtype=np.float64)
    max_skx2_neg_skx_sgn_wpxp = np.array([[0.10, 0.08, 0.30, 0.25]], dtype=np.float64)
    nz = mixt_frac.shape[1]
    ngrdcol = mixt_frac.shape[0]

    got_min_f_x, got_max_f_x = clubb_api.calc_limits_f_x_responder(
        nz=nz, ngrdcol=ngrdcol, mixt_frac=mixt_frac, skx=skx, sgn_wpxp=sgn_wpxp,
        max_skx2_pos_skx_sgn_wpxp=max_skx2_pos_skx_sgn_wpxp,
        max_skx2_neg_skx_sgn_wpxp=max_skx2_neg_skx_sgn_wpxp,
    )

    expected_min_f_x, expected_max_f_x = _expected_calc_limits_f_x_responder(
        mixt_frac, skx, sgn_wpxp, max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp)

    np.testing.assert_allclose(got_min_f_x, expected_min_f_x)
    np.testing.assert_allclose(got_max_f_x, expected_max_f_x)


def test_calc_coef_wp4_implicit_matches_formula():
    """calc_coef_wp4_implicit should match the closed-form coefficient expression."""
    mixt_frac = np.array([[0.2, 0.6, 0.8]], dtype=np.float64)
    f_w = np.array([[0.4, 0.3, 0.5]], dtype=np.float64)
    coef_sigma_w_1_sqd = np.array([[0.7, 0.4, 0.9]], dtype=np.float64)
    coef_sigma_w_2_sqd = np.array([[0.6, 0.8, 0.5]], dtype=np.float64)
    nz = mixt_frac.shape[1]
    ngrdcol = mixt_frac.shape[0]

    got = clubb_api.calc_coef_wp4_implicit(
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


def test_calc_coef_wpxp2_implicit_matches_formula_branches():
    """calc_coef_wpxp2_implicit should match both where/else branch formulas."""
    wp2 = np.array([[0.5, 0.4, 0.0, 0.9]], dtype=np.float64)
    xp2 = np.array([[0.7, 0.3, 0.8, 0.6]], dtype=np.float64)
    wpxp = np.array([[0.1, -0.05, 0.2, 0.08]], dtype=np.float64)
    sgn_wpxp = np.array([[1.0, -1.0, 1.0, 1.0]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.55, 0.35, 0.7]], dtype=np.float64)
    f_w = np.array([[0.3, 0.4, 0.2, 0.5]], dtype=np.float64)
    f_x = np.array([[0.2, 0.6, 0.3, 0.4]], dtype=np.float64)
    coef_sigma_w_1_sqd = np.array([[0.8, 0.0, 0.7, 0.9]], dtype=np.float64)
    coef_sigma_w_2_sqd = np.array([[0.5, 0.0, 0.6, 0.4]], dtype=np.float64)
    coef_sigma_x_1_sqd = np.array([[0.7, 0.0, 0.5, 0.8]], dtype=np.float64)
    coef_sigma_x_2_sqd = np.array([[0.4, 0.0, 0.3, 0.6]], dtype=np.float64)
    nz = wp2.shape[1]
    ngrdcol = wp2.shape[0]

    got = clubb_api.calc_coef_wpxp2_implicit(
        nz=nz, ngrdcol=ngrdcol, wp2=wp2, xp2=xp2, wpxp=wpxp, sgn_wpxp=sgn_wpxp, mixt_frac=mixt_frac,
        f_w=f_w, f_x=f_x, coef_sigma_w_1_sqd=coef_sigma_w_1_sqd, coef_sigma_w_2_sqd=coef_sigma_w_2_sqd,
        coef_sigma_x_1_sqd=coef_sigma_x_1_sqd, coef_sigma_x_2_sqd=coef_sigma_x_2_sqd,
    )

    expected = np.zeros_like(got)
    for j in range(got.shape[1]):
        w2 = wp2[0, j]
        x2 = xp2[0, j]
        wx = wpxp[0, j]
        sgn = sgn_wpxp[0, j]
        mf = mixt_frac[0, j]
        fw = f_w[0, j]
        fx = f_x[0, j]
        csw1 = coef_sigma_w_1_sqd[0, j]
        csw2 = coef_sigma_w_2_sqd[0, j]
        csx1 = coef_sigma_x_1_sqd[0, j]
        csx2 = coef_sigma_x_2_sqd[0, j]

        if ((csw1 * csx1 > 0.0 or csw2 * csx2 > 0.0) and (w2 * x2 > 0.0)):
            coefs_factor = (
                np.sqrt(csw1 * csx1) - np.sqrt(csw2 * csx2)
            ) / (
                mf * np.sqrt(csw1 * csx1) + (1.0 - mf) * np.sqrt(csw2 * csx2)
            )
            expected[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(w2)
                * (
                    np.sqrt(fw) * fx * ((1.0 - mf) / mf - mf / (1.0 - mf))
                    + np.sqrt(fw) * (csx1 - csx2)
                    + 2.0 * np.sqrt(fx) * coefs_factor * sgn * wx / (np.sqrt(w2) * np.sqrt(x2))
                    - 2.0 * np.sqrt(fw) * fx * coefs_factor
                )
            )
        else:
            expected[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(w2) * np.sqrt(fw)
                * (
                    fx * ((1.0 - mf) / mf - mf / (1.0 - mf))
                    + (csx1 - csx2)
                )
            )

    np.testing.assert_allclose(got, expected)


def test_calc_coefs_wp2xp_semiimpl_matches_formula_branches():
    """calc_coefs_wp2xp_semiimpl should match where/else branch formulas."""
    wp2 = np.array([[0.5, 0.4, 0.3, 0.9]], dtype=np.float64)
    xp2 = np.array([[0.7, 0.3, 0.8, 0.6]], dtype=np.float64)
    sgn_wpxp = np.array([[1.0, -1.0, 1.0, 1.0]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.55, 0.35, 0.7]], dtype=np.float64)
    f_w = np.array([[0.3, 0.4, 0.2, 0.5]], dtype=np.float64)
    f_x = np.array([[0.2, 0.6, 0.3, 0.4]], dtype=np.float64)
    coef_sigma_w_1_sqd = np.array([[0.8, 0.0, 0.7, 0.9]], dtype=np.float64)
    coef_sigma_w_2_sqd = np.array([[0.5, 0.0, 0.6, 0.4]], dtype=np.float64)
    coef_sigma_x_1_sqd = np.array([[0.7, 0.0, 0.5, 0.8]], dtype=np.float64)
    coef_sigma_x_2_sqd = np.array([[0.4, 0.0, 0.3, 0.6]], dtype=np.float64)
    nz = wp2.shape[1]
    ngrdcol = wp2.shape[0]

    got_coef, got_term = clubb_api.calc_coefs_wp2xp_semiimpl(
        nz=nz, ngrdcol=ngrdcol, wp2=wp2, xp2=xp2, sgn_wpxp=sgn_wpxp, mixt_frac=mixt_frac, f_w=f_w, f_x=f_x,
        coef_sigma_w_1_sqd=coef_sigma_w_1_sqd, coef_sigma_w_2_sqd=coef_sigma_w_2_sqd,
        coef_sigma_x_1_sqd=coef_sigma_x_1_sqd, coef_sigma_x_2_sqd=coef_sigma_x_2_sqd,
    )

    expected_coef = np.zeros_like(got_coef)
    expected_term = np.zeros_like(got_term)
    for j in range(got_coef.shape[1]):
        w2 = wp2[0, j]
        x2 = xp2[0, j]
        sgn = sgn_wpxp[0, j]
        mf = mixt_frac[0, j]
        fw = f_w[0, j]
        fx = f_x[0, j]
        csw1 = coef_sigma_w_1_sqd[0, j]
        csw2 = coef_sigma_w_2_sqd[0, j]
        csx1 = coef_sigma_x_1_sqd[0, j]
        csx2 = coef_sigma_x_2_sqd[0, j]

        if (csw1 * csx1 > 0.0 or csw2 * csx2 > 0.0):
            coefs_factor = (
                np.sqrt(csw1 * csx1) - np.sqrt(csw2 * csx2)
            ) / (
                mf * np.sqrt(csw1 * csx1) + (1.0 - mf) * np.sqrt(csw2 * csx2)
            )
            expected_coef[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * 2.0 * np.sqrt(fw) * np.sqrt(w2) * coefs_factor
            )
            expected_term[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fx) * np.sqrt(x2) * w2 * sgn
                * (
                    fw * ((1.0 - mf) / mf - mf / (1.0 - mf))
                    + (csw1 - csw2)
                    - 2.0 * fw * coefs_factor
                )
            )
        else:
            expected_coef[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fw) * np.sqrt(w2)
                * ((1.0 - mf) / mf - mf / (1.0 - mf))
            )
            expected_term[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fx) * np.sqrt(x2) * w2 * sgn
                * (csw1 - csw2)
            )

    np.testing.assert_allclose(got_coef, expected_coef)
    np.testing.assert_allclose(got_term, expected_term)


def test_calc_coefs_wpxpyp_semiimpl_matches_formula_branches():
    """calc_coefs_wpxpyp_semiimpl should match where/else branch formulas."""
    wp2 = np.array([[0.5, 0.4, 0.3, 0.9]], dtype=np.float64)
    xp2 = np.array([[0.7, 0.3, 0.8, 0.6]], dtype=np.float64)
    yp2 = np.array([[0.2, 0.9, 0.5, 0.4]], dtype=np.float64)
    wpxp = np.array([[0.1, -0.05, 0.2, 0.08]], dtype=np.float64)
    wpyp = np.array([[0.03, 0.07, -0.04, 0.02]], dtype=np.float64)
    sgn_wpxp = np.array([[1.0, -1.0, 1.0, 1.0]], dtype=np.float64)
    sgn_wpyp = np.array([[1.0, 1.0, -1.0, 1.0]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.55, 0.35, 0.7]], dtype=np.float64)
    f_w = np.array([[0.3, 0.4, 0.2, 0.5]], dtype=np.float64)
    f_x = np.array([[0.2, 0.6, 0.3, 0.4]], dtype=np.float64)
    f_y = np.array([[0.25, 0.5, 0.35, 0.45]], dtype=np.float64)
    coef_sigma_w_1_sqd = np.array([[0.8, 0.0, 0.7, 0.9]], dtype=np.float64)
    coef_sigma_w_2_sqd = np.array([[0.5, 0.0, 0.6, 0.4]], dtype=np.float64)
    coef_sigma_x_1_sqd = np.array([[0.7, 0.0, 0.5, 0.8]], dtype=np.float64)
    coef_sigma_x_2_sqd = np.array([[0.4, 0.0, 0.3, 0.6]], dtype=np.float64)
    coef_sigma_y_1_sqd = np.array([[0.6, 0.0, 0.4, 0.7]], dtype=np.float64)
    coef_sigma_y_2_sqd = np.array([[0.3, 0.0, 0.2, 0.5]], dtype=np.float64)
    nz = wp2.shape[1]
    ngrdcol = wp2.shape[0]

    got_coef, got_term = clubb_api.calc_coefs_wpxpyp_semiimpl(
        nz=nz, ngrdcol=ngrdcol, wp2=wp2, xp2=xp2, yp2=yp2, wpxp=wpxp, wpyp=wpyp, sgn_wpxp=sgn_wpxp,
        sgn_wpyp=sgn_wpyp, mixt_frac=mixt_frac, f_w=f_w, f_x=f_x, f_y=f_y,
        coef_sigma_w_1_sqd=coef_sigma_w_1_sqd, coef_sigma_w_2_sqd=coef_sigma_w_2_sqd,
        coef_sigma_x_1_sqd=coef_sigma_x_1_sqd, coef_sigma_x_2_sqd=coef_sigma_x_2_sqd,
        coef_sigma_y_1_sqd=coef_sigma_y_1_sqd, coef_sigma_y_2_sqd=coef_sigma_y_2_sqd,
    )

    expected_coef = np.zeros_like(got_coef)
    expected_term = np.zeros_like(got_term)

    for j in range(got_coef.shape[1]):
        w2 = wp2[0, j]
        x2 = xp2[0, j]
        y2 = yp2[0, j]
        wx = wpxp[0, j]
        wy = wpyp[0, j]
        sgn_wx = sgn_wpxp[0, j]
        sgn_wy = sgn_wpyp[0, j]
        mf = mixt_frac[0, j]
        fw = f_w[0, j]
        fx = f_x[0, j]
        fy = f_y[0, j]
        csw1 = coef_sigma_w_1_sqd[0, j]
        csw2 = coef_sigma_w_2_sqd[0, j]
        csx1 = coef_sigma_x_1_sqd[0, j]
        csx2 = coef_sigma_x_2_sqd[0, j]
        csy1 = coef_sigma_y_1_sqd[0, j]
        csy2 = coef_sigma_y_2_sqd[0, j]

        if (csw1 * csx1 > 0.0 or csw2 * csx2 > 0.0):
            coefs_factor_wx = (
                np.sqrt(csw1 * csx1) - np.sqrt(csw2 * csx2)
            ) / (
                mf * np.sqrt(csw1 * csx1) + (1.0 - mf) * np.sqrt(csw2 * csx2)
            )
        else:
            coefs_factor_wx = 0.0

        if (csw1 * csy1 > 0.0 or csw2 * csy2 > 0.0):
            coefs_factor_wy = (
                np.sqrt(csw1 * csy1) - np.sqrt(csw2 * csy2)
            ) / (
                mf * np.sqrt(csw1 * csy1) + (1.0 - mf) * np.sqrt(csw2 * csy2)
            )
        else:
            coefs_factor_wy = 0.0

        if (csx1 * csy1 > 0.0 or csx2 * csy2 > 0.0):
            coefs_factor_xy = (
                np.sqrt(csx1 * csy1) - np.sqrt(csx2 * csy2)
            ) / (
                mf * np.sqrt(csx1 * csy1) + (1.0 - mf) * np.sqrt(csx2 * csy2)
            )

            expected_coef[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fw) * np.sqrt(w2) * coefs_factor_xy
            )

            expected_term[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fw) * np.sqrt(w2)
                * np.sqrt(fx) * np.sqrt(x2) * sgn_wx * np.sqrt(fy) * np.sqrt(y2) * sgn_wy
                * (
                    (1.0 - mf) / mf - mf / (1.0 - mf)
                    - coefs_factor_xy - coefs_factor_wy - coefs_factor_wx
                )
                + np.sqrt(mf * (1.0 - mf)) * np.sqrt(fx) * np.sqrt(x2) * sgn_wx * coefs_factor_wy * wy
                + np.sqrt(mf * (1.0 - mf)) * np.sqrt(fy) * np.sqrt(y2) * sgn_wy * coefs_factor_wx * wx
            )
        else:
            expected_coef[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fw) * np.sqrt(w2)
                * ((1.0 - mf) / mf - mf / (1.0 - mf) - coefs_factor_wy - coefs_factor_wx)
            )
            expected_term[0, j] = (
                np.sqrt(mf * (1.0 - mf)) * np.sqrt(fx) * np.sqrt(x2) * sgn_wx * coefs_factor_wy * wy
                + np.sqrt(mf * (1.0 - mf)) * np.sqrt(fy) * np.sqrt(y2) * sgn_wy * coefs_factor_wx * wx
            )

    np.testing.assert_allclose(got_coef, expected_coef)
    np.testing.assert_allclose(got_term, expected_term)


def test_new_pdf_driver_column_consistency():
    """new_pdf_driver multi-column call should match per-column calls."""
    ncol, nz = 2, 3
    clubb_api.init_config_flags(clubb_api.get_default_config_flags())
    implicit = clubb_api.init_pdf_implicit(nz, ncol, 0)

    wm = np.array([[0.10, 0.20, 0.15], [-0.05, 0.00, 0.04]], dtype=np.float64)
    rtm = np.array([[0.012, 0.011, 0.010], [0.009, 0.010, 0.011]], dtype=np.float64)
    thlm = np.array([[300.0, 300.5, 301.0], [299.5, 300.0, 300.4]], dtype=np.float64)
    wp2 = np.array([[0.20, 0.22, 0.21], [0.18, 0.17, 0.19]], dtype=np.float64)
    rtp2 = np.array([[2e-6, 2.2e-6, 2.1e-6], [1.8e-6, 1.7e-6, 1.9e-6]], dtype=np.float64)
    thlp2 = np.array([[0.08, 0.09, 0.085], [0.07, 0.072, 0.074]], dtype=np.float64)
    skw = np.array([[0.2, -0.1, 0.15], [0.05, -0.2, 0.1]], dtype=np.float64)
    wprtp = np.array([[0.001, -0.0008, 0.0009], [0.0006, -0.0007, 0.0005]], dtype=np.float64)
    wpthlp = np.array([[0.02, -0.015, 0.018], [0.012, -0.010, 0.011]], dtype=np.float64)
    rtpthlp = np.array([[0.0002, -0.0001, 0.00015], [0.00012, -0.00008, 0.0001]], dtype=np.float64)
    clubb_params = clubb_api.init_clubb_params(ncol, iunit=10, filename="")
    skrt = np.array([[0.1, -0.05, 0.08], [0.02, -0.1, 0.04]], dtype=np.float64)
    skthl = np.array([[0.12, -0.08, 0.05], [0.03, -0.06, 0.02]], dtype=np.float64)

    got = clubb_api.new_pdf_driver(
        wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, rtpthlp, clubb_params, skrt, skthl,
        pdf_implicit_coefs_terms=implicit,
    )

    for i in range(ncol):
        implicit_i = clubb_api.init_pdf_implicit(nz, 1, 0)
        got_i = clubb_api.new_pdf_driver(
            wm[i:i+1, :], rtm[i:i+1, :], thlm[i:i+1, :], wp2[i:i+1, :],
            rtp2[i:i+1, :], thlp2[i:i+1, :], skw[i:i+1, :], wprtp[i:i+1, :],
            wpthlp[i:i+1, :], rtpthlp[i:i+1, :], clubb_params[i:i+1, :],
            skrt[i:i+1, :], skthl[i:i+1, :],
            pdf_implicit_coefs_terms=implicit_i,
        )
        for arr_multi, arr_single in zip(got, got_i):
            np.testing.assert_allclose(arr_multi[i:i+1, ...], arr_single)


def test_new_hybrid_pdf_driver_column_consistency():
    """new_hybrid_pdf_driver multi-column call should match per-column calls."""
    ncol, nz, sclr_dim = 2, 3, 1
    clubb_api.init_config_flags(clubb_api.get_default_config_flags())
    implicit = clubb_api.init_pdf_implicit(nz, ncol, sclr_dim)

    wm = np.array([[0.10, 0.20, 0.15], [-0.05, 0.00, 0.04]], dtype=np.float64)
    rtm = np.array([[0.012, 0.011, 0.010], [0.009, 0.010, 0.011]], dtype=np.float64)
    thlm = np.array([[300.0, 300.5, 301.0], [299.5, 300.0, 300.4]], dtype=np.float64)
    um = np.array([[5.0, 5.1, 5.2], [4.8, 4.9, 5.0]], dtype=np.float64)
    vm = np.array([[-1.0, -0.8, -0.9], [-1.1, -1.0, -0.95]], dtype=np.float64)
    wp2 = np.array([[0.20, 0.22, 0.21], [0.18, 0.17, 0.19]], dtype=np.float64)
    rtp2 = np.array([[2e-6, 2.2e-6, 2.1e-6], [1.8e-6, 1.7e-6, 1.9e-6]], dtype=np.float64)
    thlp2 = np.array([[0.08, 0.09, 0.085], [0.07, 0.072, 0.074]], dtype=np.float64)
    up2 = np.array([[0.15, 0.14, 0.16], [0.13, 0.12, 0.14]], dtype=np.float64)
    vp2 = np.array([[0.11, 0.10, 0.12], [0.10, 0.09, 0.11]], dtype=np.float64)
    skw = np.array([[0.2, -0.1, 0.15], [0.05, -0.2, 0.1]], dtype=np.float64)
    wprtp = np.array([[0.001, -0.0008, 0.0009], [0.0006, -0.0007, 0.0005]], dtype=np.float64)
    wpthlp = np.array([[0.02, -0.015, 0.018], [0.012, -0.010, 0.011]], dtype=np.float64)
    upwp = np.array([[0.04, 0.03, 0.035], [0.028, 0.026, 0.029]], dtype=np.float64)
    vpwp = np.array([[0.02, 0.018, 0.019], [0.016, 0.015, 0.017]], dtype=np.float64)
    sclrm = np.full((ncol, nz, sclr_dim), 0.001, dtype=np.float64)
    sclrp2 = np.full((ncol, nz, sclr_dim), 1.0e-6, dtype=np.float64)
    wpsclrp = np.full((ncol, nz, sclr_dim), 1.0e-4, dtype=np.float64)
    gamma_skw_fnc = np.array([[0.4, 0.45, 0.5], [0.35, 0.42, 0.47]], dtype=np.float64)
    slope_coef_spread_dg_means_w = np.array([1.2, 1.1], dtype=np.float64)
    pdf_component_stdev_factor_w = np.array([1.0, 1.0], dtype=np.float64)
    skrt = np.array([[0.1, -0.05, 0.08], [0.02, -0.1, 0.04]], dtype=np.float64)
    skthl = np.array([[0.12, -0.08, 0.05], [0.03, -0.06, 0.02]], dtype=np.float64)
    sku = np.array([[0.08, -0.04, 0.06], [0.01, -0.05, 0.03]], dtype=np.float64)
    skv = np.array([[0.06, -0.03, 0.05], [0.01, -0.04, 0.02]], dtype=np.float64)
    sksclr = np.full((ncol, nz, sclr_dim), 0.02, dtype=np.float64)

    got = clubb_api.new_hybrid_pdf_driver(
        wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp, upwp, vpwp,
        sclrm, sclrp2, wpsclrp, gamma_skw_fnc, slope_coef_spread_dg_means_w,
        pdf_component_stdev_factor_w, skrt, skthl, sku, skv, sksclr,
        pdf_implicit_coefs_terms=implicit,
    )

    for i in range(ncol):
        implicit_i = clubb_api.init_pdf_implicit(nz, 1, sclr_dim)
        got_i = clubb_api.new_hybrid_pdf_driver(
            wm[i:i+1, :], rtm[i:i+1, :], thlm[i:i+1, :], um[i:i+1, :], vm[i:i+1, :],
            wp2[i:i+1, :], rtp2[i:i+1, :], thlp2[i:i+1, :], up2[i:i+1, :], vp2[i:i+1, :],
            skw[i:i+1, :], wprtp[i:i+1, :], wpthlp[i:i+1, :], upwp[i:i+1, :], vpwp[i:i+1, :],
            sclrm[i:i+1, :, :], sclrp2[i:i+1, :, :], wpsclrp[i:i+1, :, :], gamma_skw_fnc[i:i+1, :],
            slope_coef_spread_dg_means_w[i:i+1], pdf_component_stdev_factor_w[i:i+1],
            skrt[i:i+1, :], skthl[i:i+1, :], sku[i:i+1, :], skv[i:i+1, :], sksclr[i:i+1, :, :],
            pdf_implicit_coefs_terms=implicit_i,
        )
        for arr_multi, arr_single in zip(got, got_i):
            np.testing.assert_allclose(arr_multi[i:i+1, ...], arr_single)

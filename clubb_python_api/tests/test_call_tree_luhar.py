"""Test wrappers for adg1_adg2_3d_luhar_pdf routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_calc_luhar_params_matches_formula():
    """calc_luhar_params should match the equations in adg1_adg2_3d_luhar_pdf.F90."""
    skx = np.array([[0.0, 0.5, -0.2, 0.9]], dtype=np.float64)
    wpxp = np.array([[0.1, -0.3, 0.0, 0.2]], dtype=np.float64)
    xp2 = np.array([[0.02, 0.3, 0.1, 0.06]], dtype=np.float64)
    x_tol_sqd = 0.05
    nz = skx.shape[1]
    ngrdcol = skx.shape[0]

    mixt_frac, big_m, small_m = clubb_api.calc_luhar_params(nz=nz, ngrdcol=ngrdcol, skx=skx, wpxp=wpxp, xp2=xp2, x_tol_sqd=x_tol_sqd)

    mask = xp2 > x_tol_sqd
    sgn_wpxp = np.where(wpxp >= 0.0, 1.0, -1.0)

    expected_small_m = np.where(
        mask,
        np.maximum((2.0 / 3.0) * np.abs(skx) ** (1.0 / 3.0), 0.05),
        0.0,
    )

    expected_big_m = np.zeros_like(skx)
    small_m_sqd = expected_small_m[mask] ** 2
    expected_big_m[mask] = (1.0 + small_m_sqd) ** 3 / ((3.0 + small_m_sqd) ** 2 * small_m_sqd)

    expected_mixt_frac = np.full_like(skx, 0.5)
    expected_mixt_frac[mask] = 0.5 * (
        1.0
        - sgn_wpxp[mask]
        * skx[mask]
        * np.sqrt(1.0 / ((4.0 / expected_big_m[mask]) + skx[mask] ** 2))
    )

    np.testing.assert_allclose(small_m, expected_small_m)
    np.testing.assert_allclose(big_m, expected_big_m)
    np.testing.assert_allclose(mixt_frac, expected_mixt_frac)


def test_close_luhar_pdf_matches_formula():
    """close_luhar_pdf should match the equations in adg1_adg2_3d_luhar_pdf.F90."""
    xm = np.array([[1.0, -0.2, 0.7, 0.1]], dtype=np.float64)
    xp2 = np.array([[0.2, 0.06, 0.4, 0.08]], dtype=np.float64)
    mixt_frac = np.array([[0.35, 0.6, 0.42, 0.5]], dtype=np.float64)
    small_m = np.array([[0.4, 0.1, 0.55, 0.25]], dtype=np.float64)
    wpxp = np.array([[0.07, -0.05, 0.0, 0.03]], dtype=np.float64)
    x_tol_sqd = 0.05
    nz = xm.shape[1]
    ngrdcol = xm.shape[0]

    (
        sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2,
        x_1_n, x_2_n, x_1, x_2,
    ) = clubb_api.close_luhar_pdf(
        nz=nz, ngrdcol=ngrdcol, xm=xm, xp2=xp2, mixt_frac=mixt_frac, small_m=small_m, wpxp=wpxp,
        x_tol_sqd=x_tol_sqd)

    mask = xp2 > x_tol_sqd
    sgn_wpxp = np.where(wpxp >= 0.0, 1.0, -1.0)
    one_plus_m2 = 1.0 + small_m**2

    expected_sigma_1 = np.where(mask, (1.0 - mixt_frac) / (mixt_frac * one_plus_m2), 1.0 / one_plus_m2)
    expected_sigma_2 = np.where(mask, mixt_frac / ((1.0 - mixt_frac) * one_plus_m2), 1.0 / one_plus_m2)

    expected_var_1 = np.where(mask, expected_sigma_1 * xp2, 0.0)
    expected_var_2 = np.where(mask, expected_sigma_2 * xp2, 0.0)

    expected_x_1_n = sgn_wpxp * small_m * np.sqrt(expected_sigma_1)
    expected_x_2_n = -sgn_wpxp * small_m * np.sqrt(expected_sigma_2)

    sqrt_xp2 = np.sqrt(np.where(mask, xp2, 0.0))
    expected_x_1 = np.where(mask, xm + sqrt_xp2 * expected_x_1_n, xm)
    expected_x_2 = np.where(mask, xm + sqrt_xp2 * expected_x_2_n, xm)

    np.testing.assert_allclose(sigma_sqd_x_1, expected_sigma_1)
    np.testing.assert_allclose(sigma_sqd_x_2, expected_sigma_2)
    np.testing.assert_allclose(varnce_x_1, expected_var_1)
    np.testing.assert_allclose(varnce_x_2, expected_var_2)
    np.testing.assert_allclose(x_1_n, expected_x_1_n)
    np.testing.assert_allclose(x_2_n, expected_x_2_n)
    np.testing.assert_allclose(x_1, expected_x_1)
    np.testing.assert_allclose(x_2, expected_x_2)


def test_adg1_w_closure_matches_formula():
    """adg1_w_closure should match the equations in adg1_adg2_3d_luhar_pdf.F90."""
    wm = np.array([[0.2, -0.1, 0.6, 0.0]], dtype=np.float64)
    wp2 = np.array([[0.5, 0.0, 0.3, 0.9]], dtype=np.float64)
    skw = np.array([[0.2, -0.4, 1.0e-7, -0.1]], dtype=np.float64)
    sigma_sqd_w = np.array([[0.4, 0.3, 0.6, 0.5]], dtype=np.float64)
    sqrt_wp2 = np.sqrt(wp2)
    mixt_frac_max_mag = 0.95
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac = clubb_api.adg1_w_closure(
        nz=nz, ngrdcol=ngrdcol, wm=wm, wp2=wp2, skw=skw, sigma_sqd_w=sigma_sqd_w, sqrt_wp2=sqrt_wp2,
        mixt_frac_max_mag=mixt_frac_max_mag
    )

    base_mixt = np.where(
        np.abs(skw) <= 1.0e-5,
        0.5,
        0.5 * (1.0 - skw / np.sqrt(4.0 * (1.0 - sigma_sqd_w) ** 3 + skw**2)),
    )
    clipped_mixt = np.minimum(np.maximum(base_mixt, 1.0 - mixt_frac_max_mag), mixt_frac_max_mag)

    expected_mixt_frac = clipped_mixt
    expected_w_1_n = np.sqrt(((1.0 - expected_mixt_frac) / expected_mixt_frac) * (1.0 - sigma_sqd_w))
    expected_w_2_n = -np.sqrt((expected_mixt_frac / (1.0 - expected_mixt_frac)) * (1.0 - sigma_sqd_w))

    expected_w_1 = wm + sqrt_wp2 * expected_w_1_n
    expected_w_2 = wm + sqrt_wp2 * expected_w_2_n
    expected_var_1 = sigma_sqd_w * wp2
    expected_var_2 = sigma_sqd_w * wp2

    np.testing.assert_allclose(mixt_frac, expected_mixt_frac)
    np.testing.assert_allclose(w_1_n, expected_w_1_n)
    np.testing.assert_allclose(w_2_n, expected_w_2_n)
    np.testing.assert_allclose(w_1, expected_w_1)
    np.testing.assert_allclose(w_2, expected_w_2)
    np.testing.assert_allclose(varnce_w_1, expected_var_1)
    np.testing.assert_allclose(varnce_w_2, expected_var_2)


def test_adg2_pdf_driver_column_consistency():
    """adg2_pdf_driver multi-column call should match per-column calls."""
    ncol, nz, sclr_dim = 2, 3, 1
    sclr_tol = np.array([1.0e-6], dtype=np.float64)
    wm = np.array([[0.10, 0.20, 0.15], [-0.05, 0.00, 0.04]], dtype=np.float64)
    rtm = np.array([[0.012, 0.011, 0.010], [0.009, 0.010, 0.011]], dtype=np.float64)
    thlm = np.array([[300.0, 300.5, 301.0], [299.5, 300.0, 300.4]], dtype=np.float64)
    wp2 = np.array([[0.20, 0.22, 0.21], [0.18, 0.17, 0.19]], dtype=np.float64)
    rtp2 = np.array([[2e-6, 2.2e-6, 2.1e-6], [1.8e-6, 1.7e-6, 1.9e-6]], dtype=np.float64)
    thlp2 = np.array([[0.08, 0.09, 0.085], [0.07, 0.072, 0.074]], dtype=np.float64)
    skw = np.array([[0.2, -0.1, 0.15], [0.05, -0.2, 0.1]], dtype=np.float64)
    wprtp = np.array([[0.001, -0.0008, 0.0009], [0.0006, -0.0007, 0.0005]], dtype=np.float64)
    wpthlp = np.array([[0.02, -0.015, 0.018], [0.012, -0.010, 0.011]], dtype=np.float64)
    sqrt_wp2 = np.sqrt(wp2)
    beta = np.array([1.5, 1.6], dtype=np.float64)
    sclrm = np.full((ncol, nz, sclr_dim), 0.001, dtype=np.float64)
    sclrp2 = np.full((ncol, nz, sclr_dim), 1.0e-6, dtype=np.float64)
    wpsclrp = np.full((ncol, nz, sclr_dim), 1.0e-4, dtype=np.float64)

    got = clubb_api.adg2_pdf_driver(
        nz=nz, ngrdcol=ncol, sclr_dim=sclr_dim, sclr_tol=sclr_tol, wm=wm, rtm=rtm, thlm=thlm, wp2=wp2,
        rtp2=rtp2, thlp2=thlp2, skw=skw, wprtp=wprtp, wpthlp=wpthlp, sqrt_wp2=sqrt_wp2, beta=beta,
        sclrm=sclrm, sclrp2=sclrp2, wpsclrp=wpsclrp, l_scalar_calc=True,
    )

    for i in range(ncol):
        got_i = clubb_api.adg2_pdf_driver(
            nz=nz, ngrdcol=1, sclr_dim=sclr_dim, sclr_tol=sclr_tol, wm=wm[i:i+1, :], rtm=rtm[i:i+1, :],
            thlm=thlm[i:i+1, :], wp2=wp2[i:i+1, :], rtp2=rtp2[i:i+1, :], thlp2=thlp2[i:i+1, :],
            skw=skw[i:i+1, :], wprtp=wprtp[i:i+1, :], wpthlp=wpthlp[i:i+1, :], sqrt_wp2=sqrt_wp2[i:i+1, :],
            beta=beta[i:i+1], sclrm=sclrm[i:i+1, :, :], sclrp2=sclrp2[i:i+1, :, :],
            wpsclrp=wpsclrp[i:i+1, :, :], l_scalar_calc=True,
        )
        for arr_multi, arr_single in zip(got, got_i):
            np.testing.assert_allclose(arr_multi[i:i+1, ...], arr_single)


def test_adg1_pdf_driver_column_consistency():
    """adg1_pdf_driver multi-column call should match per-column calls."""
    ncol, nz, sclr_dim = 2, 3, 1
    sclr_tol = np.array([1.0e-6], dtype=np.float64)
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
    sqrt_wp2 = np.sqrt(wp2)
    sigma_sqd_w = np.array([[0.4, 0.42, 0.41], [0.39, 0.38, 0.40]], dtype=np.float64)
    beta = np.array([1.5, 1.6], dtype=np.float64)
    mixt_frac_max_mag = 0.95
    sclrm = np.full((ncol, nz, sclr_dim), 0.001, dtype=np.float64)
    sclrp2 = np.full((ncol, nz, sclr_dim), 1.0e-6, dtype=np.float64)
    wpsclrp = np.full((ncol, nz, sclr_dim), 1.0e-4, dtype=np.float64)

    got = clubb_api.adg1_pdf_driver(
        nz=nz, ngrdcol=ncol, sclr_dim=sclr_dim, sclr_tol=sclr_tol, wm=wm, rtm=rtm, thlm=thlm, um=um, vm=vm,
        wp2=wp2, rtp2=rtp2, thlp2=thlp2, up2=up2, vp2=vp2, skw=skw, wprtp=wprtp, wpthlp=wpthlp, upwp=upwp,
        vpwp=vpwp, sqrt_wp2=sqrt_wp2, sigma_sqd_w=sigma_sqd_w, beta=beta, mixt_frac_max_mag=mixt_frac_max_mag,
        sclrm=sclrm, sclrp2=sclrp2, wpsclrp=wpsclrp, l_scalar_calc=True,
    )

    for i in range(ncol):
        got_i = clubb_api.adg1_pdf_driver(
            nz=nz, ngrdcol=1, sclr_dim=sclr_dim, sclr_tol=sclr_tol, wm=wm[i:i+1, :], rtm=rtm[i:i+1, :],
            thlm=thlm[i:i+1, :], um=um[i:i+1, :], vm=vm[i:i+1, :], wp2=wp2[i:i+1, :], rtp2=rtp2[i:i+1, :],
            thlp2=thlp2[i:i+1, :], up2=up2[i:i+1, :], vp2=vp2[i:i+1, :], skw=skw[i:i+1, :],
            wprtp=wprtp[i:i+1, :], wpthlp=wpthlp[i:i+1, :], upwp=upwp[i:i+1, :], vpwp=vpwp[i:i+1, :],
            sqrt_wp2=sqrt_wp2[i:i+1, :], sigma_sqd_w=sigma_sqd_w[i:i+1, :], beta=beta[i:i+1],
            mixt_frac_max_mag=mixt_frac_max_mag, sclrm=sclrm[i:i+1, :, :], sclrp2=sclrp2[i:i+1, :, :],
            wpsclrp=wpsclrp[i:i+1, :, :], l_scalar_calc=True,
        )
        for arr_multi, arr_single in zip(got, got_i):
            np.testing.assert_allclose(arr_multi[i:i+1, ...], arr_single)


def test_luhar_3d_pdf_driver_column_consistency():
    """luhar_3d_pdf_driver multi-column call should match per-column calls."""
    ncol, nz = 2, 3
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

    got = clubb_api.luhar_3d_pdf_driver(
        nz=nz, ngrdcol=ncol, wm=wm, rtm=rtm, thlm=thlm, wp2=wp2, rtp2=rtp2, thlp2=thlp2, skw=skw, skrt=skrt,
        skthl=skthl, wprtp=wprtp, wpthlp=wpthlp)

    for i in range(ncol):
        got_i = clubb_api.luhar_3d_pdf_driver(
            nz=nz, ngrdcol=1, wm=wm[i:i+1, :], rtm=rtm[i:i+1, :], thlm=thlm[i:i+1, :], wp2=wp2[i:i+1, :],
            rtp2=rtp2[i:i+1, :], thlp2=thlp2[i:i+1, :], skw=skw[i:i+1, :], skrt=skrt[i:i+1, :],
            skthl=skthl[i:i+1, :], wprtp=wprtp[i:i+1, :], wpthlp=wpthlp[i:i+1, :],
        )
        for arr_multi, arr_single in zip(got, got_i):
            np.testing.assert_allclose(arr_multi[i:i+1, :], arr_single)

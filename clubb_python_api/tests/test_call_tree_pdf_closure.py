"""Test wrappers for pdf_closure_module moment routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_calc_wp4_pdf_matches_formula():
    """calc_wp4_pdf should match the integrated two-component Gaussian formula."""
    wm = np.array([[0.2, -0.3]], dtype=np.float64)
    w_1 = np.array([[0.6, -0.1]], dtype=np.float64)
    w_2 = np.array([[-0.4, -0.7]], dtype=np.float64)
    varnce_w_1 = np.array([[0.5, 0.9]], dtype=np.float64)
    varnce_w_2 = np.array([[0.8, 0.3]], dtype=np.float64)
    mixt_frac = np.array([[0.25, 0.6]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    got = clubb_api.calc_wp4_pdf(
        nz=nz, ngrdcol=ngrdcol, wm=wm, w_1=w_1, w_2=w_2, varnce_w_1=varnce_w_1, varnce_w_2=varnce_w_2,
        mixt_frac=mixt_frac)

    expected = (
        mixt_frac * (
            3.0 * varnce_w_1**2
            + 6.0 * (w_1 - wm) ** 2 * varnce_w_1
            + (w_1 - wm) ** 4
        )
        + (1.0 - mixt_frac) * (
            3.0 * varnce_w_2**2
            + 6.0 * (w_2 - wm) ** 2 * varnce_w_2
            + (w_2 - wm) ** 4
        )
    )

    np.testing.assert_allclose(got, expected)


def test_calc_wp2xp_pdf_matches_formula():
    """calc_wp2xp_pdf should match its analytic bivariate Gaussian expression."""
    wm = np.array([[0.3, -0.2]], dtype=np.float64)
    xm = np.array([[1.0, 0.4]], dtype=np.float64)
    w_1 = np.array([[0.8, -0.1]], dtype=np.float64)
    w_2 = np.array([[-0.5, -0.4]], dtype=np.float64)
    x_1 = np.array([[1.4, 0.6]], dtype=np.float64)
    x_2 = np.array([[0.7, 0.2]], dtype=np.float64)
    varnce_w_1 = np.array([[0.5, 0.8]], dtype=np.float64)
    varnce_w_2 = np.array([[0.3, 0.7]], dtype=np.float64)
    varnce_x_1 = np.array([[0.9, 0.6]], dtype=np.float64)
    varnce_x_2 = np.array([[0.4, 0.5]], dtype=np.float64)
    corr_w_x_1 = np.array([[0.2, -0.3]], dtype=np.float64)
    corr_w_x_2 = np.array([[-0.1, 0.15]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.55]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    got = clubb_api.calc_wp2xp_pdf(
        nz=nz, ngrdcol=ngrdcol, wm=wm, xm=xm, w_1=w_1, w_2=w_2, x_1=x_1, x_2=x_2, varnce_w_1=varnce_w_1,
        varnce_w_2=varnce_w_2, varnce_x_1=varnce_x_1, varnce_x_2=varnce_x_2, corr_w_x_1=corr_w_x_1,
        corr_w_x_2=corr_w_x_2, mixt_frac=mixt_frac,
    )

    expected = (
        mixt_frac * (
            (((w_1 - wm) ** 2) + varnce_w_1) * (x_1 - xm)
            + 2.0 * corr_w_x_1 * np.sqrt(varnce_w_1 * varnce_x_1) * (w_1 - wm)
        )
        + (1.0 - mixt_frac) * (
            (((w_2 - wm) ** 2) + varnce_w_2) * (x_2 - xm)
            + 2.0 * corr_w_x_2 * np.sqrt(varnce_w_2 * varnce_x_2) * (w_2 - wm)
        )
    )

    np.testing.assert_allclose(got, expected)


def test_calc_wpxp2_pdf_matches_formula():
    """calc_wpxp2_pdf should match its analytic bivariate Gaussian expression."""
    wm = np.array([[0.3, -0.2]], dtype=np.float64)
    xm = np.array([[1.0, 0.4]], dtype=np.float64)
    w_1 = np.array([[0.8, -0.1]], dtype=np.float64)
    w_2 = np.array([[-0.5, -0.4]], dtype=np.float64)
    x_1 = np.array([[1.4, 0.6]], dtype=np.float64)
    x_2 = np.array([[0.7, 0.2]], dtype=np.float64)
    varnce_w_1 = np.array([[0.5, 0.8]], dtype=np.float64)
    varnce_w_2 = np.array([[0.3, 0.7]], dtype=np.float64)
    varnce_x_1 = np.array([[0.9, 0.6]], dtype=np.float64)
    varnce_x_2 = np.array([[0.4, 0.5]], dtype=np.float64)
    corr_w_x_1 = np.array([[0.2, -0.3]], dtype=np.float64)
    corr_w_x_2 = np.array([[-0.1, 0.15]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.55]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    got = clubb_api.calc_wpxp2_pdf(
        nz=nz, ngrdcol=ngrdcol, wm=wm, xm=xm, w_1=w_1, w_2=w_2, x_1=x_1, x_2=x_2, varnce_w_1=varnce_w_1,
        varnce_w_2=varnce_w_2, varnce_x_1=varnce_x_1, varnce_x_2=varnce_x_2, corr_w_x_1=corr_w_x_1,
        corr_w_x_2=corr_w_x_2, mixt_frac=mixt_frac,
    )

    expected = (
        mixt_frac * (
            (w_1 - wm) * (((x_1 - xm) ** 2) + varnce_x_1)
            + 2.0 * corr_w_x_1 * np.sqrt(varnce_w_1 * varnce_x_1) * (x_1 - xm)
        )
        + (1.0 - mixt_frac) * (
            (w_2 - wm) * (((x_2 - xm) ** 2) + varnce_x_2)
            + 2.0 * corr_w_x_2 * np.sqrt(varnce_w_2 * varnce_x_2) * (x_2 - xm)
        )
    )

    np.testing.assert_allclose(got, expected)


def test_calc_wpxpyp_pdf_matches_formula():
    """calc_wpxpyp_pdf should match its analytic trivariate Gaussian expression."""
    wm = np.array([[0.3, -0.2]], dtype=np.float64)
    xm = np.array([[1.0, 0.4]], dtype=np.float64)
    ym = np.array([[-0.4, 0.2]], dtype=np.float64)
    w_1 = np.array([[0.8, -0.1]], dtype=np.float64)
    w_2 = np.array([[-0.5, -0.4]], dtype=np.float64)
    x_1 = np.array([[1.4, 0.6]], dtype=np.float64)
    x_2 = np.array([[0.7, 0.2]], dtype=np.float64)
    y_1 = np.array([[-0.2, 0.5]], dtype=np.float64)
    y_2 = np.array([[-0.6, -0.1]], dtype=np.float64)
    varnce_w_1 = np.array([[0.5, 0.8]], dtype=np.float64)
    varnce_w_2 = np.array([[0.3, 0.7]], dtype=np.float64)
    varnce_x_1 = np.array([[0.9, 0.6]], dtype=np.float64)
    varnce_x_2 = np.array([[0.4, 0.5]], dtype=np.float64)
    varnce_y_1 = np.array([[0.7, 0.3]], dtype=np.float64)
    varnce_y_2 = np.array([[0.5, 0.4]], dtype=np.float64)
    corr_w_x_1 = np.array([[0.2, -0.3]], dtype=np.float64)
    corr_w_x_2 = np.array([[-0.1, 0.15]], dtype=np.float64)
    corr_w_y_1 = np.array([[0.1, 0.05]], dtype=np.float64)
    corr_w_y_2 = np.array([[-0.2, 0.1]], dtype=np.float64)
    corr_x_y_1 = np.array([[0.12, -0.08]], dtype=np.float64)
    corr_x_y_2 = np.array([[-0.04, 0.2]], dtype=np.float64)
    mixt_frac = np.array([[0.4, 0.55]], dtype=np.float64)
    nz = wm.shape[1]
    ngrdcol = wm.shape[0]

    got = clubb_api.calc_wpxpyp_pdf(
        nz=nz, ngrdcol=ngrdcol, wm=wm, xm=xm, ym=ym, w_1=w_1, w_2=w_2, x_1=x_1, x_2=x_2, y_1=y_1, y_2=y_2,
        varnce_w_1=varnce_w_1, varnce_w_2=varnce_w_2, varnce_x_1=varnce_x_1, varnce_x_2=varnce_x_2,
        varnce_y_1=varnce_y_1, varnce_y_2=varnce_y_2, corr_w_x_1=corr_w_x_1, corr_w_x_2=corr_w_x_2,
        corr_w_y_1=corr_w_y_1, corr_w_y_2=corr_w_y_2, corr_x_y_1=corr_x_y_1, corr_x_y_2=corr_x_y_2,
        mixt_frac=mixt_frac,
    )

    expected = (
        mixt_frac * (
            (w_1 - wm) * (x_1 - xm) * (y_1 - ym)
            + corr_x_y_1 * np.sqrt(varnce_x_1 * varnce_y_1) * (w_1 - wm)
            + corr_w_y_1 * np.sqrt(varnce_w_1 * varnce_y_1) * (x_1 - xm)
            + corr_w_x_1 * np.sqrt(varnce_w_1 * varnce_x_1) * (y_1 - ym)
        )
        + (1.0 - mixt_frac) * (
            (w_2 - wm) * (x_2 - xm) * (y_2 - ym)
            + corr_x_y_2 * np.sqrt(varnce_x_2 * varnce_y_2) * (w_2 - wm)
            + corr_w_y_2 * np.sqrt(varnce_w_2 * varnce_y_2) * (x_2 - xm)
            + corr_w_x_2 * np.sqrt(varnce_w_2 * varnce_x_2) * (y_2 - ym)
        )
    )

    np.testing.assert_allclose(got, expected)


def test_calc_w_up_in_cloud_zero_variance_component_split():
    """When each component is entirely updraft or downdraft, outputs should reduce exactly."""
    mixt_frac = np.array([[0.3, 0.8]], dtype=np.float64)
    cloud_frac_1 = np.array([[0.5, 0.4]], dtype=np.float64)
    cloud_frac_2 = np.array([[0.2, 0.7]], dtype=np.float64)
    w_1 = np.array([[1.2, 0.6]], dtype=np.float64)
    w_2 = np.array([[-0.7, -1.5]], dtype=np.float64)
    varnce_w_1 = np.zeros((1, 2), dtype=np.float64)
    varnce_w_2 = np.zeros((1, 2), dtype=np.float64)
    nz = mixt_frac.shape[1]
    ngrdcol = mixt_frac.shape[0]

    w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, cloudy_downdraft_frac = clubb_api.calc_w_up_in_cloud(
        nz=nz, ngrdcol=ngrdcol, mixt_frac=mixt_frac, cloud_frac_1=cloud_frac_1, cloud_frac_2=cloud_frac_2,
        w_1=w_1, w_2=w_2, varnce_w_1=varnce_w_1, varnce_w_2=varnce_w_2,
    )

    expected_cloudy_updraft_frac = mixt_frac * cloud_frac_1
    expected_cloudy_downdraft_frac = (1.0 - mixt_frac) * cloud_frac_2

    np.testing.assert_allclose(w_up_in_cloud, w_1)
    np.testing.assert_allclose(w_down_in_cloud, w_2)
    np.testing.assert_allclose(cloudy_updraft_frac, expected_cloudy_updraft_frac)
    np.testing.assert_allclose(cloudy_downdraft_frac, expected_cloudy_downdraft_frac)

"""Test wrappers for Skx_module call-tree routines."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_skx_func_matches_xp3_over_xp2_to_three_halves_when_xtol_zero():
    """With x_tol=0, Skx should reduce to xp3 / xp2^(3/2)."""
    ngrdcol, nz = 1, 4
    xp2 = np.array([[1.0, 4.0, 9.0, 16.0]], dtype=np.float64)
    xp3 = np.array([[1.0, 8.0, 27.0, 64.0]], dtype=np.float64)
    clubb_params = np.ones_like(clubb_api.init_clubb_params(ngrdcol, iunit=10, filename=""))

    skx = clubb_api.skx_func(nz=nz, ngrdcol=ngrdcol, xp2=xp2, xp3=xp3, x_tol=0.0, clubb_params=clubb_params)
    expected = xp3 / (xp2 * np.sqrt(xp2))

    assert skx.shape == (ngrdcol, nz)
    np.testing.assert_allclose(skx, expected)


def test_lg_2005_ansatz_zero_flux_gives_zero_skewness():
    """If wpxp is zero, LG05 ansatz should produce zero Skx."""
    ngrdcol, nz = 1, 5
    skw = np.full((ngrdcol, nz), 0.8, dtype=np.float64)
    wpxp = np.zeros((ngrdcol, nz), dtype=np.float64)
    wp2 = np.full((ngrdcol, nz), 1.5, dtype=np.float64)
    xp2 = np.full((ngrdcol, nz), 2.0, dtype=np.float64)
    beta = np.full(ngrdcol, 0.7, dtype=np.float64)
    sigma_sqd_w = np.full((ngrdcol, nz), 0.25, dtype=np.float64)

    skx = clubb_api.lg_2005_ansatz(
        nz=nz, ngrdcol=ngrdcol, skw=skw, wpxp=wpxp, wp2=wp2, xp2=xp2, beta=beta,
        sigma_sqd_w=sigma_sqd_w, x_tol=1.0e-6)

    assert skx.shape == (ngrdcol, nz)
    np.testing.assert_allclose(skx, 0.0)


def test_xp3_lg_2005_ansatz_matches_two_step_reconstruction():
    """xp3_LG_2005 should match LG05 Skx followed by xp3 reconstruction."""
    ngrdcol, nzt = 1, 4
    skw_zt = np.array([[0.2, 0.4, 0.6, 0.8]], dtype=np.float64)
    wpxp_zt = np.array([[0.05, 0.10, 0.12, 0.15]], dtype=np.float64)
    wp2_zt = np.array([[0.9, 1.0, 1.1, 1.2]], dtype=np.float64)
    xp2_zt = np.array([[0.7, 0.8, 0.9, 1.0]], dtype=np.float64)
    sigma_sqd_w_zt = np.full((ngrdcol, nzt), 0.2, dtype=np.float64)
    clubb_params = np.ones_like(clubb_api.init_clubb_params(ngrdcol, iunit=10, filename=""))

    skx = clubb_api.lg_2005_ansatz(
        nz=nzt, ngrdcol=ngrdcol, skw=skw_zt, wpxp=wpxp_zt, wp2=wp2_zt, xp2=xp2_zt,
        beta=np.ones(ngrdcol, dtype=np.float64),
        sigma_sqd_w=sigma_sqd_w_zt,
        x_tol=0.0,
    )
    expected_xp3 = skx * xp2_zt * np.sqrt(xp2_zt)

    xp3 = clubb_api.xp3_lg_2005_ansatz(
        nzt=nzt, ngrdcol=ngrdcol, skw_zt=skw_zt, wpxp_zt=wpxp_zt, wp2_zt=wp2_zt, xp2_zt=xp2_zt,
        sigma_sqd_w_zt=sigma_sqd_w_zt, clubb_params=clubb_params, x_tol=0.0,
    )

    assert xp3.shape == (ngrdcol, nzt)
    np.testing.assert_allclose(xp3, expected_xp3)

"""Test wrappers for transport/mixing routines in the call tree."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _make_grid(ngrdcol=1, nzmax=11, dz=500.0):
    """Create a standard test grid."""
    clubb_api.init_err_info(ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_config_flags(flags)

    zm = np.arange(nzmax, dtype=np.float64) * dz
    zt = 0.5 * (zm[:-1] + zm[1:])

    gr, _err_info = clubb_api.setup_grid(
        nzmax=nzmax, ngrdcol=ngrdcol,
        sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
        l_implemented=True, l_ascending_grid=True,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        zm_init=np.zeros(ngrdcol, dtype=np.float64),
        zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
        momentum_heights=np.tile(zm, (ngrdcol, 1)),
        thermodynamic_heights=np.tile(zt, (ngrdcol, 1)),
        err_info=ErrInfo(ngrdcol=ngrdcol),
    )
    return gr


@pytest.fixture(scope="module")
def gr():
    """Shared test grid."""
    return _make_grid()


def test_compute_current_date_one_day_and_hour():
    """86400+3600 seconds after Jan 1, 2000 should be Jan 2, 2000 + 1 hour."""
    day, month, year, sec = clubb_api.compute_current_date(1, 1, 2000, 90000.0)
    assert int(day) == 2
    assert int(month) == 1
    assert int(year) == 2000
    np.testing.assert_allclose(sec, 3600.0)


def test_set_lscale_max_behavior():
    """set_lscale_max should follow implemented/non-implemented formulas."""
    host_dx = np.array([1000.0, 2000.0], dtype=np.float64)
    host_dy = np.array([500.0, 4000.0], dtype=np.float64)

    unimplemented = clubb_api.set_lscale_max(ngrdcol=host_dx.size, host_dx=host_dx, host_dy=host_dy, l_implemented=False)
    implemented = clubb_api.set_lscale_max(ngrdcol=host_dx.size, host_dx=host_dx, host_dy=host_dy, l_implemented=True)

    np.testing.assert_allclose(unimplemented, np.array([1.0e5, 1.0e5]))
    np.testing.assert_allclose(implemented, 0.25 * np.minimum(host_dx, host_dy))


def test_term_ma_zt_lhs_zero_w(gr):
    """Zero mean vertical velocity should produce zero zt lhs contribution."""
    wm_zt = np.zeros((gr.ngrdcol, gr.nzt), dtype=np.float64)
    lhs = clubb_api.term_ma_zt_lhs(
        nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        wm_zt=wm_zt,
        weights_zt2zm=gr.weights_zt2zm,
        invrs_dzt=gr.invrs_dzt,
        invrs_dzm=gr.invrs_dzm,
        l_upwind_xm_ma=False,
        grid_dir=gr.grid_dir,
    )
    assert lhs.shape == (3, gr.ngrdcol, gr.nzt)
    np.testing.assert_allclose(lhs, 0.0)


def test_term_ma_zm_lhs_zero_w(gr):
    """Zero mean vertical velocity should produce zero zm lhs contribution."""
    wm_zm = np.zeros((gr.ngrdcol, gr.nzm), dtype=np.float64)
    lhs = clubb_api.term_ma_zm_lhs(
        nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, wm_zm=wm_zm,
        invrs_dzm=gr.invrs_dzm, weights_zm2zt=gr.weights_zm2zt,
    )
    assert lhs.shape == (3, gr.ngrdcol, gr.nzm)
    np.testing.assert_allclose(lhs, 0.0)


def test_diffusion_zt_lhs_zero_coeffs(gr):
    """Zero diffusion coefficients should yield zero zt diffusion lhs."""
    k_zm = np.zeros((gr.ngrdcol, gr.nzm), dtype=np.float64)
    k_zt = np.zeros((gr.ngrdcol, gr.nzt), dtype=np.float64)
    nu = np.zeros(gr.ngrdcol, dtype=np.float64)
    invrs_rho_ds_zt = np.ones((gr.ngrdcol, gr.nzt), dtype=np.float64)
    rho_ds_zm = np.ones((gr.ngrdcol, gr.nzm), dtype=np.float64)

    lhs = clubb_api.diffusion_zt_lhs(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        k_zm=k_zm, k_zt=k_zt, nu=nu, invrs_rho_ds_zt=invrs_rho_ds_zt, rho_ds_zm=rho_ds_zm
    )
    assert lhs.shape == (3, gr.ngrdcol, gr.nzt)
    np.testing.assert_allclose(lhs, 0.0)


def test_diffusion_zm_lhs_zero_coeffs(gr):
    """Zero diffusion coefficients should yield zero zm diffusion lhs."""
    k_zm = np.zeros((gr.ngrdcol, gr.nzm), dtype=np.float64)
    k_zt = np.zeros((gr.ngrdcol, gr.nzt), dtype=np.float64)
    nu = np.zeros(gr.ngrdcol, dtype=np.float64)
    invrs_rho_ds_zm = np.ones((gr.ngrdcol, gr.nzm), dtype=np.float64)
    rho_ds_zt = np.ones((gr.ngrdcol, gr.nzt), dtype=np.float64)

    lhs = clubb_api.diffusion_zm_lhs(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        k_zm=k_zm, k_zt=k_zt, nu=nu, invrs_rho_ds_zm=invrs_rho_ds_zm, rho_ds_zt=rho_ds_zt
    )
    assert lhs.shape == (3, gr.ngrdcol, gr.nzm)
    np.testing.assert_allclose(lhs, 0.0)

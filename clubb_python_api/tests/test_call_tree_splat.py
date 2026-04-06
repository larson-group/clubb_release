"""Test wrappers for additional advance_helper_module call-tree routines."""
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


def test_smooth_heaviside_peskin_limits_and_center():
    """Smoothed Heaviside should saturate outside range and be 0.5 at zero."""
    input_var = np.array([[-2.0, -1.0, 0.0, 1.0, 2.0]], dtype=np.float64)
    out = clubb_api.smooth_heaviside_peskin(
        nz=input_var.shape[1], ngrdcol=input_var.shape[0], input=input_var, smth_range=1.0
    )
    assert out.shape == input_var.shape
    np.testing.assert_allclose(out[0, 0], 0.0)
    np.testing.assert_allclose(out[0, -1], 1.0)
    np.testing.assert_allclose(out[0, 2], 0.5)
    assert np.all(np.diff(out[0, :]) >= -1e-15)


def test_lscale_width_vert_avg_constant_profile(gr):
    """Constant profile should remain constant under Lscale averaging."""
    shape = (gr.ngrdcol, gr.nzm)
    var_profile = np.full(shape, 7.5, dtype=np.float64)
    lscale_zm = np.full(shape, 600.0, dtype=np.float64)
    rho_ds_zm = np.ones(shape, dtype=np.float64)

    avg = clubb_api.lscale_width_vert_avg(
        gr, nzm=gr.nzm, ngrdcol=gr.ngrdcol, smth_type=1, var_profile=var_profile,
        lscale_zm=lscale_zm, rho_ds_zm=rho_ds_zm,
        var_below_ground_value=7.5,
    )
    assert avg.shape == shape
    np.testing.assert_allclose(avg, var_profile)


def test_wp2_term_splat_lhs_zero_bv(gr):
    """Zero Brunt-Vaisala squared should yield zero wp2 splat lhs."""
    c_wp2_splat = np.full(gr.ngrdcol, 0.7, dtype=np.float64)
    bv_sqd = np.zeros((gr.ngrdcol, gr.nzm), dtype=np.float64)
    lhs = clubb_api.wp2_term_splat_lhs(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        c_wp2_splat=c_wp2_splat, brunt_vaisala_freq_sqd_splat=bv_sqd
    )
    assert lhs.shape == (gr.ngrdcol, gr.nzm)
    np.testing.assert_allclose(lhs, 0.0)


def test_wp2_term_splat_lhs_constant_bv(gr):
    """Constant Brunt-Vaisala profile should preserve analytic constant result."""
    c_wp2_splat = np.full(gr.ngrdcol, 0.5, dtype=np.float64)
    bv_sqd = np.full((gr.ngrdcol, gr.nzm), 4.0, dtype=np.float64)
    lhs = clubb_api.wp2_term_splat_lhs(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        c_wp2_splat=c_wp2_splat, brunt_vaisala_freq_sqd_splat=bv_sqd
    )
    np.testing.assert_allclose(lhs, 1.0)


def test_wp3_term_splat_lhs_constant_bv(gr):
    """wp3 splat lhs should match 1.5 * C_wp2_splat * sqrt(bv_sqd) for constants."""
    c_wp2_splat = np.full(gr.ngrdcol, 0.5, dtype=np.float64)
    bv_sqd = np.full((gr.ngrdcol, gr.nzm), 4.0, dtype=np.float64)
    lhs = clubb_api.wp3_term_splat_lhs(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        c_wp2_splat=c_wp2_splat, brunt_vaisala_freq_sqd_splat=bv_sqd
    )
    assert lhs.shape == (gr.ngrdcol, gr.nzt)
    np.testing.assert_allclose(lhs, 1.5)

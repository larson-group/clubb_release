"""Test wrappers for remapping_module call-tree routines."""
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


def _rho_linear_spline(gr):
    """Return a simple constant-density linear spline covering the full grid."""
    rho_vals = np.array([[1.0, 1.0]], dtype=np.float64)
    rho_levels = np.array([[gr.zm[0, 0], gr.zm[0, -1]]], dtype=np.float64)
    return rho_vals, rho_levels


def test_remap_vals_to_target_identity_zt(gr):
    """Remapping on the same grid should preserve zt-grid data."""
    source = np.linspace(1.0, 2.0, gr.nzt, dtype=np.float64).reshape(1, -1)
    rho_vals, rho_levels = _rho_linear_spline(gr)
    p_sfc = np.array([100000.0], dtype=np.float64)

    target = clubb_api.remap_vals_to_target(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        source_values_idx=source.shape[1],
        target_values_idx=source.shape[1],
        total_idx_rho_lin_spline=rho_vals.shape[1],
        source_values=source, rho_lin_spline_vals=rho_vals, rho_lin_spline_levels=rho_levels,
        iv=1, p_sfc=p_sfc, grid_remap_method=1, l_zt_variable=True,
    )
    assert target.shape == source.shape
    np.testing.assert_allclose(target, source, rtol=1.0e-10, atol=1.0e-12)


def test_remap_vals_to_target_identity_zm(gr):
    """Remapping on the same grid should preserve zm-grid data."""
    source = np.linspace(0.2, 1.2, gr.nzm, dtype=np.float64).reshape(1, -1)
    rho_vals, rho_levels = _rho_linear_spline(gr)
    p_sfc = np.array([100000.0], dtype=np.float64)

    target = clubb_api.remap_vals_to_target(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        source_values_idx=source.shape[1],
        target_values_idx=source.shape[1],
        total_idx_rho_lin_spline=rho_vals.shape[1],
        source_values=source, rho_lin_spline_vals=rho_vals, rho_lin_spline_levels=rho_levels,
        iv=0, p_sfc=p_sfc, grid_remap_method=1, l_zt_variable=False,
    )
    assert target.shape == source.shape
    np.testing.assert_allclose(target, source, rtol=1.0e-10, atol=1.0e-12)

"""Test wrappers for sponge-layer damping call-tree routines."""
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


def test_sponge_damp_xm_simple_midpoint_case(gr):
    """With dt=tau and full-depth sponge, xm should move to midpoint with xm_ref."""
    xm_ref = np.full((gr.ngrdcol, gr.nzt), 10.0, dtype=np.float64)
    xm = np.full((gr.ngrdcol, gr.nzt), 2.0, dtype=np.float64)
    tau = np.ones((gr.ngrdcol, gr.nzt), dtype=np.float64)
    depth = np.full(gr.ngrdcol, 1.0e9, dtype=np.float64)

    damped = clubb_api.sponge_damp_xm(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        dt=1.0, zt=gr.zt, zm=gr.zm, xm_ref=xm_ref, xm=xm,
        tau_sponge_damp=tau, sponge_layer_depth=depth
    )
    np.testing.assert_allclose(damped, 6.0)


def test_sponge_damp_xp2_respects_floor(gr):
    """With dt=tau and full-depth sponge, xp2 should clip to x_tol_sqd floor."""
    xp2 = np.full((gr.ngrdcol, gr.nzm), 9.0, dtype=np.float64)
    tau = np.ones((gr.ngrdcol, gr.nzm), dtype=np.float64)
    depth = np.full(gr.ngrdcol, 1.0e9, dtype=np.float64)
    x_tol_sqd = np.full(gr.ngrdcol, 0.04, dtype=np.float64)

    damped = clubb_api.sponge_damp_xp2(
        gr, nzm=gr.nzm, ngrdcol=gr.ngrdcol,
        dt=1.0, zm=gr.zm, xp2=xp2, x_tol_sqd=x_tol_sqd,
        tau_sponge_damp=tau, sponge_layer_depth=depth,
    )
    np.testing.assert_allclose(damped, 0.04)


def test_sponge_damp_xp3_full_damping_to_zero(gr):
    """With dt=tau and full-depth sponge, xp3 damping factor should be zero."""
    xp3 = np.full((gr.ngrdcol, gr.nzt), 5.0, dtype=np.float64)
    tau = np.ones((gr.ngrdcol, gr.nzt), dtype=np.float64)
    depth = np.full(gr.ngrdcol, 1.0e9, dtype=np.float64)

    damped = clubb_api.sponge_damp_xp3(
        gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol,
        dt=1.0, z=gr.zt, zm=gr.zm, xp3=xp3,
        tau_sponge_damp=tau, sponge_layer_depth=depth,
    )
    np.testing.assert_allclose(damped, 0.0)

"""Test wrappers for routines from CLUBB_core/pos_definite_module.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _make_grid(ngrdcol=1, nzmax=7, dz=100.0):
    """Create a standard ascending grid for positive-definite tests."""
    clubb_api.init_err_info(ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_config_flags(flags)

    zm = np.arange(nzmax, dtype=np.float64) * dz
    zt = 0.5 * (zm[:-1] + zm[1:])

    gr, _err_info = clubb_api.setup_grid(
        nzmax=nzmax,
        ngrdcol=ngrdcol,
        sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
        l_implemented=True,
        l_ascending_grid=True,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        zm_init=np.zeros(ngrdcol, dtype=np.float64),
        zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
        momentum_heights=np.tile(zm, (ngrdcol, 1)),
        thermodynamic_heights=np.tile(zt, (ngrdcol, 1)),
        err_info=ErrInfo(ngrdcol=ngrdcol),
    )
    return gr


def test_pos_definite_adj_preserves_nonnegative_state_with_zero_flux():
    """Zero flux should leave a positive field unchanged and produce zero diagnostics."""
    gr = _make_grid()
    field_n = np.full((gr.ngrdcol, gr.nzt), 0.2, dtype=np.float64)
    field_np1 = field_n.copy()
    flux_np1 = np.zeros((gr.ngrdcol, gr.nzm), dtype=np.float64)

    field_np1_out, flux_np1_out, field_pd, flux_pd = clubb_api.pos_definite_adj(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        dt=10.0,
        field_np1=field_np1,
        flux_np1=flux_np1,
        field_n=field_n,
    )

    np.testing.assert_allclose(field_np1_out, field_np1)
    np.testing.assert_allclose(flux_np1_out, flux_np1)
    np.testing.assert_allclose(field_pd, 0.0)
    np.testing.assert_allclose(flux_pd, 0.0)

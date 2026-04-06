"""Test wrappers for routines from CLUBB_core/turbulent_adv_pdf.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _make_grid(ngrdcol=1, nzmax=7, dz=100.0):
    """Create a standard ascending grid for turbulent-advection tests."""
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


def test_turbulent_adv_pdf_zero_inputs_return_zero_outputs():
    """All turbulent-advection helper wrappers should return zeros for zero inputs."""
    gr = _make_grid()
    shape_zt = (gr.ngrdcol, gr.nzt)
    shape_zm = (gr.ngrdcol, gr.nzm)

    coef_zt = np.zeros(shape_zt, dtype=np.float64)
    rho_ds_zt = np.ones(shape_zt, dtype=np.float64)
    rho_ds_zm = np.ones(shape_zm, dtype=np.float64)
    invrs_rho_ds_zm = np.ones(shape_zm, dtype=np.float64)
    sgn_zm = np.zeros(shape_zm, dtype=np.float64)
    sgn_zt = np.zeros(shape_zt, dtype=np.float64)
    coef_zm = np.zeros(shape_zm, dtype=np.float64)
    term_zt = np.zeros(shape_zt, dtype=np.float64)
    term_zm = np.zeros(shape_zm, dtype=np.float64)

    lhs = clubb_api.xpyp_term_ta_pdf_lhs(
        gr, gr.nzm, gr.nzt, gr.ngrdcol, coef_zt, rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, False, sgn_zm, coef_zm
    )
    lhs_god = clubb_api.xpyp_term_ta_pdf_lhs_godunov(gr, gr.nzm, gr.nzt, gr.ngrdcol, coef_zt, invrs_rho_ds_zm, rho_ds_zm)
    rhs = clubb_api.xpyp_term_ta_pdf_rhs(
        gr, gr.nzm, gr.nzt, gr.ngrdcol, term_zt, rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, False, sgn_zm, term_zm
    )
    rhs_god = clubb_api.xpyp_term_ta_pdf_rhs_godunov(
        gr, gr.nzm, gr.nzt, gr.ngrdcol, term_zm, invrs_rho_ds_zm, sgn_zt, rho_ds_zm
    )

    np.testing.assert_allclose(lhs, 0.0)
    np.testing.assert_allclose(lhs_god, 0.0)
    np.testing.assert_allclose(rhs, 0.0)
    np.testing.assert_allclose(rhs_god, 0.0)

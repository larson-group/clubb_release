"""User-facing wrappers for routines from CLUBB_core/mono_flux_limiter.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def calc_turb_adv_range(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int,
    dt: float, w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm,
):
    set_fortran_grid(gr)
    low, high = clubb_f2py.f2py_calc_turb_adv_range(
        int(nzt), float(dt),
        f_arr(w_1_zm), f_arr(w_2_zm), f_arr(varnce_w_1_zm), f_arr(varnce_w_2_zm), f_arr(mixt_frac_zm),
        nzm=int(nzm), ngrdcol=int(ngrdcol),
    )

    """Return low and high indices with Python 0-based indexing."""
    return np.asarray(low, dtype=np.int32) - 1, np.asarray(high, dtype=np.int32) - 1


def monotonic_turbulent_flux_limit(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, solve_type: int, dt: float,
    xm_old, xp2, wm_zt, xm_forcing, rho_ds_zm, rho_ds_zt,
    invrs_rho_ds_zm, invrs_rho_ds_zt, xp2_threshold: float, xm_tol: float,
    l_implemented: bool, low_lev_effect, high_lev_effect, tridiag_solve_method: int,
    l_upwind_xm_ma: bool, l_mono_flux_lim_spikefix: bool, xm, wpxp,
):
    """Apply the monotonic turbulent flux limiter to xm and wpxp."""
    set_fortran_grid(gr)
    xm_out, wpxp_out = clubb_f2py.f2py_monotonic_turbulent_flux_limit(
        int(solve_type),
        float(dt),
        f_arr(xm_old),
        f_arr(xp2),
        f_arr(wm_zt),
        f_arr(xm_forcing),
        f_arr(rho_ds_zm),
        f_arr(rho_ds_zt),
        f_arr(invrs_rho_ds_zm),
        f_arr(invrs_rho_ds_zt),
        float(xp2_threshold),
        float(xm_tol),
        bool(l_implemented),
        f_arr(np.asarray(low_lev_effect, dtype=np.int32) + 1),
        f_arr(np.asarray(high_lev_effect, dtype=np.int32) + 1),
        int(tridiag_solve_method),
        bool(l_upwind_xm_ma),
        bool(l_mono_flux_lim_spikefix),
        f_arr(xm),
        f_arr(wpxp),
        nzm=int(nzm),
        nzt=int(nzt),
        ngrdcol=int(ngrdcol),
    )
    return xm_out, wpxp_out

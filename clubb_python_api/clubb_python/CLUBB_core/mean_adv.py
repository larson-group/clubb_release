"""User-facing wrappers for routines from CLUBB_core/mean_adv.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def term_ma_zt_lhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    wm_zt,
    weights_zt2zm,
    invrs_dzt,
    invrs_dzm,
    l_upwind_xm_ma: bool,
    grid_dir: float,
):
    """Mean-advection lhs contribution for zt-grid variables."""
    return clubb_f2py.f2py_term_ma_zt_lhs(
        f_arr(wm_zt), f_arr(weights_zt2zm), f_arr(invrs_dzt), f_arr(invrs_dzm),
        l_upwind_xm_ma, float(grid_dir),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))


def term_ma_zm_lhs(nzm: int, nzt: int, ngrdcol: int, wm_zm, invrs_dzm, weights_zm2zt):
    """Mean-advection lhs contribution for zm-grid variables."""
    return clubb_f2py.f2py_term_ma_zm_lhs(
        f_arr(wm_zm), f_arr(invrs_dzm), f_arr(weights_zm2zt),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))

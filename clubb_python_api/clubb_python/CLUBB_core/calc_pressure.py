"""User-facing wrappers for routines from CLUBB_core/calc_pressure.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def init_pressure(gr: Grid, ngrdcol: int, nzt: int, nzm: int, thvm: np.ndarray, p_sfc: np.ndarray):
    """Compute hydrostatic pressure profile."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_init_pressure(
        int(nzm), f_arr(thvm), f_arr(p_sfc),
        ngrdcol=int(ngrdcol), nzt=int(nzt))


def calculate_thvm(nzt: int, ngrdcol: int, thlm, rtm, rcm, exner, thv_ds_zt):
    """Calculate virtual potential temperature."""
    return clubb_f2py.f2py_calculate_thvm(
        f_arr(thlm), f_arr(rtm), f_arr(rcm), f_arr(exner), f_arr(thv_ds_zt),
        nzt=int(nzt), ngrdcol=int(ngrdcol))


def hydrostatic(gr: Grid, ngrdcol: int, nzt: int, nzm: int, thvm, p_sfc):
    """Compute full hydrostatic pressure, exner, and density profiles."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_hydrostatic(
        int(nzm), f_arr(thvm), f_arr(p_sfc), ngrdcol=int(ngrdcol), nzt=int(nzt))

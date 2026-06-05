"""User-facing wrappers for routines from CLUBB_core/calc_pressure.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def init_pressure(ngrdcol: int, gr: Grid, thvm: np.ndarray, p_sfc: np.ndarray, **compat_kwargs):
    """Compute hydrostatic pressure profile."""
    compat_kwargs.pop("nzt", None)
    compat_kwargs.pop("nzm", None)
    set_fortran_grid(gr)
    return clubb_f2py.f2py_init_pressure(
        int(gr.nzm), f_arr(thvm), f_arr(p_sfc),
        ngrdcol=int(ngrdcol), nzt=int(gr.nzt))


def calculate_thvm(nzt: int, ngrdcol: int, thlm, rtm, rcm, exner, thv_ds_zt):
    """Calculate virtual potential temperature."""
    return clubb_f2py.f2py_calculate_thvm(
        f_arr(thlm), f_arr(rtm), f_arr(rcm), f_arr(exner), f_arr(thv_ds_zt),
        nzt=int(nzt), ngrdcol=int(ngrdcol))


def hydrostatic(ngrdcol: int, gr: Grid, thvm, p_sfc, **compat_kwargs):
    """Compute full hydrostatic pressure, exner, and density profiles."""
    compat_kwargs.pop("nzt", None)
    compat_kwargs.pop("nzm", None)
    set_fortran_grid(gr)
    return clubb_f2py.f2py_hydrostatic(
        int(gr.nzm), f_arr(thvm), f_arr(p_sfc), ngrdcol=int(ngrdcol), nzt=int(gr.nzt))

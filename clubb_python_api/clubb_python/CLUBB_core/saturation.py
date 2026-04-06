"""User-facing wrappers for routines from CLUBB_core/saturation.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def sat_mixrat_liq(gr: Grid, nz: int, ngrdcol: int, p_in_Pa, T_in_K,
                   saturation_formula: int):
    """Compute saturation mixing ratio for liquid water."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_sat_mixrat_liq_2d(
        p_in_Pa=f_arr(p_in_Pa),
        t_in_k=f_arr(T_in_K),
        saturation_formula=int(saturation_formula),
        nz=int(nz),
        ngrdcol=int(ngrdcol),
    )


def sat_mixrat_ice(nz: int, ngrdcol: int, p_in_Pa, T_in_K,
                   saturation_formula: int):
    """Compute saturation mixing ratio for ice."""
    return clubb_f2py.f2py_sat_mixrat_ice_2d(
        p_in_Pa=f_arr(p_in_Pa),
        t_in_k=f_arr(T_in_K),
        saturation_formula=int(saturation_formula),
        nz=int(nz),
        ngrdcol=int(ngrdcol),
    )


def rcm_sat_adj(thlm: float, rtm: float, p_in_Pa: float, exner: float,
                saturation_formula: int) -> float:
    """Compute cloud water mixing ratio via saturation adjustment (scalar)."""
    return float(clubb_f2py.f2py_rcm_sat_adj(
        float(thlm), float(rtm), float(p_in_Pa), float(exner), int(saturation_formula)))

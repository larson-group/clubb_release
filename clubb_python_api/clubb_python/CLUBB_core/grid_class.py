"""User-facing wrappers for routines from CLUBB_core/grid_class.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import get_fortran_grid, set_fortran_grid
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def setup_grid(
    nzmax: int,
    ngrdcol: int,
    sfc_elevation: np.ndarray,
    l_implemented: bool,
    l_ascending_grid: bool,
    grid_type: int,
    deltaz: np.ndarray,
    zm_init: np.ndarray,
    zm_top: np.ndarray,
    momentum_heights: np.ndarray,
    thermodynamic_heights: np.ndarray,
    err_info: ErrInfo,
):
    """Set up the CLUBB vertical grid."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_setup_grid(
        f_arr(sfc_elevation),
        l_implemented,
        l_ascending_grid,
        grid_type,
        f_arr(deltaz),
        f_arr(zm_init),
        f_arr(zm_top),
        f_arr(momentum_heights),
        f_arr(thermodynamic_heights),
        nzmax=int(nzmax),
        ngrdcol=int(ngrdcol),
    )
    return get_fortran_grid(), get_fortran_err_info()


def zm2zt2zm(gr: Grid, nzm: int, nzt: int, ngrdcol: int, azm, zm_min: float = 0.0):
    """Smooth a momentum-level field by mapping zm->zt->zm."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_zm2zt2zm_2d(
        int(nzt), f_arr(azm), zm_min, nzm=int(nzm), ngrdcol=int(ngrdcol))


def zt2zm(gr: Grid, nzm: int, nzt: int, ngrdcol: int, azt):
    """Interpolate a 2D field from zt (thermo) to zm (momentum) levels."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_zt2zm_2d(
        int(nzm), f_arr(azt), nzt=int(nzt), ngrdcol=int(ngrdcol))


def zm2zt(gr: Grid, nzm: int, nzt: int, ngrdcol: int, azm):
    """Interpolate a 2D field from zm (momentum) to zt (thermo) levels."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_zm2zt_2d(
        int(nzt), f_arr(azm), nzm=int(nzm), ngrdcol=int(ngrdcol))


def zt2zm2zt(gr: Grid, nzm: int, nzt: int, ngrdcol: int, azt, zt_min: float = 0.0):
    """Smooth a thermo-level field by mapping zt->zm->zt."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_zt2zm2zt_2d(
        int(nzm), f_arr(azt), zt_min, nzt=int(nzt), ngrdcol=int(ngrdcol))


def ddzm(gr: Grid, nzm: int, nzt: int, ngrdcol: int, azm):
    """Differentiate a zm-grid field across zt levels, returning a zt-grid field."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_ddzm_2d(
        int(nzt), f_arr(azm), nzm=int(nzm), ngrdcol=int(ngrdcol))


def ddzt(gr: Grid, nzm: int, nzt: int, ngrdcol: int, azt):
    """Differentiate a zt-grid field across zm levels, returning a zm-grid field."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_ddzt_2d(
        int(nzm), f_arr(azt), nzt=int(nzt), ngrdcol=int(ngrdcol))


def cleanup_grid(gr: Grid):
    """Deallocate stored grid arrays."""
    set_fortran_grid(gr)
    clubb_f2py.f2py_cleanup_grid()


def setup_grid_heights(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int,
    l_implemented: bool, l_ascending_grid: bool, grid_type: int,
    deltaz, zm_init, momentum_heights, thermodynamic_heights,
    err_info: ErrInfo,
):
    """Update a grid's height fields and interpolation weights."""
    set_fortran_grid(gr)
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_setup_grid_heights(
        l_implemented, l_ascending_grid, int(grid_type),
        f_arr(deltaz), f_arr(zm_init), f_arr(momentum_heights), f_arr(thermodynamic_heights),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return get_fortran_grid(), get_fortran_err_info()


def read_grid_heights(
    nzmax: int, grid_type: int, zm_grid_fname: str, zt_grid_fname: str,
    file_unit: int, err_info: ErrInfo,
):
    """Read stretched-grid heights from file and return them with updated err_info."""
    set_fortran_err_info(err_info)
    momentum_heights, thermodynamic_heights = clubb_f2py.f2py_read_grid_heights(
        int(nzmax), int(grid_type), str(zm_grid_fname), str(zt_grid_fname), int(file_unit)
    )
    return get_fortran_err_info(), momentum_heights, thermodynamic_heights

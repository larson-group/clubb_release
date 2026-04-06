"""User-facing wrappers for routines from CLUBB_core/sponge_layer_damping.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def sponge_damp_xm(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int,
    dt: float, zt, zm, xm_ref, xm, tau_sponge_damp, sponge_layer_depth,
):
    """Damp a mean-field profile toward a reference profile in sponge layer."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_sponge_damp_xm(
        float(dt), f_arr(zt), f_arr(zm), f_arr(xm_ref), f_arr(xm),
        f_arr(tau_sponge_damp), f_arr(sponge_layer_depth),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))


def sponge_damp_xp2(
    gr: Grid, nzm: int, ngrdcol: int,
    dt: float, zm, xp2, x_tol_sqd, tau_sponge_damp, sponge_layer_depth,
):
    """Damp variance profile in sponge layer with a lower floor."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_sponge_damp_xp2(
        float(dt), f_arr(zm), f_arr(xp2), f_arr(x_tol_sqd),
        f_arr(tau_sponge_damp), f_arr(sponge_layer_depth),
        nzm=int(nzm), ngrdcol=int(ngrdcol))


def sponge_damp_xp3(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int,
    dt: float, z, zm, xp3, tau_sponge_damp, sponge_layer_depth,
):
    """Damp third-moment profile in sponge layer."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_sponge_damp_xp3(
        float(dt), f_arr(z), f_arr(zm), f_arr(xp3),
        f_arr(tau_sponge_damp), f_arr(sponge_layer_depth),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))

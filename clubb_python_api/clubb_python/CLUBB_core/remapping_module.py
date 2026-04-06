"""User-facing wrappers for routines from CLUBB_core/remapping_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def remap_vals_to_target(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, source_values_idx: int,
    source_values, target_values_idx: int, total_idx_rho_lin_spline: int,
    rho_lin_spline_vals, rho_lin_spline_levels,
    iv: int, p_sfc, grid_remap_method: int, l_zt_variable: bool,
):
    """Remap values using remapping_module with source and target set to the same grid."""
    set_fortran_grid(gr)
    kwargs = dict(
        source_values=f_arr(source_values),
        target_values_idx=int(target_values_idx),
        rho_lin_spline_vals=f_arr(rho_lin_spline_vals),
        rho_lin_spline_levels=f_arr(rho_lin_spline_levels),
        iv=int(iv),
        p_sfc=f_arr(p_sfc),
        grid_remap_method=int(grid_remap_method),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
        source_values_idx=int(source_values_idx),
        total_idx_rho_lin_spline=int(total_idx_rho_lin_spline),
    )
    try:
        return clubb_f2py.f2py_remap_vals_to_target_same_grid(
            l_zt_variable=l_zt_variable, **kwargs
        )
    except TypeError as exc:
        # Keep compatibility with an already-built extension exposing the old keyword.
        if "il_zt_variable" not in str(exc):
            raise
        return clubb_f2py.f2py_remap_vals_to_target_same_grid(
            il_zt_variable=l_zt_variable, **kwargs
        )

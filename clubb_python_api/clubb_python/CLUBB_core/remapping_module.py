"""User-facing wrappers for routines from CLUBB_core/remapping_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def remap_vals_to_target(
    ngrdcol: int, gr_source: Grid | None = None, gr_target: Grid | None = None, source_values_idx: int | None = None,
    source_values=None, target_values_idx: int | None = None, total_idx_rho_lin_spline: int | None = None,
    rho_lin_spline_vals=None, rho_lin_spline_levels=None,
    iv: int = 0, p_sfc=None, grid_remap_method: int = 1, l_zt_variable: bool = False,
    **compat_kwargs,
):
    """Remap values using remapping_module with source and target set to the same grid."""
    same_grid = compat_kwargs.pop("gr", None)
    if gr_source is None:
        gr_source = same_grid
    if gr_target is None:
        gr_target = same_grid if same_grid is not None else gr_source
    if gr_source is None or gr_target is None:
        raise ValueError("remap_vals_to_target requires gr_source/gr_target or legacy gr.")
    nzm = int(compat_kwargs.pop("nzm", gr_source.nzm))
    nzt = int(compat_kwargs.pop("nzt", gr_source.nzt))
    set_fortran_grid(gr_source)
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

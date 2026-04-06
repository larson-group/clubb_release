"""User-facing wrappers for routines from CLUBB_core/pos_definite_module.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def pos_definite_adj(gr: Grid, nzm: int, nzt: int, ngrdcol: int, dt: float, field_np1, flux_np1, field_n):
    """Apply the positive-definite flux limiter using the active stored grid."""
    set_fortran_grid(gr)
    field_np1_out, flux_np1_out, field_pd, flux_pd = clubb_f2py.f2py_pos_definite_adj(
        float(dt), f_arr(field_np1), f_arr(flux_np1), f_arr(field_n), nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol)
    )
    return field_np1_out, flux_np1_out, field_pd, flux_pd

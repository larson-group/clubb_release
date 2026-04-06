"""User-facing wrappers for routines from CLUBB_core/penta_bicgstab_solver.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py


def penta_bicgstab_solve(ndim: int, ngrdcol: int, lhs, rhs):
    """Solve pentadiagonal systems with CLUBB's BiCGSTAB solver."""
    return clubb_f2py.f2py_penta_bicgstab_solve(f_arr(lhs), f_arr(rhs), ndim=int(ndim), ngrdcol=int(ngrdcol))

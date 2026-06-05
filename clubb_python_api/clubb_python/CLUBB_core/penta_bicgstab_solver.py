"""User-facing wrappers for routines from CLUBB_core/penta_bicgstab_solver.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def penta_bicgstab_solve(ndim: int, ngrdcol: int, lhs, rhs, sol_init=None):
    """Solve pentadiagonal systems with CLUBB's BiCGSTAB solver."""
    if sol_init is None:
        sol_init = np.zeros_like(rhs)
    try:
        return clubb_f2py.f2py_penta_bicgstab_solve(
            f_arr(lhs), f_arr(rhs), f_arr(sol_init), ndim=int(ndim), ngrdcol=int(ngrdcol)
        )
    except TypeError:
        return clubb_f2py.f2py_penta_bicgstab_solve(
            f_arr(lhs), f_arr(rhs), ndim=int(ndim), ngrdcol=int(ngrdcol)
        )

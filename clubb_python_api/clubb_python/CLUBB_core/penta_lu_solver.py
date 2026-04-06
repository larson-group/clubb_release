"""User-facing wrappers for routines from CLUBB_core/penta_lu_solver.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def penta_lu_solve(ndim: int, ngrdcol: int, lhs, rhs):
    """Solve pentadiagonal systems using CLUBB's LU solver."""
    rhs_arr = np.asarray(rhs, dtype=np.float64)
    lhs_f = f_arr(lhs, dtype=np.float64)

    if rhs_arr.ndim == 2:
        return clubb_f2py.f2py_penta_lu_solve_single_rhs_multiple_lhs(
            lhs_f, f_arr(rhs_arr, dtype=np.float64), ndim=int(ndim), ngrdcol=int(ngrdcol)
        )

    if rhs_arr.ndim == 3:
        return clubb_f2py.f2py_penta_lu_solve_multiple_rhs_lhs(
            lhs_f, f_arr(rhs_arr, dtype=np.float64),
            ndim=int(ndim), nrhs=int(rhs_arr.shape[2]), ngrdcol=int(ngrdcol)
        )

    raise ValueError("rhs must be rank-2 or rank-3 for penta_lu_solve")

"""User-facing wrappers for routines from CLUBB_core/tridiag_lu_solver.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def tridiag_lu_solve(ndim: int, lhs, rhs):
    """Solve tridiagonal systems using CLUBB's direct LU solver."""
    lhs_arr = np.asarray(lhs, dtype=np.float64)
    rhs_arr = np.asarray(rhs, dtype=np.float64)

    if lhs_arr.ndim == 2 and rhs_arr.ndim == 1:
        return clubb_f2py.f2py_tridiag_lu_solve_single_rhs_lhs(
            f_arr(lhs_arr, dtype=np.float64), f_arr(rhs_arr, dtype=np.float64), ndim=int(ndim)
        )

    if lhs_arr.ndim == 3 and rhs_arr.ndim == 2:
        return clubb_f2py.f2py_tridiag_lu_solve_single_rhs_multiple_lhs(
            f_arr(lhs_arr, dtype=np.float64), f_arr(rhs_arr, dtype=np.float64),
            ndim=int(ndim), ngrdcol=int(lhs_arr.shape[1])
        )

    if lhs_arr.ndim == 3 and rhs_arr.ndim == 3:
        return clubb_f2py.f2py_tridiag_lu_solve_multiple_rhs_lhs(
            f_arr(lhs_arr, dtype=np.float64), f_arr(rhs_arr, dtype=np.float64),
            ndim=int(ndim), nrhs=int(rhs_arr.shape[2]), ngrdcol=int(lhs_arr.shape[1])
        )

    raise ValueError("Unsupported lhs/rhs ranks for tridiag_lu_solve")

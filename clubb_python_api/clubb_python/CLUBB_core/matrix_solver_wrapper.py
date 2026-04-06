"""User-facing wrappers for routines from CLUBB_core/matrix_solver_wrapper.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def _fa64(arr):
    """Fortran-order float64 conversion for solver buffers."""
    return f_arr(arr, dtype=np.float64)


def band_solve(solve_name: str, penta_solve_method: int, ngrdcol: int, nsup: int,
               nsub: int, ndim: int, nrhs: int, lhs, rhs, err_info: ErrInfo,
               old_soln=None, use_rcond: bool = False):
    """Solve a banded system with multiple right-hand sides."""
    lhs_f = _fa64(lhs)
    rhs_f = _fa64(rhs)
    if old_soln is None:
        old_soln_f = np.zeros((ngrdcol, ndim, nrhs), dtype=np.float64, order="F")
    else:
        old_soln_f = _fa64(old_soln)

    set_fortran_err_info(err_info)
    _, _, soln, rcond = clubb_f2py.f2py_band_solve_multiple_rhs(
        solve_name=str(solve_name),
        penta_solve_method=int(penta_solve_method),
        lhs=lhs_f,
        rhs=rhs_f,
        old_soln=old_soln_f,
        use_rcond=bool(use_rcond),
        ngrdcol=int(ngrdcol),
        nsup=int(nsup),
        nsub=int(nsub),
        ndim=int(ndim),
        nrhs=int(nrhs),
    )
    return get_fortran_err_info(), soln, rcond


def tridiag_solve(solve_name: str, tridiag_solve_method: int, ngrdcol: int,
                  ndim: int, lhs, rhs, err_info: ErrInfo, use_rcond: bool = False):
    """Solve a tridiagonal system for either single or multiple RHS inputs."""
    lhs_f = _fa64(lhs)
    set_fortran_err_info(err_info)

    rhs_arr = np.asarray(rhs, dtype=np.float64)
    if rhs_arr.ndim == 2:
        _, _, soln, rcond = clubb_f2py.f2py_tridiag_solve_single_rhs_multiple_lhs(
            solve_name=str(solve_name),
            tridiag_solve_method=int(tridiag_solve_method),
            lhs=lhs_f,
            rhs=_fa64(rhs_arr),
            use_rcond=bool(use_rcond),
            ngrdcol=int(ngrdcol),
            ndim=int(ndim),
        )
        return get_fortran_err_info(), soln, rcond

    if rhs_arr.ndim == 3:
        _, _, soln, rcond = clubb_f2py.f2py_tridiag_solve_multiple_rhs(
            solve_name=str(solve_name),
            tridiag_solve_method=int(tridiag_solve_method),
            lhs=lhs_f,
            rhs=_fa64(rhs_arr),
            use_rcond=bool(use_rcond),
            ngrdcol=int(ngrdcol),
            ndim=int(ndim),
            nrhs=int(rhs_arr.shape[2]),
        )
        return get_fortran_err_info(), soln, rcond

    raise ValueError("rhs must be rank-2 or rank-3 for tridiag_solve")

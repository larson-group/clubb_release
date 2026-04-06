"""User-facing wrappers for routines from CLUBB_core/lapack_wrap.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def lapack_tridiag_solve(solve_type: str, ndim: int, nrhs: int, ngrdcol: int, lhs, rhs, err_info: ErrInfo):
    """Solve tridiagonal systems using CLUBB's LAPACK wrapper."""
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_lapack_tridiag_solve(
        str(solve_type), f_arr(lhs), f_arr(rhs), ndim=int(ndim), nrhs=int(nrhs), ngrdcol=int(ngrdcol)
    )
    return get_fortran_err_info(), result


def lapack_tridiag_solvex(solve_type: str, ndim: int, nrhs: int, ngrdcol: int, lhs, rhs, err_info: ErrInfo):
    """Solve tridiagonal systems with CLUBB's expert LAPACK wrapper."""
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_lapack_tridiag_solvex(
        str(solve_type), f_arr(lhs), f_arr(rhs), ndim=int(ndim), nrhs=int(nrhs), ngrdcol=int(ngrdcol)
    )
    return get_fortran_err_info(), *result


def lapack_band_solve(solve_type: str, nsup: int, nsub: int, ndim: int, nrhs: int, ngrdcol: int,
                      lhs, rhs, err_info: ErrInfo):
    """Solve banded systems using CLUBB's simple LAPACK wrapper."""
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_lapack_band_solve(
        str(solve_type), int(nsup), int(nsub), f_arr(lhs), f_arr(rhs),
        ndim=int(ndim), nrhs=int(nrhs), ngrdcol=int(ngrdcol)
    )
    return get_fortran_err_info(), result


def lapack_band_solvex(solve_type: str, nsup: int, nsub: int, ndim: int, nrhs: int, ngrdcol: int,
                       lhs, rhs, err_info: ErrInfo):
    """Solve banded systems using CLUBB's expert LAPACK wrapper."""
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_lapack_band_solvex(
        str(solve_type), int(nsup), int(nsub), f_arr(lhs), f_arr(rhs),
        ndim=int(ndim), nrhs=int(nrhs), ngrdcol=int(ngrdcol)
    )
    return get_fortran_err_info(), *result

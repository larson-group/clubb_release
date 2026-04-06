"""Test wrappers for routines from CLUBB_core/lapack_wrap.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def test_lapack_tridiag_wrappers_solve_diagonal_system():
    """Simple and expert tridiagonal CLUBB LAPACK wrappers should solve a diagonal system."""
    ngrdcol = 2
    ndim = 4
    nrhs = 1
    lhs = np.zeros((3, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[1, :, :] = np.array([[2.0, 4.0, 8.0, 16.0], [3.0, 6.0, 9.0, 12.0]], dtype=np.float64)
    expected = np.array(
        [[[1.0], [2.0], [3.0], [4.0]], [[5.0], [6.0], [7.0], [8.0]]],
        dtype=np.float64,
        order="F",
    )
    rhs = (lhs[1, :, :][:, :, None] * expected).copy(order="F")
    err_info = clubb_api.init_err_info(ErrInfo(ngrdcol=ngrdcol))

    err_simple, soln_simple = clubb_api.lapack_tridiag_solve("lapack_tridiag_test", ndim, nrhs, ngrdcol, lhs, rhs, err_info)
    err_expert, soln_expert, rcond = clubb_api.lapack_tridiag_solvex("lapack_tridiag_test", ndim, nrhs, ngrdcol, lhs, rhs, err_info)

    np.testing.assert_allclose(soln_simple, expected)
    np.testing.assert_allclose(soln_expert, expected)
    assert rcond.shape == (ngrdcol,)
    assert isinstance(err_simple, ErrInfo)
    assert isinstance(err_expert, ErrInfo)


def test_lapack_band_wrappers_solve_diagonal_system():
    """Simple and expert band CLUBB LAPACK wrappers should solve a diagonal system."""
    ngrdcol = 2
    ndim = 5
    nrhs = 1
    nsup = 2
    nsub = 2
    lhs = np.zeros((nsup + nsub + 1, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[nsup, :, :] = np.array(
        [[2.0, 4.0, 8.0, 16.0, 32.0], [3.0, 6.0, 9.0, 12.0, 15.0]],
        dtype=np.float64,
    )
    expected = np.array(
        [[[1.0], [2.0], [3.0], [4.0], [5.0]], [[6.0], [7.0], [8.0], [9.0], [10.0]]],
        dtype=np.float64,
        order="F",
    )
    rhs = (lhs[nsup, :, :][:, :, None] * expected).copy(order="F")
    err_info = clubb_api.init_err_info(ErrInfo(ngrdcol=ngrdcol))

    err_simple, soln_simple = clubb_api.lapack_band_solve("lapack_band_test", nsup, nsub, ndim, nrhs, ngrdcol, lhs, rhs, err_info)
    err_expert, soln_expert, rcond = clubb_api.lapack_band_solvex("lapack_band_test", nsup, nsub, ndim, nrhs, ngrdcol, lhs, rhs, err_info)

    np.testing.assert_allclose(soln_simple, expected)
    np.testing.assert_allclose(soln_expert, expected)
    assert rcond.shape == (ngrdcol,)
    assert isinstance(err_simple, ErrInfo)
    assert isinstance(err_expert, ErrInfo)

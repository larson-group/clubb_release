"""Test wrappers for matrix solver call-tree routines."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def test_band_solve_multiple_rhs_solves_diagonal_system():
    """band_solve should solve a simple diagonal multiple-RHS system."""
    ngrdcol = 2
    ndim = 6
    nrhs = 2
    nsup = 2
    nsub = 2
    diag = np.array(
        [[2.0, 4.0, 8.0, 16.0, 32.0, 64.0],
         [3.0, 6.0, 9.0, 12.0, 15.0, 18.0]],
        dtype=np.float64,
        order="F",
    )
    expected = np.array(
        [[[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0], [5.0, 50.0], [6.0, 60.0]],
         [[7.0, 70.0], [8.0, 80.0], [9.0, 90.0], [10.0, 100.0], [11.0, 110.0], [12.0, 120.0]]],
        dtype=np.float64,
        order="F",
    )
    lhs = np.zeros((nsup + nsub + 1, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[nsup, :, :] = diag
    rhs = (diag[:, :, np.newaxis] * expected).copy(order="F")
    err_info = ErrInfo(ngrdcol=ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_err_info(err_info)

    err_info_out, soln, rcond = clubb_api.band_solve(
        solve_name="test_band",
        penta_solve_method=flags.penta_solve_method,
        ngrdcol=ngrdcol,
        nsup=nsup,
        nsub=nsub,
        ndim=ndim,
        nrhs=nrhs,
        lhs=lhs,
        rhs=rhs,
        err_info=err_info,
    )

    assert soln.shape == expected.shape
    assert np.allclose(soln, expected)
    assert rcond.shape == (ngrdcol,)
    assert isinstance(err_info_out, ErrInfo)


def test_tridiag_solve_single_rhs_multiple_lhs_solves_diagonal_system():
    """tridiag_solve should solve a simple diagonal single-RHS system."""
    ngrdcol = 2
    ndim = 4
    diag = np.array(
        [[2.0, 4.0, 8.0, 16.0],
         [3.0, 6.0, 9.0, 12.0]],
        dtype=np.float64,
        order="F",
    )
    expected = np.array(
        [[1.0, 2.0, 3.0, 4.0],
         [5.0, 6.0, 7.0, 8.0]],
        dtype=np.float64,
        order="F",
    )
    lhs = np.zeros((3, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[1, :, :] = diag
    rhs = (diag * expected).copy(order="F")
    err_info = ErrInfo(ngrdcol=ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_err_info(err_info)

    err_info_out, soln, rcond = clubb_api.tridiag_solve(
        solve_name="test_tridiag_single",
        tridiag_solve_method=flags.tridiag_solve_method,
        ngrdcol=ngrdcol,
        ndim=ndim,
        lhs=lhs,
        rhs=rhs,
        err_info=err_info,
    )

    assert soln.shape == expected.shape
    assert np.allclose(soln, expected)
    assert rcond.shape == (ngrdcol,)
    assert isinstance(err_info_out, ErrInfo)


def test_tridiag_solve_multiple_rhs_solves_diagonal_system():
    """tridiag_solve should solve a simple diagonal multiple-RHS system."""
    ngrdcol = 2
    ndim = 3
    nrhs = 2
    diag = np.array(
        [[2.0, 4.0, 8.0],
         [3.0, 6.0, 9.0]],
        dtype=np.float64,
        order="F",
    )
    expected = np.array(
        [[[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]],
         [[4.0, 40.0], [5.0, 50.0], [6.0, 60.0]]],
        dtype=np.float64,
        order="F",
    )
    lhs = np.zeros((3, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[1, :, :] = diag
    rhs = (diag[:, :, np.newaxis] * expected).copy(order="F")
    err_info = ErrInfo(ngrdcol=ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_err_info(err_info)

    err_info_out, soln, rcond = clubb_api.tridiag_solve(
        solve_name="test_tridiag_multi",
        tridiag_solve_method=flags.tridiag_solve_method,
        ngrdcol=ngrdcol,
        ndim=ndim,
        lhs=lhs,
        rhs=rhs,
        err_info=err_info,
    )

    assert soln.shape == expected.shape
    assert np.allclose(soln, expected)
    assert rcond.shape == (ngrdcol,)
    assert isinstance(err_info_out, ErrInfo)

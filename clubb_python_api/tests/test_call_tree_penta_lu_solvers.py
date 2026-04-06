"""Test wrappers for routines from CLUBB_core/penta_lu_solver.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_penta_lu_solve_single_rhs_multiple_lhs_solves_diagonal_system():
    """penta_lu_solve should solve a diagonal system with one RHS."""
    ngrdcol = 2
    ndim = 5
    lhs = np.zeros((5, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[2, :, :] = np.array(
        [[2.0, 4.0, 8.0, 16.0, 32.0], [3.0, 6.0, 9.0, 12.0, 15.0]],
        dtype=np.float64,
    )
    expected = np.array([[1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0]], dtype=np.float64, order="F")
    rhs = (lhs[2, :, :] * expected).copy(order="F")

    soln = clubb_api.penta_lu_solve(ndim=ndim, ngrdcol=ngrdcol, lhs=lhs, rhs=rhs)
    np.testing.assert_allclose(soln, expected)


def test_penta_lu_solve_multiple_rhs_lhs_solves_diagonal_system():
    """penta_lu_solve should solve a diagonal system with multiple RHS inputs."""
    ngrdcol = 2
    ndim = 4
    nrhs = 2
    lhs = np.zeros((5, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[2, :, :] = np.array([[2.0, 4.0, 8.0, 16.0], [3.0, 6.0, 9.0, 12.0]], dtype=np.float64)
    expected = np.array(
        [[[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]],
         [[5.0, 50.0], [6.0, 60.0], [7.0, 70.0], [8.0, 80.0]]],
        dtype=np.float64,
        order="F",
    )
    rhs = (lhs[2, :, :][:, :, None] * expected).copy(order="F")

    soln = clubb_api.penta_lu_solve(ndim=ndim, ngrdcol=ngrdcol, lhs=lhs, rhs=rhs)
    np.testing.assert_allclose(soln, expected)

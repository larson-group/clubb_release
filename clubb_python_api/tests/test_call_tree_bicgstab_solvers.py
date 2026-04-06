"""Test wrappers for routines from CLUBB_core/penta_bicgstab_solver.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_penta_bicgstab_solve_solves_diagonal_system():
    """penta_bicgstab_solve should preserve the zero solution for zero RHS."""
    ngrdcol = 2
    ndim = 4
    lhs = np.zeros((5, ngrdcol, ndim), dtype=np.float64, order="F")
    lhs[2, :, :] = np.array([[2.0, 4.0, 8.0, 16.0], [3.0, 6.0, 9.0, 12.0]], dtype=np.float64)
    rhs = np.zeros((ngrdcol, ndim), dtype=np.float64, order="F")

    soln, iters = clubb_api.penta_bicgstab_solve(ndim=ndim, ngrdcol=ngrdcol, lhs=lhs, rhs=rhs)

    np.testing.assert_allclose(soln, rhs)
    assert int(iters) >= 0

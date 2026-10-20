"""Check JAX tridiagonal LU solver against Fortran."""

import jax.numpy as jnp
import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_rel_diff
from clubb_jax.src.CLUBB_core.tridiag_lu_solver import tridiag_lu_solve


def call_tridiag_lu_solve(lhs, rhs):
    return tridiag_lu_solve(rhs.shape[-1], lhs, rhs)


def test_solver_f2py_oracle():
    """f2py bit-shadow vs the Fortran oracle. call_tridiag_lu_solve is a faithful port of
    `tridiag_lu_solve_single_rhs_multiple_lhs` (tridiag_lu_solver.F90) — this is THE workhorse solver every
    prognostic advance routes through, so bit-matching the SPECIFIC Fortran LU (not merely numpy.linalg.solve,
    which only confirms the system is solved) is what underpins the whole bit-faithful case suite."""
    import clubb_f2py
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(20):
        ng, nd = 2, 15
        lhs = np.zeros((3, ng, nd))
        lhs[0] = rng.uniform(-0.3, 0.3, (ng, nd))   # superdiagonal
        lhs[1] = rng.uniform(2.0, 4.0, (ng, nd))    # main (diagonally dominant → stable LU)
        lhs[2] = rng.uniform(-0.3, 0.3, (ng, nd))   # subdiagonal
        rhs = rng.uniform(-1.0, 1.0, (ng, nd))
        j = np.asarray(call_tridiag_lu_solve(jnp.asarray(lhs), jnp.asarray(rhs)))
        f = np.asarray(clubb_f2py.f2py_tridiag_lu_solve_single_rhs_multiple_lhs(lhs, rhs))
        worst = max(worst, max_rel_diff("solver", j, f, 1e-30))
    assert worst < 1e-12, f"tridiag_lu_solve f2py mismatch {worst:.2e}"

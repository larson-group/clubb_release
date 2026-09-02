"""Penta-diagonal LU solver faithfulness regression.

Establishes (Iter144) that ``penta_lu_solve_jax`` is a bit-faithful port of the
Fortran ``penta_lu_solve_single_rhs_multiple_lhs`` operation order:

  * EAGER JAX (jit disabled) reproduces a pure-numpy float64 replica of the
    Fortran elimination EXACTLY (0 ULP) — the algorithm and op-order match.
  * JIT JAX (XLA) differs from the eager/numpy result only by XLA's FMA
    contraction, bounded at ~1e-13 relative for well-conditioned systems.

This rules the solver OUT as a source of cross-case divergence: any case-level
mismatch (e.g. rico) lives in the lhs/rhs *assembly* or its inputs, not the
solve. See DESIGN.md "rico diagnosis".
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.penta_lu_solver import penta_lu_solve


def penta_lu_solve_jax(lhs, rhs):
    return penta_lu_solve(rhs.shape[-1], lhs.shape[1], lhs, rhs)


def _fortran_penta_replica(lhs, rhs):
    """Pure-numpy float64 replica of penta_lu_solve_single_rhs_multiple_lhs.

    Band map (matches penta_lu_solve_jax docstring):
      lhs[0]=super2(-2), lhs[1]=super1(-1), lhs[2]=diag(0), lhs[3]=sub1(1), lhs[4]=sub2(2)
    """
    s2, s1, dg, b1, b2 = lhs[0], lhs[1], lhs[2], lhs[3], lhs[4]
    ng, ndim = rhs.shape
    soln = np.zeros((ng, ndim))
    for i in range(ng):
        ldi = np.zeros(ndim); u1 = np.zeros(ndim); u2 = np.zeros(ndim)
        l1 = np.zeros(ndim); l2 = np.zeros(ndim)
        ldi[0] = 1.0 / dg[i, 0]
        u1[0] = ldi[0] * s1[i, 0]
        u2[0] = ldi[0] * s2[i, 0]
        l1[1] = b1[i, 1]
        ldi[1] = 1.0 / (dg[i, 1] - l1[1] * u1[0])
        u1[1] = ldi[1] * (s1[i, 1] - l1[1] * u2[0])
        u2[1] = ldi[1] * s2[i, 1]
        for k in range(2, ndim):
            l2[k] = b2[i, k]
            l1[k] = b1[i, k] - l2[k] * u1[k - 2]
            ldi[k] = 1.0 / (dg[i, k] - l2[k] * u2[k - 2] - l1[k] * u1[k - 1])
            if k < ndim - 1:
                u1[k] = ldi[k] * (s1[i, k] - l1[k] * u2[k - 1])
            if k < ndim - 2:
                u2[k] = ldi[k] * s2[i, k]
        y = np.zeros(ndim)
        y[0] = ldi[0] * rhs[i, 0]
        y[1] = ldi[1] * (rhs[i, 1] - l1[1] * y[0])
        for k in range(2, ndim):
            y[k] = ldi[k] * (rhs[i, k] - l2[k] * y[k - 2] - l1[k] * y[k - 1])
        x = y.copy()
        x[ndim - 2] = y[ndim - 2] - u1[ndim - 2] * x[ndim - 1]
        for k in range(ndim - 3, -1, -1):
            x[k] = y[k] - u1[k] * x[k + 1] - u2[k] * x[k + 2]
        soln[i] = x
    return soln


def _make_system(ndim, seed):
    """A diagonally-dominant penta system with a large-magnitude RHS (the thl
    regime where catastrophic cancellation in the coupling matters)."""
    rng = np.random.RandomState(seed)
    lhs = np.zeros((5, 1, ndim))
    lhs[2, 0] = 4.0 + rng.rand(ndim)            # diag (dominant)
    lhs[1, 0] = 0.5 * (rng.rand(ndim) - 0.5)    # super1
    lhs[0, 0] = 0.25 * (rng.rand(ndim) - 0.5)   # super2
    lhs[3, 0] = 0.5 * (rng.rand(ndim) - 0.5)    # sub1
    lhs[4, 0] = 0.25 * (rng.rand(ndim) - 0.5)   # sub2
    # boundary bands the solver ignores
    lhs[0, 0, ndim - 2:] = 0.0
    lhs[1, 0, ndim - 1:] = 0.0
    lhs[3, 0, 0] = 0.0
    lhs[4, 0, :2] = 0.0
    rhs = 300.0 + rng.rand(1, ndim)             # ~thlm magnitude
    return lhs, rhs


def test_eager_jax_matches_fortran_replica_exactly():
    """EAGER JAX == numpy Fortran-order replica to 0 ULP (faithful op order)."""
    for seed in range(4):
        lhs, rhs = _make_system(115, seed)
        with jax.disable_jit():
            jsol = np.asarray(penta_lu_solve_jax(jnp.asarray(lhs), jnp.asarray(rhs)))
        fsol = _fortran_penta_replica(lhs, rhs)
        assert np.array_equal(jsol, fsol), \
            f"seed {seed}: eager JAX != Fortran replica (max {np.abs(jsol - fsol).max():.2e})"
    print("PASS eager JAX == Fortran-order numpy replica (0 ULP)")


def test_jit_jax_matches_fortran_replica_to_fma_bound():
    """JIT JAX differs from Fortran replica only by XLA FMA (~1e-13 rel)."""
    solve = jax.jit(penta_lu_solve_jax)
    for seed in range(4):
        lhs, rhs = _make_system(115, seed)
        jsol = np.asarray(solve(jnp.asarray(lhs), jnp.asarray(rhs)))
        fsol = _fortran_penta_replica(lhs, rhs)
        rel = np.abs(jsol - fsol).max() / np.abs(fsol).max()
        assert rel < 1e-12, f"seed {seed}: jit-vs-replica rel {rel:.2e} exceeds FMA bound"
    print("PASS jit JAX == Fortran replica within XLA FMA bound (<1e-12 rel)")


if __name__ == "__main__":
    test_eager_jax_matches_fortran_replica_exactly()
    test_jit_jax_matches_fortran_replica_to_fma_bound()
    print("All penta faithfulness tests passed.")

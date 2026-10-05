"""Unit tests for JAX tridiagonal LU solver."""
from __future__ import annotations
import pytest


import numpy as np


try:
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    HAS_JAX = True
except ImportError:
    HAS_JAX = False

from clubb_jax.src.CLUBB_core.tridiag_lu_solver import tridiag_lu_solve


def call_tridiag_lu_solve(lhs, rhs):
    return tridiag_lu_solve(rhs.shape[-1], lhs, rhs)


def _make_lhs_from_bands(sup_np, mid_np, sub_np):
    """Pack numpy 1-D band arrays into (3, 1, ndim) JAX array."""
    ndim = sup_np.shape[0]
    sup = jnp.array(sup_np[None, :], dtype=jnp.float64)   # (1, ndim)
    mid = jnp.array(mid_np[None, :], dtype=jnp.float64)
    sub = jnp.array(sub_np[None, :], dtype=jnp.float64)
    return jnp.stack([sup, mid, sub], axis=0)              # (3, 1, ndim)


def _dense_tridiag(sup, mid, sub):
    """Build a dense ndim×ndim matrix from band arrays (numpy)."""
    ndim = len(mid)
    A = np.diag(mid)
    A += np.diag(sup[:-1], k=1)   # superdiagonal (above main diag)
    A += np.diag(sub[1:],  k=-1)  # subdiagonal   (below main diag)
    return A


# ──────────────────────────────────────────────────────────────────────────────


def test_solver_vs_numpy_single_col():
    """JAX solver matches numpy.linalg.solve for a random tridiagonal system."""
    if not HAS_JAX:
        pytest.skip("  SKIP")
    rng = np.random.default_rng(7)
    ndim = 10
    sup_np = rng.uniform(-2.0, -0.1, ndim)
    sub_np = rng.uniform(-2.0, -0.1, ndim)
    # Make diagonally dominant so the system is well-conditioned
    mid_np = np.abs(sup_np) + np.abs(sub_np) + 1.0
    sup_np[-1] = 0.0   # top boundary
    sub_np[0] = 0.0    # bottom boundary
    rhs_np = rng.uniform(-10.0, 10.0, ndim)

    A = _dense_tridiag(sup_np, mid_np, sub_np)
    expected = np.linalg.solve(A, rhs_np)

    lhs = _make_lhs_from_bands(sup_np, mid_np, sub_np)
    rhs = jnp.array(rhs_np[None, :], dtype=jnp.float64)
    result = np.asarray(call_tridiag_lu_solve(lhs, rhs))
    assert result.shape == rhs.shape
    soln = result[0]

    err = np.max(np.abs(soln - expected))
    print(f"  single-col vs numpy max_err = {err:.3e}",
          "  PASS" if err < 1e-12 else "  FAIL")
    assert err < 1e-12, f"single-col mismatch: {err}"


def test_solver_vs_numpy_multi_col():
    """JAX solver matches numpy for multiple independent columns."""
    if not HAS_JAX:
        pytest.skip("  SKIP")
    rng = np.random.default_rng(99)
    ndim, ngrdcol = 8, 4
    sup_np = rng.uniform(-3.0, -0.1, (ngrdcol, ndim))
    sub_np = rng.uniform(-3.0, -0.1, (ngrdcol, ndim))
    mid_np = np.abs(sup_np) + np.abs(sub_np) + 2.0
    sup_np[:, -1] = 0.0
    sub_np[:, 0] = 0.0
    rhs_np = rng.uniform(-5.0, 5.0, (ngrdcol, ndim))

    # numpy reference: solve each column independently
    expected = np.zeros_like(rhs_np)
    for i in range(ngrdcol):
        A = _dense_tridiag(sup_np[i], mid_np[i], sub_np[i])
        expected[i] = np.linalg.solve(A, rhs_np[i])

    lhs = jnp.stack([
        jnp.array(sup_np, dtype=jnp.float64),
        jnp.array(mid_np, dtype=jnp.float64),
        jnp.array(sub_np, dtype=jnp.float64),
    ], axis=0)   # (3, ngrdcol, ndim)
    rhs = jnp.array(rhs_np, dtype=jnp.float64)
    soln = np.asarray(call_tridiag_lu_solve(lhs, rhs))
    assert soln.shape == rhs.shape

    err = np.max(np.abs(soln - expected))
    print(f"  multi-col vs numpy max_err = {err:.3e}",
          "  PASS" if err < 1e-12 else "  FAIL")
    assert err < 1e-12, f"multi-col mismatch: {err}"


def test_solver_residual():
    """lhs @ soln ≈ rhs: check residual norm."""
    if not HAS_JAX:
        pytest.skip("  SKIP")
    rng = np.random.default_rng(13)
    ndim, ngrdcol = 12, 2
    sup_np = rng.uniform(-1.0, -0.01, (ngrdcol, ndim))
    sub_np = rng.uniform(-1.0, -0.01, (ngrdcol, ndim))
    mid_np = np.abs(sup_np) + np.abs(sub_np) + 1.5
    sup_np[:, -1] = 0.0
    sub_np[:, 0] = 0.0
    rhs_np = rng.uniform(-3.0, 3.0, (ngrdcol, ndim))

    lhs = jnp.stack([
        jnp.array(sup_np, dtype=jnp.float64),
        jnp.array(mid_np, dtype=jnp.float64),
        jnp.array(sub_np, dtype=jnp.float64),
    ], axis=0)
    rhs = jnp.array(rhs_np, dtype=jnp.float64)
    soln_np = np.asarray(call_tridiag_lu_solve(lhs, rhs))

    # Compute residual: A*soln - rhs for each column
    max_res = 0.0
    for i in range(ngrdcol):
        A = _dense_tridiag(sup_np[i], mid_np[i], sub_np[i])
        res = A @ soln_np[i] - rhs_np[i]
        max_res = max(max_res, np.max(np.abs(res)))

    print(f"  residual max = {max_res:.3e}",
          "  PASS" if max_res < 1e-12 else "  FAIL")
    assert max_res < 1e-12, f"residual too large: {max_res}"


# ──────────────────────────────────────────────────────────────────────────────

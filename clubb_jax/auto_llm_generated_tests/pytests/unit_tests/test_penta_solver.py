"""Compare pentadiagonal solves with NumPy for full and diagonal systems."""
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


from clubb_jax.src.CLUBB_core.penta_lu_solver import penta_lu_solve


def penta_lu_solve_jax(lhs, rhs):
    return penta_lu_solve(rhs.shape[-1], lhs.shape[1], lhs, rhs)


def _dense_penta(s2, s1, d, sb1, sb2):
    """Build dense ndim×ndim matrix from five band arrays (numpy)."""
    ndim = len(d)
    A = np.diag(d)
    A += np.diag(s1[:-1],  k=1)    # 1st superdiagonal
    A += np.diag(s2[:-2],  k=2)    # 2nd superdiagonal
    A += np.diag(sb1[1:],  k=-1)   # 1st subdiagonal
    A += np.diag(sb2[2:],  k=-2)   # 2nd subdiagonal
    return A


def _make_lhs(s2, s1, d, sb1, sb2):
    """Pack 1-D band arrays into (5, 1, ndim) JAX array."""
    return jnp.stack([
        jnp.array(s2[None, :], dtype=jnp.float64),
        jnp.array(s1[None, :], dtype=jnp.float64),
        jnp.array(d[None,  :], dtype=jnp.float64),
        jnp.array(sb1[None, :], dtype=jnp.float64),
        jnp.array(sb2[None, :], dtype=jnp.float64),
    ], axis=0)   # (5, 1, ndim)


# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize('diagonal_only', [False, True])
def test_full_penta_vs_numpy_single_col(diagonal_only):
    """Random full penta-diagonal system vs numpy.linalg.solve."""
    if not HAS_JAX:
        pytest.skip("  SKIP")
    rng = np.random.default_rng(7)
    ndim = 12
    s2  = rng.uniform(-1.5, -0.1, ndim)
    s1  = rng.uniform(-1.5, -0.1, ndim)
    sb1 = rng.uniform(-1.5, -0.1, ndim)
    sb2 = rng.uniform(-1.5, -0.1, ndim)
    # Make diagonally dominant
    d   = np.abs(s2) + np.abs(s1) + np.abs(sb1) + np.abs(sb2) + 1.0
    # Zero out bands at boundaries
    s2[-1] = s2[-2] = 0.0
    s1[-1] = 0.0
    sb1[0] = 0.0
    sb2[0] = sb2[1] = 0.0
    if diagonal_only:
        s2 = s1 = sb1 = sb2 = np.zeros(ndim)
    b = rng.uniform(-10., 10., ndim)

    A = _dense_penta(s2, s1, d, sb1, sb2)
    expected = np.linalg.solve(A, b)

    lhs = _make_lhs(s2, s1, d, sb1, sb2)
    rhs = jnp.array(b[None, :], dtype=jnp.float64)
    result = np.asarray(penta_lu_solve_jax(lhs, rhs))
    assert result.shape == rhs.shape
    soln = result[0]
    err = np.max(np.abs(soln - expected))
    print(f"  single-col vs numpy max_err = {err:.3e}  {'PASS' if err < 1e-12 else 'FAIL'}")
    assert err < 1e-12, f"single-col mismatch: {err}"


def test_multi_col_vs_numpy():
    """Multiple independent columns vs numpy."""
    if not HAS_JAX:
        pytest.skip("  SKIP")
    rng = np.random.default_rng(99)
    ndim, ngrdcol = 10, 4
    s2  = rng.uniform(-1., -0.05, (ngrdcol, ndim))
    s1  = rng.uniform(-1., -0.05, (ngrdcol, ndim))
    sb1 = rng.uniform(-1., -0.05, (ngrdcol, ndim))
    sb2 = rng.uniform(-1., -0.05, (ngrdcol, ndim))
    d   = (np.abs(s2) + np.abs(s1) + np.abs(sb1) + np.abs(sb2)) + 2.0
    s2[:, -1] = s2[:, -2] = 0.0
    s1[:, -1] = 0.0
    sb1[:, 0] = 0.0
    sb2[:, 0] = sb2[:, 1] = 0.0
    b   = rng.uniform(-5., 5., (ngrdcol, ndim))

    expected = np.zeros_like(b)
    for i in range(ngrdcol):
        A = _dense_penta(s2[i], s1[i], d[i], sb1[i], sb2[i])
        expected[i] = np.linalg.solve(A, b[i])

    lhs = jnp.stack([
        jnp.array(s2,  dtype=jnp.float64),
        jnp.array(s1,  dtype=jnp.float64),
        jnp.array(d,   dtype=jnp.float64),
        jnp.array(sb1, dtype=jnp.float64),
        jnp.array(sb2, dtype=jnp.float64),
    ], axis=0)   # (5, ngrdcol, ndim)
    rhs = jnp.array(b, dtype=jnp.float64)
    soln = np.asarray(penta_lu_solve_jax(lhs, rhs))
    assert soln.shape == rhs.shape
    err = np.max(np.abs(soln - expected))
    print(f"  multi-col vs numpy max_err = {err:.3e}  {'PASS' if err < 1e-12 else 'FAIL'}")
    assert err < 1e-12, f"multi-col mismatch: {err}"


# ──────────────────────────────────────────────────────────────────────────────

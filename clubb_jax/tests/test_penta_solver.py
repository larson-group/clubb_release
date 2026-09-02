"""Unit tests for JAX penta-diagonal LU solver.

Tests penta_lu_solve_jax against:
  1. Shape check
  2. Known 5x5 hand-constructed case
  3. numpy.linalg.solve reference on random penta-diagonal systems
  4. Multiple columns (each solved independently)
  5. Fortran penta_lu_solve reference via clubb_python_api

Run:
    python clubb_jax/tests/test_penta_solver.py
"""
from __future__ import annotations

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'clubb_python_api'))

try:
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    HAS_JAX = True
except ImportError:
    HAS_JAX = False

try:
    from clubb_python.clubb_api import penta_lu_solve as fort_penta_lu_solve
    HAS_FORTRAN = True
except Exception:
    HAS_FORTRAN = False

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

def test_shape():
    """Output shape is (ngrdcol, ndim)."""
    if not HAS_JAX:
        print("  SKIP"); return
    ndim, ngrdcol = 8, 3
    lhs = jnp.zeros((5, ngrdcol, ndim)).at[2].set(jnp.ones((ngrdcol, ndim)))
    rhs = jnp.ones((ngrdcol, ndim))
    soln = penta_lu_solve_jax(lhs, rhs)
    assert soln.shape == (ngrdcol, ndim), f"Wrong shape: {soln.shape}"
    print(f"  shape = {soln.shape}  PASS")


def test_diagonal():
    """Pure-diagonal LHS: soln[i] = rhs[i] / diag[i]."""
    if not HAS_JAX:
        print("  SKIP"); return
    ndim = 8
    rng = np.random.default_rng(42)
    d = rng.uniform(1.0, 3.0, ndim)
    b = rng.uniform(-5.0, 5.0, ndim)
    lhs = _make_lhs(np.zeros(ndim), np.zeros(ndim), d, np.zeros(ndim), np.zeros(ndim))
    rhs = jnp.array(b[None, :], dtype=jnp.float64)
    soln = np.asarray(penta_lu_solve_jax(lhs, rhs))[0]
    expected = b / d
    err = np.max(np.abs(soln - expected))
    print(f"  diagonal max_err = {err:.3e}  {'PASS' if err < 1e-14 else 'FAIL'}")
    assert err < 1e-14, f"diagonal mismatch: {err}"


def test_known_5x5():
    """5x5 tridiagonal-like system (zero 2nd off-diagonals), verified against numpy."""
    if not HAS_JAX:
        print("  SKIP"); return
    ndim = 5
    # Simple diagonally dominant tridiagonal embedded in penta storage
    d   = np.array([4., 4., 4., 4., 4.])
    s1  = np.array([-1., -1., -1., -1., 0.])
    s2  = np.zeros(ndim)
    sb1 = np.array([0., -1., -1., -1., -1.])
    sb2 = np.zeros(ndim)
    b   = np.array([1., 0., 0., 0., 1.])

    A = _dense_penta(s2, s1, d, sb1, sb2)
    expected = np.linalg.solve(A, b)

    lhs = _make_lhs(s2, s1, d, sb1, sb2)
    rhs = jnp.array(b[None, :], dtype=jnp.float64)
    soln = np.asarray(penta_lu_solve_jax(lhs, rhs))[0]
    err = np.max(np.abs(soln - expected))
    print(f"  5x5 max_err = {err:.3e}  {'PASS' if err < 1e-13 else 'FAIL'}")
    assert err < 1e-13, f"5x5 mismatch: {err}"


def test_full_penta_vs_numpy_single_col():
    """Random full penta-diagonal system vs numpy.linalg.solve."""
    if not HAS_JAX:
        print("  SKIP"); return
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
    b = rng.uniform(-10., 10., ndim)

    A = _dense_penta(s2, s1, d, sb1, sb2)
    expected = np.linalg.solve(A, b)

    lhs = _make_lhs(s2, s1, d, sb1, sb2)
    rhs = jnp.array(b[None, :], dtype=jnp.float64)
    soln = np.asarray(penta_lu_solve_jax(lhs, rhs))[0]
    err = np.max(np.abs(soln - expected))
    print(f"  single-col vs numpy max_err = {err:.3e}  {'PASS' if err < 1e-12 else 'FAIL'}")
    assert err < 1e-12, f"single-col mismatch: {err}"


def test_multi_col_vs_numpy():
    """Multiple independent columns vs numpy."""
    if not HAS_JAX:
        print("  SKIP"); return
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
    err = np.max(np.abs(soln - expected))
    print(f"  multi-col vs numpy max_err = {err:.3e}  {'PASS' if err < 1e-12 else 'FAIL'}")
    assert err < 1e-12, f"multi-col mismatch: {err}"


def test_vs_fortran_single_rhs():
    """JAX penta solver matches Fortran penta_lu_solve bit-for-bit."""
    if not (HAS_JAX and HAS_FORTRAN):
        print("  SKIP (no Fortran API)"); return

    rng = np.random.default_rng(1234)
    ndim = 20
    ngrdcol = 2
    s2  = rng.uniform(-1., -0.05, (ngrdcol, ndim)).astype(np.float64)
    s1  = rng.uniform(-1., -0.05, (ngrdcol, ndim)).astype(np.float64)
    sb1 = rng.uniform(-1., -0.05, (ngrdcol, ndim)).astype(np.float64)
    sb2 = rng.uniform(-1., -0.05, (ngrdcol, ndim)).astype(np.float64)
    d   = (np.abs(s2) + np.abs(s1) + np.abs(sb1) + np.abs(sb2) + 2.0).astype(np.float64)
    s2[:, -1] = s2[:, -2] = 0.0
    s1[:, -1] = 0.0
    sb1[:, 0] = 0.0
    sb2[:, 0] = sb2[:, 1] = 0.0
    b = rng.uniform(-5., 5., (ngrdcol, ndim)).astype(np.float64)

    # Fortran call: penta_lu_solve(ndim, ngrdcol, lhs, rhs)
    # lhs must be (-2:2, ngrdcol, ndim) → Python (5, ngrdcol, ndim)
    lhs_fort = np.stack([s2, s1, d, sb1, sb2], axis=0).copy()  # (5, ngrdcol, ndim)
    soln_fort = np.asarray(fort_penta_lu_solve(ndim, ngrdcol, lhs_fort, b.copy()))

    # JAX call
    lhs_jax = jnp.array(np.stack([s2, s1, d, sb1, sb2], axis=0), dtype=jnp.float64)
    rhs_jax = jnp.array(b, dtype=jnp.float64)
    soln_jax = np.asarray(penta_lu_solve_jax(lhs_jax, rhs_jax))

    err = np.max(np.abs(soln_jax - soln_fort))
    print(f"  vs Fortran max_err = {err:.3e}  {'PASS' if err < 1e-14 else 'FAIL'}")
    assert err < 1e-14, f"Fortran mismatch: {err}"


# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 60)
    print("CLUBB JAX penta-diagonal solver unit tests")
    print("=" * 60)

    if not HAS_JAX:
        print("ERROR: JAX not available.")
        sys.exit(1)

    tests = [
        test_shape,
        test_diagonal,
        test_known_5x5,
        test_full_penta_vs_numpy_single_col,
        test_multi_col_vs_numpy,
        test_vs_fortran_single_rhs,
    ]

    passed = failed = 0
    for t in tests:
        print(f"\n{t.__name__}:")
        try:
            t()
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {e}")
            failed += 1
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  ERROR: {type(e).__name__}: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed}/{passed+failed} passed, {failed} failed")
    print("=" * 60)
    sys.exit(0 if failed == 0 else 1)

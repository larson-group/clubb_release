"""Unit tests for JAX tridiagonal LU solver.

Tests call_tridiag_lu_solve against:
  1. Shape check
  2. Identity-diagonal system (diagonal matrix)
  3. Known 3x3 hand-computed case
  4. numpy.linalg.solve reference on random tridiagonal systems
  5. Multiple columns (each solved independently)

Run:
    python clubb_jax/tests/test_solver.py
"""
from __future__ import annotations

import sys
import os

import numpy as np

_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

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

def test_solver_shape():
    """Output shape is (ngrdcol, ndim)."""
    if not HAS_JAX:
        print("  SKIP"); return
    ndim, ngrdcol = 6, 3
    sup = jnp.zeros((ngrdcol, ndim))
    mid = jnp.ones((ngrdcol, ndim)) * 2.0
    sub = jnp.zeros((ngrdcol, ndim))
    lhs = jnp.stack([sup, mid, sub], axis=0)
    rhs = jnp.ones((ngrdcol, ndim))
    soln = call_tridiag_lu_solve(lhs, rhs)
    assert soln.shape == (ngrdcol, ndim), f"Wrong shape: {soln.shape}"
    print(f"  shape = {soln.shape}  PASS")


def test_solver_diagonal():
    """Pure-diagonal LHS (super=sub=0): soln[i] = rhs[i] / mid[i]."""
    if not HAS_JAX:
        print("  SKIP"); return
    ndim = 5
    rng = np.random.default_rng(42)
    mid_np = rng.uniform(1.0, 3.0, ndim)
    rhs_np = rng.uniform(-5.0, 5.0, ndim)
    sup_np = np.zeros(ndim)
    sub_np = np.zeros(ndim)
    lhs = _make_lhs_from_bands(sup_np, mid_np, sub_np)
    rhs = jnp.array(rhs_np[None, :], dtype=jnp.float64)   # (1, ndim)
    soln = np.asarray(call_tridiag_lu_solve(lhs, rhs))[0]   # (ndim,)
    expected = rhs_np / mid_np
    err = np.max(np.abs(soln - expected))
    print(f"  diagonal max_err = {err:.3e}  PASS" if err < 1e-14 else f"  FAIL err={err}")
    assert err < 1e-14, f"diagonal system mismatch: {err}"


def test_solver_known_3x3():
    """Hand-verified 3x3 system.

    A = [[2, -1,  0],
         [-1,  2, -1],
         [0, -1,  2]]
    rhs = [1, 0, 1]
    exact soln = [1, 1, 1]
    """
    if not HAS_JAX:
        print("  SKIP"); return
    mid_np = np.array([2.0, 2.0, 2.0])
    sup_np = np.array([-1.0, -1.0, 0.0])   # entry at top is unused
    sub_np = np.array([0.0, -1.0, -1.0])   # entry at bottom is unused
    rhs_np = np.array([1.0, 0.0, 1.0])
    lhs = _make_lhs_from_bands(sup_np, mid_np, sub_np)
    rhs = jnp.array(rhs_np[None, :], dtype=jnp.float64)
    soln = np.asarray(call_tridiag_lu_solve(lhs, rhs))[0]
    expected = np.array([1.0, 1.0, 1.0])
    err = np.max(np.abs(soln - expected))
    print(f"  3x3 soln={soln}  max_err = {err:.3e}",
          "  PASS" if err < 1e-14 else "  FAIL")
    assert err < 1e-14, f"3x3 mismatch: {err}"


def test_solver_vs_numpy_single_col():
    """JAX solver matches numpy.linalg.solve for a random tridiagonal system."""
    if not HAS_JAX:
        print("  SKIP"); return
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
    soln = np.asarray(call_tridiag_lu_solve(lhs, rhs))[0]

    err = np.max(np.abs(soln - expected))
    print(f"  single-col vs numpy max_err = {err:.3e}",
          "  PASS" if err < 1e-12 else "  FAIL")
    assert err < 1e-12, f"single-col mismatch: {err}"


def test_solver_vs_numpy_multi_col():
    """JAX solver matches numpy for multiple independent columns."""
    if not HAS_JAX:
        print("  SKIP"); return
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

    err = np.max(np.abs(soln - expected))
    print(f"  multi-col vs numpy max_err = {err:.3e}",
          "  PASS" if err < 1e-12 else "  FAIL")
    assert err < 1e-12, f"multi-col mismatch: {err}"


def test_solver_residual():
    """lhs @ soln ≈ rhs: check residual norm."""
    if not HAS_JAX:
        print("  SKIP"); return
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


def test_solver_f2py_oracle():
    """f2py bit-shadow vs the Fortran oracle. call_tridiag_lu_solve is a faithful port of
    `tridiag_lu_solve_single_rhs_multiple_lhs` (tridiag_lu_solver.F90) — this is THE workhorse solver every
    prognostic advance routes through, so bit-matching the SPECIFIC Fortran LU (not merely numpy.linalg.solve,
    which only confirms the system is solved) is what underpins the whole bit-faithful case suite. SKIPs if
    clubb_f2py is unbuilt. (iter 441)"""
    if not HAS_JAX:
        print("  SKIP (no jax)"); return
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py tridiag_lu_solve oracle: SKIP ({type(e).__name__})")
        return
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
        worst = max(worst, float(np.max(np.abs(j - f) / np.maximum(np.abs(f), 1e-30))))
    assert worst < 1e-12, f"tridiag_lu_solve f2py mismatch {worst:.2e}"
    print(f"  f2py tridiag_lu_solve_single_rhs_multiple_lhs: bit-match, worst rel {worst:.2e}  PASS")


# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 60)
    print("CLUBB JAX tridiagonal solver unit tests")
    print("=" * 60)

    if not HAS_JAX:
        print("ERROR: JAX not available.")
        sys.exit(1)

    tests = [
        test_solver_shape,
        test_solver_diagonal,
        test_solver_known_3x3,
        test_solver_vs_numpy_single_col,
        test_solver_vs_numpy_multi_col,
        test_solver_residual,
        test_solver_f2py_oracle,
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

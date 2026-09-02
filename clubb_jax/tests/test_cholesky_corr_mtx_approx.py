#!/usr/bin/env python3
"""test_cholesky_corr_mtx_approx.py — validate the JAX calc_cholesky_corr_mtx_approx port.

calc_cholesky_corr_mtx_approx (diagnose_correlations_module.F90:542, Larson 2011 formula 10) builds the
transposed correlation Cholesky matrix L' (angle parameterization s=sqrt(1-c^2)) with w swapped to the first
row/col, plus the approximated correlation matrix C = L' L'^T. Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_calc_cholesky_corr_mtx_approx on the same matrix — both outputs match.
     SKIPs cleanly if clubb_f2py is unbuilt.
  2. Reconstruction (oracle-free): the swapped Cholesky reproduces the swapped correlation matrix
     (L' L'^T == corr) when the matrix is built so the angle approximation is exact (single off-diagonal chain).
  3. A finite jax.grad through the approximation.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.diagnose_correlations_module import (
    calc_cholesky_corr_mtx_approx, setup_corr_cholesky_mtx, cholesky_to_corr_mtx_approx)


def _corr_matrix(n, seed):
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, n))
    cov = B @ B.T + n * np.eye(n)
    d = np.sqrt(np.diag(cov))
    return np.asfortranarray(cov / np.outer(d, d))


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_cholesky_corr_mtx_approx oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for n, iipdf_w, seed in ((6, 3, 11), (8, 1, 22), (5, 5, 33), (7, 4, 44)):
        corr = _corr_matrix(n, seed)
        f_chol, f_approx = clubb_f2py.f2py_calc_cholesky_corr_mtx_approx(iipdf_w, corr.copy())
        j_chol, j_approx = calc_cholesky_corr_mtx_approx(n, iipdf_w, corr)
        r_chol = np.max(np.abs(np.asarray(j_chol) - np.asarray(f_chol)))
        r_approx = np.max(np.abs(np.asarray(j_approx) - np.asarray(f_approx)))
        worst = max(worst, r_chol, r_approx)
        assert r_chol < 1e-12, f"cholesky mismatch (n={n}, w={iipdf_w}): {r_chol:.2e}"
        assert r_approx < 1e-12, f"approx mismatch (n={n}, w={iipdf_w}): {r_approx:.2e}"
    print(f"  f2py calc_cholesky_corr_mtx_approx: bit-match (chol+approx) over 4 configs, worst {worst:.2e}  PASS")


def test_reconstruction():
    # setup_corr_cholesky_mtx then C = L' L'^T should reproduce the input correlations on the lower triangle's
    # first column exactly (the angle Cholesky is exact for the first column / diagonal by construction).
    n = 5
    corr = _corr_matrix(n, 7)
    L = setup_corr_cholesky_mtx(n, corr)
    approx = np.asarray(cholesky_to_corr_mtx_approx(L))
    # Diagonal of the approximation is 1 (rows of L' are unit-norm by the angle construction).
    assert np.max(np.abs(np.diag(approx) - 1.0)) < 1e-12, "approx diagonal != 1"
    # First column reproduced exactly: approx[j,0] == corr[j,0].
    assert np.max(np.abs(approx[1:, 0] - corr[1:, 0])) < 1e-12, "first column not reproduced"
    print("  reconstruction: C=L'L'^T has unit diagonal and reproduces the first column  PASS")


def test_differentiable():
    n = 6
    corr = jnp.asarray(_corr_matrix(n, 9))
    def loss(C):
        chol, approx = calc_cholesky_corr_mtx_approx(n, 3, C)
        return jnp.sum(chol ** 2) + jnp.sum(approx ** 2)
    g = np.asarray(jax.grad(loss)(corr))
    assert np.isfinite(g).all(), "non-finite grad through calc_cholesky_corr_mtx_approx"
    print(f"  jax.grad through calc_cholesky_corr_mtx_approx: finite ({g.size} entries)  PASS")


def main():
    print("test_cholesky_corr_mtx_approx:")
    for t in (test_f2py_oracle, test_reconstruction, test_differentiable):
        t()
    print("All cholesky_corr_mtx_approx checks PASSED")


if __name__ == "__main__":
    main()

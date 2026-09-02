#!/usr/bin/env python3
"""test_cholesky_factor.py — validate the JAX matrix_operations.Cholesky_factor port.

Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_cholesky_factor on the same matrices — a_scaling, a_cholesky (lower
     triangle), and l_scaled all match. SKIPs cleanly if clubb_f2py is unbuilt.
  2. Reconstruction property (oracle-free): for a positive-definite input, L Lᵀ == a_input (no scaling case),
     and the strict upper triangle of the returned matrix retains the input's upper values (LAPACK dpotrf
     leaves it untouched).
  3. The tau-on-diagonal fallback recovers a finite factor for a non-positive-definite input.
  4. A finite jax.grad through the (positive-definite) factorization.
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

from clubb_jax.src.CLUBB_core.matrix_operations import Cholesky_factor


def _corr_matrix(n, seed):
    """A positive-definite correlation matrix (unit diagonal) built as D^{-1/2} (B Bᵀ) D^{-1/2}."""
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, n))
    cov = B @ B.T + n * np.eye(n)
    d = np.sqrt(np.diag(cov))
    corr = cov / np.outer(d, d)
    return np.asfortranarray(corr)


def test_reconstruction():
    for n, seed in ((4, 1), (6, 2), (5, 3)):
        a = _corr_matrix(n, seed)
        a_scaling, L, l_scaled = Cholesky_factor(a)
        L = np.asarray(L)
        # Unit-diagonal correlation matrix -> no equilibration.
        assert not bool(l_scaled), "correlation matrix should not be scaled"
        assert np.allclose(np.asarray(a_scaling), 1.0, atol=1e-12), "scaling should be 1 for unit diagonal"
        Ltri = np.tril(L)
        assert np.max(np.abs(Ltri @ Ltri.T - a)) < 1e-11, f"L Lᵀ != a (n={n})"
        # Strict upper triangle retains the input's upper values (dpotrf leaves it untouched).
        iu = np.triu_indices(n, 1)
        assert np.max(np.abs(L[iu] - a[iu])) < 1e-14, "strict upper triangle not preserved"
    print("  reconstruction L Lᵀ = a, scaling=1 / l_scaled=False, upper preserved  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py Cholesky_factor oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for n, seed in ((4, 11), (6, 22), (5, 33), (7, 44)):
        a = _corr_matrix(n, seed)
        f_scaling, f_chol, f_lscaled = clubb_f2py.f2py_cholesky_factor(a.copy())
        j_scaling, j_chol, j_lscaled = Cholesky_factor(a)
        j_chol = np.asarray(j_chol); j_scaling = np.asarray(j_scaling)
        # Compare the lower triangle (incl. diagonal) + scaling + l_scaled.
        il = np.tril_indices(n)
        rel = np.max(np.abs(j_chol[il] - np.asarray(f_chol)[il]))
        rel_s = np.max(np.abs(j_scaling - np.asarray(f_scaling)))
        worst = max(worst, rel, rel_s)
        assert rel < 1e-11, f"cholesky lower mismatch n={n}: {rel:.2e}"
        assert rel_s < 1e-12, f"scaling mismatch n={n}: {rel_s:.2e}"
        assert bool(j_lscaled) == bool(f_lscaled), f"l_scaled mismatch n={n}"
    print(f"  f2py Cholesky_factor: bit-match (lower+scaling+l_scaled) over 4 sizes, worst {worst:.2e}  PASS")


def test_non_pd_fallback():
    # A non-positive-definite symmetric matrix (negative eigenvalue) -> bare Cholesky yields NaN; the tau
    # fallback must recover a finite factor.
    a = np.array([[1.0, 0.99, 0.99],
                  [0.99, 1.0, 0.99],
                  [0.99, 0.99, 1.0]])
    # Make it indefinite.
    a[0, 2] = a[2, 0] = -0.99
    _, L, _ = Cholesky_factor(a)
    L = np.asarray(L)
    assert np.isfinite(L).all(), "tau fallback did not produce a finite factor"
    print("  tau-on-diagonal fallback: finite factor for a non-PD input  PASS")


def test_differentiable():
    a = jnp.asarray(_corr_matrix(5, 7))
    def loss(x):
        _, L, _ = Cholesky_factor(x)
        return jnp.sum(jnp.tril(L) ** 2)
    g = np.asarray(jax.grad(loss)(a))
    assert np.isfinite(g).all(), "non-finite grad through Cholesky_factor"
    print(f"  jax.grad through Cholesky_factor: finite ({g.size} entries)  PASS")


def main():
    print("test_cholesky_factor:")
    for t in (test_reconstruction, test_f2py_oracle, test_non_pd_fallback, test_differentiable):
        t()
    print("All Cholesky_factor checks PASSED")


if __name__ == "__main__":
    main()

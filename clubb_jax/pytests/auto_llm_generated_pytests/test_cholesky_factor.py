#!/usr/bin/env python3
"""validate the JAX matrix_operations.Cholesky_factor port."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

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

#!/usr/bin/env python3
"""validate the JAX calc_cholesky_corr_mtx_approx port."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.diagnose_correlations_module import setup_corr_cholesky_mtx, cholesky_to_corr_mtx_approx


def _corr_matrix(n, seed):
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, n))
    cov = B @ B.T + n * np.eye(n)
    d = np.sqrt(np.diag(cov))
    return np.asfortranarray(cov / np.outer(d, d))


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

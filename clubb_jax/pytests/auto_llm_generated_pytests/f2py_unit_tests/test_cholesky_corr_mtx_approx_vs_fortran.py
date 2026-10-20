"""Check the JAX calc_cholesky_corr_mtx_approx port against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.diagnose_correlations_module import calc_cholesky_corr_mtx_approx


def _corr_matrix(n, seed):
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, n))
    cov = B @ B.T + n * np.eye(n)
    d = np.sqrt(np.diag(cov))
    return np.asfortranarray(cov / np.outer(d, d))


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    for n, iipdf_w, seed in ((6, 3, 11), (8, 1, 22), (5, 5, 33), (7, 4, 44)):
        corr = _corr_matrix(n, seed)
        f_chol, f_approx = clubb_f2py.f2py_calc_cholesky_corr_mtx_approx(iipdf_w, corr.copy())
        j_chol, j_approx = calc_cholesky_corr_mtx_approx(n, iipdf_w, corr)
        r_chol = max_abs_diff("cholesky_corr_mtx_approx", np.asarray(j_chol), np.asarray(f_chol))
        r_approx = max_abs_diff("cholesky_corr_mtx_approx", np.asarray(j_approx), np.asarray(f_approx))
        worst = max(worst, r_chol, r_approx)
        assert r_chol < 1e-12, f"cholesky mismatch (n={n}, w={iipdf_w}): {r_chol:.2e}"
        assert r_approx < 1e-12, f"approx mismatch (n={n}, w={iipdf_w}): {r_approx:.2e}"

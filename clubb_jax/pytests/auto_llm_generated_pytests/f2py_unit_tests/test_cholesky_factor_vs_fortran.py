"""Check the JAX matrix_operations.Cholesky_factor port against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.matrix_operations import Cholesky_factor


def _corr_matrix(n, seed):
    """A positive-definite correlation matrix (unit diagonal) built as D^{-1/2} (B Bᵀ) D^{-1/2}."""
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((n, n))
    cov = B @ B.T + n * np.eye(n)
    d = np.sqrt(np.diag(cov))
    corr = cov / np.outer(d, d)
    return np.asfortranarray(corr)


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    for n, seed in ((4, 11), (6, 22), (5, 33), (7, 44)):
        a = _corr_matrix(n, seed)
        f_scaling, f_chol, f_lscaled = clubb_f2py.f2py_cholesky_factor(a.copy())
        j_scaling, j_chol, j_lscaled = Cholesky_factor(a)
        j_chol = np.asarray(j_chol); j_scaling = np.asarray(j_scaling)
        # Compare the lower triangle (incl. diagonal) + scaling + l_scaled.
        il = np.tril_indices(n)
        rel = max_abs_diff("cholesky_factor", j_chol[il], np.asarray(f_chol)[il])
        rel_s = max_abs_diff("cholesky_factor", j_scaling, np.asarray(f_scaling))
        worst = max(worst, rel, rel_s)
        assert rel < 1e-11, f"cholesky lower mismatch n={n}: {rel:.2e}"
        assert rel_s < 1e-12, f"scaling mismatch n={n}: {rel_s:.2e}"
        assert bool(j_lscaled) == bool(f_lscaled), f"l_scaled mismatch n={n}"

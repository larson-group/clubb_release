"""Check the JAX matrix_operations.Cholesky_factor port against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
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



def test_compiled_scaling_retry_and_lower_triangle_contracts():
    """Compare the actual Fortran/JAX routines on paths normal PD cases miss."""
    import jax
    import jax.numpy as jnp
    from clubb_python.CLUBB_core.matrix_operations import cholesky_factor

    inputs = [
        np.array([[1e-4, np.nan], [0.002, 4.0]]),
        np.diag([1e-300, 2e-300]),
        np.diag([1e300, 2e300]),
        np.array([[1.0, 17.0], [1.05, 1.0]]),
        np.array([[1.0, 17.0], [4.0, 1.0]]),
        np.array([[-1.0, 17.0], [0.0, 1.0]]),
        np.empty((0, 0)),
    ]
    compiled = jax.jit(Cholesky_factor)
    for a in inputs:
        f_scaling, f_chol, f_lscaled = cholesky_factor(a.shape[0], a.copy())
        j_scaling, j_chol, j_lscaled = compiled(jnp.asarray(a))
        assert j_chol.shape == f_chol.shape == a.shape
        assert j_scaling.shape == f_scaling.shape == (a.shape[0],)
        np.testing.assert_allclose(j_scaling, f_scaling, atol=0.0, rtol=1e-13)
        np.testing.assert_allclose(j_chol, f_chol, atol=1e-14, rtol=1e-13, equal_nan=True)
        assert bool(j_lscaled) == bool(f_lscaled)

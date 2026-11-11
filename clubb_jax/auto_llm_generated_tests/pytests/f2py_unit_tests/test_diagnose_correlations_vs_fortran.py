"""Check the JAX diagnose_correlations_module port against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_rel_diff
from clubb_jax.src.CLUBB_core.diagnose_correlations_module import diagnose_correlations


def _corr_matrix(n, seed):
    """A symmetric matrix with unit diagonal and off-diagonals in (-0.9, 0.9) — a valid 'prescribed' input."""
    rng = np.random.default_rng(seed)
    a = rng.uniform(-0.9, 0.9, size=(n, n))
    a = 0.5 * (a + a.T)
    np.fill_diagonal(a, 1.0)
    return np.asfortranarray(a)


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    for n, iipdf_w, seed in ((6, 3, 11), (8, 1, 22), (5, 5, 33), (7, 4, 44)):
        corr_pre = _corr_matrix(n, seed)
        ref = np.asarray(clubb_f2py.f2py_diagnose_correlations(iipdf_w, corr_pre.copy(), 0))
        got = np.asarray(diagnose_correlations(n, iipdf_w, corr_pre, l_calc_w_corr=False))
        rel = max_rel_diff("diagnose_correlations", got, ref, 1e-30)
        worst = max(worst, rel)
        assert rel < 1e-12, f"f2py mismatch (n={n}, iipdf_w={iipdf_w}): rel {rel:.2e}"

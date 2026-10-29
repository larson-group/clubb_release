"""Check the JAX compute_mean_binormal / compute_variance_binormal ports
against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.pdf_utilities import compute_mean_binormal, compute_variance_binormal


def test_f2py_oracle():
    import clubb_f2py
    rng = np.random.default_rng(3)
    worst_m = worst_v = 0.0
    for _ in range(50):
        mu1, mu2 = rng.uniform(-3, 3, 2)
        s1, s2 = rng.uniform(0.05, 2.0, 2)
        a = rng.uniform(0.0, 1.0)
        xm_ref = float(clubb_f2py.f2py_compute_mean_binormal(mu1, mu2, a))
        xm_j = float(compute_mean_binormal(mu1, mu2, a))
        v_ref = float(clubb_f2py.f2py_compute_variance_binormal(xm_ref, mu1, mu2, s1, s2, a))
        v_j = float(compute_variance_binormal(xm_j, mu1, mu2, s1, s2, a))
        worst_m = max(worst_m, max_abs_diff("compute_mean_binormal", xm_j, xm_ref))
        worst_v = max(worst_v, max_abs_diff("compute_variance_binormal", v_j, v_ref))
    assert worst_m < 1e-13 and worst_v < 1e-12, f"f2py mismatch mean {worst_m:.2e} var {worst_v:.2e}"

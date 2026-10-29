"""Check the JAX smooth_heaviside_peskin port (advance_helper_module.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.advance_helper as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.advance_helper_module import smooth_heaviside_peskin


def test_f2py_oracle():
    """smooth_heaviside_peskin matches Fortran, including CLUBB's truncated pi and invrs_pi constants."""
    rng = np.random.default_rng(4)
    worst = 0.0
    for smth_range in (0.6, 1.0, 0.25):
        x = rng.uniform(-3 * smth_range, 3 * smth_range, (2, 20))
        args = (x.shape[1], x.shape[0], x, smth_range)
        worst = max(worst, max_abs_diff("smooth_heaviside_peskin", smooth_heaviside_peskin(*args),
                                        fortran.smooth_heaviside_peskin(*args)))
    assert worst < 1e-15, f"smooth_heaviside_peskin Fortran mismatch {worst:.2e}"

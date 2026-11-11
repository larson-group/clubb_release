"""Check the JAX calculate_spurious_source port (advance_clubb_core_module)
against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.numerical_check import calculate_spurious_source


def test_f2py_oracle():
    import clubb_f2py
    rng = np.random.default_rng(4)
    worst = 0.0
    dt = 60.0
    for _ in range(200):
        ia, ib, ft, fs, ifc = rng.uniform(-10, 10, 5)
        ref = float(clubb_f2py.f2py_calculate_spurious_source(ia, ib, ft, fs, ifc, dt))
        got = float(calculate_spurious_source(ia, ib, ft, fs, ifc, dt))
        worst = max(worst, max_abs_diff("spurious_source", got, ref))
    assert worst < 1e-12, f"spurious_source f2py mismatch {worst:.2e}"

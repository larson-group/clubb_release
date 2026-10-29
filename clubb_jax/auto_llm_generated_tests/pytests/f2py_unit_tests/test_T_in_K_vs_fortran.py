"""Check the T_in_K_module port (thlm <-> absolute-T conversions)
against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K


def test_f2py_oracle():
    """Bit-shadow thlm2T_in_K vs f2py_thlm2t_in_k_1d (the Fortran absolute-temperature diagnosis used by the
    radiation/microphysics handoff)."""
    import clubb_f2py
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(15):
        nz = 20
        thlm = rng.uniform(280.0, 320.0, nz); exner = rng.uniform(0.6, 1.0, nz); rcm = rng.uniform(0.0, 1e-3, nz)
        j = np.asarray(thlm2T_in_K(thlm, exner, rcm))
        f = np.asarray(clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm))
        worst = max(worst, float(max_abs_diff("T_in_K", j, f)))
    assert worst < 1e-11, f"thlm2T_in_K f2py mismatch {worst:.2e}"

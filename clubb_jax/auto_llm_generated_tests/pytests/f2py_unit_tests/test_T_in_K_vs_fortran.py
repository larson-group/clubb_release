"""Check the T_in_K_module port (thlm <-> absolute-T conversions)
against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K


def test_f2py_oracle():
    """Bit-shadow thlm2T_in_K vs thlm2t_in_k (the Fortran absolute-temperature diagnosis used by the
    radiation/microphysics handoff)."""
    from clubb_python import clubb_api
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(15):
        nz = 20
        thlm = rng.uniform(280.0, 320.0, nz); exner = rng.uniform(0.6, 1.0, nz); rcm = rng.uniform(0.0, 1e-3, nz)
        j = np.asarray(thlm2T_in_K(thlm, exner, rcm))
        f = np.asarray(clubb_api.thlm2t_in_k(nz, 1, thlm[None, :], exner[None, :], rcm[None, :]))
        assert f.shape == (1, nz)
        f = f[0]
        worst = max(worst, float(max_abs_diff("T_in_K", j, f)))
    assert worst < 1e-11, f"thlm2T_in_K f2py mismatch {worst:.2e}"

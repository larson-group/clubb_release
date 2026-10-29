"""Check the JAX calculate_thvm port (calc_pressure.F90:calculate_thvm)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.calc_pressure as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.calc_pressure import calculate_thvm

NG, NZ = 2, 8


def _inputs(seed):
    rng = np.random.default_rng(seed)
    thlm = rng.uniform(280, 310, (NG, NZ))
    rtm = rng.uniform(0.0, 0.02, (NG, NZ))
    rcm = rng.uniform(0.0, 1e-3, (NG, NZ))
    exner = rng.uniform(0.8, 1.0, (NG, NZ))
    thv_ds_zt = rng.uniform(290, 305, (NG, NZ))
    return thlm, rtm, rcm, exner, thv_ds_zt


def test_f2py_oracle():
    """calculate_thvm matches Fortran."""
    worst = 0.0
    for seed in (11, 22, 33):
        thlm, rtm, rcm, exner, thv_ds_zt = _inputs(seed)
        ng, nzt = thlm.shape
        args = (nzt, ng, thlm, rtm, rcm, exner, thv_ds_zt)
        worst = max(worst, max_abs_diff("calculate_thvm", calculate_thvm(*args), fortran.calculate_thvm(*args)))
    assert worst < 1e-11, f"calculate_thvm Fortran mismatch {worst:.2e}"

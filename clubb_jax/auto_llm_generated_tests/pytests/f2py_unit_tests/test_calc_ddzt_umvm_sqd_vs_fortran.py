"""Check advance_helper_module.py:calc_ddzt_umvm_sqd against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.advance_helper as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.advance_helper_module import calc_ddzt_umvm_sqd

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_f2py_oracle():
    """Squared vertical shear of the mean wind matches Fortran."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(12):
        um = rng.uniform(-15.0, 15.0, (ng, nzt))
        vm = rng.uniform(-15.0, 15.0, (ng, nzt))
        worst = max(worst, max_abs_diff("calc_ddzt_umvm_sqd", calc_ddzt_umvm_sqd(nzm, nzt, ng, jgr, um, vm),
                                        fortran.calc_ddzt_umvm_sqd(nzm, nzt, ng, fgr, um, vm)))
    assert worst < 1e-12, f"calc_ddzt_umvm_sqd Fortran mismatch {worst:.2e}"

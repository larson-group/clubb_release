"""Check advance_helper_module.py:calc_wp3_on_wp2 against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.advance_helper as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.advance_helper_module import calc_wp3_on_wp2

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_f2py_oracle():
    """calc_wp3_on_wp2 (wp3_on_wp2 on zm and zt) matches Fortran, with the +/-1000 clip exercised."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(6)
    worst = 0.0
    for it in range(15):
        wp2 = rng.uniform(1e-4, 3.0, (ng, nzm))
        wp3 = rng.uniform(-5.0, 5.0, (ng, nzt))
        if it % 3 == 0:
            wp2[:, ::4] = 1e-7   # tiny wp2 makes |wp3/wp2| exceed the clip
        worst = max(worst, max_abs_diff("calc_wp3_on_wp2",
                                        calc_wp3_on_wp2(nzm, nzt, ng, jgr, jnp.asarray(wp2), jnp.asarray(wp3)),
                                        fortran.calc_wp3_on_wp2(nzm, nzt, ng, fgr, wp2, wp3)))
    assert worst < 1e-12, f"calc_wp3_on_wp2 Fortran mismatch {worst:.2e}"

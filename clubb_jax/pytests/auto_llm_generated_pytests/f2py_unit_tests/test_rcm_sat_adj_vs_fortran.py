"""Check the JAX rcm_sat_adj port (saturation.F90:rcm_sat_adj) against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.saturation import rcm_sat_adj, SATURATION_FLATAU, SATURATION_BOLTON


def _cases():
    rng = np.random.default_rng(3)
    cases = []
    for _ in range(40):
        thlm = rng.uniform(270.0, 305.0)
        exner = rng.uniform(0.85, 1.0)
        p = 1.0e5 * exner ** (1.0035 / 0.286)        # roughly consistent pressure
        rtm = rng.uniform(0.0, 0.03)                  # spans unsaturated and supersaturated
        cases.append((thlm, rtm, p, exner))
    return cases


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    n_sat = 0
    for formula in (SATURATION_FLATAU, SATURATION_BOLTON):
        for thlm, rtm, p, exner in _cases():
            ref = float(clubb_f2py.f2py_rcm_sat_adj(thlm, rtm, p, exner, formula))
            got = float(rcm_sat_adj(thlm, rtm, p, exner, formula))
            worst = max(worst, max_abs_diff("rcm_sat_adj", got, ref))
            if ref > 1e-8:
                n_sat += 1
    assert worst < 1e-12, f"rcm_sat_adj f2py mismatch {worst:.2e}"
    assert n_sat > 0, "no saturated cases exercised"

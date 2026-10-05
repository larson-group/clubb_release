#!/usr/bin/env python3
"""validate the JAX rcm_sat_adj port (saturation.F90:rcm_sat_adj)."""


import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.saturation import rcm_sat_adj, sat_mixrat_liq, SATURATION_FLATAU

_Cp, _Lv = 1004.67, 2.5e6


def test_self_consistency():
    # For a saturated point, rcm and the implied theta should satisfy the fixed point within ~tolerance.
    thlm, rtm, p, exner, formula = 285.0, 0.02, 95000.0, 0.95, SATURATION_FLATAU
    rcm = float(rcm_sat_adj(thlm, rtm, p, exner, formula))
    assert rcm > 0.0, "expected a saturated point"
    theta = thlm + (_Lv / (_Cp * exner)) * rcm
    rsat = float(sat_mixrat_liq(p, theta * exner, formula))
    rcm_check = max(rtm - rsat, 0.0)
    assert abs(rcm - rcm_check) < 1e-5, f"fixed point not satisfied: {abs(rcm-rcm_check):.2e}"
    print(f"  self-consistency: saturated rcm satisfies the adjustment fixed point (rcm={rcm:.2e})  PASS")

#!/usr/bin/env python3
"""test_rcm_sat_adj.py — validate the JAX rcm_sat_adj port (saturation.F90:rcm_sat_adj).

rcm_sat_adj iteratively (bisection) finds theta s.t. theta - (Lv/(Cp exner))·max(rtm - rsat(theta·exner),0) =
thlm, then returns rcm = max(rtm - rsat,0). The JAX port replicates the Fortran's exact "freeze theta once the
0.001 K tolerance is first hit" bisection, so it should bit-match. Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_rcm_sat_adj over saturated + unsaturated points, both saturation formulas.
     SKIPs if clubb_f2py is unbuilt.
  2. Self-consistency: the returned rcm satisfies the saturation-adjustment fixed point (residual ~ tolerance).
  3. A finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.saturation import (
    rcm_sat_adj, sat_mixrat_liq, SATURATION_FLATAU, SATURATION_BOLTON)

_Cp, _Lv = 1004.67, 2.5e6


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
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py rcm_sat_adj oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    n_sat = 0
    for formula in (SATURATION_FLATAU, SATURATION_BOLTON):
        for thlm, rtm, p, exner in _cases():
            ref = float(clubb_f2py.f2py_rcm_sat_adj(thlm, rtm, p, exner, formula))
            got = float(rcm_sat_adj(thlm, rtm, p, exner, formula))
            worst = max(worst, abs(got - ref))
            if ref > 1e-8:
                n_sat += 1
    assert worst < 1e-12, f"rcm_sat_adj f2py mismatch {worst:.2e}"
    assert n_sat > 0, "no saturated cases exercised"
    print(f"  f2py rcm_sat_adj: bit-match over saturated+unsaturated, both formulas, worst {worst:.2e}  PASS")


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


def test_differentiable():
    def loss(rtm):
        return rcm_sat_adj(285.0, rtm, 95000.0, 0.95, SATURATION_FLATAU) ** 2
    g = float(jax.grad(loss)(0.02))
    assert np.isfinite(g), "non-finite grad through rcm_sat_adj"
    print(f"  jax.grad through rcm_sat_adj: finite (d/drtm = {g:.3e})  PASS")


def main():
    print("test_rcm_sat_adj:")
    for t in (test_f2py_oracle, test_self_consistency, test_differentiable):
        t()
    print("All rcm_sat_adj checks PASSED")


if __name__ == "__main__":
    main()

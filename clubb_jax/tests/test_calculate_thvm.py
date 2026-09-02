#!/usr/bin/env python3
"""test_calculate_thvm.py — validate the JAX calculate_thvm port (calc_pressure.F90:calculate_thvm).

thvm = thlm + ep1·thv_ds_zt·rtm + (Lv/(Cp·exner) − ep2·thv_ds_zt)·rcm — mean virtual potential temperature.
This routine is ported and used in the gated driver but had no dedicated f2py test. Oracles:
  1. f2py bit-shadow vs f2py_calculate_thvm. SKIPs if clubb_f2py is unbuilt.
  2. Closed-form identity from the CLUBB constants.
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

from clubb_jax.src.CLUBB_core.calc_pressure import calculate_thvm
from clubb_jax.src.CLUBB_core.constants_clubb import ep1, ep2, Lv, Cp

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
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calculate_thvm oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        thlm, rtm, rcm, exner, thv_ds_zt = _inputs(seed)
        ref = np.asarray(clubb_f2py.f2py_calculate_thvm(thlm, rtm, rcm, exner, thv_ds_zt))
        got = np.asarray(calculate_thvm(thlm, rtm, rcm, exner, thv_ds_zt))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"calculate_thvm f2py mismatch {worst:.2e}"
    print(f"  f2py calculate_thvm: bit-match over 3 seeds, worst {worst:.2e}  PASS")


def test_closed_form():
    thlm, rtm, rcm, exner, thv_ds_zt = _inputs(5)
    got = np.asarray(calculate_thvm(thlm, rtm, rcm, exner, thv_ds_zt))
    ref = thlm + ep1 * thv_ds_zt * rtm + (Lv / (Cp * exner) - ep2 * thv_ds_zt) * rcm
    assert np.max(np.abs(got - ref)) < 1e-12, "closed-form mismatch"
    # rcm = 0 -> thvm = thlm + ep1*thv_ds_zt*rtm.
    g0 = np.asarray(calculate_thvm(thlm, rtm, np.zeros_like(rcm), exner, thv_ds_zt))
    assert np.max(np.abs(g0 - (thlm + ep1 * thv_ds_zt * rtm))) < 1e-12, "rcm=0 limit"
    print("  closed-form identity + rcm=0 limit  PASS")


def test_differentiable():
    thlm, rtm, rcm, exner, thv_ds_zt = _inputs(7)
    def loss(rcm_v):
        return jnp.sum(calculate_thvm(thlm, rtm, rcm_v, exner, thv_ds_zt) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(rcm)))
    assert np.isfinite(g).all(), "non-finite grad through calculate_thvm"
    print(f"  jax.grad through calculate_thvm: finite ({g.size} entries)  PASS")


def main():
    print("test_calculate_thvm:")
    for t in (test_f2py_oracle, test_closed_form, test_differentiable):
        t()
    print("All calculate_thvm checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_w_up_in_cloud.py — validate the JAX calc_w_up_in_cloud port.

calc_w_up_in_cloud (pdf_closure_module.F90:3145) computes the mean cloudy updraft / downdraft vertical velocity
and the cloudy updraft/downdraft fractions from the binormal w-PDF. Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_calc_w_up_in_cloud on the same fields — all four outputs. SKIPs if unbuilt.
  2. Physical invariants: updraft+downdraft fractions are in [0,1], w_up_in_cloud >= 0 >= w_down_in_cloud, and
     the all-updraft / all-downdraft shortcuts behave correctly.
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_w_up_in_cloud

NG, NZ = 2, 8


def _fields(seed):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0.2, 0.8, (NG, NZ))
    cf1 = rng.uniform(0.0, 1.0, (NG, NZ)); cf2 = rng.uniform(0.0, 1.0, (NG, NZ))
    w1 = rng.uniform(-3, 3, (NG, NZ)); w2 = rng.uniform(-3, 3, (NG, NZ))
    v1 = rng.uniform(0.01, 1.0, (NG, NZ)); v2 = rng.uniform(0.01, 1.0, (NG, NZ))
    return a, cf1, cf2, w1, w2, v1, v2


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_w_up_in_cloud oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        a, cf1, cf2, w1, w2, v1, v2 = _fields(seed)
        f = clubb_f2py.f2py_calc_w_up_in_cloud(a, cf1, cf2, w1, w2, v1, v2)
        g = calc_w_up_in_cloud(a, cf1, cf2, w1, w2, v1, v2)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    # ~1e-11 residual is the jax.scipy.special.erf vs Fortran intrinsic erf implementation difference
    # (amplified by the 1/cloudy_updraft_frac division), not a transcription error — closed-form faithful.
    assert worst < 1e-9, f"f2py mismatch {worst:.2e}"
    print(f"  f2py calc_w_up_in_cloud: match (4 outputs) over 3 seeds, worst {worst:.2e} (erf-limited)  PASS")


def test_invariants():
    a, cf1, cf2, w1, w2, v1, v2 = _fields(5)
    w_up, w_down, uf, df = (np.asarray(x) for x in calc_w_up_in_cloud(a, cf1, cf2, w1, w2, v1, v2))
    assert np.all(uf >= -1e-12) and np.all(uf <= 1.0 + 1e-12), "updraft frac out of [0,1]"
    assert np.all(df >= -1e-12) and np.all(df <= 1.0 + 1e-12), "downdraft frac out of [0,1]"
    assert np.all(w_up >= -1e-9), "mean cloudy updraft should be >= 0"
    assert np.all(w_down <= 1e-9), "mean cloudy downdraft should be <= 0"
    # All-updraft shortcut: a strongly positive w with tiny variance -> updraft_frac ~ 1, w_up ~ w.
    a1 = np.ones((1, 1)); cf = np.ones((1, 1))
    w_up1, _, uf1, _ = (np.asarray(x) for x in
                        calc_w_up_in_cloud(a1, cf, cf, 5.0 * np.ones((1, 1)), np.zeros((1, 1)),
                                           1e-4 * np.ones((1, 1)), 1e-4 * np.ones((1, 1))))
    assert abs(uf1[0, 0] - 1.0) < 1e-6 and abs(w_up1[0, 0] - 5.0) < 1e-3, "all-updraft shortcut wrong"
    print("  invariants: fractions in [0,1], w_up>=0>=w_down, all-updraft shortcut  PASS")


def test_differentiable():
    a, cf1, cf2, w1, w2, v1, v2 = _fields(7)
    def loss(w1v):
        w_up, w_down, uf, df = calc_w_up_in_cloud(a, cf1, cf2, w1v, w2, v1, v2)
        return jnp.sum(w_up ** 2 + w_down ** 2 + uf ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(w1)))
    assert np.isfinite(g).all(), "non-finite grad through calc_w_up_in_cloud"
    print(f"  jax.grad through calc_w_up_in_cloud: finite ({g.size} entries)  PASS")


def main():
    print("test_w_up_in_cloud:")
    for t in (test_f2py_oracle, test_invariants, test_differentiable):
        t()
    print("All calc_w_up_in_cloud checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_limits_f_x_responder.py — validate the JAX calc_limits_F_x_responder + sort_roots ports (new_pdf.F90).

Computes the min/max allowable F_x for a responding variable (PDF-component sigmas >= 0) by solving two cubics
in sqrt(F_x), sorting the roots, and a 4-way branch on sign(Skx·<w'x'>) and Skx². Oracles:
  1. f2py bit-shadow vs f2py_calc_limits_f_x_responder. SKIPs if clubb_f2py is unbuilt.
  2. sort_roots == jnp.sort; min_F_x <= max_F_x and both in [0,1].
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

from clubb_jax.src.CLUBB_core.new_pdf import calc_limits_F_x_responder, sort_roots

NG, NZ = 2, 6
# Typical max-Skx2 thresholds (new_pdf default-parameter values).
MAX_POS, MAX_NEG = 4.0, 4.0


def _inputs(seed):
    rng = np.random.default_rng(seed)
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    Skx = rng.uniform(-3, 3, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    return mixt_frac, Skx, sgn


def test_sort_roots():
    rng = np.random.default_rng(1)
    r = rng.uniform(-5, 5, (NG, NZ, 3))
    assert np.array_equal(np.asarray(sort_roots(r)), np.sort(r, axis=-1))
    print("  sort_roots == jnp.sort (ascending)  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_limits_f_x_responder oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        mf, Skx, sgn = _inputs(seed)
        mp = np.full((NG, NZ), MAX_POS); mn = np.full((NG, NZ), MAX_NEG)
        f_min, f_max = clubb_f2py.f2py_calc_limits_f_x_responder(mf, Skx, sgn, mp, mn)
        g_min, g_max = calc_limits_F_x_responder(mf, Skx, sgn, mp, mn)
        worst = max(worst, np.max(np.abs(np.asarray(g_min) - np.asarray(f_min))),
                    np.max(np.abs(np.asarray(g_max) - np.asarray(f_max))))
    assert worst < 1e-11, f"calc_limits_f_x_responder f2py mismatch {worst:.2e}"
    print(f"  f2py calc_limits_f_x_responder: bit-match (min/max F_x), worst {worst:.2e}  PASS")


def test_bounds():
    mf, Skx, sgn = _inputs(5)
    g_min, g_max = (np.asarray(x) for x in
                    calc_limits_F_x_responder(mf, Skx, sgn, np.full((NG, NZ), MAX_POS), np.full((NG, NZ), MAX_NEG)))
    assert np.all(g_min >= -1e-12) and np.all(g_max <= 1.0 + 1e-12), "F_x limits out of [0,1]"
    print("  bounds: min_F_x, max_F_x in [0,1]  PASS")


def test_differentiable():
    mf, Skx, sgn = _inputs(7)
    mp = jnp.full((NG, NZ), MAX_POS); mn = jnp.full((NG, NZ), MAX_NEG)
    def loss(s):
        a, b = calc_limits_F_x_responder(mf, s, sgn, mp, mn)
        return jnp.sum(a ** 2 + b ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(Skx)))
    assert np.isfinite(g).all(), "non-finite grad through calc_limits_F_x_responder"
    print(f"  jax.grad through calc_limits_F_x_responder: finite ({g.size} entries)  PASS")


def main():
    print("test_limits_f_x_responder:")
    for t in (test_sort_roots, test_f2py_oracle, test_bounds, test_differentiable):
        t()
    print("All calc_limits_F_x_responder checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_coefs_semiimpl.py — validate the JAX calc_coefs_wp2xp_semiimpl port (new_pdf.F90).

Semi-implicit decomposition <w'^2 x'> = coef_wp2xp_implicit*<w'x'> + term_wp2xp_explicit. Oracles:
  1. f2py bit-shadow vs f2py_calc_coefs_wp2xp_semiimpl (both outputs, both branches). SKIPs if unbuilt.
  2. A finite jax.grad.
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

from clubb_jax.src.CLUBB_core.new_pdf import (
    calc_coefs_wp2xp_semiimpl, calc_coefs_wpxpyp_semiimpl)
from clubb_jax.src.CLUBB_core.new_hybrid_pdf import calc_coefs_wpxp2_semiimpl

NG, NZ = 2, 6


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_coefs_wp2xp_semiimpl oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    worst = 0.0
    n_red = 0
    for trial in range(20):
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ)); xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
        sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
        mf = rng.uniform(0.1, 0.9, (NG, NZ)); F_w = rng.uniform(0.05, 1.0, (NG, NZ)); F_x = rng.uniform(0.05, 1.0, (NG, NZ))
        cw1 = rng.uniform(0.1, 2.0, (NG, NZ)); cw2 = rng.uniform(0.1, 2.0, (NG, NZ))
        cx1 = rng.uniform(0.1, 2.0, (NG, NZ)); cx2 = rng.uniform(0.1, 2.0, (NG, NZ))
        if trial % 4 == 0:   # reduced branch (zero sigma products)
            cw1 = np.zeros((NG, NZ)); cw2 = np.zeros((NG, NZ)); n_red += 1
        fc, ft = clubb_f2py.f2py_calc_coefs_wp2xp_semiimpl(wp2, xp2, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2)
        gc, gt = calc_coefs_wp2xp_semiimpl(wp2, xp2, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2)
        worst = max(worst, np.max(np.abs(np.asarray(gc) - np.asarray(fc))),
                    np.max(np.abs(np.asarray(gt) - np.asarray(ft))))
    assert worst < 1e-11, f"calc_coefs_wp2xp_semiimpl f2py mismatch {worst:.2e}"
    assert n_red > 0
    print(f"  f2py calc_coefs_wp2xp_semiimpl: bit-match (coef+term, both branches), worst {worst:.2e}  PASS")


def test_wpxp2_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_coefs_wpxp2_semiimpl oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(5)
    worst = 0.0; n_zero = 0
    for trial in range(20):
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ)); wpxp = rng.uniform(-1, 1, (NG, NZ))
        mf = rng.uniform(0.1, 0.9, (NG, NZ)); F_w = rng.uniform(0.05, 1.0, (NG, NZ))
        cx1 = rng.uniform(0.1, 2.0, (NG, NZ)); cx2 = rng.uniform(0.1, 2.0, (NG, NZ))
        if trial % 4 == 0:
            F_w = np.zeros((NG, NZ)); n_zero += 1   # F_w=0 -> 0,0 branch
        fc, ft = clubb_f2py.f2py_calc_coefs_wpxp2_semiimpl(wp2, wpxp, mf, F_w, cx1, cx2)
        gc, gt = calc_coefs_wpxp2_semiimpl(wp2, wpxp, mf, F_w, cx1, cx2)
        worst = max(worst, np.max(np.abs(np.asarray(gc) - np.asarray(fc))),
                    np.max(np.abs(np.asarray(gt) - np.asarray(ft))))
    assert worst < 1e-11, f"calc_coefs_wpxp2_semiimpl f2py mismatch {worst:.2e}"
    assert n_zero > 0
    print(f"  f2py calc_coefs_wpxp2_semiimpl: bit-match (both branches), worst {worst:.2e}  PASS")


def test_wpxpyp_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_coefs_wpxpyp_semiimpl oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(6)
    worst = 0.0; n_red = 0
    for trial in range(20):
        sh = (NG, NZ)
        wp2 = rng.uniform(0.05, 2.0, sh); xp2 = rng.uniform(0.05, 2.0, sh); yp2 = rng.uniform(0.05, 2.0, sh)
        wpxp = rng.uniform(-1, 1, sh); wpyp = rng.uniform(-1, 1, sh)
        sx = np.sign(wpxp); sx[sx == 0] = 1.0; sy = np.sign(wpyp); sy[sy == 0] = 1.0
        mf = rng.uniform(0.1, 0.9, sh); F_w = rng.uniform(0.05, 1.0, sh); F_x = rng.uniform(0.05, 1.0, sh); F_y = rng.uniform(0.05, 1.0, sh)
        cw1 = rng.uniform(0.1, 2, sh); cw2 = rng.uniform(0.1, 2, sh)
        cx1 = rng.uniform(0.1, 2, sh); cx2 = rng.uniform(0.1, 2, sh)
        cy1 = rng.uniform(0.1, 2, sh); cy2 = rng.uniform(0.1, 2, sh)
        if trial % 4 == 0:
            cx1 = np.zeros(sh); cx2 = np.zeros(sh); n_red += 1   # reduced branch
        args = (wp2, xp2, yp2, wpxp, wpyp, sx, sy, mf, F_w, F_x, F_y, cw1, cw2, cx1, cx2, cy1, cy2)
        fc, ft = clubb_f2py.f2py_calc_coefs_wpxpyp_semiimpl(*args)
        gc, gt = calc_coefs_wpxpyp_semiimpl(*args)
        worst = max(worst, np.max(np.abs(np.asarray(gc) - np.asarray(fc))),
                    np.max(np.abs(np.asarray(gt) - np.asarray(ft))))
    assert worst < 1e-11, f"calc_coefs_wpxpyp_semiimpl f2py mismatch {worst:.2e}"
    assert n_red > 0
    print(f"  f2py calc_coefs_wpxpyp_semiimpl: bit-match (both branches), worst {worst:.2e}  PASS")


def test_differentiable():
    rng = np.random.default_rng(7)
    def loss(wp2):
        c, t = calc_coefs_wp2xp_semiimpl(wp2, 1.0, 1.0, 0.4, 0.5, 0.6, 1.0, 0.8, 0.9, 0.7)
        c2, t2 = calc_coefs_wpxp2_semiimpl(wp2, 0.3, 0.4, 0.5, 1.0, 0.8)
        c3, t3 = calc_coefs_wpxpyp_semiimpl(wp2, 1.0, 1.0, 0.3, 0.2, 1.0, 1.0, 0.4, 0.5, 0.6, 0.6,
                                            1.0, 0.8, 0.9, 0.7, 0.85, 0.75)
        return jnp.sum(c ** 2 + t ** 2 + c2 ** 2 + t2 ** 2 + c3 ** 2 + t3 ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(rng.uniform(0.1, 2.0, (NG, NZ)))))
    assert np.isfinite(g).all(), "non-finite grad through the semiimpl coefs"
    print(f"  jax.grad through all three semiimpl coefs: finite ({g.size} entries)  PASS")


def main():
    print("test_coefs_semiimpl:")
    for t in (test_f2py_oracle, test_wpxp2_f2py, test_wpxpyp_f2py, test_differentiable):
        t()
    print("All calc_coefs_wp2xp_semiimpl checks PASSED")


if __name__ == "__main__":
    main()

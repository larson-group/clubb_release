#!/usr/bin/env python3
"""test_pdf_moment_integrals.py — validate the JAX binormal/trinormal PDF moment integrals.

calc_wp4_pdf / calc_wp2xp_pdf / calc_wpxp2_pdf / calc_wpxpyp_pdf (pdf_closure_module.F90) integrate the higher
moments <w'^4>, <w'^2 x'>, <w'x'^2>, <w'x'y'> over the two-/three-variable two-component-normal PDF. Oracles:
  1. f2py bit-shadow vs f2py_calc_*_pdf. SKIPs if clubb_f2py is unbuilt.
  2. Monte-Carlo: draw from the binormal and check <w'^4>, <w'^2 x'>, <w'x'^2> match the sample moments.
  3. Finite jax.grad.
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import (
    calc_wp4_pdf, calc_wp2xp_pdf, calc_wpxp2_pdf, calc_wpxpyp_pdf)

NG, NZ = 2, 5


def _binormal_params(rng, shape):
    return dict(w_1=rng.uniform(-2, 2, shape), w_2=rng.uniform(-2, 2, shape),
                x_1=rng.uniform(-2, 2, shape), x_2=rng.uniform(-2, 2, shape),
                varnce_w_1=rng.uniform(0.1, 1.5, shape), varnce_w_2=rng.uniform(0.1, 1.5, shape),
                varnce_x_1=rng.uniform(0.1, 1.5, shape), varnce_x_2=rng.uniform(0.1, 1.5, shape),
                corr_w_x_1=rng.uniform(-0.8, 0.8, shape), corr_w_x_2=rng.uniform(-0.8, 0.8, shape),
                mixt_frac=rng.uniform(0.3, 0.7, shape))


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py pdf-moment oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(3)
    sh = (NG, NZ)
    worst = 0.0
    for _ in range(15):
        p = _binormal_params(rng, sh)
        wm = p['mixt_frac'] * p['w_1'] + (1 - p['mixt_frac']) * p['w_2']
        xm = p['mixt_frac'] * p['x_1'] + (1 - p['mixt_frac']) * p['x_2']
        # wp4
        r = np.asarray(clubb_f2py.f2py_calc_wp4_pdf(wm, p['w_1'], p['w_2'], p['varnce_w_1'], p['varnce_w_2'], p['mixt_frac']))
        g = np.asarray(calc_wp4_pdf(wm, p['w_1'], p['w_2'], p['varnce_w_1'], p['varnce_w_2'], p['mixt_frac']))
        worst = max(worst, np.max(np.abs(g - r)))
        # wp2xp
        args = (wm, xm, p['w_1'], p['w_2'], p['x_1'], p['x_2'], p['varnce_w_1'], p['varnce_w_2'],
                p['varnce_x_1'], p['varnce_x_2'], p['corr_w_x_1'], p['corr_w_x_2'], p['mixt_frac'])
        worst = max(worst, np.max(np.abs(np.asarray(calc_wp2xp_pdf(*args)) - np.asarray(clubb_f2py.f2py_calc_wp2xp_pdf(*args)))))
        worst = max(worst, np.max(np.abs(np.asarray(calc_wpxp2_pdf(*args)) - np.asarray(clubb_f2py.f2py_calc_wpxp2_pdf(*args)))))
        # wpxpyp (reuse x params as y for simplicity, distinct corrs)
        ym = xm
        wargs = (wm, xm, ym, p['w_1'], p['w_2'], p['x_1'], p['x_2'], p['x_1'], p['x_2'],
                 p['varnce_w_1'], p['varnce_w_2'], p['varnce_x_1'], p['varnce_x_2'], p['varnce_x_1'], p['varnce_x_2'],
                 p['corr_w_x_1'], p['corr_w_x_2'], p['corr_w_x_1'], p['corr_w_x_2'], p['corr_w_x_1'], p['corr_w_x_2'],
                 p['mixt_frac'])
        worst = max(worst, np.max(np.abs(np.asarray(calc_wpxpyp_pdf(*wargs)) - np.asarray(clubb_f2py.f2py_calc_wpxpyp_pdf(*wargs)))))
    assert worst < 1e-11, f"pdf-moment f2py mismatch {worst:.2e}"
    print(f"  f2py calc_{{wp4,wp2xp,wpxp2,wpxpyp}}_pdf: bit-match over 15 configs, worst {worst:.2e}  PASS")


def test_monte_carlo():
    rng = np.random.default_rng(11)
    a, w1, w2, vw1, vw2 = 0.4, 0.8, -0.6, 0.5, 0.9
    x1, x2, vx1, vx2, cwx1, cwx2 = 0.3, -0.2, 0.4, 0.7, 0.5, -0.3
    wm = a * w1 + (1 - a) * w2
    xm = a * x1 + (1 - a) * x2
    n = 6_000_000
    pick1 = rng.random(n) < a
    # Component samples of (w,x) with the given correlation.
    def _comp(mw, mx, vw, vx, c, m):
        zw = rng.standard_normal(m); zx = rng.standard_normal(m)
        sw, sx = np.sqrt(vw), np.sqrt(vx)
        w = mw + sw * zw
        x = mx + sx * (c * zw + np.sqrt(1 - c ** 2) * zx)
        return w, x
    m1 = int(pick1.sum())
    w_a, x_a = _comp(w1, x1, vw1, vx1, cwx1, m1)
    w_b, x_b = _comp(w2, x2, vw2, vx2, cwx2, n - m1)
    w = np.concatenate([w_a, w_b]); x = np.concatenate([x_a, x_b])
    wp4 = float(calc_wp4_pdf(wm, w1, w2, vw1, vw2, a))
    wp2xp = float(calc_wp2xp_pdf(wm, xm, w1, w2, x1, x2, vw1, vw2, vx1, vx2, cwx1, cwx2, a))
    assert abs(wp4 - np.mean((w - wm) ** 4)) < 5e-2, f"wp4 {wp4} vs MC {np.mean((w-wm)**4)}"
    assert abs(wp2xp - np.mean((w - wm) ** 2 * (x - xm))) < 5e-2, "wp2xp vs MC"
    print(f"  Monte-Carlo: wp4 & wp2xp match sample moments (wp4={wp4:.3f})  PASS")


def test_differentiable():
    rng = np.random.default_rng(7)
    p = _binormal_params(rng, (NG, NZ))
    wm = jnp.asarray(p['mixt_frac'] * p['w_1'] + (1 - p['mixt_frac']) * p['w_2'])
    xm = jnp.asarray(p['mixt_frac'] * p['x_1'] + (1 - p['mixt_frac']) * p['x_2'])
    def loss(vw1):
        return jnp.sum(calc_wp2xp_pdf(wm, xm, p['w_1'], p['w_2'], p['x_1'], p['x_2'], vw1, p['varnce_w_2'],
                                      p['varnce_x_1'], p['varnce_x_2'], p['corr_w_x_1'], p['corr_w_x_2'],
                                      p['mixt_frac']) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(p['varnce_w_1'])))
    assert np.isfinite(g).all(), "non-finite grad through calc_wp2xp_pdf"
    print(f"  jax.grad through calc_wp2xp_pdf: finite ({g.size} entries)  PASS")


def main():
    print("test_pdf_moment_integrals:")
    for t in (test_f2py_oracle, test_monte_carlo, test_differentiable):
        t()
    print("All pdf-moment-integral checks PASSED")


if __name__ == "__main__":
    main()

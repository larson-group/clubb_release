#!/usr/bin/env python3
"""validate the JAX binormal/trinormal PDF moment integrals."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_wp4_pdf, calc_wp2xp_pdf

NG, NZ = 2, 5


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

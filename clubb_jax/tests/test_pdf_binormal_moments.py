#!/usr/bin/env python3
"""test_pdf_binormal_moments.py — validate the binormal mixed-moment integrals <w'^2 x'>, <w' x'^2>, <w'^2 x'^2>.

`calc_wp2xp_pdf`, `calc_wpxp2_pdf`, `calc_wp2xp2_pdf` (pdf_closure_module.py ↔ pdf_closure_module.F90) are the analytic
mixed central moments of the two-component BINORMAL PDF of (w, x) used by the ADG1 closure. They feed
calc_pdf_higher_order_moments_jax → pdf_closure_driver and were validated only end-to-end. This pins each (1) against
its closed-form transcription and (2) an INDEPENDENT Monte-Carlo — drawing (w, x) from the actual 2-component
bivariate-normal mixture (per-component means w_i/x_i, variances, correlation corr_i) and comparing the empirical
mixed central moment — so the analytic formulas (incl. the corr cross-terms and the (1+2corr²) coefficient in
<w'^2 x'^2>) are validated, not just transcribed. + finite grad. Companion to iter-559's calc_wp4_pdf. (iter 560)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.pdf_closure_module import (
    calc_wp2xp_pdf, calc_wpxp2_pdf, calc_wp2xp2_pdf)

# (w1,w2,x1,x2, vw1,vw2, vx1,vx2, c1,c2, a). O(1)-scaled so the mixed moments are well-conditioned for MC
# (the closed-form test covers the realistic small/large-mean regimes; MC needs O(1) moments for low noise).
_CASES = [
    (1.0, -0.5, 0.8, -0.6, 0.3, 0.5, 0.4, 0.6, 0.4, -0.3, 0.4),
    (1.4, -0.4, -0.7, 0.9, 0.25, 0.7, 0.5, 0.35, 0.6, 0.2, 0.55),
    (290.0, 292.0, 0.8, -0.3, 0.6, 0.9, 0.25, 0.6, 0.3, -0.4, 0.35),  # large-mean: closed-form only
]


def _refs(w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a):
    wm = a * w1 + (1 - a) * w2
    xm = a * x1 + (1 - a) * x2
    dw1, dw2, dx1, dx2 = w1 - wm, w2 - wm, x1 - xm, x2 - xm
    swx1, swx2 = (vw1 * vx1) ** 0.5, (vw2 * vx2) ** 0.5
    wp2xp = a * ((dw1 ** 2 + vw1) * dx1 + 2 * c1 * swx1 * dw1) + (1 - a) * ((dw2 ** 2 + vw2) * dx2 + 2 * c2 * swx2 * dw2)
    wpxp2 = a * (dw1 * (dx1 ** 2 + vx1) + 2 * c1 * swx1 * dx1) + (1 - a) * (dw2 * (dx2 ** 2 + vx2) + 2 * c2 * swx2 * dx2)
    wp2xp2 = (a * (dw1 ** 2 * (dx1 ** 2 + vx1) + 4 * c1 * swx1 * dx1 * dw1 + (dx1 ** 2 + (1 + 2 * c1 ** 2) * vx1) * vw1)
              + (1 - a) * (dw2 ** 2 * (dx2 ** 2 + vx2) + 4 * c2 * swx2 * dx2 * dw2 + (dx2 ** 2 + (1 + 2 * c2 ** 2) * vx2) * vw2))
    return wm, xm, wp2xp, wpxp2, wp2xp2


def test_closed_form():
    worst = 0.0
    for cs in _CASES:
        w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a = cs
        wm, xm, r_wp2xp, r_wpxp2, r_wp2xp2 = _refs(*cs)
        g_wp2xp = float(calc_wp2xp_pdf(wm, xm, w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a))
        g_wpxp2 = float(calc_wpxp2_pdf(wm, xm, w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a))
        g_wp2xp2 = float(calc_wp2xp2_pdf(wm, xm, w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a))
        for got, ref in ((g_wp2xp, r_wp2xp), (g_wpxp2, r_wpxp2), (g_wp2xp2, r_wp2xp2)):
            worst = max(worst, abs(got - ref) / (abs(ref) + 1e-300))
    assert worst < 1e-12, f"binormal moments vs closed form {worst:.2e}"
    print(f"  <w'²x'>, <w'x'²>, <w'²x'²> closed forms (worst {worst:.1e})  PASS")


def test_monte_carlo():
    rng = np.random.default_rng(560)
    N = 10_000_000
    worst = 0.0
    for cs in _CASES[:2]:   # O(1)-scaled cases only (large-mean case is closed-form-validated)
        w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a = cs
        wm, xm, r_wp2xp, r_wpxp2, r_wp2xp2 = _refs(*cs)
        pick = rng.random(N) < a
        # per-component bivariate normal via Cholesky: w = mu_w + sw*z1 ; x = mu_x + sx*(c*z1 + sqrt(1-c²)*z2)
        z1, z2 = rng.standard_normal(N), rng.standard_normal(N)
        sw = np.where(pick, np.sqrt(vw1), np.sqrt(vw2)); sx = np.where(pick, np.sqrt(vx1), np.sqrt(vx2))
        muw = np.where(pick, w1, w2); mux = np.where(pick, x1, x2); cc = np.where(pick, c1, c2)
        w = muw + sw * z1
        x = mux + sx * (cc * z1 + np.sqrt(1.0 - cc ** 2) * z2)
        dw, dx = w - wm, x - xm
        for mc, ref in ((np.mean(dw ** 2 * dx), r_wp2xp), (np.mean(dw * dx ** 2), r_wpxp2),
                        (np.mean(dw ** 2 * dx ** 2), r_wp2xp2)):
            rel = abs(ref - mc) / (abs(mc) + 1e-300)
            worst = max(worst, rel)
            assert rel < 8e-3, f"MC mismatch case {cs[:4]}: analytic {ref} vs MC {mc} rel {rel:.2e}"
    print(f"  vs {N//1_000_000}M-sample binormal mixture mixed central moments: rel <8e-3 (worst {worst:.1e})  PASS")


def test_grad_finite():
    cs = _CASES[0]; w1, w2, x1, x2, vw1, vw2, vx1, vx2, c1, c2, a = cs
    wm, xm = a * w1 + (1 - a) * w2, a * x1 + (1 - a) * x2
    g = jax.grad(lambda c: calc_wp2xp2_pdf(wm, xm, w1, w2, x1, x2, vw1, vw2, vx1, vx2, c, c2, a) ** 2)(c1)
    assert np.isfinite(float(g)), "non-finite grad wrt corr_w_x_1"
    print("  jax.grad(calc_wp2xp2_pdf) wrt corr finite  PASS")


def main():
    print("test_pdf_binormal_moments:")
    test_closed_form()
    test_monte_carlo()
    test_grad_finite()
    print("All binormal mixed-moment checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_calc_wp4_pdf.py — validate the PDF 4th-moment <w'^4> integral (pdf_closure decomposition).

`calc_wp4_pdf` (pdf_closure_module.py ↔ pdf_closure_module.F90:calc_wp4_pdf) is the analytic <w'^4> of the
two-component-normal PDF of w that the ADG1 closure uses:
    <w'^4> = Σ_i weight_i ( 3·σ_w_i^4 + 6·(μ_w_i−<w>)^2·σ_w_i^2 + (μ_w_i−<w>)^4 )
It is a core pdf_closure piece (feeds calc_pdf_higher_order_moments_jax → pdf_closure_driver) but was only validated
end-to-end. This pins it (1) against the closed-form transcription, and — crucially, an INDEPENDENT check — (2) a
Monte-Carlo: draw w from the actual 2-component Gaussian mixture and compare the empirical 4th central moment. The MC
validates the analytic moment formula itself (not a tautology). + finite grad. (iter 559)
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_wp4_pdf

# (w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac) — distinct components, the ADG1-style skewed mixture
_CASES = [
    (0.8, -0.3, 0.25, 0.6, 0.3),
    (1.5, -0.5, 0.1, 1.2, 0.45),
    (0.2, 0.2, 0.5, 0.5, 0.5),   # equal-mean symmetric case
]


def test_closed_form():
    worst = 0.0
    for w1, w2, v1, v2, a in _CASES:
        wm = a * w1 + (1.0 - a) * w2          # consistent mixture mean
        got = float(calc_wp4_pdf(wm, w1, w2, v1, v2, a))
        d1, d2 = w1 - wm, w2 - wm
        ref = a * (3 * v1 ** 2 + 6 * d1 ** 2 * v1 + d1 ** 4) + (1 - a) * (3 * v2 ** 2 + 6 * d2 ** 2 * v2 + d2 ** 4)
        worst = max(worst, abs(got - ref) / (abs(ref) + 1e-300))
    assert worst < 1e-13, f"calc_wp4_pdf vs closed form {worst:.2e}"
    print(f"  calc_wp4_pdf = Σ w_i(3σ⁴+6d²σ²+d⁴) closed form (worst {worst:.1e})  PASS")


def test_monte_carlo_central_moment():
    rng = np.random.default_rng(559)
    N = 8_000_000
    worst = 0.0
    for w1, w2, v1, v2, a in _CASES:
        wm = a * w1 + (1.0 - a) * w2
        got = float(calc_wp4_pdf(wm, w1, w2, v1, v2, a))
        # draw from the 2-component normal mixture: each sample from comp1 w.p. a, else comp2
        pick = rng.random(N) < a
        s = np.where(pick, rng.normal(w1, np.sqrt(v1), N), rng.normal(w2, np.sqrt(v2), N))
        mc = np.mean((s - wm) ** 4)
        rel = abs(got - mc) / (abs(mc) + 1e-300)
        worst = max(worst, rel)
        assert rel < 5e-3, f"wp4 case {(w1, w2, v1, v2, a)}: analytic {got} vs MC {mc} rel {rel:.2e}"
    print(f"  calc_wp4_pdf vs {N//1_000_000}M-sample mixture 4th central moment: rel <5e-3 (worst {worst:.1e})  PASS")


def test_grad_finite():
    g = jax.grad(lambda v1: calc_wp4_pdf(0.3, 0.8, -0.3, v1, 0.6, 0.3) ** 2)(0.25)
    assert np.isfinite(float(g)), "non-finite grad wrt varnce_w_1"
    print("  jax.grad(calc_wp4_pdf) wrt varnce_w_1 finite  PASS")


def main():
    print("test_calc_wp4_pdf:")
    test_closed_form()
    test_monte_carlo_central_moment()
    test_grad_finite()
    print("All calc_wp4_pdf checks PASSED")


if __name__ == "__main__":
    main()

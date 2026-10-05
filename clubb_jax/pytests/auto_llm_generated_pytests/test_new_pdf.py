#!/usr/bin/env python3
"""validate the JAX new_pdf.py ports (new-hybrid PDF helpers, Griffin & Larson 2018)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.new_pdf import calc_mixture_fraction
# calc_coef_wp2xp_implicit + the calculate_* aliases moved to new_hybrid_pdf.py (mirror-refactor iter 18)

NG, NZ = 2, 6


def _ref_mixt_frac(Skx, F, zeta, sgn):
    zp2 = zeta + 2.0
    if F > 0.0:
        sa = 4 * F ** 3 + 12 * F ** 2 * (1 - F) + 36 * F * (zeta + 1) * (1 - F) ** 2 / zp2 ** 2 + Skx ** 2
        num = (4 * F ** 3 + 18 * F * (zeta + 1) * (1 - F) / zp2 + 6 * F ** 2 * (1 - F) / zp2
               + Skx ** 2 - Skx * sgn * np.sqrt(sa))
        den = 2 * F * (F - 3) ** 2 + 2 * Skx ** 2
        return num / den
    return (zeta + 1) / zp2


def test_mixture_fraction():
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(200):
        F = rng.uniform(0.05, 1.0)
        Skx = rng.uniform(-3, 3)
        zeta = rng.uniform(0.0, 3.0)
        sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0)
        got = float(calc_mixture_fraction(Skx, F, zeta, sgn))
        ref = _ref_mixt_frac(Skx, F, zeta, sgn)
        worst = max(worst, abs(got - ref))
    assert worst < 1e-13, f"mixt_frac transcription mismatch {worst:.2e}"
    # Symmetric limit F=0, Skx=0 -> (zeta+1)/(zeta+2).
    for zeta in (0.0, 1.0, 2.5):
        v = float(calc_mixture_fraction(0.0, 0.0, zeta, 1.0))
        assert abs(v - (zeta + 1) / (zeta + 2)) < 1e-14, "symmetric limit wrong"
    print(f"  calc_mixture_fraction: literal transcription + symmetric limit, worst {worst:.2e}  PASS")

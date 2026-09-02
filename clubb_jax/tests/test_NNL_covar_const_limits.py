#!/usr/bin/env python3
"""test_NNL_covar_const_limits.py — validate the remaining trivar covariance consts as σ_x3→0 limits.

The covariance const variants `trivar_NNL_covar_const_{x3,x1x3,x2x3}` (PDF_integrals_covar.py) — the x3-deterministic
forms the dispatch selects — have no all-mixed-moment counterpart (so the iter-575/576 MM cross-check doesn't reach
them). But because x3 is LOGNORMAL, its σ_x3_n→0 limit is numerically clean (no parabolic-cylinder s_c overflow), so
each is exactly the σ_x3_n→0 limit of an already-validated routine:
    const_x3   == general covariance   at σ_x3_n→0   (general validated iter 575)
    const_x1x3 == const_x1             at σ_x3_n→0   (const_x1 validated iter 576)
    const_x2x3 == const_x2             at σ_x3_n→0   (const_x2 validated iter 576)
with the x3 mean passed as the deterministic value exp(μ_x3_n). This pins all three over random cases, completing the
trivariate covariance const-variant coverage (const_x1/x2/x1x2 via MM, 576; const_x3/x1x3/x2x3 via this limit, 577).
Oracle-independent; never SKIPs. (iter 577)
"""
import os
import sys
import math

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_covar import (
    trivar_NNL_covar, trivar_NNL_covar_const_x1, trivar_NNL_covar_const_x2,
    trivar_NNL_covar_const_x3, trivar_NNL_covar_const_x1x3, trivar_NNL_covar_const_x2x3)

_TINY = 1.0e-10


def test_const_x3_x1x3_x2x3_limits():
    rng = np.random.default_rng(577)
    w3 = w13 = w23 = 0.0
    for _ in range(150):
        mu_x1, mu_x2, mu_x3_n = rng.uniform(-.5, .8), rng.uniform(.1, .6), rng.uniform(-9, -6)
        s_x1, s_x2 = rng.uniform(.2, .7), rng.uniform(.2, .6)
        r12, r13n, r23n = rng.uniform(-.3, .3), rng.uniform(-.2, .2), rng.uniform(-.2, .2)
        x1m, Ym, al, be = rng.uniform(.3, .6), 0.0, rng.uniform(.5, 2.), rng.uniform(.3, 1.2)
        mu_x3 = math.exp(mu_x3_n)
        # const_x3 == general(σ_x3_n→0)
        cx3 = float(trivar_NNL_covar_const_x3(mu_x1, mu_x2, mu_x3, s_x1, s_x2, r12, x1m, Ym, al, be))
        gen = float(trivar_NNL_covar(mu_x1, mu_x2, mu_x3_n, s_x1, s_x2, _TINY, r12, r13n, r23n, x1m, Ym, al, be))
        w3 = max(w3, abs(cx3 - gen) / (abs(gen) + 1e-300))
        # const_x1x3 == const_x1(σ_x3_n→0)
        cx13 = float(trivar_NNL_covar_const_x1x3(mu_x1, mu_x2, mu_x3, s_x2, x1m, Ym, al, be))
        cx1 = float(trivar_NNL_covar_const_x1(mu_x1, mu_x2, mu_x3_n, s_x2, _TINY, r23n, x1m, Ym, al, be))
        w13 = max(w13, abs(cx13 - cx1) / (abs(cx1) + 1e-300))
        # const_x2x3 == const_x2(σ_x3_n→0)
        cx23 = float(trivar_NNL_covar_const_x2x3(mu_x1, mu_x2, mu_x3, x1m, Ym, al, be))
        cx2 = float(trivar_NNL_covar_const_x2(mu_x1, mu_x2, mu_x3_n, s_x1, _TINY, r13n, x1m, Ym, al, be))
        w23 = max(w23, abs(cx23 - cx2) / (abs(cx2) + 1e-300))
    assert w3 < 1e-7 and w13 < 1e-7 and w23 < 1e-7, f"limits: x3 {w3:.1e}, x1x3 {w13:.1e}, x2x3 {w23:.1e}"
    print(f"  const_x3==general, const_x1x3==const_x1, const_x2x3==const_x2 at σ_x3→0 "
          f"(worst {max(w3, w13, w23):.1e})  PASS")


def main():
    print("test_NNL_covar_const_limits:")
    test_const_x3_x1x3_x2x3_limits()
    print("All trivar covar const-limit checks PASSED")


if __name__ == "__main__":
    main()

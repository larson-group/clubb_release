#!/usr/bin/env python3
"""test_quadrivar_covar_const_limits.py — validate quadrivar covar consts via σ_x3/x4→0 limits + the all-const zero.

Extends the iter-577 limit approach to the quadrivariate covariance const variants whose x3/x4 (LOGNORMAL) directions
are degenerate — their σ_x3_n/σ_x4_n→0 limits are numerically clean (no parabolic-cylinder overflow), so each equals
the limit of an already-validated routine:
    quadrivar_NNLL_covar_const_x3   == general              at σ_x3_n→0          (general validated iter 575)
    quadrivar_NNLL_covar_const_x3x4 == general              at σ_x3_n,σ_x4_n→0
    quadrivar_NNLL_covar_const_all  == 0                    (all variables deterministic ⇒ Cov(const, Y) = 0)
with x3/x4 means passed as exp(μ_x3_n)/exp(μ_x4_n). This pins those three over random cases. The other quadrivar
consts (const_x1/x2/x1x2 vs MM, iter 576; the swap-reused forms via the dispatch test, iter 554) are already covered.
Oracle-independent; never SKIPs. (iter 578)
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
    quadrivar_NNLL_covar, quadrivar_NNLL_covar_const_x3,
    quadrivar_NNLL_covar_const_x3x4, quadrivar_NNLL_covar_const_all)

_T = 1.0e-10


def test_const_x3_x3x4_and_all():
    rng = np.random.default_rng(578)
    w3 = w34 = 0.0
    worst_all = 0.0
    for _ in range(150):
        mx1, mx2, mx3n, mx4n = rng.uniform(-.5, .8), rng.uniform(.1, .6), rng.uniform(-9, -6), rng.uniform(13, 15)
        s1, s2, s3n, s4n = rng.uniform(.2, .7), rng.uniform(.2, .6), rng.uniform(.3, .6), rng.uniform(.4, .7)
        r12, r13n, r14n = rng.uniform(-.2, .2), rng.uniform(-.15, .15), rng.uniform(-.15, .15)
        r23n, r24n, r34n = rng.uniform(-.15, .15), rng.uniform(-.1, .1), rng.uniform(-.1, .1)
        x1m, Ym, a, b, g = rng.uniform(.3, .6), 0.0, rng.uniform(.5, 1.5), rng.uniform(.3, 1.), rng.uniform(.3, 1.)
        mx3, mx4 = math.exp(mx3n), math.exp(mx4n)
        # const_x3 == general(σ_x3_n→0)
        cx3 = float(quadrivar_NNLL_covar_const_x3(mx1, mx2, mx3, mx4n, s1, s2, s4n, r12, r14n, r24n, x1m, Ym, a, b, g))
        g3 = float(quadrivar_NNLL_covar(mx1, mx2, mx3n, mx4n, s1, s2, _T, s4n,
                                        r12, r13n, r14n, r23n, r24n, r34n, x1m, Ym, a, b, g))
        w3 = max(w3, abs(cx3 - g3) / (abs(g3) + 1e-300))
        # const_x3x4 == general(σ_x3_n,σ_x4_n→0)
        cx34 = float(quadrivar_NNLL_covar_const_x3x4(mx1, mx2, mx3, mx4, s1, s2, r12, x1m, Ym, a, b, g))
        g34 = float(quadrivar_NNLL_covar(mx1, mx2, mx3n, mx4n, s1, s2, _T, _T,
                                         r12, r13n, r14n, r23n, r24n, r34n, x1m, Ym, a, b, g))
        w34 = max(w34, abs(cx34 - g34) / (abs(g34) + 1e-300))
        # const_all: all deterministic ⇒ the kernel (covariation) is 0, leaving only the linear <Y> term
        #   const_all == −<Y>·(mu_x1 − x1_mean)   (cf. the general covar's <Y> subtraction)
        Yv = rng.uniform(0, 1e-3)
        call = float(quadrivar_NNLL_covar_const_all(mx1, mx2, mx3, mx4, x1m, Yv, a, b, g))
        worst_all = max(worst_all, abs(call - (-Yv * (mx1 - x1m))))
    assert w3 < 1e-7 and w34 < 1e-7, f"limits: x3 {w3:.1e}, x3x4 {w34:.1e}"
    assert worst_all < 1e-15, f"const_all != −<Y>·(mu_x1−x1_mean): {worst_all:.2e}"
    print(f"  const_x3==general, const_x3x4==general at σ_x3/x4→0 (worst {max(w3,w34):.1e}); "
          f"const_all = −<Y>(μ_x1−x1m)  PASS")


def main():
    print("test_quadrivar_covar_const_limits:")
    test_const_x3_x3x4_and_all()
    print("All quadrivar covar const-limit checks PASSED")


if __name__ == "__main__":
    main()

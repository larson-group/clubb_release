#!/usr/bin/env python3
"""test_NNL_covar_const_vs_MM.py — validate the KK covariance CONST-variant primitives vs the MC-tested MM consts.

Extends iter-575's covar-vs-MM kernel identity to the (near-)zero-variance CONST variants that the covariance dispatch
selects. The all-mixed-moment const variants `trivar_NNL_MM_const_{x1,x2,x1x2}` / `quadrivar_NNLL_MM_const_{x1,x2,x1x2}`
(PDF_integrals_all_MM.py) are Monte-Carlo-validated (test_pdf_integrals_all_mm.py), and share the identical kernel with
the corresponding covariance const variant for a=b=1, differing only in the trailing <Y> subtraction — so at <Y>=0,
`covar_const_X == MM_const_X(a=1,b=1)`. This pins `trivar_NNL_covar_const_{x1,x2,x1x2}` and
`quadrivar_NNLL_covar_const_{x1,x2,x1x2}` (the x1-overflow / x2-degenerate const forms the dispatch uses) against those
MC-tested MM consts. The remaining covar const variants (const_x3 / const_x1x3 / const_x2x3, …) have no MM counterpart
and are covered by the dispatch-wiring tests (iters 553/554) + end-to-end. (iter 576)
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

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_covar import (
    trivar_NNL_covar_const_x1, trivar_NNL_covar_const_x2, trivar_NNL_covar_const_x1x2,
    quadrivar_NNLL_covar_const_x1, quadrivar_NNLL_covar_const_x2, quadrivar_NNLL_covar_const_x1x2)
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_all_MM import (
    trivar_NNL_MM_const_x1, trivar_NNL_MM_const_x2, trivar_NNL_MM_const_x1x2,
    quadrivar_NNLL_MM_const_x1, quadrivar_NNLL_MM_const_x2, quadrivar_NNLL_MM_const_x1x2)


def _xcheck(covar_fn, mm_fn, argf, n=150, seed=0):
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n):
        a = argf(rng)
        cov = float(covar_fn(*a)); mm = float(mm_fn(*a, 1, 1))
        worst = max(worst, abs(cov - mm) / (abs(mm) + 1e-300))
    return worst


def test_trivar_const_variants():
    # const_x1: (mu_x1,mu_x2,mu_x3_n, sigma_x2,sigma_x3_n,rho_x2x3_n, x1_mean, Y_mean=0, alpha,beta)
    f1 = lambda r: [r.uniform(-.5,.8), r.uniform(.1,.6), r.uniform(-9,-6), r.uniform(.2,.6), r.uniform(.3,.6),
                    r.uniform(-.2,.2), r.uniform(.3,.6), 0.0, r.uniform(.5,2.), r.uniform(.3,1.2)]
    # const_x2: (mu_x1,mu_x2,mu_x3_n, sigma_x1,sigma_x3_n,rho_x1x3_n, x1_mean, Y_mean=0, alpha,beta)
    f2 = lambda r: [r.uniform(-.5,.8), r.uniform(.1,.6), r.uniform(-9,-6), r.uniform(.2,.7), r.uniform(.3,.6),
                    r.uniform(-.2,.2), r.uniform(.3,.6), 0.0, r.uniform(.5,2.), r.uniform(.3,1.2)]
    # const_x1x2: (mu_x1,mu_x2,mu_x3_n, sigma_x3_n, x1_mean, Y_mean=0, alpha,beta)
    f12 = lambda r: [r.uniform(-.5,.8), r.uniform(.1,.6), r.uniform(-9,-6), r.uniform(.3,.6), r.uniform(.3,.6),
                     0.0, r.uniform(.5,2.), r.uniform(.3,1.2)]
    w1 = _xcheck(trivar_NNL_covar_const_x1, trivar_NNL_MM_const_x1, f1, seed=1)
    w2 = _xcheck(trivar_NNL_covar_const_x2, trivar_NNL_MM_const_x2, f2, seed=2)
    w12 = _xcheck(trivar_NNL_covar_const_x1x2, trivar_NNL_MM_const_x1x2, f12, seed=3)
    assert w1 < 1e-11 and w2 < 1e-11 and w12 < 1e-11, f"trivar const: {w1:.1e}/{w2:.1e}/{w12:.1e}"
    print(f"  trivar covar const_x1/x2/x1x2 == MM_const(a=1,b=1) at <Y>=0 (worst {max(w1,w2,w12):.1e})  PASS")


def test_quadrivar_const_variants():
    # const_x1: (mu_x1,mu_x2,mu_x3_n,mu_x4_n, sigma_x2,sigma_x3_n,sigma_x4_n, rho_x2x3_n,rho_x2x4_n,rho_x3x4_n, x1_mean,Y_mean=0, a,b,g)
    f1 = lambda r: [r.uniform(-.5,.8), r.uniform(.1,.6), r.uniform(-9,-6), r.uniform(13,15),
                    r.uniform(.2,.6), r.uniform(.3,.6), r.uniform(.4,.7),
                    r.uniform(-.15,.15), r.uniform(-.1,.1), r.uniform(-.1,.1),
                    r.uniform(.3,.6), 0.0, r.uniform(.5,1.5), r.uniform(.3,1.), r.uniform(.3,1.)]
    # const_x2: (mu_x1,mu_x2,mu_x3_n,mu_x4_n, sigma_x1,sigma_x3_n,sigma_x4_n, rho_x1x3_n,rho_x1x4_n,rho_x3x4_n, x1_mean,Y_mean=0, a,b,g)
    f2 = lambda r: [r.uniform(-.5,.8), r.uniform(.1,.6), r.uniform(-9,-6), r.uniform(13,15),
                    r.uniform(.2,.7), r.uniform(.3,.6), r.uniform(.4,.7),
                    r.uniform(-.15,.15), r.uniform(-.1,.1), r.uniform(-.1,.1),
                    r.uniform(.3,.6), 0.0, r.uniform(.5,1.5), r.uniform(.3,1.), r.uniform(.3,1.)]
    # const_x1x2: (mu_x1,mu_x2,mu_x3_n,mu_x4_n, sigma_x3_n,sigma_x4_n,rho_x3x4_n, x1_mean,Y_mean=0, a,b,g)
    f12 = lambda r: [r.uniform(-.5,.8), r.uniform(.1,.6), r.uniform(-9,-6), r.uniform(13,15),
                     r.uniform(.3,.6), r.uniform(.4,.7), r.uniform(-.1,.1),
                     r.uniform(.3,.6), 0.0, r.uniform(.5,1.5), r.uniform(.3,1.), r.uniform(.3,1.)]
    w1 = _xcheck(quadrivar_NNLL_covar_const_x1, quadrivar_NNLL_MM_const_x1, f1, seed=4)
    w2 = _xcheck(quadrivar_NNLL_covar_const_x2, quadrivar_NNLL_MM_const_x2, f2, seed=5)
    w12 = _xcheck(quadrivar_NNLL_covar_const_x1x2, quadrivar_NNLL_MM_const_x1x2, f12, seed=6)
    assert w1 < 1e-11 and w2 < 1e-11 and w12 < 1e-11, f"quadrivar const: {w1:.1e}/{w2:.1e}/{w12:.1e}"
    print(f"  quadrivar covar const_x1/x2/x1x2 == MM_const(a=1,b=1) at <Y>=0 (worst {max(w1,w2,w12):.1e})  PASS")


def main():
    print("test_NNL_covar_const_vs_MM:")
    test_trivar_const_variants()
    test_quadrivar_const_variants()
    print("All NNL covar-const-vs-MM cross-checks PASSED")


if __name__ == "__main__":
    main()

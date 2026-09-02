#!/usr/bin/env python3
"""test_NNL_covar_vs_MM.py — validate the KK covariance integrals against the MC-tested all-mixed-moment integrals.

The KK covariance primitives `trivar_NNL_covar` / `quadrivar_NNLL_covar` (PDF_integrals_covar.py) compute
Cov_i(x1, x2^α x3^β [x4^γ]) for x1,x2 normal (x2>0), x3[,x4] lognormal — parabolic-cylinder-based, previously
isolation-untested (deferred as "needs truncated-covariance Monte-Carlo"). The all-mixed-moment integrals
`trivar_NNL_MM`/`quadrivar_NNLL_MM` (PDF_integrals_all_MM.py) compute `<(x1−<x1>)^a (x2^α x3^β[ x4^γ]−<Y>)^b>` and are
Monte-Carlo-validated (test_pdf_integrals_all_mm.py). With the SAME args and a=b=1 they share the identical
parabolic-cylinder kernel; they differ ONLY in the trailing <Y>-subtraction term — the covariance subtracts
`<Y>·(mu_x1−x1_mean)`, the MM subtracts `<Y>·<x1−<x1>>_trunc` — so with `<Y>=0` the two coincide EXACTLY (verified to
~1e-16). This pins the covar primitives two ways, no fresh MC: (1) at <Y>=0, `covar == MM(a=1,b=1)` (cross-check
against the MC-tested kernel); (2) the <Y> dependence: `covar(<Y>) − covar(0) == −<Y>·(mu_x1−x1_mean)` (linear term).
(iter 575)
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

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_covar import trivar_NNL_covar, quadrivar_NNLL_covar
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_all_MM import trivar_NNL_MM, quadrivar_NNLL_MM


def _tri_args(rng):
    return [rng.uniform(-0.5, 0.8), rng.uniform(0.1, 0.6), rng.uniform(-9.0, -6.0),
            rng.uniform(0.2, 0.7), rng.uniform(0.2, 0.6), rng.uniform(0.3, 0.6),
            rng.uniform(-0.3, 0.3), rng.uniform(-0.2, 0.2), rng.uniform(-0.2, 0.2),
            rng.uniform(0.3, 0.6), 0.0, rng.uniform(0.5, 2.0), rng.uniform(0.3, 1.2)]   # x1_mean, Y_mean(=0), alpha, beta


def test_trivar_covar_kernel_and_Ymean_term():
    rng = np.random.default_rng(575)
    worst_k = 0.0; worst_y = 0.0
    for _ in range(200):
        a = _tri_args(rng)
        # (1) <Y>=0 cross-check vs the MC-tested MM kernel
        cov0 = float(trivar_NNL_covar(*a)); mm = float(trivar_NNL_MM(*a, 1, 1))
        worst_k = max(worst_k, abs(cov0 - mm) / (abs(mm) + 1e-300))
        # (2) the <Y> linear-subtraction term: covar(Y) − covar(0) == −Y·(mu_x1 − x1_mean)
        Y = rng.uniform(1e-4, 5e-4); aY = list(a); aY[10] = Y
        dcov = float(trivar_NNL_covar(*aY)) - cov0
        worst_y = max(worst_y, abs(dcov - (-Y * (a[0] - a[9]))))
    assert worst_k < 1e-11, f"trivar covar kernel != MM(a=1,b=1) at <Y>=0: {worst_k:.2e}"
    assert worst_y < 1e-12, f"trivar covar <Y> term wrong: {worst_y:.2e}"
    print(f"  trivar_NNL_covar: kernel==MM(1,1) at <Y>=0 (rel {worst_k:.1e}); <Y> term = −<Y>(μ_x1−x1m) ({worst_y:.1e})  PASS")


def _quad_args(rng):
    return [rng.uniform(-0.5, 0.8), rng.uniform(0.1, 0.6), rng.uniform(-9.0, -6.0), rng.uniform(13.0, 15.0),
            rng.uniform(0.2, 0.7), rng.uniform(0.2, 0.6), rng.uniform(0.3, 0.6), rng.uniform(0.4, 0.7),
            rng.uniform(-0.2, 0.2), rng.uniform(-0.15, 0.15), rng.uniform(-0.15, 0.15),
            rng.uniform(-0.15, 0.15), rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1),
            rng.uniform(0.3, 0.6), 0.0, rng.uniform(0.5, 1.5), rng.uniform(0.3, 1.0), rng.uniform(0.3, 1.0)]


def test_quadrivar_covar_kernel_and_Ymean_term():
    rng = np.random.default_rng(5751)
    worst_k = 0.0; worst_y = 0.0
    for _ in range(200):
        a = _quad_args(rng)
        cov0 = float(quadrivar_NNLL_covar(*a)); mm = float(quadrivar_NNLL_MM(*a, 1, 1))
        worst_k = max(worst_k, abs(cov0 - mm) / (abs(mm) + 1e-300))
        Y = rng.uniform(1e-4, 5e-4); aY = list(a); aY[15] = Y
        dcov = float(quadrivar_NNLL_covar(*aY)) - cov0
        worst_y = max(worst_y, abs(dcov - (-Y * (a[0] - a[14]))))
    assert worst_k < 1e-11, f"quadrivar covar kernel != MM(a=1,b=1) at <Y>=0: {worst_k:.2e}"
    assert worst_y < 1e-12, f"quadrivar covar <Y> term wrong: {worst_y:.2e}"
    print(f"  quadrivar_NNLL_covar: kernel==MM(1,1) at <Y>=0 (rel {worst_k:.1e}); <Y> term = −<Y>(μ_x1−x1m) ({worst_y:.1e})  PASS")


def main():
    print("test_NNL_covar_vs_MM:")
    test_trivar_covar_kernel_and_Ymean_term()
    test_quadrivar_covar_kernel_and_Ymean_term()
    print("All NNL covar-vs-MM cross-checks PASSED")


if __name__ == "__main__":
    main()

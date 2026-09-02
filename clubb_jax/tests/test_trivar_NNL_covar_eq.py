#!/usr/bin/env python3
"""test_trivar_NNL_covar_eq.py — pin the variance-regime DISPATCH of the KK trivariate NNL covariance integral.

`trivar_NNL_covar_eq` (KK_upscaled_covariances.py ↔ KK_upscaled_covariances.F90) computes Cov_i(x, chi^α y^β)
(x=w|eta normal, chi normal, y=Ncn|rr lognormal) by selecting among the trivariate covariance primitives per which
of σ_x1/σ_x2/σ_x3 are (near-)zero — a 7-way ladder with one subtlety: the all-const case maps to `const_x2x3`
(documented equivalence const_all == const_x2x3 for the covariance). The covar primitives (PDF_integrals_covar) are
parabolic-cylinder-based and validated only end-to-end (test_kk_rico_oracle), but the DISPATCH wiring — which
primitive + arg subset (incl. the σ_x2 denominator guard `s2g` and `x2a = mc_tndcy_mean/mc_coef`) per regime — was
never pinned. Companion to iter-549/550's mean dispatches: drives each regime and asserts `trivar_NNL_covar_eq`
returns EXACTLY the corresponding primitive with the documented args, so a mis-wired branch / wrong-arg / wrong
all-const-equivalence bug is caught directly. Oracle-independent; never SKIPs. (iter 553)
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

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_covariances import trivar_NNL_covar_eq, _CHI_TOL
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_covar import (
    trivar_NNL_covar, trivar_NNL_covar_const_x1, trivar_NNL_covar_const_x1x2,
    trivar_NNL_covar_const_x1x3, trivar_NNL_covar_const_x2, trivar_NNL_covar_const_x2x3,
    trivar_NNL_covar_const_x3)

# nominal moments (x = w/eta, chi, y = Ncn/rr). beta ≥ 0 so mu_x3 = mu_y.
_MU_X, _MU_CHI, _MU_Y, _MU_Y_N = 0.5, 0.3, 1.0e-3, -7.0
_S_Y_N, _RXC, _RXY_N, _RCY_N = 0.5, 0.2, 0.15, 0.1
_X_MEAN, _MC_TND, _MC_COEF = 0.4, 2.0, 1.5
_X_TOL, _Y_TOL, _A, _B = 1.0e-8, 1.0e-7, 1.0, 0.9
_X2A = _MC_TND / _MC_COEF


def _eq(s_x, s_chi, s_y):
    return float(trivar_NNL_covar_eq(_MU_X, _MU_CHI, _MU_Y, _MU_Y_N, s_x, s_chi, s_y, _S_Y_N,
                                     _RXC, _RXY_N, _RCY_N, _X_MEAN, _MC_TND, _MC_COEF,
                                     _X_TOL, _Y_TOL, _A, _B))


def _s2g(s_chi):
    return s_chi if s_chi > _CHI_TOL else 1.0


def _close(got, ref, name):
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"{name}: {got} vs {ref}"


def test_all_vary_base():
    _close(_eq(0.4, 0.3, 2e-4),
           float(trivar_NNL_covar(_MU_X, _MU_CHI, _MU_Y_N, 0.4, _s2g(0.3), _S_Y_N,
                                  _RXC, _RXY_N, _RCY_N, _X_MEAN, _X2A, _A, _B)), "base")
    print("  all vary -> trivar_NNL_covar (base)  PASS")


def test_c3_const_x3():
    _close(_eq(0.4, 0.3, 0.0),
           float(trivar_NNL_covar_const_x3(_MU_X, _MU_CHI, _MU_Y, 0.4, _s2g(0.3), _RXC, _X_MEAN, _X2A, _A, _B)),
           "const_x3")
    print("  sigma_y≤y_tol -> const_x3  PASS")


def test_c2_const_x2():
    _close(_eq(0.4, 0.0, 2e-4),
           float(trivar_NNL_covar_const_x2(_MU_X, _MU_CHI, _MU_Y_N, 0.4, _S_Y_N, _RXY_N, _X_MEAN, _X2A, _A, _B)),
           "const_x2")
    print("  sigma_chi≤CHI_TOL -> const_x2  PASS")


def test_c1_const_x1():
    _close(_eq(0.0, 0.3, 2e-4),
           float(trivar_NNL_covar_const_x1(_MU_X, _MU_CHI, _MU_Y_N, _s2g(0.3), _S_Y_N, _RCY_N, _X_MEAN, _X2A, _A, _B)),
           "const_x1")
    print("  sigma_x≤x_tol -> const_x1  PASS")


def test_c2c3_const_x2x3():
    _close(_eq(0.4, 0.0, 0.0),
           float(trivar_NNL_covar_const_x2x3(_MU_X, _MU_CHI, _MU_Y, _X_MEAN, _X2A, _A, _B)), "const_x2x3")
    print("  chi&y const -> const_x2x3  PASS")


def test_c1c3_const_x1x3():
    _close(_eq(0.0, 0.3, 0.0),
           float(trivar_NNL_covar_const_x1x3(_MU_X, _MU_CHI, _MU_Y, _s2g(0.3), _X_MEAN, _X2A, _A, _B)), "const_x1x3")
    print("  x&y const -> const_x1x3  PASS")


def test_c1c2_const_x1x2():
    _close(_eq(0.0, 0.0, 2e-4),
           float(trivar_NNL_covar_const_x1x2(_MU_X, _MU_CHI, _MU_Y_N, _S_Y_N, _X_MEAN, _X2A, _A, _B)), "const_x1x2")
    print("  x&chi const -> const_x1x2  PASS")


def test_all_const_maps_to_x2x3():
    # all three σ ≤ tol -> the documented equivalence const_all == const_x2x3.
    _close(_eq(0.0, 0.0, 0.0),
           float(trivar_NNL_covar_const_x2x3(_MU_X, _MU_CHI, _MU_Y, _X_MEAN, _X2A, _A, _B)), "all-const==x2x3")
    print("  all σ≤tol -> const_x2x3 (const_all equivalence)  PASS")


def main():
    print("test_trivar_NNL_covar_eq:")
    for t in (test_all_vary_base, test_c3_const_x3, test_c2_const_x2, test_c1_const_x1,
              test_c2c3_const_x2x3, test_c1c3_const_x1x3, test_c1c2_const_x1x2, test_all_const_maps_to_x2x3):
        t()
    print("All trivar_NNL_covar_eq dispatch checks PASSED")


if __name__ == "__main__":
    main()

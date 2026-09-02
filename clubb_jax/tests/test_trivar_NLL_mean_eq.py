#!/usr/bin/env python3
"""test_trivar_NLL_mean_eq.py — pin the 8-way variance DISPATCH of the KK trivariate normal-lognormal mean integral.

`trivar_NLL_mean_eq` (KK_upscaled_means.py ↔ KK_upscaled_means.F90:586) computes <chi^alpha · r_r^beta · N_r^gamma>_i
by selecting among the trivariate PDF-integral primitives per the (near-)zero-variance pattern of x1=chi, x2=r_r,
x3=N_r — an 8-way if/elseif ladder. The subtle part is ARG-SWAPPING: the chi+Nr-const case reuses `const_x1x2` with
(x2,x3) swapped, and the Nr-only-const case reuses `const_x2` with (x2,x3) swapped. The primitives are validated by
test_pdf_integrals_all_mm.py and the KK rate path is bit-exact vs rico (test_kk_rico_oracle.py), but this dispatch
wiring (incl. the swaps) was never pinned. This drives each of the 8 regimes and asserts `trivar_NLL_mean_eq` returns
EXACTLY the corresponding primitive with the documented (possibly swapped) args — companion to iter-549's bivariate
dispatch test; a mis-wired branch or wrong-swap bug is caught directly. Oracle-independent; never SKIPs. (iter 550)
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
import jax.numpy as jnp

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import trivar_NLL_mean_eq, CHI_TOL
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    trivar_NLL_mean, trivar_NLL_mean_const_x1, trivar_NLL_mean_const_x2,
    trivar_NLL_mean_const_x1x2, trivar_NLL_mean_const_x2x3, trivar_NLL_mean_const_all)

_A, _B, _G = 1.0, 0.9, 0.7      # alpha, beta, gamma — all ≥ 0 (so mu_x2=mu_rr, mu_x3=mu_Nr)
# nominal moments (chi mean, in-precip lognormal means/stdevs, correlations)
_MU_CHI, _MU_RR, _MU_NR, _MU_RR_N, _MU_NR_N = 0.4, 1.0e-4, 1.0e6, -9.2, 13.8
_S_RR_N, _S_NR_N, _R12, _R13, _R23 = 0.5, 0.6, 0.2, 0.15, 0.1


def _eq(s_chi, s_rr, s_nr, mu_chi=_MU_CHI):
    return float(trivar_NLL_mean_eq(mu_chi, _MU_RR, _MU_NR, _MU_RR_N, _MU_NR_N,
                                    s_chi, s_rr, s_nr, _S_RR_N, _S_NR_N,
                                    _R12, _R13, _R23, _A, _B, _G))


def _close(got, ref, name):
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"{name}: {got} vs {ref}"


def test_all_vary():
    _close(_eq(0.3, 2e-4, 2e5),
           float(trivar_NLL_mean(_MU_CHI, _MU_RR_N, _MU_NR_N, max(0.3, CHI_TOL), _S_RR_N, _S_NR_N,
                                 _R12, _R13, _R23, _A, _B, _G)), "all-vary->trivar_NLL_mean")
    print("  all vary -> trivar_NLL_mean (general)  PASS")


def test_all_const():
    _close(_eq(0.0, 0.0, 0.0),
           float(trivar_NLL_mean_const_all(_MU_CHI, _MU_RR, _MU_NR, _A, _B, _G)), "all-const->const_all")
    print("  all σ≤tol -> const_all  PASS")


def test_x1x2_const():
    _close(_eq(0.0, 0.0, 2e5),
           float(trivar_NLL_mean_const_x1x2(_MU_CHI, _MU_RR, _MU_NR_N, _S_NR_N, _A, _B, _G)), "x1x2->const_x1x2")
    print("  chi&rr const -> const_x1x2(mu_rr, mu_Nr_n, sig_Nr_n)  PASS")


def test_x1x3_const_uses_swapped_x1x2():
    # chi & Nr const → reuse const_x1x2 with (x2,x3)=(Nr,rr) swapped and (beta,gamma) swapped.
    _close(_eq(0.0, 2e-4, 0.0),
           float(trivar_NLL_mean_const_x1x2(_MU_CHI, _MU_NR, _MU_RR_N, _S_RR_N, _A, _G, _B)),
           "x1x3->const_x1x2 swapped")
    print("  chi&Nr const -> const_x1x2 with (x2,x3)/(beta,gamma) SWAPPED  PASS")


def test_x2x3_const():
    _close(_eq(0.3, 0.0, 0.0),
           float(trivar_NLL_mean_const_x2x3(_MU_CHI, _MU_RR, _MU_NR, max(0.3, CHI_TOL), _A, _B, _G)), "x2x3->const_x2x3")
    print("  rr&Nr const -> const_x2x3  PASS")


def test_x1_const():
    _close(_eq(0.0, 2e-4, 2e5),
           float(trivar_NLL_mean_const_x1(_MU_CHI, _MU_RR_N, _MU_NR_N, _S_RR_N, _S_NR_N, _R23, _A, _B, _G)),
           "x1->const_x1")
    print("  chi const -> const_x1(rho23)  PASS")


def test_x2_const():
    _close(_eq(0.3, 0.0, 2e5),
           float(trivar_NLL_mean_const_x2(_MU_CHI, _MU_RR, _MU_NR_N, max(0.3, CHI_TOL), _S_NR_N, _R13, _A, _B, _G)),
           "x2->const_x2")
    print("  rr const -> const_x2(mu_rr, mu_Nr_n, rho13)  PASS")


def test_x3_const_uses_swapped_x2():
    # Nr const → reuse const_x2 with (x2,x3)=(Nr,rr) swapped, (beta,gamma) swapped, rho12.
    _close(_eq(0.3, 2e-4, 0.0),
           float(trivar_NLL_mean_const_x2(_MU_CHI, _MU_NR, _MU_RR_N, max(0.3, CHI_TOL), _S_RR_N, _R12, _A, _G, _B)),
           "x3->const_x2 swapped")
    print("  Nr const -> const_x2 with (x2,x3)/(beta,gamma) SWAPPED, rho12  PASS")


def main():
    print("test_trivar_NLL_mean_eq:")
    for t in (test_all_vary, test_all_const, test_x1x2_const, test_x1x3_const_uses_swapped_x1x2,
              test_x2x3_const, test_x1_const, test_x2_const, test_x3_const_uses_swapped_x2):
        t()
    print("All trivar_NLL_mean_eq dispatch checks PASSED")


if __name__ == "__main__":
    main()

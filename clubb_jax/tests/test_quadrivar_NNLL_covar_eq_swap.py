#!/usr/bin/env python3
"""test_quadrivar_NNLL_covar_eq_swap.py — pin the SYMMETRY-SWAP logic of the KK quadrivariate covariance dispatch.

`quadrivar_NNLL_covar_eq` (KK_upscaled_covariances.py) dispatches Cov_i(x, chi^α r_r^β N_r^γ) over the full 15-way
variance-regime ladder. Its UNIQUE, most bug-prone feature (vs the trivariate dispatch pinned iter 553) is the
(x3=r_r, β) ↔ (x4=N_r, γ) SYMMETRY: an N_r-const (c4) branch reuses the corresponding r_r-const (x3) primitive with
the x3↔x4 args and β↔γ exponents SWAPPED (v_x3s, v123s, …). The covar primitives are end-to-end-covered
(test_kk_rico_oracle) but this swap wiring was never pinned — and rico may not exercise the N_r-const regimes. This
focuses on the swap branches (+ a non-swap contrast and the base/all-const sanity) and asserts each returns EXACTLY
the documented primitive with the (possibly swapped) args. A wrong-swap bug — easy to introduce, hard to catch
end-to-end — is now caught directly. Oracle-independent; never SKIPs. (iter 554)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_covariances import quadrivar_NNLL_covar_eq, _CHI_TOL
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_covar import (
    quadrivar_NNLL_covar, quadrivar_NNLL_covar_const_x3,
    quadrivar_NNLL_covar_cst_x1x2x3, quadrivar_NNLL_covar_const_all)

# nominal moments: x=w/eta, chi, r_r (lognormal), N_r (lognormal)
_MUX, _MUCHI, _MURR, _MUNR, _MURRN, _MUNRN = 0.5, 0.3, 1.0e-4, 1.0e6, -9.2, 13.8
_SRRN, _SNRN = 0.5, 0.6
_R12, _R13N, _R14N, _R23N, _R24N, _R34N = 0.2, 0.15, 0.1, 0.12, 0.08, 0.05
_XM, _MCT, _MCC = 0.4, 2.0, 1.5
_XTOL, _RRTOL, _NRTOL = 1.0e-8, 1.0e-10, 1.909859e-7
_A, _B, _G = 1.0, 0.9, 0.7
_X2A = _MCT / _MCC


def _eq(s_x, s_chi, s_rr, s_nr):
    return float(quadrivar_NNLL_covar_eq(_MUX, _MUCHI, _MURR, _MUNR, _MURRN, _MUNRN,
                                         s_x, s_chi, s_rr, s_nr, _SRRN, _SNRN,
                                         _R12, _R13N, _R14N, _R23N, _R24N, _R34N,
                                         _XM, _MCT, _MCC, _XTOL, _RRTOL, _NRTOL, _A, _B, _G))


def _s2g(s_chi):
    return s_chi if s_chi > _CHI_TOL else 1.0


def _close(got, ref, name):
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"{name}: {got} vs {ref}"


def test_base_all_vary():
    _close(_eq(0.4, 0.3, 2e-4, 2e5),
           float(quadrivar_NNLL_covar(_MUX, _MUCHI, _MURRN, _MUNRN, 0.4, _s2g(0.3), _SRRN, _SNRN,
                                      _R12, _R13N, _R14N, _R23N, _R24N, _R34N, _XM, _X2A, _A, _B, _G)), "base")
    print("  all vary -> quadrivar_NNLL_covar (base)  PASS")


def test_c3_const_x3_NONswap():
    # r_r const (c3): non-swapped x3 variant — uses mu_x3, mu_x4_n, sg4n, r14n, r24n and (a,b,g).
    _close(_eq(0.4, 0.3, 0.0, 2e5),
           float(quadrivar_NNLL_covar_const_x3(_MUX, _MUCHI, _MURR, _MUNRN, 0.4, _s2g(0.3), _SNRN,
                                               _R12, _R14N, _R24N, _XM, _X2A, _A, _B, _G)), "const_x3")
    print("  r_r const -> const_x3 (NON-swapped: mu_x4_n, sg_Nr_n, r14n, r24n, b,g)  PASS")


def test_c4_const_x3_SWAPPED():
    # N_r const (c4): reuse the x3 variant with x3↔x4 args + β↔γ swapped — mu_x4, mu_x3_n, sg3n, r13n, r23n, (a,g,b).
    _close(_eq(0.4, 0.3, 2e-4, 0.0),
           float(quadrivar_NNLL_covar_const_x3(_MUX, _MUCHI, _MUNR, _MURRN, 0.4, _s2g(0.3), _SRRN,
                                               _R12, _R13N, _R23N, _XM, _X2A, _A, _G, _B)), "const_x3 swapped")
    print("  N_r const -> const_x3 with x3↔x4 args + β↔γ SWAPPED  PASS")


def test_c1c2c4_cst_x1x2x3_SWAPPED():
    # x & chi & N_r const: reuse cst_x1x2x3 with x3↔x4 (mu_x4, mu_x3_n, sg3n) and (a,g,b).
    _close(_eq(0.0, 0.0, 2e-4, 0.0),
           float(quadrivar_NNLL_covar_cst_x1x2x3(_MUX, _MUCHI, _MUNR, _MURRN, _SRRN, _XM, _X2A, _A, _G, _B)),
           "cst_x1x2x3 swapped")
    print("  x&chi&N_r const -> cst_x1x2x3 with x3↔x4 + β↔γ SWAPPED  PASS")


def test_all_const():
    _close(_eq(0.0, 0.0, 0.0, 0.0),
           float(quadrivar_NNLL_covar_const_all(_MUX, _MUCHI, _MURR, _MUNR, _XM, _X2A, _A, _B, _G)), "const_all")
    print("  all σ≤tol -> const_all  PASS")


def main():
    print("test_quadrivar_NNLL_covar_eq_swap:")
    for t in (test_base_all_vary, test_c3_const_x3_NONswap, test_c4_const_x3_SWAPPED,
              test_c1c2c4_cst_x1x2x3_SWAPPED, test_all_const):
        t()
    print("All quadrivar_NNLL_covar_eq swap-dispatch checks PASSED")


if __name__ == "__main__":
    main()

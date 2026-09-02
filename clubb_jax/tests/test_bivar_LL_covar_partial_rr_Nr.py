#!/usr/bin/env python3
"""test_bivar_LL_covar_partial_rr_Nr.py — pin the r_r/N_r covariance-partial wrappers + the _covar_partial dispatch.

`bivar_LL_covar_partial_rr` / `bivar_LL_covar_partial_Nr` (KK_upscaled_turbulent_sed.py) are the <r_r' R_vr'> /
<N_r' R_vr'> coefA partials. Each assigns the KK mean-volume-radius exponents (rr_exp=+1/3, Nr_exp=−1/3) and routes
r_r/N_r to (x1,x2) — with r_r↔N_r SWAPPED for the N_r wrapper — then calls the private 4-way dispatch `_covar_partial`.
That dispatch selects among FOUR now-tested primitives per the (near-)zero-variance regime, with one documented
subtlety: the x1-const branch reuses the MEAN's `bivar_LL_mean_const_x1` (the extra x1' factor vanishes as σ_x1→0),
NOT a covariance partial. The primitives are tested (iters 551/555) but the wrapper routing + this dispatch were not.
This drives all 4 regimes of `_partial_rr` (asserting each == the right tested primitive with the right args/exponents)
plus `_partial_Nr` all-vary (the r_r↔N_r swap). A mis-routed var / wrong exponent / wrong-primitive-per-regime bug is
caught directly. Oracle-independent; never SKIPs. (iter 556)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import (
    bivar_LL_covar_partial_rr, bivar_LL_covar_partial_Nr, KK_MVR_RR_EXP, KK_MVR_NR_EXP)
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    bivar_LL_mean_const_x1, bivar_LL_mean_const_all)
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import (
    bivar_LL_covar_partial, bivar_LL_covar_const_x2_partial)

_MURR, _MUNR, _MURRN, _MUNRN, _SRRN, _SNRN, _RHO = 1.0e-4, 1.0e6, -9.2, 13.8, 0.5, 0.6, 0.1
_A, _B = KK_MVR_RR_EXP, KK_MVR_NR_EXP   # +1/3, −1/3


def _rr(s_rr, s_Nr):
    return float(bivar_LL_covar_partial_rr(_MURR, _MUNR, _MURRN, _MUNRN, s_rr, s_Nr, _SRRN, _SNRN, _RHO))


def _close(got, ref, name):
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"{name}: {got} vs {ref}"


def test_rr_all_vary_general():
    _close(_rr(2e-4, 2e5),
           float(bivar_LL_covar_partial(_MURRN, _MUNRN, _SRRN, _SNRN, _RHO, _A, _B)), "rr general")
    print("  partial_rr all-vary -> bivar_LL_covar_partial(rr_n, Nr_n, α=+1/3, β=−1/3)  PASS")


def test_rr_x1_const_reuses_mean_const_x1():
    # σ_rr→0 (x1 const): the documented MEAN reuse — bivar_LL_mean_const_x1(mu_rr, mu_Nr_n, σ_Nr_n, α, β).
    _close(_rr(0.0, 2e5),
           float(bivar_LL_mean_const_x1(_MURR, _MUNRN, _SNRN, _A, _B)), "rr x1-const")
    print("  partial_rr σ_rr=0 -> MEAN's bivar_LL_mean_const_x1 (extra x1' factor vanishes)  PASS")


def test_rr_x2_const_partial():
    _close(_rr(2e-4, 0.0),
           float(bivar_LL_covar_const_x2_partial(_MURRN, _MUNR, _SRRN, _A, _B)), "rr x2-const")
    print("  partial_rr σ_Nr=0 -> bivar_LL_covar_const_x2_partial  PASS")


def test_rr_both_const_all():
    _close(_rr(0.0, 0.0),
           float(bivar_LL_mean_const_all(_MURR, _MUNR, _A, _B)), "rr both-const")
    print("  partial_rr both σ=0 -> bivar_LL_mean_const_all  PASS")


def test_Nr_all_vary_swapped():
    # N_r wrapper: x1=N_r, x2=r_r swapped; α=Nr_exp, β=rr_exp.
    got = float(bivar_LL_covar_partial_Nr(_MURR, _MUNR, _MURRN, _MUNRN, 2e5, 2e-4, _SRRN, _SNRN, _RHO))
    _close(got, float(bivar_LL_covar_partial(_MUNRN, _MURRN, _SNRN, _SRRN, _RHO, _B, _A)), "Nr general swapped")
    print("  partial_Nr all-vary -> bivar_LL_covar_partial(Nr_n, rr_n, α=−1/3, β=+1/3) (r_r↔N_r SWAPPED)  PASS")


def main():
    print("test_bivar_LL_covar_partial_rr_Nr:")
    for t in (test_rr_all_vary_general, test_rr_x1_const_reuses_mean_const_x1, test_rr_x2_const_partial,
              test_rr_both_const_all, test_Nr_all_vary_swapped):
        t()
    print("All bivar_LL_covar_partial_rr/Nr checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_bivar_LL_covar_partial.py — validate the KK lognormal-lognormal covariance partials (turbulent sed).

`bivar_LL_covar_partial` and `bivar_LL_covar_const_x2_partial` (KK_upscaled_turbulent_sed.py ↔ .F90:1540/1584) are
the coefA building blocks <x1' x1^α x2^β> for two correlated LOGNORMALs, used by the KK turbulent-sedimentation
covariances. They were untested. Each is a clean closed form, and the general one is exactly the TESTED `bivar_LL_mean`
moment with one extra x1 factor, normalized by <x1>:

    bivar_LL_covar_partial(μ1n,μ2n,σ1,σ2,ρ,α,β) == bivar_LL_mean(μ1n,μ2n,σ1,σ2,ρ, α+1, β) / exp(μ1n + ½σ1²)
      (the (α²+2α) = (α+1)²−1 and (α+1) cross-term in the partial are precisely this extra-x1-factor / <x1> shift)

This validates them three ways (real math, not wiring): (1) their closed forms; (2) the general partial vs the
tested `bivar_LL_mean` via that identity (an independent cross-check); (3) the x2-const partial as the σ2→0 limit of
the general partial. + finite grad. Oracle-independent; never SKIPs. (iter 555)
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

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import (
    bivar_LL_covar_partial, bivar_LL_covar_const_x2_partial)
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import bivar_LL_mean

_CASES = [
    # (mu1_n, mu2_n, sigma1_n, sigma2_n, rho, alpha, beta)
    (-9.2, 13.8, 0.5, 0.6, 0.2, 1.0, 0.9),
    (-7.0, 12.0, 0.4, 0.45, -0.3, 1.15, 0.5),
    (-10.5, 15.0, 0.7, 0.3, 0.0, 0.83, 1.0),
]


def test_general_partial_closed_form_and_mean_identity():
    worst_cf = 0.0; worst_id = 0.0
    for m1, m2, s1, s2, rho, a, b in _CASES:
        got = float(bivar_LL_covar_partial(m1, m2, s1, s2, rho, a, b))
        cf = math.exp(m1 * a + m2 * b + 0.5 * s1 ** 2 * (a ** 2 + 2.0 * a)
                      + 0.5 * s2 ** 2 * b ** 2 + rho * s1 * (a + 1.0) * s2 * b)
        worst_cf = max(worst_cf, abs(got - cf) / (abs(cf) + 1e-300))
        # identity vs the TESTED bivar_LL_mean: == mean(α+1,β) / <x1>,  <x1> = exp(μ1n + ½σ1²)
        mean_ap1 = float(bivar_LL_mean(m1, m2, s1, s2, rho, a + 1.0, b))
        ident = mean_ap1 / math.exp(m1 + 0.5 * s1 ** 2)
        worst_id = max(worst_id, abs(got - ident) / (abs(ident) + 1e-300))
    assert worst_cf < 1e-12 and worst_id < 1e-12, f"partial cf {worst_cf:.2e}, mean-identity {worst_id:.2e}"
    print(f"  bivar_LL_covar_partial: closed form (cf {worst_cf:.1e}) == bivar_LL_mean(α+1,β)/<x1> ({worst_id:.1e})  PASS")


def test_const_x2_partial_closed_form_and_limit():
    worst_cf = 0.0; worst_lim = 0.0
    for m1, m2, s1, s2, rho, a, b in _CASES:
        mu2 = math.exp(m2 + 0.5 * s2 ** 2)        # an arbitrary positive constant x2 value
        got = float(bivar_LL_covar_const_x2_partial(m1, mu2, s1, a, b))
        cf = mu2 ** b * math.exp(m1 * a + 0.5 * s1 ** 2 * (a ** 2 + 2.0 * a))
        worst_cf = max(worst_cf, abs(got - cf) / (abs(cf) + 1e-300))
        # σ2→0 limit of the general partial with μ2n = ln(mu2): x2 deterministic at mu2.
        lim = float(bivar_LL_covar_partial(m1, math.log(mu2), s1, 0.0, rho, a, b))
        worst_lim = max(worst_lim, abs(got - lim) / (abs(lim) + 1e-300))
    assert worst_cf < 1e-12 and worst_lim < 1e-12, f"const_x2 cf {worst_cf:.2e}, limit {worst_lim:.2e}"
    print(f"  bivar_LL_covar_const_x2_partial: closed form (cf {worst_cf:.1e}) == general at σ2=0 ({worst_lim:.1e})  PASS")


def test_grad_finite():
    def loss(s1):
        return bivar_LL_covar_partial(-9.2, 13.8, s1, 0.6, 0.2, 1.0, 0.9) ** 2 \
            + bivar_LL_covar_const_x2_partial(-9.2, 1.0e-4, s1, 1.0, 0.9) ** 2
    g = jax.grad(loss)(0.5)
    assert np.isfinite(float(g)), "non-finite grad"
    print("  jax.grad of both covar partials wrt sigma_x1 finite  PASS")


def main():
    print("test_bivar_LL_covar_partial:")
    test_general_partial_closed_form_and_mean_identity()
    test_const_x2_partial_closed_form_and_limit()
    test_grad_finite()
    print("All bivar_LL_covar_partial checks PASSED")


if __name__ == "__main__":
    main()

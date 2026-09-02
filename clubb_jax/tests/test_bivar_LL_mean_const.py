#!/usr/bin/env python3
"""test_bivar_LL_mean_const.py — validate the constant-variance limits of the KK lognormal-lognormal mean integral.

`bivar_LL_mean_const_x1` and `bivar_LL_mean_const_all` (PDF_integrals_means.py ↔ PDF_integrals_means.F90:831/880) are
the σ→0 limits of `bivar_LL_mean` (the general <x1^α x2^β> for two correlated LOGNORMAL variables, tested in
test_kk_autoconversion.py). They were the only untested PDF_integrals_means LL primitives. Each has a clean closed
form (lognormal moments), and — crucially — each must coincide with the general `bivar_LL_mean` evaluated in the
corresponding zero-variance limit (with the constant variable's log-mean = ln(value)). This validates them two ways:
  1. closed form: const_all = μ_x1^α·μ_x2^β ;  const_x1 = μ_x1^α·exp(μ_x2n·β + ½σ_x2n²·β²)
  2. CONSISTENCY with the tested general form: const_x1 == bivar_LL_mean(ln μ_x1, μ_x2n, 0, σ_x2n, ρ, α, β),
     const_all == bivar_LL_mean(ln μ_x1, ln μ_x2, 0, 0, ρ, α, β)   — a genuine cross-check, not a tautology.
Oracle-independent; never SKIPs. (iter 551)
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

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    bivar_LL_mean, bivar_LL_mean_const_x1, bivar_LL_mean_const_all)

# representative KK exponents + lognormal params (x = r_r / N_r, positive hydrometeor means)
_CASES = [
    # (mu_x1, mu_x2, mu_x2_n, sigma_x2_n, rho, alpha, beta)
    (1.0e-4, 1.0e6, 13.8, 0.6, 0.2, 1.0, 0.9),
    (5.0e-4, 2.0e5, 12.2, 0.45, -0.3, 1.15, 0.5),
    (2.0e-3, 8.0e6, 15.9, 0.7, 0.0, 0.83, 1.0),
]


def test_const_all_closed_form_and_limit():
    worst_cf = 0.0; worst_lim = 0.0
    for mu1, mu2, mu2n, s2n, rho, a, b in _CASES:
        got = float(bivar_LL_mean_const_all(mu1, mu2, a, b))
        cf = mu1 ** a * mu2 ** b
        worst_cf = max(worst_cf, abs(got - cf) / (abs(cf) + 1e-300))
        # limit of the general LL mean: both σ→0, log-means = ln(value)
        lim = float(bivar_LL_mean(math.log(mu1), math.log(mu2), 0.0, 0.0, rho, a, b))
        worst_lim = max(worst_lim, abs(got - lim) / (abs(lim) + 1e-300))
    assert worst_cf < 1e-13 and worst_lim < 1e-12, f"const_all cf {worst_cf:.2e}, limit {worst_lim:.2e}"
    print(f"  const_all = μ1^α·μ2^β (cf {worst_cf:.1e}) == general LL at σ=0 (limit {worst_lim:.1e})  PASS")


def test_const_x1_closed_form_and_limit():
    worst_cf = 0.0; worst_lim = 0.0
    for mu1, mu2, mu2n, s2n, rho, a, b in _CASES:
        got = float(bivar_LL_mean_const_x1(mu1, mu2n, s2n, a, b))
        cf = mu1 ** a * math.exp(mu2n * b + 0.5 * s2n ** 2 * b ** 2)
        worst_cf = max(worst_cf, abs(got - cf) / (abs(cf) + 1e-300))
        # general LL with σ_x1=0 and μ_x1_n = ln(μ_x1) — the ρ term drops out (×σ_x1=0)
        lim = float(bivar_LL_mean(math.log(mu1), mu2n, 0.0, s2n, rho, a, b))
        worst_lim = max(worst_lim, abs(got - lim) / (abs(lim) + 1e-300))
    assert worst_cf < 1e-13 and worst_lim < 1e-12, f"const_x1 cf {worst_cf:.2e}, limit {worst_lim:.2e}"
    print(f"  const_x1 = μ1^α·exp(μ2n·β+½σ2n²β²) (cf {worst_cf:.1e}) == general LL at σ_x1=0 (limit {worst_lim:.1e})  PASS")


def test_grad_finite():
    def loss(mu1):
        return bivar_LL_mean_const_x1(mu1, 13.8, 0.6, 1.0, 0.9) ** 2 \
            + bivar_LL_mean_const_all(mu1, 1.0e6, 1.0, 0.9) ** 2
    g = jax.grad(loss)(1.0e-4)
    assert np.isfinite(float(g)), "non-finite grad"
    print("  jax.grad of both const forms wrt μ_x1 finite  PASS")


def main():
    print("test_bivar_LL_mean_const:")
    test_const_all_closed_form_and_limit()
    test_const_x1_closed_form_and_limit()
    test_grad_finite()
    print("All bivar_LL_mean const-limit checks PASSED")


if __name__ == "__main__":
    main()

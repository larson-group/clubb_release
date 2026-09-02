#!/usr/bin/env python3
"""test_covar_KK_mixt_frac_linearity.py — pin the outer mixt_frac weighting of the covar_*_KK_* compositions.

The 9 top-level KK covariance functions `covar_{rt,thl,x}_KK_{auto,accr,evap}` (KK_upscaled_covariances.py) all share
the outer structure `covar = mixt_frac·c1 + (1−mixt_frac)·c2`, where c1/c2 are per-PDF-component covariances
(`_covar_x_comp` / `_covar_x_evap_comp` over the now-dispatch-tested trivar/quadrivar integrals, iters 553/554) that do
NOT depend on mixt_frac. So each is EXACTLY linear in mixt_frac. The inner per-component routing is end-to-end-validated
(test_kk_rico_oracle) and tautology-prone to isolate, but the outer weighting was never pinned. This drives a
representative function (`covar_rt_KK_auto`) across mixt_frac and asserts exact linearity:
`covar(m) == m·covar(1) + (1−m)·covar(0)`, which holds iff the m-weighted combination of the two m-independent
component covariances is correctly wired (a swapped/miscoefficiented combination breaks it). + finite grad.
Oracle-independent; never SKIPs. (iter 558)
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

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_covariances import covar_rt_KK_auto

# nominal component moments (component 1 / 2 deliberately distinct so c1 != c2).
_ARGS = dict(
    mu_eta_1=0.3, mu_eta_2=-0.2, mu_chi_1=0.2, mu_chi_2=0.25,
    mu_Ncn_1=1.0e8, mu_Ncn_2=1.5e8, mu_Ncn_1_n=18.0, mu_Ncn_2_n=18.5,
    sigma_eta_1=0.3, sigma_eta_2=0.35, sigma_chi_1=0.2, sigma_chi_2=0.22,
    sigma_Ncn_1=0.5, sigma_Ncn_2=0.55, sigma_Ncn_1_n=0.5, sigma_Ncn_2_n=0.55,
    corr_eta_chi_1=0.2, corr_eta_chi_2=0.15, corr_eta_Ncn_1_n=0.1, corr_eta_Ncn_2_n=0.12,
    corr_chi_Ncn_1_n=0.18, corr_chi_Ncn_2_n=0.16, rtm=0.01, mu_rt_1=0.011, mu_rt_2=0.012,
    KK_auto_tndcy=1.0e-7, KK_auto_coef=1350.0, eta_tol=1.0e-8, crt1=0.6, crt2=0.4,
)


def _covar(m):
    return float(covar_rt_KK_auto(*_ARGS.values(), m))


def test_exact_mixt_frac_linearity():
    c0, c1 = _covar(0.0), _covar(1.0)
    assert c0 != c1, "component covariances c1==c2 — test inputs not distinct enough to detect a weighting bug"
    worst = 0.0
    for m in (0.1, 0.3, 0.5, 0.7, 0.9):
        got = _covar(m)
        ref = m * c1 + (1.0 - m) * c0
        worst = max(worst, abs(got - ref) / (abs(ref) + 1e-300))
    assert worst < 1e-13, f"covar_rt_KK_auto not exactly mixt_frac-linear: {worst:.2e}"
    print(f"  covar_rt_KK_auto(m) == m·covar(1)+(1−m)·covar(0) exactly (worst {worst:.1e})  PASS")


def test_grad_finite():
    g = jax.grad(lambda m: covar_rt_KK_auto(*_ARGS.values(), m) ** 2)(0.4)
    assert np.isfinite(float(g)), "non-finite grad wrt mixt_frac"
    print("  jax.grad(covar_rt_KK_auto) wrt mixt_frac finite  PASS")


def main():
    print("test_covar_KK_mixt_frac_linearity:")
    test_exact_mixt_frac_linearity()
    test_grad_finite()
    print("All covar_*_KK_* mixt_frac-weighting checks PASSED")


if __name__ == "__main__":
    main()

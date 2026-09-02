#!/usr/bin/env python3
"""test_transform_pdf_chi_eta_component.py — validate the Sommeria-Deardorff χ/η PDF-component transform.

`transform_pdf_chi_eta_component(tl, rsatl, rt, exner, varnce_rt, varnce_thl, corr_rt_thl)` (pdf_closure_module.py ↔
pdf_closure_module.F90:1699) is the extended-liquid-water transform of one PDF component from (rt, thl) to (χ, η):
    β = ep·Lv²/(Rd·Cp·tl²);  invrs = 1/(1+β·rsatl);  χ = (rt−rsatl)·invrs
    crt = invrs;  cthl = (1+β·rt)·invrs²·(Cp/Lv)·β·rsatl·exner
and the variances of the LINEAR combinations χ' = crt·rt' − cthl·thl', η' = crt·rt' + cthl·thl':
    σ_χ² = crt²σ_rt² − 2·corr·crt·cthl·σ_rt·σ_thl + cthl²σ_thl²,   σ_η² = (… + 2·corr …),
    covar_χη = crt²σ_rt² − cthl²σ_thl².
It feeds the cloud-fraction closure but was validated only end-to-end. This pins (1) the deterministic coefficients
(β, χ, crt, cthl) vs the closed form, and — crucially, INDEPENDENT — (2) the variance/covariance combination vs a
Monte-Carlo: draw (rt', thl') from the bivariate normal, form χ'/η' by the same linear maps, and compare the empirical
Var(χ'), Var(η'), Cov(χ',η') to σ_χ², σ_η², covar_χη. + the corr_chi_eta quotient and a finite grad. (iter 565)
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import transform_pdf_chi_eta_component
from clubb_jax.src.CLUBB_core.constants_clubb import ep, Lv, Rd, Cp

_TL, _RSATL, _RT, _EXNER = 290.0, 0.012, 0.011, 0.95
_VRT, _VTHL, _CORR = 1.0e-6, 1.0, 0.3


def _ref_coefs():
    beta = ep * Lv ** 2 / (Rd * Cp * _TL ** 2)
    invrs = 1.0 / (1.0 + beta * _RSATL)
    chi = (_RT - _RSATL) * invrs
    crt = invrs
    cthl = (1.0 + beta * _RT) * invrs ** 2 * (Cp / Lv) * beta * _RSATL * _EXNER
    return beta, chi, crt, cthl


def test_deterministic_coefficients():
    chi, crt, cthl, sdchi, sdeta, covar, corr = (float(x) for x in transform_pdf_chi_eta_component(
        _TL, _RSATL, _RT, _EXNER, _VRT, _VTHL, _CORR))
    beta, r_chi, r_crt, r_cthl = _ref_coefs()
    assert abs(chi - r_chi) < 1e-15 and abs(crt - r_crt) < 1e-15 and abs(cthl - r_cthl) < 1e-18, "coef mismatch"
    print(f"  β/χ/crt/cthl match closed form (χ={chi:.3e}, crt={crt:.4f}, cthl={cthl:.3e})  PASS")


def test_variance_combination_monte_carlo():
    _, _, crt, cthl = _ref_coefs()
    _, _, _, sdchi, sdeta, covar, _ = (float(x) for x in transform_pdf_chi_eta_component(
        _TL, _RSATL, _RT, _EXNER, _VRT, _VTHL, _CORR))
    rng = np.random.default_rng(565)
    N = 12_000_000
    srt, sthl = math.sqrt(_VRT), math.sqrt(_VTHL)
    z1, z2 = rng.standard_normal(N), rng.standard_normal(N)
    rtp = srt * z1
    thlp = sthl * (_CORR * z1 + math.sqrt(1 - _CORR ** 2) * z2)
    chip = crt * rtp - cthl * thlp     # χ' = crt·rt' − cthl·thl'
    etap = crt * rtp + cthl * thlp     # η' = crt·rt' + cthl·thl'
    mc_vchi, mc_veta, mc_cov = np.var(chip), np.var(etap), np.mean(chip * etap)
    r1 = abs(sdchi ** 2 - mc_vchi) / mc_vchi
    r2 = abs(sdeta ** 2 - mc_veta) / mc_veta
    r3 = abs(covar - mc_cov) / (abs(mc_cov) + 1e-300)
    assert r1 < 5e-3 and r2 < 5e-3 and r3 < 5e-3, f"variance MC: σχ² {r1:.1e}, ση² {r2:.1e}, cov {r3:.1e}"
    print(f"  σ_χ²/σ_η²/covar_χη vs {N//1_000_000}M-sample linear-transform MC: rel <5e-3 ({max(r1,r2,r3):.1e})  PASS")


def test_corr_quotient_and_grad():
    _, _, _, sdchi, sdeta, covar, corr = (float(x) for x in transform_pdf_chi_eta_component(
        _TL, _RSATL, _RT, _EXNER, _VRT, _VTHL, _CORR))
    # non-degenerate: corr_chi_eta ≈ covar / (stdev_chi·stdev_eta), within [-1,1]
    assert abs(corr - covar / (sdchi * sdeta)) < 1e-9 and -1.0 <= corr <= 1.0, "corr_chi_eta quotient off"
    g = jax.grad(lambda rt: transform_pdf_chi_eta_component(_TL, _RSATL, rt, _EXNER, _VRT, _VTHL, _CORR)[0])(_RT)
    assert np.isfinite(float(g)), "non-finite grad of chi wrt rt"
    print("  corr_chi_eta = covar/(σ_χ·σ_η) ∈ [−1,1]; finite grad  PASS")


def main():
    print("test_transform_pdf_chi_eta_component:")
    test_deterministic_coefficients()
    test_variance_combination_monte_carlo()
    test_corr_quotient_and_grad()
    print("All transform_pdf_chi_eta_component checks PASSED")


if __name__ == "__main__":
    main()

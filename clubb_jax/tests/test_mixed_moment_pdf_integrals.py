#!/usr/bin/env python3
"""test_mixed_moment_pdf_integrals.py — validate the closed-form PDF-component moment integrals.

No f2py oracle exists for these, so each closed form is checked against TWO independent references:
  1. A from-scratch raw-moment binomial expansion  E[(X-c)^n] = SUM_k C(n,k) (-c)^(n-k) E[X^k], with the raw
     moments computed independently (normal: central moments via double factorial; lognormal: exp formula).
  2. A Monte-Carlo sample of the component PDF (normal / lognormal), sample mean of (x-c)^n.
Plus a finite jax.grad.
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
import jax.numpy as jnp

from clubb_jax.src.Microphys.mixed_moment_PDF_integrals import (
    univar_N_int_PDF_comp_all_MM, univar_L_int_PDF_comp_all_MM, bivar_NL_int_PDF_comp_all_MM,
    bivar_NL_x_hm_all_MM_comp_eq, xp_a_hmpb_integrals_all_MM,
    xphmp_integral_covar, hmxphmyp_integral_covar)


def _normal_raw_central(mu, sigma, c, n):
    """Independent E[(X-c)^n], X~N(mu,sigma): expand around mu using central moments of a normal."""
    def central_normal(k):   # E[(X-mu)^k]
        if k % 2 == 1:
            return 0.0
        return sigma ** k * math.factorial(k) / (2 ** (k // 2) * math.factorial(k // 2))  # = sigma^k (k-1)!!
    d = mu - c
    return sum(math.comb(n, k) * d ** (n - k) * central_normal(k) for k in range(n + 1))


def _lognormal_raw_central(mu_n, sigma_n, c, n):
    """Independent E[(X-c)^n], ln X ~ N(mu_n,sigma_n): expand using lognormal raw moments E[X^k]."""
    raw = lambda k: math.exp(mu_n * k + 0.5 * sigma_n ** 2 * k ** 2)
    return sum(math.comb(n, k) * (-c) ** (n - k) * raw(k) for k in range(n + 1))


def test_univar_N_closed_form_and_mc():
    mu, sigma = 1.3, 0.7
    rng = np.random.default_rng(11)
    x = rng.normal(mu, sigma, 8_000_000)
    c = 0.9   # overall mean <x>, distinct from component mean
    for a in (0, 1, 2, 3, 4, 5):
        got = float(univar_N_int_PDF_comp_all_MM(mu, sigma, c, a))
        ref = _normal_raw_central(mu, sigma, c, a)
        rel = abs(got - ref) / (abs(ref) + 1e-30)
        assert rel < 1e-12, f"univar_N a={a}: closed-form rel {rel:.2e}"
        mc = np.mean((x - c) ** a)
        # Normalise by the ABSOLUTE moment so the check stays well-conditioned when the signed
        # central moment is ~0 (e.g. odd orders with mu near c).
        scale = np.mean(np.abs(x - c) ** a) + 1e-30
        rel_mc = abs(got - mc) / scale
        assert rel_mc < 5e-3, f"univar_N a={a}: MC rel {rel_mc:.2e}"
    print("  univar_N (a=0..5) vs raw-moment expansion (<1e-12) + 8M-sample MC (<5e-3)  PASS")


def test_univar_L_closed_form_and_mc():
    mu_n, sigma_n = -0.4, 0.5
    rng = np.random.default_rng(22)
    x = np.exp(rng.normal(mu_n, sigma_n, 8_000_000))
    c = float(np.exp(mu_n + 0.5 * sigma_n ** 2))   # overall lognormal mean
    for b in (0, 1, 2, 3, 4):
        got = float(univar_L_int_PDF_comp_all_MM(mu_n, sigma_n, c, b))
        ref = _lognormal_raw_central(mu_n, sigma_n, c, b)
        rel = abs(got - ref) / (abs(ref) + 1e-30)
        assert rel < 1e-12, f"univar_L b={b}: closed-form rel {rel:.2e}"
        mc = np.mean((x - c) ** b)
        scale = np.mean(np.abs(x - c) ** b) + 1e-30   # absolute moment → well-conditioned at signed ~0
        rel_mc = abs(got - mc) / scale
        assert rel_mc < 1e-2, f"univar_L b={b}: MC rel {rel_mc:.2e}"
    print("  univar_L (b=0..4) vs raw-moment expansion (<1e-12) + 8M-sample MC (<1e-2)  PASS")


def _bivar_NL_tilting(mu1, mu2n, s1, s2n, rho, c1, c2, a, b):
    """Independent E[(x1-c1)^a (x2-c2)^b] via exponential tilting (a DIFFERENT derivation than the Fortran sum).

    Expand (x2-c2)^b in raw lognormal moments; each exp(q·ln x2) tilts the joint normal, factoring out the
    lognormal raw moment exp(qμ2n+½q²σ2n²) and shifting x1's mean by q·ρ·s1·s2n. The remaining factor is the
    central a-moment of a normal with that shifted mean about c1, computed via _normal_raw_central.
    """
    out = 0.0
    for q in range(b + 1):
        raw_ln = math.exp(mu2n * q + 0.5 * s2n ** 2 * q ** 2)
        shifted_mu1 = mu1 + q * rho * s1 * s2n
        out += math.comb(b, q) * (-c2) ** (b - q) * raw_ln * _normal_raw_central(shifted_mu1, s1, c1, a)
    return out


def test_bivar_NL_closed_form_and_mc():
    mu1, mu2n, s1, s2n, rho = 0.8, -0.3, 0.6, 0.45, 0.35
    c1, c2 = 0.5, float(np.exp(mu2n + 0.5 * s2n ** 2))
    rng = np.random.default_rng(33)
    N = 12_000_000
    cov = [[s1 ** 2, rho * s1 * s2n], [rho * s1 * s2n, s2n ** 2]]
    xy = rng.multivariate_normal([mu1, mu2n], cov, size=N)
    x1 = xy[:, 0]
    x2 = np.exp(xy[:, 1])
    for a, b in ((1, 1), (2, 1), (2, 2), (3, 1), (0, 2), (3, 2)):
        got = float(bivar_NL_int_PDF_comp_all_MM(mu1, mu2n, s1, s2n, rho, c1, c2, a, b))
        ref = _bivar_NL_tilting(mu1, mu2n, s1, s2n, rho, c1, c2, a, b)
        rel = abs(got - ref) / (abs(ref) + 1e-30)
        assert rel < 1e-12, f"bivar_NL a={a},b={b}: closed-form rel {rel:.2e}"
        mc = np.mean((x1 - c1) ** a * (x2 - c2) ** b)
        scale = np.mean(np.abs((x1 - c1) ** a * (x2 - c2) ** b)) + 1e-30
        rel_mc = abs(got - mc) / scale
        assert rel_mc < 1e-2, f"bivar_NL a={a},b={b}: MC rel {rel_mc:.2e}"
    print("  bivar_NL (6 (a,b) pairs) vs tilting decomposition (<1e-12) + 12M MC (<1e-2)  PASS")


def _comp_eq_branch_ref(mu_x, mu_hm, mu_hm_n, s_x, s_hm, s_hm_n, rho, pf, xm, hmm, xtol, hmtol, a, b):
    """Independent NumPy evaluation of the 4 explicit Fortran branches of bivar_NL_x_hm_all_MM_comp_eq, using the
    raw-moment helpers (univar_N = normal central moment, univar_L = lognormal central moment)."""
    out = (1.0 - pf) * (-hmm) ** b
    if s_x <= xtol and s_hm <= hmtol:
        return (mu_x - xm) ** a * (pf * (mu_hm - hmm) ** b + out)
    if s_x <= xtol:
        return (mu_x - xm) ** a * (pf * _lognormal_raw_central(mu_hm_n, s_hm_n, hmm, b) + out)
    uN = _normal_raw_central(mu_x, s_x, xm, a)
    if s_hm <= hmtol:
        return (pf * (mu_hm - hmm) ** b + out) * uN
    return pf * _bivar_NL_tilting(mu_x, mu_hm_n, s_x, s_hm_n, rho, xm, hmm, a, b) + out * uN


def test_comp_eq_branch_selection():
    """All 4 branches of bivar_NL_x_hm_all_MM_comp_eq vs the independent per-branch NumPy formulas."""
    base = dict(mu_x=0.7, mu_hm=2.0e-4, mu_hm_n=-8.6, s_hm_n=0.5, rho=0.3,
                pf=0.7, xm=0.4, hmm=1.5e-4, xtol=1e-2, hmtol=1e-8)
    s_hm_big = 1.0e-4   # > hmtol
    regimes = {
        "both_const": dict(s_x=1e-3, s_hm=1e-9),
        "x_const":    dict(s_x=1e-3, s_hm=s_hm_big),
        "hm_const":   dict(s_x=0.6, s_hm=1e-9),
        "both_vary":  dict(s_x=0.6, s_hm=s_hm_big),
    }
    for a, b in ((1, 1), (2, 2), (3, 1)):
        for name, sig in regimes.items():
            got = float(bivar_NL_x_hm_all_MM_comp_eq(
                base['mu_x'], base['mu_hm'], base['mu_hm_n'], sig['s_x'], sig['s_hm'], base['s_hm_n'],
                base['rho'], base['pf'], base['xm'], base['hmm'], base['xtol'], base['hmtol'], a, b))
            ref = _comp_eq_branch_ref(base['mu_x'], base['mu_hm'], base['mu_hm_n'], sig['s_x'], sig['s_hm'],
                                      base['s_hm_n'], base['rho'], base['pf'], base['xm'], base['hmm'],
                                      base['xtol'], base['hmtol'], a, b)
            rel = abs(got - ref) / (abs(ref) + 1e-30)
            assert rel < 1e-12, f"comp_eq {name} a={a},b={b}: rel {rel:.2e}"
    print("  bivar_NL_x_hm_all_MM_comp_eq: all 4 branches × 3 (a,b) vs independent NumPy (<1e-12)  PASS")


def test_xp_a_hmpb_monte_carlo():
    """Full 2-component <x'^a hm'^b> vs a Monte-Carlo of the in/out-precip mixture (general both-vary regime)."""
    P = dict(mu_x_1=0.6, mu_x_2=-0.4, mu_hm_1_n=-8.4, mu_hm_2_n=-7.9,
             s_x_1=0.5, s_x_2=0.7, s_hm_1_n=0.5, s_hm_2_n=0.45, rho1=0.3, rho2=-0.2,
             mixt=0.45, pf1=0.8, pf2=0.5)
    # normal-space hm means/stds (> tol so both-vary branch is taken)
    def lm(mu_n, s_n):
        m = math.exp(mu_n + 0.5 * s_n ** 2)
        return m, max(m * math.sqrt(math.exp(s_n ** 2) - 1), 1e-9)
    P['mu_hm_1'], P['s_hm_1'] = lm(P['mu_hm_1_n'], P['s_hm_1_n'])
    P['mu_hm_2'], P['s_hm_2'] = lm(P['mu_hm_2_n'], P['s_hm_2_n'])
    xm = P['mixt'] * P['mu_x_1'] + (1 - P['mixt']) * P['mu_x_2']
    hmm = P['mixt'] * P['pf1'] * P['mu_hm_1'] + (1 - P['mixt']) * P['pf2'] * P['mu_hm_2']
    xtol, hmtol = 1e-2, 1e-12

    rng = np.random.default_rng(44)
    N = 16_000_000
    comp1 = rng.random(N) < P['mixt']
    x = np.empty(N); hm = np.zeros(N)
    for i, mask in ((1, comp1), (2, ~comp1)):
        idx = np.where(mask)[0]
        n = idx.size
        sx, shn, mxn = P[f's_x_{i}'], P[f's_hm_{i}_n'], P[f'mu_hm_{i}_n']
        mu_x_i, rho = P[f'mu_x_{i}'], P[f'rho{i}']
        in_p = rng.random(n) < P[f'pf{i}']
        # in-precip: bivariate normal (x, y) with corr rho; hm = exp(y)
        ip = idx[in_p]
        cov = [[sx ** 2, rho * sx * shn], [rho * sx * shn, shn ** 2]]
        xy = rng.multivariate_normal([mu_x_i, mxn], cov, size=ip.size)
        x[ip] = xy[:, 0]; hm[ip] = np.exp(xy[:, 1])
        # out-of-precip: x normal, hm = 0
        op = idx[~in_p]
        x[op] = rng.normal(mu_x_i, sx, op.size)
    for a, b in ((1, 1), (2, 1), (2, 2), (1, 2)):
        got = float(xp_a_hmpb_integrals_all_MM(
            P['mu_x_1'], P['mu_x_2'], P['mu_hm_1'], P['mu_hm_2'], P['mu_hm_1_n'], P['mu_hm_2_n'],
            P['s_x_1'], P['s_x_2'], P['s_hm_1'], P['s_hm_2'], P['s_hm_1_n'], P['s_hm_2_n'],
            P['rho1'], P['rho2'], P['mixt'], P['pf1'], P['pf2'], xm, hmm, xtol, hmtol, a, b))
        mc = np.mean((x - xm) ** a * (hm - hmm) ** b)
        scale = np.mean(np.abs((x - xm) ** a * (hm - hmm) ** b)) + 1e-30
        rel = abs(got - mc) / scale
        assert rel < 1e-2, f"xp_a_hmpb a={a},b={b}: MC rel {rel:.2e}"
    print("  xp_a_hmpb_integrals_all_MM (4 (a,b)) vs 16M-sample 2-component MC (<1e-2)  PASS")


def _xphmp_ref_4branch(p):
    """Literal NumPy transcription of the Fortran xphmp_integral_covar 4-way branch (independent of the port)."""
    s_x1, s_x2, s_h1, s_h2, xt, ht = p['s_x1'], p['s_x2'], p['s_h1'], p['s_h2'], p['xt'], p['ht']
    c1 = (s_x1 <= xt) or (s_h1 <= ht)
    c2 = (s_x2 <= xt) or (s_h2 <= ht)
    t1 = (p['mu_x1'] - p['xm']) * p['mu_h1'] + (0.0 if c1 else p['c_xh1'] * s_x1 * s_h1)
    t2 = (p['mu_x2'] - p['xm']) * p['mu_h2'] + (0.0 if c2 else p['c_xh2'] * s_x2 * s_h2)
    return p['mixt'] * p['pf1'] * t1 + (1 - p['mixt']) * p['pf2'] * t2


def test_xphmp_covar_branches_and_mc():
    p = dict(mu_x1=0.6, mu_x2=-0.4, mu_h1=2.0e-4, mu_h2=3.0e-4, s_x1=0.5, s_x2=0.7,
             s_h1=1.5e-4, s_h2=1.0e-4, c_xh1=0.4, c_xh2=-0.3, mixt=0.45, pf1=0.8, pf2=0.5,
             xt=1e-2, ht=1e-12)
    call = lambda pp: float(xphmp_integral_covar(
        pp['mu_x1'], pp['mu_x2'], pp['mu_h1'], pp['mu_h2'], pp['s_x1'], pp['s_x2'], pp['s_h1'], pp['s_h2'],
        pp['c_xh1'], pp['c_xh2'], pp['mixt'], pp['pf1'], pp['pf2'], pp['xm'], pp['xt'], pp['ht']))
    # xm chosen as the overall mean (so the MC covariance is well-defined); branch test independent of xm value.
    p['xm'] = p['mixt'] * p['mu_x1'] + (1 - p['mixt']) * p['mu_x2']
    # 4 branch regimes via tiny sigmas
    for s_x1, s_h1, s_x2, s_h2 in ((1e-3, 1e-13, 1e-3, 1e-13), (1e-3, 1e-13, 0.7, 1e-4),
                                   (0.5, 1.5e-4, 1e-3, 1e-13), (0.5, 1.5e-4, 0.7, 1e-4)):
        q = {**p, 's_x1': s_x1, 's_h1': s_h1, 's_x2': s_x2, 's_h2': s_h2}
        rel = abs(call(q) - _xphmp_ref_4branch(q)) / (abs(_xphmp_ref_4branch(q)) + 1e-30)
        assert rel < 1e-12, f"xphmp branch rel {rel:.2e}"
    # Monte-Carlo of the general regime: covariance is moment-based → bivariate-normal in-precip model is exact.
    rng = np.random.default_rng(55)
    N = 16_000_000
    comp1 = rng.random(N) < p['mixt']
    x = np.empty(N); hm = np.zeros(N)
    for i, (mux, muh, sx, sh, c, pf, mask) in enumerate((
            (p['mu_x1'], p['mu_h1'], p['s_x1'], p['s_h1'], p['c_xh1'], p['pf1'], comp1),
            (p['mu_x2'], p['mu_h2'], p['s_x2'], p['s_h2'], p['c_xh2'], p['pf2'], ~comp1))):
        idx = np.where(mask)[0]
        in_p = rng.random(idx.size) < pf
        ip = idx[in_p]
        cov = [[sx ** 2, c * sx * sh], [c * sx * sh, sh ** 2]]
        xy = rng.multivariate_normal([mux, muh], cov, size=ip.size)
        x[ip] = xy[:, 0]; hm[ip] = xy[:, 1]
        op = idx[~in_p]
        x[op] = rng.normal(mux, sx, op.size)   # out-of-precip hm stays 0
    hmm = np.mean(hm)
    got = call(p)
    integrand = (x - p['xm']) * (hm - hmm)
    mc = np.mean(integrand)
    # Normalize by the mean absolute integrand (well-conditioned under the inter-component cancellation that
    # makes the bare covariance small relative to the product's spread).
    scale = np.mean(np.abs(integrand)) + 1e-30
    rel = abs(got - mc) / scale
    assert rel < 3e-3, f"xphmp MC rel {rel:.2e}"
    print(f"  xphmp_integral_covar: 4 branches vs literal NumPy (<1e-12) + 16M MC covariance (rel {rel:.1e})  PASS")


def test_hmxphmyp_covar_branches_and_mc():
    # two hydrometeors hmx, hmy; covariance over the 2-component in-precip mixture.
    P = dict(mu_x1=2.0e-4, mu_x2=3.0e-4, mu_y1=5.0e5, mu_y2=7.0e5, s_x1=1.5e-4, s_x2=1.0e-4,
             s_y1=2.0e5, s_y2=3.0e5, c1=0.3, c2=-0.25, mixt=0.4, pf1=0.7, pf2=0.55,
             xt=1e-12, yt=1e-3)
    hmx_mean = P['mixt'] * P['pf1'] * P['mu_x1'] + (1 - P['mixt']) * P['pf2'] * P['mu_x2']
    hmy_mean = P['mixt'] * P['pf1'] * P['mu_y1'] + (1 - P['mixt']) * P['pf2'] * P['mu_y2']
    got = float(hmxphmyp_integral_covar(
        P['mu_x1'], P['mu_x2'], P['mu_y1'], P['mu_y2'], P['s_x1'], P['s_x2'], P['s_y1'], P['s_y2'],
        P['c1'], P['c2'], P['mixt'], P['pf1'], P['pf2'], hmx_mean, hmy_mean, P['xt'], P['yt']))
    # independent literal formula
    t1 = P['mu_x1'] * P['mu_y1'] + P['c1'] * P['s_x1'] * P['s_y1']
    t2 = P['mu_x2'] * P['mu_y2'] + P['c2'] * P['s_x2'] * P['s_y2']
    ref = P['mixt'] * P['pf1'] * t1 + (1 - P['mixt']) * P['pf2'] * t2 - hmx_mean * hmy_mean
    assert abs(got - ref) / (abs(ref) + 1e-30) < 1e-12, "hmxphmyp closed-form mismatch"
    # Monte-Carlo
    rng = np.random.default_rng(66)
    N = 16_000_000
    comp1 = rng.random(N) < P['mixt']
    hx = np.zeros(N); hy = np.zeros(N)
    for mux, muy, sx, sy, c, pf, mask in (
            (P['mu_x1'], P['mu_y1'], P['s_x1'], P['s_y1'], P['c1'], P['pf1'], comp1),
            (P['mu_x2'], P['mu_y2'], P['s_x2'], P['s_y2'], P['c2'], P['pf2'], ~comp1)):
        idx = np.where(mask)[0]
        in_p = rng.random(idx.size) < pf
        ip = idx[in_p]
        cov = [[sx ** 2, c * sx * sy], [c * sx * sy, sy ** 2]]
        xy = rng.multivariate_normal([mux, muy], cov, size=ip.size)
        hx[ip] = xy[:, 0]; hy[ip] = xy[:, 1]   # out-of-precip both stay 0
    integrand = (hx - hx.mean()) * (hy - hy.mean())
    mc = np.mean(integrand)
    scale = np.mean(np.abs(integrand)) + 1e-30
    rel = abs(got - mc) / scale
    assert rel < 3e-3, f"hmxphmyp MC rel {rel:.2e}"
    print(f"  hmxphmyp_integral_covar: closed-form (<1e-12) + 16M MC covariance (rel {rel:.1e})  PASS")


def test_differentiable():
    gN = jax.grad(lambda s: univar_N_int_PDF_comp_all_MM(1.3, s, 0.9, 4))(jnp.asarray(0.7))
    gL = jax.grad(lambda m: univar_L_int_PDF_comp_all_MM(m, 0.5, 1.0, 3))(jnp.asarray(-0.4))
    gB = jax.grad(lambda r: bivar_NL_int_PDF_comp_all_MM(0.8, -0.3, 0.6, 0.45, r, 0.5, 1.0, 2, 2))(jnp.asarray(0.35))
    # grad through the jnp.where branch selection of the full 2-component moment (both-vary regime)
    gX = jax.grad(lambda sx: xp_a_hmpb_integrals_all_MM(
        0.6, -0.4, 2.3e-4, 3.7e-4, -8.4, -7.9, sx, 0.7, 2.4e-4, 3.0e-4, 0.5, 0.45,
        0.3, -0.2, 0.45, 0.8, 0.5, 0.1, 2.0e-4, 1e-2, 1e-12, 2, 2))(jnp.asarray(0.5))
    assert all(np.isfinite(float(g)) for g in (gN, gL, gB, gX)), "non-finite grad"
    print(f"  jax.grad univar_N/dσ={float(gN):+.3e}, univar_L/dμn={float(gL):+.3e}, "
          f"bivar_NL/dρ={float(gB):+.3e}, xp_a_hmpb/dσx={float(gX):+.3e}: finite  PASS")


def main():
    print("test_mixed_moment_pdf_integrals:")
    for t in (test_univar_N_closed_form_and_mc, test_univar_L_closed_form_and_mc,
              test_bivar_NL_closed_form_and_mc, test_comp_eq_branch_selection,
              test_xp_a_hmpb_monte_carlo, test_xphmp_covar_branches_and_mc,
              test_hmxphmyp_covar_branches_and_mc, test_differentiable):
        t()
    print("All mixed_moment_PDF_integrals checks PASSED")


if __name__ == "__main__":
    main()

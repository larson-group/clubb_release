#!/usr/bin/env python3
"""test_kk_upscaled_variances.py — validate variance_KK_mvr (variance of the KK rain mean volume radius).

Two independent oracles (no f2py wrapper for this routine):
  1. Closed-form: for the non-degenerate branch, the bivariate-lognormal moment is
     `E[rr^a Nr^b] = exp(a*mu_rr_n + b*mu_Nr_n + 0.5*(a^2 s_rr^2 + b^2 s_Nr^2 + 2ab*rho*s_rr*s_Nr))`,
     so `Var(R_vr) = E[R^2]-E[R]^2` is computed from scratch (not via the ported helpers).
  2. Monte-Carlo: sample the 2-component, in-precip mixture and take the sample variance of R_vr.
Plus a finite `jax.grad`.
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

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_variances import variance_KK_mvr

ALPHA, BETA, COEF = 1.0 / 3.0, -1.0 / 3.0, 3.6e-3


def _lognormal_moment(a, b, mu_rr_n, mu_Nr_n, s_rr_n, s_Nr_n, rho):
    """Closed-form E[rr^a Nr^b] for a bivariate lognormal in log-space (mu_*_n, sigma_*_n, rho)."""
    return np.exp(a * mu_rr_n + b * mu_Nr_n
                  + 0.5 * (a ** 2 * s_rr_n ** 2 + b ** 2 * s_Nr_n ** 2 + 2 * a * b * rho * s_rr_n * s_Nr_n))


def _make_case():
    # log-space moments (component 1, 2), all non-degenerate so bivar_LL_mean_eq uses the general branch.
    p = dict(mu_rr_1_n=-9.0, mu_rr_2_n=-8.0, mu_Nr_1_n=11.0, mu_Nr_2_n=12.0,
             s_rr_1_n=0.5, s_rr_2_n=0.4, s_Nr_1_n=0.6, s_Nr_2_n=0.5, rho1=0.25, rho2=-0.15,
             mixt_frac=0.4, pf1=0.8, pf2=0.6)
    # normal-space moments (only used by bivar_LL_mean_eq's degenerate dispatch; set > tol so the general
    # branch is taken — plausible lognormal means/stds).
    def mom(mu_n, s_n):
        m = np.exp(mu_n + 0.5 * s_n ** 2)
        sd = m * np.sqrt(np.exp(s_n ** 2) - 1.0)
        return m, max(sd, 1e-6)
    p['mu_rr_1'], p['s_rr_1'] = mom(p['mu_rr_1_n'], p['s_rr_1_n'])
    p['mu_rr_2'], p['s_rr_2'] = mom(p['mu_rr_2_n'], p['s_rr_2_n'])
    p['mu_Nr_1'], p['s_Nr_1'] = mom(p['mu_Nr_1_n'], p['s_Nr_1_n'])
    p['mu_Nr_2'], p['s_Nr_2'] = mom(p['mu_Nr_2_n'], p['s_Nr_2_n'])
    return p


def _ref_mean_and_var(p):
    """Independent E[R_vr] and Var(R_vr) from the closed-form moments."""
    def comp(a, b, i):
        return _lognormal_moment(a, b, p[f'mu_rr_{i}_n'], p[f'mu_Nr_{i}_n'],
                                 p[f's_rr_{i}_n'], p[f's_Nr_{i}_n'], p[f'rho{i}'])
    w1, w2 = p['mixt_frac'] * p['pf1'], (1 - p['mixt_frac']) * p['pf2']
    e_r = COEF * (w1 * comp(ALPHA, BETA, 1) + w2 * comp(ALPHA, BETA, 2))
    e_r2 = COEF ** 2 * (w1 * comp(2 * ALPHA, 2 * BETA, 1) + w2 * comp(2 * ALPHA, 2 * BETA, 2))
    return e_r, e_r2 - e_r ** 2


def _call_jax(p, kk_mean_vol_rad):
    return float(variance_KK_mvr(
        p['mu_rr_1'], p['mu_rr_2'], p['mu_Nr_1'], p['mu_Nr_2'],
        p['mu_rr_1_n'], p['mu_rr_2_n'], p['mu_Nr_1_n'], p['mu_Nr_2_n'],
        p['s_rr_1'], p['s_rr_2'], p['s_Nr_1'], p['s_Nr_2'],
        p['s_rr_1_n'], p['s_rr_2_n'], p['s_Nr_1_n'], p['s_Nr_2_n'],
        p['rho1'], p['rho2'], kk_mean_vol_rad, COEF, p['mixt_frac'], p['pf1'], p['pf2']))


def test_closed_form():
    p = _make_case()
    e_r, var_ref = _ref_mean_and_var(p)
    got = _call_jax(p, e_r)   # pass the independently-computed mean as KK_mean_vol_rad
    rel = abs(got - var_ref) / (abs(var_ref) + 1e-30)
    assert rel < 1e-12, f"variance_KK_mvr vs closed-form rel {rel:.2e}"
    print(f"  variance_KK_mvr vs closed-form lognormal moments: rel {rel:.2e}  PASS")


def test_monte_carlo():
    p = _make_case()
    e_r, var_ref = _ref_mean_and_var(p)
    got = _call_jax(p, e_r)
    rng = np.random.default_rng(7)
    N = 4_000_000
    comp1 = rng.random(N) < p['mixt_frac']
    R = np.zeros(N)
    for i, mask in ((1, comp1), (2, ~comp1)):
        idx = np.where(mask)[0]
        n = idx.size
        prec = rng.random(n) < p[f'pf{i}']
        sub = idx[prec]
        cov = np.array([[p[f's_rr_{i}_n'] ** 2, p[f'rho{i}'] * p[f's_rr_{i}_n'] * p[f's_Nr_{i}_n']],
                        [p[f'rho{i}'] * p[f's_rr_{i}_n'] * p[f's_Nr_{i}_n'], p[f's_Nr_{i}_n'] ** 2]])
        ln = rng.multivariate_normal([p[f'mu_rr_{i}_n'], p[f'mu_Nr_{i}_n']], cov, size=sub.size)
        rr, Nr = np.exp(ln[:, 0]), np.exp(ln[:, 1])
        R[sub] = COEF * rr ** ALPHA * Nr ** BETA
    var_mc = R.var()
    rel = abs(got - var_mc) / (abs(var_mc) + 1e-30)
    assert rel < 5e-3, f"variance_KK_mvr vs Monte-Carlo rel {rel:.2e}"
    print(f"  variance_KK_mvr vs Monte-Carlo ({N//10**6}M samples): rel {rel:.2e}  PASS")


def test_differentiable():
    p = _make_case()
    e_r, _ = _ref_mean_and_var(p)
    g = jax.grad(lambda s: variance_KK_mvr(
        p['mu_rr_1'], p['mu_rr_2'], p['mu_Nr_1'], p['mu_Nr_2'],
        p['mu_rr_1_n'], p['mu_rr_2_n'], p['mu_Nr_1_n'], p['mu_Nr_2_n'],
        p['s_rr_1'], p['s_rr_2'], p['s_Nr_1'], p['s_Nr_2'],
        s, p['s_rr_2_n'], p['s_Nr_1_n'], p['s_Nr_2_n'],
        p['rho1'], p['rho2'], e_r, COEF, p['mixt_frac'], p['pf1'], p['pf2']))(jnp.asarray(p['s_rr_1_n']))
    assert np.isfinite(float(g)), "non-finite grad through variance_KK_mvr"
    print(f"  jax.grad(variance_KK_mvr)/d(sigma_rr_1_n) = {float(g):+.4e}: finite  PASS")


def main():
    print("test_kk_upscaled_variances:")
    for t in (test_closed_form, test_monte_carlo, test_differentiable):
        t()
    print("All KK_upscaled_variances checks PASSED")


if __name__ == "__main__":
    main()

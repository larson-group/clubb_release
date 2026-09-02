#!/usr/bin/env python3
"""test_pdf_integrals_all_mm.py — validate trivar_NNL_MM (KK all-mixed-moment NNL integral).

The integral is  < (x1-<x1>)^a (x2^alpha x3^beta - M)^b >  over the joint (x1,x2,x3) with x2 RESTRICTED to
x2 > 0 (the full normal PDF integrated on the half-line, not renormalized), x1 normal, x3 lognormal. No f2py
oracle exists (and the Fortran's own D_v is known low-accuracy), so two independent oracles:

  1. Base case a=b=0: the integral collapses to the x2>0 probability mass Phi(mu_x2/sigma_x2) — a closed form
     independent of the whole D_v/Gamma machinery. (Derived: (1/sqrt 2pi) exp(-s²/4) Gamma(1) D_{-1}(-s) = Phi(s).)
  2. General (a,b,alpha,beta): a truncated-domain Monte-Carlo — sample the trivariate normal (x1, x2, ln x3),
     keep the indicator [x2>0], average (x1-<x1>)^a (x2^alpha x3^beta - M)^b with the SAME M passed to the port.
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

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_all_MM import (
    trivar_NNL_MM, trivar_NNL_MM_const_x1, trivar_NNL_MM_const_x2, trivar_NNL_MM_const_x1x2,
    quadrivar_NNLL_MM, quadrivar_NNLL_MM_const_x1, quadrivar_NNLL_MM_const_x2,
    quadrivar_NNLL_MM_const_x1x2)


def _phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def test_base_case_probability_mass():
    """a=b=0 → Phi(mu_x2/sigma_x2), independent of alpha/beta/the other moments."""
    for mu_x2, sigma_x2 in ((0.3, 0.5), (-0.2, 0.8), (1.0, 0.4), (0.0, 1.0)):
        got = float(trivar_NNL_MM(
            mu_x1=0.5, mu_x2=mu_x2, mu_x3_n=-0.4, sigma_x1=0.6, sigma_x2=sigma_x2, sigma_x3_n=0.5,
            rho_x1x2=0.3, rho_x1x3_n=0.2, rho_x2x3_n=-0.1, x1_mean=0.4, x2_alpha_x3_beta_mean=1.3,
            alpha_exp=1.5, beta_exp=0.5, a_exp=0, b_exp=0))
        ref = _phi(mu_x2 / sigma_x2)
        rel = abs(got - ref) / (abs(ref) + 1e-30)
        assert rel < 1e-9, f"base case mu_x2={mu_x2}: rel {rel:.2e} (got {got}, Phi {ref})"
    print("  trivar_NNL_MM(a=0,b=0) vs Phi(mu_x2/sigma_x2) over 4 configs: rel <1e-9  PASS")


def _mc(mu_x1, mu_x2, mu_x3_n, s_x1, s_x2, s_x3n, r12, r13n, r23n, x1_mean, M, alpha, beta, a, b, N, seed):
    rng = np.random.default_rng(seed)
    cov = np.array([
        [s_x1 ** 2, r12 * s_x1 * s_x2, r13n * s_x1 * s_x3n],
        [r12 * s_x1 * s_x2, s_x2 ** 2, r23n * s_x2 * s_x3n],
        [r13n * s_x1 * s_x3n, r23n * s_x2 * s_x3n, s_x3n ** 2]])
    xyz = rng.multivariate_normal([mu_x1, mu_x2, mu_x3_n], cov, size=N)
    x1, x2, y3 = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    x3 = np.exp(y3)
    mask = x2 > 0.0
    integrand = np.where(mask, (x1 - x1_mean) ** a * (np.where(mask, x2, 1.0) ** alpha * x3 ** beta - M) ** b, 0.0)
    return integrand


def test_general_monte_carlo():
    cfg = dict(mu_x1=0.5, mu_x2=0.35, mu_x3_n=-0.4, s_x1=0.6, s_x2=0.5, s_x3n=0.45,
               r12=0.3, r13n=0.2, r23n=-0.15, x1_mean=0.45, alpha=1.5, beta=0.5)
    # Choose M = conditional mean of x2^alpha x3^beta over x2>0, from a big sample (used identically in both).
    pre = _mc(cfg['mu_x1'], cfg['mu_x2'], cfg['mu_x3_n'], cfg['s_x1'], cfg['s_x2'], cfg['s_x3n'],
              cfg['r12'], cfg['r13n'], cfg['r23n'], cfg['x1_mean'], 0.0, cfg['alpha'], cfg['beta'], 0, 1,
              20_000_000, seed=1)  # a=0,b=1,M=0 → integrand = [x2>0] x2^a x3^b ; mean = partial mean
    mass = _mc(cfg['mu_x1'], cfg['mu_x2'], cfg['mu_x3_n'], cfg['s_x1'], cfg['s_x2'], cfg['s_x3n'],
               cfg['r12'], cfg['r13n'], cfg['r23n'], cfg['x1_mean'], 0.0, cfg['alpha'], cfg['beta'], 0, 0,
               20_000_000, seed=1)
    M = float(pre.sum() / mass.sum())   # E[x2^a x3^b | x2>0]
    for a, b in ((1, 1), (2, 1), (1, 2), (2, 2)):
        got = float(trivar_NNL_MM(
            cfg['mu_x1'], cfg['mu_x2'], cfg['mu_x3_n'], cfg['s_x1'], cfg['s_x2'], cfg['s_x3n'],
            cfg['r12'], cfg['r13n'], cfg['r23n'], cfg['x1_mean'], M, cfg['alpha'], cfg['beta'], a, b))
        integ = _mc(cfg['mu_x1'], cfg['mu_x2'], cfg['mu_x3_n'], cfg['s_x1'], cfg['s_x2'], cfg['s_x3n'],
                    cfg['r12'], cfg['r13n'], cfg['r23n'], cfg['x1_mean'], M, cfg['alpha'], cfg['beta'],
                    a, b, 20_000_000, seed=2)
        mc = integ.mean()
        scale = np.mean(np.abs(integ)) + 1e-30
        rel = abs(got - mc) / scale
        assert rel < 5e-3, f"trivar_NNL_MM a={a},b={b}: MC rel {rel:.2e} (got {got}, mc {mc})"
    print("  trivar_NNL_MM (4 (a,b)) vs truncated-domain 20M MC: rel <5e-3  PASS")


def test_const_x1_limit_and_mc():
    """const_x1 vs (a) the general form with sigma_x1->0, and (b) a direct MC with x1 fixed."""
    c = dict(mu_x1=0.5, mu_x2=0.35, mu_x3_n=-0.4, s_x2=0.5, s_x3n=0.45, r23n=-0.15,
             x1_mean=0.45, M=1.2, alpha=1.5, beta=0.5)
    for a, b in ((1, 1), (2, 1), (2, 2)):
        const = float(trivar_NNL_MM_const_x1(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], c['s_x2'], c['s_x3n'],
                                             c['r23n'], c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b))
        # (a) limit consistency: general with tiny sigma_x1 (arbitrary x1 correlations vanish with it)
        gen = float(trivar_NNL_MM(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], 1e-8, c['s_x2'], c['s_x3n'],
                                  0.3, 0.2, c['r23n'], c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b))
        assert abs(const - gen) / (abs(const) + 1e-30) < 1e-5, f"const_x1 limit a={a},b={b}"
        # (b) MC with x1 fixed (sigma_x1=0)
        integ = _mc(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], 0.0, c['s_x2'], c['s_x3n'], 0.0, 0.0, c['r23n'],
                    c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b, 20_000_000, seed=3)
        rel = abs(const - integ.mean()) / (np.mean(np.abs(integ)) + 1e-30)
        assert rel < 5e-3, f"const_x1 MC a={a},b={b}: rel {rel:.2e}"
    print("  trivar_NNL_MM_const_x1: vs general(sigma_x1->0) <1e-5 + 20M MC <5e-3  PASS")


def test_const_x2_mc():
    """const_x2 (x2 fixed at mu_x2>0, x1 FULL normal, x3 lognormal) vs direct MC (no half-line truncation)."""
    c = dict(mu_x1=0.5, mu_x2=0.4, mu_x3_n=-0.4, s_x1=0.6, s_x3n=0.45, r13n=0.25,
             x1_mean=0.45, M=1.1, alpha=1.5, beta=0.5)
    for a, b in ((1, 1), (2, 1), (2, 2)):
        const = float(trivar_NNL_MM_const_x2(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], c['s_x1'], c['s_x3n'],
                                             c['r13n'], c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b))
        # MC: x2 fixed (sigma_x2=0, mu_x2>0 so [x2>0]=1 always), x1 full normal correlated w/ ln x3 via r13n
        integ = _mc(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], c['s_x1'], 0.0, c['s_x3n'], 0.0, c['r13n'], 0.0,
                    c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b, 20_000_000, seed=4)
        rel = abs(const - integ.mean()) / (np.mean(np.abs(integ)) + 1e-30)
        assert rel < 5e-3, f"const_x2 MC a={a},b={b}: rel {rel:.2e}"
    print("  trivar_NNL_MM_const_x2: vs 20M MC (x2 fixed) <5e-3  PASS")


def test_const_x1x2_limit_and_mc():
    """const_x1x2 (only x3 lognormal varies) vs (a) const_x2 with sigma_x1->0 and (b) direct MC."""
    c = dict(mu_x1=0.5, mu_x2=0.4, mu_x3_n=-0.4, s_x3n=0.45, x1_mean=0.45, M=1.1, alpha=1.5, beta=0.5)
    for a, b in ((1, 1), (2, 1), (3, 2)):
        const = float(trivar_NNL_MM_const_x1x2(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], c['s_x3n'],
                                               c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b))
        # (a) limit: const_x2 with sigma_x1 -> 0
        cx2 = float(trivar_NNL_MM_const_x2(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], 1e-9, c['s_x3n'], 0.2,
                                           c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b))
        assert abs(const - cx2) / (abs(const) + 1e-30) < 1e-7, f"const_x1x2 limit a={a},b={b}"
        # (b) MC: x1, x2 fixed; only x3 lognormal
        integ = _mc(c['mu_x1'], c['mu_x2'], c['mu_x3_n'], 0.0, 0.0, c['s_x3n'], 0.0, 0.0, 0.0,
                    c['x1_mean'], c['M'], c['alpha'], c['beta'], a, b, 8_000_000, seed=5)
        rel = abs(const - integ.mean()) / (np.mean(np.abs(integ)) + 1e-30)
        assert rel < 5e-3, f"const_x1x2 MC a={a},b={b}: rel {rel:.2e}"
    print("  trivar_NNL_MM_const_x1x2: vs const_x2(sigma_x1->0) <1e-7 + 8M MC <5e-3  PASS")


def _normal_central(mu, sigma, c, n):
    """n-th central moment E[(X-c)^n] for X~N(mu,sigma) (independent: raw double-factorial expansion)."""
    cm = lambda k: 0.0 if k % 2 else sigma ** k * math.factorial(k) / (2 ** (k // 2) * math.factorial(k // 2))
    return sum(math.comb(n, k) * (mu - c) ** (n - k) * cm(k) for k in range(n + 1))


def _Q():
    return dict(mu_x1=0.5, mu_x2=-0.3, mu_x3_n=-0.4, mu_x4_n=-0.6, s_x1=0.6, s_x2=0.5, s_x3n=0.45, s_x4n=0.5,
                r12=0.25, r13n=0.2, r14n=-0.1, r23n=-0.15, r24n=0.1, r34n=0.3,
                x1_mean=0.4, M=1.3, alpha=1.5, beta=0.5, gamma=0.4)


def _call_quad(q, a, b, **over):
    p = {**q, **over}
    return float(quadrivar_NNLL_MM(
        p['mu_x1'], p['mu_x2'], p['mu_x3_n'], p['mu_x4_n'], p['s_x1'], p['s_x2'], p['s_x3n'], p['s_x4n'],
        p['r12'], p['r13n'], p['r14n'], p['r23n'], p['r24n'], p['r34n'], p['x1_mean'], p['M'],
        p['alpha'], p['beta'], p['gamma'], a, b))


def test_quadrivar_base_and_psum():
    q = _Q()
    # (1) a=b=0 -> Phi(-mu_x2/sigma_x2) (x2<0 mass), independent of all else.
    for mx2, sx2 in ((-0.3, 0.5), (0.2, 0.8), (0.0, 1.0)):
        got = _call_quad(q, 0, 0, mu_x2=mx2, s_x2=sx2)
        ref = _phi(-mx2 / sx2)
        assert abs(got - ref) / (abs(ref) + 1e-30) < 1e-9, f"quad base mu_x2={mx2}"
    # (2) a=2,b=0,rho12=0 -> (2nd central moment of x1 about x1_mean) * Phi(-mu_x2/sigma_x2)
    for a in (1, 2, 3):
        got = _call_quad(q, a, 0, r12=0.0)
        ref = _normal_central(q['mu_x1'], q['s_x1'], q['x1_mean'], a) * _phi(-q['mu_x2'] / q['s_x2'])
        assert abs(got - ref) / (abs(ref) + 1e-30) < 1e-9, f"quad p-sum a={a}"
    print("  quadrivar_NNLL_MM: a=b=0 -> Phi(-mu_x2/sigma_x2) + (a,b=0,rho12=0) -> central-moment*Phi: <1e-9  PASS")


def _quad_complex_mc(q, a, b, N, seed, s_over=None):
    """Truncated (x2<0) complex-branch MC of <(x1-<x1>)^a (x2^α x3^β x4^γ - M)^b>; s_over overrides sigmas."""
    s = {**dict(s_x1=q['s_x1'], s_x2=q['s_x2'], s_x3n=q['s_x3n'], s_x4n=q['s_x4n']), **(s_over or {})}
    sig = [s['s_x1'], s['s_x2'], s['s_x3n'], s['s_x4n']]
    R = np.array([[1, q['r12'], q['r13n'], q['r14n']], [q['r12'], 1, q['r23n'], q['r24n']],
                  [q['r13n'], q['r23n'], 1, q['r34n']], [q['r14n'], q['r24n'], q['r34n'], 1]])
    cov = np.outer(sig, sig) * R
    rng = np.random.default_rng(seed)
    xyzw = rng.multivariate_normal([q['mu_x1'], q['mu_x2'], q['mu_x3_n'], q['mu_x4_n']], cov, size=N)
    x1, x2 = xyzw[:, 0], xyzw[:, 1]
    x3, x4 = np.exp(xyzw[:, 2]), np.exp(xyzw[:, 3])
    mask = x2 < 0.0
    x2c = np.power(np.where(mask, x2, -1.0).astype(complex), q['alpha'])
    Y = x2c * x3 ** q['beta'] * x4 ** q['gamma']
    return np.where(mask, (x1 - q['x1_mean']) ** a * np.real((Y - q['M']) ** b), 0.0)


def test_quadrivar_general_complex_mc():
    """General quadrivar_NNLL_MM vs a truncated-domain Monte-Carlo over x2<0.

    For x2<0, x2^alpha is the real principal-branch continuation (= |x2|^alpha cos(pi alpha) ...), the same
    thing `_signed_pow` encodes; computing it as Re((x2+0j)^alpha ...) makes the MC integrand exact, and the MC
    itself uses NO D_v/Gamma — so agreement independently validates the b-expansion, signed power, and the D_v
    half-line integration.
    """
    q = _Q()
    N = 24_000_000
    rng = np.random.default_rng(9)
    s = [q['s_x1'], q['s_x2'], q['s_x3n'], q['s_x4n']]
    R = np.array([[1, q['r12'], q['r13n'], q['r14n']],
                  [q['r12'], 1, q['r23n'], q['r24n']],
                  [q['r13n'], q['r23n'], 1, q['r34n']],
                  [q['r14n'], q['r24n'], q['r34n'], 1]])
    cov = np.outer(s, s) * R
    xyzw = rng.multivariate_normal([q['mu_x1'], q['mu_x2'], q['mu_x3_n'], q['mu_x4_n']], cov, size=N)
    x1, x2 = xyzw[:, 0], xyzw[:, 1]
    x3, x4 = np.exp(xyzw[:, 2]), np.exp(xyzw[:, 3])
    mask = x2 < 0.0
    x2c = np.power(np.where(mask, x2, -1.0).astype(complex), q['alpha'])   # principal branch for x2<0
    Y = x2c * x3 ** q['beta'] * x4 ** q['gamma']
    for a, b in ((1, 1), (2, 1), (1, 2), (2, 2)):
        got = _call_quad(q, a, b)
        integ = np.where(mask, (x1 - q['x1_mean']) ** a * np.real((Y - q['M']) ** b), 0.0)
        scale = np.mean(np.abs(integ)) + 1e-30
        rel = abs(got - integ.mean()) / scale
        assert rel < 5e-3, f"quadrivar MC a={a},b={b}: rel {rel:.2e} (got {got}, mc {integ.mean()})"
    print("  quadrivar_NNLL_MM (4 (a,b)) vs truncated-domain 24M complex-branch MC: rel <5e-3  PASS")


def test_quadrivar_const_x1():
    q = _Q()
    for a, b in ((1, 1), (2, 1), (2, 2)):
        const = float(quadrivar_NNLL_MM_const_x1(
            q['mu_x1'], q['mu_x2'], q['mu_x3_n'], q['mu_x4_n'], q['s_x2'], q['s_x3n'], q['s_x4n'],
            q['r23n'], q['r24n'], q['r34n'], q['x1_mean'], q['M'], q['alpha'], q['beta'], q['gamma'], a, b))
        # limit: general quadrivar with sigma_x1 -> 0 (x1 correlations vanish with it)
        gen = _call_quad(q, a, b, s_x1=1e-8)
        assert abs(const - gen) / (abs(const) + 1e-30) < 1e-5, f"const_x1 limit a={a},b={b}"
        # MC: x1 fixed (sigma_x1=0)
        integ = _quad_complex_mc(q, a, b, 24_000_000, seed=11, s_over={'s_x1': 0.0})
        rel = abs(const - integ.mean()) / (np.mean(np.abs(integ)) + 1e-30)
        assert rel < 5e-3, f"const_x1 MC a={a},b={b}: rel {rel:.2e}"
    print("  quadrivar_NNLL_MM_const_x1: vs general(sigma_x1->0) <1e-5 + 24M complex MC <5e-3  PASS")


def test_quadrivar_const_x2():
    q = _Q()   # mu_x2 = -0.3 < 0 (subsaturated, where const_x2 is nonzero)
    for a, b in ((1, 1), (2, 1), (2, 2)):
        const = float(quadrivar_NNLL_MM_const_x2(
            q['mu_x1'], q['mu_x2'], q['mu_x3_n'], q['mu_x4_n'], q['s_x1'], q['s_x3n'], q['s_x4n'],
            q['r13n'], q['r14n'], q['r34n'], q['x1_mean'], q['M'], q['alpha'], q['beta'], q['gamma'], a, b))
        # NOTE: no sigma_x2->0 limit check — the general form is singular there (s_cc=mu_x2/sigma_x2 -> inf
        # hits the D_v argument clamp; that singularity is precisely why a separate const_x2 form exists).
        # The complex-branch MC (x2 fixed at mu_x2<0) is the independent oracle.
        integ = _quad_complex_mc(q, a, b, 24_000_000, seed=12, s_over={'s_x2': 0.0})
        rel = abs(const - integ.mean()) / (np.mean(np.abs(integ)) + 1e-30)
        assert rel < 5e-3, f"const_x2 MC a={a},b={b}: rel {rel:.2e}"
    # mu_x2 > 0 -> 0 (delta outside x2<0)
    pos = float(quadrivar_NNLL_MM_const_x2(0.5, 0.3, -0.4, -0.6, 0.6, 0.45, 0.5, 0.2, -0.1, 0.3,
                                           0.4, 1.3, 1.5, 0.5, 0.4, 2, 1))
    assert pos == 0.0, "const_x2 should be 0 for mu_x2>0"
    print("  quadrivar_NNLL_MM_const_x2: 24M complex MC (x2 fixed, mu_x2<0) <5e-3 + mu_x2>0 guard  PASS")


def test_quadrivar_const_x1x2():
    q = _Q()
    for a, b in ((1, 1), (2, 1), (3, 2)):
        const = float(quadrivar_NNLL_MM_const_x1x2(
            q['mu_x1'], q['mu_x2'], q['mu_x3_n'], q['mu_x4_n'], q['s_x3n'], q['s_x4n'], q['r34n'],
            q['x1_mean'], q['M'], q['alpha'], q['beta'], q['gamma'], a, b))
        # limit: const_x2 with sigma_x1 -> 0
        cx2 = float(quadrivar_NNLL_MM_const_x2(
            q['mu_x1'], q['mu_x2'], q['mu_x3_n'], q['mu_x4_n'], 1e-9, q['s_x3n'], q['s_x4n'],
            0.2, -0.1, q['r34n'], q['x1_mean'], q['M'], q['alpha'], q['beta'], q['gamma'], a, b))
        assert abs(const - cx2) / (abs(const) + 1e-30) < 1e-7, f"const_x1x2 limit a={a},b={b}"
        # MC: x1, x2 fixed
        integ = _quad_complex_mc(q, a, b, 8_000_000, seed=13, s_over={'s_x1': 0.0, 's_x2': 0.0})
        rel = abs(const - integ.mean()) / (np.mean(np.abs(integ)) + 1e-30)
        assert rel < 5e-3, f"const_x1x2 MC a={a},b={b}: rel {rel:.2e}"
    print("  quadrivar_NNLL_MM_const_x1x2: vs const_x2(sigma_x1->0) <1e-7 + 8M MC <5e-3  PASS")


def test_differentiable():
    g = jax.grad(lambda s: trivar_NNL_MM(0.5, 0.35, -0.4, s, 0.5, 0.45, 0.3, 0.2, -0.15,
                                         0.45, 1.2, 1.5, 0.5, 2, 1))(jnp.asarray(0.6))
    g1 = jax.grad(lambda s: trivar_NNL_MM_const_x1(0.5, 0.35, -0.4, s, 0.45, -0.15, 0.45, 1.2,
                                                   1.5, 0.5, 2, 1))(jnp.asarray(0.5))
    g2 = jax.grad(lambda s: trivar_NNL_MM_const_x2(0.5, 0.4, -0.4, s, 0.45, 0.25, 0.45, 1.1,
                                                   1.5, 0.5, 2, 1))(jnp.asarray(0.6))
    assert all(np.isfinite(float(x)) for x in (g, g1, g2)), "non-finite grad"
    print(f"  jax.grad: trivar/dσx1={float(g):+.3e}, const_x1/dσx2={float(g1):+.3e}, "
          f"const_x2/dσx1={float(g2):+.3e}: finite  PASS")


def main():
    print("test_pdf_integrals_all_mm:")
    for t in (test_base_case_probability_mass, test_general_monte_carlo,
              test_const_x1_limit_and_mc, test_const_x2_mc, test_const_x1x2_limit_and_mc,
              test_quadrivar_base_and_psum, test_quadrivar_general_complex_mc,
              test_quadrivar_const_x1, test_quadrivar_const_x2, test_quadrivar_const_x1x2,
              test_differentiable):
        t()
    print("All PDF_integrals_all_MM checks PASSED")


if __name__ == "__main__":
    main()

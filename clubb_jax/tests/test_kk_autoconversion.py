"""Verification of the upscaled-KK autoconversion analytic kernel (JAX port).

Two INDEPENDENT oracles, neither of which is the Fortran binary (the upscaled KK
microphysics is not exposed by the f2py CLUBB_core API), so each closed form is
checked against first principles:

  1. dv_parabolic_cylinder  vs  scipy.special.pbdv  (an independent evaluator of the
     parabolic cylinder function D_v) over the validated argument range.

  2. The four bivar_NL_mean* closed forms  vs  brute-force numerical quadrature of
     their DEFINING integral
        < x1^a x2^b >_i = INT_{x1>0} INT x1^a x2^b P_NL,i(x1,x2) dx1 dx2,
     with x1 ~ Normal, ln(x2) ~ Normal, correlation rho. The inner (independent)
     normal direction is integrated analytically, leaving a 1-D quadrature in u
     that uses NO special functions — a fully independent check of the algebra and
     of D_v together.

  3. KK_auto_upscaled_mean composes the pieces, is finite, and is differentiable
     (jax.grad) w.r.t. the PDF moments — the project's differentiability goal.

Run: PYTHONPATH=...:. python clubb_jax/tests/test_kk_autoconversion.py
"""
import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from scipy.special import pbdv  # independent special-function oracle

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.Microphys.KK_microphys.parabolic_cylinder import (
    dv_parabolic_cylinder,
)
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    bivar_NL_mean,
    bivar_NL_mean_const_x1,
    bivar_NL_mean_const_x2,
    bivar_NL_mean_const_all,
    trivar_NLL_mean,
    bivar_LL_mean,
)
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    KK_auto_upscaled_mean,
    KK_accr_upscaled_mean,
    KK_evap_upscaled_mean,
    KK_mvr_upscaled_mean,
    KK_AUTO_RC_EXP,
    KK_AUTO_NC_EXP,
    KK_ACCR_RC_EXP,
    KK_ACCR_RR_EXP,
    KK_EVAP_SUPERSAT_EXP,
    KK_EVAP_RR_EXP,
    KK_EVAP_NR_EXP,
    KK_MVR_RR_EXP,
    KK_MVR_NR_EXP,
)
from clubb_jax.src.Microphys.KK_microphys_module import kk_auto_coef

ALPHA = KK_AUTO_RC_EXP    # 2.47
BETA = KK_AUTO_NC_EXP     # -1.79
ACCR_A = KK_ACCR_RC_EXP   # 1.15
ACCR_B = KK_ACCR_RR_EXP   # 1.15
EVAP_A = KK_EVAP_SUPERSAT_EXP  # 1
EVAP_B = KK_EVAP_RR_EXP        # 1/3
EVAP_G = KK_EVAP_NR_EXP        # 2/3
MVR_A = KK_MVR_RR_EXP     # 1/3
MVR_B = KK_MVR_NR_EXP     # -1/3


# ---------------------------------------------------------------------------
def test_dv_vs_scipy():
    """D_v matches scipy.special.pbdv on the validated zones.

    Series zone z in [-4, 4] (rel < 2e-8, ~1e-15 near 0, ~1e-8 at the z=4 edge) and
    optimally-truncated asymptotic zone
    z in [7, 48] (rel < 1e-6). The narrow z ~ [5, 6.5] handoff band carries a
    worst-case rel ~ 2e-4 for the steepest KK exponent (v=-3.47) and is checked at a
    looser tolerance — it is physically suppressed by exp(-s_c^2/4) in the integral
    that consumes D_v (see module docstring)."""
    worst_series = worst_asym = worst_band = worst_negasym = 0.0
    for v in (-3.47, -1.79, -2.0):
        zs = np.linspace(-4.0, 4.0, 49)
        mine = np.asarray(dv_parabolic_cylinder(v, jnp.asarray(zs)))
        ref = np.array([pbdv(v, z)[0] for z in zs])
        worst_series = max(worst_series, (np.abs(mine - ref) / (np.abs(ref) + 1e-300)).max())
        # Cap at z=30 (value ~1e-100): beyond that scipy.pbdv loses relative precision
        # as it approaches its own underflow boundary (our optimally-truncated
        # asymptotic stays accurate further, but scipy is no longer a clean oracle).
        za = np.linspace(7.0, 30.0, 40)
        mine_a = np.asarray(dv_parabolic_cylinder(v, jnp.asarray(za)))
        ref_a = np.array([pbdv(v, z)[0] for z in za])
        worst_asym = max(worst_asym, (np.abs(mine_a - ref_a) / np.abs(ref_a)).max())
        zb = np.linspace(5.0, 6.5, 16)
        mine_b = np.asarray(dv_parabolic_cylinder(v, jnp.asarray(zb)))
        ref_b = np.array([pbdv(v, z)[0] for z in zb])
        worst_band = max(worst_band, (np.abs(mine_b - ref_b) / (np.abs(ref_b) + 1e-300)).max())
        # Large-NEGATIVE-z growing branch (Iter130) — needed for narrow chi PDFs (e.g.
        # stratocumulus, s_c ~ 32 -> D_v arg ~ -32). pbdv stays finite to z ~ -45 here.
        zn = np.linspace(-45.0, -8.0, 40)
        mine_n = np.asarray(dv_parabolic_cylinder(v, jnp.asarray(zn)))
        ref_n = np.array([pbdv(v, z)[0] for z in zn])
        worst_negasym = max(worst_negasym, (np.abs(mine_n - ref_n) / np.abs(ref_n)).max())
    assert worst_series < 2e-8, f"D_v series vs pbdv worst rel {worst_series:.2e}"
    assert worst_asym < 1e-6, f"D_v asym vs pbdv worst rel {worst_asym:.2e}"
    assert worst_band < 1e-3, f"D_v handoff-band worst rel {worst_band:.2e}"
    assert worst_negasym < 1e-6, f"D_v neg-asym vs pbdv worst rel {worst_negasym:.2e}"
    print(f"  D_v vs scipy.pbdv: series<{worst_series:.1e}, +asym<{worst_asym:.1e}, "
          f"-asym<{worst_negasym:.1e}, handoff-band<{worst_band:.1e}  PASS")


# ---------------------------------------------------------------------------
def _quad_bivar_NL(mu_x1, mu_x2_n, sigma_x1, sigma_x2_n, rho, a, b):
    """Brute-force < x1^a x2^b > by 1-D quadrature of the defining NL integral.

    x1 = mu_x1 + sigma_x1 u, ln x2 = mu_x2_n + sigma_x2_n (rho u + sqrt(1-rho^2) v),
    u,v ~ N(0,1) indep. Integrate v analytically (E[exp(b sigma2n sqrt(1-rho^2) v)]
    = exp(1/2 b^2 sigma2n^2 (1-rho^2))); quadrature over u with x1>0 only.
    """
    u = np.linspace(-12.0, 12.0, 400001)
    phi = np.exp(-0.5 * u * u) / np.sqrt(2.0 * np.pi)
    x1 = mu_x1 + sigma_x1 * u
    integrand = np.where(x1 > 0.0, np.power(np.maximum(x1, 0.0), a)
                         * np.exp(b * sigma_x2_n * rho * u) * phi, 0.0)
    inner = np.trapezoid(integrand, u)
    return np.exp(b * mu_x2_n + 0.5 * b * b * sigma_x2_n ** 2 * (1.0 - rho * rho)) * inner


def test_bivar_general_vs_quadrature():
    """bivar_NL_mean (both vary) matches the 1-D quadrature of its defining integral."""
    # Physical-ish autoconversion moments: chi ~ small positive, Ncn ~ 1e8 num/kg.
    cases = [
        dict(mu_x1=1.0e-4, sigma_x1=4.0e-4, mu_x2_n=np.log(1.0e8),
             sigma_x2_n=0.4, rho=-0.5),
        dict(mu_x1=5.0e-4, sigma_x1=2.0e-4, mu_x2_n=np.log(5.0e7),
             sigma_x2_n=0.6, rho=-0.3),
        dict(mu_x1=-2.0e-4, sigma_x1=6.0e-4, mu_x2_n=np.log(2.0e8),
             sigma_x2_n=0.5, rho=0.2),
    ]
    worst = 0.0
    for c in cases:
        closed = float(bivar_NL_mean(c["mu_x1"], c["mu_x2_n"], c["sigma_x1"],
                                     c["sigma_x2_n"], c["rho"], ALPHA, BETA))
        quad = _quad_bivar_NL(c["mu_x1"], c["mu_x2_n"], c["sigma_x1"],
                              c["sigma_x2_n"], c["rho"], ALPHA, BETA)
        rel = abs(closed - quad) / (abs(quad) + 1e-300)
        worst = max(worst, rel)
        assert rel < 1e-5, f"bivar_NL_mean vs quad rel {rel:.2e} for {c}"
    print(f"  bivar_NL_mean (general) vs 2-D quadrature: worst rel {worst:.1e}  PASS")


def test_bivar_const_variants_vs_quadrature():
    """const_x1 / const_x2 / const_all match the limiting quadratures."""
    mu_x1, sigma_x1 = 3.0e-4, 5.0e-4
    mu_x2, mu_x2_n, sigma_x2_n = 1.0e8, np.log(1.0e8), 0.5
    # const_x1: x1 constant -> < mu_x1^a x2^b > = mu_x1^a exp(b mu2n + 1/2 b^2 s2n^2)
    f_x1 = float(bivar_NL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, ALPHA, BETA))
    ref_x1 = mu_x1 ** ALPHA * np.exp(BETA * mu_x2_n + 0.5 * BETA ** 2 * sigma_x2_n ** 2)
    assert abs(f_x1 - ref_x1) / abs(ref_x1) < 1e-12
    # const_x2: x2 constant -> mu_x2^b * INT_{x1>0} x1^a phi du
    f_x2 = float(bivar_NL_mean_const_x2(mu_x1, mu_x2, sigma_x1, ALPHA, BETA))
    u = np.linspace(-12.0, 12.0, 400001)
    x1 = mu_x1 + sigma_x1 * u
    inner = np.trapezoid(np.where(x1 > 0, np.power(np.maximum(x1, 0.0), ALPHA), 0.0)
                     * np.exp(-0.5 * u * u) / np.sqrt(2 * np.pi), u)
    ref_x2 = mu_x2 ** BETA * inner
    rel_x2 = abs(f_x2 - ref_x2) / abs(ref_x2)
    assert rel_x2 < 1e-5, f"const_x2 vs quad rel {rel_x2:.2e}"
    # const_all
    f_all = float(bivar_NL_mean_const_all(mu_x1, mu_x2, ALPHA, BETA))
    assert abs(f_all - mu_x1 ** ALPHA * mu_x2 ** BETA) / (mu_x1 ** ALPHA * mu_x2 ** BETA) < 1e-12
    # const_x1/all return 0 for negative mu_x1
    assert float(bivar_NL_mean_const_x1(-1.0, mu_x2_n, sigma_x2_n, ALPHA, BETA)) == 0.0
    assert float(bivar_NL_mean_const_all(-1.0, mu_x2, ALPHA, BETA)) == 0.0
    print(f"  const_x1/x2/all vs quadrature: PASS (const_x2 rel {rel_x2:.1e})")


# ---------------------------------------------------------------------------
def test_kk_auto_mean_composes_and_differentiable():
    """KK_auto_upscaled_mean is finite, positive, and differentiable in the moments."""
    rho_air = 1.0  # kg/m^3
    coef = kk_auto_coef(rho_air)
    args = dict(
        mu_chi_1=2.0e-4, mu_chi_2=1.0e-4, mu_Ncn_1=1.0e8, mu_Ncn_2=1.0e8,
        mu_Ncn_1_n=np.log(1.0e8), mu_Ncn_2_n=np.log(1.0e8),
        sigma_chi_1=3.0e-4, sigma_chi_2=2.0e-4,
        sigma_Ncn_1=3.0e7, sigma_Ncn_2=3.0e7,
        sigma_Ncn_1_n=0.4, sigma_Ncn_2_n=0.4,
        corr_chi_Ncn_1_n=-0.5, corr_chi_Ncn_2_n=-0.5,
        KK_auto_coef_val=coef, mixt_frac=0.3,
    )
    val = float(KK_auto_upscaled_mean(**args))
    assert np.isfinite(val) and val > 0.0, f"KK_auto mean not finite/positive: {val}"

    # differentiable w.r.t. mu_chi_1 (the cloud-water mean) — gradient finite & nonzero
    def f(mu_chi_1):
        a = dict(args)
        a["mu_chi_1"] = mu_chi_1
        return KK_auto_upscaled_mean(**a)
    g = float(jax.grad(f)(2.0e-4))
    # finite-difference check
    eps = 1.0e-8
    fd = (f(2.0e-4 + eps) - f(2.0e-4 - eps)) / (2 * eps)
    rel = abs(g - fd) / (abs(fd) + 1e-300)
    assert np.isfinite(g) and abs(g) > 0.0, "KK_auto grad not finite/nonzero"
    assert rel < 1e-4, f"KK_auto grad wrong: ad={g:.4e} fd={fd:.4e} rel={rel:.2e}"
    print(f"  KK_auto_upscaled_mean: val={val:.3e}/s, grad correct (rel {rel:.1e})  PASS")


def test_accretion_bivar_vs_quadrature():
    """bivar_NL_mean with the ACCRETION exponents (alpha=beta=1.15) matches quadrature.

    Accretion uses the same NL bivariate as autoconversion but with y=r_r and positive
    exponents; r_r ~ 1e-4 kg/kg with broad lognormal spread."""
    cases = [
        dict(mu_x1=2.0e-4, sigma_x1=3.0e-4, mu_x2_n=np.log(1.0e-4),
             sigma_x2_n=0.8, rho=0.4),
        dict(mu_x1=1.0e-4, sigma_x1=5.0e-4, mu_x2_n=np.log(5.0e-5),
             sigma_x2_n=1.0, rho=-0.2),
    ]
    worst = 0.0
    for c in cases:
        closed = float(bivar_NL_mean(c["mu_x1"], c["mu_x2_n"], c["sigma_x1"],
                                     c["sigma_x2_n"], c["rho"], ACCR_A, ACCR_B))
        quad = _quad_bivar_NL(c["mu_x1"], c["mu_x2_n"], c["sigma_x1"],
                              c["sigma_x2_n"], c["rho"], ACCR_A, ACCR_B)
        rel = abs(closed - quad) / (abs(quad) + 1e-300)
        worst = max(worst, rel)
        assert rel < 1e-5, f"accretion bivar vs quad rel {rel:.2e} for {c}"
    print(f"  bivar_NL_mean (accretion exps) vs quadrature: worst rel {worst:.1e}  PASS")


def test_kk_accr_mean_composes_and_differentiable():
    """KK_accr_upscaled_mean is finite, positive, and differentiable in the moments."""
    args = dict(
        mu_chi_1=2.0e-4, mu_chi_2=1.0e-4, mu_rr_1=1.0e-4, mu_rr_2=8.0e-5,
        mu_rr_1_n=np.log(1.0e-4), mu_rr_2_n=np.log(8.0e-5),
        sigma_chi_1=3.0e-4, sigma_chi_2=2.0e-4,
        sigma_rr_1=5.0e-5, sigma_rr_2=4.0e-5,
        sigma_rr_1_n=0.7, sigma_rr_2_n=0.7,
        corr_chi_rr_1_n=0.3, corr_chi_rr_2_n=0.3,
        mixt_frac=0.3, precip_frac_1=0.8, precip_frac_2=0.5,
    )
    val = float(KK_accr_upscaled_mean(**args))
    assert np.isfinite(val) and val > 0.0, f"KK_accr mean not finite/positive: {val}"

    def f(mu_chi_1):
        a = dict(args)
        a["mu_chi_1"] = mu_chi_1
        return KK_accr_upscaled_mean(**a)
    g = float(jax.grad(f)(2.0e-4))
    eps = 1.0e-8
    fd = (f(2.0e-4 + eps) - f(2.0e-4 - eps)) / (2 * eps)
    rel = abs(g - fd) / (abs(fd) + 1e-300)
    assert np.isfinite(g) and abs(g) > 0.0, "KK_accr grad not finite/nonzero"
    assert rel < 1e-4, f"KK_accr grad wrong: ad={g:.4e} fd={fd:.4e} rel={rel:.2e}"
    print(f"  KK_accr_upscaled_mean: val={val:.3e}/s, grad correct (rel {rel:.1e})  PASS")


def _quad_trivar_NLL(mu_x1, mu_x2_n, mu_x3_n, sigma_x1, sigma_x2_n, sigma_x3_n,
                     rho12, rho13, rho23, a, b, g):
    """Brute-force < x1^a x2^b x3^g > over x1<0 by 1-D quadrature.

    x1 = mu1 + sig1 z; ln x2, ln x3 lognormal correlated with z. The two lognormal
    directions integrate analytically conditional on z, leaving a 1-D integral in z
    over the SUBSATURATED region x1<0. V and c below are the conditional variance/mean
    coefficients (independent derivation of the closed form)."""
    V = (b**2 * sigma_x2_n**2 * (1 - rho12**2)
         + g**2 * sigma_x3_n**2 * (1 - rho13**2)
         + 2 * b * g * sigma_x2_n * sigma_x3_n * (rho23 - rho12 * rho13))
    c = b * sigma_x2_n * rho12 + g * sigma_x3_n * rho13
    z = np.linspace(-12.0, 12.0, 800001)
    x1 = mu_x1 + sigma_x1 * z
    chi_pow = np.sign(x1) * np.abs(x1) ** a          # signed power (a=1 -> x1)
    integrand = np.where(x1 < 0.0,
                         chi_pow * np.exp(c * z) * np.exp(-0.5 * z * z) / np.sqrt(2 * np.pi),
                         0.0)
    inner = np.trapezoid(integrand, z)
    return np.exp(b * mu_x2_n + g * mu_x3_n + 0.5 * V) * inner


def test_trivar_general_vs_quadrature():
    """trivar_NLL_mean (evaporation, all vary) matches the 1-D quadrature over chi<0."""
    cases = [
        dict(mu_x1=-3.0e-4, sigma_x1=4.0e-4, mu_x2_n=np.log(1.0e-4), sigma_x2_n=0.7,
             mu_x3_n=np.log(1.0e4), sigma_x3_n=0.5, rho12=0.3, rho13=0.2, rho23=0.5),
        dict(mu_x1=-1.0e-4, sigma_x1=6.0e-4, mu_x2_n=np.log(5.0e-5), sigma_x2_n=0.9,
             mu_x3_n=np.log(5.0e3), sigma_x3_n=0.6, rho12=-0.2, rho13=0.4, rho23=0.1),
    ]
    worst = 0.0
    for c in cases:
        closed = float(trivar_NLL_mean(c["mu_x1"], c["mu_x2_n"], c["mu_x3_n"],
                                       c["sigma_x1"], c["sigma_x2_n"], c["sigma_x3_n"],
                                       c["rho12"], c["rho13"], c["rho23"],
                                       EVAP_A, EVAP_B, EVAP_G))
        quad = _quad_trivar_NLL(c["mu_x1"], c["mu_x2_n"], c["mu_x3_n"],
                                c["sigma_x1"], c["sigma_x2_n"], c["sigma_x3_n"],
                                c["rho12"], c["rho13"], c["rho23"], EVAP_A, EVAP_B, EVAP_G)
        rel = abs(closed - quad) / (abs(quad) + 1e-300)
        worst = max(worst, rel)
        assert rel < 1e-5, f"trivar_NLL_mean vs quad rel {rel:.2e} for {c}"
    print(f"  trivar_NLL_mean (evaporation) vs 3-D quadrature: worst rel {worst:.1e}  PASS")


def test_kk_evap_mean_composes_and_differentiable():
    """KK_evap_upscaled_mean is finite, NEGATIVE (removes rain), and differentiable."""
    args = dict(
        mu_chi_1=-2.0e-4, mu_chi_2=-1.0e-4, mu_rr_1=1.0e-4, mu_rr_2=8.0e-5,
        mu_Nr_1=1.0e4, mu_Nr_2=8.0e3,
        mu_rr_1_n=np.log(1.0e-4), mu_rr_2_n=np.log(8.0e-5),
        mu_Nr_1_n=np.log(1.0e4), mu_Nr_2_n=np.log(8.0e3),
        sigma_chi_1=3.0e-4, sigma_chi_2=2.0e-4,
        sigma_rr_1=5.0e-5, sigma_rr_2=4.0e-5, sigma_Nr_1=5.0e3, sigma_Nr_2=4.0e3,
        sigma_rr_1_n=0.6, sigma_rr_2_n=0.6, sigma_Nr_1_n=0.5, sigma_Nr_2_n=0.5,
        corr_chi_rr_1_n=0.3, corr_chi_rr_2_n=0.3,
        corr_chi_Nr_1_n=0.2, corr_chi_Nr_2_n=0.2,
        corr_rr_Nr_1_n=0.5, corr_rr_Nr_2_n=0.5,
        KK_evap_coef_val=1.0e-3, mixt_frac=0.3, precip_frac_1=0.8, precip_frac_2=0.5,
    )
    val = float(KK_evap_upscaled_mean(**args))
    assert np.isfinite(val) and val < 0.0, f"KK_evap mean not finite/negative: {val}"

    def f(mu_chi_1):
        a = dict(args)
        a["mu_chi_1"] = mu_chi_1
        return KK_evap_upscaled_mean(**a)
    g = float(jax.grad(f)(-2.0e-4))
    eps = 1.0e-8
    fd = (f(-2.0e-4 + eps) - f(-2.0e-4 - eps)) / (2 * eps)
    rel = abs(g - fd) / (abs(fd) + 1e-300)
    assert np.isfinite(g) and abs(g) > 0.0, "KK_evap grad not finite/nonzero"
    assert rel < 1e-4, f"KK_evap grad wrong: ad={g:.4e} fd={fd:.4e} rel={rel:.2e}"
    print(f"  KK_evap_upscaled_mean: val={val:.3e}/s (<0), grad correct (rel {rel:.1e})  PASS")


def _quad_bivar_LL(mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho, a, b):
    """Brute-force < x1^a x2^b > for two lognormals by 1-D quadrature.

    ln x1 = mu1n + s1n z1, ln x2 = mu2n + s2n(rho z1 + sqrt(1-rho^2) z2). Integrate z2
    analytically; numerical 1-D integral over z1 (independent of the closed-form algebra)."""
    z = np.linspace(-12.0, 12.0, 400001)
    inner = np.trapezoid(np.exp((a * sigma_x1_n + b * sigma_x2_n * rho) * z)
                         * np.exp(-0.5 * z * z) / np.sqrt(2 * np.pi), z)
    return np.exp(a * mu_x1_n + b * mu_x2_n
                  + 0.5 * b * b * sigma_x2_n ** 2 * (1.0 - rho * rho)) * inner


def test_bivar_LL_vs_quadrature():
    """bivar_LL_mean (mvr, both lognormal) matches the 1-D quadrature of its integral."""
    cases = [
        dict(mu_x1_n=np.log(1.0e-4), mu_x2_n=np.log(1.0e4),
             sigma_x1_n=0.7, sigma_x2_n=0.5, rho=0.5),
        dict(mu_x1_n=np.log(5.0e-5), mu_x2_n=np.log(5.0e3),
             sigma_x1_n=1.0, sigma_x2_n=0.6, rho=-0.3),
    ]
    worst = 0.0
    for c in cases:
        closed = float(bivar_LL_mean(c["mu_x1_n"], c["mu_x2_n"], c["sigma_x1_n"],
                                     c["sigma_x2_n"], c["rho"], MVR_A, MVR_B))
        quad = _quad_bivar_LL(c["mu_x1_n"], c["mu_x2_n"], c["sigma_x1_n"],
                              c["sigma_x2_n"], c["rho"], MVR_A, MVR_B)
        rel = abs(closed - quad) / (abs(quad) + 1e-300)
        worst = max(worst, rel)
        assert rel < 1e-5, f"bivar_LL_mean vs quad rel {rel:.2e} for {c}"
    print(f"  bivar_LL_mean (mvr) vs quadrature: worst rel {worst:.1e}  PASS")


def test_kk_mvr_mean_composes_and_differentiable():
    """KK_mvr_upscaled_mean is finite, positive (a radius), and differentiable."""
    args = dict(
        mu_rr_1=1.0e-4, mu_rr_2=8.0e-5, mu_Nr_1=1.0e4, mu_Nr_2=8.0e3,
        mu_rr_1_n=np.log(1.0e-4), mu_rr_2_n=np.log(8.0e-5),
        mu_Nr_1_n=np.log(1.0e4), mu_Nr_2_n=np.log(8.0e3),
        sigma_rr_1=5.0e-5, sigma_rr_2=4.0e-5, sigma_Nr_1=5.0e3, sigma_Nr_2=4.0e3,
        sigma_rr_1_n=0.6, sigma_rr_2_n=0.6, sigma_Nr_1_n=0.5, sigma_Nr_2_n=0.5,
        corr_rr_Nr_1_n=0.5, corr_rr_Nr_2_n=0.5,
        mixt_frac=0.3, precip_frac_1=0.8, precip_frac_2=0.5,
    )
    val = float(KK_mvr_upscaled_mean(**args))
    # a physical rain-drop radius: ~tens of microns to ~mm
    assert np.isfinite(val) and 1e-6 < val < 1e-2, f"KK_mvr not a plausible radius: {val}"

    # Differentiate w.r.t. a LOG moment (mu_rr_1_n): the selected general bivar_LL form
    # depends on the log moments, not the linear mean (grad w.r.t. mu_rr_1 is correctly 0).
    x0 = np.log(1.0e-4)
    def f(mu_rr_1_n):
        a = dict(args)
        a["mu_rr_1_n"] = mu_rr_1_n
        return KK_mvr_upscaled_mean(**a)
    g = float(jax.grad(f)(x0))
    eps = 1.0e-6
    fd = (f(x0 + eps) - f(x0 - eps)) / (2 * eps)
    rel = abs(g - fd) / (abs(fd) + 1e-300)
    assert np.isfinite(g) and abs(g) > 0.0, "KK_mvr grad not finite/nonzero"
    assert rel < 1e-4, f"KK_mvr grad wrong: ad={g:.4e} fd={fd:.4e} rel={rel:.2e}"
    print(f"  KK_mvr_upscaled_mean: val={val:.3e} m, grad correct (rel {rel:.1e})  PASS")


if __name__ == "__main__":
    print("Upscaled-KK auto + accr + evap + mvr kernel verification:")
    test_dv_vs_scipy()
    test_bivar_general_vs_quadrature()
    test_bivar_const_variants_vs_quadrature()
    test_kk_auto_mean_composes_and_differentiable()
    test_accretion_bivar_vs_quadrature()
    test_kk_accr_mean_composes_and_differentiable()
    test_trivar_general_vs_quadrature()
    test_kk_evap_mean_composes_and_differentiable()
    test_bivar_LL_vs_quadrature()
    test_kk_mvr_mean_composes_and_differentiable()
    print("All KK auto + accr + evap + mvr tests PASSED.")

"""JAX port of KK_microphys/PDF_integrals_all_MM.F90 — the all-mixed-moment NN(LL) PDF integrals.

These evaluate the higher-order mixed moments the upscaled-KK SECOND-moment microphysics needs, of the form

    < (x1 - <x1>)^a  (x2^alpha x3^beta - <x2^alpha x3^beta>)^b >_i

over the ith PDF component, where x1 and x2 are NORMAL marginals (x2 integrated over x2 > 0, hence the
parabolic-cylinder D_v) and x3 is LOGNORMAL (ln x3 normal). `quadrivar_NNLL_MM` adds a second lognormal x4.

Reference: Griffin (2016), PhD thesis, Eq. (I22). The transcendental pieces — the parabolic-cylinder function
D_v and the Gamma function — are reused from `parabolic_cylinder` (the same accurate float64 implementations the
already-ported means/covariance integrals use). Integer orders a_exp/b_exp are static → the finite p/r/q sums
unroll at trace time; alpha_exp/beta_exp are float exponents. Pure jnp → differentiable.

FULL port (8/8): the trivariate family (`trivar_NNL_MM` + const-x1/x2/x1x2, x2>0) and the quadrivariate family
(`quadrivar_NNLL_MM` + const-x1/x2/x1x2, x2<0 subsaturated, x3,x4 lognormal). Validated in
`tests/test_pdf_integrals_all_mm.py`: the a=b=0 base cases reduce analytically to the x2>0/x2<0 probability mass
Phi(±mu_x2/sigma_x2), a=2,b=0 reduces to a normal central moment × Phi, and every general/const form matches a
truncated-domain Monte-Carlo (the quadrivar forms use a complex principal-branch power for the x2<0 region,
matching `_signed_pow`; the const-x2/x1x2 forms are guarded to mu_x2<=0, where the x2 delta lies in x2<0).
"""
import math

import jax.numpy as jnp

from clubb_jax.src.Microphys.KK_microphys.parabolic_cylinder import (
    _gamma_real,
    _dvc,          # clamped parabolic-cylinder D_v (its home, next to dv_parabolic_cylinder)
)

_INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)


def trivar_NNL_MM(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x2, sigma_x3_n,
                  rho_x1x2, rho_x1x3_n, rho_x2x3_n, x1_mean, x2_alpha_x3_beta_mean,
                  alpha_exp, beta_exp, a_exp, b_exp):
    """< (x1-<x1>)^a (x2^alpha x3^beta - <x2^alpha x3^beta>)^b >_i, x1,x2 normal (x2>0), x3 lognormal.

    PDF_integrals_all_MM.F90:trivar_NNL_MM (Griffin 2016, Eq. I22). Triple sum over p (x1 variance pairs),
    r (x1 cross terms), q (the binomial expansion of the central x2^alpha x3^beta power):

      = SUM_p SUM_r SUM_q  (1/sqrt(2 pi)) [a!/((a-2p)!p!)] [(a-2p)!/((a-2p-r)!r!)] [b!/((b-q)!q!)]
          (-<x2^a x3^b>)^(b-q) sigma_x2^(alpha q) (½(1-rho12²)sigma_x1²)^p (rho12 sigma_x1)^r
          (mu_x1 - <x1> - (mu_x2/sigma_x2) rho12 sigma_x1
             + (rho13_n - rho12 rho23_n) sigma_x1 sigma_x3_n beta q)^(a-2p-r)
          exp(mu_x3_n beta q + ½ sigma_x3_n² beta² q² - ¼ s_c²)
          Gamma(alpha q + r + 1) D_{-(alpha q + r + 1)}(-s_c),
      s_c = mu_x2/sigma_x2 + rho23_n sigma_x3_n beta q.
    """
    total = 0.0
    for p in range(a_exp // 2 + 1):
        c_p = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        rem = a_exp - 2 * p
        for r in range(rem + 1):
            c_r = math.factorial(rem) // (math.factorial(rem - r) * math.factorial(r))
            for q in range(b_exp + 1):
                c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
                s_c = (mu_x2 / sigma_x2) + rho_x2x3_n * sigma_x3_n * beta_exp * q
                shifted = (mu_x1 - x1_mean
                           - (mu_x2 / sigma_x2) * rho_x1x2 * sigma_x1
                           + (rho_x1x3_n - rho_x1x2 * rho_x2x3_n) * sigma_x1 * sigma_x3_n * beta_exp * q)
                order = alpha_exp * q + r + 1.0
                total = total + _INV_SQRT_2PI * (c_p * c_r * c_q) \
                    * (-x2_alpha_x3_beta_mean) ** (b_exp - q) \
                    * sigma_x2 ** (alpha_exp * q) \
                    * (0.5 * (1.0 - rho_x1x2 ** 2) * sigma_x1 ** 2) ** p \
                    * (rho_x1x2 * sigma_x1) ** r \
                    * shifted ** (a_exp - 2 * p - r) \
                    * jnp.exp(mu_x3_n * beta_exp * q + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2
                              - 0.25 * s_c ** 2) \
                    * _gamma_real(order) \
                    * _dvc(-order, -s_c)
    return total


def trivar_NNL_MM_const_x1(mu_x1, mu_x2, mu_x3_n, sigma_x2, sigma_x3_n, rho_x2x3_n,
                           x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, a_exp, b_exp):
    """trivar_NNL_MM with x1 constant (sigma_x1 -> 0). PDF_integrals_all_MM.F90:trivar_NNL_MM_const_x1 (Eq. I23).

    x1 collapses to (mu_x1 - <x1>)^a; the x2>0 half-line integral (D_v) and the x3 lognormal raw moment remain:
      = SUM_q (1/sqrt 2pi) (mu_x1-<x1>)^a [b!/((b-q)!q!)] (-M)^(b-q) sigma_x2^(alpha q)
          exp(mu_x3_n beta q + ½ sigma_x3_n² beta² q² - ¼ s_c²) Gamma(alpha q + 1) D_{-(alpha q + 1)}(-s_c),
      s_c = mu_x2/sigma_x2 + rho_x2x3_n sigma_x3_n beta q.
    """
    total = 0.0
    for q in range(b_exp + 1):
        c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
        s_c = (mu_x2 / sigma_x2) + rho_x2x3_n * sigma_x3_n * beta_exp * q
        order = alpha_exp * q + 1.0
        total = total + _INV_SQRT_2PI * (mu_x1 - x1_mean) ** a_exp * c_q \
            * (-x2_alpha_x3_beta_mean) ** (b_exp - q) \
            * sigma_x2 ** (alpha_exp * q) \
            * jnp.exp(mu_x3_n * beta_exp * q + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2 - 0.25 * s_c ** 2) \
            * _gamma_real(order) * _dvc(-order, -s_c)
    return total


def trivar_NNL_MM_const_x2(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x3_n, rho_x1x3_n,
                           x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, a_exp, b_exp):
    """trivar_NNL_MM with x2 constant (sigma_x2 -> 0). PDF_integrals_all_MM.F90:trivar_NNL_MM_const_x2 (Eq. I24).

    x2 collapses to mu_x2^alpha (no D_v / half-line integral); x1 stays a FULL normal (central moments via the
    p-sum, mean shifted by the x1-lnx3 correlation), x3 lognormal:
      = SUM_p SUM_q [a!/((a-2p)!p!)] [b!/((b-q)!q!)] (-M)^(b-q) mu_x2^(alpha q) (½ sigma_x1²)^p
          (mu_x1 - <x1> + rho_x1x3_n sigma_x1 sigma_x3_n beta q)^(a-2p) exp(mu_x3_n beta q + ½ sigma_x3_n² beta² q²).
    """
    total = 0.0
    for p in range(a_exp // 2 + 1):
        c_p = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        for q in range(b_exp + 1):
            c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
            shifted = mu_x1 - x1_mean + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp * q
            total = total + (c_p * c_q) \
                * (-x2_alpha_x3_beta_mean) ** (b_exp - q) \
                * mu_x2 ** (alpha_exp * q) \
                * (0.5 * sigma_x1 ** 2) ** p \
                * shifted ** (a_exp - 2 * p) \
                * jnp.exp(mu_x3_n * beta_exp * q + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2)
    return total


def trivar_NNL_MM_const_x1x2(mu_x1, mu_x2, mu_x3_n, sigma_x3_n, x1_mean, x2_alpha_x3_beta_mean,
                             alpha_exp, beta_exp, a_exp, b_exp):
    """trivar_NNL_MM with x1 AND x2 constant. PDF_integrals_all_MM.F90:trivar_NNL_MM_const_x1x2 (Eq. I26).

    Only x3 varies (lognormal):
      = SUM_q (mu_x1-<x1>)^a [b!/((b-q)!q!)] (-M)^(b-q) mu_x2^(alpha q) exp(mu_x3_n beta q + ½ sigma_x3_n² beta² q²).
    """
    total = 0.0
    for q in range(b_exp + 1):
        c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
        total = total + (mu_x1 - x1_mean) ** a_exp * c_q \
            * (-x2_alpha_x3_beta_mean) ** (b_exp - q) \
            * mu_x2 ** (alpha_exp * q) \
            * jnp.exp(mu_x3_n * beta_exp * q + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2)
    return total


# ── Quadrivariate family (x1,x2 normal with x2<0 subsaturated; x3,x4 lognormal) ──
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import _signed_pow


def quadrivar_NNLL_MM(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x1, sigma_x2, sigma_x3_n, sigma_x4_n,
                      rho_x1x2, rho_x1x3_n, rho_x1x4_n, rho_x2x3_n, rho_x2x4_n, rho_x3x4_n,
                      x1_mean, x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, a_exp, b_exp):
    """< (x1-<x1>)^a (x2^α x3^β x4^γ - M)^b >_i, x1,x2 normal (x2 over the SUBSATURATED x2<0 half-line),
    x3,x4 lognormal. PDF_integrals_all_MM.F90:quadrivar_NNLL_MM (Griffin 2016, Eq. I1).

    The x2<0 region gives `(-σ_x2)^(αq)` (the signed continuation `_signed_pow`), the parabolic-cylinder
    argument +s_cc, and the prefactor exp(¼s_cc² - (μ_x2/σ_x2)s_cc + ½μ_x2²/σ_x2²) — the same convention as the
    validated `quadrivar_NNLL_covar`. Identity: quadrivar_NNLL_MM(a=1,b=1) == that covariance.

      s_cc = μ_x2/σ_x2 + ρ_x2x3_n σ_x3_n β q + ρ_x2x4_n σ_x4_n γ q.
    """
    total = 0.0
    for p in range(a_exp // 2 + 1):
        c_p = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        rem = a_exp - 2 * p
        for r in range(rem + 1):
            c_r = math.factorial(rem) // (math.factorial(rem - r) * math.factorial(r))
            for q in range(b_exp + 1):
                c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
                s_cc = (mu_x2 / sigma_x2 + rho_x2x3_n * sigma_x3_n * beta_exp * q
                        + rho_x2x4_n * sigma_x4_n * gamma_exp * q)
                shifted = (mu_x1 - x1_mean
                           - (mu_x2 / sigma_x2) * rho_x1x2 * sigma_x1
                           + (rho_x1x3_n - rho_x1x2 * rho_x2x3_n) * sigma_x1 * sigma_x3_n * beta_exp * q
                           + (rho_x1x4_n - rho_x1x2 * rho_x2x4_n) * sigma_x1 * sigma_x4_n * gamma_exp * q)
                order = alpha_exp * q + r + 1.0
                total = total + _INV_SQRT_2PI * (c_p * c_r * c_q) \
                    * (-x2_alpha_x3_beta_x4_gamma_mean) ** (b_exp - q) \
                    * _signed_pow(-sigma_x2, alpha_exp * q) \
                    * (0.5 * (1.0 - rho_x1x2 ** 2) * sigma_x1 ** 2) ** p \
                    * (-rho_x1x2 * sigma_x1) ** r \
                    * shifted ** (rem - r) \
                    * jnp.exp(mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q
                              + 0.5 * (1.0 - rho_x2x3_n ** 2) * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2
                              + 0.5 * (1.0 - rho_x2x4_n ** 2) * sigma_x4_n ** 2 * gamma_exp ** 2 * q ** 2
                              + (rho_x3x4_n - rho_x2x3_n * rho_x2x4_n)
                              * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp * q ** 2) \
                    * jnp.exp(0.25 * s_cc ** 2 - (mu_x2 / sigma_x2) * s_cc
                              + 0.5 * (mu_x2 ** 2 / sigma_x2 ** 2)) \
                    * _gamma_real(order) \
                    * _dvc(-order, s_cc)
    return total


def quadrivar_NNLL_MM_const_x1(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x2, sigma_x3_n, sigma_x4_n,
                               rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, x1_mean,
                               x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, a_exp, b_exp):
    """quadrivar_NNLL_MM with x1 constant. PDF_integrals_all_MM.F90:quadrivar_NNLL_MM_const_x1 (Eq. I2).

    x1 factors out as (μ_x1-<x1>)^a; the x2<0 D_v half-line integral and the x3,x4 lognormal raw moments remain.
    """
    total = 0.0
    for q in range(b_exp + 1):
        c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
        s_cc = (mu_x2 / sigma_x2 + rho_x2x3_n * sigma_x3_n * beta_exp * q
                + rho_x2x4_n * sigma_x4_n * gamma_exp * q)
        order = alpha_exp * q + 1.0
        total = total + _INV_SQRT_2PI * (mu_x1 - x1_mean) ** a_exp * c_q \
            * (-x2_alpha_x3_beta_x4_gamma_mean) ** (b_exp - q) \
            * _signed_pow(-sigma_x2, alpha_exp * q) \
            * jnp.exp(mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q
                      + 0.5 * (1.0 - rho_x2x3_n ** 2) * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2
                      + 0.5 * (1.0 - rho_x2x4_n ** 2) * sigma_x4_n ** 2 * gamma_exp ** 2 * q ** 2
                      + (rho_x3x4_n - rho_x2x3_n * rho_x2x4_n)
                      * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp * q ** 2) \
            * jnp.exp(0.25 * s_cc ** 2 - (mu_x2 / sigma_x2) * s_cc + 0.5 * (mu_x2 ** 2 / sigma_x2 ** 2)) \
            * _gamma_real(order) * _dvc(-order, s_cc)
    return total


def quadrivar_NNLL_MM_const_x2(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x1, sigma_x3_n, sigma_x4_n,
                               rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, x1_mean,
                               x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, a_exp, b_exp):
    """quadrivar_NNLL_MM with x2 constant. PDF_integrals_all_MM.F90:quadrivar_NNLL_MM_const_x2 (Eq. I3).

    x2 collapses to the delta at μ_x2 (no D_v); only nonzero when that delta is inside x2<0 (μ_x2<=0), where
    x2^α = signed_pow(μ_x2, αq). x1 is a full normal (p-sum, mean shifted by x1-lnx3/x1-lnx4), x3,x4 lognormal.
    """
    total = 0.0
    for p in range(a_exp // 2 + 1):
        c_p = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        for q in range(b_exp + 1):
            c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
            shifted = (mu_x1 - x1_mean + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp * q
                       + rho_x1x4_n * sigma_x1 * sigma_x4_n * gamma_exp * q)
            total = total + (c_p * c_q) \
                * (-x2_alpha_x3_beta_x4_gamma_mean) ** (b_exp - q) \
                * _signed_pow(mu_x2, alpha_exp * q) \
                * (0.5 * sigma_x1 ** 2) ** p \
                * shifted ** (a_exp - 2 * p) \
                * jnp.exp(mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q
                          + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2
                          + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2 * q ** 2
                          + rho_x3x4_n * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp * q ** 2)
    return jnp.where(mu_x2 <= 0.0, total, 0.0)


def quadrivar_NNLL_MM_const_x1x2(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x3_n, sigma_x4_n, rho_x3x4_n,
                                 x1_mean, x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp,
                                 a_exp, b_exp):
    """quadrivar_NNLL_MM with x1 AND x2 constant. PDF_integrals_all_MM.F90:quadrivar_NNLL_MM_const_x1x2 (Eq. I6).

    Only x3,x4 (lognormal) vary; nonzero only for μ_x2<=0 (the x2 delta inside x2<0), x2^α = signed_pow(μ_x2,αq).
    """
    total = 0.0
    for q in range(b_exp + 1):
        c_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
        total = total + (mu_x1 - x1_mean) ** a_exp * c_q \
            * (-x2_alpha_x3_beta_x4_gamma_mean) ** (b_exp - q) \
            * _signed_pow(mu_x2, alpha_exp * q) \
            * jnp.exp(mu_x3_n * beta_exp * q + mu_x4_n * gamma_exp * q
                      + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 * q ** 2
                      + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2 * q ** 2
                      + rho_x3x4_n * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp * q ** 2)
    return jnp.where(mu_x2 <= 0.0, total, 0.0)

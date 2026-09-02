"""JAX port of PDF_integrals_means.F90 — the normal-lognormal (NL) bivariate means.

These evaluate, analytically, the integral that the upscaled KK microphysics needs:

    < x1^alpha x2^beta >_i  =  INT_0^inf INT_-inf^inf  x1^alpha x2^beta P_NL,i(x1,x2) dx1 dx2

for the ith PDF component, where x1 = chi is normally distributed and x2 (= N_cn for
autoconversion, r_r for accretion) is lognormally distributed, with correlation rho in
normal/ln space. (x1^alpha is taken as 0 for x1 < 0, i.e. the integral is over x1 > 0.)

Oracle: PDF_integrals_means.F90, functions bivar_NL_mean (general, lines 481-559) and the
three reduced forms used when one or both component variances vanish:
  bivar_NL_mean_const_x1  (sigma_x1 -> 0; x1 constant)        lines 562-630
  bivar_NL_mean_const_x2  (sigma_x2 -> 0; x2 constant)        lines 633-697
  bivar_NL_mean_const_all (both -> 0; both constant)          lines 700-757

References (per the Fortran): Larson & Griffin 2013 (QJRMS, Eqs 33-44); Griffin 2016
(PhD thesis, Eqs C4-C17); Griffin & Larson 2016 (GMD supplements).

The only transcendental in the general/const_x2 forms is the parabolic cylinder function
D_v, supplied by parabolic_cylinder.dv_parabolic_cylinder.
"""
import jax.numpy as jnp

from clubb_jax.src.Microphys.KK_microphys.parabolic_cylinder import (
    _gamma_real,
    _dvc,          # clamped parabolic-cylinder D_v (its home, next to dv_parabolic_cylinder)
)

_INV_SQRT_2PI = 1.0 / jnp.sqrt(2.0 * jnp.pi)


def bivar_NL_mean(mu_x1, mu_x2_n, sigma_x1, sigma_x2_n,
                  rho_x1x2_n, alpha_exp, beta_exp):
    """General NL bivariate mean — both x1 and x2 vary. PDF_integrals_means.F90:481.

    bivar_NL_mean = (1/sqrt(2 pi)) sigma_x1^alpha
                    * exp( mu_x2_n beta + 1/2 sigma_x2_n^2 beta^2 - 1/4 s_c^2 )
                    * Gamma(alpha+1) * D_{-(alpha+1)}(-s_c),
      s_c = mu_x1/sigma_x1 + rho_x1x2_n sigma_x2_n beta.
    """
    s_c = (mu_x1 / sigma_x1) + rho_x1x2_n * sigma_x2_n * beta_exp
    return (_INV_SQRT_2PI * sigma_x1 ** alpha_exp
            * jnp.exp(mu_x2_n * beta_exp
                      + 0.5 * sigma_x2_n ** 2 * beta_exp ** 2
                      - 0.25 * s_c ** 2)
            * _gamma_real(alpha_exp + 1.0)
            * _dvc(-(alpha_exp + 1.0), -s_c))


def bivar_NL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp):
    """x1 constant (sigma_x1 -> 0). PDF_integrals_means.F90:562.

    = mu_x1^alpha exp( mu_x2_n beta + 1/2 sigma_x2_n^2 beta^2 )  for mu_x1 >= 0, else 0.
    """
    # Clamp the base before the fractional power so the (unused) mu_x1<0 branch does
    # not produce a complex value that poisons the jnp.where.
    base = jnp.maximum(mu_x1, 0.0)
    val = (base ** alpha_exp
           * jnp.exp(mu_x2_n * beta_exp + 0.5 * sigma_x2_n ** 2 * beta_exp ** 2))
    return jnp.where(mu_x1 >= 0.0, val, 0.0)


def bivar_NL_mean_const_x2(mu_x1, mu_x2, sigma_x1, alpha_exp, beta_exp):
    """x2 constant (sigma_x2 -> 0). PDF_integrals_means.F90:633.

    = (1/sqrt(2 pi)) sigma_x1^alpha mu_x2^beta exp( -1/4 mu_x1^2/sigma_x1^2 )
      * Gamma(alpha+1) * D_{-(alpha+1)}( -mu_x1/sigma_x1 ).
    """
    return (_INV_SQRT_2PI * sigma_x1 ** alpha_exp * mu_x2 ** beta_exp
            * jnp.exp(-0.25 * (mu_x1 ** 2 / sigma_x1 ** 2))
            * _gamma_real(alpha_exp + 1.0)
            * _dvc(-(alpha_exp + 1.0), -(mu_x1 / sigma_x1)))


def bivar_NL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp):
    """Both constant. PDF_integrals_means.F90:700.

    = mu_x1^alpha mu_x2^beta  for mu_x1 >= 0, else 0.
    """
    base = jnp.maximum(mu_x1, 0.0)
    return jnp.where(mu_x1 >= 0.0, base ** alpha_exp * mu_x2 ** beta_exp, 0.0)


# ============================================================================
# Trivariate normal-lognormal-lognormal (NLL) means — KK rain EVAPORATION.
#
# Evaporation integrates < chi^alpha r_r^beta N_r^gamma > over the SUBSATURATED
# half-line chi < 0 (chi here is supersaturation S; evaporation occurs where S<0),
# the mirror of the autoconversion/accretion (chi>0) integrals. Hence the chi
# factor appears as (-sigma_x1)^alpha / mu_x1^alpha (mu_x1<=0 branches) and the
# parabolic cylinder argument is +s_cc (not -s_c). x2=r_r, x3=N_r are lognormal.
# Oracle: PDF_integrals_means.F90:26-478 (general + 5 reduced forms).
# The KK evaporation exponent on chi is KK_evap_Supersat_exp = 1 for the rate MEANS,
# but the second-moment COVARIANCES use alpha+1 = 2 (even); _signed_pow must give the
# parity-correct (-1)^exp continuation (Fortran `(-x)**exp`), not sign(base) (odd-only).
# ============================================================================


def _signed_pow(base, exp):
    """base^exp via the principal real continuation for base<0: (-1)^exp |base|^exp,
    with (-1)^exp = cos(pi exp) (exact for integer exp). Matches Fortran `(-x)**exp`.
    For base>=0 it is |base|^exp; for odd exp it reduces to sign(base)|base|^exp (the
    evaporation alpha=1 case — means unchanged), and for even exp the sign is +, as the
    covariance alpha+1 term requires. exp is the integer chi exponent (1 or 2)."""
    return jnp.where(base < 0.0, jnp.cos(jnp.pi * exp), 1.0) * jnp.abs(base) ** exp


def _pos_pow(base, exp):
    """base^exp for base >= 0, grad-safe at base=0 (a hydrometeor mean). For positive exp,
    base^exp = 0 at base=0 with an inf gradient (exp*0^(exp-1)); the double-where keeps the
    FORWARD exact while never evaluating the power at 0 in the autodiff graph."""
    bp = jnp.where(base > 0.0, base, 1.0)
    return jnp.where(base > 0.0, bp ** exp, 0.0)


def trivar_NLL_mean(mu_x1, mu_x2_n, mu_x3_n, sigma_x1, sigma_x2_n, sigma_x3_n,
                    rho_x1x2_n, rho_x1x3_n, rho_x2x3_n,
                    alpha_exp, beta_exp, gamma_exp):
    """General NLL trivariate mean — all three vary. PDF_integrals_means.F90:26."""
    s_cc = ((mu_x1 / sigma_x1)
            + rho_x1x2_n * sigma_x2_n * beta_exp
            + rho_x1x3_n * sigma_x3_n * gamma_exp)
    return (_INV_SQRT_2PI * _signed_pow(-sigma_x1, alpha_exp)
            * jnp.exp(mu_x2_n * beta_exp + mu_x3_n * gamma_exp)
            * jnp.exp(0.5 * ((1.0 - rho_x1x2_n ** 2) * sigma_x2_n ** 2 * beta_exp ** 2
                             + (1.0 - rho_x1x3_n ** 2) * sigma_x3_n ** 2 * gamma_exp ** 2
                             + 2.0 * (rho_x2x3_n - rho_x1x2_n * rho_x1x3_n)
                             * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp))
            * jnp.exp(0.25 * s_cc ** 2 - (mu_x1 / sigma_x1) * s_cc
                      + 0.5 * (mu_x1 ** 2 / sigma_x1 ** 2))
            * _gamma_real(alpha_exp + 1.0)
            * _dvc(-(alpha_exp + 1.0), s_cc))


def trivar_NLL_mean_const_x1(mu_x1, mu_x2_n, mu_x3_n, sigma_x2_n, sigma_x3_n,
                             rho_x2x3_n, alpha_exp, beta_exp, gamma_exp):
    """chi constant (sigma_x1->0). PDF_integrals_means.F90:123. Nonzero only mu_x1<=0."""
    val = (_signed_pow(mu_x1, alpha_exp)
           * jnp.exp(mu_x2_n * beta_exp + mu_x3_n * gamma_exp
                     + 0.5 * sigma_x2_n ** 2 * beta_exp ** 2
                     + 0.5 * sigma_x3_n ** 2 * gamma_exp ** 2
                     + rho_x2x3_n * sigma_x2_n * beta_exp * sigma_x3_n * gamma_exp))
    return jnp.where(mu_x1 <= 0.0, val, 0.0)


def trivar_NLL_mean_const_x2(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x3_n,
                             rho_x1x3_n, alpha_exp, beta_exp, gamma_exp):
    """x2 (r_r) constant (sigma_x2->0). PDF_integrals_means.F90:200."""
    s_cc = (mu_x1 / sigma_x1) + rho_x1x3_n * sigma_x3_n * gamma_exp
    return (_INV_SQRT_2PI * _signed_pow(-sigma_x1, alpha_exp)
            * _pos_pow(mu_x2, beta_exp)
            * jnp.exp(mu_x3_n * gamma_exp
                      + 0.5 * sigma_x3_n ** 2 * gamma_exp ** 2
                      - 0.25 * s_cc ** 2)
            * _gamma_real(alpha_exp + 1.0)
            * _dvc(-(alpha_exp + 1.0), s_cc))


def trivar_NLL_mean_const_x1x2(mu_x1, mu_x2, mu_x3_n, sigma_x3_n,
                               alpha_exp, beta_exp, gamma_exp):
    """chi and x2 constant. PDF_integrals_means.F90:281. Nonzero only mu_x1<=0."""
    val = (_signed_pow(mu_x1, alpha_exp) * _pos_pow(mu_x2, beta_exp)
           * jnp.exp(mu_x3_n * gamma_exp + 0.5 * sigma_x3_n ** 2 * gamma_exp ** 2))
    return jnp.where(mu_x1 <= 0.0, val, 0.0)


def trivar_NLL_mean_const_x2x3(mu_x1, mu_x2, mu_x3, sigma_x1,
                               alpha_exp, beta_exp, gamma_exp):
    """x2 and x3 constant. PDF_integrals_means.F90:349."""
    return (_INV_SQRT_2PI * _signed_pow(-sigma_x1, alpha_exp)
            * _pos_pow(mu_x2, beta_exp) * _pos_pow(mu_x3, gamma_exp)
            * jnp.exp(-0.25 * (mu_x1 ** 2 / sigma_x1 ** 2))
            * _gamma_real(alpha_exp + 1.0)
            * _dvc(-(alpha_exp + 1.0), (mu_x1 / sigma_x1)))


def trivar_NLL_mean_const_all(mu_x1, mu_x2, mu_x3, alpha_exp, beta_exp, gamma_exp):
    """All constant. PDF_integrals_means.F90:418. Nonzero only mu_x1<=0."""
    val = _signed_pow(mu_x1, alpha_exp) * _pos_pow(mu_x2, beta_exp) * _pos_pow(mu_x3, gamma_exp)
    return jnp.where(mu_x1 <= 0.0, val, 0.0)


# ============================================================================
# Bivariate lognormal-lognormal (LL) means — KK rain mean-volume-radius (mvr).
#
# Both x1=r_r and x2=N_r are lognormal, so < x1^alpha x2^beta > is the exact
# lognormal moment-generating result — a pure exponential, no special function.
# Oracle: PDF_integrals_means.F90:760-919.
# ============================================================================


def bivar_LL_mean(mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n,
                  rho_x1x2_n, alpha_exp, beta_exp):
    """General LL bivariate mean — both vary. PDF_integrals_means.F90:760.

    = exp( mu_x1_n alpha + mu_x2_n beta + 1/2 sigma_x1_n^2 alpha^2
           + 1/2 sigma_x2_n^2 beta^2 + rho sigma_x1_n alpha sigma_x2_n beta )."""
    return jnp.exp(mu_x1_n * alpha_exp + mu_x2_n * beta_exp
                   + 0.5 * sigma_x1_n ** 2 * alpha_exp ** 2
                   + 0.5 * sigma_x2_n ** 2 * beta_exp ** 2
                   + rho_x1x2_n * sigma_x1_n * alpha_exp * sigma_x2_n * beta_exp)


def bivar_LL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp):
    """x1 constant (sigma_x1->0). PDF_integrals_means.F90:831.

    = mu_x1^alpha exp( mu_x2_n beta + 1/2 sigma_x2_n^2 beta^2 ).
    (No mu_x1 sign guard in the oracle: x1=r_r/N_r is a positive hydrometeor mean.)"""
    return (mu_x1 ** alpha_exp
            * jnp.exp(mu_x2_n * beta_exp + 0.5 * sigma_x2_n ** 2 * beta_exp ** 2))


def bivar_LL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp):
    """Both constant. PDF_integrals_means.F90:880.  = mu_x1^alpha mu_x2^beta."""
    return mu_x1 ** alpha_exp * mu_x2 ** beta_exp

"""Normal/lognormal covariance integrals — JAX port of PDF_integrals_covar.F90.

These closed forms give the component covariance Cov_i(x1, x2^alpha x3^beta) (and the 4-variable
analogue) for the KK upscaled second-moment microphysics tendencies (KK_upscaled_covar_driver). They
reuse the same parabolic-cylinder `Dv` + gamma primitives as the already-ported mean integrals
(PDF_integrals_means.py). Ported incrementally with Monte-Carlo + oracle validation (see DESIGN.md
"KK microphysics"). Differentiable.

The optional ``dv`` pair supplies D_{-(alpha+1)} and D_{-(alpha+2)} from the
covariance driver's shared batch. Without it, each integral evaluates Dv locally.
"""
import math

import jax.numpy as jnp

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import _signed_pow
from clubb_jax.src.Microphys.KK_microphys.parabolic_cylinder import _dvc
from clubb_jax.src.Microphys.KK_microphys.parabolic_cylinder import _gamma_real

# Imports may occur inside JIT tracing; module constants must remain host scalars.
_SQRT_2PI = math.sqrt(2.0 * math.pi)


def quadrivar_NNLL_covar(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x1, sigma_x2, sigma_x3_n, sigma_x4_n,
                         rho_x1x2, rho_x1x3_n, rho_x1x4_n, rho_x2x3_n, rho_x2x4_n, rho_x3x4_n,
                         x1_mean, x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, *, dv=None):
    """Cov_i(x1, x2^α x3^β x4^γ) for x1,x2 ~ normal and x3,x4 ~ lognormal — the SUBSATURATED (chi<0)
    region form used by the KK evaporation / w covariances. Faithful port of
    PDF_integrals_covar.F90:quadrivar_NNLL_covar (base, all-varying); the `(-σ_x2)^α` is the signed
    continuation `_signed_pow(-σ_x2, α)` (same convention as the validated trivar_NLL_mean), and the
    parabolic-cylinder args use +s_cc (not -s_c)."""
    s_cc = (mu_x2 / sigma_x2
            + rho_x2x3_n * sigma_x3_n * beta_exp
            + rho_x2x4_n * sigma_x4_n * gamma_exp)
    d_alpha2 = _dvc(-(alpha_exp + 2.0), s_cc) if dv is None else dv[1]
    d_alpha1 = _dvc(-(alpha_exp + 1.0), s_cc) if dv is None else dv[0]
    return (1.0 / _SQRT_2PI
            * _signed_pow(-sigma_x2, alpha_exp)
            * jnp.exp(mu_x3_n * beta_exp + mu_x4_n * gamma_exp
                      + 0.5 * (1.0 - rho_x2x3_n ** 2) * sigma_x3_n ** 2 * beta_exp ** 2
                      + 0.5 * (1.0 - rho_x2x4_n ** 2) * sigma_x4_n ** 2 * gamma_exp ** 2
                      + (rho_x3x4_n - rho_x2x3_n * rho_x2x4_n)
                      * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp)
            * jnp.exp(0.25 * s_cc ** 2 - (mu_x2 / sigma_x2) * s_cc
                      + 0.5 * (mu_x2 ** 2 / sigma_x2 ** 2))
            * (-rho_x1x2 * sigma_x1 * _gamma_real(alpha_exp + 2.0) * d_alpha2
               + (mu_x1 - x1_mean
                  - (mu_x2 / sigma_x2) * rho_x1x2 * sigma_x1
                  + (rho_x1x3_n - rho_x1x2 * rho_x2x3_n) * sigma_x1 * sigma_x3_n * beta_exp
                  + (rho_x1x4_n - rho_x1x2 * rho_x2x4_n) * sigma_x1 * sigma_x4_n * gamma_exp)
               * _gamma_real(alpha_exp + 1.0) * d_alpha1)
            - x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_const_x1(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x2, sigma_x3_n, sigma_x4_n,
                                  rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, x1_mean,
                                  x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, *, dv=None):
    """x1 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x1."""
    s_cc = (mu_x2 / sigma_x2 + rho_x2x3_n * sigma_x3_n * beta_exp + rho_x2x4_n * sigma_x4_n * gamma_exp)
    d_alpha1 = _dvc(-(alpha_exp + 1.0), s_cc) if dv is None else dv[0]
    return ((1.0 / _SQRT_2PI) * (mu_x1 - x1_mean)
            * _signed_pow(-sigma_x2, alpha_exp)
            * jnp.exp(mu_x3_n * beta_exp + mu_x4_n * gamma_exp
                      + 0.5 * (1.0 - rho_x2x3_n ** 2) * sigma_x3_n ** 2 * beta_exp ** 2
                      + 0.5 * (1.0 - rho_x2x4_n ** 2) * sigma_x4_n ** 2 * gamma_exp ** 2
                      + (rho_x3x4_n - rho_x2x3_n * rho_x2x4_n) * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp)
            * jnp.exp(0.25 * s_cc ** 2 - (mu_x2 / sigma_x2) * s_cc + 0.5 * (mu_x2 ** 2 / sigma_x2 ** 2))
            * _gamma_real(alpha_exp + 1.0) * d_alpha1
            - x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_const_x2(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x1, sigma_x3_n, sigma_x4_n,
                                  rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, x1_mean,
                                  x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp):
    """x2 (chi) constant. Nonzero only for mu_x2<=0 (subsaturated). signed_pow(mu_x2,α) for the chi<0 branch.
    PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x2."""
    pos = (_signed_pow(mu_x2, alpha_exp)
           * (mu_x1 - x1_mean + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp
              + rho_x1x4_n * sigma_x1 * sigma_x4_n * gamma_exp)
           * jnp.exp(mu_x3_n * beta_exp + mu_x4_n * gamma_exp
                     + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2
                     + rho_x3x4_n * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp)
           - x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean))
    return jnp.where(mu_x2 <= 0.0, pos, -x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_const_x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x1, sigma_x2, sigma_x4_n,
                                  rho_x1x2, rho_x1x4_n, rho_x2x4_n, x1_mean,
                                  x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, *, dv=None):
    """x3 (r_r) constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x3."""
    s_cc = mu_x2 / sigma_x2 + rho_x2x4_n * sigma_x4_n * gamma_exp
    d_alpha2 = _dvc(-(alpha_exp + 2.0), s_cc) if dv is None else dv[1]
    d_alpha1 = _dvc(-(alpha_exp + 1.0), s_cc) if dv is None else dv[0]
    return (1.0 / _SQRT_2PI
            * _signed_pow(-sigma_x2, alpha_exp) * mu_x3 ** beta_exp
            * jnp.exp(mu_x4_n * gamma_exp + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2 - 0.25 * s_cc ** 2)
            * (-rho_x1x2 * sigma_x1 * _gamma_real(alpha_exp + 2.0) * d_alpha2
               + (mu_x1 - x1_mean - (mu_x2 / sigma_x2) * rho_x1x2 * sigma_x1
                  + (rho_x1x4_n - rho_x1x2 * rho_x2x4_n) * sigma_x1 * sigma_x4_n * gamma_exp)
               * _gamma_real(alpha_exp + 1.0) * d_alpha1)
            - x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean))


def _q_const_neg(pos, mu_x1, mu_x2, x1_mean, x2a):
    """Subsaturated multi-const dispatch: the `pos` form for mu_x2<=0, else the −overall-mean term."""
    return jnp.where(mu_x2 <= 0.0, pos, -x2a * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_const_x1x2(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x3_n, sigma_x4_n, rho_x3x4_n,
                                    x1_mean, x2a, alpha_exp, beta_exp, gamma_exp):
    """x1, x2 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x1x2."""
    pos = ((mu_x1 - x1_mean) * _signed_pow(mu_x2, alpha_exp)
           * jnp.exp(mu_x3_n * beta_exp + mu_x4_n * gamma_exp
                     + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2
                     + rho_x3x4_n * sigma_x3_n * beta_exp * sigma_x4_n * gamma_exp)
           - x2a * (mu_x1 - x1_mean))
    return _q_const_neg(pos, mu_x1, mu_x2, x1_mean, x2a)


def quadrivar_NNLL_covar_const_x1x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x2, sigma_x4_n, rho_x2x4_n,
                                    x1_mean, x2a, alpha_exp, beta_exp, gamma_exp, *, dv=None):
    """x1, x3 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x1x3."""
    s_cc = mu_x2 / sigma_x2 + rho_x2x4_n * sigma_x4_n * gamma_exp
    d_alpha1 = _dvc(-(alpha_exp + 1.0), s_cc) if dv is None else dv[0]
    return ((1.0 / _SQRT_2PI) * (mu_x1 - x1_mean) * _signed_pow(-sigma_x2, alpha_exp) * mu_x3 ** beta_exp
            * jnp.exp(mu_x4_n * gamma_exp + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2 - 0.25 * s_cc ** 2)
            * _gamma_real(alpha_exp + 1.0) * d_alpha1
            - x2a * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_const_x2x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x1, sigma_x4_n, rho_x1x4_n,
                                    x1_mean, x2a, alpha_exp, beta_exp, gamma_exp):
    """x2, x3 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x2x3."""
    pos = (_signed_pow(mu_x2, alpha_exp) * mu_x3 ** beta_exp
           * (mu_x1 - x1_mean + rho_x1x4_n * sigma_x1 * sigma_x4_n * gamma_exp)
           * jnp.exp(mu_x4_n * gamma_exp + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2)
           - x2a * (mu_x1 - x1_mean))
    return _q_const_neg(pos, mu_x1, mu_x2, x1_mean, x2a)


def quadrivar_NNLL_covar_const_x3x4(mu_x1, mu_x2, mu_x3, mu_x4, sigma_x1, sigma_x2, rho_x1x2,
                                    x1_mean, x2a, alpha_exp, beta_exp, gamma_exp, *, dv=None):
    """x3, x4 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_const_x3x4."""
    r = mu_x2 / sigma_x2
    d_alpha2 = _dvc(-(alpha_exp + 2.0), r) if dv is None else dv[1]
    d_alpha1 = _dvc(-(alpha_exp + 1.0), r) if dv is None else dv[0]
    return (1.0 / _SQRT_2PI * _signed_pow(-sigma_x2, alpha_exp) * mu_x3 ** beta_exp * mu_x4 ** gamma_exp
            * jnp.exp(-0.25 * (mu_x2 ** 2 / sigma_x2 ** 2))
            * (-rho_x1x2 * sigma_x1 * _gamma_real(alpha_exp + 2.0) * d_alpha2
               + (mu_x1 - x1_mean - r * rho_x1x2 * sigma_x1)
               * _gamma_real(alpha_exp + 1.0) * d_alpha1)
            - x2a * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_cst_x1x2x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x4_n,
                                    x1_mean, x2a, alpha_exp, beta_exp, gamma_exp):
    """x1, x2, x3 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_cst_x1x2x3."""
    pos = ((mu_x1 - x1_mean) * _signed_pow(mu_x2, alpha_exp) * mu_x3 ** beta_exp
           * jnp.exp(mu_x4_n * gamma_exp + 0.5 * sigma_x4_n ** 2 * gamma_exp ** 2)
           - x2a * (mu_x1 - x1_mean))
    return _q_const_neg(pos, mu_x1, mu_x2, x1_mean, x2a)


def quadrivar_NNLL_covar_cst_x1x3x4(mu_x1, mu_x2, mu_x3, mu_x4, sigma_x2,
                                    x1_mean, x2a, alpha_exp, beta_exp, gamma_exp, *, dv=None):
    """x1, x3, x4 constant. PDF_integrals_covar.F90:quadrivar_NNLL_covar_cst_x1x3x4."""
    d_alpha1 = _dvc(-(alpha_exp + 1.0), mu_x2 / sigma_x2) if dv is None else dv[0]
    return ((1.0 / _SQRT_2PI) * (mu_x1 - x1_mean) * _signed_pow(-sigma_x2, alpha_exp)
            * mu_x3 ** beta_exp * mu_x4 ** gamma_exp
            * jnp.exp(-0.25 * (mu_x2 ** 2 / sigma_x2 ** 2))
            * _gamma_real(alpha_exp + 1.0) * d_alpha1
            - x2a * (mu_x1 - x1_mean))


def quadrivar_NNLL_covar_cst_x2x3x4(mu_x1, mu_x2, mu_x3, mu_x4, x1_mean, x2a, alpha_exp, beta_exp, gamma_exp):
    """x2, x3, x4 constant (== the all-constant case). PDF_integrals_covar.F90:quadrivar_NNLL_covar_cst_x2x3x4."""
    pos = (mu_x1 - x1_mean) * (_signed_pow(mu_x2, alpha_exp) * mu_x3 ** beta_exp * mu_x4 ** gamma_exp - x2a)
    return _q_const_neg(pos, mu_x1, mu_x2, x1_mean, x2a)


# All four constant: identical formula to cst_x2x3x4.
quadrivar_NNLL_covar_const_all = quadrivar_NNLL_covar_cst_x2x3x4




def trivar_NNL_covar(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x2, sigma_x3_n,
                     rho_x1x2, rho_x1x3_n, rho_x2x3_n,
                     x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, *, dv=None):
    """Cov_i(x1, x2^alpha x3^beta) for x1,x2 ~ normal and x3 ~ lognormal (ln x3 normal), with the
    component correlations. `x1_mean` / `x2_alpha_x3_beta_mean` are the OVERALL means. Faithful port
    of PDF_integrals_covar.F90:trivar_NNL_covar (the all-varying base case; sigma_x2 > 0)."""
    s_c = (mu_x2 / sigma_x2) + rho_x2x3_n * sigma_x3_n * beta_exp
    d_alpha2 = _dvc(-(alpha_exp + 2.0), -s_c) if dv is None else dv[1]
    d_alpha1 = _dvc(-(alpha_exp + 1.0), -s_c) if dv is None else dv[0]
    return (1.0 / _SQRT_2PI
            * sigma_x2 ** alpha_exp
            * jnp.exp(mu_x3_n * beta_exp
                      + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2
                      - 0.25 * s_c ** 2)
            * (rho_x1x2 * sigma_x1 * _gamma_real(alpha_exp + 2.0)
               * d_alpha2
               + (mu_x1 - x1_mean
                  - (mu_x2 / sigma_x2) * rho_x1x2 * sigma_x1
                  + (rho_x1x3_n - rho_x1x2 * rho_x2x3_n)
                  * sigma_x1 * sigma_x3_n * beta_exp)
               * _gamma_real(alpha_exp + 1.0)
               * d_alpha1)
            - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))


def trivar_NNL_covar_const_x1(mu_x1, mu_x2, mu_x3_n, sigma_x2, sigma_x3_n,
                              rho_x2x3_n, x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, *, dv=None):
    """x1 constant within the component (sigma_x1=0). PDF_integrals_covar.F90:trivar_NNL_covar_const_x1."""
    s_c = (mu_x2 / sigma_x2) + rho_x2x3_n * sigma_x3_n * beta_exp
    d_alpha1 = _dvc(-(alpha_exp + 1.0), -s_c) if dv is None else dv[0]
    return ((1.0 / _SQRT_2PI) * (mu_x1 - x1_mean)
            * sigma_x2 ** alpha_exp
            * jnp.exp(mu_x3_n * beta_exp + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2 - 0.25 * s_c ** 2)
            * _gamma_real(alpha_exp + 1.0) * d_alpha1
            - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))


def trivar_NNL_covar_const_x2(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x3_n,
                              rho_x1x3_n, x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp):
    """x2 (chi) constant within the component (sigma_x2=0). PDF_integrals_covar.F90:trivar_NNL_covar_const_x2.
    For mu_x2<0 the x2^alpha integrand vanishes (subsaturated), leaving only the -overall-mean term."""
    pos = (mu_x2 ** alpha_exp
           * (mu_x1 - x1_mean + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp)
           * jnp.exp(mu_x3_n * beta_exp + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2)
           - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))
    neg = -x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)
    mu_x2_safe = jnp.where(mu_x2 >= 0.0, mu_x2, 1.0)   # avoid NaN^alpha in the unused branch
    pos = (mu_x2_safe ** alpha_exp
           * (mu_x1 - x1_mean + rho_x1x3_n * sigma_x1 * sigma_x3_n * beta_exp)
           * jnp.exp(mu_x3_n * beta_exp + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2)
           - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))
    return jnp.where(mu_x2 >= 0.0, pos, neg)


def trivar_NNL_covar_const_x3(mu_x1, mu_x2, mu_x3, sigma_x1, sigma_x2,
                              rho_x1x2, x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, *, dv=None):
    """x3 (N_cn or r_r) constant within the component (sigma_x3_n=0; x3=mu_x3). The l_const_Nc path.
    PDF_integrals_covar.F90:trivar_NNL_covar_const_x3."""
    r = mu_x2 / sigma_x2
    d_alpha2 = _dvc(-(alpha_exp + 2.0), -r) if dv is None else dv[1]
    d_alpha1 = _dvc(-(alpha_exp + 1.0), -r) if dv is None else dv[0]
    return (1.0 / _SQRT_2PI
            * sigma_x2 ** alpha_exp * mu_x3 ** beta_exp
            * jnp.exp(-0.25 * (mu_x2 ** 2 / sigma_x2 ** 2))
            * (rho_x1x2 * sigma_x1 * _gamma_real(alpha_exp + 2.0) * d_alpha2
               + (mu_x1 - x1_mean - r * rho_x1x2 * sigma_x1)
               * _gamma_real(alpha_exp + 1.0) * d_alpha1)
            - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))


def _x2pow(mu_x2, alpha_exp):
    """mu_x2^alpha, NaN-safe in the unused mu_x2<0 branch (the caller's jnp.where discards it)."""
    return jnp.where(mu_x2 >= 0.0, mu_x2, 1.0) ** alpha_exp


def trivar_NNL_covar_const_x1x2(mu_x1, mu_x2, mu_x3_n, sigma_x3_n,
                                x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp):
    """x1 and x2 constant. PDF_integrals_covar.F90:trivar_NNL_covar_const_x1x2."""
    pos = (_x2pow(mu_x2, alpha_exp) * (mu_x1 - x1_mean)
           * jnp.exp(mu_x3_n * beta_exp + 0.5 * sigma_x3_n ** 2 * beta_exp ** 2)
           - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))
    return jnp.where(mu_x2 >= 0.0, pos, -x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))


def trivar_NNL_covar_const_x1x3(mu_x1, mu_x2, mu_x3, sigma_x2,
                                x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, *, dv=None):
    """x1 and x3 constant. PDF_integrals_covar.F90:trivar_NNL_covar_const_x1x3."""
    d_alpha1 = _dvc(-(alpha_exp + 1.0), -(mu_x2 / sigma_x2)) if dv is None else dv[0]
    return ((1.0 / _SQRT_2PI) * (mu_x1 - x1_mean)
            * sigma_x2 ** alpha_exp * mu_x3 ** beta_exp
            * jnp.exp(-0.25 * (mu_x2 ** 2 / sigma_x2 ** 2))
            * _gamma_real(alpha_exp + 1.0) * d_alpha1
            - x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))


def trivar_NNL_covar_const_x2x3(mu_x1, mu_x2, mu_x3, x1_mean,
                                x2_alpha_x3_beta_mean, alpha_exp, beta_exp):
    """x2 and x3 constant (== the all-constant case). PDF_integrals_covar.F90:trivar_NNL_covar_const_x2x3."""
    pos = (mu_x1 - x1_mean) * (_x2pow(mu_x2, alpha_exp) * mu_x3 ** beta_exp - x2_alpha_x3_beta_mean)
    return jnp.where(mu_x2 >= 0.0, pos, -x2_alpha_x3_beta_mean * (mu_x1 - x1_mean))


# All three constant: identical formula to const_x2x3 (PDF_integrals_covar.F90:trivar_NNL_covar_const_all).
trivar_NNL_covar_const_all = trivar_NNL_covar_const_x2x3


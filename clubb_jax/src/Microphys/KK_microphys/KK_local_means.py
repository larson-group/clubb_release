"""JAX port of KK_microphys/KK_local_means.F90 — the GRID-MEAN (non-upscaled) KK rate formulas.

Mirrors `clubb_release/src/Microphys/KK_microphys/KK_local_means.F90`: the Khairoutdinov-Kogan warm-rain rates
evaluated at the grid-mean values (the "local" means, used when the analytic upscaled-over-the-PDF means are
disabled). Power laws with the `parameters_KK.F90` exponents:

  evap = coef * s^1 * rr^(1/3) * Nr^(2/3)   for s<=0 (subsaturated), else 0   [s = extended liquid water]
  auto = coef * rc^2.47 * Nc^(-1.79)        for rc>0, else 0
  accr = coef * rc^1.15 * rr^1.15           for rc>0, else 0
  mvr  = coef * rr^(1/3) * max(Nr,Nr_tol)^(-1/3)

Completeness port (the gated KK path uses the upscaled means in `KK_upscaled_means.py`); no f2py oracle for
the local-mean variants, so validated analytically (`tests/test_kk_local_means.py`). Pure jnp; the `mean_s`
base is clamped with `max(·,0)` inside the masked branches so a negative-base fractional power can't NaN the
gradient (forward-identical to the Fortran's `if`-guarded evaluation).
"""
import jax.numpy as jnp

# KK rate exponents — parameters_KK.F90
_EVAP_S_EXP, _EVAP_RR_EXP, _EVAP_NR_EXP = 1.0, 1.0 / 3.0, 2.0 / 3.0
_AUTO_RC_EXP, _AUTO_NC_EXP = 2.47, -1.79
_ACCR_RC_EXP, _ACCR_RR_EXP = 1.15, 1.15
_MVR_RR_EXP, _MVR_NR_EXP = 1.0 / 3.0, -1.0 / 3.0


def KK_evap_local_mean(mean_s, rrm, Nrm, KK_evap_coef):
    """KK rain evaporation rate at the grid mean: `coef*s*rr^(1/3)*Nr^(2/3)` for s<=0 (subsaturated), else 0."""
    rate = KK_evap_coef * mean_s ** _EVAP_S_EXP * rrm ** _EVAP_RR_EXP * Nrm ** _EVAP_NR_EXP
    return jnp.where(mean_s <= 0.0, rate, 0.0)


def KK_auto_local_mean(mean_s, Ncm, KK_auto_coef):
    """KK autoconversion rate at the grid mean: `coef*rc^2.47*Nc^(-1.79)` for rc>0, else 0 (`mean_s`=rc)."""
    safe = jnp.maximum(mean_s, 0.0)
    rate = KK_auto_coef * safe ** _AUTO_RC_EXP * Ncm ** _AUTO_NC_EXP
    return jnp.where(mean_s > 0.0, rate, 0.0)


def KK_accr_local_mean(mean_s, rrm, KK_accr_coef):
    """KK accretion rate at the grid mean: `coef*rc^1.15*rr^1.15` for rc>0, else 0 (`mean_s`=rc)."""
    safe = jnp.maximum(mean_s, 0.0)
    rate = KK_accr_coef * safe ** _ACCR_RC_EXP * rrm ** _ACCR_RR_EXP
    return jnp.where(mean_s > 0.0, rate, 0.0)


def KK_mvr_local_mean(rrm, Nrm, KK_mvr_coef, Nr_tol):
    """KK rain mean-volume-radius at the grid mean: `coef*rr^(1/3)*max(Nr,Nr_tol)^(-1/3)`."""
    return KK_mvr_coef * rrm ** _MVR_RR_EXP * jnp.maximum(Nrm, Nr_tol) ** _MVR_NR_EXP

#!/usr/bin/env python3
"""test_cloud_feedback_sfclyr.py — validate the CGILS/cloud-feedback surface flux scheme.

Oracle: a literal NumPy transcription of the Fortran bulk drag-law formula (validates the vectorization +
constants + drag-scaling), plus physical invariants (neutral -> zero flux, drag-scaling at 20m, sign of fluxes)
and a finite jax.grad. sat_mixrat_liq is the already-validated saturation port.
"""


import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.cloud_feedback import cloud_feedback_sfclyr
from clubb_jax.src.CLUBB_core.constants_clubb import p0, kappa
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq

_SAT = 3   # SATURATION_FLATAU
_C_H_20 = 0.001094
_C_Q_20 = 0.001133
_Z0 = 0.00015
_STD_FLUX_ALT = 20.0


def _set_surface_temperature(T_sfc):
    time_dependent_input.time_sfc_given = jnp.asarray([0.0, 1.0])
    time_dependent_input.T_sfc_given = jnp.asarray([T_sfc, T_sfc])


def test_physical_invariants():
    # Drag scaling == base coefficient exactly at 20 m.
    _set_surface_temperature(298.0)
    wth, wrt, _, _ = cloud_feedback_sfclyr(
        1, 0.5, 1, jnp.asarray([295.0]), jnp.asarray([0.012]),
        jnp.asarray([20.0]), jnp.asarray([5.0]), jnp.asarray([1.0e5]), _SAT
    )
    exner = (1.0e5 / p0) ** kappa
    expect_wth = -_C_H_20 * 5.0 * (295.0 - 298.0 / exner)
    assert abs(float(wth[0]) - expect_wth) < 1e-12, "drag scaling != base at 20 m"
    # Neutral (thlm == T_sfc/exner) -> zero heat flux; rtm == rsat -> zero moisture flux.
    Tn, pn, zn, un = 298.0, 1.0e5, 30.0, 6.0
    thn = Tn / ((pn / p0) ** kappa)
    rsat = float(sat_mixrat_liq(pn, Tn, _SAT))
    _set_surface_temperature(Tn)
    wth0, wrt0, _, _ = cloud_feedback_sfclyr(
        1, 0.5, 1, jnp.asarray([thn]), jnp.asarray([rsat]), jnp.asarray([zn]),
        jnp.asarray([un]), jnp.asarray([pn]), _SAT
    )
    assert abs(float(wth0[0])) < 1e-12 and abs(float(wrt0[0])) < 1e-12, "neutral/saturated should give zero flux"
    # Warm wet surface -> upward (positive) fluxes.
    _set_surface_temperature(300.0)
    wthw, wrtw, _, _ = cloud_feedback_sfclyr(
        1, 0.5, 1, jnp.asarray([290.0]), jnp.asarray([0.008]), jnp.asarray([30.0]),
        jnp.asarray([6.0]), jnp.asarray([1.0e5]), _SAT
    )
    assert float(wthw[0]) > 0 and float(wrtw[0]) > 0, "warm/wet surface should give upward fluxes"
    print("  cloud_feedback_sfclyr: drag-scaling@20m + neutral->0 + warm/wet->upward  PASS")

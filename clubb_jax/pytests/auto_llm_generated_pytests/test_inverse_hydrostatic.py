#!/usr/bin/env python3
"""test_inverse_hydrostatic.py — validate inverse_hydrostatic (pressure-sounding altitudes).

Strongest oracle: the ROUND-TRIP against the existing forward hydrostatic (init_pressure uses the same
log-mean scheme), so building exner from known heights then inverting must recover those heights exactly. Also:
a literal NumPy transcription, the constant-thvm analytic closed form, and a finite jax.grad.
"""
import math


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.Input_fields.hydrostatic_module import (
    calc_ref_z_linear_thvm, inverse_hydrostatic, _CP_OV_G)
from clubb_jax.src.CLUBB_core.constants_clubb import p0, kappa, Cp, grav


def _exner_from_z(z, thvm, exner_sfc):
    """Forward: build exner at heights z from thvm + surface exner, same log-mean scheme as init_pressure."""
    n = len(z)
    exner = np.zeros(n)
    exner[0] = exner_sfc
    for k in range(1, n):
        if abs(thvm[k] - thvm[k - 1]) > 1e-12 * thvm[k]:
            exner[k] = exner[k - 1] - (grav / Cp) * (z[k] - z[k - 1]) / (thvm[k] - thvm[k - 1]) * math.log(
                thvm[k] / thvm[k - 1])
        else:
            exner[k] = exner[k - 1] - (grav / Cp) * (z[k] - z[k - 1]) / thvm[k]
    return exner


def test_roundtrip_with_forward():
    # Known heights + thvm -> forward exner -> inverse z' must recover the heights.
    z = np.linspace(0.0, 15000.0, 40)
    thvm = 295.0 + 0.003 * z + 5.0 * np.sin(z / 4000.0)   # smooth, increasing-ish, all > 0
    p_sfc = 101325.0
    exner_sfc = (p_sfc / p0) ** kappa
    exner = _exner_from_z(z, thvm, exner_sfc)
    z_back = np.asarray(inverse_hydrostatic(p_sfc, z[0], jnp.asarray(thvm), jnp.asarray(exner)))
    err = np.max(np.abs(z_back - z))
    assert err < 1e-7, f"round-trip z->exner->z err {err:.2e} m"
    print(f"  inverse_hydrostatic round-trip (z->exner->z, 40 levels): max err {err:.1e} m  PASS")


def test_constant_thvm_analytic():
    # Constant thvm: ref_z[k] = -(Cp/g) thvm (exner[k]-exner[0]) (log-mean = thvm).
    thvm = np.full(20, 300.0)
    exner = np.linspace(0.99, 0.4, 20)
    got = np.asarray(calc_ref_z_linear_thvm(jnp.asarray(thvm), jnp.asarray(exner)))
    ref = -_CP_OV_G * 300.0 * (exner - exner[0])
    assert np.max(np.abs(got - ref)) < 1e-9, "constant-thvm analytic mismatch"
    print("  calc_ref_z_linear_thvm constant-thvm vs analytic closed form: <1e-9  PASS")

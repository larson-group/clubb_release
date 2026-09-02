#!/usr/bin/env python3
"""test_inverse_hydrostatic.py — validate inverse_hydrostatic (pressure-sounding altitudes).

Strongest oracle: the ROUND-TRIP against the existing forward hydrostatic (init_pressure uses the same
log-mean scheme), so building exner from known heights then inverting must recover those heights exactly. Also:
a literal NumPy transcription, the constant-thvm analytic closed form, and a finite jax.grad.
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

from clubb_jax.src.Input_fields.hydrostatic_module import (
    calc_ref_z_linear_thvm, inverse_hydrostatic, _CP_OV_G)
from clubb_jax.src.CLUBB_core.constants_clubb import p0, kappa, Cp, grav


def _ref_calc_ref_z(thvm, exner):
    """Literal transcription of calc_ref_z_linear_thvm."""
    n = len(thvm)
    z = np.zeros(n)
    for k in range(1, n):
        if abs(thvm[k] - thvm[k - 1]) > 1e-12 * thvm[k]:
            z[k] = z[k - 1] - (Cp / grav) * (exner[k] - exner[k - 1]) * (thvm[k] - thvm[k - 1]) / math.log(
                thvm[k] / thvm[k - 1])
        else:
            z[k] = z[k - 1] - (Cp / grav) * (exner[k] - exner[k - 1]) * thvm[k]
    return z


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


def test_vs_literal():
    rng = np.random.default_rng(1)
    thvm = np.sort(290.0 + rng.uniform(0, 40, 30))      # increasing thvm with height (stable)
    exner = np.sort(rng.uniform(0.3, 0.99, 30))[::-1]   # decreasing with level
    got = np.asarray(calc_ref_z_linear_thvm(jnp.asarray(thvm), jnp.asarray(exner)))
    ref = _ref_calc_ref_z(thvm, exner)
    rel = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-30))
    assert rel < 1e-13, f"calc_ref_z vs literal rel {rel:.2e}"
    print(f"  calc_ref_z_linear_thvm vs literal NumPy (30 levels): rel {rel:.1e}  PASS")


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


def test_differentiable():
    z = np.linspace(0.0, 12000.0, 25)
    thvm = 295.0 + 0.003 * z
    exner = _exner_from_z(z, thvm, (101325.0 / p0) ** kappa)
    g = jax.grad(lambda t: jnp.sum(inverse_hydrostatic(101325.0, 0.0, t, jnp.asarray(exner)) ** 2))(jnp.asarray(thvm))
    assert np.isfinite(np.asarray(g)).all() and np.any(np.asarray(g) != 0.0), "grad bad"
    print(f"  jax.grad(inverse_hydrostatic) wrt thvm: finite, nonzero (||g||={float(jnp.linalg.norm(g)):.3e})  PASS")


def main():
    print("test_inverse_hydrostatic:")
    for t in (test_vs_literal, test_roundtrip_with_forward, test_constant_thvm_analytic, test_differentiable):
        t()
    print("All inverse_hydrostatic checks PASSED")


if __name__ == "__main__":
    main()

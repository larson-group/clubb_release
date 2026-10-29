#!/usr/bin/env python3
"""validate the JAX calc_F_x_zeta_x_setter port (new_pdf_main.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.new_pdf_main import calc_F_x_zeta_x_setter

NG, NZ = 2, 6
_SLOPE, _STDEV_FACTOR, _LAMBDA = 1.4, 1.2, 0.5


def _ref(Skx, slope, stdev_factor, lam):
    absS = np.abs(Skx)
    min_F = np.where(absS > 0, 1e-3, 0.0)
    e = np.exp(-(absS ** lam) / slope)
    F = min_F * e + 1.0 * (1 - e)
    return F, np.full_like(Skx, stdev_factor - 1.0), min_F, np.ones_like(Skx)


def test_transcription():
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(50):
        Skx = rng.uniform(-3, 3, (NG, NZ))
        g = calc_F_x_zeta_x_setter(Skx, _SLOPE, _STDEV_FACTOR, _LAMBDA)
        r = _ref(Skx, _SLOPE, _STDEV_FACTOR, _LAMBDA)
        for gi, ri in zip(g, r):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - ri)))
    assert worst < 1e-14, f"transcription mismatch {worst:.2e}"
    print(f"  calc_F_x_zeta_x_setter: literal transcription, worst {worst:.2e}  PASS")


def test_bounds_and_limits():
    Skx = np.linspace(-3, 3, 41)
    F, zeta, minF, maxF = (np.asarray(x) for x in calc_F_x_zeta_x_setter(Skx, _SLOPE, _STDEV_FACTOR, _LAMBDA))
    assert np.all(F >= minF - 1e-12) and np.all(F <= maxF + 1e-12), "F_x out of [min,max]"
    # Skx=0 -> exp factor = 1 -> F_x = min_F_x = 0.
    F0 = float(np.asarray(calc_F_x_zeta_x_setter(np.array([0.0]), _SLOPE, _STDEV_FACTOR, _LAMBDA)[0])[0])
    assert abs(F0) < 1e-14, "Skx=0 should give F_x=0"
    # Large |Skx| -> exp(-|Skx|^lambda/slope) -> 0 -> F_x -> max_F_x = 1 (needs |Skx| large enough that
    # sqrt(|Skx|)/slope >> 1, e.g. 1e4 -> exp(-100/1.4) ~ 0).
    Fbig = float(np.asarray(calc_F_x_zeta_x_setter(np.array([1.0e4]), _SLOPE, _STDEV_FACTOR, _LAMBDA)[0])[0])
    assert abs(Fbig - 1.0) < 1e-6, "large |Skx| should give F_x->1"
    assert abs(zeta[0] - (_STDEV_FACTOR - 1.0)) < 1e-14, "zeta_x = stdev_factor - 1"
    print("  bounds F_x in [min,max] + Skx=0/large limits + zeta identity  PASS")

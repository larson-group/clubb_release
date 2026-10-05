#!/usr/bin/env python3
"""validate the JAX calc_limits_F_x_responder + sort_roots ports (new_pdf.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.new_pdf import calc_limits_F_x_responder, sort_roots

NG, NZ = 2, 6
# Typical max-Skx2 thresholds (new_pdf default-parameter values).
MAX_POS, MAX_NEG = 4.0, 4.0


def _inputs(seed):
    rng = np.random.default_rng(seed)
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    Skx = rng.uniform(-3, 3, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    return mixt_frac, Skx, sgn


def test_sort_roots():
    rng = np.random.default_rng(1)
    r = rng.uniform(-5, 5, (NG, NZ, 3))
    assert np.array_equal(np.asarray(sort_roots(r)), np.sort(r, axis=-1))
    print("  sort_roots == jnp.sort (ascending)  PASS")


def test_bounds():
    mf, Skx, sgn = _inputs(5)
    g_min, g_max = (np.asarray(x) for x in
                    calc_limits_F_x_responder(mf, Skx, sgn, np.full((NG, NZ), MAX_POS), np.full((NG, NZ), MAX_NEG)))
    assert np.all(g_min >= -1e-12) and np.all(g_max <= 1.0 + 1e-12), "F_x limits out of [0,1]"
    print("  bounds: min_F_x, max_F_x in [0,1]  PASS")

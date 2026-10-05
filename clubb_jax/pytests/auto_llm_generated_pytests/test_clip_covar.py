#!/usr/bin/env python3
"""validate the JAX clip_covar port (clip_explicit.F90:clip_covar)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.clip_explicit import clip_covar

NG, NZM = 2, 9
CLIP_RTPTHLP, CLIP_WPRTP = 3, 8
MAX_MAG = 0.99


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xp2 = rng.uniform(0.1, 2.0, (NG, NZM))
    yp2 = rng.uniform(0.1, 2.0, (NG, NZM))
    # xpyp deliberately exceeds the bound at many levels to exercise both clip branches.
    xpyp = rng.uniform(-3.0, 3.0, (NG, NZM)) * np.sqrt(xp2 * yp2)
    return xp2, yp2, xpyp


def test_realizability():
    xp2, yp2, xpyp = _inputs(5)
    g, _ = clip_covar(NZM, NG, CLIP_RTPTHLP, xp2, yp2, xpyp)
    g = np.asarray(g)
    corr = g / np.sqrt(xp2 * yp2)
    # Interior levels must satisfy |corr| <= max_mag_corr (boundaries are left as-is).
    assert np.all(np.abs(corr[:, 1:-1]) <= MAX_MAG + 1e-12), "interior correlation exceeds max_mag_corr"
    # Boundaries untouched.
    assert np.allclose(g[:, 0], xpyp[:, 0]) and np.allclose(g[:, -1], xpyp[:, -1]), "boundaries changed"
    print("  realizability: |corr|<=0.99 on interior, boundaries untouched  PASS")

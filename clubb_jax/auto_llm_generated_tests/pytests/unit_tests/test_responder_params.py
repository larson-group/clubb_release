#!/usr/bin/env python3
"""validate the JAX calc_responder_params port (new_pdf.F90, Griffin & Larson 2018)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.new_pdf import calc_responder_params

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-1.5, 1.5, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    F_x = rng.uniform(0.1, 0.9, (NG, NZ))
    mixt_frac = rng.uniform(0.3, 0.7, (NG, NZ))
    return xm, xp2, Skx, sgn, F_x, mixt_frac


def test_moment_reconstruction():
    xm, xp2, Skx, sgn, F_x, mf = _inputs(5)
    mu1, mu2, s1sq, s2sq, c1, c2 = (np.asarray(x) for x in calc_responder_params(xm, xp2, Skx, sgn, F_x, mf))
    xm_rec = mf * mu1 + (1 - mf) * mu2
    assert np.max(np.abs(xm_rec - xm)) < 1e-10, "overall mean not reproduced"
    xp2_rec = mf * ((mu1 - xm) ** 2 + s1sq) + (1 - mf) * ((mu2 - xm) ** 2 + s2sq)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-10, "overall variance not reproduced"
    print("  moment reconstruction: binormal reproduces overall mean & variance (signed comp. variances)  PASS")

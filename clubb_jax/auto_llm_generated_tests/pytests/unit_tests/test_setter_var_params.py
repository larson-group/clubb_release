#!/usr/bin/env python3
"""validate the JAX calc_setter_var_params port (new_pdf.F90, Griffin & Larson 2018)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.new_pdf import calc_setter_var_params

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-2.5, 2.5, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    F_x = rng.uniform(0.05, 0.9, (NG, NZ))
    zeta_x = rng.uniform(0.0, 2.0, (NG, NZ))
    return xm, xp2, Skx, sgn, F_x, zeta_x


def test_moment_reconstruction():
    xm, xp2, Skx, sgn, F_x, zeta = _inputs(5)
    mu1, mu2, s1, s2, mf, c1, c2 = (np.asarray(x) for x in calc_setter_var_params(xm, xp2, Skx, sgn, F_x, zeta))
    # Overall mean.
    xm_rec = mf * mu1 + (1 - mf) * mu2
    assert np.max(np.abs(xm_rec - xm)) < 1e-10, "overall mean not reproduced"
    # Overall variance = a((mu1-xm)^2 + s1^2) + (1-a)((mu2-xm)^2 + s2^2) == xp2.
    xp2_rec = mf * ((mu1 - xm) ** 2 + s1 ** 2) + (1 - mf) * ((mu2 - xm) ** 2 + s2 ** 2)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-10, "overall variance not reproduced"
    # coef * xp2 == sigma^2.
    assert np.max(np.abs(c1 * xp2 - s1 ** 2)) < 1e-12 and np.max(np.abs(c2 * xp2 - s2 ** 2)) < 1e-12
    print("  moment reconstruction: binormal reproduces overall mean & variance; coef·xp2 = sigma^2  PASS")

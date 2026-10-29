#!/usr/bin/env python3
"""validate the JAX close_Luhar_pdf port (adg1_adg2_3d_luhar_pdf.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import close_Luhar_pdf

NG, NZ = 2, 6
_X_TOL_SQD = 1.0e-8


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))      # all > x_tol_sqd
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    small_m = rng.uniform(0.05, 1.0, (NG, NZ))
    wpxp = rng.uniform(-1, 1, (NG, NZ))
    return xm, xp2, mixt_frac, small_m, wpxp


def test_moment_reconstruction():
    xm, xp2, mf, m, wpxp = _inputs(5)
    ss1, ss2, v1, v2, x1n, x2n, x1, x2 = (np.asarray(x) for x in close_Luhar_pdf(xm, xp2, mf, m, wpxp, _X_TOL_SQD))
    xm_rec = mf * x1 + (1 - mf) * x2
    assert np.max(np.abs(xm_rec - xm)) < 1e-12, "overall mean not reproduced"
    xp2_rec = mf * ((x1 - xm) ** 2 + v1) + (1 - mf) * ((x2 - xm) ** 2 + v2)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-10, "overall variance not reproduced"
    print("  moment reconstruction: binormal reproduces overall mean & variance  PASS")

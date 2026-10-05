#!/usr/bin/env python3
"""validate the JAX calc_Luhar_params port (adg1_adg2_3d_luhar_pdf.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import calc_Luhar_params

NG, NZ = 2, 6
_X_TOL_SQD = 1.0e-4


def _inputs(seed):
    rng = np.random.default_rng(seed)
    Skx = rng.uniform(-3, 3, (NG, NZ))
    wpxp = rng.uniform(-1, 1, (NG, NZ))
    xp2 = rng.uniform(0.01, 2.0, (NG, NZ))
    return Skx, wpxp, xp2


def test_invariants():
    Skx, wpxp, xp2 = _inputs(5)
    mf, big_m, small_m = (np.asarray(x) for x in calc_Luhar_params(Skx, wpxp, xp2, _X_TOL_SQD))
    assert np.all(mf >= -1e-12) and np.all(mf <= 1.0 + 1e-12), "mixt_frac out of [0,1]"
    assert np.all(small_m >= 0.05 - 1e-12), "small_m below the 0.05 floor for varying x"
    # Constant-x limit.
    mf0, M0, m0 = (np.asarray(x) for x in
                   calc_Luhar_params(np.array([[1.0]]), np.array([[0.5]]), np.array([[1e-8]]), _X_TOL_SQD))
    assert abs(mf0[0, 0] - 0.5) < 1e-14 and M0[0, 0] == 0.0 and m0[0, 0] == 0.0, "constant-x limit"
    # Unskewed -> mixt_frac = 0.5.
    mfs = np.asarray(calc_Luhar_params(np.array([[0.0]]), np.array([[0.5]]), np.array([[1.0]]), _X_TOL_SQD)[0])
    assert abs(mfs[0, 0] - 0.5) < 1e-14, "Skx=0 should give mixt_frac=0.5"
    print("  invariants: mixt_frac in [0,1], small_m>=0.05, constant-x & unskewed limits  PASS")

#!/usr/bin/env python3
"""test_wp2_wp3_lhs_terms.py — validate the implicit-LHS term builders of advance_wp2_wp3_module.

Companion to test_wp2_wp3_terms.py (RHS). `advance_wp2_wp3` assembles its banded LHS matrix from per-level term
functions, validated here vs independent per-level transcriptions of the F90 loops (advance_wp2_wp3_module.F90),
interior levels with boundaries 0:

  wp2_term_dp1_lhs (F90:3604)  C1_Skw_fnc·invrs_tau_C1
  wp2_term_pr1_lhs (F90:3703)  ⅔·C4·invrs_tau_C4
  wp3_term_pr1_lhs (F90:5130)  C8·invrs_tau_wp3·(3·C8b·Skw² + 1)
  wp2_term_ta_lhs  (F90:3388)  2-band: super = +invrs_rho_ds_zm·invrs_dzm·rho_ds_zt(k); sub = −…·rho_ds_zt(k-1)

Pure-Python (the F90 formulas are the oracle), so it never SKIPs. (iter 494)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_wp2_wp3_module import (
    wp2_term_dp1_lhs, wp2_term_pr1_lhs, wp3_term_pr1_lhs, wp2_term_ta_lhs,
    wp3_term_tp_lhs, wp3_terms_ac_pr2_lhs, wp2_terms_ac_pr2_lhs,
)
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_scalar_lhs_terms():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    rng = np.random.default_rng(41)
    C1 = rng.uniform(0.5, 2.0, (_NG, nzm)); itau1 = rng.uniform(1e-3, 1e-2, (_NG, nzm))
    C4 = rng.uniform(3.0, 6.0, _NG); itau4 = rng.uniform(1e-3, 1e-2, (_NG, nzm))
    C8 = rng.uniform(3.0, 6.0, _NG); C8b = rng.uniform(0.1, 0.5, _NG)
    itau3 = rng.uniform(1e-3, 1e-2, (_NG, nzt)); Skw = rng.uniform(-3.0, 3.0, (_NG, nzt))

    jd = np.asarray(wp2_term_dp1_lhs(jnp.asarray(C1), jnp.asarray(itau1)))
    jp = np.asarray(wp2_term_pr1_lhs(jnp.asarray(C4), jnp.asarray(itau4)))
    j3 = np.asarray(wp3_term_pr1_lhs(jnp.asarray(C8), jnp.asarray(C8b), jnp.asarray(itau3), jnp.asarray(Skw)))
    rd = np.zeros((_NG, nzm)); rp = np.zeros((_NG, nzm)); r3 = np.zeros((_NG, nzt))
    for k in range(1, nzm - 1):
        rd[:, k] = C1[:, k] * itau1[:, k]
        rp[:, k] = (2.0 * C4 * itau4[:, k]) / 3.0
    for k in range(1, nzt - 1):
        r3[:, k] = C8 * itau3[:, k] * (3.0 * C8b * Skw[:, k] ** 2 + 1.0)
    for nm, j, r in (("wp2_dp1", jd, rd), ("wp2_pr1", jp, rp), ("wp3_pr1", j3, r3)):
        assert j.shape == r.shape and float(np.max(np.abs(j - r))) == 0.0, f"{nm}_lhs != F90 loop"
    print("  wp2_term_dp1_lhs / wp2_term_pr1_lhs / wp3_term_pr1_lhs == F90 (exact)  PASS")


def test_wp2_term_ta_lhs():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    invrs_dzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(42)
    invrs_rho_ds_zm = rng.uniform(0.8, 2.0, (_NG, nzm))
    rho_ds_zt = rng.uniform(0.5, 1.2, (_NG, nzt))
    j = np.asarray(wp2_term_ta_lhs(jnp.asarray(invrs_rho_ds_zm), jnp.asarray(rho_ds_zt), gr))
    # (2, ngrdcol, nzm): band 0 = super (coeff of wp3[k]); band 1 = sub (coeff of wp3[k-1])
    ref = np.zeros((2, _NG, nzm))
    for k in range(1, nzm - 1):
        fac = invrs_rho_ds_zm[:, k] * invrs_dzm[:, k]
        ref[0, :, k] = fac * rho_ds_zt[:, k]
        ref[1, :, k] = -fac * rho_ds_zt[:, k - 1]
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) == 0.0, "wp2_term_ta_lhs != F90 loop"
    print("  wp2_term_ta_lhs (2-band super/sub) == F90 (exact)  PASS")


def test_wp3_term_tp_lhs():
    """wp3_term_tp_lhs — 2-band turbulent-production LHS (F90:4885). Band 0 = coeff of wp2[k+1], band 1 = wp2[k]."""
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    invrs_dzt = np.asarray(gr.invrs_dzt)
    rng = np.random.default_rng(43)
    coef = rng.uniform(0.5, 2.0, _NG)
    wp2 = rng.uniform(1e-3, 2.0, (_NG, nzm)); rho_zm = rng.uniform(0.5, 1.2, (_NG, nzm))
    inv_rho_zt = rng.uniform(0.8, 2.0, (_NG, nzt))
    j = np.asarray(wp3_term_tp_lhs(jnp.asarray(coef), jnp.asarray(wp2), jnp.asarray(rho_zm),
                                   jnp.asarray(inv_rho_zt), gr))
    ref = np.zeros((2, _NG, nzt))
    for k in range(1, nzt - 1):
        ref[0, :, k] = coef * (-3.0 * inv_rho_zt[:, k] * invrs_dzt[:, k] * rho_zm[:, k + 1] * wp2[:, k + 1]
                               + 1.5 * invrs_dzt[:, k] * wp2[:, k + 1])
        ref[1, :, k] = coef * (3.0 * inv_rho_zt[:, k] * invrs_dzt[:, k] * rho_zm[:, k] * wp2[:, k]
                               - 1.5 * invrs_dzt[:, k] * wp2[:, k])
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) < 1e-14, "wp3_term_tp_lhs != F90 loop"
    print("  wp3_term_tp_lhs (2-band tp) == F90 (exact/FP)  PASS")


def test_ac_pr2_lhs():
    """wp3_terms_ac_pr2_lhs (F90:5011) and wp2_terms_ac_pr2_lhs (F90:3511) — accumulation+pr2 LHS via d(wm)/dz."""
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    idzt = np.asarray(gr.invrs_dzt); idzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(44)
    C11 = rng.uniform(0.3, 0.8, (_NG, nzt)); wm_zm = rng.uniform(-0.1, 0.1, (_NG, nzm))
    C_uu_shr = rng.uniform(0.1, 0.5, _NG); wm_zt = rng.uniform(-0.1, 0.1, (_NG, nzt))
    j3 = np.asarray(wp3_terms_ac_pr2_lhs(jnp.asarray(C11), jnp.asarray(wm_zm), gr))
    j2 = np.asarray(wp2_terms_ac_pr2_lhs(jnp.asarray(C_uu_shr), jnp.asarray(wm_zt), gr))
    r3 = np.zeros((_NG, nzt)); r2 = np.zeros((_NG, nzm))
    for k in range(1, nzt - 1):
        r3[:, k] = (1.0 - C11[:, k]) * 3.0 * idzt[:, k] * (wm_zm[:, k + 1] - wm_zm[:, k])
    for k in range(1, nzm - 1):
        r2[:, k] = (1.0 - C_uu_shr) * 2.0 * idzm[:, k] * (wm_zt[:, k] - wm_zt[:, k - 1])
    assert j3.shape == r3.shape and float(np.max(np.abs(j3 - r3))) == 0.0, "wp3_terms_ac_pr2_lhs != F90 loop"
    assert j2.shape == r2.shape and float(np.max(np.abs(j2 - r2))) == 0.0, "wp2_terms_ac_pr2_lhs != F90 loop"
    print("  wp3_terms_ac_pr2_lhs / wp2_terms_ac_pr2_lhs == F90 (exact)  PASS")


def main():
    print("test_wp2_wp3_lhs_terms:")
    test_scalar_lhs_terms()
    test_wp2_term_ta_lhs()
    test_wp3_term_tp_lhs()
    test_ac_pr2_lhs()
    print("All wp2_wp3 LHS term checks PASSED")


if __name__ == "__main__":
    main()

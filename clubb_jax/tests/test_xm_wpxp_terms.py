#!/usr/bin/env python3
"""test_xm_wpxp_terms.py — validate the per-level term builders of advance_xm_wpxp_module.

`advance_xm_wpxp` (the live mean-field + turbulent-flux solve run by every case, for wprtp/wpthlp/upwp/vpwp) builds
its banded LHS/RHS from per-level term functions, previously validated only implicitly by the bit-faithful suite.
This pins them vs independent per-level transcriptions of the F90 loops (advance_xm_wpxp_module.F90):

  xm_term_ta_lhs        (F90:5300) 2-band: sup = +invrs_rho_ds_zt·invrs_dzt·rho_ds_zm(k+1); sub = −…·rho_ds_zm(k)
  wpxp_term_tp_lhs      (F90:wpxp_term_tp_lhs) 2-band: sup = +wp2·invrs_dzm (xm[k]); sub = −wp2·invrs_dzm (xm[k-1])
  wpxp_terms_ac_pr2_lhs (1 − C7)·invrs_dzm·(wm_zt(k) − wm_zt(k-1))
  wpxp_term_pr1_lhs     C6·invrs_tau_C6
  wpxp_terms_bp_pr3_rhs (F90:5740) (grav/thv_ds_zm)·(1 − C7)·x'thv'   [grav = 9.81, the CLUBB-standalone value]

Pure-Python (the F90 formulas are the oracle), so it never SKIPs. (iter 501)
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

from clubb_jax.src.CLUBB_core.advance_xm_wpxp_module import (
    xm_term_ta_lhs, wpxp_term_tp_lhs, wpxp_terms_ac_pr2_lhs, wpxp_term_pr1_lhs, wpxp_terms_bp_pr3_rhs,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _grid():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]
    return gr, nzm, nzm - 1


def test_xm_term_ta_lhs():
    gr, nzm, nzt = _grid()
    idzt = np.asarray(gr.invrs_dzt)
    rng = np.random.default_rng(51)
    inv_rho_zt = rng.uniform(0.8, 2.0, (_NG, nzt)); rho_zm = rng.uniform(0.5, 1.2, (_NG, nzm))
    j = np.asarray(xm_term_ta_lhs(jnp.asarray(inv_rho_zt), jnp.asarray(rho_zm), gr))
    ref = np.zeros((2, _NG, nzt))
    for k in range(nzt):                         # no BC zeroing — full zt array
        ref[0, :, k] = inv_rho_zt[:, k] * idzt[:, k] * rho_zm[:, k + 1]
        ref[1, :, k] = -inv_rho_zt[:, k] * idzt[:, k] * rho_zm[:, k]
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) == 0.0, "xm_term_ta_lhs != F90"
    print("  xm_term_ta_lhs (2-band) == F90 (exact)  PASS")


def test_wpxp_term_tp_lhs():
    gr, nzm, nzt = _grid()
    idzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(52)
    wp2 = rng.uniform(1e-3, 2.0, (_NG, nzm))
    j = np.asarray(wpxp_term_tp_lhs(jnp.asarray(wp2), gr))
    ref = np.zeros((2, _NG, nzm))
    for k in range(1, nzm - 1):
        ref[0, :, k] = wp2[:, k] * idzm[:, k]
        ref[1, :, k] = -wp2[:, k] * idzm[:, k]
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) == 0.0, "wpxp_term_tp_lhs != F90"
    print("  wpxp_term_tp_lhs (2-band) == F90 (exact)  PASS")


def test_wpxp_ac_pr1_bp():
    gr, nzm, nzt = _grid()
    idzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(53)
    C7 = rng.uniform(0.3, 0.8, (_NG, nzm)); wm_zt = rng.uniform(-0.1, 0.1, (_NG, nzt))
    C6 = rng.uniform(2.0, 5.0, (_NG, nzm)); itau6 = rng.uniform(1e-3, 1e-2, (_NG, nzm))
    thv = rng.uniform(290.0, 320.0, (_NG, nzm)); xpthvp = rng.uniform(-0.1, 0.3, (_NG, nzm))

    jac = np.asarray(wpxp_terms_ac_pr2_lhs(jnp.asarray(C7), jnp.asarray(wm_zt), gr))
    jp1 = np.asarray(wpxp_term_pr1_lhs(jnp.asarray(C6), jnp.asarray(itau6)))
    jbp = np.asarray(wpxp_terms_bp_pr3_rhs(jnp.asarray(C7), jnp.asarray(thv), jnp.asarray(xpthvp)))
    rac = np.zeros((_NG, nzm)); rp1 = np.zeros((_NG, nzm)); rbp = np.zeros((_NG, nzm))
    for k in range(1, nzm - 1):
        rac[:, k] = (1.0 - C7[:, k]) * idzm[:, k] * (wm_zt[:, k] - wm_zt[:, k - 1])
        rp1[:, k] = C6[:, k] * itau6[:, k]
        rbp[:, k] = (grav / thv[:, k]) * (1.0 - C7[:, k]) * xpthvp[:, k]
    for nm, j, r in (("ac_pr2", jac, rac), ("pr1", jp1, rp1), ("bp_pr3", jbp, rbp)):
        assert j.shape == r.shape and float(np.max(np.abs(j - r))) == 0.0, f"wpxp_{nm} != F90"
    print("  wpxp_terms_ac_pr2_lhs / wpxp_term_pr1_lhs / wpxp_terms_bp_pr3_rhs == F90 (exact, grav=9.81)  PASS")


def main():
    print("test_xm_wpxp_terms:")
    test_xm_term_ta_lhs()
    test_wpxp_term_tp_lhs()
    test_wpxp_ac_pr1_bp()
    print("All xm_wpxp term checks PASSED")


if __name__ == "__main__":
    main()

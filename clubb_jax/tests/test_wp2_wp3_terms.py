#!/usr/bin/env python3
"""test_wp2_wp3_terms.py — validate the explicit-RHS term builders of advance_wp2_wp3_module.

`advance_wp2_wp3` (the live wp2/wp3 closure run by every case) assembles its explicit RHS from a set of per-level
term functions previously validated only implicitly by the bit-faithful suite. This pins the buoyancy/pressure/
dissipation RHS builders as first-line guards, each vs an independent per-level transcription of the F90 loop
(advance_wp2_wp3_module.F90), interior levels k=2..nz-1 with boundaries 0:

  wp2_terms_bp_pr2_rhs (F90:3794)  (1 - C_uu_buoy)·2·(g/thv_ds_zm)·w'thv'
  wp2_term_dp1_rhs                 -(C1_Skw_fnc·invrs_tau_C1)·(u'²+v'²)
  wp2_term_pr1_rhs                 ⅓·C4·(u'²+v'²)·invrs_tau_C4
  wp2_term_pr3_rhs                 ⅔·[C_uu_buoy·(g/thv)·w'thv' + C_uu_shr·(−u'w'·dum/dz − v'w'·dvm/dz)], floored 0
  wp3_terms_bp1_pr2_rhs (F90:5340) (1 - C11_Skw_fnc)·3·(g/thv_ds_zt)·w'²thv'
  wp3_term_pr1_rhs                 C8·invrs_tau_wp3·(2·C8b·Skw²)·w'³

Pure-Python (the F90 formulas are the oracle), so it never SKIPs. (iter 493)
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
    wp2_terms_bp_pr2_rhs, wp2_term_dp1_rhs, wp2_term_pr1_rhs, wp2_term_pr3_rhs,
    wp3_terms_bp1_pr2_rhs, wp3_term_pr1_rhs,
    wp2_term_pr_dfsn_rhs, wp3_term_pr_turb_rhs, wp3_term_pr_dfsn_rhs,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _grid():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    return gr, gr.zm.shape[1]


def test_wp2_bp_pr2_dp1_pr1():
    gr, nzm = _grid()
    rng = np.random.default_rng(31)
    C_uu_buoy = rng.uniform(0.1, 0.6, _NG)
    thv = rng.uniform(290.0, 320.0, (_NG, nzm)); wpthvp = rng.uniform(-0.05, 0.2, (_NG, nzm))
    C1 = rng.uniform(0.5, 2.0, (_NG, nzm)); itau1 = rng.uniform(1e-3, 1e-2, (_NG, nzm))
    up2 = rng.uniform(1e-3, 2.0, (_NG, nzm)); vp2 = rng.uniform(1e-3, 2.0, (_NG, nzm))
    C4 = rng.uniform(3.0, 6.0, _NG); itau4 = rng.uniform(1e-3, 1e-2, (_NG, nzm))

    jb = np.asarray(wp2_terms_bp_pr2_rhs(jnp.asarray(C_uu_buoy), jnp.asarray(thv), jnp.asarray(wpthvp)))
    jd = np.asarray(wp2_term_dp1_rhs(jnp.asarray(C1), jnp.asarray(itau1), jnp.asarray(up2), jnp.asarray(vp2)))
    jp = np.asarray(wp2_term_pr1_rhs(jnp.asarray(C4), jnp.asarray(up2), jnp.asarray(vp2), jnp.asarray(itau4)))
    rb = np.zeros((_NG, nzm)); rd = np.zeros((_NG, nzm)); rp = np.zeros((_NG, nzm))
    for k in range(1, nzm - 1):
        rb[:, k] = (1.0 - C_uu_buoy) * 2.0 * (grav / thv[:, k]) * wpthvp[:, k]
        rd[:, k] = -(C1[:, k] * itau1[:, k]) * (up2[:, k] + vp2[:, k])
        rp[:, k] = (C4 * (up2[:, k] + vp2[:, k]) * itau4[:, k]) / 3.0
    for nm, j, r in (("bp_pr2", jb, rb), ("dp1", jd, rd), ("pr1", jp, rp)):
        assert j.shape == r.shape and float(np.max(np.abs(j - r))) == 0.0, f"wp2_term_{nm}_rhs != F90 loop"
    print("  wp2_terms_bp_pr2_rhs / wp2_term_dp1_rhs / wp2_term_pr1_rhs == F90 (exact)  PASS")


def test_wp2_pr3():
    gr, nzm = _grid()
    nzt = nzm - 1
    invrs_dzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(32)
    C_uu_shr = rng.uniform(0.1, 0.5, _NG); C_uu_buoy = rng.uniform(0.1, 0.5, _NG)
    thv = rng.uniform(290.0, 320.0, (_NG, nzm)); wpthvp = rng.uniform(-0.05, 0.2, (_NG, nzm))
    upwp = rng.uniform(-0.5, 0.5, (_NG, nzm)); vpwp = rng.uniform(-0.5, 0.5, (_NG, nzm))
    um = rng.uniform(-15.0, 15.0, (_NG, nzt)); vm = rng.uniform(-15.0, 15.0, (_NG, nzt))
    j = np.asarray(wp2_term_pr3_rhs(jnp.asarray(C_uu_shr), jnp.asarray(C_uu_buoy), jnp.asarray(thv),
                                    jnp.asarray(wpthvp), jnp.asarray(upwp), jnp.asarray(um),
                                    jnp.asarray(vpwp), jnp.asarray(vm), gr))
    r = np.zeros((_NG, nzm))
    for k in range(1, nzm - 1):
        du = invrs_dzm[:, k] * (um[:, k] - um[:, k - 1])
        dv = invrs_dzm[:, k] * (vm[:, k] - vm[:, k - 1])
        val = (2.0 / 3.0) * (C_uu_buoy * (grav / thv[:, k]) * wpthvp[:, k]
                             + C_uu_shr * (-upwp[:, k] * du - vpwp[:, k] * dv))
        r[:, k] = np.maximum(val, 0.0)
    # ⅔·(buoy+shear) associates differently in the vectorized JAX vs this scalar loop → FP rounding ~1e-18, not 0.0
    d = float(np.max(np.abs(j - r)))
    assert j.shape == r.shape and d < 1e-14, f"wp2_term_pr3_rhs != F90 loop ({d:.2e})"
    print(f"  wp2_term_pr3_rhs == F90 (du/dz + zero-floor; FP-rounding {d:.1e})  PASS")


def test_wp3_bp1_pr2_pr1():
    gr, nzm = _grid()
    nzt = nzm - 1
    rng = np.random.default_rng(33)
    C11 = rng.uniform(0.3, 0.8, (_NG, nzt)); thv_zt = rng.uniform(290.0, 320.0, (_NG, nzt))
    wp2thvp = rng.uniform(-0.1, 0.3, (_NG, nzt))
    C8 = rng.uniform(3.0, 6.0, _NG); C8b = rng.uniform(0.1, 0.5, _NG)
    itau3 = rng.uniform(1e-3, 1e-2, (_NG, nzt)); Skw = rng.uniform(-3.0, 3.0, (_NG, nzt))
    wp3 = rng.standard_normal((_NG, nzt))
    jb = np.asarray(wp3_terms_bp1_pr2_rhs(jnp.asarray(C11), jnp.asarray(thv_zt), jnp.asarray(wp2thvp)))
    jp = np.asarray(wp3_term_pr1_rhs(jnp.asarray(C8), jnp.asarray(C8b), jnp.asarray(itau3),
                                     jnp.asarray(Skw), jnp.asarray(wp3)))
    rb = np.zeros((_NG, nzt)); rp = np.zeros((_NG, nzt))
    for k in range(1, nzt - 1):
        rb[:, k] = (1.0 - C11[:, k]) * 3.0 * (grav / thv_zt[:, k]) * wp2thvp[:, k]
        rp[:, k] = C8 * itau3[:, k] * (2.0 * C8b * Skw[:, k] ** 2) * wp3[:, k]
    for nm, j, r in (("bp1_pr2", jb, rb), ("pr1", jp, rp)):
        assert j.shape == r.shape and float(np.max(np.abs(j - r))) == 0.0, f"wp3_term_{nm}_rhs != F90 loop"
    print("  wp3_terms_bp1_pr2_rhs / wp3_term_pr1_rhs == F90 (exact)  PASS")


def test_pr_dfsn_turb():
    """The pressure-diffusion / pressure-turbulence RHS builders (deferred from iter 494, F90-verified iter 503):
      wp2_term_pr_dfsn_rhs (F90:4229)  C·invrs_rho_ds_zm·invrs_dzm·Δ(rho_ds_zt·(wpup2+wpvp2+wp3)); lower BC = first interior
      wp3_term_pr_turb_rhs (F90:5452)  −C·invrs_rho_ds_zt·invrs_dzt·Δ(rho_ds_zm·wp2·em)   [TKE branch]
      wp3_term_pr_dfsn_rhs (F90:5564)  C·invrs_rho_ds_zt·invrs_dzt·Δ(rho_ds_zm·net), net=(wp2up2+wp2vp2+wp4)−wp2·(up2+vp2+wp2)
    """
    gr = _grid()[0]; nzm = gr.zm.shape[1]; nzt = nzm - 1
    idzm = np.asarray(gr.invrs_dzm); idzt = np.asarray(gr.invrs_dzt)
    rng = np.random.default_rng(34)
    # wp2_term_pr_dfsn_rhs
    Cw2 = rng.uniform(0.1, 0.5, _NG); rho_zt = rng.uniform(0.5, 1.2, (_NG, nzt)); irho_zm = rng.uniform(0.8, 2.0, (_NG, nzm))
    wpup2 = rng.standard_normal((_NG, nzt)); wpvp2 = rng.standard_normal((_NG, nzt)); wp3 = rng.standard_normal((_NG, nzt))
    j2 = np.asarray(wp2_term_pr_dfsn_rhs(jnp.asarray(Cw2), jnp.asarray(rho_zt), jnp.asarray(irho_zm),
                                         jnp.asarray(wpup2), jnp.asarray(wpvp2), jnp.asarray(wp3), gr))
    wpuip2 = wpup2 + wpvp2 + wp3
    r2 = np.zeros((_NG, nzm))
    for m in range(1, nzm - 1):
        r2[:, m] = Cw2 * irho_zm[:, m] * idzm[:, m] * (rho_zt[:, m] * wpuip2[:, m] - rho_zt[:, m - 1] * wpuip2[:, m - 1])
    r2[:, 0] = r2[:, 1]
    # wp3_term_pr_turb_rhs (TKE branch)
    Cw3t = rng.uniform(0.1, 0.5, _NG); rho_zm = rng.uniform(0.5, 1.2, (_NG, nzm)); irho_zt = rng.uniform(0.8, 2.0, (_NG, nzt))
    wp2s = rng.uniform(1e-3, 2.0, (_NG, nzm)); ems = rng.uniform(1e-3, 2.0, (_NG, nzm))
    jt = np.asarray(wp3_term_pr_turb_rhs(jnp.asarray(Cw3t), jnp.asarray(rho_zm), jnp.asarray(irho_zt),
                                         jnp.asarray(wp2s), jnp.asarray(ems), gr))
    rt = np.zeros((_NG, nzt))
    for k in range(1, nzt - 1):
        rt[:, k] = -Cw3t * irho_zt[:, k] * idzt[:, k] * (rho_zm[:, k + 1] * wp2s[:, k + 1] * ems[:, k + 1]
                                                         - rho_zm[:, k] * wp2s[:, k] * ems[:, k])
    # wp3_term_pr_dfsn_rhs
    Cw3d = rng.uniform(0.1, 0.5, _NG)
    wp2up2 = rng.standard_normal((_NG, nzm)); wp2vp2 = rng.standard_normal((_NG, nzm)); wp4 = rng.uniform(0, 2, (_NG, nzm))
    up2 = rng.uniform(1e-3, 2, (_NG, nzm)); vp2 = rng.uniform(1e-3, 2, (_NG, nzm)); wp2 = rng.uniform(1e-3, 2, (_NG, nzm))
    jd = np.asarray(wp3_term_pr_dfsn_rhs(jnp.asarray(Cw3d), jnp.asarray(rho_zm), jnp.asarray(irho_zt),
                                         jnp.asarray(wp2up2), jnp.asarray(wp2vp2), jnp.asarray(wp4),
                                         jnp.asarray(up2), jnp.asarray(vp2), jnp.asarray(wp2), gr))
    net = (wp2up2 + wp2vp2 + wp4) - wp2 * (up2 + vp2 + wp2)
    rd = np.zeros((_NG, nzt))
    for k in range(1, nzt - 1):
        rd[:, k] = Cw3d * irho_zt[:, k] * idzt[:, k] * (rho_zm[:, k + 1] * net[:, k + 1] - rho_zm[:, k] * net[:, k])
    for nm, j, r in (("wp2_pr_dfsn", j2, r2), ("wp3_pr_turb", jt, rt), ("wp3_pr_dfsn", jd, rd)):
        d = float(np.max(np.abs(j - r)))
        assert j.shape == r.shape and d < 1e-13, f"{nm}_rhs != F90 ({d:.2e})"
    print("  wp2_term_pr_dfsn_rhs / wp3_term_pr_turb_rhs / wp3_term_pr_dfsn_rhs == F90 (exact/FP)  PASS")


def main():
    print("test_wp2_wp3_terms:")
    test_wp2_bp_pr2_dp1_pr1()
    test_wp2_pr3()
    test_wp3_bp1_pr2_pr1()
    test_pr_dfsn_turb()
    print("All wp2_wp3 term checks PASSED")


if __name__ == "__main__":
    main()

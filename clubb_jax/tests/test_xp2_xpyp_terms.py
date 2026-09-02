#!/usr/bin/env python3
"""test_xp2_xpyp_terms.py — validate the explicit-RHS term builders of advance_xp2_xpyp_module.

`advance_xp2_xpyp` (the live <x'^2>/<x'y'> closure exercised by every case) assembles its explicit RHS from a few
per-level term functions that were previously validated only *implicitly* by the bit-faithful case suite (a single
wrong coefficient mis-tunes every variance/covariance, caught only by a slow full run). This pins the three explicit
production/pressure builders as first-line guards, each vs an INDEPENDENT per-level transcription of the F90 loop
(advance_xp2_xpyp_module.F90), interior momentum levels k=2..nzm-1:

  term_tp_rhs  (F90:5622)  rhs = -w'x_b'·invrs_dzm·(x_am(k)-x_am(k-1)) - w'x_a'·invrs_dzm·(x_bm(k)-x_bm(k-1))
  term_pr1     (F90:5937)  rhs = +⅓·C4·(x_bp2+wp2)/τ_C4 - ⅓·C14·(x_bp2+wp2)/τ_C14 + C14·w_tol²/τ_C14
  term_pr2     (F90:6053)  rhs = ⅔·[ C_uu_buoy·(g/thv_ds)·w'thv' + C_uu_shr·(-u'w'·dum/dz - v'w'·dvm/dz) ], maxed 0

Pure-Python (no JAX physics run, no Fortran oracle — the formulas are the oracle), so it never SKIPs. (iter 492)
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

from clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module import term_tp_rhs, term_pr1, term_pr2, term_dp1_lhs, term_dp1_rhs
from clubb_jax.src.CLUBB_core.constants_clubb import w_tol_sqd, grav, zero_threshold
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_term_tp_rhs():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    invrs_dzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(21)
    xam = rng.standard_normal((_NG, nzt)); xbm = rng.standard_normal((_NG, nzt))
    wpxap = rng.standard_normal((_NG, nzm)); wpxbp = rng.standard_normal((_NG, nzm))
    j = np.asarray(term_tp_rhs(jnp.asarray(xam), jnp.asarray(xbm),
                               jnp.asarray(wpxap), jnp.asarray(wpxbp), jnp.asarray(invrs_dzm)))
    ref = np.zeros((_NG, nzm - 2))
    for out, m in enumerate(range(1, nzm - 1)):       # interior momentum level m (0-based)
        ref[:, out] = (-wpxbp[:, m] * invrs_dzm[:, m] * (xam[:, m] - xam[:, m - 1])
                       - wpxap[:, m] * invrs_dzm[:, m] * (xbm[:, m] - xbm[:, m - 1]))
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) == 0.0, "term_tp_rhs != F90 loop"
    print("  term_tp_rhs == F90 per-level transcription (exact)  PASS")


def test_term_pr1():
    nzm = 11
    rng = np.random.default_rng(22)
    C4 = rng.uniform(3.0, 6.0, (_NG, 1)); C14 = rng.uniform(0.5, 2.0, (_NG, 1))
    xbp2 = rng.uniform(1e-3, 2.0, (_NG, nzm)); wp2 = rng.uniform(1e-3, 2.0, (_NG, nzm))
    itau4 = rng.uniform(1e-3, 1e-2, (_NG, nzm)); itau14 = rng.uniform(1e-3, 1e-2, (_NG, nzm))
    j = np.asarray(term_pr1(jnp.asarray(C4), jnp.asarray(C14), jnp.asarray(xbp2), jnp.asarray(wp2),
                            jnp.asarray(itau4), jnp.asarray(itau14), w_tol_sqd))
    ref = np.zeros((_NG, nzm - 2))
    for out, m in enumerate(range(1, nzm - 1)):
        ref[:, out] = ((1.0 / 3.0) * C4[:, 0] * (xbp2[:, m] + wp2[:, m]) * itau4[:, m]
                       - (1.0 / 3.0) * C14[:, 0] * (xbp2[:, m] + wp2[:, m]) * itau14[:, m]
                       + C14[:, 0] * itau14[:, m] * w_tol_sqd)
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) == 0.0, "term_pr1 != F90 loop"
    print("  term_pr1 == F90 per-level transcription (exact)  PASS")


def test_term_pr2():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    invrs_dzm = np.asarray(gr.invrs_dzm)
    rng = np.random.default_rng(23)
    C_uu_shr = rng.uniform(0.1, 0.5, (_NG, 1)); C_uu_buoy = rng.uniform(0.1, 0.5, (_NG, 1))
    thv = rng.uniform(290.0, 320.0, (_NG, nzm)); wpthvp = rng.uniform(-0.05, 0.2, (_NG, nzm))
    upwp = rng.uniform(-0.5, 0.5, (_NG, nzm)); vpwp = rng.uniform(-0.5, 0.5, (_NG, nzm))
    um = rng.uniform(-15.0, 15.0, (_NG, nzt)); vm = rng.uniform(-15.0, 15.0, (_NG, nzt))
    j = np.asarray(term_pr2(jnp.asarray(C_uu_shr), jnp.asarray(C_uu_buoy), jnp.asarray(thv),
                            jnp.asarray(wpthvp), jnp.asarray(upwp), jnp.asarray(vpwp),
                            jnp.asarray(um), jnp.asarray(vm), gr))
    ref = np.zeros((_NG, nzm - 2))
    for out, m in enumerate(range(1, nzm - 1)):
        du = invrs_dzm[:, m] * (um[:, m] - um[:, m - 1])
        dv = invrs_dzm[:, m] * (vm[:, m] - vm[:, m - 1])
        val = (2.0 / 3.0) * (C_uu_buoy[:, 0] * (grav / thv[:, m]) * wpthvp[:, m]
                             + C_uu_shr[:, 0] * (-upwp[:, m] * du - vpwp[:, m] * dv))
        ref[:, out] = np.maximum(val, zero_threshold)
    assert j.shape == ref.shape and float(np.max(np.abs(j - ref))) == 0.0, "term_pr2 != F90 loop"
    # grad through the buoyancy input is finite
    g = jax.grad(lambda w: jnp.sum(term_pr2(jnp.asarray(C_uu_shr), jnp.asarray(C_uu_buoy), jnp.asarray(thv),
                w, jnp.asarray(upwp), jnp.asarray(vpwp), jnp.asarray(um), jnp.asarray(vm), gr)))(jnp.asarray(wpthvp))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad through term_pr2"
    print("  term_pr2 == F90 per-level transcription (exact) + grad finite  PASS")


def test_term_dp1():
    """term_dp1_lhs (F90:5729) main-diagonal +(Cn/tau) [interior, BC 0]; term_dp1_rhs (F90:5809) explicit
    +(Cn/tau)·threshold [ALL levels incl. boundaries]. The DP1 dissipation terms missed by iter 492."""
    nzm = 11
    rng = np.random.default_rng(24)
    Cn = rng.uniform(0.1, 2.0, (_NG, nzm)); itau = rng.uniform(1e-3, 1e-2, (_NG, nzm)); thr = 1e-8
    jl = np.asarray(term_dp1_lhs(jnp.asarray(Cn), jnp.asarray(itau)))
    jr = np.asarray(term_dp1_rhs(jnp.asarray(Cn), jnp.asarray(itau), thr))
    rl = np.zeros((_NG, nzm)); rr = np.zeros((_NG, nzm))
    for k in range(1, nzm - 1):
        rl[:, k] = Cn[:, k] * itau[:, k]
    for k in range(nzm):                         # all levels including boundaries
        rr[:, k] = Cn[:, k] * itau[:, k] * thr
    assert jl.shape == rl.shape and float(np.max(np.abs(jl - rl))) == 0.0, "term_dp1_lhs != F90 (interior, BC 0)"
    assert jr.shape == rr.shape and float(np.max(np.abs(jr - rr))) == 0.0, "term_dp1_rhs != F90 (all levels)"
    print("  term_dp1_lhs (BC 0) / term_dp1_rhs (all levels) == F90 (exact)  PASS")


def main():
    print("test_xp2_xpyp_terms:")
    test_term_tp_rhs()
    test_term_pr1()
    test_term_pr2()
    test_term_dp1()
    print("All xp2_xpyp term checks PASSED")


if __name__ == "__main__":
    main()

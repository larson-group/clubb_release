#!/usr/bin/env python3
"""test_advance_xp3_simplified.py — validate advance_xp3_module.py:advance_xp3_simplified + term_tp_rhs/term_ac_rhs.

`advance_xp3_simplified` (advance_xp3_module.F90:576, steady-state `l_predict_xp3=.false.`, `C_xp3_dissipation=1`)
diagnoses <x'^3> as `min(tau, tau_max)·(term_tp + term_ac)` per level. It is part of the NON-ADG1 prognostic-xp3
path (gated off — the iiPDF_type init-guard rejects non-ADG1, `l_advance_xp3=False` by default), so no case
exercises it. This pins it as a first-line guard:

  (1) the per-level RHS functions `term_tp_rhs` / `term_ac_rhs` vs their literal F90 formulas (exact) + grad finite;
  (2) the whole steady-state assembly vs an independent transcription built on the f2py-bit-shadowed `zt2zm`/`zm2zt`
      regrids (so the interpolation is a real Fortran oracle), using the adjacent upper level.

The Fortran source now has `kp1 = min(k+1, nzt)` (advance_xp3_module.F90), which maps each interior k to its adjacent
upper momentum level. The JAX port follows that corrected source behavior. SKIPs the
end-to-end check if clubb_f2py/clubb_python are unbuilt; the RHS-formula checks are pure-Python and never SKIP. (iter 491)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_xp3_module import (
    term_tp_rhs, term_ac_rhs, advance_xp3_simplified as _advance_xp3_simplified,
    xp3_rtp3,
)
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.constants_clubb import zero_threshold
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def advance_xp3_simplified(
    xm,
    xp2,
    wpxp,
    wpxp2,
    rho_ds_zm,
    invrs_rho_ds_zt,
    invrs_tau_zt,
    tau_max_zt,
    x_tol,
    gr,
):
    xp3, _stats = _advance_xp3_simplified(
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        gr,
        xp3_rtp3,
        60.0,
        xm,
        xp2,
        wpxp,
        wpxp2,
        rho_ds_zm,
        invrs_rho_ds_zt,
        invrs_tau_zt,
        tau_max_zt,
        x_tol,
        False,
        JaxStats.empty(
            l_sample=False,
            names=(),
            ncol=gr.ngrdcol,
            max_nlev=max(gr.nzm, gr.nzt),
        ),
        jnp.zeros((gr.ngrdcol, gr.nzt), dtype=jnp.float64),
    )
    return xp3


def test_term_rhs_formulas():
    """term_tp_rhs / term_ac_rhs vs their literal F90 formulas (advance_xp3_module.F90:934-935 / 1004) + grad."""
    rng = np.random.default_rng(11)
    shp = (_NG, 9)
    xp2_zt = rng.uniform(1e-4, 1.0, shp); wpxpp1 = rng.standard_normal(shp); wpxp = rng.standard_normal(shp)
    rho_p1 = rng.uniform(0.5, 1.2, shp); rho = rng.uniform(0.5, 1.2, shp)
    irho = rng.uniform(0.8, 2.0, shp); idzt = rng.uniform(0.01, 0.05, shp)
    xm_p1 = rng.standard_normal(shp); xm = rng.standard_normal(shp); wpxp2 = rng.standard_normal(shp)

    tp = np.asarray(term_tp_rhs(*(jnp.asarray(a) for a in
                    (xp2_zt, wpxpp1, wpxp, rho_p1, rho, irho, idzt))))
    tp_ref = 3.0 * xp2_zt * irho * idzt * (rho_p1 * wpxpp1 - rho * wpxp)
    ac = np.asarray(term_ac_rhs(*(jnp.asarray(a) for a in (xm_p1, xm, wpxp2, idzt))))
    ac_ref = -3.0 * wpxp2 * idzt * (xm_p1 - xm)
    assert float(np.max(np.abs(tp - tp_ref))) < 1e-15, "term_tp_rhs != literal F90 formula"
    assert float(np.max(np.abs(ac - ac_ref))) < 1e-15, "term_ac_rhs != literal F90 formula"
    g = jax.grad(lambda a: jnp.sum(term_tp_rhs(a, jnp.asarray(wpxpp1), jnp.asarray(wpxp),
                jnp.asarray(rho_p1), jnp.asarray(rho), jnp.asarray(irho), jnp.asarray(idzt)) ** 2))(jnp.asarray(xp2_zt))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad through term_tp_rhs"
    print("  term_tp_rhs / term_ac_rhs == literal F90 formula (exact) + grad finite  PASS")


def test_steady_state_vs_composition():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  advance_xp3_simplified composition oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    invrs_dzt = np.asarray(jgr.invrs_dzt)
    rng = np.random.default_rng(12)
    x_tol = 1e-3
    worst = 0.0
    for _ in range(10):
        xm = rng.uniform(280.0, 320.0, (ng, nzt))
        xp2 = rng.uniform(1e-4, 1.0, (ng, nzm))
        wpxp = rng.standard_normal((ng, nzm)) * 0.3
        wpxp2 = rng.standard_normal((ng, nzt)) * 0.1
        rho_ds_zm = rng.uniform(0.5, 1.2, (ng, nzm))
        invrs_rho_ds_zt = rng.uniform(0.8, 2.0, (ng, nzt))
        invrs_tau_zt = rng.uniform(1e-3, 1e-2, (ng, nzt))
        tau_max_zt = np.full((ng, nzt), 900.0)

        xp3_j = np.asarray(advance_xp3_simplified(
            jnp.asarray(xm), jnp.asarray(xp2), jnp.asarray(wpxp), jnp.asarray(wpxp2),
            jnp.asarray(rho_ds_zm), jnp.asarray(invrs_rho_ds_zt),
            jnp.asarray(invrs_tau_zt), jnp.asarray(tau_max_zt), x_tol, jgr))

        # independent transcription on f2py-bit-shadowed regrids (zt2zm with zero_threshold floor; zm2zt with x_tol^2)
        xm_zm  = np.maximum(np.asarray(clubb_f2py.f2py_zt2zm_2d(nzm, xm)), zero_threshold)   # (ng, nzm)
        xp2_zt = np.maximum(np.asarray(clubb_f2py.f2py_zm2zt_2d(nzt, xp2)), x_tol ** 2)       # (ng, nzt)
        xp3_ref = np.zeros((ng, nzt))
        for k in range(nzt - 1):
            kp1 = k + 1
            tp = 3.0 * xp2_zt[:, k] * invrs_rho_ds_zt[:, k] * invrs_dzt[:, k] * (
                rho_ds_zm[:, kp1] * wpxp[:, kp1] - rho_ds_zm[:, k] * wpxp[:, k])
            ac = -3.0 * wpxp2[:, k] * invrs_dzt[:, k] * (xm_zm[:, kp1] - xm_zm[:, k])
            xp3_ref[:, k] = np.minimum(1.0 / invrs_tau_zt[:, k], tau_max_zt[:, k]) * (tp + ac)
        worst = max(worst, float(np.max(np.abs(xp3_j - xp3_ref))))
    assert worst < 1e-12, f"advance_xp3_simplified vs transcription mismatch {worst:.2e}"
    print(f"  advance_xp3_simplified (steady-state) vs f2py-regrid transcription: bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_advance_xp3_simplified:")
    test_term_rhs_formulas()
    test_steady_state_vs_composition()
    print("All advance_xp3_simplified checks PASSED")


if __name__ == "__main__":
    main()

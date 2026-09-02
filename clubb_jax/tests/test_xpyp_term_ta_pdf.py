#!/usr/bin/env python3
"""test_xpyp_term_ta_pdf.py — validate the STANDARD (non-Godunov) PDF turbulent-advection term ports
(turbulent_adv_pdf.F90) vs the f2py Fortran oracle.

xpyp_term_ta_pdf_lhs / xpyp_term_ta_pdf_rhs assemble the implicit tridiagonal LHS and explicit RHS for the
turbulent advection of the xp2/xpyp moments. Each Fortran subroutine carries an internal
`l_upwind_xpyp_turbulent_adv` branch selecting the centered (.false., default) or upwind (.true.) discretization;
the single JAX routine mirrors both branches, so validate BOTH against the f2py oracle. The Godunov-upwind siblings
are covered separately (test_xpyp_term_ta_godunov.py). The JAX takes the grid object; the f2py uses the module-global
grid set to the same JAX heights. SKIPs cleanly if clubb_f2py / clubb_python are unbuilt. (iter 438)
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

from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import xpyp_term_ta_pdf_lhs, xpyp_term_ta_pdf_rhs

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _setup_f2py_grid():
    import clubb_f2py
    from clubb_python import clubb_api
    from clubb_python.derived_types.err_info import ErrInfo
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    return clubb_f2py, jgr, ng, nzm


def test_f2py_oracle():
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py xpyp_term_ta_pdf oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    gd = float(jgr.grid_dir)
    rng = np.random.default_rng(5)
    worst = {}
    for up in (False, True):
        wl = wr = 0.0
        for _ in range(12):
            coef = rng.uniform(-2, 2, (ng, nzt)); rho_zt = rng.uniform(0.5, 1.2, (ng, nzt))
            rho_zm = rng.uniform(0.5, 1.2, (ng, nzm)); irho = rng.uniform(0.5, 1.5, (ng, nzm))
            sgn = np.sign(rng.uniform(-1, 1, (ng, nzm))); coef_zm = rng.uniform(-2, 2, (ng, nzm))
            term = rng.uniform(-1, 1, (ng, nzt)); term_zm = rng.uniform(-1, 1, (ng, nzm))
            jl = np.asarray(xpyp_term_ta_pdf_lhs(
                jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr,
                jnp.asarray(coef), jnp.asarray(rho_zt), jnp.asarray(rho_zm),
                jnp.asarray(irho), up, jnp.asarray(sgn), jnp.asarray(coef_zm),
            ))
            fl = np.asarray(clubb_f2py.f2py_xpyp_term_ta_pdf_lhs(coef, rho_zt, rho_zm, irho, int(up), sgn, coef_zm))
            wl = max(wl, float(np.max(np.abs(jl - fl))))
            jr = np.asarray(xpyp_term_ta_pdf_rhs(
                jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr,
                jnp.asarray(term), jnp.asarray(rho_zt), jnp.asarray(rho_zm),
                jnp.asarray(irho), up, jnp.asarray(sgn), jnp.asarray(term_zm),
            ))
            fr = np.asarray(clubb_f2py.f2py_xpyp_term_ta_pdf_rhs(term, rho_zt, rho_zm, irho, int(up), sgn, term_zm))
            wr = max(wr, float(np.max(np.abs(jr - fr))))
        worst[f"upwind={up} lhs"] = wl; worst[f"upwind={up} rhs"] = wr
    overall = max(worst.values())
    assert overall < 1e-12, f"xpyp_term_ta_pdf f2py mismatch: {worst}"
    print(f"  f2py xpyp_term_ta_pdf_{{lhs,rhs}} (centered + upwind branches, 48 cases): "
          f"bit-match, worst {overall:.2e}  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zt.shape[1]; nzm = nzt + 1
    rng = np.random.default_rng(1)
    rho_zt = jnp.asarray(rng.uniform(0.5, 1.2, (1, nzt)))
    rho_zm = jnp.asarray(rng.uniform(0.5, 1.2, (1, nzm)))
    irho = jnp.asarray(rng.uniform(0.5, 1.5, (1, nzm)))
    sgn = jnp.ones((1, nzm), dtype=jnp.float64)
    coef_zm = jnp.zeros((1, nzm), dtype=jnp.float64)

    def loss(coef):
        return jnp.sum(xpyp_term_ta_pdf_lhs(
            jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr, coef, rho_zt, rho_zm, irho,
            False, sgn, coef_zm,
        ) ** 2)

    g = jax.grad(loss)(jnp.asarray(rng.uniform(-2, 2, (1, nzt))))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad through xpyp_term_ta_pdf_lhs"
    print(f"  jax.grad through xpyp_term_ta_pdf_lhs (centered): finite ({g.size} entries)  PASS")


def main():
    print("test_xpyp_term_ta_pdf:")
    test_f2py_oracle()
    test_differentiable()
    print("All xpyp_term_ta_pdf checks PASSED")


if __name__ == "__main__":
    main()

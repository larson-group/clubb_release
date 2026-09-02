#!/usr/bin/env python3
"""test_xpyp_term_ta_godunov.py — validate the JAX Godunov-upwind TA-term ports (turbulent_adv_pdf.F90).

xpyp_term_ta_pdf_lhs_godunov / xpyp_term_ta_pdf_rhs_godunov are the Godunov-like upwind vertical
discretization of the turbulent-advection term in the xp2/xpyp equations. Oracles:
  1. f2py bit-shadow vs f2py_xpyp_term_ta_pdf_{lhs,rhs}_godunov. Both use the stored grid (gr%invrs_dzm,
     gr%grid_dir); we pull the EXACT stored invrs_dzm/grid_dir via get_grid and feed the JAX functions an
     identical grid shim, so any mismatch is the discretization itself, not the grid. SKIPs if unbuilt.
  2. A boundary invariant (k=0 and k=nzm-1 rows are zero) and a finite jax.grad.
"""
import os
import sys
from collections import namedtuple

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import (
    xpyp_term_ta_pdf_lhs_godunov, xpyp_term_ta_pdf_rhs_godunov)
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_GridShim = namedtuple("_GridShim", ["invrs_dzm"])


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py xpyp_term_ta_*_godunov oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    clubb_api.init_err_info(ng)
    cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng),
                         l_implemented=False, l_ascending_grid=True, grid_type=2,
                         deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)),
                         err_info=ErrInfo(ngrdcol=ng))
    # Pull the EXACT stored grid so the JAX side uses identical invrs_dzm / grid_dir.
    g = clubb_f2py.get_grid(ng, nzm, nzt)
    invrs_dzm = np.asarray(g[4]); grid_dir = float(g[13])
    shim = _GridShim(invrs_dzm=jnp.asarray(invrs_dzm))

    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(3):
        coef = rng.uniform(-2, 2, (ng, nzt))                 # zt-level implicit coef (sign drives upwind split)
        irho = rng.uniform(0.5, 1.5, (ng, nzm))
        rho_zm = rng.uniform(0.5, 1.2, (ng, nzm))
        term_zm = rng.uniform(-1, 1, (ng, nzm))
        sgn = np.sign(rng.uniform(-1, 1, (ng, nzt)))         # ±1 turbulent-velocity sign

        f_lhs = np.asarray(clubb_f2py.f2py_xpyp_term_ta_pdf_lhs_godunov(coef, irho, rho_zm))
        j_lhs = np.asarray(xpyp_term_ta_pdf_lhs_godunov(coef, irho, rho_zm, shim, grid_dir))
        f_rhs = np.asarray(clubb_f2py.f2py_xpyp_term_ta_pdf_rhs_godunov(term_zm, irho, sgn, rho_zm))
        j_rhs = np.asarray(xpyp_term_ta_pdf_rhs_godunov(term_zm, irho, sgn, rho_zm, shim, grid_dir))
        worst = max(worst, np.max(np.abs(j_lhs - f_lhs)), np.max(np.abs(j_rhs - f_rhs)))
    assert worst < 1e-12, f"xpyp_term_ta_*_godunov f2py mismatch {worst:.2e}"
    print(f"  f2py xpyp_term_ta_pdf_{{lhs,rhs}}_godunov: bit-match, worst {worst:.2e}  PASS")


def test_boundaries_and_grad():
    ng, nzm, nzt = 2, 10, 9
    jgr = setup_grid(ngrdcol=ng, deltaz=_DZ, zm_init=0.0, zm_top=(nzm - 1) * _DZ, grid_type=1)
    rng = np.random.default_rng(2)
    coef = rng.uniform(-2, 2, (ng, nzt)); irho = rng.uniform(0.5, 1.5, (ng, nzm))
    rho_zm = rng.uniform(0.5, 1.2, (ng, nzm)); term_zm = rng.uniform(-1, 1, (ng, nzm))
    sgn = np.sign(rng.uniform(-1, 1, (ng, nzt)))
    lhs = np.asarray(xpyp_term_ta_pdf_lhs_godunov(coef, irho, rho_zm, jgr, jgr.grid_dir))
    rhs = np.asarray(xpyp_term_ta_pdf_rhs_godunov(term_zm, irho, sgn, rho_zm, jgr, jgr.grid_dir))
    assert np.all(lhs[:, :, 0] == 0) and np.all(lhs[:, :, -1] == 0), "LHS boundaries not zero"
    assert np.all(rhs[:, 0] == 0) and np.all(rhs[:, -1] == 0), "RHS boundaries not zero"

    def loss(c):
        return jnp.sum(xpyp_term_ta_pdf_lhs_godunov(c, irho, rho_zm, jgr, jgr.grid_dir) ** 2)
    grad = np.asarray(jax.grad(loss)(jnp.asarray(coef)))
    assert np.isfinite(grad).all(), "non-finite grad through lhs_godunov"
    print(f"  zero boundaries + finite jax.grad ({grad.size} entries)  PASS")


def main():
    print("test_xpyp_term_ta_godunov:")
    for t in (test_f2py_oracle, test_boundaries_and_grad):
        t()
    print("All xpyp_term_ta_*_godunov checks PASSED")


if __name__ == "__main__":
    main()

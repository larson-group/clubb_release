"""Check the STANDARD (non-Godunov) PDF turbulent-advection term ports
against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.turbulent_adv_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import xpyp_term_ta_pdf_lhs, xpyp_term_ta_pdf_rhs

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_f2py_oracle():
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
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
            fl = np.asarray(fortran.xpyp_term_ta_pdf_lhs(nzm, nzt, ng, fgr, coef, rho_zt, rho_zm, irho, up, sgn,
                                                         coef_zm))
            wl = max(wl, float(max_abs_diff("xpyp_term_ta_pdf", jl, fl)))
            jr = np.asarray(xpyp_term_ta_pdf_rhs(
                jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr,
                jnp.asarray(term), jnp.asarray(rho_zt), jnp.asarray(rho_zm),
                jnp.asarray(irho), up, jnp.asarray(sgn), jnp.asarray(term_zm),
            ))
            fr = np.asarray(fortran.xpyp_term_ta_pdf_rhs(nzm, nzt, ng, fgr, term, rho_zt, rho_zm, irho, up, sgn,
                                                         term_zm))
            wr = max(wr, float(max_abs_diff("xpyp_term_ta_pdf", jr, fr)))
        worst[f"upwind={up} lhs"] = wl; worst[f"upwind={up} rhs"] = wr
    overall = max(worst.values())
    assert overall < 1e-12, f"xpyp_term_ta_pdf f2py mismatch: {worst}"

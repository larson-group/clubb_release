"""Check the JAX Godunov-upwind TA-term ports (turbulent_adv_pdf.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.turbulent_adv_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import (
    xpyp_term_ta_pdf_lhs_godunov,
    xpyp_term_ta_pdf_rhs_godunov,
)

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_f2py_oracle():
    """Godunov-upwind turbulent-advection lhs/rhs terms match Fortran, with mixed-sign coefficients."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(3):
        coef = rng.uniform(-2, 2, (ng, nzt))                 # zt-level implicit coef (sign drives upwind split)
        irho = rng.uniform(0.5, 1.5, (ng, nzm))
        rho_zm = rng.uniform(0.5, 1.2, (ng, nzm))
        term_zm = rng.uniform(-1, 1, (ng, nzm))
        sgn = np.sign(rng.uniform(-1, 1, (ng, nzt)))         # +/-1 turbulent-velocity sign
        worst = max(worst,
                    max_abs_diff("lhs_godunov", xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ng, jgr, coef, irho, rho_zm),
                                 fortran.xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ng, fgr, coef, irho, rho_zm)),
                    max_abs_diff("rhs_godunov",
                                 xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ng, jgr, term_zm, irho, sgn, rho_zm),
                                 fortran.xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ng, fgr, term_zm, irho, sgn, rho_zm)))
    assert worst < 1e-12, f"xpyp_term_ta_pdf_*_godunov Fortran mismatch {worst:.2e}"

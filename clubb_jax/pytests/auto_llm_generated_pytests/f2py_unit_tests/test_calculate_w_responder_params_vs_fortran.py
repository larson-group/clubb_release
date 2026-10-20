"""Check the new-hybrid PDF leaf routines (new_hybrid_pdf.F90) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_hybrid_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_hybrid_pdf import calculate_w_params, calculate_responder_params

NG, NZ = 2, 8


def test_f2py_oracle():
    """Hybrid-PDF setter (calculate_w_params) and responder (calculate_responder_params) match Fortran."""
    rng = np.random.default_rng(13)
    worst_w = worst_r = 0.0
    for _ in range(4):
        sh = (NG, NZ)
        wm = rng.uniform(-1, 1, sh); wp2 = rng.uniform(0.05, 2, sh)
        Skw = rng.uniform(-2, 2, sh)
        F_w = rng.uniform(0.05, 0.95, sh); zeta_w = rng.uniform(0.0, 3.0, sh)
        fw = fortran.calculate_w_params(NZ, NG, wm, wp2, Skw, F_w, zeta_w)
        worst_w = max(worst_w, max_abs_diff("calculate_w_params",
                                            tuple(calculate_w_params(wm, wp2, Skw, F_w, zeta_w)), tuple(fw)))
        # Responder: use the setter's mixt_frac, a nonzero covariance, and a responder variable.
        mixt_frac = np.asarray(fw[4])
        xm = rng.uniform(285, 305, sh); xp2 = rng.uniform(0.05, 1, sh)
        Skx = rng.uniform(-2, 2, sh); wpxp = rng.uniform(-0.5, 0.5, sh)
        fr = fortran.calculate_responder_params(NZ, NG, xm, xp2, Skx, wpxp, wp2, F_w, mixt_frac)
        gr = calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mixt_frac)
        worst_r = max(worst_r, max_abs_diff("calculate_responder_params", tuple(gr), tuple(fr)))
    assert worst_w < 1e-12, f"calculate_w_params Fortran mismatch {worst_w:.2e}"
    assert worst_r < 1e-12, f"calculate_responder_params Fortran mismatch {worst_r:.2e}"

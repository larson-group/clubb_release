"""Check the JAX pvertinterp port (advance_helper_module.F90) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.advance_helper as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.advance_helper_module import pvertinterp

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _p_profile(nzt):
    # Pressure decreases with height (index): ~1000 hPa at the surface to ~870 hPa at the top.
    return np.tile(np.linspace(100000.0, 87000.0, nzt)[None, :], (_NG, 1)).copy()


def test_f2py_oracle():
    """pvertinterp matches Fortran for an interior target pressure and targets beyond both ends."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(4)
    p_mid = np.asfortranarray(_p_profile(nzt))
    input_var = np.asfortranarray(rng.standard_normal((ng, nzt)))
    worst = 0.0
    # Two targets lie inside the profile, so the interpolation itself is compared.
    assert np.all((p_mid[:, -1] < 88000.0) & (95000.0 < p_mid[:, 0]))
    for p_out in (95000.0, 88000.0, 101000.0, 80000.0):
        worst = max(worst, max_abs_diff(f"pvertinterp(p_out={p_out})",
                                        pvertinterp(nzt, ng, jgr, p_mid, p_out, input_var),
                                        fortran.pvertinterp(nzt, ng, fgr, p_mid.copy(), p_out, input_var.copy())))
    assert worst < 1e-11, f"pvertinterp Fortran mismatch {worst:.2e}"

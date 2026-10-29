"""Check the JAX clip_skewness port (clip_explicit.F90) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.clip_explicit as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, no_stats, setup_matching_grids
from clubb_jax.src.CLUBB_core.clip_explicit import clip_skewness

_NG, _DZ, _ZTOP = 2, 40.0, 1500.0
_SKW_MAX = 4.5


def test_f2py_oracle():
    """clip_skewness (wp3 limiter) matches Fortran with and without the smoothed-Heaviside transition."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(4)
    sfc_elevation = np.zeros(ng)
    skw_max = np.full(ng, _SKW_MAX)
    wp2_zt = rng.uniform(0.05, 2.0, (ng, nzt))
    # wp3 large enough at many levels to trigger the skewness clip.
    wp3 = rng.uniform(-8.0, 8.0, (ng, nzt)) * wp2_zt ** 1.5
    worst = 0.0
    for l_smth in (True, False):
        ref = fortran.clip_skewness(nzt, ng, fgr, 60.0, sfc_elevation, skw_max, wp2_zt, l_smth, wp3.copy())
        got, _ = clip_skewness(nzt, ng, jgr, 60.0, sfc_elevation, skw_max, wp2_zt, l_smth, no_stats(ng), wp3)
        worst = max(worst, max_abs_diff(f"clip_skewness(l_smth={l_smth})", got, ref))
    assert worst < 1e-12, f"clip_skewness Fortran mismatch {worst:.2e}"

"""Check the JAX clip_variance port (clip_explicit.F90:clip_variance)
against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.clip_explicit as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, no_stats, setup_matching_grids
from clubb_jax.src.CLUBB_core.clip_explicit import clip_variance, clip_rcm

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
CLIP_RTP2 = 1


def test_f2py_oracle():
    """clip_variance floors x'^2 at the lower threshold exactly as Fortran does."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    rng = np.random.default_rng(4)
    # xp2 with many sub-threshold (incl. negative) values to exercise the floor.
    xp2 = rng.uniform(-0.5, 2.0, (ng, nzm))
    threshold_lo = np.full((ng, nzm), 1e-3)
    ref = fortran.clip_variance(nzm, ng, fgr, CLIP_RTP2, 60.0, threshold_lo, xp2.copy())
    got, _ = clip_variance(nzm, ng, jgr, CLIP_RTP2, 60.0, threshold_lo, no_stats(ng), xp2)
    worst = max_abs_diff("clip_variance", got, ref)
    assert worst < 1e-12, f"clip_variance Fortran mismatch {worst:.2e}"


def test_clip_rcm_f2py():
    """clip_rcm (rcm <= rtm, clipped to max(zero_threshold, rtm - eps)) matches Fortran."""
    rng = np.random.default_rng(8)
    worst = 0.0
    for _ in range(20):
        ng, nzt = 2, 10
        rcm = rng.uniform(-0.1, 2e-3, (ng, nzt)); rtm = rng.uniform(0.0, 2e-3, (ng, nzt))
        worst = max(worst, max_abs_diff("clip_rcm", clip_rcm(nzt, ng, jnp.asarray(rtm), " ", jnp.asarray(rcm)),
                                        fortran.clip_rcm(nzt, ng, rtm, " ", rcm.copy())))
    assert worst < 1e-13, f"clip_rcm Fortran mismatch {worst:.2e}"

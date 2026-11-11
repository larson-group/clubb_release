"""Check the JAX clip_covar port (clip_explicit.F90:clip_covar) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.clip_explicit as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.clip_explicit import clip_covar, clip_covars_denom
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats

NG, NZM = 2, 9
CLIP_RTPTHLP, CLIP_WPRTP = 3, 8


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xp2 = rng.uniform(0.1, 2.0, (NG, NZM))
    yp2 = rng.uniform(0.1, 2.0, (NG, NZM))
    # xpyp deliberately exceeds the bound at many levels to exercise both clip branches.
    xpyp = rng.uniform(-3.0, 3.0, (NG, NZM)) * np.sqrt(xp2 * yp2)
    return xp2, yp2, xpyp


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    for solve_type in (CLIP_RTPTHLP, CLIP_WPRTP):
        xp2, yp2, xpyp = _inputs(11 + solve_type)
        f_xpyp, f_chnge = clubb_f2py.f2py_clip_covar(
            solve_type, np.asfortranarray(xp2), np.asfortranarray(yp2), np.asfortranarray(xpyp.copy()))
        g_xpyp, g_chnge = clip_covar(NZM, NG, solve_type, xp2, yp2, xpyp)
        g_xpyp, g_chnge = np.asarray(g_xpyp), np.asarray(g_chnge)
        worst = max(worst, max_abs_diff("clip_covar", g_xpyp, np.asarray(f_xpyp)),
                    max_abs_diff("clip_covar", g_chnge, np.asarray(f_chnge)))
    assert worst < 1e-12, f"clip_covar f2py mismatch {worst:.2e}"


def test_clip_covars_denom_f2py():
    """clip_covars_denom (clip_explicit.F90) — the post-wp2/wp3-solve Cauchy-Schwarz clip driver that calls
    clip_covar on wprtp(/wp2,rtp2), wpthlp(/wp2,thlp2), upwp(/wp2,up2), and vpwp(/wp2,vp2).
    Validate the four clipped covariances against the current Fortran wrapper.
    Drive both implementations with sclr_dim=0 and a padded dummy scalar extent."""
    rng = np.random.default_rng(3)
    worst = 0.0
    stats = JaxStats.empty(l_sample=False, names=(), ncol=NG, max_nlev=NZM)
    for _ in range(15):
        wp2 = rng.uniform(0.0, 1.0, (NG, NZM)); rtp2 = rng.uniform(0.0, 1e-6, (NG, NZM))
        thlp2 = rng.uniform(0.0, 1.0, (NG, NZM)); up2 = rng.uniform(0.0, 1.0, (NG, NZM)); vp2 = rng.uniform(0.0, 1.0, (NG, NZM))
        # large covars to force clipping past the ±0.99 correlation bound
        wprtp = rng.uniform(-1, 1, (NG, NZM)) * 1e-3; wpthlp = rng.uniform(-2, 2, (NG, NZM))
        upwp = rng.uniform(-2, 2, (NG, NZM)); vpwp = rng.uniform(-2, 2, (NG, NZM))
        sclrp2 = np.zeros((NG, NZM, 1)); wpsclrp = np.zeros((NG, NZM, 1)); z = np.zeros((NG, NZM))
        j = clip_covars_denom(
            NZM, NG, 0, 60.0, rtp2, thlp2, up2, vp2, wp2, sclrp2,
            False, True, stats, 0, 0, 0, 0,
            wprtp, wpthlp, upwp, vpwp, wpsclrp, z, z,
        )
        f = fortran.clip_covars_denom(
            NZM, NG, 0, 60.0, rtp2, thlp2, up2, vp2, wp2, sclrp2,
            False, True, 0, 0, 0, 0, wprtp.copy(), wpthlp.copy(), upwp.copy(), vpwp.copy(),
            wpsclrp, z.copy(), z.copy(),
        )
        for k in range(4):
            worst = max(worst, float(max_abs_diff("clip_covar", np.asarray(j[4 + k]), np.asarray(f[4 + k]))))
    assert worst < 1e-12, f"clip_covars_denom f2py mismatch {worst:.2e}"

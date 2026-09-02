#!/usr/bin/env python3
"""test_clip_covar.py — validate the JAX clip_covar port (clip_explicit.F90:clip_covar).

Clips a covariance x'y' to ±max_mag_corr·sqrt(x'^2·y'^2) on interior momentum levels (boundaries untouched), so
|corr(x,y)| <= max_mag_corr. Ported (advance_xm_wpxp_module.py) but lacking a dedicated f2py test. Oracles:
  1. f2py bit-shadow vs f2py_clip_covar — both the clipped covariance and the net change, for a non-flux
     (clip_rtpthlp) and a flux (clip_wprtp) solve_type. SKIPs if clubb_f2py is unbuilt.
  2. Realizability: |corr| <= max_mag_corr after clipping on interior levels.
  3. A finite jax.grad.
"""
import os
import sys

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

from clubb_jax.src.CLUBB_core.clip_explicit import clip_covar, clip_covars_denom
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats

NG, NZM = 2, 9
CLIP_RTPTHLP, CLIP_WPRTP = 3, 8
MAX_MAG = 0.99


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xp2 = rng.uniform(0.1, 2.0, (NG, NZM))
    yp2 = rng.uniform(0.1, 2.0, (NG, NZM))
    # xpyp deliberately exceeds the bound at many levels to exercise both clip branches.
    xpyp = rng.uniform(-3.0, 3.0, (NG, NZM)) * np.sqrt(xp2 * yp2)
    return xp2, yp2, xpyp


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py clip_covar oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for solve_type in (CLIP_RTPTHLP, CLIP_WPRTP):
        xp2, yp2, xpyp = _inputs(11 + solve_type)
        f_xpyp, f_chnge = clubb_f2py.f2py_clip_covar(
            solve_type, np.asfortranarray(xp2), np.asfortranarray(yp2), np.asfortranarray(xpyp.copy()))
        g_xpyp, g_chnge = clip_covar(NZM, NG, solve_type, xp2, yp2, xpyp)
        g_xpyp, g_chnge = np.asarray(g_xpyp), np.asarray(g_chnge)
        worst = max(worst, np.max(np.abs(g_xpyp - np.asarray(f_xpyp))),
                    np.max(np.abs(g_chnge - np.asarray(f_chnge))))
    assert worst < 1e-12, f"clip_covar f2py mismatch {worst:.2e}"
    print(f"  f2py clip_covar: bit-match (xpyp + change, non-flux/flux), worst {worst:.2e}  PASS")


def test_clip_covars_denom_f2py():
    """clip_covars_denom (clip_explicit.F90) — the post-wp2/wp3-solve Cauchy-Schwarz clip driver that calls
    clip_covar on wprtp(/wp2,rtp2), wpthlp(/wp2,thlp2), upwp(/wp2,up2), and vpwp(/wp2,vp2).
    Validate the four clipped covariances against the current Fortran wrapper.
    Drive both implementations with sclr_dim=0 and a padded dummy scalar extent.
    SKIPs if clubb_f2py is unbuilt."""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py clip_covars_denom oracle: SKIP ({type(e).__name__})")
        return
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
        f = clubb_f2py.f2py_clip_covars_denom(
            0, 60.0, rtp2, thlp2, up2, vp2, wp2, sclrp2,
            0, 1, 0, 0, 0, 0, wprtp.copy(), wpthlp.copy(), upwp.copy(), vpwp.copy(),
            wpsclrp, z.copy(), z.copy(),
        )
        for k in range(4):
            worst = max(worst, float(np.max(np.abs(np.asarray(j[4 + k]) - np.asarray(f[4 + k])))))
    assert worst < 1e-12, f"clip_covars_denom f2py mismatch {worst:.2e}"
    print(f"  f2py clip_covars_denom (4 clipped covars, current TKE formulation, 15 cases): "
          f"bit-match, worst {worst:.2e}  PASS")


def test_realizability():
    xp2, yp2, xpyp = _inputs(5)
    g, _ = clip_covar(NZM, NG, CLIP_RTPTHLP, xp2, yp2, xpyp)
    g = np.asarray(g)
    corr = g / np.sqrt(xp2 * yp2)
    # Interior levels must satisfy |corr| <= max_mag_corr (boundaries are left as-is).
    assert np.all(np.abs(corr[:, 1:-1]) <= MAX_MAG + 1e-12), "interior correlation exceeds max_mag_corr"
    # Boundaries untouched.
    assert np.allclose(g[:, 0], xpyp[:, 0]) and np.allclose(g[:, -1], xpyp[:, -1]), "boundaries changed"
    print("  realizability: |corr|<=0.99 on interior, boundaries untouched  PASS")


def test_differentiable():
    xp2, yp2, xpyp = _inputs(7)
    def loss(v):
        clipped, _ = clip_covar(NZM, NG, CLIP_RTPTHLP, xp2, yp2, v)
        return jnp.sum(clipped ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(xpyp)))
    assert np.isfinite(g).all(), "non-finite grad through clip_covar"
    print(f"  jax.grad through clip_covar: finite ({g.size} entries)  PASS")


def main():
    print("test_clip_covar:")
    for t in (test_f2py_oracle, test_clip_covars_denom_f2py, test_realizability, test_differentiable):
        t()
    print("All clip_covar checks PASSED")


if __name__ == "__main__":
    main()

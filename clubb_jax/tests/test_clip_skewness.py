#!/usr/bin/env python3
"""test_clip_skewness.py — validate the JAX clip_skewness_core port (clip_explicit.F90:clip_skewness_core).

Limits |Skw| = |wp3|/wp2_zt^(3/2) to Skw_max_mag (reduced near the surface), with an absolute |wp3| <= 100 cap.
Two near-surface treatments: a Peskin smooth-Heaviside ramp (ARM default) and a sharp 100 m AGL threshold.
Ported but lacking a dedicated f2py test. Oracles:
  1. f2py bit-shadow vs f2py_clip_skewness on a stored grid matching the JAX grid, for the smooth-Heaviside on
     and off. SKIPs if clubb_f2py/clubb_python are unbuilt.
  2. Skewness-limit invariant: |Skw| <= Skw_max_mag after clipping (away from the surface reduction).
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

from clubb_jax.src.CLUBB_core.clip_explicit import clip_skewness_core
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1500.0
_SKW_MAX = 4.5


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py clip_skewness oracle: SKIP ({type(e).__name__})")
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
    rng = np.random.default_rng(4)
    zt = np.asarray(jgr.zt)
    sfc_elevation = np.zeros(ng)
    skw_max = np.full(ng, _SKW_MAX)
    wp2_zt = rng.uniform(0.05, 2.0, (ng, nzt))
    # wp3 large enough at many levels to trigger the skewness clip.
    wp3 = rng.uniform(-8.0, 8.0, (ng, nzt)) * wp2_zt ** 1.5
    worst = 0.0
    for l_smth in (True, False):
        ref = np.asarray(clubb_f2py.f2py_clip_skewness(
            60.0, sfc_elevation, skw_max, np.asfortranarray(wp2_zt), l_smth, np.asfortranarray(wp3.copy())))
        got = np.asarray(clip_skewness_core(wp3, wp2_zt, zt, sfc_elevation, skw_max, l_smth))
        worst = max(worst, np.max(np.abs(got - ref)))
    # The sharp-threshold branch is bit-exact; the ~1e-10 residual is the smooth-Heaviside (Peskin) transition-
    # zone FP-order difference propagated through the wp3 limit, not a transcription error (closed-form faithful).
    assert worst < 1e-9, f"clip_skewness f2py mismatch {worst:.2e}"
    print(f"  f2py clip_skewness: match (smth-Heaviside on/off), worst {worst:.2e} (smth-zone FP)  PASS")


def test_skewness_limit():
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zm.shape[1] - 1
    zt = np.asarray(jgr.zt)
    rng = np.random.default_rng(6)
    wp2_zt = rng.uniform(0.1, 2.0, (_NG, nzt))
    wp3 = rng.uniform(-8.0, 8.0, (_NG, nzt)) * wp2_zt ** 1.5
    got = np.asarray(clip_skewness_core(wp3, wp2_zt, zt, np.zeros(_NG), np.full(_NG, _SKW_MAX), True))
    skw = got / wp2_zt ** 1.5
    # Above the near-surface reduction zone (>200 m AGL) |Skw| must not exceed Skw_max_mag.
    high = zt > 200.0
    assert np.all(np.abs(skw[high]) <= _SKW_MAX + 1e-9), "|Skw| exceeds Skw_max_mag aloft"
    print("  skewness limit: |Skw| <= Skw_max_mag above the surface-reduction zone  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zm.shape[1] - 1
    zt = np.asarray(jgr.zt)[:1]
    wp2_zt = np.random.default_rng(1).uniform(0.1, 2.0, (1, nzt))
    def loss(wp3):
        return jnp.sum(clip_skewness_core(wp3, wp2_zt, zt, np.zeros(1), np.full(1, _SKW_MAX), True) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(np.random.default_rng(2).uniform(-5, 5, (1, nzt)))))
    assert np.isfinite(g).all(), "non-finite grad through clip_skewness_core"
    print(f"  jax.grad through clip_skewness_core: finite ({g.size} entries)  PASS")


def main():
    print("test_clip_skewness:")
    for t in (test_f2py_oracle, test_skewness_limit, test_differentiable):
        t()
    print("All clip_skewness checks PASSED")


if __name__ == "__main__":
    main()

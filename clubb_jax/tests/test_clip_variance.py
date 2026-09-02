#!/usr/bin/env python3
"""test_clip_variance.py — validate the JAX clip_variance port (clip_explicit.F90:clip_variance).

Clamps a variance x'^2 to >= threshold_lo (and optionally <= threshold_hi) on momentum levels 0..nzm-2 (the top
boundary is left untouched). solve_type/dt only affect STATS, so the clipped values are oracle-comparable.
Oracles:
  1. f2py bit-shadow vs f2py_clip_variance on a stored grid matching the JAX grid. SKIPs if unbuilt.
  2. Floor identity + the top-boundary-untouched invariant.
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

from clubb_jax.src.CLUBB_core.clip_explicit import clip_variance, clip_rcm
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
CLIP_RTP2 = 1


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py clip_variance oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng)
    cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng),
                         l_implemented=False, l_ascending_grid=True, grid_type=2,
                         deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)),
                         err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(4)
    # xp2 with many sub-threshold (incl. negative) values to exercise the floor.
    xp2 = rng.uniform(-0.5, 2.0, (ng, nzm))
    threshold_lo = np.full((ng, nzm), 1e-3)
    ref = np.asarray(clubb_f2py.f2py_clip_variance(CLIP_RTP2, 60.0, threshold_lo, np.asfortranarray(xp2.copy())))
    got = np.asarray(clip_variance(xp2, threshold_lo))
    worst = np.max(np.abs(got - ref))
    assert worst < 1e-12, f"clip_variance f2py mismatch {worst:.2e}"
    print(f"  f2py clip_variance: bit-match, worst {worst:.2e}  PASS")


def test_floor_and_boundary():
    nzm = 9
    xp2 = np.array([[-0.1, 0.5, 1e-4, 2.0, 0.0, 1.0, -0.2, 0.3, 5.0]])
    thr = np.full((1, nzm), 1e-3)
    got = np.asarray(clip_variance(xp2, thr))
    # Levels 0..nzm-2 floored; top level nzm-1 untouched.
    expect = np.maximum(xp2, thr)
    expect[:, -1] = xp2[:, -1]
    assert np.max(np.abs(got - expect)) < 1e-14, "floor/boundary mismatch"
    print("  floor identity (levels 0..nzm-2) + top boundary untouched  PASS")


def test_differentiable():
    xp2 = jnp.asarray(np.random.default_rng(1).uniform(-0.5, 2.0, (2, 8)))
    thr = jnp.full((2, 8), 1e-3)
    g = np.asarray(jax.grad(lambda x: jnp.sum(clip_variance(x, thr) ** 2))(xp2))
    assert np.isfinite(g).all(), "non-finite grad through clip_variance"
    print(f"  jax.grad through clip_variance: finite ({g.size} entries)  PASS")


def test_clip_rcm_f2py():
    """clip_rcm (clip_explicit.F90:clip_rcm — ensure rcm <= rtm, clipping to max(zero_threshold, rtm-eps)) vs the
    Fortran oracle. Previously untested. SKIPs if clubb_f2py is unbuilt. (iter 415)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py clip_rcm oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(8)
    worst = 0.0
    for _ in range(20):
        ng, nzt = 2, 10
        rcm = rng.uniform(-0.1, 2e-3, (ng, nzt)); rtm = rng.uniform(0.0, 2e-3, (ng, nzt))
        cj = np.asarray(clip_rcm(jnp.asarray(rcm), jnp.asarray(rtm)))
        cf = np.asarray(clubb_f2py.f2py_clip_rcm(rtm, " ", rcm.copy()))
        worst = max(worst, float(np.max(np.abs(cj - cf))))
    assert worst < 1e-13, f"clip_rcm f2py mismatch {worst:.2e}"
    print(f"  f2py clip_rcm (rcm<=rtm clamp, 20 cases): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_clip_variance:")
    for t in (test_f2py_oracle, test_floor_and_boundary, test_clip_rcm_f2py, test_differentiable):
        t()
    print("All clip_variance checks PASSED")


if __name__ == "__main__":
    main()

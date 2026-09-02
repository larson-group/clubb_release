#!/usr/bin/env python3
"""test_compute_sigma_sqd_w.py — validate the JAX compute_sigma_sqd_w port (sigma_sqd_w_module.F90).

sigma_sqd_w = gamma_Skw_fnc · (1 − min(max_x corr_wx², 1)), then smoothed zm→zt→zm. Ported and used per
timestep in the gated driver but lacking a dedicated f2py test. Oracles:
  1. f2py bit-shadow vs f2py_compute_sigma_sqd_w on a stored grid matching the JAX grid, for
     l_predict_upwp_vpwp on and off. SKIPs if clubb_f2py/clubb_python are unbuilt.
  2. Physical bounds: 0 <= sigma_sqd_w <= gamma_Skw_fnc (before smoothing the factor is in [0,1]).
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

from clubb_jax.src.CLUBB_core.sigma_sqd_w_module import compute_sigma_sqd_w
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _fields(nzm, seed):
    rng = np.random.default_rng(seed)
    sh = (_NG, nzm)
    return dict(gamma_Skw_fnc=rng.uniform(0.1, 0.4, sh),
                wp2=rng.uniform(0.05, 2.0, sh), thlp2=rng.uniform(0.01, 1.0, sh),
                rtp2=rng.uniform(1e-8, 1e-6, sh), up2=rng.uniform(0.05, 2.0, sh),
                vp2=rng.uniform(0.05, 2.0, sh), wpthlp=rng.uniform(-0.5, 0.5, sh),
                wprtp=rng.uniform(-1e-4, 1e-4, sh), upwp=rng.uniform(-0.5, 0.5, sh),
                vpwp=rng.uniform(-0.5, 0.5, sh))


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py compute_sigma_sqd_w oracle: SKIP ({type(e).__name__})")
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
    f = _fields(nzm, 5)
    order = ('gamma_Skw_fnc', 'wp2', 'thlp2', 'rtp2', 'up2', 'vp2', 'wpthlp', 'wprtp', 'upwp', 'vpwp')
    worst = 0.0
    for l_pred in (False, True):
        ref = np.asarray(clubb_f2py.f2py_compute_sigma_sqd_w(
            nzt, *[np.asfortranarray(f[k]) for k in order], l_pred))
        got = np.asarray(compute_sigma_sqd_w(*[f[k] for k in order], l_pred, jgr))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"compute_sigma_sqd_w f2py mismatch {worst:.2e}"
    print(f"  f2py compute_sigma_sqd_w: bit-match (l_predict on/off), worst {worst:.2e}  PASS")


def test_bounds():
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]
    f = _fields(nzm, 7)
    order = ('gamma_Skw_fnc', 'wp2', 'thlp2', 'rtp2', 'up2', 'vp2', 'wpthlp', 'wprtp', 'upwp', 'vpwp')
    s = np.asarray(compute_sigma_sqd_w(*[f[k] for k in order], True, jgr))
    # The factor (1 - min(max_corr,1)) is in [0,1]; gamma>0; smoothing has positive weights -> 0 <= s <= max(gamma).
    assert np.all(s >= -1e-12), "sigma_sqd_w must be >= 0"
    assert np.all(s <= np.max(f['gamma_Skw_fnc']) + 1e-9), "sigma_sqd_w must be <= gamma"
    print("  bounds: 0 <= sigma_sqd_w <= gamma_Skw_fnc  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]
    global _NG
    saved = _NG; _NG = 1
    try:
        f = _fields(nzm, 9)
        order = ('gamma_Skw_fnc', 'wp2', 'thlp2', 'rtp2', 'up2', 'vp2', 'wpthlp', 'wprtp', 'upwp', 'vpwp')
        def loss(wpthlp):
            ff = dict(f); ff['wpthlp'] = wpthlp
            return jnp.sum(compute_sigma_sqd_w(*[ff[k] for k in order], True, jgr) ** 2)
        g = np.asarray(jax.grad(loss)(jnp.asarray(f['wpthlp'])))
        assert np.isfinite(g).all(), "non-finite grad through compute_sigma_sqd_w"
    finally:
        _NG = saved
    print(f"  jax.grad through compute_sigma_sqd_w: finite ({g.size} entries)  PASS")


def main():
    print("test_compute_sigma_sqd_w:")
    for t in (test_f2py_oracle, test_bounds, test_differentiable):
        t()
    print("All compute_sigma_sqd_w checks PASSED")


if __name__ == "__main__":
    main()

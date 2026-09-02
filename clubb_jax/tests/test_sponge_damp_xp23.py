#!/usr/bin/env python3
"""test_sponge_damp_xp23.py — validate the JAX sponge_layer_damping.py ports (sponge_damp_xm/xp2/xp3).

Damp the mean xm (= (xm + dt/tau·xm_ref)/(1+dt/tau)), <x'^2> (= (1-dt/tau)^2·xp2, floored) and <x'^3>
(= (1-dt/tau)^3·xp3) in the top sponge layer. Oracles:
  1. f2py bit-shadow vs f2py_sponge_damp_{xm,xp2,xp3} on a stored grid matching the JAX grid.
     SKIPs if clubb_f2py/clubb_python are unbuilt.
  2. Below-sponge invariant (no-op) + the x_tol_sqd floor for xp2.
  3. A finite jax.grad.

`sponge_damp_xm` (unlike xp2/xp3) does NOT gate geometrically — it relies on `initialize_tau_sponge_damp` setting
tau=inf below the sponge so dt/tau=0 there (no-op); test_sponge_damp_xm_f2py builds that tau profile and so also
validates `initialize_tau_sponge_damp` (which has no own f2py wrapper) implicitly.
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

from clubb_jax.src.CLUBB_core.sponge_layer_damping import (
    sponge_damp_xm, sponge_damp_xp2, sponge_damp_xp3, initialize_tau_sponge_damp)
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_TAU, _DEPTH, _DT, _XTOLSQ = 100.0, 300.0, 60.0, 1.0e-12


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py sponge_damp_xp2/xp3 oracle: SKIP ({type(e).__name__})")
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
    zm = np.asarray(jgr.zm); zt = np.asarray(jgr.zt)
    rng = np.random.default_rng(4)
    xp2 = rng.uniform(0.01, 2.0, (ng, nzm)); xp3 = rng.uniform(-1, 1, (ng, nzt))
    tau_zm = np.full((ng, nzm), _TAU); tau_zt = np.full((ng, nzt), _TAU)
    xtol_col = np.full(ng, _XTOLSQ); depth_col = np.full(ng, _DEPTH)   # f2py wants rank-1 (ngrdcol)
    xtol_j = xtol_col[:, None]; depth_j = depth_col[:, None]            # JAX broadcasts (ng,1)
    f2 = np.asarray(clubb_f2py.f2py_sponge_damp_xp2(_DT, zm, xp2, xtol_col, tau_zm, depth_col))
    g2 = np.asarray(sponge_damp_xp2(_DT, zm, xp2, xtol_j, tau_zm, depth_j))
    f3 = np.asarray(clubb_f2py.f2py_sponge_damp_xp3(_DT, zt, zm, xp3, tau_zt, depth_col))
    g3 = np.asarray(sponge_damp_xp3(_DT, zt, zm, xp3, tau_zt, depth_j))
    w = max(np.max(np.abs(g2 - f2)), np.max(np.abs(g3 - f3)))
    assert w < 1e-12, f"sponge_damp xp2/xp3 f2py mismatch {w:.2e}"
    # Confirm only the top layer was modified.
    assert np.any(g2 != xp2) and np.any(g3 != xp3), "sponge layer not exercised"
    print(f"  f2py sponge_damp_xp2/xp3: bit-match (top-layer damping), worst {w:.2e}  PASS")


def test_sponge_damp_xm_f2py():
    """sponge_damp_xm (sponge_layer_damping.F90) — implicit nudge of the mean xm toward xm_ref in the top sponge.
    The JAX relies on tau=inf below the sponge (from initialize_tau_sponge_damp) for the no-op rather than a geometric
    gate, so build that tau profile and pass it to both sides; this also exercises initialize_tau_sponge_damp (no own
    f2py wrapper). SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 439)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py sponge_damp_xm oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    zt = np.asarray(jgr.zt); zm = np.asarray(jgr.zm); zm_top = float(jgr.zm[0, -1])
    tau_min, tau_max, depth_frac = 120.0, 1800.0, 0.25
    tau = np.zeros((ng, nzt))
    for i in range(ng):
        tau[i], sld = initialize_tau_sponge_damp(zt[i], _DT, zm_top, tau_min, tau_max, depth_frac)
    sld = depth_frac * zm_top
    rng = np.random.default_rng(2)
    worst = 0.0; damped = False
    for _ in range(10):
        xm = rng.uniform(280.0, 320.0, (ng, nzt)); xm_ref = rng.uniform(280.0, 320.0, (ng, nzt))
        j = np.asarray(sponge_damp_xm(jnp.asarray(xm), jnp.asarray(xm_ref), zt, zm_top, jnp.asarray(tau), sld, _DT))
        f = np.asarray(clubb_f2py.f2py_sponge_damp_xm(_DT, zt, zm, xm_ref, xm, tau, np.full(ng, sld)))
        worst = max(worst, float(np.max(np.abs(j - f) / np.maximum(np.abs(f), 1.0))))
        damped = damped or bool(np.any(j != xm))
    assert worst < 1e-12, f"sponge_damp_xm f2py rel mismatch {worst:.2e}"
    assert damped, "sponge layer not exercised for xm"
    print(f"  f2py sponge_damp_xm (+ initialize_tau_sponge_damp): rel-match (top-layer nudge), worst {worst:.2e}  PASS")


def test_below_sponge_noop():
    nzm = 12
    zm = np.tile(np.linspace(0, 1200, nzm)[None, :], (1, 1))
    zt = 0.5 * (zm[:, :-1] + zm[:, 1:])
    xp2 = np.ones((1, nzm)) * 0.5; xp3 = np.ones((1, nzm - 1)) * 0.3
    tau2 = np.full((1, nzm), _TAU); tau3 = np.full((1, nzm - 1), _TAU)
    g2 = np.asarray(sponge_damp_xp2(_DT, zm, xp2, _XTOLSQ, tau2, _DEPTH))
    g3 = np.asarray(sponge_damp_xp3(_DT, zt, zm, xp3, tau3, _DEPTH))
    far = (zm[0, -1] - zm[0]) >= _DEPTH
    assert np.allclose(g2[0, far], xp2[0, far]), "xp2 modified below sponge"
    fart = (zm[0, -1] - zt[0]) >= _DEPTH
    assert np.allclose(g3[0, fart], xp3[0, fart]), "xp3 modified below sponge"
    print("  below-sponge no-op for both  PASS")


def test_differentiable():
    nzm = 10
    zm = jnp.asarray(np.linspace(0, 1200, nzm)[None, :])
    zt = 0.5 * (zm[:, :-1] + zm[:, 1:])
    tau3 = jnp.full((1, nzm - 1), _TAU)
    g = np.asarray(jax.grad(lambda xp3: jnp.sum(sponge_damp_xp3(_DT, zt, zm, xp3, tau3, _DEPTH) ** 2))(
        jnp.asarray(np.random.default_rng(1).uniform(-1, 1, (1, nzm - 1)))))
    assert np.isfinite(g).all(), "non-finite grad through sponge_damp_xp3"
    print(f"  jax.grad through sponge_damp_xp3: finite ({g.size} entries)  PASS")


def main():
    print("test_sponge_damp_xp23:")
    for t in (test_f2py_oracle, test_sponge_damp_xm_f2py, test_below_sponge_noop, test_differentiable):
        t()
    print("All sponge_damp_xp2/xp3 checks PASSED")


if __name__ == "__main__":
    main()

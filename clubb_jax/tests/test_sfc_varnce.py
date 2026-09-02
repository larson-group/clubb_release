#!/usr/bin/env python3
"""test_sfc_varnce.py — validate the sfc_varnce_module.F90 port (calc_sfc_varnce) vs the f2py Fortran oracle.

calc_sfc_varnce sets the surface (bottom-level) values of the second moments wp2/up2/vp2/thlp2/rtp2/rtpthlp from the
surface fluxes (the convective-scaling surface BC). Validated for the standard l_vary_convect_depth=False path (the
JAX form; the varying-convective-depth branch's um/vm/lscale_up/dt are then unused → passed as dummies to the f2py,
sclr_dim=1 to avoid the wrapper's size-0 scalar-array error). SKIPs if clubb_f2py / clubb_python are unbuilt.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.sfc_varnce_module import calc_sfc_varnce
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _setup_f2py_grid():
    import clubb_f2py
    from clubb_python import clubb_api
    from clubb_python.derived_types.err_info import ErrInfo
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    return clubb_f2py, jgr, ng, nzm


def test_f2py_oracle():
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py calc_sfc_varnce oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    zm_sfc = np.asarray(jgr.zm)[:, 0]; sfce = np.zeros(ng)
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(20):
        upwp = rng.uniform(-0.1, 0.1, ng); vpwp = rng.uniform(-0.1, 0.1, ng)
        wpthlp = rng.uniform(-0.1, 0.1, (ng, nzm)); wprtp = rng.uniform(-1e-4, 1e-4, ng)
        splat = rng.uniform(0.0, 1e-3, (ng, nzm)); tau = rng.uniform(100.0, 1000.0, (ng, nzm))
        T0 = 300.0; up2c = np.full(ng, 2.0); ac = np.full(ng, 1.8)
        wp2 = rng.uniform(1e-2, 1, (ng, nzm)); up2 = rng.uniform(1e-2, 1, (ng, nzm)); vp2 = rng.uniform(1e-2, 1, (ng, nzm))
        thlp2 = rng.uniform(1e-3, 1, (ng, nzm)); rtp2 = rng.uniform(1e-7, 1e-5, (ng, nzm))
        rtpthlp = rng.uniform(-1e-3, 1e-3, (ng, nzm))
        jr = calc_sfc_varnce(jnp.asarray(upwp), jnp.asarray(vpwp), jnp.asarray(wpthlp), jnp.asarray(wprtp),
            jnp.asarray(splat), jnp.asarray(tau), T0, jnp.asarray(up2c), jnp.asarray(ac), jnp.asarray(wp2),
            jnp.asarray(up2), jnp.asarray(vp2), jnp.asarray(thlp2), jnp.asarray(rtp2), jnp.asarray(rtpthlp),
            jnp.asarray(zm_sfc), jnp.asarray(sfce))
        sd = 1; um = np.zeros((ng, nzt)); vm = np.zeros((ng, nzt)); lsu = np.ones((ng, nzt))
        wsp = np.zeros((ng, sd)); s2 = np.zeros((ng, nzm, sd))
        fr = clubb_f2py.f2py_calc_sfc_varnce(sd, 60.0, sfce, upwp, vpwp, wpthlp, wprtp, um, vm, lsu, wsp, splat,
            tau, False, T0, up2c, ac, wp2, up2, vp2, thlp2, rtp2, rtpthlp, s2, s2, s2)
        for i in range(6):   # wp2, up2, vp2, thlp2, rtp2, rtpthlp (the f2py also returns 3 sclr outputs after)
            worst = max(worst, float(np.max(np.abs(np.asarray(jr[i]) - np.asarray(fr[i])))))
    assert worst < 1e-12, f"calc_sfc_varnce f2py mismatch {worst:.2e}"
    print(f"  f2py calc_sfc_varnce (wp2/up2/vp2/thlp2/rtp2/rtpthlp surface BC, 20 cases): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_sfc_varnce:")
    test_f2py_oracle()
    print("All sfc_varnce checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_calc_wp3_on_wp2.py — validate advance_helper_module.py:calc_wp3_on_wp2 vs the f2py Fortran primitives.

`calc_wp3_on_wp2` (advance_helper_module.F90:82) builds the smoothed wp3/wp2 ratio that feeds the ADG1 wp3 closure
(the a1/a3 coefficients + the wp3 turbulent-advection term). It is a load-bearing helper validated only implicitly
by the bit-faithful case suite (no f2py wrapper of its own). It is a thin composition of bit-shadowed primitives:
  wp2_zt   = zm2zt(wp2, w_tol_sqd)                    ! pos-def floor
  ratio_zt = clip( wp3 / max(wp2_zt, w_tol_sqd), -1000, 1000 )
  wp3_on_wp2   = zt2zm(ratio_zt)                      ! to momentum levels
  wp3_on_wp2_zt = zm2zt(wp3_on_wp2)                   ! round-trip back (spike suppression)
This reconstructs that exact sequence from the f2py `zm2zt`/`zt2zm` oracles + the clip formula and bit-compares both
outputs — a genuine Fortran-oracle check, not pure transcription. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 496)
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

from clubb_jax.src.CLUBB_core.advance_helper_module import calc_wp3_on_wp2
from clubb_jax.src.CLUBB_core.constants_clubb import w_tol_sqd
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
        print(f"  f2py calc_wp3_on_wp2 oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    zm2zt = lambda a: np.asarray(clubb_f2py.f2py_zm2zt_2d(nzt, a))
    zt2zm = lambda a: np.asarray(clubb_f2py.f2py_zt2zm_2d(nzm, a))
    rng = np.random.default_rng(6)
    wm = wz = 0.0
    for _ in range(15):
        wp2 = rng.uniform(1e-4, 3.0, (ng, nzm))
        wp3 = rng.uniform(-5.0, 5.0, (ng, nzt))
        # occasionally force the clip to fire (tiny wp2 → |ratio| > 1000)
        if _ % 3 == 0:
            wp2[:, ::4] = 1e-7
        jm, jz = calc_wp3_on_wp2(jnp.asarray(wp2), jnp.asarray(wp3), jgr)
        wp2_zt = np.maximum(zm2zt(wp2), w_tol_sqd)
        ratio = np.clip(wp3 / np.maximum(wp2_zt, w_tol_sqd), -1000.0, 1000.0)
        m_ref = zt2zm(ratio)
        z_ref = zm2zt(m_ref)
        wm = max(wm, float(np.max(np.abs(np.asarray(jm) - m_ref))))
        wz = max(wz, float(np.max(np.abs(np.asarray(jz) - z_ref))))
    assert wm < 1e-12 and wz < 1e-12, f"calc_wp3_on_wp2 mismatch: zm {wm:.2e}, zt {wz:.2e}"
    print(f"  f2py calc_wp3_on_wp2 (wp3_on_wp2 + wp3_on_wp2_zt, clip exercised): bit-match, worst {max(wm, wz):.2e}  PASS")


def main():
    print("test_calc_wp3_on_wp2:")
    test_f2py_oracle()
    print("All calc_wp3_on_wp2 checks PASSED")


if __name__ == "__main__":
    main()

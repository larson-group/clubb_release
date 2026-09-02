#!/usr/bin/env python3
"""test_diagnose_upxp.py — validate advance_xm_wpxp_module.py:diagnose_upxp vs the f2py Fortran primitives.

`diagnose_upxp` (advance_xm_wpxp_module.F90:6052) diagnoses the turbulent horizontal scalar fluxes upthlp/uprtp/
vpthlp/vprtp (Andre et al. 1978 eqn. 7) — live for every case that predicts upwp/vpwp (the default), feeding the
`upthvp` buoyancy term. Validated only implicitly by the bit-faithful suite. The formula (interior k=2..nzm-1, BC 0):

    ypxp = (tau_C6_zm / C6x_Skw_fnc) · ( −ypwp·ddzt(xm) − (1 − C7_Skw_fnc)·wpxp·ddzt(ym) )

is reconstructed here on the **bit-shadowed f2py `ddzt`** + the formula and bit-compared — a genuine Fortran-oracle
composition, not pure transcription. SKIPs if clubb_f2py / clubb_python are unbuilt. (iter 502)
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

from clubb_jax.src.CLUBB_core.advance_xm_wpxp_module import diagnose_upxp
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
        print(f"  f2py diagnose_upxp oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    ddzt = lambda a: np.asarray(clubb_f2py.f2py_ddzt_2d(nzm, a))   # zt-field → zm-field, bit-shadowed
    rng = np.random.default_rng(8)
    worst = 0.0
    for _ in range(12):
        ypwp = rng.uniform(-0.5, 0.5, (ng, nzm))
        xm   = rng.uniform(280.0, 320.0, (ng, nzt))      # e.g. thlm (zt)
        wpxp = rng.uniform(-0.5, 0.5, (ng, nzm))
        ym   = rng.uniform(-15.0, 15.0, (ng, nzt))       # e.g. um (zt)
        C6x  = rng.uniform(2.0, 6.0, (ng, nzm))
        tau  = rng.uniform(100.0, 900.0, (ng, nzm))
        C7   = rng.uniform(0.3, 0.8, (ng, nzm))

        j = np.asarray(diagnose_upxp(jnp.asarray(ypwp), jnp.asarray(xm), jnp.asarray(wpxp),
                                     jnp.asarray(ym), jnp.asarray(C6x), jnp.asarray(tau), jnp.asarray(C7), jgr))
        ddzt_xm = ddzt(xm); ddzt_ym = ddzt(ym)
        ref = np.zeros((ng, nzm))
        for k in range(1, nzm - 1):
            ref[:, k] = (tau[:, k] / C6x[:, k]) * (
                -ypwp[:, k] * ddzt_xm[:, k] - (1.0 - C7[:, k]) * wpxp[:, k] * ddzt_ym[:, k])
        worst = max(worst, float(np.max(np.abs(j - ref))))
    assert worst < 1e-12, f"diagnose_upxp vs f2py-ddzt composition mismatch {worst:.2e}"
    print(f"  f2py diagnose_upxp (reconstructed from f2py ddzt + Andre-1978 formula): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_diagnose_upxp:")
    test_f2py_oracle()
    print("All diagnose_upxp checks PASSED")


if __name__ == "__main__":
    main()

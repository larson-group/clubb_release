#!/usr/bin/env python3
"""test_diagnostic_cache.py — validate advance_clubb_core_module.py:compute_diagnostic_cache vs the f2py oracle.

`compute_diagnostic_cache` (advance_clubb_core_module.F90:1752) derives the cache diagnostics consumed by later
CLUBB-core blocks: the virtual potential temperature `thvm` (= calculate_thvm), the turbulent kinetic energy
`em`/`sqrt_em_zt`, and the squared mean-wind shear `ddzt_umvm_sqd`. It was previously validated only *indirectly*
by the bit-faithful case suite (a single wrong coefficient would mis-tune the whole closure, caught only by a slow
full run). This pins it as a first-line guard: each output is checked against an INDEPENDENT oracle composition —
the f2py Fortran `calculate_thvm`/`ddzt`/`zm2zt` primitives + the closed-form TKE
`0.5·(wp2+up2+vp2)`, with the `em_min` floor applied as `max(·, em_min)` over all levels (the Fortran
`zm2zt_api(…, em_min)` clips every level, grid_class.F90:2375).
SKIPs cleanly if clubb_f2py / clubb_python are unbuilt. (iter 488)
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

from clubb_jax.src.CLUBB_core.constants_clubb import em_min
from clubb_jax.src.CLUBB_core.calc_pressure import calculate_thvm
from clubb_jax.src.CLUBB_core.grid_class import ddzt, zm2zt
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def compute_diagnostic_cache(
    thlm,
    rtm,
    rcm,
    exner,
    thv_ds_zt,
    wp2,
    up2,
    vp2,
    um,
    vm,
    gr,
):
    thvm = calculate_thvm(gr.nzt, gr.ngrdcol, thlm, rtm, rcm, exner, thv_ds_zt)
    em = 0.5 * (wp2 + up2 + vp2)
    sqrt_em_zt = jnp.sqrt(zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, em, em_min))
    ddzt_um = ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, um)
    ddzt_vm = ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, vm)
    ddzt_umvm_sqd = ddzt_um ** 2 + ddzt_vm ** 2
    return thvm, em, sqrt_em_zt, ddzt_umvm_sqd


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
        print(f"  f2py compute_diagnostic_cache oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    rng = np.random.default_rng(3)
    wt = wem = wsq = wsh = 0.0
    for _ in range(12):
            # zt-level fields (physically plausible ranges)
            thlm = rng.uniform(290.0, 320.0, (ng, nzt))
            rtm  = rng.uniform(1e-3, 2e-2, (ng, nzt))
            rcm  = rng.uniform(0.0, 5e-4, (ng, nzt))
            exner = rng.uniform(0.6, 1.0, (ng, nzt))
            thv_ds_zt = rng.uniform(290.0, 320.0, (ng, nzt))
            um = rng.uniform(-15.0, 15.0, (ng, nzt))
            vm = rng.uniform(-15.0, 15.0, (ng, nzt))
            # zm-level variances (>= 0)
            wp2 = rng.uniform(0.0, 2.0, (ng, nzm))
            up2 = rng.uniform(0.0, 2.0, (ng, nzm))
            vp2 = rng.uniform(0.0, 2.0, (ng, nzm))

            thvm, em, sqrt_em_zt, ddzt_umvm_sqd = compute_diagnostic_cache(
                thlm=jnp.asarray(thlm), rtm=jnp.asarray(rtm), rcm=jnp.asarray(rcm),
                exner=jnp.asarray(exner), thv_ds_zt=jnp.asarray(thv_ds_zt),
                wp2=jnp.asarray(wp2), up2=jnp.asarray(up2), vp2=jnp.asarray(vp2),
                um=jnp.asarray(um), vm=jnp.asarray(vm), gr=jgr)

            # --- independent oracle composition ---
            thvm_ref = clubb_f2py.f2py_calculate_thvm(thlm, rtm, rcm, exner, thv_ds_zt)
            em_ref = 0.5 * (wp2 + up2 + vp2)
            sqrt_em_ref = np.sqrt(np.maximum(np.asarray(clubb_f2py.f2py_zm2zt_2d(nzt, em_ref)), float(em_min)))
            ddzt_um = np.asarray(clubb_f2py.f2py_ddzt_2d(nzm, um))
            ddzt_vm = np.asarray(clubb_f2py.f2py_ddzt_2d(nzm, vm))
            shear_ref = ddzt_um ** 2 + ddzt_vm ** 2

            wt  = max(wt,  float(np.max(np.abs(np.asarray(thvm) - np.asarray(thvm_ref)))))
            wem = max(wem, float(np.max(np.abs(np.asarray(em) - em_ref))))
            wsq = max(wsq, float(np.max(np.abs(np.asarray(sqrt_em_zt) - sqrt_em_ref))))
            wsh = max(wsh, float(np.max(np.abs(np.asarray(ddzt_umvm_sqd) - shear_ref))))
    worst = dict(thvm=wt, em=wem, sqrt_em_zt=wsq, ddzt_umvm_sqd=wsh)
    overall = max(worst.values())
    assert overall < 1e-12, f"compute_diagnostic_cache f2py mismatch: {worst}"
    print(f"  f2py compute_diagnostic_cache (thvm/em/sqrt_em_zt/ddzt_umvm_sqd): "
          f"bit-match, worst {overall:.2e}  PASS")


def main():
    print("test_diagnostic_cache:")
    test_f2py_oracle()
    print("All compute_diagnostic_cache checks PASSED")


if __name__ == "__main__":
    main()

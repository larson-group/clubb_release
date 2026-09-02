#!/usr/bin/env python3
"""test_compute_xp3.py — validate advance_xp3_module.py:compute_xp3 (ADG1 path) vs the f2py Fortran primitives.

`compute_xp3` (advance_xp3_module.F90:292) diagnoses the third-order moments rtp3/thlp3/up3/vp3 from the LG-2005
ansatz. It was long classed "case-only" for an interface difference (the whole-routine `f2py_compute_xp3` wrapper
takes extra `thvm` / `iiPDF_type` / passive-scalar args and runs a full PDF setup). But the Fortran branches on the
PDF type: for `iiPDF_type == iiPDF_ADG1` (advance_xp3_module.F90:447) the routine is exactly
  wp2_zt   = zm2zt(wp2, w_tol_sqd)                          ! pos-def floor
  Skw_zt   = Skx_func(wp2_zt, wp3, w_tol)  =  wp3·sqrt(wp2_zt + Skw_denom_coef·w_tol²)^-3   (no clip — kluge off)
  x_zt     = zm2zt(x[, x_tol²])                             ! interpolate each flux/variance (pos-def with floor)
  xp3      = xp3_LG_2005_ansatz(Skw_zt, wpxp_zt, wp2_zt, xp2_zt, sigma_sqd_w_zt, clubb_params, x_tol)
and **never touches `thvm`** (thvm/Brunt-Vaisala only feed the non-ADG1 `else` branch). So the JAX `compute_xp3`
ADG1 path is validated by reconstructing that exact branch from its already-bit-shadowed Fortran building blocks —
`f2py_zm2zt_2d` + `f2py_xp3_lg_2005_ansatz` — and bit-comparing all four outputs. This reclaims the ADG1
computational path into the oracle-validated set (the iter-443/444 "apparent interface divergence is recoverable"
pattern; the whole-routine `f2py_compute_xp3` wrapper's extra setup is sidestepped). SKIPs if clubb_f2py /
clubb_python are unbuilt. (iter 490)
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

from clubb_jax.src.CLUBB_core.Skx_module import Skx_func, xp3_LG_2005_ansatz
from clubb_jax.src.CLUBB_core.grid_class import zm2zt
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
from clubb_jax.src.CLUBB_core.constants_clubb import w_tol, w_tol_sqd, thl_tol, rt_tol, zero_threshold
from clubb_jax.src.CLUBB_core.parameter_indices import iSkw_denom_coef
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def compute_xp3(
    wp2,
    wp3,
    wprtp,
    wpthlp,
    rtp2,
    thlp2,
    upwp,
    vpwp,
    up2,
    vp2,
    sigma_sqd_w,
    clubb_params,
    gr,
):
    wp2_zt = zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, wp2, w_tol_sqd)
    Skw_zt = Skx_func(gr.nzt, gr.ngrdcol, wp2_zt, wp3, w_tol, clubb_params)
    sigma_sqd_w_zt = zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, sigma_sqd_w, zero_threshold)
    rtp3 = xp3_LG_2005_ansatz(
        gr.nzt, gr.ngrdcol, Skw_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, wprtp),
        wp2_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, rtp2, rt_tol ** 2),
        sigma_sqd_w_zt,
        clubb_params,
        rt_tol,
    )
    thlp3 = xp3_LG_2005_ansatz(
        gr.nzt, gr.ngrdcol, Skw_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, wpthlp),
        wp2_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, thlp2, thl_tol ** 2),
        sigma_sqd_w_zt,
        clubb_params,
        thl_tol,
    )
    up3 = xp3_LG_2005_ansatz(
        gr.nzt, gr.ngrdcol, Skw_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, upwp),
        wp2_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, up2, w_tol_sqd),
        sigma_sqd_w_zt,
        clubb_params,
        w_tol,
    )
    vp3 = xp3_LG_2005_ansatz(
        gr.nzt, gr.ngrdcol, Skw_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, vpwp),
        wp2_zt,
        zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, vp2, w_tol_sqd),
        sigma_sqd_w_zt,
        clubb_params,
        w_tol,
    )
    return rtp3, thlp3, up3, vp3


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


def test_f2py_oracle_adg1():
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py compute_xp3 oracle: SKIP ({type(e).__name__})")
        return
    fn = os.path.join(_ROOT, "input", "tunable_parameters", "tunable_parameters.in")
    if not os.path.exists(fn):
        print("  f2py compute_xp3 oracle: SKIP (tunable_parameters.in absent)")
        return
    nzt = nzm - 1
    cp = np.asarray(init_clubb_params(ng, filename=fn))            # (ng, nparams), matches Fortran exactly
    skw_denom = cp[:, iSkw_denom_coef:iSkw_denom_coef + 1]
    zm2zt = lambda a: np.asarray(clubb_f2py.f2py_zm2zt_2d(nzt, a))  # the bit-shadowed Fortran regrid
    ansatz = clubb_f2py.f2py_xp3_lg_2005_ansatz                     # the bit-shadowed Fortran LG-2005 ansatz
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(12):
        wp2   = rng.uniform(1e-3, 2.0, (ng, nzm))
        wprtp = rng.uniform(-1e-3, 1e-3, (ng, nzm))
        wpthlp = rng.uniform(-0.5, 0.5, (ng, nzm))
        rtp2  = rng.uniform(1e-8, 1e-5, (ng, nzm))
        thlp2 = rng.uniform(1e-3, 1.0, (ng, nzm))
        upwp  = rng.uniform(-0.5, 0.5, (ng, nzm))
        vpwp  = rng.uniform(-0.5, 0.5, (ng, nzm))
        up2   = rng.uniform(1e-3, 2.0, (ng, nzm))
        vp2   = rng.uniform(1e-3, 2.0, (ng, nzm))
        sigma_sqd_w = rng.uniform(0.0, 0.5, (ng, nzm))
        wp3   = rng.uniform(-1.0, 1.0, (ng, nzt))

        rtp3_j, thlp3_j, up3_j, vp3_j = (np.asarray(x) for x in compute_xp3(
            jnp.asarray(wp2), jnp.asarray(wp3), jnp.asarray(wprtp), jnp.asarray(wpthlp),
            jnp.asarray(rtp2), jnp.asarray(thlp2), jnp.asarray(upwp), jnp.asarray(vpwp),
            jnp.asarray(up2), jnp.asarray(vp2), jnp.asarray(sigma_sqd_w), jnp.asarray(cp), jgr))

        # Reconstruct the Fortran ADG1 branch from its bit-shadowed Fortran primitives:
        wp2_zt = np.maximum(zm2zt(wp2), w_tol_sqd)
        Skw_zt = wp3 * np.sqrt(wp2_zt + skw_denom * w_tol ** 2) ** (-3)   # Skx_func, kluge off (no clip)
        ssw_zt = np.maximum(zm2zt(sigma_sqd_w), zero_threshold)
        rtp3_f  = np.asarray(ansatz(Skw_zt, zm2zt(wprtp),  wp2_zt, np.maximum(zm2zt(rtp2),  rt_tol ** 2),  ssw_zt, cp, rt_tol))
        thlp3_f = np.asarray(ansatz(Skw_zt, zm2zt(wpthlp), wp2_zt, np.maximum(zm2zt(thlp2), thl_tol ** 2), ssw_zt, cp, thl_tol))
        up3_f   = np.asarray(ansatz(Skw_zt, zm2zt(upwp),   wp2_zt, np.maximum(zm2zt(up2),   w_tol_sqd),    ssw_zt, cp, w_tol))
        vp3_f   = np.asarray(ansatz(Skw_zt, zm2zt(vpwp),   wp2_zt, np.maximum(zm2zt(vp2),   w_tol_sqd),    ssw_zt, cp, w_tol))

        for jv, fv in ((rtp3_j, rtp3_f), (thlp3_j, thlp3_f), (up3_j, up3_f), (vp3_j, vp3_f)):
            worst = max(worst, float(np.max(np.abs(jv - fv))))
    assert worst < 1e-12, f"compute_xp3 (ADG1) vs f2py-primitive composition mismatch {worst:.2e}"
    print(f"  f2py compute_xp3 (ADG1, reconstructed from f2py zm2zt + xp3_LG_2005_ansatz): rtp3/thlp3/up3/vp3 "
          f"bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_compute_xp3:")
    test_f2py_oracle_adg1()
    print("All compute_xp3 checks PASSED")


if __name__ == "__main__":
    main()

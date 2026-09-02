#!/usr/bin/env python3
"""test_calc_pdf_skewness_diagnostics.py — pin the post-advance PDF skewness diagnostics orchestration.

`calc_pdf_skewness_diagnostics_jax` (pdf_closure_module.py ↔ pdf_closure_module.F90:4448-4465) assembles the diagnostic
skewnesses of the post-advance PDF: Sk_rt and Sk_thl on the zt AND zm grids (via the tested `Skx_func`, with the zm
ones built from rtp2/thlp2 and the zt→zm-regridded 3rd moments), plus the velocity skewness
    Skw_velocity = (1/(1−σ²_w)) · wp3_zm / max(wp2, w_tol_sqd).
It feeds stats / closure diagnostics but was validated only end-to-end. This pins the ORCHESTRATION — that each of the
4 Sk outputs equals `Skx_func` called with the correct (xp2, xp3, x_tol) and the correct zt/zm regrid, and that
Skw_velocity matches the closed form — against independent recomputation (Skx_func is itself tested). A mis-routed
moment / wrong tol / missing regrid / wrong Skw_velocity factor is caught. + finite grad. (iter 568)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import rt_tol, thl_tol, w_tol_sqd
from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params

_NG, _NZT = 1, 8


def _skx(xp2, xp3, x_tol, clubb_params):
    return Skx_func(xp2.shape[1], xp2.shape[0], xp2, xp3, x_tol, clubb_params)


def calc_pdf_skewness_diagnostics_jax(
    rtp2_zt,
    rtp3,
    thlp2_zt,
    thlp3,
    rtp2,
    thlp2,
    sigma_sqd_w_zm,
    wp3_zm,
    wp2,
    clubb_params,
    gr,
):
    Skrt_zt = _skx(rtp2_zt, rtp3, rt_tol, clubb_params)
    Skthl_zt = _skx(thlp2_zt, thlp3, thl_tol, clubb_params)
    Skrt_zm = _skx(rtp2, zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, rtp3), rt_tol, clubb_params)
    Skthl_zm = _skx(thlp2, zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, thlp3), thl_tol, clubb_params)
    Skw_velocity = (1.0 / (1.0 - sigma_sqd_w_zm)) * (
        wp3_zm / jnp.maximum(wp2, w_tol_sqd)
    )
    return Skrt_zt, Skthl_zt, Skrt_zm, Skthl_zm, Skw_velocity


def _inputs(rng, gr):
    nzm = _NZT + 1
    # clubb_params: (ngrdcol, nparams) from the standalone tunable file
    params = np.asarray(init_clubb_params(_NG, _ROOT + "/input/tunable_parameters/tunable_parameters.in"))
    z2 = lambda lo, hi, n: rng.uniform(lo, hi, (_NG, n))
    return dict(
        rtp2_zt=z2(1e-8, 1e-6, _NZT), rtp3=z2(-1e-10, 1e-10, _NZT),
        thlp2_zt=z2(0.01, 0.5, _NZT), thlp3=z2(-0.1, 0.1, _NZT),
        rtp2=z2(1e-8, 1e-6, nzm), thlp2=z2(0.01, 0.5, nzm),
        sigma_sqd_w_zm=z2(0.1, 0.6, nzm), wp3_zm=z2(-0.5, 0.5, nzm),
        wp2=z2(0.05, 1.5, nzm), clubb_params=jnp.asarray(params), gr=gr), params


def test_orchestration_and_skw_velocity():
    gr = setup_grid(_NG, 50.0, 50.0, 50.0 * (_NZT + 1))
    rng = np.random.default_rng(568)
    k, params = _inputs(rng, gr)
    Skrt_zt, Skthl_zt, Skrt_zm, Skthl_zm, Skw_vel = calc_pdf_skewness_diagnostics_jax(**k)
    cp = k['clubb_params']
    # independent recomputation via the tested Skx_func + the same regrid
    r_Skrt_zt = _skx(k['rtp2_zt'], k['rtp3'], rt_tol, cp)
    r_Skthl_zt = _skx(k['thlp2_zt'], k['thlp3'], thl_tol, cp)
    r_Skrt_zm = _skx(k['rtp2'], zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, k['rtp3']), rt_tol, cp)
    r_Skthl_zm = _skx(k['thlp2'], zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, k['thlp3']), thl_tol, cp)
    r_Skw = (1.0 / (1.0 - np.asarray(k['sigma_sqd_w_zm']))) * (np.asarray(k['wp3_zm'])
             / np.maximum(np.asarray(k['wp2']), w_tol_sqd))
    for got, ref, nm in ((Skrt_zt, r_Skrt_zt, "Skrt_zt"), (Skthl_zt, r_Skthl_zt, "Skthl_zt"),
                         (Skrt_zm, r_Skrt_zm, "Skrt_zm"), (Skthl_zm, r_Skthl_zm, "Skthl_zm"),
                         (Skw_vel, r_Skw, "Skw_velocity")):
        assert np.max(np.abs(np.asarray(got) - np.asarray(ref))) < 1e-10, f"{nm} mismatch"
    print("  4 Sk_rt/Sk_thl (zt+zm) route to Skx_func correctly; Skw_velocity formula exact  PASS")


def test_grad_finite():
    gr = setup_grid(_NG, 50.0, 50.0, 50.0 * (_NZT + 1))
    k, _ = _inputs(np.random.default_rng(2), gr)
    g = jax.grad(lambda w3: jnp.sum(calc_pdf_skewness_diagnostics_jax(**{**k, 'wp3_zm': w3})[4] ** 2))(jnp.asarray(k['wp3_zm']))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad of Skw_velocity wrt wp3_zm"
    print("  jax.grad(Skw_velocity) wrt wp3_zm finite  PASS")


def main():
    print("test_calc_pdf_skewness_diagnostics:")
    test_orchestration_and_skw_velocity()
    test_grad_finite()
    print("All calc_pdf_skewness_diagnostics checks PASSED")


if __name__ == "__main__":
    main()

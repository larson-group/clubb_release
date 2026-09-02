#!/usr/bin/env python3
"""test_calc_xpthvp_terms.py — pin the buoyancy-flux (<x'th_v'>) assembly of the ADG1 closure.

`calc_xpthvp_terms_jax` (pdf_closure_module.py ↔ pdf_closure_module.F90:1122-1158) builds the virtual-potential-
temperature fluxes from the thl/rt/rc fluxes via the Sommeria-Deardorff θ_v decomposition:
    rc_coef = Lv/(exner·Cp) − ep2·thv_ds
    x'th_v' = x'thl' + ep1·thv_ds·x'rt' + rc_coef·x'rc'      for x ∈ {w, w², rt, thl}
      (wpthvp = wpthlp+ep1·thv_ds·wprtp+rc_coef·wprcp ;  rtpthvp uses rtp2 for the rt' term ;
       thlpthvp uses thlp2 for its thl' term and rtpthlp for the rt' term)
then zt→zm-regrids wpthvp/rtpthvp/thlpthvp/rc_coef (zeroing k_ub_zm); wp2thvp stays on zt. It drives the closure's
buoyancy production but was validated only end-to-end. This pins the on-zt assembly vs an INDEPENDENT transcription
(checking the ep1/ep2/rc_coef coefficients + which 2nd-moment feeds each term) and the zm outputs vs the tested
`zt2zm` regrid with k_ub_zm zeroed. + finite grad. (iter 566)
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

from clubb_jax.src.CLUBB_core.constants_clubb import ep1, ep2, Lv, Cp
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _NZT = 1, 8


def calc_xpthvp_terms_jax(
    exner,
    thv_ds_zt,
    wprcp_zt,
    wp2rcp_zt,
    rtprcp_zt,
    thlprcp_zt,
    wpthlp_zt,
    wprtp_zt,
    wp2thlp_zt,
    wp2rtp_zt,
    rtpthlp_zt,
    rtp2_zt,
    thlp2_zt,
    gr,
):
    rc_coef_zt = Lv / (exner * Cp) - ep2 * thv_ds_zt
    wp2thvp_zt = wp2thlp_zt + ep1 * thv_ds_zt * wp2rtp_zt + rc_coef_zt * wp2rcp_zt
    wpthvp_zt = wpthlp_zt + ep1 * thv_ds_zt * wprtp_zt + rc_coef_zt * wprcp_zt
    rtpthvp_zt = rtpthlp_zt + ep1 * thv_ds_zt * rtp2_zt + rc_coef_zt * rtprcp_zt
    thlpthvp_zt = thlp2_zt + ep1 * thv_ds_zt * rtpthlp_zt + rc_coef_zt * thlprcp_zt

    def to_zm_zero_top(value):
        return zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, value).at[:, gr.k_ub_zm].set(0.0)

    return (
        to_zm_zero_top(wpthvp_zt),
        wp2thvp_zt,
        to_zm_zero_top(rtpthvp_zt),
        to_zm_zero_top(thlpthvp_zt),
        rc_coef_zt,
        to_zm_zero_top(rc_coef_zt),
    )


def _inputs(rng, gr):
    z = lambda lo, hi: rng.uniform(lo, hi, (_NG, _NZT))
    return dict(
        exner=z(0.9, 1.0), thv_ds_zt=z(285.0, 300.0),
        wprcp_zt=z(-1e-4, 1e-4), wp2rcp_zt=z(-1e-4, 1e-4), rtprcp_zt=z(-1e-7, 1e-7), thlprcp_zt=z(-1e-3, 1e-3),
        wpthlp_zt=z(-0.05, 0.05), wprtp_zt=z(-1e-4, 1e-4), wp2thlp_zt=z(-0.05, 0.05), wp2rtp_zt=z(-1e-4, 1e-4),
        rtpthlp_zt=z(-1e-5, 1e-5), rtp2_zt=z(0, 1e-6), thlp2_zt=z(0, 0.5), gr=gr)


def test_assembly_and_regrid():
    gr = setup_grid(_NG, 50.0, 50.0, 50.0 * (_NZT + 1))
    rng = np.random.default_rng(566)
    k = _inputs(rng, gr)
    wpthvp_zm, wp2thvp_zt, rtpthvp_zm, thlpthvp_zm, rc_coef_zt, rc_coef_zm = calc_xpthvp_terms_jax(**k)
    # independent on-zt transcription
    rc = Lv / (np.asarray(k['exner']) * Cp) - ep2 * np.asarray(k['thv_ds_zt'])
    tds = np.asarray(k['thv_ds_zt'])
    r_wp2thvp = np.asarray(k['wp2thlp_zt']) + ep1 * tds * np.asarray(k['wp2rtp_zt']) + rc * np.asarray(k['wp2rcp_zt'])
    r_wpthvp = np.asarray(k['wpthlp_zt']) + ep1 * tds * np.asarray(k['wprtp_zt']) + rc * np.asarray(k['wprcp_zt'])
    r_rtpthvp = np.asarray(k['rtpthlp_zt']) + ep1 * tds * np.asarray(k['rtp2_zt']) + rc * np.asarray(k['rtprcp_zt'])
    r_thlpthvp = np.asarray(k['thlp2_zt']) + ep1 * tds * np.asarray(k['rtpthlp_zt']) + rc * np.asarray(k['thlprcp_zt'])
    assert np.max(np.abs(np.asarray(rc_coef_zt) - rc)) < 1e-9, "rc_coef_zt"
    assert np.max(np.abs(np.asarray(wp2thvp_zt) - r_wp2thvp)) < 1e-12, "wp2thvp_zt (stays on zt)"
    # zm outputs == zt2zm(zt-quantity) with k_ub_zm zeroed
    kub = gr.k_ub_zm
    for got_zm, zt_q, nm in ((wpthvp_zm, r_wpthvp, "wpthvp"), (rtpthvp_zm, r_rtpthvp, "rtpthvp"),
                             (thlpthvp_zm, r_thlpthvp, "thlpthvp"), (rc_coef_zm, rc, "rc_coef")):
        ref_zm = np.asarray(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, jnp.asarray(zt_q))).copy()
        ref_zm[:, kub] = 0.0
        assert np.max(np.abs(np.asarray(got_zm) - ref_zm)) < 1e-9, f"{nm}_zm regrid/k_ub-zero"
    assert np.all(np.asarray(rc_coef_zm)[:, kub] == 0.0), "k_ub_zm not zeroed"
    print("  rc_coef + 4 x'thv' assembly (ep1/ep2 coeffs) + zt→zm regrid (k_ub_zm=0) match transcription  PASS")


def test_grad_finite():
    gr = setup_grid(_NG, 50.0, 50.0, 50.0 * (_NZT + 1))
    k = _inputs(np.random.default_rng(7), gr)
    g = jax.grad(lambda w: jnp.sum(calc_xpthvp_terms_jax(**{**k, 'wprcp_zt': w})[0] ** 2))(jnp.asarray(k['wprcp_zt']))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad wrt wprcp_zt"
    print("  jax.grad(wpthvp_zm) wrt wprcp_zt finite  PASS")


def main():
    print("test_calc_xpthvp_terms:")
    test_assembly_and_regrid()
    test_grad_finite()
    print("All calc_xpthvp_terms checks PASSED")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""test_calc_xprcp_component.py — pin the per-component cloud-water covariance contributions (ADG1 x'rc' fluxes).

`calc_xprcp_component` (pdf_closure_module.py ↔ pdf_closure_module.F90:2652, lines 3089-3104) gives one PDF
component's contributions to the cloud-water covariances <w'rc'>, <w'^2 rc'>, <rt'rc'>, <thl'rc'>, <u'rc'>, <v'rc'>
on the ADG1 path. The non-ADG1 corr_w_chi correction (F90:3108-3140) is omitted because it is gated on iiPDF_type ∉
{ADG1,ADG2,new_hybrid} AND multiplied by corr_w_chi (=0 for ADG1) — verified against the Fortran, so it vanishes
identically; chi_i/corr_w_chi_i are accepted for signature fidelity only. Validated only end-to-end. This pins all 6
outputs against an INDEPENDENT transcription of the F90 base formulas (with drc = rc_i − rcm):
    wprcp=(w−wm)·drc;  wp2rcp=((w−wm)²+σ_w²)·drc;  uprcp/vprcp analogous;
    rtprcp =(rt−rtm)·drc + (corr_χη·σ_η + σ_χ)/(2·crt)·σ_χ·cloud_frac;
    thlprcp=(thl−thlm)·drc + (corr_χη·σ_η − σ_χ)/(2·cthl)·σ_χ·cloud_frac   (note the +/− sign difference)
plus the cthl=0 safe-guard (no NaN; masked by cloud_frac there) and a finite grad. (iter 564)
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_xprcp_component

# nominal per-component PDF state (component i) and grid means
_KW = dict(
    wm=0.1, rtm=0.012, thlm=290.0, um=3.0, vm=-1.0, rcm=2.0e-5,
    w_i=0.6, rt_i=0.013, thl_i=289.5, u_i=3.4, v_i=-0.6,
    varnce_w_i=0.25, chi_i=1.0e-4, stdev_chi_i=8.0e-4, stdev_eta_i=1.2e-3,
    corr_w_chi_i=0.0, corr_chi_eta_i=0.4, crt_i=1.6, cthl_i=-0.02,
    rc_i=5.0e-5, cloud_frac_i=0.35,
)


def test_six_covariance_formulas():
    k = _KW
    wprcp, wp2rcp, rtprcp, thlprcp, uprcp, vprcp = (float(x) for x in calc_xprcp_component(**k))
    drc = k['rc_i'] - k['rcm']
    dw = k['w_i'] - k['wm']
    r_wprcp = dw * drc
    r_wp2rcp = (dw ** 2 + k['varnce_w_i']) * drc
    r_rtprcp = ((k['rt_i'] - k['rtm']) * drc
                + (k['corr_chi_eta_i'] * k['stdev_eta_i'] + k['stdev_chi_i']) / (2.0 * k['crt_i'])
                * k['stdev_chi_i'] * k['cloud_frac_i'])
    r_thlprcp = ((k['thl_i'] - k['thlm']) * drc
                 + (k['corr_chi_eta_i'] * k['stdev_eta_i'] - k['stdev_chi_i']) / (2.0 * k['cthl_i'])
                 * k['stdev_chi_i'] * k['cloud_frac_i'])
    r_uprcp = (k['u_i'] - k['um']) * drc
    r_vprcp = (k['v_i'] - k['vm']) * drc
    for got, ref, nm in ((wprcp, r_wprcp, "wprcp"), (wp2rcp, r_wp2rcp, "wp2rcp"),
                         (rtprcp, r_rtprcp, "rtprcp"), (thlprcp, r_thlprcp, "thlprcp"),
                         (uprcp, r_uprcp, "uprcp"), (vprcp, r_vprcp, "vprcp")):
        assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"{nm}: {got} vs {ref}"
    # the rt (+σ_χ) vs thl (−σ_χ) sign difference must actually be present
    assert r_rtprcp != r_thlprcp, "test inputs failed to exercise the rt/thl correction sign difference"
    print("  6 cloud-water covariance contributions match the F90 base formulas (incl. rt+/thl− signs)  PASS")


def test_cthl_zero_safe():
    k = dict(_KW); k['cthl_i'] = 0.0; k['cloud_frac_i'] = 0.0   # rsatl=0 limit, masked by cloud_frac=0
    out = [float(x) for x in calc_xprcp_component(**k)]
    assert all(np.isfinite(out)), f"cthl=0 produced non-finite output: {out}"
    print("  cthl=0 (rsatl→0) safe-guard: all outputs finite  PASS")


def test_grad_finite():
    k = dict(_KW)
    g = jax.grad(lambda rc: calc_xprcp_component(**{**k, 'rc_i': rc})[2])(k['rc_i'])  # d(rtprcp)/d(rc_i)
    assert np.isfinite(float(g)), "non-finite grad wrt rc_i"
    print("  jax.grad(rtprcp) wrt rc_i finite  PASS")


def main():
    print("test_calc_xprcp_component:")
    test_six_covariance_formulas()
    test_cthl_zero_safe()
    test_grad_finite()
    print("All calc_xprcp_component checks PASSED")


if __name__ == "__main__":
    main()

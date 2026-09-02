#!/usr/bin/env python3
"""test_hydrometeor_mixed_moments.py — validate the hydrometeor_mixed_moments top driver.

The driver is pure orchestration over the already-validated integral functions (univar/bivar/covar). Its own
risk is wiring: which PDF param goes into which integral, the chi/eta->rt/thl correlation transforms, the
recomputed binormal means, and the triangular hmx/hmy loop. The oracle is therefore a LITERAL per-level,
per-hydrometeor Python transcription of the Fortran k/hm_idx/hmy_idx loops calling the same validated
integrals on scalars — the vectorized (over nzt) driver must reproduce it exactly. Plus a finite jax.grad.
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

from clubb_jax.src.CLUBB_core.pdf_utilities import compute_mean_binormal, calc_corr_rt_x, calc_corr_thl_x
from clubb_jax.src.Microphys.mixed_moment_PDF_integrals import (
    hydrometeor_mixed_moments, xphmp_integral_covar, xp_a_hmpb_integrals_all_MM, hmxphmyp_integral_covar)

NZT, HM_DIM = 12, 3


def _build_inputs(seed=7):
    rng = np.random.default_rng(seed)
    sc = lambda lo, hi: rng.uniform(lo, hi, NZT)
    col = lambda lo, hi: rng.uniform(lo, hi, (NZT, HM_DIM))
    p = dict(
        hydromet=col(1e-5, 3e-4),
        mu_w_1=sc(-0.5, 0.5), mu_w_2=sc(-0.5, 0.5),
        mu_rt_1=sc(1e-3, 1e-2), mu_rt_2=sc(1e-3, 1e-2),
        mu_thl_1=sc(290, 300), mu_thl_2=sc(290, 300),
        sigma_w_1=sc(0.2, 0.8), sigma_w_2=sc(0.2, 0.8),
        sigma_rt_1=sc(1e-4, 1e-3), sigma_rt_2=sc(1e-4, 1e-3),
        sigma_thl_1=sc(0.3, 1.0), sigma_thl_2=sc(0.3, 1.0),
        sigma_chi_1=sc(1e-4, 1e-3), sigma_chi_2=sc(1e-4, 1e-3),
        sigma_eta_1=sc(1e-4, 1e-3), sigma_eta_2=sc(1e-4, 1e-3),
        mixt_frac=sc(0.25, 0.75), precip_frac_1=sc(0.3, 0.9), precip_frac_2=sc(0.3, 0.9),
        crt_1=sc(0.5, 1.5), crt_2=sc(0.5, 1.5), cthl_1=sc(-0.02, -0.005), cthl_2=sc(-0.02, -0.005),
        mu_hm_1=col(1e-5, 3e-4), mu_hm_2=col(1e-5, 3e-4),
        sigma_hm_1=col(1e-5, 1e-4), sigma_hm_2=col(1e-5, 1e-4),
        mu_hm_1_n=col(-11, -8), mu_hm_2_n=col(-11, -8),
        sigma_hm_1_n=col(0.3, 0.8), sigma_hm_2_n=col(0.3, 0.8),
        corr_chi_hm_1=col(-0.7, 0.7), corr_chi_hm_2=col(-0.7, 0.7),
        corr_eta_hm_1=col(-0.7, 0.7), corr_eta_hm_2=col(-0.7, 0.7),
        corr_w_hm_1_n=col(-0.7, 0.7), corr_w_hm_2_n=col(-0.7, 0.7),
        corr_hmx_hmy_1=rng.uniform(-0.6, 0.6, (NZT, HM_DIM, HM_DIM)),
        corr_hmx_hmy_2=rng.uniform(-0.6, 0.6, (NZT, HM_DIM, HM_DIM)),
        hydromet_tol=np.array([1e-12, 1e-12, 1e-12]),
        rt_tol=1e-8, thl_tol=1e-2, w_tol=2e-2)
    return {k: (jnp.asarray(v) if isinstance(v, np.ndarray) else v) for k, v in p.items()}


def _ref(p):
    """Literal per-level/per-hydrometeor transcription of the Fortran loops."""
    g = lambda k: float(k)
    rt = np.zeros((NZT, HM_DIM)); th = np.zeros((NZT, HM_DIM))
    w2 = np.zeros((NZT, HM_DIM)); hh = np.zeros((NZT, HM_DIM, HM_DIM))
    P = {k: (np.asarray(v)) for k, v in p.items()}
    for k in range(NZT):
        mixt = P['mixt_frac'][k]
        rtm = float(compute_mean_binormal(P['mu_rt_1'][k], P['mu_rt_2'][k], mixt))
        thlm = float(compute_mean_binormal(P['mu_thl_1'][k], P['mu_thl_2'][k], mixt))
        wm = float(compute_mean_binormal(P['mu_w_1'][k], P['mu_w_2'][k], mixt))
        for hm in range(HM_DIM):
            crh1 = float(calc_corr_rt_x(P['crt_1'][k], P['sigma_rt_1'][k], P['sigma_chi_1'][k],
                                        P['sigma_eta_1'][k], P['corr_chi_hm_1'][k, hm], P['corr_eta_hm_1'][k, hm]))
            crh2 = float(calc_corr_rt_x(P['crt_2'][k], P['sigma_rt_2'][k], P['sigma_chi_2'][k],
                                        P['sigma_eta_2'][k], P['corr_chi_hm_2'][k, hm], P['corr_eta_hm_2'][k, hm]))
            cth1 = float(calc_corr_thl_x(P['cthl_1'][k], P['sigma_thl_1'][k], P['sigma_chi_1'][k],
                                         P['sigma_eta_1'][k], P['corr_chi_hm_1'][k, hm], P['corr_eta_hm_1'][k, hm]))
            cth2 = float(calc_corr_thl_x(P['cthl_2'][k], P['sigma_thl_2'][k], P['sigma_chi_2'][k],
                                         P['sigma_eta_2'][k], P['corr_chi_hm_2'][k, hm], P['corr_eta_hm_2'][k, hm]))
            ht = P['hydromet_tol'][hm]
            rt[k, hm] = float(xphmp_integral_covar(
                P['mu_rt_1'][k], P['mu_rt_2'][k], P['mu_hm_1'][k, hm], P['mu_hm_2'][k, hm],
                P['sigma_rt_1'][k], P['sigma_rt_2'][k], P['sigma_hm_1'][k, hm], P['sigma_hm_2'][k, hm],
                crh1, crh2, mixt, P['precip_frac_1'][k], P['precip_frac_2'][k], rtm, P['rt_tol'], ht))
            th[k, hm] = float(xphmp_integral_covar(
                P['mu_thl_1'][k], P['mu_thl_2'][k], P['mu_hm_1'][k, hm], P['mu_hm_2'][k, hm],
                P['sigma_thl_1'][k], P['sigma_thl_2'][k], P['sigma_hm_1'][k, hm], P['sigma_hm_2'][k, hm],
                cth1, cth2, mixt, P['precip_frac_1'][k], P['precip_frac_2'][k], thlm, P['thl_tol'], ht))
            w2[k, hm] = float(xp_a_hmpb_integrals_all_MM(
                P['mu_w_1'][k], P['mu_w_2'][k], P['mu_hm_1'][k, hm], P['mu_hm_2'][k, hm],
                P['mu_hm_1_n'][k, hm], P['mu_hm_2_n'][k, hm], P['sigma_w_1'][k], P['sigma_w_2'][k],
                P['sigma_hm_1'][k, hm], P['sigma_hm_2'][k, hm], P['sigma_hm_1_n'][k, hm], P['sigma_hm_2_n'][k, hm],
                P['corr_w_hm_1_n'][k, hm], P['corr_w_hm_2_n'][k, hm], mixt, P['precip_frac_1'][k],
                P['precip_frac_2'][k], wm, P['hydromet'][k, hm], P['w_tol'], ht, 2, 1))
            for hmy in range(hm + 1, HM_DIM):
                hh[k, hmy, hm] = float(hmxphmyp_integral_covar(
                    P['mu_hm_1'][k, hm], P['mu_hm_2'][k, hm], P['mu_hm_1'][k, hmy], P['mu_hm_2'][k, hmy],
                    P['sigma_hm_1'][k, hm], P['sigma_hm_2'][k, hm], P['sigma_hm_1'][k, hmy], P['sigma_hm_2'][k, hmy],
                    P['corr_hmx_hmy_1'][k, hm, hmy], P['corr_hmx_hmy_2'][k, hm, hmy], mixt,
                    P['precip_frac_1'][k], P['precip_frac_2'][k], P['hydromet'][k, hm], P['hydromet'][k, hmy],
                    ht, P['hydromet_tol'][hmy]))
    return rt, th, w2, hh


def test_driver_vs_literal_loop():
    p = _build_inputs()
    out = hydrometeor_mixed_moments(p)
    rt, th, w2, hh = _ref(p)
    for name, got, ref in (("rtphmp", out['rtphmp_zt'], rt), ("thlphmp", out['thlphmp_zt'], th),
                           ("wp2hmp", out['wp2hmp'], w2), ("hmxphmyp", out['hmxphmyp_zt'], hh)):
        got = np.asarray(got)
        rel = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-30))
        assert rel < 1e-12, f"{name} vs literal loop rel {rel:.2e}"
    print(f"  hydrometeor_mixed_moments (nzt={NZT}, hm_dim={HM_DIM}): all 4 outputs vs literal Fortran-loop "
          f"transcription rel <1e-12  PASS")


def test_differentiable():
    p = _build_inputs()
    def loss(sig_w_1):
        q = dict(p); q['sigma_w_1'] = sig_w_1
        out = hydrometeor_mixed_moments(q)
        return jnp.sum(out['wp2hmp'] ** 2) + jnp.sum(out['rtphmp_zt'] ** 2)
    g = jax.grad(loss)(p['sigma_w_1'])
    assert np.isfinite(np.asarray(g)).all(), "non-finite grad through hydrometeor_mixed_moments"
    print(f"  jax.grad(hydrometeor_mixed_moments) wrt sigma_w_1: finite (||g||={float(jnp.linalg.norm(g)):.3e})  PASS")


def test_compute_mean_binormal_f2py():
    """The literal-loop oracle above CALLS compute_mean_binormal (pdf_utilities.F90) for rtm/thlm/wm, so a bug in
    it would cancel between driver and oracle. Validate it independently against the f2py Fortran oracle to break
    that circularity. SKIPs if clubb_f2py is unbuilt. (iter 411)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py compute_mean_binormal oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(200):
        mu1, mu2, mf = float(rng.uniform(-10, 10)), float(rng.uniform(-10, 10)), float(rng.uniform(0, 1))
        j = float(compute_mean_binormal(mu1, mu2, mf))
        f = float(clubb_f2py.f2py_compute_mean_binormal(mu1, mu2, mf))
        worst = max(worst, abs(j - f))
    assert worst < 1e-13, f"compute_mean_binormal f2py mismatch {worst:.2e}"
    print(f"  compute_mean_binormal vs f2py oracle (200 cases): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_hydrometeor_mixed_moments:")
    for t in (test_driver_vs_literal_loop, test_differentiable, test_compute_mean_binormal_f2py):
        t()
    print("All hydrometeor_mixed_moments checks PASSED")


if __name__ == "__main__":
    main()

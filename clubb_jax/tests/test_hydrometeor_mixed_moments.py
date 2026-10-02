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

from types import SimpleNamespace
import pytest
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


def _run_interface(p):
    """Pack independent scalar-loop fixtures into the source interface's types."""
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
    gr = setup_grid(1, 100., 100., 100. * (NZT + 1))
    metadata = SimpleNamespace(iiPDF_w=0, iiPDF_chi=1, iiPDF_eta=2,
        iirr=0, iiNr=1, iiri=2, iiPDF_rr=3, iiPDF_Nr=4, iiPDF_ri=5,
        hydromet_tol=p['hydromet_tol'], hydromet_list=('rrm','Nrm','rim'))
    pdf = SimpleNamespace(mixt_frac=p['mixt_frac'][None,:])
    hm_pdf = SimpleNamespace()
    mu=[];sigma=[];corr=[]
    for i in (1,2):
        for name in ('rt','thl'):
            setattr(pdf,f'{name}_{i}',p[f'mu_{name}_{i}'][None,:])
            setattr(pdf,f'varnce_{name}_{i}',p[f'sigma_{name}_{i}'][None,:]**2)
        for name in ('crt','cthl'):
            setattr(pdf,f'{name}_{i}',p[f'{name}_{i}'][None,:])
        for name in ('mu_hm','sigma_hm','corr_chi_hm','corr_eta_hm','corr_hmx_hmy'):
            setattr(hm_pdf,f'{name}_{i}',p[f'{name}_{i}'][None,...])
        mu.append(jnp.concatenate((p[f'mu_w_{i}'][:,None],jnp.zeros((NZT,2)),p[f'mu_hm_{i}_n']),axis=-1)[None,...])
        sigma.append(jnp.concatenate((jnp.stack([p[f'sigma_{name}_{i}'] for name in ('w','chi','eta')],axis=-1),p[f'sigma_hm_{i}_n']),axis=-1)[None,...])
        corr.append(jnp.zeros((1,NZT,6,6)).at[0,:,3:,0].set(p[f'corr_w_hm_{i}_n']))
    frac=SimpleNamespace(**{f'precip_frac_{i}':p[f'precip_frac_{i}'][None,:] for i in (1,2)})
    stats=JaxStats.empty(l_sample=True,names=('rrpNrp','rrprip','Nrprip'),
        grids=('zm',)*3,ncol=1,max_nlev=gr.nzm)
    rt,th,w2,stats=hydrometeor_mixed_moments(gr,1,NZT,6,HM_DIM,
        p['hydromet'][None,...],metadata,mu[0],mu[1],sigma[0],sigma[1],
        corr[0],corr[1],pdf,hm_pdf,frac,stats)
    return dict(rtphmp_zt=rt[0],thlphmp_zt=th[0],wp2hmp=w2[0],stats=stats)


def test_driver_vs_literal_loop():
    p = _build_inputs()
    out = _run_interface(p)
    rt, th, w2, hh = _ref(p)
    for name, got, ref in (("rtphmp", out['rtphmp_zt'], rt), ("thlphmp", out['thlphmp_zt'], th),
                           ("wp2hmp", out['wp2hmp'], w2)):
        got = np.asarray(got)
        rel = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-30))
        assert rel < 1e-12, f"{name} vs literal loop rel {rel:.2e}"
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid, zt2zm
    gr=setup_grid(1,100.,100.,100.*(NZT+1))
    for slot,(i,j) in enumerate(((0,1),(0,2),(1,2))):
        expected=zt2zm(gr.nzm,gr.nzt,1,gr,jnp.asarray(hh[:,j,i])[None,:])
        np.testing.assert_allclose(out['stats'].buffers[1][slot],expected,rtol=1.e-12,atol=1.e-25)
    print(f"  hydrometeor_mixed_moments (nzt={NZT}, hm_dim={HM_DIM}): all 4 outputs vs literal Fortran-loop "
          f"transcription rel <1e-12  PASS")


def test_differentiable():
    p = _build_inputs()
    def loss(sig_w_1):
        q = dict(p); q['sigma_w_1'] = sig_w_1
        out = _run_interface(q)
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
        pytest.skip(f"f2py compute_mean_binormal oracle unavailable: {e}")
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

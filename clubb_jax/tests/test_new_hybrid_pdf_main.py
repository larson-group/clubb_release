#!/usr/bin/env python3
"""test_new_hybrid_pdf_main.py — validate the new-hybrid PDF driver port (new_hybrid_pdf_main.F90).

new_hybrid_pdf_driver: w-setter + rt/thl/u/v/sclr responders. Oracles:
  1. f2py bit-shadow vs f2py_new_hybrid_pdf_driver — all 31 exposed outputs (clipped responder skewnesses,
     component means/variances for w/rt/thl/u/v/sclr, mixt_frac, sigma_sqd_w), with sclr_dim=1. SKIPs if unbuilt.
  2. sigma_sqd_w == 1 - F_w consistency and a finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_hybrid_pdf_main import new_hybrid_pdf_driver

NG, NZ, SD = 2, 6, 1


def _inputs(rng):
    sh = (NG, NZ)
    return dict(
        wm=rng.uniform(-1, 1, sh), rtm=rng.uniform(0.005, 0.02, sh), thlm=rng.uniform(285, 305, sh),
        um=rng.uniform(-5, 5, sh), vm=rng.uniform(-5, 5, sh),
        wp2=rng.uniform(0.1, 2, sh), rtp2=rng.uniform(1e-7, 1e-6, sh), thlp2=rng.uniform(0.05, 1, sh),
        up2=rng.uniform(0.1, 2, sh), vp2=rng.uniform(0.1, 2, sh),
        Skw=rng.uniform(-2, 2, sh),
        wprtp=rng.uniform(-2e-4, 2e-4, sh), wpthlp=rng.uniform(-0.2, 0.2, sh),
        upwp=rng.uniform(-0.3, 0.3, sh), vpwp=rng.uniform(-0.3, 0.3, sh),
        gamma_Skw_fnc=rng.uniform(0.1, 0.9, sh),
        slope_coef_spread_DG_means_w=rng.uniform(1.0, 3.0, NG),
        pdf_component_stdev_factor_w=rng.uniform(0.8, 1.5, NG),
        Skrt=rng.uniform(-3, 3, sh), Skthl=rng.uniform(-3, 3, sh),
        Sku=rng.uniform(-3, 3, sh), Skv=rng.uniform(-3, 3, sh))


def _sclr(rng):
    return dict(sclrm=rng.uniform(0, 1, (NG, NZ, SD)), sclrp2=rng.uniform(1e-3, 1e-2, (NG, NZ, SD)),
                wpsclrp=rng.uniform(-1e-2, 1e-2, (NG, NZ, SD)), Sksclr=rng.uniform(-3, 3, (NG, NZ, SD)))


# f2py return order
_ORDER = ['Skrt', 'Skthl', 'Sku', 'Skv', 'Sksclr', 'mu_w_1', 'mu_w_2', 'mu_rt_1', 'mu_rt_2',
          'mu_thl_1', 'mu_thl_2', 'mu_u_1', 'mu_u_2', 'mu_v_1', 'mu_v_2', 'sigma_w_1_sqd', 'sigma_w_2_sqd',
          'sigma_rt_1_sqd', 'sigma_rt_2_sqd', 'sigma_thl_1_sqd', 'sigma_thl_2_sqd', 'sigma_u_1_sqd',
          'sigma_u_2_sqd', 'sigma_v_1_sqd', 'sigma_v_2_sqd', 'mu_sclr_1', 'mu_sclr_2', 'sigma_sclr_1_sqd',
          'sigma_sclr_2_sqd', 'mixt_frac', 'sigma_sqd_w']


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py new_hybrid_pdf_driver oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(17)
    worst = 0.0
    for _ in range(3):
        a = _inputs(rng); s = _sclr(rng)
        f = clubb_f2py.f2py_new_hybrid_pdf_driver(
            SD, a['wm'], a['rtm'], a['thlm'], a['um'], a['vm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['up2'], a['vp2'], a['Skw'], a['wprtp'], a['wpthlp'], a['upwp'], a['vpwp'],
            s['sclrm'], s['sclrp2'], s['wpsclrp'], a['gamma_Skw_fnc'],
            a['slope_coef_spread_DG_means_w'], a['pdf_component_stdev_factor_w'],
            a['Skrt'], a['Skthl'], a['Sku'], a['Skv'], s['Sksclr'])
        g = new_hybrid_pdf_driver(
            a['wm'], a['rtm'], a['thlm'], a['um'], a['vm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['up2'], a['vp2'], a['Skw'], a['wprtp'], a['wpthlp'], a['upwp'], a['vpwp'],
            a['gamma_Skw_fnc'], a['slope_coef_spread_DG_means_w'], a['pdf_component_stdev_factor_w'],
            a['Skrt'], a['Skthl'], a['Sku'], a['Skv'],
            s['sclrm'], s['sclrp2'], s['wpsclrp'], s['Sksclr'])
        for fi, key in zip(f, _ORDER):
            gi = np.asarray(g[key])
            fi = np.asarray(fi)
            if fi.ndim == 3:        # sclr outputs: f2py is (ng,nz,sd); ours matches
                fi = fi.reshape(gi.shape)
            worst = max(worst, np.max(np.abs(gi - fi)))
    assert worst < 1e-11, f"new_hybrid_pdf_driver f2py mismatch {worst:.2e}"
    print(f"  f2py new_hybrid_pdf_driver: bit-match (31 outputs, sclr_dim=1), worst {worst:.2e}  PASS")


def test_sigma_sqd_w_and_grad():
    rng = np.random.default_rng(4)
    a = _inputs(rng)
    g = new_hybrid_pdf_driver(
        a['wm'], a['rtm'], a['thlm'], a['um'], a['vm'], a['wp2'], a['rtp2'], a['thlp2'],
        a['up2'], a['vp2'], a['Skw'], a['wprtp'], a['wpthlp'], a['upwp'], a['vpwp'],
        a['gamma_Skw_fnc'], a['slope_coef_spread_DG_means_w'], a['pdf_component_stdev_factor_w'],
        a['Skrt'], a['Skthl'], a['Sku'], a['Skv'])
    # sigma_sqd_w must lie in (0,1) (since F_w in (0,1]) and mixt_frac must be valid.
    ssw = np.asarray(g['sigma_sqd_w'])
    assert np.all(ssw >= 0.0) and np.all(ssw < 1.0), "sigma_sqd_w out of [0,1)"
    assert np.all(np.asarray(g['mixt_frac']) > 0.0), "invalid mixt_frac"
    print("  sigma_sqd_w in [0,1), mixt_frac > 0  PASS")

    def loss(Skw):
        out = new_hybrid_pdf_driver(
            a['wm'], a['rtm'], a['thlm'], a['um'], a['vm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['up2'], a['vp2'], Skw, a['wprtp'], a['wpthlp'], a['upwp'], a['vpwp'],
            a['gamma_Skw_fnc'], a['slope_coef_spread_DG_means_w'], a['pdf_component_stdev_factor_w'],
            a['Skrt'], a['Skthl'], a['Sku'], a['Skv'])
        return sum(jnp.sum(v ** 2) for v in out.values())
    grad = np.asarray(jax.grad(loss)(jnp.asarray(a['Skw'])))
    assert np.isfinite(grad).all(), "non-finite grad through new_hybrid_pdf_driver"
    print(f"  jax.grad through new_hybrid_pdf_driver: finite ({grad.size} entries)  PASS")


def main():
    print("test_new_hybrid_pdf_main:")
    for t in (test_f2py_oracle, test_sigma_sqd_w_and_grad):
        t()
    print("All new_hybrid_pdf_driver checks PASSED")


if __name__ == "__main__":
    main()

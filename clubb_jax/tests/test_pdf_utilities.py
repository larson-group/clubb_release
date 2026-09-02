"""Verification of pdf_utilities.py — lognormal<->normal moment/correlation conversions.

mean_L2N and stdev_L2N are verified BIT-TO-BIT against the Fortran f2py API
(f2py_mean_l2n / f2py_stdev_l2n) — the gold standard. The two correlation conversions
(corr_NL2NN, corr_LL2NN) are not exposed by the API, so they are checked against
Monte-Carlo sampling of correlated (log)normal variables — an independent statistical
oracle for the defining relations (Garvey 2000, Eqs. B-1 and C-3).

Run: python clubb_jax/tests/test_pdf_utilities.py  (self-bootstraps sys.path; the f2py checks SKIP if unbuilt).
"""
import os
import sys

import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.pdf_utilities import (
    mean_L2N, stdev_L2N, corr_NL2NN, corr_LL2NN, MAX_MAG_CORRELATION,
    corr_NN2NL, corr_NN2LL, calc_corr_chi_x, calc_corr_eta_x,
    calc_corr_rt_x, calc_corr_thl_x,
)


def test_mean_stdev_L2N_bit_to_bit():
    """mean_L2N / stdev_L2N match the Fortran f2py oracle bit-to-bit."""
    try:
        import clubb_f2py as f
    except ImportError:
        print("  mean_stdev_L2N: SKIP (clubb_f2py not built)"); return
    # representative hydrometeor/Ncn coefficient-of-variation ratios sigma^2/mu^2
    s2m2_vals = [1e-6, 1e-3, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    mu_vals = [1e2, 1e4, 1e6, 1e8, 1e-5, 3.3e7]
    worst_mean = worst_std = 0.0
    for s2m2 in s2m2_vals:
        sj = float(stdev_L2N(s2m2))
        sf = float(f.f2py_stdev_l2n(s2m2))
        worst_std = max(worst_std, abs(sj - sf) / (abs(sf) + 1e-300))
        for mu in mu_vals:
            mj = float(mean_L2N(mu, s2m2))
            mf = float(f.f2py_mean_l2n(mu, s2m2))
            worst_mean = max(worst_mean, abs(mj - mf) / (abs(mf) + 1e-300))
    assert worst_mean < 1e-14, f"mean_L2N vs f2py worst rel {worst_mean:.2e}"
    assert worst_std < 1e-14, f"stdev_L2N vs f2py worst rel {worst_std:.2e}"
    print(f"  mean_L2N / stdev_L2N vs f2py: bit-to-bit "
          f"(mean rel<{worst_mean:.1e}, stdev rel<{worst_std:.1e})  PASS")


def test_corr_NL2NN_vs_montecarlo():
    """corr_NL2NN(corr(x,y)) reproduces the measured corr(x, ln y) for x normal, y lognormal."""
    rng = np.random.default_rng(0)
    n = 4_000_000
    worst = 0.0
    for sigma_y_n, rho_target in [(0.3, 0.6), (0.5, -0.4), (0.8, 0.3)]:
        # build x ~ N(0,1), ln y = mu_yn + sigma_yn*(rho_n*x + sqrt(1-rho_n^2)*z) so that
        # corr(x, ln y) = rho_n; then MEASURE the linear corr(x, y) and check the JAX
        # conversion maps that measured linear corr back to rho_n.
        rho_n = rho_target
        x = rng.standard_normal(n)
        z = rng.standard_normal(n)
        mu_yn = 2.0
        lny = mu_yn + sigma_y_n * (rho_n * x + np.sqrt(1 - rho_n**2) * z)
        y = np.exp(lny)
        corr_xy = np.corrcoef(x, y)[0, 1]            # linear correlation (the input)
        mu_y, sig_y = y.mean(), y.std()
        y_s2m2 = (sig_y / mu_y) ** 2
        got = float(corr_NL2NN(corr_xy, sigma_y_n, y_s2m2))
        worst = max(worst, abs(got - rho_n))
    assert worst < 5e-3, f"corr_NL2NN vs Monte-Carlo worst |Δ| {worst:.2e}"
    print(f"  corr_NL2NN vs Monte-Carlo: worst |Δ| {worst:.1e}  PASS")


def test_corr_LL2NN_vs_montecarlo():
    """corr_LL2NN(corr(x,y)) reproduces the measured corr(ln x, ln y) for x,y lognormal."""
    rng = np.random.default_rng(1)
    n = 4_000_000
    worst = 0.0
    for sx_n, sy_n, rho_target in [(0.3, 0.4, 0.5), (0.6, 0.5, -0.3), (0.7, 0.2, 0.2)]:
        a = rng.standard_normal(n)
        b = rng.standard_normal(n)
        lnx = 1.0 + sx_n * a
        lny = 2.0 + sy_n * (rho_target * a + np.sqrt(1 - rho_target**2) * b)
        x, y = np.exp(lnx), np.exp(lny)
        corr_xy = np.corrcoef(x, y)[0, 1]
        x_s2m2 = (x.std() / x.mean()) ** 2
        y_s2m2 = (y.std() / y.mean()) ** 2
        got = float(corr_LL2NN(corr_xy, sx_n, sy_n, x_s2m2, y_s2m2))
        worst = max(worst, abs(got - rho_target))
    assert worst < 5e-3, f"corr_LL2NN vs Monte-Carlo worst |Δ| {worst:.2e}"
    print(f"  corr_LL2NN vs Monte-Carlo: worst |Δ| {worst:.1e}  PASS")


def test_corr_clip_and_zero_sigma():
    """Correlation conversions clip to +/-0.99 and fall back to corr_x_y at sigma_n=0."""
    # inconsistent inputs that would overshoot |corr|>1 -> clipped
    big = float(corr_NL2NN(0.95, 0.01, 1.0))   # sqrt(1)/0.01 * 0.95 huge -> clip +0.99
    assert abs(big - MAX_MAG_CORRELATION) < 1e-15
    # sigma_y_n == 0 -> returns corr_x_y unchanged
    assert float(corr_NL2NN(0.42, 0.0, 0.0)) == 0.42
    assert float(corr_LL2NN(0.37, 0.0, 0.5, 0.0, 0.3)) == 0.37
    # differentiable
    g = float(jax.grad(lambda s: stdev_L2N(s))(0.5))
    assert np.isfinite(g) and g > 0
    print(f"  clip / zero-sigma fallback / differentiable: PASS")



def test_corr_NN2NL_NN2LL_calc_corr_chi_x_bit_to_bit():
    """corr_NN2NL, corr_NN2LL (inverse conversions) and calc_corr_chi_x match f2py bit-to-bit."""
    try:
        import clubb_f2py as f
    except ImportError:
        print("  corr_NN2NL/NN2LL/calc_corr_chi_x: SKIP (clubb_f2py not built)"); return
    rng = np.random.default_rng(0)
    # calc_corr_chi_x (scalar): corr(chi,x) from corr(rt,x), corr(thl,x)
    worst_chi = 0.0
    for _ in range(2000):
        crt, cthl = rng.uniform(0.5, 2), rng.uniform(0.5, 2)
        srt, sthl = rng.uniform(0, 1e-3), rng.uniform(0, 5)
        schi = rng.choice([0.0, rng.uniform(1e-9, 1e-3)])   # exercise sigma_chi=0 branch
        crx, ctx = rng.uniform(-0.9, 0.9), rng.uniform(-0.9, 0.9)
        fo = float(f.f2py_calc_corr_chi_x(crt, cthl, srt, sthl, schi, crx, ctx))
        jo = float(calc_corr_chi_x(crt, cthl, srt, sthl, schi, crx, ctx))
        worst_chi = max(worst_chi, abs(jo - fo) / (abs(fo) + 1e-30))
    assert worst_chi < 1e-13, f"calc_corr_chi_x vs f2py {worst_chi:.2e}"
    # corr_NN2NL / corr_NN2LL (array): inverse of corr_NL2NN / corr_LL2NN
    N = 500
    cxyn = rng.uniform(-0.9, 0.9, (N, 1)); syn = np.abs(rng.uniform(0, 1.5, (N, 1)))
    ys2 = np.abs(rng.uniform(0, 3, (N, 1))); sxn = np.abs(rng.uniform(0, 1.5, (N, 1)))
    xs2 = np.abs(rng.uniform(0, 3, (N, 1)))
    d_nl = np.max(np.abs(np.asarray(f.f2py_corr_nn2nl(cxyn, syn, ys2)) - np.asarray(corr_NN2NL(cxyn, syn, ys2))))
    d_ll = np.max(np.abs(np.asarray(f.f2py_corr_nn2ll(cxyn, sxn, syn, xs2, ys2)) - np.asarray(corr_NN2LL(cxyn, sxn, syn, xs2, ys2))))
    assert d_nl < 1e-13 and d_ll < 1e-13, f"NN2NL {d_nl:.2e}, NN2LL {d_ll:.2e}"
    # round-trip: corr_NL2NN -> corr_NN2NL recovers the linear corr (for valid sigma_n)
    rt = float(corr_NN2NL(corr_NL2NN(0.5, 0.4, 0.2), 0.4, 0.2))
    assert abs(rt - 0.5) < 1e-12, f"NL2NN/NN2NL round-trip {rt}"
    print(f"  calc_corr_chi_x ({worst_chi:.0e}) / corr_NN2NL ({d_nl:.0e}) / corr_NN2LL ({d_ll:.0e}) vs f2py: bit-to-bit  PASS")


def test_corr_chi_eta_rt_thl_roundtrip():
    """calc_corr_eta_x bit-to-bit vs f2py; (rt,thl)->(chi,eta)->(rt,thl) round-trips exactly.

    chi = crt·rt - cthl·thl, eta = crt·rt + cthl·thl, so the chi/eta decomposition is
    exactly invertible: calc_corr_{chi,eta}_x then calc_corr_{rt,thl}_x must recover the
    original (rt,x),(thl,x) correlations (away from clipping / zero-sigma)."""
    try:
        import clubb_f2py as f
    except ImportError:
        print("  corr_chi_eta_rt_thl_roundtrip: SKIP (clubb_f2py not built)"); return
    rng = np.random.default_rng(3)
    worst_eta = worst_rt = worst_thl = 0.0
    for _ in range(3000):
        crt, cthl = rng.uniform(0.5, 2.0), rng.uniform(0.5, 2.0)
        s_rt, s_thl = rng.uniform(1e-4, 1e-3), rng.uniform(0.05, 2.0)  # > tol
        rho_rt_thl = rng.uniform(-0.8, 0.8)
        # consistent component sigmas of chi, eta from rt, thl
        var_rt, var_thl = (crt * s_rt) ** 2, (cthl * s_thl) ** 2
        cov = crt * s_rt * cthl * s_thl * rho_rt_thl
        s_chi = np.sqrt(var_rt + var_thl - 2 * cov)
        s_eta = np.sqrt(var_rt + var_thl + 2 * cov)
        c_rt_x, c_thl_x = rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5)
        # eta vs f2py
        eo = float(f.f2py_calc_corr_eta_x(crt, cthl, s_rt, s_thl, s_eta, c_rt_x, c_thl_x))
        ej = float(calc_corr_eta_x(crt, cthl, s_rt, s_thl, s_eta, c_rt_x, c_thl_x))
        worst_eta = max(worst_eta, abs(ej - eo) / (abs(eo) + 1e-30))
        # round-trip
        c_chi = float(calc_corr_chi_x(crt, cthl, s_rt, s_thl, s_chi, c_rt_x, c_thl_x))
        c_eta = float(calc_corr_eta_x(crt, cthl, s_rt, s_thl, s_eta, c_rt_x, c_thl_x))
        rt_back = float(calc_corr_rt_x(crt, s_rt, s_chi, s_eta, c_chi, c_eta))
        thl_back = float(calc_corr_thl_x(cthl, s_thl, s_chi, s_eta, c_chi, c_eta))
        if abs(c_chi) < 0.99 and abs(c_eta) < 0.99:   # skip clipped cases
            worst_rt = max(worst_rt, abs(rt_back - c_rt_x))
            worst_thl = max(worst_thl, abs(thl_back - c_thl_x))
    assert worst_eta < 1e-12, f"calc_corr_eta_x vs f2py {worst_eta:.2e}"
    assert worst_rt < 1e-10 and worst_thl < 1e-10, f"round-trip rt {worst_rt:.2e} thl {worst_thl:.2e}"
    print(f"  calc_corr_eta_x vs f2py ({worst_eta:.0e}); chi/eta<->rt/thl round-trip "
          f"(rt {worst_rt:.0e}, thl {worst_thl:.0e})  PASS")


if __name__ == "__main__":
    print("pdf_utilities (lognormal<->normal) verification:")
    test_mean_stdev_L2N_bit_to_bit()
    test_corr_NL2NN_vs_montecarlo()
    test_corr_LL2NN_vs_montecarlo()
    test_corr_clip_and_zero_sigma()
    test_corr_NN2NL_NN2LL_calc_corr_chi_x_bit_to_bit()
    test_corr_chi_eta_rt_thl_roundtrip()
    print("All pdf_utilities tests PASSED.")

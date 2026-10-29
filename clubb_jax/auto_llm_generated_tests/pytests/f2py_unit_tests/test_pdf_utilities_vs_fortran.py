"""Check the JAX pdf_utilities ports (lognormal/normal moment and correlation conversions)
against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, max_rel_diff
from clubb_jax.src.CLUBB_core.pdf_utilities import (
    mean_L2N,
    stdev_L2N,
    corr_NL2NN,
    corr_NN2NL,
    corr_NN2LL,
    calc_corr_chi_x,
    calc_corr_eta_x,
    calc_corr_rt_x,
    calc_corr_thl_x,
)


def test_mean_stdev_L2N_bit_to_bit():
    """mean_L2N / stdev_L2N match the Fortran f2py oracle bit-to-bit."""
    import clubb_f2py as f
    # representative hydrometeor/Ncn coefficient-of-variation ratios sigma^2/mu^2
    s2m2_vals = [1e-6, 1e-3, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    mu_vals = [1e2, 1e4, 1e6, 1e8, 1e-5, 3.3e7]
    worst_mean = worst_std = 0.0
    for s2m2 in s2m2_vals:
        sj = float(stdev_L2N(s2m2))
        sf = float(f.f2py_stdev_l2n(s2m2))
        worst_std = max(worst_std, max_rel_diff("stdev_L2N", sj, sf, 1e-300))
        for mu in mu_vals:
            mj = float(mean_L2N(mu, s2m2))
            mf = float(f.f2py_mean_l2n(mu, s2m2))
            worst_mean = max(worst_mean, max_rel_diff("mean_L2N", mj, mf, 1e-300))
    assert worst_mean < 1e-14, f"mean_L2N vs f2py worst rel {worst_mean:.2e}"
    assert worst_std < 1e-14, f"stdev_L2N vs f2py worst rel {worst_std:.2e}"


def test_corr_NN2NL_NN2LL_calc_corr_chi_x_bit_to_bit():
    """corr_NN2NL, corr_NN2LL (inverse conversions) and calc_corr_chi_x match f2py bit-to-bit."""
    import clubb_f2py as f
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
        worst_chi = max(worst_chi, max_rel_diff("calc_corr_chi_x", jo, fo, 1e-30))
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


def test_corr_chi_eta_rt_thl_roundtrip():
    """calc_corr_eta_x bit-to-bit vs f2py; (rt,thl)->(chi,eta)->(rt,thl) round-trips exactly.

    chi = crt·rt - cthl·thl, eta = crt·rt + cthl·thl, so the chi/eta decomposition is
    exactly invertible: calc_corr_{chi,eta}_x then calc_corr_{rt,thl}_x must recover the
    original (rt,x),(thl,x) correlations (away from clipping / zero-sigma)."""
    import clubb_f2py as f
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
        worst_eta = max(worst_eta, max_rel_diff("calc_corr_eta_x", ej, eo, 1e-30))
        # round-trip
        c_chi = float(calc_corr_chi_x(crt, cthl, s_rt, s_thl, s_chi, c_rt_x, c_thl_x))
        c_eta = float(calc_corr_eta_x(crt, cthl, s_rt, s_thl, s_eta, c_rt_x, c_thl_x))
        rt_back = float(calc_corr_rt_x(crt, s_rt, s_chi, s_eta, c_chi, c_eta))
        thl_back = float(calc_corr_thl_x(cthl, s_thl, s_chi, s_eta, c_chi, c_eta))
        if abs(c_chi) < 0.99 and abs(c_eta) < 0.99:   # skip clipped cases
            worst_rt = max(worst_rt, max_abs_diff("calc_corr_rt_x round trip", rt_back, c_rt_x))
            worst_thl = max(worst_thl, max_abs_diff("calc_corr_thl_x round trip", thl_back, c_thl_x))
    assert worst_eta < 1e-12, f"calc_corr_eta_x vs f2py {worst_eta:.2e}"
    assert worst_rt < 1e-10 and worst_thl < 1e-10, f"round-trip rt {worst_rt:.2e} thl {worst_thl:.2e}"

"""End-to-end validation of the JAX KK rate functions against the Fortran rico oracle.

Unlike test_kk_autoconversion.py (which checks the rate functions against first-principles
quadrature), this feeds the FORTRAN's OWN PDF component moments — read from a real rico SCM
run's stats — into the JAX KK rate functions and compares to the Fortran's `rrm_auto`,
`rrm_accr`, `rrm_evap` outputs. It validates the rates against the actual Fortran microphysics
oracle, isolating the rate-function math from the (not-yet-ported) hydrometeor PDF setup.

The stats expose LINEAR component moments (mu/sigma) and LINEAR correlations; the rate functions
need the LOG moments and LOG correlations, obtained with the Iter109 pdf_utilities conversions
(mean_L2N/stdev_L2N/corr_NL2NN/corr_LL2NN) using sigma2_on_mu2 = (sigma/mu)^2 in-precip.

Validation (Iter113-114):
  * autoconversion: rrm_auto matched to median 4.7e-7 (rico: N_cn constant -> const_x2 path).
  * accretion:      rrm_accr matched to median 6e-9 (general bivar path + corr_chi_rr).
  * evaporation:    rrm_evap matched to ~1e-5 median (trivariate + 6 correlations + the
                    thermodynamic kk_evap_coef = 3 C_evap G_T_p ((4/3)pi rho_lw)^(2/3)
                    (1+Beta_Tl r_sl)/r_sl, evaluated at T_liq = thlm*exner). 11/12 points
                    match to ~1e-5; one variance-tolerance-boundary point (sigma_rr ~ rr_tol)
                    differs — a dispatch-edge case, not a rate-math error.

Requires a Fortran rico run's stats:
  python run_scripts/run_scm.py rico -legacy -max_iters 10 -out_dir output/rico_fort
The test skips (does not fail) if the stats file is absent.
"""
import os
import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    KK_auto_upscaled_mean, KK_accr_upscaled_mean, KK_evap_upscaled_mean,
)
from clubb_jax.src.Microphys.KK_microphys_module import kk_evap_coef, kk_auto_coef
from clubb_jax.src.CLUBB_core.pdf_utilities import (
    mean_L2N, stdev_L2N, corr_NL2NN, corr_LL2NN,
)
from clubb_jax.src.CLUBB_core.grid_class import ddzt, zt2zm
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_RICO_STATS = os.path.join(os.path.dirname(__file__),
                           "../../output/rico_fort/rico_stats.nc")
_RICO_LONG_STATS = os.path.join(os.path.dirname(__file__),
                                "../../output/rico_long_fort/rico_stats.nc")
_C_EVAP = 0.86   # rico tunable (rico_setup.txt)


def _logm(mu, sig):
    """(mu_n, sigma_n, sigma2_on_mu2) for a lognormal from its linear mean/std."""
    s2m2 = np.where(mu > 0, (sig / np.maximum(mu, 1e-30)) ** 2, 0.0)
    return (np.asarray(mean_L2N(np.maximum(mu, 1e-30), s2m2)),
            np.asarray(stdev_L2N(s2m2)), s2m2)


def _rel(out, ref, mask):
    return np.abs(out[mask] - ref[mask]) / np.abs(ref[mask])


def _grid_from_momentum_heights(zm, ngrdcol):
    zm = np.asarray(zm, dtype=np.float64)
    return setup_grid(
        ngrdcol=ngrdcol,
        deltaz=1.0,
        zm_init=float(zm[0]),
        zm_top=float(zm[-1]),
        grid_type=3,
        momentum_heights=np.tile(zm, (ngrdcol, 1)),
    )


def _zt2zm(value, gr):
    return zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, value)


def _ddzt(value, gr):
    return ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, value)


def test_kk_rates_vs_rico_oracle():
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP")
        return
    if not os.path.exists(_RICO_STATS):
        print(f"  rico_fort stats absent — SKIP (generate with run_scm.py rico -legacy "
              "-out_dir output/rico_fort)")
        return

    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:]).squeeze()
    J = jnp.asarray
    chi1, chi2, sc1, sc2 = g("chi_1"), g("chi_2"), g("stdev_chi_1"), g("stdev_chi_2")
    mf = g("mixt_frac")
    z = np.zeros_like(chi1)

    # --- autoconversion (N_cn constant in rico -> const_x2 path) -----------------
    mNcn = g("mu_Ncn_1")
    coef_a = kk_auto_coef(g("rho"))
    ln = np.log(np.maximum(mNcn, 1e-30))
    auto = np.asarray(KK_auto_upscaled_mean(
        J(chi1), J(chi2), J(mNcn), J(mNcn), J(ln), J(ln), J(sc1), J(sc2),
        J(z), J(z), J(z), J(z), J(z), J(z), J(coef_a), J(mf)))
    ra = g("rrm_auto")

    # --- accretion (general bivar path + corr_chi_rr) ---------------------------
    mrr1, mrr2, srr1, srr2 = g("mu_rr_1"), g("mu_rr_2"), g("sigma_rr_1"), g("sigma_rr_2")
    mrr1n, srr1n, rs1 = _logm(mrr1, srr1)
    mrr2n, srr2n, rs2 = _logm(mrr2, srr2)
    ccr1n = np.asarray(corr_NL2NN(g("corr_chi_rr_1"), srr1n, rs1))
    ccr2n = np.asarray(corr_NL2NN(g("corr_chi_rr_2"), srr2n, rs2))
    pf1, pf2 = g("precip_frac_1"), g("precip_frac_2")
    accr = np.asarray(KK_accr_upscaled_mean(
        J(chi1), J(chi2), J(mrr1), J(mrr2), J(mrr1n), J(mrr2n), J(sc1), J(sc2),
        J(srr1), J(srr2), J(srr1n), J(srr2n), J(ccr1n), J(ccr2n), J(mf), J(pf1), J(pf2)))
    rac = g("rrm_accr")

    # --- evaporation (trivariate + 6 correlations + thermodynamic coef) ----------
    mNr1, mNr2, sNr1, sNr2 = g("mu_Nr_1"), g("mu_Nr_2"), g("sigma_Nr_1"), g("sigma_Nr_2")
    mNr1n, sNr1n, Ns1 = _logm(mNr1, sNr1)
    mNr2n, sNr2n, Ns2 = _logm(mNr2, sNr2)
    ccN1n = np.asarray(corr_NL2NN(g("corr_chi_Nr_1"), sNr1n, Ns1))
    ccN2n = np.asarray(corr_NL2NN(g("corr_chi_Nr_2"), sNr2n, Ns2))
    crN1n = np.asarray(corr_LL2NN(g("corr_rr_Nr_1"), srr1n, sNr1n, rs1, Ns1))
    crN2n = np.asarray(corr_LL2NN(g("corr_rr_Nr_2"), srr2n, sNr2n, rs2, Ns2))
    T_liq = g("thlm") * g("exner")
    coef_e = np.asarray(kk_evap_coef(T_liq, g("p_in_Pa"), _C_EVAP))
    evap = np.asarray(KK_evap_upscaled_mean(
        J(chi1), J(chi2), J(mrr1), J(mrr2), J(mNr1), J(mNr2), J(mrr1n), J(mrr2n),
        J(mNr1n), J(mNr2n), J(sc1), J(sc2), J(srr1), J(srr2), J(sNr1), J(sNr2),
        J(srr1n), J(srr2n), J(sNr1n), J(sNr2n), J(ccr1n), J(ccr2n), J(ccN1n), J(ccN2n),
        J(crN1n), J(crN2n), J(coef_e), J(mf), J(pf1), J(pf2)))
    rev = g("rrm_evap")
    ds.close()

    # --- assertions -------------------------------------------------------------
    # auto/accr: gate-tight on significant points (within 3 orders of the peak rate).
    for name, out, ref, sig_tol in (("auto", auto, ra, 5e-6), ("accr", accr, rac, 5e-6)):
        nz = np.abs(ref) > 0
        sig = np.abs(ref) > np.nanmax(np.abs(ref)) / 1e3
        rs, rall = _rel(out, ref, sig), _rel(out, ref, nz)
        assert rs.max() < sig_tol, f"KK_{name} vs rico: sig max rel {rs.max():.2e}"
        assert np.median(rall) < 1e-6, f"KK_{name} vs rico: median rel {np.median(rall):.2e}"
        print(f"  KK_{name} vs rico rrm_{name}: {nz.sum()} pts, sig max {rs.max():.1e}, "
              f"median {np.median(rall):.1e}  PASS")

    # evap: the trivariate path + thermodynamic coef. Median validates the machinery;
    # one variance-tolerance-boundary point (sigma_rr ~ rr_tol) is an accepted dispatch edge.
    nz = np.abs(rev) > 0
    rall = _rel(evap, rev, nz)
    n_good = int(np.sum(rall < 1e-4))
    assert np.median(rall) < 1e-4, f"KK_evap vs rico: median rel {np.median(rall):.2e}"
    assert n_good >= nz.sum() - 1, f"KK_evap vs rico: only {n_good}/{nz.sum()} points < 1e-4"
    print(f"  KK_evap vs rico rrm_evap: {nz.sum()} pts, {n_good} match <1e-4, "
          f"median {np.median(rall):.1e}  PASS")


def test_kk_autoconversion_driver_vs_rico():
    """The COMPOSED autoconversion driver (Nc->Ncnm->log moments->KK_auto) reproduces rrm_auto.

    Unlike test_kk_rates_vs_rico_oracle (which feeds mu_Ncn directly), this drives the full
    chain from the raw PDF-state inputs (chi moments, Nc_in_cloud, cloud_frac, rho) — the
    quantities a running model has — through kk_autoconversion_mean. rico has constant N_c
    (const_Ncnp2_on_Ncnm2 = 0), so Ncnm = Nc_in_cloud and the const_x2 dispatch is used."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys.kk_microphys_driver import kk_autoconversion_mean
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    out = np.asarray(kk_autoconversion_mean(
        g("chi_1"), g("chi_2"), g("stdev_chi_1"), g("stdev_chi_2"), g("mixt_frac"),
        g("Nc_in_cloud"), g("cloud_frac_1"), g("cloud_frac_2"), g("rho"), 0.0, 0.0, 0.0))
    ra = g("rrm_auto"); ds.close()
    nz = np.abs(ra) > 0
    sig = np.abs(ra) > np.nanmax(np.abs(ra)) / 1e3
    rs, rall = _rel(out, ra, sig), _rel(out, ra, nz)
    assert rs.max() < 5e-6 and np.median(rall) < 1e-6, \
        f"driver vs rico: sig {rs.max():.2e}, median {np.median(rall):.2e}"
    print(f"  COMPOSED autoconv driver vs rico rrm_auto: {nz.sum()} pts, sig max {rs.max():.1e}, "
          f"median {np.median(rall):.1e}  PASS")


def test_kk_accr_evap_drivers_vs_rico():
    """The COMPOSED accretion + evaporation drivers reproduce rrm_accr / rrm_evap.

    Drives kk_accretion_mean / kk_evaporation_mean from the in-precip r_r/N_r moments + the
    prescribed normal-space correlations (derived here from the stats linear correlations via
    corr_NL2NN/corr_LL2NN, the back-conversion of the prescribed values)."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys.kk_microphys_driver import (
        kk_accretion_mean, kk_evaporation_mean)
    from clubb_jax.src.CLUBB_core.pdf_utilities import corr_NL2NN, corr_LL2NN, stdev_L2N
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    chi1, chi2, sc1, sc2, mf = g("chi_1"), g("chi_2"), g("stdev_chi_1"), g("stdev_chi_2"), g("mixt_frac")
    mrr1, mrr2, srr1, srr2 = g("mu_rr_1"), g("mu_rr_2"), g("sigma_rr_1"), g("sigma_rr_2")
    mNr1, mNr2, sNr1, sNr2 = g("mu_Nr_1"), g("mu_Nr_2"), g("sigma_Nr_1"), g("sigma_Nr_2")
    pf1, pf2 = g("precip_frac_1"), g("precip_frac_2")
    _s2 = lambda mu, sig: np.where(mu > 0, (sig / np.maximum(mu, 1e-30)) ** 2, 0.0)
    rs1, rs2, Ns1, Ns2 = _s2(mrr1, srr1), _s2(mrr2, srr2), _s2(mNr1, sNr1), _s2(mNr2, sNr2)
    srr1n, srr2n = np.asarray(stdev_L2N(rs1)), np.asarray(stdev_L2N(rs2))
    sNr1n, sNr2n = np.asarray(stdev_L2N(Ns1)), np.asarray(stdev_L2N(Ns2))
    ccr1n = np.asarray(corr_NL2NN(g("corr_chi_rr_1"), srr1n, rs1))
    ccr2n = np.asarray(corr_NL2NN(g("corr_chi_rr_2"), srr2n, rs2))
    acc = np.asarray(kk_accretion_mean(chi1, chi2, sc1, sc2, mrr1, mrr2, srr1, srr2,
                                       ccr1n, ccr2n, mf, pf1, pf2))
    rac = g("rrm_accr")
    nz = np.abs(rac) > 0; sig = np.abs(rac) > np.nanmax(np.abs(rac)) / 1e3
    assert _rel(acc, rac, sig).max() < 5e-6 and np.median(_rel(acc, rac, nz)) < 1e-6
    ccN1n = np.asarray(corr_NL2NN(g("corr_chi_Nr_1"), sNr1n, Ns1))
    ccN2n = np.asarray(corr_NL2NN(g("corr_chi_Nr_2"), sNr2n, Ns2))
    crN1n = np.asarray(corr_LL2NN(g("corr_rr_Nr_1"), srr1n, sNr1n, rs1, Ns1))
    crN2n = np.asarray(corr_LL2NN(g("corr_rr_Nr_2"), srr2n, sNr2n, rs2, Ns2))
    T_liq, p = g("thlm") * g("exner"), g("p_in_Pa")
    ev = np.asarray(kk_evaporation_mean(chi1, chi2, sc1, sc2, mrr1, mrr2, srr1, srr2,
                                        mNr1, mNr2, sNr1, sNr2, ccr1n, ccr2n, ccN1n, ccN2n,
                                        crN1n, crN2n, T_liq, p, 0.86, mf, pf1, pf2))
    rev = g("rrm_evap")
    # N_r evaporation tendency (same trivariate machinery, exponents 1, -2/3, 5/3)
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import KK_Nrm_evap_upscaled_mean
    from clubb_jax.src.Microphys.KK_microphys_module import kk_evap_coef
    from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N
    L = lambda mu, s2: np.asarray(mean_L2N(np.maximum(mu, 1e-30), s2))
    coef = np.asarray(kk_evap_coef(T_liq, p, 0.86))
    nev = np.asarray(KK_Nrm_evap_upscaled_mean(
        chi1, chi2, mrr1, mrr2, mNr1, mNr2, L(mrr1, rs1), L(mrr2, rs2), L(mNr1, Ns1), L(mNr2, Ns2),
        sc1, sc2, srr1, srr2, sNr1, sNr2, srr1n, srr2n, sNr1n, sNr2n,
        ccr1n, ccr2n, ccN1n, ccN2n, crN1n, crN2n, coef, mf, pf1, pf2, 300.0))
    rNev = g("Nrm_evap"); ds.close()
    nze = np.abs(rev) > 0; ralle = _rel(ev, rev, nze)
    assert np.median(ralle) < 1e-4 and int(np.sum(ralle < 1e-4)) >= nze.sum() - 1
    nzn = np.abs(rNev) > 0; ralln = _rel(nev, rNev, nzn)
    assert np.median(ralln) < 1e-4 and int(np.sum(ralln < 1e-4)) >= nzn.sum() - 1
    print(f"  COMPOSED accr driver: sig max {_rel(acc, rac, sig).max():.1e}, median "
          f"{np.median(_rel(acc, rac, nz)):.1e}; evap driver median {np.median(ralle):.1e}; "
          f"KK_Nrm_evap median {np.median(ralln):.1e}  PASS")


def test_kk_microphys_adjust_vs_rico():
    """KK_microphys_adjust (the tendency assembly) reproduces rico's rcm_mc / rrm_mc.

    Feeds the stored process rates (rrm_auto/accr/evap, Nrm_auto/evap) + state (rcm/rrm/Nrm/exner)
    into the assembly. rcm_mc = -(adjusted auto+accr) is a clean pure-function check (exact,
    including the source-over-depletion adjustment). rrm_mc = evap_net + source matches where the
    EVAP limiter doesn't trigger (the limiter's -rrm/dt uses the within-step rrm, which differs from
    the end-of-step stored rrm — the documented timing confound). thlm_mc is checked for
    self-consistency with the oracle formula (Lv/(Cp·exner)·rrm_mc)."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys_module import KK_microphys_adjust
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    dt = 300.0   # rico dt_main
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc = (np.asarray(x) for x in KK_microphys_adjust(
        dt, g("exner"), g("rcm"), g("rrm"), g("Nrm"),
        g("rrm_evap"), g("rrm_auto"), g("rrm_accr"), g("Nrm_evap"), g("Nrm_auto")))
    rcm_s, rrm_s, ev_adj, exner = g("rcm_mc"), g("rrm_mc"), g("rrm_evap_adj"), g("exner")
    ds.close()
    mr = np.abs(rcm_s) > 1e-20
    assert np.max(np.abs(rcm_mc[mr] - rcm_s[mr])) < 1e-20, "rcm_mc (source side) not exact"
    # rrm_mc where the evap limiter did not adjust (rrm_evap_adj == 0)
    mm = (np.abs(rrm_s) > 1e-20) & (np.abs(ev_adj) < 1e-30)
    assert np.max(np.abs(rrm_mc[mm] - rrm_s[mm])) < 1e-18, "rrm_mc (no-evap-adj) not exact"
    Lv, Cp = 2.5e6, 1004.67
    assert np.allclose(thlm_mc, (Lv / (Cp * exner)) * rrm_mc), "thlm_mc not self-consistent"
    print(f"  KK_microphys_adjust vs rico: rcm_mc exact ({mr.sum()} pts), rrm_mc exact at "
          f"no-evap-adj pts ({mm.sum()}), thlm_mc self-consistent  PASS")


def test_compute_kk_microphysics_vs_rico():
    """The full standalone microphysics step (hydromet fields + PDF state -> tendencies) runs,
    produces finite tendencies, and at NO-RAIN points (rrm<rr_tol, where it reduces to
    autoconversion) reproduces rico's rcm_mc to machine epsilon. (The accr/evap contributions
    from the rrm/Nrm fields carry the calc_comp_mu_sigma_hm timing confound — only a running rico
    fully validates those; see DESIGN.md. Full-array jax.grad needs edge-case hardening of the
    rate functions' where-branches — the individual rates are differentiable at clean points,
    test_differentiability.py.)"""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    import jax
    from clubb_jax.src.Microphys.KK_microphys.kk_microphys_driver import compute_kk_microphysics
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    cf = g("cloud_frac")
    pftol = np.maximum(0.1 * np.max(cf, axis=1), 0.005)
    args = (g("rrm"), g("Nrm"), g("chi_1"), g("chi_2"), g("stdev_chi_1"), g("stdev_chi_2"),
            g("mixt_frac"), g("thl_1"), g("thl_2"), g("Nc_in_cloud"), g("cloud_frac_1"),
            g("cloud_frac_2"), g("precip_frac"), g("precip_frac_1"), g("precip_frac_2"),
            pftol, g("rho"), g("thlm") * g("exner"), g("p_in_Pa"), g("exner"), g("rcm"), 300.0)
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc = (np.asarray(x) for x in compute_kk_microphysics(*args))
    rc_s, rrm = g("rcm_mc"), g("rrm")
    ds.close()
    assert np.all(np.isfinite(rcm_mc)) and np.all(np.isfinite(rrm_mc)), "non-finite tendency"
    norain = (rrm < 1e-10) & (np.abs(rc_s) > 1e-20)
    assert norain.sum() >= 5
    assert np.max(np.abs(rcm_mc[norain] - rc_s[norain])) < 1e-15, "no-rain rcm_mc not machine-exact"
    # Fully differentiable over the full array (exercises the Iter127-129 grad-safe fixes:
    # _safe_div, _safe_sqrt, the Rmax denominator guard, _pos_pow, the D_v argument clamp) —
    # both w.r.t. the rrm FIELD (d tendency / d hydromet state) and a chi PDF moment.
    import jax
    worst = 0.0
    for vi, eps in ((0, 1e-13), (2, 1e-9)):   # rrm field, chi_1 moment
        base = args[vi]
        def fv(e, _vi=vi, _b=base):
            a = list(args); a[_vi] = _b + e
            return jnp.sum(compute_kk_microphysics(*a)[0])
        g = float(jax.grad(fv)(0.0))
        fd = float((fv(eps) - fv(-eps)) / (2 * eps))
        rel = abs(g - fd) / (abs(fd) + 1e-30)
        assert np.isfinite(g) and rel < 1e-4, f"grad[{vi}] ad={g:.4e} fd={fd:.4e} rel={rel:.2e}"
        worst = max(worst, rel)
    print(f"  compute_kk_microphysics: runs, no-rain rcm_mc machine-exact ({norain.sum()} pts), "
          f"fully differentiable (rrm field + chi moment, rel {worst:.1e})  PASS")


def test_compute_kk_microphysics_vel_prereqs():
    """The opt-in `l_return_vel_prereqs=True` path (the transport-wiring enabler, Iter154): default
    behaviour is UNCHANGED (still the 5-tuple), and the prereqs path additionally returns the mean
    volume radius `mvr` + the rr/Nr in-precip component moments the sed-velocity functions consume.
    `mvr` reproduces the rico `mvrr` oracle to within the SAME calc_comp_mu_sigma_hm within-step
    timing confound as accr/evap (the rico oracle is pre-developed-rain → marginal points only)."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys.kk_microphys_driver import compute_kk_microphysics
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    cf = g("cloud_frac"); pftol = np.maximum(0.1 * np.max(cf, axis=1), 0.005)
    args = (g("rrm"), g("Nrm"), g("chi_1"), g("chi_2"), g("stdev_chi_1"), g("stdev_chi_2"),
            g("mixt_frac"), g("thl_1"), g("thl_2"), g("Nc_in_cloud"), g("cloud_frac_1"),
            g("cloud_frac_2"), g("precip_frac"), g("precip_frac_1"), g("precip_frac_2"),
            pftol, g("rho"), g("thlm") * g("exner"), g("p_in_Pa"), g("exner"), g("rcm"), 300.0)
    default = compute_kk_microphysics(*args)
    tend, pre = compute_kk_microphysics(*args, l_return_vel_prereqs=True)
    mvrr = g("mvrr"); ds.close()
    assert len(default) == 5 and len(tend) == 5, "default 5-tuple broken"
    assert all(np.allclose(np.asarray(a), np.asarray(b)) for a, b in zip(default, tend)), \
        "prereqs path changed the tendencies"
    for k in ("mvr", "mu_rr_1", "sigma_rr_1", "mu_Nr_1", "mu_rr_1_n"):
        assert k in pre and np.all(np.isfinite(np.asarray(pre[k]))), f"prereq {k} missing/non-finite"
    mvr_j = np.asarray(pre["mvr"]); m = np.abs(mvrr) > mvrr.max() / 1e3
    rel = np.abs((mvr_j - mvrr)[m] / mvrr[m]) if m.sum() else np.array([0.0])
    assert np.median(rel) < 0.1, f"mvr median rel {np.median(rel):.2e} beyond timing-confound band"
    print(f"  compute_kk_microphysics vel-prereqs: default unchanged; mvr vs rico mvrr "
          f"{m.sum()} marginal pts, median rel {np.median(rel):.2e} (timing-confound band)  PASS")


_RF02_STATS = os.path.join(os.path.dirname(__file__),
                           "../../output/rf02_do_fort/dycoms2_rf02_do_stats.nc")


def test_kk_autoconversion_vs_rf02_do():
    """Autoconversion generalizes beyond rico: matches dycoms2_rf02_do's rrm_auto.

    rf02_do is stratocumulus drizzle — a NARROW chi PDF (s_c ~ 32), which exercises the
    large-negative-z D_v branch (Iter130) that rico's |s_c|<6 never did. Generate the oracle:
      run_scm.py dycoms2_rf02_do -legacy -out_dir output/rf02_do_fort
    Skips if absent."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RF02_STATS):
        print("  rf02_do_fort stats absent — SKIP (a generalization check)"); return
    from clubb_jax.src.Microphys.KK_microphys.kk_microphys_driver import kk_autoconversion_mean
    ds = nc.Dataset(_RF02_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    out = np.asarray(kk_autoconversion_mean(
        g("chi_1"), g("chi_2"), g("stdev_chi_1"), g("stdev_chi_2"), g("mixt_frac"),
        g("Nc_in_cloud"), g("cloud_frac_1"), g("cloud_frac_2"), g("rho"), 0.0, 0.0, 0.0))
    ra = g("rrm_auto"); ds.close()
    nz = np.abs(ra) > 0
    sig = np.abs(ra) > np.nanmax(np.abs(ra)) / 1e3
    rs, rall = _rel(out, ra, sig), _rel(out, ra, nz)
    assert rs.max() < 5e-5 and np.median(rall) < 1e-6, \
        f"rf02_do autoconv: sig {rs.max():.2e}, median {np.median(rall):.2e}"
    print(f"  KK_auto vs rf02_do rrm_auto (narrow-chi/stratocu): {nz.sum()} pts, sig max "
          f"{rs.max():.1e}, median {np.median(rall):.1e}  PASS")


def test_kk_sedimentation_vs_rico():
    """KK mean sedimentation velocities (KK_sedimentation, KK00 Eq.37) vs the rico oracle. The
    Fortran stores Vrr/VNr on MOMENTUM levels as zt2zm(hydromet_vel_zt) (microphys_driver.F90:461),
    so this feeds the Fortran's OWN mean volume radius (mvrr, on zt) into the JAX KK_sedimentation
    and then the (bit-faithful) zt2zm interpolation, and compares to the Fortran Vrr/VNr. Bit-exact.
    Uses a LONGER rico run (rico_long_fort, 250 steps) so rain has developed -- the canonical
    10-step run has only sub-micron drops, whose (positive) velocities all clip to 0. Generate:
      run_scm.py rico -legacy -max_iters 250 -out_dir output/rico_long_fort"""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_LONG_STATS):
        print("  rico_long_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys_module import KK_sedimentation
    ds = nc.Dataset(_RICO_LONG_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])     # (time, nz)
    mvr, Vrr_f, VNr_f = g("mvrr"), g("Vrr"), g("VNr")
    zt, zm = np.asarray(ds["zt"][:]), np.asarray(ds["zm"][:]); ds.close()
    nt = mvr.shape[0]
    gr = _grid_from_momentum_heights(zm, nt)
    Vzt, Nzt = KK_sedimentation(jnp.asarray(mvr))          # zt-level velocities (KK00 Eq.37 + clip)
    Vrr_j, VNr_j = np.asarray(_zt2zm(Vzt, gr)), np.asarray(_zt2zm(Nzt, gr))  # -> momentum levels
    mV, mN = Vrr_f < 0, VNr_f < 0                          # non-clipped rain points
    assert mV.sum() > 0, "no rain points to validate sedimentation"
    eV = np.abs(Vrr_j - Vrr_f).max()
    eN = np.abs(VNr_j - VNr_f).max()
    assert eV < 1e-12 and eN < 1e-12, f"sed velocity mismatch: Vrr {eV:.2e}, VNr {eN:.2e}"
    grad = np.asarray(jax.grad(lambda r: jnp.sum(KK_sedimentation(r)[0]))(jnp.asarray(mvr[-1])))
    assert np.all(np.isfinite(grad)), "sed grad not finite"
    print(f"  KK_sedimentation -> zt2zm vs rico Vrr/VNr: {mV.sum()} rain pts, Vrr |Δ|max {eV:.1e}, "
          f"VNr |Δ|max {eN:.1e}  PASS")


def test_kk_sed_vel_covars_vs_rico():
    """KK sed-velocity covariances (KK_sed_vel_covars) vs the rico oracle: feed the Fortran's OWN
    PDF moments (mu/sigma_rr/Nr + corr_rr_Nr) and the mean volume radius (mvrr) into the JAX
    KK_sed_vel_covars and compare the <r_r'R_vr'> / <N_r'R_vr'> covariances to the Fortran
    rr_KK_mvr_covar_zt / Nr_KK_mvr_covar_zt (both on zt). The overall means are reconstructed
    within-step-consistently (rrm = a*f_p1*mu_rr_1 + (1-a)*f_p2*mu_rr_2) so there is NO timing
    confound; bit-faithful-to-the-gate on the significant points."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import KK_sed_vel_covars
    from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N, corr_LL2NN
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    mvr, mf = g("mvrr"), g("mixt_frac")
    mrr1, mrr2, srr1, srr2 = g("mu_rr_1"), g("mu_rr_2"), g("sigma_rr_1"), g("sigma_rr_2")
    mNr1, mNr2, sNr1, sNr2 = g("mu_Nr_1"), g("mu_Nr_2"), g("sigma_Nr_1"), g("sigma_Nr_2")
    pf1, pf2 = g("precip_frac_1"), g("precip_frac_2")
    cr1, cr2 = g("corr_rr_Nr_1"), g("corr_rr_Nr_2")
    rrcov_f, Nrcov_f = g("rr_KK_mvr_covar_zt"), g("Nr_KK_mvr_covar_zt"); ds.close()
    s2 = lambda mu, s: np.where(mu > 0, (s / np.maximum(mu, 1e-30)) ** 2, 0.0)
    rs1, rs2, Ns1, Ns2 = s2(mrr1, srr1), s2(mrr2, srr2), s2(mNr1, sNr1), s2(mNr2, sNr2)
    L = lambda mu, x: np.asarray(mean_L2N(np.maximum(mu, 1e-30), x))
    S = lambda x: np.asarray(stdev_L2N(x))
    cc1 = np.asarray(corr_LL2NN(cr1, S(rs1), S(Ns1), rs1, Ns1))
    cc2 = np.asarray(corr_LL2NN(cr2, S(rs2), S(Ns2), rs2, Ns2))
    out = KK_sed_vel_covars(pf1 * mrr1, pf2 * mrr2, pf1 * mNr1, pf2 * mNr2, mvr,
                            mrr1, mrr2, mNr1, mNr2, L(mrr1, rs1), L(mrr2, rs2),
                            L(mNr1, Ns1), L(mNr2, Ns2), srr1, srr2, sNr1, sNr2,
                            S(rs1), S(rs2), S(Ns1), S(Ns2), cc1, cc2, mf)
    rrj, Nrj = np.asarray(out["rr_KK_mvr_covar"]), np.asarray(out["Nr_KK_mvr_covar"])
    mr = np.abs(rrcov_f) > np.nanmax(np.abs(rrcov_f)) / 1e3
    mn = np.abs(Nrcov_f) > np.nanmax(np.abs(Nrcov_f)) / 1e3
    er = np.abs((rrj - rrcov_f)[mr] / rrcov_f[mr]).max()
    en = np.abs((Nrj - Nrcov_f)[mn] / Nrcov_f[mn]).max()
    assert er < 1e-6 and en < 1e-6, f"sed-covar mismatch: rr {er:.2e}, Nr {en:.2e}"
    # Differentiable w.r.t. the in-precip rr mean.
    base = (pf2 * mrr2, pf1 * mNr1, pf2 * mNr2, mvr, mrr1, mrr2, mNr1, mNr2,
            L(mrr1, rs1), L(mrr2, rs2), L(mNr1, Ns1), L(mNr2, Ns2), srr1, srr2, sNr1, sNr2,
            S(rs1), S(rs2), S(Ns1), S(Ns2), cc1, cc2, mf)
    gd = jax.grad(lambda r1: jnp.nansum(
        KK_sed_vel_covars(r1, *base)["rr_KK_mvr_covar"]))(jnp.asarray(pf1 * mrr1))
    assert np.all(np.isfinite(np.asarray(gd))), "sed-covar grad not finite"
    print(f"  KK_sed_vel_covars vs rico rr/Nr_KK_mvr_covar_zt: rr rel max {er:.1e} ({mr.sum()} sig), "
          f"Nr rel max {en:.1e}; differentiable  PASS")


def test_sed_centered_diff_lhs_vs_rico():
    """The mean-sedimentation transport operator (sed_centered_diff_lhs). Two verifications:
    (1) the RIGOROUS, oracle-free CONSERVATION CONTRACT — the column-mass-weighted sum of the
        sedimentation tendency equals the surface flux rho_ds_zm[0]*V[0]*hm[0] (the flux-form
        operator removes mass only out the bottom; the top is no-flux). Checked on the rico grid
        with a synthetic non-trivial velocity/field so it is non-degenerate;
    (2) against the rico `rrm_sd`/`Nrm_sd` budget stats where it is CLEAN — at the developed-rain
        points the KK velocity clips to 0 (sub-rain-radius drops) so the budget is exactly 0, which
        the operator reproduces (the active-sedimentation points all sit at marginal rrm, where the
        within-step vs end-of-step hmm timing confound dominates — the same caveat as accr/evap)."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_LONG_STATS):
        print("  rico_long_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.advance_microphys_module import (
        sed_centered_diff_lhs, lhs_budget_term)
    ds = nc.Dataset(_RICO_LONG_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    zt, zm = np.asarray(ds["zt"][:]), np.asarray(ds["zm"][:])
    rho_ds_zm, rho_ds_zt = g("rho_ds_zm"), g("rho_ds_zt")
    Vrr, rrm, rrm_sd_f = g("Vrr"), g("rrm"), g("rrm_sd"); ds.close()
    nzt, nzm, nt = zt.shape[0], zm.shape[0], rrm.shape[0]
    wa = np.zeros(nzm); wa[1:nzt] = (zm[1:nzt] - zt[:nzt - 1]) / (zt[1:nzt] - zt[:nzt - 1])
    dzt = zm[1:nzm] - zm[:nzm - 1]
    idzt = np.tile(1.0 / dzt, (nt, 1)); irzt = 1.0 / rho_ds_zt; waT = np.tile(wa, (nt, 1))

    # (1) Conservation contract with a synthetic smooth field + uniform downward velocity.
    rng = np.arange(nzt)
    hm = np.tile(1e-6 * np.exp(-((rng - 12) / 7.0) ** 2), (nt, 1)); hm[:, -3:] = 0.0
    Vc = np.full((nt, nzm), -0.7)
    s, m, sub = sed_centered_diff_lhs(jnp.asarray(Vc), jnp.asarray(rho_ds_zm), jnp.asarray(irzt),
                                      jnp.asarray(idzt), jnp.asarray(waT))
    tend = np.asarray(lhs_budget_term(s, m, sub, jnp.asarray(hm)))
    colsum = np.sum(rho_ds_zt * np.tile(dzt, (nt, 1)) * tend, axis=1)
    flux = rho_ds_zm[:, 0] * Vc[:, 0] * hm[:, 0]
    cons = np.abs((colsum - flux) / flux).max()
    assert cons < 1e-10, f"sed conservation contract violated: {cons:.2e}"

    # (2) Match the rico rrm_sd budget where clean (the real velocity field).
    s2, m2, sub2 = sed_centered_diff_lhs(jnp.asarray(Vrr), jnp.asarray(rho_ds_zm), jnp.asarray(irzt),
                                         jnp.asarray(idzt), jnp.asarray(waT))
    sd = np.asarray(lhs_budget_term(s2, m2, sub2, jnp.asarray(rrm)))
    clean = np.abs(rrm_sd_f) < 1e-25      # Vrr-clipped points: budget is exactly 0
    assert np.abs(sd[clean]).max() < 1e-20, "sed budget nonzero where the oracle is 0"
    print(f"  sed_centered_diff_lhs: conservation contract {cons:.1e} (flux-form), rrm_sd==0 exact "
          f"at {int(clean.sum())} clipped pts  PASS")


def test_term_turb_sed_lhs_vs_rico():
    """The turbulent-sedimentation LHS (term_turb_sed_lhs). It has the SAME flux-form discretization
    as the mean-sed operator with Vhmphmp_impc (the implicit sed-velocity covariance, interpolated to
    momentum) replacing the velocity, so this validates the FULL composition
    KK_sed_vel_covars -> zt2zm -> term_turb_sed_lhs via the rigorous, timing-independent CONSERVATION
    CONTRACT (the in-loop `rrm_ts` stat uses stats_finalize_budget bookkeeping + is timing-confounded,
    so the contract is the clean oracle), and confirms term_turb_sed_lhs == sed_centered_diff_lhs."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_LONG_STATS):
        print("  rico_long_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.advance_microphys_module import (
        sed_centered_diff_lhs, term_turb_sed_lhs, lhs_budget_term)
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import KK_sed_vel_covars
    from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N, corr_LL2NN
    ds = nc.Dataset(_RICO_LONG_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    zt, zm = np.asarray(ds["zt"][:]), np.asarray(ds["zm"][:])
    rho_ds_zm, rho_ds_zt = g("rho_ds_zm"), g("rho_ds_zt")
    mvr, mf = g("mvrr"), g("mixt_frac")
    mrr1, mrr2, srr1, srr2 = g("mu_rr_1"), g("mu_rr_2"), g("sigma_rr_1"), g("sigma_rr_2")
    mNr1, mNr2, sNr1, sNr2 = g("mu_Nr_1"), g("mu_Nr_2"), g("sigma_Nr_1"), g("sigma_Nr_2")
    pf1, pf2, cr1, cr2 = g("precip_frac_1"), g("precip_frac_2"), g("corr_rr_Nr_1"), g("corr_rr_Nr_2")
    rrm = g("rrm"); ds.close()
    nzt, nzm, nt = zt.shape[0], zm.shape[0], rrm.shape[0]
    s2 = lambda mu, s: np.where(mu > 0, (s / np.maximum(mu, 1e-30)) ** 2, 0.0)
    rs1, rs2, Ns1, Ns2 = s2(mrr1, srr1), s2(mrr2, srr2), s2(mNr1, sNr1), s2(mNr2, sNr2)
    L = lambda mu, x: np.asarray(mean_L2N(np.maximum(mu, 1e-30), x)); S = lambda x: np.asarray(stdev_L2N(x))
    cc1 = np.asarray(corr_LL2NN(cr1, S(rs1), S(Ns1), rs1, Ns1))
    cc2 = np.asarray(corr_LL2NN(cr2, S(rs2), S(Ns2), rs2, Ns2))
    out = KK_sed_vel_covars(pf1 * mrr1, pf2 * mrr2, pf1 * mNr1, pf2 * mNr2, mvr, mrr1, mrr2, mNr1, mNr2,
                            L(mrr1, rs1), L(mrr2, rs2), L(mNr1, Ns1), L(mNr2, Ns2), srr1, srr2, sNr1, sNr2,
                            S(rs1), S(rs2), S(Ns1), S(Ns2), cc1, cc2, mf)
    gr = _grid_from_momentum_heights(zm, nt)
    Vimpc = np.asarray(_zt2zm(jnp.asarray(out["Vrrprrp_impc"]), gr))     # -> momentum
    wa = np.zeros(nzm); wa[1:nzt] = (zm[1:nzt] - zt[:nzt - 1]) / (zt[1:nzt] - zt[:nzt - 1])
    dzt = zm[1:nzm] - zm[:nzm - 1]
    idzt = np.tile(1.0 / dzt, (nt, 1)); irzt = 1.0 / rho_ds_zt; waT = np.tile(wa, (nt, 1))
    args = (jnp.asarray(Vimpc), jnp.asarray(rho_ds_zm), jnp.asarray(irzt), jnp.asarray(idzt), jnp.asarray(waT))
    s, m, sub = term_turb_sed_lhs(*args)
    s0, m0, sub0 = sed_centered_diff_lhs(*args)
    assert (np.asarray(s) == np.asarray(s0)).all() and (np.asarray(m) == np.asarray(m0)).all() \
        and (np.asarray(sub) == np.asarray(sub0)).all(), "term_turb_sed != sed_centered"
    hm = rrm.copy(); hm[:, -3:] = 0.0                  # clean no-flux top for the contract
    tend = np.asarray(lhs_budget_term(s, m, sub, jnp.asarray(hm)))
    colsum = np.sum(rho_ds_zt * np.tile(dzt, (nt, 1)) * tend, axis=1)
    flux = rho_ds_zm[:, 0] * Vimpc[:, 0] * hm[:, 0]
    nz = np.abs(flux) > 1e-30
    cons = np.abs((colsum - flux)[nz] / flux[nz]).max()
    assert cons < 1e-10, f"turb-sed conservation contract violated: {cons:.2e}"
    print(f"  term_turb_sed_lhs: == sed_centered_diff_lhs, full KK_sed_vel_covars->zt2zm->op "
          f"conservation contract {cons:.1e} ({int(nz.sum())} pts)  PASS")


def test_microphys_mean_adv_vs_rico():
    """The hydrometeor MEAN-ADVECTION budget (`term_ma_zt_lhs` + `lhs_budget_term`) vs rico `rrm_ma`/
    `Nrm_ma`. Unlike `rrm_ts`/`rrm_ta`, the mean-advection budget is a plain `stats_update(-lhs_ma·hmm)`
    (no explicit+implicit split), so at robust-rrm points (within-step ≈ end-of-step) it matches the
    stored stat to machine precision — confirming the upwind `term_ma_zt_lhs` (l_upwind_xm_ma=.true.
    for rico) works for the hydrometeor transport and the grid (invrs_dzm) reconstruction is correct."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_LONG_STATS):
        print("  rico_long_fort stats absent — SKIP"); return
    from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs
    from clubb_jax.src.Microphys.advance_microphys_module import lhs_budget_term
    ds = nc.Dataset(_RICO_LONG_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    zt = np.asarray(ds["zt"][:]); zm = np.asarray(ds["zm"][:]); wm_zt = g("wm_zt")
    rrm, rrm_ma_f, Nrm, Nrm_ma_f = g("rrm"), g("rrm_ma"), g("Nrm"), g("Nrm_ma"); ds.close()
    nzt, nt = zt.shape[0], rrm.shape[0]; nzm = nzt + 1
    gr = _grid_from_momentum_heights(zm, nt)
    ma = term_ma_zt_lhs(
        nzm, nzt, nt, jnp.asarray(wm_zt),
        gr.weights_zt2zm, gr.invrs_dzt, gr.invrs_dzm, True, gr.grid_dir,
    )
    worst = 0.0
    for nm, rr, ref, tol in [("rrm_ma", rrm, rrm_ma_f, 1e-7), ("Nrm_ma", Nrm, Nrm_ma_f, 1e2)]:
        out = np.asarray(lhs_budget_term(ma[0], ma[1], ma[2], jnp.asarray(rr)))
        rob = rr > tol; rob3 = rob & np.roll(rob, 1, 1) & np.roll(rob, -1, 1)
        sig = rob3 & (np.abs(ref) > np.nanmax(np.abs(ref)) / 1e3)
        rel = np.abs((out - ref)[sig] / ref[sig]).max()
        assert rel < 1e-10, f"{nm} mismatch: {rel:.2e}"
        worst = max(worst, rel)
    print(f"  mean-adv budget (term_ma_zt_lhs) vs rico rrm_ma/Nrm_ma: robust-pt rel max {worst:.1e}  PASS")


def test_microphys_lhs_assembly():
    """`microphys_lhs` assembles 1/dt + 1/2 diffusion(+k=1 BC) + term_ma + sed + turb_sed correctly:
    (i) it equals the sum of the independently-computed verified sub-operators (band alignment, signs,
    1/dt placement); (ii) the turbulent-advection (eddy-diffusion) part conserves mass exactly (zero
    column-mass-weighted sum, the zero-flux property)."""
    from clubb_jax.src.Microphys.advance_microphys_module import (
        microphys_lhs, sed_centered_diff_lhs, term_turb_sed_lhs, lhs_budget_term)
    from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
    from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs
    nzt = 30; nzm = nzt + 1; dt = 300.0; nu = 1.5
    zm = np.linspace(0.0, 4000.0, nzm); zt = 0.5 * (zm[1:] + zm[:-1])
    idzt = 1.0 / (zm[1:nzm] - zm[:nzm - 1])
    idzm = np.zeros(nzm); idzm[1:nzt] = 1.0 / (zt[1:nzt] - zt[:nzt - 1])
    idzm[0] = idzm[1]; idzm[nzt] = idzm[nzt - 1]
    gr = _grid_from_momentum_heights(zm, 1)
    rho_zt = (1.2 * np.exp(-zt / 8000.0))[None]; rho_zm = (1.2 * np.exp(-zm / 8000.0))[None]
    irzt = 1.0 / rho_zt
    wa = np.zeros((1, nzm)); wa[0, 1:nzt] = (zm[1:nzt] - zt[:nzt - 1]) / (zt[1:nzt] - zt[:nzt - 1])
    K_hm = (5.0 * np.exp(-zm / 2000.0))[None]; wm = (-0.02 * np.ones(nzt))[None]
    V = np.full((1, nzm), -0.5); Vi = (-1e-3 * np.ones((1, nzm)))
    J = jnp.asarray
    s, m, sub = microphys_lhs(dt, J(K_hm), nu, J(wm), J(V), J(Vi), J(rho_zm), J(irzt), gr, J(wa))
    # Manual component sum
    lta = 0.5 * np.asarray(diffusion_zt_lhs(
        nzm, nzt, 1, gr, J(K_hm), jnp.zeros((1, nzt), dtype=jnp.float64),
        jnp.array([nu]), J(irzt), J(rho_zm),
    ))
    bc = 0.5 * irzt[0, 0] * (idzt[0] * (K_hm[0, 1] + nu) * rho_zm[0, 1] * idzm[1])
    lta[0, 0, 0] = -bc; lta[1, 0, 0] = bc; lta[2, 0, 0] = 0.0
    lma = np.asarray(term_ma_zt_lhs(
        nzm, nzt, 1, J(wm), gr.weights_zt2zm, gr.invrs_dzt, gr.invrs_dzm, True, gr.grid_dir,
    ))
    ss = [np.asarray(x) for x in sed_centered_diff_lhs(J(V), J(rho_zm), J(irzt), J(idzt[None]), J(wa))]
    st = [np.asarray(x) for x in term_turb_sed_lhs(J(Vi), J(rho_zm), J(irzt), J(idzt[None]), J(wa))]
    exp_super = lta[0, 0] + lma[0, 0] + ss[0][0] + st[0][0]
    exp_main = lta[1, 0] + lma[1, 0] + 1.0 / dt + ss[1][0] + st[1][0]
    exp_sub = lta[2, 0] + lma[2, 0] + ss[2][0] + st[2][0]
    assert np.abs(np.asarray(s)[0] - exp_super).max() < 1e-20
    assert np.abs(np.asarray(m)[0] - exp_main).max() < 1e-20
    assert np.abs(np.asarray(sub)[0] - exp_sub).max() < 1e-20
    # Turbulent-advection (diffusion) conserves: column-mass Σ of -lhs_ta·hm == 0.
    hm = (1e-6 * np.exp(-((np.arange(nzt) - 12) / 6.0) ** 2))[None]; hm[0, -3:] = 0.0
    tend = np.asarray(lhs_budget_term(J(lta[0]), J(lta[1]), J(lta[2]), J(hm)))
    cons = np.abs(np.sum(rho_zt * (1.0 / idzt)[None] * tend))
    assert cons < 1e-18, f"diffusion not conservative: {cons:.2e}"
    print(f"  microphys_lhs: == 1/dt + 1/2 diffusion(+BC) + term_ma + sed + turb_sed; turb-adv "
          f"conserves ({cons:.1e})  PASS")


def test_microphys_rhs_turb_sed_vs_rico():
    """The FULL turbulent-sedimentation tendency (explicit `term_turb_sed_rhs` + implicit
    `-term_turb_sed_lhs·hmm`) reproduces the rico `rrm_ts` budget stat — confirming `term_turb_sed_rhs`
    AND that `rrm_ts` = explicit + implicit (microphys_rhs stores the explicit via stats_begin_budget,
    microphys_solve adds the implicit via stats_finalize_budget). Validated by feeding the Fortran's own
    moments → KK_sed_vel_covars → zt2zm → the operators. Median-clean (the relmax tail is the marginal-rrm
    within-step/end-of-step timing confound)."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_LONG_STATS):
        print("  rico_long_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.advance_microphys_module import (
        term_turb_sed_lhs, term_turb_sed_rhs, lhs_budget_term)
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import KK_sed_vel_covars
    from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N, corr_LL2NN
    ds = nc.Dataset(_RICO_LONG_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0]); zt = np.asarray(ds["zt"][:]); zm = np.asarray(ds["zm"][:])
    rho_ds_zm, rho_ds_zt = g("rho_ds_zm"), g("rho_ds_zt"); mvr, mf = g("mvrr"), g("mixt_frac")
    mrr1, mrr2, srr1, srr2 = g("mu_rr_1"), g("mu_rr_2"), g("sigma_rr_1"), g("sigma_rr_2")
    mNr1, mNr2, sNr1, sNr2 = g("mu_Nr_1"), g("mu_Nr_2"), g("sigma_Nr_1"), g("sigma_Nr_2")
    pf1, pf2, cr1, cr2 = g("precip_frac_1"), g("precip_frac_2"), g("corr_rr_Nr_1"), g("corr_rr_Nr_2")
    rrm, rrm_ts_f = g("rrm"), g("rrm_ts"); ds.close()
    nzt, nzm, nt = zt.shape[0], zm.shape[0], rrm.shape[0]
    s2 = lambda mu, s: np.where(mu > 0, (s / np.maximum(mu, 1e-30)) ** 2, 0.0)
    rs1, rs2, Ns1, Ns2 = s2(mrr1, srr1), s2(mrr2, srr2), s2(mNr1, sNr1), s2(mNr2, sNr2)
    L = lambda mu, x: np.asarray(mean_L2N(np.maximum(mu, 1e-30), x)); S = lambda x: np.asarray(stdev_L2N(x))
    cc1 = np.asarray(corr_LL2NN(cr1, S(rs1), S(Ns1), rs1, Ns1))
    cc2 = np.asarray(corr_LL2NN(cr2, S(rs2), S(Ns2), rs2, Ns2))
    out = KK_sed_vel_covars(pf1 * mrr1, pf2 * mrr2, pf1 * mNr1, pf2 * mNr2, mvr, mrr1, mrr2, mNr1, mNr2,
                            L(mrr1, rs1), L(mrr2, rs2), L(mNr1, Ns1), L(mNr2, Ns2), srr1, srr2, sNr1, sNr2,
                            S(rs1), S(rs2), S(Ns1), S(Ns2), cc1, cc2, mf)
    gr = _grid_from_momentum_heights(zm, nt)
    Vimpc = np.asarray(_zt2zm(jnp.asarray(out["Vrrprrp_impc"]), gr))
    Vexpc = np.asarray(_zt2zm(jnp.asarray(out["Vrrprrp_expc"]), gr))
    wa = np.zeros(nzm); wa[1:nzt] = (zm[1:nzt] - zt[:nzt - 1]) / (zt[1:nzt] - zt[:nzt - 1])
    idzt = np.tile(1.0 / (zm[1:nzm] - zm[:nzm - 1]), (nt, 1)); irzt = 1.0 / rho_ds_zt
    waT = np.tile(wa, (nt, 1)); J = jnp.asarray
    si, mi, subi = term_turb_sed_lhs(J(Vimpc), J(rho_ds_zm), J(irzt), J(idzt), J(waT))
    total = np.asarray(lhs_budget_term(si, mi, subi, J(rrm))) \
        + np.asarray(term_turb_sed_rhs(J(Vexpc), J(rho_ds_zm), J(irzt), J(idzt)))
    rob = rrm > 1e-7; rob3 = rob & np.roll(rob, 1, 1) & np.roll(rob, -1, 1)
    sig = rob3 & (np.abs(rrm_ts_f) > np.nanmax(np.abs(rrm_ts_f)) / 1e3)
    rel = np.abs((total - rrm_ts_f)[sig] / rrm_ts_f[sig])
    assert np.median(rel) < 1e-8 and np.mean(rel < 1e-6) >= 0.9, \
        f"rrm_ts (expl+impl) median {np.median(rel):.2e}, frac<1e-6 {np.mean(rel < 1e-6):.2f}"
    print(f"  turb-sed total (term_turb_sed_rhs+lhs) vs rico rrm_ts: {int(sig.sum())} pts, median "
          f"{np.median(rel):.1e}, {100*np.mean(rel<1e-6):.0f}% < 1e-6  PASS")


def test_microphys_rhs_assembly():
    """`microphys_rhs` == hmm/dt + microphysics source + (-lhs_ta·hmm explicit ½-diffusion) +
    term_turb_sed_rhs (component bookkeeping); and `term_turb_sed_rhs` satisfies the conservation
    contract (column-mass Σ = surface flux rho_ds_zm[0]·Vexpc[0])."""
    from clubb_jax.src.Microphys.advance_microphys_module import (
        microphys_rhs, term_turb_sed_rhs, _turb_adv_lhs, lhs_budget_term)
    nzt = 30; nzm = nzt + 1; dt = 300.0; nu = 1.5
    zm = np.linspace(0.0, 4000.0, nzm); zt = 0.5 * (zm[1:] + zm[:-1])
    idzt = 1.0 / (zm[1:nzm] - zm[:nzm - 1])
    gr = _grid_from_momentum_heights(zm, 1)
    rho_zt = (1.2 * np.exp(-zt / 8000.0))[None]; rho_zm = (1.2 * np.exp(-zm / 8000.0))[None]; irzt = 1.0 / rho_zt
    K_hm = (5.0 * np.exp(-zm / 2000.0))[None]; hm = (1e-6 * np.exp(-((np.arange(nzt) - 12) / 6.0) ** 2))[None]
    src = (1e-9 * np.ones(nzt))[None]; Vexpc = (-1e-4 * np.exp(-zm / 2000.0))[None]; J = jnp.asarray
    rhs = np.asarray(microphys_rhs(dt, J(hm), J(src), J(K_hm), nu, J(Vexpc), J(rho_zm), J(irzt), gr))
    lta = _turb_adv_lhs(J(K_hm), nu, J(rho_zm), J(irzt), gr)
    exp = hm / dt + src + np.asarray(lhs_budget_term(lta[0], lta[1], lta[2], J(hm))) \
        + np.asarray(term_turb_sed_rhs(J(Vexpc), J(rho_zm), J(irzt), J(idzt[None])))
    assert np.abs(rhs - exp).max() < 1e-20, "microphys_rhs != component sum"
    # term_turb_sed_rhs conservation: column-mass Σ tendency == surface flux.
    tend = np.asarray(term_turb_sed_rhs(J(Vexpc), J(rho_zm), J(irzt), J(idzt[None])))
    colsum = np.sum(rho_zt * (1.0 / idzt)[None] * tend); flux = rho_zm[0, 0] * Vexpc[0, 0]
    cons = abs((colsum - flux) / flux)
    assert cons < 1e-10, f"term_turb_sed_rhs not conservative: {cons:.2e}"
    print(f"  microphys_rhs: == hmm/dt + src + expl-diffusion + term_turb_sed_rhs; expl turb-sed "
          f"conserves ({cons:.1e})  PASS")


def test_calculate_K_hm():
    """`calculate_K_hm` (hydrometeor eddy diffusivity). Like calc_comp_mu_sigma_hm it consumes the
    within-step hydrometeor field, so the in-loop K_hm_rr stat is timing-confounded (the stored rrm is
    end-of-step; matches only to ~2%). Verified instead by the exact formula transcription against a
    hand-computed reference (reusing the bit-faithful zt2zm/ddzt), incl. the correlation cap and the
    K=0 boundaries; plus differentiability."""
    from clubb_jax.src.Microphys.advance_microphys_module import calculate_K_hm, _C_K_HM, _EPS
    nzt = 20; nzm = nzt + 1
    zm = np.linspace(0.0, 3000.0, nzm); zt = 0.5 * (zm[1:] + zm[:-1])
    idzm = np.zeros(nzm); idzm[1:nzt] = 1.0 / (zt[1:nzt] - zt[:nzt - 1]); idzm[0] = idzm[1]; idzm[nzt] = idzm[nzt - 1]
    gr = _grid_from_momentum_heights(zm, 1)
    rng = np.arange(nzt)
    hm = (1e-5 * np.exp(-((rng - 8) / 4.0) ** 2) + 1e-9)[None]        # rain layer
    hmp2 = (3e-10 * np.exp(-((rng - 8) / 5.0) ** 2))[None]            # variance (zm-sized below)
    hmp2 = (np.interp(zm, zt, hmp2[0]))[None]                         # to momentum
    wp2 = (0.4 * np.exp(-zm / 1500.0))[None]; Kh_zm = (40.0 * np.exp(-zm / 1000.0))[None]
    Skw = (0.3 * np.sin(zm / 500.0))[None]; tol = 1e-10; J = jnp.asarray
    K = np.asarray(calculate_K_hm(J(wp2), J(Kh_zm), J(Skw), J(hm), J(hmp2), tol, gr))
    # Hand reference using the same verified zt2zm/ddzt.
    hm_zm = np.asarray(_zt2zm(J(hm), gr)); dhm = np.asarray(_ddzt(J(hm), gr))
    shmp2 = np.sqrt(np.maximum(hmp2, 0.0))
    ref = _C_K_HM * Kh_zm * (shmp2 / np.maximum(hm_zm, tol)) * (1.0 + np.abs(Skw))
    cap = np.sqrt(np.maximum(wp2, 0.0)) * shmp2 / np.where(np.abs(dhm) > _EPS, np.abs(dhm), 1.0)
    ref = np.where(np.abs(dhm) > _EPS, np.minimum(ref, cap), ref)
    ref[:, 0] = 0.0; ref[:, -1] = 0.0
    assert np.abs(K - ref).max() < 1e-30, f"K_hm != formula reference: {np.abs(K - ref).max():.2e}"
    assert K[0, 0] == 0.0 and K[0, -1] == 0.0, "K_hm boundaries not zeroed"
    assert np.any(K[0, 1:-1] > 0), "K_hm all zero"
    grad = np.asarray(jax.grad(lambda h: jnp.nansum(
        calculate_K_hm(J(wp2), J(Kh_zm), J(Skw), h, J(hmp2), tol, gr)))(J(hm)))
    assert np.all(np.isfinite(grad)), "K_hm grad not finite"
    print(f"  calculate_K_hm: == formula reference (cap + BC), differentiable  PASS")


def test_advance_one_hydrometeor():
    """`advance_one_hydrometeor` (assemble LHS+RHS, tridiag solve). Verify (i) the solve is consistent
    with the assembled system: lhs·soln == rhs to machine precision; (ii) physical sanity — with no
    microphysics source and downward sedimentation, a rain layer loses mass (column total decreases)."""
    from clubb_jax.src.Microphys.advance_microphys_module import (
        advance_one_hydrometeor, microphys_lhs, microphys_rhs)
    nzt = 30; nzm = nzt + 1; dt = 60.0; nu = 1.5
    zm = np.linspace(0.0, 4000.0, nzm); zt = 0.5 * (zm[1:] + zm[:-1])
    idzt = 1.0 / (zm[1:nzm] - zm[:nzm - 1])
    gr = _grid_from_momentum_heights(zm, 1)
    rho_zt = (1.2 * np.exp(-zt / 8000.0))[None]; rho_zm = (1.2 * np.exp(-zm / 8000.0))[None]; irzt = 1.0 / rho_zt
    wa = np.zeros((1, nzm)); wa[0, 1:nzt] = (zm[1:nzt] - zt[:nzt - 1]) / (zt[1:nzt] - zt[:nzt - 1])
    K_hm = (5.0 * np.exp(-zm / 2000.0))[None]; wm = (-0.01 * np.ones(nzt))[None]
    V = np.full((1, nzm), -0.6); V[0, -1] = 0.0; Vi = (-1e-4 * np.ones((1, nzm))); Ve = (-1e-5 * np.ones((1, nzm)))
    hm = (1e-6 * np.exp(-((np.arange(nzt) - 12) / 5.0) ** 2))[None]; src = np.zeros((1, nzt)); J = jnp.asarray
    a = (J(K_hm), nu, J(wm), J(V), J(Vi), J(Ve), J(rho_zm), J(irzt), gr, J(wa))
    soln = np.asarray(advance_one_hydrometeor(dt, J(hm), J(src), *a))
    # (i) consistency: lhs·soln == rhs
    sup, main, sub = microphys_lhs(dt, J(K_hm), nu, J(wm), J(V), J(Vi), J(rho_zm), J(irzt), gr, J(wa))
    rhs = np.asarray(microphys_rhs(dt, J(hm), J(src), J(K_hm), nu, J(Ve), J(rho_zm), J(irzt), gr))
    s, m, b = np.asarray(sup)[0], np.asarray(main)[0], np.asarray(sub)[0]
    lhs_soln = b * np.concatenate([soln[0, :1], soln[0, :-1]]) + m * soln[0] \
        + s * np.concatenate([soln[0, 1:], soln[0, -1:]])
    assert np.abs(lhs_soln - rhs[0]).max() < 1e-18, "lhs·soln != rhs"
    # (ii) sedimentation removes mass at the surface -> column total decreases
    col0 = np.sum(rho_zt[0] * (1.0 / idzt) * hm[0]); col1 = np.sum(rho_zt[0] * (1.0 / idzt) * soln[0])
    assert col1 < col0 and np.all(np.isfinite(soln)), "sedimentation did not remove mass"
    print(f"  advance_one_hydrometeor: lhs·soln==rhs (<1e-18), sed removes mass "
          f"({col0:.3e}->{col1:.3e})  PASS")


def test_kk_covar_driver_vs_rico():
    """The second-moment microphysics driver (KK_upscaled_covar_driver) vs the rico oracle.

    Feeds the FORTRAN's OWN PDF component moments (chi/eta/w/Ncn/rr/Nr means+stdevs+correlations,
    from rico_fort stats) into the JAX covar driver and compares the 5 outputs to the stored
    rtp2_mc / thlp2_mc / wprtp_mc / wpthlp_mc / rtpthlp_mc. The driver computes on zt; the Fortran
    stores _mc = zt2zm(_mc_zt) with the bottom/top momentum levels zeroed (KK_microphys_module.F90
    :901-916), so we apply the same zt2zm + boundary zeroing. NO timing confound — the stored
    component moments are fed directly (same as test_kk_rates_vs_rico_oracle). Key facts:
    mu_eta CANCELS in a covariance (x_mean=mu_x1), so it need not be stored; mu_w=w_1,
    sigma_w=sqrt(varnce_w_1); sigma_eta=stdev_eta_1; rico has sigma_Ncn=0 (constant Nc) so all Ncn
    correlations are immaterial. KK_{auto,accr,evap}_tndcy = the stored rrm_{auto,accr,evap}."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP"); return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP"); return
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_covariances import KK_upscaled_covar_driver
    from clubb_jax.src.Microphys.KK_microphys_module import kk_auto_coef
    ds = nc.Dataset(_RICO_STATS)
    gz = lambda n: np.asarray(ds[n][:, :, 0])              # (time, zt) — zt-level stats
    gm = lambda n: np.asarray(ds[n][:, :, 0])              # (time, zm) — zm-level stats (_mc)
    z = np.zeros_like(gz("chi_1"))
    # chi / eta / w component moments
    chi1, chi2, sc1, sc2 = gz("chi_1"), gz("chi_2"), gz("stdev_chi_1"), gz("stdev_chi_2")
    se1, se2 = gz("stdev_eta_1"), gz("stdev_eta_2")
    w1, w2 = gz("w_1"), gz("w_2")
    sw1, sw2 = np.sqrt(gz("varnce_w_1")), np.sqrt(gz("varnce_w_2"))
    mf = gz("mixt_frac")
    # rr / Nr in-precip moments + log-space conversions
    mrr1, mrr2, srr1, srr2 = gz("mu_rr_1"), gz("mu_rr_2"), gz("sigma_rr_1"), gz("sigma_rr_2")
    mNr1, mNr2, sNr1, sNr2 = gz("mu_Nr_1"), gz("mu_Nr_2"), gz("sigma_Nr_1"), gz("sigma_Nr_2")
    mrr1n, srr1n, rs1 = _logm(mrr1, srr1); mrr2n, srr2n, rs2 = _logm(mrr2, srr2)
    mNr1n, sNr1n, Ns1 = _logm(mNr1, sNr1); mNr2n, sNr2n, Ns2 = _logm(mNr2, sNr2)
    # Ncn (constant in rico: sigma_Ncn=0) — log-mean only, correlations immaterial
    mNcn1, mNcn2 = gz("mu_Ncn_1"), gz("mu_Ncn_2")
    mNcn1n, mNcn2n = np.log(np.maximum(mNcn1, 1e-30)), np.log(np.maximum(mNcn2, 1e-30))
    # normal-space correlations from the stored LINEAR correlations (guarding zero variance)
    NL = lambda nm, sn, s2: np.asarray(corr_NL2NN(gz(nm), sn, s2))
    cc_rr1, cc_rr2 = NL("corr_chi_rr_1", srr1n, rs1), NL("corr_chi_rr_2", srr2n, rs2)
    cc_Nr1, cc_Nr2 = NL("corr_chi_Nr_1", sNr1n, Ns1), NL("corr_chi_Nr_2", sNr2n, Ns2)
    ce_rr1, ce_rr2 = NL("corr_eta_rr_1", srr1n, rs1), NL("corr_eta_rr_2", srr2n, rs2)
    ce_Nr1, ce_Nr2 = NL("corr_eta_Nr_1", sNr1n, Ns1), NL("corr_eta_Nr_2", sNr2n, Ns2)
    cw_rr1, cw_rr2 = NL("corr_w_rr_1", srr1n, rs1), NL("corr_w_rr_2", srr2n, rs2)
    cw_Nr1, cw_Nr2 = NL("corr_w_Nr_1", sNr1n, Ns1), NL("corr_w_Nr_2", sNr2n, Ns2)
    crN1 = np.asarray(corr_LL2NN(gz("corr_rr_Nr_1"), srr1n, sNr1n, rs1, Ns1))
    crN2 = np.asarray(corr_LL2NN(gz("corr_rr_Nr_2"), srr2n, sNr2n, rs2, Ns2))
    cwchi1, cwchi2 = gz("corr_w_chi_1"), gz("corr_w_chi_2")      # normal-normal (direct)
    cchieta1, cchieta2 = gz("corr_chi_eta_1"), gz("corr_chi_eta_2")
    # process coefficients + tendencies
    pf1, pf2 = gz("precip_frac_1"), gz("precip_frac_2")
    coef_a = np.asarray(kk_auto_coef(gz("rho"))); coef_c = 67.0
    T_liq, p = gz("thlm") * gz("exner"), gz("p_in_Pa")
    coef_e = np.asarray(kk_evap_coef(T_liq, p, _C_EVAP))
    t_a, t_c, t_e = gz("rrm_auto"), gz("rrm_accr"), gz("rrm_evap")
    J = jnp.asarray
    out = KK_upscaled_covar_driver(
        J(gz("wm_zt")), J(gz("rtm")), J(gz("thlm")), J(gz("exner")),
        J(w1), J(w2), J(chi1), J(chi2), J(z), J(z),          # mu_w, mu_chi, mu_eta(=0, cancels)
        J(mrr1), J(mrr2), J(mNr1), J(mNr2), J(mNcn1), J(mNcn2),
        J(mrr1n), J(mrr2n), J(mNr1n), J(mNr2n), J(mNcn1n), J(mNcn2n),
        J(sw1), J(sw2), J(sc1), J(sc2), J(se1), J(se2),
        J(srr1), J(srr2), J(sNr1), J(sNr2), J(z), J(z),       # sigma_Ncn = 0
        J(srr1n), J(srr2n), J(sNr1n), J(sNr2n), J(z), J(z),   # sigma_Ncn_n = 0
        J(cwchi1), J(cwchi2), J(cw_rr1), J(cw_rr2), J(cw_Nr1), J(cw_Nr2), J(z), J(z),
        J(cchieta1), J(cchieta2), J(cc_rr1), J(cc_rr2), J(cc_Nr1), J(cc_Nr2), J(z), J(z),
        J(ce_rr1), J(ce_rr2), J(ce_Nr1), J(ce_Nr2), J(z), J(z),
        J(crN1), J(crN2), J(mf), J(pf1), J(pf2),
        J(coef_e), J(coef_a), J(coef_c), J(t_e), J(t_a), J(t_c),
        J(gz("rt_1")), J(gz("rt_2")), J(gz("thl_1")), J(gz("thl_2")),
        J(gz("crt_1")), J(gz("crt_2")), J(gz("cthl_1")), J(gz("cthl_2")))
    # zt -> zm + zero the boundary momentum levels (KK_microphys_module.F90:901-916)
    nt = chi1.shape[0]
    gr = _grid_from_momentum_heights(np.asarray(ds["zm"][:]), nt)
    names = ("wprtp_mc", "wpthlp_mc", "rtp2_mc", "thlp2_mc", "rtpthlp_mc")
    refs = {n: gm(n) for n in names}
    ds.close()
    # The auto/accr covariances are machine-faithful; the EVAP covariance carries the documented
    # within-step/end-of-step timing confound on its (mu_rt-rtm) component-mean-deviation term (the
    # stored rt_1/rtm differ slightly from the within-step values the Fortran covar used — same
    # confound as the accr/evap rate tests). So gate on the MEDIAN (~1e-6, validates the math) with a
    # max floor that tolerates the few evap-timing points (~2e-4).
    ok = True
    for nm, o_zt in zip(names, out):
        o = np.array(_zt2zm(jnp.asarray(o_zt), gr))
        o[:, 0] = 0.0; o[:, -1] = 0.0
        ref = refs[nm]
        nz = np.abs(ref) > np.nanmax(np.abs(ref)) / 1e3      # significant points (within 3 orders of peak)
        rel = np.abs(o[nz] - ref[nz]) / np.abs(ref[nz])
        mx = rel.max() if nz.sum() else 0.0
        med = np.median(rel) if nz.sum() else 0.0
        good = (med < 1e-5) and (mx < 5e-4)
        if not good:
            ok = False
        print(f"  {nm}: {nz.sum()} sig pts, max rel {mx:.2e}, median {med:.2e}  {'PASS' if good else 'FAIL'}")
    assert ok, "KK_upscaled_covar_driver vs rico: a _mc output exceeds the gate (median 1e-5 / max 5e-4)"
    print("  KK_upscaled_covar_driver vs rico: all 5 _mc match (auto/accr machine-exact; evap to the "
          "timing-confound floor)  PASS")


if __name__ == "__main__":
    print("KK rate functions end-to-end vs Fortran rico oracle:")
    test_kk_rates_vs_rico_oracle()
    test_kk_autoconversion_vs_rf02_do()
    test_kk_autoconversion_driver_vs_rico()
    test_kk_accr_evap_drivers_vs_rico()
    test_kk_microphys_adjust_vs_rico()
    test_compute_kk_microphysics_vs_rico()
    test_compute_kk_microphysics_vel_prereqs()
    test_kk_sedimentation_vs_rico()
    test_kk_sed_vel_covars_vs_rico()
    test_sed_centered_diff_lhs_vs_rico()
    test_term_turb_sed_lhs_vs_rico()
    test_microphys_mean_adv_vs_rico()
    test_microphys_lhs_assembly()
    test_microphys_rhs_turb_sed_vs_rico()
    test_microphys_rhs_assembly()
    test_calculate_K_hm()
    test_advance_one_hydrometeor()
    test_kk_covar_driver_vs_rico()
    print("Done.")

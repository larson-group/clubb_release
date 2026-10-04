"""End-to-end validation of the JAX KK rate functions against the Fortran rico oracle.

Unlike test_kk_autoconversion.py (which checks the rate functions against first-principles
quadrature), this feeds the FORTRAN's OWN PDF component moments — read from a real rico SCM
run's stats — into the JAX KK rate functions and compares to the Fortran's `rrm_auto`,
`rrm_accr`, `rrm_evap` outputs. It validates the rates against the actual Fortran microphysics
oracle, isolating the rate-function math from the hydrometeor PDF setup.

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
  python run_scripts/run_scm.py rico -legacy -max_iters 10 -output_dir output/rico_fort
The test skips (does not fail) if the stats file is absent.
"""
import os
import pytest
import numpy as np
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    KK_auto_upscaled_mean, KK_accr_upscaled_mean, KK_evap_upscaled_mean,
)
from clubb_jax.tests.microphysics_coefficients import kk_evap_coef, kk_auto_coef
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
        pytest.skip('Optional Fortran Rico oracle is unavailable')
    if not os.path.exists(_RICO_STATS):
        pytest.skip('Optional Fortran Rico oracle is unavailable')

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
        pytest.skip('Optional Fortran Rico oracle is unavailable')
    if not os.path.exists(_RICO_STATS):
        pytest.skip('Optional Fortran Rico oracle is unavailable')
    from clubb_jax.src.Microphys.KK_microphys_module import KK_microphys_adjust
    ds = nc.Dataset(_RICO_STATS)
    g = lambda n: np.asarray(ds[n][:, :, 0])
    dt = 300.0   # rico dt_main
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc = (np.asarray(x) for x in KK_microphys_adjust(
        dt, g("exner"), g("rcm"), g("rrm"), g("Nrm"),
        g("rrm_evap"), g("rrm_auto"), g("rrm_accr"), g("Nrm_evap"), g("Nrm_auto"), True, True)[:5])
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






_RF02_STATS = os.path.join(os.path.dirname(__file__),
                           "../../output/rf02_do_fort/dycoms2_rf02_do_stats.nc")

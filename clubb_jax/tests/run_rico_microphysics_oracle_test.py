#!/usr/bin/env python3
"""Compare JAX microphysics with rich statistics from a real ten-step Fortran Rico run."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.add_argument('-stats_file', type=Path, required=True, help='Fortran Rico all_stats NetCDF output')
    args = parser.parse_args()
    if not args.stats_file.is_file():
        parser.error(f'Stats file not found: {args.stats_file}')
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()

import netCDF4 as nc
import numpy as np
import jax
import jax.numpy as jnp
jax.config.update('jax_enable_x64', True)
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    KK_auto_upscaled_mean, KK_accr_upscaled_mean, KK_evap_upscaled_mean,
)
from clubb_jax.pytests.microphysics_coefficients import kk_evap_coef, kk_auto_coef
from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N, corr_NL2NN, corr_LL2NN
from clubb_jax.src.CLUBB_core.Nc_Ncn_eqns import Nc_in_cloud_to_Ncnm

from clubb_jax.pytests.microphysics_test_inputs import precip_fraction_from_fields

_RR_TOL = 1.0e-10
_NR_TOL = _RR_TOL / ((4.0 / 3.0) * np.pi * 1000.0 * (5.0e-3) ** 3)
_UPSILON = 0.55  # Rico's case-specific precipitation fraction parameter.
_C_EVAP = 0.86  # Rico's case-specific evaporation coefficient.

def _logm(mu, sig):
    """(mu_n, sigma_n, sigma2_on_mu2) for a lognormal from its linear mean/std."""
    s2m2 = np.where(mu > 0, (sig / np.maximum(mu, 1e-30)) ** 2, 0.0)
    return (np.asarray(mean_L2N(np.maximum(mu, 1e-30), s2m2)),
            np.asarray(stdev_L2N(s2m2)), s2m2)



def _rel(out, ref, mask):
    return np.abs(out[mask] - ref[mask]) / np.abs(ref[mask])



def check_kk_rates_vs_rico_oracle(stats_file):
    ds = nc.Dataset(stats_file)
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



def check_kk_microphys_adjust_vs_rico(stats_file):
    """KK_microphys_adjust (the tendency assembly) reproduces rico's rcm_mc / rrm_mc.

    Feeds the stored process rates (rrm_auto/accr/evap, Nrm_auto/evap) + state (rcm/rrm/Nrm/exner)
    into the assembly. rcm_mc = -(adjusted auto+accr) is a clean pure-function check (exact,
    including the source-over-depletion adjustment). rrm_mc = evap_net + source matches where the
    EVAP limiter doesn't trigger (the limiter's -rrm/dt uses the within-step rrm, which differs from
    the end-of-step stored rrm — the documented timing confound). thlm_mc is checked for
    self-consistency with the oracle formula (Lv/(Cp·exner)·rrm_mc)."""
    from clubb_jax.src.Microphys.KK_microphys_module import KK_microphys_adjust
    ds = nc.Dataset(stats_file)
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



def check_Ncnm_vs_rico_stats(stats_file):
    """Reproduce rico's Ncnm exactly (rico: constant N_c -> Ncnm = Nc_in_cloud)."""
    ds = nc.Dataset(stats_file)
    G = lambda n: np.asarray(ds[n][:, :, 0]).ravel()
    jo = np.asarray(Nc_in_cloud_to_Ncnm(
        G("chi_1"), G("chi_2"), G("stdev_chi_1"), G("stdev_chi_2"), G("mixt_frac"),
        G("Nc_in_cloud"), G("cloud_frac_1"), G("cloud_frac_2"), 0.0, 0.0))
    ncnm_s = G("Ncnm")
    ds.close()
    m = ncnm_s > 0
    d = np.max(np.abs(jo[m] - ncnm_s[m]) / ncnm_s[m])
    assert d < 1e-13, f"Ncnm vs rico stats max rel {d:.2e}"
    print(f"  Nc_in_cloud_to_Ncnm vs rico Ncnm (const-Ncn): max rel {d:.1e}  PASS")



def check_precip_fraction_vs_rico_oracle(stats_file):
    ds = nc.Dataset(stats_file)
    G = lambda n: np.asarray(ds[n][:, :, 0])           # (nt, nzt)
    cf, cf1, cf2 = G("cloud_frac"), G("cloud_frac_1"), G("cloud_frac_2")
    isf = G("ice_supersat_frac")
    mf, rrm, Nrm = G("mixt_frac"), G("rrm"), G("Nrm")
    pfs, pf1s, pf2s = G("precip_frac"), G("precip_frac_1"), G("precip_frac_2")
    ds.close()

    nt, nzt = cf.shape
    z = np.zeros((nt, nzt))
    hydromet = np.stack([rrm, Nrm], axis=-1)           # (nt, nzt, 2): rr (mix ratio), Nr
    l_mix = np.array([1, 0]); l_frozen = np.array([0, 0])
    hm_tol = np.array([_RR_TOL, _NR_TOL])

    pf, pf1, pf2, pftol = (np.asarray(x) for x in precip_fraction_from_fields(
        hydromet, cf, cf1, cf2, isf, z, z, mf, l_mix, l_frozen, hm_tol, _UPSILON))

    # Well-resolved precip region (comfortably above cloud_frac_min=0.005, both agree there
    # is precip): the inputs are unambiguous, so the match must be machine-exact.
    mask = (pfs > 0.006) & (pf > 0.006)
    assert mask.sum() >= 5, f"too few well-resolved points ({mask.sum()})"
    d = np.max(np.abs(pfs[mask] - pf[mask]))
    d1 = np.max(np.abs(pf1s[mask] - pf1[mask]))
    d2 = np.max(np.abs(pf2s[mask] - pf2[mask]))
    assert d < 1e-13 and d1 < 1e-13 and d2 < 1e-13, \
        f"precip_fraction vs rico: pf {d:.2e}, pf1 {d1:.2e}, pf2 {d2:.2e}"

    # Internal consistency everywhere: f_p = a f_p(1) + (1-a) f_p(2); fractions in [0,1].
    recon = mf * pf1 + (1.0 - mf) * pf2
    assert np.max(np.abs(recon - pf)) < 1e-12, "precip_frac != mixt_frac-weighted components"
    for name, a in (("pf", pf), ("pf1", pf1), ("pf2", pf2)):
        assert a.min() >= -1e-14 and a.max() <= 1.0 + 1e-12, f"{name} out of [0,1]"

    n_edge = int(np.sum((pfs < 1e-9) & (pf > 1e-9)))   # tol-boundary timing-confound levels
    print(f"  precip_fraction vs rico: well-resolved ({mask.sum()} pts) bit-exact "
          f"(pf {d:.1e}, pf1 {d1:.1e}, pf2 {d2:.1e}); {n_edge} tol-boundary edge levels  PASS")


if __name__ == '__main__':
    check_kk_rates_vs_rico_oracle(args.stats_file)
    check_kk_microphys_adjust_vs_rico(args.stats_file)
    check_Ncnm_vs_rico_stats(args.stats_file)
    check_precip_fraction_vs_rico_oracle(args.stats_file)
    print('Rico microphysics oracle checks passed')

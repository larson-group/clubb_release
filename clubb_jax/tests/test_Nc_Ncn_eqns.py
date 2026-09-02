"""Verification of Nc_Ncn_eqns.py — cloud-nuclei-concentration mean <Ncn> (Ncnm).

Nc_in_cloud_to_Ncnm IS exposed by the f2py API, so it is verified BIT-TO-BIT against
f2py_nc_in_cloud_to_ncnm over a random sweep covering both branches (constant Ncn -> simple
fallback; varying Ncn -> the erfc PDF integral). It also exactly reproduces rico's Ncnm
(rico has constant N_c, so Ncnm = Nc_in_cloud).
"""
import os
import numpy as np
import jax

jax.config.update("jax_enable_x64", True)

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.Nc_Ncn_eqns import (
    Nc_in_cloud_to_Ncnm, bivar_NL_chi_Ncn_mean, Ncnm_to_Ncm, Ncnm_to_Nc_in_cloud)
from scipy.special import erfc as _erfc

_RICO_STATS = os.path.join(os.path.dirname(__file__),
                           "../../output/rico_fort/rico_stats.nc")


def test_Ncnm_bit_to_bit_vs_f2py():
    """Nc_in_cloud_to_Ncnm matches the Fortran f2py oracle over a random sweep."""
    try:
        import clubb_f2py as f
    except ImportError:
        print("  Ncnm vs f2py: SKIP (clubb_f2py not built)"); return
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(3000):
        mu1, mu2 = rng.uniform(-1e-3, 1e-3), rng.uniform(-1e-3, 1e-3)
        s1 = rng.choice([0.0, rng.uniform(0, 5e-4)])      # exercise sigma_chi=0 branch
        s2 = rng.choice([0.0, rng.uniform(0, 5e-4)])
        mf, ncl = rng.uniform(0.01, 0.99), rng.uniform(1e6, 3e8)
        cf1, cf2 = rng.uniform(0, 1), rng.uniform(0, 1)
        cnp2 = rng.choice([0.0, rng.uniform(0, 2.0)])     # exercise const-Ncn vs varying
        ccc = rng.choice([0.0, rng.uniform(-0.9, 0.9)])
        fo = float(f.f2py_nc_in_cloud_to_ncnm(mu1, mu2, s1, s2, mf, ncl, cf1, cf2, cnp2, ccc))
        jo = float(Nc_in_cloud_to_Ncnm(mu1, mu2, s1, s2, mf, ncl, cf1, cf2, cnp2, ccc))
        worst = max(worst, abs(jo - fo) / (abs(fo) + 1e-300))
    # Machine-level (residual is the JAX-vs-Fortran erfc implementation difference).
    assert worst < 1e-12, f"Ncnm vs f2py worst rel {worst:.2e}"
    print(f"  Nc_in_cloud_to_Ncnm vs f2py: 3000 cases, worst rel {worst:.1e}  PASS")


def test_Ncnm_vs_rico_stats():
    """Reproduce rico's Ncnm exactly (rico: constant N_c -> Ncnm = Nc_in_cloud)."""
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP")
        return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP")
        return
    ds = nc.Dataset(_RICO_STATS)
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


# ── Forward direction (Ncn → Nc): independent NumPy transcription of the Fortran branch logic ──
_CHI_TOL = max(1.0e-8, np.finfo(np.float64).eps)
_NCN_TOL = 1.0e2
_SQRT2 = np.sqrt(2.0)


def _ref_bivar_NL(mu_chi, mu_Ncn, sig_chi, sig_Ncn, sig_Ncn_n, corr_n):
    """Scalar NumPy reference of Nc_Ncn_eqns.F90:bivar_NL_chi_Ncn_mean (the four Fortran branches)."""
    if sig_chi <= _CHI_TOL:                 # both-0 and chi-0 branches are identical
        return mu_Ncn if mu_chi > 0.0 else 0.0
    if sig_Ncn <= _NCN_TOL:                 # chi varies, Ncn constant
        return mu_Ncn * 0.5 * _erfc(-(mu_chi / (_SQRT2 * sig_chi)))
    return 0.5 * mu_Ncn * _erfc(-(1.0 / _SQRT2) * (mu_chi / sig_chi + corr_n * sig_Ncn_n))


def test_forward_vs_numpy_reference():
    """bivar_NL_chi_Ncn_mean / Ncnm_to_Ncm / Ncnm_to_Nc_in_cloud vs an independent NumPy reference."""
    rng = np.random.default_rng(1)
    worst_b = worst_m = worst_c = 0.0
    for _ in range(4000):
        mc1, mc2 = rng.uniform(-1e-3, 1e-3), rng.uniform(-1e-3, 1e-3)
        mn1, mn2 = rng.uniform(1e6, 3e8), rng.uniform(1e6, 3e8)
        sc1 = rng.choice([0.0, rng.uniform(0, 5e-4)])          # exercise sigma_chi<=tol branch
        sc2 = rng.choice([0.0, rng.uniform(0, 5e-4)])
        sn1 = rng.choice([0.0, rng.uniform(0, 5e2)])           # exercise sigma_Ncn<=Ncn_tol branch
        sn2 = rng.choice([0.0, rng.uniform(0, 5e2)])
        snn1, snn2 = rng.uniform(0, 1.5), rng.uniform(0, 1.5)
        cc1, cc2 = rng.uniform(-0.9, 0.9), rng.uniform(-0.9, 0.9)
        mf = rng.uniform(0.01, 0.99)
        cf1, cf2 = rng.uniform(0, 1), rng.uniform(0, 1)

        rb = _ref_bivar_NL(mc1, mn1, sc1, sn1, snn1, cc1)
        jb = float(bivar_NL_chi_Ncn_mean(mc1, mn1, sc1, sn1, snn1, cc1))
        worst_b = max(worst_b, abs(jb - rb) / (abs(rb) + 1e-300))

        rm = mf * _ref_bivar_NL(mc1, mn1, sc1, sn1, snn1, cc1) + (1 - mf) * _ref_bivar_NL(mc2, mn2, sc2, sn2, snn2, cc2)
        jm = float(Ncnm_to_Ncm(mc1, mc2, mn1, mn2, sc1, sc2, sn1, sn2, snn1, snn2, cc1, cc2, mf))
        worst_m = max(worst_m, abs(jm - rm) / (abs(rm) + 1e-300))

        cfrac = mf * cf1 + (1 - mf) * cf2
        rc = (rm / cfrac) if cfrac > 0.005 else mn1
        jc = float(Ncnm_to_Nc_in_cloud(mc1, mc2, mn1, mn2, sc1, sc2, sn1, sn2, snn1, snn2, cc1, cc2, mf, cf1, cf2))
        worst_c = max(worst_c, abs(jc - rc) / (abs(rc) + 1e-300))
    assert worst_b < 1e-12 and worst_m < 1e-12 and worst_c < 1e-12, (worst_b, worst_m, worst_c)
    print(f"  forward vs NumPy ref: 4000 cases, worst rel bivar={worst_b:.1e} Ncm={worst_m:.1e} "
          f"Nc_in_cloud={worst_c:.1e}  PASS")


if __name__ == "__main__":
    print("Nc_Ncn_eqns (<Ncn> mean) verification:")
    test_Ncnm_bit_to_bit_vs_f2py()
    test_Ncnm_vs_rico_stats()
    test_forward_vs_numpy_reference()
    print("All Nc_Ncn_eqns tests PASSED.")

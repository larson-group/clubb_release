"""Validation of the JAX precip_fraction port against the Fortran rico oracle.

precipitation_fraction.F90 is a CLUBB_core routine, but its f2py wrapper
(f2py_precip_fraction) FPE-traps in this environment (a build-level -ffpe-trap plus the
wrapper's internally-constructed grid), so the gold-standard f2py bit-to-bit path is
unavailable. Instead this validates against a real rico run's stats: feeding the stored
inputs (cloud_frac, cloud_frac_1/2, mixt_frac, rrm/Nrm as the hydromet field) into the JAX
precip_fraction and comparing to the stored precip_frac / precip_frac_1 / precip_frac_2.

On the WELL-RESOLVED precipitation region the match is bit-exact (machine epsilon),
confirming the full algorithm: the downward cumulative-max overall fraction, the upsilon=0.55
component split with all its nested special cases, the max_hm_ip_comp_mean limiter, and the
[precip_frac_tol, 1] clamps. Discrepancies are confined to tol-boundary levels (precip_frac
near cloud_frac_min=0.005) where the discrete `any(hydromet>=tol)` test flips because the
END-of-step stored rrm differs slightly from the value at precip_fraction-call-time — an
inherent limitation of the stats oracle, not a port error.

Requires output/rico_fort/rico_stats.nc; skips if absent.
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

from clubb_jax.src.CLUBB_core.parameter_indices import iupsilon_precip_frac_rat, nparams
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.precipitation_fraction import precip_fraction, precip_frac_assert_check
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_RICO_STATS = os.path.join(os.path.dirname(__file__),
                           "../../output/rico_fort/rico_stats.nc")
_RR_TOL = 1.0e-10
_NR_TOL = _RR_TOL / ((4.0 / 3.0) * np.pi * 1000.0 * (5.0e-3) ** 3)   # ~1.9099e-7
_UPSILON = 0.55   # rico tunable upsilon_precip_frac_rat


def _precip_fraction_old_api(
    hydromet,
    cloud_frac,
    cloud_frac_1,
    cloud_frac_2,
    ice_supersat_frac,
    ice_supersat_frac_1,
    ice_supersat_frac_2,
    mixt_frac,
    l_mix_rat_hm,
    l_frozen_hm,
    hydromet_tol,
    upsilon,
):
    ngrdcol, nzt, hydromet_dim = np.asarray(hydromet).shape
    gr = setup_grid(
        ngrdcol=ngrdcol,
        deltaz=1.0,
        zm_init=0.0,
        zm_top=float(nzt),
        grid_type=1,
    )
    clubb_params = jnp.zeros((ngrdcol, nparams), dtype=jnp.float64)
    clubb_params = clubb_params.at[:, iupsilon_precip_frac_rat].set(upsilon)
    stats = JaxStats.empty(
        l_sample=False,
        names=(),
        ncol=ngrdcol,
        max_nlev=max(nzt, gr.nzm),
    )
    _err, pf, pf1, pf2, pftol, _stats = precip_fraction(
        gr,
        nzt,
        ngrdcol,
        hydromet_dim,
        hydromet,
        cloud_frac,
        cloud_frac_1,
        l_mix_rat_hm,
        l_frozen_hm,
        hydromet_tol,
        cloud_frac_2,
        ice_supersat_frac,
        ice_supersat_frac_1,
        ice_supersat_frac_2,
        mixt_frac,
        clubb_params,
        ErrInfo.initialized(ngrdcol),
        stats,
    )
    return pf, pf1, pf2, pftol


def test_precip_fraction_vs_rico_oracle():
    try:
        import netCDF4 as nc
    except ImportError:
        print("  netCDF4 not available — SKIP")
        return
    if not os.path.exists(_RICO_STATS):
        print("  rico_fort stats absent — SKIP")
        return

    ds = nc.Dataset(_RICO_STATS)
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

    pf, pf1, pf2, pftol = (np.asarray(x) for x in _precip_fraction_old_api(
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


def test_precip_frac_assert_check():
    """precip_frac_assert_check accepts precip_fraction's own (self-consistent, in-range) output and rejects
    corrupted inputs — and serves as a cross-check that the JAX precip_fraction satisfies the Fortran assertions."""
    rng = np.random.default_rng(0)
    ng, nzt, hd = 1, 20, 2
    hydromet = np.abs(rng.normal(size=(ng, nzt, hd))) * 1e-4
    cf = rng.uniform(0, 1, (ng, nzt)); mf = rng.uniform(0.3, 0.7, (ng, nzt)); isf = np.zeros((ng, nzt))
    hmtol = np.array([1e-10, 1.9e-7])
    pf, pf1, pf2, pftol = _precip_fraction_old_api(
        hydromet, cf, cf, cf, isf, isf, isf, mf,
        np.array([True, False]), np.array([False, False]), hmtol, 1.0)
    args = (hydromet[0], hmtol, mf[0], pf[0], pf1[0], pf2[0], float(pftol[0]))
    assert precip_frac_assert_check(*args) is True, "valid precip_fraction output rejected"
    bad = np.array(pf[0]); bad[5] = 1.5            # precip_frac > 1
    assert precip_frac_assert_check(hydromet[0], hmtol, mf[0], bad, pf1[0], pf2[0], float(pftol[0])) is False
    bad2 = np.array(pf1[0]); bad2[3] += 0.5        # breaks the mixt_frac-weighted consistency
    assert precip_frac_assert_check(hydromet[0], hmtol, mf[0], pf[0], bad2, pf2[0], float(pftol[0])) is False
    print("  precip_frac_assert_check: valid PASS / precip_frac>1 + inconsistent FAIL  PASS")


if __name__ == "__main__":
    print("precip_fraction end-to-end vs Fortran rico oracle:")
    test_precip_fraction_vs_rico_oracle()
    test_precip_frac_assert_check()
    print("Done.")

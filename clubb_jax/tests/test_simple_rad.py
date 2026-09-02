#!/usr/bin/env python3
"""test_simple_rad.py — validate Radiation/simple_rad_module.py:simple_rad (clean l_rad_above_cloud=False path).

`simple_rad` (simple_rad_module.F90) is the simplified cloud-top-cooling longwave radiation. dycoms2_rf01 (a
bit-faithful DEFAULT_CASE) runs it with `l_rad_above_cloud=.false.` — the clean path validated here end-to-end vs
an independent transcription of the F90 formulas:

    LWP        = liq_water_path (cumulative from top)
    Frad_LW(k) = F0·exp(−κ·LWP(k)) + F1·exp(−κ·(LWP_bottom − LWP(k)))     [F1>0; else just the F0 term]   (F90:255)
    radht_LW(k)= (1/exner(k))·(−1/(Cp·rho(k)))·(Frad_LW(k+1) − Frad_LW(k))·invrs_dzt(k)                    (F90:348)

covering both the F1>0 (cloud-base flux) and F1=0 branches. The `l_rad_above_cloud=.true.` correction (dycoms2_rf02_*)
is the remaining complex piece (inversion-height + dz^(4/3) term), case-validated. Pure-Python (the F90 formulas are
the oracle), so it never SKIPs. (iter 505)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp

from clubb_jax.src.Radiation.simple_rad_module import simple_rad, simple_rad_bomex
from clubb_jax.src.Radiation.parameters_radiation import RadiationParameters
from clubb_jax.src.CLUBB_core.constants_clubb import Cp as _CP
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_EPS = np.finfo(np.float64).eps


def _radiation_parameters(f0, f1, kappa, l_rad_above_cloud):
    zeros_20 = jnp.zeros((20,), dtype=jnp.float64)
    return RadiationParameters(
        "simplified", 0.0, kappa, f0, f1, 1.0e-5, 0.85, 0.999,
        zeros_20, zeros_20, zeros_20,
        False, False, l_rad_above_cloud, 1, False,
        jnp.zeros((33,), dtype=jnp.float64), jnp.zeros((33, 36), dtype=jnp.float64),
    )


def _stats():
    return JaxStats.empty(l_sample=False, names=(), ncol=_NG, max_nlev=1)


def _ref(rho, rcm, exner, invrs_dzt, f0, f1, kappa):
    ng, nzt = rcm.shape
    nzm = nzt + 1
    contrib = rcm * rho / invrs_dzt
    lwp = np.zeros((ng, nzm))
    for k in range(nzt - 1, -1, -1):
        lwp[:, k] = lwp[:, k + 1] + contrib[:, k]
    if f1 > np.finfo(np.float64).eps:
        frad = f0 * np.exp(-kappa * lwp) + f1 * np.exp(-kappa * (lwp[:, 0:1] - lwp))
    else:
        frad = f0 * np.exp(-kappa * lwp)
    radht = (1.0 / exner) * (-1.0 / (_CP * rho)) * (frad[:, 1:] - frad[:, :-1]) * invrs_dzt
    return frad, radht


def test_simple_rad_clean():
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    invrs_dzt = np.asarray(gr.invrs_dzt)
    rng = np.random.default_rng(55)
    wf = wr = 0.0
    for f1 in (74.0, 0.0):           # F1>0 (both fluxes) and F1=0 (cloud-top only) branches
        for _ in range(8):
            rho = rng.uniform(0.5, 1.2, (_NG, nzt))
            rcm = np.zeros((_NG, nzt)); rcm[:, 6:16] = rng.uniform(1e-5, 1e-3, (_NG, 10))
            exner = rng.uniform(0.6, 1.0, (_NG, nzt))
            rho_zm = rng.uniform(0.5, 1.2, (_NG, nzm))    # unused on the clean path
            rtm = rng.uniform(1e-3, 1.5e-2, (_NG, nzt))   # unused on the clean path
            f0, kappa = 70.0, 130.0
            _, frad_j, radht_j = simple_rad(
                gr, _NG, rho, rho_zm, rtm, rcm, exner,
                _stats(), _radiation_parameters(f0, f1, kappa, False),
            )
            frad_r, radht_r = _ref(rho, rcm, exner, invrs_dzt, f0, f1, kappa)
            wf = max(wf, float(np.max(np.abs(np.asarray(frad_j) - frad_r))))
            wr = max(wr, float(np.max(np.abs(np.asarray(radht_j) - radht_r))))
    assert wf < 1e-10 and wr < 1e-12, f"simple_rad mismatch: Frad_LW {wf:.2e}, radht_LW {wr:.2e}"
    print(f"  simple_rad (l_rad_above_cloud=False, F1>0 + F1=0): Frad_LW + radht_LW match, worst {max(wf, wr):.2e}  PASS")


def test_simple_rad_above_cloud():
    """The l_rad_above_cloud=.true. correction (dycoms2_rf02_*): inversion height z_i (first level rtm<=8e-3 from
    the bottom, interpolated) then Frad_LW += rho_zm·Cp·Ls_div·Heaviside(zm-z_i)·(0.25·dz^(4/3) + z_i·dz^(1/3)).
    Uses a clear moist-below / dry-above rtm profile so k_iso is interior (the case regime)."""
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = gr.zm.shape[1]; nzt = nzm - 1
    invrs_dzt = np.asarray(gr.invrs_dzt); zt = np.asarray(gr.zt); zm = np.asarray(gr.zm)
    rng = np.random.default_rng(56)
    f0, f1, kappa = 70.0, 74.0, 130.0
    worst = 0.0
    for _ in range(8):
        rho = rng.uniform(0.5, 1.2, (_NG, nzt)); rho_zm = rng.uniform(0.5, 1.2, (_NG, nzm))
        exner = rng.uniform(0.6, 1.0, (_NG, nzt))
        rcm = np.zeros((_NG, nzt)); rcm[:, 6:14] = rng.uniform(1e-5, 1e-3, (_NG, 8))
        # rtm: well-mixed moist (~12 g/kg) up to k=14, then a sharp drop below 8 g/kg (a clear inversion)
        rtm = np.empty((_NG, nzt))
        rtm[:, :14] = rng.uniform(1.0e-2, 1.4e-2, (_NG, 14))     # > 8e-3
        rtm[:, 14:] = rng.uniform(1.0e-3, 6.0e-3, (_NG, nzt - 14))  # < 8e-3
        _, frad_j, radht_j = simple_rad(
            gr, _NG, rho, rho_zm, rtm, rcm, exner,
            _stats(), _radiation_parameters(f0, f1, kappa, True),
        )

        # --- transcription ---
        contrib = rcm * rho / invrs_dzt
        lwp = np.zeros((_NG, nzm))
        for k in range(nzt - 1, -1, -1):
            lwp[:, k] = lwp[:, k + 1] + contrib[:, k]
        frad = f0 * np.exp(-kappa * lwp) + f1 * np.exp(-kappa * (lwp[:, 0:1] - lwp))
        z_i = np.zeros(_NG)
        for i in range(_NG):
            k_iso = 0
            while k_iso < nzt and rtm[i, k_iso] > 8.0e-3:
                k_iso += 1
            assert 0 < k_iso < nzt, "test profile must put the inversion at an interior level"
            denom = rtm[i, k_iso] - rtm[i, k_iso - 1]
            z_i[i] = (0.5 * (zt[i, k_iso] + zt[i, k_iso - 1]) if abs(denom) < _EPS
                      else zt[i, k_iso - 1] + (8.0e-3 - rtm[i, k_iso - 1]) * (zt[i, k_iso] - zt[i, k_iso - 1]) / denom)
        dz = zm - z_i[:, None]
        heaviside = np.where(dz < -_EPS, 0.0, np.where(dz > _EPS, 1.0, 0.5))
        dz_pos = np.maximum(dz, 0.0)
        frad = frad + rho_zm * _CP * 3.75e-6 * heaviside * (0.25 * dz_pos ** (4.0 / 3.0) + z_i[:, None] * dz_pos ** (1.0 / 3.0))
        radht = (1.0 / exner) * (-1.0 / (_CP * rho)) * (frad[:, 1:] - frad[:, :-1]) * invrs_dzt
        worst = max(worst, float(np.max(np.abs(np.asarray(frad_j) - frad))),
                    float(np.max(np.abs(np.asarray(radht_j) - radht))))
    assert worst < 1e-9, f"simple_rad above-cloud mismatch {worst:.2e}"
    print(f"  simple_rad (l_rad_above_cloud=True): inversion-height + dz^(4/3) correction match, worst {worst:.2e}  PASS")


def test_simple_rad_bomex():
    """simple_rad_bomex (simple_rad_module.F90:360, used by the bit-faithful bomex case): the analytic piecewise
    radiative-heating profile — −2.315e-5 below 1500 m, ramping to 0 between 1500 and 2500 m, 0 above (F90:393-400)."""
    zt = np.array([[0.0, 500.0, 1499.9, 1500.0, 2000.0, 2499.9, 2500.0, 3000.0]])
    gr = setup_grid(ngrdcol=1, deltaz=500.0, zm_init=0.0, zm_top=4000.0, grid_type=1)
    gr = gr._replace(zt=jnp.asarray(zt))
    j = np.asarray(simple_rad_bomex(gr, 1))
    ref = np.where(zt < 1500.0, -2.315e-5,
                   np.where(zt < 2500.0, -2.315e-5 + 2.315e-5 * (zt - 1500.0) / 1000.0, 0.0))
    assert float(np.max(np.abs(j - ref))) < 1e-20, "simple_rad_bomex != F90 piecewise profile"
    # explicit anchors: below = -2.315e-5; midpoint(2000) = half; above = 0
    assert abs(float(j[0, 1]) - (-2.315e-5)) == 0.0 and float(j[0, 6]) == 0.0
    assert abs(float(j[0, 4]) - (-2.315e-5 + 2.315e-5 * 0.5)) < 1e-20
    print("  simple_rad_bomex == F90 piecewise analytic profile (3 regions)  PASS")


def main():
    print("test_simple_rad:")
    test_simple_rad_clean()
    test_simple_rad_above_cloud()
    test_simple_rad_bomex()
    print("All simple_rad checks PASSED")


if __name__ == "__main__":
    main()

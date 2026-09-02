#!/usr/bin/env python3
"""test_sfc_flux.py — validate Benchmark_cases/sfc_flux.py (shared surface-flux helpers, sfc_flux.F90).

Every case's `*_sfclyr` routine builds its surface fluxes from these shared helpers (mirroring `use sfc_flux`), so
they are case-active but were not directly unit-tested. Pinned here vs the literal F90 formulas:

  convert_sens_ht_to_km_s = sens_ht/(rho_sfc·Cp)              (F90:350)
  convert_latent_ht_to_m_s = latent_ht/(rho_sfc·Lv)           (F90:377)
  compute_wpthlp_sfc = −Cd·ubar·(thlm_sfc − T_sfc/exner_sfc)  (F90:218)
  compute_wprtp_sfc  = −Cd·ubar·(rtm_sfc − adjustment)        (F90:257)
  compute_ubar       = max(0.25, sqrt(um²+vm²))               (F90:104)
  compute_momentum_flux: upwp=−um·ustar²/ubar; vpwp=−vm·ustar²/ubar  (F90:compute_momentum_flux)
  compute_ht_mostr_flux: time-interp of the prescribed schedule (flat extrap outside the table)

Pure-Python (the F90 formulas are the oracle), so it never SKIPs. (iter 508)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    convert_sens_ht_to_km_s, convert_latent_ht_to_m_s, compute_wpthlp_sfc, compute_wprtp_sfc,
    compute_ubar, compute_momentum_flux, compute_ht_mostr_flux, set_sclr_sfc_rtm_thlm,
)
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx


def test_conversions_and_bulk_fluxes():
    rng = np.random.default_rng(58)
    sens = rng.uniform(-50.0, 300.0, 5); latent = rng.uniform(0.0, 400.0, 5); rho = rng.uniform(0.9, 1.25, 5)
    assert float(np.max(np.abs(np.asarray(convert_sens_ht_to_km_s(sens, rho)) - sens / (rho * Cp)))) == 0.0
    assert float(np.max(np.abs(np.asarray(convert_latent_ht_to_m_s(latent, rho)) - latent / (rho * Lv)))) == 0.0

    Cd = rng.uniform(1e-3, 2e-2, 5); ubar = rng.uniform(0.3, 12.0, 5)
    thlm = rng.uniform(285.0, 305.0, 5); T_sfc = rng.uniform(285.0, 305.0, 5); exner = rng.uniform(0.95, 1.0, 5)
    rtm = rng.uniform(1e-3, 2e-2, 5); adj = rng.uniform(1e-3, 2e-2, 5)
    assert float(np.max(np.abs(np.asarray(compute_wpthlp_sfc(2, Cd, ubar, thlm, T_sfc, exner))
                               - (-Cd * ubar * (thlm - T_sfc / exner))))) == 0.0
    assert float(np.max(np.abs(np.asarray(compute_wprtp_sfc(2, Cd, ubar, rtm, adj))
                               - (-Cd * ubar * (rtm - adj))))) == 0.0
    print("  convert_*_ht + compute_wpthlp_sfc/compute_wprtp_sfc == F90 (exact)  PASS")


def test_ubar_and_momentum():
    um = np.array([3.0, 0.1, -4.0, 0.0]); vm = np.array([4.0, 0.05, 3.0, 0.0])
    ub = np.asarray(compute_ubar(2, um, vm))
    ref = np.maximum(0.25, np.sqrt(um ** 2 + vm ** 2))
    assert float(np.max(np.abs(ub - ref))) == 0.0, "compute_ubar != F90"
    assert float(ub[1]) == 0.25 and float(ub[3]) == 0.25, "ubmin floor must apply for near-calm winds"

    ustar = np.array([0.2, 0.3, 0.5, 0.1])
    up, vp = compute_momentum_flux(2, um, vm, ub, ustar)
    assert float(np.max(np.abs(np.asarray(up) - (-um * ustar ** 2 / ub)))) == 0.0
    assert float(np.max(np.abs(np.asarray(vp) - (-vm * ustar ** 2 / ub)))) == 0.0
    print("  compute_ubar (incl. ubmin floor) + compute_momentum_flux == F90 (exact)  PASS")


def test_ht_mostr_flux_interp():
    t_given = np.array([0.0, 100.0, 200.0, 300.0])
    sens = np.array([10.0, 20.0, 0.0, -5.0]); lat = np.array([100.0, 150.0, 120.0, 80.0])
    time_dependent_input.time_sfc_given = jnp.asarray(t_given)
    time_dependent_input.sens_ht_given = jnp.asarray(sens)
    time_dependent_input.latent_ht_given = jnp.asarray(lat)
    # flat extrapolation below/above
    h, m = compute_ht_mostr_flux(-50.0, 4)
    assert (float(h), float(m)) == (10.0, 100.0)
    h, m = compute_ht_mostr_flux(500.0, 4)
    assert (float(h), float(m)) == (-5.0, 80.0)
    # exact node
    h, m = compute_ht_mostr_flux(100.0, 4)
    assert (float(h), float(m)) == (20.0, 150.0)
    # linear interior: t=150 → halfway between node 1 (100) and 2 (200)
    h, m = compute_ht_mostr_flux(150.0, 4)
    assert abs(h - 0.5 * (20.0 + 0.0)) < 1e-12 and abs(m - 0.5 * (150.0 + 120.0)) < 1e-12
    print("  compute_ht_mostr_flux: flat-extrap + linear-interior == F90  PASS")


def test_set_scalar_surface_fluxes_updates_surface_arrays_only():
    state = {
        "sclr_dim": 2,
        "edsclr_dim": 2,
        "sclr_idx": SclrIdx(
            iisclr_rt=2, iisclr_thl=1, iisclr_CO2=-1,
            iiedsclr_rt=1, iiedsclr_thl=2, iiedsclr_CO2=-1,
        ),
        "wpsclrp_sfc": jnp.full((2, 2), 99.0),
        "wpedsclrp_sfc": jnp.full((2, 2), 99.0),
        "wpsclrp": np.full((2, 4, 2), 7.0),
    }
    heat = np.array([1.0, 2.0])
    moisture = np.array([3.0, 4.0])

    state["wpsclrp_sfc"], state["wpedsclrp_sfc"] = set_sclr_sfc_rtm_thlm(
        2,
        state["sclr_dim"],
        state["edsclr_dim"],
        state["sclr_idx"],
        heat,
        moisture,
    )

    np.testing.assert_array_equal(state["wpsclrp_sfc"], np.column_stack((heat, moisture)))
    np.testing.assert_array_equal(state["wpedsclrp_sfc"], np.zeros((2, 2)))
    np.testing.assert_array_equal(state["wpsclrp"], np.full((2, 4, 2), 7.0))


def test_shared_surface_flux_path_jits_and_stays_on_device():
    def surface_kernel(um, vm, ustar, heat, moisture):
        ubar = compute_ubar(2, um, vm)
        upwp, vpwp = compute_momentum_flux(2, um, vm, ubar, ustar)
        return ubar, upwp, vpwp, compute_wprtp_sfc(2, 0.0013, ubar, moisture, heat)

    values = jax.jit(surface_kernel)(
        jnp.array([3.0, 0.1]),
        jnp.array([4.0, 0.05]),
        jnp.array([0.2, 0.3]),
        jnp.array([0.01, 0.02]),
        jnp.array([0.005, 0.006]),
    )
    assert all(isinstance(value, jax.Array) for value in values)
    assert all(bool(jnp.all(jnp.isfinite(value))) for value in values)


def main():
    print("test_sfc_flux:")
    test_conversions_and_bulk_fluxes()
    test_ubar_and_momentum()
    test_ht_mostr_flux_interp()
    test_set_scalar_surface_fluxes_updates_surface_arrays_only()
    test_shared_surface_flux_path_jits_and_stays_on_device()
    print("All sfc_flux checks PASSED")


if __name__ == "__main__":
    main()

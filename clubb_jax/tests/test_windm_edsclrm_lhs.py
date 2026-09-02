#!/usr/bin/env python3
"""test_windm_edsclrm_lhs.py — pin the windm/edsclrm tridiagonal LHS assembly.

`windm_edsclrm_lhs` (advance_windm_edsclrm_module.py ↔ ...F90:windm_edsclrm_lhs, F90:2236) assembles the banded
(3, ngrdcol, nzt) LHS for the implicit wind/edsclr eddy-diffusion solve. It is on the live path (every standalone
case advances u/v) but validated only end-to-end. This pins its four contributions against an INDEPENDENT
transcription of the Fortran:
  lhs = 0.5·lhs_diff            (Crank-Nicholson factor on all 3 bands; F90:2308-2316)
  lhs[diag] += 1/dt             (backward-Euler time accumulation on the centre band; F90:2313)
  lhs[:, :-1] += lhs_ma_zt[:, :-1]   (mean advection, INTERIOR levels only — excludes the top; F90:2334-2340)
  lhs[diag, k_lb] += invrs_rho_ds_zt·invrs_dzt·rho_ds_zm·(u_star²/wind_speed)  at the lower boundary
                                (implicit surface momentum flux, l_imp_sfc_momentum_flux=True; F90:2336-2343)
Plus a structural check that the surface term is localized to the diagonal band at the lower boundary (nowhere else),
and that the Crank-Nicholson factor scales the diffusion contribution linearly. Oracle-independent; never SKIPs. (iter 544)
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

from clubb_jax.src.CLUBB_core.advance_windm_edsclrm_module import windm_edsclrm_lhs
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _NZT = 2, 40.0, 8


def _ref(lhs_diff, lhs_ma_zt, dt, invrs_rho_ds_zt, rho_ds_zm, u_star_sqd, wind_speed, invrs_dzt):
    """Independent NumPy transcription of windm_edsclrm_lhs (the F90 assembly)."""
    lhs = 0.5 * lhs_diff.copy()
    lhs[1] += 1.0 / dt
    lhs[:, :, :-1] += lhs_ma_zt[:, :, :-1]
    sfc = invrs_rho_ds_zt[:, 0] * invrs_dzt[:, 0] * rho_ds_zm[:, 0] * (u_star_sqd / wind_speed[:, 0])
    lhs[1, :, 0] += sfc
    return lhs


def _inputs(rng, gr):
    lhs_diff = rng.uniform(-2.0, 2.0, (3, _NG, _NZT))
    lhs_ma_zt = rng.uniform(-1.0, 1.0, (3, _NG, _NZT))
    invrs_rho_ds_zt = rng.uniform(0.8, 1.2, (_NG, _NZT))
    rho_ds_zm = rng.uniform(0.8, 1.2, (_NG, _NZT + 1))
    u_star_sqd = rng.uniform(0.01, 0.5, (_NG,))
    wind_speed = rng.uniform(1.0, 15.0, (_NG, _NZT))
    invrs_dzt = np.asarray(gr.invrs_dzt)
    return lhs_diff, lhs_ma_zt, invrs_rho_ds_zt, rho_ds_zm, u_star_sqd, wind_speed, invrs_dzt


def test_lhs_matches_f90_assembly():
    gr = setup_grid(_NG, _DZ, _DZ, _DZ * (_NZT + 1))
    rng = np.random.default_rng(544)
    dt = 60.0
    lhs_diff, lhs_ma_zt, invrs_rho_ds_zt, rho_ds_zm, u_star_sqd, wind_speed, invrs_dzt = _inputs(rng, gr)
    got = np.asarray(windm_edsclrm_lhs(jnp.asarray(lhs_diff), jnp.asarray(lhs_ma_zt), dt,
                                       jnp.asarray(invrs_rho_ds_zt), jnp.asarray(rho_ds_zm),
                                       jnp.asarray(u_star_sqd), jnp.asarray(wind_speed), gr, 0, 0))
    ref = _ref(lhs_diff, lhs_ma_zt, dt, invrs_rho_ds_zt, rho_ds_zm, u_star_sqd, wind_speed, invrs_dzt)
    worst = float(np.max(np.abs(got - ref)))
    assert worst < 1e-13, f"windm_edsclrm_lhs mismatch vs F90 assembly {worst:.2e}"
    print(f"  windm_edsclrm_lhs: 0.5·diff + 1/dt + MA-interior + sfc-flux match F90 (worst {worst:.1e})  PASS")


def test_surface_term_localized_and_cn_scaling():
    gr = setup_grid(_NG, _DZ, _DZ, _DZ * (_NZT + 1))
    rng = np.random.default_rng(11)
    dt = 60.0
    lhs_diff, lhs_ma_zt, invrs_rho_ds_zt, rho_ds_zm, u_star_sqd, wind_speed, invrs_dzt = _inputs(rng, gr)
    args = (jnp.asarray(lhs_diff), jnp.asarray(lhs_ma_zt), dt, jnp.asarray(invrs_rho_ds_zt),
            jnp.asarray(rho_ds_zm), jnp.asarray(u_star_sqd), jnp.asarray(wind_speed), gr, 0, 0)
    base = np.asarray(windm_edsclrm_lhs(*args))
    # Zero u_star → surface term vanishes; the ONLY change must be band[1] at k=0.
    args_no_sfc = list(args); args_no_sfc[5] = jnp.zeros((_NG,))
    no_sfc = np.asarray(windm_edsclrm_lhs(*args_no_sfc))
    diff = base - no_sfc
    mask = np.zeros_like(diff, dtype=bool); mask[1, :, 0] = True
    assert np.max(np.abs(diff[~mask])) < 1e-18, "surface term leaked outside the diagonal/lower-boundary entry"
    assert np.all(np.abs(diff[1, :, 0]) > 0), "surface term missing at the lower boundary"
    # Crank-Nicholson: doubling lhs_diff doubles the diffusion contribution (LHS minus the diff-independent parts).
    args2 = list(args); args2[0] = jnp.asarray(2.0 * lhs_diff)
    dbl = np.asarray(windm_edsclrm_lhs(*args2))
    assert np.allclose(dbl - base, 0.5 * lhs_diff, atol=1e-13), "Crank-Nicholson 0.5 factor not linear in lhs_diff"
    print("  surface term localized to diag@lower-boundary; CN 0.5·diff scales linearly  PASS")


def main():
    print("test_windm_edsclrm_lhs:")
    test_lhs_matches_f90_assembly()
    test_surface_term_localized_and_cn_scaling()
    print("All windm_edsclrm_lhs checks PASSED")


if __name__ == "__main__":
    main()

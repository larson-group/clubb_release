"""Hydrostatic pressure / altitude integration — JAX port of hydrostatic_module.F90.

Mirrors clubb_release/src/Input_fields/hydrostatic_module.F90:
  hydrostatic                  — p/exner/rho on both grids (wraps calc_pressure.init_pressure)
  inverse_hydrostatic          — altitudes of a pressure-coordinate sounding's levels
  calc_ref_z_linear_thvm       — reference altitudes from thvm+exner (log-mean thvm integration)
  calc_ref_z_sfc_linear_thvm   — surface altitude relative to a sounding level

The forward `init_pressure` lives in CLUBB_core/calc_pressure.py (calc_pressure.F90); `hydrostatic` wraps it
(mirroring the Fortran `use calc_pressure`/init_pressure call). All pure-jnp → differentiable.
"""

import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.constants_clubb import Cp, grav, kappa, p0, Rd
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.calc_pressure import init_pressure

_CP_OV_G = Cp / grav
_THVM_TOL = 1e-12   # relative threshold for the thvm(k)=thvm(k-1) degenerate branch


def hydrostatic(thvm, p_sfc, gr):
    """Compute hydrostatic pressure, Exner, and density.

    JAX port of hydrostatic_module.F90:hydrostatic.

    Args:
        thvm:   (ngrdcol, nzt) virtual potential temperature [K]
        p_sfc:  (ngrdcol,) surface pressure [Pa]
        gr:     grid object

    Returns:
        p_in_Pa:    (ngrdcol, nzt) pressure [Pa]
        p_in_Pa_zm: (ngrdcol, nzm) pressure [Pa]
        exner:      (ngrdcol, nzt) Exner function [-]
        exner_zm:   (ngrdcol, nzm) Exner function [-]
        rho:        (ngrdcol, nzt) air density [kg/m^3]
        rho_zm:     (ngrdcol, nzm) air density [kg/m^3]
    """
    p_in_Pa, exner_zt, p_in_Pa_zm, exner_zm = init_pressure(thvm, p_sfc, gr)

    # Interpolate thvm to momentum levels (no floor, as in Fortran hydrostatic)
    thvm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, thvm)

    rho    = p_in_Pa    / (Rd * thvm    * exner_zt)
    rho_zm = p_in_Pa_zm / (Rd * thvm_zm * exner_zm)

    return p_in_Pa, p_in_Pa_zm, exner_zt, exner_zm, rho, rho_zm


def calc_ref_z_linear_thvm(thvm, exner):
    """Reference altitudes (relative to the lowest level) of a sounding from thvm + exner.

    JAX port of hydrostatic_module.F90:calc_ref_z_linear_thvm — integrates dz = -(Cp/g) thvm d(exner) assuming
    thvm varies LINEARLY IN Z between levels, which yields the LOG-MEAN of thvm:
      ref_z[k] = ref_z[k-1] - (Cp/g)(exner[k]-exner[k-1]) (thvm[k]-thvm[k-1]) / log(thvm[k]/thvm[k-1]),
    reducing to -(Cp/g)(exner[k]-exner[k-1]) thvm[k] when thvm[k]==thvm[k-1]. Differentiable (a cumulative sum).
    """
    thvm = jnp.asarray(thvm); exner = jnp.asarray(exner)
    d_exner = exner[1:] - exner[:-1]
    dthvm = thvm[1:] - thvm[:-1]
    big = jnp.abs(dthvm) > _THVM_TOL * thvm[1:]
    safe_ratio = jnp.where(big, thvm[1:] / thvm[:-1], 2.0)          # 2.0 → log(2)!=0 keeps grad finite
    logmean = jnp.where(big, dthvm / jnp.log(safe_ratio), thvm[1:])
    dz = -_CP_OV_G * d_exner * logmean
    return jnp.concatenate([jnp.zeros(1), jnp.cumsum(dz)])


def calc_ref_z_sfc_linear_thvm(thvm_km1, thvm_k, z_km1, z_k, exner_km1, exner_sfc):
    """Altitude of the surface relative to sounding level k-1 (hydrostatic_module.F90:calc_ref_z_sfc_linear_thvm)."""
    dthvm = thvm_k - thvm_km1
    big = jnp.abs(dthvm) > _THVM_TOL * jnp.abs(thvm_k)
    slope = jnp.where(big, dthvm, 1.0) / (z_k - z_km1)             # guard /dthvm via the where below
    val_big = ((z_k - z_km1) / jnp.where(big, dthvm, 1.0)) * (
        thvm_km1 * jnp.exp(-_CP_OV_G * (exner_sfc - exner_km1) * slope)
        - thvm_km1 + slope * z_km1)
    val_eq = z_km1 - _CP_OV_G * (exner_sfc - exner_km1) * thvm_k
    return jnp.where(big, val_big, val_eq)


def inverse_hydrostatic(p_sfc, zm_init, thvm, exner):
    """Heights of a pressure-coordinate sounding's levels (hydrostatic_module.F90:inverse_hydrostatic).

    Given the surface pressure/altitude and the sounding's thvm + exner profiles (exner decreasing with level),
    integrate the inverse hydrostatic equation for the altitude of each sounding level. The exact INVERSE of the
    forward `init_pressure` (same log-mean scheme), so z -> exner -> z round-trips exactly.

    Args:
        p_sfc:   surface pressure [Pa] (scalar).
        zm_init: surface altitude [m] (scalar).
        thvm:    virtual potential temperature at the sounding levels (nlevels,) [K].
        exner:   Exner function at the sounding levels (nlevels,), decreasing with level.
    Returns:
        z (nlevels,) [m].
    """
    thvm = jnp.asarray(thvm); exner = jnp.asarray(exner)
    ref_z_snd = calc_ref_z_linear_thvm(thvm, exner)
    exner_sfc = (p_sfc / p0) ** kappa
    e = np.asarray(exner)
    if float(exner_sfc) < e[-1]:
        raise ValueError("inverse_hydrostatic: the entire sounding is below the model surface.")
    # Locate the sounding interval bracketing the surface (exner decreasing → last level with exner >= exner_sfc).
    if float(exner_sfc) > e[0]:
        low, high = 0, 1                       # surface below the lowest sounding level (linear extension)
    else:
        low = int(np.nonzero(e >= float(exner_sfc))[0][-1])
        high = min(low + 1, len(e) - 1)
    ref_z_sfc = calc_ref_z_sfc_linear_thvm(thvm[low], thvm[high], ref_z_snd[low], ref_z_snd[high],
                                           exner[low], exner_sfc)
    z_snd_bottom = zm_init - ref_z_sfc
    return z_snd_bottom + ref_z_snd


__all__ = [
    "hydrostatic",
    "inverse_hydrostatic",
    "calc_ref_z_linear_thvm",
    "calc_ref_z_sfc_linear_thvm",
]

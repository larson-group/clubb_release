"""Contains subroutines for the COBRA CO2 case."""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# -----------------------------------------------------------------------
def cobra_sfclyr(
    ngrdcol,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    time,
    z,
    rho_sfc,
    thlm_sfc,
    ubar,
):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to the format used for the GCSS ARM case.

    Notes:
        The data has been altered so it can be used for the COBRA CO2 case.
    """
    M_da = 0.02897
    ntimes = 49

    # COBRA roughness height
    # z0 = 0.035  # ARM momentum roughness height
    z0 = 1.75

    # --------------------- Begin Code ---------------------
    # Compute heat and moisture fluxes from ARM data in (W/m2)

    # Use time_select to caluclate the indexes before and after time
    # and the time fraction necessary for linear_interp_factor
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        ntimes,
        time_dependent_input.time_sfc_given,
    )

    # Interpolate fluxes
    heat_flx = linear_interp_factor(
        time_frac,
        time_dependent_input.sens_ht_given[after_time],
        time_dependent_input.sens_ht_given[before_time],
    )
    moisture_flx = linear_interp_factor(
        time_frac,
        time_dependent_input.latent_ht_given[after_time],
        time_dependent_input.latent_ht_given[before_time],
    )
    CO2_flx = linear_interp_factor(
        time_frac,
        time_dependent_input.CO2_sfc_given[after_time],
        time_dependent_input.CO2_sfc_given[before_time],
    )
    T_sfc_calc = linear_interp_factor(
        time_frac,
        time_dependent_input.T_sfc_given[after_time],
        time_dependent_input.T_sfc_given[before_time],
    )

    T_sfc = jnp.full(ngrdcol, T_sfc_calc)

    # Convert heat_flx and moisture_flx to natural units
    heat_flx2 = convert_sens_ht_to_km_s(heat_flx, rho_sfc)
    moisture_flx2 = convert_latent_ht_to_m_s(moisture_flx, rho_sfc)

    # Convert CO2 surface flux to natural units.
    # The CO2 flux has been given in units of: umol/(m^2 s).
    # umol stands for micromoles. The CO2 concentration in this code is in
    # units of ppmv, which is also the molar mixing ratio times 10^6.
    # The units are: 10^6 * [ mol (CO2) / mol (dry air) ].
    # w'CO2' = (Flux) * [ M (dry air) / rho (dry air) ];
    # where M is the molecular weight of dry air.
    CO2_flx2 = CO2_flx * (M_da / rho_sfc)

    # Heat flux in units of (m2/s3) (needed by diag_ustar)
    bflx = grav / thlm_sfc * heat_flx2

    # Compute ustar
    ustar = diag_ustar(z, bflx, ubar, z0)

    # Assign fluxes
    wpthlp_sfc = heat_flx2
    wprtp_sfc = moisture_flx2

    wpsclrp_sfc = jnp.zeros((ngrdcol, sclr_dim), dtype=rho_sfc.dtype)
    wpedsclrp_sfc = jnp.zeros((ngrdcol, edsclr_dim), dtype=rho_sfc.dtype)
    if sclr_dim > 0:
        if sclr_idx.iisclr_CO2 > 0:
            wpsclrp_sfc = wpsclrp_sfc.at[:, sclr_idx.iisclr_CO2 - 1].set(CO2_flx2)
        if sclr_idx.iisclr_thl > 0:
            wpsclrp_sfc = wpsclrp_sfc.at[:, sclr_idx.iisclr_thl - 1].set(wpthlp_sfc)
        if sclr_idx.iisclr_rt > 0:
            wpsclrp_sfc = wpsclrp_sfc.at[:, sclr_idx.iisclr_rt - 1].set(wprtp_sfc)
    if edsclr_dim > 0:
        if sclr_idx.iiedsclr_CO2 > 0:
            wpedsclrp_sfc = wpedsclrp_sfc.at[:, sclr_idx.iiedsclr_CO2 - 1].set(CO2_flx2)
        if sclr_idx.iiedsclr_thl > 0:
            wpedsclrp_sfc = wpedsclrp_sfc.at[:, sclr_idx.iiedsclr_thl - 1].set(wpthlp_sfc)
        if sclr_idx.iiedsclr_rt > 0:
            wpedsclrp_sfc = wpedsclrp_sfc.at[:, sclr_idx.iiedsclr_rt - 1].set(wprtp_sfc)
    return wpthlp_sfc, wprtp_sfc, ustar, wpsclrp_sfc, wpedsclrp_sfc, T_sfc


__all__ = ["cobra_sfclyr"]

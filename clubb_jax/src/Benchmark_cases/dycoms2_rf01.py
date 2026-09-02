"""Contains subroutines for the DYCOMS II RF01 case.

References:
    http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    compute_wprtp_sfc,
    compute_wpthlp_sfc,
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# ----------------------------------------------------------------------
def dycoms2_rf01_tndcy(ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx):
    """Subroutine to set theta and water tendencies for DYCOMS RF01 case.

    References:
        http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html
    """
    thlm_forcing = jnp.zeros((ngrdcol, gr.nzt))
    rtm_forcing = jnp.zeros((ngrdcol, gr.nzt))
    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim))
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim))

    if sclr_dim > 0:
        # Test scalars with thetal and rt if desired
        if sclr_idx.iisclr_thl > 0:
            sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_thl - 1].set(thlm_forcing)
        if sclr_idx.iisclr_rt > 0:
            sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_rt - 1].set(rtm_forcing)
    if edsclr_dim > 0:
        if sclr_idx.iiedsclr_thl > 0:
            edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_thl - 1].set(
                thlm_forcing
            )
        if sclr_idx.iiedsclr_rt > 0:
            edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_rt - 1].set(
                rtm_forcing
            )
    return thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing


# ======================================================================
def dycoms2_rf01_sfclyr(
    ngrdcol,
    time,
    sfctype,
    p_sfc,
    exner_sfc,
    ubar,
    thlm_sfc,
    rtm_sfc,
    rho_sfc,
    saturation_formula,
):
    """This subroutine computes surface fluxes of heat and moisture according
    to GCSS DYCOMS II RF 01 specifications.

    References:
        http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html
    """
    # -----------------BEGIN CODE-----------------------
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        time_dependent_input.time_sfc_given.size,
        time_dependent_input.time_sfc_given,
    )
    sens_ht = linear_interp_factor(
        time_frac,
        time_dependent_input.sens_ht_given[after_time],
        time_dependent_input.sens_ht_given[before_time],
    )
    latent_ht = linear_interp_factor(
        time_frac,
        time_dependent_input.latent_ht_given[after_time],
        time_dependent_input.latent_ht_given[before_time],
    )
    T_sfc_interp = linear_interp_factor(
        time_frac,
        time_dependent_input.T_sfc_given[after_time],
        time_dependent_input.T_sfc_given[before_time],
    )

    # Compute heat and moisture fluxes
    ustar = jnp.full(ngrdcol, 0.25)
    T_sfc = jnp.full(ngrdcol, T_sfc_interp)
    if sfctype == 0:
        wpthlp_sfc = convert_sens_ht_to_km_s(sens_ht, rho_sfc)
        wprtp_sfc = convert_latent_ht_to_m_s(latent_ht, rho_sfc)
    elif sfctype == 1:
        Cd = jnp.full(ngrdcol, 0.0011)
        rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)
        wpthlp_sfc = compute_wpthlp_sfc(
            ngrdcol,
            Cd,
            ubar,
            thlm_sfc,
            T_sfc,
            exner_sfc,
        )
        wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cd, ubar, rtm_sfc, rsat)
    else:
        raise ValueError(f"Invalid sfctype value = {sfctype}")
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["dycoms2_rf01_tndcy", "dycoms2_rf01_sfclyr"]

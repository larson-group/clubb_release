"""Contains subroutines for the RICO case."""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    compute_momentum_flux,
    compute_ubar,
    compute_wprtp_sfc,
    compute_wpthlp_sfc,
)
from clubb_jax.src.Benchmark_cases.spec_hum_to_mixing_ratio import (
    force_spec_hum_to_mixing_ratio,
)
from clubb_jax.src.CLUBB_core.constants_clubb import g_per_kg
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# ----------------------------------------------------------------------
def rico_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, rtm, exner):
    """Subroutine to apply case-specific forcings to RICO case
    (Michael Falk, 13 Dec 2006).

    References:
        ATEX: http://www.atmos.ucla.edu/~bstevens/gcss/setup.html
        RICO: http://www.knmi.nl/samenw/rico/setup3d.html
    """
    # --------------------- Begin Code ---------------------
    # Compute large-scale horizontal temperature advection
    # NEW-- "And Radiation"... 15 Dec 2006, Michael Falk
    # Equations located in 1D models > Set up short composite run on reference site
    t_tendency = jnp.where(
        gr.zt < 4000.0,
        -2.51 / 86400.0 + (-2.18 + 2.51) / (86400.0 * 4000.0) * gr.zt,
        jnp.where(
            gr.zt < 5000.0,
            -2.18 / 86400.0
            + 2.18 / (86400.0 * (5000.0 - 4000.0)) * (gr.zt - 4000.0),
            0.0,
        ),
    )

    # Convert to units of [K s^-1] but potential T instead of T
    thlm_forcing = t_tendency / exner

    # Compute large-scale horizontal moisture advection [g kg^-1 s^-1]
    # Equations located in 1D models > Set up short composite run on reference site
    qtm_forcing = jnp.where(
        gr.zt < 3000.0,
        -1.0 / 86400.0 + (0.345 + 1.0) / (86400.0 * 3000.0) * gr.zt,
        jnp.where(
            gr.zt < 4000.0,
            0.345 / 86400.0,
            jnp.where(
                gr.zt < 5000.0,
                0.345 / 86400.0
                + (-0.345) / (86400.0 * (5000.0 - 4000.0)) * (gr.zt - 4000.0),
                0.0,
            ),
        ),
    )
    qtm_forcing = qtm_forcing / g_per_kg

    # Convert forcings from terms of total water specific humidity to terms of
    # total water mixing ratio.
    rtm_forcing = force_spec_hum_to_mixing_ratio(
        ngrdcol,
        gr.nzt,
        rtm,
        qtm_forcing,
    )

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=rtm.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=rtm.dtype)

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


# ----------------------------------------------------------------------
def rico_sfclyr(
    ngrdcol,
    time,
    um_sfc,
    vm_sfc,
    thlm,
    rtm,
    z_bot,
    p_sfc,
    exner_sfc,
    saturation_formula,
):
    """Surface forcing subroutine for RICO case. Written
    December 2006 by Michael Falk.

    Updated to use specific formulations for surface fluxes as specified in
    the RICO 3D LES specification, in hopes that they'll be more accurate.

    References:
        ATEX: http://www.atmos.ucla.edu/~bstevens/gcss/setup.html
        RICO: http://www.knmi.nl/samenw/rico/setup3d.html
    """
    # Constants
    C_10 = 0.0013
    C_m_20 = 0.001229
    C_h_20 = 0.001094
    C_q_20 = 0.001133
    z0 = 0.00015
    standard_flux_alt = 20.0

    # --------------------BEGIN CODE----------------------------
    # interpolate variables from time_dependent_input
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        time_dependent_input.time_sfc_given.size,
        time_dependent_input.time_sfc_given,
    )
    T_sfc_interp = linear_interp_factor(
        time_frac,
        time_dependent_input.T_sfc_given[after_time],
        time_dependent_input.T_sfc_given[before_time],
    )

    # Choose which scheme to use
    l_use_old_atex = False

    # Define variable values
    ubar = compute_ubar(ngrdcol, um_sfc, vm_sfc)
    T_sfc = jnp.full(ngrdcol, T_sfc_interp)
    ustar = jnp.full(ngrdcol, 0.3)
    rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)

    # Compute heat and moisture fluxes
    if l_use_old_atex:
        # (Stevens, et al. 2000, eq 3)
        # Modification in case lowest model level isn't at 10 m, from ATEX specification
        Cz = C_10 * (jnp.log(10.0 / z0) / jnp.log(z_bot / z0)) * (
            jnp.log(10.0 / z0) / jnp.log(z_bot / z0)
        )
        wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, Cz, ubar, thlm, T_sfc, exner_sfc)
        wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cz, ubar, rtm, rsat)
        upwp_sfc, vpwp_sfc = compute_momentum_flux(
            ngrdcol,
            um_sfc,
            vm_sfc,
            ubar,
            ustar,
        )
    else:
        # Modification in case lowest model level isn't at 10 m, from ATEX specification
        Cm = C_m_20 * (
            jnp.log(standard_flux_alt / z0) / jnp.log(z_bot / z0)
        ) * (jnp.log(standard_flux_alt / z0) / jnp.log(z_bot / z0))

        # Modification in case lowest model level isn't at 10 m, from ATEX specification
        Ch = C_h_20 * (
            jnp.log(standard_flux_alt / z0) / jnp.log(z_bot / z0)
        ) * (jnp.log(standard_flux_alt / z0) / jnp.log(z_bot / z0))

        # Modification in case lowest model level isn't at 10 m, from ATEX specification
        Cq = C_q_20 * (
            jnp.log(standard_flux_alt / z0) / jnp.log(z_bot / z0)
        ) * (jnp.log(standard_flux_alt / z0) / jnp.log(z_bot / z0))
        wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, Ch, ubar, thlm, T_sfc, exner_sfc)
        wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cq, ubar, rtm, rsat)
        upwp_sfc = -um_sfc * Cm * ubar
        vpwp_sfc = -vm_sfc * Cm * ubar

    return upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["rico_tndcy", "rico_sfclyr"]

"""Contains subroutines for the M-PACE B intercomparison.

References:
    http://science.arm.gov/wg/cpm/scm/scmic5/index.html
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Rd, g_per_kg, grav, sec_per_day
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# ----------------------------------------------------------------------
def mpace_b_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, p_in_Pa, thvm):
    """Subroutine to compute large-scale subsidence for mpace_b case (Michael
    Falk, 21 July 2006). Added ice and radiation based on Adam Smith Nov 11
    case, 27 July 2006. Comments and documentation added 31 July 2006.

    References:
        Liou, Wallace and Hobbs, Shettle and Weinman
        http://science.arm.gov/wg/cpm/scm/scmic5/index.html
    """
    # Local constants, subsidence
    D = 5.8e-6
    p_sfc = 101000.0
    pinv = 85000.0

    # --------------------- Begin Code ---------------------
    # Compute vertical motion
    velocity_omega = jnp.minimum(D * (p_sfc - p_in_Pa), D * (p_sfc - pinv))
    wm_zt = -velocity_omega * Rd * thvm / p_in_Pa / grav

    # Interpolate
    wm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, wm_zt)

    # Boundary conditions
    wm_zm = wm_zm.at[:, 0].set(0.0)
    wm_zm = wm_zm.at[:, gr.nzm - 1].set(0.0)

    # Compute large-scale tendencies
    t_tendency = jnp.minimum(-4.0, -15.0 * (1.0 - (p_sfc - p_in_Pa) / 21818.0))
    thlm_forcing = t_tendency * (p_sfc / p_in_Pa) ** (Rd / Cp) / sec_per_day
    rtm_forcing = (
        jnp.minimum(0.164, -3.0 * (1.0 - (p_sfc - p_in_Pa) / 15171.0))
        / g_per_kg
        / sec_per_day
    )

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=p_in_Pa.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=p_in_Pa.dtype)

    # Test scalars with thetal and rt if desired
    if sclr_dim > 0:
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
    return wm_zt, wm_zm, thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing


# ----------------------------------------------------------------------
def mpace_b_sfclyr(ngrdcol, time, rho_sfc):
    """Surface forcing subroutine for mpace_b case. Written July-
    November 2006 by Michael Falk.

    References:
        http://science.arm.gov/wg/cpm/scm/scmic5/index.html
    """
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        time_dependent_input.time_sfc_given.size,
        time_dependent_input.time_sfc_given,
    )

    # Get sens_ht and latent_ht from the input.
    sensible_heat_flx = linear_interp_factor(
        time_frac,
        time_dependent_input.sens_ht_given[after_time],
        time_dependent_input.sens_ht_given[before_time],
    )
    latent_heat_flx = linear_interp_factor(
        time_frac,
        time_dependent_input.latent_ht_given[after_time],
        time_dependent_input.latent_ht_given[before_time],
    )

    # Declare the value of ustar.
    ustar = jnp.full(ngrdcol, 0.25)

    # Compute heat and moisture fluxes
    wpthlp_sfc = convert_sens_ht_to_km_s(sensible_heat_flx, rho_sfc)
    wprtp_sfc = convert_latent_ht_to_m_s(latent_heat_flx, rho_sfc)
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["mpace_b_tndcy", "mpace_b_sfclyr"]

"""Contains subroutines for the GABLS2 case.

References:
    http://www.misu.su.se/~gunilla/gabls/
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Rd, grav, p0, sec_per_hr
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# -------------------------------------------------------------------------------
def gabls2_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time, time_initial):
    """Subroutine to apply case-specific forcings to GABLS2 case
    (Michael Falk, 29 Dec 2006).

    References:
        http://people.su.se/~gsven/gabls/
    """
    # --------------------- Begin Code ---------------------
    # Compute vertical motion
    # 93600 seconds = 26 hours of simulation time;
    # per GABLS2 specification
    wm_zt = jnp.where(
        time > time_initial + 93600.0,
        jnp.where(gr.zt <= 1000.0, -0.005 * gr.zt / 1000.0, -0.005),
        jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype),
    )

    wm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, wm_zt)

    # Boundary conditions on vertical motion.
    wm_zm = wm_zm.at[:, 0].set(0.0)
    wm_zm = wm_zm.at[:, gr.nzm - 1].set(0.0)

    # Compute large-scale horizontal temperature advection
    thlm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)

    # Compute large-scale horizontal moisture advection [g/kg/s]
    rtm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=gr.zt.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=gr.zt.dtype)

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


# -------------------------------------------------------------------------------
def gabls2_sfclyr(
    ngrdcol,
    time,
    time_initial,
    lowest_level,
    p_sfc,
    ubar,
    thlm,
    rtm,
    exner_sfc,
    saturation_formula,
):
    """Surface forcing subroutine for GABLS2 case. Written
    29 December 2006 by Michael Falk.

    References:
        http://people.su.se/~gsven/gabls/
    """
    # Local constants
    standard_flux_alt = 10.0
    C_10 = 0.0013
    z0 = 0.03

    # ---- Begin Code ----
    # lowest model level isn't at 10 m,
    # from ATEX specification (Stevens, et al. 2000, eq 3)
    time_in_hours_init = 14.0
    time_in_hours = (time - time_initial) / sec_per_hr + time_in_hours_init

    # at initial time,
    # time_in_hours = 14
    # (14 local; 19 UTC)

    # Compute sea surface temperature
    T_sfc_calc = jnp.where(
        time_in_hours <= 17.4,
        -10.0 - 25.0 * jnp.cos(time_in_hours * 0.22 + 0.2),
        jnp.where(
            time_in_hours <= 30.0,
            -0.54 * time_in_hours + 15.2,
            jnp.where(
                time_in_hours <= 41.9,
                -7.0 - 25.0 * jnp.cos(time_in_hours * 0.21 + 1.8),
                jnp.where(
                    time_in_hours <= 53.3,
                    -0.37 * time_in_hours + 18.0,
                    jnp.where(
                        time_in_hours <= 65.6,
                        -4.0 - 25.0 * jnp.cos(time_in_hours * 0.22 + 2.5),
                        4.4,
                    ),
                ),
            ),
        ),
    )

    Cz = C_10 * (
        jnp.log(standard_flux_alt / z0) / jnp.log(lowest_level / z0)
    ) * (jnp.log(standard_flux_alt / z0) / jnp.log(lowest_level / z0))
    T_sfc = jnp.full(ngrdcol, T_sfc_calc + 273.15)
    rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)

    # Compute heat and moisture fluxes
    wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, Cz, ubar, thlm, T_sfc, exner_sfc)
    wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cz, ubar, rtm, rsat)

    # The latent heat flux at the surface is 2.5% of its potential value
    wprtp_sfc = wprtp_sfc * 0.025

    # Compute momentum fluxes
    sstheta = T_sfc * (p0 / p_sfc) ** (Rd / Cp)
    bflx = wpthlp_sfc * grav / sstheta
    ustar = diag_ustar(lowest_level, bflx, ubar, z0)
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["gabls2_tndcy", "gabls2_sfclyr"]

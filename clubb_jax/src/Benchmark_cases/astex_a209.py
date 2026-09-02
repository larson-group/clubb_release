"""Contains subroutines for the ASTEX KK case.

References:
    http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# ----------------------------------------------------------------------
def astex_a209_tndcy(sclr_dim, edsclr_dim, sclr_idx, gr):
    """Subroutine to set theta and water tendencies for ASTEX KK case.

    References:
        http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
    """
    # --------------------- Begin Code ---------------------
    # Compute large-scale subsidence
    wm_zt = -5.0e-6 * gr.zt[0, :]

    # Interpolate to momentum levels
    wm_zm = zt2zm(gr.nzm, gr.nzt, 1, gr, wm_zt[jnp.newaxis, :])[0, :]

    # Boundary conditions on zm
    wm_zm = wm_zm.at[0].set(0.0)
    wm_zm = wm_zm.at[gr.nzm - 1].set(0.0)

    # Radiative theta-l tendency
    thlm_forcing = jnp.zeros(gr.nzt)

    # Large scale advective moisture tendency
    rtm_forcing = jnp.zeros(gr.nzt)

    sclrm_forcing = jnp.zeros((gr.nzt, sclr_dim))
    edsclrm_forcing = jnp.zeros((gr.nzt, edsclr_dim))

    # Test scalars with thetal and rt if desired
    if sclr_idx.iisclr_thl > 0:
        sclrm_forcing = sclrm_forcing.at[:, sclr_idx.iisclr_thl - 1].set(thlm_forcing)
    if sclr_idx.iisclr_rt > 0:
        sclrm_forcing = sclrm_forcing.at[:, sclr_idx.iisclr_rt - 1].set(rtm_forcing)
    if sclr_idx.iiedsclr_thl > 0:
        edsclrm_forcing = edsclrm_forcing.at[:, sclr_idx.iiedsclr_thl - 1].set(thlm_forcing)
    if sclr_idx.iiedsclr_rt > 0:
        edsclrm_forcing = edsclrm_forcing.at[:, sclr_idx.iiedsclr_rt - 1].set(rtm_forcing)
    return wm_zt, wm_zm, thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing


# ----------------------------------------------------------------------
def astex_a209_sfclyr(
    ngrdcol,
    time,
    ubar,
    rtm,
    thlm,
    lowestlevel,
    exner_sfc,
    p_sfc,
    saturation_formula,
):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to ASTEX with Khairoutdinov and Kogan
    alteration.

    References:
        http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
    """
    # Parameter Constants
    ntimes = 41
    C_h_20 = 0.001094
    C_q_20 = 0.001133
    z0 = 0.00015
    standard_flux_alt = 20.0

    # -----------------BEGIN CODE-------------------------
    # sensible_heat_flx = 10.0
    # latent_heat_flx = 25.0

    # Use time_select to determine the time indexes before and after time
    # and to calculate the time fraction necessary for linear_interp_factor
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        ntimes,
        time_dependent_input.time_sfc_given,
    )

    # Interpolate the value for T_sfc based on time.
    T_sfc_interp = linear_interp_factor(
        time_frac,
        time_dependent_input.T_sfc_given[after_time],
        time_dependent_input.T_sfc_given[before_time],
    )
    T_sfc = jnp.full(ngrdcol, T_sfc_interp)

    # We set ustar as it is set in rico
    ustar = jnp.full(ngrdcol, 0.155)
    rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)

    # (Stevens, et al. 2000, eq 3)
    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    Ch = C_h_20 * (
        (jnp.log(standard_flux_alt / z0)) / (jnp.log(lowestlevel / z0))
    ) * ((jnp.log(standard_flux_alt / z0)) / (jnp.log(lowestlevel / z0)))

    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    Cq = C_q_20 * (
        (jnp.log(standard_flux_alt / z0)) / (jnp.log(lowestlevel / z0))
    ) * ((jnp.log(standard_flux_alt / z0)) / (jnp.log(lowestlevel / z0)))

    # Compute heat and moisture fluxes
    wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, Ch, ubar, thlm, T_sfc, exner_sfc)
    wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cq, ubar, rtm, rsat)

    # wpthlp_sfc = sensible_heat_flx / ( rho_sfc * Cp )
    # wprtp_sfc = latent_heat_flx / ( rho_sfc * Lv )

    # Momentum fluxes are computed elsewhere
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["astex_a209_tndcy", "astex_a209_sfclyr"]

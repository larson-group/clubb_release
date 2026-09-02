"""Contains subroutines for the long GCSS ATEX case."""

from __future__ import annotations

import jax.numpy as jnp
from jax import lax

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    compute_wprtp_sfc,
    compute_wpthlp_sfc,
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


def calc_forcings(ngrdcol, gr):
    """Calculate the long-ATEX liquid-water and total-water tendencies."""
    # --------------------- Begin Code ---------------------
    # Theta-l tendency ! Hing - known magic numbers
    thlm_forcing = jnp.where(
        (gr.zt >= 0.0) & (gr.zt < 1400.0),
        -3.5805e-5,
        jnp.where(
            (gr.zt >= 1400.0) & (gr.zt < 1650.0),
            -3.5805e-5 + 1.1935e-5 * (gr.zt - 1400.0) * 0.004,
            jnp.where(
                (gr.zt >= 1650.0) & (gr.zt < 2990.0),
                -2.3870e-5 - 0.1155e-5 * (gr.zt - 1650.0) / 1350.0,
                0.0,
            ),
        ),
    )

    # Moisture tendency ! Hing - known magic numbers
    rtm_forcing = jnp.where(
        (gr.zt >= 0.0) & (gr.zt < 1050.0),
        -1.58e-8 * (1.0 - gr.zt / 1050.0),
        0.0,
    )
    return thlm_forcing, rtm_forcing


# ======================================================================
def atex_long_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time):
    """Subroutine to set theta-l and water tendencies and subsidence for the
    long ATEX case.

    References:
        Ong, H. The nontraditional Coriolis terms and trade-wind cumuli
    """
    # --------------------- Begin Code ---------------------
    # Large scale subsidence ! Hing - known magic numbers
    wm_zt = jnp.where(
        (gr.zt >= 0.0) & (gr.zt < 1050.0),
        -0.00636 * gr.zt / 1050.0,
        jnp.where(
            (gr.zt >= 1050.0) & (gr.zt < 1650.0),
            -0.00636 - 0.00079 * (gr.zt - 1050.0) / 600.0,
            -0.00715,
        ),
    )

    thlm_forcing, rtm_forcing = calc_forcings(ngrdcol, gr)

    # Spin up period ! Hing - known magic number
    spinup = 43200.0
    thlm_forcing, rtm_forcing, wm_zt = lax.cond(
        time < spinup,
        lambda arrays: (
            arrays[0] * time / spinup,
            arrays[1] * time / spinup,
            arrays[2] * time / spinup,
        ),
        lambda arrays: arrays,
        (thlm_forcing, rtm_forcing, wm_zt),
    )

    wm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, wm_zt)

    # Boundary conditions.
    wm_zm = wm_zm.at[:, 0].set(0.0)
    wm_zm = wm_zm.at[:, gr.nzm - 1].set(0.0)

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=wm_zt.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=wm_zt.dtype)

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


# ======================================================================
def atex_long_sfclyr(ngrdcol, time, ubar, thlm_sfc, rtm_sfc, exner_sfc, rho_sfc):
    """This subroutine computes surface fluxes of heat and moisture according
    to GCSS ATEX specifications.

    References:
        B. Stevens et al., 2000: Simulations of trade-wind cumuli
        under a strong inversion, J. Atmos. Sci, 58,1870-1891.
        http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
        ATEX_JAS2001.pdf
    """
    # -----------------BEGIN CODE-----------------------
    l_compute_flux = True

    # Interpolate T_sfc from time_dependent_input
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

    # Prescribe fluxes ! Hing - known magic number
    if not l_compute_flux:
        sens_ht = 10.70
        latent_ht = 154.0
        # Spin up period ! Hing - known magic number
        spinup = 86400.0
        if time < spinup:
            sens_ht = sens_ht - 2.20 * (spinup - time) / spinup
            latent_ht = latent_ht + 40.0 * (spinup - time) / spinup

    C_10 = jnp.full(ngrdcol, 0.0013)
    adjustment = jnp.full(ngrdcol, 0.0194664)
    ustar = jnp.full(ngrdcol, 0.3)
    T_sfc = jnp.full(ngrdcol, T_sfc_interp)

    if not l_compute_flux:
        wpthlp_sfc = convert_sens_ht_to_km_s(sens_ht, rho_sfc)
        wprtp_sfc = convert_latent_ht_to_m_s(latent_ht, rho_sfc)

    # Compute wpthlp_sfc and wprtp_sfc
    if l_compute_flux:
        wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, C_10, ubar, thlm_sfc, T_sfc, exner_sfc)
        wprtp_sfc = compute_wprtp_sfc(ngrdcol, C_10, ubar, rtm_sfc, adjustment)
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["atex_long_tndcy", "atex_long_sfclyr"]

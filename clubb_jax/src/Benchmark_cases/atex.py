"""Contains subroutines for the GCSS ATEX case.

TODO(port-scope): ``atex_tndcy`` omits the Fortran dycore-remapping arguments.
The standalone JAX driver owns only the CLUBB physics grid, so carrying those
unusable inputs would obscure the supported interface. Restore them only when
the driver also owns and tests the host-model dycore data flow.
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


def calc_forcings(ngrdcol, gr, z_inversion):
    """Calculate ATEX liquid-water and total-water tendencies."""
    # --------------------- Begin Code ---------------------
    z_inversion = z_inversion[:, None]

    # Theta-l tendency
    thlm_forcing = jnp.where(
        (gr.zt > 0.0) & (gr.zt < z_inversion),
        -1.1575e-5 * (3.0 - gr.zt / z_inversion),
        jnp.where(
            (gr.zt > z_inversion) & (gr.zt <= z_inversion + 300.0),
            -2.315e-5 * (1.0 - (gr.zt - z_inversion) / 300.0),
            0.0,
        ),
    )

    # Moisture tendency
    thlm_dtype = thlm_forcing.dtype
    rtm_forcing = jnp.where(
        (gr.zt > 0.0) & (gr.zt < z_inversion),
        -1.58e-8 * (1.0 - gr.zt / z_inversion),
        jnp.asarray(0.0, dtype=thlm_dtype),
    )
    return thlm_forcing, rtm_forcing


# ======================================================================
def atex_tndcy(
    ngrdcol,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    gr,
    time,
    time_initial,
    rtm,
    err_info,
):
    """Subroutine to set theta-l and water tendencies for ATEX case.

    References:
        B. Stevens et al., 2000: Simulations of trade-wind cumuli
        under a strong inversion, J. Atmos. Sci, 58,1870-1891.
        http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
        ATEX_JAS2001.pdf
    """
    # --------------------- Begin Code ---------------------
    # Identify height of 6.5 g/kg moisture level
    z_lev = jnp.argmax(rtm <= 6.5e-3, axis=1)
    if clubb_at_least_debug_level(2):
        invalid_level = (z_lev == 0) | ~jnp.any(rtm <= 6.5e-3, axis=1)
        err_info = err_info.set_fatal(mask=(time >= time_initial + 5400.0) & invalid_level)
    z_inversion = jnp.take_along_axis(gr.zt, (z_lev - 1)[:, None], axis=1)[:, 0]
    z_inversion_2d = z_inversion[:, None]

    # Large scale subsidence
    wm_zt_active = jnp.where(
        (gr.zt > 0.0) & (gr.zt <= z_inversion_2d),
        -0.0065 * gr.zt / z_inversion_2d,
        jnp.where(
            (gr.zt > z_inversion_2d) & (gr.zt <= z_inversion_2d + 300.0),
            -0.0065 * (1.0 - (gr.zt - z_inversion_2d) / 300.0),
            0.0,
        ),
    )
    thlm_forcing_active, rtm_forcing_active = calc_forcings(ngrdcol, gr, z_inversion)

    # Forcings are applied only after t = 5400 s
    active = time >= time_initial + 5400.0
    wm_zt = jnp.where(active, wm_zt_active, 0.0)
    thlm_forcing = jnp.where(active, thlm_forcing_active, 0.0)
    rtm_forcing = jnp.where(active, rtm_forcing_active, 0.0)

    wm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, wm_zt)

    # Boundary conditions.
    wm_zm = wm_zm.at[:, 0].set(0.0)
    wm_zm = wm_zm.at[:, gr.nzm - 1].set(0.0)

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=rtm.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=rtm.dtype)

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
    return (
        err_info,
        wm_zt,
        wm_zm,
        thlm_forcing,
        rtm_forcing,
        sclrm_forcing,
        edsclrm_forcing,
    )


# ======================================================================
def atex_sfclyr(ngrdcol, time, ubar, thlm_sfc, rtm_sfc, exner_sfc):
    """This subroutine computes surface fluxes of heat and moisture according
    to GCSS ATEX specifications.

    References:
        B. Stevens et al., 2000: Simulations of trade-wind cumuli
        under a strong inversion, J. Atmos. Sci, 58,1870-1891.
        http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
        ATEX_JAS2001.pdf
    """
    # -----------------BEGIN CODE-----------------------
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

    C_10 = jnp.full(ngrdcol, 0.0013)
    adjustment = jnp.full(ngrdcol, 0.0198293)
    ustar = jnp.full(ngrdcol, 0.3)
    T_sfc = jnp.full(ngrdcol, T_sfc_interp)

    # Compute wpthlp_sfc and wprtp_sfc
    wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, C_10, ubar, thlm_sfc, T_sfc, exner_sfc)
    wprtp_sfc = compute_wprtp_sfc(ngrdcol, C_10, ubar, rtm_sfc, adjustment)
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["atex_tndcy", "atex_sfclyr"]

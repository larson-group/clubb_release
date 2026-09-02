"""Contains subroutines for the DYCOMS II RF02 case."""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# ----------------------------------------------------------------------
def dycoms2_rf02_tndcy(
    ngrdcol,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    gr,
    wm_zt,
    wm_zm,
):
    """Compute thlm_ls, rtm_ls and adjust subsidence as needed.

    References:
        "A single-column model intercomparison of a heavily drizzling
        stratocumulus-topped boundary layer" Wyant, Matthew C., et al., (2007)
        J, Geophys. Res., 112, D24204,
        ftp://eos.atmos.washington.edu/pub/breth/papers/2007/GCSS-DYCOMS2-SCM.pdf

        ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
        Stevens, Bjorn, et al., (2003)
        Bull. Amer. Meteorol. Soc., 84, 579-593.
        http://www.atmos.ucla.edu/~bstevens/Documents/dycoms.pdf
    """
    # --------------------- Begin Code ---------------------
    # Enter the final thlm and rtm tendency

    # The driver initializes these inout arrays on the host; JAX functional
    # updates require device arrays at this boundary.
    wm_zt = jnp.asarray(wm_zt)
    wm_zm = jnp.asarray(wm_zm)

    # Imposed large-scale subsidence at the uppermost level.
    # CLUBB used a "one-sided" derivative method to compute mean advection at
    # the uppermost thermodynamic level. In order to avoid bringing in large
    # amounts of various quantities from above the top of the domain, set wm_zt
    # to 0 at level gr%nz. To stay consistent, set wm_zm to 0 at level gr%nz.
    wm_zt = wm_zt.at[:, gr.nzt - 1].set(0.0)
    wm_zm = wm_zm.at[:, gr.nzm - 1].set(0.0)

    thlm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=wm_zt.dtype)
    rtm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=wm_zt.dtype)
    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=wm_zt.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=wm_zt.dtype)

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
    return wm_zt, wm_zm, thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing


# ----------------------------------------------------------------------
def dycoms2_rf02_sfclyr(ngrdcol, time):
    """This subroutine computes surface fluxes of heat and moisture according
    to GCSS DYCOMS II RF 02 specifications.

    References:
        ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
        Stevens, Bjorn, et al., (2003)
        Bull. Amer. Meteorol. Soc., 84, 579-593.
        http://www.atmos.ucla.edu/~bstevens/Documents/dycoms.pdf
    """
    rho_sfc = 1.21
    # ------------------------BEGIN CODE-----------------------------------
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

    # Declare the value of ustar.
    ustar = jnp.full(ngrdcol, 0.25)
    wpthlp_sfc = jnp.full(ngrdcol, convert_sens_ht_to_km_s(sens_ht, rho_sfc))
    wprtp_sfc = jnp.full(ngrdcol, convert_latent_ht_to_m_s(latent_ht, rho_sfc))
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["dycoms2_rf02_tndcy", "dycoms2_rf02_sfclyr"]

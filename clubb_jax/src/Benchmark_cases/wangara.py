"""References:
    See below.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import sec_per_day


# ----------------------------------------------------------------------
def wangara_tndcy(ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx):
    """Subroutine to set theta and water tendencies for Wangara case.

    References:
        ``A PDF-Based Model for Boundary Layer Clouds. Part II:
        Model results'' Golaz, et al. (2002)
        JAS, Vol. 59, pp. 3552--3571.
    """
    wm_zm = jnp.zeros((ngrdcol, gr.nzm), dtype=gr.zt.dtype)

    # No large-scale subsidence for now
    wm_zt = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)

    # No large-scale water tendency or cooling
    rtm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)
    thlm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)

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


# ----------------------------------------------------------------------
def wangara_sfclyr(ngrdcol, time):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture for Wangara day 33.

    References:
        ``A PDF-Based Model for Boundary Layer Clouds. Part II:
        Model results'' Golaz, et al. (2002)
        JAS, Vol. 59, pp. 3552--3571.
    """
    # ---------------------BEGIN CODE-------------------------
    # Compute UTC time of the day in seconds
    time_utc = time % sec_per_day

    # Now convert UTC time to Australia EST (local time)
    est_offset = 36000.0
    time_est = (time_utc + est_offset) % sec_per_day

    # TODO(port-mirror): Compiled JAX calls cannot raise the source ERROR STOP.
    # Eager calls retain it; remove this split when the driver carries a
    # functional diagnostic status out of the forcing kernel.
    if not isinstance(time_est, jax.core.Tracer) and (
        time_est < 27000.0 or time_est > 63000.0
    ):
        raise ValueError(
            "wangara_sfclyr: error local time must be between 0730 and 1730; "
            f"time_est = {time_est}"
        )

    # Declare the value of ustar.
    ustar = jnp.full(ngrdcol, 0.13)

    # Compute heat and moisture fluxes
    wpthlp_sfc = jnp.full(ngrdcol, 0.18 * jnp.cos((time_est - 45000.0) / 36000.0 * jnp.pi))
    wprtp_sfc = 1.3e-4 * wpthlp_sfc
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["wangara_tndcy", "wangara_sfclyr"]

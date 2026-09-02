"""Contains subroutines for the GCSS BOMEX case.

References:
    http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.spec_hum_to_mixing_ratio import (
    flux_spec_hum_to_mixing_ratio,
    force_spec_hum_to_mixing_ratio,
)
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# ----------------------------------------------------------------------
def bomex_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, rtm):
    """Subroutine to set theta and water tendencies for BOMEX case.

    References:
        http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html
    """
    # --------------------- Begin Code ---------------------
    # Large scale advective moisture tendency
    # The BOMEX specifications give large-scale advective moisture tendency in
    # terms of total water specific humidity.then
    qtm_forcing = jnp.where(
        (gr.zt >= 0.0) & (gr.zt < 300.0),
        -1.2e-8,
        jnp.where(
            (gr.zt >= 300.0) & (gr.zt < 500.0),
            -1.2e-8 * (1.0 - (gr.zt - 300.0) / (500.0 - 300.0)),
            0.0,
        ),
    )

    # Radiative theta-l tendency
    thlm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=rtm.dtype)

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
    return thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing


# ----------------------------------------------------------------------
def bomex_sfclyr(ngrdcol, time, rtm_sfc):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS BOMEX specifications.

    References:
        http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html
    """
    # Compute heat and moisture fluxes
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        time_dependent_input.time_sfc_given.size,
        time_dependent_input.time_sfc_given,
    )
    wpthlp_sfc_calc = linear_interp_factor(
        time_frac,
        time_dependent_input.wpthlp_sfc_given[after_time],
        time_dependent_input.wpthlp_sfc_given[before_time],
    )

    # The BOMEX specifications give surface moisture flux in terms of total
    # water specific humidity.
    wpqtp_sfc_calc = linear_interp_factor(
        time_frac,
        time_dependent_input.wpqtp_sfc_given[after_time],
        time_dependent_input.wpqtp_sfc_given[before_time],
    )

    # Declare the value of ustar.
    ustar = jnp.full(ngrdcol, 0.28)

    wpthlp_sfc = jnp.full(ngrdcol, wpthlp_sfc_calc)
    wpqtp_sfc = jnp.full(ngrdcol, wpqtp_sfc_calc)

    # Convert flux from terms of total water specific humidity to terms of total
    # water mixing ratio.
    wprtp_sfc = flux_spec_hum_to_mixing_ratio(ngrdcol, rtm_sfc, wpqtp_sfc)
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["bomex_tndcy", "bomex_sfclyr"]

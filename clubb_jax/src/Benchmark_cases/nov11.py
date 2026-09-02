"""Contains subroutines for the Nov. 11 case.

References:
    See below.

Note: bottom point not at surface, so there is no sfc subroutine.

Note: nov11_altocu_tndcy has been marked as private, as the subroutine is now
obsolete.
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# -----------------------------------------------------------------------
def nov11_altocu_rtm_adjust(ngrdcol, gr, time, time_initial, dt, rtm):
    """This subroutine performs a one-time adjustment on the total water above
    cloud at time = 3600 seconds after the start of the simulation.

    This was moved from the nov11_altocu_tndcy subroutine, as said subroutine
    is now obsolete and no longer needs to calculate forcings.

    References:
        Larson, V. E., A. J. Smith, M. J. Falk, K. E. Kotenberg
        and J.-C. Golaz, 2006: What determines altocumulus lifetime?
        J. Geophys. Res., 111, D19207, doi:10.1029/2005JD007002.
        https://pantherfile.uwm.edu/vlarson/www/journal_articles/JGR_09_smith_clex_LES.pdf
    """
    # ---- Begin Code ----

    # -----------------------------------------------------------------------
    # SPECIAL NOV.11 CONDITION FOR TOTAL WATER ABOVE CLOUD
    # One hour after the initial time, the total water above cloud
    # is adjusted to be 0.89 of what it previously was.
    #
    # The conditional statement here is set so that if the timestep
    # is such that there is no timestep at exactly 3600.0 seconds,
    # then the operation still happnens at the first timestep and
    # only the first timestep after 3600.0 seconds.
    # -----------------------------------------------------------------------
    return jnp.where(
        (time >= time_initial + 3600.0)
        & (time < time_initial + 3600.0 + dt)
        & (gr.zt > 2900.0 + gr.zm[:, 0:1]),
        0.89 * rtm,
        rtm,
    )


# -----------------------------------------------------------------------
def nov11_altocu_read_t_dependent(time):
    """This subroutine reads in the values from the _surface.in file for this
    case.

    References:
        None.
    """

    # ---- Begin Code ----
    # interpolate T_sfc from time_dependent_input
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
    return sens_ht, latent_ht


__all__ = ["nov11_altocu_rtm_adjust", "nov11_altocu_read_t_dependent"]

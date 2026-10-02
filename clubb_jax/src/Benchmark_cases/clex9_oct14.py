"""Contains subroutines for the Clex9 Oct 14 case.

References:
    See below.

Note: bottom point not at surface, so there is no sfc subroutine.
"""

from __future__ import annotations

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# -----------------------------------------------------------------------
def clex9_oct14_read_t_dependent(time):
    """This subroutine reads in the values from the _surface.in file for
    this case.

    References:
        None.

    time is the current time [s]. Fortran output arguments sens_ht (sensible
    heat flux [W/m^2]) and latent_ht (latent heat flux [W/m^2]) are returned
    in order.
    """
    # ---- Begin Code ----

    # interpolate T_sfc from time_dependent_input
    # The shared JAX time_select returns zero-based array indices.
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time, time_dependent_input.time_sfc_given.size,
        time_dependent_input.time_sfc_given,
    )

    sens_ht = linear_interp_factor(
        time_frac, time_dependent_input.sens_ht_given[after_time],
        time_dependent_input.sens_ht_given[before_time],
    )

    latent_ht = linear_interp_factor(
        time_frac, time_dependent_input.latent_ht_given[after_time],
        time_dependent_input.latent_ht_given[before_time],
    )

    return sens_ht, latent_ht


__all__ = ["clex9_oct14_read_t_dependent"]

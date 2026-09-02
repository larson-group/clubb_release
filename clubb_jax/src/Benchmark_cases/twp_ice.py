"""Contains subroutines for the Jan. 2006 TWP_ICE case.

References:
    http://users.monash.edu.au/~ladavies/gcss.html
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# ----------------------------------------------------------------------
def twp_ice_sfclyr(
    ngrdcol,
    time,
    z,
    exner_sfc,
    thlm_sfc,
    ubar,
    rtm,
    p_sfc,
    saturation_formula,
):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ARM specifications.

    References:
        http://users.monash.edu.au/~ladavies/gcss.html
    """
    C_h_20 = 0.001094
    C_q_20 = 0.001133
    z0 = 0.00015
    standard_flux_alt = 20.0

    # ----------------------------------------------------------------------
    # interpolate T_sfc from time_dependent_input
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

    T_sfc = jnp.full(ngrdcol, T_sfc_interp)

    # Declare the value of ustar.
    ustar = jnp.full(ngrdcol, 0.3)

    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    # Cm = C_m_20 * ((log(20/z0))/(log(z/z0))) * ((log(20/z0))/(log(z/z0)))

    # (Stevens, et al. 2000, eq 3)
    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    Ch = C_h_20 * (
        (jnp.log(standard_flux_alt / z0)) / (jnp.log(z / z0))
    ) * ((jnp.log(standard_flux_alt / z0)) / (jnp.log(z / z0)))
    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    Cq = C_q_20 * (
        (jnp.log(standard_flux_alt / z0)) / (jnp.log(z / z0))
    ) * ((jnp.log(standard_flux_alt / z0)) / (jnp.log(z / z0)))
    rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)

    # Compute wpthlp_sfc and wprtp_sfc
    wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, Ch, ubar, thlm_sfc, T_sfc, exner_sfc)
    wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cq, ubar, rtm, rsat)
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["twp_ice_sfclyr"]

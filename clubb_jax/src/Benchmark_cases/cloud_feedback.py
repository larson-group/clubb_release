"""Contains subroutines for the Cloud Feedback cases.

References:
    http://atmgcm.msrc.sunysb.edu/cfmip_figs/Case_specification.html
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.constants_clubb import kappa, p0
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# ----------------------------------------------------------------------
def cloud_feedback_sfclyr(
    ngrdcol,
    time,
    sfctype,
    thlm_sfc,
    rtm_sfc,
    lowest_level,
    ubar,
    p_sfc,
    saturation_formula,
):
    """Sets up surface information for the cloud feedback case.

    References:
        http://cfmip.metoffice.com/
        http://atmgcm.msrc.sunysb.edu/cfmip_figs/Case_specification.html
    """
    # Constants
    C_h_20 = 0.001094
    C_q_20 = 0.001133
    z0 = 0.00015
    standard_flux_alt = 20.0

    # --------------BEGIN CODE---------------------
    if sfctype != 1:
        raise ValueError("Invalid value for sfctype")
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

    # Calculate exner_sfc based on p_sfc.
    exner_sfc = (p_sfc / p0) ** kappa

    # Just set ustar = 0.3
    ustar = jnp.full(ngrdcol, 0.3)

    #--------------------------------------------------------------------------------
    # Email way
    # wprtp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Lv )
    # wpthlp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Cp )

    # (Stevens, et al. 2000, eq 3)
    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    Ch = C_h_20 * (
        (jnp.log(standard_flux_alt / z0)) / (jnp.log(lowest_level / z0))
    ) * ((jnp.log(standard_flux_alt / z0)) / (jnp.log(lowest_level / z0)))
    # Modification in case lowest model level isn't at 10 m, from ATEX specification
    Cq = C_q_20 * (
        (jnp.log(standard_flux_alt / z0)) / (jnp.log(lowest_level / z0))
    ) * ((jnp.log(standard_flux_alt / z0)) / (jnp.log(lowest_level / z0)))
    rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)

    wpthlp_sfc = compute_wpthlp_sfc(
        ngrdcol,
        Ch,
        ubar,
        thlm_sfc,
        T_sfc,
        exner_sfc,
    )
    wprtp_sfc = compute_wprtp_sfc(
        ngrdcol,
        Cq,
        ubar,
        rtm_sfc,
        rsat,
    )
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["cloud_feedback_sfclyr"]

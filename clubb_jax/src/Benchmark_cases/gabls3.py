"""Contains subroutines for the GABLS3.

References:
    http://www.knmi.nl/samenw/gabls/
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.constants_clubb import grav


# -----------------------------------------------------------------------
def gabls3_sfclyr(
    ngrdcol,
    ubar,
    veg_T_in_K,
    thlm_sfc,
    rtm_sfc,
    lowest_level,
    exner_sfc,
    wpthlp_sfc,
    wprtp_sfc,
):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ATEX specifications.

    References:
        http://www.knmi.nl/samenw/gabls/
    """
    z0 = 0.15

    # Compute heat and moisture fluxes
    offset = jnp.full(ngrdcol, 9.9e-3)
    C_10 = jnp.full(ngrdcol, 0.00195)

    wpthlp_sfc = compute_wpthlp_sfc(
        ngrdcol,
        C_10,
        ubar,
        thlm_sfc,
        veg_T_in_K,
        exner_sfc,
    )
    wprtp_sfc = compute_wprtp_sfc(ngrdcol, C_10, ubar, rtm_sfc, offset)

    # Compute momentum fluxes
    wprtp_sfc = wprtp_sfc * 10.0
    veg_theta_in_K = veg_T_in_K / exner_sfc
    bflx = wpthlp_sfc * grav / veg_theta_in_K
    ustar = diag_ustar(lowest_level, bflx, ubar, z0)
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["gabls3_sfclyr"]

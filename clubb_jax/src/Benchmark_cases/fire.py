"""Contains subroutines for the GCSS FIRE case.

References:
    Moeng, C.-H., and coauthors, 1996: Simulation of a
    stratocumulus-topped PBL: Intercomparison among different
    numerical codes. Bull. Amer. Meteor. Soc., 77, 261-278.
    ftp://eos.atmos.washington.edu/pub/breth/papers/1996/GCSS1-Moeng.pdf
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_wprtp_sfc, compute_wpthlp_sfc
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


# ======================================================================
def fire_sfclyr(
    ngrdcol,
    time,
    ubar,
    p_sfc,
    thlm_sfc,
    rtm_sfc,
    exner_sfc,
    saturation_formula,
):
    """This subroutine computes surface fluxes of heat and moisture
    using aerodynamic formulas.

    References:
        Moeng, C.-H., and coauthors, 1996: Simulation of a
        stratocumulus-topped PBL: Intercomparison among different
        numerical codes. Bull. Amer. Meteor. Soc., 77, 261-278.
        ftp://eos.atmos.washington.edu/pub/breth/papers/1996/GCSS1-Moeng.pdf
    """
    # --------------BEGIN CODE---------------
    # Interpolate variables from time_dependent_input
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
    Cz = jnp.full(ngrdcol, 0.0013)
    ustar = jnp.full(ngrdcol, 0.3)
    rsat = sat_mixrat_liq(p_sfc, T_sfc, saturation_formula)

    # Compute wpthlp_sfc and wprtp_sfc
    wpthlp_sfc = compute_wpthlp_sfc(ngrdcol, Cz, ubar, thlm_sfc, T_sfc, exner_sfc)
    wprtp_sfc = compute_wprtp_sfc(ngrdcol, Cz, ubar, rtm_sfc, rsat)
    return wpthlp_sfc, wprtp_sfc, ustar, T_sfc


__all__ = ["fire_sfclyr"]

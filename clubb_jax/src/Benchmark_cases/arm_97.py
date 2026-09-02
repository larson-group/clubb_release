"""Contains subroutines for the July 26-30 1997 ARM IOP A case.

References:
    http://www.mmm.ucar.edu/gcss-wg4/gcss/case3.html
"""

from __future__ import annotations

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# ----------------------------------------------------------------------
def arm_97_sfclyr(ngrdcol, time, z, rho_sfc, thlm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ARM specifications.

    References:
        http://www.mmm.ucar.edu/gcss-wg4/gcss/case3.html
    """
    z0 = 0.035

    #----------------------------------------------------------------------
    if time_dependent_input.l_t_dependent:
        time_frac = -1.0
        before_time, after_time, time_frac = time_dependent_input.time_select(
            time,
            time_dependent_input.time_sfc_given.size,
            time_dependent_input.time_sfc_given,
        )
        heat_flx = linear_interp_factor(
            time_frac,
            time_dependent_input.sens_ht_given[after_time],
            time_dependent_input.sens_ht_given[before_time],
        )
        moisture_flx = linear_interp_factor(
            time_frac,
            time_dependent_input.latent_ht_given[after_time],
            time_dependent_input.latent_ht_given[before_time],
        )

        # Convert W/m^2 into w'thl' w'rt' units
        wpthlp_sfc = convert_sens_ht_to_km_s(heat_flx, rho_sfc)
        wprtp_sfc = convert_latent_ht_to_m_s(moisture_flx, rho_sfc)

        # Compute momentum fluxes using ARM Cu formulae
        bflx = grav / thlm_sfc * wpthlp_sfc

        # Compute ustar
        ustar = diag_ustar(z, bflx, ubar, z0)
        return wpthlp_sfc, wprtp_sfc, ustar

    return None


__all__ = ["arm_97_sfclyr"]

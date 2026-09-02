"""Contains subroutines for the March 2000 IOP ARM case.

References:
    Xie, S., et al. (2005), Simulations of midlatitude frontal
    clouds by single-column and cloud-resolving models during the
    Atmospheric Radiation Measurement March 2000 cloud intensive
    operational period, J. Geophys. Res., 110, D15S03,
    doi:10.1029/2004JD005119.
    http://www.agu.org/journals/jd/jd0506/2004JD005119/2004JD005119.pdf
"""

from __future__ import annotations

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    compute_ht_mostr_flux,
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav


# ----------------------------------------------------------------------
def arm_0003_sfclyr(ngrdcol, time, z, rho_sfc, thlm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ARM specifications.

    References:
        Xie, S., et al. (2005), Simulations of midlatitude frontal
        clouds by single-column and cloud-resolving models during the
        Atmospheric Radiation Measurement March 2000 cloud intensive
        operational period, J. Geophys. Res., 110, D15S03,
        doi:10.1029/2004JD005119.
        http://www.agu.org/journals/jd/jd0506/2004JD005119/2004JD005119.pdf
    """
    z0 = 0.035

    #----------------------------------------------------------------------

    # Compute heat and moisture fluxes from ARM data in (W/m2)
    heat_flx, moisture_flx = compute_ht_mostr_flux(
        time,
        time_dependent_input.time_sfc_given.size,
    )

    # Convert W/m^2 into w'thl' w'rt' units
    wpthlp_sfc = convert_sens_ht_to_km_s(heat_flx, rho_sfc)
    wprtp_sfc = convert_latent_ht_to_m_s(moisture_flx, rho_sfc)

    # Compute momentum fluxes using ARM Cu formulae
    bflx = grav / thlm_sfc * wpthlp_sfc

    # Compute ustar
    ustar = diag_ustar(z, bflx, ubar, z0)
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["arm_0003_sfclyr"]

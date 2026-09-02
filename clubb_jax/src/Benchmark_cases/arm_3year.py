"""Contains subroutine(s) for the ARM 3 Year case."""

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
def arm_3year_sfclyr(ngrdcol, time, z, rho_sfc, thlm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ARM specifications.

    References:
        None.
    """
    z0 = 0.035

    # -------------BEGIN CODE--------------
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


__all__ = ["arm_3year_sfclyr"]

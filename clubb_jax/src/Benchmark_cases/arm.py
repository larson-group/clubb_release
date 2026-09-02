"""Contains subroutines for the GCSS ARM case.

References:
    Brown et al., 2002:  Large-eddy simulation of the diurnal cycle
    of shallow cumulus convection over land. Quart. J. Roy.
    Meteor. Soc., 128, 1075-1093.
    http://www.atmos.washington.edu/~breth/GCSS/Brown_etal_
    DiurnalShCu_QJ2002.pdf
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    compute_ht_mostr_flux,
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav


# ----------------------------------------------------------------------
def arm_sfclyr(ngrdcol, time, z, thlm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ARM specifications.

    References:
        See module comments.
    """
    # Constants
    z0 = 0.035
    ntimes = 7

    #-----------------BEGIN CODE-------------------------

    # Compute heat and moisture fluxes from ARM data in (W/m2)
    heat_flx, moisture_flx = compute_ht_mostr_flux(time, ntimes)

    rho_sfc = 1.1

    # Convert heat_flx and moisture_flx to natural units
    heat_flx2 = convert_sens_ht_to_km_s(heat_flx, rho_sfc)
    moisture_flx2 = convert_latent_ht_to_m_s(moisture_flx, rho_sfc)

    # Compute momentum fluxes

    # Heat flux in units of (m2/s3) (needed by diag_ustar)
    bflx = grav / thlm_sfc * heat_flx2

    # Surface winds

    # Compute ustar
    ustar = diag_ustar(z, bflx, ubar, z0)

    # Assign fluxes

    wpthlp_sfc = jnp.full(ngrdcol, heat_flx2)
    wprtp_sfc = jnp.full(ngrdcol, moisture_flx2)
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["arm_sfclyr"]

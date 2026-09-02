"""Contains subroutines for the LBA case.

References:
    http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    convert_latent_ht_to_m_s,
    convert_sens_ht_to_km_s,
)
from clubb_jax.src.CLUBB_core.constants_clubb import grav, sec_per_hr


# ----------------------------------------------------------------------
def lba_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr):
    """Subroutine to set theta and water tendencies for LBA case.

    References:
        http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
    """
    # --------------------- Begin Code ---------------------
    # Large-scale temperature tendency
    thlm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)

    # Large-scale advective moisture tendency
    rtm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=gr.zt.dtype)

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=gr.zt.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=gr.zt.dtype)

    # Test scalars with thetal and rt if desired
    if sclr_dim > 0:
        if sclr_idx.iisclr_thl > 0:
            sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_thl - 1].set(thlm_forcing)
        if sclr_idx.iisclr_rt > 0:
            sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_rt - 1].set(rtm_forcing)
    if edsclr_dim > 0:
        if sclr_idx.iiedsclr_thl > 0:
            edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_thl - 1].set(
                thlm_forcing
            )
        if sclr_idx.iiedsclr_rt > 0:
            edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_rt - 1].set(
                rtm_forcing
            )
    return thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing


# ----------------------------------------------------------------------
def lba_sfclyr(ngrdcol, time_current, time_initial, z, rho_sfc, thlm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS BOMEX specifications.

    References:
        Grabowski, et al. (2005)
        http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
    """
    z0 = 0.035

    # Compute heat and moisture fluxes
    # From Table A.1.
    time = time_current - time_initial
    ft = jnp.maximum(
        0.0,
        jnp.cos(0.5 * jnp.pi * ((5.25 - time / sec_per_hr) / 5.25)),
    )

    # Known magic numbers
    wpthlp_sfc = convert_sens_ht_to_km_s(270.0 * ft**1.5, rho_sfc)
    wprtp_sfc = convert_latent_ht_to_m_s(554.0 * ft**1.3, rho_sfc)

    bflx = grav / thlm_sfc * wpthlp_sfc

    # Compute ustar
    ustar = diag_ustar(z, bflx, ubar, z0)
    return wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["lba_tndcy", "lba_sfclyr"]

"""Contains subroutines for the Moeng and Sullivan (1994) shear-driven case.

This case is modeled after George Bryan's shear-driven case in CM1.
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.sfc_flux import compute_momentum_flux


# ----------------------------------------------------------------------
def neutral_case_sfclyr(ngrdcol, time, um_sfc, vm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture for a shear-driven, neutral case.

    References:
        Moeng and Sullivan (1994).
    """
    #---- Begin code

    # Compute heat and moisture fluxes --- turn off heat flux after 3000 s
    wpthlp_sfc = jnp.full(ngrdcol, jnp.where(time > 80880.0, 0.0, 0.05))

    # No moisture in this case.
    wprtp_sfc = jnp.zeros(ngrdcol)

    # Heat flux in units of (m2/s3) (needed by diag_ustar)

    # Compute ustar
    ustar = jnp.full(ngrdcol, 0.5)

    upwp_sfc, vpwp_sfc = compute_momentum_flux(
        ngrdcol,
        um_sfc,
        vm_sfc,
        ubar,
        ustar,
    )
    return upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["neutral_case_sfclyr"]

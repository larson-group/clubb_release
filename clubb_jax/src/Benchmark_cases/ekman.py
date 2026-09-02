"""Contains subroutines for Hing Ong's Ekman layer case."""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.sfc_flux import compute_momentum_flux


# ----------------------------------------------------------------------
def ekman_sfclyr(ngrdcol, um_sfc, vm_sfc, ubar):
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture for a neutral Ekman layer case.

    References:
    """
    #---- Begin code

    # No heat and moisture fluxes
    wpthlp_sfc = jnp.zeros(ngrdcol)
    wprtp_sfc = jnp.zeros(ngrdcol)

    # Set ustar
    ustar = jnp.full(ngrdcol, 0.3)

    upwp_sfc, vpwp_sfc = compute_momentum_flux(
        ngrdcol,
        um_sfc,
        vm_sfc,
        ubar,
        ustar,
    )
    return upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar


__all__ = ["ekman_sfclyr"]

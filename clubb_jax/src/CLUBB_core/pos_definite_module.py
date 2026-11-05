"""JAX port of `src/CLUBB_core/pos_definite_module.F90`.

Porting deviations:
- The Fortran routine mutates `field_np1`, `flux_np1`, `field_pd`, and
  `flux_pd` in-place. JAX returns those four arrays as
  `(field_np1, flux_lim, field_pd, flux_pd)`.
- Loops over columns and levels are vectorized with JAX array operations.
"""

from __future__ import annotations

from functools import partial

import jax

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import eps, zero, zero_threshold
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import ddzm


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def pos_definite_adj(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    dt: float,
    field_np1,
    flux_np1,
    field_n,
):
    """Applies a  flux conservative positive definite scheme to a variable

    References:
    ``A Positive Definite Advection Scheme Obtained by
      Nonlinear Renormalization of the Advective Fluxes'' Smolarkiewicz (1989)
      Monthly Weather Review, Vol. 117, pp. 2626--2632

    Returns:
        (field_np1, flux_lim, field_pd, flux_pd)

    Fortran declarations (pos_definite_adj in src/CLUBB_core/pos_definite_module.F90):

        ! -------------------- Input variables --------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        real( kind = core_rknd ), intent(in) :: &
          dt ! Timestep    [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
          field_n ! The field (e.g. rtm) at n, prior to n+1

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
          flux_pd     ! Budget of the change in the flux term due to the scheme

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
          field_pd    ! Budget of the change in the mean term due to the scheme

        ! Output Variables

        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) :: &
          field_np1   ! Field at n+1 (e.g. rtm in [kg/kg])

        real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) :: &
          flux_np1    ! Flux applied to field
    """
    if clubb_at_least_debug_level(1):
        jax.debug.print("Correcting flux", ordered=True)

    # Def. of F+ and F- from eqn 2 Smolarkowicz
    flux_plus = jnp.maximum(zero_threshold, flux_np1)
    flux_minus = -jnp.minimum(zero_threshold, flux_np1)

    dz_over_dt = (1.0 / (gr.grid_dir * gr.invrs_dzt)) / dt

    if gr.grid_dir_indx > 0:
        # Ascending grid
        kphalf = jnp.minimum(jnp.arange(nzt) + 1, nzm - 1)  # k+1/2 flux level
        kmhalf = jnp.maximum(jnp.arange(nzt), 0)            # k-1/2 flux level
    else:
        # Descending grid
        kphalf = jnp.maximum(jnp.arange(nzt), 0)            # k+1/2 flux level
        kmhalf = jnp.minimum(jnp.arange(nzt) + 1, nzm - 1)  # k-1/2 flux level

    # Eqn A4 from Smolarkowicz
    # We place a limiter of eps to prevent a divide by zero, and
    #   after this calculation fout is on the scalar level, and
    #   fout is the total outward flux for the scalar level k.
    fout = jnp.maximum(
        jnp.take(flux_plus, kphalf, axis=1) + jnp.take(flux_minus, kmhalf, axis=1),
        eps,
    )

    interior_idx = jnp.arange(1, nzm - 1)
    # FIXME:
    # We haven't tested this for negative values at the nz level
    # -dschanen 13 June 2008
    if gr.grid_dir_indx > 0:
        # Ascending grid
        k_zt = interior_idx - 1  # k scalar level
        kp1 = interior_idx       # k+1 scalar level
    else:
        # Descending grid
        k_zt = interior_idx      # k scalar level
        kp1 = interior_idx - 1   # k+1 scalar level

    # Eqn 10 from Smolarkowicz (1989)
    upper = (
        (jnp.take(flux_plus, interior_idx, axis=1) / jnp.take(fout, k_zt, axis=1))
        * jnp.take(field_n, k_zt, axis=1)
        * jnp.take(dz_over_dt, k_zt, axis=1)
    )
    lower = -(
        (jnp.take(flux_minus, interior_idx, axis=1) / jnp.take(fout, kp1, axis=1))
        * jnp.take(field_n, kp1, axis=1)
        * jnp.take(dz_over_dt, k_zt, axis=1)
    )
    flux_lim_int = jnp.maximum(
        jnp.minimum(jnp.take(flux_np1, interior_idx, axis=1), upper),
        lower,
    )

    # Boundary conditions
    flux_lim = jnp.concatenate(
        [flux_np1[:, :1], flux_lim_int, flux_np1[:, nzm - 1:nzm]],
        axis=1,
    )

    # Only set flux_pd for a column if there is a below zero value in that column
    l_negative_before = jnp.any(field_np1 < zero, axis=1, keepdims=True)
    flux_pd = jnp.where(l_negative_before, (flux_lim - flux_np1) / dt, zero)

    field_nonlim = field_np1

    # Apply change to field at n+1
    field_np1 = -dt * ddzm(nzm, nzt, ngrdcol, gr, flux_lim - flux_np1) + field_np1

    # Determine the total time tendency in field due to this calculation
    # (for diagnostic purposes)
    l_negative_after = jnp.any(field_np1 < zero, axis=1, keepdims=True)
    # Only set flux_pd for a column if there is a below zero value in that column
    field_pd = jnp.where(l_negative_after, (field_np1 - field_nonlim) / dt, zero)

    # Replace the non-limited flux with the limited flux
    return field_np1, flux_lim, field_pd, flux_pd

"""JAX port of `src/CLUBB_core/sponge_layer_damping.F90`.

Description:
  This module is used for damping variables in upper altitudes of the grid.

References:
  None

Porting deviations:
- The Fortran module-level `sponge_damp_settings` and `sponge_damp_profile`
  derived types are not represented here. JAX callers pass the profile arrays
  and layer depth directly.
- `finalize_tau_sponge_damp_api` is omitted because JAX/Python does not
  manually deallocate the tau profile.
- `sponge_damp_xm` uses `tau = inf` below the sponge layer, produced by
  `initialize_tau_sponge_damp`, so the vectorized formula is a no-op below the
  layer instead of using the Fortran top-down loop exit.
"""
from __future__ import annotations

import numpy as np
import jax.numpy as jnp


def initialize_tau_sponge_damp(z, dt, zm_top, tau_min, tau_max, sponge_damp_depth):
    """Initialize time scale, tau_sponge_damp, used for damping.

    Description:
      Initialize time scale, tau_sponge_damp, used for damping.  The time scale
      attains its maximum value, tau_sponge_damp_max, at the bottom of the
      "sponge" damping layer, which results in minimal damping.  Likewise, the
      time scale attains its minimum value, tau_sponge_damp_min, at the top of
      the model, which results in maximum damping.  At levels in-between the top
      of the model and the base of the sponge damping layer, the value of
      tau_sponge_damp is in-between tau_sponge_damp_min and tau_sponge_damp_max,
      as calculated by an interpolation formula.

    References:
      None
    """
    nz = len(z)

    # Calculate the depth of the sponge layer.
    # The height of the model top is gr%zm(1,gr%k_ub_zm).
    sponge_layer_depth = sponge_damp_depth * zm_top

    # Check the value of tau_sponge_damp_min.
    if tau_min < 2.0 * dt:
        raise ValueError("tau_sponge_damp_min is too small (< 2*dt)")
    tau = np.full(nz, np.inf, dtype=np.float64)
    ratio = tau_max / tau_min
    # Calculate the value of the damping time scale, tau_sponge_damp, at levels
    # that are within the sponge damping layer.
    for k in range(nz - 1, -1, -1):          # top down (ascending grid)
        # The height of the model top is gr%zm(1,gr%k_ub_zm).
        d = zm_top - z[k]
        if d < sponge_layer_depth:
            # Vince Larson added code to use standard linear interpolation.
            # Brian Griffin reverted the linear interpolation in order to use code
            # that is similar to what is found in SAM LES.
            exponent = d / sponge_layer_depth
            tau[k] = tau_min * ratio ** exponent
            # End Vince Larson's change
            # End Brian Griffin's rebellious reversion.
        else:
            # Below sponge damping layer; exit loop.
            break
    return tau, sponge_layer_depth


def sponge_damp_xm(xm, xm_ref, z, zm_top, tau, sponge_layer_depth, dt):
    """Damps specified mean field toward a reference profile.

    Description:
      Damps specified mean field toward a reference profile.  The module must be
      initialized for this function to work.  Otherwise a stop is issued.

    References:
      None

    JAX note:
      `tau` is `inf` outside the sponge layer, so `dt/tau == 0` there and the
      formula is a no-op. `z`, `zm_top`, and `sponge_layer_depth` are kept only
      for signature compatibility.
    """
    del z, zm_top, sponge_layer_depth
    #  "Sponge"-layer damping at the domain top region
    dt_on_tau = dt / tau
    # Really, we should be using xm_ref at time n+1 rather than n.
    # However, for steady profiles of xm_ref, it won't matter.
    return (xm + dt_on_tau * xm_ref) / (1.0 + dt_on_tau)


def sponge_damp_xp2(dt, zm, xp2, x_tol_sqd, tau, sponge_layer_depth):
    """Calculate the effects of "sponge"-layer damping on the variance of x.

    Description:
      Calculates the effects of "sponge"-layer damping on the variance of x,
      xp2.

      Sponge damping on a local value of x is given by the equation:

      x_d = x - ( delta_t / tau ) * ( x - <x> ),

      where x is the local value prior to damping, x_d is the damped local value
      of x, <x> is the grid-level mean value of x, delta_t is the model time
      step duration, and tau is the damping time scale.  Since delta_t / tau has
      the same value everywhere at a grid level, the grid-level mean of x does
      not change as a result of damping.

      Subtracting <x> from both sides:

      x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),

      which results in:

      x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).

      Squaring both sides:

      ( x_d - <x> )^2 = ( 1 - delta_t / tau )^2 * ( x - <x> )^2.

      After averaging both sides, the damped value of xp2 is:

      < x_d'^2 > = ( 1 - delta_t / tau )^2 * < x'^2 >.

      Any sponge damping is applied to (predictive) xp2 after the field has been
      advanced in time.  This allows sponge damping to be applied in an implicit
      manner.  The damped value of xp2 must also be limited at a minimum value
      of x_tol^2.

    References:
    """
    zm = jnp.asarray(zm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    tau = jnp.asarray(tau, dtype=jnp.float64)
    zm_top = zm[:, -1:]
    in_sponge = (zm_top - zm) < sponge_layer_depth
    # Calculate the damped value of <x'^2>.
    # The damped value of <x'^2> needs to be greater than or equal to
    # x_tol^2.
    damped = jnp.maximum((1.0 - dt / tau) ** 2 * xp2, x_tol_sqd)
    return jnp.where(in_sponge, damped, xp2)


def sponge_damp_xp3(dt, z, zm, xp3, tau, sponge_layer_depth):
    """Calculate the effects of "sponge"-layer damping on xp3.

    Description:
      Calculates the effects of "sponge"-layer damping on xp3.

      Sponge damping on a local value of x is given by the equation:

      x_d = x - ( delta_t / tau ) * ( x - <x> ),

      where x is the local value prior to damping, x_d is the damped local value
      of x, <x> is the grid-level mean value of x, delta_t is the model time
      step duration, and tau is the damping time scale.  Since delta_t / tau has
      the same value everywhere at a grid level, the grid-level mean of x does
      not change as a result of damping.

      Subtracting <x> from both sides:

      x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),

      which results in:

      x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).

      Taking both sides to the third power:

      ( x_d - <x> )^3 = ( 1 - delta_t / tau )^3 * ( x - <x> )^3.

      After averaging both sides, the damped value of xp3 is:

      < x_d'^3 > = ( 1 - delta_t / tau )^3 * < x'^3 >.

      Any sponge damping is applied to (predictive) xp3 after the field has been
      advanced in time.  This allows sponge damping to be applied in an implicit
      manner.

    References:
    """
    z = jnp.asarray(z, dtype=jnp.float64); zm = jnp.asarray(zm, dtype=jnp.float64)
    xp3 = jnp.asarray(xp3, dtype=jnp.float64); tau = jnp.asarray(tau, dtype=jnp.float64)
    zm_top = zm[:, -1:]
    in_sponge = (zm_top - z) < sponge_layer_depth
    # Calculate the damped value of <x'^3>.
    damped = (1.0 - dt / tau) ** 3 * xp3
    return jnp.where(in_sponge, damped, xp3)

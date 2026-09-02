"""JAX port of exposed routines from ``src/CLUBB_core/fill_holes.F90``.

Porting deviations:
  * Fortran mutates ``field``, ``wp2``, ``up2``, and ``vp2`` in place.  JAX
    returns updated arrays.
  * Python callers pass zero-based ``lower_hf_level`` and ``upper_hf_level``,
    matching ``clubb_python.CLUBB_core.fill_holes``.
  * Only ``global_fill`` and ``sliding_window`` are implemented for
    ``fill_holes_vertical``.  The Fortran ``widening_windows``,
    ``smart_window``, ``smart_window_smooth``, and ``parallel_fill`` methods
    remain unported and deliberately raise ``NotImplementedError``.
  * Fortran diagnostic printing and debug-only conservation warnings are not
    reproduced in JAX.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    num_hf_draw_points,
    one,
    zero,
)
from clubb_jax.src.CLUBB_core.model_flags import global_fill, sliding_window

_F64_EPS = jnp.finfo(jnp.float64).eps


@partial(
    jax.jit,
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_global(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    field,
):
    """This subroutine clips values of 'field' that are below 'threshold' using
    the whole range [lower_hf_level:upper_hf_level] as the fill window. This
    maximized effectiveness, but minimized locality.

    Mass is conserved by reducing the clipped field everywhere by a constant
    multiplicative coefficient.

    This subroutine does not guarantee that the clipped field will exceed
    threshold everywhere; blunt clipping is needed for that.
    """
    del ngrdcol
    # --------------------- Begin Code ---------------------

    k_idx = jnp.arange(nz)[None, :]
    in_range = (k_idx >= lower_hf_level) & (k_idx <= upper_hf_level)
    rho_ds_dz = rho_ds * dz

    # Compute the numerator and denominator integrals
    numer_integral_global = jnp.sum(jnp.where(in_range, rho_ds_dz * field, 0.0), axis=1, keepdims=True)
    denom_integral_global = jnp.sum(jnp.where(in_range, rho_ds_dz, 0.0), axis=1, keepdims=True)

    # Find the vertical average of field, using the precomputed numerator and denominator,
    # see description of the vertical_avg function in advance_helper_module
    field_avg_global = numer_integral_global / denom_integral_global

    # Clip small or negative values from field.
    field_clipped = jnp.where(
        field_avg_global >= threshold,
        jnp.maximum(threshold, field),
        jnp.minimum(threshold, field),
    )

    # To compute the clipped field's vertical integral we only need to recompute the numerator
    numer_integral_clipped = jnp.sum(
        jnp.where(in_range, rho_ds_dz * field_clipped, 0.0),
        axis=1,
        keepdims=True,
    )
    field_clipped_avg = numer_integral_clipped / denom_integral_global

    safe_to_scale = (
        jnp.abs(field_clipped_avg - threshold)
        > jnp.abs(field_clipped_avg + threshold) * eps / 2.0
    )
    mass_fraction_global = (field_avg_global - threshold) / (
        field_clipped_avg - threshold
    )
    # Calculate normalized, filled field
    field_filled = threshold + mass_fraction_global * (field_clipped - threshold)

    # Do not complete calculations or update field values for this
    # column if there are no holes that need filling
    any_hole = jnp.any(jnp.where(in_range, field < threshold, False), axis=1, keepdims=True)
    apply_fill = in_range & any_hole & safe_to_scale
    return jnp.where(apply_fill, field_filled, field)


@partial(
    jax.jit,
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_sliding_window(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    field,
):
    """This subroutine clips values of 'field' that are below 'threshold' as much
    as possible (i.e. "fills holes"), but conserves the total integrated mass
    of 'field'.  This prevents clipping from acting as a spurious source.

    This performs a sliding window technique, modifying consecutive ranges of
    vertical levels in serial, this is computationally expensive, but highly local.
    This high locally has a tradeoff with effectiveness, and can often fail to fill
    all the holes, especially if there is more than ~5, as a result, this relies on
    the global fill if the first pass of the sliding window fails to fill all holes

    Mass is conserved by reducing the clipped field everywhere by a constant
    multiplicative coefficient.

    References:
      ``Numerical Methods for Wave Equations in Geophysical Fluid
        Dynamics'', Durran (1999), p. 292.
    """
    del nz
    # --------------------- Begin Code ---------------------

    rho_ds_dz = rho_ds * dz
    window_len = 2 * num_hf_draw_points + 1
    start_indx = lower_hf_level + num_hf_draw_points
    stop_indx = upper_hf_level - num_hf_draw_points + 1

    def fill_one_window(k, field_carry):
        k_start = k - num_hf_draw_points
        field_window = jax.lax.dynamic_slice(
            field_carry,
            (0, k_start),
            (ngrdcol, window_len),
        )
        rho_window = jax.lax.dynamic_slice(
            rho_ds_dz,
            (0, k_start),
            (ngrdcol, window_len),
        )

        invrs_denom_integral = one / jnp.sum(rho_window, axis=1, keepdims=True)
        field_avg = jnp.sum(rho_window * field_window, axis=1, keepdims=True) * invrs_denom_integral

        field_clipped = jnp.where(
            field_avg >= threshold,
            jnp.maximum(threshold, field_window),
            jnp.minimum(threshold, field_window),
        )
        # Compute the clipped field's vertical integral.
        # clipped_total_mass >= original_total_mass,
        # see description of the vertical_avg function in advance_helper_module
        field_clipped_avg = (
            jnp.sum(rho_window * field_clipped, axis=1, keepdims=True)
            * invrs_denom_integral
        )

        # Avoid divide by zero issues by doing nothing if field_clipped_avg ~= threshold
        safe_to_scale = (
            jnp.abs(field_clipped_avg - threshold)
            > jnp.abs(field_clipped_avg + threshold) * eps / 2.0
        )
        # Compute coefficient that makes the clipped field have the same mass as the
        # original field.  We should always have mass_fraction > 0.
        mass_fraction = (field_avg - threshold) / (field_clipped_avg - threshold)
        # Calculate normalized, filled field
        field_window_filled = threshold + mass_fraction * (field_clipped - threshold)

        any_hole = jnp.any(field_window < threshold, axis=1, keepdims=True)
        field_window_out = jnp.where(
            any_hole & safe_to_scale,
            field_window_filled,
            field_window,
        )
        return jax.lax.dynamic_update_slice(field_carry, field_window_out, (0, k_start))

    field = jax.lax.fori_loop(start_indx, stop_indx, fill_one_window, field)

    # Check if all holes were filled.

    # If the first sliding window pass didn't work, fallback to global fill
    return jax.lax.cond(
        jnp.any(field < threshold),
        lambda f: fill_holes_global(
            nz=field.shape[1],
            ngrdcol=ngrdcol,
            threshold=threshold,
            lower_hf_level=lower_hf_level,
            upper_hf_level=upper_hf_level,
            dz=dz,
            rho_ds=rho_ds,
            field=f,
        ),
        lambda f: f,
        field,
    )


@partial(
    jax.jit,
    static_argnames=(
        "nz",
        "ngrdcol",
        "lower_hf_level",
        "upper_hf_level",
        "grid_dir_indx",
        "fill_holes_type",
    ),
)
def fill_holes_vertical(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    grid_dir_indx: int,
    fill_holes_type: int,
    field,
):
    """This subroutine calls a hole filling method, specified by fill_holes_type.

    The lowest level (k=1) should not be included, as the hole-filling scheme
    should not alter the set value of 'field' at the surface (for momentum
    level variables), or consider the value of 'field' at a level below the
    surface (for thermodynamic level variables).

    For momentum level variables only, the hole-filling scheme should not
    alter the set value of 'field' at the upper boundary level (k=nz).
    So for momemtum level variables, call with upper_hf_level=nz-1, and
    for thermodynamic level variables, call with upper_hf_level=nz.
    """
    if grid_dir_indx not in (1, -1):
        raise ValueError(f"Unsupported grid_dir_indx={grid_dir_indx}")

    if grid_dir_indx == -1:
        field_reversed = jnp.flip(field, axis=1)
        dz_reversed = jnp.flip(dz, axis=1)
        rho_ds_reversed = jnp.flip(rho_ds, axis=1)
        lower_reversed = nz - 1 - lower_hf_level
        upper_reversed = nz - 1 - upper_hf_level
        filled = fill_holes_vertical(
            nz,
            ngrdcol,
            threshold,
            lower_reversed,
            upper_reversed,
            dz_reversed,
            rho_ds_reversed,
            1,
            fill_holes_type,
            field_reversed,
        )
        return jnp.flip(filled, axis=1)

    # Only bother will a fill call if there are values below threshold
    if fill_holes_type == global_fill:
        # This fills holes by modifying the entire range, this is maximally effective
        # and computationally cheap, but minimally local
        filled = fill_holes_global(
            nz,
            ngrdcol,
            threshold,
            lower_hf_level,
            upper_hf_level,
            dz,
            rho_ds,
            field,
        )
    elif fill_holes_type == sliding_window:
        # This performs a sliding window technique, modifying consecutive ranges of
        # vertical levels in serial, this is computationally expensive, but highly local.
        # This can also fail to fill, so this falls back to a global fill if neccesary.
        filled = fill_holes_sliding_window(
            nz,
            ngrdcol,
            threshold,
            lower_hf_level,
            upper_hf_level,
            dz,
            rho_ds,
            field,
        )
    else:
        # TODO(JAX port): port the remaining Fortran fill_holes_type options
        # rather than routing them through a Python/API fallback.
        raise NotImplementedError(
            "JAX fill_holes_vertical currently supports global_fill and "
            f"sliding_window; got fill_holes_type={fill_holes_type}."
        )

    return jnp.where(jnp.any(field < threshold), filled, field)


@partial(
    jax.jit,
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_wp2_from_horz_tke(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    wp2,
    up2,
    vp2,
):
    """This subroutine clips values of wp2 that are below 'threshold' as much
    as possible (i.e. "fills holes"), but conserves the turbulent kinetic energy
    (up2+vp2+wp2). This prevents clipping from acting as a spurious source.

    Turbulent kinetic energy at each height level is conserved by reducing up2 and vp2
    by a multiplicative coefficient.

    This subroutine does not guarantee that the clipped field will exceed
    threshold everywhere; blunt clipping is needed for that.

    The lowest level (k=1) should not be included, as the hole-filling scheme
    should not alter the set value of 'field' at the surface (for momentum
    level variables), or consider the value of 'field' at a level below the
    surface (for thermodynamic level variables).

    For momentum level variables only, the hole-filling scheme should not
    alter the set value of 'field' at the upper boundary level (k=nz).
    So for momemtum level variables, call with upper_hf_level=nz-1, and
    for thermodynamic level variables, call with upper_hf_level=nz.
    """
    del ngrdcol
    # --------------------- Begin Code ---------------------

    k_idx = jnp.arange(nz)[None, :]
    in_range = (k_idx >= lower_hf_level) & (k_idx <= upper_hf_level)

    # For each height level, fill holes in wp2 by taking tke from up2 and vp2
    missing_wp2 = threshold - wp2
    up2_avail = jnp.maximum(up2 - threshold, zero)
    vp2_avail = jnp.maximum(vp2 - threshold, zero)
    up2_vp2_avail = up2_avail + vp2_avail
    # Check if we have a hole to fill at level k and
    # there is buffer TKE in up2 and/or vp2 available
    do_fill = in_range & (wp2 < threshold) & ((up2 > threshold) | (vp2 > threshold))

    # Not enough TKE available to fill the hole.
    case_not_enough = do_fill & (missing_wp2 >= up2_vp2_avail)
    wp2_not_enough = wp2 + up2_vp2_avail
    up2_not_enough = jnp.minimum(up2, threshold)
    vp2_not_enough = jnp.minimum(vp2, threshold)

    # Enough TKE is available to fill the hole.
    case_enough = do_fill & (missing_wp2 < up2_vp2_avail)
    no_up2_avail = jnp.abs(up2_avail) < _F64_EPS * 1000.0
    no_vp2_avail = jnp.abs(vp2_avail) < _F64_EPS * 1000.0
    # Calculate portion of up2/vp2 that we want to take away
    ratio = jnp.where(up2_vp2_avail > zero, missing_wp2 / up2_vp2_avail, zero)

    up2_enough = jnp.where(
        no_up2_avail,
        up2,
        jnp.where(
            no_vp2_avail,
            up2 - missing_wp2,
            threshold + up2_avail * (one - ratio),
        ),
    )
    vp2_enough = jnp.where(
        no_up2_avail,
        vp2 - missing_wp2,
        jnp.where(
            no_vp2_avail,
            vp2,
            threshold + vp2_avail * (one - ratio),
        ),
    )

    wp2_out = jnp.where(case_not_enough, wp2_not_enough, jnp.where(case_enough, threshold, wp2))
    up2_out = jnp.where(case_not_enough, up2_not_enough, jnp.where(case_enough, up2_enough, up2))
    vp2_out = jnp.where(case_not_enough, vp2_not_enough, jnp.where(case_enough, vp2_enough, vp2))
    return wp2_out, up2_out, vp2_out


__all__ = [
    "fill_holes_global",
    "fill_holes_sliding_window",
    "fill_holes_vertical",
    "fill_holes_wp2_from_horz_tke",
]

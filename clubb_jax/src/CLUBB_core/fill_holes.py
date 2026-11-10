"""JAX port of exposed routines from ``src/CLUBB_core/fill_holes.F90``.

Porting deviations:
  * Fortran mutates ``field``, ``wp2``, ``up2``, and ``vp2`` in place.  JAX
    returns updated arrays.
  * Python callers pass zero-based ``lower_hf_level`` and ``upper_hf_level``,
    matching ``clubb_python.CLUBB_core.fill_holes``.
  * Method dispatch rejects an unknown ``fill_holes_type`` even when there are
    no holes; Fortran only dispatches after finding a below-threshold value.
"""

from __future__ import annotations

from functools import partial

import jax
import sys

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import core_rknd
from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    num_hf_draw_points,
    one,
    zero,
)
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.model_flags import (
    global_fill,
    parallel_fill,
    sliding_window,
    smart_window,
    smart_window_smooth,
    widening_windows,
)

_CORE_EPS = jnp.finfo(core_rknd).eps


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

    Fortran declarations (fill_holes_global in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          dz      ! Spacing between thermodynamic grid levels; centered over
                  ! momentum grid levels
                  ! OR
                  ! Spcaing between momentum grid levels; centered over
                  ! thermodynamic grid levels

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level, & ! Upper grid level of global hole-filling range      []
          grid_dir_indx     ! Grid direction index (+1 ascending; -1 descending)

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                     ! fall                           [Units vary; same as field]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
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
    # Guard the division itself: masking an unused 0/0 below still gives NaN
    # reverse-mode derivatives. Keep the original denominator when scaling.
    mass_fraction_global = (field_avg_global - threshold) / jnp.where(
        safe_to_scale, field_clipped_avg - threshold, one
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

    Fortran declarations (fill_holes_sliding_window in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          dz      ! Spacing between thermodynamic grid levels; centered over
                  ! momentum grid levels
                  ! OR
                  ! Spcaing between momentum grid levels; centered over
                  ! thermodynamic grid levels

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level, & ! Upper grid level of global hole-filling range      []
          grid_dir_indx     ! Grid direction index (+1 ascending; -1 descending)

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                     ! fall                           [Units vary; same as field]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
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
        # Inactive windows must have finite derivatives too; the output mask
        # alone cannot prevent an unused 0/0 from contaminating backpropagation.
        mass_fraction = (field_avg - threshold) / jnp.where(
            safe_to_scale, field_clipped_avg - threshold, one
        )
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
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_widening_windows(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    field,
):
    """This performs a multi-window technique, splitting the vertical levels into a
    series of ranges that can be filled in parallel. This is much cheaper
    (computationally) than the sliding window and still results in a high degree of
    locality, but is even less effective. If the windowing technique fails, we
    repeat the algorithm with increasingly larger window sizes rather than fallback
    to a global fill. This increases our chances of being able to fill locally, but
    if filling continuously fails, the window size grows to the entire range,
    effectively performing a global fill.

    Fortran declarations (fill_holes_widening_windows in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          dz      ! Spacing between thermodynamic grid levels; centered over
                  ! momentum grid levels
                  ! OR
                  ! Spcaing between momentum grid levels; centered over
                  ! thermodynamic grid levels

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level, & ! Upper grid level of global hole-filling range      []
          grid_dir_indx     ! Grid direction index (+1 ascending; -1 descending)

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                     ! fall                           [Units vary; same as field]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
    """
    del nz
    # --------------------- Begin Code ---------------------

    rho_ds_dz = rho_ds * dz

    # Fortran repeats the pass in a do while loop until no holes
    # remain or one window covers the whole fill range. Each window size sets the
    # window shapes, so JAX unrolls the passes over these static window sizes.
    num_levels = upper_hf_level - lower_hf_level + 1
    window_sizes = [10]
    while window_sizes[-1] < num_levels:
        window_sizes.append(2 * window_sizes[-1])

    l_field_below_threshold = jnp.asarray(True)

    for window_size in window_sizes:

        # Group the range [lower_hf_level:upper_hf_level] into windows of window_size
        # levels, starting at lower_hf_level. The last window may be shorter, so pad
        # it with levels that have no mass and no holes.
        num_windows = -(-num_levels // window_size)
        num_pad = num_windows * window_size - num_levels
        field_windows = jnp.pad(
            field[:, lower_hf_level:upper_hf_level + 1],
            ((0, 0), (0, num_pad)),
            constant_values=threshold,
        ).reshape(ngrdcol, num_windows, window_size)
        rho_ds_dz_windows = jnp.pad(
            rho_ds_dz[:, lower_hf_level:upper_hf_level + 1],
            ((0, 0), (0, num_pad)),
        ).reshape(ngrdcol, num_windows, window_size)

        l_window_below_threshold = jnp.any(field_windows < threshold, axis=2, keepdims=True)

        invrs_denom = one / jnp.sum(rho_ds_dz_windows, axis=2, keepdims=True)

        field_avg = jnp.sum(rho_ds_dz_windows * field_windows, axis=2, keepdims=True) * invrs_denom

        # field_avg < threshold
        field_clipped_avg_min = (
            jnp.sum(rho_ds_dz_windows * jnp.minimum(threshold, field_windows), axis=2, keepdims=True)
            * invrs_denom
        )

        # Both branches are evaluated, so guard the division of the unselected branch:
        # there, a window without holes has field_clipped_avg_min - threshold = 0.
        mass_fraction_min = (field_avg - threshold) / jnp.where(
            field_avg < threshold, field_clipped_avg_min - threshold, one
        )

        field_filled_min = threshold + mass_fraction_min * (
            jnp.minimum(threshold, field_windows) - threshold
        )

        # field_avg >= threshold
        field_clipped_avg_max = (
            jnp.sum(rho_ds_dz_windows * jnp.maximum(threshold, field_windows), axis=2, keepdims=True)
            * invrs_denom
        )

        l_fill_max = jnp.abs(field_clipped_avg_max - threshold) > zero

        mass_fraction_max = (field_avg - threshold) / jnp.where(
            l_fill_max, field_clipped_avg_max - threshold, one
        )

        # Calculate normalized, filled field
        field_filled_max = threshold + mass_fraction_max * (
            jnp.maximum(threshold, field_windows) - threshold
        )

        field_windows_filled = jnp.where(
            field_avg < threshold,
            field_filled_min,
            jnp.where(l_fill_max, field_filled_max, field_windows),
        )

        # Only windows with holes are filled, and only while a previous pass
        # left holes.
        field_windows = jnp.where(
            l_field_below_threshold & l_window_below_threshold,
            field_windows_filled,
            field_windows,
        )

        field = field.at[:, lower_hf_level:upper_hf_level + 1].set(
            field_windows.reshape(ngrdcol, num_windows * window_size)[:, :num_levels]
        )

        l_field_below_threshold = jnp.any(field_windows < threshold)

    return field


@partial(
    jax.jit,
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_smart_window(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    field,
):
    """This performs filling by first scanning the vertical levels to determine a
    range to fill. This is done by looping over the vertical, and tracking the start
    and end of a range of holes, then expanding that range until the field average
    is above threshold before attempting to fill. Once a range is filled, it continues
    scanning the vertical levels for more hole ranges. This is highly local when
    possible, and is kept computationally cheap by avoiding fill attempts on ranges
    whose average is below threshold. Because the range is expanded until either the
    field average of the range is above threshold, or is the whole range, this
    will effectively perform a global fill if neccesary.

    Fortran declarations (fill_holes_smart_window in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          dz      ! Spacing between thermodynamic grid levels; centered over
                  ! momentum grid levels
                  ! OR
                  ! Spcaing between momentum grid levels; centered over
                  ! thermodynamic grid levels

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level, & ! Upper grid level of global hole-filling range      []
          grid_dir_indx     ! Grid direction index (+1 ascending; -1 descending)

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                     ! fall                           [Units vary; same as field]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
    """
    del ngrdcol
    # --------------------- Begin Code ---------------------

    # Fortran marks "no level found" with 0, which is not a
    # 1-based level. JAX levels are 0-based, so it uses -1.
    no_level = -1

    # Fortran scans with do while loops and jumps k to k_end after
    # a fill. k increases by at least one per iteration, so JAX scans with a fixed
    # number of iterations (one per level in the range) and skips iterations once k
    # passes upper_hf_level. The range expansion is limited the same way: after
    # num_expansions doublings, n_steal_points covers the whole range. Fixed trip
    # counts keep the loops reverse-mode differentiable.
    num_levels = upper_hf_level - lower_hf_level + 1
    num_expansions = num_levels.bit_length() + 1

    k_levels = jnp.arange(nz)

    # Calculate the threshold normalized mass, a useful precalculation
    normalized_mass = (field - threshold) * rho_ds * jnp.abs(dz)

    def fill_holes_smart_window_column(field, normalized_mass):
        # The nested loop bodies get field and normalized_mass from their carry,
        # which holds the values updated by earlier fills.

        def expand_range(_, carry):
            (k_start, k_end, n_steal_points, stealable_mass, hole_mass,
             k_first_lte_thresh, k_last_lte_thresh, normalized_mass) = carry

            # If there isn't enough stealable mass to fill the holes, then try to repeatedly
            # expand the range until we find enough
            l_expand = (stealable_mass < jnp.abs(hole_mass)) & (
                (k_start != lower_hf_level) | (k_end != upper_hf_level)
            )

            # Double the points we search for on either side of the range
            n_steal_points_new = 2 * n_steal_points

            # Recalcualte our fill range
            k_start_new = jnp.maximum(lower_hf_level, k_first_lte_thresh - n_steal_points_new)
            k_end_new = jnp.minimum(upper_hf_level, k_last_lte_thresh + n_steal_points_new)

            # Only add to the sum what we haven't so far, this handles stuff added below
            # and above
            l_added = ((k_levels >= k_start_new) & (k_levels < k_start)) | (
                (k_levels > k_end) & (k_levels <= k_end_new)
            )
            stealable_mass_new = stealable_mass + jnp.sum(
                jnp.where(l_added & (normalized_mass >= zero), normalized_mass, zero)
            )
            hole_mass_new = hole_mass + jnp.sum(
                jnp.where(l_added & (normalized_mass < zero), normalized_mass, zero)
            )

            # Set new range
            return (
                jnp.where(l_expand, k_start_new, k_start),
                jnp.where(l_expand, k_end_new, k_end),
                jnp.where(l_expand, n_steal_points_new, n_steal_points),
                jnp.where(l_expand, stealable_mass_new, stealable_mass),
                jnp.where(l_expand, hole_mass_new, hole_mass),
                k_first_lte_thresh,
                k_last_lte_thresh,
                normalized_mass,
            )

        def scan_level(_, carry):
            k, k_first_lte_thresh, k_last_lte_thresh, n_holes, n_steal_points, field, normalized_mass = carry

            l_active = k <= upper_hf_level
            normalized_mass_k = normalized_mass[jnp.minimum(k, upper_hf_level)]

            # If normalized_mass(i,k) > 0, then field(i,k) > threshold
            # If we found field value at or below tolerance, but no holes yet, we want to
            # ignore it by reseting the k tracker for the first <= thresh
            l_above = normalized_mass_k > zero
            k_first_lte_thresh = jnp.where(l_above & (n_holes == 0), no_level, k_first_lte_thresh)

            # If normalized_mass(i,k) <= 0, then field(i,k) <= threshold

            # This is only a true hole if it's less than zero, otherwise no need to fill it
            l_hole = normalized_mass_k < zero
            n_holes = jnp.where(l_hole, n_holes + 1, n_holes)
            n_steal_points = jnp.where(l_hole, 2 * n_holes, n_steal_points)

            # Keep track of the k we found this potential hole
            k_last_lte_thresh = jnp.where(l_above, k_last_lte_thresh, k)

            # If this is our first time seeing a <= threshold value, save it
            k_first_lte_thresh = jnp.where(
                ~l_above & (k_first_lte_thresh == no_level), k, k_first_lte_thresh
            )

            # If we have holes, and we are either out of levels or k is >= n_steal_points away from
            # the first <= threshold level, then we attempt to fill the holes we've found
            l_fill = l_active & (n_holes > 0) & (
                (k == upper_hf_level) | (jnp.abs(k - k_last_lte_thresh) >= n_steal_points)
            )

            # Calcualte our fill range
            k_start = jnp.maximum(lower_hf_level, k_first_lte_thresh - n_steal_points)
            k_end = jnp.minimum(upper_hf_level, k_last_lte_thresh + n_steal_points)

            # Calculate the excess mass we can steal, and the hole mass that needs filling
            # Fortran uses max( normalized_mass, zero ) and min( normalized_mass, zero ).
            # The JAX selections give the same values, but levels exactly at threshold
            # get the derivative of the field >= threshold branch used in the fill
            # below, rather than half of it from jnp.maximum.
            l_in_range = (k_levels >= k_start) & (k_levels <= k_end)
            stealable_mass = jnp.sum(jnp.where(l_in_range & (normalized_mass >= zero), normalized_mass, zero))
            hole_mass = jnp.sum(jnp.where(l_in_range & (normalized_mass < zero), normalized_mass, zero))

            k_start, k_end, n_steal_points_fill, stealable_mass, hole_mass, _, _, _ = jax.lax.fori_loop(
                0,
                num_expansions,
                expand_range,
                (k_start, k_end, n_steal_points, stealable_mass, hole_mass,
                 k_first_lte_thresh, k_last_lte_thresh, normalized_mass),
            )

            # We can only move mass if there is some above threshold
            l_move_mass = l_fill & (stealable_mass > zero)

            # Calcualte our mass fractions
            # Both divisions are evaluated, so guard them when no mass is moved.
            mass_fraction_above = jnp.maximum(
                (stealable_mass + hole_mass) / jnp.where(l_move_mass, stealable_mass, one), zero
            )
            mass_fraction_below = jnp.maximum(
                (stealable_mass + hole_mass) / jnp.where(l_move_mass, hole_mass, one), zero
            )

            # Compute new field and update normalized mass in case we need to touch
            # these levels again in another pass
            l_in_range = (k_levels >= k_start) & (k_levels <= k_end)
            mass_fraction = jnp.where(field >= threshold, mass_fraction_above, mass_fraction_below)
            l_update = l_move_mass & l_in_range
            field = jnp.where(l_update, mass_fraction * (field - threshold) + threshold, field)
            normalized_mass = jnp.where(l_update, mass_fraction * normalized_mass, normalized_mass)

            # Filling complete, reset in case there are more holes to find
            k_first_lte_thresh = jnp.where(l_fill, no_level, k_first_lte_thresh)
            k_last_lte_thresh = jnp.where(l_fill, no_level, k_last_lte_thresh)
            n_holes = jnp.where(l_fill, 0, n_holes)
            n_steal_points = jnp.where(l_fill, n_steal_points_fill, n_steal_points)

            # set k to whatever we used for k_end, since all holes up to k_end should be filled
            k = jnp.where(l_fill, k_end, k)

            k = jnp.where(l_active, k + 1, k)

            return k, k_first_lte_thresh, k_last_lte_thresh, n_holes, n_steal_points, field, normalized_mass

        # Initialize counters
        k_first_lte_thresh = jnp.asarray(no_level)
        k_last_lte_thresh = jnp.asarray(no_level)
        n_holes = jnp.asarray(0)
        n_steal_points = jnp.asarray(0)

        # Loop from lower_hf_level to upper_hf_level
        k = jnp.asarray(lower_hf_level)
        _, _, _, _, _, field, _ = jax.lax.fori_loop(
            0,
            num_levels,
            scan_level,
            (k, k_first_lte_thresh, k_last_lte_thresh, n_holes, n_steal_points, field, normalized_mass),
        )

        return field

    # All columns are independent
    field = jax.vmap(fill_holes_smart_window_column)(field, normalized_mass)

    return field


@partial(
    jax.jit,
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_smart_window_smooth(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    field,
):
    """This is a modification of the smart_window technique, but includes two "smoothing"
    features. The first feature is psuedo-threshold, larger than the input threshold.
    This is used to maintain the "shape" of the below threshold part of the field,
    by squishing the field with values in the range [minval(field),threshold] to
    [threshold, psuedo-threshold]. As a result, this will increase field entries
    whose values are between [threshold, psuedo-threshold], even though they are
    above threshold already. The second smoothing feature is a minimum mass_fraction,
    mf_min, which prevents the above threshold field values from decreasing by
    more than a factor of mf_min. This increased smoothness requires more mass
    in general, making it slightly less local than the non-smooth smart_window
    version, and the increased complexity makes it slightly more expensive.
    Additionally if the whole field average is between [threshold, psuedo-threshold],
    then filling can fail the same way that the a global fill can fail if the average
    is below threshold, even though the field average is above threshold - this
    only occurs in situations where the hole mass is almost as large as the total excess
    (above threshold) mass.

    Fortran declarations (fill_holes_smart_window_smooth in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          dz      ! Spacing between thermodynamic grid levels; centered over
                  ! momentum grid levels
                  ! OR
                  ! Spcaing between momentum grid levels; centered over
                  ! thermodynamic grid levels

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level, & ! Upper grid level of global hole-filling range      []
          grid_dir_indx     ! Grid direction index (+1 ascending; -1 descending)

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                      ! fall                           [Units vary; same as field]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
    """
    del ngrdcol
    # --------------------- Begin Code ---------------------

    # Fortran marks "no level found" with 0, which is not a
    # 1-based level. JAX levels are 0-based, so it uses -1.
    no_level = -1

    # Fortran scans with do while loops and jumps k to k_end after
    # a fill. k increases by at least one per iteration, so JAX scans with a fixed
    # number of iterations (one per level in the range) and skips iterations once k
    # passes upper_hf_level. The range expansion is limited the same way: after
    # num_expansions doublings, n_steal_points covers the whole range. Fixed trip
    # counts keep the loops reverse-mode differentiable.
    num_levels = upper_hf_level - lower_hf_level + 1
    num_expansions = num_levels.bit_length() + 1

    k_levels = jnp.arange(nz)

    threshold_2 = 1.05 * threshold  # known magic number
    # Fortran assigns the default-real (single precision) literal 0.05, which is
    # this value in double precision.
    mf_min = 0.05000000074505806  # known magic number

    # Calculate the threshold normalized mass, a useful precalculation
    normalized_mass = (field - threshold_2) * rho_ds * jnp.abs(dz)

    def fill_holes_smart_window_smooth_column(field, normalized_mass):
        # The nested loop bodies get field and normalized_mass from their carry,
        # which holds the values updated by earlier fills.

        def expand_range(_, carry):
            (k_start, k_end, n_steal_points, stealable_mass, hole_mass, field_min,
             k_first_lte_thresh, k_last_lte_thresh, field, normalized_mass) = carry

            # If there isn't enough stealable mass to fill the holes, then try to repeatedly
            # expand the range until we find enough
            l_expand = (stealable_mass < jnp.abs(hole_mass)) & (
                (k_start != lower_hf_level) | (k_end != upper_hf_level)
            )

            # Double the points we search for on either side of the range
            n_steal_points_new = 2 * n_steal_points

            # Recalcualte our fill range
            k_start_new = jnp.maximum(lower_hf_level, k_first_lte_thresh - n_steal_points_new)
            k_end_new = jnp.minimum(upper_hf_level, k_last_lte_thresh + n_steal_points_new)

            # Only add to the sum what we haven't so far, this handles stuff added below
            # and above
            l_added = ((k_levels >= k_start_new) & (k_levels < k_start)) | (
                (k_levels > k_end) & (k_levels <= k_end_new)
            )
            stealable_mass_new = stealable_mass + jnp.sum(
                jnp.where(l_added & (normalized_mass >= zero), normalized_mass, zero)
            )
            hole_mass_new = hole_mass + jnp.sum(
                jnp.where(l_added & (normalized_mass < zero), normalized_mass, zero)
            )
            field_min_new = jnp.minimum(field_min, jnp.min(jnp.where(l_added, field, jnp.inf)))

            # Set new range
            return (
                jnp.where(l_expand, k_start_new, k_start),
                jnp.where(l_expand, k_end_new, k_end),
                jnp.where(l_expand, n_steal_points_new, n_steal_points),
                jnp.where(l_expand, stealable_mass_new, stealable_mass),
                jnp.where(l_expand, hole_mass_new, hole_mass),
                jnp.where(l_expand, field_min_new, field_min),
                k_first_lte_thresh,
                k_last_lte_thresh,
                field,
                normalized_mass,
            )

        def scan_level(_, carry):
            k, k_first_lte_thresh, k_last_lte_thresh, n_holes, n_steal_points, field, normalized_mass = carry

            l_active = k <= upper_hf_level
            normalized_mass_k = normalized_mass[jnp.minimum(k, upper_hf_level)]

            # If normalized_mass(i,k) > 0, then field(i,k) > threshold
            # If we found field value at or below tolerance, but no holes yet, we want to
            # ignore it by reseting the k tracker for the first <= thresh
            l_above = normalized_mass_k > zero
            k_first_lte_thresh = jnp.where(l_above & (n_holes == 0), no_level, k_first_lte_thresh)

            # If normalized_mass(i,k) <= 0, then field(i,k) <= threshold

            # This is only a true hole if it's less than zero, otherwise no need to fill it
            l_hole = normalized_mass_k < zero
            n_holes = jnp.where(l_hole, n_holes + 1, n_holes)
            n_steal_points = jnp.where(l_hole, 2 * n_holes, n_steal_points)

            # Keep track of the k we found this potential hole
            k_last_lte_thresh = jnp.where(l_above, k_last_lte_thresh, k)

            # If this is our first time seeing a <= threshold value, save it
            k_first_lte_thresh = jnp.where(
                ~l_above & (k_first_lte_thresh == no_level), k, k_first_lte_thresh
            )

            # If we have holes, and we are either out of levels or k is >= n_steal_points away from
            # the first <= threshold level, then we attempt to fill the holes we've found
            l_fill = l_active & (n_holes > 0) & (
                (k == upper_hf_level) | (jnp.abs(k - k_last_lte_thresh) >= n_steal_points)
            )

            # Calcualte our fill range
            k_start = jnp.maximum(lower_hf_level, k_first_lte_thresh - n_steal_points)
            k_end = jnp.minimum(upper_hf_level, k_last_lte_thresh + n_steal_points)

            # Calculate the excess mass we can steal, and the hole mass that needs filling
            # Fortran uses max( normalized_mass, zero ) and min( normalized_mass, zero ).
            # The JAX selections give the same values, but levels exactly at threshold
            # get the derivative of the field >= threshold_2 branch used in the fill
            # below, rather than half of it from jnp.maximum.
            l_in_range = (k_levels >= k_start) & (k_levels <= k_end)
            stealable_mass = jnp.sum(jnp.where(l_in_range & (normalized_mass >= zero), normalized_mass, zero))
            hole_mass = jnp.sum(jnp.where(l_in_range & (normalized_mass < zero), normalized_mass, zero))
            field_min = jnp.min(jnp.where(l_in_range, field, jnp.inf))

            (k_start, k_end, n_steal_points_fill, stealable_mass, hole_mass, field_min,
             _, _, _, _) = jax.lax.fori_loop(
                0,
                num_expansions,
                expand_range,
                (k_start, k_end, n_steal_points, stealable_mass, hole_mass, field_min,
                 k_first_lte_thresh, k_last_lte_thresh, field, normalized_mass),
            )

            # We can only move mass if there is some above threshold
            l_move_mass = l_fill & (stealable_mass > zero)

            # Both divisions are evaluated, so guard them when no mass is moved.
            T = (threshold_2 - threshold) / jnp.where(l_move_mass, threshold_2 - field_min, one)

            # Calcualte our mass fractions
            mass_fraction_above = jnp.maximum(
                (stealable_mass + (one - T) * hole_mass)
                / jnp.where(l_move_mass, stealable_mass, one), mf_min
            )
            mass_fraction_below = jnp.maximum(
                ((one - mf_min) * stealable_mass + hole_mass)
                / jnp.where(l_move_mass, hole_mass, one), T
            )

            # Compute new field and update normalized mass in case we need to touch
            # these levels again in another pass
            l_in_range = (k_levels >= k_start) & (k_levels <= k_end)
            mass_fraction = jnp.where(field >= threshold_2, mass_fraction_above, mass_fraction_below)
            l_update = l_move_mass & l_in_range
            field = jnp.where(l_update, mass_fraction * (field - threshold_2) + threshold_2, field)
            normalized_mass = jnp.where(l_update, mass_fraction * normalized_mass, normalized_mass)

            # Filling complete, reset in case there are more holes to find
            k_first_lte_thresh = jnp.where(l_fill, no_level, k_first_lte_thresh)
            k_last_lte_thresh = jnp.where(l_fill, no_level, k_last_lte_thresh)
            n_holes = jnp.where(l_fill, 0, n_holes)
            n_steal_points = jnp.where(l_fill, n_steal_points_fill, n_steal_points)

            # set k to whatever we used for k_end, since all holes up to k_end should be filled
            k = jnp.where(l_fill, k_end, k)

            k = jnp.where(l_active, k + 1, k)

            return k, k_first_lte_thresh, k_last_lte_thresh, n_holes, n_steal_points, field, normalized_mass

        # Initialize counters
        k_first_lte_thresh = jnp.asarray(no_level)
        k_last_lte_thresh = jnp.asarray(no_level)
        n_holes = jnp.asarray(0)
        n_steal_points = jnp.asarray(0)

        # Loop from lower_hf_level to upper_hf_level
        k = jnp.asarray(lower_hf_level)
        _, _, _, _, _, field, _ = jax.lax.fori_loop(
            0,
            num_levels,
            scan_level,
            (k, k_first_lte_thresh, k_last_lte_thresh, n_holes, n_steal_points, field, normalized_mass),
        )

        return field

    # All columns are independent
    field = jax.vmap(fill_holes_smart_window_smooth_column)(field, normalized_mass)

    return field


@partial(
    jax.jit,
    static_argnames=("nz", "ngrdcol", "lower_hf_level", "upper_hf_level"),
)
def fill_holes_parallel(
    nz: int,
    ngrdcol: int,
    threshold: float,
    lower_hf_level: int,
    upper_hf_level: int,
    dz,
    rho_ds,
    field,
):
    """This fills holes by considering each hole independently. It first calculates the
    number of holes (n_holes), then only allows each hole to steal at most 1/n_holes mass
    from the above threshold levels. This maximizes parallelism because each hole can
    be handled simultaneously, and maintains a high degree of locality. However,
    handling each hole by itself can require many more large passes than is required
    in theory, making this very computationally expensive if the increaseD parallelism
    cannot be exploited (e.g bad for CPUs, this is designed for GPUs).
    Note: this seems to be slower than other methods on GPUs too, perhaps worth
    a revisit. Maybe faster for a small number of holes? Was mostly tested with
    fields with a large range of holes.

    Fortran declarations (fill_holes_parallel in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          dz      ! Spacing between thermodynamic grid levels; centered over
                  ! momentum grid levels
                  ! OR
                  ! Spcaing between momentum grid levels; centered over
                  ! thermodynamic grid levels

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level, & ! Upper grid level of global hole-filling range      []
          grid_dir_indx     ! Grid direction index (+1 ascending; -1 descending)

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                     ! fall                           [Units vary; same as field]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
          rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
    """
    # --------------------- Begin Code ---------------------

    # Fortran fills the holes of a column in parallel on GPUs, with
    # atomic updates, so the order of holes is not fixed there. JAX follows the CPU
    # Fortran, which fills the holes of a column one at a time from lower_hf_level
    # to upper_hf_level, each hole seeing the field left by the previous ones.
    # The scan over levels and the range expansion use fixed trip counts, as in
    # fill_holes_smart_window, to keep the loops reverse-mode differentiable.
    num_levels = upper_hf_level - lower_hf_level + 1
    num_expansions = num_levels.bit_length() + 1

    k_levels = jnp.arange(nz)
    l_in_hf_range = (k_levels >= lower_hf_level) & (k_levels <= upper_hf_level)

    # Precalculate rho_ds_dz and invrs_rho_ds_dz
    rho_ds_dz = rho_ds * jnp.abs(dz)
    invrs_rho_ds_dz = one / rho_ds_dz

    # Calculate the number of holes in each column
    n_holes = jnp.sum(l_in_hf_range & (field < threshold), axis=1)

    # Fortran repeats the pass in a do while loop, halving mf_min
    # after each failed pass until mf_min <= 0.02. The mf_min values are fixed, so
    # JAX writes the passes out for each value and skips a pass once an earlier
    # pass has filled all holes.
    # Set mf_min, this is tunable. If filling fails on first pass, this is increased
    # and the filling is attempted again, so this is merely a starting value
    mf_min_values = [0.5]
    while mf_min_values[-1] > 0.02:
        # Allowing mf_min = 0.0 is fine in theory, but the way we decompose the stealable mass
        # can result in the sum of that being slightly greater than what we actually have.
        # So if we need to steal all of the mass from a level, it's likely that we're going to
        # steal slightly too much, and drive the above threshold value below. This is a quirk of
        # floating point calculations, and is quite unfortunate, since it forces us to rely on
        # a global fill rather than just setting it to 0. Limiting this to 0.01 is arbirary,
        # you can probably set it lower and be fine
        mf_min_values.append(max(0.01, mf_min_values[-1] / 2.0))

    def fill_holes_parallel_column(field, n_holes, stealable_mass, rho_ds_dz, invrs_rho_ds_dz):
        # The nested loop bodies get field and n_holes from their carry, which holds
        # the values updated by earlier holes.

        def expand_range(_, carry):
            k, k_start, k_end, n_steal_points, stealable_mass_total, hole_mass = carry

            # If we didn't find enough mass to fill the hole (and our bounds aren't exhausted),
            # then expand the range until we have found enough
            l_expand = (stealable_mass_total < jnp.abs(hole_mass)) & (
                (k_start != lower_hf_level) | (k_end != upper_hf_level)
            )

            # Double the number of levels to try stealing from
            n_steal_points_new = 2 * n_steal_points

            # Recalculate bounds of fill range
            k_start_new = jnp.maximum(lower_hf_level, k - n_steal_points_new)
            k_end_new = jnp.minimum(upper_hf_level, k + n_steal_points_new)

            # Sum up new mass in [k_start_new:k_start-1] and [k_end+1:k_end_new]
            l_added = ((k_levels >= k_start_new) & (k_levels < k_start)) | (
                (k_levels > k_end) & (k_levels <= k_end_new)
            )
            stealable_mass_total_new = stealable_mass_total + jnp.sum(
                jnp.where(l_added, stealable_mass, zero)
            )

            # Set new fill range
            return (
                k,
                jnp.where(l_expand, k_start_new, k_start),
                jnp.where(l_expand, k_end_new, k_end),
                jnp.where(l_expand, n_steal_points_new, n_steal_points),
                jnp.where(l_expand, stealable_mass_total_new, stealable_mass_total),
                hole_mass,
            )

        def fill_hole(k_offset, carry):
            field, n_holes, l_field_below_threshold = carry

            k = lower_hf_level + k_offset

            # Only do something if the field value were considering is below threshold
            l_hole = field[k] < threshold

            # Hole detected, calculate the mass needed to fill it
            hole_mass = (field[k] - threshold) * rho_ds_dz[k]

            # We need to steal from other levels, but how far out to reach initially?
            # 2*n_holes seems like a good guess, more holes means more "reach" needed, and
            # this scales with holes is the important part
            n_steal_points = n_holes * 2

            # Caclulate our fill range
            k_start = jnp.maximum(lower_hf_level, k - n_steal_points)
            k_end = jnp.minimum(upper_hf_level, k + n_steal_points)

            # Caclulate the total amount of stealable mass in that range
            l_in_range = (k_levels >= k_start) & (k_levels <= k_end)
            stealable_mass_total = jnp.sum(jnp.where(l_in_range, stealable_mass, zero))

            _, k_start, k_end, _, stealable_mass_total, _ = jax.lax.fori_loop(
                0,
                num_expansions,
                expand_range,
                (k, k_start, k_end, n_steal_points, stealable_mass_total, hole_mass),
            )

            # We can only fill if stealable_mass_total > 0, this is mostly to prevent divide by 0
            l_move_mass = l_hole & (stealable_mass_total > zero)

            # Calculate mass fraction
            # The division is evaluated even when no mass is moved, so guard it.
            mass_fraction_above = one - jnp.maximum(
                (stealable_mass_total + hole_mass)
                / jnp.where(l_move_mass, stealable_mass_total, one), zero
            )

            # Update field values for levels we're stealing from in the fill range
            l_in_range = (k_levels >= k_start) & (k_levels <= k_end)
            field_adj = mass_fraction_above * stealable_mass * invrs_rho_ds_dz
            field = jnp.where(l_move_mass & l_in_range & (field > threshold), field - field_adj, field)

            # Update the field value for the hole were considering
            l_not_filled = stealable_mass_total < jnp.abs(hole_mass)

            # Hole cannot be completely filled, give it all the mass we have available
            # and set the flag indicated there is still a hole
            l_field_below_threshold = l_field_below_threshold | (l_hole & l_not_filled)

            # Hole can be completely filled, which means it should be set to threshold,
            # also decrement the hole count since now there is one less
            n_holes = jnp.where(l_hole & ~l_not_filled, n_holes - 1, n_holes)
            field_k = jnp.where(
                l_not_filled,
                field[k] + stealable_mass_total * invrs_rho_ds_dz[k],
                threshold,
            )
            field = field.at[k].set(jnp.where(l_hole, field_k, field[k]))

            return field, n_holes, l_field_below_threshold

        l_field_below_threshold = jnp.asarray(False)
        field, n_holes, l_field_below_threshold = jax.lax.fori_loop(
            0,
            num_levels,
            fill_hole,
            (field, n_holes, l_field_below_threshold),
        )

        return field, n_holes, l_field_below_threshold

    l_field_below_threshold = jnp.asarray(True)

    # Loop until all holes are filled
    for mf_min in mf_min_values:

        # Calculate the steablable mass. We take stealable_scale = mf_min / n_holes(i) to
        # ensure that each hole doesn't steal too much mass from above threshold levels.
        # In theory, each hole could take up to 1 / n_holes mass, but we have mf_min as a smoothing
        # term, otherwise mass is stolen quite agressively with this approach, since holes
        # have no knowledge of how far out other holes have to search
        # A column whose holes are all filled has n_holes(i) = 0 here. It has nothing
        # to fill, so give it no stealable mass instead of dividing by zero.
        stealable_scale = jnp.where(
            n_holes > 0, (one - mf_min) / jnp.maximum(n_holes, 1), zero
        )[:, None]
        # Fortran uses max( field - threshold, 0 ); the selection gives the same value
        # and, at field == threshold, the derivative of the field > threshold branch
        # that decides below whether mass is stolen.
        stealable_mass = jnp.where(
            l_in_hf_range & (field > threshold), field - threshold, zero
        ) * rho_ds_dz * stealable_scale

        # Fill the holes
        field_pass, n_holes_pass, l_field_below_threshold_col = jax.vmap(fill_holes_parallel_column)(
            field, n_holes, stealable_mass, rho_ds_dz, invrs_rho_ds_dz
        )

        # Skip this pass if an earlier pass filled all holes
        field = jnp.where(l_field_below_threshold, field_pass, field)
        n_holes = jnp.where(l_field_below_threshold, n_holes_pass, n_holes)
        l_field_below_threshold = l_field_below_threshold & jnp.any(l_field_below_threshold_col)

    # If we still have below threshold values, just call the global fill at this point
    field = jax.lax.cond(
        l_field_below_threshold,
        lambda f: fill_holes_global(
            nz,
            ngrdcol,
            threshold,
            lower_hf_level,
            upper_hf_level,
            dz,
            rho_ds,
            f,
        ),
        lambda f: f,
        field,
    )

    return field


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
    elif fill_holes_type == widening_windows:
        # This performs a multi-window technique, splitting the vertical levels into a
        # series of ranges that can be filled in parallel. This can often fail to fill
        # and continually widens the range as neccesary, and if the window range
        # gets wide enough, it is equivalent to a global fill
        filled = fill_holes_widening_windows(
            nz,
            ngrdcol,
            threshold,
            lower_hf_level,
            upper_hf_level,
            dz,
            rho_ds,
            field,
        )
    elif fill_holes_type == smart_window:
        # This performs filling by first scanning the vertical levels to determine a
        # range to fill before actually doing the fill calculations. The window can
        # increase to the whole range, which is effectively a global fallback
        filled = fill_holes_smart_window(
            nz,
            ngrdcol,
            threshold,
            lower_hf_level,
            upper_hf_level,
            dz,
            rho_ds,
            field,
        )
    elif fill_holes_type == smart_window_smooth:
        # This is a modification of the smart_window technique, which includes two "smoothing"
        # features to preserve the field derivative a bit. The efficacy of this over the basic
        # smart_window technique is unclear, and should be investigated more.
        # WARNING: This also potentially fail to fill even if the field average is above
        # threshold (no global fallack), see the subroutine description for more detail
        filled = fill_holes_smart_window_smooth(
            nz,
            ngrdcol,
            threshold,
            lower_hf_level,
            upper_hf_level,
            dz,
            rho_ds,
            field,
        )
    elif fill_holes_type == parallel_fill:
        # This fills holes by considering each hole independently, Making this method the most
        # parallelizable one, but also the most traditionally computationally expensive.
        # This will fall back to a global fill if the locality preserving passes fail
        filled = fill_holes_parallel(
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
        print("ERROR: unknown fill_holes_type: ", fill_holes_type, file=sys.stderr)
        raise ValueError(f"unknown fill_holes_type: {fill_holes_type}")

    field = jnp.where(jnp.any(field < threshold), filled, field)
    return field


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

    Returns:
        (wp2_out, up2_out, vp2_out)

    Fortran declarations (fill_holes_wp2_from_horz_tke in src/CLUBB_core/fill_holes.F90):

        ! --------------------- Input variables ---------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        integer, intent(in) :: &
          lower_hf_level, & ! Lower grid level of global hole-filling range      []
          upper_hf_level    ! Upper grid level of global hole-filling range      []

        real( kind = core_rknd ), intent(in) :: &
          threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                     ! fall                           [Units vary; same as field]

        ! --------------------- Input/Output variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          wp2, &    ! The field (here wp2) that contains holes [m^2/s^2]
          up2, vp2  ! The fields (here up2, vp2) that will be used to fill holes [m^2/s^2]
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
    no_up2_avail = jnp.abs(up2_avail) < _CORE_EPS * 1000.0
    no_vp2_avail = jnp.abs(vp2_avail) < _CORE_EPS * 1000.0
    # Calculate portion of up2/vp2 that we want to take away
    has_donor_energy = up2_vp2_avail > zero
    ratio = jnp.where(
        has_donor_energy,
        missing_wp2 / jnp.where(has_donor_energy, up2_vp2_avail, one),
        zero,
    )

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
    if clubb_at_least_debug_level(3):
        # Host diagnostics only: retain source column/level order and old values.
        def print_tke_errors(w_old, u_old, v_old, w_new, u_new, v_new, limit):
            import numpy as np
            for i in range(w_old.shape[0]):
                for k in range(lower_hf_level, upper_hf_level + 1):
                    difference = abs(u_old[i,k] + v_old[i,k] + w_old[i,k] - u_new[i,k] - v_new[i,k] - w_new[i,k])
                    if difference > np.finfo(w_old.dtype).eps * 1000:
                        print(
                            "Warning! fill_holes_wp2_from_horz_tke is not conservative at level ",
                            i+1,
                            k+1,
                            ". Difference: ",
                            difference,
                            file=sys.stderr,
                        )
                    if w_new[i,k] > limit:
                        print(
                            "Warning in fill_holes_wp2_from_horz_tke! Too much TKE added to wp2. old: ",
                            w_old[i,k],
                            ", new: ",
                            w_new[i,k],
                            file=sys.stderr,
                        )
                    if w_new[i,k] < w_old[i,k]:
                        print(
                            "Warning in fill_holes_wp2_from_horz_tke! wp2 was reduced. old: ",
                            w_old[i,k],
                            ", new: ",
                            w_new[i,k],
                            file=sys.stderr,
                        )
                    for name, old, new in (("up2", u_old, u_new), ("vp2", v_old, v_new)):
                        if new[i,k] > old[i,k]:
                            print(
                                "Warning in fill_holes_wp2_from_horz_tke! " + name + " increased. old: ",
                                old[i,k],
                                ", new: ",
                                new[i,k],
                                file=sys.stderr,
                            )
                        if old[i,k] <= limit and new[i,k] < old[i,k]:
                            print(
                                "Warning in fill_holes_wp2_from_horz_tke! TKE was falsely taken from " + name + ". Threshold: ",
                                limit,
                                ", old: ",
                                old[i,k],
                                ", new: ",
                                new[i,k],
                                file=sys.stderr,
                            )
        jax.lax.cond(
            jnp.any(wp2 < threshold),
            lambda _: jax.debug.callback(print_tke_errors, wp2, up2, vp2, wp2_out, up2_out, vp2_out, threshold, ordered=True),
            lambda _: None, operand=None,
        )
    return wp2_out, up2_out, vp2_out


# Microphysics callers use batched columns; source vertical/species loops below
# are expressed as array operations and static species loops.
def hole_filling_hm_one_lev(num_hm_fill, hm_one_lev):
    """Fills holes between same-phase (i.e. either liquid or frozen) hydrometeors for
    one height level.

    Warning: Do not input hydrometeors of different phases, e.g. liquid and frozen.
    Otherwise heat will not be conserved.

    Returns:
        result

    Fortran declarations (hole_filling_hm_one_lev in src/CLUBB_core/fill_holes.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: num_hm_fill ! number of hydrometeors involved

        real(kind = core_rknd), dimension(num_hm_fill), intent(in) :: hm_one_lev

        !--------------------------- Output Variables ---------------------------
        real(kind = core_rknd), dimension(num_hm_fill), intent(out) :: hm_one_lev_filled
    """

    from clubb_jax.src.CLUBB_core.constants_clubb import eps

    # Total negative mass and available positive mass of the same phase.
    # Do not mix liquid and frozen species here: heat would not be conserved.
    total_hole = jnp.sum(jnp.minimum(hm_one_lev, 0.0), axis=-1, keepdims=True)
    total_mass = jnp.sum(jnp.maximum(hm_one_lev, 0.0), axis=-1, keepdims=True)

    # If positive mass is insufficient, partially fill the holes with all of it;
    # otherwise eliminate holes and reduce the positive species proportionally.
    hm_one_lev_filled = jnp.where(
        jnp.abs(total_hole) > total_mass,
        jnp.minimum(hm_one_lev, 0.0)
        * (1.0 + total_mass / jnp.where(total_hole != 0, total_hole, 1.0)),
        jnp.maximum(hm_one_lev, 0.0)
        * (1.0 + total_hole / jnp.where(total_mass != 0, total_mass, 1.0)),
    )

    # Leave the input unchanged if no material is available.
    result = jnp.where(jnp.abs(total_mass) < eps, hm_one_lev, hm_one_lev_filled)
    if clubb_at_least_debug_level(2):
        # Host adaptation of the source scalar-level/species diagnostic loops.
        def print_hm_errors(original, filled, hole, mass):
            import numpy as np
            for index in np.ndindex(original.shape[:-1]):
                if abs(mass[index][0]) < eps:
                    print("Warning: One-level hole filling was not successful! total_mass ~= 0")
                    continue
                if abs(hole[index][0]) > mass[index][0]:
                    for _ in range(num_hm_fill):
                        print(
                            'Warning: One-level hole filling was not able to fill holes completely! The holes were'
                            ' filled partially. |total_hole| > total_mass',
                        )
                before, after = np.sum(original[index]), np.sum(filled[index])
                if abs(before - after) > abs(before + after) * eps / 2:
                    print("Warning: Hole filling was not conservative!")
                if np.any(filled[index] < zero):
                    print("Warning: Hole filling failed! A hole could not be filled.")
        jax.debug.callback(print_hm_errors, hm_one_lev, result, total_hole, total_mass, ordered=True)
    return result


def fill_holes_hydromet_api(
    ngrdcol, nzt, hydromet_dim, hydromet,  # In
    l_frozen_hm, l_mix_rat_hm,    # In
):
    """Fills holes independently in each column between same-phase hydrometeors(i.e. for frozen hydrometeors).
    The hole filling conserves water substance between all same-phase (frozen or liquid)
    hydrometeors at each height level.

    Attention: The hole filling for the liquid phase hydrometeors is not yet implemented

    Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!

    Fortran declarations (fill_holes_hydromet_api in src/CLUBB_core/fill_holes.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ngrdcol, nzt, hydromet_dim ! Columns, levels and species [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hydromet ! Hydrometeor fields for all columns [units vary]

        logical, dimension(hydromet_dim), intent(in) :: &
          l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
          l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
          hydromet_filled ! Filled hydrometeor fields [units vary]
    """

    # Restrict same-level redistribution to frozen mixing ratios. Liquid
    # rain and all number concentrations bypass this phase-conserving step.
    frozen = l_frozen_hm & l_mix_rat_hm
    hydromet_frozen = jnp.where(frozen, hydromet, 0.0)
    hydromet_frozen_filled = hole_filling_hm_one_lev(hydromet_dim, hydromet_frozen)
    return jnp.where(frozen, hydromet_frozen_filled, hydromet)


def fill_holes_wv(
    ngrdcol, nzt, dt, exner, hydromet_name,  # In
    rvm_mc, thlm_mc, hydromet,      # InOut
):
    """Fills remaining hydrometeor holes from vapor in each column and adjusts
    the vapor and liquid-water potential-temperature tendencies.

    Returns:
        (rvm_mc, thlm_mc, hydromet)

    Fortran declarations (fill_holes_wv in src/CLUBB_core/fill_holes.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ngrdcol, nzt ! Number of columns and levels [-]

        real( kind = core_rknd ), intent(in) ::  &
          dt           ! Timestep         [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          exner        ! Exner function                            [-]

        character(len=10), intent(in) :: hydromet_name ! Species name [-]

        !--------------------------- Input/Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
          rvm_mc,  & ! Vapor tendency [kg/kg/s]
          thlm_mc, & ! Liquid-water potential-temperature tendency [K/s]
          hydromet   ! Hydrometeor mixing ratio [kg/kg]
    """

    from clubb_jax.src.CLUBB_core.constants_clubb import zero_threshold, Lv, Ls, Cp

    # Transfer the remaining negative mixing ratio to water vapor and
    # compensate latent heat using Lv for rain and Ls for frozen species.
    rvm_clip_tndcy = jnp.where(hydromet < zero_threshold, hydromet / dt, 0.0)
    rvm_mc = rvm_mc + rvm_clip_tndcy
    if hydromet_name == "rrm":
        thlm_mc = thlm_mc - rvm_clip_tndcy * (Lv / (Cp * exner))
    elif hydromet_name in ("rim", "rsm", "rgm"):
        thlm_mc = thlm_mc - rvm_clip_tndcy * (Ls / (Cp * exner))
    else:
        raise ValueError("Fatal error in microphys_driver: unknown hydrometeor")
    hydromet = jnp.maximum(hydromet, zero_threshold)
    return rvm_mc, thlm_mc, hydromet


def fill_holes_driver_api(
    gr, ngrdcol, nzt, dt, hydromet_dim,  # In
    hm_metadata, l_fill_holes_hm,        # In
    rho_ds_zt, exner, fill_holes_type,   # In
    stats, thlm_mc, rvm_mc, hydromet,    # InOut
):
    """Apply hydrometeor hole filling and clipping to all model columns.
    Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
    The hole filling conserves water substance between all same-phase (frozen or liquid)
    hydrometeors at each height level.

    Attention: The hole filling for the liquid phase hydrometeors is not yet implemented

    Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!

    Returns:
        (stats, thlm_mc, rvm_mc, hydromet)

    Fortran declarations (fill_holes_driver_api in src/CLUBB_core/fill_holes.F90):

        !----------------------- Input Variables -----------------------
        type(grid), intent(in) :: &
          gr

        integer, intent(in) :: &
          ngrdcol, & ! Number of model columns
          nzt, &     ! Number of vertical levels
          hydromet_dim, & ! Number of hydrometeor fields
          fill_holes_type

        real(kind=core_rknd), intent(in) :: &
          dt                   ! Timestep         [s]

        type(hm_metadata_type), intent(in) :: &
          hm_metadata

        logical, intent(in) :: &
          l_fill_holes_hm

        real(kind=core_rknd), dimension(ngrdcol,nzt), intent(in) :: &
          rho_ds_zt, & ! Dry, static density on thermo. levels    [kg/m^3]
          exner        ! Exner function                           [-]

        !----------------------- Input/Output Variables -----------------------
        type(stats_type), intent(inout) :: &
          stats

        real(kind=core_rknd), dimension(ngrdcol,nzt), intent(inout) :: &
          thlm_mc, &           ! Microphysics contributions to liquid potential temp [K/s]
          rvm_mc               ! Microphysics contributions to vapor water           [kg/kg/s]

        real(kind=core_rknd), dimension(ngrdcol,nzt,hydromet_dim), intent(inout) :: &
          hydromet    ! Mean of hydrometeor fields    [units vary]
    """

    from clubb_jax.src.CLUBB_core.constants_clubb import zero_threshold, Lv, Ls, Cp

    # Start statistics for same-phase and vertical hole filling.
    for i in range(hydromet_dim):
        _, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(
            i, hydromet_dim, hm_metadata.hydromet_list
        )
        stats = stats.begin_budget(name_hf, hydromet[..., i] / dt)
    if l_fill_holes_hm:
        hydromet = fill_holes_hydromet_api(
            ngrdcol, nzt, hydromet_dim, hydromet,                        # In
            hm_metadata.l_frozen_hm, hm_metadata.l_mix_rat_hm,  # In
        )
    for i in range(hydromet_dim):
        _, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(
            i, hydromet_dim, hm_metadata.hydromet_list
        )
        if clubb_at_least_debug_level(1):
            def print_hm_holes(values, name=hm_metadata.hydromet_list[i]):
                import numpy as np
                for k, col in np.argwhere((values < zero_threshold).T):
                    print(
                        name + " < ",
                        zero_threshold,
                        " in fill_holes_driver_api at i, k= ",
                        col+1,
                        k+1,
                        file=sys.stderr,
                    )
            jax.lax.cond(
                jnp.any(hydromet[..., i] < zero_threshold),
                lambda _: jax.debug.callback(print_hm_holes, hydromet[..., i], ordered=True),
                lambda _: None, operand=None,
            )
        hydromet_name = hm_metadata.hydromet_list[i]
        if hydromet_name.startswith("r"):
            hydromet_filled = fill_holes_vertical(
                nzt, ngrdcol, zero_threshold, 0, nzt - 1,
                gr.dzt, rho_ds_zt, 1, fill_holes_type, hydromet[..., i],
            )
            hydromet = hydromet.at[..., i].set(hydromet_filled)
        stats = stats.finalize_budget(name_hf, hydromet[..., i] / dt)
        stats = stats.begin_budget(name_wvhf, hydromet[..., i] / dt)
        if hydromet_name.startswith("r"):
            rvm_mc, thlm_mc, hydromet_filled = fill_holes_wv(
                ngrdcol, nzt, dt, exner, hydromet_name,      # In
                rvm_mc, thlm_mc, hydromet[..., i],  # InOut
            )
            hydromet = hydromet.at[..., i].set(hydromet_filled)
        stats = stats.finalize_budget(name_wvhf, hydromet[..., i] / dt)

        # Remove sub-tolerance positive mixing ratios with the same vapor and
        # latent-heat compensation, recording the clipping budget separately.
        if hydromet_name.startswith("r"):
            stats = stats.begin_budget(name_cl, hydromet[..., i] / dt)
            hydromet = hydromet.at[..., i].set(jnp.maximum(hydromet[..., i], zero_threshold))
            small = hydromet[..., i] <= hm_metadata.hydromet_tol[i]
            rvm_mc = rvm_mc + jnp.where(small, hydromet[..., i] / dt, 0.0)
            latent_heat = Lv if hydromet_name == "rrm" else Ls
            thlm_mc = thlm_mc - jnp.where(
                small, (latent_heat / (Cp * exner)) * (hydromet[..., i] / dt), 0.0
            )
            hydromet = hydromet.at[..., i].set(jnp.where(small, 0.0, hydromet[..., i]))
            stats = stats.finalize_budget(name_cl, hydromet[..., i] / dt)

    # Enforce maximum particle radii through concentration minima. Mass fields
    # remain unchanged; only number concentrations contribute to these budgets.
    hydromet_clipped = clip_hydromet_conc_mvr(ngrdcol, nzt, hydromet_dim, hm_metadata, hydromet)
    for i in range(hydromet_dim):
        if hm_metadata.hydromet_list[i].startswith("N"):
            _, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(
                i, hydromet_dim, hm_metadata.hydromet_list
            )
            stats = stats.begin_budget(name_cl, hydromet[..., i] / dt)
            hydromet = hydromet.at[..., i].set(hydromet_clipped[..., i])
            stats = stats.finalize_budget(name_cl, hydromet[..., i] / dt)
    return stats, thlm_mc, rvm_mc, hydromet


def clip_hydromet_conc_mvr(
    ngrdcol, nzt, hydromet_dim, hm_metadata,  # In
    hydromet,                        # In
):
    """Increases the value of a hydrometeor concentration when it is too small,
    according to the hydrometeor mixing ratio (which remains unchanged) and
    the maximum acceptable drop or particle mean volume radius for that
    hydrometeor species.

    Returns:
        hydromet_clipped

    Fortran declarations (clip_hydromet_conc_mvr in src/CLUBB_core/fill_holes.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol, & ! Number of columns [-]
          nzt        ! Number of levels [-]

        integer, intent(in) :: &
          hydromet_dim    ! Number of hydrometeor fields

        type (hm_metadata_type), intent(in) :: &
              hm_metadata

        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
          hydromet    ! Mean of hydrometeor fields    [units vary]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
          hydromet_clipped    ! Clipped mean of hydrometeor fields    [units vary]
    """
    from clubb_jax.src.CLUBB_core.constants_clubb import pi, rho_lw, rho_ice
    from clubb_jax.src.CLUBB_core.index_mapping import Nx2rx_hm_idx, mvr_hm_max


    # Preserve all mixing ratios; loop only over number concentrations.
    hydromet_clipped = hydromet
    for idx in range(hydromet_dim):
        if hm_metadata.hydromet_list[idx].startswith("N"):
            # N_min = r / ((4/3)*pi*rho_particle*r_max^3), using liquid
            # density for rain and ice density for frozen species.
            density = rho_lw if hm_metadata.hydromet_list[idx] == "Nrm" else rho_ice
            Nxm_min_coef = 1.0 / ((4.0 / 3.0) * pi * density * mvr_hm_max(idx, hm_metadata) ** 3)

            # Raise number where mass is present, otherwise set it to zero.
            rx = hydromet[..., Nx2rx_hm_idx(idx, hm_metadata)]
            hydromet_clipped = hydromet_clipped.at[..., idx].set(
                jnp.where(rx > 0.0, jnp.maximum(hydromet[..., idx], Nxm_min_coef * rx), 0.0)
            )
    return hydromet_clipped


def setup_stats_names(ihm, hydromet_dim, hydromet_list):
    """Determines the stats output names depending on the hydrometeor.
    Attention: hydromet_list needs to be set up before this routine is called.

    Bogus example

    Fortran declarations (setup_stats_names in src/CLUBB_core/fill_holes.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ihm, hydromet_dim

        character(len=10), dimension(hydromet_dim), intent(in) :: &
          hydromet_list

        !--------------------------- Input/Output Variables ---------------------------
        real( kind = core_rknd ), intent(inout) :: &
          max_velocity ! Maximum sedimentation velocity [m/s]

        !--------------------------- Optional Output Variables ---------------------------
        character(len=*), intent(out), optional :: &
          name_bt, name_hf, name_wvhf, name_cl, name_mc
    """

    # Negative speeds denote downward sedimentation. Ice uses the source
    # 1.2 m/s limit, snow 2 m/s, and graupel 20 m/s. Cloud number shares the
    # rain limit because Morrison has no separate cloud-water speed limit.
    name = hydromet_list[ihm]
    max_velocity = {
        "rrm": -9.1,
        "Nrm": -9.1,
        "rim": -1.2,
        "Nim": -1.2,
        "rsm": -2.0,
        "Nsm": -2.0,
        "rgm": -20.0,
        "Ngm": -20.0,
        "Ncm": -9.1,
    }.get(name, 0.0)
    if not max_velocity:
        return 0.0, "", "", "", "", ""
    return (
        max_velocity,
        name + "_bt",
        name + "_hf" if name.startswith("r") else "",
        name + "_wvhf" if name.startswith("r") else "",
        name + "_cl",
        name + "_mc",
    )

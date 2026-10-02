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
    return wp2_out, up2_out, vp2_out


__all__ = [
    "fill_holes_global",
    "fill_holes_sliding_window",
    "fill_holes_vertical",
    "fill_holes_wp2_from_horz_tke",
]


# Microphysics callers use batched columns; source vertical/species loops below
# are expressed as array operations and static species loops.
def hole_filling_hm_one_lev(num_hm_fill, hm_one_lev):
    from clubb_jax.src.CLUBB_core.constants_clubb import eps
    total_hole = jnp.sum(jnp.minimum(hm_one_lev, 0.0), axis=-1, keepdims=True)
    total_mass = jnp.sum(jnp.maximum(hm_one_lev, 0.0), axis=-1, keepdims=True)
    hm_one_lev_filled = jnp.where(
        jnp.abs(total_hole) > total_mass,
        jnp.minimum(hm_one_lev, 0.0) * (1.0 + total_mass / jnp.where(total_hole != 0, total_hole, 1.0)),
        jnp.maximum(hm_one_lev, 0.0) * (1.0 + total_hole / jnp.where(total_mass != 0, total_mass, 1.0)),
    )
    return jnp.where(jnp.abs(total_mass) < eps, hm_one_lev, hm_one_lev_filled)


def fill_holes_hydromet_api(nzt, hydromet_dim, hydromet, l_frozen_hm, l_mix_rat_hm):
    frozen = jnp.asarray(l_frozen_hm) & jnp.asarray(l_mix_rat_hm)
    hydromet_frozen = jnp.where(frozen, hydromet, 0.0)
    hydromet_frozen_filled = hole_filling_hm_one_lev(hydromet_dim, hydromet_frozen)
    return jnp.where(frozen, hydromet_frozen_filled, hydromet)


def fill_holes_wv(nzt, dt, exner, hydromet_name, rvm_mc, thlm_mc, hydromet):
    from clubb_jax.src.CLUBB_core.constants_clubb import zero_threshold, Lv, Ls, Cp
    rvm_clip_tndcy = jnp.where(hydromet < zero_threshold, hydromet / dt, 0.0)
    rvm_mc = rvm_mc + rvm_clip_tndcy
    if hydromet_name == 'rrm':
        thlm_mc = thlm_mc - rvm_clip_tndcy * (Lv / (Cp * exner))
    elif hydromet_name in ('rim', 'rsm', 'rgm'):
        thlm_mc = thlm_mc - rvm_clip_tndcy * (Ls / (Cp * exner))
    else:
        raise ValueError('Fatal error in microphys_driver: unknown hydrometeor')
    hydromet = jnp.maximum(hydromet, zero_threshold)
    return rvm_mc, thlm_mc, hydromet


def fill_holes_driver_api(gr, ngrdcol, nzt, dt, hydromet_dim, hm_metadata, l_fill_holes_hm,
                         rho_ds_zt, exner, fill_holes_type, stats,
                         thlm_mc, rvm_mc, hydromet):
    from clubb_jax.src.CLUBB_core.constants_clubb import zero_threshold, Lv, Ls, Cp
    # Start statistics for same-phase and vertical hole filling.
    for i in range(hydromet_dim):
        _, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(i, hydromet_dim, hm_metadata.hydromet_list)
        stats = stats.begin_budget(name_hf, hydromet[..., i] / dt)
    if l_fill_holes_hm:
        hydromet = fill_holes_hydromet_api(nzt, hydromet_dim, hydromet, hm_metadata.l_frozen_hm, hm_metadata.l_mix_rat_hm)
    for i in range(hydromet_dim):
        _, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(i, hydromet_dim, hm_metadata.hydromet_list)
        hydromet_name = hm_metadata.hydromet_list[i]
        if hydromet_name.startswith('r'):
            # Source calls each column with ngrdcol=1 but passes the full gr%dzt.
            # Its explicit-shape dz(1,nzt) dummy uses the first nzt elements in
            # Fortran storage order. Preserve that deliberately retained scalar
            # grid-metric mapping, including on a multicolumn grid.
            dz_scalar = gr.dzt.T.reshape(-1)[:nzt][None, :]
            hydromet_filled = jax.vmap(lambda rho_col, hm_col: fill_holes_vertical(
                nzt, 1, zero_threshold, 0, nzt - 1,
                dz_scalar, rho_col[None, :], 1, fill_holes_type,
                hm_col[None, :])[0])(rho_ds_zt, hydromet[..., i])
            hydromet = hydromet.at[..., i].set(hydromet_filled)
        stats = stats.finalize_budget(name_hf, hydromet[..., i] / dt)
        stats = stats.begin_budget(name_wvhf, hydromet[..., i] / dt)
        if hydromet_name.startswith('r'):
            rvm_mc, thlm_mc, hydromet_filled = fill_holes_wv(nzt, dt, exner, hydromet_name, rvm_mc, thlm_mc, hydromet[..., i])
            hydromet = hydromet.at[..., i].set(hydromet_filled)
        stats = stats.finalize_budget(name_wvhf, hydromet[..., i] / dt)
        if hydromet_name.startswith('r'):
            stats = stats.begin_budget(name_cl, hydromet[..., i] / dt)
            hydromet = hydromet.at[..., i].set(jnp.maximum(hydromet[..., i], zero_threshold))
            small = hydromet[..., i] <= hm_metadata.hydromet_tol[i]
            rvm_mc = rvm_mc + jnp.where(small, hydromet[..., i] / dt, 0.0)
            latent_heat = Lv if hydromet_name == 'rrm' else Ls
            thlm_mc = thlm_mc - jnp.where(small, (latent_heat / (Cp * exner)) * (hydromet[..., i] / dt), 0.0)
            hydromet = hydromet.at[..., i].set(jnp.where(small, 0.0, hydromet[..., i]))
            stats = stats.finalize_budget(name_cl, hydromet[..., i] / dt)
    hydromet_clipped = clip_hydromet_conc_mvr(nzt, hydromet_dim, hm_metadata, hydromet)
    for i in range(hydromet_dim):
        if hm_metadata.hydromet_list[i].startswith('N'):
            _, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(i, hydromet_dim, hm_metadata.hydromet_list)
            stats = stats.begin_budget(name_cl, hydromet[..., i] / dt)
            hydromet = hydromet.at[..., i].set(hydromet_clipped[..., i])
            stats = stats.finalize_budget(name_cl, hydromet[..., i] / dt)
    return stats, thlm_mc, rvm_mc, hydromet


def clip_hydromet_conc_mvr(nzt, hydromet_dim, hm_metadata, hydromet):
    from clubb_jax.src.CLUBB_core.constants_clubb import pi, rho_lw, rho_ice
    from clubb_jax.src.CLUBB_core.index_mapping import Nx2rx_hm_idx, mvr_hm_max
    hydromet_clipped = hydromet
    for idx in range(hydromet_dim):
        if hm_metadata.hydromet_list[idx].startswith('N'):
            density = rho_lw if hm_metadata.hydromet_list[idx] == 'Nrm' else rho_ice
            Nxm_min_coef = 1.0 / ((4.0 / 3.0) * pi * density * mvr_hm_max(idx, hm_metadata) ** 3)
            rx = hydromet[..., Nx2rx_hm_idx(idx, hm_metadata)]
            hydromet_clipped = hydromet_clipped.at[..., idx].set(jnp.where(rx > 0.0, jnp.maximum(hydromet[..., idx], Nxm_min_coef * rx), 0.0))
    return hydromet_clipped


def setup_stats_names(ihm, hydromet_dim, hydromet_list):
    name = hydromet_list[ihm]
    max_velocity = {'rrm': -9.1, 'Nrm': -9.1, 'rim': -1.2, 'Nim': -1.2,
                    'rsm': -2.0, 'Nsm': -2.0, 'rgm': -20.0, 'Ngm': -20.0, 'Ncm': -9.1}.get(name, 0.0)
    if not max_velocity:
        return 0.0, '', '', '', '', ''
    return (max_velocity, name + '_bt', name + '_hf' if name.startswith('r') else '',
            name + '_wvhf' if name.startswith('r') else '', name + '_cl', name + '_mc')

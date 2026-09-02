"""JAX port of ``Radiation/radiation_module.F90``.

Description:
    Advance the active radiation scheme when needed and update
    radiation statistics every sampling timestep.

JAX adaptations:
    ``radiation_variables_module`` output arrays and ``parameters_radiation``
    module variables are explicit functional inputs and returns. Their order
    follows the corresponding Fortran ``intent(inout)`` and ``intent(out)``
    values. BUGSrad and SILHS radiation are rejected by case initialization.
"""

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import Cp
from clubb_jax.src.CLUBB_core.err_info_codes import (
    ERR_RADIATION_COS_SOLAR_ZEN_INVALID_HOUR,
    ERR_RADIATION_COS_SOLAR_ZEN_TIME_NOT_FOUND,
)
from clubb_jax.src.CLUBB_core.grid_class import ddzm
from clubb_jax.src.Radiation.cos_solar_zen_module import cos_solar_zen
from clubb_jax.src.Radiation.rad_lwsw_module import sunray_sw
from clubb_jax.src.Radiation.simple_rad_module import (
    simple_rad,
    simple_rad_bomex,
    simple_rad_lba,
)
from clubb_jax.src.Radiation.soil_vegetation import advance_soil_veg

configure_jax_precision()


@partial(
    jax.jit,
    static_argnames=("ngrdcol", "hydromet_dim", "pdf_dim", "lh_num_samples", "day", "month", "year"),
)
def advance_clubb_radiation(
    gr, ngrdcol, hydromet_dim, pdf_dim, lh_num_samples,
    l_rad_itime, dt_main, day, month, year,
    lat_vals, lon_vals, time_current, time_initial,
    rho, rho_zm, p_in_Pa, exner,
    wpthlp_sfc, wprtp_sfc, p_sfc,
    cloud_frac, ice_supersat_frac,
    thlm, rtm, rcm,
    X_nl_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,
    lh_sample_point_weights, hydromet, hm_metadata,
    stats, err_info,
    deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
    # JAX adaptation: explicit ``radiation_variables_module`` state.
    radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down,
    radht_SW, radht_LW, Frad_SW, Frad_LW,
    # JAX adaptation: immutable ``parameters_radiation`` module state.
    radiation_parameters,
):
    """Advance the active radiation scheme and update radiation statistics.

    Returns the Fortran ``stats`` and ``err_info`` inout values, followed by
    soil/vegetation inout values, ``radht``, and explicit radiation-module
    arrays in the same order as their source ownership.
    """
    # ------------------------ Input Variables ------------------------

    # ------------------------ Input/Output Variables ------------------------

    # ------------------------ Output Variables ------------------------

    # ------------------------ Begin Code -----------------------

    # Soil vegetation is not neccesarily just Radiation, but we put it here to avoid
    # having to call it from the main driver, since it is only used in cases that
    # also use radiation - currently only gabls3. We could move it to the main driver
    # if we wanted to use it in more cases, but for now this is sufficient.
    if radiation_parameters.l_soil_veg:
        # This modifies the soil and vegetation temperatures at the surface, so
        # we need only _sfc arrays or surface slices, e.g (:,1)
        stats, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K = advance_soil_veg(
            ngrdcol, dt_main, rho_zm[:, 0],
            Frad_SW_up[:, 0], Frad_SW_down[:, 0], Frad_LW_down[:, 0],
            wpthlp_sfc, wprtp_sfc, p_sfc, stats,
            deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
        )

    # JAX adaptation: ``lax.cond`` requires branch callables for the source
    # ``if ( l_rad_itime )`` block and its retained module-output values.
    # Only advance radiation if l_rad_itime is true.
    def advance(_):
        # Advance a radiation scheme
        # With this call ordering, snow and ice water mixing ratio will be
        # updated by the microphysics, but thlm and rtm will not.  This
        # somewhat inconsistent, but we would need to move the call to
        # radiation before the call the microphysics to change this.
        # -dschanen 17 Aug 2009
        return radiation_driver(
            gr, time_current, time_initial, hydromet_dim,
            ngrdcol, day, month, year, lat_vals, lon_vals,
            rho, rho_zm, p_in_Pa, exner, cloud_frac, ice_supersat_frac,
            thlm, rtm, rcm, hydromet, hm_metadata, stats, err_info,
            radiation_parameters,
        )

    def retain(_):
        return (
            stats, err_info, radht, Frad, Frad_SW_up, Frad_LW_up,
            Frad_SW_down, Frad_LW_down, radht_SW, radht_LW, Frad_SW, Frad_LW,
        )

    (
        stats, err_info, radht, Frad, Frad_SW_up, Frad_LW_up,
        Frad_SW_down, Frad_LW_down, radht_SW, radht_LW, Frad_SW, Frad_LW,
    ) = jax.lax.cond(l_rad_itime, advance, retain, operand=None)

    # We update stats here each sample timestep - even if radiation is not advanced
    stats = update_radiation_variables(
        ngrdcol, gr.nzm, gr.nzt, radht, Frad, Frad_SW_up, Frad_LW_up,
        Frad_SW_down, Frad_LW_down, radht_SW, radht_LW, Frad_SW, Frad_LW,
        stats, radiation_parameters,
    )

    return (
        stats, err_info, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
        radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down,
        radht_SW, radht_LW, Frad_SW, Frad_LW,
    )


@partial(jax.jit, static_argnames=("hydromet_dim", "ngrdcol", "day", "month", "year"))
def radiation_driver(
    gr, time_current, time_initial, hydromet_dim,
    ngrdcol, day, month, year, lat_vals, lon_vals,
    rho, rho_zm, p_in_Pa, exner, cloud_frac, ice_supersat_frac,
    thlm, rtm, rcm, hydromet, hm_metadata, stats,
    err_info,
    # JAX adaptation: immutable ``parameters_radiation`` module state.
    radiation_parameters,
):
    """Compute a radiation tendency.

    References:
        None

    Returns Fortran inout values first, then ``radht`` and the explicit
    radiation-variable outputs that are module state in the source.
    """
    # ------------------------ Input Variables ------------------------

    # ------------------------ Input/Output Variables ------------------------

    # ------------------------ Output Variables ------------------------

    # ------------------------ Local Variables ------------------------

    # ------------------------ Begin Code ------------------------

    # Initialize all outputs to 0.
    Frad = jnp.zeros((ngrdcol, gr.nzm), dtype=rho.dtype)
    Frad_SW_up = jnp.zeros_like(Frad)
    Frad_LW_up = jnp.zeros_like(Frad)
    Frad_SW_down = jnp.zeros_like(Frad)
    Frad_LW_down = jnp.zeros_like(Frad)
    radht = jnp.zeros((ngrdcol, gr.nzt), dtype=rho.dtype)
    radht_SW = jnp.zeros_like(radht)
    radht_LW = jnp.zeros_like(radht)
    Frad_SW = jnp.zeros_like(Frad)
    Frad_LW = jnp.zeros_like(Frad)

    # Toggle for centered/forward differencing (in sunray_sw interpolations)
    # To use centered differencing, set the toggle to .true.
    # To use forward differencing, set the toggle to  .false.
    l_center = True

    # If l_fix_cos_solar_zen is not set in the model.in, calculate amu0
    # Otherwise, it was defined in cos_solar_zen_list file
    if radiation_parameters.l_sw_radiation:
        if radiation_parameters.l_fix_cos_solar_zen:
            if radiation_parameters.nparam > 1:
                # Find the closest time value greater than or equal to time_current
                amu0_index = jnp.searchsorted(
                    radiation_parameters.cos_solar_zen_times[:radiation_parameters.nparam],
                    time_current,
                    side="left",
                )
                l_time_not_found = amu0_index >= radiation_parameters.nparam
                # JAX indexing must remain in bounds.  Preserve the source fatal
                # path through ``err_info`` and let the host driver raise it.
                err_info = err_info.set_fatal(l_time_not_found)
                err_info = err_info.set_reason(
                    ERR_RADIATION_COS_SOLAR_ZEN_TIME_NOT_FOUND,
                    l_time_not_found,
                )
                amu0_index = jnp.minimum(amu0_index, radiation_parameters.nparam - 1)
            else:
                amu0_index = 0
            amu0 = radiation_parameters.cos_solar_zen_values[amu0_index]
        else:
            amu0 = cos_solar_zen(day, month, year, time_current, lat_vals, lon_vals)
            # The source stops when ``int(hour)`` is outside 0:23.  The JAX
            # routine carries that invalid path as NaN; preserve the source
            # fatal behavior through the functional error state.
            l_invalid_solar_time = ~jnp.isfinite(amu0)
            err_info = err_info.set_fatal(l_invalid_solar_time)
            err_info = err_info.set_reason(
                ERR_RADIATION_COS_SOLAR_ZEN_INVALID_HOUR,
                l_invalid_solar_time,
            )
    else:
        amu0 = jnp.asarray(0.0, dtype=rho.dtype)

    if radiation_parameters.rad_scheme == "simplified":
        # ----------------------------------------------------------------
        # Simplified radiation
        # ----------------------------------------------------------------
        # JAX adaptation: ``lax.cond`` requires callables for the source
        # shortwave ``if ( l_sw_radiation .and. amu0 > 0 )`` block.
        def compute_shortwave(_):
            amu0_core_rknd = amu0.astype(rho.dtype)
            if radiation_parameters.nparam > 1:
                Fs0 = jnp.interp(
                    amu0_core_rknd,
                    radiation_parameters.cos_solar_zen_values[:radiation_parameters.nparam],
                    radiation_parameters.Fs_values[:radiation_parameters.nparam],
                )
            else:
                Fs0 = radiation_parameters.Fs_values[0]
            Frad_SW = sunray_sw(
                ngrdcol, gr.nzt, rcm, rho, amu0_core_rknd, gr.dzt, gr.zm, gr.zt,
                radiation_parameters.eff_drop_radius, radiation_parameters.alvdr,
                radiation_parameters.gc, Fs0, radiation_parameters.omega, l_center,
            )
            radht_SW_ddzm = ddzm(gr.nzm, gr.nzt, ngrdcol, gr, Frad_SW)
            radht_SW = -radht_SW_ddzm / (rho * Cp)
            return Frad_SW, radht_SW

        def no_shortwave(_):
            return Frad_SW, radht_SW

        Frad_SW, radht_SW = jax.lax.cond(amu0 > 0.0, compute_shortwave, no_shortwave, operand=None)
        stats, Frad_LW, radht_LW = simple_rad(
            gr, ngrdcol, rho, rho_zm, rtm, rcm, exner, stats, radiation_parameters
        )
        Frad = Frad_SW + Frad_LW
        radht = radht_SW + radht_LW

    elif radiation_parameters.rad_scheme == "simplified_bomex":
        # ----------------------------------------------------------------
        # GCSS BOMEX specifiction radiation
        # ----------------------------------------------------------------
        radht = simple_rad_bomex(gr, ngrdcol)

    elif radiation_parameters.rad_scheme == "lba":
        radht = simple_rad_lba(gr, ngrdcol, time_current, time_initial, radiation_parameters)

    # TODO(port-mirror): The JAX lifecycle calls this driver for
    # ``rad_scheme='none'`` so sampled radiation diagnostics retain their
    # source cadence. The Fortran lifecycle normally skips this dispatch.
    elif radiation_parameters.rad_scheme != "none":
        # BUGSrad and SILHS radiation are intentionally rejected at initialization.
        raise ValueError(f"Unsupported rad_scheme={radiation_parameters.rad_scheme!r} in JAX radiation_driver.")

    # The source checks ``err_info`` after radiation.  A jitted kernel returns
    # that state; ``_advance_radiation`` performs the equivalent host-side stop.
    return (
        stats, err_info, radht, Frad, Frad_SW_up, Frad_LW_up,
        Frad_SW_down, Frad_LW_down, radht_SW, radht_LW, Frad_SW, Frad_LW,
    )


@partial(jax.jit, static_argnames=("ngrdcol", "nzm", "nzt"))
def update_radiation_variables(
    ngrdcol, nzm, nzt, radht, Frad, Frad_SW_up, Frad_LW_up,
    Frad_SW_down, Frad_LW_down, radht_SW, radht_LW, Frad_SW, Frad_LW,
    # JAX adaptation: source imports these module values rather than passing
    # them; functional statistics updates require them explicitly.
    stats, radiation_parameters,
):
    """Updates the radiation variables using the stat_var_update() subroutine.

    References:
        None

    The leading field sequence matches the source arguments. The additional
    radiation diagnostics are explicit source module state, followed by
    ``stats``.
    """
    # Input Variables

    # ---- Begin Code ----
    if not stats.l_sample:
        return stats

    stats = stats.update("Frad", Frad)
    stats = stats.update("radht", radht)
    stats = stats.update("Frad_SW_up", Frad_SW_up)
    stats = stats.update("Frad_LW_up", Frad_LW_up)
    stats = stats.update("Frad_SW_down", Frad_SW_down)
    stats = stats.update("Frad_LW_down", Frad_LW_down)

    if radiation_parameters.rad_scheme == "simplified":
        stats = stats.update("radht_LW", radht_LW)
        stats = stats.update("radht_SW", radht_SW)
        stats = stats.update("Frad_SW", Frad_SW)
        stats = stats.update("Frad_LW", Frad_LW)

    # TODO(port-mirror): The source also updates BUGSrad diagnostics here.
    # BUGSrad is rejected by the initialization boundary for this JAX driver.

    return stats


__all__ = [
    "advance_clubb_radiation", "radiation_driver", "update_radiation_variables",
]

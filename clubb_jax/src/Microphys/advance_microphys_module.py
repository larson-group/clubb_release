"""Hydrometeor transport, mirroring advance_microphys_module.F90.

JAX adaptation: column and level loops are batched; inout state is returned.
"""

from functools import partial

import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import zt2zm, zm2zt, ddzt
from clubb_jax.src.CLUBB_core.constants_clubb import Lv, Cp, rho_lw, rc_tol, ri_tol, zero_threshold
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_driver_api, setup_stats_names
from clubb_jax.src.Microphys import parameters_microphys


# -----------------------------------------------------------------------------
@partial(
    jax.jit,
    static_argnames=(
        'ngrdcol', 'time_current', 'hydromet_dim',
        'tridiag_solve_method', 'fill_holes_type', 'l_upwind_xm_ma',
    ),
)
def advance_microphys(
    gr, ngrdcol, dt, time_current,          # In
    hydromet_dim, hm_metadata,              # In
    wm_zt, wp2,                             # In
    exner, rho, rho_zm, rcm,                # In
    cloud_frac, Kh_zm, Skw_zm,              # In
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,  # In
    hydromet_mc, Ncm_mc, Lscale,            # In
    hydromet_vel_covar_zt_impc,             # In
    hydromet_vel_covar_zt_expc,             # In
    clubb_params, nu_vert_res_dep,          # In
    tridiag_solve_method,                   # In
    fill_holes_type,                        # In
    l_upwind_xm_ma,                         # In
    stats,                                  # InOut
    hydromet, hydromet_vel_zt, hydrometp2,  # InOut
    K_hm, Ncm, Nc_in_cloud, rvm_mc,         # InOut
    thlm_mc, err_info,                      # InOut
):
    """Advance mean precipitating hydrometeors and mean cloud droplet concentration one model time
    step, and calculate some statistics.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        dt: Model timestep duration [s]
        time_current: Current time [s]
        hydromet_dim: Number of precipitating hydrometeor fields.
        hm_metadata: Hydrometeor/PDF names, zero-based species indices and tolerances.
        wm_zt: w wind component on thermodynamic levels [m/s]
        wp2: Variance of vertical velocity (momentum levels) [m^2/s^2]
        exner: Exner function [-]
        rho: Density on thermodynamic levels [kg/m^3]
        rho_zm: Density on momentum levels [kg/m^3]
        rcm: Mean cloud water mixing ratio [kg/kg]
        cloud_frac: Cloud fraction [-]
        Kh_zm: Kh Eddy diffusivity on momentum grid [m^2/s]
        Skw_zm: Skewness of w on momentum levels [-]
        rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
        hydromet_mc: Microphysics tendency for mean hydrometeors [units/s]
        Ncm_mc: Microphysics tendency for Ncm [num/kg/s]
        Lscale: Length-scale [m]
        hydromet_vel_covar_zt_impc: Imp. comp. <V_hm'h_m'> t-levs [m/s]
        hydromet_vel_covar_zt_expc: Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]
        clubb_params: Column-dependent tunable CLUBB parameters.
        nu_vert_res_dep: Resolution-dependent background diffusivities.
        tridiag_solve_method: Specifier for method to solve tridiagonal systems
        fill_holes_type: Option for which type of hole filler to use
        l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
            approximation rather than a centered differencing for turbulent or mean advection
            terms. It affects rtm, thlm, sclrm, um and vm.
        stats: Immutable statistics state; return its updated value.
        hydromet: Hydrometeor mean, <h_m> (thermo. levels) [units]
        hydromet_vel_zt: Mean hydrometeor sed. vel. on thermo. levs. [m/s]
        hydrometp2: Variance of hydrometeor (overall) (m-levs.) [units^2]
        K_hm: hm eddy diffusivity on momentum grid [m^2/s]
        Ncm: Mean cloud droplet conc., <N_c> (thermo. levs.) [num/kg]
        Nc_in_cloud: Mean (in-cloud) cloud droplet concentration [num/kg]
        rvm_mc: Microphysics contributions to vapor water [kg/kg/s]
        thlm_mc: Microphysics contributions to liquid potential temp. [K/s]
        err_info: Per-column error state; return any updated fatal status.
    """

    # Description:
    # Advance mean precipitating hydrometeors and mean cloud droplet
    # concentration one model time step, and calculate some statistics.
    # References:
    # ---------------------------------------------------------------------------

    # Initialize hydrometeor and cloud-number turbulent fluxes on momentum levels.
    wphydrometp = jnp.zeros((gr.ngrdcol, gr.nzm, hydromet_dim))
    wpNcp = jnp.zeros((gr.ngrdcol, gr.nzm))

    # Return until the configured microphysics start time.
    if time_current < parameters_microphys.microphys_start_time:
        return (
            stats,
            hydromet,
            hydromet_vel_zt,
            hydrometp2,
            K_hm,
            Ncm,
            Nc_in_cloud,
            rvm_mc,
            thlm_mc,
            err_info,
            wphydrometp,
            wpNcp,
        )
    if hydromet_dim > 0:
        # Down-gradient closure: <w'hm'> = -K_hm * d<hm>/dz. The diffusivity
        # depends on the species mean and variance, as well as the turbulent state.
        K_hm = calculate_K_hm(
            gr, gr.ngrdcol, wp2, Kh_zm, Skw_zm, Lscale,  # In
            hydromet_dim, hm_metadata.hydromet_tol,      # In
            hydromet, hydrometp2,                        # InOut
            clubb_params,                                # In
            False,                                       # In
        )
        l_prevent_hm_ta_above_cloud = False  # Source local parameter.

        # Optionally suppress turbulent transport above cloud top. The momentum
        # level is above its associated thermodynamic level; exclude cloud ice
        # and ice number, as in the source. This source option is currently false.
        if l_prevent_hm_ta_above_cloud:
            cloud_top_level = get_cloud_top_level(
                gr.nzt, rcm.shape[0], rcm, hydromet,  # In
                hydromet_dim, hm_metadata.iiri,       # In
            )
            for i in range(hydromet_dim):
                if i not in (hm_metadata.iiri, hm_metadata.iiNi):
                    K_hm = K_hm.at[..., i].set(
                        jnp.where(
                            (jnp.arange(gr.nzm)[None, :] >= cloud_top_level[:, None])
                            & (cloud_top_level[:, None] > 0),
                            0.0,
                            K_hm[..., i],
                        )
                    )

    # Sample species-dependent eddy diffusivities before advancing the fields.
    for i in range(hydromet_dim):
        stats = stats.update("K_hm_" + hm_metadata.hydromet_list[i][:2], K_hm[..., i])
    if parameters_microphys.l_predict_Nc:
        # Solve for K_Nc before advancing the precipitating hydrometeors.
        from clubb_jax.src.CLUBB_core.parameter_indices import ic_K_hm

        K_Nc = clubb_params[:, ic_K_hm, None] * Kh_zm

    # Advance precipitating hydrometeors with mean/turbulent advection,
    # sedimentation and the explicit microphysics tendencies.
    if hydromet_dim > 0:
        (
            stats,
            hydromet,
            hydromet_vel_zt,
            hydrometp2,
            rvm_mc,
            thlm_mc,
            err_info,
            wphydrometp,
            hydromet_vel,
            hydromet_vel_covar,
            hydromet_vel_covar_zt,
        ) = advance_hydrometeor(
            gr, gr.ngrdcol, dt, hydromet_dim, hm_metadata,  # In
            wm_zt, exner, cloud_frac, K_hm,                 # In
            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,          # In
            hydromet_mc, hydromet_vel_covar_zt_impc,        # In
            hydromet_vel_covar_zt_expc,                     # In
            nu_vert_res_dep,                                # In
            l_upwind_xm_ma,                                 # In
            tridiag_solve_method,                           # In
            fill_holes_type,                                # In
            stats,                                          # InOut
            hydromet, hydromet_vel_zt,                      # InOut
            hydrometp2, rvm_mc, thlm_mc, err_info,          # InOut
        )

    # JAX adaptation of the source fatal RETURN after advance_hydrometeor.
    # The host emits write_adv_micro_errors after this kernel returns.
    def advance_cloud_number(carry):
        stats, Ncm, Nc_in_cloud, err_info, wpNcp = carry
        # Advance predicted cloud number; prescribed in-cloud number remains
        # constant while its grid mean follows the updated cloud fraction.
        if parameters_microphys.l_predict_Nc:
            stats, Ncm, Nc_in_cloud, err_info, wpNcp = advance_Ncm(
                gr, gr.ngrdcol, dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm,  # In
                rho_ds_zt, invrs_rho_ds_zt, Ncm_mc,                           # In
                nu_vert_res_dep,                                              # In
                l_upwind_xm_ma,                                               # In
                tridiag_solve_method,                                         # In
                stats,                                                        # InOut
                Ncm, Nc_in_cloud, err_info,                                   # InOut
            )
        else:
            # Nc is prescribed; its grid mean depends on cloud fraction.
            Ncm = Nc_in_cloud * cloud_frac
        return stats, Ncm, Nc_in_cloud, err_info, wpNcp

    carry = (stats, Ncm, Nc_in_cloud, err_info, wpNcp)
    if clubb_at_least_debug_level(0):
        stats, Ncm, Nc_in_cloud, err_info, wpNcp = jax.lax.cond(
            err_info.any_fatal(), lambda value: value, advance_cloud_number, carry
        )
    else:
        stats, Ncm, Nc_in_cloud, err_info, wpNcp = advance_cloud_number(carry)

    # Preserve the second source fatal RETURN, after advance_Ncm.
    def accumulate_statistics(stats):
        stats = stats.update("Ncm", Ncm)
        stats = stats.update("Nc_in_cloud", Nc_in_cloud)
        iirr = hm_metadata.iirr
        if iirr >= 0:
            # Rainfall rate on thermodynamic levels [mm/day], positive downward.
            # Include both the mean product <Vrr><rr> and covariance <Vrr'rr'>.
            stats = stats.update(
                "precip_rate_zt",
                jnp.maximum(
                    -(
                        hydromet[..., iirr] * hydromet_vel_zt[..., iirr]
                        + hydromet_vel_covar_zt[..., iirr]
                    ),
                    0.0,
                )
                * (rho / rho_lw)
                * 86400.0
                * 1000.0,
            )

            # Precipitation energy flux on momentum levels [W/m^2],
            # positive downward, with the same mean-plus-covariance transport.
            stats = stats.update(
                "Fprec",
                jnp.maximum(
                    -(
                        zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., iirr])
                        * hydromet_vel[..., iirr]
                        + hydromet_vel_covar[..., iirr]
                    ),
                    0.0,
                )
                * rho_zm
                * Lv,
            )

            # Morrison supplies its own surface precipitation from core fallout.
            if parameters_microphys.microphys_scheme != "morrison":
                stats = stats.update(
                    "precip_rate_sfc",
                    jnp.maximum(
                        -(
                            hydromet[:, 0, iirr] * hydromet_vel_zt[:, 0, iirr]
                            + hydromet_vel_covar_zt[:, 0, iirr]
                        ),
                        0.0,
                    )
                    * (rho[:, 0] / rho_lw)
                    * 86400.0
                    * 1000.0,
                )

            # Surface rain energy flux and interpolated mixing ratio.
            stats = stats.update(
                "rain_flux_sfc",
                jnp.maximum(
                    -(
                        zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., iirr])[:, 0]
                        * hydromet_vel[:, 0, iirr]
                        + hydromet_vel_covar[:, 0, iirr]
                    ),
                    0.0,
                )
                * rho_zm[:, 0]
                * Lv,
            )
            stats = stats.update(
                "rrm_sfc", zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., iirr])[:, 0]
            )
        from clubb_jax.src.CLUBB_core.stats_clubb_utilities import stats_accumulate_hydromet_api

        stats = stats_accumulate_hydromet_api(
            gr, ngrdcol, hydromet_dim,  # In
            hm_metadata, hydromet,      # In
            rho_ds_zt,                  # In
            stats,                      # InOut
        )
        return stats

    if clubb_at_least_debug_level(0):
        stats = jax.lax.cond(
            err_info.any_fatal(), lambda value: value, accumulate_statistics, stats
        )
    else:
        stats = accumulate_statistics(stats)
    # Additional JAX host-safety check: propagate nonfinite hydrometeors as fatal.
    err_info = err_info.set_fatal(mask=jnp.any(~jnp.isfinite(hydromet), axis=(1, 2)))
    return (
        stats,
        hydromet,
        hydromet_vel_zt,
        hydrometp2,
        K_hm,
        Ncm,
        Nc_in_cloud,
        rvm_mc,
        thlm_mc,
        err_info,
        wphydrometp,
        wpNcp,
    )


# -----------------------------------------------------------------------------
def advance_hydrometeor(
    gr, ngrdcol, dt, hydromet_dim, hm_metadata,  # In
    wm_zt, exner, cloud_frac, K_hm,              # In
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,       # In
    hydromet_mc, hydromet_vel_covar_zt_impc,     # In
    hydromet_vel_covar_zt_expc,                  # In
    nu_vert_res_dep,                             # In
    l_upwind_xm_ma,                              # In
    tridiag_solve_method,                        # In
    fill_holes_type,                             # In
    stats,                                       # InOut
    hydromet, hydromet_vel_zt,                   # InOut
    hydrometp2, rvm_mc, thlm_mc, err_info,       # InOut
):
    """Advance each hydrometeor (precipitating hydrometeor) one model time step.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        dt: Duration of one model time step [s]
        hydromet_dim: Number of precipitating hydrometeor fields.
        hm_metadata: Hydrometeor/PDF names, zero-based species indices and tolerances.
        wm_zt: mean w wind component on thermodynamic levels [m/s]
        exner: Exner function, (p/p1000mb)^(Rd/Cp) [-]
        cloud_frac: Cloud fraction [-]
        K_hm: Coefficient of diffusion (turb. adv.) for hydrometeors [m^2/s]
        rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
        hydromet_mc: Change in hydrometeors due to microphysics [units/s]
        hydromet_vel_covar_zt_impc: Imp. comp. <V_hm'h_m'> t-levs [m/s]
        hydromet_vel_covar_zt_expc: Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]
        nu_vert_res_dep: Resolution-dependent background diffusivities.
        l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
            approximation rather than a centered differencing for turbulent or mean advection
            terms. It affects rtm, thlm, sclrm, um and vm.
        tridiag_solve_method: Specifier for method to solve tridiagonal systems
        fill_holes_type: Specifier for which hole filling method to use
        stats: Immutable statistics state; return its updated value.
        hydromet: Hydrometeor mean, <h_m> (thermodynamic levs.) [units]
        hydromet_vel_zt: Mean hydrometeor sed. velocity on thermo. levs. [m/s]
        hydrometp2: Variance of hydrometeor (overall) (m-levs.) [units^2]
        rvm_mc: Microphysics contributions to vapor water [kg/kg/s]
        thlm_mc: Microphysics contributions to liquid potential temp. [K/s]
        err_info: Per-column error state; return any updated fatal status.
    """

    # Description:
    #   Advance each hydrometeor (precipitating hydrometeor) one model time step.
    # References:
    #   None
    # -----------------------------------------------------------------------

    # Allocate momentum-level velocities/fluxes and thermodynamic covariances.
    hydromet_vel = jnp.zeros((gr.ngrdcol, gr.nzm, hydromet_dim))
    ratio_hmp2_on_hmm2 = jnp.zeros_like(hydrometp2)
    wphydrometp = jnp.zeros_like(hydrometp2)
    hydromet_vel_covar = jnp.zeros_like(hydrometp2)
    hydromet_vel_covar_zt = jnp.zeros_like(hydromet)
    for i in range(hydromet_dim):
        max_velocity, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(
            i, hydromet_dim, hm_metadata.hydromet_list
        )

        # Begin the species budget before its implicit transport solve.
        stats = stats.update(name_mc, hydromet_mc[..., i])
        stats = stats.begin_budget(name_bt, hydromet[..., i] / dt)

        # Cap negative fall speeds at the species limit. Interpolate to the
        # momentum grid and impose zero sedimentation through the model top.
        # TODO(port-mirror): source per-level fall-speed warnings are omitted
        # here; restore them when ordered device diagnostics support this report.
        hydromet_vel_zt = hydromet_vel_zt.at[..., i].set(
            jnp.clip(hydromet_vel_zt[..., i], max_velocity, zero_threshold)
        )
        hydromet_vel = hydromet_vel.at[..., i].set(
            zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_zt[..., i]).at[:, -1].set(0.0)
        )

        # No variance equation is advanced here. Save <hm'^2>/<hm>^2 and
        # reconstruct the variance after the mean changes; a negligible mean
        # has no defined ratio and uses zero.
        hydromet_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., i])
        ratio_hmp2_on_hmm2 = ratio_hmp2_on_hmm2.at[..., i].set(
            jnp.where(
                hydromet_zm > hm_metadata.hydromet_tol[i],
                hydrometp2[..., i]
                / jnp.where(hydromet_zm > hm_metadata.hydromet_tol[i], hydromet_zm**2, 1.0),
                0.0,
            )
        )

        # Crank-Nicholson first half of <w'hm'>, evaluated at timestep t.
        # Impose zero turbulent flux at both domain boundaries.
        K_hm_nu_hm = K_hm[..., i] + nu_vert_res_dep.nu_hm[:, None]
        xpwp = K_hm_nu_hm * ddzt(gr.nzm, gr.nzt, ngrdcol, gr, hydromet[..., i])
        wphydrometp = wphydrometp.at[..., i].set((-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0))
        # Assemble implicit terms, explicit terms, then solve each species.
        stats, lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, lhs = microphys_lhs(
            gr, gr.ngrdcol, hm_metadata.hydromet_list[i],              # In
            parameters_microphys.l_hydromet_sed[i], dt, K_hm[..., i],  # In
            nu_vert_res_dep.nu_hm,                                     # In
            wm_zt,                                                     # In
            hydromet_vel[..., i], hydromet_vel_zt[..., i],             # In
            hydromet_vel_covar_zt_impc[..., i],                        # In
            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,                     # In
            l_upwind_xm_ma,                                            # In
            stats,                                                     # InOut
        )
        stats, rhs = microphys_rhs(
            gr, gr.ngrdcol, hm_metadata.hydromet_list[i], dt,  # In
            parameters_microphys.l_hydromet_sed[i],            # In
            hydromet[..., i], hydromet_mc[..., i],             # In
            K_hm[..., i], nu_vert_res_dep.nu_hm, cloud_frac,   # In
            hydromet_vel_covar_zt_expc[..., i],                # In
            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,             # In
            stats,                                             # InOut
        )
        stats, lhs, rhs, hmm, err_info = microphys_solve(
            gr, gr.ngrdcol, hm_metadata.hydromet_list[i],  # In
            parameters_microphys.l_hydromet_sed[i],        # In
            lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,    # In
            cloud_frac,                                    # In
            tridiag_solve_method,                          # In
            stats,                                         # InOut
            lhs, rhs, hydromet[..., i], err_info,          # InOut
        )
        hydromet = hydromet.at[..., i].set(hmm)

    # Now that all species have advanced, fill holes in the profiles.
    stats, thlm_mc, rvm_mc, hydromet = fill_holes_driver_api(
        gr, ngrdcol, gr.nzt, dt, hydromet_dim,  # In
        hm_metadata, True,                      # In
        rho_ds_zt, exner, fill_holes_type,      # In
        stats, thlm_mc, rvm_mc, hydromet,       # InOut
    )

    for i in range(hydromet_dim):
        name = hm_metadata.hydromet_list[i]

        # Sub-tolerance hydrometeors at the lower boundary sediment out of the
        # domain rather than being conserved in the atmospheric column.
        hydromet = hydromet.at[:, 0, i].set(
            jnp.where(
                hydromet[:, 0, i] < hm_metadata.hydromet_tol[i], zero_threshold, hydromet[:, 0, i]
            )
        )

        # Reconstruct the overall variance using the saved variance/mean ratio.
        hydromet_zm = jnp.maximum(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., i]), 0.0)
        hydrometp2 = hydrometp2.at[..., i].set(ratio_hmp2_on_hmm2[..., i] * hydromet_zm**2)
        # Source overwrites the first half-step flux with this second half.
        xpwp = (K_hm[..., i] + nu_vert_res_dep.nu_hm[:, None]) * ddzt(
            gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., i]
        )
        wphydrometp = wphydrometp.at[..., i].set((-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0))

        # Recover <Vhm'hm'> from its implicit coefficient and explicit part,
        # then interpolate to momentum levels with zero flux through the top.
        hydromet_vel_covar_zt = hydromet_vel_covar_zt.at[..., i].set(
            hydromet_vel_covar_zt_impc[..., i] * hydromet[..., i]
            + hydromet_vel_covar_zt_expc[..., i]
        )
        hydromet_vel_covar = hydromet_vel_covar.at[..., i].set(
            zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_covar_zt[..., i]).at[:, -1].set(0.0)
        )

        # Variance/covariance diagnostics and the completed mean-field budget.
        stats = stats.update(name[:2] + "p2", hydrometp2[..., i])
        stats = stats.update("wp" + name[:2] + "p", wphydrometp[..., i])
        if name == "rrm":
            stats = stats.update("Vrrprrp", hydromet_vel_covar[..., i])
        elif name == "Nrm":
            stats = stats.update("VNrpNrp", hydromet_vel_covar[..., i])
        stats = stats.finalize_budget(name + "_bt", hydromet[..., i] / dt)
    return (
        stats,
        hydromet,
        hydromet_vel_zt,
        hydrometp2,
        rvm_mc,
        thlm_mc,
        err_info,
        wphydrometp,
        hydromet_vel,
        hydromet_vel_covar,
        hydromet_vel_covar_zt,
    )


# -----------------------------------------------------------------------------
def advance_Ncm(
    gr, ngrdcol, dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm,  # In
    rho_ds_zt, invrs_rho_ds_zt, Ncm_mc,                        # In
    nu_vert_res_dep,                                           # In
    l_upwind_xm_ma,                                            # In
    tridiag_solve_method,                                      # In
    stats,                                                     # InOut
    Ncm, Nc_in_cloud, err_info,                                # InOut
):
    """Advance cloud droplet concentration (Ncm) one model time step.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        dt: Duration of one model time step [s]
        wm_zt: mean w wind component on thermodynamic levels [m/s]
        cloud_frac: Cloud fraction [-]
        K_Nc: Coefficient of diffusion (turb. adv.) for Nc [m^2/s]
        rcm: Mean cloud water mixing ratio [kg/kg]
        rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
        Ncm_mc: Change in Ncm due to microphysics [num/kg/s]
        nu_vert_res_dep: Resolution-dependent background diffusivities.
        l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
            approximation rather than a centered differencing for turbulent or mean advection
            terms. It affects rtm, thlm, sclrm, um and vm.
        tridiag_solve_method: Specifier for method to solve tridiagonal systems
        stats: Immutable statistics state; return its updated value.
        Ncm: Mean cloud droplet conc., <N_c> (thermo. levs.) [num/kg]
        Nc_in_cloud: Mean (in-cloud) cloud droplet concentration [num/kg]
        err_info: Per-column error state; return any updated fatal status.
    """

    # Description:
    # Advance cloud droplet concentration (Ncm) one model time step.
    # References:
    # -----------------------------------------------------------------------

    from clubb_jax.src.CLUBB_core.constants_clubb import (
        cloud_frac_min,
        pi,
        mvr_cloud_max,
        Nc_in_cloud_min,
    )

    # Cloud number does not sediment in this solve. Zero velocity/covariance
    # arrays retain the common transport interface with l_sed=False.
    Ncm_vel_covar_zt_impc = Ncm_vel_covar_zt_expc = jnp.zeros_like(Ncm)
    Ncm_vel = jnp.zeros_like(K_Nc)
    Ncm_vel_zt = jnp.zeros_like(Ncm)

    # Begin the grid-mean number budget, including the cloud fraction when
    # diffusion advances in-cloud number instead of grid-mean number.
    stats = stats.update("Ncm_mc", Ncm_mc)
    if parameters_microphys.l_in_cloud_Nc_diff:
        stats = stats.begin_budget(
            "Ncm_bt", Nc_in_cloud * jnp.maximum(cloud_frac, cloud_frac_min) / dt
        )
    else:
        stats = stats.begin_budget("Ncm_bt", Ncm / dt)

    # Crank-Nicholson first half of <w'Nc'> at timestep t; zero boundary fluxes.
    xpwp = (K_Nc + nu_vert_res_dep.nu_hm[:, None]) * ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ncm)
    wpNcp = (-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0)
    stats, lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, lhs = microphys_lhs(
        gr, gr.ngrdcol, "Ncm", False, dt, K_Nc, nu_vert_res_dep.nu_hm, wm_zt,  # In
        Ncm_vel, Ncm_vel_zt,                                                   # In
        Ncm_vel_covar_zt_impc,                                                 # In
        rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,                                 # In
        l_upwind_xm_ma,                                                        # In
        stats,                                                                 # InOut
    )
    if parameters_microphys.l_in_cloud_Nc_diff:
        stats, rhs = microphys_rhs(
            gr, gr.ngrdcol, "Ncm", dt, False,                               # In
            Nc_in_cloud, Ncm_mc / jnp.maximum(cloud_frac, cloud_frac_min),  # In
            K_Nc, nu_vert_res_dep.nu_hm, cloud_frac,                        # In
            Ncm_vel_covar_zt_expc,                                          # In
            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,                          # In
            stats,                                                          # InOut
        )
    else:
        stats, rhs = microphys_rhs(
            gr, gr.ngrdcol, "Ncm", dt, False,         # In
            Ncm, Ncm_mc,                              # In
            K_Nc, nu_vert_res_dep.nu_hm, cloud_frac,  # In
            Ncm_vel_covar_zt_expc,                    # In
            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,    # In
            stats,                                    # InOut
        )

    # Advance either Nc_in_cloud or Ncm, then recover the other using the
    # same minimum cloud fraction as the RHS and budget diagnostics.
    if parameters_microphys.l_in_cloud_Nc_diff:
        stats, lhs, rhs, Nc_in_cloud, err_info = microphys_solve(
            gr, gr.ngrdcol, "Ncm", False,                # In
            lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,  # In
            cloud_frac,                                  # In
            tridiag_solve_method,                        # In
            stats,                                       # InOut
            lhs, rhs, Nc_in_cloud, err_info,             # InOut
        )
        Ncm = Nc_in_cloud * jnp.maximum(cloud_frac, cloud_frac_min)
    else:
        stats, lhs, rhs, Ncm, err_info = microphys_solve(
            gr, gr.ngrdcol, "Ncm", False,                # In
            lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,  # In
            cloud_frac,                                  # In
            tridiag_solve_method,                        # In
            stats,                                       # InOut
            lhs, rhs, Ncm, err_info,                     # InOut
        )
        Nc_in_cloud = Ncm / jnp.maximum(cloud_frac, cloud_frac_min)

    # JAX adaptation of the fatal RETURN after microphys_solve: the callable
    # keeps clipping, fluxes, and budget updates entirely on the success path.
    def finish_cloud_number(carry):
        stats, Ncm, Nc_in_cloud, wpNcp = carry
        # Enforce the number minimum implied by the maximum cloud-drop mean
        # volume radius, as well as the specified minimum in-cloud number.
        # TODO(port-mirror): source per-level below-minimum warnings are omitted;
        # restore them alongside the ordered device clipping diagnostics.
        Ncm_mvr_min = (1.0 / ((4.0 / 3.0) * pi * rho_lw * mvr_cloud_max**3)) * rcm
        Ncm_min = jnp.maximum(
            Nc_in_cloud_min * jnp.maximum(cloud_frac, cloud_frac_min), Ncm_mvr_min
        )
        stats = stats.begin_budget("Ncm_cl", Ncm / dt)
        Ncm = jnp.maximum(Ncm, Ncm_min)
        stats = stats.finalize_budget("Ncm_cl", Ncm / dt)
        Ncic_min = jnp.maximum(
            Nc_in_cloud_min, Ncm_mvr_min / jnp.maximum(cloud_frac, cloud_frac_min)
        )
        Nc_in_cloud = jnp.maximum(Nc_in_cloud, Ncic_min)

        # Covariance half-step at t+1. As in the source, this overwrites the
        # first half-step flux; zero turbulent flux at both boundaries.
        xpwp = (K_Nc + nu_vert_res_dep.nu_hm[:, None]) * ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ncm)
        wpNcp = (-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0)
        stats = stats.update("wpNcp", wpNcp)
        stats = stats.finalize_budget("Ncm_bt", Ncm / dt)
        return stats, Ncm, Nc_in_cloud, wpNcp

    carry = (stats, Ncm, Nc_in_cloud, wpNcp)
    if clubb_at_least_debug_level(0):
        stats, Ncm, Nc_in_cloud, wpNcp = jax.lax.cond(
            err_info.any_fatal(), lambda value: value, finish_cloud_number, carry
        )
    else:
        stats, Ncm, Nc_in_cloud, wpNcp = finish_cloud_number(carry)
    return stats, Ncm, Nc_in_cloud, err_info, wpNcp


# -----------------------------------------------------------------------------
def microphys_solve(
    gr, ngrdcol, solve_type, l_sed,              # In
    lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,  # In
    cloud_frac,                                  # In
    tridiag_solve_method,                        # In
    stats,                                       # InOut
    lhs, rhs, hmm, err_info,                     # InOut
):
    """Solve the tridiagonal system for hydrometeor variable.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        solve_type: Description of which hydrometeor is being solved for.
        l_sed: Whether to add a hydrometeor sedimentation term.
        lhs_ta: LHS corresponding to contribution from turbulent adv. [1/s]
        lhs_ma: LHS corresponding to contribution from mean advection [1/s]
        sed_turb_lhs: Implicit turbulent-sedimentation super/main/subdiagonal terms [1/s].
        sed_diff_lhs: Implicit mean-sedimentation super/main/subdiagonal terms [1/s].
        cloud_frac: Cloud fraction (thermodynamic levels) [-]
        tridiag_solve_method: Specifier for method to solve tridiagonal systems
        stats: Immutable statistics state; return its updated value.
        lhs: Left hand side
        rhs: Right hand side vector
        hmm: Mean value of hydrometeor (thermodynamic levels) [units vary]
        err_info: Per-column error state; return any updated fatal status.
    """

    # Description:
    # Solve the tridiagonal system for hydrometeor variable.
    # References:
    #  None
    # ---------------------------------------------------------------------------

    from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve
    from clubb_jax.src.CLUBB_core.constants_clubb import cloud_frac_min

    # Solve system using a tridiag_solve.
    err_info, hmm, _ = tridiag_solve(
        solve_type, tridiag_solve_method, gr.ngrdcol, gr.nzt, lhs, rhs, err_info
    )

    # JAX adaptation: lax.cond needs a callable for the source's fatal RETURN.
    # No implicit statistics may be committed after a failed tridiagonal solve.
    # TODO(port-mirror): low-level source stderr call-stack banners are folded
    # into the standalone host's fatal report. Preserve individual banners when
    # shared device error diagnostics support ordered routine-context messages.
    def accumulate_implicit_statistics(stats):
        # Statistics: implicit contributions to hydrometeor hmm.
        # Tridiagonal storage: super/main/sub are JAX slots 0/1/2,
        # corresponding to Fortran 1/2/3. Clamp the neighboring boundary levels.
        km1 = jnp.maximum(jnp.arange(gr.nzt) - 1, 0)
        kp1 = jnp.minimum(jnp.arange(gr.nzt) + 1, gr.nzt - 1)
        if not stats.l_sample:
            return stats

        # In-cloud transport terms must be scaled by cloud fraction to balance
        # the grid-mean cloud-number budget.
        if solve_type == "Ncm" and parameters_microphys.l_in_cloud_Nc_diff:
            ma_term = (
                -lhs_ma[2] * hmm[:, km1] * jnp.maximum(cloud_frac, cloud_frac_min)
                - lhs_ma[1] * hmm * jnp.maximum(cloud_frac, cloud_frac_min)
                - lhs_ma[0] * hmm[:, kp1] * jnp.maximum(cloud_frac, cloud_frac_min)
            )
            ta_term = (
                -lhs_ta[2] * hmm[:, km1] * jnp.maximum(cloud_frac, cloud_frac_min)
                - lhs_ta[1] * hmm * jnp.maximum(cloud_frac, cloud_frac_min)
                - lhs_ta[0] * hmm[:, kp1] * jnp.maximum(cloud_frac, cloud_frac_min)
            )
        else:
            ma_term = -lhs_ma[2] * hmm[:, km1] - lhs_ma[1] * hmm - lhs_ma[0] * hmm[:, kp1]
            ta_term = -lhs_ta[2] * hmm[:, km1] - lhs_ta[1] * hmm - lhs_ta[0] * hmm[:, kp1]
        sd_term = (
            -sed_diff_lhs[2] * hmm[:, km1] - sed_diff_lhs[1] * hmm - sed_diff_lhs[0] * hmm[:, kp1]
        )
        ts_term = (
            -sed_turb_lhs[2] * hmm[:, km1] - sed_turb_lhs[1] * hmm - sed_turb_lhs[0] * hmm[:, kp1]
        )
        if solve_type in ("rrm", "Nrm", "rim", "rsm", "rgm", "Ncm", "Nim", "Nsm", "Ngm"):
            stats = stats.update(solve_type + "_ma", ma_term)
            if l_sed and solve_type != "Ncm":
                stats = stats.update(solve_type + "_sd", sd_term)
            if l_sed and solve_type in ("rrm", "Nrm"):
                stats = stats.finalize_budget(solve_type + "_ts", ts_term)
            stats = stats.finalize_budget(solve_type + "_ta", ta_term)
        return stats

    if clubb_at_least_debug_level(0):
        stats = jax.lax.cond(
            err_info.any_fatal(), lambda value: value, accumulate_implicit_statistics, stats
        )
    else:
        stats = accumulate_implicit_statistics(stats)
    return stats, lhs, rhs, hmm, err_info


# -----------------------------------------------------------------------------
def microphys_lhs(
    gr, ngrdcol, solve_type, l_sed, dt, K_hm, nu, wm_zt,  # In
    V_hm, V_hmt,                                          # In
    Vhmphmp_zt_impc,                                      # In
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,                # In
    l_upwind_xm_ma,                                       # In
    stats,                                                # InOut
):
    """Setup the matrix of implicit contributions to a term.

    Can include the effects of sedimentation, diffusion, and advection. The Morrison microphysics
    has an explicit sedimentation code, which is handled elsewhere.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        solve_type: Description of which hydrometeor is being solved for.
        l_sed: Whether to add a hydrometeor sedimentation term.
        dt: Model timestep [s]
        K_hm: Coefficient of diffusion (turb. adv.) for hydrometeor [m^2/s]
        nu: Background diffusion coefficient [m^2/s]
        wm_zt: w wind component on thermodynamic levels [m/s]
        V_hm: Sedimentation velocity of hydrometeor (momentum levels) [m/s]
        V_hmt: Sedimentation velocity of hydrometeor (thermo. levels) [m/s]
        Vhmphmp_zt_impc: Implicit comp. of <V_hm'h_m'> on t-levs [units(m/s)]
        rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
        l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
            approximation rather than a centered differencing for turbulent or mean advection
            terms. It affects rtm, thlm, sclrm, um and vm.
        stats: Immutable statistics state; return its updated value.
    """

    # Description:
    # Setup the matrix of implicit contributions to a term.
    # Can include the effects of sedimentation, diffusion, and advection.
    # The Morrison microphysics has an explicit sedimentation code, which is
    # handled elsewhere.
    #
    # Notes:
    # Setup for tridiagonal system and boundary conditions should be the same as
    # the original rain subroutine code.
    # -----------------------------------------------------------------------

    from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
    from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs

    # Interpolate the implicit sedimentation coefficient from thermodynamic
    # levels to momentum levels, where the vertical flux is differenced.
    Vhmphmp_impc = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Vhmphmp_zt_impc)

    # Turbulent advection uses <w'hm'> = -K_hm*d<hm>/dz. The implicit
    # diffusion operator carries half the Crank-Nicholson contribution.
    Kh_zm = K_hm
    Kh_zt = jnp.maximum(zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, K_hm), 0.0)
    lhs_ta = 0.5 * diffusion_zt_lhs(
        gr.nzm, gr.nzt, gr.ngrdcol, gr, Kh_zm, Kh_zt, nu,  # In
        invrs_rho_ds_zt, rho_ds_zm,                        # In
    )

    # Apply the lower boundary at zero-based level 0: negative superdiagonal,
    # positive main diagonal, and zero subdiagonal.
    bc = (
        0.5
        * invrs_rho_ds_zt[:, 0]
        * (gr.invrs_dzt[:, 0] * (Kh_zm[:, 1] + nu) * rho_ds_zm[:, 1] * gr.invrs_dzm[:, 1])
    )
    lhs_ta = lhs_ta.at[0, :, 0].set(-bc).at[1, :, 0].set(bc).at[2, :, 0].set(0.0)
    # LHS mean advection term.
    lhs_ma = term_ma_zt_lhs(
        gr.nzm, gr.nzt, gr.ngrdcol, wm_zt, gr.weights_zt2zm,  # In
        gr.invrs_dzt, gr.invrs_dzm,                           # In
        l_upwind_xm_ma, gr.grid_dir,                          # In
    )
    # JAX adaptation: level arguments are zero-based vectors, batching the source loop.
    k = jnp.arange(gr.nzt)
    kp1 = jnp.minimum(k + 1, gr.nzt - 1)
    # Time tendency plus implicit turbulent and mean advection.
    lhs = jnp.zeros_like(lhs_ta).at[1].set(1.0 / dt)
    lhs = lhs + lhs_ta
    lhs = lhs + lhs_ma

    # Retain separate sedimentation operators for the budget diagnostics.
    # JAX computes each operator once; the source repeats it for statistics.
    sed_diff_lhs = sed_turb_lhs = jnp.zeros_like(lhs)
    if l_sed:
        # Morrison sedimentation is handled within its core, via the _mc terms.
        if not parameters_microphys.l_upwind_diff_sed:
            sed_diff_lhs = sed_centered_diff_lhs(
                gr, V_hm[:, kp1], V_hm[:, k], rho_ds_zm[:, kp1],  # In
                rho_ds_zm[:, k], invrs_rho_ds_zt,                 # In
                gr.invrs_dzt, k,                                  # In
            )
        else:
            sed_diff_lhs = sed_upwind_diff_lhs(
                gr, V_hmt, V_hmt[:, kp1], rho_ds_zt,  # In
                rho_ds_zt[:, kp1], invrs_rho_ds_zt,   # In
                gr.invrs_dzm[:, kp1], k,              # In
            )
        lhs = lhs + sed_diff_lhs
        sed_turb_lhs = term_turb_sed_lhs(
            gr, Vhmphmp_impc[:, kp1], Vhmphmp_impc[:, k],  # In
            Vhmphmp_zt_impc[:, kp1], Vhmphmp_zt_impc,      # In
            rho_ds_zm[:, kp1], rho_ds_zm[:, k],            # In
            rho_ds_zt[:, kp1], rho_ds_zt,                  # In
            gr.invrs_dzt, gr.invrs_dzm[:, kp1],            # In
            invrs_rho_ds_zt, k,                            # In
        )
        lhs = lhs + sed_turb_lhs
    return stats, lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, lhs


# -----------------------------------------------------------------------------
def microphys_rhs(
    gr, ngrdcol, solve_type, dt, l_sed,     # In
    hmm, hmm_tndcy,                         # In
    K_hm, nu, cloud_frac,                   # In
    Vhmphmp_zt_expc,                        # In
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,  # In
    stats,                                  # InOut
):
    """Compute RHS vector for a given hydrometeor.

    This subroutine computes the explicit portion of the predictive equation for a given
    hydrometeor.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        solve_type: Description of which hydrometeor is being solved for.
        dt: Duration of model timestep [s]
        l_sed: Flag for hydrometeor sedimentation
        hmm: Mean value of hydrometeor (t-levs.) [units]
        hmm_tndcy: Microphysics tendency (thermo. levels) [units/s]
        K_hm: Coef. of diffusion for hydrometeor [m^2/s]
        nu: Background diffusion coefficient [m^2/s]
        cloud_frac: Cloud fraction [-]
        Vhmphmp_zt_expc: Explicit comp. of <V_hm'h_m'> on t-levs [units(m/s)]
        rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
        stats: Immutable statistics state; return its updated value.
    """

    # Description:
    # Compute RHS vector for a given hydrometeor.
    # This subroutine computes the explicit portion of the predictive equation
    # for a given hydrometeor.
    # References:
    # -----------------------------------------------------------------------

    from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
    from clubb_jax.src.CLUBB_core.constants_clubb import cloud_frac_min

    # Interpolate the explicit sedimentation flux to momentum levels.
    Vhmphmp_expc = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Vhmphmp_zt_expc)

    # Explicit Crank-Nicholson contribution uses the same half-weighted
    # diffusion operator and lower boundary as the implicit LHS.
    Kh_zm = K_hm
    Kh_zt = jnp.maximum(zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, K_hm), 0.0)
    lhs_ta = 0.5 * diffusion_zt_lhs(
        gr.nzm, gr.nzt, gr.ngrdcol, gr, Kh_zm, Kh_zt, nu,  # In
        invrs_rho_ds_zt, rho_ds_zm,                        # In
    )
    # The lower boundary condition needs to be applied here at level 1.
    bc = (
        0.5
        * invrs_rho_ds_zt[:, 0]
        * (gr.invrs_dzt[:, 0] * (Kh_zm[:, 1] + nu) * rho_ds_zm[:, 1] * gr.invrs_dzm[:, 1])
    )
    lhs_ta = lhs_ta.at[0, :, 0].set(-bc).at[1, :, 0].set(bc).at[2, :, 0].set(0.0)
    k = jnp.arange(gr.nzt)
    km1, kp1 = jnp.maximum(k - 1, 0), jnp.minimum(k + 1, gr.nzt - 1)

    # Explicit time tendency, microphysics rates (auto/accr/evap/etc.), and
    # the old-time turbulent-advection contribution.
    rhs = hmm / dt
    rhs = rhs + hmm_tndcy
    rhs = rhs - lhs_ta[2] * hmm[:, km1] - lhs_ta[1] * hmm - lhs_ta[0] * hmm[:, kp1]

    # Add the explicit turbulent-sedimentation divergence when enabled.
    ts_rhs = jnp.zeros_like(hmm)
    if l_sed:
        ts_rhs = term_turb_sed_rhs(
            gr, Vhmphmp_expc[:, kp1], Vhmphmp_expc[:, k],  # In
            Vhmphmp_zt_expc[:, kp1], Vhmphmp_zt_expc,      # In
            rho_ds_zm[:, kp1], rho_ds_zm[:, k],            # In
            rho_ds_zt[:, kp1], rho_ds_zt,                  # In
            gr.invrs_dzt, gr.invrs_dzm[:, kp1],            # In
            invrs_rho_ds_zt, k,                            # In
        )
        rhs = rhs + ts_rhs

    # Begin the explicit portion of the budgets; microphys_solve later
    # finalizes these with the new-time implicit contributions.
    ta_rhs = lhs_ta[2] * hmm[:, km1] + lhs_ta[1] * hmm + lhs_ta[0] * hmm[:, kp1]
    if solve_type == "Ncm" and parameters_microphys.l_in_cloud_Nc_diff:
        ta_rhs = ta_rhs * jnp.maximum(cloud_frac, cloud_frac_min)
    if solve_type in ("rrm", "Nrm", "rim", "rsm", "rgm", "Ncm", "Nim", "Nsm", "Ngm"):
        stats = stats.begin_budget(solve_type + "_ta", ta_rhs)
        if l_sed and solve_type in ("rrm", "Nrm"):
            # Reverse the RHS sign when sampling the sedimentation budget.
            stats = stats.begin_budget(solve_type + "_ts", -ts_rhs)
    return stats, rhs


# -----------------------------------------------------------------------------
def sed_centered_diff_lhs(
    gr, V_hmp1, V_hm, rho_ds_zmp1,  # In
    rho_ds_zm, invrs_rho_ds_zt,     # In
    invrs_dzt, level,               # In
):
    """Mean sedimentation of a hydrometeor:  implicit portion of the code, using the centered
    difference approximation to the vertical derivative.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        V_hmp1: Sedimentation velocity of hydrometeor (k+1) [m/s]
        V_hm: Sedimentation velocity of hydrometeor (k) [m/s]
        rho_ds_zmp1: Dry, static density at momentum level (k+1) [kg/m^3]
        rho_ds_zm: Dry, static density at momentum level (k) [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. level (k) [m^3/kg]
        invrs_dzt: Inverse of grid spacing (k) [m]
        level: Zero-based thermodynamic level index (batched over levels).
    """

    # Description:
    # Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    # the centered difference approximation to the vertical derivative.
    #
    # The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    # stands for the sedimentation velocity of the aforementioned hydrometeor.
    #
    # The d(hm)/dt equation contains a sedimentation term:
    #
    # - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    #
    # The variables hm and V_hm in the sedimentation term are divided into mean
    # and turbulent components, and the term is averaged, resulting in:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    # - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    #
    # The mean sedimentation term in the d<hm>/dt equation is:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    #
    # This term is solved for completely implicitly, such that:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    #
    # Note:  When the term is brought over to the left-hand side, the sign is
    #        reversed and the leading "-" in front of the term is changed to
    #        a "+".
    #
    # Timestep index (t) stands for the index of the current timestep, while
    # timestep index (t+1) stands for the index of the next timestep, which is
    # being advanced to in solving the d<hm>/dt equation.
    #
    # This term is discretized as follows when using the centered-difference
    # approximation:
    #
    # The values of <hm> are found on the thermodynamic levels, while the values
    # of <V_hm> are found on the momentum levels.  Additionally, the values of
    # rho_ds_zm are found on the momentum levels, and the values of
    # invrs_rho_ds_zt are found on the thermodynamic levels.  The variable <hm>
    # is interpolated to the intermediate momentum levels.  At the intermediate
    # momentum levels, the interpolated values of <hm> are multiplied by the
    # values of <V_hm> and the values of rho_ds_zm.  Then, the derivative of
    # (rho_ds*<V_hm>*<hm>) is taken over the central thermodynamic level, where
    # it is multiplied by invrs_rho_ds_zt.
    #
    # -----hmp1------------------------------------------------ t(k+1)
    #
    # =============hm(interp)=====V_hm=====rho_ds_zm=========== m(k+1)
    #
    # -----hm--------invrs_rho_ds_zt----d(rho_ds*V_hm*hm)/dz--- t(k)
    #
    # =============hm(interp)=====V_hmm1===rho_ds_zmm1========= m(k)
    #
    # -----hmm1------------------------------------------------ t(k-1)
    #
    # The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    # with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1),
    # respectively.  The letter "t" is used for thermodynamic levels and the
    # letter "m" is used for momentum levels.
    #
    # invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
    #
    #
    # Conservation Properties:
    #
    # When a hydrometeor is sedimented to the ground (or out the lower boundary
    # of the model), it is removed from the atmosphere (or from the model
    # domain).  Thus, the quantity of the hydrometeor over the entire vertical
    # domain should not be conserved due to the process of sedimentation.  Thus,
    # not all of the column totals in the left-hand side matrix should be equal
    # to 0. Instead, the sum of all the column totals should equal the flux of
    # <hm> out the bottom (zm(1) level) of the domain,
    # -rho_ds_zm(1) * V_hm(1) * hm(1), where the value of the hydrometeor at
    # the surface (zm level 1) is set equal to the value of the hydrometeor at
    # zt level 1, which is hm(1). Furthermore, most of the individual column
    # totals should sum to 0, but the 1st and 2nd (from the left) columns should
    # combine to sum to the flux out the bottom of the domain.
    #
    # To see that this modified conservation law is satisfied, compute the
    # sedimentation of hm and integrate vertically.  In discretized matrix
    # notation (where "i" stands for the matrix column and "j" stands for the
    # matrix row):
    #
    # - rho_ds_zm(1) * V_hm(1) * hm(1)
    # = Sum_j Sum_i
    #   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzt )_i
    #   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    #
    # The left-hand side matrix,
    # ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    # partially written below.  The sum over i in the above equation removes
    # invrs_rho_ds_zt and invrs_dzt everywhere from the matrix below.  The sum
    # over j leaves the column totals and the flux at zm(1) that are desired.
    #
    # Left-hand side matrix contributions from the sedimentation term (only);
    # first four vertical levels:
    #
    #     -------------------------------------------------------------------->
    # k=1 |   +invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)            0
    #    |    *invrs_dzt(k)        *invrs_dzt(k)
    #    |    *[ rho_ds_zm(k+1)    *rho_ds_zm(k+1)
    #    |       *V_hm(k+1)*B(k)   *V_hm(k+1)*A(k)
    #    |      -rho_ds_zm(k)
    #    |       *V_hm(k) ]
    #    |
    # k=2 |   -invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    #    |    *invrs_dzt(k)        *invrs_dzt(k)          *invrs_dzt(k)
    #    |    *rho_ds_zm(k)        *[ rho_ds_zm(k+1)      *rho_ds_zm(k+1)
    #    |    *V_hm(k)*D(k)           *V_hm(k+1)*B(k)     *V_hm(k+1)*A(k)
    #    |                           -rho_ds_zm(k)
    #    |                            *V_hm(k)*C(k) ]
    #    |
    # k=3 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    #    |                         *invrs_dzt(k)          *invrs_dzt(k)
    #    |                         *rho_ds_zm(k)          *[ rho_ds_zm(k+1)
    #    |                         *V_hm(k)*D(k)             *V_hm(k+1)*B(k)
    #    |                                                  -rho_ds_zm(k)
    #    |                                                   *V_hm(k)*C(k) ]
    #    |
    # k=4 |           0                     0             -invrs_rho_ds_zt(k)
    #    |                                                *invrs_dzt(k)
    #    |                                                *rho_ds_zm(k)
    #    |                                                *V_hm(k)*D(k)
    #    |
    #   \ /
    #
    # The variables A(k), B(k), C(k), and D(k) are weights of interpolation
    # around the central thermodynamic level (k), such that:
    #
    # A(k) = ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ),
    # B(k) = 1 - [ ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
    #      = 1 - A(k);
    # C(k) = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
    # D(k) = 1 - [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    #      = 1 - C(k).
    #
    # Furthermore, for all intermediate thermodynamic grid levels (as long as
    # k /= gr%nz and k /= 1), the four weighting factors have the following
    # relationships:  A(k) = C(k+1) and B(k) = D(k+1).
    #
    # Note:  The superdiagonal term from level 3 and both the main diagonal
    #        and superdiagonal terms from level 4 are not shown on this
    #        diagram.
    # References:
    # None
    # Notes:
    #   Both COAMPS Microphysics and Brian Griffin's implementation use
    #   Khairoutdinov and Kogan (2000) for the calculation of rain
    #   mixing ratio and rain droplet number concentration sedimentation
    #   velocities, but COAMPS has only the local parameterization.
    # -----------------------------------------------------------------------

    # JAX level indices are zero based. Momentum level k+1 lies between
    # thermodynamic k and k+1; momentum k lies between k-1 and k.
    mkp1, mk = level + 1, level
    # Lower boundary: surface hydrometeor equals thermodynamic level 1.
    # Upper boundary: no flux through the model top.
    # Superdiagonal: no incoming sedimentation at the upper boundary.
    super = jnp.where(
        level == gr.nzt - 1,
        0.0,
        invrs_rho_ds_zt * invrs_dzt * rho_ds_zmp1 * V_hmp1 * gr.weights_zt2zm[:, mkp1, 0],
    )

    # Main diagonal includes surface outflow and interior interpolation.
    main = (
        invrs_rho_ds_zt
        * invrs_dzt
        * (
            jnp.where(level == gr.nzt - 1, 0.0, rho_ds_zmp1 * V_hmp1 * gr.weights_zt2zm[:, mkp1, 1])
            - rho_ds_zm * V_hm * jnp.where(level == 0, 1.0, gr.weights_zt2zm[:, mk, 0])
        )
    )

    # Subdiagonal vanishes at the lower boundary.
    sub = jnp.where(
        level == 0,
        0.0,
        -invrs_rho_ds_zt * invrs_dzt * rho_ds_zm * V_hm * gr.weights_zt2zm[:, mk, 1],
    )
    return jnp.stack((super, main, sub))


# -----------------------------------------------------------------------------
def sed_upwind_diff_lhs(
    gr, V_hmt, V_hmtp1, rho_ds_zt,  # In
    rho_ds_ztp1, invrs_rho_ds_zt,   # In
    invrs_dzmp1, level,             # In
):
    """Mean sedimentation of a hydrometeor:  implicit portion of the code, using the "upwind"
    difference approximation to the vertical derivative.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        V_hmt: Sed. velocity of hydrometeor at t-lev (k) [m/s]
        V_hmtp1: Sed. velocity of hydrometeor at t-lev (k+1) [m/s]
        rho_ds_zt: Dry, static density at thermo. level (k) [kg/m^3]
        rho_ds_ztp1: Dry, static density at thermo. level (k+1) [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. level (k) [m^3/kg]
        invrs_dzmp1: Inverse of grid spacing over m-lev. (k+1) [1/m]
        level: Zero-based thermodynamic level index (batched over levels).
    """

    # Sedimentation is always a downward process, so we omit the upward case.
    # Description:
    # Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    # the "upwind" difference approximation to the vertical derivative.
    #
    # The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    # stands for the sedimentation velocity of the aforementioned hydrometeor.
    #
    # The d(hm)/dt equation contains a sedimentation term:
    #
    # - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    #
    # The variables hm and V_hm in the sedimentation term are divided into mean
    # and turbulent components, and the term is averaged, resulting in:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    # - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    #
    # The mean sedimentation term in the d<hm>/dt equation is:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    #
    # This term is solved for completely implicitly, such that:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    #
    # Note:  When the term is brought over to the left-hand side, the sign is
    #        reversed and the leading "-" in front of the term is changed to
    #        a "+".
    #
    # Timestep index (t) stands for the index of the current timestep, while
    # timestep index (t+1) stands for the index of the next timestep, which is
    # being advanced to in solving the d<hm>/dt equation.
    #
    # This term is discretized as follows when using the upwind-difference
    # approximation:
    #
    # The values of <hm> and the values of V_hmt are found on the thermodynamic
    # levels.  Additionally, the values of rho_ds_zt and the values of
    # invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    # thermodynamic levels, the values of <hm> are multiplied by the values of
    # V_hmt and the values of rho_ds_zt.  Then, the derivative of
    # (rho_ds*<V_hm>*<hm>) is taken between the thermodynamic level above the
    # central thermodynamic level and the central thermodynamic level.  The
    # derivative is multiplied by invrs_rho_ds_zt.
    #
    # --hmp1--V_hmtp1--rho_ds_ztp1--------------------------------------- t(k+1)
    #
    # =================================================================== m(k+1)
    #
    # --hm----V_hmt----rho_ds_zt--invrs_rho_ds_zt--d(rho_ds*V_hm*hm)/dz-- t(k)
    #
    # The vertical indices t(k+1), m(k+1), and t(k) correspond with altitudes
    # zt(k+1), zm(k+1), and zt(k), respectively.  The letter "t" is used for
    # thermodynamic levels and the letter "m" is used for momentum levels.
    #
    # invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) )
    #
    #
    # Conservation Properties:
    #
    # When a hydrometeor is sedimented to the ground (or out the lower boundary
    # of the model), it is removed from the atmosphere (or from the model
    # domain).  Thus, the quantity of the hydrometeor over the entire vertical
    # domain should not be conserved due to the process of sedimentation.  Thus,
    # not all of the column totals in the left-hand side matrix should be equal
    # to 0. Instead, the sum of all the column totals should equal the flux of
    # <hm> out the bottom (zm(1) level, for which the value of the hydrometeor
    # is set equal to the value of the hydrometeor at the zt(1) level) of the
    # domain, -rho_ds_zt(1) * V_hmt(1) * hm(1).  Furthermore, most of the
    # individual column totals should sum to 0, but the 2nd (from the left)
    # column should be equal to the flux out the bottom of the domain.
    #
    # To see that this modified conservation law is satisfied, compute the
    # sedimentation of hm and integrate vertically.  In discretized matrix
    # notation (where "i" stands for the matrix column and "j" stands for the
    # matrix row):
    #
    # - rho_ds_zt(1) * V_hmt(1) * hm(1)
    # = Sum_j Sum_i
    #   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzm )_i
    #   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    #
    # The left-hand side matrix,
    # ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    # partially written below.  The sum over i in the above equation removes
    # invrs_rho_ds_zt and invrs_dzm everywhere from the matrix below.  The sum
    # over j leaves the column totals and the flux at zt(1) that are desired.
    #
    # Left-hand side matrix contributions from the sedimentation term (only);
    # first three vertical levels:
    #
    #     -------------------------------------------------------------------->
    # k=1 | -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)              0
    #    |  *invrs_dzm(k+1)        *invrs_dzm(k+1)
    #    |  *rho_ds_zt(k)          *rho_ds_zt(k+1)
    #    |  *V_hmt(k)              *V_hmt(k+1)
    #    |
    # k=2 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    #    |                         *invrs_dzm(k+1)        *invrs_dzm(k+1)
    #    |                         *rho_ds_zt(k)          *rho_ds_zt(k+1)
    #    |                         *V_hmt(k)              *V_hmt(k+1)
    #    |
    # k=3 |           0                     0             -invrs_rho_ds_zt(k)
    #    |                                                *invrs_dzm(k+1)
    #    |                                                *rho_ds_zt(k)
    #    |                                                *V_hmt(k)
    #   \ /
    #
    # Note:  The superdiagonal term from level 3 is not shown on this diagram.
    # References:
    # None
    # Notes:
    # Both COAMPS Microphysics and Brian Griffin's implementation use
    # Khairoutdinov and Kogan (2000) for the calculation of rain
    # mixing ratio and rain droplet number concentration sedimentation
    # velocities, but COAMPS has only the local parameterization.
    #
    # Please note that "upwind" sedimentation is only 1st-order accurate and
    # highly diffusive.
    # -----------------------------------------------------------------------

    # Downward-only upwind transport: no incoming flux above the top, and
    # no subdiagonal contribution. The fall velocity is negative.
    super = jnp.where(
        level == gr.nzt - 1, 0.0, invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_ztp1 * V_hmtp1
    )
    main = -invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_zt * V_hmt
    return jnp.stack((super, main, jnp.zeros_like(main)))


# -----------------------------------------------------------------------------
def term_turb_sed_lhs(
    gr, Vhmphmp_impcp1, Vhmphmp_impc,    # In
    Vhmphmp_zt_impcp1, Vhmphmp_zt_impc,  # In
    rho_ds_zmp1, rho_ds_zm,              # In
    rho_ds_ztp1, rho_ds_zt,              # In
    invrs_dzt, invrs_dzmp1,              # In
    invrs_rho_ds_zt, level,              # In
):
    """Turbulent sedimentation of a hydrometeor:  implicit portion of the code.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        Vhmphmp_impcp1: Imp. comp. <V_hm'h_m'> interp. m-lev (k+1) [vary]
        Vhmphmp_impc: Imp. comp. <V_hm'h_m'> interp. m-lev (k) [vary]
        Vhmphmp_zt_impcp1: Imp. comp. <V_hm'h_m'>|_zt; t-lev (k+1) [vary]
        Vhmphmp_zt_impc: Imp. comp. <V_hm'h_m'>|_zt; t-lev (k) [vary]
        rho_ds_zmp1: Dry, static density at moment. lev (k+1) [kg/m^3]
        rho_ds_zm: Dry, static density at moment. lev (k) [kg/m^3]
        rho_ds_ztp1: Dry, static density at thermo. level (k+1) [kg/m^3]
        rho_ds_zt: Dry, static density at thermo. level (k) [kg/m^3]
        invrs_dzt: Inverse of grid spacing over t-levs. (k) [1/m]
        invrs_dzmp1: Inverse of grid spacing over m-levs. (k+1) [1/m]
        invrs_rho_ds_zt: Inv dry, static density @ thermo lev (k) [m^3/kg]
        level: Zero-based thermodynamic level index (batched over levels).
    """

    # The implicit turbulent flux has the same discretization as the mean flux.
    # JAX adaptation: share the identical operator, including boundary rows.
    # Description:
    # Turbulent sedimentation of a hydrometeor:  implicit portion of the code.
    #
    # The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    # stands for the sedimentation velocity of the aforementioned hydrometeor.
    #
    # The d(hm)/dt equation contains a sedimentation term:
    #
    # - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    #
    # The variables hm and V_hm in the sedimentation term are divided into mean
    # and turbulent components, and the term is averaged, resulting in:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    # - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    #
    # The turbulent sedimentation term in the d<hm>/dt equation is:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    #
    # This term is solved for semi-implicitly by rewriting < V_hm'hm' > based
    # on < hm > in the manner:
    #
    # < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    #
    # This term can also be solved for completely explicitly (its original
    # form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    # < V_hm'hm' >.  The equation becomes:
    #
    # - (1/rho_ds)
    #   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    #
    # where the timestep index (t+1) means that the value of < hm > being used
    # is from the next timestep, which is being advanced to in solving the
    # d<hm>/dt equation.  Implicit and explicit portions of this term are
    # produced.  The implicit portion of this term is:
    #
    # - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < hm >(t+1) ) / dz.
    #
    # Note:  When the term is brought over to the left-hand side, the sign is
    #        reversed and the leading "-" in front of the d[ ] / dz term is
    #        changed to a "+".
    #
    # This term can be discretized using the centered-difference approximation
    # (which is preferred), or else using the "upwind"-difference approximation.
    #
    # The implicit portion of this term is discretized as follows when using
    # the centered-difference approximation:
    #
    # The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    # thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    # thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    # the momentum levels, and the values of invrs_rho_ds_zt are found on the
    # thermodynamic levels.  The variables < hm > and Vhmphmp_zt_impc are both
    # interpolated to the intermediate momentum levels.  At the momentum levels,
    # the values of interpolated < hm > and interpolated Vhmphmp_zt_impc are
    # multiplied together, and their products are multiplied by the values of
    # rho_ds_zm.  The mathematical expression F is the product of these three
    # variables at momentum levels.  Then, the derivative dF/dz is taken over
    # the central thermodynamic level, where it is multiplied by
    # invrs_rho_ds_zt.  In this function, the value of F is as follows:
    #
    # F = rho_ds_zm * Vhmphmp_impc(interp) * hmm(interp).
    #
    #
    # ----hmmp1--------Vhmphmp_zt_impcp1--------------------------------- t(k+1)
    #
    # =====hmm(interp)=====Vhmphmp_impcp1(interp)=====rho_ds_zmp1======== m(k+1)
    #
    # ----hmm----------Vhmphmp_zt_impc-----invrs_rho_ds_zt-----dF/dz----- t(k)
    #
    # =====hmm(interp)=====Vhmphmp_impc(interp)=======rho_ds_zm========== m(k)
    #
    # ----hmmm1--------Vhmphmp_zt_impcm1--------------------------------- t(k-1)
    #
    # The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    # with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    # The letter "t" is used for thermodynamic levels and the letter "m" is
    # used for momentum levels.
    #
    # invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
    #
    # The implicit portion of this term is discretized as follows when using
    # the upwind-difference approximation:
    #
    # The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    # thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    # thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    # values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    # thermodynamic levels, the values of < hm > and Vhmphmp_zt_impc are
    # multiplied together, and their products are multiplied by the values of
    # rho_ds_zt.  The mathematical expression F is the product of these three
    # variables at thermodynamic levels.  Then, the derivative dF/dz is taken
    # between the thermodynamic level above the central thermodynamic level and
    # the central thermodynamic level.  The derivative is multiplied
    # by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    #
    # F = rho_ds_zt * Vhmphmp_zt_impc * hmm.
    #
    #
    # --hmmp1---Vhmphmp_zt_impcp1---rho_ds_ztp1-------------------------- t(k+1)
    #
    # =================================================================== m(k+1)
    #
    # --hmm-----Vhmphmp_zt_impc-----rho_ds_zt---invrs_rho_ds_zt---dF/dz-- t(k)
    #
    # The vertical indices t(k+1), m(k+1), and t(k) correspond with altitudes
    # zt(k+1), zm(k+1), and zt(k), respectively.  The letter "t" is used for
    # thermodynamic levels and the letter "m" is used for momentum levels.
    #
    # invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) ).
    # References:
    #  None
    #
    # Notes:
    # Please note that "upwind" sedimentation is only 1st-order accurate and
    # highly diffusive.
    # -----------------------------------------------------------------------

    if not parameters_microphys.l_upwind_diff_sed:
        return sed_centered_diff_lhs(
            gr, Vhmphmp_impcp1, Vhmphmp_impc, rho_ds_zmp1,  # In
            rho_ds_zm, invrs_rho_ds_zt,                     # In
            invrs_dzt, level,                               # In
        )
    return sed_upwind_diff_lhs(
        gr, Vhmphmp_zt_impc, Vhmphmp_zt_impcp1, rho_ds_zt,  # In
        rho_ds_ztp1, invrs_rho_ds_zt,                       # In
        invrs_dzmp1, level,                                 # In
    )


# -----------------------------------------------------------------------------
def term_turb_sed_rhs(
    gr, Vhmphmp_expcp1, Vhmphmp_expc,    # In
    Vhmphmp_zt_expcp1, Vhmphmp_zt_expc,  # In
    rho_ds_zmp1, rho_ds_zm,              # In
    rho_ds_ztp1, rho_ds_zt,              # In
    invrs_dzt, invrs_dzmp1,              # In
    invrs_rho_ds_zt, level,              # In
):
    """Turbulent sedimentation of a hydrometeor:  explicit portion of the code.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        Vhmphmp_expcp1: Exp. comp. <V_hm'h_m'> interp. m-lev (k+1) [vary]
        Vhmphmp_expc: Exp. comp. <V_hm'h_m'> interp. m-lev (k) [vary]
        Vhmphmp_zt_expcp1: Exp. comp. <V_hm'h_m'>|_zt; t-lev (k+1) [vary]
        Vhmphmp_zt_expc: Exp. comp. <V_hm'h_m'>|_zt; t-lev (k) [vary]
        rho_ds_zmp1: Dry, static density at moment. lev (k+1) [kg/m^3]
        rho_ds_zm: Dry, static density at moment. lev (k) [kg/m^3]
        rho_ds_ztp1: Dry, static density at thermo. level (k+1) [kg/m^3]
        rho_ds_zt: Dry, static density at thermo. level (k) [kg/m^3]
        invrs_dzt: Inverse of grid spacing over t-levs. (k) [1/m]
        invrs_dzmp1: Inverse of grid spacing over m-levs. (k+1) [1/m]
        invrs_rho_ds_zt: Inv dry, static density @ thermo lev (k) [m^3/kg]
        level: Zero-based thermodynamic level index (batched over levels).
    """

    # Explicit flux divergence; the flux at the model top is zero.
    # Description:
    # Turbulent sedimentation of a hydrometeor:  explicit portion of the code.
    #
    # The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    # stands for the sedimentation velocity of the aforementioned hydrometeor.
    #
    # The d(hm)/dt equation contains a sedimentation term:
    #
    # - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    #
    # The variables hm and V_hm in the sedimentation term are divided into mean
    # and turbulent components, and the term is averaged, resulting in:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    # - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    #
    # The turbulent sedimentation term in the d<hm>/dt equation is:
    #
    # - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    #
    # This term is solved for semi-implicitly by rewriting < V_hm'hm' > based
    # on < hm > in the manner:
    #
    # < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    #
    # This term can also be solved for completely explicitly (its original
    # form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    # < V_hm'hm' >.  The equation becomes:
    #
    # - (1/rho_ds)
    #   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    #
    # where the timestep index (t+1) means that the value of < hm > being used
    # is from the next timestep, which is being advanced to in solving the
    # d<hm>/dt equation.  Implicit and explicit portions of this term are
    # produced.  The explicit portion of this term is:
    #
    # - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    #
    # This term can be discretized using the centered-difference approximation
    # (which is preferred), or else using the "upwind"-difference approximation.
    #
    # The explicit portion of this term is discretized as follows when using
    # the centered-difference approximation:
    #
    # The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    # thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    # thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    # the momentum levels, and the values of invrs_rho_ds_zt are found on the
    # thermodynamic levels.  The variable Vhmphmp_zt_expc is interpolated to the
    # intermediate momentum levels.  At the momentum levels, the values of
    # interpolated Vhmphmp_zt_expc are multiplied by the values of rho_ds_zm.
    # Then, the derivative d(rho_ds*Vhmphmp_zt_expc)/dz is taken over the
    # central thermodynamic level, where it is multiplied by invrs_rho_ds_zt.
    #
    # ---Vhmphmp_zt_expcp1----------------------------------------------- t(k+1)
    #
    # ======Vhmphmp_expcp1(interp)=======rho_ds_zmp1===================== m(k+1)
    #
    # ---Vhmphmp_zt_expc--invrs_rho_ds_zt--d(rho_ds*Vhmphmp_zt_expc)/dz-- t(k)
    #
    # ======Vhmphmp_expc(interp)=========rho_ds_zm======================= m(k)
    #
    # ---Vhmphmp_zt_expcm1----------------------------------------------- t(k-1)
    #
    # The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    # with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    # The letter "t" is used for thermodynamic levels and the letter "m" is
    # used for momentum levels.
    #
    # invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
    #
    # The explicit portion of this term is discretized as follows when using
    # the upwind-difference approximation:
    #
    # The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    # thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    # thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    # values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    # thermodynamic levels, the values of Vhmphmp_zt_expc are multiplied by the
    # values of rho_ds_zt.  The mathematical expression F is the product of
    # these variables at thermodynamic levels.  Then, the derivative dF/dz is
    # taken between the thermodynamic level above the central thermodynamic
    # level and the central thermodynamic level.  The derivative is multiplied
    # by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    #
    # F = rho_ds_zt * Vhmphmp_zt_expc.
    #
    #
    # -----Vhmphmp_zt_expcp1---rho_ds_ztp1------------------------------- t(k+1)
    #
    # =================================================================== m(k+1)
    #
    # -----Vhmphmp_zt_expc-----rho_ds_zt----invrs_rho_ds_zt----dF/dz----- t(k)
    #
    # The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    # zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    # thermodynamic levels and the letter "m" is used for momentum levels.
    #
    # invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) ).
    # References:
    #  None
    #
    # Notes:
    # Please note that "upwind" sedimentation is only 1st-order accurate and
    # highly diffusive.
    # -----------------------------------------------------------------------

    if not parameters_microphys.l_upwind_diff_sed:
        return (
            -invrs_rho_ds_zt
            * invrs_dzt
            * (
                jnp.where(level == gr.nzt - 1, 0.0, rho_ds_zmp1 * Vhmphmp_expcp1)
                - rho_ds_zm * Vhmphmp_expc
            )
        )
    return (
        -invrs_rho_ds_zt
        * invrs_dzmp1
        * (
            jnp.where(level == gr.nzt - 1, 0.0, rho_ds_ztp1 * Vhmphmp_zt_expcp1)
            - rho_ds_zt * Vhmphmp_zt_expc
        )
    )


# -----------------------------------------------------------------------------
def calculate_K_hm(
    gr, ngrdcol, wp2, Kh_zm, Skw_zm, Lscale,  # In
    hydromet_dim, hydromet_tol,               # In
    hydromet, hydrometp2,                     # InOut
    clubb_params,                             # In
    l_use_non_local_diff_fac,                 # In
):
    """Calculate the diffusivity for down-gradient hydrometeor turbulent fluxes.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        wp2: Variance of vertical velocity (momentum levels) [m^2/s^2]
        Kh_zm: Kh Eddy diffusivity on momentum grid [m^2/s]
        Skw_zm: Skewness of w on momentum levels [-]
        Lscale: Length-scale [m]
        hydromet_dim: Number of precipitating hydrometeor fields.
        hydromet_tol: Species thresholds used to distinguish negligible mean hydrometeors
            [units vary].
        hydromet: Hydrometeor mean, <h_m> (thermo. levels) [units]
        hydrometp2: Variance of hydrometeor (overall) (m-levs.) [units^2]
        clubb_params: Column-dependent tunable CLUBB parameters.
        l_use_non_local_diff_fac: Flag to use a non-local factor for eddy-diffusivity applied
            to hydrometeors.
    """

    # Description:
    # The predictive equation for a hydrometeor, hm, contains a turbulent
    # advection term:
    #
    # - (1/rho_ds) * d( rho_ds * <w'hm'> )/dz.
    #
    # The value of <w'hm'> can be calculated in many ways.  For simplicity, a
    # down-gradient approximation is used here, where:
    #
    # <w'hm'> = -K_hm * d<hm>/dz.
    #
    # The coefficient of diffusion, K_hm, is variable and depends on multiple
    # factors.  It optionally includes a non-local turbulent advection factor.
    # References:
    # CLUBB ticket 651 and CLUBB ticket 739.
    # -----------------------------------------------------------------------

    from clubb_jax.src.CLUBB_core.parameter_indices import ic_K_hm, ic_K_hmb, iK_hm_min_coef
    from clubb_jax.src.CLUBB_core.constants_clubb import eps

    K_hm = jnp.zeros_like(hydrometp2)
    for h in range(hydromet_dim):
        hm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., h])
        gradient = ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., h])
        K = (
            clubb_params[:, ic_K_hm, None]
            * Kh_zm
            * (jnp.sqrt(hydrometp2[..., h]) / jnp.maximum(hm_zm, hydromet_tol[h]))
            * (1.0 + jnp.abs(Skw_zm))
        )
        if l_use_non_local_diff_fac:
            K_gamma = 1.0 - clubb_params[:, ic_K_hmb, None] * (
                jnp.maximum(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Lscale), 0.0)
                / jnp.maximum(hm_zm, hydromet_tol[h])
                * gradient
            )
            K = K * jnp.maximum(K_gamma, clubb_params[:, iK_hm_min_coef, None])
        K = jnp.where(
            jnp.abs(gradient) > eps,
            jnp.minimum(
                K,
                (jnp.sqrt(wp2) * jnp.sqrt(hydrometp2[..., h]))
                / jnp.where(jnp.abs(gradient) > eps, jnp.abs(gradient), 1.0),
            ),
            K,
        )
        K_hm = K_hm.at[..., h].set(K.at[:, 0].set(0.0).at[:, -1].set(0.0))
    return K_hm


# -----------------------------------------------------------------------------
def get_cloud_top_level(
    nzt, ngrdcol, rcm, hydromet,  # In
    hydromet_dim, iiri,           # In
):
    """Find the highest level with cloud liquid or ice above its tolerance.

    Return zero-based levels; zero is also the source no-cloud fallback.

    Arguments:
        nzt: Number of thermodynamic levels.
        ngrdcol: Number of grid columns.
        rcm: Mean cloud water mixing ratio [kg/kg]
        hydromet: Hydrometeor mean, <h_m> (thermo. levels) [units vary]
        hydromet_dim: Number of precipitating hydrometeor fields.
        iiri: Zero-based cloud-ice mixing-ratio index; negative if ice is absent.
    """

    # Description:
    # Find cloud top at a given model time step.  This function finds cloud top
    # by looping downward from the top of the model and returning the index of
    # the first vertical level that has a mean cloud water mixing ratio (or a
    # mean cloud ice mixing ratio, when ice is included in the microphysics
    # scheme) that is greater than the tolerance amount.  In a scenario that
    # there is not any cloud found, the function returns a value of 1 (for
    # vertical level 1, which is below the model surface). In JAX this
    # no-cloud fallback is zero, as are all returned level indices.
    # References:
    # -----------------------------------------------------------------------

    rim = hydromet[..., iiri] if iiri >= 0 else jnp.zeros_like(rcm)
    return jnp.max(jnp.where((rcm > rc_tol) | (rim > ri_tol), jnp.arange(nzt)[None, :], 0), axis=-1)


# -----------------------------------------------------------------------------
def write_adv_micro_errors(
    gr, ngrdcol, dt, time_current, hydromet_dim,  # In
    wm_zt, wp2,                                   # In
    exner, rho, rho_zm, rcm,                      # In
    cloud_frac, Kh_zm, Skw_zm,                    # In
    rho_ds_zm, rho_ds_zt,                         # In
    invrs_rho_ds_zt,                              # In
    hydromet_mc, Ncm_mc, Lscale,                  # In
    hydromet_vel_covar_zt_impc,                   # In
    hydromet_vel_covar_zt_expc,                   # In
    clubb_params, nu_vert_res_dep,                # In
    l_upwind_xm_ma,                               # In
    hydromet, hydromet_vel_zt,                    # In
    hydrometp2, K_hm, Ncm,                        # In
    Nc_in_cloud, rvm_mc, thlm_mc,                 # In
    wphydrometp, wpNcp, err_info,                 # In
):
    """Host I/O boundary for the source fatal-error variable dump.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        dt: Model timestep duration [s]
        time_current: Current time [s]
        hydromet_dim: Number of precipitating hydrometeor fields.
        wm_zt: w wind component on thermodynamic levels [m/s]
        wp2: Variance of vertical velocity (momentum levels) [m^2/s^2]
        exner: Exner function [-]
        rho: Density on thermodynamic levels [kg/m^3]
        rho_zm: Density on momentum levels [kg/m^3]
        rcm: Mean cloud water mixing ratio [kg/kg]
        cloud_frac: Cloud fraction [-]
        Kh_zm: Kh Eddy diffusivity on momentum grid [m^2/s]
        Skw_zm: Skewness of w on momentum levels [-]
        rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
        invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
        hydromet_mc: Microphysics tendency for mean hydrometeors [units/s]
        Ncm_mc: Microphysics tendency for Ncm [num/kg/s]
        Lscale: Length-scale [m]
        hydromet_vel_covar_zt_impc: Imp. comp. <V_hm'h_m'> t-levs [m/s]
        hydromet_vel_covar_zt_expc: Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]
        clubb_params: Column-dependent tunable CLUBB parameters.
        nu_vert_res_dep: Resolution-dependent background diffusivities.
        l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
            approximation rather than a centered differencing for turbulent or mean advection
            terms. It affects rtm, thlm, sclrm, um and vm.
        hydromet: Hydrometeor mean, <h_m> (thermo. levels) [units]
        hydromet_vel_zt: Mean hydrometeor sed. vel. on thermo. levs. [m/s]
        hydrometp2: Variance of hydrometeor (overall) (m-levs.) [units^2]
        K_hm: hm eddy diffusivity on momentum grid [m^2/s]
        Ncm: Mean cloud droplet conc., <N_c> (thermo. levs.) [num/kg]
        Nc_in_cloud: Mean (in-cloud) cloud droplet concentration [num/kg]
        rvm_mc: Microphysics contributions to vapor water [kg/kg/s]
        thlm_mc: Microphysics contributions to liquid potential temp. [K/s]
        wphydrometp: Covariance < w'h_m' > (momentum levels) [(m/s)units]
        wpNcp: Covariance < w'N_c' > (momentum levels) [(m/s)(num/kg)]
        err_info: Per-column error state; return any updated fatal status.
    """
    # Description:
    # Writes to screen the values of all variables that are passed into and out
    # of subroutine advance_microphys if a fatal error has been detected.
    # JAX adaptation: the standalone host invokes this after the compiled fatal
    # return, before raising; no subsequent physics/statistics is committed.
    import sys
    import numpy as np  # Host diagnostic formatting only; no physics math.

    if not clubb_at_least_debug_level(0) or not err_info.is_fatal():
        return
    print("Error in advance_microphys", file=sys.stderr)
    # Source field order; device transfer is confined to this host I/O boundary.
    for section, values in (
        (
            "Intent(in)",
            (
                ("dt", dt),
                ("time_current", time_current),
                ("wm_zt", wm_zt),
                ("wp2", wp2),
                ("exner", exner),
                ("rho", rho),
                ("rho_zm", rho_zm),
                ("rcm", rcm),
                ("cloud_frac", cloud_frac),
                ("Kh_zm", Kh_zm),
                ("Skw_zm", Skw_zm),
                ("rho_ds_zm", rho_ds_zm),
                ("rho_ds_zt", rho_ds_zt),
                ("invrs_rho_ds_zt", invrs_rho_ds_zt),
                ("hydromet_mc", hydromet_mc),
                ("Ncm_mc", Ncm_mc),
                ("Lscale", Lscale),
                ("hydromet_vel_covar_zt_impc", hydromet_vel_covar_zt_impc),
                ("hydromet_vel_covar_zt_expc", hydromet_vel_covar_zt_expc),
                ("clubb_params", clubb_params),
                ("nu_hm", nu_vert_res_dep.nu_hm),
                ("l_upwind_xm_ma", l_upwind_xm_ma),
            ),
        ),
        (
            "Intent(inout)",
            (
                ("hydromet", hydromet),
                ("hydromet_vel_zt", hydromet_vel_zt),
                ("hydrometp2", hydrometp2),
                ("K_hm", K_hm),
                ("Ncm", Ncm),
                ("Nc_in_cloud", Nc_in_cloud),
                ("rvm_mc", rvm_mc),
                ("thlm_mc", thlm_mc),
            ),
        ),
        ("Intent(out)", (("wphydrometp", wphydrometp), ("wpNcp", wpNcp))),
    ):
        print(section, file=sys.stderr)
        for name, value in values:
            with np.printoptions(threshold=np.inf):
                print(name, "=", jax.device_get(value), file=sys.stderr)

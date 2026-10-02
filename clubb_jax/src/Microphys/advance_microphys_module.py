"""Hydrometeor transport, mirroring advance_microphys_module.F90.

JAX adaptation: column and level loops are batched; inout state is returned.
"""
import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import zt2zm, zm2zt, ddzt
from clubb_jax.src.CLUBB_core.constants_clubb import Lv, Cp, rho_lw, rc_tol, ri_tol, zero_threshold
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_driver_api, setup_stats_names
from clubb_jax.src.Microphys import parameters_microphys


def advance_microphys(gr, ngrdcol, dt, time_current, hydromet_dim,
    hm_metadata, wm_zt, wp2, exner, rho,
    rho_zm, rcm, cloud_frac, Kh_zm, Skw_zm,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, hydromet_mc, Ncm_mc,
    Lscale, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, clubb_params, nu_vert_res_dep,
    tridiag_solve_method, fill_holes_type, l_upwind_xm_ma, stats, hydromet,
    hydromet_vel_zt, hydrometp2, K_hm, Ncm, Nc_in_cloud,
    rvm_mc, thlm_mc, err_info):
    # Description:
    # Advance mean precipitating hydrometeors and mean cloud droplet
    # concentration one model time step, and calculate some statistics.
    # References:
    #---------------------------------------------------------------------------
    wphydrometp = jnp.zeros((gr.ngrdcol, gr.nzm, hydromet_dim))
    wpNcp = jnp.zeros((gr.ngrdcol, gr.nzm))
    if time_current < parameters_microphys.microphys_start_time:
        return (stats, hydromet, hydromet_vel_zt, hydrometp2, K_hm, Ncm, Nc_in_cloud,
                rvm_mc, thlm_mc, err_info, wphydrometp, wpNcp)
    if hydromet_dim > 0:
        # Solve for the coefficient of diffusion for hydrometeors.
        K_hm = calculate_K_hm(gr, gr.ngrdcol, wp2, Kh_zm, Skw_zm,
                   Lscale, hydromet_dim, hm_metadata.hydromet_tol, hydromet, hydrometp2,
                   clubb_params, False)
        l_prevent_hm_ta_above_cloud = False  # Source local parameter.
        if l_prevent_hm_ta_above_cloud:
            cloud_top_level = get_cloud_top_level(gr.nzt, rcm.shape[0], rcm, hydromet, hydromet_dim,
                hm_metadata.iiri)
            for i in range(hydromet_dim):
                if i not in (hm_metadata.iiri, hm_metadata.iiNi):
                    K_hm = K_hm.at[..., i].set(jnp.where(
                        (jnp.arange(gr.nzm)[None, :] >= cloud_top_level[:, None]) & (cloud_top_level[:, None] > 0),
                        0.0, K_hm[..., i]))
    for i in range(hydromet_dim):
        stats = stats.update('K_hm_' + hm_metadata.hydromet_list[i][:2], K_hm[..., i])
    if parameters_microphys.l_predict_Nc:
        # Solve for K_Nc before advancing the precipitating hydrometeors.
        from clubb_jax.src.CLUBB_core.parameter_indices import ic_K_hm
        K_Nc = clubb_params[:, ic_K_hm, None] * Kh_zm
    if hydromet_dim > 0:
        (stats, hydromet, hydromet_vel_zt, hydrometp2, rvm_mc, thlm_mc,
         err_info, wphydrometp, hydromet_vel, hydromet_vel_covar,
         hydromet_vel_covar_zt) = advance_hydrometeor(gr, gr.ngrdcol, dt, hydromet_dim, hm_metadata,
                wm_zt, exner, cloud_frac, K_hm, rho_ds_zm,
                rho_ds_zt, invrs_rho_ds_zt, hydromet_mc, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc,
                nu_vert_res_dep, l_upwind_xm_ma, tridiag_solve_method, fill_holes_type, stats,
                hydromet, hydromet_vel_zt, hydrometp2, rvm_mc, thlm_mc,
                err_info)
    # JAX adaptation of the source fatal RETURN after advance_hydrometeor.
    # The host emits write_adv_micro_errors after this kernel returns.
    def advance_cloud_number(carry):
        stats, Ncm, Nc_in_cloud, err_info, wpNcp = carry
        if parameters_microphys.l_predict_Nc:
            stats, Ncm, Nc_in_cloud, err_info, wpNcp = advance_Ncm(gr, gr.ngrdcol, dt, wm_zt, cloud_frac,
                K_Nc, rcm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,
                Ncm_mc, nu_vert_res_dep, l_upwind_xm_ma, tridiag_solve_method, stats,
                Ncm, Nc_in_cloud, err_info)
        else:
            # Nc is prescribed; its grid mean depends on cloud fraction.
            Ncm = Nc_in_cloud * cloud_frac
        return stats, Ncm, Nc_in_cloud, err_info, wpNcp

    carry = (stats, Ncm, Nc_in_cloud, err_info, wpNcp)
    if clubb_at_least_debug_level(0):
        stats, Ncm, Nc_in_cloud, err_info, wpNcp = jax.lax.cond(
            err_info.any_fatal(), lambda value: value, advance_cloud_number, carry)
    else:
        stats, Ncm, Nc_in_cloud, err_info, wpNcp = advance_cloud_number(carry)
    # Preserve the second source fatal RETURN, after advance_Ncm.
    def accumulate_statistics(stats):
        stats = stats.update('Ncm', Ncm)
        stats = stats.update('Nc_in_cloud', Nc_in_cloud)
        iirr = hm_metadata.iirr
        if iirr >= 0:
            # Rainfall rate (positive downward) and precipitation flux.
            stats = stats.update('precip_rate_zt', jnp.maximum(-(hydromet[..., iirr] * hydromet_vel_zt[..., iirr] + hydromet_vel_covar_zt[..., iirr]), 0.0) * (rho / rho_lw) * 86400.0 * 1000.0)
            stats = stats.update('Fprec', jnp.maximum(-(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., iirr]) * hydromet_vel[..., iirr] + hydromet_vel_covar[..., iirr]), 0.0) * rho_zm * Lv)
            if parameters_microphys.microphys_scheme != 'morrison':
                stats = stats.update('precip_rate_sfc', jnp.maximum(-(hydromet[:, 0, iirr] * hydromet_vel_zt[:, 0, iirr] + hydromet_vel_covar_zt[:, 0, iirr]), 0.0) * (rho[:, 0] / rho_lw) * 86400.0 * 1000.0)
            stats = stats.update('rain_flux_sfc', jnp.maximum(-(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., iirr])[:, 0] * hydromet_vel[:, 0, iirr] + hydromet_vel_covar[:, 0, iirr]), 0.0) * rho_zm[:, 0] * Lv)
            stats = stats.update('rrm_sfc', zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., iirr])[:, 0])
        from clubb_jax.src.CLUBB_core.stats_clubb_utilities import stats_accumulate_hydromet_api
        stats = stats_accumulate_hydromet_api(gr, ngrdcol, hydromet_dim, hm_metadata, hydromet, rho_ds_zt, stats)
        return stats

    if clubb_at_least_debug_level(0):
        stats = jax.lax.cond(err_info.any_fatal(), lambda value: value,
                             accumulate_statistics, stats)
    else:
        stats = accumulate_statistics(stats)
    # Additional JAX host-safety check: propagate nonfinite hydrometeors as fatal.
    err_info = err_info.set_fatal(mask=jnp.any(~jnp.isfinite(hydromet), axis=(1, 2)))
    return (stats, hydromet, hydromet_vel_zt, hydrometp2, K_hm, Ncm, Nc_in_cloud,
            rvm_mc, thlm_mc, err_info, wphydrometp, wpNcp)


def advance_hydrometeor(gr, ngrdcol, dt, hydromet_dim, hm_metadata,
    wm_zt, exner, cloud_frac, K_hm, rho_ds_zm,
    rho_ds_zt, invrs_rho_ds_zt, hydromet_mc, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc,
    nu_vert_res_dep, l_upwind_xm_ma, tridiag_solve_method, fill_holes_type, stats,
    hydromet, hydromet_vel_zt, hydrometp2, rvm_mc, thlm_mc,
    err_info):
    # Description:
    #   Advance each hydrometeor (precipitating hydrometeor) one model time step.
    # References:
    #   None
    #-----------------------------------------------------------------------
    hydromet_vel = jnp.zeros((gr.ngrdcol, gr.nzm, hydromet_dim))
    ratio_hmp2_on_hmm2 = jnp.zeros_like(hydrometp2)
    wphydrometp = jnp.zeros_like(hydrometp2)
    hydromet_vel_covar = jnp.zeros_like(hydrometp2)
    hydromet_vel_covar_zt = jnp.zeros_like(hydromet)
    for i in range(hydromet_dim):
        max_velocity, name_bt, name_hf, name_wvhf, name_cl, name_mc = setup_stats_names(i, hydromet_dim, hm_metadata.hydromet_list)
        stats = stats.update(name_mc, hydromet_mc[..., i])
        stats = stats.begin_budget(name_bt, hydromet[..., i] / dt)
        hydromet_vel_zt = hydromet_vel_zt.at[..., i].set(jnp.clip(hydromet_vel_zt[..., i], max_velocity, zero_threshold))
        hydromet_vel = hydromet_vel.at[..., i].set(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_zt[..., i]).at[:, -1].set(0.0))
        hydromet_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., i])
        ratio_hmp2_on_hmm2 = ratio_hmp2_on_hmm2.at[..., i].set(jnp.where(
            hydromet_zm > hm_metadata.hydromet_tol[i],
            hydrometp2[..., i] / jnp.where(hydromet_zm > hm_metadata.hydromet_tol[i], hydromet_zm ** 2, 1.0), 0.0))
        # Portion of <w'hm'> using hydromet from timestep t.
        K_hm_nu_hm = K_hm[..., i] + nu_vert_res_dep.nu_hm[:, None]
        xpwp = K_hm_nu_hm * ddzt(gr.nzm, gr.nzt, ngrdcol, gr, hydromet[..., i])
        wphydrometp = wphydrometp.at[..., i].set(
            (-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0))
        # Assemble implicit terms, explicit terms, then solve each species.
        stats, lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, lhs = microphys_lhs(gr, gr.ngrdcol, hm_metadata.hydromet_list[i], parameters_microphys.l_hydromet_sed[i], dt,
            K_hm[..., i], nu_vert_res_dep.nu_hm, wm_zt, hydromet_vel[..., i], hydromet_vel_zt[..., i],
            hydromet_vel_covar_zt_impc[..., i], rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, l_upwind_xm_ma,
            stats)
        stats, rhs = microphys_rhs(gr, gr.ngrdcol, hm_metadata.hydromet_list[i], dt, parameters_microphys.l_hydromet_sed[i],
            hydromet[..., i], hydromet_mc[..., i], K_hm[..., i], nu_vert_res_dep.nu_hm, cloud_frac,
            hydromet_vel_covar_zt_expc[..., i], rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, stats)
        stats, lhs, rhs, hmm, err_info = microphys_solve(gr, gr.ngrdcol, hm_metadata.hydromet_list[i], parameters_microphys.l_hydromet_sed[i], lhs_ta,
            lhs_ma, sed_turb_lhs, sed_diff_lhs, cloud_frac, tridiag_solve_method,
            stats, lhs, rhs, hydromet[..., i], err_info)
        hydromet = hydromet.at[..., i].set(hmm)
    # Now that all species have advanced, fill holes in the profiles.
    stats, thlm_mc, rvm_mc, hydromet = fill_holes_driver_api(
        gr, ngrdcol, gr.nzt, dt, hydromet_dim, hm_metadata, True, rho_ds_zt, exner,
        fill_holes_type, stats, thlm_mc, rvm_mc, hydromet)
    for i in range(hydromet_dim):
        name = hm_metadata.hydromet_list[i]
        hydromet = hydromet.at[:, 0, i].set(jnp.where(hydromet[:, 0, i] < hm_metadata.hydromet_tol[i], zero_threshold, hydromet[:, 0, i]))
        hydromet_zm = jnp.maximum(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., i]), 0.0)
        hydrometp2 = hydrometp2.at[..., i].set(ratio_hmp2_on_hmm2[..., i] * hydromet_zm ** 2)
        # Source overwrites the first half-step flux with this second half.
        xpwp = (K_hm[..., i] + nu_vert_res_dep.nu_hm[:, None]) * ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., i])
        wphydrometp = wphydrometp.at[..., i].set((-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0))
        hydromet_vel_covar_zt = hydromet_vel_covar_zt.at[..., i].set(hydromet_vel_covar_zt_impc[..., i] * hydromet[..., i] + hydromet_vel_covar_zt_expc[..., i])
        hydromet_vel_covar = hydromet_vel_covar.at[..., i].set(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_covar_zt[..., i]).at[:, -1].set(0.0))
        stats = stats.update(name[:2] + 'p2', hydrometp2[..., i])
        stats = stats.update('wp' + name[:2] + 'p', wphydrometp[..., i])
        if name == 'rrm':
            stats = stats.update('Vrrprrp', hydromet_vel_covar[..., i])
        elif name == 'Nrm':
            stats = stats.update('VNrpNrp', hydromet_vel_covar[..., i])
        stats = stats.finalize_budget(name + '_bt', hydromet[..., i] / dt)
    return (stats, hydromet, hydromet_vel_zt, hydrometp2, rvm_mc, thlm_mc, err_info,
            wphydrometp, hydromet_vel, hydromet_vel_covar, hydromet_vel_covar_zt)




def advance_Ncm(gr, ngrdcol, dt, wm_zt, cloud_frac,
    K_Nc, rcm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,
    Ncm_mc, nu_vert_res_dep, l_upwind_xm_ma, tridiag_solve_method, stats,
    Ncm, Nc_in_cloud, err_info):
    # Description:
    # Advance cloud droplet concentration (Ncm) one model time step.
    # References:
    #-----------------------------------------------------------------------
    from clubb_jax.src.CLUBB_core.constants_clubb import cloud_frac_min, pi, mvr_cloud_max, Nc_in_cloud_min
    Ncm_vel_covar_zt_impc = Ncm_vel_covar_zt_expc = jnp.zeros_like(Ncm)
    Ncm_vel = jnp.zeros_like(K_Nc)
    Ncm_vel_zt = jnp.zeros_like(Ncm)
    stats = stats.update('Ncm_mc', Ncm_mc)
    if parameters_microphys.l_in_cloud_Nc_diff:
        stats = stats.begin_budget('Ncm_bt', Nc_in_cloud * jnp.maximum(cloud_frac, cloud_frac_min) / dt)
    else:
        stats = stats.begin_budget('Ncm_bt', Ncm / dt)
    # Solve for <w'Nc'> using the Crank-Nicholson down-gradient approximation.
    xpwp = (K_Nc + nu_vert_res_dep.nu_hm[:, None]) * ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ncm)
    wpNcp = (-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0)
    stats, lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, lhs = microphys_lhs(gr, gr.ngrdcol, 'Ncm', False, dt,
            K_Nc, nu_vert_res_dep.nu_hm, wm_zt, Ncm_vel, Ncm_vel_zt,
            Ncm_vel_covar_zt_impc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, l_upwind_xm_ma,
            stats)
    if parameters_microphys.l_in_cloud_Nc_diff:
        stats, rhs = microphys_rhs(gr, gr.ngrdcol, 'Ncm', dt, False,
            Nc_in_cloud, Ncm_mc / jnp.maximum(cloud_frac, cloud_frac_min), K_Nc, nu_vert_res_dep.nu_hm, cloud_frac,
            Ncm_vel_covar_zt_expc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, stats)
    else:
        stats, rhs = microphys_rhs(gr, gr.ngrdcol, 'Ncm', dt, False,
            Ncm, Ncm_mc, K_Nc, nu_vert_res_dep.nu_hm, cloud_frac,
            Ncm_vel_covar_zt_expc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, stats)
    # Advance Ncm one time step.
    if parameters_microphys.l_in_cloud_Nc_diff:
        stats, lhs, rhs, Nc_in_cloud, err_info = microphys_solve(gr, gr.ngrdcol, 'Ncm', False, lhs_ta,
            lhs_ma, sed_turb_lhs, sed_diff_lhs, cloud_frac, tridiag_solve_method,
            stats, lhs, rhs, Nc_in_cloud, err_info)
        Ncm = Nc_in_cloud * jnp.maximum(cloud_frac, cloud_frac_min)
    else:
        stats, lhs, rhs, Ncm, err_info = microphys_solve(gr, gr.ngrdcol, 'Ncm', False, lhs_ta,
            lhs_ma, sed_turb_lhs, sed_diff_lhs, cloud_frac, tridiag_solve_method,
            stats, lhs, rhs, Ncm, err_info)
        Nc_in_cloud = Ncm / jnp.maximum(cloud_frac, cloud_frac_min)
    # JAX adaptation of the fatal RETURN after microphys_solve: the callable
    # keeps clipping, fluxes, and budget updates entirely on the success path.
    def finish_cloud_number(carry):
        stats, Ncm, Nc_in_cloud, wpNcp = carry
        # Clipping for mean cloud droplet concentration, <Nc>.
        Ncm_mvr_min = (1.0 / ((4.0 / 3.0) * pi * rho_lw * mvr_cloud_max ** 3)) * rcm
        Ncm_min = jnp.maximum(Nc_in_cloud_min * jnp.maximum(cloud_frac, cloud_frac_min), Ncm_mvr_min)
        stats = stats.begin_budget('Ncm_cl', Ncm / dt)
        Ncm = jnp.maximum(Ncm, Ncm_min)
        stats = stats.finalize_budget('Ncm_cl', Ncm / dt)
        Ncic_min = jnp.maximum(Nc_in_cloud_min, Ncm_mvr_min / jnp.maximum(cloud_frac, cloud_frac_min))
        Nc_in_cloud = jnp.maximum(Nc_in_cloud, Ncic_min)
        # Portion of the covariance calculation using <Nc> from timestep t+1.
        xpwp = (K_Nc + nu_vert_res_dep.nu_hm[:, None]) * ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ncm)
        wpNcp = (-0.5 * xpwp).at[:, 0].set(0.0).at[:, -1].set(0.0)
        stats = stats.update('wpNcp', wpNcp)
        stats = stats.finalize_budget('Ncm_bt', Ncm / dt)
        return stats, Ncm, Nc_in_cloud, wpNcp

    carry = (stats, Ncm, Nc_in_cloud, wpNcp)
    if clubb_at_least_debug_level(0):
        stats, Ncm, Nc_in_cloud, wpNcp = jax.lax.cond(
            err_info.any_fatal(), lambda value: value, finish_cloud_number, carry)
    else:
        stats, Ncm, Nc_in_cloud, wpNcp = finish_cloud_number(carry)
    return stats, Ncm, Nc_in_cloud, err_info, wpNcp


def microphys_solve(gr, ngrdcol, solve_type, l_sed, lhs_ta,
    lhs_ma, sed_turb_lhs, sed_diff_lhs, cloud_frac, tridiag_solve_method,
    stats, lhs, rhs, hmm, err_info):
    # Description:
    # Solve the tridiagonal system for hydrometeor variable.
    # References:
    #  None
    #---------------------------------------------------------------------------
    from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve
    from clubb_jax.src.CLUBB_core.constants_clubb import cloud_frac_min
    # Solve system using a tridiag_solve.
    err_info, hmm, _ = tridiag_solve(solve_type, tridiag_solve_method,
                                    gr.ngrdcol, gr.nzt, lhs, rhs, err_info)
    # JAX adaptation: lax.cond needs a callable for the source's fatal RETURN.
    # No implicit statistics may be committed after a failed tridiagonal solve.
    # TODO(port-mirror): low-level source stderr call-stack banners are folded
    # into the standalone host's fatal report. Preserve individual banners when
    # shared device error diagnostics support ordered routine-context messages.
    def accumulate_implicit_statistics(stats):
        # Statistics: implicit contributions to hydrometeor hmm.
        km1 = jnp.maximum(jnp.arange(gr.nzt) - 1, 0)
        kp1 = jnp.minimum(jnp.arange(gr.nzt) + 1, gr.nzt - 1)
        if not stats.l_sample:
            return stats
        if solve_type == 'Ncm' and parameters_microphys.l_in_cloud_Nc_diff:
            ma_term = (-lhs_ma[2] * hmm[:, km1] * jnp.maximum(cloud_frac, cloud_frac_min)
                       - lhs_ma[1] * hmm * jnp.maximum(cloud_frac, cloud_frac_min)
                       - lhs_ma[0] * hmm[:, kp1] * jnp.maximum(cloud_frac, cloud_frac_min))
            ta_term = (-lhs_ta[2] * hmm[:, km1] * jnp.maximum(cloud_frac, cloud_frac_min)
                       - lhs_ta[1] * hmm * jnp.maximum(cloud_frac, cloud_frac_min)
                       - lhs_ta[0] * hmm[:, kp1] * jnp.maximum(cloud_frac, cloud_frac_min))
        else:
            ma_term = -lhs_ma[2] * hmm[:, km1] - lhs_ma[1] * hmm - lhs_ma[0] * hmm[:, kp1]
            ta_term = -lhs_ta[2] * hmm[:, km1] - lhs_ta[1] * hmm - lhs_ta[0] * hmm[:, kp1]
        sd_term = -sed_diff_lhs[2] * hmm[:, km1] - sed_diff_lhs[1] * hmm - sed_diff_lhs[0] * hmm[:, kp1]
        ts_term = -sed_turb_lhs[2] * hmm[:, km1] - sed_turb_lhs[1] * hmm - sed_turb_lhs[0] * hmm[:, kp1]
        if solve_type in ('rrm', 'Nrm', 'rim', 'rsm', 'rgm', 'Ncm', 'Nim', 'Nsm', 'Ngm'):
            stats = stats.update(solve_type + '_ma', ma_term)
            if l_sed and solve_type != 'Ncm':
                stats = stats.update(solve_type + '_sd', sd_term)
            if l_sed and solve_type in ('rrm', 'Nrm'):
                stats = stats.finalize_budget(solve_type + '_ts', ts_term)
            stats = stats.finalize_budget(solve_type + '_ta', ta_term)
        return stats

    if clubb_at_least_debug_level(0):
        stats = jax.lax.cond(err_info.any_fatal(), lambda value: value,
                             accumulate_implicit_statistics, stats)
    else:
        stats = accumulate_implicit_statistics(stats)
    return stats, lhs, rhs, hmm, err_info


def microphys_lhs(gr, ngrdcol, solve_type, l_sed, dt,
    K_hm, nu, wm_zt, V_hm, V_hmt,
    Vhmphmp_zt_impc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, l_upwind_xm_ma,
    stats):
    # Description:
    # Setup the matrix of implicit contributions to a term.
    # Can include the effects of sedimentation, diffusion, and advection.
    # The Morrison microphysics has an explicit sedimentation code, which is
    # handled elsewhere.
    #
    # Notes:
    # Setup for tridiagonal system and boundary conditions should be the same as
    # the original rain subroutine code.
    #-----------------------------------------------------------------------
    from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
    from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs
    Vhmphmp_impc = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Vhmphmp_zt_impc)
    # LHS turbulent advection term. A Crank-Nicholson time-stepping scheme is used.
    Kh_zm = K_hm
    Kh_zt = jnp.maximum(zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, K_hm), 0.0)
    lhs_ta = 0.5 * diffusion_zt_lhs(gr.nzm, gr.nzt, gr.ngrdcol, gr,
                                   Kh_zm, Kh_zt, nu, invrs_rho_ds_zt, rho_ds_zm)
    # The lower boundary condition needs to be applied here at level 1.
    bc = 0.5 * invrs_rho_ds_zt[:, 0] * (gr.invrs_dzt[:, 0] * (Kh_zm[:, 1] + nu) * rho_ds_zm[:, 1] * gr.invrs_dzm[:, 1])
    lhs_ta = lhs_ta.at[0, :, 0].set(-bc).at[1, :, 0].set(bc).at[2, :, 0].set(0.0)
    # LHS mean advection term.
    lhs_ma = term_ma_zt_lhs(gr.nzm, gr.nzt, gr.ngrdcol, wm_zt, gr.weights_zt2zm,
                            gr.invrs_dzt, gr.invrs_dzm, l_upwind_xm_ma, gr.grid_dir)
    # JAX adaptation: level arguments are zero-based vectors, batching the source loop.
    k = jnp.arange(gr.nzt)
    kp1 = jnp.minimum(k + 1, gr.nzt - 1)
    lhs = jnp.zeros_like(lhs_ta).at[1].set(1.0 / dt)
    lhs = lhs + lhs_ta
    lhs = lhs + lhs_ma
    sed_diff_lhs = sed_turb_lhs = jnp.zeros_like(lhs)
    if l_sed:
        # Morrison sedimentation is handled within its core, via the _mc terms.
        if not parameters_microphys.l_upwind_diff_sed:
            sed_diff_lhs = sed_centered_diff_lhs(gr, V_hm[:, kp1], V_hm[:, k],
                rho_ds_zm[:, kp1], rho_ds_zm[:, k], invrs_rho_ds_zt, gr.invrs_dzt, k)
        else:
            sed_diff_lhs = sed_upwind_diff_lhs(gr, V_hmt, V_hmt[:, kp1], rho_ds_zt,
                rho_ds_zt[:, kp1], invrs_rho_ds_zt, gr.invrs_dzm[:, kp1], k)
        lhs = lhs + sed_diff_lhs
        sed_turb_lhs = term_turb_sed_lhs(gr, Vhmphmp_impc[:, kp1], Vhmphmp_impc[:, k],
            Vhmphmp_zt_impc[:, kp1], Vhmphmp_zt_impc, rho_ds_zm[:, kp1], rho_ds_zm[:, k],
            rho_ds_zt[:, kp1], rho_ds_zt, gr.invrs_dzt, gr.invrs_dzm[:, kp1], invrs_rho_ds_zt, k)
        lhs = lhs + sed_turb_lhs
    return stats, lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, lhs


def microphys_rhs(gr, ngrdcol, solve_type, dt, l_sed,
    hmm, hmm_tndcy, K_hm, nu, cloud_frac,
    Vhmphmp_zt_expc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, stats):
    # Description:
    # Compute RHS vector for a given hydrometeor.
    # This subroutine computes the explicit portion of the predictive equation
    # for a given hydrometeor.
    # References:
    #-----------------------------------------------------------------------
    from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
    from clubb_jax.src.CLUBB_core.constants_clubb import cloud_frac_min
    Vhmphmp_expc = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Vhmphmp_zt_expc)
    # LHS turbulent advection term.
    Kh_zm = K_hm
    Kh_zt = jnp.maximum(zm2zt(gr.nzm, gr.nzt, gr.ngrdcol, gr, K_hm), 0.0)
    lhs_ta = 0.5 * diffusion_zt_lhs(gr.nzm, gr.nzt, gr.ngrdcol, gr,
                                   Kh_zm, Kh_zt, nu, invrs_rho_ds_zt, rho_ds_zm)
    # The lower boundary condition needs to be applied here at level 1.
    bc = 0.5 * invrs_rho_ds_zt[:, 0] * (gr.invrs_dzt[:, 0] * (Kh_zm[:, 1] + nu) * rho_ds_zm[:, 1] * gr.invrs_dzm[:, 1])
    lhs_ta = lhs_ta.at[0, :, 0].set(-bc).at[1, :, 0].set(bc).at[2, :, 0].set(0.0)
    k = jnp.arange(gr.nzt)
    km1, kp1 = jnp.maximum(k - 1, 0), jnp.minimum(k + 1, gr.nzt - 1)
    rhs = hmm / dt
    rhs = rhs + hmm_tndcy
    rhs = rhs - lhs_ta[2] * hmm[:, km1] - lhs_ta[1] * hmm - lhs_ta[0] * hmm[:, kp1]
    ts_rhs = jnp.zeros_like(hmm)
    if l_sed:
        ts_rhs = term_turb_sed_rhs(gr, Vhmphmp_expc[:, kp1], Vhmphmp_expc[:, k],
            Vhmphmp_zt_expc[:, kp1], Vhmphmp_zt_expc, rho_ds_zm[:, kp1], rho_ds_zm[:, k],
            rho_ds_zt[:, kp1], rho_ds_zt, gr.invrs_dzt, gr.invrs_dzm[:, kp1], invrs_rho_ds_zt, k)
        rhs = rhs + ts_rhs
    ta_rhs = lhs_ta[2] * hmm[:, km1] + lhs_ta[1] * hmm + lhs_ta[0] * hmm[:, kp1]
    if solve_type == 'Ncm' and parameters_microphys.l_in_cloud_Nc_diff:
        ta_rhs = ta_rhs * jnp.maximum(cloud_frac, cloud_frac_min)
    if solve_type in ('rrm', 'Nrm', 'rim', 'rsm', 'rgm', 'Ncm', 'Nim', 'Nsm', 'Ngm'):
        stats = stats.begin_budget(solve_type + '_ta', ta_rhs)
        if l_sed and solve_type in ('rrm', 'Nrm'):
            stats = stats.begin_budget(solve_type + '_ts', -ts_rhs)
    return stats, rhs


def sed_centered_diff_lhs(gr, V_hmp1, V_hm, rho_ds_zmp1, rho_ds_zm,
                          invrs_rho_ds_zt, invrs_dzt, level):
    # Momentum level (k+1) is between thermodynamic level (k+1) and level (k).
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
    #k=1 |   +invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)            0
    #    |    *invrs_dzt(k)        *invrs_dzt(k)
    #    |    *[ rho_ds_zm(k+1)    *rho_ds_zm(k+1)
    #    |       *V_hm(k+1)*B(k)   *V_hm(k+1)*A(k)
    #    |      -rho_ds_zm(k)
    #    |       *V_hm(k) ]
    #    |
    #k=2 |   -invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    #    |    *invrs_dzt(k)        *invrs_dzt(k)          *invrs_dzt(k)
    #    |    *rho_ds_zm(k)        *[ rho_ds_zm(k+1)      *rho_ds_zm(k+1)
    #    |    *V_hm(k)*D(k)           *V_hm(k+1)*B(k)     *V_hm(k+1)*A(k)
    #    |                           -rho_ds_zm(k)
    #    |                            *V_hm(k)*C(k) ]
    #    |
    #k=3 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    #    |                         *invrs_dzt(k)          *invrs_dzt(k)
    #    |                         *rho_ds_zm(k)          *[ rho_ds_zm(k+1)
    #    |                         *V_hm(k)*D(k)             *V_hm(k+1)*B(k)
    #    |                                                  -rho_ds_zm(k)
    #    |                                                   *V_hm(k)*C(k) ]
    #    |
    #k=4 |           0                     0             -invrs_rho_ds_zt(k)
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
    #-----------------------------------------------------------------------
    mkp1, mk = level + 1, level
    # Lower boundary: surface hydrometeor equals thermodynamic level 1.
    # Upper boundary: no flux through the model top.
    super = jnp.where(level == gr.nzt - 1, 0.0,
        invrs_rho_ds_zt * invrs_dzt * rho_ds_zmp1 * V_hmp1 * gr.weights_zt2zm[:, mkp1, 0])
    main = invrs_rho_ds_zt * invrs_dzt * (
        jnp.where(level == gr.nzt - 1, 0.0, rho_ds_zmp1 * V_hmp1 * gr.weights_zt2zm[:, mkp1, 1])
        - rho_ds_zm * V_hm * jnp.where(level == 0, 1.0, gr.weights_zt2zm[:, mk, 0]))
    sub = jnp.where(level == 0, 0.0,
        -invrs_rho_ds_zt * invrs_dzt * rho_ds_zm * V_hm * gr.weights_zt2zm[:, mk, 1])
    return jnp.stack((super, main, sub))


def sed_upwind_diff_lhs(gr, V_hmt, V_hmtp1, rho_ds_zt, rho_ds_ztp1,
                        invrs_rho_ds_zt, invrs_dzmp1, level):
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
    #k=1 | -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)              0
    #    |  *invrs_dzm(k+1)        *invrs_dzm(k+1)
    #    |  *rho_ds_zt(k)          *rho_ds_zt(k+1)
    #    |  *V_hmt(k)              *V_hmt(k+1)
    #    |
    #k=2 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    #    |                         *invrs_dzm(k+1)        *invrs_dzm(k+1)
    #    |                         *rho_ds_zt(k)          *rho_ds_zt(k+1)
    #    |                         *V_hmt(k)              *V_hmt(k+1)
    #    |
    #k=3 |           0                     0             -invrs_rho_ds_zt(k)
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
    #-----------------------------------------------------------------------
    super = jnp.where(level == gr.nzt - 1, 0.0,
        invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_ztp1 * V_hmtp1)
    main = -invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_zt * V_hmt
    return jnp.stack((super, main, jnp.zeros_like(main)))


def term_turb_sed_lhs(gr, Vhmphmp_impcp1, Vhmphmp_impc, Vhmphmp_zt_impcp1,
                      Vhmphmp_zt_impc, rho_ds_zmp1, rho_ds_zm, rho_ds_ztp1,
                      rho_ds_zt, invrs_dzt, invrs_dzmp1, invrs_rho_ds_zt, level):
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
    # This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    # on < hm > in the manner:
    #
    # < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    #
    # This term can also be solved for completely explicitly (it's original
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
    #-----------------------------------------------------------------------
    if not parameters_microphys.l_upwind_diff_sed:
        return sed_centered_diff_lhs(gr, Vhmphmp_impcp1, Vhmphmp_impc, rho_ds_zmp1,
                                     rho_ds_zm, invrs_rho_ds_zt, invrs_dzt, level)
    return sed_upwind_diff_lhs(gr, Vhmphmp_zt_impc, Vhmphmp_zt_impcp1, rho_ds_zt,
                               rho_ds_ztp1, invrs_rho_ds_zt, invrs_dzmp1, level)


def term_turb_sed_rhs(gr, Vhmphmp_expcp1, Vhmphmp_expc, Vhmphmp_zt_expcp1,
                      Vhmphmp_zt_expc, rho_ds_zmp1, rho_ds_zm, rho_ds_ztp1,
                      rho_ds_zt, invrs_dzt, invrs_dzmp1, invrs_rho_ds_zt, level):
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
    # This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    # on < hm > in the manner:
    #
    # < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    #
    # This term can also be solved for completely explicitly (it's original
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
    #-----------------------------------------------------------------------
    if not parameters_microphys.l_upwind_diff_sed:
        return -invrs_rho_ds_zt * invrs_dzt * (
            jnp.where(level == gr.nzt - 1, 0.0, rho_ds_zmp1 * Vhmphmp_expcp1)
            - rho_ds_zm * Vhmphmp_expc)
    return -invrs_rho_ds_zt * invrs_dzmp1 * (
        jnp.where(level == gr.nzt - 1, 0.0, rho_ds_ztp1 * Vhmphmp_zt_expcp1)
        - rho_ds_zt * Vhmphmp_zt_expc)


def calculate_K_hm(gr, ngrdcol, wp2, Kh_zm, Skw_zm,
    Lscale, hydromet_dim, hydromet_tol, hydromet, hydrometp2,
    clubb_params, l_use_non_local_diff_fac):
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
    #-----------------------------------------------------------------------
    from clubb_jax.src.CLUBB_core.parameter_indices import ic_K_hm, ic_K_hmb, iK_hm_min_coef
    from clubb_jax.src.CLUBB_core.constants_clubb import eps
    K_hm = jnp.zeros_like(hydrometp2)
    for h in range(hydromet_dim):
        hm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., h])
        gradient = ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet[..., h])
        K = (clubb_params[:, ic_K_hm, None] * Kh_zm
             * (jnp.sqrt(hydrometp2[..., h]) / jnp.maximum(hm_zm, hydromet_tol[h]))
             * (1.0 + jnp.abs(Skw_zm)))
        if l_use_non_local_diff_fac:
            K_gamma = 1.0 - clubb_params[:, ic_K_hmb, None] * (
                jnp.maximum(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Lscale), 0.0)
                / jnp.maximum(hm_zm, hydromet_tol[h]) * gradient)
            K = K * jnp.maximum(K_gamma, clubb_params[:, iK_hm_min_coef, None])
        K = jnp.where(jnp.abs(gradient) > eps,
                      jnp.minimum(K, (jnp.sqrt(wp2) * jnp.sqrt(hydrometp2[..., h]))
                                  / jnp.where(jnp.abs(gradient) > eps, jnp.abs(gradient), 1.0)), K)
        K_hm = K_hm.at[..., h].set(K.at[:, 0].set(0.0).at[:, -1].set(0.0))
    return K_hm


def get_cloud_top_level(nzt, ngrdcol, rcm, hydromet, hydromet_dim,
    iiri):
    # Description:
    # Find cloud top at a given model time step.  This function finds cloud top
    # by looping downward from the top of the model and returning the index of
    # the first vertical level that has a mean cloud water mixing ratio (or a
    # mean cloud ice mixing ratio, when ice is included in the microphysics
    # scheme) that is greater than the tolerance amount.  In a scenario that
    # there is not any cloud found, the function returns a value of 1 (for
    # vertical level 1, which is below the model surface).
    # References:
    #-----------------------------------------------------------------------
    rim = hydromet[..., iiri] if iiri >= 0 else jnp.zeros_like(rcm)
    return jnp.max(jnp.where((rcm > rc_tol) | (rim > ri_tol), jnp.arange(nzt)[None, :], 0), axis=-1)


def write_adv_micro_errors(gr, ngrdcol, dt, time_current, hydromet_dim,
    wm_zt, wp2, exner, rho, rho_zm,
    rcm, cloud_frac, Kh_zm, Skw_zm, rho_ds_zm,
    rho_ds_zt, invrs_rho_ds_zt, hydromet_mc, Ncm_mc, Lscale,
    hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, clubb_params, nu_vert_res_dep, l_upwind_xm_ma,
    hydromet, hydromet_vel_zt, hydrometp2, K_hm, Ncm,
    Nc_in_cloud, rvm_mc, thlm_mc, wphydrometp, wpNcp,
    err_info):
    """Host I/O boundary for the source fatal-error variable dump."""
    # Description:
    # Writes to screen the values of all variables that are passed into and out
    # of subroutine advance_microphys if a fatal error has been detected.
    # JAX adaptation: the standalone host invokes this after the compiled fatal
    # return, before raising; no subsequent physics/statistics is committed.
    import sys
    import numpy as np  # Host diagnostic formatting only; no physics math.
    if not clubb_at_least_debug_level(0) or not err_info.is_fatal():
        return
    print('Error in advance_microphys', file=sys.stderr)
    # Source field order; device transfer is confined to this host I/O boundary.
    for section, values in (
        ('Intent(in)', (
            ('dt', dt), ('time_current', time_current), ('wm_zt', wm_zt),
            ('wp2', wp2), ('exner', exner), ('rho', rho), ('rho_zm', rho_zm),
            ('rcm', rcm), ('cloud_frac', cloud_frac), ('Kh_zm', Kh_zm),
            ('Skw_zm', Skw_zm), ('rho_ds_zm', rho_ds_zm), ('rho_ds_zt', rho_ds_zt),
            ('invrs_rho_ds_zt', invrs_rho_ds_zt), ('hydromet_mc', hydromet_mc),
            ('Ncm_mc', Ncm_mc), ('Lscale', Lscale),
            ('hydromet_vel_covar_zt_impc', hydromet_vel_covar_zt_impc),
            ('hydromet_vel_covar_zt_expc', hydromet_vel_covar_zt_expc),
            ('clubb_params', clubb_params), ('nu_hm', nu_vert_res_dep.nu_hm),
            ('l_upwind_xm_ma', l_upwind_xm_ma))),
        ('Intent(inout)', (
            ('hydromet', hydromet), ('hydromet_vel_zt', hydromet_vel_zt),
            ('hydrometp2', hydrometp2), ('K_hm', K_hm), ('Ncm', Ncm),
            ('Nc_in_cloud', Nc_in_cloud), ('rvm_mc', rvm_mc), ('thlm_mc', thlm_mc))),
        ('Intent(out)', (('wphydrometp', wphydrometp), ('wpNcp', wpNcp))),
    ):
        print(section, file=sys.stderr)
        for name, value in values:
            with np.printoptions(threshold=np.inf):
                print(name, '=', jax.device_get(value), file=sys.stderr)

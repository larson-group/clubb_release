"""Scheme tendency dispatch from microphys_driver.F90.

JAX adaptation: columns are batched, Fortran inout/out values are returned.
"""
import jax.numpy as jnp
from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.src.Microphys.KK_microphys_module import KK_local_microphys_driver, KK_upscaled_microphys
from clubb_jax.src.Microphys.cloud_sed_module import cloud_drop_sed
from clubb_jax.src.CLUBB_core.grid_class import zt2zm, zm2zt, zm2zt2zm
from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
from clubb_jax.src.CLUBB_core.constants_clubb import w_tol, w_tol_sqd


def calc_microphys_scheme_tendcies(gr, ngrdcol, dt, time_current, pdf_dim, hydromet_dim, runtype,
                                  thlm, p_in_Pa, exner, rho, rho_zm, rtm,
                                  rcm, cloud_frac, wm_zt, wm_zm, wp2, wp3, clubb_params,
                                  hydromet, Nc_in_cloud, hm_metadata,
                                  pdf_params, hydromet_pdf_params, precip_fracs,
                                  X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights,
                                  mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,
                                  corr_array_1_n, corr_array_2_n,
                                  lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped,
                                  l_lh_importance_sampling, l_lh_instant_var_covar_src,
                                  saturation_formula, stats, Nccnm):
    # We must initialize intent(out) variables. If they do not otherwise get
    # assigned, they must still have defined values.
    # Description:
    # Call a microphysics scheme and output microphysics tendencies for the
    # predictive variables.
    # References:
    # H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
    # moment microphysics scheme for application in cloud and
    # climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.
    #
    # Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
    # large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
    # 229-243, 2000.
    #-----------------------------------------------------------------------
    unit_sample_weight_2d = jnp.ones_like(rcm)
    hydromet_mc = jnp.zeros_like(hydromet)
    Ncm_mc = rcm_mc = rvm_mc = thlm_mc = jnp.zeros_like(rcm)
    hydromet_vel_zt = hydromet_vel_covar_zt_impc = hydromet_vel_covar_zt_expc = jnp.zeros_like(hydromet)
    wprtp_mc = wpthlp_mc = rtp2_mc = thlp2_mc = rtpthlp_mc = jnp.zeros((gr.ngrdcol, gr.nzm))
    # Calculate Skw_zm for use in advance_microphys.
    wp3_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, wp3)
    Skw_zm = Skx_func(gr.nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params)
    # Smooth by interpolating to thermodynamic levels and back.
    Skw_zm_smooth = zm2zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, Skw_zm)
    wp2_zt = zm2zt(gr.nzm, gr.nzt, ngrdcol, gr, wp2, w_tol_sqd)
    # Return if there is delay between model start and microphysics start.
    if time_current < parameters.microphys_start_time:
        return (stats, Nccnm, hydromet_mc, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
                hydromet_vel_zt, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc,
                wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc, Skw_zm_smooth)
    if not runtype:
        raise ValueError('Runtype is null, which should not happen')
    # Calculate the updated mean cloud droplet concentration from Nc_in_cloud
    # and the updated cloud fraction.
    Ncm_microphys = Nc_in_cloud * cloud_frac
    # Determine 's' from Mellor (1977).
    chi = pdf_params.mixt_frac * pdf_params.chi_1 + (1.0 - pdf_params.mixt_frac) * pdf_params.chi_2
    # Compute standard deviation of vertical velocity in the grid column.
    wtmp = jnp.sqrt(wp2_zt)
    # Morrison delta_zt(k) is zt(k+1)-zt(k), which is CLUBB dzm(k+1).
    delta_zt = gr.dzm[:, 1:]
    # TODO(port-mirror): SILHS sampling and COAMPS dispatch remain disabled;
    # their source-signature dummy modules fail explicitly until cores exist.
    # GFDL-only aerosol mass and temperature preparation is likewise disabled.
    if parameters.lh_microphys_type != parameters.lh_microphys_disabled:
        raise NotImplementedError('SILHS microphysics is disabled')
    if parameters.microphys_scheme == 'morrison':
        from clubb_jax.src.Microphys.morrison_microphys_module import morrison_microphys_driver
        (stats, hydromet_mc, hydromet_vel_zt, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
         rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag) = morrison_microphys_driver(
            gr, ngrdcol, dt, gr.nzt, hydromet_dim, hm_metadata, False, thlm, wm_zt, p_in_Pa,
            exner, rho, cloud_frac, wtmp, delta_zt, rcm, Ncm_microphys, chi,
            rtm - rcm, hydromet, saturation_formula, unit_sample_weight_2d, stats)
        # Output rain sedimentation velocity.
        stats = stats.update('Vrr', zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_zt[..., hm_metadata.iirr]))
    elif parameters.microphys_scheme == 'khairoutdinov_kogan':
        if parameters.l_local_kk:
            (stats, hydromet_mc, hydromet_vel_zt, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
             rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag) = KK_local_microphys_driver(
                gr, ngrdcol, dt, gr.nzt, hydromet_dim, hm_metadata, False, thlm, wm_zt,
                p_in_Pa, exner, rho, cloud_frac, wtmp, delta_zt, rcm, Ncm_microphys,
                chi, rtm - rcm, hydromet, saturation_formula, unit_sample_weight_2d, stats)
        else:
            (stats, hydromet_mc, hydromet_vel_zt, rcm_mc, rvm_mc, thlm_mc,
             hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc,
             wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc) = KK_upscaled_microphys(
                gr, ngrdcol, dt, gr.nzt, gr.nzm, pdf_dim, hydromet_dim, hm_metadata,
                wm_zt, rtm, thlm, p_in_Pa, exner, rho, rcm, pdf_params, hydromet_pdf_params,
                precip_fracs, hydromet, mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,
                corr_array_1_n, corr_array_2_n, saturation_formula, stats)
            if parameters.l_silhs_KK_convergence_adj_mean:
                hydromet_vel_covar_zt_impc = jnp.zeros_like(hydromet)
                hydromet_vel_covar_zt_expc = jnp.zeros_like(hydromet)
        stats = stats.update('Vrr', zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_zt[..., hm_metadata.iirr]))
        stats = stats.update('VNr', zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet_vel_zt[..., hm_metadata.iiNr]))
    elif parameters.microphys_scheme != 'none':
        raise ValueError(f'Unsupported microphysics scheme: {parameters.microphys_scheme}')
    stats = stats.update('Nccnm', Nccnm)
    if parameters.l_gfdl_activation:
        raise NotImplementedError('GFDL activation core is disabled')
    # Cloud water sedimentation.
    if parameters.l_cloud_sed:
        stats, rcm_mc, thlm_mc = cloud_drop_sed(gr, gr.ngrdcol, rcm, Ncm_microphys, rho_zm,
            rho, exner, parameters.sigma_g, stats, rcm_mc,
            thlm_mc)
    return (stats, Nccnm, hydromet_mc, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
            hydromet_vel_zt, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc,
            wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc, Skw_zm_smooth)

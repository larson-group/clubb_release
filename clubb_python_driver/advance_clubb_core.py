"""Python port of src/CLUBB_core/advance_clubb_core_module.F90.

Translates the Fortran advance_clubb_core subroutine into Python/NumPy,
calling individual Fortran subroutines via the F2PY API for the complex
advance/closure routines, and using NumPy for the inline array math.
"""

import numpy as np

from clubb_python import clubb_api
from clubb_python_driver.clubb_constants import (
    em_min,
    three_halves,
    two,
    zero,
    zero_threshold,
    # Parameter indices
    ia_const,
    ibv_efold,
    ic_K,
    iC_wp2_splat,
    ilambda0_stability_coef,
    iup2_sfc_coef,
    # Model flag constants
    iiPDF_ADG1,
    ipdf_post_advance_fields,
    ipdf_pre_advance_fields,
    ipdf_pre_post_advance_fields,
    l_advance_xp3,
    l_use_invrs_tau_N2_iso,
    order_windm,
    order_wp2_wp3,
    order_xm_wpxp,
    order_xp2_xpyp,
    tau_const,
    ufmin,
)

CLUBB_FATAL_ERROR = 99

def advance_clubb_core(
    gr, nzm, nzt, ngrdcol,
    l_implemented, dt, fcor, fcor_y, sfc_elevation,
    hydromet_dim,
    sclr_dim, sclr_tol, edsclr_dim, sclr_idx,
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing,
    sclrm_forcing, edsclrm_forcing, wprtp_forcing,
    wpthlp_forcing, rtp2_forcing, thlp2_forcing,
    rtpthlp_forcing, wm_zm, wm_zt,
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc,
    wpsclrp_sfc, wpedsclrp_sfc,
    upwp_sfc_pert, vpwp_sfc_pert,
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg,
    p_in_Pa, rho_zm, rho, exner,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
    invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt,
    l_mix_rat_hm,
    rfrzm,
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt,
    host_dx, host_dy,
    clubb_params, nu_vert_res_dep, lmin,
    mixt_frac_max_mag, t0, ts_nudge,
    rtm_min, rtm_nudge_max_altitude,
    clubb_config_flags,
    l_sample,
    um, vm, upwp, vpwp, up2, vp2, up3, vp3,
    thlm, rtm, wprtp, wpthlp,
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp,
    sclrm,
    sclrp2, sclrp3, sclrprtp, sclrpthlp,
    wpsclrp, edsclrm,
    rcm, cloud_frac,
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
    sclrpthvp,
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4,
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac,
    um_pert, vm_pert, upwp_pert, vpwp_pert,
    pdf_params, pdf_params_zm,
    pdf_implicit_coefs_terms,
    err_info,
):
    """Advance CLUBB one timestep with an explicit argument surface."""
    Kh_zm = np.zeros((ngrdcol, nzm))
    Kh_zt = np.zeros((ngrdcol, nzt))
    Lscale = np.zeros((ngrdcol, nzt))
    invrs_tau_zm = np.zeros((ngrdcol, nzm))
    cloud_frac_zm = np.zeros((ngrdcol, nzm))
    ice_supersat_frac_zm = np.zeros((ngrdcol, nzm))
    rc_coef = np.zeros((ngrdcol, nzt))
    rcm_supersat_adj = np.zeros((ngrdcol, nzt))
    rcm_zm = np.zeros((ngrdcol, nzm))
    rcp2 = np.zeros((ngrdcol, nzm))
    rcp2_zt = np.zeros((ngrdcol, nzt))
    rtm_zm = np.zeros((ngrdcol, nzm))
    rtprcp = np.zeros((ngrdcol, nzm))
    sclrprcp = np.zeros((ngrdcol, nzm, max(sclr_dim, 1)))
    skw_velocity = np.zeros((ngrdcol, nzm))
    thlm_zm = np.zeros((ngrdcol, nzm))
    wp2rcp = np.zeros((ngrdcol, nzt))
    wp2sclrp = np.zeros((ngrdcol, nzt, max(sclr_dim, 1)))
    wprtp2 = np.zeros((ngrdcol, nzt))
    wprtpthlp = np.zeros((ngrdcol, nzt))
    wpsclrp2 = np.zeros((ngrdcol, nzt, max(sclr_dim, 1)))
    wpsclrprtp = np.zeros((ngrdcol, nzt, max(sclr_dim, 1)))
    wpsclrpthlp = np.zeros((ngrdcol, nzt, max(sclr_dim, 1)))
    wpthlp2 = np.zeros((ngrdcol, nzt))
    cloud_cover = np.zeros((ngrdcol, nzt))
    cloudy_downdraft_frac = np.zeros((ngrdcol, nzt))
    cloudy_updraft_frac = np.zeros((ngrdcol, nzt))
    rcm_in_layer = np.zeros((ngrdcol, nzt))
    thlprcp = np.zeros((ngrdcol, nzm))
    w_down_in_cloud = np.zeros((ngrdcol, nzt))
    w_up_in_cloud = np.zeros((ngrdcol, nzt))
    wprcp_out = np.zeros((ngrdcol, nzm))
    rsat = np.zeros((ngrdcol, nzt))
    l_sample = bool(l_sample)

    dt_advance = two * dt if clubb_config_flags.l_lmm_stepping else dt

    #----------------------------------------------------------------
    # Test input variables
    #----------------------------------------------------------------
    if clubb_api.clubb_at_least_debug_level(2):

        err_info = clubb_api.parameterization_check(
            nzm, nzt, ngrdcol, sclr_dim, edsclr_dim,
            thlm_forcing, rtm_forcing, um_forcing,
            vm_forcing, wm_zm, wm_zt, p_in_Pa,
            rho_zm, rho, exner, rho_ds_zm,
            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
            thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,
            vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,
            rtm, wprtp, thlm, wpthlp, wp2, wp3,
            rtp2, thlp2, rtpthlp,
            "beginning of ",
            wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,
            sclrp2,
            sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,
            edsclrm_forcing,
            err_info,
        )

        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

    if l_sample:
        rtm_before = rtm.copy()
        thlm_before = thlm.copy()

        clubb_api.stats_update("rfrzm", rfrzm)

        # Set up budget stats variables.
        clubb_api.stats_begin_budget("wp2_bt", wp2 / dt)
        clubb_api.stats_begin_budget("vp2_bt", vp2 / dt)
        clubb_api.stats_begin_budget("up2_bt", up2 / dt)
        clubb_api.stats_begin_budget("wprtp_bt", wprtp / dt)
        clubb_api.stats_begin_budget("wpthlp_bt", wpthlp / dt)
        if clubb_config_flags.l_predict_upwp_vpwp:
            clubb_api.stats_begin_budget("upwp_bt", upwp / dt)
            clubb_api.stats_begin_budget("vpwp_bt", vpwp / dt)
        clubb_api.stats_begin_budget("rtp2_bt", rtp2 / dt)
        clubb_api.stats_begin_budget("thlp2_bt", thlp2 / dt)
        clubb_api.stats_begin_budget("rtpthlp_bt", rtpthlp / dt)
        clubb_api.stats_begin_budget("rtm_bt", rtm / dt)
        clubb_api.stats_begin_budget("thlm_bt", thlm / dt)
        clubb_api.stats_begin_budget("um_bt", um / dt)
        clubb_api.stats_begin_budget("vm_bt", vm / dt)
        clubb_api.stats_begin_budget("wp3_bt", wp3 / dt)

    wprtp_cl_num = 0
    wpthlp_cl_num = 0
    upwp_cl_num = 0
    vpwp_cl_num = 0
        
    wpedsclrp = np.zeros((ngrdcol, nzm, max(edsclr_dim, 1)))
    (wpthlp, wprtp, upwp, vpwp,
     upwp_pert, vpwp_pert,
     wpsclrp, wpedsclrp) = set_sfc_value_of_flux_profiles(
        nzm, ngrdcol, sclr_dim, edsclr_dim, gr,
        clubb_config_flags.l_host_applies_sfc_fluxes, clubb_config_flags.l_linearize_pbl_winds,
        wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
        upwp_sfc_pert, vpwp_sfc_pert,
        wpsclrp_sfc, wpedsclrp_sfc,
        wpthlp, wprtp, upwp, vpwp,
        upwp_pert, vpwp_pert,
        wpsclrp, wpedsclrp,
    )

    sigma_sqd_w = clubb_api.compute_sigma_sqd_w(
        nzm, nzt, ngrdcol, gr,
        wp3, wp2, thlp2, rtp2,
        up2, vp2, wpthlp, wprtp, upwp, vpwp,
        clubb_params,
        clubb_config_flags.l_predict_upwp_vpwp,
    )

    if (clubb_config_flags.ipdf_call_placement == ipdf_pre_advance_fields
            or clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields):

        # Sample stats in this call to subroutine pdf_closure_driver for
        # both of these options (ipdf_pre_advance_fields and
        # ipdf_pre_post_advance_fields).
        if clubb_config_flags.ipdf_call_placement == ipdf_pre_advance_fields:
            l_samp_stats = True
        elif clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields:
            l_samp_stats = True

        #########################################################################
        ########                     CALL CLUBB's PDF                     #######
        ########   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
        #########################################################################
        pdf_result = clubb_api.pdf_closure_driver(
            gr, nzm, nzt, ngrdcol,
            dt, hydromet_dim, sclr_dim, sclr_tol,
            wprtp, thlm, wpthlp, rtp2, rtp3,
            thlp2, thlp3, rtpthlp, wp2,
            wp3, wm_zm, wm_zt,
            um, up2, upwp, up3,
            vm, vp2, vpwp, vp3,
            p_in_Pa, exner,
            thv_ds_zm, thv_ds_zt, rtm_ref,
            wphydrometp,
            wp2hmp, rtphmp_zt, thlphmp_zt,
            sclrm, wpsclrp, sclrp2,
            sclrprtp, sclrpthlp, sclrp3,
            p_sfc, l_samp_stats,
            mixt_frac_max_mag, ts_nudge,
            rtm_min, rtm_nudge_max_altitude,
            clubb_params,
            clubb_config_flags.iiPDF_type,
            clubb_config_flags.saturation_formula,
            clubb_config_flags.l_rtm_nudge,
            clubb_config_flags.l_trapezoidal_rule_zt,
            clubb_config_flags.l_trapezoidal_rule_zm,
            clubb_config_flags.l_call_pdf_closure_twice,
            clubb_config_flags.l_use_cloud_cover,
            clubb_config_flags.l_rcm_supersat_adj,
            l_mix_rat_hm,
            rtm, sigma_sqd_w,
            pdf_implicit_coefs_terms,
            pdf_params, pdf_params_zm, err_info,
        )
        (rtm, pdf_implicit_coefs_terms, pdf_params, pdf_params_zm, err_info,
         rcm, cloud_frac, ice_supersat_frac, wprcp_out, sigma_sqd_w, wpthvp, wp2thvp,
         wp2up, rtpthvp, thlpthvp, rc_coef, rcm_in_layer, cloud_cover,
         rcp2_zt, thlprcp, rc_coef_zm, sclrpthvp, wpup2, wpvp2, wp2up2,
         wp2vp2, wp4, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, wp2rcp,
         rtprcp, rcp2, uprcp, vprcp, w_up_in_cloud, w_down_in_cloud,
         cloudy_updraft_frac, cloudy_downdraft_frac, skw_velocity,
         cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, thlm_zm, rcm_zm,
         rcm_supersat_adj, wp2sclrp, wpsclrp2, sclrprcp,
         wpsclrprtp, wpsclrpthlp) = pdf_result

        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

    # This feels like an awkward place to have this, but this keeps things BFB
    # with where it was before
    if l_sample:
        if clubb_config_flags.iiPDF_type == iiPDF_ADG1:
            sigma_sqd_w_zt = clubb_api.zm2zt(
                nzm, nzt, ngrdcol, gr, sigma_sqd_w,
            )
            sigma_sqd_w_zt = np.maximum(sigma_sqd_w_zt, zero_threshold)
        else:
            sigma_sqd_w_zt = np.zeros((ngrdcol, nzt))
        clubb_api.stats_update("sigma_sqd_w_zt", sigma_sqd_w_zt)

    thvm = clubb_api.calculate_thvm(
        nzt, ngrdcol,
        thlm, rtm, rcm, exner, thv_ds_zt,
    )

    # Compute tke (turbulent kinetic energy).
    if not clubb_config_flags.l_tke_aniso:
        # tke is assumed to be 3/2 of wp2.
        em = three_halves * wp2
    else:
        em = 0.5 * (wp2 + vp2 + up2)

    if clubb_config_flags.l_call_pdf_closure_twice:
        w_1_zm = pdf_params_zm.w_1.copy()
        w_2_zm = pdf_params_zm.w_2.copy()
        varnce_w_1_zm = pdf_params_zm.varnce_w_1.copy()
        varnce_w_2_zm = pdf_params_zm.varnce_w_2.copy()
        mixt_frac_zm = pdf_params_zm.mixt_frac.copy()
    else:
        w_1_zm = clubb_api.zt2zm(nzm, nzt, ngrdcol, gr, pdf_params.w_1)
        w_2_zm = clubb_api.zt2zm(nzm, nzt, ngrdcol, gr, pdf_params.w_2)
        varnce_w_1_zm = clubb_api.zt2zm(
            nzm, nzt, ngrdcol, gr, pdf_params.varnce_w_1,
        )
        varnce_w_2_zm = clubb_api.zt2zm(
            nzm, nzt, ngrdcol, gr, pdf_params.varnce_w_2,
        )
        mixt_frac_zm = clubb_api.zt2zm(
            nzm, nzt, ngrdcol, gr, pdf_params.mixt_frac,
        )

    if l_sample:
        clubb_api.stats_update("rvm", rtm - rcm)
        # Output relative humidity (q/q* where q* is the saturation mixing ratio
        # over liquid).
        # Added an extra check for rel_humidity or rsat output; otherwise,
        # if both are omitted, rsat may not be computed, leading to a
        # floating-point exception when rel_humidity is written.
        if clubb_api.var_on_stats_list("rel_humidity") or clubb_api.var_on_stats_list("rsat"):
            T_in_K = np.empty((ngrdcol, nzt))
            for i in range(ngrdcol):
                T_in_K[i, :] = clubb_api.thlm2t_in_k(
                    nzt, thlm[i, :], exner[i, :], rcm[i, :]
                )
            rsat = clubb_api.sat_mixrat_liq(
                nzt, ngrdcol, gr, p_in_Pa,
                T_in_K,
                clubb_config_flags.saturation_formula,
            )
            # Recompute rsat and rel_humidity. They might have changed.
            rel_humidity = (rtm - rcm) / rsat
            # Write the var to stats
            clubb_api.stats_update("rel_humidity", rel_humidity)

    bv_result = clubb_api.calc_brunt_vaisala_freq_sqd(
        nzm, nzt, ngrdcol, gr, thlm,
        exner, rtm, rcm, p_in_Pa, thvm,
        ice_supersat_frac,
        clubb_config_flags.saturation_formula,
        clubb_config_flags.l_brunt_vaisala_freq_moist,
        clubb_config_flags.l_use_thvm_in_bv_freq,
        clubb_config_flags.l_modify_limiters_for_cnvg_test,
        clubb_params[:, ibv_efold - 1], t0,
    )
    (brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed,
     brunt_vaisala_freq_sqd_smth) = bv_result

    #----------------------------------------------------------------
    # Compute mixing length and dissipation time
    #----------------------------------------------------------------
    ddzt_umvm_sqd = clubb_api.calc_ddzt_umvm_sqd(
        nzm, nzt, ngrdcol, gr, um, vm,
    )

    lscale_result = clubb_api.calc_lscale(
        nzm, nzt, ngrdcol, gr, l_implemented, host_dx, host_dy,
        p_in_Pa, exner, rtm, thlm, thvm,
        thlp2, rtp2, rtpthlp,
        pdf_params, em, thv_ds_zt, lmin,
        upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac,
        ufmin, tau_const, sfc_elevation, clubb_params,
        clubb_config_flags.saturation_formula,
        clubb_config_flags.l_Lscale_plume_centered,
        clubb_config_flags.l_diag_Lscale_from_tau,
        clubb_config_flags.l_e3sm_config,
        clubb_config_flags.l_smooth_Heaviside_tau_wpxp,
        clubb_config_flags.l_modify_limiters_for_cnvg_test,
        l_use_invrs_tau_N2_iso,
        brunt_vaisala_freq_sqd_smth,
        err_info,
    )
    (err_info,
     invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm, invrs_tau_wp3_zt,
     invrs_tau_C1_zm, invrs_tau_C4_zm, invrs_tau_C6_zm, invrs_tau_C14_zm,
     tau_max_zm, tau_max_zt, tau_zm,
     Lscale, Lscale_zm, Lscale_up, Lscale_down) = lscale_result
    err_code = err_info.err_code
    if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
        return

    if clubb_config_flags.l_stability_correct_tau_zm:
        stability_correction = clubb_api.calc_stability_correction(
            nzm, ngrdcol, brunt_vaisala_freq_sqd,
            Lscale_zm, em,
            clubb_params[:, ilambda0_stability_coef - 1],
        )
        invrs_tau_C6_zm = invrs_tau_zm * stability_correction
        invrs_tau_C1_zm = invrs_tau_C6_zm.copy()
    else:
        stability_correction = np.zeros((ngrdcol, nzm))

    #----------------------------------------------------------------
    # Eddy diffusivity coefficient
    #----------------------------------------------------------------
    # c_K is 0.548 usually (Duynkerke and Driedonks 1987)
    # CLUBB uses a smaller value to better fit empirical data.
    #
    # Calculate CLUBB's eddy diffusivity as
    #   CLUBB's length scale times a velocity scale.
    em_zt = np.maximum(
        clubb_api.zm2zt(nzm, nzt, ngrdcol, gr, em),
        em_min,
    )
    c_K = clubb_params[:, ic_K - 1]
    Kh_zt = c_K[:, None] * Lscale * np.sqrt(em_zt)

    Lscale_zm = clubb_api.zt2zm(
        nzm, nzt, ngrdcol, gr, Lscale,
    )
    Kh_zm = (c_K[:, None]
             * np.maximum(Lscale_zm, zero_threshold)
             * np.sqrt(np.maximum(em, em_min)))

    lhs_splat_wp2, lhs_splat_wp3 = clubb_api.wp23_term_splat_lhs(
        nzm, nzt, ngrdcol, gr, clubb_params[:, iC_wp2_splat - 1],
        brunt_vaisala_freq_sqd_mixed, Lscale_zm, rho_ds_zm,
    )

    #----------------------------------------------------------------
    # Set Surface variances
    #----------------------------------------------------------------
    # Surface variances should be set here, before the call to either
    # advance_xp2_xpyp or advance_wp2_wp3.
    # Surface effects should not be included with any case where the lowest
    # level is not the ground level.  Brian Griffin.  December 22, 2005.
    #
    # Diagnose surface variances based on surface fluxes.
    sfc_result = clubb_api.calc_sfc_varnce(
        nzm, nzt, ngrdcol, sclr_dim, sclr_idx,
        gr, dt, sfc_elevation,
        upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc,
        um, vm, Lscale_up, wpsclrp_sfc,
        lhs_splat_wp2, tau_zm,
        clubb_config_flags.l_vary_convect_depth,
        t0,
        clubb_params[:, iup2_sfc_coef - 1],
        clubb_params[:, ia_const - 1],
        wp2, up2, vp2,
        thlp2, rtp2, rtpthlp,
        sclrp2, sclrprtp, sclrpthlp,
        err_info,
    )
    (wp2, up2, vp2,
     thlp2, rtp2, rtpthlp,
     sclrp2, sclrprtp, sclrpthlp, err_info) = sfc_result
    err_code = err_info.err_code
    if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
        return

    ########################################################################
    ############### ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
    ########################################################################

    # Cx_fnc_Richardson is only used if one of these clubb_config_flags is true,
    # otherwise its value is irrelevant, set it to 0 to avoid NaN problems
    if clubb_config_flags.l_use_C7_Richardson or clubb_config_flags.l_use_C11_Richardson:
        Cx_fnc_Richardson = clubb_api.compute_cx_fnc_richardson(
            nzm, nzt, ngrdcol, gr,
            Lscale_zm, ddzt_umvm_sqd, rho_ds_zm,
            brunt_vaisala_freq_sqd,
            brunt_vaisala_freq_sqd_mixed,
            clubb_params,
            clubb_config_flags.l_use_shear_Richardson,
            clubb_config_flags.l_modify_limiters_for_cnvg_test,
        )
    else:
        Cx_fnc_Richardson = np.zeros((ngrdcol, nzm))

    # Loop over the 4 main advance subroutines -- advance_xm_wpxp,
    # advance_wp2_wp3, advance_xp2_xpyp, and advance_windm_edsclrm -- in the
    # order determined by order_xm_wpxp, order_wp2_wp3, order_xp2_xpyp, and
    # order_windm.
    for advance_iter in range(1, 5):

        if advance_iter == order_xm_wpxp:

            #----------------------------------------------------------------
            # Advance rtm/wprtp and thlm/wpthlp one time step
            #----------------------------------------------------------------
            #
            # Advance the prognostic equations for
            #   the scalar grid means (rtm, thlm, sclrm) and
            #   scalar turbulent fluxes (wprtp, wpthlp, and wpsclrp)
            #   by one time step.
            result = clubb_api.advance_xm_wpxp(
                nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt_advance,
                sigma_sqd_w, wm_zm, wm_zt, wp2,
                Lscale_zm, wp3,
                Kh_zt, Kh_zm, stability_correction,
                invrs_tau_C6_zm, tau_max_zm,
                wp2rtp, rtpthvp, rtm_forcing,
                wprtp_forcing, rtm_ref, wp2thlp,
                thlpthvp, thlm_forcing,
                wpthlp_forcing, thlm_ref, rho_ds_zm,
                rho_ds_zt, invrs_rho_ds_zm,
                invrs_rho_ds_zt, thv_ds_zm, rtp2,
                thlp2, w_1_zm, w_2_zm,
                varnce_w_1_zm, varnce_w_2_zm,
                mixt_frac_zm, l_implemented, em,
                wp2sclrp, sclrpthvp, sclrm_forcing,
                sclrp2, Cx_fnc_Richardson,
                pdf_implicit_coefs_terms,
                um_forcing, vm_forcing, ug, vg,
                wpthvp, fcor, fcor_y, um_ref, vm_ref,
                up2, vp2, uprcp, vprcp, rc_coef_zm,
                clubb_params, nu_vert_res_dep, ts_nudge,
                clubb_config_flags.iiPDF_type,
                clubb_config_flags.penta_solve_method,
                clubb_config_flags.tridiag_solve_method,
                clubb_config_flags.fill_holes_type,
                clubb_config_flags.l_predict_upwp_vpwp,
                clubb_config_flags.l_ho_nontrad_coriolis,
                clubb_config_flags.l_ho_trad_coriolis,
                clubb_config_flags.l_diffuse_rtm_and_thlm,
                clubb_config_flags.l_stability_correct_Kh_N2_zm,
                clubb_config_flags.l_godunov_upwind_wpxp_ta,
                clubb_config_flags.l_upwind_xm_ma,
                clubb_config_flags.l_uv_nudge,
                clubb_config_flags.l_tke_aniso,
                clubb_config_flags.l_diag_Lscale_from_tau,
                clubb_config_flags.l_use_C7_Richardson,
                clubb_config_flags.l_lmm_stepping,
                clubb_config_flags.l_enable_relaxed_clipping,
                clubb_config_flags.l_linearize_pbl_winds,
                clubb_config_flags.l_mono_flux_lim_thlm,
                clubb_config_flags.l_mono_flux_lim_rtm,
                clubb_config_flags.l_mono_flux_lim_um,
                clubb_config_flags.l_mono_flux_lim_vm,
                clubb_config_flags.l_mono_flux_lim_spikefix,
                wprtp_cl_num, wpthlp_cl_num,
                upwp_cl_num, vpwp_cl_num,
                rtm, wprtp, thlm,
                wpthlp, sclrm, wpsclrp, um, upwp,
                vm, vpwp, um_pert, vm_pert,
                upwp_pert, vpwp_pert,
                err_info,
            )
            (wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
             rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
             um_pert, vm_pert, upwp_pert, vpwp_pert, err_info) = result
            wprtp_cl_num = int(wprtp_cl_num)
            wpthlp_cl_num = int(wpthlp_cl_num)
            upwp_cl_num = int(upwp_cl_num)
            vpwp_cl_num = int(vpwp_cl_num)

            # Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
            # This code won't work unless rtm >= 0 !!!
            # We do not clip rcm_in_layer because rcm_in_layer only influences
            # radiation, and we do not want to bother recomputing it.  6 Aug 2009
            rcm = clubb_api.clip_rcm(
                nzt, ngrdcol, rtm,
                'rtm < rcm in advance_xm_wpxp',
                rcm,
            )
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

        elif advance_iter == order_xp2_xpyp:
            #----------------------------------------------------------------
            # Compute some of the variances and covariances.  These include the
            # variance of total water (rtp2), liquid water potential temperature
            # (thlp2), their covariance (rtpthlp), and the variance of horizontal
            # wind (up2 and vp2).  The variance of vertical velocity is computed
            # in a different section, which will come either earlier or later
            # depending on the chosen call order.
            #----------------------------------------------------------------
            #
            # We found that certain cases require a time tendency to run
            # at shorter timesteps so these are prognosed now.
            #
            # We found that if we call advance_xp2_xpyp first, we can use a longer timestep.
            #
            # Advance the prognostic equations
            #   for scalar variances and covariances,
            #   plus the horizontal wind variances by one time step, by one time step.
            result = clubb_api.advance_xp2_xpyp(
                nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx,
                invrs_tau_xp2_zm, invrs_tau_C4_zm,
                invrs_tau_C14_zm, wm_zm,
                rtm, wprtp, thlm, wpthlp, wpthvp, um, vm,
                wp2, wp3, upwp, vpwp,
                sigma_sqd_w, wprtp2, wpthlp2,
                wprtpthlp, Kh_zt, rtp2_forcing,
                thlp2_forcing, rtpthlp_forcing,
                rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
                thv_ds_zm, cloud_frac,
                pdf_implicit_coefs_terms,
                dt_advance, fcor_y,
                sclrm, wpsclrp,
                wpsclrp2, wpsclrprtp, wpsclrpthlp,
                lhs_splat_wp2,
                clubb_params, nu_vert_res_dep,
                clubb_config_flags.iiPDF_type,
                clubb_config_flags.tridiag_solve_method,
                clubb_config_flags.fill_holes_type,
                clubb_config_flags.l_ho_nontrad_coriolis,
                clubb_config_flags.l_min_xp2_from_corr_wx,
                clubb_config_flags.l_C2_cloud_frac,
                clubb_config_flags.l_upwind_xpyp_ta,
                clubb_config_flags.l_godunov_upwind_xpyp_ta,
                clubb_config_flags.l_lmm_stepping,
                l_implemented,
                rtp2, thlp2, rtpthlp, up2, vp2,
                sclrp2, sclrprtp, sclrpthlp, err_info,
            )
            rtp2, thlp2, rtpthlp, up2, vp2, sclrp2, sclrprtp, sclrpthlp, err_info = result

            #----------------------------------------------------------------
            # Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
            # after subroutine advance_xp2_xpyp updated xp2.
            #----------------------------------------------------------------
            (wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
             wprtp, wpthlp, upwp, vpwp,
             wpsclrp, upwp_pert, vpwp_pert) = clubb_api.clip_covars_denom(
                nzm, ngrdcol, sclr_dim,
                dt,
                rtp2, thlp2, up2, vp2, wp2,
                sclrp2, clubb_config_flags.l_tke_aniso,
                clubb_config_flags.l_linearize_pbl_winds,
                clubb_config_flags.l_predict_upwp_vpwp,
                wprtp_cl_num, wpthlp_cl_num,
                upwp_cl_num, vpwp_cl_num,
                wprtp, wpthlp, upwp, vpwp, wpsclrp,
                upwp_pert, vpwp_pert,
            )
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

        elif advance_iter == order_wp2_wp3:
            #----------------------------------------------------------------
            # Advance the 2nd- and 3rd-order moments
            #   of vertical velocity (wp2, wp3) by one timestep.
            #----------------------------------------------------------------
            result = clubb_api.advance_wp2_wp3(
                nzm, nzt, ngrdcol, gr, dt_advance,
                sfc_elevation, fcor_y, sigma_sqd_w, wm_zm,
                wm_zt,
                wpup2, wpvp2, wp2up2, wp2vp2, wp4,
                wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp,
                em, Kh_zm, Kh_zt, invrs_tau_C4_zm,
                invrs_tau_wp3_zt, invrs_tau_C1_zm,
                rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,
                invrs_rho_ds_zt, thv_ds_zm,
                thv_ds_zt, pdf_params.mixt_frac, Cx_fnc_Richardson,
                lhs_splat_wp2, lhs_splat_wp3,
                pdf_implicit_coefs_terms,
                wprtp, wpthlp, rtp2, thlp2,
                clubb_params, nu_vert_res_dep,
                clubb_config_flags.iiPDF_type,
                clubb_config_flags.penta_solve_method,
                clubb_config_flags.fill_holes_type,
                clubb_config_flags.l_min_wp2_from_corr_wx,
                clubb_config_flags.l_upwind_xm_ma,
                clubb_config_flags.l_tke_aniso,
                clubb_config_flags.l_standard_term_ta,
                clubb_config_flags.l_partial_upwind_wp3,
                clubb_config_flags.l_damp_wp2_using_em,
                clubb_config_flags.l_use_C11_Richardson,
                clubb_config_flags.l_damp_wp3_Skw_squared,
                clubb_config_flags.l_lmm_stepping,
                clubb_config_flags.l_use_tke_in_wp3_pr_turb_term,
                clubb_config_flags.l_use_tke_in_wp2_wp3_K_dfsn,
                clubb_config_flags.l_use_wp3_lim_with_smth_Heaviside,
                clubb_config_flags.l_wp2_fill_holes_tke,
                clubb_config_flags.l_ho_nontrad_coriolis,
                l_implemented,
                up2, vp2, wp2, wp3, err_info,
            )
            up2, vp2, wp2, wp3, err_info = result

            #----------------------------------------------------------------
            # Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
            # after subroutine advance_wp2_wp3 updated wp2.
            #----------------------------------------------------------------
            (wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
             wprtp, wpthlp, upwp, vpwp,
             wpsclrp, upwp_pert, vpwp_pert) = clubb_api.clip_covars_denom(
                nzm, ngrdcol, sclr_dim,
                dt,
                rtp2, thlp2, up2, vp2, wp2,
                sclrp2, clubb_config_flags.l_tke_aniso,
                clubb_config_flags.l_linearize_pbl_winds,
                clubb_config_flags.l_predict_upwp_vpwp,
                wprtp_cl_num, wpthlp_cl_num,
                upwp_cl_num, vpwp_cl_num,
                wprtp, wpthlp, upwp, vpwp, wpsclrp,
                upwp_pert, vpwp_pert,
            )
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

        elif advance_iter == order_windm:
            #----------------------------------------------------------------
            # Advance the horizontal mean winds (um, vm),
            #   the mean of the eddy-diffusivity scalars (i.e. edsclrm),
            #   and their fluxes (upwp, vpwp, wpedsclrp) by one time step.
            #----------------------------------------------------------------
            result = clubb_api.advance_windm_edsclrm(
                nzm, nzt, ngrdcol, edsclr_dim, gr, dt,
                wm_zt, Kh_zm, clubb_params,
                ug, vg, um_ref, vm_ref,
                wp2, up2, vp2, um_forcing, vm_forcing,
                edsclrm_forcing, p_in_Pa,
                rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,
                fcor, l_implemented,
                nu_vert_res_dep, ts_nudge,
                clubb_config_flags.tridiag_solve_method,
                clubb_config_flags.l_predict_upwp_vpwp,
                clubb_config_flags.l_upwind_xm_ma,
                clubb_config_flags.l_uv_nudge,
                clubb_config_flags.l_tke_aniso,
                clubb_config_flags.l_lmm_stepping,
                clubb_config_flags.l_linearize_pbl_winds,
                clubb_config_flags.l_do_expldiff_rtm_thlm,
                clubb_config_flags.fill_holes_type,
                upwp_cl_num, vpwp_cl_num,
                um, vm, thlm, rtm, edsclrm,
                upwp, vpwp, wpedsclrp,
                um_pert, vm_pert, upwp_pert,
                vpwp_pert, err_info,
            )
            (upwp_cl_num, vpwp_cl_num,
             um, vm, thlm, rtm, edsclrm, upwp, vpwp, wpedsclrp,
             um_pert, vm_pert, upwp_pert, vpwp_pert, err_info) = result
            upwp_cl_num = int(upwp_cl_num)
            vpwp_cl_num = int(vpwp_cl_num)
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

    #----------------------------------------------------------------
    # Advance or otherwise calculate <thl'^3>, <rt'^3>, and
    # <sclr'^3>.
    #----------------------------------------------------------------
    if l_advance_xp3 and clubb_config_flags.iiPDF_type != iiPDF_ADG1:
        # Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
        # simplified form of the <x'^3> predictive equation.  The simplified
        # <x'^3> equation can either be advanced from its previous value or
        # calculated using a steady-state approximation.
        rtp3, thlp3, sclrp3, up3, vp3 = clubb_api.advance_xp3(
            nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,
            rtm, thlm, rtp2, thlp2, wprtp,
            wpthlp, wprtp2, wpthlp2, rho_ds_zm,
            invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,
            sclrm, sclrp2, wpsclrp, wpsclrp2,
            wp2, wp3, upwp, vpwp, up2, vp2,
            thvm, clubb_params,
            clubb_config_flags.l_lmm_stepping,
            rtp3, thlp3, sclrp3, up3, vp3,
        )
    else:
        rtp3, thlp3, up3, vp3, sclrp3 = clubb_api.diagnose_xp3(
            nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr,
            clubb_config_flags.iiPDF_type, clubb_params,
            wp2, wp3, thvm,
            wprtp, wpthlp, rtp2, thlp2, upwp, vpwp, up2, vp2,
            sigma_sqd_w, wpsclrp, sclrp2,
            rtp3, thlp3, up3, vp3,
            sclrp3,
        )

    if (clubb_config_flags.ipdf_call_placement == ipdf_post_advance_fields
            or clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields):

        # Sample stats in this call to subroutine pdf_closure_driver for
        # ipdf_post_advance_fields, but not for ipdf_pre_post_advance_fields
        # because stats were sampled during the first call to subroutine
        # pdf_closure_driver.
        if clubb_config_flags.ipdf_call_placement == ipdf_post_advance_fields:
            l_samp_stats = True
        elif clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields:
            l_samp_stats = False

        #########################################################################
        ########                     CALL CLUBB's PDF                     #######
        ########   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
        #########################################################################
        # Given CLUBB's prognosed moments, diagnose CLUBB's PDF parameters
        #   and quantities integrated over that PDF, including
        #   quantities related to clouds, buoyancy, and turbulent advection.
        sigma_sqd_w = clubb_api.compute_sigma_sqd_w(
            nzm, nzt, ngrdcol, gr,
            wp3, wp2, thlp2, rtp2,
            up2, vp2, wpthlp, wprtp, upwp, vpwp,
            clubb_params,
            clubb_config_flags.l_predict_upwp_vpwp,
        )

        pdf_result = clubb_api.pdf_closure_driver(
            gr, nzm, nzt, ngrdcol,
            dt, hydromet_dim, sclr_dim, sclr_tol,
            wprtp, thlm, wpthlp, rtp2, rtp3,
            thlp2, thlp3, rtpthlp, wp2,
            wp3, wm_zm, wm_zt,
            um, up2, upwp, up3,
            vm, vp2, vpwp, vp3,
            p_in_Pa, exner,
            thv_ds_zm, thv_ds_zt, rtm_ref,
            wphydrometp,
            wp2hmp, rtphmp_zt, thlphmp_zt,
            sclrm, wpsclrp, sclrp2,
            sclrprtp, sclrpthlp, sclrp3,
            p_sfc, l_samp_stats,
            mixt_frac_max_mag, ts_nudge,
            rtm_min, rtm_nudge_max_altitude,
            clubb_params,
            clubb_config_flags.iiPDF_type,
            clubb_config_flags.saturation_formula,
            clubb_config_flags.l_rtm_nudge,
            clubb_config_flags.l_trapezoidal_rule_zt,
            clubb_config_flags.l_trapezoidal_rule_zm,
            clubb_config_flags.l_call_pdf_closure_twice,
            clubb_config_flags.l_use_cloud_cover,
            clubb_config_flags.l_rcm_supersat_adj,
            l_mix_rat_hm,
            rtm, sigma_sqd_w,
            pdf_implicit_coefs_terms,
            pdf_params, pdf_params_zm, err_info,
        )
        (rtm, pdf_implicit_coefs_terms, pdf_params, pdf_params_zm, err_info,
         rcm, cloud_frac, ice_supersat_frac, wprcp_out, sigma_sqd_w, wpthvp, wp2thvp,
         wp2up, rtpthvp, thlpthvp, rc_coef, rcm_in_layer, cloud_cover,
         rcp2_zt, thlprcp, rc_coef_zm, sclrpthvp, wpup2, wpvp2, wp2up2,
         wp2vp2, wp4, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, wp2rcp,
         rtprcp, rcp2, uprcp, vprcp, w_up_in_cloud, w_down_in_cloud,
         cloudy_updraft_frac, cloudy_downdraft_frac, skw_velocity,
         cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, thlm_zm, rcm_zm,
         rcm_supersat_adj, wp2sclrp, wpsclrp2, sclrprcp,
         wpsclrprtp, wpsclrpthlp) = pdf_result

        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return


    ########################################################################
    #############            ACCUMULATE STATISTICS            #############
    ########################################################################
    if l_sample:
        clubb_api.stats_accumulate(
            nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, gr, dt,
            l_implemented, clubb_config_flags.l_host_applies_sfc_fluxes,
            clubb_config_flags.l_stability_correct_tau_zm,
            clubb_params,
            um, vm, upwp, vpwp, up2, vp2,
            thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing,
            wpthlp_sfc, wprtp_sfc, wprtp, wpthlp,
            wp2, wp3, rtp2, rtp3, thlp2, thlp3,
            rtpthlp,
            p_in_Pa, exner, rho, rho_zm,
            rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt,
            wm_zt, wm_zm, rcm,
            cloud_frac,
            thvm, ug, vg,
            ddzt_umvm_sqd, stability_correction,
            Kh_zt,
            rsat,
            Kh_zm,
            em,
            sclrm, sclrp2,
            sclrprtp, sclrpthlp, sclrm_forcing,
            wpsclrp, wpedsclrp, edsclrm,
            edsclrm_forcing,
            clubb_config_flags.saturation_formula,
        )

        clubb_api.stats_finalize_budget("wp2_bt", wp2 / dt)
        clubb_api.stats_finalize_budget("vp2_bt", vp2 / dt)
        clubb_api.stats_finalize_budget("up2_bt", up2 / dt)
        clubb_api.stats_finalize_budget("wprtp_bt", wprtp / dt)
        clubb_api.stats_finalize_budget("wpthlp_bt", wpthlp / dt)
        if clubb_config_flags.l_predict_upwp_vpwp:
            clubb_api.stats_finalize_budget("upwp_bt", upwp / dt)
            clubb_api.stats_finalize_budget("vpwp_bt", vpwp / dt)
        clubb_api.stats_finalize_budget("rtp2_bt", rtp2 / dt)
        clubb_api.stats_finalize_budget("thlp2_bt", thlp2 / dt)
        clubb_api.stats_finalize_budget("rtpthlp_bt", rtpthlp / dt)
        clubb_api.stats_finalize_budget("rtm_bt", rtm / dt)
        clubb_api.stats_finalize_budget("thlm_bt", thlm / dt)
        clubb_api.stats_finalize_budget("um_bt", um / dt)
        clubb_api.stats_finalize_budget("vm_bt", vm / dt)
        clubb_api.stats_finalize_budget("wp3_bt", wp3 / dt)

    if clubb_api.clubb_at_least_debug_level(2):

        err_info = clubb_api.parameterization_check(
            nzm, nzt, ngrdcol, sclr_dim, edsclr_dim,
            thlm_forcing, rtm_forcing, um_forcing,
            vm_forcing, wm_zm, wm_zt, p_in_Pa,
            rho_zm, rho, exner, rho_ds_zm,
            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
            thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,
            vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,
            rtm, wprtp, thlm, wpthlp, wp2, wp3,
            rtp2, thlp2, rtpthlp,
            "end of ",
            wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,
            sclrp2,
            sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,
            edsclrm_forcing,
            err_info,
        )

        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

    # ================================================================== #
    return (
        um, vm, up3, vp3, thlm, rtm, rtp3, thlp3, wp3,
        upwp, vpwp, up2, vp2, wprtp, wpthlp, rtp2, thlp2, rtpthlp, wp2,
        sclrm, sclrp3, wpsclrp, sclrp2, sclrprtp, sclrpthlp,
        p_in_Pa, exner, rcm, cloud_frac, wp2thvp, wp2up,
        wpthvp, rtpthvp, thlpthvp, sclrpthvp,
        wp2rtp, wp2thlp, wpup2, wpvp2, ice_supersat_frac,
        uprcp, vprcp, rc_coef_zm, wp4, wp2up2, wp2vp2,
        um_pert, vm_pert, upwp_pert, vpwp_pert,
        edsclrm,
        pdf_params, pdf_params_zm,
        pdf_implicit_coefs_terms, err_info,
        rcm_in_layer, cloud_cover, w_up_in_cloud, w_down_in_cloud,
        cloudy_updraft_frac, cloudy_downdraft_frac, wprcp_out, invrs_tau_zm,
        Kh_zt, Kh_zm, thlprcp, Lscale,
        sigma_sqd_w,
        rc_coef, rcp2_zt, wprtp2,
        wpthlp2, wprtpthlp, wp2rcp,
        rtprcp, rcp2, skw_velocity,
        cloud_frac_zm, ice_supersat_frac_zm,
        rtm_zm, thlm_zm, rcm_zm,
        rcm_supersat_adj, wp2sclrp,
        wpsclrp2, sclrprcp,
        wpsclrprtp, wpsclrpthlp,
    )


def set_sfc_value_of_flux_profiles(
    nzm,
    ngrdcol,
    sclr_dim,
    edsclr_dim,
    gr,
    l_host_applies_sfc_fluxes,
    l_linearize_pbl_winds,
    wpthlp_sfc,
    wprtp_sfc,
    upwp_sfc,
    vpwp_sfc,
    upwp_sfc_pert,
    vpwp_sfc_pert,
    wpsclrp_sfc,
    wpedsclrp_sfc,
    wpthlp,
    wprtp,
    upwp,
    vpwp,
    upwp_pert,
    vpwp_pert,
    wpsclrp,
    wpedsclrp,
):
    k_lb = gr.k_lb_zm

    # SET SURFACE VALUES OF FLUXES (BROUGHT IN)
    # We only do this for host models that do not apply the flux
    # elsewhere in the code (e.g. WRF).  In other cases the _sfc variables will
    # only be used to compute the variance at the surface. -dschanen 8 Sept 2009
    if not l_host_applies_sfc_fluxes:
        wpthlp[:, k_lb] = wpthlp_sfc
        wprtp[:, k_lb] = wprtp_sfc
        upwp[:, k_lb] = upwp_sfc
        vpwp[:, k_lb] = vpwp_sfc

        if l_linearize_pbl_winds:
            upwp_pert[:, k_lb] = upwp_sfc_pert
            vpwp_pert[:, k_lb] = vpwp_sfc_pert

        # Set fluxes for passive scalars (if enabled)
        if sclr_dim > 0:
            for sclr in range(sclr_dim):
                wpsclrp[:, k_lb, sclr] = wpsclrp_sfc[:, sclr]

        if edsclr_dim > 0:
            # wpedsclrp is a local variable that is not saved from
            # timestep to timestep. Set it to 0 and overwrite the surface.
            wpedsclrp[:, :, :] = zero
            for edsclr in range(edsclr_dim):
                wpedsclrp[:, k_lb, edsclr] = wpedsclrp_sfc[:, edsclr]

    else:
        wpthlp[:, k_lb] = zero
        wprtp[:, k_lb] = zero
        upwp[:, k_lb] = zero
        vpwp[:, k_lb] = zero

        # Set fluxes for passive scalars (if enabled)
        if sclr_dim > 0:
            for sclr in range(sclr_dim):
                wpsclrp[:, k_lb, sclr] = zero

        if edsclr_dim > 0:
            # wpedsclrp is a local variable that is not saved from
            # timestep to timestep. Set it to 0 and overwrite the surface.
            wpedsclrp[:, :, :] = zero
            for edsclr in range(edsclr_dim):
                wpedsclrp[:, k_lb, edsclr] = zero

    return wpthlp, wprtp, upwp, vpwp, upwp_pert, vpwp_pert, wpsclrp, wpedsclrp

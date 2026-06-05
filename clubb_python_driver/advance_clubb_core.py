"""Python port of src/CLUBB_core/advance_clubb_core_module.F90.

Translates the Fortran advance_clubb_core subroutine into Python/NumPy,
calling individual Fortran subroutines via the F2PY API for the complex
advance/closure routines, and using NumPy for the inline array math.
"""

import numpy as np

from clubb_python import clubb_api
from clubb_python_driver.clubb_constants import (
    em_min,
    eps,
    grav,
    l_smooth_min_max,
    min_max_smth_mag,
    one,
    rt_tol,
    thl_tol,
    three_halves,
    two,
    unused_var,
    w_tol,
    w_tol_sqd,
    zero,
    zero_threshold,
    # Parameter indices
    ia3_coef_min,
    ia_const,
    ibv_efold,
    ic_K,
    ic_K10,
    ic_K10h,
    iC_wp2_splat,
    igamma_coef,
    igamma_coefb,
    igamma_coefc,
    ilambda0_stability_coef,
    imu,
    itaumax,
    iup2_sfc_coef,
    ixp3_coef_base,
    ixp3_coef_slope,
    # Model flag constants
    iiPDF_ADG1,
    ipdf_post_advance_fields,
    ipdf_pre_advance_fields,
    ipdf_pre_post_advance_fields,
    smth_type,
    tau_const,
    ufmin,
    below_grnd_val,
)

CLUBB_FATAL_ERROR = 99


def advance_clubb_core(
    *,
    gr,
    nzm,
    nzt,
    ngrdcol,
    dt_main,
    flags,
    sclr_dim,
    edsclr_dim,
    hydromet_dim,
    clubb_params,
    fcor,
    fcor_y,
    host_dx,
    host_dy,
    wm_zm,
    wm_zt,
    rho_ds_zt,
    rtm,
    thlm,
    rho,
    rfrzm,
    sfc_elevation,
    upwp_sfc,
    vpwp_sfc,
    wpthlp,
    wprtp_sfc,
    upwp,
    vpwp,
    upwp_sfc_pert,
    vpwp_sfc_pert,
    wpsclrp,
    wpedsclrp_sfc,
    p_sfc,
    thv_ds_zm,
    thv_ds_zt,
    wp2,
    wp3,
    thlp2,
    rtp2,
    rtpthlp,
    um,
    vm,
    p_in_Pa,
    exner,
    rcm,
    ice_supersat_frac,
    up2,
    vp2,
    wprtp,
    wpthlp_sfc,
    wp2thvp,
    wp2up,
    rtpthvp,
    thlpthvp,
    wpthvp,
    wphydrometp,
    wp2hmp,
    rtphmp_zt,
    thlphmp_zt,
    lmin,
    mixt_frac_max_mag,
    T0,
    ts_nudge,
    rtm_min,
    rtm_nudge_max_altitude,
    um_forcing,
    vm_forcing,
    thlm_forcing,
    rtm_forcing,
    wprtp_forcing,
    wpthlp_forcing,
    rtp2_forcing,
    thlp2_forcing,
    rtpthlp_forcing,
    err_info,
    sclr_tol,
    thlm_ref,
    rtm_ref,
    um_ref,
    vm_ref,
    ug,
    vg,
    sclrm_forcing,
    edsclrm_forcing,
    sclrp2,
    sclrprtp,
    sclrpthlp,
    sclr_idx,
    pdf_params,
    pdf_params_zm,
    pdf_implicit_coefs_terms,
    nu_vert_res_dep,
    sclrm,
    sclrpthvp,
    up3,
    vp3,
    um_pert,
    vm_pert,
    uprcp,
    vprcp,
    rc_coef_zm,
    wp4,
    wpup2,
    wpvp2,
    wp2up2,
    wp2vp2,
    wp2rtp,
    wp2thlp,
    upwp_pert,
    vpwp_pert,
    sclrp3,
    cloud_frac,
    thlp3,
    rtp3,
    edsclrm,
    wpsclrp_sfc,
    l_mix_rat_hm,
    rho_ds_zm,
    invrs_rho_ds_zm,
    invrs_rho_ds_zt,
    rho_zm,
    l_sample=False,
    l_gamma_Skw=True,
    l_advance_xp3=False,
    l_use_invrs_tau_N2_iso=False,
    order_xm_wpxp=1,
    order_xp2_xpyp=2,
    order_wp2_wp3=3,
    order_windm=4,
    debug_level=0,
):
    """Advance CLUBB one timestep with an explicit argument surface."""
    shzt = (ngrdcol, nzt)
    shzm = (ngrdcol, nzm)
    shzts = (ngrdcol, nzt, max(sclr_dim, 1))
    shzms = (ngrdcol, nzm, max(sclr_dim, 1))
    Kh_zm = np.zeros(shzm)
    Kh_zt = np.zeros(shzt)
    Lscale = np.zeros(shzt)
    invrs_tau_zm = np.zeros(shzm)
    _cloud_frac_zm = np.zeros(shzm)
    _ice_supersat_frac_zm = np.zeros(shzm)
    _rc_coef = np.zeros(shzt)
    _rcm_supersat_adj = np.zeros(shzt)
    _rcm_zm = np.zeros(shzm)
    _rcp2 = np.zeros(shzm)
    _rcp2_zt = np.zeros(shzt)
    _rtm_zm = np.zeros(shzm)
    _rtprcp = np.zeros(shzm)
    _sclrprcp = np.zeros(shzms)
    _sigma_sqd_w = np.zeros(shzm)
    _skw_velocity = np.zeros(shzm)
    _thlm_zm = np.zeros(shzm)
    _wp2rcp = np.zeros(shzt)
    _wp2sclrp = np.zeros(shzts)
    _wprtp2 = np.zeros(shzt)
    _wprtpthlp = np.zeros(shzt)
    _wpsclrp2 = np.zeros(shzts)
    _wpsclrprtp = np.zeros(shzts)
    _wpsclrpthlp = np.zeros(shzts)
    _wpthlp2 = np.zeros(shzt)
    cloud_cover = np.zeros(shzt)
    cloudy_downdraft_frac = np.zeros(shzt)
    cloudy_updraft_frac = np.zeros(shzt)
    rcm_in_layer = np.zeros(shzt)
    thlprcp = np.zeros(shzm)
    w_down_in_cloud = np.zeros(shzt)
    w_up_in_cloud = np.zeros(shzt)
    wprcp_out = np.zeros(shzm)
    rsat = np.zeros(shzt)
    dt = dt_main
    l_gamma_skw = l_gamma_Skw
    l_advance_xp3_flag = l_advance_xp3
    l_use_invrs_tau_n2_iso = l_use_invrs_tau_N2_iso
    order_xm_wpxp_val = order_xm_wpxp
    order_xp2_xpyp_val = order_xp2_xpyp
    order_wp2_wp3_val = order_wp2_wp3
    order_windm_val = order_windm

    l_implemented = False  # standalone mode

    # ================================================================== #
    # Block A: Setup
    # ================================================================== #
    dt_advance = two * dt if flags.l_lmm_stepping else dt

    Lscale_max = clubb_api.set_lscale_max(
        ngrdcol=ngrdcol,
        l_implemented=l_implemented,
        host_dx=host_dx,
        host_dy=host_dy,
    )

    # ================================================================== #
    # Block B: Stats — spurious source pre-integration
    # ================================================================== #
    thlm_before = thlm.copy()
    rtm_before = rtm.copy()

    if debug_level >= 2:
        err_info = clubb_api.parameterization_check(
            nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
            thlm_forcing=thlm_forcing, rtm_forcing=rtm_forcing,
            um_forcing=um_forcing, vm_forcing=vm_forcing,
            wm_zm=wm_zm, wm_zt=wm_zt, p_in_pa=p_in_Pa,
            rho_zm=rho_zm, rho=rho, exner=exner,
            rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
            invrs_rho_ds_zm=invrs_rho_ds_zm, invrs_rho_ds_zt=invrs_rho_ds_zt,
            thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt,
            wpthlp_sfc=wpthlp_sfc, wprtp_sfc=wprtp_sfc,
            upwp_sfc=upwp_sfc, vpwp_sfc=vpwp_sfc, p_sfc=p_sfc,
            um=um, upwp=upwp, vm=vm, vpwp=vpwp,
            up2=up2, vp2=vp2, rtm=rtm, wprtp=wprtp,
            thlm=thlm, wpthlp=wpthlp, wp2=wp2, wp3=wp3,
            rtp2=rtp2, thlp2=thlp2, rtpthlp=rtpthlp,
            prefix="beginning of ", wpsclrp_sfc=wpsclrp_sfc, wpedsclrp_sfc=wpedsclrp_sfc,
            sclrm=sclrm, wpsclrp=wpsclrp, sclrp2=sclrp2,
            sclrprtp=sclrprtp, sclrpthlp=sclrpthlp,
            sclrm_forcing=sclrm_forcing, edsclrm=edsclrm,
            edsclrm_forcing=edsclrm_forcing, err_info=err_info,
        )
        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

    # ================================================================== #
    # Block D: Stats — begin budget tracking
    # ================================================================== #
    if l_sample:
        clubb_api.stats_update("rfrzm", rfrzm)
        clubb_api.stats_begin_budget("wp2_bt", wp2 / dt)
        clubb_api.stats_begin_budget("vp2_bt", vp2 / dt)
        clubb_api.stats_begin_budget("up2_bt", up2 / dt)
        clubb_api.stats_begin_budget("wprtp_bt", wprtp / dt)
        clubb_api.stats_begin_budget("wpthlp_bt", wpthlp / dt)
        if flags.l_predict_upwp_vpwp:
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

    # ================================================================== #
    # Block E: Set surface boundary conditions
    # ================================================================== #
    k_lb = gr.k_lb_zm  # 0-based lower boundary for momentum levels

    if not flags.l_host_applies_sfc_fluxes:
        wpthlp[:, k_lb] = wpthlp_sfc
        wprtp[:, k_lb] = wprtp_sfc
        upwp[:, k_lb] = upwp_sfc
        vpwp[:, k_lb] = vpwp_sfc

        if flags.l_linearize_pbl_winds:
            upwp_pert[:, k_lb] = upwp_sfc_pert
            vpwp_pert[:, k_lb] = vpwp_sfc_pert

        if sclr_dim > 0:
            for sclr in range(sclr_dim):
                wpsclrp[:, k_lb, sclr] = wpsclrp_sfc[:, sclr]

        if edsclr_dim > 0:
            wpedsclrp = np.zeros((ngrdcol, nzm, edsclr_dim))
            for edsclr in range(edsclr_dim):
                wpedsclrp[:, k_lb, edsclr] = wpedsclrp_sfc[:, edsclr]
        else:
            wpedsclrp = np.zeros((ngrdcol, nzm, max(edsclr_dim, 1)))
    else:
        wpthlp[:, k_lb] = 0.0
        wprtp[:, k_lb] = 0.0
        upwp[:, k_lb] = 0.0
        vpwp[:, k_lb] = 0.0

        if sclr_dim > 0:
            for sclr in range(sclr_dim):
                wpsclrp[:, k_lb, sclr] = 0.0

        if edsclr_dim > 0:
            wpedsclrp = np.zeros((ngrdcol, nzm, edsclr_dim))
        else:
            wpedsclrp = np.zeros((ngrdcol, nzm, max(edsclr_dim, 1)))

    # ================================================================== #
    # Block F: Set mu
    # ================================================================== #
    # Standalone mode: mu from tunable parameters (no CLUBBND_CAM)
    # Parameter indices are 1-based into clubb_params(:, nparams)
    mu = clubb_params[:, imu - 1].copy()

    # ================================================================== #
    # Block G: Pre-advance PDF closure (conditional)
    # ================================================================== #
    if (flags.ipdf_call_placement == ipdf_pre_advance_fields
            or flags.ipdf_call_placement == ipdf_pre_post_advance_fields):

        l_samp_stats = True  # always sample for pre or pre_post

        pdf_result = clubb_api.pdf_closure_driver(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            dt=dt, hydromet_dim=hydromet_dim, sclr_dim=sclr_dim,
            sclr_tol=sclr_tol,
            wprtp=wprtp, thlm=thlm, wpthlp=wpthlp,
            rtp2=rtp2, rtp3=rtp3,
            thlp2=thlp2, thlp3=thlp3,
            rtpthlp=rtpthlp, wp2=wp2, wp3=wp3,
            wm_zm=wm_zm, wm_zt=wm_zt,
            um=um, up2=up2, upwp=upwp, up3=up3,
            vm=vm, vp2=vp2, vpwp=vpwp, vp3=vp3,
            p_in_pa=p_in_Pa, exner=exner,
            thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt,
            rtm_ref=rtm_ref,
            wphydrometp=wphydrometp,
            wp2hmp=wp2hmp,
            rtphmp_zt=rtphmp_zt,
            thlphmp_zt=thlphmp_zt,
            sclrm=sclrm, wpsclrp=wpsclrp,
            sclrp2=sclrp2, sclrprtp=sclrprtp,
            sclrpthlp=sclrpthlp, sclrp3=sclrp3,
            p_sfc=p_sfc,
            l_samp_stats_in_pdf_call=l_samp_stats,
            mixt_frac_max_mag=mixt_frac_max_mag,
            ts_nudge=ts_nudge,
            rtm_min=rtm_min,
            rtm_nudge_max_altitude=rtm_nudge_max_altitude,
            clubb_params=clubb_params,
            iiPDF_type=flags.iiPDF_type,
            saturation_formula=flags.saturation_formula,
            l_rtm_nudge=flags.l_rtm_nudge,
            l_trapezoidal_rule_zt=flags.l_trapezoidal_rule_zt,
            l_trapezoidal_rule_zm=flags.l_trapezoidal_rule_zm,
            l_call_pdf_closure_twice=flags.l_call_pdf_closure_twice,
            l_use_cloud_cover=flags.l_use_cloud_cover,
            l_rcm_supersat_adj=flags.l_rcm_supersat_adj,
            l_mix_rat_hm=l_mix_rat_hm,
            pdf_params=pdf_params,
            pdf_params_zm=pdf_params_zm,
            pdf_implicit_coefs_terms=pdf_implicit_coefs_terms,
            err_info=err_info,
            rtm=rtm,
        )
        (rtm, pdf_implicit_coefs_terms, pdf_params, pdf_params_zm, err_info,
         rcm, cloud_frac, ice_supersat_frac, wprcp_out, _sigma_sqd_w, wpthvp, wp2thvp,
         wp2up, rtpthvp, thlpthvp, _rc_coef, rcm_in_layer, cloud_cover,
         _rcp2_zt, thlprcp, rc_coef_zm, sclrpthvp, wpup2, wpvp2, wp2up2,
         wp2vp2, wp4, wp2rtp, _wprtp2, wp2thlp, _wpthlp2, _wprtpthlp, _wp2rcp,
         _rtprcp, _rcp2, uprcp, vprcp, w_up_in_cloud, w_down_in_cloud,
         cloudy_updraft_frac, cloudy_downdraft_frac, _skw_velocity,
         _cloud_frac_zm, _ice_supersat_frac_zm, _rtm_zm, _thlm_zm, _rcm_zm,
         _rcm_supersat_adj, _wp2sclrp, _wpsclrp2, _sclrprcp,
         _wpsclrprtp, _wpsclrpthlp) = pdf_result

        # Check for fatal error
        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

    # ================================================================== #
    # Block H: Interpolations — wp2_zt, sigma_sqd_w, a3_coef
    # ================================================================== #
    wp2 = wp2
    wp3 = wp3

    wp2_zt = np.maximum(
        clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=wp2),
        w_tol_sqd,
    )

    sigma_sqd_w = _sigma_sqd_w  # may be set by PDF closure above

    if flags.ipdf_call_placement == ipdf_post_advance_fields:
        # Calculate sigma_sqd_w here
        sigma_sqd_w = clubb_api.compute_sigma_sqd_w(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            wp3=wp3,
            wp2=wp2, thlp2=thlp2, rtp2=rtp2,
            up2=up2, vp2=vp2,
            wpthlp=wpthlp, wprtp=wprtp,
            upwp=upwp, vpwp=vpwp,
            clubb_params=clubb_params,
            l_predict_upwp_vpwp=flags.l_predict_upwp_vpwp,
        )

    if sigma_sqd_w is None:
        # If PDF closure was called pre-advance, sigma_sqd_w was set there
        # It should already be available from the pdf_closure unpack.
        sigma_sqd_w = np.zeros((ngrdcol, nzm))

    # a3 coefficient
    a3_coef = -two * (one - sigma_sqd_w) ** 2 + 3.0
    a3_min = clubb_params[:, ia3_coef_min - 1]
    for k in range(nzm):
        a3_coef[:, k] = np.maximum(a3_coef[:, k], a3_min)

    # Interpolate variances/covariances to thermodynamic levels
    thlp2_zt = np.maximum(
        clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=thlp2),
        thl_tol ** 2,
    )
    rtp2_zt = np.maximum(
        clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=rtp2),
        rt_tol ** 2,
    )
    rtpthlp_zt = clubb_api.zm2zt(
        gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=rtpthlp
    )
    if l_sample:
        clubb_api.stats_update("rtpthlp_zt", rtpthlp_zt)

    # wp3_on_wp2 on zt levels
    wp3_on_wp2_zt = wp3 / np.maximum(wp2_zt, w_tol_sqd)
    wp3_on_wp2_zt = np.clip(wp3_on_wp2_zt, -1000.0, 1000.0)

    wp3_on_wp2 = clubb_api.zt2zm(
        gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=wp3_on_wp2_zt
    )
    wp3_on_wp2_zt = clubb_api.zm2zt(
        gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=wp3_on_wp2
    )

    # ================================================================== #
    # Block I: Compute thvm
    # ================================================================== #
    thvm = clubb_api.calculate_thvm(
        nzt=nzt, ngrdcol=ngrdcol,
        thlm=thlm, rtm=rtm, rcm=rcm,
        exner=exner, thv_ds_zt=thv_ds_zt,
    )

    # ================================================================== #
    # Block J: TKE computation
    # ================================================================== #
    if not flags.l_tke_aniso:
        em = three_halves * wp2
    else:
        em = 0.5 * (wp2 + vp2 + up2)

    sqrt_em_zt = np.maximum(
        clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=em),
        em_min,
    )
    sqrt_em_zt = np.sqrt(sqrt_em_zt)

    # ================================================================== #
    # Block K: Brunt-Vaisala, wind shear, Richardson number
    # ================================================================== #
    bv_result = clubb_api.calc_brunt_vaisala_freq_sqd(
        gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
        thlm=thlm, exner=exner, rtm=rtm,
        rcm=rcm, p_in_Pa=p_in_Pa, thvm=thvm,
        ice_supersat_frac=ice_supersat_frac,
        saturation_formula=flags.saturation_formula,
        l_brunt_vaisala_freq_moist=flags.l_brunt_vaisala_freq_moist,
        l_use_thvm_in_bv_freq=flags.l_use_thvm_in_bv_freq,
        l_modify_limiters_for_cnvg_test=flags.l_modify_limiters_for_cnvg_test,
        bv_efold=clubb_params[:, ibv_efold - 1], T0=T0,
    )
    (brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed,
     brunt_vaisala_freq_sqd_smth) = bv_result

    ddzt_um = clubb_api.ddzt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=um)
    ddzt_vm = clubb_api.ddzt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=vm)
    ddzt_umvm_sqd = ddzt_um ** 2 + ddzt_vm ** 2

    if l_sample:
        clubb_api.stats_update("ddzt_umvm_sqd", ddzt_umvm_sqd)

    # Richardson number
    if flags.l_modify_limiters_for_cnvg_test:
        Ri_zm = clubb_api.calc_ri_zm(
            nzm=nzm, ngrdcol=ngrdcol,
            bv_freq_sqd=brunt_vaisala_freq_sqd_smth,
            shear=ddzt_umvm_sqd,
            lim_bv=0.0, lim_shear=1.0e-12,
        )
        Ri_zm = clubb_api.zm2zt2zm(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=Ri_zm, zm_min=0.0
        )
    else:
        if l_smooth_min_max:
            brunt_vaisala_freq_clipped = clubb_api.smooth_max(
                nz=nzm, ngrdcol=ngrdcol,
                input_var1=1.0e-7,
                input_var2=brunt_vaisala_freq_sqd_smth,
                smth_coef=1.0e-4 * min_max_smth_mag,
            )
            ddzt_umvm_sqd_clipped = clubb_api.smooth_max(
                nz=nzm, ngrdcol=ngrdcol,
                input_var1=ddzt_umvm_sqd,
                input_var2=1.0e-7,
                smth_coef=1.0e-6 * min_max_smth_mag,
            )
            Ri_zm = clubb_api.calc_ri_zm(
                nzm=nzm, ngrdcol=ngrdcol,
                bv_freq_sqd=brunt_vaisala_freq_clipped,
                shear=ddzt_umvm_sqd_clipped,
                lim_bv=0.0, lim_shear=0.0,
            )
        else:
            Ri_zm = clubb_api.calc_ri_zm(
                nzm=nzm, ngrdcol=ngrdcol,
                bv_freq_sqd=brunt_vaisala_freq_sqd_smth,
                shear=ddzt_umvm_sqd,
                lim_bv=1.0e-7, lim_shear=1.0e-7,
            )

    # ================================================================== #
    # Block L: Mixing length / dissipation time scale
    # ================================================================== #
    if not flags.l_diag_Lscale_from_tau:
        lscale_result = clubb_api.calc_lscale_directly(
            gr=gr, ngrdcol=ngrdcol, nzm=nzm, nzt=nzt,
            l_implemented=l_implemented, p_in_pa=p_in_Pa,
            exner=exner, rtm=rtm, thlm=thlm,
            thvm=thvm, newmu=mu,
            rtp2_zt=rtp2_zt, thlp2_zt=thlp2_zt, rtpthlp_zt=rtpthlp_zt,
            em=em, thv_ds_zt=thv_ds_zt,
            lscale_max=Lscale_max, lmin=lmin,
            clubb_params=clubb_params,
            saturation_formula=flags.saturation_formula,
            l_lscale_plume_centered=flags.l_Lscale_plume_centered,
            pdf_params=pdf_params,
            err_info=err_info,
        )
        err_info, Lscale, Lscale_up, Lscale_down = lscale_result
        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

        # tau from Lscale
        tau_zt = np.minimum(Lscale / sqrt_em_zt,
                            clubb_params[:, itaumax - 1:itaumax])

        Lscale_zm = np.maximum(
            clubb_api.zt2zm(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=Lscale),
            zero_threshold,
        )

        tau_zm = np.minimum(
            Lscale_zm / np.sqrt(np.maximum(em_min, em)),
            clubb_params[:, itaumax - 1:itaumax],
        )

        invrs_tau_zm = one / tau_zm
        invrs_tau_wp2_zm = invrs_tau_zm.copy()
        invrs_tau_xp2_zm = invrs_tau_zm.copy()
        invrs_tau_wpxp_zm = invrs_tau_zm.copy()
        invrs_tau_wp3_zm = invrs_tau_zm.copy()
        tau_max_zm = np.broadcast_to(
            clubb_params[:, itaumax - 1:itaumax], (ngrdcol, nzm)).copy()

        invrs_tau_zt = one / tau_zt
        invrs_tau_wp3_zt = invrs_tau_zt.copy()
        tau_max_zt = np.broadcast_to(
            clubb_params[:, itaumax - 1:itaumax], (ngrdcol, nzt)).copy()

        # Placeholder variables not computed in this branch
        invrs_tau_no_N2_zm = np.zeros((ngrdcol, nzm))
        invrs_tau_bkgnd = np.zeros((ngrdcol, nzm))
        invrs_tau_shear = np.zeros((ngrdcol, nzm))
        invrs_tau_sfc = np.zeros((ngrdcol, nzm))
        invrs_tau_N2_iso = np.zeros((ngrdcol, nzm))
    else:
        lscale_result = clubb_api.diagnose_lscale_from_tau(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            upwp_sfc=upwp_sfc, vpwp_sfc=vpwp_sfc,
            ddzt_umvm_sqd=ddzt_umvm_sqd,
            ice_supersat_frac=ice_supersat_frac,
            em=em, sqrt_em_zt=sqrt_em_zt,
            ufmin=ufmin, tau_const=tau_const,
            sfc_elevation=sfc_elevation,
            lscale_max=Lscale_max,
            clubb_params=clubb_params,
            l_e3sm_config=flags.l_e3sm_config,
            l_smooth_heaviside_tau_wpxp=flags.l_smooth_Heaviside_tau_wpxp,
            brunt_vaisala_freq_sqd_smth=brunt_vaisala_freq_sqd_smth,
            ri_zm=Ri_zm,
            err_info=err_info,
        )
        (err_info,
         invrs_tau_zt, invrs_tau_zm,
         invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd,
         invrs_tau_shear, invrs_tau_N2_iso,
         invrs_tau_wp2_zm, invrs_tau_xp2_zm,
         invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm,
         tau_max_zm, tau_max_zt, tau_zm, tau_zt,
         Lscale, Lscale_up, Lscale_down) = lscale_result
        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

    if l_sample:
        clubb_api.stats_update("tau_zt", tau_zt)

    # ================================================================== #
    # Block M: Eddy diffusivity
    # ================================================================== #
    # Kh_zt = c_K * Lscale * sqrt_em_zt
    c_K = clubb_params[:, ic_K - 1]
    Kh_zt = c_K[:, None] * Lscale * sqrt_em_zt

    Lscale_zm = clubb_api.zt2zm(
        gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=Lscale
    )
    Kh_zm = (c_K[:, None]
             * np.maximum(Lscale_zm, zero_threshold)
             * np.sqrt(np.maximum(em, em_min)))

    lhs_splat_wp2, lhs_splat_wp3 = clubb_api.wp23_term_splat_lhs(
        nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, gr=gr,
        c_wp2_splat=clubb_params[:, iC_wp2_splat - 1],
        brunt_vaisala_freq_sqd_mixed=brunt_vaisala_freq_sqd_mixed,
        lscale_zm=Lscale_zm,
        rho_ds_zm=rho_ds_zm,
    )

    # ================================================================== #
    # Block N: Surface variances
    # ================================================================== #
    sfc_result = clubb_api.calc_sfc_varnce(
        gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, sclr_dim=sclr_dim,
        dt=dt, sfc_elevation=sfc_elevation,
        upwp_sfc=upwp_sfc, vpwp_sfc=vpwp_sfc,
        wpthlp=wpthlp, wprtp_sfc=wprtp_sfc,
        um=um, vm=vm,
        lscale_up=Lscale_up, wpsclrp_sfc=wpsclrp_sfc,
        lhs_splat_wp2=lhs_splat_wp2, tau_zm=tau_zm,
        l_vary_convect_depth=flags.l_vary_convect_depth,
        t0=T0,
        up2_sfc_coef=clubb_params[:, iup2_sfc_coef - 1],
        a_const=clubb_params[:, ia_const - 1],
        wp2=wp2, up2=up2, vp2=vp2,
        thlp2=thlp2, rtp2=rtp2, rtpthlp=rtpthlp,
        sclrp2=sclrp2, sclrprtp=sclrprtp,
        sclrpthlp=sclrpthlp,
        sclr_idx=sclr_idx,
        err_info=err_info,
    )
    # F2PY returns: wp2, up2, vp2, thlp2, rtp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp
    (wp2, up2, vp2,
     thlp2, rtp2, rtpthlp,
     sclrp2, sclrprtp, sclrpthlp, err_info) = sfc_result
    err_code = err_info.err_code
    if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
        return

    # Update local aliases after sfc_varnce modified them
    wp2 = wp2

    # ================================================================== #
    # Block O: Stats — pre-advance outputs (rvm, rel_humidity)
    # ================================================================== #
    if l_sample:
        clubb_api.stats_update("rvm", rtm - rcm)
        if clubb_api.var_on_stats_list("rel_humidity") or clubb_api.var_on_stats_list("rsat"):
            T_in_K = np.empty((ngrdcol, nzt))
            for i in range(ngrdcol):
                T_in_K[i, :] = clubb_api.thlm2t_in_k(
                    nzt, thlm[i, :], exner[i, :], rcm[i, :]
                )
            rsat = clubb_api.sat_mixrat_liq(
                gr=gr,
                nz=nzt,
                ngrdcol=ngrdcol,
                p_in_Pa=p_in_Pa,
                T_in_K=T_in_K,
                saturation_formula=flags.saturation_formula,
            )
            rel_humidity = (rtm - rcm) / rsat
            clubb_api.stats_update("rel_humidity", rel_humidity)

    # ================================================================== #
    # Block P: Extract PDF params for zm grid
    # ================================================================== #
    pdf_params = pdf_params

    if flags.l_call_pdf_closure_twice:
        pdf_params_zm = pdf_params_zm
        w_1_zm = pdf_params_zm.w_1.copy()
        w_2_zm = pdf_params_zm.w_2.copy()
        varnce_w_1_zm = pdf_params_zm.varnce_w_1.copy()
        varnce_w_2_zm = pdf_params_zm.varnce_w_2.copy()
        mixt_frac_zm = pdf_params_zm.mixt_frac.copy()
    else:
        w_1_zm = clubb_api.zt2zm(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=pdf_params.w_1)
        w_2_zm = clubb_api.zt2zm(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=pdf_params.w_2)
        varnce_w_1_zm = clubb_api.zt2zm(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=pdf_params.varnce_w_1
        )
        varnce_w_2_zm = clubb_api.zt2zm(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=pdf_params.varnce_w_2
        )
        mixt_frac_zm = clubb_api.zt2zm(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azt=pdf_params.mixt_frac
        )

    # ================================================================== #
    # Block Q: Stability correction, invrs_tau
    # ================================================================== #
    if flags.l_stability_correct_tau_zm:
        stability_correction = clubb_api.calc_stability_correction(
            nzm=nzm, ngrdcol=ngrdcol,
            brunt_vaisala_freq_sqd=brunt_vaisala_freq_sqd,
            lscale_zm=Lscale_zm, em=em,
            lambda0_stability_coef=clubb_params[:, ilambda0_stability_coef - 1],
        )
        if l_sample:
            clubb_api.stats_update("stability_correction", stability_correction)

        invrs_tau_N2_zm = invrs_tau_zm * stability_correction
        invrs_tau_C6_zm = invrs_tau_N2_zm.copy()
        invrs_tau_C1_zm = invrs_tau_N2_zm.copy()
    else:
        stability_correction = np.zeros((ngrdcol, nzm))
        invrs_tau_N2_zm = np.full((ngrdcol, nzm), unused_var)
        invrs_tau_C6_zm = invrs_tau_wpxp_zm.copy()
        invrs_tau_C1_zm = invrs_tau_wp2_zm.copy()

    # C14 always uses wp2 tau
    invrs_tau_C14_zm = invrs_tau_wp2_zm.copy()

    # C4 tau
    if (not flags.l_diag_Lscale_from_tau) and l_use_invrs_tau_n2_iso:
        raise RuntimeError(
            "Error! l_use_invrs_tau_N2_iso is not used when "
            "l_diag_Lscale_from_tau=false."
            "If you want to use Lscale code, go to file "
            "src/CLUBB_core/advance_clubb_core_module.F90 and "
            "change l_use_invrs_tau_N2_iso to false"
        )
    if not l_use_invrs_tau_n2_iso:
        invrs_tau_C4_zm = invrs_tau_wp2_zm.copy()
    else:
        invrs_tau_C4_zm = invrs_tau_N2_iso.copy()

    if l_sample:
        clubb_api.stats_update("invrs_tau_zm", invrs_tau_zm)
        clubb_api.stats_update("invrs_tau_xp2_zm", invrs_tau_xp2_zm)
        clubb_api.stats_update("invrs_tau_wp2_zm", invrs_tau_wp2_zm)
        clubb_api.stats_update("invrs_tau_wpxp_zm", invrs_tau_wpxp_zm)
        clubb_api.stats_update("Ri_zm", Ri_zm)
        clubb_api.stats_update("invrs_tau_wp3_zm", invrs_tau_wp3_zm)
        if flags.l_diag_Lscale_from_tau:
            clubb_api.stats_update("invrs_tau_no_N2_zm", invrs_tau_no_N2_zm)
            clubb_api.stats_update("invrs_tau_bkgnd", invrs_tau_bkgnd)
            clubb_api.stats_update("invrs_tau_sfc", invrs_tau_sfc)
            clubb_api.stats_update("invrs_tau_shear", invrs_tau_shear)
    # ================================================================== #
    # Block R: Cx_fnc_Richardson
    # ================================================================== #
    if flags.l_use_C7_Richardson or flags.l_use_C11_Richardson:
        Cx_fnc_Richardson = clubb_api.compute_cx_fnc_richardson(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            Lscale_zm=Lscale_zm, ddzt_umvm_sqd=ddzt_umvm_sqd,
            rho_ds_zm=rho_ds_zm,
            brunt_vaisala_freq_sqd=brunt_vaisala_freq_sqd,
            brunt_vaisala_freq_sqd_mixed=brunt_vaisala_freq_sqd_mixed,
            clubb_params=clubb_params,
            l_use_shear_Richardson=flags.l_use_shear_Richardson,
            l_modify_limiters_for_cnvg_test=flags.l_modify_limiters_for_cnvg_test,
        )
    else:
        Cx_fnc_Richardson = np.zeros((ngrdcol, nzm))

    # ================================================================== #
    # Block S: Main advance loop (4 iterations)
    # ================================================================== #
    # Locally-needed PDF closure outputs for advance routines
    wprtp2 = _wprtp2
    wpthlp2 = _wpthlp2
    wprtpthlp = _wprtpthlp
    wp2sclrp = _wp2sclrp
    wpsclrp2 = _wpsclrp2
    wpsclrprtp = _wpsclrprtp
    wpsclrpthlp = _wpsclrpthlp

    wprtp_cl_num = 0
    wpthlp_cl_num = 0
    upwp_cl_num = 0
    vpwp_cl_num = 0
    wprtp_cl_max = 3
    wpthlp_cl_max = 3
    upwp_cl_max = 3
    vpwp_cl_max = 3

    for advance_iter in range(1, 5):

        if advance_iter == order_xm_wpxp_val:
            result = clubb_api.advance_xm_wpxp(
                gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
                sclr_dim=sclr_dim, sclr_tol=sclr_tol, dt=dt_advance,
                sigma_sqd_w=sigma_sqd_w, wm_zm=wm_zm, wm_zt=wm_zt, wp2=wp2,
                lscale_zm=Lscale_zm, wp3=wp3,
                kh_zt=Kh_zt, kh_zm=Kh_zm, stability_correction=stability_correction,
                invrs_tau_c6_zm=invrs_tau_C6_zm, tau_max_zm=tau_max_zm,
                wp2rtp=wp2rtp, rtpthvp=rtpthvp, rtm_forcing=rtm_forcing,
                wprtp_forcing=wprtp_forcing, rtm_ref=rtm_ref, wp2thlp=wp2thlp,
                thlpthvp=thlpthvp, thlm_forcing=thlm_forcing,
                wpthlp_forcing=wpthlp_forcing, thlm_ref=thlm_ref, rho_ds_zm=rho_ds_zm,
                rho_ds_zt=rho_ds_zt, invrs_rho_ds_zm=invrs_rho_ds_zm,
                invrs_rho_ds_zt=invrs_rho_ds_zt, thv_ds_zm=thv_ds_zm, rtp2=rtp2,
                thlp2=thlp2, w_1_zm=w_1_zm, w_2_zm=w_2_zm,
                varnce_w_1_zm=varnce_w_1_zm, varnce_w_2_zm=varnce_w_2_zm,
                mixt_frac_zm=mixt_frac_zm, l_implemented=l_implemented, em=em,
                wp2sclrp=wp2sclrp, sclrpthvp=sclrpthvp, sclrm_forcing=sclrm_forcing,
                sclrp2=sclrp2, cx_fnc_richardson=Cx_fnc_Richardson,
                um_forcing=um_forcing, vm_forcing=vm_forcing, ug=ug, vg=vg,
                wpthvp=wpthvp, fcor=fcor, fcor_y=fcor_y, um_ref=um_ref, vm_ref=vm_ref,
                up2=up2, vp2=vp2, uprcp=uprcp, vprcp=vprcp, rc_coef_zm=rc_coef_zm,
                clubb_params=clubb_params, ts_nudge=ts_nudge,
                iipdf_type=flags.iiPDF_type, penta_solve_method=flags.penta_solve_method,
                tridiag_solve_method=flags.tridiag_solve_method,
                fill_holes_type=flags.fill_holes_type,
                l_predict_upwp_vpwp=flags.l_predict_upwp_vpwp,
                l_ho_nontrad_coriolis=flags.l_ho_nontrad_coriolis,
                l_ho_trad_coriolis=flags.l_ho_trad_coriolis,
                l_diffuse_rtm_and_thlm=flags.l_diffuse_rtm_and_thlm,
                l_stability_correct_kh_n2_zm=flags.l_stability_correct_Kh_N2_zm,
                l_godunov_upwind_wpxp_ta=flags.l_godunov_upwind_wpxp_ta,
                l_upwind_xm_ma=flags.l_upwind_xm_ma, l_uv_nudge=flags.l_uv_nudge,
                l_tke_aniso=flags.l_tke_aniso,
                l_diag_lscale_from_tau=flags.l_diag_Lscale_from_tau,
                l_use_c7_richardson=flags.l_use_C7_Richardson,
                l_lmm_stepping=flags.l_lmm_stepping,
                l_enable_relaxed_clipping=flags.l_enable_relaxed_clipping,
                l_linearize_pbl_winds=flags.l_linearize_pbl_winds,
                l_mono_flux_lim_thlm=flags.l_mono_flux_lim_thlm,
                l_mono_flux_lim_rtm=flags.l_mono_flux_lim_rtm,
                l_mono_flux_lim_um=flags.l_mono_flux_lim_um,
                l_mono_flux_lim_vm=flags.l_mono_flux_lim_vm,
                l_mono_flux_lim_spikefix=flags.l_mono_flux_lim_spikefix,
                wprtp_cl_num=wprtp_cl_num, wpthlp_cl_num=wpthlp_cl_num,
                upwp_cl_num=upwp_cl_num, vpwp_cl_num=vpwp_cl_num,
                rtm=rtm, wprtp=wprtp, thlm=thlm,
                wpthlp=wpthlp, sclrm=sclrm, wpsclrp=wpsclrp, um=um, upwp=upwp,
                vm=vm, vpwp=vpwp, um_pert=um_pert, vm_pert=vm_pert,
                upwp_pert=upwp_pert, vpwp_pert=vpwp_pert,
                nu_vert_res_dep=nu_vert_res_dep,
                pdf_implicit_coefs_terms=pdf_implicit_coefs_terms, err_info=err_info,
            )
            (wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
             rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
             um_pert, vm_pert, upwp_pert, vpwp_pert, err_info) = result
            wprtp_cl_num = int(wprtp_cl_num)
            wpthlp_cl_num = int(wpthlp_cl_num)
            upwp_cl_num = int(upwp_cl_num)
            vpwp_cl_num = int(vpwp_cl_num)
            rcm = clubb_api.clip_rcm(
                nzt=nzt, ngrdcol=ngrdcol, rtm=rtm,
                message='rtm < rcm in advance_xm_wpxp', rcm=rcm,
            )
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

        elif advance_iter == order_xp2_xpyp_val:
            result = clubb_api.advance_xp2_xpyp(
                gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
                sclr_dim=sclr_dim, sclr_tol=sclr_tol,
                invrs_tau_xp2_zm=invrs_tau_xp2_zm, invrs_tau_c4_zm=invrs_tau_C4_zm,
                invrs_tau_c14_zm=invrs_tau_C14_zm, wm_zm=wm_zm, rtm=rtm, wprtp=wprtp,
                thlm=thlm, wpthlp=wpthlp, wpthvp=wpthvp, um=um, vm=vm,
                wp2=wp2, wp2_zt=wp2_zt, wp3=wp3, upwp=upwp, vpwp=vpwp,
                sigma_sqd_w=sigma_sqd_w, wprtp2=wprtp2, wpthlp2=wpthlp2,
                wprtpthlp=wprtpthlp, kh_zt=Kh_zt, rtp2_forcing=rtp2_forcing,
                thlp2_forcing=thlp2_forcing, rtpthlp_forcing=rtpthlp_forcing,
                rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
                invrs_rho_ds_zm=invrs_rho_ds_zm, thv_ds_zm=thv_ds_zm,
                cloud_frac=cloud_frac, dt=dt_advance, fcor_y=fcor_y,
                sclrm=sclrm, wpsclrp=wpsclrp, wpsclrp2=wpsclrp2,
                wpsclrprtp=wpsclrprtp, wpsclrpthlp=wpsclrpthlp,
                lhs_splat_wp2=lhs_splat_wp2, clubb_params=clubb_params,
                iipdf_type=flags.iiPDF_type,
                tridiag_solve_method=flags.tridiag_solve_method,
                fill_holes_type=flags.fill_holes_type,
                l_ho_nontrad_coriolis=flags.l_ho_nontrad_coriolis,
                l_min_xp2_from_corr_wx=flags.l_min_xp2_from_corr_wx,
                l_c2_cloud_frac=flags.l_C2_cloud_frac,
                l_upwind_xpyp_ta=flags.l_upwind_xpyp_ta,
                l_godunov_upwind_xpyp_ta=flags.l_godunov_upwind_xpyp_ta,
                l_lmm_stepping=flags.l_lmm_stepping, rtp2=rtp2, thlp2=thlp2,
                rtpthlp=rtpthlp, up2=up2, vp2=vp2, sclrp2=sclrp2,
                sclrprtp=sclrprtp, sclrpthlp=sclrpthlp, sclr_idx=sclr_idx,
                nu_vert_res_dep=nu_vert_res_dep,
                pdf_implicit_coefs_terms=pdf_implicit_coefs_terms, err_info=err_info,
            )
            rtp2, thlp2, rtpthlp, up2, vp2, sclrp2, sclrprtp, sclrpthlp, err_info = result
            (wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
             wprtp, wpthlp, upwp, vpwp,
             wpsclrp, upwp_pert, vpwp_pert) = clubb_api.clip_covars_denom(
                nzm=nzm, ngrdcol=ngrdcol, sclr_dim=sclr_dim, dt=dt,
                rtp2=rtp2, thlp2=thlp2, up2=up2, vp2=vp2, wp2=wp2, sclrp2=sclrp2,
                l_tke_aniso=flags.l_tke_aniso,
                l_linearize_pbl_winds=flags.l_linearize_pbl_winds,
                l_predict_upwp_vpwp=flags.l_predict_upwp_vpwp,
                wprtp_cl_num=wprtp_cl_num, wpthlp_cl_num=wpthlp_cl_num,
                upwp_cl_num=upwp_cl_num, vpwp_cl_num=vpwp_cl_num,
                wprtp=wprtp, wpthlp=wpthlp, upwp=upwp, vpwp=vpwp,
                wpsclrp=wpsclrp, upwp_pert=upwp_pert, vpwp_pert=vpwp_pert,
            )
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

        elif advance_iter == order_wp2_wp3_val:
            result = clubb_api.advance_wp2_wp3(
                gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, dt=dt_advance,
                sfc_elevation=sfc_elevation, fcor_y=fcor_y, sigma_sqd_w=sigma_sqd_w,
                wm_zm=wm_zm, wm_zt=wm_zt,
                wpup2=wpup2, wpvp2=wpvp2, wp2up2=wp2up2,
                wp2vp2=wp2vp2, wp4=wp4, wpthvp=wpthvp, wp2thvp=wp2thvp, wp2up=wp2up,
                um=um, vm=vm, upwp=upwp, vpwp=vpwp, em=em, kh_zm=Kh_zm, kh_zt=Kh_zt,
                invrs_tau_c4_zm=invrs_tau_C4_zm, invrs_tau_wp3_zt=invrs_tau_wp3_zt,
                invrs_tau_c1_zm=invrs_tau_C1_zm,
                rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
                invrs_rho_ds_zm=invrs_rho_ds_zm, invrs_rho_ds_zt=invrs_rho_ds_zt,
                thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt,
                mixt_frac=pdf_params.mixt_frac, cx_fnc_richardson=Cx_fnc_Richardson,
                lhs_splat_wp2=lhs_splat_wp2, lhs_splat_wp3=lhs_splat_wp3,
                wprtp=wprtp, wpthlp=wpthlp, rtp2=rtp2, thlp2=thlp2,
                clubb_params=clubb_params, iipdf_type=flags.iiPDF_type,
                penta_solve_method=flags.penta_solve_method,
                fill_holes_type=flags.fill_holes_type,
                l_min_wp2_from_corr_wx=flags.l_min_wp2_from_corr_wx,
                l_upwind_xm_ma=flags.l_upwind_xm_ma, l_tke_aniso=flags.l_tke_aniso,
                l_standard_term_ta=flags.l_standard_term_ta,
                l_partial_upwind_wp3=flags.l_partial_upwind_wp3,
                l_damp_wp2_using_em=flags.l_damp_wp2_using_em,
                l_use_c11_richardson=flags.l_use_C11_Richardson,
                l_damp_wp3_skw_squared=flags.l_damp_wp3_Skw_squared,
                l_lmm_stepping=flags.l_lmm_stepping,
                l_use_tke_in_wp3_pr_turb_term=flags.l_use_tke_in_wp3_pr_turb_term,
                l_use_tke_in_wp2_wp3_k_dfsn=flags.l_use_tke_in_wp2_wp3_K_dfsn,
                l_use_wp3_lim_with_smth_heaviside=flags.l_use_wp3_lim_with_smth_Heaviside,
                l_wp2_fill_holes_tke=flags.l_wp2_fill_holes_tke,
                l_ho_nontrad_coriolis=flags.l_ho_nontrad_coriolis,
                up2=up2, vp2=vp2, wp2=wp2, wp3=wp3, wp2_zt=wp2_zt,
                nu_vert_res_dep=nu_vert_res_dep,
                pdf_implicit_coefs_terms=pdf_implicit_coefs_terms, err_info=err_info,
            )
            up2, vp2, wp2, wp3, wp2_zt, err_info = result
            (wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,
             wprtp, wpthlp, upwp, vpwp,
             wpsclrp, upwp_pert, vpwp_pert) = clubb_api.clip_covars_denom(
                nzm=nzm, ngrdcol=ngrdcol, sclr_dim=sclr_dim, dt=dt,
                rtp2=rtp2, thlp2=thlp2, up2=up2, vp2=vp2, wp2=wp2, sclrp2=sclrp2,
                l_tke_aniso=flags.l_tke_aniso,
                l_linearize_pbl_winds=flags.l_linearize_pbl_winds,
                l_predict_upwp_vpwp=flags.l_predict_upwp_vpwp,
                wprtp_cl_num=wprtp_cl_num, wpthlp_cl_num=wpthlp_cl_num,
                upwp_cl_num=upwp_cl_num, vpwp_cl_num=vpwp_cl_num,
                wprtp=wprtp, wpthlp=wpthlp, upwp=upwp, vpwp=vpwp,
                wpsclrp=wpsclrp, upwp_pert=upwp_pert, vpwp_pert=vpwp_pert,
            )
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

        elif advance_iter == order_windm_val:
            result = clubb_api.advance_windm_edsclrm(
                gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, edsclr_dim=edsclr_dim, dt=dt,
                wm_zt=wm_zt, kh_zm=Kh_zm, clubb_params=clubb_params, ug=ug, vg=vg,
                um_ref=um_ref, vm_ref=vm_ref, wp2=wp2, up2=up2, vp2=vp2,
                um_forcing=um_forcing, vm_forcing=vm_forcing,
                edsclrm_forcing=edsclrm_forcing, p_in_pa=p_in_Pa,
                rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
                invrs_rho_ds_zt=invrs_rho_ds_zt, fcor=fcor, l_implemented=l_implemented,
                ts_nudge=ts_nudge, tridiag_solve_method=flags.tridiag_solve_method,
                l_predict_upwp_vpwp=flags.l_predict_upwp_vpwp,
                l_upwind_xm_ma=flags.l_upwind_xm_ma, l_uv_nudge=flags.l_uv_nudge,
                l_tke_aniso=flags.l_tke_aniso, l_lmm_stepping=flags.l_lmm_stepping,
                l_linearize_pbl_winds=flags.l_linearize_pbl_winds,
                l_do_expldiff_rtm_thlm=flags.l_do_expldiff_rtm_thlm,
                fill_holes_type=flags.fill_holes_type,
                upwp_cl_num=upwp_cl_num, vpwp_cl_num=vpwp_cl_num,
                um=um, vm=vm, thlm=thlm, rtm=rtm, edsclrm=edsclrm,
                upwp=upwp, vpwp=vpwp, wpedsclrp=wpedsclrp, um_pert=um_pert,
                vm_pert=vm_pert, upwp_pert=upwp_pert, vpwp_pert=vpwp_pert,
                nu_vert_res_dep=nu_vert_res_dep, err_info=err_info,
            )
            (upwp_cl_num, vpwp_cl_num,
             um, vm, thlm, rtm, edsclrm, upwp, vpwp, _wpedsclrp_out,
             um_pert, vm_pert, upwp_pert, vpwp_pert, err_info) = result
            upwp_cl_num = int(upwp_cl_num)
            vpwp_cl_num = int(vpwp_cl_num)
            wpedsclrp = _wpedsclrp_out
            err_code = err_info.err_code
            if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
                return

    # Update local aliases after advance loop
    wp2 = wp2
    wp2_zt = np.maximum(
        clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=wp2),
        w_tol_sqd,
    )

    # ================================================================== #
    # Block T: Advance or diagnose third-order moments (xp3)
    # ================================================================== #
    if l_advance_xp3_flag and flags.iiPDF_type != iiPDF_ADG1:
        rtp3, thlp3, sclrp3, up3, vp3 = clubb_api.advance_xp3(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, sclr_dim=sclr_dim,
            sclr_tol=sclr_tol, dt=dt, rtm=rtm, thlm=thlm, rtp2=rtp2,
            thlp2=thlp2, wprtp=wprtp, wpthlp=wpthlp, wprtp2=wprtp2,
            wpthlp2=wpthlp2, rho_ds_zm=rho_ds_zm,
            invrs_rho_ds_zt=invrs_rho_ds_zt, invrs_tau_zt=invrs_tau_zt,
            tau_max_zt=tau_max_zt, sclrm=sclrm, sclrp2=sclrp2,
            wpsclrp=wpsclrp, wpsclrp2=_wpsclrp2,
            wp2=wp2, wp3=wp3, upwp=upwp, vpwp=vpwp, up2=up2, vp2=vp2,
            thvm=thvm, clubb_params=clubb_params,
            l_lmm_stepping=flags.l_lmm_stepping, rtp3=rtp3, thlp3=thlp3,
            sclrp3=sclrp3, up3=up3, vp3=vp3,
        )
    else:
        rtp3, thlp3, up3, vp3, sclrp3 = clubb_api.diagnose_xp3(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, sclr_dim=sclr_dim,
            sclr_tol=sclr_tol, iipdf_type=flags.iiPDF_type,
            clubb_params=clubb_params, wp2=wp2, wp3=wp3, thvm=thvm,
            wprtp=wprtp, wpthlp=wpthlp, rtp2=rtp2, thlp2=thlp2,
            upwp=upwp, vpwp=vpwp, up2=up2, vp2=vp2,
            sigma_sqd_w=sigma_sqd_w, wpsclrp=wpsclrp, sclrp2=sclrp2,
            rtp3=rtp3, thlp3=thlp3, up3=up3, vp3=vp3, sclrp3=sclrp3,
        )

    # ================================================================== #
    # Block U: Post-advance PDF closure (conditional)
    # ================================================================== #
    if (flags.ipdf_call_placement == ipdf_post_advance_fields
            or flags.ipdf_call_placement == ipdf_pre_post_advance_fields):

        if flags.ipdf_call_placement == ipdf_post_advance_fields:
            l_samp_stats = True
        else:
            l_samp_stats = False  # already sampled in pre-advance call

        pdf_result = clubb_api.pdf_closure_driver(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            dt=dt, hydromet_dim=hydromet_dim, sclr_dim=sclr_dim,
            sclr_tol=sclr_tol,
            wprtp=wprtp, thlm=thlm, wpthlp=wpthlp,
            rtp2=rtp2, rtp3=rtp3,
            thlp2=thlp2, thlp3=thlp3,
            rtpthlp=rtpthlp, wp2=wp2, wp3=wp3,
            wm_zm=wm_zm, wm_zt=wm_zt,
            um=um, up2=up2, upwp=upwp, up3=up3,
            vm=vm, vp2=vp2, vpwp=vpwp, vp3=vp3,
            p_in_pa=p_in_Pa, exner=exner,
            thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt,
            rtm_ref=rtm_ref,
            wphydrometp=wphydrometp,
            wp2hmp=wp2hmp,
            rtphmp_zt=rtphmp_zt,
            thlphmp_zt=thlphmp_zt,
            sclrm=sclrm, wpsclrp=wpsclrp,
            sclrp2=sclrp2, sclrprtp=sclrprtp,
            sclrpthlp=sclrpthlp, sclrp3=sclrp3,
            p_sfc=p_sfc,
            l_samp_stats_in_pdf_call=l_samp_stats,
            mixt_frac_max_mag=mixt_frac_max_mag,
            ts_nudge=ts_nudge,
            rtm_min=rtm_min,
            rtm_nudge_max_altitude=rtm_nudge_max_altitude,
            clubb_params=clubb_params,
            iiPDF_type=flags.iiPDF_type,
            saturation_formula=flags.saturation_formula,
            l_rtm_nudge=flags.l_rtm_nudge,
            l_trapezoidal_rule_zt=flags.l_trapezoidal_rule_zt,
            l_trapezoidal_rule_zm=flags.l_trapezoidal_rule_zm,
            l_call_pdf_closure_twice=flags.l_call_pdf_closure_twice,
            l_use_cloud_cover=flags.l_use_cloud_cover,
            l_rcm_supersat_adj=flags.l_rcm_supersat_adj,
            l_mix_rat_hm=l_mix_rat_hm,
            pdf_params=pdf_params,
            pdf_params_zm=pdf_params_zm,
            pdf_implicit_coefs_terms=pdf_implicit_coefs_terms,
            err_info=err_info,
            rtm=rtm,
        )
        (rtm, pdf_implicit_coefs_terms, pdf_params, pdf_params_zm, err_info,
         rcm, cloud_frac, ice_supersat_frac, wprcp_out, _sigma_sqd_w, wpthvp, wp2thvp,
         wp2up, rtpthvp, thlpthvp, _rc_coef, rcm_in_layer, cloud_cover,
         _rcp2_zt, thlprcp, rc_coef_zm, sclrpthvp, wpup2, wpvp2, wp2up2,
         wp2vp2, wp4, wp2rtp, _wprtp2, wp2thlp, _wpthlp2, _wprtpthlp, _wp2rcp,
         _rtprcp, _rcp2, uprcp, vprcp, w_up_in_cloud, w_down_in_cloud,
         cloudy_updraft_frac, cloudy_downdraft_frac, _skw_velocity,
         _cloud_frac_zm, _ice_supersat_frac_zm, _rtm_zm, _thlm_zm, _rcm_zm,
         _rcm_supersat_adj, _wp2sclrp, _wpsclrp2, _sclrprcp,
         _wpsclrprtp, _wpsclrpthlp) = pdf_result

        err_code = err_info.err_code
        if err_code is not None and np.any(np.asarray(err_code) == CLUBB_FATAL_ERROR):
            return

        # The post-advance PDF closure recomputes sigma_sqd_w and related
        # moments. Stats and later diagnostics should use the updated value.
        sigma_sqd_w = _sigma_sqd_w
        wprtp2 = _wprtp2
        wpthlp2 = _wpthlp2
        wprtpthlp = _wprtpthlp
        wp2rcp = _wp2rcp
        wp2sclrp = _wp2sclrp
        wpsclrp2 = _wpsclrp2
        wpsclrprtp = _wpsclrprtp
        wpsclrpthlp = _wpsclrpthlp

    # ================================================================== #
    # Block V: Stats — accumulate and finalize budgets
    # ================================================================== #
    if l_sample:
        # Recompute interpolated variables for stats
        wpthlp_zt = clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=wpthlp)
        wprtp_zt = clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=wprtp)
        up2_zt = np.maximum(
            clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=up2),
            w_tol_sqd,
        )
        vp2_zt = np.maximum(
            clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=vp2),
            w_tol_sqd,
        )
        upwp_zt = clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=upwp)
        vpwp_zt = clubb_api.zm2zt(gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, azm=vpwp)

        clubb_api.stats_accumulate(
            gr=gr, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
            dt=dt,
            l_implemented=l_implemented,
            l_host_applies_sfc_fluxes=flags.l_host_applies_sfc_fluxes,
            l_stability_correct_tau_zm=flags.l_stability_correct_tau_zm,
            clubb_params=clubb_params,
            um=um, vm=vm,
            upwp=upwp, vpwp=vpwp,
            up2=up2, vp2=vp2,
            thlm=thlm, rtm=rtm,
            thlm_before=thlm_before, rtm_before=rtm_before,
            thlm_forcing=thlm_forcing, rtm_forcing=rtm_forcing,
            wpthlp_sfc=wpthlp_sfc, wprtp_sfc=wprtp_sfc,
            wprtp=wprtp, wpthlp=wpthlp,
            wp2=wp2, wp3=wp3,
            rtp2=rtp2, rtp3=rtp3,
            thlp2=thlp2, thlp3=thlp3,
            rtpthlp=rtpthlp,
            p_in_pa=p_in_Pa, exner=exner,
            rho=rho, rho_zm=rho_zm,
            rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
            thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt,
            wm_zt=wm_zt, wm_zm=wm_zm,
            rcm=rcm,
            cloud_frac=cloud_frac,
            thvm=thvm,
            ug=ug, vg=vg,
            ddzt_umvm_sqd=ddzt_umvm_sqd,
            stability_correction=stability_correction,
            kh_zt=Kh_zt,
            rsat=rsat,
            kh_zm=Kh_zm,
            em=em,
            sclrm=sclrm,
            sclrp2=sclrp2,
            sclrprtp=sclrprtp,
            sclrpthlp=sclrpthlp,
            sclrm_forcing=sclrm_forcing,
            wpsclrp=wpsclrp,
            wpedsclrp=wpedsclrp,
            edsclrm=edsclrm,
            edsclrm_forcing=edsclrm_forcing,
            saturation_formula=flags.saturation_formula,
        )

        clubb_api.stats_finalize_budget("wp2_bt", wp2 / dt)
        clubb_api.stats_finalize_budget("vp2_bt", vp2 / dt)
        clubb_api.stats_finalize_budget("up2_bt", up2 / dt)
        clubb_api.stats_finalize_budget("wprtp_bt", wprtp / dt)
        clubb_api.stats_finalize_budget("wpthlp_bt", wpthlp / dt)
        if flags.l_predict_upwp_vpwp:
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

    if debug_level >= 2:
        err_info = clubb_api.parameterization_check(
            nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
            thlm_forcing=thlm_forcing, rtm_forcing=rtm_forcing,
            um_forcing=um_forcing, vm_forcing=vm_forcing,
            wm_zm=wm_zm, wm_zt=wm_zt, p_in_pa=p_in_Pa,
            rho_zm=rho_zm, rho=rho, exner=exner,
            rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
            invrs_rho_ds_zm=invrs_rho_ds_zm, invrs_rho_ds_zt=invrs_rho_ds_zt,
            thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt,
            wpthlp_sfc=wpthlp_sfc, wprtp_sfc=wprtp_sfc,
            upwp_sfc=upwp_sfc, vpwp_sfc=vpwp_sfc, p_sfc=p_sfc,
            um=um, upwp=upwp, vm=vm, vpwp=vpwp,
            up2=up2, vp2=vp2, rtm=rtm, wprtp=wprtp,
            thlm=thlm, wpthlp=wpthlp, wp2=wp2, wp3=wp3,
            rtp2=rtp2, thlp2=thlp2, rtpthlp=rtpthlp,
            prefix="end of ", wpsclrp_sfc=wpsclrp_sfc, wpedsclrp_sfc=wpedsclrp_sfc,
            sclrm=sclrm, wpsclrp=wpsclrp, sclrp2=sclrp2,
            sclrprtp=sclrprtp, sclrpthlp=sclrpthlp,
            sclrm_forcing=sclrm_forcing, edsclrm=edsclrm,
            edsclrm_forcing=edsclrm_forcing, err_info=err_info,
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
        _sigma_sqd_w,
        _rc_coef, _rcp2_zt, _wprtp2,
        _wpthlp2, _wprtpthlp, _wp2rcp,
        _rtprcp, _rcp2, _skw_velocity,
        _cloud_frac_zm, _ice_supersat_frac_zm,
        _rtm_zm, _thlm_zm, _rcm_zm,
        _rcm_supersat_adj, _wp2sclrp,
        _wpsclrp2, _sclrprcp,
        _wpsclrprtp, _wpsclrpthlp,
    )

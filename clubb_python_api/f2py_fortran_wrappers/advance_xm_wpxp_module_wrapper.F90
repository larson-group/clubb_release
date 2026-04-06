! advance_xm_wpxp_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_advance_xm_wpxp(nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport, sclr_tol, dt, &
    sigma_sqd_w, wm_zm, wm_zt, wp2, lscale_zm, &
    wp3_on_wp2, wp3_on_wp2_zt, kh_zt, kh_zm, &
    stability_correction, &
    invrs_tau_c6_zm, tau_max_zm, skw_zm, wp2rtp, rtpthvp, &
    rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp, &
    thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
    invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, &
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
    mixt_frac_zm, l_implemented, em, wp2sclrp, &
    sclrpthvp, sclrm_forcing, sclrp2, cx_fnc_richardson, &
    um_forcing, vm_forcing, ug, vg, wpthvp, &
    fcor, fcor_y, um_ref, vm_ref, up2, vp2, &
    uprcp, vprcp, rc_coef_zm, &
    clubb_params, ts_nudge, &
    iipdf_type, &
    penta_solve_method, &
    tridiag_solve_method, &
    fill_holes_type, &
    l_predict_upwp_vpwp, &
    l_ho_nontrad_coriolis, &
    l_ho_trad_coriolis, &
    l_diffuse_rtm_and_thlm, &
    l_stability_correct_kh_n2_zm, &
    l_godunov_upwind_wpxp_ta, &
    l_upwind_xm_ma, &
    l_uv_nudge, &
    l_tke_aniso, &
    l_diag_lscale_from_tau, &
    l_use_c7_richardson, &
    l_lmm_stepping, &
    l_enable_relaxed_clipping, &
    l_linearize_pbl_winds, &
    l_mono_flux_lim_thlm, &
    l_mono_flux_lim_rtm, &
    l_mono_flux_lim_um, &
    l_mono_flux_lim_vm, &
    l_mono_flux_lim_spikefix, &
    order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, &
    rtm, wprtp, thlm, wpthlp, &
    sclrm, wpsclrp, um, upwp, vm, vpwp, &
    um_pert, vm_pert, upwp_pert, vpwp_pert)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_pdf_implicit_coefs_terms, stored_nu_vert_res_dep, &
    stored_stats, stored_err_info
  use advance_xm_wpxp_module, only: advance_xm_wpxp

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclr_tol
  real(core_rknd), intent(in) :: dt

  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    wm_zt, wp3_on_wp2_zt, kh_zt, wp2rtp, rtm_forcing, rtm_ref, wp2thlp, &
    thlm_forcing, thlm_ref, rho_ds_zt, invrs_rho_ds_zt, &
    um_forcing, vm_forcing, ug, vg, um_ref, vm_ref

  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    sigma_sqd_w, wm_zm, wp2, lscale_zm, wp3_on_wp2, kh_zm, stability_correction, &
    invrs_tau_c6_zm, tau_max_zm, skw_zm, rtpthvp, wprtp_forcing, thlpthvp, &
    wpthlp_forcing, rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, rtp2, thlp2, &
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, em, &
    cx_fnc_richardson, wpthvp, up2, vp2, uprcp, vprcp, rc_coef_zm

  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(in) :: &
    wp2sclrp, sclrm_forcing
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(in) :: &
    sclrpthvp, sclrp2

  real(core_rknd), dimension(ngrdcol), intent(in) :: fcor, fcor_y
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), intent(in) :: ts_nudge

  integer, intent(in) :: iipdf_type, penta_solve_method, tridiag_solve_method, fill_holes_type
  logical, intent(in) :: &
    l_implemented, l_predict_upwp_vpwp, l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
    l_diffuse_rtm_and_thlm, l_stability_correct_kh_n2_zm, l_godunov_upwind_wpxp_ta, &
    l_upwind_xm_ma, l_uv_nudge, l_tke_aniso, l_diag_lscale_from_tau, l_use_c7_richardson, &
    l_lmm_stepping, l_enable_relaxed_clipping, l_linearize_pbl_winds, l_mono_flux_lim_thlm, &
    l_mono_flux_lim_rtm, l_mono_flux_lim_um, l_mono_flux_lim_vm, l_mono_flux_lim_spikefix
  integer, intent(in) :: order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: &
    rtm, thlm, um, vm, um_pert, vm_pert
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    wprtp, wpthlp, upwp, vpwp, upwp_pert, vpwp_pert
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(inout) :: sclrm
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: wpsclrp

  call advance_xm_wpxp(nzm, nzt, ngrdcol, sclr_dim, sclr_tol, stored_grid, dt, &
    sigma_sqd_w, wm_zm, wm_zt, wp2, lscale_zm, &
    wp3_on_wp2, wp3_on_wp2_zt, kh_zt, kh_zm, &
    stability_correction, &
    invrs_tau_c6_zm, tau_max_zm, skw_zm, wp2rtp, rtpthvp, &
    rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp, &
    thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
    invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, &
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
    mixt_frac_zm, l_implemented, em, wp2sclrp, &
    sclrpthvp, sclrm_forcing, sclrp2, cx_fnc_richardson, &
    stored_pdf_implicit_coefs_terms, &
    um_forcing, vm_forcing, ug, vg, wpthvp, &
    fcor, fcor_y, um_ref, vm_ref, up2, vp2, &
    uprcp, vprcp, rc_coef_zm, &
    clubb_params, stored_nu_vert_res_dep, ts_nudge, &
    iipdf_type, &
    penta_solve_method, &
    tridiag_solve_method, &
    fill_holes_type, &
    l_predict_upwp_vpwp, &
    l_ho_nontrad_coriolis, &
    l_ho_trad_coriolis, &
    l_diffuse_rtm_and_thlm, &
    l_stability_correct_kh_n2_zm, &
    l_godunov_upwind_wpxp_ta, &
    l_upwind_xm_ma, &
    l_uv_nudge, &
    l_tke_aniso, &
    l_diag_lscale_from_tau, &
    l_use_c7_richardson, &
    l_lmm_stepping, &
    l_enable_relaxed_clipping, &
    l_linearize_pbl_winds, &
    l_mono_flux_lim_thlm, &
    l_mono_flux_lim_rtm, &
    l_mono_flux_lim_um, &
    l_mono_flux_lim_vm, &
    l_mono_flux_lim_spikefix, &
    order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, &
    stored_stats, &
    rtm, wprtp, thlm, wpthlp, &
    sclrm, wpsclrp, um, upwp, vm, vpwp, &
    um_pert, vm_pert, upwp_pert, vpwp_pert, stored_err_info)

end subroutine f2py_advance_xm_wpxp

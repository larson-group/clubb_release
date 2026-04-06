! advance_clubb_core_module_wrapper.F90 — wrappers for advance_clubb_core_module APIs

subroutine f2py_advance_clubb_core( &
    nzm, nzt, ngrdcol, &
    l_implemented, dt, &
    fcor, fcor_y, sfc_elevation, &
    hydromet_dim, sclr_dim, edsclr_dim, &
    hydromet_dim_transport, sclr_dim_transport, edsclr_dim_transport, &
    sclr_tol, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
    sclrm_forcing, edsclrm_forcing, &
    wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    rtpthlp_forcing, wm_zm, wm_zt, &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
    wpsclrp_sfc, wpedsclrp_sfc, &
    upwp_sfc_pert, vpwp_sfc_pert, &
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
    rho_zm, rho, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    l_mix_rat_hm, &
    rfrzm, &
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
    host_dx, host_dy, &
    clubb_params, lmin, mixt_frac_max_mag, t0_val, ts_nudge, &
    rtm_min, rtm_nudge_max_altitude, &
    um, vm, up3, vp3, thlm, rtm, rtp3, thlp3, wp3, &
    upwp, vpwp, up2, vp2, wprtp, wpthlp, rtp2, thlp2, rtpthlp, wp2, &
    sclrm, sclrp3, &
    wpsclrp, sclrp2, sclrprtp, sclrpthlp, &
    p_in_Pa, exner, rcm, cloud_frac, wp2thvp, wp2up, &
    wpthvp, rtpthvp, thlpthvp, &
    sclrpthvp, &
    wp2rtp, wp2thlp, wpup2, wpvp2, ice_supersat_frac, &
    uprcp, vprcp, rc_coef_zm, wp4, wp2up2, wp2vp2, &
    um_pert, vm_pert, upwp_pert, vpwp_pert, &
    edsclrm, &
    rcm_in_layer, cloud_cover, w_up_in_cloud, w_down_in_cloud, &
    cloudy_updraft_frac, cloudy_downdraft_frac, wprcp, invrs_tau_zm, &
    Kh_zt, Kh_zm, thlprcp, Lscale)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_grid, stored_sclr_idx, &
      stored_config_flags, stored_nu_vert_res_dep, &
      stored_stats, stored_pdf_params, stored_pdf_params_zm, &
      stored_pdf_implicit_coefs_terms, stored_err_info
  use advance_clubb_core_module, only: advance_clubb_core

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  integer, intent(in) :: hydromet_dim, sclr_dim, edsclr_dim
  integer, intent(in) :: hydromet_dim_transport, sclr_dim_transport, edsclr_dim_transport

  logical, intent(in) :: l_implemented
  real(core_rknd), intent(in) :: dt
  real(core_rknd), intent(in) :: lmin, mixt_frac_max_mag, t0_val, ts_nudge
  real(core_rknd), intent(in) :: rtm_min, rtm_nudge_max_altitude

  real(core_rknd), dimension(ngrdcol), intent(in) :: &
    fcor, fcor_y, sfc_elevation, &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
    upwp_sfc_pert, vpwp_sfc_pert, host_dx, host_dy

  real(core_rknd), dimension(sclr_dim), intent(in) :: sclr_tol
  logical, dimension(hydromet_dim_transport), intent(in) :: l_mix_rat_hm

  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
    wm_zt, rho, rho_ds_zt, invrs_rho_ds_zt, thv_ds_zt, rfrzm, &
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg

  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    rtpthlp_forcing, wm_zm, rho_zm, rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm

  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params

  real(core_rknd), dimension(ngrdcol, nzm, hydromet_dim_transport), intent(in) :: wphydrometp
  real(core_rknd), dimension(ngrdcol, nzt, hydromet_dim_transport), intent(in) :: &
    wp2hmp, rtphmp_zt, thlphmp_zt

  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(in) :: sclrm_forcing
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(in) :: edsclrm_forcing
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(in) :: wpsclrp_sfc
  real(core_rknd), dimension(ngrdcol, edsclr_dim_transport), intent(in) :: wpedsclrp_sfc

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: &
    um, vm, up3, vp3, thlm, rtm, rtp3, thlp3, wp3

  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    upwp, vpwp, up2, vp2, wprtp, wpthlp, rtp2, thlp2, rtpthlp, wp2

  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(inout) :: &
    sclrm, sclrp3

  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: &
    wpsclrp, sclrp2, sclrprtp, sclrpthlp

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: &
    p_in_Pa, exner, rcm, cloud_frac, wp2thvp, wp2up

  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    wpthvp, rtpthvp, thlpthvp

  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: &
    sclrpthvp

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: &
    wp2rtp, wp2thlp, wpup2, wpvp2, ice_supersat_frac

  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    uprcp, vprcp, rc_coef_zm, wp4, wp2up2, wp2vp2

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: &
    um_pert, vm_pert

  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    upwp_pert, vpwp_pert

  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(inout) :: &
    edsclrm

  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: &
    rcm_in_layer, cloud_cover, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
    cloudy_downdraft_frac, Kh_zt, Lscale

  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: &
    wprcp, invrs_tau_zm, Kh_zm, thlprcp

  call advance_clubb_core( &
    stored_grid, nzm, nzt, ngrdcol, &
    l_implemented, dt, fcor, fcor_y, sfc_elevation, &
    hydromet_dim, &
    sclr_dim, sclr_tol, edsclr_dim, stored_sclr_idx, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
    sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
    wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    rtpthlp_forcing, wm_zm, wm_zt, &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
    wpsclrp_sfc, wpedsclrp_sfc, &
    upwp_sfc_pert, vpwp_sfc_pert, &
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
    p_in_Pa, rho_zm, rho, exner, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
    invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    l_mix_rat_hm(1:hydromet_dim), &
    rfrzm, &
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
    host_dx, host_dy, &
    clubb_params, stored_nu_vert_res_dep, lmin, &
    mixt_frac_max_mag, t0_val, ts_nudge, &
    rtm_min, rtm_nudge_max_altitude, &
    stored_config_flags, &
    stored_stats, &
    um, vm, upwp, vpwp, up2, vp2, up3, vp3, &
    thlm, rtm, wprtp, wpthlp, &
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
    sclrm, &
    sclrp2, sclrp3, sclrprtp, sclrpthlp, &
    wpsclrp, edsclrm, &
    rcm, cloud_frac, &
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
    sclrpthvp, &
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &
    um_pert, vm_pert, upwp_pert, vpwp_pert, &
    stored_pdf_params, stored_pdf_params_zm, &
    stored_pdf_implicit_coefs_terms, &
    stored_err_info, &
    Kh_zm, Kh_zt, &
    thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &
    cloudy_updraft_frac, cloudy_downdraft_frac, &
    rcm_in_layer, cloud_cover, invrs_tau_zm, &
    Lscale)

end subroutine f2py_advance_clubb_core

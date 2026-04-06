! stats_clubb_utilities_wrapper.F90 — wrapper for module stats_clubb_utilities

subroutine f2py_stats_accumulate( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport, &
    invrs_dzm, zt, dzm, dzt, dt, &
    um, vm, upwp, vpwp, up2, vp2, &
    thlm, rtm, wprtp, wpthlp, &
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
    rtpthlp, &
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
    p_in_Pa, exner, rho, rho_zm, &
    rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    wm_zt, wm_zm, rcm, wprcp, rc_coef, &
    rc_coef_zm, &
    rcm_zm, rtm_zm, thlm_zm, cloud_frac, &
    ice_supersat_frac, &
    cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
    cloud_cover, rcm_supersat_adj, sigma_sqd_w, &
    thvm, ug, vg, Lscale, wpthlp2, wp2thlp, &
    wprtp2, wp2rtp, &
    Lscale_up, Lscale_down, tau_zt, Kh_zt, wp2rcp, &
    wprtpthlp, sigma_sqd_w_zt, rsat, wp2_zt, &
    thlp2_zt, &
    wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, &
    up2_zt, &
    vp2_zt, upwp_zt, vpwp_zt, wpup2, wpvp2, &
    wp2up2, wp2vp2, wp4, &
    tau_zm, Kh_zm, thlprcp, &
    rtprcp, rcp2, em, a3_coef, a3_coef_zt, &
    wp3_zm, wp3_on_wp2, wp3_on_wp2_zt, Skw_velocity, &
    w_up_in_cloud, w_down_in_cloud, &
    cloudy_updraft_frac, cloudy_downdraft_frac, &
    sclrm, sclrp2, &
    sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp, &
    wpsclrp, sclrprcp, wp2sclrp, wpsclrp2, &
    wpsclrprtp, &
    wpsclrpthlp, wpedsclrp, edsclrm, &
    edsclrm_forcing, &
    saturation_formula, &
    l_call_pdf_closure_twice)

  use clubb_precision, only: core_rknd
  use stats_clubb_utilities, only: stats_accumulate
  use derived_type_storage, only: stored_pdf_params, stored_pdf_params_zm, stored_stats

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: invrs_dzm, dzm
  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: zt, dzt
  real(core_rknd), intent(in) :: dt

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    um, vm, thlm, rtm, wp3, rtp3, thlp3, wp2thvp, wp2up

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    upwp, vpwp, up2, vp2, wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp, &
    wpthvp, rtpthvp, thlpthvp

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    p_in_Pa, exner, rho, rho_ds_zt, thv_ds_zt, wm_zt

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    rho_zm, rho_ds_zm, thv_ds_zm, wm_zm

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    rcm, rc_coef, cloud_frac, ice_supersat_frac, rcm_in_layer, cloud_cover, rcm_supersat_adj

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    rcm_zm, rtm_zm, thlm_zm, wprcp, rc_coef_zm, cloud_frac_zm, ice_supersat_frac_zm

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: sigma_sqd_w

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    thvm, ug, vg, Lscale, wpthlp2, wp2thlp, wprtp2, wp2rtp, Lscale_up, Lscale_down, tau_zt, &
    Kh_zt, wp2rcp, wprtpthlp, sigma_sqd_w_zt, rsat

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    wp2_zt, thlp2_zt, wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, up2_zt, vp2_zt, upwp_zt, vpwp_zt, &
    wpup2, wpvp2, a3_coef_zt, wp3_on_wp2_zt, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
    cloudy_downdraft_frac

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    wp2up2, wp2vp2, wp4, tau_zm, Kh_zm, thlprcp, rtprcp, rcp2, em, a3_coef, wp3_zm, wp3_on_wp2, &
    Skw_velocity

  real(core_rknd), intent(in), dimension(ngrdcol, nzt, sclr_dim_transport) :: &
    sclrm, sclrm_forcing, wp2sclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp

  real(core_rknd), intent(in), dimension(ngrdcol, nzm, sclr_dim_transport) :: &
    sclrprcp, sclrp2, sclrprtp, sclrpthlp, sclrpthvp, wpsclrp

  real(core_rknd), intent(in), dimension(ngrdcol, nzt, edsclr_dim_transport) :: &
    edsclrm, edsclrm_forcing

  real(core_rknd), intent(in), dimension(ngrdcol, nzm, edsclr_dim_transport) :: &
    wpedsclrp

  integer, intent(in) :: saturation_formula
  logical, intent(in) :: l_call_pdf_closure_twice

  call stats_accumulate( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, &
    invrs_dzm, zt, dzm, dzt, dt, &
    um, vm, upwp, vpwp, up2, vp2, &
    thlm, rtm, wprtp, wpthlp, &
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
    rtpthlp, &
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
    p_in_Pa, exner, rho, rho_zm, &
    rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    wm_zt, wm_zm, rcm, wprcp, rc_coef, &
    rc_coef_zm, &
    rcm_zm, rtm_zm, thlm_zm, cloud_frac, &
    ice_supersat_frac, &
    cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
    cloud_cover, rcm_supersat_adj, sigma_sqd_w, &
    thvm, ug, vg, Lscale, wpthlp2, wp2thlp, &
    wprtp2, wp2rtp, &
    Lscale_up, Lscale_down, tau_zt, Kh_zt, wp2rcp, &
    wprtpthlp, sigma_sqd_w_zt, rsat, wp2_zt, &
    thlp2_zt, &
    wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, &
    up2_zt, &
    vp2_zt, upwp_zt, vpwp_zt, wpup2, wpvp2, &
    wp2up2, wp2vp2, wp4, &
    tau_zm, Kh_zm, thlprcp, &
    rtprcp, rcp2, em, a3_coef, a3_coef_zt, &
    wp3_zm, wp3_on_wp2, wp3_on_wp2_zt, Skw_velocity, &
    w_up_in_cloud, w_down_in_cloud, &
    cloudy_updraft_frac, cloudy_downdraft_frac, &
    stored_pdf_params, stored_pdf_params_zm, &
    sclrm, sclrp2, &
    sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp, &
    wpsclrp, sclrprcp, wp2sclrp, wpsclrp2, &
    wpsclrprtp, &
    wpsclrpthlp, wpedsclrp, edsclrm, &
    edsclrm_forcing, &
    saturation_formula, &
    l_call_pdf_closure_twice, &
    stored_stats )

end subroutine f2py_stats_accumulate

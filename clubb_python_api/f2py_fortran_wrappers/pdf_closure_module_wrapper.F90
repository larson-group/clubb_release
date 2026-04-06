! pdf_closure_module_wrapper.F90 — wrappers for pdf_closure_module APIs

subroutine f2py_calc_wp4_pdf(nz, ngrdcol, wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac, wp4)

  use clubb_precision, only: core_rknd
  use pdf_closure_module, only: calc_wp4_pdf

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: wp4

  call calc_wp4_pdf(nz, ngrdcol, wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac, wp4)

end subroutine f2py_calc_wp4_pdf

subroutine f2py_calc_wp2xp_pdf(nz, ngrdcol, wm, xm, w_1, w_2, x_1, x_2, &
    varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac, wp2xp)

  use clubb_precision, only: core_rknd
  use pdf_closure_module, only: calc_wp2xp_pdf

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, &
    corr_w_x_1, corr_w_x_2, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: wp2xp

  call calc_wp2xp_pdf(nz, ngrdcol, wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2, &
    varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac, wp2xp)

end subroutine f2py_calc_wp2xp_pdf

subroutine f2py_calc_wpxp2_pdf(nz, ngrdcol, wm, xm, w_1, w_2, x_1, x_2, &
    varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac, wpxp2)

  use clubb_precision, only: core_rknd
  use pdf_closure_module, only: calc_wpxp2_pdf

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, &
    corr_w_x_1, corr_w_x_2, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: wpxp2

  call calc_wpxp2_pdf(nz, ngrdcol, wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2, &
    varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac, wpxp2)

end subroutine f2py_calc_wpxp2_pdf

subroutine f2py_calc_wpxpyp_pdf(nz, ngrdcol, wm, xm, ym, w_1, w_2, x_1, x_2, y_1, y_2, &
    varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, varnce_y_1, varnce_y_2, &
    corr_w_x_1, corr_w_x_2, corr_w_y_1, corr_w_y_2, corr_x_y_1, corr_x_y_2, mixt_frac, wpxpyp)

  use clubb_precision, only: core_rknd
  use pdf_closure_module, only: calc_wpxpyp_pdf

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, xm, ym, w_1, w_2, x_1, x_2, y_1, y_2, varnce_w_1, varnce_w_2, &
    varnce_x_1, varnce_x_2, varnce_y_1, varnce_y_2, corr_w_x_1, corr_w_x_2, &
    corr_w_y_1, corr_w_y_2, corr_x_y_1, corr_x_y_2, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: wpxpyp

  call calc_wpxpyp_pdf(nz, ngrdcol, wm, xm, ym, w_1, w_2, x_1, x_2, y_1, y_2, &
    varnce_w_1, varnce_w_2, varnce_x_1, varnce_x_2, varnce_y_1, varnce_y_2, &
    corr_w_x_1, corr_w_x_2, corr_w_y_1, corr_w_y_2, corr_x_y_1, corr_x_y_2, mixt_frac, wpxpyp)

end subroutine f2py_calc_wpxpyp_pdf

subroutine f2py_calc_w_up_in_cloud(nz, ngrdcol, mixt_frac, cloud_frac_1, cloud_frac_2, &
    w_1, w_2, varnce_w_1, varnce_w_2, w_up_in_cloud, w_down_in_cloud, &
    cloudy_updraft_frac, cloudy_downdraft_frac)

  use clubb_precision, only: core_rknd
  use pdf_closure_module, only: calc_w_up_in_cloud

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    mixt_frac, cloud_frac_1, cloud_frac_2, w_1, w_2, varnce_w_1, varnce_w_2
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, cloudy_downdraft_frac

  call calc_w_up_in_cloud(nz, ngrdcol, mixt_frac, cloud_frac_1, cloud_frac_2, w_1, w_2, &
    varnce_w_1, varnce_w_2, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
    cloudy_downdraft_frac)

end subroutine f2py_calc_w_up_in_cloud

subroutine f2py_pdf_closure_driver(nzm, nzt, ngrdcol, dt, hydromet_dim, sclr_dim, hydromet_dim_transport, sclr_dim_transport, &
    sclr_tol, wprtp, thlm, wpthlp, rtp2, rtp3, thlp2, thlp3, rtpthlp, wp2, wp3, &
    wm_zm, wm_zt, um, up2, upwp, up3, vm, vp2, vpwp, vp3, p_in_pa, exner, thv_ds_zm, &
    thv_ds_zt, rtm_ref, wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
    sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp, sclrp3, p_sfc, &
    l_samp_stats_in_pdf_call, mixt_frac_max_mag, ts_nudge, rtm_min, rtm_nudge_max_altitude, &
    clubb_params, iipdf_type, saturation_formula, l_predict_upwp_vpwp, l_rtm_nudge, &
    l_trapezoidal_rule_zt, l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
    l_use_cloud_cover, l_rcm_supersat_adj, l_mix_rat_hm, rtm, &
    rcm, cloud_frac, ice_supersat_frac, wprcp, sigma_sqd_w, wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
    rc_coef, rcm_in_layer, cloud_cover, rcp2_zt, thlprcp, rc_coef_zm, sclrpthvp, &
    wpup2, wpvp2, wp2up2, wp2vp2, wp4, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, wp2rcp, &
    rtprcp, rcp2, uprcp, vprcp, w_up_in_cloud, w_down_in_cloud, &
    cloudy_updraft_frac, cloudy_downdraft_frac, skw_velocity, cloud_frac_zm, &
    ice_supersat_frac_zm, rtm_zm, thlm_zm, rcm_zm, rcm_supersat_adj, &
    wp2sclrp, wpsclrp2, sclrprcp, wpsclrprtp, wpsclrpthlp)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: &
    stored_grid, stored_err_info, stored_pdf_params, stored_pdf_params_zm, &
    stored_pdf_implicit_coefs_terms, stored_stats
  use pdf_closure_module, only: pdf_closure_driver
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, hydromet_dim
  integer, intent(in) :: hydromet_dim_transport, sclr_dim_transport
  real(core_rknd), intent(in) :: dt
  integer, intent(in) :: saturation_formula, iipdf_type
  logical, intent(in) :: l_predict_upwp_vpwp, l_call_pdf_closure_twice
  logical, intent(in) :: l_samp_stats_in_pdf_call, l_rtm_nudge
  logical, intent(in) :: l_trapezoidal_rule_zt, l_trapezoidal_rule_zm
  logical, intent(in) :: l_use_cloud_cover, l_rcm_supersat_adj
  logical, dimension(hydromet_dim_transport), intent(in) :: l_mix_rat_hm

  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclr_tol
  real(core_rknd), intent(in) :: mixt_frac_max_mag, ts_nudge, rtm_min, rtm_nudge_max_altitude
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params

  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    thlm, rtp3, thlp3, wp3, wm_zt, p_in_pa, exner, thv_ds_zt, rtm_ref, &
    um, up3, vm, vp3

  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    wprtp, wpthlp, rtp2, thlp2, rtpthlp, wp2, wm_zm, thv_ds_zm, &
    up2, upwp, vp2, vpwp

  real(core_rknd), dimension(ngrdcol, nzm, hydromet_dim_transport), intent(in) :: wphydrometp
  real(core_rknd), dimension(ngrdcol, nzt, hydromet_dim_transport), intent(in) :: &
    wp2hmp, rtphmp_zt, thlphmp_zt

  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(in) :: sclrm, sclrp3
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(in) :: &
    wpsclrp, sclrp2, sclrprtp, sclrpthlp

  real(core_rknd), dimension(ngrdcol), intent(in) :: p_sfc
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: rtm

  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: &
    rcm, cloud_frac, ice_supersat_frac, wp2thvp, wp2up, rc_coef, rcm_in_layer, cloud_cover, &
    rcp2_zt, wpup2, wpvp2, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, wp2rcp, &
    w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, cloudy_downdraft_frac, rcm_supersat_adj
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: &
    wprcp, sigma_sqd_w, wpthvp, rtpthvp, thlpthvp, thlprcp, rc_coef_zm, &
    wp2up2, wp2vp2, wp4, rtprcp, rcp2, uprcp, vprcp, &
    skw_velocity, cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, thlm_zm, rcm_zm
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(out) :: sclrpthvp, sclrprcp
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(out) :: &
    wp2sclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp

  call pdf_closure_driver( &
    gr=stored_grid, nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, &
    dt=dt, hydromet_dim=hydromet_dim, sclr_dim=sclr_dim, sclr_tol=sclr_tol, &
    wprtp=wprtp, thlm=thlm, wpthlp=wpthlp, rtp2=rtp2, rtp3=rtp3, &
    thlp2=thlp2, thlp3=thlp3, rtpthlp=rtpthlp, wp2=wp2, wp3=wp3, wm_zm=wm_zm, wm_zt=wm_zt, &
    um=um, up2=up2, upwp=upwp, up3=up3, vm=vm, vp2=vp2, vpwp=vpwp, vp3=vp3, &
    p_in_pa=p_in_pa, exner=exner, thv_ds_zm=thv_ds_zm, thv_ds_zt=thv_ds_zt, rtm_ref=rtm_ref, &
    wphydrometp=wphydrometp, wp2hmp=wp2hmp, rtphmp_zt=rtphmp_zt, thlphmp_zt=thlphmp_zt, &
    sclrm=sclrm, wpsclrp=wpsclrp, sclrp2=sclrp2, sclrprtp=sclrprtp, sclrpthlp=sclrpthlp, sclrp3=sclrp3, &
    p_sfc=p_sfc, l_samp_stats_in_pdf_call=l_samp_stats_in_pdf_call, mixt_frac_max_mag=mixt_frac_max_mag, &
    ts_nudge=ts_nudge, rtm_min=rtm_min, rtm_nudge_max_altitude=rtm_nudge_max_altitude, &
    clubb_params=clubb_params, iiPDF_type=iipdf_type, saturation_formula=saturation_formula, &
    l_predict_upwp_vpwp=l_predict_upwp_vpwp, l_rtm_nudge=l_rtm_nudge, &
    l_trapezoidal_rule_zt=l_trapezoidal_rule_zt, l_trapezoidal_rule_zm=l_trapezoidal_rule_zm, &
    l_call_pdf_closure_twice=l_call_pdf_closure_twice, l_use_cloud_cover=l_use_cloud_cover, &
    l_rcm_supersat_adj=l_rcm_supersat_adj, l_mix_rat_hm=l_mix_rat_hm, &
    stats=stored_stats, rtm=rtm, pdf_implicit_coefs_terms=stored_pdf_implicit_coefs_terms, &
    pdf_params=stored_pdf_params, pdf_params_zm=stored_pdf_params_zm, err_info=stored_err_info, &
    rcm=rcm, cloud_frac=cloud_frac, ice_supersat_frac=ice_supersat_frac, wprcp=wprcp, &
    sigma_sqd_w=sigma_sqd_w, wpthvp=wpthvp, wp2thvp=wp2thvp, wp2up=wp2up, &
    rtpthvp=rtpthvp, thlpthvp=thlpthvp, rc_coef=rc_coef, rcm_in_layer=rcm_in_layer, cloud_cover=cloud_cover, &
    rcp2_zt=rcp2_zt, thlprcp=thlprcp, rc_coef_zm=rc_coef_zm, sclrpthvp=sclrpthvp, &
    wpup2=wpup2, wpvp2=wpvp2, wp2up2=wp2up2, wp2vp2=wp2vp2, wp4=wp4, wp2rtp=wp2rtp, &
    wprtp2=wprtp2, wp2thlp=wp2thlp, wpthlp2=wpthlp2, wprtpthlp=wprtpthlp, wp2rcp=wp2rcp, &
    rtprcp=rtprcp, rcp2=rcp2, uprcp=uprcp, vprcp=vprcp, &
    w_up_in_cloud=w_up_in_cloud, w_down_in_cloud=w_down_in_cloud, &
    cloudy_updraft_frac=cloudy_updraft_frac, cloudy_downdraft_frac=cloudy_downdraft_frac, &
    skw_velocity=skw_velocity, cloud_frac_zm=cloud_frac_zm, ice_supersat_frac_zm=ice_supersat_frac_zm, &
    rtm_zm=rtm_zm, thlm_zm=thlm_zm, rcm_zm=rcm_zm, rcm_supersat_adj=rcm_supersat_adj, &
    wp2sclrp=wp2sclrp, wpsclrp2=wpsclrp2, sclrprcp=sclrprcp, wpsclrprtp=wpsclrprtp, wpsclrpthlp=wpsclrpthlp )

end subroutine f2py_pdf_closure_driver

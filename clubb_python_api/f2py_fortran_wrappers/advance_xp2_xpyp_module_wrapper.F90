! advance_xp2_xpyp_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_advance_xp2_xpyp(nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport, sclr_tol, &
    invrs_tau_xp2_zm, invrs_tau_c4_zm, invrs_tau_c14_zm, wm_zm, &
    rtm, wprtp, thlm, wpthlp, wpthvp, um, vm, &
    wp2, wp2_zt, wp3, upwp, vpwp, &
    sigma_sqd_w, wprtp2, wpthlp2, wprtpthlp, kh_zt, &
    rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, thv_ds_zm, cloud_frac, &
    wp3_on_wp2, wp3_on_wp2_zt, dt, fcor_y, &
    sclrm, wpsclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp, lhs_splat_wp2, &
    clubb_params, iipdf_type, tridiag_solve_method, fill_holes_type, &
    l_predict_upwp_vpwp, l_ho_nontrad_coriolis, l_min_xp2_from_corr_wx, &
    l_c2_cloud_frac, l_upwind_xpyp_ta, l_godunov_upwind_xpyp_ta, l_lmm_stepping, &
    rtp2, thlp2, rtpthlp, up2, vp2, sclrp2, sclrprtp, sclrpthlp)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_sclr_idx, stored_pdf_implicit_coefs_terms, &
    stored_nu_vert_res_dep, stored_stats, stored_err_info
  use advance_xp2_xpyp_module, only: advance_xp2_xpyp

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclr_tol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    invrs_tau_xp2_zm, invrs_tau_c4_zm, invrs_tau_c14_zm, wm_zm, wprtp, wpthlp, wpthvp, &
    wp2, upwp, vpwp, sigma_sqd_w, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
    rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, wp3_on_wp2, lhs_splat_wp2
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    rtm, thlm, um, vm, wp2_zt, wp3, wprtp2, wpthlp2, wprtpthlp, kh_zt, &
    rho_ds_zt, cloud_frac, wp3_on_wp2_zt
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol), intent(in) :: fcor_y
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(in) :: &
    sclrm, wpsclrp2, wpsclrprtp, wpsclrpthlp
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(in) :: wpsclrp
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  integer, intent(in) :: iipdf_type, tridiag_solve_method, fill_holes_type
  logical, intent(in) :: &
    l_predict_upwp_vpwp, l_ho_nontrad_coriolis, l_min_xp2_from_corr_wx, l_c2_cloud_frac, &
    l_upwind_xpyp_ta, l_godunov_upwind_xpyp_ta, l_lmm_stepping

  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: rtp2, thlp2, rtpthlp, up2, vp2
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: &
    sclrp2, sclrprtp, sclrpthlp

  call advance_xp2_xpyp(nzm, nzt, ngrdcol, sclr_dim, sclr_tol, stored_grid, stored_sclr_idx, &
    invrs_tau_xp2_zm, invrs_tau_c4_zm, invrs_tau_c14_zm, wm_zm, &
    rtm, wprtp, thlm, wpthlp, wpthvp, um, vm, &
    wp2, wp2_zt, wp3, upwp, vpwp, &
    sigma_sqd_w, wprtp2, wpthlp2, wprtpthlp, kh_zt, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, thv_ds_zm, cloud_frac, &
    wp3_on_wp2, wp3_on_wp2_zt, stored_pdf_implicit_coefs_terms, dt, fcor_y, &
    sclrm, wpsclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp, lhs_splat_wp2, &
    clubb_params, stored_nu_vert_res_dep, iipdf_type, tridiag_solve_method, fill_holes_type, &
    l_predict_upwp_vpwp, l_ho_nontrad_coriolis, &
    l_min_xp2_from_corr_wx, l_c2_cloud_frac, &
    l_upwind_xpyp_ta, l_godunov_upwind_xpyp_ta, l_lmm_stepping, &
    stored_stats, rtp2, thlp2, rtpthlp, up2, vp2, sclrp2, sclrprtp, sclrpthlp, stored_err_info)

end subroutine f2py_advance_xp2_xpyp

subroutine f2py_update_xp2_mc(nzm, nzt, ngrdcol, dt, cloud_frac, rcm, rvm, thlm, &
    wm, exner, rrm_evap, rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, rtpthlp_mc)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_pdf_params
  use advance_xp2_xpyp_module, only: update_xp2_mc

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, rtpthlp_mc

  call update_xp2_mc(stored_grid, nzm, nzt, ngrdcol, dt, cloud_frac, rcm, rvm, thlm, &
    wm, exner, rrm_evap, stored_pdf_params, rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, &
    rtpthlp_mc)

end subroutine f2py_update_xp2_mc

! advance_wp2_wp3_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_advance_wp2_wp3(nzm, nzt, ngrdcol, dt, &
    sfc_elevation, fcor_y, sigma_sqd_w, wm_zm, wm_zt, &
    a3_coef, a3_coef_zt, wp3_on_wp2, &
    wpup2, wpvp2, wp2up2, wp2vp2, wp4, &
    wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp, &
    em, kh_zm, kh_zt, invrs_tau_c4_zm, invrs_tau_wp3_zt, invrs_tau_c1_zm, &
    skw_zm, skw_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    thv_ds_zm, thv_ds_zt, mixt_frac, cx_fnc_richardson, &
    lhs_splat_wp2, lhs_splat_wp3, &
    wprtp, wpthlp, rtp2, thlp2, clubb_params, &
    iipdf_type, penta_solve_method, fill_holes_type, &
    l_min_wp2_from_corr_wx, l_upwind_xm_ma, l_tke_aniso, l_standard_term_ta, &
    l_partial_upwind_wp3, l_damp_wp2_using_em, l_use_c11_richardson, &
    l_damp_wp3_skw_squared, l_lmm_stepping, &
    l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_k_dfsn, &
    l_use_wp3_lim_with_smth_heaviside, l_wp2_fill_holes_tke, l_ho_nontrad_coriolis, &
    up2, vp2, wp2, wp3, wp3_zm, wp2_zt)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_stats, stored_pdf_implicit_coefs_terms, &
    stored_nu_vert_res_dep, stored_err_info
  use advance_wp2_wp3_module, only: advance_wp2_wp3

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol), intent(in) :: sfc_elevation, fcor_y
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    sigma_sqd_w, wm_zm, a3_coef, wp3_on_wp2, wp2up2, wp2vp2, wp4, &
    wpthvp, upwp, vpwp, em, kh_zm, invrs_tau_c4_zm, invrs_tau_c1_zm, &
    skw_zm, rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, cx_fnc_richardson, &
    lhs_splat_wp2, wprtp, wpthlp, rtp2, thlp2
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    wm_zt, a3_coef_zt, wpup2, wpvp2, wp2thvp, wp2up, um, vm, kh_zt, &
    invrs_tau_wp3_zt, skw_zt, rho_ds_zt, invrs_rho_ds_zt, thv_ds_zt, &
    mixt_frac, lhs_splat_wp3
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  integer, intent(in) :: iipdf_type, penta_solve_method, fill_holes_type
  logical, intent(in) :: &
    l_min_wp2_from_corr_wx, l_upwind_xm_ma, l_tke_aniso, l_standard_term_ta, &
    l_partial_upwind_wp3, l_damp_wp2_using_em, l_use_c11_richardson, l_damp_wp3_skw_squared, &
    l_lmm_stepping, l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_k_dfsn, &
    l_use_wp3_lim_with_smth_heaviside, l_wp2_fill_holes_tke, l_ho_nontrad_coriolis
  logical, parameter :: l_implemented = .false.

  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: up2, vp2, wp2, wp3_zm
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: wp3, wp2_zt

  call advance_wp2_wp3(nzm, nzt, ngrdcol, stored_grid, dt, &
    sfc_elevation, fcor_y, sigma_sqd_w, wm_zm, wm_zt, &
    a3_coef, a3_coef_zt, wp3_on_wp2, &
    wpup2, wpvp2, wp2up2, wp2vp2, wp4, &
    wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp, &
    em, kh_zm, kh_zt, invrs_tau_c4_zm, invrs_tau_wp3_zt, invrs_tau_c1_zm, &
    skw_zm, skw_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    thv_ds_zm, thv_ds_zt, mixt_frac, cx_fnc_richardson, &
    lhs_splat_wp2, lhs_splat_wp3, &
    stored_pdf_implicit_coefs_terms, &
    wprtp, wpthlp, rtp2, thlp2, &
    clubb_params, stored_nu_vert_res_dep, &
    iipdf_type, penta_solve_method, fill_holes_type, &
    l_min_wp2_from_corr_wx, l_upwind_xm_ma, &
    l_tke_aniso, l_standard_term_ta, &
    l_partial_upwind_wp3, l_damp_wp2_using_em, &
    l_use_c11_richardson, l_damp_wp3_skw_squared, &
    l_lmm_stepping, l_use_tke_in_wp3_pr_turb_term, &
    l_use_tke_in_wp2_wp3_k_dfsn, &
    l_use_wp3_lim_with_smth_heaviside, &
    l_wp2_fill_holes_tke, l_ho_nontrad_coriolis, &
    l_implemented, stored_stats, up2, vp2, wp2, wp3, wp3_zm, wp2_zt, stored_err_info)

end subroutine f2py_advance_wp2_wp3

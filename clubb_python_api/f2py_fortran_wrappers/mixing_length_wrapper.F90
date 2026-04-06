! mixing_length_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module mixing_length

subroutine f2py_set_lscale_max(ngrdcol, l_implemented, host_dx, host_dy, lscale_max)

  use clubb_precision, only: core_rknd
  use mixing_length, only: set_Lscale_max

  implicit none

  integer, intent(in) :: ngrdcol
  logical, intent(in) :: l_implemented
  real(core_rknd), dimension(ngrdcol), intent(in) :: host_dx, host_dy
  real(core_rknd), dimension(ngrdcol), intent(out) :: lscale_max

  call set_Lscale_max(ngrdcol, l_implemented, host_dx, host_dy, lscale_max)

end subroutine f2py_set_lscale_max

subroutine f2py_diagnose_lscale_from_tau(nzm, nzt, ngrdcol, &
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, &
    ice_supersat_frac, em, sqrt_em_zt, &
    ufmin, tau_const, sfc_elevation, lscale_max, clubb_params, &
    l_e3sm_config, l_smooth_heaviside_tau_wpxp, &
    brunt_vaisala_freq_sqd_smth, ri_zm, &
    invrs_tau_zt, invrs_tau_zm, &
    invrs_tau_sfc, invrs_tau_no_n2_zm, invrs_tau_bkgnd, &
    invrs_tau_shear, invrs_tau_n2_iso, &
    invrs_tau_wp2_zm, invrs_tau_xp2_zm, &
    invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm, &
    tau_max_zm, tau_max_zt, tau_zm, tau_zt, &
    lscale, lscale_up, lscale_down)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_grid, stored_stats, stored_err_info
  use mixing_length, only: diagnose_lscale_from_tau

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol), intent(in) :: upwp_sfc, vpwp_sfc
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: ice_supersat_frac, sqrt_em_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    ddzt_umvm_sqd, em, brunt_vaisala_freq_sqd_smth, ri_zm
  real(core_rknd), intent(in) :: ufmin, tau_const
  real(core_rknd), dimension(ngrdcol), intent(in) :: sfc_elevation, lscale_max
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  logical, intent(in) :: l_e3sm_config, l_smooth_heaviside_tau_wpxp

  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: &
    invrs_tau_zt, invrs_tau_wp3_zt, tau_max_zt, tau_zt, lscale, lscale_up, lscale_down
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: &
    invrs_tau_zm, invrs_tau_sfc, invrs_tau_no_n2_zm, invrs_tau_bkgnd, &
    invrs_tau_shear, invrs_tau_n2_iso, invrs_tau_wp2_zm, invrs_tau_xp2_zm, &
    invrs_tau_wp3_zm, invrs_tau_wpxp_zm, tau_max_zm, tau_zm

  call diagnose_lscale_from_tau(nzm, nzt, ngrdcol, stored_grid, &
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, &
    ice_supersat_frac, em, sqrt_em_zt, &
    ufmin, tau_const, sfc_elevation, lscale_max, clubb_params, &
    stored_stats, l_e3sm_config, l_smooth_heaviside_tau_wpxp, &
    brunt_vaisala_freq_sqd_smth, ri_zm, stored_err_info, &
    invrs_tau_zt, invrs_tau_zm, &
    invrs_tau_sfc, invrs_tau_no_n2_zm, invrs_tau_bkgnd, &
    invrs_tau_shear, invrs_tau_n2_iso, &
    invrs_tau_wp2_zm, invrs_tau_xp2_zm, &
    invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm, &
    tau_max_zm, tau_max_zt, tau_zm, tau_zt, &
    lscale, lscale_up, lscale_down)

end subroutine f2py_diagnose_lscale_from_tau

subroutine f2py_calc_lscale_directly(ngrdcol, nzm, nzt, &
    l_implemented, p_in_pa, exner, rtm, thlm, thvm, newmu, &
    rtp2_zt, thlp2_zt, rtpthlp_zt, em, thv_ds_zt, lscale_max, lmin, &
    clubb_params, saturation_formula, l_lscale_plume_centered, &
    lscale, lscale_up, lscale_down)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_pdf_params, stored_stats, stored_err_info
  use mixing_length, only: calc_lscale_directly

  implicit none

  integer, intent(in) :: ngrdcol, nzm, nzt
  integer, intent(in) :: saturation_formula
  logical, intent(in) :: l_implemented, l_lscale_plume_centered
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    p_in_pa, exner, rtm, thlm, thvm, rtp2_zt, thlp2_zt, rtpthlp_zt, thv_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: em
  real(core_rknd), dimension(ngrdcol), intent(in) :: newmu, lscale_max
  real(core_rknd), intent(in) :: lmin
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: lscale, lscale_up, lscale_down

  call calc_lscale_directly(ngrdcol, nzm, nzt, stored_grid, &
    l_implemented, p_in_pa, exner, rtm, thlm, thvm, newmu, &
    rtp2_zt, thlp2_zt, rtpthlp_zt, stored_pdf_params, em, thv_ds_zt, lscale_max, lmin, &
    clubb_params, saturation_formula, l_lscale_plume_centered, &
    stored_stats, stored_err_info, lscale, lscale_up, lscale_down)

end subroutine f2py_calc_lscale_directly


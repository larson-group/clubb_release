! mixing_length_wrapper.F90
!
! Expose the Fortran mixing-length routines to the compiled Python API.
!
! Usage:
!   f2py_calc_lscale(): complete length/time diagnosis for Lscale_method 1--3.
!   f2py_set_lscale_max(): standalone or host-spacing length cap.
!   f2py_diagnose_lscale_from_tau(): method 2 and its inverse-time contributions.
!   f2py_calc_lscale_directly(): legacy/coherent parcel methods 1 or 3.
!
! Understanding:
!   Python's public wrappers populate derived_type_storage before these calls.
!   Numeric fields retain (ngrdcol,nzt) thermodynamic or (ngrdcol,nzm) momentum
!   layout. Grid, PDF, statistics, and error objects use shared Fortran storage;
!   callers must initialize them and retrieve the updated error object afterward.
!   Lengths are in m, time scales s, and inverse times 1/s. No files are written
!   by these wrappers. src/CLUBB_core/mixing_length.F90 defines the calculations.

subroutine f2py_set_lscale_max(ngrdcol, l_implemented, host_dx, host_dy, lscale_max)

  ! Description:
  !   Return per-column length caps [m] from host_dx/host_dy [m], or 100000 m
  !   when l_implemented is false. No stored derived-type state is needed.

  use clubb_precision, only: core_rknd
  use mixing_length, only: set_Lscale_max

  implicit none

  integer, intent(in) :: ngrdcol
  logical, intent(in) :: l_implemented
  real(core_rknd), dimension(ngrdcol), intent(in) :: host_dx, host_dy
  real(core_rknd), dimension(ngrdcol), intent(out) :: lscale_max

  call set_Lscale_max(ngrdcol, l_implemented, host_dx, host_dy, lscale_max)

end subroutine f2py_set_lscale_max

subroutine f2py_calc_lscale(nzm, nzt, ngrdcol, l_implemented, host_dx, host_dy, &
    p_in_pa, exner, rtm, thlm, thvm, thlp2, rtp2, rtpthlp, em, &
    thv_ds_zt, lmin, &
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac, &
    ufmin, tau_const, sfc_elevation, clubb_params, &
    saturation_formula, l_lscale_plume_centered, lscale_method, &
    l_e3sm_config, &
    l_smooth_heaviside_tau_wpxp, l_modify_limiters_for_cnvg_test, l_use_invrs_tau_n2_iso, &
    brunt_vaisala_freq_sqd_smth, &
    invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm, invrs_tau_wp3_zt, &
    invrs_tau_c1_zm, invrs_tau_c4_zm, invrs_tau_c6_zm, invrs_tau_c14_zm, &
    tau_max_zm, tau_max_zt, tau_zm, lscale, lscale_zm, lscale_up, lscale_down)

  ! Description:
  !   Read thermodynamic fields, momentum-level moments/TKE, per-column
  !   parameters, and method selector 1--3. Return length/time profiles and
  !   update stored_stats/stored_err_info using stored_grid/stored_pdf_params.
  !   Output suffixes identify grids; lengths are m, times s, inverse times 1/s.

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_pdf_params, stored_stats, stored_err_info
  use mixing_length, only: calc_Lscale

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  integer, intent(in) :: saturation_formula, lscale_method
  logical, intent(in) :: &
    l_implemented, l_lscale_plume_centered, l_e3sm_config, &
    l_smooth_heaviside_tau_wpxp, l_modify_limiters_for_cnvg_test, l_use_invrs_tau_n2_iso
  real(core_rknd), intent(in) :: lmin, ufmin, tau_const
  real(core_rknd), dimension(ngrdcol), intent(in) :: &
    host_dx, host_dy, upwp_sfc, vpwp_sfc, sfc_elevation
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    p_in_pa, exner, rtm, thlm, thvm, thv_ds_zt, ice_supersat_frac
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    thlp2, rtp2, rtpthlp, em, ddzt_umvm_sqd, brunt_vaisala_freq_sqd_smth
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params

  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: &
    invrs_tau_zt, invrs_tau_wp3_zt, tau_max_zt, lscale, lscale_up, lscale_down
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: &
    invrs_tau_zm, invrs_tau_xp2_zm, invrs_tau_c1_zm, invrs_tau_c4_zm, &
    invrs_tau_c6_zm, invrs_tau_c14_zm, tau_max_zm, tau_zm, lscale_zm

  call calc_Lscale(nzm, nzt, ngrdcol, stored_grid, l_implemented, host_dx, host_dy, &
    p_in_pa, exner, rtm, thlm, thvm, thlp2, rtp2, rtpthlp, stored_pdf_params, &
    em, &
    thv_ds_zt, lmin, upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac, &
    ufmin, tau_const, sfc_elevation, clubb_params, saturation_formula, &
    l_lscale_plume_centered, lscale_method, &
    l_e3sm_config, &
    l_smooth_heaviside_tau_wpxp, l_modify_limiters_for_cnvg_test, &
    l_use_invrs_tau_n2_iso, brunt_vaisala_freq_sqd_smth, stored_stats, stored_err_info, &
    invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm, invrs_tau_wp3_zt, &
    invrs_tau_c1_zm, invrs_tau_c4_zm, invrs_tau_c6_zm, invrs_tau_c14_zm, &
    tau_max_zm, tau_max_zt, tau_zm, lscale, lscale_zm, lscale_up, lscale_down)

end subroutine f2py_calc_lscale

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

  ! Description:
  !   Diagnose method-2 inverse-time contributions [1/s], bounded times [s],
  !   and lengths [m] from surface stress, TKE, shear, and stability. Use
  !   stored_grid and update stored_stats/stored_err_info. Directional lengths
  !   are zero placeholders for this method.

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_grid, stored_stats, stored_err_info
  use mixing_length, only: diagnose_lscale_from_tau

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol), intent(in) :: upwp_sfc, vpwp_sfc
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: ice_supersat_frac, sqrt_em_zt
  ! TODO: Remove unused sqrt_em_zt with a coordinated binding/public API update;
  ! the core computes its own thermodynamic-level TKE interpolation.
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
    ice_supersat_frac, em, &
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
    l_implemented, p_in_pa, exner, rtm, thlm, thvm, &
    rtp2_zt, thlp2_zt, rtpthlp_zt, em, thv_ds_zt, lscale_max, lmin, &
    clubb_params, saturation_formula, l_lscale_plume_centered, &
    lscale_method, &
    lscale, lscale_up, lscale_down)

  ! Description:
  !   Select parcel method 1 or 3 from environmental/PDF properties, TKE,
  !   and per-column parameters. Return three thermodynamic-level lengths [m].
  !   Read stored_grid/stored_pdf_params and update stored_stats/stored_err_info.

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_pdf_params, stored_stats, stored_err_info
  use mixing_length, only: calc_lscale_directly

  implicit none

  integer, intent(in) :: ngrdcol, nzm, nzt
  integer, intent(in) :: saturation_formula, lscale_method
  logical, intent(in) :: &
    l_implemented, l_lscale_plume_centered
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    p_in_pa, exner, rtm, thlm, thvm, rtp2_zt, thlp2_zt, rtpthlp_zt, thv_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: em
  real(core_rknd), dimension(ngrdcol), intent(in) :: lscale_max
  real(core_rknd), intent(in) :: lmin
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: lscale, lscale_up, lscale_down

  call calc_lscale_directly(ngrdcol, nzm, nzt, stored_grid, &
    l_implemented, p_in_pa, exner, rtm, thlm, thvm, &
    rtp2_zt, thlp2_zt, rtpthlp_zt, stored_pdf_params, em, &
    thv_ds_zt, lscale_max, lmin, &
    clubb_params, saturation_formula, l_lscale_plume_centered, &
    lscale_method, &
    stored_stats, stored_err_info, lscale, lscale_up, lscale_down)

end subroutine f2py_calc_lscale_directly

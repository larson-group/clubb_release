! mono_flux_limiter_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module mono_flux_limiter

subroutine f2py_calc_turb_adv_range(nzm, nzt, ngrdcol, dt, &
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, &
    low_lev_effect, high_lev_effect)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_stats
  use mono_flux_limiter, only: calc_turb_adv_range

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm
  integer, dimension(ngrdcol, nzt), intent(out) :: low_lev_effect, high_lev_effect

  call calc_turb_adv_range(nzm, nzt, ngrdcol, stored_grid, dt, &
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm, &
    stored_stats, low_lev_effect, high_lev_effect)

end subroutine f2py_calc_turb_adv_range

subroutine f2py_monotonic_turbulent_flux_limit(nzm, nzt, ngrdcol, solve_type, dt, &
    xm_old, xp2, wm_zt, xm_forcing, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    xp2_threshold, xm_tol, l_implemented, low_lev_effect, high_lev_effect, &
    tridiag_solve_method, l_upwind_xm_ma, l_mono_flux_lim_spikefix, xm, wpxp)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_stats, stored_err_info
  use mono_flux_limiter, only: monotonic_turbulent_flux_limit

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, solve_type, tridiag_solve_method
  real(core_rknd), intent(in) :: dt, xp2_threshold, xm_tol
  logical, intent(in) :: l_implemented, l_upwind_xm_ma, l_mono_flux_lim_spikefix
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: xm_old, wm_zt, xm_forcing, rho_ds_zt, invrs_rho_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: xp2, rho_ds_zm, invrs_rho_ds_zm
  integer, dimension(ngrdcol, nzt), intent(in) :: low_lev_effect, high_lev_effect
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: xm
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: wpxp

  call monotonic_turbulent_flux_limit( &
    nzm, nzt, ngrdcol, stored_grid, solve_type, dt, &
    xm_old, xp2, wm_zt, xm_forcing, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
    xp2_threshold, xm_tol, l_implemented, &
    low_lev_effect, high_lev_effect, tridiag_solve_method, &
    l_upwind_xm_ma, l_mono_flux_lim_spikefix, stored_stats, xm, wpxp, stored_err_info)

end subroutine f2py_monotonic_turbulent_flux_limit

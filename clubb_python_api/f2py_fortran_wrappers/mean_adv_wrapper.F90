! mean_adv_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module mean_adv

subroutine f2py_term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt, weights_zt2zm, &
    invrs_dzt, invrs_dzm, l_upwind_xm_ma, grid_dir, lhs_ma)

  use clubb_precision, only: core_rknd
  use mean_adv, only: term_ma_zt_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: wm_zt
  real(core_rknd), dimension(ngrdcol, nzm, 2), intent(in) :: weights_zt2zm
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: invrs_dzt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: invrs_dzm
  logical, intent(in) :: l_upwind_xm_ma
  real(core_rknd), intent(in) :: grid_dir
  real(core_rknd), dimension(3, ngrdcol, nzt), intent(out) :: lhs_ma

  call term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt, weights_zt2zm, invrs_dzt, &
    invrs_dzm, l_upwind_xm_ma, grid_dir, lhs_ma)

end subroutine f2py_term_ma_zt_lhs

subroutine f2py_term_ma_zm_lhs(nzm, nzt, ngrdcol, wm_zm, invrs_dzm, &
    weights_zm2zt, lhs_ma)

  use clubb_precision, only: core_rknd
  use mean_adv, only: term_ma_zm_lhs

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: wm_zm
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: invrs_dzm
  real(core_rknd), dimension(ngrdcol, nzt, 2), intent(in) :: weights_zm2zt
  real(core_rknd), dimension(3, ngrdcol, nzm), intent(out) :: lhs_ma

  call term_ma_zm_lhs(nzm, nzt, ngrdcol, wm_zm, invrs_dzm, weights_zm2zt, lhs_ma)

end subroutine f2py_term_ma_zm_lhs

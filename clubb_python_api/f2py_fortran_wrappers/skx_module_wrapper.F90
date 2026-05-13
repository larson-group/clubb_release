! skx_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module skx_module

subroutine f2py_skx_func(nz, ngrdcol, xp2, xp3, x_tol, clubb_params, skx)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use skx_module, only: skx_func

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xp2, xp3
  real(core_rknd), intent(in) :: x_tol
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: skx

  call skx_func(nz, ngrdcol, xp2, xp3, x_tol, clubb_params, skx)

end subroutine f2py_skx_func

subroutine f2py_lg_2005_ansatz(nz, ngrdcol, skw, wpxp, wp2, xp2, beta, &
    sigma_sqd_w, x_tol, skx)

  use clubb_precision, only: core_rknd
  use skx_module, only: lg_2005_ansatz

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    skw, wpxp, wp2, xp2, sigma_sqd_w
  real(core_rknd), dimension(ngrdcol), intent(in) :: beta
  real(core_rknd), intent(in) :: x_tol
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: skx

  call lg_2005_ansatz(nz, ngrdcol, skw, wpxp, wp2, xp2, beta, sigma_sqd_w, x_tol, skx)

end subroutine f2py_lg_2005_ansatz

subroutine f2py_xp3_lg_2005_ansatz(nzt, ngrdcol, skw_zt, wpxp_zt, wp2_zt, &
    xp2_zt, sigma_sqd_w_zt, clubb_params, x_tol, xp3)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use skx_module, only: xp3_lg_2005_ansatz

  implicit none

  integer, intent(in) :: nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    skw_zt, wpxp_zt, wp2_zt, xp2_zt, sigma_sqd_w_zt
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), intent(in) :: x_tol
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: xp3

  call xp3_lg_2005_ansatz(nzt, ngrdcol, skw_zt, wpxp_zt, wp2_zt, xp2_zt, &
    sigma_sqd_w_zt, clubb_params, x_tol, xp3)

end subroutine f2py_xp3_lg_2005_ansatz

subroutine f2py_compute_gamma_skw(nzm, nzt, ngrdcol, l_gamma_skw, &
    skw_zm, clubb_params, skw_zt, gamma_skw_fnc, gamma_skw_fnc_zt)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use skx_module, only: compute_gamma_Skw

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  logical, intent(in) :: l_gamma_skw
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: skw_zm
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: skw_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: gamma_skw_fnc
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: gamma_skw_fnc_zt

  call compute_gamma_Skw(nzm, nzt, ngrdcol, l_gamma_skw, &
    skw_zm, clubb_params, gamma_skw_fnc, skw_zt, gamma_skw_fnc_zt)

end subroutine f2py_compute_gamma_skw

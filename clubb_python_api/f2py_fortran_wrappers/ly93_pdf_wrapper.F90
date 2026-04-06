! ly93_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module ly93_pdf

subroutine f2py_calc_params_ly93(nz, ngrdcol, xm, xp2, skx, mixt_frac, &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use LY93_pdf, only: calc_params_LY93

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm, xp2, skx, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd

  integer :: i

  do i = 1, ngrdcol
    call calc_params_LY93(nz, xm(i,:), xp2(i,:), skx(i,:), mixt_frac(i,:), &
      mu_x_1(i,:), mu_x_2(i,:), sigma_x_1_sqd(i,:), sigma_x_2_sqd(i,:))
  end do

end subroutine f2py_calc_params_ly93

subroutine f2py_ly93_driver(nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, &
    skw, skrt, skthl, mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, &
    sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, &
    sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac)

  use clubb_precision, only: core_rknd
  use LY93_pdf, only: LY93_driver

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, &
    sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, &
    sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac

  integer :: i

  do i = 1, ngrdcol
    call LY93_driver(nz, wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), &
      thlp2(i,:), skw(i,:), skrt(i,:), skthl(i,:), &
      mu_w_1(i,:), mu_w_2(i,:), mu_rt_1(i,:), mu_rt_2(i,:), &
      mu_thl_1(i,:), mu_thl_2(i,:), &
      sigma_w_1_sqd(i,:), sigma_w_2_sqd(i,:), &
      sigma_rt_1_sqd(i,:), sigma_rt_2_sqd(i,:), &
      sigma_thl_1_sqd(i,:), sigma_thl_2_sqd(i,:), mixt_frac(i,:))
  end do

end subroutine f2py_ly93_driver


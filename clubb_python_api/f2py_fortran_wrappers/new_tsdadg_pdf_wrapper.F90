! new_tsdadg_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_tsdadg_pdf

subroutine f2py_calc_l_x_skx_fnc(skx, sgn_wpxp, small_l_x_1, small_l_x_2, &
    big_l_x_1, big_l_x_2)

  use clubb_precision, only: core_rknd
  use new_tsdadg_pdf, only: calc_L_x_Skx_fnc

  implicit none

  real(core_rknd), intent(in) :: skx, sgn_wpxp, small_l_x_1, small_l_x_2
  real(core_rknd), intent(out) :: big_l_x_1, big_l_x_2

  call calc_L_x_Skx_fnc(skx, sgn_wpxp, small_l_x_1, small_l_x_2, big_l_x_1, big_l_x_2)

end subroutine f2py_calc_l_x_skx_fnc

subroutine f2py_calc_setter_parameters_tsdadg(xm, xp2, skx, sgn_wpxp, &
    big_l_x_1, big_l_x_2, mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    mixt_frac, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use new_tsdadg_pdf, only: calc_setter_parameters

  implicit none

  real(core_rknd), intent(in) :: xm, xp2, skx, sgn_wpxp, big_l_x_1, big_l_x_2
  real(core_rknd), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, mixt_frac, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd

  call calc_setter_parameters(xm, xp2, skx, sgn_wpxp, big_l_x_1, big_l_x_2, &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, mixt_frac, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

end subroutine f2py_calc_setter_parameters_tsdadg

subroutine f2py_tsdadg_pdf_driver(nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, &
    skw, skrt, skthl, wprtp, wpthlp, mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, &
    mu_thl_2, sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, &
    sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac)

  use clubb_precision, only: core_rknd
  use new_tsdadg_pdf, only: tsdadg_pdf_driver

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, &
    sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, &
    sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac

  integer :: i

  do i = 1, ngrdcol
    call tsdadg_pdf_driver(nz, wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:), &
      skw(i,:), skrt(i,:), skthl(i,:), wprtp(i,:), wpthlp(i,:), &
      mu_w_1(i,:), mu_w_2(i,:), mu_rt_1(i,:), mu_rt_2(i,:), mu_thl_1(i,:), mu_thl_2(i,:), &
      sigma_w_1_sqd(i,:), sigma_w_2_sqd(i,:), sigma_rt_1_sqd(i,:), sigma_rt_2_sqd(i,:), &
      sigma_thl_1_sqd(i,:), sigma_thl_2_sqd(i,:), mixt_frac(i,:))
  end do

end subroutine f2py_tsdadg_pdf_driver


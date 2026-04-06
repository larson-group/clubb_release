! new_pdf_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_pdf

subroutine f2py_calc_setter_var_params(nz, ngrdcol, xm, xp2, skx, sgn_wpxp, &
    f_x, zeta_x, mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, mixt_frac, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_setter_var_params

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: xm, xp2, skx, sgn_wpxp
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: f_x, zeta_x
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, mixt_frac, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd

  integer :: i

  do i = 1, ngrdcol
    call calc_setter_var_params(nz, xm(i,:), xp2(i,:), skx(i,:), sgn_wpxp(i,:), &
      f_x(i,:), zeta_x(i,:), mu_x_1(i,:), mu_x_2(i,:), sigma_x_1(i,:), sigma_x_2(i,:), &
      mixt_frac(i,:), coef_sigma_x_1_sqd(i,:), coef_sigma_x_2_sqd(i,:))
  end do

end subroutine f2py_calc_setter_var_params

subroutine f2py_calc_responder_params(nz, ngrdcol, xm, xp2, skx, sgn_wpxp, &
    f_x, mixt_frac, mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_responder_params

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    xm, xp2, skx, sgn_wpxp, f_x, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd

  integer :: i

  do i = 1, ngrdcol
    call calc_responder_params(nz, xm(i,:), xp2(i,:), skx(i,:), sgn_wpxp(i,:), &
      f_x(i,:), mixt_frac(i,:), mu_x_1(i,:), mu_x_2(i,:), sigma_x_1_sqd(i,:), &
      sigma_x_2_sqd(i,:), coef_sigma_x_1_sqd(i,:), coef_sigma_x_2_sqd(i,:))
  end do

end subroutine f2py_calc_responder_params

subroutine f2py_calc_limits_f_x_responder(nz, ngrdcol, mixt_frac, skx, sgn_wpxp, &
    max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp, min_f_x, max_f_x)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_limits_f_x_responder

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    mixt_frac, skx, sgn_wpxp, max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: min_f_x, max_f_x

  integer :: i

  do i = 1, ngrdcol
    call calc_limits_f_x_responder(nz, mixt_frac(i,:), skx(i,:), sgn_wpxp(i,:), &
      max_skx2_pos_skx_sgn_wpxp(i,:), max_skx2_neg_skx_sgn_wpxp(i,:), &
      min_f_x(i,:), max_f_x(i,:))
  end do

end subroutine f2py_calc_limits_f_x_responder

subroutine f2py_calc_coef_wp4_implicit(nz, ngrdcol, mixt_frac, f_w, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_wp4_implicit)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coef_wp4_implicit

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wp4_implicit

  integer :: i

  do i = 1, ngrdcol
    coef_wp4_implicit(i,:) = calc_coef_wp4_implicit(nz, mixt_frac(i,:), f_w(i,:), &
      coef_sigma_w_1_sqd(i,:), coef_sigma_w_2_sqd(i,:))
  end do

end subroutine f2py_calc_coef_wp4_implicit

subroutine f2py_calc_coef_wpxp2_implicit(nz, ngrdcol, wp2, xp2, wpxp, sgn_wpxp, &
    mixt_frac, f_w, f_x, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, coef_wpxp2_implicit)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coef_wpxp2_implicit

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wp2, xp2, wpxp, sgn_wpxp, mixt_frac, f_w, f_x, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wpxp2_implicit

  integer :: i

  do i = 1, ngrdcol
    coef_wpxp2_implicit(i,:) = calc_coef_wpxp2_implicit(nz, wp2(i,:), xp2(i,:), &
      wpxp(i,:), sgn_wpxp(i,:), mixt_frac(i,:), f_w(i,:), f_x(i,:), &
      coef_sigma_w_1_sqd(i,:), coef_sigma_w_2_sqd(i,:), &
      coef_sigma_x_1_sqd(i,:), coef_sigma_x_2_sqd(i,:))
  end do

end subroutine f2py_calc_coef_wpxp2_implicit

subroutine f2py_calc_coefs_wp2xp_semiimpl(nz, ngrdcol, wp2, xp2, sgn_wpxp, &
    mixt_frac, f_w, f_x, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, coef_wp2xp_implicit, &
    term_wp2xp_explicit)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coefs_wp2xp_semiimpl

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wp2, xp2, sgn_wpxp, mixt_frac, f_w, f_x, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    coef_wp2xp_implicit, term_wp2xp_explicit

  integer :: i

  do i = 1, ngrdcol
    call calc_coefs_wp2xp_semiimpl(nz, wp2(i,:), xp2(i,:), sgn_wpxp(i,:), &
      mixt_frac(i,:), f_w(i,:), f_x(i,:), coef_sigma_w_1_sqd(i,:), &
      coef_sigma_w_2_sqd(i,:), coef_sigma_x_1_sqd(i,:), coef_sigma_x_2_sqd(i,:), &
      coef_wp2xp_implicit(i,:), term_wp2xp_explicit(i,:))
  end do

end subroutine f2py_calc_coefs_wp2xp_semiimpl

subroutine f2py_calc_coefs_wpxpyp_semiimpl(nz, ngrdcol, wp2, xp2, yp2, wpxp, &
    wpyp, sgn_wpxp, sgn_wpyp, mixt_frac, f_w, f_x, f_y, coef_sigma_w_1_sqd, &
    coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, &
    coef_sigma_y_1_sqd, coef_sigma_y_2_sqd, coef_wpxpyp_implicit, &
    term_wpxpyp_explicit)

  use clubb_precision, only: core_rknd
  use new_pdf, only: calc_coefs_wpxpyp_semiimpl

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wp2, xp2, yp2, wpxp, wpyp, sgn_wpxp, sgn_wpyp, mixt_frac, f_w, f_x, f_y, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, &
    coef_sigma_y_1_sqd, coef_sigma_y_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    coef_wpxpyp_implicit, term_wpxpyp_explicit

  integer :: i

  do i = 1, ngrdcol
    call calc_coefs_wpxpyp_semiimpl(nz, wp2(i,:), xp2(i,:), yp2(i,:), wpxp(i,:), &
      wpyp(i,:), sgn_wpxp(i,:), sgn_wpyp(i,:), mixt_frac(i,:), f_w(i,:), f_x(i,:), &
      f_y(i,:), coef_sigma_w_1_sqd(i,:), coef_sigma_w_2_sqd(i,:), &
      coef_sigma_x_1_sqd(i,:), coef_sigma_x_2_sqd(i,:), coef_sigma_y_1_sqd(i,:), &
      coef_sigma_y_2_sqd(i,:), coef_wpxpyp_implicit(i,:), term_wpxpyp_explicit(i,:))
  end do

end subroutine f2py_calc_coefs_wpxpyp_semiimpl

! new_hybrid_pdf_wrapper.F90 — wrappers for helper routines from module new_hybrid_pdf

subroutine f2py_calculate_w_params(nz, ngrdcol, wm, wp2, skw, f_w, zeta_w, &
    mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac, coef_sigma_w_1_sqd, &
    coef_sigma_w_2_sqd)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calculate_w_params

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: wm, wp2, skw, f_w, zeta_w
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac, coef_sigma_w_1_sqd, &
    coef_sigma_w_2_sqd

  integer :: i, k

  do i = 1, ngrdcol
    do k = 1, nz
      call calculate_w_params(wm(i,k), wp2(i,k), skw(i,k), f_w(i,k), zeta_w(i,k), &
        mu_w_1(i,k), mu_w_2(i,k), sigma_w_1(i,k), sigma_w_2(i,k), mixt_frac(i,k), &
        coef_sigma_w_1_sqd(i,k), coef_sigma_w_2_sqd(i,k))
    end do
  end do

end subroutine f2py_calculate_w_params

subroutine f2py_calculate_responder_params(nz, ngrdcol, xm, xp2, skx, wpxp, wp2, &
    f_w, mixt_frac, mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calculate_responder_params

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    xm, xp2, skx, wpxp, wp2, f_w, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, coef_sigma_x_1_sqd, &
    coef_sigma_x_2_sqd

  integer :: i, k

  do i = 1, ngrdcol
    do k = 1, nz
      call calculate_responder_params(xm(i,k), xp2(i,k), skx(i,k), wpxp(i,k), &
        wp2(i,k), f_w(i,k), mixt_frac(i,k), mu_x_1(i,k), mu_x_2(i,k), &
        sigma_x_1_sqd(i,k), sigma_x_2_sqd(i,k), coef_sigma_x_1_sqd(i,k), &
        coef_sigma_x_2_sqd(i,k))
    end do
  end do

end subroutine f2py_calculate_responder_params

subroutine f2py_calculate_coef_wp4_implicit(nz, ngrdcol, mixt_frac, f_w, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_wp4_implicit)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calculate_coef_wp4_implicit

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wp4_implicit

  integer :: i, k

  do i = 1, ngrdcol
    do k = 1, nz
      coef_wp4_implicit(i,k) = calculate_coef_wp4_implicit(mixt_frac(i,k), f_w(i,k), &
        coef_sigma_w_1_sqd(i,k), coef_sigma_w_2_sqd(i,k))
    end do
  end do

end subroutine f2py_calculate_coef_wp4_implicit

subroutine f2py_calc_coef_wp2xp_implicit(nz, ngrdcol, wp2, mixt_frac, f_w, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_wp2xp_implicit)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calc_coef_wp2xp_implicit

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wp2, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: coef_wp2xp_implicit

  integer :: i, k

  do i = 1, ngrdcol
    do k = 1, nz
      coef_wp2xp_implicit(i,k) = calc_coef_wp2xp_implicit(wp2(i,k), mixt_frac(i,k), &
        f_w(i,k), coef_sigma_w_1_sqd(i,k), coef_sigma_w_2_sqd(i,k))
    end do
  end do

end subroutine f2py_calc_coef_wp2xp_implicit

subroutine f2py_calc_coefs_wpxp2_semiimpl(nz, ngrdcol, wp2, wpxp, mixt_frac, f_w, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, coef_wpxp2_implicit, &
    term_wpxp2_explicit)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calc_coefs_wpxp2_semiimpl

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wp2, wpxp, mixt_frac, f_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    coef_wpxp2_implicit, term_wpxp2_explicit

  integer :: i, k

  do i = 1, ngrdcol
    do k = 1, nz
      call calc_coefs_wpxp2_semiimpl(wp2(i,k), wpxp(i,k), mixt_frac(i,k), f_w(i,k), &
        coef_sigma_x_1_sqd(i,k), coef_sigma_x_2_sqd(i,k), coef_wpxp2_implicit(i,k), &
        term_wpxp2_explicit(i,k))
    end do
  end do

end subroutine f2py_calc_coefs_wpxp2_semiimpl

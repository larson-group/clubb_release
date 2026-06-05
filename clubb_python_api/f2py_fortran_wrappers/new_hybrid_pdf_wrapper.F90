! new_hybrid_pdf_wrapper.F90 — wrappers for helper routines from module new_hybrid_pdf

subroutine f2py_calculate_w_params(wm, wp2, skw, f_w, zeta_w, &
    mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac, coef_sigma_w_1_sqd, &
    coef_sigma_w_2_sqd)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calculate_w_params

  implicit none

  real(core_rknd), intent(in) :: wm, wp2, skw, f_w, zeta_w
  real(core_rknd), intent(out) :: &
    mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac, coef_sigma_w_1_sqd, &
    coef_sigma_w_2_sqd

  call calculate_w_params(wm, wp2, skw, f_w, zeta_w, mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, &
    mixt_frac, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)

end subroutine f2py_calculate_w_params

subroutine f2py_calculate_responder_params(xm, xp2, skx, wpxp, wp2, &
    f_w, mixt_frac, mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calculate_responder_params

  implicit none

  real(core_rknd), intent(in) :: &
    xm, xp2, skx, wpxp, wp2, f_w, mixt_frac
  real(core_rknd), intent(out) :: &
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, coef_sigma_x_1_sqd, &
    coef_sigma_x_2_sqd

  call calculate_responder_params(xm, xp2, skx, wpxp, wp2, f_w, mixt_frac, mu_x_1, mu_x_2, &
    sigma_x_1_sqd, sigma_x_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)

end subroutine f2py_calculate_responder_params

subroutine f2py_calculate_coef_wp4_implicit(mixt_frac, f_w, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_wp4_implicit)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calculate_coef_wp4_implicit

  implicit none

  real(core_rknd), intent(in) :: &
    mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd
  real(core_rknd), intent(out) :: coef_wp4_implicit

  coef_wp4_implicit = calculate_coef_wp4_implicit(mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)

end subroutine f2py_calculate_coef_wp4_implicit

subroutine f2py_calc_coef_wp2xp_implicit(wp2, mixt_frac, f_w, &
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_wp2xp_implicit)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calc_coef_wp2xp_implicit

  implicit none

  real(core_rknd), intent(in) :: &
    wp2, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd
  real(core_rknd), intent(out) :: coef_wp2xp_implicit

  coef_wp2xp_implicit = calc_coef_wp2xp_implicit(wp2, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)

end subroutine f2py_calc_coef_wp2xp_implicit

subroutine f2py_calc_coefs_wpxp2_semiimpl(wp2, wpxp, mixt_frac, f_w, &
    coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, coef_wpxp2_implicit, &
    term_wpxp2_explicit)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf, only: calc_coefs_wpxp2_semiimpl

  implicit none

  real(core_rknd), intent(in) :: &
    wp2, wpxp, mixt_frac, f_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd
  real(core_rknd), intent(out) :: &
    coef_wpxp2_implicit, term_wpxp2_explicit

  call calc_coefs_wpxp2_semiimpl(wp2, wpxp, mixt_frac, f_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd, &
    coef_wpxp2_implicit, term_wpxp2_explicit)

end subroutine f2py_calc_coefs_wpxp2_semiimpl

! pdf_utilities_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module pdf_utilities

subroutine f2py_smooth_corr_quotient(ngrdcol, nz, numerator, denominator, &
    denom_thresh, quotient)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: smooth_corr_quotient

  implicit none

  integer, intent(in) :: ngrdcol, nz
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: numerator, denominator
  real(core_rknd), intent(in) :: denom_thresh
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: quotient

  call smooth_corr_quotient(ngrdcol, nz, numerator, denominator, denom_thresh, quotient)

end subroutine f2py_smooth_corr_quotient

subroutine f2py_calc_comp_corrs_binormal(nz, ngrdcol, xpyp, xm, ym, &
    mu_x_1, mu_x_2, mu_y_1, mu_y_2, sigma_x_1_sqd, sigma_x_2_sqd, &
    sigma_y_1_sqd, sigma_y_2_sqd, mixt_frac, corr_x_y_1, corr_x_y_2)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: calc_comp_corrs_binormal

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    xpyp, xm, ym, mu_x_1, mu_x_2, mu_y_1, mu_y_2, &
    sigma_x_1_sqd, sigma_x_2_sqd, sigma_y_1_sqd, sigma_y_2_sqd, mixt_frac
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: corr_x_y_1, corr_x_y_2

  call calc_comp_corrs_binormal(nz, ngrdcol, xpyp, xm, ym, mu_x_1, mu_x_2, mu_y_1, mu_y_2, &
    sigma_x_1_sqd, sigma_x_2_sqd, sigma_y_1_sqd, sigma_y_2_sqd, mixt_frac, corr_x_y_1, corr_x_y_2)

end subroutine f2py_calc_comp_corrs_binormal

subroutine f2py_mean_l2n(mu_x, sigma2_on_mu2, mu_x_n)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: mean_L2N

  implicit none

  real(core_rknd), intent(in) :: mu_x, sigma2_on_mu2
  real(core_rknd), intent(out) :: mu_x_n

  mu_x_n = mean_L2N(mu_x, sigma2_on_mu2)

end subroutine f2py_mean_l2n

subroutine f2py_stdev_l2n(sigma2_on_mu2, sigma_x_n)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: stdev_L2N

  implicit none

  real(core_rknd), intent(in) :: sigma2_on_mu2
  real(core_rknd), intent(out) :: sigma_x_n

  sigma_x_n = stdev_L2N(sigma2_on_mu2)

end subroutine f2py_stdev_l2n

subroutine f2py_corr_nn2nl(nz, ngrdcol, corr_x_y_n, sigma_y_n, y_sigma2_on_mu2, corr_x_y)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: corr_NN2NL

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    corr_x_y_n, sigma_y_n, y_sigma2_on_mu2
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: corr_x_y

  call corr_NN2NL(nz, ngrdcol, corr_x_y_n, sigma_y_n, y_sigma2_on_mu2, corr_x_y)

end subroutine f2py_corr_nn2nl

subroutine f2py_corr_nn2ll(nz, ngrdcol, corr_x_y_n, sigma_x_n, sigma_y_n, &
    x_sigma2_on_mu2, y_sigma2_on_mu2, corr_x_y)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: corr_NN2LL

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    corr_x_y_n, sigma_x_n, sigma_y_n, x_sigma2_on_mu2, y_sigma2_on_mu2
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: corr_x_y

  call corr_NN2LL(nz, ngrdcol, corr_x_y_n, sigma_x_n, sigma_y_n, &
    x_sigma2_on_mu2, y_sigma2_on_mu2, corr_x_y)

end subroutine f2py_corr_nn2ll

subroutine f2py_compute_mean_binormal(mu_x_1, mu_x_2, mixt_frac, xm)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: compute_mean_binormal

  implicit none

  real(core_rknd), intent(in) :: mu_x_1, mu_x_2, mixt_frac
  real(core_rknd), intent(out) :: xm

  xm = compute_mean_binormal(mu_x_1, mu_x_2, mixt_frac)

end subroutine f2py_compute_mean_binormal

subroutine f2py_compute_variance_binormal(xm, mu_x_1, mu_x_2, stdev_x_1, stdev_x_2, &
    mixt_frac, xp2)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: compute_variance_binormal

  implicit none

  real(core_rknd), intent(in) :: xm, mu_x_1, mu_x_2, stdev_x_1, stdev_x_2, mixt_frac
  real(core_rknd), intent(out) :: xp2

  xp2 = compute_variance_binormal(xm, mu_x_1, mu_x_2, stdev_x_1, stdev_x_2, mixt_frac)

end subroutine f2py_compute_variance_binormal

subroutine f2py_calc_corr_chi_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_chi_i, &
    corr_rt_x_i, corr_thl_x_i, corr_chi_x_i)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: calc_corr_chi_x

  implicit none

  real(core_rknd), intent(in) :: &
    crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_chi_i, corr_rt_x_i, corr_thl_x_i
  real(core_rknd), intent(out) :: corr_chi_x_i

  corr_chi_x_i = calc_corr_chi_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_chi_i, &
    corr_rt_x_i, corr_thl_x_i)

end subroutine f2py_calc_corr_chi_x

subroutine f2py_calc_corr_eta_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_eta_i, &
    corr_rt_x_i, corr_thl_x_i, corr_eta_x_i)

  use clubb_precision, only: core_rknd
  use pdf_utilities, only: calc_corr_eta_x

  implicit none

  real(core_rknd), intent(in) :: &
    crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_eta_i, corr_rt_x_i, corr_thl_x_i
  real(core_rknd), intent(out) :: corr_eta_x_i

  corr_eta_x_i = calc_corr_eta_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_eta_i, &
    corr_rt_x_i, corr_thl_x_i)

end subroutine f2py_calc_corr_eta_x

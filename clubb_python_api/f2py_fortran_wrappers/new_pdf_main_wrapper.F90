! new_pdf_main_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_pdf_main

subroutine f2py_new_pdf_driver(nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, &
    wprtp, wpthlp, rtpthlp, clubb_params, skrt, skthl, mu_w_1, mu_w_2, mu_rt_1, &
    mu_rt_2, mu_thl_1, mu_thl_2, sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, &
    sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use new_pdf_main, only: new_pdf_driver
  use derived_type_storage, only: stored_pdf_implicit_coefs_terms

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, rtpthlp
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nz), intent(inout) :: skrt, skthl
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, &
    sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, &
    sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac

  call new_pdf_driver(nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, &
    rtpthlp, clubb_params, skrt, skthl, mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, &
    mu_thl_2, sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, &
    sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac, stored_pdf_implicit_coefs_terms)

end subroutine f2py_new_pdf_driver


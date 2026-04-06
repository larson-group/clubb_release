! new_hybrid_pdf_main_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module new_hybrid_pdf_main

subroutine f2py_new_hybrid_pdf_driver(nz, ngrdcol, sclr_dim, sclr_dim_transport, wm, rtm, thlm, um, vm, wp2, &
    rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp, upwp, vpwp, sclrm, sclrp2, wpsclrp, &
    gamma_skw_fnc, slope_coef_spread_dg_means_w, pdf_component_stdev_factor_w, skrt, &
    skthl, sku, skv, sksclr, mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, &
    mu_u_1, mu_u_2, mu_v_1, mu_v_2, sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, &
    sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, sigma_u_1_sqd, sigma_u_2_sqd, &
    sigma_v_1_sqd, sigma_v_2_sqd, mu_sclr_1, mu_sclr_2, sigma_sclr_1_sqd, &
    sigma_sclr_2_sqd, mixt_frac, sigma_sqd_w)

  use clubb_precision, only: core_rknd
  use new_hybrid_pdf_main, only: new_hybrid_pdf_driver
  use derived_type_storage, only: stored_pdf_implicit_coefs_terms

  implicit none

  integer, intent(in) :: nz, ngrdcol, sclr_dim, sclr_dim_transport
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: &
    wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp, upwp, vpwp, gamma_skw_fnc
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(in) :: sclrm, sclrp2, wpsclrp
  real(core_rknd), dimension(ngrdcol), intent(in) :: &
    slope_coef_spread_dg_means_w, pdf_component_stdev_factor_w
  real(core_rknd), dimension(ngrdcol, nz), intent(inout) :: skrt, skthl, sku, skv
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(inout) :: sksclr
  real(core_rknd), dimension(ngrdcol, nz), intent(out) :: &
    mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, mu_u_1, mu_u_2, mu_v_1, mu_v_2, &
    sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, &
    sigma_u_1_sqd, sigma_u_2_sqd, sigma_v_1_sqd, sigma_v_2_sqd, mixt_frac, sigma_sqd_w
  real(core_rknd), dimension(ngrdcol, nz, sclr_dim_transport), intent(out) :: &
    mu_sclr_1, mu_sclr_2, sigma_sclr_1_sqd, sigma_sclr_2_sqd

  call new_hybrid_pdf_driver(nz, ngrdcol, sclr_dim, wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, &
    skw, wprtp, wpthlp, upwp, vpwp, sclrm, sclrp2, wpsclrp, gamma_skw_fnc, &
    slope_coef_spread_dg_means_w, pdf_component_stdev_factor_w, skrt, skthl, sku, skv, sksclr, &
    mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, mu_u_1, mu_u_2, mu_v_1, mu_v_2, &
    sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd, sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, &
    sigma_u_1_sqd, sigma_u_2_sqd, sigma_v_1_sqd, sigma_v_2_sqd, mu_sclr_1, mu_sclr_2, &
    sigma_sclr_1_sqd, sigma_sclr_2_sqd, mixt_frac, sigma_sqd_w, stored_pdf_implicit_coefs_terms)

end subroutine f2py_new_hybrid_pdf_driver

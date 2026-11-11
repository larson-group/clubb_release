! Export compiled native constants for the C++ core boundary test.
subroutine f2k_core_constants_c(values, indices, flags) bind(C)
  ! Description: Return constant bits, one-based parameter indices and core options.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use constants_clubb, only: zero, two, w_tol, w_tol_sqd, em_min, zero_threshold
  use parameter_indices, only: iC_wp2_splat, ic_K, ilambda0_stability_coef, &
                              iup2_sfc_coef, ia_const, ibv_efold
  use model_flags, only: order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, order_windm, &
                        iiPDF_ADG1, ipdf_pre_advance_fields, ipdf_post_advance_fields, &
                        ipdf_pre_post_advance_fields, l_advance_xp3, Lscale_tau
  implicit none

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: values(6) ! Scalar constants [units vary]
  integer(c_int), intent(out) :: indices(6) ! Native parameter-column indices [-]
  integer(c_int), intent(out) :: flags(10) ! Native compile-time choices [-]

  !--------------------------- Begin Code ---------------------------
  values = [zero, two, w_tol, w_tol_sqd, em_min, zero_threshold]
  indices = [iC_wp2_splat, ic_K, ilambda0_stability_coef, iup2_sfc_coef, ia_const, ibv_efold]
  flags = [order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, order_windm, iiPDF_ADG1, &
           ipdf_pre_advance_fields, ipdf_post_advance_fields, ipdf_pre_post_advance_fields, &
           merge(1,0,l_advance_xp3), Lscale_tau]
end subroutine f2k_core_constants_c

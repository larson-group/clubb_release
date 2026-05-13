! sigma_sqd_w_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_compute_sigma_sqd_w(nzm, nzt, ngrdcol, &
    gamma_Skw_fnc, wp2, thlp2, rtp2, up2, vp2, &
    wpthlp, wprtp, upwp, vpwp, &
    l_predict_upwp_vpwp, sigma_sqd_w)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use sigma_sqd_w_module, only: compute_sigma_sqd_w

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    gamma_Skw_fnc, wp2, thlp2, rtp2, up2, vp2, &
    wpthlp, wprtp, upwp, vpwp
  logical, intent(in) :: l_predict_upwp_vpwp
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: sigma_sqd_w

  call compute_sigma_sqd_w(nzm, nzt, ngrdcol, stored_grid, &
    gamma_Skw_fnc, wp2, thlp2, rtp2, up2, vp2, &
    wpthlp, wprtp, upwp, vpwp, &
    l_predict_upwp_vpwp, &
    sigma_sqd_w)

end subroutine f2py_compute_sigma_sqd_w

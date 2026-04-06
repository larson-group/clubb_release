! sfc_varnce_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module sfc_varnce_module

subroutine f2py_calc_sfc_varnce(nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport, &
    dt, sfc_elevation, upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc, um, vm, lscale_up, &
    wpsclrp_sfc, lhs_splat_wp2, tau_zm, l_vary_convect_depth, t0, up2_sfc_coef, a_const, &
    wp2, up2, vp2, thlp2, rtp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_sclr_idx, stored_stats, stored_err_info
  use sfc_varnce_module, only: calc_sfc_varnce

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport
  logical, intent(in) :: l_vary_convect_depth
  real(core_rknd), intent(in) :: dt, t0
  real(core_rknd), dimension(ngrdcol), intent(in) :: sfc_elevation, upwp_sfc, vpwp_sfc, wprtp_sfc
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: wpthlp, lhs_splat_wp2, tau_zm
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: um, vm, lscale_up
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(in) :: wpsclrp_sfc
  real(core_rknd), dimension(ngrdcol), intent(in) :: up2_sfc_coef, a_const
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: wp2, up2, vp2, thlp2, rtp2, rtpthlp
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: sclrp2, sclrprtp, sclrpthlp

  call calc_sfc_varnce(nzm, nzt, ngrdcol, sclr_dim, stored_sclr_idx, stored_grid, dt, sfc_elevation, &
    upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc, um, vm, lscale_up, wpsclrp_sfc, lhs_splat_wp2, tau_zm, &
    l_vary_convect_depth, t0, up2_sfc_coef, a_const, stored_stats, wp2, up2, vp2, thlp2, &
    rtp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, stored_err_info)

end subroutine f2py_calc_sfc_varnce

! advance_xp3_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_advance_xp3(nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport, &
    sclr_tol, dt, &
    rtm, thlm, rtp2, thlp2, wprtp, wpthlp, wprtp2, wpthlp2, &
    rho_ds_zm, invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt, &
    sclrm, sclrp2, wpsclrp, wpsclrp2, &
    l_lmm_stepping, &
    rtp3, thlp3, sclrp3)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_stats
  use advance_xp3_module, only: advance_xp3

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, sclr_dim_transport
  real(core_rknd), dimension(sclr_dim_transport), intent(in) :: sclr_tol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    rtm, thlm, wprtp2, wpthlp2, invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    rtp2, thlp2, wprtp, wpthlp, rho_ds_zm
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(in) :: &
    sclrm, wpsclrp2
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(in) :: &
    sclrp2, wpsclrp
  logical, intent(in) :: l_lmm_stepping
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: rtp3, thlp3
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(inout) :: sclrp3

  call advance_xp3(nzm, nzt, ngrdcol, sclr_dim, sclr_tol, stored_grid, dt, &
    rtm, thlm, rtp2, thlp2, wprtp, wpthlp, wprtp2, wpthlp2, rho_ds_zm, &
    invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt, &
    sclrm, sclrp2, wpsclrp, wpsclrp2, l_lmm_stepping, &
    stored_stats, rtp3, thlp3, sclrp3)

end subroutine f2py_advance_xp3

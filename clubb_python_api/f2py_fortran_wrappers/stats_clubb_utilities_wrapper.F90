! stats_clubb_utilities_wrapper.F90 — wrapper for module stats_clubb_utilities

subroutine f2py_stats_accumulate( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, dt, &
    l_implemented, l_host_applies_sfc_fluxes, l_stability_correct_tau_zm, &
    um, vm, upwp, vpwp, up2, vp2, &
    thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, &
    wpthlp_sfc, wprtp_sfc, wprtp, wpthlp, &
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
    rtpthlp, &
    p_in_Pa, exner, rho, rho_zm, &
    rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    wm_zt, wm_zm, rcm, cloud_frac, &
    thvm, ug, vg, ddzt_umvm_sqd, stability_correction, &
    Kh_zt, rsat, Kh_zm, em, wp3_on_wp2, wp3_on_wp2_zt, &
    sclrm, sclrp2, &
    sclrprtp, sclrpthlp, sclrm_forcing, &
    wpsclrp, wpedsclrp, edsclrm, &
    edsclrm_forcing, &
    saturation_formula)

  use clubb_precision, only: core_rknd
  use stats_clubb_utilities, only: stats_accumulate
  use derived_type_storage, only: stored_grid, stored_stats

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, edsclr_dim
  integer, intent(in) :: saturation_formula
  logical, intent(in) :: l_implemented, l_host_applies_sfc_fluxes, l_stability_correct_tau_zm
  real(core_rknd), intent(in) :: dt

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    um, vm, thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, wp3, rtp3, thlp3

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    upwp, vpwp, up2, vp2, wprtp, wpthlp, wp2, rtp2, thlp2, rtpthlp

  real(core_rknd), intent(in), dimension(ngrdcol) :: wpthlp_sfc, wprtp_sfc

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    p_in_Pa, exner, rho, rho_ds_zt, thv_ds_zt, wm_zt

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    rho_zm, rho_ds_zm, thv_ds_zm, wm_zm

  real(core_rknd), intent(in), dimension(ngrdcol, nzt) :: &
    rcm, cloud_frac, thvm, ug, vg, Kh_zt, rsat, wp3_on_wp2_zt

  real(core_rknd), intent(in), dimension(ngrdcol, nzm) :: &
    ddzt_umvm_sqd, stability_correction, Kh_zm, em, wp3_on_wp2

  real(core_rknd), intent(in), dimension(ngrdcol, nzt, sclr_dim) :: &
    sclrm, sclrm_forcing

  real(core_rknd), intent(in), dimension(ngrdcol, nzm, sclr_dim) :: &
    sclrp2, sclrprtp, sclrpthlp, wpsclrp

  real(core_rknd), intent(in), dimension(ngrdcol, nzt, edsclr_dim) :: &
    edsclrm, edsclrm_forcing

  real(core_rknd), intent(in), dimension(ngrdcol, nzm, edsclr_dim) :: &
    wpedsclrp

  call stats_accumulate( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, &
    stored_grid, dt, &
    l_implemented, l_host_applies_sfc_fluxes, l_stability_correct_tau_zm, &
    um, vm, upwp, vpwp, up2, vp2, &
    thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, &
    wpthlp_sfc, wprtp_sfc, wprtp, wpthlp, &
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
    rtpthlp, &
    p_in_Pa, exner, rho, rho_zm, &
    rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
    wm_zt, wm_zm, rcm, cloud_frac, &
    thvm, ug, vg, ddzt_umvm_sqd, stability_correction, &
    Kh_zt, rsat, Kh_zm, em, wp3_on_wp2, wp3_on_wp2_zt, &
    sclrm, sclrp2, &
    sclrprtp, sclrpthlp, sclrm_forcing, &
    wpsclrp, wpedsclrp, edsclrm, &
    edsclrm_forcing, &
    saturation_formula, &
    stored_stats )

end subroutine f2py_stats_accumulate

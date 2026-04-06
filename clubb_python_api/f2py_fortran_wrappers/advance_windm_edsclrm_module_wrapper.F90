! advance_windm_edsclrm_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_advance_windm_edsclrm(nzm, nzt, ngrdcol, edsclr_dim, edsclr_dim_transport, dt, &
    wm_zt, km_zm, kmh_zm, &
    ug, vg, um_ref, vm_ref, &
    wp2, up2, vp2, um_forcing, vm_forcing, &
    edsclrm_forcing, &
    rho_ds_zm, invrs_rho_ds_zt, &
    fcor, l_implemented, ts_nudge, tridiag_solve_method, &
    l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge, l_tke_aniso, &
    l_lmm_stepping, l_linearize_pbl_winds, &
    order_xp2_xpyp, order_wp2_wp3, order_windm, &
    um, vm, edsclrm, upwp, vpwp, wpedsclrp, &
    um_pert, vm_pert, upwp_pert, vpwp_pert)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: &
    stored_grid, stored_nu_vert_res_dep, stored_stats, stored_err_info
  use advance_windm_edsclrm_module, only: advance_windm_edsclrm

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, edsclr_dim, edsclr_dim_transport
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    wm_zt, ug, vg, um_ref, vm_ref, um_forcing, vm_forcing, invrs_rho_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    km_zm, kmh_zm, wp2, up2, vp2, rho_ds_zm
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(in) :: &
    edsclrm_forcing
  real(core_rknd), dimension(ngrdcol), intent(in) :: fcor
  logical, intent(in) :: l_implemented
  real(core_rknd), intent(in) :: ts_nudge
  integer, intent(in) :: tridiag_solve_method
  logical, intent(in) :: &
    l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge, l_tke_aniso, l_lmm_stepping, &
    l_linearize_pbl_winds
  integer, intent(in) :: order_xp2_xpyp, order_wp2_wp3, order_windm

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: um, vm, um_pert, vm_pert
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: upwp, vpwp, upwp_pert, vpwp_pert
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(inout) :: edsclrm
  real(core_rknd), dimension(ngrdcol, nzm, edsclr_dim_transport), intent(inout) :: wpedsclrp

  call advance_windm_edsclrm(nzm, nzt, ngrdcol, edsclr_dim, stored_grid, dt, &
    wm_zt, km_zm, kmh_zm, ug, vg, um_ref, vm_ref, &
    wp2, up2, vp2, um_forcing, vm_forcing, edsclrm_forcing, &
    rho_ds_zm, invrs_rho_ds_zt, fcor, l_implemented, &
    stored_nu_vert_res_dep, ts_nudge, tridiag_solve_method, &
    l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge, &
    l_tke_aniso, l_lmm_stepping, l_linearize_pbl_winds, &
    order_xp2_xpyp, order_wp2_wp3, order_windm, &
    stored_stats, um, vm, edsclrm, upwp, vpwp, wpedsclrp, &
    um_pert, vm_pert, upwp_pert, vpwp_pert, stored_err_info)

end subroutine f2py_advance_windm_edsclrm

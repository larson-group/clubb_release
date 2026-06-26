! advance_windm_edsclrm_module_wrapper.F90 — wrappers organized by source module

subroutine f2py_advance_windm_edsclrm(nzm, nzt, ngrdcol, edsclr_dim, edsclr_dim_transport, dt, &
    wm_zt, kh_zm, clubb_params, &
    ug, vg, um_ref, vm_ref, &
    wp2, up2, vp2, um_forcing, vm_forcing, &
    edsclrm_forcing, p_in_pa, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
    fcor, l_implemented, ts_nudge, tridiag_solve_method, &
    l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge, l_tke_aniso, &
    l_lmm_stepping, l_linearize_pbl_winds, l_do_expldiff_rtm_thlm, &
    fill_holes_type, &
    upwp_cl_num, vpwp_cl_num, &
    um, vm, thlm, rtm, edsclrm, upwp, vpwp, wpedsclrp, &
    um_pert, vm_pert, upwp_pert, vpwp_pert)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: &
    stored_grid, stored_nu_vert_res_dep, stored_stats, stored_err_info
  use advance_windm_edsclrm_module, only: advance_windm_edsclrm

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, edsclr_dim, edsclr_dim_transport
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    wm_zt, ug, vg, um_ref, vm_ref, um_forcing, vm_forcing, &
    p_in_pa, rho_ds_zt, invrs_rho_ds_zt
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: &
    kh_zm, wp2, up2, vp2, rho_ds_zm
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(in) :: &
    edsclrm_forcing
  real(core_rknd), dimension(ngrdcol), intent(in) :: fcor
  logical, intent(in) :: l_implemented
  real(core_rknd), intent(in) :: ts_nudge
  integer, intent(in) :: tridiag_solve_method
  logical, intent(in) :: &
    l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge, l_tke_aniso, l_lmm_stepping, &
    l_linearize_pbl_winds, l_do_expldiff_rtm_thlm
  integer, intent(in) :: fill_holes_type
  integer, intent(inout) :: upwp_cl_num, vpwp_cl_num

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: um, vm, thlm, rtm, um_pert, vm_pert
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: upwp, vpwp, upwp_pert, vpwp_pert
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(inout) :: edsclrm
  real(core_rknd), dimension(ngrdcol, nzm, edsclr_dim_transport), intent(inout) :: wpedsclrp

  call advance_windm_edsclrm(nzm, nzt, ngrdcol, edsclr_dim, stored_grid, dt, &
    wm_zt, kh_zm, clubb_params, ug, vg, um_ref, vm_ref, &
    wp2, up2, vp2, um_forcing, vm_forcing, edsclrm_forcing, p_in_pa, &
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, fcor, l_implemented, &
    stored_nu_vert_res_dep, ts_nudge, tridiag_solve_method, &
    l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge, &
    l_tke_aniso, l_lmm_stepping, l_linearize_pbl_winds, &
    l_do_expldiff_rtm_thlm, fill_holes_type, &
    upwp_cl_num, vpwp_cl_num, &
    stored_stats, um, vm, thlm, rtm, edsclrm, upwp, vpwp, wpedsclrp, &
    um_pert, vm_pert, upwp_pert, vpwp_pert, stored_err_info)

end subroutine f2py_advance_windm_edsclrm

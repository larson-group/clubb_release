! prescribe_forcings_module_wrapper.F90 — F2PY wrapper for prescribe_forcings_module
!
! Standalone subroutine (NOT in a module) so it gets standard Fortran
! name mangling (f2py_prescribe_forcings_) which F2PY expects.
! UDT arguments come from derived_type_storage.

subroutine f2py_prescribe_forcings( &
    nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, sclr_dim_transport, edsclr_dim_transport, &
    runtype, sfctype, time_current, time_initial, dt, &
    um, vm, thlm, p_in_pa, exner, rho, rho_zm, thvm, zt_in, &
    l_t_dependent, l_ignore_forcings, l_input_xpwp_sfc, &
    l_modify_bc_for_cnvg_test, &
    saturation_formula, l_add_dycore_grid, grid_remap_method, &
    grid_adapt_in_time_method, &
    rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
    wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
    wpsclrp, sclrm_forcing, edsclrm_forcing, &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
    t_sfc, p_sfc, sens_ht, latent_ht, &
    wpsclrp_sfc, wpedsclrp_sfc)

  use clubb_precision, only: core_rknd, time_precision
  use grid_class, only: grid
  use error_code, only: clubb_fatal_error
  use constants_clubb, only: fstderr
  use prescribe_forcings_module, only: prescribe_forcings
  use time_dependent_input, only: &
    initialize_t_dependent_input, &
    l_t_dependent_mod => l_t_dependent, &
    l_ignore_forcings_mod => l_ignore_forcings, &
    l_input_xpwp_sfc_mod => l_input_xpwp_sfc
  use derived_type_storage, only: stored_grid, stored_sclr_idx, stored_stats, stored_err_info

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, sclr_dim, edsclr_dim
  integer, intent(in) :: sclr_dim_transport, edsclr_dim_transport
  character(len=50), intent(in) :: runtype
  integer, intent(in) :: sfctype
  real(time_precision), intent(in) :: time_current, time_initial
  real(core_rknd), intent(in) :: dt
  logical, intent(in) :: l_t_dependent, l_ignore_forcings, l_input_xpwp_sfc
  logical, intent(in) :: l_modify_bc_for_cnvg_test
  integer, intent(in) :: saturation_formula
  logical, intent(in) :: l_add_dycore_grid
  integer, intent(in) :: grid_remap_method, grid_adapt_in_time_method

  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    um, vm, thlm, p_in_pa, exner, rho, thvm, zt_in
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: rho_zm

  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: &
    rtm, wm_zt, ug, vg, um_ref, vm_ref, &
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: &
    wm_zm, wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, rtpthlp_forcing
  real(core_rknd), dimension(ngrdcol, nzm, sclr_dim_transport), intent(inout) :: wpsclrp
  real(core_rknd), dimension(ngrdcol, nzt, sclr_dim_transport), intent(inout) :: sclrm_forcing
  real(core_rknd), dimension(ngrdcol, nzt, edsclr_dim_transport), intent(inout) :: edsclrm_forcing

  real(core_rknd), dimension(ngrdcol), intent(inout) :: &
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, t_sfc, p_sfc
  real(core_rknd), intent(inout) :: sens_ht, latent_ht
  real(core_rknd), dimension(ngrdcol, sclr_dim_transport), intent(inout) :: wpsclrp_sfc
  real(core_rknd), dimension(ngrdcol, edsclr_dim_transport), intent(inout) :: wpedsclrp_sfc
  real(core_rknd), dimension(ngrdcol) :: veg_t_in_k
  logical, save :: l_tdep_initialized = .false.
  integer, parameter :: iunit = 10
  integer :: total_idx_rho_lin_spline
  type(grid) :: gr_dycore
  real(core_rknd), dimension(ngrdcol, 1) :: rho_lin_spline_vals, rho_lin_spline_levels

  total_idx_rho_lin_spline = 1
  rho_lin_spline_vals = 0.0_core_rknd
  rho_lin_spline_levels = 0.0_core_rknd
  ! No vegetation-temperature input is currently exposed through the F2PY API.
  ! Use surface temperature as a placeholder to satisfy the updated interface.
  veg_t_in_k = t_sfc

  l_t_dependent_mod = l_t_dependent
  l_ignore_forcings_mod = l_ignore_forcings
  l_input_xpwp_sfc_mod = l_input_xpwp_sfc

  if ( l_t_dependent .and. (.not. l_tdep_initialized) ) then
    call initialize_t_dependent_input( &
        iunit, trim(runtype), nzt, zt_in(1,:), p_in_pa(1,:), &
        l_add_dycore_grid, grid_adapt_in_time_method )
    l_tdep_initialized = .true.
  end if

  if (l_add_dycore_grid) then
    write(fstderr, *) "f2py_prescribe_forcings: l_add_dycore_grid is not supported yet."
    stored_err_info%err_code = clubb_fatal_error
    return
  end if

  call prescribe_forcings( &
      stored_grid, nzm, nzt, ngrdcol, &
      sclr_dim, edsclr_dim, stored_sclr_idx, &
      runtype, sfctype, &
      time_current, time_initial, dt, &
      um, vm, thlm, &
      p_in_pa, exner, rho, rho_zm, thvm, &
      veg_t_in_k, &
      l_modify_bc_for_cnvg_test, &
      saturation_formula, &
      stored_stats, &
      l_add_dycore_grid, &
      grid_remap_method, &
      total_idx_rho_lin_spline, rho_lin_spline_vals, rho_lin_spline_levels, &
      gr_dycore, &
      rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref, &
      thlm_forcing, rtm_forcing, um_forcing, &
      vm_forcing, wprtp_forcing, wpthlp_forcing, &
      rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
      wpsclrp, &
      sclrm_forcing, &
      edsclrm_forcing, &
      wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
      t_sfc, p_sfc, sens_ht, latent_ht, &
      wpsclrp_sfc, &
      wpedsclrp_sfc, &
      stored_err_info )

end subroutine f2py_prescribe_forcings

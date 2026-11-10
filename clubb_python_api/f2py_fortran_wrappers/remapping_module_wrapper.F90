! remapping_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module remapping_module

subroutine f2py_remap_vals_to_target_same_grid( nzm, nzt, ngrdcol, source_values_idx, & ! In
                                                source_values, target_values_idx, & ! In
                                                total_idx_rho_lin_spline, & ! In
                                                rho_lin_spline_vals, & ! In
                                                rho_lin_spline_levels, iv, p_sfc, & ! In
                                                grid_remap_method, l_zt_variable, & ! In
                                                l_apply_remap_mask, l_remap, & ! In
                                                target_values ) ! Out

  ! Description:
  ! Reconstruct the stored grid and apply optional column selection at the F2PY boundary.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: grid
  use remapping_module, only: remap_vals_to_target

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, source_values_idx, target_values_idx, &
    total_idx_rho_lin_spline
  real(core_rknd), dimension(ngrdcol, source_values_idx), intent(in) :: source_values
  real(core_rknd), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
    rho_lin_spline_vals, rho_lin_spline_levels
  integer, intent(in) :: iv, grid_remap_method
  logical, intent(in) :: l_zt_variable
  logical, intent(in) :: l_apply_remap_mask ! Whether the optional selection was supplied [-]
  logical, intent(in) :: l_remap(ngrdcol) ! Columns selected for remapping [-]
  real(core_rknd), dimension(ngrdcol), intent(in) :: p_sfc
  real(core_rknd), dimension(ngrdcol, target_values_idx), intent(out) :: target_values

  type(grid) :: gr_source, gr_target

  gr_source = stored_grid
  gr_target = stored_grid

  if ( l_apply_remap_mask ) then

    target_values = remap_vals_to_target( ngrdcol, & ! In
                                          gr_source, gr_target, & ! In
                                          source_values_idx, & ! In
                                          source_values, & ! In
                                          target_values_idx, & ! In
                                          total_idx_rho_lin_spline, & ! In
                                          rho_lin_spline_vals, & ! In
                                          rho_lin_spline_levels, & ! In
                                          iv, p_sfc, & ! In
                                          grid_remap_method, & ! In
                                          l_zt_variable, & ! In
                                          l_remap ) ! Optional In

  else

    target_values = remap_vals_to_target( ngrdcol, & ! In
                                          gr_source, gr_target, & ! In
                                          source_values_idx, & ! In
                                          source_values, & ! In
                                          target_values_idx, & ! In
                                          total_idx_rho_lin_spline, & ! In
                                          rho_lin_spline_vals, & ! In
                                          rho_lin_spline_levels, & ! In
                                          iv, p_sfc, & ! In
                                          grid_remap_method, & ! In
                                          l_zt_variable ) ! In

  end if

end subroutine f2py_remap_vals_to_target_same_grid

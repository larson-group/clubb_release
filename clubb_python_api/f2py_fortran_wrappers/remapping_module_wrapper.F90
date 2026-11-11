! remapping_module_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module remapping_module

subroutine f2py_remap_vals_to_target_same_grid(nzm, nzt, ngrdcol, source_values_idx, &
    source_values, target_values_idx, total_idx_rho_lin_spline, &
    rho_lin_spline_vals, rho_lin_spline_levels, iv, p_sfc, grid_remap_method, &
    l_zt_variable, target_values)

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
  real(core_rknd), dimension(ngrdcol), intent(in) :: p_sfc
  real(core_rknd), dimension(ngrdcol, target_values_idx), intent(out) :: target_values

  type(grid) :: gr_source, gr_target

  gr_source = stored_grid
  gr_target = stored_grid

  target_values = remap_vals_to_target(ngrdcol, gr_source, gr_target, source_values_idx, &
    source_values, target_values_idx, total_idx_rho_lin_spline, rho_lin_spline_vals, &
    rho_lin_spline_levels, iv, p_sfc, grid_remap_method, l_zt_variable)

end subroutine f2py_remap_vals_to_target_same_grid

! fill_holes_wrapper.F90 — wrappers organized by source module

subroutine f2py_fill_holes_vertical(nz, ngrdcol, threshold, &
    lower_hf_level, upper_hf_level, dz, rho_ds, &
    grid_dir_indx, fill_holes_type, field)

  use clubb_precision, only: core_rknd
  use fill_holes, only: fill_holes_vertical_api

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), intent(in) :: threshold
  integer, intent(in) :: lower_hf_level, upper_hf_level
  real(core_rknd), dimension(ngrdcol, nz), intent(in) :: dz, rho_ds
  integer, intent(in) :: grid_dir_indx, fill_holes_type
  real(core_rknd), dimension(ngrdcol, nz), intent(inout) :: field

  call fill_holes_vertical_api(nz, ngrdcol, threshold, &
    lower_hf_level, upper_hf_level, dz, rho_ds, &
    grid_dir_indx, fill_holes_type, field)

end subroutine f2py_fill_holes_vertical

subroutine f2py_fill_holes_wp2_from_horz_tke(nz, ngrdcol, threshold, &
    lower_hf_level, upper_hf_level, wp2, up2, vp2)

  use clubb_precision, only: core_rknd
  use fill_holes, only: fill_holes_wp2_from_horz_tke

  implicit none

  integer, intent(in) :: nz, ngrdcol
  real(core_rknd), intent(in) :: threshold
  integer, intent(in) :: lower_hf_level, upper_hf_level
  real(core_rknd), dimension(ngrdcol, nz), intent(inout) :: wp2, up2, vp2

  call fill_holes_wp2_from_horz_tke(nz, ngrdcol, threshold, &
    lower_hf_level, upper_hf_level, wp2, up2, vp2)

end subroutine f2py_fill_holes_wp2_from_horz_tke

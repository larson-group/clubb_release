! grid_class_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module grid_class

subroutine f2py_cleanup_grid()

  use derived_type_storage, only: stored_grid
  use grid_class, only: cleanup_grid_api

  implicit none

  call cleanup_grid_api(stored_grid)

end subroutine f2py_cleanup_grid

subroutine f2py_zm2zt2zm_2d(nzm, nzt, ngrdcol, azm, zm_min, azm_smoothed)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: zm2zt2zm

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: azm
  real(core_rknd), intent(in) :: zm_min
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: azm_smoothed

  azm_smoothed = zm2zt2zm(nzm, nzt, ngrdcol, stored_grid, azm, zm_min)

end subroutine f2py_zm2zt2zm_2d

subroutine f2py_zt2zm2zt_2d(nzm, nzt, ngrdcol, azt, zt_min, azt_smoothed)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: zt2zm2zt

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: azt
  real(core_rknd), intent(in) :: zt_min
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: azt_smoothed

  azt_smoothed = zt2zm2zt(nzm, nzt, ngrdcol, stored_grid, azt, zt_min)

end subroutine f2py_zt2zm2zt_2d

subroutine f2py_setup_grid(nzmax, ngrdcol, sfc_elevation, &
    l_implemented, l_ascending_grid, grid_type, &
    deltaz, zm_init, zm_top, &
    momentum_heights, thermodynamic_heights) &
  bind(C, name="f2py_setup_grid_")

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_err_info
  use grid_class, only: setup_grid, cleanup_grid_api

  implicit none

  integer, intent(in) :: nzmax, ngrdcol
  real(core_rknd), dimension(ngrdcol), intent(in) :: &
    sfc_elevation, deltaz, zm_init, zm_top
  integer, intent(in) :: grid_type
  logical, intent(in) :: l_implemented, l_ascending_grid
  real(core_rknd), dimension(ngrdcol, nzmax), intent(in) :: momentum_heights
  real(core_rknd), dimension(ngrdcol, nzmax-1), intent(in) :: thermodynamic_heights

  ! Cleanup previously allocated grid if any
  if (allocated(stored_grid%zm)) then
    call cleanup_grid_api(stored_grid)
  end if

  call setup_grid(nzmax, ngrdcol, sfc_elevation, &
                  l_implemented, l_ascending_grid, &
                  grid_type, deltaz, zm_init, zm_top, &
                  momentum_heights, thermodynamic_heights, &
                  stored_grid, stored_err_info)

end subroutine f2py_setup_grid

subroutine f2py_get_grid_dims(nzm_out, nzt_out) &
  bind(C, name="f2py_get_grid_dims_")

  use derived_type_storage, only: stored_grid

  implicit none

  integer, intent(out) :: nzm_out, nzt_out
  nzm_out = stored_grid%nzm
  nzt_out = stored_grid%nzt

end subroutine f2py_get_grid_dims

subroutine f2py_get_grid_shape(ngrdcol_out, nzm_out, nzt_out) &
  bind(C, name="f2py_get_grid_shape_")

  use derived_type_storage, only: stored_grid

  implicit none

  integer, intent(out) :: ngrdcol_out, nzm_out, nzt_out

  if (allocated(stored_grid%zm)) then
    ngrdcol_out = size(stored_grid%zm, 1)
    nzm_out = stored_grid%nzm
    nzt_out = stored_grid%nzt
  else
    ngrdcol_out = 0
    nzm_out = 0
    nzt_out = 0
  end if

end subroutine f2py_get_grid_shape

subroutine f2py_zt2zm_2d(nzm, nzt, ngrdcol, azt, azm)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: zt2zm_api

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: azt
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: azm

  azm = zt2zm_api(nzm, nzt, ngrdcol, stored_grid, azt)

end subroutine f2py_zt2zm_2d

subroutine f2py_zm2zt_2d(nzm, nzt, ngrdcol, azm, azt)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: zm2zt_api

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: azm
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: azt

  azt = zm2zt_api(nzm, nzt, ngrdcol, stored_grid, azm)

end subroutine f2py_zm2zt_2d

subroutine f2py_ddzm_2d(nzm, nzt, ngrdcol, azm, d_azm_dzt)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: ddzm

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: azm
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: d_azm_dzt

  d_azm_dzt = ddzm(nzm, nzt, ngrdcol, stored_grid, azm)

end subroutine f2py_ddzm_2d

subroutine f2py_ddzt_2d(nzm, nzt, ngrdcol, azt, d_azt_dzm)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use grid_class, only: ddzt

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: azt
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: d_azt_dzm

  d_azt_dzm = ddzt(nzm, nzt, ngrdcol, stored_grid, azt)

end subroutine f2py_ddzt_2d

subroutine f2py_setup_grid_heights(nzm, nzt, ngrdcol, &
    l_implemented, l_ascending_grid, grid_type, deltaz, zm_init, &
    momentum_heights, thermodynamic_heights)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_err_info
  use grid_class, only: setup_grid_heights

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol, grid_type
  logical, intent(in) :: l_implemented, l_ascending_grid
  real(core_rknd), dimension(ngrdcol), intent(in) :: deltaz, zm_init
  real(core_rknd), dimension(ngrdcol, nzm), intent(in) :: momentum_heights
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: thermodynamic_heights

  call setup_grid_heights(nzm, nzt, ngrdcol, l_implemented, l_ascending_grid, &
    grid_type, deltaz, zm_init, momentum_heights, thermodynamic_heights, &
    stored_grid, stored_err_info)

end subroutine f2py_setup_grid_heights

subroutine f2py_read_grid_heights(nzmax, grid_type, zm_grid_fname, zt_grid_fname, &
    file_unit, momentum_heights, thermodynamic_heights)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use grid_class, only: read_grid_heights

  implicit none

  integer, intent(in) :: nzmax, grid_type, file_unit
  character(len=*), intent(in) :: zm_grid_fname, zt_grid_fname
  real(core_rknd), dimension(nzmax), intent(out) :: momentum_heights
  real(core_rknd), dimension(nzmax-1), intent(out) :: thermodynamic_heights

  call read_grid_heights(nzmax, grid_type, zm_grid_fname, zt_grid_fname, file_unit, &
    momentum_heights, thermodynamic_heights, stored_err_info)

end subroutine f2py_read_grid_heights

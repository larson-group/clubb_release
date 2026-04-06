subroutine f2py_pos_definite_adj(nzm, nzt, ngrdcol, dt, field_np1, flux_np1, field_n, &
    field_pd, flux_pd)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid
  use pos_definite_module, only: pos_definite_adj

  implicit none

  integer, intent(in) :: nzm, nzt, ngrdcol
  real(core_rknd), intent(in) :: dt
  real(core_rknd), dimension(ngrdcol, nzt), intent(inout) :: field_np1
  real(core_rknd), dimension(ngrdcol, nzm), intent(inout) :: flux_np1
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: field_n
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: field_pd
  real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: flux_pd

  call pos_definite_adj(nzm, nzt, ngrdcol, stored_grid, dt, field_np1, flux_np1, field_n, &
    field_pd, flux_pd)

end subroutine f2py_pos_definite_adj

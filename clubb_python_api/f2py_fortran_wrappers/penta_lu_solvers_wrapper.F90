subroutine f2py_penta_lu_solve_single_rhs_multiple_lhs(ndim, ngrdcol, lhs, rhs, soln)

  use clubb_precision, only: core_rknd
  use penta_lu_solvers, only: penta_lu_solve

  implicit none

  integer, intent(in) :: ndim, ngrdcol
  real(core_rknd), dimension(5, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim), intent(out) :: soln

  real(core_rknd), dimension(-2:2, ngrdcol, ndim) :: lhs_local

  lhs_local(-2,:,:) = lhs(1,:,:)
  lhs_local(-1,:,:) = lhs(2,:,:)
  lhs_local( 0,:,:) = lhs(3,:,:)
  lhs_local( 1,:,:) = lhs(4,:,:)
  lhs_local( 2,:,:) = lhs(5,:,:)

  call penta_lu_solve(ndim, ngrdcol, lhs_local, rhs, soln)

end subroutine f2py_penta_lu_solve_single_rhs_multiple_lhs

subroutine f2py_penta_lu_solve_multiple_rhs_lhs(ndim, nrhs, ngrdcol, lhs, rhs, soln)

  use clubb_precision, only: core_rknd
  use penta_lu_solvers, only: penta_lu_solve

  implicit none

  integer, intent(in) :: ndim, nrhs, ngrdcol
  real(core_rknd), dimension(5, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln

  real(core_rknd), dimension(-2:2, ngrdcol, ndim) :: lhs_local

  lhs_local(-2,:,:) = lhs(1,:,:)
  lhs_local(-1,:,:) = lhs(2,:,:)
  lhs_local( 0,:,:) = lhs(3,:,:)
  lhs_local( 1,:,:) = lhs(4,:,:)
  lhs_local( 2,:,:) = lhs(5,:,:)

  call penta_lu_solve(ndim, nrhs, ngrdcol, lhs_local, rhs, soln)

end subroutine f2py_penta_lu_solve_multiple_rhs_lhs

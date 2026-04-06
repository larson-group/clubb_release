subroutine f2py_penta_bicgstab_solve(ndim, ngrdcol, lhs, rhs, soln, iters)

  use clubb_precision, only: core_rknd
  use bicgstab_solvers, only: penta_bicgstab_solve

  implicit none

  integer, intent(in) :: ndim, ngrdcol
  real(core_rknd), dimension(5, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim), intent(out) :: soln
  integer, intent(out) :: iters

  real(core_rknd), dimension(-2:2, ngrdcol, ndim) :: lhs_local

  lhs_local(-2,:,:) = lhs(1,:,:)
  lhs_local(-1,:,:) = lhs(2,:,:)
  lhs_local( 0,:,:) = lhs(3,:,:)
  lhs_local( 1,:,:) = lhs(4,:,:)
  lhs_local( 2,:,:) = lhs(5,:,:)

  call penta_bicgstab_solve(ndim, ngrdcol, lhs_local, rhs, soln, iters=iters)

end subroutine f2py_penta_bicgstab_solve

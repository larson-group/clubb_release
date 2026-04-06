! matrix_solver_wrapper_wrapper.F90 - wrappers for matrix solver interfaces

subroutine f2py_band_solve_multiple_rhs(solve_name, penta_solve_method, &
    ngrdcol, nsup, nsub, ndim, nrhs, lhs, rhs, old_soln, use_rcond, soln, rcond)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use matrix_solver_wrapper, only: band_solve

  implicit none

  character(len=*), intent(in) :: solve_name
  integer, intent(in) :: penta_solve_method, ngrdcol, nsup, nsub, ndim, nrhs
  logical, intent(in) :: use_rcond
  real(core_rknd), dimension(nsup+nsub+1, ngrdcol, ndim), intent(inout) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(inout) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(in) :: old_soln
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln
  real(core_rknd), dimension(ngrdcol), intent(out) :: rcond

  if (use_rcond) then
    call band_solve(solve_name, penta_solve_method, ngrdcol, nsup, nsub, ndim, nrhs, &
                    lhs, rhs, stored_err_info, soln, old_soln=old_soln, rcond=rcond)
  else
    call band_solve(solve_name, penta_solve_method, ngrdcol, nsup, nsub, ndim, nrhs, &
                    lhs, rhs, stored_err_info, soln, old_soln=old_soln)
    rcond = 0.0_core_rknd
  end if

end subroutine f2py_band_solve_multiple_rhs

subroutine f2py_tridiag_solve_single_rhs_multiple_lhs(solve_name, tridiag_solve_method, &
    ngrdcol, ndim, lhs, rhs, use_rcond, soln, rcond)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use matrix_solver_wrapper, only: tridiag_solve

  implicit none

  character(len=*), intent(in) :: solve_name
  integer, intent(in) :: tridiag_solve_method, ngrdcol, ndim
  logical, intent(in) :: use_rcond
  real(core_rknd), dimension(3, ngrdcol, ndim), intent(inout) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim), intent(inout) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim), intent(out) :: soln
  real(core_rknd), dimension(ngrdcol), intent(out) :: rcond

  if (use_rcond) then
    call tridiag_solve(solve_name, tridiag_solve_method, ngrdcol, ndim, &
                       lhs, rhs, stored_err_info, soln, rcond)
  else
    call tridiag_solve(solve_name, tridiag_solve_method, ngrdcol, ndim, &
                       lhs, rhs, stored_err_info, soln)
    rcond = 0.0_core_rknd
  end if

end subroutine f2py_tridiag_solve_single_rhs_multiple_lhs

subroutine f2py_tridiag_solve_multiple_rhs(solve_name, tridiag_solve_method, &
    ngrdcol, ndim, nrhs, lhs, rhs, use_rcond, soln, rcond)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use matrix_solver_wrapper, only: tridiag_solve

  implicit none

  character(len=*), intent(in) :: solve_name
  integer, intent(in) :: tridiag_solve_method, ngrdcol, ndim, nrhs
  logical, intent(in) :: use_rcond
  real(core_rknd), dimension(3, ngrdcol, ndim), intent(inout) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(inout) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln
  real(core_rknd), dimension(ngrdcol), intent(out) :: rcond

  if (use_rcond) then
    call tridiag_solve(solve_name, tridiag_solve_method, ngrdcol, ndim, nrhs, &
                       lhs, rhs, stored_err_info, soln, rcond)
  else
    call tridiag_solve(solve_name, tridiag_solve_method, ngrdcol, ndim, nrhs, &
                       lhs, rhs, stored_err_info, soln)
    rcond = 0.0_core_rknd
  end if

end subroutine f2py_tridiag_solve_multiple_rhs

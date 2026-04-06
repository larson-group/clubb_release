subroutine f2py_lapack_tridiag_solve(solve_type, ndim, nrhs, ngrdcol, lhs, rhs, soln)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use lapack_wrap, only: lapack_tridiag_solve

  implicit none

  character(len=*), intent(in) :: solve_type
  integer, intent(in) :: ndim, nrhs, ngrdcol
  real(core_rknd), dimension(3, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln

  real(core_rknd), dimension(3, ngrdcol, ndim) :: lhs_local
  real(core_rknd), dimension(ngrdcol, ndim, nrhs) :: rhs_local

  lhs_local = lhs
  rhs_local = rhs

  call lapack_tridiag_solve(solve_type, ndim, nrhs, ngrdcol, lhs_local, rhs_local, stored_err_info, soln)

end subroutine f2py_lapack_tridiag_solve

subroutine f2py_lapack_tridiag_solvex(solve_type, ndim, nrhs, ngrdcol, lhs, rhs, soln, rcond)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use lapack_wrap, only: lapack_tridiag_solvex

  implicit none

  character(len=*), intent(in) :: solve_type
  integer, intent(in) :: ndim, nrhs, ngrdcol
  real(core_rknd), dimension(3, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln
  real(core_rknd), dimension(ngrdcol), intent(out) :: rcond

  real(core_rknd), dimension(3, ngrdcol, ndim) :: lhs_local
  real(core_rknd), dimension(ngrdcol, ndim, nrhs) :: rhs_local

  lhs_local = lhs
  rhs_local = rhs

  call lapack_tridiag_solvex(solve_type, ndim, nrhs, ngrdcol, lhs_local, rhs_local, stored_err_info, soln, rcond)

end subroutine f2py_lapack_tridiag_solvex

subroutine f2py_lapack_band_solve(solve_type, nsup, nsub, ndim, nrhs, ngrdcol, lhs, rhs, soln)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use lapack_wrap, only: lapack_band_solve

  implicit none

  character(len=*), intent(in) :: solve_type
  integer, intent(in) :: nsup, nsub, ndim, nrhs, ngrdcol
  real(core_rknd), dimension(nsup+nsub+1, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln

  real(core_rknd), dimension(nsup+nsub+1, ngrdcol, ndim) :: lhs_local
  real(core_rknd), dimension(ngrdcol, ndim, nrhs) :: rhs_local

  lhs_local = lhs
  rhs_local = rhs

  call lapack_band_solve(solve_type, nsup, nsub, ndim, nrhs, ngrdcol, lhs_local, rhs_local, stored_err_info, soln)

end subroutine f2py_lapack_band_solve

subroutine f2py_lapack_band_solvex(solve_type, nsup, nsub, ndim, nrhs, ngrdcol, lhs, rhs, soln, rcond)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_err_info
  use lapack_wrap, only: lapack_band_solvex

  implicit none

  character(len=*), intent(in) :: solve_type
  integer, intent(in) :: nsup, nsub, ndim, nrhs, ngrdcol
  real(core_rknd), dimension(nsup+nsub+1, ngrdcol, ndim), intent(in) :: lhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(in) :: rhs
  real(core_rknd), dimension(ngrdcol, ndim, nrhs), intent(out) :: soln
  real(core_rknd), dimension(ngrdcol), intent(out) :: rcond

  real(core_rknd), dimension(nsup+nsub+1, ngrdcol, ndim) :: lhs_local
  real(core_rknd), dimension(ngrdcol, ndim, nrhs) :: rhs_local

  lhs_local = lhs
  rhs_local = rhs

  call lapack_band_solvex(solve_type, nsup, nsub, ndim, nrhs, ngrdcol, lhs_local, rhs_local, stored_err_info, soln, rcond)

end subroutine f2py_lapack_band_solvex

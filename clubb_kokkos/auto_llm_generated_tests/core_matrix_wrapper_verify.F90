! Exercise the native generic matrix wrapper, including optional conditions and error gates.
subroutine f2k_matrix_reference_c(variant, solve_method, ndim, nrhs, ngrdcol, use_rcond, &
                                   debug_level, implemented, lhs, rhs, soln, rcond, codes) bind(C)
  ! Description: Reference all5 native solve surfaces with matching debug/host gates.
  ! Seed an existing nonfatal code to check that successful solves preserve status.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use err_info_type_module, only: err_info_type
  use error_code, only: set_clubb_debug_level_api
  use matrix_solver_wrapper, only: band_solve, tridiag_solve
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: variant, solve_method, ndim, nrhs, ngrdcol ! Solve shape
  integer(c_int), value, intent(in) :: use_rcond, debug_level, implemented ! Diagnostic gates

  !--------------------------- Input/Output Variables ---------------------------
  real( kind = core_rknd ), intent(inout) :: lhs(merge(5,3,variant<=2),ngrdcol,ndim) ! Bands
  real( kind = core_rknd ), intent(inout) :: rhs(ngrdcol,ndim,nrhs) ! Right-hand sides

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: soln(ngrdcol,ndim,nrhs) ! Solutions
  real( kind = core_rknd ), intent(out) :: rcond(ngrdcol) ! Optional diagnostic or sentinel [-]
  integer(c_int), intent(out) :: codes(ngrdcol) ! Native status flags [-]

  !--------------------------- Local Variables ---------------------------
  type(err_info_type) :: err_info ! Local status context

  !--------------------------- Begin Code ---------------------------
  allocate(err_info%err_code(ngrdcol))
  err_info%err_code = 0
  err_info%err_code(1) = 50
  err_info%err_header_global = "matrix wrapper test"
  rcond = 1234._core_rknd
  call set_clubb_debug_level_api(debug_level) ! In

  select case(variant)
  case(1)
    if (use_rcond /= 0) then
      call band_solve("matrix test", solve_method, ngrdcol, 2, 2, ndim, implemented /= 0, & ! In
                        lhs, rhs(:,:,1), err_info, & ! InOut
                        soln(:,:,1), & ! Out
                        rcond=rcond) ! Optional Out
    else
      call band_solve("matrix test", solve_method, ngrdcol, 2, 2, ndim, implemented /= 0, & ! In
                        lhs, rhs(:,:,1), err_info, & ! InOut
                        soln(:,:,1)) ! Out
    end if
  case(2)
    if (use_rcond /= 0) then
      call band_solve("matrix test", solve_method, & ! In
                        ngrdcol, 2, 2, ndim, nrhs, implemented /= 0, & ! In
                        lhs, rhs, err_info, & ! InOut
                        soln, & ! Out
                        rcond=rcond) ! Optional Out
    else
      call band_solve("matrix test", solve_method, & ! In
                        ngrdcol, 2, 2, ndim, nrhs, implemented /= 0, & ! In
                        lhs, rhs, err_info, & ! InOut
                        soln) ! Out
    end if
  case(3)
    if (use_rcond /= 0) then
      call tridiag_solve("matrix test", solve_method, ndim, & ! In
                        lhs(:,1,:), rhs(1,:,1), err_info, & ! InOut
                        soln(1,:,1), & ! Out
                        rcond=rcond) ! Optional Out
    else
      call tridiag_solve("matrix test", solve_method, ndim, & ! In
                        lhs(:,1,:), rhs(1,:,1), err_info, & ! InOut
                        soln(1,:,1)) ! Out
    end if
  case(4)
    if (use_rcond /= 0) then
      call tridiag_solve("matrix test", solve_method, ngrdcol, ndim, implemented /= 0, & ! In
                        lhs, rhs(:,:,1), err_info, & ! InOut
                        soln(:,:,1), & ! Out
                        rcond=rcond) ! Optional Out
    else
      call tridiag_solve("matrix test", solve_method, ngrdcol, ndim, implemented /= 0, & ! In
                        lhs, rhs(:,:,1), err_info, & ! InOut
                        soln(:,:,1)) ! Out
    end if
  case(5)
    if (use_rcond /= 0) then
      call tridiag_solve("matrix test", solve_method, & ! In
                        ngrdcol, ndim, nrhs, implemented /= 0, & ! In
                        lhs, rhs, err_info, & ! InOut
                        soln, & ! Out
                        rcond=rcond) ! Optional Out
    else
      call tridiag_solve("matrix test", solve_method, & ! In
                        ngrdcol, ndim, nrhs, implemented /= 0, & ! In
                        lhs, rhs, err_info, & ! InOut
                        soln) ! Out
    end if
  end select
  codes = err_info%err_code
end subroutine f2k_matrix_reference_c

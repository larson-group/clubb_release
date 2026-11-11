! Compare every native tridiagonal/pentadiagonal LU overload with the C++ port.
subroutine f2k_lu_reference_c(ndim, nrhs, ngrdcol, nbands, rhs, lhs, soln, single, legacy) bind(C)
  ! Description: Solve column-varying band systems through the compiled native
  ! generic interfaces. legacy also exercises the single-matrix tridiagonal form.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use tridiag_lu_solvers, only: tridiag_lu_solve
  use penta_lu_solvers, only: penta_lu_solve
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: ndim, nrhs, ngrdcol, nbands ! Problem sizes [-]
  real( kind = core_rknd ), intent(in) :: rhs(ngrdcol,ndim,nrhs) ! Right-hand sides [units vary]

  !--------------------------- Input/Output Variables ---------------------------
  real( kind = core_rknd ), intent(inout) :: lhs(nbands,ngrdcol,ndim) ! Native band storage

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: soln(ngrdcol,ndim,nrhs) ! Multiple-RHS solution
  real( kind = core_rknd ), intent(out) :: single(ngrdcol,ndim) ! Single-RHS solution
  real( kind = core_rknd ), intent(out) :: legacy(ngrdcol,ndim) ! Single-matrix solution

  !--------------------------- Local Variables ---------------------------
  integer :: i ! Column index [-]

  !--------------------------- Begin Code ---------------------------
  if (nbands == 3) then
    call tridiag_lu_solve(ndim, nrhs, ngrdcol, & ! In
                           lhs, & ! InOut
                           rhs, & ! In
                           soln) ! Out
    if (nrhs > 0) then
      call tridiag_lu_solve(ndim, ngrdcol, & ! In
                             lhs, & ! InOut
                             rhs(:,:,1), & ! In
                             single) ! Out
      do i = 1, ngrdcol
        call tridiag_lu_solve(ndim, & ! In
                               lhs(:,i,:), & ! InOut
                               rhs(i,:,1), & ! In
                               legacy(i,:)) ! Out
      end do
    end if
  else
    call penta_lu_solve(ndim, nrhs, ngrdcol, & ! In
                         lhs, & ! InOut
                         rhs, & ! In
                         soln) ! Out
    if (nrhs > 0) then
      call penta_lu_solve(ndim, ngrdcol, & ! In
                           lhs, & ! InOut
                           rhs(:,:,1), & ! In
                           single) ! Out
    end if
  end if
end subroutine f2k_lu_reference_c

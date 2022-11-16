!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module matrix_solver_wrapper

  use clubb_precision, only: &
        core_rknd ! Variable(s)

  implicit none
  
  public :: band_solve, &
            tridiag_solve

  interface band_solve
      module procedure band_solve_single_rhs_multiple_lhs
      module procedure band_solve_multiple_rhs_lhs
  end interface

  interface tridiag_solve
      module procedure tridiag_solve_single_rhs_lhs
      module procedure tridiag_solve_single_rhs_multiple_lhs
      module procedure tridiag_solve_multiple_rhs_lhs
  end interface
            
  contains


  !-------------------------------------------------------------------
  !                     Band Solvers Procedures
  !-------------------------------------------------------------------

  subroutine band_solve_single_rhs_multiple_lhs( &
                solve_name, & !solver_method, & ! Intent(in)
                ngrdcol, nsup, nsub, ndim,    & ! Intent(in)
                lhs, rhs,                     & ! Intent(inout)
                solution, rcond )               ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_band_solve,  & ! Procedure(s)
          lapack_band_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    !integer, intent(in) :: &
    !  solver_method

    integer, intent(in) :: & 
      ngrdcol,  & ! Number of grid columns
      nsup,     & ! Number of superdiagonals
      nsub,     & ! Number of subdiagonals
      ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(inout) ::  & 
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(out) :: &
      solution

    ! ----------------------- Optional Out -----------------------

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) ::  & 
      rcond

    ! ----------------------- Local Variables -----------------------
    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call lapack_band_solvex( "xm_wpxp", nsup, nsub, ndim, 1,  & ! Intent(in) 
                          lhs(:,i,:), rhs(i,:),                   & ! Intent(inout)
                          solution(i,:), rcond(i) )                 ! Intent(out)
      end do

    else

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_band_solve( "xm_wpxp", nsup, nsub, ndim, 1, & ! Intent(in) 
                                lhs(:,i,:), rhs(i,:),           & ! Intent(inout)
                                solution(i,:) )                   ! Intent(out)
      end do

    end if

    return

  end subroutine band_solve_single_rhs_multiple_lhs

  subroutine band_solve_multiple_rhs_lhs( &
                solve_name, & !solver_method,     & ! Intent(in)
                ngrdcol, nsup, nsub, ndim, nrhs,  & ! Intent(in)
                lhs, rhs,                         & ! Intent(inout)
                solution, rcond )                   ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_band_solve,  & ! Procedure(s)
          lapack_band_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    !integer, intent(in) :: &
    !  solver_method

    integer, intent(in) :: & 
      ngrdcol,  & ! Number of grid columns
      nsup,     & ! Number of superdiagonals
      nsub,     & ! Number of subdiagonals
      ndim,     & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs        ! Number of RHS's to back substitute for

    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(inout) ::  & 
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
      solution

    ! ----------------------- Optional Out -----------------------

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) ::  & 
      rcond

    ! ----------------------- Local Variables -----------------------
    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call lapack_band_solvex( "xm_wpxp", nsup, nsub, ndim, nrhs, & ! Intent(in) 
                          lhs(:,i,:), rhs(i,:,:),                   & ! Intent(inout)
                          solution(i,:,:), rcond(i) )                 ! Intent(out)
      end do

    else

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_band_solve( "xm_wpxp", nsup, nsub, ndim, nrhs,  & ! Intent(in) 
                                lhs(:,i,:), rhs(i,:,:),             & ! Intent(inout)
                                solution(i,:,:) )                     ! Intent(out)
      end do

    end if

    return

  end subroutine band_solve_multiple_rhs_lhs



  !-------------------------------------------------------------------
  !                    Tridiag Solver Procedures
  !-------------------------------------------------------------------

  subroutine tridiag_solve_single_rhs_lhs( &
                solve_name, & !solver_method, & ! Intent(in)
                ndim,                         & ! Intent(in)
                lhs, rhs,                     & ! Intent(inout)
                solution, rcond )               ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_tridiag_solve,  & ! Procedure(s)
          lapack_tridiag_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    !integer, intent(in) :: &
    !  solver_method

    integer, intent(in) :: & 
      ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

    real( kind = core_rknd ), dimension(3,ndim), intent(inout) ::  & 
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ndim), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), dimension(ndim), intent(out) :: &
      solution

    ! ----------------------- Optional Out -----------------------

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), optional, intent(out) ::  & 
      rcond

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      call lapack_tridiag_solvex( & 
             solve_name, ndim, 1,                   & ! Intent(in) 
             lhs(1,:), lhs(2,:), lhs(3,:), rhs(:),  & ! Intent(inout)
             solution(:), rcond )                     ! Intent(out)

    else

      ! Perform LU decomp and solve system (LAPACK)
      call lapack_tridiag_solve( & 
             solve_name, ndim, 1,                   & ! Intent(in)
             lhs(1,:), lhs(2,:), lhs(3,:), rhs(:),  & ! Intent(inout)
             solution(:) )                            ! Intent(out)

    end if

    return

  end subroutine tridiag_solve_single_rhs_lhs


  subroutine tridiag_solve_single_rhs_multiple_lhs( &
                solve_name, & !solver_method, & ! Intent(in)
                ngrdcol, ndim,                & ! Intent(in)
                lhs, rhs,                     & ! Intent(inout)
                solution, rcond )               ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_tridiag_solve,  & ! Procedure(s)
          lapack_tridiag_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    !integer, intent(in) :: &
    !  solver_method

    integer, intent(in) :: & 
      ngrdcol,  & ! Number of grid columns
      ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

    real( kind = core_rknd ), dimension(3,ngrdcol,ndim), intent(inout) ::  & 
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(out) :: &
      solution

    ! ----------------------- Optional Out -----------------------

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) ::  & 
      rcond

    ! ----------------------- Local Variables -----------------------
    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call lapack_tridiag_solvex( & 
               solve_name, ndim, 1,                           & ! Intent(in) 
               lhs(1,i,:), lhs(2,i,:), lhs(3,i,:), rhs(i,:),  & ! Intent(inout)
               solution(i,:), rcond(i) )                        ! Intent(out)
      end do

    else

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_tridiag_solve( & 
               solve_name, ndim, 1,                           & ! Intent(in)
               lhs(1,i,:), lhs(2,i,:), lhs(3,i,:), rhs(i,:),  & ! Intent(inout)
               solution(i,:) )                                  ! Intent(out)
      end do

    end if

    return

  end subroutine tridiag_solve_single_rhs_multiple_lhs

  subroutine tridiag_solve_multiple_rhs_lhs( &
                solve_name, & !solver_method, & ! Intent(in)
                ngrdcol, ndim, nrhs,          & ! Intent(in)
                lhs, rhs,                     & ! Intent(inout)
                solution, rcond )               ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_tridiag_solve,  & ! Procedure(s)
          lapack_tridiag_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    !integer, intent(in) :: &
    !  solver_method

    integer, intent(in) :: & 
      ngrdcol,  & ! Number of grid columns
      ndim,     & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs        ! Number of RHS's to back substitute for

    real( kind = core_rknd ), dimension(3,ngrdcol,ndim), intent(inout) ::  & 
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
      solution

    ! ----------------------- Optional Out -----------------------

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) ::  & 
      rcond

    ! ----------------------- Local Variables -----------------------
    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call lapack_tridiag_solvex( & 
               solve_name, ndim, nrhs,                          & ! Intent(in) 
               lhs(1,i,:), lhs(2,i,:), lhs(3,i,:), rhs(i,:,:),  & ! Intent(inout)
               solution(i,:,:), rcond(i) )                        ! Intent(out)
      end do

    else

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_tridiag_solve( & 
               solve_name, ndim, nrhs,                          & ! Intent(in)
               lhs(1,i,:), lhs(2,i,:), lhs(3,i,:), rhs(i,:,:),  & ! Intent(inout)
               solution(i,:,:) )                                  ! Intent(out)
      end do

    end if

    return

  end subroutine tridiag_solve_multiple_rhs_lhs


end module matrix_solver_wrapper
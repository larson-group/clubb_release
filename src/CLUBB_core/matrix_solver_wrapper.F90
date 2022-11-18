!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module matrix_solver_wrapper

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  use error_code, only: &
    clubb_at_least_debug_level, & ! Procedure
    err_code,                   & ! Error indicator
    clubb_no_error,             & ! Constant
    clubb_fatal_error

  use constants_clubb, only: &
    fstderr     ! Constant(s)

  use model_flags, only: &
    lapack ,        & ! Variable(s)
    penta_lu ,      & 
    penta_bicgstab    

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
                solve_name, penta_solve_method, & ! Intent(in)
                ngrdcol, nsup, nsub, ndim,      & ! Intent(in)
                old_solution,                   & ! Intent(in)
                lhs, rhs,                       & ! Intent(inout)
                solution, rcond )                 ! Intent(out)

    use lapack_wrap, only:  & 
      lapack_band_solve,  & ! Procedure(s)
      lapack_band_solvex

    use penta_lu_solvers, only: &
      penta_lu_solve    ! Procedure(s)

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    integer, intent(in) :: &
      penta_solve_method

    integer, intent(in) :: & 
      ngrdcol,  & ! Number of grid columns
      nsup,     & ! Number of superdiagonals
      nsub,     & ! Number of subdiagonals
      ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

    real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(in) ::  & 
      old_solution ! Old solution, used as an initial guess in the bicgstab method

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
    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim) ::  & 
      lhs_copy ! Copy of left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  & 
      rhs_copy ! Copy of right hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim) :: &
      dummy_solution

    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Lapack overwrites lhs and rhs, so we'll give it copies of them.
      lhs_copy = lhs
      rhs_copy = rhs

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      ! Using dummy_solution, since we only want this routine for diagnostics
      do i = 1, ngrdcol
        call lapack_band_solvex( "xm_wpxp", nsup, nsub, ndim, 1,  & ! Intent(in) 
                          lhs(:,i,:), rhs(i,:),                   & ! Intent(inout)
                          dummy_solution(i,:), rcond(i) )           ! Intent(out)
      end do

    end if


    if ( penta_solve_method == lapack ) then

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_band_solve( "xm_wpxp", nsup, nsub, ndim, 1, & ! Intent(in) 
                                lhs(:,i,:), rhs(i,:),           & ! Intent(inout)
                                solution(i,:) )                   ! Intent(out)
      end do

    else if ( penta_solve_method == penta_lu ) then 

      ! Solve the system with a penta-diagonal specific LU decomp
      call penta_lu_solve( ndim, ngrdcol,         & ! Intent(in)
                           lhs(:,:,:), rhs(:,:),  & ! Intent(in)
                           solution(:,:) )          ! Intent(out)

    else

      ! The solve method should match one of the above
      if ( clubb_at_least_debug_level( 0 ) ) then
        write(fstderr,*) "Error in band_solve_single_rhs_multiple_lhs: "
        write(fstderr,*) "  no case for penta_solve_method = ", penta_solve_method
        err_code = clubb_fatal_error
      end if

    end if

    return

  end subroutine band_solve_single_rhs_multiple_lhs

  subroutine band_solve_multiple_rhs_lhs( &
                solve_name, penta_solve_method,   & ! Intent(in)
                ngrdcol, nsup, nsub, ndim, nrhs,  & ! Intent(in)
                old_solution,                     & ! Intent(in)
                lhs, rhs,                         & ! Intent(inout)
                solution, rcond )                   ! Intent(out)

    use lapack_wrap, only:  & 
      lapack_band_solve,  & ! Procedure(s)
      lapack_band_solvex

    use penta_lu_solvers, only: &
      penta_lu_solve    ! Procedure(s)

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    integer, intent(in) :: &
      penta_solve_method

    integer, intent(in) :: & 
      ngrdcol,  & ! Number of grid columns
      nsup,     & ! Number of superdiagonals
      nsub,     & ! Number of subdiagonals
      ndim,     & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs        ! Number of RHS's to back substitute for

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(in) ::  & 
      old_solution ! Old solution, used as an initial guess in the bicgstab method

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
    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim) ::  & 
      lhs_copy ! Copy of left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs) ::  & 
      rhs_copy ! Copy of right hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs) :: &
      dummy_solution

    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Lapack overwrites lhs and rhs, so we'll give it copies of them.
      lhs_copy = lhs
      rhs_copy = rhs

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      ! Using dummy_solution, since we only want this routine for diagnostics
      do i = 1, ngrdcol
        call lapack_band_solvex( "xm_wpxp", nsup, nsub, ndim, nrhs, & ! Intent(in) 
                          lhs(:,i,:), rhs(i,:,:),                   & ! Intent(inout)
                          dummy_solution(i,:,:), rcond(i) )           ! Intent(out)
      end do

    end if


    if ( penta_solve_method == lapack ) then

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_band_solve( "xm_wpxp", nsup, nsub, ndim, nrhs,  & ! Intent(in) 
                                lhs(:,i,:), rhs(i,:,:),             & ! Intent(inout)
                                solution(i,:,:) )                     ! Intent(out)
      end do

    else if ( penta_solve_method == penta_lu ) then 

      ! Solve the system with a penta-diagonal specific LU decomp
      call penta_lu_solve( ndim, nrhs, ngrdcol,     & ! Intent(in)
                           lhs(:,:,:), rhs(:,:,:),  & ! Intent(in)
                           solution(:,:,:) )          ! Intent(out)

    else

      ! The solve method should match one of the above
      if ( clubb_at_least_debug_level( 0 ) ) then
        write(fstderr,*) "Error in band_solve_multiple_rhs_lhs: "
        write(fstderr,*) "  no case for penta_solve_method = ", penta_solve_method
        err_code = clubb_fatal_error
      end if

    end if

    return

  end subroutine band_solve_multiple_rhs_lhs



  !-------------------------------------------------------------------
  !                    Tridiag Solver Procedures
  !-------------------------------------------------------------------

  subroutine tridiag_solve_single_rhs_lhs( &
                solve_name, tridiag_solve_method, & ! Intent(in)
                ndim,                             & ! Intent(in)
                lhs, rhs,                         & ! Intent(inout)
                solution, rcond )                   ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_tridiag_solve,  & ! Procedure(s)
          lapack_tridiag_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    integer, intent(in) :: &
      tridiag_solve_method

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

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(3,ndim) ::  & 
      lhs_copy ! Copy of left hand side

    real( kind = core_rknd ), dimension(ndim) ::  & 
      rhs_copy ! Copy of right hand side

    real( kind = core_rknd ), dimension(ndim) :: &
      dummy_solution

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Lapack overwrites lhs and rhs, so we'll give it copies of them.
      lhs_copy = lhs
      rhs_copy = rhs

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      call lapack_tridiag_solvex( & 
             solve_name, ndim, 1,           & ! Intent(in) 
             lhs_copy(1,:), lhs_copy(2,:),  & ! Intent(in) 
             lhs_copy(3,:), rhs_copy(:),    & ! Intent(inout)
             dummy_solution(:), rcond )       ! Intent(out)

    end if


    if ( tridiag_solve_method == lapack ) then

      ! Perform LU decomp and solve system (LAPACK)
      call lapack_tridiag_solve( & 
             solve_name, ndim, 1,                   & ! Intent(in)
             lhs(1,:), lhs(2,:), lhs(3,:), rhs(:),  & ! Intent(inout)
             solution(:) )                            ! Intent(out)

    else

      ! The solve method should match one of the above
      if ( clubb_at_least_debug_level( 0 ) ) then
        write(fstderr,*) "Error in tridiag_solve_single_rhs_lhs: "
        write(fstderr,*) "  no case for tridiag_solve_method = ", tridiag_solve_method
        err_code = clubb_fatal_error
      end if

    end if

    return

  end subroutine tridiag_solve_single_rhs_lhs


  subroutine tridiag_solve_single_rhs_multiple_lhs( &
                solve_name, tridiag_solve_method, & ! Intent(in)
                ngrdcol, ndim,                    & ! Intent(in)
                lhs, rhs,                         & ! Intent(inout)
                solution, rcond )                   ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_tridiag_solve,  & ! Procedure(s)
          lapack_tridiag_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    integer, intent(in) :: &
      tridiag_solve_method

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
    real( kind = core_rknd ), dimension(3,ngrdcol,ndim) ::  & 
      lhs_copy ! Copy of left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  & 
      rhs_copy ! Copy of right hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim) :: &
      dummy_solution

    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Lapack overwrites lhs and rhs, so we'll give it copies of them.
      lhs_copy = lhs
      rhs_copy = rhs

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call lapack_tridiag_solvex( & 
               solve_name, ndim, 1,               & ! Intent(in) 
               lhs_copy(1,i,:), lhs_copy(2,i,:),  & ! Intent(in) 
               lhs_copy(3,i,:), rhs_copy(i,:),    & ! Intent(inout)
               dummy_solution(i,:), rcond(i) )      ! Intent(out)
      end do

    end if


    if ( tridiag_solve_method == lapack ) then

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_tridiag_solve( & 
               solve_name, ndim, 1,                           & ! Intent(in)
               lhs(1,i,:), lhs(2,i,:), lhs(3,i,:), rhs(i,:),  & ! Intent(inout)
               solution(i,:) )                                  ! Intent(out)
      end do

    else

      ! The solve method should match one of the above
      if ( clubb_at_least_debug_level( 0 ) ) then
        write(fstderr,*) "Error in tridiag_solve_single_rhs_multiple_lhs: "
        write(fstderr,*) "  no case for tridiag_solve_method = ", tridiag_solve_method
        err_code = clubb_fatal_error
      end if

    end if

    return

  end subroutine tridiag_solve_single_rhs_multiple_lhs

  subroutine tridiag_solve_multiple_rhs_lhs( &
                solve_name, tridiag_solve_method, & ! Intent(in)
                ngrdcol, ndim, nrhs,              & ! Intent(in)
                lhs, rhs,                         & ! Intent(inout)
                solution, rcond )                   ! Intent(out)

    use lapack_wrap, only:  & 
          lapack_tridiag_solve,  & ! Procedure(s)
          lapack_tridiag_solvex

    implicit none

    ! ----------------------- Input Variables -----------------------
    character(len=*), intent(in) :: &
      solve_name

    integer, intent(in) :: &
      tridiag_solve_method

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
    real( kind = core_rknd ), dimension(3,ngrdcol,ndim) ::  & 
      lhs_copy ! Copy of left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs) ::  & 
      rhs_copy ! Copy of right hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs) :: &
      dummy_solution

    integer :: i

    ! ----------------------- Begin Code -----------------------

    if ( present(rcond) ) then

      ! Lapack overwrites lhs and rhs, so we'll give it copies of them.
      lhs_copy = lhs
      rhs_copy = rhs

      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call lapack_tridiag_solvex( & 
               solve_name, ndim, nrhs,            & ! Intent(in) 
               rhs_copy(1,i,:), rhs_copy(2,i,:),  & ! Intent(in) 
               rhs_copy(3,i,:), rhs_copy(i,:,:),  & ! Intent(inout)
               dummy_solution(i,:,:), rcond(i) )    ! Intent(out)
      end do

    end if


    if ( tridiag_solve_method == lapack ) then

      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call lapack_tridiag_solve( & 
               solve_name, ndim, nrhs,                          & ! Intent(in)
               lhs(1,i,:), lhs(2,i,:), lhs(3,i,:), rhs(i,:,:),  & ! Intent(inout)
               solution(i,:,:) )                                  ! Intent(out)
      end do

    else

      ! The solve method should match one of the above
      if ( clubb_at_least_debug_level( 0 ) ) then
        write(fstderr,*) "Error in tridiag_solve_multiple_rhs_lhs: "
        write(fstderr,*) "  no case for tridiag_solve_method = ", tridiag_solve_method
        err_code = clubb_fatal_error
      end if

    end if

    return

  end subroutine tridiag_solve_multiple_rhs_lhs


end module matrix_solver_wrapper
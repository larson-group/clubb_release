!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module lapack_wrap

! Description:
!   Wrappers for the band diagonal and tridiagonal direct matrix
!   solvers contained in the LAPACK library.

! References:
!   LAPACK--Linear Algebra PACKage
!   URL: <http://www.netlib.org/lapack/>
!-----------------------------------------------------------------------
  use constants_clubb, only: &
    fstderr ! Variable(s)

  use clubb_precision, only: &
    core_rknd ! Variable

  implicit none

  ! Simple routines
  public :: lapack_tridiag_solve, &
            lapack_band_solve

  ! Expert routines
  public :: lapack_tridiag_solvex, &
            lapack_band_solvex

  private ! Set Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine lapack_tridiag_solvex( solve_type, ndim, nrhs, ngrdcol, &
                                    lhs, rhs, err_info, &
                                    soln, rcond )

! Description:
!   Solves a tridiagonal system of equations (expert routine).

! References:
!   <http://www.netlib.org/lapack/single/sgtsvx.f>
!   <http://www.netlib.org/lapack/double/dgtsvx.f>

! Notes:
!   More expensive than the simple routine, but tridiagonal
!   decomposition is still relatively cheap.
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error                 ! Constants

    use lapack_interfaces, only: &
        lapack_gtsvx, &      ! Procedure
        lapack_isnan

    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    intrinsic :: kind

    ! ----------------------- Input variables -----------------------
    character(len=*), intent(in) :: &
      solve_type ! Used to write a message if this fails

    integer, intent(in) :: &
      ndim,     & ! N-dimension of matrix
      ngrdcol,  & ! Number of grid columns
      nrhs        ! # of right hand sides to back subst. after LU-decomp.

    ! ----------------------- Input/Output variables -----------------------
    real( kind = core_rknd ), intent(inout), dimension(3,ngrdcol,ndim) :: &
      lhs ! Tridiagonal LHS

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,ndim,nrhs) :: &
      rhs ! RHS input

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! The estimate of the reciprocal of the condition number on the LHS matrix.
    ! If rcond is < machine precision the matrix is singular to working
    ! precision, and info == ndim+1.  If rcond == 0, then the LHS matrix
    ! is singular.  This condition is indicated by a return code of info > 0
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: rcond

    ! ----------------------- Output variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim,nrhs) :: &
      soln ! solution

    ! ----------------------- Local Variables -----------------------
    ! These contain the decomposition of the matrix
    real( kind = core_rknd ), dimension(ndim-1) :: dlf, duf
    real( kind = core_rknd ), dimension(ndim)   :: df
    real( kind = core_rknd ), dimension(ndim-2) :: du2

    integer, dimension(ndim) :: &
      ipivot  ! Index of pivots done during decomposition

    integer, dimension(ndim) :: &
      iwork   ! `scrap' array


    real( kind = core_rknd ), dimension(ngrdcol,nrhs) :: &
      ferr,  & ! Forward error estimate
      berr     ! Backward error estimate

    real( kind = core_rknd ), dimension(3*ndim) :: &
      work  ! `Scrap' array

    integer :: info ! Diagnostic output

    integer :: i, n  ! Array index

!-----------------------------------------------------------------------
!     *** The LAPACK Routine ***
!     SUBROUTINE SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,
!    $                   DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,
!    $                   WORK, IWORK, INFO )
!-----------------------------------------------------------------------

    ! Lapack tridiagonal matrix solver, expert version, sgtsvx for single
    ! or dgtsvx for double precision
    do i = 1, ngrdcol
      call lapack_gtsvx( "Not Factored", "No Transpose lhs", ndim, nrhs,  &
                         lhs(3,i,2:ndim), lhs(2,i,:), lhs(1,i,1:ndim-1),  &
                         dlf, df, duf, du2, ipivot,  &
                         rhs(i,:,:), ndim, soln(i,:,:), ndim, rcond(i), &
                         ferr(i,:), berr(i,:), work, iwork, info )
    end do

    ! Print diagnostics for when ferr is large
    if ( clubb_at_least_debug_level_api( 2 ) .and. any( ferr > 1.e-3_core_rknd ) ) then

      write(fstderr,*) "Warning, large error est. for: " // trim( solve_type )

      do n = 1, nrhs
        do i = 1, ngrdcol
          write(fstderr,*) "grdcol #", i, "rhs # ", i, "tridiag forward error est. =", ferr(i,n)
          write(fstderr,*) "grdcol #", i, "rhs # ", i, "tridiag backward error est. =", berr(i,n)
        end do
      end do

      write(fstderr,'(2(a20,e15.6))') "rcond est. = ", rcond, &
        "machine epsilon = ", epsilon( lhs(1,1,1) )
    end if

    select case( info )
    case( :-1 )
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "lapack_tridiag_solvex"
      write(fstderr, *) trim( solve_type )// &
        "illegal value in argument", -info
      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error

    case( 0 )
      ! Success
      do i = 1, ngrdcol
        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( lapack_isnan( ndim, nrhs, soln(i,:,:) ) ) then
            write(fstderr, *) err_info%err_header(i)
            write(fstderr, *) "lapack_tridiag_solvex"
            write(fstderr, *) trim( solve_type )// &
              "NaNs in solution", info
            ! Error in grid column i -> set ith entry to clubb_fatal_error
            err_info%err_code(i) = clubb_fatal_error
          end if
        end if
      end do

    case( 1: )
      if ( info == ndim+1 ) then
        write(fstderr, *) "lapack_tridiag_solvex"
        write(fstderr,*) trim( solve_type) // &
          " Warning: matrix is singular to working precision."
        write(fstderr,'(a,e12.5)')  &
          "Estimate of the reciprocal of the condition number: ", rcond
      else
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) solve_type// &
          " singular matrix."
        ! General error -> set all entries to clubb_fatal_error
        err_info%err_code = clubb_fatal_error
      end if

    end select

    return
  end subroutine lapack_tridiag_solvex

!-----------------------------------------------------------------------
  subroutine lapack_tridiag_solve( solve_type, ndim, nrhs, ngrdcol, &
                                   lhs, rhs, err_info, &
                                   soln )

! Description:
!   Solves a tridiagonal system of equations (simple routine)

! References:
!   <http://www.netlib.org/lapack/single/sgtsv.f>
!   <http://www.netlib.org/lapack/double/dgtsv.f>
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

#ifdef E3SM
#ifndef NDEBUG
#if defined(ARCH_MIC_KNL) && defined(CPRINTEL)
    use, intrinsic :: ieee_exceptions
#endif
#endif
#endif /*E3SM*/

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error                 ! Constants

    use lapack_interfaces, only: &
        lapack_gtsv, &       ! Procedure
        lapack_isnan

    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    intrinsic :: kind

    ! ----------------------- Input variables -----------------------
    character(len=*), intent(in) :: &
      solve_type ! Used to write a message if this fails

    integer, intent(in) :: &
      ndim,     & ! N-dimension of matrix
      ngrdcol,  & ! Number of grid columns
      nrhs        ! # of right hand sides to back subst. after LU-decomp.

    ! ----------------------- Input/Output Variables -----------------------
    real( kind = core_rknd ), intent(inout), dimension(3,ngrdcol,ndim) :: &
      lhs ! Tridiagonal LHS input

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,ndim,nrhs) :: &
      rhs ! RHS input

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! ----------------------- Output variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim,nrhs) :: &
      soln ! solution

    ! ----------------------- Local Variables -----------------------
    integer :: &
      info, & ! Diagnostic output
      i       ! Loop var

!-----------------------------------------------------------------------
!       *** The LAPACK Routine ***
!       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!-----------------------------------------------------------------------

#ifdef E3SM
#ifndef NDEBUG
#if defined(ARCH_MIC_KNL) && defined(CPRINTEL)
    ! when floating-point exceptions are turned on, this call was failing with
    ! a div-by-zero on KNL with Intel/MKL. Solution was to turn off exceptions
    ! only here at this call (and only for machine with ARCH_MIC_KNL defined)
    ! (github 1183)
    call ieee_set_halting_mode(IEEE_DIVIDE_BY_ZERO, .false.) ! Turn off stopping on div-by-zero only
#endif
#endif
#endif /*E3SM*/

    ! Interface for Lapack tridiagonal matrix solver, sgtsv for single
    ! or dgtsv for double precision
    do i = 1, ngrdcol
      call lapack_gtsv( ndim, nrhs, lhs(3,i,2:ndim), lhs(2,i,:), lhs(1,i,1:ndim-1),  &
                        rhs(i,:,:), ndim, info )
    end do

#ifdef E3SM
#ifndef NDEBUG
#if defined(ARCH_MIC_KNL) && defined(CPRINTEL)
    ! Turn back on stopping on div-by-zero only
    call ieee_set_halting_mode(IEEE_DIVIDE_BY_ZERO, .true.)
#endif
#endif
#endif /*E3SM*/

    select case( info )
    case( :-1 )
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "lapack_tridiag_solve"
      write(fstderr,*) trim( solve_type )// &
        " illegal value in argument", -info
      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error

      soln = -999._core_rknd

    case( 0 )
      ! Success
      do i = 1, ngrdcol
        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( lapack_isnan( ndim, nrhs, rhs(i,:,:) ) ) then
            write(fstderr, *) err_info%err_header(i)
            write(fstderr, *) "lapack_tridiag_solve"
            write(fstderr, *) trim( solve_type )// &
              "NaNs in solution", info
            ! Error in grid column i -> set ith entry to clubb_fatal_error
            err_info%err_code(i) = clubb_fatal_error
          end if
        end if
      end do

      soln = rhs

    case( 1: )
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "lapack_tridiag_solve"
      write(fstderr,*) trim( solve_type )//" singular matrix."
      ! General error -> set all entries to clubb_fatal_error
      err_info%err_code = clubb_fatal_error

      soln = -999._core_rknd

    end select

    return
  end subroutine lapack_tridiag_solve

!-----------------------------------------------------------------------
  subroutine lapack_band_solvex( solve_type, nsup, nsub, &
                                 ndim, nrhs, ngrdcol, &
                                 lhs, rhs, err_info, &
                                 soln, rcond )

! Description:
!   Restructure and then solve a band diagonal system, with
!   diagnostic output

! References:
!   <http://www.netlib.org/lapack/single/sgbsvx.f>
!   <http://www.netlib.org/lapack/double/dgbsvx.f>

! Notes:
!   I found that due to the use of sgbcon/dgbcon it is much
!   more expensive to use this on most systems than the simple
!   driver. Use this version only if you don't case about compute time.
!   Also note that this version equilibrates the lhs and does an iterative
!   refinement of the solutions, which results in a slightly different answer
!   than the simple driver does. -dschanen 24 Sep 2008
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error                 ! Constants

    use lapack_interfaces, only: &
        lapack_gbsvx, &      ! Procedures
        lapack_isnan

    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    ! ------------------------------ Input Variables ------------------------------
    character(len=*), intent(in) :: solve_type

    integer, intent(in) :: &
      nsup,    & ! Number of superdiagonals
      nsub,    & ! Number of subdiagonals
      ngrdcol, & ! Number of grid columns
      ndim,    & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs       ! Number of RHS's to solve for

    ! ------------------------------ InOut Variables ------------------------------
    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(inout) :: &
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) :: &
      rhs ! Right hand side(s)

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! ------------------------------ Output Variables ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
      soln

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: &
      rcond

    ! ------------------------------ Local Variables ------------------------------

    ! Workspaces
    real( kind = core_rknd ), dimension(3*ndim) :: work
    integer, dimension(ndim) :: iwork

    real( kind = core_rknd ), dimension(2*nsub+nsup+1,ndim) :: &
      lulhs ! LU Decomposition of the LHS

    integer, dimension(ndim) :: &
      ipivot

    real( kind = core_rknd ), dimension(ngrdcol,nrhs) :: &
      ferr, berr ! Forward and backward error estimate

    real( kind = core_rknd ), dimension(ndim) :: &
      rscale, cscale ! Row and column scale factors for the LHS

    integer :: &
      info,   & ! If this doesn't come back as 0, something went wrong
      offset, & ! Loop iterator
      imain,  & ! Main diagonal of the matrix
      i,n       ! Loop iterator

    character :: &
      equed ! Row equilibration status


!-----------------------------------------------------------------------
!       Reorder Matrix to use LAPACK band matrix format (5x6)

!       Shift example:

!       [    *        *     lhs(1,1) lhs(1,2) lhs(1,3) lhs(1,4) ] (2)=>
!       [    *     lhs(2,1) lhs(2,2) lhs(2,3) lhs(2,4) lhs(2,5) ] (1)=>
!       [ lhs(3,1) lhs(3,2) lhs(3,3) lhs(3,4) lhs(3,5) lhs(3,6) ]
! <=(1) [ lhs(4,2) lhs(4,3) lhs(4,4) lhs(4,5) lhs(4,6)    *     ]
! <=(2) [ lhs(5,3) lhs(5,4) lhs(5,5) lhs(5,6)    *        *     ]

!       The '*' indicates unreferenced elements.
!       For additional bands above and below the main diagonal, the
!       shifts to the left or right increases by the distance from the
!       main diagonal of the matrix.
!-----------------------------------------------------------------------

    imain = nsup + 1

    ! For the offset, (+) is left, and (-) is right

    ! Sub diagonals
    do i = 1, ngrdcol
      do offset = 1, nsub, 1
        lhs(imain+offset,i,1:ndim) = eoshift( lhs(imain+offset,i,1:ndim), offset )
      end do
    end do

    ! Super diagonals
    do i = 1, ngrdcol
      do offset = 1, nsup, 1
        lhs(imain-offset,i,1:ndim) = eoshift( lhs(imain-offset,i,1:ndim), -offset )
      end do
    end do

!-----------------------------------------------------------------------
!     *** The LAPACK Routine ***
!     SUBROUTINE SGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,
!    $                   LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,
!    $                   RCOND, FERR, BERR, WORK, IWORK, INFO )
!-----------------------------------------------------------------------

    ! Lapack general band solver, expert version, sgbsvx for single
    ! or dgbsvx for double precision
    do i = 1, ngrdcol
      call lapack_gbsvx( 'Equilibrate lhs', 'No Transpose lhs', &
                         ndim, nsub, nsup, nrhs, &
                         lhs(:,i,:), nsup+nsub+1, lulhs, 2*nsub+nsup+1, &
                         ipivot, equed, rscale, cscale, &
                         rhs(i,:,:), ndim, soln(i,:,:), ndim, &
                         rcond(i), ferr(i,:), berr(i,:), work, iwork, info )
    end do


! %% debug
!       select case ( equed )
!       case ('N')
!         print *, "No equilib. was required for lhs."
!       case ('R')
!         print *, "Row equilib. was done on lhs."
!       case ('C')
!         print *, "Column equilib. was done on lhs."
!       case ('B')
!         print *, "Row and column equilib. was done on lhs."
!       end select

!       write(*,'(a,e12.5)') "Row scale : ", rscale
!       write(*,'(a,e12.5)') "Column scale: ", cscale
!       write(*,'(a,e12.5)') "Estimate of the reciprocal of the "//
!                            "condition number: ", rcond
!       write(*,'(a,e12.5)') "Forward Error Estimate: ", ferr
!       write(*,'(a,e12.5)') "Backward Error Estimate: ", berr
! %% end debug

    ! Diagnostic information
    if ( clubb_at_least_debug_level_api( 2 ) .and. any( ferr > 1.e-3_core_rknd ) ) then

      write(fstderr,*) "Warning, large error est. for: " // trim( solve_type )

      do n = 1, nrhs
        do i = 1, ngrdcol
          write(fstderr,*) "grdcol #", i, "rhs # ", n, "band_solvex forward error est. =", &
                           ferr(i,n)
          write(fstderr,*) "grdcol #", i, "rhs # ", n, "band_solvex backward error est. =", &
                           berr(i,n)
        end do
      end do

      write(fstderr,'(2(a20,e15.6))') "rcond est. = ", rcond, &
        "machine epsilon = ", epsilon( lhs(1,1,1) )
    end if

    select case( info )

    case( :-1 )
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "in band_solvex for ", trim( solve_type ), &
            ": illegal value for argument", -info
        ! General error -> set all entries to clubb_fatal_error
        err_info%err_code = clubb_fatal_error

    case( 0 )
      ! Success!
      do i = 1, ngrdcol
        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( lapack_isnan( ndim, nrhs, soln(i,:,:) ) ) then
            write(fstderr, *) err_info%err_header(i)
            write(fstderr, *) "lapack_band_solvex"
            write(fstderr, *) trim( solve_type )// &
              "NaNs in solution", info
            ! Error in grid column i -> set ith entry to clubb_fatal_error
            err_info%err_code(i) = clubb_fatal_error
          end if
        end if
      end do

    case( 1: )
      if ( info == ndim+1 ) then

        write(fstderr,*) "in band_solvex for", trim( solve_type ), &
          " Warning: matrix singular to working precision."
        write(fstderr,'(a,e12.5)') &
          "Estimate of the reciprocal of the"// &
          " condition number: ", rcond
      else
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "in band_solvex for", trim( solve_type ), &
          ": singular matrix, solution not computed"
        ! General error -> set all entries to clubb_fatal_error
        err_info%err_code = clubb_fatal_error
      end if

    end select

    return
  end subroutine lapack_band_solvex

!-----------------------------------------------------------------------
  subroutine lapack_band_solve( solve_type, nsup, nsub, &
                                ndim, nrhs, ngrdcol, &
                                lhs, rhs, err_info, &
                                soln )
! Description:
!   Restructure and then solve a band diagonal system

! References:
!   <http://www.netlib.org/lapack/single/sgbsv.f>
!   <http://www.netlib.org/lapack/double/dgbsv.f>
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error                 ! Constants

    use lapack_interfaces, only: &
        lapack_gbsv, &       ! Procedures
        lapack_isnan

    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    ! ------------------------------ Input Variables ------------------------------
    character(len=*), intent(in) :: solve_type

    integer, intent(in) :: &
      nsup,    & ! Number of superdiagonals
      nsub,    & ! Number of subdiagonals
      ngrdcol, & ! Number of grid columns
      ndim,    & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs       ! Number of RHS's to solve for

    ! Note: matrix lhs is intent(in), not intent(inout)
    ! as in the subroutine band_solvex( )
    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(in) :: &
      lhs ! Left hand side

    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) :: &
      rhs ! Right hand side(s)

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! ------------------------------ Output Variables ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
      soln

    ! ------------------------------ Local Variables ------------------------------

    ! Workspaces
    real( kind = core_rknd ), dimension(2*nsub+nsup+1,ndim,ngrdcol) :: &
      lulhs ! LU Decomposition of the LHS

    integer, dimension(ndim) :: &
      ipivot

    integer :: &
      info,   & ! If this doesn't come back as 0, something went wrong
      imain  ! Main diagonal of the matrix

    integer :: i, j, d

    !-----------------------------------------------------------------------
    !       Reorder LU Matrix to use LAPACK band matrix format
    ! 
    !       Shift example for lulhs matrix given a 5x5 lhs matrix
    !
    !
    !  lulhs =
    !                            Columns
    !         1  2       3          4          5         6           7
    ! Rows
    ! 1     [ 0  0       0          0      lhs(3,1)   lhs(4,2)   lhs(5,3) ]
    ! 2     [ 0  0       0      lhs(2,1)   lhs(3,2)   lhs(4,3)   lhs(5,4) ]
    ! 3     [ 0  0   lhs(1,1)   lhs(2,2)   lhs(3,3)   lhs(4,4)   lhs(5,5) ]
    ! 4     [ 0  0   lhs(1,2)   lhs(2,3)   lhs(3,4)   lhs(4,5)       0    ]
    ! 5     [ 0  0   lhs(1,3)   lhs(2,4)   lhs(3,5)       0          0    ]
    !
    !         all       lhs        lhs        lhs        lhs        lhs
    !        set to   shifted     shifted      no      shifted    shifted
    !          0       down 2      down 1    shift      up 1        up 2
    !
    !   The first nsup columns of lulhs are always set to 0;
    !   the rest of the columns are set to shifted
    !   columns of lhs. This can be thought of as taking lhs, never touching the middle column, but
    !   shifting the columns that are n columns to the left of the middle down by n rows, and then
    !   shifting the columns that are n columns to the right of the middle up by n rows, finally
    !   adding nsup columns of zeros onto the left of the array. This results in lulhs.
    !-----------------------------------------------------------------------

    ! Reorder lulhs, omitting the additional 2*nsub bands
    ! that are used for the LU decomposition of the matrix.

    imain = nsub + nsup + 1


    ! The first nsup rows of lulhs will contain 0s that are end-shifted lhs values. This needs
    ! to be handled differently so the algorithm to access lhs will not try to use out of bound
    ! values.
    !             ...   nsup     nsup+1    ...       imain   ...
    !                        \ /                 \ /
    !                always   |  begins with nsup | all lhs values
    !                   0       0s, and decreases
    !                           by one 0 each row
    !
    !              ...  nsup     nsup+1    ...       imain   ...
    ! lulhs(:,1) =  0     0       0         0         lhs    lhs
    ! lulhs(:,2) =  0     0       0        lhs        lhs    lhs
    !
    ! Since the first nsup rows are the first rows in lulhs, we're going to access them first to
    ! avoid out of order memory accesses.
    do i = 1, ngrdcol
      do d = 1, nsup

        ! Add 0s to first nsup columns, and decreasing number of end-shift affected columns
        do j = 1, imain-d
          lulhs(j,d,i) = 0.0_core_rknd
        end do

        ! Copy lhs values into appropriate lulhs spots
        do j = imain-d+1, imain+nsub
          lulhs(j,d,i) = lhs(j-nsub,i,d+j-imain)
        end do

      end do
    end do

    ! After the first nsup rows are dealt with, the offset lhs values can be copied into lulhs
    ! until the last nsup rows are reached. This is because the last nsup rows also contain
    ! end-shifted values, set to 0 in the next loop.
    !
    !                      ...  nsup     nsup+1    ...
    !                                \ /
    !                       always    |    all lhs values
    !
    !                      ...  nsup     nsup+1    ...
    ! lulhs(:,nsup+1)    =  0     0       lhs      lhs
    ! lulhs(:,ndim-nsub) =  0     0       lhs      lhs
    !
    ! For all values not affected by end-shifting
    do i = 1, ngrdcol
      do d = nsup+1, ndim-nsub

        ! Set first nsup columns to 0
        do j = 1, nsub
          lulhs(j,d,i) = 0.0_core_rknd
        end do

        ! Copy lhs values into appropriate lulhs spots
        do j = imain-nsub, imain+nsub
          lulhs(j,d,i) = lhs(j-nsub,i, d+j-imain)
        end do

      end do
    end do


    ! The last nsup rows of lulhs will contain 0s that are end-shifted lhs values. This needs
    ! to be handled differently so the algorithm to access lhs will not try to use out of bound
    ! values.
    !
    !
    !                        ...  nsup     nsup+1    ...      imain+1    ...
    ! lulhs(:,ndim-nsub+1) =  0     0       lhs      lhs        lhs       0
    ! lulhs(:,ndim)        =  0     0       lhs      lhs         0        0
    !
    !                                   |                   |   starts with one 0, then
    !                          always   |   all lhs values  |  then increases to nsup 0s
    !                                   |                   |       towards ndim
    !                                  / \                 / \
    !                        ...  nsup     nsup+1    ...          ndim-nsub+1 ... ndim
    !
    ! Finish the lulhs setup by accessing the last values last, keeping memory access ordered
    do i = 1, ngrdcol
      do d = ndim-nsub+1, ndim

        ! Set first nsup columns to 0
        do j = 1, nsub
          lulhs(j,d,i) = 0.0_core_rknd
        end do

        ! Copy lhs values into appropriate lulhs spots
        do j = imain-nsup, imain-(d-ndim)
          lulhs(j,d,i) = lhs(j-nsub,i, d+j-imain)
        end do

        ! Set increasing number of end-shift affected columns to 0
        do j = imain-(d-ndim)+1, imain+nsub
          lulhs(j,d,i) = 0.0_core_rknd
        end do

      end do
    end do

!-----------------------------------------------------------------------
!       *** LAPACK routine ***
!       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!-----------------------------------------------------------------------

    ! Lapack general band solver, sgbsv for single
    ! or dgbsv for double precision
    do i = 1, ngrdcol
      call lapack_gbsv( ndim, nsub, nsup, nrhs, lulhs(:,:,i), nsub*2+nsup+1, &
                        ipivot, rhs(i,:,:), ndim, info )
    end do


    select case( info )

    case( :-1 )
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "in band_solve for ", trim( solve_type ), &
            ": illegal value for argument", -info
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
    case( 0 )
          ! Success!
          if ( clubb_at_least_debug_level_api( 0 ) ) then
            do i = 1, ngrdcol
              if ( lapack_isnan( ndim, nrhs, rhs(i,:,:) ) ) then
                write(fstderr, *) err_info%err_header(i)
                write(fstderr, *) "lapack_band_solve"
                write(fstderr, *) trim( solve_type )// &
                "NaNs in solution", info
                ! Error in grid column i -> set ith entry to clubb_fatal_error
                err_info%err_code(i) = clubb_fatal_error
              end if
            end do
          end if

          soln = rhs

    case( 1: )
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "in band_solve for ", trim( solve_type ), &
                       ": singular matrix, solution not computed"
        ! General error -> set all entries to clubb_fatal_error
        err_info%err_code = clubb_fatal_error
    end select

    return
  end subroutine lapack_band_solve

end module lapack_wrap

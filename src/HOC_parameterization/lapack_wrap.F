!-----------------------------------------------------------------------
! $Id: lapack_wrap.F,v 1.4 2008-03-31 18:12:30 faschinj Exp $
        module lapack_wrap

!       Description:
!       Wrappers for the band diagonal and tridiagonal direct matrix
!       solvers contained in the LAPACK library.

!       References:
!       LAPACK--Linear Algebra PACKage version 3.0
!       URL: <http://netlib2.cs.utk.edu/lapack/>
!-----------------------------------------------------------------------
        use constants, only: fstderr
        use error_code
        
        implicit none

        ! Simple routines
        public :: tridag_solve, band_solve

        ! Expert routines
        public :: tridag_solvex, band_solvex

        private

        contains

!-----------------------------------------------------------------------
        subroutine tridag_solvex
     .             ( solve_type, ndim, nrhs, supd, diag, subd, 
     .               rhs, solution, rcond, code )

!       Description:
!       Solves a tridiagonal system of equations (expert routine).

!       References:
!       <http://www.netlib.org/lapack/single/sgtsvx.f>
!       <http://www.netlib.org/lapack/double/dgtsvx.f>

!       Notes:
!       More expensive than the simple routine, but tridiagonal
!       decomposition is still relatively cheap.
!-----------------------------------------------------------------------

        implicit none

        ! External
        external :: 
     .  sgtsvx, ! Single-prec. General Tridiagonal Solver eXpert
     .  dgtsvx  ! Double-prec. General Tridiagonal Solver eXpert

        intrinsic :: kind

        ! Input variables
        character(len=*), intent(in) :: 
     .  solve_type ! Used to write a message if this fails

        integer :: 
     .  ndim, ! N-dimension of matrix
     .  nrhs  ! # of right hand sides to back subst. after LU-decomp.

        ! Input/Output variables
        real, intent(inout), dimension(ndim) :: 
     .  diag,      ! Main diagonal
     .  subd, supd ! Sub and super diagonal

        real, intent(inout), dimension(ndim,nrhs) :: 
     .  rhs ! RHS input

        ! The estimate of the reciprocal condition number of the matrix.
        ! If rcond is less than the machine precision (in
        ! particular, if rcond = 0), the matrix is singular to working
        ! precision.  This condition is indicated by a return code of
        ! info > 0.
        real, intent(out) :: rcond 

        integer, intent(out) :: 
     .  code ! Used to determine when a decomp. failed

        ! Output variables
        real, intent(out), dimension(ndim,nrhs) :: 
     .  solution ! Solution

        ! Local Variables
        ! These contain the decomposition of the matrix
        real, dimension(ndim-1) :: dlf, duf
        real, dimension(ndim)   :: df
        real, dimension(ndim-2) :: du2
         
        integer, dimension(ndim) :: 
     .  ipivot  ! Index of pivots done during decomposition

        integer, dimension(ndim) :: 
     .  iwork   ! `scrap' array


        real, dimension(nrhs) :: 
     .  ferr, ! Forward error estimate
     .  berr  ! Backward error estimate

        real, dimension(3*ndim) :: 
     .  work  ! `Scrap' array

        integer :: info ! Diagnostic output

!-----------------------------------------------------------------------
!     *** The LAPACK Routine ***
!     SUBROUTINE SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,
!    $                   DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,
!    $                   WORK, IWORK, INFO )
!-----------------------------------------------------------------------

        if ( kind( diag(1) ) == 4 ) then
          call sgtsvx( "Not Factored", "No Transpose lhs", ndim, nrhs, 
     .                 subd(2:ndim), diag, supd(1:ndim-1), 
     .                 dlf, df, duf, du2, ipivot, 
     .                 rhs, ndim, solution, ndim, rcond,
     .                 ferr, berr, work, iwork, info )
        else if ( kind( diag(1) ) == 8 ) then
          call dgtsvx( "Not Factored", "No Transpose lhs", ndim, nrhs, 
     .                 subd(2:ndim), diag, supd(1:ndim-1), 
     .                 dlf, df, duf, du2, ipivot, 
     .                 rhs, ndim, solution, ndim, rcond,
     .                 ferr, berr, work, iwork, info )
        else
          stop "Tridag: Cannot resolve precision of real datatype"
        end if

        ! %% Debug
        ! where( ferr  >= 1.0e-3 ) print *, "tridag ferr=", ferr
        ! where( berr  >= 1.0e-3 ) print *, "tridag berr=", berr
        ! print *, "rcond", rcond
        ! %% End Debug

        select case( info )
        case( :-1 )
          write(fstderr,*) trim( solve_type )//
     .      "illegal value in argument", -info
          code = CLUBB_bad_lapack_arg

        case( 0 )
        ! Success!
        code = CLUBB_no_error

        case( 1: )
          if ( info == ndim+1 ) then
            write(fstderr,*) trim( solve_type) //
     .        " Warning: matrix is singular to working precision."
            write(fstderr,'(a,e12.5)') 
     .        "Estimate of the reciprocal of the"//
     .        " condition number: ", rcond
            code = CLUBB_no_error
          else
            write(fstderr,*) solve_type//
     .        " singular matrix."
            code = CLUBB_singular_matrix
          end if

        end select

        return
        end subroutine tridag_solvex

!-----------------------------------------------------------------------
        subroutine tridag_solve
     .             ( solve_type, ndim, nrhs, supd, diag, subd, 
     .               rhs, solution, code )

!       Description:
!       Solves a tridiagonal system of equations (simple routine)

!       References:
!       <http://www.netlib.org/lapack/single/sgtsv.f>
!       <http://www.netlib.org/lapack/double/dgtsv.f>
!-----------------------------------------------------------------------
        implicit none

        ! External
        external :: 
     .  sgtsv, ! Single-prec. General Tridiagonal Solver eXpert
     .  dgtsv  ! Double-prec. General Tridiagonal Solver eXpert

        intrinsic :: kind

        ! Input variables
        character(len=*), intent(in) :: 
     .  solve_type ! Used to write a message if this fails

        integer :: 
     .  ndim, ! N-dimension of matrix
     .  nrhs  ! # of right hand sides to back subst. after LU-decomp.

        ! Input/Output variables
        real, intent(inout), dimension(ndim) :: 
     .  diag,      ! Main diagonal
     .  subd, supd ! Sub and super diagonal

        real, intent(inout), dimension(ndim,nrhs) :: 
     .  rhs ! RHS input

        integer, intent(out) :: 
     .  code ! Used to determine when a decomp. failed

        ! Output variables
        real, intent(out), dimension(ndim,nrhs) :: 
     .  solution ! Solution

        ! Local Variables
         
        integer :: info ! Diagnostic output

!-----------------------------------------------------------------------
!       *** The LAPACK Routine ***
!       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!-----------------------------------------------------------------------

        if ( kind( diag(1) ) == 4 ) then
          call sgtsv( ndim, nrhs, subd(2:ndim), diag, supd(1:ndim-1), 
     .                rhs, ndim, info )
        else if ( kind( diag(1) ) == 8 ) then
          call dgtsv( ndim, nrhs, subd(2:ndim), diag, supd(1:ndim-1), 
     .                rhs, ndim, info )
        else
          stop "Tridag: Cannot resolve precision of real datatype"
        end if

        select case( info )
        case( :-1 )
          write(fstderr,*) trim( solve_type )//
     .      " illegal value in argument", -info
          code = CLUBB_bad_lapack_arg

        case( 0 )
          ! Success!
          solution = rhs
          code = CLUBB_no_error

        case( 1: )
          if ( info == ndim+1 ) then
            write(fstderr,*) trim( solve_type )//
     .        " matrix is singular to working precision."
            code = CLUBB_no_error
          else
            write(fstderr,*) trim( solve_type )//
     .        " singular matrix."
            code = CLUBB_singular_matrix
          end if

        end select

        return
        end subroutine tridag_solve

!-----------------------------------------------------------------------
        subroutine band_solvex( solve_type, nsup, nsub, ndim, nrhs, 
     .                          lhs, rhs, solution, rcond, code )
!       Description:
!       Restructure and then solve a band diagonal system, with
!       diagnostic output

!       References:
!       <http://www.netlib.org/lapack/single/sgbsvx.f>
!       <http://www.netlib.org/lapack/double/dgbsvx.f>

!       Notes:
!       I found that due to the use of sgbcon/dgbcon it is much
!       more expensive to use this on most systems than the simple 
!       driver. -dschanen 3 Sept 2006
!-----------------------------------------------------------------------

        implicit none

        ! External
        external :: 
     .  sgbsvx, ! Single-prec. General Band Solver eXpert
     .  dgbsvx  ! Double-prec. General Band Solver eXpert

        intrinsic :: eoshift, kind, trim

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        integer, intent(in) ::
     .  nsup, ! Number of superdiagonals
     .  nsub, ! Number of subdiagonals
     .  ndim, ! The order of the LHS Matrix, i.e. the # of linear equations
     .  nrhs  ! Number of RHS's to solve for

        real, dimension(nsup+nsub+1,ndim), intent(inout) :: 
     .  lhs ! Left hand side
        real, dimension(ndim,nrhs), intent(inout) :: 
     .  rhs ! Right hand side(s)

        ! Output Variables
        real, dimension(ndim,nrhs), intent(out) :: solution

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real, intent(out) :: 
     .  rcond

        integer, intent(out) :: code ! Valid calculation?

        ! Local Variables

        ! Workspaces
        real, dimension(3*ndim)  :: work
        integer, dimension(ndim) :: iwork

        real, dimension(2*nsub+nsup+1,ndim) ::
     .  lulhs ! LU Decomposition of the LHS

        integer, dimension(ndim) :: 
     .  ipivot

        real, dimension(nrhs) :: 
     .  ferr, berr ! Forward and backward error estimate

        real, dimension(ndim) :: 
     .  rscale, cscale ! Row and column scale factors for the LHS

        integer :: 
     .  info,  ! If this doesn't come back as 0, something went wrong
     .  offset,! Loop iterator
     .  imain  ! Main diagonal of the matrix

        character :: 
     .  equed ! Row equilibration status


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
        do offset = 1, nsub, 1
          lhs(imain+offset, 1:ndim)
     .    = eoshift( lhs(imain+offset, 1:ndim), offset )
        end do 

        ! Super diagonals
        do offset = 1, nsup, 1
          lhs(imain-offset, 1:ndim)
     .    = eoshift( lhs(imain-offset, 1:ndim), -offset )
        end do

!-----------------------------------------------------------------------
!     *** The LAPACK Routine ***
!     SUBROUTINE SGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,
!    $                   LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,
!    $                   RCOND, FERR, BERR, WORK, IWORK, INFO ) 
!-----------------------------------------------------------------------

        if ( kind( lhs(1,1) ) == 4 ) then
          call sgbsvx( 'Equilibrate lhs', 'No Transpose lhs',
     .                 ndim, nsub, nsup, nrhs,
     .                 lhs, nsup+nsub+1, lulhs, 2*nsub+nsup+1, 
     .                 ipivot, equed, rscale, cscale,
     .                 rhs, ndim, solution, ndim,
     .                 rcond, ferr, berr, work, iwork, info )
        else if ( kind( lhs(1,1) ) == 8 ) then
          call dgbsvx( 'Equilibrate lhs', 'No Transpose lhs',
     .                 ndim, nsub, nsup, nrhs,
     .                 lhs, nsup+nsub+1, lulhs, 2*nsub+nsup+1,
     .                 ipivot, equed, rscale, cscale,
     .                 rhs, ndim, solution, ndim,
     .                 rcond, ferr, berr, work, iwork, info )
        else
          stop "Cannot handle this precision" 
          ! One implication of this is that HOC cannot be used with quad
          ! precision variables without a quad precision band diagonal solver
        end if 

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

        select case( info )

        case( :-1 )
          write(fstderr,*) trim( solve_type )//
     .      " illegal value for argument", -info
          code = CLUBB_bad_lapack_arg

        case( 0 )
        ! Success!
        code = CLUBB_no_error

        case( 1: )
          if ( info == ndim+1 ) then
            write(fstderr,*) trim( solve_type )//
     .        " Warning: matrix singular to working precision."
            write(fstderr,'(a,e12.5)') 
     .        "Estimate of the reciprocal of the"//
     .        " condition number: ", rcond
            code = CLUBB_no_error
          else
            write(fstderr,*) trim( solve_type )//
     .        " band solver: singular matrix"
            code = CLUBB_singular_matrix
          end if

        end select

        return
        end subroutine band_solvex

!-----------------------------------------------------------------------
        subroutine band_solve( solve_type, nsup, nsub, ndim, nrhs, 
     .                          lhs, rhs, solution, code )
!       Description:
!       Restructure and then solve a band diagonal system

!       References:
!       <http://www.netlib.org/lapack/single/sgbsv.f>
!       <http://www.netlib.org/lapack/double/dgbsv.f>
!-----------------------------------------------------------------------

        implicit none

        ! External
        external :: 
     .  sgbsv, ! Single-prec. General Band Solver
     .  dgbsv  ! Double-prec. General Band Solver

        intrinsic :: eoshift, kind, trim

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        integer, intent(in) ::
     .  nsup, ! Number of superdiagonals
     .  nsub, ! Number of subdiagonals
     .  ndim, ! The order of the LHS Matrix, i.e. the # of linear equations
     .  nrhs  ! Number of RHS's to solve for

        ! Note: matrix lhs is intent(in), not intent(inout) 
        ! as in the subroutine band_solvex( )
        real, dimension(nsup+nsub+1,ndim), intent(in) :: 
     .  lhs ! Left hand side
        real, dimension(ndim,nrhs), intent(inout) :: 
     .  rhs ! Right hand side(s)

        ! Output Variables
        real, dimension(ndim,nrhs), intent(out) :: solution

        integer, intent(out) :: code ! Valid calculation?

        ! Local Variables

        ! Workspaces
        real, dimension(2*nsub+nsup+1,ndim) ::
     .  lulhs ! LU Decomposition of the LHS

        integer, dimension(ndim) :: 
     .  ipivot

        integer :: 
     .  info,  ! If this doesn't come back as 0, something went wrong
     .  offset,! Loop iterator
     .  imain  ! Main diagonal of the matrix

        ! Copy LHS into Decomposition scratch space

        lulhs(nsub+1:2*nsub+nsup+1, 1:ndim) = lhs(1:nsub+nsup+1, 1:ndim)

!-----------------------------------------------------------------------
!       Reorder LU Matrix to use LAPACK band matrix format

!       Shift example for lulhs matrix (note the extra bands):

!       [    +        +        +        +        +        +     ]    
!       [    +        +        +        +        +        +     ]    
!       [    *        *     lhs(1,1) lhs(1,2) lhs(1,3) lhs(1,4) ] (2)=>
!       [    *     lhs(2,1) lhs(2,2) lhs(2,3) lhs(2,4) lhs(2,5) ] (1)=>
!       [ lhs(3,1) lhs(3,2) lhs(3,3) lhs(3,4) lhs(3,5) lhs(3,6) ]
! <=(1) [ lhs(4,2) lhs(4,3) lhs(4,4) lhs(4,5) lhs(4,6)    *     ]
! <=(2) [ lhs(5,3) lhs(5,4) lhs(5,5) lhs(5,6)    *        *     ]
!       [    +        +        +        +        +        +     ]    
!       [    +        +        +        +        +        +     ]    

!       The '*' indicates unreferenced elements.
!       The '+' indicates an element overwritten during decomposition.
!       For additional bands above and below the main diagonal, the
!       shifts to the left or right increases by the distance from the
!       main diagonal of the matrix.
!-----------------------------------------------------------------------

        ! Reorder lulhs, omitting the additional 2*nsub bands
        ! that are used for the LU decomposition of the matrix.

        imain = nsub + nsup + 1

        ! For the offset, (+) is left, and (-) is right

        ! Sub diagonals
        do offset = 1, nsub, 1
          lulhs(imain+offset, 1:ndim)
     .    = eoshift( lulhs(imain+offset, 1:ndim), offset )
        end do 

        ! Super diagonals
        do offset = 1, nsup, 1
          lulhs(imain-offset, 1:ndim)
     .    = eoshift( lulhs(imain-offset, 1:ndim), -offset )
        end do

!-----------------------------------------------------------------------
!       *** LAPACK routine ***
!       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!-----------------------------------------------------------------------

        if ( kind( lhs(1,1) ) == 4 ) then
          call sgbsv( ndim, nsub, nsup, nrhs, lulhs, nsub*2+nsup+1, 
     .                ipivot, rhs, ndim, info )
        else if ( kind( lhs(1,1) ) == 8 ) then
          call dgbsv( ndim, nsub, nsup, nrhs, lulhs, nsub*2+nsup+1, 
     .                ipivot, rhs, ndim, info )
        else
          stop "Cannot handle this precision" 
          ! One implication of this is that HOC cannot be used with quad 
          ! precision variables without a quad precision band diagonal solver
        end if 

        select case( info )

        case( :-1 )
          write(fstderr,*) trim( solve_type )//
     .      " illegal value for argument ", -info
          code = CLUBB_bad_lapack_arg

        case( 0 )
        ! Success!
        solution = rhs
        code = CLUBB_no_error

        case( 1: )
          if ( info == ndim+1 ) then
            write(fstderr,*) trim( solve_type )//
     .        " Warning: matrix singular to working precision"
            code = CLUBB_no_error
          else
            write(fstderr,*) trim( solve_type )//
     .        " band solver: singular matrix"
            code = CLUBB_singular_matrix
          end if

        end select

        return
        end subroutine band_solve

!-----------------------------------------------------------------------

        end module lapack_wrap

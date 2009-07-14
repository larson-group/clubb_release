! $Id$
module matrix_operations

  implicit none

  public :: linear_eqn_solve, linear_symm_upper_eqn_solve, band_mult

  private ! Default scope

  contains

!----------------------------------------------------------------------
  subroutine linear_eqn_solve( n, a, x, b )
!   Description:
!    Solve for A * X = B for a general matrix.
!   References:
!     <http://www.netlib.org/lapack/double/dgesv.f> 
!-----------------------------------------------------------------------
    implicit none

    ! External
    external :: dgesv   ! LAPACK subroutine

    ! Parameters  
    integer, parameter :: nrhs = 1

    ! Input Variables
    integer, intent(in) :: n

    double precision, dimension(n,n), intent(in) :: a

    double precision, dimension(n), intent(in) :: b

    ! Output Variables
    double precision, dimension(n), intent(out) :: x

    ! Local Variables
    double precision, dimension(n,n) :: a_decomp

    integer, dimension(n) :: &
      ipiv ! Pivot indices from the permutation matrix

    integer :: info

    ! ---- Begin code ----

    a_decomp = a
    x = b

    call dgesv( n, nrhs, a_decomp, n, ipiv, x, n, info )

    select case( info )
    case( :-1 )
      write(0,*) "linear_eqn_solve" // & 
        " illegal value for argument ", -info
      stop
    case( 0 )
      ! Success!

    case( 1: )
      write(0,*) "linear_eqn_solve: singular matrix"
      stop
    end select

    return
  end subroutine linear_eqn_solve

!----------------------------------------------------------------------
  subroutine linear_symm_upper_eqn_solve( n, a, x, b )
!   Description:
!    Solve for A * X = B for a symmetric matrix, using the upper diagonals.
!   References:
!     <http://www.netlib.org/lapack/double/dsysv.f> 
!-----------------------------------------------------------------------
    implicit none

    ! External
    external :: dsysv   ! LAPACK subroutine

    ! Parameters  
    integer, parameter :: nrhs = 1

    ! Input Variables
    integer, intent(in) :: n

    double precision, dimension(n,n), intent(in) :: a

    double precision, dimension(n), intent(in) :: b

    ! Output Variables
    double precision, dimension(n), intent(out) :: x

    ! Local Variables
    double precision, dimension(n,n) :: a_decomp

    double precision, allocatable, dimension(:) :: work

    integer, dimension(n) :: &
      ipiv ! Pivot indices from the permutation matrix

    integer :: info, work_dim
!   integer :: i, j
    ! ---- Begin code ----

    work_dim = n * 128 ! Best guess for an optimal blocksize

    allocate( work(work_dim) )

    a_decomp = a
    x = b

!   do i = 1, n
!     do j = 1, n
!       write(6,'(e10.3)',advance='no') a(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause

    call dsysv( 'Upper', n, nrhs, a_decomp, n, ipiv, x, n, work, work_dim, info )

    select case( info )
    case( :-1 )
      write(0,*) "linear_symm_upper_eqn_solve" // & 
        " illegal value for argument ", -info
      stop
    case( 0 )
      ! Success!

    case( 1: )
      write(0,*) "linear_symm_upper_eqn_solve: singular matrix"
      stop
    end select

    deallocate( work )

    return
  end subroutine linear_symm_upper_eqn_solve

!-----------------------------------------------------------------------
  subroutine band_mult( trans, ndim, mdim, nsup, nsub, yinc, xinc, & 
                        alpha, beta, lhs, xvec, yvec )
!       Description:
!       Wrapper subroutine for banded matrix by vector multiplication in
!       the level 2 BLAS library.

!       References:
!       <http://www.netlib.org/blas/>
!-----------------------------------------------------------------------

    implicit none

    ! External
    ! Level 2 BLAS to multiply a vector by a band diagonal matrix
    external :: sgbmv, dgbmv

    intrinsic :: kind

    character(len=1), intent(in) ::  & 
      trans  ! Whether to use the transposition of the lhs matrix

    integer, intent(in) :: & 
      ndim, mdim,   & ! Dimensions of the matrix when not compact
      nsup, nsub,   & ! Super and Sub diagonals
      yinc, xinc   ! Increments of y and x vector

    real, intent(in) :: & 
      alpha,  & ! Coefficient of the matrix lhs
      beta      ! Coefficient of Y vector

    real, dimension(nsup+nsub+1, ndim), intent(in) :: & 
      lhs  ! The matrix 'A' in the blas subroutine

    real, dimension(ndim), intent(in) :: & 
      xvec ! The vector X

    real, dimension(ndim), intent(inout) :: & 
      yvec ! The vector Y


!-----------------------------------------------------------------------
!       *** BLAS 2 routine ***
!       SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!-----------------------------------------------------------------------
    ! Multiply so that Y := alpha*A*X + beta*Y

    if ( kind( lhs(1,1) ) == 4 ) then
      call sgbmv( trans, ndim, mdim, nsub, nsup,  & 
                  alpha, lhs, nsup+nsub+1,  & 
                  xvec, xinc, beta, yvec, yinc )

    else if ( kind( lhs(1,1) ) == 8 ) then
      call dgbmv( trans, ndim, mdim, nsub, nsup,  & 
                  alpha, lhs, nsup+nsub+1,  & 
                  xvec, xinc, beta, yvec, yinc )
    else
      stop "Cannot multiply this precision"
    end if

    return
  end subroutine band_mult
!-----------------------------------------------------------------------

end module matrix_operations

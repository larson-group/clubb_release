! $Id$
module matrix_operations

  implicit none

  interface linear_symm_upper_eqn_solve
    module procedure sp_lin_symm_upper_eqn_solve, dp_lin_symm_upper_eqn_solve
  end interface linear_symm_upper_eqn_solve

  public :: linear_eqn_solve, linear_symm_upper_eqn_solve, band_mult, &
    covar_matrix_2_corr_matrix, Cholesky_factor

  private :: sp_lin_symm_upper_eqn_solve, dp_lin_symm_upper_eqn_solve

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
  subroutine dp_lin_symm_upper_eqn_solve( n, a, x, b )
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
      write(0,*) "dp_lin_symm_upper_eqn_solve" // & 
        " illegal value for argument ", -info
      stop
    case( 0 )
      ! Success!

    case( 1: )
      write(0,*) "dp_lin_symm_upper_eqn_solve: singular matrix"
      stop
    end select

    deallocate( work )

    return
  end subroutine dp_lin_symm_upper_eqn_solve
!----------------------------------------------------------------------
  subroutine sp_lin_symm_upper_eqn_solve( n, a, x, b )
!   Description:
!    Solve for A * X = B for a symmetric matrix, using the upper diagonals.
!   References:
!     <http://www.netlib.org/lapack/single/ssysv.f>
!-----------------------------------------------------------------------
    implicit none

    ! External
    external :: ssysv   ! LAPACK subroutine

    ! Parameters
    integer, parameter :: nrhs = 1

    ! Input Variables
    integer, intent(in) :: n

    real(kind=4), dimension(n,n), intent(in) :: a

    real(kind=4), dimension(n), intent(in) :: b

    ! Output Variables
    real(kind=4), dimension(n), intent(out) :: x

    ! Local Variables
    real(kind=4), dimension(n,n) :: a_decomp

    real(kind=4), allocatable, dimension(:) :: work

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

    call ssysv( 'Upper', n, nrhs, a_decomp, n, ipiv, x, n, work, work_dim, info )

    select case( info )
    case( :-1 )
      write(0,*) "sp_lin_symm_upper_eqn_solve" // & 
        " illegal value for argument ", -info
      stop
    case( 0 )
      ! Success!

    case( 1: )
      write(0,*) "sp_lin_symm_upper_eqn_solve: singular matrix"
      stop
    end select

    deallocate( work )

    return
  end subroutine sp_lin_symm_upper_eqn_solve

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
  subroutine covar_matrix_2_corr_matrix( ndim, cov, corr )

! Description:
!   Convert a matrix of covariances in to a matrix of correlations
! References:
!   None
!-----------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: sqrt

    ! Input Variables
    integer, intent(in) :: ndim

    double precision, dimension(ndim,ndim), intent(in) :: &
      cov ! Covariance Matrix [units vary]

    ! Output Variables
    double precision, dimension(ndim,ndim), intent(out) :: & 
      corr ! Correlation Matrix [-]

    ! Local Variables
    integer :: i, j

    ! ---- Begin Code ----
    i = 1   ! These 4 lines eliminate a g95 compiler warning of uninitialized variables.
    i = i   ! Simply initializing them to 1 is not sufficient as it generates two new
    j = 1   ! warnings.  -meyern
    j = j

    forall( i = 1:ndim, j= 1:ndim )
    corr(i,j) = cov(i,j) / sqrt( cov(i,i) * cov(j,j) )
    end forall

    return
  end subroutine covar_matrix_2_corr_matrix

!----------------------------------------------------------------------
  subroutine Cholesky_factor( n, a, a_Cholesky )
!  Description:
!    Create a Cholesky factorization of a.
!  References:
!    <http://www.netlib.org/lapack/explore-html-old/dpotrf.f.html>
!-----------------------------------------------------------------------
    implicit none

    ! External
    external :: dpotrf ! LAPACK subroutine

    ! Constant Parameters
    integer, parameter :: itermax = 5

    ! Input Variables
    integer, intent(in) :: n

    double precision, dimension(n,n), intent(in) :: a

    ! Output Variables
    double precision, dimension(n,n), intent(out) :: a_Cholesky

    ! Local Variables
    double precision :: tau

    integer :: info
    integer :: i, j, iter

    ! ---- Begin code ----

    a_Cholesky = a ! Copy input array into output array
!   do i = 1, n
!     do j = 1, n
!       write(6,'(e10.3)',advance='no') a(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause

    do iter = 1, itermax
      call dpotrf( 'Lower', n, a_Cholesky, n, info )

      select case( info )
      case( :-1 )
        write(0,*) "Cholesky_factor " // & 
          " illegal value for argument ", -info
        stop
      case( 0 )
        ! Success!
        exit
      case( 1: )
        write(0,*) "Cholesky_factor: leading minor of order ", info, " is not positive definite."
        write(0,*) "factorization failed."
!       write(6,*) "a="
!       do i = 1, n
!         do j = 1, n
!           write(6,'(e10.3)',advance='no') a(i,j)
!         end do
!         write(6,*) ""
!       end do
!       write(6,*) "a_factor="
!       do i = 1, n
!         do j = 1, n
!           write(6,'(e10.3)',advance='no') a_Cholesky(i,j)
!         end do
!         write(6,*) ""
!       end do
        if ( iter == itermax ) then
          write(0,*) "iteration =", iter, "itermax =", itermax
          stop
        else
          write(0,*) "Attempting to modify matrix to allow factorization."
        end if

        ! The number used for tau here is case specific to the Sigma covariance
        ! matrix in the latin hypercube code and is not at all general.  
        ! Tau should be number that is small relative to the other diagonal 
        ! elements of the matrix to have keep the error caused by modifying 'a' low.
        ! -dschanen 30 Aug 2010
        tau = a(1,1) * iter ! Use the s_mellor element * iteration for now

        do i = 1, n
          do j = 1, n
            if ( i == j ) then
              a_Cholesky(i,j) = a(i,j) + tau ! Add tau to the diagonal 
            else
              a_Cholesky(i,j) = a(i,j)
            end if
          end do
        end do

      end select ! info

    end do ! 1..itermax

    return
  end subroutine Cholesky_factor

end module matrix_operations

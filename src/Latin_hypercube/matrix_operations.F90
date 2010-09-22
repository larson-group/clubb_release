! $Id$
module matrix_operations

  implicit none


  public :: symm_covar_matrix_2_corr_matrix, Cholesky_factor

  private :: Symm_matrix_eigenvalues

  private ! Default scope

  contains
!
!-----------------------------------------------------------------------
  subroutine symm_covar_matrix_2_corr_matrix( ndim, cov, corr )

! Description:
!   Convert a matrix of covariances in to a matrix of correlations.
!   This only does the computation the lower triangular portion of the 
!   matrix.
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
    integer :: row, col

    ! ---- Begin Code ----

    corr = 0. ! Initialize to 0

    do col = 1, ndim
      do row = col, ndim
        corr(row,col) = cov(row,col) / sqrt( cov(row,col) * cov(row,col) )
      end do
    end do

    return
  end subroutine symm_covar_matrix_2_corr_matrix

!----------------------------------------------------------------------
  subroutine Cholesky_factor( ndim, a_input, a_scaling, a_Cholesky, l_scaled )
!  Description:
!    Create a Cholesky factorization of a_input.
!  References:
!    <http://www.netlib.org/lapack/explore-html/a00868.html> dpotrf
!    <http://www.netlib.org/lapack/explore-html/a00860.html> dpoequ
!    <http://www.netlib.org/lapack/explore-html/a00753.html> dlaqsy
!-----------------------------------------------------------------------
    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use constants_clubb, only: &
      fstderr ! Constant

    implicit none

    ! External
    external :: dpotrf, dpoequ, dlaqsy ! LAPACK subroutines

    ! Constant Parameters
    integer, parameter :: itermax = 5

    ! Input Variables
    integer, intent(in) :: ndim

    double precision, dimension(ndim,ndim), intent(in) :: a_input

    ! Output Variables
    double precision, dimension(ndim), intent(out) :: a_scaling

    double precision, dimension(ndim,ndim), intent(out) :: a_Cholesky

    logical, intent(out) :: l_scaled

    ! Local Variables
    double precision, dimension(ndim) :: a_eigenvalues
    double precision, dimension(ndim,ndim) ::  a_corr, a_scaled

    double precision :: tau, d_smallest

    double precision :: amax, scond
    integer :: info
    integer :: i, j, iter

    character :: equed

    ! ---- Begin code ----

    a_scaled = a_input ! Copy input array into output array

!   do i = 1, n
!     do j = 1, n
!       write(6,'(e10.3)',advance='no') a(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause

    equed = 'N'

    ! Compute scaling for a_input
    call dpoequ( ndim, a_input, ndim, a_scaling, scond, amax, info )

    if ( info == 0 ) then
      ! Apply scaling to a_input
      call dlaqsy( 'Lower', ndim, a_scaled, ndim, a_scaling, scond, amax, equed )
    end if

    ! Determine if scaling was necessary
    if ( equed == 'Y' ) then
      l_scaled = .true.
      a_Cholesky = a_scaled
    else
      l_scaled = .false.
      a_Cholesky = a_input
    end if

    do iter = 1, itermax
      call dpotrf( 'Lower', ndim, a_Cholesky, ndim, info )

      select case( info )
      case( :-1 )
        write(fstderr,*) "Cholesky_factor " // & 
          " illegal value for argument ", -info
        stop
      case( 0 )
        ! Success!
        if ( clubb_at_least_debug_level( 1 ) .and. iter > 1 ) then
          write(fstderr,*) "a_factored (worked)="
          do i = 1, ndim
            do j = 1, ndim
              write(fstderr,'(e10.3)',advance='no') a_Cholesky(i,j)
            end do
            write(fstderr,*) ""
          end do
        end if
        exit
      case( 1: )
        if ( clubb_at_least_debug_level( 1 ) ) then
          ! This shouldn't happen now that the s and t Mellor elements have been
          ! modified to never be perfectly correlated, but it's here just in case.
          ! -dschanen 10 Sept 2010
          write(fstderr,*) "Cholesky_factor: leading minor of order ", &
            info, " is not positive definite."
          write(fstderr,*) "factorization failed."
          write(fstderr,*) "a_input="
          do i = 1, ndim
            do j = 1, ndim
              write(fstderr,'(e10.3)',advance='no') a_input(i,j)
            end do
            write(fstderr,*) ""
          end do
          write(fstderr,*) "a_Cholesky="
          do i = 1, ndim
            do j = 1, ndim
              write(fstderr,'(e10.3)',advance='no') a_Cholesky(i,j)
            end do
            write(fstderr,*) ""
          end do
        end if

        if ( clubb_at_least_debug_level( 2 ) ) then
          call Symm_matrix_eigenvalues( ndim, a_input, a_eigenvalues )
          write(fstderr,*) "a_eigenvalues="
          do i = 1, ndim
            write(fstderr,'(e10.3)',advance='no') a_eigenvalues(i)
          end do
          write(fstderr,*) ""

          call symm_covar_matrix_2_corr_matrix( ndim, a_input, a_corr )
          write(fstderr,*) "a_correlations="
          do i = 1, ndim
            do j = 1, ndim
              write(fstderr,'(g10.3)',advance='no') a_corr(i,j)
            end do
            write(fstderr,*) ""
          end do
        end if

        if ( iter == itermax ) then
          write(fstderr,*) "iteration =", iter, "itermax =", itermax
          stop "Fatal error in Cholesky_factor"
        else if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) "Attempting to modify matrix to allow factorization."
        end if

        if ( l_scaled ) then
          a_Cholesky = a_scaled
        else
          a_Cholesky = a_input
        end if
        ! The number used for tau here is case specific to the Sigma covariance
        ! matrix in the latin hypercube code and is not at all general.
        ! Tau should be number that is small relative to the other diagonal
        ! elements of the matrix to have keep the error caused by modifying 'a' low.
        ! -dschanen 30 Aug 2010
        d_smallest = a_Cholesky(1,1)
        do i = 2, ndim
          if ( d_smallest > a_Cholesky(i,i) ) d_smallest = a_Cholesky(i,i)
        end do
        tau = d_smallest * 0.01 * dble( iter ) ! Use the smallest element * 0.01 * iteration
!       print *, "tau =", tau, "d_smallest = ", d_smallest

        do i = 1, ndim
          do j = 1, ndim
            if ( i == j ) then
              a_Cholesky(i,j) = a_Cholesky(i,j) + tau ! Add tau to the diagonal
            else
              a_Cholesky(i,j) = a_Cholesky(i,j)
            end if
          end do
        end do

        if ( clubb_at_least_debug_level( 2 ) ) then
          call Symm_matrix_eigenvalues( ndim, a_Cholesky, a_eigenvalues )
          write(fstderr,*) "a_modified eigenvalues="
          do i = 1, ndim
            write(fstderr,'(e10.3)',advance='no') a_eigenvalues(i)
          end do
          write(fstderr,*) ""
        end if

      end select ! info

    end do ! 1..itermax

    return
  end subroutine Cholesky_factor

!----------------------------------------------------------------------
  subroutine Symm_matrix_eigenvalues( ndim, a_input, a_eigenvalues )
!   Description:
!   References:
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Constant

    implicit none

    ! External
    external :: dsyev ! LAPACK subroutine

    ! Parameters
    integer, parameter :: &
      lwork = 180 ! This is the optimal value I obtained for an n of 5 -dschanen 31 Aug 2010

    ! Input Variables
    integer, intent(in) :: ndim

    double precision, dimension(ndim,ndim), intent(in) :: a_input

    ! Output Variables
    double precision, dimension(ndim), intent(out) :: a_eigenvalues

    ! Local Variables
    double precision, dimension(ndim,ndim) :: a_scratch

    double precision, dimension(lwork) :: work

    integer :: info
!   integer :: i, j
    ! ---- Begin code ----

    a_scratch = a_input

!   do i = 1, ndim
!     do j = 1, ndim
!       write(6,'(e10.3)',advance='no') a(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause

    call dsyev( 'No eigenvectors', 'Lower', ndim, a_scratch, ndim, &
                a_eigenvalues, work, lwork, info )

    select case( info )
    case( :-1 )
      write(fstderr,*) "Symm_matrix_eigenvalues:" // & 
        " illegal value for argument ", -info
      stop
    case( 0 )
      ! Success!

    case( 1: )
      write(fstderr,*) "Symm_matrix_eigenvalues: Algorithm failed to converge."
      stop
    end select

    return
  end subroutine Symm_matrix_eigenvalues

end module matrix_operations

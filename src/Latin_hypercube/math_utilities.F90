!$Id$
module math_utilities
!-----------------------------------------------------------------------
! Various mathematical utilities
!-----------------------------------------------------------------------
  implicit none

  public :: corrcoef, std, cov, mean, compute_variance

  private

  contains

!-----------------------------------------------------------------------
  function corrcoef( vect1, vect2, n )

! Description:
!   Correlation coefficient of two vectors

! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! Input
    integer, intent(in) :: n

    double precision, dimension(n), intent(in) :: &
      vect1, vect2

    ! Return type
    double precision :: corrcoef

    corrcoef = cov( vect1, vect2, n ) / & 
           sqrt( cov( vect1, vect1, n ) * cov( vect2, vect2, n ) )

    return
  end function corrcoef

!-----------------------------------------------------------------------
  function std( vector, n )
! Description:
!   Compute standard deviation of vector
! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: n

    double precision, dimension(n), intent(in) :: &
      vector

    ! Return type
    double precision std

    std = sqrt( cov( vector, vector, n )*( n/(n-1) ) )

    return
  end function std


!-----------------------------------------------------------------------
  function cov( vect1, vect2, n )

! Description:
!   Covariance of two vectors
! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: n

    double precision, dimension(n), intent(in) :: &
      vect1, vect2

    ! Return type
    double precision :: cov

    ! Internal
    double precision :: sum, avg1, avg2
    integer :: j

    avg1 = mean( vect1, n )
    avg2 = mean( vect2, n )

    sum = 0.d0
    do j = 1, n
      sum = sum + (vect1(j) - avg1) * (vect2(j) - avg2)
    enddo

    cov = sum / n

    return
  end function cov

!-----------------------------------------------------------------------
  function mean( vector, n )
! Description:
!   Find the mean of the vector

! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: dble, sum

    ! Input Varibles
    integer, intent(in) :: n

    double precision, dimension(n), intent(in) :: &
      vector

    ! Return type
    double precision :: mean

    ! ---- Begin Code ----

    mean = sum( vector ) / dble( n )

    return
  end function mean
!-----------------------------------------------------------------------
  function compute_variance( n_pts, n_samples, x_sample, x_mean ) result( variance )

! Description:
!   Compute variance of a set of sample points

! References:
!   None
!-----------------------------------------------------------------------
    implicit none

    integer, intent(in) :: &
      n_pts, &   ! Number of data points in the mean / variance
      n_samples  ! Number of sample points compute the variance of

    real,dimension(n_pts,n_samples) :: &
      x_sample ! Collection of sample points

    real,dimension(n_pts) :: &
      x_mean  ! Mean sample points

    real,dimension(n_pts) :: &
      variance ! Variance of x

    integer :: sample ! Loop iterator

    ! ---- Begin Code ----

    variance(1:n_pts) = 0.0

    do sample=1, n_samples
      variance(1:n_pts) = variance(1:n_pts) + ( x_sample(1:n_pts,sample) - x_mean(1:n_pts) )**2
    end do

    variance(1:n_pts) = variance(1:n_pts) / real( n_samples )

    return
  end function compute_variance

end module math_utilities

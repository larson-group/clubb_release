!$Id$
module math_utilities
!-----------------------------------------------------------------------
! Various mathematical utilities
!-----------------------------------------------------------------------
  implicit none

  public :: corrcoef, std, cov, compute_mean, compute_sample_variance

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

    ! External
    intrinsic :: sqrt

    ! Input
    integer, intent(in) :: n

    real, dimension(n), intent(in) :: &
      vect1, vect2

    ! Return type
    real :: corrcoef

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

    ! External
    intrinsic :: sqrt

    ! Input Variables
    integer, intent(in) :: n

    real, dimension(n), intent(in) :: &
      vector

    ! Return type
    real std

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

    real, dimension(n), intent(in) :: &
      vect1, vect2

    ! Return type
    real :: cov

    ! Internal
    real :: sum, avg1, avg2
    integer :: j

    avg1 = compute_mean( n, vect1 )
    avg2 = compute_mean( n, vect2 )

    sum = 0.
    do j = 1, n
      sum = sum + (vect1(j) - avg1) * (vect2(j) - avg2)
    enddo

    cov = sum / real( n )

    return
  end function cov

!-----------------------------------------------------------------------
  pure function compute_mean( n_dim, vector )
! Description:
!   Find the mean of the vector

! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: real, sum

    ! Input Varibles
    integer, intent(in) :: n_dim

    real, dimension(n_dim), intent(in) :: &
      vector

    ! Return type
    real :: compute_mean

    ! ---- Begin Code ----

    compute_mean = sum( vector ) / real( n_dim )

    return
  end function compute_mean
!-----------------------------------------------------------------------
  pure function compute_sample_variance( n_levels, n_samples, x_sample, weight, x_mean ) &
    result( variance )

! Description:
!   Compute the variance of a set of sample points

! References:
!   None
!-----------------------------------------------------------------------
    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples   ! Number of sample points compute the variance of

    real,dimension(n_levels,n_samples), intent(in) :: &
      x_sample ! Collection of sample points    [units vary]

    real,dimension(n_samples), intent(in) :: &
      weight ! Coefficient to weight the nth sample point by [-]

    real,dimension(n_levels), intent(in) :: &
      x_mean ! Mean sample points [units vary]

    ! Output Variable
    real,dimension(n_levels) :: &
      variance ! Variance of x [(units vary)^2]

    ! Local Variable(s)
    integer :: sample ! Loop iterator

    ! ---- Begin Code ----

    variance(1:n_levels) = 0.0

    do sample=1, n_samples
      variance(1:n_levels) = variance(1:n_levels) &
        + weight(sample) * ( x_sample(1:n_levels,sample) - x_mean(1:n_levels) )**2
    end do

    variance(1:n_levels) = variance(1:n_levels) / real( n_samples )

    return
  end function compute_sample_variance

end module math_utilities

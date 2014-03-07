!$Id$
module math_utilities
!-----------------------------------------------------------------------
! Various mathematical utilities
!-----------------------------------------------------------------------
  implicit none

  public :: corrcoef, std, covar, compute_sample_mean, & 
            compute_sample_variance, compute_sample_covariance

  private

  contains

!-----------------------------------------------------------------------
  function corrcoef( vect1, vect2, n )

! Description:
!   Correlation coefficient of two vectors

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt

    ! Input
    integer, intent(in) :: n

    real( kind = core_rknd ), dimension(n), intent(in) :: &
      vect1, vect2

    ! Return type
    real(kind=core_rknd) :: corrcoef

    corrcoef = covar( vect1, vect2, n ) / & 
           sqrt( covar( vect1, vect1, n ) * covar( vect2, vect2, n ) )

    return
  end function corrcoef

!-----------------------------------------------------------------------
  function std( vector, n )
! Description:
!   Compute standard deviation of vector
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt

    ! Input Variables
    integer, intent(in) :: n

    real( kind = core_rknd ), dimension(n), intent(in) :: &
      vector

    ! Return type
    real(kind=core_rknd) :: std

    std = sqrt( covar( vector, vector, n )*real( n/(n-1), kind=core_rknd ) )

    return
  end function std


!-----------------------------------------------------------------------
  function covar( vect1, vect2, n )

! Description:
!   Covariance of two vectors
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: n

    real( kind = core_rknd ), dimension(n), intent(in) :: &
      vect1, vect2

    ! Return type
    real(kind=core_rknd) :: covar

    ! Local variables
    real(kind=core_rknd) :: total, avg1, avg2
    integer :: j

    ! ---- Begin Code ----

    avg1 = sum( vect1 ) / real( n, kind = core_rknd ) 
    avg2 = sum( vect2 ) / real( n, kind = core_rknd ) 

    total = 0._core_rknd
    do j = 1, n
      total = total + (vect1(j) - avg1) * (vect2(j) - avg2)
    enddo

    covar = total / real( n, kind = core_rknd )

    return
  end function covar

!-----------------------------------------------------------------------
  pure function compute_sample_mean( n_levels, n_samples, weight, x_sample ) &
    result( mean )
! Description:
!   Find the mean of a set of sample points

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: real, sum

    ! Input Varibles
    integer, intent(in) :: &
      n_levels, &
      n_samples

    real( kind = core_rknd ),dimension(n_levels,n_samples), intent(in) :: &
      x_sample ! Collection of sample points    [units vary]

    real( kind = core_rknd ), dimension(n_samples), intent(in) :: &
      weight   ! Weights for individual points of the vector

    ! Return type
    real( kind = core_rknd ), dimension(n_levels) :: mean

    integer :: k

    ! ---- Begin Code ----

    ! Get rid of an annoying compiler warning.
    k = 1
    k = k

    forall( k = 1:n_levels )
      mean(k) = sum( weight(1:n_samples) * x_sample(k,1:n_samples) ) &
              / real( n_samples, kind=core_rknd )
    end forall


    return

  end function compute_sample_mean

!-----------------------------------------------------------------------
  pure function compute_sample_variance( n_levels, n_samples, x_sample, weight, x_mean ) &
    result( variance )

! Description:
!   Compute the variance of a set of sample points

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples   ! Number of sample points compute the variance of

    real( kind = core_rknd ),dimension(n_levels,n_samples), intent(in) :: &
      x_sample ! Collection of sample points    [units vary]

    real( kind = core_rknd ),dimension(n_samples), intent(in) :: &
      weight ! Coefficient to weight the nth sample point by [-]

    real( kind = core_rknd ),dimension(n_levels), intent(in) :: &
      x_mean ! Mean sample points [units vary]

    ! Output Variable
    real( kind = core_rknd ),dimension(n_levels) :: &
      variance ! Variance of x [(units vary)^2]

    ! Local Variable(s)
    integer :: sample ! Loop iterator

    ! ---- Begin Code ----

    variance(1:n_levels) = 0.0_core_rknd

    do sample=1, n_samples
      variance(1:n_levels) = variance(1:n_levels) &
        + weight(sample) * ( x_sample(1:n_levels,sample) - x_mean(1:n_levels) )**2
    end do

    variance(1:n_levels) = variance(1:n_levels) / real( n_samples, kind=core_rknd )

    return
  end function compute_sample_variance

!-----------------------------------------------------------------------
  pure function compute_sample_covariance( n_levels, n_samples, weight, &
                   x_sample, x_mean, y_sample, y_mean ) &
    result( covariance )

! Description:
!   Compute the covariance of a set of sample points of 2 variables
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples   ! Number of sample points compute the variance of

    real( kind = core_rknd ),dimension(n_levels,n_samples), intent(in) :: &
      x_sample, & ! Collection of sample points    [units vary]
      y_sample

    real( kind = core_rknd ),dimension(n_samples), intent(in) :: &
      weight ! Coefficient to weight the nth sample point by [-]

    real( kind = core_rknd ),dimension(n_levels), intent(in) :: &
      x_mean, & ! Mean sample points [units vary]
      y_mean

    ! Output Variable
    real( kind = core_rknd ),dimension(n_levels) :: &
      covariance ! Coariance of x and y [(units vary)^2]

    ! Local Variable(s)
    integer :: sample ! Loop iterator

    ! ---- Begin Code ----

    covariance(1:n_levels) = 0.0_core_rknd

    do sample=1, n_samples
      covariance(1:n_levels) = covariance(1:n_levels) &
        + weight(sample) * ( x_sample(1:n_levels,sample) - x_mean(1:n_levels) ) &
           * ( y_sample(1:n_levels,sample) - y_mean(1:n_levels) )
    end do

    covariance(1:n_levels) = covariance(1:n_levels) / real( n_samples, kind=core_rknd )

    return
  end function compute_sample_covariance

end module math_utilities

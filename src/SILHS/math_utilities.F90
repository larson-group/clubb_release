!----------------------------------------------------------------------- 
!$Id$
!===============================================================================
module math_utilities
!-----------------------------------------------------------------------
! Various mathematical utilities
!-----------------------------------------------------------------------
  implicit none

  public :: compute_sample_mean, compute_sample_variance, compute_sample_covariance, &
            rand_integer_in_range

  private

  contains

!-----------------------------------------------------------------------
  pure function compute_sample_mean( n_levels, n_samples, ngrdcol, &
                                     weight, x_sample ) result( mean )
! Description:
!   Find the mean of a set of sample points in every model column.
!
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, &  ! Number of sample levels
      n_samples, & ! Number of sample points
      ngrdcol      ! Number of model columns

    real( kind = core_rknd ), dimension(ngrdcol,n_samples,n_levels), intent(in) :: &
      weight, & ! Weights for individual points of the vector
      x_sample ! Collection of sample points    [units vary]

    ! Return type
    real( kind = core_rknd ), dimension(ngrdcol,n_levels) :: mean ! Mean sample points [units vary]

    ! Local Variables
    integer :: i, k, sample ! Column, level, and sample loop iterators

    ! ---- Begin Code ----

    mean = 0._core_rknd

    ! Preserve the sample accumulation order for each column while making the
    ! column index the contiguous, innermost loop.
    do sample = 1, n_samples
      do k = 1, n_levels
        do i = 1, ngrdcol
          mean(i,k) = mean(i,k) + weight(i,sample,k) * x_sample(i,sample,k)
        enddo
      enddo
    enddo

    mean = mean / real( n_samples, kind=core_rknd )

  end function compute_sample_mean


!-----------------------------------------------------------------------
  pure function compute_sample_variance( n_levels, n_samples, ngrdcol, &
                                         x_sample, weight, x_mean ) result( variance )
! Description:
!   Compute the variance of a set of sample points in every model column.
!
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples, & ! Number of sample points to compute the variance of
      ngrdcol      ! Number of model columns

    real( kind = core_rknd ), dimension(ngrdcol,n_samples,n_levels), intent(in) :: &
      x_sample, &          ! Collection of sample points    [units vary]
      weight    ! Coefficient to weight the nth sample point by [-]

    real( kind = core_rknd ), dimension(ngrdcol,n_levels), intent(in) :: &
      x_mean ! Mean sample points [units vary]

    ! Return type
    real( kind = core_rknd ), dimension(ngrdcol,n_levels) :: variance ! Variance of x [(units vary)^2]

    ! Local Variables
    integer :: i, k, sample ! Column, level, and sample loop iterators

    ! ---- Begin Code ----

    variance = 0._core_rknd

    do sample = 1, n_samples
      do k = 1, n_levels
        do i = 1, ngrdcol
          variance(i,k) = variance(i,k) + weight(i,sample,k) &
            * ( x_sample(i,sample,k) - x_mean(i,k) )**2
        enddo
      enddo
    enddo

    variance = variance / real( n_samples, kind=core_rknd )

  end function compute_sample_variance

!-----------------------------------------------------------------------
  pure function compute_sample_covariance( n_levels, n_samples, ngrdcol, &
                                           weight, x_sample, x_mean, &
                                           y_sample, y_mean ) result( covariance )
! Description:
!   Compute the covariance of a set of sample points of 2 variables
!   in every model column.
!
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples, & ! Number of sample points to compute the covariance of
      ngrdcol      ! Number of model columns

    real( kind = core_rknd ), dimension(ngrdcol,n_samples,n_levels), intent(in) :: &
      weight, &   ! Coefficient to weight the nth sample point by [-]
      x_sample, & ! Collection of sample points    [units vary]
      y_sample    ! Collection of sample points    [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,n_levels), intent(in) :: &
      x_mean, & ! Mean sample points [units vary]
      y_mean    ! Mean sample points [units vary]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,n_levels) :: covariance ! Covariance of x and y [(units vary)^2]

    ! Local Variables
    integer :: i, k, sample ! Column, level, and sample loop iterators

    ! ---- Begin Code ----

    covariance = 0._core_rknd

    do sample = 1, n_samples
      do k = 1, n_levels
        do i = 1, ngrdcol
          covariance(i,k) = covariance(i,k) + weight(i,sample,k) &
            * ( x_sample(i,sample,k) - x_mean(i,k) ) &
            * ( y_sample(i,sample,k) - y_mean(i,k) )
        enddo
      enddo
    enddo

    covariance = covariance / real( n_samples, kind=core_rknd )

  end function compute_sample_covariance

  !-----------------------------------------------------------------------
  function rand_integer_in_range(low, high)

  ! Description:
  !   Returns a uniformly distributed integer in the range [low,high]
  !   using the Mersenne Twister PRNG library.
  !
  !   The integers returned from this function are actually not quite
  !   evenly distributed because of the use of MOD. Smaller numbers are
  !   slightly more likely than larger ones. This could be fixed someday.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use mt95, only: &
      genrand_intg, & ! Constant
      genrand_int32   ! Procedure

    implicit none

    ! Local Constants

    ! Input Variables
    integer, intent(in) :: &
      low,   &      ! Lowest possible returned value
      high          ! Highest possible returned value

    ! Output Variable
    integer :: &
      rand_integer_in_range  ! Random integer in range [low,high]

    ! Local Variables
    integer( kind = genrand_intg ) :: &
      rand_32                ! Random integer in range[-2^31, +2^31-1]

    integer :: &
      range_width

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    range_width = high - low + 1
    call genrand_int32( rand_32 )
    rand_integer_in_range = abs( mod( rand_32, range_width ) ) + low

    return
  end function rand_integer_in_range
  !-----------------------------------------------------------------------


end module math_utilities

! $Id$
!===============================================================================
module KK_utilities

  private ! Set default scope to private

  public :: mean_L2N,   &
            stdev_L2N,  &
            corr_NL2NN, &
            corr_LL2NN, &
            factorial,  &
            Dv_fnc        ! Parabolic Cylinder Function, D.

  contains

  !=============================================================================
  pure function mean_L2N( mu_x, sigma_sqd_x )  &
  result( mu_x_n )
  
    ! Description:
    ! For a lognormally-distributed variable x, this function finds the mean of
    ! ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
    ! and the variance of x (sigma_sqd_x) for the ith component of the PDF.
    ! The value ln x is distributed normally when x is distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) ::  &
      mu_x,        & ! Mean of x (ith PDF component)       [-]
      sigma_sqd_x    ! Variance of x (ith PDF component)   [-]

    ! Return Variable
    real ::  &
      mu_x_n  ! Mean of ln x (ith PDF component)           [-]

    ! Find the mean of ln x for the ith component of the PDF.
    mu_x_n = log( mu_x * ( 1.0 + sigma_sqd_x / mu_x**2.0 )**(-0.5) )

    return
  end function mean_L2N

  !=============================================================================
  pure function stdev_L2N( mu_x, sigma_sqd_x )  &
  result( sigma_x_n )

    ! Description:
    ! For a lognormally-distributed variable x, this function finds the standard
    ! deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
    ! mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
    ! of the PDF.  The value ln x is distributed normally when x is distributed
    ! lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) ::  &
      mu_x,        & ! Mean of x (ith PDF component)       [-]
      sigma_sqd_x    ! Variance of x (ith PDF component)   [-]

    ! Return Variable
    real ::  &
      sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]

    ! Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = sqrt( log( 1.0 + sigma_sqd_x / mu_x**2.0 ) )

    return
  end function stdev_L2N

  !=============================================================================
  pure function corr_NL2NN( corr_xy, sigma_y_n )  &
  result( corr_xy_n )

    ! Description:
    ! For a normally-distributed variable x and a lognormally-distributed
    ! variable y, this function finds the correlation between x and ln y
    ! (corr_xy_n) for the ith component of the PDF, given the correlation
    ! between x and y (corr_xy) and the standard deviation of ln y (sigma_y_n)
    ! for the ith component of the PDF.  The value ln y is distributed normally
    ! when y is distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. B-1.
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: &
      corr_xy,   & ! Correlation between x and y (ith PDF component)  [-]
      sigma_y_n    ! Standard deviation of ln y (ith PDF component)   [-]

    ! Return Variable
    real ::  &
      corr_xy_n  ! Correlation between x and ln y (ith PDF component) [-]

    ! Find the correlation between x and ln y for the ith component of the PDF.
    corr_xy_n = corr_xy * sqrt( exp( sigma_y_n**2.0 ) - 1.0 ) / sigma_y_n 

    return
  end function corr_NL2NN

  !=============================================================================
  pure function corr_LL2NN( corr_xy, sigma_x_n, sigma_y_n )  &
  result( corr_xy_n )

    ! Description:
    ! For lognormally-distributed variables x and y, this function finds the
    ! correlation between ln x and ln y (corr_xy_n) for the ith component of the
    ! PDF, given the correlation between x and y (corr_xy), the standard
    ! deviation of ln x (sigma_x_n), and the standard deviation of ln y
    ! (sigma_y_n) for the ith component of the PDF.  The value of ln x (or ln y)
    ! is distributed normally when x (or y) is distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. C-3.
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) ::  &
      corr_xy,   & ! Correlation between x and y (ith PDF component)  [-]
      sigma_x_n, & ! Standard deviation of ln x (ith PDF component)   [-]
      sigma_y_n    ! Standard deviation of ln y (ith PDF component)   [-]

    ! Return Variable
    real ::  &
      corr_xy_n  ! Correlation between ln x and ln y (ith PDF component)  [-]

    ! Find the correlation between ln x and ln y for the ith component of the
    ! PDF.
    corr_xy_n = log( 1.0 + corr_xy * sqrt( exp( sigma_x_n**2.0 ) - 1.0 )  &
                                   * sqrt( exp( sigma_y_n**2.0 ) - 1.0 )  )  &
                / ( sigma_x_n * sigma_y_n )

    return
  end function corr_LL2NN

  !=============================================================================
  recursive function factorial( num )  &
  result( fact )

    ! Description:
    ! Calculates the factorial of an integer (num).
    !
    ! Note:  Performing this operation on an integer data type means that
    !        overflow occurs at any integer higher than 12.
    !        In the future, this function may be better replaced by a simple
    !        call to the gamma function.

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      num  ! Integer of which to take the factorial.

    ! Return Variable
    integer ::  &
      fact  ! Factorial of num.

    if ( num == 0 ) then
       fact = 1
    else
       fact = num * factorial( num - 1 )
    endif

    return
  end function factorial

  !=============================================================================
  double precision function Dv_fnc( order, argument )

    ! Description:
    ! Compute the parabolic cylinder function in terms of Dv
    ! using an Algorithm from ACM TOMS.  Replaces the more expensive
    ! D_fnc function used previously.

    ! References:
    !  Algorithm 850, collected algorithms from ACM.
    !  ACM Transactions on Mathematical Software,
    !  Vol. 32, No. 1, March 2006 pp. 102--112

    !-----------------------------------------------------------------------

    use parabolic, only:  & 
        gamma,  & ! Procedure(s) 
        parab

    use constants_clubb, only:  & 
        pi_dp ! Variable(s)

    implicit none

    ! External
    intrinsic :: sin

    ! Parameter constants
    integer, parameter :: & 
      scaling = 0 ! 0 = Unscaled functions, 1 = scaled functions

    double precision, parameter :: limit = 10.0d0**308

    ! Input Variables
    double precision, intent(in) :: & 
      order,    & ! Order 'a' of Dv(a,x)         [-]
      argument    ! Argument 'x' of Dv(a,x)      [-]

    ! Local Variables
    double precision, dimension(2) :: & 
      uaxx, & ! U(a,x), U'(a,x)                [-]
      vaxx    ! V(a,x), V'(a,x)                [-]
    ! Where a is the order and x is the argument

    integer :: ierr ! Error condition

    if ( argument <= 0.0d0 ) then
      call parab( -order-0.5, -argument, scaling, uaxx, vaxx, ierr )
      Dv_fnc = vaxx(1) / ( (1.0d0/pi_dp) * gamma( -order ) ) & 
             - sin( pi_dp * ( -order-0.5 ) ) * uaxx(1)
    else
      call parab( -order-0.5, argument, scaling, uaxx, vaxx, ierr )
      Dv_fnc = uaxx(1)
    end if

    ! Handle the overflow condition
    if ( ierr /= 0 ) then
      Dv_fnc = limit
    end if

    return
  end function Dv_fnc

!===============================================================================

end module KK_utilities

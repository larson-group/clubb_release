! $Id$
!===============================================================================
module KK_utilities

  implicit none

  private ! Set default scope to private

  public :: mean_L2N,      &
            mean_L2N_dp,   &
            stdev_L2N,     &
            stdev_L2N_dp,  &
            corr_NL2NN,    &
            corr_NL2NN_dp, &
            corr_LL2NN,    &
            corr_LL2NN_dp, &
            factorial,     &
            Dv_fnc,        & ! Parabolic Cylinder Function, D.
            calc_corr_sx,  &
            G_T_p

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

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      mu_x,        & ! Mean of x (ith PDF component)       [-]
      sigma_sqd_x    ! Variance of x (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      mu_x_n  ! Mean of ln x (ith PDF component)           [-]

    ! Find the mean of ln x for the ith component of the PDF.
    mu_x_n = log( mu_x / sqrt( 1.0 + sigma_sqd_x / mu_x**2 ) )

    return
  end function mean_L2N

  !=============================================================================
  pure function mean_L2N_dp( mu_x, sigma_sqd_x )  &
  result( mu_x_n )
  
    ! Description:
    ! For a lognormally-distributed variable x, this function finds the mean of
    ! ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
    ! and the variance of x (sigma_sqd_x) for the ith component of the PDF.
    ! The value ln x is distributed normally when x is distributed lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) ::  &
      mu_x,        & ! Mean of x (ith PDF component)       [-]
      sigma_sqd_x    ! Variance of x (ith PDF component)   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      mu_x_n  ! Mean of ln x (ith PDF component)           [-]

    ! Find the mean of ln x for the ith component of the PDF.
    mu_x_n = log( mu_x / sqrt( 1.0_dp + sigma_sqd_x / mu_x**2 ) )

    return
  end function mean_L2N_dp

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

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      mu_x,        & ! Mean of x (ith PDF component)       [-]
      sigma_sqd_x    ! Variance of x (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]

    ! Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = sqrt( log( 1.0 + sigma_sqd_x / mu_x**2 ) )

    return
  end function stdev_L2N

  !=============================================================================
  pure function stdev_L2N_dp( mu_x, sigma_sqd_x )  &
  result( sigma_x_n )

    ! Description:
    ! For a lognormally-distributed variable x, this function finds the standard
    ! deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
    ! mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
    ! of the PDF.  The value ln x is distributed normally when x is distributed
    ! lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) ::  &
      mu_x,        & ! Mean of x (ith PDF component)       [-]
      sigma_sqd_x    ! Variance of x (ith PDF component)   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]

    ! Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = sqrt( log( 1.0_dp + sigma_sqd_x / mu_x**2 ) )

    return
  end function stdev_L2N_dp

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

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_xy,   & ! Correlation between x and y (ith PDF component)  [-]
      sigma_y_n    ! Standard deviation of ln y (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      corr_xy_n  ! Correlation between x and ln y (ith PDF component) [-]

    ! Find the correlation between x and ln y for the ith component of the PDF.
    corr_xy_n = corr_xy * sqrt( exp( sigma_y_n**2 ) - 1.0 ) / sigma_y_n 

    return
  end function corr_NL2NN

  !=============================================================================
  pure function corr_NL2NN_dp( corr_xy, sigma_y_n )  &
  result( corr_xy_n )

    ! Description:
    ! For a normally-distributed variable x and a lognormally-distributed
    ! variable y, this function finds the correlation between x and ln y
    ! (corr_xy_n) for the ith component of the PDF, given the correlation
    ! between x and y (corr_xy) and the standard deviation of ln y (sigma_y_n)
    ! for the ith component of the PDF.  The value ln y is distributed normally
    ! when y is distributed lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. B-1.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      corr_xy,   & ! Correlation between x and y (ith PDF component)  [-]
      sigma_y_n    ! Standard deviation of ln y (ith PDF component)   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      corr_xy_n  ! Correlation between x and ln y (ith PDF component) [-]

    ! Find the correlation between x and ln y for the ith component of the PDF.
    corr_xy_n = corr_xy * sqrt( exp( sigma_y_n**2 ) - 1.0_dp ) / sigma_y_n 

    return
  end function corr_NL2NN_dp

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

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      corr_xy,   & ! Correlation between x and y (ith PDF component)  [-]
      sigma_x_n, & ! Standard deviation of ln x (ith PDF component)   [-]
      sigma_y_n    ! Standard deviation of ln y (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      corr_xy_n  ! Correlation between ln x and ln y (ith PDF component)  [-]

    ! Find the correlation between ln x and ln y for the ith component of the
    ! PDF.
    corr_xy_n = log( 1.0 + corr_xy * sqrt( exp( sigma_x_n**2 ) - 1.0 )  &
                                   * sqrt( exp( sigma_y_n**2 ) - 1.0 )  )  &
                / ( sigma_x_n * sigma_y_n )

    return
  end function corr_LL2NN

  !=============================================================================
  pure function corr_LL2NN_dp( corr_xy, sigma_x_n, sigma_y_n )  &
  result( corr_xy_n )

    ! Description:
    ! For lognormally-distributed variables x and y, this function finds the
    ! correlation between ln x and ln y (corr_xy_n) for the ith component of the
    ! PDF, given the correlation between x and y (corr_xy), the standard
    ! deviation of ln x (sigma_x_n), and the standard deviation of ln y
    ! (sigma_y_n) for the ith component of the PDF.  The value of ln x (or ln y)
    ! is distributed normally when x (or y) is distributed lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. C-3.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) ::  &
      corr_xy,   & ! Correlation between x and y (ith PDF component)  [-]
      sigma_x_n, & ! Standard deviation of ln x (ith PDF component)   [-]
      sigma_y_n    ! Standard deviation of ln y (ith PDF component)   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      corr_xy_n  ! Correlation between ln x and ln y (ith PDF component)  [-]

    ! Find the correlation between ln x and ln y for the ith component of the
    ! PDF.
    corr_xy_n = log( 1.0_dp + corr_xy * sqrt( exp( sigma_x_n**2 ) - 1.0_dp )  &
                                     * sqrt( exp( sigma_y_n**2 ) - 1.0_dp )  )  &
                / ( sigma_x_n * sigma_y_n )

    return
  end function corr_LL2NN_dp

!  !=============================================================================
!  recursive function factorial( num )  &
!  result( fact )
!
!    ! Description:
!    ! Calculates the factorial of an integer (num).
!    !
!    ! Note:  Performing this operation on an integer data type means that
!    !        overflow occurs at any integer higher than 12.
!    !        In the future, this function may be better replaced by a simple
!    !        call to the gamma function.
!
!    ! References:
!    !-----------------------------------------------------------------------
!
!    implicit none
!
!    ! Input Variable
!    integer, intent(in) :: &
!      num  ! Integer of which to take the factorial.
!
!    ! Return Variable
!    integer ::  &
!      fact  ! Factorial of num.
!
!    if ( num == 0 ) then
!       fact = 1
!    else
!       fact = num * factorial( num - 1 )
!    endif
!
!    return
!  end function factorial
!
  !=============================================================================
  function factorial( num )

    ! Description:
    ! Calculates the factorial of an integer (num).

    ! References:
    !-----------------------------------------------------------------------

    use parabolic, only:  &
        gamma  ! Procedure(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      num  ! Integer of which to take the factorial.

    ! Return Variable
    real( kind = dp ) :: &
      factorial  ! Factorial of num.


    ! Calculate the factorial of integer num using the gamma function.
    factorial = gamma( dble( num+1 ) )


    return

  end function factorial

  !=============================================================================
  function Dv_fnc( order, argument )

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

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! External
    intrinsic :: sin

    ! Parameter constants
    integer, parameter :: & 
      scaling = 0 ! 0 = Unscaled functions, 1 = scaled functions

    real( kind = dp ), parameter :: limit = 10.0_dp**308

    ! Input Variables
    real( kind = dp ), intent(in) :: & 
      order,    & ! Order 'a' of Dv(a,x)         [-]
      argument    ! Argument 'x' of Dv(a,x)      [-]

    ! Return Variable
    real( kind = dp ) :: &
      Dv_fnc

    ! Local Variables
    real( kind = dp ), dimension(2) :: & 
      uaxx, & ! U(a,x), U'(a,x)                [-]
      vaxx    ! V(a,x), V'(a,x)                [-]
    ! Where a is the order and x is the argument

    integer :: ierr ! Error condition

    if ( argument <= 0.0_dp ) then
      call parab( -order-0.5_dp, -argument, scaling, uaxx, vaxx, ierr )
      Dv_fnc = vaxx(1) / ( (1.0_dp/pi_dp) * gamma( -order ) ) & 
             - sin( pi_dp * ( -order-0.5_dp ) ) * uaxx(1)
    else
      call parab( -order-0.5_dp, argument, scaling, uaxx, vaxx, ierr )
      Dv_fnc = uaxx(1)
    end if

    ! Handle the overflow condition
    if ( ierr /= 0 ) then
      Dv_fnc = limit
    end if

    return
  end function Dv_fnc

  !=============================================================================
  pure function calc_corr_sx( crt_i, cthl_i, sigma_rt_i, sigma_thl_i,  &
                              sigma_s_i, corr_rtx_i, corr_thlx_i )  &
  result( corr_sx_i )

    ! Description:
    ! This function calculates the correlation between extended liquid water
    ! mixing ratio, s, and a generic variable x, within the ith component of the
    ! PDF.  The variable s can be split into mean and turbulent components, such
    ! that:
    !
    ! s = <s> + s';
    !
    ! where < > denotes a mean field an ' denotes a turbulent component.
    !
    ! The linearized equation for s' is given in Larson et al. (2001), where
    ! within the ith component of the PDF:
    !
    ! s_(i)' = Coef_rt(i) * r_t(i)' - Coef_thl(i) * th_l(i)'.
    !
    ! The equation for s' can be multiplied by x'.  The equation becomes:
    !
    ! s'x'_(i) = Coef_rt(i) * r_t'x'_(i) - Coef_thl(i) * th_l'x'_(i).
    !
    ! Averaging both sides, the covariance <s'x'> is given by the equation:
    !
    ! <s'x'_(i)> = Coef_rt(i) * <r_t'x'_(i)> - Coef_thl(i) * <th_l'x'_(i)>.
    !
    ! This equation can be rewritten as:
    !
    ! sigma_s(i) * sigma_x(i) * corr_sx(i)
    !   = Coef_rt(i) * sigma_rt(i) * sigma_x(i) * corr_rtx(i)
    !     - Coef_thl(i) * sigma_thl(i) * sigma_x(i) * corr_thlx(i).
    !
    ! This equation can be solved for corr_sx(i):
    !
    ! corr_sx(i) = Coef_rt(i) * ( sigma_rt(i) / sigma_s(i) ) * corr_rtx(i)
    !              - Coef_thl(i) * ( sigma_thl(i) / sigma_s(i) ) * corr_thlx(i).
    !
    ! The correlation between s and x within the ith component of the PDF is
    ! calculated.

    ! References:
    !  Larson, V. E., R. Wood, P. R. Field, J.-C. Golaz, T. H. Vonder Haar,
    !    W. R. Cotton, 2001: Systematic Biases in the Microphysics and
    !    Thermodynamics of Numerical Models That Ignore Subgrid-Scale
    !    Variability. J. Atmos. Sci., 58, 1117--1128.
    !  -- Eq. 13 and 14.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      crt_i,       & ! Coefficient of r_t for s' (ith PDF component)         [-]
      cthl_i,      & ! Coefficient of th_l for s' (ith PDF component)      [1/K]
      sigma_rt_i,  & ! Standard deviation of r_t (ith PDF component)     [kg/kg]
      sigma_thl_i, & ! Standard deviation of th_l (ith PDF component)        [K]
      sigma_s_i,   & ! Standard deviation of s (ith PDF component)       [kg/kg]
      corr_rtx_i,  & ! Correlation between r_t and x (ith PDF component)     [-]
      corr_thlx_i    ! Correlation between th_l and x (ith PDF component)    [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_sx_i  ! Correlation of s and x (ith PDF component)   [-]


    ! Calculate the correlation of s and x in the ith PDF component.
    if ( sigma_s_i > 0.0 ) then

       corr_sx_i = crt_i * ( sigma_rt_i / sigma_s_i ) * corr_rtx_i  &
                   - cthl_i * ( sigma_thl_i / sigma_s_i ) * corr_thlx_i

    else  ! sigma_s_i = 0

       ! The variance of s_(i) is 0.  This means that s is constant within the
       ! ith PDF component and covariance <s'x'_(i)> is also 0.  The correlation
       ! between s and x is 0 in the ith PDF component.
       corr_sx_i = 0.0

    endif


    return

  end function calc_corr_sx

  !=============================================================================
  function G_T_p( T_in_K, p_in_Pa )

    ! Description:
    ! Calculate G(T,p) from the drop radius growth equation.  This is used as a
    ! coefficient in the KK evaporation equation.  The equation for G(T,p) is
    ! taken from Rogers and Yau (1989).

    ! References:
    !  Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    !    Parameterization in a Large-Eddy Simulation Model of Marine
    !    Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
    !  -- Eq. 22.
    !  Rogers, R. R. and M. K. Yau, 1989:  A Short Course in Cloud Physics. 3rd
    !    edition, Butterworth-Heinemann, 290 pp.
    !  -- Eq. 7.17 and 7.18.
    !-----------------------------------------------------------------------

    use saturation, only: & 
        sat_vapor_press_liq ! Procedure(s)

    use constants_clubb, only: & 
        T_freeze_K, & ! Constant(s)
        ep,         &
        rho_lw,     & 
        Lv,         & 
        Rv

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      T_in_K,  & ! Temperature  [K]
      p_in_Pa    ! Pressure     [Pa]

    ! Return Variable
    real( kind = core_rknd ) :: &
      G_T_p  ! Function G(T,p) as found in the evaporation equation  [m^2/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      Ka,      & ! Coefficient of thermal conductivity of air          [J/(msK)]
      Dv,      & ! Coefficient of diffusion of water vapor in air      [m^2/s]
      Fk,      & ! Term in denominator associated with heat conduction [s/m^2]
      Fd,      & ! Term in denominator associated with vapor diffusion [s/m^2]
      esatv,   & ! Saturation vapor pressure                           [Pa]
      Celsius    ! Temperature in Celsius                              [deg C]


    ! Temperature in degrees Celsius.
    Celsius = T_in_K - T_freeze_K

    ! Coefficient of thermal conductivity of air.
    Ka = ( 5.69_core_rknd + 0.017_core_rknd * Celsius )  &
         * 0.00001_core_rknd                       ! Ka in cal./(cm.*sec.*C)
    Ka = 4.1868_core_rknd * 100.0_core_rknd * Ka   ! Ka in J./(m.*sec.*K)

    ! Coefficient of diffusion of water vapor in air.
    Dv = 0.221_core_rknd * ( (T_in_K/273.16_core_rknd)**1.94_core_rknd )  &
         * ( 101325.0_core_rknd / p_in_Pa )   ! Dv in (cm.^2)/sec.
    Dv = Dv / 10000.0_core_rknd               ! Dv in (m.^2)/sec.

    ! Calculate saturation mixing ratio and saturation vapor pressure.
    esatv = sat_vapor_press_liq( T_in_K )

    ! The values of F_k and F_d are found in Rogers and Yau (1989);
    ! Eq. 7.17 and 7.18.
    Fk = ( Lv / ( Rv * T_in_K ) - 1.0 ) * ( Lv * rho_lw ) / ( Ka * T_in_K )
    Fd = ( rho_lw * Rv * T_in_K ) / ( Dv * esatv )

    ! Calculate G(T,p).
    G_T_p = 1.0 / (Fk + Fd)


    return

  end function G_T_p

!===============================================================================

end module KK_utilities

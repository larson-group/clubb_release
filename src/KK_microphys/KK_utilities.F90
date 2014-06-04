! $Id$
!===============================================================================
module KK_utilities

  implicit none

  private ! Set default scope to private

  public :: factorial,           &
            Dv_fnc,              & ! Parabolic Cylinder Function, D.
            get_cloud_top_level, &
            G_T_p

  contains

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
    factorial = gamma( real( num+1, kind=dp ) )


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
        pi_dp,       & ! Constant(s)
        one_dp,      &
        one_half_dp, &
        zero_dp

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

    if ( argument <= zero_dp ) then
      call parab( -order-one_half_dp, -argument, scaling, uaxx, vaxx, ierr )
      Dv_fnc = vaxx(1) / ( (one_dp/pi_dp) * gamma( -order ) ) & 
             - sin( pi_dp * ( -order-one_half_dp ) ) * uaxx(1)
    else
      call parab( -order-one_half_dp, argument, scaling, uaxx, vaxx, ierr )
      Dv_fnc = uaxx(1)
    end if

    ! Handle the overflow condition
    if ( ierr /= 0 ) then
      Dv_fnc = limit
    end if

    return
  end function Dv_fnc

  !=============================================================================
  function get_cloud_top_level( nz, mean_rc ) &
  result( cloud_top_level )

    ! Description:
    ! Find cloud top at a given model time step.  This function can be used to
    ! find overall cloud top or the cloud top within a given PDF component.
    ! This function finds cloud top by looping downward from the top of the
    ! model and returning the index of the first vertical level that has a mean
    ! cloud water mixing ratio greater than the tolerance amount.  In a scenario
    ! that there is not any cloud found, the function returns a value of 1 (for
    ! vertical level 1, which is below the model surface).

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        rc_tol    ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      mean_rc    ! Mean cloud water mixing ratio                [kg/kg]
 
    ! Return Variable
    integer :: &
      cloud_top_level    ! Vertical level index of cloud top

    ! Local Variable
    integer :: k    ! Vertical level index


    ! Start at the model upper boundary and loop downwards until cloud top is
    ! found or the model lower boundary is reached.
    k = nz
    do
       if ( mean_rc(k) > rc_tol ) then
          ! A level with mean cloud water mixing ratio greater than the
          ! tolerance amount has been found.  Cloud top has been found.
          cloud_top_level = k
          exit
       elseif ( k == 1 ) then
          ! There was not any cloud found in the model vertical domain.
          ! Return a value of 1.
          cloud_top_level = 1
          exit
       else
          ! Continue down another vertical level.
          k = k - 1
       endif
    enddo


    return

  end function get_cloud_top_level

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

    use constants_clubb, only: &
        one_hundred, & ! Constant(s)
        one

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
    Ka = 4.1868_core_rknd * one_hundred * Ka       ! Ka in J./(m.*sec.*K)

    ! Coefficient of diffusion of water vapor in air.
    Dv = 0.221_core_rknd * ( (T_in_K/273.16_core_rknd)**1.94_core_rknd )  &
         * ( 101325.0_core_rknd / p_in_Pa )   ! Dv in (cm.^2)/sec.
    Dv = Dv / 10000.0_core_rknd               ! Dv in (m.^2)/sec.

    ! Calculate saturation mixing ratio and saturation vapor pressure.
    esatv = sat_vapor_press_liq( T_in_K )

    ! The values of F_k and F_d are found in Rogers and Yau (1989);
    ! Eq. 7.17 and 7.18.
    Fk = ( Lv / ( Rv * T_in_K ) - one ) * ( Lv * rho_lw ) / ( Ka * T_in_K )
    Fd = ( rho_lw * Rv * T_in_K ) / ( Dv * esatv )

    ! Calculate G(T,p).
    G_T_p = one / (Fk + Fd)


    return

  end function G_T_p

!===============================================================================

end module KK_utilities
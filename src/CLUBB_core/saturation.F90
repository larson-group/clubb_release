!$Id$

!-----------------------------------------------------------------------
module saturation
!       Description: 
!         Contains functions that compute saturation with respect
!         to liquid or ice.
!-----------------------------------------------------------------------
  implicit none

  private  ! Change default so all items private

  public   :: sat_mixrat_liq, sat_mixrat_ice, sat_rcm

  private  :: sat_vapor_press_liq, sat_vapor_press_ice

contains

!-------------------------------------------------------------------------
  elemental real function sat_mixrat_liq( p_in_Pa, T_in_K )

!       Description:
!       Used to compute the saturation mixing ratio of liquid water.

!       References:
!       Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

    use constants, only: & 
        ep ! Variable

    implicit none

    ! Input Variables
    real, intent(in) ::  & 
    p_in_Pa,  & ! Pressure    [Pa]
    T_in_K      ! Temperature [K]

    ! Local Variables
    real :: esatv

    ! Saturation Vapor Pressure, esat, can be found to be approximated
    ! in many different ways.

    esatv = sat_vapor_press_liq( T_in_K )

    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v
    sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - esatv ) )

    return

  end function sat_mixrat_liq

!------------------------------------------------------------------------
  pure function sat_vapor_press_liq( T_in_K ) result ( esat )

!       Description:
!       Computes SVP for water vapor.

!       References:
!       ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!         and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!         pp. 1507--1513
!------------------------------------------------------------------------

    use constants, only: T_freeze_K

    implicit none

    ! External
    intrinsic :: exp

    ! Parameter Constants
    logical, parameter :: l_Flatau = .true.

    ! Relative error norm expansion (-50 to 50 deg_C) from
    ! Table 3 of pp. 1510 of Flatau et al. 1992 (Water Vapor)
    ! (The 100 coefficient converts from mb to Pa)
    real, dimension(7), parameter :: a = & 
    100.* (/ 6.11176750,      0.443986062,     0.143053301E-01, & 
             0.265027242E-03, 0.302246994E-05, 0.203886313E-07, & 
             0.638780966E-10 /)

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esat  ! Saturation vapor pressure over water [Pa]

    ! Local Variables
    real :: T_in_C
!   integer :: i ! Loop index

    if ( l_Flatau ) then
      ! Determine deg K - 273.15
      T_in_C = T_in_K - T_freeze_K 

      ! Polynomial approx. (Flatau, et al. 1992)

     ! Formulation 1 
     ! This is the most generalized and the least efficient.
     ! Based on Wexler's expressions(2.1)-(2.4) (See Flatau et al. p 1508)
     ! e_{sat} = a_1 + a_2 ( T - T_0 ) + ... + a_{n+1} ( T - T_0 )^n

!     esat = a(1)

!     do i = 2, size( a ) , 1
!       esat = esat + a(i) * ( T_in_C )**(i-1)
!     end do

     ! Formulation 2 (no loop)
     ! This would be efficient if GNU's libm pow function wasn't slow
!    esat = 100.*( a(1) + a(2) * T_in_C + a(3) * T_in_C**2 + a(4) * T_in_C**3 &
!       + a(5) * T_in_C**4 + a(6) * T_in_C**5 + a(7) * T_in_C**6 )

     ! Formulation 3 (no pow, no loop)
     esat = ( a(1) + T_in_C *( a(2) + T_in_C *( a(3) + T_in_C *( a(4) + T_in_C &
        *( a(5) + T_in_C*( a(6) + T_in_C*( a(7) ) ) ) ) ) ) )


    else
      ! (Bolton 1980) approx.
      ! Generally more expensive, but not on Sun Fortran with -xlibmopt
      esat = 611.2 * exp( (17.67*(T_in_K-T_freeze_K)) / (T_in_K-29.65) )
    end if

    return

  end function sat_vapor_press_liq

!------------------------------------------------------------------------
  real function sat_mixrat_ice( p_in_Pa, T_in_K )

!       Description:
!       Used to compute the saturation mixing ratio of ice. 

!       References:
!       Formula from Emanuel 1994, 4.4.15
!-------------------------------------------------------------------------

    use constants, only: & 
        ep ! Variable(s)

    implicit none

    ! Input Variables

    real, intent(in) :: &
    p_in_Pa, &          ! Pressure [Pa]
    T_in_K              ! Temperature [K]

    ! Local Variables

    real :: esat_ice

    ! Compute SVP for ice

    esat_ice = sat_vapor_press_ice( T_in_K )

    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v

    sat_mixrat_ice = ep * ( esat_ice / ( p_in_Pa - esat_ice ) )

    return

  end function sat_mixrat_ice

!------------------------------------------------------------------------
  real pure function sat_vapor_press_ice( T_in_K ) result ( esati )
!
!       Description:
!       Computes SVP for ice.
!
!       References:
!       ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!         and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!         pp. 1507--1513
!------------------------------------------------------------------------
    use constants, only: T_freeze_K


    implicit none

    ! External
    intrinsic :: exp, log

    ! Parameter Constants
    logical, parameter :: l_Flatau = .false.

    ! Relative error norm expansion (-50 to 0 deg_C) from
    ! Table 3 of pp. 1510 of Flatau et al. 1992 (Ice)
    real, dimension(7), parameter :: a = & 
    100. * (/ 6.10952665,      0.501948366,     0.18628899E-01, & 
              0.403488906E-03, 0.539797852E-05, 0.420713632E-07, & 
              0.147271071E-09 /)

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Local Variables
    integer :: i

    if ( l_Flatau ) then
      ! Polynomial approx. (Flatau, et al. 1992)
      esati = a(1)

      do i = 2, size( a ), 1
        esati = esati + a(i) * ( T_in_K-T_freeze_K )**(i-1)
      end do

    else
      ! Exponential approx. (Bolton?)
      esati = 100.0 * exp( 23.33086 - (6111.72784/T_in_K) + (0.15215*log( T_in_K )) )
    end if

    return

  end function sat_vapor_press_ice

!-------------------------------------------------------------------------
  FUNCTION sat_rcm( thlm, rtm, p_in_Pa, exner )

    ! Description:
    !
    ! This function uses an iterative method to find the value of rcm 
    ! from an initial profile that has saturation at some point.
    !-------------------------------------------------------------------------

    USE constants, only: & 
        Cp,            & ! Variable(s)
        Lv,            &
        zero_threshold

    implicit none
    
    ! Input Variable(s)
    REAL, INTENT(IN):: thlm         ! Liquid Water Potential Temperature [K]
    REAL, INTENT(IN):: rtm          ! Total Water Mixing Ratio       [kg/kg]
    REAL, INTENT(IN):: p_in_Pa      ! Pressure                          [Pa]
    REAL, INTENT(IN):: exner        ! Exner function                     [-]

    REAL:: sat_rcm

    REAL:: theta
    REAL:: answer, too_low, too_high

    INTEGER:: iteration

    REAL, PARAMETER:: tolerance = 0.001

    ! Default initialization
    theta = thlm
    iteration = 0
    too_high = 0.0
    too_low = 0.0

    DO

      iteration = iteration + 1

      answer = &
      theta - (Lv/(Cp*exner)) &
             *(MAX( rtm - sat_mixrat_liq(p_in_Pa,theta*exner), zero_threshold ))

      IF ( ABS(answer - thlm) <= tolerance ) THEN
        EXIT
      ELSEIF ( answer - thlm > tolerance ) THEN
        too_high = theta
      ELSEIF ( thlm - answer > tolerance ) THEN
        too_low = theta
      ENDIF

      ! For the first timestep, be sure to set a "too_high"
      ! that is "way too high."
      IF ( iteration == 1 ) THEN
        too_high = theta + 20.0
      ENDIF

      theta = (too_low + too_high)/2.0

    ENDDO

    sat_rcm = MAX( rtm - sat_mixrat_liq( p_in_Pa, theta*exner), zero_threshold )

  END FUNCTION sat_rcm
  
end module saturation

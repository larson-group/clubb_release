!$Id: saturation.F90,v 1.7 2008-08-08 15:13:18 faschinj Exp $

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
    logical, parameter :: lFlatau = .false.

    ! Relative error norm expansion (-50 to 50 deg_C) from
    ! Table 3 of pp. 1511 of Flatau et al. 1992 (Water Vapor)
    real, dimension(7), parameter :: a = & 
    (/ 6.11176750,      0.443986062,     0.143053301E-01, & 
       0.265027242E-03, 0.302246994E-05, 0.203886313E-07, & 
       0.638780966E-10 /)

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esat

    ! Local Variables
    integer :: i

    if ( lFlatau ) then
      ! Polynomial approx. (Flatau, et al. 1992)
      esat = a(1)

      do i = 2, 7, 1
        esat = esat + a(i) * ( T_in_K-T_freeze_K )**(i-1)
      end do

      esat = 100.0 * esat ! Convert units

    else
      ! (Bolton 1980) approx.
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
    logical, parameter :: lFlatau = .false.

    ! Relative error norm expansion (-50 to 0 deg_C) from
    ! Table 3 of pp. 1511 of Flatau et al. 1992 (Ice)
    real, dimension(7), parameter :: a = & 
    (/ 6.10952665,      0.501948366,     0.18628899E-01, & 
       0.403488906E-03, 0.539797852E-05, 0.420713632E-07, & 
       0.147271071E-09 /)

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Local Variables
    integer :: i

    if ( lFlatau ) then
      ! Polynomial approx. (Flatau, et al. 1992)
      esati = a(1)

      do i = 2, 7, 1
        esati = esati + a(i) * ( T_in_K-T_freeze_K )**(i-1)
      end do

      esati = 100.0 * esati ! Convert units

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
        Cp,  & ! Variable(s)
        Lv

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

      answer = theta - (Lv/(Cp*exner))*(MAX( rtm - sat_mixrat_liq(p_in_Pa,theta*exner), 0.0 ))

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

    sat_rcm = MAX( rtm - sat_mixrat_liq( p_in_Pa, theta*exner), 0.0 )

  END FUNCTION sat_rcm
  
end module saturation

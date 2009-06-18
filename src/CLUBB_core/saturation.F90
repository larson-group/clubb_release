!$Id$
!-----------------------------------------------------------------------
module saturation

! Description:
!   Contains functions that compute saturation with respect
!   to liquid or ice.
!-----------------------------------------------------------------------
  implicit none

  private  ! Change default so all items private

  public   :: sat_mixrat_liq, sat_mixrat_ice, sat_rcm

  private  :: sat_vapor_press_liq_flatau, sat_vapor_press_liq_bolton
  private  :: sat_vapor_press_ice_flatau, sat_vapor_press_ice_bolton

  contains

!-------------------------------------------------------------------------
  elemental real function sat_mixrat_liq( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of liquid water.

! References:
!   Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

    use constants, only: & 
      ep, & ! Variable
      fstderr

    use model_flags, only: &
      saturation_formula ! Variable

    implicit none

    ! Input Variables
    real, intent(in) ::  & 
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

    ! Local Variables
    real :: esatv

    ! --- Begin Code ---

    ! Undefined approximation
    esatv = -99999.999

    ! Saturation Vapor Pressure, esat, can be found to be approximated
    ! in many different ways.

    !**************************************************************************
    !                         ***** IMPORTANT *****
    ! If new saturation_formula's are added with a different string length 
    ! then this select case must be modified accordingly:
    ! 
    ! (1) If you add to the length of saturation_formula, spaces must be 
    ! added at the end of current cases. E.g. case ( "bolton  " ) for len=8. 
    ! 
    ! (2) If you use as name for a saturation_formula with less than 6 
    ! characters, you must add spaces at the end of the name.  
    ! E.g. case ( "mine  " ) for len=6.
    !                         ***** IMPORTANT *****
    !**************************************************************************
    select case ( saturation_formula )
    case ( "bolton", "Bolton" )
      ! Using the Bolton 1980 approximations for SVP over vapor
      esatv = sat_vapor_press_liq_bolton( T_in_K )

    case ( "flatau", "Flatau" )
      ! Using the Flatau, et al. polynomial approximation for SVP over vapor
      esatv = sat_vapor_press_liq_flatau( T_in_K )

    ! Add new cases after this

    end select

    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v
    sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - esatv ) )

    return
  end function sat_mixrat_liq

!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_flatau( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor.

! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------

    use constants, only: T_freeze_K

    implicit none

    ! Relative error norm expansion (-50 to 50 deg_C) from
    ! Table 3 of pp. 1510 of Flatau et al. 1992 (Water Vapor)
    ! (The 100 coefficient converts from mb to Pa)
!   real, dimension(7), parameter :: a = & 
!   100.* (/ 6.11176750,      0.443986062,     0.143053301E-01, & 
!            0.265027242E-03, 0.302246994E-05, 0.203886313E-07, & 
!            0.638780966E-10 /)
    ! Relative error norm expansion (-85 to 70 deg_C) from
    ! Table 4 of pp. 1511 of Flatau et al.
    real, dimension(9), parameter :: a = & 
    100.* (/ 6.11583699,      0.444606896,     0.143177157E-01, & 
             0.264224321E-03, 0.299291081E-05, 0.203154182E-07, & 
             0.702620698E-10, 0.379534310E-13,-0.321582393E-15 /)

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esat  ! Saturation vapor pressure over water [Pa]

    ! Local Variables
    real :: T_in_C
!   integer :: i ! Loop index

    ! Determine deg K - 273.15
    T_in_C = T_in_K - T_freeze_K

    ! Polynomial approx. (Flatau, et al. 1992)

    ! This is the generalized formula but is not computationally efficient.
    ! Based on Wexler's expressions(2.1)-(2.4) (See Flatau et al. p 1508)
    ! e_{sat} = a_1 + a_2 ( T - T_0 ) + ... + a_{n+1} ( T - T_0 )^n

!   esat = a(1)

!   do i = 2, size( a ) , 1
!     esat = esat + a(i) * ( T_in_C )**(i-1)
!   end do

    ! The 8th order polynomial fit.  When running deep 
    ! convective cases I noted that absolute temperature often dips below
    ! -50 deg_C at higher altitudes, where the 6th order approximation is
    ! not accurate.  -dschanen 20 Nov 2008
    esat = a(1) + T_in_C*( a(2) + T_in_C*( a(3) + T_in_C*( a(4) + T_in_C &
    *( a(5) + T_in_C*( a(6) + T_in_C*( a(7) + T_in_C*( a(8) + T_in_C*( a(9) ) ) ) ) ) ) ) )

    return
  end function sat_vapor_press_liq_flatau


!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_bolton( T_in_K ) result ( esat )
! Description:
!   Computes SVP for water vapor.
! References:
!   Bolton 1980
!------------------------------------------------------------------------

    use constants, only: T_freeze_K

    implicit none

    ! External
    intrinsic :: exp

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esat  ! Saturation vapor pressure over water [Pa]

    ! (Bolton 1980) approx.
    ! Generally this more computationally expensive than the Flatau polnomial expansion
    esat = 611.2 * exp( (17.67*(T_in_K-T_freeze_K)) / (T_in_K-29.65) )

    return
  end function sat_vapor_press_liq_bolton

!------------------------------------------------------------------------
  elemental real function sat_mixrat_ice( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of ice.

! References:
!   Formula from Emanuel 1994, 4.4.15
!-------------------------------------------------------------------------

    use constants, only: & 
        ep ! Variable(s)
    use model_flags, only: &
      saturation_formula ! Variable(s)

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables

    real, intent(in) :: &
      p_in_Pa, &          ! Pressure [Pa]
      T_in_K              ! Temperature [K]

    ! Local Variables

    real :: esat_ice

    ! --- Begin Code ---

    ! Undefined approximation
    esat_ice = -99999.999

    select case ( trim( saturation_formula ) )
    case ( "bolton", "Bolton" )
      ! Using the Bolton 1980 approximations for SVP over ice
      esat_ice = sat_vapor_press_ice_bolton( T_in_K )

    case ( "flatau", "Flatau" )
      ! Using the Flatau, et al. polynomial approximation for SVP over ice
      esat_ice = sat_vapor_press_ice_flatau( T_in_K )

    ! Add new cases after this

    end select

    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v

    sat_mixrat_ice = ep * ( esat_ice / ( p_in_Pa - esat_ice ) )

    return

  end function sat_mixrat_ice

!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice_flatau( T_in_K ) result ( esati )
!
! Description:
!   Computes SVP for ice.
!
! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------
    use constants, only: T_freeze_K


    implicit none

    ! Relative error norm expansion (-90 to 0 deg_C) from
    ! Table 4 of pp. 1511 of Flatau et al. 1992 (Ice)
    real, dimension(9), parameter :: a = & 
    100. * (/ 6.09868993,      0.499320233,     0.184672631E-01, &
              0.402737184E-03, 0.565392987E-05, 0.521693933E-07, &
              0.307839583E-09, 0.105785160E-11, 0.161444444E-14 /)

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [deg_K]

    ! Output Variables
    real :: esati  ! Saturation vapor pressure over ice [Pa]

    ! Local Variables
    real :: T_in_C ! Temperature [deg_C]
!   integer :: i

    ! ---- Begin Code ----

    ! Determine deg K - 273.15
    T_in_C = T_in_K - T_freeze_K

    ! Polynomial approx. (Flatau, et al. 1992)
!   esati = a(1)

!   do i = 2, size( a ), 1
!     esati = esati + a(i) * ( T_in_C )**(i-1)
!   end do

    esati = a(1) + T_in_C*( a(2) + T_in_C*( a(3) + T_in_C*( a(4) + T_in_C &
    *( a(5) + T_in_C*( a(6) + T_in_C*( a(7) + T_in_C*( a(8) + T_in_C*( a(9) ) ) ) ) ) ) ) )

    return

  end function sat_vapor_press_ice_flatau

!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice_bolton( T_in_K ) result ( esati )
!
! Description:
!   Computes SVP for ice.
!
! References:
!   Bolton 1980
!------------------------------------------------------------------------
    use constants, only: T_freeze_K

    implicit none

    ! External
    intrinsic :: exp, log

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esati  ! Saturation vapor pressure over ice [Pa]

    ! Exponential approx.
    esati = 100.0 * exp( 23.33086 - (6111.72784/T_in_K) + (0.15215*log( T_in_K )) )

    return

  end function sat_vapor_press_ice_bolton

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

!$Id$
!-----------------------------------------------------------------------
module saturation

! Description:
!   Contains functions that compute saturation with respect
!   to liquid or ice.
!-----------------------------------------------------------------------

#ifdef GFDL
    use model_flags, only: &  ! h1g, 2010-06-18
       I_sat_sphum
#endif


  implicit none

  private  ! Change default so all items private

  public   :: sat_mixrat_liq, sat_mixrat_liq_lookup, sat_mixrat_ice, sat_rcm

  private  :: sat_vapor_press_liq_flatau, sat_vapor_press_liq_bolton
  private  :: sat_vapor_press_ice_flatau, sat_vapor_press_ice_bolton

  ! Lookup table of values for saturation 
  real, private, dimension(188:343) :: &
    svp_liq_lookup_table

  data svp_liq_lookup_table(188:343) / &
    0.049560547, 0.059753418, 0.070129395, 0.083618164, 0.09814453, &
    0.11444092, 0.13446045, 0.15686035, 0.18218994, 0.21240234, &
    0.24725342, 0.28668213, 0.33184814, 0.3826294, 0.4416504, &
    0.50775146, 0.58343506, 0.6694946, 0.7668457, 0.87750244, &
    1.0023804, 1.1434937, 1.3028564, 1.482544, 1.6847534, &
    1.9118042, 2.1671143, 2.4535522, 2.774231, 3.1330566, &
    3.5343628, 3.9819336, 4.480713, 5.036072, 5.6540527, &
    6.340088, 7.1015015, 7.9450684, 8.8793335, 9.91217, &
    11.053528, 12.313049, 13.70166, 15.231018, 16.91394, &
    18.764038, 20.795898, 23.025574, 25.470093, 28.147766, &
    31.078003, 34.282043, 37.782593, 41.60382, 45.771606, &
    50.31366, 55.259644, 60.641174, 66.492004, 72.84802, &
    79.74756, 87.23126, 95.34259, 104.12747, 113.634796, &
    123.91641, 135.02725, 147.02563, 159.97308, 173.93488, &
    188.97995, 205.18109, 222.61517, 241.36334, 261.51108, &
    283.14853, 306.37054, 331.27698, 357.97278, 386.56842, &
    417.17978, 449.9286, 484.94254, 522.3556, 562.30804, &
    604.947, 650.42645, 698.9074, 750.55835, 805.55554, &
    864.0828, 926.3325, 992.5052, 1062.8102, 1137.4657, &
    1216.6995, 1300.7483, 1389.8594, 1484.2896, 1584.3064, &
    1690.1881, 1802.224, 1920.7146, 2045.9724, 2178.3218, &
    2318.099, 2465.654, 2621.3489, 2785.5596, 2958.6758, &
    3141.101, 3333.2534, 3535.5657, 3748.4863, 3972.4792, &
    4208.024, 4455.616, 4715.7686, 4989.0127, 5275.8945, &
    5576.9795, 5892.8535, 6224.116, 6571.3926, 6935.3213, &
    7316.5674, 7715.8105, 8133.755, 8571.125, 9028.667, &
    9507.15, 10007.367, 10530.132, 11076.282, 11646.683, &
    12242.221, 12863.808, 13512.384, 14188.913, 14894.385, &
    15629.823, 16396.268, 17194.799, 18026.516, 18892.55, &
    19794.07, 20732.262, 21708.352, 22723.592, 23779.273, &
    24876.709, 26017.258, 27202.3, 28433.256, 29711.578, &
    31038.766 /

  contains

!-------------------------------------------------------------------------
  elemental real function sat_mixrat_liq( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of liquid water.

! References:
!   Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

    use constants_clubb, only: & 
      ep, & ! Variable
      fstderr

    implicit none

    ! Input Variables
    real, intent(in) ::  & 
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

   ! Local Variables
    real :: esatv

    ! --- Begin Code ---

    ! Calculate the SVP for water vapor.
    esatv = sat_vapor_press_liq( T_in_K )

#ifdef GFDL

    ! GFDL uses specific humidity
    ! Formula for Saturation Specific Humidity
     if( I_sat_sphum )  then   ! h1g, 2010-06-18 begin mod
           sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - (1.0-ep) * esatv ) )
     else
           sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - esatv ) )
     endif                     ! h1g, 2010-06-18 end mod
#else
    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v
    sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - esatv ) )
#endif

    return
  end function sat_mixrat_liq

!-------------------------------------------------------------------------
  elemental real function sat_mixrat_liq_lookup( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of liquid water.
!   This function utilizes sat_vapor_press_liq_lookup; the SVP is found
!   using a lookup table rather than calculating it using various
!   approximations.

! References:
!   Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

    use constants_clubb, only: & 
      ep, & ! Variable
      fstderr

    implicit none

    ! Input Variables
    real, intent(in) ::  & 
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

   ! Local Variables
    real :: esatv

    ! --- Begin Code ---

    ! Calculate the SVP for water vapor using a lookup table.
    esatv = sat_vapor_press_liq_lookup( T_in_K )

#ifdef GFDL

    ! GFDL uses specific humidity
    ! Formula for Saturation Specific Humidity
     if( I_sat_sphum )  then   ! h1g, 2010-06-18 begin mod
           sat_mixrat_liq_lookup = ep * ( esatv / ( p_in_Pa - (1.0-ep) * esatv ) )
     else
           sat_mixrat_liq_lookup = ep * ( esatv / ( p_in_Pa - esatv ) )
     endif                     ! h1g, 2010-06-18 end mod
#else
    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v
    sat_mixrat_liq_lookup = ep * ( esatv / ( p_in_Pa - esatv ) )
#endif

    return
  end function sat_mixrat_liq_lookup

!-----------------------------------------------------------------
  elemental function sat_vapor_press_liq( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor. Calls one of the other functions
!   that calculate an approximation to SVP.

! References:
!   None

    use model_flags, only: &
      saturation_formula, & ! Variable
      saturation_bolton, &
#ifdef GFDL
      saturation_gfdl, &
#endif
      saturation_flatau

    implicit none

    ! Input Variables
    real, intent(in) :: T_in_K     ! Temperature                          [K]

    ! Output Variables
    real :: esat      ! Saturation Vapor Pressure over Water [Pa]

    ! Undefined approximation
    esat = -99999.999

    ! Saturation Vapor Pressure, esat, can be found to be approximated
    ! in many different ways.
    select case ( saturation_formula )
    case ( saturation_bolton )
      ! Using the Bolton 1980 approximations for SVP over vapor
      esat = sat_vapor_press_liq_bolton( T_in_K )

    case ( saturation_flatau )
      ! Using the Flatau, et al. polynomial approximation for SVP over vapor
      esat = sat_vapor_press_liq_flatau( T_in_K )

    ! Add new cases after this
! ---> h1g
#ifdef GFDL
    case ( saturation_gfdl )
      ! Using GFDL polynomial approximation for SVP with respect to liquid
      esat = sat_vapor_press_liq_gfdl( T_in_K )
#endif
! <--- h1g

    end select

    return

  end function sat_vapor_press_liq

!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_lookup( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor, using a lookup table.
!
!   The lookup table was constructed using the Flatau approximation.

! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------

  implicit none

  ! Input Variables
  real, intent(in) :: T_in_K   ! Temperature   [K]

  ! Output Variables
  real :: esat  ! Saturation vapor pressure over water [Pa]

  ! Local Variables
  integer :: T_in_K_int

  T_in_K_int = int( anint( T_in_K ) )

  if ( T_in_K_int >= 188 .and. T_in_K_int <= 343 ) then
    ! Use the lookup table to determine the saturation vapor pressure.
    esat = svp_liq_lookup_table( T_in_K_int )
  else
    ! If we're outside the bounds of the lookup table, fall back to
    ! using a Flatau approximation.
    esat = sat_vapor_press_liq_flatau( T_in_K )
  end if

  end function sat_vapor_press_liq_lookup

!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_flatau( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor.

! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------

    use constants_clubb, only: T_freeze_K

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

    use constants_clubb, only: T_freeze_K

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


! ---> h1g, 2010-06-16
#ifdef GFDL
!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_gfdl( T_in_K ) result ( esat )
! Description:
! copy from "GFDL polysvp.F90" 
!  Compute saturation vapor pressure with respect to liquid  by using 
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T_in_K  is input in units of K.
!------------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esat  ! Saturation vapor pressure over water [Pa]

! Goff Gatch equation, uncertain below -70 C
      
         esat = 10.**(-7.90298*(373.16/T_in_K-1.)+ &
             5.02808*log10(373.16/T_in_K)- &
             1.3816e-7*(10.**(11.344*(1.-T_in_K/373.16))-1.)+ &
             8.1328e-3*(10.**(-3.49149*(373.16/T_in_K-1.))-1.)+ &
             log10(1013.246))*100.

    return
  end function sat_vapor_press_liq_gfdl
#endif
! <--- h1g, 2010-06-16

!------------------------------------------------------------------------
  elemental real function sat_mixrat_ice( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of ice.

! References:
!   Formula from Emanuel 1994, 4.4.15
!-------------------------------------------------------------------------

    use constants_clubb, only: & 
        ep ! Variable(s)

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

    ! Determine the SVP for the given temperature
    esat_ice = sat_vapor_press_ice( T_in_K )

#ifdef GFDL
    ! GFDL uses specific humidity
    ! Formula for Saturation Specific Humidity
     if( I_sat_sphum )  then   ! h1g, 2010-06-18 begin mod
           sat_mixrat_ice = ep * ( esat_ice / ( p_in_Pa - (1.0-ep) * esat_ice ) )
     else
           sat_mixrat_ice  = ep * ( esat_ice / ( p_in_Pa - esat_ice ) )
     endif                     ! h1g, 2010-06-18 end mod
#else
    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v

    sat_mixrat_ice = ep * ( esat_ice / ( p_in_Pa - esat_ice ) )
#endif

    return

  end function sat_mixrat_ice

!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice( T_in_K ) result ( esat_ice )
!
! Description:
!   Computes SVP for ice, using one of the various approximations.
!
! References:
!   None
!------------------------------------------------------------------------
 
    use model_flags, only: &
      saturation_formula, & ! Variable(s)
      saturation_bolton, &
#ifdef GFDL
      saturation_gfdl, &
#endif
      saturation_flatau

    implicit none

    ! Input Variable
    real, intent(in) :: &
      T_in_K      ! Temperature     [K]

    ! Output Variable
    real :: esat_ice    ! Saturation Vapor Pressure over Ice [Pa]

    ! Undefined approximation
    esat_ice = -99999.999

    select case ( saturation_formula )
    case ( saturation_bolton )
      ! Using the Bolton 1980 approximations for SVP over ice
      esat_ice = sat_vapor_press_ice_bolton( T_in_K )

    case ( saturation_flatau )
      ! Using the Flatau, et al. polynomial approximation for SVP over ice
      esat_ice = sat_vapor_press_ice_flatau( T_in_K )

    ! Add new cases after this
! ---> h1g, 2010-06-16
#ifdef GFDL
    case ( saturation_gfdl )
      ! Using GFDL polynomial approximation for SVP with respect to ice
      esat_ice = sat_vapor_press_ice_gfdl( T_in_K )
#endif
! <--- h1g, 2010-06-16
    end select

    return

  end function sat_vapor_press_ice

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
    use constants_clubb, only: T_freeze_K


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
    use constants_clubb, only: T_freeze_K

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


! ---> h1g, 2010-06-16
#ifdef GFDL
!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice_gfdl( T_in_K ) result ( esati )
! Description:
! copy from "GFDL polysvp.F90" 
!  Compute saturation vapor pressure with respect to liquid  by using 
! function from Goff and Gatch (1946)
! 
!  Polysvp returned in units of pa.
!  T_in_K is input in units of K.
!------------------------------------------------------------------------
 
    implicit none

    ! Input Variables
    real, intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real :: esati  ! Saturation vapor pressure over ice [Pa]

! Goff Gatch equation (good down to -100 C)

          esati = 10.**(-9.09718*(273.16/T_in_k-1.)-3.56654* &
          log10(273.16/T_in_k)+0.876793*(1.-T_in_k/273.16)+ &
          log10(6.1071))*100.

    return

  end function sat_vapor_press_ice_gfdl
#endif
! <--- h1g, 2010-06-16

!-------------------------------------------------------------------------
  FUNCTION sat_rcm( thlm, rtm, p_in_Pa, exner )

    ! Description:
    !
    ! This function uses an iterative method to find the value of rcm
    ! from an initial profile that has saturation at some point.
    !-------------------------------------------------------------------------

    use constants_clubb, only: & 
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

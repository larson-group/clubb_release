!$Id: saturation.F,v 1.4 2008-05-30 16:24:25 faschinj Exp $

        module saturation

!       Description: 
!         Contains functions that compute saturation with respect
!       to liquid or ice.

        implicit none

        private        ! Change default so all items private

        public   :: rsat, rsati, sat_rcm

        private  :: svp_liq, svp_ice

        contains

!-------------------------------------------------------------------------
        real function rsat( p, T )

!       Description:
!       Used to compute the saturation mixing ratio.

!       References:
!       Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

        use constants, only:
     .      ep ! Variable

        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  p, ! Pressure    [Pa]
     .  T  ! Temperature [K]

        ! Local Variables
        real :: esatv

        ! Saturation Vapor Pressure, esat, can be found to be approximated
        ! in many different ways.

        esatv = svp_liq( T )

        ! Formula for Saturation Mixing Ratio:
        !
        ! rs = (epsilon) * [ es / ( p - es ) ];
        ! where epsilon = R_d / R_v
       
        rsat = ep * ( esatv / ( p - esatv ) )

        return
        end function rsat

!------------------------------------------------------------------------
        pure function svp_liq( T ) result ( esat )

!       Description:
!       Computes SVP for water vapor.

!       References:
!       ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!         and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!         pp. 1507--1513
!------------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: exp
        
        ! Parameter Constants
        logical, parameter :: lFlatau = .false.

        ! Relative error norm expansion (-50 to 50 deg_C) from
        ! Table 3 of pp. 1511 of Flatau et al. 1992 (Water Vapor)
        real, dimension(7), parameter :: a =
     .  (/ 6.11176750,      0.443986062,     0.143053301E-01,
     .     0.265027242E-03, 0.302246994E-05, 0.203886313E-07,
     .     0.638780966E-10 /)

        ! Input Variables
        real, intent(in) :: T   ! Temperature   [K]

        ! Output Variables
        real :: esat

        ! Local Variables
        integer :: i

        if ( lFlatau ) then
          ! Polynomial approx. (Flatau, et al. 1992)
          esat = a(1)
          do i = 2, 7, 1
            esat = esat + a(i) * ( T-273.15 )**(i-1)
          end do
          esat = 100.0 * esat ! Convert units

        else
          ! (Bolton 1980) approx.
          esat = 611.2 * exp( (17.67*(T-273.15)) / (T-29.65) )
        end if

        return
        end function svp_liq

!------------------------------------------------------------------------
        real function rsati( p, T )

!       Description:
!       Used to compute the saturation mixing ratio. 

!       References:
!       Formula from Emanuel 1994, 4.4.15
!-------------------------------------------------------------------------

        use constants, only:
     .      ep ! Variable(s)

        implicit none

        ! Input Variables
        real, intent(in) :: p, T

        ! Local Variables
        real :: esat_ice

        ! Compute SVP for ice

        esat_ice = svp_ice( T )

        ! Formula for Saturation Mixing Ratio:
        !
        ! rs = (epsilon) * [ es / ( p - es ) ];
        ! where epsilon = R_d / R_v

        rsati = ep * ( esat_ice / ( p - esat_ice ) )

        return
        end function rsati

!------------------------------------------------------------------------
        real pure function svp_ice( T ) result ( esati )

!       Description:

!       References:
!       ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!         and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!         pp. 1507--1513
!------------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: exp, log
        
        ! Parameter Constants
        logical, parameter :: lFlatau = .false.

        ! Relative error norm expansion (-50 to 0 deg_C) from
        ! Table 3 of pp. 1511 of Flatau et al. 1992 (Ice)
        real, dimension(7), parameter :: a =
     .  (/ 6.10952665,      0.501948366,     0.18628899E-01,
     .     0.403488906E-03, 0.539797852E-05, 0.420713632E-07,
     .     0.147271071E-09 /)

        ! Input Variables
        real, intent(in) :: T   ! Temperature   [K]

        ! Local Variables
        integer :: i

        if ( lFlatau ) then
          ! Polynomial approx. (Flatau, et al. 1992)
          esati = a(1)
          do i = 2, 7, 1
            esati = esati + a(i) * ( T-273.15 )**(i-1)
          end do
          esati = 100.0 * esati ! Convert units

        else
          ! Exponential approx. (Bolton?)
          esati = 100.0 *
     .            exp( 23.33086 - (6111.72784/T) + (0.15215*log( T )) )
        end if

        return
        end function svp_ice

!------------------------------------------------------------------------
        ! This function was created by Brian Griffin.
        ! If you have information about your initial profile in thlm and rtm,
        ! and your profile is saturated at some point, you need to use an
        ! iterative method in order to accurately find the value of rcm.

        FUNCTION sat_rcm( thlm, rtm, p, exner )

        USE constants, only:
     .      Cp, ! Variable(s)
     .      Lv

        REAL, INTENT(IN):: thlm
        REAL, INTENT(IN):: rtm
        REAL, INTENT(IN):: p
        REAL, INTENT(IN):: exner

        REAL:: sat_rcm

        REAL:: theta
        REAL:: answer, too_low, too_high

        INTEGER:: iteration

        REAL, PARAMETER:: tolerance = 0.001


        theta = thlm
        iteration = 0

        DO

           iteration = iteration + 1

           answer = theta - (Lv/(Cp*exner))
     .                     *(MAX( rtm - rsat(p,theta*exner), 0.0 ))

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

        sat_rcm = MAX( rtm - rsat(p,theta*exner), 0.0 )

        END FUNCTION sat_rcm
        
        end module saturation

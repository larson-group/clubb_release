! $Id: rsat.F,v 1.6 2007-06-06 15:06:45 dschanen Exp $

!-------------------------------------------------------------------------
        real function rsat( p, T )

!       Description:
!       Used to compute the saturation mixing ratio.

!       References:
!       Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

        use constants

        implicit none

        ! External
        real, external :: svp_liq

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

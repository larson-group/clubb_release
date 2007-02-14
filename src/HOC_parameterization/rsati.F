!-------------------------------------------------------------------------
! $Id: rsati.F,v 1.4 2007-02-14 21:08:07 dschanen Exp $
        real function rsati( p, T )

!       Description:
!       Used to compute the saturation mixing ratio. 

!       References:
!       Formula from Emanuel 1994, 4.4.15
!-------------------------------------------------------------------------

        use constants

        implicit none

        ! External
        real, external :: esati

        ! Input Variables
        real, intent(in) :: p, T

        ! Local Variables
        real :: esat_ice

        ! Compute SVP for ice
        esat_ice = esati( T )

        ! Formula for Saturation Mixing Ratio:
        !
        ! rs = (epsilon) * [ es / ( p - es ) ];
        ! where epsilon = R_d / R_v
        rsati = ep * ( esat_ice / ( p - esat_ice ) )

        return
        end function rsati

!------------------------------------------------------------------------
        real pure function esati( T )

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
        end function esati
!------------------------------------------------------------------------

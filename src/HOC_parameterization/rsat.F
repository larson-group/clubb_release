! $Id: rsat.F,v 1.3 2006-11-18 17:26:27 griffinb Exp $

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
        intrinsic :: exp

        ! Input Variables
        real, intent(in) :: 
     .  p, ! Pressure    [Pa]
     .  T  ! Temperature [K]

        ! Local Variables
        real :: esat

        ! Saturation Vapor Pressure, esat, can be found to be approximated
        ! many differenet ways.

        ! rs = (epsilon) * [ es / ( p - es ) ];
        ! where epsilon = R_d / R_v

        esat = 611.2 * exp( (17.67*(T-273.15)) / (T-29.65) )
        rsat = ep * ( esat / ( p - esat ) )

        return
        end function rsat
!------------------------------------------------------------------------

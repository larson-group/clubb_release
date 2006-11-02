! $Id: rsat.F,v 1.2 2006-11-02 01:45:02 dschanen Exp $

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

        esat = 611.2 * exp( (17.67*(T-273.15)) / (T-29.65) )
        rsat = ep*esat/(p-(1-ep)*esat)

        return
        end function rsat
!------------------------------------------------------------------------

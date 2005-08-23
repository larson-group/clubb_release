!-------------------------------------------------------------------------
! $Id: rsat.F,v 1.3 2005-08-23 21:11:10 dschanen Exp $

! real function rstat( )

! Used to compute the saturation mixing ratio. Formula
! from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------
        real function rsat( p, T )

        use constants

        implicit none

!       Input
        real, intent(in) :: p, T

!       Internal
        real esat

        esat = 611.2 * exp( (17.67*(T-273.15)) / (T-29.65) )
        rsat = ep*esat/(p-(1-ep)*esat)

        return
        end function rsat
!------------------------------------------------------------------------

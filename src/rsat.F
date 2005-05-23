!-------------------------------------------------------------------------
! Function to compute the saturation mixing ratio. Formula
! from Emanuel 1994, 4.4.14

        real function rsat(p,T)

        use constants

        implicit none

        real p,T
        real esat

        esat = 611.2 * exp( (17.67*(T-273.15)) / (T-29.65) )
        rsat = ep*esat/(p-(1-ep)*esat)

        return
        end

!------------------------------------------------------------------------

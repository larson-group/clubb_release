!-------------------------------------------------------------------------
! $Id: rsati.F,v 1.3 2006-11-18 18:53:09 griffinb Exp $

! real function rstati( )

! Used to compute the saturation mixing ratio. Formula
! from Emanuel 1994, 4.4.15
! Michael Falk, 25 Aug 2006
!-------------------------------------------------------------------------
        real function rsati( p, T )

        use constants

        implicit none

!       Input
        real, intent(in) :: p, T

!       Internal
        real esati

        ! Formula for Saturation Mixing Ratio:
        !
        ! rs = (epsilon) * [ es / ( p - es ) ];
        ! where epsilon = R_d / R_v

        esati = 100 * exp (23.33086 - (6111.72784/T) + (0.15215*log(T)))
        rsati = ep * ( esati / ( p - esati ) )

        return
        end function rsati
!------------------------------------------------------------------------

!------------------------------------------------------------------------
! $Id: hydrostatic.F,v 1.5 2005-08-23 21:11:10 dschanen Exp $

! SUBROUTINE hydrostatic( )

! Subprogram to integrate hydrostatic equation
!------------------------------------------------------------------------
        subroutine hydrostatic(thvm, psfc, p, exner, rhot, rhom)

        use constants
        use grid_class

        implicit none

! Input

        real, intent(in)        :: thvm(1:gr%nnzp)
        real, intent(in)        :: psfc

! Output

        real, intent(out)       :: p(1:gr%nnzp)
        real, intent(out)       :: exner(1:gr%nnzp)
        real, intent(out)       :: rhot(1:gr%nnzp)   ! thermodynamics points
        real, intent(out)       :: rhom(1:gr%nnzp)   ! momentum points

! Internal

        integer k

! Integrate hydrostatic equation: we first compute Exner function
! on the momentum grid

        exner(1) = ( psfc/p0 ) ** kappa
        do k=2,gr%nnzp
           exner(k) = exner(k-1) - grav/( Cp * thvm(k) * gr%dzt(k) )
        end do

! Now interpolate Exner to the thermodynamic grid points

        exner = zm2zt( exner )

! Exner is defined on the thermodynamic grid point except for the first
! element which corresponds to surface value

        exner(1) = ( psfc/p0 ) ** kappa

! Compute pressure on thermodynamic points

        do k=1,gr%nnzp
           p(k) = p0 * exner(k)**(1./kappa)
        end do

! Compute density on thermodynamic grid

        do k=1,gr%nnzp
           rhot(k) = p(k) / ( R * thvm(k) * exner(k) )
        end do

! Interpolate density back to momentum grid

        rhom = zt2zm( rhot )
        rhom(1) = p(1) / ( R * thvm(1) * exner(1) )

        return
        end subroutine hydrostatic

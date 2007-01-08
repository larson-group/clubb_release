!------------------------------------------------------------------------
! $Id: hydrostatic.F,v 1.3 2007-01-08 20:37:45 dschanen Exp $

        subroutine hydrostatic( thvm, psfc, pinpa, exner, rhot, rhom )

!       Description:
!       Subprogram to integrate hydrostatic equation

!       References:
!
!------------------------------------------------------------------------

        use constants
        use grid_class

        implicit none

        ! Input Variables
        real, intent(in) :: psfc ! Pressure at the surface      [Pa]

        real, intent(in), dimension(gr%nnzp) :: 
     .  thvm  ! Virtual potential temperature   [K]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  pinpa, ! Pressure                       [Pa]
     .  exner, ! Exner function                 [-]
     .  rhot,  ! Density on thermo. points      [kg/m^3]
     .  rhom   ! Density on moment. points      [kg/m^3]

        !  Local Variables

        integer :: k

        ! Integrate hydrostatic equation: we first compute Exner function
        ! on the momentum grid

        exner(1) = ( psfc/p0 )**kappa
        do k=2,gr%nnzp
          exner(k) = exner(k-1) - grav/( Cp * thvm(k) * gr%dzt(k) )
        end do

        ! Now interpolate Exner to the thermodynamic grid points

        exner = zm2zt( exner )

        ! Exner is defined on the thermodynamic grid point except for the first
        ! element which corresponds to surface value

        ! Note: kappa = Rd / Cp

        exner(1) = ( psfc/p0 )**kappa

        ! Compute pressure on thermodynamic points

        do k=1,gr%nnzp
          pinpa(k) = p0 * exner(k)**( 1./kappa )
        end do

        ! Compute density on thermodynamic grid

        do k=1,gr%nnzp
          rhot(k) = pinpa(k) / ( Rd * thvm(k) * exner(k) )
        end do

        ! Interpolate density back to momentum grid

        rhom = zt2zm( rhot )
        rhom(1) = pinpa(1) / ( Rd * thvm(1) * exner(1) )

        return
        end subroutine hydrostatic

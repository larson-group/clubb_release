!------------------------------------------------------------------------
! $Id: hydrostatic.F,v 1.7 2008-07-21 20:31:54 faschinj Exp $

        module hydrostatic_mod

        implicit none

        private ! Default Scope

        public :: hydrostatic

        contains
        
        subroutine hydrostatic( thvm, psfc, p_in_Pa, exner, rhot, rhom )

!       Description:
!       Subprogram to integrate hydrostatic equation

!       References:
!
!------------------------------------------------------------------------

        use constants, only:
     .      kappa, ! Variable(s)
     .      p0,
     .      Cp,
     .      grav,
     .      Rd
        use grid_class, only:
     .      gr, ! Variable(s)
     .      zm2zt, ! Procedure(s)
     .      zt2zm


        implicit none

        ! Input Variables
        real, intent(in) :: psfc ! Pressure at the surface      [Pa]

        real, intent(in), dimension(gr%nnzp) :: 
     .  thvm  ! Virtual potential temperature   [K]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  p_in_Pa, ! Pressure                       [Pa]
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
          p_in_Pa(k) = p0 * exner(k)**( 1./kappa )
        end do

        ! Compute density on thermodynamic grid

        do k=1,gr%nnzp
          rhot(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
        end do

        ! Interpolate density back to momentum grid

        rhom = max( zt2zm( rhot ), 0.0 )   ! Positive definite quantity
        rhom(1) = p_in_Pa(1) / ( Rd * thvm(1) * exner(1) )

        return
        end subroutine hydrostatic

        end module hydrostatic_mod

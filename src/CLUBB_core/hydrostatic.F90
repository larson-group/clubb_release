!------------------------------------------------------------------------
! $Id$

module hydrostatic_mod

  implicit none

  private ! Default Scope

  public :: hydrostatic, inverse_hydrostatic

  contains

  subroutine hydrostatic( thvm, psfc, &
                          p_in_Pa, exner, rho, rho_zm )
    !       Description:
    !       Subprogram to integrate hydrostatic equation

    !       References:
    !
    !------------------------------------------------------------------------

    use constants, only: & 
        kappa,  & ! Variable(s)
        p0, & 
        Cp, & 
        grav, & 
        Rd, &
        zero_threshold

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zm2zt,  & ! Procedure(s)
        zt2zm


    implicit none

    ! Input Variables
    real, intent(in) :: psfc ! Pressure at the surface      [Pa]

    real, intent(in), dimension(gr%nnzp) ::  & 
      thvm  ! Virtual potential temperature   [K]

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      p_in_Pa,  & ! Pressure                       [Pa]
      exner,    & ! Exner function                 [-]
      rho,      & ! Density on thermo. points      [kg/m^3]
      rho_zm      ! Density on moment. points      [kg/m^3]

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
      rho(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
    end do

    ! Interpolate density back to momentum grid

    rho_zm = max( zt2zm( rho ), zero_threshold )   ! Positive definite quantity
    rho_zm(1) = p_in_Pa(1) / ( Rd * thvm(1) * exner(1) )

    return
  end subroutine hydrostatic

  subroutine inverse_hydrostatic( thvm, zm_init, exner, nVar, &
                                   z )
    !       Description:
    !       Subprogram to integrate the inverse of hydrostatic equation

    !       References:
    !
    !------------------------------------------------------------------------

    use constants, only: & 
        kappa,  & ! Variable(s)
        p0, & 
        Cp, & 
        grav, & 
        Rd, &
        zero_threshold

    implicit none

    ! Input Variables
    real, intent(in) :: zm_init ! Pressure at the surface      [Pa]

    integer, intent(in) :: nVar ! Number of points in the profile

    real, intent(in), dimension(nVar) ::  & 
      thvm, &  ! Virtual potential temperature   [K]
      exner    ! Exner function [-]

    ! Output Variables
    real, intent(out), dimension(nVar) ::  & 
      z        ! Height                    [m]

    !  Local Variables
    integer :: k

    real, dimension(nVar) :: zm_snd, exner_zm, d_exner_zm

    do k=1, nVar-1
      exner_zm(k) = 0.5 * ( exner( k ) + exner( k+1 ) )
    end do

    exner_zm(nVar) = exner(nVar) + 0.5 * ( exner(nVar) - exner(nVar -1) )

    zm_snd(1) = zm_init

    do k=2, nVar
      d_exner_zm(k) = exner_zm(k) - exner_zm(k-1)
      zm_snd(k) = zm_snd(k-1) - ( Cp / grav ) * thvm(k) * d_exner_zm(k)
    end do

    z(1) = 0

    do k = 2, nVar, 1
      z(k) = 0.5 * ( zm_snd(k) + zm_snd(k-1) )
    enddo

    return

  end subroutine inverse_hydrostatic

end module hydrostatic_mod

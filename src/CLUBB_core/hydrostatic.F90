!------------------------------------------------------------------------
!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module hydrostatic_mod

  implicit none

  private ! Default Scope

  public :: hydrostatic, inverse_hydrostatic

  private :: calc_exner_const_theta, &
             calc_exner_linear_theta

  contains

!===============================================================================
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

    ! Variables for test code for new formulation of integrating over the
    ! hydrostatic equation.  Brian; 12/28/2009.
    real, dimension(gr%nnzp) ::  &
      thvm_zm,    & ! Theta_v interpolated to momentum levels           [K]
      exner_zm,   & ! Exner function computed on momentum levels        [-]
      p_in_Pa_zm    ! Pressure computed on momentum levels              [Pa]

    real :: &
      dthv_dz       ! Constant d(theta_v)/dz between successive levels  [K/m]
    ! Brian.

    integer :: k

    ! Integrate hydrostatic equation: we first compute Exner function
    ! on the momentum grid

    exner(1) = ( psfc/p0 )**kappa
    do k = 2, gr%nnzp
       exner(k) &
       = calc_exner_const_theta( thvm(k), gr%zm(k), gr%zm(k-1), exner(k-1) )
    enddo

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

!    ! Test code for new formulation of integrating over the hydrostatic
!    ! equation.  Brian; 12/28/2009.
!    thvm_zm = zt2zm( thvm )
!
!    exner(1) = ( psfc/p0 )**kappa
!    exner_zm(1) = ( psfc/p0 )**kappa
!
!    ! Consider the value of exner at thermodynamic level (2) to be based on
!    ! a constant theta between thermodynamic level (2) and momentum level (1),
!    ! which is the surface or model lower boundary.  Since thlm(1) is set equal
!    ! to thlm(2), the values of theta are considered to be basically constant
!    ! near the ground.
!    exner(2) &
!    = calc_exner_const_theta( thvm(2), gr%zt(2), gr%zm(1), exner(1) )
!
!    ! Given the value of exner at thermodynamic level k-1, and considering
!    ! theta to vary linearly between its values at thermodynamic levels k
!    ! and k-1, the value of exner can be found at thermodynamic level k,
!    ! as well as at intermediate momentum level k-1.
!    do k = 3, gr%nnzp
!
!       dthv_dz = gr%dzm(k-1) * ( thvm(k) - thvm(k-1) )
!
!       if ( dthv_dz /= 0.0 ) then
!
!          exner(k) &
!          = calc_exner_linear_theta( thvm(k-1), dthv_dz, &
!                                     gr%zt(k-1), gr%zt(k), exner(k-1) )
!
!          exner_zm(k-1) &
!          = calc_exner_linear_theta( thvm(k-1), dthv_dz, &
!                                     gr%zt(k-1), gr%zm(k-1), exner(k-1) )
!
!       else ! dthv_dz = 0
!
!          exner(k) &
!          = calc_exner_const_theta &
!               ( thvm(k), gr%zt(k), gr%zt(k-1), exner(k-1) )
!
!          exner_zm(k-1) &
!          = calc_exner_const_theta &
!               ( thvm(k), gr%zm(k-1), gr%zt(k-1), exner(k-1) )
!
!       endif
!
!    enddo ! k = 3, gr%nnzp
!
!    ! Find the value of exner_zm at momentum level gr%nnzp by using a linear
!    ! extension of theta from the two thermodynamic level immediately below
!    ! momentum level gr%nnzp.
!    dthv_dz = ( thvm_zm(gr%nnzp) - thvm(gr%nnzp) ) &
!              / ( gr%zm(gr%nnzp) - gr%zt(gr%nnzp) )
!
!    if ( dthv_dz /= 0.0 ) then
!
!       exner_zm(gr%nnzp) &
!       = calc_exner_linear_theta &
!            ( thvm(gr%nnzp), dthv_dz, &
!              gr%zt(gr%nnzp), gr%zm(gr%nnzp), exner(gr%nnzp) )
!
!    else ! dthv_dz = 0
!
!       exner_zm(gr%nnzp) &
!       = calc_exner_const_theta &
!            ( thvm(gr%nnzp), gr%zm(gr%nnzp), gr%zt(gr%nnzp), exner(gr%nnzp) )
!
!    endif
!
!    do k = 1, gr%nnzp
!       p_in_Pa(k) = p0 * exner(k)**( 1./kappa )
!       p_in_Pa_zm(k) = p0 * exner_zm(k)**( 1./kappa )
!    enddo
!
!    do k = 1, gr%nnzp
!       rho(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
!       rho_zm(k) = p_in_Pa_zm(k) / ( Rd * thvm_zm(k) * exner_zm(k) )
!    enddo

    return
  end subroutine hydrostatic

!===============================================================================
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

!===============================================================================
  pure function calc_exner_const_theta( theta, z_2, z_1, exner_1 ) &
  result( exner_2 )

    ! Description:
    ! This function solves for exner at a level, given exner at another level,
    ! the altitudes of both levels, and a constant theta over the depth of the
    ! level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * theta).
    !
    ! This equation is integrated to solve for exner, such that:
    !
    ! INT(exner_1:exner_2) d(exner) = - ( grav / Cp ) INT(z_1:z_2) (1/theta) dz.
    !
    ! Since theta is considered to be a constant over the depth of the layer
    ! between z_1 and z_2, the equation can be written as:
    !
    ! INT(exner_1:exner_2) d(exner) = - grav / ( Cp * theta ) INT(z_1:z_2) dz.
    !
    ! Solving the integral:
    !
    ! exner_2 = exner_1 - [ grav / ( Cp * theta ) ] * ( z_2 - z_1 ).

    ! References:
    !-------------------------------------------------------------------

    use constants, only: &
        grav,  & ! Gravitational acceleration                  [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure [J/(kg*K)]

    implicit none

    ! Input Variables
    real, intent(in) :: &
      theta,   & ! Constant value of theta over the layer  [K]
      z_2,     & ! Altitude at the top of the layer        [m]
      z_1,     & ! Altitude at the bottom of the layer     [m]
      exner_1    ! Exner at the bottom of the layer        [-]

    ! Return Variable
    real :: exner_2  ! Exner at the top of the layer       [-]

    ! Calculate exner at top of the layer.
    exner_2 = exner_1 - ( grav / ( Cp * theta ) ) * ( z_2 - z_1 )

    return
  end function calc_exner_const_theta

!===============================================================================
  pure function calc_exner_linear_theta( theta_k, dtheta_dz, &
                                         z_k, z_2, exner_k  ) &
  result( exner_2 )

    ! Description:
    ! This function solves for exner at a level, given exner at another level,
    ! the altitudes of both levels, and a value of theta that is considered to
    ! vary linearly over the depth of the level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * theta).
    !
    ! This equation is integrated to solve for exner, such that:
    !
    ! INT(exner_1:exner_2) d(exner) = - ( grav / Cp ) INT(z_1:z_2) (1/theta) dz.
    !
    ! The value of theta is considered to vary linearly over the depth of the
    ! level (resulting in a constant d(theta)/dz over the depth of the level).
    ! The entire level must be encompassed between two levels with two known
    ! values of theta.  The value of theta at the upper level (z_up) is called
    ! theta_up, and the value of theta at the lower level (z_low) is called
    ! theta_low.  Again, the values of theta at all interior altitudes,
    ! z_low <= z <= z_up, behave linearly between theta_low and theta_up, such
    ! that:
    !
    ! theta(z) = [ ( theta_up - theta_low ) / ( z_up - z_low ) ] * ( z - z_low)
    !            + theta_low
    !          = [ d(theta)/dz ] * ( z - z_low ) + theta_low
    !          = C_a*z + C_b;
    !
    ! where:
    !
    ! C_a = ( theta_up - theta_low ) / ( z_up - z_low )
    !     = d(theta)/dz; and
    ! C_b = theta_low - [ ( theta_up - theta_low ) / ( z_up - z_low ) ] * z_low
    !     = theta_low - [ d(theta)/dz ] * z_low.
    !
    ! The integral becomes:
    !
    ! INT(exner_1:exner_2) d(exner) 
    ! = - ( grav / Cp ) INT(z_1:z_2) [ 1 / ( C_a*z + C_b ) ] dz.
    !
    ! Performing a u-substitution ( u = C_a*z + C_b ), the equation becomes:
    !
    ! INT(exner_1:exner_2) d(exner)
    ! = - ( grav / Cp ) * ( 1 / C_a ) INT(z=z_1:z=z_2) (1/u) du.
    !
    ! Solving the integral, and then re-substituting for u:
    !
    ! exner_2 = exner_1 
    !           - ( grav / Cp ) * ( 1 / C_a )
    !             * ln [ ( C_a*z_2 + C_b ) / ( C_a*z_1 + C_b ) ].
    !
    ! Re-substituting for C_a and C_b:
    !
    ! exner_2 
    ! = exner_1 
    !   - ( grav / Cp ) * ( 1 / {d(theta)/dz} )
    !     * ln [   ( {d(theta)/dz}*z_2 + theta_low - {d(theta)/dz}*z_low )
    !            / ( {d(theta)/dz}*z_1 + theta_low - {d(theta)/dz}*z_low ) ].
    !
    ! This equation is used to calculate exner_2 using exner_1, which is at the
    ! same level as z_1.  Furthermore, theta_low and z_low are taken from the
    ! same level as z_1 and exner_1.  Thus, z_1 = z_low.  Therefore:
    !
    ! exner_2 
    ! = exner_1 
    !   - ( grav / Cp ) * ( 1 / {d(theta)/dz} )
    !     * ln [ ( theta_low + {d(theta)/dz} * ( z_2 - z_low ) ) / theta_low ].
    !
    ! Considering either a thermodynamic or sounding level k as the low level
    ! in the integration, and that theta varies linearly between level k and
    ! level k+1:
    !
    ! exner_2
    ! = exner(k)
    !   - ( grav / Cp ) * ( 1 / {d(theta)/dz} )
    !     * ln [ ( theta(k) + {d(theta)/dz} * ( z_2 - z(k) ) ) / theta(k) ];
    !
    ! where:
    !
    ! d(theta)/dz = ( theta(k+1) - theta(k) ) / ( z(k+1) - z(k) );
    !
    ! and where z(k) < z_2 <= z(k+1); and {d(theta)/dz} /= 0.  If the value of
    ! {d(theta)/dz} is 0, then theta is considered to be a constant over the
    ! depth of the level.  The appropriate equation is found in pure function
    ! calc_exner_const_theta.

    ! References:
    !-------------------------------------------------------------------

    use constants, only: &
        grav,  & ! Gravitational acceleration                  [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure [J/(kg*K)]

    implicit none

    ! Input Variables
    real, intent(in) :: &
      theta_k,   & ! Value of theta at level k                     [K]
      dtheta_dz, & ! Constant d(theta)/dz between levels k and k+1 [K/m]
      z_k,       & ! Altitude at level k                           [m]
      z_2,       & ! Altitude at the top of the layer              [m]
      exner_k      ! Exner at level k                              [-]

    ! Return Variable
    real :: exner_2 ! Exner at the top of the layer                [-]

    ! Calculate exner at the top of the layer.
    exner_2 = exner_k &
              - ( grav / Cp ) * ( 1.0 / dtheta_dz ) &
                * log(  ( theta_k + dtheta_dz * ( z_2 - z_k ) )  /  theta_k  )

    return
  end function calc_exner_linear_theta

!===============================================================================

end module hydrostatic_mod

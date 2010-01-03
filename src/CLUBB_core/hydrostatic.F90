!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module hydrostatic_mod

  implicit none

  private ! Default Scope

  public :: hydrostatic, inverse_hydrostatic

  private :: calc_exner_const_th_var, &
             calc_exner_linear_th_var

  contains

!===============================================================================
  subroutine hydrostatic( th_var, psfc, &
                          p_in_Pa, exner, rho, rho_zm )

    ! Description:
    ! This subroutine integrates the hydrostatic equation of the form:
    !
    ! d(exner)/dz = - grav / ( Cp * th_var );
    !
    ! where th_var can be either virtual potential temperature (thvm), which
    ! is used to calculate total exner, pressure, and density, or potential
    ! temperature (thm), which is used to calculate dry, static, base-state
    ! exner, pressure, and density.
    !
    ! Sometimes, this subroutine is called with liquid water potential
    ! temperature (thlm) as the th_var.  This allows for a calculation of
    ! approximate pressure and exner, which allows initial cloud water mixing
    ! ratio (rcm) to be solved through an iterative method.  Then, thm is
    ! calculated from thlm and rcm.  Then, this subroutine is called again with
    ! thm as th_var in order to calculate dry, static, base-state variables.
    !
    ! After dry, static, base-state values have been determined, this subroutine
    ! is called with thvm as th_var in order to determined the total (moist)
    ! values of exner, pressure, and density.

    ! References:
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
    real, intent(in) :: &
      psfc    ! Pressure at the surface                     [Pa]

    real, intent(in), dimension(gr%nnzp) ::  & 
      th_var  ! Theta variable:  virtual potential temperature 
              !                  or potential temperature   [K]

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      p_in_Pa,  & ! Pressure (thermodynamic levels)         [Pa]
      exner,    & ! Exner function (thermodynamic levels)   [-]
      rho,      & ! Density (thermodynamic levels)          [kg/m^3]
      rho_zm      ! Density on momentum levels              [kg/m^3]

    !  Local Variables
    real, dimension(gr%nnzp) ::  &
      th_var_zm,  & ! Theta varariable interpolated to momentum levels  [K]
      exner_zm,   & ! Exner function on momentum levels                 [-]
      p_in_Pa_zm    ! Pressure on momentum levels                       [Pa]

    real :: &
      dth_var_dz    ! Constant d(th_var_v)/dz between successive levels [K/m]

    integer :: k

    ! Interpolate th_var from thermodynamic to momentum levels.  Linear
    ! interpolation is used, except for the uppermost momentum level, where a
    ! linear extension is used.  Since th_var is considered to either be
    ! constant or vary linearly over the depth of a grid level, this
    ! interpolation is consistent with the rest of this code.
    th_var_zm = zt2zm( th_var )

    ! Exner is defined on thermodynamic grid levels except for the value at
    ! index 1.  Since thermodynamic level 1 is below the surface, it is
    ! disregarded, and the value of exner(1) corresponds to surface value, which
    ! is actually at momentum level 1.
    exner(1) = ( psfc/p0 )**kappa
    exner_zm(1) = ( psfc/p0 )**kappa

    ! Consider the value of exner at thermodynamic level (2) to be based on
    ! a constant th_var between thermodynamic level (2) and momentum level (1),
    ! which is the surface or model lower boundary.  Since thlm(1) is set equal
    ! to thlm(2), the values of th_var are considered to be basically constant
    ! near the ground.
    exner(2) &
    = calc_exner_const_th_var( th_var(2), gr%zt(2), gr%zm(1), exner(1) )

    ! Given the value of exner at thermodynamic level k-1, and considering
    ! th_var to vary linearly between its values at thermodynamic levels k
    ! and k-1, the value of exner can be found at thermodynamic level k,
    ! as well as at intermediate momentum level k-1.
    do k = 3, gr%nnzp

       dth_var_dz = gr%dzm(k-1) * ( th_var(k) - th_var(k-1) )

       if ( dth_var_dz /= 0.0 ) then

          exner(k) &
          = calc_exner_linear_th_var( th_var(k-1), dth_var_dz, &
                                      gr%zt(k-1), gr%zt(k), exner(k-1) )

          exner_zm(k-1) &
          = calc_exner_linear_th_var( th_var(k-1), dth_var_dz, &
                                      gr%zt(k-1), gr%zm(k-1), exner(k-1) )

       else ! dth_var_dz = 0

          exner(k) &
          = calc_exner_const_th_var &
               ( th_var(k), gr%zt(k), gr%zt(k-1), exner(k-1) )

          exner_zm(k-1) &
          = calc_exner_const_th_var &
               ( th_var(k), gr%zm(k-1), gr%zt(k-1), exner(k-1) )

       endif

    enddo ! k = 3, gr%nnzp

    ! Find the value of exner_zm at momentum level gr%nnzp by using a linear
    ! extension of th_var from the two thermodynamic level immediately below
    ! momentum level gr%nnzp.
    dth_var_dz = ( th_var_zm(gr%nnzp) - th_var(gr%nnzp) ) &
                 / ( gr%zm(gr%nnzp) - gr%zt(gr%nnzp) )

    if ( dth_var_dz /= 0.0 ) then

       exner_zm(gr%nnzp) &
       = calc_exner_linear_th_var &
            ( th_var(gr%nnzp), dth_var_dz, &
              gr%zt(gr%nnzp), gr%zm(gr%nnzp), exner(gr%nnzp) )

    else ! dth_var_dz = 0

       exner_zm(gr%nnzp) &
       = calc_exner_const_th_var &
            ( th_var(gr%nnzp), gr%zm(gr%nnzp), gr%zt(gr%nnzp), exner(gr%nnzp) )

    endif

    ! Calculate pressure based on the values of exner.

    do k = 1, gr%nnzp
       p_in_Pa(k) = p0 * exner(k)**( 1./kappa )
       p_in_Pa_zm(k) = p0 * exner_zm(k)**( 1./kappa )
    enddo

    ! Calculate density based on pressure, exner, and th_var.

    do k = 1, gr%nnzp
       rho(k) = p_in_Pa(k) / ( Rd * th_var(k) * exner(k) )
       rho_zm(k) = p_in_Pa_zm(k) / ( Rd * th_var_zm(k) * exner_zm(k) )
    enddo


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
  pure function calc_exner_const_th_var( th_var, z_2, z_1, exner_1 ) &
  result( exner_2 )

    ! Description:
    ! This function solves for exner at a level, given exner at another level,
    ! the altitudes of both levels, and a constant th_var over the depth of the
    ! level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * th_var).
    !
    ! This equation is integrated to solve for exner, such that:
    !
    ! INT(exner_1:exner_2) d(exner) 
    ! = - ( grav / Cp ) INT(z_1:z_2) (1/th_var) dz.
    !
    ! Since th_var is considered to be a constant over the depth of the layer
    ! between z_1 and z_2, the equation can be written as:
    !
    ! INT(exner_1:exner_2) d(exner) = - grav / ( Cp * th_var ) INT(z_1:z_2) dz.
    !
    ! Solving the integral:
    !
    ! exner_2 = exner_1 - [ grav / ( Cp * th_var ) ] * ( z_2 - z_1 ).

    ! References:
    !-------------------------------------------------------------------

    use constants, only: &
        grav,  & ! Gravitational acceleration                  [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure [J/(kg*K)]

    implicit none

    ! Input Variables
    real, intent(in) :: &
      th_var,  & ! Constant value of th_var over the layer  [K]
      z_2,     & ! Altitude at the top of the layer         [m]
      z_1,     & ! Altitude at the bottom of the layer      [m]
      exner_1    ! Exner at the bottom of the layer         [-]

    ! Return Variable
    real :: exner_2  ! Exner at the top of the layer        [-]

    ! Calculate exner at top of the layer.
    exner_2 = exner_1 - ( grav / ( Cp * th_var ) ) * ( z_2 - z_1 )

    return
  end function calc_exner_const_th_var

!===============================================================================
  pure function calc_exner_linear_th_var( th_var_km1, dth_var_dz, &
                                          z_km1, z_2, exner_km1  ) &
  result( exner_2 )

    ! Description:
    ! This function solves for exner at a level, given exner at another level,
    ! the altitudes of both levels, and a value of th_var that is considered to
    ! vary linearly over the depth of the level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * th_var).
    !
    ! This equation is integrated to solve for exner, such that:
    !
    ! INT(exner_1:exner_2) d(exner)
    ! = - ( grav / Cp ) INT(z_1:z_2) (1/th_var) dz.
    !
    ! The value of th_var is considered to vary linearly over the depth of the
    ! level (resulting in a constant d(th_var)/dz over the depth of the level).
    ! The entire level between z_1 and z_2 must be encompassed between two
    ! levels with two known values of th_var.  The value of th_var at the upper
    ! level (z_up) is called th_var_up, and the value of th_var at the lower
    ! level (z_low) is called th_var_low.  Again, the values of th_var at all
    ! interior altitudes, z_low <= z_1 < z <= z_2 <= z_up, behave linearly
    ! between th_var_low and th_var_up, such that:
    !
    ! th_var(z)
    ! = [ ( th_var_up - th_var_low ) / ( z_up - z_low ) ] * ( z - z_low)
    !   + th_var_low
    ! = [ d(th_var)/dz ] * ( z - z_low ) + th_var_low
    ! = C_a*z + C_b;
    !
    ! where:
    !
    ! C_a 
    ! = ( th_var_up - th_var_low ) / ( z_up - z_low )
    ! = d(th_var)/dz;
    !
    ! and:
    !
    ! C_b
    ! = th_var_low - [ ( th_var_up - th_var_low ) / ( z_up - z_low ) ] * z_low
    ! = th_var_low - [ d(th_var)/dz ] * z_low.
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
    !   - ( grav / Cp ) * ( 1 / {d(th_var)/dz} )
    !     * ln [   ( {d(th_var)/dz}*z_2 + th_var_low - {d(th_var)/dz}*z_low )
    !            / ( {d(th_var)/dz}*z_1 + th_var_low - {d(th_var)/dz}*z_low ) ].
    !
    ! This equation is used to calculate exner_2 using exner_1, which is at the
    ! same level as z_1.  Furthermore, th_var_low and z_low are taken from the
    ! same level as z_1 and exner_1.  Thus, z_1 = z_low.  Therefore:
    !
    ! exner_2 
    ! = exner_low
    !   - ( grav / Cp ) * ( 1 / {d(th_var)/dz} )
    !     * ln [ ( th_var_low + {d(th_var)/dz}*(z_2-z_low) ) / th_var_low ].
    !
    ! Considering either a thermodynamic or sounding level k-1 as the low level
    ! in the integration, and that th_var varies linearly between level k-1 and
    ! level k:
    !
    ! exner_2
    ! = exner(k-1)
    !   - ( grav / Cp ) * ( 1 / {d(th_var)/dz} )
    !     * ln [ ( th_var(k-1) + {d(th_var)/dz}*(z_2-z(k-1)) ) / th_var(k-1) ];
    !
    ! where:
    !
    ! d(th_var)/dz = ( th_var(k) - th_var(k-1) ) / ( z(k) - z(k-1) );
    !
    ! and where z(k-1) < z_2 <= z(k); and {d(th_var)/dz} /= 0.  If the value of
    ! {d(th_var)/dz} is 0, then th_var is considered to be a constant over the
    ! depth of the level.  The appropriate equation is found in pure function
    ! calc_exner_const_th_var.

    ! References:
    !-------------------------------------------------------------------

    use constants, only: &
        grav,  & ! Gravitational acceleration                   [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure  [J/(kg*K)]

    implicit none

    ! Input Variables
    real, intent(in) :: &
      th_var_km1, & ! Value of th_var at level k-1                    [K]
      dth_var_dz, & ! Constant d(th_var)/dz between levels k-1 and k  [K/m]
      z_km1,      & ! Altitude at level k-1                           [m]
      z_2,        & ! Altitude at the top of the layer                [m]
      exner_km1     ! Exner at level k-1                              [-]

    ! Return Variable
    real :: exner_2 ! Exner at the top of the layer                [-]

    ! Calculate exner at the top of the layer.
    exner_2  &
    = exner_km1 &
      - ( grav / Cp ) * ( 1.0 / dth_var_dz )  &
        * log(  ( th_var_km1 + dth_var_dz * ( z_2 - z_km1 ) )  /  th_var_km1  )

    return
  end function calc_exner_linear_th_var

!===============================================================================

end module hydrostatic_mod

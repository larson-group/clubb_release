!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module hydrostatic_module

  implicit none

  private ! Default Scope

  public :: hydrostatic, &
            inverse_hydrostatic

  private :: calc_z_linear_thvm

  contains

!===============================================================================
  subroutine hydrostatic( thvm, p_sfc, &
                          p_in_Pa, p_in_Pa_zm, &
                          exner, exner_zm, &
                          rho, rho_zm )

    ! Description:
    ! This subroutine integrates the hydrostatic equation.
    !
    ! The hydrostatic equation is of the form:
    !
    ! dp/dz = - rho * grav.
    !
    ! This equation can be re-written in terms of d(exner)/dz, such that:
    !
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) } / { R_d * rho } ] * d(exner)/dz
    ! = - grav / C_p;
    !
    ! which can also be expressed as:
    !
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) } / { R_d * rho_d * ( 1 + r_v + r_c ) } ]
    !    * d(exner)/dz
    ! = - grav / C_p.
    !
    ! Furthermore, the moist equation of state can be written as:
    !
    ! theta =
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) }
    !   / { R_d * rho_d * ( 1 + (R_v/R_d)*r_v ) } ].
    !
    ! The relationship between theta and theta_v (including water vapor and
    ! cloud water) is:
    !
    ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
    !
    ! which, when substituted into the above equation, changes the equation of
    ! state to:
    !
    ! theta_v =
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) }
    !   / { R_d * rho_d * ( 1 + r_v + r_c ) } ].
    !
    ! This equation is substituted into the d(exner)/dz form of the hydrostatic
    ! equation, resulting in:
    !
    ! theta_v * d(exner)/dz = - grav / C_p;
    !
    ! which can be re-written as:
    !
    ! d(exner)/dz = - grav / ( C_p * theta_v ).
    !
    ! This subroutine integrates the above equation to solve for exner, such
    ! that:
    !
    ! INT(exner_1:exner_2) d(exner) =
    ! - ( grav / C_p ) * INT(z_1:z_2) ( 1 / theta_v ) dz.
    !
    !
    ! The resulting value of exner is used to calculate pressure.  Then, the
    ! values of pressure, exner, and theta_v can be used to calculate density.

    ! References:
    !
    !------------------------------------------------------------------------

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zt2zm

    use constants_clubb, only: & 
        Rd    ! Gas Constant for Dry Air  [J/(kg K)]

    use calc_pressure, only: &
        init_pressure    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      p_sfc    ! Pressure at the surface                     [Pa]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      thvm    ! Virtual potential temperature               [K]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  & 
      p_in_Pa,    & ! Pressure (thermodynamic levels)         [Pa]
      p_in_Pa_zm, & ! Pressure on momentum levels             [Pa]
      exner,      & ! Exner function (thermodynamic levels)   [-]
      exner_zm,   & ! Exner function on momentum levels       [-]
      rho,        & ! Density (thermodynamic levels)          [kg/m^3]
      rho_zm        ! Density on momentum levels              [kg/m^3]

    !  Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      thvm_zm       ! Theta_v interpolated to momentum levels  [K]

    integer :: k

    ! Calculate pressure and exner on both thermodynamic and momentum levels.
    call init_pressure( thvm, p_sfc, &
                        p_in_Pa, exner, p_in_Pa_zm, exner_zm )

    ! Interpolate thvm from thermodynamic to momentum levels.  Linear
    ! interpolation is used, except for the uppermost momentum level, where a
    ! linear extension is used.  Since thvm is considered to either be constant
    ! or vary linearly over the depth of a grid level, this interpolation is
    ! consistent with the rest of this code.
    thvm_zm = zt2zm( thvm )

    ! Calculate density based on pressure, exner, and thvm.
    do k = 1, gr%nz
      rho(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
      rho_zm(k) = p_in_Pa_zm(k) / ( Rd * thvm_zm(k) * exner_zm(k) )
    enddo


    return

  end subroutine hydrostatic

!===============================================================================
  subroutine inverse_hydrostatic( p_sfc, zm_init, nlevels, thvm, exner, &
                                  z )

    ! Description:
    ! Subprogram to integrate the inverse of hydrostatic equation

    ! References:
    !
    !------------------------------------------------------------------------

    use constants_clubb, only: &
        p0,     & ! Constant(s)
        kappa,  &
        fstderr

    use interpolation, only: &
        binary_search ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      p_sfc,    & ! Pressure at the surface      [Pa]
      zm_init    ! Altitude at the surface      [m]

    integer, intent(in) ::  &
      nlevels  ! Number of levels in the sounding [-]

    real( kind = core_rknd ), intent(in), dimension(nlevels) ::  & 
      thvm,  & ! Virtual potential temperature   [K]
      exner    ! Exner function                  [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(nlevels) ::  & 
      z        ! Height                    [m]

    !  Local Variables
    integer :: k

    real( kind = core_rknd ), dimension(nlevels) ::  &
      ref_z_snd  ! Altitude minus altitude of the lowest sounding level  [m]

    real( kind = core_rknd ), dimension(nlevels) ::  &
      exner_reverse_array  ! Array of exner snd. values in reverse order [-]

    real( kind = core_rknd ) ::  &
      exner_sfc,    & ! Value of exner at the surface                    [-]
      ref_z_sfc,    & ! Alt. diff between surface and lowest snd. level  [m]
      z_snd_bottom, & ! Altitude of the bottom of the input sounding     [m]
      dthvm_dexner    ! Constant rate of change of thvm with respect to
    ! exner between sounding levels k-1 and k          [K]

    integer ::  &
      rev_low_idx, &
      low_idx, &
      high_idx


    ! Variable ref_z_sfc is initialized to 0.0 to avoid a compiler warning.
    ref_z_sfc = 0.0_core_rknd

    ! The variable ref_z_snd is the altitude of each sounding level compared to
    ! the altitude of the lowest sounding level.  Thus, the value of ref_z_snd
    ! at sounding level 1 is 0.  The lowest sounding level may or may not be
    ! right at the surface, and therefore an adjustment may be required to find
    ! the actual altitude above ground.
    ref_z_snd(1) = 0.0_core_rknd

    do k = 2, nlevels

      ! The value of thvm is given at two successive sounding levels.  For
      ! purposes of achieving a quality estimate of altitude at each pressure
      ! sounding level, the value of thvm is considered to vary linearly
      ! with respect to exner between two successive sounding levels.  Thus,
      ! there is a constant d(thvm)/d(exner) between the two successive
      ! sounding levels.  If thvm is constant, then d(thvm)/d(exner) is 0.
      dthvm_dexner = ( thvm(k) - thvm(k-1) ) / ( exner(k) - exner(k-1) )

      ! Calculate the value of the reference height at sounding level k, based
      ! the value of thvm at sounding level k-1, the constant value of
      ! d(thvm)/d(exner), the value of exner at sounding levels k-1 and k, and
      ! the reference altitude at sounding level k-1.
      ref_z_snd(k) &
      = calc_z_linear_thvm( thvm(k-1), dthvm_dexner, &
                            exner(k-1), exner(k), ref_z_snd(k-1) )

    enddo

    ! Find the actual (above ground) altitude of the sounding levels from the
    ! reference altitudes.

    ! The pressure at the surface (or model lower boundary), p_sfc, is found at
    ! the altitude of the surface (or model lower boundary), zm_init.

    ! Find the value of exner at the surface from the pressure at the surface.
    exner_sfc = ( p_sfc / p0 )**kappa

    ! Find the value of exner_sfc compared to the values of exner in the exner
    ! sounding profile.

    if ( exner_sfc < exner(nlevels) ) then

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is less than all the
      ! values of exner in the sounding (and thus the surface is located above
      ! all the levels of the sounding), then there is insufficient information
      ! to run the model.  Stop the run.

      write(fstderr,*) "The entire sounding is below the model surface."
      stop

    elseif ( exner_sfc > exner(1) ) then

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is greater than all the
      ! values of exner in the sounding (and thus the surface is located below
      ! all the levels of the sounding), use a linear extension of thvm to find
      ! thvm at the surface.  Thus, d(thvm)/d(exner) is the same as its value
      ! between sounding levels 1 and 2.  If the surface is so far below the
      ! sounding that gr%zt(2) is below the first sounding level, the code in
      ! subroutine read_sounding (found in sounding.F90) will stop the run.

      ! Calculate the appropriate d(thvm)/d(exner).
      dthvm_dexner = ( thvm(2) - thvm(1) ) / ( exner(2) - exner(1) )

      ! Calculate the difference between the altitude of the surface (or model
      ! lower boundary) and the altitude of the lowest level of the sounding.
      ref_z_sfc  &
      = calc_z_linear_thvm( thvm(1), dthvm_dexner, &
                            exner(1), exner_sfc, ref_z_snd(1) )

    else  ! exner(nlevels) < exner_sfc < exner(1)

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is between two values of
      ! exner (at some levels k-1 and k) in the sounding, and the value of
      ! d(thvm)/d(exner) is the same as between those two levels in the above
      ! calculation.

      ! The value of exner_sfc is between two levels of the exner sounding.
      ! Find the index of the lower level.

      ! In order to use the binary search, the array must be sorted from least
      ! value to greatest value.  Since exner decreases with altitude (and
      ! vertical level), the array that is sent to function binary_search must
      ! be the exact reverse of exner.
      ! Thus, exner(1) becomes exner_reverse_array(nlevels), exner(nlevels)
      ! becomes exner_reverse_array(1), etc.
      do k = 1, nlevels, 1
        exner_reverse_array(k) = exner(nlevels-k+1)
      enddo
      ! The output from the binary search yields the first value in the
      ! exner_reverse_array that is greater than or equal to exner_sfc.  Thus,
      ! in regards to the regular exner array, this is the reverse index of
      ! the lower sounding level for exner_sfc.  For example, if exner_sfc
      ! is found between exner(1) and exner(2), the binary search for exner_sfc
      ! in regards to exner_reverse_index will return a value of nlevels.
      ! Once the actual lower level index is calculated, the result will be 1.
      rev_low_idx = binary_search( nlevels, exner_reverse_array, exner_sfc )

      ! Find the lower level index for the regular exner profile from the
      ! lower level index for the reverse exner profile.
      low_idx = nlevels - rev_low_idx + 1

      ! Find the index of the upper level.
      high_idx = low_idx + 1

      ! Calculate the appropriate d(thvm)/d(exner).
      dthvm_dexner = ( thvm(high_idx) - thvm(low_idx) )  &
                       /  ( exner(high_idx) - exner(low_idx) )

      ! Calculate the difference between the altitude of the surface (or model
      ! lower boundary) and the altitude of the lowest level of the sounding.
      ref_z_sfc  &
      = calc_z_linear_thvm( thvm(low_idx), dthvm_dexner, &
                            exner(low_idx), exner_sfc, ref_z_snd(low_idx) )

    endif  ! exner_sfc

    ! Find the altitude of the bottom of the sounding.
    z_snd_bottom = zm_init - ref_z_sfc

    ! Calculate the sounding altitude profile based
    ! on z_snd_bottom and ref_z_snd.
    do k = 1, nlevels, 1
      z(k) = z_snd_bottom + ref_z_snd(k)
    enddo


    return
  end subroutine inverse_hydrostatic

!===============================================================================
  pure function calc_z_linear_thvm( thvm_km1, dthvm_dexner, &
                                    exner_km1, exner_2, z_km1 ) &
  result( z_2 )

    ! Description:
    ! This function solves for z (altitude) at a level, given altitude at
    ! another level, the values of exner at both levels, and a value of thvm
    ! that is considered to vary linearly over the depth of the level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * thvm).
    !
    ! This equation is integrated to solve for z, such that:
    !
    ! INT(exner_1:exner_2) thvm d(exner) = - ( grav / Cp ) INT(z_1:z_2) dz.
    !
    ! The value of thvm is considered to vary linearly (with respect to exner)
    ! over the depth of the level (resulting in a constant d(thvm)/d(exner) over
    ! the depth of the level).  The entire level between exner_1 and exner_2
    ! must be encompassed between two levels with two known values of thvm.  The
    ! value of thvm at the upper level (exner_up) is called thvm_up, and the
    ! value of thvm at the lower level (exner_low) is called thvm_low.  Again,
    ! the values of thvm at all interior exner levels,
    ! exner_low >= exner_1 > exner >= exner_2 >= exner_up, behave linearly
    ! between thvm_low and thvm_up, such that:
    !
    ! thvm(exner)
    ! = [ ( thvm_up - thvm_low ) / ( exner_up - exner_low ) ]
    !     * ( exner - exner_low )
    !   + thvm_low
    ! = [ d(thvm)/d(exner) ] * ( exner - exner_low ) + thvm_low
    ! = C_a*z + C_b;
    !
    ! where:
    !
    ! C_a
    ! = ( thvm_up - thvm_low ) / ( exner_up - exner_low )
    ! = d(thvm)/d(exner);
    !
    ! and:
    !
    ! C_b
    ! = thvm_low
    !   - [ ( thvm_up - thvm_low ) / ( exner_up - exner_low ) ] * exner_low
    ! = thvm_low - [ d(thvm)/d(exner) ] * exner_low.
    !
    ! The integral becomes:
    !
    ! INT(exner_1:exner_2) ( C_a*exner + C_b ) d(exner)
    ! = - ( grav / Cp ) INT(z_1:z_2) dz.
    !
    ! Solving the integral:
    !
    ! z_2
    ! = z_1
    !   - ( Cp / grav )
    !     * [    (1/2) * {d(thvm)/d(exner)} * ( {exner_2}^2 - {exner_1}^2 )
    !          + ( thvm_low - {d(thvm)/d(exner)} * exner_low )
    !            * ( exner_2 - exner_1 )  ].
    !
    ! This equation is used to calculate z_2 using z_1, which is at the same
    ! level as exner_1.  Furthermore, thvm_low and exner_low are taken from the
    ! same level as exner_1 and z_1.  Thus, exner_1 = exner_low.  Therefore:
    !
    ! z_2
    ! = z_low
    !   - ( Cp / grav )
    !     * [    (1/2) * {d(thvm)/d(exner)} * ( {exner_2}^2 - {exner_low}^2 )
    !          + ( thvm_low - {d(thvm)/d(exner)} * exner_low )
    !            * ( exner_2 - exner_low )  ].
    !
    ! Considering a sounding level k-1 as the low level in the integration, and
    ! that thvm varies linearly (with respect to exner) between level k-1 and
    ! level k:
    !
    ! z_2
    ! = z(k-1)
    !   - ( Cp / grav )
    !     * [    (1/2) * {d(thvm)/d(exner)} * ( {exner_2}^2 - {exner(k-1)}^2 )
    !          + ( thvm(k-1) - {d(thvm)/d(exner)} * exner(k-1) )
    !            * ( exner_2 - exner(k-1) )  ];
    !
    ! where:
    !
    ! d(thvm)/d(exner)
    ! = ( thvm(k) - thvm(k-1) ) / ( exner(k) - exner(k-1) );
    !
    ! and where exner(k-1) > exner_2 >= exner(k).  If the value of
    ! d(thvm)/d(exner) is 0, then thvm is considered to be a constant over the
    ! depth of the level, and the equation will reduce to:
    !
    ! z_2 = z(k-1) - ( Cp / grav ) * thvm(k-1) * ( exner_2 - exner(k-1) ).
    !
    !
    ! IMPORTANT NOTE:
    !
    ! CLUBB is an altitude-based model.  All linear interpolations (and
    ! extensions) are based on considering a variable to change linearly with
    ! respect to altitude, rather than with respect to exner.  An exception is
    ! made here to calculate the altitude of a sounding level based on a
    ! sounding given in terms of a pressure coordinate rather than a height
    ! coordinate.  After the altitude of the sounding level has been calculated,
    ! the values of the sounding variables are interpolated onto the model grid
    ! linearly with respect to altitude.  Therefore, considering a variable to
    ! change linearly with respect to exner is not consistent with the rest of
    ! the model code, but provides for a better estimation of the altitude of
    ! the sounding levels (than simply considering thvm to be constant over the
    ! depth of the sounding level).

    ! References:
    !-------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        grav,  & ! Gravitational acceleration                   [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure  [J/(kg*K)]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      thvm_km1,     & ! Value of thvm at sounding level k-1                 [K]
      dthvm_dexner, & ! Constant d(thvm)/d(exner) between levels k-1 and k  [K]
      exner_km1,    & ! Value of exner at sounding level k-1                [-]
      exner_2,      & ! Value of exner at the top of the layer              [-]
      z_km1           ! Altitude at sounding level k-1                      [m]

    ! Return Variable
    real( kind = core_rknd ) :: z_2       ! Altitude at the top of the layer                    [m]

    ! Calculate z_2 at the top of the layer.
    z_2  &
    = z_km1  &
      - ( Cp / grav )  &
        * (   0.5_core_rknd * dthvm_dexner * ( exner_2**2 - exner_km1**2 )  &
            + ( thvm_km1 - dthvm_dexner * exner_km1 )  &
              * ( exner_2 - exner_km1 )  &
          )

    return
  end function calc_z_linear_thvm

!===============================================================================

end module hydrostatic_module

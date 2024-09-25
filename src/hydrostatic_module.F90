!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module hydrostatic_module

  implicit none

  private ! Default Scope

  public :: hydrostatic, &
            inverse_hydrostatic

  private :: calc_ref_z_linear_thvm,     &
             calc_ref_z_sfc_linear_thvm

  contains

!===============================================================================
  subroutine hydrostatic( gr, thvm, p_sfc, &
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

    use grid_class, only: grid ! Type

    use grid_class, only: & 
        zt2zm

    use constants_clubb, only: & 
        Rd    ! Gas Constant for Dry Air  [J/(kg K)]

    use calc_pressure, only: &
        init_pressure    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      p_sfc    ! Pressure at the surface                     [Pa]

    real( kind = core_rknd ), intent(in), dimension(gr%nzt) ::  & 
      thvm    ! Virtual potential temperature               [K]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nzt) ::  & 
      p_in_Pa,    & ! Pressure (thermodynamic levels)         [Pa]
      exner,      & ! Exner function (thermodynamic levels)   [-]
      rho           ! Density (thermodynamic levels)          [kg/m^3]

    real( kind = core_rknd ), intent(out), dimension(gr%nzm) ::  & 
      p_in_Pa_zm, & ! Pressure on momentum levels             [Pa]
      exner_zm,   & ! Exner function on momentum levels       [-]
      rho_zm        ! Density on momentum levels              [kg/m^3]

    !  Local Variables
    real( kind = core_rknd ), dimension(gr%nzm) ::  &
      thvm_zm       ! Theta_v interpolated to momentum levels  [K]

    integer :: k

    ! Calculate pressure and exner on both thermodynamic and momentum levels.
    call init_pressure( gr, thvm, p_sfc, &
                        p_in_Pa, exner, p_in_Pa_zm, exner_zm )

    ! Interpolate thvm from thermodynamic to momentum levels.  Linear
    ! interpolation is used, except for the uppermost momentum level, where a
    ! linear extension is used.  Since thvm is considered to either be constant
    ! or vary linearly over the depth of a grid level, this interpolation is
    ! consistent with the rest of this code.
    thvm_zm = zt2zm( gr, thvm )

    ! Calculate density based on pressure, exner, and thvm.
    do k = 1, gr%nzt
      rho(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
    enddo
    do k = 1, gr%nzm
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
      z_snd_bottom    ! Altitude of the bottom of the input sounding     [m]

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
    ref_z_snd = calc_ref_z_linear_thvm( nlevels, thvm, exner )

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
      error stop

    elseif ( exner_sfc > exner(1) ) then

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is greater than all the
      ! values of exner in the sounding (and thus the surface is located below
      ! all the levels of the sounding), use a linear extension of thvm to find
      ! thvm at the surface.  Thus, d(thvm)/d(exner) is the same as its value
      ! between sounding levels 1 and 2.  If the surface is so far below the
      ! sounding that gr%zt(1,2) is below the first sounding level, the code in
      ! subroutine read_sounding (found in sounding.F90) will stop the run.

      ! Calculate the difference between the altitude of the surface (or model
      ! lower boundary) and the altitude of the lowest level of the sounding.
      ref_z_sfc  &
      = calc_ref_z_sfc_linear_thvm( thvm(1), thvm(2), z(1), z(2), &
                                 exner(1), exner_sfc )

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

      ! Calculate the difference between the altitude of the surface (or model
      ! lower boundary) and the altitude of the lowest level of the sounding.
      ref_z_sfc  &
      = calc_ref_z_sfc_linear_thvm( thvm(low_idx), thvm(high_idx), &
                                    ref_z_snd(low_idx), ref_z_snd(high_idx), &
                                    exner(low_idx), exner_sfc )

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
  pure function calc_ref_z_linear_thvm( nlevels, thvm, exner ) &
  result( z )

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
    ! INT(exner_1:exner_2) d(exner) = - ( grav / Cp ) INT(z_1:z_2) (1/thvm) dz.
    !
    ! The value of mean theta_v (thvm) is calculated at each level, and linear
    ! interpolation is used in the integral equation for all altitudes
    ! in-between successive levels, such that:
    !
    ! thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
    !
    ! The integrals are solved, and the equation can be solved for the thickness
    ! z2 - z1, such that:
    !
    ! z2 - z1
    !   | - ( Cp / grav )
    !   |   * ( exner2 - exner1 ) * ( thvm2 - thvm1 ) / ln( thvm2 / thvm1 );
    ! = | where thvm2 /= thvm1;
    !   |
    !   | - ( Cp / grav ) * ( exner2 - exner1 ) * thvm;
    !   | where thvm2 = thvm1 (= thvm).
    ! 
    ! The value of height (z) can be calculated using the above equation at all
    ! levels once the value of height is known at one level.

    ! References:
    !-------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        grav,  & ! Gravitational acceleration                   [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure  [J/(kg*K)]

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nlevels

    real( kind = core_rknd ), dimension(nlevels), intent(in) :: &
      thvm,  & ! Virtual potential temperature   [K]
      exner    ! Exner function                  [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(nlevels) :: &
      z        ! Height                    [m]

    integer :: k

    z(1) = 0.0_core_rknd

    ! Calculate the value of the reference height at sounding level k, based
    ! the value of thvm at sounding levels k-1 and k, the value of exner at
    ! sounding levels k-1 and k, and the reference altitude at sounding level
    ! k-1.
    do k = 2, nlevels

      if ( abs( thvm(k) - thvm(k-1) ) > epsilon( thvm(k) ) * thvm(k) ) then

        z(k) &
        = z(k-1) &
          - ( Cp / grav ) * ( exner(k) - exner(k-1) ) * ( thvm(k) - thvm(k-1) ) &
            / log( thvm(k) / thvm(k-1) )

      else ! thvm(k) = thvm(k-1)

        z(k) = z(k-1) - ( Cp / grav ) * ( exner(k) - exner(k-1) ) * thvm(k)

      endif
      
    end do


    return

  end function calc_ref_z_linear_thvm

!===============================================================================
  pure function calc_ref_z_sfc_linear_thvm( thvm_km1, thvm_k, z_km1, z_k, &
                                            exner_km1, exner_sfc ) &
  result( ref_z_sfc )

    ! Description:
    ! This function solves for the difference between the altitude of the model
    ! lower boundary and the altitude of a sounding level.  Since the altitude
    ! of the model lower boundary is known, this routine allows the actual
    ! altitude of a sounding level to be calculated.  Once a sounding level
    ! altitude is known, the altitudes of the remaining sounding levels can also
    ! be calculated.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * thvm).
    !
    ! This equation is integrated to solve for ref_z_sfc, such that:
    !
    ! INT(exner1:exner_sfc) d(exner)
    ! = - ( grav / Cp ) INT(z1:ref_z_sfc) (1/thvm) dz;
    !
    ! where exner1 and z1 are exner and altitude, respectively, at a sounding
    ! level, exner_sfc is surface exner (calculated from surface pressure), and
    ! ref_z_sfc is the height of the surface relative to sounding level z1.  
    !
    ! The value of mean theta_v (thvm) is calculated at each level, and linear
    ! interpolation is used in the integral equation for all altitudes
    ! in-between successive levels, such that:
    !
    ! thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
    !
    ! In the above equation, z2 - z1 is the thickness between two sounding
    ! levels (which has already been calculated), and thvm1 and thvm2 are the
    ! values of thvm at those two sounding levels.  The thickness z - z1 is
    ! the difference between altitude z and the sounding level z1.
    !
    ! The integrals are solved, and the equation can be solved for ref_z_sfc,
    ! such that:
    !
    ! ref_z_sfc
    !   | ( ( z2 - z1 ) / ( thvm2 - thvm1 ) )
    !   | * ( thvm1 * exp{ - ( Cp / grav ) * ( exner_sfc - exner1 )
    !   |                    * ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) }
    ! = |     - thvm1 + ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * z1 );
    !   | where thvm2 /= thvm1;
    !   |
    !   | z1 - ( Cp / grav ) * ( exner2 - exner1 ) * thvm;
    !   | where thvm2 = thvm1 (= thvm).
    ! 
    ! Prior to calculating ref_z_sfc, the thickness between all sounding levels
    ! were already calculated.  The calculation of ref_z_sfc allows for the
    ! difference between the surface and a sounding level to be calculated
    ! (ref_z_sfc - z1).  Since the altitude of the surface is already provided
    ! (zm_init), the ref_z_sfc - z1 calculation allows for the calculation of
    ! the altitude of a sounding level.  Once the altitude of a sounding level
    ! is calculated, the altitude of all other sounding levels are known, since
    ! the thicknesses between all sounding levels were already calculated.

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
      thvm_km1,  & ! Value of thvm at sounding level k-1          ` [K]
      thvm_k,    & ! Value of thvm at sounding level k              [K]
      z_km1,     & ! Reference altitude at sounding level k-1       [m]
      z_k,       & ! Reference altitude at sounding level k         [m]
      exner_km1, & ! Value of exner at sounding level k-1           [-]
      exner_sfc    ! Value of exner at the surface                  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      ref_z_sfc   ! Altitude at of surface compared to sounding level k-1    [m]

    ! Calculate z at sounding level k.
    if ( abs( thvm_k - thvm_km1 ) > epsilon( thvm_k ) * thvm_k ) then

       ref_z_sfc &
       = ( ( z_k - z_km1 ) / ( thvm_k - thvm_km1 ) ) &
         * ( thvm_km1 * exp( - ( Cp / grav ) * ( exner_sfc - exner_km1 ) &
                               * ( ( thvm_k - thvm_km1 ) / ( z_k - z_km1 ) ) ) &
             - thvm_km1 + ( ( thvm_k - thvm_km1 ) / ( z_k - z_km1 ) ) * z_km1 )

    else ! thvm_k = thvm_km1

       ref_z_sfc = z_km1 - ( Cp / grav ) * ( exner_sfc - exner_km1 ) * thvm_k

    endif


    return

  end function calc_ref_z_sfc_linear_thvm

!===============================================================================

end module hydrostatic_module

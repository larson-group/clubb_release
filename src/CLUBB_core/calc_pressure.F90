!===============================================================================
module calc_pressure

  implicit none

  public :: update_pressure, & ! Procedure(s)
            init_pressure,   &
            calculate_thvm

  private ! default scope

  contains

  !=============================================================================
  subroutine update_pressure( nz, ngrdcol, gr, &
                              thlm, rtm, rcm, rho_ds_zt, thv_ds_zt, &
                              p_in_Pa, exner, &
                              p_in_Pa_zm, exner_zm )

    ! Description:
    ! Updates pressure according to the hydrostatic approximation.  Combining
    ! the moist equation of state and the hydrostatic approximation, the change
    ! of pressure with respect to height can be calculated based on theta_v,
    ! such that:
    !
    ! dp/dz = - p * grav / ( Rd * theta_v * exner );
    !
    ! where exner = ( p / p0 )^(Rd/Cp);
    !
    ! and where p0 is a reference pressure of 100000 Pa.
    !
    ! The integral equation is set up to integrate over p on the left-hand side
    ! and integrate over z on the right-hand side.  The equation is:
    !
    ! INT(p1:p2) p^(Rd/Cp-1) dp
    ! = - p0^(Rd/Cp) * ( grav / Rd ) * INT(z1:z2) (1/thvm) dz.
    !
    ! The value of mean theta_v (thvm) is calculated at each thermodynamic grid
    ! level, and linear interpolation is used in the integral equation for all
    ! altitudes in-between successive thermodynamic levels, such that:
    !
    ! thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
    !
    ! The integrals are solved, and the results for pressure can be rewritten
    ! in terms of exner, such that:
    !
    ! exner2 - exner1
    !   | - ( grav / Cp )
    !   |   * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    ! = | where thvm2 /= thvm1;
    !   |
    !   | - ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    ! 
    ! The value of pressure (exner) can be calculated using the above equation
    ! at all vertical levels once the value of pressure is known at one level.
    !
    ! The model fields are usually least variant over time at the top of the
    ! model, so the value of pressure is calculated at the uppermost
    ! thermodynamic level, and then the values of pressure (exner) can be
    ! calculated at all other vertical levels (both thermodynamic levels and
    ! momentum levels).  The following equation is used to calculate pressure
    ! at the uppermost thermodynamic level:
    !
    ! p = ( ( rho_dry * Rd / p0^(Rd/Cp) )
    !       * thm * ( 1 + (Rv/Rd) * rvm ) )^(1/(1-Rd/Cp));
    !
    ! where rho_dry is the density of dry air, where rvm = rtm - rcm, and where:
    !
    ! thm = thlm + ( Lv / ( Cp * exner|_old ) ) * rcm;
    !
    ! and exner|_old is exner from the previous timestep.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid, & ! Type
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one,   & ! 1
        Rd,    & ! Gas constant for dry air                    [J/(kg K)]
        Lv,    & ! Latent heat of vaporization                 [J/kg]
        Cp,    & ! Specific heat of dry air                    [J/(kg K)]
        kappa, & ! Rd/Cp                                       [-]
        ep2,   & ! Rv/Rd (Rv is gas constant for water vapor)  [-]
        p0,    & ! Reference pressure of 100000 Pa             [Pa]
        grav     ! Acceleration of gravity (9.81 m/s^2)        [m/s^2]

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! ---------------------------- Input Variables ----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      thlm,      & ! Mean liquid water potential temperature         [K]
      rtm,       & ! Mean total water mixing ratio                   [kg/kg]
      rcm,       & ! Mean cloud water mixing ratio                   [kg/kg]
      rho_ds_zt, & ! Dry, state base-state density (thermo. levels)  [kg/m^3]
      thv_ds_zt    ! Reference theta_v on thermodynamic levels       [K]

    ! ---------------------------- Input/Output Variables ---------------------------- 
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      p_in_Pa,    & ! Pressure (thermodynamic levels)        [Pa]
      exner         ! Exner function (thermodynamic levels)  [-]

    ! ---------------------------- Output Variables ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      p_in_Pa_zm, & ! Pressure on momentum levels            [Pa]
      exner_zm      ! Exner function on momentum levels      [-]

    ! ---------------------------- Local Variables ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      thvm,    & ! Mean theta_v (thermodynamic levels)                [K]
      thvm_zm    ! Mean theta_v interpolated to momentum grid levels  [K]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      thm_nz, & ! Theta at the uppermost thermodynamic grid level    [K]
      rvm_nz    ! Water vapor mixing ratio; uppermost thermo. level  [kg/kg]

    real( kind = core_rknd ), parameter :: &
      g_ov_Cp = grav / Cp, &     ! g / Cp  [K/m]
      invrs_kappa = one / kappa  ! 1 / kappa

    ! Flag to calculate pressure and exner on momentum levels.  Otherwise,
    ! linear interpolation of exner will be used.
    logical, parameter :: &
      l_calc_p_exner_m_levs = .true.

    integer :: i, k  ! Vertical level index


    ! ---------------------------- Begin Code ----------------------------

    ! Calculate thvm on thermodynamic grid levels.
    call calculate_thvm( nz, ngrdcol, &
                         thlm, rtm, rcm, exner, thv_ds_zt, &
                         thvm )

    ! Interpolate thvm to momentum grid levels.
    thvm_zm = zt2zm( nz, ngrdcol, gr, thvm )

    do i = 1, ngrdcol
      ! Calculate mean theta (thm) at the uppermost thermodynamic vertical grid
      ! level.
      thm_nz(i) = thlm(i,nz) + ( Lv / ( Cp * exner(i,nz) ) ) * rcm(i,nz)

      ! Calculate mean water vapor mixing ratio (rvm) at the uppermost
      ! thermodynamic vertical grid level.
      rvm_nz(i) = rtm(i,nz) - rcm(i,nz)

      ! Update pressure at the uppermost thermodynamic grid level.
      p_in_Pa(i,nz) &
      = ( ( rho_ds_zt(i,nz) * Rd / p0**kappa ) &
          * thm_nz(i) * ( one + ep2 * rvm_nz(i) ) )**(one/(one-kappa))

      ! Calculate exner at the uppermost thermodynamic grid level.
      exner(i,nz) = ( p_in_Pa(i,nz) / p0 )**kappa

    end do

    ! Calculate exner at the uppermost momentum grid level, which is located
    ! above the uppermost thermodynamic grid level.
    ! exner2
    ! = exner1
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   - | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    if ( l_calc_p_exner_m_levs ) then

      do i = 1, ngrdcol

        if ( abs( thvm(i,nz) - thvm_zm(i,nz) ) &
            > epsilon( thvm ) * thvm(i,nz) ) then

          exner_zm(i,nz) &
          = exner(i,nz) &
            - g_ov_Cp * ( gr%zm(i,nz) - gr%zt(i,nz) ) &
              / ( thvm_zm(i,nz) - thvm(i,nz) ) &
              * log( thvm_zm(i,nz) / thvm(i,nz) )

        else ! thvm(k+1) = thvm_zm(k)

          exner_zm(i,nz) = exner(i,nz) &
                            - g_ov_Cp * ( gr%zm(i,nz) - gr%zt(i,nz) ) / thvm_zm(i,nz)
        endif
      end do
    else ! .not. l_calc_p_exner_m_levs

      ! Interpolate exner to momentum levels
      do i = 1, ngrdcol
        exner_zm(i,nz) = zt2zm( gr, exner(i,:), nz )
      end do

    endif ! l_calc_p_exner_m_levs

    ! Calculate pressure on the uppermost momentum level.
    do i = 1, ngrdcol
      p_in_Pa_zm(i,nz) = p0 * exner_zm(i,nz)**invrs_kappa
    end do

    ! Calculate exner at all other thermodynamic and momentum grid levels,
    ! which are all located below the uppermost thermodynamic grid level.
    ! exner1
    ! = exner2
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   + | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    do k = nz-1, 2, -1
      do i = 1, ngrdcol

        ! Calculate exner on thermodynamic levels.
        if ( abs( thvm(i,k+1) - thvm(i,k) ) > epsilon( thvm ) * thvm(i,k+1) ) then

          exner(i,k) &
          = exner(i,k+1) &
            + g_ov_Cp * ( gr%zt(i,k+1) - gr%zt(i,k) ) / ( thvm(i,k+1) - thvm(i,k) ) &
              * log( thvm(i,k+1) / thvm(i,k) )

        else ! thvm(k+1) = thvm(k)

          exner(i,k) = exner(i,k+1) + g_ov_Cp * ( gr%zt(i,k+1) - gr%zt(i,k) ) / thvm(i,k)

        endif

        if ( l_calc_p_exner_m_levs ) then

          ! Calculate exner on momentum levels.
          if ( abs( thvm(i,k+1) - thvm_zm(i,k) ) &
               > epsilon( thvm ) * thvm(i,k+1) ) then

             exner_zm(i,k) &
             = exner(i,k+1) &
               + g_ov_Cp * ( gr%zt(i,k+1) - gr%zm(i,k) ) &
                 / ( thvm(i,k+1) - thvm_zm(i,k) ) &
                 * log( thvm(i,k+1) / thvm_zm(i,k) )

          else ! thvm(k+1) = thvm_zm(k)

             exner_zm(i,k) &
             = exner(i,k+1) + g_ov_Cp * ( gr%zt(i,k+1) - gr%zm(i,k) ) / thvm_zm(i,k)

          endif

        else ! .not. l_calc_p_exner_m_levs

          ! Interpolate exner to momentum levels
          exner_zm(i,k) = zt2zm( gr, exner(i,:), k )

        endif ! l_calc_p_exner_m_levs

      enddo ! k = nz-1, 2, -1
    end do

    do k = 2, nz-1
      do i = 1, ngrdcol
        ! Calculate pressure on thermodynamic levels.
        p_in_Pa(i,k) = p0 * exner(i,k)**invrs_kappa

        ! Calculate pressure on momentum levels.
        p_in_Pa_zm(i,k) = p0 * exner_zm(i,k)**invrs_kappa
      end do
    end do

    do i = 1, ngrdcol
      ! Calculate exner the model lower boundary or surface.
      ! exner1
      ! = exner2
      !     | ( grav / Cp )
      !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
      !   + | where thvm2 /= thvm1;
      !     |
      !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
      if ( abs( thvm(i,2) - thvm_zm(i,1) ) > epsilon( thvm ) * thvm(i,2) ) then

         exner_zm(i,1) &
         = exner(i,2) &
           + g_ov_Cp * ( gr%zt(i,2) - gr%zm(i,1) ) / ( thvm(i,2) - thvm_zm(i,1) ) &
             * log( thvm(i,2) / thvm_zm(i,1) )

       else ! thvm(k+1) = thvm_zm(k)

         exner_zm(i,1) &
         = exner(i,2) + g_ov_Cp * ( gr%zt(i,2) - gr%zm(i,1) ) / thvm_zm(i,1)

      endif

      ! Calculate pressure at the model lower boundary or surface.
      p_in_Pa_zm(i,1) = p0 * exner_zm(i,1)**invrs_kappa

      ! For the lowest thermodynamic level, which is below the model lower
      ! boundary, set pressure and exner to the pressure and exner found at the
      ! model lower boundary.
      p_in_Pa(i,1) = p_in_Pa_zm(i,1)
      exner(i,1) = exner_zm(i,1)
    end do

    return

  end subroutine update_pressure

  !=============================================================================
  subroutine init_pressure( gr, thvm, p_sfc, &
                            p_in_Pa, exner, p_in_Pa_zm, exner_zm )

    ! Description:
    ! Calculates the initial pressure according to the hydrostatic
    ! approximation.  Combining the moist equation of state and the hydrostatic
    ! approximation, the change of pressure with respect to height can be
    ! calculated based on theta_v, such that:
    !
    ! dp/dz = - p * grav / ( Rd * theta_v * exner );
    !
    ! where exner = ( p / p0 )^(Rd/Cp);
    !
    ! and where p0 is a reference pressure of 100000 Pa.
    !
    ! The integral equation is set up to integrate over p on the left-hand side
    ! and integrate over z on the right-hand side.  The equation is:
    !
    ! INT(p1:p2) p^(Rd/Cp-1) dp
    ! = - p0^(Rd/Cp) * ( grav / Rd ) * INT(z1:z2) (1/thvm) dz.
    !
    ! The value of mean theta_v (thvm) is calculated at each thermodynamic grid
    ! level, and linear interpolation is used in the integral equation for all
    ! altitude in-between successive thermodynamic levels, such that:
    !
    ! thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
    !
    ! The integrals are solved, and the results for pressure can be rewritten
    ! in terms of exner, such that:
    !
    ! exner2 - exner1
    !   | - ( grav / Cp )
    !   |   * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    ! = | where thvm2 /= thvm1;
    !   |
    !   | - ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    ! 
    ! The value of pressure (exner) can be calculated using the above equation
    ! at all vertical levels once the value of pressure is known at one level.
    ! Since the surface pressure is known at the initial time, that allows
    ! pressure to be calculated for the rest of the vertical profile.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid, & ! Type
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one,   & ! 1
        Cp,    & ! Specific heat of dry air                    [J/(kg K)]
        kappa, & ! Rd/Cp                                       [-]
        p0,    & ! Reference pressure of 100000 Pa             [Pa]
        grav     ! Acceleration of gravity (9.81 m/s^2)        [m/s^2]

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      thvm    ! Mean theta_v (thermodynamic levels)                [K]

    real( kind = core_rknd ), intent(in) :: &
      p_sfc    ! Surface pressure                                   [Pa]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      p_in_Pa,    & ! Pressure (thermodynamic levels)        [Pa]
      exner,      & ! Exner function (thermodynamic levels)  [-]
      p_in_Pa_zm, & ! Pressure on momentum levels            [Pa]
      exner_zm      ! Exner function on momentum levels      [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      thvm_zm    ! Mean theta_v interpolated to momentum grid levels  [K]

    real( kind = core_rknd ), parameter :: &
      g_ov_Cp = grav / Cp  ! g / Cp  [K/m]

    integer :: k  ! Vertical level index


    ! The pressure (and exner) at the lowest momentum level is the surface
    ! pressure (and exner based on the surface pressure).
    p_in_Pa_zm(1) = p_sfc
    exner_zm(1) = ( p_sfc / p0 )**kappa

    ! Set the pressure (and exner) at the lowest thermodynamic level, which is
    ! below the model lower boundary, to their values at the model lower
    ! boundary or surface.
    p_in_Pa(1) = p_in_Pa_zm(1)
    exner(1) = exner_zm(1)

    ! Interpolate theta_v to momentum levels.
    thvm_zm = zt2zm( gr, thvm )

    ! Calculate exner at all other thermodynamic and momentum grid levels.
    ! exner2
    ! = exner1
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   - | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).

    ! Calculate exner at thermodynamic level 2 (first thermodynamic grid level
    ! that is above the lower boundary).
    if ( abs( thvm(2) - thvm_zm(1) ) > epsilon( thvm ) * thvm(2) ) then

       exner(2) &
       = exner_zm(1) &
         - g_ov_Cp * ( gr%zt(1,2) - gr%zm(1,1) ) / ( thvm(2) - thvm_zm(1) ) &
           * log( thvm(2) / thvm_zm(1) )

    else ! thvm(2) = thvm_zm(1)

       exner(2) = exner_zm(1) - g_ov_Cp * ( gr%zt(1,2) - gr%zm(1,1) ) / thvm(2)

    endif

    ! Calculate pressure on thermodynamic levels.
    p_in_Pa(2) = p0 * exner(2)**(one/kappa)

    ! Loop over all other thermodynamic vertical grid levels.
    do k = 3, gr%nz, 1

       ! Calculate exner on thermodynamic levels.
       if ( abs( thvm(k) - thvm(k-1) ) > epsilon( thvm ) * thvm(k) ) then

          exner(k) &
          = exner(k-1) &
            - g_ov_Cp * ( gr%zt(1,k) - gr%zt(1,k-1) ) / ( thvm(k) - thvm(k-1) ) &
              * log( thvm(k) / thvm(k-1) )

       else ! thvm(k+1) = thvm(k)

          exner(k) = exner(k-1) - g_ov_Cp * ( gr%zt(1,k) - gr%zt(1,k-1) ) / thvm(k)

       endif

       ! Calculate pressure on thermodynamic levels.
       p_in_Pa(k) = p0 * exner(k)**(one/kappa)

    enddo ! k = 2, gr%nz, 1

    ! Loop over all momentum grid levels.
    do k = 2, gr%nz, 1

       ! Calculate exner on momentum levels.
       if ( abs( thvm(k) - thvm_zm(k) ) > epsilon( thvm ) * thvm(k) ) then

          exner_zm(k) &
          = exner(k) &
            - g_ov_Cp * ( gr%zm(1,k) - gr%zt(1,k) ) / ( thvm_zm(k) - thvm(k) ) &
              * log( thvm_zm(k) / thvm(k) )

       else ! thvm(k) = thvm_zm(k)

          exner_zm(k) &
          = exner(k) - g_ov_Cp * ( gr%zm(1,k) - gr%zt(1,k) ) / thvm_zm(k)

       endif

       ! Calculate pressure on momentum levels.
       p_in_Pa_zm(k) = p0 * exner_zm(k)**(one/kappa)

    enddo ! k = 2, gr%nz, 1


    return

  end subroutine init_pressure

  !=============================================================================
  subroutine calculate_thvm( nz, ngrdcol, &
                             thlm, rtm, rcm, exner, thv_ds_zt, &
                             thvm) 
  ! Description:
  ! Calculates mean theta_v based on a linearized approximation to the theta_v
  ! equation.  This equation also includes liquid water loading.

  ! References:
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
        Lv,  & ! Latent Heat of Vaporizaion    [J/kg]
        Cp,  & ! Specific Heat of Dry Air      [J/(kg K)]
        ep1, & ! Rv/Rd - 1                     [-]
        ep2    ! Rv/Rd                         [-]

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! ---------------------------- Input Variables ----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      thlm,      & ! Mean theta_l (thermodynamic levels)          [K]
      rtm,       & ! Mean total water (thermodynamic levels)      [kg/kg]
      rcm,       & ! Mean cloud water (thermodynamic levels)      [kg/kg]
      exner,     & ! Exner function (thermodynamic levels)        [-]
      thv_ds_zt    ! Reference theta_v on thermodynamic levels    [K]

    ! ---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      thvm    ! Mean theta_v (thermodynamic levels)    [K]

    ! ---------------------------- Local Variables ----------------------------
    integer :: i, k

    ! ---------------------------- Begin Code ----------------------------

    ! Calculate mean theta_v
    do k = 1, nz
      do i = 1, ngrdcol
        thvm(i,k) = thlm(i,k) + ep1 * thv_ds_zt(i,k) * rtm(i,k) &
               + ( Lv / ( Cp * exner(i,k) ) - ep2 * thv_ds_zt(i,k) ) * rcm(i,k)
      end do
    end do

    return

  end subroutine calculate_thvm
  !=============================================================================

end module calc_pressure

!===============================================================================
module calc_pressure

  implicit none

  public :: init_pressure,   & ! Procedure(s)
            calculate_thvm

  private ! default scope

  contains

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
        grav,  & ! Acceleration of gravity (9.81 m/s^2)        [m/s^2]
        zero_threshold

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      thvm    ! Mean theta_v (thermodynamic levels)                [K]

    real( kind = core_rknd ), intent(in) :: &
      p_sfc    ! Surface pressure                                   [Pa]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: &
      p_in_Pa,    & ! Pressure (thermodynamic levels)        [Pa]
      exner         ! Exner function (thermodynamic levels)  [-]

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      p_in_Pa_zm, & ! Pressure on momentum levels            [Pa]
      exner_zm      ! Exner function on momentum levels      [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nzm) :: &
      thvm_zm    ! Mean theta_v interpolated to momentum grid levels  [K]

    real( kind = core_rknd ), parameter :: &
      g_ov_Cp = grav / Cp  ! g / Cp  [K/m]

    integer :: k  ! Vertical level index


    ! The pressure (and exner) at the lowest momentum level is the surface
    ! pressure (and exner based on the surface pressure).
    p_in_Pa_zm(1) = p_sfc
    exner_zm(1) = ( p_sfc / p0 )**kappa

    ! Interpolate theta_v to momentum levels.
    thvm_zm = zt2zm( gr, thvm, zero_threshold )

    ! Calculate exner at all other thermodynamic and momentum grid levels.
    ! exner2
    ! = exner1
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   - | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).

    ! Calculate exner at thermodynamic level 1 (first thermodynamic grid level
    ! that is above the lower boundary).
    if ( abs( thvm(1) - thvm_zm(1) ) > epsilon( thvm ) * thvm(1) ) then

       exner(1) &
       = exner_zm(1) &
         - g_ov_Cp * ( gr%zt(1,1) - gr%zm(1,1) ) / ( thvm(1) - thvm_zm(1) ) &
           * log( thvm(1) / thvm_zm(1) )

    else ! thvm(1) = thvm_zm(1)

       exner(1) = exner_zm(1) - g_ov_Cp * ( gr%zt(1,1) - gr%zm(1,1) ) / thvm(1)

    endif

    ! Calculate pressure on thermodynamic levels.
    p_in_Pa(1) = p0 * exner(1)**(one/kappa)

    ! Loop over all other thermodynamic vertical grid levels.
    do k = 2, gr%nzt, 1

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

    enddo ! k = 2, gr%nzt, 1

    ! Loop over all momentum grid levels.
    do k = 2, gr%nzm, 1

       ! Calculate exner on momentum levels.
       if ( abs( thvm(k-1) - thvm_zm(k) ) > epsilon( thvm ) * thvm(k-1) ) then

          exner_zm(k) &
          = exner(k-1) &
            - g_ov_Cp * ( gr%zm(1,k) - gr%zt(1,k-1) ) / ( thvm_zm(k) - thvm(k-1) ) &
              * log( thvm_zm(k) / thvm(k-1) )

       else ! thvm(k) = thvm_zm(k)

          exner_zm(k) &
          = exner(k-1) - g_ov_Cp * ( gr%zm(1,k) - gr%zt(1,k-1) ) / thvm_zm(k)

       endif

       ! Calculate pressure on momentum levels.
       p_in_Pa_zm(k) = p0 * exner_zm(k)**(one/kappa)

    enddo ! k = 2, gr%nzm, 1


    return

  end subroutine init_pressure

  !=============================================================================
  subroutine calculate_thvm( nzt, ngrdcol, &
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
      nzt, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      thlm,      & ! Mean theta_l (thermodynamic levels)          [K]
      rtm,       & ! Mean total water (thermodynamic levels)      [kg/kg]
      rcm,       & ! Mean cloud water (thermodynamic levels)      [kg/kg]
      exner,     & ! Exner function (thermodynamic levels)        [-]
      thv_ds_zt    ! Reference theta_v on thermodynamic levels    [K]

    ! ---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      thvm    ! Mean theta_v (thermodynamic levels)    [K]

    ! ---------------------------- Local Variables ----------------------------
    integer :: i, k

    ! ---------------------------- Begin Code ----------------------------

    !$acc data copyin( thlm, rtm, rcm, exner, thv_ds_zt ) & 
    !$acc     copyout( thvm )

    ! Calculate mean theta_v
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        thvm(i,k) = thlm(i,k) + ep1 * thv_ds_zt(i,k) * rtm(i,k) &
               + ( Lv / ( Cp * exner(i,k) ) - ep2 * thv_ds_zt(i,k) ) * rcm(i,k)
      end do
    end do
    !$acc end parallel loop
    
    !$acc end data

    return

  end subroutine calculate_thvm
  !=============================================================================

end module calc_pressure 

! $Id$
!===============================================================================
module microphys_utilities

  ! Description:
  ! Utilities used in relation to the general microphysics code.  This file
  ! contains the microphysics term limiters, such as the turbulent sedimentation
  ! flux limiter and the evaporation limiter, which limits excessive
  ! evaporation.

  implicit none

  private ! Default scope

  public :: turb_sed_flux_limiter

  contains

  !=============================================================================
  subroutine turb_sed_flux_limiter( solve_type, dt, &
                                    rho_ds_zm, invrs_rho_ds_zt, &
                                    hmm, Vhmphmp, Vhmphmp_zt )

    ! Description:
    ! Adjusts the turbulent sedimentation flux of a hydrometeor, as well as the
    ! value of the hydrometeor, when the turbulent sedimentation flux has
    ! caused the value of the hydrometeor to become negative.  Wherever a
    ! turbulent sedimentation flux (momentum levels) needs to be adjusted, the
    ! value of the hydrometeor needs to be adjusted at each of the thermodynamic
    ! levels that are adjacent to the momentum level where the turbulent
    ! sedimentation flux adjustment took place.
    !
    ! In this section, the hydrometeor (either a mixing ratio or a
    ! concentration) will be denoted hm, and the sedimentation velocity (fall
    ! velocity) of that hydrometeor will be denoted V_hm.
    !
    ! Many times, near the top of the precipitation region (usually near cloud
    ! top), the grid level mean value of a hydrometeor, < hm >, decreases very
    ! sharply with height.  Small values of < hm > are often found near the top
    ! of the precipitation region.  At the same altitudes, the turbulent
    ! sedimentation flux of the hydrometeor, < V_hm'hm' >, often increases
    ! sharply with height.  The predictive equation for < hm > contains the
    ! turbulent sedimentation term:
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz;
    ! where rho_ds is the dry, static base-state density.  When < V_hm'hm' >
    ! increases sharply with height, the rate of change of < hm > with respect
    ! to time, d<hm>/dt, becomes more negative.  When the duration of a model
    ! time step is long enough, the large, negative contribution of the explicit
    ! turbulent sedimentation term can cause the value of < hm > to become
    ! negative at the given grid level.  When this occurs, the magnitude of
    ! d<hm>/dt can be restricted (so that < hm > has a value of 0 at the end of
    ! the model time step instead of a negative value) by altering the
    ! derivative involving < V_hm'hm' >.
    !
    ! The predictive equation for < hm > is (Equation 1):
    ! d<hm>/dt = - < w > * d< hm > / dz
    !            - (1/rho_ds) * d( rho_ds * < w'hm' > ) / dz
    !            - < V_hm > * d< hm > / dz
    !            - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz
    !            + nu * d^2<hm>/dz^2
    !            + d<hm>/dt|_mc;
    !
    ! where d<hm>/dt|_mc is the rate of change of <hm> over time due to
    ! microphysics processes and w is the vertical velocity (vertical wind
    ! component).  When d<hm>/dt and < V_hm'hm' > need to be adjusted due to an
    ! excessive contribution from the turbulent sedimentation flux, the adjusted
    ! predictive equation for < hm >, resulting in the net time tendency for
    ! < hm >, d<hm>/dt|_net, is (Equation 2):
    !
    ! d<hm>/dt|_net = - < w > * d< hm > / dz
    !                 - (1/rho_ds) * d( rho_ds * < w'hm' > ) / dz
    !                 - < V_hm > * d< hm > / dz
    !                 - ( (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz )|_net
    !                 + nu * d^2<hm>/dz^2
    !                 + d<hm>/dt|_mc.
    !
    ! The net < hm > time tendency can also be written as:
    !
    ! d<hm>/dt|_net = d<hm>/dt + d<hm>/dt|_chg;
    !
    ! where d<hm>/dt|_chg is the change (or difference) between the net < hm >
    ! time tendency and the original < hm > time tendency.  This can be
    ! rewritten as:
    !
    ! d<hm>/dt|_chg = d<hm>/dt|_net - d<hm>/dt.
    !
    ! Subtracting Equation (2) from Equation (1), d<hm>/dt_chg can be written
    ! as:
    !
    ! d<hm>/dt|_chg = - ( (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz )|_net
    !                 + (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The only term on the right-hand side that can be altered is < V_hm'hm' >
    ! (rho_ds cannot be changed).  The equation can be written as:
    !
    ! d<hm>/dt|_chg = - (1/rho_ds) * ( d( rho_ds * < V_hm'hm' >|_net ) / dz
    !                                  - d( rho_ds * < V_hm'hm' > ) / dz );
    !
    ! and can also be written as (Equation 3):
    !
    ! d<hm>/dt|_chg =
    !    - (1/rho_ds) * d( rho_ds * ( < V_hm'hm' >|_net - < V_hm'hm' > ) ) / dz.
    !
    ! The above equation is discretized as follows:
    !
    ! The values of < hm > are found on the thermodynamic levels.  The left-hand
    ! side of the equation, d<hm>/dt_chg, is based on d<hm>/dt and
    ! d<hm>/dt|_net.  The original time tendency, d<hm>/dt, is discretized as:
    ! ( hm(t+1,k) - hm(t,k) ) / delta_t, where t is the time index, k is the
    ! vertical level index, and delta_t is the model time step duration.  The
    ! net time tendency, d<hm>/dt|_net, is discretized as:
    ! ( hm(t+1,k)|_net - hm(t,k) ) / delta_t.  The amount of change in the
    ! < hm > time tendency, d<hm>/dt|_chg, is discretized as:
    !
    ! ( hm(t+1,k)|_net - hm(t,k) ) / delta_t - ( hm(t+1,k) - hm(t,k) ) / delta_t
    ! = ( hm(t+1,k)|_net - hm(t,k) - hm(t+1,k) + hm(t,k) ) / delta_t
    ! = ( hm(t+1,k)|_net - hm(t+1,k) ) / delta_t.
    !
    ! The values of < V_hm'hm' > are found on the momenum levels.  Additionally,
    ! The values of rho_ds_zm are found on the momentum levels, and the values
    ! of invrs_rho_ds_zt are found on the thermodynamic levels.  On the momentum
    ! levels, the values of rho_ds_zm are multiplied by the values of
    ! < V_hm'hm' >.  The derivative of (rho_ds_zm * < V_hm'hm' >) is taken over
    ! the intermediate (central) thermodynamic level, where it is multiplied by
    ! invrs_rho_ds_zt.
    !
    ! =====rho_ds_zm(k)=====Vhmphmp(k)============================= m(k)
    !
    ! ------invrs_rho_ds_zt(k)--------d(rho_ds*Vhmphmp)/dz--------- t(k)
    !
    ! =====rho_ds_zm(k-1)===Vhmphmp(k-1)=========================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )
    !
    ! Based on the above discretization for -(1/rho_ds)*d(rho_ds*<V_hm'hm'>)/dz,
    ! the right-hand side of Equation (3) is discretized as:
    !
    ! - invrs_rho_ds_zt(k)
    !   * ( ( rho_ds_zm(k) * Vhmphmp(t,k)|_net
    !         - rho_ds_zm(k-1) * Vhmphmp(t,k-1)|_net )
    !       - ( rho_ds_zm(k) * Vhmphmp(t,k)
    !           - rho_ds_zm(k-1) * Vhmphmp(t,k-1) ) ) / delta_z;
    !
    ! where delta_z is the altitude difference between momentum levels k and
    ! k-1.  This can also can be written as:
    !
    ! - invrs_rho_ds_zt(k)
    !   * ( rho_ds_zm(k)
    !       * ( Vhmphmp(t,k)|_net - Vhmphmp(t,k) )
    !       - rho_ds_zm(k-1)
    !         * ( Vhmphmp(t,k-1)|_net - Vhmphmp(t,k-1) ) ) / delta_z.
    !
    ! The full discretized equation is:
    !
    ! ( hm(t+1,k)|_net - hm(t+1,k) ) / delta_t
    ! = - invrs_rho_ds_zt(k)
    !     * ( rho_ds_zm(k)
    !         * ( Vhmphmp(t,k)|_net - Vhmphmp(t,k) )
    !         - rho_ds_zm(k-1)
    !           * ( Vhmphmp(t,k-1)|_net - Vhmphmp(t,k-1) ) ) / delta_z.
    !
    ! The full discretized equation can also be written as (Equation 4):
    !
    ! ( hm(t+1,k)|_net - hm(t+1,k) ) / delta_t
    ! = - invrs_rho_ds_zt(k)
    !     * invrs_dzt(k)
    !     * ( rho_ds_zm(k) * ( Vhmphmp(t,k)|_net - Vhmphmp(t,k) )
    !         - rho_ds_zm(k-1) * ( Vhmphmp(t,k-1)|_net - Vhmphmp(t,k-1) ) ).
    !
    !
    ! The procedure for adjusting the values of the hydrometeor and the values
    ! of the hydrometeor turbulent sedimentation flux is as follows:
    !
    ! The values are adjusted one grid level at a time, starting with the top
    ! grid level and going downwards.  The value of < V_hm'hm' > at the highest
    ! momentum grid level (model top) is fixed with a value of 0, meaning that
    ! there is not any precipitation flux through the top of the model.  The
    ! value of adjusted < V_hm'hm' > (or < V_hm'hm' >|_net) at the top grid
    ! level is always the same.  The difference between < V_hm'hm' >|_net and
    ! the original < V_hm'hm' > is always 0 at the top of the model.
    !
    ! Working downwards, the highest thermodynamic level is below the highest
    ! momentum level.  The value of < hm > is checked to find out if it is less
    ! than 0.  In the scenario that the value of < hm > is less than 0, < hm >
    ! is reset to 0, and the value of < V_hm'hm' > at the momentum level below
    ! that thermodynamic level is reset accordingly.  When the value of < hm >
    ! is greater than or equal to 0, the value of < hm > is not reset and the
    ! value of < V_hm'hm' > at the momentum level below that thermodynamic level
    ! remains the same.  When the value of < V_hm'hm' > needs to be adjusted at
    ! a momentum level, the value of < hm > also needs to be adjusted at the
    ! thermodynamic level below the momentum level where the < V_hm'hm' >
    ! adjustment took place.  The process repeats at the next vertical level
    ! down and continues until the model lower boundary is reached.
    !
    ! When the mean value of a hydrometeor is below 0 and needs to be returned
    ! to 0 at thermodynamic level k, Equation (4) is used by setting
    ! hm(t+1,k)|_net = 0.  The value of < V_hm'hm' > at momentum level k, which
    ! is above thermodynamic level k, is held fixed, as it was already looped
    ! over and was either originally okay or was reset to an appropriate value.
    ! In Equation (4), set Vhmphmp(t,k)|_net = Vhmphmp(t,k).  The resulting
    ! equation is:
    !
    ! - hm(t+1,k) / delta_t = invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k-1)
    !                         * ( Vhmphmp(t,k-1)|_net - Vhmphmp(t,k-1) ).
    !
    ! Solving for the amount of adjustment to < V_hm'hm' > at momentum
    ! level (k-1), Vhmphmp_adj_amt(k-1), where:
    !
    ! Vhmphmp_adj_amt(k-1) = Vhmphmp(t,k-1)|_net - Vhmphmp(t,k-1);
    !
    ! the amount of necessary adjustment is:
    !
    ! Vhmphmp_adj_amt(k-1) = ( - hm(t+1,k) / delta_t )
    !                        / ( invrs_rho_ds_zt(k) * invrs_dzt(k)
    !                            * rho_ds_zm(k-1) ).
    !
    ! The new value of < V_hm'hm' > at momentum level k-1 is:
    !
    ! Vhmphmp(t,k-1)|_net = Vhmphmp(t,k-1) + Vhmphmp_adj_amt(k-1).
    !
    ! In the scenario when the value of < hm > needs to be adjusted at
    ! thermodynamic level k because the value < V_hm'hm' > was previously
    ! adjusted at momentum level k (the momentum level above thermodynamic
    ! level k), Equation 4 is used to solve for hm(t+1,k)|_net.  The values of
    ! < V_hm'hm' > are held fixed (for now) at momentum level k-1, such that
    ! Vhmphmp(t,k-1)|_net = Vhmphmp(t,k-1).  The values of Vhmphmp(t,k)|_net
    ! and Vhmphmp(t,k) are used that were calculated at the above level
    ! (previous loop iteration).  The resulting equation is:
    !
    ! ( hm(t+1,k)|_net - hm(t+1,k) ) / delta_t =
    !    - invrs_rho_ds_zt(k) * invrs_dzt(k)
    !      * rho_ds_zm(k) * ( Vhmphmp(t,k)|_net - Vhmphmp(t,k) ).
    !
    ! The new value of < hm >, hm(t+1,k)|_net, is given by the equation:
    !
    ! hm(t+1,k)|_net = hm(t+1,k)
    !                  - invrs_rho_ds_zt(k) * delta_t * invrs_dzt(k)
    !                    * rho_ds_zm(k) * ( Vhmphmp(t,k)|_net - Vhmphmp(t,k) ).
    !
    ! In some scenarios, the new value of < hm > may be pushed below 0 as a
    ! result of the adjustment.  The next adjustment will bring the value of
    ! < hm > back to 0, while adjusting the value of < V_hm'hm' > at momentum
    ! level (k-1).  The process will continue cycling down through the vertical
    ! levels, as all values of < hm > on all vertical levels are eventually
    ! checked and adjusted when appropriate.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr,    & ! Variable(s)
        zm2zt    ! Procedure(s)

    use constants_clubb, only: &
        zero    ! Constant(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type, only: &
        stat_update_var    ! Procedure(s)

    use stats_variables, only: &
        iVrrprrp_net, & ! Variable(s)
        iVNrpNrp_net, &
        zm,           &
        l_stats_samp

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type    ! Description of which hydrometeor is being solved for.

    real( kind = time_precision ), intent(in) :: &
      dt          ! Model time step duration    [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      hmm,        & ! Mean value of hydrometeor (thermo. levels)   [units vary]
      Vhmphmp,    & ! Covariance of V_hm and hm (momentum levels)  [units(m/s)]
      Vhmphmp_zt    ! Covariance of V_hm and hm (thermo. levels)   [units(m/s)]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Vhmphmp_adj_amt  ! Amount adjustment necessary to <V_hm'hm'>  [units(m/s)]

    integer :: k  ! Loop index

    integer :: iVhmphmp_net  ! Stat index


    ! Initialize stats variable iVhmphmp_net to avoid compiler warnings.
    iVhmphmp_net = 0
       
    select case( solve_type )
    case( "rrainm" )
      iVhmphmp_net = iVrrprrp_net
    case( "Nrm" )
      iVhmphmp_net = iVNrpNrp_net
    case default
      iVhmphmp_net = 0
    end select

    ! The value of < V_hm'hm' > is set to 0 at the top of the model, as
    ! precipitation doesn't flux out the top of the model.  Likewise, the new
    ! value of < V_hm'hm' > is also 0 at the top of the model.  The difference
    ! between the values is also 0.
    Vhmphmp_adj_amt(gr%nz) = zero

    ! Loop downward through all vertical levels and make adjustments
    ! as necessary.
    do k = gr%nz, 2, -1

       if ( Vhmphmp_adj_amt(k) /= zero ) then

          ! The value of < V_hm'hm' > has been adjusted at momentum level k,
          ! which is above the current thermodynamic level k, due to < hm >
          ! having a value less than 0 at thermodynamic level k+1.  The value
          ! of < hm > needs to be adjusted due to a new value of < V_hm'hm' >,
          ! which changes the turbulent sedimentation term at thermodynamic
          ! level k.
          hmm(k) = hmm(k) &
                   - invrs_rho_ds_zt(k) * real( dt, kind = core_rknd ) &
                     * gr%invrs_dzt(k) * rho_ds_zm(k) * Vhmphmp_adj_amt(k)

       endif

       if ( hmm(k) < zero ) then

          ! The value of < hm > is less than 0.  It was either less than 0
          ! originally, after the predictive equation was solved for < hm >, or
          ! else it became less than 0 due to the adjustment above, which was
          ! due to having < hm > less than 0 at thermodynamic level k+1 and
          ! having to adjust < V_hm'hm' > at momentum level k.  Adjust
          ! < V_hm'hm' > at momentum level k-1 and reset < hm > to 0 at
          ! thermodynamic level k.
          Vhmphmp_adj_amt(k-1) &
          = - hmm(k) / ( invrs_rho_ds_zt(k) * real( dt, kind = core_rknd ) &
                         * gr%invrs_dzt(k) * rho_ds_zm(k-1) )

          hmm(k) = zero

       else  ! hmm(k) >= 0

          ! The value of < hm > is greater than or equal to 0.  The value of
          ! < V_hm'hm' > at momentum level k-1 doesn't need to be adjusted.
          Vhmphmp_adj_amt(k-1) = zero
          
       endif ! hmm(k) < 0

    enddo ! k = gr%nz, 2, -1

    ! Calculate the new value of < V_hm'hm' >.
    Vhmphmp = Vhmphmp + Vhmphmp_adj_amt

    ! Output < V_hm'hm' > on thermodynamic levels.
    ! In the scenario that < V_hm'hm' > was adjusted at either momentum level
    ! adjacent to thermodynamic level k, interpolate the adjusted < V_hm'hm' >
    ! output to thermodynamic level k.  When < V_hm'hm' > hasn't been adjusted
    ! at either momentum level, don't interpolate, as the values of Vhmphmp_zt
    ! shouldn't be unnecessarily smoothed out.  They are used later for
    ! statistical values of rain rate.
    do k = 2, gr%nz, 1
       if ( Vhmphmp_adj_amt(k) /= zero .or. Vhmphmp_adj_amt(k-1) /= zero ) then
          Vhmphmp_zt(k) = zm2zt( Vhmphmp, k )
       else ! Vhmphmp_adj_amt(k) = 0 .and. Vhmphmp_adj_amt(k-1) = 0
          Vhmphmp_zt(k) = Vhmphmp_zt(k)
       endif
    enddo
    ! The value of < V_hm'hm' > at thermodynamic level 1, which is below the
    ! model lower boundary, is set equal to the value of < V_hm'hm' > at
    ! momentum level 1, which is the model lower boundary.
    Vhmphmp_zt(1) = Vhmphmp(1)

    ! Statistics
    if ( l_stats_samp ) then

       ! The orignal value of < V_hm'hm' > is stored for statistics as Vhmphmp.
       ! That was done before this subroutine was called.  The updated or
       ! adjusted value of < V_hm'hm' > is stored for statistics as Vhmphmp_net,
       ! which is done here.  The value of Vhmphmp_net overwrites Vhmphmp so
       ! that the updated value can be used in the code to produce the
       ! appropriate statisical values for precipitation flux and rain rate.
       call stat_update_var( iVhmphmp_net, Vhmphmp, zm )

     endif

 
  return

  end subroutine turb_sed_flux_limiter

!===============================================================================

end module microphys_utilities

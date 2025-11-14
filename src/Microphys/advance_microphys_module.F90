! $Id$
!===============================================================================
module advance_microphys_module

  ! Description:
  ! Advance mean precipitating hydrometeors and mean cloud droplet concentration
  ! one time step.

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: advance_microphys,   &
            get_cloud_top_level, &
            calculate_K_hm

  ! Subroutines
  private :: advance_hydrometeor, &
             advance_Ncm, &
             microphys_solve, &
             microphys_lhs, &
             microphys_rhs, &
             write_adv_micro_errors

  ! Functions
  private :: sed_centered_diff_lhs, &
             sed_upwind_diff_lhs, &
             term_turb_sed_lhs, &
             term_turb_sed_rhs

  private ! Default Scope

  contains

  !=============================================================================
  subroutine advance_microphys( gr, dt, time_current, &                    ! In
                                hydromet_dim, hm_metadata, &               ! In
                                wm_zt, wp2, &                              ! In
                                exner, rho, rho_zm, rcm, &                 ! In
                                cloud_frac, Kh_zm, Skw_zm, &               ! In
                                rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &   ! In
                                hydromet_mc, Ncm_mc, Lscale, &             ! In
                                hydromet_vel_covar_zt_impc, &              ! In
                                hydromet_vel_covar_zt_expc, &              ! In
                                clubb_params, nu_vert_res_dep, &           ! In
                                tridiag_solve_method, &                    ! In
                                fill_holes_type, &                         ! In
                                l_upwind_xm_ma, &                          ! In
                                stats_metadata, &                          ! In
                                stats_zt, stats_zm, stats_sfc, &           ! Inout
                                hydromet, hydromet_vel_zt, hydrometp2, &   ! Inout
                                K_hm, Ncm, Nc_in_cloud, rvm_mc, &          ! Inout
                                thlm_mc, err_info, &                       ! Inout
                                wphydrometp, wpNcp )                       ! Out

    ! Description:
    ! Advance mean precipitating hydrometeors and mean cloud droplet
    ! concentration one model time step, and calculate some statistics.

    ! References:
    !---------------------------------------------------------------------------

    use grid_class, only: & 
        zt2zm_api    ! Procedure(s)

    use grid_class, only: grid ! Type

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        ic_K_hm

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use constants_clubb, only: & 
        zero,        & ! Constant(s)
        Lv,          &
        rho_lw,      & 
        fstderr,     &
        sec_per_day, &
        mm_per_m

    use parameters_microphys, only: &
        l_predict_Nc,          & ! Predict cloud droplet number conc (Morrison)
        microphys_scheme,      & ! The microphysical scheme in use
        microphys_start_time     ! When to start the microphysics [s]

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error             ! Constant

    use corr_varnce_module, only: &
        hm_metadata_type

    use stats_variables, only: &
        stats_metadata_type

    use stats_type_utilities, only: & 
        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_clubb_utilities, only: &
        stats_accumulate_hydromet_api  ! Procedure(s)

    use stats_type, only: &
        stats ! Type

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    !---------------------- Input Variables ----------------------
    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) ::  & 
      dt           ! Model timestep duration         [s]

    real( kind = time_precision ), intent(in) ::  & 
      time_current   ! Current time     [s]

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      wm_zt,           & ! w wind component on thermodynamic levels  [m/s]
      exner,           & ! Exner function                            [-]
      rho,             & ! Density on thermodynamic levels           [kg/m^3]
      rcm,             & ! Mean cloud water mixing ratio             [kg/kg]
      cloud_frac,      & ! Cloud fraction                            [-]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: & 
      wp2,        & ! Variance of vertical velocity (momentum levels) [m^2/s^2]
      rho_zm,     & ! Density on momentum levels                      [kg/m^3]
      rho_ds_zm,  & ! Dry, static density on momentum levels          [kg/m^3]
      Kh_zm,      & ! Kh Eddy diffusivity on momentum grid            [m^2/s]
      Skw_zm        ! Skewness of w on momentum levels                [-]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: & 
      hydromet_mc    ! Microphysics tendency for mean hydrometeors  [units/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      Ncm_mc, & ! Microphysics tendency for Ncm                     [num/kg/s]
      Lscale    ! Length-scale                                      [m]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(in) :: &
      tridiag_solve_method, & ! Specifier for method to solve tridiagonal systems
      fill_holes_type         ! Option for which type of hole filler to use

    logical, intent(in) :: &
      l_upwind_xm_ma ! This flag determines whether we want to use an upwind differencing
                     ! approximation rather than a centered differencing for turbulent or
                     ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.

    type (stats_metadata_type), intent(in) :: &
          stats_metadata

    !---------------------- Input/Output Variables ----------------------
    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(inout) :: &
      hydromet,        & ! Hydrometeor mean, <h_m> (thermo. levels)    [units]
      hydromet_vel_zt    ! Mean hydrometeor sed. vel. on thermo. levs. [m/s]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(inout) :: &
      hydrometp2,      & ! Variance of hydrometeor (overall) (m-levs.) [units^2]
      K_hm               ! hm eddy diffusivity on momentum grid        [m^2/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water          [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp.   [K/s]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !---------------------- Output Variables ----------------------
    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(out) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)   [(m/s)units]

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      wpNcp          ! Covariance < w'N_c' > (momentum levels)   [(m/s)(num/kg)]

    !---------------------- Local Variables ----------------------
    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim) :: &
      hydromet_vel    ! Mean hydrometeor sedimentation velocity, <V_xx> [m/s]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim) :: &
      hydromet_vel_covar       ! Covariance of V_xx & x_x (m-levs)  [units(m/s)]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim) :: &
      hydromet_vel_covar_zt    ! Covariance of V_xx & x_x (t-levs)  [units(m/s)]

    ! Turbulent advection for hydrometeors -- down-gradient approximation for
    ! covariances <w'hm'> (for any hydrometeor, hm) and <w'Nc'>:
    ! <w'hm'> = - K_hm * d<hm>/dz; and
    ! <w'Nc'> = - K_Nc & d<Nc>/dz;
    ! where the coefficients of diffusion, K_hm and K_Nc are variable and depend
    ! on multiple factors.
    real( kind = core_rknd ), dimension(gr%nzm) :: &
      K_Nc    ! Coefficient of diffusion (turb. adv.) for Nc             [m^2/s]

    integer :: k, i    ! Loop indices

    integer :: &
      cloud_top_level    ! Vertical level index of cloud top    [-]

    logical, parameter :: &
      l_use_non_local_diff_fac = .false. ! Use a non-local factor for
                                         ! eddy-diffusivity applied to
                                         ! hydrometeors

    logical, parameter :: &
      l_prevent_hm_ta_above_cloud = .false. ! Set K_hm to 0 above cloud top to
                                            ! prevent turbulent advection of
                                            ! hydrometeors to those levels.
    integer :: &
      iirr, & ! Variable(s)
      iiri, &
      iiNi

    !---------------------- Begin Code ----------------------

    iirr = hm_metadata%iirr
    iiri = hm_metadata%iiri
    iiNi = hm_metadata%iiNi

    ! Initialize intent(out) variables -- covariances <w'hm'> (for any
    ! hydrometeor, hm) and <w'Nc'>.
    if ( hydromet_dim > 0 ) then
       wphydrometp = zero
    endif
    wpNcp = zero

    ! Return if there is delay between the model start time and start of the
    ! microphysics
    if ( time_current < microphys_start_time ) return

    ! Turbulent advection for hydrometeors -- down-gradient approximation for
    ! covariances <w'hm'> (for any hydrometeor, hm) and <w'Nc'>:
    ! <w'hm'> = - K_hm * d<hm>/dz; and
    ! <w'Nc'> = - K_Nc & d<Nc>/dz;
    ! where the coefficients of diffusion, K_hm and K_Nc are variable and depend
    ! on multiple factors.
    if ( hydromet_dim > 0 ) then

       ! Prevent the turbulent advection of hydrometeors to altitudes above
       ! cloud top.
       if ( l_prevent_hm_ta_above_cloud ) then

          ! Find the vertical level index of cloud top.
          cloud_top_level = get_cloud_top_level( gr%nzt, rcm, hydromet, &
                                                 hydromet_dim, hm_metadata%iiri )

       endif ! l_prevent_hm_ta_above_cloud

       ! Solve for the value of K_hm, the coefficient of diffusion for
       ! hydrometeors.
       K_hm = calculate_K_hm( gr, wp2, Kh_zm, Skw_zm, Lscale, &
                              hydromet_dim, hm_metadata%hydromet_tol, &
                              hydromet, hydrometp2, &
                              clubb_params, &
                              l_use_non_local_diff_fac )

       do i = 1, hydromet_dim, 1

          ! Prevent the turbulent advection of hydrometeors to altitudes above
          ! cloud top.
          if ( l_prevent_hm_ta_above_cloud ) then

             ! Set K_hm to 0 above cloud top.  Since K_hm is a momentum-level
             ! variable, and since momentum grid levels are located above their
             ! corresponding thermodynamic grid levels, set K_hm to 0 starting
             ! at the vertical-level index of cloud_top_level.
             if ( cloud_top_level > 1 &
                  .and. i /= iiri .and. i /= iiNi ) then
                K_hm(cloud_top_level:gr%nzm,i) = zero
             endif ! cloud_top_level > 1
             
          endif ! l_prevent_hm_ta_above_cloud

       enddo ! i = 1, hydromet_dim, 1

    endif ! hydromet_dim > 0

    if ( stats_metadata%l_stats_samp ) then
      do i = 1, hydromet_dim, 1
        do k = 1, gr%nzm, 1
          call stat_update_var_pt( stats_metadata%iK_hm(i), k, K_hm(k,i), stats_zm )
        end do
      end do ! i = hydromet_dim, 1
    end if 


    if ( l_predict_Nc ) then

       ! Solve for the value of K_Nc, the coefficient of diffusion for cloud
       ! droplet concentration.
       do k = 1, gr%nzm, 1
          K_Nc(k) = clubb_params(ic_K_hm) * Kh_zm(k)
       enddo ! k = 1, gr%nzm, 1

    endif ! l_predict_Nc

    !------------------------------------------------------------------------
    ! Advance predictive mean precipitating hydrometeors one model timestep.
    ! Loop over all hydrometeor species, apply mean advection, turbulent
    ! advection (diffusion), mean sedimentation, turbulent sedimentation, and
    ! microphysics tendency terms.
    !------------------------------------------------------------------------

    if ( hydromet_dim > 0 ) then

       call advance_hydrometeor( gr, dt, hydromet_dim, hm_metadata, &
                                 wm_zt, exner, cloud_frac, K_hm, &
                                 rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                                 hydromet_mc, hydromet_vel_covar_zt_impc, &
                                 hydromet_vel_covar_zt_expc, &
                                 nu_vert_res_dep, &
                                 l_upwind_xm_ma, &
                                 tridiag_solve_method, &
                                 fill_holes_type, &
                                 stats_metadata, &
                                 stats_zt, stats_zm, &
                                 hydromet, hydromet_vel_zt, &
                                 hydrometp2, rvm_mc, thlm_mc, err_info, &
                                 wphydrometp, hydromet_vel, &
                                 hydromet_vel_covar, hydromet_vel_covar_zt )

    endif ! hydromet_dim > 0


    if ( clubb_at_least_debug_level_api( 0 ) ) then
       if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr,*) "calling advance_hydrometeor"
          call write_adv_micro_errors( gr, dt, time_current, hydromet_dim,  & ! In
                                       wm_zt, wp2, &                          ! In
                                       exner, rho, rho_zm, rcm, &             ! In
                                       cloud_frac, Kh_zm, Skw_zm, &           ! In
                                       rho_ds_zm, rho_ds_zt, &                ! In
                                       invrs_rho_ds_zt, &                     ! In
                                       hydromet_mc, Ncm_mc, Lscale, &         ! In
                                       hydromet_vel_covar_zt_impc, &          ! In
                                       hydromet_vel_covar_zt_expc, &          ! In
                                       clubb_params, nu_vert_res_dep, &       ! In
                                       l_upwind_xm_ma, &                      ! In
                                       hydromet, hydromet_vel_zt, &           ! In
                                       hydrometp2, K_hm, Ncm, &               ! In
                                       Nc_in_cloud, rvm_mc, thlm_mc, &        ! In
                                       wphydrometp, wpNcp, err_info )         ! In
          return
       endif !  err_info%err_code == clubb_fatal_error
    endif

    !-----------------------------------------------------------------------
    ! When mean cloud droplet concentration, Ncm, is predicted, apply
    ! sedimentation, advection, and diffusion, and advance Ncm one model
    ! timestep.
    !-----------------------------------------------------------------------

    if ( l_predict_Nc ) then

       ! Nc is predicted.
       call advance_Ncm( gr, dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm, &
                         rho_ds_zt, invrs_rho_ds_zt, Ncm_mc, &
                         nu_vert_res_dep, &
                         l_upwind_xm_ma, &
                         tridiag_solve_method, &
                         stats_metadata, &
                         stats_zt, stats_zm, &
                         Ncm, Nc_in_cloud, err_info, &
                         wpNcp )

       if ( clubb_at_least_debug_level_api( 0 ) ) then
         if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
           write(fstderr,*) "in advance_Ncm"
           call write_adv_micro_errors( gr, dt, time_current, hydromet_dim, & ! In
                                        wm_zt, wp2, &                         ! In
                                        exner, rho, rho_zm, rcm, &            ! In
                                        cloud_frac, Kh_zm, Skw_zm, &          ! In
                                        rho_ds_zm, rho_ds_zt, &               ! In
                                        invrs_rho_ds_zt, &                    ! In
                                        hydromet_mc, Ncm_mc, Lscale, &        ! In
                                        hydromet_vel_covar_zt_impc, &         ! In
                                        hydromet_vel_covar_zt_expc, &         ! In
                                        clubb_params, nu_vert_res_dep, &      ! In
                                        l_upwind_xm_ma, &                     ! In
                                        hydromet, hydromet_vel_zt, &          ! In
                                        hydrometp2, K_hm, Ncm, &              ! In
                                        Nc_in_cloud, rvm_mc, thlm_mc, &       ! In
                                        wphydrometp, wpNcp, err_info )        ! In
           return
         endif
       endif

    else

       ! Nc is prescribed.
       ! The in-cloud mean of cloud droplet concentration, Nc_in_cloud, is
       ! constant in this scenario.  The overall mean of cloud droplet
       ! concentration, Ncm, depends of the in-cloud mean and cloud fraction.

       Ncm = Nc_in_cloud * cloud_frac

    endif ! l_predict_Nc

    if ( stats_metadata%l_stats_samp ) then
      call stat_update_var( stats_metadata%iNcm, Ncm, stats_zt )
      call stat_update_var( stats_metadata%iNc_in_cloud, Nc_in_cloud, stats_zt )
    endif ! stats_metadata%l_stats_samp

    if ( stats_metadata%l_stats_samp .and. iirr > 0 ) then

      ! Rainfall rate (mm/day) is defined on thermodynamic levels.  The rainfall
      ! rate is given by the equation:
      ! Rainfall rate (mm/day) = - V_rr * r_r * ( rho_a / rho_lw )
      !                            * ( 86400 s/day ) * ( 1000 mm/m ).
      ! The level average rainfall rate is given by:
      ! < Rainfall rate (mm/day) > = - < V_rr * r_r > * ( rho_a / rho_lw )
      !                                * ( 86400 s/day ) * ( 1000 mm/m );
      ! which can also be written as:
      ! < Rainfall rate (mm/day) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
      !                                * ( rho_a / rho_lw )
      !                                * ( 86400 s/day ) * ( 1000 mm/m ).
      ! Rainfall rate is defined as positive.  Since V_rr is always negative,
      ! the minus (-) sign is necessary.
      call stat_update_var( stats_metadata%iprecip_rate_zt,  & 
                            max( - ( hydromet(:,iirr) &
                                     * hydromet_vel_zt(:,iirr) &
                                     + hydromet_vel_covar_zt(:,iirr) ), &
                                 zero ) &
                            * ( rho / rho_lw ) & 
                            * sec_per_day &
                            * mm_per_m, stats_zt )

      ! Precipitation Flux (W/m^2) is defined on momentum levels.  The
      ! precipitation flux is given by the equation:
      ! Precip. flux (W/m^2) = - V_rr * r_r * rho_a * L_v.
      ! The level average precipitation flux is given by:
      ! < Precip. flux (W/m^2) > = - < V_rr * r_r > * rho_a * L_v;
      ! which can also be written as:
      ! < Precip. flux (W/m^2) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
      !                              * rho_a * L_v.
      ! It is generally a convention in meteorology to show Precipitation Flux
      ! as a positive downward quantity, so the minus (-) sign is necessary.
      call stat_update_var( stats_metadata%iFprec,  & 
                            max( - ( zt2zm_api( gr, hydromet(:,iirr) )  & 
                                     * hydromet_vel(:,iirr) &
                                     + hydromet_vel_covar(:,iirr) ), &
                                 zero ) &
                            * rho_zm * Lv, stats_zm )

      ! Store values of surface fluxes for statistics
      ! See notes above.

      if ( trim( microphys_scheme ) /= "morrison" ) then
        call stat_update_var_pt( stats_metadata%iprecip_rate_sfc, 1,  & 
                                 max( - ( hydromet(1,iirr) &
                                          * hydromet_vel_zt(1,iirr) &
                                          + hydromet_vel_covar_zt(1,iirr) ), &
                                      zero ) &
                                  * ( rho(1) / rho_lw ) & 
                                  * sec_per_day &
                                  * mm_per_m, stats_sfc )
      endif ! microphys_scheme /= "morrison"

      call stat_update_var_pt( stats_metadata%irain_flux_sfc, 1, & 
                               max( - ( zt2zm_api( gr, hydromet(:,iirr), 1 )  & 
                                        * hydromet_vel(1,iirr) &
                                        + hydromet_vel_covar(1,iirr) ), &
                                    zero ) &
                               * rho_zm(1) * Lv, stats_sfc )

      ! Also store the value of surface rain water mixing ratio.
      call stat_update_var_pt( stats_metadata%irrm_sfc, 1,  & 
                               ( zt2zm_api( gr, hydromet(:,iirr), 1 ) ), stats_sfc )

    endif ! stats_metadata%l_stats_samp and iirr > 0

    call stats_accumulate_hydromet_api( gr, hydromet_dim, hm_metadata, &
                                        hydromet, rho_ds_zt, &
                                        stats_metadata, &
                                        stats_zt, stats_sfc )

    call write_adv_micro_errors( gr, dt, time_current, hydromet_dim, & ! In
                                 wm_zt, wp2, &                         ! In
                                 exner, rho, rho_zm, rcm, &            ! In
                                 cloud_frac, Kh_zm, Skw_zm, &          ! In
                                 rho_ds_zm, rho_ds_zt, &               ! In
                                 invrs_rho_ds_zt, &                    ! In
                                 hydromet_mc, Ncm_mc, Lscale, &        ! In
                                 hydromet_vel_covar_zt_impc, &         ! In
                                 hydromet_vel_covar_zt_expc, &         ! In
                                 clubb_params, nu_vert_res_dep, &      ! In
                                 l_upwind_xm_ma, &                     ! In
                                 hydromet, hydromet_vel_zt, &          ! In
                                 hydrometp2, K_hm, Ncm, &              ! In
                                 Nc_in_cloud, rvm_mc, thlm_mc, &       ! In
                                 wphydrometp, wpNcp, err_info )        ! In

    return

  end subroutine advance_microphys

  !=============================================================================
  subroutine advance_hydrometeor( gr, dt, hydromet_dim, hm_metadata, &
                                  wm_zt, exner, cloud_frac, K_hm, &
                                  rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                                  hydromet_mc, hydromet_vel_covar_zt_impc, &
                                  hydromet_vel_covar_zt_expc, &
                                  nu_vert_res_dep, &
                                  l_upwind_xm_ma, &
                                  tridiag_solve_method, &
                                  fill_holes_type, &
                                  stats_metadata, &
                                  stats_zt, stats_zm, &
                                  hydromet, hydromet_vel_zt, &
                                  hydrometp2, rvm_mc, thlm_mc, err_info, &
                                  wphydrometp, hydromet_vel, &
                                  hydromet_vel_covar, hydromet_vel_covar_zt )

    ! Description:
    !   Advance each hydrometeor (precipitating hydrometeor) one model time step.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        zt2zm_api    ! Procedure(s)

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: & 
        one_half,       & ! Constant(s)
        zero,           &
        zero_threshold, &
        fstderr

    use advance_helper_module, only : &
        calc_xpwp  ! Procedure(s)

    use fill_holes, only: &
        fill_holes_driver_api,   & ! Procedure(s)
        setup_stats_indices

    use parameters_tunable, only: & 
        nu_vertical_res_dep  ! Type(s)

    use parameters_microphys, only: &
        l_hydromet_sed,    & ! Variable(s)
        l_upwind_diff_sed

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error             ! Constant

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_update_var,      & ! Procedure(s)
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use stats_variables, only: &
        stats_metadata_type

    use corr_varnce_module, only: &
        hm_metadata_type

    use stats_type, only: stats ! Type

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    !---------------------------- Input Variables ----------------------------
    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) ::  & 
      dt    ! Duration of one model time step         [s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      wm_zt,      & ! mean w wind component on thermodynamic levels  [m/s]
      exner,      & ! Exner function, (p/p1000mb)^(Rd/Cp)            [-]
      cloud_frac    ! Cloud fraction                                 [-]

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(in) :: &
      K_hm,      & ! Coefficient of diffusion (turb. adv.) for hydrometeors [m^2/s]
      rho_ds_zm    ! Dry, static density on momentum levels   [kg/m^3]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: & 
      hydromet_mc     ! Change in hydrometeors due to microphysics  [units/s]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values


    logical, intent(in) :: &
      l_upwind_xm_ma ! This flag determines whether we want to use an upwind differencing
                     ! approximation rather than a centered differencing for turbulent or
                     ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.

    integer, intent(in) :: &
      tridiag_solve_method, & ! Specifier for method to solve tridiagonal systems
      fill_holes_type         ! Specifier for which hole filling method to use

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !---------------------------- Input/Output Variables ----------------------------
    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(inout) :: &
      hydromet,        & ! Hydrometeor mean, <h_m> (thermodynamic levs.) [units]
      hydromet_vel_zt    ! Mean hydrometeor sed. velocity on thermo. levs. [m/s]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(inout) :: &
      hydrometp2         ! Variance of hydrometeor (overall) (m-levs.) [units^2]

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water          [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp.   [K/s]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !---------------------------- Output Variables ----------------------------
    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(out) :: &
      wphydrometp,  & ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]
      hydromet_vel    ! Mean hydrometeor sedimentation velocity, <V_hm>    [m/s]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(out) :: &
      hydromet_vel_covar       ! Covariance of V_hm & h_m (m-levs)  [units(m/s)]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(out) :: &
      hydromet_vel_covar_zt    ! Covariance of V_hm & h_m (t-levs)  [units(m/s)]

    !---------------------------- Local Variables ----------------------------
    real( kind = core_rknd ), dimension(3,1,gr%nzt) :: & 
      lhs_ta    ! LHS corresponding to contribution from turbulent adv.  [1/s]

    real( kind = core_rknd ), dimension(3,gr%nzt) :: &
      sed_diff_lhs,       & ! Variables use to save stats results
      sed_turb_lhs
      
    real( kind = core_rknd ), dimension(3,gr%nzt) :: & 
      lhs_ma    ! LHS corresponding to contribution from mean advection  [1/s]

    real( kind = core_rknd ), dimension(3,gr%nzt) :: & 
      lhs    ! Left hand side array

    real( kind = core_rknd ), dimension(gr%nzt) :: & 
      rhs    ! Right hand side vector

    character(len=10) :: hydromet_name

    real( kind = core_rknd ) :: &
      max_velocity    ! Maximum sedimentation velocity         [m/s]

    real( kind = core_rknd ), dimension(gr%nzm) :: &
      hydromet_zm     ! Mean of hydrometeor interp. to m-levs. [units vary]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim) :: &
      ratio_hmp2_on_hmm2    ! Value of <hm'^2> / <hm>^2    [-]
      
    real( kind = core_rknd ), dimension(gr%nzm) :: &
      xpwp,       & ! x'w' for arbitrary x
      K_hm_nu_hm    ! K_hm + nu_vert_res_dep%nu_hm

    logical :: l_fill_holes_hm = .true.

    integer :: i, k    ! Loop indices

    ! Stat indices
    integer :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
               ixrm_bt, ixrm_mc

    !---------------------------- Begin Code ----------------------------

    ! Loop over each type of precipitating hydrometeor and advance one model
    ! time step.
    do i = 1, hydromet_dim

       ! Set up the stats indices for hydrometeor index i
       call setup_stats_indices( i, stats_metadata, hydromet_dim, & ! Intent(in)
                                 hm_metadata%hydromet_list,    & ! Intent(in)
                                 ixrm_bt, ixrm_hf, ixrm_wvhf,     & ! Intent(inout)
                                 ixrm_cl, ixrm_mc,                & ! Intent(inout)
                                 max_velocity )                     ! Intent(inout)

       if ( stats_metadata%l_stats_samp ) then

          ! Update explicit contributions to the hydrometeor species
          call stat_update_var( ixrm_mc, hydromet_mc(:,i), stats_zt )

          ! Save prior value of the hydrometeors for determining total time
          ! tendency.
          if ( l_hydromet_sed(i) .and. l_upwind_diff_sed ) then
             call stat_begin_update( gr%nzt, ixrm_bt, &
                                     hydromet(:,i) &
                                     / dt, stats_zt )
          else
             do k = 1, gr%nzt, 1
                call stat_begin_update_pt( ixrm_bt, k, &
                                           hydromet(k,i) &
                                           / dt, stats_zt )
             enddo
          endif

       endif ! stats_metadata%l_stats_samp

       ! Set realistic limits on sedimentation velocities, following the
       ! numbers in the Morrison microphysics.
       do k = 1, gr%nzt
          if ( clubb_at_least_debug_level_api( 1 ) ) then
            ! Print a warning if the velocity has a large magnitude or the
            ! velocity is in the wrong direction.
             if ( hydromet_vel_zt(k,i) < max_velocity .or. &
                  hydromet_vel_zt(k,i) > zero_threshold ) then

                write(fstderr,*) trim( hm_metadata%hydromet_list(i) )// &
                                 " velocity at k = ", k, " = ", &
                                 hydromet_vel_zt(k,i), "m/s"
             endif
          endif
          hydromet_vel_zt(k,i) &
          = min( max( hydromet_vel_zt(k,i), max_velocity ), zero_threshold )
       enddo ! k = 1, gr%nzt, 1

       ! Interpolate velocity to the momentum grid for a centered difference
       ! approximation of the sedimenation term.
       hydromet_vel(:,i) = zt2zm_api( gr, hydromet_vel_zt(:,i) )
       hydromet_vel(gr%nzm,i) = zero ! Upper boundary condition

       ! Calculate the value of <hm'^2> / <hm>^2.  This will be used to update
       ! <hm'^2> after <hm> has been advanced one model timestep.  This method
       ! is being used because CLUBB does not currently have a predictive
       ! equation for <hm'^2> (hydrometp2).
       hydromet_zm = zt2zm_api( gr, hydromet(:,i) )

       do k = 1, gr%nzm, 1

          if ( hydromet_zm(k) > hm_metadata%hydromet_tol(i) ) then

             ! Calculate the ratio of the overall variance of the hydrometeor
             ! to the overall mean of the hydrometeor squared.
             ratio_hmp2_on_hmm2(k,i) = hydrometp2(k,i) / hydromet_zm(k)**2

          else  ! hydromet_zm(k) <= hydromet_tol(i)

             ! The overall mean of the hydrometeor is 0 or is treated as 0.  The
             ! overall variance of the hydrometeor is also 0 in this scenario.
             ! The ratio is undefined, but will be assigned a value of 0.
             ratio_hmp2_on_hmm2(k,i) = zero

          endif ! hydromet_zm(k) > hydromet_tol(i)

       enddo ! k = 1, gr%nzm, 1

       ! Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
       ! a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
       ! A Crank-Nicholson time-stepping scheme is used for this variable.
       ! This is the portion of the calculation using < h_m > from timestep t. 
       K_hm_nu_hm(:) = K_hm(:,i) + nu_vert_res_dep%nu_hm(1)
       
       call calc_xpwp( gr, K_hm_nu_hm, hydromet(:,i), &
                       xpwp )
        
       wphydrometp(2:gr%nzm-1,i) = - one_half * xpwp(2:gr%nzm-1)

       ! A zero-flux boundary condition is used for hydrometeors.
       wphydrometp(1,i) = zero
       wphydrometp(gr%nzm,i) = zero

       ! Add implicit terms to the LHS matrix
       call microphys_lhs( gr, trim( hm_metadata%hydromet_list(i) ), l_hydromet_sed(i), & ! In
                           dt, K_hm(:,i), nu_vert_res_dep%nu_hm(1), wm_zt,  & ! In
                           hydromet_vel(:,i), hydromet_vel_zt(:,i),         & ! In
                           hydromet_vel_covar_zt_impc(:,i),                 & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,           & ! In
                           l_upwind_xm_ma,                                  & ! In
                           stats_metadata,                                  & ! In
                           lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,      & ! Out
                           lhs )                                              ! Out

       ! Set up explicit term in the RHS vector
       call microphys_rhs( gr, trim( hm_metadata%hydromet_list(i) ), dt, l_hydromet_sed(i), & ! In
                           hydromet(:,i), hydromet_mc(:,i),                     & ! In
                           K_hm(:,i), nu_vert_res_dep%nu_hm(1), cloud_frac,     & ! In
                           hydromet_vel_covar_zt_expc(:,i),                     & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,               & ! In
                           stats_metadata,                                      & ! In
                           stats_zt,                                            & ! In
                           rhs )                                                  ! Out

       !!!!! Advance hydrometeor one time step.
       call microphys_solve( gr, trim( hm_metadata%hydromet_list(i) ), l_hydromet_sed(i), & ! In
                             lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,      & ! In
                             cloud_frac,                                      & ! In
                             tridiag_solve_method,                            & ! In
                             stats_metadata,                                  & ! In
                             stats_zt,                                        & ! InOut
                             lhs, rhs, hydromet(:,i), err_info )                ! InOut

       if ( clubb_at_least_debug_level_api( 0 ) ) then
           if ( any(err_info%err_code == clubb_fatal_error) ) then
                write(fstderr,*) err_info%err_header_global
                write(fstderr,*) "Error in hydrometeor field " &
                                  // trim( hm_metadata%hydromet_list(i) )

                write(fstderr,*) trim( hm_metadata%hydromet_list(i) ) // " = ", hydromet(:,i)
     
           endif !  err_info%err_code == clubb_fatal_error 
        end if

    enddo ! i = 1, hydromet_dim, 1

    ! Now that all precipitating hydrometeors have been advanced, fill holes in
    ! hydromet profiles.
    call fill_holes_driver_api( gr, gr%nzt, dt, hydromet_dim,         & ! Intent(in)
                                hm_metadata, l_fill_holes_hm,         & ! Intent(in)
                                rho_ds_zt, exner,                     & ! Intent(in)
                                fill_holes_type,                      & ! Intent(in)
                                stats_metadata,                       & ! Intent(in)
                                stats_zt,                             & ! intent(inout)
                                thlm_mc, rvm_mc, hydromet )             ! Intent(inout)

    ! Loop over each type of precipitating hydrometeor and calculate hydrometeor
    ! covariances (<w'hm'> and <V_hm'hm'>) and other quantities requiring the
    ! value of hydromet (<hm>) from the (t+1) timestep.
    do i = 1, hydromet_dim

       ! Set up the stats indices for hydrometeor at index i
       call setup_stats_indices( i, stats_metadata, hydromet_dim, & ! Intent(in)
                                 hm_metadata%hydromet_list,       & ! Intent(in)
                                 ixrm_bt, ixrm_hf, ixrm_wvhf,     & ! Intent(inout)
                                 ixrm_cl, ixrm_mc,                & ! Intent(inout)
                                 max_velocity )                     ! Intent(inout)

       ! Print warning message if any hydrometeor species has a value < 0.
       if ( any( hydromet(:,i) < zero_threshold ) ) then

          hydromet_name = hm_metadata%hydromet_list(i)

          if ( clubb_at_least_debug_level_api( 1 ) ) then
             do k = 1, gr%nzt
                if ( hydromet(k,i) < zero_threshold ) then
                   write(fstderr,*) trim( hydromet_name ) //" < ", &
                                    zero_threshold, &
                                    " in advance_microphys at k= ", k
                endif
             enddo
          endif

       endif ! hydromet(:,i) < 0

       ! Lower boundary condition
       ! Hydrometeors that are below the model lower boundary level have
       ! sedimented out of the model domain, and is not conserved.
       if ( hydromet(1,i) < hm_metadata%hydromet_tol(i) ) then
          hydromet(1,i) = zero_threshold
       endif

       ! Calculate the value of <hm'^2> (hydrometp2).  This is based on the
       ! ratio of <hm'^2> / <hm>^2 (ratio_hmp2_on_hmm2) that was saved before
       ! hydrometeors were updated.  This method is being used because CLUBB
       ! does not currently have a predictive equation for <hm'^2>.
       hydromet_zm = max( zt2zm_api( gr, hydromet(:,i) ), 0.0_core_rknd )

       do k = 1, gr%nzm, 1
          hydrometp2(k,i) = ratio_hmp2_on_hmm2(k,i) * hydromet_zm(k)**2
       enddo ! k = 1, gr%nz, 1

       ! Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
       ! a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
       ! A Crank-Nicholson time-stepping scheme is used for this variable.
       ! This is the portion of the calculation using < h_m > from timestep t+1.
       K_hm_nu_hm(:) = K_hm(:,i) + nu_vert_res_dep%nu_hm(1)
       
       call calc_xpwp( gr, K_hm_nu_hm, hydromet(:,i), &
                       xpwp )
        
       wphydrometp(2:gr%nzm-1,i) = - one_half * xpwp(2:gr%nzm-1)
       
       ! A zero-flux boundary condition is used for hydrometeors.
       wphydrometp(1,i) = zero
       wphydrometp(gr%nzm,i) = zero

       !!! Calculate the covariance of hydrometeor sedimentation velocity and
       !!! the hydrometeor, which is solved semi-implicitly on thermodynamic
       !!! levels.
       hydromet_vel_covar_zt(:,i) &
       = hydromet_vel_covar_zt_impc(:,i) * hydromet(:,i) &
         + hydromet_vel_covar_zt_expc(:,i)

       !!! Calculate the covariance of hydrometeor sedimentation velocity and
       !!! the hydrometeor, < V_hm'h_m' >, by interpolating the thermodynamic
       !!! level results to momentum levels.
       hydromet_vel_covar(:,i) = zt2zm_api( gr, hydromet_vel_covar_zt(:,i) )

       ! Boundary conditions for < V_hm'hm' >.
       hydromet_vel_covar(gr%nzm,i) = zero

       ! Statistics for all covariances involving hydrometeors:  < w'h_m' >,
       ! <V_rr'r_r'>, and <V_Nr'N_r'>.
       if ( stats_metadata%l_stats_samp ) then

          if ( stats_metadata%ihydrometp2(i) > 0 ) then

             ! Covariance of vertical velocity and the hydrometeor.
             call stat_update_var( stats_metadata%ihydrometp2(i), &
                                   hydrometp2(:,i), stats_zm )

          endif

          if ( stats_metadata%iwphydrometp(i) > 0 ) then

             ! Covariance of vertical velocity and the hydrometeor.
             call stat_update_var( stats_metadata%iwphydrometp(i), &
                                   wphydrometp(:,i), stats_zm )

          endif

          if ( trim( hm_metadata%hydromet_list(i) ) == "rrm" .and. &
                     stats_metadata%iVrrprrp > 0 ) then

             ! Covariance of sedimentation velocity of r_r and r_r.
             call stat_update_var( stats_metadata%iVrrprrp, &
                                   hydromet_vel_covar(:,hm_metadata%iirr), stats_zm )

          elseif ( trim( hm_metadata%hydromet_list(i) ) == "Nrm" .and. &
                   stats_metadata%iVNrpNrp > 0 ) then

             ! Covariance of sedimentation velocity of N_r and N_r.
             call stat_update_var( stats_metadata%iVNrpNrp, &
                                   hydromet_vel_covar(:,hm_metadata%iiNr), stats_zm )

          endif

       endif ! stats_metadata%l_stats_samp

       if ( stats_metadata%l_stats_samp ) then

          ! Total time tendency
          if ( l_hydromet_sed(i) .and. l_upwind_diff_sed ) then
             call stat_end_update( gr%nzt, ixrm_bt, &
                                   hydromet(:,i) &
                                   / dt, stats_zt )
          else
             do k = 1, gr%nzt, 1
                call stat_end_update_pt( ixrm_bt, k, &
                                         hydromet(k,i) &
                                         / dt, stats_zt )
             enddo
          endif

       endif ! stats_metadata%l_stats_samp

    enddo ! i = 1, hydromet_dim, 1


    return

  end subroutine advance_hydrometeor

  !=============================================================================
  subroutine advance_Ncm( gr, dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm, &
                          rho_ds_zt, invrs_rho_ds_zt, Ncm_mc, &
                          nu_vert_res_dep, &
                          l_upwind_xm_ma, &
                          tridiag_solve_method, &
                          stats_metadata, &
                          stats_zt, stats_zm, &
                          Ncm, Nc_in_cloud, err_info, &
                          wpNcp )

    ! Description:
    ! Advance cloud droplet concentration (Ncm) one model time step.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        zt2zm_api    ! Procedure(s)

    use grid_class, only: grid ! Type

    use constants_clubb, only: &
        pi,              & ! Constant(s)
        four_thirds,     &
        one,             &
        one_half,        &
        zero,            &
        rho_lw,          &
        mvr_cloud_max,   &
        cloud_frac_min,  &
        Nc_in_cloud_min, &
        fstderr

    use advance_helper_module, only : &
        calc_xpwp  ! Procedure(s)

    use parameters_tunable, only: & 
        nu_vertical_res_dep  ! Type(s)

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff  ! Use in cloud values of Nc for diffusion

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error                 ! Constant

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var,      & ! Procedure(s)
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use stats_variables, only: &
        stats_metadata_type

    use stats_type, only: stats ! Type

    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    !------------------------------- Input Variables -------------------------------
    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) ::  & 
      dt    ! Duration of one model time step         [s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      wm_zt,      & ! mean w wind component on thermodynamic levels  [m/s]
      cloud_frac, & ! Cloud fraction                                 [-]
      rcm           ! Mean cloud water mixing ratio                  [kg/kg]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: & 
      K_Nc,       & ! Coefficient of diffusion (turb. adv.) for Nc   [m^2/s]
      rho_ds_zm     ! Dry, static density on momentum levels   [kg/m^3]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      Ncm_mc     ! Change in Ncm due to microphysics  [num/kg/s]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    logical, intent(in) :: &
      l_upwind_xm_ma ! This flag determines whether we want to use an upwind differencing
                     ! approximation rather than a centered differencing for turbulent or
                     ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.

    integer, intent(in) :: &
      tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------------- InOut Variables -------------------------------
    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !------------------------------- Output Variables -------------------------------
    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      wpNcp    ! Covariance < w'N_c' > (momentum levels)    [(m/s)(num/kg)]

    !------------------------------- Locak Variables -------------------------------
    real( kind = core_rknd ), dimension(gr%nzt) :: &
      Ncm_vel_covar_zt_impc, & ! Imp. comp. <V_Nc'N_c'> t-levs [m/s]
      Ncm_vel_covar_zt_expc    ! Exp. comp. <V_Nc'N_c'> t-levs [(num/kg)(m/s)]

    real( kind = core_rknd ), dimension(gr%nzm) :: &
      Ncm_vel       ! Mean N_c sedimentation velocity, <V_Nc> (m-levs.) [m/s]

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      Ncm_vel_zt    ! Mean N_c sedimenation velocity on thermo. levs.   [m/s]

    real( kind = core_rknd ), dimension(3,gr%nzt) :: & 
      lhs    ! Left hand side array

    real( kind = core_rknd ), dimension(gr%nzt) :: & 
      rhs    ! Right hand side vector

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      Ncm_mvr_min, & ! Min. allowable Ncm due to max. cloud droplet mvr [num/kg]
      Ncm_min,     & ! Minimum allowable mean cloud droplet conc.       [num/kg]
      Ncic_min       ! Min. allowable in-cloud mean cloud droplet conc. [num/kg]

    real( kind = core_rknd ), dimension(gr%nzm) :: &
      xpwp,       & ! x'w' for arbitrary x
      K_Nc_nu_hm    ! K_Nc + nu_vert_res_dep%nu_hm

    real( kind = core_rknd ), dimension(3,1,gr%nzt) :: & 
      lhs_ta    ! LHS corresponding to contribution from turbulent adv.  [1/s]

    real( kind = core_rknd ), dimension(3,gr%nzt) :: &
      sed_diff_lhs,       & ! Variables use to save stats results
      sed_turb_lhs
      
    real( kind = core_rknd ), dimension(3,gr%nzt) :: & 
      lhs_ma    ! LHS corresponding to contribution from mean advection  [1/s]

    ! Sedimentation velocity of cloud droplet concentration is not considered in
    ! the solution of N_c.
    logical, parameter :: &
      l_Ncm_sed = .false. ! Flag to sediment Ncm with sedimentation vel. Ncm_vel

    integer :: k    ! Loop index

    !------------------------------- Begin Code -------------------------------

    ! The mean sedimentation velocity of cloud droplet concentration, < V_Nc >,
    ! and the covariance of N_c sedimentation velocity with N_c, < V_Nc'Nc' >,
    ! are not used in the advancement of < N_c >.  Sedimentation velocity of N_c
    ! is not considered in the solution.  However, these variables need to be
    ! initialized to 0 and passed into microphys_lhs and microphys_rhs.
    Ncm_vel_covar_zt_impc = zero
    Ncm_vel_covar_zt_expc = zero
    Ncm_vel    = zero
    Ncm_vel_zt = zero

    if ( stats_metadata%l_stats_samp ) then

       ! Update explicit contributions to cloud droplet concentration.
       call stat_update_var( stats_metadata%iNcm_mc, &
                             Ncm_mc, stats_zt )

       ! Save prior value of Ncm for determining total time tendency.
       do k = 1, gr%nzt, 1
          if ( l_in_cloud_Nc_diff ) then
             call stat_begin_update_pt( stats_metadata%iNcm_bt, k, &
                                        ( Nc_in_cloud(k) &
                                         * max(cloud_frac(k),cloud_frac_min) ) &
                                        / dt, stats_zt )
          else
             call stat_begin_update_pt( stats_metadata%iNcm_bt, k, &
                                        Ncm(k) / dt, &
                                        stats_zt )
          endif
       enddo

    endif ! stats_metadata%l_stats_samp

    ! Solve for < w'N_c' > at all intermediate (momentum) grid levels, using
    ! a down-gradient approximation:  < w'N_c' > = - K * d< N_c >/dz.
    ! A Crank-Nicholson time-stepping scheme is used for this variable.
    ! This is the portion of the calculation using < N_c > from timestep t.
    K_Nc_nu_hm(:) = K_Nc(:) + nu_vert_res_dep%nu_hm(1)
    
    call calc_xpwp( gr, K_Nc_nu_hm, Ncm, &
                    xpwp )
     
    wpNcp(2:gr%nzm-1) = - one_half * xpwp(2:gr%nzm-1)

    ! A zero-flux boundary condition is used for N_c.
    wpNcp(1) = zero
    wpNcp(gr%nzm) = zero

    ! Add implicit terms to the LHS array
    call microphys_lhs( gr, "Ncm", l_Ncm_sed, & ! In
                        dt, K_Nc, nu_vert_res_dep%nu_hm(1), wm_zt, &  ! In
                        Ncm_vel, Ncm_vel_zt, & ! In
                        Ncm_vel_covar_zt_impc, & ! In
                        rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, & ! In
                        l_upwind_xm_ma, & ! In
                        stats_metadata, & ! In
                        lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, &
                        lhs ) ! Out

    ! Set up explicit term in the RHS vector
    if ( l_in_cloud_Nc_diff ) then

       call microphys_rhs( gr, "Ncm", dt, l_Ncm_sed, &
                           Nc_in_cloud, &
                           Ncm_mc / max( cloud_frac, cloud_frac_min ), &
                           K_Nc, nu_vert_res_dep%nu_hm(1), cloud_frac, &
                           Ncm_vel_covar_zt_expc, &
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                           stats_metadata, &
                           stats_zt, &
                           rhs )

    else

       call microphys_rhs( gr, "Ncm", dt, l_Ncm_sed, &
                           Ncm, Ncm_mc, &
                           K_Nc, nu_vert_res_dep%nu_hm(1), cloud_frac, &
                           Ncm_vel_covar_zt_expc, &
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                           stats_metadata, &
                           stats_zt, &
                           rhs )

    endif


    !!!!! Advance Ncm one time step.
    if ( l_in_cloud_Nc_diff ) then

       call microphys_solve( gr, "Ncm", l_Ncm_sed,                       & ! In
                             lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, & ! In
                             cloud_frac,                                 & ! In
                             tridiag_solve_method,                       & ! In
                             stats_metadata,                             & ! In
                             stats_zt,                                   & ! InOut
                             lhs, rhs, Nc_in_cloud, err_info )             ! InOut

       Ncm = Nc_in_cloud * max( cloud_frac, cloud_frac_min )

    else

       call microphys_solve( gr, "Ncm", l_Ncm_sed,                       & ! In
                             lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, & ! In
                             cloud_frac,                                 & ! In
                             tridiag_solve_method,                       & ! In
                             stats_metadata,                             & ! In
                             stats_zt,                                   & ! InOut
                             lhs, rhs, Ncm, err_info )                     ! InOut

       Nc_in_cloud = Ncm / max( cloud_frac, cloud_frac_min )

    endif

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "calling microphys_solve in advance_Ncm"
        return
      endif
    endif

    ! Clipping for mean cloud droplet concentration, <Nc>.

    ! When mean cloud water mixing ratio, <rc>, is found at a grid level, mean
    ! cloud droplet concentration must be at least a minimum value so that
    ! average cloud droplet mean volume radius stays within an upper bound.

    ! The minimum mean rain drop concentration is given by:
    !
    ! <Nc> = <rc> / ( (4/3) * pi * rho_lw * mvr_cloud_max^3 ).

    Ncm_mvr_min &
    = ( one / ( four_thirds * pi * rho_lw * mvr_cloud_max**3 ) ) * rcm

    ! Calculate the minimum threshold value for mean cloud droplet
    ! concentration.  This is the greater of the minimum allowable Ncm based on
    ! the maximum cloud droplet mean volume radius and Nc_in_cloud_min
    ! multiplied by max( cloud_frac, cloud_frac_min ).
    Ncm_min &
    = max( Nc_in_cloud_min * max( cloud_frac, cloud_frac_min ), Ncm_mvr_min )

    ! Print warning message if Ncm (at any level) has a value < Ncm_min.
    if ( any( Ncm < Ncm_min ) ) then

       if ( clubb_at_least_debug_level_api( 1 ) ) then
          do k = 1, gr%nzt
             if ( Ncm(k) < Ncm_min(k) ) then
                write(fstderr,*) "Ncm < ", Ncm_min(k), &
                                 " in advance_microphys at k = ", k
             endif ! Ncm(k) < Ncm_min(k)
          enddo ! k = 1, gr%nzt, 1
       endif ! clubb_at_least_debug_level_api( 1 )

    endif ! Ncm < Ncm_min

    ! Store the previous value of Ncm for the effect of clipping.
    if ( stats_metadata%l_stats_samp ) then
       call stat_begin_update( gr%nzt, stats_metadata%iNcm_cl, &
                               Ncm / dt, stats_zt )
    endif

    if ( any( Ncm < Ncm_min ) ) then

       ! Clip any values of Ncm below the minimum threshold value of Ncm_min to
       ! the minimum threshold value.
       where ( Ncm < Ncm_min )
          Ncm = Ncm_min
       end where

    endif ! Ncm < Ncm_min

    ! Enter the new value of Ncm for the effect of clipping.
    if ( stats_metadata%l_stats_samp ) then
       call stat_end_update( gr%nzt, stats_metadata%iNcm_cl, &
                             Ncm / dt, stats_zt )
    endif

    ! Calculate the minimum threshold value for in-cloud mean cloud droplet
    ! concentration.  This is the greater of the minimum allowable Nc_in_cloud
    ! based on the maximum cloud droplet mean volume radius (Ncm_mvr_min)
    ! divided by max( cloud_frac, cloud_frac_min ) and Nc_in_cloud_min.
    Ncic_min &
    = max( Nc_in_cloud_min, Ncm_mvr_min /  max( cloud_frac, cloud_frac_min ) )

    ! Clip any negative values of Nc_in_cloud to the minimum threshold value of
    ! Ncic_min.
    where ( Nc_in_cloud < Ncic_min )
       Nc_in_cloud = Ncic_min
    end where

    ! Solve for < w'N_c' > at all intermediate (momentum) grid levels, using
    ! a down-gradient approximation:  < w'N_c' > = - K * d< N_c >/dz.
    ! A Crank-Nicholson time-stepping scheme is used for this variable.
    ! This is the portion of the calculation using < N_c > from timestep t+1. 
    K_Nc_nu_hm(:) = K_Nc(:) + nu_vert_res_dep%nu_hm(1)
    
    call calc_xpwp( gr, K_Nc_nu_hm, Ncm, &
                    xpwp )
     
    wpNcp(2:gr%nzm-1) = - one_half * xpwp(2:gr%nzm-1)

    ! A zero-flux boundary condition is used for N_c.
    wpNcp(1) = zero
    wpNcp(gr%nzm) = zero

    ! Statistics for all covariances involving N_c:  < w'N_c' >.
    if ( stats_metadata%l_stats_samp ) then

       if ( stats_metadata%iwpNcp > 0 ) then

          ! Covariance of vertical velocity and N_c.
          call stat_update_var( stats_metadata%iwpNcp, &
                                wpNcp, stats_zm )

       endif

    endif ! stats_metadata%l_stats_samp

    if ( stats_metadata%l_stats_samp ) then

       ! Total time tendency
       do k = 1, gr%nzt, 1
          call stat_end_update_pt( stats_metadata%iNcm_bt, k, &
                                   Ncm(k) / dt, stats_zt )
       enddo ! k = 1, gr%nzt, 1

    endif ! stats_metadata%l_stats_samp

    return

  end subroutine advance_Ncm

  !=============================================================================
  subroutine microphys_solve( gr, solve_type, l_sed, &
                              lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, &
                              cloud_frac, &
                              tridiag_solve_method, &
                              stats_metadata, &
                              stats_zt, &
                              lhs, rhs, hmm, err_info )

    ! Description:
    ! Solve the tridiagonal system for hydrometeor variable.

    ! References:
    !  None
    !---------------------------------------------------------------------------

    use grid_class, only: grid

    use constants_clubb, only: &
        cloud_frac_min, & ! Constant(s)
        fstderr

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error                 ! Constant

    use matrix_solver_wrapper, only: &
        tridiag_solve ! Procedure(s)

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff  ! Use in cloud values of Nc for diffusion

    use stats_variables, only: &
        stats_metadata_type

    use stats_type_utilities, only: &
        stat_update_var_pt, & ! Procedure(s)
        stat_end_update_pt

    use stats_type, only: &
        stats ! Type

    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    !------------------------ Input Variables ------------------------
    type (grid), intent(in) :: gr

    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      cloud_frac    ! Cloud fraction (thermodynamic levels)        [-]

    integer, intent(in) :: &
      tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

    real( kind = core_rknd ), intent(in), dimension(3,1,gr%nzt) :: & 
      lhs_ta    ! LHS corresponding to contribution from turbulent adv.  [1/s]

    real( kind = core_rknd ), intent(in), dimension(3,gr%nzt) :: &
      sed_diff_lhs,       & ! Variables use to save stats results
      sed_turb_lhs
      
    real( kind = core_rknd ), intent(in), dimension(3,gr%nzt) :: & 
      lhs_ma    ! LHS corresponding to contribution from mean advection  [1/s]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------ Input/Output Variables ------------------------
    type(stats), intent(inout) :: &
      stats_zt

    real( kind = core_rknd ), intent(inout), dimension(3,gr%nzt) :: &
      lhs    ! Left hand side

    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: & 
      rhs    ! Right hand side vector

    real( kind = core_rknd ), intent(inout), dimension(gr%nzt) :: & 
      hmm    ! Mean value of hydrometeor (thermodynamic levels)    [units vary]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !------------------------ Local Variables ------------------------
    integer :: k, km1, kp1  ! Array indices

    integer :: & 
      ihmm_ma,  & ! Mean advection budget stats toggle
      ihmm_ta,  & ! Turbulent advection budget stats toggle
      ihmm_sd,  & ! Mean sedimentation budget stats toggle
      ihmm_ts     ! Turbulent sedimentation budget stats toggle

    !---------------------------- Begin Code ----------------------------

    ! Initializing ihmm_ma, ihmm_ta, ihmm_sd, ihmm_ts, and in order to avoid
    ! compiler warnings.
    ihmm_ma = 0
    ihmm_ta = 0
    ihmm_sd = 0
    ihmm_ts = 0

    select case( solve_type )
    case( "rrm" )
      ihmm_ma = stats_metadata%irrm_ma
      ihmm_ta = stats_metadata%irrm_ta
      ihmm_sd = stats_metadata%irrm_sd
      ihmm_ts = stats_metadata%irrm_ts
    case( "rim" )
      ihmm_ma = stats_metadata%irim_ma
      ihmm_ta = stats_metadata%irim_ta
      ihmm_sd = stats_metadata%irim_sd
      ihmm_ts = 0
    case( "rsm" )
      ihmm_ma = stats_metadata%irsm_ma
      ihmm_ta = stats_metadata%irsm_ta
      ihmm_sd = stats_metadata%irsm_sd
      ihmm_ts = 0
    case( "rgm" )
      ihmm_ma = stats_metadata%irgm_ma
      ihmm_ta = stats_metadata%irgm_ta
      ihmm_sd = stats_metadata%irgm_sd
      ihmm_ts = 0
    case( "Ncm" )
      ihmm_ma = stats_metadata%iNcm_ma
      ihmm_ta = stats_metadata%iNcm_ta
      ihmm_sd = 0
      ihmm_ts = 0
    case( "Nrm" )
      ihmm_ma = stats_metadata%iNrm_ma
      ihmm_ta = stats_metadata%iNrm_ta
      ihmm_sd = stats_metadata%iNrm_sd
      ihmm_ts = stats_metadata%iNrm_ts
    case( "Nim" )
      ihmm_ma = stats_metadata%iNim_ma
      ihmm_ta = stats_metadata%iNim_ta
      ihmm_sd = stats_metadata%iNim_sd
      ihmm_ts = 0
    case( "Nsm" )
      ihmm_ma = stats_metadata%iNsm_ma
      ihmm_ta = stats_metadata%iNsm_ta
      ihmm_sd = stats_metadata%iNsm_sd
      ihmm_ts = 0
    case( "Ngm" )
      ihmm_ma = stats_metadata%iNgm_ma
      ihmm_ta = stats_metadata%iNgm_ta
      ihmm_sd = stats_metadata%iNgm_sd
      ihmm_ts = 0
    case default
      ihmm_ma = 0
      ihmm_ta = 0
      ihmm_sd = 0
      ihmm_ts = 0
    end select

    ! Solve system using a tridiag_solve.
    call tridiag_solve( solve_type, tridiag_solve_method, & ! Intent(in)
                        gr%nzt,                           & ! Intent(in)
                        lhs, rhs, err_info,               & ! Intent(inout)
                        hmm )                               ! Intent(out)

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "calling tridiag_solve in microphys_solve"
        return
      endif
    endif

    ! Statistics
    if ( stats_metadata%l_stats_samp ) then

       do k = 1, gr%nzt, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nzt )

          ! Finalize implicit contributions

          ! hmm term ma is completely implicit; call stat_update_var_pt.
          if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

             ! For Ncm, we divide by cloud_frac when entering the subroutine,
             ! but do not multiply until we return from the subroutine, so we
             ! must account for this here for the budget to balance.
             call stat_update_var_pt( ihmm_ma, k, &
                     -lhs_ma(3,k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) &
                     -lhs_ma(2,k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) &
                     -lhs_ma(1,k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                      stats_zt )

          else

             call stat_update_var_pt( ihmm_ma, k, &
                                      -lhs_ma(3,k) * hmm(km1) &
                                      -lhs_ma(2,k) * hmm(k) &
                                      -lhs_ma(1,k) * hmm(kp1), stats_zt)

          endif

          ! hmm term sd is completely implicit; call stat_update_var_pt.
          if ( l_sed ) then
             call stat_update_var_pt( ihmm_sd, k, & 
                                      - sed_diff_lhs(3,k) * hmm(km1) & 
                                      - sed_diff_lhs(2,k) * hmm(k) & 
                                      - sed_diff_lhs(1,k) * hmm(kp1), stats_zt )
          endif

          ! hmm term ts has both implicit and explicit components; call
          ! stat_end_update_pt.
          if ( l_sed ) then
             call stat_end_update_pt( ihmm_ts, k, & 
                                      - sed_turb_lhs(3,k) * hmm(km1) & 
                                      - sed_turb_lhs(2,k) * hmm(k) & 
                                      - sed_turb_lhs(1,k) * hmm(kp1), stats_zt )
          endif

          ! hmm term ta has both implicit and explicit components; call
          ! stat_end_update_pt.
          if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

             ! For Ncm, we divide by cloud_frac when entering the subroutine,
             ! but do not multiply until we return from the subroutine, so we
             ! must account for this here for the budget to balance.
             call stat_end_update_pt( ihmm_ta, k, & 
                   -lhs_ta(3,1,k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
                   -lhs_ta(2,1,k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
                   -lhs_ta(1,1,k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                      stats_zt )

          else

             call stat_end_update_pt( ihmm_ta, k, & 
                                      -lhs_ta(3,1,k) * hmm(km1) & 
                                      -lhs_ta(2,1,k) * hmm(k) & 
                                      -lhs_ta(1,1,k) * hmm(kp1), stats_zt )

          endif

       enddo ! 1..gr%nzt

    endif ! stats_metadata%l_stats_samp


    return

  end subroutine microphys_solve

  !=============================================================================
  subroutine microphys_lhs( gr, solve_type, l_sed, dt, K_hm, nu, wm_zt, &
                            V_hm, V_hmt, &
                            Vhmphmp_zt_impc, &
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                            l_upwind_xm_ma, &
                            stats_metadata, &
                            lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs, &
                            lhs )

    ! Description:
    ! Setup the matrix of implicit contributions to a term.
    ! Can include the effects of sedimentation, diffusion, and advection.
    ! The Morrison microphysics has an explicit sedimentation code, which is
    ! handled elsewhere.
    !
    ! Notes:
    ! Setup for tridiagonal system and boundary conditions should be the same as
    ! the original rain subroutine code.
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        zm2zt_api, & ! Procedure(s)
        zt2zm_api    ! Procedure(s)

    use grid_class, only: grid ! Type

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use mean_adv, only:  & 
        term_ma_zt_lhs ! Procedure(s)

    use constants_clubb, only: &
        one,         & ! Constant(s)
        one_half,    &
        zero

    use parameters_microphys, only: &
        l_upwind_diff_sed  ! Use "upwind" differencing approx. for sedimentation

    use stats_variables, only: &
        stats_metadata_type

    implicit none

    !------------------------------ Input Variables ------------------------------
    type (grid), intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real( kind = core_rknd ), intent(in) ::  & 
      dt       ! Model timestep                                          [s]

    real( kind = core_rknd ), intent(in) ::  & 
      nu       ! Background diffusion coefficient                        [m^2/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nzt) ::  & 
      wm_zt, & ! w wind component on thermodynamic levels                [m/s]
      V_hmt    ! Sedimentation velocity of hydrometeor (thermo. levels)  [m/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nzm) ::  & 
      V_hm,  & ! Sedimentation velocity of hydrometeor (momentum levels) [m/s]
      K_hm     ! Coefficient of diffusion (turb. adv.) for hydrometeor   [m^2/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in)  :: & 
      Vhmphmp_zt_impc, & ! Implicit comp. of <V_hm'h_m'> on t-levs  [units(m/s)]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in)  :: & 
      rho_ds_zm          ! Dry, static density on momentum levels   [kg/m^3]

    logical, intent(in) :: &
      l_upwind_xm_ma ! This flag determines whether we want to use an upwind differencing
                     ! approximation rather than a centered differencing for turbulent or
                     ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------------ Output Variables ------------------------------
    real( kind = core_rknd ), intent(out), dimension(3,1,gr%nzt) :: & 
      lhs_ta    ! LHS corresponding to contribution from turbulent adv.  [1/s]

    real( kind = core_rknd ), intent(out), dimension(3,gr%nzt) :: &
      sed_diff_lhs,       & ! Variables use to save stats results
      sed_turb_lhs
      
    real( kind = core_rknd ), intent(out), dimension(3,gr%nzt) :: & 
      lhs_ma    ! LHS corresponding to contribution from mean advection  [1/s]

    real( kind = core_rknd ), intent(out), dimension(3,gr%nzt) :: & 
      lhs    ! Left hand side of tridiagonal matrix

    !------------------------------ Local Variables ------------------------------
    real( kind = core_rknd ), dimension(gr%nzm) :: &
      Vhmphmp_impc ! Implicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]

    real( kind = core_rknd ), dimension(1,gr%nzt) :: &
      Kh_zt,                & ! Eddy diffusivity coefficient, thermo levels [m2/s]
      invrs_rho_ds_zt_col     ! Inv. dry, static density @ thermo. levs.    [m^3/kg]

    real( kind = core_rknd ), dimension(1,gr%nzm) :: &
      Kh_zm,                & ! Eddy diffusivity coefficient, momentum levels [m2/s]
      rho_ds_zm_col           ! Dry, static density on momentum levels        [kg/m^3]
      
    real( kind = core_rknd ), dimension(1) :: &
      nu_col
      
    real( kind = core_rknd ), dimension(1,gr%nzt) ::  & 
      wm_zt_col

    integer :: k, km1, kp1, i  ! Array indices

    integer :: & 
      ihmm_ma, & ! Mean advection budget stats toggle
      ihmm_ta, & ! Turbulent advection budget stats toggle
      ihmm_sd, & ! Mean sedimentation budget stats toggle
      ihmm_ts    ! Turbulent sedimentation budget stats toggle

    !---------------------------- Begin Code ----------------------------

    ! Initializing ihmm_ma, ihmm_sd, ihmm_ts, and ihmm_ta in order to avoid
    ! compiler warnings.
    ihmm_ma = 0
    ihmm_ta = 0
    ihmm_sd = 0
    ihmm_ts = 0

    select case( solve_type )
    case ( "rrm" )
      ihmm_ma = stats_metadata%irrm_ma
      ihmm_ta = stats_metadata%irrm_ta
      ihmm_sd = stats_metadata%irrm_sd
      ihmm_ts = stats_metadata%irrm_ts
    case ( "Nrm" )
      ihmm_ma = stats_metadata%iNrm_ma
      ihmm_ta = stats_metadata%iNrm_ta
      ihmm_sd = stats_metadata%iNrm_sd
      ihmm_ts = stats_metadata%iNrm_ts
    case ( "rim" )
      ihmm_ma = stats_metadata%irim_ma
      ihmm_ta = stats_metadata%irim_ta
      ihmm_sd = stats_metadata%irim_sd
      ihmm_ts = 0
    case ( "rsm" )
      ihmm_ma = stats_metadata%irsm_ma
      ihmm_ta = stats_metadata%irsm_ta
      ihmm_sd = stats_metadata%irsm_sd
      ihmm_ts = 0
    case ( "rgm" )
      ihmm_ma = stats_metadata%irgm_ma
      ihmm_ta = stats_metadata%irgm_ta
      ihmm_sd = stats_metadata%irgm_sd
      ihmm_ts = 0
    case ( "Ncm" )
      ihmm_ma = stats_metadata%iNcm_ma
      ihmm_ta = stats_metadata%iNcm_ta
      ihmm_sd = 0
      ihmm_ts = 0
    case( "Nim" )
      ihmm_ma = stats_metadata%iNim_ma
      ihmm_ta = stats_metadata%iNim_ta
      ihmm_sd = stats_metadata%iNim_sd
      ihmm_ts = 0
    case( "Nsm" )
      ihmm_ma = stats_metadata%iNsm_ma
      ihmm_ta = stats_metadata%iNsm_ta
      ihmm_sd = stats_metadata%iNsm_sd
      ihmm_ts = 0
    case( "Ngm" )
      ihmm_ma = stats_metadata%iNgm_ma
      ihmm_ta = stats_metadata%iNgm_ta
      ihmm_sd = stats_metadata%iNgm_sd
      ihmm_ts = 0
    case default
      ihmm_ma = 0
      ihmm_ta = 0
      ihmm_sd = 0
      ihmm_ts = 0
    end select
    
    ! Some procedures expect arrays that are dimension(ngrdcol,nz), so the aruments 
    ! to be passed to those procedures have a dummy dimension hardcoded to 1.
    ! For visual clarity and consistency with other routines, we use i to index
    ! those dummy dimensions, but we always want i=1.
    i = 1

    ! Interpolate the implicit component of < V_hm'h_m' >, a momentum-level
    ! variable that is calculated on thermodynamic levels, from thermodynamic
    ! levels to momentum levels.
    Vhmphmp_impc = zt2zm_api( gr, Vhmphmp_zt_impc )

    ! Reset LHS Matrix for current timestep.
    lhs = zero

    ! LHS turbulent advection term.
    ! - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
    ! Note:  a down gradient closure approximation is made for < w'h_m' >, so
    !        the turbulent advection term is solved as an eddy-diffusion
    !        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
    ! A Crank-Nicholson time-stepping scheme is used for this term.
    Kh_zm(i,:) = K_hm
    Kh_zt(i,:) = max( zm2zt_api(gr,K_hm), 0._core_rknd )
    invrs_rho_ds_zt_col(i,:) = invrs_rho_ds_zt
    rho_ds_zm_col(i,:) = rho_ds_zm
    nu_col(i) = nu

    !$acc data copyin( gr, gr%invrs_dzm, gr%invrs_dzt, Kh_zm, Kh_zt, nu_col, &
    !$acc              invrs_rho_ds_zt_col, rho_ds_zm_col ) &
    !$acc       copyout( lhs_ta )
    call diffusion_zt_lhs( gr%nzm, gr%nzt, i, gr, Kh_zm, Kh_zt, nu_col, &
                           invrs_rho_ds_zt_col, rho_ds_zm_col, &
                           lhs_ta )
    !$acc end data

    do k = 1, gr%nzt, 1
       lhs_ta(:,i,k) = one_half * lhs_ta(:,i,k)
    enddo ! k = 1, gr%nzt, 1

    ! The lower boundary condition needs to be applied here at level 1.

    ! Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
    lhs_ta(kp1_tdiag,i,1) &
    = - one_half * invrs_rho_ds_zt(1) &
        * ( gr%invrs_dzt(i,1) * ( Kh_zm(i,2) + nu ) * rho_ds_zm(2) * gr%invrs_dzm(i,2) )

    ! Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
    lhs_ta(k_tdiag,i,1) &
    = + one_half * invrs_rho_ds_zt(1) &
        * ( gr%invrs_dzt(i,1) * ( Kh_zm(i,2) + nu ) * rho_ds_zm(2) * gr%invrs_dzm(i,2) )

    ! Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
    lhs_ta(km1_tdiag,i,1) = zero
    
    wm_zt_col(i,:) = wm_zt

    ! LHS mean advection term.
    !$acc data copyin( gr, gr%weights_zt2zm, gr%invrs_dzt, gr%invrs_dzm, wm_zt_col ) &
    !$acc      copyout( lhs_ma )
    call term_ma_zt_lhs( gr%nzm, gr%nzt, 1, wm_zt_col, gr%weights_zt2zm, & ! intent(in)
                         gr%invrs_dzt, gr%invrs_dzm,  & ! intent(in)
                         l_upwind_xm_ma, gr%grid_dir, & ! intent(in)
                         lhs_ma(:,:) )
    !$acc end data

    ! Setup LHS Matrix
    do k = 1, gr%nzt, 1

       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nzt )

       ! LHS time tendency.
       lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( one / dt )

       lhs(:,k) = lhs(:,k) + lhs_ta(:,i,k)

       ! LHS mean advection term.
       lhs(kp1_tdiag:km1_tdiag,k) &
       = lhs(kp1_tdiag:km1_tdiag,k) + lhs_ma(kp1_tdiag:km1_tdiag,k)

       if ( l_sed ) then

          ! LHS mean sedimentation term.
          ! Note: the Morrison microphysics has its own sedimentation code,
          ! which is applied through the _mc terms to each hydrometeor species.
          ! Therefore, l_sed will always be false when the Morrison microphysics
          ! is enabled.  -dschanen 24 Jan 2011
          if ( .not. l_upwind_diff_sed ) then

             ! Sedimentation (both mean and turbulent) uses centered
             ! differencing.  This is the default method.
             lhs(kp1_tdiag:km1_tdiag,k) & 
             = lhs(kp1_tdiag:km1_tdiag,k) & 
               + sed_centered_diff_lhs( gr, V_hm(kp1), V_hm(k), rho_ds_zm(kp1), &
                                        rho_ds_zm(k), invrs_rho_ds_zt(k), &
                                        gr%invrs_dzt(i,k), k )

          else

             ! Sedimentation (both mean and turbulent) uses "upwind"
             ! differencing.
             lhs(kp1_tdiag:km1_tdiag,k) & 
             = lhs(kp1_tdiag:km1_tdiag,k) & 
               + sed_upwind_diff_lhs( gr, V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                      rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                      gr%invrs_dzm(i,kp1), k )

          endif

          ! LHS turbulent sedimentation term.
          lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) & 
             + term_turb_sed_lhs( gr, Vhmphmp_impc(kp1), Vhmphmp_impc(k), &
                                  Vhmphmp_zt_impc(kp1), Vhmphmp_zt_impc(k), &
                                  rho_ds_zm(kp1), rho_ds_zm(k), &
                                  rho_ds_zt(kp1), rho_ds_zt(k), &
                                  gr%invrs_dzt(i,k), gr%invrs_dzm(i,kp1), &
                                  invrs_rho_ds_zt(k), k )

       endif ! l_sed

       if ( stats_metadata%l_stats_samp ) then

          ! Statistics:  implicit contributions to hydrometeor hmm.
          if ( ihmm_sd > 0 .and. l_sed ) then

             if ( .not. l_upwind_diff_sed ) then
                sed_diff_lhs(1:3,k) &
                = sed_centered_diff_lhs( gr, V_hm(kp1), V_hm(k), rho_ds_zm(kp1), &
                                         rho_ds_zm(k), invrs_rho_ds_zt(k), &
                                         gr%invrs_dzt(i,k), k )
             else
                sed_diff_lhs(1:3,k) &
                = sed_upwind_diff_lhs( gr, V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                       rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                       gr%invrs_dzm(i,kp1), k )
             endif

          endif

          if ( ihmm_ts > 0 .and. l_sed ) then
             sed_turb_lhs(1:3,k) &
             = term_turb_sed_lhs( gr, Vhmphmp_impc(kp1), Vhmphmp_impc(k), &
                                  Vhmphmp_zt_impc(kp1), Vhmphmp_zt_impc(k), &
                                  rho_ds_zm(kp1), rho_ds_zm(k), &
                                  rho_ds_zt(kp1), rho_ds_zt(k), &
                                  gr%invrs_dzt(i,k), gr%invrs_dzm(i,kp1), &
                                  invrs_rho_ds_zt(k), k )
          endif

       endif ! stats_metadata%l_stats_samp

    enddo ! 1..gr%nzt


    return

  end subroutine microphys_lhs

  !=============================================================================
  subroutine microphys_rhs( gr, solve_type, dt, l_sed, &
                            hmm, hmm_tndcy, &
                            K_hm, nu, cloud_frac, &
                            Vhmphmp_zt_expc, &
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                            stats_metadata, &
                            stats_zt, &
                            rhs )

    ! Description:
    ! Compute RHS vector for a given hydrometeor.
    ! This subroutine computes the explicit portion of the predictive equation
    ! for a given hydrometeor.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        zt2zm_api, &    ! Procedure(s)
        zm2zt_api

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        one_half,       & ! Constant(s)
        zero,           &
        cloud_frac_min

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff  ! Use in cloud values of Nc for diffusion

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use stats_variables, only: &
        stats_metadata_type

    use stats_type_utilities, only: &
        stat_begin_update_pt ! Procedure(s)

    use stats_type, only: &
        stats ! Type

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    !-------------------------- Input Variables --------------------------
    type (grid), intent(in) :: &
      gr

    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Duration of model timestep     [s]

    logical, intent(in) :: &
      l_sed    ! Flag for hydrometeor sedimentation

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      hmm,             & ! Mean value of hydrometeor (t-levs.)      [units]
      hmm_tndcy,       & ! Microphysics tendency (thermo. levels)   [units/s]
      cloud_frac,      & ! Cloud fraction                           [-]
      Vhmphmp_zt_expc, & ! Explicit comp. of <V_hm'h_m'> on t-levs  [units(m/s)]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: &
      K_hm,            & ! Coef. of diffusion for hydrometeor       [m^2/s]
      rho_ds_zm          ! Dry, static density on momentum levels   [kg/m^3]

    real( kind = core_rknd ), intent(in) :: &
      nu                 ! Background diffusion coefficient         [m^2/s]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !-------------------------- Input Variables --------------------------
    type(stats), intent(inout) :: &
      stats_zt

    !-------------------------- Output Variable --------------------------
    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: & 
      rhs   ! Right hand side

    !-------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(gr%nzm) :: &
      Vhmphmp_expc  ! Explicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]

    real( kind = core_rknd ), dimension(1,gr%nzt) :: &
      Kh_zt,                & ! Eddy diffusivity coefficient, thermo levels [m2/s]
      invrs_rho_ds_zt_col     ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    real( kind = core_rknd ), dimension(1,gr%nzm) :: &
      Kh_zm,                & ! Eddy diffusivity coefficient, momentum levels [m2/s]
      rho_ds_zm_col           ! Dry, static density on momentum levels      [kg/m^3]
      
    real( kind = core_rknd ), dimension(1) :: &
      nu_col

    real( kind = core_rknd ), dimension(3,gr%nzt) :: & 
      lhs_ta    ! LHS corresponding to contribution from turbulent advection


    integer :: i, k, kp1, km1  ! Array indices

    integer :: ihmm_ta, & ! Turbulent advection budget toggle.
               ihmm_ts    ! Turbulent sedimentation budget toggle.

    !-------------------------- Begin Code --------------------------

    ! Initializing ihmm_ta and ihmm_ts in order to avoid compiler warnings.
    ihmm_ta = 0
    ihmm_ts = 0

    select case( solve_type )
    case ( "rrm" )
      ihmm_ta = stats_metadata%irrm_ta
      ihmm_ts = stats_metadata%irrm_ts
    case ( "Nrm" )
      ihmm_ta = stats_metadata%iNrm_ta
      ihmm_ts = stats_metadata%iNrm_ts
    case( "rim" )
      ihmm_ta = stats_metadata%irim_ta
      ihmm_ts = 0
    case( "rsm" )
      ihmm_ta = stats_metadata%irsm_ta
      ihmm_ts = 0
    case( "rgm" )
      ihmm_ta = stats_metadata%irgm_ta
      ihmm_ts = 0
    case( "Ncm" )
      ihmm_ta = stats_metadata%iNcm_ta
      ihmm_ts = 0
    case( "Nim" )
      ihmm_ta = stats_metadata%iNim_ta
      ihmm_ts = 0
    case( "Nsm" )
      ihmm_ta = stats_metadata%iNsm_ta
      ihmm_ts = 0
    case( "Ngm" )
      ihmm_ta = stats_metadata%iNgm_ta
      ihmm_ts = 0
    case default
      ihmm_ta = 0
      ihmm_ts = 0
    end select
    
    ! Some procedures expect arrays that are dimension(ngrdcol,nz), so the aruments 
    ! to be passed to those procedures have a dummy dimension hardcoded to 1.
    ! For visual clarity and consistency with other routines, we use i to index
    ! those dummy dimensions, but we always want i=1.
    i = 1

    ! Interpolate the explicit component of < V_hm'h_m' >, a momentum-level
    ! variable that is calculated on thermodynamic levels, from thermodynamic
    ! levels to momentum levels.
    Vhmphmp_expc = zt2zm_api( gr, Vhmphmp_zt_expc )

    ! Initialize right-hand side vector to 0.
    rhs = zero

    ! LHS turbulent advection term.
    ! - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
    ! Note:  a down gradient closure approximation is made for < w'h_m' >, so
    !        the turbulent advection term is solved as an eddy-diffusion
    !        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
    ! A Crank-Nicholson time-stepping scheme is used for this term.
    Kh_zm(i,:) = K_hm
    Kh_zt(i,:) = max( zm2zt_api(gr,K_hm), zero ) 
    invrs_rho_ds_zt_col(i,:) = invrs_rho_ds_zt
    rho_ds_zm_col(i,:) = rho_ds_zm
    nu_col(i) = nu

    !$acc data copyin( gr, gr%invrs_dzm, gr%invrs_dzt, Kh_zm, Kh_zt, nu_col, &
    !$acc              invrs_rho_ds_zt_col, rho_ds_zm_col ) &
    !$acc       copyout( lhs_ta )
    call diffusion_zt_lhs( gr%nzm, gr%nzt, i, gr, Kh_zm, Kh_zt, nu_col, &
                           invrs_rho_ds_zt_col, rho_ds_zm_col, &
                           lhs_ta )
    !$acc end data

    do k = 1, gr%nzt, 1
       lhs_ta(:,k) = one_half * lhs_ta(:,k)
    enddo ! k = 1, gr%nz, 1

    ! The lower boundary condition needs to be applied here at level 1.

    ! Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
    lhs_ta(kp1_tdiag,1) &
    = - one_half * invrs_rho_ds_zt(1) &
        * ( gr%invrs_dzt(i,1) * ( Kh_zm(i,2) + nu ) * rho_ds_zm(2) * gr%invrs_dzm(i,2) )

    ! Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
    lhs_ta(k_tdiag,1) &
    = + one_half * invrs_rho_ds_zt(1) &
        * ( gr%invrs_dzt(i,1) * ( Kh_zm(i,2) + nu ) * rho_ds_zm(2) * gr%invrs_dzm(i,2) )

    ! Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
    lhs_ta(km1_tdiag,1) = zero

    ! Hydrometeor right-hand side (explicit portion of the code).
    do k = 1, gr%nzt, 1

       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nzt )

       ! RHS time tendency.
       rhs(k) = hmm(k) / dt

       ! RHS microphysics tendency term (autoconversion, accretion, evaporation,
       ! etc.).
       rhs(k) = rhs(k) + hmm_tndcy(k)

       rhs(k) &
       = rhs(k) &
         - lhs_ta(km1_tdiag,k) * hmm(km1) &
         - lhs_ta(k_tdiag,k) * hmm(k) &
         - lhs_ta(kp1_tdiag,k) * hmm(kp1)

       ! RHS turbulent sedimentation term.
       if ( l_sed ) then
          rhs(k) &
          = rhs(k) &
            + term_turb_sed_rhs( gr, Vhmphmp_expc(kp1), Vhmphmp_expc(k), &
                                 Vhmphmp_zt_expc(kp1), Vhmphmp_zt_expc(k), &
                                 rho_ds_zm(kp1), rho_ds_zm(k), &
                                 rho_ds_zt(kp1), rho_ds_zt(k), &
                                 gr%invrs_dzt(i,k), gr%invrs_dzm(i,kp1), &
                                 invrs_rho_ds_zt(k), k )
       endif

       
       if ( stats_metadata%l_stats_samp ) then

          ! Statistics: explicit contributions for the hydrometeor.

          ! hmm term ta has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on the right-hand side
          ! turbulent advection component.
          if ( ihmm_ta > 0 ) then

             if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

                ! For Ncm, we divide by cloud_frac when entering the subroutine,
                ! but do not multiply until we return from the subroutine, so we
                ! must account for this here for the budget to balance.
                call stat_begin_update_pt( ihmm_ta, k, & 
              lhs_ta(3,k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
              + lhs_ta(2,k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
              + lhs_ta(1,k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                           stats_zt )

             else

                call stat_begin_update_pt( ihmm_ta, k, & 
                                           lhs_ta(3,k) * hmm(km1) &
                                           + lhs_ta(2,k) * hmm(k)   &
                                           + lhs_ta(1,k) * hmm(kp1), stats_zt )

             endif

          endif

          ! hmm term ts has both implicit and explicit components; call
          ! stat_update_var_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_turb_sed_rhs.
          if ( ihmm_ts > 0 .and. l_sed ) then
             call stat_begin_update_pt( ihmm_ts, k, &
                 -term_turb_sed_rhs( gr, Vhmphmp_expc(kp1), Vhmphmp_expc(k), &
                                     Vhmphmp_zt_expc(kp1), Vhmphmp_zt_expc(k), &
                                     rho_ds_zm(kp1), rho_ds_zm(k), &
                                     rho_ds_zt(kp1), rho_ds_zt(k), &
                                     gr%invrs_dzt(i,k), gr%invrs_dzm(i,kp1), &
                                     invrs_rho_ds_zt(k), k ), &
                                        stats_zt )
          endif ! ihmm_ts > 0 and l_sed

       endif ! stats_metadata%l_stats_samp

    enddo ! k = 1, gr%nzt, 1
    

    return

  end subroutine microphys_rhs

  !=============================================================================
  pure function sed_centered_diff_lhs( gr, V_hmp1, V_hm, rho_ds_zmp1, &
                                       rho_ds_zm, invrs_rho_ds_zt, &
                                       invrs_dzt, level ) &
    result( lhs )

    ! Description:
    ! Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    ! the centered difference approximation to the vertical derivative.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The mean sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<hm>/dt equation.
    !
    ! This term is discretized as follows when using the centered-difference
    ! approximation:
    !
    ! The values of <hm> are found on the thermodynamic levels, while the values
    ! of <V_hm> are found on the momentum levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The variable <hm>
    ! is interpolated to the intermediate momentum levels.  At the intermediate
    ! momentum levels, the interpolated values of <hm> are multiplied by the
    ! values of <V_hm> and the values of rho_ds_zm.  Then, the derivative of
    ! (rho_ds*<V_hm>*<hm>) is taken over the central thermodynamic level, where
    ! it is multiplied by invrs_rho_ds_zt.
    !
    ! -----hmp1------------------------------------------------ t(k+1)
    !
    ! =============hm(interp)=====V_hm=====rho_ds_zm=========== m(k+1)
    !
    ! -----hm--------invrs_rho_ds_zt----d(rho_ds*V_hm*hm)/dz--- t(k)
    !
    ! =============hm(interp)=====V_hmm1===rho_ds_zmm1========= m(k)
    !
    ! -----hmm1------------------------------------------------ t(k-1)
    !
    ! The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1),
    ! respectively.  The letter "t" is used for thermodynamic levels and the
    ! letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
    !
    !
    ! Conservation Properties:
    !
    ! When a hydrometeor is sedimented to the ground (or out the lower boundary
    ! of the model), it is removed from the atmosphere (or from the model
    ! domain).  Thus, the quantity of the hydrometeor over the entire vertical
    ! domain should not be conserved due to the process of sedimentation.  Thus,
    ! not all of the column totals in the left-hand side matrix should be equal
    ! to 0. Instead, the sum of all the column totals should equal the flux of
    ! <hm> out the bottom (zm(1) level) of the domain,
    ! -rho_ds_zm(1) * V_hm(1) * hm(1), where the value of the hydrometeor at
    ! the surface (zm level 1) is set equal to the value of the hydrometeor at
    ! zt level 1, which is hm(1). Furthermore, most of the individual column
    ! totals should sum to 0, but the 1st and 2nd (from the left) columns should
    ! combine to sum to the flux out the bottom of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    ! - rho_ds_zm(1) * V_hm(1) * hm(1)
    ! = Sum_j Sum_i
    !   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzt )_i
    !   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    !
    ! The left-hand side matrix,
    ! ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_rho_ds_zt and invrs_dzt everywhere from the matrix below.  The sum
    ! over j leaves the column totals and the flux at zm(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first four vertical levels:
    !
    !     -------------------------------------------------------------------->
    !k=1 |   +invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)            0
    !    |    *invrs_dzt(k)        *invrs_dzt(k)
    !    |    *[ rho_ds_zm(k+1)    *rho_ds_zm(k+1)
    !    |       *V_hm(k+1)*B(k)   *V_hm(k+1)*A(k)
    !    |      -rho_ds_zm(k)
    !    |       *V_hm(k) ]
    !    |
    !k=2 |   -invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |    *invrs_dzt(k)        *invrs_dzt(k)          *invrs_dzt(k)
    !    |    *rho_ds_zm(k)        *[ rho_ds_zm(k+1)      *rho_ds_zm(k+1)
    !    |    *V_hm(k)*D(k)           *V_hm(k+1)*B(k)     *V_hm(k+1)*A(k)
    !    |                           -rho_ds_zm(k)
    !    |                            *V_hm(k)*C(k) ]
    !    |
    !k=3 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |                         *invrs_dzt(k)          *invrs_dzt(k)
    !    |                         *rho_ds_zm(k)          *[ rho_ds_zm(k+1)
    !    |                         *V_hm(k)*D(k)             *V_hm(k+1)*B(k)
    !    |                                                  -rho_ds_zm(k)
    !    |                                                   *V_hm(k)*C(k) ]
    !    |
    !k=4 |           0                     0             -invrs_rho_ds_zt(k)
    !    |                                                *invrs_dzt(k)
    !    |                                                *rho_ds_zm(k)
    !    |                                                *V_hm(k)*D(k)
    !    |
    !   \ /
    !
    ! The variables A(k), B(k), C(k), and D(k) are weights of interpolation
    ! around the central thermodynamic level (k), such that:
    !
    ! A(k) = ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ),
    ! B(k) = 1 - [ ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
    !      = 1 - A(k);
    ! C(k) = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
    ! D(k) = 1 - [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    !      = 1 - C(k).
    !
    ! Furthermore, for all intermediate thermodynamic grid levels (as long as
    ! k /= gr%nz and k /= 1), the four weighting factors have the following
    ! relationships:  A(k) = C(k+1) and B(k) = D(k+1).
    !
    ! Note:  The superdiagonal term from level 3 and both the main diagonal
    !        and superdiagonal terms from level 4 are not shown on this
    !        diagram.

    ! References:
    ! None

    ! Notes:
    !   Both COAMPS Microphysics and Brian Griffin's implementation use
    !   Khairoutdinov and Kogan (2000) for the calculation of rain
    !   mixing ratio and rain droplet number concentration sedimentation
    !   velocities, but COAMPS has only the local parameterization.
    !-----------------------------------------------------------------------

    use grid_class, only: grid

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      V_hmp1,          & ! Sedimentation velocity of hydrometeor (k+1)  [m/s]
      V_hm,            & ! Sedimentation velocity of hydrometeor (k)    [m/s]
      rho_ds_zmp1,     & ! Dry, static density at momentum level (k+1)  [kg/m^3]
      rho_ds_zm,       & ! Dry, static density at momentum level (k)    [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. level (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                  [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mkp1, & ! Momentum level directly above central thermodynamic level.
      mk      ! Momentum level directly below central thermodynamic level.
      
    integer :: i

    ! ---- Begin Code ----
    
    ! Some procedures expect arrays that are dimension(ngrdcol,nzt), so the aruments 
    ! to be passed to those procedures have a dummy dimension hardcoded to 1.
    ! For visual clarity and consistency with other routines, we use i to index
    ! those dummy dimensions, but we always want i=1.
    i = 1

    ! Momentum level (k+1) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mkp1 = level + 1
    ! Momentum level (k) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mk   = level

    if ( level == 1 ) then

       ! k = 1 (bottom level); lower boundary level.
       ! Special discretization where the value of the hydrometeor at the lower
       ! boundary or surface (momentum level 1) is implicitly set equal to the
       ! value of the hydrometeor at thermodynamic level 1.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag)  & 
       = + invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zmp1 * V_hmp1 * gr%weights_zt2zm(i,mkp1,t_above)

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)  & 
       = + invrs_rho_ds_zt &
           * invrs_dzt &
           * ( rho_ds_zmp1 * V_hmp1 * gr%weights_zt2zm(i,mkp1,t_below) &
               - rho_ds_zm * V_hm )

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    elseif ( level > 1 .and. level < gr%nzt ) then

       ! Most of the interior model; normal conditions.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag)  & 
       = + invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zmp1 * V_hmp1 * gr%weights_zt2zm(i,mkp1,t_above)

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)  & 
       = + invrs_rho_ds_zt &
           * invrs_dzt &
           * ( rho_ds_zmp1 * V_hmp1 * gr%weights_zt2zm(i,mkp1,t_below) & 
               - rho_ds_zm * V_hm * gr%weights_zt2zm(i,mk,t_above)  )

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag)  & 
       = - invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zm * V_hm * gr%weights_zt2zm(i,mk,t_below)


    elseif ( level == gr%nzt ) then

       ! k = gr%nzt (top level); upper boundary level; no flux.
       ! Special discretization where the value of the hydrometeor at the upper 
       ! boundary (momentum level gr%nzm) is implicitly set equal to 0. The
       ! process of sedimentation can only decrease the value of the hydrometeor
       ! at level gr%nzt. 

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)  & 
       = - invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zm * V_hm * gr%weights_zt2zm(i,mk,t_above)

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag)  & 
       = - invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zm * V_hm * gr%weights_zt2zm(i,mk,t_below)


    endif


    return

  end function sed_centered_diff_lhs

  !=============================================================================
  pure function sed_upwind_diff_lhs( gr, V_hmt, V_hmtp1, rho_ds_zt, &
                                     rho_ds_ztp1, invrs_rho_ds_zt, &
                                     invrs_dzmp1, level ) &
    result( lhs )

    ! Description:
    ! Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    ! the "upwind" difference approximation to the vertical derivative.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The mean sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<hm>/dt equation.
    !
    ! This term is discretized as follows when using the upwind-difference
    ! approximation:
    !
    ! The values of <hm> and the values of V_hmt are found on the thermodynamic
    ! levels.  Additionally, the values of rho_ds_zt and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of <hm> are multiplied by the values of
    ! V_hmt and the values of rho_ds_zt.  Then, the derivative of
    ! (rho_ds*<V_hm>*<hm>) is taken between the thermodynamic level above the
    ! central thermodynamic level and the central thermodynamic level.  The
    ! derivative is multiplied by invrs_rho_ds_zt.
    !
    ! --hmp1--V_hmtp1--rho_ds_ztp1--------------------------------------- t(k+1)
    !
    ! =================================================================== m(k+1)
    !
    ! --hm----V_hmt----rho_ds_zt--invrs_rho_ds_zt--d(rho_ds*V_hm*hm)/dz-- t(k)
    !
    ! The vertical indices t(k+1), m(k+1), and t(k) correspond with altitudes
    ! zt(k+1), zm(k+1), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) )
    !
    !
    ! Conservation Properties:
    !
    ! When a hydrometeor is sedimented to the ground (or out the lower boundary
    ! of the model), it is removed from the atmosphere (or from the model
    ! domain).  Thus, the quantity of the hydrometeor over the entire vertical
    ! domain should not be conserved due to the process of sedimentation.  Thus,
    ! not all of the column totals in the left-hand side matrix should be equal
    ! to 0. Instead, the sum of all the column totals should equal the flux of
    ! <hm> out the bottom (zm(1) level, for which the value of the hydrometeor
    ! is set equal to the value of the hydrometeor at the zt(1) level) of the
    ! domain, -rho_ds_zt(1) * V_hmt(1) * hm(1).  Furthermore, most of the
    ! individual column totals should sum to 0, but the 2nd (from the left)
    ! column should be equal to the flux out the bottom of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    ! - rho_ds_zt(1) * V_hmt(1) * hm(1)
    ! = Sum_j Sum_i
    !   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzm )_i
    !   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    !
    ! The left-hand side matrix,
    ! ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_rho_ds_zt and invrs_dzm everywhere from the matrix below.  The sum
    ! over j leaves the column totals and the flux at zt(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first three vertical levels:
    !
    !     -------------------------------------------------------------------->
    !k=1 | -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)              0
    !    |  *invrs_dzm(k+1)        *invrs_dzm(k+1)
    !    |  *rho_ds_zt(k)          *rho_ds_zt(k+1)
    !    |  *V_hmt(k)              *V_hmt(k+1)
    !    |
    !k=2 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |                         *invrs_dzm(k+1)        *invrs_dzm(k+1)
    !    |                         *rho_ds_zt(k)          *rho_ds_zt(k+1)
    !    |                         *V_hmt(k)              *V_hmt(k+1)
    !    |
    !k=3 |           0                     0             -invrs_rho_ds_zt(k)
    !    |                                                *invrs_dzm(k+1)
    !    |                                                *rho_ds_zt(k)
    !    |                                                *V_hmt(k)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 3 is not shown on this diagram.

    ! References:
    ! None

    ! Notes:
    ! Both COAMPS Microphysics and Brian Griffin's implementation use
    ! Khairoutdinov and Kogan (2000) for the calculation of rain
    ! mixing ratio and rain droplet number concentration sedimentation
    ! velocities, but COAMPS has only the local parameterization.
    !
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only: grid

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      V_hmtp1,         & ! Sed. velocity of hydrometeor at t-lev (k+1)  [m/s]
      V_hmt,           & ! Sed. velocity of hydrometeor at t-lev (k)    [m/s]
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1)   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. level (k) [m^3/kg]
      invrs_dzmp1        ! Inverse of grid spacing over m-lev. (k+1)    [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! ---- Begin Code ----

    ! Sedimention is always a downward process, so we omit the upward case
    ! (i.e. the V_hmt variable will always be negative).
    if ( level == gr%nzt ) then

       ! k = gr%nzt (top level); upper boundary level; no flux.
       ! Special discretization where the value of the hydrometeor above the
       ! upper boundary is implicitly set equal to 0. The process of
       ! sedimentation can only decrease the value of the hydrometeor
       ! at level gr%nzt. 

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = - invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_zt * V_hmt

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    else

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = + invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_ztp1 * V_hmtp1

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = - invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_zt * V_hmt

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero

    endif

    return

  end function sed_upwind_diff_lhs

  !=============================================================================
  pure function term_turb_sed_lhs( gr, Vhmphmp_impcp1, Vhmphmp_impc, &
                                   Vhmphmp_zt_impcp1, Vhmphmp_zt_impc, &
                                   rho_ds_zmp1, rho_ds_zm, &
                                   rho_ds_ztp1, rho_ds_zt, &
                                   invrs_dzt, invrs_dzmp1, &
                                   invrs_rho_ds_zt, level ) &
    result( lhs )

    ! Description:
    ! Turbulent sedimentation of a hydrometeor:  implicit portion of the code.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The turbulent sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    ! on < hm > in the manner:
    !
    ! < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    !
    ! This term can also be solved for completely explicitly (it's original
    ! form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    ! < V_hm'hm' >.  The equation becomes:
    !
    ! - (1/rho_ds) 
    !   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    !
    ! where the timestep index (t+1) means that the value of < hm > being used
    ! is from the next timestep, which is being advanced to in solving the
    ! d<hm>/dt equation.  Implicit and explicit portions of this term are
    ! produced.  The implicit portion of this term is:
    !
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < hm >(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the d[ ] / dz term is
    !        changed to a "+".
    !
    ! This term can be discretized using the centered-difference approximation
    ! (which is preferred), or else using the "upwind"-difference approximation.
    !
    ! The implicit portion of this term is discretized as follows when using
    ! the centered-difference approximation:
    ! 
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    ! the momentum levels, and the values of invrs_rho_ds_zt are found on the
    ! thermodynamic levels.  The variables < hm > and Vhmphmp_zt_impc are both
    ! interpolated to the intermediate momentum levels.  At the momentum levels,
    ! the values of interpolated < hm > and interpolated Vhmphmp_zt_impc are
    ! multiplied together, and their products are multiplied by the values of
    ! rho_ds_zm.  The mathematical expression F is the product of these three
    ! variables at momentum levels.  Then, the derivative dF/dz is taken over
    ! the central thermodynamic level, where it is multiplied by
    ! invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zm * Vhmphmp_impc(interp) * hmm(interp).
    !
    !
    ! ----hmmp1--------Vhmphmp_zt_impcp1--------------------------------- t(k+1)
    !
    ! =====hmm(interp)=====Vhmphmp_impcp1(interp)=====rho_ds_zmp1======== m(k+1)
    !
    ! ----hmm----------Vhmphmp_zt_impc-----invrs_rho_ds_zt-----dF/dz----- t(k)
    !
    ! =====hmm(interp)=====Vhmphmp_impc(interp)=======rho_ds_zm========== m(k)
    !
    ! ----hmmm1--------Vhmphmp_zt_impcm1--------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
    !
    ! The implicit portion of this term is discretized as follows when using
    ! the upwind-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    ! values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of < hm > and Vhmphmp_zt_impc are
    ! multiplied together, and their products are multiplied by the values of
    ! rho_ds_zt.  The mathematical expression F is the product of these three
    ! variables at thermodynamic levels.  Then, the derivative dF/dz is taken
    ! between the thermodynamic level above the central thermodynamic level and
    ! the central thermodynamic level.  The derivative is multiplied
    ! by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zt * Vhmphmp_zt_impc * hmm.
    !
    !
    ! --hmmp1---Vhmphmp_zt_impcp1---rho_ds_ztp1-------------------------- t(k+1)
    !
    ! =================================================================== m(k+1)
    !
    ! --hmm-----Vhmphmp_zt_impc-----rho_ds_zt---invrs_rho_ds_zt---dF/dz-- t(k)
    !
    ! The vertical indices t(k+1), m(k+1), and t(k) correspond with altitudes
    ! zt(k+1), zm(k+1), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) ).

    ! References:
    !  None
    !
    ! Notes:
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only: grid

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_microphys, only: &
        l_upwind_diff_sed  ! Use "upwind" differencing approx. for sedimentation

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    ! Constant parameters
    integer, parameter :: &
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: &
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Vhmphmp_impcp1,    & ! Imp. comp. <V_hm'h_m'> interp. m-lev (k+1) [vary]
      Vhmphmp_impc,      & ! Imp. comp. <V_hm'h_m'> interp. m-lev (k)   [vary]
      Vhmphmp_zt_impcp1, & ! Imp. comp. <V_hm'h_m'>|_zt; t-lev (k+1)    [vary]
      Vhmphmp_zt_impc,   & ! Imp. comp. <V_hm'h_m'>|_zt; t-lev (k)      [vary]
      rho_ds_zmp1,       & ! Dry, static density at moment. lev (k+1)   [kg/m^3]
      rho_ds_zm,         & ! Dry, static density at moment. lev (k)     [kg/m^3]
      rho_ds_ztp1,       & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,         & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_dzt,         & ! Inverse of grid spacing over t-levs. (k)   [1/m]
      invrs_dzmp1,       & ! Inverse of grid spacing over m-levs. (k+1) [1/m]
      invrs_rho_ds_zt      ! Inv dry, static density @ thermo lev (k)   [m^3/kg]

    integer, intent(in) :: &
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mkp1, & ! Momentum level directly above central thermodynamic level.
      mk      ! Momentum level directly below central thermodynamic level.

    integer :: i

    !---------------------------- Begin Code ----------------------------

    ! Some procedures expect arrays that are dimension(ngrdcol,nz), so the aruments
    ! to be passed to those procedures have a dummy dimension hardcoded to 1.
    ! For visual clarity and consistency with other routines, we use i to index
    ! those dummy dimensions, but we always want i=1.
    i = 1

    ! Momentum level (k+1) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mkp1 = level + 1
    ! Momentum level (k) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mk   = level


    ! LHS (implicit component) of turbulent sedimentation term,
    ! - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
    ! = - (1/rho_ds)
    !     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
    ! Implicit component:
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < h_m > ) / dz.
    if ( .not. l_upwind_diff_sed ) then

       ! Sedimentation (both mean and turbulent) uses centered differencing.

       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level.
          ! Special discretization where the value of the hydrometeor at the
          ! lower boundary or surface (momentum level 1) is implicitly set equal
          ! to the value of the hydrometeor at thermodynamic level 1.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zmp1 * Vhmphmp_impcp1 * gr%weights_zt2zm(i,mkp1,t_above)

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * ( rho_ds_zmp1 * Vhmphmp_impcp1 &
                              * gr%weights_zt2zm(i,mkp1,t_below) & 
                  - rho_ds_zm * Vhmphmp_impc )

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level > 1 .and. level < gr%nzt ) then

          ! Most of the interior model; normal conditions.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zmp1 * Vhmphmp_impcp1 * gr%weights_zt2zm(i,mkp1,t_above)

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * ( rho_ds_zmp1 * Vhmphmp_impcp1 &
                              * gr%weights_zt2zm(i,mkp1,t_below) & 
                  - rho_ds_zm * Vhmphmp_impc &
                              * gr%weights_zt2zm(i,mk,t_above)  )

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag)  & 
          = - invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zm * Vhmphmp_impc * gr%weights_zt2zm(i,mk,t_below)


       elseif ( level == gr%nzt ) then

          ! k = gr%nzt (top level); upper boundary level; no flux.
          ! Special discretization where the values of both the hydrometeor and
          ! V_hm'hm' at the upper boundary (momentum level gr%nzm) are set equal
          ! to 0.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag) &
          = - invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zm * Vhmphmp_impc * gr%weights_zt2zm(i,mk,t_above)

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) &
          = - invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zm * Vhmphmp_impc * gr%weights_zt2zm(i,mk,t_below)


       endif  ! level


    else ! l_upwind_diff_sed

       ! Sedimentation (both mean and turbulent) uses "upwind" differencing.
       if ( level < gr%nzt ) then

          ! Most of the interior model; normal conditions.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) &
          = + invrs_rho_ds_zt &
              * invrs_dzmp1 * rho_ds_ztp1 * Vhmphmp_zt_impcp1

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag) &
          = - invrs_rho_ds_zt &
              * invrs_dzmp1 * rho_ds_zt * Vhmphmp_zt_impc

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level == gr%nzt ) then

          ! k = gr%nzt (top level); upper boundary level; no flux.
          ! Special discretization where the values of both the hydrometeor and
          ! V_hm'hm' at the upper boundary (momentum level gr%nzm) are set equal
          ! to 0.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag) &
          = - invrs_rho_ds_zt &
              * invrs_dzmp1 * rho_ds_zt * Vhmphmp_zt_impc

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero

       endif  ! level

    endif

    return

  end function term_turb_sed_lhs

  !=============================================================================
  pure function term_turb_sed_rhs( gr, Vhmphmp_expcp1, Vhmphmp_expc, &
                                   Vhmphmp_zt_expcp1, Vhmphmp_zt_expc, &
                                   rho_ds_zmp1, rho_ds_zm, &
                                   rho_ds_ztp1, rho_ds_zt, &
                                   invrs_dzt, invrs_dzmp1, &
                                   invrs_rho_ds_zt, level ) &
    result( rhs )

    ! Description:
    ! Turbulent sedimentation of a hydrometeor:  explicit portion of the code.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The turbulent sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    ! on < hm > in the manner:
    !
    ! < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    !
    ! This term can also be solved for completely explicitly (it's original
    ! form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    ! < V_hm'hm' >.  The equation becomes:
    !
    ! - (1/rho_ds) 
    !   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    !
    ! where the timestep index (t+1) means that the value of < hm > being used
    ! is from the next timestep, which is being advanced to in solving the
    ! d<hm>/dt equation.  Implicit and explicit portions of this term are
    ! produced.  The explicit portion of this term is:
    !
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    !
    ! This term can be discretized using the centered-difference approximation
    ! (which is preferred), or else using the "upwind"-difference approximation.
    !
    ! The explicit portion of this term is discretized as follows when using
    ! the centered-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    ! the momentum levels, and the values of invrs_rho_ds_zt are found on the
    ! thermodynamic levels.  The variable Vhmphmp_zt_expc is interpolated to the
    ! intermediate momentum levels.  At the momentum levels, the values of
    ! interpolated Vhmphmp_zt_expc are multiplied by the values of rho_ds_zm.
    ! Then, the derivative d(rho_ds*Vhmphmp_zt_expc)/dz is taken over the
    ! central thermodynamic level, where it is multiplied by invrs_rho_ds_zt.
    !
    ! ---Vhmphmp_zt_expcp1----------------------------------------------- t(k+1)
    !
    ! ======Vhmphmp_expcp1(interp)=======rho_ds_zmp1===================== m(k+1)
    !
    ! ---Vhmphmp_zt_expc--invrs_rho_ds_zt--d(rho_ds*Vhmphmp_zt_expc)/dz-- t(k)
    !
    ! ======Vhmphmp_expc(interp)=========rho_ds_zm======================= m(k)
    !
    ! ---Vhmphmp_zt_expcm1----------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
    !
    ! The explicit portion of this term is discretized as follows when using
    ! the upwind-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    ! values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of Vhmphmp_zt_expc are multiplied by the
    ! values of rho_ds_zt.  The mathematical expression F is the product of
    ! these variables at thermodynamic levels.  Then, the derivative dF/dz is
    ! taken between the thermodynamic level above the central thermodynamic
    ! level and the central thermodynamic level.  The derivative is multiplied
    ! by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zt * Vhmphmp_zt_expc.
    !
    !
    ! -----Vhmphmp_zt_expcp1---rho_ds_ztp1------------------------------- t(k+1)
    !
    ! =================================================================== m(k+1)
    !
    ! -----Vhmphmp_zt_expc-----rho_ds_zt----invrs_rho_ds_zt----dF/dz----- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) ).
    
    ! References:
    !  None
    !
    ! Notes:
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only: grid

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_microphys, only: &
        l_upwind_diff_sed  ! Use "upwind" differencing approx. for sedimentation

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Vhmphmp_expcp1,    & ! Exp. comp. <V_hm'h_m'> interp. m-lev (k+1) [vary]
      Vhmphmp_expc,      & ! Exp. comp. <V_hm'h_m'> interp. m-lev (k)   [vary]
      Vhmphmp_zt_expcp1, & ! Exp. comp. <V_hm'h_m'>|_zt; t-lev (k+1)    [vary]
      Vhmphmp_zt_expc,   & ! Exp. comp. <V_hm'h_m'>|_zt; t-lev (k)      [vary]
      rho_ds_zmp1,       & ! Dry, static density at moment. lev (k+1)   [kg/m^3]
      rho_ds_zm,         & ! Dry, static density at moment. lev (k)     [kg/m^3]
      rho_ds_ztp1,       & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,         & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_dzt,         & ! Inverse of grid spacing over t-levs. (k)   [1/m]
      invrs_dzmp1,       & ! Inverse of grid spacing over m-levs. (k+1) [1/m]
      invrs_rho_ds_zt      ! Inv dry, static density @ thermo lev (k)   [m^3/kg]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    !---------------------------- Begin Code ----------------------------

    ! Initialize Right-hand side term.
    rhs = zero

    ! RHS (explicit component) of turbulent sedimentation term,
    ! - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
    ! = - (1/rho_ds)
    !     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
    ! Explicit component:  - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    if ( .not. l_upwind_diff_sed ) then

       ! Sedimentation (both mean and turbulent) uses centered differencing.
       if ( level < gr%nzt ) then

          ! Most of the interior model; normal conditions.
          rhs &
          = - invrs_rho_ds_zt &
              * invrs_dzt * ( rho_ds_zmp1 * Vhmphmp_expcp1 &
                              - rho_ds_zm * Vhmphmp_expc )


       elseif ( level == gr%nzt ) then

          ! k = gr%nzt (top level); upper boundary level; no flux.
          ! Special discretization where the value of V_hm'hm' at the upper
          ! boundary (momentum level gr%nzm) is set equal to 0.
          ! This results in an equation that is the same as setting
          ! Vhmphmp_expcp1 to 0 as found in the equation above for most of the
          ! interior model.
          rhs = + invrs_rho_ds_zt * invrs_dzt * rho_ds_zm * Vhmphmp_expc


       endif


    else ! l_upwind_diff_sed

       ! Sedimentation (both mean and turbulent) uses "upwind" differencing.
       if ( level < gr%nzt ) then

          ! Most of the interior model; normal conditions.
          rhs &
          = - invrs_rho_ds_zt &
              * invrs_dzmp1 * ( rho_ds_ztp1 * Vhmphmp_zt_expcp1 &
                                - rho_ds_zt * Vhmphmp_zt_expc )


       elseif ( level == gr%nzt ) then

          ! k = gr%nzt (top level); upper boundary level; no flux.
          ! Special discretization where the value of V_hm'hm' at the upper
          ! boundary (momentum level gr%nzm) is set equal to 0.
          ! This results in an equation that is the same as setting
          ! Vhmphmp_zt_expcp1 to 0 as found in the equation above for most of the
          ! interior model.
          rhs = + invrs_rho_ds_zt * invrs_dzmp1 * rho_ds_zt * Vhmphmp_zt_expc


       endif

    endif

    return

  end function term_turb_sed_rhs

  !=============================================================================
  function calculate_K_hm( gr, wp2, Kh_zm, Skw_zm, Lscale, &
                           hydromet_dim, hydromet_tol, &
                           hydromet, hydrometp2, &
                           clubb_params, &
                           l_use_non_local_diff_fac ) &
  result( K_hm )

    ! Description:
    ! The predictive equation for a hydrometeor, hm, contains a turbulent
    ! advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'hm'> )/dz.
    !
    ! The value of <w'hm'> can be calculated in many ways.  For simplicity, a
    ! down-gradient approximation is used here, where:
    !
    ! <w'hm'> = -K_hm * d<hm>/dz.
    !
    ! The coefficient of diffusion, K_hm, is variable and depends on multiple
    ! factors.  It optionally includes a non-local turbulent advection factor.

    ! References:
    ! CLUBB ticket 651 and CLUBB ticket 739.
    !-----------------------------------------------------------------------

    use grid_class, only: grid ! Type

    use grid_class, only: & 
        zt2zm_api    ! Procedure(s)

    use parameter_indices, only: &
        nparams,        & ! Variable(s)
        ic_K_hm,        & 
        ic_K_hmb,       &  
        iK_hm_min_coef

    use constants_clubb, only: & 
        one,  & ! Constant(s)
        zero, &
        eps

    use clubb_precision, only:  & 
        core_rknd    ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr

    integer, intent(in) :: &
      hydromet_dim

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: & 
      hydromet_tol

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: & 
      wp2,    & ! Variance of vertical velocity (momentum levels) [m^2/s^2]
      Kh_zm,  & ! Kh Eddy diffusivity on momentum grid            [m^2/s]
      Skw_zm    ! Skewness of w on momentum levels                [-]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      Lscale    ! Length-scale                                    [m]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(inout) :: &
      hydromet      ! Hydrometeor mean, <h_m> (thermo. levels)    [units]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of hydrometeor (overall) (m-levs.) [units^2]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_use_non_local_diff_fac    ! Flag to use a non-local factor for
                                  ! eddy-diffusivity applied to hydrometeors.

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim) :: &
      K_hm    ! Hydrometeor eddy diffusivity on momentum grid     [m^2/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim) :: &
      K_gamma ! Non-local factor of diffusion (t. adv.) for hydrometeors [m^2/s]

    integer :: k, km1, i, h    ! Loop indices

    !---------------------------- Begin Code ----------------------------

    ! Some procedures expect arrays that are dimension(ngrdcol,nz), so the aruments
    ! to be passed to those procedures have a dummy dimension hardcoded to 1.
    ! For visual clarity and consistency with other routines, we use i to index
    ! those dummy dimensions, but we always want i=1.
    i = 1

    ! Loop over all hydrometeors.
    do h = 1, hydromet_dim, 1

       ! Loop over all vertical levels for each hydrometeor.
       do k = 2, gr%nzm-1, 1

          km1 = max( k-1, 1 )

          K_hm(k,h) &
          = clubb_params(ic_K_hm) * Kh_zm(k) &
            * ( sqrt( hydrometp2(k,h) ) &
                / max( zt2zm_api( gr, hydromet(:,h), k ), hydromet_tol(h) ) ) &
            * ( one + abs( Skw_zm(k) ) ) 

          if ( l_use_non_local_diff_fac ) then
             K_gamma(k,h) &
             = one &
               - clubb_params(ic_K_hmb) &
                 * ( ( max( zt2zm_api( gr, Lscale(:), k ), 0.0_core_rknd ) &
                       / max( zt2zm_api( gr, hydromet(:,h), k ), hydromet_tol(h) ) ) &
                     * ( gr%invrs_dzm(i,k) &
                         * ( hydromet(k,h) - hydromet(km1,h) ) ) )

             K_hm(k,h) &
             = K_hm(k,h) * max( K_gamma(k,h), clubb_params(iK_hm_min_coef) )
          endif

          if ( abs( gr%invrs_dzm(i,k) &
                    * ( hydromet(k,h) - hydromet(km1,h) ) ) > eps ) then

             ! Ensure the abs( correlation ) between w and hydromet does not
             ! have a value greater than one.
             K_hm(k,h) &
             = min( K_hm(k,h), &
                    ( sqrt( wp2(k) ) * sqrt( hydrometp2(k,h) ) ) &
                    / abs( gr%invrs_dzm(i,k) &
                           * ( hydromet(k,h) - hydromet(km1,h) ) ) )

          endif ! | d<hm>/dz | > 0

       enddo ! k = 1, gr%nzm-1, 1

       ! Set K_hm at the lower and upper boundaries (not used in calculations).
       K_hm(1,h) = zero
       K_hm(gr%nzm,h) = zero

    enddo ! i = 1, hydromet_dim, 1


    return

  end function calculate_K_hm

  !=============================================================================
  function get_cloud_top_level( nzt, rcm, hydromet, &
                                hydromet_dim, iiri ) &
  result( cloud_top_level )

    ! Description:
    ! Find cloud top at a given model time step.  This function finds cloud top
    ! by looping downward from the top of the model and returning the index of
    ! the first vertical level that has a mean cloud water mixing ratio (or a
    ! mean cloud ice mixing ratio, when ice is included in the microphysics
    ! scheme) that is greater than the tolerance amount.  In a scenario that
    ! there is not any cloud found, the function returns a value of 1 (for
    ! vertical level 1, which is below the model surface).

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        rc_tol, & ! Constant(s)
        ri_tol, &
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzt, &    ! Number of model thermodynamic vertical grid levels
      hydromet_dim, &
      iiri

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      rcm    ! Mean cloud water mixing ratio                [kg/kg]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet    ! Hydrometeor mean, <h_m> (thermo. levels)    [units vary]

    ! Return Variable
    integer :: &
      cloud_top_level    ! Vertical level index of cloud top

    ! Local Variable
    real( kind = core_rknd ), dimension(nzt) :: &
      rim    ! Mean cloud ice mixing ratio                [kg/kg]

    integer :: k    ! Vertical level index

    !---------------------------- Begin Code ----------------------------

    ! Include cloud ice when the microphysics scheme includes ice mixing ratio.
    if ( iiri > 0 ) then
       rim = hydromet(:, iiri)
    else
       rim = zero
    endif ! iiri > 0

    ! Start at the model upper boundary and loop downwards until cloud top is
    ! found or the model lower boundary is reached.
    k = nzt
    do
       if ( rcm(k) > rc_tol .or. rim(k) > ri_tol ) then
          ! A level with mean cloud water mixing ratio (or mean cloud ice mixing
          ! ratio, when it is included in the microphysics scheme) that is
          ! greater than the tolerance amount has been found.  Cloud top has
          ! been found.
          cloud_top_level = k
          exit
       elseif ( k == 1 ) then
          ! There was not any cloud found in the model vertical domain.
          ! Return a value of 1.
          cloud_top_level = 1
          exit
       else
          ! Continue down another vertical level.
          k = k - 1
       endif
    enddo


    return

  end function get_cloud_top_level

  !=============================================================================
  subroutine write_adv_micro_errors( gr, dt, time_current, hydromet_dim,  & ! In
                                     wm_zt, wp2, &                          ! In
                                     exner, rho, rho_zm, rcm, &             ! In
                                     cloud_frac, Kh_zm, Skw_zm, &           ! In
                                     rho_ds_zm, rho_ds_zt, &                ! In
                                     invrs_rho_ds_zt, &                     ! In
                                     hydromet_mc, Ncm_mc, Lscale, &         ! In
                                     hydromet_vel_covar_zt_impc, &          ! In
                                     hydromet_vel_covar_zt_expc, &          ! In
                                     clubb_params, nu_vert_res_dep, &       ! In
                                     l_upwind_xm_ma, &                      ! In
                                     hydromet, hydromet_vel_zt, &           ! In
                                     hydrometp2, K_hm, Ncm, &               ! In
                                     Nc_in_cloud, rvm_mc, thlm_mc, &        ! In
                                     wphydrometp, wpNcp, err_info )         ! In

    ! Description:
    ! Writes to screen the values of all variables that are passed into and out
    ! of subroutine advance_microphys if a fatal error has been detected.

    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid    ! Type

    use constants_clubb, only: &
        fstderr    ! Variable(s)

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use error_code, only: &
        clubb_at_least_debug_level_api, & ! Procedure
        clubb_fatal_error             ! Constant

    use clubb_precision, only:  &
        time_precision, & ! Variable(s)
        core_rknd

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    ! Input Variables
    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: &
      dt           ! Model timestep duration         [s]

    real( kind = time_precision ), intent(in) :: &
      time_current   ! Current time     [s]

    integer, intent(in) :: &
      hydromet_dim

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      exner,           & ! Exner function                           [-]
      rho,             & ! Density on thermodynamic levels          [kg/m^3]
      rcm,             & ! Mean cloud water mixing ratio            [kg/kg]
      cloud_frac,      & ! Cloud fraction                           [-]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: &
      wp2,        & ! Variance of vertical velocity (momentum levels) [m^2/s^2]
      rho_zm,     & ! Density on momentum levels                      [kg/m^3]
      rho_ds_zm,  & ! Dry, static density on momentum levels   [kg/m^3]
      Kh_zm,      & ! Kh Eddy diffusivity on momentum grid            [m^2/s]
      Skw_zm        ! Skewness of w on momentum levels                [-]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet_mc    ! Microphysics tendency for mean hydrometeors  [units/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      Ncm_mc, & ! Microphysics tendency for Ncm                     [num/kg/s]
      Lscale    ! Length-scale                                      [m]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    logical, intent(in) :: &
      l_upwind_xm_ma ! This flag determines whether we want to use an upwind differencing
                     ! approximation rather than a centered differencing for turbulent or
                     ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.

    ! Input/Output Variables for advance_microphys
    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet,        & ! Hydrometeor mean, <h_m> (thermo. levels)    [units]
      hydromet_vel_zt    ! Mean hydrometeor sed. vel. on thermo. levs. [m/s]

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(in) :: &
      hydrometp2,      & ! Variance of hydrometeor (overall) (m-levs.) [units^2]
      K_hm               ! hm eddy diffusivity on momentum grid        [m^2/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water          [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp.   [K/s]

    ! Output Variables for advance_microphys
    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)   [(m/s)units]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: &
      wpNcp          ! Covariance < w'N_c' > (momentum levels)   [(m/s)(num/kg)]

    type(err_info_type), intent(in) :: &
      err_info        ! err_info struct containing err_code and err_header

    !---------------------------- Begin Code ----------------------------

    if ( clubb_at_least_debug_level_api( 0 ) ) then
       if ( any(err_info%err_code == clubb_fatal_error) ) then

          write(fstderr,*) "Error in advance_microphys"

          write(fstderr,*) "Intent(in)"

          write(fstderr,*) "dt = ", dt
          write(fstderr,*) "time_current = ", time_current

          write(fstderr,*) "wm_zt = ", wm_zt
          write(fstderr,*) "wp2 = ", wp2
          write(fstderr,*) "exner = ", exner
          write(fstderr,*) "rho = ", rho
          write(fstderr,*) "rho_zm = ", rho_zm
          write(fstderr,*) "rcm = ", rcm
          write(fstderr,*) "cloud_frac = ", cloud_frac
          write(fstderr,*) "Kh_zm = ", Kh_zm
          write(fstderr,*) "Skw_zm = ", Skw_zm
          write(fstderr,*) "rho_ds_zm = ", rho_ds_zm
          write(fstderr,*) "rho_ds_zt = ", rho_ds_zt
          write(fstderr,*) "invrs_rho_ds_zt = ", invrs_rho_ds_zt

          write(fstderr,*) "hydromet_mc = ", hydromet_mc

          write(fstderr,*) "Ncm_mc = ", Ncm_mc
          write(fstderr,*) "Lscale = ", Lscale

          write(fstderr,*) "hydromet_vel_covar_zt_impc = ", &
                           hydromet_vel_covar_zt_impc
          write(fstderr,*) "hydromet_vel_covar_zt_expc = ", &
                           hydromet_vel_covar_zt_expc

          write(fstderr,*) "clubb_params = ", clubb_params
          write(fstderr,*) "nu_vert_res_dep%nu_hm = ", nu_vert_res_dep%nu_hm(1)
          write(fstderr,*) "l_upwind_xm_ma = ", l_upwind_xm_ma

          write(fstderr,*) "Intent(inout)"

          write(fstderr,*) "hydromet = ", hydromet
          write(fstderr,*) "hydromet_vel_zt = ", hydromet_vel_zt
          write(fstderr,*) "hydrometp2 = ", hydrometp2
          write(fstderr,*) "K_hm = ", K_hm

          write(fstderr,*) "Ncm = ", Ncm
          write(fstderr,*) "Nc_in_cloud = ", Nc_in_cloud

          write(fstderr,*) "rvm_mc = ", rvm_mc
          write(fstderr,*) "thlm_mc = ", thlm_mc

          write(fstderr,*) "Intent(out)" 

          write(fstderr,*) "wphydrometp = ", wphydrometp

          write(fstderr,*) "wpNcp = ", wpNcp

        endif
    endif


    return

  end subroutine write_adv_micro_errors

!===============================================================================

end module advance_microphys_module

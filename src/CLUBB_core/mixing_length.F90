!-----------------------------------------------------------------------
! $Id: mixing_length.F90 8664 2018-05-10 20:21:35Z huebler@uwm.edu $
!
! Mixing distances and turbulence dissipation time scales for CLUBB.
!
! Description:
!   The mixing length estimates how far turbulent air motions exchange heat,
!   water, and momentum. This module computes that distance and the damping
!   time scales used by the moment equations. Lscale_method selects the
!   legacy buoyant-parcel method (1), length diagnosed from a time scale
!   (2), or the coherence-weighted parcel method (3).
!
! Usage:
!   Normal timestep:
!     - calc_Lscale(): select the method, compute lengths and inverse damping
!       times, and accumulate registered mixing-length/time-scale statistics.
!
!   Components, when only part of the calculation is needed:
!     - set_Lscale_max(): set the per-column length cap from host grid spacing,
!       or use the standalone cap of 100000 m.
!     - calc_Lscale_directly(): compute method 1 or 3 from thermodynamic state
!       and PDF (probability density function) properties. The source constant
!       l_avg_Lscale controls the optional perturbed-launch calculation.
!     - diagnose_Lscale_from_tau(): compute method 2 from background, surface,
!       shear, and stability contributions to inverse damping times.
!
! Understanding:
!   Fields are shaped (ngrdcol,nzt) on thermodynamic levels or (ngrdcol,nzm)
!   on momentum levels. Lengths are in m, time scales in s, inverse time scales
!   in 1/s, temperatures in K, and water mixing ratios in kg/kg. Column
!   parameters are constant with height. Callers supply initialized grid, PDF,
!   statistics, and error structures; no parcel state is kept between calls.
!   stats and err_info are updated in place by the public calculation routines.
!
!   compute_mixing_length() and compute_coherent_mixing_length() are private
!   parcel kernels. Both account for entrainment and moist buoyancy along
!   upward and downward paths. Method 1 stops at the first energy exhaustion.
!   Method 3 integrates a survival weight, includes the first layer with mean
!   coherence <= 0.01, and then stops. Its source constant overwrite_mode
!   chooses local reaches, maximum-endpoint extension, or arrival-weighted
!   averaging; arrival-weighted averaging is selected. This choice is distinct
!   from the runtime Lscale_method flag.
!
!   Method 3 uses sat_mixrat_liq_flatau() directly for its internal parcel
!   saturation; method 1 uses saturation_formula. In the parcel methods only
!   the combined length is capped and copied to the upper boundary; directional
!   lengths retain separate boundary values. Method 2 returns zero directional
!   lengths. See doc/coherence_mixing_length.md
!   for equations, parameter meanings, and configuration examples.
!
! References:
!   Golaz et al. (2002), J. Atmos. Sci., 59, 3540--3551: parcel-method background.
!   Guo et al. (2021, JAMES): length diagnosed from the dissipation time scale.
!===============================================================================
module mixing_length

  implicit none

  private ! Default Scope

  public :: calc_Lscale, &
            calc_Lscale_directly,  &
            diagnose_Lscale_from_tau, &
            set_Lscale_max

  contains

!===============================================================================

  subroutine calc_Lscale( nzm, nzt, ngrdcol, gr, l_implemented, host_dx, host_dy, &
                          p_in_Pa, exner, rtm, thlm, thvm, thlp2, rtp2, rtpthlp, &
                          pdf_params, em, thv_ds_zt, lmin, &
                          upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac, &
                          ufmin, tau_const, sfc_elevation, clubb_params, &
                          saturation_formula, l_Lscale_plume_centered, &
                          Lscale_method, &
                          l_e3sm_config, &
                          l_smooth_Heaviside_tau_wpxp, &
                          l_modify_limiters_for_cnvg_test, &
                          l_use_invrs_tau_N2_iso, &
                          brunt_vaisala_freq_sqd_smth, &
                          stats, err_info, &
                          invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm, &
                          invrs_tau_wp3_zt, &
                          invrs_tau_C1_zm, invrs_tau_C4_zm, &
                          invrs_tau_C6_zm, invrs_tau_C14_zm, &
                          tau_max_zm, tau_max_zt, tau_zm, &
                          Lscale, Lscale_zm, Lscale_up, Lscale_down )

    ! Description:
    !   Compute CLUBB's mixing length and dissipation time scales, including
    !   Richardson-number-dependent damping terms and Lscale-related stats.
    !   Read thermodynamic fields on nzt levels, moments on nzm levels, and
    !   per-column tunable parameters; Lscale_method selects methods 1--3.
    !   Return lengths [m], time scales [s], and inverse times [1/s], and update
    !   stats and err_info. Grid and PDF state must be initialized by the caller.
    !
    ! References:
    !   Guo et al. (2021, JAMES) for the diagnosed Lscale-from-tau option.
    !-----------------------------------------------------------------------

    use advance_helper_module, only: &
        calc_Ri_zm, &
        smooth_max

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        em_min, &
        fstderr, &
        min_max_smth_mag, &
        one, &
        rt_tol, &
        thl_tol, &
        zero_threshold

    use parameter_indices, only: &
        nparams, &
        itaumax

    use model_flags, only: &
        Lscale_tau

    use grid_class, only: &
        grid, &
        zm2zt_api, &
        zm2zt2zm, &
        zt2zm_api

    use pdf_parameter_module, only: &
        pdf_parameter

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    use error_code, only: &
        clubb_at_least_debug_level_api, &
        clubb_fatal_error

    use err_info_type_module, only: &
        err_info_type

    use code_timer_module, only: &
        timer_start, &
        timer_stop

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm,    & ! Number of momentum levels
      nzt,    & ! Number of thermodynamic levels
      ngrdcol   ! Number of grid columns

    type(grid), intent(in) :: &
      gr   ! Grid structure

    logical, intent(in) :: &
      l_implemented,                    & ! True when current configuration is implemented
      l_Lscale_plume_centered,          & ! Whether to use plume-centered Lscale calculation
      l_e3sm_config,                    & ! Whether E3SM configuration options are active
      l_smooth_Heaviside_tau_wpxp,      & ! Whether to smooth Heaviside functions in wpxp tau
      l_modify_limiters_for_cnvg_test,  & ! Whether to modify limiters for convergence tests
      l_use_invrs_tau_N2_iso              ! Whether to use isotropic N2 inverse tau

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      host_dx,        & ! Host-model grid spacing in x direction [m]
      host_dy,        & ! Host-model grid spacing in y direction [m]
      upwp_sfc,       & ! Surface value of u'w' [m^2/s^2]
      vpwp_sfc,       & ! Surface value of v'w' [m^2/s^2]
      sfc_elevation     ! Surface elevation [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      p_in_Pa,           & ! Pressure [Pa]
      exner,             & ! Exner function [-]
      rtm,               & ! Total water mixing ratio [kg/kg]
      thlm,              & ! Liquid water potential temperature [K]
      thvm,              & ! Virtual potential temperature [K]
      thv_ds_zt,         & ! Reference-state virtual potential temperature [K]
      ice_supersat_frac    ! Ice supersaturation fraction [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      thlp2,                       & ! Variance of thl [K^2]
      rtp2,                        & ! Variance of rt [(kg/kg)^2]
      rtpthlp,                     & ! Covariance of rt and thl [K kg/kg]
      em,                          & ! Turbulent kinetic energy [m^2/s^2]
      ddzt_umvm_sqd,               & ! Squared vertical shear of horizontal wind [s^-2]
      brunt_vaisala_freq_sqd_smth    ! Smoothed Brunt-Vaisala frequency squared [s^-2]

    real( kind = core_rknd ), intent(in) :: &
      lmin,      & ! Minimum mixing length [m]
      ufmin,     & ! Minimum velocity scale [m/s]
      tau_const    ! Constant eddy dissipation time scale [s]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params   ! Tunable CLUBB parameters

    integer, intent(in) :: &
      saturation_formula, & ! Saturation vapor pressure formula selector
      Lscale_method         ! Mixing-length formulation

    type(pdf_parameter), intent(in) :: &
      pdf_params   ! PDF parameters

    !----------------------- Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: &
      stats   ! Statistics accumulator

    type(err_info_type), intent(inout) :: &
      err_info   ! Error information

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      invrs_tau_zm,      & ! Inverse time scale on momentum levels [1/s]
      invrs_tau_xp2_zm,  & ! Inverse time scale for xp2 terms on momentum levels [1/s]
      invrs_tau_C1_zm,   & ! Inverse time scale used for C1 term on momentum levels [1/s]
      invrs_tau_C4_zm,   & ! Inverse time scale used for C4 terms on momentum levels [1/s]
      invrs_tau_C6_zm,   & ! Inverse time scale used for C6 term on momentum levels [1/s]
      invrs_tau_C14_zm,  & ! Inverse time scale used for C14 terms on momentum levels [1/s]
      tau_max_zm,        & ! Maximum allowable eddy dissipation time scale on momentum levels [s]
      tau_zm,            & ! Eddy dissipation time scale on momentum levels [s]
      Lscale_zm            ! Turbulent mixing length on momentum levels [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      invrs_tau_zt,      & ! Inverse time scale on thermodynamic levels [1/s]
      invrs_tau_wp3_zt,  & ! Inverse time scale for wp3 terms on thermodynamic levels [1/s]
      tau_max_zt,        & ! Maximum allowable eddy dissipation time scale on thermodynamic levels [s]
      Lscale,            & ! Turbulent mixing length on thermodynamic levels [m]
      Lscale_up,         & ! Upward component of turbulent mixing length [m]
      Lscale_down          ! Downward component of turbulent mixing length [m]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      Lscale_max   ! Maximum allowable mixing length [m]

    logical, parameter :: l_smooth_min_max = .false.  ! whether to apply smooth min
                                                      ! and max functions

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      em_zt,       & ! TKE interpolated to thermodynamic levels [m^2/s^2]
      thlp2_zt,    & ! Variance of thl on thermodynamic levels [K^2]
      rtp2_zt,     & ! Variance of rt on thermodynamic levels [(kg/kg)^2]
      rtpthlp_zt,  & ! Covariance of rt and thl on thermodynamic levels [K kg/kg]
      tau_zt         ! Eddy dissipation time scale on thermodynamic levels [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      invrs_tau_sfc,              & ! Surface contribution to inverse tau [1/s]
      invrs_tau_no_N2_zm,         & ! Inverse tau without N2 contribution on momentum levels [1/s]
      invrs_tau_bkgnd,            & ! Background inverse tau [1/s]
      invrs_tau_N2_iso,           & ! Isotropic N2 contribution to inverse tau [1/s]
      invrs_tau_wp2_zm,           & ! Inverse time scale for wp2 terms on momentum levels [1/s]
      invrs_tau_wp3_zm,           & ! Inverse time scale for wp3 terms on momentum levels [1/s]
      invrs_tau_wpxp_zm,          & ! Inverse time scale for wpxp terms on momentum levels [1/s]
      invrs_tau_shear,            & ! Shear contribution to inverse tau [1/s]
      Ri_zm,                      & ! Gradient Richardson number on momentum levels [-]
      ddzt_umvm_sqd_clipped,      & ! Clipped squared vertical shear of horizontal wind [s^-2]
      brunt_vaisala_freq_clipped    ! Clipped Brunt-Vaisala frequency [1/s]

    integer :: &
      i,  & ! Grid-column loop index
      k     ! Vertical-level loop index

    !----------------------------- Begin Code ------------------------------

    call timer_start( "calc_Lscale" )

    !$acc enter data create( Lscale_max, em_zt, thlp2_zt, rtp2_zt, rtpthlp_zt, &
    !$acc                    tau_zt, &
    !$acc                    invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
    !$acc                    invrs_tau_N2_iso, invrs_tau_wp2_zm, &
    !$acc                    invrs_tau_wp3_zm, invrs_tau_wpxp_zm, &
    !$acc                    invrs_tau_shear, Ri_zm, &
    !$acc                    ddzt_umvm_sqd_clipped, brunt_vaisala_freq_clipped )

    ! Determine the maximum allowable value for Lscale (in meters).
    call set_Lscale_max( ngrdcol, l_implemented, host_dx, host_dy, & ! In
                         Lscale_max )                                ! Out

    ! Interpolate rtpthlp here so the diagnostic stats value is generated at
    ! the same point as the Lscale calculation that consumes it.
    rtpthlp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, rtpthlp(:,:) )

    if ( stats%l_sample ) then
      !$acc update host( rtpthlp_zt )
      call stats_update( "rtpthlp_zt", rtpthlp_zt, stats )
    end if

    ! Calculate Richardson number Ri_zm.
    if ( l_modify_limiters_for_cnvg_test ) then

      call calc_Ri_zm( nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd, &
                       0.0_core_rknd, 1.0e-12_core_rknd, Ri_zm )

      Ri_zm = zm2zt2zm( nzm, nzt, ngrdcol, gr, Ri_zm )

    else

      if ( l_smooth_min_max ) then

        brunt_vaisala_freq_clipped = smooth_max( nzm, ngrdcol, 1.0e-7_core_rknd, &
                                                 brunt_vaisala_freq_sqd_smth, &
                                                 1.0e-4_core_rknd * min_max_smth_mag )

        ddzt_umvm_sqd_clipped = smooth_max( nzm, ngrdcol, ddzt_umvm_sqd, &
                                            1.0e-7_core_rknd, &
                                            1.0e-6_core_rknd * min_max_smth_mag )

        call calc_Ri_zm( nzm, ngrdcol, brunt_vaisala_freq_clipped, &
                         ddzt_umvm_sqd_clipped, 0.0_core_rknd, 0.0_core_rknd, Ri_zm )

      else

        call calc_Ri_zm( nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd, &
                         1.0e-7_core_rknd, 1.0e-7_core_rknd, Ri_zm )

      end if

    end if

    if ( stats%l_sample ) then
      !$acc update host( Ri_zm )
      call stats_update( "Ri_zm", Ri_zm, stats )
    end if

    if ( Lscale_method /= Lscale_tau ) then
      ! Compute Lscale first using the buoyant parcel calculation.

      ! Interpolate thlp2 and rtp2 to thermodynamic levels.
      thlp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, thlp2(:,:), &
                                  thl_tol**2 )  ! Positive def. quantity
      rtp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, rtp2(:,:), &
                                  rt_tol**2 )   ! Positive def. quantity

      call calc_Lscale_directly ( ngrdcol, nzm, nzt, gr,                             & ! In
                                  l_implemented, p_in_Pa,                            & ! In
                                  exner, rtm, thlm, thvm,                            & ! In
                                  rtp2_zt, thlp2_zt, rtpthlp_zt, pdf_params,         & ! In
                                  em, thv_ds_zt, Lscale_max, lmin,                   & ! In
                                  clubb_params,                                      & ! In
                                  saturation_formula,                                & ! In
                                  l_Lscale_plume_centered,                           & ! In
                                  Lscale_method,                                     & ! In
                                  stats, err_info,                                   & ! InOut
                                  Lscale, Lscale_up, Lscale_down )                     ! Out

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr,*) "Error calling calc_Lscale_directly in calc_Lscale"
          !$acc exit data delete( Lscale_max, em_zt, thlp2_zt, rtp2_zt, rtpthlp_zt, &
          !$acc                   tau_zt, &
          !$acc                   invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
          !$acc                   invrs_tau_N2_iso, invrs_tau_wp2_zm, &
          !$acc                   invrs_tau_wp3_zm, invrs_tau_wpxp_zm, &
          !$acc                   invrs_tau_shear, Ri_zm, &
          !$acc                   ddzt_umvm_sqd_clipped, brunt_vaisala_freq_clipped )
          return
        end if
      end if

      ! Calculate CLUBB's turbulent eddy-turnover time scale as CLUBB's
      ! length scale divided by a velocity scale.
      em_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, em(:,:), em_min )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          tau_zt(i,k) = min( Lscale(i,k) / sqrt( em_zt(i,k) ), clubb_params(i,itaumax) )
        end do
      end do
      !$acc end parallel loop

      Lscale_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, Lscale(:,:), zero_threshold )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_zm(i,k) = min( Lscale_zm(i,k) / sqrt( max( em_min, em(i,k) ) ),  &
                             clubb_params(i,itaumax) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_zm(i,k)      = one / tau_zm(i,k)
          invrs_tau_wp2_zm(i,k)  = invrs_tau_zm(i,k)
          invrs_tau_xp2_zm(i,k)  = invrs_tau_zm(i,k)
          invrs_tau_wpxp_zm(i,k) = invrs_tau_zm(i,k)
          invrs_tau_wp3_zm(i,k)  = invrs_tau_zm(i,k)
          tau_max_zm(i,k) = clubb_params(i,itaumax)
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          invrs_tau_zt(i,k)     = one / tau_zt(i,k)
          invrs_tau_wp3_zt(i,k) = invrs_tau_zt(i,k)
          tau_max_zt(i,k) = clubb_params(i,itaumax)
        end do
      end do
      !$acc end parallel loop

    else
      ! Diagnose simple tau and Lscale.

      call diagnose_Lscale_from_tau( nzm, nzt, ngrdcol, gr,                       & ! In
                        upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,                        & ! In
                        ice_supersat_frac,                                        & ! In
                        em,                                                       & ! In
                        ufmin, tau_const,                                         & ! In
                        sfc_elevation, Lscale_max,                                & ! In
                        clubb_params,                                             & ! In
                        stats,                                                    & ! InOut
                        l_e3sm_config,                                            & ! In
                        l_smooth_Heaviside_tau_wpxp,                              & ! In
                        brunt_vaisala_freq_sqd_smth, Ri_zm,                       & ! In
                        err_info,                                                 & ! InOut
                        invrs_tau_zt, invrs_tau_zm,                               & ! Out
                        invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd,       & ! Out
                        invrs_tau_shear, invrs_tau_N2_iso,                        & ! Out
                        invrs_tau_wp2_zm, invrs_tau_xp2_zm,                       & ! Out
                        invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm,    & ! Out
                        tau_max_zm, tau_max_zt, tau_zm, tau_zt,                   & ! Out
                        Lscale, Lscale_up, Lscale_down )                            ! Out

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr,*) "Error calling diagnose_Lscale_from_tau in calc_Lscale"
          !$acc exit data delete( Lscale_max, em_zt, thlp2_zt, rtp2_zt, rtpthlp_zt, &
          !$acc                   tau_zt, &
          !$acc                   invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
          !$acc                   invrs_tau_N2_iso, invrs_tau_wp2_zm, &
          !$acc                   invrs_tau_wp3_zm, invrs_tau_wpxp_zm, &
          !$acc                   invrs_tau_shear, Ri_zm, &
          !$acc                   ddzt_umvm_sqd_clipped, brunt_vaisala_freq_clipped )
          return
        end if
      end if

    end if

    Lscale_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, Lscale(:,:) )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_C6_zm(i,k) = invrs_tau_wpxp_zm(i,k)
        invrs_tau_C1_zm(i,k) = invrs_tau_wp2_zm(i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_C14_zm(i,k) = invrs_tau_wp2_zm(i,k)
      end do
    end do
    !$acc end parallel loop

    if ( Lscale_method /= Lscale_tau .and. l_use_invrs_tau_N2_iso) then
      write(fstderr,*) "Error! l_use_invrs_tau_N2_iso is not used when "// &
                       "Lscale_method is not the tau method."// &
                       "If you want to use Lscale code, go to file "// &
                       "src/CLUBB_core/mixing_length.F90 and "// &
                       "change l_use_invrs_tau_N2_iso to false"
      error stop
    end if

    if ( .not. l_use_invrs_tau_N2_iso ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_C4_zm(i,k) = invrs_tau_wp2_zm(i,k)
        end do
      end do
      !$acc end parallel loop
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_C4_zm(i,k) = invrs_tau_N2_iso(i,k)
        end do
      end do
      !$acc end parallel loop
    end if

    if ( stats%l_sample ) then
      !$acc update host( Lscale, Lscale_up, Lscale_down, tau_zm, &
      !$acc              invrs_tau_zm, invrs_tau_xp2_zm, invrs_tau_wp2_zm, &
      !$acc              invrs_tau_wpxp_zm, invrs_tau_wp3_zm )
      call stats_update( "Lscale", Lscale, stats )
      call stats_update( "Lscale_up", Lscale_up, stats )
      call stats_update( "Lscale_down", Lscale_down, stats )
      call stats_update( "tau_zm", tau_zm, stats )
      call stats_update( "invrs_tau_zm", invrs_tau_zm, stats )
      call stats_update( "invrs_tau_xp2_zm", invrs_tau_xp2_zm, stats )
      call stats_update( "invrs_tau_wp2_zm", invrs_tau_wp2_zm, stats )
      call stats_update( "invrs_tau_wpxp_zm", invrs_tau_wpxp_zm, stats )
      call stats_update( "invrs_tau_wp3_zm", invrs_tau_wp3_zm, stats )
      !$acc update host( tau_zt )
      call stats_update( "tau_zt", tau_zt, stats )

      if ( Lscale_method == Lscale_tau ) then
        !$acc update host( invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
        !$acc              invrs_tau_sfc, invrs_tau_shear )
        call stats_update( "invrs_tau_no_N2_zm", invrs_tau_no_N2_zm, stats )
        call stats_update( "invrs_tau_bkgnd", invrs_tau_bkgnd, stats )
        call stats_update( "invrs_tau_sfc", invrs_tau_sfc, stats )
        call stats_update( "invrs_tau_shear", invrs_tau_shear, stats )
      end if
    end if

    !$acc exit data delete( Lscale_max, em_zt, thlp2_zt, rtp2_zt, rtpthlp_zt, &
    !$acc                   tau_zt, &
    !$acc                   invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
    !$acc                   invrs_tau_N2_iso, invrs_tau_wp2_zm, &
    !$acc                   invrs_tau_wp3_zm, invrs_tau_wpxp_zm, &
    !$acc                   invrs_tau_shear, Ri_zm, &
    !$acc                   ddzt_umvm_sqd_clipped, brunt_vaisala_freq_clipped )

    call timer_stop( "calc_Lscale" )

    return

  end subroutine calc_Lscale

  !=============================================================================

  subroutine set_Lscale_max( ngrdcol, l_implemented, host_dx, host_dy, &
                              Lscale_max )

    ! Description:
    !   This subroutine sets the value of Lscale_max, which is the maximum
    !   allowable value of Lscale.  For standard CLUBB, it is set to a very large
    !   value so that Lscale will not be limited.  However, when CLUBB is running
    !   as part of a host model, the value of Lscale_max is dependent on the size
    !   of the host model's horizontal grid spacing.  The smaller the host model's
    !   horizontal grid spacing, the smaller the value of Lscale_max.  When Lscale
    !   is limited to a small value, the value of time-scale Tau is reduced, which
    !   in turn produces greater damping on CLUBB's turbulent parameters.  This
    !   is the desired effect on turbulent parameters for a host model with small
    !   horizontal grid spacing, for small areas usually contain much less
    !   variation in meteorological quantities than large areas.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    !----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      ngrdcol

    logical, intent(in) :: &
      l_implemented     ! True if CLUBB is being implemented and run in a host model

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      host_dx, & ! Host model's east-west horizontal grid spacing     [m]
      host_dy    ! Host model's north-south horizontal grid spacing   [m]

    !----------------------- Output Variable -----------------------
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: &
      Lscale_max    ! Maximum allowable value for Lscale   [m]

    !----------------------- Local Variable -----------------------
    integer :: i

    !----------------------- Begin Code-----------------------

    ! Determine the maximum allowable value for Lscale (in meters).
    if ( l_implemented ) then
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        Lscale_max(i) = 0.25_core_rknd * min( host_dx(i), host_dy(i) )
      end do
      !$acc end parallel loop
    else
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        Lscale_max(i) = 1.0e5_core_rknd
      end do
      !$acc end parallel loop
    end if

    return
  end subroutine set_Lscale_max

  !=============================================================================

  subroutine compute_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm, &
                                    rtm, em, Lscale_max, p_in_Pa, &
                                    exner, thv_ds, Lscale_mu_entrain, lmin, &
                                    saturation_formula, &
                                    l_implemented, &
                                    err_info, &
                                    Lscale, Lscale_up, Lscale_down )

    ! Description:
    !   Larson's 5th moist, nonlocal length scale
    !
    ! References:
    !   Section 3b ( /Eddy length formulation/ ) of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !
    ! Notes:
    !
    !   The equation for the rate of change of theta_l and r_t of the parcel with
    !   respect to height, due to entrainment, is:
    !
    !           d(thl_par)/dz = - Lscale_mu_entrain * ( thl_parcel - thl_environment );
    !
    !           d(rt_par)/dz = - Lscale_mu_entrain * ( rt_parcel - rt_environment );
    !
    !   where Lscale_mu_entrain is the entrainment rate,
    !   such that:
    !
    !           Lscale_mu_entrain = (1/m)*(dm/dz);
    !
    !   where m is the mass of the parcel.  The value of Lscale_mu_entrain is set to be a
    !   constant.
    !
    !   The differential equations are solved for given the boundary condition
    !   and given the fact that the value of thl_environment and rt_environment
    !   are treated as changing linearly for a parcel of air from one grid level
    !   to the next.
    !
    !   For the special case where entrainment rate, Lscale_mu_entrain, is set to 0,
    !   thl_parcel and rt_parcel remain constant
    !
    !
    !   The equation for Lscale_up is:
    !
    !       INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);
    !
    !   and for Lscale_down
    !
    !       INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);
    !
    !   where thv_par is theta_v of the parcel, thvm is the mean
    !   environmental value of theta_v, z_i is the altitude that the parcel
    !   started from, and em is the mean value of TKE at
    !   altitude z_i (which gives the parcel its initial boost)
    !
    !   The increment of CAPE (convective air potential energy) for any two
    !   successive vertical levels is:
    !
    !       Upwards:
    !           CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
    !
    !       Downwards:
    !           CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz
    !
    !   Thus, the derivative of CAPE with respect to height is:
    !
    !           dCAPE/dz = g * ( thv_par - thvm ) / thvm.
    !
    !   A purely trapezoidal rule is used between levels, and is considered
    !   to vary linearly at all altitudes.  Thus, dCAPE/dz is considered to be
    !   of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
    !   where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 )
    !
    !   The integral is evaluated to find the CAPE increment between two
    !   successive vertical levels.  The result either adds to or depletes
    !   from the total amount of energy that keeps the parcel ascending/descending.
    !
    !
    ! IMPORTANT NOTE:
    !   This subroutine has been optimized by adding precalculations, rearranging
    !   equations to avoid divides, and modifying the algorithm entirely.
    !       -Gunther Huebler
    !
    !   The algorithm previously used looped over every grid level, following a
    !   a parcel up from its initial grid level to its max. The very nature of this
    !   algorithm is an N^2
    !--------------------------------------------------------------------------------

    ! Lscale_mu_entrain = (1/M) dM/dz > 0.  Lscale_mu_entrain=0 for no entrainment.
    ! Siebesma recommends Lscale_mu_entrain=2e-3, although most schemes use Lscale_mu_entrain=1e-4
    ! When Lscale_mu_entrain was fixed, we used the value Lscale_mu_entrain = 6.e-4

    use constants_clubb, only:  &  ! Variable(s)
        Cp,             & ! Dry air specific heat at constant pressure [J/kg/K]
        Rd,             & ! Dry air gas constant                       [J/kg/K]
        ep,             & ! Rd / Rv                                    [-]
        ep1,            & ! (1-ep)/ep                                  [-]
        ep2,            & ! 1/ep                                       [-]
        Lv,             & ! Latent heat of vaporiztion                 [J/kg/K]
        grav,           & ! Gravitational acceleration                 [m/s^2]
        fstderr,        &
        zero_threshold, &
        eps,            &
        one_half,       &
        one,            &
        two,            &
        zero

    use grid_class, only:  &
        grid, & ! Type
        zm2zt_api ! Procedure(s)

    use numerical_check, only:  &
        length_check ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error              ! Constant

    use saturation, only:  &
        sat_mixrat_liq_api ! Procedure(s)

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    ! Constant Parameters
    real( kind = core_rknd ), parameter ::  &
      zlmin = 0.1_core_rknd, & ! Minimum value for Lscale [m]
      Lscale_sfclyr_depth = 500._core_rknd ! [m]

    !--------------------------------- Input Variables ---------------------------------
    integer, intent(in) :: &
      nzm, & ! Number of momentum levels
      nzt, & ! Number of thermodynamic levels
      ngrdcol ! Number of grid columns

    type (grid), intent(in) :: &
      gr ! Grid heights, interpolation weights, and orientation

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      thvm,    & ! Virtual potential temp. on themodynamic level  [K]
      thlm,    & ! Liquid potential temp. on themodynamic level   [K]
      rtm,     & ! Total water mixing ratio on themodynamic level [kg/kg]
      exner,   & ! Exner function on thermodynamic level          [-]
      p_in_Pa, & ! Pressure on thermodynamic level                [Pa]
      thv_ds     ! Dry, base-state theta_v on thermodynamic level [K]
    ! Note:  thv_ds used as a reference theta_l here

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      em         ! em = 3/2 * w'^2; on momentum level             [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      Lscale_max ! Maximum allowable value for Lscale             [m]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      Lscale_mu_entrain      ! Fractional entrainment rate per unit altitude [1/m]

    real( kind = core_rknd ), intent(in) :: &
      lmin    ! CLUBB tunable parameter lmin

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    logical, intent(in) :: &
      l_implemented ! True if CLUBB is being implemented and run in a host model

    !--------------------------------- InOut Variables ---------------------------------
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    !--------------------------------- Output Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    !--------------------------------- Local Variables ---------------------------------

    integer :: i, j, k, start_index, j_zm, jp1_zm, k_zm, kp1_zm

    real( kind = core_rknd ) :: tke, CAPE_incr

    real( kind = core_rknd ) :: dCAPE_dz_j, dCAPE_dz_j_minus_1, dCAPE_dz_j_plus_1

    ! Temporary 2D arrays to store calculations to speed runtime
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
        grav_on_thvm, &
        Lv_coef, &
        thl_par_j_precalc, &
        rt_par_j_precalc, &
        tl_par_1, &
        rt_par_1, &
        rsatl_par_1, &
        thl_par_1, &
        dCAPE_dz_1, &
        s_par_1, &
        rc_par_1, &
        CAPE_incr_1, &
        thv_par_1

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
        exp_mu_dzm, &
        invrs_dzm_on_mu, &
        entrain_coef

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
        tke_i

    ! Minimum value for Lscale that will taper off with height
    real( kind = core_rknd ) :: lminh

    ! Parcel quantities at grid level j
    real( kind = core_rknd ) :: thl_par_j, rt_par_j, rc_par_j, thv_par_j

    ! Used in latent heating calculation
    real( kind = core_rknd ) :: tl_par_j, rsatl_par_j, s_par_j

    ! Variables to make L nonlocal
    real( kind = core_rknd ) :: Lscale_up_max_alt, Lscale_down_min_alt

    ! Variables used to precalculate values
    real( kind = core_rknd ) :: &
        Lv2_coef, &
        tl_par_j_sqd, &
        invrs_dCAPE_diff, &
        invrs_Lscale_sfclyr_depth

    ! --------------------------------- Begin Code ---------------------------------

    !$acc enter data create( exp_mu_dzm, invrs_dzm_on_mu, grav_on_thvm, Lv_coef, &
    !$acc                    entrain_coef, thl_par_j_precalc, rt_par_j_precalc, &
    !$acc                    tl_par_1, rt_par_1, rsatl_par_1, thl_par_1, dCAPE_dz_1, &
    !$acc                    s_par_1, rc_par_1, CAPE_incr_1, thv_par_1, tke_i )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      if ( abs(Lscale_mu_entrain(i)) < eps ) then
        ! Error in grid column i -> set ith entry to clubb_fatal_error
        err_info%err_code(i) = clubb_fatal_error
        write(fstderr,*) err_info%err_header(i)
        write(fstderr,*) "Lscale_mu_entrain = ", Lscale_mu_entrain(i)
      end if
    end do
    !$acc end parallel loop

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      !$acc update host( err_info%err_code )
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) "Entrainment rate Lscale_mu_entrain cannot be 0"
        write(fstderr, *) "Fatal error in subroutine compute_mixing_length"
        return
      end if
    end if

    ! Calculate initial turbulent kinetic energy for each grid level
    tke_i = zm2zt_api( nzm, nzt, ngrdcol, gr, em )

    ! Initialize arrays and precalculate values for computational efficiency
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        ! Initialize up and down arrays
        Lscale_up(i,k) = zlmin
        Lscale_down(i,k) = zlmin

        ! Precalculate values to avoid unnecessary calculations later
        grav_on_thvm(i,k) = grav / thvm(i,k)
        Lv_coef(i,k) = Lv / ( exner(i,k) * cp ) - ep2 * thv_ds(i,k)

      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do k = 1, nzm

        ! Precalculate values to avoid unnecessary calculations later
        exp_mu_dzm(i,k) = exp( -Lscale_mu_entrain(i) * gr%grid_dir * gr%dzm(i,k) )
        invrs_dzm_on_mu(i,k) = ( gr%grid_dir * gr%invrs_dzm(i,k) ) / Lscale_mu_entrain(i)
        entrain_coef(i,k) = ( one - exp_mu_dzm(i,k) ) * invrs_dzm_on_mu(i,k)

      end do
    end do
    !$acc end parallel loop

    ! Precalculations of single values to avoid unnecessary calculations later
    Lv2_coef = ep * Lv**2 / ( Rd * cp )
    invrs_Lscale_sfclyr_depth = one / Lscale_sfclyr_depth


    ! ---------------- Upwards Length Scale Calculation ----------------

    ! Precalculate values for upward Lscale, these are useful only if a parcel can rise
    ! more than one level. They are used in the equations that calculate thl and rt
    ! recursively for a parcel as it ascends

    !$acc parallel loop gang vector collapse(2) default(present)
    do j = gr%k_lb_zt+gr%grid_dir_indx, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          j_zm = j
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          j_zm = j+1
        endif

        thl_par_j_precalc(i,j) = thlm(i,j) - thlm(i,j-gr%grid_dir_indx) * exp_mu_dzm(i,j_zm)  &
                                 - ( thlm(i,j) - thlm(i,j-gr%grid_dir_indx) ) * entrain_coef(i,j_zm)

        rt_par_j_precalc(i,j) = rtm(i,j) - rtm(i,j-gr%grid_dir_indx) * exp_mu_dzm(i,j_zm)  &
                                - ( rtm(i,j) - rtm(i,j-gr%grid_dir_indx) ) * entrain_coef(i,j_zm)

      end do
    end do
    !$acc end parallel loop

    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations.

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    thl_par_1 = zero
    tl_par_1 = zero
    rt_par_1 = zero
    !$acc parallel loop gang vector collapse(2) default(present)
    do j = gr%k_lb_zt+gr%grid_dir_indx, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          j_zm = j
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          j_zm = j+1
        endif

        thl_par_1(i,j) &
        = thlm(i,j) &
          - ( thlm(i,j) - thlm(i,j-gr%grid_dir_indx) ) * entrain_coef(i,j_zm)

        tl_par_1(i,j) = thl_par_1(i,j) * exner(i,j)

        rt_par_1(i,j) &
        = rtm(i,j) &
          - ( rtm(i,j) - rtm(i,j-gr%grid_dir_indx) ) * entrain_coef(i,j_zm)

      end do
    end do
    !$acc end parallel loop


    ! Caclculate initial rsatl for parcels at each grid level

    ! The entire pressure and temperature arrays are passed as
    ! argument and the sub-arrays are choosen using
    ! start_index. This workaround is used to solve
    ! subarray issues with OpenACC.
    ! rsatl_par_1(i,3:) = sat_mixrat_liq_acc( nz-2, ngrdcol, p_in_Pa(i,3:), tl_par_1(i,3:) )
    ! since subarray 3:, the start_index is 3 and it is an optional argument
    start_index = gr%k_lb_zt+gr%grid_dir_indx
    rsatl_par_1 = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, tl_par_1, saturation_formula, &
                                      start_index )

    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do j = gr%k_lb_zt+gr%grid_dir_indx, gr%k_ub_zt, gr%grid_dir_indx

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          j_zm = j
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          j_zm = j+1
        endif

        tl_par_j_sqd = tl_par_1(i,j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
        !
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )
        s_par_1(i,j) = ( rt_par_1(i,j) - rsatl_par_1(i,j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(i,j) )

        rc_par_1(i,j) = max( s_par_1(i,j), zero_threshold )

        ! theta_v of entraining parcel at grid level j
        thv_par_1(i,j) = thl_par_1(i,j) + ep1 * thv_ds(i,j) * rt_par_1(i,j) + Lv_coef(i,j) &
                         * rc_par_1(i,j)


        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(i,j) = grav_on_thvm(i,j) * ( thv_par_1(i,j) - thvm(i,j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(i,j) = one_half * dCAPE_dz_1(i,j) * gr%grid_dir * gr%dzm(i,j_zm)

      end do


      ! Calculate Lscale_up for each grid level. If the TKE from a parcel has not been completely
      ! exhausted by the initial change then continue the exhaustion calculations here for a single
      ! grid level at a time until the TKE is exhausted.

      Lscale_up_max_alt = zero    ! Set initial max value for Lscale_up to 0
      do k = gr%k_lb_zt, gr%k_ub_zt-2*gr%grid_dir_indx, gr%grid_dir_indx

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i,k) + CAPE_incr_1(i,k+gr%grid_dir_indx) > zero ) then

          ! Calculate new TKE for parcel
          tke = tke_i(i,k) + CAPE_incr_1(i,k+gr%grid_dir_indx)

          ! Set j to 2 levels above current Lscale_up level, this is because we've already
          ! determined that the parcel can rise at least 1 full level
          j = k + 2*gr%grid_dir_indx

          ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
          thl_par_j = thl_par_1(i,k+gr%grid_dir_indx)
          rt_par_j  = rt_par_1(i,k+gr%grid_dir_indx)
          dCAPE_dz_j_minus_1 = dCAPE_dz_1(i,k+gr%grid_dir_indx)


          ! Continue change in TKE calculations until it is exhausted or the max grid
          ! level has been reached. j is the next grid level above the level that can
          ! be reached for a parcel starting at level k. If TKE is exhausted in this loop
          ! that means the parcel starting at k cannot reach level j, but has reached j-1
          !
          ! Update for generalized gridding:
          !
          ! Ascending grid: do while ( j < nzt ) -- j is increasing
          ! Descending grid: do while ( j > 1 ) -- j is decreasing
          !
          ! Either way, j is going upwards in altitude.
          !
          ! Reconciling the grids:
          !
          ! Step 1: Use gr%k_ub_zt for the upper boundary:
          ! Ascending grid: do while ( j < gr%k_ub_zt ) -- j is increasing
          ! Descending grid: do while ( j > gr%k_ub_zt ) -- j is decreasing
          !
          ! Step 2: Reconciling the greater than / less than:
          ! Ascending grid: do while ( j < nzt ) -- j is increasing
          ! Descending grid: do while ( -j < -1 ) -- j is increasing
          !
          ! Generalized Expression:
          ! do while ( gr%grid_dir_indx * j < gr%grid_dir_indx * gr%k_ub_zt );
          ! where for an ascending grid, gr%grid_dir_indx = 1 and gr%k_ub_zt = nzt;
          ! and for a descending grid, gr%grid_dir_indx = -1 and gr%k_ub_zt = 1.
          do while ( gr%grid_dir_indx * j < gr%grid_dir_indx * gr%k_ub_zt )

            if ( gr%grid_dir_indx > 0 ) then
              ! Ascending grid
              j_zm = j
            else ! gr%grid_dir_indx < 0
              ! Descending grid
              j_zm = j+1
            endif

            ! thl, rt of parcel are conserved except for entrainment
            !
            ! The values of thl_env and rt_env are treated as changing linearly for a parcel
            ! of air ascending from level j-1 to level j

            ! theta_l of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(thl_par)/dz = - Lscale_mu_entrain * ( thl_par - thl_env )
            thl_par_j = thl_par_j_precalc(i,j) + thl_par_j * exp_mu_dzm(i,j_zm)


            ! r_t of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(rt_par)/dz = - Lscale_mu_entrain * ( rt_par - rt_env )
            rt_par_j = rt_par_j_precalc(i,j) + rt_par_j * exp_mu_dzm(i,j_zm)


            ! Include effects of latent heating on Lscale_up 6/12/00
            ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
            ! Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(i,j)

            rsatl_par_j = sat_mixrat_liq_api( p_in_Pa(i,j), tl_par_j, saturation_formula )

            tl_par_j_sqd = tl_par_j**2

            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
            !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )
            ! and SD's beta (eqn. 8),
            !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
            !
            ! Simplified by multiplying top and bottom by tl^2 to avoid a
            ! divide and precalculating ep * Lv**2 / ( Rd * cp )
            s_par_j = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
                      / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

            rc_par_j = max( s_par_j, zero_threshold )

            ! theta_v of entraining parcel at grid level j
            thv_par_j = thl_par_j + ep1 * thv_ds(i,j) * rt_par_j  &
                        + Lv_coef(i,j) * rc_par_j

            ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
            dCAPE_dz_j = grav_on_thvm(i,j) * ( thv_par_j - thvm(i,j) )

            ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
            ! Trapezoidal estimate between grid levels j and j-1
            CAPE_incr = one_half * ( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) &
                        * gr%grid_dir * gr%dzm(i,j_zm)

            ! Exit loop early if tke has been exhaused between level j and j+1
            if ( tke + CAPE_incr <= zero ) then
                exit
            end if

            ! Save previous dCAPE value for next cycle
            dCAPE_dz_j_minus_1 = dCAPE_dz_j

            ! Caclulate new TKE and increment j
            tke = tke + CAPE_incr
            j = j + gr%grid_dir_indx

          enddo


          ! Add full grid level thickness for each grid level that was passed without the TKE
          ! being exhausted, difference between starting level (k) and last level passed (j-1)
          Lscale_up(i,k) = Lscale_up(i,k) + gr%zt(i,j-gr%grid_dir_indx) - gr%zt(i,k)


          if ( gr%grid_dir_indx * j < gr%grid_dir_indx * gr%k_ub_zt ) then

            ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
            ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.

            if ( abs( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) * 2 <= &
                 abs( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * eps ) then

              ! Special case where dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) = 0
              ! Find the remaining distance z - z_0 that it takes to
              ! exhaust the remaining TKE

              Lscale_up(i,k) = Lscale_up(i,k) + ( - tke / dCAPE_dz_j )

            else

              if ( gr%grid_dir_indx > 0 ) then
                ! Ascending grid
                j_zm = j
              else ! gr%grid_dir_indx < 0
                ! Descending grid
                j_zm = j+1
              endif

              ! Case used for most scenarios where dCAPE/dz|_(z_1) /= dCAPE/dz|_(z_0)
              ! Find the remaining distance z - z_0 that it takes to exhaust the
              ! remaining TKE (tke_i), using the quadratic formula (only the
              ! negative (-) root works in this scenario).
              invrs_dCAPE_diff = one / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )

              Lscale_up(i,k) &
              = Lscale_up(i,k) &
                - dCAPE_dz_j_minus_1 * invrs_dCAPE_diff &
                  * gr%grid_dir * gr%dzm(i,j_zm) &
                - sqrt( dCAPE_dz_j_minus_1**2 &
                        - two * tke * gr%grid_dir * gr%invrs_dzm(i,j_zm) &
                          * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) ) &
                  * invrs_dCAPE_diff  * gr%grid_dir * gr%dzm(i,j_zm)

            endif

          end if

        else    ! TKE for parcel at level (k) was exhaused before one full grid level

          if ( gr%grid_dir_indx > 0 ) then
            ! Ascending grid
            kp1_zm = k+1
          else ! gr%grid_dir_indx < 0
            ! Descending grid
            kp1_zm = k
          endif

          ! Find the remaining distance z - z_0 that it takes to exhaust the
          ! remaining TKE (tke_i), using the quadratic formula. Simplified
          ! since dCAPE_dz_j_minus_1 = 0.0
          Lscale_up(i,k) &
          = Lscale_up(i,k) &
            - sqrt( - two * tke_i(i,k) &
                      * gr%grid_dir * gr%dzm(i,kp1_zm) * dCAPE_dz_1(i,k+gr%grid_dir_indx) ) &
              / dCAPE_dz_1(i,k+gr%grid_dir_indx)

        endif


        ! If a parcel at a previous grid level can rise past the parcel at this grid level
        ! then this one should also be able to rise up to that height. This feature insures
        ! that the profile of Lscale_up will be smooth, thus reducing numerical instability.
        if ( gr%zt(i,k) + Lscale_up(i,k) < Lscale_up_max_alt ) then

            ! A lower starting parcel can ascend higher than this one, set height to the max
            ! that any lower starting parcel can ascend to
            Lscale_up(i,k) = Lscale_up_max_alt - gr%zt(i,k)
        else

            ! This parcel can ascend higher than any below it, save final height
            Lscale_up_max_alt = Lscale_up(i,k) + gr%zt(i,k)
        end if


      end do
    end do
    !$acc end parallel loop

    ! ---------------- Downwards Length Scale Calculation ----------------

    ! Precalculate values for downward Lscale, these are useful only if a parcel can descend
    ! more than one level. They are used in the equations that calculate thl and rt
    ! recursively for a parcel as it descends
    !$acc parallel loop gang vector collapse(2) default(present)
    do j = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          jp1_zm = j+1
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          jp1_zm = j
        endif

        thl_par_j_precalc(i,j) = thlm(i,j) - thlm(i,j+gr%grid_dir_indx) * exp_mu_dzm(i,jp1_zm)  &
                               - ( thlm(i,j) - thlm(i,j+gr%grid_dir_indx) ) * entrain_coef(i,jp1_zm)

        rt_par_j_precalc(i,j) = rtm(i,j) - rtm(i,j+gr%grid_dir_indx) * exp_mu_dzm(i,jp1_zm)  &
                              - ( rtm(i,j) - rtm(i,j+gr%grid_dir_indx) ) * entrain_coef(i,jp1_zm)

      end do
    end do
    !$acc end parallel loop

    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations.

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    !$acc parallel loop gang vector collapse(2) default(present)
    do j = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          jp1_zm = j+1
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          jp1_zm = j
        endif

        thl_par_1(i,j) &
        = thlm(i,j) &
          - ( thlm(i,j) - thlm(i,j+gr%grid_dir_indx) ) * entrain_coef(i,jp1_zm)

        tl_par_1(i,j) = thl_par_1(i,j) * exner(i,j)

        rt_par_1(i,j) &
        = rtm(i,j) &
          - ( rtm(i,j) - rtm(i,j+gr%grid_dir_indx) ) * entrain_coef(i,jp1_zm)

      end do
    end do
    !$acc end parallel loop

    ! Caclculate initial rsatl for parcels at each grid level, this function is elemental

    ! The entire pressure and temperature arrays are passed as
    ! argument and the sub-arrays are choosen using
    ! start_index. This workaround is used to solve
    ! subarray issues with OpenACC.
    ! rsatl_par_1(i,2:) = sat_mixrat_liq_acc( nz-1, p_in_Pa(i,2:), tl_par_1(i,2:) )
    ! since subarray 2:, the start_index is 2 and it is an optional argument
    start_index = gr%k_lb_zt
    rsatl_par_1 = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, tl_par_1, saturation_formula, &
                                      start_index )

    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do j = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          jp1_zm = j+1
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          jp1_zm = j
        endif

        tl_par_j_sqd = tl_par_1(i,j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
        !
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )
        s_par_1(i,j) = ( rt_par_1(i,j) - rsatl_par_1(i,j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(i,j) )

        rc_par_1(i,j) = max( s_par_1(i,j), zero_threshold )

        ! theta_v of entraining parcel at grid level j
        thv_par_1(i,j) = thl_par_1(i,j) + ep1 * thv_ds(i,j) * rt_par_1(i,j) + Lv_coef(i,j) &
                         * rc_par_1(i,j)

        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(i,j) = grav_on_thvm(i,j) * ( thv_par_1(i,j) - thvm(i,j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(i,j) = one_half * dCAPE_dz_1(i,j) * gr%grid_dir * gr%dzm(i,jp1_zm)

      end do


      ! Calculate Lscale_down for each grid level. If the TKE from a parcel has not been completely
      ! exhausted by the initial change then continue the exhaustion calculations here for a single
      ! grid level at a time until the TKE is exhausted.

      Lscale_down_min_alt = gr%zt(i,gr%k_ub_zt)  ! Set initial min value for Lscale_down to max zt
      do k = gr%k_ub_zt, gr%k_lb_zt+gr%grid_dir_indx, -gr%grid_dir_indx

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i,k) - CAPE_incr_1(i,k-gr%grid_dir_indx) > zero ) then

          ! Calculate new TKE for parcel
          tke = tke_i(i,k) - CAPE_incr_1(i,k-gr%grid_dir_indx)

          ! Set j to 2 levels below current Lscale_down level, this is because we've already
          ! determined that the parcel can descend at least 1 full level
          j = k - 2*gr%grid_dir_indx

          ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
          thl_par_j = thl_par_1(i,k-gr%grid_dir_indx)
          rt_par_j = rt_par_1(i,k-gr%grid_dir_indx)
          dCAPE_dz_j_plus_1 = dCAPE_dz_1(i,k-gr%grid_dir_indx)


          ! Continue change in TKE calculations until it is exhausted or the min grid
          ! level has been reached. j is the next grid level below the level that can
          ! be reached for a parcel starting at level k. If TKE is exhausted in this loop
          ! that means the parcel starting at k cannot sink to level j, but can sink to j+1
          !
          ! Update for generalized gridding:
          !
          ! Ascending grid: do while ( j >= 1 ) -- j is decreasing
          ! Descending grid: do while ( j <= nzt ) -- j is increasing
          !
          ! Either way, j is going downwards in altitude.
          !
          ! Reconciling the grids:
          !
          ! Step 1: Use gr%k_lb_zt for the lower boundary:
          ! Ascending grid: do while ( j >= gr%k_lb_zt ) -- j is decreasing
          ! Descending grid: do while ( j <= gr%k_lb_zt ) -- j is increasing
          !
          ! Step 2: Reconciling the greater than / less than:
          ! Ascending grid: do while ( j >= 1 ) -- j is decreasing
          ! Descending grid: do while ( -j >= -nzt ) -- j is decreasing
          !
          ! Generalized Expression:
          ! do while ( gr%grid_dir_indx * j >= gr%grid_dir_indx * gr%k_lb_zt );
          ! where for an ascending grid, gr%grid_dir_indx = 1 and gr%k_lb_zt = 1;
          ! and for a descending grid, gr%grid_dir_indx = -1 and gr%k_lb_zt = nzt.
          do while ( gr%grid_dir_indx * j >= gr%grid_dir_indx * gr%k_lb_zt )

            if ( gr%grid_dir_indx > 0 ) then
              ! Ascending grid
              jp1_zm = j+1
            else ! gr%grid_dir_indx < 0
              ! Descending grid
              jp1_zm = j
            endif

            ! thl, rt of parcel are conserved except for entrainment
            !
            ! The values of thl_env and rt_env are treated as changing linearly for a parcel
            ! of air descending from level j to level j-1

            ! theta_l of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(thl_par)/dz = - Lscale_mu_entrain * ( thl_par - thl_env )
            thl_par_j = thl_par_j_precalc(i,j) + thl_par_j * exp_mu_dzm(i,jp1_zm)


            ! r_t of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(rt_par)/dz = - Lscale_mu_entrain * ( rt_par - rt_env )
            rt_par_j = rt_par_j_precalc(i,j) + rt_par_j * exp_mu_dzm(i,jp1_zm)


            ! Include effects of latent heating on Lscale_up 6/12/00
            ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
            ! Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(i,j)

            rsatl_par_j = sat_mixrat_liq_api( p_in_Pa(i,j), tl_par_j, saturation_formula )

            tl_par_j_sqd = tl_par_j**2

            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
            !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )
            ! and SD's beta (eqn. 8),
            !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
            !
            ! Simplified by multiplying top and bottom by tl^2 to avoid a
            ! divide and precalculating ep * Lv**2 / ( Rd * cp )
            s_par_j = (rt_par_j - rsatl_par_j) * tl_par_j_sqd &
                      / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

            rc_par_j = max( s_par_j, zero_threshold )

            ! theta_v of entraining parcel at grid level j
            thv_par_j = thl_par_j + ep1 * thv_ds(i,j) * rt_par_j + Lv_coef(i,j) * rc_par_j

            ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
            dCAPE_dz_j = grav_on_thvm(i,j) * ( thv_par_j - thvm(i,j) )

            ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
            ! Trapezoidal estimate between grid levels j+1 and j
            CAPE_incr = one_half * ( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) &
                        * gr%grid_dir * gr%dzm(i,jp1_zm)

            ! Exit loop early if tke has been exhaused between level j+1 and j
            if ( tke - CAPE_incr <= zero ) then
              exit
            endif

            ! Save previous dCAPE value for next cycle
            dCAPE_dz_j_plus_1 = dCAPE_dz_j

            ! Caclulate new TKE and increment j
            tke = tke - CAPE_incr
            j = j - gr%grid_dir_indx

          enddo

          ! Add full grid level thickness for each grid level that was passed without the TKE
          ! being exhausted, difference between starting level (k) and last level passed (j+1)
          Lscale_down(i,k) = Lscale_down(i,k) + gr%zt(i,k) - gr%zt(i,j+gr%grid_dir_indx)


          if ( gr%grid_dir_indx * j >= gr%grid_dir_indx * gr%k_lb_zt ) then

            ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
            ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.

            if ( abs( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) * 2 <= &
                 abs( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * eps ) then

              ! Special case where dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) = 0
              ! Find the remaining distance z_0 - z that it takes to
              ! exhaust the remaining TKE

              Lscale_down(i,k) = Lscale_down(i,k) + ( tke / dCAPE_dz_j )

            else

              if ( gr%grid_dir_indx > 0 ) then
                ! Ascending grid
                jp1_zm = j+1
              else ! gr%grid_dir_indx < 0
                ! Descending grid
                jp1_zm = j
              endif

              ! Case used for most scenarios where dCAPE/dz|_(z_(-1)) /= dCAPE/dz|_(z_0)
              ! Find the remaining distance z_0 - z that it takes to exhaust the
              ! remaining TKE (tke_i), using the quadratic formula (only the
              ! negative (-) root works in this scenario) -- however, the
              ! negative (-) root is divided by another negative (-) factor,
              ! which results in an overall plus (+) sign in front of the
              ! square root term in the equation below).
              invrs_dCAPE_diff = one / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )

              Lscale_down(i,k) &
              = Lscale_down(i,k) &
                - dCAPE_dz_j_plus_1 * invrs_dCAPE_diff &
                  * gr%grid_dir * gr%dzm(i,jp1_zm)  &
                + sqrt( dCAPE_dz_j_plus_1**2 &
                        + two * tke * gr%grid_dir * gr%invrs_dzm(i,jp1_zm)  &
                          * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) )  &
                  * invrs_dCAPE_diff * gr%grid_dir * gr%dzm(i,jp1_zm)

            endif

          end if

        else    ! TKE for parcel at level (k) was exhaused before one full grid level

          if ( gr%grid_dir_indx > 0 ) then
            ! Ascending grid
            k_zm = k
          else ! gr%grid_dir_indx < 0
            ! Descending grid
            k_zm = k+1
          endif

          ! Find the remaining distance z_0 - z that it takes to exhaust the
          ! remaining TKE (tke_i), using the quadratic formula. Simplified
          ! since dCAPE_dz_j_plus_1 = 0.0
          Lscale_down(i,k) &
          = Lscale_down(i,k) &
            + sqrt( two * tke_i(i,k) &
                    * gr%grid_dir * gr%dzm(i,k_zm) * dCAPE_dz_1(i,k-gr%grid_dir_indx) ) &
              / dCAPE_dz_1(i,k-gr%grid_dir_indx)

        endif

        ! If a parcel at a previous grid level can descend past the parcel at this grid level
        ! then this one should also be able to descend down to that height. This feature insures
        ! that the profile of Lscale_down will be smooth, thus reducing numerical instability.
        if ( gr%zt(i,k) - Lscale_down(i,k) > Lscale_down_min_alt ) then
          Lscale_down(i,k) = gr%zt(i,k) - Lscale_down_min_alt
        else
          Lscale_down_min_alt = gr%zt(i,k) - Lscale_down(i,k)
        end if

      end do
    end do
    !$acc end parallel loop

      ! ---------------- Final Lscale Calculation ----------------

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do k = 1, nzt, 1

        ! Make lminh a linear function starting at value lmin at the bottom
        ! and going to zero at 500 meters in altitude.
        if( l_implemented ) then

          ! Within a host model, increase mixing length in 500 m layer above *ground*
          lminh = max( zero_threshold, &
                       Lscale_sfclyr_depth - ( gr%zt(i,k) - gr%zm(i,gr%k_lb_zm) ) ) &
                  * lmin * invrs_Lscale_sfclyr_depth
        else

          ! In standalone mode, increase mixing length in 500 m layer above *mean sea level*
          lminh = max( zero_threshold, Lscale_sfclyr_depth - gr%zt(i,k) ) &
                  * lmin * invrs_Lscale_sfclyr_depth
        end if

        Lscale_up(i,k)    = max( lminh, Lscale_up(i,k) )
        Lscale_down(i,k)  = max( lminh, Lscale_down(i,k) )

        ! When L is large, turbulence is weakly damped
        ! When L is small, turbulence is strongly damped
        ! Use a geometric mean to determine final Lscale so that L tends to become small
        ! if either Lscale_up or Lscale_down becomes small.
        Lscale(i,k) = sqrt( Lscale_up(i,k) * Lscale_down(i,k) )

      enddo

      ! Set the value of Lscale at the upper and lower boundaries.
      Lscale(i,gr%k_ub_zt) = Lscale(i,gr%k_ub_zt-gr%grid_dir_indx)

      ! Vince Larson limited Lscale to allow host
      ! model to take over deep convection.  13 Feb 2008.
      Lscale(i,:) = min( Lscale(i,:), Lscale_max(i) )

    end do
    !$acc end parallel loop

    ! Ensure that no Lscale values are NaN
    if ( clubb_at_least_debug_level_api( 1 ) ) then

      !$acc update host( Lscale, Lscale_up, Lscale_down, &
      !$acc              thvm, thlm, rtm, em, exner, p_in_Pa, thv_ds )

      do i = 1, ngrdcol
        call length_check( nzt, Lscale(i,:), Lscale_up(i,:), Lscale_down(i,:), & ! intent(in)
                           err_info ) ! intent(inout)

        if ( err_info%err_code(i) == clubb_fatal_error ) then

          write(fstderr,*) err_info%err_header(i)
          write(fstderr,*) "Errors in compute_mixing_length subroutine in grid column ", i

          write(fstderr,*) "Intent(in)"

          write(fstderr,*) "thvm = ", thvm(i,:)
          write(fstderr,*) "thlm = ", thlm(i,:)
          write(fstderr,*) "rtm = ", rtm(i,:)
          write(fstderr,*) "em = ", em(i,:)
          write(fstderr,*) "exner = ", exner(i,:)
          write(fstderr,*) "p_in_Pa = ", p_in_Pa(i,:)
          write(fstderr,*) "thv_ds = ", thv_ds(i,:)

          write(fstderr,*) "Intent(out)"

          write(fstderr,*) "Lscale = ", Lscale(i,:)
          write(fstderr,*) "Lscale_up = ", Lscale_up(i,:)
          write(fstderr,*) "Lscale_down = ", Lscale_down(i,:)

        endif ! Fatal error
      end do

    end if

    !$acc exit data delete( exp_mu_dzm, invrs_dzm_on_mu, grav_on_thvm, Lv_coef, &
    !$acc                   entrain_coef, thl_par_j_precalc, rt_par_j_precalc, &
    !$acc                   tl_par_1, rt_par_1, rsatl_par_1, thl_par_1, dCAPE_dz_1, &
    !$acc                   s_par_1, rc_par_1, CAPE_incr_1, thv_par_1, tke_i )

    return

  end subroutine compute_mixing_length

  !=============================================================================

  subroutine calc_Lscale_directly ( ngrdcol, nzm, nzt, gr, &
                                    l_implemented, p_in_Pa, exner, rtm,    &
                                    thlm, thvm, rtp2_zt, thlp2_zt, rtpthlp_zt, &
                                    pdf_params, em, thv_ds_zt, Lscale_max, lmin, &
                                    clubb_params, &
                                    saturation_formula, &
                                    l_Lscale_plume_centered, &
                                    Lscale_method, &
                                    stats, err_info,         &
                                    Lscale, Lscale_up, Lscale_down)

    ! Description:
    !   Select the legacy or coherent parcel calculation (methods 1 or 3)
    !   using environmental profiles, momentum-level TKE, PDF properties, and
    !   per-column parameters. Return thermodynamic-level combined, upward,
    !   and downward lengths [m], updating stats and err_info.
    !   If source constant l_avg_Lscale is enabled, compute perturbed launches
    !   before the mean launch so the directional outputs describe the mean.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        thl_tol,      &
        rt_tol,       &
        one_half,     &
        one_third,    &
        one,          &
        unused_var

    use parameter_indices, only: &
        nparams, &
        iLscale_mu_entrain, &
        iLscale_mu_coef, &
        iLscale_pert_coef, &
        iLscale_tke_multiplier, &
        iLscale_pdf_launch_frac, &
        iLscale_up_weight

    use model_flags, only: &
        Lscale_direct, &
        Lscale_coherence

    use grid_class, only: &
        grid ! Type

    use clubb_precision, only: &
        core_rknd

    use pdf_parameter_module, only: &
        pdf_parameter

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error              ! Constant

    use constants_clubb, only:  &
        fstderr  ! Variable(s)

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    !--------------------------------- Input Variables ---------------------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(in) :: &
      gr

    logical, intent(in) ::  &
      l_implemented ! True if CLUBB is being implemented and run in a host model

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      rtp2_zt,    & ! Total-water variance [(kg/kg)^2]
      thlp2_zt,   & ! Liquid-water potential-temperature variance [K^2]
      rtpthlp_zt, & ! Total-water/temperature covariance [K kg/kg]
      thlm,       & ! Mean liquid-water potential temperature [K]
      thvm,       & ! Mean virtual potential temperature [K]
      rtm,        & ! Mean total-water mixing ratio [kg/kg]
      p_in_Pa,    & ! Air pressure (thermodynamic levels)       [Pa]
      exner,      & ! Exner function [-]
      thv_ds_zt      ! Reference-state virtual potential temperature [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      em   ! Turbulent kinetic energy [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      Lscale_max ! Per-column cap on combined mixing length [m]

    real( kind = core_rknd ), intent(in) ::  &
      lmin ! Surface value of the height-dependent length floor [m]

    type (pdf_parameter), intent(in) :: &
      pdf_params    ! PDF Parameters  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      saturation_formula, & ! Integer that stores the saturation formula to be used
      Lscale_method         ! Mixing-length formulation

    logical, intent(in) :: &
      l_Lscale_plume_centered ! Alternate that uses the PDF to compute the perturbed values

    !----------------------- Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: &
      stats ! Statistics accumulator

    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    !--------------------------------- Output Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    !--------------------------------- Local Variables ---------------------------------
    integer :: k, i ! Vertical-level and column loop indices

    logical, parameter :: &
      l_avg_Lscale = .false.  ! Lscale is calculated in subroutine compute_mixing_length
                              ! if l_avg_Lscale is true, compute_mixing_length is called
                              ! two additional times with
                              ! perturbed values of rtm and thlm.  An average value of Lscale
                              ! from the three calls to compute_mixing_length is then calculated.
                              ! This reduces temporal noise in RICO, BOMEX, LBA, and other cases.

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      sign_rtpthlp_zt,              & ! Sign of the covariance rtpthlp       [-]
      Lscale_pert_1, Lscale_pert_2, & ! For avg. calculation of Lscale  [m]
      thlm_pert_1, thlm_pert_2,     & ! For avg. calculation of Lscale  [K]
      rtm_pert_1, rtm_pert_2,       & ! For avg. calculation of Lscale  [kg/kg]
      thlm_pert_pos_rt, thlm_pert_neg_rt, & ! For avg. calculation of Lscale [K]
      rtm_pert_pos_rt, rtm_pert_neg_rt      ! For avg. calculation of Lscale [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      mu_pert_1, mu_pert_2, &
      mu_pert_pos_rt, mu_pert_neg_rt  ! For l_Lscale_plume_centered

    real( kind = core_rknd ) :: &
      Lscale_pert_coef

    !Lscale_weight Uncomment this if you need to use this vairable at some
    !point.

    !--------------------------------- Begin Code ---------------------------------

    if ( l_avg_Lscale .and. Lscale_method /= Lscale_direct ) then
      write(fstderr,*) err_info%err_header_global
      write(fstderr,*) "l_avg_Lscale requires Lscale_method=1"
      write(fstderr,*) "Fatal error in calc_Lscale_directly"
      err_info%err_code = clubb_fatal_error
      return
    end if

    !$acc enter data create( sign_rtpthlp_zt, Lscale_pert_1, Lscale_pert_2, &
    !$acc                    thlm_pert_1, thlm_pert_2, rtm_pert_1, rtm_pert_2, &
    !$acc                    thlm_pert_pos_rt, thlm_pert_neg_rt, rtm_pert_pos_rt, &
    !$acc                    rtm_pert_neg_rt, &
    !$acc                    mu_pert_1, mu_pert_2, mu_pert_pos_rt, mu_pert_neg_rt )

    if ( clubb_at_least_debug_level_api( 0 ) ) then

      if ( l_Lscale_plume_centered .and. .not. l_avg_Lscale ) then
        write(fstderr,*) err_info%err_header_global
        write(fstderr,*) "l_Lscale_plume_centered requires l_avg_Lscale"
        write(fstderr,*) "Fatal error in calc_Lscale_directly"
        ! General error -> set all entries to clubb_fatal_error
        err_info%err_code = clubb_fatal_error
        return
      end if

    end if

    if ( l_avg_Lscale .and. .not. l_Lscale_plume_centered ) then

      ! Call compute length two additional times with perturbed values
      ! of rtm and thlm so that an average value of Lscale may be calculated.

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt, 1
        do  i = 1, ngrdcol
          sign_rtpthlp_zt(i,k) = sign( one, rtpthlp_zt(i,k) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt, 1
        do  i = 1, ngrdcol
          rtm_pert_1(i,k)  = rtm(i,k) + clubb_params(i,iLscale_pert_coef) &
                                        * sqrt( max( rtp2_zt(i,k), rt_tol**2 ) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt, 1
        do  i = 1, ngrdcol
          thlm_pert_1(i,k) = thlm(i,k) + sign_rtpthlp_zt(i,k) * clubb_params(i,iLscale_pert_coef) &
                                         * sqrt( max( thlp2_zt(i,k), thl_tol**2 ) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do  i = 1, ngrdcol
        mu_pert_1(i) = clubb_params(i,iLscale_mu_entrain) / clubb_params(i,iLscale_mu_coef)
      end do
      !$acc end parallel loop

      call compute_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm_pert_1, & ! In
                                  rtm_pert_1, em, Lscale_max, p_in_Pa,      & ! In
                                  exner, thv_ds_zt, mu_pert_1, lmin,        & ! In
                                  saturation_formula, l_implemented,        & ! In
                                  err_info,                                 & ! InOut
                                  Lscale_pert_1, Lscale_up, Lscale_down )     ! Out

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt, 1
        do  i = 1, ngrdcol
          rtm_pert_2(i,k)  = rtm(i,k) - clubb_params(i,iLscale_pert_coef) &
                                        * sqrt( max( rtp2_zt(i,k), rt_tol**2 ) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt, 1
        do  i = 1, ngrdcol
          thlm_pert_2(i,k) = thlm(i,k) - sign_rtpthlp_zt(i,k) * clubb_params(i,iLscale_pert_coef) &
                               * sqrt( max( thlp2_zt(i,k), thl_tol**2 ) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do  i = 1, ngrdcol
        mu_pert_2(i) = clubb_params(i,iLscale_mu_entrain) * clubb_params(i,iLscale_mu_coef)
      end do
      !$acc end parallel loop

      call compute_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm_pert_2, & ! In
                                  rtm_pert_2, em, Lscale_max, p_in_Pa,      & ! In
                                  exner, thv_ds_zt, mu_pert_2, lmin,        & ! In
                                  saturation_formula, l_implemented,        & ! In
                                  err_info,                                 & ! InOut
                                  Lscale_pert_2, Lscale_up, Lscale_down )     ! Out

    else if ( l_avg_Lscale .and. l_Lscale_plume_centered ) then

      ! Take the values of thl and rt based one 1st or 2nd plume

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          sign_rtpthlp_zt(i,k) = sign( one, rtpthlp_zt(i,k) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol

          Lscale_pert_coef = clubb_params(i,iLscale_pert_coef)

          if ( pdf_params%rt_1(i,k) > pdf_params%rt_2(i,k) ) then

            rtm_pert_pos_rt(i,k) = pdf_params%rt_1(i,k) &
                       + Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_1(i,k), rt_tol**2 ) )

            thlm_pert_pos_rt(i,k) = pdf_params%thl_1(i,k) + ( sign_rtpthlp_zt(i,k) &
                     * Lscale_pert_coef * sqrt( max( pdf_params%varnce_thl_1(i,k), thl_tol**2 ) ) )

            thlm_pert_neg_rt(i,k) = pdf_params%thl_2(i,k) - ( sign_rtpthlp_zt(i,k) &
                     * Lscale_pert_coef * sqrt( max( pdf_params%varnce_thl_2(i,k), thl_tol**2 ) ) )

            rtm_pert_neg_rt(i,k) = pdf_params%rt_2(i,k) &
                       - Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_2(i,k), rt_tol**2 ) )

            !Lscale_weight = pdf_params%mixt_frac(i,k)

          else

            rtm_pert_pos_rt(i,k) = pdf_params%rt_2(i,k) &
                       + Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_2(i,k), rt_tol**2 ) )

            thlm_pert_pos_rt(i,k) = pdf_params%thl_2(i,k) + ( sign_rtpthlp_zt(i,k) &
                     * Lscale_pert_coef * sqrt( max( pdf_params%varnce_thl_2(i,k), thl_tol**2 ) ) )

            thlm_pert_neg_rt(i,k) = pdf_params%thl_1(i,k) - ( sign_rtpthlp_zt(i,k) &
                     * Lscale_pert_coef * sqrt( max( pdf_params%varnce_thl_1(i,k), thl_tol**2 ) ) )

            rtm_pert_neg_rt(i,k) = pdf_params%rt_1(i,k) &
                       - Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_1(i,k), rt_tol**2 ) )

            !Lscale_weight = 1.0_core_rknd - pdf_params%mixt_frac(i,k)

          end if

        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        mu_pert_pos_rt(i) = clubb_params(i,iLscale_mu_entrain) / clubb_params(i,iLscale_mu_coef)
        mu_pert_neg_rt(i) = clubb_params(i,iLscale_mu_entrain) * clubb_params(i,iLscale_mu_coef)
      end do
      !$acc end parallel loop

      ! Call length with perturbed values of thl and rt
      call compute_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm_pert_pos_rt, & ! In
                                  rtm_pert_pos_rt, em, Lscale_max, p_in_Pa,      & ! In
                                  exner, thv_ds_zt, mu_pert_pos_rt, lmin,        & ! In
                                  saturation_formula, l_implemented,             & ! In
                                  err_info,                                      & ! InOut
                                  Lscale_pert_1, Lscale_up, Lscale_down )          ! Out

      call compute_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm_pert_neg_rt, & ! In
                                  rtm_pert_neg_rt, em, Lscale_max, p_in_Pa,      & ! In
                                  exner, thv_ds_zt, mu_pert_neg_rt, lmin,        & ! In
                                  saturation_formula, l_implemented,             & ! In
                                  err_info,                                      & ! InOut
                                  Lscale_pert_2, Lscale_up, Lscale_down )          ! Out

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          Lscale_pert_1(i,k) = unused_var ! Undefined
          Lscale_pert_2(i,k) = unused_var ! Undefined
        end do
      end do
      !$acc end parallel loop

    end if ! l_avg_Lscale

    if ( stats%l_sample ) then
      !$acc update host( Lscale_pert_1, Lscale_pert_2 )
      call stats_update( "Lscale_pert_1", Lscale_pert_1, stats )
      call stats_update( "Lscale_pert_2", Lscale_pert_2, stats )
    end if


    ! ********** NOTE: **********
    ! This call to compute_mixing_length must be last.  Otherwise, the values
    ! of
    ! Lscale_up and Lscale_down in stats will be based on perturbation length
    ! scales
    ! rather than the mean length scale.

    ! Diagnose CLUBB's turbulent mixing length scale.
    select case ( Lscale_method )

    case ( Lscale_coherence )

      call compute_coherent_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm, & ! In
                                            rtm, pdf_params, em, Lscale_max, p_in_Pa, & ! In
                                            exner, thv_ds_zt, clubb_params(:,iLscale_mu_entrain), & ! In
                                            clubb_params(:,iLscale_tke_multiplier), & ! In
                                            clubb_params(:,iLscale_pdf_launch_frac), & ! In
                                            clubb_params(:,iLscale_up_weight), & ! In
                                            lmin, l_implemented, & ! In
                                            err_info, & ! InOut
                                            Lscale, Lscale_up, Lscale_down ) ! Out

    case ( Lscale_direct )

      call compute_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm, & ! In
                                  rtm, em, Lscale_max, p_in_Pa,      & ! In
                                  exner, thv_ds_zt, clubb_params(:,iLscale_mu_entrain), lmin, & ! In
                                  saturation_formula, l_implemented, & ! In
                                  err_info,                          & ! InOut
                                  Lscale, Lscale_up, Lscale_down )     ! Out

    case default

      write(fstderr,*) err_info%err_header_global
      write(fstderr,*) "calc_Lscale_directly cannot diagnose Lscale_method=", Lscale_method
      err_info%err_code = clubb_fatal_error
      return

    end select

    if ( l_avg_Lscale ) then
      if ( l_Lscale_plume_centered ) then
        ! Weighted average of mean, pert_1, & pert_2
        !       Lscale = 0.5_core_rknd * ( Lscale + Lscale_weight*Lscale_pert_1 &
        !                                  + (1.0_core_rknd-Lscale_weight)*Lscale_pert_2
        !                                  )
        ! Weighted average of just the perturbed values
        !       Lscale = Lscale_weight*Lscale_pert_1 +
        !       (1.0_core_rknd-Lscale_weight)*Lscale_pert_2

        ! Un-weighted average of just the perturbed values
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            Lscale(i,k) = one_half *( Lscale_pert_1(i,k) + Lscale_pert_2(i,k) )
          end do
        end do
        !$acc end parallel loop
      else
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            Lscale(i,k) = one_third * ( Lscale(i,k) + Lscale_pert_1(i,k) + Lscale_pert_2(i,k) )
          end do
        end do
        !$acc end parallel loop
      end if
    end if

    !$acc exit data delete( sign_rtpthlp_zt, Lscale_pert_1, Lscale_pert_2, &
    !$acc                   thlm_pert_1, thlm_pert_2, rtm_pert_1, rtm_pert_2, &
    !$acc                   thlm_pert_pos_rt, thlm_pert_neg_rt, rtm_pert_pos_rt, &
    !$acc                   rtm_pert_neg_rt, &
    !$acc                   mu_pert_1, mu_pert_2, mu_pert_pos_rt, mu_pert_neg_rt )

   return

 end subroutine calc_Lscale_directly

  !===============================================================================

 subroutine diagnose_Lscale_from_tau( nzm, nzt, ngrdcol, gr, & ! intent in
                        upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, & !intent in
                        ice_supersat_frac, &! intent in
                        em, & ! intent in
                        ufmin, tau_const, & ! intent in
                        sfc_elevation, Lscale_max, & ! intent in
                        clubb_params, & ! intent in
                        stats,         & ! intent in
                        l_e3sm_config, & ! intent in
                        l_smooth_Heaviside_tau_wpxp, & ! intent in
                        brunt_vaisala_freq_sqd_smth, Ri_zm, & ! intent in
                        err_info, & ! intent inout
                        invrs_tau_zt, invrs_tau_zm, & ! intent out
                        invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd, & ! intent out
                        invrs_tau_shear, invrs_tau_N2_iso, & ! intent out
                        invrs_tau_wp2_zm, invrs_tau_xp2_zm, & ! intent out
                        invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm, & ! intent out
                        tau_max_zm, tau_max_zt, tau_zm, tau_zt, & !intent out
                        Lscale, Lscale_up, Lscale_down)! intent out

! Description:
!     Diagnose inverse damping times and mixing length from surface stress,
!     momentum-level TKE, shear, and stability profiles (method 2). Return the
!     separate inverse-time contributions [1/s], bounded time scales [s], and
!     thermodynamic-level lengths [m]; update statistics and error information.
! References:
!     Guo et al.(2021, JAMES)
!--------------------------------------------------------------------------------------------------

    use advance_helper_module, only: &
        smooth_heaviside_peskin, &
        smooth_min, smooth_max

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    use constants_clubb, only: &
        one_fourth,         &
        one_half,           &
        vonk,               &
        zero,               &
        one,                &
        two,                &
        em_min,             &
        zero_threshold,     &
        eps,                &
        min_max_smth_mag,   &
        fstderr

    use grid_class, only: &
        grid, & ! Type
        zt2zm_api, &
        zm2zt_api, &
        zm2zt2zm, &
        zm2zt2zm

    use clubb_precision, only: &
        core_rknd

    use parameter_indices, only: &
        nparams,                     & ! Variable(s)
        iC_invrs_tau_bkgnd,          &
        iC_invrs_tau_shear,          &
        iC_invrs_tau_sfc,            &
        iC_invrs_tau_N2,             &
        iC_invrs_tau_N2_wp2 ,        &
        iC_invrs_tau_N2_wpxp,        &
        iC_invrs_tau_N2_xp2,         &
        iC_invrs_tau_wpxp_N2_thresh, &
        iC_invrs_tau_N2_clear_wp3,   &
        iC_invrs_tau_wpxp_Ri,        &
        ialtitude_threshold,         &
        iwpxp_Ri_exp,                &
        iz_displace

    use error_code, only: &
      clubb_fatal_error, &
      clubb_at_least_debug_level_api

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    !--------------------------------- Input Variables ---------------------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      upwp_sfc,      &
      vpwp_sfc

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      ice_supersat_frac

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      em,                           &
      ddzt_umvm_sqd,                &
      brunt_vaisala_freq_sqd_smth,  &     ! smoothed Buoyancy frequency squared, N^2     [s^-2]
      Ri_zm                               ! Richardson number

    real(kind = core_rknd), intent(in) :: &
      ufmin,         &
      tau_const

    real(kind = core_rknd), dimension(ngrdcol), intent(in) :: &
      sfc_elevation, &
      Lscale_max

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(stats_type), intent(inout) :: &
      stats

    logical, intent(in) :: &
      l_e3sm_config,              &
      l_smooth_Heaviside_tau_wpxp   ! Use the smoothed Heaviside 'Peskin' function
                                    ! to compute invrs_tau_wpxp_zm

    !--------------------------- Input/Output Variables ---------------------------
    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    !--------------------------------- Output Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      invrs_tau_zm,                 &
      invrs_tau_sfc,                &
      invrs_tau_no_N2_zm,           &
      invrs_tau_bkgnd,              &
      invrs_tau_shear,              &
      invrs_tau_N2_iso,             &
      invrs_tau_wp2_zm,             &
      invrs_tau_xp2_zm,             &
      invrs_tau_wp3_zm,             &
      invrs_tau_wpxp_zm,            &
      tau_max_zm,                   &
      tau_zm

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      invrs_tau_zt,                 &
      invrs_tau_wp3_zt,             &
      tau_max_zt,                   &
      tau_zt,                       &
      Lscale,                       &
      Lscale_up,                    &
      Lscale_down

    !--------------------------------- Local Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      brunt_freq_pos,               &
      brunt_freq_out_cloud,         &
      bvf_thresh,                   & ! temporatory array
      H_invrs_tau_wpxp_N2             ! Heaviside function for clippings of invrs_tau_wpxp_N2

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      ustar

    real( kind = core_rknd ) :: &
      C_invrs_tau_N2,             &
      C_invrs_tau_N2_wp2

    real( kind = core_rknd ), parameter :: &
      heaviside_smth_range = 1.0e-0_core_rknd ! range where Heaviside function is smoothed

    logical, parameter :: l_smooth_min_max = .false. ! whether to apply smooth min and max function

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      norm_ddzt_umvm, &
      smooth_norm_ddzt_umvm, &
      brunt_vaisala_freq_clipped, &
      ice_supersat_frac_zm, &
      invrs_tau_shear_smooth, &
      Ri_zm_smooth, &
      em_clipped, &
      tau_zm_unclipped, &
      tmp_calc, &
      tmp_calc_max, &
      tmp_calc_min_max

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      em_zt, &
      tau_zt_unclipped

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      tmp_calc_ngrdcol

    integer :: i, k

    !--------------------------------- Begin Code ---------------------------------

    !$acc enter data create( brunt_freq_pos, brunt_freq_out_cloud, &
    !$acc                    bvf_thresh, H_invrs_tau_wpxp_N2, ustar, &
    !$acc                    norm_ddzt_umvm, smooth_norm_ddzt_umvm, &
    !$acc                    brunt_vaisala_freq_clipped, &
    !$acc                    ice_supersat_frac_zm, invrs_tau_shear_smooth, &
    !$acc                    em_zt, tmp_calc_ngrdcol )

    !$acc enter data if( l_smooth_min_max ) &
    !$acc            create( tau_zm_unclipped, tau_zt_unclipped, Ri_zm_smooth, em_clipped, &
    !$acc                    tmp_calc, tmp_calc_max, tmp_calc_min_max )

    em_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, em(:,:), em_min )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      if ( gr%zm(i,gr%k_lb_zm) - sfc_elevation(i) + clubb_params(i,iz_displace) < eps ) then
        ! Error in grid column i -> set ith entry to clubb_fatal_error
        err_info%err_code(i) = clubb_fatal_error
        write(fstderr, *) err_info%err_header(i)
      end if
    end do
    !$acc end parallel loop

    if ( clubb_at_least_debug_level_api(0) ) then
      !$acc update host( err_info%err_code )
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) "Lowest zm grid level is below ground in CLUBB."
        return
      end if
    end if

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      tmp_calc_ngrdcol(i) = ( upwp_sfc(i)**2 + vpwp_sfc(i)**2 )**one_fourth
    end do
    !$acc end parallel loop

    if ( l_smooth_min_max ) then

      ! tmp_calc_ngrdcol used here because openacc has an issue with
      !  a variable being both input and output of a function
      ustar = smooth_max( ngrdcol, tmp_calc_ngrdcol, &
                          ufmin, &
                          min_max_smth_mag )

    else

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        ustar(i) = max( tmp_calc_ngrdcol(i), ufmin )
      end do
      !$acc end parallel loop

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_bkgnd(i,k) = clubb_params(i,iC_invrs_tau_bkgnd) / tau_const
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        norm_ddzt_umvm(i,k) = sqrt( ddzt_umvm_sqd(i,k) )
      end do
    end do
    !$acc end parallel loop

    smooth_norm_ddzt_umvm = zm2zt2zm( nzm, nzt, ngrdcol, gr, norm_ddzt_umvm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_shear_smooth(i,k) = clubb_params(i,iC_invrs_tau_shear) &
                                      * smooth_norm_ddzt_umvm(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Enforce that invrs_tau_shear is positive
    invrs_tau_shear = smooth_max( nzm, ngrdcol, invrs_tau_shear_smooth, &
                                  zero_threshold, min_max_smth_mag )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_sfc(i,k) = clubb_params(i,iC_invrs_tau_sfc) &
                             * ( ustar(i) / vonk ) &
                               / ( gr%zm(i,k) - sfc_elevation(i) + clubb_params(i,iz_displace) )
         !C_invrs_tau_sfc * ( wp2 / vonk /ustar ) / ( gr%zm(1,:) -sfc_elevation + z_displace )
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_no_N2_zm(i,k) = invrs_tau_bkgnd(i,k) + invrs_tau_sfc(i,k) + invrs_tau_shear(i,k)
      end do
    end do
    !$acc end parallel loop

    if ( l_smooth_min_max ) then

      brunt_vaisala_freq_clipped = smooth_max( nzm, ngrdcol, zero_threshold, &
                                               brunt_vaisala_freq_sqd_smth, &
                                               1.0e-4_core_rknd * min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          brunt_freq_pos(i,k) = sqrt( brunt_vaisala_freq_clipped(i,k) )
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          brunt_freq_pos(i,k) = sqrt( max( zero_threshold, brunt_vaisala_freq_sqd_smth(i,k) ) )
        end do
      end do
      !$acc end parallel loop

    end if

    ice_supersat_frac_zm = zt2zm_api( nzm, nzt, ngrdcol, gr, ice_supersat_frac, zero_threshold )

    if ( l_smooth_min_max ) then

      ! roll this back as well once checks have passed
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tmp_calc(i,k) = one - ice_supersat_frac_zm(i,k) / 0.001_core_rknd
        end do
      end do
      !$acc end parallel loop

      tmp_calc_max = smooth_max( nzm, ngrdcol, zero_threshold, tmp_calc, &
                                 min_max_smth_mag)

      tmp_calc_min_max = smooth_min( nzm, ngrdcol, one, &
                                     tmp_calc_max, min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          brunt_freq_out_cloud(i,k) =  brunt_freq_pos(i,k) * tmp_calc_min_max(i,k)
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          brunt_freq_out_cloud(i,k) &
            = brunt_freq_pos(i,k) &
                * min(one, max(zero_threshold, &
                               one - ( ( ice_supersat_frac_zm(i,k) / 0.001_core_rknd) )))
        end do
      end do
      !$acc end parallel loop

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        if ( gr%zm(i,k) < clubb_params(i,ialtitude_threshold) ) then
          brunt_freq_out_cloud(i,k) = zero
        end if
      end do
    end do
    !$acc end parallel loop

    ! Write both bv extra terms to invrs_taus to disk
    if ( stats%l_sample ) then
      !$acc update host( brunt_freq_pos, brunt_freq_out_cloud )
      call stats_update( "bv_freq_pos", brunt_freq_pos, stats )
      call stats_update( "bv_freq_out_cloud", brunt_freq_out_cloud, stats )
    end if

    ! This time scale is used optionally for the return-to-isotropy term. It
    ! omits invrs_tau_sfc based on the rationale that the isotropization
    ! rate shouldn't be enhanced near the ground.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol

        C_invrs_tau_N2     = clubb_params(i,iC_invrs_tau_N2)
        C_invrs_tau_N2_wp2 = clubb_params(i,iC_invrs_tau_N2_wp2)

        invrs_tau_N2_iso(i,k) = invrs_tau_bkgnd(i,k) + invrs_tau_shear(i,k) &
                                + C_invrs_tau_N2_wp2 * brunt_freq_pos(i,k)

        invrs_tau_wp2_zm(i,k) = invrs_tau_no_N2_zm(i,k) + &
                                C_invrs_tau_N2 * brunt_freq_pos(i,k) + &
                                C_invrs_tau_N2_wp2 * brunt_freq_out_cloud(i,k)

        invrs_tau_zm(i,k) = invrs_tau_no_N2_zm(i,k) + C_invrs_tau_N2 * brunt_freq_pos(i,k)

      end do
    end do
    !$acc end parallel loop


    if ( l_e3sm_config ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_zm(i,k) = one_half * invrs_tau_zm(i,k)
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_xp2_zm(i,k) = invrs_tau_no_N2_zm(i,k) &
                                  + clubb_params(i,iC_invrs_tau_N2_xp2) * brunt_freq_pos(i,k) & ! 0
                                  + clubb_params(i,iC_invrs_tau_sfc) * two &
                                  * sqrt(em(i,k)) &
                                  / ( gr%zm(i,k) - sfc_elevation(i) &
                                      + clubb_params(i,iz_displace) ) ! small
        end do
      end do
      !$acc end parallel loop

      if ( l_smooth_min_max ) then

        brunt_vaisala_freq_clipped = smooth_max( nzm, ngrdcol, 1.0e-7_core_rknd, &
                                                 brunt_vaisala_freq_sqd_smth, &
                                                 1.0e-4_core_rknd * min_max_smth_mag )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            tmp_calc(i,k) = sqrt( ddzt_umvm_sqd(i,k) / brunt_vaisala_freq_clipped(i,k) )
          end do
        end do
        !$acc end parallel loop

        tmp_calc_max = smooth_max( nzm, ngrdcol, tmp_calc, &
                                   0.3_core_rknd, 0.3_core_rknd * min_max_smth_mag )

        tmp_calc_min_max = smooth_min( nzm, ngrdcol, tmp_calc_max, &
                                       one, min_max_smth_mag )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            invrs_tau_xp2_zm(i,k) =  tmp_calc_min_max(i,k) * invrs_tau_xp2_zm(i,k)
          end do
        end do
        !$acc end parallel loop

      else

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            invrs_tau_xp2_zm(i,k) &
              = min( max( sqrt( ddzt_umvm_sqd(i,k) &
                                / max( 1.0e-7_core_rknd, brunt_vaisala_freq_sqd_smth(i,k) ) ), &
                          0.3_core_rknd ), one ) &
                     * invrs_tau_xp2_zm(i,k)
          end do
        end do
        !$acc end parallel loop

      end if

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_wpxp_zm(i,k) = two * invrs_tau_zm(i,k) &
                                   + clubb_params(i,iC_invrs_tau_N2_wpxp) &
                                   * brunt_freq_out_cloud(i,k)
        end do
      end do
      !$acc end parallel loop

    else ! l_e3sm_config = false

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_xp2_zm(i,k) = invrs_tau_no_N2_zm(i,k) + &
                                  clubb_params(i,iC_invrs_tau_N2) * brunt_freq_pos(i,k) + &
                                  clubb_params(i,iC_invrs_tau_N2_xp2) * brunt_freq_out_cloud(i,k)
        end do
      end do
      !$acc end parallel loop

      ice_supersat_frac_zm = zt2zm_api( nzm, nzt, ngrdcol, gr, ice_supersat_frac, zero_threshold )

!      !$acc parallel loop gang vector collapse(2) default(present)
!      do k = 1, nzm
!        do i = 1, ngrdcol
!          if ( ice_supersat_frac_zm(i,k) <= 0.01_core_rknd &
!               .and. invrs_tau_xp2_zm(i,k)  >= 0.003_core_rknd ) then
!            invrs_tau_xp2_zm(i,k) = 0.003_core_rknd
!          end if
!        end do
!      end do
!      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_wpxp_zm(i,k) = invrs_tau_no_N2_zm(i,k) + &
                                   clubb_params(i,iC_invrs_tau_N2) * brunt_freq_pos(i,k) + &
                                   clubb_params(i,iC_invrs_tau_N2_wpxp) * brunt_freq_out_cloud(i,k)
        end do
      end do
      !$acc end parallel loop

    end if ! l_e3sm_config

    if ( l_smooth_Heaviside_tau_wpxp ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          bvf_thresh(i,k) = brunt_vaisala_freq_sqd_smth(i,k) &
                            / clubb_params(i,iC_invrs_tau_wpxp_N2_thresh) - one
          end do
      end do
      !$acc end parallel loop

      H_invrs_tau_wpxp_N2 = smooth_heaviside_peskin( nzm, ngrdcol, bvf_thresh, &
                                                     heaviside_smth_range )

    else ! l_smooth_Heaviside_tau_wpxp = .false.

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          if ( brunt_vaisala_freq_sqd_smth(i,k) > clubb_params(i,iC_invrs_tau_wpxp_N2_thresh) ) then
            H_invrs_tau_wpxp_N2(i,k) = one
          else
            H_invrs_tau_wpxp_N2(i,k) = zero
          end if
        end do
      end do
      !$acc end parallel loop

    end if ! l_smooth_Heaviside_tau_wpxp

    if ( l_smooth_min_max ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tmp_calc(i,k) = clubb_params(i,iC_invrs_tau_wpxp_Ri) * &
                          Ri_zm(i,k)**clubb_params(i,iwpxp_Ri_exp)
        end do
      end do
      !$acc end parallel loop

      Ri_zm_smooth = smooth_min( nzm, ngrdcol, tmp_calc, &
                                 12.0_core_rknd, 12.0_core_rknd * min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol

          if ( gr%zm(i,k) > clubb_params(i,ialtitude_threshold) ) then
             invrs_tau_wpxp_zm(i,k) = invrs_tau_wpxp_zm(i,k) &
                                      * ( one + H_invrs_tau_wpxp_N2(i,k) &
                                          * Ri_zm_smooth(i,k) )

          end if
        end do
      end do
      !$acc end parallel loop

    else ! l_smooth_min_max

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          if ( gr%zm(i,k) > clubb_params(i,ialtitude_threshold) ) then
             invrs_tau_wpxp_zm(i,k) = invrs_tau_wpxp_zm(i,k) &
                         * ( one  + H_invrs_tau_wpxp_N2(i,k) &
                         * min( clubb_params(i,iC_invrs_tau_wpxp_Ri) &
                         * max( Ri_zm(i,k), zero)**clubb_params(i,iwpxp_Ri_exp), 12.0_core_rknd ) )
          end if
        end do
      end do
      !$acc end parallel loop

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_wp3_zm(i,k) = invrs_tau_wp2_zm(i,k) &
                                + clubb_params(i,iC_invrs_tau_N2_clear_wp3) &
                                * brunt_freq_out_cloud(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Calculate the maximum allowable value of time-scale tau,
    ! which depends of the value of Lscale_max.
    if ( l_smooth_min_max ) then

      em_clipped = smooth_max( nzm, ngrdcol, em, em_min, min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          tau_max_zt(i,k) = Lscale_max(i) / sqrt( em_zt(i,k) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_max_zm(i,k) = Lscale_max(i) / sqrt( em_clipped(i,k) )
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          tau_max_zt(i,k) = Lscale_max(i) / sqrt( em_zt(i,k) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_max_zm(i,k) = Lscale_max(i) / sqrt( max( em(i,k), em_min ) )
        end do
      end do
      !$acc end parallel loop

    end if

    if ( l_smooth_min_max ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_zm_unclipped(i,k) = one / invrs_tau_zm(i,k)
        end do
      end do
      !$acc end parallel loop

      tau_zm = smooth_min( nzm, ngrdcol, tau_zm_unclipped, &
                           tau_max_zm, 1.0e3_core_rknd * min_max_smth_mag )

      tau_zt_unclipped = zm2zt_api( nzm, nzt, ngrdcol, gr, tau_zm )

      tau_zt = smooth_min( nzt, ngrdcol, tau_zt_unclipped, &
                           tau_max_zt, 1.0e3_core_rknd * min_max_smth_mag )

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_zm(i,k) = min( one / invrs_tau_zm(i,k), tau_max_zm(i,k) )
        end do
      end do
      !$acc end parallel loop

      tau_zt = zm2zt_api( nzm, nzt, ngrdcol, gr, tau_zm )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          tau_zt(i,k) = min( tau_zt(i,k), tau_max_zt(i,k) )
        end do
      end do
      !$acc end parallel loop

    end if

    invrs_tau_zt     = zm2zt_api( nzm, nzt, ngrdcol, gr, invrs_tau_zm )
    invrs_tau_wp3_zt = zm2zt_api( nzm, nzt, ngrdcol, gr, invrs_tau_wp3_zm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        Lscale(i,k) = tau_zt(i,k) * sqrt( em_zt(i,k) )

        ! Lscale_up and Lscale_down aren't calculated with this option.
        ! They are set to 0 for stats output.
        Lscale_up(i,k) = zero
        Lscale_down(i,k) = zero

      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete( brunt_freq_pos, brunt_freq_out_cloud, &
    !$acc                   bvf_thresh, H_invrs_tau_wpxp_N2, ustar, &
    !$acc                   norm_ddzt_umvm, smooth_norm_ddzt_umvm, &
    !$acc                   brunt_vaisala_freq_clipped, &
    !$acc                   ice_supersat_frac_zm, invrs_tau_shear_smooth, &
    !$acc                   em_zt, tmp_calc_ngrdcol )

    !$acc exit data if( l_smooth_min_max ) &
    !$acc           delete( tau_zm_unclipped, tau_zt_unclipped, Ri_zm_smooth, em_clipped, &
    !$acc                   tmp_calc, tmp_calc_max, tmp_calc_min_max )

    return

  end subroutine diagnose_Lscale_from_tau

  !=============================================================================

  subroutine compute_coherent_mixing_length( nzm, nzt, ngrdcol, gr, thvm, thlm, &
                                             rtm, pdf_params, em, Lscale_max, p_in_Pa, &
                                             exner, thv_ds, Lscale_mu_entrain, &
                                             Lscale_tke_multiplier, &
                                             Lscale_pdf_launch_frac, &
                                             Lscale_up_weight, lmin, &
                                             l_implemented, &
                                             err_info, &
                                             Lscale, Lscale_up, Lscale_down )

    ! Description:
    !   Diagnose mixing length by weighting entraining moist-parcel travel with
    !   the fraction of assumed launch energies able to cross each energy barrier.
    !   Read thermodynamic-level environmental and PDF properties, momentum-level
    !   TKE, positive entrainment rates [1/m], and per-column launch/mixing weights.
    !   Return combined, upward, and downward lengths on thermodynamic levels [m];
    !   err_info is passed through unchanged.
    !
    !   Analytic layer updates advance each trajectory in path order. With launch
    !   energy scale E and running adverse-work barrier H [m^2/s^2], endpoint
    !   coherence is (E/(E+H))^2. Integrate this weight assuming H is linear across
    !   the layer, include the first layer at mean coherence <= 0.01, then stop.
    !   A source constant selects local, maximum-endpoint, or arrival-weighted
    !   directional lengths before floors, geometric mixing, and the final cap.
    !   Parcel saturation always uses the local Flatau polynomial.
    !
    ! References:
    !   doc/coherence_mixing_length.md: survival law and aggregation equations.
    !   Golaz et al. (2002), J. Atmos. Sci., 59, 3540--3551: parcel background.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        Cp, &
        Rd, &
        ep, &
        ep1, &
        ep2, &
        Lv, &
        grav, &
        zero_threshold, &
        one, &
        one_half, &
        zero

    use grid_class, only: &
        grid, &
        zm2zt_api

    use pdf_parameter_module, only: &
        pdf_parameter

    use clubb_precision, only: &
        core_rknd

    use err_info_type_module, only: &
        err_info_type

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm,     & ! Number of momentum levels
      nzt,     & ! Number of thermodynamic levels
      ngrdcol    ! Number of grid columns

    type(grid), intent(in) :: &
      gr   ! Grid structure

    type(pdf_parameter), intent(in) :: &
      pdf_params   ! PDF parameters used for directional launch properties

    logical, intent(in) :: &
      l_implemented   ! True when CLUBB is run within a host model

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      thvm,     & ! Environmental virtual potential temperature [K]
      thlm,     & ! Environmental liquid-water potential temperature [K]
      rtm,      & ! Environmental total-water mixing ratio [kg/kg]
      p_in_Pa,  & ! Pressure [Pa]
      exner,    & ! Exner function [-]
      thv_ds      ! Reference-state virtual potential temperature [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      em   ! Turbulent kinetic energy on momentum levels [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      Lscale_max,  & ! Maximum allowable mixing length [m]
      Lscale_mu_entrain,                    & ! Fractional entrainment rate [1/m]
      Lscale_tke_multiplier, & ! Mean launch energy divided by TKE [-]
      Lscale_pdf_launch_frac,    & ! PDF departure retained in launch properties [-]
      Lscale_up_weight           ! Upward weight in directional geometric mixing [-]

    real( kind = core_rknd ), intent(in) :: &
      lmin   ! Surface minimum mixing length [m]

    !----------------------- Input/Output Variables ------------------------
    type(err_info_type), intent(inout) :: &
      err_info   ! Error information

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Lscale,      & ! Turbulent mixing length [m]
      Lscale_up,   & ! Upward component of turbulent mixing length [m]
      Lscale_down    ! Downward component of turbulent mixing length [m]

    !--------------------------- Local Variables ---------------------------
    integer, parameter :: &
      no_overwrite = 0,                & ! Retain each locally launched reach
      legacy_max_overwrite = 1,        & ! Extend the most distant parcel endpoint
      coherence_average_overwrite = 2, & ! Average remaining reach by arrival coherence
      overwrite_mode = coherence_average_overwrite ! Selected overwrite method

    real( kind = core_rknd ), parameter :: &
      zlmin = 0.1_core_rknd,                      & ! Directional floor after aggregation [m]
      surface_layer_depth = 500._core_rknd,       & ! Depth of surface floor [m]
      invs_sfclyr_depth = one / surface_layer_depth ! Inverse surface-floor depth [1/m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      coherence_energy,       & ! Coherent parcel energy scale [m^2/s^2]
      crossed_length_up,      & ! Upward coherence integral crossing each layer [m]
      crossed_length_down,    & ! Downward coherence integral crossing each layer [m]
      arrival_coherence_up,   & ! Upward parcel coherence at each target [-]
      arrival_coherence_down, & ! Downward parcel coherence at each target [-]
      layer_retention,        & ! Incoming parcel-property retention factor [-]
      thl_entrainment,        & ! Environmental contribution to parcel thl [K]
      rt_entrainment,         & ! Environmental contribution to parcel rt [kg/kg]
      thl_launch,             & ! Directional parcel launch thl [K]
      rt_launch,              & ! Directional parcel launch rt [kg/kg]
      launch_buoyancy,        & ! Directional parcel launch buoyancy [m/s^2]
      grav_coef,              & ! Gravity divided by environmental thv [m/s^2/K]
      vapor_coef,             & ! Water-vapor contribution to parcel thv [K]
      latent_coef               ! Liquid-water contribution to parcel thv [K]

    integer :: &
      i, j, k, prev_j, & ! Grid indices
      step,            & ! Paired directional postprocessing step
      up_k,            & ! Bottom-up target level
      down_k             ! Top-down target level

    real( kind = core_rknd ) :: &
      layer_depth,          & ! Distance across the current layer [m]
      retention,            & ! Incoming parcel-property retention factor [-]
      response,             & ! Response to a layer-linear environment [-]
      thl_parcel,           & ! Parcel liquid-water potential temperature [K]
      rt_parcel,            & ! Parcel total-water mixing ratio [kg/kg]
      liquid_temperature,   & ! Parcel liquid temperature, then its square [K, K^2]
      saturation_rt,        & ! Parcel saturation mixing ratio [kg/kg]
      rc_parcel,            & ! Parcel liquid-water mixing ratio [kg/kg]
      oriented_buoyancy,    & ! Buoyancy in the path direction [m/s^2]
      prev_buoyancy,        & ! Oriented buoyancy at the prior point [m/s^2]
      prev_work,            & ! Cumulative buoyancy work [m^2/s^2]
      max_barrier,          & ! Largest adverse-work barrier encountered [m^2/s^2]
      prev_barrier,         & ! Adverse-work barrier at the layer start [m^2/s^2]
      layer_coherence,      & ! Mean survival weight in the layer [-]
      endpoint_coherence,   & ! Survival weight at the layer endpoint [-]
      layer_length,         & ! Coherence integral through the current layer [m]
      Lscale_tmp,           & ! Coherence-weighted parcel travel distance [m]
      remaining_length_up,  & ! Total remaining upward reach [m]
      remaining_length_down, & ! Total remaining downward reach [m]
      maximum_endpoint,     & ! Highest upward parcel endpoint [m]
      minimum_endpoint,     & ! Lowest downward parcel endpoint [m]
      surface_floor,        & ! Height-dependent minimum mixing length [m]
      liquid_coef             ! Constant in the saturation-excess formula [K^2]

    !----------------------------- Begin Code ------------------------------

    !$acc enter data create( coherence_energy, layer_retention, &
    !$acc                    thl_entrainment, rt_entrainment, &
    !$acc                    thl_launch, rt_launch, launch_buoyancy, &
    !$acc                    grav_coef, vapor_coef, latent_coef )

    if ( overwrite_mode == coherence_average_overwrite ) then
      !$acc enter data create( crossed_length_up, crossed_length_down, &
      !$acc                    arrival_coherence_up, arrival_coherence_down )
    end if

    ! Interpolate launch TKE to thermodynamic levels.
    coherence_energy = zm2zt_api( nzm, nzt, ngrdcol, gr, em )
    liquid_coef = ep * Lv**2 / ( Rd * Cp )

    ! Precompute the inverse-square energy scale and buoyancy coefficients.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        ! This is the value we use for the mean of the assumed parcel energy distribution.
        ! It is a somewhat arbitrary choice, so it was given a scaling parameter
        coherence_energy(i,k) = coherence_energy(i,k) * Lscale_tke_multiplier(i)

        ! These coefficients reconstruct parcel thv and scale its buoyancy.
        grav_coef(i,k)   = grav / thvm(i,k)
        vapor_coef(i,k)  = ep1  * thv_ds(i,k)
        latent_coef(i,k) = Lv   / ( exner(i,k) * Cp ) - ep2 * thv_ds(i,k)

      end do
    end do
    !$acc end parallel loop

    ! The coherence-average mode needs directional crossing sums. At each target,
    ! its own local launch contributes unit arrival coherence, even at zero TKE.
    if ( overwrite_mode == coherence_average_overwrite ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          crossed_length_up(i,k)      = zero
          crossed_length_down(i,k)    = zero
          arrival_coherence_up(i,k)   = one
          arrival_coherence_down(i,k) = one
        end do
      end do
      !$acc end parallel loop
    end if

    !-------------------- Upward calculation ---------------------

    ! Precompute the upward entrainment update for every layer.
    !$acc parallel loop gang vector collapse(2) default(present)
    do j = gr%k_lb_zt + gr%grid_dir_indx, gr%k_ub_zt, gr%grid_dir_indx
      do i = 1, ngrdcol
        prev_j = j - gr%grid_dir_indx
        layer_depth = abs( gr%zt(i,j) - gr%zt(i,prev_j) )

        ! Integrate entrainment analytically through a layer-linear environment.
        retention = exp( -Lscale_mu_entrain(i) * layer_depth )
        response = ( one - retention ) / ( Lscale_mu_entrain(i) * layer_depth )

        layer_retention(i,j) = retention
        thl_entrainment(i,j) = thlm(i,j) - retention * thlm(i,prev_j) &
                               - response * ( thlm(i,j) - thlm(i,prev_j) )
        rt_entrainment(i,j)  = rtm(i,j) - retention * rtm(i,prev_j) &
                               - response * ( rtm(i,j) - rtm(i,prev_j) )
      end do
    end do
    !$acc end parallel loop

    ! Diagnose every upward launch state and buoyancy in parallel.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        ! Blend toward PDF component 1 to define the upward launch state.
        ! Using Lscale_pdf_launch_frac > 0 is experimental, and seems to introduce
        ! a tight feedback loop that can cause instability
        thl_launch(i,k) = thlm(i,k) + Lscale_pdf_launch_frac(i) &
                                    * ( pdf_params%thl_1(i,k) - thlm(i,k) )
        rt_launch(i,k) = rtm(i,k) + Lscale_pdf_launch_frac(i) &
                                  * ( pdf_params%rt_1(i,k) - rtm(i,k) )

        ! Saturation-adjust the launch state before diagnosing its buoyancy.
        liquid_temperature = thl_launch(i,k) * exner(i,k)
        saturation_rt      = sat_mixrat_liq_flatau( p_in_Pa(i,k), liquid_temperature )
        liquid_temperature = liquid_temperature**2

        rc_parcel = max( ( rt_launch(i,k) - saturation_rt ) * liquid_temperature &
                       / ( liquid_temperature + liquid_coef * saturation_rt ), &
                         zero_threshold )

        launch_buoyancy(i,k) = grav_coef(i,k) &
                             * ( thl_launch(i,k) + vapor_coef(i,k) * rt_launch(i,k) &
                               + latent_coef(i,k) * rc_parcel ) - grav
      end do
    end do
    !$acc end parallel loop

    ! Compute the coherence-weighted upward length for each launch level.
    !$acc parallel loop collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        Lscale_tmp = zero

        ! Initialize the parcel state and path-integral trackers at launch.
        thl_parcel  = thl_launch(i,k)
        rt_parcel   =  rt_launch(i,k)

        prev_buoyancy = launch_buoyancy(i,k)
        prev_work     = zero
        prev_barrier  = zero

        ! The upper boundary is outside the active mixing-length domain.
        do j = k + gr%grid_dir_indx, gr%k_ub_zt - gr%grid_dir_indx, gr%grid_dir_indx

          layer_depth = abs( gr%zt(i,j) - gr%zt(i,j - gr%grid_dir_indx) )

          ! Entrain layer-mean air into the parcel as it rises.
          thl_parcel  = thl_entrainment(i,j) + layer_retention(i,j) * thl_parcel
          rt_parcel   = rt_entrainment(i,j) + layer_retention(i,j) * rt_parcel

          ! Diagnose cloud liquid water for the parcel buoyancy calculation.
          liquid_temperature = thl_parcel * exner(i,j)
          saturation_rt = sat_mixrat_liq_flatau( p_in_Pa(i,j), liquid_temperature )
          liquid_temperature = liquid_temperature**2
          rc_parcel = ( rt_parcel - saturation_rt ) * liquid_temperature &
                      / ( liquid_temperature + liquid_coef * saturation_rt )
          rc_parcel = max( rc_parcel, zero_threshold )

          ! Diagnose upward-oriented buoyancy: B = g * ( thv_parcel / thvm - 1 ).
          oriented_buoyancy = grav_coef(i,j) &
                            * ( thl_parcel + vapor_coef(i,j) * rt_parcel &
                              + latent_coef(i,j) * rc_parcel ) - grav

          ! Integrate buoyancy over the layer to update cumulative parcel work.
          prev_work = prev_work + one_half * ( prev_buoyancy + oriented_buoyancy ) &
                    * layer_depth

          ! The barrier is the greatest energy deficit encountered along the path.
          max_barrier = max( prev_barrier, -prev_work )

          ! Integrate survival exactly for a barrier linear between these endpoints.
          layer_coherence = coherence_energy(i,k)**2 &
                          / ( ( coherence_energy(i,k) + prev_barrier ) &
                            * ( coherence_energy(i,k) + max_barrier ) )

          layer_length = layer_coherence * layer_depth

          if ( overwrite_mode == coherence_average_overwrite ) then

            endpoint_coherence = coherence_energy(i,k)**2 &
                                / ( coherence_energy(i,k) + max_barrier )**2

            ! These sums reconstruct remaining reach without storing trajectories.
            !$acc atomic update
            crossed_length_up(i,j) = crossed_length_up(i,j) + layer_length
            !$acc atomic update
            arrival_coherence_up(i,j) = arrival_coherence_up(i,j) &
                                      + endpoint_coherence
          end if

          ! Accumulate coherence-weighted distance. Coherence is one without a barrier
          ! and cannot increase because the running barrier cannot decrease.
          Lscale_tmp = Lscale_tmp + layer_length

          ! Include this layer, then neglect the remaining low-coherence tail.
          ! This has almost no effect on the result, and reduces the runtime cost
          if ( layer_coherence <= 0.01 ) exit

          prev_buoyancy = oriented_buoyancy
          prev_barrier  = max_barrier

        end do

        Lscale_up(i,k) = Lscale_tmp

      end do
    end do
    !$acc end parallel loop

    !------------------- Downward calculation --------------------

    ! Reuse the same scratch arrays for the downward layer updates.
    !$acc parallel loop gang vector collapse(2) default(present)
    do j = gr%k_ub_zt - gr%grid_dir_indx, gr%k_lb_zt, -gr%grid_dir_indx
      do i = 1, ngrdcol
        prev_j = j + gr%grid_dir_indx
        layer_depth = abs( gr%zt(i,j) - gr%zt(i,prev_j) )

        ! Integrate entrainment analytically through a layer-linear environment.
        retention = exp( -Lscale_mu_entrain(i) * layer_depth )
        response = ( one - retention ) / ( Lscale_mu_entrain(i) * layer_depth )

        layer_retention(i,j) = retention
        thl_entrainment(i,j) = thlm(i,j) - retention * thlm(i,prev_j) &
                             - response * ( thlm(i,j) - thlm(i,prev_j) )
        rt_entrainment(i,j) = rtm(i,j) - retention * rtm(i,prev_j) &
                            - response * ( rtm(i,j) - rtm(i,prev_j) )
      end do
    end do
    !$acc end parallel loop

    ! Diagnose every downward launch state and buoyancy in parallel.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        ! Blend toward PDF component 2 to define the downward launch state.
        ! Using Lscale_pdf_launch_frac > 0 is experimental, and seems to introduce
        ! a tight feedback loop that can cause instability
        thl_launch(i,k) = thlm(i,k) + Lscale_pdf_launch_frac(i) &
                                    * ( pdf_params%thl_2(i,k) - thlm(i,k) )
        rt_launch(i,k) = rtm(i,k) + Lscale_pdf_launch_frac(i) &
                                  * ( pdf_params%rt_2(i,k) - rtm(i,k) )

        ! Saturation-adjust the launch state before diagnosing its buoyancy.
        liquid_temperature = thl_launch(i,k) * exner(i,k)
        saturation_rt = sat_mixrat_liq_flatau( p_in_Pa(i,k), liquid_temperature )
        liquid_temperature = liquid_temperature**2

        rc_parcel = max( ( rt_launch(i,k) - saturation_rt ) * liquid_temperature &
                       / ( liquid_temperature + liquid_coef * saturation_rt ), &
                         zero_threshold )

        launch_buoyancy(i,k) = grav - grav_coef(i,k) &
                             * ( thl_launch(i,k) + vapor_coef(i,k) * rt_launch(i,k) &
                               + latent_coef(i,k) * rc_parcel )
      end do
    end do
    !$acc end parallel loop

    ! Compute the coherence-weighted downward length for each launch level.
    !$acc parallel loop collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        Lscale_tmp = zero

        ! Initialize the parcel state and path-integral trackers at launch.
        thl_parcel  = thl_launch(i,k)
        rt_parcel   =  rt_launch(i,k)

        prev_buoyancy = launch_buoyancy(i,k)
        prev_work     = zero
        prev_barrier  = zero

        ! Track the parcel from its launch level down to the lower boundary.
        do j = k - gr%grid_dir_indx, gr%k_lb_zt, -gr%grid_dir_indx

          layer_depth = abs( gr%zt(i,j) - gr%zt(i,j + gr%grid_dir_indx) )

          ! Entrain layer-mean air into the parcel as it descends.
          thl_parcel = thl_entrainment(i,j) + layer_retention(i,j) * thl_parcel
          rt_parcel = rt_entrainment(i,j) + layer_retention(i,j) * rt_parcel

          ! Diagnose cloud liquid water for the parcel buoyancy calculation.
          liquid_temperature = thl_parcel * exner(i,j)
          saturation_rt = sat_mixrat_liq_flatau( p_in_Pa(i,j), liquid_temperature )
          liquid_temperature = liquid_temperature**2
          rc_parcel = max( ( rt_parcel - saturation_rt ) * liquid_temperature &
                        / ( liquid_temperature + liquid_coef * saturation_rt ), &
                          zero_threshold )

          ! Reverse the physical buoyancy sign to orient it along the downward path.
          oriented_buoyancy = grav - grav_coef(i,j) &
                            * ( thl_parcel + vapor_coef(i,j) * rt_parcel &
                              + latent_coef(i,j) * rc_parcel )

          ! Integrate oriented buoyancy to update cumulative parcel work.
          prev_work = prev_work + one_half * ( prev_buoyancy + oriented_buoyancy ) &
                    * layer_depth

          ! The barrier is the greatest energy deficit encountered along the path.
          max_barrier = max( prev_barrier, -prev_work )

          ! Integrate survival for a barrier linear between the layer endpoints.
          layer_coherence = coherence_energy(i,k)**2 &
                          / ( ( coherence_energy(i,k) + prev_barrier ) &
                            * ( coherence_energy(i,k) + max_barrier ) )

          layer_length = layer_coherence * layer_depth

          ! Accumulate coherence-weighted distance before testing the tail cutoff.
          if ( overwrite_mode == coherence_average_overwrite &
                .and. k /= gr%k_ub_zt ) then

            endpoint_coherence = coherence_energy(i,k)**2 &
                                / ( coherence_energy(i,k) + max_barrier )**2
            !$acc atomic update
            crossed_length_down(i,j) = crossed_length_down(i,j) + layer_length
            !$acc atomic update
            arrival_coherence_down(i,j) = arrival_coherence_down(i,j) &
                                        + endpoint_coherence
          end if

          Lscale_tmp = Lscale_tmp + layer_length

          ! Include this layer, then neglect the remaining low-coherence tail.
          ! This has almost no effect on the result, and reduces the runtime cost
          if ( layer_coherence <= 0.01 ) exit

          prev_buoyancy = oriented_buoyancy
          prev_barrier  = max_barrier

        end do

        Lscale_down(i,k) = Lscale_tmp

      end do
    end do
    !$acc end parallel loop

    !-------------------- Aggregation ---------------------

    ! Do some final smoothing on the result
    if ( overwrite_mode == legacy_max_overwrite ) then

      ! The legacy_max_overwrite method replaces Lscale values with the
      ! maximum of all parcels that crossed it. For example, if a parcel
      ! that starts on k=10 rises to k=20, but one from k=11 only rises
      ! to k=15, this overwriting method will boost the k=11 result as
      ! if it rose to k=20.

      ! This was created primarily to boost the final Lscale result,
      ! but also has a smoothing effect.

      !$acc parallel loop gang default(present)
      do i = 1, ngrdcol
        maximum_endpoint = gr%zt(i,gr%k_lb_zt) &
                         + max( zlmin, Lscale_up(i,gr%k_lb_zt) )
        minimum_endpoint = gr%zt(i,gr%k_ub_zt) &
                         - max( zlmin, Lscale_down(i,gr%k_ub_zt) )

        !$acc loop seq
        do step = 0, nzt - 2
          up_k = gr%k_lb_zt + step * gr%grid_dir_indx
          down_k = gr%k_ub_zt - ( step + 1 ) * gr%grid_dir_indx

          Lscale_up(i,up_k) = max( zlmin, Lscale_up(i,up_k) )
          Lscale_down(i,down_k) = max( zlmin, Lscale_down(i,down_k) )

          ! Match the legacy active bounds at the surface and upper boundary.
          if ( step > 0 .and. step < nzt - 2 ) then
            Lscale_up(i,up_k) = max( Lscale_up(i,up_k), &
                                    maximum_endpoint - gr%zt(i,up_k) )
            maximum_endpoint = max( maximum_endpoint, &
                                    gr%zt(i,up_k) + Lscale_up(i,up_k) )
          end if
          if ( step < nzt - 2 ) then
            Lscale_down(i,down_k) = max( Lscale_down(i,down_k), &
                                        gr%zt(i,down_k) - minimum_endpoint )
            minimum_endpoint = min( minimum_endpoint, &
                                    gr%zt(i,down_k) - Lscale_down(i,down_k) )
          end if
        end do
      end do
      !$acc end parallel loop

    else if ( overwrite_mode == coherence_average_overwrite ) then

      ! The coherence_average_overwrite averages how far arriving parcels can still
      ! travel, weighted by their coherence at the level.
      ! For an upward example on an evenly spaced grid, suppose the mixing lengths
      ! from k=10,11,12,13 would reach k=20, and those from k=14,15 would reach k=18.
      ! If all arriving parcels have coherence 1 and none from below k=10 arrive,
      ! the averaged lengths at k=14,15 correspond to reaching
      ! k=(4/5)*20+(1/5)*18 and k=(4/6)*20+(2/6)*18 respectively.

      !$acc parallel loop gang default(present)
      do i = 1, ngrdcol
        remaining_length_up = Lscale_up(i,gr%k_lb_zt)
        remaining_length_down = Lscale_down(i,gr%k_ub_zt - gr%grid_dir_indx)

        !$acc loop seq
        do step = 0, nzt - 2
          up_k = gr%k_lb_zt + step * gr%grid_dir_indx
          down_k = gr%k_ub_zt - ( step + 1 ) * gr%grid_dir_indx

          if ( step > 0 ) then
            remaining_length_up = remaining_length_up - crossed_length_up(i,up_k) &
                                + Lscale_up(i,up_k)
            remaining_length_down = remaining_length_down &
                                  - crossed_length_down(i,down_k) &
                                  + Lscale_down(i,down_k)
          end if
          Lscale_up(i,up_k) = max( zlmin, remaining_length_up &
                                 / arrival_coherence_up(i,up_k) )
          Lscale_down(i,down_k) = max( zlmin, remaining_length_down &
                                     / arrival_coherence_down(i,down_k) )
        end do
      end do
      !$acc end parallel loop

    else if ( overwrite_mode == no_overwrite ) then

      ! No overwrite, just ensure zlmin is respected

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          Lscale_up(i,k)   = max( zlmin, Lscale_up(i,k) )
          Lscale_down(i,k) = max( zlmin, Lscale_down(i,k) )
        end do
      end do
      !$acc end parallel loop

    end if

    ! Apply surface constraints and combine the two directional lengths.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = gr%k_lb_zt, gr%k_ub_zt - gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        ! The minimum length decays linearly over the lowest 500 m.
        if ( l_implemented ) then
          surface_floor = ( surface_layer_depth - ( gr%zt(i,k) - gr%zm(i,gr%k_lb_zm) ) ) &
                          * lmin * invs_sfclyr_depth
        else
          surface_floor = ( surface_layer_depth - gr%zt(i,k) ) &
                          * lmin * invs_sfclyr_depth
        end if

        ! Apply the same surface floor before combining the directional lengths.
        Lscale_up(i,k) = max( surface_floor, Lscale_up(i,k) )
        Lscale_down(i,k) = max( surface_floor, Lscale_down(i,k) )

        ! Weighted geometric mean between Lscale_up and Lscale_down.
        Lscale(i,k) = Lscale_up(i,k)**Lscale_up_weight(i) &
                    * Lscale_down(i,k)**( one - Lscale_up_weight(i) )

        ! Enforce the configured upper bound on the combined mixing length.
        Lscale(i,k) = min( Lscale(i,k), Lscale_max(i) )

        ! The highest active level also supplies the zero-gradient lid value.
        if ( k == gr%k_ub_zt - gr%grid_dir_indx ) then
          Lscale(i,gr%k_ub_zt) = Lscale(i,k)
        end if

      end do
    end do
    !$acc end parallel loop

    if ( overwrite_mode == coherence_average_overwrite ) then
      !$acc exit data delete( crossed_length_up, crossed_length_down, &
      !$acc                   arrival_coherence_up, arrival_coherence_down )
    end if

    !$acc exit data delete( coherence_energy, &
    !$acc                   layer_retention, &
    !$acc                   thl_entrainment, rt_entrainment, &
    !$acc                   thl_launch, rt_launch, launch_buoyancy, &
    !$acc                   grav_coef, vapor_coef, latent_coef )

    return

  end subroutine compute_coherent_mixing_length


  !=============================================================================

  pure function sat_mixrat_liq_flatau( pressure, temperature ) result( saturation_rt )
  !$acc routine seq

    ! Description:
    !   Compute liquid saturation mixing ratio using only the Flatau polynomial.
    !   This compact specialization is colocated with the mixing-length closure
    !   so the compiler can inline it into parcel trajectories.
    !
    !   It is also intentionally the serial method as opposed to the
    !   complex pipeline / multi-fma unit friendlier version. This was
    !   to minimize register usage in the complex loop it is used in.
    !   See saturation.F90 for more details.
    !
    ! Note:
    !   It is usually frowned upon to make custom saturation equations like
    !   this because it can cause different parts of code to detect cloud
    !   transitions differently, leading to inconsistent calculations.
    !   We believe that there is no similar concern in this instance
    !   because this is used only to compute the behavior of "virtual" parcels
    !   in the above methods. Virtual in the sense that their lifetime is
    !   only inside the Lscale calculation, which is already somewhat
    !   arbitrary, and none of that internal information escapes.
    !
    ! Inputs and output:
    !   pressure [Pa] and temperature [K] are scalar inputs; saturation_rt [kg/kg]
    !   is returned. Evaluate the polynomial no lower than -85 deg C. If pressure
    !   minus saturation vapor pressure is less than 1 Pa, return ep (the ratio
    !   of dry-air to water-vapor gas constants) instead of dividing.
    !
    ! References:
    !   Flatau et al. (1992), J. Appl. Meteor., 31, 1507--1513.
    !   saturation.F90: shared saturation formulas and polynomial evaluations.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        ep, &
        T_freeze_K

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), intent(in) :: &
      pressure,     & ! Pressure [Pa]
      temperature    ! Temperature [K]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ) :: &
      saturation_rt   ! Saturation mixing ratio [kg/kg]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), parameter :: &
      minimum_temperature = -85._core_rknd   ! Polynomial lower bound [deg C]

    real( kind = core_rknd ) :: &
      temperature_celsius, & ! Clipped temperature [deg C]
      saturation_pressure    ! Saturation vapor pressure [Pa]

    !----------------------------- Begin Code ------------------------------
    temperature_celsius = max( temperature - T_freeze_K, minimum_temperature )

    saturation_pressure = 100._core_rknd &
      * ( 6.11583699_core_rknd &
        + temperature_celsius * ( 0.444606896_core_rknd &
        + temperature_celsius * ( 0.143177157e-1_core_rknd &
        + temperature_celsius * ( 0.264224321e-3_core_rknd &
        + temperature_celsius * ( 0.299291081e-5_core_rknd &
        + temperature_celsius * ( 0.203154182e-7_core_rknd &
        + temperature_celsius * ( 0.702620698e-10_core_rknd &
        + temperature_celsius * ( 0.379534310e-13_core_rknd &
        - temperature_celsius *   0.321582393e-15_core_rknd ) ) ) ) ) ) ) )

    if ( pressure - saturation_pressure < 1._core_rknd ) then
      saturation_rt = ep
    else
      saturation_rt = ep * saturation_pressure / ( pressure - saturation_pressure )
    end if

    return

  end function sat_mixrat_liq_flatau

end module mixing_length

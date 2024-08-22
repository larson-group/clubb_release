! $Id$
!===============================================================================
module microphys_driver

  ! Description:
  ! Call a microphysical scheme and output tendencies for mean precipitating
  ! hydrometeors, mean cloud droplet concentration, predictive mean, variance,
  ! and covariance fields, as well as sedimentation velocity variables for
  ! precipitating hydrometeors.

  ! References:
  !  H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
  !  moment microphysics scheme for application in cloud and
  !  climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.

  !  Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
  !  large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
  !  229-243, 2000.
  !-------------------------------------------------------------------------

  implicit none

  ! Subroutines
  public :: calc_microphys_scheme_tendcies

  private ! Default Scope

  contains

  !=============================================================================
  subroutine calc_microphys_scheme_tendcies( gr, dt, time_current,      & ! In
                                pdf_dim, hydromet_dim, runtype,         & ! In
                                thlm, p_in_Pa, exner, rho, rho_zm, rtm, & ! In
                                rcm, cloud_frac, wm_zt, wm_zm, wp2_zt, &  ! In
                                hydromet, Nc_in_cloud, &                  ! In
                                hm_metadata, &                            ! In
                                pdf_params, hydromet_pdf_params, &        ! In
                                precip_fracs, &                           ! In
                                X_nl_all_levs, X_mixt_comp_all_levs, &    ! In
                                lh_sample_point_weights, &                ! In
                                mu_x_1_n, mu_x_2_n, &                     ! In
                                sigma_x_1_n, sigma_x_2_n, &               ! In
                                corr_array_1_n, corr_array_2_n, &         ! In
                                lh_rt_clipped, lh_thl_clipped, &          ! In
                                lh_rc_clipped, lh_rv_clipped, &           ! In
                                lh_Nc_clipped, &                          ! In
                                l_lh_importance_sampling, &               ! In
                                l_lh_instant_var_covar_src, &             ! In
                                saturation_formula, &                     ! In
                                stats_metadata, &                         ! In
                                stats_zt, stats_zm, &                     ! Inout
                                stats_sfc, stats_lh_zt, &                 ! Inout
                                Nccnm, &                                  ! Inout
                                hydromet_mc, Ncm_mc, rcm_mc, rvm_mc, &    ! Out
                                thlm_mc, hydromet_vel_zt, &               ! Out
                                hydromet_vel_covar_zt_impc, &             ! Out
                                hydromet_vel_covar_zt_expc, &             ! Out
                                wprtp_mc, wpthlp_mc,  rtp2_mc, &          ! Out
                                thlp2_mc, rtpthlp_mc )                    ! Out

    ! Description:
    ! Call a microphysics scheme and output microphysics tendencies for the
    ! predictive variables.

    ! References:
    ! H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
    ! moment microphysics scheme for application in cloud and
    ! climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.
    !
    ! Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
    ! large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
    ! 229-243, 2000.
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        zt2zm    ! Procedure(s)

    use grid_class, only: grid ! Type

    use constants_clubb, only: & 
        one,            & ! Constant(s)
        zero,           &
        cm3_per_m3,     &
        cloud_frac_min

    use KK_microphys_module, only: & 
        KK_local_microphys,    & ! Procedure(s)
        KK_upscaled_microphys

    use cloud_sed_module, only: &
        cloud_drop_sed  ! Procedure(s)

    use morrison_microphys_module, only: &
        morrison_microphys_driver  ! Procedure(s)

#ifdef COAMPS_MICRO
    use coamps_microphys_driver_module, only:  & 
        coamps_microphys_driver  ! Procedure(s)
#endif

#ifdef SILHS
    use lh_microphys_driver_module, only: &
        lh_microphys_driver  ! Procedure(s)
#endif /* SILHS */

    use ice_dfsn_module, only: & 
        ice_dfsn  ! Procedure(s)

    use T_in_K_module, only: &
        thlm2T_in_K  ! Procedure(s)

    use gfdl_activation, only: &
        aer_act_clubb_quadrature_Gauss, & ! Procedure(s)
        aeromass_value                    ! Variable(s)

    use clubb_api_module, only: &
      update_xp2_mc_api

    use model_flags, only: &
        l_morr_xp2_mc  ! Flag(s)

    use parameters_microphys, only: &
        l_cloud_sed,          & ! Cloud water sedimentation (K&K or no microphysics)
        l_predict_Nc,         & ! Predict cloud droplet number conc (Morrison)
        lh_num_samples,   & ! # of SILHS samples for which call the microphysics
        l_local_kk,           & ! Use local formula for K&K
        microphys_scheme,     & ! The microphysical scheme in use
        l_gfdl_activation,    & ! Flag to use GFDL activation scheme
        microphys_start_time, & ! When to start the microphysics [s]
        sigma_g,              & ! Parameter used in the cloud droplet sedimentation code
        l_silhs_KK_convergence_adj_mean ! Flag used to adjust a run of KK upscaled or SILHS

    use pdf_parameter_module, only:  &
        pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only:  &
        hydromet_pdf_parameter, &  ! Type
        precipitation_fractions

    use stats_type_utilities, only: & 
        stat_update_var,   & ! Procedure(s)
        stat_begin_update, &
        stat_end_update

    use corr_varnce_module, only: &
        hm_metadata_type

    use stats_variables, only: &
        stats_metadata_type

    use stats_clubb_utilities, only: & 
        stats_accumulate_lh_tend ! Procedure(s)

    use parameters_microphys, only: &
        lh_microphys_type,        & ! Determines how the LH samples are used
        lh_microphys_interactive, & ! Feed the subcols into microphys and allow feedback
        lh_microphys_disabled       ! Disable latin hypercube entirely

    use clubb_precision, only: &
        time_precision, & ! Constant(s)
        core_rknd

    use microphys_stats_vars_module, only: &
        microphys_stats_vars_type,  & ! Type
        microphys_stats_accumulate, & ! Procedure(s)
        microphys_get_var,          &
        microphys_stats_cleanup

    use stats_type, only: stats ! Type

    implicit none

    ! Input Variables
    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in) ::  & 
      dt           ! Timestep         [s]

    real( kind = time_precision ), intent(in) ::  & 
      time_current ! Current time     [s]

    integer, intent(in) :: &
      pdf_dim,      & ! Number of variables in the correlation arrays
      hydromet_dim    ! Number of hydrometeor species

    character(len=*), intent(in) :: & 
      runtype ! Name of the run, for case specific effects.

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: & 
      thlm,       & ! Liquid potential temp.                    [K]
      p_in_Pa,    & ! Pressure                                  [Pa]
      exner,      & ! Exner function                            [-]
      rho,        & ! Density on thermodynamic levels           [kg/m^3]
      rtm,        & ! Total water mixing ratio                  [kg/kg]
      rcm,        & ! Liquid water mixing ratio                 [kg/kg]
      cloud_frac, & ! Cloud fraction                            [-]
      wm_zt,      & ! w wind component on thermodynamic levels  [m/s]
      wp2_zt        ! w'^2 on the thermo. grid                  [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: & 
      rho_zm,     & ! Density on momentum levels                [kg/m^3]
      wm_zm         ! w wind component on momentum levels       [m/s]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet       ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(pdf_parameter), intent(in) :: & 
      pdf_params     ! PDF parameters

    type(hydromet_pdf_parameter), dimension(gr%nzt), intent(in) :: &
      hydromet_pdf_params     ! PDF parameters
      
    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(lh_num_samples,gr%nzt,pdf_dim), &
    intent(in) :: &
      X_nl_all_levs ! Normally and lognormally distributed hydrometeors and other variables

    integer, dimension(lh_num_samples,gr%nzt), intent(in) :: &
      X_mixt_comp_all_levs ! Which mixture component the sample is in

    real( kind = core_rknd ), dimension(lh_num_samples,gr%nzt), intent(in) :: &
      lh_sample_point_weights ! Weights for cloud weighted sampling

    real( kind = core_rknd ), dimension(gr%nzt,pdf_dim), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(gr%nzt, pdf_dim, pdf_dim), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]
      
    real( kind = core_rknd ), dimension(lh_num_samples,gr%nzt), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    logical, intent(in) :: &
      l_lh_importance_sampling, & ! Do importance sampling (SILHS) [-]
      l_lh_instant_var_covar_src  ! Produce instantaneous var/covar tendencies [-]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output Variables
    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt

    ! Note:
    ! For COAMPS Nccnm is initialized and Nim & Ncm are computed within
    ! subroutine adjtg.
    real( kind = core_rknd ), dimension(gr%nzt), intent(inout) :: & 
      Nccnm    ! Cloud condensation nuclei concentration (COAMPS)  [num/kg]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(out) :: & 
      hydromet_mc     ! Change in hydrometeors due to microphysics  [units/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: & 
      Ncm_mc     ! Change in Ncm due to microphysics  [num/kg/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: & 
      rcm_mc,  & ! Microphysics contributions to liquid water          [kg/kg/s]
      rvm_mc,  & ! Microphysics contributions to vapor water           [kg/kg/s]
      thlm_mc    ! Microphysics contributions to theta-l               [K/s]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(out) :: &
      hydromet_vel_zt   ! Mean hydrometeor sed. velocity on thermo. levs. [m/s]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(out) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_xx'x_x'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      wprtp_mc,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nzt) :: &
      delta_zt  ! Difference in thermo. height levels     [m]

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      T_in_K,        & ! Temperature                                [K]
      rvm,           & ! Vapor water mixing ratio                   [kg/kg]
      thlm_morr,     & ! Thlm fed into morrison microphysics        [K]
      Ncm_microphys    ! Mean cloud droplet concentration, <N_c>    [num/kg]

    real( kind = core_rknd ), dimension(gr%nzt) :: & 
      rrm_evap    ! Evaporation rate for rrm         [kg/kg/s]

    real( kind = core_rknd ), dimension(1,1,gr%nzt) :: & 
      cond ! COAMPS stat for condesation/evap of rcm

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      wtmp,    & ! Standard dev. of w                   [m/s]
      chi   ! The variable 's' in Mellor (1977)    [kg/kg]

    integer :: k    ! Loop index

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      Ndrop_max  ! GFDL droplet activation concentration [#/kg]

    !Input aerosol mass concentration: the unit is 10^12 ug/m3.
    !For example, aeromass=2.25e-12 means that the aerosol mass
    !concentration is 2.25 ug/m3.
    !This value of aeromass was recommended by Huan Guo
    real( kind = core_rknd ), dimension(gr%nzt, 4) :: &
      aeromass ! ug/m^3

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      lh_AKm,     & ! Kessler ac estimate                 [kg/kg/s]
      AKm,        & ! Exact Kessler ac                    [kg/kg/s]
      AKstd,      & ! St dev of exact Kessler ac          [kg/kg/s]
      AKstd_cld,  & ! Stdev of exact w/in cloud ac        [kg/kg/s]
      lh_rcm_avg, & ! Monte Carlo rcm estimate            [kg/kg]
      AKm_rcm,    & ! Kessler ac based on rcm             [kg/kg/s]
      AKm_rcc       ! Kessler ac based on rcm/cloud_frac  [kg/kg/s]

    type(microphys_stats_vars_type) :: &
      microphys_stats_zt, &   ! Stats. vars. from microphys. on zt and sfc grids
      microphys_stats_sfc

    logical :: l_latin_hypercube_input

    ! Just to avoid typing hm_metadata%iixx everywhere
    integer :: &
      iirr, &
      iirs, &
      iiri, &
      iirg, &
      iiNr, &
      iiNi

    !------------------------ Begin Code ------------------------

    iirr = hm_metadata%iirr
    iirs = hm_metadata%iirs
    iiri = hm_metadata%iiri
    iirg = hm_metadata%iirg
    iiNr = hm_metadata%iiNr
    iiNi = hm_metadata%iiNi

    ! We must initialize intent(out) variables. If we do not, and they do not
    ! otherwise get assigned (e.g. microphys_scheme=="none"), then we are not
    ! within the realms of Fortran standards compliance.
    rvm_mc = zero
    rcm_mc = zero
    thlm_mc = zero
    wprtp_mc = zero
    wpthlp_mc = zero
    rtp2_mc = zero
    thlp2_mc = zero
    rtpthlp_mc = zero

    ! Initialize predicitive hydrometeor tendencies and sedimentation
    ! velocities.
    if ( hydromet_dim > 0 ) then
      hydromet_mc = zero
      hydromet_vel_zt = zero
    endif
    ! Initialize predicitive Ncm tendency.
    Ncm_mc = zero

    ! Initialize the values of the implicit and explicit components used to
    ! calculate the covariances of hydrometeor sedimentation velocities and
    ! their associated hydrometeors (for example, <V_rr'r_r'> and <V_Nr'N_r'>)
    ! to 0.
    if ( hydromet_dim > 0 ) then
      hydromet_vel_covar_zt_impc = zero
      hydromet_vel_covar_zt_expc = zero
    endif

    ! Return if there is delay between the model start time and start of the
    ! microphysics
    if ( time_current < microphys_start_time ) return

    ! Make some compiler warnings go away for external users
    if ( runtype == "" ) then
      error stop "Runtype is null, which should not happen."
    endif

    ! Calculate the updated value of mean cloud droplet concentration, based on
    ! Nc_in_cloud and the updated value of cloud_frac.
    Ncm_microphys = Nc_in_cloud * cloud_frac

    Ndrop_max = zero

    ! Set the value of the aerosol mass array to a constant
    aeromass = aeromass_value

    ! Determine 's' from Mellor (1977)
    chi(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
                  + ( one - pdf_params%mixt_frac(1,:) ) * pdf_params%chi_2(1,:)

    ! Compute standard deviation of vertical velocity in the grid column
    wtmp(:) = sqrt( wp2_zt(:) )

    ! Compute difference in thermodynamic height levels.
    ! The variable delta_zt is on the Morrison Microphysics grid, which
    ! identifies delta_zt at level k as being the difference between zt(k+1) and
    ! zt(k). On CLUBB's grid, the difference between zt(k+1) and zt(k) is
    ! identified as gr%dzm at level k+1 (since momentum level k+1 is between
    ! thermodynamic levels k and k+1).
    do k = 2, gr%nzm
       delta_zt(k-1) = gr%dzm(1,k)
    enddo

    ! Calculate T_in_K
    T_in_K = thlm2T_in_K( gr%nzt, thlm, exner, rcm )

    ! Begin by calling Brian Griffin's implementation of the
    ! Khairoutdinov and Kogan microphysics (analytic or local formulas),
    ! the COAMPS implementation of Rutlege and Hobbes, or the Morrison
    ! microphysics.
    ! Note: COAMPS appears to have some K&K elements to it as well.

    select case ( trim( microphys_scheme ) )

    case ( "coamps" )

#ifdef COAMPS_MICRO
      call coamps_microphys_driver & 
           ( gr, runtype, time_current, dt, & ! In
             rtm, wm_zm, p_in_Pa, exner, rho, & ! In
             thlm, hydromet(:,iiri), hydromet(:,iirr),  &  ! In
             hydromet(:,iirg), hydromet(:,iirs), & ! In
             rcm, Ncm_microphys, hydromet(:,iiNr), hydromet(:,iiNi), & !In
             saturation_formula, &
             stats_metadata, &
             stats_zt, &
             Nccnm, cond, & ! Inout
             hydromet_vel_zt(:,iirs), hydromet_vel_zt(:,iiri), & ! Out
             hydromet_vel_zt(:,iirr), hydromet_vel_zt(:,iiNr),  &  ! Out
             hydromet_vel_zt(:,iirg), &  ! Out
             hydromet_mc(:,iiri), hydromet_mc(:,iirr), & ! Out
             hydromet_mc(:,iirg), hydromet_mc(:,iirs), & ! Out
             hydromet_mc(:,iiNr), & ! Out
             Ncm_mc, hydromet_mc(:,iiNi), & ! Out
             rvm_mc, rcm_mc, thlm_mc )
#else
      error stop "Not compiled with COAMPS microphysics"
      cond = -999._core_rknd
      if ( cond(1,1,1) /= cond(1,1,1) ) error stop
#endif

      if ( stats_metadata%l_stats_samp ) then

        ! Sedimentation velocity for rrm
        call stat_update_var(stats_metadata%iVrr, zt2zm( gr, hydromet_vel_zt(:,iirr) ), stats_zm)

        ! Sedimentation velocity for Nrm
        call stat_update_var(stats_metadata%iVNr, zt2zm( gr, hydromet_vel_zt(:,iiNr) ), stats_zm )

        ! Sedimentation velocity for snow
        call stat_update_var(stats_metadata%iVrs, zt2zm( gr, hydromet_vel_zt(:,iirs) ), stats_zm )

        ! Sedimentation velocity for pristine ice
        call stat_update_var( stats_metadata%iVri, zt2zm( gr, hydromet_vel_zt(:,iiri) ), stats_zm )

        ! Sedimentation velocity for graupel
        call stat_update_var( stats_metadata%iVrg, &
                            zt2zm( gr, hydromet_vel_zt(:,iirg) ), stats_zm )
      endif ! stats_metadata%l_stats_samp

    case ( "morrison" )

      if ( lh_microphys_type /= lh_microphys_disabled ) then
#ifdef SILHS
        call lh_microphys_driver( &
               gr, dt, gr%nzt, gr%nzm, lh_num_samples, & ! In
               pdf_dim, hydromet_dim, hm_metadata, & ! In
               X_nl_all_levs, lh_sample_point_weights, & ! In
               pdf_params, precip_fracs, p_in_Pa, exner, rho, & ! In
               rcm, delta_zt(1:gr%nzt), cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs,  & ! In
               lh_rt_clipped, lh_thl_clipped, & ! In
               lh_rc_clipped, lh_rv_clipped, & ! In
               lh_Nc_clipped, & ! In
               l_lh_importance_sampling, & ! In
               l_lh_instant_var_covar_src, & ! In
               saturation_formula, & ! In
               stats_metadata, & ! In
               stats_zt, stats_zm, & ! InOut
               stats_sfc, stats_lh_zt, & ! InOut
               hydromet_mc, hydromet_vel_zt, Ncm_mc, & ! Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               rtp2_mc, thlp2_mc, wprtp_mc, & ! Out
               wpthlp_mc, rtpthlp_mc, & ! Out
               lh_AKm, AKm, AKstd, AKstd_cld, &
               lh_rcm_avg, AKm_rcm, AKm_rcc, &
               morrison_microphys_driver )  ! Procedure
#else
        error stop "Latin hypercube was not enabled at compile time"
        ! Get rid of compiler warnings
        if ( .false. .and. size( X_nl_all_levs ) < 1 ) then
          rcm_mc(1) = &
            + lh_sample_point_weights(1,1) + real( X_mixt_comp_all_levs(1,1) )
        endif
#endif /* SILHS */
        call stats_accumulate_lh_tend( gr, hydromet_dim, hm_metadata, &
                                       hydromet_mc, Ncm_mc, &
                                       thlm_mc, rvm_mc, rcm_mc, &
                                       lh_AKm, AKm, AKstd, AKstd_cld, &
                                       lh_rcm_avg, AKm_rcm, AKm_rcc, &
                                       stats_metadata, &
                                       stats_lh_zt ) ! intent(inout)

      endif ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( lh_microphys_type /= lh_microphys_interactive ) then
        l_latin_hypercube_input = .false.

        if ( l_morr_xp2_mc ) then
          !Use the moister rt_1/rt_2 rather than rtm in morrison microphys
          !Also use the colder of thl_1/thl_2
          where ( pdf_params%rt_1(1,:) > pdf_params%rt_2(1,:) )
            rvm = pdf_params%rt_1(1,:) - pdf_params%rc_1(1,:)
          else where
            rvm = pdf_params%rt_2(1,:) - pdf_params%rc_2(1,:)
          end where

          where ( pdf_params%thl_1(1,:) < pdf_params%thl_2(1,:) )
            thlm_morr = pdf_params%thl_1(1,:)
          else where
            thlm_morr = pdf_params%thl_2(1,:)
          end where
        else
          rvm = rtm - rcm
          thlm_morr = thlm
        endif

        call morrison_microphys_driver( &
               gr, dt, gr%nzt, &
               hydromet_dim, hm_metadata, &
               l_latin_hypercube_input, thlm_morr, wm_zt, p_in_Pa, &
               exner, rho, cloud_frac, wtmp, &
               delta_zt, rcm, Ncm_microphys, chi, rvm, hydromet, &
               saturation_formula, &
               stats_metadata, &
               hydromet_mc, hydromet_vel_zt, Ncm_mc, &
               rcm_mc, rvm_mc, thlm_mc, &
               microphys_stats_zt, microphys_stats_sfc )

        if ( l_morr_xp2_mc) then

          rrm_evap = microphys_get_var( stats_metadata%irrm_evap, microphys_stats_zt )

          call update_xp2_mc_api( gr, gr%nzm, gr%nzt, dt, &
                                  cloud_frac, rcm, rvm, thlm_morr,            & ! Intent(in)  
                                  wm_zt, exner, rrm_evap, pdf_params,         & ! Intent(in)
                                  rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,     & ! Intent(out)
                                  rtpthlp_mc )                                  ! Intent(out)

        else

          ! Set microphysics tendencies for model variances to 0.
          rtp2_mc  = zero
          thlp2_mc = zero
          wprtp_mc = zero
          wpthlp_mc = zero
          rtpthlp_mc = zero

        endif

      endif ! lh_microphys_type /= interactive

      ! Output rain sedimentation velocity
      if ( stats_metadata%l_stats_samp ) then
        call stat_update_var(stats_metadata%iVrr, zt2zm( gr, hydromet_vel_zt(:,iirr) ), stats_zm)
      endif

    case ( "khairoutdinov_kogan" )

      if ( lh_microphys_type /= lh_microphys_disabled ) then

#ifdef SILHS
        call lh_microphys_driver( &
               gr, dt, gr%nzt, gr%nzm, lh_num_samples, & ! In
               pdf_dim, hydromet_dim, hm_metadata, & ! In
               X_nl_all_levs, lh_sample_point_weights, & ! In
               pdf_params, precip_fracs, p_in_Pa, exner, rho, & ! In
               rcm, delta_zt, cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs,  & ! In
               lh_rt_clipped, lh_thl_clipped, & ! In
               lh_rc_clipped, lh_rv_clipped, & ! In
               lh_Nc_clipped, & ! In
               l_lh_importance_sampling, & ! In
               l_lh_instant_var_covar_src, & ! In
               saturation_formula, & ! In
               stats_metadata, & ! In
               stats_zt, stats_zm, & ! InOut
               stats_sfc, stats_lh_zt, & ! InOut
               hydromet_mc, hydromet_vel_zt, Ncm_mc, & ! Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               rtp2_mc, thlp2_mc, wprtp_mc, & ! Out
               wpthlp_mc, rtpthlp_mc, & ! Out
               lh_AKm, AKm, AKstd, AKstd_cld, &
               lh_rcm_avg, AKm_rcm, AKm_rcc, &
               KK_local_microphys ) ! Procedure
#else
        error stop "Subgrid Importance Latin Hypercube was not enabled at compile time"
#endif /* SILHS */

        call stats_accumulate_lh_tend( gr, hydromet_dim, hm_metadata, &
                                       hydromet_mc, Ncm_mc, &
                                       thlm_mc, rvm_mc, rcm_mc, &
                                       lh_AKm, AKm, AKstd, AKstd_cld, &
                                       lh_rcm_avg, AKm_rcm, AKm_rcc, &
                                       stats_metadata, &
                                       stats_lh_zt ) ! intent(inout)

        if ( stats_metadata%l_stats_samp ) then
          ! Latin hypercube estimate for sedimentation velocities
          call stat_update_var( stats_metadata%ilh_Vrr, hydromet_vel_zt(:,iirr), stats_lh_zt )

          call stat_update_var( stats_metadata%ilh_VNr, hydromet_vel_zt(:,iiNr), stats_lh_zt )

        endif

      endif ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( lh_microphys_type /= lh_microphys_interactive ) then

        l_latin_hypercube_input = .false.
        rvm = rtm - rcm

        if ( l_local_kk ) then

          call KK_local_microphys( gr, dt, gr%nzt,                         & ! In
                                   hydromet_dim, hm_metadata,              & ! In
                                   l_latin_hypercube_input,                & ! In
                                   thlm, wm_zt, p_in_Pa, exner, rho,       & ! In
                                   cloud_frac, wtmp, delta_zt, rcm,        & ! In
                                   Ncm_microphys, chi, rvm, hydromet,      & ! In
                                   saturation_formula,                     & ! In
                                   stats_metadata,                         & ! In
                                   hydromet_mc, hydromet_vel_zt,           & ! Out
                                   Ncm_mc, rcm_mc, rvm_mc, thlm_mc,        & ! Out
                                   microphys_stats_zt,                     & ! Out
                                   microphys_stats_sfc                     ) ! Out

        else

          call KK_upscaled_microphys( gr, dt, gr%nzt, gr%nzm,                 & ! In
                                      pdf_dim, hydromet_dim, hm_metadata,     & ! In
                                      wm_zt, rtm, thlm, p_in_Pa,              & ! In
                                      exner, rho, rcm,                        & ! In
                                      pdf_params, hydromet_pdf_params,        & ! In
                                      precip_fracs,                           & ! In
                                      hydromet,                               & ! In
                                      mu_x_1_n, mu_x_2_n,                     & ! In
                                      sigma_x_1_n, sigma_x_2_n,               & ! In
                                      corr_array_1_n, corr_array_2_n,         & ! In
                                      saturation_formula,                     & ! In
                                      stats_metadata,                         & ! In
                                      stats_zt, stats_zm,                     & ! InOut
                                      hydromet_mc, hydromet_vel_zt,           & ! Out
                                      rcm_mc, rvm_mc, thlm_mc,                & ! Out
                                      hydromet_vel_covar_zt_impc,             & ! Out
                                      hydromet_vel_covar_zt_expc,             & ! Out
                                      wprtp_mc, wpthlp_mc, rtp2_mc,           & ! Out
                                      thlp2_mc, rtpthlp_mc                    ) ! Out

          if ( l_silhs_KK_convergence_adj_mean ) then
            ! SILHS cannot currently predict the turbulent sedimentation
            ! terms, so zero them out for a convergence test.
            hydromet_vel_covar_zt_impc = zero
            hydromet_vel_covar_zt_expc = zero
          end if

        endif

      endif ! lh_microphys_type /= interactive

      if ( stats_metadata%l_stats_samp ) then

        ! Sedimentation velocity for rrm
        call stat_update_var( stats_metadata%iVrr, zt2zm( gr, hydromet_vel_zt(:,iirr) ), stats_zm )

        ! Sedimentation velocity for Nrm
        call stat_update_var( stats_metadata%iVNr, zt2zm( gr, hydromet_vel_zt(:,iiNr) ), stats_zm )

      endif ! stats_metadata%l_stats_samp

    case ( "simplified_ice" )

      ! Call the simplified ice diffusion scheme
      call ice_dfsn( gr, dt, thlm, rcm, exner, p_in_Pa, rho, &
                     saturation_formula, &
                     stats_metadata, &
                     stats_zt, &
                     rcm_mc, thlm_mc )

    case default
      ! Do nothing
    end select ! microphys_scheme

    if ( stats_metadata%l_stats_samp ) then
      call stat_update_var( stats_metadata%iNccnm, Nccnm, stats_zt )
    endif ! stats_metadata%l_stats_samp

    ! Call GFDL activation code
    if ( l_gfdl_activation ) then
      ! Ensure a microphysics that has Ncm is being used
      if ( l_predict_Nc ) then

        ! Save the initial Ncm mc value for the Ncm_act term
        if ( stats_metadata%l_stats_samp ) then
          call stat_begin_update( gr%nzt, stats_metadata%iNcm_act, Ncm_mc, stats_zt )
        endif

        call aer_act_clubb_quadrature_Gauss( gr, pdf_params, p_in_Pa, &
                                             aeromass, T_in_K, &
                                             Ndrop_max)

        ! Convert to #/kg
        Ndrop_max = Ndrop_max * cm3_per_m3 / rho

        if( stats_metadata%l_stats_samp ) then
          call stat_update_var( stats_metadata%iNc_activated, Ndrop_max, stats_zt)
        endif

        ! Clip Ncm values that are outside of cloud by CLUBB standards
        do k = 1, gr%nzt

! ---> h1g, 2011-04-20,   no liquid drop nucleation if T < -40 C
          if ( T_in_K(k) <= 233.15_core_rknd )  Ndrop_max(k) = &
             zero  ! if T<-40C, no liquid drop nucleation
! <--- h1g, 2011-04-20

          ! Apply "clipped" Ncm to the Ncm tendency, Ncm_mc.
          if( cloud_frac(k) > cloud_frac_min ) then
            Ncm_mc(k) &
            = Ncm_mc(k) &
              + ( ( max( Ndrop_max(k), Ncm_microphys(k) ) &
                    - Ncm_microphys(k) ) &
                  / dt )
          else
            Ncm_mc(k) &
            = Ncm_mc(k) &
              - ( Ncm_microphys(k) / dt )
          endif

        enddo

        ! Update the Ncm_act term
        if( stats_metadata%l_stats_samp ) then
          call stat_end_update( gr%nzt, stats_metadata%iNcm_act, Ncm_mc, stats_zt )
        endif

      else

        error stop "Unsupported microphysics scheme for GFDL activation."

      endif  ! l_predict_Nc
    endif ! l_gfdl_activation

    ! Cloud water sedimentation.
    if ( l_cloud_sed ) then

      ! Note:  it would be very easy to upscale the cloud water sedimentation
      !        flux, so we should look into adding an upscaled option.

      call cloud_drop_sed( gr, rcm, Ncm_microphys,      & ! Intent(in)
                           rho_zm, rho, exner, sigma_g, & ! Intent(in)
                           stats_metadata,              & ! Intent(in)
                           stats_zt, stats_zm,          & ! Intent(inout)
                           rcm_mc, thlm_mc )              ! Intent(inout)

    endif ! l_cloud_sed

    ! Sample microphysics variables if necessary
    if ( microphys_stats_zt%l_allocated ) then
      call microphys_stats_accumulate( microphys_stats_zt, stats_metadata, stats_zt )
      call microphys_stats_cleanup( microphys_stats_zt )
    endif
    if ( microphys_stats_sfc%l_allocated ) then
      call microphys_stats_accumulate( microphys_stats_sfc, stats_metadata, stats_sfc )
      call microphys_stats_cleanup( microphys_stats_sfc )
    endif


    return

  end subroutine calc_microphys_scheme_tendcies

!===============================================================================

end module microphys_driver

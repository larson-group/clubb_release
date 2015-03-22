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
  public :: microphys_schemes

  private ! Default Scope

  contains

  !=============================================================================
  subroutine microphys_schemes( dt, time_current, n_variables, runtype, & ! In
                                thlm, p_in_Pa, exner, rho, rho_zm, rtm, & ! In
                                rcm, cloud_frac, wm_zt, wm_zm, wp2_zt, &  ! In
                                hydromet, Nc_in_cloud, &                  ! In
                                pdf_params, hydromet_pdf_params, &        ! In
                                X_nl_all_levs, X_mixt_comp_all_levs, &    ! In
                                lh_sample_point_weights, &                ! In
                                mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &   ! In
                                corr_array_1, corr_array_2, &             ! In
                                lh_clipped_vars, &                        ! In
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
        gr,    & ! Variable(s)
        zt2zm    ! Procedure(s)

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

    use mg_microphys_driver_module, only: &
        mg_microphys_driver  ! Procedure(s)

    use latin_hypercube_driver_module, only: &
        lh_clipped_variables_type  ! Type

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

    use advance_xp2_xpyp_module, only: &
        update_xp2_mc  ! Procedure(s)

    use parameters_model, only: & 
        hydromet_dim   ! Variable(s)

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
        hydromet_pdf_parameter  ! Type

    use stats_type_utilities, only: & 
        stat_update_var,   & ! Procedure(s)
        stat_begin_update, &
        stat_end_update

    use array_index, only:  & 
        iirrm, & ! Variable(s)
        iirsm, &
        iirim, &
        iirgm, &
        iiNrm, &
        iiNim

    use stats_variables, only: & 
        stats_zt,  & ! Variable(s)
        stats_zm,  & 
        stats_sfc, & 
        stats_lh_zt, &
        l_stats_samp

    use stats_variables, only: & 
        iVrr,  & ! Variable(s)
        iVNr, & 
        iVrs, & 
        iVri, & 
        iVrg, &
        ilh_Vrr, &
        ilh_VNr

    use stats_variables, only: & 
        iNcm_act,      & ! Variable(s)
        iNc_activated, &
        iNccnm,        &
        irrm_cond

    use stats_clubb_utilities, only: & 
        stats_accumulate_lh_tend ! Procedure(s)

    use phys_buffer, only: & ! Used for placing wp2_zt in morrison_gettelman microphysics
        pbuf_add,      &
        pbuf_allocate, &
        pbuf_setval

    use shr_kind_mod, only: &
        shr_kind_r8  ! Variable(s)

    use parameters_microphys, only: &
        lh_microphys_type,        & ! Determines how the LH samples are used
        lh_microphys_interactive, & ! Feed the subcols into microphys and allow feedback
        lh_microphys_disabled       ! Disable latin hypercube entirely

    use clubb_precision, only: &
        time_precision, & ! Variable(s)
        dp,             &
        core_rknd

    use corr_varnce_module, only: &
        d_variables ! Variable(s)

    use microphys_stats_vars_module, only: &
        microphys_stats_vars_type,  & ! Type
        microphys_stats_accumulate, & ! Procedure(s)
        microphys_get_var,          &
        microphys_stats_cleanup

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt           ! Timestep         [s]

    real( kind = time_precision ), intent(in) ::  & 
      time_current ! Current time     [s]

    integer, intent(in) :: &
      n_variables   ! Number of variables in the correlation arrays

    character(len=*), intent(in) :: & 
      runtype ! Name of the run, for case specific effects.

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      thlm,       & ! Liquid potential temp.                    [K]
      p_in_Pa,    & ! Pressure                                  [Pa]
      exner,      & ! Exner function                            [-]
      rho,        & ! Density on thermodynamic levels           [kg/m^3]
      rho_zm,     & ! Density on momentum levels                [kg/m^3]
      rtm,        & ! Total water mixing ratio                  [kg/kg]
      rcm,        & ! Liquid water mixing ratio                 [kg/kg]
      cloud_frac, & ! Cloud fraction                            [-]
      wm_zt,      & ! w wind component on thermodynamic levels  [m/s]
      wm_zm,      & ! w wind component on momentum levels       [m/s]
      wp2_zt        ! w'^2 on the thermo. grid                  [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet       ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    type(pdf_parameter), dimension(gr%nz), intent(in) :: & 
      pdf_params     ! PDF parameters

    type(hydromet_pdf_parameter), dimension(gr%nz), intent(in) :: &
      hydromet_pdf_params     ! PDF parameters

    real( kind = dp ), dimension(gr%nz,lh_num_samples,d_variables), &
    intent(in) :: &
      X_nl_all_levs ! Normally and lognormally distributed hydrometeors and other variables

    integer, dimension(gr%nz,lh_num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Which mixture component the sample is in

    real( kind = core_rknd ), dimension(lh_num_samples), intent(in) :: &
      lh_sample_point_weights ! Weights for cloud weighted sampling

    real( kind = core_rknd ), dimension(n_variables, gr%nz), intent(in) :: &
      mu_x_1,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(n_variables, n_variables, gr%nz), &
    intent(in) :: &
      corr_array_1, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    type(lh_clipped_variables_type), dimension(gr%nz,lh_num_samples) :: &
      lh_clipped_vars ! SILHS sample variables

    ! Input/Output Variables
    ! Note:
    ! For COAMPS Nccnm is initialized and Nim & Ncm are computed within
    ! subroutine adjtg.
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      Nccnm    ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(out) :: & 
      hydromet_mc     ! Change in hydrometeors due to microphysics  [units/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      Ncm_mc     ! Change in Ncm due to microphysics  [num/kg/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      rcm_mc,  & ! Microphysics contributions to liquid water          [kg/kg/s]
      rvm_mc,  & ! Microphysics contributions to vapor water           [kg/kg/s]
      thlm_mc    ! Microphysics contributions to theta-l               [K/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(out) :: &
      hydromet_vel_zt   ! Mean hydrometeor sed. velocity on thermo. levs. [m/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(out) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_xx'x_x'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      wprtp_mc,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      delta_zt  ! Difference in thermo. height levels     [m]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      T_in_K,        & ! Temperature                                [K]
      rvm,           & ! Vapor water mixing ratio                   [kg/kg]
      thlm_morr,     & ! Thlm fed into morrison microphysics        [K]
      Ncm_microphys    ! Mean cloud droplet concentration, <N_c>    [num/kg]

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      rrm_evap    ! Evaporation rate for rrm         [kg/kg/s]

    real( kind = core_rknd ), dimension(1,1,gr%nz) :: & 
      cond ! COAMPS stat for condesation/evap of rcm

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wtmp,    & ! Standard dev. of w                   [m/s]
      chi   ! The variable 's' in Mellor (1977)    [kg/kg]

    integer :: k    ! Loop index

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Ndrop_max  ! GFDL droplet activation concentration [#/kg]

    !Input aerosol mass concentration: the unit is 10^12 ug/m3.
    !For example, aeromass=2.25e-12 means that the aerosol mass
    !concentration is 2.25 ug/m3.
    !This value of aeromass was recommended by Huan Guo
    real( kind = core_rknd ), dimension(gr%nz, 4) :: &
      aeromass ! ug/m^3

    type(microphys_stats_vars_type) :: &
      microphys_stats_zt, &   ! Stats. vars. from microphys. on zt and sfc grids
      microphys_stats_sfc

    logical :: l_latin_hypercube_input


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
      stop "Runtype is null, which should not happen."
    endif

    ! Calculate the updated value of mean cloud droplet concentration, based on
    ! Nc_in_cloud and the updated value of cloud_frac.
    Ncm_microphys = Nc_in_cloud * cloud_frac

    Ndrop_max = zero

    ! Set the value of the aerosol mass array to a constant
    aeromass = aeromass_value

    ! Determine 's' from Mellor (1977)
    chi(:) = pdf_params(:)%mixt_frac * pdf_params(:)%chi_1 &
                  + ( one - pdf_params(:)%mixt_frac ) * pdf_params(:)%chi_2

    ! Compute standard deviation of vertical velocity in the grid column
    wtmp(:) = sqrt( wp2_zt(:) )

    ! Compute difference in thermodynamic height levels
    delta_zt(1:gr%nz) = one / gr%invrs_dzm(1:gr%nz)

    ! Calculate T_in_K
    T_in_K = thlm2T_in_K( thlm, exner, rcm )

    ! Begin by calling Brian Griffin's implementation of the
    ! Khairoutdinov and Kogan microphysics (analytic or local formulas),
    ! the COAMPS implementation of Rutlege and Hobbes, or the Morrison
    ! microphysics.
    ! Note: COAMPS appears to have some K&K elements to it as well.

    select case ( trim( microphys_scheme ) )

    case ( "coamps" )

#ifdef COAMPS_MICRO
      call coamps_microphys_driver & 
           ( runtype, time_current, dt, & ! In
             rtm, wm_zm, p_in_Pa, exner, rho, & ! In
             thlm, hydromet(:,iirim), hydromet(:,iirrm),  &  ! In
             hydromet(:,iirgm), hydromet(:,iirsm), & ! In
             rcm, Ncm_microphys, hydromet(:,iiNrm), hydromet(:,iiNim), & !In
             Nccnm, cond, & ! Inout
             hydromet_vel_zt(:,iirsm), hydromet_vel_zt(:,iirim), & ! Out
             hydromet_vel_zt(:,iirrm), hydromet_vel_zt(:,iiNrm),  &  ! Out
             hydromet_vel_zt(:,iirgm), &  ! Out
             hydromet_mc(:,iirim), hydromet_mc(:,iirrm), & ! Out
             hydromet_mc(:,iirgm), hydromet_mc(:,iirsm), & ! Out
             hydromet_mc(:,iiNrm), & ! Out
             Ncm_mc, hydromet_mc(:,iiNim), & ! Out
             rvm_mc, rcm_mc, thlm_mc )
#else
      stop "Not compiled with COAMPS microphysics"
      cond = -999._core_rknd
      if ( cond(1,1,1) /= cond(1,1,1) ) stop
#endif

      if ( l_stats_samp ) then

        ! Sedimentation velocity for rrm
        call stat_update_var(iVrr, zt2zm( hydromet_vel_zt(:,iirrm) ), stats_zm)

        ! Sedimentation velocity for Nrm
        call stat_update_var(iVNr, zt2zm( hydromet_vel_zt(:,iiNrm) ), stats_zm )

        ! Sedimentation velocity for snow
        call stat_update_var(iVrs, zt2zm( hydromet_vel_zt(:,iirsm) ), stats_zm )

        ! Sedimentation velocity for pristine ice
        call stat_update_var( iVri, zt2zm( hydromet_vel_zt(:,iirim) ), stats_zm )

        ! Sedimentation velocity for graupel
        call stat_update_var( iVrg, &
                            zt2zm( hydromet_vel_zt(:,iirgm) ), stats_zm )
      endif ! l_stats_samp

    case ( "morrison" )

      if ( lh_microphys_type /= lh_microphys_disabled ) then
#ifdef SILHS
        call lh_microphys_driver &
             ( dt, gr%nz, lh_num_samples, d_variables, & ! In
               X_nl_all_levs, lh_sample_point_weights, & ! In
               pdf_params, hydromet_pdf_params, p_in_Pa, exner, rho, & ! In
               rcm, delta_zt, cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs,  & !In
               lh_clipped_vars, & ! In
               hydromet_mc, hydromet_vel_zt, Ncm_mc, & ! Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               rtp2_mc, thlp2_mc, wprtp_mc, & ! Out
               wpthlp_mc, rtpthlp_mc, & ! Out
               morrison_microphys_driver )  ! Procedure
#else
        stop "Latin hypercube was not enabled at compile time"
        ! Get rid of compiler warnings
        if ( .false. .and. size( X_nl_all_levs ) < 1 ) then
          rcm_mc(1) = &
            + lh_sample_point_weights(1) + real( X_mixt_comp_all_levs(1,1) )
        endif
#endif /* SILHS */
        call stats_accumulate_lh_tend( hydromet_mc, Ncm_mc, &
                                     thlm_mc, rvm_mc, rcm_mc )

      endif ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( lh_microphys_type /= lh_microphys_interactive ) then
        l_latin_hypercube_input = .false.

        if ( l_morr_xp2_mc ) then
          !Use the moister rt_1/rt_2 rather than rtm in morrison microphys
          !Also use the colder of thl_1/thl_2
          where ( pdf_params%rt_1 > pdf_params%rt_2 )
            rvm = pdf_params%rt_1 - pdf_params%rc_1
          else where
            rvm = pdf_params%rt_2 - pdf_params%rc_2
          end where

          where ( pdf_params%thl_1 < pdf_params%thl_2 )
            thlm_morr = pdf_params%thl_1
          else where
            thlm_morr = pdf_params%thl_2
          end where
        else
          rvm = rtm - rcm
          thlm_morr = thlm
        endif

        call morrison_microphys_driver &
             ( dt, gr%nz, &
               l_latin_hypercube_input, thlm_morr, wm_zt, p_in_Pa, &
               exner, rho, cloud_frac, wtmp, &
               delta_zt, rcm, Ncm_microphys, chi, rvm, hydromet, &
               hydromet_mc, hydromet_vel_zt, Ncm_mc, &
               rcm_mc, rvm_mc, thlm_mc, &
               microphys_stats_zt, microphys_stats_sfc )

        if ( l_morr_xp2_mc) then

          rrm_evap = microphys_get_var( irrm_cond, microphys_stats_zt )

          call update_xp2_mc( gr%nz, dt, cloud_frac, rcm, rvm, thlm_morr, & ! Intent(in)  
                              wm_zt, exner, rrm_evap, pdf_params,      & ! Intent(in)
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


        ! Output rain sedimentation velocity
        if ( l_stats_samp ) then
          call stat_update_var(iVrr, zt2zm( hydromet_vel_zt(:,iirrm) ), stats_zm)
        endif

      endif ! lh_microphys_type /= interactive

    case ( "morrison_gettelman" )

      ! Place wp2 into the dummy phys_buffer module to import it into microp_aero_ts.
      ! Placed here because parameters cannot be changed on mg_microphys_driver with
      ! the way LH is currently set up.
      call pbuf_add( 'WP2', 1, gr%nz, 1 )
      call pbuf_allocate()
      call pbuf_setval( 'WP2', real( wp2_zt, kind=shr_kind_r8 ) )

      rvm = rtm - rcm
      call mg_microphys_driver &
           ( dt, gr%nz, l_stats_samp, gr%invrs_dzt, thlm, p_in_Pa, exner, &
             rho, cloud_frac, rcm, Ncm_microphys, rvm, Nccnm, pdf_params, hydromet, &
             hydromet_mc, hydromet_vel_zt, rcm_mc, rvm_mc, thlm_mc )

    case ( "khairoutdinov_kogan" )

      if ( lh_microphys_type /= lh_microphys_disabled ) then

#ifdef SILHS
        call lh_microphys_driver &
             ( dt, gr%nz, lh_num_samples, d_variables, & ! In
               X_nl_all_levs, lh_sample_point_weights, & ! In
               pdf_params, hydromet_pdf_params, p_in_Pa, exner, rho, & ! In
               rcm, delta_zt, cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs,  & !In
               lh_clipped_vars, & ! In
               hydromet_mc, hydromet_vel_zt, Ncm_mc, & ! Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               rtp2_mc, thlp2_mc, wprtp_mc, & ! Out
               wpthlp_mc, rtpthlp_mc, & ! Out               
               KK_local_microphys ) ! Procedure
#else
        stop "Subgrid Importance Latin Hypercube was not enabled at compile time"
#endif /* SILHS */

        call stats_accumulate_lh_tend( hydromet_mc, Ncm_mc, &
                                       thlm_mc, rvm_mc, rcm_mc )

        if ( l_stats_samp ) then
          ! Latin hypercube estimate for sedimentation velocities
          call stat_update_var( ilh_Vrr, hydromet_vel_zt(:,iirrm), stats_lh_zt )

          call stat_update_var( ilh_VNr, hydromet_vel_zt(:,iiNrm), stats_lh_zt )

        endif

      endif ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( lh_microphys_type /= lh_microphys_interactive ) then

        l_latin_hypercube_input = .false.
        rvm = rtm - rcm

        if ( l_local_kk ) then

          call KK_local_microphys( dt, gr%nz, l_latin_hypercube_input,     & ! In
                                   thlm, wm_zt, p_in_Pa, exner, rho,       & ! In
                                   cloud_frac, wtmp, delta_zt, rcm,        & ! In
                                   Ncm_microphys, chi, rvm, hydromet, & ! In
                                   hydromet_mc, hydromet_vel_zt,           & ! Out
                                   Ncm_mc, rcm_mc, rvm_mc, thlm_mc,        & ! Out
                                   microphys_stats_zt,                     & ! Out
                                   microphys_stats_sfc                     ) ! Out

        else

          call KK_upscaled_microphys( dt, gr%nz, n_variables, l_stats_samp, & ! In
                                      wm_zt, rtm, thlm, p_in_Pa,            & ! In
                                      exner, rho, rcm, Nc_in_cloud,         & ! In
                                      pdf_params, hydromet_pdf_params,      & ! In
                                      hydromet,                             & ! In
                                      mu_x_1, mu_x_2,                       & ! In
                                      sigma_x_1, sigma_x_2,                 & ! In
                                      corr_array_1, corr_array_2,           & ! In
                                      hydromet_mc, hydromet_vel_zt,         & ! Out
                                      rcm_mc, rvm_mc, thlm_mc,              & ! Out
                                      hydromet_vel_covar_zt_impc,           & ! Out
                                      hydromet_vel_covar_zt_expc,           & ! Out
                                      wprtp_mc, wpthlp_mc, rtp2_mc,         & ! Out
                                      thlp2_mc, rtpthlp_mc                  ) ! Out

          if ( l_silhs_KK_convergence_adj_mean ) then
            ! SILHS cannot currently predict the turbulent sedimentation
            ! terms, so zero them out for a convergence test.
            hydromet_vel_covar_zt_impc = zero
            hydromet_vel_covar_zt_expc = zero
          end if

        endif

      endif ! lh_microphys_type /= interactive

      if ( l_stats_samp ) then

        ! Sedimentation velocity for rrm
        call stat_update_var( iVrr, zt2zm( hydromet_vel_zt(:,iirrm) ), stats_zm )

        ! Sedimentation velocity for Nrm
        call stat_update_var( iVNr, zt2zm( hydromet_vel_zt(:,iiNrm) ), stats_zm )

      endif ! l_stats_samp

    case ( "simplified_ice" )

      ! Call the simplified ice diffusion scheme
      call ice_dfsn( dt, thlm, rcm, exner, p_in_Pa, rho, rcm_mc, thlm_mc )

    case default
      ! Do nothing
    end select ! microphys_scheme

    if ( l_stats_samp ) then
      call stat_update_var( iNccnm, Nccnm, stats_zt )
    endif ! l_stats_samp

    ! Call GFDL activation code
    if ( l_gfdl_activation ) then
      ! Ensure a microphysics that has Ncm is being used
      if ( l_predict_Nc ) then

        ! Save the initial Ncm mc value for the Ncm_act term
        if ( l_stats_samp ) then
          call stat_begin_update( iNcm_act, Ncm_mc, stats_zt )
        endif

        call aer_act_clubb_quadrature_Gauss( aeromass, T_in_K, Ndrop_max )

        ! Convert to #/kg
        Ndrop_max = Ndrop_max * cm3_per_m3 / rho

        if( l_stats_samp ) then
          call stat_update_var( iNc_activated, Ndrop_max, stats_zt)
        endif

        ! Clip Ncm values that are outside of cloud by CLUBB standards
        do k = 1, gr%nz

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
        if( l_stats_samp ) then
          call stat_end_update( iNcm_act, Ncm_mc, stats_zt )
        endif

      else

        stop "Unsupported microphysics scheme for GFDL activation."

      endif  ! l_predict_Nc
    endif ! l_gfdl_activation

    ! Cloud water sedimentation.
    if ( l_cloud_sed ) then

      ! Note:  it would be very easy to upscale the cloud water sedimentation
      !        flux, so we should look into adding an upscaled option.

      call cloud_drop_sed( rcm, Ncm_microphys,          & ! Intent(in)
                           rho_zm, rho, exner, sigma_g, & ! Intent(in)
                           rcm_mc, thlm_mc )              ! Intent(inout)

    endif ! l_cloud_sed

    ! Sample microphysics variables if necessary
    if ( microphys_stats_zt%l_allocated ) then
      call microphys_stats_accumulate( microphys_stats_zt, l_stats_samp, stats_zt )
      call microphys_stats_cleanup( microphys_stats_zt )
    endif
    if ( microphys_stats_sfc%l_allocated ) then
      call microphys_stats_accumulate( microphys_stats_sfc, l_stats_samp, stats_sfc )
      call microphys_stats_cleanup( microphys_stats_sfc )
    endif


    return

  end subroutine microphys_schemes

!===============================================================================

end module microphys_driver

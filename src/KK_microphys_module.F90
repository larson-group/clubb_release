! $Id$
!===============================================================================
module KK_microphys_module

  implicit none

  private

  public :: KK_local_micro_driver, &
            KK_upscaled_micro_driver, &
            KK_in_precip_values, &
            KK_upscaled_setup

  private :: KK_micro_init, &
             component_means_rain, &
             precip_fraction, &
             KK_tendency_coefs, &
             KK_upscaled_means_driver, &
             KK_upscaled_covar_driver, &
             KK_microphys_adjust, &
             KK_upscaled_stats, &
             KK_stats_output, &
             KK_sedimentation, &
             variance_KK_mvr

  contains

  !=============================================================================
  subroutine KK_local_micro_driver( dt, nz, l_stats_samp, &
                                    l_latin_hypercube, thlm, wm_zt, &
                                    p_in_Pa, exner, rho, cloud_frac, &
                                    pdf_params, w_std_dev, dzq, rcm, &
                                    Ncm, s_mellor, rvm, Nc0_in_cloud, &
                                    hydromet, &
                                    hydromet_mc, hydromet_vel, &
                                    rcm_mc, rvm_mc, thlm_mc, &
                                    wprtp_mc_tndcy, wpthlp_mc_tndcy, &
                                    rtp2_mc_tndcy, thlp2_mc_tndcy, &
                                    rtpthlp_mc_tndcy, &
                                    KK_auto_tndcy, KK_accr_tndcy )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr

    use constants_clubb, only: &
        one,    & ! Constant(s)
        zero,   &
        rho_lw, &
        rr_tol, &
        Nr_tol, &
        Nc_tol

    use KK_local_means, only: &
        KK_mvr_local_mean,  & ! Procedure(s)
        KK_evap_local_mean, &
        KK_auto_local_mean, &
        KK_accr_local_mean

    use KK_Nrm_tendencies, only: &
        KK_Nrm_evap_local_mean, & ! Procedure(s)
        KK_Nrm_auto_mean

    use KK_utilities, only: &
        get_cloud_top_level    ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    implicit none

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    logical, intent(in) :: &
      l_stats_samp,      & ! Flag to sample statistics
      l_latin_hypercube    ! Flag to use Latin Hypercube interface

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      thlm,       & ! Mean liquid water potential temperature         [K]
      wm_zt,      & ! Mean vertical velocity on thermodynamic levels  [m/s]
      p_in_Pa,    & ! Pressure                                        [Pa]
      exner,      & ! Exner function                                  [-]
      rho,        & ! Density                                         [kg/m^3]
      cloud_frac, & ! Cloud fraction                                  [-]
      rcm,        & ! Mean cloud water mixing ratio                   [kg/kg]
      Ncm,        & ! Mean cloud droplet conc., < N_c >               [num/kg]
      s_mellor      ! Mean extended liquid water mixing ratio         [kg/kg]

    type(pdf_parameter), dimension(nz), target, intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      w_std_dev, & ! Standard deviation of w (for LH interface)          [m/s]
      dzq,       & ! Thickness between thermo. levels (for LH interface) [m]
      rvm          ! Mean water vapor mixing ratio (for LH interface)    [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Nc0_in_cloud    ! Constant in-cloud value of cloud droplet conc.  [num/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet    ! Hydrometeor species                      [units vary]

    ! Input / Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(inout) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      wprtp_mc_tndcy,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc_tndcy,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc_tndcy,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc_tndcy, & ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]
      KK_auto_tndcy,    & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy       ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(:), pointer ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
      Vrr,             & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,             & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) :: &
      KK_evap_tndcy,   & ! Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
      KK_mean_vol_rad    ! Mean KK rain drop mean volume radius     [m]

    real( kind = core_rknd ), dimension(nz) :: &
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    real( kind = core_rknd ) ::  &
      rrainm_source,     & ! Total source term rate for rrainm       [(kg/kg)/s]
      Nrm_source           ! Total source term rate for Nrm         [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) ::  &
      rrainm_src_adj,  & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj,     & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrainm_evap_net, & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net       ! Net evaporation rate of <N_r>            [(num/kg)/s]

    logical :: &
      l_src_adj_enabled,  & ! Flag to enable rrainm/Nrm source adjustment
      l_evap_adj_enabled, & ! Flag to enable rrainm/Nrm evaporation adjustment
      l_stats_samp_in_sub   ! Used to disable stats when SILHS is enabled

    integer :: &
      cloud_top_level, & ! Vertical level index of cloud top 
      k                  ! Loop index


    ! Remove compiler warnings
    if ( .false. ) then
      rrainm_src_adj = dzq
      rrainm_src_adj = rvm
      rrainm_src_adj = w_std_dev
      rrainm_src_adj = cloud_frac
    end if

    ! Disable stats when latin hypercube is enabled
    if ( .not. l_latin_hypercube .and. l_stats_samp ) then
     l_stats_samp_in_sub = .true.
    else
     l_stats_samp_in_sub = .false.
    end if

    !!! Initialize microphysics fields.
    call KK_micro_init( nz, hydromet, hydromet_mc, hydromet_vel, &
                        rrainm, Nrm, Vrr, VNr, &
                        rrainm_mc_tndcy, Nrm_mc_tndcy, &
                        KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                        KK_mean_vol_rad, KK_Nrm_evap_tndcy, &
                        KK_Nrm_auto_tndcy, &
                        l_src_adj_enabled, l_evap_adj_enabled )

    !!! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

      !!! Calculate the coefficients for the KK microphysics tendencies.
      call KK_tendency_coefs( thlm(k), exner(k), p_in_Pa(k), rho(k), &
                              KK_evap_coef, KK_auto_coef, &
                              KK_accr_coef, KK_mvr_coef )


      !!! Calculate the local KK rain drop mean volume radius.
      if ( rrainm(k) > rr_tol ) then

         KK_mean_vol_rad(k)  &
         = KK_mvr_local_mean( rrainm(k), Nrm(k), KK_mvr_coef )

      else  ! r_r or N_r = 0.

         KK_mean_vol_rad(k) = zero

      endif

      !!! Calculate the values of the local KK microphysics tendencies.

      !!! Calculate the local KK evaporation tendency.
      if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

         KK_evap_tndcy(k)  &
         = KK_evap_local_mean( s_mellor(k), rrainm(k), Nrm(k), &
                               KK_evap_coef )

      else  ! r_r or N_r = 0.

         KK_evap_tndcy(k) = zero

      endif

      !!! Calculate the local KK autoconversion tendency.
      if ( Ncm(k) > Nc_tol ) then

         KK_auto_tndcy(k)  &
         = KK_auto_local_mean( s_mellor(k), Ncm(k), KK_auto_coef )

      else  ! N_c = 0.

         KK_auto_tndcy(k) = zero

      endif

      !!! Calculate the local KK accretion tendency.
      if ( rrainm(k) > rr_tol ) then

         KK_accr_tndcy(k)  &
         = KK_accr_local_mean( s_mellor(k), rrainm(k), KK_accr_coef )

      else  ! r_r = 0.

         KK_accr_tndcy(k) = zero

      endif


      !!! KK rain drop concentration microphysics tendencies.

      !!! Calculate the KK N_r evaporation tendency.
      if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

         KK_Nrm_evap_tndcy(k)  &
         = KK_Nrm_evap_local_mean( KK_evap_tndcy(k), Nrm(k), rrainm(k), dt )

      else  ! r_r or N_r = 0.

          KK_Nrm_evap_tndcy(k) = zero

      endif

      !!! Calculate the KK N_r autoconversion tendency.
      KK_Nrm_auto_tndcy(k) = KK_Nrm_auto_mean( KK_auto_tndcy(k) )


      !!! Calculate any necessary adjustments to KK microphysics tendencies.
      call KK_microphys_adjust( dt, exner(k), rcm(k), rrainm(k), Nrm(k), &
                                KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                KK_accr_tndcy(k), KK_Nrm_evap_tndcy(k), &
                                KK_Nrm_auto_tndcy(k), l_src_adj_enabled, &
                                l_evap_adj_enabled, l_stats_samp_in_sub, k, &
                                rrainm_mc_tndcy(k), Nrm_mc_tndcy(k), &
                                rvm_mc(k), rcm_mc(k), thlm_mc(k) )

      !!! Statistical output for mean microphysics tendenices.
      call KK_stats_output( KK_evap_tndcy(k), KK_auto_tndcy(k), &
                            KK_accr_tndcy(k), KK_mean_vol_rad(k), &
                            KK_Nrm_evap_tndcy(k), KK_Nrm_auto_tndcy(k), &
                            l_stats_samp_in_sub, k )


    enddo  ! Microphysics tendency loop: k = 2, nz, 1


    ! Microphysics tendency terms for model variances and covariances
    ! are set to 0.
    wprtp_mc_tndcy   = zero
    wpthlp_mc_tndcy  = zero
    rtp2_mc_tndcy    = zero
    thlp2_mc_tndcy   = zero
    rtpthlp_mc_tndcy = zero

    !!! Boundary conditions for microphysics tendencies.

    ! Explicit contributions to rrainm and Nrm from microphysics are not set at
    ! thermodynamic level k = 1 because it is below the model lower boundary.
    rrainm_mc_tndcy(1) = zero
    Nrm_mc_tndcy(1)    = zero

    rrainm_mc_tndcy(nz) = zero
    Nrm_mc_tndcy(nz)    = zero

    ! Boundary conditions
    KK_mean_vol_rad(1)  = zero
    KK_mean_vol_rad(nz) = zero

    rvm_mc(1)  = zero
    rvm_mc(nz) = zero

    rcm_mc(1)  = zero
    rcm_mc(nz) = zero

    thlm_mc(1)  = zero
    thlm_mc(nz) = zero

    ! Find the vertical level index of cloud top.
    cloud_top_level = get_cloud_top_level( nz, rcm )

    !!! Microphysics sedimentation velocities.
    call KK_sedimentation( nz, cloud_top_level, KK_mean_vol_rad, Vrr, VNr )


    return

  end subroutine KK_local_micro_driver

  !=============================================================================
  subroutine KK_upscaled_micro_driver( dt, nz, l_stats_samp, thlm, wm_zt, &
                                       p_in_Pa, exner, rho, cloud_frac, &
                                       pdf_params, w_std_dev, rcm, Ncm, &
                                       s_mellor, Nc0_in_cloud, &
                                       hydromet, wphydrometp, &
                                       hydromet_mc, hydromet_vel, &
                                       rcm_mc, rvm_mc, thlm_mc, &
                                       hydromet_vel_covar_zt_impc, &
                                       hydromet_vel_covar_zt_expc, &
                                       wprtp_mc_tndcy, wpthlp_mc_tndcy, &
                                       rtp2_mc_tndcy, thlp2_mc_tndcy, &
                                       rtpthlp_mc_tndcy, &
                                       KK_auto_tndcy, KK_accr_tndcy )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm, &  ! Procedure(s)
        zm2zt, &
        gr

    use constants_clubb, only: &
        one,    & ! Constant(s)
        zero,   &
        rr_tol, &
        Nr_tol, &
        Nc_tol, &
        eps

    use parameters_microphys, only: &
        l_var_covar_src,     & ! Flag for using variance/covariance src terms
        l_const_Nc_in_cloud    ! Flag to use a const. value of N_c within cloud

    use KK_Nrm_tendencies, only: &
        KK_Nrm_evap_upscaled_mean, & ! Procedure(s)
        KK_Nrm_auto_mean

    use KK_upscaled_turbulent_sed, only: &
        KK_sed_vel_covars  ! Procedure(s)

    use KK_utilities, only: &
        get_cloud_top_level    ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        iirrainm, & ! Constant(s)
        iiNrm

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        irr1,             & ! Variable(s)
        irr2,             &
        iNr1,             &
        iNr2,             &
        iprecip_frac,     &
        iprecip_frac_1,   &
        iprecip_frac_2,   &
        iVrrprrp_expcalc, &
        iVNrpNrp_expcalc, &
        zt,               &
        zm

    use model_flags, only: &
        l_use_precip_frac, & ! Flag(s)
        l_calc_w_corr

    use advance_windm_edsclrm_module, only: &
        xpwp_fnc

    use variables_diagnostic_module, only: &
        Kh_zm

    use parameters_tunable, only: &
        c_Krrainm

    implicit none

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      thlm,       & ! Mean liquid water potential temperature         [K]
      wm_zt,      & ! Mean vertical velocity on thermodynamic levels  [m/s]
      p_in_Pa,    & ! Pressure                                        [Pa]
      exner,      & ! Exner function                                  [-]
      rho,        & ! Density                                         [kg/m^3]
      cloud_frac, & ! Cloud fraction                                  [-]
      rcm,        & ! Mean cloud water mixing ratio                   [kg/kg]
      Ncm,        & ! Mean cloud droplet conc., < N_c >               [num/kg]
      s_mellor      ! Mean extended liquid water mixing ratio         [kg/kg]

    type(pdf_parameter), dimension(nz), target, intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      w_std_dev    ! Standard deviation of w (for LH interface)          [m/s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Nc0_in_cloud    ! Constant in-cloud value of cloud droplet conc.  [num/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet,    & ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)      [(m/s)units]

    ! Input / Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(inout) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(out) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      wprtp_mc_tndcy,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc_tndcy,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc_tndcy,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc_tndcy, & ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]
      KK_auto_tndcy,    & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy       ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(:), pointer ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r > [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r > [num/kg]
      wprrp,           & ! Covariance of w and r_r, < w'r_r' >   [(m/s)(kg/kg)]
      wpNrp,           & ! Covariance of w and N_r, < w'N_r' >   [(m/s)(num/kg)]
      Vrr,             & ! Mean sedimentation velocity of r_r    [m/s]
      VNr,             & ! Mean sedimentation velocity of N_r    [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics    [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics    [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) :: &
      KK_evap_tndcy,   & ! Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
      KK_mean_vol_rad    ! Mean KK rain drop mean volume radius     [m]

    real( kind = core_rknd ), dimension(nz) :: &
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    real( kind = core_rknd ), dimension(nz) :: &
      rc1,         & ! Mean of r_c (1st PDF component)              [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)              [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)           [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)           [-]
      mixt_frac      ! Mixture fraction                             [-]

    real( kind = core_rknd ), dimension(nz) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    real( kind = core_rknd ) :: &
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Nc_1,       & ! Mean of Nc (1st PDF component)                 [num/kg]
      mu_Nc_2,       & ! Mean of Nc (2nd PDF component)                 [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Nc_1_n,     & ! Mean of ln Nc (1st PDF component)          [ln(num/kg)]
      mu_Nc_2_n,     & ! Mean of ln Nc (2nd PDF component)          [ln(num/kg)]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Nc_1,    & ! Standard deviation of Nc (1st PDF component)   [num/kg]
      sigma_Nc_2,    & ! Standard deviation of Nc (2nd PDF component)   [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Nc_1_n,  & ! Standard dev. of ln Nc (1st PDF comp.)     [ln(num/kg)]
      sigma_Nc_2_n,  & ! Standard dev. of ln Nc (2nd PDF comp.)     [ln(num/kg)]
      corr_srr_1,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_srr_2,    & ! Correlation between s and rr (2nd PDF component) ip [-]
      corr_sNr_1,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNr_2,    & ! Correlation between s and Nr (2nd PDF component) ip [-]
      corr_sNc_1,    & ! Correlation between s and Nc (1st PDF component)    [-]
      corr_sNc_2,    & ! Correlation between s and Nc (2nd PDF component)    [-]
      corr_rrNr_1,   & ! Correlation between rr & Nr (1st PDF component) ip  [-]
      corr_rrNr_2,   & ! Correlation between rr & Nr (2nd PDF component) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF comp.)     [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF comp.)     [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      corr_sw,       & ! Correlation between s & w (both components)         [-]
      corr_wrr,      & ! Correlation between rr & w (both components)        [-]
      corr_wNr,      & ! Correlation between Nr & w (both components)        [-]
      corr_wNc         ! Correlation between Nc & w (both components)        [-]

    real( kind = core_rknd ), dimension(:), pointer :: &
      Vrrprrp_zt_impc, & ! Imp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)]
      Vrrprrp_zt_expc, & ! Exp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)(kg/kg)]
      VNrpNrp_zt_impc, & ! Imp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)]
      VNrpNrp_zt_expc    ! Exp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz) :: &
      wprtp_mc_tndcy_zt,   & ! Micro. tend. for <w'rt'>; t-lev   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy_zt,  & ! Micro. tend. for <w'thl'>; t-lev  [m*K/s^2]
      rtp2_mc_tndcy_zt,    & ! Micro. tend. for <rt'^2>; t-lev   [(kg/kg)^2/s]
      thlp2_mc_tndcy_zt,   & ! Micro. tend. for <thl'^2>; t-lev  [K^2/s]
      rtpthlp_mc_tndcy_zt    ! Micro. tend. for <rt'thl'>; t-lev [K*(kg/kg)/s]

    real( kind = core_rknd ) ::  &
      rrainm_source,     & ! Total source term rate for rrainm       [(kg/kg)/s]
      Nrm_source           ! Total source term rate for Nrm         [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) ::  &
      rrainm_src_adj,  & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj,     & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrainm_evap_net, & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net       ! Net evaporation rate of <N_r>            [(num/kg)/s]

    ! changes by janhft 10/04/12
    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zm,     & ! Covariance of s and w (momentum levels)   [(m/s)(kg/kg)]
      wpNcp_zm,    & ! Covariance of N_c and w (momentum levels) [(m/s)(num/kg)]
      wpsp_zt,     & ! Covariance of s and w on t-levs           [(m/s)(kg/kg)]
      wprrp_ip_zt, & ! Covar. of r_r and w (in-precip) on t-levs [(m/s)(kg/kg)]
      wpNrp_ip_zt, & ! Covar. of N_r and w (in-precip) on t-levs [(m/s)(num/kg)]
      wpNcp_zt       ! Covariance of N_c and w on t-levs         [(m/s)(num/kg)]
    ! end changes by janhft 10/04/12

    logical :: &
      l_src_adj_enabled,  & ! Flag to enable rrainm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrainm/Nrm evaporation adjustment

    integer :: &
      cloud_top_level, & ! Vertical level index of cloud top 
      k                  ! Loop index


    !!! Initialize microphysics fields.
    call KK_micro_init( nz, hydromet, hydromet_mc, hydromet_vel, &
                        rrainm, Nrm, Vrr, VNr, &
                        rrainm_mc_tndcy, Nrm_mc_tndcy, &
                        KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                        KK_mean_vol_rad, KK_Nrm_evap_tndcy, &
                        KK_Nrm_auto_tndcy, &
                        l_src_adj_enabled, l_evap_adj_enabled )

    ! Covariance of vertical velocity and a hydrometeor
    ! (< w'r_r' > and < w'N_r' >).
    wprrp => wphydrometp(:,iirrainm)
    wpNrp => wphydrometp(:,iiNrm)

    ! The implicit and explicit components used to calculate the covariances of
    ! hydrometeor sedimentation velocities and their associated hydrometeors
    ! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp_zt_impc => hydromet_vel_covar_zt_impc(:,iirrainm)
    Vrrprrp_zt_expc => hydromet_vel_covar_zt_expc(:,iirrainm)
    VNrpNrp_zt_impc => hydromet_vel_covar_zt_impc(:,iiNrm)
    VNrpNrp_zt_expc => hydromet_vel_covar_zt_expc(:,iiNrm)

    ! Setup some of the PDF parameters
    rc1         = pdf_params%rc1
    rc2         = pdf_params%rc2
    cloud_frac1 = pdf_params%cloud_frac1
    cloud_frac2 = pdf_params%cloud_frac2
    mixt_frac   = pdf_params%mixt_frac

    ! Precipitation fraction
    if ( l_use_precip_frac ) then

       call component_means_rain( nz, rrainm, Nrm, rho, rc1, rc2, &
                                  mixt_frac, l_stats_samp, &
                                  rr1, rr2, Nr1, Nr2 )

       call precip_fraction( nz, rrainm, rr1, rr2, Nrm, Nr1, Nr2, &
                             cloud_frac, cloud_frac1, mixt_frac, &
                             precip_frac, precip_frac_1, precip_frac_2 )

    else

       rr1 = rrainm
       rr2 = rrainm
       Nr1 = Nrm
       Nr2 = Nrm

       precip_frac   = one
       precip_frac_1 = one
       precip_frac_2 = one

    endif

    ! Statistics
    if ( l_stats_samp ) then

       if ( irr1 > 0 ) then
          ! Mean rain water mixing ratio in PDF component 1.
          call stat_update_var( irr1, rr1, zt )
       endif

       if ( irr2 > 0 ) then
          ! Mean rain water mixing ratio in PDF component 2.
          call stat_update_var( irr2, rr2, zt )
       endif

       if ( iNr1 > 0 ) then
          ! Mean rain drop concentration in PDF component 1.
          call stat_update_var( iNr1, Nr1, zt )
       endif

       if ( iNr2 > 0 ) then
          ! Mean rain drop concentration in PDF component 2.
          call stat_update_var( iNr2, Nr2, zt )
       endif

       if ( iprecip_frac > 0 ) then
          ! Overall precipitation fraction.
          call stat_update_var( iprecip_frac, precip_frac, zt )
       endif

       if ( iprecip_frac_1 > 0 ) then
          ! Precipitation fraction in PDF component 1.
          call stat_update_var( iprecip_frac_1, precip_frac_1, zt )
       endif

       if ( iprecip_frac_2 > 0 ) then
          ! Precipitation fraction in PDF component 2.
          call stat_update_var( iprecip_frac_2, precip_frac_2, zt )
       endif

    endif

    ! calculate the covariances of w with the hydrometeors
    if ( l_calc_w_corr ) then

       ! calculate the covariances of w with the hydrometeors
       do k = 1, nz
          wpsp_zm(k) = pdf_params(k)%mixt_frac &
                       * ( one - pdf_params(k)%mixt_frac ) &
                       * ( pdf_params(k)%s1 - pdf_params(k)%s2 ) &
                       * ( pdf_params(k)%w1 - pdf_params(k)%w2 )
       enddo

       wpNcp_zm(1:nz-1) = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), Ncm(1:nz-1), &
                                    Ncm(2:nz), gr%invrs_dzm(1:nz-1) )

       ! Boundary conditions; We are assuming zero flux at the top.
       wpNcp_zm(nz) = zero

       ! interpolate back to zt-grid
       wpsp_zt     = zm2zt(wpsp_zm)
       wprrp_ip_zt = zm2zt(wprrp) / max( precip_frac, eps )
       wpNrp_ip_zt = zm2zt(wpNrp) / max( precip_frac, eps )
       wpNcp_zt    = zm2zt(wpNcp_zm)

       do k = 1, nz, 1
          if ( rrainm(k) <= rr_tol ) then
             wprrp_ip_zt(k) = zero
          endif
          if ( Nrm(k) <= Nr_tol ) then
             wpNrp_ip_zt(k) = zero
          endif
          if ( Ncm(k) <= Nc_tol ) then
             wpNcp_zt(k) = zero
          endif
       enddo

    endif


    !!! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

      !!! Calculate the coefficients for the KK microphysics tendencies.
      call KK_tendency_coefs( thlm(k), exner(k), p_in_Pa(k), rho(k), &
                              KK_evap_coef, KK_auto_coef, &
                              KK_accr_coef, KK_mvr_coef )

      !!! Calculate the means, standard deviations, and correlations involving
      !!! rain water mixing ratio and rain drop concentration for each PDF
      !!! component.
      call KK_in_precip_values( rr1(k), rr2(k), Nr1(k), Nr2(k), rc1(k), &
                                rc2(k), cloud_frac1(k), cloud_frac2(k), &
                                precip_frac_1(k), precip_frac_2(k), &
                                mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                sigma_Nr_2, corr_srr_1, corr_srr_2, &
                                corr_sNr_1, corr_sNr_2, corr_rrNr_1, &
                                corr_rrNr_2 )

      !!! Calculate the mean, standard deviations, and correlations involving
      !!! ln r_r, ln N_r, and ln N_c for each PDF component.
      call KK_upscaled_setup( rcm(k), rrainm(k), Nrm(k), Ncm(k), &
                              rr1(k), rr2(k), Nr1(k), Nr2(k), &
                              mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                              sigma_rr_1, sigma_rr_2, &
                              sigma_Nr_1, sigma_Nr_2, &
                              wpsp_zt(k), wprrp_ip_zt(k), wpNrp_ip_zt(k), &
                              wpNcp_zt(k), w_std_dev(k), mixt_frac(k), &
                              pdf_params(k), &
                              corr_srr_1, corr_srr_2, corr_sNr_1, &
                              corr_sNr_2, corr_rrNr_1, corr_rrNr_2, &
                              mu_s_1, mu_s_2, mu_Nc_1, mu_Nc_2, &
                              mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                              mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                              sigma_s_1, sigma_s_2, &
                              sigma_Nc_1, sigma_Nc_2, &
                              sigma_rr_1_n, sigma_rr_2_n, &
                              sigma_Nr_1_n, sigma_Nr_2_n, &
                              sigma_Nc_1_n, sigma_Nc_2_n, &
                              corr_sNc_1, corr_sNc_2, &
                              corr_srr_1_n, corr_srr_2_n, &
                              corr_sNr_1_n, corr_sNr_2_n, &
                              corr_sNc_1_n, corr_sNc_2_n, &
                              corr_rrNr_1_n, corr_rrNr_2_n, &
                              corr_sw, corr_wrr, corr_wNr, corr_wNc )

      !!! Calculate the values of the upscaled KK microphysics tendencies.
      call KK_upscaled_means_driver( rrainm(k), Nrm(k), Ncm(k), &
                                     mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                     mu_Nr_1, mu_Nr_2, mu_Nr_1_n, &
                                     mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                                     sigma_s_1, sigma_s_2, &
                                     sigma_rr_1_n, sigma_rr_2_n, &
                                     sigma_Nr_1_n, sigma_Nr_2_n, &
                                     sigma_Nc_1_n, sigma_Nc_2_n, &
                                     corr_srr_1_n, corr_srr_2_n, &
                                     corr_sNr_1_n, corr_sNr_2_n, &
                                     corr_sNc_1_n, corr_sNc_2_n, &
                                     corr_rrNr_1_n, corr_rrNr_2_n, &
                                     mixt_frac(k), precip_frac_1(k), &
                                     precip_frac_2(k), Nc0_in_cloud(k), &
                                     l_const_Nc_in_cloud, &
                                     KK_evap_coef, KK_auto_coef, &
                                     KK_accr_coef, KK_mvr_coef, &
                                     KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                     KK_accr_tndcy(k), KK_mean_vol_rad(k) )

      call KK_sed_vel_covars( rrainm(k), rr1(k), rr2(k), Nrm(k), &
                              Nr1(k), Nr2(k), KK_mean_vol_rad(k), &
                              mu_rr_1_n, mu_rr_2_n, mu_Nr_1, mu_Nr_2, &
                              mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1_n, &
                              sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                              corr_rrNr_1_n, corr_rrNr_2_n, KK_mvr_coef, &
                              mixt_frac(k), precip_frac_1(k), &
                              precip_frac_2(k), k, l_stats_samp, &
                              Vrrprrp_zt_impc(k), Vrrprrp_zt_expc(k), &
                              VNrpNrp_zt_impc(k), VNrpNrp_zt_expc(k) )


      if ( l_var_covar_src ) then

        call KK_upscaled_covar_driver( wm_zt(k), exner(k), rcm(k),  &
                                       rrainm(k), Nrm(k), Ncm(k), &
                                       mu_s_1, mu_s_2, mu_rr_1_n, mu_Nr_1_n, &
                                       mu_Nc_1_n, sigma_s_1, sigma_s_2, &
                                       sigma_rr_1_n, sigma_Nr_1_n, sigma_Nc_1_n, &
                                       corr_srr_1_n, corr_srr_2_n, &
                                       corr_sNr_1_n, corr_sNr_2_n, &
                                       corr_sNc_1_n, corr_sNc_2_n, &
                                       corr_rrNr_1_n,  mixt_frac(k), &
                                       precip_frac(k), Nc0_in_cloud(k), &
                                       l_const_Nc_in_cloud, &
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_evap_tndcy(k), &
                                       KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                       pdf_params(k), k, &
                                       l_stats_samp, &
                                       wprtp_mc_tndcy_zt(k), &
                                       wpthlp_mc_tndcy_zt(k), &
                                       rtp2_mc_tndcy_zt(k), &
                                       thlp2_mc_tndcy_zt(k), &
                                       rtpthlp_mc_tndcy_zt(k) )

      endif

      !!! KK rain drop concentration microphysics tendencies.

      !!! Calculate the KK N_r evaporation tendency.
      if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

         KK_Nrm_evap_tndcy(k)  &
         = KK_Nrm_evap_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                      mu_Nr_1_n, mu_Nr_2_n, sigma_s_1, &
                                      sigma_s_2, sigma_rr_1_n, sigma_rr_2_n, &
                                      sigma_Nr_1_n, sigma_Nr_2_n, &
                                      corr_srr_1_n, corr_srr_2_n, &
                                      corr_sNr_1_n, corr_sNr_2_n, &
                                      corr_rrNr_1_n, corr_rrNr_2_n, &
                                      KK_evap_coef, mixt_frac(k), &
                                      precip_frac_1(k), precip_frac_2(k), dt )

      else  ! r_r or N_r = 0.

          KK_Nrm_evap_tndcy(k) = zero

      endif

      !!! Calculate the KK N_r autoconversion tendency.
      KK_Nrm_auto_tndcy(k) = KK_Nrm_auto_mean( KK_auto_tndcy(k) )


      !!! Calculate any necessary adjustments to KK microphysics tendencies.
      call KK_microphys_adjust( dt, exner(k), rcm(k), rrainm(k), Nrm(k), &
                                KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                KK_accr_tndcy(k), KK_Nrm_evap_tndcy(k), &
                                KK_Nrm_auto_tndcy(k), l_src_adj_enabled, &
                                l_evap_adj_enabled, l_stats_samp, k, &
                                rrainm_mc_tndcy(k), Nrm_mc_tndcy(k), &
                                rvm_mc(k), rcm_mc(k), thlm_mc(k) )

      !!! Statistical output for upscaled KK.
      call KK_upscaled_stats( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                              mu_Nc_1, mu_Nc_2, mu_rr_1_n, mu_rr_2_n, &
                              mu_Nr_1_n, mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                              sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                              sigma_Nr_2, sigma_Nc_1, sigma_Nc_2, &
                              sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                              sigma_Nr_2_n, sigma_Nc_1_n, sigma_Nc_2_n, &
                              corr_srr_1, corr_srr_2, corr_sNr_1, &
                              corr_sNr_2, corr_sNc_1, corr_sNc_2, &
                              corr_rrNr_1, corr_rrNr_2, corr_srr_1_n, &
                              corr_srr_2_n, corr_sNr_1_n, corr_sNr_2_n, &
                              corr_sNc_1_n, corr_sNc_2_n, corr_rrNr_1_n, &
                              corr_rrNr_2_n, mixt_frac(k), precip_frac_1(k), &
                              precip_frac_2(k), KK_mvr_coef, &
                              KK_mean_vol_rad(k), rrainm(k), Nrm(k), k, &
                              l_stats_samp )

      !!! Statistical output for mean microphysics tendenices.
      call KK_stats_output( KK_evap_tndcy(k), KK_auto_tndcy(k), &
                            KK_accr_tndcy(k), KK_mean_vol_rad(k), &
                            KK_Nrm_evap_tndcy(k), KK_Nrm_auto_tndcy(k), &
                            l_stats_samp, k )

    enddo  ! Microphysics tendency loop: k = 2, nz, 1


    if ( l_var_covar_src ) then

       ! Output microphysics tendency terms for
       ! model variances and covariances on momentum levels.
       wprtp_mc_tndcy   = zt2zm( wprtp_mc_tndcy_zt )
       wpthlp_mc_tndcy  = zt2zm( wpthlp_mc_tndcy_zt )
       rtp2_mc_tndcy    = zt2zm( rtp2_mc_tndcy_zt )
       thlp2_mc_tndcy   = zt2zm( thlp2_mc_tndcy_zt )
       rtpthlp_mc_tndcy = zt2zm( rtpthlp_mc_tndcy_zt )

       ! Set values of microphysics tendency terms to 0 at model lower boundary.
       wprtp_mc_tndcy(1)   = zero
       wpthlp_mc_tndcy(1)  = zero
       rtp2_mc_tndcy(1)    = zero
       thlp2_mc_tndcy(1)   = zero
       rtpthlp_mc_tndcy(1) = zero
       ! Set values of microphysics tendency terms to 0 at model upper boundary.
       wprtp_mc_tndcy(nz)   = zero
       wpthlp_mc_tndcy(nz)  = zero
       rtp2_mc_tndcy(nz)    = zero
       thlp2_mc_tndcy(nz)   = zero
       rtpthlp_mc_tndcy(nz) = zero

    else

       ! Microphysics tendency terms for model variances and covariances
       ! are set to 0.
       wprtp_mc_tndcy   = zero
       wpthlp_mc_tndcy  = zero
       rtp2_mc_tndcy    = zero
       thlp2_mc_tndcy   = zero
       rtpthlp_mc_tndcy = zero

    endif

    !!! Boundary conditions for microphysics tendencies.

    ! Explicit contributions to rrainm and Nrm from microphysics are not set at
    ! thermodynamic level k = 1 because it is below the model lower boundary.
    rrainm_mc_tndcy(1) = zero
    Nrm_mc_tndcy(1)    = zero

    rrainm_mc_tndcy(nz) = zero
    Nrm_mc_tndcy(nz)    = zero

    ! Boundary conditions
    KK_mean_vol_rad(1)  = zero
    KK_mean_vol_rad(nz) = zero

    rvm_mc(1)  = zero
    rvm_mc(nz) = zero

    rcm_mc(1)  = zero
    rcm_mc(nz) = zero

    thlm_mc(1)  = zero
    thlm_mc(nz) = zero

    ! Find the vertical level index of cloud top.
    cloud_top_level = get_cloud_top_level( nz, rcm )

    !!! Microphysics sedimentation velocities.
    call KK_sedimentation( nz, cloud_top_level, KK_mean_vol_rad, Vrr, VNr )

    !!! Turbulent sedimentation above cloud top should have a value of 0.
    if ( cloud_top_level > 1 ) then
       Vrrprrp_zt_impc(cloud_top_level+1:nz-1) = zero
       Vrrprrp_zt_expc(cloud_top_level+1:nz-1) = zero
       VNrpNrp_zt_impc(cloud_top_level+1:nz-1) = zero
       VNrpNrp_zt_expc(cloud_top_level+1:nz-1) = zero
    endif

    !!! Boundary conditions (lower) for the covariances of hydrometeor
    !!! sedimentation velocities and their associated hydrometeors
    !!! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp_zt_impc(1) = Vrrprrp_zt_impc(2)
    Vrrprrp_zt_expc(1) = Vrrprrp_zt_expc(2)
    VNrpNrp_zt_impc(1) = VNrpNrp_zt_impc(2)
    VNrpNrp_zt_expc(1) = VNrpNrp_zt_expc(2)

    !!! Boundary conditions (upper) for the covariances of hydrometeor
    !!! sedimentation velocities and their associated hydrometeors
    !!! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp_zt_impc(nz) = zero
    Vrrprrp_zt_expc(nz) = zero
    VNrpNrp_zt_impc(nz) = zero
    VNrpNrp_zt_expc(nz) = zero

    ! Statistics
    if ( l_stats_samp ) then

       if ( iVrrprrp_expcalc > 0 ) then

          ! The covariance < V_rr'r_r' > calculated completely explicitly.
          ! When semi-implicit turbulent advection is used, this result can be
          ! compared to the < V_rr'r_r' > results used in the code, which are
          ! calculated semi-implicitly.
          call stat_update_var( iVrrprrp_expcalc, &
                                zt2zm( Vrrprrp_zt_impc * rrainm &
                                       + Vrrprrp_zt_expc ), zm )

       endif

       if ( iVNrpNrp_expcalc > 0 ) then

          ! The covariance < V_Nr'N_r' > calculated completely explicitly.
          ! When semi-implicit turbulent advection is used, this result can be
          ! compared to the < V_Nr'N_r' > results used in the code, which are
          ! calculated semi-implicitly.
          call stat_update_var( iVNrpNrp_expcalc, &
                                zt2zm( VNrpNrp_zt_impc * Nrm &
                                       + VNrpNrp_zt_expc ), zm )

       endif

    endif ! l_stats_samp


    return

  end subroutine KK_upscaled_micro_driver

  !=============================================================================
  subroutine KK_micro_init( nz, hydromet, hydromet_mc, hydromet_vel, &
                            rrainm, Nrm, Vrr, VNr, &
                            rrainm_mc_tndcy, Nrm_mc_tndcy, &
                            KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                            KK_mean_vol_rad, KK_Nrm_evap_tndcy, &
                            KK_Nrm_auto_tndcy, &
                            l_src_adj_enabled, l_evap_adj_enabled )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        iirrainm, & ! Constant(s)
        iiNrm

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet    ! Hydrometeor species                      [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(:), pointer, intent(out) ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
      Vrr,             & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,             & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      KK_evap_tndcy,     & ! Mean KK (dr_r/dt) due to evaporation    [(kg/kg)/s]
      KK_auto_tndcy,     & ! Mean KK (dr_r/dt) due to autoconversion [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK (dr_r/dt) due to accretion      [(kg/kg)/s]
      KK_mean_vol_rad,   & ! Mean KK rain drop mean volume radius            [m]
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation   [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.     [(num/kg)/s]

    logical, intent(out) :: &
      l_src_adj_enabled,  & ! Flag to enable rrainm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrainm/Nrm evaporation adjustment


    !!! Initialize microphysics tendencies and rain drop mean volume radius.
    KK_evap_tndcy = zero
    KK_auto_tndcy = zero
    KK_accr_tndcy = zero

    KK_mean_vol_rad = zero

    KK_Nrm_evap_tndcy = zero
    KK_Nrm_auto_tndcy = zero

    !!! Assign pointers for hydrometeor variables.

    ! Mean fields.
    rrainm => hydromet(:,iirrainm)
    Nrm    => hydromet(:,iiNrm)

    ! Sedimentation Velocities.
    Vrr => hydromet_vel(:,iirrainm)
    VNr => hydromet_vel(:,iiNrm)

    ! Mean field tendencies.
    rrainm_mc_tndcy => hydromet_mc(:,iirrainm)
    Nrm_mc_tndcy    => hydromet_mc(:,iiNrm)

    !!! Set KK microphysics tendency adjustment flags
    l_src_adj_enabled  = .true.
    l_evap_adj_enabled = .true.


    return

  end subroutine KK_micro_init

  !=============================================================================
  subroutine component_means_rain( nz, rrainm, Nrm, rho, rc1, rc2, &
                                   mixt_frac, l_stats_samp, &
                                   rr1, rr2, Nr1, Nr2 )

    ! Description:
    ! The values of grid-level mean rain water mixing ratio, <r_r>, and
    ! grid-level mean rain drop concentration, <N_r>, are solved as part of the
    ! CLUBB model's predictive equation set.  However, CLUBB has a two component
    ! PDF.  The grid-level means of r_r and N_r must be subdivided into
    ! component means for each PDF component.  The equations relating the
    ! overall means to the component means are:
    !
    ! <r_r> = a * rr1 + (1-a) * rr2, and
    ! <N_r> = a * Nr1 + (1-a) * Nr2;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), rr1
    ! is the mean rain water mixing ratio in PDF component 1, rr2 is the mean
    ! rain water mixing ratio in PDF component 2, Nr1 is the mean rain drop
    ! concentration in PDF component 1, and Nr2 is the mean rain drop
    ! concentration in PDF component 2.  These equations can be rewritten as:
    !
    ! <r_r> = rr1 * ( a + (1-a) * rr2/rr1 ), and
    ! <N_r> = Nr1 * ( a + (1-a) * Nr2/Nr1 ).
    !
    ! One way to solve for a component mean is to relate the ratios rr2/rr1 and
    ! Nr2/Nr1 to other factors.  For now, these ratios based on other factors
    ! will be called rr2_rr1_ratio (for rr2/rr1) and Nr2_Nr1_ratio
    ! (for Nr2/Nr1).  These ratios are entered into the above equations,
    ! allowing the equations to be solved for rr1 and Nr1:
    !
    ! rr1 = <r_r> / ( a + (1-a) * rr2_rr1_ratio ), and
    ! Nr1 = <N_r> / ( a + (1-a) * Nr2_Nr1_ratio ).
    !
    ! Once that rr1 and Nr1 have been solved, rr2 and Nr2 can be solved by:
    !
    ! rr2 = ( <r_r> - a * rr1 ) / (1-a); and
    ! Nr2 = ( <N_r> - a * Nr1 ) / (1-a).
    !
    ! At a grid level that is at least mostly cloudy, the simplest way to handle
    ! the ratios rr2/rr1 and Nr2/Nr1 is to set them equal to the ratio rc2/rc1,
    ! where rc1 is the mean cloud water mixing ratio in PDF component 1 and rc2
    ! is the mean cloud water mixing ratio in PDF component 2.  However, rain
    ! sediments, falling from higher altitudes downwards.  The values of cloud
    ! water mixing ratio at a given grid level are not necessarily indicative
    ! of the amount of cloud water at higher levels, which has already produced
    ! rain which has fallen downwards to the given grid level.  Additionally,
    ! using grid-level cloud water mixing ratio especially does not work for
    ! rain below cloud base (near the ground).
    !
    ! However, an alternative to component cloud water mixing ratio is component
    ! liquid water path.  Liquid water path accounts for the cloud water mixing
    ! ratio at the given grid level and at all grid levels higher in altitude.
    !
    ! In a stratocumulus case, the cloud water is spread out over all or almost
    ! all of the horizontal domain over a group of vertical levels.  At a given
    ! vertical level, the component mean cloud water mixing ratios should be
    ! almost equal, although usually slightly larger in the component with the
    ! larger component mean extended liquid water mixing ratio, s.  Likewise,
    ! the component liquid water paths should be nearly equal, with one
    ! component having a slightly larger liquid water path than the other
    ! component.
    ! 
    ! In a case of cumulus rising into stratocumulus, the upper portion of the
    ! cloudy domain will be very similar to the stratocumulus case described
    ! above, with similar cloud water mixing ratio and liquid water path
    ! results.  However, below the base of the stratocumulus clouds, where the
    ! cumulus clouds are found, the horizontal domain at each vertical level is
    ! only partially cloudy.  At these levels, rain produced in the
    ! stratocumulus clouds above is evaporating in the clear-air portions, while
    ! rain is not evaporating in the cloudy portions.  Additionally, more rain
    ! is being produced in the cloudy portions.  The rain in the cloudy portions
    ! becomes significantly larger than the rain in the clear portions.  The
    ! partiallly cloudy levels usually have a PDF where one component is
    ! significantly more saturated than the other component.  By the time the
    ! cloud base of the cumulus clouds is reached, the liquid water path for one
    ! PDF component should be significantly greater than the liquid water path
    ! for the other PDF component.
    ! 
    ! In a cumulus case, the horizontal domain at each level is usually partly
    ! cloudy.  Throughout the entire vertical domain, at every vertical level,
    ! one component usually is much more saturated than the other component.
    ! The liquid water path for one component is much greater than the liquid
    ! water path in the other component.  Likewise, rain that is formed in cloud
    ! and falls preferentially through cloud will have large values in a portion
    ! of the horizontal domain and very small or 0 values over the rest of the 
    ! horizontal domain.
    !
    ! In order to estimate the amount of rain in each PDF component, the ratios
    ! rr2/rr1 and Nr2/Nr1 are going to be set equal to the ratio LWP2/LWP1,
    ! where LWP1 is the liquid water path in PDF component 1 and LWP2 is the
    ! liquid water path in PDF component 2.  LWP1 will be computed by taking the
    ! vertical integral of cloud water (see equation below) through the 1st PDF
    ! component from the given vertical level all the way to the top of the
    ! model.  LWP2 will be computed in the same manner.   It should be noted
    ! that this method makes the poor assumption that PDF component 1 always
    ! overlaps PDF component 1 between vertical levels, and likewise for PDF
    ! component 2.
    !
    ! Total liquid water path, LWP, is given by the following equation:
    !
    ! LWP(z) = INT(z:z_top) rho_a <r_c> dz';
    !
    ! where z is the altitude of the vertical level for which LWP is desired,
    ! z_top is the altitude at the top of the model domain, and z' is the
    ! dummy variable of integration.  Mean cloud water mixing ratio can be
    ! written as:
    !
    ! <r_c> = a * rc1 + (1-a) * rc2.
    !
    ! The equation for liquid water path is rewritten as:
    !
    ! LWP(z) = INT(z:z_top) rho_a ( a rc1 + (1-a) rc2 ) dz'; or
    !
    ! LWP(z) = INT(z:z_top) a rho_a rc1 dz'
    !          + INT(z:z_top) (1-a) rho_a rc2 dz'.
    !
    ! This can be rewritten as:
    !
    ! LWP(z) = LWP1(z) + LWP2(z);
    !
    ! where:
    !
    ! LWP1(z) = INT(z:z_top) a rho_a rc1 dz'; and
    ! LWP2(z) = INT(z:z_top) (1-a) rho_a rc2 dz'.
    !
    ! The trapezoidal rule will be used to numerically integrate for LWP1
    ! and LWP2.
    

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one,      & ! Constant(s)
        one_half, &
        zero,     &
        rr_tol,   &
        Nr_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type, only: &
        stat_update_var  ! Procedure(s)

    use stats_variables, only : &
        iLWP1, & ! Variable(s)
        iLWP2, &
        zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz           ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rrainm,    & ! Overall mean rain water mixing ratio               [kg/kg]
      Nrm,       & ! Overall mean rain drop concentration               [num/kg]
      rho,       & ! Air density                                        [kg/m^3]
      rc1,       & ! Mean cloud water mixing ratio (1st PDF component)  [kg/kg]
      rc2,       & ! Mean cloud water mixing ratio (2nd PDF component)  [kg/kg]
      mixt_frac    ! Mixture fraction                                   [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: &
      LWP1, & ! Liquid water path (1st PDF component) on thermo. levs.  [kg/m^2]
      LWP2    ! Liquid water path (2nd PDF component) on thermo. levs.  [kg/m^2]

    integer :: k  ! Array index

    real( kind = core_rknd ), parameter :: &
      LWP_tol = 5.0e-7_core_rknd  ! Tolerance value for component LWP


    !!! Compute component liquid water paths using trapezoidal rule for
    !!! numerical integration.

    ! At the uppermost thermodynamic level (k = nz), use the trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k = nz,
    ! integrand_b is the integrand at momentum level k = nz (model upper
    ! boundary), and delta_z = zm(nz) - zt(nz).  At the upper boundary, r_c is
    ! set to 0, and the form of the trapezoidal rule is simply:
    !
    ! 0.5 * integrand_a * delta_z.

    ! Liquid water path in PDF component 1.
    LWP1(nz) &
    = one_half * mixt_frac(nz) * rho(nz) * rc1(nz) * ( gr%zm(nz) - gr%zt(nz) )

    ! Liquid water path in PDF component 2.
    LWP2(nz) &
    = one_half * ( one - mixt_frac(nz) ) * rho(nz) * rc2(nz) &
      * ( gr%zm(nz) - gr%zt(nz) )

    ! At all other thermodynamic levels, compute liquid water path using the
    ! trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k, integrand_b
    ! is the integrand at thermodynamic level k+1, and
    ! delta_z = zt(k+1) - zt(k), or 1/invrs_dzm(k).  The total for the segment
    ! is added to the sum total of all higher vertical segments to compute the
    ! total vertical integral.
    do k = nz-1, 1, -1

       ! Liquid water path in PDF component 1.
       LWP1(k) &
       = LWP1(k+1) &
         + one_half * ( mixt_frac(k+1) * rho(k+1) * rc1(k+1) &
                        + mixt_frac(k) * rho(k) * rc1(k) ) / gr%invrs_dzm(k)

       ! Liquid water path in PDF component 2.
       LWP2(k) &
       = LWP2(k+1) &
         + one_half * ( ( one - mixt_frac(k+1) ) * rho(k+1) * rc2(k+1) &
                        + ( one - mixt_frac(k) ) * rho(k) * rc2(k) ) &
           / gr%invrs_dzm(k)

    enddo ! k = nz-1, 1, -1


    !!! Find rr1, rr2, Nr1, and Nr2 based on the ratio of LWP2/LWP1, such that:
    !!! rr2/rr1 = Nr2/Nr1 = LWP2/LWP1.
    do k = 1, nz, 1

       !!! Calculate the component means for rain water mixing ratio.
       if ( rrainm(k) > rr_tol ) then

          if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
             !
             ! There is rain at this level, yet no cloud at or above the
             ! current level.  This is usually due to a numerical artifact.
             ! For example, rain is diffused above cloud top.  Simply set
             ! each component mean equal to the overall mean.
             rr1(k) = rrainm(k)
             rr2(k) = rrainm(k)

          elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 1st PDF component.  All rain water
             ! mixing ratio is found in the 1st PDF component.
             rr1(k) = rrainm(k) / mixt_frac(k)
             rr2(k) = zero

          elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

             ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 2nd PDF component.  All rain water
             ! mixing ratio is found in the 2nd PDF component.
             rr1(k) = zero
             rr2(k) = rrainm(k) / ( one - mixt_frac(k) )

          else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

             ! Both LWP1 and LWP2 are (significantly) greater than 0.
             !
             ! There is rain at this level, and there is sufficient cloud water
             ! at or above this level in both PDF components to find rain in
             ! both PDF components.  Delegate rain water mixing ratio between
             ! the 1st and 2nd PDF components according to the above equations.
             rr1(k) &
             = rrainm(k) &
               / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

             rr2(k) &
             = ( rrainm(k) - mixt_frac(k) * rr1(k) ) / ( one - mixt_frac(k) )

          endif


       else ! rrainm(k) <= rr_tol

          ! Overall mean rain water mixing ratio is either 0 or below tolerance
          ! value (any postive value is considered to be a numerical artifact).
          ! Simply set each component mean equal to the overall mean.
           rr1(k) = rrainm(k)
           rr2(k) = rrainm(k)

       endif


       !!! Calculate the component means for rain drop concentration.
       if ( Nrm(k) > Nr_tol ) then

          if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
             !
             ! There is rain at this level, yet no cloud at or above the
             ! current level.  This is usually due to a numerical artifact.
             ! For example, rain is diffused above cloud top.  Simply set
             ! each component mean equal to the overall mean.
             Nr1(k) = Nrm(k)
             Nr2(k) = Nrm(k)

          elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 1st PDF component.  All rain drop
             ! concentration is found in the 1st PDF component.
             Nr1(k) = Nrm(k) / mixt_frac(k)
             Nr2(k) = zero

          elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

             ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 2nd PDF component.  All rain drop
             ! concentration is found in the 2nd PDF component.
             Nr1(k) = zero
             Nr2(k) = Nrm(k) / ( one - mixt_frac(k) )

          else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

             ! Both LWP1 and LWP2 are (significantly) greater than 0.
             !
             ! There is rain at this level, and there is sufficient cloud water
             ! at or above this level in both PDF components to find rain in
             ! both PDF components.  Delegate rain drop concentration between
             ! the 1st and 2nd PDF components according to the above equations.
             Nr1(k) &
             = Nrm(k) &
               / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

             Nr2(k) &
             = ( Nrm(k) - mixt_frac(k) * Nr1(k) ) / ( one - mixt_frac(k) )

          endif


       else ! Nrm(k) <= Nr_tol

          ! Overall mean rain drop concentration is either 0 or below tolerance
          ! value (any postive value is considered to be a numerical artifact).
          ! Simply set each component mean equal to the overall mean.
           Nr1(k) = Nrm(k)
           Nr2(k) = Nrm(k)

       endif


    enddo ! k = 1, nz, 1

    ! Statistics
    if ( l_stats_samp ) then

       if ( iLWP1 > 0 ) then
          ! Liquid water path in PDF component 1.
          call stat_update_var( iLWP1, LWP1, zt )
       endif

       if ( iLWP2 > 0 ) then
          ! Liquid water path in PDF component 2.
          call stat_update_var( iLWP2, LWP2, zt )
       endif
       
    endif


    return

  end subroutine component_means_rain

  !=============================================================================
  subroutine precip_fraction( nz, rrainm, rr1, rr2, Nrm, Nr1, Nr2, &
                              cloud_frac, cloud_frac1, mixt_frac, &
                              precip_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! Determines (overall) precipitation fraction over the horizontal domain, as
    ! well as the precipitation fraction within each PDF component, at every
    ! vertical grid level.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        rr_tol,         &
        Nr_tol,         &
        cloud_frac_min

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rrainm,      & ! Mean rain water mixing ratio (overall)           [kg/kg]
      rr1,         & ! Mean rain water mixing ratio (1st PDF component) [kg/kg]
      rr2,         & ! Mean rain water mixing ratio (2nd PDF component) [kg/kg]
      Nrm,         & ! Mean rain drop concentration (overall)           [num/kg]
      Nr1,         & ! Mean rain drop concentration (1st PDF component) [num/kg]
      Nr2,         & ! Mean rain drop concentration (2nd PDF component) [num/kg]
      cloud_frac,  & ! Cloud fraction (overall)                         [-] 
      cloud_frac1, & ! Cloud fraction (1st PDF component)               [-]
      mixt_frac      ! Mixture fraction                                 [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      weighted_pfrac1    ! Product of mixt_frac and precip_frac_1       [-]

    real( kind = core_rknd ), parameter :: &
      precip_frac_tol = cloud_frac_min  ! Minimum precip. frac.         [-]
    
    integer :: &
      k   ! Loop index


    !!! Find overall precipitation fraction.
    do k = nz, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       if ( k < nz ) then
          precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )
       else  ! k = nz
          precip_frac(k) = cloud_frac(k)
       endif

       if ( ( rrainm(k) > rr_tol .or. Nrm(k) > Nr_tol ) &
            .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find rain at this grid level, but no cloud at
          ! or above this grid level, set precipitation fraction to a minimum
          ! threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( ( rrainm(k) < rr_tol .and. Nrm(k) < Nr_tol ) &
                .and. precip_frac(k) < precip_frac_tol ) then

          ! Mean (overall) rain water mixing ratio and mean (overall) rain drop
          ! concentration are both less than their respective tolerance amounts.
          ! They are both considered to have values of 0.  There is not any rain
          ! at this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo ! Overall precipitation fraction loop: k = nz, 1, -1.


    !!! Find precipitation fraction within each PDF component.
    !
    ! The overall precipitation fraction, f_p, is given by the equation:
    !
    ! f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    !
    ! where a is the mixture fraction (weight of PDF component 1), f_p(1) is
    ! the precipitation fraction within PDF component 1, and f_p(2) is the
    ! precipitation fraction within PDF component 2.  Overall precipitation
    ! fraction is found according the method above, and mixture fraction is
    ! already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    ! values for f_p(1) and f_p(2) must satisfy the above equation.
    if ( .false. ) then

       !!! Find precipitation fraction within PDF component 1.
       ! The method used to find overall precipitation fraction will also be to
       ! find precipitation fraction within PDF component 1.  In order to do so,
       ! it is assumed (poorly) that PDF component 1 overlaps PDF component 1 at
       ! every vertical level in the vertical profile.
       do k = nz, 1, -1

          ! The weighted precipitation fraction (PDF component 1) is the
          ! greatest value of the product of mixture fraction and cloud fraction
          ! (PDF component 1) at or above a vertical level.
          if ( k < nz ) then
             weighted_pfrac1(k) = max( weighted_pfrac1(k+1), &
                                       mixt_frac(k) * cloud_frac1(k) )
          else  ! k = nz
             weighted_pfrac1(k) = mixt_frac(k) * cloud_frac1(k)
          endif

          precip_frac_1(k) = weighted_pfrac1(k) / mixt_frac(k)

          ! Special cases for precip_frac_1.
          if ( precip_frac_1(k) > one ) then

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  For example, the mixture fraction at level k+1
             ! is 0.10 and the cloud_frac1 at level k+1 is 1, resulting in a
             ! weighted_pfrac1 of 0.10.  This product is greater than the
             ! product of mixt_frac and cloud_frac1 at level k.  The mixture
             ! fraction at level k is 0.05, resulting in a precip_frac_1 of 2.
             ! The value of precip_frac_1 is limited at 1.  The leftover
             ! precipitation fraction (a result of the decreasing weight of PDF
             ! component 1 between the levels) is applied to PDF component 2.
             precip_frac_1(k) = one

          elseif ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 1st PDF component at this
             ! grid level, but no cloud in the 1st PDF component at or above
             ! this grid level, set precipitation fraction (in the 1st PDF
             ! component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 1st PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif

       enddo ! Precipitation fraction (1st PDF component) loop: k = nz, 1, -1.


       !!! Find precipitation fraction within PDF component 2.
       ! The equation for precipitation fraction within PDF component 2 is:
       !
       ! f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
       !
       ! given the overall precipitation fraction, f_p (calculated above), the
       ! precipitation fraction within PDF component 1, f_p(1) (calculated
       ! above), and mixture fraction, a.  Any leftover precipitation fraction
       ! from precip_frac_1 will be included in this calculation of
       ! precip_frac_2.
       do k = 1, nz, 1

          precip_frac_2(k) &
          = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
            / ( one - mixt_frac(k) )

          ! Special cases for precip_frac_2.
          if ( precip_frac_2(k) > one ) then

             ! Again, it is possible for precip_frac_2 to be greater than 1.
             ! For example, the mixture fraction at level k+1 is 0.10 and the
             ! cloud_frac1 at level k+1 is 1, resulting in a weighted_pfrac1 of
             ! 0.10.  This product is greater than the product of mixt_frac and
             ! cloud_frac1 at level k.  Additionally, precip_frac (overall) is 1
             ! for level k.  The mixture fraction at level k is 0.5, resulting
             ! in a precip_frac_1 of 0.2.  Using the above equation,
             ! precip_frac_2 is calculated to be 1.8.  The value of
             ! precip_frac_2 is limited at 1.  The leftover precipitation
             ! fraction (as a result of the increasing weight of component 1
             ! between the levels) is applied to PDF component 1.
             precip_frac_2(k) = one

             ! Recalculate the precipitation fraction in PDF component 1.
             precip_frac_1(k) &
             = ( precip_frac(k) - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
               / mixt_frac(k)

             ! Double check for errors in PDF component 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             elseif ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = precip_frac_tol
             elseif ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = zero
             endif

          elseif ( ( rr2(k) > rr_tol .or. Nr2(k) > Nr_tol ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 2nd PDF component at this
             ! grid level, but no cloud in the 2nd PDF component at or above
             ! this grid level, set precipitation fraction (in the 2nd PDF
             ! component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 2nd PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif

       enddo ! Precipitation fraction (2nd PDF component) loop: k = 1, nz, 1.


    else  ! .true.

       ! Precipitation fraction in each PDF component is based on mean rain
       ! water mixing ratio in each PDF component.  The ratio it is based on is:
       !
       ! rr1/f_p(1) = rr2/f_p(2);
       !
       ! which can be rewritten as:
       !
       ! f_p(2)/f_p(1) = rr2/rr1.
       !
       ! Since overall precipitation fraction is given by the equation:
       !
       ! f_p = a f_p(1) + (1-a) f_p(2);
       !
       ! it can be rewritten as:
       !
       ! f_p = f_p(1) ( a + (1-a) f_p(2)/f_p(1) ).
       !
       ! Substituting the ratio rr2/rr1 for the ratio f_p(2)/f_p(1), the above
       ! equation can be solved for f_p(1):
       !
       ! f_p(1) = f_p / ( a + (1-a) rr2/rr1 ).
       !
       ! Then, f_p(2) can be solved for according to the equation:
       !
       ! f_p(2) = ( f_p - a f_p(1) ) / (1-a).
       do k = 1, nz, 1

          if ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
               .and. ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) ) then

             ! There is no rain in each PDF component.  Precipitation fraction
             ! within each component is set to 0.
             precip_frac_1(k) = zero
             precip_frac_2(k) = zero

          elseif ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
                   .and. ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) ) then

             ! All the rain is within the 1st PDF component.
             precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
             precip_frac_2(k) = zero

          elseif ( ( rr2(k) > rr_tol .or. Nr2(k) > Nr_tol ) &
                   .and. ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) ) then

             ! All the rain is within the 2nd PDF component.
             precip_frac_1(k) = zero
             precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )

          else ! rr1(k) > rr_tol or Nr1(k) > Nr_tol
               ! AND rr2(k) > rr_tol or Nr2(k) > Nr_tol

             ! Rain within both PDF components.

             !!! Find precipitation fraction within PDF component 1.
             precip_frac_1(k) &
             = precip_frac(k) &
               / ( mixt_frac(k) + ( one - mixt_frac(k) ) * rr2(k)/rr1(k) )

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  The value of precip_frac_1 is limited at 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             endif

             !!! Find precipitation fraction within PDF component 2.
             precip_frac_2(k) &
             = ( precip_frac(k) - mixt_frac(k) *  precip_frac_1(k) ) &
               / ( one - mixt_frac(k) )

             ! Using the above method, it is possible for precip_frac_2 to be
             ! greater than 1.  The value of precip_frac_2 is limited at 1.
             if ( precip_frac_2(k) > one ) then

                precip_frac_2(k) = one

                ! Recalculate the precipitation fraction in PDF component 1.
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

             endif

          endif


          ! Special cases for PDF component 1.
          if ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
               .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 1st PDF component at this
             ! grid level, but no cloud in the 1st PDF component at or above
             ! this grid level, set precipitation fraction (in the 1st PDF
             ! component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 1st PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif


          ! Special cases for PDF component 2.
          if ( ( rr2(k) > rr_tol .or. Nr2(k) > Nr_tol ) &
               .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 2nd PDF component at this
             ! grid level, but no cloud in the 2nd PDF component at or above
             ! this grid level, set precipitation fraction (in the 2nd PDF
             ! component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 2nd PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif


       enddo ! Component precipitation fraction loop: k = 1, nz, 1.


    endif ! Select component precipitation fraction method.


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine KK_tendency_coefs( thlm, exner, p_in_Pa, rho, &
                                KK_evap_coef, KK_auto_coef, &
                                KK_accr_coef, KK_mvr_coef )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use saturation, only: &
        sat_mixrat_liq  ! Procedure(s)

    use constants_clubb, only: &
        Rd,           & ! Constant(s)
        Rv,           &
        Lv,           &
        Cp,           &
        pi,           &
        three,        &
        four_thirds,  &
        one,          &
        two_thirds,   &
        one_third,    &
        rho_lw,       &
        cm3_per_m3

    use parameters_microphys, only: &
        KK_auto_Nc_exp,      & ! Constant(s)
        C_evap

    use KK_utilities, only: &
        G_T_p  ! Procedure(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      thlm,    & ! Mean liquid water potential temperature         [K]
      p_in_Pa, & ! Pressure                                        [Pa]
      exner,   & ! Exner function                                  [-]
      rho        ! Density                                         [kg/m^3]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      KK_evap_coef, & ! KK evaporation coefficient                 [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient              [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                   [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient          [m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      T_liq_in_K, & ! Mean liquid water temperature, T_l           [K]
      r_sl,       & ! Liquid water sat. mixing ratio, r_s(T_l,p)   [kg/kg]
      Beta_Tl       ! Parameter Beta, Beta(T_l)                    [1/(kg/kg)]


    ! Liquid water temperature.
    T_liq_in_K = thlm * exner

    ! Saturation mixing ratio (based on liquid water temperature and
    ! pressure), r_sl = r_s(T_l,p).
    r_sl = sat_mixrat_liq( p_in_Pa, T_liq_in_K )

    ! Beta(T_l).
    Beta_Tl = (Rd/Rv) * ( Lv / ( Rd * T_liq_in_K ) )  &
                      * ( Lv / ( Cp * T_liq_in_K ) )

    ! Coefficient for KK evaporation.
    KK_evap_coef = three * C_evap * G_T_p( T_liq_in_K, p_in_Pa )   &
                         * ( four_thirds * pi * rho_lw )**two_thirds  &
                         * ( ( one + Beta_Tl * r_sl ) / r_sl )

    ! Coefficient for KK autoconversion.
    KK_auto_coef = 1350.0_core_rknd * ( rho / cm3_per_m3 )**KK_auto_Nc_exp

    ! Coefficient for KK accretion.
    KK_accr_coef = 67.0_core_rknd

    ! Coefficient for KK rain drop mean volume radius.
    KK_mvr_coef = ( four_thirds * pi * rho_lw )**(-one_third)


    return

  end subroutine KK_tendency_coefs

  !=============================================================================
  subroutine KK_in_precip_values( rr1, rr2, Nr1, Nr2, rc1, &
                                  rc2, cloud_frac1, cloud_frac2, &
                                  precip_frac_1, precip_frac_2, &
                                  mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                  sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                  sigma_Nr_2, corr_srr_1, corr_srr_2, &
                                  corr_sNr_1, corr_sNr_2, corr_rrNr_1, &
                                  corr_rrNr_2 )
       
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol, & ! Constant(s)
        rr_tol, & 
        Nr_tol, &
        one,    & 
        zero

    use parameters_microphys, only: &
        rrp2_on_rrm2_cloud, & ! Variable(s)
        rrp2_on_rrm2_below, &
        Nrp2_on_Nrm2_cloud, &
        Nrp2_on_Nrm2_below

    use KK_fixed_correlations, only: &
        corr_srr_NL_cloud,  & ! Variable(s) 
        corr_srr_NL_below,  &
        corr_sNr_NL_cloud,  &
        corr_sNr_NL_below,  &
        corr_rrNr_LL_cloud, &
        corr_rrNr_LL_below

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rr1,           & ! Mean rain water mixing ratio (1st PDF comp.)  [kg/kg]
      rr2,           & ! Mean rain water mixing ratio (2nd PDF comp.)  [kg/kg]
      Nr1,           & ! Mean rain drop concentration (1st PDF comp.)  [num/kg]
      Nr2,           & ! Mean rain drop concentration (2nd PDF comp.)  [num/kg]
      rc1,           & ! Mean of r_c (1st PDF component)               [kg/kg]
      rc2,           & ! Mean of r_c (2nd PDF component)               [kg/kg]
      cloud_frac1,   & ! Cloud fraction (1st PDF component)            [-]
      cloud_frac2,   & ! Cloud fraction (2nd PDF component)            [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)    [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)    [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_rr_1,     & ! Mean of rr (1st PDF component) in-precip (ip)     [kg/kg]
      mu_rr_2,     & ! Mean of rr (2nd PDF component) ip                 [kg/kg]
      mu_Nr_1,     & ! Mean of Nr (1st PDF component) ip                [num/kg]
      mu_Nr_2,     & ! Mean of Nr (2nd PDF component) ip                [num/kg]
      sigma_rr_1,  & ! Standard deviation of rr (1st PDF component) ip   [kg/kg]
      sigma_rr_2,  & ! Standard deviation of rr (2nd PDF component) ip   [kg/kg]
      sigma_Nr_1,  & ! Standard deviation of Nr (1st PDF component) ip  [num/kg]
      sigma_Nr_2,  & ! Standard deviation of Nr (2nd PDF component) ip  [num/kg]
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component) ip   [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component) ip   [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component) ip   [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component) ip   [-]
      corr_rrNr_1, & ! Correlation between rr and Nr (1st PDF component) ip  [-]
      corr_rrNr_2    ! Correlation between rr and Nr (2nd PDF component) ip  [-]

    ! Local Variable

    ! Prescribed parameters are set to in-cloud or outside-cloud (below-cloud)
    ! values based on whether or not cloud water mixing ratio has a value of at
    ! least rc_tol.  However, this does not take into account the amount of
    ! cloudiness in a component, just whether or not there is any cloud in the
    ! component.  The option l_interp_prescribed_params allows for an
    ! interpolated value between the in-cloud and below-cloud parameter value
    ! based on the component cloud fraction.
    logical, parameter :: &
      l_interp_prescribed_params = .false.


    ! Mean of in-precip rain water mixing ratio in PDF component 1.
    if ( rr1 > rr_tol ) then
       mu_rr_1 = rr1 / precip_frac_1
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.
       mu_rr_1 = zero
    endif

    ! Mean of in-precip rain water mixing ratio in PDF component 2.
    if ( rr2 > rr_tol ) then
       mu_rr_2 = rr2 / precip_frac_2
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.
       mu_rr_2 = zero
    endif

    ! Mean of in-precip rain drop concentration in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       mu_Nr_1 = Nr1 / precip_frac_1
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.
       mu_Nr_1 = zero
    endif

    ! Mean of in-precip rain drop concentration in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       mu_Nr_2 = Nr2 / precip_frac_2
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.
       mu_Nr_2 = zero
    endif

    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    ! Brian Griffin; February 3, 2007.

    ! Standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_rr_1 = sqrt( cloud_frac1 * rrp2_on_rrm2_cloud &
                             + ( one - cloud_frac1 ) * rrp2_on_rrm2_below ) &
                       * mu_rr_1
       else
          if ( rc1 > rc_tol ) then
             sigma_rr_1 = sqrt( rrp2_on_rrm2_cloud ) * mu_rr_1
          else
             sigma_rr_1 = sqrt( rrp2_on_rrm2_below ) * mu_rr_1
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain water mixing ratio does not vary in
       ! this component at this grid level.
       sigma_rr_1 = zero
    endif

    ! Standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_rr_2 = sqrt( cloud_frac2 * rrp2_on_rrm2_cloud &
                             + ( one - cloud_frac2 ) * rrp2_on_rrm2_below ) &
                       * mu_rr_2
       else
          if ( rc2 > rc_tol ) then
             sigma_rr_2 = sqrt( rrp2_on_rrm2_cloud ) * mu_rr_2
          else
             sigma_rr_2 = sqrt( rrp2_on_rrm2_below ) * mu_rr_2
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain water mixing ratio does not vary in
       ! this component at this grid level.
       sigma_rr_2 = zero
    endif

    ! Standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_Nr_1 = sqrt( cloud_frac1 * Nrp2_on_Nrm2_cloud &
                             + ( one - cloud_frac1 ) * Nrp2_on_Nrm2_below ) &
                       * mu_Nr_1
       else
          if ( rc1 > rc_tol ) then
             sigma_Nr_1 = sqrt( Nrp2_on_Nrm2_cloud ) * mu_Nr_1
          else
             sigma_Nr_1 = sqrt( Nrp2_on_Nrm2_below ) * mu_Nr_1
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain drop concentration does not vary in
       ! this component at this grid level.
       sigma_Nr_1 = zero
    endif

    ! Standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_Nr_2 = sqrt( cloud_frac2 * Nrp2_on_Nrm2_cloud &
                             + ( one - cloud_frac2 ) * Nrp2_on_Nrm2_below ) &
                       * mu_Nr_2
       else
          if ( rc2 > rc_tol ) then
             sigma_Nr_2 = sqrt( Nrp2_on_Nrm2_cloud ) * mu_Nr_2
          else
             sigma_Nr_2 = sqrt( Nrp2_on_Nrm2_below ) * mu_Nr_2
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain drop concentration does not vary in
       ! this component at this grid level.
       sigma_Nr_2 = zero
    endif

    ! Correlation (in-precip) between s and r_r in PDF component 1.
    if ( rr1 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_srr_1 = cloud_frac1 * corr_srr_NL_cloud &
                       + ( one - cloud_frac1 ) * corr_srr_NL_below
       else
          if ( rc1 > rc_tol ) then
             corr_srr_1 = corr_srr_NL_cloud
          else
             corr_srr_1 = corr_srr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 1st PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_srr_1 = zero
    endif

    ! Correlation (in-precip) between s and r_r in PDF component 2.
    if ( rr2 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_srr_2 = cloud_frac2 * corr_srr_NL_cloud &
                       + ( one - cloud_frac2 ) * corr_srr_NL_below
       else
          if ( rc2 > rc_tol ) then
             corr_srr_2 = corr_srr_NL_cloud
          else
             corr_srr_2 = corr_srr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 2nd PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_srr_2 = zero
    endif

    ! Correlation (in-precip) between s and N_r in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_sNr_1 = cloud_frac1 * corr_sNr_NL_cloud &
                       + ( one - cloud_frac1 ) * corr_sNr_NL_below
       else
          if ( rc1 > rc_tol ) then
             corr_sNr_1 = corr_sNr_NL_cloud
          else
             corr_sNr_1 = corr_sNr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 1st PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_sNr_1 = zero
    endif

    ! Correlation (in-precip) between s and N_r in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_sNr_2 = cloud_frac2 * corr_sNr_NL_cloud &
                       + ( one - cloud_frac2 ) * corr_sNr_NL_below
       else
          if ( rc2 > rc_tol ) then
             corr_sNr_2 = corr_sNr_NL_cloud
          else
             corr_sNr_2 = corr_sNr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 2nd PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_sNr_2 = zero
    endif

    ! Correlation (in-precip) between r_r and N_r in PDF component 1.
    if ( rr1 > rr_tol .and. Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_rrNr_1 = cloud_frac1 * corr_rrNr_LL_cloud &
                        + ( one - cloud_frac1 ) * corr_rrNr_LL_below
       else
          if ( rc1 > rc_tol ) then
             corr_rrNr_1 = corr_rrNr_LL_cloud
          else
             corr_rrNr_1 = corr_rrNr_LL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 and (or) mean
       ! in-precip rain drop concentration in PDF component 1 are (is) less than
       ! their (its) respective tolerance amount(s), and are (is) considered to
       ! have a value of 0.  There is not any rain in the 1st PDF component at
       ! this grid level.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_1 = zero
    endif

    ! Correlation (in-precip) between r_r and N_r in PDF component 2.
    if ( rr2 > rr_tol .and. Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_rrNr_2 = cloud_frac2 * corr_rrNr_LL_cloud &
                        + ( one - cloud_frac2 ) * corr_rrNr_LL_below
       else
          if ( rc2 > rc_tol ) then
             corr_rrNr_2 = corr_rrNr_LL_cloud
          else
             corr_rrNr_2 = corr_rrNr_LL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 and (or) mean
       ! in-precip rain drop concentration in PDF component 2 are (is) less than
       ! their (its) respective tolerance amount(s), and are (is) considered to
       ! have a value of 0.  There is not any rain in the 2nd PDF component at
       ! this grid level.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_2 = zero
    endif


    return    

  end subroutine KK_in_precip_values

  !=============================================================================
  subroutine KK_upscaled_setup( rcm, rrainm, Nrm, Ncm, &
                                rr1, rr2, Nr1, Nr2, &
                                mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                sigma_rr_1, sigma_rr_2, &
                                sigma_Nr_1, sigma_Nr_2, &
                                wpsp, wprrp_ip, wpNrp_ip, &
                                wpNcp, stdev_w, mixt_frac, &
                                pdf_params, &
                                corr_srr_1, corr_srr_2, corr_sNr_1, &
                                corr_sNr_2, corr_rrNr_1, corr_rrNr_2, &
                                mu_s_1, mu_s_2, mu_Nc_1, mu_Nc_2, &
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                                sigma_s_1, sigma_s_2, &
                                sigma_Nc_1, sigma_Nc_2, &
                                sigma_rr_1_n, sigma_rr_2_n, &
                                sigma_Nr_1_n, sigma_Nr_2_n, &
                                sigma_Nc_1_n, sigma_Nc_2_n, &
                                corr_sNc_1, corr_sNc_2, &
                                corr_srr_1_n, corr_srr_2_n, &
                                corr_sNr_1_n, corr_sNr_2_n, &
                                corr_sNc_1_n, corr_sNc_2_n, &
                                corr_rrNr_1_n, corr_rrNr_2_n, &
                                corr_sw, corr_wrr, corr_wNr, corr_wNc )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol, & ! Constant(s)
        rr_tol, & 
        Nr_tol, & 
        Nc_tol, &
        zero

    use KK_utilities, only: &
        mean_L2N,   & ! Procedure(s)
        stdev_L2N,  &
        corr_NL2NN, &
        corr_LL2NN

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use parameters_microphys, only: &
        Ncp2_on_Ncm2_cloud, & ! Variable(s)
        Ncp2_on_Ncm2_below

    use KK_fixed_correlations, only: &
        corr_wrr_NL_cloud,  & ! Variable(s) 
        corr_wNr_NL_cloud,  & 
        corr_wNc_NL_cloud,  & 
        corr_sw_NN_cloud,   & 
        corr_sNc_NL_cloud,  &
        corr_sNr_NL_cloud,  &
        corr_srr_NL_cloud,  &
        corr_rrNr_LL_cloud, &
        corr_wrr_NL_below,  & 
        corr_wNr_NL_below,  & 
        corr_wNc_NL_below,  & 
        corr_sw_NN_below,   & 
        corr_sNc_NL_below,  &
        corr_sNr_NL_below,  &
        corr_srr_NL_below,  &
        corr_rrNr_LL_below

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_diagnose_correlations, & ! Variable(s)
        l_calc_w_corr

    use diagnose_correlations_module, only: &
        diagnose_KK_corr, & ! Procedure(s)
        calc_mean, &
        calc_w_corr

    use constants_clubb, only: &
        w_tol,        & ! [m/s]
        s_mellor_tol, & ! [kg/kg]
        Nc_tol,       & ! [#/kg]
        rr_tol,       & ! [kg/kg] 
        Nr_tol          ! [#/kg]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rcm,        & ! Mean cloud water mixing ratio                      [kg/kg]
      rrainm,     & ! Mean rain water mixing ratio                       [kg/kg]
      Nrm,        & ! Mean rain drop concentration                      [num/kg]
      Ncm,        & ! Mean cloud droplet concentration                  [num/kg]
      rr1,        & ! Mean rain water mixing ratio (1st PDF component)   [kg/kg]
      rr2,        & ! Mean rain water mixing ratio (2nd PDF component)   [kg/kg]
      Nr1,        & ! Mean rain drop concentration (1st PDF component)  [num/kg]
      Nr2,        & ! Mean rain drop concentration (2nd PDF component)  [num/kg]
      mu_rr_1,    & ! Mean of rr (1st PDF component) in-precip (ip)      [kg/kg]
      mu_rr_2,    & ! Mean of rr (2nd PDF component) ip                  [kg/kg]
      mu_Nr_1,    & ! Mean of Nr (1st PDF component) ip                 [num/kg]
      mu_Nr_2,    & ! Mean of Nr (2nd PDF component) ip                 [num/kg]
      sigma_rr_1, & ! Standard deviation of rr (1st PDF component) ip    [kg/kg]
      sigma_rr_2, & ! Standard deviation of rr (2nd PDF component) ip    [kg/kg]
      sigma_Nr_1, & ! Standard deviation of Nr (1st PDF component) ip   [num/kg]
      sigma_Nr_2, & ! Standard deviation of Nr (2nd PDF component) ip   [num/kg]
      wpsp,       & ! Covariance of w and s                         [(m/s)kg/kg]
      wprrp_ip,   & ! Covariance of w and r_r (overall) ip          [(m/s)kg/kg]
      wpNrp_ip,   & ! Covariance of w and N_r (overall) ip         [(m/s)num/kg]
      wpNcp,      & ! Covariance of w and N_c                      [(m/s)num/kg]
      stdev_w,    & ! Standard deviation of w                              [m/s]
      mixt_frac     ! Mixture fraction                                       [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout) :: &
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component) ip   [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component) ip   [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component) ip   [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component) ip   [-]
      corr_rrNr_1, & ! Correlation between rr and Nr (1st PDF component) ip  [-]
      corr_rrNr_2    ! Correlation between rr and Nr (2nd PDF component) ip  [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_Nc_1,       & ! Mean of Nc (1st PDF component)                 [num/kg]
      mu_Nc_2,       & ! Mean of Nc (2nd PDF component)                 [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Nc_1_n,     & ! Mean of ln Nc (1st PDF component)          [ln(num/kg)]
      mu_Nc_2_n,     & ! Mean of ln Nc (2nd PDF component)          [ln(num/kg)]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_Nc_1,    & ! Standard deviation of Nc (1st PDF component)   [num/kg]
      sigma_Nc_2,    & ! Standard deviation of Nc (2nd PDF component)   [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Nc_1_n,  & ! Standard dev. of ln Nc (1st PDF comp.)     [ln(num/kg)]
      sigma_Nc_2_n,  & ! Standard dev. of ln Nc (2nd PDF comp.)     [ln(num/kg)]
      corr_sNc_1,    & ! Correlation between s and Nc (1st PDF component)    [-]
      corr_sNc_2,    & ! Correlation between s and Nc (2nd PDF component)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF comp.)     [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF comp.)     [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      corr_sw,       & ! Correlation between s & w (both components)         [-]
      corr_wrr,      & ! Correlation between rr & w (both components)        [-]
      corr_wNr,      & ! Correlation between Nr & w (both components)        [-]
      corr_wNc         ! Correlation between Nc & w (both components)        [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      rrp2_on_rrm2,    & ! Ratio of < r_r'^2 > to < r_r >^2              [-]
      Nrp2_on_Nrm2,    & ! Ratio of < N_r'^2 > to < N_r >^2              [-]
      Ncp2_on_Ncm2,    & ! Ratio of < N_c'^2 > to < N_c >^2              [-]
      s_mellor_m,      & ! Mean of s_mellor                              [kg/kg]
      stdev_s_mellor     ! Standard deviation of s_mellor                [kg/kg]


    ! --- Begin Code ---

    ! Enter the PDF parameters.
    mu_s_1    = pdf_params%s1
    mu_s_2    = pdf_params%s2
    sigma_s_1 = pdf_params%stdev_s1
    sigma_s_2 = pdf_params%stdev_s2

    ! Mean of cloud droplet concentration in PDF component 1.
    if ( Ncm > Nc_tol ) then
       mu_Nc_1 = Ncm
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.
       mu_Nc_1 = zero
    endif

    ! Mean of cloud droplet concentration in PDF component 2.
    if ( Ncm > Nc_tol ) then
       mu_Nc_2 = Ncm
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.
       mu_Nc_2 = zero
    endif

    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    ! Brian Griffin; February 3, 2007.
    !
    ! Set the value of the parameters based on whether the altitude is above or
    ! below cloud base.  Determine whether there is cloud at any given vertical
    ! level.  In order for a vertical level to have cloud, the amount of cloud
    ! water (rcm) must be greater than or equal to the tolerance level (rc_tol).
    ! If there is cloud at a given vertical level, then the ###_cloud value is
    ! used.  Otherwise, the ###_below value is used.

    ! Standard deviation of cloud droplet concentration in PDF component 1.
    if ( Ncm > Nc_tol ) then
       if ( rcm > rc_tol ) then
          sigma_Nc_1 = sqrt( Ncp2_on_Ncm2_cloud ) * mu_Nc_1
       else
          sigma_Nc_1 = sqrt( Ncp2_on_Ncm2_below ) * mu_Nc_1
       endif
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The standard deviation is simply 0 since cloud droplet
       ! concentration does not vary at this grid level.
       sigma_Nc_1 = zero
    endif

    ! Standard deviation of cloud droplet concentration in PDF component 2.
    if ( Ncm > Nc_tol ) then
       if ( rcm > rc_tol ) then
          sigma_Nc_2 = sqrt( Ncp2_on_Ncm2_cloud ) * mu_Nc_2
       else
          sigma_Nc_2 = sqrt( Ncp2_on_Ncm2_below ) * mu_Nc_2
       endif
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The standard deviation is simply 0 since cloud droplet
       ! concentration does not vary at this grid level.
       sigma_Nc_2 = zero
    endif

    ! Correlation between s and N_c in PDF component 1.
    if ( Ncm > Nc_tol ) then
       if ( rcm > rc_tol ) then
          corr_sNc_1 = corr_sNc_NL_cloud
       else
          corr_sNc_1 = corr_sNc_NL_below
       endif
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The correlations involving cloud droplet concentration
       ! are 0 since cloud droplet concentration does not vary at this grid
       ! level.
       corr_sNc_1 = zero
    endif

    ! Correlation between s and N_c in PDF component 2.
    if ( Ncm > Nc_tol ) then
       if ( rcm > rc_tol ) then
          corr_sNc_2 = corr_sNc_NL_cloud
       else
          corr_sNc_2 = corr_sNc_NL_below
       endif
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The correlations involving cloud droplet concentration
       ! are 0 since cloud droplet concentration does not vary at this grid
       ! level.
       corr_sNc_2 = zero
    endif

    if ( l_calc_w_corr ) then

      s_mellor_m &
      = calc_mean( pdf_params%mixt_frac, pdf_params%s1, pdf_params%s2 )

      stdev_s_mellor &
        = sqrt( pdf_params%mixt_frac &
                * ( ( pdf_params%s1 - s_mellor_m )**2 &
                    + pdf_params%stdev_s1**2 ) &
              + ( 1 - pdf_params%mixt_frac ) &
                * ( ( pdf_params%s2 - s_mellor_m )**2 &
                    + pdf_params%stdev_s2**2 ) )

          corr_sw &
          = calc_w_corr( wpsp, stdev_w, stdev_s_mellor, w_tol, s_mellor_tol )
          corr_wrr = calc_w_corr( wprrp_ip, stdev_w, sigma_rr_1, w_tol, rr_tol )
          corr_wNr = calc_w_corr( wpNrp_ip, stdev_w, sigma_Nr_1, w_tol, Nr_tol )
          corr_wNc = calc_w_corr( wpNcp, stdev_w, sigma_Nc_1, w_tol, Nc_tol )
      
    else ! .not. l_calc_w_corr

       if ( rcm > rc_tol ) then
  
          corr_sw  = corr_sw_NN_cloud
          corr_wrr = corr_wrr_NL_cloud
          corr_wNr = corr_wNr_NL_cloud
          corr_wNc = corr_wNc_NL_cloud

       else
  
          corr_sw  = corr_sw_NN_below
          corr_wrr = corr_wrr_NL_below
          corr_wNr = corr_wNr_NL_below
          corr_wNc = corr_wNc_NL_below

       endif 

    endif ! l_calc_w_corr


    if ( l_diagnose_correlations ) then

       if ( rrainm > rr_tol ) then
          rrp2_on_rrm2 = (sigma_rr_1/mu_rr_1)**2
       else
          ! The ratio is undefined; set it equal to 0.
          rrp2_on_rrm2 = zero
       endif

       if ( Nrm > Nr_tol ) then
          Nrp2_on_Nrm2 = (sigma_Nr_1/mu_Nr_1)**2
       else
          ! The ratio is undefined; set it equal to 0.
          Nrp2_on_Nrm2 = zero
       endif

       if ( Ncm > Nc_tol ) then
          Ncp2_on_Ncm2 = (sigma_Nc_1/mu_Nc_1)**2
       else
          ! The ratio is undefined; set it equal to 0.
          Ncp2_on_Ncm2 = zero
       endif

       if ( rcm > rc_tol ) then
         call diagnose_KK_corr( Ncm, rrainm, Nrm, &
                                Ncp2_on_Ncm2, rrp2_on_rrm2, Nrp2_on_Nrm2, &
                                corr_sw, corr_wrr, corr_wNr, corr_wNc,  &
                                pdf_params, &
                                corr_rrNr_LL_cloud, corr_srr_NL_cloud, &
                                corr_sNr_NL_cloud, corr_sNc_NL_cloud, & 
                                corr_rrNr_1, corr_srr_1, &
                                corr_sNr_1, corr_sNc_1 )
       else
         call diagnose_KK_corr( Ncm, rrainm, Nrm, &
                                Ncp2_on_Ncm2, rrp2_on_rrm2, Nrp2_on_Nrm2, &
                                corr_sw, corr_wrr, corr_wNr, corr_wNc,  &
                                pdf_params, &
                                corr_rrNr_LL_below, corr_srr_NL_below, &
                                corr_sNr_NL_below, corr_sNc_NL_below, & 
                                corr_rrNr_1, corr_srr_1, &
                                corr_sNr_1, corr_sNc_1 )
       endif

       corr_srr_2  = corr_srr_1
       corr_sNr_2  = corr_sNr_1
       corr_sNc_2  = corr_sNc_1
       corr_rrNr_2 = corr_rrNr_1

    endif


    !!! Calculate the normalized mean of variables that have an assumed
    !!! lognormal distribution, given the mean and variance of those variables.

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 1.
    if ( rr1 > rr_tol ) then
       mu_rr_1_n = mean_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! value of mu_rr_1_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_rr_1_n = -huge( mu_rr_1_n )
    endif

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 2.
    if ( rr2 > rr_tol ) then
       mu_rr_2_n = mean_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! value of mu_rr_2_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_rr_2_n = -huge( mu_rr_2_n )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       mu_Nr_1_n = mean_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! value of mu_Nr_1_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_Nr_1_n = -huge( mu_Nr_1_n )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       mu_Nr_2_n = mean_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! value of mu_Nr_2_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_Nr_2_n = -huge( mu_Nr_2_n )
    endif

    ! Normalized mean of cloud droplet concentration in PDF component 1.
    if ( Ncm > Nc_tol ) then
       mu_Nc_1_n = mean_L2N( mu_Nc_1, sigma_Nc_1**2 )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There isn't any cloud at this
       ! grid level.  The value of mu_Nc_n should be -inf.  It will be set to
       ! -huge for purposes of assigning it a value.
       mu_Nc_1_n = -huge( mu_Nc_1_n )
    endif

    ! Normalized mean of cloud droplet concentration in PDF component 2.
    if ( Ncm > Nc_tol ) then
       mu_Nc_2_n = mean_L2N( mu_Nc_2, sigma_Nc_2**2 )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There isn't any cloud at this
       ! grid level.  The value of mu_Nc_n should be -inf.  It will be set to
       ! -huge for purposes of assigning it a value.
       mu_Nc_2_n = -huge( mu_Nc_2_n )
    endif

    !!! Calculate the normalized standard deviation of variables that have
    !!! an assumed lognormal distribution, given the mean and variance of
    !!! those variables.

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       sigma_rr_1_n = stdev_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_rr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       sigma_rr_2_n = stdev_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_rr_2_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       sigma_Nr_1_n = stdev_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_Nr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       sigma_Nr_2_n = stdev_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_Nr_2_n = zero
    endif

    ! Normalized standard deviation of cloud droplet concentration
    ! in PDF component 1.
    if ( Ncm > Nc_tol ) then
       sigma_Nc_1_n = stdev_L2N( mu_Nc_1, sigma_Nc_1**2 )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The standard deviation is simply 0 since cloud droplet
       ! concentration does not vary at this grid level.
       sigma_Nc_1_n = zero
    endif

    ! Normalized standard deviation of cloud droplet concentration
    ! in PDF component 2.
    if ( Ncm > Nc_tol ) then
       sigma_Nc_2_n = stdev_L2N( mu_Nc_2, sigma_Nc_2**2 )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The standard deviation is simply 0 since cloud droplet
       ! concentration does not vary at this grid level.
       sigma_Nc_2_n = zero
    endif

    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_srr_1_n = corr_NL2NN( corr_srr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_srr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_srr_2_n = corr_NL2NN( corr_srr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_srr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_sNr_1_n = corr_NL2NN( corr_sNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_sNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_sNr_2_n = corr_NL2NN( corr_sNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_sNr_2_n = zero
    endif

    ! Normalize the correlation between s and N_c in PDF component 1.
    if ( Ncm > Nc_tol ) then
       corr_sNc_1_n = corr_NL2NN( corr_sNc_1, sigma_Nc_1_n )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The correlations involving cloud droplet concentration are
       ! 0 since cloud droplet concentration does not vary at this grid level.
       corr_sNc_1_n = zero
    endif

    ! Normalize the correlation between s and N_c in PDF component 2.
    if ( Ncm > Nc_tol ) then
       corr_sNc_2_n = corr_NL2NN( corr_sNc_2, sigma_Nc_2_n )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The correlations involving cloud droplet concentration are
       ! 0 since cloud droplet concentration does not vary at this grid level.
       corr_sNc_2_n = zero
    endif

    !!! Calculate the normalized correlation between two variables that both
    !!! have an assumed lognormal distribution, given their correlation and both
    !!! of their normalized standard deviations.

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 1.
    if ( rr1 > rr_tol .and. Nr1 > Nr_tol ) then
       corr_rrNr_1_n = corr_LL2NN( corr_rrNr_1, sigma_rr_1_n, sigma_Nr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 and (or) mean rain drop
       ! concentration in PDF component 1 are (is) less than their (its)
       ! respective tolerance amount(s), and are (is) considered to have a value
       ! of 0.  There is not any rain at this grid level.  The mean in-precip
       ! rain water mixing ratio (1st PDF component) and (or) mean in-precip
       ! rain drop concentration (1st PDF component) are also considered to have
       ! a value of 0.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 2.
    if ( rr2 > rr_tol .and. Nr2 > Nr_tol ) then
       corr_rrNr_2_n = corr_LL2NN( corr_rrNr_2, sigma_rr_2_n, sigma_Nr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 and (or) mean rain drop
       ! concentration in PDF component 2 are (is) less than their (its)
       ! respective tolerance amount(s), and are (is) considered to have a value
       ! of 0.  There is not any rain at this grid level.  The mean in-precip
       ! rain water mixing ratio (2nd PDF component) and (or) mean in-precip
       ! rain drop concentration (2nd PDF component) are also considered to have
       ! a value of 0.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_2_n = zero
    endif


    return

  end subroutine KK_upscaled_setup

  !=============================================================================
  subroutine KK_upscaled_means_driver( rrainm, Nrm, Ncm, &
                                       mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                       mu_Nr_1, mu_Nr_2, mu_Nr_1_n, &
                                       mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                                       sigma_s_1, sigma_s_2, &
                                       sigma_rr_1_n, sigma_rr_2_n, &
                                       sigma_Nr_1_n, sigma_Nr_2_n, &
                                       sigma_Nc_1_n, sigma_Nc_2_n, &
                                       corr_srr_1_n, corr_srr_2_n, &
                                       corr_sNr_1_n, corr_sNr_2_n, &
                                       corr_sNc_1_n, corr_sNc_2_n, &
                                       corr_rrNr_1_n, corr_rrNr_2_n, &
                                       mixt_frac, precip_frac_1, &
                                       precip_frac_2, Nc0_in_cloud, &
                                       l_const_Nc_in_cloud, &
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_mvr_coef, &
                                       KK_evap_tndcy, KK_auto_tndcy, &
                                       KK_accr_tndcy, KK_mean_vol_rad )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rr_tol, & ! Constant(s)
        Nr_tol, & 
        Nc_tol, &
        zero

    use KK_upscaled_means, only: & 
        KK_evap_upscaled_mean, & ! Procedure(s)
        KK_auto_upscaled_mean, &
        KK_accr_upscaled_mean, &
        KK_mvr_upscaled_mean

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rrainm, & ! Mean rain water mixing ratio        [kg/kg]
      Nrm,    & ! Mean rain drop concentration        [num/kg]
      Ncm       ! Mean cloud droplet concentration    [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF comp.) in-precip (ip)[ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF comp.) ip            [ln(kg/kg)]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Nc_1_n,     & ! Mean of ln Nc (1st PDF component)          [ln(num/kg)]
      mu_Nc_2_n,     & ! Mean of ln Nc (2nd PDF component)          [ln(num/kg)]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Nc_1_n,  & ! Standard dev. of ln Nc (1st PDF comp.)     [ln(num/kg)]
      sigma_Nc_2_n,  & ! Standard dev. of ln Nc (2nd PDF comp.)     [ln(num/kg)]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF comp.) ip  [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF comp.) ip  [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      Nc0_in_cloud     ! Constant in-cloud value of cloud droplet conc. [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_coef, & ! KK evaporation coefficient          [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient       [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient            [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient   [m]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      KK_evap_tndcy,   & ! KK evaporation tendency          [(kg/kg)/s]
      KK_auto_tndcy,   & ! KK autoconversion tendency       [(kg/kg)/s]
      KK_accr_tndcy,   & ! KK accretion tendency            [(kg/kg)/s]
      KK_mean_vol_rad    ! KK rain drop mean volume radius  [m]


    !!! Calculate the upscaled KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       KK_evap_tndcy  &
       = KK_evap_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                mu_Nr_1_n, mu_Nr_2_n, sigma_s_1, &
                                sigma_s_2, sigma_rr_1_n, sigma_rr_2_n, &
                                sigma_Nr_1_n, sigma_Nr_2_n, corr_srr_1_n, &
                                corr_srr_2_n, corr_sNr_1_n, corr_sNr_2_n, &
                                corr_rrNr_1_n, corr_rrNr_2_n, KK_evap_coef, &
                                mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       KK_evap_tndcy = zero

    endif

    !!! Calculate the upscaled KK autoconversion tendency.
    if ( Ncm > Nc_tol ) then

       KK_auto_tndcy  &
       = KK_auto_upscaled_mean( mu_s_1, mu_s_2, mu_Nc_1_n, mu_Nc_2_n, &
                                sigma_s_1, sigma_s_2, sigma_Nc_1_n, &
                                sigma_Nc_2_n, corr_sNc_1_n, corr_sNc_2_n, &
                                KK_auto_coef, mixt_frac, &
                                Nc0_in_cloud, l_const_Nc_in_cloud )

    else  ! N_c = 0.

       KK_auto_tndcy = zero

    endif

    !!! Calculate the upscaled KK accretion tendency.
    if ( rrainm > rr_tol ) then

       KK_accr_tndcy  &
       = KK_accr_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                sigma_s_1, sigma_s_2, sigma_rr_1_n, &
                                sigma_rr_2_n, corr_srr_1_n, corr_srr_2_n, &
                                KK_accr_coef, mixt_frac, precip_frac_1, &
                                precip_frac_2 )

    else  ! r_r = 0.

       KK_accr_tndcy = zero

    endif

    !!! Calculate the upscaled KK rain drop mean volume radius.
    if ( rrainm > rr_tol ) then

       KK_mean_vol_rad &
       = KK_mvr_upscaled_mean( mu_rr_1_n, mu_rr_2_n, mu_Nr_1, mu_Nr_2, &
                               mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1_n, &
                               sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                               corr_rrNr_1_n, corr_rrNr_2_n, KK_mvr_coef, &
                               mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r = 0.

       KK_mean_vol_rad = zero

    endif


    return

  end subroutine KK_upscaled_means_driver

  !=============================================================================
  subroutine KK_upscaled_covar_driver( w_mean, exner, rcm, &
                                       rrainm, Nrm, Ncm, &
                                       mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                       mu_Nc_n, sigma_s_1, sigma_s_2, &
                                       sigma_rr_n, sigma_Nr_n, sigma_Nc_n, &
                                       corr_srr_1_n, corr_srr_2_n, &
                                       corr_sNr_1_n, corr_sNr_2_n, &
                                       corr_sNc_1_n, corr_sNc_2_n, &
                                       corr_rrNr_n,  mixt_frac, &
                                       precip_frac, Nc0_in_cloud, &
                                       l_const_Nc_in_cloud, &
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_evap_tndcy, &
                                       KK_auto_tndcy, KK_accr_tndcy, &
                                       pdf_params, level, &
                                       l_stats_samp, &
                                       wprtp_mc_src_tndcy, &
                                       wpthlp_mc_src_tndcy, &
                                       rtp2_mc_src_tndcy, &
                                       thlp2_mc_src_tndcy, &
                                       rtpthlp_mc_src_tndcy )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        two,    & ! Constant(s)
        zero,   &
        Lv,     &
        Cp,     &
        w_tol,  &
        rc_tol, &
        rr_tol, & 
        Nr_tol, & 
        Nc_tol

    use constants_clubb, only:  &
        t_tol => t_mellor_tol  ! Constant

    use KK_upscaled_covariances, only: &
        covar_x_KK_evap,   & ! Procedure(s)
        covar_rt_KK_evap,  &
        covar_thl_KK_evap, &
        covar_x_KK_auto,   &
        covar_rt_KK_auto,  &
        covar_thl_KK_auto, &
        covar_x_KK_accr,   &
        covar_rt_KK_accr,  &
        covar_thl_KK_accr

    use KK_utilities, only: &
        corr_NL2NN  ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use KK_fixed_correlations, only: &
        corr_wrr_NL_cloud, & ! Variable(s)
        corr_wNr_NL_cloud, &
        corr_wNc_NL_cloud, &
        corr_wrr_NL_below, &
        corr_wNr_NL_below, &
        corr_wNc_NL_below, &
        corr_trr_NL_cloud, &
        corr_tNr_NL_cloud, &
        corr_tNc_NL_cloud, &
        corr_trr_NL_below, &
        corr_tNr_NL_below, &
        corr_tNc_NL_below, &
        corr_st_NN_cloud,  &
        corr_st_NN_below

    use parameters_microphys, only: &
        l_fix_s_t_correlations ! Variable(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type, only: & 
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        zt,                    & ! Variable(s)
        iw_KK_evap_covar_zt,   &
        irt_KK_evap_covar_zt,  &
        ithl_KK_evap_covar_zt, &
        iw_KK_auto_covar_zt,   &
        irt_KK_auto_covar_zt,  &
        ithl_KK_auto_covar_zt, &
        iw_KK_accr_covar_zt,   &
        irt_KK_accr_covar_zt,  &
        ithl_KK_accr_covar_zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      w_mean, & ! Mean vertical velocity              [m/s]
      exner,  & ! Exner function                      [-]
      rcm,    & ! Mean cloud water mixing ratio       [kg/kg]
      rrainm, & ! Mean rain water mixing ratio        [kg/kg]
      Nrm,    & ! Mean rain drop concentration        [num/kg]
      Ncm       ! Mean cloud droplet concentration    [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      mu_s_1,       & ! Mean of s (1st PDF component)                    [kg/kg]
      mu_s_2,       & ! Mean of s (2nd PDF component)                    [kg/kg]
      mu_rr_n,      & ! Mean of ln rr (both components)              [ln(kg/kg)]
      mu_Nr_n,      & ! Mean of ln Nr (both components)             [ln(num/kg)]
      mu_Nc_n,      & ! Mean of ln Nc (both components)             [ln(num/kg)]
      sigma_s_1,    & ! Standard deviation of s (1st PDF component)      [kg/kg]
      sigma_s_2,    & ! Standard deviation of s (2nd PDF component)      [kg/kg]
      sigma_rr_n,   & ! Standard deviation of ln rr (both comps.)    [ln(kg/kg)]
      sigma_Nr_n,   & ! Standard deviation of ln Nr (both comps.)   [ln(num/kg)]
      sigma_Nc_n,   & ! Standard deviation of ln Nc (both comps.)   [ln(num/kg)]
      corr_srr_1_n, & ! Correlation between s and ln rr (1st PDF component)  [-]
      corr_srr_2_n, & ! Correlation between s and ln rr (2nd PDF component)  [-]
      corr_sNr_1_n, & ! Correlation between s and ln Nr (1st PDF component)  [-]
      corr_sNr_2_n, & ! Correlation between s and ln Nr (2nd PDF component)  [-]
      corr_sNc_1_n, & ! Correlation between s and ln Nc (1st PDF component)  [-]
      corr_sNc_2_n, & ! Correlation between s and ln Nc (2nd PDF component)  [-]
      corr_rrNr_n,  & ! Correlation between ln rr & ln Nr (both components)  [-]
      mixt_frac,    & ! Mixture fraction                                     [-]
      precip_frac,  & ! Precipitation fraction                               [-]
      Nc0_in_cloud    ! Constant in-cloud value of cloud droplet conc.  [num/kg]

    logical, intent(in) :: &
      l_const_Nc_in_cloud  ! Flag to use a constant value of N_c within cloud

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_coef, & ! KK evaporation coefficient          [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient       [(kg/kg)/s]
      KK_accr_coef    ! KK accretion coefficient            [(kg/kg)/s]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy, & ! KK evaporation tendency            [(kg/kg)/s]
      KK_auto_tndcy, & ! KK autoconversion tendency         [(kg/kg)/s]
      KK_accr_tndcy    ! KK accretion tendency              [(kg/kg)/s]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                        [units vary]

    integer, intent(in) :: &
      level         ! Vertical level index                  [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      wprtp_mc_src_tndcy,   & ! Microphysics tendency for w'rt'  [m*(kg/kg)/s^2]
      wpthlp_mc_src_tndcy,  & ! Microphysics tendency for w'thl' [m*K/s^2]
      rtp2_mc_src_tndcy,    & ! Microphysics tendency for rt'^2  [(kg/kg)^2/s]
      thlp2_mc_src_tndcy,   & ! Microphysics tendency for thl'^2 [K^2/s]
      rtpthlp_mc_src_tndcy    ! Microphysics tend. for rt'thl'   [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      mu_w_1,     & ! Mean of w (1st PDF component)                    [m/s]
      mu_w_2,     & ! Mean of w (2nd PDF component)                    [m/s]
      sigma_w_1,  & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2,  & ! Standard deviation of w (2nd PDF component)      [m/s]
      sigma_t_1,  & ! Standard deviation of t (1st PDF component)      [kg/kg]
      sigma_t_2,  & ! Standard deviation of t (2nd PDF component)      [kg/kg]
      corr_wrr_1, & ! Correlation between w and rr (1st PDF component) [-]
      corr_wrr_2, & ! Correlation between w and rr (2nd PDF component) [-] 
      corr_wNr_1, & ! Correlation between w and Nr (1st PDF component) [-] 
      corr_wNr_2, & ! Correlation between w and Nr (2nd PDF component) [-] 
      corr_wNc_1, & ! Correlation between w and Nc (1st PDF component) [-] 
      corr_wNc_2, & ! Correlation between w and Nc (2nd PDF component) [-] 
      corr_ts_1,  & ! Correlation between t and s (1st PDF component)  [-]
      corr_ts_2,  & ! Correlation between t and s (2nd PDF component)  [-]
      corr_trr_1, & ! Correlation between t and rr (1st PDF component) [-] 
      corr_trr_2, & ! Correlation between t and rr (2nd PDF component) [-]
      corr_tNr_1, & ! Correlation between t and Nr (1st PDF component) [-]
      corr_tNr_2, & ! Correlation between t and Nr (2nd PDF component) [-]
      corr_tNc_1, & ! Correlation between t and Nc (1st PDF component) [-]
      corr_tNc_2    ! Correlation between t and Nc (2nd PDF component) [-]

    real( kind = core_rknd ) :: &
      corr_wrr_1_n, & ! Correlation between w and ln rr (1st PDF component) [-]
      corr_wrr_2_n, & ! Correlation between w and ln rr (2nd PDF component) [-]
      corr_wNr_1_n, & ! Correlation between w and ln Nr (1st PDF component) [-]
      corr_wNr_2_n, & ! Correlation between w and ln Nr (2nd PDF component) [-]
      corr_wNc_1_n, & ! Correlation between w and ln Nc (1st PDF component) [-]
      corr_wNc_2_n, & ! Correlation between w and ln Nc (2nd PDF component) [-]
      corr_trr_1_n, & ! Correlation between t and ln rr (1st PDF component) [-]
      corr_trr_2_n, & ! Correlation between t and ln rr (2nd PDF component) [-]
      corr_tNr_1_n, & ! Correlation between t and ln Nr (1st PDF component) [-]
      corr_tNr_2_n, & ! Correlation between t and ln Nr (2nd PDF component) [-]
      corr_tNc_1_n, & ! Correlation between t and ln Nc (1st PDF component) [-]
      corr_tNc_2_n    ! Correlation between t and ln Nc (2nd PDF component) [-]

    real( kind = core_rknd ) :: &
      crt1,  & ! Coefficient c_rt (1st PDF component)    [-]
      crt2,  & ! Coefficient c_rt (2nd PDF component)    [-]
      cthl1, & ! Coefficient c_thl (1st PDF component)   [(kg/kg)/K]
      cthl2    ! Coefficient c_thl (2nd PDF component)   [(kg/kg)/K]

    real( kind = core_rknd ) :: &
      w_KK_evap_covar,   & ! Covar. btw. w and KK evap. tend.    [m*(kg/kg)/s^2]
      rt_KK_evap_covar,  & ! Covar. btw. rt and KK evap. tend.   [(kg/kg)^2/s]
      thl_KK_evap_covar, & ! Covar. btw. thl and KK evap. tend.  [K*(kg/kg)/s]
      w_KK_auto_covar,   & ! Covar. btw. w and KK auto. tend.    [m*(kg/kg)/s^2]
      rt_KK_auto_covar,  & ! Covar. btw. rt and KK auto. tend.   [(kg/kg)^2/s]
      thl_KK_auto_covar, & ! Covar. btw. thl and KK auto. tend.  [K*(kg/kg)/s]
      w_KK_accr_covar,   & ! Covar. btw. w and KK accr. tend.    [m*(kg/kg)/s^2]
      rt_KK_accr_covar,  & ! Covar. btw. rt and KK accr. tend.   [(kg/kg)^2/s]
      thl_KK_accr_covar    ! Covar. btw. thl and KK accr. tend.  [K*(kg/kg)/s]

    ! Constant Parameters
    !
    ! Set the correlations between vertical velocity (w) and extended liquid
    ! water mixing ratio (s_*) to 0.
    ! Note:  The component correlations of w and r_t and the component
    !        correlations of w and theta_l are both set to be 0 within the CLUBB
    !        model code.  In other words, w and r_t (theta_l) have overall
    !        covariance w'r_t' (w'theta_l'), but the single component covariance
    !        and correlation are defined to be 0.  Likewise, the single
    !        component correlation and covariance of w and s, as well as
    !        w and t, are defined to be 0.
    real( kind = core_rknd ), parameter :: &
      corr_ws_1 = zero, & ! Correlation between w and s (1st PDF component) [-]
      corr_ws_2 = zero    ! Correlation between w and s (2nd PDF component) [-]
    !
    ! Set the component mean values of t_* to 0.
    ! Note:  The component mean values of t_* are not important.  They can be
    !        set to anything.  They cancel out in the model code.  However, the
    !        best thing to do is to set them to 0 and avoid any kind of
    !        numerical error.
    real( kind = core_rknd ), parameter :: &
      mu_t_1 = zero, & ! Mean of t (1st PDF component)  [kg/kg]
      mu_t_2 = zero    ! Mean of t (2nd PDF component)  [kg/kg]


    ! Enter the PDF parameters.
    mu_w_1    = pdf_params%w1
    mu_w_2    = pdf_params%w2
    sigma_w_1 = sqrt( pdf_params%varnce_w1 )
    sigma_w_2 = sqrt( pdf_params%varnce_w2 )
    sigma_t_1 = pdf_params%stdev_t1
    sigma_t_2 = pdf_params%stdev_t2

    if ( l_fix_s_t_correlations ) then
       if ( mu_s_1 > zero ) then
          corr_ts_1 = corr_st_NN_cloud
       else
          corr_ts_1 = corr_st_NN_below
       endif
       if ( mu_s_2 > zero ) then
          corr_ts_2 = corr_st_NN_cloud
       else
          corr_ts_2 = corr_st_NN_below
       endif
    else
       corr_ts_1 = pdf_params%corr_st_1
       corr_ts_2 = pdf_params%corr_st_2
    endif

    crt1      = pdf_params%crt1
    crt2      = pdf_params%crt2
    cthl1     = pdf_params%cthl1
    cthl2     = pdf_params%cthl2

    
    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    !
    ! Set the value of the parameters based on whether the altitude is above or
    ! below cloud base.  Determine whether there is cloud at any given vertical
    ! level.  In order for a vertical level to have cloud, the amount of cloud
    ! water (rcm) must be greater than or equal to the tolerance level (rc_tol).
    ! If there is cloud at a given vertical level, then the ###_cloud value is
    ! used.  Otherwise, the ###_below value is used.
    if ( rcm > rc_tol ) then
       corr_wrr_1 = corr_wrr_NL_cloud
       corr_wrr_2 = corr_wrr_NL_cloud
       corr_wNr_1 = corr_wNr_NL_cloud
       corr_wNr_2 = corr_wNr_NL_cloud
       corr_wNc_1 = corr_wNc_NL_cloud
       corr_wNc_2 = corr_wNc_NL_cloud
       corr_trr_1 = corr_trr_NL_cloud
       corr_trr_2 = corr_trr_NL_cloud
       corr_tNr_1 = corr_tNr_NL_cloud
       corr_tNr_2 = corr_tNr_NL_cloud
       corr_tNc_1 = corr_tNc_NL_cloud
       corr_tNc_2 = corr_tNc_NL_cloud
    else
       corr_wrr_1 = corr_wrr_NL_below
       corr_wrr_2 = corr_wrr_NL_below
       corr_wNr_1 = corr_wNr_NL_below
       corr_wNr_2 = corr_wNr_NL_below
       corr_wNc_1 = corr_wNc_NL_below
       corr_wNc_2 = corr_wNc_NL_below
       corr_trr_1 = corr_trr_NL_below
       corr_trr_2 = corr_trr_NL_below
       corr_tNr_1 = corr_tNr_NL_below
       corr_tNr_2 = corr_tNr_NL_below
       corr_tNc_1 = corr_tNc_NL_below
       corr_tNc_2 = corr_tNc_NL_below
    endif

    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! (single) lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    if ( rrainm > rr_tol ) then

       ! Normalize the correlation between w and r_r in PDF component 1.
       corr_wrr_1_n = corr_NL2NN( corr_wrr_1, sigma_rr_n )

       ! Normalize the correlation between w and r_r in PDF component 2.
       corr_wrr_2_n = corr_NL2NN( corr_wrr_2, sigma_rr_n )

       ! Normalize the correlation between t and r_r in PDF component 1.
       corr_trr_1_n = corr_NL2NN( corr_trr_1, sigma_rr_n )

       ! Normalize the correlation between t and r_r in PDF component 2.
       corr_trr_2_n = corr_NL2NN( corr_trr_2, sigma_rr_n )

    else

       ! Mean rain water mixing ratio is less than the tolerance amount.  It is
       ! considered to have a value of 0.  There is not any rain at this grid
       ! level.  The correlations involving rain water mixing ratio are 0 since
       ! rain water mixing ratio does not vary at this grid level.
       corr_wrr_1_n = zero
       corr_wrr_2_n = zero
       corr_trr_1_n = zero
       corr_trr_2_n = zero

    endif

    if ( Nrm > Nr_tol ) then

       ! Normalize the correlation between w and N_r in PDF component 1.
       corr_wNr_1_n = corr_NL2NN( corr_wNr_1, sigma_Nr_n )

       ! Normalize the correlation between w and N_r in PDF component 2.
       corr_wNr_2_n = corr_NL2NN( corr_wNr_2, sigma_Nr_n )

       ! Normalize the correlation between t and N_r in PDF component 1.
       corr_tNr_1_n = corr_NL2NN( corr_tNr_1, sigma_Nr_n )

       ! Normalize the correlation between t and N_r in PDF component 2.
       corr_tNr_2_n = corr_NL2NN( corr_tNr_2, sigma_Nr_n )

    else

       ! Mean rain drop concentration is less than the tolerance amount.  It is
       ! considered to have a value of 0.  There is not any rain at this grid
       ! level.  The correlations involving rain drop concentration are 0 since
       ! rain drop concentration does not vary at this grid level.
       corr_wNr_1_n = zero
       corr_wNr_2_n = zero
       corr_tNr_1_n = zero
       corr_tNr_2_n = zero

    endif

    if ( Ncm > Nc_tol ) then

       ! Normalize the correlation between s and N_c in PDF component 1.
       corr_wNc_1_n = corr_NL2NN( corr_wNc_1, sigma_Nc_n )

       ! Normalize the correlation between s and N_c in PDF component 2.
       corr_wNc_2_n = corr_NL2NN( corr_wNc_2, sigma_Nc_n )

       ! Normalize the correlation between t and N_c in PDF component 1.
       corr_tNc_1_n = corr_NL2NN( corr_tNc_1, sigma_Nc_n )

       ! Normalize the correlation between t and N_c in PDF component 2.
       corr_tNc_2_n = corr_NL2NN( corr_tNc_2, sigma_Nc_n )

    else

       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There is not any cloud at this
       ! grid level.  The correlations involving cloud droplet concentration are
       ! 0 since cloud droplet concentration does not vary at this grid level.
       corr_wNc_1_n = zero
       corr_wNc_2_n = zero
       corr_tNc_1_n = zero
       corr_tNc_2_n = zero

    endif


    ! Calculate the covariance of vertical velocity and KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       w_KK_evap_covar  &
       = covar_x_KK_evap( mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_rr_n, &
                          mu_Nr_n, sigma_w_1, sigma_w_2, sigma_s_1, &
                          sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                          corr_ws_1, corr_ws_2, corr_wrr_1_n, &
                          corr_wrr_2_n, corr_wNr_1_n, corr_wNr_2_n, &
                          corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                          corr_sNr_2_n, corr_rrNr_n, w_mean, &
                          KK_evap_tndcy, KK_evap_coef, w_tol, &
                          mixt_frac, precip_frac )

    else  ! r_r or N_r = 0.

       w_KK_evap_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK evaporation
    ! tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       rt_KK_evap_covar  &
       = covar_rt_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                           mu_Nr_n, sigma_t_1, sigma_t_2, sigma_s_1, &
                           sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                           corr_ts_1, corr_ts_2, corr_trr_1_n, &
                           corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                           corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                           corr_sNr_2_n, corr_rrNr_n, KK_evap_tndcy, &
                           KK_evap_coef, t_tol, crt1, crt2, &
                           mixt_frac, precip_frac )

    else  ! r_r or N_r = 0.

       rt_KK_evap_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       thl_KK_evap_covar  &
       = covar_thl_KK_evap( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                            mu_Nr_n, sigma_t_1, sigma_t_2, sigma_s_1, &
                            sigma_s_2, sigma_rr_n, sigma_Nr_n, &
                            corr_ts_1, corr_ts_2, corr_trr_1_n, &
                            corr_trr_2_n, corr_tNr_1_n, corr_tNr_2_n, &
                            corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                            corr_sNr_2_n, corr_rrNr_n, KK_evap_tndcy, &
                            KK_evap_coef, t_tol, cthl1, cthl2, &
                            mixt_frac, precip_frac )

    else  ! r_r or N_r = 0.

       thl_KK_evap_covar = zero

    endif

    ! Calculate the covariance of vertical velocity and KK autoconversion
    ! tendency.
    if ( Ncm > Nc_tol ) then

       w_KK_auto_covar  &
       = covar_x_KK_auto( mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_Nc_n, &
                          sigma_w_1, sigma_w_2, sigma_s_1, sigma_s_2, &
                          sigma_Nc_n, corr_ws_1, corr_ws_2, &
                          corr_wNc_1_n, corr_wNc_2_n, corr_sNc_1_n, &
                          corr_sNc_2_n, w_mean, KK_auto_tndcy, &
                          KK_auto_coef, w_tol, mixt_frac, &
                          Nc0_in_cloud, l_const_Nc_in_cloud )

    else  ! N_c = 0.

       w_KK_auto_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK autoconversion
    ! tendency.
    if ( Ncm > Nc_tol ) then

       rt_KK_auto_covar  &
       = covar_rt_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Nc_n, &
                           sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                           sigma_Nc_n, corr_ts_1, corr_ts_2, &
                           corr_tNc_1_n, corr_tNc_2_n, corr_sNc_1_n, &
                           corr_sNc_2_n, KK_auto_tndcy, KK_auto_coef, &
                           t_tol, crt1, crt2, mixt_frac, &
                           Nc0_in_cloud, l_const_Nc_in_cloud )

    else  ! N_c = 0.

       rt_KK_auto_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK autoconversion tendency.
    if ( Ncm > Nc_tol ) then

       thl_KK_auto_covar  &
       = covar_thl_KK_auto( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_Nc_n, &
                            sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                            sigma_Nc_n, corr_ts_1, corr_ts_2, &
                            corr_tNc_1_n, corr_tNc_2_n, corr_sNc_1_n, &
                            corr_sNc_2_n, KK_auto_tndcy, KK_auto_coef, &
                            t_tol, cthl1, cthl2, mixt_frac, &
                            Nc0_in_cloud, l_const_Nc_in_cloud )

    else  ! N_c = 0.

       thl_KK_auto_covar = zero

    endif

    ! Calculate the covariance of vertical velocity and KK accretion tendency.
    if ( rrainm > rr_tol ) then

       w_KK_accr_covar  &
       = covar_x_KK_accr( mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_rr_n, &
                          sigma_w_1, sigma_w_2, sigma_s_1, sigma_s_2, &
                          sigma_rr_n, corr_ws_1, corr_ws_2, &
                          corr_wrr_1_n, corr_wrr_2_n, corr_srr_1_n, &
                          corr_srr_2_n, w_mean, KK_accr_tndcy, &
                          KK_accr_coef, w_tol, mixt_frac, precip_frac )

    else  ! r_r = 0.

       w_KK_accr_covar = zero

    endif

    ! Calculate the covariance of total water mixing ratio and KK accretion
    ! tendency.
    if ( rrainm > rr_tol ) then

       rt_KK_accr_covar  &
       = covar_rt_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                           sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                           sigma_rr_n, corr_ts_1, corr_ts_2, &
                           corr_trr_1_n, corr_trr_2_n, corr_srr_1_n, &
                           corr_srr_2_n, KK_accr_tndcy, KK_accr_coef, &
                           t_tol, crt1, crt2, mixt_frac, precip_frac )

    else  ! r_r = 0.

       rt_KK_accr_covar = zero

    endif

    ! Calculate the covariance of liquid water potential temperature and
    ! KK accretion tendency.
    if ( rrainm > rr_tol ) then

       thl_KK_accr_covar  &
       = covar_thl_KK_accr( mu_t_1, mu_t_2, mu_s_1, mu_s_2, mu_rr_n, &
                            sigma_t_1, sigma_t_2, sigma_s_1, sigma_s_2, &
                            sigma_rr_n, corr_ts_1, corr_ts_2, &
                            corr_trr_1_n, corr_trr_2_n, corr_srr_1_n, &
                            corr_srr_2_n, KK_accr_tndcy, KK_accr_coef, &
                            t_tol, cthl1, cthl2, mixt_frac, precip_frac )

    else  ! r_r = 0.

       thl_KK_accr_covar = zero

    endif


    ! Statistics
    if ( l_stats_samp ) then
       ! All of these covariance variables are being calculated on thermodynamic
       ! grid levels (all inputs are on thermodynamic grid levels, so the output
       ! is also on thermodynamic grid levels).  These covariances will be
       ! combined in various ways to produce the microphysics tendency terms for
       ! various model predictive variances and covariances.  These source
       ! tendency terms will be interpolated to momentum grid levels.

       ! Covariance of w and KK evaporation tendency.
       if ( iw_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_evap_covar_zt, level, &
                                   w_KK_evap_covar, zt )
       endif

       ! Covariance of r_t and KK evaporation tendency.
       if ( irt_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_evap_covar_zt, level, &
                                   rt_KK_evap_covar, zt )
       endif

       ! Covariance of theta_l and KK evaporation tendency.
       if ( ithl_KK_evap_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_evap_covar_zt, level, &
                                   thl_KK_evap_covar, zt )
       endif

       ! Covariance of w and KK autoconversion tendency.
       if ( iw_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_auto_covar_zt, level, &
                                   w_KK_auto_covar, zt )
       endif

       ! Covariance of r_t and KK autoconversion tendency.
       if ( irt_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_auto_covar_zt, level, &
                                   rt_KK_auto_covar, zt )
       endif

       ! Covariance of theta_l and KK autoconversion tendency.
       if ( ithl_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_auto_covar_zt, level, &
                                   thl_KK_auto_covar, zt )
       endif

       ! Covariance of w and KK accretion tendency.
       if ( iw_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( iw_KK_accr_covar_zt, level, &
                                   w_KK_accr_covar, zt )
       endif

       ! Covariance of r_t and KK accretion tendency.
       if ( irt_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( irt_KK_accr_covar_zt, level, &
                                   rt_KK_accr_covar, zt )
       endif

       ! Covariance of theta_l and KK accretion tendency.
       if ( ithl_KK_auto_covar_zt > 0 ) then
          call stat_update_var_pt( ithl_KK_accr_covar_zt, level, &
                                   thl_KK_accr_covar, zt )
       endif

    endif ! l_stats_samp


    ! Calculate the microphysics tendency for <w'r_t'>.
    wprtp_mc_src_tndcy  &
    = - ( w_KK_auto_covar + w_KK_accr_covar + w_KK_evap_covar )

    ! Calculate the microphysics tendency for <w'th_l'>.
    wpthlp_mc_src_tndcy  &
    = ( Lv / ( Cp * exner ) )  &
      * ( w_KK_auto_covar + w_KK_accr_covar + w_KK_evap_covar )

    ! Calculate the microphysics tendency for <r_t'^2>.
    rtp2_mc_src_tndcy  &
    = - two * ( rt_KK_auto_covar + rt_KK_accr_covar + rt_KK_evap_covar )

    ! Calculate the microphysics tendency for <th_l'^2>.
    thlp2_mc_src_tndcy  &
    = two * ( Lv / ( Cp * exner ) )  &
          * ( thl_KK_auto_covar + thl_KK_accr_covar + thl_KK_evap_covar )

    ! Calculate the microphysics tendency for <r_t'th_l'>.
    rtpthlp_mc_src_tndcy  &
    = ( Lv / ( Cp * exner ) )  &
      * ( rt_KK_auto_covar + rt_KK_accr_covar + rt_KK_evap_covar )  &
      - ( thl_KK_auto_covar + thl_KK_accr_covar + thl_KK_evap_covar )


    return

  end subroutine KK_upscaled_covar_driver

  !=============================================================================
  subroutine KK_microphys_adjust( dt, exner, rcm, rrainm, Nrm, &
                                  KK_evap_tndcy, KK_auto_tndcy, &
                                  KK_accr_tndcy, KK_Nrm_evap_tndcy, &
                                  KK_Nrm_auto_tndcy, l_src_adj_enabled, &
                                  l_evap_adj_enabled, l_stats_samp, level, &
                                  rrainm_mc_tndcy, Nrm_mc_tndcy, &
                                  rvm_mc, rcm_mc, thlm_mc )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd, &  ! Variable(s)
        time_precision

    use constants_clubb, only: &
        Lv,           & ! Constant(s)
        Cp,           &
        zero,         &
        rr_tol,       &
        Nr_tol

    use KK_Nrm_tendencies, only: &
        KK_Nrm_evap_local_mean, & ! Procedure(s)
        KK_Nrm_auto_mean

    use stats_variables, only: &
        irrainm_src_adj,  & ! Variable(s)
        irrainm_cond_adj, &
        iNrm_src_adj,     &
        iNrm_cond_adj,    &
        zt

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    implicit none

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt        ! Model time step duration                 [s]

    real( kind = core_rknd ), intent(in) :: &
      exner,  & ! Exner function                           [-]
      rcm,    & ! Mean cloud water mixing ratio            [kg/kg]
      rrainm, & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm       ! Mean rain drop concentration, < N_r >    [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy,     & ! Mean KK <r_r> evaporation tendency     [(kg/kg)/s]
      KK_auto_tndcy,     & ! Mean KK <r_r> autoconversion tendency  [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK <r_r> accretion tendency       [(kg/kg)/s]
      KK_Nrm_evap_tndcy, & ! Mean KK <N_r> evaporation tendency     [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK <N_r> autoconversion tendency  [(num/kg)/s]

    logical, intent(in) :: &
      l_src_adj_enabled,  & ! Flag to enable rrainm/Nrm source adjustment
      l_evap_adj_enabled, & ! Flag to enable rrainm/Nrm evaporation adjustment
      l_stats_samp          ! Flag to sample statistical output

    integer, intent(in) :: & 
      level    ! Vertical level index

    ! Output Variables
    real( kind = core_rknd ), intent(out) ::  &
      rrainm_mc_tndcy, & ! Mean <dr_r/dt> due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean <dN_r/dt> due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio       [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio        [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature    [K/s]

    ! Local Variables
    real( kind = core_rknd ) ::  &
      rrainm_source,   & ! Total source term rate for rrainm        [(kg/kg)/s]
      Nrm_source,      & ! Total source term rate for Nrm           [(num/kg)/s]
      rrainm_src_adj,  & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj,     & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrainm_evap_net, & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net       ! Net evaporation rate of <N_r>            [(num/kg)/s]

    real( kind = core_rknd ) :: &
      rrainm_src_max,    & ! Maximum allowable rrainm source rate    [(kg/kg)/s]
      rrainm_auto_ratio, & ! Ratio of rrainm autoconv to overall source term [-]
      total_rc_needed      ! Amount of r_c needed to over the timestep
                           ! for rain source terms                       [kg/kg]


    !!! Source-adjustment code for rrainm and Nrm.
    rrainm_source = KK_auto_tndcy + KK_accr_tndcy
    Nrm_source = KK_Nrm_auto_tndcy

    if ( l_src_adj_enabled ) then

       ! The increase of rain due to autoconversion and accretion both draw
       ! water from the available cloud water.  Over a long time step, these
       ! rates may over-deplete cloud water.  In other words, these processes
       ! may draw more cloud water than there is available.  Thus, the total
       ! source rate multiplied by the duration of the time step cannot exceed
       ! the total amount of cloud water available.  If it does, then the rate
       ! must be adjusted.
       total_rc_needed = rrainm_source * real( dt, kind = core_rknd )

       if ( total_rc_needed > rcm ) then

          ! The maximum allowable rate of the source terms is rcm/dt.
          rrainm_src_max = rcm / real( dt, kind = core_rknd )

          ! The amount of adjustment to the source terms.
          ! This value should always be negative.
          rrainm_src_adj = rrainm_src_max - rrainm_source

          ! Reset the value of the source terms to the maximum allowable value
          ! of the source terms.
          rrainm_source = rrainm_src_max

          ! The rrainm source terms are made up of autoconversion and accretion.
          ! Only the sum of those two terms is corrected.  However, Nrm has only
          ! an autoconversion term for a source term.  Assume that change in the
          ! rrainm autoconversion term is proportional to to the total rrainm
          ! adjustment rate by the ratio of rrainm autoconversion to the overall
          ! source term.  Then, plug the rrainm autoconversion adjustment into
          ! the equation for Nrm autoconversion to determine the effect on the
          ! Nrm source term.
          rrainm_auto_ratio = KK_auto_tndcy /  &
                              ( KK_auto_tndcy + KK_accr_tndcy )
          Nrm_src_adj = KK_Nrm_auto_mean( rrainm_auto_ratio * rrainm_src_adj )

          ! Change Nrm by Nrm_src_adj.  Nrm_src_adj will always be negative.
          Nrm_source = Nrm_source + Nrm_src_adj

       else ! total_rc_needed <= rcm: enough available cloud water.

          rrainm_src_adj = zero
          Nrm_src_adj    = zero

       endif

    else ! l_src_adj_enabled is false: source adjustment is disabled.

       rrainm_src_adj = zero
       Nrm_src_adj    = zero

    endif


    !!! Evaporation-adjustment code for rrainm and Nrm.
    if ( l_evap_adj_enabled ) then

       ! Prevent over-evaporation of rain over a long model time step.
       ! Limit the evaporation rate of both rain water mixing ratio and rain
       ! drop concentration.  The total amount of rain lost due to evaporation
       ! cannot be so great as to result in negative rain.

       ! Calculate net evaporation rate of <r_r>.
       rrainm_evap_net = max( KK_evap_tndcy, &
                              - rrainm / real( dt, kind = core_rknd ) )

       ! Recalcuate the net evaporation rate of <N_r> based on the net
       ! evaporation rate of <r_r>.
       if ( KK_evap_tndcy /= rrainm_evap_net .and. &
            rrainm > rr_tol .and. Nrm > Nr_tol ) then
          Nrm_evap_net &
          = KK_Nrm_evap_local_mean( rrainm_evap_net, Nrm, rrainm, dt )
       else
          Nrm_evap_net = KK_Nrm_evap_tndcy
       endif

       Nrm_evap_net = max( Nrm_evap_net, &
                           - Nrm / real( dt, kind = core_rknd ) )

    else ! l_evap_adj_enabled is false: evaporation adjustment is disabled.

       rrainm_evap_net = KK_evap_tndcy
       Nrm_evap_net    = KK_Nrm_evap_tndcy

    endif


    !!! Calculate overall KK microphysics tendencies.
    rrainm_mc_tndcy = rrainm_evap_net + rrainm_source
    Nrm_mc_tndcy    = Nrm_evap_net + Nrm_source

    !!! Explicit contributions to thlm and rtm from the microphysics
    rvm_mc  = -rrainm_evap_net
    rcm_mc  = -rrainm_source  ! Accretion + Autoconversion
    thlm_mc = ( Lv / ( Cp * exner ) ) * rrainm_mc_tndcy


    ! Statistics
    if ( l_stats_samp ) then

       if ( irrainm_src_adj > 0 ) then
          call stat_update_var_pt( irrainm_src_adj, level, rrainm_src_adj, zt )
       endif

       if ( iNrm_src_adj > 0 ) then
          call stat_update_var_pt( iNrm_src_adj, level, Nrm_src_adj, zt )
       endif

       if ( irrainm_cond_adj > 0 ) then
          call stat_update_var_pt( irrainm_cond_adj, level, &
                                   rrainm_evap_net - KK_evap_tndcy, zt )
       endif

       if ( irrainm_cond_adj > 0 ) then
          call stat_update_var_pt( iNrm_cond_adj, level, &
                                   Nrm_evap_net - KK_Nrm_evap_tndcy, zt )
       endif

    endif


    return

  end subroutine KK_microphys_adjust

  !=============================================================================
  subroutine KK_upscaled_stats( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                mu_Nc_1, mu_Nc_2, mu_rr_1_n, mu_rr_2_n, &
                                mu_Nr_1_n, mu_Nr_2_n, mu_Nc_1_n, mu_Nc_2_n, &
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                sigma_Nr_2, sigma_Nc_1, sigma_Nc_2, &
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                sigma_Nr_2_n, sigma_Nc_1_n, sigma_Nc_2_n, &
                                corr_srr_1, corr_srr_2, corr_sNr_1, &
                                corr_sNr_2, corr_sNc_1, corr_sNc_2, &
                                corr_rrNr_1, corr_rrNr_2, corr_srr_1_n, &
                                corr_srr_2_n, corr_sNr_1_n, corr_sNr_2_n, &
                                corr_sNc_1_n, corr_sNc_2_n, corr_rrNr_1_n, &
                                corr_rrNr_2_n, mixt_frac, precip_frac_1, &
                                precip_frac_2, KK_mvr_coef, &
                                KK_mean_vol_rad, rrainm, Nrm, level, &
                                l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        rr_tol, & ! Constant(s)
        Nr_tol, &
        zero

    use KK_utilities, only: &
        calc_xp2    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        imu_rr_1,       & ! Variable(s)
        imu_rr_2,       &
        imu_Nr_1,       &
        imu_Nr_2,       &
        imu_Nc_1,       &
        imu_Nc_2,       &
        imu_rr_1_n,     &
        imu_rr_2_n,     &
        imu_Nr_1_n,     &
        imu_Nr_2_n,     &
        imu_Nc_1_n,     &
        imu_Nc_2_n,     &
        isigma_rr_1,    &
        isigma_rr_2,    &
        isigma_Nr_1,    &
        isigma_Nr_2,    &
        isigma_Nc_1,    &
        isigma_Nc_2,    &
        isigma_rr_1_n,  &
        isigma_rr_2_n,  &
        isigma_Nr_1_n,  &
        isigma_Nr_2_n,  &
        isigma_Nc_1_n,  &
        isigma_Nc_2_n,  &
        icorr_srr_1,    &
        icorr_srr_2,    &
        icorr_sNr_1,    &
        icorr_sNr_2,    &
        icorr_sNc_1,    &
        icorr_sNc_2,    &
        icorr_rrNr_1,   &
        icorr_rrNr_2,   &
        icorr_srr_1_n,  &
        icorr_srr_2_n,  &
        icorr_sNr_1_n,  &
        icorr_sNr_2_n,  &
        icorr_sNc_1_n,  &
        icorr_sNc_2_n,  &
        icorr_rrNr_1_n, &
        icorr_rrNr_2_n, &
        zt

    use stats_variables, only : &
        iKK_mvr_variance_zt, &
        irrp2_zt,            &
        iNrp2_zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Nc_1,       & ! Mean of Nc (1st PDF component)                 [num/kg]
      mu_Nc_2,       & ! Mean of Nc (2nd PDF component)                 [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Nc_1_n,     & ! Mean of ln Nc (1st PDF component)          [ln(num/kg)]
      mu_Nc_2_n,     & ! Mean of ln Nc (2nd PDF component)          [ln(num/kg)]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Nc_1,    & ! Standard deviation of Nc (1st PDF component)   [num/kg]
      sigma_Nc_2,    & ! Standard deviation of Nc (2nd PDF component)   [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Nc_1_n,  & ! Standard dev. of ln Nc (1st PDF comp.)     [ln(num/kg)]
      sigma_Nc_2_n,  & ! Standard dev. of ln Nc (2nd PDF comp.)     [ln(num/kg)]
      corr_srr_1,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_srr_2,    & ! Correlation between s and rr (2nd PDF component) ip [-]
      corr_sNr_1,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNr_2,    & ! Correlation between s and Nr (2nd PDF component) ip [-]
      corr_sNc_1,    & ! Correlation between s and Nc (1st PDF component)    [-]
      corr_sNc_2,    & ! Correlation between s and Nc (2nd PDF component)    [-]
      corr_rrNr_1,   & ! Correlation between rr & Nr (1st PDF component) ip  [-]
      corr_rrNr_2,   & ! Correlation between rr & Nr (2nd PDF component) ip  [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNc_1_n,  & ! Correlation between s and ln Nc (1st PDF comp.)     [-]
      corr_sNc_2_n,  & ! Correlation between s and ln Nc (2nd PDF comp.)     [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac,       & ! Mixture fraction                             [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)   [-]
      precip_frac_2,   & ! Precipitation fraction (2nd PDF component)   [-]
      KK_mvr_coef,     & ! KK mean volume radius coefficient            [m]
      KK_mean_vol_rad, & ! KK rain drop mean volume radius              [m]
      rrainm,          & ! Mean rain water mixing ratio                 [kg/kg]
      Nrm                ! Mean rain drop concentration                 [num/kg]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variables
    real( kind = core_rknd ) :: &
      KK_mvr_variance, & ! Variance of KK rain drop mean vol rad   [m^2]
      rrp2,            & ! Overall variance of r_r                 [(kg/kg)^2]
      Nrp2               ! Overall variance of N_r                 [(num/kg)^2]


    !!! Output the statistics for upscaled KK.

    ! Statistics
    if ( l_stats_samp ) then

       ! Mean of in-precip rain drop concentration in PDF component 1.
       if ( imu_rr_1 > 0 ) then
          call stat_update_var_pt( imu_rr_1, level, mu_rr_1, zt )
       endif

       ! Mean of in-precip rain drop concentration in PDF component 2.
       if ( imu_rr_2 > 0 ) then
          call stat_update_var_pt( imu_rr_2, level, mu_rr_2, zt )
       endif

       ! Mean of in-precip rain drop concentration in PDF component 1.
       if ( imu_Nr_1 > 0 ) then
          call stat_update_var_pt( imu_Nr_1, level, mu_Nr_1, zt )
       endif

       ! Mean of in-precip rain drop concentration in PDF component 2.
       if ( imu_Nr_2 > 0 ) then
          call stat_update_var_pt( imu_Nr_2, level, mu_Nr_2, zt )
       endif

       ! Mean of cloud droplet concentration in PDF component 1.
       if ( imu_Nc_1 > 0 ) then
          call stat_update_var_pt( imu_Nc_1, level, mu_Nc_1, zt )
       endif

       ! Mean of cloud droplet concentration in PDF component 2.
       if ( imu_Nc_2 > 0 ) then
          call stat_update_var_pt( imu_Nc_2, level, mu_Nc_2, zt )
       endif

       ! Mean (in-precip) of ln r_r in PDF component 1.
       if ( imu_rr_1_n > 0 ) then
          if ( mu_rr_1_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_rr_1_n, level, mu_rr_1_n, zt )
          else
             ! When rr1 is 0 (or below tolerance value), mu_rr_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_rr_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln r_r in PDF component 2.
       if ( imu_rr_2_n > 0 ) then
          if ( mu_rr_2_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_rr_2_n, level, mu_rr_2_n, zt )
          else
             ! When rr2 is 0 (or below tolerance value), mu_rr_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_rr_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln N_r in PDF component 1.
       if ( imu_Nr_1_n > 0 ) then
          if ( mu_Nr_1_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Nr_1_n, level, mu_Nr_1_n, zt )
          else
             ! When Nr1 is 0 (or below tolerance value), mu_Nr_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nr_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln N_r in PDF component 2.
       if ( imu_Nr_2_n > 0 ) then
          if ( mu_Nr_2_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Nr_2_n, level, mu_Nr_2_n, zt )
          else
             ! When Nr2 is 0 (or below tolerance value), mu_Nr_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nr_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean of ln N_c in PDF component 1.
       if ( imu_Nc_1_n > 0 ) then
          if ( mu_Nc_1_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Nc_1_n, level, mu_Nc_1_n, zt )
          else
             ! When Ncm is 0 (or below tolerance value), mu_Nc_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nc_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean of ln N_c in PDF component 2.
       if ( imu_Nc_2_n > 0 ) then
          if ( mu_Nc_2_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Nc_2_n, level, mu_Nc_2_n, zt )
          else
             ! When Ncm is 0 (or below tolerance value), mu_Nc_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nc_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Standard deviation of in-precip rain water mixing ratio
       ! in PDF component 1.
       if ( isigma_rr_1 > 0 ) then
          call stat_update_var_pt( isigma_rr_1, level, sigma_rr_1, zt )
       endif

       ! Standard deviation of in-precip rain water mixing ratio
       ! in PDF component 2.
       if ( isigma_rr_2 > 0 ) then
          call stat_update_var_pt( isigma_rr_2, level, sigma_rr_2, zt )
       endif

       ! Standard deviation of in-precip rain drop concentration
       ! in PDF component 1.
       if ( isigma_Nr_1 > 0 ) then
          call stat_update_var_pt( isigma_Nr_1, level, sigma_Nr_1, zt )
       endif

       ! Standard deviation of in-precip rain drop concentration
       ! in PDF component 2.
       if ( isigma_Nr_2 > 0 ) then
          call stat_update_var_pt( isigma_Nr_2, level, sigma_Nr_2, zt )
       endif

       ! Standard deviation of cloud droplet concentration in PDF component 1.
       if ( isigma_Nc_1 > 0 ) then
          call stat_update_var_pt( isigma_Nc_1, level, sigma_Nc_1, zt )
       endif

       ! Standard deviation of cloud droplet concentration in PDF component 2.
       if ( isigma_Nc_2 > 0 ) then
          call stat_update_var_pt( isigma_Nc_2, level, sigma_Nc_2, zt )
       endif

       ! Standard deviation (in-precip) of ln r_r in PDF component 1.
       if ( isigma_rr_1_n > 0 ) then
          call stat_update_var_pt( isigma_rr_1_n, level, sigma_rr_1_n, zt )
       endif

       ! Standard deviation (in-precip) of ln r_r in PDF component 2.
       if ( isigma_rr_2_n > 0 ) then
          call stat_update_var_pt( isigma_rr_2_n, level, sigma_rr_2_n, zt )
       endif

       ! Standard deviation (in-precip) of ln N_r in PDF component 1.
       if ( isigma_Nr_1_n > 0 ) then
          call stat_update_var_pt( isigma_Nr_1_n, level, sigma_Nr_1_n, zt )
       endif

       ! Standard deviation (in-precip) of ln N_r in PDF component 2.
       if ( isigma_Nr_2_n > 0 ) then
          call stat_update_var_pt( isigma_Nr_2_n, level, sigma_Nr_2_n, zt )
       endif

       ! Standard deviation of ln N_c in PDF component 1.
       if ( isigma_Nc_1_n > 0 ) then
          call stat_update_var_pt( isigma_Nc_1_n, level, sigma_Nc_1_n, zt )
       endif

       ! Standard deviation of ln N_c in PDF component 2.
       if ( isigma_Nc_2_n > 0 ) then
          call stat_update_var_pt( isigma_Nc_2_n, level, sigma_Nc_2_n, zt )
       endif

       ! Correlation (in-precip) between s and r_r in PDF component 1.
       if ( icorr_srr_1 > 0 ) then
          call stat_update_var_pt( icorr_srr_1, level, corr_srr_1, zt )
       endif

       ! Correlation (in-precip) between s and r_r in PDF component 2.
       if ( icorr_srr_2 > 0 ) then
          call stat_update_var_pt( icorr_srr_2, level, corr_srr_2, zt )
       endif

       ! Correlation (in-precip) between s and N_r in PDF component 1.
       if ( icorr_sNr_1 > 0 ) then
          call stat_update_var_pt( icorr_sNr_1, level, corr_sNr_1, zt )
       endif

       ! Correlation (in-precip) between s and N_r in PDF component 2.
       if ( icorr_sNr_2 > 0 ) then
          call stat_update_var_pt( icorr_sNr_2, level, corr_sNr_2, zt )
       endif

       ! Correlation between s and N_c in PDF component 1.
       if ( icorr_sNc_1 > 0 ) then
          call stat_update_var_pt( icorr_sNc_1, level, corr_sNc_1, zt )
       endif

       ! Correlation between s and N_c in PDF component 2.
       if ( icorr_sNc_2 > 0 ) then
          call stat_update_var_pt( icorr_sNc_2, level, corr_sNc_2, zt )
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 1.
       if ( icorr_rrNr_1 > 0 ) then
          call stat_update_var_pt( icorr_rrNr_1, level, corr_rrNr_1, zt )
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 2.
       if ( icorr_rrNr_2 > 0 ) then
          call stat_update_var_pt( icorr_rrNr_2, level, corr_rrNr_2, zt )
       endif

       ! Correlation (in-precip) between s and ln r_r in PDF component 1.
       if ( icorr_srr_1_n > 0 ) then
          call stat_update_var_pt( icorr_srr_1_n, level, corr_srr_1_n, zt )
       endif

       ! Correlation (in-precip) between s and ln r_r in PDF component 2.
       if ( icorr_srr_2_n > 0 ) then
          call stat_update_var_pt( icorr_srr_2_n, level, corr_srr_2_n, zt )
       endif

       ! Correlation (in-precip) between s and ln N_r in PDF component 1.
       if ( icorr_sNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_sNr_1_n, level, corr_sNr_1_n, zt )
       endif

       ! Correlation (in-precip) between s and ln N_r in PDF component 2.
       if ( icorr_sNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_sNr_2_n, level, corr_sNr_2_n, zt )
       endif

       ! Correlation between s and ln N_c in PDF component 1.
       if ( icorr_sNc_1_n > 0 ) then
          call stat_update_var_pt( icorr_sNc_1_n, level, corr_sNc_1_n, zt )
       endif

       ! Correlation between s and ln N_c in PDF component 2.
       if ( icorr_sNc_2_n > 0 ) then
          call stat_update_var_pt( icorr_sNc_2_n, level, corr_sNc_2_n, zt )
       endif

       ! Correlation (in-precip) between ln r_r and ln N_r in PDF component 1.
       if ( icorr_rrNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_rrNr_1_n, level, corr_rrNr_1_n, zt )
       endif

       ! Correlation (in-precip) between ln r_r and ln N_r in PDF component 2.
       if ( icorr_rrNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_rrNr_2_n, level, corr_rrNr_2_n, zt )
       endif

       ! Overall variance of r_r.
       if ( irrp2_zt > 0 ) then

          ! Calculate the overall variance of r_r, <r_r'^2>.
          if ( rrainm > rr_tol ) then
             rrp2 &
             = calc_xp2( mu_rr_1_n, mu_rr_2_n, sigma_rr_1_n, sigma_rr_2_n, &
                         mixt_frac, precip_frac_1, precip_frac_2, rrainm )
          else ! r_r = 0.
             rrp2 = zero
          endif

          call stat_update_var_pt( irrp2_zt, level, rrp2, zt )

       endif ! irrp2_zt > 0

       ! Overall variance of N_r.
       if ( iNrp2_zt > 0 ) then

          ! Calculate the overall variance of N_r, <N_r'^2>.
          if ( Nrm > Nr_tol ) then
             Nrp2 &
             = calc_xp2( mu_Nr_1_n, mu_Nr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                         mixt_frac, precip_frac_1, precip_frac_2, Nrm )
          else ! N_r = 0.
             Nrp2 = zero
          endif

          call stat_update_var_pt( iNrp2_zt, level, Nrp2, zt )

       endif ! iNrp2_zt > 0

       ! Variance of KK rain drop mean volume radius.
       if ( iKK_mvr_variance_zt > 0 ) then

          ! Calculate the variance of KK rain drop mean volume radius,
          ! < R_vr'^2 >.
          if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then
             KK_mvr_variance &
             = variance_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1, mu_Nr_2, &
                                mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1_n, &
                                sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                                corr_rrNr_1_n, corr_rrNr_2_n, &
                                KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                                precip_frac_1, precip_frac_2 )
          else  ! r_r or N_r = 0.
             KK_mvr_variance = zero
          endif

          call stat_update_var_pt( iKK_mvr_variance_zt, level, &
                                   KK_mvr_variance, zt )

       endif ! iKK_mvr_variance_zt > 0

    endif ! l_stats_samp


    return

  end subroutine KK_upscaled_stats

  !=============================================================================
  subroutine KK_stats_output( KK_evap_tndcy, KK_auto_tndcy, &
                              KK_accr_tndcy, KK_mean_vol_rad, &
                              KK_Nrm_evap_tndcy, KK_Nrm_auto_tndcy, &
                              l_stats_samp, level )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use stats_variables, only: &
        zt,              & ! Variable(s)
        im_vol_rad_rain, &
        irrainm_cond,    &
        irrainm_auto,    &
        irrainm_accr,    &
        iNrm_cond,       &
        iNrm_auto

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy,     & ! Mean KK (dr_r/dt) due to evaporation    [(kg/kg)/s]
      KK_auto_tndcy,     & ! Mean KK (dr_r/dt) due to autoconversion [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK (dr_r/dt) due to accretion      [(kg/kg)/s]
      KK_mean_vol_rad,   & ! Mean KK rain drop mean volume radius            [m]
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation   [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.     [(num/kg)/s]

    logical, intent(in) :: &
      l_stats_samp         ! Flag to sample statistical output

    integer, intent(in) :: & 
      level    ! Vertical level index


    ! Statistics
    if ( l_stats_samp ) then

       ! Mean rain water mixing ratio microphysics tendencies.
       if ( irrainm_cond > 0 ) then
          call stat_update_var_pt( irrainm_cond, level, KK_evap_tndcy, zt )
       endif

       if ( irrainm_auto > 0 ) then
          call stat_update_var_pt( irrainm_auto, level, KK_auto_tndcy, zt )
       endif

       if ( irrainm_accr > 0 ) then
          call stat_update_var_pt( irrainm_accr, level, KK_accr_tndcy, zt )
       endif

       ! Rain drop mean volume radius.
       if ( im_vol_rad_rain > 0 ) then
          call stat_update_var_pt( im_vol_rad_rain, level, KK_mean_vol_rad, zt )
       endif

       ! Mean rain drop concentration microphysics tendencies.
       if ( iNrm_cond > 0 ) then
          call stat_update_var_pt( iNrm_cond, level, KK_Nrm_evap_tndcy, zt )
       endif

       if ( iNrm_auto > 0 ) then
          call stat_update_var_pt( iNrm_auto, level, KK_Nrm_auto_tndcy, zt )
       endif

    endif ! l_stats_samp


    return

  end subroutine KK_stats_output

  !=============================================================================
  subroutine KK_sedimentation( nz, cloud_top_level, KK_mean_vol_rad, Vrr, VNr )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        micron_per_m    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use constants_clubb, only: &
        zero   ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,              & ! Number of model vertical grid levels
      cloud_top_level    ! Vertical level index of cloud top

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      KK_mean_vol_rad    ! KK rain drop mean volume radius       [m]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(:), pointer, intent(inout) ::  &
      Vrr, & ! Mean sedimentation velocity of < r_r >            [m/s]
      VNr    ! Mean sedimentation velocity of < N_r >            [m/s]

    ! Local Variables
    integer :: k ! Loop iterator


    !!! Mean sedimentation velocities
    do k = 1, nz-1, 1

       ! Mean sedimentation velocity of rain water mixing ratio.
       Vrr(k) = - ( 0.012_core_rknd * ( micron_per_m * KK_mean_vol_rad(k) ) &
                    - 0.2_core_rknd )

       ! Mean sedimentation velocity of rain water mixing ratio cannot have a
       ! positive value.
       if ( Vrr(k) > zero ) then
          Vrr(k) = zero
       endif

       ! Mean sedimentation velocity of rain drop concentration.
       VNr(k) = - ( 0.007_core_rknd * ( micron_per_m * KK_mean_vol_rad(k) )  &
                    - 0.1_core_rknd )

       ! Mean sedimentation velocity of rain drop concentration cannot have a
       ! positive value.
       if ( VNr(k) > zero ) then
          VNr(k) = zero
       endif

    enddo ! Sedimentation velocity loop: k = 1, nz-1, 1

    !!! Mean sedimentation above cloud top should have a value of 0.
    if ( cloud_top_level > 1 ) then
       Vrr(cloud_top_level+1:nz-1) = zero
       VNr(cloud_top_level+1:nz-1) = zero
    endif

    !!! Boundary conditions for sedimentation velocities.

    ! The flux of rain water through the model top is 0.
    ! Vrr and VNr are set to 0 at the highest model level.
    Vrr(nz) = zero
    VNr(nz) = zero


    return

  end subroutine KK_sedimentation

  !=============================================================================
  function variance_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1, mu_Nr_2, &
                            mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1_n, &
                            sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, &
                            corr_rrNr_1_n, corr_rrNr_2_n, &
                            KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                            precip_frac_1, precip_frac_2 )

    ! Description:
    ! This function calculates the variance of KK mean volume radius of rain
    ! drops (R_vr), which can be written as < R_vr'^2 >.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, & ! Constant(s)
        two

    use KK_upscaled_means, only:  &
        bivar_LL_mean_eq  ! Procedure(s)

    use parameters_microphys, only: &
        KK_mvr_rr_exp, & ! Variable(s)
        KK_mvr_Nr_exp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) in-precip (ip)  [-]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip              [-]
      mu_Nr_1,         & ! Mean of Nr (1st PDF component) ip                 [-]
      mu_Nr_2,         & ! Mean of Nr (2nd PDF component) ip                 [-]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip              [-]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip              [-]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_rrNr_1_n,   & ! Corr. betw. ln rr & ln Nr (1st PDF comp.) ip      [-]
      corr_rrNr_2_n,   & ! Corr. betw. ln rr & ln Nr (2nd PDF comp.) ip      [-]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      variance_KK_mvr    ! Variance of KK rain drop mean volume radius     [m^2]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha_exp, & ! Exponent on r_r
      beta_exp     ! Exponent on N_r


    ! Values of the KK exponents.
    alpha_exp = KK_mvr_rr_exp
    beta_exp  = KK_mvr_Nr_exp

    ! Calculate the covariance of r_r and KK mean volume radius of rain drops.
    variance_KK_mvr  &
    = KK_mvr_coef**2 &
      * ( mixt_frac &
          * precip_frac_1 &
          * bivar_LL_mean_eq( mu_rr_1_n, mu_Nr_1, mu_Nr_1_n, &
                              sigma_rr_1_n, sigma_Nr_1_n, corr_rrNr_1_n, &
                              two * alpha_exp, two * beta_exp ) &
        + ( one - mixt_frac ) &
          * precip_frac_2 &
          * bivar_LL_mean_eq( mu_rr_2_n, mu_Nr_2, mu_Nr_2_n, &
                              sigma_rr_2_n, sigma_Nr_2_n, corr_rrNr_2_n, &
                              two * alpha_exp, two * beta_exp ) &
        ) &
      - KK_mean_vol_rad**2


    return

  end function variance_KK_mvr

!===============================================================================

end module KK_microphys_module

! $Id$
!===============================================================================
module KK_microphys_module

  implicit none

  private

  public :: KK_micro_driver, &
            precip_fraction, &
            KK_in_precip_values, &
            KK_upscaled_setup, &
            KK_stat_output

  private :: KK_upscaled_means_driver, &
             KK_sed_vel_covars, &
             KK_upscaled_covar_driver

  contains

  !=============================================================================
  subroutine KK_micro_driver( dt, nz, l_stats_samp, l_local_kk, &
                              l_latin_hypercube, thlm, wm_zt, p_in_Pa, &
                              exner, rho, cloud_frac, pdf_params, w_std_dev, &
                              dzq, rcm, Ncm, s_mellor, rvm, Nc0_in_cloud, &
                              hydromet, &
                              hydromet_mc, hydromet_vel, &
                              rcm_mc, rvm_mc, thlm_mc, &
                              hydromet_vel_covar, hydromet_vel_covar_zt,  &
                              wprtp_mc_tndcy, wpthlp_mc_tndcy, &
                              rtp2_mc_tndcy, thlp2_mc_tndcy, rtpthlp_mc_tndcy, &
                              KK_auto_tndcy, KK_accr_tndcy )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm, &  ! Procedure(s)
        zm2zt, &
        gr

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
        zero,         &
        rho_lw,       &
        rc_tol,       &
        rr_tol,       &
        Nr_tol,       &
        Nc_tol,       &
        cm3_per_m3,   &
        micron_per_m, &
        eps

    use parameters_microphys, only: &
        KK_auto_Nc_exp,      & ! Constant(s)
        C_evap,              &
        l_var_covar_src,     & ! Flag for using variance/covariance src terms
        l_const_Nc_in_cloud    ! Flag to use a constant value of N_c within cloud

    use KK_utilities, only: &
        G_T_p  ! Procedure(s)

    use KK_local_means, only: &
        KK_mvr_local_mean,  & ! Procedure(s)
        KK_evap_local_mean, &
        KK_auto_local_mean, &
        KK_accr_local_mean

    use KK_Nrm_tendencies, only: &
        KK_Nrm_evap, & ! Procedure(s)
        KK_Nrm_auto

    use saturation, only: & 
        sat_mixrat_liq  ! Procedure(s)

    use array_index, only: &
        iirrainm, & ! Constant(s)
        iiNrm

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type, only: & 
        stat_update_var_pt, & ! Procedure(s)
        stat_update_var

    use stats_variables, only: & 
        zt,               & ! Variable(s)
        im_vol_rad_rain,  &
        irrainm_cond,     &
        irrainm_auto,     &
        irrainm_accr,     &
        irrainm_src_adj,  &
        irrainm_cond_adj, &
        iNrm_cond,        &
        iNrm_auto,        &
        iNrm_src_adj,     &
        iNrm_cond_adj,    &
        iprecip_frac

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
      l_stats_samp,      & ! Flag to sample statistics
      l_local_kk,        & ! Flag to use the local form of KK microphysics
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

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(out) :: &
      hydromet_vel_covar,    & ! Covariance of V_xx & x_x (m-levs)  [units(m/s)]
      hydromet_vel_covar_zt    ! Covariance of V_xx & x_x (t-levs)  [units(m/s)]

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
      T_liq_in_K, & ! Mean liquid water temperature, T_l           [K]
      r_sl,       & ! Liquid water sat. mixing ratio, r_s(T_l,p)   [kg/kg]
      Beta_Tl       ! Parameter Beta, Beta(T_l)                    [1/(kg/kg)]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    real( kind = core_rknd ), dimension(nz) :: &
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component) [-]
      precip_frac      ! Precipitation fraction (overall)           [-]

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
      corr_wNc,      & ! Correlation between Nc & w (both components)        [-]
      mixt_frac        ! Mixture fraction                                    [-]

    real( kind = core_rknd ), dimension(:), pointer :: &
      Vrrprrp, & ! Covariance of V_rr and r_r (momentum levels)  [(m/s)(kg/kg)]
      VNrpNrp    ! Covariance of V_Nr and N_r (momentum levels)  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(:), pointer :: &
      Vrrprrp_zt, & ! Covariance of V_rr and r_r; thermo. levs.  [(m/s)(kg/kg)]
      VNrpNrp_zt    ! Covariance of V_Nr and N_r; thermo. levs.  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz) :: &
      wprtp_mc_tndcy_zt,   & ! Micro. tend. for <w'rt'>; t-lev   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy_zt,  & ! Micro. tend. for <w'thl'>; t-lev  [m*K/s^2]
      rtp2_mc_tndcy_zt,    & ! Micro. tend. for <rt'^2>; t-lev   [(kg/kg)^2/s]
      thlp2_mc_tndcy_zt,   & ! Micro. tend. for <thl'^2>; t-lev  [K^2/s]
      rtpthlp_mc_tndcy_zt    ! Micro. tend. for <rt'thl'>; t-lev [K*(kg/kg)/s]

    real( kind = core_rknd ) ::  &
      rrainm_source,     & ! Total source term rate for rrainm       [(kg/kg)/s]
      Nrm_source,        & ! Total source term rate for Nrm         [(num/kg)/s]
      rrainm_src_max,    & ! Maximum allowable rrainm source rate    [(kg/kg)/s]
      rrainm_auto_ratio, & ! Ratio of rrainm autoconv to overall source term [-]
      total_rc_needed      ! Amount of r_c needed to over the timestep
                           ! for rain source terms                       [kg/kg]

    real( kind = core_rknd ), dimension(nz) ::  &
      rrainm_src_adj,  & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj,     & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrainm_evap_net, & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net       ! Net evaporation rate of <N_r>            [(num/kg)/s]

    ! changes by janhft 10/04/12
    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zm,  & ! Covariance of s and w on the zm-grid    [(m/s)(kg/kg)]
      wprrp_zm, & ! Covariance of r_r and w on the zm-grid  [(m/s)(kg/kg)]
      wpNrp_zm, & ! Covariance of N_r and w on the zm-grid  [(m/s)(#/kg)]
      wpNcp_zm, & ! Covariance of N_c and w on the zm-grid  [(m/s)(#/kg)]
      wpsp_zt,  & ! Covariance of s and w on the zt-grid    [(m/s)(kg/kg)]
      wprrp_zt, & ! Covariance of r_r and w on the zt-grid  [(m/s)(kg/kg)]
      wpNrp_zt, & ! Covariance of N_r and w on the zt-grid  [(m/s)(#/kg)]
      wpNcp_zt    ! Covariance of N_c and w on the zt-grid  [(m/s)(#/kg)]
    ! end changes by janhft 10/04/12

    logical :: &
      l_upscaled,        & ! Flag for using upscaled KK microphysics.
      l_src_adj_enabled, & ! Flag to enable rrainm/Nrm source adjustment
      l_stats_samp_in_sub  ! Used to disable stats when SILHS is enabled

    integer :: &
      k   ! Loop index

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

    call KK_init_micro_driver( nz, hydromet, hydromet_mc, hydromet_vel, & ! Intent(in)
                                   hydromet_vel_covar, hydromet_vel_covar_zt, &
                                   rrainm, Nrm, Vrr, VNr, & ! Intent(out)
                                   rrainm_mc_tndcy, Nrm_mc_tndcy, &
                                   KK_auto_tndcy, KK_accr_tndcy, &
                                   Vrrprrp, VNrpNrp, Vrrprrp_zt, &
                                   VNrpNrp_zt, l_src_adj_enabled )

    if ( .not. l_local_kk ) then
       l_upscaled = .true.
    else
       l_upscaled = .false.
    endif

    ! Precipitation fraction
    if ( l_use_precip_frac ) then

       call precip_fraction( nz, rrainm, cloud_frac, precip_frac )

    else

       precip_frac = one

    endif

    ! Set precipitation fraction for each PDF component.
    precip_frac_1 = precip_frac
    precip_frac_2 = precip_frac

    ! Statistics
    if ( l_stats_samp_in_sub ) then
       call stat_update_var( iprecip_frac, precip_frac, zt )
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

       wprrp_zm(1:nz-1) &
       = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), &
                   rrainm(1:nz-1) / max( precip_frac(1:nz-1), eps ), &
                   rrainm(2:nz) / max( precip_frac(2:nz), eps ), &
                   gr%invrs_dzm(1:nz-1) )

       wpNrp_zm(1:nz-1) &
       = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), &
                   Nrm(1:nz-1) / max( precip_frac(1:nz-1), eps ), & 
                   Nrm(2:nz) / max( precip_frac(2:nz), eps ), &
                   gr%invrs_dzm(1:nz-1) )

       wpNcp_zm(1:nz-1) = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), Ncm(1:nz-1), & 
                                    Ncm(2:nz), gr%invrs_dzm(1:nz-1) )

       ! Boundary conditions; We are assuming constant flux at the top.
       wprrp_zm(nz) = wprrp_zm(nz-1)
       wpNrp_zm(nz) = wpNrp_zm(nz-1)
       wpNcp_zm(nz) = wpNcp_zm(nz-1)

       ! interpolate back to zt-grid
       wpsp_zt  = zm2zt(wpsp_zm)
       wprrp_zt = zm2zt(wprrp_zm)
       wpNrp_zt = zm2zt(wpNrp_zm)
       wpNcp_zt = zm2zt(wpNcp_zm)

    endif


    ! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

      ! Compute supersaturation via s1, s2.
      !     Larson et al 2002, JAS, Vol 59, p 3534.
      ! This allows a more direct comparison of local, upscaled formulas.

      call KK_supersaturation( thlm(k), exner(k), p_in_Pa(k), rho(k), & ! Intent(in)
                               KK_evap_coef, KK_auto_coef, & ! Intent(out)
                               KK_accr_coef, KK_mvr_coef )

       !!! KK rain water mixing ratio microphysics tendencies.
       if ( l_upscaled ) then

          call KK_in_precip_values( rrainm(k), Nrm(k), rcm(k), &
                                    precip_frac_1(k), precip_frac_2(k), &
                                    mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                    sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                    sigma_Nr_2, corr_srr_1, corr_srr_2, &
                                    corr_sNr_1, corr_sNr_2, corr_rrNr_1, &
                                    corr_rrNr_2 )

          call KK_upscaled_setup( rcm(k), rrainm(k), Nrm(k), Ncm(k), &
                                  mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                  sigma_rr_1, sigma_rr_2, &
                                  sigma_Nr_1, sigma_Nr_2, &
                                  wpsp_zt(k), wprrp_zt(k), wpNrp_zt(k), &
                                  wpNcp_zt(k), w_std_dev(k), pdf_params(k), &
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
                                  corr_sw, corr_wrr, corr_wNr, corr_wNc, &
                                  mixt_frac )

          call KK_stat_output( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
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
                               corr_rrNr_2_n, corr_sw, corr_wrr, corr_wNr, &
                               corr_wNc, k, l_stats_samp_in_sub )

          !!! Calculate the values of the upscaled KK microphysics tendencies.
          call KK_upscaled_means_driver( rrainm(k), Nrm(k), Ncm(k), &
                                         mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                         mu_Nr_1_n, mu_Nr_2_n, mu_Nc_1_n, &
                                         mu_Nc_2_n, sigma_s_1, sigma_s_2, &
                                         sigma_rr_1_n, sigma_rr_2_n, &
                                         sigma_Nr_1_n, sigma_Nr_2_n, &
                                         sigma_Nc_1_n, sigma_Nc_2_n, &
                                         corr_srr_1_n, corr_srr_2_n, &
                                         corr_sNr_1_n, corr_sNr_2_n, &
                                         corr_sNc_1_n, corr_sNc_2_n, &
                                         corr_rrNr_1_n, corr_rrNr_2_n, &
                                         mixt_frac, precip_frac_1(k), &
                                         precip_frac_2(k), Nc0_in_cloud(k), &
                                         l_const_Nc_in_cloud, &
                                         KK_evap_coef, KK_auto_coef, &
                                         KK_accr_coef, KK_mvr_coef, &
                                         KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                         KK_accr_tndcy(k), KK_mean_vol_rad(k) )

          call KK_sed_vel_covars( rrainm(k), Nrm(k), KK_mean_vol_rad(k), &
                                  mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                  sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                  sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                                  KK_mvr_coef, mixt_frac, precip_frac_1(k), &
                                  precip_frac_2(k), k, l_stats_samp_in_sub, &
                                  Vrrprrp_zt(k), VNrpNrp_zt(k) )
          

          if ( l_var_covar_src ) then

            call KK_upscaled_covar_driver( wm_zt(k), exner(k), rcm(k),  &
                                           rrainm(k), Nrm(k), Ncm(k), &
                                           mu_s_1, mu_s_2, mu_rr_1_n, mu_Nr_1_n, &
                                           mu_Nc_1_n, sigma_s_1, sigma_s_2, &
                                           sigma_rr_1_n, sigma_Nr_1_n, sigma_Nc_1_n, &
                                           corr_srr_1_n, corr_srr_2_n, &
                                           corr_sNr_1_n, corr_sNr_2_n, &
                                           corr_sNc_1_n, corr_sNc_2_n, &
                                           corr_rrNr_1_n,  mixt_frac, &
                                           precip_frac(k), Nc0_in_cloud(k), &
                                           l_const_Nc_in_cloud, &
                                           KK_evap_coef, KK_auto_coef, &
                                           KK_accr_coef, KK_evap_tndcy(k), &
                                           KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                           pdf_params(k), k, &
                                           l_stats_samp_in_sub, &
                                           wprtp_mc_tndcy_zt(k), &
                                           wpthlp_mc_tndcy_zt(k), &
                                           rtp2_mc_tndcy_zt(k), &
                                           thlp2_mc_tndcy_zt(k), &
                                           rtpthlp_mc_tndcy_zt(k) )

          endif


       else  ! local KK

          !!! Calculate the local KK rain drop mean volume radius.
          if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

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

          ! Set the covariances of hydrometeor sedimentation velocities and
          ! their associated hydrometeors (<V_rr'r_r'> and <V_Nr'N_r'>) to 0.
          Vrrprrp_zt(k) = zero
          VNrpNrp_zt(k) = zero

       endif ! l_upscaled


      call KK_microphys_adjust( rrainm(k), Nrm(k), rcm(k), exner(k), & ! Intent(in)
                                 KK_evap_tndcy(k), KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                 dt, l_src_adj_enabled, &
                                 rrainm_source, Nrm_source, & ! Intent(out)
                                 KK_Nrm_auto_tndcy(k), KK_Nrm_evap_tndcy(k), &
                                 rrainm_src_adj(k), Nrm_src_adj(k), &
                                 rrainm_evap_net(k), Nrm_evap_net(k), &
                                 rrainm_mc_tndcy(k), Nrm_mc_tndcy(k), &
                                 rvm_mc(k), rcm_mc(k), thlm_mc(k) )


      if ( l_stats_samp_in_sub ) then

        call KK_stats_output_samp_in_sub( KK_mean_vol_rad(k), KK_evap_tndcy(k), & ! Intent(in)
                                        KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                        KK_Nrm_evap_tndcy(k), KK_Nrm_auto_tndcy(k), &
                                        rrainm_src_adj(k), Nrm_src_adj(k), &
                                        rrainm_evap_net(k), Nrm_evap_net(k), &
                                        k )

      endif

    enddo  ! Microphysics tendency loop: k = 2, nz, 1


    if ( l_upscaled .and. l_var_covar_src ) then

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

       ! Set values to 0.
       wprtp_mc_tndcy   = zero
       wpthlp_mc_tndcy  = zero
       rtp2_mc_tndcy    = zero
       thlp2_mc_tndcy   = zero
       rtpthlp_mc_tndcy = zero

    endif


    call KK_sedimentation( nz, & ! Intent(in)
                           Vrr, VNr, & ! Intent(InOut)
                           rrainm_mc_tndcy, Nrm_mc_tndcy, &
                           rcm_mc, rvm_mc, thlm_mc, &
                           Vrrprrp, VNrpNrp, &
                           Vrrprrp_zt, VNrpNrp_zt, &
                           KK_mean_vol_rad )
    return

  end subroutine KK_micro_driver

  !=============================================================================
  subroutine precip_fraction( nz, rrainm, cloud_frac, precip_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero,           & ! Constant(s)
        rr_tol,         &
        cloud_frac_min

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rrainm,     & ! Mean rain water mixing ratio     [kg/kg]
      cloud_frac    ! Cloud fraction                   [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac    ! Precipitation fraction          [-]

    ! Local Variables
    real( kind = core_rknd ), parameter :: &
      precip_frac_tol = cloud_frac_min  ! Minimum precip. frac. [-]
    
    integer :: &
      k   ! Loop index


    precip_frac(nz) = zero

    do k = nz-1, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )

       if ( rrainm(k) > rr_tol .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find rain at this grid level, but no cloud at
          ! or above this grid level, set precipitation fraction to a minimum
          ! threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( rrainm(k) < rr_tol .and. precip_frac(k) < precip_frac_tol ) then

          ! Mean rain water mixing ratio is less than the tolerance amount.  It
          ! is considered to have a value of 0.  There is not any rain at this
          ! grid level.  There is also no cloud at or above this grid level, so
          ! set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine KK_in_precip_values( rrainm, Nrm, rcm, &
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
      rrainm,        & ! Mean rain water mixing ratio, < r_r >       [kg/kg]
      Nrm,           & ! Mean rain drop concentration, < N_r >       [num/kg]
      rcm,           & ! Mean cloud water mixing ratio, < r_c >      [kg/kg]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)  [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)  [-]

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


    ! Mean of in-precip rain water mixing ratio in PDF component 1.
    if ( rrainm > rr_tol ) then
       mu_rr_1 = rrainm / precip_frac_1
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.
       mu_rr_1 = zero
    endif

    ! Mean of in-precip rain water mixing ratio in PDF component 2.
    if ( rrainm > rr_tol ) then
       mu_rr_2 = rrainm / precip_frac_2
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.
       mu_rr_2 = zero
    endif

    ! Mean of in-precip rain drop concentration in PDF component 1.
    if ( Nrm > Nr_tol ) then
       mu_Nr_1 = Nrm / precip_frac_1
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.
       mu_Nr_1 = zero
    endif

    ! Mean of in-precip rain drop concentration in PDF component 2.
    if ( Nrm > Nr_tol ) then
       mu_Nr_2 = Nrm / precip_frac_2
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.
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
    !
    ! Set the value of the parameters based on whether the altitude is above or
    ! below cloud base.  Determine whether there is cloud at any given vertical
    ! level.  In order for a vertical level to have cloud, the amount of cloud
    ! water (rcm) must be greater than or equal to the tolerance level (rc_tol).
    ! If there is cloud at a given vertical level, then the ###_cloud value is
    ! used.  Otherwise, the ###_below value is used.

    ! Standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( rrainm > rr_tol ) then
       if ( rcm > rc_tol ) then
          sigma_rr_1 = sqrt( rrp2_on_rrm2_cloud ) * mu_rr_1
       else
          sigma_rr_1 = sqrt( rrp2_on_rrm2_below ) * mu_rr_1
       endif
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain
       ! water mixing ratio does not vary at this grid level.
       sigma_rr_1 = zero
    endif

    ! Standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( rrainm > rr_tol ) then
       if ( rcm > rc_tol ) then
          sigma_rr_2 = sqrt( rrp2_on_rrm2_cloud ) * mu_rr_2
       else
          sigma_rr_2 = sqrt( rrp2_on_rrm2_below ) * mu_rr_2
       endif
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain
       ! water mixing ratio does not vary at this grid level.
       sigma_rr_2 = zero
    endif

    ! Standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( Nrm > Nr_tol ) then
       if ( rcm > rc_tol ) then
          sigma_Nr_1 = sqrt( Nrp2_on_Nrm2_cloud ) * mu_Nr_1
       else
          sigma_Nr_1 = sqrt( Nrp2_on_Nrm2_below ) * mu_Nr_1
       endif
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain drop
       ! concentration does not vary at this grid level.
       sigma_Nr_1 = zero
    endif

    ! Standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( Nrm > Nr_tol ) then
       if ( rcm > rc_tol ) then
          sigma_Nr_2 = sqrt( Nrp2_on_Nrm2_cloud ) * mu_Nr_2
       else
          sigma_Nr_2 = sqrt( Nrp2_on_Nrm2_below ) * mu_Nr_2
       endif
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain drop
       ! concentration does not vary at this grid level.
       sigma_Nr_2 = zero
    endif

    ! Correlation (in-precip) between s and r_r.
    if ( rrainm > rr_tol ) then

       ! Correlation (in-precip) between s and r_r in PDF component 1.
       if ( rcm > rc_tol ) then
          corr_srr_1 = corr_srr_NL_cloud
       else
          corr_srr_1 = corr_srr_NL_below
       endif

       ! Correlation (in-precip) between s and r_r in PDF component 2.
       if ( rcm > rc_tol ) then
          corr_srr_2 = corr_srr_NL_cloud
       else
          corr_srr_2 = corr_srr_NL_below
       endif

    else

       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The correlations involving rain water mixing ratio
       ! are 0 since rain water mixing ratio does not vary at this grid level.
       corr_srr_1 = zero
       corr_srr_2 = zero

    endif

    ! Correlation (in-precip) between s and N_r.
    if ( Nrm > Nr_tol ) then

       ! Correlation (in-precip) between s and N_r in PDF component 1.
       if ( rcm > rc_tol ) then
          corr_sNr_1 = corr_sNr_NL_cloud
       else
          corr_sNr_1 = corr_sNr_NL_below
       endif

       ! Correlation (in-precip) between s and N_r in PDF component 2.
       if ( rcm > rc_tol ) then
          corr_sNr_2 = corr_sNr_NL_cloud
       else
          corr_sNr_2 = corr_sNr_NL_below
       endif

    else

       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The correlations involving rain drop concentration
       ! are 0 since rain water mixing ratio does not vary at this grid level.
       corr_sNr_1 = zero
       corr_sNr_2 = zero

    endif

    ! Correlation (in-precip) between r_r and N_r.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       ! Correlation (in-precip) between r_r and N_r in PDF component 1.
       if ( rcm > rc_tol ) then
          corr_rrNr_1 = corr_rrNr_LL_cloud
       else
          corr_rrNr_1 = corr_rrNr_LL_below
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 2.
       if ( rcm > rc_tol ) then
          corr_rrNr_2 = corr_rrNr_LL_cloud
       else
          corr_rrNr_2 = corr_rrNr_LL_below
       endif

    else

       ! Mean in-precip rain water mixing ratio and (or) mean in-precip rain
       ! drop concentration are (is) less than their (its) respective tolerance
       ! amount(s), and are (is) considered to have a value of 0.  There is not
       ! any rain at this grid level.  The correlation is 0 since rain does not
       ! vary at this grid level.
       corr_rrNr_1 = zero
       corr_rrNr_2 = zero

    endif


    return    

  end subroutine KK_in_precip_values

  !=============================================================================
  subroutine KK_upscaled_setup( rcm, rrainm, Nrm, Ncm, &
                                mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                sigma_rr_1, sigma_rr_2, &
                                sigma_Nr_1, sigma_Nr_2, &
                                wpsp, wprrp, wpNrp, &
                                wpNcp, stdev_w, pdf_params, &
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
                                corr_sw, corr_wrr, corr_wNr, corr_wNc, &
                                mixt_frac )

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
      mu_rr_1,    & ! Mean of rr (1st PDF component) in-precip (ip)      [kg/kg]
      mu_rr_2,    & ! Mean of rr (2nd PDF component) ip                  [kg/kg]
      mu_Nr_1,    & ! Mean of Nr (1st PDF component) ip                 [num/kg]
      mu_Nr_2,    & ! Mean of Nr (2nd PDF component) ip                 [num/kg]
      sigma_rr_1, & ! Standard deviation of rr (1st PDF component) ip    [kg/kg]
      sigma_rr_2, & ! Standard deviation of rr (2nd PDF component) ip    [kg/kg]
      sigma_Nr_1, & ! Standard deviation of Nr (1st PDF component) ip   [num/kg]
      sigma_Nr_2, & ! Standard deviation of Nr (2nd PDF component) ip   [num/kg]
      wpsp,       & ! Covariance of w and s                         [(m/s)kg/kg]
      wprrp,      & ! Covariance of w and rrain                     [(m/s)kg/kg]
      wpNrp,      & ! Covariance of w and Nr                       [(m/s)num/kg]
      wpNcp,      & ! Covariance of w and Nc                       [(m/s)num/kg]
      stdev_w       ! Standard deviation of w                              [m/s]

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
      corr_wNc,      & ! Correlation between Nc & w (both components)        [-]
      mixt_frac        ! Mixture fraction                                    [-]

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
    mixt_frac = pdf_params%mixt_frac

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

    ! Correlation between s and N_c.
    if ( Ncm > Nc_tol ) then

       ! Correlation between s and N_c in PDF component 1.
       if ( rcm > rc_tol ) then
          corr_sNc_1 = corr_sNc_NL_cloud
       else
          corr_sNc_1 = corr_sNc_NL_below
       endif

       ! Correlation between s and N_c in PDF component 2.
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
       corr_sNc_1 = zero
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
          corr_wrr = calc_w_corr( wprrp, stdev_w, sigma_rr_1, w_tol, rr_tol )
          corr_wNr = calc_w_corr( wpNrp, stdev_w, sigma_Nr_1, w_tol, Nr_tol )
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
    if ( mu_rr_1 > rr_tol ) then
       mu_rr_1_n = mean_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The value of mu_rr_n should be -inf.  It will be
       ! set to -huge for purposes of assigning it a value.  This value will not
       ! be used again in the CLUBB code.
       !mu_rr_1_n = -huge( mu_rr_1_n )
       ! Some compilers have issues outputting to stats files (in single
       ! precision) when the default CLUBB kind is in double precision.
       ! Set to -huge for single precision.
       mu_rr_1_n = -huge( 0.0 )
    endif

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 2.
    if ( mu_rr_2 > rr_tol ) then
       mu_rr_2_n = mean_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The value of mu_rr_n should be -inf.  It will be
       ! set to -huge for purposes of assigning it a value.  This value will not
       ! be used again in the CLUBB code.
       !mu_rr_2_n = -huge( mu_rr_2_n )
       ! Some compilers have issues outputting to stats files (in single
       ! precision) when the default CLUBB kind is in double precision.
       ! Set to -huge for single precision.
       mu_rr_2_n = -huge( 0.0 )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 1.
    if ( mu_Nr_1 > Nr_tol ) then
       mu_Nr_1_n = mean_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The value of mu_Nr_n should be -inf.  It will be
       ! set to -huge for purposes of assigning it a value.  This value will not
       ! be used again in the CLUBB code.
       !mu_Nr_1_n = -huge( mu_Nr_1_n )
       ! Some compilers have issues outputting to stats files (in single
       ! precision) when the default CLUBB kind is in double precision.
       ! Set to -huge for single precision.
       mu_Nr_1_n = -huge( 0.0 )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 2.
    if ( mu_Nr_2 > Nr_tol ) then
       mu_Nr_2_n = mean_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The value of mu_Nr_n should be -inf.  It will be
       ! set to -huge for purposes of assigning it a value.  This value will not
       ! be used again in the CLUBB code.
       !mu_Nr_2_n = -huge( mu_Nr_2_n )
       ! Some compilers have issues outputting to stats files (in single
       ! precision) when the default CLUBB kind is in double precision.
       ! Set to -huge for single precision.
       mu_Nr_2_n = -huge( 0.0 )
    endif

    ! Normalized mean of cloud droplet concentration in PDF component 1.
    if ( Ncm > Nc_tol ) then
       mu_Nc_1_n = mean_L2N( mu_Nc_1, sigma_Nc_1**2 )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There isn't any cloud at this
       ! grid level.  The value of mu_Nc_n should be -inf.  It will be set to
       ! -huge for purposes of assigning it a value.  This value will not be
       ! used again in the CLUBB code.
       !mu_Nc_1_n = -huge( mu_Nc_1_n )
       ! Some compilers have issues outputting to stats files (in single
       ! precision) when the default CLUBB kind is in double precision.
       ! Set to -huge for single precision.
       mu_Nc_1_n = -huge( 0.0 )
    endif

    ! Normalized mean of cloud droplet concentration in PDF component 2.
    if ( Ncm > Nc_tol ) then
       mu_Nc_2_n = mean_L2N( mu_Nc_2, sigma_Nc_2**2 )
    else
       ! Mean cloud droplet concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There isn't any cloud at this
       ! grid level.  The value of mu_Nc_n should be -inf.  It will be set to
       ! -huge for purposes of assigning it a value.  This value will not be
       ! used again in the CLUBB code.
       !mu_Nc_2_n = -huge( mu_Nc_2_n )
       ! Some compilers have issues outputting to stats files (in single
       ! precision) when the default CLUBB kind is in double precision.
       ! Set to -huge for single precision.
       mu_Nc_2_n = -huge( 0.0 )
    endif

    !!! Calculate the normalized standard deviation of variables that have
    !!! an assumed lognormal distribution, given the mean and variance of
    !!! those variables.

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( mu_rr_1 > rr_tol ) then
       sigma_rr_1_n = stdev_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain
       ! water mixing ratio does not vary at this grid level.
       sigma_rr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( mu_rr_2 > rr_tol ) then
       sigma_rr_2_n = stdev_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain
       ! water mixing ratio does not vary at this grid level.
       sigma_rr_2_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( mu_Nr_1 > Nr_tol ) then
       sigma_Nr_1_n = stdev_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain drop
       ! concentration does not vary at this grid level.
       sigma_Nr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( mu_Nr_2 > Nr_tol ) then
       sigma_Nr_2_n = stdev_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The standard deviation is simply 0 since rain drop
       ! concentration does not vary at this grid level.
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
    if ( mu_rr_1 > rr_tol ) then
       corr_srr_1_n = corr_NL2NN( corr_srr_1, sigma_rr_1_n )
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The correlations involving rain water mixing ratio
       ! are 0 since rain water mixing ratio does not vary at this grid level.
       corr_srr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 2.
    if ( mu_rr_2 > rr_tol ) then
       corr_srr_2_n = corr_NL2NN( corr_srr_2, sigma_rr_2_n )
    else
       ! Mean in-precip rain water mixing ratio is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The correlations involving rain water mixing ratio
       ! are 0 since rain water mixing ratio does not vary at this grid level.
       corr_srr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 1.
    if ( mu_Nr_1 > Nr_tol ) then
       corr_sNr_1_n = corr_NL2NN( corr_sNr_1, sigma_Nr_1_n )
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The correlations involving rain drop concentration
       ! are 0 since rain drop concentration does not vary at this grid level.
       corr_sNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 2.
    if ( mu_Nr_2 > Nr_tol ) then
       corr_sNr_2_n = corr_NL2NN( corr_sNr_2, sigma_Nr_2_n )
    else
       ! Mean in-precip rain drop concentration is less than the tolerance
       ! amount.  It is considered to have a value of 0.  There is not any rain
       ! at this grid level.  The correlations involving rain drop concentration
       ! are 0 since rain drop concentration does not vary at this grid level.
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
    if ( mu_rr_1 > rr_tol .and. mu_Nr_1 > Nr_tol ) then
       corr_rrNr_1_n = corr_LL2NN( corr_rrNr_1, sigma_rr_1_n, sigma_Nr_1_n )
    else
       ! Mean in-precip rain water mixing ratio and (or) mean in-precip rain
       ! drop concentration are (is) less than their (its) respective tolerance
       ! amount(s), and are (is) considered to have a value of 0.  There is not
       ! any rain at this grid level.  The correlation is 0 since rain does not
       ! vary at this grid level.
       corr_rrNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 2.
    if ( mu_rr_2 > rr_tol .and. mu_Nr_2 > Nr_tol ) then
       corr_rrNr_2_n = corr_LL2NN( corr_rrNr_2, sigma_rr_2_n, sigma_Nr_2_n )
    else
       ! Mean in-precip rain water mixing ratio and (or) mean in-precip rain
       ! drop concentration are (is) less than their (its) respective tolerance
       ! amount(s), and are (is) considered to have a value of 0.  There is not
       ! any rain at this grid level.  The correlation is 0 since rain does not
       ! vary at this grid level.
       corr_rrNr_2_n = zero
    endif


    return

  end subroutine KK_upscaled_setup

  !=============================================================================
  subroutine KK_stat_output( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
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
                             corr_rrNr_2_n, corr_sw, corr_wrr, corr_wNr, &
                             corr_wNc, level, l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

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
        icorr_sw,       &
        icorr_wrr,      &
        icorr_wNr,      &
        icorr_wNc,      &
        zt

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
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      corr_sw,       & ! Correlation between s & w (both components)         [-]
      corr_wrr,      & ! Correlation between rr & w (both components)        [-]
      corr_wNr,      & ! Correlation between Nr & w (both components)        [-]
      corr_wNc         ! Correlation between Nc & w (both components)        [-]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.


    !!! Output the correlations

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
          call stat_update_var_pt( imu_rr_1_n, level, mu_rr_1_n, zt )
       endif

       ! Mean (in-precip) of ln r_r in PDF component 2.
       if ( imu_rr_2_n > 0 ) then
          call stat_update_var_pt( imu_rr_2_n, level, mu_rr_2_n, zt )
       endif

       ! Mean (in-precip) of ln N_r in PDF component 1.
       if ( imu_Nr_1_n > 0 ) then
          call stat_update_var_pt( imu_Nr_1_n, level, mu_Nr_1_n, zt )
       endif

       ! Mean (in-precip) of ln N_r in PDF component 2.
       if ( imu_Nr_2_n > 0 ) then
          call stat_update_var_pt( imu_Nr_2_n, level, mu_Nr_2_n, zt )
       endif

       ! Mean of ln N_c in PDF component 1.
       if ( imu_Nc_1_n > 0 ) then
          call stat_update_var_pt( imu_Nc_1_n, level, mu_Nc_1_n, zt )
       endif

       ! Mean of ln N_c in PDF component 2.
       if ( imu_Nc_2_n > 0 ) then
          call stat_update_var_pt( imu_Nc_2_n, level, mu_Nc_2_n, zt )
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

       ! Correlation between s and w.
       if ( icorr_sw > 0 ) then
          call stat_update_var_pt( icorr_sw, level, corr_sw, zt )
       endif

       ! Correlation between w and r_r.
       if ( icorr_wrr > 0 ) then
          call stat_update_var_pt( icorr_wrr, level, corr_wrr, zt )
       endif

       ! Correlation between w and N_r.
       if ( icorr_wNr > 0 ) then
          call stat_update_var_pt( icorr_wNr, level, corr_wNr, zt )
       endif

       ! Correlation between w and N_c.
       if ( icorr_wNc > 0 ) then
          call stat_update_var_pt( icorr_wNc, level, corr_wNc, zt )
       endif

    endif


    return

  end subroutine KK_stat_output

  !=============================================================================
  subroutine KK_upscaled_means_driver( rrainm, Nrm, Ncm, &
                                       mu_s_1, mu_s_2, mu_rr_1_n, mu_rr_2_n, &
                                       mu_Nr_1_n, mu_Nr_2_n, mu_Nc_1_n, &
                                       mu_Nc_2_n, sigma_s_1, sigma_s_2, &
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
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       KK_mean_vol_rad &
       = KK_mvr_upscaled_mean( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                               sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                               sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                               KK_mvr_coef, mixt_frac, precip_frac_1, & 
                               precip_frac_2 )

    else  ! r_r or N_r = 0.

       KK_mean_vol_rad = zero

    endif


    return

  end subroutine KK_upscaled_means_driver

  !=============================================================================
  subroutine KK_sed_vel_covars( rrainm, Nrm, KK_mean_vol_rad, &
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                                KK_mvr_coef, mixt_frac, precip_frac_1, &
                                precip_frac_2, level, l_stats_samp, &
                                Vrrprrp, VNrpNrp )

    ! Description:
    ! 

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        micron_per_m, & ! Constant(s)
        rr_tol, &
        Nr_tol, &
        zero

    use KK_upscaled_covariances, only: &
        covar_rr_KK_mvr, & ! Procedure(s)
        covar_Nr_KK_mvr

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use stats_type, only: & 
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        zt,                  & ! Variable(s)
        irr_KK_mvr_covar_zt, &
        iNr_KK_mvr_covar_zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rrainm,          & ! Mean rain water mixing ratio                  [kg/kg]
      Nrm,             & ! Mean rain drop concentration                 [num/kg]
      KK_mean_vol_rad, & ! KK mean volume radius of rain drops               [m]
      mu_rr_1_n,       & ! Mean of ln rr (1st PDF component) in-precip (ip)  [-]
      mu_rr_2_n,       & ! Mean of ln rr (2nd PDF component) ip              [-]
      mu_Nr_1_n,       & ! Mean of ln Nr (1st PDF component) ip              [-]
      mu_Nr_2_n,       & ! Mean of ln Nr (2nd PDF component) ip              [-]
      sigma_rr_1_n,    & ! Standard deviation of ln rr (1st PDF comp.) ip    [-]
      sigma_rr_2_n,    & ! Standard deviation of ln rr (2nd PDF comp.) ip    [-]
      sigma_Nr_1_n,    & ! Standard deviation of ln Nr (1st PDF comp.) ip    [-]
      sigma_Nr_2_n,    & ! Standard deviation of ln Nr (2nd PDF comp.) ip    [-]
      corr_rrNr_1_n,   & ! Corr. betw. ln rr & ln Nr (1st PDF comp.) ip      [-]
      corr_rrNr_2_n,   & ! Corr. betw. ln rr & ln Nr (2nd PDF comp.) ip      [-]
      KK_mvr_coef,     & ! KK mean volume radius coefficient                 [m]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2      ! Precipitation fraction (2nd PDF component)        [-]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      Vrrprrp, & ! Covariance between V_rr and r_r, <V_rr'r_r'>   [(m/s)(kg/kg)]
      VNrpNrp    ! Covariance between V_Nr and N_r, <V_Nr'N_r'>  [(m/s)(num/kg)]

    ! Local Variables
    real( kind = core_rknd ) :: &
      rr_KK_mvr_covar, & ! Covariance of r_r and KK mean vol rad   [(kg/kg)m]
      Nr_KK_mvr_covar    ! Covariance of N_r and KK mean vol rad  [(num/kg)m]


    ! Calculate the covariance between the sedimentation velocity of r_r
    ! and r_r, < V_rr'r_r' >.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       rr_KK_mvr_covar  &
       = covar_rr_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                          sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                          sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                          rrainm, KK_mean_vol_rad, KK_mvr_coef, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       rr_KK_mvr_covar = zero

    endif

    Vrrprrp = - 0.012_core_rknd * micron_per_m * rr_KK_mvr_covar


    ! Calculate the covariance between the sedimentation velocity of N_r
    ! and N_r, < V_Nr'N_r' >.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       Nr_KK_mvr_covar  &
       = covar_Nr_KK_mvr( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                          sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                          sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                          Nrm, KK_mean_vol_rad, KK_mvr_coef, &
                          mixt_frac, precip_frac_1, precip_frac_2 )

    else  ! r_r or N_r = 0.

       Nr_KK_mvr_covar = zero

    endif

    VNrpNrp = - 0.007_core_rknd * micron_per_m * Nr_KK_mvr_covar


    if ( l_stats_samp ) then

       ! Covariance between r_r and KK rain drop mean volume radius.
       if ( irr_KK_mvr_covar_zt > 0 ) then
          call stat_update_var_pt( irr_KK_mvr_covar_zt, level, &
                                   rr_KK_mvr_covar, zt )
       endif

       ! Covariance between N_r and KK rain drop mean volume radius.
       if ( iNr_KK_mvr_covar_zt > 0 ) then
          call stat_update_var_pt( iNr_KK_mvr_covar_zt, level, &
                                   Nr_KK_mvr_covar, zt )
       endif

    endif ! l_stats_samp


    return

  end subroutine KK_sed_vel_covars

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


!===============================================================================
  subroutine KK_init_micro_driver( nz, hydromet, hydromet_mc, hydromet_vel, & ! Intent(in)
                                   hydromet_vel_covar, hydromet_vel_covar_zt, &
                                   rrainm, Nrm, Vrr, VNr, & ! Intent(out)
                                   rrainm_mc_tndcy, Nrm_mc_tndcy, &
                                   KK_auto_tndcy, KK_accr_tndcy, &
                                   Vrrprrp, VNrpNrp, Vrrprrp_zt, &
                                   VNrpNrp_zt, l_src_adj_enabled )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use array_index, only: &
        iirrainm, & ! Constant(s)
        iiNrm

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use constants_clubb, only: &
        zero

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

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet_vel_covar,    & ! Covariance of V_xx & x_x (m-levs)  [units(m/s)]
      hydromet_vel_covar_zt    ! Covariance of V_xx & x_x (t-levs)  [units(m/s)]

    ! Output Variables
    real( kind = core_rknd ), dimension(:), pointer, intent(out) ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
      Vrr,             & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,             & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(:), pointer, intent(out) :: &
      Vrrprrp, & ! Covariance of V_rr and r_r (momentum levels)  [(m/s)(kg/kg)]
      VNrpNrp    ! Covariance of V_Nr and N_r (momentum levels)  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(:), pointer, intent(out) :: &
      Vrrprrp_zt, & ! Covariance of V_rr and r_r; thermo. levs.  [(m/s)(kg/kg)]
      VNrpNrp_zt    ! Covariance of V_Nr and N_r; thermo. levs.  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      KK_auto_tndcy,    & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy       ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]

    logical, intent(out) :: &
      l_src_adj_enabled ! Flag to enable rrainm/Nrm source adjustment

    KK_auto_tndcy = zero
    KK_accr_tndcy = zero

    ! Assign pointers for hydrometeor variables.

    ! Mean fields.
    rrainm => hydromet(:,iirrainm)
    Nrm    => hydromet(:,iiNrm)

    ! Sedimentation Velocities.
    Vrr => hydromet_vel(:,iirrainm)
    VNr => hydromet_vel(:,iiNrm)

    ! Mean field tendencies.
    rrainm_mc_tndcy => hydromet_mc(:,iirrainm)
    Nrm_mc_tndcy    => hydromet_mc(:,iiNrm)

    ! Covariances of hydrometeor sedimentation velocities and their
    ! associated hydrometeors (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp => hydromet_vel_covar(:,iirrainm)
    VNrpNrp => hydromet_vel_covar(:,iiNrm)

    Vrrprrp_zt => hydromet_vel_covar_zt(:,iirrainm)
    VNrpNrp_zt => hydromet_vel_covar_zt(:,iiNrm)

    l_src_adj_enabled = .true.

  end subroutine KK_init_micro_driver
!===============================================================================

!===============================================================================
  subroutine KK_supersaturation( thlm, exner, p_in_Pa, rho, & ! Intent(in)
                                 KK_evap_coef, KK_auto_coef, & ! Intent(out)
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
      thlm,       & ! Mean liquid water potential temperature         [K]
      p_in_Pa,    & ! Pressure                                        [Pa]
      exner,      & ! Exner function                                  [-]
      rho           ! Density                                         [kg/m^3]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    ! Local Variables
    real( kind = core_rknd ) :: &
      T_liq_in_K, & ! Mean liquid water temperature, T_l           [K]
      r_sl,       & ! Liquid water sat. mixing ratio, r_s(T_l,p)   [kg/kg]
      Beta_Tl       ! Parameter Beta, Beta(T_l)                    [1/(kg/kg)]

    ! Compute supersaturation via s1, s2.
    !     Larson et al 2002, JAS, Vol 59, p 3534.
    ! This allows a more direct comparison of local, upscaled formulas.

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

  end subroutine KK_supersaturation
!===============================================================================

!===============================================================================
  subroutine KK_microphys_adjust( rrainm, Nrm, rcm, exner, & ! Intent(in)
                                 KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                                 dt, l_src_adj_enabled, &
                                 rrainm_source, Nrm_source, & ! Intent(out)
                                 KK_Nrm_auto_tndcy, KK_Nrm_evap_tndcy, &
                                 rrainm_src_adj, Nrm_src_adj, &
                                 rrainm_evap_net, Nrm_evap_net, &
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
        KK_Nrm_evap, & ! Procedure(s)
        KK_Nrm_auto

    implicit none

    ! Input Variables

    real( kind = time_precision ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    real( kind = core_rknd ), intent(in) :: &
      exner,      & ! Exner function                                  [-]
      rcm           ! Mean cloud water mixing ratio                   [kg/kg]

    real( kind = core_rknd ), intent(in) ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm                ! Mean rain drop concentration, < N_r >    [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_tndcy, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_tndcy    ! KK accretion coefficient                    [(kg/kg)/s]

    logical, intent(in) :: &
      l_src_adj_enabled   ! Flag to enable rrainm/Nrm source adjustment

    ! Output Variables

    real( kind = core_rknd ), intent(out) ::  &
      rrainm_source,     & ! Total source term rate for rrainm       [(kg/kg)/s]
      Nrm_source,        & ! Total source term rate for Nrm         [(num/kg)/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), intent(out) :: &
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]

    real( kind = core_rknd ), intent(out) ::  &
      rrainm_src_adj,  & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj,     & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrainm_evap_net, & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net       ! Net evaporation rate of <N_r>            [(num/kg)/s]

    real( kind = core_rknd ), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    ! Local Variables
    real( kind = core_rknd ) :: &
      rrainm_src_max,    & ! Maximum allowable rrainm source rate    [(kg/kg)/s]
      rrainm_auto_ratio, & ! Ratio of rrainm autoconv to overall source term [-]
      total_rc_needed      ! Amount of r_c needed to over the timestep
                           ! for rain source terms                       [kg/kg]
    ! ---- Begin Code ----

    !!! KK rain drop concentration microphysics tendencies.

    !!! Calculate the KK N_r evaporation tendency.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       KK_Nrm_evap_tndcy  &
       = KK_Nrm_evap( KK_evap_tndcy, Nrm, rrainm )

    else  ! r_r or N_r = 0.

       KK_Nrm_evap_tndcy = zero

    endif

    !!! Calculate the KK N_r autoconversion tendency.
    KK_Nrm_auto_tndcy = KK_Nrm_auto( KK_auto_tndcy )

    !!! Source-adjustment code for rrainm and Nrm.

    rrainm_source = KK_auto_tndcy + KK_accr_tndcy
    Nrm_source = KK_Nrm_auto_tndcy

    ! The increase of rain due to autoconversion and accretion both draw
    ! their water from the available cloud water.  Over a long time step
    ! these rates may over-deplete cloud water.  In other words, these
    ! processes may draw more cloud water than there is available.  Thus,
    ! the total source rate multiplied by the time step length cannot exceed
    ! the total amount of cloud water available.  If it does, then the rate
    ! must be adjusted.
    total_rc_needed = rrainm_source * real( dt, kind = core_rknd )

    if ( total_rc_needed > rcm .and. l_src_adj_enabled ) then

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
       ! an autoconversion term for a source term.  Figure that change in the
       ! rrainm autoconversion term is proportional to to the total rrainm
       ! adjustment rate by the ratio of rrainm autoconversion to the overall
       ! source term.  Then, plug the rrainm autoconversion adjustment into
       ! the equation for Nrm autoconversion to determine the effect on the
       ! Nrm source term.
       rrainm_auto_ratio = KK_auto_tndcy /  &
                           ( KK_auto_tndcy + KK_accr_tndcy )
       Nrm_src_adj = KK_Nrm_auto( rrainm_auto_ratio * rrainm_src_adj )

       ! Change Nrm by Nrm_src_adj.  Nrm_src_adj will always be negative.
       Nrm_source = Nrm_source + Nrm_src_adj

    else

       rrainm_src_adj = zero
       Nrm_src_adj    = zero

    endif

    ! Prevent over-evaporation of rain over a long model time step.
    ! Limit the evaporation rate.  The total amount of rain lost due to
    ! evaporation cannot be so great as to result in negative rain.
    ! Calculate net evaporation rate of <r_r>.
    rrainm_evap_net = max( KK_evap_tndcy, &
                              - rrainm / real( dt, kind = core_rknd ) )

    ! Recalcuate the net evaporation rate of <N_r> based on the net
    ! evaporation rate of <r_r>.
    if ( KK_evap_tndcy /= rrainm_evap_net .and. &
         rrainm > rr_tol .and. Nrm > Nr_tol ) then
       Nrm_evap_net = KK_Nrm_evap( rrainm_evap_net, Nrm, rrainm )
    else
       Nrm_evap_net = KK_Nrm_evap_tndcy
    endif

    Nrm_evap_net = max( Nrm_evap_net, &
                           - Nrm / real( dt, kind = core_rknd ) )


    !!! Calculate overall KK microphysics tendencies.
    rrainm_mc_tndcy = rrainm_evap_net + rrainm_source
    Nrm_mc_tndcy    = Nrm_evap_net + Nrm_source

    !!! Explicit contributions to thlm and rtm from the microphysics
    rvm_mc  = -rrainm_evap_net
    rcm_mc  = -rrainm_source  ! Accretion + Autoconversion
    thlm_mc = ( Lv / ( Cp * exner ) ) * rrainm_mc_tndcy

  end subroutine KK_microphys_adjust
!===============================================================================

!===============================================================================
  subroutine KK_stats_output_samp_in_sub( KK_mean_vol_rad, KK_evap_tndcy, & ! Intent(in)
                                        KK_auto_tndcy, KK_accr_tndcy, &
                                        KK_Nrm_evap_tndcy, KK_Nrm_auto_tndcy, &
                                        rrainm_src_adj, Nrm_src_adj, &
                                        rrainm_evap_net, Nrm_evap_net, &
                                        k )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use stats_variables, only: &
        zt,               & ! Variable(s)
        im_vol_rad_rain,  &
        irrainm_cond,     &
        irrainm_auto,     &
        irrainm_accr,     &
        irrainm_src_adj,  &
        irrainm_cond_adj, &
        iNrm_cond,        &
        iNrm_auto,        &
        iNrm_src_adj,     &
        iNrm_cond_adj,    &
        iprecip_frac

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)


    implicit none

    ! Input Variables

    integer, intent(in) :: k ! height level

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy,     & ! Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
      KK_mean_vol_rad,   & ! Mean KK rain drop mean volume radius     [m]
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
      KK_Nrm_auto_tndcy, & ! Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]
      KK_auto_tndcy,     & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]
      rrainm_src_adj,    & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj,       & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrainm_evap_net,   & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net         ! Net evaporation rate of <N_r>            [(num/kg)/s]


    ! Rain drop mean volume radius.
    call stat_update_var_pt( im_vol_rad_rain, k, KK_mean_vol_rad, zt )

    ! Explicit contributions to rrainm.
    call stat_update_var_pt( irrainm_cond, k, KK_evap_tndcy, zt )

    call stat_update_var_pt( irrainm_auto, k, KK_auto_tndcy, zt )

    call stat_update_var_pt( irrainm_accr, k, KK_accr_tndcy, zt )

    ! Explicit contributions to Nrm.
    call stat_update_var_pt( iNrm_cond, k, KK_Nrm_evap_tndcy, zt )

    call stat_update_var_pt( iNrm_auto, k, KK_Nrm_auto_tndcy, zt )


    call stat_update_var_pt( irrainm_src_adj, k, rrainm_src_adj, zt )

    call stat_update_var_pt( iNrm_src_adj, k, Nrm_src_adj, zt )

    call stat_update_var_pt( irrainm_cond_adj, k, &
                             rrainm_evap_net - KK_evap_tndcy, zt )

    call stat_update_var_pt( iNrm_cond_adj, k, &
                             Nrm_evap_net - KK_Nrm_evap_tndcy, zt )

  end subroutine KK_stats_output_samp_in_sub
!===============================================================================

subroutine KK_sedimentation( nz, & ! Intent(in)
                           Vrr, VNr, & ! Intent(InOut)
                           rrainm_mc_tndcy, Nrm_mc_tndcy, &
                           rcm_mc, rvm_mc, thlm_mc, &
                           Vrrprrp, VNrpNrp, &
                           Vrrprrp_zt, VNrpNrp_zt, &
                           KK_mean_vol_rad )

  ! Description:
  !
  ! Bogus example
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
        micron_per_m  ! Constant(s)

    use clubb_precision, only: &
        core_rknd       ! Variable(s)

    use constants_clubb, only: &
        zero ! Constant(s)

    use grid_class, only: &
        zt2zm ! Procedure(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(:), pointer, intent(inout) ::  &
      Vrr,             & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,             & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(:), pointer, intent(inout) :: &
      Vrrprrp, & ! Covariance of V_rr and r_r (momentum levels)  [(m/s)(kg/kg)]
      VNrpNrp    ! Covariance of V_Nr and N_r (momentum levels)  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(:), pointer, intent(inout) :: &
      Vrrprrp_zt, & ! Covariance of V_rr and r_r; thermo. levs.  [(m/s)(kg/kg)]
      VNrpNrp_zt    ! Covariance of V_Nr and N_r; thermo. levs.  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      KK_mean_vol_rad    ! Mean KK rain drop mean volume radius     [m]

    ! Local Variables

    integer :: k ! Loop iterator


  !-----------------------------------------------------------------------

    !----- Begin Code -----


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

    !!! Sedimentation velocities
    forall ( k = 1:nz-1 )

       ! Sedimentation velocity of rrainm.
!       Vrr(k) = 0.012_core_rknd * ( micron_per_m * zt2zm(KK_mean_vol_rad,k) ) &
!                - 0.2_core_rknd
       Vrr(k) = 0.012_core_rknd * ( micron_per_m * KK_mean_vol_rad(k) )  &
                - 0.2_core_rknd

       ! Sedimentation velocity is positive upwards.
       Vrr(k) = -max( Vrr(k), zero )

       ! Sedimentation velocity of Nrm.
!       VNr(k) = 0.007_core_rknd * ( micron_per_m * zt2zm(KK_mean_vol_rad,k) ) &
!                - 0.1_core_rknd
       VNr(k) = 0.007_core_rknd * ( micron_per_m * KK_mean_vol_rad(k) )  &
                - 0.1_core_rknd

       ! Sedimentation velocity is positive upwards.
       VNr(k) = -max( VNr(k), zero )

    end forall ! 1..nz-1

    !!! Boundary conditions for sedimentation velocities.

    ! The flux of rain water through the model top is 0.
    ! Vrr and VNr are set to 0 at the highest model level.
    Vrr(nz) = zero
    VNr(nz) = zero

    !!! Boundary conditions (lower) for the covariances of hydrometeor
    !!! sedimentation velocities and their associated hydrometeors
    !!! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp_zt(1) = Vrrprrp_zt(2)
    VNrpNrp_zt(1) = VNrpNrp_zt(2)

    !!! Interpolate the covariances of hydrometeor sedimentation velocities
    !!! and their associated hydrometeors (<V_rr'r_r'> and <V_Nr'N_r'>) to
    !!! momentum levels.
    Vrrprrp = zt2zm( Vrrprrp_zt )
    VNrpNrp = zt2zm( VNrpNrp_zt )

    !!! Boundary conditions (upper) for the covariances of hydrometeor
    !!! sedimentation velocities and their associated hydrometeors
    !!! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp(nz) = zero
    VNrpNrp(nz) = zero


  return
end subroutine KK_sedimentation
!-----------------------------------------------------------------------

end module KK_microphys_module

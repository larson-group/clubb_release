! $Id$
!===============================================================================
module KK_microphys_module

  implicit none

  private

  public :: KK_micro_driver

  contains

  !=============================================================================
!  subroutine KK_micro_driver( dt, l_local_kk, thlm, rho, p_in_Pa, &
!                              exner, s_mellor, rcm, hydromet, &
!                              pdf_params, hydromet_mc, hydromet_vel, &
!                              rcm_mc, rvm_mc, thlm_mc )
  subroutine KK_micro_driver( dt, nz, l_stats_samp, l_local_kk, &
                              l_latin_hypercube, thlm, p_in_Pa, exner, rho, &
                              cloud_frac, pdf_params, wm, w_std_dev, dzq, rcm, s_mellor, &
                              rvm, hydromet, hydromet_mc, hydromet_vel, &
                              rcm_mc, rvm_mc, thlm_mc )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

!    use grid_class, only: &
!        gr  ! Variable(s)

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
        micron_per_m

    use parameters_microphys, only: &
        rrp2_on_rrainm2_cloud, & ! Constant(s)
        rrp2_on_rrainm2_below, &
        Nrp2_on_Nrm2_cloud,    &
        Nrp2_on_Nrm2_below,    &
        Ncp2_on_Ncm2_cloud,    &
        Ncp2_on_Ncm2_below,    &
!        corr_srr_NL_cloud,     &
!        corr_srr_NL_below,     &
!        corr_sNr_NL_cloud,     &
!        corr_sNr_NL_below,     &
!        corr_sNc_NL_cloud,     &
!        corr_sNc_NL_below,     &
        corr_rrNr_LL_cloud,    &
        corr_rrNr_LL_below,    &
        corr_rtrr_NL_cloud,    &
        corr_rtrr_NL_below,    &
        corr_rtNr_NL_cloud,    &
        corr_rtNr_NL_below,    &
        corr_rtNc_NL_cloud,    &
        corr_rtNc_NL_below,    &
        corr_thlrr_NL_cloud,   &
        corr_thlrr_NL_below,   &
        corr_thlNr_NL_cloud,   &
        corr_thlNr_NL_below,   &
        corr_thlNc_NL_cloud,   &
        corr_thlNc_NL_below,   &
        KK_auto_Nc_exp,        &
        C_evap

    use KK_utilities, only: &
        mean_L2N,     & ! Procedure(s)
        stdev_L2N,    &
        corr_NL2NN,   &
        corr_LL2NN,   &
        calc_corr_sx, &
        G_T_p

    use KK_upscaled_means, only: &
        KK_mvr_upscaled_mean,  & ! Procedure(s)
        KK_evap_upscaled_mean, &
        KK_auto_upscaled_mean, &
        KK_accr_upscaled_mean

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
        iiNcm,    &
        iiNrm

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type, only: & 
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        zt,                 & ! Variable(s)
!        l_stats_samp,       &
        im_vol_rad_rain, &
        irrainm_cond,       &
        irrainm_auto,       &
        irrainm_accr,       &
        irrainm_src_adj,    &
        iNrm_cond,          &
        iNrm_auto,          &
        iNrm_src_adj

    implicit none

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      nz        ! Number of model vertical grid levels

    logical, intent(in) :: &
      l_stats_samp,      & ! Flag to sample statistics
      l_local_kk,        & ! Flag to use the local form of KK microphysics
      l_latin_hypercube    ! Flag to use Latin Hypercube interface

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      thlm,       & ! Mean liquid water potential temperature  [K]
      rho,        & ! Density                                  [kg/m^3]
      cloud_frac, & ! Cloud fraction                           [-]
      p_in_Pa,    & ! Pressure                                 [Pa]
      exner,      & ! Exner function                           [-]
      s_mellor,   & ! Mean extended liquid water mixing ratio  [kg/kg]
      rcm           ! Mean cloud water mixing ratio            [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wm,        & ! Mean vertical velocity, w (for LH interface)        [m/s] 
      w_std_dev, & ! Standard deviation of w (for LH interface)          [m/s]
      dzq,       & ! Thickness between thermo. levels (for LH interface) [m]
      rvm          ! Mean water vapor mixing ratio (for LH interface)    [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet    ! Hydrometeor species                      [units vary]

    type(pdf_parameter), dimension(nz), target, intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

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

    ! Local Variables
    real( kind = core_rknd ), dimension(:), pointer ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
      Ncm,             & ! Mean cloud droplet conc., < N_c >        [num/kg]
      Vrr,             & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,             & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) :: &
      rrainm_cond,  & ! Mean KK (dr_r/dt) due to evaporation      [(kg/kg)/s]
      rrainm_auto,  & ! Mean KK (dr_r/dt) due to autoconversion   [(kg/kg)/s]
      rrainm_accr,  & ! Mean KK (dr_r/dt) due to accretion        [(kg/kg)/s]
      Nrm_cond,     & ! Mean KK (dN_r/dt) due to condensation     [(num/kg)/s]
      Nrm_auto,     & ! Mean KK (dN_r/dt) due to autoconversion   [(num/kg)/s]
      mean_vol_rad    ! Mean KK rain drop mean volume radius      [m]

    real( kind = core_rknd ) :: &
      T_liq_in_K, & ! Mean liquid water temperature, T_l                [K]
      r_sl,       & ! Liquid water saturation mixing ratio, r_s(T_l,p)  [kg/kg]
      Beta_Tl       ! Parameter Beta, Beta(T_l)                         [-]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    real( kind = core_rknd ), dimension(:), pointer :: &
      varnce_rt1,  & ! Variance of r_t (1st PDF component)           [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t (2nd PDF component)           [kg^2/kg^2]
      varnce_thl1, & ! Variance of theta_l (1st PDF component)       [K^2]
      varnce_thl2, & ! Variance of theta_l (2nd PDF component)       [K^2]
      crt1,        & ! Coefficient of r_t for s' (1st PDF component) [-]
      crt2,        & ! Coefficient of r_t for s' (2nd PDF component) [-]
      cthl1,       & ! Coefficient of th_l for s' (1st PDF comp.)    [1/K]
      cthl2,       & ! Coefficient of th_l for s' (2nd PDF comp.)    [1/K]
      s1,          & ! Mean of s (1st PDF component)                 [kg/kg]
      s2,          & ! Mean of s (2nd PDF component)                 [kg/kg]
      stdev_s1,    & ! Standard deviation of s (1st PDF component)   [kg/kg]
      stdev_s2,    & ! Standard deviation of s (2nd PDF component)   [kg/kg]
      mixt_frac      ! Mixture fraction                              [-]

    real( kind = core_rknd ) :: &
      mu_rr_n,      & ! Mean of ln rr (both components)                     [-]
      mu_Nr_n,      & ! Mean of ln Nr (both components)                     [-]
      mu_Nc_n,      & ! Mean of ln Nc (both components)                     [-]
      sigma_rr_n,   & ! Standard deviation of ln rr (both components)       [-]
      sigma_Nr_n,   & ! Standard deviation of ln Nr (both components)       [-]
      sigma_Nc_n,   & ! Standard deviation of ln Nc (both components)       [-]
      corr_srr_1_n, & ! Correlation between s and ln rr (1st PDF component) [-]
      corr_srr_2_n, & ! Correlation between s and ln rr (2nd PDF component) [-]
      corr_sNr_1_n, & ! Correlation between s and ln Nr (1st PDF component) [-]
      corr_sNr_2_n, & ! Correlation between s and ln Nr (2nd PDF component) [-]
      corr_sNc_1_n, & ! Correlation between s and ln Nc (1st PDF component) [-]
      corr_sNc_2_n, & ! Correlation between s and ln Nc (2nd PDF component) [-]
      corr_rrNr_n     ! Correlation between ln rr & ln Nr (both components) [-]

    real( kind = core_rknd ) :: &
      sigma_rt_1,  & ! Standard deviation of r_t (1st PDF component)     [kg/kg]
      sigma_rt_2,  & ! Standard deviation of r_t (2nd PDF component)     [kg/kg]
      sigma_thl_1, & ! Standard deviation of theta_l (1st PDF component) [K]
      sigma_thl_2, & ! Standard deviation of theta_l (2nd PDF component) [K]
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component)  [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component)  [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component)  [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component)  [-]
      corr_sNc_1,  & ! Correlation between s and Nc (1st PDF component)  [-]
      corr_sNc_2     ! Correlation between s and Nc (2nd PDF component)  [-]

    real( kind = core_rknd ), dimension(nz) :: &
      rrp2_on_rrainm2, & ! Specified ratio of < r_r >^2 to < r_r'^2 >       [-]
      Nrp2_on_Nrm2,    & ! Specified ratio of < N_r >^2 to < N_r'^2 >       [-]
      Ncp2_on_Ncm2,    & ! Specified ratio of < N_c >^2 to < N_c'^2 >       [-]
      corr_rrNr_LL,    & ! Specified correlation between r_r and N_r        [-]
      corr_rtrr_NL,    & ! Specified correlation between r_t and r_r        [-]
      corr_rtNr_NL,    & ! Specified correlation between r_t and N_r        [-]
      corr_rtNc_NL,    & ! Specified correlation between r_t and N_c        [-]
      corr_thlrr_NL,   & ! Specified correlation between theta_l and r_r    [-]
      corr_thlNr_NL,   & ! Specified correlation between theta_l and N_r    [-]
      corr_thlNc_NL      ! Specified correlation between theta_l and N_c    [-]
!      corr_srr_NL,     & ! Specified correlation between s and r_r          [-]
!      corr_sNr_NL,     & ! Specified correlation between s and N_r          [-]
!      corr_sNc_NL        ! Specified correlation between s and N_c          [-]

!      ! Normalized Correlations
!      corr_rtrr_n
!      corr_thlrr_n
!      corr_rtNr_n
!      corr_thlNr_n
!      corr_rtNc_n
!      corr_thlNc_n

    real( kind = core_rknd ) ::  &
      rrainm_source,     & ! Total source term rate for rrainm       [(kg/kg)/s]
      Nrm_source,        & ! Total source term rate for Nrm         [(num/kg)/s]
      rrainm_src_max,    & ! Maximum allowable rrainm source rate    [(kg/kg)/s]
      rrainm_auto_ratio, & ! Ratio of rrainm autoconv to overall source term [-]
      total_rc_needed      ! Amount of r_c needed to over the timestep
                           ! for rain source terms                       [kg/kg]

    real( kind = core_rknd ), dimension(nz) ::  &
      rrainm_src_adj, & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj       ! Total adjustment to Nrm source terms     [{num/kg)/s]

    logical :: &
      l_upscaled,        & ! Flag for using upscaled KK microphysics.
      l_src_adj_enabled    ! Flag to enable rrainm/Nrm source adjustment

    integer :: &
      k   ! Loop index


    ! Remove compiler warnings
    if ( .false. .and. l_latin_hypercube ) then
      rrainm_src_adj = dzq
      rrainm_src_adj = rvm
      rrainm_src_adj = w_std_dev
      rrainm_src_adj = wm
      rrainm_src_adj = cloud_frac
    end if

    ! Assign pointers for hydrometeor variables.

    ! Mean fields.
    rrainm => hydromet(:,iirrainm)
    Nrm    => hydromet(:,iiNrm)
    Ncm    => hydromet(:,iiNcm)

    ! Sedimentation Velocities.
    Vrr => hydromet_vel(:,iirrainm)
    VNr => hydromet_vel(:,iiNrm)
    hydromet_vel(:,iiNcm) = zero

    ! Mean field tendencies.
    rrainm_mc_tndcy => hydromet_mc(:,iirrainm)
    Nrm_mc_tndcy    => hydromet_mc(:,iiNrm)
    
    ! Assign pointers for PDF parameters.
    ! Note:  these are only necessary for upscaled KK microphysics; however,
    !        they need to be initialized here to avoid a compiler warning.
    varnce_rt1  => pdf_params(:)%varnce_rt1
    varnce_rt2  => pdf_params(:)%varnce_rt2
    varnce_thl1 => pdf_params(:)%varnce_thl1
    varnce_thl2 => pdf_params(:)%varnce_thl2
    crt1        => pdf_params(:)%crt1
    crt2        => pdf_params(:)%crt2
    cthl1       => pdf_params(:)%cthl1
    cthl2       => pdf_params(:)%cthl2
    mixt_frac   => pdf_params(:)%mixt_frac
    s1          => pdf_params(:)%s1
    s2          => pdf_params(:)%s2
    stdev_s1    => pdf_params(:)%stdev_s1
    stdev_s2    => pdf_params(:)%stdev_s2


    if ( .not. l_local_kk ) then
       l_upscaled = .true.
    else
       l_upscaled = .false.
    endif

    l_src_adj_enabled = .true.


    if ( l_upscaled ) then

       ! Set up the values of the statistical correlations and variances.  Since
       ! we currently do not have enough variables to compute the correlations
       ! and variances directly, we have obtained these values by analyzing LES
       ! runs of certain cases.  We have divided those results into an
       ! inside-cloud average and an outside-cloud (or below-cloud) average.
       ! This coding leaves the software architecture in place in case we ever
       ! have the variables in place to compute these values directly.  It also
       ! allows us to use separate inside-cloud and outside-cloud parameter
       ! values.
       ! Brian Griffin; February 3, 2007.
       !
       ! Set the value of the parameters based on whether the altitude is above
       ! or below cloud base.  Determine whether there is cloud at any given
       ! vertical level.  In order for a vertical level to have cloud, the
       ! amount of cloud water (rcm) must be greater than or equal to the
       ! tolerance level (rc_tol).  If there is cloud at a given vertical level,
       ! then the ###_cloud value is used.  Otherwise, the ###_below value is
       ! used.
       where ( rcm >= rc_tol )
          rrp2_on_rrainm2 = rrp2_on_rrainm2_cloud
          Nrp2_on_Nrm2    = Nrp2_on_Nrm2_cloud
          Ncp2_on_Ncm2    = Ncp2_on_Ncm2_cloud
          corr_rrNr_LL    = corr_rrNr_LL_cloud
          corr_rtrr_NL    = corr_rtrr_NL_cloud
          corr_rtNr_NL    = corr_rtNr_NL_cloud
          corr_rtNc_NL    = corr_rtNc_NL_cloud
          corr_thlrr_NL   = corr_thlrr_NL_cloud
          corr_thlNr_NL   = corr_thlNr_NL_cloud
          corr_thlNc_NL   = corr_thlNc_NL_cloud
!          corr_srr_NL     = corr_srr_NL_cloud
!          corr_sNr_NL     = corr_sNr_NL_cloud
!          corr_sNc_NL     = corr_sNc_NL_cloud
       else where
          rrp2_on_rrainm2 = rrp2_on_rrainm2_below
          Nrp2_on_Nrm2    = Nrp2_on_Nrm2_below
          Ncp2_on_Ncm2    = Ncp2_on_Ncm2_below
          corr_rrNr_LL    = corr_rrNr_LL_below
          corr_rtrr_NL    = corr_rtrr_NL_below
          corr_rtNr_NL    = corr_rtNr_NL_below
          corr_rtNc_NL    = corr_rtNc_NL_below
          corr_thlrr_NL   = corr_thlrr_NL_below
          corr_thlNr_NL   = corr_thlNr_NL_below
          corr_thlNc_NL   = corr_thlNc_NL_below
!          corr_srr_NL     = corr_srr_NL_below
!          corr_sNr_NL     = corr_sNr_NL_below
!          corr_sNc_NL     = corr_sNc_NL_below
       end where

    endif  ! l_upscaled


    ! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz-1, 1

       ! Compute supersaturation via s1, s2.
       !     Larson et al 2002, JAS, Vol 59, p 3534.
       ! This allows a more direct comparison of local, upscaled formulas.

       ! Liquid water temperature.
       T_liq_in_K = thlm(k) * exner(k)

       ! Saturation mixing ratio (based on liquid water temperature and
       ! pressure), r_sl = r_s(T_l,p).
       r_sl = sat_mixrat_liq( p_in_Pa(k), T_liq_in_K )

       ! Beta(T_l).
       Beta_Tl = (Rd/Rv) * ( Lv / ( Rd * T_liq_in_K ) )  &
                         * ( Lv / ( Cp * T_liq_in_K ) )

       ! Coefficient for KK evaporation.
       KK_evap_coef = three * C_evap * G_T_p( T_liq_in_K, p_in_Pa(k) )   &
                            * ( four_thirds * pi * rho_lw )**two_thirds  &
                            * ( ( one + Beta_Tl * r_sl ) / r_sl )

       ! Coefficient for KK autoconversion.
       KK_auto_coef = 1350.0_core_rknd * ( rho(k) / cm3_per_m3 )**KK_auto_Nc_exp

       ! Coefficient for KK accretion.
       KK_accr_coef = 67.0_core_rknd

       ! Coefficient for KK rain drop mean volume radius.
       KK_mvr_coef = ( four_thirds * pi * rho_lw )**(-one_third)


       !!! KK rain water mixing ratio microphysics tendencies.
       if ( l_upscaled ) then


          !!! Calculate the normalized mean of variables that have an assumed
          !!! (single) lognormal distribution, given the mean and variance of
          !!! those variables.

          ! Normalized mean of rain water mixing ratio.
          if ( rrainm(k) > rr_tol ) then
             mu_rr_n = mean_L2N( rrainm(k), &
                                 rrp2_on_rrainm2(k) * rrainm(k)**2 )
          endif

          ! Normalized mean of rain drop concentration.
          if ( Nrm(k) > Nr_tol ) then
             mu_Nr_n = mean_L2N( Nrm(k), Nrp2_on_Nrm2(k) * Nrm(k)**2 )
          endif

          ! Normalized mean of cloud droplet concentration.
          if ( Ncm(k) > Nc_tol ) then
             mu_Nc_n = mean_L2N( Ncm(k), Ncp2_on_Ncm2(k) * Ncm(k)**2 )
          endif

          !!! Calculate the standard deviation of variables that have an assumed
          !!! normal distribution for the ith PDF component.

          ! Standard deviation of liquid water potential temperature in PDF
          ! component 1.
          sigma_thl_1 = sqrt( varnce_thl1(k) )

          ! Standard deviation of liquid water potential temperature in PDF
          ! component 2.
          sigma_thl_2 = sqrt( varnce_thl2(k) )

          ! Standard deviation of total water mixing ratio in PDF component 1.
          sigma_rt_1 = sqrt( varnce_rt1(k) )

          ! Standard deviation of total water mixing ratio in PDF component 2.
          sigma_rt_2 = sqrt( varnce_rt2(k) )

          !!! Calculate the normalized standard deviation of variables that have
          !!! an assumed (single) lognormal distribution, given the mean and
          !!! variance of those variables.

          ! Normalized standard deviation of rain water mixing ratio.
          if ( rrainm(k) > rr_tol ) then
             sigma_rr_n = stdev_L2N( rrainm(k), &
                                     rrp2_on_rrainm2(k) * rrainm(k)**2 )
          endif

          ! Normalized standard deviation of rain drop concentration.
          if ( Nrm(k) > Nr_tol ) then
             sigma_Nr_n = stdev_L2N( Nrm(k), Nrp2_on_Nrm2(k) * Nrm(k)**2 )
          endif

          ! Normalized standard deviation of cloud droplet concentration.
          if ( Ncm(k) > Nc_tol ) then
             sigma_Nc_n = stdev_L2N( Ncm(k), Ncp2_on_Ncm2(k) * Ncm(k)**2 )
          endif

          ! Note:  the standard deviation of extended liquid water mixing ratio,
          !        s, is given by stdev_s1 for PDF component 1 and stdev_s2 for
          !        PDF component 2.

          !!! Calculate correlations between extended liquid water mixing ratio
          !!! and another variable.

          ! Calculate the correlation between s and r_r in PDF component 1.
          corr_srr_1  &
          = calc_corr_sx( crt1(k), cthl1(k), sigma_rt_1, sigma_thl_1,  &
                          stdev_s1(k), corr_rtrr_NL(k), corr_thlrr_NL(k) )

          ! Calculate the correlation between s and r_r in PDF component 2.
          corr_srr_2  &
          = calc_corr_sx( crt2(k), cthl2(k), sigma_rt_2, sigma_thl_2,  &
                          stdev_s2(k), corr_rtrr_NL(k), corr_thlrr_NL(k) )

          ! Calculate the correlation between s and N_r in PDF component 1.
          corr_sNr_1  &
          = calc_corr_sx( crt1(k), cthl1(k), sigma_rt_1, sigma_thl_1,  &
                          stdev_s1(k), corr_rtNr_NL(k), corr_thlNr_NL(k) )

          ! Calculate the correlation between s and N_r in PDF component 2.
          corr_sNr_2  &
          = calc_corr_sx( crt2(k), cthl2(k), sigma_rt_2, sigma_thl_2,  &
                          stdev_s2(k), corr_rtNr_NL(k), corr_thlNr_NL(k) )

          ! Calculate the correlation between s and N_c in PDF component 1.
          corr_sNc_1  &
          = calc_corr_sx( crt1(k), cthl1(k), sigma_rt_1, sigma_thl_1,  &
                          stdev_s1(k), corr_rtNc_NL(k), corr_thlNc_NL(k) )

          ! Calculate the correlation between s and N_c in PDF component 2.
          corr_sNc_2  &
          = calc_corr_sx( crt2(k), cthl2(k), sigma_rt_2, sigma_thl_2,  &
                          stdev_s2(k), corr_rtNc_NL(k), corr_thlNc_NL(k) )

          !!! Calculate the normalized correlation between variables that have
          !!! an assumed normal distribution and variables that have an assumed
          !!! (single) lognormal distribution for the ith PDF component, given
          !!! their correlation and the normalized standard deviation of the
          !!! variable with the assumed lognormal distribution.

          if ( rrainm(k) > rr_tol ) then

             ! Normalize the correlation between s and r_r in PDF component 1.
             corr_srr_1_n = corr_NL2NN( corr_srr_1, sigma_rr_n )

             ! Normalize the correlation between s and r_r in PDF component 2.
             corr_srr_2_n = corr_NL2NN( corr_srr_2, sigma_rr_n )

           endif

          if ( Nrm(k) > Nr_tol ) then

             ! Normalize the correlation between s and N_r in PDF component 1.
             corr_sNr_1_n = corr_NL2NN( corr_sNr_1, sigma_Nr_n )

             ! Normalize the correlation between s and N_r in PDF component 2.
             corr_sNr_2_n = corr_NL2NN( corr_sNr_2, sigma_Nr_n )

          endif

          if ( Ncm(k) > Nc_tol ) then

             ! Normalize the correlation between s and N_c in PDF component 1.
             corr_sNc_1_n = corr_NL2NN( corr_sNc_1, sigma_Nc_n )

             ! Normalize the correlation between s and N_c in PDF component 2.
             corr_sNc_2_n = corr_NL2NN( corr_sNc_2, sigma_Nc_n )

          endif

!          if ( rrainm(k) > rr_tol ) then
!
!             ! Normalize the correlation between s and r_r in PDF component 1.
!             corr_srr_1_n = corr_NL2NN( corr_srr_NL(k), sigma_rr_n )
!
!             ! Normalize the correlation between s and r_r in PDF component 2.
!             corr_srr_2_n = corr_srr_1_n
!
!          endif
!
!          if ( Nrm(k) > Nr_tol ) then
!
!             ! Normalize the correlation between s and N_r in PDF component 1.
!             corr_sNr_1_n = corr_NL2NN( corr_sNr_NL(k), sigma_Nr_n )
!
!             ! Normalize the correlation between s and N_r in PDF component 2.
!             corr_sNr_2_n = corr_sNr_1_n
!
!          endif
!
!          if ( Ncm(k) > Nc_tol ) then
!
!             ! Normalize the correlation between s and N_c in PDF component 1.
!             corr_sNc_1_n = corr_NL2NN( corr_sNc_NL(k), sigma_Nc_n )
!
!             ! Normalize the correlation between s and N_c in PDF component 2.
!             corr_sNc_2_n = corr_sNc_1_n
!
!          endif

          !!! Calculate the normalized correlation between two variables that
          !!! both have an assumed lognormal distribution, given their
          !!! correlation and both of their normalized standard deviations.

          ! Normalize the correlation between rr and Nr (this is the same for
          ! both PDF components).
          if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then
             corr_rrNr_n = corr_LL2NN( corr_rrNr_LL(k), sigma_rr_n, sigma_Nr_n )
          endif


          !!! Calculate the upscaled KK rain drop mean volume radius.
          if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

             mean_vol_rad(k) &
             = KK_mvr_upscaled_mean( mu_rr_n, mu_Nr_n, sigma_rr_n, &
                                     sigma_Nr_n, corr_rrNr_n, KK_mvr_coef )

          else  ! r_r or N_r = 0.

             mean_vol_rad(k) = zero

          endif

          !!! Calculate the values of the upscaled KK microphysics tendencies.

          !!! Calculate the upscaled KK evaporation tendency.
          if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

             rrainm_cond(k)  &
             = KK_evap_upscaled_mean( s1(k), s2(k), mu_rr_n, mu_Nr_n, &
                                      stdev_s1(k), stdev_s2(k), sigma_rr_n, &
                                      sigma_Nr_n, corr_srr_1_n, corr_srr_2_n, &
                                      corr_sNr_1_n, corr_sNr_2_n, corr_rrNr_n, &
                                      KK_evap_coef, mixt_frac(k) )

          else  ! r_r or N_r = 0.

             rrainm_cond(k) = zero

          endif

          !!! Calculate the upscaled KK autoconversion tendency.
          if ( Ncm(k) > Nc_tol ) then

             rrainm_auto(k)  &
             = KK_auto_upscaled_mean( s1(k), s2(k), mu_Nc_n, stdev_s1(k), &
                                      stdev_s2(k), sigma_Nc_n, corr_sNc_1_n, &
                                      corr_sNc_2_n, KK_auto_coef, mixt_frac(k) )

          else  ! N_c = 0.

             rrainm_auto(k) = zero

          endif

          !!! Calculate the upscaled KK accretion tendency.
          if ( rrainm(k) > rr_tol ) then

             rrainm_accr(k)  &
             = KK_accr_upscaled_mean( s1(k), s2(k), mu_rr_n, stdev_s1(k), &
                                      stdev_s2(k), sigma_rr_n, corr_srr_1_n, &
                                      corr_srr_2_n, KK_accr_coef, mixt_frac(k) )

          else  ! r_r = 0.

             rrainm_accr(k) = zero

          endif

!          if ( l_var_covar_src ) then
!
!             ! Calculate the correlation between s and r_t in PDF component 1.
!
!             ! Calculate the correlation between s and r_t in PDF component 2.
!
!             ! Calculate the correlation between s and th_l in PDF component 1.
!
!             ! Calculate the correlation between s and th_l in PDF component 2.
!
!          endif


       else  ! local KK

          !!! Calculate the local KK rain drop mean volume radius.
          if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

             mean_vol_rad(k)  &
             = KK_mvr_local_mean( rrainm(k), Nrm(k), KK_mvr_coef )

          else  ! r_r or N_r = 0.

             mean_vol_rad(k) = zero

          endif

          !!! Calculate the values of the local KK microphysics tendencies.

          !!! Calculate the local KK evaporation tendency.
          if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

             rrainm_cond(k)  &
             = KK_evap_local_mean( s_mellor(k), rrainm(k), Nrm(k), &
                                   KK_evap_coef )

          else  ! r_r or N_r = 0.

             rrainm_cond(k) = zero

          endif

          !!! Calculate the local KK autoconversion tendency.
          if ( Ncm(k) > Nc_tol ) then

             rrainm_auto(k)  &
             = KK_auto_local_mean( s_mellor(k), Ncm(k), KK_auto_coef )

          else  ! N_c = 0.

             rrainm_auto(k) = zero

          endif

          !!! Calculate the local KK accretion tendency.
          if ( rrainm(k) > rr_tol ) then

             rrainm_accr(k)  &
             = KK_accr_local_mean( s_mellor(k), rrainm(k), KK_accr_coef )

          else  ! r_r = 0.

             rrainm_accr(k) = zero

          endif


       endif ! l_upscaled


       !!! KK rain drop concentration microphysics tendencies.

       !!! Calculate the KK N_r evaporation tendency.
       if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

          Nrm_cond(k) = KK_Nrm_evap( rrainm_cond(k), Nrm(k), rrainm(k) )

       else  ! r_r or N_r = 0.

          Nrm_cond(k) = zero

       endif

       !!! Calculate the KK N_r autoconversion tendency.
       Nrm_auto(k) = KK_Nrm_auto( rrainm_auto(k) )


       ! Statistics
       if ( l_stats_samp ) then

          ! Rain drop mean volume radius.
          call stat_update_var_pt( im_vol_rad_rain, k, mean_vol_rad(k), zt )

          ! Explicit contributions to rrainm.
          call stat_update_var_pt( irrainm_cond, k, rrainm_cond(k), zt )

          call stat_update_var_pt( irrainm_auto, k, rrainm_auto(k), zt )

          call stat_update_var_pt( irrainm_accr, k, rrainm_accr(k), zt )

          ! Explicit contributions to Nrm.
          call stat_update_var_pt( iNrm_cond, k, Nrm_cond(k), zt )

          call stat_update_var_pt( iNrm_auto, k, Nrm_auto(k), zt )

       endif  ! l_stats_samp


       !!! Source-adjustment code for rrainm and Nrm.

       rrainm_source = rrainm_auto(k) + rrainm_accr(k)
       Nrm_source = Nrm_auto(k)

       ! The increase of rain due to autoconversion and accretion both draw
       ! their water from the available cloud water.  Over a long time step
       ! these rates may over-deplete cloud water.  In other words, these
       ! processes may draw more cloud water than there is available.  Thus,
       ! the total source rate multiplied by the time step length cannot exceed
       ! the total amount of cloud water available.  If it does, then the rate
       ! must be adjusted.
       total_rc_needed = rrainm_source * real( dt, kind = core_rknd )

       if ( total_rc_needed > rcm(k) .and. l_src_adj_enabled ) then

          ! The maximum allowable rate of the source terms is rcm/dt.
          rrainm_src_max = rcm(k) / real( dt, kind = core_rknd )

          ! The amount of adjustment to the source terms.
          ! This value should always be negative.
          rrainm_src_adj(k) = rrainm_src_max - rrainm_source

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
          rrainm_auto_ratio = rrainm_auto(k) /  &
                              ( rrainm_auto(k) + rrainm_accr(k) )
          Nrm_src_adj(k) = KK_Nrm_auto( rrainm_auto_ratio * rrainm_src_adj(k) )

          ! Change Nrm by Nrm_src_adj.  Nrm_src_adj will always be negative.
          Nrm_source = Nrm_source + Nrm_src_adj(k)

       else

          rrainm_src_adj(k) = zero
          Nrm_src_adj(k)    = zero

       endif

       if ( l_stats_samp ) then

          call stat_update_var_pt( irrainm_src_adj, k, rrainm_src_adj(k), zt )

          call stat_update_var_pt( iNrm_src_adj, k, Nrm_src_adj(k), zt )

       endif ! l_stats_samp


       !!! Calculate overall KK microphysics tendencies.
       rrainm_mc_tndcy(k) = rrainm_cond(k) + rrainm_source
       Nrm_mc_tndcy(k)    = Nrm_cond(k) + Nrm_source

       !!! Explicit contributions to thlm and rtm from the microphysics
       rvm_mc(k)  = -rrainm_cond(k)
       rcm_mc(k)  = -rrainm_source  ! Accretion + Autoconversion
       thlm_mc(k) = ( Lv / ( Cp * exner(k) ) ) * rrainm_mc_tndcy(k)


    enddo  ! Microphysics tendency loop: k = 2, nz-1, 1


    !!! Boundary conditions for microphysics tendencies.

    ! Explicit contributions to rrainm and Nrm from microphysics are not set at
    ! thermodynamic level k = 1 because it is below the model lower boundary.
    rrainm_mc_tndcy(1) = zero
    Nrm_mc_tndcy(1)    = zero

    rrainm_mc_tndcy(nz) = zero
    Nrm_mc_tndcy(nz)    = zero

    ! Boundary conditions
    mean_vol_rad(1)     = zero
    mean_vol_rad(nz) = zero

    rvm_mc(1)     = zero
    rvm_mc(nz) = zero

    rcm_mc(1)     = zero
    rcm_mc(nz) = zero

    thlm_mc(1)     = zero
    thlm_mc(nz) = zero

    !!! Sedimentation velocities
    forall ( k = 1:nz-1 )

       ! Sedimentation velocity of rrainm.
!       Vrr(k) = 0.012_core_rknd * ( micron_per_m * zt2zm(mean_vol_rad,k) ) &
!                - 0.2_core_rknd
       Vrr(k) = 0.012_core_rknd * ( micron_per_m * mean_vol_rad(k) )  &
                - 0.2_core_rknd

       ! Sedimentation velocity is positive upwards.
       Vrr(k) = -max( Vrr(k), zero )

       ! Sedimentation velocity of Nrm.
!       VNr(k) = 0.007_core_rknd * ( micron_per_m * zt2zm(mean_vol_rad,k) ) &
!                - 0.1_core_rknd
       VNr(k) = 0.007_core_rknd * ( micron_per_m * mean_vol_rad(k) )  &
                - 0.1_core_rknd

       ! Sedimentation velocity is positive upwards.
       VNr(k) = -max( VNr(k), zero )

    end forall ! 1..nz-1

    !!! Boundary conditions for sedimentation velocities.

    ! The flux of rain water through the model top is 0.
    ! Vrr and VNr are set to 0 at the highest model level.
    Vrr(nz) = zero
    VNr(nz) = zero


    return

  end subroutine KK_micro_driver

!===============================================================================

end module KK_microphys_module

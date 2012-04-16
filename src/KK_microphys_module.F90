! $Id$
!===============================================================================
module KK_microphys_module

  implicit none

  private

  public :: KK_micro_driver

  private :: KK_upscaled_setup, &
             KK_upscaled_means_driver, &
             KK_upscaled_covar_driver

  contains

  !=============================================================================
  subroutine KK_micro_driver( dt, nz, l_stats_samp, l_local_kk, &
                              l_latin_hypercube, thlm, wm_zt, p_in_Pa, &
                              exner, rho, cloud_frac, pdf_params, w_std_dev, &
                              dzq, rcm, Ncm, s_mellor, rvm, hydromet, &
                              hydromet_mc, hydromet_vel, &
                              rcm_mc, rvm_mc, thlm_mc, &
                              wprtp_mc_tndcy, wpthlp_mc_tndcy, &
                              rtp2_mc_tndcy, thlp2_mc_tndcy, rtpthlp_mc_tndcy, &
                              KK_auto_tndcy, KK_accr_tndcy )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm  ! Procedure(s)

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
        KK_auto_Nc_exp,        &
        C_evap

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
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        zt,                 & ! Variable(s)
        im_vol_rad_rain,    &
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
      T_liq_in_K, & ! Mean liquid water temperature, T_l           [K]
      r_sl,       & ! Liquid water sat. mixing ratio, r_s(T_l,p)   [kg/kg]
      Beta_Tl       ! Parameter Beta, Beta(T_l)                    [1/(kg/kg)]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    real( kind = core_rknd ) :: &
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
      mixt_frac       ! Mixture fraction                                     [-]

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
      rrainm_src_adj, & ! Total adjustment to rrainm source terms  [(kg/kg)/s]
      Nrm_src_adj       ! Total adjustment to Nrm source terms     [{num/kg)/s]

    logical :: &
      l_upscaled,        & ! Flag for using upscaled KK microphysics.
      l_var_covar_src,   & ! Flag for using upscaled covariances.
      l_src_adj_enabled    ! Flag to enable rrainm/Nrm source adjustment

    integer :: &
      k   ! Loop index

    ! Remove compiler warnings
    if ( .false. .and. l_latin_hypercube ) then
      rrainm_src_adj = dzq
      rrainm_src_adj = rvm
      rrainm_src_adj = w_std_dev
      rrainm_src_adj = cloud_frac
    end if

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


    if ( .not. l_local_kk ) then
       l_upscaled = .true.
    else
       l_upscaled = .false.
    endif

    l_var_covar_src = .false.

    l_src_adj_enabled = .true.


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

          call KK_upscaled_setup( rcm(k), rrainm(k), Nrm(k), &
                                  Ncm(k), pdf_params(k), &
                                  mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                  mu_Nc_n, sigma_s_1, sigma_s_2, &
                                  sigma_rr_n, sigma_Nr_n, sigma_Nc_n, &
                                  corr_srr_1_n, corr_srr_2_n, &
                                  corr_sNr_1_n, corr_sNr_2_n, &
                                  corr_sNc_1_n, corr_sNc_2_n, &
                                  corr_rrNr_n, mixt_frac )

          !!! Calculate the values of the upscaled KK microphysics tendencies.
          call KK_upscaled_means_driver( rrainm(k), Nrm(k), Ncm(k), &
                                         mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                         mu_Nc_n, sigma_s_1, sigma_s_2, &
                                         sigma_rr_n, sigma_Nr_n, sigma_Nc_n, &
                                         corr_srr_1_n, corr_srr_2_n, &
                                         corr_sNr_1_n, corr_sNr_2_n, &
                                         corr_sNc_1_n, corr_sNc_2_n, &
                                         corr_rrNr_n,  mixt_frac, &
                                         KK_evap_coef, KK_auto_coef, &
                                         KK_accr_coef, KK_mvr_coef, &
                                         KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                         KK_accr_tndcy(k), KK_mean_vol_rad(k) )
          

          if ( l_var_covar_src ) then

            call KK_upscaled_covar_driver( wm_zt(k), exner(k), rcm(k),  &
                                           rrainm(k), Nrm(k), Ncm(k), &
                                           mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                           mu_Nc_n, sigma_s_1, sigma_s_2, &
                                           sigma_rr_n, sigma_Nr_n, sigma_Nc_n, &
                                           corr_srr_1_n, corr_srr_2_n, &
                                           corr_sNr_1_n, corr_sNr_2_n, &
                                           corr_sNc_1_n, corr_sNc_2_n, &
                                           corr_rrNr_n,  mixt_frac, &
                                           KK_evap_coef, KK_auto_coef, &
                                           KK_accr_coef, KK_evap_tndcy(k), &
                                           KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                           pdf_params(k), k, l_stats_samp, &
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


       endif ! l_upscaled


       !!! KK rain drop concentration microphysics tendencies.

       !!! Calculate the KK N_r evaporation tendency.
       if ( rrainm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

          KK_Nrm_evap_tndcy(k)  &
          = KK_Nrm_evap( KK_evap_tndcy(k), Nrm(k), rrainm(k) )

       else  ! r_r or N_r = 0.

          KK_Nrm_evap_tndcy(k) = zero

       endif

       !!! Calculate the KK N_r autoconversion tendency.
       KK_Nrm_auto_tndcy(k) = KK_Nrm_auto( KK_auto_tndcy(k) )


       ! Statistics
       if ( l_stats_samp ) then

          ! Rain drop mean volume radius.
          call stat_update_var_pt( im_vol_rad_rain, k, KK_mean_vol_rad(k), zt )

          ! Explicit contributions to rrainm.
          call stat_update_var_pt( irrainm_cond, k, KK_evap_tndcy(k), zt )

          call stat_update_var_pt( irrainm_auto, k, KK_auto_tndcy(k), zt )

          call stat_update_var_pt( irrainm_accr, k, KK_accr_tndcy(k), zt )

          ! Explicit contributions to Nrm.
          call stat_update_var_pt( iNrm_cond, k, KK_Nrm_evap_tndcy(k), zt )

          call stat_update_var_pt( iNrm_auto, k, KK_Nrm_auto_tndcy(k), zt )

       endif  ! l_stats_samp


       !!! Source-adjustment code for rrainm and Nrm.

       rrainm_source = KK_auto_tndcy(k) + KK_accr_tndcy(k)
       Nrm_source = KK_Nrm_auto_tndcy(k)

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
          rrainm_auto_ratio = KK_auto_tndcy(k) /  &
                              ( KK_auto_tndcy(k) + KK_accr_tndcy(k) )
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
       rrainm_mc_tndcy(k) = KK_evap_tndcy(k) + rrainm_source
       Nrm_mc_tndcy(k)    = KK_Nrm_evap_tndcy(k) + Nrm_source

       !!! Explicit contributions to thlm and rtm from the microphysics
       rvm_mc(k)  = -KK_evap_tndcy(k)
       rcm_mc(k)  = -rrainm_source  ! Accretion + Autoconversion
       thlm_mc(k) = ( Lv / ( Cp * exner(k) ) ) * rrainm_mc_tndcy(k)


    enddo  ! Microphysics tendency loop: k = 2, nz-1, 1


    if ( l_upscaled .and. l_var_covar_src ) then

       ! Output microphysics tendency terms for
       ! model variances and covariances on momentum levels.
       wprtp_mc_tndcy   = zt2zm( wprtp_mc_tndcy_zt )
       wpthlp_mc_tndcy  = zt2zm( wpthlp_mc_tndcy_zt )
       rtp2_mc_tndcy    = zt2zm( rtp2_mc_tndcy_zt )
       thlp2_mc_tndcy   = zt2zm( thlp2_mc_tndcy_zt )
       rtpthlp_mc_tndcy = zt2zm( rtpthlp_mc_tndcy_zt )

    else

       ! Set values to 0.
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


    return

  end subroutine KK_micro_driver

  !=============================================================================
  subroutine KK_upscaled_setup( rcm, rrainm, Nrm, & 
                                Ncm, pdf_params, &
                                mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                mu_Nc_n, sigma_s_1, sigma_s_2, &
                                sigma_rr_n, sigma_Nr_n, sigma_Nc_n, &
                                corr_srr_1_n, corr_srr_2_n, &
                                corr_sNr_1_n, corr_sNr_2_n, &
                                corr_sNc_1_n, corr_sNc_2_n, &
                                corr_rrNr_n, mixt_frac )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol, & ! Constant(s)
        rr_tol, & 
        Nr_tol, & 
        Nc_tol

    use KK_utilities, only: &
        mean_L2N,   & ! Procedure(s)
        stdev_L2N,  &
        corr_NL2NN, &
        corr_LL2NN

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use parameters_microphys, only: &
        rrp2_on_rrainm2_cloud, & ! Variable(s)
        rrp2_on_rrainm2_below, &
        Nrp2_on_Nrm2_cloud,    &
        Nrp2_on_Nrm2_below,    &
        Ncp2_on_Ncm2_cloud,    &
        Ncp2_on_Ncm2_below

    use KK_fixed_correlations, only: &
        corr_srr_NL_cloud,  & ! Variable(s)
        corr_sNr_NL_cloud,  &
        corr_sNc_NL_cloud,  &
        corr_srr_NL_below,  &
        corr_sNr_NL_below,  &
        corr_sNc_NL_below,  &
        corr_rrNr_LL_cloud, &
        corr_rrNr_LL_below

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rcm,    & ! Mean cloud water mixing ratio       [kg/kg]
      rrainm, & ! Mean rain water mixing ratio        [kg/kg]
      Nrm,    & ! Mean rain drop concentration        [num/kg]
      Ncm       ! Mean cloud droplet concentration    [num/kg]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                        [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
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
      mixt_frac       ! Mixture fraction                                     [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      rrp2_on_rrainm2, & ! Ratio of < r_r >^2 to < r_r'^2 >                 [-]
      Nrp2_on_Nrm2,    & ! Ratio of < N_r >^2 to < N_r'^2 >                 [-]
      Ncp2_on_Ncm2,    & ! Ratio of < N_c >^2 to < N_c'^2 >                 [-]
      corr_srr_1,      & ! Correlation between s and rr (1st PDF component) [-] 
      corr_srr_2,      & ! Correlation between s and rr (2nd PDF component) [-]
      corr_sNr_1,      & ! Correlation between s and Nr (1st PDF component) [-]
      corr_sNr_2,      & ! Correlation between s and Nr (2nd PDF component) [-]
      corr_sNc_1,      & ! Correlation between s and Nc (1st PDF component) [-]
      corr_sNc_2,      & ! Correlation between s and Nc (2nd PDF component) [-]
      corr_rrNr          ! Correlation between rr and Nr (both components)  [-]


    ! Enter the PDF parameters.
    mu_s_1    = pdf_params%s1
    mu_s_2    = pdf_params%s2
    sigma_s_1 = pdf_params%stdev_s1
    sigma_s_2 = pdf_params%stdev_s2
    mixt_frac = pdf_params%mixt_frac
    
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
    if ( rcm > rc_tol ) then
       rrp2_on_rrainm2 = rrp2_on_rrainm2_cloud
       Nrp2_on_Nrm2    = Nrp2_on_Nrm2_cloud
       Ncp2_on_Ncm2    = Ncp2_on_Ncm2_cloud
       corr_rrNr       = corr_rrNr_LL_cloud
       corr_srr_1      = corr_srr_NL_cloud
       corr_srr_2      = corr_srr_NL_cloud
       corr_sNr_1      = corr_sNr_NL_cloud
       corr_sNr_2      = corr_sNr_NL_cloud
       corr_sNc_1      = corr_sNc_NL_cloud
       corr_sNc_2      = corr_sNc_NL_cloud
    else
       rrp2_on_rrainm2 = rrp2_on_rrainm2_below
       Nrp2_on_Nrm2    = Nrp2_on_Nrm2_below
       Ncp2_on_Ncm2    = Ncp2_on_Ncm2_below
       corr_rrNr       = corr_rrNr_LL_below
       corr_srr_1      = corr_srr_NL_below
       corr_srr_2      = corr_srr_NL_below
       corr_sNr_1      = corr_sNr_NL_below
       corr_sNr_2      = corr_sNr_NL_below
       corr_sNc_1      = corr_sNc_NL_below
       corr_sNc_2      = corr_sNc_NL_below
    endif

    !!! Calculate the normalized mean of variables that have an assumed (single)
    !!! lognormal distribution, given the mean and variance of those variables.

    ! Normalized mean of rain water mixing ratio.
    if ( rrainm > rr_tol ) then
       mu_rr_n = mean_L2N( rrainm, rrp2_on_rrainm2 * rrainm**2 )
    endif

    ! Normalized mean of rain drop concentration.
    if ( Nrm > Nr_tol ) then
       mu_Nr_n = mean_L2N( Nrm, Nrp2_on_Nrm2 * Nrm**2 )
    endif

    ! Normalized mean of cloud droplet concentration.
    if ( Ncm > Nc_tol ) then
       mu_Nc_n = mean_L2N( Ncm, Ncp2_on_Ncm2 * Ncm**2 )
    endif

    !!! Calculate the normalized standard deviation of variables that have
    !!! an assumed (single) lognormal distribution, given the mean and
    !!! variance of those variables.

    ! Normalized standard deviation of rain water mixing ratio.
    if ( rrainm > rr_tol ) then
       sigma_rr_n = stdev_L2N( rrainm, rrp2_on_rrainm2 * rrainm**2 )
    endif

    ! Normalized standard deviation of rain drop concentration.
    if ( Nrm > Nr_tol ) then
       sigma_Nr_n = stdev_L2N( Nrm, Nrp2_on_Nrm2 * Nrm**2 )
    endif

    ! Normalized standard deviation of cloud droplet concentration.
    if ( Ncm > Nc_tol ) then
       sigma_Nc_n = stdev_L2N( Ncm, Ncp2_on_Ncm2 * Ncm**2 )
    endif

    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! (single) lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    if ( rrainm > rr_tol ) then

       ! Normalize the correlation between s and r_r in PDF component 1.
       corr_srr_1_n = corr_NL2NN( corr_srr_1, sigma_rr_n )

       ! Normalize the correlation between s and r_r in PDF component 2.
       corr_srr_2_n = corr_NL2NN( corr_srr_2, sigma_rr_n )

    endif

    if ( Nrm > Nr_tol ) then

       ! Normalize the correlation between s and N_r in PDF component 1.
       corr_sNr_1_n = corr_NL2NN( corr_sNr_1, sigma_Nr_n )

       ! Normalize the correlation between s and N_r in PDF component 2.
       corr_sNr_2_n = corr_NL2NN( corr_sNr_2, sigma_Nr_n )

    endif

    if ( Ncm > Nc_tol ) then

       ! Normalize the correlation between s and N_c in PDF component 1.
       corr_sNc_1_n = corr_NL2NN( corr_sNc_1, sigma_Nc_n )

       ! Normalize the correlation between s and N_c in PDF component 2.
       corr_sNc_2_n = corr_NL2NN( corr_sNc_2, sigma_Nc_n )

    endif

    !!! Calculate the normalized correlation between two variables that both
    !!! have an assumed lognormal distribution, given their correlation and both
    !!! of their normalized standard deviations.

    ! Normalize the correlation between rr and Nr (this is the same for
    ! both PDF components).
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then
       corr_rrNr_n = corr_LL2NN( corr_rrNr, sigma_rr_n, sigma_Nr_n )
    endif


    return

  end subroutine KK_upscaled_setup

  !=============================================================================
  subroutine KK_upscaled_means_driver( rrainm, Nrm, Ncm, &
                                       mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                       mu_Nc_n, sigma_s_1, sigma_s_2, &
                                       sigma_rr_n, sigma_Nr_n, sigma_Nc_n, &
                                       corr_srr_1_n, corr_srr_2_n, &
                                       corr_sNr_1_n, corr_sNr_2_n, &
                                       corr_sNc_1_n, corr_sNc_2_n, &
                                       corr_rrNr_n,  mixt_frac, &
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
      mixt_frac       ! Mixture fraction                                     [-]

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
       = KK_evap_upscaled_mean( mu_s_1, mu_s_2, mu_rr_n, mu_Nr_n, &
                                sigma_s_1, sigma_s_2, sigma_rr_n, &
                                sigma_Nr_n, corr_srr_1_n, corr_srr_2_n, &
                                corr_sNr_1_n, corr_sNr_2_n, corr_rrNr_n, &
                                KK_evap_coef, mixt_frac )

    else  ! r_r or N_r = 0.

       KK_evap_tndcy = zero

    endif

    !!! Calculate the upscaled KK autoconversion tendency.
    if ( Ncm > Nc_tol ) then

       KK_auto_tndcy  &
       = KK_auto_upscaled_mean( mu_s_1, mu_s_2, mu_Nc_n, sigma_s_1, &
                                sigma_s_2, sigma_Nc_n, corr_sNc_1_n, &
                                corr_sNc_2_n, KK_auto_coef, mixt_frac )

    else  ! N_c = 0.

       KK_auto_tndcy = zero

    endif

    !!! Calculate the upscaled KK accretion tendency.
    if ( rrainm > rr_tol ) then

       KK_accr_tndcy  &
       = KK_accr_upscaled_mean( mu_s_1, mu_s_2, mu_rr_n, sigma_s_1, &
                                sigma_s_2, sigma_rr_n, corr_srr_1_n, &
                                corr_srr_2_n, KK_accr_coef, mixt_frac )

    else  ! r_r = 0.

       KK_accr_tndcy = zero

    endif

    !!! Calculate the upscaled KK rain drop mean volume radius.
    if ( rrainm > rr_tol .and. Nrm > Nr_tol ) then

       KK_mean_vol_rad &
       = KK_mvr_upscaled_mean( mu_rr_n, mu_Nr_n, sigma_rr_n, &
                               sigma_Nr_n, corr_rrNr_n, KK_mvr_coef )

    else  ! r_r or N_r = 0.

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
                                       KK_evap_coef, KK_auto_coef, &
                                       KK_accr_coef, KK_evap_tndcy, &
                                       KK_auto_tndcy, KK_accr_tndcy, &
                                       pdf_params, level, l_stats_samp, &
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
        corr_trr_NL_cloud, & ! Variable(s)
        corr_tNr_NL_cloud, &
        corr_tNc_NL_cloud, &
        corr_trr_NL_below, &
        corr_tNr_NL_below, &
        corr_tNc_NL_below

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
      mixt_frac       ! Mixture fraction                                     [-]

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
    corr_ts_1 = pdf_params%corr_st_1
    corr_ts_2 = pdf_params%corr_st_2
    crt1      = pdf_params%crt1
    crt2      = pdf_params%crt2
    cthl1     = pdf_params%cthl1
    cthl2     = pdf_params%cthl2

    ! Set all the correlations between vertical velocity (w) and hydrometeor
    ! variables to 0 for now.
    corr_wrr_1 = zero
    corr_wrr_2 = zero
    corr_wNr_1 = zero
    corr_wNr_2 = zero
    corr_wNc_1 = zero
    corr_wNc_2 = zero

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

    endif

    if ( Nrm > Nr_tol ) then

       ! Normalize the correlation between w and N_r in PDF component 1.
       corr_wNr_1_n = corr_NL2NN( corr_wNr_1, sigma_Nr_n )

       ! Normalize the correlation between w and N_r in PDF component 2.
       corr_wNr_2_n = corr_NL2NN( corr_wNr_2, sigma_Nr_n )

    endif

    if ( Ncm > Nc_tol ) then

       ! Normalize the correlation between s and N_c in PDF component 1.
       corr_wNc_1_n = corr_NL2NN( corr_wNc_1, sigma_Nc_n )

       ! Normalize the correlation between s and N_c in PDF component 2.
       corr_wNc_2_n = corr_NL2NN( corr_wNc_2, sigma_Nc_n )

    endif

    
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
       corr_trr_1 = corr_trr_NL_cloud
       corr_trr_2 = corr_trr_NL_cloud
       corr_tNr_1 = corr_tNr_NL_cloud
       corr_tNr_2 = corr_tNr_NL_cloud
       corr_tNc_1 = corr_tNc_NL_cloud
       corr_tNc_2 = corr_tNc_NL_cloud
    else
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

       ! Normalize the correlation between t and r_r in PDF component 1.
       corr_trr_1_n = corr_NL2NN( corr_trr_1, sigma_rr_n )

       ! Normalize the correlation between t and r_r in PDF component 2.
       corr_trr_2_n = corr_NL2NN( corr_trr_2, sigma_rr_n )

    endif

    if ( Nrm > Nr_tol ) then

       ! Normalize the correlation between t and N_r in PDF component 1.
       corr_tNr_1_n = corr_NL2NN( corr_tNr_1, sigma_Nr_n )

       ! Normalize the correlation between t and N_r in PDF component 2.
       corr_tNr_2_n = corr_NL2NN( corr_tNr_2, sigma_Nr_n )

    endif

    if ( Ncm > Nc_tol ) then

       ! Normalize the correlation between t and N_c in PDF component 1.
       corr_tNc_1_n = corr_NL2NN( corr_tNc_1, sigma_Nc_n )

       ! Normalize the correlation between t and N_c in PDF component 2.
       corr_tNc_2_n = corr_NL2NN( corr_tNc_2, sigma_Nc_n )

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
                          KK_evap_tndcy, KK_evap_coef, w_tol, mixt_frac )

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
                           KK_evap_coef, t_tol, crt1, crt2, mixt_frac )

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
                            KK_evap_coef, t_tol, cthl1, cthl2, mixt_frac )

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
                          KK_auto_coef, w_tol, mixt_frac )

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
                           t_tol, crt1, crt2, mixt_frac )

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
                            t_tol, cthl1, cthl2, mixt_frac )

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
                          KK_accr_coef, w_tol, mixt_frac )

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
                           t_tol, crt1, crt2, mixt_frac )

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
                            t_tol, cthl1, cthl2, mixt_frac )

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

end module KK_microphys_module

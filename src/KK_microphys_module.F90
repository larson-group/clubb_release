! $Id$
!===============================================================================
module KK_microphys_module

  implicit none

  private

  public :: KK_local_micro_driver, &
            KK_upscaled_micro_driver

  private :: KK_micro_init, &
             KK_tendency_coefs, &
             KK_microphys_adjust, &
             KK_upscaled_stats, &
             KK_stats_output, &
             KK_sedimentation

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
                                    rtp2_mc_tndcy, thlp2_mc_tndcy, &
                                    KK_auto_tndcy, KK_accr_tndcy, KK_evap_tndcy )

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
      rtp2_mc_tndcy,  & ! Microphysics tendency for <rt'^2>        [(kg/kg)^2/s]
      thlp2_mc_tndcy, & ! Microphysics tendency for <thl'^2>       [K^2/s]
      KK_auto_tndcy,  & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy,  & ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]
      KK_evap_tndcy     ! Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(:), pointer ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
      Vrr,             & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,             & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) :: &
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


    ! Microphysics tendency terms for model variances are set to 0.
    rtp2_mc_tndcy  = zero
    thlp2_mc_tndcy = zero

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
  subroutine KK_upscaled_micro_driver( dt, nz, l_stats_samp, thlm, wm_zt,       & ! Intent(in)
                                       p_in_Pa, exner, rho, cloud_frac,         & ! Intent(in)
                                       pdf_params, w_std_dev, rcm, Ncnm,        & ! Intent(in)
                                       s_mellor, Nc0_in_cloud,                  & ! Intent(in)
                                       hydromet, wphydrometp,                   & ! Intent(in)
                                       d_variables, corr_array_1, corr_array_2, & ! Intent(in)
                                       mu_x_1, mu_x_2, sigma_x_1, sigma_x_2,    & ! Intent(in)
                                       hydromet_pdf_params,                     & ! Intent(in)
                                       hydromet_mc, hydromet_vel,               & ! Intent(inout)
                                       rcm_mc, rvm_mc, thlm_mc,                 & ! Intent(out)
                                       hydromet_vel_covar_zt_impc,              & ! Intent(out)
                                       hydromet_vel_covar_zt_expc,              & ! Intent(out)
                                       wprtp_mc_tndcy, wpthlp_mc_tndcy,         & ! Intent(out)
                                       rtp2_mc_tndcy, thlp2_mc_tndcy,           & ! Intent(out)
                                       rtpthlp_mc_tndcy )                         ! Intent(out)
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        zero    ! Constant(s)

    use parameters_microphys, only: &
        l_var_covar_src,     & ! Flag for using variance/covariance src terms
        l_const_Nc_in_cloud    ! Flag to use a const. value of N_c within cloud

    use setup_clubb_pdf_params, only: &
        setup_pdf_parameters, & ! Procedure(s)
        unpack_pdf_params, &
        normalize_pdf_params, &
        pdf_param_hm_stats, &
        pdf_param_log_hm_stats

    use KK_upscaled_means, only: &
        KK_upscaled_means_driver ! Procedure(s)

    use KK_Nrm_tendencies, only: &
        KK_Nrm_evap_upscaled_mean, & ! Procedure(s)
        KK_Nrm_auto_mean

    use KK_upscaled_turbulent_sed, only: &
        KK_sed_vel_covars  ! Procedure(s)

    use KK_upscaled_covariances, only: &
        KK_upscaled_covar_driver    ! Procedure(s)

    use KK_utilities, only: &
        get_cloud_top_level    ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

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
        iVrrprrp_expcalc, &
        iVNrpNrp_expcalc, &
        zm

    implicit none

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      nz,         & ! Number of model vertical grid levels
      d_variables   ! Number of variables in the correlation arrays

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
      Ncnm,       & ! Mean cloud nuclei concentration < N_cn >        [num/kg]
      s_mellor      ! Mean extended liquid water mixing ratio         [kg/kg]

    type(pdf_parameter), dimension(nz), target, intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      w_std_dev    ! Standard deviation of w (for LH interface)          [m/s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Nc0_in_cloud    ! Constant in-cloud value of cloud droplet conc.  [num/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet,    & ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)      [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), intent(in) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    real( kind = core_rknd ), dimension(d_variables, nz), intent(in) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

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
      rtpthlp_mc_tndcy    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(:), pointer ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r > [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r > [num/kg]
      Vrr,             & ! Mean sedimentation velocity of r_r    [m/s]
      VNr,             & ! Mean sedimentation velocity of N_r    [m/s]
      rrainm_mc_tndcy, & ! Mean (dr_r/dt) due to microphysics    [(kg/kg)/s]
      Nrm_mc_tndcy       ! Mean (dN_r/dt) due to microphysics    [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) :: &
      KK_evap_tndcy,   & ! Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
      KK_auto_tndcy,   & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy,   & ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]
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
      mixt_frac   ! Mixture fraction                                [-]

    real( kind = core_rknd ) :: &
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n    ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]

    real( kind = core_rknd ) :: &
      corr_ws_1,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_ws_2,     & ! Correlation between w and s (2nd PDF component)     [-]
      !corr_wrr_1,    & ! Correlation between w and rr (1st PDF component) ip [-]
      !corr_wrr_2,    & ! Correlation between w and rr (2nd PDF component) ip [-]
      !corr_wNr_1,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      !corr_wNr_2,    & ! Correlation between w and Nr (2nd PDF component) ip [-]
      !corr_wNcn_1,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      !corr_wNcn_2,   & ! Correlation between w and Ncn (2nd PDF component)   [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2!,     & ! Correlation between s and t (2nd PDF component)     [-]
      !corr_srr_1,    & ! Correlation between s and rr (1st PDF component) ip [-]
      !corr_srr_2,    & ! Correlation between s and rr (2nd PDF component) ip [-]
      !corr_sNr_1,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      !corr_sNr_2,    & ! Correlation between s and Nr (2nd PDF component) ip [-]
      !corr_sNcn_1,   & ! Correlation between s and Ncn (1st PDF component)   [-]
      !corr_sNcn_2,   & ! Correlation between s and Ncn (2nd PDF component)   [-]
      !corr_trr_1,    & ! Correlation between t and rr (1st PDF component) ip [-]
      !corr_trr_2,    & ! Correlation between t and rr (2nd PDF component) ip [-]
      !corr_tNr_1,    & ! Correlation between t and Nr (1st PDF component) ip [-]
      !corr_tNr_2,    & ! Correlation between t and Nr (2nd PDF component) ip [-]
      !corr_tNcn_1,   & ! Correlation between t and Ncn (1st PDF component)   [-]
      !corr_tNcn_2,   & ! Correlation between t and Ncn (2nd PDF component)   [-]
      !corr_rrNr_1,   & ! Correlation between rr & Nr (1st PDF component) ip  [-]
      !corr_rrNr_2      ! Correlation between rr & Nr (2nd PDF component) ip  [-]

    real( kind = core_rknd ) :: &
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

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

    logical :: &
      l_src_adj_enabled,  & ! Flag to enable rrainm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrainm/Nrm evaporation adjustment

    integer :: &
      cloud_top_level, & ! Vertical level index of cloud top 
      k                  ! Loop index

    ! ---- Begin Code ----

    !!! Initialize microphysics fields.
    call KK_micro_init( nz, hydromet, hydromet_mc, hydromet_vel, &
                        rrainm, Nrm, Vrr, VNr, &
                        rrainm_mc_tndcy, Nrm_mc_tndcy, &
                        KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                        KK_mean_vol_rad, KK_Nrm_evap_tndcy, &
                        KK_Nrm_auto_tndcy, &
                        l_src_adj_enabled, l_evap_adj_enabled )

    ! The implicit and explicit components used to calculate the covariances of
    ! hydrometeor sedimentation velocities and their associated hydrometeors
    ! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp_zt_impc => hydromet_vel_covar_zt_impc(:,iirrainm)
    Vrrprrp_zt_expc => hydromet_vel_covar_zt_expc(:,iirrainm)
    VNrpNrp_zt_impc => hydromet_vel_covar_zt_impc(:,iiNrm)
    VNrpNrp_zt_expc => hydromet_vel_covar_zt_expc(:,iiNrm)

    ! Setup mixture fraction.
    mixt_frac   = pdf_params%mixt_frac    

    !!! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

       !!! Calculate the coefficients for the KK microphysics tendencies.
       call KK_tendency_coefs( thlm(k), exner(k), p_in_Pa(k), rho(k), &
                               KK_evap_coef, KK_auto_coef, &
                               KK_accr_coef, KK_mvr_coef )


       !!! Unpack the PDF parameters.
       call unpack_pdf_params( d_variables, corr_array_1(:,:,k), corr_array_2(:,:,k), & ! In
                               mu_x_1(:,k), mu_x_2(:,k), sigma_x_1(:,k), sigma_x_2(:,k), &
                               hydromet_pdf_params(k), &
                               mu_w_1, mu_w_2, mu_s_1, mu_s_2, & ! Out
                               mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &
                               mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                               mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                               mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                               sigma_w_1, sigma_w_2, sigma_s_1, &
                               sigma_s_2, sigma_t_1, sigma_t_2, &
                               sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                               sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                               sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                               sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                               corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, &
                               corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &
                               corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &
                               corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                               corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                               corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &
                               corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &
                               corr_rrNr_1_n, corr_rrNr_2_n, &
                               rr1, rr2, Nr1, Nr2, &
                               precip_frac, precip_frac_1, precip_frac_2 )

       !!! Calculate the values of the upscaled KK microphysics tendencies.
       call KK_upscaled_means_driver( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, &
                                      mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                      mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                      mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                      sigma_s_1, sigma_s_2, &
                                      sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                      sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                      sigma_rr_1_n, sigma_rr_2_n, &
                                      sigma_Nr_1_n, sigma_Nr_2_n, &
                                      sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                      corr_srr_1_n, corr_srr_2_n, &
                                      corr_sNr_1_n, corr_sNr_2_n, &
                                      corr_sNcn_1_n, corr_sNcn_2_n, &
                                      corr_rrNr_1_n, corr_rrNr_2_n, &
                                      mixt_frac(k), precip_frac_1, &
                                      precip_frac_2, Nc0_in_cloud(k), &
                                      l_const_Nc_in_cloud, &
                                      KK_evap_coef, KK_auto_coef, &
                                      KK_accr_coef, KK_mvr_coef, &
                                      KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                      KK_accr_tndcy(k), KK_mean_vol_rad(k) )

       call KK_sed_vel_covars( rrainm(k), rr1, rr2, Nrm(k), &
                               Nr1, Nr2, KK_mean_vol_rad(k), &
                               mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                               mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                               sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                               sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                               sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                               KK_mvr_coef, mixt_frac(k), precip_frac_1, &
                               precip_frac_2, k, l_stats_samp, &
                               Vrrprrp_zt_impc(k), Vrrprrp_zt_expc(k), &
                               VNrpNrp_zt_impc(k), VNrpNrp_zt_expc(k) )


       if ( l_var_covar_src ) then

          call KK_upscaled_covar_driver( wm_zt(k), exner(k), &
                                         rrainm(k), Nrm(k), Ncnm(k), &
                                         mu_w_1, mu_w_2, mu_s_1, mu_s_2, &
                                         mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &
                                         mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                         mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                         mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                         sigma_w_1, sigma_w_2, sigma_s_1, &
                                         sigma_s_2, sigma_t_1, sigma_t_2, &
                                         sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                         sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                         sigma_rr_1_n, sigma_rr_2_n, &
                                         sigma_Nr_1_n, sigma_Nr_2_n, &
                                         sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                         corr_ws_1, corr_ws_2, corr_wrr_1_n, &
                                         corr_wrr_2_n, corr_wNr_1_n, &
                                         corr_wNr_2_n, corr_wNcn_1_n, &
                                         corr_wNcn_2_n, corr_st_1, corr_st_2, &
                                         corr_srr_1_n, corr_srr_2_n, &
                                         corr_sNr_1_n, corr_sNr_2_n, &
                                         corr_sNcn_1_n, corr_sNcn_2_n, &
                                         corr_trr_1_n, corr_trr_2_n, &
                                         corr_tNr_1_n, corr_tNr_2_n, &
                                         corr_tNcn_1_n, corr_tNcn_2_n, &
                                         corr_rrNr_1_n, corr_rrNr_2_n, &
                                         mixt_frac(k), precip_frac_1, &
                                         precip_frac_2, Nc0_in_cloud(k), &
                                         KK_evap_coef, KK_auto_coef, &
                                         KK_accr_coef, KK_evap_tndcy(k), &
                                         KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                         pdf_params(k), k, &
                                         l_const_Nc_in_cloud, l_stats_samp, &
                                         wprtp_mc_tndcy_zt(k), &
                                         wpthlp_mc_tndcy_zt(k), &
                                         rtp2_mc_tndcy_zt(k), &
                                         thlp2_mc_tndcy_zt(k), &
                                         rtpthlp_mc_tndcy_zt(k) )

       endif

       !!! KK rain drop concentration microphysics tendencies.

       !!! Calculate the KK N_r evaporation tendency.
       KK_Nrm_evap_tndcy(k) &
       = KK_Nrm_evap_upscaled_mean( mu_s_1, mu_s_2, mu_rr_1, mu_rr_2, &
                                    mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                                    mu_Nr_1_n, mu_Nr_2_n, sigma_s_1, &
                                    sigma_s_2, sigma_rr_1, sigma_rr_2, &
                                    sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                    sigma_rr_2_n, sigma_Nr_1_n, &
                                    sigma_Nr_2_n, corr_srr_1_n, &
                                    corr_srr_2_n, corr_sNr_1_n, &
                                    corr_sNr_2_n, corr_rrNr_1_n, &
                                    corr_rrNr_2_n, KK_evap_coef, mixt_frac(k), &
                                    precip_frac_1, precip_frac_2, dt )


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
                               mu_Ncn_1, mu_Ncn_2, mu_rr_1_n, mu_rr_2_n, &
                               mu_Nr_1_n, mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                               sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                               sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                               sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                               sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                               corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                               corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                               corr_rrNr_1_n, corr_rrNr_2_n, mixt_frac(k), &
                               precip_frac_1, precip_frac_2, &
                               KK_mvr_coef, KK_mean_vol_rad(k), rrainm(k), &
                               Nrm(k), k, l_stats_samp )

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
                                mu_Ncn_1, mu_Ncn_2, mu_rr_1_n, mu_rr_2_n, &
                                mu_Nr_1_n, mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                                corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                                corr_rrNr_1_n, corr_rrNr_2_n, mixt_frac, &
                                precip_frac_1, precip_frac_2, &
                                KK_mvr_coef, KK_mean_vol_rad, rrainm, &
                                Nrm, level, l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        rr_tol, & ! Constant(s)
        Nr_tol, &
        zero

    use KK_upscaled_variances, only: &
        variance_KK_mvr    ! Procedure(s)

    use KK_utilities, only: &
        calc_xp2    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        iKK_mvr_variance_zt, &
        irrp2_zt,            &
        iNrp2_zt,            &
        zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n, & ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n, & ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), intent(in) :: &
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

       ! Overall variance of r_r.
       if ( irrp2_zt > 0 ) then

          ! Calculate the overall variance of r_r, <r_r'^2>.
          if ( rrainm > rr_tol ) then
             rrp2 &
             = calc_xp2( mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_rr_1, &
                         sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, mixt_frac, &
                         precip_frac_1, precip_frac_2, rrainm, rr_tol )
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
             = calc_xp2( mu_Nr_1, mu_Nr_2, mu_Nr_1_n, mu_Nr_2_n, sigma_Nr_1, &
                         sigma_Nr_2, sigma_Nr_1_n, sigma_Nr_2_n, mixt_frac, &
                         precip_frac_1, precip_frac_2, Nrm, Nr_tol )
          else ! N_r = 0.
             Nrp2 = zero
          endif

          call stat_update_var_pt( iNrp2_zt, level, Nrp2, zt )

       endif ! iNrp2_zt > 0

       ! Variance of KK rain drop mean volume radius.
       if ( iKK_mvr_variance_zt > 0 ) then

          ! Calculate the variance of KK rain drop mean volume radius,
          ! < R_vr'^2 >.
          KK_mvr_variance &
          = variance_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                             mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                             sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                             sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                             sigma_Nr_2_n, corr_rrNr_1_n, corr_rrNr_2_n, &
                             KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                             precip_frac_1, precip_frac_2 )

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

!===============================================================================

end module KK_microphys_module

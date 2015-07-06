! $Id$
!===============================================================================
module KK_microphys_module

  use constants_clubb, only: &
      zero  ! Constant(s)

  use clubb_precision, only: &
      core_rknd  ! Variable(s)

  implicit none

  private

  public :: KK_local_microphys,          &
            KK_upscaled_microphys,       &
            KK_microphys_adjust,         &
            KK_microphys_adj_terms_type

  private :: KK_microphys_init,    &
             KK_tendency_coefs,    &
             KK_upscaled_stats,    &
             KK_stats_output,      &
             KK_sedimentation,     &
             KK_microphys_output,  &
             unpack_pdf_params_KK

  ! This derived type stores information about adjustment to microphysics terms
  ! in KK microphysics.
  type KK_microphys_adj_terms_type

    real( kind = core_rknd ) :: &
      rrm_src_adj  = zero, & ! Total adj. to rrm source terms    [(kg/kg)/s]
      rrm_cond_adj = zero, & ! Total adj. to rrm evap terms      [(kg/kg)/s]
      Nrm_src_adj  = zero, & ! Total adj. to Nrm source terms    [(num/kg)/s]
      Nrm_cond_adj = zero    ! Total adj. to Nrm evap terms      [(num/kg)/s]

  end type KK_microphys_adj_terms_type

  contains

  !=============================================================================
  subroutine KK_local_microphys( dt, nz, l_latin_hypercube,             & ! In
                                 thlm, wm_zt, p_in_Pa, exner, rho,      & ! In
                                 cloud_frac, w_std_dev, dzq, rcm,       & ! In
                                 Ncm, chi, rvm, hydromet,               & ! In
                                 hydromet_mc, hydromet_vel,             & ! Out
                                 Ncm_mc, rcm_mc, rvm_mc, thlm_mc,       & ! Out
                                 microphys_stats_zt,                    & ! Out
                                 microphys_stats_sfc )                    ! Out

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero,   &
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

    use parameters_microphys, only: &
        l_silhs_KK_convergence_adj_mean ! Variable

    use microphys_stats_vars_module, only: &
        microphys_stats_vars_type, &
        microphys_stats_alloc, &
        microphys_put_var

    use stats_variables, only: &
        im_vol_rad_rain, & ! Variables
        iNrm_auto, &
        iNrm_cond, &
        iNrm_cond_adj, &
        iNrm_src_adj, &
        irrm_accr, &
        irrm_auto, &
        irrm_cond, &
        irrm_cond_adj, &
        irrm_src_adj, &
        irrm_mc_nonadj

    implicit none

    ! Local Constants
    integer, parameter :: &
      num_stats_zt = 20,  &       ! Overestimate of the number of statistics
                                  ! variables sampled in this subroutine

      num_stats_sfc = 0           ! No sfc variables sampled in this routine

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    logical, intent(in) :: &
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
      chi           ! Mean extended liquid water mixing ratio         [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      w_std_dev, & ! Standard deviation of w (for LH interface)          [m/s]
      dzq,       & ! Thickness between thermo. levels (for LH interface) [m]
      rvm          ! Mean water vapor mixing ratio (for LH interface)    [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Hydrometeor species                      [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      Ncm_mc,  & ! Time tendency of cloud droplet concentration  [num/kg/s]
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    type(microphys_stats_vars_type), intent(out) :: &
      microphys_stats_zt, & ! Variables output for statistical sampling (zt grid)
      microphys_stats_sfc   ! Variables output for statistical sampling (sfc grid)

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      KK_auto_tndcy,     & ! Mean KK (dr_r/dt) due to autoconversion   [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK (dr_r/dt) due to accretion        [(kg/kg)/s]
      KK_evap_tndcy,     & ! Mean KK (dr_r/dt) due to evaporation      [(kg/kg)/s]
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation      [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.        [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) ::  &
      rrm,    & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,    & ! Mean rain drop concentration, < N_r >    [num/kg]
      Vrr,    & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,    & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrm_mc, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc    ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(nz) :: &
      KK_mean_vol_rad    ! Mean KK rain drop mean volume radius     [m]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    type(KK_microphys_adj_terms_type), dimension(nz) :: &
      adj_terms    ! Adjustments to microphysics terms

    logical :: &
      l_src_adj_enabled,   & ! Flag to enable rrm/Nrm source adjustment
      l_evap_adj_enabled,  & ! Flag to enable rrm/Nrm evaporation adjustment
      l_clip_positive_sed    ! Flag to enable Vrr/VNr positive clipping

    integer :: &
      cloud_top_level, & ! Vertical level index of cloud top 
      k                  ! Loop index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Eliminate some compiler warnings
    if ( .false. ) then
      KK_mean_vol_rad = dzq
      KK_mean_vol_rad = rvm
      KK_mean_vol_rad = w_std_dev
      KK_mean_vol_rad = wm_zt
      KK_mean_vol_rad = cloud_frac
    end if
    ! Initialize microphys_stats_vars for statistics sampling
    call microphys_stats_alloc( nz, num_stats_zt, microphys_stats_zt )
    call microphys_stats_alloc( 1, num_stats_sfc, microphys_stats_sfc )

    !!! Initialize microphysics fields.
    call KK_microphys_init( nz, hydromet, &
                            hydromet_mc, hydromet_vel, rrm, Nrm, &
                            KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                            KK_mean_vol_rad, KK_Nrm_evap_tndcy, &
                            KK_Nrm_auto_tndcy, &
                            l_src_adj_enabled, l_evap_adj_enabled )
    l_clip_positive_sed = .true.

    ! The time tendency of cloud droplet concentration is only present because
    ! of Latin Hypercube interface.  Ncm_mc is needed for other microphysics
    ! schemes, but not KK.  Simply set Ncm_mc to 0.
    Ncm_mc = zero

    ! Do not use these adjustments when l_silhs_KK_convergence_adj_mean is
    ! true. We need to do this to get Latin hypercube to converge to KK
    ! analytic. The adjustment code will be called by Latin hypercube for the
    ! means only.
    if ( l_silhs_KK_convergence_adj_mean .and. l_latin_hypercube ) then
      l_src_adj_enabled = .false.
      l_evap_adj_enabled = .false.
      l_clip_positive_sed = .false.
    end if

    !!! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

       !!! Calculate the coefficients for the KK microphysics tendencies.
       call KK_tendency_coefs( thlm(k), exner(k), p_in_Pa(k), rho(k), &
                               KK_evap_coef, KK_auto_coef, &
                               KK_accr_coef, KK_mvr_coef )


       !!! Calculate the local KK rain drop mean volume radius.
       if ( rrm(k) > rr_tol ) then

          KK_mean_vol_rad(k)  &
          = KK_mvr_local_mean( rrm(k), Nrm(k), KK_mvr_coef )

       else  ! r_r or N_r = 0.

          KK_mean_vol_rad(k) = zero

       endif

       !!! Calculate the values of the local KK microphysics tendencies.

       !!! Calculate the local KK evaporation tendency.
       if ( rrm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

          KK_evap_tndcy(k)  &
          = KK_evap_local_mean( chi(k), rrm(k), Nrm(k), &
                                KK_evap_coef )

       else  ! r_r or N_r = 0.

          KK_evap_tndcy(k) = zero

       endif

       !!! Calculate the local KK autoconversion tendency.
       if ( Ncm(k) > Nc_tol ) then

          KK_auto_tndcy(k)  &
          = KK_auto_local_mean( chi(k), Ncm(k), KK_auto_coef )

       else  ! N_c = 0.

          KK_auto_tndcy(k) = zero

       endif

       !!! Calculate the local KK accretion tendency.
       if ( rrm(k) > rr_tol ) then

          KK_accr_tndcy(k)  &
          = KK_accr_local_mean( chi(k), rrm(k), KK_accr_coef )

       else  ! r_r = 0.

          KK_accr_tndcy(k) = zero

       endif


       !!! KK rain drop concentration microphysics tendencies.

       !!! Calculate the KK N_r evaporation tendency.
       if ( rrm(k) > rr_tol .and. Nrm(k) > Nr_tol ) then

          KK_Nrm_evap_tndcy(k)  &
          = KK_Nrm_evap_local_mean( KK_evap_tndcy(k), Nrm(k), rrm(k), dt )

       else  ! r_r or N_r = 0.

          KK_Nrm_evap_tndcy(k) = zero

       endif

       !!! Calculate the KK N_r autoconversion tendency.
       KK_Nrm_auto_tndcy(k) = KK_Nrm_auto_mean( KK_auto_tndcy(k) )


       !!! Calculate any necessary adjustments to KK microphysics tendencies.
       call KK_microphys_adjust( dt, exner(k), rcm(k), rrm(k), Nrm(k), &
                                 KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                 KK_accr_tndcy(k), KK_Nrm_evap_tndcy(k), &
                                 KK_Nrm_auto_tndcy(k), l_src_adj_enabled, &
                                 l_evap_adj_enabled, &
                                 rrm_mc(k), Nrm_mc(k), &
                                 rvm_mc(k), rcm_mc(k), thlm_mc(k), &
                                 adj_terms(k) )

    enddo  ! Microphysics tendency loop: k = 2, nz, 1


    !!! Boundary conditions for microphysics tendencies.

    ! Explicit contributions to rrm and Nrm from microphysics are not set at
    ! thermodynamic level k = 1 because it is below the model lower boundary.
    rrm_mc(1) = zero
    Nrm_mc(1) = zero

    rrm_mc(nz) = zero
    Nrm_mc(nz) = zero

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
    call KK_sedimentation( nz, cloud_top_level, KK_mean_vol_rad, Vrr, VNr, &
                           l_clip_positive_sed )

    !!! Output hydrometeor mean tendencies and mean sedimentation velocities
    !!! in output arrays.
    call KK_microphys_output( nz, Vrr, VNr, rrm_mc, Nrm_mc, &
                              hydromet_mc, hydromet_vel )

    !!! Output values for statistics
    call microphys_put_var( irrm_cond, KK_evap_tndcy, microphys_stats_zt )
    call microphys_put_var( irrm_auto, KK_auto_tndcy, microphys_stats_zt )
    call microphys_put_var( irrm_accr, KK_accr_tndcy, microphys_stats_zt )
    call microphys_put_var( im_vol_rad_rain, KK_mean_vol_rad, microphys_stats_zt )
    call microphys_put_var( iNrm_cond, KK_Nrm_evap_tndcy, microphys_stats_zt )
    call microphys_put_var( iNrm_auto, KK_Nrm_auto_tndcy, microphys_stats_zt )
    if ( l_src_adj_enabled ) then
      call microphys_put_var( irrm_src_adj, adj_terms%rrm_src_adj, microphys_stats_zt )
      call microphys_put_var( iNrm_src_adj, adj_terms%Nrm_src_adj, microphys_stats_zt )
    end if
    if ( l_evap_adj_enabled ) then
      call microphys_put_var( irrm_cond_adj, adj_terms%rrm_cond_adj, microphys_stats_zt )
      call microphys_put_var( iNrm_cond_adj, adj_terms%Nrm_cond_adj, microphys_stats_zt )
    end if
    call microphys_put_var( irrm_mc_nonadj, KK_auto_tndcy+KK_accr_tndcy+KK_evap_tndcy, &
                            microphys_stats_zt )

    return

  end subroutine KK_local_microphys

  !=============================================================================
  subroutine KK_upscaled_microphys( dt, nz, d_variables, l_stats_samp, & ! In
                                    wm_zt, rtm, thlm, p_in_Pa,         & ! In
                                    exner, rho, rcm, Nc_in_cloud,      & ! In
                                    pdf_params, hydromet_pdf_params,   & ! In
                                    hydromet,                          & ! In
                                    mu_x_1_n, mu_x_2_n,                & ! In
                                    sigma_x_1_n, sigma_x_2_n,          & ! In
                                    corr_array_1_n, corr_array_2_n,    & ! In
                                    hydromet_mc, hydromet_vel,         & ! Out
                                    rcm_mc, rvm_mc, thlm_mc,           & ! Out
                                    hydromet_vel_covar_zt_impc,        & ! Out
                                    hydromet_vel_covar_zt_expc,        & ! Out
                                    wprtp_mc, wpthlp_mc, rtp2_mc,      & ! Out
                                    thlp2_mc, rtpthlp_mc )               ! Out

    ! Description:
    ! Version of KK microphysics scheme that is analytically upscaled by
    ! integrating over the product of the microphysics tendency and the
    ! functional form of the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        zero    ! Constant(s)

    use model_flags, only: &
        l_const_Nc_in_cloud

    use parameters_microphys, only: &
        l_var_covar_src        ! Flag for using variance/covariance src terms

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
        iirrm, & ! Constant(s)
        iiNrm

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision

    use stats_type_utilities, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        iVrrprrp_expcalc, &
        iVNrpNrp_expcalc, &
        stats_zm, &
        stats_zt, &
        irrm_src_adj, &
        iNrm_src_adj, &
        irrm_cond_adj, &
        iNrm_cond_adj, &
        irrm_mc_nonadj

    implicit none

    ! Local Constants
    logical, parameter :: &
      l_clip_positive_sed = .true.  ! Clip positive Vrr and VNr terms to zero

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation arrays

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wm_zt,       & ! Mean vertical velocity on thermodynamic levels  [m/s]
      rtm,         & ! Mean total water mixing ratio                   [kg/kg]
      thlm,        & ! Mean liquid water potential temperature         [K]
      p_in_Pa,     & ! Pressure                                        [Pa]
      exner,       & ! Exner function                                  [-]
      rho,         & ! Density                                         [kg/m^3]
      rcm,         & ! Mean cloud water mixing ratio                   [kg/kg]
      Nc_in_cloud    ! Constant in-cloud value of cloud droplet conc.  [num/kg]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet       ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]

    real( kind = core_rknd ), dimension(d_variables, nz), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      wprtp_mc,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      rrm,    & ! Mean rain water mixing ratio, < r_r > [kg/kg]
      Nrm,    & ! Mean rain drop concentration, < N_r > [num/kg]
      Vrr,    & ! Mean sedimentation velocity of r_r    [m/s]
      VNr,    & ! Mean sedimentation velocity of N_r    [m/s]
      rrm_mc, & ! Mean (dr_r/dt) due to microphysics    [(kg/kg)/s]
      Nrm_mc    ! Mean (dN_r/dt) due to microphysics    [(num/kg)/s]

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
      mu_chi_1,      & ! Mean of chi (1st PDF component)                 [kg/kg]
      mu_chi_2,      & ! Mean of chi (2nd PDF component)                 [kg/kg]
      mu_eta_1,      & ! Mean of eta (1st PDF component)                 [kg/kg]
      mu_eta_2,      & ! Mean of eta (2nd PDF component)                 [kg/kg]
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
      sigma_chi_1,   & ! Standard deviation of chi (1st PDF component)   [kg/kg]
      sigma_chi_2,   & ! Standard deviation of chi (2nd PDF component)   [kg/kg]
      sigma_eta_1,   & ! Standard deviation of eta (1st PDF component)   [kg/kg]
      sigma_eta_2,   & ! Standard deviation of eta (2nd PDF component)   [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      sigma_Ncn_1_n, & ! Standard deviation of ln Ncn (1st PDF component)    [-]
      sigma_Ncn_2_n    ! Standard deviation of ln Ncn (2nd PDF component)    [-]

    real( kind = core_rknd ) :: &
      corr_w_chi_1,   & ! Correlation of w and chi (1st PDF component)     [-]
      corr_w_chi_2,   & ! Correlation of w and chi (2nd PDF component)     [-]
     !corr_w_rr_1,    & ! Correlation of w and rr (1st PDF component) ip   [-]
     !corr_w_rr_2,    & ! Correlation of w and rr (2nd PDF component) ip   [-]
     !corr_w_Nr_1,    & ! Correlation of w and Nr (1st PDF component) ip   [-]
     !corr_w_Nr_2,    & ! Correlation of w and Nr (2nd PDF component) ip   [-]
     !corr_w_Ncn_1,   & ! Correlation of w and Ncn (1st PDF component)     [-]
     !corr_w_Ncn_2,   & ! Correlation of w and Ncn (2nd PDF component)     [-]
      corr_chi_eta_1, & ! Correlation of chi and eta (1st PDF component)   [-]
      corr_chi_eta_2!,& ! Correlation of chi and eta (2nd PDF component)   [-]
     !corr_chi_rr_1,  & ! Correlation of chi and rr (1st PDF component) ip [-]
     !corr_chi_rr_2,  & ! Correlation of chi and rr (2nd PDF component) ip [-]
     !corr_chi_Nr_1,  & ! Correlation of chi and Nr (1st PDF component) ip [-]
     !corr_chi_Nr_2,  & ! Correlation of chi and Nr (2nd PDF component) ip [-]
     !corr_chi_Ncn_1, & ! Correlation of chi and Ncn (1st PDF component)   [-]
     !corr_chi_Ncn_2, & ! Correlation of chi and Ncn (2nd PDF component)   [-]
     !corr_eta_rr_1,  & ! Correlation of eta and rr (1st PDF component) ip [-]
     !corr_eta_rr_2,  & ! Correlation of eta and rr (2nd PDF component) ip [-]
     !corr_eta_Nr_1,  & ! Correlation of eta and Nr (1st PDF component) ip [-]
     !corr_eta_Nr_2,  & ! Correlation of eta and Nr (2nd PDF component) ip [-]
     !corr_eta_Ncn_1, & ! Correlation of eta and Ncn (1st PDF component)   [-]
     !corr_eta_Ncn_2, & ! Correlation of eta and Ncn (2nd PDF component)   [-]
     !corr_rr_Nr_1,   & ! Correlation of rr and Nr (1st PDF component) ip  [-]
     !corr_rr_Nr_2      ! Correlation of rr and Nr (2nd PDF component) ip  [-]

    real( kind = core_rknd ) :: &
      corr_w_rr_1_n,    & ! Correlation of w and ln rr (1st PDF comp.) ip    [-]
      corr_w_rr_2_n,    & ! Correlation of w and ln rr (2nd PDF comp.) ip    [-]
      corr_w_Nr_1_n,    & ! Correlation of w and ln Nr (1st PDF comp.) ip    [-]
      corr_w_Nr_2_n,    & ! Correlation of w and ln Nr (2nd PDF comp.) ip    [-]
      corr_w_Ncn_1_n,   & ! Correlation of w and ln Ncn (1st PDF component)  [-]
      corr_w_Ncn_2_n,   & ! Correlation of w and ln Ncn (2nd PDF component)  [-]
      corr_chi_rr_1_n,  & ! Correlation of chi and ln rr (1st PDF comp.) ip  [-]
      corr_chi_rr_2_n,  & ! Correlation of chi and ln rr (2nd PDF comp.) ip  [-]
      corr_chi_Nr_1_n,  & ! Correlation of chi and ln Nr (1st PDF comp.) ip  [-]
      corr_chi_Nr_2_n,  & ! Correlation of chi and ln Nr (2nd PDF comp.) ip  [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      corr_eta_rr_1_n,  & ! Correlation of eta and ln rr (1st PDF comp.) ip  [-]
      corr_eta_rr_2_n,  & ! Correlation of eta and ln rr (2nd PDF comp.) ip  [-]
      corr_eta_Nr_1_n,  & ! Correlation of eta and ln Nr (1st PDF comp.) ip  [-]
      corr_eta_Nr_2_n,  & ! Correlation of eta and ln Nr (2nd PDF comp.) ip  [-]
      corr_eta_Ncn_1_n, & ! Correlation of eta and ln Ncn (1st PDF comp.)    [-]
      corr_eta_Ncn_2_n, & ! Correlation of eta and ln Ncn (2nd PDF comp.)    [-]
      corr_rr_Nr_1_n,   & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rr_Nr_2_n      ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ) :: &
      rr_1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr_2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr_1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr_2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    real( kind = core_rknd ), dimension(nz) :: &
      Vrrprrp_zt_impc, & ! Imp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)]
      Vrrprrp_zt_expc, & ! Exp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)(kg/kg)]
      VNrpNrp_zt_impc, & ! Imp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)]
      VNrpNrp_zt_expc    ! Exp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz) :: &
      wprtp_mc_zt,   & ! Micro. tend. for <w'rt'>; t-lev   [m*(kg/kg)/s^2]
      wpthlp_mc_zt,  & ! Micro. tend. for <w'thl'>; t-lev  [m*K/s^2]
      rtp2_mc_zt,    & ! Micro. tend. for <rt'^2>; t-lev   [(kg/kg)^2/s]
      thlp2_mc_zt,   & ! Micro. tend. for <thl'^2>; t-lev  [K^2/s]
      rtpthlp_mc_zt    ! Micro. tend. for <rt'thl'>; t-lev [K*(kg/kg)/s]

    type(KK_microphys_adj_terms_type), dimension(nz) :: &
      adj_terms           ! Adjustments to microphysics terms

    logical :: &
      l_src_adj_enabled,  & ! Flag to enable rrm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrm/Nrm evaporation adjustment

    integer :: &
      cloud_top_level, & ! Vertical level index of cloud top 
      k                  ! Loop index

    ! ---- Begin Code ----

    !!! Initialize microphysics fields.
    call KK_microphys_init( nz, hydromet, &
                            hydromet_mc, hydromet_vel, rrm, Nrm, &
                            KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, &
                            KK_mean_vol_rad, KK_Nrm_evap_tndcy, &
                            KK_Nrm_auto_tndcy, &
                            l_src_adj_enabled, l_evap_adj_enabled )

    ! Setup mixture fraction.
    mixt_frac = pdf_params%mixt_frac    

    !!! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

       !!! Calculate the coefficients for the KK microphysics tendencies.
       call KK_tendency_coefs( thlm(k), exner(k), p_in_Pa(k), rho(k), &
                               KK_evap_coef, KK_auto_coef, &
                               KK_accr_coef, KK_mvr_coef )


       !!! Unpack the PDF parameters.
       call unpack_pdf_params_KK( d_variables, mu_x_1_n(:,k), mu_x_2_n(:,k), &
                                  sigma_x_1_n(:,k), sigma_x_2_n(:,k), &
                                  corr_array_1_n(:,:,k), &
                                  corr_array_2_n(:,:,k), &
                                  hydromet_pdf_params(k), &
                                  mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, &
                                  mu_eta_1, mu_eta_2, mu_rr_1, mu_rr_2, &
                                  mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                  mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                  mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                  sigma_w_1, sigma_w_2, sigma_chi_1, &
                                  sigma_chi_2, sigma_eta_1, sigma_eta_2, &
                                  sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                  sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                  sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                  sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                  corr_w_chi_1, corr_w_chi_2, &
                                  corr_chi_eta_1, corr_chi_eta_2, &
                                  corr_w_rr_1_n, corr_w_rr_2_n, &
                                  corr_w_Nr_1_n, corr_w_Nr_2_n, &
                                  corr_w_Ncn_1_n, corr_w_Ncn_2_n, &
                                  corr_chi_rr_1_n, corr_chi_rr_2_n, &
                                  corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                                  corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                  corr_eta_rr_1_n, corr_eta_rr_2_n, &
                                  corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                                  corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                                  corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                  rr_1, rr_2, Nr_1, Nr_2, &
                                  precip_frac, precip_frac_1, precip_frac_2 )

       !!! Calculate the values of the upscaled KK microphysics tendencies.
       call KK_upscaled_means_driver( mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, &
                                      mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                      mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                      mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                      sigma_chi_1, sigma_chi_2, &
                                      sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                      sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                      sigma_rr_1_n, sigma_rr_2_n, &
                                      sigma_Nr_1_n, sigma_Nr_2_n, &
                                      sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                      corr_chi_rr_1_n, corr_chi_rr_2_n, &
                                      corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                                      corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                      corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                      mixt_frac(k), precip_frac_1, &
                                      precip_frac_2, KK_evap_coef, &
                                      KK_auto_coef, KK_accr_coef, &
                                      KK_mvr_coef, KK_evap_tndcy(k), &
                                      KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                      KK_mean_vol_rad(k) )

       !!! Calculate the values of the turbulent sedimentation terms,
       !!! <V_rr'rr'> and <V_Nr'Nr'>.  Each turbulent sedimentation term has an
       !!! implicit component and an explicit component to be used in the
       !!! predictive equation for the hydrometeor.
       call KK_sed_vel_covars( rrm(k), rr_1, rr_2, Nrm(k), &
                               Nr_1, Nr_2, KK_mean_vol_rad(k), &
                               mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                               mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                               sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                               sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                               sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                               KK_mvr_coef, mixt_frac(k), precip_frac_1, &
                               precip_frac_2, k, l_stats_samp, &
                               Vrrprrp_zt_impc(k), Vrrprrp_zt_expc(k), &
                               VNrpNrp_zt_impc(k), VNrpNrp_zt_expc(k) )

       if ( l_var_covar_src ) then

          call KK_upscaled_covar_driver( wm_zt(k), rtm(k), thlm(k), &
                                         exner(k), rrm(k), Nrm(k), &
                                         mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, &
                                         mu_eta_1, mu_eta_2, mu_rr_1, mu_rr_2, &
                                         mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                         mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                         mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                         sigma_w_1, sigma_w_2, sigma_chi_1, &
                                         sigma_chi_2, sigma_eta_1, sigma_eta_2, &
                                         sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                         sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                         sigma_rr_1_n, sigma_rr_2_n, &
                                         sigma_Nr_1_n, sigma_Nr_2_n, &
                                         sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                         corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, &
                                         corr_w_rr_2_n, corr_w_Nr_1_n, &
                                         corr_w_Nr_2_n, corr_w_Ncn_1_n, &
                                         corr_w_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, &
                                         corr_chi_rr_1_n, corr_chi_rr_2_n, &
                                         corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                                         corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                         corr_eta_rr_1_n, corr_eta_rr_2_n, &
                                         corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                                         corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                                         corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                         mixt_frac(k), precip_frac_1, &
                                         precip_frac_2, Nc_in_cloud(k), &
                                         KK_evap_coef, KK_auto_coef, &
                                         KK_accr_coef, KK_evap_tndcy(k), &
                                         KK_auto_tndcy(k), KK_accr_tndcy(k), &
                                         pdf_params(k), k, &
                                         l_const_Nc_in_cloud, l_stats_samp, &
                                         wprtp_mc_zt(k), &
                                         wpthlp_mc_zt(k), &
                                         rtp2_mc_zt(k), &
                                         thlp2_mc_zt(k), &
                                         rtpthlp_mc_zt(k) )

       endif

       !!! KK rain drop concentration microphysics tendencies.

       !!! Calculate the KK N_r evaporation tendency.
       KK_Nrm_evap_tndcy(k) &
       = KK_Nrm_evap_upscaled_mean( mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, &
                                    mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                                    mu_Nr_1_n, mu_Nr_2_n, sigma_chi_1, &
                                    sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                                    sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                    sigma_rr_2_n, sigma_Nr_1_n, &
                                    sigma_Nr_2_n, corr_chi_rr_1_n, &
                                    corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                                    corr_chi_Nr_2_n, corr_rr_Nr_1_n, &
                                    corr_rr_Nr_2_n, KK_evap_coef, mixt_frac(k), &
                                    precip_frac_1, precip_frac_2, dt )


       !!! Calculate the KK N_r autoconversion tendency.
       KK_Nrm_auto_tndcy(k) = KK_Nrm_auto_mean( KK_auto_tndcy(k) )


       !!! Calculate any necessary adjustments to KK microphysics tendencies.
       call KK_microphys_adjust( dt, exner(k), rcm(k), rrm(k), Nrm(k), &
                                 KK_evap_tndcy(k), KK_auto_tndcy(k), &
                                 KK_accr_tndcy(k), KK_Nrm_evap_tndcy(k), &
                                 KK_Nrm_auto_tndcy(k), l_src_adj_enabled, &
                                 l_evap_adj_enabled, &
                                 rrm_mc(k), Nrm_mc(k), &
                                 rvm_mc(k), rcm_mc(k), thlm_mc(k), &
                                 adj_terms(k) )


       !!! Statistical output for upscaled KK.
       call KK_upscaled_stats( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                               mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                               sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                               sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                               sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, &
                               corr_rr_Nr_2_n, mixt_frac(k), precip_frac_1, &
                               precip_frac_2, KK_mvr_coef, KK_mean_vol_rad(k), &
                               k, l_stats_samp )

       !!! Statistical output for mean microphysics tendenices.
       call KK_stats_output( KK_evap_tndcy(k), KK_auto_tndcy(k), &
                             KK_accr_tndcy(k), KK_mean_vol_rad(k), &
                             KK_Nrm_evap_tndcy(k), KK_Nrm_auto_tndcy(k), &
                             l_stats_samp, k )


    enddo  ! Microphysics tendency loop: k = 2, nz, 1


    if ( l_var_covar_src ) then

       ! Set values of wprtp_mc_zt, wpthlp_mc_zt, rtp2_mc_zt, thlp2_mc_zt, and
       ! rtpthlp_mc_zt to 0 at the lowest thermodynamic grid level
       ! (thermodynamic level 1), which is below the model lower boundary.  This
       ! will prevent an unset value of these variables (from thermodynamic
       ! level 1, which is not part of the above microphysics tendency loop)
       ! from being used in the interpolation of these variables to momentum
       ! levels (in this case, the model lower boundary at momentum level 1).
       ! The interpolated value of these variables at momentum level 1 is not
       ! used in the model code.
       wprtp_mc_zt(1)   = zero
       wpthlp_mc_zt(1)  = zero
       rtp2_mc_zt(1)    = zero
       thlp2_mc_zt(1)   = zero
       rtpthlp_mc_zt(1) = zero

       ! Output microphysics tendency terms for
       ! model variances and covariances on momentum levels.
       wprtp_mc   = zt2zm( wprtp_mc_zt )
       wpthlp_mc  = zt2zm( wpthlp_mc_zt )
       rtp2_mc    = zt2zm( rtp2_mc_zt )
       thlp2_mc   = zt2zm( thlp2_mc_zt )
       rtpthlp_mc = zt2zm( rtpthlp_mc_zt )

       ! Set values of microphysics tendency terms to 0 at model lower boundary.
       wprtp_mc(1)   = zero
       wpthlp_mc(1)  = zero
       rtp2_mc(1)    = zero
       thlp2_mc(1)   = zero
       rtpthlp_mc(1) = zero
       ! Set values of microphysics tendency terms to 0 at model upper boundary.
       wprtp_mc(nz)   = zero
       wpthlp_mc(nz)  = zero
       rtp2_mc(nz)    = zero
       thlp2_mc(nz)   = zero
       rtpthlp_mc(nz) = zero

    else

       ! Microphysics tendency terms for model variances and covariances
       ! are set to 0.
       wprtp_mc   = zero
       wpthlp_mc  = zero
       rtp2_mc    = zero
       thlp2_mc   = zero
       rtpthlp_mc = zero

    endif

    !!! Boundary conditions for microphysics tendencies.

    ! Explicit contributions to rrm and Nrm from microphysics are not set at
    ! thermodynamic level k = 1 because it is below the model lower boundary.
    rrm_mc(1) = zero
    Nrm_mc(1) = zero

    rrm_mc(nz) = zero
    Nrm_mc(nz) = zero

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
    call KK_sedimentation( nz, cloud_top_level, KK_mean_vol_rad, Vrr, VNr, &
                           l_clip_positive_sed )

    !!! Output hydrometeor mean tendencies and mean sedimentation velocities
    !!! in output arrays.
    call KK_microphys_output( nz, Vrr, VNr, rrm_mc, Nrm_mc, &
                              hydromet_mc, hydromet_vel )

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

    ! The implicit and explicit components used to calculate the covariances of
    ! hydrometeor sedimentation velocities and their associated hydrometeors
    ! (<V_rr'r_r'> and <V_Nr'N_r'>) are fed into the output arrays.
    hydromet_vel_covar_zt_impc(:,iirrm) = Vrrprrp_zt_impc
    hydromet_vel_covar_zt_expc(:,iirrm) = Vrrprrp_zt_expc
    hydromet_vel_covar_zt_impc(:,iiNrm) = VNrpNrp_zt_impc
    hydromet_vel_covar_zt_expc(:,iiNrm) = VNrpNrp_zt_expc

    ! Statistics
    if ( l_stats_samp ) then

       if ( iVrrprrp_expcalc > 0 ) then

          ! The covariance < V_rr'r_r' > calculated completely explicitly.
          ! When semi-implicit turbulent advection is used, this result can be
          ! compared to the < V_rr'r_r' > results used in the code, which are
          ! calculated semi-implicitly.
          call stat_update_var( iVrrprrp_expcalc, &
                                zt2zm( Vrrprrp_zt_impc * rrm &
                                       + Vrrprrp_zt_expc ), stats_zm )

       endif

       if ( iVNrpNrp_expcalc > 0 ) then

          ! The covariance < V_Nr'N_r' > calculated completely explicitly.
          ! When semi-implicit turbulent advection is used, this result can be
          ! compared to the < V_Nr'N_r' > results used in the code, which are
          ! calculated semi-implicitly.
          call stat_update_var( iVNrpNrp_expcalc, &
                                zt2zm( VNrpNrp_zt_impc * Nrm &
                                       + VNrpNrp_zt_expc ), stats_zm )

       endif

       ! Stats output for microphysics adjustment terms
       if ( irrm_src_adj > 0 ) then
         call stat_update_var( irrm_src_adj, adj_terms%rrm_src_adj, stats_zt )
       end if

       if ( iNrm_src_adj > 0 ) then
         call stat_update_var( iNrm_src_adj, adj_terms%Nrm_src_adj, stats_zt )
       end if

       if ( irrm_cond_adj > 0 ) then
         call stat_update_var( irrm_cond_adj, adj_terms%rrm_cond_adj, stats_zt )
       end if

       if ( iNrm_cond_adj > 0 ) then
         call stat_update_var( iNrm_cond_adj, adj_terms%Nrm_cond_adj, stats_zt )
       end if

       if ( irrm_mc_nonadj > 0 ) then
         call stat_update_var( irrm_mc_nonadj, KK_auto_tndcy+KK_accr_tndcy+KK_evap_tndcy, &
                               stats_zt )
       end if

    endif ! l_stats_samp


    return

  end subroutine KK_upscaled_microphys

  !=============================================================================
  subroutine KK_microphys_init( nz, hydromet, &
                                hydromet_mc, hydromet_vel, rrm, Nrm, &
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
        iirrm, & ! Constant(s)
        iiNrm

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Hydrometeor species                      [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz), intent(out) ::  &
      rrm, & ! Mean rain water mixing ratio (overall), < r_r >    [kg/kg]
      Nrm    ! Mean rain drop concentration (overall), < N_r >    [num/kg]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      KK_evap_tndcy,     & ! Mean KK (dr_r/dt) due to evaporation    [(kg/kg)/s]
      KK_auto_tndcy,     & ! Mean KK (dr_r/dt) due to autoconversion [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK (dr_r/dt) due to accretion      [(kg/kg)/s]
      KK_mean_vol_rad,   & ! Mean KK rain drop mean volume radius            [m]
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation   [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.     [(num/kg)/s]

    logical, intent(out) :: &
      l_src_adj_enabled,  & ! Flag to enable rrm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrm/Nrm evaporation adjustment


    !!! Initialize microphysics tendencies and rain drop mean volume radius.
    KK_evap_tndcy = zero
    KK_auto_tndcy = zero
    KK_accr_tndcy = zero

    KK_mean_vol_rad = zero

    KK_Nrm_evap_tndcy = zero
    KK_Nrm_auto_tndcy = zero

    hydromet_mc(:,:) = zero
    hydromet_vel(:,:) = zero

    ! Set up mean field variables for <r_r> and <N_r>.
    rrm = hydromet(:,iirrm)
    Nrm = hydromet(:,iiNrm)

    ! Set KK microphysics tendency adjustment flags
    l_src_adj_enabled  = .true.
    l_evap_adj_enabled = .true.

    return

  end subroutine KK_microphys_init

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

    use parameters_KK, only: &
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
  subroutine KK_microphys_adjust( dt, exner, rcm, rrm, Nrm, &
                                  KK_evap_tndcy, KK_auto_tndcy, &
                                  KK_accr_tndcy, KK_Nrm_evap_tndcy, &
                                  KK_Nrm_auto_tndcy, l_src_adj_enabled, &
                                  l_evap_adj_enabled, &
                                  rrm_mc, Nrm_mc, &
                                  rvm_mc, rcm_mc, thlm_mc, &
                                  adj_terms )

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

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt        ! Model time step duration                 [s]

    real( kind = core_rknd ), intent(in) :: &
      exner, & ! Exner function                           [-]
      rcm,   & ! Mean cloud water mixing ratio            [kg/kg]
      rrm,   & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm      ! Mean rain drop concentration, < N_r >    [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      KK_evap_tndcy,     & ! Mean KK <r_r> evaporation tendency     [(kg/kg)/s]
      KK_auto_tndcy,     & ! Mean KK <r_r> autoconversion tendency  [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK <r_r> accretion tendency       [(kg/kg)/s]
      KK_Nrm_evap_tndcy, & ! Mean KK <N_r> evaporation tendency     [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK <N_r> autoconversion tendency  [(num/kg)/s]

    logical, intent(in) :: &
      l_src_adj_enabled,  & ! Flag to enable rrm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrm/Nrm evaporation adjustment

    ! Output Variables
    real( kind = core_rknd ), intent(out) ::  &
      rrm_mc, & ! Mean <dr_r/dt> due to microphysics       [(kg/kg)/s]
      Nrm_mc    ! Mean <dN_r/dt> due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio       [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio        [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature    [K/s]

    type(KK_microphys_adj_terms_type), intent(out) :: &
      adj_terms    ! Structure containing information about adjustment of
                      ! microphysics terms

    ! Local Variables
    real( kind = core_rknd ) ::  &
      rrm_source,   & ! Total source term rate for rrm           [(kg/kg)/s]
      Nrm_source,   & ! Total source term rate for Nrm           [(num/kg)/s]
      rrm_src_adj,  & ! Total adjustment to rrm source terms     [(kg/kg)/s]
      Nrm_src_adj,  & ! Total adjustment to Nrm source terms     [(num/kg)/s]
      rrm_evap_net, & ! Net evaporation rate of <r_r>            [(kg/kg)/s]
      Nrm_evap_net    ! Net evaporation rate of <N_r>            [(num/kg)/s]

    real( kind = core_rknd ) :: &
      rrm_src_max,    & ! Maximum allowable rrm source rate    [(kg/kg)/s]
      rrm_auto_ratio, & ! Ratio of rrm autoconv to overall source term [-]
      total_rc_needed   ! Amount of r_c needed to over the timestep
                        ! for rain source terms                    [kg/kg]
  !-----------------------------------------------------------------------

    !---- Begin Code -----

    !!! Source-adjustment code for rrm and Nrm.
    rrm_source = KK_auto_tndcy + KK_accr_tndcy
    Nrm_source = KK_Nrm_auto_tndcy

    if ( l_src_adj_enabled ) then

       ! The increase of rain due to autoconversion and accretion both draw
       ! water from the available cloud water.  Over a long time step, these
       ! rates may over-deplete cloud water.  In other words, these processes
       ! may draw more cloud water than there is available.  Thus, the total
       ! source rate multiplied by the duration of the time step cannot exceed
       ! the total amount of cloud water available.  If it does, then the rate
       ! must be adjusted.
       total_rc_needed = rrm_source * dt

       if ( total_rc_needed > rcm ) then

          ! The maximum allowable rate of the source terms is rcm/dt.
          rrm_src_max = rcm / dt

          ! The amount of adjustment to the source terms.
          ! This value should always be negative.
          rrm_src_adj = rrm_src_max - rrm_source

          ! Reset the value of the source terms to the maximum allowable value
          ! of the source terms.
          rrm_source = rrm_src_max

          ! The rrm source terms are made up of autoconversion and accretion.
          ! Only the sum of those two terms is corrected.  However, Nrm has only
          ! an autoconversion term for a source term.  Assume that change in the
          ! rrm autoconversion term is proportional to the total rrm adjustment
          ! rate by the ratio of rrm autoconversion to the overall source term.
          ! Then, plug the rrm autoconversion adjustment into the equation for
          ! Nrm autoconversion to determine the effect on the Nrm source term.
          rrm_auto_ratio = KK_auto_tndcy /  &
                              ( KK_auto_tndcy + KK_accr_tndcy )
          Nrm_src_adj = KK_Nrm_auto_mean( rrm_auto_ratio * rrm_src_adj )

          ! Change Nrm by Nrm_src_adj.  Nrm_src_adj will always be negative.
          Nrm_source = Nrm_source + Nrm_src_adj

       else ! total_rc_needed <= rcm: enough available cloud water.

          rrm_src_adj = zero
          Nrm_src_adj = zero

       endif

    else ! l_src_adj_enabled is false: source adjustment is disabled.

       rrm_src_adj = zero
       Nrm_src_adj = zero

    endif


    !!! Evaporation-adjustment code for rrm and Nrm.
    if ( l_evap_adj_enabled ) then

       ! Prevent over-evaporation of rain over a long model time step.
       ! Limit the evaporation rate of both rain water mixing ratio and rain
       ! drop concentration.  The total amount of rain lost due to evaporation
       ! cannot be so great as to result in negative rain.

       ! Calculate net evaporation rate of <r_r>.
       rrm_evap_net = max( KK_evap_tndcy, &
                              - rrm / dt )

       ! Recalcuate the net evaporation rate of <N_r> based on the net
       ! evaporation rate of <r_r>.
       if ( KK_evap_tndcy /= rrm_evap_net .and. &
            rrm > rr_tol .and. Nrm > Nr_tol ) then
          Nrm_evap_net &
          = KK_Nrm_evap_local_mean( rrm_evap_net, Nrm, rrm, dt )
       else
          Nrm_evap_net = KK_Nrm_evap_tndcy
       endif

       Nrm_evap_net = max( Nrm_evap_net, &
                           - Nrm / dt )

    else ! l_evap_adj_enabled is false: evaporation adjustment is disabled.

       rrm_evap_net = KK_evap_tndcy
       Nrm_evap_net = KK_Nrm_evap_tndcy

    endif


    !!! Calculate overall KK microphysics tendencies.
    rrm_mc = rrm_evap_net + rrm_source
    Nrm_mc = Nrm_evap_net + Nrm_source

    !!! Explicit contributions to thlm and rtm from the microphysics
    rvm_mc  = -rrm_evap_net
    rcm_mc  = -rrm_source  ! Accretion + Autoconversion
    thlm_mc = ( Lv / ( Cp * exner ) ) * rrm_mc

    ! Return adjustment terms
    adj_terms%rrm_src_adj = rrm_src_adj
    adj_terms%Nrm_src_adj = Nrm_src_adj
    adj_terms%rrm_cond_adj = rrm_evap_net - KK_evap_tndcy
    adj_terms%Nrm_cond_adj = Nrm_evap_net - KK_Nrm_evap_tndcy

    return

  end subroutine KK_microphys_adjust

  !=============================================================================
  subroutine KK_upscaled_stats( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, &
                                sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, &
                                corr_rr_Nr_2_n, mixt_frac, precip_frac_1, &
                                precip_frac_2, KK_mvr_coef, KK_mean_vol_rad, &
                                level, l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use KK_upscaled_variances, only: &
        variance_KK_mvr    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        iKK_mvr_variance_zt, &
        stats_zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1,        & ! Mean of rr (1st PDF component) in-precip (ip)  [kg/kg]
      mu_rr_2,        & ! Mean of rr (2nd PDF component) ip              [kg/kg]
      mu_Nr_1,        & ! Mean of Nr (1st PDF component) ip             [num/kg]
      mu_Nr_2,        & ! Mean of Nr (2nd PDF component) ip             [num/kg]
      mu_rr_1_n,      & ! Mean of ln rr (1st PDF component) ip       [ln(kg/kg)]
      mu_rr_2_n,      & ! Mean of ln rr (2nd PDF component) ip       [ln(kg/kg)]
      mu_Nr_1_n,      & ! Mean of ln Nr (1st PDF component) ip      [ln(num/kg)]
      mu_Nr_2_n,      & ! Mean of ln Nr (2nd PDF component) ip      [ln(num/kg)]
      sigma_rr_1,     & ! Standard deviation of rr (1st PDF comp.) ip    [kg/kg]
      sigma_rr_2,     & ! Standard deviation of rr (2nd PDF comp.) ip    [kg/kg]
      sigma_Nr_1,     & ! Standard deviation of Nr (1st PDF comp.) ip   [num/kg]
      sigma_Nr_2,     & ! Standard deviation of Nr (2nd PDF comp.) ip   [num/kg]
      sigma_rr_1_n,   & ! Standard deviation of ln rr (1st PDF component) ip [-]
      sigma_rr_2_n,   & ! Standard deviation of ln rr (2nd PDF component) ip [-]
      sigma_Nr_1_n,   & ! Standard deviation of ln Nr (1st PDF component) ip [-]
      sigma_Nr_2_n,   & ! Standard deviation of ln Nr (2nd PDF component) ip [-]
      corr_rr_Nr_1_n, & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip    [-]
      corr_rr_Nr_2_n, & ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip    [-]
      mixt_frac,      & ! Mixture fraction                                   [-]
      precip_frac_1,  & ! Precipitation fraction (1st PDF component)         [-]
      precip_frac_2     ! Precipitation fraction (2nd PDF component)         [-]

    real( kind = core_rknd ), intent(in) :: &
      KK_mvr_coef,     & ! KK mean volume radius coefficient            [m]
      KK_mean_vol_rad    ! KK rain drop mean volume radius              [m]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variables
    real( kind = core_rknd ) :: &
      KK_mvr_variance    ! Variance of KK rain drop mean vol rad   [m^2]


    !!! Output the statistics for upscaled KK.

    ! Statistics
    if ( l_stats_samp ) then

       ! Variance of KK rain drop mean volume radius.
       if ( iKK_mvr_variance_zt > 0 ) then

          ! Calculate the variance of KK rain drop mean volume radius,
          ! < R_vr'^2 >.
          KK_mvr_variance &
          = variance_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                             mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                             sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                             sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                             sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                             KK_mean_vol_rad, KK_mvr_coef, mixt_frac, &
                             precip_frac_1, precip_frac_2 )

          call stat_update_var_pt( iKK_mvr_variance_zt, level, &
                                   KK_mvr_variance, stats_zt )

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
        stats_zt,              & ! Variable(s)
        im_vol_rad_rain, &
        irrm_cond,       &
        irrm_auto,       &
        irrm_accr,       &
        iNrm_cond,       &
        iNrm_auto

    use stats_type_utilities, only: &
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
       if ( irrm_cond > 0 ) then
          call stat_update_var_pt( irrm_cond, level, KK_evap_tndcy, stats_zt )
       endif

       if ( irrm_auto > 0 ) then
          call stat_update_var_pt( irrm_auto, level, KK_auto_tndcy, stats_zt )
       endif

       if ( irrm_accr > 0 ) then
          call stat_update_var_pt( irrm_accr, level, KK_accr_tndcy, stats_zt )
       endif

       ! Rain drop mean volume radius.
       if ( im_vol_rad_rain > 0 ) then
          call stat_update_var_pt( im_vol_rad_rain, level, KK_mean_vol_rad, stats_zt )
       endif

       ! Mean rain drop concentration microphysics tendencies.
       if ( iNrm_cond > 0 ) then
          call stat_update_var_pt( iNrm_cond, level, KK_Nrm_evap_tndcy, stats_zt )
       endif

       if ( iNrm_auto > 0 ) then
          call stat_update_var_pt( iNrm_auto, level, KK_Nrm_auto_tndcy, stats_zt )
       endif

    endif ! l_stats_samp


    return

  end subroutine KK_stats_output

  !=============================================================================
  subroutine KK_sedimentation( nz, cloud_top_level, KK_mean_vol_rad, Vrr, VNr, &
                               l_clip_positive_sed )

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
      nz,              &  ! Number of model vertical grid levels
      cloud_top_level     ! Vertical level index of cloud top

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      KK_mean_vol_rad     ! KK rain drop mean volume radius       [m]

    logical, intent(in) :: &
      l_clip_positive_sed ! Clip positive values of Vrr and VNr   [T/F]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(:), intent(inout) ::  &
      Vrr, & ! Mean sedimentation velocity of < r_r >             [m/s]
      VNr    ! Mean sedimentation velocity of < N_r >             [m/s]

    ! Local Variables
    integer :: k ! Loop iterator


    !!! Mean sedimentation velocities
    do k = 1, nz-1, 1

       ! Mean sedimentation velocity of rain water mixing ratio.
       Vrr(k) = - ( 0.012_core_rknd * ( micron_per_m * KK_mean_vol_rad(k) ) &
                    - 0.2_core_rknd )

       ! Mean sedimentation velocity of rain drop concentration.
       VNr(k) = - ( 0.007_core_rknd * ( micron_per_m * KK_mean_vol_rad(k) )  &
                    - 0.1_core_rknd )

       if ( l_clip_positive_sed ) then

         ! Mean sedimentation velocity of rain water mixing ratio cannot
         ! have a
         ! positive value.
         if ( Vrr(k) > zero ) then
            Vrr(k) = zero
         endif

         ! Mean sedimentation velocity of rain drop concentration cannot have a
         ! positive value.
         if ( VNr(k) > zero ) then
            VNr(k) = zero
         endif

       end if

    enddo ! Sedimentation velocity loop: k = 1, nz-1, 1

    if ( l_clip_positive_sed ) then
      !!! Mean sedimentation above cloud top should have a value of 0.
      if ( cloud_top_level > 1 ) then
         Vrr(cloud_top_level+1:nz-1) = zero
         VNr(cloud_top_level+1:nz-1) = zero
      endif
    end if

    !!! Boundary conditions for sedimentation velocities.

    ! The flux of rain water through the model top is 0.
    ! Vrr and VNr are set to 0 at the highest model level.
    Vrr(nz) = zero
    VNr(nz) = zero


    return

  end subroutine KK_sedimentation

  !=============================================================================
  subroutine KK_microphys_output( nz, Vrr, VNr, rrm_mc, Nrm_mc, &
                                  hydromet_mc, hydromet_vel )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        iirrm, & ! Constant(s)
        iiNrm

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      Vrr,    & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,    & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrm_mc, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc    ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]


    !!! Output mean hydrometeor tendencies and mean sedimentation velocities.

    ! Mean field tendencies.
    hydromet_mc(:,iirrm) = rrm_mc
    hydromet_mc(:,iiNrm) = Nrm_mc

    ! Sedimentation Velocities.
    hydromet_vel(:,iirrm) = Vrr
    hydromet_vel(:,iiNrm) = VNr


    return

  end subroutine KK_microphys_output

  !=============================================================================
  subroutine unpack_pdf_params_KK( d_variables, mu_x_1_n, mu_x_2_n, &
                                   sigma_x_1_n, sigma_x_2_n, &
                                   corr_array_1_n, &
                                   corr_array_2_n, &
                                   hydromet_pdf_params, &
                                   mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, &
                                   mu_eta_1, mu_eta_2, mu_rr_1, mu_rr_2, &
                                   mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                   mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                   mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                   sigma_w_1, sigma_w_2, sigma_chi_1, &
                                   sigma_chi_2, sigma_eta_1, sigma_eta_2, &
                                   sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                   sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                   sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                   sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                   corr_w_chi_1, corr_w_chi_2, &
                                   corr_chi_eta_1, corr_chi_eta_2, &
                                   corr_w_rr_1_n, corr_w_rr_2_n, &
                                   corr_w_Nr_1_n, corr_w_Nr_2_n, &
                                   corr_w_Ncn_1_n, corr_w_Ncn_2_n, &
                                   corr_chi_rr_1_n, corr_chi_rr_2_n, &
                                   corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                                   corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                   corr_eta_rr_1_n, corr_eta_rr_2_n, &
                                   corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                                   corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                                   corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                   rr_1, rr_2, Nr_1, Nr_2, &
                                   precip_frac, precip_frac_1, precip_frac_2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use corr_varnce_module, only: &
        iiPDF_w,        & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_rr,    &
        iiPDF_Nr,       &
        iiPDF_Ncn

    use array_index, only: &
        iirrm, & ! Variable(s)
        iiNrm

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array.

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_chi_1,      & ! Mean of chi (1st PDF component)                 [kg/kg]
      mu_chi_2,      & ! Mean of chi (2nd PDF component)                 [kg/kg]
      mu_eta_1,      & ! Mean of eta (1st PDF component)                 [kg/kg]
      mu_eta_2,      & ! Mean of eta (2nd PDF component)                 [kg/kg]
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
      sigma_chi_1,   & ! Standard deviation of chi (1st PDF component)   [kg/kg]
      sigma_chi_2,   & ! Standard deviation of chi (2nd PDF component)   [kg/kg]
      sigma_eta_1,   & ! Standard deviation of eta (1st PDF component)   [kg/kg]
      sigma_eta_2,   & ! Standard deviation of eta (2nd PDF component)   [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard deviation of ln rr (1st PDF component) ip  [-]
      sigma_rr_2_n,  & ! Standard deviation of ln rr (2nd PDF component) ip  [-]
      sigma_Nr_1_n,  & ! Standard deviation of ln Nr (1st PDF component) ip  [-]
      sigma_Nr_2_n,  & ! Standard deviation of ln Nr (2nd PDF component) ip  [-]
      sigma_Ncn_1_n, & ! Standard deviation of ln Ncn (1st PDF component)    [-]
      sigma_Ncn_2_n    ! Standard deviation of ln Ncn (2nd PDF component)    [-]

    real( kind = core_rknd ), intent(out) :: &
      corr_w_chi_1,     & ! Correlation of w and chi (1st PDF component)     [-]
      corr_w_chi_2,     & ! Correlation of w and chi (2nd PDF component)     [-]
      corr_chi_eta_1,   & ! Correlation of chi and eta (1st PDF component)   [-]
      corr_chi_eta_2      ! Correlation of chi and eta (2nd PDF component)   [-]

   real( kind = core_rknd ), intent(out) :: &
      corr_w_rr_1_n,    & ! Correlation of w and ln rr (1st PDF comp.) ip    [-]
      corr_w_rr_2_n,    & ! Correlation of w and ln rr (2nd PDF comp.) ip    [-]
      corr_w_Nr_1_n,    & ! Correlation of w and ln Nr (1st PDF comp.) ip    [-]
      corr_w_Nr_2_n,    & ! Correlation of w and ln Nr (2nd PDF comp.) ip    [-]
      corr_w_Ncn_1_n,   & ! Correlation of w and ln Ncn (1st PDF component)  [-]
      corr_w_Ncn_2_n,   & ! Correlation of w and ln Ncn (2nd PDF component)  [-]
      corr_chi_rr_1_n,  & ! Correlation of chi and ln rr (1st PDF comp.) ip  [-]
      corr_chi_rr_2_n,  & ! Correlation of chi and ln rr (2nd PDF comp.) ip  [-]
      corr_chi_Nr_1_n,  & ! Correlation of chi and ln Nr (1st PDF comp.) ip  [-]
      corr_chi_Nr_2_n,  & ! Correlation of chi and ln Nr (2nd PDF comp.) ip  [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      corr_eta_rr_1_n,  & ! Correlation of eta and ln rr (1st PDF comp.) ip  [-]
      corr_eta_rr_2_n,  & ! Correlation of eta and ln rr (2nd PDF comp.) ip  [-]
      corr_eta_Nr_1_n,  & ! Correlation of eta and ln Nr (1st PDF comp.) ip  [-]
      corr_eta_Nr_2_n,  & ! Correlation of eta and ln Nr (2nd PDF comp.) ip  [-]
      corr_eta_Ncn_1_n, & ! Correlation of eta and ln Ncn (1st PDF comp.)    [-]
      corr_eta_Ncn_2_n, & ! Correlation of eta and ln Ncn (2nd PDF comp.)    [-]
      corr_rr_Nr_1_n,   & ! Correlation of ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rr_Nr_2_n      ! Correlation of ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(out) :: &
      rr_1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr_2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr_1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr_2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]


    ! Unpack mu_x_i and sigma_x_i into Means and Standard Deviations.
    mu_w_1        = mu_x_1_n(iiPDF_w)
    mu_w_2        = mu_x_2_n(iiPDF_w)
    mu_chi_1      = mu_x_1_n(iiPDF_chi)
    mu_chi_2      = mu_x_2_n(iiPDF_chi)
    mu_eta_1      = mu_x_1_n(iiPDF_eta)
    mu_eta_2      = mu_x_2_n(iiPDF_eta)
    mu_rr_1_n     = mu_x_1_n(iiPDF_rr)
    mu_rr_2_n     = mu_x_2_n(iiPDF_rr)
    mu_Nr_1_n     = mu_x_1_n(iiPDF_Nr)
    mu_Nr_2_n     = mu_x_2_n(iiPDF_Nr)
    mu_Ncn_1_n    = mu_x_1_n(iiPDF_Ncn)
    mu_Ncn_2_n    = mu_x_2_n(iiPDF_Ncn)
    sigma_w_1     = sigma_x_1_n(iiPDF_w)
    sigma_w_2     = sigma_x_2_n(iiPDF_w)
    sigma_chi_1   = sigma_x_1_n(iiPDF_chi)
    sigma_chi_2   = sigma_x_2_n(iiPDF_chi)
    sigma_eta_1   = sigma_x_1_n(iiPDF_eta)
    sigma_eta_2   = sigma_x_2_n(iiPDF_eta)
    sigma_rr_1_n  = sigma_x_1_n(iiPDF_rr)
    sigma_rr_2_n  = sigma_x_2_n(iiPDF_rr)
    sigma_Nr_1_n  = sigma_x_1_n(iiPDF_Nr)
    sigma_Nr_2_n  = sigma_x_2_n(iiPDF_Nr)
    sigma_Ncn_1_n = sigma_x_1_n(iiPDF_Ncn)
    sigma_Ncn_2_n = sigma_x_2_n(iiPDF_Ncn)

    ! Unpack variables from hydromet_pdf_params
    mu_rr_1     = hydromet_pdf_params%mu_hm_1(iirrm)
    mu_rr_2     = hydromet_pdf_params%mu_hm_2(iirrm)
    mu_Nr_1     = hydromet_pdf_params%mu_hm_1(iiNrm)
    mu_Nr_2     = hydromet_pdf_params%mu_hm_2(iiNrm)
    mu_Ncn_1    = hydromet_pdf_params%mu_Ncn_1
    mu_Ncn_2    = hydromet_pdf_params%mu_Ncn_2
    sigma_rr_1  = hydromet_pdf_params%sigma_hm_1(iirrm)
    sigma_rr_2  = hydromet_pdf_params%sigma_hm_2(iirrm)
    sigma_Nr_1  = hydromet_pdf_params%sigma_hm_1(iiNrm)
    sigma_Nr_2  = hydromet_pdf_params%sigma_hm_2(iiNrm)
    sigma_Ncn_1 = hydromet_pdf_params%sigma_Ncn_1
    sigma_Ncn_2 = hydromet_pdf_params%sigma_Ncn_2

    rr_1          = hydromet_pdf_params%hm_1(iirrm)
    rr_2          = hydromet_pdf_params%hm_2(iirrm)
    Nr_1          = hydromet_pdf_params%hm_1(iiNrm)
    Nr_2          = hydromet_pdf_params%hm_2(iiNrm)
    precip_frac   = hydromet_pdf_params%precip_frac
    precip_frac_1 = hydromet_pdf_params%precip_frac_1
    precip_frac_2 = hydromet_pdf_params%precip_frac_2

    ! Unpack corr_array_1_n into correlations (1st PDF component).
    corr_chi_eta_1   = corr_array_1_n(iiPDF_eta, iiPDF_chi)
    corr_w_chi_1     = corr_array_1_n(iiPDF_w,iiPDF_chi)
    corr_chi_rr_1_n  = corr_array_1_n(iiPDF_rr, iiPDF_chi)
    corr_chi_Nr_1_n  = corr_array_1_n(iiPDF_Nr, iiPDF_chi)
    corr_chi_Ncn_1_n = corr_array_1_n(iiPDF_Ncn, iiPDF_chi)
    corr_eta_rr_1_n  = corr_array_1_n(iiPDF_rr, iiPDF_eta)
    corr_eta_Nr_1_n  = corr_array_1_n(iiPDF_Nr, iiPDF_eta)
    corr_eta_Ncn_1_n = corr_array_1_n(iiPDF_Ncn, iiPDF_eta)
    corr_w_rr_1_n    = corr_array_1_n(iiPDF_rr, iiPDF_w)
    corr_w_Nr_1_n    = corr_array_1_n(iiPDF_Nr, iiPDF_w)
    corr_w_Ncn_1_n   = corr_array_1_n(iiPDF_Ncn, iiPDF_w)
    corr_rr_Nr_1_n   = corr_array_1_n(iiPDF_Nr, iiPDF_rr)

    ! Unpack corr_array_2_n into correlations (2nd PDF component).
    corr_chi_eta_2   = corr_array_2_n(iiPDF_eta, iiPDF_chi)
    corr_w_chi_2     = corr_array_2_n(iiPDF_w,iiPDF_chi)
    corr_chi_rr_2_n  = corr_array_2_n(iiPDF_rr, iiPDF_chi)
    corr_chi_Nr_2_n  = corr_array_2_n(iiPDF_Nr, iiPDF_chi)
    corr_chi_Ncn_2_n = corr_array_2_n(iiPDF_Ncn, iiPDF_chi)
    corr_eta_rr_2_n  = corr_array_2_n(iiPDF_rr, iiPDF_eta)
    corr_eta_Nr_2_n  = corr_array_2_n(iiPDF_Nr, iiPDF_eta)
    corr_eta_Ncn_2_n = corr_array_2_n(iiPDF_Ncn, iiPDF_eta)
    corr_w_rr_2_n    = corr_array_2_n(iiPDF_rr, iiPDF_w)
    corr_w_Nr_2_n    = corr_array_2_n(iiPDF_Nr, iiPDF_w)
    corr_w_Ncn_2_n   = corr_array_2_n(iiPDF_Ncn, iiPDF_w)
    corr_rr_Nr_2_n   = corr_array_2_n(iiPDF_Nr, iiPDF_rr)


    return

  end subroutine unpack_pdf_params_KK

!===============================================================================

end module KK_microphys_module

! $Id$
!===============================================================================
module KK_microphys_module

  use constants_clubb, only: &
      zero  ! Constant(s)

  use clubb_precision, only: &
      core_rknd  ! Variable(s)

  implicit none

  private

  public :: KK_local_microphys_driver,   &
            KK_upscaled_microphys,       &
            KK_microphys_adjust,         &
            KK_microphys_adj_terms_type

  private :: KK_local_microphys_core, &
             KK_tendency_coefs,    &
             KK_upscaled_sedimentation, &
             KK_microphys_output

  ! This derived type stores information about adjustment to microphysics terms
  ! in KK microphysics.
  type KK_microphys_adj_terms_type

    real( kind = core_rknd ) :: &
      rrm_src_adj  = zero, & ! Total adj. to rrm source terms    [(kg/kg)/s]
      rrm_evap_adj = zero, & ! Total adj. to rrm evap terms      [(kg/kg)/s]
      Nrm_src_adj  = zero, & ! Total adj. to Nrm source terms    [(num/kg)/s]
      Nrm_evap_adj = zero    ! Total adj. to Nrm evap terms      [(num/kg)/s]

  end type KK_microphys_adj_terms_type

  contains

  !=============================================================================

  !=============================================================================
  subroutine KK_local_microphys_driver( gr, ngrdcol, dt, nzt,           & ! In
                                        hydromet_dim, hm_metadata,      & ! In
                                        l_latin_hypercube,              & ! In
                                        thlm, wm_zt, p_in_Pa, exner, rho, & ! In
                                        cloud_frac, w_std_dev, dzq, rcm, & ! In
                                        Ncm, chi, rvm, hydromet,        & ! In
                                        saturation_formula, sample_weight, & ! In
                                        stats,                          & ! InOut
                                        hydromet_mc, hydromet_vel,      & ! Out
                                        Ncm_mc, rcm_mc, rvm_mc, thlm_mc, & ! Out
                                        rrm_auto_diag, rrm_accr_diag,   & ! Out
                                        rrm_evap_diag, Nrm_auto_diag,   & ! Out
                                        Nrm_evap_diag )                   ! Out

    ! Description:
    !   Run local K-K microphysics for all model columns and update statistics
    !   with full-column fields after the microphysics calculations finish.
    !-----------------------------------------------------------------------

    use grid_class, only: grid

    use corr_varnce_module, only: &
      hm_metadata_type

    use parameters_microphys, only: &
      lh_microphys_type, &
      lh_microphys_non_interactive

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    implicit none

    type(grid), intent(in) :: gr

    integer, intent(in) :: &
      ngrdcol, &
      nzt, &
      hydromet_dim

    real( kind = core_rknd ), intent(in) :: dt

    type(hm_metadata_type), intent(in) :: hm_metadata

    logical, intent(in) :: l_latin_hypercube

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      thlm,        & ! Liquid potential temperature                       [K]
      wm_zt,       & ! Mean vertical velocity on thermodynamic levels      [m/s]
      p_in_Pa,     & ! Pressure                                           [Pa]
      exner,       & ! Exner function                                     [-]
      rho,         & ! Density on thermo. grid                            [kg/m^3]
      cloud_frac,  & ! Cloud fraction                                     [-]
      w_std_dev,   & ! Standard deviation of w (for LH interface)          [m/s]
      dzq,         & ! Thickness between thermo. levels (for LH interface) [m]
      rcm,         & ! Cloud water mixing ratio                           [kg/kg]
      Ncm,         & ! Grid mean value for cloud droplet conc.             [num/kg]
      chi,         & ! The variable 's' from Mellor                        [kg/kg]
      rvm,         & ! Mean water vapor mixing ratio (for LH interface)    [kg/kg]
      sample_weight  ! SILHS sample weight; unity for an ordinary call     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet

    integer, intent(in) :: saturation_formula

    type(stats_type), intent(inout) :: stats

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_mc, &       ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel         ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Ncm_mc, rcm_mc, rvm_mc, thlm_mc, &
      rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, &
      Nrm_auto_diag, Nrm_evap_diag

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      mvrr_diag, &
      rrm_src_adj_diag, &
      Nrm_src_adj_diag, &
      rrm_evap_adj_diag, &
      Nrm_evap_adj_diag, &
      rrm_mc_nonadj_diag

    call KK_local_microphys_core( &
           gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, l_latin_hypercube, &
           thlm, p_in_Pa, exner, rho, rcm, Ncm, chi, hydromet, &
           saturation_formula, hydromet_mc, hydromet_vel, &
           Ncm_mc, rcm_mc, rvm_mc, thlm_mc, &
           rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, &
           Nrm_auto_diag, Nrm_evap_diag, mvrr_diag, &
           rrm_src_adj_diag, Nrm_src_adj_diag, &
           rrm_evap_adj_diag, Nrm_evap_adj_diag, rrm_mc_nonadj_diag )

    !!! Output values for statistics
    call update_microphys_stat( "rrm_evap", rrm_evap_diag )
    call update_microphys_stat( "rrm_auto", rrm_auto_diag )
    call update_microphys_stat( "rrm_accr", rrm_accr_diag )
    call update_microphys_stat( "mvrr", mvrr_diag )
    call update_microphys_stat( "Nrm_evap", Nrm_evap_diag )
    call update_microphys_stat( "Nrm_auto", Nrm_auto_diag )
    call update_microphys_stat( "rrm_src_adj", rrm_src_adj_diag )
    call update_microphys_stat( "Nrm_src_adj", Nrm_src_adj_diag )
    call update_microphys_stat( "rrm_evap_adj", rrm_evap_adj_diag )
    call update_microphys_stat( "Nrm_evap_adj", Nrm_evap_adj_diag )
    call update_microphys_stat( "rrm_mc_nonadj", rrm_mc_nonadj_diag )

  contains

    subroutine update_microphys_stat( name, value )

      implicit none

      character(len=*), intent(in) :: name
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: value

      character(len=64) :: output_name

      output_name = name
      ! When SILHS is non-interactive, output sample-mean diagnostics under
      ! their SILHS names (e.g., lh_rrm_auto), leaving ordinary stats unchanged.
      if ( l_latin_hypercube .and. lh_microphys_type == lh_microphys_non_interactive ) then
        select case ( trim(name) )
        case ( "rrm_auto", "rrm_accr", "rrm_evap", "Nrm_auto", "Nrm_evap" )
          output_name = "lh_" // trim(name)
        case ( "mvrr" )
          output_name = "lh_m_vol_rad_rain"
        case ( "rrm_mc_nonadj" )
          output_name = "lh_rrm_mc_nonadj"
        case default
          return
        end select
      endif

      if ( l_latin_hypercube ) then
        call stats_update( trim(output_name), value * sample_weight, stats, &
                           sub_timestep_average=.true. )
      else
        call stats_update( trim(output_name), value, stats )
      endif

    end subroutine update_microphys_stat

  end subroutine KK_local_microphys_driver

  !=============================================================================
  subroutine KK_local_microphys_core( gr, ngrdcol, dt, nzt,             & ! In
                                      hydromet_dim, hm_metadata,        & ! In
                                      l_latin_hypercube,                & ! In
                                      thlm, p_in_Pa, exner, rho,        & ! In
                                      rcm, Ncm, chi, hydromet,          & ! In
                                      saturation_formula,               & ! In
                                      hydromet_mc, hydromet_vel,        & ! Out
                                      Ncm_mc, rcm_mc, rvm_mc, thlm_mc,  & ! Out
                                      rrm_auto_diag, rrm_accr_diag,     & ! Out
                                      rrm_evap_diag, Nrm_auto_diag,     & ! Out
                                      Nrm_evap_diag, mvrr_diag,         & ! Out
                                      rrm_src_adj_diag, Nrm_src_adj_diag, & ! Out
                                      rrm_evap_adj_diag, Nrm_evap_adj_diag, & ! Out
                                      rrm_mc_nonadj_diag )                ! Out

    ! Description:

    ! References:
    ! Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    !    Parameterization in a Large-Eddy Simulation Model of Marine
    !    Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
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

    use advance_microphys_module, only: &
        get_cloud_top_level    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd       ! Variable(s)

    use parameters_microphys, only: &
        l_silhs_KK_convergence_adj_mean ! Variable

    use corr_varnce_module, only: &
      hm_metadata_type

    use grid_class, only: &
        grid ! Type

    implicit none

    ! Input Variables
    type(grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      ngrdcol, &      ! Number of model columns
      nzt, &          ! Number of model vertical grid levels
      hydromet_dim  

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    logical, intent(in) :: &
      l_latin_hypercube    ! Flag to use Latin Hypercube interface

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      thlm,       & ! Mean liquid water potential temperature         [K]
      p_in_Pa,    & ! Pressure                                        [Pa]
      exner,      & ! Exner function                                  [-]
      rho,        & ! Density                                         [kg/m^3]
      rcm,        & ! Mean cloud water mixing ratio                   [kg/kg]
      Ncm,        & ! Mean cloud droplet conc., < N_c >               [num/kg]
      chi           ! Mean extended liquid water mixing ratio         [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet    ! Hydrometeor species                      [units vary]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Ncm_mc,  & ! Time tendency of cloud droplet concentration  [num/kg/s]
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      rrm_auto_diag, &
      rrm_accr_diag, &
      rrm_evap_diag, &
      Nrm_auto_diag, &
      Nrm_evap_diag, &
      mvrr_diag, &
      rrm_src_adj_diag, &
      Nrm_src_adj_diag, &
      rrm_evap_adj_diag, &
      Nrm_evap_adj_diag, &
      rrm_mc_nonadj_diag

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      KK_auto_tndcy,     & ! Mean KK (dr_r/dt) due to autoconversion   [(kg/kg)/s]
      KK_accr_tndcy,     & ! Mean KK (dr_r/dt) due to accretion        [(kg/kg)/s]
      KK_evap_tndcy,     & ! Mean KK (dr_r/dt) due to evaporation      [(kg/kg)/s]
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation      [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.        [(num/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      rrm,    & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,    & ! Mean rain drop concentration, < N_r >    [num/kg]
      Vrr,    & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,    & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrm_mc, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc    ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      KK_mean_vol_rad    ! Mean KK rain drop mean volume radius     [m]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

    type(KK_microphys_adj_terms_type), dimension(ngrdcol,nzt) :: &
      adj_terms    ! Adjustments to microphysics terms

    logical :: &
      l_src_adj_enabled,   & ! Flag to enable rrm/Nrm source adjustment
      l_evap_adj_enabled,  & ! Flag to enable rrm/Nrm evaporation adjustment
      l_clip_positive_sed    ! Flag to enable Vrr/VNr positive clipping

    integer :: &
      i, & ! Column loop index
      k    ! Vertical loop index

    integer, dimension(ngrdcol) :: &
      cloud_top_level ! Vertical level index of cloud top

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    ! Set up mean field variables for <r_r> and <N_r>.
    rrm = hydromet(:,:,hm_metadata%iirr)
    Nrm = hydromet(:,:,hm_metadata%iiNr)
    ! Set KK microphysics tendency adjustment flags
    l_src_adj_enabled = .true.
    l_evap_adj_enabled = .true.
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
    do k = 1, nzt, 1
      do i = 1, ngrdcol

        !!! Calculate the coefficients for the KK microphysics tendencies.
        call KK_tendency_coefs( thlm(i,k), exner(i,k), p_in_Pa(i,k), rho(i,k), &
                                KK_evap_coef, KK_auto_coef, &
                                KK_accr_coef, KK_mvr_coef, &
                                saturation_formula )

        !!! Calculate the local KK rain drop mean volume radius.
        if ( rrm(i,k) > rr_tol ) then
          KK_mean_vol_rad(i,k) = KK_mvr_local_mean( rrm(i,k), Nrm(i,k), KK_mvr_coef )
        else
          KK_mean_vol_rad(i,k) = zero
        endif

        !!! Calculate the local KK evaporation tendency.
        if ( rrm(i,k) > rr_tol .and. Nrm(i,k) > Nr_tol ) then
          KK_evap_tndcy(i,k) = KK_evap_local_mean( &
                                  chi(i,k), rrm(i,k), Nrm(i,k), KK_evap_coef )
        else
          KK_evap_tndcy(i,k) = zero
        endif

        !!! Calculate the local KK autoconversion tendency.
        if ( Ncm(i,k) > Nc_tol ) then
          KK_auto_tndcy(i,k) = KK_auto_local_mean( chi(i,k), Ncm(i,k), KK_auto_coef )
        else
          KK_auto_tndcy(i,k) = zero
        endif

        !!! Calculate the local KK accretion tendency.
        if ( rrm(i,k) > rr_tol ) then
          KK_accr_tndcy(i,k) = KK_accr_local_mean( chi(i,k), rrm(i,k), KK_accr_coef )
        else
          KK_accr_tndcy(i,k) = zero
        endif

        !!! Calculate the KK rain-number tendencies.
        if ( rrm(i,k) > rr_tol .and. Nrm(i,k) > Nr_tol ) then
          KK_Nrm_evap_tndcy(i,k) = KK_Nrm_evap_local_mean( &
                                      KK_evap_tndcy(i,k), Nrm(i,k), rrm(i,k), dt )
        else
          KK_Nrm_evap_tndcy(i,k) = zero
        endif
        KK_Nrm_auto_tndcy(i,k) = KK_Nrm_auto_mean( KK_auto_tndcy(i,k) )

        !!! Calculate any necessary adjustments to KK microphysics tendencies.
        call KK_microphys_adjust( &
          dt, exner(i,k), rcm(i,k), rrm(i,k), Nrm(i,k), &
          KK_evap_tndcy(i,k), KK_auto_tndcy(i,k), KK_accr_tndcy(i,k), &
          KK_Nrm_evap_tndcy(i,k), KK_Nrm_auto_tndcy(i,k), &
          l_src_adj_enabled, l_evap_adj_enabled, &
          rrm_mc(i,k), Nrm_mc(i,k), rvm_mc(i,k), rcm_mc(i,k), thlm_mc(i,k), &
          adj_terms(i,k) )

      enddo

    enddo  ! Microphysics tendency loop: k = 1, nzt, 1


    !!! Boundary conditions for microphysics tendencies.
    do i = 1, ngrdcol
      rrm_mc(i,nzt) = zero
      Nrm_mc(i,nzt) = zero
      KK_mean_vol_rad(i,nzt) = zero
      rvm_mc(i,nzt) = zero
      rcm_mc(i,nzt) = zero
      thlm_mc(i,nzt) = zero
    enddo

    ! Find the vertical level index of cloud top.
    cloud_top_level = get_cloud_top_level( nzt, ngrdcol, rcm, hydromet, &
                                           hydromet_dim, hm_metadata%iiri )

    !!! Microphysics sedimentation velocities.
    call KK_upscaled_sedimentation( ngrdcol, nzt, cloud_top_level, &
                                    KK_mean_vol_rad, l_clip_positive_sed, Vrr, VNr )

    !!! Output hydrometeor mean tendencies and mean sedimentation velocities
    !!! in output arrays.
    call KK_microphys_output( ngrdcol, nzt, hydromet_dim, hm_metadata, &
                              Vrr, VNr, rrm_mc, Nrm_mc, &
                              hydromet_mc, hydromet_vel )

    rrm_auto_diag = KK_auto_tndcy
    rrm_accr_diag = KK_accr_tndcy
    rrm_evap_diag = KK_evap_tndcy
    Nrm_auto_diag = KK_Nrm_auto_tndcy
    Nrm_evap_diag = KK_Nrm_evap_tndcy
    mvrr_diag = KK_mean_vol_rad
    rrm_src_adj_diag = adj_terms%rrm_src_adj
    Nrm_src_adj_diag = adj_terms%Nrm_src_adj
    rrm_evap_adj_diag = adj_terms%rrm_evap_adj
    Nrm_evap_adj_diag = adj_terms%Nrm_evap_adj
    rrm_mc_nonadj_diag = KK_auto_tndcy + KK_accr_tndcy + KK_evap_tndcy

    return

  end subroutine KK_local_microphys_core

  !=============================================================================
  subroutine KK_upscaled_microphys( gr, ngrdcol, dt, nzt, nzm,            & ! In
                                    pdf_dim, hydromet_dim, hm_metadata,   & ! In
                                    wm_zt, rtm, thlm, p_in_Pa,            & ! In
                                    exner, rho, rcm,                      & ! In
                                    pdf_params, hydromet_pdf_params,      & ! In
                                    precip_fracs,                         & ! In
                                    hydromet,                             & ! In
                                    mu_x_1_n, mu_x_2_n,                   & ! In
                                    sigma_x_1_n, sigma_x_2_n,             & ! In
                                    corr_array_1_n, corr_array_2_n,       & ! In
                                    saturation_formula,                   & ! In
                                    stats,                                & ! InOut
                                    hydromet_mc, hydromet_vel,            & ! Out
                                    rcm_mc, rvm_mc, thlm_mc,              & ! Out
                                    hydromet_vel_covar_zt_impc,           & ! Out
                                    hydromet_vel_covar_zt_expc,           & ! Out
                                    wprtp_mc, wpthlp_mc, rtp2_mc,         & ! Out
                                    thlp2_mc, rtpthlp_mc )                  ! Out

    ! Description:
    ! Version of KK microphysics scheme that is analytically upscaled by
    ! integrating over the product of the microphysics tendency and the
    ! functional form of the PDF.

    ! References:
    ! Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of a local
    !    microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol. Soc.,
    !    139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Griffin, B. M. and V. E. Larson, 2013:  Analytic upscaling of a local
    !    microphysics scheme. Part II: Simulations.  Q. J. Roy. Meteorol. Soc.,
    !    139, 670, 58--69, doi:http://dx.doi.org/10.1002/qj.1966.
    !
    ! Griffin, B. M., 2016:  Improving the Subgrid-Scale Representation of
    !    Hydrometeors and Microphysical Feedback Effects Using a Multivariate
    !    PDF.  Doctoral dissertation, University of Wisconsin -- Milwaukee,
    !    Milwaukee, WI, Paper 1144, 165 pp., URL
    !    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Griffin, B. M. and V. E. Larson, 2016:  Supplement of A new subgrid-scale
    !    representation of hydrometeor fields using a multivariate PDF.
    !    Geosci. Model Dev., 9, 6,
    !    doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Griffin, B. M. and V. E. Larson, 2016:  Parameterizing microphysical
    !    effects on variances and covariances of moisture and heat content using
    !    a multivariate probability density function: a study with CLUBB (tag
    !    MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    !    doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm_api,  & ! Procedure(s)
        grid          ! Type

    use constants_clubb, only: &
        zero

    use advance_microphys_module, only: &
        get_cloud_top_level ! Procedure(s)

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

    use KK_upscaled_variances, only: &
        variance_KK_mvr

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter, &  ! Variable(s)
        precipitation_fractions

    use clubb_precision, only: &
        core_rknd        ! Variable(s)

    use stats_netcdf, only: &
        stats_type, &
        stats_update, &
        var_on_stats_list

    use corr_varnce_module, only: &
        hm_metadata_type

    implicit none

    ! Local Constants
    logical, parameter :: &
      l_clip_positive_sed = .true.  ! Clip positive Vrr and VNr terms to zero

    !-------------------------- Input Variables --------------------------
    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: &
      dt          ! Model time step duration                 [s]

    integer, intent(in) :: &
      ngrdcol,    & ! Number of model columns
      nzt,        & ! Number of model thermodynamic vertical grid levels
      nzm,        & ! Number of model momentum vertical grid levels
      pdf_dim,    & ! Number of variables in the correlation arrays
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      wm_zt,   & ! Mean vertical velocity on thermodynamic levels  [m/s]
      rtm,     & ! Mean total water mixing ratio                   [kg/kg]
      thlm,    & ! Mean liquid water potential temperature         [K]
      p_in_Pa, & ! Pressure                                        [Pa]
      exner,   & ! Exner function                                  [-]
      rho,     & ! Density                                         [kg/m^3]
      rcm        ! Mean cloud water mixing ratio                   [kg/kg]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    type(hydromet_pdf_parameter), dimension(ngrdcol,nzt), intent(in) :: &
      hydromet_pdf_params
       
    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet       ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type(stats_type), intent(inout) :: &
      stats

    !-------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      rcm_mc,  & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc,  & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc    ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      wprtp_mc,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    !-------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rrm,    & ! Mean rain water mixing ratio, < r_r > [kg/kg]
      Nrm,    & ! Mean rain drop concentration, < N_r > [num/kg]
      Vrr,    & ! Mean sedimentation velocity of r_r    [m/s]
      VNr,    & ! Mean sedimentation velocity of N_r    [m/s]
      rrm_mc, & ! Mean (dr_r/dt) due to microphysics    [(kg/kg)/s]
      Nrm_mc    ! Mean (dN_r/dt) due to microphysics    [(num/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      KK_evap_tndcy,   & ! Mean KK (dr_r/dt) due to evaporation     [(kg/kg)/s]
      KK_auto_tndcy,   & ! Mean KK (dr_r/dt) due to autoconversion  [(kg/kg)/s]
      KK_accr_tndcy,   & ! Mean KK (dr_r/dt) due to accretion       [(kg/kg)/s]
      KK_mean_vol_rad    ! Mean KK rain drop mean volume radius     [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      KK_Nrm_evap_tndcy, & ! Mean KK (dN_r/dt) due to evaporation  [(num/kg)/s]
      KK_Nrm_auto_tndcy    ! Mean KK (dN_r/dt) due to autoconv.    [(num/kg)/s]

    real( kind = core_rknd ) :: &
      KK_evap_coef, & ! KK evaporation coefficient                  [(kg/kg)/s]
      KK_auto_coef, & ! KK autoconversion coefficient               [(kg/kg)/s]
      KK_accr_coef, & ! KK accretion coefficient                    [(kg/kg)/s]
      KK_mvr_coef     ! KK mean volume radius coefficient           [m]

     real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
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
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Vrrprrp_zt_impc, & ! Imp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)]
      Vrrprrp_zt_expc, & ! Exp. comp. of <V_rr'r_r'>: <r_r> eq.  [(m/s)(kg/kg)]
      VNrpNrp_zt_impc, & ! Imp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)]
      VNrpNrp_zt_expc    ! Exp. comp. of <V_Nr'N_r'>: <N_r> eq.  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rr_KK_mvr_covar, & ! Covariance of r_r and KK mean vol. radius [(kg/kg)m]
      Nr_KK_mvr_covar, & ! Covariance of N_r and KK mean vol. radius [(num/kg)m]
      KK_mvr_variance    ! Variance of KK rain drop mean vol. radius [m^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      w_KK_evap_covar,   & ! Covar. of w and KK evap. tend.    [m*(kg/kg)/s^2]
      rt_KK_evap_covar,  & ! Covar. of rt and KK evap. tend.   [(kg/kg)^2/s]
      thl_KK_evap_covar, & ! Covar. of thl and KK evap. tend.  [K*(kg/kg)/s]
      w_KK_auto_covar,   & ! Covar. of w and KK auto. tend.    [m*(kg/kg)/s^2]
      rt_KK_auto_covar,  & ! Covar. of rt and KK auto. tend.   [(kg/kg)^2/s]
      thl_KK_auto_covar, & ! Covar. of thl and KK auto. tend.  [K*(kg/kg)/s]
      w_KK_accr_covar,   & ! Covar. of w and KK accr. tend.    [m*(kg/kg)/s^2]
      rt_KK_accr_covar,  & ! Covar. of rt and KK accr. tend.   [(kg/kg)^2/s]
      thl_KK_accr_covar    ! Covar. of thl and KK accr. tend.  [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wprtp_mc_zt,   & ! Micro. tend. for <w'rt'>; t-lev   [m*(kg/kg)/s^2]
      wpthlp_mc_zt,  & ! Micro. tend. for <w'thl'>; t-lev  [m*K/s^2]
      rtp2_mc_zt,    & ! Micro. tend. for <rt'^2>; t-lev   [(kg/kg)^2/s]
      thlp2_mc_zt,   & ! Micro. tend. for <thl'^2>; t-lev  [K^2/s]
      rtpthlp_mc_zt    ! Micro. tend. for <rt'thl'>; t-lev [K*(kg/kg)/s]

    type(KK_microphys_adj_terms_type), dimension(ngrdcol,nzt) :: &
      adj_terms           ! Adjustments to microphysics terms

    logical :: &
      l_src_adj_enabled,  & ! Flag to enable rrm/Nrm source adjustment
      l_evap_adj_enabled    ! Flag to enable rrm/Nrm evaporation adjustment

    integer, dimension(ngrdcol) :: &
      cloud_top_level    ! Vertical level index of cloud top

    integer :: &
      i, & ! Column loop index
      k    ! Vertical loop index

    !-------------------------- Begin Code --------------------------

    rrm = hydromet(:,:,hm_metadata%iirr)
    Nrm = hydromet(:,:,hm_metadata%iiNr)
    l_src_adj_enabled = .true.
    l_evap_adj_enabled = .true.

    ! Setup mixture fraction.
    mixt_frac = pdf_params%mixt_frac

    !!! Microphysics tendency loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 1, nzt, 1
      do i = 1, ngrdcol

       !!! Calculate the coefficients for the KK microphysics tendencies.
       call KK_tendency_coefs( thlm(i,k), exner(i,k), p_in_Pa(i,k), rho(i,k), &
                               KK_evap_coef, KK_auto_coef, &
                               KK_accr_coef, KK_mvr_coef, &
                               saturation_formula )


       !!! Unpack the PDF parameters directly from the multicolumn arrays.
       ! Unpack mu_x_i and sigma_x_i into Means and Standard Deviations.
       mu_w_1 = mu_x_1_n(i,k,hm_metadata%iiPDF_w)
       mu_w_2 = mu_x_2_n(i,k,hm_metadata%iiPDF_w)
       mu_chi_1 = mu_x_1_n(i,k,hm_metadata%iiPDF_chi)
       mu_chi_2 = mu_x_2_n(i,k,hm_metadata%iiPDF_chi)
       mu_eta_1 = mu_x_1_n(i,k,hm_metadata%iiPDF_eta)
       mu_eta_2 = mu_x_2_n(i,k,hm_metadata%iiPDF_eta)
       mu_rr_1_n = mu_x_1_n(i,k,hm_metadata%iiPDF_rr)
       mu_rr_2_n = mu_x_2_n(i,k,hm_metadata%iiPDF_rr)
       mu_Nr_1_n = mu_x_1_n(i,k,hm_metadata%iiPDF_Nr)
       mu_Nr_2_n = mu_x_2_n(i,k,hm_metadata%iiPDF_Nr)
       mu_Ncn_1_n = mu_x_1_n(i,k,hm_metadata%iiPDF_Ncn)
       mu_Ncn_2_n = mu_x_2_n(i,k,hm_metadata%iiPDF_Ncn)
       sigma_w_1 = sigma_x_1_n(i,k,hm_metadata%iiPDF_w)
       sigma_w_2 = sigma_x_2_n(i,k,hm_metadata%iiPDF_w)
       sigma_chi_1 = sigma_x_1_n(i,k,hm_metadata%iiPDF_chi)
       sigma_chi_2 = sigma_x_2_n(i,k,hm_metadata%iiPDF_chi)
       sigma_eta_1 = sigma_x_1_n(i,k,hm_metadata%iiPDF_eta)
       sigma_eta_2 = sigma_x_2_n(i,k,hm_metadata%iiPDF_eta)
       sigma_rr_1_n = sigma_x_1_n(i,k,hm_metadata%iiPDF_rr)
       sigma_rr_2_n = sigma_x_2_n(i,k,hm_metadata%iiPDF_rr)
       sigma_Nr_1_n = sigma_x_1_n(i,k,hm_metadata%iiPDF_Nr)
       sigma_Nr_2_n = sigma_x_2_n(i,k,hm_metadata%iiPDF_Nr)
       sigma_Ncn_1_n = sigma_x_1_n(i,k,hm_metadata%iiPDF_Ncn)
       sigma_Ncn_2_n = sigma_x_2_n(i,k,hm_metadata%iiPDF_Ncn)

       ! Unpack variables from hydromet_pdf_params
       mu_rr_1 = hydromet_pdf_params(i,k)%mu_hm_1(hm_metadata%iirr)
       mu_rr_2 = hydromet_pdf_params(i,k)%mu_hm_2(hm_metadata%iirr)
       mu_Nr_1 = hydromet_pdf_params(i,k)%mu_hm_1(hm_metadata%iiNr)
       mu_Nr_2 = hydromet_pdf_params(i,k)%mu_hm_2(hm_metadata%iiNr)
       mu_Ncn_1 = hydromet_pdf_params(i,k)%mu_Ncn_1
       mu_Ncn_2 = hydromet_pdf_params(i,k)%mu_Ncn_2
       sigma_rr_1 = hydromet_pdf_params(i,k)%sigma_hm_1(hm_metadata%iirr)
       sigma_rr_2 = hydromet_pdf_params(i,k)%sigma_hm_2(hm_metadata%iirr)
       sigma_Nr_1 = hydromet_pdf_params(i,k)%sigma_hm_1(hm_metadata%iiNr)
       sigma_Nr_2 = hydromet_pdf_params(i,k)%sigma_hm_2(hm_metadata%iiNr)
       sigma_Ncn_1 = hydromet_pdf_params(i,k)%sigma_Ncn_1
       sigma_Ncn_2 = hydromet_pdf_params(i,k)%sigma_Ncn_2
       rr_1 = hydromet_pdf_params(i,k)%hm_1(hm_metadata%iirr)
       rr_2 = hydromet_pdf_params(i,k)%hm_2(hm_metadata%iirr)
       Nr_1 = hydromet_pdf_params(i,k)%hm_1(hm_metadata%iiNr)
       Nr_2 = hydromet_pdf_params(i,k)%hm_2(hm_metadata%iiNr)

       ! Unpack corr_array_1_n into correlations (1st PDF component).
       corr_chi_eta_1 = corr_array_1_n(i,k,hm_metadata%iiPDF_eta,hm_metadata%iiPDF_chi)
       corr_w_chi_1 = corr_array_1_n(i,k,hm_metadata%iiPDF_w,hm_metadata%iiPDF_chi)
       corr_chi_rr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_rr,hm_metadata%iiPDF_chi)
       corr_chi_Nr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_chi)
       corr_chi_Ncn_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Ncn,hm_metadata%iiPDF_chi)
       corr_eta_rr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_rr,hm_metadata%iiPDF_eta)
       corr_eta_Nr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_eta)
       corr_eta_Ncn_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Ncn,hm_metadata%iiPDF_eta)
       corr_w_rr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_rr,hm_metadata%iiPDF_w)
       corr_w_Nr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_w)
       corr_w_Ncn_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Ncn,hm_metadata%iiPDF_w)
       corr_rr_Nr_1_n = corr_array_1_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_rr)

       ! Unpack corr_array_2_n into correlations (2nd PDF component).
       corr_chi_eta_2 = corr_array_2_n(i,k,hm_metadata%iiPDF_eta,hm_metadata%iiPDF_chi)
       corr_w_chi_2 = corr_array_2_n(i,k,hm_metadata%iiPDF_w,hm_metadata%iiPDF_chi)
       corr_chi_rr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_rr,hm_metadata%iiPDF_chi)
       corr_chi_Nr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_chi)
       corr_chi_Ncn_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Ncn,hm_metadata%iiPDF_chi)
       corr_eta_rr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_rr,hm_metadata%iiPDF_eta)
       corr_eta_Nr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_eta)
       corr_eta_Ncn_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Ncn,hm_metadata%iiPDF_eta)
       corr_w_rr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_rr,hm_metadata%iiPDF_w)
       corr_w_Nr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_w)
       corr_w_Ncn_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Ncn,hm_metadata%iiPDF_w)
       corr_rr_Nr_2_n = corr_array_2_n(i,k,hm_metadata%iiPDF_Nr,hm_metadata%iiPDF_rr)
                                  
       precip_frac_1 = precip_fracs%precip_frac_1(i,k)
       precip_frac_2 = precip_fracs%precip_frac_2(i,k)

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
                                      mixt_frac(i,k), precip_frac_1, &
                                      precip_frac_2, KK_evap_coef, &
                                      KK_auto_coef, KK_accr_coef, &
                                      KK_mvr_coef, KK_evap_tndcy(i,k), &
                                      KK_auto_tndcy(i,k), KK_accr_tndcy(i,k), &
                                      KK_mean_vol_rad(i,k) )

       !!! Calculate the values of the turbulent sedimentation terms,
       !!! <V_rr'rr'> and <V_Nr'Nr'>.  Each turbulent sedimentation term has an
       !!! implicit component and an explicit component to be used in the
       !!! predictive equation for the hydrometeor.
       call KK_sed_vel_covars( rrm(i,k), rr_1, rr_2, Nrm(i,k), &
                               Nr_1, Nr_2, KK_mean_vol_rad(i,k), &
                               mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, &
                               mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, &
                               sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                               sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                               sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                               KK_mvr_coef, mixt_frac(i,k), precip_frac_1, &
                               precip_frac_2, &
                               Vrrprrp_zt_impc(i,k), Vrrprrp_zt_expc(i,k), &
                               VNrpNrp_zt_impc(i,k), VNrpNrp_zt_expc(i,k), &
                               rr_KK_mvr_covar(i,k), Nr_KK_mvr_covar(i,k) )

       if ( l_var_covar_src ) then

          call KK_upscaled_covar_driver( wm_zt(i,k), rtm(i,k), thlm(i,k), &
                                         exner(i,k), rrm(i,k), Nrm(i,k), &
                                         mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, &
                                         mu_eta_1, mu_eta_2, mu_rr_1, mu_rr_2, &
                                         mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                         mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                         mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                         sigma_w_1, sigma_w_2, sigma_chi_1, &
                                         sigma_chi_2, sigma_eta_1, sigma_eta_2,&
                                         sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                         sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                         sigma_rr_1_n, sigma_rr_2_n, &
                                         sigma_Nr_1_n, sigma_Nr_2_n, &
                                         sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                         corr_w_chi_1, corr_w_chi_2, &
                                         corr_w_rr_1_n, corr_w_rr_2_n, &
                                         corr_w_Nr_1_n, corr_w_Nr_2_n, &
                                         corr_w_Ncn_1_n, corr_w_Ncn_2_n, &
                                         corr_chi_eta_1, corr_chi_eta_2, &
                                         corr_chi_rr_1_n, corr_chi_rr_2_n, &
                                         corr_chi_Nr_1_n, corr_chi_Nr_2_n, &
                                         corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                         corr_eta_rr_1_n, corr_eta_rr_2_n, &
                                         corr_eta_Nr_1_n, corr_eta_Nr_2_n, &
                                         corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, &
                                         corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                                         mixt_frac(i,k), precip_frac_1, &
                                         precip_frac_2, &
                                         KK_evap_coef, KK_auto_coef, &
                                         KK_accr_coef, KK_evap_tndcy(i,k), &
                                         KK_auto_tndcy(i,k), KK_accr_tndcy(i,k), &
                                         pdf_params%rt_1(i,k), pdf_params%rt_2(i,k), &
                                         pdf_params%thl_1(i,k), pdf_params%thl_2(i,k), &
                                         pdf_params%crt_1(i,k), pdf_params%crt_2(i,k), &
                                         pdf_params%cthl_1(i,k), pdf_params%cthl_2(i,k), &
                                         wprtp_mc_zt(i,k), &
                                         wpthlp_mc_zt(i,k), &
                                         rtp2_mc_zt(i,k), &
                                         thlp2_mc_zt(i,k), &
                                         rtpthlp_mc_zt(i,k), &
                                         w_KK_evap_covar(i,k), rt_KK_evap_covar(i,k), &
                                         thl_KK_evap_covar(i,k), w_KK_auto_covar(i,k), &
                                         rt_KK_auto_covar(i,k), thl_KK_auto_covar(i,k), &
                                         w_KK_accr_covar(i,k), rt_KK_accr_covar(i,k), &
                                         thl_KK_accr_covar(i,k) )

       endif

       !!! KK rain drop concentration microphysics tendencies.

       !!! Calculate the KK N_r evaporation tendency.
       KK_Nrm_evap_tndcy(i,k) &
       = KK_Nrm_evap_upscaled_mean( mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, &
                                    mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, &
                                    mu_Nr_1_n, mu_Nr_2_n, sigma_chi_1, &
                                    sigma_chi_2, sigma_rr_1, sigma_rr_2, &
                                    sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, &
                                    sigma_rr_2_n, sigma_Nr_1_n, &
                                    sigma_Nr_2_n, corr_chi_rr_1_n, &
                                    corr_chi_rr_2_n, corr_chi_Nr_1_n, &
                                    corr_chi_Nr_2_n, corr_rr_Nr_1_n, &
                                    corr_rr_Nr_2_n, KK_evap_coef, mixt_frac(i,k), &
                                    precip_frac_1, precip_frac_2, dt )


       !!! Calculate the KK N_r autoconversion tendency.
       KK_Nrm_auto_tndcy(i,k) = KK_Nrm_auto_mean( KK_auto_tndcy(i,k) )


       !!! Calculate any necessary adjustments to KK microphysics tendencies.
       call KK_microphys_adjust( dt, exner(i,k), rcm(i,k), rrm(i,k), Nrm(i,k), &
                                 KK_evap_tndcy(i,k), KK_auto_tndcy(i,k), &
                                 KK_accr_tndcy(i,k), KK_Nrm_evap_tndcy(i,k), &
                                 KK_Nrm_auto_tndcy(i,k), l_src_adj_enabled, &
                                 l_evap_adj_enabled, &
                                 rrm_mc(i,k), Nrm_mc(i,k), &
                                 rvm_mc(i,k), rcm_mc(i,k), thlm_mc(i,k), &
                                 adj_terms(i,k) )


       if ( stats%l_sample .and. var_on_stats_list( stats, "KK_mvr_variance_zt" ) ) then
          ! Calculate the variance of KK rain drop mean volume radius,
          ! < R_vr'^2 >.
          KK_mvr_variance(i,k) &
          = variance_KK_mvr( mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &
                             mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, &
                             sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, &
                             sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                             sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, &
                             KK_mean_vol_rad(i,k), KK_mvr_coef, mixt_frac(i,k), &
                             precip_frac_1, precip_frac_2 )
       endif


      enddo
    enddo  ! Microphysics tendency loops

    ! Statistics
    if ( stats%l_sample ) then
      ! Covariance of r_r and KK rain drop mean volume radius.
      call stats_update( "rr_KK_mvr_covar_zt", rr_KK_mvr_covar, stats )
      ! Covariance of N_r and KK rain drop mean volume radius.
      call stats_update( "Nr_KK_mvr_covar_zt", Nr_KK_mvr_covar, stats )
      if ( l_var_covar_src ) then
        ! All of these covariance variables are being calculated on thermodynamic
        ! grid levels (all inputs are on thermodynamic grid levels, so the output
        ! is also on thermodynamic grid levels).  These covariances will be
        ! combined in various ways to produce the microphysics tendency terms for
        ! various model predictive variances and covariances.  These source
        ! tendency terms will be interpolated to momentum grid levels.
        ! Covariance of w and KK evaporation tendency.
        call stats_update( "w_KK_evap_covar_zt", w_KK_evap_covar, stats )
        ! Covariance of r_t and KK evaporation tendency.
        call stats_update( "rt_KK_evap_covar_zt", rt_KK_evap_covar, stats )
        ! Covariance of theta_l and KK evaporation tendency.
        call stats_update( "thl_KK_evap_covar_zt", thl_KK_evap_covar, stats )
        ! Covariance of w and KK autoconversion tendency.
        call stats_update( "w_KK_auto_covar_zt", w_KK_auto_covar, stats )
        ! Covariance of r_t and KK autoconversion tendency.
        call stats_update( "rt_KK_auto_covar_zt", rt_KK_auto_covar, stats )
        ! Covariance of theta_l and KK autoconversion tendency.
        call stats_update( "thl_KK_auto_covar_zt", thl_KK_auto_covar, stats )
        ! Covariance of w and KK accretion tendency.
        call stats_update( "w_KK_accr_covar_zt", w_KK_accr_covar, stats )
        ! Covariance of r_t and KK accretion tendency.
        call stats_update( "rt_KK_accr_covar_zt", rt_KK_accr_covar, stats )
        ! Covariance of theta_l and KK accretion tendency.
        call stats_update( "thl_KK_accr_covar_zt", thl_KK_accr_covar, stats )
      endif
      call stats_update( "rrm_evap", KK_evap_tndcy, stats )
      call stats_update( "rrm_auto", KK_auto_tndcy, stats )
      call stats_update( "rrm_accr", KK_accr_tndcy, stats )
      call stats_update( "mvrr", KK_mean_vol_rad, stats )
      call stats_update( "Nrm_evap", KK_Nrm_evap_tndcy, stats )
      call stats_update( "Nrm_auto", KK_Nrm_auto_tndcy, stats )
      if ( var_on_stats_list( stats, "KK_mvr_variance_zt" ) ) then
        call stats_update( "KK_mvr_variance_zt", KK_mvr_variance, stats )
      endif
    endif


    if ( l_var_covar_src ) then

       ! Output microphysics tendency terms for
       ! model variances and covariances on momentum levels.
       wprtp_mc   = zt2zm_api( nzm, nzt, ngrdcol, gr, wprtp_mc_zt )
       wpthlp_mc  = zt2zm_api( nzm, nzt, ngrdcol, gr, wpthlp_mc_zt )
       rtp2_mc    = zt2zm_api( nzm, nzt, ngrdcol, gr, rtp2_mc_zt )
       thlp2_mc   = zt2zm_api( nzm, nzt, ngrdcol, gr, thlp2_mc_zt )
       rtpthlp_mc = zt2zm_api( nzm, nzt, ngrdcol, gr, rtpthlp_mc_zt )

       ! Set values of microphysics tendency terms to 0 at model lower boundary.
       wprtp_mc(:,1)   = zero
       wpthlp_mc(:,1)  = zero
       rtp2_mc(:,1)    = zero
       thlp2_mc(:,1)   = zero
       rtpthlp_mc(:,1) = zero
       ! Set values of microphysics tendency terms to 0 at model upper boundary.
       wprtp_mc(:,nzm)   = zero
       wpthlp_mc(:,nzm)  = zero
       rtp2_mc(:,nzm)    = zero
       thlp2_mc(:,nzm)   = zero
       rtpthlp_mc(:,nzm) = zero

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
    rrm_mc(:,nzt) = zero
    Nrm_mc(:,nzt) = zero

    ! Boundary conditions
    KK_mean_vol_rad(:,nzt) = zero
    rvm_mc(:,nzt) = zero
    rcm_mc(:,nzt) = zero
    thlm_mc(:,nzt) = zero

    ! Find the vertical level index of cloud top.
    cloud_top_level = get_cloud_top_level( nzt, ngrdcol, rcm, hydromet, &
                                           hydromet_dim, hm_metadata%iiri )

    !!! Microphysics sedimentation velocities.
    call KK_upscaled_sedimentation( &
           ngrdcol, nzt, cloud_top_level, KK_mean_vol_rad, & ! In
           l_clip_positive_sed, & ! In
           Vrr, VNr ) ! Out

    !!! Output hydrometeor mean tendencies and mean sedimentation velocities
    !!! in output arrays.
    hydromet_mc = zero
    hydromet_vel = zero
    hydromet_mc(:,:,hm_metadata%iirr) = rrm_mc
    hydromet_mc(:,:,hm_metadata%iiNr) = Nrm_mc
    hydromet_vel(:,:,hm_metadata%iirr) = Vrr
    hydromet_vel(:,:,hm_metadata%iiNr) = VNr

    !!! Turbulent sedimentation above cloud top should have a value of 0.
    do k = 1, nzt-1
      do i = 1, ngrdcol
        if ( cloud_top_level(i) > 1 .and. k > cloud_top_level(i) ) then
          Vrrprrp_zt_impc(i,k) = zero
          Vrrprrp_zt_expc(i,k) = zero
          VNrpNrp_zt_impc(i,k) = zero
          VNrpNrp_zt_expc(i,k) = zero
        endif
      enddo
    enddo

    !!! Boundary conditions (upper) for the covariances of hydrometeor
    !!! sedimentation velocities and their associated hydrometeors
    !!! (<V_rr'r_r'> and <V_Nr'N_r'>).
    Vrrprrp_zt_impc(:,nzt) = zero
    Vrrprrp_zt_expc(:,nzt) = zero
    VNrpNrp_zt_impc(:,nzt) = zero
    VNrpNrp_zt_expc(:,nzt) = zero

    ! The implicit and explicit components used to calculate the covariances of
    ! hydrometeor sedimentation velocities and their associated hydrometeors
    ! (<V_rr'r_r'> and <V_Nr'N_r'>) are fed into the output arrays.
    hydromet_vel_covar_zt_impc = zero
    hydromet_vel_covar_zt_expc = zero
    hydromet_vel_covar_zt_impc(:,:,hm_metadata%iirr) = Vrrprrp_zt_impc
    hydromet_vel_covar_zt_expc(:,:,hm_metadata%iirr) = Vrrprrp_zt_expc
    hydromet_vel_covar_zt_impc(:,:,hm_metadata%iiNr) = VNrpNrp_zt_impc
    hydromet_vel_covar_zt_expc(:,:,hm_metadata%iiNr) = VNrpNrp_zt_expc

    if ( stats%l_sample ) then
      !!! Statistical output for upscaled KK.
      ! The covariance < V_rr'r_r' > calculated completely explicitly.
      ! When semi-implicit turbulent advection is used, this result can be
      ! compared to the < V_rr'r_r' > results used in the code, which are
      ! calculated semi-implicitly.
      call stats_update( "Vrrprrp_expcalc", zt2zm_api( nzm, nzt, ngrdcol, gr, &
                         Vrrprrp_zt_impc * rrm + Vrrprrp_zt_expc ), stats )
      ! The covariance < V_Nr'N_r' > calculated completely explicitly.
      ! When semi-implicit turbulent advection is used, this result can be
      ! compared to the < V_Nr'N_r' > results used in the code, which are
      ! calculated semi-implicitly.
      call stats_update( "VNrpNrp_expcalc", zt2zm_api( nzm, nzt, ngrdcol, gr, &
                         VNrpNrp_zt_impc * Nrm + VNrpNrp_zt_expc ), stats )
      !!! Statistical output for mean microphysics tendencies.
      ! Stats output for microphysics adjustment terms.
      call stats_update( "rrm_src_adj", adj_terms%rrm_src_adj, stats )
      call stats_update( "Nrm_src_adj", adj_terms%Nrm_src_adj, stats )
      call stats_update( "rrm_evap_adj", adj_terms%rrm_evap_adj, stats )
      call stats_update( "Nrm_evap_adj", adj_terms%Nrm_evap_adj, stats )
      call stats_update( "rrm_mc_nonadj", &
                         KK_auto_tndcy + KK_accr_tndcy + KK_evap_tndcy, stats )
    end if

    return

  end subroutine KK_upscaled_microphys

  !=============================================================================
  subroutine KK_tendency_coefs( thlm, exner, p_in_Pa, rho, &
                                KK_evap_coef, KK_auto_coef, &
                                KK_accr_coef, KK_mvr_coef, &
                                saturation_formula )

    ! Description:

    ! References:
    ! Eq. (3), Eq. (22), Eq. (29), and Eq. (33) of Khairoutdinov, M. and
    ! Y. Kogan, 2000:  A New Cloud Physics Parameterization in a Large-Eddy
    ! Simulation Model of Marine Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
    !
    ! Eq. (22), Eq. (28), Eq. (38), and Eq. (51) of Larson, V. E. and
    ! B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
    ! Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
    ! doi:http://dx.doi.org/10.1002/qj.1967.
    !
    ! Eq. (C21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !
    ! Eq. (S21) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    ! A new subgrid-scale representation of hydrometeor fields using a
    ! multivariate PDF.  Geosci. Model Dev., 9, 6,
    ! doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    !
    ! Eq. (A27) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    ! microphysical effects on variances and covariances of moisture and heat
    ! content using a multivariate probability density function: a study with
    ! CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295, 
    ! doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use saturation, only: &
        sat_mixrat_liq_api  ! Procedure(s)

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

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

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
    r_sl = sat_mixrat_liq_api( p_in_Pa, T_liq_in_K, saturation_formula )

    ! Beta(T_l).
    Beta_Tl = (Rd/Rv) * ( Lv / ( Rd * T_liq_in_K ) )  &
                      * ( Lv / ( Cp * T_liq_in_K ) )

    ! Coefficient for KK evaporation.
    KK_evap_coef = three * C_evap * G_T_p( T_liq_in_K, p_in_Pa, saturation_formula )   &
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
        core_rknd  ! Variable(s)

    use constants_clubb, only: &
        Lv,           & ! Constant(s)
        Cp,           &
        zero,         &
        rr_tol,       &
        Nr_tol,       &
        eps

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
       if ( abs(KK_evap_tndcy-rrm_evap_net) > abs(KK_evap_tndcy+rrm_evap_net)*eps/2 &
            .and. rrm > rr_tol &
            .and. Nrm > Nr_tol ) then
          Nrm_evap_net = KK_Nrm_evap_local_mean( rrm_evap_net, Nrm, rrm, dt )
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
    adj_terms%rrm_evap_adj = rrm_evap_net - KK_evap_tndcy
    adj_terms%Nrm_evap_adj = Nrm_evap_net - KK_Nrm_evap_tndcy

    return

  end subroutine KK_microphys_adjust

  !=============================================================================
  subroutine KK_upscaled_sedimentation( ngrdcol, nzt, cloud_top_level, &
                                        KK_mean_vol_rad, &
                                        l_clip_positive_sed, Vrr, VNr )

    ! Description:
    ! Calculate sedimentation velocities for every upscaled K-K model column.

    ! References:
    ! Eq. (37) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    ! Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    ! Mon. Wea. Rev., 128, 229--243.
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
      ngrdcol, & ! Number of model columns
      nzt        ! Number of model thermodynamic vertical grid levels

    integer, dimension(ngrdcol), intent(in) :: &
      cloud_top_level     ! Vertical level index of cloud top

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      KK_mean_vol_rad     ! KK rain drop mean volume radius       [m]

    logical, intent(in) :: &
      l_clip_positive_sed ! Clip positive values of Vrr and VNr   [T/F]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      Vrr, & ! Mean sedimentation velocity of < r_r >             [m/s]
      VNr    ! Mean sedimentation velocity of < N_r >             [m/s]

    ! Local Variables
    integer :: &
      i, & ! Column loop iterator
      k    ! Vertical loop iterator


    !!! Mean sedimentation velocities
    do k = 1, nzt-1
      do i = 1, ngrdcol

       ! Mean sedimentation velocity of rain water mixing ratio.
       Vrr(i,k) = - ( 0.012_core_rknd * ( micron_per_m * KK_mean_vol_rad(i,k) ) &
                      - 0.2_core_rknd )

       ! Mean sedimentation velocity of rain drop concentration.
       VNr(i,k) = - ( 0.007_core_rknd * ( micron_per_m * KK_mean_vol_rad(i,k) ) &
                      - 0.1_core_rknd )

       if ( l_clip_positive_sed ) then

         ! Mean sedimentation velocity of rain water mixing ratio cannot
         ! have a
         ! positive value.
         if ( Vrr(i,k) > zero ) then
            Vrr(i,k) = zero
         endif

         ! Mean sedimentation velocity of rain drop concentration cannot have a
         ! positive value.
         if ( VNr(i,k) > zero ) then
            VNr(i,k) = zero
         endif

         if ( cloud_top_level(i) > 1 .and. k > cloud_top_level(i) ) then
            Vrr(i,k) = zero
            VNr(i,k) = zero
         endif

       end if

      enddo
    enddo ! Sedimentation velocity loop

    !!! Boundary conditions for sedimentation velocities.

    ! The flux of rain water through the model top is 0.
    ! Vrr and VNr are set to 0 at the highest model level.
    Vrr(:,nzt) = zero
    VNr(:,nzt) = zero


    return

  end subroutine KK_upscaled_sedimentation

  !=============================================================================
  subroutine KK_microphys_output( ngrdcol, nzt, hydromet_dim, hm_metadata, &
                                  Vrr, VNr, rrm_mc, Nrm_mc, &
                                  hydromet_mc, hydromet_vel )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use corr_varnce_module, only: &
      hm_metadata_type

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol, &      ! Number of model columns
      nzt, &          ! Number of model thermodynamic vertical grid levels
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      Vrr,    & ! Mean sedimentation velocity of < r_r >   [m/s]
      VNr,    & ! Mean sedimentation velocity of < N_r >   [m/s]
      rrm_mc, & ! Mean (dr_r/dt) due to microphysics       [(kg/kg)/s]
      Nrm_mc    ! Mean (dN_r/dt) due to microphysics       [(num/kg)/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]


    !!! Output mean hydrometeor tendencies and mean sedimentation velocities.

    hydromet_mc = zero
    hydromet_vel = zero

    ! Mean field tendencies.
    hydromet_mc(:,:,hm_metadata%iirr) = rrm_mc
    hydromet_mc(:,:,hm_metadata%iiNr) = Nrm_mc

    ! Sedimentation Velocities.
    hydromet_vel(:,:,hm_metadata%iirr) = Vrr
    hydromet_vel(:,:,hm_metadata%iiNr) = VNr


    return

  end subroutine KK_microphys_output

!===============================================================================

end module KK_microphys_module

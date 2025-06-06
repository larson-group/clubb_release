!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module estimate_scm_microphys_module

  implicit none

  public :: est_single_column_tndcy

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine est_single_column_tndcy( &
               gr, dt, nzt, nzm, num_samples, &
               pdf_dim, hydromet_dim, hm_metadata, &
               X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights, &
               pdf_params, precip_fracs, p_in_Pa, exner, rho, &
               dzq, hydromet, rcm, &
               lh_rt_clipped, lh_thl_clipped, &
               lh_rc_clipped, lh_rv_clipped, &
               lh_Nc_clipped, &
               l_lh_instant_var_covar_src, &
               saturation_formula, &
               stats_metadata, &
               stats_zt, stats_zm, stats_sfc, stats_lh_zt, &
               lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
               lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc, &
               microphys_sub )
! Description:
!   Estimate the tendency of a microphysics scheme via latin hypercube sampling
!
! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      zero, & ! Constant(s)
      unused_var

    use grid_class, only: grid ! Type


    use parameters_microphys, only: &
      l_var_covar_src

    use clubb_precision, only: &
      core_rknd

    use grid_class, only: &
      zt2zm_api

    use stats_variables, only: &
        stats_metadata_type

    use stats_type_utilities, only: &
      stat_update_var ! Procedure

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, &
      microphys_stats_accumulate, &
      microphys_stats_cleanup

    use math_utilities, only: &
      compute_sample_mean             ! Procedure

    use parameters_microphys, only: &
      lh_microphys_type,            &           ! Variable
      lh_microphys_non_interactive              ! Constant

    use latin_hypercube_driver_module, only: &
      copy_X_nl_into_hydromet_all_pts    ! Procedure(s)


    use parameters_microphys, only: &
      l_silhs_KK_convergence_adj_mean ! Variable(s)

    use corr_varnce_module, only: &
      hm_metadata_type

    use silhs_api_module, only: &
      lh_microphys_var_covar_driver_api   ! Procedure

    use silhs_category_variance_module, only: &
      silhs_category_variance_driver  ! Procedure

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions

    use stats_type, only: stats ! Type

    implicit none

    ! External
#include "microphys_interface.inc"

    intrinsic :: real, all, any

    ! Constant parameters
    logical, parameter :: &
      l_latin_hypercube = .true. ! We are the Latin hypercube!

    ! Input Variables
    type (grid), intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      nzt,           & ! Number of thermodynamic vertical levels
      nzm,           & ! Number of momentum vertical levels
      num_samples,   & ! Number of calls to microphysics
      pdf_dim,       & ! Number of variates
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs    ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component of each sample

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight for cloud weighted sampling

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! The PDF parameters

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho,        & ! Density on thermo. grid  [kg/m^3]
      dzq,        & ! Difference in height per gridbox   [m]
      rcm           ! Mean liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    logical, intent(in) :: &
      l_lh_instant_var_covar_src ! Produce instantaneous var/covar tendencies [-]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! InOut Variables
    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt

    ! Output Variables

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nzt), intent(out) :: &
      lh_Ncm_mc,     & ! LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
      lh_rcm_mc,     & ! LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
      lh_rvm_mc,     & ! LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
      lh_thlm_mc       ! LH estimate of time tndcy. of liquid potential temp.  [K/s]

    real( kind = core_rknd ), dimension(nzm), intent(out) :: &
      lh_rtp2_mc,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]


    ! Local Variables
    real( kind = core_rknd ), dimension(num_samples,nzt,hydromet_dim) :: &
      lh_hydromet_mc_all, & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_all   ! LH est of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(num_samples,nzt) :: &
      lh_Ncm_mc_all,       & ! LH est of time tendency of cloud droplet concentration  [#/kg/s]
      lh_rcm_mc_all,       & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_all,       & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_all         ! LH est of time tendency of liquid potential temperature     [K/s]

    real( kind = core_rknd ), dimension(num_samples,nzt,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species                    [units vary]

    real( kind = core_rknd ), dimension(num_samples,nzt) :: &
      Ncn_all_points, &    ! Cloud Nuclei conc. (simplified); Nc=Ncn*H(chi) [#/kg]
      chi_all_points, &    ! 's' (Mellor 1977)                              [kg/kg]
      w_all_points         ! Vertical velocity                              [m/s]

    real( kind = core_rknd ), dimension(nzt) :: &
      lh_rtp2_mc_zt,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc_zt,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc_zt,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc_zt,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc_zt    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]

    ! These parameters are not used by the microphysics scheme when SILHS is
    ! turned on.
    real( kind = core_rknd ), dimension(nzt) :: &
      cloud_frac_unused, &
      w_std_dev_unused

    type(microphys_stats_vars_type), dimension(num_samples) :: &
      microphys_stats_zt_all, &   ! Statistics variables output from microphysics on zt grid
      microphys_stats_sfc_all     ! Statistics variables on sfc grid

    type(microphys_stats_vars_type) :: &
      microphys_stats_zt_avg, &   ! Average of zt statistics over all sample points
      microphys_stats_sfc_avg     ! Average of sfc statistics over all sample points

    integer :: ivar, sample

    ! Just to avoid typing hm_metadata%iiPDF_x everywhere
    integer :: &
      iiNr, &
      iirr, &
      iiPDF_chi, &
      iiPDF_w

    ! ---- Begin Code ----

    iiNr      = hm_metadata%iiNr
    iirr      = hm_metadata%iirr
    iiPDF_chi = hm_metadata%iiPDF_chi
    iiPDF_w   = hm_metadata%iiPDF_w

    w_all_points   = real( X_nl_all_levs(:,:,iiPDF_w), kind=core_rknd )
    chi_all_points = real( X_nl_all_levs(:,:,iiPDF_chi), kind=core_rknd )

    call copy_X_nl_into_hydromet_all_pts &
         ( nzt, pdf_dim, num_samples, & ! Intent(in)
           X_nl_all_levs,                & ! Intent(in)
           hydromet_dim, hm_metadata, & ! Intent(in)
           hydromet,                     & ! Intent(in)
           hydromet_all_points,          & ! Intent(out)
           Ncn_all_points )                ! Intent(out)

    do sample = 1, num_samples

      cloud_frac_unused = unused_var
      w_std_dev_unused  = unused_var
      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub( &
             gr, dt, nzt, & ! In
             hydromet_dim, hm_metadata, & ! In
             l_latin_hypercube, lh_thl_clipped(sample,:), w_all_points(sample,:), p_in_Pa, & ! In
             exner, rho, cloud_frac_unused, w_std_dev_unused, & ! In
             dzq, lh_rc_clipped(sample,:), lh_Nc_clipped(sample,:), & ! In
             chi_all_points(sample,:), lh_rv_clipped(sample,:), & ! In
             hydromet_all_points(sample,:,:), & ! In
             saturation_formula, &
             stats_metadata, & ! In
             lh_hydromet_mc_all(sample,:,:), lh_hydromet_vel_all(sample,:,:), & ! Out
             lh_Ncm_mc_all(sample,:), & ! Out
             lh_rcm_mc_all(sample,:), lh_rvm_mc_all(sample,:), lh_thlm_mc_all(sample,:), & ! Out
             microphys_stats_zt_all(sample), microphys_stats_sfc_all(sample) ) ! Out

      ! Loop to get new sample
    end do ! sample = 1, num_samples

    if ( l_var_covar_src ) then

      call lh_microphys_var_covar_driver_api &
           ( nzt, num_samples, dt, lh_sample_point_weights,  &  ! Intent(in)
             pdf_params, lh_rt_clipped, lh_thl_clipped, w_all_points,  &  ! Intent(in)
             lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all,  &  ! Intent(in)
             l_lh_instant_var_covar_src, &                     ! Intent(in)
             lh_rtp2_mc_zt, lh_thlp2_mc_zt, lh_wprtp_mc_zt, &  ! Intent(out)
             lh_wpthlp_mc_zt, lh_rtpthlp_mc_zt )               ! Intent(out)

      ! Convert from the zt grid to the zm grid.
      lh_rtp2_mc    = zt2zm_api( gr, lh_rtp2_mc_zt )
      lh_thlp2_mc   = zt2zm_api( gr, lh_thlp2_mc_zt )
      lh_wprtp_mc   = zt2zm_api( gr, lh_wprtp_mc_zt )
      lh_wpthlp_mc  = zt2zm_api( gr, lh_wpthlp_mc_zt )
      lh_rtpthlp_mc = zt2zm_api( gr, lh_rtpthlp_mc_zt )

      ! Stats sampling
      if ( stats_metadata%l_stats_samp ) then
        call stat_update_var( stats_metadata%ilh_rtp2_mc, lh_rtp2_mc, stats_zm )
        call stat_update_var( stats_metadata%ilh_thlp2_mc, lh_thlp2_mc, stats_zm )
        call stat_update_var( stats_metadata%ilh_wprtp_mc, lh_wprtp_mc, stats_zm )
        call stat_update_var( stats_metadata%ilh_wpthlp_mc, lh_wpthlp_mc, stats_zm )
        call stat_update_var( stats_metadata%ilh_rtpthlp_mc, lh_rtpthlp_mc, stats_zm )
      end if

    else ! .not. l_var_covar_src

      lh_rtp2_mc     = zero
      lh_thlp2_mc    = zero
      lh_wprtp_mc    = zero
      lh_wpthlp_mc   = zero
      lh_rtpthlp_mc  = zero

    end if

    ! Grid box average.
    forall( ivar = 1:hydromet_dim )

      lh_hydromet_vel(:,ivar) &
        = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                               lh_hydromet_vel_all(:,:,ivar) )

      lh_hydromet_mc(:,ivar) &
        = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                               lh_hydromet_mc_all(:,:,ivar) )

    end forall

    lh_Ncm_mc = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, lh_Ncm_mc_all(:,:) )
    lh_rcm_mc = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, lh_rcm_mc_all(:,:) )
    lh_rvm_mc = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, lh_rvm_mc_all(:,:) )
    lh_thlm_mc= compute_sample_mean( nzt, num_samples, lh_sample_point_weights, lh_thlm_mc_all(:,:))

    ! Sample variables from microphys_stats_vars objects for statistics
    microphys_stats_zt_avg =  silhs_microphys_stats_avg( num_samples, microphys_stats_zt_all, &
                                                         lh_sample_point_weights )
    microphys_stats_sfc_avg = silhs_microphys_stats_avg( num_samples, microphys_stats_sfc_all, &
                                                         lh_sample_point_weights )

    if ( lh_microphys_type /= lh_microphys_non_interactive ) then
      call microphys_stats_accumulate( microphys_stats_zt_avg, stats_metadata, stats_zt )
      call microphys_stats_accumulate( microphys_stats_sfc_avg, stats_metadata, stats_sfc )
    else
      call silhs_noninteractive_stats( microphys_stats_zt_avg, &
                                       stats_metadata, &
                                       stats_lh_zt )
    end if

    ! Adjust the mean if l_silhs_KK_convergence_adj_mean is true
    if ( l_silhs_KK_convergence_adj_mean ) then
      call adjust_KK_src_means( dt, nzt, exner, rcm, hydromet(:,iirr),          & ! intent(in)
                                hydromet(:,iiNr), hydromet,                     & ! intent(in)
                                hydromet_dim, hm_metadata%iiri,                 & ! intent(in)
                                microphys_stats_zt_avg,                         & ! intent(in)
                                stats_metadata,                                 & ! intent(in)
                                stats_lh_zt,                                    & ! intent(inout)
                                lh_hydromet_vel(:,iirr),                        & ! intent(inout)
                                lh_hydromet_vel(:,iiNr),                        & ! intent(inout)
                                lh_hydromet_mc(:,iirr), lh_hydromet_mc(:,iiNr), & ! intent(out)
                                lh_rvm_mc, lh_rcm_mc, lh_thlm_mc )                ! intent(out)
    end if

    ! Invoke the SILHS category variance sampler (if desired by user)!!

    if ( stats_metadata%l_stats_samp ) then

      if ( allocated( stats_metadata%isilhs_variance_category ) ) then

        if( stats_metadata%isilhs_variance_category(1) > 0 ) then
          call silhs_category_variance_driver( &
                 nzt, num_samples, pdf_dim, hydromet_dim, hm_metadata,      & ! Intent(in)
                 X_nl_all_levs,                                             & ! Intent(in)
                 X_mixt_comp_all_levs, microphys_stats_zt_all,              & ! Intent(in)
                 lh_hydromet_mc_all, lh_sample_point_weights, pdf_params,   & ! Intent(in)
                 precip_fracs,                                              & ! intent(in)
                 stats_metadata,                                            & ! intent(in)
                 stats_lh_zt )                                                ! intent(inout)
        end if ! stats_metadata%isilhs_variance_category(1) > 0

      end if ! allocated( stats_metadata%isilhs_variance_category ) 

    end if ! stats_metadata%l_stats_samp

    ! Cleanup microphys_stats_vars objects
    do ivar=1, num_samples
      call microphys_stats_cleanup( microphys_stats_zt_all(ivar) )
      call microphys_stats_cleanup( microphys_stats_sfc_all(ivar) )
    end do
    call microphys_stats_cleanup( microphys_stats_zt_avg )
    call microphys_stats_cleanup( microphys_stats_sfc_avg )

    return
  end subroutine est_single_column_tndcy

  !-----------------------------------------------------------------------
  function silhs_microphys_stats_avg &
           ( num_samples, microphys_stats_all, lh_sample_point_weights ) &
           result( microphys_stats_avg )

  ! Description:
  !   Computes a weighted average of all statistical variables 

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd            ! Compile-time constant

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, & ! Type
      microphys_put_var,         & ! Procedure(s)
      microphys_stats_alloc

    use math_utilities, only: &
      compute_sample_mean

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples           ! Number of SILHS sample points

    type(microphys_stats_vars_type), dimension(num_samples), intent(in) :: &
      microphys_stats_all   ! Statistics variables

    real( kind = core_rknd ), dimension(num_samples,microphys_stats_all(1)%nz), intent(in) :: &
      lh_sample_point_weights ! Weight of each SILHS sample point

    ! Output Variables
    type(microphys_stats_vars_type) :: &
      microphys_stats_avg   ! Average of statistical variables

    ! Local Variables
    integer :: ivar, svar, nz, num_vars, stats_index

    real( kind = core_rknd ), dimension(num_samples,microphys_stats_all(1)%nz) :: &
      stats_var_all

    real( kind = core_rknd ), dimension(microphys_stats_all(1)%nz) :: &
      stats_var_avg

  !-----------------------------------------------------------------------
    !----- Begin Code-----

    nz = microphys_stats_all(1)%nz
    num_vars = microphys_stats_all(1)%num_vars

    call microphys_stats_alloc( nz, num_vars, microphys_stats_avg )

    do ivar=1, num_vars

      stats_index = microphys_stats_all(1)%stats_indices(ivar)

      do svar=1, num_samples
        stats_var_all(svar,:) = microphys_stats_all(svar)%output_values(:,ivar)
      end do

      stats_var_avg = compute_sample_mean( nz, num_samples, lh_sample_point_weights, stats_var_all )

      call microphys_put_var( stats_index, stats_var_avg, microphys_stats_avg )

    end do

    return
  end function silhs_microphys_stats_avg
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine silhs_noninteractive_stats( microphys_stats_zt_avg, &
                                         stats_metadata, &
                                         stats_lh_zt )

  ! Description:
  !   When SILHS is run in non-interactive mode in CLUBB standalone, this
  !   subroutine outputs some averages from microphysics to the corresponding
  !   SILHS output variables (e.g., lh_rrm_auto)

  ! References:
  !   none
  !-----------------------------------------------------------------------

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, &  ! Type
      microphys_get_var ! Procedure

    use stats_type_utilities, only: &
      stat_update_var  ! Procedure

    use parameters_microphys, only: &
      microphys_scheme    ! Variable

    use stats_type, only: &
      stats ! Type

    use stats_variables, only: &
      stats_metadata_type

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_zt_avg  ! Statistics structure from the microphysics scheme

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !--------------------------- InOut Variables ---------------------------
    type(stats), intent(inout) :: &
      stats_lh_zt

    !--------------------------- Begin Code ---------------------------

    ! Statistical sampling
    if ( stats_metadata%l_stats_samp ) then

      if ( stats_metadata%ilh_rrm_auto > 0 ) then
        call stat_update_var( stats_metadata%ilh_rrm_auto, microphys_get_var( &
             stats_metadata%irrm_auto, microphys_stats_zt_avg ), stats_lh_zt )
      end if

      if ( stats_metadata%ilh_rrm_accr > 0 ) then
        call stat_update_var( stats_metadata%ilh_rrm_accr, microphys_get_var( &
             stats_metadata%irrm_accr, microphys_stats_zt_avg ), stats_lh_zt )
      end if

      if ( stats_metadata%ilh_rrm_evap > 0 ) then
        call stat_update_var( stats_metadata%ilh_rrm_evap, microphys_get_var( &
             stats_metadata%irrm_evap, microphys_stats_zt_avg ), stats_lh_zt )
      end if

      if ( stats_metadata%ilh_Nrm_auto > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nrm_auto, microphys_get_var( &
             stats_metadata%iNrm_auto, microphys_stats_zt_avg ), stats_lh_zt )
      end if

      if ( stats_metadata%ilh_Nrm_evap > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nrm_evap, microphys_get_var( &
             stats_metadata%iNrm_evap, microphys_stats_zt_avg ), stats_lh_zt )
      end if

      if ( trim( microphys_scheme ) == "khairoutdinov_kogan" ) then
        ! These variables are output only from KK microphysics.
        if ( stats_metadata%ilh_m_vol_rad_rain > 0 ) then
          call stat_update_var( stats_metadata%ilh_m_vol_rad_rain, microphys_get_var( &
               stats_metadata%im_vol_rad_rain, microphys_stats_zt_avg ), stats_lh_zt )
        end if

        if ( stats_metadata%ilh_rrm_mc_nonadj > 0 ) then
          call stat_update_var( stats_metadata%ilh_rrm_mc_nonadj, microphys_get_var( &
               stats_metadata%irrm_mc_nonadj, microphys_stats_zt_avg ), stats_lh_zt )
        end if
      end if ! trim( microphys_scheme ) == "khairoutdinov_kogan"

    end if ! stats_metadata%l_stats_samp

    return
  end subroutine silhs_noninteractive_stats
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine adjust_KK_src_means( dt, nzt, exner, rcm, rrm, Nrm, hydromet, &
                                  hydromet_dim, iiri,                      &
                                  microphys_stats_zt,                      &
                                  stats_metadata,                          &
                                  stats_lh_zt,                             &
                                  lh_Vrr, lh_VNr,                          &
                                  rrm_mc, Nrm_mc,                          &
                                  rvm_mc, rcm_mc, thlm_mc )

  ! Description:
  !   Adjusts the means of microphysics terms for KK microphysics by calling the
  !   KK microphysics adjustment subroutine

  ! References:
  !   clubb:ticket:558
  !-----------------------------------------------------------------------------

    use KK_Nrm_tendencies, only: &
        KK_Nrm_auto_mean, & ! Procedure(s)
        KK_Nrm_evap_local_mean

    use KK_microphys_module, only: &
        KK_microphys_adjust, &      ! Procedure
        KK_microphys_adj_terms_type ! Type

    use advance_microphys_module, only: &
        get_cloud_top_level    ! Procedure

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        zero    ! Constant

    use stats_type_utilities, only: &
        stat_update_var ! Procedure

    use microphys_stats_vars_module, only: &
        microphys_stats_vars_type, &     ! Type
        microphys_get_var                ! Procedure

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    implicit none

    type(stats), intent(inout) :: &
      stats_lh_zt

    ! Local Constants
    logical, parameter :: &
      ! Whether to adjust rrm_source to not over-deplete cloud water
      l_src_adj_enabled = .true.,  &
      ! Whether to adjust rrm_evap to not over-evaporate rain
      l_evap_adj_enabled = .true.

    !-------------------------- Input variables --------------------------
    real( kind = core_rknd ), intent(in) :: &
      dt   ! Model timestep

    integer, intent(in) :: &
      nzt, &   ! Number of thermodynamic vertical grid levels
      hydromet_dim, &
      iiri

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      exner, & ! Exner function                            [-]
      rcm,   & ! Mean liquid water mixing ratio            [kg/kg]
      rrm,   & ! Rain water mixing ration                  [kg/kg]
      Nrm      ! Rain drop concentration                   [num/kg]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet    ! Mean value of hydrometeor              [units vary]

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_zt     ! Statistics variables        [units vary]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !-------------------------- InOut Variables --------------------------
    real( kind = core_rknd ), dimension(nzt), intent(inout) :: &
      lh_Vrr, &         ! Mean sedimentation velocity of < r_r > [m/s]
      lh_VNr            ! Mean sedimentation velocity of < N_r > [m/s]

    !-------------------------- Output variables --------------------------
    real( kind = core_rknd ), dimension(nzt), intent(out) :: &
      rrm_mc, & ! Mean change in rain due to microphysics [(kg/kg)/s] 
      Nrm_mc,    & ! Mean change in Nrm due to microphysics  [(kg/kg)/s]
      rvm_mc,    & ! Time tendency of rvm                    [(kg/kg)/s]
      rcm_mc,    & ! Time tendency of rcm                    [(kg/kg)/s]
      thlm_mc      ! Time tendency of thlm                   [(kg/kg)/s]

    !-------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(nzt) :: &
      rrm_evap, & ! Mean change in rain due to evap           [(kg/kg)/s]
      rrm_auto, & ! Mean change in rain due to autoconversion [(kg/kg)/s]
      rrm_accr, & ! Mean change in rain due to accretion      [(kg/kg)/s]
      Nrm_auto,    & ! Mean change in Nrm due to autoconversion  [(num/kg)/s]
      Nrm_evap       ! Mean change in Nrm due to evaporation     [(num/kg)/s]

    type(KK_microphys_adj_terms_type), dimension(nzt) :: &
      adj_terms    ! Adjustment terms returned from the adjustment routine

    integer :: k, cloud_top_level

    !-------------------------- Begin code --------------------------

    ! Initialize output
    rrm_mc = zero
    Nrm_mc = zero
    rvm_mc = zero
    rcm_mc = zero
    thlm_mc = zero

    rrm_auto = microphys_get_var( stats_metadata%irrm_auto, microphys_stats_zt )
    rrm_accr = microphys_get_var( stats_metadata%irrm_accr, microphys_stats_zt )
    rrm_evap = microphys_get_var( stats_metadata%irrm_evap, microphys_stats_zt )

    Nrm_auto    = microphys_get_var( stats_metadata%iNrm_auto,    microphys_stats_zt )
    Nrm_evap    = microphys_get_var( stats_metadata%iNrm_evap,    microphys_stats_zt )

    ! Loop over each vertical level above the lower boundary
    do k = 1, nzt, 1

      ! We call KK_microphys_adjust to adjust the means of the mc terms
      call KK_microphys_adjust( dt, exner(k), rcm(k), rrm(k), Nrm(k),      & !intent(in)
                                rrm_evap(k), rrm_auto(k),                  & !intent(in)
                                rrm_accr(k), Nrm_evap(k),                  & !intent(in)
                                Nrm_auto(k), l_src_adj_enabled,            & !intent(in)
                                l_evap_adj_enabled,                        & !intent(in)
                                rrm_mc(k), Nrm_mc(k),                      & !intent(out)
                                rvm_mc(k), rcm_mc(k), thlm_mc(k),          & !intent(out)
                                adj_terms(k) )                               !intent(out)
    end do ! k = 1, nzt, 1

    ! Clip positive values of Vrr and VNr
    do k = 1, nzt-1, 1

      if ( lh_Vrr(k) > zero ) then
        lh_Vrr(k) = zero
      end if

      if ( lh_VNr(k) > zero ) then
        lh_VNr(k) = zero
      end if

    end do

    cloud_top_level = get_cloud_top_level( nzt, rcm, hydromet, &
                                           hydromet_dim, iiri )

    !!! Mean sedimentation above cloud top should have a value of 0.
    if ( cloud_top_level > 1 ) then
       lh_Vrr(cloud_top_level+1:nzt-1) = zero
       lh_VNr(cloud_top_level+1:nzt-1) = zero
    endif

    ! Set boundary conditions
    rrm_mc(nzt) = zero
    Nrm_mc(nzt) = zero
    rvm_mc(nzt) = zero
    rcm_mc(nzt) = zero
    thlm_mc(nzt) = zero

    ! Statistical sampling
    if ( stats_metadata%l_stats_samp ) then

      if ( stats_metadata%ilh_rrm_src_adj > 0 ) then
        call stat_update_var( stats_metadata%ilh_rrm_src_adj, adj_terms%rrm_src_adj, stats_lh_zt )
      end if

      if ( stats_metadata%ilh_Nrm_src_adj > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nrm_src_adj, adj_terms%Nrm_src_adj, stats_lh_zt )
      end if

      if ( stats_metadata%ilh_rrm_evap_adj > 0 ) then
        call stat_update_var( stats_metadata%ilh_rrm_evap_adj, adj_terms%rrm_evap_adj, stats_lh_zt )
      end if

      if ( stats_metadata%ilh_Nrm_evap_adj > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nrm_evap_adj, adj_terms%Nrm_evap_adj, stats_lh_zt )
      end if

    end if ! stats_metadata%l_stats_samp

  end subroutine adjust_KK_src_means
  !-----------------------------------------------------------------------

end module estimate_scm_microphys_module

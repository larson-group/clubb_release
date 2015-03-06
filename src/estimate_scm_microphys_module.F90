!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module estimate_scm_microphys_module

  implicit none

  public :: est_single_column_tndcy

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine est_single_column_tndcy &
             ( dt, nz, num_samples, d_variables, &
               X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights, &
               pdf_params, hydromet_pdf_params, p_in_Pa, exner, rho, &
               dzq, hydromet, rcm, &
               lh_clipped_vars, &
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
      fstderr, &  ! Constant(s)
      zero, &
      unused_var

    use parameters_model, only: &
      hydromet_dim ! Variable

    use parameters_microphys, only: &
      l_var_covar_src

    use model_flags, only: &
      l_const_Nc_in_cloud

    use corr_varnce_module, only: &
      iiPDF_chi, &
      iiPDF_w

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    use stats_variables, only: &
      stats_zt,           & ! Variable(s)
      stats_sfc,          &
      isilhs_variance_category, &
      irrm_auto,          &
      l_stats_samp

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
      copy_X_nl_into_hydromet_all_pts, &    ! Procedure(s)
      clip_transform_silhs_output,     &
      lh_clipped_variables_type             ! Type

    use parameters_microphys, only: &
      l_silhs_KK_convergence_adj_mean ! Variable(s)

    use array_index, only: &
      iiNrm, & ! Variable(s)
      iirrm

    use lh_microphys_var_covar_module, only: &
      lh_microphys_var_covar_driver   ! Procedure

    use silhs_category_variance_module, only: &
      silhs_category_variance_driver  ! Procedure

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    implicit none

    ! External
#include "microphys_interface.inc"

    intrinsic :: real, all, any

    ! Constant parameters
    logical, parameter :: &
      l_latin_hypercube = .true. ! We are the Latin hypercube!

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      num_samples,   & ! Number of calls to microphysics
      d_variables      ! Number of variates

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs    ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component of each sample

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight for cloud weighted sampling

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! The PDF parameters

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho,        & ! Density on thermo. grid  [kg/m^3]
      dzq,        & ! Difference in height per gridbox   [m]
      rcm           ! Mean liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    type(lh_clipped_variables_type), dimension(nz,num_samples), intent(in) :: &
      lh_clipped_vars   ! Variables from SILHS sample

    ! Output Variables

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_Ncm_mc,     & ! LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
      lh_rcm_mc,     & ! LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
      lh_rvm_mc,     & ! LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
      lh_thlm_mc,    & ! LH estimate of time tndcy. of liquid potential temp.  [K/s]
      lh_rtp2_mc,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]


    ! Local Variables
    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim) :: &
      lh_hydromet_mc_all, & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_all   ! LH est of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      lh_Ncm_mc_all,       & ! LH est of time tendency of cloud droplet concentration  [#/kg/s]
      lh_rcm_mc_all,       & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_all,       & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_all         ! LH est of time tendency of liquid potential temperature     [K/s]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species                    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      Ncn_all_points, &    ! Cloud Nuclei conc. (simplified); Nc=Ncn*H(chi) [#/kg]
      rt_all_points,  &    ! Total water mixing ratio                       [kg/kg]
      thl_all_points, &    ! Liquid potential temperature                   [K]
      chi_all_points, &    ! 's' (Mellor 1977)                              [kg/kg]
      rv_all_points,  &    ! Vapor water                                    [kg/kg]
      rc_all_points,  &    ! Liquid water                                   [kg/kg]
      w_all_points,   &    ! Vertical velocity                              [m/s]
      Nc_all_points        ! Cloud droplet concentration                    [#/kg]

    ! These parameters are not used by the microphysics scheme when SILHS is
    ! turned on.
    real( kind = core_rknd ), dimension(nz) :: &
      cloud_frac_unused, &
      w_std_dev_unused

    type(microphys_stats_vars_type), dimension(num_samples) :: &
      microphys_stats_zt_all, &   ! Statistics variables output from microphysics on zt grid
      microphys_stats_sfc_all     ! Statistics variables on sfc grid

    type(microphys_stats_vars_type) :: &
      microphys_stats_zt_avg, &   ! Average of zt statistics over all sample points
      microphys_stats_sfc_avg     ! Average of sfc statistics over all sample points

    integer :: ivar, sample

    ! ---- Begin Code ----

    w_all_points   = real( X_nl_all_levs(:,:,iiPDF_w), kind=core_rknd )
    chi_all_points = real( X_nl_all_levs(:,:,iiPDF_chi), kind=core_rknd )

    call copy_X_nl_into_hydromet_all_pts &
         ( nz, d_variables, num_samples, & ! Intent(in)
           X_nl_all_levs,                & ! Intent(in)
           hydromet,                     & ! Intent(in)
           hydromet_all_points,          & ! Intent(out)
           Ncn_all_points )                ! Intent(out)

    ! Unpack the lh_clipped_vars structure
    rt_all_points  = lh_clipped_vars%rt
    thl_all_points = lh_clipped_vars%thl
    rc_all_points  = lh_clipped_vars%rc
    rv_all_points  = lh_clipped_vars%rv
    Nc_all_points  = lh_clipped_vars%Nc

    do sample = 1, num_samples

      cloud_frac_unused = unused_var
      w_std_dev_unused  = unused_var
      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub &
           ( dt, nz, & ! In
             l_latin_hypercube, thl_all_points(:,sample), w_all_points(:,sample), p_in_Pa, & ! In
             exner, rho, cloud_frac_unused, w_std_dev_unused, & ! In
             dzq, rc_all_points(:,sample), Nc_all_points(:,sample), & ! In
             chi_all_points(:,sample), rv_all_points(:,sample), & ! In
             hydromet_all_points(:,sample,:), & ! In
             lh_hydromet_mc_all(:,sample,:), lh_hydromet_vel_all(:,sample,:), & ! Out
             lh_Ncm_mc_all(:,sample), & ! Out
             lh_rcm_mc_all(:,sample), lh_rvm_mc_all(:,sample), lh_thlm_mc_all(:,sample), & ! Out
             microphys_stats_zt_all(sample), microphys_stats_sfc_all(sample) ) ! Out

      ! Loop to get new sample
    end do ! sample = 1, num_samples

    if ( l_var_covar_src ) then

      call lh_microphys_var_covar_driver &
           ( nz, num_samples, dt, lh_sample_point_weights, &   ! Intent(in)
             rt_all_points, thl_all_points, w_all_points,  &   ! Intent(in)
             lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, &   ! Intent(in)
             lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc,         &   ! Intent(out)
             lh_wpthlp_mc, lh_rtpthlp_mc )                     ! Intent(out)

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
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               lh_hydromet_vel_all(:,:,ivar) )

      lh_hydromet_mc(:,ivar) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               lh_hydromet_mc_all(:,:,ivar) )

    end forall

    lh_Ncm_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_Ncm_mc_all(:,:) )
    lh_rcm_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_rcm_mc_all(:,:) )
    lh_rvm_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_rvm_mc_all(:,:) )
    lh_thlm_mc= compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_thlm_mc_all(:,:))

    ! Sample variables from microphys_stats_vars objects for statistics
    microphys_stats_zt_avg =  silhs_microphys_stats_avg( num_samples, microphys_stats_zt_all, &
                                                         lh_sample_point_weights )
    microphys_stats_sfc_avg = silhs_microphys_stats_avg( num_samples, microphys_stats_sfc_all, &
                                                         lh_sample_point_weights )

    if ( lh_microphys_type /= lh_microphys_non_interactive ) then
      call microphys_stats_accumulate( microphys_stats_zt_avg, l_stats_samp, stats_zt )
      call microphys_stats_accumulate( microphys_stats_sfc_avg, l_stats_samp, stats_sfc )
    end if

    ! Adjust the mean if l_silhs_KK_convergence_adj_mean is true
    if ( l_silhs_KK_convergence_adj_mean ) then
      call adjust_KK_src_means( dt, nz, exner, rcm, hydromet(:,iirrm),           & ! intent(in)
                                hydromet(:,iiNrm),                               & ! intent(in)
                                microphys_stats_zt_avg, l_stats_samp,            & ! intent(in)
                                lh_hydromet_vel(:,iirrm),                        & ! intent(inout)
                                lh_hydromet_vel(:,iiNrm),                        & ! intent(inout)
                                lh_hydromet_mc(:,iirrm), lh_hydromet_mc(:,iiNrm),& ! intent(out)
                                lh_rvm_mc, lh_rcm_mc, lh_thlm_mc )                 ! intent(out)
    end if

    ! Invoke the SILHS category variance sampler (if desired by user)!!

    if ( l_stats_samp ) then

      if ( allocated( isilhs_variance_category ) ) then

        if( isilhs_variance_category(1) > 0 ) then
          call silhs_category_variance_driver &
               ( nz, num_samples, d_variables, hydromet_dim, X_nl_all_levs, & ! Intent(in)
                 X_mixt_comp_all_levs, microphys_stats_zt_all,              & ! Intent(in)
                 lh_hydromet_mc_all, lh_sample_point_weights, pdf_params,   & ! Intent(in)
                 hydromet_pdf_params )                                        ! Intent(in)
        end if ! isilhs_variance_category(1) > 0

      end if ! allocated( isilhs_variance_category ) 

    end if ! l_stats_samp

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

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight of each SILHS sample point

    ! Output Variables
    type(microphys_stats_vars_type) :: &
      microphys_stats_avg   ! Average of statistical variables

    ! Local Variables
    integer :: ivar, svar, nz, num_vars, stats_index

    real( kind = core_rknd ), dimension(microphys_stats_all(1)%nz,num_samples) :: &
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
        stats_var_all(:,svar) = microphys_stats_all(svar)%output_values(:,ivar)
      end do

      stats_var_avg = compute_sample_mean( nz, num_samples, lh_sample_point_weights, stats_var_all )

      call microphys_put_var( stats_index, stats_var_avg, microphys_stats_avg )

    end do

    return
  end function silhs_microphys_stats_avg
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine adjust_KK_src_means( dt, nz, exner, rcm, rrm, Nrm,         &
                                  microphys_stats_zt, l_stats_samp,     &
                                  lh_Vrr, lh_VNr,                       &
                                  rrm_mc, Nrm_mc,                       &
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

    use KK_utilities, only: &
      get_cloud_top_level         ! Procedure

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: &
      rr_tol, & ! Constant(s)
      Nr_tol, &
      zero

    use stats_type_utilities, only: &
      stat_update_var ! Procedure

    use stats_variables, only: &
      stats_zt,        &
      stats_lh_zt,     &
      irrm_auto, &
      irrm_accr, &
      irrm_cond, &
      iNrm_auto, &
      iNrm_cond, &
      ilh_rrm_auto, &
      ilh_rrm_accr, &
      ilh_rrm_evap, &
      ilh_Nrm_auto, &
      ilh_Nrm_cond, &
      ilh_rrm_src_adj, &
      ilh_Nrm_src_adj, &
      ilh_rrm_cond_adj, &
      ilh_Nrm_cond_adj, &
      im_vol_rad_rain, &
      ilh_m_vol_rad_rain

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, &     ! Type
      microphys_get_var                ! Procedure

    implicit none

    ! Local Constants
    logical, parameter :: &
      ! Whether to adjust rrm_source to not over-deplete cloud water
      l_src_adj_enabled = .true.,  &
      ! Whether to adjust rrm_evap to not over-evaporate rain
      l_evap_adj_enabled = .true.

    ! Input variables
    real( kind = core_rknd ), intent(in) :: &
      dt   ! Model timestep

    integer, intent(in) :: &
      nz   ! Number of vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner,       & ! Exner function                            [-]
      rcm,         & ! Mean liquid water mixing ratio            [kg/kg]
      rrm,      & ! Rain water mixing ration                  [kg/kg]
      Nrm            ! Rain drop concentration                   [num/kg]

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_zt ! Statistics variables                  [units vary]

    logical, intent(in) :: &
      l_stats_samp   ! Whether to sample this timestep

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      lh_Vrr, &         ! Mean sedimentation velocity of < r_r > [m/s]
      lh_VNr            ! Mean sedimentation velocity of < N_r > [m/s]

    ! Output variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rrm_mc, & ! Mean change in rain due to microphysics [(kg/kg)/s] 
      Nrm_mc,    & ! Mean change in Nrm due to microphysics  [(kg/kg)/s]
      rvm_mc,    & ! Time tendency of rvm                    [(kg/kg)/s]
      rcm_mc,    & ! Time tendency of rcm                    [(kg/kg)/s]
      thlm_mc      ! Time tendency of thlm                   [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      rrm_evap, & ! Mean change in rain due to evap           [(kg/kg)/s]
      rrm_auto, & ! Mean change in rain due to autoconversion [(kg/kg)/s]
      rrm_accr, & ! Mean change in rain due to accretion      [(kg/kg)/s]
      Nrm_auto,    & ! Mean change in Nrm due to autoconversion  [(num/kg)/s]
      Nrm_evap       ! Mean change in Nrm due to evaporation     [(num/kg)/s]

    type(KK_microphys_adj_terms_type), dimension(nz) :: &
      adj_terms    ! Adjustment terms returned from the adjustment routine

    integer :: k, cloud_top_level

    !----- Begin code -----

    ! Initialize output
    rrm_mc = zero
    Nrm_mc = zero
    rvm_mc = zero
    rcm_mc = zero
    thlm_mc = zero

    rrm_auto = microphys_get_var( irrm_auto, microphys_stats_zt )
    rrm_accr = microphys_get_var( irrm_accr, microphys_stats_zt )
    rrm_evap = microphys_get_var( irrm_cond, microphys_stats_zt )

    Nrm_auto    = microphys_get_var( iNrm_auto,    microphys_stats_zt )
    Nrm_evap    = microphys_get_var( iNrm_cond,    microphys_stats_zt )

    ! Loop over each vertical level above the lower boundary
    do k = 2, nz, 1

      ! We call KK_microphys_adjust to adjust the means of the mc terms
      call KK_microphys_adjust( dt, exner(k), rcm(k), rrm(k), Nrm(k),      & !intent(in)
                                rrm_evap(k), rrm_auto(k),                  & !intent(in)
                                rrm_accr(k), Nrm_evap(k),                  & !intent(in)
                                Nrm_auto(k), l_src_adj_enabled,            & !intent(in)
                                l_evap_adj_enabled,                        & !intent(in)
                                rrm_mc(k), Nrm_mc(k),                      & !intent(out)
                                rvm_mc(k), rcm_mc(k), thlm_mc(k),          & !intent(out)
                                adj_terms(k) )                               !intent(out)
    end do ! k = 2, nz, 1

    ! Clip positive values of Vrr and VNr
    do k = 1, nz-1, 1

      if ( lh_Vrr(k) > zero ) then
        lh_Vrr(k) = zero
      end if

      if ( lh_VNr(k) > zero ) then
        lh_VNr(k) = zero
      end if

    end do

    cloud_top_level = get_cloud_top_level( nz, rcm )

    !!! Mean sedimentation above cloud top should have a value of 0.
    if ( cloud_top_level > 1 ) then
       lh_Vrr(cloud_top_level+1:nz-1) = zero
       lh_VNr(cloud_top_level+1:nz-1) = zero
    endif

    ! Set boundary conditions
    rrm_mc(1) = zero
    rrm_mc(nz) = zero

    Nrm_mc(1) = zero
    Nrm_mc(nz) = zero

    rvm_mc(1) = zero
    rvm_mc(nz) = zero

    rcm_mc(1) = zero
    rcm_mc(nz) = zero

    thlm_mc(1) = zero
    thlm_mc(nz) = zero

    ! Statistical sampling
    if ( l_stats_samp ) then

      if ( ilh_rrm_src_adj > 0 ) then
        call stat_update_var( ilh_rrm_src_adj, adj_terms%rrm_src_adj, stats_lh_zt )
      end if

      if ( ilh_Nrm_src_adj > 0 ) then
        call stat_update_var( ilh_Nrm_src_adj, adj_terms%Nrm_src_adj, stats_lh_zt )
      end if

      if ( ilh_rrm_cond_adj > 0 ) then
        call stat_update_var( ilh_rrm_cond_adj, adj_terms%rrm_cond_adj, stats_lh_zt )
      end if

      if ( ilh_Nrm_cond_adj > 0 ) then
        call stat_update_var( ilh_Nrm_cond_adj, adj_terms%Nrm_cond_adj, stats_lh_zt )
      end if

      if ( ilh_m_vol_rad_rain > 0 ) then
        call stat_update_var( ilh_m_vol_rad_rain, microphys_get_var( &
             im_vol_rad_rain, microphys_stats_zt ), stats_lh_zt )
      end if

      if ( ilh_rrm_auto > 0 ) then
        call stat_update_var( ilh_rrm_auto, microphys_get_var( &
             irrm_auto, microphys_stats_zt ), stats_lh_zt )
      end if

      if ( ilh_rrm_accr > 0 ) then
        call stat_update_var( ilh_rrm_accr, microphys_get_var( &
             irrm_accr, microphys_stats_zt ), stats_lh_zt )
      end if

      if ( ilh_rrm_evap > 0 ) then
        call stat_update_var( ilh_rrm_evap, microphys_get_var( &
             irrm_cond, microphys_stats_zt ), stats_lh_zt )
      end if

      if ( ilh_Nrm_auto > 0 ) then
        call stat_update_var( ilh_Nrm_auto, microphys_get_var( &
             iNrm_auto, microphys_stats_zt ), stats_lh_zt )
      end if

      if ( ilh_Nrm_cond > 0 ) then
        call stat_update_var( ilh_Nrm_cond, microphys_get_var( &
             iNrm_cond, microphys_stats_zt ), stats_lh_zt )
      end if

    end if ! l_stats_samp

  end subroutine adjust_KK_src_means

  !-----------------------------------------------------------------------
  subroutine lh_moments ( n_samples, lh_weights, nz, & 
                 rt_all_samples, thl_all_samples, w_all_samples, &
                 lh_rtp2, lh_thlp2, &
                 lh_wprtp, lh_wpthlp, &
                 lh_rtpthlp )

  ! Description:
  !   Calculates variances and covariances using LH sample columns
  !
  ! References:
  !   None
  !
  ! TODO: This code assumes nz == gr%nnzp since it references zt2zm;  this is
  ! not necessarily the case.
  !-----------------------------------------------------------------------------

    use grid_class, only: &
      zt2zm    ! Procedures

    use math_utilities, only: &
      compute_sample_mean, & ! functions
      compute_sample_variance, &
      compute_sample_covariance

    use clubb_precision, only: &
      core_rknd ! Constants

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      n_samples        ! Number of sample columns from latin hypercube

    real( kind = core_rknd ), dimension(n_samples), intent(in) :: &
      lh_weights   ! Sample  point weights                   [-]

    real( kind = core_rknd ), dimension(nz,n_samples), intent(in) :: &
      rt_all_samples, &  ! rt columns from latin hypercube   [kg/kg]
      thl_all_samples, & ! thl columns from latin hypercube  [K]
      w_all_samples      ! w columns from latin hypercube    [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rtp2, &    ! Latin hypercube estimate of <rt'^2>     [(kg/kg)^2]
      lh_thlp2, &   ! Latin hypercube estimate of <thl'^2>    [K^2]
      lh_wprtp, &   ! Latin hypercube estimate of <w'rt'>     [m*(kg/kg)/s]
      lh_wpthlp, &  ! Latin hypercube estimate of <w'thl'>    [m*K/s]
      lh_rtpthlp    ! Latin hypercube estimate of <rt'thl'>   [K*(kg/kg)]

    ! Local variables
    real( kind = core_rknd ), dimension(nz) :: &  
      rt_mean, &    ! Latin hypercube estimate of rtm         [kg/kg]
      thl_mean, &   ! Latin hypercube estimate of thlm        [K]
      w_mean, &     ! Latin hypercube estimate of wm          [m/s]
      rtp2_zt, &    ! Estimate of <rt'^2> on the zt grid      [(kg/kg)^2]
      thlp2_zt, &   ! Estimate of <thl'^2> on the zt grid     [K^2]
      wprtp_zt, &   ! Estimate of <w'rt'> on the zt grid      [m*(kg/kg)/s]
      wpthlp_zt, &  ! Estimate of <w'thl'> on the zt grid     [m*K/s]
      rtpthlp_zt    ! Estimate of <rt'thl'> on the zt grid    [K*(kg/kg)]


    ! ---- Begin code ----

    rt_mean = compute_sample_mean( nz, n_samples, lh_weights, rt_all_samples )
    thl_mean = compute_sample_mean( nz, n_samples, lh_weights, thl_all_samples )
    w_mean = compute_sample_mean( nz, n_samples, lh_weights, w_all_samples )

    rtp2_zt = compute_sample_variance( nz, n_samples, rt_all_samples, lh_weights, rt_mean )
    thlp2_zt = compute_sample_variance( nz, n_samples, thl_all_samples, lh_weights, thl_mean )
  
    wprtp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   w_all_samples, w_mean, rt_all_samples, rt_mean ) 
    wpthlp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   w_all_samples, w_mean, thl_all_samples, thl_mean )
    rtpthlp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   rt_all_samples, rt_mean, thl_all_samples, thl_mean ) 

    lh_rtp2 = zt2zm( rtp2_zt )
    lh_thlp2 = zt2zm( thlp2_zt )
    lh_wprtp = zt2zm( wprtp_zt )
    lh_wpthlp = zt2zm( wpthlp_zt )
    lh_rtpthlp = zt2zm( rtpthlp_zt )

    return
  end subroutine lh_moments
  !-----------------------------------------------------------------------------
end module estimate_scm_microphys_module

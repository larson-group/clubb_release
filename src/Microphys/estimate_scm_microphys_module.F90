!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module estimate_scm_microphys_module

  implicit none

  public :: est_silhs_tndcy

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine est_silhs_tndcy( &
               gr, ngrdcol, dt, nzt, nzm, num_samples, &
               pdf_dim, hydromet_dim, hm_metadata, &
               X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights, &
               pdf_params, precip_fracs, p_in_Pa, exner, rho, &
               dzq, hydromet, rcm, &
               lh_rt_clipped, lh_thl_clipped, &
               lh_rc_clipped, lh_rv_clipped, &
               lh_Nc_clipped, &
               l_lh_instant_var_covar_src, &
               saturation_formula, &
               stats,               &
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

    use stats_netcdf, only: &
      stats_type, &
      stats_update, &
      var_on_stats_list

    use math_utilities, only: &
      compute_sample_mean             ! Procedure

    use latin_hypercube_driver_module, only: &
      copy_X_nl_into_hydromet_all_pts    ! Procedure(s)


    use parameters_microphys, only: &
      l_silhs_KK_convergence_adj_mean ! Variable(s)

    use corr_varnce_module, only: &
      hm_metadata_type

    use lh_microphys_var_covar_module, only: &
      lh_microphys_var_covar_driver_api   ! Procedure

    use silhs_category_variance_module, only: &
      silhs_category_variance_driver  ! Procedure

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions

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
      ngrdcol,       & ! Number of model columns
      num_samples,   & ! Number of calls to microphysics
      pdf_dim,       & ! Number of variates
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs    ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component of each sample

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight for cloud weighted sampling

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! The PDF parameters

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho,        & ! Density on thermo. grid  [kg/m^3]
      dzq,        & ! Difference in height per gridbox   [m]
      rcm           ! Mean liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    logical, intent(in) :: &
      l_lh_instant_var_covar_src ! Produce instantaneous var/covar tendencies [-]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type(stats_type), intent(inout) :: &
      stats

    ! Output Variables

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      lh_Ncm_mc,     & ! LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
      lh_rcm_mc,     & ! LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
      lh_rvm_mc,     & ! LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
      lh_thlm_mc       ! LH estimate of time tndcy. of liquid potential temp.  [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      lh_rtp2_mc,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]


    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,hydromet_dim) :: &
      lh_hydromet_mc_all, & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_all   ! LH est of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt) :: &
      lh_Ncm_mc_all,       & ! LH est of time tendency of cloud droplet concentration  [#/kg/s]
      lh_rcm_mc_all,       & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_all,       & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_all         ! LH est of time tendency of liquid potential temperature     [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species                    [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt) :: &
      Ncn_all_points, &    ! Cloud Nuclei conc. (simplified); Nc=Ncn*H(chi) [#/kg]
      chi_all_points, &    ! 's' (Mellor 1977)                              [kg/kg]
      w_all_points         ! Vertical velocity                              [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      lh_rtp2_mc_zt,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc_zt,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc_zt,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc_zt,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc_zt    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]

    ! These parameters are not used by the microphysics scheme when SILHS is
    ! turned on.
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      cloud_frac_unused, &
      w_std_dev_unused

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt) :: &
      rrm_auto_diag_all, &
      rrm_accr_diag_all, &
      rrm_evap_diag_all, &
      Nrm_auto_diag_all, &
      Nrm_evap_diag_all

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rrm_auto_diag_avg, &
      rrm_accr_diag_avg, &
      rrm_evap_diag_avg, &
      Nrm_auto_diag_avg, &
      Nrm_evap_diag_avg

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      lh_rrm_src_adj, lh_Nrm_src_adj, lh_rrm_evap_adj, lh_Nrm_evap_adj

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

    w_all_points   = real( X_nl_all_levs(:,:,:,iiPDF_w), kind=core_rknd )
    chi_all_points = real( X_nl_all_levs(:,:,:,iiPDF_chi), kind=core_rknd )

    call copy_X_nl_into_hydromet_all_pts( &
           nzt, pdf_dim, num_samples, ngrdcol, X_nl_all_levs, & ! In
           hydromet_dim, hm_metadata, hydromet, & ! In
           hydromet_all_points, Ncn_all_points ) ! Out

    cloud_frac_unused = unused_var
    w_std_dev_unused  = unused_var

    do sample = 1, num_samples
      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub( &
             gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, &
             l_latin_hypercube, lh_thl_clipped(:,sample,:), &
             w_all_points(:,sample,:), p_in_Pa, exner, rho, &
             cloud_frac_unused, w_std_dev_unused, dzq, &
             lh_rc_clipped(:,sample,:), lh_Nc_clipped(:,sample,:), &
             chi_all_points(:,sample,:), lh_rv_clipped(:,sample,:), &
             hydromet_all_points(:,sample,:,:), saturation_formula, &
             lh_sample_point_weights(:,sample,:), stats, &
             lh_hydromet_mc_all(:,sample,:,:), &
             lh_hydromet_vel_all(:,sample,:,:), &
             lh_Ncm_mc_all(:,sample,:), lh_rcm_mc_all(:,sample,:), &
             lh_rvm_mc_all(:,sample,:), lh_thlm_mc_all(:,sample,:), &
             rrm_auto_diag_all(:,sample,:), rrm_accr_diag_all(:,sample,:), &
             rrm_evap_diag_all(:,sample,:), Nrm_auto_diag_all(:,sample,:), &
             Nrm_evap_diag_all(:,sample,:) )
    enddo

    if ( l_var_covar_src ) then
      call lh_microphys_var_covar_driver_api( &
             nzt, num_samples, ngrdcol, dt, lh_sample_point_weights, & ! In
             pdf_params, lh_rt_clipped, lh_thl_clipped, w_all_points, & ! In
             lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, & ! In
             l_lh_instant_var_covar_src, & ! In
             lh_rtp2_mc_zt, lh_thlp2_mc_zt, lh_wprtp_mc_zt, & ! Out
             lh_wpthlp_mc_zt, lh_rtpthlp_mc_zt ) ! Out

      ! Convert from the zt grid to the zm grid.
      lh_rtp2_mc = zt2zm_api( nzm, nzt, ngrdcol, gr, lh_rtp2_mc_zt )
      lh_thlp2_mc = zt2zm_api( nzm, nzt, ngrdcol, gr, lh_thlp2_mc_zt )
      lh_wprtp_mc = zt2zm_api( nzm, nzt, ngrdcol, gr, lh_wprtp_mc_zt )
      lh_wpthlp_mc = zt2zm_api( nzm, nzt, ngrdcol, gr, lh_wpthlp_mc_zt )
      lh_rtpthlp_mc = zt2zm_api( nzm, nzt, ngrdcol, gr, lh_rtpthlp_mc_zt )

      ! Stats sampling for LH variance/covariance tendencies.
      if ( stats%l_sample ) then
        call stats_update( "lh_rtp2_mc", lh_rtp2_mc, stats )
        call stats_update( "lh_thlp2_mc", lh_thlp2_mc, stats )
        call stats_update( "lh_wprtp_mc", lh_wprtp_mc, stats )
        call stats_update( "lh_wpthlp_mc", lh_wpthlp_mc, stats )
        call stats_update( "lh_rtpthlp_mc", lh_rtpthlp_mc, stats )
      endif
    else
      lh_rtp2_mc     = zero
      lh_thlp2_mc    = zero
      lh_wprtp_mc    = zero
      lh_wpthlp_mc   = zero
      lh_rtpthlp_mc  = zero
    endif

    ! Grid box average.
    do ivar = 1, hydromet_dim
      lh_hydromet_vel(:,:,ivar) = compute_sample_mean( &
        nzt, num_samples, ngrdcol, lh_sample_point_weights, &
        lh_hydromet_vel_all(:,:,:,ivar) )
      lh_hydromet_mc(:,:,ivar) = compute_sample_mean( &
        nzt, num_samples, ngrdcol, lh_sample_point_weights, &
        lh_hydromet_mc_all(:,:,:,ivar) )
    enddo

    lh_Ncm_mc = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, lh_Ncm_mc_all )
    lh_rcm_mc = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, lh_rcm_mc_all )
    lh_rvm_mc = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, lh_rvm_mc_all )
    lh_thlm_mc = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, lh_thlm_mc_all )
    rrm_auto_diag_avg = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, rrm_auto_diag_all )
    rrm_accr_diag_avg = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, rrm_accr_diag_all )
    rrm_evap_diag_avg = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, rrm_evap_diag_all )
    Nrm_auto_diag_avg = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, Nrm_auto_diag_all )
    Nrm_evap_diag_avg = compute_sample_mean( &
      nzt, num_samples, ngrdcol, lh_sample_point_weights, Nrm_evap_diag_all )

    ! Adjust the mean if l_silhs_KK_convergence_adj_mean is true
    if ( l_silhs_KK_convergence_adj_mean ) then
      call adjust_KK_src_means( &
             dt, nzt, ngrdcol, exner, rcm, hydromet(:,:,iirr), & ! In
             hydromet(:,:,iiNr), hydromet, hydromet_dim, hm_metadata%iiri, & ! In
             rrm_auto_diag_avg, rrm_accr_diag_avg, rrm_evap_diag_avg, & ! In
             Nrm_auto_diag_avg, Nrm_evap_diag_avg, & ! In
             lh_hydromet_vel(:,:,iirr), lh_hydromet_vel(:,:,iiNr), & ! InOut
             lh_hydromet_mc(:,:,iirr), lh_hydromet_mc(:,:,iiNr), & ! Out
             lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, & ! Out
             lh_rrm_src_adj, lh_Nrm_src_adj, & ! Out
             lh_rrm_evap_adj, lh_Nrm_evap_adj ) ! Out

      ! Statistical sampling (moved here with the adjustment diagnostics).
      if ( stats%l_sample ) then
        call stats_update( "lh_rrm_src_adj", lh_rrm_src_adj, stats )
        call stats_update( "lh_Nrm_src_adj", lh_Nrm_src_adj, stats )
        call stats_update( "lh_rrm_evap_adj", lh_rrm_evap_adj, stats )
        call stats_update( "lh_Nrm_evap_adj", lh_Nrm_evap_adj, stats )
      endif
    endif

    if ( stats%l_sample ) then
      ! Invoke the SILHS category variance sampler (if desired by user)!!
      if ( var_on_stats_list( stats, "silhs_var_cat_1" ) ) then
        call silhs_category_variance_driver( &
               ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim, hm_metadata, &
               X_nl_all_levs, X_mixt_comp_all_levs, lh_hydromet_mc_all, &
               lh_sample_point_weights, pdf_params, precip_fracs, stats )
      endif
    endif

  end subroutine est_silhs_tndcy

  !-----------------------------------------------------------------------
  subroutine adjust_KK_src_means( &
               dt, nzt, ngrdcol, exner, rcm, rrm, Nrm, hydromet, & ! In
               hydromet_dim, iiri, rrm_auto, rrm_accr, rrm_evap, & ! In
               Nrm_auto, Nrm_evap, & ! In
               lh_Vrr, lh_VNr, & ! InOut
               rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc, & ! Out
               rrm_src_adj, Nrm_src_adj, rrm_evap_adj, Nrm_evap_adj ) ! Out

  ! Description:
  !   Adjusts the means of microphysics terms for KK microphysics by calling the
  !   KK microphysics adjustment subroutine for every model column.
  !
  ! References:
  !   clubb:ticket:558
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: &
      rc_tol, &
      ri_tol, &
      zero

    use KK_microphys_module, only: &
      KK_microphys_adjust, &
      KK_microphys_adj_terms_type

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      dt                   ! Model timestep [s]

    integer, intent(in) :: &
      nzt, &
      ngrdcol, &
      hydromet_dim, &
      iiri

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      exner, & ! Exner function                            [-]
      rcm,   & ! Mean liquid water mixing ratio            [kg/kg]
      rrm, &               ! Rain water mixing ratio                   [kg/kg]
      Nrm, &               ! Rain drop concentration                   [num/kg]
      rrm_auto, &          ! Mean change in rain due to autoconversion [(kg/kg)/s]
      rrm_accr, &          ! Mean change in rain due to accretion      [(kg/kg)/s]
      rrm_evap, &          ! Mean change in rain due to evap           [(kg/kg)/s]
      Nrm_auto, &          ! Mean change in Nrm due to autoconversion  [(num/kg)/s]
      Nrm_evap             ! Mean change in Nrm due to evaporation     [(num/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      lh_Vrr, &         ! Mean sedimentation velocity of < r_r > [m/s]
      lh_VNr            ! Mean sedimentation velocity of < N_r > [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      rrm_mc, & ! Mean change in rain due to microphysics [(kg/kg)/s]
      Nrm_mc,    & ! Mean change in Nrm due to microphysics  [(kg/kg)/s]
      rvm_mc,    & ! Time tendency of rvm                    [(kg/kg)/s]
      rcm_mc,    & ! Time tendency of rcm                    [(kg/kg)/s]
      thlm_mc, &           ! Time tendency of thlm                   [K/s]
      rrm_src_adj, &
      Nrm_src_adj, &
      rrm_evap_adj, &
      Nrm_evap_adj

    logical, parameter :: &
      l_src_adj_enabled = .true., & ! Whether to adjust rrm_source to not over-deplete cloud water
      l_evap_adj_enabled = .true.   ! Whether to adjust rrm_evap to not over-evaporate rain

    type(KK_microphys_adj_terms_type), dimension(ngrdcol,nzt) :: &
      adj_terms    ! Adjustment terms returned from the adjustment routine

    integer, dimension(ngrdcol) :: &
      cloud_top_level

    integer :: i, k

    !-------------------------- Begin Code --------------------------

    ! Initialize output
    rrm_mc = zero
    Nrm_mc = zero
    rvm_mc = zero
    rcm_mc = zero
    thlm_mc = zero

    do k = 1, nzt
      do i = 1, ngrdcol
        ! We call KK_microphys_adjust to adjust the means of the mc terms
        call KK_microphys_adjust( &
          dt, exner(i,k), rcm(i,k), rrm(i,k), Nrm(i,k), & ! In
          rrm_evap(i,k), rrm_auto(i,k), rrm_accr(i,k), Nrm_evap(i,k), & ! In
          Nrm_auto(i,k), l_src_adj_enabled, l_evap_adj_enabled, & ! In
          rrm_mc(i,k), Nrm_mc(i,k), rvm_mc(i,k), rcm_mc(i,k), thlm_mc(i,k), & ! Out
          adj_terms(i,k) ) ! Out
      enddo
    enddo

    ! Clip positive values of Vrr and VNr
    do k = 1, nzt-1
      do i = 1, ngrdcol
        if ( lh_Vrr(i,k) > zero ) lh_Vrr(i,k) = zero
        if ( lh_VNr(i,k) > zero ) lh_VNr(i,k) = zero
      enddo
    enddo

    cloud_top_level = 1
    if ( iiri > 0 ) then
      do k = nzt, 1, -1
        do i = 1, ngrdcol
          if ( cloud_top_level(i) == 1 ) then
            if ( rcm(i,k) > rc_tol .or. hydromet(i,k,iiri) > ri_tol ) then
              cloud_top_level(i) = k
            endif
          endif
        enddo
      enddo
    else
      do k = nzt, 1, -1
        do i = 1, ngrdcol
          if ( cloud_top_level(i) == 1 .and. rcm(i,k) > rc_tol ) then
            cloud_top_level(i) = k
          endif
        enddo
      enddo
    endif

    !!! Mean sedimentation above cloud top should have a value of 0.
    do k = 2, nzt-1
      do i = 1, ngrdcol
        if ( k > cloud_top_level(i) ) then
          lh_Vrr(i,k) = zero
          lh_VNr(i,k) = zero
        endif
      enddo
    enddo

    ! Set boundary conditions
    rrm_mc(:,nzt) = zero
    Nrm_mc(:,nzt) = zero
    rvm_mc(:,nzt) = zero
    rcm_mc(:,nzt) = zero
    thlm_mc(:,nzt) = zero

    rrm_src_adj = adj_terms%rrm_src_adj
    Nrm_src_adj = adj_terms%Nrm_src_adj
    rrm_evap_adj = adj_terms%rrm_evap_adj
    Nrm_evap_adj = adj_terms%Nrm_evap_adj

  end subroutine adjust_KK_src_means

end module estimate_scm_microphys_module

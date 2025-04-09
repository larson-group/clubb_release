!===============================================================================
module pdf_hydromet_microphys_wrapper

  implicit none

  public :: pdf_hydromet_microphys_prep    ! Procedure(s)

  contains

  !=============================================================================
  subroutine pdf_hydromet_microphys_prep &
             ( gr, ngrdcol, pdf_dim, hydromet_dim,              & ! In
               itime, dt_main, vert_decorr_coef,                & ! In
               Nc_in_cloud, cloud_frac, ice_supersat_frac,      & ! In
               rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp, & ! In
               corr_array_n_cloud, corr_array_n_below,          & ! In
               hm_metadata, pdf_params, clubb_params,           & ! In
               clubb_config_flags, silhs_config_flags,          & ! In
               l_rad_itime, stats_metadata,                     & ! In
               stats_zt, stats_zm, stats_sfc,                   & ! In/Out
               stats_lh_zt, stats_lh_sfc, err_code,             & ! In/Out
               time_clubb_pdf, time_stop, time_start,           & ! In/Out
               hydrometp2,                                      & ! Out
               mu_x_1_n, mu_x_2_n,                              & ! Out
               sigma_x_1_n, sigma_x_2_n,                        & ! Out
               corr_array_1_n, corr_array_2_n,                  & ! Out
               corr_cholesky_mtx_1, corr_cholesky_mtx_2,        & ! Out
               precip_fracs,                                    & ! In/Out
               rtphmp_zt, thlphmp_zt, wp2hmp,                   & ! Out
               X_nl_all_levs, X_mixt_comp_all_levs,             & ! Out
               lh_sample_point_weights,                         & ! Out
               lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,    & ! Out
               lh_rv_clipped, lh_Nc_clipped,                    & ! Out
               hydromet_pdf_params )                              ! Optional(out)

    ! Description:
    ! Prepares the hydrometeor PDF and, when SILHS is in use, sample points,
    ! for the call to microphysics.

    use constants_clubb, only: &
        fstderr    ! Variable(s)

    use grid_class, only: &
        grid    ! Type(s)

    use setup_clubb_pdf_params, only: &
        setup_pdf_parameters    ! Procedure(s)

    use mixed_moment_PDF_integrals, only: &
        hydrometeor_mixed_moments    ! Procedure(s)

#ifdef SILHS
    use silhs_api_module, only: &
        generate_silhs_sample_api, & ! Procedure(s)
        clip_transform_silhs_output_api

    use latin_hypercube_driver_module, only: &
        stats_accumulate_lh
#endif /*SILHS*/

    use pdf_parameter_module, only: &
        pdf_parameter    ! Type(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter,   &  ! Type(s)
        precipitation_fractions

    use corr_varnce_module, only: &
        hm_metadata_type    ! Type(s)

    use model_flags, only: &
        clubb_config_flags_type, & ! Types(s)
        l_silhs_rad                ! Variable(s)

    use parameters_silhs, only: &
        silhs_config_flags_type    ! Types(s)

    use stats_type, only: &
        stats ! Type(s)

    use stats_variables, only: &
        stats_metadata_type    ! Type(s)

    use parameters_microphys, only: &
        microphys_scheme    ! Variable(s)

#ifdef SILHS
    use parameters_microphys, only: &
        lh_microphys_type,     & ! Variable(s)
        lh_microphys_disabled, &
        lh_seed,               &
        lh_num_samples,        &
        lh_sequence_length
#endif /*SILHS*/

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level, & ! Procedure(s)
        clubb_fatal_error             ! Constant

    use clubb_precision, only: &
        core_rknd    ! Variable(s) 

    implicit none

    ! Input Variables
    type (grid), intent(in) :: gr

    integer, intent(in) :: &
      ngrdcol,      & ! Number of grid columns
      pdf_dim,      & ! Number of variables in the correlation array
      hydromet_dim, & ! Number of hydrometeor species
      itime

    real( kind = core_rknd ), intent(in) :: &
      dt_main,          & ! Model timestep                              [s]
      vert_decorr_coef    ! Empirically defined de-correlation constant [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac, & ! Ice supersaturation fraction              [-]
      rho_ds_zt,         & ! Dry, base-state density on thermo. levs.  [kg/m^3]
      Lscale               ! Turbulent Mixing Length                   [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      Kh_zm                ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(in) :: &
      hydromet       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(in) :: &
      corr_array_n_cloud, & ! Prescribed normal space corr. array in cloud  [-]
      corr_array_n_below    ! Prescribed normal space corr. array below cl. [-]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags

    logical, intent(in) :: &
      l_rad_itime

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output Variables
    type (stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc

    integer, intent(inout) :: &
      err_code      ! Error code catching and relaying any errors

    real( kind = core_rknd ), intent(inout) :: &
      time_clubb_pdf, & ! time spent in setup_pdf_parameters and hydrometeor_mixed_moments [s]
      time_start,     & ! help variables to measure the time [s]
      time_stop         ! help variables to measure the time [s]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim), intent(out) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim,pdf_dim), intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim,pdf_dim), intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(out) :: &
      wp2hmp,     & ! Higher-order mixed moment:  < w'^2 hm' > [(m/s)^2<hm un.>]
      rtphmp_zt,  & ! Covariance of rt and hm (on t-levs.)     [(kg/kg)<hm un.>]
      thlphmp_zt    ! Covariance of thl and hm (on t-levs.)    [K<hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt,pdf_dim), intent(out) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(ngrdcol,lh_num_samples,gr%nzt), intent(out) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt), intent(out) :: &
      lh_sample_point_weights

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt), intent(out) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    type(hydromet_pdf_parameter), dimension(ngrdcol,gr%nzt), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables
    logical :: &
      l_calc_weights_all_levs_itime ! determines if vertically correlated sample points are needed

    integer :: i

    logical, parameter :: &
      l_calc_weights_all_levs = .false. ! .false. if all time steps use the same weights at all grid
                                        ! levels
    
    logical, parameter :: &
      l_use_Ncn_to_Nc = .true. ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                               ! Ncn_to_Nc might cause problems with the MG microphysics 
                               ! since the changes made here (Nc-tendency) are not fed into 
                               ! the microphysics



    if ( .not. trim( microphys_scheme ) == "none" ) then

      !$acc update host( cloud_frac, Kh_zm, ice_supersat_frac )

      !$acc if( hydromet_dim > 0 ) update host( wphydrometp )

      !!! Setup the PDF parameters.
      call setup_pdf_parameters( &
              gr, gr%nzm, gr%nzt, ngrdcol, pdf_dim, hydromet_dim, dt_main, & ! In
              Nc_in_cloud, cloud_frac, Kh_zm,                              & ! In
              ice_supersat_frac, hydromet, wphydrometp,                    & ! In
              corr_array_n_cloud, corr_array_n_below,                      & ! In
              hm_metadata,                                                 & ! In
              pdf_params,                                                  & ! In
              clubb_params,                                                & ! In
              clubb_config_flags%iiPDF_type,                               & ! In
              clubb_config_flags%l_use_precip_frac,                        & ! In
              clubb_config_flags%l_predict_upwp_vpwp,                      & ! In
              clubb_config_flags%l_diagnose_correlations,                  & ! In
              clubb_config_flags%l_calc_w_corr,                            & ! In
              clubb_config_flags%l_const_Nc_in_cloud,                      & ! In
              clubb_config_flags%l_fix_w_chi_eta_correlations,             & ! In
              stats_metadata,                                              & ! In
              stats_zt, stats_zm, stats_sfc, err_code,                     & ! In/Out
              hydrometp2,                                                  & ! Out
              mu_x_1_n, mu_x_2_n,                                          & ! Out
              sigma_x_1_n, sigma_x_2_n,                                    & ! Out
              corr_array_1_n, corr_array_2_n,                              & ! Out
              corr_cholesky_mtx_1, corr_cholesky_mtx_2,                    & ! Out
              precip_fracs,                                                & ! In/Out
              hydromet_pdf_params )                                          ! Optional(out)

      ! Error check after setup_pdf_parameters
      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error after setup_pdf_parameters_api"
          return
        end if
      end if

      ! Calculate < rt'hm' >, < thl'hm' >, and < w'^2 hm' >.
      do i = 1, ngrdcol
        call hydrometeor_mixed_moments( gr, gr%nzt, pdf_dim, hydromet_dim,                 & ! In
                                        hydromet(i,:,:), hm_metadata,                      & ! In
                                        mu_x_1_n(i,:,:), mu_x_2_n(i,:,:),                  & ! In
                                        sigma_x_1_n(i,:,:), sigma_x_2_n(i,:,:),            & ! In
                                        corr_array_1_n(i,:,:,:), corr_array_2_n(i,:,:,:),  & ! In
                                        pdf_params, hydromet_pdf_params(i,:),              & ! In
                                        precip_fracs,                                      & ! In
                                        stats_metadata,                                    & ! In
                                        stats_zt(i), stats_zm(i),                          & ! In/Out
                                        rtphmp_zt(i,:,:), thlphmp_zt(i,:,:), wp2hmp(i,:,:) ) ! Out
       end do

      !$acc update device( mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n, corr_array_1_n, corr_array_2_n, &
      !$acc                corr_cholesky_mtx_1, corr_cholesky_mtx_2 )

      !$acc if( hydromet_dim > 0 ) update device( wp2hmp, rtphmp_zt, thlphmp_zt )

    endif ! not microphys_scheme == "none"
      
    ! Measure time in setup_pdf_parameters and hydrometeor_mixed_moments
    call cpu_time(time_stop)
    time_clubb_pdf = time_clubb_pdf + time_stop - time_start
    call cpu_time(time_start) ! initialize timer for SILHS
      
#ifdef SILHS
    !----------------------------------------------------------------
    ! Compute subcolumns if enabled
    !----------------------------------------------------------------

    if ( lh_microphys_type /= lh_microphys_disabled .or. l_silhs_rad ) then
      
      ! Calculate sample weights separately at all grid levels when radiation is not called
      l_calc_weights_all_levs_itime = l_calc_weights_all_levs .and. .not. l_rad_itime

      ! The profile of CLUBB's mixing length, Lscale, is passed into
      ! subroutine generate_silhs_sample for use in calculating the vertical
      ! correlation coefficient.  Rather than output Lscale directly, its
      ! value can be calculated from other fields that are already output from
      ! subroutine advance_clubb_core.  The equation relating Lscale to eddy
      ! diffusivity is:
      !
      ! Kh = c_K * Lscale * sqrt( TKE ).
      !
      ! Kh is available, TKE can be calculated from wp2, up2, and vp2, all of
      ! which are available, and c_K is easily extracted from CLUBB's tunable
      ! parameters.  The equation for Lscale is:
      !
      ! Lscale = Kh / ( c_K * sqrt( TKE ) ).
      !
      ! Since Kh and TKE are output on momentum grid levels, the resulting
      ! calculation of Lscale is also found on momentum levels.  It needs to
      ! be interpolated back to thermodynamic (midpoint) grid levels for
      ! further use.

      ! Calculate TKE
      ! em is calculated but never used??
      ! if ( .not. clubb_config_flags%l_tke_aniso ) then
      !   ! tke is assumed to be 3/2 of wp2
      !   do k = 1, gr%nzm
      !     do i = 1, ngrdcol
      !       em(i,k) = three_halves * wp2(i,k)
      !     end do
      !   end do
      ! else
      !   do k = 1, gr%nzm
      !     do i = 1, ngrdcol
      !       em(i,k) = one_half * ( wp2(i,k) + vp2(i,k) + up2(i,k) )
      !     end do
      !   end do
      ! endif

      call generate_silhs_sample_api( &
             itime, pdf_dim, lh_num_samples, lh_sequence_length, gr%nzt, ngrdcol, & ! In
             l_calc_weights_all_levs_itime,                                & ! In
             gr, pdf_params, gr%dzt, Lscale,                               & ! In
             lh_seed, hm_metadata,                                         & ! In
             !rho_ds_zt,                                                    & ! In
             mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,                 & ! In
             corr_cholesky_mtx_1, corr_cholesky_mtx_2,                     & ! In
             precip_fracs, silhs_config_flags,                             & ! In
             vert_decorr_coef,                                             & ! In
             stats_metadata,                                               & ! In
             stats_lh_zt, stats_lh_sfc,                                    & ! InOut
             X_nl_all_levs, X_mixt_comp_all_levs,                          & ! Out
             lh_sample_point_weights ) ! Out


      call clip_transform_silhs_output_api( gr%nzt, ngrdcol, lh_num_samples,        & ! In
                                            pdf_dim, hydromet_dim, hm_metadata,     & ! In
                                            X_mixt_comp_all_levs,                   & ! In
                                            X_nl_all_levs,                          & ! Inout
                                            pdf_params, l_use_Ncn_to_Nc,            & ! In
                                            lh_rt_clipped, lh_thl_clipped,          & ! Out
                                            lh_rc_clipped, lh_rv_clipped,           & ! Out
                                            lh_Nc_clipped                           ) ! Out

      if ( stats_metadata%l_stats_samp ) then     

        !$acc update host( rho_ds_zt, lh_sample_point_weights, X_nl_all_levs, &
        !$acc              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped )

        do i = 1, ngrdcol
          call stats_accumulate_lh( &
                gr, gr%nzt, lh_num_samples, pdf_dim, rho_ds_zt(i,:),     & ! In
                hydromet_dim, hm_metadata,                               & ! In
                lh_sample_point_weights(i,:,:),  X_nl_all_levs(i,:,:,:), & ! In
                lh_rt_clipped(i,:,:), lh_thl_clipped(i,:,:),             & ! In
                lh_rc_clipped(i,:,:), lh_rv_clipped(i,:,:),              & ! In
                lh_Nc_clipped(i,:,:),                                    & ! In
                stats_metadata,                                          & ! In
                stats_lh_zt(i), stats_lh_sfc(i) )                          ! intent(inout)
        end do
      end if

    end if ! lh_microphys_enabled

#endif /* SILHS */

    return

  end subroutine pdf_hydromet_microphys_prep

!===============================================================================

end module pdf_hydromet_microphys_wrapper

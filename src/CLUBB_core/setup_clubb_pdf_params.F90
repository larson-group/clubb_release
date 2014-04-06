! $Id$
!===============================================================================
module setup_clubb_pdf_params

  implicit none

  private

  public :: setup_pdf_parameters, &
            unpack_pdf_params,    &
            compute_mean_stdev,   &
            normalize_mean_stdev, &
            compute_corr,         &
            normalize_corr

  private :: component_means_hydromet, &
             precip_fraction,          &
             component_mean_hm_ip,     &
             component_stdev_hm_ip,    &
             component_corr_wx,        &
             component_corr_st,        &
             component_corr_whm_ip,    &
             component_corr_xhm_ip,    &
             component_corr_hmxhmy_ip, &
             calc_corr_whm,            &
             pdf_param_hm_stats,       &
             pdf_param_log_hm_stats,   &
             pack_pdf_params

  ! Prescribed parameters are set to in-cloud or outside-cloud (below-cloud)
  ! values based on whether or not cloud water mixing ratio has a value of at
  ! least rc_tol.  However, this does not take into account the amount of
  ! cloudiness in a component, just whether or not there is any cloud in the
  ! component.  The option l_interp_prescribed_params allows for an interpolated
  ! value between the in-cloud and below-cloud parameter value based on the
  ! component cloud fraction.
  logical, parameter, private :: &
    l_interp_prescribed_params = .false.

  contains

  !=============================================================================
  subroutine setup_pdf_parameters( nz, d_variables, dt, wm_zt, rho, &           ! Intent(in)
                                   wp2_zt, Ncm, Nc_in_cloud, rcm, cloud_frac, & ! Intent(in)
                                   ice_supersat_frac, hydromet, wphydrometp, &  ! Intent(in)
                                   corr_array_cloud, corr_array_below, &        ! Intent(in)
                                   pdf_params, l_stats_samp, &                  ! Intent(in)
                                   mu_x_1_n, mu_x_2_n, &                        ! Intent(out)
                                   sigma_x_1_n, sigma_x_2_n, &                  ! Intent(out)
                                   corr_array_1_n, corr_array_2_n, &            ! Intent(out)
                                   corr_cholesky_mtx_1, corr_cholesky_mtx_2, &  ! Intent(out)
                                   hydromet_pdf_params, hydrometp2 )            ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr,    & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        rc_tol,         &
        Ncn_tol,        &
        cloud_frac_min, &
        fstderr,        &
        zero_threshold

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use model_flags, only: &
        l_use_precip_frac,   & ! Flag(s)
        l_calc_w_corr,       &
        l_use_modified_corr

    use parameters_microphys, only: &
        l_const_Nc_in_cloud, & ! Flag(s)
        l_predictnc,         &
        hydromet_list,       & ! Variable(s)
        hydromet_tol

    use advance_windm_edsclrm_module, only: &
        xpwp_fnc

    use variables_diagnostic_module, only: &
        Kh_zm

    use parameters_tunable, only: &
        c_K_hm

    use PDF_utilities, only: &
        calc_xp2  ! Procedure(s)

    use clip_explicit, only: &
        clip_covar_level, & ! Procedure(s)
        clip_wphydrometp    ! Variables(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        time_precision, &
        dp

    use matrix_operations, only: &
        Cholesky_factor ! Procedure(s)

    use stats_type, only: &
        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: &
        irr1,           & ! Variable(s)
        irr2,           &
        iNr1,           &
        iNr2,           &
        iprecip_frac,   &
        iprecip_frac_1, &
        iprecip_frac_2, &
        iNcnm,          &
        irrp2_zt,       &
        iNrp2_zt,       &
        zt

    use model_flags, only: &
        l_diagnose_correlations ! Variable(s)

    use diagnose_correlations_module, only: &
        diagnose_correlations, & ! Procedure(s)
        calc_cholesky_corr_mtx_approx

    use corr_matrix_module, only: &
        sigma2_on_mu2_ip_array_cloud, & ! Variable(s)
        sigma2_on_mu2_ip_array_below, &
        iiPDF_rrain, &
        iiPDF_Nr, &
        iiPDF_Ncn, &
        iiPDF_s_mellor, &
        iiPDF_t_mellor, &
        iiPDF_w

    use index_mapping, only: &
        hydromet2pdf_idx    ! Procedure(s)

    use array_index, only: &
        iirrainm, & ! Variable(s)
        iiNrm

    use error_code, only : &
        clubb_at_least_debug_level   ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = time_precision ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wm_zt,       & ! Mean vertical velocity, <w>, on thermo. levels  [m/s]
      rho,         & ! Density                                         [kg/m^3]
      wp2_zt,      & ! Variance of w, <w'^2> (interp. to t-levs.)      [m^2/s^2]
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.) [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration     [num/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_cloud, & ! Prescribed correlation array in cloud      [-]
      corr_array_below    ! Prescribed correlation array below cloud   [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2_n    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Local Variables
    real( kind = dp ), dimension(d_variables,d_variables,nz) :: &
      corr_cholesky_mtx_1_dp, & ! Used for call to Cholesky_factor, requires dp
      corr_cholesky_mtx_2_dp

    real( kind = core_rknd ), dimension(d_variables,d_variables) :: &
      corr_mtx_approx_1,   & ! Approximated corr. matrix (C = LL'), 1st comp. [-]
      corr_mtx_approx_2      ! Approximated corr. matrix (C = LL'), 2nd comp. [-]

    real( kind = core_rknd ), dimension(nz) :: &
      rc1,         & ! Mean of r_c (1st PDF component)              [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)              [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)           [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)           [-]
      mixt_frac      ! Mixture fraction                             [-]

    real( kind = core_rknd ), dimension(nz) :: &
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >        [num/kg]

    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zm,     & ! Covariance of s and w (momentum levels)   [(m/s)(kg/kg)]
      wpNcnp_zm,   & ! Covariance of N_cn and w (momentum levs.) [(m/s)(num/kg)]
      wpsp_zt,     & ! Covariance of s and w on t-levs           [(m/s)(kg/kg)]
      wpNcnp_zt      ! Covariance of N_cn and w on t-levs        [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)    [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)    [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hydrometp2_zt,  & ! Variance of a hydrometeor (overall); t-lev   [units^2]
      wphydrometp_zt    ! Covariance of w and hm interp. to t-levs. [(m/s)units]

    real( kind = core_rknd ), dimension(nz) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    real( kind = core_rknd ), dimension(d_variables) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = dp ), dimension(d_variables) :: &
      corr_array_scaling

    real( kind = core_rknd ), dimension(d_variables) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      wphydrometp_chnge    ! Change in wphydrometp_zt: covar. clip. [(m/s)units]

    logical :: l_corr_array_scaling

    ! Flags used for covariance clipping of <w'hm'>.
    logical, parameter :: &
      l_first_clip_ts = .true., & ! First instance of clipping in a timestep.
      l_last_clip_ts  = .true.    ! Last instance of clipping in a timestep.

    character(len=10) :: &
      hydromet_name    ! Name of a hydrometeor

    integer :: pdf_idx  ! Index of precipitating hydrometeor in PDF array.

    integer :: k, i  ! Loop indices

    ! ---- Begin Code ----


    ! Assertion check
    ! Check that all hydrometeors are positive otherwise exit the program
    if ( clubb_at_least_debug_level( 2 ) ) then
       do i = 1, hydromet_dim
          if ( any( hydromet(:,i) < zero_threshold ) ) then
             hydromet_name = hydromet_list(i)
             do k = 1, nz
                if ( hydromet(k,i) < zero_threshold ) then

                   ! Write error message
                   write(fstderr,*) trim( hydromet_name )// " = ", hydromet(k,i)," < ", &
                                   zero_threshold, &
                                   " at beginning of setup_pdf_parameters at k= ", k

                   ! Exit program
                   stop "Exiting..."

                endif ! hydromet(k,i) < 0
             enddo ! k = 1, nz
          endif ! hydromet(:,i) < 0
       enddo ! i = 1, hydromet_dim

    endif !clubb_at_least_debug_level( 2 )

    ! Interpolate the covariances of w and precipitating hydrometeors to
    ! thermodynamic grid levels.
    do i = 1, hydromet_dim, 1

       wphydrometp_zt(:,i) = zm2zt( wphydrometp(:,i) )

       ! When the mean value of a precipitating hydrometeor is below tolerance
       ! value, it is considered to have a value of 0, and the precipitating
       ! hydrometeor does not vary over the grid level.  Any covariance
       ! involving that precipitating hydrometeor also has a value of 0 at that
       ! grid level. 
       do k = 1, nz, 1
          if ( hydromet(k,i) <= hydromet_tol(i) ) then
             wphydrometp_zt(k,i) = zero
          endif
       enddo ! k = 1, nz, 1

    enddo ! i = 1, hydromet_dim, 1

    ! Setup some of the PDF parameters
    rc1         = pdf_params%rc1
    rc2         = pdf_params%rc2
    cloud_frac1 = pdf_params%cloud_frac1
    cloud_frac2 = pdf_params%cloud_frac2
    mixt_frac   = pdf_params%mixt_frac

    ! Component mean values for r_r and N_r, and precipitation fraction.
    if ( l_use_precip_frac ) then

       call component_means_hydromet( nz, hydromet, rho, rc1, rc2, &
                                      mixt_frac, l_stats_samp, &
                                      hm1, hm2 )

       call precip_fraction( nz, hydromet, hm1, hm2, &
                             cloud_frac, cloud_frac1, mixt_frac, &
                             ice_supersat_frac, &
                             precip_frac, precip_frac_1, precip_frac_2 )

    else

       hm1 = hydromet
       hm2 = hydromet

       precip_frac   = one
       precip_frac_1 = one
       precip_frac_2 = one

    endif

    ! Calculate <N_cn> from <N_c>, whether <N_c> is predicted or based on a
    ! prescribed value, and whether the in-cloud value is constant or varying.
    if ( l_predictnc .and. .not. l_const_Nc_in_cloud ) then
       ! The mean of N_c is predicted and the value of N_c in-cloud may vary.
       ! I will call the full equation here later.  Brian; 1/28/2014.
       do k = 1, nz
          if ( cloud_frac(k) > cloud_frac_min ) then
             Ncnm(k) = Ncm(k) / cloud_frac(k)
          else
             ! The model is producing a positive mean cloud droplet
             ! concentration, but is not producing any cloud.  Set Ncnm to the
             ! constant, prescribed value of Nc (which is still set to a default
             ! value even when Ncm is predicted).  This will avoid a huge value
             ! of Ncnm.
             Ncnm(k) = Nc_in_cloud(k)
          endif
       enddo ! k = 1, nz
    elseif ( l_predictnc .and. l_const_Nc_in_cloud ) then
       ! The mean of N_c is predicted and the value of N_c in-cloud is constant.
       do k = 1, nz
          if ( cloud_frac(k) > cloud_frac_min ) then
             Ncnm(k) = Ncm(k) / cloud_frac(k)
          else
             ! The model is producing a positive mean cloud droplet
             ! concentration, but is not producing any cloud.  Set Ncnm to the
             ! constant, prescribed value of Nc (which is still set to a default
             ! value even when Ncm is predicted).  This will avoid a huge value
             ! of Ncnm.
             Ncnm(k) = Nc_in_cloud(k)
          endif
       enddo ! k = 1, nz
    elseif ( .not. l_const_Nc_in_cloud .and. .not. l_predictNc ) then
       ! The value of N_c in-cloud is based on a prescribed parameter, which is
       ! used as the in-cloud mean of N_c.  However, N_c in-cloud is allowed to
       ! vary around this prescribed mean value.  The value of N_cn also varies.
       ! Find the mean of N_cn, <N_cn>.
       ! I will call the full equation here later.  Brian; 1/28/2014.
       Ncnm = Nc_in_cloud
    elseif ( l_const_Nc_in_cloud .and. .not. l_predictNc ) then
       ! The value of N_c in-cloud is constant, prescribed parameter.  The value
       ! of N_cn is also constant at any grid level.
       Ncnm = Nc_in_cloud
    endif

    ! Calculate correlations involving w by first calculating total covariances
    ! involving w (<w'r_r'>, etc.) using the down-gradient approximation.
    if ( l_calc_w_corr ) then

       ! Calculate the covariances of w with the hydrometeors
       do k = 1, nz
          wpsp_zm(k) = pdf_params(k)%mixt_frac &
                       * ( one - pdf_params(k)%mixt_frac ) &
                       * ( pdf_params(k)%s1 - pdf_params(k)%s2 ) &
                       * ( pdf_params(k)%w1 - pdf_params(k)%w2 )
       enddo

       wpNcnp_zm(1:nz-1) = xpwp_fnc( -c_K_hm * Kh_zm(1:nz-1), Ncnm(1:nz-1), &
                                     Ncnm(2:nz), gr%invrs_dzm(1:nz-1) )

       ! Boundary conditions; We are assuming zero flux at the top.
       wpNcnp_zm(nz) = zero

       ! Interpolate the covariances to thermodynamic grid levels.
       wpsp_zt = zm2zt( wpsp_zm )
       wpNcnp_zt = zm2zt( wpNcnp_zm )

       ! When the mean value of Ncn is below tolerance value, it is considered
       ! to have a value of 0, and Ncn does not vary over the grid level.  Any
       ! covariance involving Ncn also has a value of 0 at that grid level.
       do k = 1, nz, 1
          if ( Ncnm(k) <= Ncn_tol ) then
             wpNcnp_zt(k) = zero
          endif
       enddo ! k = 1, nz, 1

    endif ! l_calc_w_corr

    ! Initialize the correlation Cholesky matrices
    corr_cholesky_mtx_1 = zero
    corr_cholesky_mtx_2 = zero


    ! Statistics
    if ( l_stats_samp ) then

       if ( iirrainm > 0 ) then

          if ( irr1 > 0 ) then
             ! Mean rain water mixing ratio in PDF component 1.
             call stat_update_var( irr1, hm1(:,iirrainm), zt )
          endif

          if ( irr2 > 0 ) then
             ! Mean rain water mixing ratio in PDF component 2.
             call stat_update_var( irr2, hm2(:,iirrainm), zt )
          endif

       endif ! iirrainm > 0

       if ( iiNrm > 0 ) then

          if ( iNr1 > 0 ) then
             ! Mean rain drop concentration in PDF component 1.
             call stat_update_var( iNr1, hm1(:,iiNrm), zt )
          endif

          if ( iNr2 > 0 ) then
             ! Mean rain drop concentration in PDF component 2.
             call stat_update_var( iNr2, hm2(:,iiNrm), zt )
          endif

       endif ! iiNrm > 0

       if ( iprecip_frac > 0 ) then
          ! Overall precipitation fraction.
          call stat_update_var( iprecip_frac, precip_frac, zt )
       endif

       if ( iprecip_frac_1 > 0 ) then
          ! Precipitation fraction in PDF component 1.
          call stat_update_var( iprecip_frac_1, precip_frac_1, zt )
       endif

       if ( iprecip_frac_2 > 0 ) then
          ! Precipitation fraction in PDF component 2.
          call stat_update_var( iprecip_frac_2, precip_frac_2, zt )
       endif

       if ( iNcnm > 0 ) then
          ! Mean simplified cloud nuclei concentration (overall).
          call stat_update_var( iNcnm, Ncnm, zt )
       endif

    endif


    !!! Setup PDF parameters loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    ! Now also including "model lower boundary" -- Eric Raut Aug 2013
    do k = 1, nz, 1

       if ( rc1(k) > rc_tol ) then
          sigma2_on_mu2_ip_1 = sigma2_on_mu2_ip_array_cloud
       else
          sigma2_on_mu2_ip_1 = sigma2_on_mu2_ip_array_below
       endif

       if ( rc2(k) > rc_tol ) then
          sigma2_on_mu2_ip_2 = sigma2_on_mu2_ip_array_cloud
       else
          sigma2_on_mu2_ip_2 = sigma2_on_mu2_ip_array_below
       endif

       !!! Calculate the means and standard deviations involving PDF variables
       !!! -- w, s, t, N_cn, and any precipitating hydrometeors (hm in-precip)
       !!! -- for each PDF component.
       call compute_mean_stdev( Ncnm(k), rc1(k), rc2(k), &            ! Intent(in)
                                cloud_frac1(k), cloud_frac2(k), &     ! Intent(in)
                                hm1(k,:), hm2(k,:), &                 ! Intent(in)
                                precip_frac_1(k), precip_frac_2(k), & ! Intent(in)
                                sigma2_on_mu2_ip_array_cloud, &       ! Intent(in)
                                sigma2_on_mu2_ip_array_below, &       ! Intent(in)
                                pdf_params(k), d_variables, &         ! Intent(in)
                                mu_x_1, mu_x_2, &                     ! Intent(out)
                                sigma_x_1, sigma_x_2 )                ! Intent(out)


       !!! Calculate the normalized means and normalized standard deviations
       !!! involving precipitating hydrometeors (hm in-precip) and N_cn --
       !!! ln hm and ln N_cn -- for each PDF component.
       call normalize_mean_stdev( hm1(k,:), hm2(k,:), Ncnm(k), d_variables, &
                                  mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                                  sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                                  mu_x_1_n(:,k), mu_x_2_n(:,k), &
                                  sigma_x_1_n(:,k), sigma_x_2_n(:,k) )

       ! Calculate the overall variance of a precipitating hydrometeor (hm),
       ! <hm'^2>.
       do i = 1, hydromet_dim, 1

          if ( hydromet(k,i) > hydromet_tol(i) ) then

             ! There is some of the hydrometeor species found at level k.
             ! Calculate the variance (overall) of the hydrometeor.

             pdf_idx = hydromet2pdf_idx(i)

             hydrometp2_zt(k,i) &
             = calc_xp2( mu_x_1(pdf_idx), mu_x_2(pdf_idx), &
                         mu_x_1_n(pdf_idx,k), mu_x_2_n(pdf_idx,k), &
                         sigma_x_1(pdf_idx), sigma_x_2(pdf_idx), &
                         sigma_x_1_n(pdf_idx,k), sigma_x_2_n(pdf_idx,k), &
                         mixt_frac(k), precip_frac_1(k), precip_frac_2(k), &
                         hydromet(k,i) )

          else ! hydromet(k,i) = 0.

             hydrometp2_zt(k,i) = zero

          endif

          ! Statistics
          if ( l_stats_samp ) then

             if ( i == iirrainm ) then

                ! Variance of rain water mixing ratio, <r_r'^2>.
                call stat_update_var_pt( irrp2_zt, k, hydrometp2_zt(k,i), zt )

             elseif ( i == iiNrm ) then

                ! Variance of rain drop concentration, <N_r'^2>.
                call stat_update_var_pt( iNrp2_zt, k, hydrometp2_zt(k,i), zt )

             endif

          endif ! l_stats_samp

          ! Clip the value of covariance <w'hm'> on thermodynamic levels.
          call clip_covar_level( clip_wphydrometp, k, l_first_clip_ts, &
                                 l_last_clip_ts, dt, wp2_zt(k), &
                                 hydrometp2_zt(k,i), &
                                 wphydrometp_zt(k,i), wphydrometp_chnge(k,i) )

       enddo ! i = 1, hydromet_dim, 1

       if ( l_diagnose_correlations ) then

          if ( rcm(k) > rc_tol ) then

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_array_1 )                   ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_array_2 )                   ! Intent(out)

          else

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_array_1 )                   ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_array_2 )                   ! Intent(out)

          endif

       else ! if .not. l_diagnose_correlations

          call compute_corr( wm_zt(k), rc1(k), rc2(k), cloud_frac1(k), &
                             cloud_frac2(k), wpsp_zt(k), wpNcnp_zt(k), &
                             sqrt(wp2_zt(k)), mixt_frac(k), precip_frac_1(k), &
                             precip_frac_2(k), wphydrometp_zt(k,:), &
                             mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                             corr_array_cloud, corr_array_below, &
                             pdf_params(k), d_variables, &
                             corr_array_1, corr_array_2 )

       endif ! l_diagnose_correlations

       !!! Statistics
       !!! We should generalize the statistics output to write all the
       !!! hydrometeor species to disk.
       call pdf_param_hm_stats(  mu_x_1(iiPDF_rrain), mu_x_2(iiPDF_rrain), &
                                 mu_x_1(iiPDF_Nr), mu_x_2(iiPDF_Nr), &
                                 mu_x_1(iiPDF_Ncn), mu_x_2(iiPDF_Ncn), &
                                 sigma_x_1(iiPDF_rrain), sigma_x_2(iiPDF_rrain), &
                                 sigma_x_1(iiPDF_Nr), sigma_x_2(iiPDF_Nr), &
                                 sigma_x_1(iiPDF_Ncn), sigma_x_2(iiPDF_Ncn), &
                                 corr_array_1(iiPDF_rrain, iiPDF_w), &
                                 corr_array_2(iiPDF_rrain, iiPDF_w), &
                                 corr_array_1(iiPDF_Nr, iiPDF_w), &
                                 corr_array_2(iiPDF_Nr, iiPDF_w), &
                                 corr_array_1(iiPDF_Ncn, iiPDF_w), &
                                 corr_array_2(iiPDF_Ncn, iiPDF_w), &
                                 corr_array_1(iiPDF_rrain, iiPDF_s_mellor), &
                                 corr_array_2(iiPDF_rrain, iiPDF_s_mellor), &
                                 corr_array_1(iiPDF_Nr, iiPDF_s_mellor), &
                                 corr_array_2(iiPDF_Nr, iiPDF_s_mellor), &
                                 corr_array_1(iiPDF_Ncn, iiPDF_s_mellor), &
                                 corr_array_2(iiPDF_Ncn, iiPDF_s_mellor), &
                                 corr_array_1(iiPDF_rrain, iiPDF_t_mellor), &
                                 corr_array_2(iiPDF_rrain, iiPDF_t_mellor), &
                                 corr_array_1(iiPDF_Nr, iiPDF_t_mellor), &
                                 corr_array_2(iiPDF_Nr, iiPDF_t_mellor), &
                                 corr_array_1(iiPDF_Ncn, iiPDF_t_mellor), &
                                 corr_array_2(iiPDF_Ncn, iiPDF_t_mellor), &
                                 corr_array_1(iiPDF_Nr, iiPDF_rrain), &
                                 corr_array_2(iiPDF_Nr, iiPDF_rrain), &
                                 k, l_stats_samp )

       !!! Calculate the correlations involving the natural logarithm of
       !!! precipitating hydrometeors, ln hm (for example, ln r_r and ln N_r),
       !!! and ln N_cn for each PDF component.
       call normalize_corr( d_variables, sigma_x_1_n(:,k), sigma_x_2_n(:,k), &
                            sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                            corr_array_1, corr_array_2, &
                            corr_array_1_n(:,:,k), corr_array_2_n(:,:,k) )


       !!! Statistics
       !!! We should generalize the statistics output to write all the
       !!! hydrometeor species to disk.
       call pdf_param_log_hm_stats(  mu_x_1_n(iiPDF_rrain, k), mu_x_2_n(iiPDF_rrain, k), &
                                     mu_x_1_n(iiPDF_Nr, k), mu_x_2_n(iiPDF_Nr, k), &
                                     mu_x_1_n(iiPDF_Ncn, k), mu_x_2_n(iiPDF_Ncn, k), &
                                     sigma_x_1_n(iiPDF_rrain, k), sigma_x_2_n(iiPDF_rrain, k), &
                                     sigma_x_1_n(iiPDF_Nr, k), sigma_x_2_n(iiPDF_Nr, k), &
                                     sigma_x_1_n(iiPDF_Ncn, k), sigma_x_2_n(iiPDF_Ncn, k), &
                                     corr_array_1_n(iiPDF_rrain, iiPDF_w, k), &
                                     corr_array_2_n(iiPDF_rrain, iiPDF_w, k), &
                                     corr_array_1_n(iiPDF_Nr, iiPDF_w, k), &
                                     corr_array_2_n(iiPDF_Nr, iiPDF_w, k), &
                                     corr_array_1_n(iiPDF_Ncn, iiPDF_w, k), &
                                     corr_array_2_n(iiPDF_Ncn, iiPDF_w, k), &
                                     corr_array_1_n(iiPDF_rrain, iiPDF_s_mellor, k), &
                                     corr_array_2_n(iiPDF_rrain, iiPDF_s_mellor, k), &
                                     corr_array_1_n(iiPDF_Nr, iiPDF_s_mellor, k), &
                                     corr_array_2_n(iiPDF_Nr, iiPDF_s_mellor, k), &
                                     corr_array_1_n(iiPDF_Ncn, iiPDF_s_mellor, k), &
                                     corr_array_2_n(iiPDF_Ncn, iiPDF_s_mellor, k), &
                                     corr_array_1_n(iiPDF_rrain, iiPDF_t_mellor, k), &
                                     corr_array_2_n(iiPDF_rrain, iiPDF_t_mellor, k), &
                                     corr_array_1_n(iiPDF_Nr, iiPDF_t_mellor, k), &
                                     corr_array_2_n(iiPDF_Nr, iiPDF_t_mellor, k), &
                                     corr_array_1_n(iiPDF_Ncn, iiPDF_t_mellor, k), &
                                     corr_array_2_n(iiPDF_Ncn, iiPDF_t_mellor, k), &
                                     corr_array_1_n(iiPDF_Nr, iiPDF_rrain, k), &
                                     corr_array_2_n(iiPDF_Nr, iiPDF_rrain, k), &
                                     k, l_stats_samp )

       !!! Has to be generalized
       call pack_pdf_params( hm1(k,:), hm2(k,:), mu_x_1, mu_x_2, &  ! Intent(in)
                             sigma_x_1, sigma_x_2, d_variables, &   ! Intent(in)
                             precip_frac(k), precip_frac_1(k), precip_frac_2(k), & ! Intent(in)
                             hydromet_pdf_params(k) )               ! Intent(out)

       if ( l_use_modified_corr ) then

          if ( l_diagnose_correlations ) then

             call calc_cholesky_corr_mtx_approx &
                            ( d_variables, corr_array_1_n(:,:,k), &           ! intent(in)
                              corr_cholesky_mtx_1(:,:,k), corr_mtx_approx_1 ) ! intent(out)

             call calc_cholesky_corr_mtx_approx &
                            ( d_variables, corr_array_2_n(:,:,k), &           ! intent(in)
                              corr_cholesky_mtx_2(:,:,k), corr_mtx_approx_2 ) ! intent(out)

             corr_array_1_n(:,:,k) = corr_mtx_approx_1
             corr_array_2_n(:,:,k) = corr_mtx_approx_2

          else

             ! Compute choleksy factorization for the correlation matrix (out of
             ! cloud)
             call Cholesky_factor( d_variables, real(corr_array_1_n(:,:,k), kind = dp), & ! In
                                   corr_array_scaling, corr_cholesky_mtx_1_dp(:,:,k), &  ! Out
                                   l_corr_array_scaling ) ! Out

             call Cholesky_factor( d_variables, real(corr_array_2_n(:,:,k), kind = dp), & ! In
                                   corr_array_scaling, corr_cholesky_mtx_2_dp(:,:,k), &  ! Out
                                   l_corr_array_scaling ) ! Out
             corr_cholesky_mtx_1(:,:,k) = real( corr_cholesky_mtx_1_dp(:,:,k), kind = core_rknd )
             corr_cholesky_mtx_2(:,:,k) = real( corr_cholesky_mtx_2_dp(:,:,k), kind = core_rknd )
          endif

       endif

    enddo  ! Setup PDF parameters loop: k = 2, nz, 1

    ! Boundary condition for the variance (overall) of a hydrometeor, <hm'^2>,
    ! on thermodynamic grid levels at the lowest thermodynamic grid level, k = 1
    ! (which is below the model lower boundary).
    hydrometp2_zt(1,:) = hydrometp2_zt(2,:)

    ! Interpolate the overall variance of a hydrometeor, <hm'^2>, to its home on
    ! momentum grid levels.
    do i = 1, hydromet_dim, 1
       hydrometp2(:,i)  = zt2zm( hydrometp2_zt(:,i) )
       hydrometp2(nz,i) = zero
    enddo


    return

  end subroutine setup_pdf_parameters

  !=============================================================================
  subroutine component_means_hydromet( nz, hydromet, rho, rc1, rc2, &
                                       mixt_frac, l_stats_samp, &
                                       hm1, hm2 )

    ! Description:
    ! The values of grid-level mean hydrometeor fields, <hm>, (for example,
    ! grid-level mean rain water mixing ratio, <r_r>, and grid-level mean rain
    ! drop concentration, <N_r>) are solved as part of the predictive equation
    ! set, based on the microphysics scheme.  However, CLUBB has a two component
    ! PDF.  The grid-level means of all hydrometeors must be subdivided into
    ! component means for each PDF component.  The equation relating the overall
    ! mean to the component means (for any hydrometeor, hm) is:
    !
    ! <hm> = a * hm1 + (1-a) * hm2;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), hm1
    ! is the mean of the hydrometeor in PDF component 1, and hm2 is the mean of
    ! the hydrometeor in PDF component 2.  This equation can be rewritten as:
    !
    ! <hm> = hm1 * ( a + (1-a) * hm2/hm1 ).
    !
    ! One way to solve for a component mean is to relate the ratio hm2/hm1 to
    ! other factors.  For now, this ratio based on other factors will be called
    ! hm2_hm1_ratio.  This ratio is entered into the above equation, allowing
    ! the equation to be solved for hm1:
    !
    ! hm1 = <hm> / ( a + (1-a) * hm2_hm1_ratio ).
    !
    ! Once hm1 has been solved for, hm2 can be solved by:
    !
    ! hm2 = ( <hm> - a * hm1 ) / (1-a).
    !
    ! At a grid level that is at least mostly cloudy, the simplest way to handle
    ! the ratio hm2/hm1 is to set it equal to the ratio rc2/rc1, where rc1 is
    ! the mean cloud water mixing ratio in PDF component 1 and rc2 is the mean
    ! cloud water mixing ratio in PDF component 2.  However, a precipitating
    ! hydrometeor sediments, falling from higher altitudes downwards.  The
    ! values of cloud water mixing ratio at a given grid level are not
    ! necessarily indicative of the amount of cloud water at higher levels.  A
    ! precipitating hydrometeor may have been already produced from cloud water
    ! at a higher altitude (vertical level) and fallen downwards to the given
    ! grid level.  Additionally, using grid-level cloud water mixing ratio
    ! especially does not work for a precipitating hydrometeor below cloud base
    ! (near the ground).
    !
    ! However, an alternative to component cloud water mixing ratio is component
    ! liquid water path.  Liquid water path accounts for the cloud water mixing
    ! ratio at the given grid level and at all grid levels higher in altitude.
    !
    ! In a stratocumulus case, the cloud water is spread out over all or almost
    ! all of the horizontal domain over a group of vertical levels.  At a given
    ! vertical level, the component mean cloud water mixing ratios should be
    ! almost equal, although usually slightly larger in the component with the
    ! larger component mean extended liquid water mixing ratio, s.  Likewise,
    ! the component liquid water paths should be nearly equal, with one
    ! component having a slightly larger liquid water path than the other
    ! component.
    !
    ! In a case of cumulus rising into stratocumulus, the upper portion of the
    ! cloudy domain will be very similar to the stratocumulus case described
    ! above, with similar cloud water mixing ratio and liquid water path
    ! results.  However, below the base of the stratocumulus clouds, where the
    ! cumulus clouds are found, the horizontal domain at each vertical level is
    ! only partially cloudy.  At these levels, any precipitating hydrometeor
    ! that was produced in the stratocumulus clouds above and fallen downwards
    ! is evaporating in the clear-air portions, while not evaporating in the
    ! cloudy portions.  Additionally, new amounts of a hydrometeor are being
    ! produced in the cloudy portions.  The amount of a hydrometeor in the
    ! cloudy portions becomes significantly larger than the amount of a
    ! hydrometeor in the clear portions.  The partially cloudy levels usually
    ! have a PDF where one component is significantly more saturated than the
    ! other component.  By the time the cloud base of the cumulus clouds is
    ! reached, the liquid water path for one PDF component should be
    ! significantly greater than the liquid water path for the other PDF
    ! component.
    !
    ! In a cumulus case, the horizontal domain at each level is usually partly
    ! cloudy.  Throughout the entire vertical domain, at every vertical level,
    ! one component usually is much more saturated than the other component.
    ! The liquid water path for one component is much greater than the liquid
    ! water path in the other component.  Likewise, a precipitating hydrometeor
    ! that is formed in cloud and falls preferentially through cloud will have
    ! large values in a portion of the horizontal domain and very small or 0
    ! values over the rest of the horizontal domain.
    !
    ! In order to estimate the amount of a hydrometeor in each PDF component,
    ! the ratio hm2/hm1 is going to be set equal to the ratio LWP2/LWP1, where
    ! LWP1 is the liquid water path in PDF component 1 and LWP2 is the liquid
    ! water path in PDF component 2.  LWP1 will be computed by taking the
    ! vertical integral of cloud water (see equation below) through the 1st PDF
    ! component from the given vertical level all the way to the top of the
    ! model.  LWP2 will be computed in the same manner.   It should be noted
    ! that this method makes the poor assumption that PDF component 1 always
    ! overlaps PDF component 1 between vertical levels, and likewise for PDF
    ! component 2.
    !
    ! Total liquid water path, LWP, is given by the following equation:
    !
    ! LWP(z) = INT(z:z_top) rho_a <r_c> dz';
    !
    ! where z is the altitude of the vertical level for which LWP is desired,
    ! z_top is the altitude at the top of the model domain, and z' is the
    ! dummy variable of integration.  Mean cloud water mixing ratio can be
    ! written as:
    !
    ! <r_c> = a * rc1 + (1-a) * rc2.
    !
    ! The equation for liquid water path is rewritten as:
    !
    ! LWP(z) = INT(z:z_top) rho_a ( a rc1 + (1-a) rc2 ) dz'; or
    !
    ! LWP(z) = INT(z:z_top) a rho_a rc1 dz'
    !          + INT(z:z_top) (1-a) rho_a rc2 dz'.
    !
    ! This can be rewritten as:
    !
    ! LWP(z) = LWP1(z) + LWP2(z);
    !
    ! where:
    !
    ! LWP1(z) = INT(z:z_top) a rho_a rc1 dz'; and
    ! LWP2(z) = INT(z:z_top) (1-a) rho_a rc2 dz'.
    !
    ! The trapezoidal rule will be used to numerically integrate for LWP1
    ! and LWP2.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one,      & ! Constant(s)
        one_half, &
        zero

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use parameters_microphys, only: &
        hydromet_tol  ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type, only: &
        stat_update_var  ! Procedure(s)

    use stats_variables, only : &
        iLWP1, & ! Variable(s)
        iLWP2, &
        zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho,       & ! Air density                                        [kg/m^3]
      rc1,       & ! Mean cloud water mixing ratio (1st PDF component)  [kg/kg]
      rc2,       & ! Mean cloud water mixing ratio (2nd PDF component)  [kg/kg]
      mixt_frac    ! Mixture fraction                                   [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hm1, & ! Mean of hydrometeor (1st PDF component)          [units vary]
      hm2    ! Mean of hydrometeor (2nd PDF component)          [units vary]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: &
      LWP1, & ! Liquid water path (1st PDF component) on thermo. levs.  [kg/m^2]
      LWP2    ! Liquid water path (2nd PDF component) on thermo. levs.  [kg/m^2]

    integer :: k, i  ! Array index

    real( kind = core_rknd ), parameter :: &
      LWP_tol = 5.0e-7_core_rknd  ! Tolerance value for component LWP


    !!! Compute component liquid water paths using trapezoidal rule for
    !!! numerical integration.

    ! At the uppermost thermodynamic level (k = nz), use the trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k = nz,
    ! integrand_b is the integrand at momentum level k = nz (model upper
    ! boundary), and delta_z = zm(nz) - zt(nz).  At the upper boundary, r_c is
    ! set to 0, and the form of the trapezoidal rule is simply:
    !
    ! 0.5 * integrand_a * delta_z.

    ! Liquid water path in PDF component 1.
    LWP1(nz) &
    = one_half * mixt_frac(nz) * rho(nz) * rc1(nz) * ( gr%zm(nz) - gr%zt(nz) )

    ! Liquid water path in PDF component 2.
    LWP2(nz) &
    = one_half * ( one - mixt_frac(nz) ) * rho(nz) * rc2(nz) &
      * ( gr%zm(nz) - gr%zt(nz) )

    ! At all other thermodynamic levels, compute liquid water path using the
    ! trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k, integrand_b
    ! is the integrand at thermodynamic level k+1, and
    ! delta_z = zt(k+1) - zt(k), or 1/invrs_dzm(k).  The total for the segment
    ! is added to the sum total of all higher vertical segments to compute the
    ! total vertical integral.
    do k = nz-1, 1, -1

       ! Liquid water path in PDF component 1.
       LWP1(k) &
       = LWP1(k+1) &
         + one_half * ( mixt_frac(k+1) * rho(k+1) * rc1(k+1) &
                        + mixt_frac(k) * rho(k) * rc1(k) ) / gr%invrs_dzm(k)

       ! Liquid water path in PDF component 2.
       LWP2(k) &
       = LWP2(k+1) &
         + one_half * ( ( one - mixt_frac(k+1) ) * rho(k+1) * rc2(k+1) &
                        + ( one - mixt_frac(k) ) * rho(k) * rc2(k) ) &
           / gr%invrs_dzm(k)

    enddo ! k = nz-1, 1, -1


    !!! Find hm1 and hm2 based on the ratio of LWP2/LWP1, such that:
    !!! hm2/hm1 ( = rr2/rr1 = nr2/nr1, etc. ) = LWP2/LWP1.
    do i = 1, hydromet_dim, 1

       do k = 1, nz, 1

          !!! Calculate the component means for the hydrometeor.
          if ( hydromet(k,i) > hydromet_tol(i) ) then

             if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

                ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
                !
                ! The hydrometeor is found at this level, yet there is no cloud
                ! at or above the current level.  This is usually due to a
                ! numerical artifact.  For example, the hydrometeor is diffused
                ! above cloud top.  Simply set each component mean equal to the
                ! overall mean.
                hm1(k,i) = hydromet(k,i)
                hm2(k,i) = hydromet(k,i)

             elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

                ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
                ! insignificant amount).
                !
                ! The hydrometeor is found at this level, and all cloud water at
                ! or above this level is found in the 1st PDF component.  All of
                ! the hydrometeor is found in the 1st PDF component.
                hm1(k,i) = hydromet(k,i) / mixt_frac(k)
                hm2(k,i) = zero

             elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

                ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
                ! insignificant amount).
                !
                ! The hydrometeor is found at this level, and all cloud water at
                ! or above this level is found in the 2nd PDF component.  All of
                ! the hydrometeor is found in the 2nd PDF component.
                hm1(k,i) = zero
                hm2(k,i) = hydromet(k,i) / ( one - mixt_frac(k) )

             else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

                ! Both LWP1 and LWP2 are (significantly) greater than 0.
                !
                ! The hydrometeor is found at this level, and there is
                ! sufficient cloud water at or above this level in both PDF
                ! components to find the hydrometeor in both PDF components.
                ! Delegate the hydrometeor between the 1st and 2nd PDF
                ! components according to the above equations.
                hm1(k,i) &
                = hydromet(k,i) &
                  / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

                hm2(k,i) &
                = ( hydromet(k,i) - mixt_frac(k) * hm1(k,i) ) &
                  / ( one - mixt_frac(k) )

                if ( hm1(k,i) <= hydromet_tol(i) ) then

                   ! The mean value of the hydrometeor within the 1st PDF
                   ! component is below the tolerance value for the hydrometeor.
                   ! It is considered to have a value of 0.  All the the
                   ! hydrometeor is found within the 2nd PDF component.
                   hm1(k,i) = zero
                   hm2(k,i) = hydromet(k,i) / ( one - mixt_frac(k) )

                elseif ( hm2(k,i) <= hydromet_tol(i) ) then

                   ! The mean value of the hydrometeor within the 2nd PDF
                   ! component is below the tolerance value for the hydrometeor.
                   ! It is considered to have a value of 0.  All the the
                   ! hydrometeor is found within the 1st PDF component.
                   hm1(k,i) = hydromet(k,i) / mixt_frac(k)
                   hm2(k,i) = zero

                endif

             endif


          else ! hydromet(k,i) <= hydromet_tol(i)

             ! The overall hydrometeor is either 0 or below tolerance value (any
             ! postive value is considered to be a numerical artifact).  Simply
             ! set each pdf component mean equal to 0.
             hm1(k,i) = zero
             hm2(k,i) = zero

          endif

       enddo ! k = 1, nz, 1

    enddo ! i = 1, hydromet_dim, 1


    ! Statistics
    if ( l_stats_samp ) then

       if ( iLWP1 > 0 ) then
          ! Liquid water path in PDF component 1.
          call stat_update_var( iLWP1, LWP1, zt )
       endif

       if ( iLWP2 > 0 ) then
          ! Liquid water path in PDF component 2.
          call stat_update_var( iLWP2, LWP2, zt )
       endif
       
    endif


    return

  end subroutine component_means_hydromet

  !=============================================================================
  subroutine precip_fraction( nz, hydromet, hm1, hm2, &
                              cloud_frac, cloud_frac1, mixt_frac, &
                              ice_supersat_frac, &
                              precip_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! Determines (overall) precipitation fraction over the horizontal domain, as
    ! well as the precipitation fraction within each PDF component, at every
    ! vertical grid level.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        cloud_frac_min

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use parameters_microphys, only: &
        hydromet_tol  ! Variable(s)

    use array_index, only: &
        l_mix_rat_hm  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet, & ! Mean of hydrometeor, hm (overall)           [units vary]
      hm1,      & ! Mean of hydrometeor (1st PDF component)     [units vary]
      hm2         ! Mean of hydrometeor (2nd PDF component)     [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,  &     ! Cloud fraction (overall)                     [-] 
      cloud_frac1, &     ! Cloud fraction (1st PDF component)           [-]
      mixt_frac, &       ! Mixture fraction                             [-]
      ice_supersat_frac  ! Ice cloud fraction                           [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      weighted_pfrac1    ! Product of mixt_frac and precip_frac_1       [-]

    real( kind = core_rknd ) :: &
      r_tot_hm_1, & ! Mean total hydromet mixing ratio (1st PDF comp.)  [kg/kg]
      r_tot_hm_2, & ! Mean total hydromet mixing ratio (2nd PDF comp.)  [kg/kg]
      N_tot_hm_1, & ! Mean total hydromet concentration (1st PDF comp.) [num/kg]
      N_tot_hm_2    ! Mean total hydromet concentration (2nd PDF comp.) [num/kg]

    real( kind = core_rknd ), parameter :: &
      precip_frac_tol = cloud_frac_min  ! Minimum precip. frac.         [-]
    
    ! "Maximum allowable" hydrometeor mixing ratio in-precip component mean.
    real( kind = core_rknd ), parameter :: &
      max_hm_ip_comp_mean = 0.0025_core_rknd  ! [kg/kg]

    integer :: &
      k, i   ! Loop indices



    ! Initialize the precipitation fraction variables (precip_frac,
    ! precip_frac_1, and precip_frac_2) to 0.
    precip_frac   = zero
    precip_frac_1 = zero
    precip_frac_2 = zero

    !!! Find overall precipitation fraction.
    do k = nz, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       if ( k < nz ) then
          precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )
       else  ! k = nz
          precip_frac(k) = cloud_frac(k)
       endif

       if ( any( hydromet(k,:) > hydromet_tol(:) ) &
            .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find any hydrometeor at this grid level, but
          ! no cloud at or above this grid level, set precipitation fraction to
          ! a minimum threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( all( hydromet(k,:) <= hydromet_tol(:) ) &
                .and. precip_frac(k) < precip_frac_tol ) then

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There are not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo ! Overall precipitation fraction loop: k = nz, 1, -1.

    !!! Account for ice cloud fraction
    do k = nz, 1, -1
      precip_frac(k) = max( precip_frac(k), ice_supersat_frac(k) )
    enddo


    !!! Find precipitation fraction within each PDF component.
    !
    ! The overall precipitation fraction, f_p, is given by the equation:
    !
    ! f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    !
    ! where a is the mixture fraction (weight of PDF component 1), f_p(1) is
    ! the precipitation fraction within PDF component 1, and f_p(2) is the
    ! precipitation fraction within PDF component 2.  Overall precipitation
    ! fraction is found according the method above, and mixture fraction is
    ! already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    ! values for f_p(1) and f_p(2) must satisfy the above equation.
    if ( .false. ) then

       !!! Find precipitation fraction within PDF component 1.
       ! The method used to find overall precipitation fraction will also be to
       ! find precipitation fraction within PDF component 1.  In order to do so,
       ! it is assumed (poorly) that PDF component 1 overlaps PDF component 1 at
       ! every vertical level in the vertical profile.
       do k = nz, 1, -1

          ! The weighted precipitation fraction (PDF component 1) is the
          ! greatest value of the product of mixture fraction and cloud fraction
          ! (PDF component 1) at or above a vertical level.
          if ( k < nz ) then
             weighted_pfrac1(k) = max( weighted_pfrac1(k+1), &
                                       mixt_frac(k) * cloud_frac1(k) )
          else  ! k = nz
             weighted_pfrac1(k) = mixt_frac(k) * cloud_frac1(k)
          endif

          precip_frac_1(k) = weighted_pfrac1(k) / mixt_frac(k)

          ! Special cases for precip_frac_1.
          if ( precip_frac_1(k) > one ) then

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  For example, the mixture fraction at level k+1
             ! is 0.10 and the cloud_frac1 at level k+1 is 1, resulting in a
             ! weighted_pfrac1 of 0.10.  This product is greater than the
             ! product of mixt_frac and cloud_frac1 at level k.  The mixture
             ! fraction at level k is 0.05, resulting in a precip_frac_1 of 2.
             ! The value of precip_frac_1 is limited at 1.  The leftover
             ! precipitation fraction (a result of the decreasing weight of PDF
             ! component 1 between the levels) is applied to PDF component 2.
             precip_frac_1(k) = one

          elseif ( any( hm1(k,:) > hydromet_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 1st PDF
             ! component at this grid level, but no cloud in the 1st PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 1st PDF component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( all( hm1(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 1st PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There are not any
             ! hydrometeor species found in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif

       enddo ! Precipitation fraction (1st PDF component) loop: k = nz, 1, -1.


       !!! Find precipitation fraction within PDF component 2.
       ! The equation for precipitation fraction within PDF component 2 is:
       !
       ! f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
       !
       ! given the overall precipitation fraction, f_p (calculated above), the
       ! precipitation fraction within PDF component 1, f_p(1) (calculated
       ! above), and mixture fraction, a.  Any leftover precipitation fraction
       ! from precip_frac_1 will be included in this calculation of
       ! precip_frac_2.
       do k = 1, nz, 1

          precip_frac_2(k) &
          = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
            / ( one - mixt_frac(k) )

          ! Special cases for precip_frac_2.
          if ( precip_frac_2(k) > one ) then

             ! Again, it is possible for precip_frac_2 to be greater than 1.
             ! For example, the mixture fraction at level k+1 is 0.10 and the
             ! cloud_frac1 at level k+1 is 1, resulting in a weighted_pfrac1 of
             ! 0.10.  This product is greater than the product of mixt_frac and
             ! cloud_frac1 at level k.  Additionally, precip_frac (overall) is 1
             ! for level k.  The mixture fraction at level k is 0.5, resulting
             ! in a precip_frac_1 of 0.2.  Using the above equation,
             ! precip_frac_2 is calculated to be 1.8.  The value of
             ! precip_frac_2 is limited at 1.  The leftover precipitation
             ! fraction (as a result of the increasing weight of component 1
             ! between the levels) is applied to PDF component 1.
             precip_frac_2(k) = one

             ! Recalculate the precipitation fraction in PDF component 1.
             precip_frac_1(k) &
             = ( precip_frac(k) - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
               / mixt_frac(k)

             ! Double check for errors in PDF component 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             elseif ( any( hm1(k,:) > hydromet_tol(:) ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = precip_frac_tol
             elseif ( all( hm1(k,:) <= hydromet_tol(:) ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = zero
             endif

          elseif ( any( hm2(k,:) > hydromet_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 2nd PDF
             ! component at this grid level, but no cloud in the 2nd PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 2nd PDF component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( all( hm2(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 2nd PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There are not any
             ! hydrometeor species found in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif

       enddo ! Precipitation fraction (2nd PDF component) loop: k = 1, nz, 1.


    else  ! .true.

       ! Precipitation fraction in each PDF component is based on the mean total
       ! hydrometeor mixing ratio in each PDF component, where total hydrometeor
       ! mixing ratio, r_Thm, is the sum of all precipitating hydrometeor
       ! species mixing ratios (which doesn't include cloud water), such that:
       !
       ! r_Thm = r_r + r_i + r_s + r_g;
       !
       ! where r_r is rain water mixing ratio, r_i is ice mixing ratio, r_s is
       ! snow mixing ratio, and r_g is graupel mixing ratio.
       !
       ! Precipitation fraction in each PDF component is based on the ratio:
       !
       ! r_Thm_1/f_p(1) = r_Thm_2/f_p(2);
       !
       ! where r_Thm_1 is mean total hydrometeor mixing ratio is the 1st PDF
       ! component and r_Thm_2 is mean total hydrometeor mixing ratio in the 2nd
       ! PDF component.  The equation can be rewritten as:
       !
       ! f_p(2)/f_p(1) = r_Thm_2/r_Thm_1.
       !
       ! Since overall precipitation fraction is given by the equation:
       !
       ! f_p = a f_p(1) + (1-a) f_p(2);
       !
       ! it can be rewritten as:
       !
       ! f_p = f_p(1) ( a + (1-a) f_p(2)/f_p(1) ).
       !
       ! Substituting the ratio r_Thm_2/r_Thm_1 for the ratio f_p(2)/f_p(1), the
       ! above equation can be solved for f_p(1):
       !
       ! f_p(1) = f_p / ( a + (1-a) r_Thm_2/r_Thm_1 ).
       !
       ! Then, f_p(2) can be solved for according to the equation:
       !
       ! f_p(2) = ( f_p - a f_p(1) ) / (1-a).
       !
       ! In the event where hydrometeor concentrations are found at a given
       ! vertical level, but not hydrometeor mixing ratios (due to numerical
       ! artifacts), the mean total hydrometeor concentrations in each PDF
       ! component will be used in place of mean total hydrometeor mixing ratios
       ! in the above equations to solve for component precipitation fractions.
       do k = 1, nz, 1

          if ( all( hm1(k,:) <= hydromet_tol(:) ) &
               .and. all( hm2(k,:) <= hydromet_tol(:) ) ) then

             ! There are no hydrometeors found in each PDF component.
             ! Precipitation fraction within each component is set to 0.
             precip_frac_1(k) = zero
             precip_frac_2(k) = zero

          elseif ( any( hm1(k,:) > hydromet_tol(:) ) &
                   .and. all( hm2(k,:) <= hydromet_tol(:) ) ) then

             ! All the hydrometeors are found within the 1st PDF component.
             precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
             precip_frac_2(k) = zero

          elseif ( any( hm2(k,:) > hydromet_tol(:) ) &
                   .and. all( hm1(k,:) <= hydromet_tol(:) ) ) then

             ! All the hydrometeors are found within the 2nd PDF component.
             precip_frac_1(k) = zero
             precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )

          else

             ! any( hm1(k,:) > hydromet_tol(:) )
             ! AND any( hm2(k,:) > hydromet_tol(:) )

             ! Hydrometeors are found within both PDF components.
             r_tot_hm_1 = zero
             r_tot_hm_2 = zero
             N_tot_hm_1 = zero
             N_tot_hm_2 = zero
             do i = 1, hydromet_dim, 1

                if ( l_mix_rat_hm(i) ) then

                   ! The hydrometeor is a mixing ratio.
                   ! Find total hydrometeor mixing ratio in each PDF component.
                   if ( hm1(k,i) > hydromet_tol(i) ) then
                      r_tot_hm_1 = r_tot_hm_1 + hm1(k,i)
                   endif
                   if ( hm2(k,i) > hydromet_tol(i) ) then
                      r_tot_hm_2 = r_tot_hm_2 + hm2(k,i)
                   endif

                else ! l_mix_rat_hm(i) is false

                   ! The hydrometeor is a concentration.
                   ! Find total hydrometeor concentration in each PDF component.
                   if ( hm1(k,i) > hydromet_tol(i) ) then
                      N_tot_hm_1 = N_tot_hm_1 + hm1(k,i)
                   endif
                   if ( hm2(k,i) > hydromet_tol(i) ) then
                      N_tot_hm_2 = N_tot_hm_2 + hm2(k,i)
                   endif

                endif ! l_mix_rat_hm(i)

             enddo ! i = 1, hydromet_dim, 1

             !!! Find precipitation fraction within PDF component 1.
             if ( r_tot_hm_1 > zero ) then
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) &
                      + ( one - mixt_frac(k) ) * r_tot_hm_2/r_tot_hm_1 )
             else ! N_tot_hm_1 > zero 
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) &
                      + ( one - mixt_frac(k) ) * N_tot_hm_2/N_tot_hm_1 )
             endif

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  The value of precip_frac_1 is limited at 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             endif

             !!! Find precipitation fraction within PDF component 2.
             precip_frac_2(k) &
             = ( precip_frac(k) - mixt_frac(k) *  precip_frac_1(k) ) &
               / ( one - mixt_frac(k) )

             ! Using the above method, it is possible for precip_frac_2 to be
             ! greater than 1.  The value of precip_frac_2 is limited at 1.
             if ( precip_frac_2(k) > one ) then

                precip_frac_2(k) = one

                ! Recalculate the precipitation fraction in PDF component 1.
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

             endif

          endif


          ! Special cases for PDF component 1.
          if ( any( hm1(k,:) > hydromet_tol(:) ) &
               .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 1st PDF
             ! component at this grid level, but no cloud in the 1st PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 1st PDF component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( all( hm1(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 1st PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There is not any
             ! hydrometeor species found in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif


          ! Special cases for PDF component 2.
          if ( any( hm2(k,:) > hydromet_tol(:) ) &
               .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 2nd PDF
             ! component at this grid level, but no cloud in the 2nd PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 2nd PDF component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( all( hm2(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 2nd PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There is not any
             ! hydrometeor species found in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif


       enddo ! Component precipitation fraction loop: k = 1, nz, 1.


    endif ! Select component precipitation fraction method.


    ! Increase Precipiation Fraction under special conditions.
    !
    ! There are scenarios that sometimes occur that require precipitation
    ! fraction to be boosted.  Precipitation fraction is calculated from cloud
    ! fraction and ice supersaturation fraction.  For numerical reasons, CLUBB's
    ! PDF may become entirely subsaturated with respect to liquid and ice,
    ! resulting in both a cloud fraction of 0 and an ice supersaturation
    ! fraction of 0.  When this happens, precipitation fraction drops to 0 when
    ! there aren't any hydrometeors present at that grid level, or to
    ! precip_frac_tol when there is at least one hydrometeor present at that
    ! grid level.  However, sometimes there are large values of hydrometeors
    ! found at that grid level.  When this occurs, the PDF component in-precip
    ! mean of a hydrometeor can become ridiculously large.  This is because the
    ! ith PDF component in-precip mean of a hydrometeor, mu_hm_i,  is given by
    ! the equation:
    !
    ! mu_hm_i = hmi / precip_frac_i;
    !
    ! where hmi is the overall ith PDF component mean of the hydrometeor, and
    ! precip_frac_i is the ith PDF component precipitation fraction.  When
    ! precip_frac_i has a value of precip_frac_tol and hmi is large, mu_hm_i can
    ! be huge.  This can cause enormous microphysical process rates and result
    ! in numerical instability.  It is also very inaccurate.
    !
    ! In order to limit this problem, the ith PDF component precipitation
    ! fraction is increased in order to decrease mu_hm_i.  First, an "upper
    ! limit" is set for mu_hm_i when the hydrometeor is a mixing ratio.  This is
    ! called max_hm_ip_comp_mean.  At every vertical level and for every
    ! hydrometeor mixing ratio, a check is made to try to prevent mu_hm_i from
    ! exceeding the "upper limit".  The check is:
    ! hmi / precip_frac_i ( which = mu_hm_i ) > max_hm_ip_comp_mean, which can
    ! be rewritten:  hmi > precip_frac_i * max_hm_ip_comp_mean.  When this
    ! occurs, precip_frac_i is increased to hmi/max_hm_ip_comp_mean.  Of course,
    ! precip_frac_i is not allowed to exceed 1, so when hmi is already greater
    ! than max_hm_ip_comp_mean, mu_hm_i will also have to be greater than
    ! max_hm_ip_comp_mean.  However, the value of mu_hm_i is still reduced when
    ! compared to what it would have been using precip_frac_tol.  In the event
    ! that multiple hydrometeor mixing ratios violate the check, the code is set
    ! up so that precip_frac_i is increased based on the highest hmi.
    do k = 1, nz, 1

       do i = 1, hydromet_dim, 1

          if ( l_mix_rat_hm(i) ) then

             ! The hydrometeor is a mixing ratio.

             if ( hm1(k,i) > precip_frac_1(k) * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 1st PDF component.
                precip_frac_1(k) = min( hm1(k,i)/max_hm_ip_comp_mean, one )

                ! Recalculate overall precipitation fraction.
                precip_frac(k) = mixt_frac(k) * precip_frac_1(k) &
                                 + ( one - mixt_frac(k) ) * precip_frac_2(k)

             endif ! mu_hm_1 = hm1/precip_frac_1 > max_hm_ip_comp_mean

             if ( hm2(k,i) > precip_frac_2(k) * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 2nd PDF component.
                precip_frac_2(k) = min( hm2(k,i)/max_hm_ip_comp_mean, one )

                ! Recalculate overall precipitation fraction.
                precip_frac(k) = mixt_frac(k) * precip_frac_1(k) &
                                 + ( one - mixt_frac(k) ) * precip_frac_2(k)

             endif ! mu_hm_2 = hm2/precip_frac_2 > max_hm_ip_comp_mean

          endif ! l_mix_rat_hm(i)

       enddo ! i = 1, hydromet_dim, 1

    enddo ! k = 1, nz, 1


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine compute_mean_stdev( Ncnm, rc1, rc2, &                      ! Intent(in)
                                 cloud_frac1, cloud_frac2, &            ! Intent(in)
                                 hm1, hm2, &                            ! Intent(in)
                                 precip_frac_1, precip_frac_2, &        ! Intent(in)
                                 sigma2_on_mu2_ip_array_cloud, &        ! Intent(in)
                                 sigma2_on_mu2_ip_array_below, &        ! Intent(in)
                                 pdf_params, d_variables, &             ! Intent(in)
                                 mu_x_1, mu_x_2, sigma_x_1, sigma_x_2 ) ! Intent(out)
       
    ! Description:
    ! Calculates the means and standard deviations (for each PDF component) of
    ! s, t, w, Ncn, and the precipitating hydrometeors.  For the precipitating
    ! hydrometeors, the component means and standard deviations are in-precip. 

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,  & ! Constant(s)
        zero

    use parameters_microphys, only: &
        Ncnp2_on_Ncnm2, & ! Variable(s)
        hydromet_tol

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use corr_matrix_module, only: &
        iiPDF_s_mellor, & ! Variable(s)
        iiPDF_t_mellor, &
        iiPDF_w,        &
        iiPDF_Ncn

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: d_variables ! Number of PDF variables

    real( kind = core_rknd ), intent(in) :: &
      Ncnm,          & ! Mean cloud nuclei concentration                [num/kg]
      rc1,           & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc2,           & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac1,   & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac2,   & ! Cloud fraction (2nd PDF component)                  [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array_cloud, & ! Prescribed ratio array: cloudy levs. [-]
      sigma2_on_mu2_ip_array_below    ! Prescribed ratio array: clear levs.  [-]

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)    [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)    [units vary]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    ! Note:  This code assumes to be these arrays in the same order as the
    ! correlation arrays, etc., which is determined by the iiPDF indices.
    ! The order should be as follows:  s, t, w, Ncn, <precip. hydrometeors>
    ! (indices increasing from left to right).
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    ! Local Variables
    integer :: ivar ! Loop iterator


    !!! Enter the PDF parameters.

    !!! Means.

    ! Mean of vertical velocity, w, in PDF component 1.
    mu_x_1(iiPDF_w) = pdf_params%w1

    ! Mean of vertical velocity, w, in PDF component 2.
    mu_x_2(iiPDF_w) = pdf_params%w2

    ! Mean of extended liquid water mixing ratio, s, in PDF component 1.
    mu_x_1(iiPDF_s_mellor) = pdf_params%s1

    ! Mean of extended liquid water mixing ratio, s, in PDF component 2.
    mu_x_2(iiPDF_s_mellor) = pdf_params%s2

    ! Mean of t in PDF component 1.
    ! Set the component mean values of t to 0.
    ! The component mean values of t are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_x_1(iiPDF_t_mellor) = zero

    ! Mean of t in PDF component 2.
    ! Set the component mean values of t to 0.
    ! The component mean values of t are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_x_2(iiPDF_t_mellor) = zero

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 1.
    mu_x_1(iiPDF_Ncn) = Ncnm

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 2.
    mu_x_2(iiPDF_Ncn) = Ncnm

    ! Mean of the hydrometeor species
    do ivar = iiPDF_Ncn+1, d_variables

       ! Mean of hydrometeor, hm, in PDF component 1.
       mu_x_1(ivar) &
       = component_mean_hm_ip( hm1(pdf2hydromet_idx(ivar)), precip_frac_1, &
                               hydromet_tol(pdf2hydromet_idx(ivar)) )

       ! Mean of hydrometeor, hm, in PDF component 2.
       mu_x_2(ivar) &
       = component_mean_hm_ip( hm2(pdf2hydromet_idx(ivar)), precip_frac_2, &
                               hydromet_tol(pdf2hydromet_idx(ivar)) )

    enddo


    !!! Standard deviations.

    ! Standard deviation of vertical velocity, w, in PDF component 1.
    sigma_x_1(iiPDF_w) = sqrt( pdf_params%varnce_w1 )

    ! Standard deviation of vertical velocity, w, in PDF component 2.
    sigma_x_2(iiPDF_w) = sqrt( pdf_params%varnce_w2 )

    ! Standard deviation of extended liquid water mixing ratio, s,
    ! in PDF component 1.
    sigma_x_1(iiPDF_s_mellor) = pdf_params%stdev_s1

    ! Standard deviation of extended liquid water mixing ratio, s,
    ! in PDF component 2.
    sigma_x_2(iiPDF_s_mellor) = pdf_params%stdev_s2

    ! Standard deviation of t in PDF component 1.
    sigma_x_1(iiPDF_t_mellor) = pdf_params%stdev_t1

    ! Standard deviation of t in PDF component 2.
    sigma_x_2(iiPDF_t_mellor) = pdf_params%stdev_t2

    ! Standard deviation of simplified cloud nuclei concentration, Ncn,
    ! in PDF component 1.
    sigma_x_1(iiPDF_Ncn) &
    = component_stdev_hm_ip( mu_x_1(iiPDF_Ncn), rc1, one, &
                             Ncnp2_on_Ncnm2, Ncnp2_on_Ncnm2 )

    ! Standard deviation of simplified cloud nuclei concentration, Ncn,
    ! in PDF component 2.
    sigma_x_2(iiPDF_Ncn) &
    = component_stdev_hm_ip( mu_x_2(iiPDF_Ncn), rc2, one, &
                             Ncnp2_on_Ncnm2, Ncnp2_on_Ncnm2 )

    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    ! Brian Griffin; February 3, 2007.

    do ivar = iiPDF_Ncn+1, d_variables

       ! Standard deviation of hydrometeor, hm, in PDF component 1.
       sigma_x_1(ivar) &
       =  component_stdev_hm_ip( mu_x_1(ivar), &
                                 rc1, cloud_frac1, &
                                 sigma2_on_mu2_ip_array_cloud(ivar), &
                                 sigma2_on_mu2_ip_array_below(ivar) )

       ! Standard deviation of hydrometeor, hm, in PDF component 2.
       sigma_x_2(ivar) &
       =  component_stdev_hm_ip( mu_x_2(ivar), &
                                 rc2, cloud_frac2, &
                                 sigma2_on_mu2_ip_array_cloud(ivar), &
                                 sigma2_on_mu2_ip_array_below(ivar) )

    enddo


    return

  end subroutine compute_mean_stdev

  !=============================================================================
  subroutine compute_corr( wm_zt, rc1, rc2, cloud_frac1, &
                           cloud_frac2, wpsp, wpNcnp, &
                           stdev_w, mixt_frac, precip_frac_1, &
                           precip_frac_2, wphydrometp_zt, &
                           mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                           corr_array_cloud, corr_array_below, &
                           pdf_params, d_variables, &
                           corr_array_1, corr_array_2 )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        Ncn_tol,      &
        w_tol,        & ! [m/s]
        s_mellor_tol, & ! [kg/kg]
        one,          &
        zero

    use model_flags, only: &
        l_calc_w_corr, &
        l_use_modified_corr

    use diagnose_correlations_module, only: &
        calc_mean,        & ! Procedure(s)
        calc_w_corr

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use corr_matrix_module, only: &
        iiPDF_s_mellor, & ! Variable(s)
        iiPDF_t_mellor, &
        iiPDF_w,        &
        iiPDF_Ncn

    implicit none

    ! Input Variables
    integer, intent(in) :: d_variables ! Number of variables in the corr/mean/stdev arrays

    real( kind = core_rknd ), intent(in) :: &
      wm_zt,         & ! Mean vertical velocity, <w>, on thermo. levels    [m/s]
      rc1,           & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc2,           & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac1,   & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac2,   & ! Cloud fraction (2nd PDF component)                  [-]
      wpsp,          & ! Covariance of w and s                      [(m/s)kg/kg]
      wpNcnp,        & ! Covariance of w and N_cn (overall)       [(m/s) num/kg]
      stdev_w,       & ! Standard deviation of w                           [m/s]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      wphydrometp_zt    ! Covariance of w and hm interp. to t-levs.  [(m/s)u.v.]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean of x array (1st PDF component)            [units vary]
      mu_x_2,    & ! Mean of x array (2nd PDF component)            [units vary]
      sigma_x_1, & ! Standard deviation of x array (1st PDF comp.)  [units vary]
      sigma_x_2    ! Standard deviation of x array (2nd PDF comp.)  [units vary]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_cloud, & ! Prescribed correlation array in cloud        [-]
      corr_array_below    ! Prescribed correlation array below cloud     [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(out) :: &
      corr_array_1, & ! Correlation array 1st component (order: s,t,w,Ncn <hydrometeors>)   [-]
      corr_array_2    ! Correlation array 2nd component (order: s,t,w,Ncn <hydrometeors>)   [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      sigma_Ncn_1

    real( kind = core_rknd ), dimension(d_variables)  :: &
      corr_whm_1, & ! Correlation between w and hm (1st PDF component) ip    [-]
      corr_whm_2    ! Correlation between w and hm (2nd PDF component) ip    [-]

    real( kind = core_rknd ) :: &
      s_mellor_m,     & ! Mean of s_mellor                               [kg/kg]
      stdev_s_mellor, & ! Standard deviation of s_mellor                 [kg/kg]
      corr_ws,        & ! Correlation between w and s (overall)              [-]
      corr_wNcn         ! Correlation between w and Ncn (overall)            [-]

    logical :: &
      l_limit_corr_st    ! If true, we limit the correlation between s and t [-]

    integer :: ivar, jvar ! Loop iterators

    ! ---- Begin Code ----

    !!! Enter the PDF parameters.
    sigma_Ncn_1 = sigma_x_1(iiPDF_Ncn)

    !!! Correlations

    ! Initialize corr_whm_1 and corr_whm_2 arrays to 0.
    corr_whm_1 = zero
    corr_whm_2 = zero

    ! Calculate correlations involving w by first calculating total covariances
    ! involving w (<w'r_r'>, etc.) using the down-gradient approximation.
    if ( l_calc_w_corr ) then

       s_mellor_m &
       = calc_mean( pdf_params%mixt_frac, pdf_params%s1, pdf_params%s2 )

       stdev_s_mellor &
       = sqrt( pdf_params%mixt_frac &
               * ( ( pdf_params%s1 - s_mellor_m )**2 &
                   + pdf_params%stdev_s1**2 ) &
             + ( one - pdf_params%mixt_frac ) &
               * ( ( pdf_params%s2 - s_mellor_m )**2 &
                   + pdf_params%stdev_s2**2 ) &
             )

       corr_ws &
       = calc_w_corr( wpsp, stdev_w, stdev_s_mellor, w_tol, s_mellor_tol )

       corr_wNcn = calc_w_corr( wpNcnp, stdev_w, sigma_Ncn_1, w_tol, Ncn_tol )

       do jvar = iiPDF_Ncn+1, d_variables

          call calc_corr_whm( wm_zt, wphydrometp_zt(pdf2hydromet_idx(jvar)), &
                              mu_x_1(iiPDF_w), mu_x_2(iiPDF_w), &
                              mu_x_1(jvar), mu_x_2(jvar), &
                              sigma_x_1(iiPDF_w), sigma_x_2(iiPDF_w), &
                              sigma_x_1(jvar), sigma_x_2(jvar), &
                              mixt_frac, precip_frac_1, precip_frac_2, &
                              corr_whm_1(jvar), corr_whm_2(jvar) )

       enddo ! jvar = iiPDF_Ncn+1, d_variables


    endif

    ! In order to decompose the correlation matrix that we may or may not make
    ! (l_use_modified_corr), we must not have a perfect correlation between s
    ! and t. Thus, we impose a limitation.
    if ( l_use_modified_corr ) then
      l_limit_corr_st = .true.
    else
      l_limit_corr_st = .false.
    end if

    ! Initialize the correlation arrays
    corr_array_1 = zero
    corr_array_2 = zero

    !!! The corr_arrays are assumed to be lower triangular matrices
    ! Set diagonal elements to 1
    do ivar=1, d_variables
      corr_array_1(ivar, ivar) = one
      corr_array_2(ivar, ivar) = one
    end do


    !!! This code assumes the following order in the prescribed correlation arrays (iiPDF indices):
    !!! s, t, w, Ncn, <hydrometeors> (indices increasing from left to right)

    ! Correlation between s and t
    corr_array_1(iiPDF_t_mellor, iiPDF_s_mellor) &
    = component_corr_st( pdf_params%corr_st_1, rc1, cloud_frac1, &
                         corr_array_cloud(iiPDF_t_mellor, iiPDF_s_mellor), &
                         corr_array_below(iiPDF_t_mellor, iiPDF_s_mellor), &
                         l_limit_corr_st )

    corr_array_2(iiPDF_t_mellor, iiPDF_s_mellor) &
    = component_corr_st( pdf_params%corr_st_2, rc2, cloud_frac2, &
                         corr_array_cloud(iiPDF_t_mellor, iiPDF_s_mellor), &
                         corr_array_below(iiPDF_t_mellor, iiPDF_s_mellor), &
                         l_limit_corr_st )

    ! Correlation between s and w
    corr_array_1(iiPDF_w, iiPDF_s_mellor) &
    = component_corr_wx( corr_ws, rc1, cloud_frac1, &
                         corr_array_cloud(iiPDF_w, iiPDF_s_mellor), &
                         corr_array_below(iiPDF_w, iiPDF_s_mellor) )

    corr_array_2(iiPDF_w, iiPDF_s_mellor) &
    = component_corr_wx( corr_ws, rc2, cloud_frac2, &
                         corr_array_cloud(iiPDF_w, iiPDF_s_mellor), &
                         corr_array_below(iiPDF_w, iiPDF_s_mellor) )


    ! Correlation between s and Ncn
    corr_array_1(iiPDF_Ncn, iiPDF_s_mellor) &
    = component_corr_xhm_ip( rc1, one, &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_s_mellor), &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_s_mellor) )

    corr_array_2(iiPDF_Ncn, iiPDF_s_mellor) &
    = component_corr_xhm_ip( rc2, one, &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_s_mellor), &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_s_mellor) )

    ! Correlation between s and the hydrometeors
    ivar = iiPDF_s_mellor
    do jvar = iiPDF_Ncn+1, d_variables
       corr_array_1(jvar, ivar) &
       = component_corr_xhm_ip( rc1, cloud_frac1,&
                                corr_array_cloud(jvar, ivar), corr_array_below(jvar, ivar) )

       corr_array_2(jvar, ivar) &
       = component_corr_xhm_ip( rc2, cloud_frac2,&
                                corr_array_cloud(jvar, ivar), corr_array_below(jvar, ivar) )
    enddo

    ! Correlation between t and w
    corr_array_1(iiPDF_w, iiPDF_t_mellor) = zero
    corr_array_2(iiPDF_w, iiPDF_t_mellor) = zero

    ! Correlation between t and Ncn
    corr_array_1(iiPDF_Ncn, iiPDF_t_mellor) &
    = component_corr_xhm_ip( rc1, one, &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_t_mellor), &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_t_mellor) )

    corr_array_2(iiPDF_Ncn, iiPDF_t_mellor) &
    = component_corr_xhm_ip( rc2, one, &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_t_mellor), &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_t_mellor) )

    ! Correlation between t and the hydrometeors
    ivar = iiPDF_t_mellor
    do jvar = iiPDF_Ncn+1, d_variables

       if ( l_use_modified_corr ) then

          corr_array_1(jvar, ivar) &
          = component_corr_thm_ip( corr_array_1( iiPDF_t_mellor, iiPDF_s_mellor), &
                                   corr_array_1( jvar, iiPDF_s_mellor) )

          corr_array_2(jvar, ivar) &
          = component_corr_thm_ip( corr_array_2( iiPDF_t_mellor, iiPDF_s_mellor), &
                                   corr_array_2( jvar, iiPDF_s_mellor) )

       else ! .not. l_use_modified_corr

          corr_array_1(jvar, ivar) &
          = component_corr_xhm_ip( rc1, cloud_frac1, &
                                   corr_array_cloud(jvar, ivar), &
                                   corr_array_below(jvar, ivar) )

          corr_array_2(jvar, ivar) &
          = component_corr_xhm_ip( rc2, cloud_frac2, &
                                   corr_array_cloud(jvar, ivar), &
                                   corr_array_below(jvar, ivar) )

       endif ! l_use_modified_corr

    enddo


    ! Correlation between w and Ncn
    corr_array_1(iiPDF_Ncn, iiPDF_w) &
    = component_corr_whm_ip( corr_wNcn, rc1, one, &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_w), &
                             corr_array_below(iiPDF_Ncn, iiPDF_w) )

    corr_array_2(iiPDF_Ncn, iiPDF_w) &
    = component_corr_whm_ip( corr_wNcn, rc2, one, &
                             corr_array_cloud(iiPDF_Ncn, iiPDF_w), &
                             corr_array_below(iiPDF_Ncn, iiPDF_w) )

    ! Correlation between w and the hydrometeors
    ivar = iiPDF_w
    do jvar = iiPDF_Ncn+1, d_variables

       corr_array_1(jvar, ivar) &
       = component_corr_whm_ip( corr_whm_1(jvar), rc1, cloud_frac1, &
                                corr_array_cloud(jvar, ivar), corr_array_below(jvar, ivar) )

       corr_array_2(jvar, ivar) &
       = component_corr_whm_ip( corr_whm_2(jvar), rc2, cloud_frac2, &
                                corr_array_cloud(jvar, ivar), corr_array_below(jvar, ivar) )

    enddo

    ! Correlation between Ncn and the hydrometeors
    ivar = iiPDF_Ncn
    do jvar = iiPDF_Ncn+1, d_variables
       corr_array_1(jvar, ivar) &
       = component_corr_hmxhmy_ip( rc1, cloud_frac1, &
                                   corr_array_cloud(jvar, ivar), &
                                   corr_array_below(jvar, ivar) )

       corr_array_2(jvar, ivar) &
       = component_corr_hmxhmy_ip( rc2, cloud_frac2, &
                                   corr_array_cloud(jvar, ivar), &
                                   corr_array_below(jvar, ivar) )
    enddo

    ! Correlation between two hydrometeors
    do ivar = iiPDF_Ncn+1, d_variables-1
       do jvar = ivar+1, d_variables

          corr_array_1(jvar, ivar) &
          = component_corr_hmxhmy_ip( rc1, cloud_frac1, &
                                      corr_array_cloud(jvar, ivar), &
                                      corr_array_below(jvar, ivar) )

          corr_array_2(jvar, ivar) &
          = component_corr_hmxhmy_ip( rc2, cloud_frac2, &
                                      corr_array_cloud(jvar, ivar), &
                                      corr_array_below(jvar, ivar) )

       enddo ! jvar
    enddo ! ivar

    return

  end subroutine compute_corr

  !=============================================================================
  function component_mean_hm_ip( hmi, precip_frac_i, hydromet_tol )  &
  result( mu_hm_i )

    ! Description:
    ! Calculates the in-precip mean of a hydrometeor species within the ith
    ! PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmi,           & ! Mean of hydrometeor, hm (ith PDF component) [hm units]
      precip_frac_i, & ! Precipitation fraction (ith PDF component)  [-]
      hydromet_tol     ! Tolerance value for the hydrometeor         [hm units]

    ! Return Variable
    real( kind = core_rknd ) :: &
      mu_hm_i    ! Mean of hm (ith PDF component) in-precip (ip)     [hm units]


    ! Mean of the hydrometeor (in-precip) in the ith PDF component.
    if ( hmi > hydromet_tol ) then
       mu_hm_i = hmi / precip_frac_i
    else
       ! The mean of the hydrometeor in the ith PDF component is less than the
       ! tolerance amount for the particular hydrometeor.  It is considered to
       ! have a value of 0.  There is not any of this hydrometeor species in the
       ! ith PDF component at this grid level.
       mu_hm_i = zero
    endif


    return

  end function component_mean_hm_ip

  !=============================================================================
  function component_stdev_hm_ip( mu_hm_i, rci, cloud_fraci, &
                                  hm_sigma2_on_mu2_cloud, &
                                  hm_sigma2_on_mu2_below )  &
  result( sigma_hm_i )

    ! Description:
    ! Calculates the in-precip standard deviation of a hydrometeor species
    ! within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_hm_i,     & ! Mean of hm (ith PDF component) in-precip (ip) [hm units]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci    ! Cloud fraction (ith PDF component)            [-]

    real( kind = core_rknd ), intent(in) :: &
      hm_sigma2_on_mu2_cloud, & ! Ratio sigma_hm_1^2/mu_hm_1^2; cloudy levs. [-]
      hm_sigma2_on_mu2_below    ! Ratio sigma_hm_2^2/mu_hm_2^2; clear levs.  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      sigma_hm_i    ! Standard deviation of hm (ith PDF component) ip [hm units]


    ! Standard deviation of the hydrometeor (in-precip) in the
    ! ith PDF component.
    if ( l_interp_prescribed_params ) then
       sigma_hm_i = sqrt( cloud_fraci * hm_sigma2_on_mu2_cloud &
                        + ( one - cloud_fraci ) * hm_sigma2_on_mu2_below ) &
                    * mu_hm_i
    else
       if ( rci > rc_tol ) then
          sigma_hm_i = sqrt( hm_sigma2_on_mu2_cloud ) * mu_hm_i
       else
          sigma_hm_i = sqrt( hm_sigma2_on_mu2_below ) * mu_hm_i
       endif
    endif

    return

  end function component_stdev_hm_ip

  !=============================================================================
  function component_corr_wx( corr_wx, rci, cloud_fraci, &
                              corr_wx_NN_cloud, corr_wx_NN_below ) &
  result( corr_wx_i )

    ! Description:
    ! Calculates the correlation between w and x within the ith PDF component.
    ! Here, x is a variable with a normally distributed individual marginal PDF,
    ! such as s or t.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        zero,   &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_calc_w_corr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_wx,     & ! Correlation between w and x (overall)          [-]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.)  [kg/kg]
      cloud_fraci    ! Cloud fraction (ith PDF component)             [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_wx_NN_cloud, & ! Corr. btwn. w and x (ith PDF comp.); cloudy levs [-]
      corr_wx_NN_below    ! Corr. btwn. w and x (ith PDF comp.); clear levs  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_wx_i    ! Correlation between w and x (ith PDF component)  [-]

    ! Local Variables

    ! The component correlations of w and r_t and the component correlations of
    ! w and theta_l are both set to be 0 within the CLUBB model code.  In other
    ! words, w and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'),
    ! but the single component covariance and correlation are defined to be 0.
    ! Since the component covariances (or correlations) between w and s and
    ! between w and t are based on the covariances (or correlations) between w
    ! and r_t and between w and theta_l, the single component correlation and
    ! covariance of w and s, as well as w and t, are defined to be 0.
    logical, parameter :: &
      l_follow_CLUBB_PDF_standards = .true.


    ! Correlation between w and x in the ith PDF component.
    if ( l_follow_CLUBB_PDF_standards ) then

       ! The component correlations of w and r_t and the component correlations
       ! of w and theta_l are both set to be 0 within the CLUBB model code.  In
       ! other words, w and r_t (theta_l) have overall covariance w'r_t'
       ! (w'theta_l'), but the single component covariance and correlation are
       ! defined to be 0.  Since the component covariances (or correlations)
       ! between w and s and between w and t are based on the covariances (or
       ! correlations) between w and r_t and between w and theta_l, the single
       ! component correlation and covariance of w and s, as well as w and t,
       ! are defined to be 0.
       corr_wx_i = zero

    else ! not following CLUBB PDF standards

       ! WARNING:  the standards used in the generation of the two-component
       !           CLUBB PDF are not being obeyed.  The use of this code is
       !           inconsistent with the rest of CLUBB's PDF.
       if ( l_calc_w_corr ) then
          corr_wx_i = corr_wx
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_wx_i = cloud_fraci * corr_wx_NN_cloud &
                         + ( one - cloud_fraci ) * corr_wx_NN_below
          else
             if ( rci > rc_tol ) then
                corr_wx_i = corr_wx_NN_cloud
             else
                corr_wx_i = corr_wx_NN_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr

    endif ! l_follow_CLUBB_PDF_standards


    return

  end function component_corr_wx

  !=============================================================================
  function component_corr_st( pdf_corr_st_i, rci, cloud_fraci, &
                              corr_st_NN_cloud, corr_st_NN_below, &
                              l_limit_corr_st ) &
  result( corr_st_i )

    ! Description:
    ! Calculates the correlation between s and t within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol, &
        max_mag_correlation

    use parameters_microphys, only: &
        l_fix_s_t_correlations  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      pdf_corr_st_i, & ! Correlation between s and t (ith PDF component)     [-]
      rci,           & ! Mean cloud water mixing ratio (ith PDF comp.)   [kg/kg]
      cloud_fraci      ! Cloud fraction (ith PDF component)                  [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_st_NN_cloud, & ! Corr. btwn. s and t (ith PDF comp.); cloudy levs [-]
      corr_st_NN_below    ! Corr. btwn. s and t (ith PDF comp.); clear levs  [-]

    logical, intent(in) :: &
      l_limit_corr_st    ! We must limit the correlation between s and t if we
                          ! are to take the Cholesky decomposition of the
                          ! resulting correlation matrix. This is because a
                          ! perfect correlation between s and t was found to be
                          ! unrealizable.

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_st_i    ! Correlation between s and t (ith PDF component)         [-]


    ! Correlation between s and t in the ith PDF component.

    ! The PDF variables s and t result from a transformation of the PDF
    ! involving r_t and theta_l.  The correlation between s and t depends on the
    ! correlation between r_t and theta_l, as well as the variances of r_t and
    ! theta_l, and other factors.  The correlation between s and t is subject to
    ! change at every vertical level and model time step, and is calculated as
    ! part of the CLUBB PDF parameters.
    if ( .not. l_fix_s_t_correlations ) then

       ! Preferred, more accurate version.
       corr_st_i = pdf_corr_st_i

    else ! fix the correlation between s and t.

       ! WARNING:  this code is inconsistent with the rest of CLUBB's PDF.  This
       !           code is necessary because SILHS is lazy and wussy, and only
       !           wants to declare correlation arrays at the start of the model
       !           run, rather than updating them throughout the model run.
       if ( l_interp_prescribed_params ) then
          corr_st_i = cloud_fraci * corr_st_NN_cloud &
                      + ( one - cloud_fraci ) * corr_st_NN_below
       else
          if ( rci > rc_tol ) then
             corr_st_i = corr_st_NN_cloud
          else
             corr_st_i = corr_st_NN_below
          endif
       endif

    endif

    ! We cannot have a perfect correlation between s and t if we plan to
    ! decompose this matrix and we don't want the Cholesky_factor code to
    ! throw a fit.
    if ( l_limit_corr_st ) then

       corr_st_i = max( min( corr_st_i, max_mag_correlation ), -max_mag_correlation )

    end if


    return

  end function component_corr_st

  !=============================================================================
  function component_corr_whm_ip( corr_whm_i_in, rci, cloud_fraci, &
                                  corr_whm_NL_cloud, corr_whm_NL_below ) &
  result( corr_whm_i )

    ! Description:
    ! Calculates the in-precip correlation between w and a hydrometeor species
    ! within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_calc_w_corr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_whm_i_in, & ! Correlation between w and hm (ith PDF comp.) ip     [-]
      rci,           & ! Mean cloud water mixing ratio (ith PDF comp.)   [kg/kg]
      cloud_fraci      ! Cloud fraction (ith PDF component)                  [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_whm_NL_cloud, & ! Corr. btwn. w and hm (ith PDF comp.) ip; cloudy [-]
      corr_whm_NL_below    ! Corr. btwn. w and hm (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_whm_i    ! Correlation between w and hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) between w and the hydrometeor in the ith
    ! PDF component.
    if ( l_calc_w_corr ) then
       corr_whm_i = corr_whm_i_in
    else ! use prescribed parameter values
       if ( l_interp_prescribed_params ) then
          corr_whm_i = cloud_fraci * corr_whm_NL_cloud &
                       + ( one - cloud_fraci ) * corr_whm_NL_below
       else
          if ( rci > rc_tol ) then
             corr_whm_i = corr_whm_NL_cloud
          else
             corr_whm_i = corr_whm_NL_below
          endif
       endif ! l_interp_prescribed_params
    endif ! l_calc_w_corr

    return

  end function component_corr_whm_ip

  !=============================================================================
  function component_corr_xhm_ip( rci, cloud_fraci, &
                                  corr_xhm_NL_cloud, corr_xhm_NL_below ) &
  result( corr_xhm_i )

    ! Description:
    ! Calculates the in-precip correlation between x and a hydrometeor species
    ! within the ith PDF component.  Here, x is a variable with a normally
    ! distributed individual marginal PDF, such as s or t.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci    ! Cloud fraction (ith PDF component)            [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_xhm_NL_cloud, & ! Corr. btwn. x and hm (ith PDF comp.) ip; cloudy [-]
      corr_xhm_NL_below    ! Corr. btwn. x and hm (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_xhm_i    ! Correlation between x and hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) between x and the hydrometeor in the ith
    ! PDF component.
    if ( l_interp_prescribed_params ) then
       corr_xhm_i = cloud_fraci * corr_xhm_NL_cloud &
                    + ( one - cloud_fraci ) * corr_xhm_NL_below
    else
       if ( rci > rc_tol ) then
          corr_xhm_i = corr_xhm_NL_cloud
       else
          corr_xhm_i = corr_xhm_NL_below
       endif
    endif

    return

  end function component_corr_xhm_ip

  !=============================================================================
  function component_corr_hmxhmy_ip( rci, cloud_fraci, &
                                     corr_hmxhmy_LL_cloud, &
                                     corr_hmxhmy_LL_below ) &
  result( corr_hmxhmy_i )

    ! Description:
    ! Calculates the in-precip correlation between hydrometeor x and
    ! hydrometeor y within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci    ! Cloud fraction (ith PDF component)            [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_hmxhmy_LL_cloud, & ! Corr.: hmx & hmy (ith PDF comp.) ip; cloudy [-]
      corr_hmxhmy_LL_below    ! Corr.: hmx & hmy (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_hmxhmy_i   ! Correlation between hmx & hmy (ith PDF component) ip [-]


    ! Correlation (in-precip) between hydrometeor x and hydrometeor y in the
    ! ith PDF component.
    if ( l_interp_prescribed_params ) then
       corr_hmxhmy_i = cloud_fraci * corr_hmxhmy_LL_cloud &
                       + ( one - cloud_fraci ) * corr_hmxhmy_LL_below
    else
       if ( rci > rc_tol ) then
          corr_hmxhmy_i = corr_hmxhmy_LL_cloud
       else
          corr_hmxhmy_i = corr_hmxhmy_LL_below
       endif
    endif

    return

  end function component_corr_hmxhmy_ip

  !=============================================================================
  pure function component_corr_thm_ip( corr_st_i, corr_shm_i ) result( corr_thm_i )

    ! Description:
    !   Estimates the correlation between t_mellor and a hydrometeor species
    !   using the correlation between s_mellor and t_mellor, and between
    !   s_mellor and the hydrometeor. This facilities the Cholesky
    !   decomposability of the correlation array that will inevitably be
    !   decomposed for SILHS purposes. Without this estimation, we have found
    !   that the resulting correlation matrix cannot be decomposed.

    ! References:
    !   See clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd       ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_st_i,    & ! Component correlation of s and t                [-]
      corr_shm_i      ! Component correlation of s and the hydrometeor  [-]

    ! Output Variables
    real( kind = core_rknd ) :: &
      corr_thm_i      ! Component correlation of t and the hydrometeor  [-]

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    corr_thm_i = corr_st_i * corr_shm_i

    return
  end function component_corr_thm_ip

  !=============================================================================
  subroutine normalize_mean_stdev( hm1, hm2, Ncnm, d_variables, &
                                   mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                                   sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                                   mu_x_1_n, mu_x_2_n, &
                                   sigma_x_1_n, sigma_x_2_n )

    ! Description:
    ! Calculates the normalized means and the normalized standard deviations
    ! of PDF variables that have assumed lognormal distributions -- which are
    ! precipitating hydrometeors (in precipitation) and N_cn.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        Ncn_tol  ! Constant(s)

    use PDF_utilities, only: &
        mean_L2N,  & ! Procedure(s)
        stdev_L2N

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use corr_matrix_module, only: &
        iiPDF_Ncn  ! Variable(s)

    use parameters_microphys, only: &
        Ncnp2_on_Ncnm2, & ! Variable(s)
        hydromet_tol

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)    [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)    [units vary]

    real( kind = core_rknd ), intent(in) :: &
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >               [num/kg]

    integer, intent(in) :: &
      d_variables    ! Number of variables in CLUBB's PDF

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    ! Local Variable
    integer :: ivar  ! Loop index


    ! The means and standard deviations in each PDF component of w, s, and t do
    ! not need to be normalized, since w, s, and t already follow assumed normal
    ! distributions in each PDF component.  The normalized means and standard
    ! deviations are the same as the actual means and standard deviations.    
    mu_x_1_n = mu_x_1
    mu_x_2_n = mu_x_2
    sigma_x_1_n = sigma_x_1
    sigma_x_2_n = sigma_x_2

    !!! Calculate the normalized mean and standard deviation in each PDF
    !!! component for variables that have an assumed lognormal distribution,
    !!! given the mean and standard deviation in each PDF component for those
    !!! variables.  A precipitating hydrometeor has an assumed lognormal
    !!! distribution in precipitation in each PDF component.  Simplified cloud
    !!! nuclei concentration, N_cn, has an assumed lognormal distribution in
    !!! each PDF component, and furthermore, mu_Ncn_1 = mu_Ncn_2 and
    !!! sigma_Ncn_1 = sigma_Ncn_2, so N_cn has an assumed single lognormal
    !!! distribution over the entire domain.

    ! Normalized mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 1.
    if ( Ncnm > Ncn_tol ) then

       mu_x_1_n(iiPDF_Ncn) = mean_L2N( mu_x_1(iiPDF_Ncn), Ncnp2_on_Ncnm2 )

    else

       ! Mean simplified cloud nuclei concentration in PDF component 1 is less
       ! than the tolerance amount.  It is considered to have a value of 0.
       ! There are not any cloud nuclei or cloud at this grid level.  The value
       ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
       ! assigning it a value.
       mu_x_1_n(iiPDF_Ncn) = -huge( mu_x_1(iiPDF_Ncn) )

    endif

    ! Normalized standard deviation of simplified cloud nuclei concentration,
    ! N_cn, in PDF component 1.
    sigma_x_1_n(iiPDF_Ncn) = stdev_L2N( Ncnp2_on_Ncnm2 )

    ! Normalized mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 2.
    if ( Ncnm > Ncn_tol ) then

       mu_x_2_n(iiPDF_Ncn) = mean_L2N( mu_x_2(iiPDF_Ncn), Ncnp2_on_Ncnm2 )

    else

       ! Mean simplified cloud nuclei concentration in PDF component 1 is less
       ! than the tolerance amount.  It is considered to have a value of 0.
       ! There are not any cloud nuclei or cloud at this grid level.  The value
       ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
       ! assigning it a value.
       mu_x_2_n(iiPDF_Ncn) = -huge( mu_x_2(iiPDF_Ncn) )

    endif

    ! Normalized standard deviation of simplified cloud nuclei concentration,
    ! N_cn, in PDF component 2.
    sigma_x_2_n(iiPDF_Ncn) = stdev_L2N( Ncnp2_on_Ncnm2 )


    ! Normalize precipitating hydrometeor means and standard deviations.
    do ivar = iiPDF_Ncn+1, d_variables, 1

       ! Normalized mean of a precipitating hydrometeor, hm, in PDF component 1.
       if ( hm1(pdf2hydromet_idx(ivar)) &
            > hydromet_tol(pdf2hydromet_idx(ivar)) ) then

          mu_x_1_n(ivar) = mean_L2N( mu_x_1(ivar), sigma2_on_mu2_ip_1(ivar) )

       else

          ! The mean of a precipitating hydrometeor in PDF component 1 is less
          ! than its tolerance amount.  It is considered to have a value of 0.
          ! There is not any of this precipitating hydrometeor in the 1st PDF
          ! component at this grid level.  The in-precip mean of this
          ! precipitating hydrometeor (1st PDF component) is also 0.  The value
          ! of mu_hm_1_n should be -inf.  It will be set to -huge for purposes
          ! of assigning it a value.
          mu_x_1_n(ivar) = -huge( mu_x_1(ivar) )

       endif

       ! Normalized standard deviation of a precipitating hydrometeor, hm, in
       ! PDF component 1.
       sigma_x_1_n(ivar) = stdev_L2N( sigma2_on_mu2_ip_1(ivar) )

       ! Normalized mean of a precipitating hydrometeor, hm, in PDF component 2.
       if ( hm2(pdf2hydromet_idx(ivar)) &
            > hydromet_tol(pdf2hydromet_idx(ivar)) ) then

          mu_x_2_n(ivar) = mean_L2N( mu_x_2(ivar), sigma2_on_mu2_ip_2(ivar) )

       else

          ! The mean of a precipitating hydrometeor in PDF component 2 is less
          ! than its tolerance amount.  It is considered to have a value of 0.
          ! There is not any of this precipitating hydrometeor in the 2nd PDF
          ! component at this grid level.  The in-precip mean of this
          ! precipitating hydrometeor (2nd PDF component) is also 0.  The value
          ! of mu_hm_2_n should be -inf.  It will be set to -huge for purposes
          ! of assigning it a value.
          mu_x_2_n(ivar) = -huge( mu_x_2(ivar) )

       endif

       ! Normalized standard deviation of a precipitating hydrometeor, hm, in
       ! PDF component 2.
       sigma_x_2_n(ivar) = stdev_L2N( sigma2_on_mu2_ip_2(ivar) )

    enddo ! ivar = iiPDF_Ncn+1, d_variables, 1


    return

  end subroutine normalize_mean_stdev

  !=============================================================================
  subroutine normalize_corr( d_variables, sigma_x_1_n, sigma_x_2_n, &
                             sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                             corr_array_1, corr_array_2, &
                             corr_array_1_n, corr_array_2_n )

    ! Description:
    ! Calculates the normalized correlations between PDF variables, where at
    ! least one of the variables that is part of a correlation has an assumed
    ! lognormal distribution -- which are the precipitating hydrometeors (in
    ! precipitation) and N_cn.

    ! References:
    !-----------------------------------------------------------------------

    use PDF_utilities, only: &
        corr_NL2NN, & ! Procedure(s)
        corr_LL2NN

    use corr_matrix_module, only: &
        iiPDF_s_mellor, & ! Variable(s)
        iiPDF_t_mellor, &
        iiPDF_w,        &
        iiPDF_Ncn

    use parameters_microphys, only: &
        Ncnp2_on_Ncnm2  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables ! Number of PDF variables

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    real ( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2_n    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    ! Local Variables
    integer :: ivar, jvar ! Loop indices


    ! The correlations in each PDF component between two of w, s, and t do not
    ! need to be normalized, since w, s, and t already follow assumed normal
    ! distributions in each PDF component.  The normalized correlations between
    ! any two of these variables are the same as the actual correlations.    
    corr_array_1_n = corr_array_1
    corr_array_2_n = corr_array_2


    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    ! Normalize the correlations between s/t/w and N_cn.

    ! Normalize the correlation between w and N_cn in PDF component 1.
    corr_array_1_n(iiPDF_Ncn, iiPDF_w) &
    = corr_NL2NN( corr_array_1(iiPDF_Ncn, iiPDF_w), sigma_x_1_n(iiPDF_Ncn), &
                  Ncnp2_on_Ncnm2 )

    ! Normalize the correlation between w and N_cn in PDF component 2.
    corr_array_2_n(iiPDF_Ncn, iiPDF_w) &
    = corr_NL2NN( corr_array_2(iiPDF_Ncn, iiPDF_w), sigma_x_2_n(iiPDF_Ncn), &
                  Ncnp2_on_Ncnm2 )

    ! Normalize the correlation between s and N_cn in PDF component 1.
    corr_array_1_n(iiPDF_Ncn, iiPDF_s_mellor) &
    = corr_NL2NN( corr_array_1(iiPDF_Ncn, iiPDF_s_mellor), &
                  sigma_x_1_n(iiPDF_Ncn), Ncnp2_on_Ncnm2 )

    ! Normalize the correlation between s and N_cn in PDF component 2.
    corr_array_2_n(iiPDF_Ncn, iiPDF_s_mellor) &
    = corr_NL2NN( corr_array_2(iiPDF_Ncn, iiPDF_s_mellor), &
                  sigma_x_2_n(iiPDF_Ncn), Ncnp2_on_Ncnm2 )

    ! Normalize the correlation between t and N_cn in PDF component 1.
    corr_array_1_n(iiPDF_Ncn, iiPDF_t_mellor) &
    = corr_NL2NN( corr_array_1(iiPDF_Ncn, iiPDF_t_mellor), &
                  sigma_x_1_n(iiPDF_Ncn), Ncnp2_on_Ncnm2 )

    ! Normalize the correlation between t and N_cn in PDF component 2.
    corr_array_2_n(iiPDF_Ncn, iiPDF_t_mellor) &
    = corr_NL2NN( corr_array_2(iiPDF_Ncn, iiPDF_t_mellor), &
                  sigma_x_2_n(iiPDF_Ncn), Ncnp2_on_Ncnm2 )

    ! Normalize the correlations (in-precip) between s/t/w and the precipitating
    ! hydrometeors.
    do ivar = iiPDF_s_mellor, iiPDF_w
       do jvar = iiPDF_Ncn+1, d_variables

          ! Normalize the correlation (in-precip) between w, s, or t and a
          ! precipitating hydrometeor, hm, in PDF component 1.
          corr_array_1_n(jvar, ivar) &
          = corr_NL2NN( corr_array_1(jvar, ivar), sigma_x_1_n(jvar), &
                        sigma2_on_mu2_ip_1(jvar) )

          ! Normalize the correlation (in-precip) between w, s, or t and a
          ! precipitating hydrometeor, hm, in PDF component 2.
          corr_array_2_n(jvar, ivar) &
          = corr_NL2NN( corr_array_2(jvar, ivar), sigma_x_2_n(jvar), &
                        sigma2_on_mu2_ip_2(jvar) )

       enddo ! jvar = iiPDF_Ncn+1, d_variables
    enddo ! ivar = iiPDF_s_mellor, iiPDF_w


    !!! Calculate the normalized correlation between two variables that both
    !!! have an assumed lognormal distribution for the ith PDF component, given
    !!! their correlation and both of their normalized standard deviations.

    ! Normalize the correlations (in-precip) between N_cn and the precipitating
    ! hydrometeors.
    ivar = iiPDF_Ncn
    do jvar = ivar+1, d_variables

       ! Normalize the correlation (in-precip) between N_cn and a precipitating
       ! hydrometeor, hm, in PDF component 1.
       corr_array_1_n(jvar, ivar) &
       = corr_LL2NN( corr_array_1(jvar, ivar), &
                     sigma_x_1_n(ivar), sigma_x_1_n(jvar), &
                     Ncnp2_on_Ncnm2, sigma2_on_mu2_ip_1(jvar) )

       ! Normalize the correlation (in-precip) between N_cn and a precipitating
       ! hydrometeor, hm, in PDF component 2.
       corr_array_2_n(jvar, ivar) &
       = corr_LL2NN( corr_array_2(jvar, ivar), &
                     sigma_x_2_n(ivar), sigma_x_2_n(jvar), &
                     Ncnp2_on_Ncnm2, sigma2_on_mu2_ip_2(jvar) )

    enddo ! jvar = ivar+1, d_variables

    ! Normalize the correlations (in-precip) between two precipitating
    ! hydrometeors.
    do ivar = iiPDF_Ncn+1, d_variables-1
       do jvar = ivar+1, d_variables

          ! Normalize the correlation (in-precip) between two precipitating
          ! hydrometeors (for example, r_r and N_r) in PDF component 1.
          corr_array_1_n(jvar, ivar) &
          = corr_LL2NN( corr_array_1(jvar, ivar), &
                        sigma_x_1_n(ivar), sigma_x_1_n(jvar), &
                        sigma2_on_mu2_ip_1(ivar), sigma2_on_mu2_ip_1(jvar) )

          ! Normalize the correlation (in-precip) between two precipitating
          ! hydrometeors (for example, r_r and N_r) in PDF component 2.
          corr_array_2_n(jvar, ivar) &
          = corr_LL2NN( corr_array_2(jvar, ivar), &
                        sigma_x_2_n(ivar), sigma_x_2_n(jvar), &
                        sigma2_on_mu2_ip_2(ivar), sigma2_on_mu2_ip_2(jvar) )

       enddo ! jvar = ivar+1, d_variables
    enddo ! ivar = iiPDF_Ncn+1, d_variables-1


    return

  end subroutine normalize_corr

  !=============================================================================
  subroutine calc_corr_whm( wm, wphydrometp, &
                            mu_w_1, mu_w_2, &
                            mu_hm_1, mu_hm_2, &
                            sigma_w_1, sigma_w_2, &
                            sigma_hm_1, sigma_hm_2, &
                            mixt_frac, precip_frac_1, precip_frac_2, &
                            corr_whm_1, corr_whm_2 )

    ! Description:
    ! Calculates the PDF component correlation (in-precip) between vertical
    ! velocity, w, and a hydrometeor, hm.  The overall covariance of w and hm,
    ! <w'hm'> can be written in terms of the PDF parameters.  When both w and hm
    ! vary in both PDF components, the equation is written as:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( ( mu_w_1 - <w> ) * mu_hm_1
    !               + corr_wrr_1 * sigma_w_1 * sigma_rr_1 )
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( ( mu_w_2 - <w> ) * mu_hm_2
    !                 + corr_wrr_2 * sigma_w_2 * sigma_rr_2 ).
    !
    ! The overall covariance is provided, so the component correlation is solved
    ! by setting corr_wrr_1 = corr_wrr_2 ( = corr_wrr ).  The equation is:
    !
    ! corr_wrr
    ! = ( <w'hm'>
    !     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
    !     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
    !   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1
    !       + ( 1 - mixt_frac ) * precip_frac_2 * sigma_w_2 * sigma_hm_2 );
    !
    ! again, where corr_wrr_1 = corr_wrr_2 = corr_wrr.  When either w or hm is
    ! constant in one PDF component, but both w and hm vary in the other PDF
    ! component, the equation for <w'hm'> is written as:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( ( mu_w_1 - <w> ) * mu_hm_1
    !               + corr_wrr_1 * sigma_w_1 * sigma_rr_1 )
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w> ) * mu_hm_2.
    !
    ! In the above equation, either w or hm (or both) is (are) constant in PDF
    ! component 2, but both w and hm vary in PDF component 1.  When both w and
    ! hm vary in PDF component 2, but at least one of w or hm is constant in PDF
    ! component 1, the equation is similar.  The above equation can be rewritten
    ! to solve for corr_wrr_1, such that:
    !
    ! corr_wrr_1
    ! = ( <w'hm'>
    !     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
    !     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
    !   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1 ).
    !
    ! Since either w or hm is constant in PDF component 2, corr_wrr_2 is
    ! undefined.  When both w and hm vary in PDF component 2, but at least one
    ! of w or hm is constant in PDF component 1, the equation is similar, but
    ! is in terms of corr_wrr_2, while corr_wrr_1 is undefined.  When either w
    ! or hm is constant in both PDF components, the equation for <w'hm'> is:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( mu_w_1 - <w> ) * mu_hm_1
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w> ) * mu_hm_2.
    !
    ! When this is the case, both corr_wrr_1 and corr_wrr_2 are undefined.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,                 & ! Constant(s)
        zero,                &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      wm,            & ! Mean vertical velocity (overall), <w>             [m/s]
      wphydrometp,   & ! Covariance of w and hm (overall), <w'hm'>  [m/s(hm un)]
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_hm_1,       & ! Mean of hm (1st PDF component) in-precip (ip)   [hm un]
      mu_hm_2,       & ! Mean of hm (2nd PDF component) ip               [hm un]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_hm_1,    & ! Standard deviation of hm (1st PDF component) ip [hm un]
      sigma_hm_2,    & ! Standard deviation of hm (2nd PDF component) ip [hm un]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      corr_whm_1, & ! Correlation between w and hm (1st PDF component) ip    [-]
      corr_whm_2    ! Correlation between w and hm (2nd PDF component) ip    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_whm    ! Correlation between w and hm (both PDF components) ip    [-]


    ! Calculate the PDF component correlation between vertical velocity, w, and
    ! a hydrometeor, hm, in precipitation.
    if ( sigma_w_1 * sigma_hm_1 > zero .and. &
         sigma_w_2 * sigma_hm_2 > zero ) then

       ! Both w and hm vary in both PDF components.
       ! Calculate corr_whm (where corr_whm_1 = corr_whm_2 = corr_whm).
       corr_whm &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1 &
             + ( one - mixt_frac ) * precip_frac_2 * sigma_w_2 * sigma_hm_2 )

       ! Check that the PDF component correlations have reasonable values.
       if ( corr_whm > max_mag_correlation ) then
          corr_whm = max_mag_correlation
       elseif ( corr_whm < -max_mag_correlation ) then
          corr_whm = -max_mag_correlation
       endif

       ! The PDF component correlations between w and hm (in-precip) are equal.
       corr_whm_1 = corr_whm
       corr_whm_2 = corr_whm


    elseif ( sigma_w_1 * sigma_hm_1 > zero ) then

       ! Both w and hm vary in PDF component 1, but at least one of w and hm is
       ! constant in PDF component 2.
       ! Calculate the PDF component 1 correlation between w and hm (in-precip).
       corr_whm_1 &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1 )

       ! Check that the PDF component 1 correlation has a reasonable value.
       if ( corr_whm_1 > max_mag_correlation ) then
          corr_whm_1 = max_mag_correlation
       elseif ( corr_whm_1 < -max_mag_correlation ) then
          corr_whm_1 = -max_mag_correlation
       endif

       ! The PDF component 2 correlation is undefined.
       corr_whm_2 = zero
       

    elseif ( sigma_w_2 * sigma_hm_2 > zero ) then

       ! Both w and hm vary in PDF component 2, but at least one of w and hm is
       ! constant in PDF component 1.
       ! Calculate the PDF component 2 correlation between w and hm (in-precip).
       corr_whm_2 &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( ( one - mixt_frac ) * precip_frac_2 * sigma_w_2 * sigma_hm_2 )

       ! Check that the PDF component 2 correlation has a reasonable value.
       if ( corr_whm_2 > max_mag_correlation ) then
          corr_whm_2 = max_mag_correlation
       elseif ( corr_whm_2 < -max_mag_correlation ) then
          corr_whm_2 = -max_mag_correlation
       endif

       ! The PDF component 1 correlation is undefined.
       corr_whm_1 = zero
       

    else    ! sigma_w_1 * sigma_hm_1 = 0 .and. sigma_w_2 * sigma_hm_2 = 0.

       ! At least one of w and hm is constant in both PDF components.

       ! The PDF component 1 and component 2 correlations are both undefined.
       corr_whm_1 = zero
       corr_whm_2 = zero


    endif


    return

  end subroutine calc_corr_whm

  !=============================================================================
  subroutine pdf_param_hm_stats( mu_rr_1, mu_rr_2, mu_Nr_1, &
                                 mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                 sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                 sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                 corr_wrr_1, corr_wrr_2, corr_wNr_1, &
                                 corr_wNr_2, corr_wNcn_1, corr_wNcn_2, &
                                 corr_srr_1, corr_srr_2, corr_sNr_1, &
                                 corr_sNr_2, corr_sNcn_1, corr_sNcn_2, &
                                 corr_trr_1, corr_trr_2, corr_tNr_1, &
                                 corr_tNr_2, corr_tNcn_1, corr_tNcn_2, &
                                 corr_rrNr_1, corr_rrNr_2, &
                                 level, l_stats_samp )

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
        imu_Ncn_1,      &
        imu_Ncn_2,      &
        isigma_rr_1,    &
        isigma_rr_2,    &
        isigma_Nr_1,    &
        isigma_Nr_2,    &
        isigma_Ncn_1,   &
        isigma_Ncn_2


    use stats_variables, only : &
        icorr_wrr_1,    & ! Variable(s)
        icorr_wrr_2,    &
        icorr_wNr_1,    &
        icorr_wNr_2,    &
        icorr_wNcn_1,   &
        icorr_wNcn_2,   &
        icorr_srr_1,    &
        icorr_srr_2,    &
        icorr_sNr_1,    &
        icorr_sNr_2,    &
        icorr_sNcn_1,   &
        icorr_sNcn_2,   &
        icorr_trr_1,    &
        icorr_trr_2,    &
        icorr_tNr_1,    &
        icorr_tNr_2,    &
        icorr_tNcn_1,   &
        icorr_tNcn_2,   &
        icorr_rrNr_1,   &
        icorr_rrNr_2, &
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
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2      ! Standard deviation of Ncn (2nd PDF component)  [num/kg]


    real( kind = core_rknd ), intent(in) :: &
      corr_wrr_1,    & ! Correlation between w and rr (1st PDF component) ip [-]
      corr_wrr_2,    & ! Correlation between w and rr (2nd PDF component) ip [-]
      corr_wNr_1,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      corr_wNr_2,    & ! Correlation between w and Nr (2nd PDF component) ip [-]
      corr_wNcn_1,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      corr_wNcn_2,   & ! Correlation between w and Ncn (2nd PDF component)   [-]
      corr_srr_1,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_srr_2,    & ! Correlation between s and rr (2nd PDF component) ip [-]
      corr_sNr_1,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNr_2,    & ! Correlation between s and Nr (2nd PDF component) ip [-]
      corr_sNcn_1,   & ! Correlation between s and Ncn (1st PDF component)   [-]
      corr_sNcn_2,   & ! Correlation between s and Ncn (2nd PDF component)   [-]
      corr_trr_1,    & ! Correlation between t and rr (1st PDF component) ip [-]
      corr_trr_2,    & ! Correlation between t and rr (2nd PDF component) ip [-]
      corr_tNr_1,    & ! Correlation between t and Nr (1st PDF component) ip [-]
      corr_tNr_2,    & ! Correlation between t and Nr (2nd PDF component) ip [-]
      corr_tNcn_1,   & ! Correlation between t and Ncn (1st PDF component)   [-]
      corr_tNcn_2,   & ! Correlation between t and Ncn (2nd PDF component)   [-]
      corr_rrNr_1,   & ! Correlation between rr & Nr (1st PDF component) ip  [-]
      corr_rrNr_2      ! Correlation between rr & Nr (2nd PDF component) ip  [-]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.


    !!! Output the statistics for upscaled KK.

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

       ! Mean of cloud nuclei concentration in PDF component 1.
       if ( imu_Ncn_1 > 0 ) then
          call stat_update_var_pt( imu_Ncn_1, level, mu_Ncn_1, zt )
       endif

       ! Mean of cloud nuclei concentration in PDF component 2.
       if ( imu_Ncn_2 > 0 ) then
          call stat_update_var_pt( imu_Ncn_2, level, mu_Ncn_2, zt )
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

       ! Standard deviation of cloud nuclei concentration in PDF component 1.
       if ( isigma_Ncn_1 > 0 ) then
          call stat_update_var_pt( isigma_Ncn_1, level, sigma_Ncn_1, zt )
       endif

       ! Standard deviation of cloud nuclei concentration in PDF component 2.
       if ( isigma_Ncn_2 > 0 ) then
          call stat_update_var_pt( isigma_Ncn_2, level, sigma_Ncn_2, zt )
       endif

       ! Correlation (in-precip) between w and r_r in PDF component 1.
       if ( icorr_wrr_1 > 0 ) then
          call stat_update_var_pt( icorr_wrr_1, level, corr_wrr_1, zt )
       endif

       ! Correlation (in-precip) between w and r_r in PDF component 2.
       if ( icorr_wrr_2 > 0 ) then
          call stat_update_var_pt( icorr_wrr_2, level, corr_wrr_2, zt )
       endif

       ! Correlation (in-precip) between w and N_r in PDF component 1.
       if ( icorr_wNr_1 > 0 ) then
          call stat_update_var_pt( icorr_wNr_1, level, corr_wNr_1, zt )
       endif

       ! Correlation (in-precip) between w and N_r in PDF component 2.
       if ( icorr_wNr_2 > 0 ) then
          call stat_update_var_pt( icorr_wNr_2, level, corr_wNr_2, zt )
       endif

       ! Correlation between w and N_cn in PDF component 1.
       if ( icorr_wNcn_1 > 0 ) then
          call stat_update_var_pt( icorr_wNcn_1, level, corr_wNcn_1, zt )
       endif

       ! Correlation between w and N_cn in PDF component 2.
       if ( icorr_wNcn_2 > 0 ) then
          call stat_update_var_pt( icorr_wNcn_2, level, corr_wNcn_2, zt )
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

       ! Correlation between s and N_cn in PDF component 1.
       if ( icorr_sNcn_1 > 0 ) then
          call stat_update_var_pt( icorr_sNcn_1, level, corr_sNcn_1, zt )
       endif

       ! Correlation between s and N_cn in PDF component 2.
       if ( icorr_sNcn_2 > 0 ) then
          call stat_update_var_pt( icorr_sNcn_2, level, corr_sNcn_2, zt )
       endif

       ! Correlation (in-precip) between t and r_r in PDF component 1.
       if ( icorr_trr_1 > 0 ) then
          call stat_update_var_pt( icorr_trr_1, level, corr_trr_1, zt )
       endif

       ! Correlation (in-precip) between t and r_r in PDF component 2.
       if ( icorr_trr_2 > 0 ) then
          call stat_update_var_pt( icorr_trr_2, level, corr_trr_2, zt )
       endif

       ! Correlation (in-precip) between t and N_r in PDF component 1.
       if ( icorr_tNr_1 > 0 ) then
          call stat_update_var_pt( icorr_tNr_1, level, corr_tNr_1, zt )
       endif

       ! Correlation (in-precip) between t and N_r in PDF component 2.
       if ( icorr_tNr_2 > 0 ) then
          call stat_update_var_pt( icorr_tNr_2, level, corr_tNr_2, zt )
       endif

       ! Correlation between t and N_cn in PDF component 1.
       if ( icorr_tNcn_1 > 0 ) then
          call stat_update_var_pt( icorr_tNcn_1, level, corr_tNcn_1, zt )
       endif

       ! Correlation between t and N_cn in PDF component 2.
       if ( icorr_tNcn_2 > 0 ) then
          call stat_update_var_pt( icorr_tNcn_2, level, corr_tNcn_2, zt )
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 1.
       if ( icorr_rrNr_1 > 0 ) then
          call stat_update_var_pt( icorr_rrNr_1, level, corr_rrNr_1, zt )
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 2.
       if ( icorr_rrNr_2 > 0 ) then
          call stat_update_var_pt( icorr_rrNr_2, level, corr_rrNr_2, zt )
       endif

    endif ! l_stats_samp

    return


  end subroutine pdf_param_hm_stats

  !=============================================================================
  subroutine pdf_param_log_hm_stats( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                     mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                     sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                     sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                     corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &
                                     corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &
                                     corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                                     corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                                     corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &
                                     corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &
                                     corr_rrNr_1_n, corr_rrNr_2_n, level, &
                                     l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        imu_rr_1_n,     & ! Variable(s)
        imu_rr_2_n,     &
        imu_Nr_1_n,     &
        imu_Nr_2_n,     &
        imu_Ncn_1_n,    &
        imu_Ncn_2_n,    &
        isigma_rr_1_n,  &
        isigma_rr_2_n,  &
        isigma_Nr_1_n,  &
        isigma_Nr_2_n,  &
        isigma_Ncn_1_n, &
        isigma_Ncn_2_n

    use stats_variables, only : &
        icorr_wrr_1_n,  & ! Variables
        icorr_wrr_2_n,  &
        icorr_wNr_1_n,  &
        icorr_wNr_2_n,  &
        icorr_wNcn_1_n, &
        icorr_wNcn_2_n, &
        icorr_srr_1_n,  &
        icorr_srr_2_n,  &
        icorr_sNr_1_n,  &
        icorr_sNr_2_n,  &
        icorr_sNcn_1_n, &
        icorr_sNcn_2_n, &
        icorr_trr_1_n,  &
        icorr_trr_2_n,  &
        icorr_tNr_1_n,  &
        icorr_tNr_2_n,  &
        icorr_tNcn_1_n, &
        icorr_tNcn_2_n, &
        icorr_rrNr_1_n, &
        icorr_rrNr_2_n, &
        zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n    ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]

    real( kind = core_rknd ), intent(in) :: &
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

    integer, intent(in) :: &
      level   ! Vertical level index

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.


    !!! Output the statistics for upscaled KK.

    ! Statistics
    if ( l_stats_samp ) then

       ! Mean (in-precip) of ln r_r in PDF component 1.
       if ( imu_rr_1_n > 0 ) then
          if ( mu_rr_1_n > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_rr_1_n, level, mu_rr_1_n, zt )
          else
             ! When rr1 is 0 (or below tolerance value), mu_rr_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_rr_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln r_r in PDF component 2.
       if ( imu_rr_2_n > 0 ) then
          if ( mu_rr_2_n > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_rr_2_n, level, mu_rr_2_n, zt )
          else
             ! When rr2 is 0 (or below tolerance value), mu_rr_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_rr_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln N_r in PDF component 1.
       if ( imu_Nr_1_n > 0 ) then
          if ( mu_Nr_1_n > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_Nr_1_n, level, mu_Nr_1_n, zt )
          else
             ! When Nr1 is 0 (or below tolerance value), mu_Nr_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nr_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln N_r in PDF component 2.
       if ( imu_Nr_2_n > 0 ) then
          if ( mu_Nr_2_n > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_Nr_2_n, level, mu_Nr_2_n, zt )
          else
             ! When Nr2 is 0 (or below tolerance value), mu_Nr_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nr_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean of ln N_cn in PDF component 1.
       if ( imu_Ncn_1_n > 0 ) then
          if ( mu_Ncn_1_n > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_Ncn_1_n, level, mu_Ncn_1_n, zt )
          else
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_1_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Ncn_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean of ln N_cn in PDF component 2.
       if ( imu_Ncn_2_n > 0 ) then
          if ( mu_Ncn_2_n > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_Ncn_2_n, level, mu_Ncn_2_n, zt )
          else
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_2_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Ncn_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
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

       ! Standard deviation of ln N_cn in PDF component 1.
       if ( isigma_Ncn_1_n > 0 ) then
          call stat_update_var_pt( isigma_Ncn_1_n, level, sigma_Ncn_1_n, zt )
       endif

       ! Standard deviation of ln N_cn in PDF component 2.
       if ( isigma_Ncn_2_n > 0 ) then
          call stat_update_var_pt( isigma_Ncn_2_n, level, sigma_Ncn_2_n, zt )
       endif

       ! Correlation (in-precip) between w and ln r_r in PDF component 1.
       if ( icorr_wrr_1_n > 0 ) then
          call stat_update_var_pt( icorr_wrr_1_n, level, corr_wrr_1_n, zt )
       endif

       ! Correlation (in-precip) between w and ln r_r in PDF component 2.
       if ( icorr_wrr_2_n > 0 ) then
          call stat_update_var_pt( icorr_wrr_2_n, level, corr_wrr_2_n, zt )
       endif

       ! Correlation (in-precip) between w and ln N_r in PDF component 1.
       if ( icorr_wNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_wNr_1_n, level, corr_wNr_1_n, zt )
       endif

       ! Correlation (in-precip) between w and ln N_r in PDF component 2.
       if ( icorr_wNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_wNr_2_n, level, corr_wNr_2_n, zt )
       endif

       ! Correlation between w and ln N_cn in PDF component 1.
       if ( icorr_wNcn_1_n > 0 ) then
          call stat_update_var_pt( icorr_wNcn_1_n, level, corr_wNcn_1_n, zt )
       endif

       ! Correlation between w and ln N_cn in PDF component 2.
       if ( icorr_wNcn_2_n > 0 ) then
          call stat_update_var_pt( icorr_wNcn_2_n, level, corr_wNcn_2_n, zt )
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

       ! Correlation between s and ln N_cn in PDF component 1.
       if ( icorr_sNcn_1_n > 0 ) then
          call stat_update_var_pt( icorr_sNcn_1_n, level, corr_sNcn_1_n, zt )
       endif

       ! Correlation between s and ln N_cn in PDF component 2.
       if ( icorr_sNcn_2_n > 0 ) then
          call stat_update_var_pt( icorr_sNcn_2_n, level, corr_sNcn_2_n, zt )
       endif

       ! Correlation (in-precip) between t and ln r_r in PDF component 1.
       if ( icorr_trr_1_n > 0 ) then
          call stat_update_var_pt( icorr_trr_1_n, level, corr_trr_1_n, zt )
       endif

       ! Correlation (in-precip) between t and ln r_r in PDF component 2.
       if ( icorr_trr_2_n > 0 ) then
          call stat_update_var_pt( icorr_trr_2_n, level, corr_trr_2_n, zt )
       endif

       ! Correlation (in-precip) between t and ln N_r in PDF component 1.
       if ( icorr_tNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_tNr_1_n, level, corr_tNr_1_n, zt )
       endif

       ! Correlation (in-precip) between t and ln N_r in PDF component 2.
       if ( icorr_tNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_tNr_2_n, level, corr_tNr_2_n, zt )
       endif

       ! Correlation between t and ln N_cn in PDF component 1.
       if ( icorr_tNcn_1_n > 0 ) then
          call stat_update_var_pt( icorr_tNcn_1_n, level, corr_tNcn_1_n, zt )
       endif

       ! Correlation between t and ln N_cn in PDF component 2.
       if ( icorr_tNcn_2_n > 0 ) then
          call stat_update_var_pt( icorr_tNcn_2_n, level, corr_tNcn_2_n, zt )
       endif

       ! Correlation (in-precip) between ln r_r and ln N_r in PDF component 1.
       if ( icorr_rrNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_rrNr_1_n, level, corr_rrNr_1_n, zt )
       endif

       ! Correlation (in-precip) between ln r_r and ln N_r in PDF component 2.
       if ( icorr_rrNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_rrNr_2_n, level, corr_rrNr_2_n, zt )
       endif

    endif ! l_stats_samp

    return


  end subroutine pdf_param_log_hm_stats

  !=============================================================================
  subroutine pack_pdf_params( hm1, hm2, mu_x_1, mu_x_2, &                   ! Intent(in)
                              sigma_x_1, sigma_x_2, d_variables, &          ! Intent(in)
                              precip_frac, precip_frac_1, precip_frac_2, &  ! Intent(in)
                              hydromet_pdf_params )                         ! Intent(out)

    ! Description: Setup the structure hydromet_pdf_params.

    !!! This code still has to be generalized to deal with an arbitrary set of
    !!! hydrometeors.

    ! References:
    !-----------------------------------------------------------------------

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Constant

    use corr_matrix_module, only: &
        iiPDF_rrain, &
        iiPDF_Nr, &
        iiPDF_Ncn

    use array_index, only: &
        iirrainm, & ! Variable(s)
        iiNrm

    implicit none

    ! Input Variables
    integer, intent(in) :: d_variables ! Numeber of variables in the mean/stdev arrays

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,       & ! Mean array (1st PDF component) in-precip (ip)   [units vary]
      mu_x_2,       & ! Mean array (2nd PDF component) ip               [units vary]
      sigma_x_1,    & ! Standard deviation array (1st PDF component) ip [units vary]
      sigma_x_2       ! Standard deviation array (2nd PDF component) ip [units vary]

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)  [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)  [units vary]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    ! Output Variables
    type(hydromet_pdf_parameter), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! ---- Begin Code ----

    ! Pack remaining variables into hydromet_pdf_params
    if ( iirrainm > 0 ) then
       hydromet_pdf_params%rr1           = hm1(iirrainm)
       hydromet_pdf_params%rr2           = hm2(iirrainm)

       hydromet_pdf_params%mu_rr_1     = mu_x_1(iiPDF_rrain)
       hydromet_pdf_params%mu_rr_2     = mu_x_2(iiPDF_rrain)
       hydromet_pdf_params%sigma_rr_1  = sigma_x_1(iiPDF_rrain)
       hydromet_pdf_params%sigma_rr_2  = sigma_x_2(iiPDF_rrain)
    endif

    if ( iiNrm > 0 ) then
       hydromet_pdf_params%Nr1           = hm1(iiNrm)
       hydromet_pdf_params%Nr2           = hm2(iiNrm)

       hydromet_pdf_params%mu_Nr_1     = mu_x_1(iiPDF_Nr)
       hydromet_pdf_params%mu_Nr_2     = mu_x_2(iiPDF_Nr)
       hydromet_pdf_params%sigma_Nr_1  = sigma_x_1(iiPDF_Nr)
       hydromet_pdf_params%sigma_Nr_2  = sigma_x_2(iiPDF_Nr)
    endif

    hydromet_pdf_params%mu_Ncn_1    = mu_x_1(iiPDF_Ncn)
    hydromet_pdf_params%mu_Ncn_2    = mu_x_2(iiPDF_Ncn)
    hydromet_pdf_params%sigma_Ncn_1 = sigma_x_1(iiPDF_Ncn)
    hydromet_pdf_params%sigma_Ncn_2 = sigma_x_2(iiPDF_Ncn)

    hydromet_pdf_params%precip_frac   = precip_frac
    hydromet_pdf_params%precip_frac_1 = precip_frac_1
    hydromet_pdf_params%precip_frac_2 = precip_frac_2

    return

  end subroutine pack_pdf_params

  !=============================================================================
  subroutine unpack_pdf_params( d_variables, corr_array_1, corr_array_2, &    ! Intent(in)
                                mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &       ! Intent(in)
                                hydromet_pdf_params, &                        ! Intent(in)
                                mu_w_1, mu_w_2, mu_s_1, mu_s_2, &             ! Intent(out)
                                mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &           ! Intent(out)
                                mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &       ! Intent(out)
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(out)
                                mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(out)
                                sigma_w_1, sigma_w_2, sigma_s_1, &            ! Intent(out)
                                sigma_s_2, sigma_t_1, sigma_t_2, &            ! Intent(out)
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &         ! Intent(out)
                                sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &       ! Intent(out)
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(out)
                                sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(out)
                                corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, & ! Intent(out)
                                corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(out)
                                corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(out)
                                corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(out)
                                corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(out)
                                corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(out)
                                corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(out)
                                corr_rrNr_1_n, corr_rrNr_2_n , &              ! Intent(out)
                                rr1, rr2, Nr1, Nr2, &                         ! Intent(out)
                                precip_frac, precip_frac_1, precip_frac_2 )   ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use corr_matrix_module, only: &
        iiPDF_w,        & ! Variable(s)
        iiPDF_s_mellor, &
        iiPDF_t_mellor, &
        iiPDF_rrain,    &
        iiPDF_Nr,       &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array.

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
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

    real( kind = core_rknd ), intent(out) :: &
      corr_ws_1,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_ws_2,     & ! Correlation between w and s (2nd PDF component)     [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2        ! Correlation between s and t (2nd PDF component)     [-]

   real( kind = core_rknd ), intent(out) :: &
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

    real( kind = core_rknd ), intent(out) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]


    ! Unpack mu_x_i and sigma_x_i into Means and Standard Deviations.
    mu_w_1        = mu_x_1(iiPDF_w)
    mu_w_2        = mu_x_2(iiPDF_w)
    mu_s_1        = mu_x_1(iiPDF_s_mellor)
    mu_s_2        = mu_x_2(iiPDF_s_mellor)
    mu_t_1        = mu_x_1(iiPDF_t_mellor)
    mu_t_2        = mu_x_2(iiPDF_t_mellor)
    mu_rr_1_n     = mu_x_1(iiPDF_rrain)
    mu_rr_2_n     = mu_x_2(iiPDF_rrain)
    mu_Nr_1_n     = mu_x_1(iiPDF_Nr)
    mu_Nr_2_n     = mu_x_2(iiPDF_Nr)
    mu_Ncn_1_n    = mu_x_1(iiPDF_Ncn)
    mu_Ncn_2_n    = mu_x_2(iiPDF_Ncn)
    sigma_w_1     = sigma_x_1(iiPDF_w)
    sigma_w_2     = sigma_x_2(iiPDF_w)
    sigma_s_1     = sigma_x_1(iiPDF_s_mellor)
    sigma_s_2     = sigma_x_2(iiPDF_s_mellor)
    sigma_t_1     = sigma_x_1(iiPDF_t_mellor)
    sigma_t_2     = sigma_x_2(iiPDF_t_mellor)
    sigma_rr_1_n  = sigma_x_1(iiPDF_rrain)
    sigma_rr_2_n  = sigma_x_2(iiPDF_rrain)
    sigma_Nr_1_n  = sigma_x_1(iiPDF_Nr)
    sigma_Nr_2_n  = sigma_x_2(iiPDF_Nr)
    sigma_Ncn_1_n = sigma_x_1(iiPDF_Ncn)
    sigma_Ncn_2_n = sigma_x_2(iiPDF_Ncn)

    ! Unpack variables from hydromet_pdf_params
    mu_rr_1     = hydromet_pdf_params%mu_rr_1
    mu_rr_2     = hydromet_pdf_params%mu_rr_2
    mu_Nr_1     = hydromet_pdf_params%mu_Nr_1
    mu_Nr_2     = hydromet_pdf_params%mu_Nr_2
    mu_Ncn_1    = hydromet_pdf_params%mu_Ncn_1
    mu_Ncn_2    = hydromet_pdf_params%mu_Ncn_2
    sigma_rr_1  = hydromet_pdf_params%sigma_rr_1
    sigma_rr_2  = hydromet_pdf_params%sigma_rr_2
    sigma_Nr_1  = hydromet_pdf_params%sigma_Nr_1
    sigma_Nr_2  = hydromet_pdf_params%sigma_Nr_2
    sigma_Ncn_1 = hydromet_pdf_params%sigma_Ncn_1
    sigma_Ncn_2 = hydromet_pdf_params%sigma_Ncn_2

    rr1           = hydromet_pdf_params%rr1
    rr2           = hydromet_pdf_params%rr2
    Nr1           = hydromet_pdf_params%Nr1
    Nr2           = hydromet_pdf_params%Nr2
    precip_frac   = hydromet_pdf_params%precip_frac
    precip_frac_1 = hydromet_pdf_params%precip_frac_1
    precip_frac_2 = hydromet_pdf_params%precip_frac_2

    ! Unpack corr_array_1 into correlations (1st PDF component).
    corr_st_1     = corr_array_1(iiPDF_t_mellor, iiPDF_s_mellor)
    corr_ws_1     = corr_array_1(iiPDF_w,iiPDF_s_mellor)
    corr_srr_1_n  = corr_array_1(iiPDF_rrain, iiPDF_s_mellor)
    corr_sNr_1_n  = corr_array_1(iiPDF_Nr, iiPDF_s_mellor)
    corr_sNcn_1_n = corr_array_1(iiPDF_Ncn, iiPDF_s_mellor)
    corr_trr_1_n  = corr_array_1(iiPDF_rrain, iiPDF_t_mellor)
    corr_tNr_1_n  = corr_array_1(iiPDF_Nr, iiPDF_t_mellor)
    corr_tNcn_1_n = corr_array_1(iiPDF_Ncn, iiPDF_t_mellor)
    corr_wrr_1_n  = corr_array_1(iiPDF_rrain, iiPDF_w)
    corr_wNr_1_n  = corr_array_1(iiPDF_Nr, iiPDF_w)
    corr_wNcn_1_n = corr_array_1(iiPDF_Ncn, iiPDF_w)
    corr_rrNr_1_n = corr_array_1(iiPDF_Nr, iiPDF_rrain)

    ! Unpack corr_array_2 into correlations (2nd PDF component).
    corr_st_2     = corr_array_2(iiPDF_t_mellor, iiPDF_s_mellor)
    corr_ws_2     = corr_array_2(iiPDF_w,iiPDF_s_mellor)
    corr_srr_2_n  = corr_array_2(iiPDF_rrain, iiPDF_s_mellor)
    corr_sNr_2_n  = corr_array_2(iiPDF_Nr, iiPDF_s_mellor)
    corr_sNcn_2_n = corr_array_2(iiPDF_Ncn, iiPDF_s_mellor)
    corr_trr_2_n  = corr_array_2(iiPDF_rrain, iiPDF_t_mellor)
    corr_tNr_2_n  = corr_array_2(iiPDF_Nr, iiPDF_t_mellor)
    corr_tNcn_2_n = corr_array_2(iiPDF_Ncn, iiPDF_t_mellor)
    corr_wrr_2_n  = corr_array_2(iiPDF_rrain, iiPDF_w)
    corr_wNr_2_n  = corr_array_2(iiPDF_Nr, iiPDF_w)
    corr_wNcn_2_n = corr_array_2(iiPDF_Ncn, iiPDF_w)
    corr_rrNr_2_n = corr_array_2(iiPDF_Nr, iiPDF_rrain)


    return

  end subroutine unpack_pdf_params

!===============================================================================

end module setup_clubb_pdf_params

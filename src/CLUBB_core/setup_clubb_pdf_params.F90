! $Id$
!===============================================================================
module setup_clubb_pdf_params

  implicit none

  private

  public :: setup_pdf_parameters,   &
            unpack_pdf_params,      &
            normalize_pdf_params,   &
            compute_mean_stdev,     &
            compute_corr,           &
            init_precip_hm_arrays

  private :: component_means_hydromet, &
             precip_fraction,          &
             component_mean_hm_ip,     &
             component_stdev_hm_ip,    &
             component_corr_wx,        &
             component_corr_st,        &
             component_corr_whm_ip,    &
             component_corr_xhm_ip,    &
             component_corr_hmxhmy_ip, &
             pdf_param_hm_stats,       &
             pdf_param_log_hm_stats,   &
             precip_hm_arrays,         &
             pack_pdf_params

  ! New hydrometeor arrays are needed which are made out of only the
  ! precipitating hydrometeors, which do not include N_c or N_cn.  CLUBB's PDF
  ! does not include N_c, but rather N_cn.  However, N_cn is handled in the PDF
  ! in a different manner than the rest of the hydrometeors.  The remaining
  ! hydrometeors (the precipitating hydrometeors -- rain, ice, snow, etc.) need
  ! their own array for use in the PDF parameter setup.
  integer, private :: &
    iirr, & ! Index of rain water mixing ratio: precipitating hydrometeor arrays
    iiNr, & ! Index of rain drop concentration: precipitating hydrometeor arrays
    iirs, & ! Index of snow mixing ratio: precipitating hydrometeor arrays
    iiri, & ! Index of ice mixing ratio: precipitating hydrometeor arrays
    iirg, & ! Index of graupel mixing ratio: precipitating hydrometeor arrays
    iiNs, & ! Index of snow concentration: precipitating hydrometeor arrays
    iiNi, & ! Index of ice concentration: precipitating hydrometeor arrays
    iiNg    ! Index of graupel concentration: precipitating hydrometeor arrays
!$omp threadprivate( iirr, iiNr, iirs, iiri, iirg, iiNs, iiNi, iiNg )

  integer, private :: &
    num_hm    ! Number of precipitating hydrometeors
!$omp threadprivate( num_hm )

  ! Prescribed parameters are set to in-cloud or outside-cloud (below-cloud)
  ! values based on whether or not cloud water mixing ratio has a value of at
  ! least rc_tol.  However, this does not take into account the amount of
  ! cloudiness in a component, just whether or not there is any cloud in the
  ! component.  The option l_interp_prescribed_params allows for an interpolated
  ! value between the in-cloud and below-cloud parameter value based on the
  ! component cloud fraction.
  logical, parameter :: &
    l_interp_prescribed_params = .false.

  contains

  !=============================================================================
  subroutine setup_pdf_parameters( nz, hydromet, Ncnm, rho, rcm, &             ! Intent(in)
                                   cloud_frac, w_std_dev, wphydrometp, &       ! Intent(in)
                                   corr_array_cloud, corr_array_below, &       ! Intent(in)
                                   pdf_params, l_stats_samp, d_variables, &    ! Intent(in)
                                   corr_array_1, corr_array_2, &               ! Intent(out)
                                   mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &     ! Intent(out)
                                   corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
                                   hydromet_pdf_params )                       ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zm2zt, & ! Procedure(s)
        gr

    use constants_clubb, only: &
        one,     & ! Constant(s)
        zero,    &
        rc_tol,  &
        rr_tol,  &
        Nr_tol,  &
        Ncn_tol, &
        eps

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use model_flags, only: &
        l_use_precip_frac, & ! Flag(s)
        l_calc_w_corr, &
        l_use_modified_corr

    use advance_windm_edsclrm_module, only: &
        xpwp_fnc

    use variables_diagnostic_module, only: &
        Kh_zm

    use parameters_tunable, only: &
        c_Krrainm

    use clubb_precision, only: &
        core_rknd, &    ! Variable(s)
        dp

    use matrix_operations, only: &
        Cholesky_factor ! Procedure(s)

    use parameters_microphys, only: &
        rrp2_on_rrm2_cloud,       & ! Variable(s)
        rrp2_on_rrm2_below,       &
        Nrp2_on_Nrm2_cloud,       &
        Nrp2_on_Nrm2_below,       &
        ricep2_on_ricem2_cloud,   &
        ricep2_on_ricem2_below,   &
        Nicep2_on_Nicem2_cloud,   &
        Nicep2_on_Nicem2_below,   &
        rsnowp2_on_rsnowm2_cloud, &
        rsnowp2_on_rsnowm2_below, &
        Nsnowp2_on_Nsnowm2_cloud, &
        Nsnowp2_on_Nsnowm2_below

    use stats_type, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        irr1,             & ! Variable(s)
        irr2,             &
        iNr1,             &
        iNr2,             &
        iprecip_frac,     &
        iprecip_frac_1,   &
        iprecip_frac_2,   &
        zt

    use model_flags, only: &
        l_diagnose_correlations ! Variable(s)

    use diagnose_correlations_module, only: &
        diagnose_correlations, & ! Procedure(s)
        calc_cholesky_corr_mtx_approx

    use corr_matrix_module, only: &
        iiPDF_w,        & ! Variable(s)
        iiPDF_s_mellor, &
        iiPDF_t_mellor, &
        iiPDF_rrain,    &
        iiPDF_Nr,       &
        iiPDF_Ncn

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Ncnm,       & ! Mean cloud nuclei concentration, < N_cn >       [num/kg]
      rho,        & ! Density                                         [kg/m^3]
      rcm,        & ! Mean cloud water mixing ratio, < r_c >          [kg/kg]
      cloud_frac, & ! Cloud fraction                                  [-]
      w_std_dev     ! Standard deviation of vertical velocity, w      [m/s]

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

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed correlation cholesky matrix, 1st comp.     [-]
      corr_cholesky_mtx_2    ! Transposed correlation cholesky matrix, 2nd comp.     [-]

    ! Local Variables
    real( kind = dp ), dimension(d_variables,d_variables,nz) :: &
      corr_cholesky_mtx_1_dp, & ! Used for call to Cholesky_factor, which requires dp
      corr_cholesky_mtx_2_dp

    real( kind = core_rknd ), dimension(d_variables,d_variables) :: &
      corr_mtx_approx_1,   & ! Approximated correlation matrix (C = LL'), 1st comp.  [-]
      corr_mtx_approx_2      ! Approximated correlation matrix (C = LL'), 2nd comp.  [-]

    real( kind = core_rknd ), dimension(nz) ::  &
      rrainm, & ! Mean rain water mixing ratio, <r_r> (t-levs.)  [kg/kg]
      Nrm,    & ! Mean rain drop concentration, <N_r> (t-levs.)  [num/kg]
      wprrp,  & ! Covariance of w and r_r, <w'r_r'> (m-levs.)    [(m/s)(kg/kg)]
      wpNrp     ! Covariance of w and N_r, <w'N_r'> (m-levs.)    [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz) :: &
      rc1,         & ! Mean of r_c (1st PDF component)              [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)              [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)           [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)           [-]
      mixt_frac      ! Mixture fraction                             [-]

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
      corr_wrr_1,    & ! Correlation between w and rr (1st PDF component) ip [-]
      corr_wrr_2,    & ! Correlation between w and rr (2nd PDF component) ip [-]
      corr_wNr_1,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      corr_wNr_2,    & ! Correlation between w and Nr (2nd PDF component) ip [-]
      corr_wNcn_1,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      corr_wNcn_2,   & ! Correlation between w and Ncn (2nd PDF component)   [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2,     & ! Correlation between s and t (2nd PDF component)     [-]
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

    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zm,     & ! Covariance of s and w (momentum levels)   [(m/s)(kg/kg)]
      wpNcnp_zm,   & ! Covariance of N_cn and w (momentum levs.) [(m/s)(num/kg)]
      wpsp_zt,     & ! Covariance of s and w on t-levs           [(m/s)(kg/kg)]
      wprrp_ip_zt, & ! Covar. of r_r and w (in-precip) on t-levs [(m/s)(kg/kg)]
      wpNrp_ip_zt, & ! Covar. of N_r and w (in-precip) on t-levs [(m/s)(num/kg)]
      wpNcnp_zt      ! Covariance of N_cn and w on t-levs        [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz,num_hm) :: &
      hmm,   & ! Mean of a precipitating hydrometeor, hm (overall)  [units vary]
      hm1,   & ! Mean of a precip. hydrometeor (1st PDF component)  [units vary]
      hm2,   & ! Mean of a precip. hydrometeor (2nd PDF component)  [units vary]
      wphmp    ! Covariance of w and a precipitating hydrometeor    [(m/s)units]

    real( kind = core_rknd ), dimension(nz,num_hm) :: &
      wphmp_ip_zt    ! Covar. of hm and w (in-precip) on t-levs   [(m/s)(kg/kg)]

    character(len=10), dimension(num_hm) :: &
      hm_list    ! Names of all precipitating hydrometeors

    real( kind = core_rknd ), dimension(num_hm) :: &
      hm_tol    ! Tolerance value for the hydrometeor           [units vary]

    real( kind = core_rknd ), dimension(num_hm) :: &
      hmp2_on_hmm2_cloud, & ! Ratio of < hm'^2 > to < hm >^2 at cloudy levs. [-]
      hmp2_on_hmm2_below    ! Ratio of < hm'^2 > to < hm >^2 at clear levs.  [-]

     real( kind = core_rknd ), dimension(nz) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    real( kind = dp ), dimension(d_variables) :: &
      corr_array_scaling

    logical :: l_corr_array_scaling

    integer :: k, i  ! Loop indices

    ! ---- Begin Code ----

    ! Setup hydrometeor arrays out of precipitating hydrometeors, which do not
    ! include N_c or N_cn.
    call precip_hm_arrays( nz, hydromet, wphydrometp, &
                           hmm, wphmp, hm_list, hm_tol )

    ! Initialize arrays for < hm'^2 > / < hm >^2.
    hmp2_on_hmm2_cloud = one
    hmp2_on_hmm2_below = one

    ! Setup within cloud and outside cloud arrays for < hm'^2 > / < hm >^2.
    if ( iirr > 0 ) then
       hmp2_on_hmm2_cloud(iirr) = rrp2_on_rrm2_cloud
       hmp2_on_hmm2_below(iirr) = rrp2_on_rrm2_below
    endif
    if ( iiNr > 0 ) then
       hmp2_on_hmm2_cloud(iiNr) = Nrp2_on_Nrm2_cloud
       hmp2_on_hmm2_below(iiNr) = Nrp2_on_Nrm2_below
    endif
    if ( iiri > 0 ) then
       hmp2_on_hmm2_cloud(iiri) = ricep2_on_ricem2_cloud
       hmp2_on_hmm2_below(iiri) = ricep2_on_ricem2_below
    endif
    if ( iiNi > 0 ) then
       hmp2_on_hmm2_cloud(iiNi) = Nicep2_on_Nicem2_cloud
       hmp2_on_hmm2_below(iiNi) = Nicep2_on_Nicem2_below
    endif
    if ( iirs > 0 ) then
       hmp2_on_hmm2_cloud(iirs) = rsnowp2_on_rsnowm2_cloud
       hmp2_on_hmm2_below(iirs) = rsnowp2_on_rsnowm2_below
    endif
    if ( iiNs > 0 ) then
       hmp2_on_hmm2_cloud(iiNs) = Nsnowp2_on_Nsnowm2_cloud
       hmp2_on_hmm2_below(iiNs) = Nsnowp2_on_Nsnowm2_below
    endif

    ! Setup some of the PDF parameters
    rc1         = pdf_params%rc1
    rc2         = pdf_params%rc2
    cloud_frac1 = pdf_params%cloud_frac1
    cloud_frac2 = pdf_params%cloud_frac2
    mixt_frac   = pdf_params%mixt_frac

    ! Component mean values for r_r and N_r, and precipitation fraction.
    if ( l_use_precip_frac ) then

       call component_means_hydromet( nz, hmm, rho, rc1, rc2, &
                                      mixt_frac, hm_tol, l_stats_samp, &
                                      hm1, hm2 )

       call precip_fraction( nz, hmm, hm1, hm2, hm_list, hm_tol, &
                             cloud_frac, cloud_frac1, mixt_frac, &
                             precip_frac, precip_frac_1, precip_frac_2 )

    else

       hm1 = hmm
       hm2 = hmm

       precip_frac   = one
       precip_frac_1 = one
       precip_frac_2 = one

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

       wpNcnp_zm(1:nz-1) = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), Ncnm(1:nz-1), &
                                     Ncnm(2:nz), gr%invrs_dzm(1:nz-1) )

       ! Boundary conditions; We are assuming zero flux at the top.
       wpNcnp_zm(nz) = zero

       ! Interpolate the covariances to thermodynamic grid levels.
       wpsp_zt = zm2zt( wpsp_zm )
       do i = 1, num_hm, 1
          wphmp_ip_zt(:,i) = zm2zt( wphmp(:,i) ) / max( precip_frac, eps )
       enddo ! i = 1, num_hm, 1
       wpNcnp_zt = zm2zt( wpNcnp_zm )

       ! When the mean value of a hydrometeor is below tolerance value, it is
       ! considered to have a value of 0, and does not vary over the grid level.
       ! Any covariance involving that hydrometeor also has a value of 0 at that
       ! grid level.
       do k = 1, nz, 1
          do i = 1, num_hm, 1
             if ( hmm(k,i) <= hm_tol(i) ) then
                wphmp_ip_zt(k,i) = zero
             endif
          enddo ! i = 1, num_hm, 1
          if ( Ncnm(k) <= Ncn_tol ) then
             wpNcnp_zt(k) = zero
          endif
       enddo ! k = 1, nz, 1

    endif ! l_calc_w_corr

    ! Initialize the correlation Cholesky matrices
    corr_cholesky_mtx_1 = zero
    corr_cholesky_mtx_2 = zero


    ! Mean of hydrometeor (overall)
    rrainm = hmm(:,iirr)
    Nrm    = hmm(:,iiNr)

    ! Mean of hydrometeor (each PDF component)
    rr1 = hm1(:,iirr)
    rr2 = hm2(:,iirr)
    Nr1 = hm1(:,iiNr)
    Nr2 = hm2(:,iiNr)

    ! Covariance of vertical velocity and a hydrometeor
    ! (< w'r_r' > and < w'N_r' >).
    wprrp = wphmp(:,iirr)
    wpNrp = wphmp(:,iiNr)

    wprrp_ip_zt = wphmp_ip_zt(:,iirr)
    wpNrp_ip_zt = wphmp_ip_zt(:,iiNr)


    ! Statistics
    if ( l_stats_samp ) then

       if ( iirr > 0 ) then

          if ( irr1 > 0 ) then
             ! Mean rain water mixing ratio in PDF component 1.
             call stat_update_var( irr1, rr1, zt )
          endif

          if ( irr2 > 0 ) then
             ! Mean rain water mixing ratio in PDF component 2.
             call stat_update_var( irr2, rr2, zt )
          endif

       endif ! iirr > 0

       if ( iiNr > 0 ) then

          if ( iNr1 > 0 ) then
             ! Mean rain drop concentration in PDF component 1.
             call stat_update_var( iNr1, Nr1, zt )
          endif

          if ( iNr2 > 0 ) then
             ! Mean rain drop concentration in PDF component 2.
             call stat_update_var( iNr2, Nr2, zt )
          endif

       endif ! iiNr > 0

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

    endif


    !!! Setup PDF parameters loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    ! Now also including "model lower boundary" -- Eric Raut Aug 2013
    do k = 1, nz, 1

       !!! Calculate the means and standard deviations involving w, s, t, N_cn,
       !!! and all precipitating hydrometeors (in-precip) in each PDF component.
       call compute_mean_stdev( pdf_params(k), hm1(k,:), hm2(k,:), Ncnm(k), &
                                hmp2_on_hmm2_cloud, hmp2_on_hmm2_below, &
                                rc1(k), rc2(k), cloud_frac1(k), &
                                cloud_frac2(k), precip_frac_1(k), &
                                precip_frac_2(k), hm_tol, d_variables, &
                                mu_x_1(:,k), mu_x_2(:,k), &
                                sigma_x_1(:,k), sigma_x_2(:,k) )

       mu_w_1      = mu_x_1(iiPDF_w,k)
       mu_w_2      = mu_x_2(iiPDF_w,k)
       mu_s_1      = mu_x_1(iiPDF_s_mellor,k)
       mu_s_2      = mu_x_2(iiPDF_s_mellor,k)
       mu_t_1      = mu_x_1(iiPDF_t_mellor,k)
       mu_t_2      = mu_x_2(iiPDF_t_mellor,k)
       mu_rr_1     = mu_x_1(iiPDF_rrain,k)
       mu_rr_2     = mu_x_2(iiPDF_rrain,k)
       mu_Nr_1     = mu_x_1(iiPDF_Nr,k)
       mu_Nr_2     = mu_x_2(iiPDF_Nr,k)
       mu_Ncn_1    = mu_x_1(iiPDF_Ncn,k)
       mu_Ncn_2    = mu_x_2(iiPDF_Ncn,k)
       sigma_w_1   = sigma_x_1(iiPDF_w,k)
       sigma_w_2   = sigma_x_2(iiPDF_w,k)
       sigma_s_1   = sigma_x_1(iiPDF_s_mellor,k)
       sigma_s_2   = sigma_x_2(iiPDF_s_mellor,k)
       sigma_t_1   = sigma_x_1(iiPDF_t_mellor,k)
       sigma_t_2   = sigma_x_2(iiPDF_t_mellor,k)
       sigma_rr_1  = sigma_x_1(iiPDF_rrain,k)
       sigma_rr_2  = sigma_x_2(iiPDF_rrain,k)
       sigma_Nr_1  = sigma_x_1(iiPDF_Nr,k)
       sigma_Nr_2  = sigma_x_2(iiPDF_Nr,k)
       sigma_Ncn_1 = sigma_x_1(iiPDF_Ncn,k)
       sigma_Ncn_2 = sigma_x_2(iiPDF_Ncn,k)

       if ( l_diagnose_correlations ) then

          if ( rcm(k) > rc_tol ) then

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_ws_1, corr_wrr_1, &         ! Intent(out)
                                         corr_wNr_1, corr_wNcn_1, &       ! Intent(out)
                                         corr_st_1, corr_srr_1, &         ! Intent(out)
                                         corr_sNr_1, corr_sNcn_1, &       ! Intent(out)
                                         corr_trr_1, corr_tNr_1, &        ! Intent(out)
                                         corr_tNcn_1, corr_rrNr_1 )       ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_ws_2, corr_wrr_2, &         ! Intent(out)
                                         corr_wNr_2, corr_wNcn_2, &       ! Intent(out)
                                         corr_st_2, corr_srr_2, &         ! Intent(out)
                                         corr_sNr_2, corr_sNcn_2, &       ! Intent(out)
                                         corr_trr_2, corr_tNr_2, &        ! Intent(out)
                                         corr_tNcn_2, corr_rrNr_2 )       ! Intent(out)

          else

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_ws_1, corr_wrr_1, &         ! Intent(out)
                                         corr_wNr_1, corr_wNcn_1, &       ! Intent(out)
                                         corr_st_1, corr_srr_1, &         ! Intent(out)
                                         corr_sNr_1, corr_sNcn_1, &       ! Intent(out)
                                         corr_trr_1, corr_tNr_1, &        ! Intent(out)
                                         corr_tNcn_1, corr_rrNr_1 )       ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_ws_2, corr_wrr_2, &         ! Intent(out)
                                         corr_wNr_2, corr_wNcn_2, &       ! Intent(out)
                                         corr_st_2, corr_srr_2, &         ! Intent(out)
                                         corr_sNr_2, corr_sNcn_2, &       ! Intent(out)
                                         corr_trr_2, corr_tNr_2, &        ! Intent(out)
                                         corr_tNcn_2, corr_rrNr_2 )       ! Intent(out)

          endif

       else ! if .not. l_diagnose_correlations

          call compute_corr( Ncnm(k), rr1(k), rr2(k), Nr1(k), Nr2(k), rc1(k), & ! Intent(in)
                             rc2(k), cloud_frac1(k), cloud_frac2(k), &          ! Intent(in)
                             wpsp_zt(k), wprrp_ip_zt(k), wpNrp_ip_zt(k), &      ! Intent(in)
                             wpNcnp_zt(k), w_std_dev(k), &                      ! Intent(in)
                             sigma_rr_1, sigma_Nr_1, sigma_Ncn_1, &             ! Intent(in)
                             pdf_params(k), &                                   ! Intent(in)
                             corr_ws_1, corr_ws_2, corr_wrr_1, &                ! Intent(out)
                             corr_wrr_2, corr_wNr_1, corr_wNr_2, &              ! Intent(out)
                             corr_wNcn_1, corr_wNcn_2, corr_st_1, &             ! Intent(out)
                             corr_st_2, corr_srr_1, corr_srr_2, &               ! Intent(out)
                             corr_sNr_1, corr_sNr_2, corr_sNcn_1, &             ! Intent(out)
                             corr_sNcn_2, corr_trr_1, corr_trr_2, &             ! Intent(out)
                             corr_tNr_1, corr_tNr_2, corr_tNcn_1, &             ! Intent(out)
                             corr_tNcn_2, corr_rrNr_1, corr_rrNr_2 )            ! Intent(out)

       endif ! l_diagnose_correlations

       !!! Calculate the mean, standard deviations, and correlations involving
       !!! ln r_r, ln N_r, and ln N_cn for each PDF component.
       call normalize_pdf_params( rr1(k), rr2(k), Nr1(k), Nr2(k), Ncnm(k), &     ! Intent(in)
                                   mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &         ! Intent(in)
                                   mu_Ncn_1, mu_Ncn_2, sigma_rr_1, sigma_rr_2, & ! Intent(in)
                                   sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, &        ! Intent(in)
                                   sigma_Ncn_2, corr_wrr_1, corr_wrr_2, &        ! Intent(in)
                                   corr_wNr_1, corr_wNr_2, corr_wNcn_1, &        ! Intent(in)
                                   corr_wNcn_2, corr_srr_1, corr_srr_2, &        ! Intent(in)
                                   corr_sNr_1, corr_sNr_2, corr_sNcn_1, &        ! Intent(in)
                                   corr_sNcn_2, corr_trr_1, corr_trr_2, &        ! Intent(in)
                                   corr_tNr_1, corr_tNr_2, corr_tNcn_1, &        ! Intent(in)
                                   corr_tNcn_2, corr_rrNr_1, corr_rrNr_2, &      ! Intent(in)
                                   mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(out)
                                   mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(out)
                                   sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(out)
                                   sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(out)
                                   corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(out)
                                   corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(out)
                                   corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(out)
                                   corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(out)
                                   corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(out)
                                   corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(out)
                                   corr_rrNr_1_n, corr_rrNr_2_n )                ! Intent(out)

       !!! Statistics
       call pdf_param_hm_stats(  mu_rr_1, mu_rr_2, mu_Nr_1, &
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
                                 k, l_stats_samp )

       call pdf_param_log_hm_stats( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                    mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                    sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                    sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                    corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &
                                    corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &
                                    corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                                    corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                                    corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &
                                    corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &
                                    corr_rrNr_1_n, corr_rrNr_2_n, k, &
                                    l_stats_samp )

       !!! Pack the PDF parameters
       call pack_pdf_params( mu_w_1, mu_w_2, mu_s_1, mu_s_2, &                     ! Intent(in)
                             mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &                   ! Intent(in)
                             mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &               ! Intent(in)
                             mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &                    ! Intent(in)
                             mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &                  ! Intent(in)
                             sigma_w_1, sigma_w_2, sigma_s_1, &                    ! Intent(in)
                             sigma_s_2, sigma_t_1, sigma_t_2, &                    ! Intent(in)
                             sigma_rr_1, sigma_rr_2, sigma_Nr_1, &                 ! Intent(in)
                             sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &               ! Intent(in)
                             sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &           ! Intent(in)
                             sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &         ! Intent(in)
                             corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, &         ! Intent(in)
                             corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &           ! Intent(in)
                             corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &         ! Intent(in)
                             corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &           ! Intent(in)
                             corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &         ! Intent(in)
                             corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &           ! Intent(in)
                             corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &         ! Intent(in)
                             corr_rrNr_1_n, corr_rrNr_2_n, &                       ! Intent(in)
                             rr1(k), rr2(k), Nr1(k), Nr2(k), &                     ! Intent(in)
                             precip_frac(k), precip_frac_1(k), precip_frac_2(k), & ! Intent(in)
                             d_variables, &                                        ! Intent(in)
                             corr_array_1(:,:,k), corr_array_2(:,:,k), &           ! Intent(out)
                             mu_x_1(:,k), mu_x_2(:,k), &                           ! Intent(out)
                             sigma_x_1(:,k), sigma_x_2(:,k), &                     ! Intent(out)
                             hydromet_pdf_params(k) )                              ! Intent(out)

       if ( l_use_modified_corr ) then

          if ( l_diagnose_correlations ) then

             call calc_cholesky_corr_mtx_approx &
                            ( d_variables, corr_array_1(:,:,k), &                    ! intent(in)
                              corr_cholesky_mtx_1(:,:,k), corr_mtx_approx_1 ) ! intent(out)

             call calc_cholesky_corr_mtx_approx &
                            ( d_variables, corr_array_2(:,:,k), &                    ! intent(in)
                              corr_cholesky_mtx_2(:,:,k), corr_mtx_approx_2 ) ! intent(out)

             corr_array_1(:,:,k) = corr_mtx_approx_1
             corr_array_2(:,:,k) = corr_mtx_approx_2

          else

             ! Compute choleksy factorization for the correlation matrix (out of cloud)
             call Cholesky_factor( d_variables, real(corr_array_1(:,:,k), kind = dp), & ! In
                                   corr_array_scaling, corr_cholesky_mtx_1_dp(:,:,k), &  ! Out
                                   l_corr_array_scaling ) ! Out

             call Cholesky_factor( d_variables, real(corr_array_2(:,:,k), kind = dp), & ! In
                                   corr_array_scaling, corr_cholesky_mtx_2_dp(:,:,k), &  ! Out
                                   l_corr_array_scaling ) ! Out
             corr_cholesky_mtx_1(:,:,k) = real( corr_cholesky_mtx_1_dp(:,:,k), kind = core_rknd )
             corr_cholesky_mtx_2(:,:,k) = real( corr_cholesky_mtx_2_dp(:,:,k), kind = core_rknd )
          endif

       endif

    enddo  ! Setup PDF parameters loop: k = 2, nz, 1


    return

  end subroutine setup_pdf_parameters

  !=============================================================================
  subroutine component_means_hydromet( nz, hmm, rho, rc1, rc2, &
                                       mixt_frac, hm_tol, l_stats_samp, &
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

    real( kind = core_rknd ), dimension(nz,num_hm), intent(in) :: &
      hmm    ! Mean of hydrometeor, hm (overall)                [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho,       & ! Air density                                        [kg/m^3]
      rc1,       & ! Mean cloud water mixing ratio (1st PDF component)  [kg/kg]
      rc2,       & ! Mean cloud water mixing ratio (2nd PDF component)  [kg/kg]
      mixt_frac    ! Mixture fraction                                   [-]

    real( kind = core_rknd ), dimension(num_hm), intent(in) :: &
      hm_tol    ! Tolerance value for the hydrometeor           [units vary]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,num_hm), intent(out) :: &
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
    !!! hm2/hm1 ( = rr2/rr1 = Nr2/Nr1, etc. ) = LWP2/LWP1.
    do i = 1, num_hm, 1

       do k = 1, nz, 1

          !!! Calculate the component means for the hydrometeor.
          if ( hmm(k,i) > hm_tol(i) ) then

             if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

                ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
                !
                ! The hydrometeor is found at this level, yet there is no cloud
                ! at or above the current level.  This is usually due to a
                ! numerical artifact.  For example, the hydrometeor is diffused
                ! above cloud top.  Simply set each component mean equal to the
                ! overall mean.
                hm1(k,i) = hmm(k,i)
                hm2(k,i) = hmm(k,i)

             elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

                ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
                ! insignificant amount).
                !
                ! The hydrometeor is found at this level, and all cloud water at
                ! or above this level is found in the 1st PDF component.  All of
                ! the hydrometeor is found in the 1st PDF component.
                hm1(k,i) = hmm(k,i) / mixt_frac(k)
                hm2(k,i) = zero

             elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

                ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
                ! insignificant amount).
                !
                ! The hydrometeor is found at this level, and all cloud water at
                ! or above this level is found in the 2nd PDF component.  All of
                ! the hydrometeor is found in the 2nd PDF component.
                hm1(k,i) = zero
                hm2(k,i) = hmm(k,i) / ( one - mixt_frac(k) )

             else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

                ! Both LWP1 and LWP2 are (significantly) greater than 0.
                !
                ! The hydrometeor is found at this level, and there is
                ! sufficient cloud water at or above this level in both PDF
                ! components to find the hydrometeor in both PDF components.
                ! Delegate the hydrometeor between the 1st and 2nd PDF
                ! components according to the above equations.
                hm1(k,i) &
                = hmm(k,i) &
                  / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

                hm2(k,i) &
                = ( hmm(k,i) - mixt_frac(k) * hm1(k,i) ) &
                  / ( one - mixt_frac(k) )

             endif


          else ! hmm(k) <= hm_tol

             ! The overall hydrometeor is either 0 or below tolerance value (any
             ! postive value is considered to be a numerical artifact).  Simply
             ! set each component mean equal to the overall mean.
             hm1(k,i) = hmm(k,i)
             hm2(k,i) = hmm(k,i)

          endif

       enddo ! k = 1, nz, 1

    enddo ! i = 1, num_hm, 1


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
  subroutine precip_fraction( nz, hmm, hm1, hm2, hm_list, hm_tol, &
                              cloud_frac, cloud_frac1, mixt_frac, &
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

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,num_hm), intent(in) :: &
      hmm, & ! Mean of hydrometeor, hm (overall)                [units vary]
      hm1, & ! Mean of hydrometeor (1st PDF component)          [units vary]
      hm2    ! Mean of hydrometeor (2nd PDF component)          [units vary]

    character(len=10), dimension(num_hm), intent(in) :: &
      hm_list    ! Names of all hydrometeors

    real( kind = core_rknd ), dimension(num_hm), intent(in) :: &
      hm_tol    ! Tolerance value for the hydrometeor           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,  & ! Cloud fraction (overall)                         [-] 
      cloud_frac1, & ! Cloud fraction (1st PDF component)               [-]
      mixt_frac      ! Mixture fraction                                 [-]

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

    real( kind = core_rknd ) :: &
      r_tot_hm_tol    ! Tolerance value for total hydromet mixing ratio [kg/kg]

    character(len=10) :: &
      hydromet_name    ! Name of a hydrometeor

    real( kind = core_rknd ), parameter :: &
      precip_frac_tol = cloud_frac_min  ! Minimum precip. frac.         [-]
    
    integer :: &
      k, i   ! Loop indices


    !!! Find overall precipitation fraction.
    do k = nz, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       if ( k < nz ) then
          precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )
       else  ! k = nz
          precip_frac(k) = cloud_frac(k)
       endif

       if ( any( hmm(k,:) > hm_tol(:) ) &
            .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find any hydrometeor at this grid level, but
          ! no cloud at or above this grid level, set precipitation fraction to
          ! a minimum threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( all( hmm(k,:) <= hm_tol(:) ) &
                .and. precip_frac(k) < precip_frac_tol ) then

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There is not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo ! Overall precipitation fraction loop: k = nz, 1, -1.


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

          elseif ( any( hm1(k,:) > hm_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 1st PDF
             ! component at this grid level, but no cloud in the 1st PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 1st PDF component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( all( hm1(k,:) <= hm_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 1st PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There is not any
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
             elseif ( any( hm1(k,:) > hm_tol(:) ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = precip_frac_tol
             elseif ( all( hm1(k,:) <= hm_tol(:) ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = zero
             endif

          elseif ( any( hm2(k,:) > hm_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 2nd PDF
             ! component at this grid level, but no cloud in the 2nd PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 2nd PDF component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( all( hm2(k,:) <= hm_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 2nd PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There is not any
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

          if ( all( hm1(k,:) <= hm_tol(:) ) &
               .and. all( hm2(k,:) <= hm_tol(:) ) ) then

             ! There are no hydrometeors found in each PDF component.
             ! Precipitation fraction within each component is set to 0.
             precip_frac_1(k) = zero
             precip_frac_2(k) = zero

          elseif ( any( hm1(k,:) > hm_tol(:) ) &
                   .and. all( hm2(k,:) <= hm_tol(:) ) ) then

             ! All the hydrometeors are found within the 1st PDF component.
             precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
             precip_frac_2(k) = zero

          elseif ( any( hm2(k,:) > hm_tol(:) ) &
                   .and. all( hm1(k,:) <= hm_tol(:) ) ) then

             ! All the hydrometeors are found within the 2nd PDF component.
             precip_frac_1(k) = zero
             precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )

          else ! any( hm1(k,:) > hm_tol(:) ) AND any( hm2(k,:) > hm_tol(:) )

             ! Hydrometeors are found within both PDF components.
             r_tot_hm_1 = zero
             r_tot_hm_2 = zero
             N_tot_hm_1 = zero
             N_tot_hm_2 = zero
             r_tot_hm_tol = one
             do i = 1, num_hm, 1

                hydromet_name = hm_list(i)

                if ( hydromet_name(1:1) == "r" ) then

                   ! The hydrometeor is a mixing ratio.
                   ! Find total hydrometeor mixing ratio in each PDF component.
                   r_tot_hm_1 = r_tot_hm_1 + hm1(k,i)
                   r_tot_hm_2 = r_tot_hm_2 + hm2(k,i)

                   if ( hm1(k,i) > hm_tol(i) ) then
                      r_tot_hm_tol = min( r_tot_hm_tol, hm_tol(i) )
                   endif

                elseif ( hydromet_name(1:1) == "N" ) then

                   ! The hydrometeor is a concentration.
                   ! Find total hydrometeor concentration in each PDF component.
                   N_tot_hm_1 = N_tot_hm_1 + hm1(k,i)
                   N_tot_hm_2 = N_tot_hm_2 + hm2(k,i)

                endif

             enddo ! i = 1, num_hm, 1

             !!! Find precipitation fraction within PDF component 1.
             if ( r_tot_hm_1 > r_tot_hm_tol ) then
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) &
                      + ( one - mixt_frac(k) ) * r_tot_hm_2/r_tot_hm_1 )
             else ! N_tot_hm_1 > N_tot_hm_tol
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
          if ( any( hm1(k,:) > hm_tol(:) ) &
               .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 1st PDF
             ! component at this grid level, but no cloud in the 1st PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 1st PDF component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( all( hm1(k,:) <= hm_tol(:) ) &
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
          if ( any( hm2(k,:) > hm_tol(:) ) &
               .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 2nd PDF
             ! component at this grid level, but no cloud in the 2nd PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 2nd PDF component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( all( hm2(k,:) <= hm_tol(:) ) &
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


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine compute_mean_stdev( pdf_params, hm1, hm2, Ncnm, &
                                 hmp2_on_hmm2_cloud, hmp2_on_hmm2_below, &
                                 rc1, rc2, cloud_frac1, &
                                 cloud_frac2, precip_frac_1, &
                                 precip_frac_2, hm_tol, d_variables, &
                                 mu_x_1, mu_x_2, &
                                 sigma_x_1, sigma_x_2 )
       
    ! Description:
    ! Setup the PDF mean vector and the PDF standard deviation vector for each
    ! PDF component.  PDF variables that are precipitating hydrometeors use an
    ! in-precip mean and an in-precip standard deviation.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,     & ! Constant(s)
        zero,    &
        Ncn_tol

    use corr_matrix_module, only: &
        iiPDF_w, & ! Constant(s)
        iiPDF_s => iiPDF_s_mellor, &
        iiPDF_t => iiPDF_t_mellor, &
        iiPDF_Ncn

    use parameters_microphys, only: &
        Ncnp2_on_Ncnm2_cloud  ! Variable(s) 

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    implicit none

    ! Input Variables
    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    real( kind = core_rknd ), dimension(num_hm), intent(in) :: &
      hm1, & ! Mean of a precipitating hydrometeor (1st PDF comp.)  [units vary]
      hm2    ! Mean of a precipitating hydrometeor (2nd PDF comp.)  [units vary]

    real( kind = core_rknd ), intent(in) :: &
      Ncnm    ! Mean cloud nuclei concentration                     [num/kg]

    real( kind = core_rknd ), dimension(num_hm), intent(in) :: &
      hmp2_on_hmm2_cloud, & ! Ratio of < hm'^2 > to < hm >^2 at cloudy levs. [-]
      hmp2_on_hmm2_below    ! Ratio of < hm'^2 > to < hm >^2 at clear levs.  [-]

    real( kind = core_rknd ), intent(in) :: &
      rc1,           & ! Mean of r_c (1st PDF component)                [kg/kg]
      rc2,           & ! Mean of r_c (2nd PDF component)                [kg/kg]
      cloud_frac1,   & ! Cloud fraction (1st PDF component)             [-]
      cloud_frac2,   & ! Cloud fraction (2nd PDF component)             [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    real( kind = core_rknd ), dimension(num_hm), intent(in) :: &
      hm_tol    ! Tolerance value for the hydrometeor               [units vary]

    integer, intent(in) :: &
      d_variables    ! Number of PDF variables

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1,    & ! Means of PDF-variables (1st PDF component)     [units vary]
      mu_x_2,    & ! Means of PDF-variables (2nd PDF component)     [units vary]
      sigma_x_1, & ! Standard devs. of PDF-vars. (1st PDF comp.)    [units vary]
      sigma_x_2    ! Standard devs. of PDF-vars. (2nd PDF comp.)    [units vary]

    ! Local Variable
    integer :: i    ! Loop index


    !!! Enter the PDF parameters.

    !!! Means.

    ! Vertical Velocity, w.

    ! Mean of vertical velocity, w, in PDF component 1.
    ! Note:  mu_x_1(iiPDF_w) and mu_w_1 denote the same thing.
    mu_x_1(iiPDF_w) = pdf_params%w1

    ! Mean of vertical velocity, w, in PDF component 2.
    ! Note:  mu_x_2(iiPDF_w) and mu_w_2 denote the same thing.
    mu_x_2(iiPDF_w) = pdf_params%w2


    ! Extended liquid water mixing ratio, s.

    ! Mean of extended liquid water mixing ratio, s, in PDF component 1.
    ! Note:  mu_x_1(iiPDF_s) and mu_s_1 denote the same thing.
    mu_x_1(iiPDF_s) = pdf_params%s1

    ! Mean of extended liquid water mixing ratio, s, in PDF component 2.
    ! Note:  mu_x_2(iiPDF_s) and mu_s_2 denote the same thing.
    mu_x_2(iiPDF_s) = pdf_params%s2


    ! PDF variable t (where s and t replace r_t and theta_l in the PDF).

    ! Mean of t in PDF component 1.
    ! Set the component mean values of t to 0.
    ! The component mean values of t are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    ! Note:  mu_x_1(iiPDF_t) and mu_t_1 denote the same thing.
    mu_x_1(iiPDF_t) = zero

    ! Mean of t in PDF component 2.
    ! Set the component mean values of t to 0.
    ! The component mean values of t are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    ! Note:  mu_x_2(iiPDF_t) and mu_t_2 denote the same thing.
    mu_x_2(iiPDF_t) = zero


    ! Cloud nuclei concentration, N_cn.

    ! Mean of cloud nuclei concentration in PDF component 1.
    ! Note:  mu_x_1(iiPDF_Ncn) and mu_Ncn_1 denote the same thing.
    ! Note:  Ncn does not use an in-precip mean, but rather an overall mean,
    !        so 1 is sent in for precip_frac in this function call.
    mu_x_1(iiPDF_Ncn) = component_mean_hm_ip( Ncnm, one, Ncn_tol )

    ! Mean of cloud nuclei concentration in PDF component 2.
    ! Note:  mu_x_2(iiPDF_Ncn) and mu_Ncn_2 denote the same thing.
    ! Note:  Ncn does not use an in-precip mean, but rather an overall mean,
    !        so 1 is sent in for precip_frac in this function call.
    mu_x_2(iiPDF_Ncn) = component_mean_hm_ip( Ncnm, one, Ncn_tol )


    ! Precipitating hydrometeors, hm.
    if ( num_hm > 0 ) then

       do i = 1, num_hm, 1

          ! Mean of hydrometeor (in-precip) in PDF component 1.
          ! Note:  mu_x_1(iiPDF_rr) and mu_rr_1 denote the same thing;
          !        mu_x_1(iiPDF_Nr) and mu_Nr_1 denote the same thing; etc.
          mu_x_1(i+4) &
          = component_mean_hm_ip( hm1(i), precip_frac_1, hm_tol(i) )

          ! Mean of hydrometeor (in-precip) in PDF component 2.
          ! Note:  mu_x_2(iiPDF_rr) and mu_rr_2 denote the same thing;
          !        mu_x_2(iiPDF_Nr) and mu_Nr_2 denote the same thing; etc.
          mu_x_2(i+4) &
          = component_mean_hm_ip( hm2(i), precip_frac_2, hm_tol(i) )

       enddo ! i = 1, num_hm, 1

    endif ! num_hm > 0


    !!! Standard deviations.

    ! Vertical Velocity, w.

    ! Standard deviation of vertical velocity, w, in PDF component 1.
    ! Note:  sigma_x_1(iiPDF_w) and sigma_w_1 denote the same thing.
    sigma_x_1(iiPDF_w) = sqrt( pdf_params%varnce_w1 )

    ! Standard deviation of vertical velocity, w, in PDF component 2.
    ! Note:  sigma_x_2(iiPDF_w) and sigma_w_2 denote the same thing.
    sigma_x_2(iiPDF_w) = sqrt( pdf_params%varnce_w2 )


    ! Extended liquid water mixing ratio, s.

    ! Standard deviation of extended liquid water mixing ratio, s,
    ! in PDF component 1.
    ! Note:  sigma_x_1(iiPDF_s) and sigma_s_1 denote the same thing.
    sigma_x_1(iiPDF_s) = pdf_params%stdev_s1

    ! Standard deviation of extended liquid water mixing ratio, s,
    ! in PDF component 2.
    ! Note:  sigma_x_2(iiPDF_s) and sigma_s_2 denote the same thing.
    sigma_x_2(iiPDF_s) = pdf_params%stdev_s2


    ! PDF variable t (where s and t replace r_t and theta_l in the PDF).

    ! Standard deviation of t in PDF component 1.
    ! Note:  sigma_x_1(iiPDF_t) and sigma_t_1 denote the same thing.
    sigma_x_1(iiPDF_t) = pdf_params%stdev_t1

    ! Standard deviation of t in PDF component 2.
    ! Note:  sigma_x_2(iiPDF_t) and sigma_t_2 denote the same thing.
    sigma_x_2(iiPDF_t) = pdf_params%stdev_t2

    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    ! Brian Griffin; February 3, 2007.

    ! Cloud nuclei concentration, N_cn.

    ! Standard deviation of cloud nuclei concentration in PDF component 1.
    ! Note:  sigma_x_1(iiPDF_Ncn) and sigma_Ncn_1 denote the same thing.
    ! Note:  Ncn does not use an in-precip mean, but rather an overall mean,
    !        so 1 is sent in for cloud_frac in this function call.
    sigma_x_1(iiPDF_Ncn) &
    = component_stdev_hm_ip( Ncnm, mu_x_1(iiPDF_Ncn), rc1, one, Ncn_tol, &
                             Ncnp2_on_Ncnm2_cloud, Ncnp2_on_Ncnm2_cloud )

    ! Standard deviation of cloud nuclei concentration in PDF component 2.
    ! Note:  sigma_x_2(iiPDF_Ncn) and sigma_Ncn_2 denote the same thing.
    ! Note:  Ncn does not use an in-precip mean, but rather an overall mean,
    !        so 1 is sent in for cloud_frac in this function call.
    sigma_x_2(iiPDF_Ncn) &
    = component_stdev_hm_ip( Ncnm, mu_x_2(iiPDF_Ncn), rc2, one, Ncn_tol, &
                             Ncnp2_on_Ncnm2_cloud, Ncnp2_on_Ncnm2_cloud )


    ! Precipitating hydrometeors, hm.
    if ( num_hm > 0 ) then

       do i = 1, num_hm, 1

          ! Standard deviation of hydrometeor (in-precip) in PDF component 1.
          ! Note:  sigma_x_1(iiPDF_rr) and sigma_rr_1 denote the same thing;
          !        sigma_x_1(iiPDF_Nr) and sigma_Nr_1 denote the same thing;
          !        etc.
          sigma_x_1(i+4) &
          = component_stdev_hm_ip( hm1(i), mu_x_1(i+4), rc1, cloud_frac1, hm_tol(i), &
                                   hmp2_on_hmm2_cloud(i), hmp2_on_hmm2_below(i) )

          ! Standard deviation of hydrometeor (in-precip) in PDF component 2.
          ! Note:  sigma_x_2(iiPDF_rr) and sigma_rr_2 denote the same thing;
          !        sigma_x_2(iiPDF_Nr) and sigma_Nr_2 denote the same thing;
          !        etc.
          sigma_x_2(i+4) &
          = component_stdev_hm_ip( hm2(i), mu_x_2(i+4), rc2, cloud_frac2, hm_tol(i), &
                                   hmp2_on_hmm2_cloud(i), hmp2_on_hmm2_below(i) )

       enddo ! i = 1, num_hm, 1

    endif ! num_hm > 0


    return

  end subroutine compute_mean_stdev

!=============================================================================
  subroutine compute_corr( Ncnm, rr1, rr2, Nr1, Nr2, rc1, &        ! Intent(in)
                           rc2, cloud_frac1, cloud_frac2, &        ! Intent(in)
                           wpsp, wprrp_ip, wpNrp_ip, &             ! Intent(in)
                           wpNcnp, stdev_w, &                      ! Intent(in)
                           sigma_rr_1, sigma_Nr_1, sigma_Ncn_1, &  ! Intent(in)
                           pdf_params, &                           ! Intent(in)
                           corr_ws_1, corr_ws_2, corr_wrr_1, &     ! Intent(out)
                           corr_wrr_2, corr_wNr_1, corr_wNr_2, &   ! Intent(out)
                           corr_wNcn_1, corr_wNcn_2, corr_st_1, &  ! Intent(out)
                           corr_st_2, corr_srr_1, corr_srr_2, &    ! Intent(out)
                           corr_sNr_1, corr_sNr_2, corr_sNcn_1, &  ! Intent(out)
                           corr_sNcn_2, corr_trr_1, corr_trr_2, &  ! Intent(out)
                           corr_tNr_1, corr_tNr_2, corr_tNcn_1, &  ! Intent(out)
                           corr_tNcn_2, corr_rrNr_1, corr_rrNr_2 ) ! Intent(out)
       
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol,       & ! Constant(s)
        rr_tol,       &
        Nr_tol,       &
        Ncn_tol,      &
        w_tol,        & ! [m/s]
        s_mellor_tol, & ! [kg/kg]
        one,          &
        zero

    use parameters_microphys, only: &
        l_fix_s_t_correlations ! Variable(s)

    use fixed_correlations, only: &
        corr_sw_NN_cloud,   & ! Variable(s)
        corr_wrr_NL_cloud,  &
        corr_wNr_NL_cloud,  &
        corr_wNcn_NL_cloud, &
        corr_st_NN_cloud,   &
        corr_srr_NL_cloud,  &
        corr_sNr_NL_cloud,  &
        corr_sNcn_NL_cloud, &
        corr_trr_NL_cloud,  &
        corr_tNr_NL_cloud,  &
        corr_tNcn_NL_cloud, &
        corr_rrNr_LL_cloud

    use fixed_correlations, only: &
        corr_sw_NN_below,   & ! Variable(s)
        corr_wrr_NL_below,  &
        corr_wNr_NL_below,  &
        corr_wNcn_NL_below, &
        corr_st_NN_below,   &
        corr_srr_NL_below,  &
        corr_sNr_NL_below,  &
        corr_sNcn_NL_below, &
        corr_trr_NL_below,  &
        corr_tNr_NL_below,  &
        corr_tNcn_NL_below, &
        corr_rrNr_LL_below

    use model_flags, only: &
        l_calc_w_corr, &
        l_use_hydromet_tolerance, &
        l_use_modified_corr

    use diagnose_correlations_module, only: &
        calc_mean,        & ! Procedure(s)
        calc_w_corr

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Ncnm,        & ! Mean cloud nuclei concentration                [num/kg]
      rr1,         & ! Mean rain water mixing ratio (1st PDF comp.)    [kg/kg]
      rr2,         & ! Mean rain water mixing ratio (2nd PDF comp.)    [kg/kg]
      Nr1,         & ! Mean rain drop concentration (1st PDF comp.)   [num/kg]
      Nr2,         & ! Mean rain drop concentration (2nd PDF comp.)   [num/kg]
      rc1,         & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)                  [-]
      wpsp,        & ! Covariance of w and s                      [(m/s)kg/kg]
      wprrp_ip,    & ! Covariance of w and r_r (overall) ip       [(m/s)kg/kg]
      wpNrp_ip,    & ! Covariance of w and N_r (overall) ip      [(m/s)num/kg]
      wpNcnp,      & ! Covariance of w and N_cn                  [(m/s)num/kg]
      stdev_w        ! Standard deviation of w                           [m/s]

    real( kind = core_rknd ), intent(in) :: &
      sigma_rr_1,  &
      sigma_Nr_1,  &
      sigma_Ncn_1

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      corr_ws_1,   & ! Correlation between w and s (1st PDF component)       [-]
      corr_ws_2,   & ! Correlation between w and s (2nd PDF component)       [-]
      corr_wrr_1,  & ! Correlation between w and rr (1st PDF component) ip   [-]
      corr_wrr_2,  & ! Correlation between w and rr (2nd PDF component) ip   [-]
      corr_wNr_1,  & ! Correlation between w and Nr (1st PDF component) ip   [-]
      corr_wNr_2,  & ! Correlation between w and Nr (2nd PDF component) ip   [-]
      corr_wNcn_1, & ! Correlation between w and Ncn (1st PDF component)     [-]
      corr_wNcn_2, & ! Correlation between w and Ncn (2nd PDF component)     [-]
      corr_st_1,   & ! Correlation between s and t (1st PDF component)       [-]
      corr_st_2,   & ! Correlation between s and t (2nd PDF component)       [-]
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component) ip   [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component) ip   [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component) ip   [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component) ip   [-]
      corr_sNcn_1, & ! Correlation between s and Ncn (1st PDF component)     [-]
      corr_sNcn_2, & ! Correlation between s and Ncn (2nd PDF component)     [-]
      corr_trr_1,  & ! Correlation between t and rr (1st PDF component) ip   [-]
      corr_trr_2,  & ! Correlation between t and rr (2nd PDF component) ip   [-]
      corr_tNr_1,  & ! Correlation between t and Nr (1st PDF component) ip   [-]
      corr_tNr_2,  & ! Correlation between t and Nr (2nd PDF component) ip   [-]
      corr_tNcn_1, & ! Correlation between t and Ncn (1st PDF component)     [-]
      corr_tNcn_2, & ! Correlation between t and Ncn (2nd PDF component)     [-]
      corr_rrNr_1, & ! Correlation between rr and Nr (1st PDF component) ip  [-]
      corr_rrNr_2    ! Correlation between rr and Nr (2nd PDF component) ip  [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      s_mellor_m,     & ! Mean of s_mellor                            [kg/kg]
      stdev_s_mellor, & ! Standard deviation of s_mellor              [kg/kg]
      corr_ws,        & ! Correlation between w and s (overall)       [-]
      corr_wrr,       & ! Correlation between w and rr (overall) ip   [-]
      corr_wNr,       & ! Correlation between w and Nr (overall) ip   [-]
      corr_wNcn         ! Correlation between w and Ncn (overall)     [-]

    logical :: &
      l_limit_corr_st    ! If true, we limit the correlation between s and t [-]

    !!! Enter the PDF parameters.

    !!! Correlations

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

       corr_wrr = calc_w_corr( wprrp_ip, stdev_w, sigma_rr_1, w_tol, rr_tol )

       corr_wNr = calc_w_corr( wpNrp_ip, stdev_w, sigma_Nr_1, w_tol, Nr_tol )

       corr_wNcn = calc_w_corr( wpNcnp, stdev_w, sigma_Ncn_1, w_tol, Ncn_tol )

    endif

    ! Correlation between w and s in PDF component 1.
    corr_ws_1 &
    = component_corr_wx( corr_ws, rc1, cloud_frac1, &
                         corr_sw_NN_cloud, corr_sw_NN_below )

    ! Correlation between w and s in PDF component 2.
    corr_ws_2 &
    = component_corr_wx( corr_ws, rc2, cloud_frac2, &
                         corr_sw_NN_cloud, corr_sw_NN_below )

    ! Correlation (in-precip) between w and r_r in PDF component 1.
    corr_wrr_1 &
    = component_corr_whm_ip( rr1, corr_wrr, rc1, cloud_frac1, rr_tol, &
                             corr_wrr_NL_cloud, corr_wrr_NL_below )

    ! Correlation (in-precip) between w and r_r in PDF component 2.
    corr_wrr_2 &
    = component_corr_whm_ip( rr2, corr_wrr, rc2, cloud_frac2, rr_tol, &
                             corr_wrr_NL_cloud, corr_wrr_NL_below )

    ! Correlation (in-precip) between w and N_r in PDF component 1.
    corr_wNr_1 &
    = component_corr_whm_ip( Nr1, corr_wNr, rc1, cloud_frac1, Nr_tol, &
                             corr_wNr_NL_cloud, corr_wNr_NL_below )

    ! Correlation (in-precip) between w and N_r in PDF component 2.
    corr_wNr_2 &
    = component_corr_whm_ip( Nr2, corr_wNr, rc2, cloud_frac2, Nr_tol, &
                             corr_wNr_NL_cloud, corr_wNr_NL_below )

    ! Correlation between w and N_cn in PDF component 1.
    corr_wNcn_1 &
    = component_corr_whm_ip( Ncnm, corr_wNcn, rc1, one, Ncn_tol, &
                             corr_wNcn_NL_cloud, corr_wNcn_NL_cloud )

    ! Correlation between w and N_cn in PDF component 2.
    corr_wNcn_2 &
    = component_corr_whm_ip( Ncnm, corr_wNcn, rc2, one, Ncn_tol, &
                             corr_wNcn_NL_cloud, corr_wNcn_NL_cloud )

    ! In order to decompose the correlation matrix that we may or may not make
    ! (l_use_modified_corr), we must not have a perfect correlation between s
    ! and t. Thus, we impose a limitation.
    if ( l_use_modified_corr ) then
      l_limit_corr_st = .true.
    else
      l_limit_corr_st = .false.
    end if

    ! Correlation between s and t in PDF component 1.
    corr_st_1 &
    = component_corr_st( pdf_params%corr_st_1, rc1, cloud_frac1, &
                         corr_st_NN_cloud, corr_st_NN_below, &
                         l_limit_corr_st )

    ! Correlation between s and t in PDF component 2.
    corr_st_2 &
    = component_corr_st( pdf_params%corr_st_2, rc2, cloud_frac2, &
                         corr_st_NN_cloud, corr_st_NN_below, &
                         l_limit_corr_st )

    ! Correlation (in-precip) between s and r_r in PDF component 1.
    corr_srr_1 &
    = component_corr_xhm_ip( rr1, rc1, cloud_frac1, rr_tol, &
                             corr_srr_NL_cloud, corr_srr_NL_below )

    ! Correlation (in-precip) between s and r_r in PDF component 2.
    corr_srr_2 &
    = component_corr_xhm_ip( rr2, rc2, cloud_frac2, rr_tol, &
                             corr_srr_NL_cloud, corr_srr_NL_below )

    ! Correlation (in-precip) between s and N_r in PDF component 1.
    corr_sNr_1 &
    = component_corr_xhm_ip( Nr1, rc1, cloud_frac1, Nr_tol, &
                             corr_sNr_NL_cloud, corr_sNr_NL_below )

    ! Correlation (in-precip) between s and N_r in PDF component 2.
    corr_sNr_2 &
    = component_corr_xhm_ip( Nr2, rc2, cloud_frac2, Nr_tol, &
                             corr_sNr_NL_cloud, corr_sNr_NL_below )

    ! Correlation between s and N_cn in PDF component 1.
    corr_sNcn_1 &
    = component_corr_xhm_ip( Ncnm, rc1, one, Ncn_tol, &
                             corr_sNcn_NL_cloud, corr_sNcn_NL_cloud )

    ! Correlation between s and N_cn in PDF component 2.
    corr_sNcn_2 &
    = component_corr_xhm_ip( Ncnm, rc2, one, Ncn_tol, &
                             corr_sNcn_NL_cloud, corr_sNcn_NL_cloud )

    if ( l_use_modified_corr ) then

      ! Correlation (in-precip) between t and r_r in PDF component 1.
      corr_trr_1 &
      = component_corr_thm_ip( corr_st_1, corr_srr_1 )

      ! Correlation (in-precip) between t and r_r in PDF component 2.
      corr_trr_2 &
      = component_corr_thm_ip( corr_st_2, corr_srr_2 )

      ! Correlation (in-precip) between t and N_r in PDF component 1.
      corr_tNr_1 &
      = component_corr_thm_ip( corr_st_1, corr_sNr_1 )

      ! Correlation (in-precip) between t and N_r in PDF component 2.
      corr_tNr_2 &
      = component_corr_thm_ip( corr_st_2, corr_sNr_2 )

    else ! .not. l_use_modified_corr

      ! Correlation (in-precip) between t and r_r in PDF component 1.
      corr_trr_1 &
      = component_corr_xhm_ip( rr1, rc1, cloud_frac1, rr_tol, &
                               corr_trr_NL_cloud, corr_trr_NL_below )

      ! Correlation (in-precip) between t and r_r in PDF component 2.
      corr_trr_2 &
      = component_corr_xhm_ip( rr2, rc2, cloud_frac2, rr_tol, &
                               corr_trr_NL_cloud, corr_trr_NL_below )

      ! Correlation (in-precip) between t and N_r in PDF component 1.
      corr_tNr_1 &
      = component_corr_xhm_ip( Nr1, rc1, cloud_frac1, Nr_tol, &
                               corr_tNr_NL_cloud, corr_tNr_NL_below )

      ! Correlation (in-precip) between t and N_r in PDF component 2.
      corr_tNr_2 &
      = component_corr_xhm_ip( Nr2, rc2, cloud_frac2, Nr_tol, &
                               corr_tNr_NL_cloud, corr_tNr_NL_below )

    end if ! l_use_modified_corr

    ! Correlation between t and N_cn in PDF component 1.
    corr_tNcn_1 &
    = component_corr_xhm_ip( Ncnm, rc1, one, Ncn_tol, &
                             corr_tNcn_NL_cloud, corr_tNcn_NL_cloud )

    ! Correlation between t and N_cn in PDF component 2.
    corr_tNcn_2 &
    = component_corr_xhm_ip( Ncnm, rc2, one, Ncn_tol, &
                             corr_tNcn_NL_cloud, corr_tNcn_NL_cloud )

    ! Correlation (in-precip) between r_r and N_r in PDF component 1.
    corr_rrNr_1 &
    = component_corr_hmxhmy_ip( rr1, Nr1, rc1, cloud_frac1, &
                                rr_tol, Nr_tol, &
                                corr_rrNr_LL_cloud, &
                                corr_rrNr_LL_below )

    ! Correlation (in-precip) between r_r and N_r in PDF component 2.
    corr_rrNr_2 &
    = component_corr_hmxhmy_ip( rr2, Nr2, rc2, cloud_frac2, &
                                rr_tol, Nr_tol, &
                                corr_rrNr_LL_cloud, &
                                corr_rrNr_LL_below )


    return

  end subroutine compute_corr

  !=============================================================================
  function component_mean_hm_ip( hmi, precip_frac_i, hm_tol )  &
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

    use model_flags, only: &
        l_use_hydromet_tolerance

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmi,           & ! Mean of hydrometeor, hm (ith PDF component) [hm units]
      precip_frac_i, & ! Precipitation fraction (ith PDF component)  [-]
      hm_tol           ! Tolerance value for hydrometeor             [hm units]

    ! Return Variable
    real( kind = core_rknd ) :: &
      mu_hm_i    ! Mean of hm (ith PDF component) in-precip (ip)     [hm units]


    ! Mean of the hydrometeor (in-precip) in the ith PDF component.
    if ( ( hmi > hm_tol ) .or. ( .not. l_use_hydromet_tolerance ) ) then
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
  function component_stdev_hm_ip( hmi, mu_hm_i, rci, cloud_fraci, hm_tol, &
                                  hmp2_on_hmm2_cloud, hmp2_on_hmm2_below )  &
  result( sigma_hm_i )

    ! Description:
    ! Calculates the in-precip standard deviation of a hydrometeor species
    ! within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        zero,   &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_use_hydromet_tolerance

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmi,         & ! Mean of hydrometeor, hm (ith PDF component)   [hm units]
      mu_hm_i,     & ! Mean of hm (ith PDF component) in-precip (ip) [hm units]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci, & ! Cloud fraction (ith PDF component)            [-]
      hm_tol         ! Tolerance value for hydrometeor               [hm units]

    real( kind = core_rknd ), intent(in) :: &
      hmp2_on_hmm2_cloud, & ! Ratio of <hm'^2> to <hm>^2: cloudy grid levels [-]
      hmp2_on_hmm2_below    ! Ratio of <hm'^2> to <hm>^2: clear grid levels  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      sigma_hm_i    ! Standard deviation of hm (ith PDF component) ip [hm units]


    ! Standard deviation of the hydrometeor (in-precip) in the
    ! ith PDF component.
    if ( ( hmi > hm_tol ) .or. ( .not. l_use_hydromet_tolerance ) ) then
       if ( l_interp_prescribed_params ) then
          sigma_hm_i = sqrt( cloud_fraci * hmp2_on_hmm2_cloud &
                             + ( one - cloud_fraci ) * hmp2_on_hmm2_below ) &
                       * mu_hm_i
       else
          if ( rci > rc_tol ) then
             sigma_hm_i = sqrt( hmp2_on_hmm2_cloud ) * mu_hm_i
          else
             sigma_hm_i = sqrt( hmp2_on_hmm2_below ) * mu_hm_i
          endif
       endif
    else
       ! The mean of the hydrometeor in the ith PDF component is less than the
       ! tolerance amount for the particular hydrometeor.  It is considered to
       ! have a value of 0.  There is not any of this hydrometeor species in the
       ! ith PDF component at this grid level.  The standard deviation of the
       ! hydrometeor is simply 0 since the hydrometeor does not vary in this
       ! component at this grid level.
       sigma_hm_i = zero
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
  function component_corr_whm_ip( hmi, corr_whm, rci, cloud_fraci, hm_tol, &
                                  corr_whm_NL_cloud, corr_whm_NL_below ) &
  result( corr_whm_i )

    ! Description:
    ! Calculates the in-precip correlation between w and a hydrometeor species
    ! within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        zero,   &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_calc_w_corr, &
        l_use_hydromet_tolerance

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmi,         & ! Mean of hydrometeor, hm (ith PDF component)   [hm units]
      corr_whm,    & ! Correlation between w and hm (overall) ip     [-]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci, & ! Cloud fraction (ith PDF component)            [-]
      hm_tol         ! Tolerance value for hydrometeor               [hm units]

    real( kind = core_rknd ), intent(in) :: &
      corr_whm_NL_cloud, & ! Corr. btwn. w and hm (ith PDF comp.) ip; cloudy [-]
      corr_whm_NL_below    ! Corr. btwn. w and hm (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_whm_i    ! Correlation between w and hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) between w and the hydrometeor in the ith
    ! PDF component.
    if ( ( hmi > hm_tol ) .or. ( .not. l_use_hydromet_tolerance ) ) then
       if ( l_calc_w_corr ) then
          corr_whm_i = corr_whm
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
    else
       ! The mean of the hydrometeor in the ith PDF component is less than the
       ! tolerance amount for the particular hydrometeor.  It is considered to
       ! have a value of 0.  There is not any of this hydrometeor species in the
       ! ith PDF component at this grid level.  Both the standard deviation of
       ! the hydrometeor and the covariance of w and the hydrometeor are simply
       ! 0 since the hydrometeor does not vary in this component at this grid
       ! level.  The correlation involving w and the hydrometeor in the ith PDF
       ! component is undefined (and will be set to 0) since the hydrometeor
       ! does not vary in this component at this grid level.
       corr_whm_i = zero
    endif


    return

  end function component_corr_whm_ip

  !=============================================================================
  function component_corr_xhm_ip( hmi, rci, cloud_fraci, hm_tol, &
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
        zero,   &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_use_hydromet_tolerance

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmi,         & ! Mean of hydrometeor, hm (ith PDF component)   [hm units]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci, & ! Cloud fraction (ith PDF component)            [-]
      hm_tol         ! Tolerance value for hydrometeor               [hm units]

    real( kind = core_rknd ), intent(in) :: &
      corr_xhm_NL_cloud, & ! Corr. btwn. x and hm (ith PDF comp.) ip; cloudy [-]
      corr_xhm_NL_below    ! Corr. btwn. x and hm (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_xhm_i    ! Correlation between x and hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) between x and the hydrometeor in the ith
    ! PDF component.
    if ( ( hmi > hm_tol ) .or. ( .not. l_use_hydromet_tolerance ) ) then
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
    else
       ! The mean of the hydrometeor in the ith PDF component is less than the
       ! tolerance amount for the particular hydrometeor.  It is considered to
       ! have a value of 0.  There is not any of this hydrometeor species in the
       ! ith PDF component at this grid level.  Both the standard deviation of
       ! the hydrometeor and the covariance of x and the hydrometeor are simply
       ! 0 since the hydrometeor does not vary in this component at this grid
       ! level.  The correlation involving x and the hydrometeor in the ith PDF
       ! component is undefined (and will be set to 0) since the hydrometeor
       ! does not vary in this component at this grid level.
       corr_xhm_i = zero
    endif


    return

  end function component_corr_xhm_ip

  !=============================================================================
  function component_corr_hmxhmy_ip( hmxi, hmyi, rci, cloud_fraci, &
                                     hmx_tol, hmy_tol, &
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
        zero,   &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_use_hydromet_tolerance

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmxi,        & ! Mean of hydrometeor x, hmx (ith PDF comp.)    [hmx units]
      hmyi,        & ! Mean of hydrometeor y, hmy (ith PDF comp.)    [hmy units]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci, & ! Cloud fraction (ith PDF component)            [-]
      hmx_tol,     & ! Tolerance value for hydrometeor x             [hmx units]
      hmy_tol        ! Tolerance value for hydrometeor y             [hmy units]

    real( kind = core_rknd ), intent(in) :: &
      corr_hmxhmy_LL_cloud, & ! Corr.: hmx & hmy (ith PDF comp.) ip; cloudy [-]
      corr_hmxhmy_LL_below    ! Corr.: hmx & hmy (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_hmxhmy_i   ! Correlation between hmx & hmy (ith PDF component) ip [-]


    ! Correlation (in-precip) between hydrometeor x and hydrometeor y in the
    ! ith PDF component.
    if ( ( ( hmxi > hmx_tol ) .and. ( hmyi > hmy_tol ) ) &
         .or. ( .not. l_use_hydromet_tolerance ) ) then
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
    else
       ! The mean(s) of both/either hydrometeor x and/or hydrometeor y in the
       ! ith PDF component are/is less than the tolerance amount for the
       ! particular hydrometeor(s), and are/is considered to have a value of 0.
       ! There are/is not any of these/this hydrometeor species in the ith PDF
       ! component at this grid level.  Both the standard deviation of the
       ! hydrometeor(s) and the covariance of hydrometeor x and hydrometeor y
       ! are simply 0 since the hydrometeor(s) do/does not vary in this
       ! component at this grid level.  The correlation involving hydrometeor x
       ! hydrometeor y in the ith PDF component is undefined (and will be set
       ! to 0) since the hydrometeor(s) do/does not vary in this component at
       ! this grid level.
       corr_hmxhmy_i = zero
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
  subroutine normalize_pdf_params( rr1, rr2, Nr1, Nr2, Ncnm, &                   ! Intent(in)
                                   mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &         ! Intent(in)
                                   mu_Ncn_1, mu_Ncn_2, sigma_rr_1, sigma_rr_2, & ! Intent(in)
                                   sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, &        ! Intent(in)
                                   sigma_Ncn_2, corr_wrr_1, corr_wrr_2, &        ! Intent(in)
                                   corr_wNr_1, corr_wNr_2, corr_wNcn_1, &        ! Intent(in)
                                   corr_wNcn_2, corr_srr_1, corr_srr_2, &        ! Intent(in)
                                   corr_sNr_1, corr_sNr_2, corr_sNcn_1, &        ! Intent(in)
                                   corr_sNcn_2, corr_trr_1, corr_trr_2, &        ! Intent(in)
                                   corr_tNr_1, corr_tNr_2, corr_tNcn_1, &        ! Intent(in)
                                   corr_tNcn_2, corr_rrNr_1, corr_rrNr_2, &      ! Intent(in)
                                   mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(out)
                                   mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(out)
                                   sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(out)
                                   sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(out)
                                   corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(out)
                                   corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(out)
                                   corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(out)
                                   corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(out)
                                   corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(out)
                                   corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(out)
                                   corr_rrNr_1_n, corr_rrNr_2_n )                ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rr_tol, & ! Constant(s)
        Nr_tol, &
        Ncn_tol, &
        zero

    use PDF_utilities, only: &
        mean_L2N,   & ! Procedure(s)
        stdev_L2N,  &
        corr_NL2NN, &
        corr_LL2NN

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rr1,         & ! Mean rain water mixing ratio (1st PDF component)  [kg/kg]
      rr2,         & ! Mean rain water mixing ratio (2nd PDF component)  [kg/kg]
      Nr1,         & ! Mean rain drop concentration (1st PDF component) [num/kg]
      Nr2,         & ! Mean rain drop concentration (2nd PDF component) [num/kg]
      Ncnm,        & ! Mean cloud nuclei concentration (overall)        [num/kg]
      mu_rr_1,     & ! Mean of rr (1st PDF component) in-precip (ip)     [kg/kg]
      mu_rr_2,     & ! Mean of rr (2nd PDF component) ip                 [kg/kg]
      mu_Nr_1,     & ! Mean of Nr (1st PDF component) ip                [num/kg]
      mu_Nr_2,     & ! Mean of Nr (2nd PDF component) ip                [num/kg]
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_rr_1,  & ! Standard deviation of rr (1st PDF component) ip   [kg/kg]
      sigma_rr_2,  & ! Standard deviation of rr (2nd PDF component) ip   [kg/kg]
      sigma_Nr_1,  & ! Standard deviation of Nr (1st PDF component) ip  [num/kg]
      sigma_Nr_2,  & ! Standard deviation of Nr (2nd PDF component) ip  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2, & ! Standard deviation of Ncn (2nd PDF component)    [num/kg]
      corr_wrr_1,  & ! Correlation between w and rr (1st PDF component) ip   [-]
      corr_wrr_2,  & ! Correlation between w and rr (2nd PDF component) ip   [-]
      corr_wNr_1,  & ! Correlation between w and Nr (1st PDF component) ip   [-]
      corr_wNr_2,  & ! Correlation between w and Nr (2nd PDF component) ip   [-]
      corr_wNcn_1, & ! Correlation between w and Ncn (1st PDF component)     [-]
      corr_wNcn_2, & ! Correlation between w and Ncn (2nd PDF component)     [-]
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component) ip   [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component) ip   [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component) ip   [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component) ip   [-]
      corr_sNcn_1, & ! Correlation between s and Ncn (1st PDF component)     [-]
      corr_sNcn_2, & ! Correlation between s and Ncn (2nd PDF component)     [-]
      corr_trr_1,  & ! Correlation between t and rr (1st PDF component) ip   [-]
      corr_trr_2,  & ! Correlation between t and rr (2nd PDF component) ip   [-]
      corr_tNr_1,  & ! Correlation between t and Nr (1st PDF component) ip   [-]
      corr_tNr_2,  & ! Correlation between t and Nr (2nd PDF component) ip   [-]
      corr_tNcn_1, & ! Correlation between t and Ncn (1st PDF component)     [-]
      corr_tNcn_2, & ! Correlation between t and Ncn (2nd PDF component)     [-]
      corr_rrNr_1, & ! Correlation between rr & Nr (1st PDF component) ip    [-]
      corr_rrNr_2    ! Correlation between rr & Nr (2nd PDF component) ip    [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
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
      sigma_Ncn_2_n, & ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]
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


    ! --- Begin Code ---

    !!! Calculate the normalized mean of variables that have an assumed
    !!! lognormal distribution, given the mean and variance of those variables.

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 1.
    if ( rr1 > rr_tol ) then
       mu_rr_1_n = mean_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! value of mu_rr_1_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_rr_1_n = -huge( mu_rr_1_n )
    endif

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 2.
    if ( rr2 > rr_tol ) then
       mu_rr_2_n = mean_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! value of mu_rr_2_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_rr_2_n = -huge( mu_rr_2_n )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       mu_Nr_1_n = mean_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! value of mu_Nr_1_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_Nr_1_n = -huge( mu_Nr_1_n )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       mu_Nr_2_n = mean_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! value of mu_Nr_2_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_Nr_2_n = -huge( mu_Nr_2_n )
    endif

    ! Normalized mean of cloud nuclei concentration in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       mu_Ncn_1_n = mean_L2N( mu_Ncn_1, sigma_Ncn_1**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The value of mu_Ncn_1_n should be -inf.  It
       ! will be set to -huge for purposes of assigning it a value.
       mu_Ncn_1_n = -huge( mu_Ncn_1_n )
    endif

    ! Normalized mean of cloud nuclei concentration in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       mu_Ncn_2_n = mean_L2N( mu_Ncn_2, sigma_Ncn_2**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The value of mu_Ncn_2_n should be -inf.  It
       ! will be set to -huge for purposes of assigning it a value.
       mu_Ncn_2_n = -huge( mu_Ncn_2_n )
    endif

    !!! Calculate the normalized standard deviation of variables that have
    !!! an assumed lognormal distribution, given the mean and variance of
    !!! those variables.

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       sigma_rr_1_n = stdev_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_rr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       sigma_rr_2_n = stdev_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_rr_2_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       sigma_Nr_1_n = stdev_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_Nr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       sigma_Nr_2_n = stdev_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_Nr_2_n = zero
    endif

    ! Normalized standard deviation of cloud nuclei concentration
    ! in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       sigma_Ncn_1_n = stdev_L2N( mu_Ncn_1, sigma_Ncn_1**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The standard deviation is simply 0 since
       ! cloud nuclei concentration does not vary at this grid level.
       sigma_Ncn_1_n = zero
    endif

    ! Normalized standard deviation of cloud nuclei concentration
    ! in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       sigma_Ncn_2_n = stdev_L2N( mu_Ncn_2, sigma_Ncn_2**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The standard deviation is simply 0 since
       ! cloud nuclei concentration does not vary at this grid level.
       sigma_Ncn_2_n = zero
    endif

    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    ! Normalize the correlation (in-precip) between w and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_wrr_1_n = corr_NL2NN( corr_wrr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_wrr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between w and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_wrr_2_n = corr_NL2NN( corr_wrr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_wrr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between w and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_wNr_1_n = corr_NL2NN( corr_wNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_wNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between w and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_wNr_2_n = corr_NL2NN( corr_wNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_wNr_2_n = zero
    endif

    ! Normalize the correlation between w and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_wNcn_1_n = corr_NL2NN( corr_wNcn_1, sigma_Ncn_1_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_1_n = zero
    endif

    ! Normalize the correlation between w and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_wNcn_2_n = corr_NL2NN( corr_wNcn_2, sigma_Ncn_2_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_srr_1_n = corr_NL2NN( corr_srr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_srr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_srr_2_n = corr_NL2NN( corr_srr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_srr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_sNr_1_n = corr_NL2NN( corr_sNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_sNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_sNr_2_n = corr_NL2NN( corr_sNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_sNr_2_n = zero
    endif

    ! Normalize the correlation between s and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_sNcn_1_n = corr_NL2NN( corr_sNcn_1, sigma_Ncn_1_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_sNcn_1_n = zero
    endif

    ! Normalize the correlation between s and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_sNcn_2_n = corr_NL2NN( corr_sNcn_2, sigma_Ncn_2_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_sNcn_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_trr_1_n = corr_NL2NN( corr_trr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_trr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_trr_2_n = corr_NL2NN( corr_trr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_trr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_tNr_1_n = corr_NL2NN( corr_tNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_tNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_tNr_2_n = corr_NL2NN( corr_tNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_tNr_2_n = zero
    endif

    ! Normalize the correlation between t and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_tNcn_1_n = corr_NL2NN( corr_tNcn_1, sigma_Ncn_1_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_tNcn_1_n = zero
    endif

    ! Normalize the correlation between t and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_tNcn_2_n = corr_NL2NN( corr_tNcn_2, sigma_Ncn_2_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_tNcn_2_n = zero
    endif

    !!! Calculate the normalized correlation between two variables that both
    !!! have an assumed lognormal distribution, given their correlation and both
    !!! of their normalized standard deviations.

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 1.
    if ( rr1 > rr_tol .and. Nr1 > Nr_tol ) then
       corr_rrNr_1_n = corr_LL2NN( corr_rrNr_1, sigma_rr_1_n, sigma_Nr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 and (or) mean rain drop
       ! concentration in PDF component 1 are (is) less than their (its)
       ! respective tolerance amount(s), and are (is) considered to have a value
       ! of 0.  There is not any rain at this grid level.  The mean in-precip
       ! rain water mixing ratio (1st PDF component) and (or) mean in-precip
       ! rain drop concentration (1st PDF component) are also considered to have
       ! a value of 0.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 2.
    if ( rr2 > rr_tol .and. Nr2 > Nr_tol ) then
       corr_rrNr_2_n = corr_LL2NN( corr_rrNr_2, sigma_rr_2_n, sigma_Nr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 and (or) mean rain drop
       ! concentration in PDF component 2 are (is) less than their (its)
       ! respective tolerance amount(s), and are (is) considered to have a value
       ! of 0.  There is not any rain at this grid level.  The mean in-precip
       ! rain water mixing ratio (2nd PDF component) and (or) mean in-precip
       ! rain drop concentration (2nd PDF component) are also considered to have
       ! a value of 0.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_2_n = zero
    endif


    return

  end subroutine normalize_pdf_params

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
  subroutine init_precip_hm_arrays( hydromet_dim )

    ! Description:
    ! Precipitating hydrometeor arrays are hydrometeor arrays made out of only
    ! the precipitating hydrometeors, which do not include N_c or N_cn.  CLUBB's
    ! PDF does not include N_c, but rather N_cn.  However, N_cn is handled in
    ! the PDF in a different manner than the rest of the hydrometeors.  The
    ! remaining hydrometeors (the precipitating hydrometeors -- rain, ice, snow,
    ! etc.) need their own array for use in the PDF parameter setup.

    ! References:
    !-----------------------------------------------------------------------

    use array_index, only: & 
        iirrainm,    & ! Variables
        iiNrm,       &
        iirsnowm,    &
        iiricem,     &
        iirgraupelm, &
        iiNcnm,      &
        iiNsnowm,    & 
        iiNim,       &
        iiNgraupelm, &
        iiNcm          ! Note: Ncm is not part of CLUBB's PDF.

    use clubb_precision, only: &
        core_rknd     ! Variable(s)

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      hydromet_dim    ! Number of hydrometeor fields total

    ! Local Variables
    integer :: &
      idx, & ! Precipitating hydrometeor index
      i      ! Hydrometeor index


    ! Initialize indices of precipitating hydrometeors to -1.
    iirr = -1
    iiNr = -1
    iirs = -1
    iiri = -1
    iirg = -1
    iiNs = -1
    iiNi = -1
    iiNg = -1

    ! Initialize precipitating hydrometeor index to 0.
    idx = 0

    ! Count the number of precipitating hydrometeors and set indices for
    ! precipitating hydrometeors.
    do i = 1, hydromet_dim, 1
       if ( ( i /= iiNcm ) .and. ( i /= iiNcnm ) &
            .and. ( i /= iirgraupelm ) .and. ( i /= iiNgraupelm ) ) then

          ! The hydrometeor is a precipitating hydrometeor, which does not
          ! include N_c or N_cn.
          ! This currently, for use in setting up the PDF arrays, does not
          ! include graupel.

          ! Iterate the precipitating hydrometeor index.
          idx = idx + 1

          ! Assign indices to each hydrometeor in the precipitating hydrometeor
          ! arrays.
          if ( i == iirrainm ) then
             iirr = idx
          elseif ( i == iiNrm ) then
             iiNr = idx
          elseif ( i == iirsnowm ) then
             iirs = idx
          elseif ( i == iiricem ) then
             iiri = idx
!          elseif ( i == iirgraupelm ) then
!             iirg = idx
          elseif ( i == iiNsnowm ) then
             iiNs = idx
          elseif ( i == iiNim ) then
             iiNi = idx
!          elseif ( i == iiNgraupelm ) then
!             iiNg = idx
          endif

       endif ! i /= iiNcm and i /= iiNcnm

    enddo ! i = 1, hydromet_dim, 1

    ! Total number of precipitating hydrometeors
    num_hm = idx


    return

  end subroutine init_precip_hm_arrays

  !=============================================================================
  subroutine precip_hm_arrays( nz, hydromet, wphydrometp, &
                               hmm, wphmp, hm_list, hm_tol )

    ! Description:
    ! Makes new hydrometeor arrays out of only the precipitating hydrometeors,
    ! which do not include N_c or N_cn.  CLUBB's PDF does not include N_c, but
    ! rather N_cn.  However, N_cn is handled in the PDF in a different manner
    ! than the rest of the hydrometeors.  The remaining hydrometeors (the
    ! precipitating hydrometeors -- rain, ice, snow, etc.) need their own array
    ! for use in the PDF parameter setup.

    ! References:
    !-----------------------------------------------------------------------

    use array_index, only: & 
        iiNcnm, & ! Variable(s)
        iiNcm     ! Note: Ncm is not part of CLUBB's PDF.

    use array_index, only: & 
        iirgraupelm, & ! Variable(s)
        iiNgraupelm    ! Note: graupel is not currently part of CLUBB's PDF.

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use parameters_microphys, only: &
        hydromet_list, & ! Variable(s)
        hydromet_tol

    use clubb_precision, only: &
        core_rknd     ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,num_hm), intent(out) :: &
      hmm,   & ! Mean of a precipitating hydrometeor, hm (overall) [units vary]
      wphmp    ! Covariance of w and a precipitating hydrometeor   [(m/s)units]

    character(len=10), dimension(num_hm), intent(out) :: &
      hm_list    ! Names of all precipitating hydrometeors

    real( kind = core_rknd ), dimension(num_hm), intent(out) :: &
      hm_tol    ! Tolerance value for a precipitating hydrometeor  [units vary]

    ! Local Variables
    integer :: &
      idx, & ! Precipitating hydrometeor index
      i      ! Hydrometeor index

    ! Initialize precipitating hydrometeor index to 0.
    idx = 0

    ! Loop over all hydrometeors and set up precipitating hydrometeor arrays.
    do i = 1, hydromet_dim, 1

       if ( ( i /= iiNcm ) .and. ( i /= iiNcnm ) &
            .and. ( i /= iirgraupelm ) .and. ( i /= iiNgraupelm ) ) then

          ! The hydrometeor is a precipitating hydrometeor, which does not
          ! include N_c or N_cn.
          ! This currently, for use in setting up the PDF arrays, does not
          ! include graupel.

          ! Iterate the precipitating hydrometeor index.
          idx = idx + 1

          ! Mean of hydrometeor (overall)
          hmm(:,idx) = hydromet(:,i)

          ! Covariance of vertical velocity and a hydrometeor
          wphmp(:,idx) = wphydrometp(:,i)

          ! Hydrometeor name from hydrometeor list
          hm_list(idx) = hydromet_list(i)

          ! Hydrometeor tolerance values
          hm_tol(idx) = hydromet_tol(i)

       endif ! i /= iiNcm and i /= iiNcnm

    enddo ! i = 1, hydromet_dim, 1


    return

  end subroutine precip_hm_arrays

  !=============================================================================
  subroutine pack_pdf_params( mu_w_1, mu_w_2, mu_s_1, mu_s_2, &             ! Intent(in)
                              mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &           ! Intent(in)
                              mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &       ! Intent(in)
                              mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(in)
                              mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(in)
                              sigma_w_1, sigma_w_2, sigma_s_1, &            ! Intent(in)
                              sigma_s_2, sigma_t_1, sigma_t_2, &            ! Intent(in)
                              sigma_rr_1, sigma_rr_2, sigma_Nr_1, &         ! Intent(in)
                              sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &       ! Intent(in)
                              sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(in)
                              sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(in)
                              corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, & ! Intent(in)
                              corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(in)
                              corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(in)
                              corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(in)
                              corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(in)
                              corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(in)
                              corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(in)
                              corr_rrNr_1_n, corr_rrNr_2_n, &               ! Intent(in)
                              rr1, rr2, Nr1, Nr2, &                         ! Intent(in)
                              precip_frac, precip_frac_1, precip_frac_2, &  ! Intent(in)
                              d_variables, &                                ! Intent(in)
                              corr_array_1, corr_array_2, &                 ! Intent(out)
                              mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &       ! Intent(out)
                              hydromet_pdf_params )                         ! Intent(out)

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
        core_rknd    ! Constant

    use constants_clubb, only: &
        zero, &      ! Constant(s)
        one

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
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

    real( kind = core_rknd ), intent(in) :: &
      corr_ws_1,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_ws_2,     & ! Correlation between w and s (2nd PDF component)     [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2        ! Correlation between s and t (2nd PDF component)     [-]

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

    real( kind = core_rknd ), intent(in) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]


    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array.

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(out) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(hydromet_pdf_parameter), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local variables
    integer :: i ! Loop variable

    ! ---- Begin Code ----
    ! Initialize output variables
    mu_x_1 = zero
    mu_x_2 = zero
    sigma_x_1 = zero
    sigma_x_2 = zero
    corr_array_1 = zero
    corr_array_2 = zero

    ! Pack Means and Standard Deviations into arrays
    mu_x_1(iiPDF_w)        = mu_w_1
    mu_x_2(iiPDF_w)        = mu_w_2
    mu_x_1(iiPDF_s_mellor) = mu_s_1
    mu_x_2(iiPDF_s_mellor) = mu_s_2
    mu_x_1(iiPDF_t_mellor) = mu_t_1
    mu_x_2(iiPDF_t_mellor) = mu_t_2
    mu_x_1(iiPDF_rrain)    = mu_rr_1_n
    mu_x_2(iiPDF_rrain)    = mu_rr_2_n
    mu_x_1(iiPDF_Nr)       = mu_Nr_1_n
    mu_x_2(iiPDF_Nr)       = mu_Nr_2_n
    mu_x_1(iiPDF_Ncn)      = mu_Ncn_1_n
    mu_x_2(iiPDF_Ncn)      = mu_Ncn_2_n
    sigma_x_1(iiPDF_w)        = sigma_w_1
    sigma_x_2(iiPDF_w)        = sigma_w_2
    sigma_x_1(iiPDF_s_mellor) = sigma_s_1
    sigma_x_2(iiPDF_s_mellor) = sigma_s_2
    sigma_x_1(iiPDF_t_mellor) = sigma_t_1
    sigma_x_2(iiPDF_t_mellor) = sigma_t_2
    sigma_x_1(iiPDF_rrain)    = sigma_rr_1_n
    sigma_x_2(iiPDF_rrain)    = sigma_rr_2_n
    sigma_x_1(iiPDF_Nr)       = sigma_Nr_1_n
    sigma_x_2(iiPDF_Nr)       = sigma_Nr_2_n
    sigma_x_1(iiPDF_Ncn)      = sigma_Ncn_1_n
    sigma_x_2(iiPDF_Ncn)      = sigma_Ncn_2_n

    ! Pack remaining variables into hydromet_pdf_params
    hydromet_pdf_params%mu_rr_1     = mu_rr_1
    hydromet_pdf_params%mu_rr_2     = mu_rr_2
    hydromet_pdf_params%mu_Nr_1     = mu_Nr_1
    hydromet_pdf_params%mu_Nr_2     = mu_Nr_2
    hydromet_pdf_params%mu_Ncn_1    = mu_Ncn_1
    hydromet_pdf_params%mu_Ncn_2    = mu_Ncn_2
    hydromet_pdf_params%sigma_rr_1  = sigma_rr_1
    hydromet_pdf_params%sigma_rr_2  = sigma_rr_2
    hydromet_pdf_params%sigma_Nr_1  = sigma_Nr_1
    hydromet_pdf_params%sigma_Nr_2  = sigma_Nr_2
    hydromet_pdf_params%sigma_Ncn_1 = sigma_Ncn_1
    hydromet_pdf_params%sigma_Ncn_2 = sigma_Ncn_2

    hydromet_pdf_params%rr1           = rr1
    hydromet_pdf_params%rr2           = rr2
    hydromet_pdf_params%Nr1           = Nr1
    hydromet_pdf_params%Nr2           = Nr2
    hydromet_pdf_params%precip_frac   = precip_frac
    hydromet_pdf_params%precip_frac_1 = precip_frac_1
    hydromet_pdf_params%precip_frac_2 = precip_frac_2


    !!! The corr_arrays are assumed to be lower triangular matrices
    ! Set diagonal elements to 1
    do i=1, d_variables
      corr_array_1(i, i) = one
      corr_array_2(i, i) = one
    end do

    ! Pack correlations (1st PDF component) into corr_array_1.
    corr_array_1(iiPDF_t_mellor, iiPDF_s_mellor) = corr_st_1
    corr_array_1(iiPDF_w,iiPDF_s_mellor)         = corr_ws_1
    corr_array_1(iiPDF_rrain, iiPDF_s_mellor)    = corr_srr_1_n
    corr_array_1(iiPDF_Nr, iiPDF_s_mellor)       = corr_sNr_1_n
    corr_array_1(iiPDF_Ncn, iiPDF_s_mellor)      = corr_sNcn_1_n
    corr_array_1(iiPDF_rrain, iiPDF_t_mellor)    = corr_trr_1_n
    corr_array_1(iiPDF_Nr, iiPDF_t_mellor)       = corr_tNr_1_n
    corr_array_1(iiPDF_Ncn, iiPDF_t_mellor)      = corr_tNcn_1_n
    corr_array_1(iiPDF_rrain, iiPDF_w)           = corr_wrr_1_n
    corr_array_1(iiPDF_Nr, iiPDF_w)              = corr_wNr_1_n
    corr_array_1(iiPDF_Ncn, iiPDF_w)             = corr_wNcn_1_n
    corr_array_1(iiPDF_Nr,iiPDF_rrain)           = corr_rrNr_1_n

    ! Pack correlations (2nd PDF component) into corr_array_2.
    corr_array_2(iiPDF_t_mellor, iiPDF_s_mellor) = corr_st_2
    corr_array_2(iiPDF_w,iiPDF_s_mellor)         = corr_ws_2
    corr_array_2(iiPDF_rrain, iiPDF_s_mellor)    = corr_srr_2_n
    corr_array_2(iiPDF_Nr, iiPDF_s_mellor)       = corr_sNr_2_n
    corr_array_2(iiPDF_Ncn, iiPDF_s_mellor)      = corr_sNcn_2_n
    corr_array_2(iiPDF_rrain, iiPDF_t_mellor)    = corr_trr_2_n
    corr_array_2(iiPDF_Nr, iiPDF_t_mellor)       = corr_tNr_2_n
    corr_array_2(iiPDF_Ncn, iiPDF_t_mellor)      = corr_tNcn_2_n
    corr_array_2(iiPDF_rrain, iiPDF_w)           = corr_wrr_2_n
    corr_array_2(iiPDF_Nr, iiPDF_w)              = corr_wNr_2_n
    corr_array_2(iiPDF_Ncn, iiPDF_w)             = corr_wNcn_2_n
    corr_array_2(iiPDF_Nr,iiPDF_rrain)           = corr_rrNr_2_n

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
